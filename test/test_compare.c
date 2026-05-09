/********************************************************************
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject
* to the following conditions:
*
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*
* Description: Black-box comparison of micro-segment preprocessors.
*              Runs three strategies on the same CAM NGC files:
*                  baseline    — pass-through G1/G2/G3 (no compression)
*                  dp-xyz      — Douglas-Peucker polyline simplification
*                                (Douglas & Peucker, 1973; GIS / CAM staple)
*                  liscio      — full liscio pipeline (arc+bezier+bspline)
*              Reports per strategy: primitive count, max deviation,
*              fit-path length error, elapsed time, throughput.
*              Usage: test_compare [--tol_xyz=N] <ngc>...
* Author:      杨阳 (Yang Yang) <mika-net@outlook.com>
* License:     MIT (SPDX-License-Identifier: MIT)
* Copyright (c) 2026 杨阳 (Yang Yang)
********************************************************************/

#define _USE_MATH_DEFINES
#define _POSIX_C_SOURCE 199309L
#include "liscio/liscio.h"
#include "ngc_parser.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAX_PTS 200000

/* ---------- polyline (output of NGC parse, arcs pre-sampled) ---------- */
typedef struct { double x, y, z; } vec3_t;

typedef struct {
    vec3_t pts[MAX_PTS];
    int    n;
    double raw_length;
} polyline_t;

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* NGC parse bridge: turn each line/arc move into polyline points
 * (arcs sampled into arc_samples sub-segments for fair comparison). */
typedef struct {
    polyline_t *pl;
    int arc_samples;
} pl_bridge_t;

static void pl_on_move(const ngc_move_t *m, void *user)
{
    pl_bridge_t *b = (pl_bridge_t *)user;
    polyline_t *pl = b->pl;
    if (m->is_rapid) return;    /* rapids excluded */

    /* Seed first point. */
    if (pl->n == 0) {
        if (pl->n < MAX_PTS)
            pl->pts[pl->n++] = (vec3_t){m->start.x, m->start.y, m->start.z};
    }
    if (m->type == NGC_MOVE_LINE) {
        double dx = m->end.x - pl->pts[pl->n-1].x;
        double dy = m->end.y - pl->pts[pl->n-1].y;
        double dz = m->end.z - pl->pts[pl->n-1].z;
        pl->raw_length += sqrt(dx*dx + dy*dy + dz*dz);
        if (pl->n < MAX_PTS) pl->pts[pl->n++] = (vec3_t){m->end.x, m->end.y, m->end.z};
        return;
    }
    /* Arc sample — plane-aware. */
    for (int s = 1; s <= b->arc_samples; s++) {
        double t = (double)s / b->arc_samples;
        double ang = t * m->arc_angle;
        /* Build plane-local start-radial u and perp v = n×u. */
        double vx = m->start.x - m->cx;
        double vy = m->start.y - m->cy;
        double vz = m->start.z - m->cz;
        double un = sqrt(vx*vx + vy*vy + vz*vz);
        if (un < 1e-15) continue;
        double ux = vx/un, uy = vy/un, uz = vz/un;
        double wx = m->ny*uz - m->nz*uy;
        double wy = m->nz*ux - m->nx*uz;
        double wz = m->nx*uy - m->ny*ux;
        double sx = m->cx + m->radius * (cos(ang)*ux + sin(ang)*wx);
        double sy = m->cy + m->radius * (cos(ang)*uy + sin(ang)*wy);
        double sz = m->cz + m->radius * (cos(ang)*uz + sin(ang)*wz);
        /* Helical off-plane linear interp (plane basis is planar). */
        if (m->plane == 0 && fabs(uz)+fabs(wz) < 1e-12)
            sz = m->start.z + t * (m->end.z - m->start.z);
        else if (m->plane == 1 && fabs(uy)+fabs(wy) < 1e-12)
            sy = m->start.y + t * (m->end.y - m->start.y);
        else if (m->plane == 2 && fabs(ux)+fabs(wx) < 1e-12)
            sx = m->start.x + t * (m->end.x - m->start.x);
        if (pl->n < MAX_PTS) {
            double dx = sx - pl->pts[pl->n-1].x;
            double dy = sy - pl->pts[pl->n-1].y;
            double dz = sz - pl->pts[pl->n-1].z;
            pl->raw_length += sqrt(dx*dx+dy*dy+dz*dz);
            pl->pts[pl->n++] = (vec3_t){sx,sy,sz};
        }
    }
}

static int load_ngc(const char *path, polyline_t *pl, int arc_samples)
{
    pl->n = 0; pl->raw_length = 0;
    pl_bridge_t br = { pl, arc_samples };
    ngc_stats_t st;
    int rc = ngc_parse_file(path, pl_on_move, &br, &st);
    if (rc != 0) return -1;
    return 0;
}

/* ---------- Strategy A: baseline pass-through ----------
 * Output = input polyline.  Always max_dev=0, len_err=0, ratio=1.0.
 * Serves as the reference point. */
typedef struct {
    int    n_emit;
    double max_dev;
    double out_length;
    double elapsed_sec;
} strat_result_t;

static void strat_baseline(const polyline_t *pl, strat_result_t *r)
{
    double t0 = now_sec();
    r->n_emit = pl->n - 1;   /* each segment is a primitive */
    r->max_dev = 0.0;
    r->out_length = pl->raw_length;
    r->elapsed_sec = now_sec() - t0;
}

/* ---------- Strategy B: Douglas-Peucker 3D ----------
 * Classic 1973 polyline simplification.  Recursively retains points
 * whose perpendicular distance from the current chord exceeds the
 * tolerance.  Output is a LINE-only polyline — no curves — but the
 * compression ratio vs liscio on curved CAM paths is a telling
 * benchmark (polyline simplification is the standard naive approach). */
static double point_to_segment_dist(vec3_t p, vec3_t a, vec3_t b)
{
    double abx = b.x-a.x, aby = b.y-a.y, abz = b.z-a.z;
    double ab2 = abx*abx + aby*aby + abz*abz;
    if (ab2 < 1e-30) {
        double dx=p.x-a.x, dy=p.y-a.y, dz=p.z-a.z;
        return sqrt(dx*dx + dy*dy + dz*dz);
    }
    double apx = p.x-a.x, apy = p.y-a.y, apz = p.z-a.z;
    double t = (apx*abx + apy*aby + apz*abz) / ab2;
    if (t < 0) t = 0; if (t > 1) t = 1;
    double cx = a.x + t*abx, cy = a.y + t*aby, cz = a.z + t*abz;
    double dx = p.x-cx, dy = p.y-cy, dz = p.z-cz;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

static void dp_recurse(const vec3_t *pts, int i0, int i1, double tol, char *keep)
{
    if (i1 <= i0 + 1) return;
    double maxd = 0; int maxi = i0;
    for (int k = i0+1; k < i1; k++) {
        double d = point_to_segment_dist(pts[k], pts[i0], pts[i1]);
        if (d > maxd) { maxd = d; maxi = k; }
    }
    if (maxd > tol) {
        keep[maxi] = 1;
        dp_recurse(pts, i0, maxi, tol, keep);
        dp_recurse(pts, maxi, i1, tol, keep);
    }
}

static void strat_dp(const polyline_t *pl, double tol, strat_result_t *r)
{
    double t0 = now_sec();
    static char keep[MAX_PTS];
    memset(keep, 0, (size_t)pl->n);
    keep[0] = 1; keep[pl->n-1] = 1;
    dp_recurse(pl->pts, 0, pl->n-1, tol, keep);
    int n_kept = 0;
    vec3_t prev = pl->pts[0];
    double out_len = 0, max_d = 0;
    int last_kept = 0;
    for (int i = 0; i < pl->n; i++) {
        if (!keep[i]) continue;
        if (i > 0) {
            /* Evaluate max deviation on the preceding subsegment
             * (last_kept..i) against the straight chord. */
            for (int k = last_kept+1; k < i; k++) {
                double d = point_to_segment_dist(pl->pts[k],
                    pl->pts[last_kept], pl->pts[i]);
                if (d > max_d) max_d = d;
            }
        }
        if (n_kept > 0) {
            double dx = pl->pts[i].x-prev.x, dy = pl->pts[i].y-prev.y, dz = pl->pts[i].z-prev.z;
            out_len += sqrt(dx*dx+dy*dy+dz*dz);
        }
        prev = pl->pts[i];
        last_kept = i;
        n_kept++;
    }
    r->n_emit = n_kept - 1;
    r->max_dev = max_d;
    r->out_length = out_len;
    r->elapsed_sec = now_sec() - t0;
}

/* ---------- Strategy C: liscio ---------- */
typedef struct {
    long n_prim;
    double max_dev;
    double out_length;
} liscio_ctx_local_t;

static liscio_ctx_local_t g_ll;

static void sample_prim_length(const liscio_primitive_t *p, double *plen)
{
    const int N = 40;
    double xs[64], ys[64], zs[64];
    for (int i = 0; i <= N; i++) {
        double t = (double)i / N;
        switch (p->type) {
            case LISCIO_PRIM_LINE:
                xs[i] = p->start.x + t*(p->end.x - p->start.x);
                ys[i] = p->start.y + t*(p->end.y - p->start.y);
                zs[i] = p->start.z + t*(p->end.z - p->start.z);
                break;
            case LISCIO_PRIM_HELIX:    /* same data layout as ARC */
            case LISCIO_PRIM_ARC: {
                double vx=p->start.x-p->cx, vy=p->start.y-p->cy, vz=p->start.z-p->cz;
                double m=sqrt(vx*vx+vy*vy+vz*vz);
                if (m<1e-15) { xs[i]=p->start.x; ys[i]=p->start.y; zs[i]=p->start.z; break; }
                double ux=vx/m, uy=vy/m, uz=vz/m;
                double rx=p->ny*uz - p->nz*uy;
                double ry=p->nz*ux - p->nx*uz;
                double rz=p->nx*uy - p->ny*ux;
                double ang = t * p->arc_angle;
                xs[i] = p->cx + p->radius*(cos(ang)*ux + sin(ang)*rx);
                ys[i] = p->cy + p->radius*(cos(ang)*uy + sin(ang)*ry);
                zs[i] = p->cz + p->radius*(cos(ang)*uz + sin(ang)*rz);
                if (fabs(uz)+fabs(rz) < 1e-12)
                    zs[i] = p->start.z + t*(p->end.z - p->start.z);
                break;
            }
            case LISCIO_PRIM_BEZIER:
                liscio_rbezier_eval(p, t, &xs[i], &ys[i], &zs[i]);
                break;
            case LISCIO_PRIM_SPLINE:
                liscio_bspline_eval(p, t, &xs[i], &ys[i], &zs[i]);
                break;
        }
    }
    *plen = 0;
    for (int i = 1; i <= N; i++) {
        double dx=xs[i]-xs[i-1], dy=ys[i]-ys[i-1], dz=zs[i]-zs[i-1];
        *plen += sqrt(dx*dx+dy*dy+dz*dz);
    }
}

static void ll_cb(const liscio_primitive_t *p, void *u)
{
    (void)u;
    g_ll.n_prim++;
    double plen = 0; sample_prim_length(p, &plen);
    g_ll.out_length += plen;
    double d = 0;
    if (p->type == LISCIO_PRIM_BEZIER) d = p->bezier_max_deviation;
    else if (p->type == LISCIO_PRIM_SPLINE) d = p->spline_max_deviation;
    if (d > g_ll.max_dev) g_ll.max_dev = d;
}

static void strat_liscio(const polyline_t *pl, double tol, strat_result_t *r)
{
    double t0 = now_sec();
    liscio_cfg_t cfg; liscio_cfg_default(&cfg);
    cfg.tol_xyz = tol;
    liscio_ctx_t *ctx = liscio_create(&cfg);
    memset(&g_ll, 0, sizeof(g_ll));
    liscio_set_callback(ctx, ll_cb, NULL);
    for (int i = 1; i < pl->n; i++) {
        liscio_pose_t s = {pl->pts[i-1].x, pl->pts[i-1].y, pl->pts[i-1].z,0,0,0,0,0,0};
        liscio_pose_t e = {pl->pts[i].x,   pl->pts[i].y,   pl->pts[i].z,  0,0,0,0,0,0};
        liscio_add_line(ctx, &s, &e, 1000.0, i);
    }
    liscio_flush(ctx);
    liscio_destroy(ctx);
    r->n_emit = (int)g_ll.n_prim;
    r->max_dev = g_ll.max_dev;
    r->out_length = g_ll.out_length;
    r->elapsed_sec = now_sec() - t0;
}

/* ---------- driver ---------- */
static void print_row(const char *strat, const polyline_t *pl,
                      const strat_result_t *r)
{
    double ratio = r->n_emit > 0 ? (double)(pl->n-1) / r->n_emit : 0;
    double lenerr = pl->raw_length > 1e-12
                    ? (r->out_length - pl->raw_length) / pl->raw_length * 100.0
                    : 0;
    double us_per = pl->n > 0 ? r->elapsed_sec * 1e6 / pl->n : 0;
    printf("  %-10s  segs=%-6d  ratio=%6.2fx  max_dev=%.4f  len_err=%+7.3f%%  %.3f ms  (%.2f µs/pt)\n",
        strat, r->n_emit, ratio, r->max_dev, lenerr,
        r->elapsed_sec * 1000.0, us_per);
}

int main(int argc, char **argv)
{
    double tol = 0.05;
    int arc_samp = 8;
    int i0 = 1;
    for (int i = 1; i < argc; i++) {
        if (strncmp(argv[i], "--tol_xyz=", 10) == 0) {
            tol = atof(argv[i]+10); i0 = i+1;
        } else if (strncmp(argv[i], "--arc_samples=", 14) == 0) {
            arc_samp = atoi(argv[i]+14); i0 = i+1;
        } else break;
    }
    if (i0 >= argc) {
        fprintf(stderr, "usage: %s [--tol_xyz=N] [--arc_samples=N] <ngc>...\n", argv[0]);
        return 1;
    }

    printf("# liscio black-box comparator (tol_xyz=%.4f, arc_samples=%d)\n", tol, arc_samp);
    printf("# Reference: Douglas & Peucker (1973), classic polyline-simplification baseline.\n");
    printf("# baseline = identity (each raw segment emitted as-is; deviation=0)\n");
    printf("# dp-xyz   = Douglas-Peucker 3D with same tolerance\n");
    printf("# liscio   = full pipeline (LINE+ARC+cubic Bezier+B-spline)\n\n");

    static polyline_t pl;
    double grand_t_base = 0, grand_t_dp = 0, grand_t_lc = 0;
    long grand_in = 0, grand_out_base = 0, grand_out_dp = 0, grand_out_lc = 0;

    for (int i = i0; i < argc; i++) {
        if (load_ngc(argv[i], &pl, arc_samp) != 0) {
            fprintf(stderr, "skip %s (read fail)\n", argv[i]);
            continue;
        }
        const char *name = strrchr(argv[i], '/');
        name = name ? name+1 : argv[i];
        printf("%s  (waypoints=%d, raw_len=%.3f)\n", name, pl.n, pl.raw_length);

        strat_result_t rb, rd, rl;
        strat_baseline(&pl, &rb);  print_row("baseline", &pl, &rb);
        strat_dp(&pl, tol, &rd);   print_row("dp-xyz",   &pl, &rd);
        strat_liscio(&pl, tol, &rl); print_row("liscio",  &pl, &rl);
        printf("\n");

        grand_t_base += rb.elapsed_sec;
        grand_t_dp   += rd.elapsed_sec;
        grand_t_lc   += rl.elapsed_sec;
        grand_in     += pl.n - 1;
        grand_out_base += rb.n_emit;
        grand_out_dp   += rd.n_emit;
        grand_out_lc   += rl.n_emit;
    }

    printf("=== GRAND TOTAL ===\n");
    printf("  baseline : segs=%-8ld  ratio=%6.2fx  time=%.3f s\n",
        grand_out_base, grand_out_base ? (double)grand_in/grand_out_base : 0, grand_t_base);
    printf("  dp-xyz   : segs=%-8ld  ratio=%6.2fx  time=%.3f s\n",
        grand_out_dp, grand_out_dp ? (double)grand_in/grand_out_dp : 0, grand_t_dp);
    printf("  liscio   : segs=%-8ld  ratio=%6.2fx  time=%.3f s\n",
        grand_out_lc, grand_out_lc ? (double)grand_in/grand_out_lc : 0, grand_t_lc);
    double speedup_dp = grand_t_dp > 1e-9 ? grand_t_lc / grand_t_dp : 0;
    printf("\n  liscio vs dp-xyz: %.2fx compression advantage, %.2fx time\n",
        grand_out_dp > 0 && grand_out_lc > 0
         ? (double)grand_out_dp / grand_out_lc : 0,
        speedup_dp);
    return 0;
}
