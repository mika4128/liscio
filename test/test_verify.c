/********************************************************************
* SPDX-License-Identifier: MIT
* Copyright (c) 2026 杨阳 (Yang Yang) <mika-net@outlook.com>
*
* Independent fit-quality verifier.  Does NOT trust liscio's self-
* reported max_deviation fields.  For each emitted primitive, takes the
* raw waypoints it absorbed (tracked via n_absorbed), samples the
* primitive densely, computes actual max perpendicular distance.
*
* Reports per file:
*   max_obs_dev  observed max waypoint→fit_curve distance
*   avg_obs_dev  mean
*   OVER_TOL     # waypoints that exceed tol_xyz (should be 0)
*
* Usage: test_verify [--tol_xyz=N] [--samples=N] <ngc>...
********************************************************************/

#define _USE_MATH_DEFINES
#define _POSIX_C_SOURCE 199309L
#include "liscio/liscio.h"
#include "ngc_parser.h"
#include "test_timing.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAX_RAW   300000

typedef struct {
    double x, y, z;
} vec3_t;

/* Accumulated per-file state. */
static vec3_t g_raw[MAX_RAW];
static int    g_raw_n;
static long   g_abs_cursor;      /* raw index consumed by primitives so far */
static double g_obs_max;
static double g_obs_sum;
static long   g_obs_count;
static long   g_over_tol;
static double g_tol;
static int    g_samples;

/* ---------- NGC → raw waypoint array ----------
 * Record one point per NGC move: start-seed on very first feed move, then
 * each move's endpoint.  Matches liscio's n_absorbed accounting (each
 * liscio_add_line / liscio_add_arc call absorbs one segment = 1 unit). */
static void raw_on_move(const ngc_move_t *m, void *user)
{
    (void)user;
    if (m->is_rapid) return;
    if (g_raw_n == 0 && g_raw_n < MAX_RAW) {
        g_raw[g_raw_n++] = (vec3_t){m->start.x, m->start.y, m->start.z};
    }
    if (g_raw_n < MAX_RAW)
        g_raw[g_raw_n++] = (vec3_t){m->end.x, m->end.y, m->end.z};
}

/* Sample primitive at t ∈ [0,1]. */
static void prim_eval(const liscio_primitive_t *p, double t,
                      double *ox, double *oy, double *oz)
{
    switch (p->type) {
        case LISCIO_PRIM_LINE:
            *ox = p->start.x + t*(p->end.x - p->start.x);
            *oy = p->start.y + t*(p->end.y - p->start.y);
            *oz = p->start.z + t*(p->end.z - p->start.z);
            return;
        case LISCIO_PRIM_HELIX:    /* same data layout as ARC + pitch */
        case LISCIO_PRIM_ARC: {
            double vx = p->start.x - p->cx;
            double vy = p->start.y - p->cy;
            double vz = p->start.z - p->cz;
            double un = sqrt(vx*vx + vy*vy + vz*vz);
            if (un < 1e-15) { *ox = p->start.x; *oy = p->start.y; *oz = p->start.z; return; }
            double ux = vx/un, uy = vy/un, uz = vz/un;
            double wx = p->ny*uz - p->nz*uy;
            double wy = p->nz*ux - p->nx*uz;
            double wz = p->nx*uy - p->ny*ux;
            double ang = t * p->arc_angle;
            double axial = ang * p->pitch;     /* HELIX advance; 0 for ARC */
            *ox = p->cx + p->radius*(cos(ang)*ux + sin(ang)*wx) + axial*p->nx;
            *oy = p->cy + p->radius*(cos(ang)*uy + sin(ang)*wy) + axial*p->ny;
            *oz = p->cz + p->radius*(cos(ang)*uz + sin(ang)*wz) + axial*p->nz;
            if (p->type == LISCIO_PRIM_ARC && fabs(uz)+fabs(wz) < 1e-12)
                *oz = p->start.z + t*(p->end.z - p->start.z);
            return;
        }
        case LISCIO_PRIM_BEZIER:
            liscio_rbezier_eval(p, t, ox, oy, oz);
            return;
        case LISCIO_PRIM_SPLINE:
            liscio_bspline_eval(p, t, ox, oy, oz);
            return;
    }
}

/* Feed one NGC move into liscio. */
static void feed_liscio(const ngc_move_t *m, void *user)
{
    liscio_ctx_t *c = (liscio_ctx_t *)user;
    if (m->is_rapid) {
        liscio_pose_t s = {m->start.x,m->start.y,m->start.z,0,0,0,0,0,0};
        liscio_pose_t e = {m->end.x,m->end.y,m->end.z,0,0,0,0,0,0};
        liscio_add_rapid(c, &s, &e, m->line_no);
        return;
    }
    if (m->type == NGC_MOVE_LINE) {
        liscio_pose_t s = {m->start.x,m->start.y,m->start.z,0,0,0,0,0,0};
        liscio_pose_t e = {m->end.x,m->end.y,m->end.z,0,0,0,0,0,0};
        liscio_add_line(c, &s, &e, m->feedrate, m->line_no);
    } else {
        liscio_primitive_t a; memset(&a,0,sizeof(a));
        a.type  = LISCIO_PRIM_ARC;
        a.start = (liscio_pose_t){m->start.x,m->start.y,m->start.z,0,0,0,0,0,0};
        a.end   = (liscio_pose_t){m->end.x,m->end.y,m->end.z,0,0,0,0,0,0};
        a.cx=m->cx; a.cy=m->cy; a.cz=m->cz;
        a.nx=m->nx; a.ny=m->ny; a.nz=m->nz;
        a.radius=m->radius; a.arc_angle=m->arc_angle;
        liscio_add_arc(c, &a, m->feedrate, m->line_no);
    }
}

/* Per-primitive verification: sample densely, then for each raw point
 * absorbed by this primitive, find nearest sample, accumulate max dist. */
static int    g_prim_idx_dbg = 0;
static long   g_drift_count = 0;
static double g_drift_max = 0;
static void prim_on_prim(const liscio_primitive_t *p, void *user)
{
    (void)user;
    /* RAPID / STOP are passthrough events, not fits — skip verification. */
    if (p->type == LISCIO_PRIM_RAPID || p->type == LISCIO_PRIM_STOP) return;
    int absorbed = p->n_absorbed;
    if (absorbed <= 0) absorbed = 1;

    /* Sample primitive at 200 evenly-spaced points. */
    const int N = 200;
    static double sx[201], sy[201], sz[201];
    for (int i = 0; i <= N; i++) {
        double t = (double)i / N;
        prim_eval(p, t, &sx[i], &sy[i], &sz[i]);
    }

    /* Raw range absorbed by THIS primitive: [cursor, cursor+absorbed]. */
    int i0 = (int)g_abs_cursor;
    int i1 = i0 + absorbed;
    if (i1 > g_raw_n - 1) i1 = g_raw_n - 1;

    /* Chord-length sanity (BEZIER/SPLINE only).  Skip:
     *  - ARC: ratio = arc_length/chord > 1 by construction (π/2 for half).
     *  - LINE: fit_len = chord(start, end) by construction, but raw[i0]
     *    may differ from prim.start when liscio flushes on small G0
     *    retracts (drilling peck cycle).  Then "fit_len > raw_len" is
     *    just the rapid offset and not a real LSQ drift. */
    if ((p->type == LISCIO_PRIM_BEZIER || p->type == LISCIO_PRIM_SPLINE)
        && absorbed >= 1 && i0 >= 0 && i1 < g_raw_n) {
        double raw_len = 0;
        for (int k = i0; k < i1; k++) {
            double dx = g_raw[k+1].x - g_raw[k].x;
            double dy = g_raw[k+1].y - g_raw[k].y;
            double dz = g_raw[k+1].z - g_raw[k].z;
            raw_len += sqrt(dx*dx + dy*dy + dz*dz);
        }
        double fit_len = 0;
        for (int i = 1; i <= N; i++) {
            double dx = sx[i]-sx[i-1], dy = sy[i]-sy[i-1], dz = sz[i]-sz[i-1];
            fit_len += sqrt(dx*dx + dy*dy + dz*dz);
        }
        double ratio = raw_len > 1e-9 ? fit_len / raw_len : 1.0;
        if (ratio > 1.5) {
            g_drift_count++;
            if (ratio > g_drift_max) g_drift_max = ratio;
            if (getenv("LISCIO_VDBG")) {
                const char *tn = (p->type==LISCIO_PRIM_LINE)?"LINE":(p->type==LISCIO_PRIM_ARC)?"ARC":(p->type==LISCIO_PRIM_BEZIER)?"BEZ":"SPL";
                fprintf(stderr, "  DRIFT prim[%d] %s n_abs=%d range=[%d..%d] raw_len=%.3f fit_len=%.3f ratio=%.2f\n",
                    g_prim_idx_dbg, tn, absorbed, i0, i1, raw_len, fit_len, ratio);
                if (g_drift_count <= 3) {
                    fprintf(stderr, "    p.start=(%.3f, %.3f, %.3f) p.end=(%.3f, %.3f, %.3f)\n",
                        p->start.x, p->start.y, p->start.z,
                        p->end.x, p->end.y, p->end.z);
                    fprintf(stderr, "    raw waypoints (incl 2 around):\n");
                    int k0 = i0 > 1 ? i0 - 1 : i0;
                    int k1 = i1 + 1 < g_raw_n ? i1 + 1 : i1;
                    for (int k = k0; k <= k1; k++)
                        fprintf(stderr, "      [%d] (%.3f, %.3f, %.3f)\n",
                            k, g_raw[k].x, g_raw[k].y, g_raw[k].z);
                }
            }
        }
    }

    double this_max = 0;
    int    this_worst = i0;
    /* Skip the shared start endpoint — it belongs to the previous primitive
     * geometrically, and after G0 rapids or discontinuities the "raw[i0]"
     * may not equal this primitive's start.  End point (raw[i1]) and all
     * interior absorbed waypoints must lie on the fit curve. */
    for (int k = i0 + 1; k <= i1; k++) {
        if (k < 0 || k >= g_raw_n) continue;
        double rx = g_raw[k].x, ry = g_raw[k].y, rz = g_raw[k].z;
        /* Min distance to the sampled polyline (segment, not vertex). */
        double best = 1e30;
        for (int i = 0; i < N; i++) {
            double ax = sx[i],   ay = sy[i],   az = sz[i];
            double bx = sx[i+1], by = sy[i+1], bz = sz[i+1];
            double ab_x = bx-ax, ab_y = by-ay, ab_z = bz-az;
            double ab2 = ab_x*ab_x + ab_y*ab_y + ab_z*ab_z;
            double tt;
            double cx_, cy_, cz_;
            if (ab2 < 1e-30) { cx_ = ax; cy_ = ay; cz_ = az; }
            else {
                tt = ((rx-ax)*ab_x + (ry-ay)*ab_y + (rz-az)*ab_z) / ab2;
                if (tt < 0) tt = 0; else if (tt > 1) tt = 1;
                cx_ = ax + tt*ab_x; cy_ = ay + tt*ab_y; cz_ = az + tt*ab_z;
            }
            double dx = cx_-rx, dy = cy_-ry, dz = cz_-rz;
            double d2 = dx*dx + dy*dy + dz*dz;
            if (d2 < best) best = d2;
        }
        double d = sqrt(best);
        if (d > g_obs_max) g_obs_max = d;
        if (d > this_max) { this_max = d; this_worst = k; }
        g_obs_sum += d;
        g_obs_count++;
        if (d > g_tol) g_over_tol++;
    }

    if (getenv("LISCIO_VDBG2")) {
        const char *tn2 = (p->type==LISCIO_PRIM_LINE)?"LINE":(p->type==LISCIO_PRIM_ARC)?"ARC":(p->type==LISCIO_PRIM_BEZIER)?"BEZ":"SPL";
        fprintf(stderr, "  prim[%d] %s n_abs=%d cur=%ld→%ld start=(%.3f,%.3f,%.3f) end=(%.3f,%.3f,%.3f)\n",
            g_prim_idx_dbg, tn2, absorbed, g_abs_cursor, g_abs_cursor+absorbed,
            p->start.x, p->start.y, p->start.z, p->end.x, p->end.y, p->end.z);
    }
    if (getenv("LISCIO_VDBG") && this_max > g_tol) {
        const char *tn = "?";
        if (p->type==LISCIO_PRIM_LINE) tn="LINE";
        else if (p->type==LISCIO_PRIM_ARC) tn="ARC";
        else if (p->type==LISCIO_PRIM_BEZIER) tn="BEZ";
        else if (p->type==LISCIO_PRIM_SPLINE) tn="SPL";
        fprintf(stderr, "  prim[%d] %-4s n_abs=%d range=[%d..%d] worst_pt=%d max=%.4f\n",
            g_prim_idx_dbg, tn, absorbed, i0, i1, this_worst, this_max);
        fprintf(stderr, "    p.start=(%.4f,%.4f,%.4f) p.end=(%.4f,%.4f,%.4f)\n",
            p->start.x, p->start.y, p->start.z, p->end.x, p->end.y, p->end.z);
        fprintf(stderr, "    raw[%d]=(%.4f,%.4f,%.4f)\n",
            this_worst, g_raw[this_worst].x, g_raw[this_worst].y, g_raw[this_worst].z);
        fprintf(stderr, "    raw[i0]=(%.4f,%.4f,%.4f) raw[i1]=(%.4f,%.4f,%.4f)\n",
            g_raw[i0].x, g_raw[i0].y, g_raw[i0].z, g_raw[i1].x, g_raw[i1].y, g_raw[i1].z);
        if (p->type==LISCIO_PRIM_BEZIER) {
            fprintf(stderr, "    P1=(%.4f,%.4f,%.4f) P2=(%.4f,%.4f,%.4f) w=(%.3f,%.3f,%.3f,%.3f) liscio_dev=%.4f\n",
                p->P1.x, p->P1.y, p->P1.z, p->P2.x, p->P2.y, p->P2.z,
                p->w0, p->w1, p->w2, p->w3, p->bezier_max_deviation);
            /* Find t that minimizes distance to worst_pt. */
            double best_t = 0, best_d2 = 1e30;
            for (int i = 0; i <= 1000; i++) {
                double tt = (double)i/1000;
                double ex, ey, ez;
                liscio_rbezier_eval(p, tt, &ex, &ey, &ez);
                double dx = ex - g_raw[this_worst].x;
                double dy = ey - g_raw[this_worst].y;
                double dz = ez - g_raw[this_worst].z;
                double d2 = dx*dx+dy*dy+dz*dz;
                if (d2 < best_d2) { best_d2 = d2; best_t = tt; }
            }
            fprintf(stderr, "    fine-grained min dist=%.4f at t=%.3f\n",
                sqrt(best_d2), best_t);
        }
    }
    g_prim_idx_dbg++;
    g_abs_cursor += absorbed;
}

static int run_file(const char *path, const liscio_cfg_t *cfg)
{
    /* Reset state. */
    g_raw_n = 0; g_abs_cursor = 0;
    g_obs_max = 0; g_obs_sum = 0; g_obs_count = 0; g_over_tol = 0;
    g_prim_idx_dbg = 0;
    g_drift_count = 0; g_drift_max = 0;

    LISCIO_TEST_TIMER_START(t0);

    /* Pass 1: parse NGC → fill g_raw[] with raw waypoints. */
    ngc_stats_t st;
    if (ngc_parse_file(path, raw_on_move, NULL, &st) != 0) return -1;

    /* Pass 2: feed same stream through liscio; on emit, compare against raw. */
    liscio_ctx_t *ctx = liscio_create(cfg);
    liscio_set_callback(ctx, prim_on_prim, NULL);

    ngc_parse_file(path, feed_liscio, ctx, &st);
    liscio_flush(ctx);
    liscio_stats_t ls;
    liscio_get_stats(ctx, &ls);
    liscio_destroy(ctx);
    long liscio_absorbed = ls.absorbed_total;
    long liscio_input = ls.input_count;
    (void)liscio_absorbed; (void)liscio_input;

    double elapsed_ms = LISCIO_TEST_TIMER_MS(t0);
    const char *name = strrchr(path, '/');
    name = name ? name+1 : path;
    double avg = g_obs_count ? g_obs_sum / g_obs_count : 0;
    long cursor_final = g_abs_cursor;
    printf("%-30s  raw=%-6d  abs=%-6ld  max_obs=%-7.4f  OVER_TOL=%-4ld  "
           "DRIFT=%-3ld (max %.2fx)  %.1f ms\n",
        name, g_raw_n, liscio_absorbed, g_obs_max, g_over_tol,
        g_drift_count, g_drift_max, elapsed_ms);
    (void)cursor_final; (void)g_obs_count; (void)avg; (void)liscio_input;
    return 0;
}

int main(int argc, char **argv)
{
    liscio_cfg_t cfg; liscio_cfg_default(&cfg);
    g_tol = 0.05;
    g_samples = 8;
    int i0 = 1;
    for (int i = 1; i < argc; i++) {
        if (strncmp(argv[i], "--tol_xyz=", 10) == 0) {
            cfg.tol_xyz = atof(argv[i]+10); g_tol = cfg.tol_xyz; i0 = i+1;
        } else if (strncmp(argv[i], "--samples=", 10) == 0) {
            g_samples = atoi(argv[i]+10); i0 = i+1;
        } else if (strncmp(argv[i], "--arc_to_helix=", 15) == 0) {
            cfg.arc_to_helix_max = atoi(argv[i]+15); i0 = i+1;
        } else if (strncmp(argv[i], "--arc_subseg=", 13) == 0) {
            cfg.arc_subseg_samples = atoi(argv[i]+13); i0 = i+1;
        } else break;
    }
    if (i0 >= argc) {
        fprintf(stderr, "usage: %s [--tol_xyz=N] [--samples=N] <ngc>...\n", argv[0]);
        return 1;
    }

    printf("# liscio independent fit-quality verifier\n");
    printf("# tol_xyz=%.4f (OVER_TOL counts waypoints > tol)\n#\n",
        cfg.tol_xyz);
    printf("%-30s  %-10s %-13s %-14s %-14s %s\n",
        "file", "raw", "checked", "max_obs", "avg", "status");

    LISCIO_TEST_TIMER_START(__g);
    int fails = 0;
    long total_over_tol = 0, total_drift = 0;
    for (int i = i0; i < argc; i++) {
        if (run_file(argv[i], &cfg) != 0) fails++;
        total_over_tol += g_over_tol;
        total_drift    += g_drift_count;
    }
    printf("---- liscio verify: %d file%s, total OVER_TOL=%ld DRIFT=%ld, %.3f ms total ----\n",
        argc - i0, (argc - i0 == 1 ? "" : "s"),
        total_over_tol, total_drift, LISCIO_TEST_TIMER_MS(__g));
    /* Exit non-zero if any waypoint exceeds tolerance OR any primitive
     * drifts > 1.5× its absorbed chord — gates regressions in CTest. */
    if (total_over_tol > 0 || total_drift > 0) return 2;
    return fails ? 1 : 0;
}
