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
* Description: Quantitative fit-quality report for CAM NGCs.
*              Measures per file:
*                - compression ratio (waypoints / primitives)
*                - max fit deviation (primitive sample vs raw)
*                - path-length ratio (fit vs raw polyline)
*                - tangent mismatch at primitive boundaries (max, mean)
*                - liscio pass time (µs/waypoint)
*              Produces stdout table + CSV report suitable for
*              comparison across config changes.
*              Usage: test_metrics <ngc>... [--tol_xyz=N] [--tol_abc=N]
* Author:      杨阳 (Yang Yang) <mika-net@outlook.com>
* License:     MIT (SPDX-License-Identifier: MIT)
* Copyright (c) 2026 杨阳 (Yang Yang)
********************************************************************/

#define _USE_MATH_DEFINES
#define _POSIX_C_SOURCE 199309L
#include "liscio/liscio.h"
#include "test_timing.h"
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

/* ---------- accumulated metrics ---------- */
typedef struct {
    long   n_input_g01;
    long   n_input_g23;
    long   n_emit_line, n_emit_arc, n_emit_bezier, n_emit_spline, n_emit_helix;
    double raw_length;
    double fit_length;
    double max_dev;          /* max reported primitive max_deviation */
    double tan_max_deg;      /* max tangent mismatch at primitive boundaries */
    double tan_sum_deg;      /* sum for mean */
    long   tan_count;
    /* Per-bucket counts (for diagnosis): smooth vs corner-induced tan
     * mismatch at primitive boundaries.  Buckets: <5°, 5-30°, 30-90°,
     * 90-150°, 150-180°. */
    long   tan_bucket[5];
    long   n_emit_rapid;     /* # LISCIO_PRIM_RAPID passthrough */
    long   n_emit_stop;      /* # LISCIO_PRIM_STOP events */
    double prev_end_tan_x, prev_end_tan_y, prev_end_tan_z;
    int    have_prev_tan;
    double t_elapsed;
    /* debug breakdown */
    double fit_len_line, fit_len_arc, fit_len_bez, fit_len_spl;
    int    dbg_dumped;
} metrics_t;

static metrics_t g_m;

static double rad2deg(double r) { return r * 180.0 / M_PI; }

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* Sample a primitive into N points to measure fit length + deviation. */
static void sample_prim(const liscio_primitive_t *p, int N,
                        double *xs, double *ys, double *zs)
{
    for (int i = 0; i <= N; i++) {
        double t = (double)i / N;
        switch (p->type) {
            case LISCIO_PRIM_LINE:
                xs[i] = p->start.x + t*(p->end.x - p->start.x);
                ys[i] = p->start.y + t*(p->end.y - p->start.y);
                zs[i] = p->start.z + t*(p->end.z - p->start.z);
                break;
            case LISCIO_PRIM_HELIX:    /* arc fields + axial pitch term */
            case LISCIO_PRIM_ARC: {
                double vx = p->start.x - p->cx;
                double vy = p->start.y - p->cy;
                double vz = p->start.z - p->cz;
                double m = sqrt(vx*vx + vy*vy + vz*vz);
                if (m < 1e-15) { xs[i]=p->start.x; ys[i]=p->start.y; zs[i]=p->start.z; break; }
                double ux = vx/m, uy = vy/m, uz = vz/m;
                double rx = p->ny*uz - p->nz*uy;
                double ry = p->nz*ux - p->nx*uz;
                double rz = p->nx*uy - p->ny*ux;
                double ang = t * p->arc_angle;
                double axial = ang * p->pitch;     /* HELIX advance, 0 for ARC */
                xs[i] = p->cx + p->radius * (cos(ang)*ux + sin(ang)*rx) + axial * p->nx;
                ys[i] = p->cy + p->radius * (cos(ang)*uy + sin(ang)*ry) + axial * p->ny;
                zs[i] = p->cz + p->radius * (cos(ang)*uz + sin(ang)*rz) + axial * p->nz;
                /* Pure planar ARC basis with z-mismatch → linear Z. */
                if (p->type == LISCIO_PRIM_ARC && fabs(uz) + fabs(rz) < 1e-12)
                    zs[i] = p->start.z + t * (p->end.z - p->start.z);
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
}

static void on_prim(const liscio_primitive_t *p, void *user)
{
    (void)user;
    switch (p->type) {
        case LISCIO_PRIM_LINE:   g_m.n_emit_line++; break;
        case LISCIO_PRIM_HELIX:  g_m.n_emit_helix++; g_m.n_emit_arc++; break;
        case LISCIO_PRIM_ARC:    g_m.n_emit_arc++; break;
        case LISCIO_PRIM_BEZIER: g_m.n_emit_bezier++; break;
        case LISCIO_PRIM_SPLINE: g_m.n_emit_spline++; break;
        case LISCIO_PRIM_RAPID:
            g_m.n_emit_rapid++;
            g_m.have_prev_tan = 0;   /* break tangent chain across G0 */
            return;
        case LISCIO_PRIM_STOP:
            g_m.n_emit_stop++;
            g_m.have_prev_tan = 0;   /* break across stop event */
            return;
    }

    /* fit length via 40-sample polyline */
    const int N = 40;
    double xs[64], ys[64], zs[64];
    sample_prim(p, N, xs, ys, zs);
    double plen = 0;
    for (int i = 1; i <= N; i++) {
        double dx=xs[i]-xs[i-1], dy=ys[i]-ys[i-1], dz=zs[i]-zs[i-1];
        plen += sqrt(dx*dx + dy*dy + dz*dz);
    }
    g_m.fit_length += plen;
    switch (p->type) {
        case LISCIO_PRIM_LINE:   g_m.fit_len_line += plen; break;
        case LISCIO_PRIM_ARC:    g_m.fit_len_arc  += plen; break;
        case LISCIO_PRIM_BEZIER: g_m.fit_len_bez  += plen; break;
        case LISCIO_PRIM_SPLINE: g_m.fit_len_spl  += plen; break;
    }
    if (getenv("LISCIO_DBG") && g_m.dbg_dumped < 30) {
        const char *tn = "?";
        if (p->type==LISCIO_PRIM_LINE) tn="LINE";
        else if (p->type==LISCIO_PRIM_ARC) tn="ARC";
        else if (p->type==LISCIO_PRIM_BEZIER) tn="BEZ";
        else if (p->type==LISCIO_PRIM_SPLINE) tn="SPL";
        fprintf(stderr, "  prim[%d] %-4s n_abs=%d plen=%.6f  s=(%.4f,%.4f,%.4f) e=(%.4f,%.4f,%.4f)",
            g_m.dbg_dumped, tn, p->n_absorbed, plen,
            p->start.x, p->start.y, p->start.z,
            p->end.x, p->end.y, p->end.z);
        if (p->type==LISCIO_PRIM_ARC)
            fprintf(stderr, "  c=(%.4f,%.4f,%.4f) r=%.4f ang=%.4f",
                p->cx, p->cy, p->cz, p->radius, p->arc_angle);
        fprintf(stderr, "\n");
        g_m.dbg_dumped++;
    }

    /* track max fit deviation */
    double d = 0;
    if (p->type == LISCIO_PRIM_ARC)    d = 0.0; /* arc exact within fitted pts */
    if (p->type == LISCIO_PRIM_BEZIER) d = p->bezier_max_deviation;
    if (p->type == LISCIO_PRIM_SPLINE) d = p->spline_max_deviation;
    if (d > g_m.max_dev) g_m.max_dev = d;

    /* (RAPID / STOP already returned above; tangent accounting only
     * runs for fitted primitives.  RAPID/STOP also break the cross-
     * primitive tangent chain — handled by the early-returns since
     * have_prev_tan stays false until next fit emit.) */

    /* tangent-continuity check at start of primitive */
    if (g_m.have_prev_tan) {
        double dot = g_m.prev_end_tan_x * p->tan_start_x
                   + g_m.prev_end_tan_y * p->tan_start_y
                   + g_m.prev_end_tan_z * p->tan_start_z;
        if (dot > 1.0) dot = 1.0;
        if (dot < -1.0) dot = -1.0;
        double ang = rad2deg(acos(dot));
        if (ang > g_m.tan_max_deg) g_m.tan_max_deg = ang;
        g_m.tan_sum_deg += ang;
        g_m.tan_count++;
        if      (ang <   5.0) g_m.tan_bucket[0]++;
        else if (ang <  30.0) g_m.tan_bucket[1]++;
        else if (ang <  90.0) g_m.tan_bucket[2]++;
        else if (ang < 150.0) g_m.tan_bucket[3]++;
        else                  g_m.tan_bucket[4]++;
    }
    g_m.prev_end_tan_x = p->tan_end_x;
    g_m.prev_end_tan_y = p->tan_end_y;
    g_m.prev_end_tan_z = p->tan_end_z;
    g_m.have_prev_tan = 1;
}

/* ---------- NGC parser bridge (uses shared ngc_parser) ---------- */
typedef struct {
    liscio_ctx_t *ctx;
    double        feed;
} metrics_ctx_t;

static void on_move(const ngc_move_t *m, void *user)
{
    metrics_ctx_t *mc = (metrics_ctx_t *)user;
    if (m->is_rapid) {
        /* G0 → emit RAPID passthrough event */
        liscio_pose_t s = {m->start.x,m->start.y,m->start.z,
                           m->start.a,m->start.b,m->start.c,
                           m->start.u,m->start.v,m->start.w};
        liscio_pose_t e = {m->end.x,m->end.y,m->end.z,
                           m->end.a,m->end.b,m->end.c,
                           m->end.u,m->end.v,m->end.w};
        liscio_add_rapid(mc->ctx, &s, &e, m->line_no);
        return;
    }

    if (m->type == NGC_MOVE_LINE) {
        g_m.n_input_g01++;
        double dx = m->end.x - m->start.x;
        double dy = m->end.y - m->start.y;
        double dz = m->end.z - m->start.z;
        g_m.raw_length += sqrt(dx*dx + dy*dy + dz*dz);
        liscio_pose_t s = {m->start.x, m->start.y, m->start.z,
                            m->start.a, m->start.b, m->start.c,
                            m->start.u, m->start.v, m->start.w};
        liscio_pose_t e = {m->end.x, m->end.y, m->end.z,
                            m->end.a, m->end.b, m->end.c,
                            m->end.u, m->end.v, m->end.w};
        liscio_add_line(mc->ctx, &s, &e, m->feedrate, m->line_no);
    } else {
        g_m.n_input_g23++;
        liscio_primitive_t a; memset(&a, 0, sizeof(a));
        a.type  = LISCIO_PRIM_ARC;
        a.start = (liscio_pose_t){m->start.x, m->start.y, m->start.z,
                                  m->start.a, m->start.b, m->start.c,
                                  m->start.u, m->start.v, m->start.w};
        a.end   = (liscio_pose_t){m->end.x, m->end.y, m->end.z,
                                  m->end.a, m->end.b, m->end.c,
                                  m->end.u, m->end.v, m->end.w};
        a.cx = m->cx; a.cy = m->cy; a.cz = m->cz;
        a.nx = m->nx; a.ny = m->ny; a.nz = m->nz;
        a.radius = m->radius;
        a.arc_angle = m->arc_angle;
        /* raw arc length including helical off-plane component. */
        double planar = fabs(m->arc_angle) * m->radius;
        double dzh, dzo;
        if (m->plane == 0)      { dzh = m->end.z - m->start.z; dzo = 0; }
        else if (m->plane == 1) { dzh = m->end.y - m->start.y; dzo = 0; }
        else                    { dzh = m->end.x - m->start.x; dzo = 0; }
        (void)dzo;
        g_m.raw_length += sqrt(planar*planar + dzh*dzh);
        liscio_add_arc(mc->ctx, &a, m->feedrate, m->line_no);
    }
}

static int run_file(const char *path, const liscio_cfg_t *cfg, FILE *csv)
{
    memset(&g_m, 0, sizeof(g_m));
    liscio_ctx_t *ctx = liscio_create(cfg);
    liscio_set_callback(ctx, on_prim, NULL);

    metrics_ctx_t mc = { ctx, 1000.0 };
    double t0 = now_sec();
    ngc_stats_t st;
    int rc = ngc_parse_file(path, on_move, &mc, &st);
    if (rc != 0) { liscio_destroy(ctx); return -1; }
    liscio_flush(ctx);
    g_m.t_elapsed = now_sec() - t0;
    liscio_destroy(ctx);

    long n_in = g_m.n_input_g01 + g_m.n_input_g23;
    long n_em = g_m.n_emit_line + g_m.n_emit_arc + g_m.n_emit_bezier + g_m.n_emit_spline;
    double ratio  = n_em ? (double)n_in / n_em : 0.0;
    double len_rel = g_m.raw_length > 1e-12
                    ? (g_m.fit_length - g_m.raw_length) / g_m.raw_length
                    : 0.0;
    double tan_mean = g_m.tan_count ? g_m.tan_sum_deg / g_m.tan_count : 0.0;
    double us_per = n_in ? g_m.t_elapsed * 1e6 / n_in : 0.0;

    const char *name = strrchr(path, '/');
    name = name ? name + 1 : path;

    printf("%-30s  in=%-6ld  out=%-5ld  ratio=%6.2fx  "
           "max_dev=%-7.4f  len_err=%+7.3f%%  tan_max=%-5.1f°  tan_mean=%-4.1f°  "
           "rapid=%-3ld stop=%-2ld helix=%-3ld  %.1f µs/wp\n",
        name, n_in, n_em, ratio,
        g_m.max_dev, len_rel * 100.0,
        g_m.tan_max_deg, tan_mean,
        g_m.n_emit_rapid, g_m.n_emit_stop, g_m.n_emit_helix, us_per);
    if (getenv("LISCIO_DBG")) {
        fprintf(stderr, "  raw=%.4f fit=%.4f (line=%.4f arc=%.4f bez=%.4f spl=%.4f)\n",
            g_m.raw_length, g_m.fit_length,
            g_m.fit_len_line, g_m.fit_len_arc, g_m.fit_len_bez, g_m.fit_len_spl);
    }
#ifdef LISCIO_DEBUG_G1
    if (getenv("LISCIO_G1_STATS")) {
        extern long liscio_dbg_g1_ok, liscio_dbg_g1_fall;
        extern long liscio_dbg_g1_fast_ok, liscio_dbg_g1_fast_fall;
        long c_total = liscio_dbg_g1_ok + liscio_dbg_g1_fall;
        long f_total = liscio_dbg_g1_fast_ok + liscio_dbg_g1_fast_fall;
        fprintf(stderr,
            "  g1-stats: composite OK=%ld FALL=%ld (%.1f%% G1)  "
            "fast-path OK=%ld FALL=%ld (%.1f%% G1)\n",
            liscio_dbg_g1_ok, liscio_dbg_g1_fall,
            c_total ? 100.0 * liscio_dbg_g1_ok / c_total : 0.0,
            liscio_dbg_g1_fast_ok, liscio_dbg_g1_fast_fall,
            f_total ? 100.0 * liscio_dbg_g1_fast_ok / f_total : 0.0);
        liscio_dbg_g1_ok = liscio_dbg_g1_fall = 0;
        liscio_dbg_g1_fast_ok = liscio_dbg_g1_fast_fall = 0;
    }
#endif
    if (getenv("LISCIO_TYPE_HIST")) {
        long e = g_m.n_emit_line + g_m.n_emit_arc + g_m.n_emit_bezier + g_m.n_emit_spline;
        if (e > 0) {
            fprintf(stderr,
                "  type-hist: LINE=%ld(%.0f%%) ARC=%ld(%.0f%%) "
                "BEZ=%ld(%.0f%%) SPL=%ld(%.0f%%)\n",
                g_m.n_emit_line,   100.0*g_m.n_emit_line/e,
                g_m.n_emit_arc,    100.0*g_m.n_emit_arc/e,
                g_m.n_emit_bezier, 100.0*g_m.n_emit_bezier/e,
                g_m.n_emit_spline, 100.0*g_m.n_emit_spline/e);
        }
    }
    if (getenv("LISCIO_TAN_HIST") && g_m.tan_count > 0) {
        long c = g_m.tan_count;
        fprintf(stderr,
            "  tan-hist: <5°=%ld(%.0f%%) 5-30°=%ld(%.0f%%) 30-90°=%ld(%.0f%%) "
            "90-150°=%ld(%.0f%%) >150°=%ld(%.0f%%)\n",
            g_m.tan_bucket[0], 100.0*g_m.tan_bucket[0]/c,
            g_m.tan_bucket[1], 100.0*g_m.tan_bucket[1]/c,
            g_m.tan_bucket[2], 100.0*g_m.tan_bucket[2]/c,
            g_m.tan_bucket[3], 100.0*g_m.tan_bucket[3]/c,
            g_m.tan_bucket[4], 100.0*g_m.tan_bucket[4]/c);
    }

    if (csv) {
        /* Per-type breakdown: ARC counter already includes HELIX (tracked
         * in n_emit_helix), so n_arc_only = n_arc - n_helix. */
        long n_arc_only = g_m.n_emit_arc - g_m.n_emit_helix;
        fprintf(csv,
            "%s,%ld,%ld,%.6f,%.6f,%.6f,%.2f,%.4f,%.6f,%.2f,%.2f,%.3f,"
            "%ld,%ld,%ld,%ld,%ld,%ld,%ld\n",
            name, n_in, n_em, ratio, g_m.raw_length, g_m.fit_length,
            len_rel * 100.0, g_m.max_dev, g_m.t_elapsed * 1e6 / (n_in > 0 ? n_in : 1),
            g_m.tan_max_deg, tan_mean, us_per,
            g_m.n_emit_line, n_arc_only, g_m.n_emit_helix,
            g_m.n_emit_bezier, g_m.n_emit_spline,
            g_m.n_emit_rapid, g_m.n_emit_stop);
    }
    return 0;
}

int main(int argc, char **argv)
{
    liscio_cfg_t cfg; liscio_cfg_default(&cfg);
    const char *csv_path = NULL;
    int i0 = 1;
    for (int i = 1; i < argc; i++) {
        if (strncmp(argv[i], "--tol_xyz=", 10) == 0) {
            cfg.tol_xyz = atof(argv[i] + 10); i0 = i + 1;
        } else if (strncmp(argv[i], "--tol_abc=", 10) == 0) {
            cfg.tol_abc = atof(argv[i] + 10); i0 = i + 1;
        } else if (strncmp(argv[i], "--csv=", 6) == 0) {
            csv_path = argv[i] + 6; i0 = i + 1;
        } else if (strncmp(argv[i], "--arc_subseg=", 13) == 0) {
            cfg.arc_subseg_samples = atoi(argv[i] + 13); i0 = i + 1;
        } else if (strncmp(argv[i], "--arc_merge=", 12) == 0) {
            int m = atoi(argv[i] + 12);
            liscio_cfg_set_arc_merge(&cfg,
                (m == 2) ? LISCIO_ARC_MERGE_HELIX :
                (m == 1) ? LISCIO_ARC_MERGE_ARC :
                           LISCIO_ARC_MERGE_OFF, 0, 0);
            i0 = i + 1;
        } else if (strncmp(argv[i], "--arc_merge_xy_tol=", 19) == 0) {
            cfg.arc_merge_xy_tol = atof(argv[i] + 19); i0 = i + 1;
        } else if (strncmp(argv[i], "--arc_merge_pitch_tol=", 22) == 0) {
            cfg.arc_merge_pitch_tol = atof(argv[i] + 22); i0 = i + 1;
        } else if (strncmp(argv[i], "--arc_to_helix=", 15) == 0) {
            cfg.arc_to_helix_max = atoi(argv[i] + 15); i0 = i + 1;
        } else if (strncmp(argv[i], "--corner_mode=", 14) == 0) {
            int m = atoi(argv[i] + 14);
            liscio_cfg_set_corner_detection(&cfg,
                m ? LISCIO_CORNER_LOOKAHEAD : LISCIO_CORNER_IMMEDIATE,
                0, 0);
            i0 = i + 1;
        } else if (strncmp(argv[i], "--corner_soft=", 14) == 0) {
            cfg.corner_angle_deg = atof(argv[i] + 14); i0 = i + 1;
        } else if (strncmp(argv[i], "--corner_hard=", 14) == 0) {
            cfg.corner_hard_angle_deg = atof(argv[i] + 14); i0 = i + 1;
        } else break;
    }
    if (i0 >= argc) {
        fprintf(stderr,
            "usage: %s [--tol_xyz=N] [--tol_abc=N] [--arc_subseg=N] "
            "[--arc_merge=0|1|2] [--csv=out.csv] <ngc>...\n",
            argv[0]);
        return 1;
    }

    printf("# liscio fit-quality metrics (tol_xyz=%.4g tol_abc=%.4f "
           "arc_subseg=%d arc_merge=%d)\n",
        cfg.tol_xyz, cfg.tol_abc,
        cfg.arc_subseg_samples, cfg.arc_merge_mode);
    printf("#\n");
    printf("# in        = input motion blocks (G1+G2+G3)\n");
    printf("# out       = emitted primitives\n");
    printf("# ratio     = compression (in/out), higher = better\n");
    printf("# max_dev   = worst reported primitive fit deviation (mm)\n");
    printf("# len_err   = (fit_length - raw_length) / raw_length  (should be ≈0)\n");
    printf("# tan_max   = max tangent-mismatch angle at primitive boundaries (deg)\n");
    printf("# tan_mean  = mean tangent-mismatch\n");
    printf("# µs/wp     = CPU time per input waypoint\n#\n");

    FILE *csv = NULL;
    if (csv_path) {
        csv = fopen(csv_path, "w");
        if (csv) fprintf(csv,
            "file,input,output,ratio,raw_len,fit_len,len_err_pct,"
            "max_dev,elapsed_us,tan_max_deg,tan_mean_deg,us_per_wp,"
            "n_line,n_arc,n_helix,n_bezier,n_spline,n_rapid,n_stop\n");
    }

    LISCIO_TEST_TIMER_START(__t0);
    int fails = 0;
    for (int i = i0; i < argc; i++) {
        if (run_file(argv[i], &cfg, csv) != 0) fails++;
    }
    if (csv) fclose(csv);
    printf("---- liscio metrics: %d file%s, %.3f ms total ----\n",
        argc - i0, (argc - i0 == 1 ? "" : "s"), LISCIO_TEST_TIMER_MS(__t0));
    return fails ? 1 : 0;
}
