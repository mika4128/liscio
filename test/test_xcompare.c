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
* Description: Cross-library comparison harness.  Same input polyline
*              (extracted from CAM NGCs by rs274) is fed to each of:
*                baseline   — passthrough (segs == waypoints-1)
*                dp-xyz     — Douglas-Peucker 3D simplification (1973)
*                schneider  — Schneider FitCurves (Graphics Gems I 1990,
*                             3D ported, recursive cubic Bezier LSQ)
*                sisl       — SINTEF SISL s1961 LSQ B-spline approx
*                geomdl     — Piegl & Tiller §9.4 LSQ B-spline (Python)
*                liscio     — full pipeline (LINE+ARC+Bezier+B-spline)
*              Reports per file: # primitives, max_dev (XYZ Euclidean),
*              elapsed time.  Final Markdown table aggregates totals.
* Usage: test_xcompare [--tol_xyz=N] [--no-geomdl] [--no-sisl]
*                      [--no-schneider] <ngc>...
* Author:      杨阳 (Yang Yang) <mika-net@outlook.com>
* License:     MIT (SPDX-License-Identifier: MIT)
* Copyright (c) 2026 杨阳 (Yang Yang)
********************************************************************/

#define _USE_MATH_DEFINES
#define _POSIX_C_SOURCE 199309L

#include "liscio/liscio.h"
#include "ngc_parser.h"
#include "third_party/schneider3d.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

#ifdef HAVE_SISL
#include <sisl.h>
#endif

#define MAX_PTS 200000
#define ARC_SAMPLES 8

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

/* ============================================================ */
/* NGC → polyline bridge (arcs sampled, rapids excluded)        */
/* ============================================================ */

typedef struct { polyline_t *pl; int arc_samples; } pl_bridge_t;

static void pl_on_move(const ngc_move_t *m, void *user)
{
    pl_bridge_t *b = (pl_bridge_t *)user;
    polyline_t *pl = b->pl;
    if (m->is_rapid) return;
    if (pl->n == 0 && pl->n < MAX_PTS)
        pl->pts[pl->n++] = (vec3_t){m->start.x, m->start.y, m->start.z};
    if (m->type == NGC_MOVE_LINE) {
        double dx = m->end.x - pl->pts[pl->n-1].x;
        double dy = m->end.y - pl->pts[pl->n-1].y;
        double dz = m->end.z - pl->pts[pl->n-1].z;
        pl->raw_length += sqrt(dx*dx+dy*dy+dz*dz);
        if (pl->n < MAX_PTS) pl->pts[pl->n++] = (vec3_t){m->end.x, m->end.y, m->end.z};
        return;
    }
    for (int s = 1; s <= b->arc_samples; s++) {
        double t = (double)s / b->arc_samples;
        double ang = t * m->arc_angle;
        double vx = m->start.x - m->cx, vy = m->start.y - m->cy, vz = m->start.z - m->cz;
        double un = sqrt(vx*vx+vy*vy+vz*vz);
        if (un < 1e-15) continue;
        double ux=vx/un, uy=vy/un, uz=vz/un;
        double wx=m->ny*uz - m->nz*uy;
        double wy=m->nz*ux - m->nx*uz;
        double wz=m->nx*uy - m->ny*ux;
        double sx = m->cx + m->radius*(cos(ang)*ux + sin(ang)*wx);
        double sy = m->cy + m->radius*(cos(ang)*uy + sin(ang)*wy);
        double sz = m->cz + m->radius*(cos(ang)*uz + sin(ang)*wz);
        if (m->plane == 0 && fabs(uz)+fabs(wz) < 1e-12)
            sz = m->start.z + t * (m->end.z - m->start.z);
        else if (m->plane == 1 && fabs(uy)+fabs(wy) < 1e-12)
            sy = m->start.y + t * (m->end.y - m->start.y);
        else if (m->plane == 2 && fabs(ux)+fabs(wx) < 1e-12)
            sx = m->start.x + t * (m->end.x - m->start.x);
        if (pl->n < MAX_PTS) {
            double dx=sx-pl->pts[pl->n-1].x, dy=sy-pl->pts[pl->n-1].y, dz=sz-pl->pts[pl->n-1].z;
            pl->raw_length += sqrt(dx*dx+dy*dy+dz*dz);
            pl->pts[pl->n++] = (vec3_t){sx, sy, sz};
        }
    }
}

static int load_ngc(const char *path, polyline_t *pl)
{
    pl->n = 0; pl->raw_length = 0;
    pl_bridge_t br = { pl, ARC_SAMPLES };
    ngc_stats_t st;
    return ngc_parse_file(path, pl_on_move, &br, &st);
}

/* ============================================================ */
/* Result struct                                                */
/* ============================================================ */

typedef struct {
    int    n_emit;       /* output primitive / control-point count */
    double max_dev;      /* worst observed XYZ deviation (mm)      */
    double elapsed_sec;
    int    available;    /* 1 = ran, 0 = skipped (lib missing)     */
} strat_result_t;

/* ============================================================ */
/* Strategy: baseline                                           */
/* ============================================================ */

static void strat_baseline(const polyline_t *pl, strat_result_t *r)
{
    double t0 = now_sec();
    r->n_emit = pl->n - 1;
    r->max_dev = 0.0;
    r->elapsed_sec = now_sec() - t0;
    r->available = 1;
}

/* ============================================================ */
/* Strategy: Douglas-Peucker 3D                                 */
/* ============================================================ */

static double pt_seg_dist(vec3_t p, vec3_t a, vec3_t b)
{
    double abx=b.x-a.x, aby=b.y-a.y, abz=b.z-a.z;
    double ab2 = abx*abx+aby*aby+abz*abz;
    if (ab2 < 1e-30) {
        double dx=p.x-a.x, dy=p.y-a.y, dz=p.z-a.z;
        return sqrt(dx*dx+dy*dy+dz*dz);
    }
    double apx=p.x-a.x, apy=p.y-a.y, apz=p.z-a.z;
    double t = (apx*abx+apy*aby+apz*abz)/ab2;
    if (t<0) t=0; if (t>1) t=1;
    double cx=a.x+t*abx, cy=a.y+t*aby, cz=a.z+t*abz;
    double dx=p.x-cx, dy=p.y-cy, dz=p.z-cz;
    return sqrt(dx*dx+dy*dy+dz*dz);
}

static void dp_recurse(const vec3_t *pts, int i0, int i1, double tol, char *keep)
{
    if (i1 <= i0+1) return;
    double maxd=0; int maxi=i0;
    for (int k=i0+1; k<i1; k++) {
        double d = pt_seg_dist(pts[k], pts[i0], pts[i1]);
        if (d > maxd) { maxd=d; maxi=k; }
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
    int last_kept = 0;
    double max_d = 0;
    for (int i = 0; i < pl->n; i++) {
        if (!keep[i]) continue;
        if (i > 0) {
            for (int k = last_kept+1; k < i; k++) {
                double d = pt_seg_dist(pl->pts[k], pl->pts[last_kept], pl->pts[i]);
                if (d > max_d) max_d = d;
            }
        }
        last_kept = i;
        n_kept++;
    }
    r->n_emit = n_kept - 1;
    r->max_dev = max_d;
    r->elapsed_sec = now_sec() - t0;
    r->available = 1;
}

/* ============================================================ */
/* Strategy: Schneider FitCurves (3D port)                      */
/* ============================================================ */

typedef struct {
    int n_bez;
    schneider3d_bezier_t bez[MAX_PTS];   /* upper bound */
} schneider_collect_t;

static void sn_emit(const schneider3d_bezier_t *b, void *user)
{
    schneider_collect_t *c = (schneider_collect_t *)user;
    if (c->n_bez < MAX_PTS) c->bez[c->n_bez++] = *b;
}

static void schneider_eval(const schneider3d_bezier_t *b, double t,
                           double *x, double *y, double *z)
{
    double u = 1.0 - t;
    double b0 = u*u*u, b1 = 3*t*u*u, b2 = 3*t*t*u, b3 = t*t*t;
    *x = b0*b->P0[0] + b1*b->P1[0] + b2*b->P2[0] + b3*b->P3[0];
    *y = b0*b->P0[1] + b1*b->P1[1] + b2*b->P2[1] + b3*b->P3[1];
    *z = b0*b->P0[2] + b1*b->P1[2] + b2*b->P2[2] + b3*b->P3[2];
}

static void strat_schneider(const polyline_t *pl, double tol, strat_result_t *r)
{
    static double xyz_flat[MAX_PTS*3];
    for (int i = 0; i < pl->n; i++) {
        xyz_flat[3*i+0] = pl->pts[i].x;
        xyz_flat[3*i+1] = pl->pts[i].y;
        xyz_flat[3*i+2] = pl->pts[i].z;
    }
    static schneider_collect_t cc;
    cc.n_bez = 0;

    double t0 = now_sec();
    schneider3d_fit(xyz_flat, pl->n, tol*tol, sn_emit, &cc);
    double elapsed = now_sec() - t0;

    /* Measure actual max XYZ deviation by projecting input points
     * onto nearest sampled curve point (32 samples per Bezier). */
    double max_dev = 0.0;
    if (cc.n_bez > 0) {
        const int SAMP = 32;
        static double sx[MAX_PTS*32], sy[MAX_PTS*32], sz[MAX_PTS*32];
        int total = 0;
        for (int b = 0; b < cc.n_bez; b++) {
            for (int s = 0; s <= SAMP; s++) {
                double t = (double)s / SAMP;
                schneider_eval(&cc.bez[b], t, &sx[total], &sy[total], &sz[total]);
                total++;
                if (total >= MAX_PTS*32 - 1) break;
            }
            if (total >= MAX_PTS*32 - 1) break;
        }
        int hint = 0;
        for (int i = 0; i < pl->n; i++) {
            double best = 1e18;
            int lo = hint > 32 ? hint - 32 : 0;
            int hi = hint + 32 < total ? hint + 32 : total;
            int ji = hint;
            for (int j = lo; j < hi; j++) {
                double dx=pl->pts[i].x-sx[j], dy=pl->pts[i].y-sy[j], dz=pl->pts[i].z-sz[j];
                double d2 = dx*dx+dy*dy+dz*dz;
                if (d2 < best) { best = d2; ji = j; }
            }
            if (ji == lo || ji == hi-1) {
                /* Fallback global scan */
                for (int j = 0; j < total; j++) {
                    double dx=pl->pts[i].x-sx[j], dy=pl->pts[i].y-sy[j], dz=pl->pts[i].z-sz[j];
                    double d2 = dx*dx+dy*dy+dz*dz;
                    if (d2 < best) { best = d2; ji = j; }
                }
            }
            hint = ji;
            double d = sqrt(best);
            if (d > max_dev) max_dev = d;
        }
    }

    r->n_emit = cc.n_bez;
    r->max_dev = max_dev;
    r->elapsed_sec = elapsed;
    r->available = 1;
}

/* ============================================================ */
/* Strategy: SISL s1961 (LSQ B-spline)                          */
/* ============================================================ */

#ifdef HAVE_SISL
static void strat_sisl(const polyline_t *pl, double tol, strat_result_t *r)
{
    int n = pl->n;
    if (n < 4) { r->available = 0; return; }
    static double ep[MAX_PTS*3];
    for (int i = 0; i < n; i++) {
        ep[3*i+0] = pl->pts[i].x;
        ep[3*i+1] = pl->pts[i].y;
        ep[3*i+2] = pl->pts[i].z;
    }
    int    idim   = 3;
    int    ipar   = 1;       /* chord-length parameterization */
    double dummy_par[1] = {0};
    double eeps[3]= { tol, tol, tol };
    int    ilend  = 1;       /* keep position at left end */
    int    irend  = 1;       /* keep position at right end */
    int    iopen  = 1;       /* open curve */
    double afctol = 0.0;     /* default reduction */
    int    itmax  = 5;
    int    ik     = 4;       /* cubic */
    SISLCurve *rc = NULL;
    double emxerr[3] = {0,0,0};
    int    jstat = 0;

    double t0 = now_sec();
    s1961(ep, n, idim, ipar, dummy_par, eeps,
          ilend, irend, iopen, afctol, itmax, ik,
          &rc, emxerr, &jstat);
    double elapsed = now_sec() - t0;

    if (jstat < 0 || rc == NULL) {
        r->available = 1;
        r->n_emit = -1;
        r->max_dev = NAN;
        r->elapsed_sec = elapsed;
        return;
    }

    /* Per-dim max errors → Euclidean upper bound (sqrt of sum). */
    double dev3 = sqrt(emxerr[0]*emxerr[0] +
                       emxerr[1]*emxerr[1] +
                       emxerr[2]*emxerr[2]);

    r->n_emit = rc->in;       /* number of B-spline control points */
    r->max_dev = dev3;
    r->elapsed_sec = elapsed;
    r->available = 1;
    freeCurve(rc);
}
#else
static void strat_sisl(const polyline_t *pl, double tol, strat_result_t *r)
{
    (void)pl; (void)tol;
    r->available = 0;
}
#endif

/* ============================================================ */
/* Strategy: geomdl (Python subprocess)                          */
/* ============================================================ */

static const char *g_geomdl_python  = NULL;   /* venv python interpreter */
static const char *g_geomdl_script  = NULL;   /* path to geomdl_fit.py   */

static void strat_geomdl(const polyline_t *pl, double tol, strat_result_t *r)
{
    (void)tol;
    if (!g_geomdl_python || !g_geomdl_script) { r->available = 0; return; }

    int p_in[2], p_out[2];
    if (pipe(p_in) || pipe(p_out)) { r->available = 0; return; }

    double t0 = now_sec();
    pid_t pid = fork();
    if (pid < 0) {
        close(p_in[0]); close(p_in[1]); close(p_out[0]); close(p_out[1]);
        r->available = 0; return;
    }
    if (pid == 0) {
        /* child: stdin from p_in[0], stdout to p_out[1] */
        dup2(p_in[0], 0);
        dup2(p_out[1], 1);
        close(p_in[0]);  close(p_in[1]);
        close(p_out[0]); close(p_out[1]);
        execl(g_geomdl_python, g_geomdl_python, g_geomdl_script, (char*)NULL);
        _exit(127);
    }
    /* parent */
    close(p_in[0]);
    close(p_out[1]);

    FILE *out = fdopen(p_in[1], "w");
    fprintf(out, "%d\n3\n0\n", pl->n);
    for (int i = 0; i < pl->n; i++)
        fprintf(out, "%.6f %.6f %.6f\n", pl->pts[i].x, pl->pts[i].y, pl->pts[i].z);
    fclose(out);

    FILE *in = fdopen(p_out[0], "r");
    int n_ctrl = 0; double max_dev = 0.0; double elapsed_ms = 0.0;
    if (fscanf(in, "%d %lf %lf", &n_ctrl, &max_dev, &elapsed_ms) != 3) {
        n_ctrl = -1; max_dev = NAN; elapsed_ms = 0;
    }
    fclose(in);

    int status = 0;
    waitpid(pid, &status, 0);
    double total = now_sec() - t0;

    /* Use elapsed reported by Python (excludes fork/IPC overhead). */
    r->n_emit = n_ctrl;
    r->max_dev = max_dev;
    r->elapsed_sec = elapsed_ms > 0 ? elapsed_ms / 1000.0 : total;
    r->available = (n_ctrl > 0) ? 1 : 0;
}

/* ============================================================ */
/* Strategy: liscio                                              */
/* ============================================================ */

typedef struct { long n_prim; double max_dev; } liscio_collect_t;
static liscio_collect_t g_ll;

static void ll_cb(const liscio_primitive_t *p, void *u)
{
    (void)u;
    g_ll.n_prim++;
    double d = 0;
    if (p->type == LISCIO_PRIM_BEZIER)      d = p->bezier_max_deviation;
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
        liscio_pose_t s = {pl->pts[i-1].x, pl->pts[i-1].y, pl->pts[i-1].z, 0,0,0,0,0,0};
        liscio_pose_t e = {pl->pts[i].x,   pl->pts[i].y,   pl->pts[i].z,   0,0,0,0,0,0};
        liscio_add_line(ctx, &s, &e, 1000.0, i);
    }
    liscio_flush(ctx);
    liscio_destroy(ctx);
    r->n_emit = (int)g_ll.n_prim;
    r->max_dev = g_ll.max_dev;
    r->elapsed_sec = now_sec() - t0;
    r->available = 1;
}

/* ============================================================ */
/* Driver                                                        */
/* ============================================================ */

static const char *strat_names[] = {
    "baseline", "dp-xyz", "schneider", "sisl", "geomdl", "liscio"
};
#define N_STRAT 6
enum { S_BASE=0, S_DP=1, S_SCHN=2, S_SISL=3, S_GEOMDL=4, S_LISCIO=5 };

static void print_row(int strat, int n_input, const strat_result_t *r)
{
    if (!r->available) {
        printf("  %-10s  N/A\n", strat_names[strat]);
        return;
    }
    if (r->n_emit < 0) {
        printf("  %-10s  FAIL\n", strat_names[strat]);
        return;
    }
    double ratio = r->n_emit > 0 ? (double)n_input / r->n_emit : 0;
    double us_per = n_input > 0 ? r->elapsed_sec * 1e6 / n_input : 0;
    printf("  %-10s  segs=%-6d  ratio=%6.2fx  max_dev=%.4f  %.3f ms  (%.2f µs/pt)\n",
        strat_names[strat], r->n_emit, ratio, r->max_dev,
        r->elapsed_sec * 1000.0, us_per);
}

int main(int argc, char **argv)
{
    double tol = 0.05;
    int do_schn = 1, do_sisl = 1, do_geomdl = 1;
    int i0 = 1;
    for (int i = 1; i < argc; i++) {
        if (strncmp(argv[i], "--tol_xyz=", 10) == 0) {
            tol = atof(argv[i]+10); i0 = i+1;
        } else if (strcmp(argv[i], "--no-schneider") == 0) {
            do_schn = 0; i0 = i+1;
        } else if (strcmp(argv[i], "--no-sisl") == 0) {
            do_sisl = 0; i0 = i+1;
        } else if (strcmp(argv[i], "--no-geomdl") == 0) {
            do_geomdl = 0; i0 = i+1;
        } else break;
    }
    if (i0 >= argc) {
        fprintf(stderr,
            "usage: %s [--tol_xyz=N] [--no-schneider] [--no-sisl] "
            "[--no-geomdl] <ngc>...\n", argv[0]);
        return 1;
    }

    /* geomdl wiring: env LISCIO_GEOMDL_PYTHON / LISCIO_GEOMDL_SCRIPT */
    if (do_geomdl) {
        g_geomdl_python = getenv("LISCIO_GEOMDL_PYTHON");
        g_geomdl_script = getenv("LISCIO_GEOMDL_SCRIPT");
        if (!g_geomdl_python || !g_geomdl_script) {
            fprintf(stderr, "[info] geomdl skipped (LISCIO_GEOMDL_PYTHON/SCRIPT not set)\n");
            do_geomdl = 0;
        }
    }

    printf("# liscio cross-library comparator (tol_xyz=%.4f mm)\n", tol);
    printf("# strategies: baseline / dp-xyz / schneider(3D port) / sisl / geomdl / liscio\n#\n");

    long grand_in = 0;
    long grand_out[N_STRAT] = {0};
    double grand_t[N_STRAT] = {0};
    double grand_max_dev[N_STRAT] = {0};
    int    grand_avail[N_STRAT] = {0};
    int    files_ok = 0;

    for (int i = i0; i < argc; i++) {
        static polyline_t pl;
        if (load_ngc(argv[i], &pl) != 0 || pl.n < 4) {
            fprintf(stderr, "skip %s\n", argv[i]);
            continue;
        }
        const char *name = strrchr(argv[i], '/');
        name = name ? name+1 : argv[i];
        printf("%s  (waypoints=%d, raw_len=%.3f)\n", name, pl.n, pl.raw_length);

        strat_result_t res[N_STRAT] = {0};
        strat_baseline(&pl,  &res[S_BASE]);
        strat_dp(&pl, tol,   &res[S_DP]);
        if (do_schn)   strat_schneider(&pl, tol, &res[S_SCHN]);
        if (do_sisl)   strat_sisl(&pl, tol, &res[S_SISL]);
        if (do_geomdl) strat_geomdl(&pl, tol, &res[S_GEOMDL]);
        strat_liscio(&pl, tol, &res[S_LISCIO]);

        for (int s = 0; s < N_STRAT; s++) print_row(s, pl.n - 1, &res[s]);
        printf("\n");

        grand_in += pl.n - 1;
        for (int s = 0; s < N_STRAT; s++) {
            if (!res[s].available) continue;
            grand_avail[s] = 1;
            if (res[s].n_emit > 0) grand_out[s] += res[s].n_emit;
            grand_t[s]  += res[s].elapsed_sec;
            if (res[s].max_dev > grand_max_dev[s]) grand_max_dev[s] = res[s].max_dev;
        }
        files_ok++;
    }

    /* Markdown summary table */
    printf("\n## Cross-library summary  (tol_xyz=%.4f, %d files)\n\n", tol, files_ok);
    printf("| strategy   | total segs | compression | worst max_dev | total time | µs/wp |\n");
    printf("|------------|-----------:|------------:|--------------:|-----------:|------:|\n");
    for (int s = 0; s < N_STRAT; s++) {
        if (!grand_avail[s]) {
            printf("| %-10s | _N/A_ | _N/A_ | _N/A_ | _N/A_ | _N/A_ |\n", strat_names[s]);
            continue;
        }
        double ratio   = grand_out[s] > 0 ? (double)grand_in / grand_out[s] : 0;
        double us_per  = grand_in > 0 ? grand_t[s] * 1e6 / grand_in : 0;
        printf("| %-10s | %10ld | %10.2fx | %13.4f | %8.3f s | %5.2f |\n",
               strat_names[s], grand_out[s], ratio, grand_max_dev[s],
               grand_t[s], us_per);
    }
    return 0;
}
