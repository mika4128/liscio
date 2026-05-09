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
* Description: Dump waypoint input + emitted primitives to CSV for
*              plotting.  Reads an NGC file, feeds through liscio,
*              writes two CSVs:
*                - <ngc>.input.csv   raw G-code waypoints
*                - <ngc>.prims.csv   sampled points along emitted
*                                    primitives (30 samples each)
*              Usage: test_dump_csv <ngc> [--tol_xyz=N] [--tol_abc=N]
* Author:      杨阳 (Yang Yang) <mika-net@outlook.com>
* License:     MIT (SPDX-License-Identifier: MIT)
* Copyright (c) 2026 杨阳 (Yang Yang)
********************************************************************/

#define _USE_MATH_DEFINES
#include "liscio/liscio.h"
#include "test_timing.h"
#include "ngc_parser.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static FILE *g_prim_csv = NULL;
static int   g_prim_idx = 0;
static int   g_samples  = 30;

/* NGC parsing delegated to shared ngc_parser (handles G17/G18/G19, G20/G21,
 * G90/G91, O-word while, named variables). */

/* --- primitive sampling → CSV row --- */
static void sample_line(const liscio_primitive_t *p)
{
    for (int i = 0; i <= g_samples; i++) {
        double t = (double)i / g_samples;
        double x = p->start.x + t*(p->end.x - p->start.x);
        double y = p->start.y + t*(p->end.y - p->start.y);
        double z = p->start.z + t*(p->end.z - p->start.z);
        fprintf(g_prim_csv, "%d,LINE,%d,%.6f,%.6f,%.6f\n",
            g_prim_idx, i, x, y, z);
    }
}

static void sample_arc(const liscio_primitive_t *p)
{
    /* Arc param: start = vec from center at t=0, rotate about normal
     * by t·arc_angle to reach end. */
    double vx = p->start.x - p->cx;
    double vy = p->start.y - p->cy;
    double vz = p->start.z - p->cz;

    /* Build orthonormal (u, v) basis in plane: u = start-radial; v = normal × u */
    double un = sqrt(vx*vx + vy*vy + vz*vz);
    if (un < 1e-15) return;
    double ux = vx/un, uy = vy/un, uz = vz/un;
    double vx2 = p->ny*uz - p->nz*uy;
    double vy2 = p->nz*ux - p->nx*uz;
    double vz2 = p->nx*uy - p->ny*ux;

    const char *label = (p->type == LISCIO_PRIM_HELIX) ? "HELIX" : "ARC";
    for (int i = 0; i <= g_samples; i++) {
        double t = (double)i / g_samples;
        double ang = t * p->arc_angle;
        double cx = p->cx + p->radius * (cos(ang)*ux + sin(ang)*vx2);
        double cy = p->cy + p->radius * (cos(ang)*uy + sin(ang)*vy2);
        double cz = p->cz + p->radius * (cos(ang)*uz + sin(ang)*vz2);
        /* HELIX: add axial advance along normal (pitch·θ).  For planar
         * ARC pitch is 0; the legacy uz≈vz2≈0 fallback for ARC keeps
         * single-arc passthrough behaviour. */
        if (p->type == LISCIO_PRIM_HELIX) {
            double axial = ang * p->pitch;
            cx += axial * p->nx;
            cy += axial * p->ny;
            cz += axial * p->nz;
        } else {
            double plane_z_var = fabs(uz) + fabs(vz2);
            if (plane_z_var < 1e-12) {
                cz = p->start.z + t * (p->end.z - p->start.z);
            }
        }
        fprintf(g_prim_csv, "%d,%s,%d,%.6f,%.6f,%.6f\n",
            g_prim_idx, label, i, cx, cy, cz);
    }
}

static void sample_bezier(const liscio_primitive_t *p)
{
    for (int i = 0; i <= g_samples; i++) {
        double t = (double)i / g_samples;
        double x, y, z;
        liscio_rbezier_eval(p, t, &x, &y, &z);
        fprintf(g_prim_csv, "%d,BEZIER,%d,%.6f,%.6f,%.6f\n",
            g_prim_idx, i, x, y, z);
    }
}

static void sample_spline(const liscio_primitive_t *p)
{
    for (int i = 0; i <= g_samples; i++) {
        double t = (double)i / g_samples;
        double x, y, z;
        liscio_bspline_eval(p, t, &x, &y, &z);
        fprintf(g_prim_csv, "%d,SPLINE,%d,%.6f,%.6f,%.6f\n",
            g_prim_idx, i, x, y, z);
    }
}

static void on_prim(const liscio_primitive_t *p, void *user)
{
    (void)user;
    switch (p->type) {
        case LISCIO_PRIM_LINE:   sample_line(p); break;
        case LISCIO_PRIM_HELIX:  /* same fields as ARC */
        case LISCIO_PRIM_ARC:    sample_arc(p); break;
        case LISCIO_PRIM_BEZIER: sample_bezier(p); break;
        case LISCIO_PRIM_SPLINE: sample_spline(p); break;
    }
    g_prim_idx++;
}

/* --- parser bridge --- */
typedef struct {
    liscio_ctx_t *ctx;
    FILE *fi;      /* input waypoints CSV */
    FILE *fr;      /* G0 rapids CSV */
    int   in_seq;
} dump_bridge_t;

static void on_move(const ngc_move_t *m, void *user)
{
    dump_bridge_t *db = (dump_bridge_t *)user;
    if (m->is_rapid) {
        fprintf(db->fr, "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
            db->in_seq, m->start.x, m->start.y, m->start.z,
            m->end.x, m->end.y, m->end.z);
        db->in_seq++;
        liscio_flush(db->ctx);
        return;
    }
    fprintf(db->fi, "%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
        db->in_seq, m->end.x, m->end.y, m->end.z,
        m->end.a, m->end.b, m->end.c,
        m->end.u, m->end.v, m->end.w);
    db->in_seq++;
    if (m->type == NGC_MOVE_LINE) {
        liscio_pose_t s = {m->start.x, m->start.y, m->start.z,
                            m->start.a, m->start.b, m->start.c,
                            m->start.u, m->start.v, m->start.w};
        liscio_pose_t e = {m->end.x, m->end.y, m->end.z,
                            m->end.a, m->end.b, m->end.c,
                            m->end.u, m->end.v, m->end.w};
        liscio_add_line(db->ctx, &s, &e, m->feedrate, m->line_no);
    } else {
        liscio_primitive_t a; memset(&a, 0, sizeof(a));
        a.type  = LISCIO_PRIM_ARC;
        a.start = (liscio_pose_t){m->start.x, m->start.y, m->start.z,
                                  m->start.a, m->start.b, m->start.c,
                                  m->start.u, m->start.v, m->start.w};
        a.end   = (liscio_pose_t){m->end.x, m->end.y, m->end.z,
                                  m->end.a, m->end.b, m->end.c,
                                  m->end.u, m->end.v, m->end.w};
        a.cx=m->cx; a.cy=m->cy; a.cz=m->cz;
        a.nx=m->nx; a.ny=m->ny; a.nz=m->nz;
        a.radius=m->radius;
        a.arc_angle=m->arc_angle;
        liscio_add_arc(db->ctx, &a, m->feedrate, m->line_no);
    }
}

/* --- main --- */
static int dump_file(const char *path, const liscio_cfg_t *cfg)
{
    char in_csv[1024], pr_csv[1024], rp_csv[1024];
    snprintf(in_csv, sizeof(in_csv), "%s.input.csv",  path);
    snprintf(pr_csv, sizeof(pr_csv), "%s.prims.csv",  path);
    snprintf(rp_csv, sizeof(rp_csv), "%s.rapids.csv", path);

    FILE *fi = fopen(in_csv, "w");
    FILE *fr = fopen(rp_csv, "w");
    g_prim_csv = fopen(pr_csv, "w");
    if (!fi || !g_prim_csv || !fr) {
        fprintf(stderr, "cannot open output CSVs\n"); return -1;
    }
    fprintf(fi,         "seq,x,y,z,a,b,c,u,v,w\n");
    fprintf(g_prim_csv, "prim,type,seq,x,y,z\n");
    fprintf(fr,         "seg,sx,sy,sz,ex,ey,ez\n");

    liscio_ctx_t *ctx = liscio_create(cfg);
    liscio_set_callback(ctx, on_prim, NULL);
    g_prim_idx = 0;

    dump_bridge_t db = { ctx, fi, fr, 0 };
    ngc_stats_t st;
    int rc = ngc_parse_file(path, on_move, &db, &st);
    liscio_flush(ctx);
    fclose(fi); fclose(fr); fclose(g_prim_csv);
    liscio_destroy(ctx);
    if (rc != 0) return -1;

    printf("wrote %s (%d waypoints), %s (%d primitives), %s\n",
        in_csv, db.in_seq, pr_csv, g_prim_idx, rp_csv);
    return 0;
}

int main(int argc, char **argv)
{
    liscio_cfg_t cfg; liscio_cfg_default(&cfg);
    int i0 = 1;
    for (int i = 1; i < argc; i++) {
        if (strncmp(argv[i], "--tol_xyz=", 10) == 0) {
            cfg.tol_xyz = atof(argv[i] + 10); i0 = i + 1;
        } else if (strncmp(argv[i], "--tol_abc=", 10) == 0) {
            cfg.tol_abc = atof(argv[i] + 10); i0 = i + 1;
        } else if (strncmp(argv[i], "--samples=", 10) == 0) {
            g_samples = atoi(argv[i] + 10); i0 = i + 1;
        } else if (strncmp(argv[i], "--arc_subseg=", 13) == 0) {
            cfg.arc_subseg_samples = atoi(argv[i] + 13); i0 = i + 1;
        } else if (strncmp(argv[i], "--arc_merge=", 12) == 0) {
            int m = atoi(argv[i] + 12);
            liscio_cfg_set_arc_merge(&cfg,
                (m == 2) ? LISCIO_ARC_MERGE_HELIX :
                (m == 1) ? LISCIO_ARC_MERGE_ARC :
                           LISCIO_ARC_MERGE_OFF, 0, 0);
            i0 = i + 1;
        } else if (strncmp(argv[i], "--arc_to_helix=", 15) == 0) {
            cfg.arc_to_helix_max = atoi(argv[i] + 15); i0 = i + 1;
        } else break;
    }
    if (i0 >= argc) {
        fprintf(stderr,
            "usage: %s [--tol_xyz=N] [--tol_abc=N] [--samples=N] "
            "[--arc_subseg=N] [--arc_merge=0|1|2] <ngc>...\n",
            argv[0]);
        return 1;
    }
    LISCIO_TEST_TIMER_START(__t0);
    int rc = 0;
    for (int i = i0; i < argc; i++) {
        if (dump_file(argv[i], &cfg) != 0) rc = 1;
    }
    printf("---- liscio dump_csv: %d file%s, %.3f ms total ----\n",
        argc - i0, (argc - i0 == 1 ? "" : "s"), LISCIO_TEST_TIMER_MS(__t0));
    return rc;
}
