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
* Description: Replay real CAM NGC files through liscio.  Ships a
*              minimal subset G-code parser (G0/G1/G2/G3 with
*              X/Y/Z/A/B/C/U/V/W/I/J/K/R/F), skips comments, O-words,
*              M-codes, #vars.  Reports per file: input G01/G23 count,
*              emit (line/arc/bezier), absorb ratio, max-run length.
*              Usage: test_ngc_replay <ngc>... [--tol_xyz=N] [--tol_abc=N]
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

typedef struct {
    long input_g01;
    long input_g23;
    long emit_line;
    long emit_arc;
    long emit_bezier;
    long emit_spline;
    long absorbed_total;
    double max_ratio;
} stats_t;

static stats_t g_stats;

static void on_prim(const liscio_primitive_t *p, void *user)
{
    (void)user;
    if (p->type == LISCIO_PRIM_LINE) g_stats.emit_line++;
    else if (p->type == LISCIO_PRIM_ARC) g_stats.emit_arc++;
    else if (p->type == LISCIO_PRIM_BEZIER) g_stats.emit_bezier++;
    else if (p->type == LISCIO_PRIM_SPLINE) g_stats.emit_spline++;
    g_stats.absorbed_total += p->n_absorbed;
    double r = (double)p->n_absorbed;
    if (r > g_stats.max_ratio) g_stats.max_ratio = r;
}

/* ---------- NGC parser bridge ---------- */
typedef struct { liscio_ctx_t *ctx; } replay_bridge_t;

static void on_move(const ngc_move_t *m, void *user)
{
    replay_bridge_t *rb = (replay_bridge_t *)user;
    if (m->type == NGC_MOVE_LINE) {
        /* G0 and G1 both fed as lines (liscio decides how to absorb). */
        liscio_pose_t s = {m->start.x, m->start.y, m->start.z,
                            m->start.a, m->start.b, m->start.c,
                            m->start.u, m->start.v, m->start.w};
        liscio_pose_t e = {m->end.x, m->end.y, m->end.z,
                            m->end.a, m->end.b, m->end.c,
                            m->end.u, m->end.v, m->end.w};
        liscio_add_line(rb->ctx, &s, &e, m->feedrate, m->line_no);
        g_stats.input_g01++;
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
        liscio_add_arc(rb->ctx, &a, m->feedrate, m->line_no);
        g_stats.input_g23++;
    }
}

/* ---------- main replay ---------- */
static int replay_file(const char *path, const liscio_cfg_t *cfg)
{
    liscio_ctx_t *ctx = liscio_create(cfg);
    if (!ctx) return -1;
    memset(&g_stats, 0, sizeof(g_stats));
    liscio_set_callback(ctx, on_prim, NULL);

    replay_bridge_t rb = { ctx };
    ngc_stats_t pst;
    int rc = ngc_parse_file(path, on_move, &rb, &pst);
    if (rc != 0) { liscio_destroy(ctx); return -1; }
    liscio_flush(ctx);

    /* Report */
    long total_input = g_stats.input_g01 + g_stats.input_g23;
    long total_emit = g_stats.emit_line + g_stats.emit_arc;
    double absorb_ratio = total_input > 0 ?
        (double)total_input / (double)total_emit : 0.0;

    printf("%-40s  G01=%-6ld G23=%-4ld  L=%-5ld A=%-5ld B=%-5ld S=%-5ld  "
           "ratio=%.2fx  max_absorb=%.0f\n",
        path, g_stats.input_g01, g_stats.input_g23,
        g_stats.emit_line, g_stats.emit_arc, g_stats.emit_bezier,
        g_stats.emit_spline, absorb_ratio, g_stats.max_ratio);

    liscio_destroy(ctx);
    return 0;
}

int main(int argc, char **argv)
{
    liscio_cfg_t cfg;
    liscio_cfg_default(&cfg);

    int i0 = 1;
    for (int i = 1; i < argc; i++) {
        if (strncmp(argv[i], "--tol_xyz=", 10) == 0) {
            cfg.tol_xyz = atof(argv[i] + 10);
            i0 = i + 1;
        } else if (strncmp(argv[i], "--tol_abc=", 10) == 0) {
            cfg.tol_abc = atof(argv[i] + 10);
            i0 = i + 1;
        } else {
            break;
        }
    }

    if (i0 >= argc) {
        fprintf(stderr, "usage: %s [--tol_xyz=N] [--tol_abc=N] <ngc_file>...\n",
            argv[0]);
        return 1;
    }

    printf("liscio replay (tol_xyz=%.4f tol_abc=%.4f)\n",
        cfg.tol_xyz, cfg.tol_abc);
    printf("%-40s  %-11s %-9s  %-13s %-14s  %-11s %s\n",
        "file", "G01", "G23", "emit_line", "emit_arc",
        "ratio", "max_absorb");

    LISCIO_TEST_TIMER_START(__t0);
    int fails = 0;
    for (int i = i0; i < argc; i++) {
        if (replay_file(argv[i], &cfg) != 0) fails++;
    }
    printf("---- liscio replay: %d file%s processed in %.3f ms ----\n",
        argc - i0, (argc - i0 == 1 ? "" : "s"), LISCIO_TEST_TIMER_MS(__t0));
    return fails ? 1 : 0;
}
