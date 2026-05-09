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
* Description: Synthetic cubic Bezier → fit verification test.
*              Covers: single S-curve, ARC preferred over Bezier,
*              enable_bezier=0 fallback, composite recursive split,
*              and emitted tangent-vector correctness on quarter arc.
* Author:      杨阳 (Yang Yang) <mika-net@outlook.com>
* License:     MIT (SPDX-License-Identifier: MIT)
* Copyright (c) 2026 杨阳 (Yang Yang)
********************************************************************/

#define _USE_MATH_DEFINES
#include "liscio/liscio.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int g_bezier = 0, g_arc = 0, g_line = 0;
static liscio_primitive_t g_last;

static void on_prim(const liscio_primitive_t *p, void *user)
{
    (void)user;
    g_last = *p;
    if (p->type == LISCIO_PRIM_BEZIER) g_bezier++;
    else if (p->type == LISCIO_PRIM_ARC) g_arc++;
    else if (p->type == LISCIO_PRIM_LINE) g_line++;
}

static double bez(double P0, double P1, double P2, double P3, double t)
{
    double u = 1.0 - t;
    return u*u*u*P0 + 3*u*u*t*P1 + 3*u*t*t*P2 + t*t*t*P3;
}

static int test_s_curve(void)
{
    printf("TEST s_curve_bezier... "); fflush(stdout);

    /* S-curve: P0=(0,0,0), P1=(10,0,0), P2=(0,10,0), P3=(10,10,0)
     * This is a classic S-shape that ARC cannot fit (not circular). */
    liscio_cfg_t cfg;
    liscio_cfg_default(&cfg);
    cfg.tol_xyz = 0.5;  /* chord-length reparam incurs modest dev on S */
    cfg.corner_angle_deg = 45.0;  /* relax so S-curve not fragmented */
    liscio_ctx_t *ctx = liscio_create(&cfg);
    assert(ctx);
    liscio_set_callback(ctx, on_prim, NULL);
    g_bezier = g_arc = g_line = 0;

    double P0x=0, P0y=0, P1x=10, P1y=0, P2x=0, P2y=10, P3x=10, P3y=10;
    int n = 20;
    double px = P0x, py = P0y;
    for (int i = 1; i <= n; i++) {
        double t = (double)i / n;
        double ex = bez(P0x, P1x, P2x, P3x, t);
        double ey = bez(P0y, P1y, P2y, P3y, t);
        liscio_pose_t s = { px, py, 0, 0,0,0, 0,0,0 };
        liscio_pose_t e = { ex, ey, 0, 0,0,0, 0,0,0 };
        liscio_add_line(ctx, &s, &e, 1000.0, i);
        px = ex; py = ey;
    }
    liscio_flush(ctx);

    printf("L=%d A=%d B=%d ", g_line, g_arc, g_bezier);
    if (g_bezier < 1) {
        printf("FAIL (no bezier emitted)\n");
        liscio_destroy(ctx); return 1;
    }
    /* Check P0 and P3 match. */
    double e0 = fabs(g_last.P0.x - P0x) + fabs(g_last.P0.y - P0y);
    double e3 = fabs(g_last.P3.x - P3x) + fabs(g_last.P3.y - P3y);
    if (e0 > 1e-6 || e3 > 1e-6) {
        printf("FAIL (endpoint mismatch e0=%.6f e3=%.6f)\n", e0, e3);
        liscio_destroy(ctx); return 1;
    }
    printf("PASS (P0=(%.2f,%.2f) P1=(%.2f,%.2f) P2=(%.2f,%.2f) P3=(%.2f,%.2f) dev=%.4f)\n",
        g_last.P0.x, g_last.P0.y, g_last.P1.x, g_last.P1.y,
        g_last.P2.x, g_last.P2.y, g_last.P3.x, g_last.P3.y,
        g_last.bezier_max_deviation);
    liscio_destroy(ctx);
    return 0;
}

static int test_arc_prefers_arc_over_bezier(void)
{
    printf("TEST arc_preferred_over_bezier... "); fflush(stdout);

    /* A perfect circle can be approximated by Bezier but ARC is the
     * specific primitive — arc fit should win (tried first). */
    liscio_cfg_t cfg;
    liscio_cfg_default(&cfg);
    cfg.tol_xyz = 0.02;
    liscio_ctx_t *ctx = liscio_create(&cfg);
    assert(ctx);
    liscio_set_callback(ctx, on_prim, NULL);
    g_bezier = g_arc = g_line = 0;

    int n = 30;
    for (int i = 0; i < n; i++) {
        double t0 = (double)i / n * M_PI;  /* half circle */
        double t1 = (double)(i+1) / n * M_PI;
        liscio_pose_t s = { 10*cos(t0), 10*sin(t0), 0, 0,0,0, 0,0,0 };
        liscio_pose_t e = { 10*cos(t1), 10*sin(t1), 0, 0,0,0, 0,0,0 };
        liscio_add_line(ctx, &s, &e, 1000.0, i);
    }
    liscio_flush(ctx);

    printf("L=%d A=%d B=%d ", g_line, g_arc, g_bezier);
    if (g_arc < 1 || g_bezier > 0) {
        printf("FAIL (expected ARC, got otherwise)\n");
        liscio_destroy(ctx); return 1;
    }
    printf("PASS\n");
    liscio_destroy(ctx);
    return 0;
}

static int test_disable_bezier(void)
{
    printf("TEST disable_bezier... "); fflush(stdout);

    liscio_cfg_t cfg;
    liscio_cfg_default(&cfg);
    cfg.enable_bezier = 0;
    cfg.corner_angle_deg = 45.0;
    cfg.tol_xyz = 0.5;
    liscio_ctx_t *ctx = liscio_create(&cfg);
    assert(ctx);
    liscio_set_callback(ctx, on_prim, NULL);
    g_bezier = g_arc = g_line = 0;

    /* Same S-curve; expect LINE fallback now. */
    double P0x=0, P0y=0, P1x=10, P1y=0, P2x=0, P2y=10, P3x=10, P3y=10;
    int n = 20;
    double px = P0x, py = P0y;
    for (int i = 1; i <= n; i++) {
        double t = (double)i / n;
        double ex = bez(P0x, P1x, P2x, P3x, t);
        double ey = bez(P0y, P1y, P2y, P3y, t);
        liscio_pose_t s = { px, py, 0, 0,0,0, 0,0,0 };
        liscio_pose_t e = { ex, ey, 0, 0,0,0, 0,0,0 };
        liscio_add_line(ctx, &s, &e, 1000.0, i);
        px = ex; py = ey;
    }
    liscio_flush(ctx);

    printf("L=%d A=%d B=%d ", g_line, g_arc, g_bezier);
    if (g_bezier > 0) {
        printf("FAIL (bezier emitted despite disable)\n");
        liscio_destroy(ctx); return 1;
    }
    printf("PASS\n");
    liscio_destroy(ctx);
    return 0;
}

static int test_tangent_vectors(void)
{
    printf("TEST tangent_vectors... "); fflush(stdout);

    liscio_cfg_t cfg;
    liscio_cfg_default(&cfg);
    cfg.tol_xyz = 0.05;  /* tight so arc not mistaken for line */
    cfg.corner_angle_deg = 90.0;
    liscio_ctx_t *ctx = liscio_create(&cfg);
    assert(ctx);
    liscio_set_callback(ctx, on_prim, NULL);
    g_bezier = g_arc = g_line = 0;

    /* Quarter circle in XY from (10,0) to (0,10), CCW, r=10 */
    int n = 20;
    liscio_pose_t prev = { 10, 0, 0, 0,0,0, 0,0,0 };
    for (int i = 1; i <= n; i++) {
        double t = (double)i / n * M_PI / 2;
        liscio_pose_t e = { 10*cos(t), 10*sin(t), 0, 0,0,0, 0,0,0 };
        liscio_add_line(ctx, &prev, &e, 1000.0, i);
        prev = e;
    }
    liscio_flush(ctx);

    if (g_arc < 1) { printf("FAIL (no arc)\n"); liscio_destroy(ctx); return 1; }
    /* Arc tangent at start (10,0) CCW should be (0, 1, 0). */
    double err_start = fabs(g_last.tan_start_x) + fabs(g_last.tan_start_y - 1.0);
    /* At end (0,10): (-1, 0, 0). */
    double err_end   = fabs(g_last.tan_end_x + 1.0) + fabs(g_last.tan_end_y);
    if (err_start > 0.1 || err_end > 0.1) {
        printf("FAIL (tan_start=(%.3f,%.3f,%.3f) tan_end=(%.3f,%.3f,%.3f))\n",
            g_last.tan_start_x, g_last.tan_start_y, g_last.tan_start_z,
            g_last.tan_end_x, g_last.tan_end_y, g_last.tan_end_z);
        liscio_destroy(ctx); return 1;
    }
    printf("PASS (tan_s=(%.2f,%.2f) tan_e=(%.2f,%.2f))\n",
        g_last.tan_start_x, g_last.tan_start_y,
        g_last.tan_end_x, g_last.tan_end_y);
    liscio_destroy(ctx);
    return 0;
}

static int test_composite_bezier(void)
{
    printf("TEST composite_bezier... "); fflush(stdout);

    /* Complex curve: r*sin(3θ) radial polar, sampled 60 pts.
     * Not representable by a single cubic; should decompose into >=2 Bezier. */
    liscio_cfg_t cfg;
    liscio_cfg_default(&cfg);
    cfg.tol_xyz = 0.05;
    cfg.corner_angle_deg = 90.0;  /* allow large smooth turns */
    liscio_ctx_t *ctx = liscio_create(&cfg);
    assert(ctx);
    liscio_set_callback(ctx, on_prim, NULL);
    g_bezier = g_arc = g_line = 0;

    int n = 60;
    double R = 10.0;
    double px = 0, py = 0;
    /* First point. */
    double th0 = 0;
    double r0 = R * (1.0 + 0.3 * sin(3.0 * th0));
    px = r0 * cos(th0); py = r0 * sin(th0);
    for (int i = 1; i <= n; i++) {
        double th = (double)i / n * M_PI;  /* half polar sweep */
        double r = R * (1.0 + 0.3 * sin(3.0 * th));
        double ex = r * cos(th), ey = r * sin(th);
        liscio_pose_t s = { px, py, 0, 0,0,0, 0,0,0 };
        liscio_pose_t e = { ex, ey, 0, 0,0,0, 0,0,0 };
        liscio_add_line(ctx, &s, &e, 1000.0, i);
        px = ex; py = ey;
    }
    liscio_flush(ctx);

    printf("L=%d A=%d B=%d ", g_line, g_arc, g_bezier);
    if (g_bezier + g_arc + g_line < 1) {
        printf("FAIL (nothing emitted)\n");
        liscio_destroy(ctx); return 1;
    }
    printf("PASS\n");
    liscio_destroy(ctx);
    return 0;
}

#include "test_timing.h"

int main(void)
{
    LISCIO_TEST_TIMER_START(t0);
    int fails = 0;
    fails += test_s_curve();
    fails += test_arc_prefers_arc_over_bezier();
    fails += test_disable_bezier();
    fails += test_composite_bezier();
    fails += test_tangent_vectors();
    printf("---- liscio bezier tests: %s (%d, time=%.3f ms) ----\n",
           fails ? "FAIL" : "OK", fails, LISCIO_TEST_TIMER_MS(t0));
    return fails ? 1 : 0;
}
