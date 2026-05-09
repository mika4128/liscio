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
* Description: B-spline fit + evaluator verification test.
*              Covers:
*                - complex polar curve → single SPLINE primitive
*                  (verifies fit + eval round-trip within tol)
*                - enable_bspline=0 fallback to composite Bezier
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

static int g_line = 0, g_arc = 0, g_bezier = 0, g_spline = 0;
static liscio_primitive_t g_last_spline;

static void on_prim(const liscio_primitive_t *p, void *user)
{
    (void)user;
    switch (p->type) {
        case LISCIO_PRIM_LINE:   g_line++; break;
        case LISCIO_PRIM_ARC:    g_arc++; break;
        case LISCIO_PRIM_BEZIER: g_bezier++; break;
        case LISCIO_PRIM_SPLINE: g_spline++; g_last_spline = *p; break;
    }
}

static int test_polar_single_spline(void)
{
    printf("TEST polar_single_spline... "); fflush(stdout);

    /* 40-pt polar curve r·(1+0.2·sin(2θ)) over θ ∈ [0, π].
     * Smoother than the test_bezier_fit composite case — single
     * B-spline with 6-10 ctrl pts should fit within moderate tol. */
    liscio_cfg_t cfg;
    liscio_cfg_default(&cfg);
    cfg.tol_xyz = 0.08;
    cfg.corner_angle_deg = 90.0;
    cfg.enable_bezier = 0;     /* force B-spline path */
    cfg.enable_composite_bezier = 0;
    cfg.bspline_init_ctrl = 6;
    liscio_ctx_t *ctx = liscio_create(&cfg);
    assert(ctx);
    liscio_set_callback(ctx, on_prim, NULL);
    g_line = g_arc = g_bezier = g_spline = 0;

    int n = 40;
    double R = 10.0;
    double th0 = 0;
    double r0 = R * (1.0 + 0.2 * sin(2.0 * th0));
    double px = r0 * cos(th0), py = r0 * sin(th0);
    for (int i = 1; i <= n; i++) {
        double th = (double)i / n * M_PI;
        double r = R * (1.0 + 0.2 * sin(2.0 * th));
        double ex = r * cos(th), ey = r * sin(th);
        liscio_pose_t s = { px, py, 0, 0,0,0, 0,0,0 };
        liscio_pose_t e = { ex, ey, 0, 0,0,0, 0,0,0 };
        liscio_add_line(ctx, &s, &e, 1000.0, i);
        px = ex; py = ey;
    }
    liscio_flush(ctx);

    printf("L=%d A=%d B=%d S=%d ", g_line, g_arc, g_bezier, g_spline);
    if (g_spline < 1) {
        printf("FAIL (no spline)\n"); liscio_destroy(ctx); return 1;
    }

    /* Verify eval matches endpoints. */
    double x, y, z;
    liscio_bspline_eval(&g_last_spline, 0.0, &x, &y, &z);
    double e_start = fabs(x - g_last_spline.start.x) + fabs(y - g_last_spline.start.y);
    liscio_bspline_eval(&g_last_spline, 1.0, &x, &y, &z);
    double e_end = fabs(x - g_last_spline.end.x) + fabs(y - g_last_spline.end.y);
    if (e_start > 1e-6 || e_end > 1e-6) {
        printf("FAIL (endpoints e_s=%.2e e_e=%.2e)\n", e_start, e_end);
        liscio_destroy(ctx); return 1;
    }
    printf("PASS (M=%d dev=%.4f)\n",
        g_last_spline.n_ctrl, g_last_spline.spline_max_deviation);
    liscio_destroy(ctx);
    return 0;
}

static int test_bspline_disabled_fallback(void)
{
    printf("TEST bspline_disabled_fallback... "); fflush(stdout);

    liscio_cfg_t cfg;
    liscio_cfg_default(&cfg);
    cfg.tol_xyz = 0.08;
    cfg.corner_angle_deg = 90.0;
    cfg.enable_bspline = 0;
    cfg.enable_bezier = 1;
    cfg.enable_composite_bezier = 1;
    liscio_ctx_t *ctx = liscio_create(&cfg);
    assert(ctx);
    liscio_set_callback(ctx, on_prim, NULL);
    g_line = g_arc = g_bezier = g_spline = 0;

    int n = 40;
    double R = 10.0;
    double px = R, py = 0;
    for (int i = 1; i <= n; i++) {
        double th = (double)i / n * M_PI;
        double r = R * (1.0 + 0.2 * sin(2.0 * th));
        double ex = r * cos(th), ey = r * sin(th);
        liscio_pose_t s = { px, py, 0, 0,0,0, 0,0,0 };
        liscio_pose_t e = { ex, ey, 0, 0,0,0, 0,0,0 };
        liscio_add_line(ctx, &s, &e, 1000.0, i);
        px = ex; py = ey;
    }
    liscio_flush(ctx);

    printf("L=%d A=%d B=%d S=%d ", g_line, g_arc, g_bezier, g_spline);
    if (g_spline > 0) {
        printf("FAIL (spline emitted despite disable)\n");
        liscio_destroy(ctx); return 1;
    }
    if (g_bezier < 1 && g_arc < 1) {
        printf("FAIL (no curve primitive)\n");
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
    fails += test_polar_single_spline();
    fails += test_bspline_disabled_fallback();
    printf("---- liscio bspline tests: %s (%d, time=%.3f ms) ----\n",
           fails ? "FAIL" : "OK", fails, LISCIO_TEST_TIMER_MS(t0));
    return fails ? 1 : 0;
}
