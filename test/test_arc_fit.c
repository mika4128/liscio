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
* Description: Synthetic circle → arc-fit verification test.
*              Generates N points along a known circle in arbitrary
*              plane, with optional noise < tol.  Feeds to
*              liscio_add_line; expects a single ARC primitive
*              whose center matches within tolerance.
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

static int g_line_count = 0;
static int g_arc_count = 0;
static liscio_primitive_t g_last_arc;

static void on_prim(const liscio_primitive_t *p, void *user)
{
    (void)user;
    if (p->type == LISCIO_PRIM_LINE) {
        g_line_count++;
    } else if (p->type == LISCIO_PRIM_ARC) {
        g_arc_count++;
        g_last_arc = *p;
    }
}

/* Generate N points along circle centered at (cx,cy,cz), radius r, in
 * plane defined by basis vectors (ux,uy,uz) and (vx,vy,vz).  angles
 * in [a0, a1].  abc/uvw vary linearly from 0 to end values. */
static void gen_arc(liscio_ctx_t *ctx,
                    double cx, double cy, double cz, double r,
                    double ux, double uy, double uz,
                    double vx, double vy, double vz,
                    double a0, double a1, int n,
                    double abc_end, double uvw_end)
{
    for (int i = 0; i < n - 1; i++) {
        double t0 = (double)i / (n - 1);
        double t1 = (double)(i + 1) / (n - 1);
        double ang0 = a0 + t0 * (a1 - a0);
        double ang1 = a0 + t1 * (a1 - a0);
        liscio_pose_t s = {
            cx + r * (cos(ang0)*ux + sin(ang0)*vx),
            cy + r * (cos(ang0)*uy + sin(ang0)*vy),
            cz + r * (cos(ang0)*uz + sin(ang0)*vz),
            t0 * abc_end, 0, 0, t0 * uvw_end, 0, 0
        };
        liscio_pose_t e = {
            cx + r * (cos(ang1)*ux + sin(ang1)*vx),
            cy + r * (cos(ang1)*uy + sin(ang1)*vy),
            cz + r * (cos(ang1)*uz + sin(ang1)*vz),
            t1 * abc_end, 0, 0, t1 * uvw_end, 0, 0
        };
        liscio_add_line(ctx, &s, &e, 1000.0, i);
    }
}

static int test_xy_full_circle(void)
{
    printf("TEST xy_full_circle... ");
    fflush(stdout);

    liscio_cfg_t cfg;
    liscio_cfg_default(&cfg);
    cfg.tol_xyz = 0.01;
    liscio_ctx_t *ctx = liscio_create(&cfg);
    assert(ctx);
    liscio_set_callback(ctx, on_prim, NULL);
    g_line_count = g_arc_count = 0;

    /* 3/4 circle in XY plane, radius 10, center (5,5,0). */
    gen_arc(ctx, 5, 5, 0, 10,
            1, 0, 0,   0, 1, 0,
            0, 1.5 * M_PI, 30, 0, 0);
    liscio_flush(ctx);

    printf("lines=%d arcs=%d ", g_line_count, g_arc_count);
    if (g_arc_count < 1 || g_line_count > 0) {
        printf("FAIL\n"); liscio_destroy(ctx); return 1;
    }
    double ec = fabs(g_last_arc.cx - 5) + fabs(g_last_arc.cy - 5) +
                fabs(g_last_arc.cz - 0);
    if (ec > 0.05) {
        printf("center err=%.4f FAIL\n", ec); liscio_destroy(ctx);
        return 1;
    }
    if (fabs(g_last_arc.radius - 10) > 0.05) {
        printf("radius err=%.4f FAIL\n", fabs(g_last_arc.radius - 10));
        liscio_destroy(ctx); return 1;
    }
    printf("PASS (c=(%.3f,%.3f,%.3f) r=%.3f)\n",
        g_last_arc.cx, g_last_arc.cy, g_last_arc.cz, g_last_arc.radius);
    liscio_destroy(ctx);
    return 0;
}

static int test_xz_arc_with_rotary(void)
{
    printf("TEST xz_arc_with_rotary... ");
    fflush(stdout);

    liscio_cfg_t cfg;
    liscio_cfg_default(&cfg);
    cfg.tol_xyz = 0.005;
    cfg.tol_abc = 0.1;
    liscio_ctx_t *ctx = liscio_create(&cfg);
    assert(ctx);
    liscio_set_callback(ctx, on_prim, NULL);
    g_line_count = g_arc_count = 0;

    /* Half circle in XZ plane at Y=3, with C varying linearly 0→90 */
    gen_arc(ctx, 0, 3, 0, 5,
            1, 0, 0,   0, 0, 1,
            -M_PI/2, M_PI/2, 20, 90.0, 0);
    liscio_flush(ctx);

    printf("lines=%d arcs=%d ", g_line_count, g_arc_count);
    if (g_arc_count < 1 || g_line_count > 0) {
        printf("FAIL\n"); liscio_destroy(ctx); return 1;
    }
    printf("PASS (c=(%.3f,%.3f,%.3f) r=%.3f)\n",
        g_last_arc.cx, g_last_arc.cy, g_last_arc.cz, g_last_arc.radius);
    liscio_destroy(ctx);
    return 0;
}

static int test_collinear_line(void)
{
    printf("TEST collinear_line... ");
    fflush(stdout);

    liscio_ctx_t *ctx = liscio_create(NULL);
    assert(ctx);
    liscio_set_callback(ctx, on_prim, NULL);
    g_line_count = g_arc_count = 0;

    /* 10 colinear points along X. Should NOT fit as arc (degenerate) →
     * fall back to LINE primitive. */
    for (int i = 0; i < 10; i++) {
        liscio_pose_t s = { i*1.0, 0, 0, 0,0,0, 0,0,0 };
        liscio_pose_t e = { (i+1)*1.0, 0, 0, 0,0,0, 0,0,0 };
        liscio_add_line(ctx, &s, &e, 1000.0, i);
    }
    liscio_flush(ctx);

    printf("lines=%d arcs=%d ", g_line_count, g_arc_count);
    if (g_arc_count > 0) {
        printf("FAIL (unexpected arc fit on colinear)\n");
        liscio_destroy(ctx); return 1;
    }
    if (g_line_count < 1) {
        printf("FAIL (no primitive)\n");
        liscio_destroy(ctx); return 1;
    }
    printf("PASS\n");
    liscio_destroy(ctx);
    return 0;
}

static int test_noisy_circle_within_tol(void)
{
    printf("TEST noisy_circle_within_tol... ");
    fflush(stdout);

    liscio_cfg_t cfg;
    liscio_cfg_default(&cfg);
    cfg.tol_xyz = 0.02;
    liscio_ctx_t *ctx = liscio_create(&cfg);
    assert(ctx);
    liscio_set_callback(ctx, on_prim, NULL);
    g_line_count = g_arc_count = 0;

    /* Circle + small noise < tol/3. Must fit. */
    srand(42);
    double r = 50.0;
    int n = 40;
    double prev_x = r, prev_y = 0, prev_z = 0;
    for (int i = 0; i < n; i++) {
        double ang = (double)(i + 1) / n * 2.0 * M_PI * 0.5;  /* half circle */
        double jx = ((rand() & 0x7fff) / 32767.0 - 0.5) * (cfg.tol_xyz * 0.5);
        double jy = ((rand() & 0x7fff) / 32767.0 - 0.5) * (cfg.tol_xyz * 0.5);
        double ex = r * cos(ang) + jx;
        double ey = r * sin(ang) + jy;
        liscio_pose_t s = { prev_x, prev_y, prev_z, 0,0,0, 0,0,0 };
        liscio_pose_t e = { ex, ey, 0, 0,0,0, 0,0,0 };
        liscio_add_line(ctx, &s, &e, 1000.0, i);
        prev_x = ex; prev_y = ey;
    }
    liscio_flush(ctx);

    printf("lines=%d arcs=%d ", g_line_count, g_arc_count);
    if (g_arc_count < 1) {
        printf("FAIL (no arc fit)\n");
        liscio_destroy(ctx); return 1;
    }
    printf("PASS\n");
    liscio_destroy(ctx);
    return 0;
}

static int test_sharp_corner_breaks_arc(void)
{
    printf("TEST sharp_corner_breaks_arc... ");
    fflush(stdout);

    liscio_ctx_t *ctx = liscio_create(NULL);
    assert(ctx);
    liscio_set_callback(ctx, on_prim, NULL);
    g_line_count = g_arc_count = 0;

    /* 5 points along X, then sharp 90° corner, then 5 along Y.
     * Expect two line primitives (or a line + flush), no arc. */
    liscio_pose_t p0 = {0,0,0,0,0,0,0,0,0};
    liscio_add_line(ctx, &p0, &(liscio_pose_t){1,0,0,0,0,0,0,0,0}, 1000, 1);
    liscio_add_line(ctx, &(liscio_pose_t){1,0,0,0,0,0,0,0,0},
                         &(liscio_pose_t){2,0,0,0,0,0,0,0,0}, 1000, 2);
    liscio_add_line(ctx, &(liscio_pose_t){2,0,0,0,0,0,0,0,0},
                         &(liscio_pose_t){3,0,0,0,0,0,0,0,0}, 1000, 3);
    /* Sharp corner: same end = next start, but direction changes 90°. */
    liscio_add_line(ctx, &(liscio_pose_t){3,0,0,0,0,0,0,0,0},
                         &(liscio_pose_t){3,1,0,0,0,0,0,0,0}, 1000, 4);
    liscio_add_line(ctx, &(liscio_pose_t){3,1,0,0,0,0,0,0,0},
                         &(liscio_pose_t){3,2,0,0,0,0,0,0,0}, 1000, 5);
    liscio_add_line(ctx, &(liscio_pose_t){3,2,0,0,0,0,0,0,0},
                         &(liscio_pose_t){3,3,0,0,0,0,0,0,0}, 1000, 6);
    liscio_flush(ctx);

    printf("lines=%d arcs=%d ", g_line_count, g_arc_count);
    /* Acceptable: arc=0, lines>=1 (the corner may fragment into >1). */
    if (g_arc_count > 0) {
        printf("FAIL (unexpected arc at corner)\n");
        liscio_destroy(ctx); return 1;
    }
    if (g_line_count < 1) {
        printf("FAIL (no primitive)\n");
        liscio_destroy(ctx); return 1;
    }
    printf("PASS\n");
    liscio_destroy(ctx);
    return 0;
}

static int test_helical_arc(void)
{
    printf("TEST helical_arc (non-planar 3D helix)... "); fflush(stdout);
    /* Helix: X=r·cos, Y=r·sin, Z=pitch·θ.  Non-planar → arc fit
     * should REJECT (plane verification fails). */
    liscio_cfg_t cfg; liscio_cfg_default(&cfg);
    cfg.tol_xyz = 0.01;
    cfg.corner_angle_deg = 90.0;
    liscio_ctx_t *ctx = liscio_create(&cfg); assert(ctx);
    liscio_set_callback(ctx, on_prim, NULL);
    g_line_count = g_arc_count = 0;

    int n = 24;
    double r = 10.0;
    double px = r, py = 0, pz = 0;
    for (int i = 1; i <= n; i++) {
        double t = (double)i / n * M_PI;
        double ex = r * cos(t);
        double ey = r * sin(t);
        double ez = 5.0 * t / (2.0 * M_PI);
        liscio_pose_t s = { px, py, pz, 0,0,0, 0,0,0 };
        liscio_pose_t e = { ex, ey, ez, 0,0,0, 0,0,0 };
        liscio_add_line(ctx, &s, &e, 1000.0, i);
        px = ex; py = ey; pz = ez;
    }
    liscio_flush(ctx);

    printf("lines=%d arcs=%d ", g_line_count, g_arc_count);
    if (g_arc_count > 0) {
        printf("FAIL (helix mistakenly fit as arc)\n");
        liscio_destroy(ctx); return 1;
    }
    printf("PASS (no false-positive arc)\n");
    liscio_destroy(ctx);
    return 0;
}

static int test_very_small_circle(void)
{
    printf("TEST very_small_circle (r=0.5mm)... "); fflush(stdout);
    liscio_cfg_t cfg; liscio_cfg_default(&cfg);
    cfg.tol_xyz = 0.005;
    liscio_ctx_t *ctx = liscio_create(&cfg); assert(ctx);
    liscio_set_callback(ctx, on_prim, NULL);
    g_line_count = g_arc_count = 0;

    double r = 0.5;
    int n = 20;
    double px = r, py = 0;
    for (int i = 1; i <= n; i++) {
        double t = (double)i / n * M_PI;
        double ex = r * cos(t), ey = r * sin(t);
        liscio_pose_t s = { px, py, 0, 0,0,0, 0,0,0 };
        liscio_pose_t e = { ex, ey, 0, 0,0,0, 0,0,0 };
        liscio_add_line(ctx, &s, &e, 1000.0, i);
        px = ex; py = ey;
    }
    liscio_flush(ctx);

    printf("lines=%d arcs=%d ", g_line_count, g_arc_count);
    if (g_arc_count < 1) {
        printf("FAIL\n"); liscio_destroy(ctx); return 1;
    }
    printf("PASS (r=%.4f)\n", g_last_arc.radius);
    liscio_destroy(ctx);
    return 0;
}

static int test_long_arc_semicircle(void)
{
    printf("TEST long_arc_semicircle (100 pts, r=50)... "); fflush(stdout);
    liscio_cfg_t cfg; liscio_cfg_default(&cfg);
    cfg.tol_xyz = 0.01;
    cfg.max_window = 128;
    liscio_ctx_t *ctx = liscio_create(&cfg); assert(ctx);
    liscio_set_callback(ctx, on_prim, NULL);
    g_line_count = g_arc_count = 0;

    int n = 100;
    double r = 50.0;
    double px = r, py = 0;
    for (int i = 1; i <= n; i++) {
        double t = (double)i / n * M_PI;
        double ex = r * cos(t), ey = r * sin(t);
        liscio_pose_t s = { px, py, 0, 0,0,0, 0,0,0 };
        liscio_pose_t e = { ex, ey, 0, 0,0,0, 0,0,0 };
        liscio_add_line(ctx, &s, &e, 1000.0, i);
        px = ex; py = ey;
    }
    liscio_flush(ctx);

    printf("lines=%d arcs=%d ", g_line_count, g_arc_count);
    if (g_arc_count < 1) {
        printf("FAIL\n"); liscio_destroy(ctx); return 1;
    }
    printf("PASS (r=%.3f angle=%.4f)\n",
        g_last_arc.radius, g_last_arc.arc_angle);
    liscio_destroy(ctx);
    return 0;
}

#include "test_timing.h"

int main(void)
{
    LISCIO_TEST_TIMER_START(t0);
    int fails = 0;
    fails += test_collinear_line();
    fails += test_sharp_corner_breaks_arc();
    fails += test_xy_full_circle();
    fails += test_xz_arc_with_rotary();
    fails += test_noisy_circle_within_tol();
    fails += test_helical_arc();
    fails += test_very_small_circle();
    fails += test_long_arc_semicircle();
    printf("---- liscio arc tests: %s (%d failure%s, time=%.3f ms) ----\n",
           fails ? "FAIL" : "OK", fails, fails == 1 ? "" : "s",
           LISCIO_TEST_TIMER_MS(t0));
    return fails ? 1 : 0;
}
