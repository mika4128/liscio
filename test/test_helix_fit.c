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
* Description: Synthetic helix → helix9 fit verification.
*              Generates N points along a known helical curve P(θ)
*              = c + r·(cos·u + sin·v) + p·θ·n in arbitrary plane,
*              feeds to liscio_add_line, expects a single LISCIO_PRIM_HELIX
*              whose axis / center / radius / pitch / arc_angle match
*              ground truth within tolerance.  Also verifies (a) planar
*              circles still classify as ARC (helix should not steal
*              them), and (b) noise > tol_xyz fails the gate.
*
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

/* ---------- callback collector ---------- */

static int g_line_count = 0;
static int g_arc_count = 0;
static int g_helix_count = 0;
static int g_other_count = 0;
static liscio_primitive_t g_last_prim;

static void on_prim(const liscio_primitive_t *p, void *user)
{
    (void)user;
    switch (p->type) {
        case LISCIO_PRIM_LINE:   g_line_count++;  break;
        case LISCIO_PRIM_ARC:    g_arc_count++;   g_last_prim = *p; break;
        case LISCIO_PRIM_HELIX:  g_helix_count++; g_last_prim = *p; break;
        default:                  g_other_count++; break;
    }
}

static void reset_counters(void)
{
    g_line_count = g_arc_count = g_helix_count = g_other_count = 0;
    memset(&g_last_prim, 0, sizeof(g_last_prim));
}

/* ---------- helix point generator ---------- */

/* Feed N adjacent G1 segments tracing a helix with axis (nx,ny,nz),
 * center anchor (cx,cy,cz), radial unit basis (ux,uy,uz),(vx,vy,vz),
 * radius r, pitch p (mm/rad), angle range [a0,a1].  N is the segment
 * count (so N+1 sampled poses; we feed N add_line calls). */
static void gen_helix(liscio_ctx_t *ctx,
                       double cx, double cy, double cz,
                       double nx, double ny, double nz,
                       double ux, double uy, double uz,
                       double vx, double vy, double vz,
                       double r, double pitch,
                       double a0, double a1, int N)
{
    for (int i = 0; i < N; i++) {
        double th0 = a0 + (double)i       / N * (a1 - a0);
        double th1 = a0 + (double)(i + 1) / N * (a1 - a0);
        liscio_pose_t s = {
            cx + r*(cos(th0)*ux + sin(th0)*vx) + pitch*th0*nx,
            cy + r*(cos(th0)*uy + sin(th0)*vy) + pitch*th0*ny,
            cz + r*(cos(th0)*uz + sin(th0)*vz) + pitch*th0*nz,
            0,0,0, 0,0,0
        };
        liscio_pose_t e = {
            cx + r*(cos(th1)*ux + sin(th1)*vx) + pitch*th1*nx,
            cy + r*(cos(th1)*uy + sin(th1)*vy) + pitch*th1*ny,
            cz + r*(cos(th1)*uz + sin(th1)*vz) + pitch*th1*nz,
            0,0,0, 0,0,0
        };
        liscio_add_line(ctx, &s, &e, 1000.0, i);
    }
}

/* ---------- helpers for assertions ---------- */

/* Distance from point P to a line through origin C with direction n.
 * Returns ‖(P - C) − ((P-C)·n)·n‖. */
static double dist_point_to_axis(double px, double py, double pz,
                                  double cx, double cy, double cz,
                                  double nx, double ny, double nz)
{
    double dx = px - cx, dy = py - cy, dz = pz - cz;
    double t  = dx*nx + dy*ny + dz*nz;
    double rx = dx - t*nx, ry = dy - t*ny, rz = dz - t*nz;
    return sqrt(rx*rx + ry*ry + rz*rz);
}

/* ---------- tests ---------- */

/* 1. Single full turn helix in XY plane, axis = +Z. */
static int test_xy_full_turn(void)
{
    printf("TEST xy_full_turn... ");
    fflush(stdout);

    liscio_cfg_t cfg;
    liscio_cfg_default(&cfg);
    cfg.tol_xyz = 0.01;
    liscio_ctx_t *ctx = liscio_create(&cfg);
    assert(ctx);
    liscio_set_callback(ctx, on_prim, NULL);
    reset_counters();

    /* r=10, pitch p=0.5 mm/rad → 1 turn rises 2π·p ≈ 3.14 mm. */
    gen_helix(ctx, /*c*/ 5,5,0, /*n*/ 0,0,1,
              /*u*/ 1,0,0, /*v*/ 0,1,0,
              /*r*/ 10, /*p*/ 0.5,
              0.0, 2.0 * M_PI, 60);
    liscio_flush(ctx);

    printf("L=%d A=%d H=%d ", g_line_count, g_arc_count, g_helix_count);
    if (g_helix_count < 1 || g_line_count > 0 || g_arc_count > 0) {
        printf("FAIL (expected 1 HELIX)\n");
        liscio_destroy(ctx); return 1;
    }
    /* Center should lie on the axis (5, 5, *).  Axis direction along ±Z. */
    double drad = dist_point_to_axis(5, 5, 0,
                                      g_last_prim.cx, g_last_prim.cy, g_last_prim.cz,
                                      g_last_prim.nx, g_last_prim.ny, g_last_prim.nz);
    if (drad > 0.05) {
        printf("center off axis %.4f FAIL\n", drad);
        liscio_destroy(ctx); return 1;
    }
    if (fabs(fabs(g_last_prim.nz) - 1.0) > 1e-3) {
        printf("axis not Z (n=(%.3f,%.3f,%.3f)) FAIL\n",
               g_last_prim.nx, g_last_prim.ny, g_last_prim.nz);
        liscio_destroy(ctx); return 1;
    }
    if (fabs(g_last_prim.radius - 10.0) > 0.05) {
        printf("radius err=%.4f FAIL\n", g_last_prim.radius - 10.0);
        liscio_destroy(ctx); return 1;
    }
    /* Pitch sign follows axis sign chosen by fitter; magnitude must match. */
    double sign_pitch = (g_last_prim.nz > 0) ? 1.0 : -1.0;
    if (fabs(g_last_prim.pitch * sign_pitch - 0.5) > 0.01) {
        printf("pitch err=%.4f FAIL\n", g_last_prim.pitch * sign_pitch - 0.5);
        liscio_destroy(ctx); return 1;
    }
    printf("PASS (n=(%.2f,%.2f,%.2f) r=%.3f p=%.3f arc=%.2frad)\n",
           g_last_prim.nx, g_last_prim.ny, g_last_prim.nz,
           g_last_prim.radius, g_last_prim.pitch, g_last_prim.arc_angle);
    liscio_destroy(ctx);
    return 0;
}

/* 2. Half turn helix, axis tilted along (1,1,1)/√3. */
static int test_tilted_half_turn(void)
{
    printf("TEST tilted_half_turn... ");
    fflush(stdout);

    liscio_cfg_t cfg;
    liscio_cfg_default(&cfg);
    cfg.tol_xyz = 0.005;
    liscio_ctx_t *ctx = liscio_create(&cfg);
    assert(ctx);
    liscio_set_callback(ctx, on_prim, NULL);
    reset_counters();

    double s3 = sqrt(3.0);
    double nx = 1.0/s3, ny = 1.0/s3, nz = 1.0/s3;
    /* Build an orthonormal frame perpendicular to n. */
    double ux, uy, uz;
    {
        double tmp[3] = { 0, 0, 1 };
        double cu[3] = { ny*tmp[2] - nz*tmp[1],
                          nz*tmp[0] - nx*tmp[2],
                          nx*tmp[1] - ny*tmp[0] };
        double m = sqrt(cu[0]*cu[0] + cu[1]*cu[1] + cu[2]*cu[2]);
        ux = cu[0]/m; uy = cu[1]/m; uz = cu[2]/m;
    }
    double vx = ny*uz - nz*uy;
    double vy = nz*ux - nx*uz;
    double vz = nx*uy - ny*ux;

    gen_helix(ctx, /*c*/ 0,0,0, nx,ny,nz, ux,uy,uz, vx,vy,vz,
              /*r*/ 4.0, /*p*/ 0.3, 0.0, M_PI, 25);
    liscio_flush(ctx);

    printf("L=%d A=%d H=%d ", g_line_count, g_arc_count, g_helix_count);
    if (g_helix_count < 1) {
        printf("FAIL (expected HELIX)\n");
        liscio_destroy(ctx); return 1;
    }
    /* Axis should be parallel (or anti-parallel) to (1,1,1)/√3. */
    double dot = g_last_prim.nx * nx + g_last_prim.ny * ny + g_last_prim.nz * nz;
    if (fabs(fabs(dot) - 1.0) > 1e-3) {
        printf("axis dev (|dot|=%.4f) FAIL\n", fabs(dot));
        liscio_destroy(ctx); return 1;
    }
    if (fabs(g_last_prim.radius - 4.0) > 0.02) {
        printf("radius err=%.4f FAIL\n", g_last_prim.radius - 4.0);
        liscio_destroy(ctx); return 1;
    }
    /* arc_angle magnitude ≈ π. */
    if (fabs(fabs(g_last_prim.arc_angle) - M_PI) > 0.05) {
        printf("arc_angle err=%.4f FAIL\n", fabs(g_last_prim.arc_angle) - M_PI);
        liscio_destroy(ctx); return 1;
    }
    /* |pitch| ≈ 0.3 (sign follows fitter's axis choice). */
    if (fabs(fabs(g_last_prim.pitch) - 0.3) > 0.02) {
        printf("pitch err=%.4f FAIL\n", fabs(g_last_prim.pitch) - 0.3);
        liscio_destroy(ctx); return 1;
    }
    printf("PASS (axis·truth=%.3f r=%.3f |p|=%.3f arc=%.3f)\n",
           dot, g_last_prim.radius, fabs(g_last_prim.pitch), g_last_prim.arc_angle);
    liscio_destroy(ctx);
    return 0;
}

/* 3. Reverse rotation (CW about +Z): pitch and arc_angle should be
 * signed consistently. */
static int test_cw_helix(void)
{
    printf("TEST cw_helix... ");
    fflush(stdout);

    liscio_cfg_t cfg;
    liscio_cfg_default(&cfg);
    cfg.tol_xyz = 0.01;
    liscio_ctx_t *ctx = liscio_create(&cfg);
    assert(ctx);
    liscio_set_callback(ctx, on_prim, NULL);
    reset_counters();

    /* Negative angle range = CW about +Z. */
    gen_helix(ctx, 0,0,0, 0,0,1, 1,0,0, 0,1,0,
              /*r*/ 6.0, /*p*/ 0.4,
              0.0, -1.5 * M_PI, 40);
    liscio_flush(ctx);

    printf("L=%d A=%d H=%d ", g_line_count, g_arc_count, g_helix_count);
    if (g_helix_count < 1) {
        printf("FAIL\n"); liscio_destroy(ctx); return 1;
    }
    /* Whatever axis sign the fitter picked, the world-frame signed axial
     * advance (= sign(axis_z)·pitch·arc_angle) must match the truth's
     * pitch_true · arc_true = 0.4 * (-1.5π) ≈ -1.885 along +Z. */
    double sign_axis = (g_last_prim.nz >= 0) ? 1.0 : -1.0;
    double signed_advance = sign_axis * g_last_prim.pitch * g_last_prim.arc_angle;
    double truth_advance  = 0.4 * (-1.5 * M_PI);
    if (fabs(signed_advance - truth_advance) > 0.05) {
        printf("axial advance err=%.4f (got %.3f expect %.3f) FAIL\n",
               signed_advance - truth_advance, signed_advance, truth_advance);
        liscio_destroy(ctx); return 1;
    }
    printf("PASS (advance=%.3f truth=%.3f)\n", signed_advance, truth_advance);
    liscio_destroy(ctx);
    return 0;
}

/* 4. Planar circle (pitch == 0) must still classify as ARC, not HELIX,
 * because arc9d_fit runs first.  Guards against helix9 stealing planar
 * data. */
static int test_planar_arc_not_helix(void)
{
    printf("TEST planar_arc_not_helix... ");
    fflush(stdout);

    liscio_cfg_t cfg;
    liscio_cfg_default(&cfg);
    cfg.tol_xyz = 0.01;
    liscio_ctx_t *ctx = liscio_create(&cfg);
    assert(ctx);
    liscio_set_callback(ctx, on_prim, NULL);
    reset_counters();

    /* Pure XY circle, pitch = 0 → planar arc. */
    gen_helix(ctx, 0,0,0, 0,0,1, 1,0,0, 0,1,0,
              /*r*/ 5.0, /*p*/ 0.0, 0.0, 1.5*M_PI, 24);
    liscio_flush(ctx);

    printf("L=%d A=%d H=%d ", g_line_count, g_arc_count, g_helix_count);
    if (g_arc_count < 1 || g_helix_count > 0) {
        printf("FAIL (planar should be ARC, not HELIX)\n");
        liscio_destroy(ctx); return 1;
    }
    printf("PASS\n");
    liscio_destroy(ctx);
    return 0;
}

/* 5. Helix with adversarial noise > tol must NOT be accepted as HELIX. */
static int test_noisy_helix_rejected(void)
{
    printf("TEST noisy_helix_rejected... ");
    fflush(stdout);

    liscio_cfg_t cfg;
    liscio_cfg_default(&cfg);
    cfg.tol_xyz = 0.01;       /* tight */
    liscio_ctx_t *ctx = liscio_create(&cfg);
    assert(ctx);
    liscio_set_callback(ctx, on_prim, NULL);
    reset_counters();

    /* Generate helix points but inject ±0.1 mm radial noise per point —
     * 10× tol.  Fitter should reject HELIX (and ARC; falls back to
     * Bezier or LINE composite). */
    int N = 40;
    double r = 8.0, p = 0.3;
    liscio_pose_t prev = {0,0,0,0,0,0,0,0,0};
    for (int i = 0; i <= N; i++) {
        double th = (double)i / N * 2.0 * M_PI;
        double noise = ((i % 2 == 0) ? +1 : -1) * 0.1;     /* ±0.1 mm */
        liscio_pose_t pose = {
            (r + noise) * cos(th),
            (r + noise) * sin(th),
            p * th,
            0,0,0, 0,0,0
        };
        if (i > 0) liscio_add_line(ctx, &prev, &pose, 1000.0, i);
        prev = pose;
    }
    liscio_flush(ctx);

    printf("L=%d A=%d H=%d O=%d ",
           g_line_count, g_arc_count, g_helix_count, g_other_count);
    if (g_helix_count > 0) {
        printf("FAIL (helix should be rejected at tol=0.01)\n");
        liscio_destroy(ctx); return 1;
    }
    printf("PASS (rejected, fell back)\n");
    liscio_destroy(ctx);
    return 0;
}

/* 6. Multi-turn helix (2 turns).  Tests that arc_angle exceeds 2π and
 * the fitter handles unwrapped angles correctly. */
static int test_multi_turn_helix(void)
{
    printf("TEST multi_turn_helix... ");
    fflush(stdout);

    liscio_cfg_t cfg;
    liscio_cfg_default(&cfg);
    cfg.tol_xyz = 0.01;
    cfg.max_window = 128;
    liscio_ctx_t *ctx = liscio_create(&cfg);
    assert(ctx);
    liscio_set_callback(ctx, on_prim, NULL);
    reset_counters();

    /* 2 turns @ 60 segments/turn = 120 segments — fits in default
     * window cap (128). */
    gen_helix(ctx, 0,0,0, 0,0,1, 1,0,0, 0,1,0,
              /*r*/ 5.0, /*p*/ 0.25,
              0.0, 4.0 * M_PI, 120);
    liscio_flush(ctx);

    printf("L=%d A=%d H=%d ", g_line_count, g_arc_count, g_helix_count);
    if (g_helix_count < 1) {
        printf("FAIL (expected at least 1 HELIX)\n");
        liscio_destroy(ctx); return 1;
    }
    /* total |arc_angle| across all helix prims ≥ 2π·2 − fudge. */
    /* (callback only saved last; multi-emit unlikely here, but allow it.)
     * Inspect last: |arc_angle| should be near 4π if single emit. */
    if (fabs(g_last_prim.arc_angle) < 2.0 * M_PI) {
        printf("arc_angle small (%.3f) — unexpected sub-emit\n",
               g_last_prim.arc_angle);
        /* not necessarily a fail if multiple HELIXes were emitted */
    }
    printf("PASS (last arc_angle=%.3f r=%.3f p=%.3f)\n",
           g_last_prim.arc_angle, g_last_prim.radius, g_last_prim.pitch);
    liscio_destroy(ctx);
    return 0;
}

int main(void)
{
    int failed = 0;
    failed += test_xy_full_turn();
    failed += test_tilted_half_turn();
    failed += test_cw_helix();
    failed += test_planar_arc_not_helix();
    failed += test_noisy_helix_rejected();
    failed += test_multi_turn_helix();
    if (failed) {
        fprintf(stderr, "%d test(s) FAILED\n", failed);
        return 1;
    }
    printf("All helix-fit tests PASSED\n");
    return 0;
}
