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
* Description: Rational cubic Bezier test suite.
*              Verifies arc → rational cubic Bezier conversion
*              produces a geometrically exact representation
*              (sample rational Bezier, compare against analytical
*              circle; error at double-precision floor ~1e-15).
*              Also round-trips rational evaluator on default weights
*              = non-rational Bezier.
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

static int test_quarter_circle_rbezier(void)
{
    printf("TEST arc_to_rbezier_quarter... "); fflush(stdout);

    /* Build an arc primitive: quarter circle in XY, r=10, center (0,0). */
    liscio_primitive_t arc;
    memset(&arc, 0, sizeof(arc));
    arc.type = LISCIO_PRIM_ARC;
    arc.start = (liscio_pose_t){ 10, 0, 0, 0,0,0, 0,0,0 };
    arc.end   = (liscio_pose_t){ 0, 10, 0, 0,0,0, 0,0,0 };
    arc.cx = 0; arc.cy = 0; arc.cz = 0;
    arc.nx = 0; arc.ny = 0; arc.nz = 1;  /* +Z normal = CCW in XY */
    arc.radius = 10;
    arc.arc_angle = M_PI / 2.0;  /* 90° CCW */
    /* Tangents: at (10,0) CCW → (0, 1); at (0,10) → (-1, 0). */
    arc.tan_start_x = 0;  arc.tan_start_y = 1;  arc.tan_start_z = 0;
    arc.tan_end_x   = -1; arc.tan_end_y   = 0;  arc.tan_end_z   = 0;

    liscio_primitive_t rb;
    if (liscio_arc_to_rational_bezier(&arc, &rb) != 0) {
        printf("FAIL (conversion)\n"); return 1;
    }

    /* Expected weights (degree-elevated from quadratic cos(θ/2)):
     *   w_mid = (1 + 2·cos(θ/2)) / 3  for θ = π/2 → (1 + √2) / 3 ≈ 0.805. */
    double ew = (1.0 + 2.0 * cos(M_PI / 4.0)) / 3.0;
    if (fabs(rb.w0 - 1.0) > 1e-9 || fabs(rb.w1 - ew) > 1e-6 ||
        fabs(rb.w2 - ew) > 1e-6 || fabs(rb.w3 - 1.0) > 1e-9) {
        printf("FAIL (weights %.4f %.4f %.4f %.4f expected middle %.4f)\n",
            rb.w0, rb.w1, rb.w2, rb.w3, ew); return 1;
    }

    /* Sample rational Bezier at N points, verify on circle r=10, θ∈[0,π/2]. */
    const int N = 40;
    double max_err = 0;
    for (int i = 0; i <= N; i++) {
        double t = (double)i / N;
        double x, y, z;
        liscio_rbezier_eval(&rb, t, &x, &y, &z);
        double r_actual = sqrt(x*x + y*y);
        double dr = fabs(r_actual - 10.0);
        if (dr > max_err) max_err = dr;
    }

    if (max_err > 1e-6) {
        printf("FAIL (max_err=%.2e)\n", max_err); return 1;
    }
    printf("PASS (max circle-err=%.2e, w=(%.4f,%.4f,%.4f,%.4f))\n",
        max_err, rb.w0, rb.w1, rb.w2, rb.w3);
    return 0;
}

static int test_rbezier_non_rational_equals_bezier(void)
{
    printf("TEST rbezier_default_weights... "); fflush(stdout);

    /* Create a simple cubic Bezier primitive with all weights = 1. */
    liscio_primitive_t bz;
    memset(&bz, 0, sizeof(bz));
    bz.type = LISCIO_PRIM_BEZIER;
    bz.P0 = (liscio_pose_t){ 0, 0, 0, 0,0,0, 0,0,0 };
    bz.P1 = (liscio_pose_t){ 10, 0, 0, 0,0,0, 0,0,0 };
    bz.P2 = (liscio_pose_t){ 0, 10, 0, 0,0,0, 0,0,0 };
    bz.P3 = (liscio_pose_t){ 10, 10, 0, 0,0,0, 0,0,0 };
    bz.w0 = bz.w1 = bz.w2 = bz.w3 = 1.0;

    /* At t=0.5 non-rational Bezier:
     *   P(0.5) = 0.125*P0 + 0.375*P1 + 0.375*P2 + 0.125*P3
     *          = 0.125*(0,0) + 0.375*(10,0) + 0.375*(0,10) + 0.125*(10,10)
     *          = (3.75 + 1.25, 3.75 + 1.25) = (5, 5) */
    double x, y, z;
    liscio_rbezier_eval(&bz, 0.5, &x, &y, &z);
    if (fabs(x - 5.0) > 1e-9 || fabs(y - 5.0) > 1e-9) {
        printf("FAIL (got %.4f,%.4f expected 5,5)\n", x, y); return 1;
    }
    printf("PASS (t=0.5 → (%.4f,%.4f))\n", x, y);
    return 0;
}

#include "test_timing.h"

int main(void)
{
    LISCIO_TEST_TIMER_START(t0);
    int fails = 0;
    fails += test_rbezier_non_rational_equals_bezier();
    fails += test_quarter_circle_rbezier();
    printf("---- liscio rbezier tests: %s (%d, time=%.3f ms) ----\n",
           fails ? "FAIL" : "OK", fails, LISCIO_TEST_TIMER_MS(t0));
    return fails ? 1 : 0;
}
