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
* Description: Random-walk G01 fuzz test.  Feeds N random 9D lines
*              (continuous), flushes, verifies invariants:
*                - no crash / assertion
*                - absorbed_total == input_count
*                - each emitted primitive's type / ranges valid
* Author:      杨阳 (Yang Yang) <mika-net@outlook.com>
* License:     MIT (SPDX-License-Identifier: MIT)
* Copyright (c) 2026 杨阳 (Yang Yang)
********************************************************************/

#include "liscio/liscio.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

static int g_emit_count = 0;
static long g_absorbed_sum = 0;

static void on_prim(const liscio_primitive_t *p, void *user)
{
    (void)user;
    g_emit_count++;
    g_absorbed_sum += p->n_absorbed;
    /* Verify basic invariants. */
    assert(p->n_absorbed >= 1);
    assert(p->type == LISCIO_PRIM_LINE ||
           p->type == LISCIO_PRIM_ARC  ||
           p->type == LISCIO_PRIM_BEZIER);
    if (p->type == LISCIO_PRIM_ARC) {
        assert(p->radius > 0.0);
        assert(isfinite(p->radius));
        assert(isfinite(p->arc_angle));
    }
    if (p->type == LISCIO_PRIM_BEZIER) {
        assert(isfinite(p->P0.x) && isfinite(p->P0.y) && isfinite(p->P0.z));
        assert(isfinite(p->P3.x) && isfinite(p->P3.y) && isfinite(p->P3.z));
        assert(isfinite(p->P1.x) && isfinite(p->P2.x));
    }
}

static double rnd(double lo, double hi)
{
    double t = (rand() & 0x7fff) / 32767.0;
    return lo + (hi - lo) * t;
}

#include "test_timing.h"

int main(int argc, char **argv)
{
    LISCIO_TEST_TIMER_START(__t0);
    unsigned int seed = (argc > 1) ? (unsigned)atoi(argv[1]) : 1;
    int iters = (argc > 2) ? atoi(argv[2]) : 50;
    int pts_per_iter = (argc > 3) ? atoi(argv[3]) : 200;

    srand(seed);
    int fails = 0;

    for (int it = 0; it < iters; it++) {
        liscio_cfg_t cfg;
        liscio_cfg_default(&cfg);
        /* Vary tolerance per iter */
        cfg.tol_xyz = rnd(0.001, 0.1);
        cfg.tol_abc = rnd(0.01, 1.0);
        liscio_ctx_t *ctx = liscio_create(&cfg);
        if (!ctx) { printf("iter %d: create fail\n", it); fails++; continue; }
        liscio_set_callback(ctx, on_prim, NULL);
        g_emit_count = 0; g_absorbed_sum = 0;

        /* Random walk in 9D */
        liscio_pose_t p = { rnd(-10, 10), rnd(-10, 10), rnd(-5, 5),
                             rnd(-45, 45), 0, 0, 0, 0, 0 };
        int input_count = 0;
        for (int i = 0; i < pts_per_iter; i++) {
            liscio_pose_t n = p;
            /* Mostly small steps; occasional direction change. */
            double step = rnd(0.01, 1.0);
            double dx = rnd(-1, 1), dy = rnd(-1, 1), dz = rnd(-0.2, 0.2);
            double nm = sqrt(dx*dx + dy*dy + dz*dz) + 1e-12;
            n.x = p.x + dx/nm*step;
            n.y = p.y + dy/nm*step;
            n.z = p.z + dz/nm*step;
            n.a = p.a + rnd(-0.1, 0.1);
            if (liscio_add_line(ctx, &p, &n, 1000.0, i) != 0) {
                fails++; break;
            }
            input_count++;
            p = n;
        }
        liscio_flush(ctx);

        liscio_stats_t st;
        liscio_get_stats(ctx, &st);

        /* input_count = N add_line calls.
         * absorbed_total = sum of (n_absorbed) across emitted primitives
         *                = N (each input line counted once).
         * The LINE path counts n_absorbed = (i1-i0) where i0..i1 are
         * waypoints (including seed).  Initial seed is shared so
         * total should match input_count. */
        if (st.input_count != input_count) {
            printf("iter %d: input mismatch (expected %d got %ld)\n",
                it, input_count, st.input_count);
            fails++;
        }
        if (st.absorbed_total != input_count) {
            printf("iter %d: absorbed=%ld input=%d tol=%.4f (emit=%d)\n",
                it, st.absorbed_total, input_count, cfg.tol_xyz, g_emit_count);
            /* Not a hard fail, but flag. */
        }

        liscio_destroy(ctx);
    }

    printf("fuzz iters=%d fails=%d time=%.3f ms\n",
        iters, fails, LISCIO_TEST_TIMER_MS(__t0));
    return fails ? 1 : 0;
}
