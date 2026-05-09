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
* Description: Throughput benchmark.  Feeds N random 9D waypoints
*              (noisy helical walk), measures wall-clock time for
*              the complete preprocess pipeline.  Reports
*              waypoints/sec and per-primitive absorb ratio.
* Author:      杨阳 (Yang Yang) <mika-net@outlook.com>
* License:     MIT (SPDX-License-Identifier: MIT)
* Copyright (c) 2026 杨阳 (Yang Yang)
********************************************************************/

#define _USE_MATH_DEFINES
#define _POSIX_C_SOURCE 199309L
#include "liscio/liscio.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static long g_emit = 0;
static long g_absorb = 0;
static void on_prim(const liscio_primitive_t *p, void *user)
{
    (void)user;
    g_emit++;
    g_absorb += p->n_absorbed;
}

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv)
{
    int N = (argc > 1) ? atoi(argv[1]) : 100000;
    int seed = (argc > 2) ? atoi(argv[2]) : 42;

    liscio_ctx_t *ctx = liscio_create(NULL);
    liscio_set_callback(ctx, on_prim, NULL);
    srand(seed);

    /* Generate smooth curve: helical mix. */
    double t0 = now_sec();
    liscio_pose_t prev = {0,0,0,0,0,0,0,0,0};
    for (int i = 0; i < N; i++) {
        double t = (double)i * 0.01;
        liscio_pose_t e = {
            10.0 * cos(t), 10.0 * sin(t), 0.05 * t,
            0.1 * t, 0, 0, 0, 0, 0
        };
        /* Add small noise. */
        e.x += ((rand() & 0xff) / 255.0 - 0.5) * 0.005;
        e.y += ((rand() & 0xff) / 255.0 - 0.5) * 0.005;
        liscio_add_line(ctx, &prev, &e, 1000.0, i);
        prev = e;
    }
    liscio_flush(ctx);
    double dt = now_sec() - t0;

    liscio_stats_t st;
    liscio_get_stats(ctx, &st);

    printf("BENCH %d waypoints in %.4f s = %.1f kpts/s\n",
        N, dt, N / dt / 1000.0);
    printf("  emit: line=%ld arc=%ld spline(bezier)=%ld\n",
        st.emitted_line, st.emitted_arc, st.emitted_spline);
    printf("  absorbed total: %ld (ratio %.2fx)\n",
        st.absorbed_total,
        (st.emitted_line + st.emitted_arc + st.emitted_spline) > 0
            ? (double)st.absorbed_total /
              (st.emitted_line + st.emitted_arc + st.emitted_spline)
            : 0.0);

    liscio_destroy(ctx);
    return 0;
}
