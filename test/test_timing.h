/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 杨阳 (Yang Yang) <mika-net@outlook.com>
 *
 * Tiny timing helper shared by liscio test binaries. */

#ifndef LISCIO_TEST_TIMING_H
#define LISCIO_TEST_TIMING_H

#define _POSIX_C_SOURCE 199309L
#include <time.h>

static inline double liscio_test_now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

#define LISCIO_TEST_TIMER_START(var)  double var = liscio_test_now_sec()
#define LISCIO_TEST_TIMER_MS(var)     ((liscio_test_now_sec() - (var)) * 1000.0)

#endif
