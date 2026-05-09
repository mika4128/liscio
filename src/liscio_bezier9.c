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
* Description: liscio 9D cubic Bezier fit.
*                Cubic Bezier: B(t) = b0·P0 + b1·P1 + b2·P2 + b3·P3,
*                              t ∈ [0,1]
*                  b0=(1-t)³, b1=3t(1-t)², b2=3t²(1-t), b3=t³
*                P0, P3 pinned to first/last waypoint.
*                Per-dimension 2x2 LSQ (decoupled; shared Bernstein
*                basis). 9D: xyz + abc + uvw each independent.
*                Chord-length parameterization, then Hoschek Newton
*                reparameterization (3 iter) for tight fit.
*                Composite recursion: if single-span exceeds tol,
*                split at mid-chord and fit each half (max depth 8).
* References:  L. Piegl & W. Tiller, *The NURBS Book*, 2nd ed.,
*              Springer 1997 — §9.4 (LSQ curve fitting).
*              J. Hoschek, *Intrinsic parametrization for
*              approximation*, CAGD 5(1), 1988.
* Author:      杨阳 (Yang Yang) <mika-net@outlook.com>
* License:     MIT (SPDX-License-Identifier: MIT)
* Copyright (c) 2026 杨阳 (Yang Yang)
********************************************************************/

#define _USE_MATH_DEFINES
#include "liscio/liscio.h"
#include "liscio_internal.h"

#include <stdio.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Debug counters: G1 successes vs unconstrained fallbacks.  Compiled
 * in only with -DLISCIO_DEBUG_G1; tests print under env LISCIO_G1_STATS=1. */
#ifdef LISCIO_DEBUG_G1
long liscio_dbg_g1_ok   = 0;
long liscio_dbg_g1_fall = 0;
long liscio_dbg_g1_fast_ok   = 0;
long liscio_dbg_g1_fast_fall = 0;
#define LISCIO_DBG_G1_INC(x) ((x)++)
#else
#define LISCIO_DBG_G1_INC(x) ((void)0)
#endif

/* Bernstein basis. */
static inline double B0(double t) { double u = 1.0 - t; return u*u*u; }
static inline double B1(double t) { double u = 1.0 - t; return 3.0*t*u*u; }
static inline double B2(double t) { double u = 1.0 - t; return 3.0*t*t*u; }
static inline double B3(double t) { return t*t*t; }

/* Solve 2x2 linear system  [a b; c d] x = [e; f]  via Cramer's rule.
 * Returns 0 on success, -1 on singular. */
static int solve2x2(double a, double b, double c, double d,
                    double e, double f,
                    double *x1, double *x2)
{
    double det = a*d - b*c;
    if (fabs(det) < 1e-20) return -1;
    *x1 = (e*d - b*f) / det;
    *x2 = (a*f - e*c) / det;
    return 0;
}

/* LSQ fit cubic Bezier 1D: given P0, P3, and N points (xs[i] at t[i]),
 * find P1, P2 minimizing sum(residual²).  Returns 0 on success.
 * If dimension is nearly constant (variance < 1e-18), sets
 * P1 = P0 + (P3-P0)/3, P2 = P0 + 2*(P3-P0)/3 (linear interpolation
 * Bezier control points). */
static int bezier_lsq_1d(const double *xs, const double *ts, int n,
                         double P0, double P3,
                         double *P1_out, double *P2_out)
{
    if (n < 4) return -1;

    /* Degenerate (near-constant) dimension: normal equations are singular
     * because all residuals are ≈ 0.  Use linear interpolation CPs. */
    double mean = 0;
    for (int i = 0; i < n; i++) mean += xs[i];
    mean /= n;
    double var = 0;
    for (int i = 0; i < n; i++) {
        double d = xs[i] - mean;
        var += d*d;
    }
    if (var < 1e-18) {
        *P1_out = P0 + (P3 - P0) / 3.0;
        *P2_out = P0 + 2.0 * (P3 - P0) / 3.0;
        return 0;
    }

    double sum_b1b1 = 0.0, sum_b1b2 = 0.0, sum_b2b2 = 0.0;
    double sum_b1c = 0.0, sum_b2c = 0.0;
    for (int i = 0; i < n; i++) {
        double t = ts[i];
        double b1 = B1(t), b2 = B2(t);
        double c  = xs[i] - B0(t)*P0 - B3(t)*P3;
        sum_b1b1 += b1*b1;
        sum_b1b2 += b1*b2;
        sum_b2b2 += b2*b2;
        sum_b1c  += b1*c;
        sum_b2c  += b2*c;
    }
    return solve2x2(sum_b1b1, sum_b1b2, sum_b1b2, sum_b2b2,
                    sum_b1c, sum_b2c, P1_out, P2_out);
}

/* Evaluate Bezier at t: B(t) = b0·P0 + b1·P1 + b2·P2 + b3·P3. */
static inline double bezier_eval(double P0, double P1, double P2, double P3,
                                  double t)
{
    return B0(t)*P0 + B1(t)*P1 + B2(t)*P2 + B3(t)*P3;
}

/* First derivative: C'(t) = 3(1-t)²(P1-P0) + 6t(1-t)(P2-P1) + 3t²(P3-P2). */
static inline double bezier_eval_d1(double P0, double P1, double P2, double P3,
                                      double t)
{
    double u = 1.0 - t;
    return 3.0*u*u*(P1-P0) + 6.0*t*u*(P2-P1) + 3.0*t*t*(P3-P2);
}

/* Second derivative: C''(t) = 6(1-t)(P2-2P1+P0) + 6t(P3-2P2+P1). */
static inline double bezier_eval_d2(double P0, double P1, double P2, double P3,
                                      double t)
{
    double u = 1.0 - t;
    return 6.0*u*(P2 - 2.0*P1 + P0) + 6.0*t*(P3 - 2.0*P2 + P1);
}

/* Hoschek single Newton iteration on t_i for point p_i.
 * Minimizes |C(t) - p|² using one Newton step. Returns refined t ∈ [0,1]. */
static double hoschek_reparam(double p_x, double p_y, double p_z,
                              double P0x, double P1x, double P2x, double P3x,
                              double P0y, double P1y, double P2y, double P3y,
                              double P0z, double P1z, double P2z, double P3z,
                              double t)
{
    double cx = bezier_eval(P0x, P1x, P2x, P3x, t);
    double cy = bezier_eval(P0y, P1y, P2y, P3y, t);
    double cz = bezier_eval(P0z, P1z, P2z, P3z, t);

    double d1x = bezier_eval_d1(P0x, P1x, P2x, P3x, t);
    double d1y = bezier_eval_d1(P0y, P1y, P2y, P3y, t);
    double d1z = bezier_eval_d1(P0z, P1z, P2z, P3z, t);

    double d2x = bezier_eval_d2(P0x, P1x, P2x, P3x, t);
    double d2y = bezier_eval_d2(P0y, P1y, P2y, P3y, t);
    double d2z = bezier_eval_d2(P0z, P1z, P2z, P3z, t);

    double diff_x = cx - p_x, diff_y = cy - p_y, diff_z = cz - p_z;
    double num = diff_x*d1x + diff_y*d1y + diff_z*d1z;
    double den = d1x*d1x + d1y*d1y + d1z*d1z
               + diff_x*d2x + diff_y*d2y + diff_z*d2z;
    if (fabs(den) < 1e-15) return t;

    double t_new = t - num / den;
    if (t_new < 0.0) t_new = 0.0;
    if (t_new > 1.0) t_new = 1.0;
    return t_new;
}

/* ---------- Public entry ---------- */
int liscio_bezier9_fit(const struct liscio_ctx *ctx, int i0, int i1,
                        liscio_bezier9_fit_t *out)
{
    if (!ctx || !out) return -1;
    int n = i1 - i0 + 1;
    if (n < 4) return -1;   /* need >=4 pts for cubic */
    if (n > LISCIO_MAX_WINDOW) return -1;

    const liscio_pose_t *pts = &ctx->pts[i0];

    /* Chord-length parameterization. */
    double ts[LISCIO_MAX_WINDOW];
    double cum_len[LISCIO_MAX_WINDOW];
    cum_len[0] = 0.0;
    for (int i = 1; i < n; i++) {
        double dx = pts[i].x - pts[i-1].x;
        double dy = pts[i].y - pts[i-1].y;
        double dz = pts[i].z - pts[i-1].z;
        cum_len[i] = cum_len[i-1] + sqrt(dx*dx + dy*dy + dz*dz);
    }
    double total = cum_len[n-1];
    if (total < 1e-12) return -1;
    for (int i = 0; i < n; i++) ts[i] = cum_len[i] / total;

    /* Fit each of 9 dimensions independently with shared chord-length ts. */
    double coords[9][LISCIO_MAX_WINDOW];
    double P0c[9], P3c[9], P1c[9], P2c[9];
    for (int i = 0; i < n; i++) {
        coords[0][i] = pts[i].x; coords[1][i] = pts[i].y; coords[2][i] = pts[i].z;
        coords[3][i] = pts[i].a; coords[4][i] = pts[i].b; coords[5][i] = pts[i].c;
        coords[6][i] = pts[i].u; coords[7][i] = pts[i].v; coords[8][i] = pts[i].w;
    }
    for (int d = 0; d < 9; d++) {
        P0c[d] = coords[d][0];
        P3c[d] = coords[d][n-1];
        if (bezier_lsq_1d(coords[d], ts, n, P0c[d], P3c[d], &P1c[d], &P2c[d]) != 0)
            return -1;
    }
    /* Legacy xyz shortcuts for Hoschek below. */
    double P0x=P0c[0], P3x=P3c[0], p1x=P1c[0], p2x=P2c[0];
    double P0y=P0c[1], P3y=P3c[1], p1y=P1c[1], p2y=P2c[1];
    double P0z=P0c[2], P3z=P3c[2], p1z=P1c[2], p2z=P2c[2];
    double *xs = coords[0], *ys = coords[1], *zs = coords[2];

    /* Hoschek iterative reparameterization: refine t_i toward closest-
     * point-on-curve, refit, repeat.  Stops when max Δt < eps or
     * max_iter reached.  Industry standard for tight-tolerance fits. */
    const int N_HOSCHEK_ITERS = 3;
    double max_dev = 0.0;
    for (int iter = 0; iter < N_HOSCHEK_ITERS; iter++) {
        /* Evaluate current fit; find max deviation. */
        max_dev = 0.0;
        for (int i = 0; i < n; i++) {
            double bx = bezier_eval(P0x, p1x, p2x, P3x, ts[i]);
            double by = bezier_eval(P0y, p1y, p2y, P3y, ts[i]);
            double bz = bezier_eval(P0z, p1z, p2z, P3z, ts[i]);
            double dx = bx - xs[i], dy = by - ys[i], dz = bz - zs[i];
            double d = sqrt(dx*dx + dy*dy + dz*dz);
            if (d > max_dev) max_dev = d;
        }
        if (max_dev <= ctx->cfg.tol_xyz) break;

        /* Newton step on each t_i (fix endpoints). */
        double max_dt = 0.0;
        for (int i = 1; i < n - 1; i++) {
            double t_new = hoschek_reparam(xs[i], ys[i], zs[i],
                P0x, p1x, p2x, P3x,
                P0y, p1y, p2y, P3y,
                P0z, p1z, p2z, P3z,
                ts[i]);
            double dt = fabs(t_new - ts[i]);
            if (dt > max_dt) max_dt = dt;
            ts[i] = t_new;
        }
        if (max_dt < 1e-8) break;

        /* Ensure ts monotone (safety; Newton may locally overshoot). */
        for (int i = 1; i < n; i++) {
            if (ts[i] <= ts[i-1]) ts[i] = ts[i-1] + 1e-9;
        }

        /* Refit with refined ts. */
        if (bezier_lsq_1d(xs, ts, n, P0x, P3x, &p1x, &p2x) != 0) return -1;
        if (bezier_lsq_1d(ys, ts, n, P0y, P3y, &p1y, &p2y) != 0) return -1;
        if (bezier_lsq_1d(zs, ts, n, P0z, P3z, &p1z, &p2z) != 0) return -1;
    }
    if (max_dev > ctx->cfg.tol_xyz) return -1;

    /* After Hoschek, xyz refined.  Refit rotary/uvw with updated ts so
     * their Bezier coeffs match the refined param.  Then verify each
     * within its own tolerance. */
    for (int d = 3; d < 9; d++) {
        if (bezier_lsq_1d(coords[d], ts, n, P0c[d], P3c[d],
                          &P1c[d], &P2c[d]) != 0) return -1;
    }
    for (int d = 3; d < 9; d++) {
        double tol = (d < 6) ? ctx->cfg.tol_abc : ctx->cfg.tol_uvw;
        for (int i = 0; i < n; i++) {
            double bv = bezier_eval(P0c[d], P1c[d], P2c[d], P3c[d], ts[i]);
            if (fabs(bv - coords[d][i]) > tol) return -1;
        }
    }

    /* Fill output: 9D control points. */
    out->P0 = (liscio_pose_t){ P0c[0], P0c[1], P0c[2],
                                P0c[3], P0c[4], P0c[5],
                                P0c[6], P0c[7], P0c[8] };
    out->P1 = (liscio_pose_t){ p1x, p1y, p1z,
                                P1c[3], P1c[4], P1c[5],
                                P1c[6], P1c[7], P1c[8] };
    out->P2 = (liscio_pose_t){ p2x, p2y, p2z,
                                P2c[3], P2c[4], P2c[5],
                                P2c[6], P2c[7], P2c[8] };
    out->P3 = (liscio_pose_t){ P3c[0], P3c[1], P3c[2],
                                P3c[3], P3c[4], P3c[5],
                                P3c[6], P3c[7], P3c[8] };
    out->total_chord = total;
    out->max_deviation = max_dev;
    return 0;
}

/* ---------- G1-constrained 9D fit (Schneider-style αL/αR LSQ) ----- */

static int bezier_lsq_1d_pinned(const double *xs, const double *ts, int n,
                                double P0, double P3,
                                double *P1_out, double *P2_out)
{
    return bezier_lsq_1d(xs, ts, n, P0, P3, P1_out, P2_out);
}

int liscio_bezier9_fit_g1(const struct liscio_ctx *ctx, int i0, int i1,
                          double t1x, double t1y, double t1z,
                          double t2x, double t2y, double t2z,
                          liscio_bezier9_fit_t *out)
{
    if (!ctx || !out) return -1;
    int n = i1 - i0 + 1;
    if (n < 4) return -1;
    if (n > LISCIO_MAX_WINDOW) return -1;

    /* Normalize tHats defensively. */
    double n1 = sqrt(t1x*t1x + t1y*t1y + t1z*t1z);
    double n2 = sqrt(t2x*t2x + t2y*t2y + t2z*t2z);
    if (n1 < 1e-15 || n2 < 1e-15) return -1;
    t1x /= n1; t1y /= n1; t1z /= n1;
    t2x /= n2; t2y /= n2; t2z /= n2;

    const liscio_pose_t *pts = &ctx->pts[i0];

    /* Chord-length parameterization. */
    double ts[LISCIO_MAX_WINDOW];
    double cum_len[LISCIO_MAX_WINDOW];
    cum_len[0] = 0.0;
    for (int i = 1; i < n; i++) {
        double dx = pts[i].x - pts[i-1].x;
        double dy = pts[i].y - pts[i-1].y;
        double dz = pts[i].z - pts[i-1].z;
        cum_len[i] = cum_len[i-1] + sqrt(dx*dx + dy*dy + dz*dz);
    }
    double seg_total = cum_len[n-1];
    if (seg_total < 1e-12) return -1;
    for (int i = 0; i < n; i++) ts[i] = cum_len[i] / seg_total;

    double P0x = pts[0].x,   P0y = pts[0].y,   P0z = pts[0].z;
    double P3x = pts[n-1].x, P3y = pts[n-1].y, P3z = pts[n-1].z;
    double xs[LISCIO_MAX_WINDOW], ys[LISCIO_MAX_WINDOW], zs[LISCIO_MAX_WINDOW];
    for (int i = 0; i < n; i++) {
        xs[i] = pts[i].x; ys[i] = pts[i].y; zs[i] = pts[i].z;
    }

    double t12 = t1x*t2x + t1y*t2y + t1z*t2z;  /* tHat1 · tHat2 */

    double alpha_l = 0.0, alpha_r = 0.0;
    double p1x = 0, p1y = 0, p1z = 0;
    double p2x = 0, p2y = 0, p2z = 0;
    double max_dev = 0.0;

    const int N_HOSCHEK_ITERS = 4;
    int converged = 0;
    for (int iter = 0; iter < N_HOSCHEK_ITERS; iter++) {
        /* 2x2 constrained-LSQ on αL, αR. */
        double sAA = 0, sBB = 0, sBcoeff = 0;
        double sAr = 0, sBr = 0;
        for (int i = 0; i < n; i++) {
            double t = ts[i];
            double b0 = B0(t), b1 = B1(t), b2 = B2(t), b3 = B3(t);
            sAA     += b1 * b1;
            sBB     += b2 * b2;
            sBcoeff += b1 * b2;
            /* r_i = q_i − (b0+b1)P0 − (b2+b3)P3 */
            double rx = xs[i] - (b0+b1)*P0x - (b2+b3)*P3x;
            double ry = ys[i] - (b0+b1)*P0y - (b2+b3)*P3y;
            double rz = zs[i] - (b0+b1)*P0z - (b2+b3)*P3z;
            sAr += b1 * (t1x*rx + t1y*ry + t1z*rz);
            sBr += b2 * (t2x*rx + t2y*ry + t2z*rz);
        }
        double sAB = sBcoeff * t12;
        if (solve2x2(sAA, sAB, sAB, sBB, sAr, sBr, &alpha_l, &alpha_r) != 0)
            return -1;

        /* Wu/Barsky fallback (Schneider 1990): when LSQ produces non-
         * positive α (P1 behind P0 or P2 past P3), clamp to seg_len/3.
         * Tangent direction (= tHat1, tHat2) is preserved either way, so
         * G1 continuity at the join is retained.  Caller-visible fit
         * quality is then validated below by the tol_xyz check. */
        double eps = 1e-6 * seg_total;
        if (alpha_l < eps) alpha_l = seg_total / 3.0;
        if (alpha_r < eps) alpha_r = seg_total / 3.0;

        p1x = P0x + alpha_l * t1x;
        p1y = P0y + alpha_l * t1y;
        p1z = P0z + alpha_l * t1z;
        p2x = P3x + alpha_r * t2x;
        p2y = P3y + alpha_r * t2y;
        p2z = P3z + alpha_r * t2z;

        /* Evaluate fit; track max XYZ deviation. */
        max_dev = 0.0;
        for (int i = 0; i < n; i++) {
            double bx = bezier_eval(P0x, p1x, p2x, P3x, ts[i]);
            double by = bezier_eval(P0y, p1y, p2y, P3y, ts[i]);
            double bz = bezier_eval(P0z, p1z, p2z, P3z, ts[i]);
            double dx = bx - xs[i], dy = by - ys[i], dz = bz - zs[i];
            double d = sqrt(dx*dx + dy*dy + dz*dz);
            if (d > max_dev) max_dev = d;
        }
        if (max_dev <= ctx->cfg.tol_xyz) { converged = 1; break; }

        /* Newton step on each interior t_i (Hoschek). */
        double max_dt = 0.0;
        for (int i = 1; i < n - 1; i++) {
            double t_new = hoschek_reparam(xs[i], ys[i], zs[i],
                P0x, p1x, p2x, P3x,
                P0y, p1y, p2y, P3y,
                P0z, p1z, p2z, P3z,
                ts[i]);
            double dt = fabs(t_new - ts[i]);
            if (dt > max_dt) max_dt = dt;
            ts[i] = t_new;
        }
        if (max_dt < 1e-8) break;
        for (int i = 1; i < n; i++)
            if (ts[i] <= ts[i-1]) ts[i] = ts[i-1] + 1e-9;
    }
    if (!converged && max_dev > ctx->cfg.tol_xyz) return -1;

    /* ABC/UVW free per-axis fit with the refined ts. */
    double P1abc[6], P2abc[6];
    double abc[6][LISCIO_MAX_WINDOW];
    for (int i = 0; i < n; i++) {
        abc[0][i] = pts[i].a; abc[1][i] = pts[i].b; abc[2][i] = pts[i].c;
        abc[3][i] = pts[i].u; abc[4][i] = pts[i].v; abc[5][i] = pts[i].w;
    }
    double P0abc[6] = { pts[0].a, pts[0].b, pts[0].c,
                        pts[0].u, pts[0].v, pts[0].w };
    double P3abc[6] = { pts[n-1].a, pts[n-1].b, pts[n-1].c,
                        pts[n-1].u, pts[n-1].v, pts[n-1].w };
    for (int d = 0; d < 6; d++) {
        if (bezier_lsq_1d_pinned(abc[d], ts, n, P0abc[d], P3abc[d],
                                 &P1abc[d], &P2abc[d]) != 0) return -1;
    }
    for (int d = 0; d < 6; d++) {
        double tol = (d < 3) ? ctx->cfg.tol_abc : ctx->cfg.tol_uvw;
        for (int i = 0; i < n; i++) {
            double bv = bezier_eval(P0abc[d], P1abc[d], P2abc[d], P3abc[d], ts[i]);
            if (fabs(bv - abc[d][i]) > tol) return -1;
        }
    }

    /* Pack 9D output. */
    out->P0 = pts[0];
    out->P3 = pts[n-1];
    out->P1 = (liscio_pose_t){ p1x, p1y, p1z,
                                P1abc[0], P1abc[1], P1abc[2],
                                P1abc[3], P1abc[4], P1abc[5] };
    out->P2 = (liscio_pose_t){ p2x, p2y, p2z,
                                P2abc[0], P2abc[1], P2abc[2],
                                P2abc[3], P2abc[4], P2abc[5] };
    out->total_chord = seg_total;
    out->max_deviation = max_dev;
    return 0;
}

/* ---------- Composite Bezier: recursive split until each segment fits.
 *
 * If the single-Bezier fit over [i0..i1] exceeds tol, find the index of
 * maximum deviation and recursively fit the two halves.  Emit each
 * successful sub-fit via callback.  Terminates when sub-range has
 * fewer than 4 points (falls back to LINE).
 *
 * Returns total number of emitted Bezier primitives (>=0), or -1 on
 * allocation or config error.  emit_cb is called with user-data for
 * each emitted sub-fit.
 */
/* Emit one linear cubic bezier (P1 = P0+(P3-P0)/3, P2 = 2/3 way) that
 * exactly coincides with the chord P0→P3.  Helper for per-segment
 * fallback emission. */
static void emit_single_linear_bezier(const struct liscio_ctx *ctx,
                                      int i0, int i1,
                                      void (*emit)(int i0, int i1,
                                                   const liscio_bezier9_fit_t *,
                                                   void *user),
                                      void *user)
{
    liscio_bezier9_fit_t fit;
    memset(&fit, 0, sizeof(fit));
    fit.P0 = ctx->pts[i0];
    fit.P3 = ctx->pts[i1];
#define LERP(field, alpha) \
    fit.P1.field = fit.P0.field + (alpha)*(fit.P3.field - fit.P0.field)
    LERP(x, 1.0/3.0); LERP(y, 1.0/3.0); LERP(z, 1.0/3.0);
    LERP(a, 1.0/3.0); LERP(b, 1.0/3.0); LERP(c, 1.0/3.0);
    LERP(u, 1.0/3.0); LERP(v, 1.0/3.0); LERP(w, 1.0/3.0);
#undef LERP
#define LERP2(field) \
    fit.P2.field = fit.P0.field + (2.0/3.0)*(fit.P3.field - fit.P0.field)
    LERP2(x); LERP2(y); LERP2(z);
    LERP2(a); LERP2(b); LERP2(c);
    LERP2(u); LERP2(v); LERP2(w);
#undef LERP2
    fit.max_deviation = 0.0;
    emit(i0, i1, &fit, user);
}

/* Emit the range [i0..i1] as one linear bezier PER SEGMENT, so the
 * emitted polyline exactly traces the waypoints (no chord shortcut). */
static void emit_polyline_fallback(const struct liscio_ctx *ctx,
                                   int i0, int i1,
                                   void (*emit)(int i0, int i1,
                                                const liscio_bezier9_fit_t *,
                                                void *user),
                                   void *user)
{
    for (int k = i0; k < i1; k++) {
        emit_single_linear_bezier(ctx, k, k + 1, emit, user);
    }
}

/* Compute chord-direction unit vector pts[a] → pts[b] (returns 0 vec if
 * coincident). */
static void chord_unit(const liscio_pose_t *pts, int a, int b,
                       double *ux, double *uy, double *uz)
{
    double dx = pts[b].x - pts[a].x;
    double dy = pts[b].y - pts[a].y;
    double dz = pts[b].z - pts[a].z;
    double n  = sqrt(dx*dx + dy*dy + dz*dz);
    if (n < 1e-15) { *ux=0; *uy=0; *uz=0; return; }
    *ux = dx/n; *uy = dy/n; *uz = dz/n;
}

/* Forward unit tangent at interior split point s — average of incoming
 * and outgoing chord directions, normalized.  G1 join uses ±this vector
 * for left/right sub-curves. */
static void center_tangent_forward(const liscio_pose_t *pts, int s,
                                   double *cx, double *cy, double *cz)
{
    double ux, uy, uz, vx, vy, vz;
    chord_unit(pts, s-1, s,   &ux, &uy, &uz);   /* incoming */
    chord_unit(pts, s,   s+1, &vx, &vy, &vz);   /* outgoing */
    double sx = ux + vx, sy = uy + vy, sz = uz + vz;
    double n = sqrt(sx*sx + sy*sy + sz*sz);
    if (n < 1e-15) {  /* opposing dirs (cusp) — fall back to outgoing */
        *cx = vx; *cy = vy; *cz = vz; return;
    }
    *cx = sx/n; *cy = sy/n; *cz = sz/n;
}

static int fit_recursive(const struct liscio_ctx *ctx,
                         int i0, int i1,
                         double t1x, double t1y, double t1z,
                         double t2x, double t2y, double t2z,
                         void (*emit)(int i0, int i1,
                                      const liscio_bezier9_fit_t *,
                                      void *user),
                         void *user,
                         int depth)
{
    const int MAX_DEPTH = 8;  /* 2^8 = 256 splits max */
    int n = i1 - i0 + 1;
    if (n < 2) return 0;                /* degenerate */
    if (n < 4) {                        /* 2–3 points: emit as linear bezier */
        emit_polyline_fallback(ctx, i0, i1, emit, user);
        return 1;
    }
    if (depth > MAX_DEPTH) {            /* recursion cap: emit linear to avoid loss */
        emit_polyline_fallback(ctx, i0, i1, emit, user);
        return 1;
    }

    liscio_bezier9_fit_t fit;
    int rc = liscio_bezier9_fit_g1(ctx, i0, i1,
                                   t1x, t1y, t1z, t2x, t2y, t2z, &fit);
    if (rc == 0) {
        LISCIO_DBG_G1_INC(liscio_dbg_g1_ok);
        emit(i0, i1, &fit, user);
        return 1;
    }

    /* G1 fit missed tol.  Prefer to split (smaller halves more likely to
     * satisfy G1 + tol) over an unconstrained fit that breaks continuity.
     * Only fall back to unconstrained when the range is too small to
     * usefully split or recursion depth is exhausted — in those last-
     * resort cases, the C0 join is acceptable since it's bounded by
     * either the recursion cap or the 4-point Bezier minimum. */
    int can_split = (n >= 8) && (depth < MAX_DEPTH);
    if (!can_split) {
        rc = liscio_bezier9_fit(ctx, i0, i1, &fit);
        if (rc == 0) {
            LISCIO_DBG_G1_INC(liscio_dbg_g1_fall);
            emit(i0, i1, &fit, user);
            return 1;
        }
    }

    /* Fit failed — find the point of maximum deviation and split. */
    const liscio_pose_t *pts = &ctx->pts[i0];

    /* Recompute ts to find max-dev index. Cheap redo. */
    double cum[LISCIO_MAX_WINDOW];
    cum[0] = 0.0;
    for (int i = 1; i < n; i++) {
        double dx = pts[i].x - pts[i-1].x;
        double dy = pts[i].y - pts[i-1].y;
        double dz = pts[i].z - pts[i-1].z;
        cum[i] = cum[i-1] + sqrt(dx*dx + dy*dy + dz*dz);
    }
    double total = cum[n-1];
    if (total < 1e-12) return 0;

    /* Split at midpoint by chord length (robust default). */
    int split = 1;
    for (int i = 1; i < n - 1; i++) {
        if (cum[i] >= total * 0.5) { split = i; break; }
    }
    if (split < 2) split = 2;
    if (split > n - 3) split = n - 3;
    if (split < 2 || split > n - 3) {
        /* Too short to split cubic/cubic: emit linear bezier covering all. */
        emit_polyline_fallback(ctx, i0, i1, emit, user);
        return 1;
    }

    /* G1-shared tangent at the split point: average chord dirs in/out of
     * pts[i0+split].  Left sub-curve receives -tHat (Schneider inward
     * convention at right end); right sub-curve receives +tHat (forward
     * at left end).  Both sub-fits then evaluate to the same tangent
     * direction at the join, eliminating composite-recursion tan_mismatch. */
    double cx, cy, cz;
    center_tangent_forward(&ctx->pts[0], i0 + split, &cx, &cy, &cz);

    int cnt = 0;
    cnt += fit_recursive(ctx, i0, i0 + split,
                         t1x, t1y, t1z, -cx, -cy, -cz,
                         emit, user, depth + 1);
    cnt += fit_recursive(ctx, i0 + split, i1,
                         cx, cy, cz, t2x, t2y, t2z,
                         emit, user, depth + 1);
    return cnt;
}

/* Wrapper: try a single fit; if it fails, attempt composite (recursive
 * split).  Calls back once per emitted sub-bezier.  Outer-end tangents
 * are derived from the first/last raw chord directions (Schneider
 * convention: tHat1 forward at start, tHat2 inward at end). */
int liscio_bezier9_composite_fit(const struct liscio_ctx *ctx, int i0, int i1,
                                  liscio_bezier9_emit_fn emit, void *user)
{
    return liscio_bezier9_composite_fit_g1(ctx, i0, i1, NULL, emit, user);
}

int liscio_bezier9_composite_fit_g1(const struct liscio_ctx *ctx,
                                     int i0, int i1,
                                     const double *left_tan_xyz,
                                     liscio_bezier9_emit_fn emit, void *user)
{
    if (!ctx || !emit) return -1;
    if (i1 - i0 < 1) return 0;
    const liscio_pose_t *pts = &ctx->pts[0];
    double t1x, t1y, t1z, t2x, t2y, t2z;
    if (left_tan_xyz) {
        t1x = left_tan_xyz[0]; t1y = left_tan_xyz[1]; t1z = left_tan_xyz[2];
    } else {
        chord_unit(pts, i0, i0+1, &t1x, &t1y, &t1z);
    }
    chord_unit(pts, i1, i1-1, &t2x, &t2y, &t2z);   /* inward at right */
    return fit_recursive(ctx, i0, i1,
                         t1x, t1y, t1z, t2x, t2y, t2z,
                         emit, user, 0);
}
