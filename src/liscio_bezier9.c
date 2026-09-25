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
#include <stdlib.h>   /* getenv: TP2_G2DBG diagnostic */
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Fit-input densification capacity.  The fit acceptance gate (max_dev <=
 * tol_xyz) samples the curve only AT the input points; a long straight
 * segment inside a window contributes just 2 waypoints, so an unchecked
 * bulge can grow inside its chord (3_1001 N4755 stay-down: 72 mm G1
 * between scallop micro-lines → 17.7 mm phantom bulge, smooth and
 * jerk-invisible).  densify_window() inserts collinear interior samples
 * on over-long spans so BOTH the LSQ objective and the gate see the
 * chord interior: the fit is then provably within tol everywhere (span
 * shorter than 8·tol bulges at most ~(8·tol/3)·sin15° ≈ 0.69·tol by the
 * control-polygon convex-hull bound, i.e. sub-tolerance even unsampled).
 * Pure-micro windows (all spans < step) get zero added points. */
#define LISCIO_FIT_DENSIFY_CAP 256

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

/* Copy src[0..n-1] into dst, inserting evenly spaced collinear interior
 * samples on every span whose chord exceeds
 *     step = max(8·tol_xyz, total_chord / (cap − n))
 * so the output never exceeds cap points and spans shorter than ~8·tol
 * (geometrically incapable of hiding a >tol bulge) stay untouched.
 * All 9 pose dims are linearly interpolated along the chord — exact
 * ground truth for a G1 segment.  First/last dst points are bit-copies
 * of the original endpoints (endpoint semantics unchanged for callers).
 * Returns the new count, or -1 on defensive misuse (caller falls back
 * to the original points). */
static int densify_window(const liscio_pose_t *src, int n, double tol_xyz,
                          liscio_pose_t *dst, int cap)
{
    if (!src || !dst || n < 2 || cap < n) return -1;

    double total = 0.0;
    for (int i = 1; i < n; i++) {
        double dx = src[i].x - src[i-1].x;
        double dy = src[i].y - src[i-1].y;
        double dz = src[i].z - src[i-1].z;
        total += sqrt(dx*dx + dy*dy + dz*dz);
    }
    if (total < 1e-12) {                 /* zero-XYZ (W-only) window */
        for (int i = 0; i < n; i++) dst[i] = src[i];
        return n;
    }

    double step = total / (double)(cap - n);
    if (tol_xyz > 0.0 && 8.0 * tol_xyz > step) step = 8.0 * tol_xyz;

    int out = 0;
    dst[out++] = src[0];
    for (int k = 0; k < n - 1; k++) {
        const liscio_pose_t *a = &src[k];
        const liscio_pose_t *b = &src[k + 1];
        double dx = b->x - a->x, dy = b->y - a->y, dz = b->z - a->z;
        double L = sqrt(dx*dx + dy*dy + dz*dz);
        int m = (int)(L / step);         /* interior samples on this span */
        if (m >= 1) {
            double inv = 1.0 / (double)(m + 1);
            for (int j = 1; j <= m; j++) {
                double f = (double)j * inv;
                liscio_pose_t p;
                p.x = a->x + f * dx;
                p.y = a->y + f * dy;
                p.z = a->z + f * dz;
                p.a = a->a + f * (b->a - a->a);
                p.b = a->b + f * (b->b - a->b);
                p.c = a->c + f * (b->c - a->c);
                p.u = a->u + f * (b->u - a->u);
                p.v = a->v + f * (b->v - a->v);
                p.w = a->w + f * (b->w - a->w);
                dst[out++] = p;
            }
        }
        dst[out++] = *b;
    }
    return out;                          /* <= cap by construction of step */
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

    /* Densify long spans so the LSQ and the tol gate see chord interiors
     * (see LISCIO_FIT_DENSIFY_CAP above).  Fallback: original points. */
    liscio_pose_t dpts[LISCIO_FIT_DENSIFY_CAP];
    int dn = densify_window(pts, n, ctx->cfg.tol_xyz,
                            dpts, LISCIO_FIT_DENSIFY_CAP);
    if (dn >= n) { pts = dpts; n = dn; }

    /* Chord-length parameterization. */
    double ts[LISCIO_FIT_DENSIFY_CAP];
    double cum_len[LISCIO_FIT_DENSIFY_CAP];
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
    double coords[9][LISCIO_FIT_DENSIFY_CAP];
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

    /* Densify long spans (see LISCIO_FIT_DENSIFY_CAP).  Endpoints are
     * bit-copies of the originals, so the given G1 tangents keep their
     * meaning.  Fallback: original points. */
    liscio_pose_t dpts[LISCIO_FIT_DENSIFY_CAP];
    int dn = densify_window(pts, n, ctx->cfg.tol_xyz,
                            dpts, LISCIO_FIT_DENSIFY_CAP);
    if (dn >= n) { pts = dpts; n = dn; }

    /* Chord-length parameterization. */
    double ts[LISCIO_FIT_DENSIFY_CAP];
    double cum_len[LISCIO_FIT_DENSIFY_CAP];
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
    double xs[LISCIO_FIT_DENSIFY_CAP], ys[LISCIO_FIT_DENSIFY_CAP], zs[LISCIO_FIT_DENSIFY_CAP];
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
    double abc[6][LISCIO_FIT_DENSIFY_CAP];
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

/* ============================================================
 *  Quintic G2 fit — endpoint curvature pinned to 0 (line seams)
 * ============================================================
 * Degree-5 Bézier whose first three / last three control points lie ON the
 * start / end tangent lines, so κ(0)=κ(5)=0: the curve leaves and re-enters its
 * neighbour straight segments with ZERO curvature → no centripetal-jerk step at
 * the bezier↔line seam. Four scalar DOF (a,e along the start tangent; f,b along
 * the end tangent) → a 4×4 linear least squares (the natural extension of the
 * cubic G1 αL/αR fit). xyz carry the curvature; abc/uvw are interpolated linearly
 * in the control polygon (exact for the constant/linear rotary the corpus uses).
 * Emitted as a single-span clamped degree-5 B-spline. Returns 0 on success; -1 if
 * a κ=0 quintic cannot meet tolerance (region genuinely curved at a seam → caller
 * keeps the cubic, and the residual κ-step is handled by the normal-jerk cap). */
static inline double Q50(double t){double u=1.0-t;return u*u*u*u*u;}
static inline double Q51(double t){double u=1.0-t;return 5.0*t*u*u*u*u;}
static inline double Q52(double t){double u=1.0-t;return 10.0*t*t*u*u*u;}
static inline double Q53(double t){double u=1.0-t;return 10.0*t*t*t*u*u;}
static inline double Q54(double t){double u=1.0-t;return 5.0*t*t*t*t*u;}
static inline double Q55(double t){return t*t*t*t*t;}

static int solve4x4(double A[4][4], double rhs[4], double x[4])
{
    for (int c = 0; c < 4; c++) {
        int piv = c; double mx = fabs(A[c][c]);
        for (int r = c+1; r < 4; r++) if (fabs(A[r][c]) > mx) { mx = fabs(A[r][c]); piv = r; }
        if (mx < 1e-18) return -1;
        if (piv != c) {
            for (int k = 0; k < 4; k++) { double t = A[c][k]; A[c][k] = A[piv][k]; A[piv][k] = t; }
            double t = rhs[c]; rhs[c] = rhs[piv]; rhs[piv] = t;
        }
        for (int r = c+1; r < 4; r++) {
            double fct = A[r][c] / A[c][c];
            for (int k = c; k < 4; k++) A[r][k] -= fct * A[c][k];
            rhs[r] -= fct * rhs[c];
        }
    }
    for (int r = 3; r >= 0; r--) {
        double s = rhs[r];
        for (int k = r+1; k < 4; k++) s -= A[r][k] * x[k];
        x[r] = s / A[r][r];
    }
    return 0;
}

int liscio_quintic9_g2zero_fit(const struct liscio_ctx *ctx, int i0, int i1,
    double t1x, double t1y, double t1z,
    double t2x, double t2y, double t2z,
    liscio_bspline9_fit_t *out)
{
    /* TODO(densify): same waypoint-only residual hole as the cubic
     * fitters had — apply densify_window() here before enabling
     * (env-gated OFF in TP2 production: enable_g2_quintic/curvature). */
    if (!ctx || !out) return -1;
    int n = i1 - i0 + 1;
    if (n < 4 || n > LISCIO_MAX_WINDOW) return -1;
    double m1 = sqrt(t1x*t1x+t1y*t1y+t1z*t1z), m2 = sqrt(t2x*t2x+t2y*t2y+t2z*t2z);
    if (m1 < 1e-15 || m2 < 1e-15) return -1;
    t1x/=m1; t1y/=m1; t1z/=m1;  t2x/=m2; t2y/=m2; t2z/=m2;

    const liscio_pose_t *pts = &ctx->pts[i0];
    double ts[LISCIO_MAX_WINDOW], cl[LISCIO_MAX_WINDOW]; cl[0] = 0.0;
    for (int i = 1; i < n; i++) {
        double dx=pts[i].x-pts[i-1].x, dy=pts[i].y-pts[i-1].y, dz=pts[i].z-pts[i-1].z;
        cl[i] = cl[i-1] + sqrt(dx*dx+dy*dy+dz*dz);
    }
    double total = cl[n-1]; if (total < 1e-12) return -1;
    for (int i = 0; i < n; i++) ts[i] = cl[i]/total;

    double P0x=pts[0].x,P0y=pts[0].y,P0z=pts[0].z;
    double P5x=pts[n-1].x,P5y=pts[n-1].y,P5z=pts[n-1].z;

    /* 4×4 normal equations for (a,e,f,b). */
    double M[4][4] = {{0}}, rhs[4] = {0};
    for (int i = 0; i < n; i++) {
        double t=ts[i];
        double b0=Q50(t),b1=Q51(t),b2=Q52(t),b3=Q53(t),b4=Q54(t),b5=Q55(t);
        double cP0=b0+b1+b2, cP5=b3+b4+b5;
        double rx=pts[i].x-cP0*P0x-cP5*P5x;
        double ry=pts[i].y-cP0*P0y-cP5*P5y;
        double rz=pts[i].z-cP0*P0z-cP5*P5z;
        double co[4]={b1,b2,b3,b4};
        double dx[4]={t1x,t1x,t2x,t2x}, dy[4]={t1y,t1y,t2y,t2y}, dz[4]={t1z,t1z,t2z,t2z};
        for (int k=0;k<4;k++) {
            double gkx=co[k]*dx[k], gky=co[k]*dy[k], gkz=co[k]*dz[k];
            rhs[k]+=gkx*rx+gky*ry+gkz*rz;
            for (int l=0;l<4;l++) {
                double glx=co[l]*dx[l], gly=co[l]*dy[l], glz=co[l]*dz[l];
                M[k][l]+=gkx*glx+gky*gly+gkz*glz;
            }
        }
    }
    double sol[4];
    if (solve4x4(M, rhs, sol) != 0) return -1;
    double a=sol[0], e=sol[1], f=sol[2], b=sol[3];
    double eps = 1e-6 * total;
    /* Control polygon must be FORWARD-MONOTONE on each tangent line: Q0→Q1→Q2 is
     * 0<a<e and Q3→Q4→Q5 is f>b>0. A fold (e<a or b>f) places a control point
     * behind its neighbour → a tiny in-tolerance wiggle but a huge spurious
     * curvature spike (κ→hundreds). Reject → keep the cubic. Cap handles to the
     * chord so a runaway LSQ solution can't balloon. */
    if (a < eps || b < eps || e < a + eps || f < b + eps) return -1;
    if (e > 0.9*total || f > 0.9*total) return -1;

    double Qx[6]={P0x, P0x+a*t1x, P0x+e*t1x, P5x+f*t2x, P5x+b*t2x, P5x};
    double Qy[6]={P0y, P0y+a*t1y, P0y+e*t1y, P5y+f*t2y, P5y+b*t2y, P5y};
    double Qz[6]={P0z, P0z+a*t1z, P0z+e*t1z, P5z+f*t2z, P5z+b*t2z, P5z};

    double maxdev=0;
    for (int i=0;i<n;i++) {
        double t=ts[i]; double bb[6]={Q50(t),Q51(t),Q52(t),Q53(t),Q54(t),Q55(t)};
        double bx=0,by=0,bz=0;
        for (int j=0;j<6;j++){bx+=bb[j]*Qx[j];by+=bb[j]*Qy[j];bz+=bb[j]*Qz[j];}
        double dx=bx-pts[i].x,dy=by-pts[i].y,dz=bz-pts[i].z;
        double d=sqrt(dx*dx+dy*dy+dz*dz); if(d>maxdev)maxdev=d;
    }
    if (maxdev > ctx->cfg.tol_xyz) return -1;

    /* Interior-curvature guard: a monotone control polygon can still bend the
     * MIDDLE (Q2→Q3 transition) into a near-loop — a small in-tolerance wiggle
     * that spikes κ to hundreds and would blow joint accel/jerk in RT. Sample κ
     * via the analytic 1st/2nd derivatives and reject if it exceeds a sane bound
     * (real corpus curvature is ≲ a few /mm); the cubic then stands. */
    {
        double maxk = 0.0;
        for (int s = 1; s < 64; s++) {
            double t = s/64.0, u = 1.0-t;
            /* B'(t), B''(t) of a quintic from control points (deg-4/deg-3 diffs). */
            double d1[3], d2[3];
            double Qc[3][6] = {{Qx[0],Qx[1],Qx[2],Qx[3],Qx[4],Qx[5]},
                               {Qy[0],Qy[1],Qy[2],Qy[3],Qy[4],Qy[5]},
                               {Qz[0],Qz[1],Qz[2],Qz[3],Qz[4],Qz[5]}};
            /* deg-4 Bernstein for B', deg-3 for B'' */
            double e4[5]={u*u*u*u, 4*t*u*u*u, 6*t*t*u*u, 4*t*t*t*u, t*t*t*t};
            double e3[4]={u*u*u, 3*t*u*u, 3*t*t*u, t*t*t};
            for (int c=0;c<3;c++){
                double s1=0,s2=0;
                for (int j=0;j<5;j++) s1 += e4[j]*5.0*(Qc[c][j+1]-Qc[c][j]);
                for (int j=0;j<4;j++) s2 += e3[j]*20.0*(Qc[c][j+2]-2*Qc[c][j+1]+Qc[c][j]);
                d1[c]=s1; d2[c]=s2;
            }
            double cx=d1[1]*d2[2]-d1[2]*d2[1], cy=d1[2]*d2[0]-d1[0]*d2[2], cz=d1[0]*d2[1]-d1[1]*d2[0];
            double sp1=sqrt(d1[0]*d1[0]+d1[1]*d1[1]+d1[2]*d1[2]);
            if (sp1 > 1e-9) {
                double k = sqrt(cx*cx+cy*cy+cz*cz)/(sp1*sp1*sp1);
                if (k > maxk) maxk = k;
            }
        }
        if (maxk > 20.0) return -1;   /* fold/loop → keep the cubic */
    }

    /* abc/uvw: linear-in-polygon, verify tolerance. */
    const double *sp=(const double*)&pts[0], *ep=(const double*)&pts[n-1];
    double st_r[6], en_r[6];
    for (int d=0;d<6;d++){ st_r[d]=sp[3+d]; en_r[d]=ep[3+d]; }
    for (int i=0;i<n;i++){
        double t=ts[i]; const double *pi=(const double*)&pts[i];
        for (int d=0;d<6;d++){
            double v=st_r[d]+(en_r[d]-st_r[d])*t;
            double tol=(d<3)?ctx->cfg.tol_abc:ctx->cfg.tol_uvw;
            if (fabs(v-pi[3+d])>tol) return -1;
        }
    }

    memset(out, 0, sizeof(*out));
    out->degree = 5; out->n_ctrl = 6;
    for (int j=0;j<6;j++){
        double *cj=(double*)&out->ctrl[j];
        cj[0]=Qx[j]; cj[1]=Qy[j]; cj[2]=Qz[j];
        for (int d=0;d<6;d++) cj[3+d]=st_r[d]+(en_r[d]-st_r[d])*(j/5.0);
        out->weights[j]=1.0;
    }
    for (int j=0;j<6;j++){ out->knots[j]=0.0; out->knots[6+j]=1.0; }
    out->max_deviation = maxdev;
    return 0;
}

/* ============================================================
 *  Generalized quintic G2 fit — pins ARBITRARY endpoint curvature (κ≠0)  [MARK]
 * ============================================================
 * Brings the bezier9 curvature-matching (blend/bezier9.cc δ=5α²κ/4 normal
 * offset) INTO liscio so fitted curves are C2 curvature-continuous at seams
 * (each curve continues the neighbour's κ) instead of unbending to κ=0.  For an
 * ARC run (subseg points) this recovers κ=1/R at the seams → no curvature step →
 * no ripple → constant velocity → low jerk.  Endpoint κ + curve-normal are
 * estimated from the first/last 3 points (osculating circle) so the fit is
 * self-contained (arcs AND curved micro-line runs). */
static double est_endpoint_kappa(const double pa[3], const double pb[3],
                                 const double pc[3], double n[3])
{
    double ux=pb[0]-pa[0], uy=pb[1]-pa[1], uz=pb[2]-pa[2];
    double vx=pc[0]-pa[0], vy=pc[1]-pa[1], vz=pc[2]-pa[2];
    double cx=uy*vz-uz*vy, cy=uz*vx-ux*vz, cz=ux*vy-uy*vx;
    double cm=sqrt(cx*cx+cy*cy+cz*cz);
    double um=sqrt(ux*ux+uy*uy+uz*uz);
    double vm=sqrt(vx*vx+vy*vy+vz*vz);
    double wx=pc[0]-pb[0], wy=pc[1]-pb[1], wz=pc[2]-pb[2];
    double wm=sqrt(wx*wx+wy*wy+wz*wz);
    n[0]=n[1]=n[2]=0.0;
    if (cm < 1e-18 || um < 1e-15 || vm < 1e-15 || wm < 1e-15) return 0.0;
    /* Circumradius R = abc/(4A), A = triangle area = ½·cm (cm=|AB×AC|=2A).
     * So κ = 1/R = 4A/(abc) = 2·cm/(abc).  (Prior code dropped the ×2 → κ/2,
     * which under-bent the g2 quintic construction → mid-curve bulge →
     * BULGE-reject on sharp arcs → bisect-to-native crawl.  See arc-chain
     * G2 fusion.) */
    double kappa = 2.0 * cm / (um * vm * wm);   /* κ = 1/R */
    double uhx=ux/um, uhy=uy/um, uhz=uz/um;
    double vdotu = vx*uhx+vy*uhy+vz*uhz;
    double px=vx-vdotu*uhx, py=vy-vdotu*uhy, pz=vz-vdotu*uhz;
    double pm=sqrt(px*px+py*py+pz*pz);
    if (pm < 1e-15) return 0.0;
    n[0]=px/pm; n[1]=py/pm; n[2]=pz/pm;
    return kappa;
}

int liscio_quintic9_g2_fit(const struct liscio_ctx *ctx, int i0, int i1,
    double t1x, double t1y, double t1z,
    double t2x, double t2y, double t2z,
    liscio_bspline9_fit_t *out)
{
    /* TODO(densify): same waypoint-only residual hole as the cubic
     * fitters had — apply densify_window() here before enabling
     * (env-gated OFF in TP2 production: enable_g2_quintic/curvature). */
    if (!ctx || !out) return -1;
    if (getenv("TP2_G2REJECT")) return -1;   /* diag: force-reject to isolate wedge */
    int n = i1 - i0 + 1;
    if (n < 4 || n > LISCIO_MAX_WINDOW) return -1;
    double m1 = sqrt(t1x*t1x+t1y*t1y+t1z*t1z), m2 = sqrt(t2x*t2x+t2y*t2y+t2z*t2z);
    if (m1 < 1e-15 || m2 < 1e-15) return -1;
    t1x/=m1; t1y/=m1; t1z/=m1;  t2x/=m2; t2y/=m2; t2z/=m2;

    const liscio_pose_t *pts = &ctx->pts[i0];
    double ts[LISCIO_MAX_WINDOW], cl[LISCIO_MAX_WINDOW]; cl[0] = 0.0;
    for (int i = 1; i < n; i++) {
        double dx=pts[i].x-pts[i-1].x, dy=pts[i].y-pts[i-1].y, dz=pts[i].z-pts[i-1].z;
        cl[i] = cl[i-1] + sqrt(dx*dx+dy*dy+dz*dz);
    }
    double total = cl[n-1]; if (total < 1e-12) return -1;
    for (int i = 0; i < n; i++) ts[i] = cl[i]/total;

    double P0x=pts[0].x,P0y=pts[0].y,P0z=pts[0].z;
    double P5x=pts[n-1].x,P5y=pts[n-1].y,P5z=pts[n-1].z;

    double ns[3], ne[3];
    double pa0[3]={pts[0].x,pts[0].y,pts[0].z};
    double pb0[3]={pts[1].x,pts[1].y,pts[1].z};
    double pc0[3]={pts[2].x,pts[2].y,pts[2].z};
    double k_start = est_endpoint_kappa(pa0, pb0, pc0, ns);
    double pa1[3]={pts[n-1].x,pts[n-1].y,pts[n-1].z};
    double pb1[3]={pts[n-2].x,pts[n-2].y,pts[n-2].z};
    double pc1[3]={pts[n-3].x,pts[n-3].y,pts[n-3].z};
    double k_end = est_endpoint_kappa(pa1, pb1, pc1, ne);

    /* G2 CHAIN: if the previous emitted primitive was a curvature-matched
     * quintic, PIN this fit's start curvature+normal to its end → the two
     * curves share the seam curvature (C2), no κ step, no velocity hunting.
     * Requires the G1 tangent chain too (G2 ⇒ G1); a corner/event that broke
     * the tangent (have_prev_emit_tan=0) also breaks the G2 chain here, so the
     * seam falls to a blend join instead — the user's "split where you can't
     * hold G2, let the blend connect" design. */
    if (ctx->have_prev_emit_kappa && ctx->have_prev_emit_tan) {
        k_start = ctx->prev_emit_kappa;
        ns[0] = ctx->prev_emit_norm_x;
        ns[1] = ctx->prev_emit_norm_y;
        ns[2] = ctx->prev_emit_norm_z;
    }

    double a, e, f, b;
    if (getenv("TP2_G2LSQ")) {
        /* Legacy facet+LSQ leg solve — DISPROVEN (bumps → jerk 94k/173k vs 81k
         * native); kept for A/B only.  See git history / memory. */
        double M[4][4] = {{0}}, rhs[4] = {0};
        for (int i = 0; i < n; i++) {
            double t=ts[i];
            double b1=Q51(t),b2=Q52(t),b3=Q53(t),b4=Q54(t);
            double b0=Q50(t),b5=Q55(t);
            double cP0=b0+b1+b2, cP5=b3+b4+b5;
            double rx=pts[i].x-cP0*P0x-cP5*P5x;
            double ry=pts[i].y-cP0*P0y-cP5*P5y;
            double rz=pts[i].z-cP0*P0z-cP5*P5z;
            double co[4]={b1,b2,b3,b4};
            double dx[4]={t1x,t1x,t2x,t2x}, dy[4]={t1y,t1y,t2y,t2y}, dz[4]={t1z,t1z,t2z,t2z};
            for (int k=0;k<4;k++) {
                double gkx=co[k]*dx[k], gky=co[k]*dy[k], gkz=co[k]*dz[k];
                rhs[k]+=gkx*rx+gky*ry+gkz*rz;
                for (int l=0;l<4;l++) {
                    double glx=co[l]*dx[l], gly=co[l]*dy[l], glz=co[l]*dz[l];
                    M[k][l]+=gkx*glx+gky*gly+gkz*glz;
                }
            }
        }
        double sol[4];
        if (solve4x4(M, rhs, sol) != 0) return -1;
        a=sol[0]; e=sol[1]; f=sol[2]; b=sol[3];
    } else {
        /* DIRECT construction (no LSQ, no facet noise): evenly-spaced control
         * legs along the chord; the δ curvature offset (below) bends P2/P3 to
         * the pinned endpoint κ.  Reconstructs a gentle arc / smooth curved run
         * cleanly.  A run too bent for ONE quintic busts the deviation check
         * below → return -1 → caller splits (composite) → blend joins the
         * split.  This is the direct-geometry path the old LSQ note called for. */
        a = 0.2 * total;  e = 0.4 * total;
        b = 0.2 * total;  f = 0.4 * total;
    }
    double eps = 1e-6 * total;
    if (a < eps || b < eps || e < a + eps || f < b + eps) return -1;
    if (e > 0.9*total || f > 0.9*total) return -1;

    double d0 = 1.25 * a * a * k_start;   /* Quintic δ = 5·a²·κ/4 */
    double d1 = 1.25 * b * b * k_end;
    double Qx[6]={P0x, P0x+a*t1x, P0x+e*t1x + d0*ns[0], P5x+f*t2x + d1*ne[0], P5x+b*t2x, P5x};
    double Qy[6]={P0y, P0y+a*t1y, P0y+e*t1y + d0*ns[1], P5y+f*t2y + d1*ne[1], P5y+b*t2y, P5y};
    double Qz[6]={P0z, P0z+a*t1z, P0z+e*t1z + d0*ns[2], P5z+f*t2z + d1*ne[2], P5z+b*t2z, P5z};

    if (getenv("TP2_G2DBG")) {
        static int c=0;
        if (c++ < 40) { FILE*gf=fopen("/tmp/tp2_g2dbg.log","a");
        if(gf){fprintf(gf,"g2fit n=%d tot=%.3f a=%.3f b=%.3f e=%.3f f=%.3f k0=%.4f k1=%.4f d0=%.4f d1=%.4f ns=(%.2f,%.2f,%.2f)\n",
            n,total,a,b,e,f,k_start,k_end,d0,d1,ns[0],ns[1],ns[2]); fclose(gf);} }
    }

    double maxdev=0;
    for (int i=0;i<n;i++) {
        double t=ts[i]; double bb[6]={Q50(t),Q51(t),Q52(t),Q53(t),Q54(t),Q55(t)};
        double bx=0,by=0,bz=0;
        for (int j=0;j<6;j++){bx+=bb[j]*Qx[j];by+=bb[j]*Qy[j];bz+=bb[j]*Qz[j];}
        double dx=bx-pts[i].x,dy=by-pts[i].y,dz=bz-pts[i].z;
        double d=sqrt(dx*dx+dy*dy+dz*dz); if(d>maxdev)maxdev=d;
    }
    {
        /* Reconstruction gate: tol_arcfit (fusion, decoupled from G64 P —
         * P must never coarsen rebuilt arc geometry); fallback tol_xyz. */
        double fit_tol = (ctx->cfg.tol_arcfit > 0.0)
                         ? ctx->cfg.tol_arcfit : ctx->cfg.tol_xyz;
        if (maxdev > fit_tol) return -1;
    }

    {
        double maxk = 0.0;
        for (int s = 1; s < 64; s++) {
            double t = s/64.0, u = 1.0-t;
            double d1v[3], d2v[3];
            double Qc[3][6] = {{Qx[0],Qx[1],Qx[2],Qx[3],Qx[4],Qx[5]},
                               {Qy[0],Qy[1],Qy[2],Qy[3],Qy[4],Qy[5]},
                               {Qz[0],Qz[1],Qz[2],Qz[3],Qz[4],Qz[5]}};
            double e4[5]={u*u*u*u, 4*t*u*u*u, 6*t*t*u*u, 4*t*t*t*u, t*t*t*t};
            double e3[4]={u*u*u, 3*t*u*u, 3*t*t*u, t*t*t};
            for (int c=0;c<3;c++){
                double s1=0,s2=0;
                for (int j=0;j<5;j++) s1 += e4[j]*5.0*(Qc[c][j+1]-Qc[c][j]);
                for (int j=0;j<4;j++) s2 += e3[j]*20.0*(Qc[c][j+2]-2*Qc[c][j+1]+Qc[c][j]);
                d1v[c]=s1; d2v[c]=s2;
            }
            double cx=d1v[1]*d2v[2]-d1v[2]*d2v[1], cy=d1v[2]*d2v[0]-d1v[0]*d2v[2], cz=d1v[0]*d2v[1]-d1v[1]*d2v[0];
            double sp1=sqrt(d1v[0]*d1v[0]+d1v[1]*d1v[1]+d1v[2]*d1v[2]);
            if (sp1 > 1e-9) {
                double k = sqrt(cx*cx+cy*cy+cz*cz)/(sp1*sp1*sp1);
                if (k > maxk) maxk = k;
            }
        }
        if (maxk > 20.0) return -1;
        /* BULGE REJECT: a curvature-matched fit of an ARC-like run (endpoints
         * share a real curvature) should stay near that curvature.  A mid-curve
         * bulge > BULGE_MAX × endpoint κ means one quintic cannot hold this span
         * (thomam's ~150° arc window bulged maxk to 2.3× endpoint κ → the RT
         * centripetal cap throttled it → 38s crawl).  Reject so the composite
         * splitter makes SMALLER pieces (each a near-arc, no bulge), which the
         * ctx κ-chain then stitches C2 → continuous curvature, no seam step.
         * Only for arc-like runs (endpoint κ significant); low-κ gentle runs are
         * exempt so we don't over-split them.  Tunable / off via env. */
        {
            double kref = fmax(k_start, k_end);
            double bmax = getenv("TP2_G2BULGE") ? atof(getenv("TP2_G2BULGE")) : 1.4;
            if (bmax > 1.0 && kref > 0.02 && maxk > bmax * kref) return -1;
        }
        /* 2026-09-05: G1 DROP REJECT (pinned k_start vs 8–15% avg < 0.5)
         * → cubic/blend, 2_1001.sweep ±306k @ ~6 mm/s.  Do not reopen.
         * 2026-09-05: G1 C2-split on that drop (keep quintic if split fails)
         * → 2_1001.sweep ±1.34M @ Z=−42.53.  Do not reopen.
         * 2026-09-05: one-TC Boehm-insert + snap Q[3] at the dump
         * → take1 ±158k @ 70 mm/s same Z (FO-up), take2 pass was FO luck.
         * Do not reopen. */
    }

    const double *sp=(const double*)&pts[0], *ep=(const double*)&pts[n-1];
    double st_r[6], en_r[6];
    for (int d=0;d<6;d++){ st_r[d]=sp[3+d]; en_r[d]=ep[3+d]; }
    for (int i=0;i<n;i++){
        double t=ts[i]; const double *pi=(const double*)&pts[i];
        for (int d=0;d<6;d++){
            double v=st_r[d]+(en_r[d]-st_r[d])*t;
            double tol=(d<3)?ctx->cfg.tol_abc:ctx->cfg.tol_uvw;
            if (fabs(v-pi[3+d])>tol) return -1;
        }
    }

    memset(out, 0, sizeof(*out));
    out->degree = 5; out->n_ctrl = 6;
    for (int j=0;j<6;j++){
        double *cj=(double*)&out->ctrl[j];
        cj[0]=Qx[j]; cj[1]=Qy[j]; cj[2]=Qz[j];
        for (int d=0;d<6;d++) cj[3+d]=st_r[d]+(en_r[d]-st_r[d])*(j/5.0);
        out->weights[j]=1.0;
    }
    for (int j=0;j<6;j++){ out->knots[j]=0.0; out->knots[6+j]=1.0; }
    out->max_deviation = maxdev;
    /* Publish end curvature+normal for the ctx G2 chain (next curve's k_start). */
    out->end_kappa  = k_end;
    out->end_norm_x = ne[0];
    out->end_norm_y = ne[1];
    out->end_norm_z = ne[2];
    return 0;
}
