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
* Description: liscio 9D arc detection + Kasa least-squares fit.
*              Algorithm:
*                1. Pick 3 most-spread points to establish XYZ plane.
*                2. Verify all N points within tol_xyz of that plane.
*                3. Project to 2D (u,v); fit circle by Kasa linearization
*                     xi² + yi² = 2a·xi + 2b·yi + c   (3x3 normal eq.)
*                4. Verify all projected points within tol_xyz of fit.
*                5. Verify ABC/UVW linear in arc length.
*                6. Compute signed arc_angle (start → end).
*              Pure C99, no allocations (stack buffers only).
* References:  J. Kasa, *A circle fitting procedure and its error
*              analysis*, IEEE Trans. Instrum. Meas., 1976.
*              N. Chernov & C. Lesort, *Least squares fitting of
*              circles*, J. Math. Imaging Vis., 2005.
* Author:      杨阳 (Yang Yang) <mika-net@outlook.com>
* License:     MIT (SPDX-License-Identifier: MIT)
* Copyright (c) 2026 杨阳 (Yang Yang)
********************************************************************/

#define _USE_MATH_DEFINES
#include "liscio/liscio.h"
#include "liscio_internal.h"

#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static inline double sqr(double v) { return v * v; }

/* 3D vector helpers */
static inline void v_sub(double *o, const double *a, const double *b) {
    o[0] = a[0] - b[0];
    o[1] = a[1] - b[1];
    o[2] = a[2] - b[2];
}
static inline void v_cross(double *o, const double *a, const double *b) {
    o[0] = a[1]*b[2] - a[2]*b[1];
    o[1] = a[2]*b[0] - a[0]*b[2];
    o[2] = a[0]*b[1] - a[1]*b[0];
}
static inline double v_dot(const double *a, const double *b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
static inline double v_norm(const double *a) {
    return sqrt(v_dot(a, a));
}
static inline int v_normalize(double *a) {
    double n = v_norm(a);
    if (n < 1e-15) return -1;
    a[0] /= n; a[1] /= n; a[2] /= n;
    return 0;
}

/* 3x3 linear solve Ax = b via Gaussian elimination with partial pivot.
 * A is row-major. Returns 0 on success, -1 on singular. */
static int solve3x3(double A[3][3], double b[3], double x[3])
{
    /* Augment */
    double M[3][4];
    for (int i = 0; i < 3; i++) {
        M[i][0] = A[i][0]; M[i][1] = A[i][1]; M[i][2] = A[i][2];
        M[i][3] = b[i];
    }
    for (int k = 0; k < 3; k++) {
        /* Pivot */
        int piv = k;
        double best = fabs(M[k][k]);
        for (int r = k + 1; r < 3; r++) {
            if (fabs(M[r][k]) > best) { best = fabs(M[r][k]); piv = r; }
        }
        if (best < 1e-15) return -1;
        if (piv != k) {
            for (int c = 0; c < 4; c++) {
                double t = M[k][c]; M[k][c] = M[piv][c]; M[piv][c] = t;
            }
        }
        /* Eliminate below */
        for (int r = k + 1; r < 3; r++) {
            double f = M[r][k] / M[k][k];
            for (int c = k; c < 4; c++)
                M[r][c] -= f * M[k][c];
        }
    }
    /* Back-sub */
    for (int i = 2; i >= 0; i--) {
        double s = M[i][3];
        for (int j = i + 1; j < 3; j++) s -= M[i][j] * x[j];
        x[i] = s / M[i][i];
    }
    return 0;
}

/* Find plane normal from span of window points.  Uses first point as
 * origin; picks two longest-span pairs for stability. */
static int plane_normal(const liscio_pose_t *pts, int n, double n_out[3])
{
    if (n < 3) return -1;

    double p0[3] = { pts[0].x, pts[0].y, pts[0].z };
    /* Find point furthest from p0 */
    int iA = -1;
    double max1 = -1.0;
    for (int i = 1; i < n; i++) {
        double d[3] = { pts[i].x - p0[0], pts[i].y - p0[1], pts[i].z - p0[2] };
        double m = v_norm(d);
        if (m > max1) { max1 = m; iA = i; }
    }
    if (iA < 0 || max1 < 1e-12) return -1;

    double dA[3] = { pts[iA].x - p0[0], pts[iA].y - p0[1], pts[iA].z - p0[2] };
    v_normalize(dA);

    /* Find point maximizing |cross(dA, dB)| = point furthest from line p0-pA */
    int iB = -1;
    double max2 = -1.0;
    double best_n[3] = { 0, 0, 0 };
    for (int i = 1; i < n; i++) {
        if (i == iA) continue;
        double dB[3] = { pts[i].x - p0[0], pts[i].y - p0[1], pts[i].z - p0[2] };
        double cr[3];
        v_cross(cr, dA, dB);
        double m = v_norm(cr);
        if (m > max2) {
            max2 = m; iB = i;
            best_n[0] = cr[0]; best_n[1] = cr[1]; best_n[2] = cr[2];
        }
    }
    if (iB < 0 || max2 < 1e-9) return -1;  /* degenerate (all collinear) */

    if (v_normalize(best_n) != 0) return -1;
    n_out[0] = best_n[0]; n_out[1] = best_n[1]; n_out[2] = best_n[2];
    return 0;
}

/* Build orthonormal plane basis (u,v) given normal n.  u perpendicular
 * to n and to world Z if possible; else to world X. */
static void plane_basis(const double n[3], double u[3], double v[3])
{
    double z[3] = { 0, 0, 1 };
    double cu[3];
    v_cross(cu, n, z);
    if (v_norm(cu) < 1e-6) {
        double x[3] = { 1, 0, 0 };
        v_cross(cu, n, x);
    }
    v_normalize(cu);
    u[0] = cu[0]; u[1] = cu[1]; u[2] = cu[2];
    double cv[3];
    v_cross(cv, n, u);
    v_normalize(cv);
    v[0] = cv[0]; v[1] = cv[1]; v[2] = cv[2];
}

/* Fit circle (Kasa) to 2D points.  Returns 0 on success with center
 * (cx,cy) and radius r.  */
static int fit_circle_kasa(const double *xs, const double *ys, int n,
                           double *cx, double *cy, double *r)
{
    if (n < 3) return -1;

    double Sx = 0, Sy = 0, Sxx = 0, Syy = 0, Sxy = 0;
    double Sxxx = 0, Syyy = 0, Sxxy = 0, Sxyy = 0;
    for (int i = 0; i < n; i++) {
        double x = xs[i], y = ys[i];
        double x2 = x*x, y2 = y*y;
        Sx += x; Sy += y;
        Sxx += x2; Syy += y2; Sxy += x*y;
        Sxxx += x2*x; Syyy += y2*y;
        Sxxy += x2*y; Sxyy += x*y2;
    }

    /* System:
     *   [Sxx  Sxy  Sx ] [a]   [Sxxx + Sxyy]
     *   [Sxy  Syy  Sy ] [b] = [Sxxy + Syyy]  where a=2*cx, b=2*cy, c=r²-cx²-cy²
     *   [Sx   Sy   n  ] [c]   [Sxx + Syy  ]
     */
    double A[3][3] = {
        { Sxx, Sxy, Sx },
        { Sxy, Syy, Sy },
        { Sx,  Sy,  (double)n }
    };
    double B[3] = { Sxxx + Sxyy, Sxxy + Syyy, Sxx + Syy };
    double X[3];
    if (solve3x3(A, B, X) != 0) return -1;

    *cx = X[0] * 0.5;
    *cy = X[1] * 0.5;
    double r2 = X[2] + (*cx)*(*cx) + (*cy)*(*cy);
    if (r2 < 0.0) return -1;
    *r = sqrt(r2);
    return 0;
}

/* ---------- Public entry ---------- */
int liscio_arc9d_fit(const struct liscio_ctx *ctx, int i0, int i1,
                      liscio_arc9d_fit_t *out)
{
    if (!ctx || !out) return -1;
    int n = i1 - i0 + 1;
    if (n < ctx->cfg.min_arc_pts) return -1;
    if (n > LISCIO_MAX_WINDOW) return -1;

    const liscio_pose_t *pts = &ctx->pts[i0];

    /* 1. Plane. */
    double nrm[3];
    if (plane_normal(pts, n, nrm) != 0) return -1;

    /* 2. Reject if any point too far from plane (flat check). */
    double p0[3] = { pts[0].x, pts[0].y, pts[0].z };
    double max_plane_dev = 0.0;
    for (int i = 0; i < n; i++) {
        double d[3] = { pts[i].x - p0[0], pts[i].y - p0[1], pts[i].z - p0[2] };
        double off = fabs(v_dot(d, nrm));
        if (off > max_plane_dev) max_plane_dev = off;
    }
    if (max_plane_dev > ctx->cfg.tol_xyz) return -1;

    /* 3. Project to 2D. */
    double u[3], v[3];
    plane_basis(nrm, u, v);

    double xs[LISCIO_MAX_WINDOW], ys[LISCIO_MAX_WINDOW];
    for (int i = 0; i < n; i++) {
        double d[3] = { pts[i].x - p0[0], pts[i].y - p0[1], pts[i].z - p0[2] };
        xs[i] = v_dot(d, u);
        ys[i] = v_dot(d, v);
    }

    /* 4. Fit circle. */
    double cx2, cy2, radius;
    if (fit_circle_kasa(xs, ys, n, &cx2, &cy2, &radius) != 0) return -1;
    if (radius < ctx->cfg.min_arc_radius) return -1;
    if (radius > ctx->cfg.max_arc_radius) return -1;

    /* 5. Deviation check on projected points. */
    double max_dev = 0.0;
    for (int i = 0; i < n; i++) {
        double dx = xs[i] - cx2, dy = ys[i] - cy2;
        double dr = fabs(sqrt(dx*dx + dy*dy) - radius);
        if (dr > max_dev) max_dev = dr;
    }
    if (max_dev > ctx->cfg.tol_xyz) return -1;

    /* 6. ABC/UVW linearity in arc length. */
    /* Cumulative angular position. */
    double angs[LISCIO_MAX_WINDOW];
    angs[0] = atan2(ys[0] - cy2, xs[0] - cx2);
    /* Unwrap to monotone. */
    for (int i = 1; i < n; i++) {
        double a = atan2(ys[i] - cy2, xs[i] - cx2);
        double da = a - angs[i-1];
        while (da >  M_PI) { da -= 2*M_PI; }
        while (da < -M_PI) { da += 2*M_PI; }
        angs[i] = angs[i-1] + da;
    }
    double arc_total = angs[n-1] - angs[0];
    double arc_len_total = fabs(arc_total) * radius;
    if (arc_len_total < 1e-9) return -1;

    /* Check linearity for each of abc, uvw.  Linear regression slope
     * check: max deviation from linear interpolation < tol. */
    #define LINCHECK(FIELD, TOL)                                             \
    do {                                                                     \
        double v0 = pts[0].FIELD, vN = pts[n-1].FIELD;                      \
        for (int i = 1; i < n - 1; i++) {                                   \
            double t = (angs[i] - angs[0]) / arc_total;                     \
            double predicted = v0 + t * (vN - v0);                          \
            if (fabs(pts[i].FIELD - predicted) > (TOL))                     \
                return -1;                                                  \
        }                                                                    \
    } while (0)

    LINCHECK(a, ctx->cfg.tol_abc);
    LINCHECK(b, ctx->cfg.tol_abc);
    LINCHECK(c, ctx->cfg.tol_abc);
    LINCHECK(u, ctx->cfg.tol_uvw);
    LINCHECK(v, ctx->cfg.tol_uvw);
    LINCHECK(w, ctx->cfg.tol_uvw);
    #undef LINCHECK

    /* 7. Fill output. Center in world coords: p0 + cx2*u + cy2*v. */
    out->cx = p0[0] + cx2 * u[0] + cy2 * v[0];
    out->cy = p0[1] + cx2 * u[1] + cy2 * v[1];
    out->cz = p0[2] + cx2 * u[2] + cy2 * v[2];
    out->nx = nrm[0]; out->ny = nrm[1]; out->nz = nrm[2];
    out->radius = radius;
    out->arc_angle = arc_total;
    out->max_deviation = max_dev;
    return 0;
}
