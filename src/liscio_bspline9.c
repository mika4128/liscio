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
* Description: liscio 9D cubic B-spline least-squares fit.
*              Cubic (degree 3) B-spline with clamped uniform knot
*              vector and M control points (M ≥ 4).  Non-rational
*              (weights = 1).  Each XYZ / ABC / UVW dimension fit
*              independently sharing chord-length parameterization;
*              normal equations (B^T·B) c = B^T·q solved via 3x3-band
*              LDL^T Cholesky (symmetric positive-definite).
*              Pins endpoints P_0 = first waypoint, P_{M-1} = last.
*              C² continuous across interior knots automatically.
*              Reduces to composite cubic Bezier at knot spans but
*              stores as single unified primitive (fewer object
*              boundaries for downstream consumers).
* References:  L. Piegl & W. Tiller, *The NURBS Book*, 2nd ed.,
*              Springer 1997 — §9.4 (global LSQ curve fitting).
*              C. de Boor, *A Practical Guide to Splines*, Springer
*              1978 — basis evaluation via Cox-de Boor recursion.
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

/* Cox-de Boor: evaluate non-zero cubic B-spline basis functions
 * N_{i-3..i,3}(t) at parameter t, given span index i and clamped
 * uniform knot vector U[]. Writes 4 values to N_out. */
static void basis_funcs_cubic(int i, double t, const double *U, double *N_out)
{
    double left[4], right[4];
    N_out[0] = 1.0;
    for (int j = 1; j <= 3; j++) {
        left[j]  = t - U[i + 1 - j];
        right[j] = U[i + j] - t;
        double saved = 0.0;
        for (int r = 0; r < j; r++) {
            double denom = right[r + 1] + left[j - r];
            if (denom < 1e-15) { N_out[r] = 0.0; continue; }
            double temp = N_out[r] / denom;
            N_out[r] = saved + right[r + 1] * temp;
            saved    = left[j - r] * temp;
        }
        N_out[j] = saved;
    }
}

/* Find knot span index for t in clamped knot vector U[] of M+4 knots
 * (M = #ctrl pts, degree 3). Returns i such that U[i] ≤ t < U[i+1]. */
static int find_span(int M, double t, const double *U)
{
    if (t >= U[M]) return M - 1;
    if (t <= U[3]) return 3;
    int lo = 3, hi = M, mid = (lo + hi) / 2;
    while (t < U[mid] || t >= U[mid + 1]) {
        if (t < U[mid]) hi = mid; else lo = mid;
        mid = (lo + hi) / 2;
    }
    return mid;
}

/* Solve symmetric positive-definite system A·x = b (size n×n) via
 * Cholesky LDL^T decomposition.  A destroyed in place.  Returns 0 on
 * success, -1 on non-positive pivot. */
static int solve_spd(double *A, double *b, double *x, int n)
{
    /* Cholesky: A = L·L^T (lower triangular). */
    for (int i = 0; i < n; i++) {
        double sum = A[i * n + i];
        for (int k = 0; k < i; k++) sum -= A[i * n + k] * A[i * n + k];
        if (sum <= 1e-15) return -1;
        A[i * n + i] = sqrt(sum);
        for (int j = i + 1; j < n; j++) {
            double s = A[j * n + i];
            for (int k = 0; k < i; k++) s -= A[j * n + k] * A[i * n + k];
            A[j * n + i] = s / A[i * n + i];
        }
    }
    /* Forward: L·y = b */
    double y[LISCIO_BSPLINE_MAX_CTRL];
    for (int i = 0; i < n; i++) {
        double s = b[i];
        for (int k = 0; k < i; k++) s -= A[i * n + k] * y[k];
        y[i] = s / A[i * n + i];
    }
    /* Backward: L^T·x = y */
    for (int i = n - 1; i >= 0; i--) {
        double s = y[i];
        for (int k = i + 1; k < n; k++) s -= A[k * n + i] * x[k];
        x[i] = s / A[i * n + i];
    }
    return 0;
}

/* LSQ fit cubic B-spline to points along one 1D dimension.
 * Returns 0 on success; M control points written to ctrl_out[0..M-1].
 * P_0 = points[0], P_{M-1} = points[n-1] pinned; unknowns are M-2
 * interior control points. */
static int bspline_lsq_1d(const double *q, const double *ts, int n,
                          const double *U, int M,
                          double *ctrl_out)
{
    if (M < 4 || n < M) return -1;

    ctrl_out[0]     = q[0];
    ctrl_out[M - 1] = q[n - 1];

    int F = M - 2;  /* # free interior control points */

    /* Build A = B_f^T · B_f, b = B_f^T · R where
     *   R_i = q_i − N_0(t_i)·P_0 − N_{M-1}(t_i)·P_{M-1}
     *   B_f: n × F matrix, B_f[i][j] = N_{j+1}(t_i)  (j = 0..F-1) */
    double A[LISCIO_BSPLINE_MAX_CTRL * LISCIO_BSPLINE_MAX_CTRL];
    double b[LISCIO_BSPLINE_MAX_CTRL];
    memset(A, 0, sizeof(A));
    memset(b, 0, sizeof(b));

    for (int i = 0; i < n; i++) {
        int span = find_span(M, ts[i], U);
        double N[4];
        basis_funcs_cubic(span, ts[i], U, N);
        /* Non-zero basis indices: span-3 .. span (4 values). */
        double N_end[4];  /* just the 4 weights mapped to CP index */
        int cp_idx[4];
        for (int k = 0; k < 4; k++) {
            cp_idx[k] = span - 3 + k;
            N_end[k]  = N[k];
        }

        /* Residual target after subtracting pinned endpoints. */
        double N0 = 0.0, NM = 0.0;
        for (int k = 0; k < 4; k++) {
            if (cp_idx[k] == 0)       N0 = N_end[k];
            if (cp_idx[k] == M - 1)   NM = N_end[k];
        }
        double R = q[i] - N0 * ctrl_out[0] - NM * ctrl_out[M - 1];

        /* Accumulate A and b over interior CPs only. */
        for (int a = 0; a < 4; a++) {
            int ia = cp_idx[a] - 1;
            if (ia < 0 || ia >= F) continue;
            b[ia] += N_end[a] * R;
            for (int c = 0; c < 4; c++) {
                int ic = cp_idx[c] - 1;
                if (ic < 0 || ic >= F) continue;
                A[ia * F + ic] += N_end[a] * N_end[c];
            }
        }
    }

    double sol[LISCIO_BSPLINE_MAX_CTRL];
    if (solve_spd(A, b, sol, F) != 0) return -1;
    for (int j = 0; j < F; j++) ctrl_out[j + 1] = sol[j];
    return 0;
}

/* Evaluate 1D cubic B-spline at t. */
static double bspline_eval_1d(const double *ctrl, int M,
                              const double *U, double t)
{
    int span = find_span(M, t, U);
    double N[4];
    basis_funcs_cubic(span, t, U, N);
    double v = 0.0;
    for (int k = 0; k < 4; k++) v += N[k] * ctrl[span - 3 + k];
    return v;
}

/* Build clamped uniform knot vector: first + last 4 clamped at 0/1,
 * M − 4 interior knots evenly spaced in (0, 1). Total size M + 4. */
static void build_clamped_uniform_knots(double *U, int M)
{
    for (int i = 0; i < 4; i++) U[i] = 0.0;
    for (int i = M; i < M + 4; i++) U[i] = 1.0;
    int n_interior = M - 4;
    for (int i = 0; i < n_interior; i++) {
        U[4 + i] = (double)(i + 1) / (double)(n_interior + 1);
    }
}

/* ---------- Public entry ---------- */
int liscio_bspline9_fit(const struct liscio_ctx *ctx, int i0, int i1,
                        int M, liscio_bspline9_fit_t *out)
{
    if (!ctx || !out) return -1;
    int n = i1 - i0 + 1;
    if (n < M || M < 4 || M > LISCIO_BSPLINE_MAX_CTRL) return -1;

    const liscio_pose_t *pts = &ctx->pts[i0];

    /* Chord-length parameterization (ts in [0,1]). */
    double ts[LISCIO_MAX_WINDOW];
    double cum[LISCIO_MAX_WINDOW];
    cum[0] = 0.0;
    for (int i = 1; i < n; i++) {
        double dx = pts[i].x - pts[i-1].x;
        double dy = pts[i].y - pts[i-1].y;
        double dz = pts[i].z - pts[i-1].z;
        cum[i] = cum[i-1] + sqrt(dx*dx + dy*dy + dz*dz);
    }
    double total = cum[n-1];
    if (total < 1e-12) return -1;
    for (int i = 0; i < n; i++) ts[i] = cum[i] / total;

    /* Clamped uniform knot vector, size M + 4. */
    double U[LISCIO_BSPLINE_MAX_CTRL + 4];
    build_clamped_uniform_knots(U, M);

    /* Fit each 9D dimension independently. */
    double coords[9][LISCIO_MAX_WINDOW];
    double ctrl[9][LISCIO_BSPLINE_MAX_CTRL];
    for (int i = 0; i < n; i++) {
        coords[0][i] = pts[i].x; coords[1][i] = pts[i].y; coords[2][i] = pts[i].z;
        coords[3][i] = pts[i].a; coords[4][i] = pts[i].b; coords[5][i] = pts[i].c;
        coords[6][i] = pts[i].u; coords[7][i] = pts[i].v; coords[8][i] = pts[i].w;
    }
    for (int d = 0; d < 9; d++) {
        /* Skip near-constant dims: pin endpoints, linear interior. */
        double mean = 0;
        for (int i = 0; i < n; i++) mean += coords[d][i];
        mean /= n;
        double var = 0;
        for (int i = 0; i < n; i++) {
            double df = coords[d][i] - mean;
            var += df*df;
        }
        if (var < 1e-18) {
            double v0 = coords[d][0], vN = coords[d][n-1];
            for (int j = 0; j < M; j++) {
                ctrl[d][j] = v0 + (vN - v0) * (double)j / (double)(M - 1);
            }
            continue;
        }
        if (bspline_lsq_1d(coords[d], ts, n, U, M, ctrl[d]) != 0)
            return -1;
    }

    /* Deviation check in XYZ at waypoints. */
    double max_dev = 0.0;
    for (int i = 0; i < n; i++) {
        double bx = bspline_eval_1d(ctrl[0], M, U, ts[i]);
        double by = bspline_eval_1d(ctrl[1], M, U, ts[i]);
        double bz = bspline_eval_1d(ctrl[2], M, U, ts[i]);
        double dx = bx - coords[0][i];
        double dy = by - coords[1][i];
        double dz = bz - coords[2][i];
        double d = sqrt(dx*dx + dy*dy + dz*dz);
        if (d > max_dev) max_dev = d;
    }
    if (max_dev > ctx->cfg.tol_xyz) return -1;

    /* Length-reasonableness check: reject global fits that oscillate
     * wildly between waypoints (LSQ with too few CPs on data with sharp
     * direction reversals).  Waypoint-tolerance alone doesn't catch
     * inter-waypoint excursions; comparing curve length to chord length
     * does.  Factor 1.5 allows genuinely curvy data. */
    double fit_len = 0.0;
    const int SAMP = 4 * n;
    double prev_bx = coords[0][0], prev_by = coords[1][0], prev_bz = coords[2][0];
    for (int s = 1; s <= SAMP; s++) {
        double tt = (double)s / (double)SAMP;
        double bx = bspline_eval_1d(ctrl[0], M, U, tt);
        double by = bspline_eval_1d(ctrl[1], M, U, tt);
        double bz = bspline_eval_1d(ctrl[2], M, U, tt);
        double dx = bx - prev_bx, dy = by - prev_by, dz = bz - prev_bz;
        fit_len += sqrt(dx*dx + dy*dy + dz*dz);
        prev_bx = bx; prev_by = by; prev_bz = bz;
    }
    if (fit_len > total * 1.5) return -1;

    /* Rotary / UVW per-dim tol check. */
    for (int d = 3; d < 9; d++) {
        double tol = (d < 6) ? ctx->cfg.tol_abc : ctx->cfg.tol_uvw;
        for (int i = 0; i < n; i++) {
            double v = bspline_eval_1d(ctrl[d], M, U, ts[i]);
            if (fabs(v - coords[d][i]) > tol) return -1;
        }
    }

    /* Fill output. */
    out->n_ctrl = M;
    for (int j = 0; j < M; j++) {
        out->ctrl[j].x = ctrl[0][j]; out->ctrl[j].y = ctrl[1][j]; out->ctrl[j].z = ctrl[2][j];
        out->ctrl[j].a = ctrl[3][j]; out->ctrl[j].b = ctrl[4][j]; out->ctrl[j].c = ctrl[5][j];
        out->ctrl[j].u = ctrl[6][j]; out->ctrl[j].v = ctrl[7][j]; out->ctrl[j].w = ctrl[8][j];
        out->weights[j] = 1.0;
    }
    memcpy(out->knots, U, (size_t)(M + 4) * sizeof(double));
    out->degree = 3;
    out->total_chord = total;
    out->max_deviation = max_dev;
    return 0;
}

/* Evaluate B-spline primitive at t ∈ [0,1]. */
void liscio_bspline_eval(const liscio_primitive_t *prim, double t,
                         double *out_x, double *out_y, double *out_z)
{
    if (!prim || !out_x || !out_y || !out_z) return;
    if (prim->type != LISCIO_PRIM_SPLINE) return;
    if (prim->n_ctrl < 4) return;

    int M = prim->n_ctrl;
    int span = find_span(M, t, prim->knots);
    double N[4];
    basis_funcs_cubic(span, t, prim->knots, N);

    double x = 0, y = 0, z = 0, w = 0;
    for (int k = 0; k < 4; k++) {
        int idx = span - 3 + k;
        double wk = prim->weights[idx];
        x += N[k] * wk * prim->ctrl_pts[idx].x;
        y += N[k] * wk * prim->ctrl_pts[idx].y;
        z += N[k] * wk * prim->ctrl_pts[idx].z;
        w += N[k] * wk;
    }
    if (fabs(w) < 1e-15) { *out_x = *out_y = *out_z = 0; return; }
    *out_x = x / w;
    *out_y = y / w;
    *out_z = z / w;
}
