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
* Description: liscio 9D helical-arc least-squares fit.  Detects when
*              a CAM-cut helix (G-code chord-approximation of a screw
*              or lathe thread) can be folded into a single helical
*              primitive instead of many short LINE/BEZIER segments.
*
*              Algorithm:
*                1. Estimate axis n via average of adjacent chord
*                   cross products (chord_i × chord_{i+1} ≈ axis on
*                   a helix).
*                2. Build orthonormal plane basis (u, v ⟂ n).
*                3. Project waypoints onto n-perpendicular plane;
*                   Kasa-fit the planar circle (cx, cy, r).
*                4. Compute axial coordinate z_i = (P_i - P0) · n and
*                   unwrapped angle θ_i; LSQ-solve z = z0 + p·θ.
*                5. Reconstruct each waypoint and gate against tol_xyz.
*
*              Reduces to planar ARC when p ≈ 0; in that case arc9d
*              fitter would have already accepted the run upstream so
*              this fitter is the fallback for helical climb data.
*
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

/* ---------- 3D vector helpers (kept private, mirror arc9d) ---------- */
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

/* Kasa circle fit on 2D points (same algorithm as arc9d).  Returns 0
 * on success with center (cx, cy) and radius r in the local frame. */
static int kasa_circle(const double *xs, const double *ys, int n,
                       double *cx_out, double *cy_out, double *r_out)
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
    /* Normal equations:
     *   [Sxx  Sxy  Sx ] [a]   [Sxxx + Sxyy]
     *   [Sxy  Syy  Sy ] [b] = [Sxxy + Syyy]   a = 2cx, b = 2cy, c = r²-cx²-cy²
     *   [Sx   Sy   N  ] [c]   [Sxx  + Syy ] */
    double A[3][3] = {
        { Sxx, Sxy, Sx },
        { Sxy, Syy, Sy },
        { Sx,  Sy,  (double)n }
    };
    double b[3] = { Sxxx + Sxyy, Sxxy + Syyy, Sxx + Syy };
    /* Gauss with pivot. */
    double M[3][4];
    for (int i = 0; i < 3; i++) {
        M[i][0] = A[i][0]; M[i][1] = A[i][1]; M[i][2] = A[i][2];
        M[i][3] = b[i];
    }
    for (int k = 0; k < 3; k++) {
        int piv = k;
        double best = fabs(M[k][k]);
        for (int r = k + 1; r < 3; r++)
            if (fabs(M[r][k]) > best) { best = fabs(M[r][k]); piv = r; }
        if (best < 1e-15) return -1;
        if (piv != k) {
            for (int c = 0; c < 4; c++) {
                double t = M[k][c]; M[k][c] = M[piv][c]; M[piv][c] = t;
            }
        }
        for (int r = k + 1; r < 3; r++) {
            double f = M[r][k] / M[k][k];
            for (int c = k; c < 4; c++) M[r][c] -= f * M[k][c];
        }
    }
    double X[3];
    for (int i = 2; i >= 0; i--) {
        double s = M[i][3];
        for (int j = i + 1; j < 3; j++) s -= M[i][j] * X[j];
        X[i] = s / M[i][i];
    }
    double cx = X[0] * 0.5;
    double cy = X[1] * 0.5;
    double r2 = X[2] + cx*cx + cy*cy;
    if (r2 < 0) return -1;
    *cx_out = cx; *cy_out = cy; *r_out = sqrt(r2);
    return 0;
}

/* Estimate helix axis as the mean of adjacent chord cross products.
 * For a helix P(θ) = c + r·u(θ) + p·θ·n, c_i × c_{i+1} is dominated by
 * r²·sin(Δθ)·n + lower-order axial-coupling terms; averaging over a
 * dense run cancels the noise.  Pure PCA gives a biased axis on helices
 * with multiple turns (sin(θ)·θ correlation is non-zero for finite
 * sample), so cross-product is preferred for this geometry. */
static int helix_axis(const liscio_pose_t *pts, int n, double n_out[3])
{
    if (n < 4) return -1;
    double sum[3] = {0, 0, 0};
    int valid = 0;
    for (int i = 0; i + 2 < n; i++) {
        double c1[3] = {
            pts[i+1].x - pts[i].x,
            pts[i+1].y - pts[i].y,
            pts[i+1].z - pts[i].z
        };
        double c2[3] = {
            pts[i+2].x - pts[i+1].x,
            pts[i+2].y - pts[i+1].y,
            pts[i+2].z - pts[i+1].z
        };
        double cr[3];
        v_cross(cr, c1, c2);
        double m = v_norm(cr);
        if (m < 1e-15) continue;          /* skip collinear triplets */
        sum[0] += cr[0]; sum[1] += cr[1]; sum[2] += cr[2];
        valid++;
    }
    if (valid == 0) return -1;
    if (v_normalize(sum) != 0) return -1;
    n_out[0] = sum[0]; n_out[1] = sum[1]; n_out[2] = sum[2];
    return 0;
}

/* Refine an existing axis estimate by minimizing planar roundness
 * residual on a few finite-difference axis perturbations.  One pass of
 * coordinate descent with two perpendicular tangent vectors (b1, b2)
 * is usually enough to recover from chord-cross-product's angular
 * error on sparse helices.  Returns refined axis in n_io. */
static double helix_round_dev(const liscio_pose_t *pts, int n,
                               const double n_axis[3])
{
    /* Build local frame perpendicular to axis. */
    double u[3], v[3];
    double z[3] = { 0, 0, 1 };
    double cu[3] = { n_axis[1]*z[2] - n_axis[2]*z[1],
                     n_axis[2]*z[0] - n_axis[0]*z[2],
                     n_axis[0]*z[1] - n_axis[1]*z[0] };
    if (v_norm(cu) < 1e-6) {
        double xax[3] = { 1, 0, 0 };
        cu[0] = n_axis[1]*xax[2] - n_axis[2]*xax[1];
        cu[1] = n_axis[2]*xax[0] - n_axis[0]*xax[2];
        cu[2] = n_axis[0]*xax[1] - n_axis[1]*xax[0];
    }
    v_normalize(cu);
    u[0] = cu[0]; u[1] = cu[1]; u[2] = cu[2];
    v[0] = n_axis[1]*u[2] - n_axis[2]*u[1];
    v[1] = n_axis[2]*u[0] - n_axis[0]*u[2];
    v[2] = n_axis[0]*u[1] - n_axis[1]*u[0];

    /* Project + Kasa-fit + measure max roundness deviation. */
    if (n > LISCIO_MAX_WINDOW) return 1e9;
    double xs[LISCIO_MAX_WINDOW], ys[LISCIO_MAX_WINDOW];
    for (int i = 0; i < n; i++) {
        double dx = pts[i].x - pts[0].x;
        double dy = pts[i].y - pts[0].y;
        double dz = pts[i].z - pts[0].z;
        xs[i] = dx*u[0] + dy*u[1] + dz*u[2];
        ys[i] = dx*v[0] + dy*v[1] + dz*v[2];
    }
    double cx2, cy2, r;
    if (kasa_circle(xs, ys, n, &cx2, &cy2, &r) != 0) return 1e9;
    double mrd = 0;
    for (int i = 0; i < n; i++) {
        double dx = xs[i] - cx2, dy = ys[i] - cy2;
        double dr = fabs(sqrt(dx*dx + dy*dy) - r);
        if (dr > mrd) mrd = dr;
    }
    return mrd;
}

/* Coordinate-descent axis refinement.  Tries small angular perturbations
 * of the axis along two perpendicular directions; keeps any reduction
 * in planar roundness deviation.  3 iterations is enough to recover
 * the chord-cross seed's ~5° bias on sparse helices. */
static void refine_axis(const liscio_pose_t *pts, int n, double axis[3])
{
    /* Build perpendicular tangent basis (b1, b2). */
    double b1[3], b2[3];
    double zw[3] = { 0, 0, 1 };
    double cu[3] = { axis[1]*zw[2] - axis[2]*zw[1],
                     axis[2]*zw[0] - axis[0]*zw[2],
                     axis[0]*zw[1] - axis[1]*zw[0] };
    if (v_norm(cu) < 1e-6) {
        double xw[3] = { 1, 0, 0 };
        cu[0] = axis[1]*xw[2] - axis[2]*xw[1];
        cu[1] = axis[2]*xw[0] - axis[0]*xw[2];
        cu[2] = axis[0]*xw[1] - axis[1]*xw[0];
    }
    v_normalize(cu);
    b1[0]=cu[0]; b1[1]=cu[1]; b1[2]=cu[2];
    b2[0] = axis[1]*b1[2] - axis[2]*b1[1];
    b2[1] = axis[2]*b1[0] - axis[0]*b1[2];
    b2[2] = axis[0]*b1[1] - axis[1]*b1[0];

    double step = 0.05;     /* radians */
    for (int iter = 0; iter < 12; iter++) {
        double base = helix_round_dev(pts, n, axis);
        int improved = 0;
        for (int dir = 0; dir < 4; dir++) {
            double sign = (dir & 1) ? -1.0 : 1.0;
            const double *b = (dir < 2) ? b1 : b2;
            double trial[3] = {
                axis[0] + sign * step * b[0],
                axis[1] + sign * step * b[1],
                axis[2] + sign * step * b[2]
            };
            v_normalize(trial);
            double dev = helix_round_dev(pts, n, trial);
            if (dev < base * 0.999) {
                axis[0] = trial[0]; axis[1] = trial[1]; axis[2] = trial[2];
                /* Re-orthogonalize tangents after axis update. */
                cu[0] = axis[1]*zw[2] - axis[2]*zw[1];
                cu[1] = axis[2]*zw[0] - axis[0]*zw[2];
                cu[2] = axis[0]*zw[1] - axis[1]*zw[0];
                if (v_norm(cu) < 1e-6) {
                    double xw[3] = { 1, 0, 0 };
                    cu[0] = axis[1]*xw[2] - axis[2]*xw[1];
                    cu[1] = axis[2]*xw[0] - axis[0]*xw[2];
                    cu[2] = axis[0]*xw[1] - axis[1]*xw[0];
                }
                v_normalize(cu);
                b1[0]=cu[0]; b1[1]=cu[1]; b1[2]=cu[2];
                b2[0] = axis[1]*b1[2] - axis[2]*b1[1];
                b2[1] = axis[2]*b1[0] - axis[0]*b1[2];
                b2[2] = axis[0]*b1[1] - axis[1]*b1[0];
                base = dev;
                improved = 1;
                break;
            }
        }
        if (!improved) step *= 0.5;     /* shrink step */
        if (step < 1e-6) break;
    }
}

/* Build orthonormal plane basis (u, v) given axis n.  Picks u ⟂ world
 * Z when possible; falls back to world X when n ≈ ±Z. */
static void plane_basis(const double n[3], double u[3], double v[3])
{
    double z[3] = { 0, 0, 1 };
    double cu[3];
    v_cross(cu, n, z);
    if (v_norm(cu) < 1e-6) {
        double xax[3] = { 1, 0, 0 };
        v_cross(cu, n, xax);
    }
    v_normalize(cu);
    u[0] = cu[0]; u[1] = cu[1]; u[2] = cu[2];
    double cv[3];
    v_cross(cv, n, u);
    v_normalize(cv);
    v[0] = cv[0]; v[1] = cv[1]; v[2] = cv[2];
}

/* ------------------------------------------------------------------ */
int liscio_helix9_fit(const struct liscio_ctx *ctx, int i0, int i1,
                       liscio_helix9_fit_t *out)
{
    if (!ctx || !out) return -1;
    int n = i1 - i0 + 1;
    /* Need at least 4 points: 3 to define a circle plus 1 for axial
     * regression.  Reuse arc-fit's min_arc_pts (typically 5). */
    if (n < ctx->cfg.min_arc_pts) return -1;
    if (n < 4) return -1;
    if (n > LISCIO_MAX_WINDOW) return -1;

    const liscio_pose_t *pts = &ctx->pts[i0];

    /* 1. Helix axis (chord-cross seed + coordinate-descent refinement). */
    double nrm[3];
    if (helix_axis(pts, n, nrm) != 0) return -1;
    refine_axis(pts, n, nrm);

    /* 2. Plane basis (u, v ⟂ n). */
    double u[3], v[3];
    plane_basis(nrm, u, v);

    /* 3. Project to local frame: q = (P - P0) · {u, v}, z = (P - P0) · n. */
    double xs[LISCIO_MAX_WINDOW], ys[LISCIO_MAX_WINDOW];
    double zs[LISCIO_MAX_WINDOW];
    double p0[3] = { pts[0].x, pts[0].y, pts[0].z };
    for (int i = 0; i < n; i++) {
        double d[3] = { pts[i].x - p0[0], pts[i].y - p0[1], pts[i].z - p0[2] };
        xs[i] = v_dot(d, u);
        ys[i] = v_dot(d, v);
        zs[i] = v_dot(d, nrm);
    }

    /* 4. Kasa-fit the planar circle. */
    double cx2, cy2, radius;
    if (kasa_circle(xs, ys, n, &cx2, &cy2, &radius) != 0) return -1;
    if (radius < ctx->cfg.min_arc_radius) return -1;
    if (radius > ctx->cfg.max_arc_radius) return -1;

    /* 5. Roundness check (planar).  Half-tolerance gate gives the
     * downstream verifier (which uses continuous-curve closest-point
     * distance) headroom against discrete-sample reconstruction
     * residual.  Axis refinement keeps this strict gate satisfiable
     * even on sparse helices. */
    const double helix_tol = ctx->cfg.tol_xyz * 0.5;
    for (int i = 0; i < n; i++) {
        double dx = xs[i] - cx2, dy = ys[i] - cy2;
        double dr = fabs(sqrt(dx*dx + dy*dy) - radius);
        if (dr > helix_tol) return -1;
    }

    /* 6. Unwrapped angles around the planar center. */
    double angs[LISCIO_MAX_WINDOW];
    angs[0] = atan2(ys[0] - cy2, xs[0] - cx2);
    for (int i = 1; i < n; i++) {
        double a = atan2(ys[i] - cy2, xs[i] - cx2);
        double da = a - angs[i-1];
        while (da >  M_PI) da -= 2*M_PI;
        while (da < -M_PI) da += 2*M_PI;
        angs[i] = angs[i-1] + da;
    }
    double arc_total = angs[n-1] - angs[0];
    double arc_len_total = fabs(arc_total) * radius;
    if (arc_len_total < 1e-9) return -1;

    /* 7. Linear regression z = z0 + p·θ.  Standard 1D LSQ. */
    double t_mean = 0, z_mean = 0;
    for (int i = 0; i < n; i++) { t_mean += angs[i]; z_mean += zs[i]; }
    t_mean /= n; z_mean /= n;
    double Stt = 0, Stz = 0;
    for (int i = 0; i < n; i++) {
        double dt = angs[i] - t_mean;
        Stt += dt * dt;
        Stz += dt * (zs[i] - z_mean);
    }
    if (Stt < 1e-15) return -1;
    double pitch = Stz / Stt;
    double z0 = z_mean - pitch * t_mean;

    /* 8. Reconstruct + 3D max-deviation gate. */
    double max_dev = 0.0;
    for (int i = 0; i < n; i++) {
        double cs = cos(angs[i]), sn = sin(angs[i]);
        double zr = z0 + pitch * angs[i];
        double pred[3] = {
            p0[0] + (cx2 + radius * cs) * u[0]
                  + (cy2 + radius * sn) * v[0]
                  + zr * nrm[0],
            p0[1] + (cx2 + radius * cs) * u[1]
                  + (cy2 + radius * sn) * v[1]
                  + zr * nrm[1],
            p0[2] + (cx2 + radius * cs) * u[2]
                  + (cy2 + radius * sn) * v[2]
                  + zr * nrm[2]
        };
        double dx = pts[i].x - pred[0];
        double dy = pts[i].y - pred[1];
        double dz = pts[i].z - pred[2];
        double dev = sqrt(dx*dx + dy*dy + dz*dz);
        if (dev > max_dev) max_dev = dev;
        if (max_dev > helix_tol) return -1;
    }

    /* 9. ABC/UVW linearity along arc length.  Same idea as arc9d. */
    #define LINCHECK(FIELD, TOL)                                            \
    do {                                                                    \
        double v0 = pts[0].FIELD, vN = pts[n-1].FIELD;                     \
        for (int i = 1; i < n - 1; i++) {                                  \
            double t = (angs[i] - angs[0]) / arc_total;                    \
            double predicted = v0 + t * (vN - v0);                         \
            if (fabs(pts[i].FIELD - predicted) > (TOL)) return -1;         \
        }                                                                   \
    } while (0)

    LINCHECK(a, ctx->cfg.tol_abc);
    LINCHECK(b, ctx->cfg.tol_abc);
    LINCHECK(c, ctx->cfg.tol_abc);
    LINCHECK(u, ctx->cfg.tol_uvw);
    LINCHECK(v, ctx->cfg.tol_uvw);
    LINCHECK(w, ctx->cfg.tol_uvw);
    #undef LINCHECK

    /* 10. Center in world coords: P0 + cx2·u + cy2·v + z0·n.  This is
     *     the point on the axis at the start's plane (z0 axial offset
     *     puts the start exactly on the rim). */
    out->cx = p0[0] + cx2 * u[0] + cy2 * v[0] + z0 * nrm[0];
    out->cy = p0[1] + cx2 * u[1] + cy2 * v[1] + z0 * nrm[1];
    out->cz = p0[2] + cx2 * u[2] + cy2 * v[2] + z0 * nrm[2];
    out->nx = nrm[0]; out->ny = nrm[1]; out->nz = nrm[2];
    out->radius        = radius;
    out->arc_angle     = arc_total;
    out->pitch         = pitch;
    out->max_deviation = max_dev;
    return 0;
}
