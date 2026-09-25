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
* Description: liscio private internal types and fitter declarations.
*              NOT public API; do not include from consumers.
* Author:      杨阳 (Yang Yang) <mika-net@outlook.com>
* License:     MIT (SPDX-License-Identifier: MIT)
* Copyright (c) 2026 杨阳 (Yang Yang)
********************************************************************/

#ifndef LISCIO_INTERNAL_H
#define LISCIO_INTERNAL_H

#include "liscio/liscio.h"

#define LISCIO_MAX_WINDOW 128

struct liscio_ctx {
    liscio_cfg_t    cfg;
    liscio_emit_cb  cb;
    void            *cb_user;

    /* Waypoint window (ring? no — linear, flushed when full). */
    liscio_pose_t   pts[LISCIO_MAX_WINDOW];
    double           feed[LISCIO_MAX_WINDOW];
    int              tag[LISCIO_MAX_WINDOW];
    int              n_pts;

    /* Cross-flush G1 continuity: tangent at the most recently emitted
     * primitive's end (forward direction).  Used to seed the next
     * window's first Bezier fit so the join is G1.  Cleared on
     * corner-induced flushes (where discontinuity is desired). */
    double           prev_emit_tan_x, prev_emit_tan_y, prev_emit_tan_z;
    int              have_prev_emit_tan;

    /* Cross-flush G2 continuity: curvature (kappa) + curve-normal at the most
     * recently emitted CURVATURE-MATCHED quintic's END.  Seeds the next
     * window's g2 quintic k_start so consecutive curved primitives share the
     * seam curvature (C2) instead of each self-estimating κ → no curvature step
     * → no velocity hunting.  Set ONLY right after a successful g2 quintic;
     * cleared by every other emit (via cache_exit_tangent) and at corners /
     * events — so the G2 chain runs through a smooth curved run and BREAKS at
     * transitions, where the RT blend engine (bezier9, curvature-matched) makes
     * the join G2 instead.  This is the user's "keep G2, split where you can't,
     * let blend connect" design. */
    double           prev_emit_kappa;
    double           prev_emit_norm_x, prev_emit_norm_y, prev_emit_norm_z;
    int              have_prev_emit_kappa;

    /* Micro-cusp cluster collapse state (cfg.collapse_cusp_enable).  When a
     * hard corner between two sub-tolerance segments is detected, enter
     * collapse mode: subsequent points within collapse_cusp_tol of
     * cluster_anchor are decimated (dropped) so a 2 µm thrash cluster emits
     * ONE stop junction instead of N.  Cleared when a point escapes the ball
     * or on flush/reset. */
    int              in_cluster;
    liscio_pose_t    cluster_anchor;

    /* Pending arc buffer (arc-merge).  When cfg.arc_merge_mode
     * != OFF, liscio_add_arc holds the most recent arc here and only
     * emits when the next call cannot extend it.  pending_feed/tag
     * carry the source-line metadata of the FIRST arc in the run. */
    liscio_primitive_t pending_arc;
    double             pending_feed;
    int                pending_tag;          /* first absorbed arc's line */
    int                pending_tag_last;     /* last absorbed arc's line  */
    int                has_pending_arc;
    int                pending_n_absorbed;

    /* Arc-helix buffer: when cfg.arc_to_helix_max > 0, liscio_add_arc
     * accumulates G2/G3 arcs whose centers, axes and radii are
     * approximately consistent, then flushes them as a single HELIX
     * primitive after testing the buffered set against the helix
     * geometric model (axis line through all centers, constant pitch,
     * constant radius).  Falls back to per-arc emission on failure. */
#define LISCIO_ARC_HELIX_BUF 64
    liscio_primitive_t arc_buf[LISCIO_ARC_HELIX_BUF];
    int                arc_buf_n;
    int                arc_buf_tag_first;
    int                arc_buf_tag_last;
    double             arc_buf_feed;

    /* Arc-chain G2 fusion pending buffer (cfg.arc_g2fuse).  liscio_add_arc
     * holds the most recent arc so the NEXT arc's start tangent is known when
     * it is emitted: a soft kink (seam tangent break ≤ corner_angle_deg —
     * thomam has 26 of 1-19°) is then ABSORBED into the quintics by averaging
     * the seam tangent on both sides (the fit tolerance gate rejects kinks too
     * big to absorb → bisect → native).  Without this each kink seam needs a
     * corner blend whose dκ/ds-dominated velocity cap (~21 mm/s) wrecks feed
     * holding.  arcg2_pend_t1 carries the averaged tangent override for the
     * PENDING arc's own start (set when the PREVIOUS seam was soft-kink). */
    liscio_primitive_t arcg2_pend;
    double             arcg2_pend_feed;
    int                arcg2_pend_tag;
    int                has_arcg2_pend;
    double             arcg2_pend_t1[3];
    int                arcg2_pend_has_t1;
    int                emitting_arc_g2;      /* emit_bspline is an arc-chain fuse */

    liscio_stats_t  stats;
};

/* Arc fit result — filled by liscio_arc9d_fit on success. */
typedef struct {
    double cx, cy, cz;    /* center */
    double nx, ny, nz;    /* unit normal of plane */
    double radius;
    double arc_angle;     /* signed */
    double max_deviation; /* xyz (ref plane) */
} liscio_arc9d_fit_t;

/* Cubic Bezier 9D fit result. */
typedef struct {
    liscio_pose_t P0, P1, P2, P3;
    double total_chord;
    double max_deviation;
} liscio_bezier9_fit_t;

/* Attempt to fit a 9D arc through pts[i0..i1] (inclusive). Returns 0
 * on success (all points within tol_xyz of arc), -1 on failure. */
int liscio_arc9d_fit(const struct liscio_ctx *ctx, int i0, int i1,
                      liscio_arc9d_fit_t *out);

/* Attempt to fit a 9D cubic Bezier through pts[i0..i1] via LSQ. */
int liscio_bezier9_fit(const struct liscio_ctx *ctx, int i0, int i1,
                        liscio_bezier9_fit_t *out);

/* G1-continuous variant of liscio_bezier9_fit.
 *   tHat1: unit tangent at left end (P0), pointing FORWARD along curve.
 *   tHat2: unit tangent at right end (P3), pointing INWARD (toward P2) —
 *          Schneider Graphics-Gems convention.  This means at t=1 the
 *          curve tangent direction is -tHat2.
 * XYZ control points are constrained: P1 = P0 + αL·tHat1,
 *                                      P2 = P3 + αR·tHat2.
 * The 2-DoF system (αL, αR) is solved via 2×2 normal equations,
 * iterated with Hoschek reparameterization.  ABC/UVW dimensions are
 * fitted unconstrained per-axis with the refined parameter values.
 *
 * Use this when both endpoint tangents are dictated by neighbouring
 * primitives (composite recursion sharing a split-point tangent, or
 * matching a previous primitive's exit tangent across a flush).
 *
 * Returns 0 on success (fit within tol_xyz / tol_abc / tol_uvw),
 * -1 on tolerance miss or singular system. */
int liscio_bezier9_fit_g1(const struct liscio_ctx *ctx, int i0, int i1,
                          double tHat1_x, double tHat1_y, double tHat1_z,
                          double tHat2_x, double tHat2_y, double tHat2_z,
                          liscio_bezier9_fit_t *out);

/* Composite Bezier: recursively split until each sub-range fits.
 * emit_fn is called once per successful sub-fit.
 * Returns number of emitted primitives (>=0), -1 on error. */
typedef void (*liscio_bezier9_emit_fn)(int i0, int i1,
    const liscio_bezier9_fit_t *fit, void *user);

int liscio_bezier9_composite_fit(const struct liscio_ctx *ctx,
                                  int i0, int i1,
                                  liscio_bezier9_emit_fn emit,
                                  void *user);

/* Same as composite_fit, but takes an explicit forward-direction tangent
 * for the LEFT outermost endpoint.  Used to enforce G1 across flush
 * boundaries when the previous primitive's exit tangent is known.
 * Pass NULL for left_tan_xyz to fall back to chord-derived tangent. */
int liscio_bezier9_composite_fit_g1(const struct liscio_ctx *ctx,
                                     int i0, int i1,
                                     const double *left_tan_xyz,
                                     liscio_bezier9_emit_fn emit,
                                     void *user);

/* Helix 9D fit result.  A helix is parameterized as
 *   P(θ) = c + r·(cos(θ)·u + sin(θ)·v) + (z0 + p·θ)·n
 * where (u,v,n) is an orthonormal frame with n the helix axis, r the
 * planar radius, p the pitch (mm/rad; 2π·p = mm per turn), z0 axial
 * offset of the start point.  Reduces to a planar arc when p == 0.
 *
 * Center (cx,cy,cz) = c + z0·n is the projection of the start point
 * onto a plane perpendicular to the axis through the helix center —
 * matches arc9d emit semantics (start point is on the rim, center on
 * the axis at start's plane). */
typedef struct {
    double cx, cy, cz;        /* center on axis at start plane */
    double nx, ny, nz;        /* axis direction (unit) */
    double radius;
    double arc_angle;         /* signed total angle (rad) */
    double pitch;             /* mm/rad along +n */
    double max_deviation;     /* worst |raw - reconstructed| (mm) */
} liscio_helix9_fit_t;

/* Attempt to fit a 9D helix through pts[i0..i1] (inclusive).  Returns 0
 * on success (max_dev ≤ tol_xyz) and fills out; -1 on geometry failure
 * or tolerance miss.  Caller should try arc9d_fit first; this is the
 * fallback when the points don't lie on a plane (i.e. helical climb). */
int liscio_helix9_fit(const struct liscio_ctx *ctx, int i0, int i1,
                       liscio_helix9_fit_t *out);

/* Cubic B-spline 9D fit result. */
typedef struct {
    int           degree;   /* 3 */
    int           n_ctrl;
    liscio_pose_t ctrl[LISCIO_BSPLINE_MAX_CTRL];
    double        weights[LISCIO_BSPLINE_MAX_CTRL];
    double        knots[LISCIO_BSPLINE_MAX_CTRL + 4];
    double        total_chord;
    double        max_deviation;
    /* Endpoint curvature at s=1 (filled by the g2 quintic fit) — fed into the
     * ctx G2 chain so the NEXT curve's k_start matches this curve's k_end. */
    double        end_kappa;
    double        end_norm_x, end_norm_y, end_norm_z;
} liscio_bspline9_fit_t;

/* κ=0 quintic (G2 at bezier↔line seams) — degree-5 spline out. */
int liscio_quintic9_g2zero_fit(const struct liscio_ctx *ctx, int i0, int i1,
    double t1x, double t1y, double t1z, double t2x, double t2y, double t2z,
    liscio_bspline9_fit_t *out);

/* Generalized quintic G2 fit: pins ARBITRARY endpoint curvature (κ≠0), estimated
 * from the endpoint osculating circles — C2 curvature-continuous seams for arc /
 * curved-micro-line runs (vs the κ=0 g2zero which is only correct at line seams). */
int liscio_quintic9_g2_fit(const struct liscio_ctx *ctx, int i0, int i1,
    double t1x, double t1y, double t1z, double t2x, double t2y, double t2z,
    liscio_bspline9_fit_t *out);

/* Attempt to fit a 9D cubic B-spline through pts[i0..i1] with M
 * control points (4 ≤ M ≤ LISCIO_BSPLINE_MAX_CTRL).  Clamped uniform
 * knot vector; chord-length parameterization; LSQ normal equations
 * via SPD Cholesky. */
int liscio_bspline9_fit(const struct liscio_ctx *ctx, int i0, int i1,
                         int M, liscio_bspline9_fit_t *out);

#endif /* LISCIO_INTERNAL_H */
