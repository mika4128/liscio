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
* Description: liscio coordinator + window state machine.  Drives
*              the waypoint-accumulation → classify → fit → emit
*              pipeline.  Strategy (most-specific first):
*                1. LINE if all points colinear
*                2. ARC (9D Kasa LSQ, ≥ min_arc_pts)
*                3. BEZIER cubic (single LSQ, then composite split)
*                4. LINE fallback
*              Streaming corner detection flushes window at tangent
*              discontinuity (Siemens G641/642 equivalent).
* Author:      杨阳 (Yang Yang) <mika-net@outlook.com>
* License:     MIT (SPDX-License-Identifier: MIT)
* Copyright (c) 2026 杨阳 (Yang Yang)
********************************************************************/

#define _USE_MATH_DEFINES
#include "liscio/liscio.h"
#include "liscio_internal.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- DIAGNOSTIC instrumentation (#2 dense-crawl investigation) ----
 * Unconditional globals (negligible cost); read by test_fit_stats to see
 * WHY windows flush + which emit paths fire.  STRIP before final commit. */
struct liscio_dbg_counters {
    long flush_corner, flush_disc, flush_maxwin, flush_final;
    long flush_feedchg;                /* same-F gate: F change split the window */
    long emit_line_collin, emit_line_fallback;
    long bez_g1_used, bez_g1_miss;   /* cross-window G1 pin on first bezier */
    long stream_absorb;              /* add_line collinear extend (no new pt) */
};
struct liscio_dbg_counters liscio_dbg = {0};
void liscio_dbg_reset(void) { memset(&liscio_dbg, 0, sizeof(liscio_dbg)); }

/* "Same F" for merge gating (cfg.merge_same_feed_only): non-positive is a
 * wildcard (unset); otherwise require agreement up to fp noise (rel 1e-9,
 * same yardstick as the streaming collinear absorb's feed match).
 * TP2_LISCIO_NOFGATE=1 bypasses the gate entirely — diagnostic A/B switch
 * for F-change window splits (process-cached: relaunch to flip). */
static int feed_match(double fa, double fb)
{
    static int nogate = -1;
    if (nogate < 0) nogate = getenv("TP2_LISCIO_NOFGATE") ? 1 : 0;
    if (nogate) return 1;
    if (fa <= 0.0 || fb <= 0.0) return 1;
    return fabs(fa - fb) < 1e-9 * fmax(1.0, fmax(fabs(fa), fabs(fb)));
}

/* ------------------------------------------------------------------ */
void liscio_cfg_default(liscio_cfg_t *cfg)
{
    if (!cfg) return;
    cfg->tol_xyz      = 0.01;      /* 10 μm */
    cfg->tol_abc      = 0.1;       /* 0.1 deg */
    cfg->tol_uvw      = -1.0;      /* follow xyz */
    cfg->tol_arcfit   = 0.0;       /* <=0 follows tol_xyz (fusion sets its own cap) */
    cfg->max_window   = LISCIO_MAX_WINDOW;
    cfg->min_arc_pts  = 5;
    cfg->min_arc_radius = 1e-3;    /* 1 μm */
    cfg->max_arc_radius = 1e6;     /* 1 km — beyond is effectively linear */
    cfg->corner_angle_deg = 15.0;  /* 15° sharp-corner threshold */
    cfg->collinear_tol   = -1.0;   /* <=0 follows tol_xyz */
    cfg->enable_bezier   = 1;
    cfg->enable_composite_bezier = 1;
    cfg->enable_bspline  = 1;
    cfg->enable_g2_quintic = 0;
    cfg->enable_g2_curvature = 0;
    cfg->bspline_init_ctrl = 6;
    cfg->arc_subseg_samples = 0;   /* 0 = passthrough (default lossless) */
    cfg->arc_merge_mode     = 0;   /* LISCIO_ARC_MERGE_OFF (default) */
    cfg->arc_merge_xy_tol   = 1e-6;
    cfg->arc_merge_pitch_tol = 1e-6;
    cfg->arc_to_helix_max    = 16;      /* on by default; set 0 to disable.
                                         * Auto-detects helical G2/G3 runs
                                         * (e.g. lathe threads) and folds
                                         * them into LISCIO_PRIM_HELIX with
                                         * tol_xyz residual gate. */
    cfg->corner_detection_mode = 0;     /* LISCIO_CORNER_IMMEDIATE (default) */
    cfg->corner_hard_angle_deg = 60.0;
    cfg->enable_g1_chain       = 0;     /* legacy off; #2 dense-crawl opt-in */
    cfg->prefer_bspline        = 0;     /* legacy off; #2 dense-crawl opt-in */
    cfg->collapse_cusp_enable  = 1;     /* on: fold 2µm hard-corner thrash → 1 stop */
    cfg->collapse_cusp_tol     = -1.0;  /* <=0 follows tol_xyz */
    cfg->arc_g2fuse            = 0;     /* off: native arc passthrough (default) */
    cfg->merge_same_feed_only  = 1;     /* on: an F change splits the window */
}

void liscio_cfg_set_corner_detection(liscio_cfg_t *cfg,
                                       liscio_corner_mode_t mode,
                                       double soft_angle_deg,
                                       double hard_angle_deg)
{
    if (!cfg) return;
    cfg->corner_detection_mode = (int)mode;
    if (soft_angle_deg > 0) cfg->corner_angle_deg = soft_angle_deg;
    if (hard_angle_deg > 0) cfg->corner_hard_angle_deg = hard_angle_deg;
    /* Sanity: hard ≥ soft. */
    if (cfg->corner_hard_angle_deg < cfg->corner_angle_deg)
        cfg->corner_hard_angle_deg = cfg->corner_angle_deg;
}

/* Public API setter — single entry to opt into arc merge with sensible
 * tolerance defaults.  Pass <=0 for either tol to keep its default. */
void liscio_cfg_set_arc_merge(liscio_cfg_t *cfg,
                                liscio_arc_merge_mode_t mode,
                                double xy_tol,
                                double pitch_tol)
{
    if (!cfg) return;
    cfg->arc_merge_mode = (int)mode;
    if (xy_tol    > 0) cfg->arc_merge_xy_tol    = xy_tol;
    else if (cfg->arc_merge_xy_tol <= 0) cfg->arc_merge_xy_tol = 1e-6;
    if (pitch_tol > 0) cfg->arc_merge_pitch_tol = pitch_tol;
    else if (cfg->arc_merge_pitch_tol <= 0) cfg->arc_merge_pitch_tol = 1e-6;
}

/* ------------------------------------------------------------------ */
liscio_ctx_t *liscio_create(const liscio_cfg_t *cfg)
{
    liscio_ctx_t *ctx = (liscio_ctx_t *)calloc(1, sizeof(*ctx));
    if (!ctx) return NULL;

    if (cfg) {
        ctx->cfg = *cfg;
    } else {
        liscio_cfg_default(&ctx->cfg);
    }
    if (ctx->cfg.max_window > LISCIO_MAX_WINDOW)
        ctx->cfg.max_window = LISCIO_MAX_WINDOW;
    if (ctx->cfg.max_window < 2)
        ctx->cfg.max_window = 2;
    if (ctx->cfg.tol_uvw <= 0.0)
        ctx->cfg.tol_uvw = ctx->cfg.tol_xyz;
    if (ctx->cfg.collinear_tol <= 0.0)
        ctx->cfg.collinear_tol = ctx->cfg.tol_xyz;

    ctx->n_pts = 0;
    ctx->cb = NULL;
    ctx->cb_user = NULL;
    return ctx;
}

/* ------------------------------------------------------------------ */
void liscio_destroy(liscio_ctx_t *ctx)
{
    if (!ctx) return;
    free(ctx);
}

/* ------------------------------------------------------------------ */
void liscio_set_callback(liscio_ctx_t *ctx, liscio_emit_cb cb, void *user)
{
    if (!ctx) return;
    ctx->cb = cb;
    ctx->cb_user = user;
}

/* ------------------------------------------------------------------ */
void liscio_get_stats(const liscio_ctx_t *ctx, liscio_stats_t *out)
{
    if (!ctx || !out) return;
    *out = ctx->stats;
}

/* Fill tangent vectors + default weights/tolerance on a primitive. */
static void fill_common_fields(liscio_ctx_t *ctx, liscio_primitive_t *prim)
{
    prim->w0 = prim->w1 = prim->w2 = prim->w3 = 1.0;
    prim->tol_xyz = ctx->cfg.tol_xyz;
    prim->tol_abc = ctx->cfg.tol_abc;
}

/* Cache the most recent primitive's exit tangent so the next window's
 * first Bezier fit can pin its start tangent for cross-flush G1
 * continuity.  Must be called by every emit_* path after ctx->cb. */
static void cache_exit_tangent(liscio_ctx_t *ctx, const liscio_primitive_t *p)
{
    /* Break the G2 curvature chain by default: it is re-armed ONLY by the g2
     * quintic branch right after this call.  So any non-quintic emit (line,
     * cubic, arc, spline) ends the chain → its seam becomes a blend join. */
    ctx->have_prev_emit_kappa = 0;
    double m = sqrt(p->tan_end_x*p->tan_end_x +
                    p->tan_end_y*p->tan_end_y +
                    p->tan_end_z*p->tan_end_z);
    if (m < 1e-12) { ctx->have_prev_emit_tan = 0; return; }
    ctx->prev_emit_tan_x = p->tan_end_x / m;
    ctx->prev_emit_tan_y = p->tan_end_y / m;
    ctx->prev_emit_tan_z = p->tan_end_z / m;
    ctx->have_prev_emit_tan = 1;
}

/* No-op shim retained while internal call sites are migrated to the
 * G-code-style RAPID/STOP event primitives.  Kept tiny so it's free at
 * -O2; will be removed when bridges + tests have been updated. */
static inline void consume_run_marker(liscio_ctx_t *ctx,
                                      liscio_primitive_t *p)
{
    (void)ctx; (void)p;
}

static void fill_tangent_from_points(liscio_primitive_t *prim,
                                      const liscio_pose_t *s,
                                      const liscio_pose_t *e)
{
    double dx = e->x - s->x, dy = e->y - s->y, dz = e->z - s->z;
    double m = sqrt(dx*dx + dy*dy + dz*dz);
    if (m < 1e-15) {
        prim->tan_start_x = prim->tan_end_x = 1.0;
        prim->tan_start_y = prim->tan_end_y = 0.0;
        prim->tan_start_z = prim->tan_end_z = 0.0;
    } else {
        prim->tan_start_x = prim->tan_end_x = dx / m;
        prim->tan_start_y = prim->tan_end_y = dy / m;
        prim->tan_start_z = prim->tan_end_z = dz / m;
    }
}

/* ------------------------------------------------------------------ */
/* Emit a LINE primitive spanning waypoints[i0..i1] (inclusive).      */
static void emit_line(liscio_ctx_t *ctx, int i0, int i1)
{
    if (!ctx->cb) return;
    if (i0 < 0 || i1 < i0 || i1 >= ctx->n_pts) return;

    liscio_primitive_t prim;
    memset(&prim, 0, sizeof(prim));
    prim.type         = LISCIO_PRIM_LINE;
    prim.start        = ctx->pts[i0];
    prim.end          = ctx->pts[i1];
    prim.feedrate       = ctx->feed[i0];
    prim.tag_line_first = ctx->tag[i0];
    prim.tag_line_last  = ctx->tag[i1];
    prim.tag_line_id    = prim.tag_line_first;
    prim.n_absorbed     = i1 - i0;  /* # segments merged (= # waypoints − 1) */
    fill_common_fields(ctx, &prim);
    fill_tangent_from_points(&prim, &prim.start, &prim.end);

    consume_run_marker(ctx, &prim);
    ctx->cb(&prim, ctx->cb_user);
    cache_exit_tangent(ctx, &prim);
    ctx->stats.emitted_line++;
    ctx->stats.absorbed_total += prim.n_absorbed;
}

/* Emit a BEZIER primitive. */
static void emit_bezier(liscio_ctx_t *ctx, int i0, int i1,
                        const liscio_bezier9_fit_t *fit);

/* Emit a SPLINE (cubic B-spline) primitive. */
static void emit_bspline(liscio_ctx_t *ctx, int i0, int i1,
                         const liscio_bspline9_fit_t *fit)
{
    if (!ctx->cb) return;
    if (!fit || i0 < 0 || i1 < i0 || i1 >= ctx->n_pts) return;

    liscio_primitive_t prim;
    memset(&prim, 0, sizeof(prim));
    prim.type                 = LISCIO_PRIM_SPLINE;
    prim.start                = ctx->pts[i0];
    prim.end                  = ctx->pts[i1];
    prim.degree               = fit->degree;
    prim.n_ctrl               = fit->n_ctrl;
    for (int j = 0; j < fit->n_ctrl; j++) {
        prim.ctrl_pts[j] = fit->ctrl[j];
        prim.weights[j]  = fit->weights[j];
    }
    memcpy(prim.knots, fit->knots,
           (size_t)(fit->n_ctrl + fit->degree + 1) * sizeof(double));
    prim.spline_max_deviation = fit->max_deviation;
    prim.feedrate             = ctx->feed[i0];
    prim.tag_line_first       = ctx->tag[i0];
    prim.tag_line_last        = ctx->tag[i1];
    prim.tag_line_id          = prim.tag_line_first;
    prim.n_absorbed           = i1 - i0;
    prim.from_arc_g2          = ctx->emitting_arc_g2 ? 1 : 0;

    /* Defaults: rational weights = 1, tol from cfg. */
    prim.w0 = prim.w1 = prim.w2 = prim.w3 = 1.0;
    prim.tol_xyz = ctx->cfg.tol_xyz;
    prim.tol_abc = ctx->cfg.tol_abc;

    /* Tangent from first/last control polygon segment. */
    double sx = fit->ctrl[1].x - fit->ctrl[0].x;
    double sy = fit->ctrl[1].y - fit->ctrl[0].y;
    double sz = fit->ctrl[1].z - fit->ctrl[0].z;
    double sm = sqrt(sx*sx + sy*sy + sz*sz);
    if (sm > 1e-15) {
        prim.tan_start_x = sx/sm;
        prim.tan_start_y = sy/sm;
        prim.tan_start_z = sz/sm;
    }
    int L = fit->n_ctrl - 1;
    double ex = fit->ctrl[L].x - fit->ctrl[L-1].x;
    double ey = fit->ctrl[L].y - fit->ctrl[L-1].y;
    double ez = fit->ctrl[L].z - fit->ctrl[L-1].z;
    double em = sqrt(ex*ex + ey*ey + ez*ez);
    if (em > 1e-15) {
        prim.tan_end_x = ex/em;
        prim.tan_end_y = ey/em;
        prim.tan_end_z = ez/em;
    }

    consume_run_marker(ctx, &prim);
    ctx->cb(&prim, ctx->cb_user);
    cache_exit_tangent(ctx, &prim);
    ctx->stats.emitted_spline++;
    ctx->stats.absorbed_total += prim.n_absorbed;
}

/* Adapter for composite recursion: translates (i0,i1,fit) → emit_bezier. */
static void bezier_composite_emit_adapter(int i0, int i1,
                                           const liscio_bezier9_fit_t *fit,
                                           void *user)
{
    liscio_ctx_t *ctx = (liscio_ctx_t *)user;
    emit_bezier(ctx, i0, i1, fit);
}

static void emit_bezier(liscio_ctx_t *ctx, int i0, int i1,
                        const liscio_bezier9_fit_t *fit)
{
    if (!ctx->cb) return;
    if (!fit || i0 < 0 || i1 < i0 || i1 >= ctx->n_pts) return;

    /* G2 upgrade: try a κ=0 quintic that leaves/re-enters the neighbour straight
     * segments with ZERO curvature, killing the centripetal-jerk step at the
     * bezier↔line seam. The κ=0 constraint self-selects: it only meets tolerance
     * where the data genuinely flattens at both ends (line-adjacent); a curve↔
     * curve seam (κ≠0 at the end) busts tolerance and falls through to the cubic.
     * Tangent directions come from the cubic fit, so G1 is preserved either way. */
    /* C2 upgrade: a curvature-MATCHED quintic (κ≠0) — for curve/arc seams where
     * the κ=0 quintic busts tolerance.  Pins the actual endpoint curvature
     * (estimated from the osculating circles) so the fit follows a genuinely
     * curved run without the κ=0 unbend-rebend ripple.  Tried FIRST (when
     * enabled) since it is the general case; g2zero remains for pure line seams. */
    if (ctx->cfg.enable_g2_curvature) {
        double t1x=fit->P1.x-fit->P0.x, t1y=fit->P1.y-fit->P0.y, t1z=fit->P1.z-fit->P0.z;
        double t2x=fit->P2.x-fit->P3.x, t2y=fit->P2.y-fit->P3.y, t2z=fit->P2.z-fit->P3.z;
        liscio_bspline9_fit_t q;
        if (liscio_quintic9_g2_fit(ctx, i0, i1, t1x,t1y,t1z, t2x,t2y,t2z, &q) == 0) {
            emit_bspline(ctx, i0, i1, &q);   /* clears the G2 chain via cache_exit_tangent */
            /* Re-arm the G2 chain: the NEXT window's g2 fit will pin its
             * k_start to THIS curve's k_end → shared seam curvature = C2. */
            ctx->prev_emit_kappa   = q.end_kappa;
            ctx->prev_emit_norm_x  = q.end_norm_x;
            ctx->prev_emit_norm_y  = q.end_norm_y;
            ctx->prev_emit_norm_z  = q.end_norm_z;
            ctx->have_prev_emit_kappa = 1;
            ctx->stats.emitted_g2_quintic++;
            return;
        }
    }

    if (ctx->cfg.enable_g2_quintic) {
        double t1x=fit->P1.x-fit->P0.x, t1y=fit->P1.y-fit->P0.y, t1z=fit->P1.z-fit->P0.z;
        double t2x=fit->P2.x-fit->P3.x, t2y=fit->P2.y-fit->P3.y, t2z=fit->P2.z-fit->P3.z;
        liscio_bspline9_fit_t q;
        if (liscio_quintic9_g2zero_fit(ctx, i0, i1, t1x,t1y,t1z, t2x,t2y,t2z, &q) == 0) {
            emit_bspline(ctx, i0, i1, &q);
            ctx->stats.emitted_g2_quintic++;
            return;
        }
    }

    liscio_primitive_t prim;
    memset(&prim, 0, sizeof(prim));
    prim.type         = LISCIO_PRIM_BEZIER;
    prim.start        = ctx->pts[i0];
    prim.end          = ctx->pts[i1];
    prim.P0           = fit->P0;
    prim.P1           = fit->P1;
    prim.P2           = fit->P2;
    prim.P3           = fit->P3;
    prim.bezier_max_deviation = fit->max_deviation;
    prim.feedrate       = ctx->feed[i0];
    prim.tag_line_first = ctx->tag[i0];
    prim.tag_line_last  = ctx->tag[i1];
    prim.tag_line_id    = prim.tag_line_first;
    prim.n_absorbed     = i1 - i0;
    fill_common_fields(ctx, &prim);

    /* Tangent = C'(0) direction = 3·(P1-P0); end = 3·(P3-P2). */
    {
        double dx = fit->P1.x - fit->P0.x;
        double dy = fit->P1.y - fit->P0.y;
        double dz = fit->P1.z - fit->P0.z;
        double m = sqrt(dx*dx + dy*dy + dz*dz);
        if (m < 1e-15) {
            /* Degenerate: fall back to chord. */
            fill_tangent_from_points(&prim, &prim.start, &prim.end);
        } else {
            prim.tan_start_x = dx / m;
            prim.tan_start_y = dy / m;
            prim.tan_start_z = dz / m;
        }
        double ex = fit->P3.x - fit->P2.x;
        double ey = fit->P3.y - fit->P2.y;
        double ez = fit->P3.z - fit->P2.z;
        double em = sqrt(ex*ex + ey*ey + ez*ez);
        if (em < 1e-15) {
            /* Use chord direction. */
            double cx = prim.end.x - prim.start.x;
            double cy = prim.end.y - prim.start.y;
            double cz = prim.end.z - prim.start.z;
            double cm = sqrt(cx*cx + cy*cy + cz*cz);
            if (cm > 1e-15) {
                prim.tan_end_x = cx / cm;
                prim.tan_end_y = cy / cm;
                prim.tan_end_z = cz / cm;
            }
        } else {
            prim.tan_end_x = ex / em;
            prim.tan_end_y = ey / em;
            prim.tan_end_z = ez / em;
        }
    }

    consume_run_marker(ctx, &prim);
    ctx->cb(&prim, ctx->cb_user);
    cache_exit_tangent(ctx, &prim);
    ctx->stats.emitted_spline++;  /* reuse spline counter for now */
    ctx->stats.absorbed_total += prim.n_absorbed;
}

/* Emit an ARC primitive. Caller supplies fitted parameters. */
static void emit_arc(liscio_ctx_t *ctx, int i0, int i1,
                     const liscio_arc9d_fit_t *fit)
{
    if (!ctx->cb) return;
    if (!fit || i0 < 0 || i1 < i0 || i1 >= ctx->n_pts) return;

    liscio_primitive_t prim;
    memset(&prim, 0, sizeof(prim));
    prim.type         = LISCIO_PRIM_ARC;
    prim.start        = ctx->pts[i0];
    prim.end          = ctx->pts[i1];
    prim.cx           = fit->cx;
    prim.cy           = fit->cy;
    prim.cz           = fit->cz;
    prim.nx           = fit->nx;
    prim.ny           = fit->ny;
    prim.nz           = fit->nz;
    prim.radius         = fit->radius;
    prim.arc_angle      = fit->arc_angle;
    prim.feedrate       = ctx->feed[i0];
    prim.tag_line_first = ctx->tag[i0];
    prim.tag_line_last  = ctx->tag[i1];
    prim.tag_line_id    = prim.tag_line_first;
    prim.n_absorbed     = i1 - i0;
    fill_common_fields(ctx, &prim);

    /* Tangent at arc endpoints: cross(normal, radius_vector) with
     * direction sign from arc_angle. */
    {
        double rs_x = prim.start.x - fit->cx;
        double rs_y = prim.start.y - fit->cy;
        double rs_z = prim.start.z - fit->cz;
        double re_x = prim.end.x - fit->cx;
        double re_y = prim.end.y - fit->cy;
        double re_z = prim.end.z - fit->cz;
        double sign = (fit->arc_angle >= 0.0) ? 1.0 : -1.0;
        double ts_x = sign * (fit->ny*rs_z - fit->nz*rs_y);
        double ts_y = sign * (fit->nz*rs_x - fit->nx*rs_z);
        double ts_z = sign * (fit->nx*rs_y - fit->ny*rs_x);
        double te_x = sign * (fit->ny*re_z - fit->nz*re_y);
        double te_y = sign * (fit->nz*re_x - fit->nx*re_z);
        double te_z = sign * (fit->nx*re_y - fit->ny*re_x);
        double sm = sqrt(ts_x*ts_x + ts_y*ts_y + ts_z*ts_z);
        double em = sqrt(te_x*te_x + te_y*te_y + te_z*te_z);
        if (sm > 1e-15) {
            prim.tan_start_x = ts_x/sm;
            prim.tan_start_y = ts_y/sm;
            prim.tan_start_z = ts_z/sm;
        }
        if (em > 1e-15) {
            prim.tan_end_x = te_x/em;
            prim.tan_end_y = te_y/em;
            prim.tan_end_z = te_z/em;
        }
    }

    consume_run_marker(ctx, &prim);
    ctx->cb(&prim, ctx->cb_user);
    cache_exit_tangent(ctx, &prim);
    ctx->stats.emitted_arc++;
    ctx->stats.absorbed_total += prim.n_absorbed;
    if (fit->max_deviation > ctx->stats.max_arc_deviation)
        ctx->stats.max_arc_deviation = fit->max_deviation;
}

/* Emit a HELIX primitive (planar arc + axial pitch).  Mirrors emit_arc
 * but writes pitch and uses the helical tangent (radial × axis + pitch
 * along axis). */
static void emit_helix(liscio_ctx_t *ctx, int i0, int i1,
                       const liscio_helix9_fit_t *fit)
{
    if (!ctx->cb) return;
    if (!fit || i0 < 0 || i1 < i0 || i1 >= ctx->n_pts) return;

    liscio_primitive_t prim;
    memset(&prim, 0, sizeof(prim));
    prim.type           = LISCIO_PRIM_HELIX;
    prim.start          = ctx->pts[i0];
    prim.end            = ctx->pts[i1];
    /* Shift center along axis so (start - c) is purely radial.  This
     * lets eval code (and verify) reconstruct as
     *   P(t) = c + r·radial(t·arc_angle) + t·arc_angle·pitch·n
     * with the radial frame anchored at start (matches ARC convention
     * where start lies on the rim at θ=0). */
    {
        double rs_x = ctx->pts[i0].x - fit->cx;
        double rs_y = ctx->pts[i0].y - fit->cy;
        double rs_z = ctx->pts[i0].z - fit->cz;
        double axial = rs_x*fit->nx + rs_y*fit->ny + rs_z*fit->nz;
        prim.cx = fit->cx + axial * fit->nx;
        prim.cy = fit->cy + axial * fit->ny;
        prim.cz = fit->cz + axial * fit->nz;
    }
    prim.nx             = fit->nx;
    prim.ny             = fit->ny;
    prim.nz             = fit->nz;
    prim.radius         = fit->radius;
    prim.arc_angle      = fit->arc_angle;
    prim.pitch          = fit->pitch;
    prim.feedrate       = ctx->feed[i0];
    prim.tag_line_first = ctx->tag[i0];
    prim.tag_line_last  = ctx->tag[i1];
    prim.tag_line_id    = prim.tag_line_first;
    prim.n_absorbed     = i1 - i0;
    fill_common_fields(ctx, &prim);

    /* Tangent = sign·(n × radial) + pitch·n; normalized.
     * For a helix P(θ) = c + r·u(θ) + (z0+p·θ)·n, dP/dθ has magnitude
     * sqrt(r² + p²) and direction = radial-perpendicular + p·axis. */
    {
        double rs_x = prim.start.x - fit->cx;
        double rs_y = prim.start.y - fit->cy;
        double rs_z = prim.start.z - fit->cz;
        double re_x = prim.end.x - fit->cx;
        double re_y = prim.end.y - fit->cy;
        double re_z = prim.end.z - fit->cz;
        double sign = (fit->arc_angle >= 0.0) ? 1.0 : -1.0;
        double ts_x = sign * (fit->ny*rs_z - fit->nz*rs_y) + fit->pitch * fit->nx;
        double ts_y = sign * (fit->nz*rs_x - fit->nx*rs_z) + fit->pitch * fit->ny;
        double ts_z = sign * (fit->nx*rs_y - fit->ny*rs_x) + fit->pitch * fit->nz;
        double te_x = sign * (fit->ny*re_z - fit->nz*re_y) + fit->pitch * fit->nx;
        double te_y = sign * (fit->nz*re_x - fit->nx*re_z) + fit->pitch * fit->ny;
        double te_z = sign * (fit->nx*re_y - fit->ny*re_x) + fit->pitch * fit->nz;
        double sm = sqrt(ts_x*ts_x + ts_y*ts_y + ts_z*ts_z);
        double em = sqrt(te_x*te_x + te_y*te_y + te_z*te_z);
        if (sm > 1e-15) {
            prim.tan_start_x = ts_x/sm;
            prim.tan_start_y = ts_y/sm;
            prim.tan_start_z = ts_z/sm;
        }
        if (em > 1e-15) {
            prim.tan_end_x = te_x/em;
            prim.tan_end_y = te_y/em;
            prim.tan_end_z = te_z/em;
        }
    }

    consume_run_marker(ctx, &prim);
    ctx->cb(&prim, ctx->cb_user);
    cache_exit_tangent(ctx, &prim);
    ctx->stats.emitted_arc++;          /* HELIX counts as arc-family    */
    ctx->stats.absorbed_total += prim.n_absorbed;
    if (fit->max_deviation > ctx->stats.max_arc_deviation)
        ctx->stats.max_arc_deviation = fit->max_deviation;
}

/* ------------------------------------------------------------------ */
/* Try to emit window[0..n-1] as the largest possible primitive.
 * Strategy:
 *   1. If n >= min_arc_pts, try 9D arc fit on [0..n-1].
 *   2. Fall back: [0..n-1] as collinear line (passed tol checks on add).
 * After emission, window collapses to the last point (to seed next run).
 */
static void flush_window_as_primitive(liscio_ctx_t *ctx)
{
    if (ctx->n_pts < 2) {
        /* Nothing to emit. */
        ctx->n_pts = 0;
        return;
    }

    int emitted = 0;

    /* G1-CHAIN (opt-in): when a previous primitive's exit tangent is
     * known, try a G1-constrained cubic Bezier across the WHOLE window
     * first, so this primitive leaves the previous one tangentially
     * (sub-degree boundary → no jerk-kink corner cap).  Gated on tol_xyz
     * via fit_g1; a window that bends too hard for a single G1 cubic
     * falls through to the normal ladder below (which still benefits from
     * the preserved tangent in its own bezier path).  This is the core of
     * the #2 dense-path smoothing fix. */
    if (ctx->cfg.enable_g1_chain && ctx->have_prev_emit_tan
        && ctx->n_pts >= 4 && ctx->cfg.enable_bezier)
    {
        int n = ctx->n_pts;
        double dx = ctx->pts[n-2].x - ctx->pts[n-1].x;
        double dy = ctx->pts[n-2].y - ctx->pts[n-1].y;
        double dz = ctx->pts[n-2].z - ctx->pts[n-1].z;
        double m = sqrt(dx*dx + dy*dy + dz*dz);
        if (m > 1e-15) {
            liscio_bezier9_fit_t bf;
            int rc = liscio_bezier9_fit_g1(ctx, 0, ctx->n_pts - 1,
                ctx->prev_emit_tan_x, ctx->prev_emit_tan_y, ctx->prev_emit_tan_z,
                dx/m, dy/m, dz/m, &bf);
            if (rc == 0) {
                emit_bezier(ctx, 0, ctx->n_pts - 1, &bf);
                emitted = 1;
                liscio_dbg.bez_g1_used++;
            }
        }
    }

    /* Strategy: try most-specific first (preserves CAM intent).
     *   0. LINE if all points colinear (degenerate cases)
     *   1. ARC (3-point circle; exact match for CAM-output G01 approx of arcs)
     *   2. BEZIER cubic (smooth curves, CAM tool-path splines)
     *   3. LINE fallback
     */
    if (!emitted) {
        /* Linearity check: max perpendicular deviation from line pts[0]→pts[n-1]. */
        const liscio_pose_t *p0 = &ctx->pts[0];
        const liscio_pose_t *pN = &ctx->pts[ctx->n_pts - 1];
        double dx = pN->x - p0->x, dy = pN->y - p0->y, dz = pN->z - p0->z;
        double l2 = dx*dx + dy*dy + dz*dz;
        if (l2 > 1e-18) {
            double max_perp = 0.0;
            for (int i = 1; i < ctx->n_pts - 1; i++) {
                double vx = ctx->pts[i].x - p0->x;
                double vy = ctx->pts[i].y - p0->y;
                double vz = ctx->pts[i].z - p0->z;
                double t = (vx*dx + vy*dy + vz*dz) / l2;
                double ex = vx - t*dx, ey = vy - t*dy, ez = vz - t*dz;
                double d = sqrt(ex*ex + ey*ey + ez*ez);
                if (d > max_perp) max_perp = d;
            }
            if (max_perp < ctx->cfg.collinear_tol) {
                emit_line(ctx, 0, ctx->n_pts - 1);
                emitted = 1;
                liscio_dbg.emit_line_collin++;
            }
        }
    }

    if (!emitted && ctx->n_pts >= ctx->cfg.min_arc_pts) {
        liscio_arc9d_fit_t fit;
        int rc = liscio_arc9d_fit(ctx, 0, ctx->n_pts - 1, &fit);
        if (rc == 0) {
            emit_arc(ctx, 0, ctx->n_pts - 1, &fit);
            emitted = 1;
        }
    }
    /* HELIX: planar-arc fit failed (typically because waypoints climb
     * along an axis).  Try fitting a helix before falling back to
     * Bezier; CAM helical-cut and lathe-thread runs benefit. */
    if (!emitted && ctx->n_pts >= ctx->cfg.min_arc_pts) {
        liscio_helix9_fit_t hf;
        if (liscio_helix9_fit(ctx, 0, ctx->n_pts - 1, &hf) == 0) {
            emit_helix(ctx, 0, ctx->n_pts - 1, &hf);
            emitted = 1;
        }
    }
    /* B-spline FIRST (opt-in via prefer_bspline): a C2 B-spline spreads a
     * dense window's curvature so the RT centripetal-jerk cap does not
     * throttle it to a crawl (the #2 dense-path fix).  Tried before Bezier;
     * Bezier still catches windows the B-spline fitter rejects. */
    if (!emitted && ctx->cfg.prefer_bspline && ctx->cfg.enable_bspline
        && ctx->n_pts >= (ctx->cfg.bspline_init_ctrl > 4
                           ? ctx->cfg.bspline_init_ctrl : 4))
    {
        int M_init = ctx->cfg.bspline_init_ctrl;
        if (M_init < 4) M_init = 4;
        if (M_init > LISCIO_BSPLINE_MAX_CTRL) M_init = LISCIO_BSPLINE_MAX_CTRL;
        int M_cap = M_init + 4;
        if (M_cap > LISCIO_BSPLINE_MAX_CTRL) M_cap = LISCIO_BSPLINE_MAX_CTRL;
        for (int M = M_init; M <= M_cap && M <= ctx->n_pts; M += 4) {
            liscio_bspline9_fit_t bs;
            if (liscio_bspline9_fit(ctx, 0, ctx->n_pts - 1, M, &bs) == 0) {
                emit_bspline(ctx, 0, ctx->n_pts - 1, &bs);
                emitted = 1;
                break;
            }
        }
    }

    if (!emitted && ctx->n_pts >= 4 && ctx->cfg.enable_bezier) {
        /* Try simple single-fit first (fast path).  When a previous
         * primitive's exit tangent is available (cross-window G1), use
         * the constrained variant so the join is smooth. */
#ifdef LISCIO_DEBUG_G1
        extern long liscio_dbg_g1_fast_ok, liscio_dbg_g1_fast_fall;
#endif
        liscio_bezier9_fit_t bf;
        int rc = -1;
        int g1_tried = 0;
        if (ctx->have_prev_emit_tan) {
            int n = ctx->n_pts;
            double dx = ctx->pts[n-2].x - ctx->pts[n-1].x;
            double dy = ctx->pts[n-2].y - ctx->pts[n-1].y;
            double dz = ctx->pts[n-2].z - ctx->pts[n-1].z;
            double m = sqrt(dx*dx + dy*dy + dz*dz);
            if (m > 1e-15) {
                rc = liscio_bezier9_fit_g1(ctx, 0, ctx->n_pts - 1,
                    ctx->prev_emit_tan_x, ctx->prev_emit_tan_y, ctx->prev_emit_tan_z,
                    dx/m, dy/m, dz/m, &bf);
                g1_tried = 1;
            }
        }
#ifdef LISCIO_DEBUG_G1
        if (rc == 0) liscio_dbg_g1_fast_ok++;
        else if (g1_tried) liscio_dbg_g1_fast_fall++;
#else
        (void)g1_tried;
#endif
        if (g1_tried) { if (rc == 0) liscio_dbg.bez_g1_used++; else liscio_dbg.bez_g1_miss++; }
        if (rc != 0) rc = liscio_bezier9_fit(ctx, 0, ctx->n_pts - 1, &bf);
        if (rc == 0) {
            emit_bezier(ctx, 0, ctx->n_pts - 1, &bf);
            emitted = 1;
        }
    }

    /* B-spline: try when Bezier alone fails (one multi-CP primitive
     * vs composite-N Beziers).  Progressive M escalation is costly on
     * wild CAM data — the length-reasonableness check rejects bad fits,
     * so cap escalation to two attempts. */
    if (!emitted && ctx->cfg.enable_bspline
        && ctx->n_pts >= (ctx->cfg.bspline_init_ctrl > 4
                           ? ctx->cfg.bspline_init_ctrl : 4))
    {
        int M_init = ctx->cfg.bspline_init_ctrl;
        if (M_init < 4) M_init = 4;
        if (M_init > LISCIO_BSPLINE_MAX_CTRL) M_init = LISCIO_BSPLINE_MAX_CTRL;
        int M_cap = M_init + 4;
        if (M_cap > LISCIO_BSPLINE_MAX_CTRL) M_cap = LISCIO_BSPLINE_MAX_CTRL;
        for (int M = M_init; M <= M_cap && M <= ctx->n_pts; M += 4) {
            liscio_bspline9_fit_t bs;
            if (liscio_bspline9_fit(ctx, 0, ctx->n_pts - 1, M, &bs) == 0) {
                emit_bspline(ctx, 0, ctx->n_pts - 1, &bs);
                emitted = 1;
                break;
            }
        }
    }

    /* Composite Bezier fallback (after B-spline) — splits window when
     * single-span fit too loose.  Cross-window G1: pass the previous
     * primitive's exit tangent so the first sub-Bezier joins smoothly. */
    if (!emitted && ctx->n_pts >= 4 && ctx->cfg.enable_bezier
        && ctx->cfg.enable_composite_bezier)
    {
        const double *left_tan = NULL;
        double left_tan_buf[3];
        if (ctx->have_prev_emit_tan) {
            left_tan_buf[0] = ctx->prev_emit_tan_x;
            left_tan_buf[1] = ctx->prev_emit_tan_y;
            left_tan_buf[2] = ctx->prev_emit_tan_z;
            left_tan = left_tan_buf;
        }
        int n = liscio_bezier9_composite_fit_g1(ctx, 0, ctx->n_pts - 1,
                left_tan, bezier_composite_emit_adapter, ctx);
        if (n > 0) emitted = 1;
    }

    if (!emitted) {
        /* Fallback: emit each window segment as its own LINE so we never
         * silently shortcut a non-collinear run with a chord that exceeds
         * tol_xyz at interior waypoints. */
        for (int k = 0; k < ctx->n_pts - 1; k++) {
            emit_line(ctx, k, k + 1);
            liscio_dbg.emit_line_fallback++;
        }
    }

    /* TEMP debug: per-flush cumulative counter snapshot (flower bisect). */
    if (getenv("TP2_FLUSHDBG")) {
        FILE *ff = fopen("/tmp/tp2_flushdbg.log", "a");
        if (ff) {
            fprintf(ff, "flush npts=%d emitted=%d | corner=%ld disc=%ld maxwin=%ld feedchg=%ld final=%ld | coll=%ld fall=%ld g1u=%ld g1m=%ld absorb=%ld\n",
                ctx->n_pts, emitted,
                liscio_dbg.flush_corner, liscio_dbg.flush_disc,
                liscio_dbg.flush_maxwin, liscio_dbg.flush_feedchg,
                liscio_dbg.flush_final,
                liscio_dbg.emit_line_collin, liscio_dbg.emit_line_fallback,
                liscio_dbg.bez_g1_used, liscio_dbg.bez_g1_miss,
                liscio_dbg.stream_absorb);
            fclose(ff);
        }
    }

    /* Keep last waypoint as seed for next run. */
    liscio_pose_t last = ctx->pts[ctx->n_pts - 1];
    double last_f = ctx->feed[ctx->n_pts - 1];
    int last_tag = ctx->tag[ctx->n_pts - 1];
    ctx->pts[0] = last;
    ctx->feed[0] = last_f;
    ctx->tag[0] = last_tag;
    ctx->n_pts = 1;
}

/* ------------------------------------------------------------------ */
/* Categorise the angular jump from segment pts[n-2]→pts[n-1] to
 * pts[n-1]→end:
 *   0 = below corner_angle_deg (smooth)
 *   1 = soft corner (corner_angle_deg ≤ angle ≤ corner_hard_angle_deg)
 *   2 = hard corner (above corner_hard_angle_deg)
 * Soft + hard distinction matters only in LOOKAHEAD mode; in IMMEDIATE
 * mode the caller treats 1 and 2 identically. */
static int corner_level(const liscio_ctx_t *ctx,
                        const liscio_pose_t *end)
{
    if (ctx->n_pts < 2 || ctx->cfg.corner_angle_deg <= 0.0) return 0;
    const liscio_pose_t *pA = &ctx->pts[ctx->n_pts - 2];
    const liscio_pose_t *pB = &ctx->pts[ctx->n_pts - 1];

    double ax = pB->x - pA->x, ay = pB->y - pA->y, az = pB->z - pA->z;
    double bx = end->x - pB->x, by = end->y - pB->y, bz = end->z - pB->z;
    double la = sqrt(ax*ax + ay*ay + az*az);
    double lb = sqrt(bx*bx + by*by + bz*bz);
    if (la < 1e-12 || lb < 1e-12) return 0;

    double cosang = (ax*bx + ay*by + az*bz) / (la * lb);
    if (cosang > 1.0) cosang = 1.0;
    if (cosang < -1.0) cosang = -1.0;
    double cos_soft = cos(ctx->cfg.corner_angle_deg * M_PI / 180.0);
    if (cosang >= cos_soft) return 0;            /* below soft threshold */
    double hard = ctx->cfg.corner_hard_angle_deg;
    if (hard <= 0) hard = ctx->cfg.corner_angle_deg; /* fallback */
    double cos_hard = cos(hard * M_PI / 180.0);
    return (cosang < cos_hard) ? 2 : 1;
}

static void flush_pending_arc(liscio_ctx_t *ctx);
static void flush_arc_helix_buffer(liscio_ctx_t *ctx);
static void arcg2_flush_pending(liscio_ctx_t *ctx);

int liscio_add_line(liscio_ctx_t *ctx,
                     const liscio_pose_t *start,
                     const liscio_pose_t *end,
                     double feedrate,
                     int tag_line_id)
{
    if (!ctx || !start || !end) return -1;
    ctx->stats.input_count++;

    /* G1 interrupts any pending arc-merge / arc-helix / fusion run. */
    if (ctx->has_pending_arc) flush_pending_arc(ctx);
    if (ctx->arc_buf_n > 0)   flush_arc_helix_buffer(ctx);
    if (ctx->has_arcg2_pend)  arcg2_flush_pending(ctx);

    /* Micro-cusp cluster collapse: while inside a hard-corner thrash cluster,
     * decimate points that stay within collapse_cusp_tol of the anchor — they
     * are geometrically indistinct at that tolerance and each would otherwise
     * force its own v=0 junction (N stops in a few µm → millions of jerk).
     * Resume normal fitting when the path escapes the ball. */
    if (ctx->in_cluster) {
        double _ct = (ctx->cfg.collapse_cusp_tol > 0.0)
                     ? ctx->cfg.collapse_cusp_tol : ctx->cfg.tol_xyz;
        double dx = end->x - ctx->cluster_anchor.x;
        double dy = end->y - ctx->cluster_anchor.y;
        double dz = end->z - ctx->cluster_anchor.z;
        if (dx*dx + dy*dy + dz*dz <= _ct * _ct)
            return 0;                       /* decimate: drop this point */
        /* Escaped the ball → cluster ends.  Append `end` to the anchored
         * window ([anchor]) as the single exit segment: the anchor becomes
         * ONE junction (v=0 when the net turn across the cluster is sharp)
         * for the whole collapsed thrash region.  Bypass the continuity
         * check — `start` no longer matches the anchor (intermediate points
         * were decimated). */
        ctx->in_cluster = 0;
        if (ctx->n_pts < LISCIO_MAX_WINDOW) {
            ctx->pts[ctx->n_pts]  = *end;
            ctx->feed[ctx->n_pts] = feedrate;
            ctx->tag[ctx->n_pts]  = tag_line_id;
            ctx->n_pts++;
        }
        return 0;
    }

    /* Seed window on first call or after reset/flush. */
    if (ctx->n_pts == 0) {
        ctx->pts[0]  = *start;
        ctx->feed[0] = feedrate;
        ctx->tag[0]  = tag_line_id;
        ctx->n_pts   = 1;
    } else {
        /* Verify continuity: start must match previous end. Drop +
         * re-seed if not (upstream should handle but be defensive). */
        const liscio_pose_t *prev_end = &ctx->pts[ctx->n_pts - 1];
        double dx = start->x - prev_end->x;
        double dy = start->y - prev_end->y;
        double dz = start->z - prev_end->z;
        double d2 = dx*dx + dy*dy + dz*dz;
        if (d2 > (ctx->cfg.tol_xyz * ctx->cfg.tol_xyz) * 4.0) {
            /* Discontinuity — flush, then seed.  Drop cross-flush
             * tangent and mark the next emission as a new run: the
             * caller probably skipped a G0 rapid without calling
             * liscio_mark_stop. */
            liscio_dbg.flush_disc++;
            flush_window_as_primitive(ctx);
            ctx->have_prev_emit_tan  = 0;
            ctx->pts[0]  = *start;
            ctx->feed[0] = feedrate;
            ctx->tag[0]  = tag_line_id;
            ctx->n_pts   = 1;
        }
    }

    /* Corner detection.
     *
     * In IMMEDIATE mode any soft-or-hard corner triggers an instant flush.
     * In LOOKAHEAD mode only HARD corners flush; soft corners are deferred
     * — composite Bezier fit downstream often absorbs them with G1
     * tangent continuity, dropping tan_mean substantially. */
    {
        int level = corner_level(ctx, end);
        int do_flush = 0;
        if (level >= 2) {
            do_flush = 1;                                /* hard: always */
        } else if (level == 1) {
            do_flush = (ctx->cfg.corner_detection_mode
                       != LISCIO_CORNER_LOOKAHEAD);     /* soft in legacy mode */
        }
        /* Micro-cusp cluster entry: a HARD corner between two sub-tolerance
         * segments is direction thrashing at the resolution floor.  Flush the
         * clean pre-cluster geometry once, then enter collapse mode so the
         * whole thrash region folds into a single stop junction instead of
         * one v=0 junction per corner. */
        int enter_cluster = 0;
        if (do_flush && level >= 2 && ctx->cfg.collapse_cusp_enable
            && ctx->n_pts >= 2) {
            double _ct = (ctx->cfg.collapse_cusp_tol > 0.0)
                         ? ctx->cfg.collapse_cusp_tol : ctx->cfg.tol_xyz;
            int n = ctx->n_pts;
            double ax = ctx->pts[n-1].x - ctx->pts[n-2].x;
            double ay = ctx->pts[n-1].y - ctx->pts[n-2].y;
            double az = ctx->pts[n-1].z - ctx->pts[n-2].z;
            double bx = end->x - ctx->pts[n-1].x;
            double by = end->y - ctx->pts[n-1].y;
            double bz = end->z - ctx->pts[n-1].z;
            if (ax*ax + ay*ay + az*az < _ct * _ct &&
                bx*bx + by*by + bz*bz < _ct * _ct)
                enter_cluster = 1;
        }
        if (do_flush) {
            liscio_dbg.flush_corner++;
            flush_window_as_primitive(ctx);
            ctx->have_prev_emit_tan = 0;   /* corner is G1 discontinuity */
            /* After flush, n_pts = 1 with last waypoint = corner pt. */
        }
        if (enter_cluster) {
            ctx->in_cluster     = 1;
            ctx->cluster_anchor = ctx->pts[0];   /* = the corner point */
            return 0;                            /* decimate current `end` */
        }
    }

    /* Same-F merge gate (cfg.merge_same_feed_only, default on).  A fitted
     * primitive carries ONE feedrate (prim.feedrate = feed of the window's
     * first segment), so a window may only span segments that share F: an
     * F change flushes the accumulated (all-old-F) window and the incoming
     * segment starts a fresh window under the new F.  Runs AFTER corner
     * detection so a corner+F-change junction keeps its corner handling
     * (flush + tangent clear); a pure F change is not a geometric break,
     * so the cross-flush tangent survives (same as a maxwin flush). */
    if (ctx->cfg.merge_same_feed_only && ctx->n_pts >= 1
        && !feed_match(ctx->feed[ctx->n_pts - 1], feedrate)) {
        if (getenv("TP2_FLUSHDBG")) {   /* TEMP: flower bisect */
            FILE *ff = fopen("/tmp/tp2_feedchg.log", "a");
            if (ff) { fprintf(ff, "feedchg stored=%.9f incoming=%.9f npts=%d\n",
                              ctx->feed[ctx->n_pts - 1], feedrate, ctx->n_pts); fclose(ff); }
        }
        if (ctx->n_pts >= 2) {
            liscio_dbg.flush_feedchg++;
            flush_window_as_primitive(ctx);   /* keeps last pt as seed */
        }
        /* Re-label the seed waypoint: it only anchors the incoming segment,
         * whose feedrate rules the window now starting (the seed's stored
         * feed is the PREVIOUS run's F — emitting it would mis-feed the
         * first segment after every F-change or flush boundary). */
        ctx->feed[0] = feedrate;
    }

    /* Streaming collinear absorption (header: collinear_tol).  Extend the
     * last chord pts[n-2]→pts[n-1] to pts[n-2]→end when the abandoned
     * waypoint stays inside tol of the new chord.  Cuts µm G1 junction
     * count (130207LZW engraving) without touching arc_g2fuse (add_arc).
     * Gates: same sense (no reverse), Z/W peck not reversed, feed match. */
    if (ctx->n_pts >= 2) {
        const liscio_pose_t *p0 = &ctx->pts[ctx->n_pts - 2];
        const liscio_pose_t *p1 = &ctx->pts[ctx->n_pts - 1];
        double ax = p1->x - p0->x, ay = p1->y - p0->y, az = p1->z - p0->z;
        double aw = p1->w - p0->w;
        double bx = end->x - p1->x, by = end->y - p1->y, bz = end->z - p1->z;
        double bw = end->w - p1->w;
        double la = sqrt(ax * ax + ay * ay + az * az);
        double lb = sqrt(bx * bx + by * by + bz * bz);
        double ctol = ctx->cfg.collinear_tol;
        if (ctol <= 0.0) ctol = ctx->cfg.tol_xyz;
        /* 6.6.135: a 0.01 cap left flower Y=121 12 µm XZ bows as
         * LINEAR wraps (rec-129-fls +81k) even though G64 P=0.05 and
         * g1_chain already trust that band.  LZW P=0.005 stays 5 µm.
         * Ceiling 0.05 so looser-than-flower G64 cannot eat visible
         * bends.  same_dir still ~8°.  Do not reopen 128–134. */
        if (ctol > 0.05) ctol = 0.05;
        if (la > 1e-15 && lb > 1e-15 && ctol > 0.0) {
            double dot = (ax * bx + ay * by + az * bz) / (la * lb);
            /* ~8° — below IMMEDIATE soft-corner flush (15°); path still
             * gated by collinear_tol / ctol so G64 P is respected. */
            int same_dir = (dot > 0.99);
            const double peck_eps = 1e-6;
            int z_rev = (fabs(az) > peck_eps && fabs(bz) > peck_eps && az * bz < 0.0);
            int w_rev = (fabs(aw) > peck_eps && fabs(bw) > peck_eps && aw * bw < 0.0);
            double fprev = ctx->feed[ctx->n_pts - 1];
            int feed_ok = (fprev <= 0.0 || feedrate <= 0.0
                           || fabs(fprev - feedrate) < 1e-9 * fmax(1.0, fabs(fprev)));
            if (same_dir && !z_rev && !w_rev && feed_ok) {
                double cx = end->x - p0->x, cy = end->y - p0->y, cz = end->z - p0->z;
                double lc2 = cx * cx + cy * cy + cz * cz;
                double vx = p1->x - p0->x, vy = p1->y - p0->y, vz = p1->z - p0->z;
                double t = (lc2 > 1e-30) ? (vx * cx + vy * cy + vz * cz) / lc2 : 0.0;
                if (t < 0.0) t = 0.0;
                else if (t > 1.0) t = 1.0;
                double qx = p0->x + t * cx - p1->x;
                double qy = p0->y + t * cy - p1->y;
                double qz = p0->z + t * cz - p1->z;
                double dev = sqrt(qx * qx + qy * qy + qz * qz);
                if (dev <= ctol) {
                    ctx->pts[ctx->n_pts - 1]  = *end;
                    ctx->feed[ctx->n_pts - 1] = feedrate;
                    ctx->tag[ctx->n_pts - 1]  = tag_line_id;
                    liscio_dbg.stream_absorb++;
                    return 0;
                }
            }
        }
    }

    /* Append new endpoint. */
    if (ctx->n_pts >= ctx->cfg.max_window) {
        /* Window full — flush and re-seed. */
        liscio_dbg.flush_maxwin++;
        flush_window_as_primitive(ctx);
    }

    ctx->pts[ctx->n_pts]  = *end;
    ctx->feed[ctx->n_pts] = feedrate;
    ctx->tag[ctx->n_pts]  = tag_line_id;
    ctx->n_pts++;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Compute arc tangents + invoke ctx->cb.  Shared by direct-passthrough
 * and pending-arc-flush paths. */
static void emit_arc_primitive(liscio_ctx_t *ctx,
                               const liscio_primitive_t *arc,
                               double feedrate, int tag_line_first,
                               int tag_line_last,
                               int n_absorbed)
{
    if (!ctx->cb) return;
    liscio_primitive_t pass = *arc;
    if (feedrate > 0.0) pass.feedrate = feedrate;
    if (tag_line_first > 0) {
        pass.tag_line_first = tag_line_first;
        pass.tag_line_id    = tag_line_first;
    }
    if (tag_line_last > 0) pass.tag_line_last = tag_line_last;
    else                   pass.tag_line_last = pass.tag_line_first;
    pass.n_absorbed = n_absorbed;

    /* Tangent at start = (normal × start_radial), end = (normal × end_radial),
     * sign-flipped for CW (negative arc_angle).  For HELIX (multi-loop)
     * the end tangent is identical in direction to the single-loop case
     * since arc_angle mod 2π determines orientation. */
    double rsx = arc->start.x - arc->cx;
    double rsy = arc->start.y - arc->cy;
    double rsz = arc->start.z - arc->cz;
    double rex = arc->end.x - arc->cx;
    double rey = arc->end.y - arc->cy;
    double rez = arc->end.z - arc->cz;
    double tsx = arc->ny*rsz - arc->nz*rsy;
    double tsy = arc->nz*rsx - arc->nx*rsz;
    double tsz = arc->nx*rsy - arc->ny*rsx;
    double tex = arc->ny*rez - arc->nz*rey;
    double tey = arc->nz*rex - arc->nx*rez;
    double tez = arc->nx*rey - arc->ny*rex;
    if (arc->arc_angle < 0) {
        tsx = -tsx; tsy = -tsy; tsz = -tsz;
        tex = -tex; tey = -tey; tez = -tez;
    }
    double tsm = sqrt(tsx*tsx + tsy*tsy + tsz*tsz);
    double tem = sqrt(tex*tex + tey*tey + tez*tez);
    if (tsm > 1e-15) { pass.tan_start_x=tsx/tsm; pass.tan_start_y=tsy/tsm; pass.tan_start_z=tsz/tsm; }
    if (tem > 1e-15) { pass.tan_end_x=tex/tem; pass.tan_end_y=tey/tem; pass.tan_end_z=tez/tem; }

    /* Promote to HELIX type when configured + |arc_angle| > 2π. */
    if (ctx->cfg.arc_merge_mode == 2 /*HELIX*/ &&
        fabs(arc->arc_angle) > 2.0 * M_PI + 1e-9)
        pass.type = LISCIO_PRIM_HELIX;

    consume_run_marker(ctx, &pass);
    ctx->cb(&pass, ctx->cb_user);
    cache_exit_tangent(ctx, &pass);
    ctx->stats.emitted_arc++;
    ctx->stats.absorbed_total += pass.n_absorbed;
}

/* Flush any buffered pending arc (arc-merge).  Called on:
 *   - add_arc with non-mergeable next arc
 *   - add_line (G1 interrupts arc run)
 *   - liscio_flush
 * Resets has_pending_arc afterward. */
static void flush_pending_arc(liscio_ctx_t *ctx)
{
    if (!ctx->has_pending_arc) return;
    emit_arc_primitive(ctx, &ctx->pending_arc,
                       ctx->pending_feed,
                       ctx->pending_tag,
                       ctx->pending_tag_last,
                       ctx->pending_n_absorbed);
    ctx->has_pending_arc = 0;
}

/* ------------------------------------------------------------------ */
/* Arc-helix LSQ detection: takes N buffered ARCs and tests whether
 * they jointly fit a helical primitive (analytic, no sampling).
 * On success emits one HELIX with merged params; on failure emits
 * each buffered arc individually.  Always clears the buffer. */
static void flush_arc_helix_buffer(liscio_ctx_t *ctx)
{
    if (ctx->arc_buf_n == 0) return;
    int N = ctx->arc_buf_n;

    if (N < 3) {
        /* Too few to call a helix run; emit individually. */
        for (int i = 0; i < N; i++) {
            emit_arc_primitive(ctx, &ctx->arc_buf[i], ctx->arc_buf_feed,
                               ctx->arc_buf[i].tag_line_first,
                               ctx->arc_buf[i].tag_line_last, 1);
        }
        ctx->arc_buf_n = 0;
        return;
    }

    /* 1. Sign reference (must agree across the run). */
    double sign_ref = (ctx->arc_buf[0].arc_angle >= 0) ? 1.0 : -1.0;
    for (int i = 1; i < N; i++) {
        double s = (ctx->arc_buf[i].arc_angle >= 0) ? 1.0 : -1.0;
        if (s != sign_ref) goto fail;
    }

    /* 2. Mean axis (assume signs already aligned by buffer admission). */
    double nx = 0, ny = 0, nz = 0;
    for (int i = 0; i < N; i++) {
        nx += ctx->arc_buf[i].nx;
        ny += ctx->arc_buf[i].ny;
        nz += ctx->arc_buf[i].nz;
    }
    double nm = sqrt(nx*nx + ny*ny + nz*nz);
    if (nm < 1e-12) goto fail;
    nx /= nm; ny /= nm; nz /= nm;

    /* 3. Mean radius. */
    double r_avg = 0;
    for (int i = 0; i < N; i++) r_avg += ctx->arc_buf[i].radius;
    r_avg /= N;

    /* 4. Mean center (anchor for axis line). */
    double mcx = 0, mcy = 0, mcz = 0;
    for (int i = 0; i < N; i++) {
        mcx += ctx->arc_buf[i].cx;
        mcy += ctx->arc_buf[i].cy;
        mcz += ctx->arc_buf[i].cz;
    }
    mcx /= N; mcy /= N; mcz /= N;

    /* 5. Project each center onto axis line through (mcx,mcy,mcz);
     *    record perpendicular offset (residual 1) and axial position. */
    double axial[LISCIO_ARC_HELIX_BUF];
    double max_off = 0, max_rerr = 0;
    for (int i = 0; i < N; i++) {
        double dx = ctx->arc_buf[i].cx - mcx;
        double dy = ctx->arc_buf[i].cy - mcy;
        double dz = ctx->arc_buf[i].cz - mcz;
        axial[i] = dx*nx + dy*ny + dz*nz;
        double px = dx - axial[i]*nx;
        double py = dy - axial[i]*ny;
        double pz = dz - axial[i]*nz;
        double off = sqrt(px*px + py*py + pz*pz);
        if (off > max_off) max_off = off;
        double rerr = fabs(ctx->arc_buf[i].radius - r_avg);
        if (rerr > max_rerr) max_rerr = rerr;
    }

    /* 6. LSQ fit pitch: axial[i] vs cumulative angle theta_cum[i] (at
     *    the START of each arc, with theta_cum[0] = 0). */
    double theta_cum[LISCIO_ARC_HELIX_BUF];
    theta_cum[0] = 0;
    for (int i = 1; i < N; i++)
        theta_cum[i] = theta_cum[i-1] + ctx->arc_buf[i-1].arc_angle;

    double tmean = 0, zmean = 0;
    for (int i = 0; i < N; i++) { tmean += theta_cum[i]; zmean += axial[i]; }
    tmean /= N; zmean /= N;
    double Stt = 0, Stz = 0;
    for (int i = 0; i < N; i++) {
        double dt = theta_cum[i] - tmean;
        Stt += dt * dt;
        Stz += dt * (axial[i] - zmean);
    }
    if (Stt < 1e-15) goto fail;
    double pitch = Stz / Stt;
    double z0    = zmean - pitch * tmean;

    /* 7. Pitch residual. */
    double max_pres = 0;
    for (int i = 0; i < N; i++) {
        double pred = z0 + pitch * theta_cum[i];
        double res = fabs(axial[i] - pred);
        if (res > max_pres) max_pres = res;
    }

    /* 8. Combined max-deviation gate. */
    double max_dev = max_off;
    if (max_rerr > max_dev) max_dev = max_rerr;
    if (max_pres > max_dev) max_dev = max_pres;
    if (max_dev > ctx->cfg.tol_xyz) goto fail;

    /* 9. Build merged HELIX primitive.  Center on axis at first arc's
     *    start plane: foot of arc[0].cx on the axis line. */
    double a0_dx = ctx->arc_buf[0].cx - mcx;
    double a0_dy = ctx->arc_buf[0].cy - mcy;
    double a0_dz = ctx->arc_buf[0].cz - mcz;
    double a0_axial = a0_dx*nx + a0_dy*ny + a0_dz*nz;
    double cx0 = mcx + a0_axial * nx;
    double cy0 = mcy + a0_axial * ny;
    double cz0 = mcz + a0_axial * nz;

    /* Cumulative arc_angle across the run. */
    double total_angle = 0;
    for (int i = 0; i < N; i++) total_angle += ctx->arc_buf[i].arc_angle;

    liscio_primitive_t merged = ctx->arc_buf[0];
    merged.cx        = cx0;
    merged.cy        = cy0;
    merged.cz        = cz0;
    merged.nx        = nx;
    merged.ny        = ny;
    merged.nz        = nz;
    merged.radius    = r_avg;
    merged.arc_angle = total_angle;
    merged.pitch     = pitch;
    merged.start     = ctx->arc_buf[0].start;
    merged.end       = ctx->arc_buf[N-1].end;
    merged.type      = LISCIO_PRIM_HELIX;

    if (ctx->cb) {
        liscio_primitive_t pass = merged;
        pass.feedrate       = ctx->arc_buf_feed;
        pass.tag_line_first = ctx->arc_buf_tag_first;
        pass.tag_line_last  = ctx->arc_buf_tag_last;
        pass.tag_line_id    = ctx->arc_buf_tag_first;
        pass.n_absorbed     = N;

        /* Tangent: sign·(n × radial_start) + pitch·n, normalized. */
        double rsx = pass.start.x - pass.cx;
        double rsy = pass.start.y - pass.cy;
        double rsz = pass.start.z - pass.cz;
        double rex = pass.end.x - pass.cx;
        double rey = pass.end.y - pass.cy;
        double rez = pass.end.z - pass.cz;
        double sign = sign_ref;
        double tsx = sign*(ny*rsz - nz*rsy) + pitch*nx;
        double tsy = sign*(nz*rsx - nx*rsz) + pitch*ny;
        double tsz = sign*(nx*rsy - ny*rsx) + pitch*nz;
        double tex = sign*(ny*rez - nz*rey) + pitch*nx;
        double tey = sign*(nz*rex - nx*rez) + pitch*ny;
        double tez = sign*(nx*rey - ny*rex) + pitch*nz;
        double sm = sqrt(tsx*tsx + tsy*tsy + tsz*tsz);
        double em = sqrt(tex*tex + tey*tey + tez*tez);
        if (sm > 1e-15) { pass.tan_start_x = tsx/sm; pass.tan_start_y = tsy/sm; pass.tan_start_z = tsz/sm; }
        if (em > 1e-15) { pass.tan_end_x = tex/em; pass.tan_end_y = tey/em; pass.tan_end_z = tez/em; }

        consume_run_marker(ctx, &pass);
        ctx->cb(&pass, ctx->cb_user);
        cache_exit_tangent(ctx, &pass);
        ctx->stats.emitted_arc++;
        ctx->stats.absorbed_total += N;
        if (max_dev > ctx->stats.max_arc_deviation)
            ctx->stats.max_arc_deviation = max_dev;
    }
    ctx->arc_buf_n = 0;
    return;

fail:
    /* Helix detection failed — emit each buffered arc individually. */
    for (int i = 0; i < N; i++) {
        emit_arc_primitive(ctx, &ctx->arc_buf[i], ctx->arc_buf_feed,
                           ctx->arc_buf[i].tag_line_first,
                           ctx->arc_buf[i].tag_line_last, 1);
    }
    ctx->arc_buf_n = 0;
}

/* Lightweight admission check: is the new arc plausibly part of the
 * helix run buffered so far?  Strict rules (axis ≈ same, radius ≈ same,
 * sign same, start ≈ previous end) act as a coarse pre-filter; the
 * final accept/reject decision is the LSQ residual gate in the flush
 * routine. */
static int arc_helix_admissible(const liscio_primitive_t *prev,
                                 const liscio_primitive_t *next,
                                 double radius_tol, double xy_tol)
{
    /* Same rotation direction. */
    if ((prev->arc_angle >= 0) != (next->arc_angle >= 0)) return 0;
    /* Axes within ~3°. */
    double dot = prev->nx*next->nx + prev->ny*next->ny + prev->nz*next->nz;
    if (dot < 0.998) return 0;
    /* Radius within tolerance. */
    if (fabs(prev->radius - next->radius) > radius_tol) return 0;
    /* Endpoint continuity. */
    double dx = next->start.x - prev->end.x;
    double dy = next->start.y - prev->end.y;
    double dz = next->start.z - prev->end.z;
    if (sqrt(dx*dx + dy*dy + dz*dz) > xy_tol) return 0;
    return 1;
}

/* Compatibility check: can `b` follow `a` as part of the same merged
 * helical/planar arc?  Requires:
 *   - shared center (cx, cy, cz within xy_tol)
 *   - shared normal (parallel within ~1e-9 cross-product magnitude)
 *   - shared radius (within xy_tol)
 *   - same rotation direction (arc_angle sign matches)
 *   - b.start ≈ a.end (within xy_tol)
 *   - helical pitch consistent (within pitch_tol mm/rad)
 */
static int arcs_mergeable(const liscio_primitive_t *a,
                          const liscio_primitive_t *b,
                          double xy_tol, double pitch_tol)
{
    if (xy_tol < 1e-15) xy_tol = 1e-6;
    if (pitch_tol < 1e-15) pitch_tol = 1e-6;
    /* Center match. */
    if (fabs(a->cx - b->cx) > xy_tol) return 0;
    if (fabs(a->cy - b->cy) > xy_tol) return 0;
    if (fabs(a->cz - b->cz) > xy_tol) return 0;
    /* Normal parallel (cross product near 0). */
    double cx = a->ny*b->nz - a->nz*b->ny;
    double cy = a->nz*b->nx - a->nx*b->nz;
    double cz = a->nx*b->ny - a->ny*b->nx;
    if (sqrt(cx*cx + cy*cy + cz*cz) > 1e-9) return 0;
    /* Same direction (positive dot). */
    if (a->nx*b->nx + a->ny*b->ny + a->nz*b->nz < 0.999) return 0;
    /* Radius. */
    if (fabs(a->radius - b->radius) > xy_tol) return 0;
    /* Same rotation direction. */
    if ((a->arc_angle >= 0) != (b->arc_angle >= 0)) return 0;
    /* Continuity: b.start ≈ a.end. */
    double dxs = b->start.x - a->end.x;
    double dys = b->start.y - a->end.y;
    double dzs = b->start.z - a->end.z;
    if (sqrt(dxs*dxs + dys*dys + dzs*dzs) > xy_tol) return 0;
    /* Pitch consistency for helical motion.
     * Pitch = (axial_displacement) / (arc_angle), where axial_displ. =
     * (end - start) · normal. */
    double a_axial = (a->end.x - a->start.x)*a->nx
                   + (a->end.y - a->start.y)*a->ny
                   + (a->end.z - a->start.z)*a->nz;
    double b_axial = (b->end.x - b->start.x)*b->nx
                   + (b->end.y - b->start.y)*b->ny
                   + (b->end.z - b->start.z)*b->nz;
    if (fabs(a->arc_angle) < 1e-15 || fabs(b->arc_angle) < 1e-15) return 0;
    double pitch_a = a_axial / a->arc_angle;
    double pitch_b = b_axial / b->arc_angle;
    if (fabs(pitch_a - pitch_b) > pitch_tol) return 0;
    return 1;
}

/* ==================================================================
 *  Arc-chain G2 spline fusion  (cfg.arc_g2fuse)
 * ==================================================================
 * Emit ONE G2/G3 arc as an ANALYTIC C2 quintic (or, on a wide sweep or a
 * curvature-mismatched seam, a bisected run of quintics), reusing the
 * validated liscio_quintic9_g2_fit direct-construction + tolerance gate.
 *
 * Unlike arc_subseg_samples (which samples arcs into LINE chords pushed
 * through the G1 LSQ pipeline — injecting curvature noise → ripple), this
 * samples the EXACT circle into the fit's check window and relies on the
 * fit pinning the EXACT κ=1/R (est_endpoint_kappa on 3 exact circle points
 * = 1/R; a chained arc's start κ is pinned to the neighbour via the ctx G2
 * chain).  The quintic therefore follows the TRUE arc, and consecutive
 * tangent-continuous arcs share seam curvature (C2) → no curvature step →
 * no velocity hunting.  A real corner or turn-sign flip breaks the chain so
 * the seam falls to a blend join. */

/* Reconstructed arc geometry, mirroring PmCircle / pmCirclePoint so SPIRAL
 * (radius growth) and HELIX (axial advance) are reproduced exactly — thomam's
 * "arcs" are expanding spirals, so a constant-R circle is badly wrong. */
typedef struct {
    double C[3];        /* center                                   */
    double rtan[3];     /* start radial (= start − center), |·| = R0 */
    double rperp[3];    /* n × rtan (same magnitude R0, CCW dir)     */
    double R0;          /* start radius                             */
    double spiral;      /* Rend − R0 (total radial growth)          */
    double hvec[3];     /* total axial (helix) vector = ((E−S)·n)·n */
    double phi;         /* signed total sweep (arc_angle)           */
    /* endpoint ABC/UVW for linear interpolation */
    liscio_pose_t s, e;
} arc_geom_t;

/* Exact pose on the arc at fraction `frac` ∈ [0,1], mirroring pmCirclePoint:
 * base = rtan·cos(a) + rperp·sin(a) (|·|=R0), then add unit(base)·scale·spiral
 * (spiral growth) + hvec·scale (helix), scale = frac.  ABC/UVW interpolate
 * linearly. */
static void arc_pose_at(const arc_geom_t *g, double frac, liscio_pose_t *out)
{
    double a  = frac * g->phi;
    double ca = cos(a), sa = sin(a);
    double base[3];
    base[0] = g->rtan[0]*ca + g->rperp[0]*sa;
    base[1] = g->rtan[1]*ca + g->rperp[1]*sa;
    base[2] = g->rtan[2]*ca + g->rperp[2]*sa;
    double bm = sqrt(base[0]*base[0] + base[1]*base[1] + base[2]*base[2]);
    double grow = (bm > 1e-15) ? (frac * g->spiral / bm) : 0.0;  /* unit(base)·scale·spiral */
    out->x = g->C[0] + base[0]*(1.0 + grow) + frac*g->hvec[0];
    out->y = g->C[1] + base[1]*(1.0 + grow) + frac*g->hvec[1];
    out->z = g->C[2] + base[2]*(1.0 + grow) + frac*g->hvec[2];
    out->a = g->s.a + frac*(g->e.a - g->s.a);
    out->b = g->s.b + frac*(g->e.b - g->s.b);
    out->c = g->s.c + frac*(g->e.c - g->s.c);
    out->u = g->s.u + frac*(g->e.u - g->s.u);
    out->v = g->s.v + frac*(g->e.v - g->s.v);
    out->w = g->s.w + frac*(g->e.w - g->s.w);
}

/* Unit FORWARD tangent at fraction `frac`, by central finite difference of
 * arc_pose_at — robust for spiral + helix without hand-derived derivatives. */
static void arc_tan_at(const arc_geom_t *g, double frac, double t[3])
{
    double h = 1e-5;
    double f0 = frac - h, f1 = frac + h;
    if (f0 < 0.0) { f0 = 0.0; f1 = 2.0*h; }
    if (f1 > 1.0) { f1 = 1.0; f0 = 1.0 - 2.0*h; }
    liscio_pose_t p0, p1;
    arc_pose_at(g, f0, &p0);
    arc_pose_at(g, f1, &p1);
    t[0] = p1.x - p0.x; t[1] = p1.y - p0.y; t[2] = p1.z - p0.z;
    double m = sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);
    if (m > 1e-15) { t[0]/=m; t[1]/=m; t[2]/=m; }
}

#define ARCG2_M         16   /* exact-arc check points per quintic fit */
#define ARCG2_MAXDEPTH   7   /* bisection floor → native sub-arc fallback */

int liscio_arc_g2_burst_bound(void)
{
    /* Full bisection tree of depth ARCG2_MAXDEPTH = 2^MAXDEPTH leaves, each
     * emitting one primitive (quintic fit or native sub-arc floor). */
    return 1 << ARCG2_MAXDEPTH;
}

/* Recursively emit the arc range [f0,f1] as one quintic; bisect on fit miss.
 * t1o/teo (nullable): forward-tangent OVERRIDES applied at the arc's own outer
 * endpoints (f=0 / f=1) — used to absorb a soft seam kink by handing both
 * sides of the seam the AVERAGED tangent.  Interior bisection endpoints always
 * use the exact arc tangent. */
static void arc_g2_emit_range(liscio_ctx_t *ctx, const liscio_primitive_t *arc,
                              const arc_geom_t *g,
                              double feed, int tag, double f0, double f1, int depth,
                              const double *t1o, const double *teo)
{
    for (int i = 0; i < ARCG2_M; i++) {
        double fr = f0 + (f1 - f0) * ((double)i / (ARCG2_M - 1));
        arc_pose_at(g, fr, &ctx->pts[i]);
        ctx->feed[i] = feed;
        ctx->tag[i]  = tag;
    }
    ctx->n_pts = ARCG2_M;

    double t1[3], te[3];
    if (t1o && f0 <= 0.0) { t1[0]=t1o[0]; t1[1]=t1o[1]; t1[2]=t1o[2]; }
    else arc_tan_at(g, f0, t1);
    if (teo && f1 >= 1.0) { te[0]=teo[0]; te[1]=teo[1]; te[2]=teo[2]; }
    else arc_tan_at(g, f1, te);   /* forward end tangent */

    liscio_bspline9_fit_t q;
    /* quintic t2 convention = INWARD at the end (= −forward). */
    if (liscio_quintic9_g2_fit(ctx, 0, ARCG2_M - 1,
                               t1[0], t1[1], t1[2],
                               -te[0], -te[1], -te[2], &q) == 0) {
        emit_bspline(ctx, 0, ARCG2_M - 1, &q);   /* clears G2 chain */
        ctx->prev_emit_kappa   = q.end_kappa;    /* re-arm: next fit pins to this κ */
        ctx->prev_emit_norm_x  = q.end_norm_x;
        ctx->prev_emit_norm_y  = q.end_norm_y;
        ctx->prev_emit_norm_z  = q.end_norm_z;
        ctx->have_prev_emit_kappa = 1;
        ctx->stats.emitted_g2_quintic++;
        return;
    }

    if (depth >= ARCG2_MAXDEPTH) {
        /* Pinned C2 from a struggling neighbor can reject a circular leaf
         * that would otherwise fit.  Native fallback is a few-cycle
         * TC_CIRCULAR whose independent scatti then sign-flips through the
         * wrap (thomam.random X −55k; const same seam already ~48k).  Unpin
         * and retry: the leaf stays a quintic so the G2 blend can join. */
        ctx->have_prev_emit_kappa = 0;
        if (liscio_quintic9_g2_fit(ctx, 0, ARCG2_M - 1,
                                   t1[0], t1[1], t1[2],
                                   -te[0], -te[1], -te[2], &q) == 0) {
            emit_bspline(ctx, 0, ARCG2_M - 1, &q);
            ctx->prev_emit_kappa   = q.end_kappa;
            ctx->prev_emit_norm_x  = q.end_norm_x;
            ctx->prev_emit_norm_y  = q.end_norm_y;
            ctx->prev_emit_norm_z  = q.end_norm_z;
            ctx->have_prev_emit_kappa = 1;
            ctx->stats.emitted_g2_quintic++;
            return;
        }
        ctx->have_prev_emit_kappa = 0;
        /* Floor: a native TC_CIRCULAR of this leftover is a 4-cycle
         * independent scatti whose wrap sign-flips axis jerk (thomam
         * X −55k).  The span was already judged too small/bent for a
         * quintic — emit the chord.  cache_exit_tangent breaks C2 so
         * curve_corner can G2-blend both ends (circular neighbors are
         * excluded from that path). */
        emit_line(ctx, 0, ARCG2_M - 1);
        return;
    }

    double mid = 0.5 * (f0 + f1);
    arc_g2_emit_range(ctx, arc, g, feed, tag, f0, mid, depth + 1, t1o, teo);
    arc_g2_emit_range(ctx, arc, g, feed, tag, mid, f1, depth + 1, t1o, teo);
}

/* Reconstruct exact PmCircle geometry (spiral + helix) from an ARC primitive.
 * Returns 0 on success, -1 degenerate (caller should fall back to a line). */
static int arc_g2_geom_init(const liscio_primitive_t *arc, arc_geom_t *g)
{
    if (fabs(arc->arc_angle) < 1e-9) return -1;
    g->C[0] = arc->cx; g->C[1] = arc->cy; g->C[2] = arc->cz;
    g->rtan[0] = arc->start.x - arc->cx;
    g->rtan[1] = arc->start.y - arc->cy;
    g->rtan[2] = arc->start.z - arc->cz;
    g->R0 = sqrt(g->rtan[0]*g->rtan[0] + g->rtan[1]*g->rtan[1] + g->rtan[2]*g->rtan[2]);
    if (g->R0 < 1e-12) return -1;
    /* rperp = n × rtan  (magnitude R0, CCW about n) */
    g->rperp[0] = arc->ny*g->rtan[2] - arc->nz*g->rtan[1];
    g->rperp[1] = arc->nz*g->rtan[0] - arc->nx*g->rtan[2];
    g->rperp[2] = arc->nx*g->rtan[1] - arc->ny*g->rtan[0];
    /* axial (helix) total = ((E−S)·n)·n */
    double h = (arc->end.x - arc->start.x)*arc->nx
             + (arc->end.y - arc->start.y)*arc->ny
             + (arc->end.z - arc->start.z)*arc->nz;
    g->hvec[0] = h*arc->nx; g->hvec[1] = h*arc->ny; g->hvec[2] = h*arc->nz;
    /* end radius (in-plane) → spiral growth */
    double ex = arc->end.x - arc->cx, ey = arc->end.y - arc->cy, ez = arc->end.z - arc->cz;
    double edotn = ex*arc->nx + ey*arc->ny + ez*arc->nz;
    double eip[3] = { ex - edotn*arc->nx, ey - edotn*arc->ny, ez - edotn*arc->nz };
    double Rend = sqrt(eip[0]*eip[0] + eip[1]*eip[1] + eip[2]*eip[2]);
    g->spiral = Rend - g->R0;
    g->phi = arc->arc_angle;
    g->s = arc->start;
    g->e = arc->end;
    return 0;
}

/* Emit one arc as fused C2 quintic(s).  t1o/teo (nullable) are forward-tangent
 * overrides at the arc's outer endpoints — the seam-kink absorption handshake
 * (both sides of a soft-kink seam receive the same averaged tangent).  Breaks
 * the G2 chain at a real corner / turn-sign flip before emitting so a wrong
 * pin never crosses a seam. */
static void arc_g2_emit_one(liscio_ctx_t *ctx,
                            const liscio_primitive_t *arc,
                            double feed, int tag,
                            const double *t1o, const double *teo)
{
    arc_geom_t g;
    if (arc_g2_geom_init(arc, &g) != 0) {
        liscio_add_line(ctx, &arc->start, &arc->end, feed, tag);
        return;
    }

    if (getenv("TP2_ARCG2_DBG")) {
        static int c = 0;
        if (c++ < 30) {
            liscio_pose_t pe; arc_pose_at(&g, 1.0, &pe);
            double dx = pe.x - arc->end.x, dy = pe.y - arc->end.y, dz = pe.z - arc->end.z;
            FILE *df = fopen("/tmp/tp2_arcg2.log", "a");
            if (df) { fprintf(df, "arc R0=%.4f spiral=%.5f phi=%.4f  reconErr(end)=%.3e kink1=%d\n",
                        g.R0, g.spiral, g.phi, sqrt(dx*dx+dy*dy+dz*dz), t1o!=NULL); fclose(df); }
        }
    }

    /* Chain-break at the arc START: a tangent discontinuity (real corner) or a
     * turn-sign flip (S-curve seam — start principal normal reverses) makes the
     * C2 pin invalid → break so the seam gets a blend join.  Uses the EFFECTIVE
     * start tangent (averaged override when absorbing a soft kink). */
    if (ctx->have_prev_emit_kappa && ctx->have_prev_emit_tan) {
        double t1[3];
        if (t1o) { t1[0]=t1o[0]; t1[1]=t1o[1]; t1[2]=t1o[2]; }
        else arc_tan_at(&g, 0.0, t1);
        double dot  = t1[0]*ctx->prev_emit_tan_x
                    + t1[1]*ctx->prev_emit_tan_y
                    + t1[2]*ctx->prev_emit_tan_z;
        double cmin = cos(ctx->cfg.corner_angle_deg * M_PI / 180.0);
        /* start principal normal (toward centre) = −unit(rtan) */
        double ndot = -(g.rtan[0]*ctx->prev_emit_norm_x
                      + g.rtan[1]*ctx->prev_emit_norm_y
                      + g.rtan[2]*ctx->prev_emit_norm_z) / g.R0;
        if (dot < cmin || ndot < 0.0) ctx->have_prev_emit_kappa = 0;
    }

    ctx->emitting_arc_g2 = 1;
    arc_g2_emit_range(ctx, arc, &g, feed, tag, 0.0, 1.0, 0, t1o, teo);
    ctx->emitting_arc_g2 = 0;
}

/* Flush the arc-chain fusion pending arc (if any): emit it with its stored
 * start-tangent override (from a previous soft-kink seam) and its OWN exact
 * end tangent (no successor to average with). */
static void arcg2_flush_pending(liscio_ctx_t *ctx)
{
    if (!ctx->has_arcg2_pend) return;
    ctx->has_arcg2_pend = 0;   /* clear FIRST: emit path may re-enter flush */
    arc_g2_emit_one(ctx, &ctx->arcg2_pend,
                    ctx->arcg2_pend_feed, ctx->arcg2_pend_tag,
                    ctx->arcg2_pend_has_t1 ? ctx->arcg2_pend_t1 : NULL, NULL);
    ctx->arcg2_pend_has_t1 = 0;
    /* The emit SCRIBBLED the shared pts[] window (quintic fit check points).
     * Callers like liscio_add_line flush the pending arc BEFORE appending
     * their own point — without this reset those 16 arc samples would be
     * mistaken for window waypoints (measured: pathlen +200mm, ghost fits,
     * 4e10 jerk).  Empty window: the next G1 line starts a fresh run. */
    ctx->n_pts = 0;
}

/* Arc-chain fusion entry: hold the arc until its successor is known, so a soft
 * seam kink (tangent break ≤ corner_angle_deg — thomam: 26 seams of 1-19°) is
 * ABSORBED by handing both sides the averaged seam tangent (fit tolerance gate
 * rejects kinks too big to absorb → bisect → native).  Without this every kink
 * seam falls to a corner blend whose dκ/ds-dominated cap (~21 mm/s) wrecks
 * feed holding (thomam fused: 51% of time stuck at 20-45 mm/s).  A feed change
 * or a hard corner emits the pending arc as-is (seam falls to blend/kink cap). */
static void liscio_arc_g2_add(liscio_ctx_t *ctx,
                              const liscio_primitive_t *arc,
                              double feed, int tag)
{
    if (ctx->has_arcg2_pend) {
        arc_geom_t gp, gn;
        int soft = 0;
        double tm[3] = {0,0,0};
        if (fabs(feed - ctx->arcg2_pend_feed) < 1e-9
            && arc_g2_geom_init(&ctx->arcg2_pend, &gp) == 0
            && arc_g2_geom_init(arc, &gn) == 0) {
            double ta[3], tb[3];
            arc_tan_at(&gp, 1.0, ta);   /* pending exit tangent  */
            arc_tan_at(&gn, 0.0, tb);   /* new     entry tangent */
            double dot = ta[0]*tb[0] + ta[1]*tb[1] + ta[2]*tb[2];
            double cmin = cos(ctx->cfg.corner_angle_deg * M_PI / 180.0);
            if (dot > cmin) {
                tm[0] = ta[0] + tb[0]; tm[1] = ta[1] + tb[1]; tm[2] = ta[2] + tb[2];
                double m = sqrt(tm[0]*tm[0] + tm[1]*tm[1] + tm[2]*tm[2]);
                if (m > 1e-12) { tm[0]/=m; tm[1]/=m; tm[2]/=m; soft = 1; }
            }
        }
        /* Emit the pending arc: averaged exit tangent on a soft seam, own exact
         * tangent otherwise. */
        int had_t1 = ctx->arcg2_pend_has_t1;
        double t1s[3] = { ctx->arcg2_pend_t1[0], ctx->arcg2_pend_t1[1], ctx->arcg2_pend_t1[2] };
        ctx->has_arcg2_pend = 0;
        ctx->arcg2_pend_has_t1 = 0;
        arc_g2_emit_one(ctx, &ctx->arcg2_pend,
                        ctx->arcg2_pend_feed, ctx->arcg2_pend_tag,
                        had_t1 ? t1s : NULL, soft ? tm : NULL);
        /* The new arc inherits the averaged tangent as ITS start override. */
        if (soft) {
            ctx->arcg2_pend_t1[0] = tm[0];
            ctx->arcg2_pend_t1[1] = tm[1];
            ctx->arcg2_pend_t1[2] = tm[2];
            ctx->arcg2_pend_has_t1 = 1;
        }
    }
    ctx->arcg2_pend      = *arc;
    ctx->arcg2_pend_feed = feed;
    ctx->arcg2_pend_tag  = tag;
    ctx->has_arcg2_pend  = 1;
}

int liscio_add_arc(liscio_ctx_t *ctx,
                    const liscio_primitive_t *arc,
                    double feedrate,
                    int tag_line_id)
{
    if (!ctx || !arc) return -1;
    ctx->stats.input_count++;

    /* Subdivision mode: split arc into N straight sub-segments and
     * push them through the G1 pipeline.  Lets adjacent G1+arc runs
     * fold into longer Bezier/Spline primitives at the cost of
     * arc-precision.  Default cfg.arc_subseg_samples = 0 keeps the
     * lossless passthrough below. */
    if (ctx->cfg.arc_subseg_samples > 0) {
        int N = ctx->cfg.arc_subseg_samples;
        /* Plane basis: u = unit start-radial, v = normal × u (so [u,v,n]
         * right-handed; arc CCW about n traces cos(a)·u + sin(a)·v). */
        double ux = arc->start.x - arc->cx;
        double uy = arc->start.y - arc->cy;
        double uz = arc->start.z - arc->cz;
        double um = sqrt(ux*ux + uy*uy + uz*uz);
        if (um < 1e-15)
            return liscio_add_line(ctx, &arc->start, &arc->end,
                                    feedrate, tag_line_id);
        ux /= um; uy /= um; uz /= um;
        double vx = arc->ny*uz - arc->nz*uy;
        double vy = arc->nz*ux - arc->nx*uz;
        double vz = arc->nx*uy - arc->ny*ux;
        /* Helical off-plane displacement = (end - start) projected on n. */
        double dx = arc->end.x - arc->start.x;
        double dy = arc->end.y - arc->start.y;
        double dz = arc->end.z - arc->start.z;
        double off_delta = dx*arc->nx + dy*arc->ny + dz*arc->nz;
        double ang = arc->arc_angle;
        double r   = arc->radius;
        liscio_pose_t prev = arc->start;
        for (int s = 1; s <= N; s++) {
            double t = (double)s / N;
            liscio_pose_t p;
            if (s == N) {
                p = arc->end;   /* exact endpoint, no float drift */
            } else {
                double a = t * ang;
                double off = t * off_delta;
                p.x = arc->cx + r*(cos(a)*ux + sin(a)*vx) + off*arc->nx;
                p.y = arc->cy + r*(cos(a)*uy + sin(a)*vy) + off*arc->ny;
                p.z = arc->cz + r*(cos(a)*uz + sin(a)*vz) + off*arc->nz;
                p.a = arc->start.a + t * (arc->end.a - arc->start.a);
                p.b = arc->start.b + t * (arc->end.b - arc->start.b);
                p.c = arc->start.c + t * (arc->end.c - arc->start.c);
                p.u = arc->start.u + t * (arc->end.u - arc->start.u);
                p.v = arc->start.v + t * (arc->end.v - arc->start.v);
                p.w = arc->start.w + t * (arc->end.w - arc->start.w);
            }
            liscio_add_line(ctx, &prev, &p, feedrate, tag_line_id);
            prev = p;
        }
        return 0;
    }

    /* Arc-chain G2 fusion: emit the arc as analytic C2 quintic(s) so a run of
     * tangent-continuous arcs becomes one curvature-continuous spline chain
     * (no per-seam curvature step); soft seam kinks are absorbed via the
     * one-arc pending buffer (averaged seam tangent).  Drain any pending G1
     * line window first. */
    if (ctx->cfg.arc_g2fuse) {
        flush_window_as_primitive(ctx);
        liscio_arc_g2_add(ctx, arc, feedrate, tag_line_id);
        ctx->pts[0]  = arc->end;   /* seed window for a following G1 run */
        ctx->feed[0] = feedrate;
        ctx->tag[0]  = tag_line_id;
        ctx->n_pts   = 1;
        return 0;
    }

    /* G2/G3 passes through: flush pending G1 window first. */
    flush_window_as_primitive(ctx);

    /* Arc-to-helix LSQ mode: buffer multiple G2/G3 arcs, then test
     * whether they jointly fit a HELIX primitive (analytic, no point
     * sampling).  Coarse admission: same axis / radius / direction /
     * endpoint continuity.  Final decision made on flush by
     * flush_arc_helix_buffer() (LSQ residual gate). */
    if (ctx->cfg.arc_to_helix_max > 0) {
        int cap = ctx->cfg.arc_to_helix_max;
        if (cap > LISCIO_ARC_HELIX_BUF) cap = LISCIO_ARC_HELIX_BUF;
        if (ctx->arc_buf_n == 0) {
            /* Start a new run. */
            ctx->arc_buf[0]          = *arc;
            ctx->arc_buf[0].tag_line_first = tag_line_id;
            ctx->arc_buf[0].tag_line_last  = tag_line_id;
            ctx->arc_buf_n           = 1;
            ctx->arc_buf_tag_first   = tag_line_id;
            ctx->arc_buf_tag_last    = tag_line_id;
            ctx->arc_buf_feed        = feedrate;
        } else {
            const liscio_primitive_t *prev = &ctx->arc_buf[ctx->arc_buf_n - 1];
            double radius_tol = ctx->cfg.tol_xyz;
            double xy_tol     = ctx->cfg.tol_xyz;
            /* Same-F gate: a merged primitive carries one feedrate. */
            int feed_ok = !ctx->cfg.merge_same_feed_only
                          || feed_match(ctx->arc_buf_feed, feedrate);
            if (feed_ok
                && arc_helix_admissible(prev, arc, radius_tol, xy_tol)
                && ctx->arc_buf_n < cap)
            {
                ctx->arc_buf[ctx->arc_buf_n]                 = *arc;
                ctx->arc_buf[ctx->arc_buf_n].tag_line_first  = tag_line_id;
                ctx->arc_buf[ctx->arc_buf_n].tag_line_last   = tag_line_id;
                ctx->arc_buf_n++;
                ctx->arc_buf_tag_last = tag_line_id;
            } else {
                /* Run break or buffer full: flush + restart with new arc. */
                flush_arc_helix_buffer(ctx);
                ctx->arc_buf[0]                 = *arc;
                ctx->arc_buf[0].tag_line_first  = tag_line_id;
                ctx->arc_buf[0].tag_line_last   = tag_line_id;
                ctx->arc_buf_n                  = 1;
                ctx->arc_buf_tag_first          = tag_line_id;
                ctx->arc_buf_tag_last           = tag_line_id;
                ctx->arc_buf_feed               = feedrate;
            }
        }
        /* Seed G1 window with arc's end for any following G1 run. */
        ctx->pts[0]  = arc->end;
        ctx->feed[0] = feedrate;
        ctx->tag[0]  = tag_line_id;
        ctx->n_pts   = 1;
        return 0;
    }

    /* Arc-merge mode: try to extend a buffered prior arc.  Same-F gate: a
     * merged primitive carries one feedrate, so an F change breaks the run
     * (flush pending, buffer the new arc) exactly like a geometry break. */
    if (ctx->cfg.arc_merge_mode != LISCIO_ARC_MERGE_OFF) {
        if (ctx->has_pending_arc &&
            (!ctx->cfg.merge_same_feed_only
             || feed_match(ctx->pending_feed, feedrate)) &&
            arcs_mergeable(&ctx->pending_arc, arc,
                           ctx->cfg.arc_merge_xy_tol,
                           ctx->cfg.arc_merge_pitch_tol))
        {
            /* Extend pending: keep start/center/normal/radius, accumulate
             * arc_angle, replace end with new arc's end. */
            ctx->pending_arc.arc_angle += arc->arc_angle;
            ctx->pending_arc.end        = arc->end;
            ctx->pending_tag_last       = tag_line_id;
            ctx->pending_n_absorbed++;
            /* Seed window for next G1 run after eventual flush. */
            ctx->pts[0]  = arc->end;
            ctx->feed[0] = feedrate;
            ctx->tag[0]  = tag_line_id;
            ctx->n_pts   = 1;
            return 0;
        }
        /* Cannot extend — flush whatever was pending, buffer the new one. */
        flush_pending_arc(ctx);
        ctx->pending_arc        = *arc;
        ctx->pending_feed       = feedrate;
        ctx->pending_tag        = tag_line_id;
        ctx->pending_tag_last   = tag_line_id;
        ctx->pending_n_absorbed = 1;
        ctx->has_pending_arc    = 1;
        ctx->pts[0]  = arc->end;
        ctx->feed[0] = feedrate;
        ctx->tag[0]  = tag_line_id;
        ctx->n_pts   = 1;
        return 0;
    }

    /* OFF mode: legacy direct passthrough — 1 arc → 1 ARC primitive. */
    emit_arc_primitive(ctx, arc, feedrate, tag_line_id, tag_line_id, 1);

    /* Seed window with arc's end for next line run. */
    ctx->pts[0]  = arc->end;
    ctx->feed[0] = feedrate;
    ctx->tag[0]  = tag_line_id;
    ctx->n_pts   = 1;
    return 0;
}

/* ------------------------------------------------------------------ */
void liscio_flush(liscio_ctx_t *ctx)
{
    if (!ctx) return;
    ctx->in_cluster = 0;
    if (ctx->n_pts >= 2)
        flush_window_as_primitive(ctx);
    ctx->n_pts = 0;
    flush_pending_arc(ctx);
    flush_arc_helix_buffer(ctx);
    arcg2_flush_pending(ctx);
}

/* Flush ONLY the pending G1 line window (arc-chain fusion helper).  The tp2
 * bridge calls this before switching modal from lines to an arc so the fitted
 * lines emit under the correct (line) modal — WITHOUT emitting the fusion
 * pending arc, whose seam-tangent handshake with the NEXT arc must survive
 * across consecutive liscio_add_arc calls. */
void liscio_flush_lines(liscio_ctx_t *ctx)
{
    if (!ctx) return;
    ctx->in_cluster = 0;
    if (ctx->n_pts >= 2)
        flush_window_as_primitive(ctx);
    ctx->n_pts = 0;
}

/* ------------------------------------------------------------------ */
/* Internal: drain pending state before injecting a passthrough event. */
static void drain_pending(liscio_ctx_t *ctx)
{
    ctx->in_cluster = 0;           /* an injected event breaks the cluster */
    if (ctx->n_pts >= 2)
        flush_window_as_primitive(ctx);
    ctx->n_pts = 0;
    flush_pending_arc(ctx);
    flush_arc_helix_buffer(ctx);
    arcg2_flush_pending(ctx);
    ctx->have_prev_emit_tan = 0;   /* tangent broken across event */
}

int liscio_add_rapid(liscio_ctx_t *ctx,
                      const liscio_pose_t *start,
                      const liscio_pose_t *end,
                      int tag_line_id)
{
    if (!ctx || !start || !end) return -1;
    ctx->stats.input_count++;
    drain_pending(ctx);
    if (!ctx->cb) return 0;

    liscio_primitive_t prim;
    memset(&prim, 0, sizeof(prim));
    prim.type           = LISCIO_PRIM_RAPID;
    prim.start          = *start;
    prim.end            = *end;
    prim.feedrate       = 0.0;            /* rapid: TP uses max-velocity */
    prim.tag_line_first = tag_line_id;
    prim.tag_line_last  = tag_line_id;
    prim.tag_line_id    = tag_line_id;
    prim.n_absorbed     = 1;
    /* Tangent = chord direction (informational, TP usually ignores). */
    fill_tangent_from_points(&prim, &prim.start, &prim.end);
    prim.tol_xyz = ctx->cfg.tol_xyz;
    prim.tol_abc = ctx->cfg.tol_abc;
    ctx->cb(&prim, ctx->cb_user);
    /* RAPID emits don't seed cross-window G1 — tangent stays cleared. */
    return 0;
}

int liscio_emit_stop(liscio_ctx_t *ctx,
                      liscio_stop_reason_t reason,
                      double dwell_seconds,
                      const liscio_pose_t *position,
                      int tag_line_id)
{
    if (!ctx) return -1;
    ctx->stats.input_count++;
    drain_pending(ctx);
    if (!ctx->cb) return 0;

    liscio_primitive_t prim;
    memset(&prim, 0, sizeof(prim));
    prim.type               = LISCIO_PRIM_STOP;
    prim.stop_reason        = (int)reason;
    prim.stop_dwell_seconds = (reason == LISCIO_STOP_DWELL) ? dwell_seconds : 0.0;
    if (position) { prim.start = *position; prim.end = *position; }
    prim.tag_line_first     = tag_line_id;
    prim.tag_line_last      = tag_line_id;
    prim.tag_line_id        = tag_line_id;
    prim.n_absorbed         = 1;
    prim.tol_xyz            = ctx->cfg.tol_xyz;
    prim.tol_abc            = ctx->cfg.tol_abc;
    ctx->cb(&prim, ctx->cb_user);
    return 0;
}

/* ------------------------------------------------------------------ */
int liscio_primitive_line_at(const liscio_primitive_t *prim, double t)
{
    if (!prim) return 0;
    int lo = prim->tag_line_first;
    int hi = prim->tag_line_last;
    if (lo == 0 && hi == 0) lo = hi = prim->tag_line_id;
    if (lo == hi) return lo;
    if (t <= 0.0) return lo;
    if (t >= 1.0) return hi;
    /* Linear interpolation; rounds to nearest line.  G-code line numbers
     * are typically contiguous so this matches actual execution order;
     * for non-contiguous numbering (O-word jumps) it is still monotone
     * within an absorbed run. */
    double f = (double)lo + t * (double)(hi - lo);
    return (int)(f + 0.5);
}

/* ------------------------------------------------------------------ */
void liscio_reset(liscio_ctx_t *ctx)
{
    if (!ctx) return;
    ctx->n_pts = 0;
    ctx->have_prev_emit_tan = 0;
    ctx->have_prev_emit_kappa = 0;
    ctx->prev_emit_kappa = 0.0;
    ctx->prev_emit_norm_x = ctx->prev_emit_norm_y = ctx->prev_emit_norm_z = 0.0;
    ctx->has_pending_arc = 0;
    ctx->pending_n_absorbed = 0;
    ctx->in_cluster = 0;
    /* Drop fusion pendings without emit — program boundary / abort. */
    ctx->arc_buf_n = 0;
    ctx->has_arcg2_pend = 0;
    ctx->arcg2_pend_has_t1 = 0;
    ctx->emitting_arc_g2 = 0;
    memset(&ctx->stats, 0, sizeof(ctx->stats));
}

void liscio_set_path_tol(liscio_ctx_t *ctx, double tol_xyz, double tol_arcfit)
{
    if (!ctx) return;
    if (tol_xyz < 0.0) tol_xyz = 0.0;
    ctx->cfg.tol_xyz = tol_xyz;
    if (tol_arcfit > 0.0)
        ctx->cfg.tol_arcfit = tol_arcfit;
    else
        ctx->cfg.tol_arcfit = tol_xyz;
    ctx->cfg.collinear_tol = tol_xyz;
    if (ctx->cfg.tol_uvw <= 0.0)
        ctx->cfg.tol_uvw = tol_xyz;
    if (ctx->cfg.collapse_cusp_tol <= 0.0 ||
        ctx->cfg.collapse_cusp_tol > tol_xyz)
        ctx->cfg.collapse_cusp_tol = tol_xyz;
    /* Mid-tight + flower-loose G64: G1-chain kills residual θ.  6.6.126
     * fused_c2 |Δκ| samples did not touch flower Y=73.2 (P=0.05 was above
     * the old 0.03 HI → g2_quintic κ=0 ends, 3.6° left on LINE seams).
     * 6.6.136: LO 0.0025→0.001 so 1_1001/2_1001 same-dir leftover
     * wraps can G1-chain (n_pts>=4).  Do not stack with g2_quintic
     * (P=0.005 + both → SER-40 +203M). */
    ctx->cfg.enable_g1_chain = (tol_xyz >= 0.001 && tol_xyz <= 0.05)
                               && !getenv("TP2_NO_G1CHAIN");
    /* Only looser than flower P=0.05.  TP2_G2ZERO forces on (A/B). */
    ctx->cfg.enable_g2_quintic = !getenv("TP2_NO_G2ZERO")
        && (getenv("TP2_G2ZERO") || (tol_xyz > 0.05));
}
