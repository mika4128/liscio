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
* Description: liscio - pure C99 geometry preprocessor library for
*              CNC trajectory planners.  Absorbs CAM-generated
*              micro-line-segment sequences into fitted primitives
*              (LINE / ARC / cubic Bezier / rational Bezier) before
*              look-ahead.  9D (XYZ + ABC + UVW), RT-friendly
*              (bounded memory, no alloc after create()).
*              Industry-standard 2-stage pipeline:
*                G-code → [liscio] → look-ahead → interpolator
*              Comparable to Fanuc Nano Smoothing / Siemens COMPCURV.
* Author:      杨阳 (Yang Yang) <mika-net@outlook.com>
* License:     MIT (SPDX-License-Identifier: MIT)
* Copyright (c) 2026 杨阳 (Yang Yang)
********************************************************************/

#ifndef LISCIO_H
#define LISCIO_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* 9D pose: LinuxCNC EmcPose-compatible layout but no LinuxCNC header
 * dependency (adapter in user code). */
typedef struct {
    double x, y, z;  /* tran */
    double a, b, c;  /* rotary */
    double u, v, w;  /* additional linear */
} liscio_pose_t;

typedef enum {
    LISCIO_PRIM_LINE   = 0,  /* fitted G1 group → straight line       */
    LISCIO_PRIM_ARC    = 1,  /* fitted G1 group → planar arc (Kasa)   */
    LISCIO_PRIM_BEZIER = 2,  /* cubic 9D Bezier             */
    LISCIO_PRIM_SPLINE = 3,  /* 9D cubic B-spline / NURBS   */
    LISCIO_PRIM_HELIX  = 4,  /* multi-loop helical arc      */
    /* G-code-style event primitives.  These are NOT fits; they
     * are passthrough events that the caller (rs274 bridge) injects so
     * downstream TP receives the full G-code intent in order.  Geometry
     * fields (start, end) are valid; max_dev is 0 by definition. */
    LISCIO_PRIM_RAPID  = 5,  /* G0 rapid traversal (start, end)       */
    LISCIO_PRIM_STOP   = 6,  /* M0/M1/M2/M30/M6/G4/G92/...
                                  reason in stop_reason; G4 dwell sec
                                  in stop_dwell_seconds; position in
                                  start (= end). */
} liscio_prim_type_t;

/* Reasons for LISCIO_PRIM_STOP — caller maps from canon events.  See
 * docs/architecture.md "STOP / run boundary" + results/events.md for the
 * full LinuxCNC canon → liscio mapping table. */
typedef enum {
    LISCIO_STOP_NONE            = 0,
    LISCIO_STOP_PROGRAM         = 1,  /* M0 mandatory pause            */
    LISCIO_STOP_OPTIONAL        = 2,  /* M1 optional pause             */
    LISCIO_STOP_PROGRAM_END     = 3,  /* M2                            */
    LISCIO_STOP_PROGRAM_END_REWIND = 4, /* M30                         */
    LISCIO_STOP_TOOL_CHANGE     = 5,  /* M6                            */
    LISCIO_STOP_DWELL           = 6,  /* G4 P<sec>                     */
    LISCIO_STOP_SET_OFFSETS     = 7,  /* G92 / G92.1 / G10             */
    LISCIO_STOP_TOOL_LENGTH     = 8,  /* G43 / G49                     */
    LISCIO_STOP_COORD_SYSTEM    = 9,  /* G54-G59 P switch              */
    LISCIO_STOP_USER            = 10, /* explicit caller mark          */
} liscio_stop_reason_t;

/* Max # control points for B-spline primitive.  Inline arrays
 * keep the primitive RT-safe (no heap allocation).  Typical composite
 * CAM curve fits within 8-12 CPs; 16 gives headroom. */
#ifndef LISCIO_BSPLINE_MAX_CTRL
#define LISCIO_BSPLINE_MAX_CTRL 16
#endif

typedef struct {
    liscio_prim_type_t type;
    liscio_pose_t start;
    liscio_pose_t end;

    /* Arc-specific (valid when type == LISCIO_PRIM_ARC).
     * Center + normal define the plane. turn encodes direction + wrap:
     * positive = CCW, negative = CW, |turn| = full-turns + 1. */
    double cx, cy, cz;        /* center */
    double nx, ny, nz;        /* plane normal (unit) */
    double radius;            /* xyz-plane radius */
    double arc_angle;         /* signed traversed angle (rad) */
    double pitch;             /* HELIX: mm/rad along +n; 0 for planar ARC */

    /* Cubic Bezier control points (when type == LISCIO_PRIM_BEZIER).
     * P0 == start, P3 == end (redundant with above, kept for clarity). */
    liscio_pose_t P0, P1, P2, P3;
    double bezier_max_deviation;

    /* Rational weights for NURBS / rational Bezier.
     * For non-rational curve (default): all 1.0.
     * For arc-as-rational-Bezier: (1, cos(θ/2), cos(θ/2), 1). */
    double w0, w1, w2, w3;

    /* Tangent vectors at endpoints (XYZ only; unit-length).
     * Computed from control points; used by downstream blend code to
     * establish tangent-continuity between adjacent primitives.
     * tangent_start = 3·(P1 - P0)/|3·(P1-P0)|; mirror at end. */
    double tan_start_x, tan_start_y, tan_start_z;
    double tan_end_x,   tan_end_y,   tan_end_z;

    /* Per-primitive tolerance (G64 P/Q hint forwarded from G-code).
     * 0 or negative means "use planner default". */
    double tol_xyz, tol_abc;

    /* NURBS / B-spline primitive fields (when type == LISCIO_PRIM_SPLINE).
     * Fixed capacity LISCIO_BSPLINE_MAX_CTRL; unused slots zeroed. */
    int degree;                                /* 3 for cubic */
    int n_ctrl;                                /* active control points */
    liscio_pose_t ctrl_pts[LISCIO_BSPLINE_MAX_CTRL];
    double        weights [LISCIO_BSPLINE_MAX_CTRL];
    double        knots   [LISCIO_BSPLINE_MAX_CTRL + 4];  /* clamped, size n_ctrl+degree+1 */
    double        spline_max_deviation;

    /* Source metadata (passed through).
     *
     * Line tracking — when many G1 segments fold into one primitive,
     * `tag_line_first` / `tag_line_last` carry the line range absorbed.
     * Use `liscio_primitive_line_at(prim, t)` for parameter→line lookup
     * (UI highlight, run-from-line, error reporting).
     *
     * `tag_line_id` is preserved as an alias of `tag_line_first` for
     * backward-compat with single-line consumers. */
    double feedrate;          /* requested F */
    int    tag_line_id;       /* alias of tag_line_first (back-compat) */
    int    tag_line_first;    /* first absorbed G-code line */
    int    tag_line_last;     /* last  absorbed G-code line  */
    int    n_absorbed;        /* # input segs merged into this primitive */

    /* STOP-event metadata (valid when type == LISCIO_PRIM_STOP).
     * For other primitive types these are zero. */
    int    stop_reason;             /* see liscio_stop_reason_t */
    double stop_dwell_seconds;      /* > 0 when stop_reason == DWELL */
} liscio_primitive_t;

typedef void (*liscio_emit_cb)(const liscio_primitive_t *prim, void *user);

/* Tolerance configuration: G64 P applies to xyz, G64 Q applies to abc.
 * uvw shares xyz tol by default. */
typedef struct {
    double tol_xyz;       /* mm */
    double tol_abc;       /* deg */
    double tol_uvw;       /* mm (<=0 uses tol_xyz) */
    int    max_window;    /* max points buffered before forced emit (default 128) */
    int    min_arc_pts;   /* min consecutive pts to trigger arc fit (default 5) */
    double min_arc_radius;  /* lower bound on detected radius (mm, default 0.001) */
    double max_arc_radius;  /* upper bound; beyond this treat as line (default 1e6) */

    /* Corner detection: when the XYZ direction between two consecutive
     * input segments rotates by more than corner_angle_deg, treat as a
     * tangent discontinuity and split the window there (flush before
     * accepting the new segment).  Industry default: G641-equivalent
     * ~15-45°. 0 disables corner splitting. */
    double corner_angle_deg;

    /* Collinear tolerance: segments whose endpoints' perpendicular
     * offset from the running line direction < collinear_tol are
     * considered same-line.  Used by streaming line absorption. */
    double collinear_tol;

    /* Enable cubic Bezier LSQ fitting (default 1). Try after arc fit
     * fails; captures smooth CAM tool-path splines.  Disable for strict
     * geometric preservation. */
    int    enable_bezier;

    /* Enable composite Bezier recursive split (default 1). When single
     * cubic Bezier fit exceeds tol, split the window at mid-chord and
     * fit each half (recursive up to 8 levels / 256 sub-curves). */
    int    enable_composite_bezier;

    /* Enable free-form cubic B-spline LSQ fit (default 1).
     * Tried after composite Bezier; captures complex smooth curves
     * more compactly (one primitive, any M ctrl points). */
    int    enable_bspline;

    /* Initial # control points for B-spline fit attempt (default 6).
     * If fit fails, window is split and each half attempts smaller
     * M via the composite/split code path. */
    int    bspline_init_ctrl;

    /* G2/G3 input handling.  liscio_add_arc() default behaviour is
     * geometric passthrough (max_dev=0, one ARC primitive per input).
     *
     * Set arc_subseg_samples > 0 to instead subdivide each arc into N
     * straight sub-segments and feed them through the G1 LSQ pipeline.
     * Useful when many small arcs (e.g. CAM biarc fillets) could fold
     * into longer Bezier/Spline primitives, OR when CAM helical
     * chord-approximation (multiple G2/G3 with slightly drifting
     * radius) should fold into a single HELIX primitive — the helix9
     * fitter only runs on the G1 path, so non-zero subseg lets it see
     * G2/G3 helical runs.
     *
     *   0  = passthrough (default, lossless, 1 arc → 1 ARC primitive)
     *   ≥1 = subdivide each arc into N sub-segments, fed via add_line.
     * Recommended values: 4–16; aim for sub-segment chord length ≤ tol.
     * For helical CAM data (e.g. thomam.ngc lathe-thread cuts), N=8
     * with tol_xyz=0.05 typically yields 4-10× compression with several
     * HELIX primitives covering full multi-turn runs. */
    int    arc_subseg_samples;

    /* Adjacent-arc merge.  When successive G2/G3 inputs share
     * the same plane (center XY, normal, radius, helical pitch within
     * arc_merge_xy_tol / arc_merge_pitch_tol), they can be folded into
     * a single primitive whose arc_angle accumulates.  CAM helix output
     * (e.g. thomam.ngc) compresses ~5-10× under this mode.
     *
     * Configure via liscio_cfg_set_arc_merge() (preferred, see below)
     * rather than touching these fields directly.  Geometric precision
     * is preserved (max_dev = 0). */
    int    arc_merge_mode;            /* see liscio_arc_merge_mode_t */
    double arc_merge_xy_tol;          /* mm; center / radius / start point */
    double arc_merge_pitch_tol;       /* mm/rad; helical pitch consistency */

    /* Direct ARC → HELIX least-squares detection (analytic, no point
     * sampling).  liscio_add_arc accumulates up to N consecutive G2/G3
     * arcs whose axis / radius / center align approximately, then
     * runs a LSQ helix fit on the buffered ARC parameters (centers
     * along common axis, axial advance ∝ cumulative angle).  When the
     * fit's max_dev (radius + center-to-axis offset) is within tol_xyz,
     * the whole run emits as one HELIX primitive; otherwise the buffer
     * flushes as individual ARCs.
     *
     * Unlike `arc_subseg_samples`, this preserves the original arc
     * representation (no sampling) and only writes a HELIX when the
     * geometric error is genuinely below tolerance.  Best for CAM
     * lathe-thread / chord-approximated helical cuts (e.g. thomam.ngc).
     *
     *   0   = disable (force per-arc passthrough)
     *   N>0 = max arcs to buffer per helix detection (default 16) */
    int    arc_to_helix_max;

    /* Corner detection mode (lookahead).
     *
     *   IMMEDIATE = legacy: any segment-direction change > corner_angle_deg
     *               flushes the window AT ONCE.  Hard corners and soft
     *               curvature transitions get the same hard split.  Higher
     *               tan_max but predictable.
     *
     *   LOOKAHEAD = two-tier:
     *               soft (corner_angle_deg ≤ angle ≤ corner_hard_angle_deg)
     *                  → defer flush, let composite Bezier fit absorb it
     *               hard (angle > corner_hard_angle_deg)
     *                  → flush immediately (preserve real-corner intent)
     *               Soft curves get fitted with G1 tangent continuity →
     *               substantial tan_mean drop on flower-style CAM data.
     *
     * Configure via liscio_cfg_set_corner_detection(); fields below are
     * for read-only inspection.  Defaults via liscio_cfg_default():
     *   corner_detection_mode  = LISCIO_CORNER_IMMEDIATE
     *   corner_hard_angle_deg  = 60.0  (only consulted in LOOKAHEAD) */
    int    corner_detection_mode;
    double corner_hard_angle_deg;
} liscio_cfg_t;

typedef enum {
    LISCIO_CORNER_IMMEDIATE = 0,   /* default — back-compatible */
    LISCIO_CORNER_LOOKAHEAD = 1,
} liscio_corner_mode_t;

/* Configure corner detection.  Pass <=0 to keep current/default values for
 * either threshold.  In IMMEDIATE mode only soft_angle_deg is consulted;
 * hard_angle_deg becomes meaningful only in LOOKAHEAD mode. */
void liscio_cfg_set_corner_detection(liscio_cfg_t *cfg,
                                       liscio_corner_mode_t mode,
                                       double soft_angle_deg,
                                       double hard_angle_deg);

/* Arc-merge modes for cfg.arc_merge_mode.  Use the setter below rather
 * than assigning the enum directly so default tols stay coherent. */
typedef enum {
    LISCIO_ARC_MERGE_OFF   = 0,   /* default — 1 input arc = 1 ARC out */
    LISCIO_ARC_MERGE_ARC   = 1,   /* fold contiguous arcs into single ARC
                                     primitive (arc_angle may exceed 2π) */
    LISCIO_ARC_MERGE_HELIX = 2,   /* same merge logic, but emit HELIX
                                     primitive type when |arc_angle| > 2π
                                     (multi-loop wrap), ARC otherwise */
} liscio_arc_merge_mode_t;

/* Configure arc-merge mode + sensible default tolerances.
 *   mode      = LISCIO_ARC_MERGE_{OFF,ARC,HELIX}
 *   xy_tol    = mm.  ≤ 0 → 1e-6.  Center XY, radius, helical Z step.
 *   pitch_tol = mm/rad.  ≤ 0 → 1e-6.  Helical pitch consistency.
 * Call after liscio_cfg_default(). */
void liscio_cfg_set_arc_merge(liscio_cfg_t *cfg,
                                liscio_arc_merge_mode_t mode,
                                double xy_tol,
                                double pitch_tol);

/* Run-time statistics (read-only to caller). */
typedef struct {
    long input_count;         /* # liscio_add_line calls */
    long emitted_line;        /* # LINE primitives emitted */
    long emitted_arc;         /* # ARC primitives emitted */
    long emitted_spline;      /* # SPLINE primitives emitted */
    long absorbed_total;      /* sum of all primitive.n_absorbed */
    double max_arc_deviation; /* max deviation observed during successful arc fits */
} liscio_stats_t;

/* Opaque context (heap-allocated in liscio_create). */
typedef struct liscio_ctx liscio_ctx_t;

/* Default-initialize config. */
void liscio_cfg_default(liscio_cfg_t *cfg);

/* Create context. Returns NULL on alloc failure.
 * ctx owns its internal waypoint buffer — no further allocation afterward. */
liscio_ctx_t *liscio_create(const liscio_cfg_t *cfg);
void           liscio_destroy(liscio_ctx_t *ctx);

/* Register emit callback. Called once per emitted primitive. */
void liscio_set_callback(liscio_ctx_t *ctx, liscio_emit_cb cb, void *user);

/* Feed one line segment (start→end, 9D). feedrate/tag are passed through.
 * Returns 0 on success, -1 on error. Context may emit 0 or more primitives
 * during this call (via callback). */
int liscio_add_line(liscio_ctx_t *ctx,
                     const liscio_pose_t *start,
                     const liscio_pose_t *end,
                     double feedrate,
                     int tag_line_id);

/* Pass-through: caller has a G2/G3 arc primitive; emit it directly after
 * flushing any buffered line accumulation. */
int liscio_add_arc(liscio_ctx_t *ctx,
                    const liscio_primitive_t *arc,
                    double feedrate,
                    int tag_line_id);

/* Inject a G0 rapid traversal into the primitive stream.  Flushes any
 * pending G1 window and arc-merge buffer first (so the prior cutting
 * run completes), then emits a LISCIO_PRIM_RAPID with the given endpoints.
 *
 * liscio does not "compress" rapids — they're passthrough events for
 * the downstream TP to schedule traversal motion (max-velocity, no
 * precision constraint). */
int liscio_add_rapid(liscio_ctx_t *ctx,
                      const liscio_pose_t *start,
                      const liscio_pose_t *end,
                      int tag_line_id);

/* Inject a STOP event (M0/M1/M2/M30/M6/G4/G92/...) into the stream.
 * Flushes pending state, then emits LISCIO_PRIM_STOP carrying the
 * reason + dwell.  `position` is optional (NULL = use last known
 * tool position from the prior emit); pass when the caller knows the
 * exact pose (e.g. some events are tied to current cmd-position).
 *
 *   reason          = which canon event triggered the stop
 *   dwell_seconds   = > 0 only for LISCIO_STOP_DWELL (G4 P<sec>)
 *   tag_line_id     = original G-code line number for traceability */
int liscio_emit_stop(liscio_ctx_t *ctx,
                      liscio_stop_reason_t reason,
                      double dwell_seconds,
                      const liscio_pose_t *position,
                      int tag_line_id);

/* Flush: emit any buffered points as primitive(s). Call at program end
 * or when upstream context breaks (tool change, dwell, rapid). */
void liscio_flush(liscio_ctx_t *ctx);

/* Map a primitive parameter t∈[0,1] to the original G-code line number.
 * Used by UI/debug:  "what line is currently executing?"  Linear
 * interpolation between tag_line_first..tag_line_last; for single-line
 * primitives both are equal and the answer is constant.
 * Out-of-range t is clamped. */
int liscio_primitive_line_at(const liscio_primitive_t *prim, double t);

/* Reset internal state (drop buffered points, zero stats) without
 * emitting.  Use for program restart. */
void liscio_reset(liscio_ctx_t *ctx);

/* Get runtime stats snapshot. */
void liscio_get_stats(const liscio_ctx_t *ctx, liscio_stats_t *out);

/* ---------- NURBS / rational Bezier utilities ---------- */

/* Convert an ARC primitive to exact rational cubic Bezier with weights
 * (1, cos(θ/2), cos(θ/2), 1). For arc_angle ≤ 2π/3 (120°), conversion
 * is geometrically exact. Wider arcs should be split upstream.
 * out receives the new primitive (type LISCIO_PRIM_BEZIER). */
int liscio_arc_to_rational_bezier(const liscio_primitive_t *arc,
                                    liscio_primitive_t *out);

/* Evaluate a (rational) cubic Bezier at t ∈ [0,1]. Uses primitive's
 * weights; reduces to non-rational if all are 1.0. */
void liscio_rbezier_eval(const liscio_primitive_t *prim, double t,
                          double *out_x, double *out_y, double *out_z);

/* ---------- B-spline utilities ---------- */

/* Evaluate a (rational) cubic B-spline primitive at t ∈ [0,1].
 * Reduces to non-rational B-spline when all weights are 1.0. */
void liscio_bspline_eval(const liscio_primitive_t *prim, double t,
                         double *out_x, double *out_y, double *out_z);

#ifdef __cplusplus
}
#endif

#endif /* LISCIO_H */
