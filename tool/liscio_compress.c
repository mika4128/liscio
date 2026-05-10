/********************************************************************
 * liscio_compress — G-code → compressed G-code tool.
 *
 *   input.ngc  →  rs274 解析 (含单位/平面/O-word 展开)
 *               →  liscio (微段 → ARC / Bezier / Spline / Helix 拟合)
 *               →  primitive emitter (按类型回写 G-code)
 *               →  output.ngc
 *
 * 输出标准 G-code: G0/G1/G2/G3/G4/G17/G18/G19/M0/M2/M6 等, 跨控制器
 * 兼容.  BEZIER/SPLINE 因 G-code 无原生表示, 当前实现回退为 G1 polyline
 * (该段无压缩, 但保几何精度).  HELIX 多圈用 LinuxCNC P<turns> 形式;
 * --portable 模式下拆成多个单圈 G2/G3.
 *
 * Usage:
 *   liscio_compress [--tol_xyz=N] [--output=out.ngc] [--precision=N]
 *                    [--arc_merge=N] [--corner_mode=N] [--portable]
 *                    [--bezier_samples=N] <input.ngc>
 *
 * License: MIT
 ********************************************************************/

#define _USE_MATH_DEFINES
#define _POSIX_C_SOURCE 199309L
#include "liscio/liscio.h"
#include "ngc_parser.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ */
/* Modal state for output dedup.                                        */

typedef enum {
    LISCIO_TARGET_PORTABLE = 0,   /* Standard G-code only (Fanuc/Siemens/...) */
    LISCIO_TARGET_LINUXCNC = 1,   /* Allow G5 cubic Bezier (XY plane only) */
} liscio_emit_target_t;

typedef struct {
    FILE  *out;
    int    precision;        /* decimal digits, default 5 */
    int    portable;         /* 1 = no LinuxCNC P<turns>, split helices */
    int    bezier_samples;   /* G1 sub-segments per Bezier/Spline fallback */
    int    target;           /* see liscio_emit_target_t */

    /* Modal output state */
    int    plane;            /* 17 / 18 / 19 */
    int    motion_mode;      /* 0/1/2/3, -1 = unset */
    double last_feed;
    double last_x, last_y, last_z;
    double last_a, last_b, last_c;
    double last_u, last_v, last_w;
    int    have_last_pos;

    /* Stats */
    long   n_lines, n_arcs, n_helix, n_rapids, n_stops, n_g1_resampled;
    long   n_g5_bezier;     /* LinuxCNC G5 cubic Bezier emissions */
    long   bytes_out;
} emit_ctx_t;

/* Test whether a BEZIER primitive lies on one of LinuxCNC's three
 * cardinal planes (G17 XY / G18 XZ / G19 YZ).  G5 is emitted in the
 * currently-active plane, so any of the three is acceptable.  Returns:
 *   17 → XY-planar (constant Z), 18 → XZ-planar (constant Y),
 *   19 → YZ-planar (constant X), -1 → tilted/3D (must G1-resample).
 *
 * Also rejects rotary / additional-linear motion across the curve —
 * G5 has no syntax for those axes. */
static int bezier_planar_g5(const liscio_primitive_t *p)
{
    const double zt = 1e-6;
    const liscio_pose_t *pp[4] = { &p->P0, &p->P1, &p->P2, &p->P3 };

    /* No rotary / additional-linear motion across the curve. */
    for (int i = 1; i < 4; i++) {
        if (fabs(pp[i]->a - pp[0]->a) > zt) return -1;
        if (fabs(pp[i]->b - pp[0]->b) > zt) return -1;
        if (fabs(pp[i]->c - pp[0]->c) > zt) return -1;
        if (fabs(pp[i]->u - pp[0]->u) > zt) return -1;
        if (fabs(pp[i]->v - pp[0]->v) > zt) return -1;
        if (fabs(pp[i]->w - pp[0]->w) > zt) return -1;
    }

    /* Plane detection — which axis is invariant across all 4 control points? */
    int z_const = 1, y_const = 1, x_const = 1;
    double x = pp[0]->x, y = pp[0]->y, z = pp[0]->z;
    for (int i = 1; i < 4; i++) {
        if (fabs(pp[i]->z - z) > zt) z_const = 0;
        if (fabs(pp[i]->y - y) > zt) y_const = 0;
        if (fabs(pp[i]->x - x) > zt) x_const = 0;
    }
    if (z_const) return 17;   /* XY plane → G17 */
    if (y_const) return 18;   /* XZ plane → G18 */
    if (x_const) return 19;   /* YZ plane → G19 */
    return -1;                /* tilted */
}

static void fmt_axis(char *buf, size_t n, char letter, double val,
                     double last, int have_last, int prec)
{
    if (have_last && fabs(val - last) < 5e-7) {
        buf[0] = '\0';
        return;
    }
    snprintf(buf, n, " %c%.*g", letter, prec, val);
}

static void emit_pose_xyz(emit_ctx_t *e, double x, double y, double z)
{
    char bx[32], by[32], bz[32];
    fmt_axis(bx, sizeof(bx), 'X', x, e->last_x, e->have_last_pos, e->precision);
    fmt_axis(by, sizeof(by), 'Y', y, e->last_y, e->have_last_pos, e->precision);
    fmt_axis(bz, sizeof(bz), 'Z', z, e->last_z, e->have_last_pos, e->precision);
    fputs(bx, e->out);
    fputs(by, e->out);
    fputs(bz, e->out);
    e->last_x = x; e->last_y = y; e->last_z = z;
    e->have_last_pos = 1;
}

static void emit_rotary_uvw(emit_ctx_t *e, const liscio_pose_t *p)
{
    char b[10][32];
    int  k = 0;
    fmt_axis(b[k++], sizeof(b[0]), 'A', p->a, e->last_a, e->have_last_pos, e->precision);
    fmt_axis(b[k++], sizeof(b[0]), 'B', p->b, e->last_b, e->have_last_pos, e->precision);
    fmt_axis(b[k++], sizeof(b[0]), 'C', p->c, e->last_c, e->have_last_pos, e->precision);
    fmt_axis(b[k++], sizeof(b[0]), 'U', p->u, e->last_u, e->have_last_pos, e->precision);
    fmt_axis(b[k++], sizeof(b[0]), 'V', p->v, e->last_v, e->have_last_pos, e->precision);
    fmt_axis(b[k++], sizeof(b[0]), 'W', p->w, e->last_w, e->have_last_pos, e->precision);
    for (int i = 0; i < k; i++) fputs(b[i], e->out);
    e->last_a = p->a; e->last_b = p->b; e->last_c = p->c;
    e->last_u = p->u; e->last_v = p->v; e->last_w = p->w;
}

static void emit_feed(emit_ctx_t *e, double feed)
{
    if (feed <= 0.0) return;
    if (e->have_last_pos && fabs(feed - e->last_feed) < 1e-6) return;
    fprintf(e->out, " F%.*g", e->precision, feed);
    e->last_feed = feed;
}

static void emit_motion_mode(emit_ctx_t *e, int mode)
{
    if (e->motion_mode == mode) return;
    fprintf(e->out, "G%d", mode);
    e->motion_mode = mode;
}

static void emit_plane(emit_ctx_t *e, int plane)
{
    if (e->plane == plane) return;
    fprintf(e->out, "G%d\n", plane);
    e->plane = plane;
    e->motion_mode = -1;   /* plane change forces re-emit of motion mode */
}

/* Pick output plane for arc based on its normal.  Returns 17/18/19;
 * 0 means "tilted plane, can't use G2/G3 — caller must resample". */
static int arc_plane_from_normal(double nx, double ny, double nz)
{
    double anx = fabs(nx), any = fabs(ny), anz = fabs(nz);
    if (anz > 0.999 && anx + any < 1e-3) return 17;
    if (any > 0.999 && anx + anz < 1e-3) return 18;
    if (anx > 0.999 && any + anz < 1e-3) return 19;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Per-primitive emitters                                               */

static void emit_line(emit_ctx_t *e, const liscio_primitive_t *p)
{
    /* Initial state may not have last_pos: emit absolute G1 */
    fputs("\n", e->out);
    emit_motion_mode(e, 1);
    emit_pose_xyz(e, p->end.x, p->end.y, p->end.z);
    emit_rotary_uvw(e, &p->end);
    emit_feed(e, p->feedrate);
    e->n_lines++;
}

static void emit_arc(emit_ctx_t *e, const liscio_primitive_t *p,
                     int n_full_turns)
{
    int plane = arc_plane_from_normal(p->nx, p->ny, p->nz);
    if (plane == 0) {
        /* Tilted plane: caller should have resampled.  Emit warning. */
        fprintf(stderr, "[warn] tilted-plane arc, resampling to G1 polyline\n");
        /* Resample arc here into G1 segments — fall back. */
        const int N = 32;
        for (int i = 1; i <= N; i++) {
            double t = (double)i / N;
            double ang = t * p->arc_angle;
            double ux = p->start.x - p->cx;
            double uy = p->start.y - p->cy;
            double uz = p->start.z - p->cz;
            double um = sqrt(ux*ux + uy*uy + uz*uz);
            if (um < 1e-15) break;
            ux /= um; uy /= um; uz /= um;
            double vx = p->ny*uz - p->nz*uy;
            double vy = p->nz*ux - p->nx*uz;
            double vz = p->nx*uy - p->ny*ux;
            double off = t * ((p->end.x - p->start.x)*p->nx
                            + (p->end.y - p->start.y)*p->ny
                            + (p->end.z - p->start.z)*p->nz);
            double x = p->cx + p->radius*(cos(ang)*ux + sin(ang)*vx) + off*p->nx;
            double y = p->cy + p->radius*(cos(ang)*uy + sin(ang)*vy) + off*p->ny;
            double z = p->cz + p->radius*(cos(ang)*uz + sin(ang)*vz) + off*p->nz;
            fputs("\n", e->out);
            emit_motion_mode(e, 1);
            emit_pose_xyz(e, x, y, z);
            emit_feed(e, p->feedrate);
            e->n_g1_resampled++;
        }
        return;
    }

    emit_plane(e, plane);
    /* Sign of arc_angle + sign of normal-dot-plane-axis decides G2/G3.
     * Convention: arc_angle > 0 = CCW about +n; XY plane n=+Z → G3; n=-Z → G2.
     * Simpler: rotation direction in 2D plane = sign(arc_angle * n_dot_plane). */
    double n_axis = (plane == 17) ? p->nz :
                    (plane == 18) ? p->ny : p->nx;
    int g_code = ((p->arc_angle * n_axis) > 0) ? 3 : 2;

    fputs("\n", e->out);
    emit_motion_mode(e, g_code);
    emit_pose_xyz(e, p->end.x, p->end.y, p->end.z);
    /* Plane-relevant IJK only (drop the orthogonal-to-plane offset which
     * is structurally zero).  Saves bytes vs always-3-axis form. */
    double di = p->cx - p->start.x;
    double dj = p->cy - p->start.y;
    double dk = p->cz - p->start.z;
    if (plane == 17) {
        fprintf(e->out, " I%.*g J%.*g", e->precision, di, e->precision, dj);
    } else if (plane == 18) {
        fprintf(e->out, " I%.*g K%.*g", e->precision, di, e->precision, dk);
    } else /* 19 */ {
        fprintf(e->out, " J%.*g K%.*g", e->precision, dj, e->precision, dk);
    }
    if (n_full_turns > 0)
        fprintf(e->out, " P%d", n_full_turns);
    emit_feed(e, p->feedrate);
    if (n_full_turns > 0) e->n_helix++;
    else                  e->n_arcs++;
}

static void emit_helix(emit_ctx_t *e, const liscio_primitive_t *p)
{
    /* HELIX: |arc_angle| > 2π. */
    double turns_total = fabs(p->arc_angle) / (2.0 * M_PI);
    int    extra = (int)floor(turns_total - 1.0 + 1e-9);  /* P<n> = additional full turns */
    if (extra < 0) extra = 0;

    int plane = arc_plane_from_normal(p->nx, p->ny, p->nz);
    if (e->portable || plane == 0) {
        /* Split into single-loop G2/G3 segments.  Each one rotates 2π
         * (or less for the final) in the same direction. */
        int n_segs = (int)ceil(fabs(p->arc_angle) / (2.0 * M_PI));
        if (n_segs < 1) n_segs = 1;
        liscio_primitive_t seg = *p;
        for (int i = 0; i < n_segs; i++) {
            double t0 = (double)i     / n_segs;
            double t1 = (double)(i+1) / n_segs;
            double a0 = t0 * p->arc_angle;
            double a1 = t1 * p->arc_angle;
            double ux = p->start.x - p->cx;
            double uy = p->start.y - p->cy;
            double uz = p->start.z - p->cz;
            double um = sqrt(ux*ux + uy*uy + uz*uz);
            if (um < 1e-15) break;
            ux /= um; uy /= um; uz /= um;
            double vx = p->ny*uz - p->nz*uy;
            double vy = p->nz*ux - p->nx*uz;
            double vz = p->nx*uy - p->ny*ux;
            double off_total = (p->end.x - p->start.x)*p->nx
                             + (p->end.y - p->start.y)*p->ny
                             + (p->end.z - p->start.z)*p->nz;
            double off0 = t0 * off_total;
            double off1 = t1 * off_total;
            seg.start.x = p->cx + p->radius*(cos(a0)*ux + sin(a0)*vx) + off0*p->nx;
            seg.start.y = p->cy + p->radius*(cos(a0)*uy + sin(a0)*vy) + off0*p->ny;
            seg.start.z = p->cz + p->radius*(cos(a0)*uz + sin(a0)*vz) + off0*p->nz;
            seg.end.x   = p->cx + p->radius*(cos(a1)*ux + sin(a1)*vx) + off1*p->nx;
            seg.end.y   = p->cy + p->radius*(cos(a1)*uy + sin(a1)*vy) + off1*p->ny;
            seg.end.z   = p->cz + p->radius*(cos(a1)*uz + sin(a1)*vz) + off1*p->nz;
            seg.arc_angle = a1 - a0;
            emit_arc(e, &seg, 0);
        }
    } else {
        /* LinuxCNC compact form: single G2/G3 with P<extra_turns>. */
        emit_arc(e, p, extra);
    }
}

static void emit_bezier(emit_ctx_t *e, const liscio_primitive_t *p)
{
    /* LinuxCNC G5: cubic Bezier in any of XY/XZ/YZ planes (active plane
     * decided by G17/G18/G19 mode).  Available only when target=LINUXCNC
     * AND the curve is planar in one cardinal plane AND no ABCUVW motion.
     *
     * In each plane the X/Y G5 letters map to the in-plane axes:
     *   G17 (XY): X=X Y=Y, I=ΔX_first J=ΔY_first, P=ΔX_last Q=ΔY_last
     *   G18 (XZ): X=X Z=Z, I=ΔX_first K=ΔZ_first, P=ΔX_last R=ΔZ_last  *
     *   G19 (YZ): Y=Y Z=Z, J=ΔY_first K=ΔZ_first, Q=ΔY_last R=ΔZ_last  *
     *
     * (* LinuxCNC source still names letters I/J/P/Q regardless of plane;
     * the active plane chooses which 2D axes the letters refer to.) */
    int g5_plane = (e->target == LISCIO_TARGET_LINUXCNC) ? bezier_planar_g5(p) : -1;
    if (g5_plane > 0) {
        /* LinuxCNC's interp_convert.cc convert_spline() reads I/J/P/Q
         * with plane-dependent letter semantics:
         *   XY: I=ΔX_first  J=ΔY_first  P=ΔX_last  Q=ΔY_last
         *   YZ: I=ΔY_first  J=ΔZ_first  P=ΔY_last  Q=ΔZ_last
         *   XZ: I=ΔX_first  J=ΔZ_first  P=ΔZ_last  Q=ΔX_last  ← P/Q reversed!
         * Endpoint coordinates are written with their actual axis letters
         * (X/Y, X/Z, Y/Z) regardless of plane. */
        double a3=0, b3=0;                 /* in-plane endpoint coords   */
        double da1=0, db1=0, da2=0, db2=0; /* (P1-P0), (P2-P3) in-plane  */
        char la = 'X', lb = 'Y';
        int xz_swap = 0;
        switch (g5_plane) {
        case 17:
            la='X'; lb='Y';
            a3=p->P3.x; b3=p->P3.y;
            da1=p->P1.x-p->P0.x; db1=p->P1.y-p->P0.y;
            da2=p->P2.x-p->P3.x; db2=p->P2.y-p->P3.y;
            break;
        case 18:
            la='X'; lb='Z';
            a3=p->P3.x; b3=p->P3.z;
            da1=p->P1.x-p->P0.x; db1=p->P1.z-p->P0.z;
            da2=p->P2.x-p->P3.x; db2=p->P2.z-p->P3.z;
            xz_swap = 1;       /* P/Q letters swap for XZ plane only */
            break;
        case 19:
            la='Y'; lb='Z';
            a3=p->P3.y; b3=p->P3.z;
            da1=p->P1.y-p->P0.y; db1=p->P1.z-p->P0.z;
            da2=p->P2.y-p->P3.y; db2=p->P2.z-p->P3.z;
            break;
        }
        emit_plane(e, g5_plane);
        fputs("\n", e->out);
        e->motion_mode = -1;   /* G5 is non-modal — force re-emit next */
        if (xz_swap) {
            /* XZ: P=ΔZ_last, Q=ΔX_last */
            fprintf(e->out,
                "G5 %c%.*g %c%.*g I%.*g J%.*g P%.*g Q%.*g",
                la, e->precision, a3,
                lb, e->precision, b3,
                e->precision, da1, e->precision, db1,
                e->precision, db2, e->precision, da2);
        } else {
            fprintf(e->out,
                "G5 %c%.*g %c%.*g I%.*g J%.*g P%.*g Q%.*g",
                la, e->precision, a3,
                lb, e->precision, b3,
                e->precision, da1, e->precision, db1,
                e->precision, da2, e->precision, db2);
        }
        emit_feed(e, p->feedrate);
        e->last_x = p->P3.x; e->last_y = p->P3.y; e->last_z = p->P3.z;
        e->have_last_pos = 1;
        e->n_g5_bezier++;
        return;
    }

    /* Fallback: resample to G1 polyline (portable G-code). */
    int N = e->bezier_samples > 0 ? e->bezier_samples : 16;
    for (int i = 1; i <= N; i++) {
        double t = (double)i / N;
        double u = 1.0 - t;
        double b0 = u*u*u, b1 = 3*t*u*u, b2 = 3*t*t*u, b3 = t*t*t;
        double x = b0*p->P0.x + b1*p->P1.x + b2*p->P2.x + b3*p->P3.x;
        double y = b0*p->P0.y + b1*p->P1.y + b2*p->P2.y + b3*p->P3.y;
        double z = b0*p->P0.z + b1*p->P1.z + b2*p->P2.z + b3*p->P3.z;
        fputs("\n", e->out);
        emit_motion_mode(e, 1);
        emit_pose_xyz(e, x, y, z);
        emit_feed(e, p->feedrate);
        e->n_g1_resampled++;
    }
}

static void emit_spline(emit_ctx_t *e, const liscio_primitive_t *p)
{
    /* Resample via liscio_bspline_eval (declared in liscio.h). */
    int N = e->bezier_samples > 0 ? e->bezier_samples : 16;
    for (int i = 1; i <= N; i++) {
        double t = (double)i / N;
        double x, y, z;
        liscio_bspline_eval(p, t, &x, &y, &z);
        fputs("\n", e->out);
        emit_motion_mode(e, 1);
        emit_pose_xyz(e, x, y, z);
        emit_feed(e, p->feedrate);
        e->n_g1_resampled++;
    }
}

static void emit_rapid(emit_ctx_t *e, const liscio_primitive_t *p)
{
    fputs("\n", e->out);
    emit_motion_mode(e, 0);
    emit_pose_xyz(e, p->end.x, p->end.y, p->end.z);
    emit_rotary_uvw(e, &p->end);
    e->n_rapids++;
}

static void emit_stop(emit_ctx_t *e, const liscio_primitive_t *p)
{
    fputs("\n", e->out);
    e->motion_mode = -1;   /* M-code/G4 break modal motion */
    switch (p->stop_reason) {
    case LISCIO_STOP_DWELL:
        fprintf(e->out, "G4 P%.*g", e->precision, p->stop_dwell_seconds);
        break;
    case LISCIO_STOP_PROGRAM:           fputs("M0", e->out); break;
    case LISCIO_STOP_OPTIONAL:          fputs("M1", e->out); break;
    case LISCIO_STOP_PROGRAM_END:       fputs("M2", e->out); break;
    case LISCIO_STOP_PROGRAM_END_REWIND:fputs("M30", e->out); break;
    case LISCIO_STOP_TOOL_CHANGE:       fputs("M6", e->out); break;
    case LISCIO_STOP_SET_OFFSETS:
        fprintf(e->out, "( G92/G10 set-offsets event from line %d )", p->tag_line_id);
        break;
    case LISCIO_STOP_TOOL_LENGTH:
        fprintf(e->out, "( G43/G49 tool-length event from line %d )", p->tag_line_id);
        break;
    case LISCIO_STOP_COORD_SYSTEM:
        fprintf(e->out, "( coord-system change from line %d )", p->tag_line_id);
        break;
    default:
        fprintf(e->out, "( STOP reason=%d line=%d )", p->stop_reason, p->tag_line_id);
        break;
    }
    e->n_stops++;
}

/* ------------------------------------------------------------------ */
/* Liscio primitive callback                                            */

static void on_prim(const liscio_primitive_t *p, void *user)
{
    emit_ctx_t *e = (emit_ctx_t *)user;
    switch (p->type) {
    case LISCIO_PRIM_LINE:    emit_line(e, p);    break;
    case LISCIO_PRIM_ARC:     emit_arc(e, p, 0);  break;
    case LISCIO_PRIM_HELIX:   emit_helix(e, p);   break;
    case LISCIO_PRIM_BEZIER:  emit_bezier(e, p);  break;
    case LISCIO_PRIM_SPLINE:  emit_spline(e, p);  break;
    case LISCIO_PRIM_RAPID:   emit_rapid(e, p);   break;
    case LISCIO_PRIM_STOP:    emit_stop(e, p);    break;
    }
}

/* ------------------------------------------------------------------ */
/* NGC parser bridge                                                    */

typedef struct {
    liscio_ctx_t *ctx;
    double        feed;
} bridge_ctx_t;

static void on_move(const ngc_move_t *m, void *user)
{
    bridge_ctx_t *b = (bridge_ctx_t *)user;
    if (m->is_rapid) {
        liscio_pose_t s = {m->start.x,m->start.y,m->start.z,
                           m->start.a,m->start.b,m->start.c,
                           m->start.u,m->start.v,m->start.w};
        liscio_pose_t e = {m->end.x,m->end.y,m->end.z,
                           m->end.a,m->end.b,m->end.c,
                           m->end.u,m->end.v,m->end.w};
        liscio_add_rapid(b->ctx, &s, &e, m->line_no);
        return;
    }
    if (m->feedrate > 0.0) b->feed = m->feedrate;
    if (m->type == NGC_MOVE_LINE) {
        liscio_pose_t s = {m->start.x,m->start.y,m->start.z,
                           m->start.a,m->start.b,m->start.c,
                           m->start.u,m->start.v,m->start.w};
        liscio_pose_t e = {m->end.x,m->end.y,m->end.z,
                           m->end.a,m->end.b,m->end.c,
                           m->end.u,m->end.v,m->end.w};
        liscio_add_line(b->ctx, &s, &e, b->feed, m->line_no);
    } else {
        liscio_primitive_t a; memset(&a, 0, sizeof(a));
        a.type     = LISCIO_PRIM_ARC;
        a.start    = (liscio_pose_t){m->start.x,m->start.y,m->start.z,
                                     m->start.a,m->start.b,m->start.c,
                                     m->start.u,m->start.v,m->start.w};
        a.end      = (liscio_pose_t){m->end.x,m->end.y,m->end.z,
                                     m->end.a,m->end.b,m->end.c,
                                     m->end.u,m->end.v,m->end.w};
        a.cx = m->cx; a.cy = m->cy; a.cz = m->cz;
        a.nx = m->nx; a.ny = m->ny; a.nz = m->nz;
        a.radius = m->radius;
        a.arc_angle = m->arc_angle;
        liscio_add_arc(b->ctx, &a, b->feed, m->line_no);
    }
}

/* ------------------------------------------------------------------ */
/* Main                                                                  */

int main(int argc, char **argv)
{
    liscio_cfg_t cfg; liscio_cfg_default(&cfg);
    cfg.tol_xyz = 0.01;        /* default 10µm — industrial sweet spot */
    const char *out_path = NULL;
    int precision = 5;
    int portable  = 0;
    /* Default 8: 18-file round-trip stays within tol_xyz (verified).
     * 16 is safer for sub-µm work but inflates BEZIER-dense 3D files
     * (flower.ngc 0.59× → 1.12× compression at N=8 vs N=16). */
    int bezier_samples = 8;
    int i0 = 1;
    for (int i = 1; i < argc; i++) {
        if (strncmp(argv[i], "--tol_xyz=", 10) == 0) {
            cfg.tol_xyz = atof(argv[i] + 10); i0 = i + 1;
        } else if (strncmp(argv[i], "--output=", 9) == 0) {
            out_path = argv[i] + 9; i0 = i + 1;
        } else if (strncmp(argv[i], "--precision=", 12) == 0) {
            precision = atoi(argv[i] + 12); i0 = i + 1;
        } else if (strncmp(argv[i], "--bezier_samples=", 17) == 0) {
            bezier_samples = atoi(argv[i] + 17); i0 = i + 1;
        } else if (strncmp(argv[i], "--arc_merge=", 12) == 0) {
            int m = atoi(argv[i] + 12);
            liscio_cfg_set_arc_merge(&cfg,
                m == 2 ? LISCIO_ARC_MERGE_HELIX :
                m == 1 ? LISCIO_ARC_MERGE_ARC : LISCIO_ARC_MERGE_OFF, 0, 0);
            i0 = i + 1;
        } else if (strncmp(argv[i], "--corner_mode=", 14) == 0) {
            int m = atoi(argv[i] + 14);
            liscio_cfg_set_corner_detection(&cfg,
                m ? LISCIO_CORNER_LOOKAHEAD : LISCIO_CORNER_IMMEDIATE, 0, 0);
            i0 = i + 1;
        } else if (strcmp(argv[i], "--portable") == 0) {
            portable = 1; i0 = i + 1;
        } else if (strncmp(argv[i], "--target=", 9) == 0) {
            const char *t = argv[i] + 9;
            if      (!strcmp(t, "linuxcnc")) /* enables G5 */;
            else if (!strcmp(t, "portable")) /* default */;
            else { fprintf(stderr, "unknown --target=%s\n", t); return 1; }
            /* set later once emit_ctx exists; for now stash */
            i0 = i + 1;
        } else break;
    }
    int target = LISCIO_TARGET_PORTABLE;
    for (int i = 1; i < i0; i++) {
        if (strncmp(argv[i], "--target=", 9) == 0 &&
            !strcmp(argv[i]+9, "linuxcnc"))
            target = LISCIO_TARGET_LINUXCNC;
    }
    if (i0 >= argc) {
        fprintf(stderr,
            "usage: %s [--tol_xyz=N] [--output=out.ngc] [--precision=N]\n"
            "          [--arc_merge=0|1|2] [--corner_mode=0|1] [--portable]\n"
            "          [--bezier_samples=N] <input.ngc>\n",
            argv[0]);
        return 1;
    }
    /* arc_merge defaults to HELIX so multi-loop spirals get compact P-form */
    if (cfg.arc_merge_mode == 0)
        liscio_cfg_set_arc_merge(&cfg, LISCIO_ARC_MERGE_HELIX, 0, 0);

    FILE *out = stdout;
    if (out_path) {
        out = fopen(out_path, "w");
        if (!out) { perror(out_path); return 1; }
    }

    emit_ctx_t e = {0};
    e.out = out;
    e.precision = precision;
    e.portable  = portable;
    e.target    = target;
    e.bezier_samples = bezier_samples;
    e.plane = 0;          /* unset → first arc emits its plane */
    e.motion_mode = -1;
    e.last_feed = -1.0;

    /* Preamble */
    const char *in_name = argv[i0];
    const char *base = strrchr(in_name, '/');
    base = base ? base + 1 : in_name;
    fprintf(out, "%%\n");
    fprintf(out, "( liscio-compressed: tol=%.*g mm; in=%s )\n",
            precision, cfg.tol_xyz, base);
    fprintf(out, "G21\n");                      /* mm */
    fprintf(out, "G17\n");  e.plane = 17;       /* default XY */
    fprintf(out, "G90\n");                      /* absolute */
    fprintf(out, "G64 P%.*g\n", precision, cfg.tol_xyz); /* path tolerance */

    liscio_ctx_t *ctx = liscio_create(&cfg);
    liscio_set_callback(ctx, on_prim, &e);
    bridge_ctx_t br = { ctx, 0.0 };
    ngc_stats_t st;
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    int rc = ngc_parse_file(in_name, on_move, &br, &st);
    if (rc != 0) {
        fprintf(stderr, "ngc_parse_file %s failed\n", in_name);
        liscio_destroy(ctx);
        if (out_path) fclose(out);
        return 1;
    }
    liscio_flush(ctx);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed_ms = (t1.tv_sec - t0.tv_sec)*1000.0
                      + (t1.tv_nsec - t0.tv_nsec)*1e-6;
    liscio_destroy(ctx);

    fprintf(out, "\nM2\n%%\n");

    long out_emitted = e.n_lines + e.n_arcs + e.n_helix + e.n_g5_bezier
                     + e.n_g1_resampled;
    long input_motions = st.n_g1 + st.n_g2g3;
    fprintf(stderr,
        "%-30s  in_g1=%ld  in_g2g3=%ld  in_g0=%ld  target=%s\n"
        "  out: line=%ld arc=%ld helix=%ld G5=%ld rapid=%ld stop=%ld resampled=%ld\n"
        "  compression=%.2fx (motion blocks); elapsed=%.2f ms\n",
        base, st.n_g1, st.n_g2g3, st.n_g0,
        target == LISCIO_TARGET_LINUXCNC ? "linuxcnc" : "portable",
        e.n_lines, e.n_arcs, e.n_helix, e.n_g5_bezier,
        e.n_rapids, e.n_stops, e.n_g1_resampled,
        out_emitted ? (double)input_motions / out_emitted : 0,
        elapsed_ms);

    if (out_path) fclose(out);

    /* Inflation check.  Standard G-code has no native cubic Bezier, so
     * BEZIER-dense inputs may inflate when resampled to G1.  When the
     * compressed file is bigger than the original, warn the user and
     * suggest mitigation. */
    if (out_path) {
        struct stat in_st, out_st;
        if (stat(in_name, &in_st) == 0 && stat(out_path, &out_st) == 0
            && out_st.st_size > in_st.st_size)
        {
            double bloat = (double)out_st.st_size / (double)in_st.st_size;
            fprintf(stderr,
                "warning: compressed output %ld B > input %ld B (%.2f×).\n"
                "  This file resampled %ld BEZIER segments → %ld G1 lines.\n"
                "  Try a smaller --bezier_samples (current default 8); or, on\n"
                "  LinuxCNC with 2D BEZIER, --target=linuxcnc to emit native G5.\n"
                "  3D free-form data (BEZIER + Z variation) cannot be expressed\n"
                "  compactly in standard G-code.\n",
                (long)out_st.st_size, (long)in_st.st_size, bloat,
                e.n_g1_resampled / 8 /*approx*/, e.n_g1_resampled);
        }
    }
    return 0;
}
