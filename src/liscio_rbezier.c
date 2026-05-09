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
* Description: liscio rational cubic Bezier utilities (NURBS).
*                (1) Exact arc → rational cubic Bezier conversion.
*                    Degree-elevate quadratic rational form (P0, Q,
*                    P2, weights (1, cos(θ/2), 1)) to cubic:
*                      P1' = (P0 + 2·wq·Q) / (1 + 2·wq)
*                      P2' = (2·wq·Q + P3) / (2·wq + 1)
*                      w_mid = (1 + 2·wq) / 3
*                    Q = tangent-intersection apex.
*                    Geometrically exact for arc_angle ≤ 2π/3.
*                (2) Rational cubic Bezier point evaluator:
*                      C(t) = Σ(w_i·B_i(t)·P_i) / Σ(w_i·B_i(t))
*              Reduces to non-rational Bezier when all w_i = 1.
* References:  L. Piegl & W. Tiller, *The NURBS Book*, 2nd ed.,
*              Springer 1997 — App. A.1 (arc to rational Bezier).
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

static inline double B0(double t) { double u = 1.0 - t; return u*u*u; }
static inline double B1(double t) { double u = 1.0 - t; return 3.0*t*u*u; }
static inline double B2(double t) { double u = 1.0 - t; return 3.0*t*t*u; }
static inline double B3(double t) { return t*t*t; }

/* Convert an ARC primitive to an equivalent rational cubic Bezier.
 * out must point to a primitive; on return out->type = LISCIO_PRIM_BEZIER
 * with weights set to (1, cos(θ/2), cos(θ/2), 1).
 *
 * Works exactly for |arc_angle| ≤ 2π/3 (~120°); for wider arcs, fit
 * is still continuous but no longer G² exact. Caller should split
 * wide arcs upstream for best precision.
 *
 * Returns 0 on success. */
int liscio_arc_to_rational_bezier(const liscio_primitive_t *arc,
                                   liscio_primitive_t *out)
{
    if (!arc || !out) return -1;
    if (arc->type != LISCIO_PRIM_ARC) return -1;

    double theta = arc->arc_angle;
    double half = theta * 0.5;
    double tan_half = tan(fabs(half));
    double cos_half = cos(fabs(half));
    double r = arc->radius;

    double tsx = arc->tan_start_x, tsy = arc->tan_start_y, tsz = arc->tan_start_z;
    double tex = arc->tan_end_x,   tey = arc->tan_end_y,   tez = arc->tan_end_z;
    double tsmag = sqrt(tsx*tsx + tsy*tsy + tsz*tsz);
    double temag = sqrt(tex*tex + tey*tey + tez*tez);
    if (tsmag < 1e-15 || temag < 1e-15) return -1;

    /* Quadratic rational Bezier is the canonical exact-arc form:
     *   P0, Q (tangent intersection = apex), P2  with weights (1, cos(θ/2), 1)
     * Q = P0 + (r·tan(θ/2))·T_start = P2 − (r·tan(θ/2))·T_end
     * Degree-elevate to cubic for unified primitive-type storage:
     *   P1' = (P0 + 2·w_q·Q) / (1 + 2·w_q)
     *   P2' = (2·w_q·Q + P3) / (2·w_q + 1)
     *   w1' = (1 + 2·w_q) / 3
     *   w2' = (2·w_q + 1) / 3
     *   w0' = w3' = 1 */
    double Qx = arc->start.x + r * tan_half * tsx / tsmag;
    double Qy = arc->start.y + r * tan_half * tsy / tsmag;
    double Qz = arc->start.z + r * tan_half * tsz / tsmag;

    double wq = cos_half;
    double w_mid = (1.0 + 2.0 * wq) / 3.0;

    memset(out, 0, sizeof(*out));
    out->type  = LISCIO_PRIM_BEZIER;
    out->start = arc->start;
    out->end   = arc->end;
    out->P0    = arc->start;
    out->P3    = arc->end;

    out->P1.x = (arc->start.x + 2.0 * wq * Qx) / (1.0 + 2.0 * wq);
    out->P1.y = (arc->start.y + 2.0 * wq * Qy) / (1.0 + 2.0 * wq);
    out->P1.z = (arc->start.z + 2.0 * wq * Qz) / (1.0 + 2.0 * wq);
    out->P2.x = (2.0 * wq * Qx + arc->end.x) / (2.0 * wq + 1.0);
    out->P2.y = (2.0 * wq * Qy + arc->end.y) / (2.0 * wq + 1.0);
    out->P2.z = (2.0 * wq * Qz + arc->end.z) / (2.0 * wq + 1.0);

    out->w0 = 1.0;
    out->w1 = w_mid;
    out->w2 = w_mid;
    out->w3 = 1.0;

    /* Rotary/uvw: linear (Bezier CP with 1/3-offset). */
    out->P0.a = arc->start.a; out->P3.a = arc->end.a;
    out->P1.a = arc->start.a + (arc->end.a - arc->start.a) / 3.0;
    out->P2.a = arc->start.a + 2.0 * (arc->end.a - arc->start.a) / 3.0;
    out->P0.b = arc->start.b; out->P3.b = arc->end.b;
    out->P1.b = arc->start.b + (arc->end.b - arc->start.b) / 3.0;
    out->P2.b = arc->start.b + 2.0 * (arc->end.b - arc->start.b) / 3.0;
    out->P0.c = arc->start.c; out->P3.c = arc->end.c;
    out->P1.c = arc->start.c + (arc->end.c - arc->start.c) / 3.0;
    out->P2.c = arc->start.c + 2.0 * (arc->end.c - arc->start.c) / 3.0;
    out->P0.u = arc->start.u; out->P3.u = arc->end.u;
    out->P1.u = arc->start.u + (arc->end.u - arc->start.u) / 3.0;
    out->P2.u = arc->start.u + 2.0 * (arc->end.u - arc->start.u) / 3.0;
    out->P0.v = arc->start.v; out->P3.v = arc->end.v;
    out->P1.v = arc->start.v + (arc->end.v - arc->start.v) / 3.0;
    out->P2.v = arc->start.v + 2.0 * (arc->end.v - arc->start.v) / 3.0;
    out->P0.w = arc->start.w; out->P3.w = arc->end.w;
    out->P1.w = arc->start.w + (arc->end.w - arc->start.w) / 3.0;
    out->P2.w = arc->start.w + 2.0 * (arc->end.w - arc->start.w) / 3.0;

    out->tan_start_x = arc->tan_start_x;
    out->tan_start_y = arc->tan_start_y;
    out->tan_start_z = arc->tan_start_z;
    out->tan_end_x   = arc->tan_end_x;
    out->tan_end_y   = arc->tan_end_y;
    out->tan_end_z   = arc->tan_end_z;

    out->feedrate       = arc->feedrate;
    out->tag_line_id    = arc->tag_line_id;
    out->tag_line_first = arc->tag_line_first;
    out->tag_line_last  = arc->tag_line_last;
    out->n_absorbed     = arc->n_absorbed;
    out->tol_xyz     = arc->tol_xyz;
    out->tol_abc     = arc->tol_abc;
    return 0;
}

/* Evaluate rational cubic Bezier at t ∈ [0,1]. Writes 3D point to out.
 * For non-rational (w_i all 1.0) this reduces to Bezier evaluation. */
void liscio_rbezier_eval(const liscio_primitive_t *prim, double t,
                          double *ox, double *oy, double *oz)
{
    if (!prim || !ox || !oy || !oz) return;
    double b0 = B0(t), b1 = B1(t), b2 = B2(t), b3 = B3(t);
    double den = prim->w0*b0 + prim->w1*b1 + prim->w2*b2 + prim->w3*b3;
    if (fabs(den) < 1e-15) { *ox = *oy = *oz = 0; return; }

    double num_x = prim->w0*b0*prim->P0.x + prim->w1*b1*prim->P1.x
                 + prim->w2*b2*prim->P2.x + prim->w3*b3*prim->P3.x;
    double num_y = prim->w0*b0*prim->P0.y + prim->w1*b1*prim->P1.y
                 + prim->w2*b2*prim->P2.y + prim->w3*b3*prim->P3.y;
    double num_z = prim->w0*b0*prim->P0.z + prim->w1*b1*prim->P1.z
                 + prim->w2*b2*prim->P2.z + prim->w3*b3*prim->P3.z;
    *ox = num_x / den;
    *oy = num_y / den;
    *oz = num_z / den;
}
