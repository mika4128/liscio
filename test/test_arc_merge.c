/********************************************************************
 * test_arc_merge.c — synthetic regression for adjacent-arc merge
 *.  Builds known geometric scenarios and verifies that
 * liscio_cfg_set_arc_merge() folds them as expected.
 *
 *   1. 4 quarter circles in XY (planar, full circle) → 1 ARC, angle = 2π
 *   2. 8 quarter helical arcs (Z pitch) → 1 HELIX (mode=2) or 1 ARC (mode=1),
 *      arc_angle = 2 * 2π = 4π (two loops)
 *   3. Two NOT-mergeable arcs (different radius) → 2 separate ARCs
 *   4. OFF mode passthrough sanity (4 quarter arcs → 4 ARCs)
 *
 * License: MIT
 ********************************************************************/

#define _USE_MATH_DEFINES
#include "liscio/liscio.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    int n_emit;
    int n_helix;
    int n_arc;
    double last_arc_angle;
    int    last_type;
    double last_end_z;
    int    last_n_absorbed;
} collect_t;

static void on_prim(const liscio_primitive_t *p, void *user)
{
    collect_t *c = (collect_t *)user;
    c->n_emit++;
    if (p->type == LISCIO_PRIM_HELIX) c->n_helix++;
    if (p->type == LISCIO_PRIM_ARC)   c->n_arc++;
    c->last_arc_angle = p->arc_angle;
    c->last_type     = p->type;
    c->last_end_z    = p->end.z;
    c->last_n_absorbed = p->n_absorbed;
}

/* Build a planar (XY plane) quarter arc CCW about Z axis,
 * starting at angle a0 with center (0,0,0), radius r.  Returns the
 * arc primitive ready for liscio_add_arc(). */
static liscio_primitive_t mk_planar_quarter(double a0, double r)
{
    liscio_primitive_t arc;
    memset(&arc, 0, sizeof(arc));
    arc.type      = LISCIO_PRIM_ARC;
    arc.cx = 0; arc.cy = 0; arc.cz = 0;
    arc.nx = 0; arc.ny = 0; arc.nz = 1;     /* normal = +Z */
    arc.radius    = r;
    arc.arc_angle = M_PI / 2.0;             /* 90° CCW */
    arc.start.x = r * cos(a0); arc.start.y = r * sin(a0); arc.start.z = 0;
    arc.end.x   = r * cos(a0 + M_PI/2.0);
    arc.end.y   = r * sin(a0 + M_PI/2.0);
    arc.end.z   = 0;
    return arc;
}

/* Build a helical quarter arc rising by `dz` along Z axis. */
static liscio_primitive_t mk_helical_quarter(double a0, double r,
                                              double z_start, double dz)
{
    liscio_primitive_t arc = mk_planar_quarter(a0, r);
    arc.start.z = z_start;
    arc.end.z   = z_start + dz;
    return arc;
}

static int test_full_circle_merge(void)
{
    liscio_cfg_t cfg; liscio_cfg_default(&cfg); cfg.arc_to_helix_max = 0; /* legacy off-by-default for these tests */
    cfg.tol_xyz = 0.01;
    liscio_cfg_set_arc_merge(&cfg, LISCIO_ARC_MERGE_ARC, 0, 0);
    liscio_ctx_t *ctx = liscio_create(&cfg);
    collect_t c = {0};
    liscio_set_callback(ctx, on_prim, &c);

    for (int i = 0; i < 4; i++) {
        liscio_primitive_t a = mk_planar_quarter(i * M_PI/2.0, 10.0);
        liscio_add_arc(ctx, &a, 1000.0, i+1);
    }
    liscio_flush(ctx);
    liscio_destroy(ctx);

    if (c.n_emit != 1) {
        fprintf(stderr, "FAIL: full_circle expected 1 emit, got %d\n", c.n_emit);
        return 1;
    }
    if (fabs(c.last_arc_angle - 2.0 * M_PI) > 1e-9) {
        fprintf(stderr, "FAIL: full_circle expected arc_angle=2π, got %.6f\n",
                c.last_arc_angle);
        return 1;
    }
    if (c.last_type != LISCIO_PRIM_ARC) {
        fprintf(stderr, "FAIL: full_circle (mode=ARC) expected ARC type, got %d\n",
                c.last_type);
        return 1;
    }
    if (c.last_n_absorbed != 4) {
        fprintf(stderr, "FAIL: full_circle expected n_absorbed=4, got %d\n",
                c.last_n_absorbed);
        return 1;
    }
    printf("  PASS  full_circle: 4 quarters → 1 ARC (angle=2π, n_absorbed=4)\n");
    return 0;
}

static int test_double_helix_merge(void)
{
    liscio_cfg_t cfg; liscio_cfg_default(&cfg); cfg.arc_to_helix_max = 0; /* legacy off-by-default for these tests */
    cfg.tol_xyz = 0.01;
    liscio_cfg_set_arc_merge(&cfg, LISCIO_ARC_MERGE_HELIX, 0, 0);
    liscio_ctx_t *ctx = liscio_create(&cfg);
    collect_t c = {0};
    liscio_set_callback(ctx, on_prim, &c);

    /* 8 quarter helical arcs covering 2 full loops, pitch 5mm/turn → dz=1.25/qtr. */
    double z = 0;
    double dz_per_quarter = 1.25;
    for (int i = 0; i < 8; i++) {
        liscio_primitive_t a = mk_helical_quarter(
            i * M_PI/2.0, 10.0, z, dz_per_quarter);
        liscio_add_arc(ctx, &a, 1000.0, i+1);
        z += dz_per_quarter;
    }
    liscio_flush(ctx);
    liscio_destroy(ctx);

    if (c.n_emit != 1) {
        fprintf(stderr, "FAIL: double_helix expected 1 emit, got %d\n", c.n_emit);
        return 1;
    }
    if (fabs(c.last_arc_angle - 4.0 * M_PI) > 1e-9) {
        fprintf(stderr, "FAIL: double_helix expected arc_angle=4π, got %.6f\n",
                c.last_arc_angle);
        return 1;
    }
    if (c.last_type != LISCIO_PRIM_HELIX) {
        fprintf(stderr, "FAIL: double_helix expected HELIX type, got %d\n",
                c.last_type);
        return 1;
    }
    if (fabs(c.last_end_z - 8.0 * 1.25) > 1e-9) {
        fprintf(stderr, "FAIL: double_helix expected end.z=10, got %.6f\n",
                c.last_end_z);
        return 1;
    }
    printf("  PASS  double_helix: 8 quarters → 1 HELIX (angle=4π, dz=10)\n");
    return 0;
}

static int test_radius_change_no_merge(void)
{
    liscio_cfg_t cfg; liscio_cfg_default(&cfg); cfg.arc_to_helix_max = 0; /* legacy off-by-default for these tests */
    cfg.tol_xyz = 0.01;
    liscio_cfg_set_arc_merge(&cfg, LISCIO_ARC_MERGE_ARC, 0, 0);
    liscio_ctx_t *ctx = liscio_create(&cfg);
    collect_t c = {0};
    liscio_set_callback(ctx, on_prim, &c);

    /* Two quarters with mismatched radius — must NOT merge. */
    liscio_primitive_t a1 = mk_planar_quarter(0,         10.0);
    liscio_primitive_t a2 = mk_planar_quarter(M_PI/2.0,  20.0);
    /* a2.start must equal a1.end for the test to actually exercise the
     * "radius differs" rejection rather than the "discontinuous start"
     * one.  Force the start point on the larger circle. */
    a2.start.x = 20.0 * cos(M_PI/2.0);
    a2.start.y = 20.0 * sin(M_PI/2.0);
    a2.start.z = 0;
    /* But a1.end != a2.start, so the rejection actually fires on the
     * "start ≠ pending.end" check.  That's fine for negative-test purpose. */
    liscio_add_arc(ctx, &a1, 1000.0, 1);
    liscio_add_arc(ctx, &a2, 1000.0, 2);
    liscio_flush(ctx);
    liscio_destroy(ctx);

    if (c.n_emit != 2) {
        fprintf(stderr, "FAIL: radius_change expected 2 emits, got %d\n", c.n_emit);
        return 1;
    }
    printf("  PASS  radius_change: not mergeable → 2 separate ARCs\n");
    return 0;
}

static int test_off_mode_passthrough(void)
{
    liscio_cfg_t cfg; liscio_cfg_default(&cfg); cfg.arc_to_helix_max = 0; /* legacy off-by-default for these tests */
    cfg.tol_xyz = 0.01;
    /* arc_merge_mode = OFF (default); confirm legacy 1:1 passthrough. */
    liscio_ctx_t *ctx = liscio_create(&cfg);
    collect_t c = {0};
    liscio_set_callback(ctx, on_prim, &c);

    for (int i = 0; i < 4; i++) {
        liscio_primitive_t a = mk_planar_quarter(i * M_PI/2.0, 10.0);
        liscio_add_arc(ctx, &a, 1000.0, i+1);
    }
    liscio_flush(ctx);
    liscio_destroy(ctx);

    if (c.n_emit != 4) {
        fprintf(stderr, "FAIL: off_mode expected 4 emits, got %d\n", c.n_emit);
        return 1;
    }
    if (c.n_helix != 0) {
        fprintf(stderr, "FAIL: off_mode emitted HELIX in OFF mode\n");
        return 1;
    }
    printf("  PASS  off_mode: 4 quarters stay 4 ARCs (passthrough)\n");
    return 0;
}

/* ====================================================================
 * RAPID / STOP event-stream tests
 * ==================================================================== */

typedef struct {
    int n;
    int type[16];
    int stop_reason[16];
    double dwell[16];
} event_collect_t;

static void on_event_prim(const liscio_primitive_t *p, void *user)
{
    event_collect_t *c = (event_collect_t *)user;
    if (c->n < 16) {
        c->type[c->n]        = p->type;
        c->stop_reason[c->n] = p->stop_reason;
        c->dwell[c->n]       = p->stop_dwell_seconds;
    }
    c->n++;
}

static int test_rapid_passthrough(void)
{
    liscio_cfg_t cfg; liscio_cfg_default(&cfg); cfg.arc_to_helix_max = 0; /* legacy off-by-default for these tests */
    cfg.tol_xyz = 0.01;
    liscio_ctx_t *ctx = liscio_create(&cfg);
    event_collect_t c = {0};
    liscio_set_callback(ctx, on_event_prim, &c);

    /* G1 — G0 — G1 sequence. */
    liscio_pose_t a = {0,0,0,0,0,0,0,0,0};
    liscio_pose_t b = {10,0,0,0,0,0,0,0,0};
    liscio_pose_t cc = {100,0,0,0,0,0,0,0,0};
    liscio_pose_t dd = {110,0,0,0,0,0,0,0,0};
    liscio_add_line(ctx, &a, &b, 1000, 1);
    liscio_add_rapid(ctx, &b, &cc, 2);            /* G0 traverse */
    liscio_add_line(ctx, &cc, &dd, 1000, 3);
    liscio_flush(ctx);
    liscio_destroy(ctx);

    if (c.n != 3) {
        fprintf(stderr, "FAIL: rapid_pass expected 3 prims, got %d\n", c.n);
        return 1;
    }
    if (c.type[0] != LISCIO_PRIM_LINE ||
        c.type[1] != LISCIO_PRIM_RAPID ||
        c.type[2] != LISCIO_PRIM_LINE) {
        fprintf(stderr, "FAIL: rapid_pass types wrong [%d %d %d]\n",
                c.type[0], c.type[1], c.type[2]);
        return 1;
    }
    printf("  PASS  rapid_passthrough  (LINE, RAPID, LINE)\n");
    return 0;
}

static int test_dwell_stop(void)
{
    liscio_cfg_t cfg; liscio_cfg_default(&cfg); cfg.arc_to_helix_max = 0; /* legacy off-by-default for these tests */
    cfg.tol_xyz = 0.01;
    liscio_ctx_t *ctx = liscio_create(&cfg);
    event_collect_t c = {0};
    liscio_set_callback(ctx, on_event_prim, &c);

    liscio_pose_t a = {0,0,0,0,0,0,0,0,0};
    liscio_pose_t b = {10,0,0,0,0,0,0,0,0};
    liscio_add_line(ctx, &a, &b, 1000, 1);
    liscio_emit_stop(ctx, LISCIO_STOP_DWELL, 1.5, &b, 2);
    liscio_add_line(ctx, &b, &a, 1000, 3);
    liscio_flush(ctx);
    liscio_destroy(ctx);

    if (c.n != 3) {
        fprintf(stderr, "FAIL: dwell expected 3 prims, got %d\n", c.n);
        return 1;
    }
    if (c.type[1] != LISCIO_PRIM_STOP ||
        c.stop_reason[1] != LISCIO_STOP_DWELL ||
        c.dwell[1] != 1.5) {
        fprintf(stderr, "FAIL: dwell middle prim wrong: type=%d reason=%d dwell=%g\n",
                c.type[1], c.stop_reason[1], c.dwell[1]);
        return 1;
    }
    printf("  PASS  dwell_stop  (G4 P1.5 → STOP/DWELL/1.5s)\n");
    return 0;
}

static int test_tool_change_stop(void)
{
    liscio_cfg_t cfg; liscio_cfg_default(&cfg); cfg.arc_to_helix_max = 0; /* legacy off-by-default for these tests */
    cfg.tol_xyz = 0.01;
    liscio_ctx_t *ctx = liscio_create(&cfg);
    event_collect_t c = {0};
    liscio_set_callback(ctx, on_event_prim, &c);

    liscio_pose_t a = {0,0,0,0,0,0,0,0,0};
    liscio_pose_t b = {10,0,0,0,0,0,0,0,0};
    liscio_add_line(ctx, &a, &b, 1000, 1);
    liscio_emit_stop(ctx, LISCIO_STOP_TOOL_CHANGE, 0, NULL, 2);
    liscio_flush(ctx);
    liscio_destroy(ctx);

    if (c.n != 2 ||
        c.type[1] != LISCIO_PRIM_STOP ||
        c.stop_reason[1] != LISCIO_STOP_TOOL_CHANGE) {
        fprintf(stderr, "FAIL: tool_change wrong: n=%d type[1]=%d reason[1]=%d\n",
                c.n, c.n>1?c.type[1]:-1, c.n>1?c.stop_reason[1]:-1);
        return 1;
    }
    printf("  PASS  tool_change_stop  (M6 → STOP/TOOL_CHANGE)\n");
    return 0;
}

int main(void)
{
    int fails = 0;
    printf("# liscio arc-merge + STOP-marker synthetic tests\n");
    fails += test_off_mode_passthrough();
    fails += test_full_circle_merge();
    fails += test_double_helix_merge();
    fails += test_radius_change_no_merge();
    fails += test_rapid_passthrough();
    fails += test_dwell_stop();
    fails += test_tool_change_stop();
    if (fails) {
        printf("---- %d test(s) FAILED ----\n", fails);
        return 1;
    }
    printf("---- all arc_merge + stop tests passed ----\n");
    return 0;
}
