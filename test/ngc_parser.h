/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 杨阳 (Yang Yang) <mika-net@outlook.com>
 *
 * Shared NGC parser for liscio test tools.  Handles:
 *   - G0 rapids, G1 feeds, G2/G3 arcs (IJK or R form)
 *   - G17 / G18 / G19 plane select (affects arc center + normal)
 *   - G20 (inch) / G21 (mm) — all coords normalized to MM in output
 *   - G90 (abs) / G91 (incremental) position
 *   - G90.1 (abs IJK) / G91.1 (incremental IJK)
 *   - O-word while/endwhile loops with numeric #var counter
 *   - Comments: (…) and ;…
 *   - Blank/M/T/S lines skipped
 *
 * Output: stream of ngc_move_t (line or arc) via callback.  Coords in mm.
 * G0 rapids are reported with is_rapid=1; consumer decides how to treat.
 */

#ifndef LISCIO_NGC_PARSER_H
#define LISCIO_NGC_PARSER_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double x, y, z, a, b, c, u, v, w;
} ngc_pose_t;

typedef enum {
    NGC_MOVE_LINE = 0,
    NGC_MOVE_ARC  = 1,
} ngc_move_type_t;

typedef struct {
    ngc_move_type_t type;
    ngc_pose_t      start;
    ngc_pose_t      end;
    double          feedrate;    /* mm/min */
    int             line_no;
    int             is_rapid;    /* G0 rapid vs G1/G2/G3 */
    int             plane;       /* 0=G17 XY, 1=G18 XZ, 2=G19 YZ */
    /* Arc-only fields: */
    double          cx, cy, cz;
    double          radius;
    double          arc_angle;   /* signed */
    double          nx, ny, nz;  /* plane normal */
} ngc_move_t;

typedef void (*ngc_move_cb)(const ngc_move_t *m, void *user);

typedef struct {
    /* Returned after parse: statistics. */
    long n_g0;
    long n_g1;
    long n_g2g3;
    long n_lines_parsed;
    int  unit_scale_used;    /* 1 = inches → mm applied somewhere */
    int  plane_g18_used;
    int  plane_g19_used;
    int  g91_used;
    int  owhile_iters;
} ngc_stats_t;

int ngc_parse_file(const char *path, ngc_move_cb cb, void *user,
                   ngc_stats_t *stats_out);

#ifdef __cplusplus
}
#endif

#endif
