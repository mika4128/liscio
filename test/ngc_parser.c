/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 杨阳 (Yang Yang) <mika-net@outlook.com>
 *
 * NGC parser = thin wrapper around LinuxCNC's own `rs274` interpreter.
 *
 * Why rs274:
 *   Writing a correct G-code parser (G20/G21 units, G17/G18/G19 planes,
 *   G90/G91 incremental, G90.1/G91.1 IJK modes, O-word flow control,
 *   named variables, expression evaluation, tool offsets, work coord
 *   systems…) is a several-thousand-line project.  LinuxCNC already ships
 *   a validated standalone interpreter at bin/rs274 that expands NGC into
 *   a canonical move stream.  We popen() it and translate its output.
 *
 * Canonical output lines we care about:
 *   USE_LENGTH_UNITS(CANON_UNITS_MM | CANON_UNITS_INCHES)
 *   SELECT_PLANE(CANON_PLANE_XY | CANON_PLANE_XZ | CANON_PLANE_YZ)
 *   SET_FEED_RATE(f)
 *   STRAIGHT_TRAVERSE(x, y, z, a, b, c)
 *   STRAIGHT_FEED     (x, y, z, a, b, c[, u, v, w])
 *   ARC_FEED(first_end, second_end, first_axis, second_axis,
 *            rotation, axis_end_point, a, b, c)
 *
 * Coordinates are emitted in the currently-active units — if
 * USE_LENGTH_UNITS(CANON_UNITS_INCHES) is in effect we multiply linear
 * axes by 25.4 so the ngc_move_t stream is always in mm.
 *
 * ARC_FEED plane conventions (LinuxCNC canon.hh):
 *   G17 (XY): first=X_end, second=Y_end, first_axis=X_center,
 *             second_axis=Y_center, axis_end_point=Z_end
 *   G18 (XZ): first=Z_end, second=X_end, first_axis=Z_center,
 *             second_axis=X_center, axis_end_point=Y_end
 *   G19 (YZ): first=Y_end, second=Z_end, first_axis=Y_center,
 *             second_axis=Z_center, axis_end_point=X_end
 *   rotation: signed turns.  |rot| = full revolutions, sign = direction
 *             (positive = CCW in plane as viewed from + axis). */

#define _USE_MATH_DEFINES
#define _DEFAULT_SOURCE
#define _GNU_SOURCE
#include "ngc_parser.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------- locate rs274 ---------- */

static const char *find_rs274(void)
{
    /* 1. LISCIO_RS274 env var (absolute path). */
    const char *e = getenv("LISCIO_RS274");
    if (e && *e) return e;
    /* 2. Hard-coded path that matches this source tree layout. */
    static const char *candidates[] = {
        "/home/aitalmac/PROGRAM/CNC/linuxcnc-current/bin/rs274",
        "../../../bin/rs274",
        "../../../../bin/rs274",
        "/usr/bin/rs274",
        "rs274",   /* PATH fallback */
        NULL
    };
    struct stat st;
    for (int i = 0; candidates[i]; i++) {
        if (stat(candidates[i], &st) == 0) return candidates[i];
    }
    return "rs274";
}

/* Write a dummy tool table covering any T-number a CAM post might use.
 * rs274 refuses to run without one for programs that call T/M6. */
static const char *ensure_tool_table(void)
{
    static char path[256];
    snprintf(path, sizeof(path), "/tmp/liscio_rs274_tools.tbl");
    struct stat st;
    if (stat(path, &st) == 0) return path;
    FILE *f = fopen(path, "w");
    if (!f) return "";
    for (int t = 1; t <= 200; t++) {
        fprintf(f, "T%d P%d D0.0 ;\n", t, t);
    }
    fclose(f);
    return path;
}

/* ---------- canonical-output line parser ---------- */

/* Extract the command name + paren arg list from a line like:
 *     "   42 N1234  STRAIGHT_FEED(1.0, 2.0, ...)"
 * Returns pointers into the buffer; the name is NUL-terminated in-place. */
static int split_cmd(char *line, char **out_name, char **out_args)
{
    char *paren = strchr(line, '(');
    if (!paren) return 0;
    char *e = paren;
    /* Walk backwards over the command name (alphabet / digits / underscore). */
    char *s = paren;
    while (s > line && (isalnum((unsigned char)s[-1]) || s[-1] == '_')) s--;
    if (s == paren) return 0;
    *paren = 0;     /* end-of-name */
    *out_name = s;
    /* Args: strip closing ')'. */
    char *close = strrchr(paren+1, ')');
    if (close) *close = 0;
    *out_args = paren + 1;
    *out_args = e + 1;  /* after the now-NULd '(' */
    return 1;
}

static int parse_args_doubles(char *args, double *out, int max)
{
    int n = 0;
    char *p = args;
    while (*p && n < max) {
        while (*p && (*p == ' ' || *p == ',')) p++;
        if (!*p) break;
        char *ep;
        double v = strtod(p, &ep);
        if (ep == p) break;
        out[n++] = v;
        p = ep;
    }
    return n;
}

/* ---------- main entry ---------- */

int ngc_parse_file(const char *path, ngc_move_cb cb, void *user,
                   ngc_stats_t *stats_out)
{
    ngc_stats_t stats; memset(&stats, 0, sizeof(stats));

    const char *rs274 = find_rs274();
    const char *tools = ensure_tool_table();

    /* rs274 standalone hardcodes axis_mask=0x3f (XYZABC only, no UVW) in
     * src/emc/sai/saicanon.cc:GET_EXTERNAL_AXIS_MASK().  Any G-line with a
     * W (or U/V) field triggers "Bad character 'w'" and the entire line is
     * dropped — drilling/engraving programs that mirror Z onto W lose 90%+
     * of their motion blocks.  Workaround: pre-process input through sed
     * to strip U/V/W fields and comment out unknown M-codes (M111 etc.) so
     * rs274 sees a 6-axis-only stream.  liscio's own pipeline is 9D, but
     * here we only feed the rs274-extractable subset; in CAM files that
     * mirror W onto Z (the common case) Z motion is preserved as-is. */
    char cmd[4096];
    snprintf(cmd, sizeof(cmd),
        "TMPF=$(mktemp /tmp/liscio_XXXXXX.ngc); "
        "sed -E "
        "-e 's/[UuVvWw][[:space:]]*[-+]?[0-9]*\\.?[0-9]+//g' "
        "-e 's/^[[:space:]]*[Mm]111\\b/(M111)/' "
        "\"%s\" > \"$TMPF\" && "
        "\"%s\" -g -n 0 -t \"%s\" \"$TMPF\" 2>/dev/null; "
        "ret=$?; rm -f \"$TMPF\"; exit $ret",
        path, rs274, tools);
    FILE *f = popen(cmd, "r");
    if (!f) return -1;

    /* Modal state derived from canonical stream. */
    int    units_inch = 0;     /* default CANON_UNITS_MM per rs274 init line */
    int    plane = 0;          /* 0=XY, 1=XZ, 2=YZ */
    double feed = 1000.0;
    ngc_pose_t pos = {0};
    int    first_move = 1;

    char line[2048];
    int  lineno = 0;
    while (fgets(line, sizeof(line), f)) {
        lineno++;
        /* Trim trailing \n / \r. */
        size_t L = strlen(line);
        while (L > 0 && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = 0;
        stats.n_lines_parsed++;

        char *name, *args;
        if (!split_cmd(line, &name, &args)) continue;

        double s = units_inch ? 25.4 : 1.0;

        if (!strcmp(name, "USE_LENGTH_UNITS")) {
            if (strstr(args, "INCHES")) { units_inch = 1; stats.unit_scale_used = 1; }
            else                        { units_inch = 0; }
            continue;
        }
        if (!strcmp(name, "SELECT_PLANE")) {
            if (strstr(args, "XZ"))      { plane = 1; stats.plane_g18_used = 1; }
            else if (strstr(args, "YZ")) { plane = 2; stats.plane_g19_used = 1; }
            else                         { plane = 0; }
            continue;
        }
        if (!strcmp(name, "SET_FEED_RATE")) {
            double a[1];
            if (parse_args_doubles(args, a, 1) == 1) feed = a[0] * s;
            continue;
        }

        if (!strcmp(name, "STRAIGHT_TRAVERSE") || !strcmp(name, "STRAIGHT_FEED")) {
            double a[9] = {0};
            int got = parse_args_doubles(args, a, 9);
            if (got < 6) continue;
            ngc_pose_t end = pos;
            end.x = a[0] * s; end.y = a[1] * s; end.z = a[2] * s;
            end.a = a[3];     end.b = a[4];     end.c = a[5];
            if (got >= 9) { end.u = a[6]*s; end.v = a[7]*s; end.w = a[8]*s; }

            if (first_move) { pos = end; first_move = 0; continue; }

            ngc_move_t m; memset(&m, 0, sizeof(m));
            m.type = NGC_MOVE_LINE;
            m.start = pos; m.end = end;
            m.feedrate = feed;
            m.is_rapid = (name[9] == 'T');   /* STRAIGHT_TRAVERSE vs STRAIGHT_FEED */
            m.line_no = lineno;
            m.plane = plane;
            if (m.is_rapid) stats.n_g0++; else stats.n_g1++;
            if (cb) cb(&m, user);
            pos = end;
            continue;
        }

        if (!strcmp(name, "ARC_FEED")) {
            double a[9] = {0};
            int got = parse_args_doubles(args, a, 9);
            if (got < 6) continue;
            /* Unpack plane-dependent canonical ARC_FEED args. */
            double first_end = a[0] * s;
            double second_end = a[1] * s;
            double first_axis = a[2] * s;
            double second_axis = a[3] * s;
            int    rotation = (int)a[4];
            double axis_end = a[5] * s;
            double ea = got > 6 ? a[6] : pos.a;
            double eb = got > 7 ? a[7] : pos.b;
            double ec = got > 8 ? a[8] : pos.c;

            ngc_pose_t end = pos;
            end.a = ea; end.b = eb; end.c = ec;

            ngc_move_t m; memset(&m, 0, sizeof(m));
            m.type = NGC_MOVE_ARC;
            m.start = pos;
            m.feedrate = feed;
            m.line_no = lineno;
            m.plane = plane;

            /* Plane-specific mapping → m.end and (cx, cy, cz) in 3D.
             *
             * Convention: normal always +z (G17), +y (G18), +x (G19) — the
             * rotation direction is carried by the SIGN of arc_angle alone,
             * never by flipping the normal.  Downstream arc-sampling uses
             * Rodrigues rotation about the normal by arc_angle, so flipping
             * both would double-negate. */
            double sx_pl, sy_pl, ex_pl, ey_pl;    /* coords in plane-local 2D */
            double cx_pl, cy_pl;                   /* center in 2D */
            if (plane == 0) {         /* G17: first=X second=Y, axis_end=Z */
                end.x = first_end; end.y = second_end; end.z = axis_end;
                m.cx = first_axis; m.cy = second_axis; m.cz = pos.z;
                sx_pl = pos.x; sy_pl = pos.y; ex_pl = end.x; ey_pl = end.y;
                cx_pl = m.cx; cy_pl = m.cy;
                m.nx = 0; m.ny = 0; m.nz = 1.0;
            } else if (plane == 1) {  /* G18: first=Z second=X, axis_end=Y */
                end.z = first_end; end.x = second_end; end.y = axis_end;
                m.cz = first_axis; m.cx = second_axis; m.cy = pos.y;
                /* Plane-local basis for G18 in LinuxCNC canon: first axis Z,
                 * second axis X.  Keep that orientation so the signed
                 * arc_angle derived below matches canon rotation sign. */
                sx_pl = pos.z; sy_pl = pos.x; ex_pl = end.z; ey_pl = end.x;
                cx_pl = m.cz; cy_pl = m.cx;
                m.nx = 0; m.ny = 1.0; m.nz = 0;
            } else {                   /* G19: first=Y second=Z, axis_end=X */
                end.y = first_end; end.z = second_end; end.x = axis_end;
                m.cy = first_axis; m.cz = second_axis; m.cx = pos.x;
                sx_pl = pos.y; sy_pl = pos.z; ex_pl = end.y; ey_pl = end.z;
                cx_pl = m.cy; cy_pl = m.cz;
                m.nx = 1.0; m.ny = 0; m.nz = 0;
            }

            /* Radius + signed arc_angle: atan2 gives principal angle (-π,π],
             * then use rs274's signed `rotation` to disambiguate >π arcs. */
            double dsx = sx_pl - cx_pl, dsy = sy_pl - cy_pl;
            double dex = ex_pl - cx_pl, dey = ey_pl - cy_pl;
            m.radius = sqrt(dsx*dsx + dsy*dsy);
            double ang = atan2(dsx*dey - dsy*dex, dsx*dex + dsy*dey);
            if (rotation > 0 && ang < 0) ang += 2.0 * M_PI;
            if (rotation < 0 && ang > 0) ang -= 2.0 * M_PI;
            if (abs(rotation) > 1) {
                int extra = abs(rotation) - 1;
                ang += (rotation > 0 ? 1.0 : -1.0) * 2.0 * M_PI * extra;
            }
            m.arc_angle = ang;
            m.end = end;

            if (first_move) { pos = end; first_move = 0; continue; }
            stats.n_g2g3++;
            if (cb) cb(&m, user);
            pos = end;
            continue;
        }

        /* Other canon ops (COMMENT, SELECT_TOOL, START_SPINDLE_*, etc.):
         * ignored. */
    }

    int rc = pclose(f);
    (void)rc;   /* rs274 may non-zero-exit on some files even after emit */
    if (stats_out) *stats_out = stats;
    return 0;
}
