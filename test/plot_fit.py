#!/usr/bin/env python3
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction...
#
# Plot liscio fit result: raw waypoints (gray) + fitted primitives
# color-coded by type.
#
# Usage:
#   test_dump_csv foo.ngc                  # writes foo.ngc.{input,prims}.csv
#   python3 plot_fit.py foo.ngc            # plots both CSVs
#
# License: MIT  Author: 杨阳 (Yang Yang) <mika-net@outlook.com>

import sys
import csv
from pathlib import Path

try:
    import matplotlib.pyplot as plt
except ImportError:
    sys.stderr.write("matplotlib not installed; try: pip install matplotlib\n")
    sys.exit(2)

COLOR = {
    "LINE":   "#2196F3",   # blue
    "ARC":    "#FF9800",   # orange
    "BEZIER": "#4CAF50",   # green
    "SPLINE": "#E91E63",   # pink
    "HELIX":  "#9C27B0",   # purple — multi-turn arcs folded by helix9 LSQ
}


def load_input(path):
    rows = []
    with open(path) as f:
        r = csv.DictReader(f)
        for row in r:
            rows.append((float(row["x"]), float(row["y"]), float(row["z"])))
    return rows


def load_prims(path):
    by_prim = {}
    with open(path) as f:
        r = csv.DictReader(f)
        for row in r:
            pid = int(row["prim"])
            t = row["type"]
            by_prim.setdefault(pid, (t, []))
            by_prim[pid][1].append(
                (float(row["x"]), float(row["y"]), float(row["z"])))
    return by_prim


def load_rapids(path):
    rows = []
    if not path.exists():
        return rows
    with open(path) as f:
        r = csv.DictReader(f)
        for row in r:
            rows.append((
                float(row["sx"]), float(row["sy"]), float(row["sz"]),
                float(row["ex"]), float(row["ey"]), float(row["ez"])))
    return rows


def plot_one(ngc_path, threed=False):
    base = Path(ngc_path)
    in_csv = base.with_suffix(base.suffix + ".input.csv")
    pr_csv = base.with_suffix(base.suffix + ".prims.csv")
    rp_csv = base.with_suffix(base.suffix + ".rapids.csv")
    if not in_csv.exists() or not pr_csv.exists():
        sys.stderr.write(f"missing CSVs next to {ngc_path}\n")
        return

    raw = load_input(in_csv)
    prims = load_prims(pr_csv)
    rapids = load_rapids(rp_csv)

    # Always render 3D — easier to read CAM retract/plunge moves.
    if not threed and raw:
        threed = True

    if threed:
        from mpl_toolkits.mplot3d import Axes3D  # noqa

        xs_all = [p[0] for p in raw]
        ys_all = [p[1] for p in raw]
        zs_all = [p[2] for p in raw]
        if xs_all and ys_all and zs_all:
            rx = max(xs_all) - min(xs_all)
            ry = max(ys_all) - min(ys_all)
            rz = max(zs_all) - min(zs_all)
        else:
            rx = ry = rz = 1
        m = max(rx, ry, rz) or 1
        floor = m * 0.01
        aspect = (max(rx, floor) / m,
                  max(ry, floor) / m,
                  max(rz, floor) / m)

        fig = plt.figure(figsize=(20, 17))

        # Three orthographic projections as 2D — avoids mpl-3D box_aspect
        # axis-limit warping when one dimension is much smaller than others.
        # Index 0,1,2: Top (XY), Front (XZ), Side (YZ).
        ORTHO = [
            ("Top  (XY)",   0, 1, 'X', 'Y'),
            ("Front (XZ)",  0, 2, 'X', 'Z'),
            ("Side (YZ)",   1, 2, 'Y', 'Z'),
        ]
        for i, (name, ai, bi, alab, blab) in enumerate(ORTHO):
            ax = fig.add_subplot(2, 2, i + 1)
            # NOTE: skip G0 rapids in orthographic views — they pull axis
            # range far beyond cutting region and shrink the part of interest.
            # Raw polyline (light grey).
            ra = [p[ai] for p in raw]
            rb = [p[bi] for p in raw]
            ax.plot(ra, rb, color='#bbbbbb', linewidth=0.6, alpha=0.7,
                    label='raw G1' if i == 0 else None)
            types_seen = set()
            for pid, (t, pts) in prims.items():
                xa = [p[ai] for p in pts]
                xb = [p[bi] for p in pts]
                lbl = t if (i == 0 and t not in types_seen) else None
                types_seen.add(t)
                ax.plot(xa, xb, color=COLOR.get(t, 'k'),
                        linewidth=0.8, alpha=0.85, label=lbl)
            ax.set_xlabel(alab); ax.set_ylabel(blab)
            # 'box' adjusts subplot shape, not axis limits — preserves the
            # real Z range on flat parts (datalim mode would inflate Z to
            # match X span and make a 18 mm-deep cut look 150 mm deep).
            ax.set_aspect('equal', adjustable='box')
            ax.grid(True, alpha=0.3)
            ax.set_title(name, fontsize=10)

        # Iso 3D in subplot 4.  Skip G0 rapids here too — they pull axis
        # range far beyond the cutting region.
        ax3d = fig.add_subplot(2, 2, 4, projection='3d')
        ax3d.scatter(xs_all, ys_all, zs_all, c='#cccccc', s=1, alpha=0.25,
                     label='raw G1 waypoints')
        types_seen = set()
        for pid, (t, pts) in prims.items():
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            zs = [p[2] for p in pts]
            lbl = t if t not in types_seen else None
            types_seen.add(t)
            ax3d.plot(xs, ys, zs, color=COLOR.get(t, 'k'),
                      linewidth=2.2, label=lbl)
        ax3d.set_xlabel('X'); ax3d.set_ylabel('Y'); ax3d.set_zlabel('Z')
        ax3d.view_init(elev=30, azim=-60)
        ax3d.set_title('Iso', fontsize=10)
        # Tight axis limits to raw waypoint extent (excluding G0 rapids
        # is already implicit since rapids aren't fed to scatter/plot).
        # Then set box_aspect to the real (rx, ry, rz) ratio so X:Y:Z
        # render at 1:1:1 physical scale — a 18 mm cut on a 150 mm part
        # looks like 12% of the box height, matching reality.  Floor any
        # near-zero axis to 5% of the largest so flat parts still show.
        if xs_all and ys_all and zs_all:
            ax3d.set_xlim(min(xs_all), max(xs_all))
            ax3d.set_ylim(min(ys_all), max(ys_all))
            zlo, zhi = min(zs_all), max(zs_all)
            if zhi - zlo < 1e-9:
                zlo -= 0.5; zhi += 0.5
            ax3d.set_zlim(zlo, zhi)
            floor = m * 0.05
            try:
                ax3d.set_box_aspect((max(rx, floor),
                                     max(ry, floor),
                                     max(rz, floor)))
            except AttributeError:
                pass
        ax = ax3d   # legend + suptitle anchor
    else:
        fig, ax = plt.subplots(figsize=(14, 11))
        # Rapids first (background)
        r_lbl = 'G0 rapid'
        for s in rapids:
            ax.plot([s[0], s[3]], [s[1], s[4]],
                    color='#e0e0e0', linestyle='--', linewidth=0.6, label=r_lbl)
            r_lbl = None
        # Raw G1 waypoints as a connected polyline — the TRUE input shape.
        # Solid black, thicker, slightly transparent so fit overlays show.
        rx = [p[0] for p in raw]; ry = [p[1] for p in raw]
        ax.plot(rx, ry, color='black', linewidth=1.0, alpha=0.85,
                label='raw G1 polyline')
        types_seen = set()
        for pid, (t, pts) in prims.items():
            xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
            lbl = t if t not in types_seen else None
            types_seen.add(t)
            ax.plot(xs, ys, color=COLOR.get(t, 'k'),
                    linewidth=1.4, alpha=0.55, label=lbl)
        ax.set_xlabel('X'); ax.set_ylabel('Y')
        ax.set_aspect('equal', adjustable='datalim')
        ax.grid(True, alpha=0.3)

    stats = {}
    for pid, (t, _) in prims.items():
        stats[t] = stats.get(t, 0) + 1
    stat_str = " ".join(f"{t}:{n}" for t, n in sorted(stats.items()))
    title = f"{base.name}   waypoints:{len(raw)}   {stat_str}"
    if threed:
        fig.suptitle(title, fontsize=12)
        ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0), fontsize=8)
    else:
        ax.set_title(title, fontsize=11)
        ax.legend(loc='best', fontsize=9)

    out = base.with_suffix(base.suffix + ".png")
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close(fig)
    print(f"saved {out}")


def main():
    args = sys.argv[1:]
    threed = False
    files = []
    for a in args:
        if a in ('-3', '--3d'):
            threed = True
        else:
            files.append(a)
    if not files:
        sys.stderr.write("usage: plot_fit.py [--3d] <ngc>...\n")
        sys.exit(1)
    for p in files:
        plot_one(p, threed=threed)


if __name__ == "__main__":
    main()
