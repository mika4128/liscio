#!/usr/bin/env python3
"""
geomdl_fit.py — subprocess helper for liscio cross-comparison benchmark.

Reads a polyline from stdin and runs geomdl's LSQ B-spline approximation
(geomdl.fitting.approximate_curve, Piegl & Tiller §9.4 reference impl).

Stdin protocol (text, one number per line, in order):
    n_points
    degree
    ctrlpts_size_target   (0 = auto, else explicit M)
    x0  y0  z0
    x1  y1  z1
    ...

Stdout protocol (text):
    line 1:  N_CTRLPTS  MAX_DEV  ELAPSED_MS
    lines 2..N_CTRLPTS+1:  cpx cpy cpz  (one control point per line)
"""

import sys
import time
import math

def read_input():
    f = sys.stdin
    n   = int(f.readline().strip())
    deg = int(f.readline().strip())
    sz  = int(f.readline().strip())
    pts = []
    for _ in range(n):
        parts = f.readline().split()
        pts.append((float(parts[0]), float(parts[1]), float(parts[2])))
    return pts, deg, sz

def main():
    try:
        from geomdl import fitting
        from geomdl import operations
    except Exception as e:
        sys.stderr.write(f"geomdl import failed: {e}\n")
        sys.exit(2)

    pts, deg, sz_target = read_input()
    n = len(pts)
    if n < 4:
        sys.stderr.write("geomdl: need >= 4 points\n")
        print("0 0.0 0.0")
        return

    # Heuristic ctrlpts size: max( deg+1 , min( N//4, 200 ) ) for fairness
    if sz_target <= 0:
        sz = max(deg + 1, min(max(4, n // 4), 200))
    else:
        sz = max(deg + 1, min(sz_target, n))

    t0 = time.perf_counter()
    try:
        curve = fitting.approximate_curve(pts, degree=deg, ctrlpts_size=sz)
    except Exception as e:
        sys.stderr.write(f"geomdl fit failed: {e}\n")
        print("0 0.0 0.0")
        return
    elapsed_ms = (time.perf_counter() - t0) * 1000.0

    # Sample curve at moderate resolution to measure max deviation
    # against original polyline (point-to-curve closest distance).
    # Approximate: project each input point to nearest curve sample.
    n_eval = max(200, n)
    curve.delta = 1.0 / n_eval
    eval_pts = curve.evalpts  # list of [x,y,z]

    max_dev = 0.0
    # Brute-force NN search vs sampled curve. O(N * n_eval), fine for
    # benchmark sizes <= 50k.
    j_hint = 0
    for p in pts:
        # Local search: start near previous hint to keep O(N+n_eval).
        best = None
        # Bounded scan around hint (window 64), wraps to full scan if missing.
        lo = max(0, j_hint - 32)
        hi = min(n_eval, j_hint + 32)
        for j in range(lo, hi):
            q = eval_pts[j]
            d2 = (p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2
            if best is None or d2 < best[0]:
                best = (d2, j)
        # Fallback global scan if local hit window edge.
        if best[1] == lo or best[1] == hi - 1:
            for j in range(0, n_eval):
                q = eval_pts[j]
                d2 = (p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2
                if d2 < best[0]:
                    best = (d2, j)
        j_hint = best[1]
        d = math.sqrt(best[0])
        if d > max_dev: max_dev = d

    n_ctrl = len(curve.ctrlpts)
    print(f"{n_ctrl} {max_dev:.6f} {elapsed_ms:.3f}")
    for cp in curve.ctrlpts:
        print(f"{cp[0]:.6f} {cp[1]:.6f} {cp[2]:.6f}")

if __name__ == "__main__":
    main()
