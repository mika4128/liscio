**English** | [中文](README.zh.md)

# liscio test results

Auto-generated regression archive — refresh with:

```bash
cmake --build build --target liscio_results
```

This rebuilds CSV/PNGs from `test/data/*.ngc` and copies them here.
The CMake `metrics_csv` target also runs as part of `ALL`, so
`data/metrics.csv` and `data/metrics_g1.csv` stay current on every
build.

## Layout

```
results/
├── README.md            ← this file
├── images/              ← PNG fit visualization, one per ngc
│   ├── thomam.ngc.png   ← helix detection demo (purple = HELIX prims)
│   └── ...
├── data/
│   ├── metrics.csv      ← tol_xyz=0.01, default helix detection ON
│   ├── metrics_g1.csv   ← tol_xyz=0.01, arc_subseg=8 (chord-sampled)
│   ├── metrics.txt      ← human-readable snapshot (older run)
│   └── xcompare.txt     ← 5-library cross-comparison log
├── sweep/               ← tol-sweep artefacts (8 tol × 2 modes)
                         ← phase-feature deep-dives
```

## Configuration used

`metrics.csv` and the PNG dumps run with default tol_xyz, helix
auto-detection on (`cfg.arc_to_helix_max = 16`).  No file-specific
flags — every `.ngc` (or `.NGC`) in `test/data/` is fitted with the
same pipeline:

```
G1 → window LSQ → ARC | HELIX | BEZIER | SPLINE | LINE
G2/G3 → arc-merge / arc-to-helix LSQ / passthrough
```

Helix detection is automatic: the fitter buffers up to 16 consecutive
G2/G3 arcs, runs an analytic LSQ helix fit on the buffered geometry
(centers along common axis, constant pitch, constant radius), and
emits a single `LISCIO_PRIM_HELIX` when the residual ≤ tol_xyz.
Otherwise the buffered arcs flush as individual `LISCIO_PRIM_ARC`s
(passthrough fallback, max_dev=0).

## Highlights from the current run

(See `data/metrics.csv` for the full 18-file table.)

### Helix detection (thomam.ngc — lathe-thread cut)

Without arc-to-helix this file passes through as 254 ARCs.  With
auto-detection it folds the same geometry into:

```
in=276  out=80  ratio=3.45×  HELIX=14  max_dev=0.0000  tan_mean=2.6°
```

`max_dev=0` and len_err < 0.5% prove the helix LSQ is geometrically
consistent with the original arc set.  The PNG (`images/thomam.ngc.png`)
shows the multi-turn purple HELIX primitives; orange ARCs at the rim
are single-turn boundaries the LSQ rejected as non-helical.

### Compression at tol_xyz=0.01

| File | input | output | ratio | max_dev | tan_mean |
|---|---:|---:|---:|---:|---:|
| sine_test | 501 | 18 | **27.83×** | 0.005 | 2.6° |
| 5_1001 | 10 150 | 906 | 11.20× | 0.010 | 22.0° |
| 3_1001 | 11 391 | 1 034 | 11.02× | 0.010 | 15.4° |
| 1_1001 | 19 879 | 2 881 | 6.90× | 0.010 | 1.9° |
| flower | 43 997 | 14 660 | 3.00× | 0.010 | 27.5° |
| thomam | 276 | 254 | 1.09× | 0.000 | 0.7° |

(`thomam` row above is at tol=0.01 strict — at the recommended
0.05 it folds to 80 prims with 14 HELIX as shown earlier.)

## Quality gates (ctest verify_geometry)

`OVER_TOL` (waypoints farther than tol_xyz from any fitted curve)
and `DRIFT` (primitives whose fit length > 1.5× raw chord) must both
be 0 for every file.  CI runs this gate on every push.

```bash
ctest --test-dir build -R verify_geometry --output-on-failure
```

Expected output: 18 files, total `OVER_TOL=0 DRIFT=0`.

## Reproducing manually

```bash
# 1. Build
cmake -B build && cmake --build build -j

# 2. Run all unit + regression tests
ctest --test-dir build --output-on-failure

# 3. Single-file probe
./build/test_metrics --tol_xyz=0.05 test/data/thomam.ngc

# 4. Refresh PNGs + CSVs into results/
cmake --build build --target liscio_results
```

## Adding a new test sample

Drop the `.ngc` (or `.NGC`) into `test/data/`.  CMake's
`CONFIGURE_DEPENDS` glob picks it up automatically on the next build —
no list edits needed.  All targets above (verify_geometry,
dump_ngcs, plot_ngcs, metrics_csv, liscio_results) include it.

