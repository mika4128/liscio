**English** | [中文](README.zh.md)

# liscio

<p align="center">
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-green.svg" alt="MIT"></a>
  <img src="https://img.shields.io/badge/lang-C99-blue.svg" alt="C99">
  <img src="https://img.shields.io/badge/deps-none-brightgreen.svg" alt="no deps">
</p>

> **Pronunciation**: Italian /ˈliʃʃo/ — `liscio` means "smooth", which fits the trajectory-smoothing preprocessor theme.

**Pure C99 geometry preprocessor for CNC trajectory planners — micro-segments → arc / cubic Bezier / rational Bezier / B-spline fit.**

## Features

- Pure C99, zero third-party deps, portable to bare metal / RT kernel (LinuxCNC RTAPI).
- **9D (XYZ + ABC + UVW)**, per-axis independent fit, 5-axis CAM friendly.
- Primitives: **LINE**, **ARC** (Kasa LSQ), **cubic Bezier** (9D per-axis LSQ),
  **rational Bezier** (NURBS), **B-spline** (clamped cubic).
- Composite recursive subdivision + Hoschek reparameterization, approaching Fanuc AICC II / Siemens COMPCURV accuracy.
- Automatic corner split, tangent-vector output, streaming API, fixed memory (no heap allocation).
- Throughput 14M waypoints / second (x86 single core), typical CAM compression 10–25×.

## Industry comparison

| Feature | Fanuc AICC II | Siemens 840D | liscio |
|---|---|---|---|
| Micro-segment → geometry compression | ✅ Nano Smoothing | ✅ COMPCURV | ✅ |
| 9D (incl. rotary/uvw) | ✅ | ✅ | ✅ |
| Cubic Bezier LSQ | ✅ | ✅ | ✅ |
| Composite recursion | ✅ | ✅ | ✅ |
| Hoschek reparameterization | ✅ | ✅ | ✅ |
| Auto corner split | ✅ | G641/642 | ✅ |
| Rational weights | ✅ | ✅ | ✅ |
| Exact arc → NURBS | ✅ | ✅ | ✅ |
| Free-form B-spline LSQ | ✅ | ✅ | ✅ |

## Recommended configuration

For most users, **only `tol_xyz` needs tuning**; keep other fields at defaults.

```c
liscio_cfg_t cfg;
liscio_cfg_default(&cfg);
cfg.tol_xyz            = 0.01;     /* 10 µm — industry sweet spot */
cfg.tol_abc            = 0.1;      /* 0.1° — rotary axes */
cfg.arc_subseg_samples = 0;        /* G2/G3 passthrough, lossless (default) */
liscio_ctx_t *ctx = liscio_create(&cfg);
```

### How to pick `tol_xyz`

Match LinuxCNC `[TRAJ] DEFAULT_LINEAR_TOLERANCE` / G-code `G64 P`:

| Scenario | Suggested tol | Measured compression | Note |
|---|---:|---:|---|
| Wood / plastic / 3D printing | 0.05 mm | 7.49× | CAM lower bound |
| **Metal cutting mainstream** | **0.01 mm** | **4.43×** | **Recommended default ✅** |
| Mold / optics / semiconductor | 1e-3 mm | 2.02× | Ultra-precision |
| Form-tolerance zero | ≤1e-4 mm | ≈ 1.4× | Near passthrough |

> liscio runs down to double-precision floor `1e-8 mm` with no performance penalty.
> At tight tolerance, compression ratio approaches the 1.40× asymptote (≈ geometric passthrough).

### G2/G3 handling: **passthrough** (`arc_subseg_samples = 0`)

CAM G2/G3 arcs go through `liscio_add_arc()` directly → emit `LISCIO_PRIM_ARC`,
**max_dev = 0, lossless geometry**. Downstream TP using analytic arc formulas
for jerk planning is much faster than numerical integration along Bezier.

`arc_subseg_samples > 0` (split mode) only when **all three** conditions hold:
1. CAM emits many short arcs (e.g. biarc fillet)
2. Workpiece tolerance ≥ 10 µm
3. Downstream TP buffer is tight, must reduce primitive count

Empirical: at tight tolerance, split mode **inflates** up to 5×. Default 0 is
the industrial optimum.

### Advanced: adjacent arc merge + HELIX (opt-in)

```c
liscio_cfg_set_arc_merge(&cfg, LISCIO_ARC_MERGE_HELIX, 0, 0);
```

Fold adjacent concentric + same-radius + same-pitch G2/G3 arcs into a single ARC
(or HELIX when arc_angle > 2π). Lossless. Use cases:

| Applicable | Not applicable |
|---|---|
| Hand-written G-code with concentric arcs | CAM chord-approximation |
| Thread cutting (lathe G33) | Free-form surface cuts |
| Multi-layer concentric drilling | Spiral-wall machining (R diminishing) |

### Advanced: Look-ahead corner decision (opt-in)

```c
liscio_cfg_set_corner_detection(&cfg, LISCIO_CORNER_LOOKAHEAD,
                                 15.0 /*soft°*/, 60.0 /*hard°*/);
```

Decouple "angle judgment" from "split-or-not", two thresholds:
- **Soft corner** (15-60°): don't split immediately, let Bezier G1 fit absorb it → tangent continuity improves significantly
- **Hard corner** (>60°): split immediately (preserve corner intent)

Empirical (16 files, tol=0.05): weighted tan_mean **20.4° → 14.4° (-30%)**, cost
compression ratio 7.49× → 6.05× (-19%). Best gain on flower (-40%) / flower-one-line (-55%).

Applicable: when downstream does **G64 blend**, jerk-limited speed planning becomes smoother.
Not applicable: G61 exact-stop + tight prim-count budget. Default still IMMEDIATE.

### G-code event stream: RAPID + STOP passthrough

liscio output = full G-code event stream, not only cutting primitives. Non-cutting
events get **dedicated APIs**:

```c
liscio_add_rapid(ctx, &start, &end, line);                   /* G0 */
liscio_emit_stop(ctx, LISCIO_STOP_DWELL, 1.5, &pos, line);   /* G4 P1.5 */
liscio_emit_stop(ctx, LISCIO_STOP_TOOL_CHANGE, 0, NULL, line); /* M6 */
liscio_emit_stop(ctx, LISCIO_STOP_PROGRAM_END, 0, NULL, line); /* M2/M30 */
```

Output `LISCIO_PRIM_RAPID` / `LISCIO_PRIM_STOP` primitives; downstream TP
single entry (`on_prim`) dispatches by `p->type`: cutting primitives go to
the jerk-limited queue, RAPID goes max-velocity, STOP handled per `p->stop_reason`.

Full LinuxCNC canon (STRAIGHT_FEED/STRAIGHT_TRAVERSE/DWELL/CHANGE_TOOL/
SET_ORIGIN_OFFSETS/...) → liscio API mapping table at

### Configuration flowchart

```
What is your G64 P?
  ├─ unset ─────────────→ tol_xyz = 0.01 (LinuxCNC default sweet spot)
  ├─ 0.05 (wood/plastic) → tol_xyz = 0.05
  ├─ 0.01 (metal) ──────→ tol_xyz = 0.01
  └─ 0.001 (mold) ──────→ tol_xyz = 0.001

CAM file with many G2/G3 arcs?
  ├─ no (95% of files) ──→ arc_subseg_samples = 0  ✅
  └─ yes + tol≥0.01 ────→ try arc_subseg_samples = 8 (compare empirically)
```

## Quick start

```bash
cmake -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Produces `libliscio.a` + `libliscio.so`, all tests PASS.

### Minimal usage

```c
#include <liscio/liscio.h>

static void on_prim(const liscio_primitive_t *p, void *user) {
    /* p->type = LINE / ARC / BEZIER / RBEZIER / SPLINE
     * Consume p->start / p->end / P0..P3 / cx,cy,... */
}

liscio_cfg_t cfg;
liscio_cfg_default(&cfg);           /* tol_xyz=0.01, arc passthrough */
liscio_ctx_t *ctx = liscio_create(&cfg);
liscio_set_callback(ctx, on_prim, NULL);

for (/* each G01 segment */) {
    liscio_pose_t s = {...}, e = {...};
    liscio_add_line(ctx, &s, &e, feedrate, line_no);
}
/* G2/G3 arcs passthrough via add_arc: */
for (/* each G02/G03 segment */) {
    liscio_primitive_t arc = {...};      /* cx,cy,radius,arc_angle,normal */
    liscio_add_arc(ctx, &arc, feedrate, line_no);
}
liscio_flush(ctx);
liscio_destroy(ctx);
```

### CMake integration

```cmake
add_subdirectory(liscio)
target_link_libraries(your_target PRIVATE liscio m)
```

Or `find_package(liscio CONFIG REQUIRED)` after install.

### Visualizing fit results

```bash
cmake --build build --target dump_ngcs    # test/data/*.ngc → CSV
cmake --build build --target plot_ngcs    # CSV → PNG (needs matplotlib)
ls build/ngc_out/*.png
# grey dots = raw G01 / blue = LINE / orange = ARC / green = BEZIER / pink = SPLINE
```

Single file:
```bash
./build/test_dump_csv --tol_xyz=0.05 --samples=30 foo.ngc
python3 test/plot_fit.py --3d foo.ngc
```

## Layout

```
liscio/
├── include/liscio/liscio.h   public API
├── src/                       fit algorithm impl (arc/bezier/rbezier/bspline)
├── test/                      unit tests + CAM regression + visualization
└── results/                   test results (PNG / CSV / sweep / cross-lib)
```

## Tools

`liscio_compress` — **G-code → G-code compressor**, cross-controller compatible.
Compresses CAM-output NGC files; typical machining files achieve **1.6-1.9× file-size
compression**:

```bash
./build/liscio_compress --tol_xyz=0.01 --output=part.ngc.gz part.ngc
```

See [tool/README.md](tool/README.md).

## Adding new test samples

Drop `.ngc` into [test/data/](test/data/) and re-run `cmake -B build` (CMake
`file(GLOB)` re-scans). Auto-included:

- `ctest verify_geometry` — fit-quality gate (OVER_TOL=0, DRIFT=0)
- `make dump_ngcs` / `make plot_ngcs` — CSV + PNG visualization
- `./test_metrics --csv=report.csv test/data/*.ngc` — quantitative metrics

## License

MIT (SPDX `MIT`). See [LICENSE](LICENSE).
