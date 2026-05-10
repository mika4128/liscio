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
cfg.tol_xyz            = 0.01;     /* 10 µm — typical metal-cutting setting */
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
| Form-tolerance zero | ≤1e-4 mm | ≈ 1.4× | Almost no compression — near-passthrough |

> liscio still runs at the double-precision floor `1e-8 mm` with no
> slowdown.  As you tighten the tolerance, the compression ratio
> approaches its lower bound (~1.40×), which is essentially "pass
> the input through unchanged".  

### G2/G3 arc handling: **passthrough** (`arc_subseg_samples = 0`)

A CAM-emitted G2/G3 arc fed through `liscio_add_arc()` is preserved
**bit-for-bit geometrically** as a `LISCIO_PRIM_ARC` (max_dev = 0).
Downstream the trajectory planner uses the analytic arc formulas to
compute jerk-limited velocity profiles — much faster than numerically
integrating along a Bezier curve.

Set `arc_subseg_samples > 0` to instead chop each G2/G3 into N short
G1 lines and feed them to the G1 LSQ pipeline; this lets adjacent
arcs fold with neighbouring G1 runs into longer BEZIER primitives.
Only useful when **all three** are true:

1. CAM produces many short arcs (e.g. fillet biarcs)
2. Workpiece tolerance ≥ 10 µm (chord error is tolerable)
3. Downstream queue is tight on primitive count

Empirical: at tight tolerance the split mode **inflates** up to 5×
(one arc → 8 lines).  Default `0` (passthrough) is the right industrial
choice.

### Advanced: fold adjacent concentric arcs into HELIX (opt-in)

```c
liscio_cfg_set_arc_merge(&cfg, LISCIO_ARC_MERGE_HELIX, 0, 0);
```

Fold adjacent concentric + same-radius + same-pitch G2/G3 arcs into a single ARC
(or HELIX when arc_angle > 2π). Lossless. Use cases:

| Applicable | Not applicable |
|---|---|
| Hand-written G-code with several truly concentric arcs | CAM that uses many small arcs to approximate one curve |
| Thread cutting on a lathe (G33) | Free-form surface milling |
| Multi-layer concentric drilling | Spiral wall (radius shrinks per turn) |

> For CAM chord-approximated helices use the default-on `arc_to_helix_max`
> LSQ path instead.  
> It handles centres that drift in the last few floating-point digits.

### Advanced: don't split on "soft" corners (opt-in)

```c
liscio_cfg_set_corner_detection(&cfg, LISCIO_CORNER_LOOKAHEAD,
                                 15.0 /*soft°*/, 60.0 /*hard°*/);
```

Two angle thresholds, separating "is this a corner?" from "do we split?":
- **Soft corner** (15-60°): don't split — let the BEZIER fit absorb
  it.  The output curve has a continuous tangent through this region.
- **Hard corner** (>60°): split immediately, this is a real corner.

Empirical (16 files, tol=0.05): mean tangent mismatch at primitive
joins drops **20.4° → 14.4° (-30%)**, cost is compression ratio
7.49× → 6.05× (-19%).  flower (-40%) / flower-one-line (-55%) benefit
the most.

When to enable: downstream uses G64 blend (smooth corner pass-through)
and wants smoother jerk-limited velocity profiles.  When **not** to
enable: G61 exact-stop machining or a tight primitive-count budget —
keep the default of splitting at every corner.

### G-code event stream: G0, dwell, tool-change pass through too

liscio's output is **the full G-code event stream**, not only cutting
primitives.  Non-cutting events have dedicated APIs:

```c
liscio_add_rapid(ctx, &start, &end, line);                   /* G0 */
liscio_emit_stop(ctx, LISCIO_STOP_DWELL, 1.5, &pos, line);   /* G4 P1.5 */
liscio_emit_stop(ctx, LISCIO_STOP_TOOL_CHANGE, 0, NULL, line); /* M6 */
liscio_emit_stop(ctx, LISCIO_STOP_PROGRAM_END, 0, NULL, line); /* M2/M30 */
```

These emit as `LISCIO_PRIM_RAPID` / `LISCIO_PRIM_STOP` primitives.
The downstream planner takes a single `on_prim` callback and
dispatches on `p->type`:

- cutting primitives → jerk-limited velocity queue
- RAPID (G0) → max-velocity traverse
- STOP → handled by `p->stop_reason` (dwell / tool change / program
  end / coordinate-system change)

### Configuration flowchart

```
What is your G64 P?
  ├─ unset ─────────────→ tol_xyz = 0.01 (LinuxCNC default — good first guess)
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
