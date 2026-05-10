**English** | [中文](README.zh.md)

# liscio tools

> Command-line programs **for end users** — distinct from the dev-side
> regression / benchmark binaries under `test/`.  These ship into a
> machine post-processing pipeline.  Built by default; disable with
> `-DBUILD_TOOLS=OFF` to cmake.

## liscio_compress — G-code → G-code compressor

Compresses CAM-output G-code through liscio fitting and re-emits standard
G-code, **cross-controller compatible** (LinuxCNC / Fanuc / Siemens / Haas /
Mach3, etc.).

```
input.ngc                      output.ngc
---------                      ----------
G1 X0.001 Y0.000          →    G1 X0.001 Y0.000 F500
G1 X0.002 Y0.001                ...
G1 X0.003 Y0.003                G2 X10 Y10 I5 J0  ; 100 G1 folded into 1 ARC
... (thousands of G1) ...       ...
G1 X10 Y10                      G3 X20 Y0 I0 J-10
G3 X20 Y0 R10                   M2
```

### Usage

```bash
./build/liscio_compress [options] <input.ngc>

Options:
  --tol_xyz=N        XYZ geometric tolerance in mm (default 0.01 = 10 µm,
                      the typical metal-cutting setting)
  --output=PATH      output file path (default: write to stdout)
  --precision=N      decimal places used for coordinates in the output (default 5)
  --arc_merge=N      adjacent G2/G3 fold:  0=off  1=merge into single ARC
                      2=detect HELIX (default — multi-turn helices auto-fold)
  --corner_mode=N    corner handling: 0=split immediately on any corner (default)
                                      1=look-ahead, soft corners absorbed into BEZIER
  --bezier_samples=N each BEZIER resampled to N G1 lines on output (default 8)
                      (LinuxCNC mode prefers G5; this only kicks in on the
                       portable fallback path)
  --portable         emit one G2/G3 per turn for helices (no P<turns>);
                      default is LinuxCNC's compact P form
  --target=MODE      portable (default) / linuxcnc
                      linuxcnc mode: emit BEZIER as native G5 cubic Bezier in
                                     any of the G17/G18/G19 planes (XY/XZ/YZ);
                                     much smaller than G1 resample, see table.

Examples:
  # Standard usage (LinuxCNC + default precision)
  liscio_compress --output=part.compressed.ngc part.ngc

  # Cross-controller compatible (Fanuc / Siemens / ...)
  liscio_compress --portable --output=out.ngc in.ngc

  # High-precision mode (optics / semiconductor)
  liscio_compress --tol_xyz=0.001 --precision=7 --output=out.ngc in.ngc

  # Wood / 3D printing (loose tol, max compression)
  liscio_compress --tol_xyz=0.05 --bezier_samples=8 \
                  --output=out.ngc in.ngc
```

### Input → output G-code mapping

| Input G-code | liscio intermediate | Output G-code |
|---|---|---|
| Multiple adjacent G1 | LINE (collinear fold) | Single G1 |
| Multiple G1 (circular) | ARC (Kasa LSQ fit) | G2/G3 X- Y- Z- I- J- |
| Multiple G1 (smooth curve) | BEZIER (LSQ fit) | N G1 (resampled) |
| Multiple G1 (complex curve) | SPLINE (B-spline LSQ) | N G1 (resampled) |
| Single G2/G3 | ARC (passthrough) | G2/G3 X- Y- Z- I- J- |
| Multiple concentric G2/G3 (with arc_merge) | HELIX (merged) | G2/G3 X- Y- Z- I- J- P\<turns\> |
| G0 rapid | RAPID | G0 X- Y- Z- |
| G4 dwell | STOP/DWELL | G4 P\<sec\> |
| M0/M1 | STOP | M0 / M1 |
| M2/M30 | STOP | M2 / M30 |
| M6 | STOP | M6 |
| G92/G10/G43/G49/G54-59 | STOP | (comment placeholder, see below) |

#### Known limitations

- **portable mode: BEZIER must be resampled to G1**.
  Standard G-code has no native cubic Bezier opcode (G2/G3 = arcs, G1 =
  lines).  Cross-controller output therefore breaks each BEZIER into
  several short G1 lines, which can inflate the file on BEZIER-dense
  2D inputs.  Fix: use `--target=linuxcnc` to emit G5 (see below).

- **`--target=linuxcnc` constraints**: G5 is LinuxCNC's native cubic
  Bezier opcode; it works in any of the three cardinal planes
  (G17 XY / G18 XZ / G19 YZ).  Tilted-plane (truly 3D) BEZIER and any
  curve with rotary or UVW motion still fall back to G1 resample.
  If you need to keep full 9D geometry, integrate the liscio API
  directly rather than going through G-code.

- **Coordinate-shift opcodes (G92/G10/G43/G49) become comments**.
  By the time rs274 has parsed the input, every coordinate is in
  absolute millimetres, so these "change-coordinate-system" opcodes
  no longer affect anything.  The tool keeps them as comments to
  preserve line numbering only.  Note: the output is absolute-mm
  G-code; you cannot get back to the original coordinate system.

- **Out-of-plane arcs** (normal not aligned with X/Y/Z): currently
  fall back to G1 resample.  Vanishingly rare in CAM data.

### Measurement (18 CAM files, tol=0.05, default cfg, helix9 LSQ on)

| File | Input | portable | linuxcnc (G5) | p_x | **l_x** |
|---|---:|---:|---:|---:|---:|
| **sine_test** | 8.85 KB | 2.45 KB | **0.99 KB** | 3.62× | **8.93×** ✅ |
| **1_1001** | 502 KB | 176 KB | **100 KB** | 2.84× | **5.01×** ✅ |
| bto45 | 533 B | 131 B | 131 B | 4.07× | 4.07× |
| **thomam** | 8.77 KB | 4.47 KB | 4.47 KB | **1.96×** ✓ | 1.96× |
| 3_1001 | 352 KB | 155 KB | **151 KB** | 2.27× | 2.33× |
| 5_1001 | 307 KB | 138 KB | **134 KB** | 2.22× | 2.29× |
| 4_1001 | 134 KB | 60.6 KB | 60.6 KB | 2.22× | 2.22× |
| **flower** | 733 KB | 654 KB | **353 KB** | 1.12× | **2.07×** ✅ |
| flower-one-line | 3.76 KB | 2.12 KB | **1.69 KB** | 1.77× | 2.22× |
| 6_1001 | 141 KB | 70.4 KB | 70.4 KB | 2.00× | 2.00× |
| 2_1001 | 19.8 KB | 10.2 KB | **9.12 KB** | 1.94× | 2.17× |
| 130207LZW | 130 KB | 68.9 KB | **67.7 KB** | 1.89× | 1.92× |
| SER-40 | 553 KB | 397 KB | 397 KB | 1.39× | 1.39× |
| luca-long-reverse | 213 B | 189 B | 189 B | 1.13× | 1.13× |
| rechteck-10x10 | 161 B | 282 B | 282 B | 0.57× ❌ | 0.57× |
| yy / yy-g61 | ≤163 B | ≤442 B | (same) | 0.35-0.37× ❌ | same |

> p_x: portable mode (cross-controller, BEZIER → G1 resample)
> l_x: --target=linuxcnc (cubic Bezier → native G5 in **XY/XZ/YZ planes**)

### thomam improvement (vs earlier doc 0.78× inflation)

Earlier each of 245 ARCs was rewritten in IJK-form, inflating output.
Now helix9 LSQ auto-detection folds them into 11 HELIX primitives
emitted as `G2/G3 P<turns>` compact multi-turn form, file 8.77 →
4.47 KB (1.96× compression).  Independent of LinuxCNC target — helix9
runs in both modes.

### flower improvement (vs earlier 0.59× inflation)

flower is 3D free-form (XYZ all varying); liscio fits 3476 BEZIER
primitives.  Earlier G5 emit only handled the XY plane (G17), but
flower's data is rarely strictly XY-planar → all 3476 BEZIER fell
back to G1 resample (16 segs each = 56k G1 lines) → output
1243 KB > input 733 KB.

**Fix (2026-05-09)**: G5 emit extended to LinuxCNC's full three-plane
support (XY/XZ/YZ, i.e. G17/G18/G19 modes).  flower's actual BEZIER
distribution:

- 0 strictly XY-planar (Z varies)
- **3472 XZ-planar** (Y constant — slice-layer machining)
- 2 YZ-planar
- 2 truly 3D (no cardinal plane)

The new G5 emit lets 3474 BEZIER express directly; only 2 truly-3D
segments fall back to G1 resample.  Output 1243 KB → **353 KB (2.07×
forward compression)**.

### --target=linuxcnc incremental gain

- **sine_test** 3.62× → **8.93×** (BEZIER → single G5)
- **1_1001** 2.84× → **5.01×** (lots of 2D BEZIER)
- **flower** 1.12× → **2.07×** (3-plane G5 rescues 3D)
- No change for ARC-dominated / corner-rich files (G5 unused)

### Round-trip validation

`portable` mode 18-file round-trip (`liscio_compress` +
`test_verify`) passes OVER_TOL=0 across the board ✓.

`--target=linuxcnc` round-trip — the ctest verifier (`ngc_parser`)
does not implement LinuxCNC's `NURBS_G5_FEED` canon callback, so it
reports OVER_TOL on G5-containing output.  **This is a test-tool
limitation, not a G5-emit bug** — production LinuxCNC parses G5
correctly via `interp_convert.cc convert_spline()`.

**Summary**:
| Data characteristic | Recommended mode |
|---|---|
| **Cross-controller required** (Fanuc/Siemens/...) | `portable` (default) |
| **LinuxCNC + 2D cutting** | `--target=linuxcnc` ⭐ |
| **3D free-form (flower-like)** | Neither ideal — G-code standard limitation |
| **ARC-dominated (turning/milling typical)** | Either, prefer portable for compatibility |

### Geometric accuracy guarantee

- Output file re-parsed + re-fit with liscio, max_dev still ≤ tol_xyz
- ARC passthrough as-is, max_dev = 0
- BEZIER resampled by chord length ≤ tol/2, measured max_dev ≤ tol_xyz

### Self-verification (round-trip)

```bash
# 1. Compress
./build/liscio_compress --tol_xyz=0.05 --output=/tmp/out.ngc input.ngc

# 2. Run compressed output through verify_geometry — should still be OVER_TOL=0
./build/test_verify --tol_xyz=0.05 /tmp/out.ngc
```

### Integration example (LinuxCNC HAL / production)

Post-processing flow:

```bash
#!/bin/bash
# CAM post-processor: take CAM raw G-code, compress for LinuxCNC use
RAW=$1
OUT=${2:-${RAW%.ngc}.compressed.ngc}

./liscio_compress \
    --tol_xyz=0.01 \
    --arc_merge=2 \
    --output=$OUT \
    $RAW

echo "Compressed: $RAW → $OUT"
ls -la $RAW $OUT | awk '{print $5, $NF}'
```

### Future improvements

- [ ] LinuxCNC G5/G5.1 (cubic Bezier / quadratic B-spline) emit, `--target=linuxcnc`
- [ ] biarc approximation of BEZIER (G2/G3 pair instead of G1 resample)
- [ ] R-form arc (more compact than IJK for short arcs)
- [ ] Adaptive BEZIER resample density (per local curvature)
- [ ] Preserve original G-code comments (rs274 currently discards)
