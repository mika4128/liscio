**English** | [中文](README.zh.md)

# liscio tools

> User-facing tools (`BUILD_TOOLS=ON`). Different from `test/` regression /
> benchmarks; these are for production.

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
  --tol_xyz=N        XYZ fit tolerance mm (default 0.01 = 10µm sweet spot)
  --output=PATH      output file (default stdout)
  --precision=N      G-code numeric decimals (default 5)
  --arc_merge=N      G2/G3 merge: 0=off 1=ARC 2=HELIX (default 2)
  --corner_mode=N    Corner detection: 0=IMMEDIATE 1=LOOKAHEAD (default 0)
  --bezier_samples=N BEZIER/SPLINE resampled to N G1 (default 16)
  --portable         Split helical multi-turn into multiple single-turn G2/G3
                      (no P<turns>); default uses LinuxCNC P form (more compact)
  --target=MODE      portable (default) / linuxcnc
                      linuxcnc mode: XY-planar BEZIER emitted as G5 cubic Bezier
                                     (significantly reduces BEZIER-dense files; see table)

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

- **portable mode BEZIER → G1 resample**: standard G-code has no native Bezier.
  For BEZIER-dense 2D inputs, output inflates. Solution:
  → `--target=linuxcnc` enables G5 emit (XY-plane Bezier directly G5 output;
    measured sine file 1.83× → **8.93×**).
- **`--target=linuxcnc` only for XY-plane Bezier**: G5 G-code standard supports
  only X/Y axes; 3D Bezier with Z variation (e.g. flower) still falls back to
  G1 resampling. 9D rotary/UVW also fall back. For full 9D, the better path is
  to integrate liscio API directly (bypassing G-code).
- **G92/G10/G43/G49 offset events**: currently emitted as comment placeholders,
  no coordinate-system restoration (because ngc_parser already converted all
  coordinates to absolute mm via rs274). Note: compressed NGC is in absolute
  mm form, cannot return to original coordinate system.
- **Out-of-plane arc (tilted plane)**: currently falls back to G1 resampling.
  Rare.

### Measurement (16 CAM files, tol=0.05) — portable vs linuxcnc

| File | Input | portable | linuxcnc (G5) | p_x | **l_x** |
|---|---:|---:|---:|---:|---:|
| **1_1001** | 502 KB | 308 KB | **107 KB** | 1.63× | **4.68×** ✅ |
| 2_1001 | 19.8 KB | 11.9 KB | **9.06 KB** | 1.66× | 2.18× |
| 3_1001 | 352 KB | 188 KB | **177 KB** | 1.87× | 1.99× |
| 4_1001 | 134 KB | 71.6 KB | 71.6 KB | 1.88× | 1.88× |
| 5_1001 | 307 KB | 178 KB | **168 KB** | 1.72× | 1.82× |
| 6_1001 | 141 KB | 82.1 KB | 82.1 KB | 1.71× | 1.71× |
| bto45 | 533 B | 131 B | 131 B | 4.07× | 4.07× |
| **sine_test** | 8.85 KB | 4.83 KB | **0.99 KB** | 1.83× | **8.93×** ✅ |
| **flower** | 733 KB | 1243 KB | 1243 KB | 0.59× | 0.59× ❌ (3D) |
| thomam | 8.77 KB | 11.3 KB | 11.3 KB | 0.78× | 0.78× |
| rechteck/yy/luca | ≤213 B | ≤442 B | ≤442 B | <1× | <1× |

**linuxcnc mode incremental gain**:
- sine_test 1.83× → **8.93×** (5× more compression)
- 1_1001 1.63× → **4.68×** (3× more compression)
- 2D BEZIER-dense files (sine, 1_1001, 2_1001, ...) benefit greatly
- 3D BEZIER (flower, 4_1001, 6_1001) no benefit — G5 is XY-plane only
- ARC-dominated files (thomam) unchanged

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
