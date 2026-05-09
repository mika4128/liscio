#!/usr/bin/env bash
# sweep_tol.sh — run liscio test_metrics across tolerance × arc-mode grid.
#
# Output layout (relative to project root):
#   results/sweep/
#   ├── metrics_<mode>_<tol>.csv         per-config raw metrics (16 files)
#   ├── summary.csv                      aggregate (mode, tol) → totals
#   └── images/<mode>_<tol>/<file>.png   PNG for representative configs
#
# Usage:  bash test/sweep_tol.sh
set -euo pipefail
cd "$(dirname "$0")/.."

OUT=results/sweep
mkdir -p "$OUT/images"

TOLS=(0.05 0.01 0.001 0.0001 0.00001 0.000001 0.0000001 0.00000001)
MODES=(0 8)                         # arc_subseg: 0=passthrough, 8=split
PNG_TOLS=(0.05 0.001 0.000001)      # subset for PNG to avoid bloat
PNG_FILES=(sine_test.ngc flower.ngc thomam.ngc 1_1001.ngc)

shopt -s nullglob
NGC=(test/data/*.ngc test/data/*.NGC)

# ---------- metrics sweep ----------
echo "mode,tol,total_input,total_output,total_ratio,worst_max_dev,worst_len_err_pct,total_elapsed_ms" \
    > "$OUT/summary.csv"

for mode in "${MODES[@]}"; do
    for tol in "${TOLS[@]}"; do
        csv="$OUT/metrics_arc${mode}_tol${tol}.csv"
        ./build/test_metrics --tol_xyz="$tol" --arc_subseg="$mode" \
            --csv="$csv" "${NGC[@]}" > /dev/null 2>&1 || \
            { echo "WARN: failed mode=$mode tol=$tol"; continue; }

        # Aggregate from CSV (skip header).
        awk -F, -v mode="$mode" -v tol="$tol" '
            NR>1 {
                ti += $2; to += $3
                if ($8+0 > md) md = $8+0
                if (($7+0 < 0 ? -$7 : $7) > le) le = ($7+0 < 0 ? -$7 : $7)
                te += $9
            }
            END {
                ratio = (to>0) ? ti/to : 0
                printf "%d,%s,%d,%d,%.3f,%.6g,%.4f,%.3f\n",
                       mode, tol, ti, to, ratio, md, le, te/1000.0
            }' "$csv" >> "$OUT/summary.csv"
        echo "  done: arc_subseg=$mode tol_xyz=$tol"
    done
done

echo
echo "## Sweep summary"
column -t -s, "$OUT/summary.csv"
echo

# ---------- representative PNG generation ----------
echo "Rendering representative PNGs..."
for mode in "${MODES[@]}"; do
    for tol in "${PNG_TOLS[@]}"; do
        outdir="$OUT/images/arc${mode}_tol${tol}"
        mkdir -p "$outdir"
        for f in "${PNG_FILES[@]}"; do
            base="$outdir/${f}"
            (cd "$outdir" && cp "../../../../test/data/$f" . && \
                "../../../../build/test_dump_csv" --tol_xyz="$tol" \
                    --arc_subseg="$mode" --samples=30 "$f" > /dev/null && \
                python3 ../../../../test/plot_fit.py --3d "$f" > /dev/null 2>&1) || \
                { echo "  PNG fail mode=$mode tol=$tol file=$f"; continue; }
        done
        echo "  PNGs ready: arc${mode}_tol${tol}"
    done
done

echo "Done. See $OUT/"
