#!/usr/bin/env bash
# Copy build artifacts into results/.
# Usage: sync_results.sh BUILD_DIR SRC_DIR
set -e
BUILD="$1"
SRC="$2"
mkdir -p "$SRC/results/images" "$SRC/results/data"
shopt -s nullglob
png=("$BUILD/ngc_out/"*.png)
if [ ${#png[@]} -gt 0 ]; then
    cp "${png[@]}" "$SRC/results/images/"
fi
[ -f "$BUILD/metrics.csv" ]    && cp "$BUILD/metrics.csv"    "$SRC/results/data/"
[ -f "$BUILD/metrics_g1.csv" ] && cp "$BUILD/metrics_g1.csv" "$SRC/results/data/"
echo "synced ${#png[@]} png + 2 csv → results/"
