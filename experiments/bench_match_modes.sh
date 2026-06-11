#!/bin/bash
# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0
#
# Compare the matching modes (COLMAP exhaustive / sequential vs the in-house
# cluster matcher) end to end on every bootstrapped workspace: time `sfm match`,
# solve from the resulting .matches with a fixed seed, then report solve quality
# (registered images, points, mean reprojection error) and the `sfm compare`
# verdict against the workspace's stored baseline solve.
#
# All three modes run through the same `sfm match` harness (DB populate -> match
# -> geometric verification -> .matches), so the times are directly comparable.
# Cluster uses the production defaults (d=10, alpha=0.8). CPU only (pycolmap
# without CUDA); the cluster matcher is CPU Rust, so this is a CPU-to-CPU read.
#
# Usage:  pixi run -e test bash experiments/bench_match_modes.sh [dataset ...]
#         (datasets default to all four; names: seoul seattle kerry dino)
set -u
export PATH="$HOME/.pixi/bin:$PATH"
cd "$(dirname "$0")/.."

# dataset -> "workspace:total_images:image_globs"
spec_for () {
  case "$1" in
    seoul)   echo "seoul_bull_ws:17:seoul_bull_ws/images/*.jpg" ;;
    seattle) echo "seattle_backyard_ws:26:seattle_backyard_ws/images/*.jpg" ;;
    kerry)   echo "kerry_park_ws:48:kerry_park_ws/fisheye_left/*.jpg kerry_park_ws/fisheye_right/*.jpg" ;;
    dino)    echo "dino_dog_toy_ws:85:dino_dog_toy_ws/images/*.jpg" ;;
    *) echo "" ;;
  esac
}

DATASETS=("${@:-}")
[ -z "${DATASETS[*]}" ] && DATASETS=(seoul seattle kerry dino)

TMP="${TMPDIR:-/tmp}"
printf "%-9s %-11s %8s  %-7s %-8s %-9s  %s\n" \
  dataset mode match reg points reproj "verdict vs baseline"
echo "---------------------------------------------------------------------------------"

for ds in "${DATASETS[@]}"; do
  spec="$(spec_for "$ds")"; [ -z "$spec" ] && { echo "unknown dataset: $ds" >&2; continue; }
  ws="${spec%%:*}"; rest="${spec#*:}"; total="${rest%%:*}"; globs="${rest#*:}"
  base="$(ls "$ws"/sfmr/*solve*.sfmr 2>/dev/null | head -1)"

  for mode in exhaustive sequential cluster; do
    mf="$TMP/bench_${ds}_${mode}.matches"
    sf="$TMP/solve_${ds}_${mode}.sfmr"

    t0=$(date +%s.%N)
    pixi run sfm match "--$mode" $globs -o "$mf" >"$TMP/m.log" 2>&1
    if [ $? -ne 0 ]; then printf "%-9s %-11s   FAILED (match)\n" "$ds" "$mode"; continue; fi
    t1=$(date +%s.%N)
    mt=$(echo "$t1 - $t0" | bc)

    pixi run sfm solve -i "$mf" -s 42 -o "$sf" >"$TMP/s.log" 2>&1
    if [ $? -ne 0 ]; then printf "%-9s %-11s %7.1fs FAILED (solve)\n" "$ds" "$mode" "$mt"; continue; fi
    reg=$(grep -oE "Images: [0-9]+" "$TMP/s.log" | grep -oE "[0-9]+" | tail -1)
    pts=$(grep -oE "Points: [0-9]+" "$TMP/s.log" | grep -oE "[0-9]+" | tail -1)

    reproj=$(pixi run sfm analyze "$sf" --metrics 2>&1 \
      | grep -oE "Mean reprojection error: [0-9.]+" | grep -oE "[0-9.]+" | head -1)

    verdict="(no baseline)"
    if [ -n "$base" ]; then
      verdict=$(pixi run sfm compare "$base" "$sf" 2>&1 \
        | grep -oE "IDENTICAL|VERY SIMILAR|SIGNIFICANT DIFFERENCES" | head -1)
    fi

    printf "%-9s %-11s %7.1fs  %s/%-3s %-8s %-7s px  %s\n" \
      "$ds" "$mode" "$mt" "${reg:-?}" "$total" "${pts:-?}" "${reproj:-?}" "${verdict:-?}"
  done
done
