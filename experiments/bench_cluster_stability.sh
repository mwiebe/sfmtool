#!/bin/bash
# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0
#
# Stability sweep for the cluster matcher, built to chase down the seattle
# "no cluster solves" failure. Holds the matcher deterministic and varies, one
# axis at a time, the knobs that turned out to matter:
#
#   1. --max-features      feature budget per image (match + incremental solve)
#   2. solve seed          re-solve the SAME matches across seeds (is it random?)
#   3. --cluster-d         background rank
#   4. --cluster-alpha     radius multiplier
#   5. incremental vs global at sparse feature budgets
#
# Findings on seattle (default ~3072 feat/image): full registration and quality
# are stable across seed / d / alpha; the only failure mode is feature
# starvation, and it is cluster-specific — at 512 feat the cluster matcher drops
# an image while exhaustive keeps all 26. Incremental collapses below ~96
# feat/image (64 = no reconstruction at all); the global solver degrades far more
# gracefully (full registration down to ~128). alpha=1.0 is worse than the 0.8
# default (over-merged clusters get dropped).
#
# Usage:  pixi run -e test bash experiments/bench_cluster_stability.sh [dataset]
#         dataset in {seoul, seattle (default), kerry, dino}
set -u
export PATH="$HOME/.pixi/bin:$PATH"
cd "$(dirname "$0")/.."

DS="${1:-seattle}"
case "$DS" in
  seoul)   WS=seoul_bull_ws;       TOTAL=17; GLOBS="seoul_bull_ws/images/*.jpg" ;;
  seattle) WS=seattle_backyard_ws; TOTAL=26; GLOBS="seattle_backyard_ws/images/*.jpg" ;;
  kerry)   WS=kerry_park_ws;       TOTAL=48; GLOBS="kerry_park_ws/fisheye_left/*.jpg kerry_park_ws/fisheye_right/*.jpg" ;;
  dino)    WS=dino_dog_toy_ws;     TOTAL=85; GLOBS="dino_dog_toy_ws/images/*.jpg" ;;
  *) echo "unknown dataset: $DS" >&2; exit 1 ;;
esac
TMP="${TMPDIR:-/tmp}"

# match with the cluster matcher; extra args are passed straight to `sfm match`.
cluster_match () { # $1 out.matches ; rest = extra flags
  local out="$1"; shift
  pixi run sfm match --cluster "$@" $GLOBS -o "$out" >"$TMP/m.log" 2>&1
}

# solve a matches file and echo "reg/pts" (or FAIL / NOMATCH). $2 = i|g.
solve_report () { # $1 matches  $2 mode(i|g)
  local m="$1" mode="${2:-i}" flag="-i"
  [ "$mode" = "g" ] && flag="-g"
  [ -s "$m" ] || { echo "NOMATCH"; return; }
  pixi run sfm solve $flag "$m" -s 42 -o "$TMP/r.sfmr" >"$TMP/s.log" 2>&1
  local rc=$? reg pts
  reg=$(grep -oE "Images: [0-9]+" "$TMP/s.log" | grep -oE "[0-9]+" | tail -1)
  pts=$(grep -oE "Points: [0-9]+" "$TMP/s.log" | grep -oE "[0-9]+" | tail -1)
  { [ $rc -ne 0 ] || [ -z "$reg" ]; } && echo "SOLVE-FAIL" || echo "${reg}/${TOTAL} ${pts}p"
}

reproj_of () { # $1 sfmr
  pixi run sfm analyze "$1" --metrics 2>&1 \
    | grep -oE "Mean reprojection error: [0-9.]+" | grep -oE "[0-9.]+" | head -1
}

echo "###############  cluster stability sweep: $DS ($TOTAL images)  ###############"

echo
echo "== 1. --max-features sweep (cluster vs exhaustive, incremental seed 42) =="
printf "  %-9s %-14s %-14s\n" maxfeat cluster exhaustive
for mf in 64 96 128 256 512 1024 2048 0; do
  flag=""; tag="$mf"; [ "$mf" != "0" ] && flag="--max-features $mf" || tag="default"
  cluster_match "$TMP/sw_${DS}_c_${tag}.matches" $flag
  cres="$(solve_report "$TMP/sw_${DS}_c_${tag}.matches" i)"
  pixi run sfm match --exhaustive $flag $GLOBS -o "$TMP/sw_${DS}_e_${tag}.matches" >"$TMP/m.log" 2>&1
  eres="$(solve_report "$TMP/sw_${DS}_e_${tag}.matches" i)"
  printf "  %-9s %-14s %-14s\n" "$tag" "$cres" "$eres"
done

echo
echo "== 2. solve-seed sensitivity (same matches, default + a marginal budget) =="
cluster_match "$TMP/seed_def.matches"
cluster_match "$TMP/seed_512.matches" --max-features 512
for lab in def 512; do
  printf "  maxfeat=%-7s:" "$lab"
  for seed in 0 1 2 7 42 99 2024; do
    pixi run sfm solve -i "$TMP/seed_${lab}.matches" -s "$seed" -o "$TMP/r.sfmr" >"$TMP/s.log" 2>&1
    reg=$(grep -oE "Images: [0-9]+" "$TMP/s.log" | grep -oE "[0-9]+" | tail -1)
    pts=$(grep -oE "Points: [0-9]+" "$TMP/s.log" | grep -oE "[0-9]+" | tail -1)
    printf " s%s=%s/%sp" "$seed" "${reg:-FAIL}" "${pts:-?}"
  done; echo
done

echo
echo "== 3. --cluster-d sweep (default features, incremental seed 42) =="
for d in 4 6 10 20 28 40; do
  cluster_match "$TMP/d_$d.matches" --cluster-d "$d"
  printf "  d=%-3s %s\n" "$d" "$(solve_report "$TMP/d_$d.matches" i)"
done

echo
echo "== 4. --cluster-alpha sweep (default features, incremental seed 42) =="
for a in 0.5 0.6 0.7 0.8 0.9 1.0; do
  cluster_match "$TMP/a_$a.matches" --cluster-alpha "$a"
  printf "  alpha=%-4s %s\n" "$a" "$(solve_report "$TMP/a_$a.matches" i)"
done

echo
echo "== 5. incremental vs global at sparse feature budgets (cluster) =="
printf "  %-9s %-14s %-14s\n" maxfeat incremental global
for mf in 128 256 512; do
  cluster_match "$TMP/sp_$mf.matches" --max-features "$mf"
  printf "  %-9s %-14s %-14s\n" "$mf" \
    "$(solve_report "$TMP/sp_$mf.matches" i)" "$(solve_report "$TMP/sp_$mf.matches" g)"
done
