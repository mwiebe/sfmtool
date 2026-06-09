#!/bin/bash
# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0
#
# End-to-end reconstruction from the purely-local background-floor matcher
# (exp05 --mode bgfloor): build matches, verify (seeded), incremental SfM
# (seeded), compare to the workspace's baseline solve. Distinct bgf* output
# names so it never clobbers the nbr_* / pct_* solves.
#
# Usage:  bash experiments/run_bgfloor.sh <workspace> [alpha] [b0]
set -e
export PATH="$HOME/.pixi/bin:$PATH"

WS="$1"
ALPHA="${2:-0.8}"
B0="${3:-8}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)/$WS"
BASE=$(ls "$ROOT"/sfmr/*solve*.sfmr | head -1)
echo "============================== $WS (bgfloor a=$ALPHA b0=$B0) =============================="

cd "$(dirname "$0")"
pixi run -e experiments python -u exp05_cluster_match.py "$ROOT/sfmr/"'*solve*.sfmr' \
    --out "$ROOT/matches/bgf.matches" --mode bgfloor --bg-alpha "$ALPHA" --bg-b0 "$B0" \
    2>&1 | grep -vE "WARN" | grep -E "after one"

cd "$ROOT/.."
rm -f "$ROOT/bgf.db"
pixi run sfm to-colmap-db "$ROOT/matches/bgf.matches" --out-db "$ROOT/bgf.db" 2>&1 \
    | grep -vE "WARN" >/dev/null
pixi run python - "$ROOT" <<'PY' 2>&1 | grep -vE "WARN" >/dev/null
import sys, pycolmap
pycolmap.set_random_seed(42)
pycolmap.geometric_verification(f"{sys.argv[1]}/bgf.db")
PY

OUT="$ROOT/bgf_inc"; rm -rf "$OUT"
REC=$(pixi run python - "$ROOT" "$OUT" <<'PY' 2>&1 | grep -vE "WARN" | grep -E "^recon:"
import os, sys, pycolmap
root, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
pycolmap.set_random_seed(42)
o = pycolmap.IncrementalPipelineOptions(); o.random_seed = 42
recs = pycolmap.incremental_mapping(f"{root}/bgf.db", root, out, o)
best = max(recs.values(), key=lambda r: r.num_reg_images()) if recs else None
if best is None:
    print("recon: FAILED 0 0 0")
else:
    os.makedirs(f"{out}/best", exist_ok=True); best.write(f"{out}/best")
    print(f"recon: OK {best.num_reg_images()} {best.num_points3D()} "
          f"{best.compute_mean_reprojection_error():.2f}")
PY
)
read -r _ ok reg pts reproj <<<"$REC"
verdict="-"
if [ "$ok" = "OK" ]; then
    pixi run sfm from-colmap-bin "$OUT/best" --image-dir "$ROOT" \
        -o "$ROOT/sfmr/bgf_incremental.sfmr" --tool-name colmap >/dev/null 2>&1
    verdict=$(pixi run sfm compare "$BASE" "$ROOT/sfmr/bgf_incremental.sfmr" 2>&1 \
        | grep -vE "WARN" | grep -oE "VERY SIMILAR|SIGNIFICANT DIFFERENCES" | head -1)
fi
printf "  bgfloor  reg=%s points=%s reproj=%s  %s\n" "$reg" "$pts" "$reproj" "$verdict"
