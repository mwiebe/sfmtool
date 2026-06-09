#!/bin/bash
# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0
#
# Global SfM from the purely-local background-floor matcher. Reuses the verified
# bgf.db that run_bgfloor.sh already built, runs COLMAP's global mapper (seeded),
# and compares to the workspace's baseline solve. Distinct bgf_global output so it
# never clobbers the incremental bgf_inc / bgf_incremental.sfmr results.
#
# Usage:  bash experiments/run_bgfloor_global.sh <workspace>
set -e
export PATH="$HOME/.pixi/bin:$PATH"

WS="$1"
ROOT="$(cd "$(dirname "$0")/.." && pwd)/$WS"
BASE=$(ls "$ROOT"/sfmr/*solve*.sfmr | head -1)
echo "============================== $WS (bgfloor global) =============================="

cd "$ROOT/.."
OUT="$ROOT/bgf_global"; rm -rf "$OUT" "$ROOT/sfmr/bgf_global.sfmr"
REC=$(pixi run python - "$ROOT" "$OUT" <<'PY' 2>&1 | grep -vE "WARN" | grep -E "^recon:"
import os, sys, pycolmap
root, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
pycolmap.set_random_seed(42)
o = pycolmap.GlobalPipelineOptions()
recs = pycolmap.global_mapping(f"{root}/bgf.db", root, out, o)
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
        -o "$ROOT/sfmr/bgf_global.sfmr" --tool-name colmap >/dev/null 2>&1
    verdict=$(pixi run sfm compare "$BASE" "$ROOT/sfmr/bgf_global.sfmr" 2>&1 \
        | grep -vE "WARN" | grep -oE "VERY SIMILAR|SIGNIFICANT DIFFERENCES" | head -1)
fi
printf "  bgfloor-global  reg=%s points=%s reproj=%s  %s\n" "$reg" "$pts" "$reproj" "$verdict"
