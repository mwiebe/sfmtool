#!/bin/bash
# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0
#
# Global SfM from the background-floor matcher. Reuses the verified per-mode db
# that run_bgfloor.sh already built (bgf = edges, bgc = materialized clusters),
# runs COLMAP's global mapper (seeded), and compares to the workspace's baseline
# solve. Distinct <pfx>_global output so it never clobbers the incremental
# <pfx>_inc / <pfx>_incremental.sfmr results.
#
# Usage:  bash experiments/run_bgfloor_global.sh <workspace> [pfx]   (pfx: bgf|bgc)
set -e
export PATH="$HOME/.pixi/bin:$PATH"

WS="$1"
PFX="${2:-bgf}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)/$WS"
BASE=$(ls "$ROOT"/sfmr/*solve*.sfmr | head -1)
echo "============================== $WS ($PFX global) =============================="

cd "$ROOT/.."
OUT="$ROOT/${PFX}_global"; rm -rf "$OUT" "$ROOT/sfmr/${PFX}_global.sfmr"
REC=$(pixi run python - "$ROOT" "$PFX" "$OUT" <<'PY' 2>&1 | grep -vE "WARN" | grep -E "^recon:"
import os, sys, pycolmap
root, pfx, out = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(out, exist_ok=True)
pycolmap.set_random_seed(42)
o = pycolmap.GlobalPipelineOptions()
recs = pycolmap.global_mapping(f"{root}/{pfx}.db", root, out, o)
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
        -o "$ROOT/sfmr/${PFX}_global.sfmr" --tool-name colmap >/dev/null 2>&1
    verdict=$(pixi run sfm compare "$BASE" "$ROOT/sfmr/${PFX}_global.sfmr" 2>&1 \
        | grep -vE "WARN" | grep -oE "VERY SIMILAR|SIGNIFICANT DIFFERENCES" | head -1)
fi
printf "  %s-global  reg=%s points=%s reproj=%s  %s\n" "$PFX" "$reg" "$pts" "$reproj" "$verdict"
