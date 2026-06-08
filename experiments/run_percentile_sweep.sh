#!/bin/bash
# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0
#
# Full incremental reconstruction at a variety of global radii, each derived as a
# percentile of the per-point "first-excluded" distance d_(r+1) (the neighbour
# just past each point's cliff; see exp14). Plus the current Otsu x1.25 baseline
# for reference. For each radius: build cluster matches, verify (seeded), run
# incremental SfM (seeded), and compare to the workspace's baseline solve.
#
# Usage:  bash experiments/run_percentile_sweep.sh <workspace> [pcts]
#         pcts default "10,20,35,50"
set -e
export PATH="$HOME/.pixi/bin:$PATH"

WS="$1"
PCTS="${2:-10,20,35,50}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)/$WS"
BASE=$(ls "$ROOT"/sfmr/*solve*.sfmr | head -1)
cd "$(dirname "$0")"

echo "============================== $WS =============================="
echo "baseline=$(basename "$BASE")  percentiles=$PCTS"

# Percentile -> absolute radius T, computed once from d_(r+1) over all points.
mapfile -t SPECS < <(pixi run -e experiments python - "$ROOT" "$PCTS" <<'PY' 2>/dev/null
import sys, glob, numpy as np
from exp03_radius_clusters import K
from exp05_cluster_match import derive_threshold
from sfm_descriptors import load_descriptor_bank
from sfmtool import KdForest
root, pcts = sys.argv[1], [float(x) for x in sys.argv[2].split(",")]
bank = load_descriptor_bank(sorted(glob.glob(f"{root}/sfmr/*solve*.sfmr"))[0])
desc = np.ascontiguousarray(bank.descriptors)
_, dst = KdForest(desc, preset="accurate").query(desc, k=K)
dn = dst[:, 1:8]
g = (dn[:, 1:] / np.maximum(dn[:, :-1], 1e-6)).argmax(1)
outer = dst[np.arange(len(dst)), 2 + g]
print(f"otsu125 {derive_threshold(dst[:, 1], 'otsu') * 1.25:.1f}")
for p in pcts:
    print(f"p{int(p)} {np.percentile(outer, p):.1f}")
PY
)

printf "  %-8s %6s %5s %7s %8s  %s\n" tag T reg points reproj compare
for spec in "${SPECS[@]}"; do
    tag="${spec%% *}"; T="${spec##* }"
    M="$ROOT/matches/pct_${tag}.matches"
    pixi run -e experiments python -u exp05_cluster_match.py \
        "$ROOT/sfmr/"'*solve*.sfmr' --out "$M" --mode clusters \
        --threshold "$T" --t-scale 1.0 2>&1 | grep -vE "WARN" >/dev/null

    cd "$ROOT/.."
    rm -f "$ROOT/pct.db"
    pixi run sfm to-colmap-db "$M" --out-db "$ROOT/pct.db" 2>&1 | grep -vE "WARN" >/dev/null
    pixi run python - "$ROOT" <<'PY' 2>&1 | grep -vE "WARN" >/dev/null
import sys, pycolmap
pycolmap.set_random_seed(42)
pycolmap.geometric_verification(f"{sys.argv[1]}/pct.db")
PY
    OUT="$ROOT/pct_inc"; rm -rf "$OUT"
    REC=$(pixi run python - "$ROOT" "$OUT" <<'PY' 2>&1 | grep -vE "WARN" | grep -E "^recon:"
import os, sys, pycolmap
root, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
pycolmap.set_random_seed(42)
o = pycolmap.IncrementalPipelineOptions(); o.random_seed = 42
recs = pycolmap.incremental_mapping(f"{root}/pct.db", root, out, o)
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
            -o "$ROOT/sfmr/pct_${tag}.sfmr" --tool-name colmap >/dev/null 2>&1
        verdict=$(pixi run sfm compare "$BASE" "$ROOT/sfmr/pct_${tag}.sfmr" 2>&1 \
            | grep -vE "WARN" | grep -oE "VERY SIMILAR|SIGNIFICANT DIFFERENCES" | head -1)
    fi
    printf "  %-8s %6s %5s %7s %8s  %s\n" "$tag" "$T" "$reg" "$pts" "$reproj" "$verdict"
    cd "$(dirname "$0")"
done
