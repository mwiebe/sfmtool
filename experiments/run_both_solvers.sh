#!/bin/bash
# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0
#
# Run BOTH incremental and global SfM on the cluster matcher's matches.
# Usage: run_both_solvers.sh <workspace> [mode]   (mode default: clusters)
# Builds + verifies the COLMAP db once, then runs each mapper on it and
# compares the result to the workspace's baseline solve.
#
# Usage:  bash experiments/run_both_solvers.sh <workspace_dir_name>
set -e
export PATH="$HOME/.pixi/bin:$PATH"

WS="$1"
MODE="${2:-clusters}"
REFINE="${3:-0}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)/$WS"
BASE=$(ls "$ROOT"/sfmr/*solve*.sfmr | head -1)
echo "============================== $WS (mode=$MODE refine=$REFINE) =============================="
echo "baseline=$(basename "$BASE")"

# (Re)generate the cluster matches.
cd "$(dirname "$0")"
pixi run -e experiments python -u exp05_cluster_match.py "$ROOT/sfmr/"'*solve*.sfmr' \
    --out "$ROOT/matches/nbr.matches" --mode "$MODE" --refine "$REFINE" 2>&1 \
    | grep -vE "WARN" | grep -E "after one|wrote"

# Build + verify the database once.
cd "$ROOT/.."
rm -f "$ROOT/nbr.db"
pixi run sfm to-colmap-db "$ROOT/matches/nbr.matches" --out-db "$ROOT/nbr.db" 2>&1 | grep -vE "WARN" | tail -1
pixi run python - "$ROOT" <<'PY' 2>&1 | grep -vE "WARN" | grep -E "verify:"
import sys, sqlite3, pycolmap
db = f"{sys.argv[1]}/nbr.db"
pycolmap.set_random_seed(42)  # RANSAC verification is randomized; pin it
pycolmap.geometric_verification(db)
con = sqlite3.connect(db)
nm = con.execute("SELECT COALESCE(SUM(rows),0) FROM matches").fetchone()[0]
ni = con.execute("SELECT COALESCE(SUM(rows),0) FROM two_view_geometries").fetchone()[0]
nvp = con.execute("SELECT COUNT(*) FROM two_view_geometries WHERE rows>0").fetchone()[0]
con.close()
print(f"verify: in_matches={nm} inliers={ni} verified_pairs={nvp} "
      f"inlier_rate={ni / max(nm, 1):.2f}")
PY

for SOLVER in incremental global; do
    OUT="$ROOT/nbr_${SOLVER}"
    rm -rf "$OUT" "$ROOT/sfmr/nbr_${SOLVER}.sfmr"
    pixi run python - "$ROOT" "$SOLVER" "$OUT" <<'PY' 2>&1 | grep -vE "WARN" | grep -E "recon:"
import os, sys, pycolmap
root, solver, out = sys.argv[1], sys.argv[2], sys.argv[3]
db = f"{root}/nbr.db"
os.makedirs(out, exist_ok=True)
pycolmap.set_random_seed(42)
if solver == "incremental":
    o = pycolmap.IncrementalPipelineOptions(); o.random_seed = 42
    recs = pycolmap.incremental_mapping(db, root, out, o)
else:
    o = pycolmap.GlobalPipelineOptions()
    recs = pycolmap.global_mapping(db, root, out, o)
best = max(recs.values(), key=lambda r: r.num_reg_images()) if recs else None
if best is None:
    print(f"recon: {solver}: FAILED (no model)")
else:
    os.makedirs(f"{out}/best", exist_ok=True); best.write(f"{out}/best")
    print(f"recon: {solver:>11}: reg={best.num_reg_images()} "
          f"points={best.num_points3D()} obs={best.compute_num_observations()} "
          f"reproj={best.compute_mean_reprojection_error():.3f}px")
PY
    if [ -d "$OUT/best" ]; then
        pixi run sfm from-colmap-bin "$OUT/best" --image-dir "$ROOT" \
            -o "$ROOT/sfmr/nbr_${SOLVER}.sfmr" --tool-name colmap 2>&1 >/dev/null
        echo -n "  compare ($SOLVER): "
        pixi run sfm compare "$BASE" "$ROOT/sfmr/nbr_${SOLVER}.sfmr" 2>&1 | grep -vE "WARN" | \
            grep -E "VERY SIMILAR|SIGNIFICANT DIFFERENCES|Scale:|Mean:" | tr '\n' ' '; echo
    fi
done
