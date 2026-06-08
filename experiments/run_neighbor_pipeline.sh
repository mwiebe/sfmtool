#!/bin/bash
# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0
#
# End-to-end driver for the neighbours-mode cluster matcher (exp05):
#   .matches -> COLMAP db -> geometric verification -> incremental mapping
#   -> .sfmr -> `sfm compare` against the workspace's baseline solve.
#
# Usage:  bash experiments/run_neighbor_pipeline.sh <workspace_dir_name>
#   e.g.  bash experiments/run_neighbor_pipeline.sh seattle_backyard_ws
set -e
export PATH="$HOME/.pixi/bin:$PATH"

WS="$1"
MODE="${2:-neighbors}"
MINSIZE="${3:-2}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)/$WS"
BASE=$(ls "$ROOT"/sfmr/*solve*.sfmr | head -1)
echo "=== $WS (mode=$MODE min_size=$MINSIZE) ===  baseline=$(basename "$BASE")"

cd "$(dirname "$0")"
pixi run -e experiments python -u exp05_cluster_match.py "$ROOT/sfmr/"'*solve*.sfmr' \
    --out "$ROOT/matches/nbr.matches" --mode "$MODE" \
    --min-cluster-size "$MINSIZE" 2>&1 | grep -vE "WARN"

cd "$ROOT/.."
rm -f "$ROOT/nbr.db" "$ROOT/sfmr/nbr_match.sfmr"; rm -rf "$ROOT/nbr_recon"
pixi run sfm to-colmap-db "$ROOT/matches/nbr.matches" --out-db "$ROOT/nbr.db" 2>&1 | grep -vE "WARN" | tail -1

pixi run python - "$ROOT" <<'PY' 2>&1 | grep -vE "WARN" | grep -E "recon:|verify:" | tail -3
import os, sys, sqlite3, pycolmap
root = sys.argv[1]; db = f"{root}/nbr.db"; out = f"{root}/nbr_recon"
os.makedirs(out, exist_ok=True)
pycolmap.set_random_seed(42)                 # RANSAC verification is randomized
pycolmap.geometric_verification(db)          # step 4: per-pair RANSAC
# read verification stats straight from the COLMAP SQLite db
con = sqlite3.connect(db)
nm = con.execute("SELECT COALESCE(SUM(rows),0) FROM matches").fetchone()[0]
npair = con.execute("SELECT COUNT(*) FROM matches WHERE rows>0").fetchone()[0]
ni = con.execute("SELECT COALESCE(SUM(rows),0) FROM two_view_geometries").fetchone()[0]
nvp = con.execute("SELECT COUNT(*) FROM two_view_geometries WHERE rows>0").fetchone()[0]
con.close()
print(f"verify: in_matches={nm} in_pairs={npair} inliers={ni} "
      f"verified_pairs={nvp} inlier_rate={ni / max(nm, 1):.2f}")
pycolmap.set_random_seed(42)
o = pycolmap.IncrementalPipelineOptions(); o.random_seed = 42
recs = pycolmap.incremental_mapping(db, root, out, o)
best = max(recs.values(), key=lambda r: r.num_reg_images()) if recs else None
if best is not None:
    os.makedirs(f"{out}/best", exist_ok=True)
    best.write(f"{out}/best")
    print(f"recon: reg={best.num_reg_images()} points={best.num_points3D()} "
          f"obs={best.compute_num_observations()} "
          f"reproj={best.compute_mean_reprojection_error():.3f}px")
PY

pixi run sfm from-colmap-bin "$ROOT/nbr_recon/best" --image-dir "$ROOT" \
    -o "$ROOT/sfmr/nbr_match.sfmr" --tool-name colmap 2>&1 | grep -vE "WARN" | tail -1
echo "--- compare (baseline vs matcher) ---"
pixi run sfm compare "$BASE" "$ROOT/sfmr/nbr_match.sfmr" 2>&1 | grep -vE "WARN" | \
    grep -E "Scale:|VERY SIMILAR|SIGNIFICANT DIFFERENCES"
