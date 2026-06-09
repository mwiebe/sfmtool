#!/bin/bash
# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0
#
# End-to-end reconstruction from the background-floor matcher (exp05): build
# matches, verify (seeded), incremental SfM (seeded), compare to the workspace's
# baseline solve. MODE selects the matcher form: bgfloor (per-descriptor edges)
# or bgclusters (materialized clusters). Distinct per-mode output names so runs
# never clobber each other or the nbr_* / pct_* solves.
#
# Usage:  bash experiments/run_bgfloor.sh <workspace> [alpha] [d] [mode]
set -e
export PATH="$HOME/.pixi/bin:$PATH"

WS="$1"
ALPHA="${2:-0.8}"
D="${3:-28}"
MODE="${4:-bgfloor}"
case "$MODE" in
  bgfloor) PFX=bgf ;;
  bgclusters) PFX=bgc ;;
  *) echo "unknown mode: $MODE" >&2; exit 1 ;;
esac
ROOT="$(cd "$(dirname "$0")/.." && pwd)/$WS"
BASE=$(ls "$ROOT"/sfmr/*solve*.sfmr | head -1)
echo "============================== $WS ($MODE a=$ALPHA d=$D) =============================="

cd "$(dirname "$0")"
pixi run -e experiments python -u exp05_cluster_match.py "$ROOT/sfmr/"'*solve*.sfmr' \
    --out "$ROOT/matches/$PFX.matches" --mode "$MODE" --bg-alpha "$ALPHA" --bg-d "$D" \
    2>&1 | grep -vE "WARN" | grep -E "after one"

cd "$ROOT/.."
rm -f "$ROOT/$PFX.db"
pixi run sfm to-colmap-db "$ROOT/matches/$PFX.matches" --out-db "$ROOT/$PFX.db" 2>&1 \
    | grep -vE "WARN" >/dev/null
pixi run python - "$ROOT" "$PFX" <<'PY' 2>&1 | grep -vE "WARN" >/dev/null
import sys, pycolmap
pycolmap.set_random_seed(42)
pycolmap.geometric_verification(f"{sys.argv[1]}/{sys.argv[2]}.db")
PY

OUT="$ROOT/${PFX}_inc"; rm -rf "$OUT"
REC=$(pixi run python - "$ROOT" "$PFX" "$OUT" <<'PY' 2>&1 | grep -vE "WARN" | grep -E "^recon:"
import os, sys, pycolmap
root, pfx, out = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(out, exist_ok=True)
pycolmap.set_random_seed(42)
o = pycolmap.IncrementalPipelineOptions(); o.random_seed = 42
recs = pycolmap.incremental_mapping(f"{root}/{pfx}.db", root, out, o)
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
        -o "$ROOT/sfmr/${PFX}_incremental.sfmr" --tool-name colmap >/dev/null 2>&1
    verdict=$(pixi run sfm compare "$BASE" "$ROOT/sfmr/${PFX}_incremental.sfmr" 2>&1 \
        | grep -vE "WARN" | grep -oE "VERY SIMILAR|SIGNIFICANT DIFFERENCES" | head -1)
fi
printf "  %s  reg=%s points=%s reproj=%s  %s\n" "$MODE" "$reg" "$pts" "$reproj" "$verdict"
