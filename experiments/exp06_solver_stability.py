# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 06 — solver stability/flakiness on the dino cluster matches.

The cluster matches and verified COLMAP database are fixed, so any variation
here is the *mapper's* RNG.  Re-run incremental and global SfM N times with
different random seeds on the same verified db, and report, per run:

  * registered images, 3D points, mean reprojection error,
  * agreement with the baseline solve (position RMS and mean rotation error
    after a similarity alignment — the same metric `sfm compare` reports).

Usage:
    pixi run python experiments/exp06_solver_stability.py \
        ../dino_dog_toy_ws --runs 10
"""

from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

import numpy as np
import pycolmap

from sfmtool._sfmtool import RotQuaternion, SfmrReconstruction
from sfmtool.align.core import ImageMatch, estimate_pairwise_alignment


def pose_errors_vs_baseline(base: SfmrReconstruction, target: SfmrReconstruction):
    """Position RMS (after similarity alignment) and per-image rotation errors,
    matching images by workspace-relative name. Returns (rms, pos_errs, rot_errs).
    """
    name_to_t = {n.replace("\\", "/"): i for i, n in enumerate(target.image_names)}
    matches = []
    for i1, n in enumerate(base.image_names):
        j = name_to_t.get(n.replace("\\", "/"))
        if j is None:
            continue
        q1 = RotQuaternion.from_wxyz_array(base.quaternions_wxyz[i1])
        c1 = q1.camera_center(base.translations[i1])
        q2 = RotQuaternion.from_wxyz_array(target.quaternions_wxyz[j])
        c2 = q2.camera_center(target.translations[j])
        matches.append(
            ImageMatch(
                image_name=Path(n).name,
                source_index=j,
                target_index=i1,
                source_quat=q2,
                source_camera_center=c2,
                target_quat=q1,
                target_camera_center=c1,
                quality=1.0,
            )
        )
    res = estimate_pairwise_alignment(
        matches=matches, confidence_threshold=0.0, source_id="t", target_id="b"
    )
    transform = res.transform
    pos, rot = [], []
    for m in matches:
        c2t = transform @ m.source_camera_center
        pos.append(float(np.linalg.norm(c2t - m.target_camera_center)))
        q2t = m.source_quat * transform.rotation.conjugate()
        dot = np.clip(
            abs(np.dot(np.asarray(q2t.to_wxyz_array()),
                       np.asarray(m.target_quat.to_wxyz_array()))),
            0, 1,
        )
        rot.append(float(np.degrees(2 * np.arccos(dot))))
    return res.total_rms_error, np.array(pos), np.array(rot), len(matches)


def run_once(db, root, solver, seed, baseline, tmp):
    pycolmap.set_random_seed(seed)
    out = Path(tmp) / f"{solver}_{seed}"
    out.mkdir(parents=True, exist_ok=True)
    if solver == "incremental":
        o = pycolmap.IncrementalPipelineOptions()
        o.random_seed = seed
        recs = pycolmap.incremental_mapping(db, root, str(out), o)
    else:
        o = pycolmap.GlobalPipelineOptions()
        recs = pycolmap.global_mapping(db, root, str(out), o)
    if not recs:
        return None
    best = max(recs.values(), key=lambda r: r.num_reg_images())
    bdir = out / "best"
    bdir.mkdir(exist_ok=True)
    best.write(str(bdir))
    # convert to sfmr (CLI) to reuse the workspace's .sift/poses
    sfmr_path = out / "m.sfmr"
    import subprocess
    subprocess.run(
        ["pixi", "run", "sfm", "from-colmap-bin", str(bdir),
         "--image-dir", root, "-o", str(sfmr_path), "--tool-name", "colmap"],
        check=True, capture_output=True,
    )
    tgt = SfmrReconstruction.load(str(sfmr_path))
    rms, pos, rot, nmatch = pose_errors_vs_baseline(baseline, tgt)
    return {
        "reg": best.num_reg_images(),
        "points": best.num_points3D(),
        "reproj": best.compute_mean_reprojection_error(),
        "rms": rms,
        "pos_mean": float(pos.mean()),
        "rot_mean": float(rot.mean()),
        "rot_med": float(np.median(rot)),
        "matched": nmatch,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("workspace")
    ap.add_argument("--runs", type=int, default=10)
    args = ap.parse_args()

    root = str(Path(args.workspace).resolve())
    db = f"{root}/nbr.db"
    base_path = sorted(Path(root, "sfmr").glob("*solve*.sfmr"))[0]
    baseline = SfmrReconstruction.load(str(base_path))
    print(f"workspace={root}\nbaseline={base_path.name} db={db}")

    with tempfile.TemporaryDirectory() as tmp:
        for solver in ("incremental", "global"):
            print(f"\n=== {solver} x{args.runs} (seeds 1..{args.runs}) ===")
            print(f"  {'seed':>4} {'reg':>4} {'points':>7} {'reproj':>7} "
                  f"{'posRMS':>7} {'rotMean':>8} {'rotMed':>7}")
            rows = []
            for seed in range(1, args.runs + 1):
                r = run_once(db, root, solver, seed, baseline, tmp)
                if r is None:
                    print(f"  {seed:>4}  FAILED (no model)")
                    continue
                rows.append(r)
                print(f"  {seed:>4} {r['reg']:>4} {r['points']:>7} "
                      f"{r['reproj']:>7.3f} {r['rms']:>7.3f} "
                      f"{r['rot_mean']:>8.3f} {r['rot_med']:>7.3f}")
            if rows:
                def col(k):
                    v = np.array([x[k] for x in rows])
                    return f"{v.mean():.3f}±{v.std():.3f} [{v.min():.3f},{v.max():.3f}]"
                print(f"  --- {len(rows)}/{args.runs} succeeded ---")
                for k in ("reg", "points", "reproj", "rms", "rot_mean"):
                    print(f"    {k:>8}: {col(k)}")


if __name__ == "__main__":
    main()
