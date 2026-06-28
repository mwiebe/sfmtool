#!/usr/bin/env python
"""Compare exhaustive vs "+" descent on the patch-normal grid search.

Prototype for the experiment described on the
``claude/quirky-heisenberg-2vlpmp`` branch:

  At each coarse-to-fine level the production refiner scores every cell
  of an ``init_steps × init_steps`` grid (49 evals for the default
  ``init_steps = 7``, ~135 over 3 levels). On reconstructions of any
  size this dominates the wall time.

  The "+" descent alternative walks the integer grid from the level
  center: at each cell it evaluates the 4 axis neighbors (skipping
  out-of-bounds / out-of-disk / already-visited) and moves to the best
  improver. Stops when no axis neighbor beats the current cell. Each
  cell is scored at most once (visited cache). Adaptive: fewer evals
  when the seed is already near the level argmax, more when it has to
  walk further.

This script runs both strategies on the same reconstruction with
otherwise-identical params and reports:

  * total wall time and per-call overhead
  * total ``Φ`` evaluations (across all rayon threads) and the per-patch
    average
  * the distribution of disagreement angles between the two strategies'
    final normals (median, p90, max, fraction > 1°, fraction > 5°)
  * the per-patch ``Φ`` delta (plus minus exhaustive) — negative means
    "+"-descent landed on a worse cell

Example:

  pixi run -e test python scripts/exp_plus_descent_normal.py \\
      /tmp/pxv/seoul/sfmr/seoul.sfmr --n-points 300 --seed 0

The ``SFMTOOL_PROFILE=1`` environment variable is set automatically so
the Rust ``N_EVAL`` counter is live; pass ``--phase-report`` to also
print the per-phase timing breakdown to stderr.

For per-cell visibility into the "+"-descent walk (which `(i, j)` cells
were evaluated, which were cache hits, which were skipped, the per-level
walk length), set ``SFMTOOL_PLUS_DESCENT_TRACE=1`` and pass
``--n-points 1`` so the rayon parallelism doesn't interleave the traces.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np


def angle_deg(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Per-row angle between unit-ish vectors, in degrees."""
    an = a / np.linalg.norm(a, axis=1, keepdims=True)
    bn = b / np.linalg.norm(b, axis=1, keepdims=True)
    dot = np.clip((an * bn).sum(axis=1), -1.0, 1.0)
    return np.degrees(np.arccos(dot))


def run_one(s, recon, images, ids, params, strategy):
    """Refine `ids` with the given strategy and return timings + per-patch arrays."""
    cloud = s.PatchCloud.from_reconstruction(
        recon,
        normal="mean_viewing",
        extent_value=5.0,
        exclude_points_at_infinity=True,
    )
    # The eval counter resets at the start of `refine_patch_cloud`, so reading
    # it right after the call gives this call's count.
    t0 = time.perf_counter()
    res = cloud.refine_normals(
        recon,
        images,
        point_indexes=ids.tolist(),
        search_strategy=strategy,
        **params,
    )
    wall = time.perf_counter() - t0
    evals = s.normal_refine_eval_count()
    return {
        "strategy": strategy,
        "wall_s": wall,
        "evals": evals,
        "normal": np.asarray(res["normal"], dtype=np.float64),
        "photoconsistency": np.asarray(res["photoconsistency"], dtype=np.float64),
        "init_photoconsistency": np.asarray(
            res["init_photoconsistency"], dtype=np.float64
        ),
    }


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("sfmr", help="path to the .sfmr reconstruction")
    ap.add_argument(
        "--n-points",
        type=int,
        default=300,
        help="random subset of finite points to refine (0 = all)",
    )
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--resolution", type=int, default=16)
    ap.add_argument("--init-steps", type=int, default=7)
    ap.add_argument("--refine-levels", type=int, default=3)
    ap.add_argument("--angular-range-deg", type=float, default=25.0)
    ap.add_argument(
        "--sampler", default="anisotropic", choices=["bilinear", "anisotropic"]
    )
    ap.add_argument("--cache", default="fronto", choices=["off", "fronto"])
    ap.add_argument(
        "--phase-report",
        action="store_true",
        help="also print the Rust per-phase timing summary to stderr "
        "(SFMTOOL_PROFILE=1)",
    )
    ap.add_argument(
        "--out",
        default=None,
        help="append a JSON record to this path (besides the stdout summary)",
    )
    args = ap.parse_args()

    # The Rust eval counter is gated on SFMTOOL_PROFILE; always-on here so the
    # script's measurements are honest. `--phase-report` additionally lets the
    # per-phase timing summary through to stderr.
    os.environ["SFMTOOL_PROFILE"] = "1"

    import cv2  # type: ignore[import-not-found]

    import sfmtool._sfmtool as s

    if s.build_profile() != "release":
        sys.exit(
            "sfmtool._sfmtool was built in debug; rebuild with "
            "`pixi run maturin develop --release` before measuring."
        )

    recon = s.SfmrReconstruction.load(args.sfmr)
    ws = recon.workspace_dir
    t0 = time.perf_counter()
    images = [
        np.ascontiguousarray(cv2.imread(os.path.join(ws, name), cv2.IMREAD_COLOR))
        for name in recon.image_names
    ]
    load_s = time.perf_counter() - t0

    # Sample the finite-point subset once so both strategies refine identical patches.
    cloud0 = s.PatchCloud.from_reconstruction(
        recon,
        normal="mean_viewing",
        extent_value=5.0,
        exclude_points_at_infinity=True,
    )
    all_ids = np.asarray(cloud0.point_indexes, dtype=np.uint32)
    rng = np.random.default_rng(args.seed)
    if args.n_points and args.n_points < len(all_ids):
        ids = np.sort(rng.choice(all_ids, size=args.n_points, replace=False))
    else:
        ids = np.sort(all_ids)
    # Map cloud indices to keep only the refined rows when comparing.
    index_of = {pid: i for i, pid in enumerate(all_ids.tolist())}
    rows = np.asarray([index_of[pid] for pid in ids.tolist()])

    dataset = os.path.splitext(os.path.basename(args.sfmr))[0]
    print(
        f"# {dataset}: {len(recon.image_names)} images (loaded in {load_s:.1f}s), "
        f"refining {len(ids)} points; init_steps={args.init_steps}, "
        f"refine_levels={args.refine_levels}, angular_range={args.angular_range_deg}°, "
        f"resolution={args.resolution}, sampler={args.sampler}, cache={args.cache}",
        file=sys.stderr,
    )

    params = dict(
        resolution=args.resolution,
        angular_range_deg=args.angular_range_deg,
        init_steps=args.init_steps,
        refine_levels=args.refine_levels,
        sampler=args.sampler,
        cache=args.cache,
        objective="robust",
        robust_iters=3,
    )

    runs = {}
    for strategy in ("exhaustive", "plus_descent"):
        runs[strategy] = run_one(s, recon, images, ids, params, strategy)
        r = runs[strategy]
        print(
            f"[{strategy:>13}] wall {r['wall_s']:.3f}s  evals {r['evals']:>9}  "
            f"evals/patch {r['evals'] / max(len(ids), 1):.1f}",
            file=sys.stderr,
        )

    # Per-patch comparison on the refined rows.
    a = runs["exhaustive"]
    b = runs["plus_descent"]
    na = a["normal"][rows]
    nb = b["normal"][rows]
    pa = a["photoconsistency"][rows]
    pb = b["photoconsistency"][rows]
    scored = np.isfinite(pa) & np.isfinite(pb)
    angles = angle_deg(na[scored], nb[scored])
    dphi = pb[scored] - pa[scored]

    n_scored = int(scored.sum())
    summary = {
        "dataset": dataset,
        "n_points_requested": int(len(ids)),
        "n_scored": n_scored,
        "params": params,
        "exhaustive": {
            "wall_s": round(a["wall_s"], 4),
            "evals": int(a["evals"]),
            "evals_per_patch": round(a["evals"] / max(len(ids), 1), 2),
            "median_phi": round(float(np.median(pa[scored])), 5) if n_scored else None,
        },
        "plus_descent": {
            "wall_s": round(b["wall_s"], 4),
            "evals": int(b["evals"]),
            "evals_per_patch": round(b["evals"] / max(len(ids), 1), 2),
            "median_phi": round(float(np.median(pb[scored])), 5) if n_scored else None,
        },
        "speedup_wall": round(a["wall_s"] / max(b["wall_s"], 1e-9), 2),
        "eval_ratio": round(a["evals"] / max(b["evals"], 1), 2),
        "agreement": {
            "median_disagree_deg": round(float(np.median(angles)), 4)
            if n_scored
            else None,
            "p90_disagree_deg": round(float(np.percentile(angles, 90)), 4)
            if n_scored
            else None,
            "max_disagree_deg": round(float(angles.max()), 4) if n_scored else None,
            "frac_within_1deg": round(float((angles <= 1.0).mean()), 4)
            if n_scored
            else None,
            "frac_within_5deg": round(float((angles <= 5.0).mean()), 4)
            if n_scored
            else None,
        },
        "phi_delta": {
            "mean": round(float(dphi.mean()), 6) if n_scored else None,
            "median": round(float(np.median(dphi)), 6) if n_scored else None,
            "p10": round(float(np.percentile(dphi, 10)), 6) if n_scored else None,
            "p90": round(float(np.percentile(dphi, 90)), 6) if n_scored else None,
            "frac_worse": round(float((dphi < -1e-6).mean()), 4) if n_scored else None,
        },
    }

    print(json.dumps(summary, indent=2))
    if args.out:
        with open(args.out, "a", encoding="utf-8") as f:
            f.write(json.dumps(summary) + "\n")


if __name__ == "__main__":
    main()
