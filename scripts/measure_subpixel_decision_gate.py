#!/usr/bin/env python3
# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Decision-gate measurement: compare every available keypoint-localization /
sub-pixel-refinement path on the four checked-in datasets and dump a structured
result.

For each (dataset, variant), this script runs ``embed_patches(...)`` end-to-end,
captures wall time, point count, mean LOO ZNCC / ECC score (depending on the
variant — these proxy different photometric agreement metrics, see Methodology
in the report), and the mean keypoint magnitude shift vs. the baseline (the
no-refinement path).

Output is written to ``reports/<date>-decision-gate-grid-vs-lk-data.json`` (a
machine-readable dump) and used to fill the prose tables in
``reports/<date>-decision-gate-grid-vs-lk.md``.

The variants are:

| Variant            | search_resolution_multiplier | subpixel       |
|--------------------|------------------------------|----------------|
| baseline           | 1.0                          | "none"         |
| grid_m2            | 2.0                          | "none"         |
| grid_m3            | 3.0                          | "none"         |
| lk_per_sweep       | 1.0                          | "lk"           |
| lk_per_move        | 1.0                          | "lk_per_move"  |
| grid_m2_then_lk    | 2.0                          | "lk"           |

Datasets: seoul_bull_sculpture (17), dino_dog_toy (85), seattle_backyard (26),
kerry_park (48; two-sensor rig). The dino_dog_toy run is the slowest by far
(85 high-res images × ~600 patches), so it is opt-in via the env var
``SFMTOOL_DECISION_GATE_FAST=1`` to skip it; the report records the dino number
from a separate slower run.

Run with::

    pixi run -e test python scripts/measure_subpixel_decision_gate.py [DATASET ...]

When no datasets are named, every dataset is run.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

# Add repo src/ so the script runs as a standalone python invocation if needed
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from sfmtool._embed_patches import embed_patches  # noqa: E402
from sfmtool._sfmtool import SfmrReconstruction  # noqa: E402
from sfmtool._workspace_image import read_workspace_image  # noqa: E402

DATASETS: dict[str, str] = {
    "seoul_bull": "seoul_bull_ws/sfmr/20260621-00-solve-seoul_bull_sculpture_1-17.sfmr",
    "dino_dog_toy": "dino_dog_toy_ws/sfmr/20260621-00-solve-dino_dog_toy_1-85.sfmr",
    "seattle_backyard": "seattle_backyard_ws/sfmr/20260621-00-solve-seattle_backyard_1-26.sfmr",
    "kerry_park": "kerry_park_ws/sfmr/20260621-00-solve-frame_1-24.sfmr",
}

# (label, kwargs for embed_patches)
VARIANTS: list[tuple[str, dict]] = [
    ("baseline", dict(search_resolution_multiplier=1.0, subpixel="none")),
    ("grid_m2", dict(search_resolution_multiplier=2.0, subpixel="none")),
    ("grid_m3", dict(search_resolution_multiplier=3.0, subpixel="none")),
    ("lk_per_sweep", dict(search_resolution_multiplier=1.0, subpixel="lk")),
    ("lk_per_move", dict(search_resolution_multiplier=1.0, subpixel="lk_per_move")),
    ("grid_m2_then_lk", dict(search_resolution_multiplier=2.0, subpixel="lk")),
]


def load_recon_and_images(
    sfmr_path: Path,
) -> tuple[SfmrReconstruction, list[np.ndarray]]:
    recon = SfmrReconstruction.load(str(sfmr_path))
    images = [
        read_workspace_image(recon.workspace_dir, name) for name in recon.image_names
    ]
    return recon, images


def per_obs_keypoints_by_world(
    recon: SfmrReconstruction,
) -> dict[tuple[bytes, int], np.ndarray]:
    """Map ``(rounded world position bytes, image_index) -> [x, y]`` for an
    embedded_patches recon's inline keypoints. Used to compare per-observation
    locations across variants where the *output* point index renumbering can
    differ (compaction renumbers survivors into a dense set, so a single dropped
    point shifts every later index by one).

    The world position is rounded to 6 decimal places and packed to bytes for
    use as a hash key — the same source-point will round-trip to identical bits
    across variants (compaction never re-bundles, only culls + renumbers), so
    keying on the world position is the invariant join key.
    """
    kxy = np.asarray(recon.keypoints_xy, dtype=np.float64)
    tpid = np.asarray(recon.track_point_ids, dtype=np.int64)
    timg = np.asarray(recon.track_image_indexes, dtype=np.int64)
    positions = np.asarray(recon.positions, dtype=np.float64)
    out: dict[tuple[bytes, int], np.ndarray] = {}
    for k, (p, i) in enumerate(zip(tpid, timg)):
        key_pos = np.round(positions[int(p)], 6).tobytes()
        out[(key_pos, int(i))] = kxy[k]
    return out


def keypoint_shift_vs_baseline(
    baseline: SfmrReconstruction, variant: SfmrReconstruction
) -> dict[str, float]:
    """Mean / median per-observation keypoint shift, in source-image px, over the
    observations present in BOTH the baseline and the variant. We can't compare a
    point that one variant culled and the other kept; we report ``n_overlap``
    so the magnitude is in context. See :func:`per_obs_keypoints_by_world` for
    the join key (world position, not the renumbered output index).
    """
    b = per_obs_keypoints_by_world(baseline)
    v = per_obs_keypoints_by_world(variant)
    shared = b.keys() & v.keys()
    if not shared:
        return {
            "mean_shift_px": float("nan"),
            "median_shift_px": float("nan"),
            "n_overlap": 0,
        }
    diffs = np.array([np.linalg.norm(b[k] - v[k]) for k in shared])
    return {
        "mean_shift_px": float(diffs.mean()),
        "median_shift_px": float(np.median(diffs)),
        "n_overlap": int(len(shared)),
    }


def measure_photometric_agreement(
    recon: SfmrReconstruction, images: list[np.ndarray], resolution: int = 12
) -> dict[str, float]:
    """Build a patch cloud over the embedded_patches recon's stored frames, run
    ``refine_keypoints`` with **zero** GN steps (seed-only) against the stored
    inline keypoints, and report the mean post-refinement ECC.

    This is the cross-variant comparable measurement: every variant's stored
    keypoints are scored against the same ECC objective (channel-averaged
    windowed ZNCC against the IRLS consensus rendered from the stored
    keypoints). It is computed AFTER ``embed_patches`` writes the result, so it
    sees only the inline keypoints that survived compaction.

    The metric inherits its limitations from the refiner's ``score`` field:
    NaN for any point left with < 2 views (no consensus), so we only average
    finite values and report the count.
    """
    cloud = recon.patches
    if cloud is None or len(cloud) == 0:
        return {"mean_ecc": float("nan"), "n_scored": 0}

    # Per-point inline view set + keypoints from the compacted recon (the data
    # the metric is measuring).
    tpid = np.asarray(recon.track_point_ids, dtype=np.int64)
    timg = np.asarray(recon.track_image_indexes, dtype=np.int64)
    kxy = np.asarray(recon.keypoints_xy, dtype=np.float64)
    cloud_pids = set(int(p) for p in cloud.point_ids)
    view_sets: dict[int, list[int]] = {}
    seeds: dict[int, list[list[float]]] = {}
    for k, (p, i) in enumerate(zip(tpid.tolist(), timg.tolist())):
        if p not in cloud_pids:
            continue
        view_sets.setdefault(p, []).append(i)
        seeds.setdefault(p, []).append([float(kxy[k, 0]), float(kxy[k, 1])])

    if not view_sets:
        return {"mean_ecc": float("nan"), "n_scored": 0}

    results = cloud.refine_keypoints(
        recon,
        images,
        view_sets=view_sets,
        starting_keypoints=seeds,
        point_ids=list(view_sets.keys()),
        resolution=resolution,
        max_gn_steps=0,  # seed-only: just score the seed against the consensus.
    )
    scores = []
    for r in results:
        s = np.asarray(r["scores"], dtype=np.float64)
        scores.extend(s[np.isfinite(s)].tolist())
    if not scores:
        return {"mean_ecc": float("nan"), "n_scored": 0}
    return {"mean_ecc": float(np.mean(scores)), "n_scored": len(scores)}


def run_variant(
    recon: SfmrReconstruction,
    images: list[np.ndarray],
    variant: dict,
) -> tuple[SfmrReconstruction, float]:
    t0 = time.perf_counter()
    out = embed_patches(recon, images, patch_size=10.0, **variant)
    dt = time.perf_counter() - t0
    return out, dt


def measure_dataset(dataset_label: str, sfmr_path: Path) -> dict:
    print(f"\n=== {dataset_label} ({sfmr_path.name}) ===", flush=True)
    recon, images = load_recon_and_images(sfmr_path)
    print(
        f"  images={recon.image_count} input_points={recon.point_count}",
        flush=True,
    )

    rows: list[dict] = []
    baseline_recon: SfmrReconstruction | None = None
    for label, kw in VARIANTS:
        print(f"  - {label} ... ", end="", flush=True)
        try:
            out, dt = run_variant(recon, images, kw)
        except Exception as e:  # surface, don't crash the run
            print(f"FAILED: {e}", flush=True)
            rows.append(
                {
                    "variant": label,
                    "wall_secs": float("nan"),
                    "out_points": 0,
                    "mean_ecc": float("nan"),
                    "n_scored": 0,
                    "mean_shift_px": float("nan"),
                    "median_shift_px": float("nan"),
                    "n_overlap": 0,
                    "kwargs": kw,
                    "error": str(e),
                }
            )
            continue
        agree = measure_photometric_agreement(out, images)
        if label == "baseline":
            baseline_recon = out
            shift = {"mean_shift_px": 0.0, "median_shift_px": 0.0, "n_overlap": -1}
        else:
            shift = keypoint_shift_vs_baseline(baseline_recon, out)
        row = {
            "variant": label,
            "wall_secs": dt,
            "out_points": int(out.point_count),
            "mean_ecc": agree["mean_ecc"],
            "n_scored": agree["n_scored"],
            **shift,
            "kwargs": kw,
        }
        rows.append(row)
        print(
            f"{dt:5.2f}s  points={out.point_count:4d}  "
            f"mean_ecc={agree['mean_ecc']:.4f}  "
            f"shift={shift['mean_shift_px']:.4f}px (median {shift['median_shift_px']:.4f}px)",
            flush=True,
        )
    return {
        "dataset": dataset_label,
        "sfmr_path": str(sfmr_path),
        "image_count": recon.image_count,
        "input_point_count": recon.point_count,
        "rows": rows,
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "datasets",
        nargs="*",
        choices=sorted(DATASETS.keys()) + [],
        help="Subset of datasets to measure (default: all)",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=ROOT / "reports" / "2026-06-27-decision-gate-grid-vs-lk-data.json",
        help="Where to write the structured measurement results",
    )
    args = p.parse_args()

    selected = args.datasets or list(DATASETS.keys())
    fast = os.environ.get("SFMTOOL_DECISION_GATE_FAST")
    if fast and not args.datasets:
        selected = [d for d in selected if d != "dino_dog_toy"]
        print(f"[fast mode] skipping dino_dog_toy; running {selected}", flush=True)

    results = []
    for d in selected:
        sfmr = ROOT / DATASETS[d]
        if not sfmr.exists():
            print(f"  SKIP {d}: {sfmr} missing", flush=True)
            continue
        results.append(measure_dataset(d, sfmr))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(results, indent=2))
    print(f"\nWrote {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
