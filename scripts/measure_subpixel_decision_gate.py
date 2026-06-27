#!/usr/bin/env python3
# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Decision-gate measurement: compare every available keypoint-localization /
sub-pixel-refinement path on the four checked-in datasets and dump a structured
result.

For each (dataset, variant), this script runs ``embed_patches(...)`` end-to-end
(or, for variants that exercise LK parameters the production API does not
expose, runs the pipeline by hand and slips a custom :meth:`refine_keypoints`
call into step 3.5), captures wall time, point count, mean ECC, and
per-observation keypoint shift relative to:

- the **baseline** variant (the no-refinement production path), and
- the **SIFT seeds** that ``to_embedded_patches`` copied inline at step 0
  (before the normal-refine → select-views → localize → optional LK pipeline
  ever moved them).

Output is written to ``reports/<date>-decision-gate-grid-vs-lk-data.json`` (a
machine-readable dump) and used to fill the prose tables in
``reports/<date>-decision-gate-grid-vs-lk.md``.

The variants are:

| Variant                    | search_resolution_multiplier | step 3.5                              |
|----------------------------|------------------------------|---------------------------------------|
| baseline                   | 1.0                          | none                                  |
| grid_m2                    | 2.0                          | none                                  |
| grid_m3                    | 3.0                          | none                                  |
| lk_per_sweep               | 1.0                          | LK, ``max_outer_sweeps=1``            |
| lk_per_move                | 1.0                          | LK, per-move, ``max_outer_sweeps=5``  |
| grid_m2_then_lk            | 2.0                          | LK, ``max_outer_sweeps=1``            |
| lk_per_sweep_aniso         | 1.0                          | LK, ``sampler='anisotropic'``         |
| lk_per_sweep_tight_offset  | 1.0                          | LK, ``max_offset_px=1.0``             |
| lk_per_move_10sweeps       | 1.0                          | LK per-move, ``max_outer_sweeps=10``  |

The last three exercise LK parameters the production ``embed_patches`` API
does not expose (the production decision was deliberately to keep the surface
at ``subpixel = none|lk|lk_per_move``); the script runs those by replaying
the full pipeline and calling ``cloud.refine_keypoints(...)`` with custom
kwargs.

Datasets: seoul_bull_sculpture (17), dino_dog_toy (85), seattle_backyard (26),
kerry_park (48; two-sensor rig). The dino_dog_toy run is the slowest by far
(85 high-res images × ~600 patches); previous rounds skipped it. This round
includes it, but only for the not-prohibitively-slow variants (baseline,
grid_m2, lk_per_sweep, lk_per_move, grid_m2_then_lk plus the three custom-LK
variants — i.e. everything except ``grid_m3``, the ~21× baseline supersampled
grid). If even that is too slow, set ``SFMTOOL_DINO_SUBSAMPLE_STRIDE`` to take
every Nth point in the reconstruction's point order; the sampling is recorded
in the result.

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
from typing import Any

import numpy as np

# Add repo src/ so the script runs as a standalone python invocation if needed
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from sfmtool._embed_patches import compact_to_embedded_patches  # noqa: E402
from sfmtool._embed_patches import embed_patches  # noqa: E402
from sfmtool._sfmtool import PatchCloud, SfmrReconstruction  # noqa: E402
from sfmtool._workspace_image import read_workspace_image  # noqa: E402

DATASETS: dict[str, str] = {
    "seoul_bull": "seoul_bull_ws/sfmr/20260621-00-solve-seoul_bull_sculpture_1-17.sfmr",
    "dino_dog_toy": "dino_dog_toy_ws/sfmr/20260621-00-solve-dino_dog_toy_1-85.sfmr",
    "seattle_backyard": "seattle_backyard_ws/sfmr/20260621-00-solve-seattle_backyard_1-26.sfmr",
    "kerry_park": "kerry_park_ws/sfmr/20260621-00-solve-frame_1-24.sfmr",
}

# Variants: each entry is (label, dict) where the dict has either
#   "embed_kwargs": dict          -- pass straight to embed_patches(...)
# or
#   "refine_kwargs": dict         -- run the pipeline by hand and call
#                                    cloud.refine_keypoints(**refine_kwargs)
#                                    as step 3.5, seeded at the localizer's
#                                    keypoints. ``search_resolution_multiplier``
#                                    may be set alongside for the localize step.
VARIANTS: list[tuple[str, dict]] = [
    # Production-API variants — drive embed_patches directly.
    (
        "baseline",
        {"embed_kwargs": dict(search_resolution_multiplier=1.0, subpixel="none")},
    ),
    (
        "grid_m2",
        {"embed_kwargs": dict(search_resolution_multiplier=2.0, subpixel="none")},
    ),
    (
        "grid_m3",
        {"embed_kwargs": dict(search_resolution_multiplier=3.0, subpixel="none")},
    ),
    (
        "lk_per_sweep",
        {"embed_kwargs": dict(search_resolution_multiplier=1.0, subpixel="lk")},
    ),
    (
        "lk_per_move",
        {
            "embed_kwargs": dict(
                search_resolution_multiplier=1.0, subpixel="lk_per_move"
            )
        },
    ),
    (
        "grid_m2_then_lk",
        {"embed_kwargs": dict(search_resolution_multiplier=2.0, subpixel="lk")},
    ),
    # Custom-refine variants — exercise LK params that the production
    # embed_patches API does not expose. Each runs the pipeline by hand
    # (same as `subpixel="none"`) and then calls refine_keypoints with the
    # given kwargs. The plausible fix-candidates for the kerry_park regression
    # come first.
    (
        "lk_per_sweep_aniso",
        {
            "search_resolution_multiplier": 1.0,
            "refine_kwargs": dict(
                max_outer_sweeps=1,
                consensus_refresh="per_sweep",
                sampler="anisotropic",
            ),
        },
    ),
    (
        "lk_per_sweep_tight_offset",
        {
            "search_resolution_multiplier": 1.0,
            "refine_kwargs": dict(
                max_outer_sweeps=1,
                consensus_refresh="per_sweep",
                max_offset_px=1.0,
            ),
        },
    ),
    (
        "lk_per_move_10sweeps",
        {
            "search_resolution_multiplier": 1.0,
            "refine_kwargs": dict(
                max_outer_sweeps=10,
                consensus_refresh="per_move",
            ),
        },
    ),
]

# Per-dataset variant overrides — drop the prohibitively-slow ones from dino.
DATASET_VARIANT_SKIPS: dict[str, set[str]] = {
    # dino is 85 @ 2040x1536; grid_m3 is ~21× baseline cost (from prior data on
    # the smaller datasets), which would push the dino run over an hour by
    # itself. Every other variant still gets measured.
    "dino_dog_toy": {"grid_m3"},
}


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


def keypoint_shift_summary(
    reference: dict[tuple[bytes, int], np.ndarray],
    variant: dict[tuple[bytes, int], np.ndarray],
) -> dict[str, float]:
    """Mean / median / p95 per-observation keypoint shift, in source-image px,
    over the observations present in BOTH ``reference`` and ``variant``. We
    can't compare an observation that one variant culled and the other kept;
    we report ``n_overlap`` so the magnitude is in context.
    """
    shared = reference.keys() & variant.keys()
    if not shared:
        return {
            "mean_shift_px": float("nan"),
            "median_shift_px": float("nan"),
            "p95_shift_px": float("nan"),
            "n_overlap": 0,
        }
    diffs = np.array([np.linalg.norm(reference[k] - variant[k]) for k in shared])
    return {
        "mean_shift_px": float(diffs.mean()),
        "median_shift_px": float(np.median(diffs)),
        "p95_shift_px": float(np.percentile(diffs, 95)),
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
        max_gn_steps=0,
    )
    scores = []
    for r in results:
        s = np.asarray(r["scores"], dtype=np.float64)
        scores.extend(s[np.isfinite(s)].tolist())
    if not scores:
        return {"mean_ecc": float("nan"), "n_scored": 0}
    return {"mean_ecc": float(np.mean(scores)), "n_scored": len(scores)}


def sift_seed_keypoints(
    recon: SfmrReconstruction,
    *,
    patch_size: float,
) -> dict[tuple[bytes, int], np.ndarray]:
    """Snapshot the per-observation SIFT detection keypoints copied inline by
    ``to_embedded_patches``, keyed by ``(rounded world position bytes,
    image_index)``. This is the input the rest of the pipeline operates on —
    every variant's per-observation shift relative to this is the *cumulative*
    photometric move the pipeline applied to that observation (localizer +
    optional LK).

    The conversion is the same one ``embed_patches`` uses at step 0, with
    identical normal/extent policy, so the resulting per-observation keypoints
    are bit-equivalent to what every variant starts from.
    """
    half_extent = patch_size / 2.0
    embedded = recon.to_embedded_patches(
        normal="mean_viewing", extent="feature_size", extent_value=half_extent
    )
    return per_obs_keypoints_by_world(embedded)


def _refine_keypoints_custom(
    cloud: PatchCloud,
    embedded: SfmrReconstruction,
    images: list[np.ndarray],
    localizations: list[dict[str, Any]],
    *,
    refine_kwargs: dict[str, Any],
    resolution: int,
) -> list[dict[str, Any]]:
    """Like ``sfmtool._embed_patches._refine_subpixel`` but takes the LK kwargs
    verbatim instead of mapping a ``mode`` to a fixed kwarg set. The splice
    semantics are identical: per-point view sets + per-view seeds are derived
    from the localizer's output, the refiner runs, and each point's per-view
    keypoints are spliced back into its localization dict while every other
    field stays put.
    """
    view_sets: dict[int, list[int]] = {}
    seeds: dict[int, list[list[float]]] = {}
    for loc in localizations:
        pid = int(loc["point_id"])
        views = np.asarray(loc["views"], dtype=np.uint32).tolist()
        kpts = np.asarray(loc["keypoints"], dtype=np.float64).reshape(-1, 2)
        if not views:
            continue
        view_sets[pid] = views
        seeds[pid] = [[float(p[0]), float(p[1])] for p in kpts]

    if not view_sets:
        return localizations

    refined = cloud.refine_keypoints(
        embedded,
        images,
        view_sets=view_sets,
        starting_keypoints=seeds,
        point_ids=list(view_sets.keys()),
        resolution=resolution,
        **refine_kwargs,
    )

    refined_by_pid = {int(r["point_id"]): r for r in refined}
    out: list[dict[str, Any]] = []
    for loc in localizations:
        pid = int(loc["point_id"])
        r = refined_by_pid.get(pid)
        if r is None:
            out.append(loc)
            continue
        r_views = np.asarray(r["views"], dtype=np.uint32)
        r_kpts = np.asarray(r["keypoints"], dtype=np.float64).reshape(-1, 2)
        l_views = np.asarray(loc["views"], dtype=np.uint32)
        l_kpts = np.asarray(loc["keypoints"], dtype=np.float64).reshape(-1, 2)
        r_map = {int(v): r_kpts[i] for i, v in enumerate(r_views.tolist())}
        new_kpts = np.array(
            [r_map.get(int(v), l_kpts[i]) for i, v in enumerate(l_views.tolist())],
            dtype=np.float64,
        ).reshape(-1, 2)
        new_loc = dict(loc)
        new_loc["keypoints"] = new_kpts
        out.append(new_loc)
    return out


def run_custom_refine_variant(
    recon: SfmrReconstruction,
    images: list[np.ndarray],
    *,
    search_resolution_multiplier: float,
    refine_kwargs: dict[str, Any],
    patch_size: float = 10.0,
    min_relative_zncc: float = 0.7,
    max_shift_px: float = 3.0,
    min_views: int = 2,
    max_iters: int = 5,
    search: float = 6.0,
    resolution: int = 24,
) -> SfmrReconstruction:
    """Reproduce ``embed_patches`` end to end, but at step 3.5 call
    ``cloud.refine_keypoints(**refine_kwargs)`` directly (rather than going
    through the production ``subpixel="lk"|"lk_per_move"`` enum, which hard-
    codes its LK params). Used by the measurement script to exercise LK
    knobs (anisotropic sampler, tighter ``max_offset_px``, more outer sweeps)
    that the production API does not expose by design.
    """
    half_extent = patch_size / 2.0
    embedded = recon.to_embedded_patches(
        normal="mean_viewing", extent="feature_size", extent_value=half_extent
    )
    cloud = embedded.patches
    if cloud is None:
        raise ValueError("to_embedded_patches produced no patch frames to refine")
    refine = cloud.refine_normals(
        embedded,
        images,
        resolution=resolution,
        use_stored_keypoints=True,
        render_bitmaps=True,
    )
    selections = cloud.select_views(
        embedded, images, min_relative_zncc=min_relative_zncc, resolution=resolution
    )
    view_sets = {
        int(s["point_id"]): np.asarray(s["admitted"]).tolist() for s in selections
    }
    localizations = cloud.localize_keypoints(
        embedded,
        images,
        view_sets=view_sets,
        max_iters=max_iters,
        search=search,
        max_shift_px=max_shift_px,
        min_relative_zncc=min_relative_zncc,
        resolution=resolution,
        search_resolution_multiplier=search_resolution_multiplier,
    )
    localizations = _refine_keypoints_custom(
        cloud,
        embedded,
        images,
        localizations,
        refine_kwargs=refine_kwargs,
        resolution=resolution,
    )
    hashes = embedded.image_file_hashes
    return compact_to_embedded_patches(
        recon,
        cloud,
        localizations,
        hashes,
        patch_bitmaps=refine.get("bitmaps"),
        min_views=min_views,
    )


def run_variant(
    recon: SfmrReconstruction,
    images: list[np.ndarray],
    variant: dict,
) -> tuple[SfmrReconstruction, float]:
    t0 = time.perf_counter()
    if "embed_kwargs" in variant:
        out = embed_patches(recon, images, patch_size=10.0, **variant["embed_kwargs"])
    elif "refine_kwargs" in variant:
        out = run_custom_refine_variant(
            recon,
            images,
            search_resolution_multiplier=variant.get(
                "search_resolution_multiplier", 1.0
            ),
            refine_kwargs=variant["refine_kwargs"],
        )
    else:
        raise ValueError(
            f"variant dict must carry 'embed_kwargs' or 'refine_kwargs', got: {variant!r}"
        )
    dt = time.perf_counter() - t0
    return out, dt


def variant_kwargs_for_record(variant: dict) -> dict:
    """A JSON-friendly representation of the variant's knobs for the dump."""
    if "embed_kwargs" in variant:
        return {"embed_kwargs": variant["embed_kwargs"]}
    return {
        "search_resolution_multiplier": variant.get(
            "search_resolution_multiplier", 1.0
        ),
        "refine_kwargs": variant["refine_kwargs"],
    }


def subsampled_recon(recon: SfmrReconstruction, stride: int) -> SfmrReconstruction:
    """Take every ``stride``-th point in the reconstruction's existing point
    order. Used to fit the dino sweep into a reasonable budget when the
    full-point run is too slow.

    The implementation re-exports the recon as a ``sift_files`` recon (the
    pipeline's only valid input) keeping only the stride-selected points and
    the observations that reference them; the resulting recon is functionally
    a smaller version of the same dataset (same images, same cameras, fewer
    points and fewer per-point observations).
    """
    if stride <= 1:
        return recon
    n = recon.point_count
    keep = np.arange(0, n, stride, dtype=np.int64)
    if len(keep) == 0:
        raise ValueError(f"stride {stride} drops every point")
    # Slice per-point arrays to the survivors.
    positions = np.asarray(recon.positions_xyzw, dtype=np.float64)[keep]
    colors = np.asarray(recon.colors, dtype=np.uint8)[keep]
    errors = np.asarray(recon.errors, dtype=np.float32)[keep]
    if recon.has_normals:
        normals = np.asarray(recon.normals, dtype=np.float32)[keep]
    else:
        normals = None
    # Renumber observations: keep only obs whose point id survived; remap
    # surviving point ids to dense 0..len(keep)-1.
    old_to_new = -np.ones(n, dtype=np.int64)
    old_to_new[keep] = np.arange(len(keep), dtype=np.int64)
    tpid = np.asarray(recon.track_point_ids, dtype=np.int64)
    timg = np.asarray(recon.track_image_indexes, dtype=np.uint32)
    tfid = np.asarray(recon.track_feature_indexes, dtype=np.uint32)
    keep_mask = old_to_new[tpid] >= 0
    new_tpid = old_to_new[tpid[keep_mask]].astype(np.uint32)
    new_timg = timg[keep_mask]
    new_tfid = tfid[keep_mask]
    kwargs: dict = {
        "positions": positions,
        "colors": colors,
        "errors": errors,
        "track_image_indexes": new_timg,
        "track_feature_indexes": new_tfid,
        "track_point_ids": new_tpid,
    }
    if normals is not None:
        kwargs["normals"] = normals
    return recon.clone_with_changes(**kwargs)


def measure_dataset(
    dataset_label: str,
    sfmr_path: Path,
    *,
    subsample_stride: int = 1,
) -> dict:
    print(f"\n=== {dataset_label} ({sfmr_path.name}) ===", flush=True)
    recon, images = load_recon_and_images(sfmr_path)
    print(
        f"  images={recon.image_count} input_points={recon.point_count}",
        flush=True,
    )
    if subsample_stride > 1:
        recon = subsampled_recon(recon, subsample_stride)
        print(
            f"  sub-sampled (stride={subsample_stride}): {recon.point_count} points",
            flush=True,
        )

    # Snapshot the SIFT seed keypoints once — the reference for every variant's
    # `shift_vs_sift` metric. The conversion is identical to the one
    # `embed_patches` runs as step 0, so each variant's per-observation
    # keypoint shares its starting position with this map.
    sift_seeds = sift_seed_keypoints(recon, patch_size=10.0)
    print(f"  SIFT seed observations: {len(sift_seeds)}", flush=True)

    skips = DATASET_VARIANT_SKIPS.get(dataset_label, set())

    rows: list[dict] = []
    baseline_obs: dict[tuple[bytes, int], np.ndarray] | None = None
    for label, variant_def in VARIANTS:
        if label in skips:
            print(f"  - {label} ... skipped (dataset opt-out)", flush=True)
            rows.append(
                {
                    "variant": label,
                    "skipped": "dataset opt-out",
                    "kwargs": variant_kwargs_for_record(variant_def),
                }
            )
            continue
        print(f"  - {label} ... ", end="", flush=True)
        try:
            out, dt = run_variant(recon, images, variant_def)
        except Exception as e:
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
                    "p95_shift_px": float("nan"),
                    "n_overlap": 0,
                    "mean_shift_vs_sift_px": float("nan"),
                    "median_shift_vs_sift_px": float("nan"),
                    "p95_shift_vs_sift_px": float("nan"),
                    "n_overlap_vs_sift": 0,
                    "kwargs": variant_kwargs_for_record(variant_def),
                    "error": str(e),
                }
            )
            continue
        agree = measure_photometric_agreement(out, images)
        variant_obs = per_obs_keypoints_by_world(out)
        if label == "baseline":
            baseline_obs = variant_obs
            shift_vs_baseline = {
                "mean_shift_px": 0.0,
                "median_shift_px": 0.0,
                "p95_shift_px": 0.0,
                "n_overlap": -1,
            }
        else:
            assert baseline_obs is not None, (
                "baseline must come first to define the join reference"
            )
            shift_vs_baseline = keypoint_shift_summary(baseline_obs, variant_obs)
        sift_shift = keypoint_shift_summary(sift_seeds, variant_obs)
        sift_shift_record = {
            "mean_shift_vs_sift_px": sift_shift["mean_shift_px"],
            "median_shift_vs_sift_px": sift_shift["median_shift_px"],
            "p95_shift_vs_sift_px": sift_shift["p95_shift_px"],
            "n_overlap_vs_sift": sift_shift["n_overlap"],
        }
        row = {
            "variant": label,
            "wall_secs": dt,
            "out_points": int(out.point_count),
            "mean_ecc": agree["mean_ecc"],
            "n_scored": agree["n_scored"],
            **shift_vs_baseline,
            **sift_shift_record,
            "kwargs": variant_kwargs_for_record(variant_def),
        }
        rows.append(row)
        print(
            f"{dt:6.2f}s  points={out.point_count:4d}  "
            f"mean_ecc={agree['mean_ecc']:.4f}  "
            f"shift_vs_base={shift_vs_baseline['mean_shift_px']:.4f}px  "
            f"shift_vs_sift={sift_shift['mean_shift_px']:.4f}px",
            flush=True,
        )
    return {
        "dataset": dataset_label,
        "sfmr_path": str(sfmr_path),
        "image_count": recon.image_count,
        "input_point_count": recon.point_count,
        "subsample_stride": subsample_stride,
        "sift_seed_obs_count": len(sift_seeds),
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
    p.add_argument(
        "--dino-stride",
        type=int,
        default=int(os.environ.get("SFMTOOL_DINO_SUBSAMPLE_STRIDE", "1")),
        help=(
            "Stride for dino_dog_toy point sub-sampling (default 1 = full set; "
            "use a larger N if the full sweep takes too long, e.g. 4 keeps every "
            "4th point in source order)"
        ),
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
        stride = args.dino_stride if d == "dino_dog_toy" else 1
        results.append(measure_dataset(d, sfmr, subsample_stride=stride))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(results, indent=2))
    print(f"\nWrote {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
