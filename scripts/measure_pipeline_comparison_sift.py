#!/usr/bin/env python3
# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Pipeline comparison: how much do four SfM refinement pipelines diverge from
the raw SIFT input, and from each other?

All four pipelines start by running ``recon.to_embedded_patches(...)``, which
copies the raw SIFT-detection keypoints inline and seeds each point's patch
frame from the mean-viewing direction. They then differ in what (if anything)
moves the keypoints, before always re-refining the normals to be consistent
with the final keypoint locations:

- **Pipeline A** ("normals-only"): ``to_embedded_patches`` → ``refine_normals``.
  Keypoints stay at SIFT positions; only the per-point normal is photometrically
  optimized.
- **Pipeline B** ("grid"): ``to_embedded_patches`` → ``refine_normals`` →
  ``localize_keypoints`` (the supersampled grid keypoint search) →
  ``refine_normals`` again, with the moved keypoints.
- **Pipeline C** ("LK-bilinear"): ``to_embedded_patches`` → ``refine_normals`` →
  ``refine_keypoints`` (the LK subpixel refiner, default ``sampler="bilinear"``)
  → ``refine_normals`` again.
- **Pipeline D** ("LK-anisotropic"): identical to C except the LK refiner is
  invoked with ``sampler="anisotropic"`` — the anti-aliased oblique-view
  sampler. The decision-gate report already measured that this sampler does
  not meaningfully change ECC scores; this run answers a different question:
  how much does it move the persisted keypoints, normals, and bitmaps?

For each (dataset × pipeline), three artifacts are captured:

1. Per-observation **keypoint** positions, source-image px — joined across
   pipelines by ``(rounded world-position bytes, image_index)``, the same join
   the existing :file:`measure_subpixel_decision_gate.py` uses.
2. Per-point **normals** (3D unit vectors), compared via the inter-normal angle
   in degrees: ``arccos(clip(dot(n1, n2), -1, 1))``.
3. Per-point **patch bitmaps** — the RGBA reference texture each pipeline
   persists via ``refine_normals(render_bitmaps=True)``. Comparison is a
   mean per-pixel L1 distance over the RGB channels (alpha ignored), plus the
   count of points whose mean-L1 exceeds a "substantially different" threshold
   (default 16/255).

Results are dumped to a JSON sidecar and (in the companion prose report)
rendered as per-dataset tables.

Run with::

    pixi run -e test python scripts/measure_pipeline_comparison_sift.py [DATASET ...]

When no datasets are named, every dataset is measured.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

# Add repo src/ so the script runs as a standalone python invocation.
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from sfmtool._embed_patches import compact_to_embedded_patches  # noqa: E402
from sfmtool._sfmtool import SfmrReconstruction  # noqa: E402
from sfmtool._workspace_image import read_workspace_image  # noqa: E402

# Reuse the join helper from the decision-gate script so the join key cannot
# drift across the two measurements.
from measure_subpixel_decision_gate import (  # noqa: E402
    per_obs_keypoints_by_world,
    subsampled_recon,
)

DATASETS: dict[str, str] = {
    "seoul_bull": "seoul_bull_ws/sfmr/20260621-00-solve-seoul_bull_sculpture_1-17.sfmr",
    "dino_dog_toy": "dino_dog_toy_ws/sfmr/20260621-00-solve-dino_dog_toy_1-85.sfmr",
    "seattle_backyard": "seattle_backyard_ws/sfmr/20260621-00-solve-seattle_backyard_1-26.sfmr",
    "kerry_park": "kerry_park_ws/sfmr/20260621-00-solve-frame_1-24.sfmr",
}

# Stride at which to sub-sample the dino_dog_toy points. dino is 85
# high-resolution images × ~19k points; with the full point set, three pipelines
# × two refine_normals calls each would blow the run-time budget. The same
# stride-5 sub-sample that measure_subpixel_decision_gate.py defaults to keeps
# the per-pipeline relative comparison intact while fitting comfortably.
DEFAULT_DINO_STRIDE = 5

# Pipeline knobs — match the defaults the production `embed_patches` / `xform
# refine-normals` paths use. Pick a single resolution for everything: it's the
# patch grid the refiner samples on AND the resolution of the rendered RGBA
# bitmap (every patch contributes an RxRx4 texture). 24 matches the production
# default; the storage cost per point is 24*24*4 = 2.3 KiB, which is fine.
RESOLUTION = 24

# `localize_keypoints` knobs that match `embed_patches`'s defaults (these
# inputs are the same ones the production `subpixel="none"` / "lk" /
# "lk_per_move" paths feed `localize_keypoints`).
LOCALIZE_KWARGS = dict(
    max_iters=5,
    search=6.0,
    max_shift_px=3.0,
    min_relative_zncc=0.7,
    # The "grid" pipeline is the supersampled-grid keypoint search. The
    # production CLI exposes this behind `--search-resolution-multiplier`; the
    # decision-gate report identified m=2 as the sweet spot (~5× baseline cost,
    # but a meaningful per-observation move vs m=1 which is degenerate at one
    # grid cell per image), so this run uses m=2.
    search_resolution_multiplier=2.0,
)

# `refine_keypoints` knobs that match the production `subpixel="lk"` path
# (per-sweep consensus, one outer sweep — the cheapest LK setting that still
# moves keypoints to sub-pixel).
REFINE_KEYPOINTS_KWARGS = dict(
    max_outer_sweeps=1,
    consensus_refresh="per_sweep",
)

# Bitmap-difference threshold: a point is counted as "substantially different"
# between two pipelines if its mean per-pixel L1 distance over the RGB channels
# (computed over the union of "either pipeline has any alpha here" pixels)
# exceeds this threshold. 16/255 ~= 6.3% is a fairly conservative cut — much
# smaller than e.g. JPEG-noise levels — so a point that crosses it really has a
# different reference appearance, not just a per-texel sub-pixel jitter.
BITMAP_DIFF_THRESHOLD_L1 = 16.0 / 255.0


# ---------------------------------------------------------------------------
# Pipeline drivers
# ---------------------------------------------------------------------------


def _load_recon_and_images(
    sfmr_path: Path,
) -> tuple[SfmrReconstruction, list[np.ndarray]]:
    recon = SfmrReconstruction.load(str(sfmr_path))
    images = [
        read_workspace_image(recon.workspace_dir, name) for name in recon.image_names
    ]
    return recon, images


def _to_embedded(
    recon: SfmrReconstruction, patch_size: float = 10.0
) -> SfmrReconstruction:
    """Step 0 — convert sift_files → embedded_patches with the production
    knobs the rest of the pipeline assumes (mean-viewing normal seed, frame
    sized by SIFT feature scale, inline-copied SIFT keypoints)."""
    half_extent = patch_size / 2.0
    return recon.to_embedded_patches(
        normal="mean_viewing", extent="feature_size", extent_value=half_extent
    )


def _refine_normals_on(
    embedded: SfmrReconstruction, images: list[np.ndarray]
) -> SfmrReconstruction:
    """Run a single ``refine_normals`` pass on an embedded_patches recon and
    return a new embedded_patches recon that carries the refined normals AND the
    rendered per-point RGBA bitmap (so the next pipeline step can read either
    back from the recon, just like the production path does).

    Reuses the production defaults (use_stored_keypoints=True — anchor every
    view on the recon's stored keypoint; render_bitmaps=True — persist the
    fused RGBA reference appearance per point).
    """
    cloud = embedded.patches
    if cloud is None:
        raise ValueError("embedded recon has no patch cloud to refine")
    result = cloud.refine_normals(
        embedded,
        images,
        resolution=RESOLUTION,
        use_stored_keypoints=True,
        render_bitmaps=True,
    )
    # Scatter the refined normals back to per-point rows (cloud.point_ids
    # indexes the recon's 3D points; finite points only — infinity points keep
    # their stored normal).
    refined = np.asarray(result["normal"], dtype=np.float32)
    point_ids = np.asarray(cloud.point_ids, dtype=np.int64)
    finite = ~np.asarray(embedded.point_is_at_infinity)[point_ids]
    normals = np.asarray(embedded.normals, dtype=np.float32).copy()
    finite_pids = point_ids[finite]
    normals[finite_pids] = refined[finite]
    bitmaps = np.asarray(result["bitmaps"], dtype=np.uint8)
    return embedded.clone_with_changes(
        normals=normals, patches=cloud, patch_bitmaps=bitmaps
    )


def _localize_and_recompact(
    embedded: SfmrReconstruction, images: list[np.ndarray]
) -> SfmrReconstruction:
    """Run ``localize_keypoints`` and roll the result into a new
    embedded_patches recon. Reuses ``compact_to_embedded_patches`` so the join
    semantics line up with the production pipeline; ``min_views=1`` keeps every
    point the localizer admitted, so the join-vs-SIFT overlap is maximal."""
    cloud = embedded.patches
    if cloud is None:
        raise ValueError("embedded recon has no patch cloud to localize")
    # Step 2 (view selection) — re-run before localize so the view set is built
    # to the same shape `embed_patches` uses.
    selections = cloud.select_views(
        embedded,
        images,
        min_relative_zncc=LOCALIZE_KWARGS["min_relative_zncc"],
        resolution=RESOLUTION,
    )
    view_sets = {
        int(s["point_id"]): np.asarray(s["admitted"]).tolist() for s in selections
    }
    localizations = cloud.localize_keypoints(
        embedded,
        images,
        view_sets=view_sets,
        resolution=RESOLUTION,
        **LOCALIZE_KWARGS,
    )
    return compact_to_embedded_patches(
        embedded,
        cloud,
        localizations,
        embedded.image_file_hashes,
        # Carry the just-rendered bitmaps through compaction, then they'll be
        # overwritten by the second refine_normals call.
        patch_bitmaps=np.asarray(embedded.patch_bitmaps, dtype=np.uint8)
        if embedded.patch_bitmaps is not None
        else None,
        min_views=1,
    )


def _refine_keypoints_and_recompact(
    embedded: SfmrReconstruction,
    images: list[np.ndarray],
    *,
    sampler: str = "bilinear",
) -> SfmrReconstruction:
    """LK counterpart of ``_localize_and_recompact``: build per-point view sets
    + per-view seeds from the recon's stored keypoints, run ``refine_keypoints``
    seeded there, and compact the result into a new embedded_patches recon.

    ``sampler`` is passed straight through to ``refine_keypoints`` —
    ``"bilinear"`` (default; production ``subpixel="lk"`` path) or
    ``"anisotropic"`` (anti-aliased oblique-view sampler; pipeline D).
    """
    cloud = embedded.patches
    if cloud is None:
        raise ValueError("embedded recon has no patch cloud for LK refinement")
    # Use the stored keypoints as the seed: build per-point {views, keypoints}
    # from `track_*` arrays (the same logic `_refine_subpixel` uses on the
    # localizer's output).
    tpid = np.asarray(embedded.track_point_ids, dtype=np.int64)
    timg = np.asarray(embedded.track_image_indexes, dtype=np.int64)
    kxy = np.asarray(embedded.keypoints_xy, dtype=np.float64)
    cloud_pids = set(int(p) for p in cloud.point_ids)
    view_sets: dict[int, list[int]] = {}
    seeds: dict[int, list[list[float]]] = {}
    for k, (p, i) in enumerate(zip(tpid.tolist(), timg.tolist())):
        if p not in cloud_pids:
            continue
        view_sets.setdefault(p, []).append(i)
        seeds.setdefault(p, []).append([float(kxy[k, 0]), float(kxy[k, 1])])
    refined = cloud.refine_keypoints(
        embedded,
        images,
        view_sets=view_sets,
        starting_keypoints=seeds,
        point_ids=list(view_sets.keys()),
        resolution=RESOLUTION,
        sampler=sampler,
        **REFINE_KEYPOINTS_KWARGS,
    )
    # Splice refined keypoints into a localizations-shaped list. The refiner
    # returns the same view set it was handed (per-view membership is
    # preserved).
    localizations: list[dict[str, Any]] = []
    refined_by_pid = {int(r["point_id"]): r for r in refined}
    for pid, views in view_sets.items():
        r = refined_by_pid.get(pid)
        if r is None:
            # Fell out of the refiner — fall back to the seed keypoints so the
            # observation is still emitted by compaction (this matches the
            # behavior `_refine_subpixel` falls back to).
            localizations.append(
                {
                    "point_id": pid,
                    "views": views,
                    "keypoints": np.asarray(seeds[pid], dtype=np.float64),
                }
            )
            continue
        localizations.append(
            {
                "point_id": pid,
                "views": np.asarray(r["views"], dtype=np.uint32).tolist(),
                "keypoints": np.asarray(r["keypoints"], dtype=np.float64).reshape(
                    -1, 2
                ),
            }
        )
    return compact_to_embedded_patches(
        embedded,
        cloud,
        localizations,
        embedded.image_file_hashes,
        patch_bitmaps=np.asarray(embedded.patch_bitmaps, dtype=np.uint8)
        if embedded.patch_bitmaps is not None
        else None,
        min_views=1,
    )


def run_pipeline(
    pipeline: str,
    recon: SfmrReconstruction,
    images: list[np.ndarray],
) -> tuple[SfmrReconstruction, float]:
    """Run one pipeline end to end. Returns the final embedded_patches recon
    (carrying stored keypoints, refined normals, and rendered bitmaps) and the
    wall time."""
    t0 = time.perf_counter()
    embedded = _to_embedded(recon)
    # Step 1 — always refine normals once over the SIFT-seeded keypoints.
    embedded = _refine_normals_on(embedded, images)
    if pipeline == "A":
        out = embedded
    elif pipeline == "B":
        # Move keypoints via the grid localizer, then refine normals again
        # against the moved keypoints.
        moved = _localize_and_recompact(embedded, images)
        out = _refine_normals_on(moved, images)
    elif pipeline == "C":
        moved = _refine_keypoints_and_recompact(embedded, images, sampler="bilinear")
        out = _refine_normals_on(moved, images)
    elif pipeline == "D":
        moved = _refine_keypoints_and_recompact(embedded, images, sampler="anisotropic")
        out = _refine_normals_on(moved, images)
    else:
        raise ValueError(f"unknown pipeline {pipeline!r}")
    dt = time.perf_counter() - t0
    return out, dt


# ---------------------------------------------------------------------------
# Comparison helpers (keypoints, normals, bitmaps)
# ---------------------------------------------------------------------------


def keypoint_shift_summary(
    reference: dict[tuple[bytes, int], np.ndarray],
    variant: dict[tuple[bytes, int], np.ndarray],
) -> dict[str, float]:
    """Mean / median / p95 per-observation keypoint shift in source-image px,
    over the observations present in BOTH dicts."""
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


def normals_by_world(recon: SfmrReconstruction) -> dict[bytes, np.ndarray]:
    """Map ``rounded world position bytes -> normal (3,)`` over the recon's
    points. The same world-position rounding the keypoint join uses, so the
    two comparisons line up on identical join keys."""
    positions = np.asarray(recon.positions, dtype=np.float64)
    normals = np.asarray(recon.normals, dtype=np.float64)
    out: dict[bytes, np.ndarray] = {}
    for i in range(positions.shape[0]):
        key = np.round(positions[i], 6).tobytes()
        out[key] = normals[i]
    return out


def normal_angle_summary(
    reference: dict[bytes, np.ndarray],
    variant: dict[bytes, np.ndarray],
) -> dict[str, float]:
    """Mean / median / p95 inter-normal angle in degrees, over the points
    present in both dicts and whose normals are finite, non-zero unit vectors
    on both sides. Points with degenerate normals (the all-zero seed reserved
    for points at infinity, or a NaN row produced by a refine_normals failure)
    are excluded — they would otherwise dominate the percentiles."""
    shared = reference.keys() & variant.keys()
    if not shared:
        return {
            "mean_angle_deg": float("nan"),
            "median_angle_deg": float("nan"),
            "p95_angle_deg": float("nan"),
            "n_overlap": 0,
            "n_skipped": 0,
        }
    a = np.array([reference[k] for k in shared], dtype=np.float64)
    b = np.array([variant[k] for k in shared], dtype=np.float64)
    na = np.linalg.norm(a, axis=1)
    nb = np.linalg.norm(b, axis=1)
    valid = (
        np.isfinite(a).all(axis=1)
        & np.isfinite(b).all(axis=1)
        & (na > 1e-6)
        & (nb > 1e-6)
    )
    n_skipped = int((~valid).sum())
    if not valid.any():
        return {
            "mean_angle_deg": float("nan"),
            "median_angle_deg": float("nan"),
            "p95_angle_deg": float("nan"),
            "n_overlap": 0,
            "n_skipped": n_skipped,
        }
    a_unit = a[valid] / na[valid, None]
    b_unit = b[valid] / nb[valid, None]
    dots = np.clip((a_unit * b_unit).sum(axis=1), -1.0, 1.0)
    angles_deg = np.degrees(np.arccos(dots))
    return {
        "mean_angle_deg": float(angles_deg.mean()),
        "median_angle_deg": float(np.median(angles_deg)),
        "p95_angle_deg": float(np.percentile(angles_deg, 95)),
        "n_overlap": int(valid.sum()),
        "n_skipped": n_skipped,
    }


def bitmaps_by_world(recon: SfmrReconstruction) -> dict[bytes, np.ndarray]:
    """Map ``rounded world position bytes -> RGBA bitmap (R, R, 4) uint8``
    over the recon's points. Same join key as the keypoint / normal maps."""
    positions = np.asarray(recon.positions, dtype=np.float64)
    bitmaps = recon.patch_bitmaps
    if bitmaps is None:
        return {}
    bitmaps = np.asarray(bitmaps, dtype=np.uint8)
    if bitmaps.shape[0] != positions.shape[0]:
        # Defensive — the scatter is per-3D-point and should always be parallel
        # to `positions`. If it isn't, something is upstream-wrong.
        raise RuntimeError(
            f"patch_bitmaps shape {bitmaps.shape} not parallel to {positions.shape[0]} points"
        )
    out: dict[bytes, np.ndarray] = {}
    for i in range(positions.shape[0]):
        key = np.round(positions[i], 6).tobytes()
        out[key] = bitmaps[i]
    return out


def _bitmap_sharpness_per_point(
    bitmap_dict: dict[bytes, np.ndarray],
) -> list[tuple[float, float]]:
    """Per-point ``(laplacian_var, gradient_mag_mean)`` over each bitmap's
    luminance channel. Caller-side aggregate convenience: returns just the
    survivors with no keys. See ``_bitmap_sharpness_by_world`` for the keyed
    variant used to pair across pipelines.

    Per-bitmap recipe is documented on ``_bitmap_sharpness_by_world``; this
    function shares it.
    """
    return [v for v in _bitmap_sharpness_by_world(bitmap_dict).values()]


def _bitmap_sharpness_by_world(
    bitmap_dict: dict[bytes, np.ndarray],
) -> dict[bytes, tuple[float, float]]:
    """Per-point ``world_key -> (laplacian_var, gradient_mag_mean)`` over each
    bitmap's luminance channel. Same join key as
    :func:`bitmaps_by_world` / :func:`normals_by_world`, so the result can be
    intersected across pipelines for paired per-point sharpness deltas.

    Inputs:
      ``bitmap_dict`` — ``{world_key: RGBA uint8 (H, W, 4)}`` (the dict
      returned by :func:`bitmaps_by_world`).

    Per bitmap:
      * Luminance ``I_gray = mean(rgb, axis=-1)`` in ``float32`` normalized to
        ``[0, 1]`` (so the metric is on the same scale as `bitmap_l1_summary`'s
        normalized output).
      * Alpha mask matches `bitmap_l1_summary` semantics: we use ``alpha > 0``
        per bitmap so all-transparent border texels don't deflate sharpness.
      * Laplacian: 3x3 kernel ``[[0,1,0],[1,-4,1],[0,1,0]]`` applied via
        numpy slicing with ``'reflect'`` boundary handling (one-row/column
        reflection — the standard symmetric padding so interior texels aren't
        biased by zero borders). Variance is taken over masked texels only.
      * Gradient magnitude: forward differences ``dx = diff(I, axis=1)``,
        ``dy = diff(I, axis=0)``, padded to original shape with zeros on the
        trailing column/row. ``mag = sqrt(dx**2 + dy**2)``. Mean is taken over
        masked texels (same mask).
      * If the masked region is empty OR has < 2 masked texels (the Laplacian
        variance is undefined for n<2) the bitmap is skipped (key omitted from
        the returned dict, NOT NaN — caller aggregates / intersects survivors).
    """
    out: dict[bytes, tuple[float, float]] = {}
    for key, rgba in bitmap_dict.items():
        if rgba.ndim != 3 or rgba.shape[-1] != 4:
            continue
        alpha = rgba[:, :, 3]
        mask = alpha > 0
        if mask.sum() < 2:
            continue
        rgb = rgba[:, :, :3].astype(np.float32) / 255.0
        gray = rgb.mean(axis=2)  # (H, W) in [0, 1]
        # 3x3 Laplacian [[0,1,0],[1,-4,1],[0,1,0]] via numpy with reflect
        # padding. np.pad(..., 'reflect') reflects across the edge sample
        # without duplicating it, matching the standard image-processing
        # convention used by scipy.ndimage.convolve(mode='reflect').
        padded = np.pad(gray, 1, mode="reflect")
        lap = (
            padded[:-2, 1:-1]
            + padded[2:, 1:-1]
            + padded[1:-1, :-2]
            + padded[1:-1, 2:]
            - 4.0 * padded[1:-1, 1:-1]
        )
        lap_masked = lap[mask]
        lap_var = float(lap_masked.var())
        # Forward differences padded with zeros to original shape.
        dx = np.zeros_like(gray)
        dy = np.zeros_like(gray)
        dx[:, :-1] = np.diff(gray, axis=1)
        dy[:-1, :] = np.diff(gray, axis=0)
        mag = np.sqrt(dx * dx + dy * dy)
        mag_mean = float(mag[mask].mean())
        out[key] = (lap_var, mag_mean)
    return out


def bitmap_sharpness_summary(
    bitmap_dict: dict[bytes, np.ndarray],
) -> dict[str, float | int]:
    """Aggregate per-pipeline sharpness from ``_bitmap_sharpness_per_point``.

    Returns mean / median / p95 of Laplacian variance and of gradient-mag mean
    across the bitmaps that survived the alpha-mask check, plus ``n_compared``
    (survivors) and ``n_skipped`` (bitmaps with <2 masked texels — fully
    transparent or single-pixel coverage)."""
    per_point = _bitmap_sharpness_per_point(bitmap_dict)
    n_total = len(bitmap_dict)
    n_compared = len(per_point)
    n_skipped = n_total - n_compared
    if not per_point:
        return {
            "lap_var_mean": float("nan"),
            "lap_var_median": float("nan"),
            "lap_var_p95": float("nan"),
            "grad_mag_mean_mean": float("nan"),
            "grad_mag_mean_median": float("nan"),
            "grad_mag_mean_p95": float("nan"),
            "n_compared": n_compared,
            "n_skipped": n_skipped,
        }
    arr = np.array(per_point, dtype=np.float64)
    lap = arr[:, 0]
    grad = arr[:, 1]
    return {
        "lap_var_mean": float(lap.mean()),
        "lap_var_median": float(np.median(lap)),
        "lap_var_p95": float(np.percentile(lap, 95)),
        "grad_mag_mean_mean": float(grad.mean()),
        "grad_mag_mean_median": float(np.median(grad)),
        "grad_mag_mean_p95": float(np.percentile(grad, 95)),
        "n_compared": n_compared,
        "n_skipped": n_skipped,
    }


def common_subset_sharpness(
    sharpness_by_pipeline: dict[str, dict[bytes, tuple[float, float]]],
    *,
    pairs: list[tuple[str, str]],
) -> dict[str, Any]:
    """Apples-to-apples sharpness aggregates restricted to the per-dataset
    intersection of points covered by every pipeline.

    Inputs:
      ``sharpness_by_pipeline`` — ``{pipeline_label: {world_key: (lap_var,
        grad_mag)}}`` as returned by :func:`_bitmap_sharpness_by_world` per
        pipeline. Only pipelines present here (i.e. those whose sharpness was
        successfully computed) participate in the intersection.
      ``pairs`` — list of ``(X, Y)`` pipeline labels (e.g. ``("B", "A")``);
        the function emits per-pair paired-delta aggregates for each. Pairs
        whose either side is missing from ``sharpness_by_pipeline`` are
        skipped (so the call site can pass the full 6-pair list and have
        partial-pipeline runs Just Work).

    Returns a dict with three keys:
      * ``n_common`` — size of the intersection (0 if any pipeline has no
        bitmaps or the intersection is empty).
      * ``per_pipeline`` — ``{pipeline_label: {lap_var_*, grad_mag_mean_*,
        n_common}}`` with mean/median/p95 of each metric restricted to the
        intersection, in the same shape as :func:`bitmap_sharpness_summary`'s
        return (minus ``n_skipped`` — every pipeline has exactly ``n_common``
        survivors on the intersection by construction).
      * ``paired_deltas`` — ``{"X_vs_Y": {delta_lap_var_*,
        delta_grad_mag_mean_*, n_common}}`` for each pair, computed as
        ``metric_X - metric_Y`` per point on the intersection. Paired deltas
        eliminate per-point variance, so they give a cleaner signal than
        comparing the per-pipeline aggregates.
    """
    pipelines_present = [p for p, d in sharpness_by_pipeline.items() if d]
    if not pipelines_present:
        return {"n_common": 0, "per_pipeline": {}, "paired_deltas": {}}
    common: set[bytes] = set(sharpness_by_pipeline[pipelines_present[0]].keys())
    for p in pipelines_present[1:]:
        common &= sharpness_by_pipeline[p].keys()
    n_common = len(common)
    if n_common == 0:
        return {
            "n_common": 0,
            "per_pipeline": {p: None for p in pipelines_present},
            "paired_deltas": {},
        }
    # Stable key order so paired-delta indexing is deterministic across
    # pipelines (numpy arrays are zipped index-wise).
    keys = sorted(common)
    per_pipeline_arrays: dict[str, np.ndarray] = {}
    per_pipeline: dict[str, dict[str, float | int]] = {}
    for p in pipelines_present:
        arr = np.array([sharpness_by_pipeline[p][k] for k in keys], dtype=np.float64)
        per_pipeline_arrays[p] = arr
        lap = arr[:, 0]
        grad = arr[:, 1]
        per_pipeline[p] = {
            "lap_var_mean": float(lap.mean()),
            "lap_var_median": float(np.median(lap)),
            "lap_var_p95": float(np.percentile(lap, 95)),
            "grad_mag_mean_mean": float(grad.mean()),
            "grad_mag_mean_median": float(np.median(grad)),
            "grad_mag_mean_p95": float(np.percentile(grad, 95)),
            "n_common": n_common,
        }
    paired_deltas: dict[str, dict[str, float | int]] = {}
    for x, y in pairs:
        if x not in per_pipeline_arrays or y not in per_pipeline_arrays:
            continue
        dlap = per_pipeline_arrays[x][:, 0] - per_pipeline_arrays[y][:, 0]
        dgrad = per_pipeline_arrays[x][:, 1] - per_pipeline_arrays[y][:, 1]
        paired_deltas[f"{x}_vs_{y}"] = {
            "delta_lap_var_mean": float(dlap.mean()),
            "delta_lap_var_median": float(np.median(dlap)),
            "delta_lap_var_p95": float(np.percentile(dlap, 95)),
            "delta_grad_mag_mean_mean": float(dgrad.mean()),
            "delta_grad_mag_mean_median": float(np.median(dgrad)),
            "delta_grad_mag_mean_p95": float(np.percentile(dgrad, 95)),
            "n_common": n_common,
        }
    return {
        "n_common": n_common,
        "per_pipeline": per_pipeline,
        "paired_deltas": paired_deltas,
    }


def bitmap_l1_summary(
    reference: dict[bytes, np.ndarray],
    variant: dict[bytes, np.ndarray],
    *,
    diff_threshold: float = BITMAP_DIFF_THRESHOLD_L1,
) -> dict[str, float]:
    """Per-point mean per-pixel L1 distance over the RGB channels (alpha
    ignored; mean is taken over the union of texels where either pipeline has
    any alpha — i.e. either bitmap is non-zero). Points with no covered texel
    on either side are skipped; points where every texel is zero in BOTH
    pipelines (an unrefined / no-bitmap point) are also skipped.

    Returns mean / median / p95 of per-point mean-L1 (in normalized [0, 1]
    units), the count of points whose per-point mean-L1 exceeds
    ``diff_threshold``, and the population sizes."""
    shared = reference.keys() & variant.keys()
    per_point_mean_l1: list[float] = []
    n_skipped = 0
    for key in shared:
        a = reference[key]
        b = variant[key]
        if a.shape != b.shape:
            n_skipped += 1
            continue
        # Mask = either pipeline rendered to this texel. Defined as
        # `alpha > 0` on either side.
        a_alpha = a[:, :, 3]
        b_alpha = b[:, :, 3]
        mask = (a_alpha > 0) | (b_alpha > 0)
        if not mask.any():
            n_skipped += 1
            continue
        a_rgb = a[:, :, :3].astype(np.float32)
        b_rgb = b[:, :, :3].astype(np.float32)
        # Mean per-pixel L1 over RGB, taken across the covered texels.
        diff = np.abs(a_rgb - b_rgb).mean(axis=2)  # (R, R)
        masked_mean = float(diff[mask].mean()) / 255.0
        per_point_mean_l1.append(masked_mean)
    if not per_point_mean_l1:
        return {
            "mean_l1": float("nan"),
            "median_l1": float("nan"),
            "p95_l1": float("nan"),
            "n_substantially_different": 0,
            "diff_threshold_l1": diff_threshold,
            "n_compared": 0,
            "n_skipped": n_skipped,
        }
    arr = np.array(per_point_mean_l1, dtype=np.float64)
    return {
        "mean_l1": float(arr.mean()),
        "median_l1": float(np.median(arr)),
        "p95_l1": float(np.percentile(arr, 95)),
        "n_substantially_different": int((arr > diff_threshold).sum()),
        "diff_threshold_l1": diff_threshold,
        "n_compared": int(arr.size),
        "n_skipped": n_skipped,
    }


# ---------------------------------------------------------------------------
# Per-dataset driver
# ---------------------------------------------------------------------------


PIPELINES = ("A", "B", "C", "D")
PIPELINE_LABELS = {
    "A": "normals-only",
    "B": "grid (localize_keypoints)",
    "C": "LK (refine_keypoints, bilinear)",
    "D": "LK (refine_keypoints, anisotropic)",
}


def measure_dataset(
    dataset_label: str,
    sfmr_path: Path,
    *,
    subsample_stride: int = 1,
    pipelines: tuple[str, ...] = PIPELINES,
) -> dict:
    print(f"\n=== {dataset_label} ({sfmr_path.name}) ===", flush=True)
    recon, images = _load_recon_and_images(sfmr_path)
    print(
        f"  images={recon.image_count} input_points={recon.point_count}",
        flush=True,
    )
    full_point_count = int(recon.point_count)
    if subsample_stride > 1:
        recon = subsampled_recon(recon, subsample_stride)
        print(
            f"  sub-sampled (stride={subsample_stride}): {recon.point_count} points",
            flush=True,
        )

    # SIFT baseline: just `to_embedded_patches`. The recon has the raw inline
    # SIFT keypoints and the mean-viewing seed normals; we use it as the
    # reference for keypoint and normal comparisons. (No bitmaps yet — those
    # only get rendered by `refine_normals(render_bitmaps=True)`. The
    # `*_vs_SIFT` bitmap row is therefore reported as N/A and only the three
    # pairwise comparisons are populated for bitmaps.)
    sift = _to_embedded(recon)
    sift_keypoints = per_obs_keypoints_by_world(sift)
    sift_normals = normals_by_world(sift)
    print(
        f"  SIFT seeds: {len(sift_keypoints)} observations, {len(sift_normals)} points",
        flush=True,
    )

    pipeline_outputs: dict[str, dict] = {}
    for pipeline in pipelines:
        print(
            f"  - pipeline {pipeline} ({PIPELINE_LABELS[pipeline]}) ... ",
            end="",
            flush=True,
        )
        try:
            out, dt = run_pipeline(pipeline, recon, images)
        except Exception as e:
            print(f"FAILED: {e}", flush=True)
            pipeline_outputs[pipeline] = {"error": str(e), "wall_secs": float("nan")}
            continue
        pipeline_outputs[pipeline] = {
            "wall_secs": dt,
            "out_points": int(out.point_count),
            "keypoints": per_obs_keypoints_by_world(out),
            "normals": normals_by_world(out),
            "bitmaps": bitmaps_by_world(out),
        }
        print(
            f"{dt:6.2f}s  points={out.point_count:4d}  "
            f"obs={len(pipeline_outputs[pipeline]['keypoints'])}",
            flush=True,
        )

    # Build the pair comparisons. Order kept stable across A/B/C from the
    # original report; D-vs-* rows appended at the end so JSON readers that
    # index by position on the legacy rows are unaffected.
    pairs = [
        ("A_vs_SIFT", "A", "_SIFT"),
        ("B_vs_SIFT", "B", "_SIFT"),
        ("C_vs_SIFT", "C", "_SIFT"),
        ("B_vs_A", "B", "A"),
        ("C_vs_A", "C", "A"),
        ("C_vs_B", "C", "B"),
        ("D_vs_SIFT", "D", "_SIFT"),
        ("D_vs_A", "D", "A"),
        ("D_vs_B", "D", "B"),
        ("D_vs_C", "D", "C"),
    ]

    def _pipe_kp(p: str) -> dict[tuple[bytes, int], np.ndarray]:
        if p == "_SIFT":
            return sift_keypoints
        return pipeline_outputs[p].get("keypoints", {})

    def _pipe_n(p: str) -> dict[bytes, np.ndarray]:
        if p == "_SIFT":
            return sift_normals
        return pipeline_outputs[p].get("normals", {})

    def _pipe_b(p: str) -> dict[bytes, np.ndarray]:
        if p == "_SIFT":
            return {}  # no bitmaps for SIFT baseline
        return pipeline_outputs[p].get("bitmaps", {})

    keypoint_rows: list[dict] = []
    normal_rows: list[dict] = []
    bitmap_rows: list[dict] = []
    for label, lhs, rhs in pairs:
        kp = keypoint_shift_summary(_pipe_kp(rhs), _pipe_kp(lhs))
        n = normal_angle_summary(_pipe_n(rhs), _pipe_n(lhs))
        keypoint_rows.append({"pair": label, **kp})
        normal_rows.append({"pair": label, **n})
        if rhs == "_SIFT":
            # No SIFT-baseline bitmap — report a sentinel row so the JSON shape
            # stays uniform across pairs (and the report can mark it N/A).
            bitmap_rows.append(
                {
                    "pair": label,
                    "mean_l1": float("nan"),
                    "median_l1": float("nan"),
                    "p95_l1": float("nan"),
                    "n_substantially_different": 0,
                    "diff_threshold_l1": BITMAP_DIFF_THRESHOLD_L1,
                    "n_compared": 0,
                    "n_skipped": 0,
                    "note": "no SIFT-baseline bitmap (rendering only happens during refine_normals)",
                }
            )
        else:
            bitmap_rows.append(
                {"pair": label, **bitmap_l1_summary(_pipe_b(rhs), _pipe_b(lhs))}
            )

    # Per-pipeline bitmap sharpness (Laplacian variance + gradient-mag mean
    # over the rendered RGBA bitmaps, aggregated to mean/median/p95 per
    # pipeline per dataset). Computed against the same alpha>0 mask
    # `bitmap_l1_summary` uses so sharpness isn't deflated by transparent
    # border texels. Bitmaps were already rendered in memory — no re-render.
    sharpness_by_pipeline: dict[str, dict] = {}
    sharpness_per_point: dict[str, dict[bytes, tuple[float, float]]] = {}
    for p in pipelines:
        info = pipeline_outputs.get(p, {})
        if "error" in info:
            continue
        sharpness_by_pipeline[p] = bitmap_sharpness_summary(info.get("bitmaps", {}))
        sharpness_per_point[p] = _bitmap_sharpness_by_world(info.get("bitmaps", {}))

    # Apples-to-apples sharpness: restrict to the per-dataset intersection of
    # points covered by all pipelines (vs the per-pipeline `sharpness` view
    # above, where A/B/C/D each aggregate over different populations because
    # `alpha > 0` coverage differs per pipeline — B in particular admits extra
    # views via `localize_keypoints` and so has more covered points). The
    # restricted view emits per-pipeline aggregates on the same `n_common` set
    # AND paired per-point deltas (six pairs, mean/median/p95). Paired deltas
    # eliminate per-point variance and give a cleaner blurriness signal than
    # comparing the unrestricted aggregates.
    sharpness_common_subset = common_subset_sharpness(
        sharpness_per_point,
        pairs=[
            ("B", "A"),
            ("C", "A"),
            ("D", "A"),
            ("C", "B"),
            ("D", "B"),
            ("D", "C"),
        ],
    )

    # Strip the heavy in-memory dicts before returning so we can serialize.
    pipeline_summary = {}
    for p, info in pipeline_outputs.items():
        if "error" in info:
            pipeline_summary[p] = {
                "label": PIPELINE_LABELS[p],
                "error": info["error"],
                "wall_secs": info["wall_secs"],
            }
        else:
            pipeline_summary[p] = {
                "label": PIPELINE_LABELS[p],
                "wall_secs": info["wall_secs"],
                "out_points": info["out_points"],
                "n_observations": len(info["keypoints"]),
            }

    return {
        "dataset": dataset_label,
        "sfmr_path": str(sfmr_path),
        "image_count": recon.image_count,
        "input_point_count": full_point_count,
        "measured_point_count": int(recon.point_count),
        "subsample_stride": subsample_stride,
        "sift_seed_obs_count": len(sift_keypoints),
        "sift_seed_point_count": len(sift_normals),
        "resolution": RESOLUTION,
        "bitmap_diff_threshold_l1": BITMAP_DIFF_THRESHOLD_L1,
        "pipelines": pipeline_summary,
        "keypoint_comparisons": keypoint_rows,
        "normal_comparisons": normal_rows,
        "bitmap_comparisons": bitmap_rows,
        "sharpness": sharpness_by_pipeline,
        "sharpness_common_subset": sharpness_common_subset,
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
        default=ROOT / "reports" / "2026-06-27-pipeline-comparison-sift-data.json",
        help="Where to write the structured measurement results",
    )
    p.add_argument(
        "--dino-stride",
        type=int,
        default=DEFAULT_DINO_STRIDE,
        help=(
            "Stride for dino_dog_toy point sub-sampling "
            f"(default {DEFAULT_DINO_STRIDE}; set to 1 for the full point set)"
        ),
    )
    p.add_argument(
        "--pipelines",
        type=str,
        default=",".join(PIPELINES),
        help=(
            "Comma-separated subset of pipelines to run (default: all four — "
            f"{','.join(PIPELINES)}). Pairs involving a non-run pipeline are "
            "still emitted but with NaN summary stats and 0 overlap, since the "
            "comparison needs both sides' in-memory keypoint/normal/bitmap maps."
        ),
    )
    args = p.parse_args()

    pipelines_subset: tuple[str, ...] = tuple(
        x.strip() for x in args.pipelines.split(",") if x.strip()
    )
    unknown = [p for p in pipelines_subset if p not in PIPELINES]
    if unknown:
        raise SystemExit(
            f"unknown pipeline(s) {unknown!r}; choose from {list(PIPELINES)}"
        )

    selected = args.datasets or list(DATASETS.keys())
    results = []
    for d in selected:
        sfmr = ROOT / DATASETS[d]
        if not sfmr.exists():
            print(f"  SKIP {d}: {sfmr} missing", flush=True)
            continue
        stride = args.dino_stride if d == "dino_dog_toy" else 1
        results.append(
            measure_dataset(
                d, sfmr, subsample_stride=stride, pipelines=pipelines_subset
            )
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(results, indent=2))
    print(f"\nWrote {args.out}", flush=True)

    # Brief stdout summary.
    print("\n=== Cross-dataset summary ===", flush=True)
    for r in results:
        print(f"\n  {r['dataset']} (stride={r['subsample_stride']}):", flush=True)
        for kp in r["keypoint_comparisons"]:
            print(
                f"    kp {kp['pair']:>11s}: mean={kp['mean_shift_px']:.3f}px  "
                f"median={kp['median_shift_px']:.3f}px  p95={kp['p95_shift_px']:.3f}px  "
                f"n={kp['n_overlap']}",
                flush=True,
            )
        for n in r["normal_comparisons"]:
            print(
                f"    n  {n['pair']:>11s}: mean={n['mean_angle_deg']:.3f}deg  "
                f"median={n['median_angle_deg']:.3f}deg  p95={n['p95_angle_deg']:.3f}deg  "
                f"n={n['n_overlap']}",
                flush=True,
            )
        for p, s in r.get("sharpness", {}).items():
            print(
                f"    sh {p:>11s}: lap_var mean={s['lap_var_mean']:.4f} "
                f"median={s['lap_var_median']:.4f} p95={s['lap_var_p95']:.4f}  "
                f"grad_mag mean={s['grad_mag_mean_mean']:.4f} "
                f"n={s['n_compared']} skipped={s['n_skipped']}",
                flush=True,
            )
        scs = r.get("sharpness_common_subset", {})
        if scs.get("n_common"):
            print(f"    common-subset n_common={scs['n_common']}", flush=True)
            for p, s in scs.get("per_pipeline", {}).items():
                if s is None:
                    continue
                print(
                    f"    sh* {p:>10s}: lap_var mean={s['lap_var_mean']:.4f} "
                    f"median={s['lap_var_median']:.4f} p95={s['lap_var_p95']:.4f}",
                    flush=True,
                )
            for pair, d in scs.get("paired_deltas", {}).items():
                print(
                    f"    dsh {pair:>10s}: dlap mean={d['delta_lap_var_mean']:+.4f} "
                    f"median={d['delta_lap_var_median']:+.4f} p95={d['delta_lap_var_p95']:+.4f}  "
                    f"dgrad mean={d['delta_grad_mag_mean_mean']:+.4f}",
                    flush=True,
                )
        for b in r["bitmap_comparisons"]:
            if "note" in b:
                print(f"    bm {b['pair']:>11s}: (N/A) {b['note']}", flush=True)
            else:
                print(
                    f"    bm {b['pair']:>11s}: mean={b['mean_l1']:.4f}  "
                    f"median={b['median_l1']:.4f}  p95={b['p95_l1']:.4f}  "
                    f"n_diff={b['n_substantially_different']}/{b['n_compared']}",
                    flush=True,
                )

    return 0


if __name__ == "__main__":
    sys.exit(main())
