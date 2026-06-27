# Pipeline Comparison: SIFT Input vs Three Refinement Pipelines

_Date: 2026-06-27._

Three SfM refinement pipelines, all starting from the raw SIFT-detection
keypoints copied inline by `recon.to_embedded_patches(...)`, compared on
the four checked-in datasets:

- **Pipeline A** (normals-only): `to_embedded_patches` → `refine_normals`.
  Keypoints stay at SIFT positions; only the per-point normal is
  photometrically optimized.
- **Pipeline B** (grid): `to_embedded_patches` → `refine_normals` →
  `localize_keypoints` (supersampled grid, `search_resolution_multiplier=2`) →
  `refine_normals` (second time, against the moved keypoints).
- **Pipeline C** (LK): `to_embedded_patches` → `refine_normals` →
  `refine_keypoints` (LK, `subpixel='lk'` defaults) → `refine_normals`
  (second time, against the moved keypoints).

Three artifacts captured per (dataset × pipeline):

1. **Per-observation keypoints** in source-image px, joined across pipelines
   by `(rounded world position bytes, image_index)`. Compared by per-obs L2
   distance; mean / median / p95 reported.
2. **Per-point normals**, compared via the angle between corresponding unit
   vectors: `acos(clip(dot(n1, n2), -1, 1))`. Degenerate (all-zero / NaN)
   normals are excluded; mean / median / p95 angle in degrees reported.
3. **Per-point RGBA patch bitmaps** (24×24, rendered by
   `refine_normals(render_bitmaps=True)`). Compared per point as the mean
   per-pixel L1 over the RGB channels over the union of covered texels
   (`alpha > 0` on either side); mean / median / p95 of per-point mean-L1
   reported in normalized `[0, 1]` units, plus the count of points whose
   per-point mean-L1 exceeds `16/255 ≈ 0.0627` ("substantially different").

## TL;DR

**Pipeline B (grid) diverges from the SIFT input significantly more than
Pipeline C (LK)** on every dataset. The per-observation keypoint shifts
from SIFT range **0.48–1.09 px mean** for B vs **0.19–0.36 px mean** for
C — roughly a **2.5–3× ratio** in favour of C being closer to SIFT. The
biggest single-dataset divergences from SIFT come from dino_dog_toy: B
moves keypoints 1.09 px mean, C 0.36 px mean. Pipeline A's keypoints are
**exactly** the SIFT positions (0 px shift by definition, since A never
runs a keypoint mover) — its only delta vs SIFT is in normals and bitmaps.

**Normals follow the keypoints, and the second `refine_normals` call
actually does new work.** On every dataset, the normal-angle deltas
order the same way as the keypoint shifts: A < C < B vs SIFT
(seoul_bull: 17.5° / 24.5° / 27.4°; dino: 21.8° / 27.9° / 31.0°;
seattle: 21.1° / 34.0° / 36.5°; kerry: 15.1° / 23.4° / 24.5°). The
pairwise normal angles (B_vs_A, C_vs_A) are 7.9–20.4° mean — *not*
zero — confirming that re-running `refine_normals` with the moved
keypoints lands on a meaningfully different normal than A's normal.

**Bitmaps differ visibly between B and C** on every dataset. The
fraction of points whose bitmap mean-L1 exceeds the
substantially-different threshold (`16/255 ≈ 6.3%`):

| Dataset | B-vs-A (n_diff/n) | C-vs-A (n_diff/n) | C-vs-B (n_diff/n) |
|---|---|---|---|
| seoul_bull | 229 / 836 (27.4%) | 29 / 677 (4.3%) | 214 / 837 (25.6%) |
| dino_dog_toy (stride=5) | 1109 / 3576 (31.0%) | 129 / 3242 (4.0%) | 1088 / 3574 (30.4%) |
| seattle_backyard | 813 / 3226 (25.2%) | 493 / 3142 (15.7%) | 572 / 3225 (17.7%) |
| kerry_park | 265 / 708 (37.4%) | 98 / 551 (17.8%) | 284 / 708 (40.1%) |

The headline pattern: **B and C produce visibly different bitmaps for
roughly a quarter to a third of points across the catalog** (C-vs-B
fractions: 17.7–40.1%). C (LK) stays much closer to A (only ~4–18% of
points cross the threshold), which is consistent with C moving keypoints
far less than B.

**Headline finding for the lead's question** — "how much do the three
pipelines diverge in keypoints/normals/bitmaps relative to the SIFT
input?": **Pipeline B diverges substantially on all three artifacts;
Pipeline C diverges modestly on keypoints but still non-trivially on
normals (because the second `refine_normals` re-fits to the moved
keypoints); Pipeline A leaves keypoints untouched but still rotates
normals 15–22° from the mean-viewing seed.** B and C are *not* small
perturbations of the same answer — the keypoint-mover choice is a real
signal that propagates through the second `refine_normals` and into
the persisted reference appearance.

**One mild surprise.** On kerry_park, the *median* `B_vs_A` and `C_vs_A`
normal angles are essentially 0° (0.000° and 1.309°) despite mean angles
of ~12°. The mean and p95 (~36°) confirm there's real per-point movement
on the tail, but a large fraction of finite points end up with the same
post-refinement normal across A vs B and A vs C. That suggests the
kerry_park (fisheye-rig) view geometry leaves many points with a
well-determined normal that the second refine reproduces regardless of
which keypoint-mover ran upstream. The other three datasets show no
such bimodality (their medians are 4–22°, in line with their means).


## seoul_bull

`20260621-00-solve-seoul_bull_sculpture_1-17.sfmr` — 17 images, 1054 input points, 3137 SIFT-seed observations, 1046 SIFT-seed points.

### Pipeline wall times + outputs

| Pipeline | Label | Wall (s) | Out points | Observations |
|---|---|---:|---:|---:|
| A | normals-only | 0.71 | 1054 | 3137 |
| B | grid (localize_keypoints) | 16.88 | 1043 | 4715 |
| C | LK (refine_keypoints) | 2.87 | 1054 | 3137 |

### Keypoint divergence (per-observation L2 in source-image px)

| Pair | Mean | Median | P95 | n_overlap |
|---|---:|---:|---:|---:|
| A_vs_SIFT | 0.000 | 0.000 | 0.000 | 3137 |
| B_vs_SIFT | 0.850 | 0.465 | 2.926 | 2904 |
| C_vs_SIFT | 0.234 | 0.150 | 0.786 | 3137 |
| B_vs_A | 0.850 | 0.465 | 2.926 | 2904 |
| C_vs_A | 0.234 | 0.150 | 0.786 | 3137 |
| C_vs_B | 0.720 | 0.349 | 2.687 | 2904 |

### Normal divergence (inter-normal angle in degrees)

| Pair | Mean | Median | P95 | n_overlap | n_skipped |
|---|---:|---:|---:|---:|---:|
| A_vs_SIFT | 17.547 | 19.039 | 36.020 | 1046 | 0 |
| B_vs_SIFT | 27.354 | 29.120 | 61.754 | 1035 | 0 |
| C_vs_SIFT | 24.461 | 21.950 | 64.373 | 1046 | 0 |
| B_vs_A | 13.633 | 9.551 | 36.080 | 1035 | 0 |
| C_vs_A | 9.386 | 2.036 | 34.015 | 1046 | 0 |
| C_vs_B | 11.656 | 5.986 | 36.638 | 1035 | 0 |

### Bitmap divergence (per-point mean-L1 over RGB, normalized `[0,1]`)

Threshold for the "substantially different" count: per-point mean-L1 > 0.0627 (= 16/255).

| Pair | Mean | Median | P95 | n_substantially_different / n_compared |
|---|---:|---:|---:|---:|
| A_vs_SIFT | — | — | — | N/A (no SIFT-baseline bitmap (rendering only happens during refine_normals)) |
| B_vs_SIFT | — | — | — | N/A (no SIFT-baseline bitmap (rendering only happens during refine_normals)) |
| C_vs_SIFT | — | — | — | N/A (no SIFT-baseline bitmap (rendering only happens during refine_normals)) |
| B_vs_A | 0.1139 | 0.0307 | 0.4769 | 229 / 836 |
| C_vs_A | 0.0241 | 0.0147 | 0.0584 | 29 / 677 |
| C_vs_B | 0.1081 | 0.0249 | 0.4748 | 214 / 837 |

## dino_dog_toy

`20260621-00-solve-dino_dog_toy_1-85.sfmr` — 85 images, 19024 input points (stride-5 sub-sample: 3805 points), 17581 SIFT-seed observations, 3803 SIFT-seed points.

### Pipeline wall times + outputs

| Pipeline | Label | Wall (s) | Out points | Observations |
|---|---|---:|---:|---:|
| A | normals-only | 11.86 | 3805 | 17581 |
| B | grid (localize_keypoints) | 242.43 | 3802 | 62989 |
| C | LK (refine_keypoints) | 22.20 | 3805 | 17581 |

### Keypoint divergence (per-observation L2 in source-image px)

| Pair | Mean | Median | P95 | n_overlap |
|---|---:|---:|---:|---:|
| A_vs_SIFT | 0.000 | 0.000 | 0.000 | 17581 |
| B_vs_SIFT | 1.086 | 0.840 | 2.527 | 14378 |
| C_vs_SIFT | 0.363 | 0.222 | 1.203 | 17581 |
| B_vs_A | 1.086 | 0.840 | 2.527 | 14378 |
| C_vs_A | 0.363 | 0.222 | 1.203 | 17581 |
| C_vs_B | 0.965 | 0.721 | 2.300 | 14378 |

### Normal divergence (inter-normal angle in degrees)

| Pair | Mean | Median | P95 | n_overlap | n_skipped |
|---|---:|---:|---:|---:|---:|
| A_vs_SIFT | 21.824 | 27.928 | 35.977 | 3803 | 0 |
| B_vs_SIFT | 30.976 | 32.829 | 55.767 | 3800 | 0 |
| C_vs_SIFT | 27.866 | 30.490 | 57.112 | 3803 | 0 |
| B_vs_A | 11.986 | 7.857 | 35.207 | 3800 | 0 |
| C_vs_A | 7.909 | 3.339 | 28.763 | 3803 | 0 |
| C_vs_B | 8.945 | 3.689 | 34.554 | 3800 | 0 |

### Bitmap divergence (per-point mean-L1 over RGB, normalized `[0,1]`)

Threshold for the "substantially different" count: per-point mean-L1 > 0.0627 (= 16/255).

| Pair | Mean | Median | P95 | n_substantially_different / n_compared |
|---|---:|---:|---:|---:|
| A_vs_SIFT | — | — | — | N/A (no SIFT-baseline bitmap (rendering only happens during refine_normals)) |
| B_vs_SIFT | — | — | — | N/A (no SIFT-baseline bitmap (rendering only happens during refine_normals)) |
| C_vs_SIFT | — | — | — | N/A (no SIFT-baseline bitmap (rendering only happens during refine_normals)) |
| B_vs_A | 0.1134 | 0.0431 | 0.5784 | 1109 / 3576 |
| C_vs_A | 0.0279 | 0.0116 | 0.0580 | 129 / 3242 |
| C_vs_B | 0.1134 | 0.0405 | 0.5834 | 1088 / 3574 |

## seattle_backyard

`20260621-00-solve-seattle_backyard_1-26.sfmr` — 26 images, 3343 input points, 14341 SIFT-seed observations, 3281 SIFT-seed points.

### Pipeline wall times + outputs

| Pipeline | Label | Wall (s) | Out points | Observations |
|---|---|---:|---:|---:|
| A | normals-only | 3.99 | 3343 | 14341 |
| B | grid (localize_keypoints) | 109.01 | 3314 | 31040 |
| C | LK (refine_keypoints) | 14.71 | 3343 | 14341 |

### Keypoint divergence (per-observation L2 in source-image px)

| Pair | Mean | Median | P95 | n_overlap |
|---|---:|---:|---:|---:|
| A_vs_SIFT | 0.000 | 0.000 | 0.000 | 14341 |
| B_vs_SIFT | 0.501 | 0.360 | 1.442 | 13280 |
| C_vs_SIFT | 0.228 | 0.141 | 0.770 | 14341 |
| B_vs_A | 0.501 | 0.360 | 1.442 | 13280 |
| C_vs_A | 0.228 | 0.141 | 0.770 | 14341 |
| C_vs_B | 0.378 | 0.253 | 1.149 | 13280 |

### Normal divergence (inter-normal angle in degrees)

| Pair | Mean | Median | P95 | n_overlap | n_skipped |
|---|---:|---:|---:|---:|---:|
| A_vs_SIFT | 21.116 | 28.018 | 36.055 | 3278 | 3 |
| B_vs_SIFT | 36.457 | 37.776 | 66.386 | 3249 | 3 |
| C_vs_SIFT | 34.013 | 35.690 | 66.647 | 3278 | 3 |
| B_vs_A | 20.397 | 22.443 | 36.094 | 3249 | 3 |
| C_vs_A | 16.171 | 15.850 | 35.256 | 3278 | 3 |
| C_vs_B | 13.824 | 8.469 | 39.367 | 3249 | 3 |

### Bitmap divergence (per-point mean-L1 over RGB, normalized `[0,1]`)

Threshold for the "substantially different" count: per-point mean-L1 > 0.0627 (= 16/255).

| Pair | Mean | Median | P95 | n_substantially_different / n_compared |
|---|---:|---:|---:|---:|
| A_vs_SIFT | — | — | — | N/A (no SIFT-baseline bitmap (rendering only happens during refine_normals)) |
| B_vs_SIFT | — | — | — | N/A (no SIFT-baseline bitmap (rendering only happens during refine_normals)) |
| C_vs_SIFT | — | — | — | N/A (no SIFT-baseline bitmap (rendering only happens during refine_normals)) |
| B_vs_A | 0.0603 | 0.0407 | 0.2191 | 813 / 3226 |
| C_vs_A | 0.0384 | 0.0265 | 0.0948 | 493 / 3142 |
| C_vs_B | 0.0511 | 0.0319 | 0.1945 | 572 / 3225 |

## kerry_park

`20260621-00-solve-frame_1-24.sfmr` — 48 images, 786 input points, 2728 SIFT-seed observations, 770 SIFT-seed points.

### Pipeline wall times + outputs

| Pipeline | Label | Wall (s) | Out points | Observations |
|---|---|---:|---:|---:|
| A | normals-only | 0.87 | 786 | 2728 |
| B | grid (localize_keypoints) | 29.64 | 786 | 8319 |
| C | LK (refine_keypoints) | 3.10 | 786 | 2728 |

### Keypoint divergence (per-observation L2 in source-image px)

| Pair | Mean | Median | P95 | n_overlap |
|---|---:|---:|---:|---:|
| A_vs_SIFT | 0.000 | 0.000 | 0.000 | 2728 |
| B_vs_SIFT | 0.483 | 0.342 | 1.379 | 2675 |
| C_vs_SIFT | 0.194 | 0.120 | 0.631 | 2728 |
| B_vs_A | 0.483 | 0.342 | 1.379 | 2675 |
| C_vs_A | 0.194 | 0.120 | 0.631 | 2728 |
| C_vs_B | 0.400 | 0.257 | 1.178 | 2675 |

### Normal divergence (inter-normal angle in degrees)

| Pair | Mean | Median | P95 | n_overlap | n_skipped |
|---|---:|---:|---:|---:|---:|
| A_vs_SIFT | 15.129 | 11.453 | 35.510 | 756 | 14 |
| B_vs_SIFT | 24.470 | 27.038 | 61.833 | 756 | 14 |
| C_vs_SIFT | 23.366 | 19.823 | 64.598 | 756 | 14 |
| B_vs_A | 12.951 | 0.000 | 36.087 | 756 | 14 |
| C_vs_A | 11.454 | 1.309 | 34.108 | 756 | 14 |
| C_vs_B | 15.704 | 13.494 | 39.155 | 756 | 14 |

### Bitmap divergence (per-point mean-L1 over RGB, normalized `[0,1]`)

Threshold for the "substantially different" count: per-point mean-L1 > 0.0627 (= 16/255).

| Pair | Mean | Median | P95 | n_substantially_different / n_compared |
|---|---:|---:|---:|---:|
| A_vs_SIFT | — | — | — | N/A (no SIFT-baseline bitmap (rendering only happens during refine_normals)) |
| B_vs_SIFT | — | — | — | N/A (no SIFT-baseline bitmap (rendering only happens during refine_normals)) |
| C_vs_SIFT | — | — | — | N/A (no SIFT-baseline bitmap (rendering only happens during refine_normals)) |
| B_vs_A | 0.1405 | 0.0472 | 0.5645 | 265 / 708 |
| C_vs_A | 0.0562 | 0.0272 | 0.2743 | 98 / 551 |
| C_vs_B | 0.1386 | 0.0498 | 0.5509 | 284 / 708 |

## Cross-dataset summary

### Mean per-observation keypoint shift vs SIFT, by pipeline (px)

| Dataset | A | B | C |
|---|---:|---:|---:|
| seoul_bull | 0.000 | 0.850 | 0.234 |
| dino_dog_toy (stride=5) | 0.000 | 1.086 | 0.363 |
| seattle_backyard | 0.000 | 0.501 | 0.228 |
| kerry_park | 0.000 | 0.483 | 0.194 |

### Mean normal angle vs SIFT, by pipeline (degrees)

| Dataset | A | B | C |
|---|---:|---:|---:|
| seoul_bull | 17.547 | 27.354 | 24.461 |
| dino_dog_toy (stride=5) | 21.824 | 30.976 | 27.866 |
| seattle_backyard | 21.116 | 36.457 | 34.013 |
| kerry_park | 15.129 | 24.470 | 23.366 |

### Bitmap divergence between B and C (per-point mean-L1 + substantially-different count)

These are the directly informative B-vs-C bitmap comparisons (the *vs_SIFT bitmap rows are N/A — there is no rendered bitmap before the first refine_normals).

| Dataset | mean-L1 | median-L1 | p95-L1 | n_diff / n_compared |
|---|---:|---:|---:|---:|
| seoul_bull | 0.1081 | 0.0249 | 0.4748 | 214 / 837 |
| dino_dog_toy (stride=5) | 0.1134 | 0.0405 | 0.5834 | 1088 / 3574 |
| seattle_backyard | 0.0511 | 0.0319 | 0.1945 | 572 / 3225 |
| kerry_park | 0.1386 | 0.0498 | 0.5509 | 284 / 708 |

## Methodology

### Pipelines

All three pipelines share step 0 (`to_embedded_patches`, `normal='mean_viewing'`, `extent='feature_size'`, `extent_value=5.0`) and step 1 (`refine_normals` with `use_stored_keypoints=True`, `render_bitmaps=True`, `resolution=24` — matching the production `embed_patches` defaults). Pipeline A stops there. Pipelines B and C add a keypoint-moving step and then a second `refine_normals` pass (same knobs) against the moved keypoints.

- Pipeline B's keypoint mover is `PatchCloud.localize_keypoints(search_resolution_multiplier=2.0, max_iters=5, search=6.0, max_shift_px=3.0, min_relative_zncc=0.7)` — the supersampled-grid keypoint search the production `embed_patches` exposes behind `--search-resolution-multiplier 2`. Before localize, `select_views` is run to produce the per-point view sets, matching the production pipeline's step 2.
- Pipeline C's keypoint mover is `PatchCloud.refine_keypoints(max_outer_sweeps=1, consensus_refresh='per_sweep')` — the production `subpixel='lk'` LK refiner (per-sweep consensus, one outer sweep), seeded at the stored SIFT keypoints with the view sets derived from the recon's tracks.

After the keypoint-moving step, both B and C call `compact_to_embedded_patches(min_views=1)` to rebuild an `embedded_patches` recon whose stored `keypoints_xy` reflect the moved keypoints; that recon is what the second `refine_normals` then sees. `min_views=1` keeps every point the mover admits with at least one view, preserving the join-vs-SIFT population.

### Join key

Per-observation keypoint comparisons are joined on `(rounded world position bytes, image_index)` — the same key the existing [`measure_subpixel_decision_gate.py`](../scripts/measure_subpixel_decision_gate.py) uses (the helper `per_obs_keypoints_by_world` is imported verbatim). Per-point normal and bitmap comparisons are joined on the same rounded world position alone (no image_index). Compaction can drop a point between pipelines (B is the one that culls — `localize_keypoints` rejects views and `min_views=1` can leave a point with zero views), so we report the join-overlap size alongside each statistic to keep magnitudes in context.

### Bitmap distance

Per point we have a `24×24×4` `uint8` RGBA texture for each pipeline (the rendered reference appearance the second `refine_normals` pass uses for its consensus). The per-point distance is `mean(|a_rgb - b_rgb|.mean(-1))` over the texels where either pipeline has any alpha (so a texel covered by neither pipeline doesn't artificially deflate the mean). Per-point means are normalized into `[0, 1]` (`/= 255`). The "substantially-different" threshold is `16/255 ≈ 0.0627` — much smaller than typical JPEG-noise levels but well above bilinear-resample jitter, so a point that crosses it really has a different reference appearance, not just a per-texel sub-pixel wobble.

The `*_vs_SIFT` bitmap rows are N/A: there is no rendered bitmap before the first `refine_normals` call, and the SIFT-baseline is `to_embedded_patches` alone (no rendering). Only the three pipeline-to-pipeline bitmap rows (B-vs-A, C-vs-A, C-vs-B) are populated.

### Dino sub-sampling

`dino_dog_toy` is sub-sampled at stride-5 (every 5th point in the reconstruction's existing order) — the same approach `measure_subpixel_decision_gate.py` defaults to. Even at stride-5 it's ~85 high-res images × ~3800 points across three pipelines with two `refine_normals` calls each; the full-point run would push the wall-clock past the budget for a measurement of this kind. The per-pipeline relative comparison the report is built around is invariant to sub-sampling (the per-point and per-observation join keys still match across pipelines on the survivors).

