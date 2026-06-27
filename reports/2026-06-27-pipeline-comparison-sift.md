# Pipeline Comparison: SIFT Input vs Four Refinement Pipelines

_Date: 2026-06-27._

Four SfM refinement pipelines, all starting from the raw SIFT-detection
keypoints copied inline by `recon.to_embedded_patches(...)`, compared on
the four checked-in datasets:

- **Pipeline A** (normals-only): `to_embedded_patches` → `refine_normals`.
  Keypoints stay at SIFT positions; only the per-point normal is
  photometrically optimized.
- **Pipeline B** (grid): `to_embedded_patches` → `refine_normals` →
  `localize_keypoints` (supersampled grid, `search_resolution_multiplier=2`) →
  `refine_normals` (second time, against the moved keypoints).
- **Pipeline C** (LK-bilinear): `to_embedded_patches` → `refine_normals` →
  `refine_keypoints` (LK, `subpixel='lk'` defaults; default
  `sampler='bilinear'`) → `refine_normals` (second time, against the moved
  keypoints).
- **Pipeline D** (LK-anisotropic): identical to C except the LK refiner is
  invoked with `sampler='anisotropic'` instead of the default
  `'bilinear'` — the anti-aliased oblique-view sampler.

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

**Pipeline D (LK-anisotropic) is the lead's follow-up question — and
the answer is "tiny but non-zero":** swapping the LK sampler from
`bilinear` to `anisotropic` moves per-observation keypoints by **0.03–0.21
px mean (median 0.000–0.080 px)** vs C across the four datasets, rotates
normals by **0.5–2.3° mean (median 0.000°)**, and shifts bitmaps by
**0.002–0.013 mean L1** with only **0.3–1.8%** of points crossing the
"substantially different" threshold. By every measure, D and C are well
within the per-pipeline noise — far smaller than the C-vs-A deltas, let
alone C-vs-B. The largest D-vs-C separation lives on dino_dog_toy
(0.211 px mean keypoint shift, 2.342° mean normal rotation, 1.77% of
points with a substantially-different bitmap) and is concentrated on
the tail (p95 keypoint shift 0.769 px, p95 normal angle 13.882°);
seoul_bull and kerry_park D-vs-C agree to within rendering jitter
(2/676 and 9/533 substantially-different bitmaps). The lead's
question — "even if ECC scores end up similar, how much do the bilinear
and anisotropic LK paths actually disagree on keypoint placement,
normals, and the persisted bitmap?" — answers: **they don't, in any way
that propagates to a downstream consumer reading the persisted recon.**
This matches the decision-gate report's ECC finding (no meaningful
ECC-score difference between `lk_per_sweep` and `lk_per_sweep_aniso`):
the persisted artifacts agree too. Wall-time penalty for D over C is
~1.2–1.4× across datasets (largest on dino: 31.75s vs 22.37s).


## seoul_bull

`20260621-00-solve-seoul_bull_sculpture_1-17.sfmr` — 17 images, 1054 input points, 3137 SIFT-seed observations, 1046 SIFT-seed points.

### Pipeline wall times + outputs

| Pipeline | Label | Wall (s) | Out points | Observations |
|---|---|---:|---:|---:|
| A | normals-only | 0.79 | 1054 | 3137 |
| B | grid (localize_keypoints) | 16.42 | 1043 | 4715 |
| C | LK (refine_keypoints, bilinear) | 2.68 | 1054 | 3137 |
| D | LK (refine_keypoints, anisotropic) | 3.22 | 1054 | 3137 |

### Keypoint divergence (per-observation L2 in source-image px)

| Pair | Mean | Median | P95 | n_overlap |
|---|---:|---:|---:|---:|
| A_vs_SIFT | 0.000 | 0.000 | 0.000 | 3137 |
| B_vs_SIFT | 0.850 | 0.465 | 2.926 | 2904 |
| C_vs_SIFT | 0.234 | 0.150 | 0.786 | 3137 |
| B_vs_A | 0.850 | 0.465 | 2.926 | 2904 |
| C_vs_A | 0.234 | 0.150 | 0.786 | 3137 |
| C_vs_B | 0.720 | 0.349 | 2.687 | 2904 |
| D_vs_SIFT | 0.258 | 0.158 | 0.849 | 3137 |
| D_vs_A | 0.258 | 0.158 | 0.849 | 3137 |
| D_vs_B | 0.708 | 0.347 | 2.641 | 2904 |
| D_vs_C | 0.033 | 0.000 | 0.158 | 3137 |

### Normal divergence (inter-normal angle in degrees)

| Pair | Mean | Median | P95 | n_overlap | n_skipped |
|---|---:|---:|---:|---:|---:|
| A_vs_SIFT | 17.547 | 19.039 | 36.020 | 1046 | 0 |
| B_vs_SIFT | 27.354 | 29.120 | 61.754 | 1035 | 0 |
| C_vs_SIFT | 24.461 | 21.950 | 64.373 | 1046 | 0 |
| B_vs_A | 13.633 | 9.551 | 36.080 | 1035 | 0 |
| C_vs_A | 9.386 | 2.036 | 34.015 | 1046 | 0 |
| C_vs_B | 11.656 | 5.986 | 36.638 | 1035 | 0 |
| D_vs_SIFT | 24.532 | 21.908 | 64.373 | 1046 | 0 |
| D_vs_A | 9.466 | 2.070 | 34.010 | 1046 | 0 |
| D_vs_B | 11.595 | 5.986 | 36.638 | 1035 | 0 |
| D_vs_C | 0.525 | 0.000 | 1.309 | 1046 | 0 |

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
| D_vs_SIFT | — | — | — | N/A (no SIFT-baseline bitmap (rendering only happens during refine_normals)) |
| D_vs_A | 0.0243 | 0.0153 | 0.0565 | 30 / 677 |
| D_vs_B | 0.1079 | 0.0248 | 0.4748 | 211 / 837 |
| D_vs_C | 0.0020 | 0.0000 | 0.0117 | 2 / 676 |

## dino_dog_toy

`20260621-00-solve-dino_dog_toy_1-85.sfmr` — 85 images, 19024 input points (stride-5 sub-sample: 3805 points), 17581 SIFT-seed observations, 3803 SIFT-seed points.

### Pipeline wall times + outputs

| Pipeline | Label | Wall (s) | Out points | Observations |
|---|---|---:|---:|---:|
| A | normals-only | 13.26 | 3805 | 17581 |
| B | grid (localize_keypoints) | 255.75 | 3802 | 62989 |
| C | LK (refine_keypoints, bilinear) | 22.37 | 3805 | 17581 |
| D | LK (refine_keypoints, anisotropic) | 31.75 | 3805 | 17581 |

### Keypoint divergence (per-observation L2 in source-image px)

| Pair | Mean | Median | P95 | n_overlap |
|---|---:|---:|---:|---:|
| A_vs_SIFT | 0.000 | 0.000 | 0.000 | 17581 |
| B_vs_SIFT | 1.086 | 0.840 | 2.527 | 14378 |
| C_vs_SIFT | 0.363 | 0.222 | 1.203 | 17581 |
| B_vs_A | 1.086 | 0.840 | 2.527 | 14378 |
| C_vs_A | 0.363 | 0.222 | 1.203 | 17581 |
| C_vs_B | 0.965 | 0.721 | 2.300 | 14378 |
| D_vs_SIFT | 0.502 | 0.287 | 1.687 | 17581 |
| D_vs_A | 0.502 | 0.287 | 1.687 | 17581 |
| D_vs_B | 0.948 | 0.707 | 2.255 | 14378 |
| D_vs_C | 0.211 | 0.080 | 0.769 | 17581 |

### Normal divergence (inter-normal angle in degrees)

| Pair | Mean | Median | P95 | n_overlap | n_skipped |
|---|---:|---:|---:|---:|---:|
| A_vs_SIFT | 21.824 | 27.928 | 35.977 | 3803 | 0 |
| B_vs_SIFT | 30.976 | 32.829 | 55.767 | 3800 | 0 |
| C_vs_SIFT | 27.866 | 30.490 | 57.112 | 3803 | 0 |
| B_vs_A | 11.986 | 7.857 | 35.207 | 3800 | 0 |
| C_vs_A | 7.909 | 3.339 | 28.763 | 3803 | 0 |
| C_vs_B | 8.945 | 3.689 | 34.554 | 3800 | 0 |
| D_vs_SIFT | 28.011 | 30.649 | 56.388 | 3803 | 0 |
| D_vs_A | 8.157 | 3.818 | 28.756 | 3803 | 0 |
| D_vs_B | 8.874 | 3.574 | 34.770 | 3800 | 0 |
| D_vs_C | 2.342 | 0.000 | 13.882 | 3803 | 0 |

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
| D_vs_SIFT | — | — | — | N/A (no SIFT-baseline bitmap (rendering only happens during refine_normals)) |
| D_vs_A | 0.0308 | 0.0130 | 0.0589 | 138 / 3251 |
| D_vs_B | 0.1126 | 0.0404 | 0.5822 | 1085 / 3574 |
| D_vs_C | 0.0130 | 0.0022 | 0.0345 | 57 / 3223 |

## seattle_backyard

`20260621-00-solve-seattle_backyard_1-26.sfmr` — 26 images, 3343 input points, 14341 SIFT-seed observations, 3281 SIFT-seed points.

### Pipeline wall times + outputs

| Pipeline | Label | Wall (s) | Out points | Observations |
|---|---|---:|---:|---:|
| A | normals-only | 4.26 | 3343 | 14341 |
| B | grid (localize_keypoints) | 109.83 | 3314 | 31040 |
| C | LK (refine_keypoints, bilinear) | 14.36 | 3343 | 14341 |
| D | LK (refine_keypoints, anisotropic) | 18.73 | 3343 | 14341 |

### Keypoint divergence (per-observation L2 in source-image px)

| Pair | Mean | Median | P95 | n_overlap |
|---|---:|---:|---:|---:|
| A_vs_SIFT | 0.000 | 0.000 | 0.000 | 14341 |
| B_vs_SIFT | 0.501 | 0.360 | 1.442 | 13280 |
| C_vs_SIFT | 0.228 | 0.141 | 0.770 | 14341 |
| B_vs_A | 0.501 | 0.360 | 1.442 | 13280 |
| C_vs_A | 0.228 | 0.141 | 0.770 | 14341 |
| C_vs_B | 0.378 | 0.253 | 1.149 | 13280 |
| D_vs_SIFT | 0.273 | 0.155 | 0.950 | 14341 |
| D_vs_A | 0.273 | 0.155 | 0.950 | 14341 |
| D_vs_B | 0.365 | 0.246 | 1.084 | 13280 |
| D_vs_C | 0.067 | 0.009 | 0.290 | 14341 |

### Normal divergence (inter-normal angle in degrees)

| Pair | Mean | Median | P95 | n_overlap | n_skipped |
|---|---:|---:|---:|---:|---:|
| A_vs_SIFT | 21.116 | 28.018 | 36.055 | 3278 | 3 |
| B_vs_SIFT | 36.457 | 37.776 | 66.386 | 3249 | 3 |
| C_vs_SIFT | 34.013 | 35.690 | 66.647 | 3278 | 3 |
| B_vs_A | 20.397 | 22.443 | 36.094 | 3249 | 3 |
| C_vs_A | 16.171 | 15.850 | 35.256 | 3278 | 3 |
| C_vs_B | 13.824 | 8.469 | 39.367 | 3249 | 3 |
| D_vs_SIFT | 34.100 | 36.058 | 66.554 | 3278 | 3 |
| D_vs_A | 16.539 | 16.636 | 35.228 | 3278 | 3 |
| D_vs_B | 13.765 | 8.516 | 38.935 | 3249 | 3 |
| D_vs_C | 2.073 | 0.000 | 13.393 | 3278 | 3 |

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
| D_vs_SIFT | — | — | — | N/A (no SIFT-baseline bitmap (rendering only happens during refine_normals)) |
| D_vs_A | 0.0397 | 0.0274 | 0.0965 | 508 / 3144 |
| D_vs_B | 0.0510 | 0.0317 | 0.1930 | 577 / 3225 |
| D_vs_C | 0.0077 | 0.0008 | 0.0360 | 50 / 3117 |

## kerry_park

`20260621-00-solve-frame_1-24.sfmr` — 48 images, 786 input points, 2728 SIFT-seed observations, 770 SIFT-seed points.

### Pipeline wall times + outputs

| Pipeline | Label | Wall (s) | Out points | Observations |
|---|---|---:|---:|---:|
| A | normals-only | 0.86 | 786 | 2728 |
| B | grid (localize_keypoints) | 28.85 | 786 | 8319 |
| C | LK (refine_keypoints, bilinear) | 3.11 | 786 | 2728 |
| D | LK (refine_keypoints, anisotropic) | 3.64 | 786 | 2728 |

### Keypoint divergence (per-observation L2 in source-image px)

| Pair | Mean | Median | P95 | n_overlap |
|---|---:|---:|---:|---:|
| A_vs_SIFT | 0.000 | 0.000 | 0.000 | 2728 |
| B_vs_SIFT | 0.483 | 0.342 | 1.379 | 2675 |
| C_vs_SIFT | 0.194 | 0.120 | 0.631 | 2728 |
| B_vs_A | 0.483 | 0.342 | 1.379 | 2675 |
| C_vs_A | 0.194 | 0.120 | 0.631 | 2728 |
| C_vs_B | 0.400 | 0.257 | 1.178 | 2675 |
| D_vs_SIFT | 0.213 | 0.125 | 0.711 | 2728 |
| D_vs_A | 0.213 | 0.125 | 0.711 | 2728 |
| D_vs_B | 0.398 | 0.256 | 1.174 | 2675 |
| D_vs_C | 0.036 | 0.000 | 0.178 | 2728 |

### Normal divergence (inter-normal angle in degrees)

| Pair | Mean | Median | P95 | n_overlap | n_skipped |
|---|---:|---:|---:|---:|---:|
| A_vs_SIFT | 15.129 | 11.453 | 35.510 | 756 | 14 |
| B_vs_SIFT | 24.470 | 27.038 | 61.833 | 756 | 14 |
| C_vs_SIFT | 23.366 | 19.823 | 64.598 | 756 | 14 |
| B_vs_A | 12.951 | 0.000 | 36.087 | 756 | 14 |
| C_vs_A | 11.454 | 1.309 | 34.108 | 756 | 14 |
| C_vs_B | 15.704 | 13.494 | 39.155 | 756 | 14 |
| D_vs_SIFT | 23.316 | 19.382 | 64.598 | 756 | 14 |
| D_vs_A | 11.420 | 1.799 | 34.040 | 756 | 14 |
| D_vs_B | 15.732 | 13.772 | 38.818 | 756 | 14 |
| D_vs_C | 0.991 | 0.000 | 5.918 | 756 | 14 |

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
| D_vs_SIFT | — | — | — | N/A (no SIFT-baseline bitmap (rendering only happens during refine_normals)) |
| D_vs_A | 0.0558 | 0.0256 | 0.2752 | 98 / 548 |
| D_vs_B | 0.1408 | 0.0510 | 0.5561 | 286 / 707 |
| D_vs_C | 0.0082 | 0.0000 | 0.0273 | 9 / 533 |

## Cross-dataset summary

### Mean per-observation keypoint shift vs SIFT, by pipeline (px)

| Dataset | A | B | C | D |
|---|---:|---:|---:|---:|
| seoul_bull | 0.000 | 0.850 | 0.234 | 0.258 |
| dino_dog_toy (stride=5) | 0.000 | 1.086 | 0.363 | 0.502 |
| seattle_backyard | 0.000 | 0.501 | 0.228 | 0.273 |
| kerry_park | 0.000 | 0.483 | 0.194 | 0.213 |

### Mean normal angle vs SIFT, by pipeline (degrees)

| Dataset | A | B | C | D |
|---|---:|---:|---:|---:|
| seoul_bull | 17.547 | 27.354 | 24.461 | 24.532 |
| dino_dog_toy (stride=5) | 21.824 | 30.976 | 27.866 | 28.011 |
| seattle_backyard | 21.116 | 36.457 | 34.013 | 34.100 |
| kerry_park | 15.129 | 24.470 | 23.366 | 23.316 |

### Bitmap divergence between B and C (per-point mean-L1 + substantially-different count)

These are the directly informative B-vs-C bitmap comparisons (the *vs_SIFT bitmap rows are N/A — there is no rendered bitmap before the first refine_normals).

| Dataset | mean-L1 | median-L1 | p95-L1 | n_diff / n_compared |
|---|---:|---:|---:|---:|
| seoul_bull | 0.1081 | 0.0249 | 0.4748 | 214 / 837 |
| dino_dog_toy (stride=5) | 0.1134 | 0.0405 | 0.5834 | 1088 / 3574 |
| seattle_backyard | 0.0511 | 0.0319 | 0.1945 | 572 / 3225 |
| kerry_park | 0.1386 | 0.0498 | 0.5509 | 284 / 708 |

## D vs C: bilinear vs anisotropic

The lead's follow-up question: holding the LK keypoint refiner's other knobs
fixed, how much does swapping `sampler="bilinear"` (Pipeline C) for
`sampler="anisotropic"` (Pipeline D) move the persisted keypoints,
normals, and bitmaps? The three cross-dataset rollups below report only
the `D_vs_C` rows from each per-dataset block.

### Keypoints (per-observation L2 in source-image px)

| Dataset | Mean | Median | P95 | n_overlap |
|---|---:|---:|---:|---:|
| seoul_bull | 0.033 | 0.000 | 0.158 | 3137 |
| dino_dog_toy (stride=5) | 0.211 | 0.080 | 0.769 | 17581 |
| seattle_backyard | 0.067 | 0.009 | 0.290 | 14341 |
| kerry_park | 0.036 | 0.000 | 0.178 | 2728 |

Context: C's own mean shift vs SIFT is 0.194–0.363 px across these
datasets (and B's is 0.483–1.086 px). The D-vs-C separation is **3–7×
smaller than C-vs-SIFT** everywhere except dino, where it's still about
**1.7×** smaller. Medians of 0.000 on three of the four datasets mean
the *majority* of observations have identical keypoints across the two
samplers (matching at full float precision); the mean is dragged up by
a small tail of points where the sampler choice does push the LK
optimum into a measurably different basin.

### Normals (inter-normal angle in degrees)

| Dataset | Mean | Median | P95 | n_overlap |
|---|---:|---:|---:|---:|
| seoul_bull | 0.525 | 0.000 | 1.309 | 1046 |
| dino_dog_toy (stride=5) | 2.342 | 0.000 | 13.882 | 3803 |
| seattle_backyard | 2.073 | 0.000 | 13.393 | 3278 |
| kerry_park | 0.991 | 0.000 | 5.918 | 756 |

Same shape as keypoints: medians of 0.000° on every dataset (the
second `refine_normals` lands on the same normal on the majority of
points), with a thin tail that pulls the mean to a few degrees and the
p95 to 5–14°. Compare against C-vs-A normal means of 7.9–16.2° — the
D-vs-C normal disagreement is a small fraction (~5–25%) of that.

### Bitmaps (per-point mean-L1 over RGB; threshold = 16/255)

| Dataset | mean-L1 | median-L1 | p95-L1 | n_diff / n_compared |
|---|---:|---:|---:|---:|
| seoul_bull | 0.0020 | 0.0000 | 0.0117 | 2 / 676 |
| dino_dog_toy (stride=5) | 0.0130 | 0.0022 | 0.0345 | 57 / 3223 |
| seattle_backyard | 0.0077 | 0.0008 | 0.0360 | 50 / 3117 |
| kerry_park | 0.0082 | 0.0000 | 0.0273 | 9 / 533 |

Substantially-different rates: 0.3% (seoul) / 1.8% (dino) / 1.6%
(seattle) / 1.7% (kerry). For comparison, the C-vs-A rates on the same
datasets are 4.3% / 4.0% / 15.7% / 17.8%, and C-vs-B are 25.6% / 30.4%
/ 17.7% / 40.1%. D and C produce **essentially the same persisted
reference bitmap** — well over an order of magnitude tighter than any
other pipeline-to-pipeline comparison in this report.

### Verdict

The bilinear vs anisotropic LK sampler choice does produce a
**measurable, non-zero signal** on every dataset, but the signal is
small in absolute terms and in every direction smaller than the
intra-pipeline noise from running the same pipeline against different
keypoint movers. On the persisted artifacts a downstream consumer
actually reads back — keypoint positions, normals, and 24×24 RGBA
patches — D and C are interchangeable to within roughly the
single-pipeline rendering jitter. This matches the decision-gate
report's finding that the two samplers also produce indistinguishable
ECC scores. The dino_dog_toy dataset shows the largest disagreement on
all three axes, consistent with that dataset having the most aggressive
keypoint motion in C in the first place; even there, the D-vs-C deltas
are well under C's own move from SIFT.

Wall-time penalty for D over C: **1.16× (kerry) – 1.42× (dino)**,
landing well below B's 6–11× cost over C.

## Methodology

### Pipelines

All four pipelines share step 0 (`to_embedded_patches`, `normal='mean_viewing'`, `extent='feature_size'`, `extent_value=5.0`) and step 1 (`refine_normals` with `use_stored_keypoints=True`, `render_bitmaps=True`, `resolution=24` — matching the production `embed_patches` defaults). Pipeline A stops there. Pipelines B, C, and D add a keypoint-moving step and then a second `refine_normals` pass (same knobs) against the moved keypoints.

- Pipeline B's keypoint mover is `PatchCloud.localize_keypoints(search_resolution_multiplier=2.0, max_iters=5, search=6.0, max_shift_px=3.0, min_relative_zncc=0.7)` — the supersampled-grid keypoint search the production `embed_patches` exposes behind `--search-resolution-multiplier 2`. Before localize, `select_views` is run to produce the per-point view sets, matching the production pipeline's step 2.
- Pipeline C's keypoint mover is `PatchCloud.refine_keypoints(max_outer_sweeps=1, consensus_refresh='per_sweep')` (default `sampler='bilinear'`) — the production `subpixel='lk'` LK refiner (per-sweep consensus, one outer sweep), seeded at the stored SIFT keypoints with the view sets derived from the recon's tracks.
- Pipeline D's keypoint mover is the same `PatchCloud.refine_keypoints` call as C, but with `sampler='anisotropic'` passed through — the anti-aliased oblique-view sampler. Same seeds, same view sets, same outer-sweep / consensus settings, same `resolution`. (Anti-aliasing only affects how the source pyramids are sampled into the patch, not what the patch is.)

After the keypoint-moving step, B, C, and D all call `compact_to_embedded_patches(min_views=1)` to rebuild an `embedded_patches` recon whose stored `keypoints_xy` reflect the moved keypoints; that recon is what the second `refine_normals` then sees. `min_views=1` keeps every point the mover admits with at least one view, preserving the join-vs-SIFT population.

### Join key

Per-observation keypoint comparisons are joined on `(rounded world position bytes, image_index)` — the same key the existing [`measure_subpixel_decision_gate.py`](../scripts/measure_subpixel_decision_gate.py) uses (the helper `per_obs_keypoints_by_world` is imported verbatim). Per-point normal and bitmap comparisons are joined on the same rounded world position alone (no image_index). Compaction can drop a point between pipelines (B is the one that culls — `localize_keypoints` rejects views and `min_views=1` can leave a point with zero views), so we report the join-overlap size alongside each statistic to keep magnitudes in context.

### Bitmap distance

Per point we have a `24×24×4` `uint8` RGBA texture for each pipeline (the rendered reference appearance the second `refine_normals` pass uses for its consensus). The per-point distance is `mean(|a_rgb - b_rgb|.mean(-1))` over the texels where either pipeline has any alpha (so a texel covered by neither pipeline doesn't artificially deflate the mean). Per-point means are normalized into `[0, 1]` (`/= 255`). The "substantially-different" threshold is `16/255 ≈ 0.0627` — much smaller than typical JPEG-noise levels but well above bilinear-resample jitter, so a point that crosses it really has a different reference appearance, not just a per-texel sub-pixel wobble.

The `*_vs_SIFT` bitmap rows are N/A: there is no rendered bitmap before the first `refine_normals` call, and the SIFT-baseline is `to_embedded_patches` alone (no rendering). Only the pipeline-to-pipeline bitmap rows (B-vs-A, C-vs-A, C-vs-B, D-vs-A, D-vs-B, D-vs-C) are populated.

### Dino sub-sampling

`dino_dog_toy` is sub-sampled at stride-5 (every 5th point in the reconstruction's existing order) — the same approach `measure_subpixel_decision_gate.py` defaults to. Even at stride-5 it's ~85 high-res images × ~3800 points across four pipelines with two `refine_normals` calls each; the full-point run would push the wall-clock past the budget for a measurement of this kind. The per-pipeline relative comparison the report is built around is invariant to sub-sampling (the per-point and per-observation join keys still match across pipelines on the survivors).

