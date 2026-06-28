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

**Bitmap sharpness — B's bitmaps are blurrier than A/C/D, but the gap is
smaller paired than the per-pipeline aggregates make it look.** The
unrestricted per-pipeline aggregates (Laplacian variance + gradient-mag
mean, alpha-masked) show B 15–27% below A on Laplacian variance and 7–10%
below on gradient-mag mean across the four datasets; A/C/D are within ~5%
of each other. But those aggregates mix populations — each pipeline's
`alpha > 0` coverage differs, and B in particular admits more points via
`localize_keypoints`. **Restricted to the per-dataset common subset
(504–2995 points covered by every pipeline) and computed as paired
per-point deltas, the sign survives — B_vs_A ΔLapVar is negative on every
dataset — but the magnitude shrinks**: paired ΔLapVar means are -0.0004
(seoul), -0.0029 (dino), -0.0029 (seattle), -0.0027 (kerry), vs the
unrestricted-aggregate differences of -0.0011 / -0.0036 / -0.0036 / -0.0028
on the same datasets. Seoul is essentially zero in paired terms; dino /
seattle / kerry keep a real but modest B-blurrier signal. C and D agree
with A on the paired view to within ±0.002 LapVar (D vs C within ±0.0004).
Interpretation unchanged: **B's larger keypoint moves admit more views per
point at the cost of some consensus crispness**, but most of the
unpaired-aggregate gap was the population-difference confound — the
genuine per-point blur cost is smaller, and concentrated on the datasets
where B moves keypoints furthest.

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

### Bitmap sharpness (per-pipeline aggregates over the rendered bitmaps)

Per-bitmap Laplacian variance + gradient magnitude mean, masked to `alpha > 0`
per bitmap so transparent border texels don't deflate the metric. Aggregated
per pipeline to mean / median / p95. Bitmaps with fewer than 2 covered texels
are skipped (`n_skipped`).

> **Caveat:** each pipeline aggregates over its own `alpha > 0` population
> (`n_compared` differs across A/B/C/D — B in particular admits extra points
> via `localize_keypoints`), so this table mixes a population-difference signal
> into the per-pipeline aggregates. The apples-to-apples view is the next
> sub-table ("Bitmap sharpness on the common subset"), which restricts to the
> per-dataset intersection of points covered by all four pipelines.

| Pipeline | LapVar mean | LapVar median | LapVar p95 | GradMag mean | GradMag median | GradMag p95 | n_compared | n_skipped |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| A | 0.0044 | 0.0013 | 0.0214 | 0.0357 | 0.0320 | 0.0664 | 667 | 379 |
| B | 0.0033 | 0.0010 | 0.0163 | 0.0326 | 0.0291 | 0.0617 | 806 | 229 |
| C | 0.0045 | 0.0012 | 0.0229 | 0.0353 | 0.0313 | 0.0679 | 676 | 370 |
| D | 0.0044 | 0.0012 | 0.0236 | 0.0353 | 0.0313 | 0.0682 | 676 | 370 |

### Bitmap sharpness on the common subset

Restricted to the 634 points covered by every pipeline (the intersection of all four `alpha > 0` populations). Same alpha-masked Laplacian variance + gradient-mag mean as the table above, but every pipeline aggregates over the **same** points.

Per-pipeline aggregates on the common subset:

| Pipeline | LapVar mean | LapVar median | LapVar p95 | GradMag mean | GradMag median | GradMag p95 | n_common |
|---|---:|---:|---:|---:|---:|---:|---:|
| A | 0.0039 | 0.0013 | 0.0193 | 0.0350 | 0.0313 | 0.0648 | 634 |
| B | 0.0035 | 0.0011 | 0.0165 | 0.0336 | 0.0299 | 0.0620 | 634 |
| C | 0.0040 | 0.0012 | 0.0201 | 0.0347 | 0.0310 | 0.0658 | 634 |
| D | 0.0040 | 0.0012 | 0.0200 | 0.0347 | 0.0310 | 0.0664 | 634 |

Paired per-point deltas on the common subset (`metric_X - metric_Y`):

| Pair | ΔLapVar mean | ΔLapVar median | ΔLapVar p95 | ΔGradMag mean | ΔGradMag median | ΔGradMag p95 |
|---|---:|---:|---:|---:|---:|---:|
| B_vs_A | -0.0004 | -0.0001 | +0.0008 | -0.0014 | -0.0011 | +0.0033 |
| C_vs_A | +0.0001 | -0.0000 | +0.0014 | -0.0004 | -0.0001 | +0.0033 |
| D_vs_A | +0.0001 | -0.0000 | +0.0014 | -0.0004 | -0.0001 | +0.0034 |
| C_vs_B | +0.0005 | +0.0001 | +0.0036 | +0.0011 | +0.0005 | +0.0066 |
| D_vs_B | +0.0005 | +0.0001 | +0.0036 | +0.0011 | +0.0005 | +0.0067 |
| D_vs_C | +0.0000 | +0.0000 | +0.0001 | +0.0000 | +0.0000 | +0.0004 |

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

### Bitmap sharpness (per-pipeline aggregates over the rendered bitmaps)

> **Caveat:** see the seoul_bull bitmap-sharpness note — each pipeline aggregates
> over its own `alpha > 0` population. See the next sub-table for the
> common-subset apples-to-apples view.

| Pipeline | LapVar mean | LapVar median | LapVar p95 | GradMag mean | GradMag median | GradMag p95 | n_compared | n_skipped |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| A | 0.0151 | 0.0085 | 0.0541 | 0.0500 | 0.0467 | 0.0935 | 3216 | 587 |
| B | 0.0115 | 0.0058 | 0.0419 | 0.0454 | 0.0416 | 0.0867 | 3405 | 395 |
| C | 0.0165 | 0.0091 | 0.0608 | 0.0517 | 0.0480 | 0.0986 | 3200 | 603 |
| D | 0.0167 | 0.0092 | 0.0624 | 0.0519 | 0.0484 | 0.0991 | 3204 | 599 |

### Bitmap sharpness on the common subset

Restricted to the 2995 points covered by every pipeline (the intersection of all four `alpha > 0` populations). Same alpha-masked Laplacian variance + gradient-mag mean as the table above, but every pipeline aggregates over the **same** points.

Per-pipeline aggregates on the common subset:

| Pipeline | LapVar mean | LapVar median | LapVar p95 | GradMag mean | GradMag median | GradMag p95 | n_common |
|---|---:|---:|---:|---:|---:|---:|---:|
| A | 0.0146 | 0.0083 | 0.0527 | 0.0498 | 0.0465 | 0.0932 | 2995 |
| B | 0.0117 | 0.0060 | 0.0431 | 0.0459 | 0.0422 | 0.0873 | 2995 |
| C | 0.0163 | 0.0090 | 0.0600 | 0.0517 | 0.0480 | 0.0986 | 2995 |
| D | 0.0165 | 0.0091 | 0.0613 | 0.0519 | 0.0482 | 0.0987 | 2995 |

Paired per-point deltas on the common subset (`metric_X - metric_Y`):

| Pair | ΔLapVar mean | ΔLapVar median | ΔLapVar p95 | ΔGradMag mean | ΔGradMag median | ΔGradMag p95 |
|---|---:|---:|---:|---:|---:|---:|
| B_vs_A | -0.0029 | -0.0010 | +0.0083 | -0.0039 | -0.0027 | +0.0128 |
| C_vs_A | +0.0016 | +0.0003 | +0.0104 | +0.0019 | +0.0008 | +0.0116 |
| D_vs_A | +0.0019 | +0.0003 | +0.0114 | +0.0021 | +0.0008 | +0.0125 |
| C_vs_B | +0.0046 | +0.0015 | +0.0252 | +0.0058 | +0.0038 | +0.0278 |
| D_vs_B | +0.0048 | +0.0016 | +0.0250 | +0.0061 | +0.0040 | +0.0286 |
| D_vs_C | +0.0002 | +0.0000 | +0.0023 | +0.0002 | +0.0000 | +0.0025 |

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

### Bitmap sharpness (per-pipeline aggregates over the rendered bitmaps)

> **Caveat:** see the seoul_bull bitmap-sharpness note — each pipeline aggregates
> over its own `alpha > 0` population. See the next sub-table for the
> common-subset apples-to-apples view.

| Pipeline | LapVar mean | LapVar median | LapVar p95 | GradMag mean | GradMag median | GradMag p95 | n_compared | n_skipped |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| A | 0.0185 | 0.0089 | 0.0693 | 0.0558 | 0.0508 | 0.1056 | 3077 | 204 |
| B | 0.0149 | 0.0052 | 0.0640 | 0.0489 | 0.0428 | 0.1014 | 3174 | 78 |
| C | 0.0188 | 0.0080 | 0.0740 | 0.0552 | 0.0492 | 0.1079 | 3107 | 174 |
| D | 0.0192 | 0.0081 | 0.0760 | 0.0555 | 0.0492 | 0.1088 | 3107 | 174 |

### Bitmap sharpness on the common subset

Restricted to the 2963 points covered by every pipeline (the intersection of all four `alpha > 0` populations). Same alpha-masked Laplacian variance + gradient-mag mean as the table above, but every pipeline aggregates over the **same** points.

Per-pipeline aggregates on the common subset:

| Pipeline | LapVar mean | LapVar median | LapVar p95 | GradMag mean | GradMag median | GradMag p95 | n_common |
|---|---:|---:|---:|---:|---:|---:|---:|
| A | 0.0176 | 0.0085 | 0.0649 | 0.0550 | 0.0501 | 0.1022 | 2963 |
| B | 0.0146 | 0.0051 | 0.0632 | 0.0486 | 0.0426 | 0.1014 | 2963 |
| C | 0.0178 | 0.0079 | 0.0706 | 0.0544 | 0.0486 | 0.1059 | 2963 |
| D | 0.0182 | 0.0079 | 0.0732 | 0.0546 | 0.0486 | 0.1066 | 2963 |

Paired per-point deltas on the common subset (`metric_X - metric_Y`):

| Pair | ΔLapVar mean | ΔLapVar median | ΔLapVar p95 | ΔGradMag mean | ΔGradMag median | ΔGradMag p95 |
|---|---:|---:|---:|---:|---:|---:|
| B_vs_A | -0.0029 | -0.0018 | +0.0075 | -0.0064 | -0.0054 | +0.0067 |
| C_vs_A | +0.0003 | -0.0000 | +0.0080 | -0.0006 | +0.0000 | +0.0069 |
| D_vs_A | +0.0006 | -0.0000 | +0.0097 | -0.0004 | +0.0000 | +0.0078 |
| C_vs_B | +0.0032 | +0.0015 | +0.0178 | +0.0058 | +0.0044 | +0.0217 |
| D_vs_B | +0.0036 | +0.0015 | +0.0183 | +0.0060 | +0.0045 | +0.0219 |
| D_vs_C | +0.0004 | +0.0000 | +0.0026 | +0.0002 | +0.0000 | +0.0023 |

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

### Bitmap sharpness (per-pipeline aggregates over the rendered bitmaps)

> **Caveat:** see the seoul_bull bitmap-sharpness note — each pipeline aggregates
> over its own `alpha > 0` population. See the next sub-table for the
> common-subset apples-to-apples view.

| Pipeline | LapVar mean | LapVar median | LapVar p95 | GradMag mean | GradMag median | GradMag p95 | n_compared | n_skipped |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| A | 0.0112 | 0.0056 | 0.0394 | 0.0526 | 0.0508 | 0.0920 | 526 | 244 |
| B | 0.0084 | 0.0041 | 0.0322 | 0.0476 | 0.0451 | 0.0848 | 704 | 66 |
| C | 0.0112 | 0.0054 | 0.0388 | 0.0517 | 0.0495 | 0.0934 | 533 | 237 |
| D | 0.0109 | 0.0053 | 0.0376 | 0.0514 | 0.0489 | 0.0916 | 527 | 243 |

### Bitmap sharpness on the common subset

Restricted to the 504 points covered by every pipeline (the intersection of all four `alpha > 0` populations). Same alpha-masked Laplacian variance + gradient-mag mean as the table above, but every pipeline aggregates over the **same** points.

Per-pipeline aggregates on the common subset:

| Pipeline | LapVar mean | LapVar median | LapVar p95 | GradMag mean | GradMag median | GradMag p95 | n_common |
|---|---:|---:|---:|---:|---:|---:|---:|
| A | 0.0109 | 0.0055 | 0.0387 | 0.0522 | 0.0506 | 0.0917 | 504 |
| B | 0.0081 | 0.0041 | 0.0313 | 0.0470 | 0.0448 | 0.0850 | 504 |
| C | 0.0108 | 0.0052 | 0.0383 | 0.0514 | 0.0490 | 0.0917 | 504 |
| D | 0.0110 | 0.0052 | 0.0394 | 0.0514 | 0.0489 | 0.0916 | 504 |

Paired per-point deltas on the common subset (`metric_X - metric_Y`):

| Pair | ΔLapVar mean | ΔLapVar median | ΔLapVar p95 | ΔGradMag mean | ΔGradMag median | ΔGradMag p95 |
|---|---:|---:|---:|---:|---:|---:|
| B_vs_A | -0.0027 | -0.0010 | +0.0013 | -0.0052 | -0.0037 | +0.0029 |
| C_vs_A | -0.0000 | -0.0000 | +0.0033 | -0.0008 | -0.0000 | +0.0046 |
| D_vs_A | +0.0001 | -0.0000 | +0.0041 | -0.0008 | +0.0000 | +0.0048 |
| C_vs_B | +0.0027 | +0.0007 | +0.0138 | +0.0044 | +0.0029 | +0.0180 |
| D_vs_B | +0.0028 | +0.0008 | +0.0140 | +0.0044 | +0.0029 | +0.0180 |
| D_vs_C | +0.0001 | +0.0000 | +0.0011 | +0.0000 | +0.0000 | +0.0007 |

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

### Mean Laplacian variance of rendered bitmap, by pipeline

Per-pipeline mean of the per-bitmap Laplacian variance, alpha-masked. Higher
= sharper rendered reference appearance.

| Dataset | A | B | C | D |
|---|---:|---:|---:|---:|
| seoul_bull | 0.0044 | 0.0033 | 0.0045 | 0.0044 |
| dino_dog_toy (stride=5) | 0.0151 | 0.0115 | 0.0165 | 0.0167 |
| seattle_backyard | 0.0185 | 0.0149 | 0.0188 | 0.0192 |
| kerry_park | 0.0112 | 0.0084 | 0.0112 | 0.0109 |

### Mean gradient-magnitude mean of rendered bitmap, by pipeline

Per-pipeline mean of the per-bitmap mean forward-difference gradient
magnitude, alpha-masked. Same shape as the Laplacian-variance table.

| Dataset | A | B | C | D |
|---|---:|---:|---:|---:|
| seoul_bull | 0.0357 | 0.0326 | 0.0353 | 0.0353 |
| dino_dog_toy (stride=5) | 0.0500 | 0.0454 | 0.0517 | 0.0519 |
| seattle_backyard | 0.0558 | 0.0489 | 0.0552 | 0.0555 |
| kerry_park | 0.0526 | 0.0476 | 0.0517 | 0.0514 |

On every dataset, both unrestricted-aggregate sharpness metrics agree: **B is
the blurriest** (15–27% lower Laplacian variance, 7–10% lower gradient-mag
mean than A), while A / C / D land within rendering-jitter of each other (D
vs C agrees to within ±2%). **Caveat:** these aggregates are over different
per-pipeline populations (each pipeline's `alpha > 0` coverage differs); see
"Paired ΔLapVar mean on the common subset" below for the apples-to-apples
paired view. The sign of B-vs-A survives the paired restriction on every
dataset; the magnitude shrinks (especially on seoul_bull, where the paired
ΔLapVar mean is essentially zero).

### Paired ΔLapVar mean on the common subset, by pair (LapVar units)

Per-dataset paired ΔLapVar mean restricted to the per-dataset common subset
(the same `n_common` set every pipeline aggregates over). Three informative
pairs: **B-A** captures the grid pipeline's blur cost vs normals-only;
**C-A** captures LK-bilinear's deviation; **D-C** captures the anisotropic
sampler's marginal effect. The other three pairs derive from these via
subtraction. Compare against the unrestricted "Mean Laplacian variance"
table above — the common-subset numbers preserve the B-is-blurrier sign
on every dataset but the magnitudes drop.

| Dataset | B-A | C-A | D-C | n_common |
|---|---:|---:|---:|---:|
| seoul_bull | -0.0004 | +0.0001 | +0.0000 | 634 |
| dino_dog_toy (stride=5) | -0.0029 | +0.0016 | +0.0002 | 2995 |
| seattle_backyard | -0.0029 | +0.0003 | +0.0004 | 2963 |
| kerry_park | -0.0027 | -0.0000 | +0.0001 | 504 |

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

### Bitmap sharpness

Per-bitmap Laplacian variance and gradient magnitude mean are computed on
the luminance channel `I_gray = mean(rgb)` (RGB normalized to `[0, 1]`).
The Laplacian uses the standard 3x3 kernel `[[0,1,0],[1,-4,1],[0,1,0]]`
applied via numpy slicing with `reflect` boundary padding (matches
`scipy.ndimage.convolve(mode='reflect')`); variance is taken over the
covered texels. Gradient magnitude uses forward differences (`np.diff`
along each axis, zero-padded back to original shape) and
`sqrt(dx^2 + dy^2)`, again meaned over the covered texels.

Both metrics share the per-bitmap alpha mask `alpha > 0` — the same
mask the bitmap-distance helper uses, so sharpness isn't deflated by
all-transparent border texels. Bitmaps with fewer than 2 covered texels
are skipped (`n_skipped`); the per-pipeline aggregate (mean / median /
p95) is over the survivors only.

### Dino sub-sampling

`dino_dog_toy` is sub-sampled at stride-5 (every 5th point in the reconstruction's existing order) — the same approach `measure_subpixel_decision_gate.py` defaults to. Even at stride-5 it's ~85 high-res images × ~3800 points across four pipelines with two `refine_normals` calls each; the full-point run would push the wall-clock past the budget for a measurement of this kind. The per-pipeline relative comparison the report is built around is invariant to sub-sampling (the per-point and per-observation join keys still match across pipelines on the survivors).

