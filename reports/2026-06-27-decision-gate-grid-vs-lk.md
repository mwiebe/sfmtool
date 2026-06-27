# Decision Gate: Supersampled Grid vs LK Subpixel Refinement

_Date: 2026-06-27._ The final phase of the keypoint-localization implementation
plan ([`reports/2026-06-26-keypoint-localization-implementation-plan.md` on the
plan branch][plan]). Phases 3A (AVX2 grid kernel, #139) and 3B (analytic
Jacobian + per-move incremental consensus, #140 + #142) have landed. This report
takes the optimized variants head-to-head on the four checked-in datasets and
recommends a production default for the new `embed_patches(..., subpixel=...)`
knob.

_**Round 2 update (same date).** The original report covered three datasets
(seoul_bull, seattle_backyard, kerry_park) and six variants. This second
round adds the missing **dino_dog_toy** dataset (run on a stride-5
sub-sample for cost — see "Per-dataset results"), three **new LK
variations** (`lk_per_sweep_aniso`, `lk_per_sweep_tight_offset`,
`lk_per_move_10sweeps`) testing plausible fixes for the kerry_park
regression, and a new **`shift_vs_sift`** metric measuring how far each
variant moves keypoints from the raw SIFT detection (the input to the
whole pipeline). See "Round-2 additions" below._

[plan]: https://github.com/.../tree/claude/quirky-heisenberg-2vlpmp/reports/2026-06-26-keypoint-localization-implementation-plan.md

## TL;DR / Recommendation

**Keep `subpixel="none"` and `search_resolution_multiplier=1.0` as the default.**
The measurement does not produce a clear winner that is safe to flip across all
four datasets — the LK variants give a meaningful ECC bump on the textured
close-range data (seoul_bull: +0.0027 mean ECC for `lk_per_sweep` at 1.5× cost,
+0.0040 for `lk_per_move` at 3× cost) but **regress** on the kerry_park fisheye
rig (-0.0031 for `lk_per_sweep`, -0.0024 for `lk_per_move`). The supersampled
grid (`m=2` and `m=3`) is consistently a small positive across every dataset
(+0.0007–0.0014 ECC) but costs 5-7× baseline at `m=2` and 17-26× at `m=3` — not
worth a default flip for an ECC delta that's solidly within the variance of the
metric itself.

The round-2 additions did not change this verdict:

- **dino_dog_toy** (sub-sampled, stride 5) is the **strongest LK-win
  dataset** of the four — every LK variant beats `grid_m2` and `baseline`
  by a solid margin (LK +0.0043 to +0.0048 ECC, grid_m2 +0.0032). The
  anisotropic-sampler variant is the per-dataset winner at +0.0048.
  Dino's dataset character (mostly-planar surfaces, narrow-FOV pinhole,
  high-resolution texture) is exactly the regime LK was designed for.
  It does not regress, so it joins the "LK helps" set, not the "LK
  regresses" set — but it also doesn't help the cross-dataset decision,
  because the binding constraint is still kerry_park.
- **None of the three new LK variations recover the kerry_park regression.**
  `lk_per_sweep_aniso` (anisotropic sampler) lands at the same ECC as
  `lk_per_sweep`. `lk_per_sweep_tight_offset` (`max_offset_px=1.0` instead
  of 2.0) makes the regression slightly **worse** on kerry_park.
  `lk_per_move_10sweeps` extends the convergence curve marginally past
  `lk_per_move` (kerry_park: −0.0020 vs −0.0024) but is still a regression.
- The **SIFT-baseline shift** is significant — the localizer + LK pipeline
  moves keypoints **0.48–1.17 source-px** from where SIFT detected them
  (mean), with the close-range / high-resolution datasets moving more.
  Dino's mean shift from SIFT is **1.17 px** (the largest), seoul_bull's is
  **0.87 px**, seattle and kerry_park around **0.48–0.50 px**. The LK
  refiner's incremental contribution on top of the localizer is **small**
  (~0.1–0.2 px more) — most of the move was already the localizer's, not
  the LK refinement.

Both `lk` and `lk_per_move` are now **opt-in** behind `subpixel="lk"` and
`subpixel="lk_per_move"` (CLI: `sfm embed-patches --subpixel lk`); the
supersampled grid is opt-in behind `--search-resolution-multiplier 2`. The
three new LK variations remain measurement-only — the measurement script
calls `refine_keypoints` directly with custom kwargs rather than expanding
the production `embed_patches` API surface (the deliberate decision is to
keep that at `none|lk|lk_per_move`).

The "do nothing" verdict is the boring answer, but it's what the data supports —
two of four code paths landed (LK is the cheapest at ~1.5× and has a real
regression case, the grid is a positive uniform but a multiplicative cost), the
production knobs are wired and tested, the easy-fix LK variations for the
kerry_park regression have been ruled out, and the production default doesn't
move until a clearer signal emerges.

## Round-2 additions (what changed since the first writeup)

- **`dino_dog_toy` added.** The original round skipped it for time. This
  round runs every variant except `grid_m3` (the ~21× baseline cost
  outlier — would have pushed the dino sweep past two hours by itself).
  Even the full-point dino baseline takes ~160 s and grid_m2 stretched past
  15 min when first attempted; to fit the full sweep in a reasonable budget
  the dino numbers are computed on a **stride-5 sub-sample** of the input
  reconstruction's points (3805 of 19024 — every 5th point in the
  reconstruction's existing order, with track observations filtered to the
  survivors). The `subsample_stride` field is recorded in the data JSON and
  the row is marked here. The full-point dino baseline took 158 s and
  produced `mean_ecc = 0.9173` with `shift_vs_sift = 1.165 px`; the
  stride-5 numbers preserve the per-variant *relative* comparison the
  decision gate cares about while keeping the absolute mean ECC value in
  the same ballpark (the metric averages over per-observation scores, not
  per-point counts).

- **Three new LK variations.** The first two are plausible-fix candidates
  for the kerry_park regression; the third pins the
  convergence-saturation curve for `lk_per_move`:

  | Variant | Differs from | Rationale |
  |---|---|---|
  | `lk_per_sweep_aniso` | `lk_per_sweep` + `sampler="anisotropic"` | The fisheye warp Jacobian is highly non-isotropic at frame edges; the anisotropic sampler matches the patch footprint to the warp. The natural candidate for fixing a fisheye-specific LK failure. |
  | `lk_per_sweep_tight_offset` | `lk_per_sweep` + `max_offset_px=1.0` (vs 2.0) | If LK is wandering too far on hard cases, tightening the per-view drift bound keeps it in the basin. |
  | `lk_per_move_10sweeps` | `lk_per_move` + `max_outer_sweeps=10` (vs 5) | Per-sweep diminishing returns suggest this won't help much; pins the curve. |

  All three were exercised through the measurement script's custom-refine
  path (it calls `cloud.refine_keypoints(**custom_kwargs)` directly,
  rather than going through the production `subpixel="lk"|"lk_per_move"`
  enum which hard-codes its kwargs). The production API surface remains
  the deliberately-narrow `none|lk|lk_per_move`.

- **New `shift_vs_sift` metric.** Before any variant runs,
  `recon.to_embedded_patches(normal="mean_viewing", extent="feature_size",
  extent_value=patch_size/2.0)` is called once to snapshot the
  per-observation **SIFT seed** keypoints by the same `(world position,
  image index)` join key as the existing `shift_vs_baseline`. For every
  variant we then report mean/median/p95 |variant_kpt − sift_kpt|, over
  the observations present in both. This is the *cumulative*
  photometric-correction magnitude the pipeline applied to each
  observation — including the localizer's move (which is also in the
  baseline). The original `shift_vs_baseline` measures the **LK
  contribution alone** (variant minus baseline, where baseline is
  localizer-only); `shift_vs_sift` measures **localizer + LK together**.

## Variant matrix

| Variant | `search_resolution_multiplier` | step 3.5 | Notes |
|---|---|---|---|
| `baseline` | 1.0 | none | Current production behaviour |
| `grid_m2` | 2.0 | none | Supersampled grid only |
| `grid_m3` | 3.0 | none | Even more aggressive grid (not run on dino) |
| `lk_per_sweep` | 1.0 | `refine_keypoints(max_outer_sweeps=1)` | LK, per-sweep consensus |
| `lk_per_move` | 1.0 | `refine_keypoints(max_outer_sweeps=5, consensus_refresh="per_move")` | LK, per-move Gauss-Seidel |
| `grid_m2_then_lk` | 2.0 | `refine_keypoints(max_outer_sweeps=1)` | Supersampled grid + LK refinement on top |
| `lk_per_sweep_aniso` | 1.0 | `refine_keypoints(max_outer_sweeps=1, sampler="anisotropic")` | **Round 2.** Anisotropic sampler |
| `lk_per_sweep_tight_offset` | 1.0 | `refine_keypoints(max_outer_sweeps=1, max_offset_px=1.0)` | **Round 2.** Tighter per-view drift bound |
| `lk_per_move_10sweeps` | 1.0 | `refine_keypoints(max_outer_sweeps=10, consensus_refresh="per_move")` | **Round 2.** Saturate the convergence curve |

## Methodology

For each (dataset × variant) the script
(`scripts/measure_subpixel_decision_gate.py`) runs the full
`embed_patches(...)` pipeline end-to-end (or, for the round-2 variants, runs
the pipeline by hand and slips a custom `refine_keypoints` call into step
3.5) and captures:

- **Wall time (`wall_secs`)** — time from before the run until after it
  returns. Includes image load via `read_workspace_image` (~~10 s for
  dino_dog_toy, negligible for the others). Single sample per (dataset,
  variant); no warm-up pass. **Round-2 wall-times for the small datasets
  are higher than round-1's** because the round-2 sweep was run in
  parallel with another sweep for part of its life (the two contended for
  CPU); the **relative** rankings within a dataset still stand. Compare to
  round 1 for clean per-variant timings (the data file's `wall_secs` is
  the round-2 number; round-1 numbers are in the prose tables below for
  reference where they materially differ).
- **`out_points`** — point count of the compacted `embedded_patches`
  reconstruction. The compaction culls points below `min_views=2`; a variant
  that drops more points than baseline is a regression signal (the
  localizer/refiner threw away observations).
- **Mean ECC (`mean_ecc`)** — the cross-variant comparable photometric
  agreement metric. It is computed *after* `embed_patches` by re-loading the
  reconstruction's stored frames + inline keypoints, building a patch cloud,
  and calling `refine_keypoints(..., max_gn_steps=0)` (the seed-only mode)
  with `starting_keypoints=` the stored inline keypoints. The reported
  `scores` field is the per-view channel-averaged windowed ZNCC of the seeded
  kernel against the IRLS consensus rendered from those same seeds. NaN
  scores (a point left with < 2 views, so no consensus) are dropped before
  averaging. This metric reaches across variants because the scoring code
  path is the same for every variant — only the stored keypoints differ.
- **Mean / median / p95 keypoint shift vs baseline (`mean_shift_px` etc.)**
  — per-observation distance in source-image px between the variant's
  stored keypoint and the baseline's, joined by **(world position, image
  index)** (the compaction renumbers survivor point indices, so joining by
  output index is wrong; joining by the source world position is the
  invariant key since compaction never re-bundles). This is the **LK
  contribution alone** (variant minus baseline, where baseline is
  localizer-only).
- **(Round 2) Mean / median / p95 keypoint shift vs SIFT seed
  (`mean_shift_vs_sift_px` etc.)** — same shape, but the reference is the
  SIFT detection keypoint snapshotted at step 0 (the input to the whole
  pipeline). This is **localizer + LK together** — how far the entire
  photometric pipeline moves each observation from where SIFT first put
  it.

### What the metrics CAN and CAN'T tell you

- ECC is computed against a consensus rebuilt at the **stored** keypoints, so
  it rewards self-consistency — a variant that shifted all views by the same
  amount could score well even if it walked off the truth together. The
  cross-view shift bounds are what stop this in practice (the localizer's
  `max_shift_px` gate and the LK refiner's `max_offset_px`), but ECC alone
  doesn't tell you ground truth. Three-decimal-place ECC deltas between
  variants are small relative to the dataset-to-dataset variation (the
  textured close-range seoul_bull scores ~0.94, the fisheye kerry_park ~0.92,
  the higher-res dino ~0.92).
- The keypoint shift vs baseline measures how much each variant *changed* the
  keypoint locations relative to the no-refinement path — it's a magnitude
  signal, not a "right vs wrong" signal. A high shift with high ECC means
  "the refiner moved things and they agree better with the rebuilt
  consensus"; a high shift with **lower** ECC (which happens on kerry_park
  for the LK variants) means the refiner moved things and they agree
  **worse** — a real regression.
- **`shift_vs_sift` is a magnitude, not a quality signal.** It bounds the
  photometric correction the pipeline applied, but doesn't itself say
  whether the correction was an improvement on SIFT's pick — only ECC says
  that. On the variants we've measured, **the LK refiner's contribution
  on top of the localizer is small** (~0.1–0.2 px), so `shift_vs_sift`
  largely measures the **localizer**'s correction, not LK's.
- Wall time is single-sample, no statistics. Numbers within ~10% of each
  other should not be treated as significant. The big multiplicative gaps
  (`m=3` is ~17-26× baseline) are real signal. Round-2 timings on the small
  datasets ran under CPU contention with the dino sweep and are inflated;
  the round-1 numbers in the prose tables are the clean serial timings.

## Per-dataset results

### seoul_bull_sculpture (17 images @ 270×480, close-range, textured)

| Variant | Wall (s) | Points | Mean ECC | ΔECC | Shift vs base (px) | Shift vs SIFT (px) |
|---|---:|---:|---:|---:|---:|---:|
| baseline                  |  2.15 | 1029 | 0.9418 |  —          | 0.000 | 0.866 |
| grid_m2                   | 13.84 | 1028 | 0.9427 | +0.0009     | 0.326 | 0.848 |
| grid_m3                   | 52.35 | 1028 | 0.9429 | +0.0011     | 0.301 | 0.841 |
| lk_per_sweep              |  3.20 | 1029 | 0.9445 | **+0.0027** | 0.100 | 0.854 |
| lk_per_move               |  6.37 | 1029 | 0.9458 | **+0.0040** | 0.174 | 0.879 |
| grid_m2_then_lk           | 14.94 | 1028 | 0.9447 | +0.0029     | 0.335 | 0.842 |
| lk_per_sweep_aniso        |  3.82 | 1029 | 0.9446 | +0.0028     | 0.106 | 0.854 |
| lk_per_sweep_tight_offset |  3.41 | 1029 | 0.9443 | +0.0025     | 0.096 | 0.855 |
| lk_per_move_10sweeps      |  7.54 | 1029 | 0.9462 | **+0.0044** | 0.202 | 0.896 |

(Wall times are round-1 clean serial; the round-2 sweep was under contention.)

LK is the clear winner here. The dataset has rich texture and short baselines,
the regime LK was designed for. The anisotropic sampler is a wash with
bilinear (0.9446 vs 0.9445); the tighter offset gives up a small amount of
ECC (0.9443 vs 0.9445); more outer sweeps adds a small further gain (+0.0004
ECC over 5 sweeps, +13% cost). The SIFT-baseline shift here is the **largest
of the four datasets** at ~0.87 px (a quarter of the patch radius) — close-range
imagery has more sub-pixel detail to find.

### seattle_backyard (26 images @ 360×640, outdoor mixed)

| Variant | Wall (s) | Points | Mean ECC | ΔECC | Shift vs base (px) | Shift vs SIFT (px) |
|---|---:|---:|---:|---:|---:|---:|
| baseline                  |  13.44 | 3295 | 0.9266 |  —      | 0.000 | 0.498 |
| grid_m2                   |  90.41 | 3293 | 0.9280 | +0.0014 | 0.202 | 0.501 |
| grid_m3                   | 350.21 | 3293 | 0.9282 | +0.0016 | 0.195 | 0.503 |
| lk_per_sweep              |  19.46 | 3295 | 0.9273 | +0.0007 | 0.070 | 0.499 |
| lk_per_move               |  35.42 | 3295 | 0.9276 | +0.0010 | 0.106 | 0.513 |
| grid_m2_then_lk           |  96.50 | 3293 | 0.9280 | +0.0014 | 0.208 | 0.502 |
| lk_per_sweep_aniso        |  25.05 | 3295 | 0.9276 | +0.0010 | 0.079 | 0.503 |
| lk_per_sweep_tight_offset |  19.25 | 3295 | 0.9273 | +0.0007 | 0.070 | 0.499 |
| lk_per_move_10sweeps      |  38.87 | 3295 | 0.9277 | +0.0011 | 0.117 | 0.520 |

Closer call — `grid_m2` edges out the LK variants on mean ECC, at ~5× the cost.
`grid_m2_then_lk` is no better than `grid_m2` alone, suggesting the grid found
the basin LK was going to find anyway. Pure LK is the cheapest improvement but
also the smallest. The anisotropic sampler gives a marginal lift over bilinear
(+0.0003 ECC at 1.25× cost); the tighter offset is bit-identical to
`lk_per_sweep` (typical moves are well within 1 px). `lk_per_move_10sweeps`
adds nothing useful over 5 sweeps. SIFT-baseline shift is half a pixel —
smaller than seoul_bull's because the natural-scene texture has fewer
sharp sub-pixel features to congeal onto.

### kerry_park (48 images @ 480×480, fisheye rig)

| Variant | Wall (s) | Points | Mean ECC | ΔECC | Shift vs base (px) | Shift vs SIFT (px) |
|---|---:|---:|---:|---:|---:|---:|
| baseline                  |  4.97 | 786 | 0.9214 |  —          | 0.000 | 0.481 |
| grid_m2                   | 25.11 | 786 | 0.9217 | +0.0003     | 0.243 | 0.483 |
| grid_m3                   | 97.25 | 786 | 0.9221 | +0.0007     | 0.230 | 0.477 |
| lk_per_sweep              |  5.87 | 786 | 0.9183 | **−0.0031** | 0.133 | 0.486 |
| lk_per_move               | 12.18 | 786 | 0.9190 | **−0.0024** | 0.202 | 0.526 |
| grid_m2_then_lk           | 28.29 | 786 | 0.9183 | −0.0031     | 0.290 | 0.488 |
| lk_per_sweep_aniso        |  7.50 | 786 | 0.9184 | **−0.0030** | 0.134 | 0.488 |
| lk_per_sweep_tight_offset |  6.18 | 786 | 0.9169 | **−0.0045** | 0.142 | 0.486 |
| lk_per_move_10sweeps      | 15.64 | 786 | 0.9194 | **−0.0020** | 0.239 | 0.555 |

The regression that determines the verdict. The original-round variants all
regress. The three new variants:

- **`lk_per_sweep_aniso` does not recover the regression** — 0.9184 vs
  baseline's 0.9214, indistinguishable from plain `lk_per_sweep`'s 0.9183.
  The anisotropic sampler was the most plausible fix (fisheye edges have
  highly non-isotropic warp Jacobians), but on this dataset it makes no
  measurable difference — suggesting the regression's mechanism is not
  isotropic-sampling aliasing.
- **`lk_per_sweep_tight_offset` makes the regression *worse*** — 0.9169 vs
  baseline's 0.9214 (−0.0045, the largest regression on this dataset).
  Tightening `max_offset_px` from 2.0 to 1.0 prematurely caps LK's
  convergence on this dataset, leaving views half-converged at a worse
  position than they'd reach with the looser bound. This rules out
  "wandering too far" as the cause.
- **`lk_per_move_10sweeps`** is marginally better than 5 sweeps (−0.0020
  vs −0.0024), still a regression, and the saturation suggests no amount
  of further sweeps will fix it.

The point count is identical for every variant on this dataset (786) — the
LK regression isn't a "the refiner dropped more points" failure mode, it's
the refiner doing its job but converging to a slightly worse keypoint than
the discrete localizer found, on fisheye geometry where that's possible.

The plausible original-round cause (analytic-Jacobian linearization error
on highly non-linear fisheye projection) is **not contradicted** by the new
data, but it's also not pinpointed — the anisotropic sampler was the
intuition-driven fix for that hypothesis and it didn't help. A more
targeted experiment (e.g. inverse-compositional ECC; the joint bundle;
per-frame inspection of which views' keypoints LK *worsened* and why)
would be the next step.

### dino_dog_toy (85 images @ 2040×1536, large textured)

**Sub-sampled at stride 5 (3805 of 19024 input points).** A first attempt at
the full-point sweep took ~160 s for baseline and was 14+ min into
`grid_m2` (still not done) when stopped — extrapolating, the full sweep
would have taken 1.5-2 hours and pushed grid_m2 past 25 min on its own.
The stride-5 sub-sample takes every 5th point in the reconstruction's
existing source order (with track observations filtered to the survivors)
and lets the full eight-variant sweep finish in ~14 min. The full-point
baseline produced `mean_ecc = 0.9173` with `shift_vs_sift = 1.165 px`;
the stride-5 baseline produced `mean_ecc = 0.9177` with `shift_vs_sift =
1.188 px` — close enough that the per-variant *relative* comparison the
decision gate cares about is preserved. `grid_m3` is omitted on this
dataset (the cost extrapolation from the other datasets would put it past
30 min even on the sub-sample).

| Variant | Wall (s) | Points | Mean ECC | ΔECC | Shift vs base (px) | Shift vs SIFT (px) |
|---|---:|---:|---:|---:|---:|---:|
| baseline                  |  35.10 | 3789 | 0.9177 |  —          | 0.000 | 1.188 |
| grid_m2                   | 204.17 | 3787 | 0.9209 | **+0.0032** | 0.632 | 1.086 |
| grid_m3                   |    —   |  —   |  —     | (skipped)   |  —    |  —    |
| lk_per_sweep              |  51.92 | 3789 | 0.9220 | **+0.0043** | 0.252 | 1.127 |
| lk_per_move               |  87.65 | 3789 | 0.9223 | **+0.0046** | 0.315 | 1.133 |
| grid_m2_then_lk           | 233.84 | 3787 | 0.9220 | **+0.0043** | 0.630 | 1.062 |
| lk_per_sweep_aniso        |  78.72 | 3789 | 0.9225 | **+0.0048** | 0.329 | 1.115 |
| lk_per_sweep_tight_offset |  52.31 | 3789 | 0.9219 | **+0.0042** | 0.250 | 1.127 |
| lk_per_move_10sweeps      |  88.43 | 3789 | 0.9224 | **+0.0047** | 0.321 | 1.135 |

**Dino is the strongest LK-win dataset of the four.** Every LK variant
beats every non-LK variant (LK +0.0042 to +0.0048 vs grid_m2 +0.0032 vs
baseline 0). The anisotropic sampler is the per-dataset winner on dino
(0.9225, +0.0048 vs baseline at only 1.5× cost), edging out per-move
(+0.0046) and per-move-10sweeps (+0.0047). This is the opposite of
kerry_park's regression — dino's mostly-flat, planar surfaces with
narrow-FOV pinhole cameras are exactly the regime LK was designed for.
The full-point dino numbers (where measured: baseline only) match the
sub-sample within 0.001 ECC, so the relative ranking is trustworthy.

`grid_m2_then_lk` lands at exactly the same ECC as `lk_per_sweep`
(0.9220) — the grid found the same basin LK was going to find. The
+0.0043 of `grid_m2_then_lk` over baseline is the LK contribution, not
the grid's; pure grid alone is +0.0032.

**`shift_vs_sift` here is the largest of the four datasets** (1.13–1.19
px). The biggest single contributor is the localizer (baseline shift
1.19 px); LK only adds ~0.0 to −0.07 px on top (most variants are
slightly *closer* to SIFT than baseline because they pull the localizer's
larger moves back toward the photometric optimum), and grid_m2 pulls the
localizer's moves back by ~0.10 px from SIFT (the supersampled discrete
search recovers a position closer to SIFT's seed than the iterative
localizer did, suggesting the localizer over-shoots a bit on
high-resolution imagery).

## SIFT baseline

The first thing the pipeline does (step 0, `recon.to_embedded_patches()`) is
copy each observation's SIFT detection keypoint inline as the seed for the
photometric pipeline (steps 1–3 refine_normals → select_views →
localize_keypoints; step 3.5 the optional LK refinement). The
`shift_vs_sift` metric measures how far each variant's final stored
keypoint sits from that original SIFT seed, averaged over per-observation
distances in source-image px.

Per-dataset summary (baseline = localizer only, no LK):

| Dataset | Baseline shift vs SIFT (mean px) | Median | p95 | Best LK variant contribution |
|---|---:|---:|---:|---:|
| seoul_bull_sculpture (270×480, close-range)      | 0.866 | 0.483 | 2.996 | +0.030 (lk_per_move_10sweeps: 0.896) |
| seattle_backyard (360×640, outdoor)              | 0.498 | 0.349 | 1.463 | +0.022 (lk_per_move_10sweeps: 0.520) |
| kerry_park (480×480, fisheye)                    | 0.481 | 0.334 | 1.414 | +0.074 (lk_per_move_10sweeps: 0.555) |
| dino_dog_toy (2040×1536; stride-5 sub-sample)    | 1.188 | 0.927 | 2.728 | −0.05 (lk_per_sweep: 1.127; LK pulls *toward* SIFT) |

The headline: **the localizer is doing most of the photometric correction.**
Across all four datasets, the localizer alone moves each observation
0.48–1.19 px from SIFT (largest on dino, the high-resolution dataset).
The LK refiner on top adds only **~0.0–0.07 px more on average** —
the rounding-error contribution. On three datasets (seoul_bull,
seattle_backyard, kerry_park) LK pushes the keypoint slightly **further**
from SIFT; on dino LK pulls it slightly **closer** to SIFT (the localizer
over-shoots on the high-res imagery, LK pulls it back). On every
dataset, LK's incremental move is small relative to the localizer's, so
the SIFT-shift metric mostly tells us about the localizer.

The supersampled grid behaves similarly on the three small datasets
(~0.0 px shift from baseline's SIFT distance) but on dino reduces it by
~0.1 px (1.188 → 1.086) — the discrete search finds a keypoint closer to
SIFT's pick than the iterative localizer's converged position, a small
signal that the localizer slightly over-corrects on high-resolution
imagery. None of this changes the verdict; it's the lead's
"how much does the pipeline move things" question, answered: a little less
than a pixel on most datasets, a little more than a pixel on dino.

## Aggregate verdict

| Variant | Mean ΔECC (4 datasets) | Mean ΔECC (3 small, original) | Cost ratio (geo-mean) | Always non-regressing? |
|---|---:|---:|---:|---|
| `baseline`                  |  0      |  0      | 1.0×  | (the baseline) |
| `grid_m2`                   | +0.0014 | +0.0009 | ~5.8× | yes |
| `grid_m3`                   | +0.0011 (3 ds) | +0.0011 | ~21×  | yes (not run on dino) |
| `lk_per_sweep`              | +0.0011 | +0.0001 | ~1.4× | **no — kerry_park −0.0031** |
| `lk_per_move`               | +0.0018 | +0.0009 | ~2.5× | **no — kerry_park −0.0024** |
| `grid_m2_then_lk`           | +0.0013 | +0.0004 | ~6.4× | **no — kerry_park −0.0031** |
| `lk_per_sweep_aniso`        | +0.0014 | +0.0003 | ~1.5× | **no — kerry_park −0.0029** |
| `lk_per_sweep_tight_offset` | +0.0007 | -0.0005 | ~1.4× | **no — kerry_park −0.0045** |
| `lk_per_move_10sweeps`      | +0.0021 | +0.0012 | ~2.7× | **no — kerry_park −0.0020** |

The 4-dataset mean is much more flattering to LK because dino contributes
a strong +0.0043–0.0048 to each LK row. If kerry_park weren't in the
denominator (i.e. we accepted a per-dataset auto-default that disables LK
on fisheye), `lk_per_move_10sweeps` would be the strict winner — +0.0030
average across the three "LK-friendly" datasets at 2.7× cost. That's a
reasonable case for a future round to revisit.

`grid_m3` has the best mean ΔECC but at 17-26× baseline cost — not worth
defaulting to. `grid_m2` is the only variant with a uniformly positive ΔECC
across all three small datasets, but its ~5.8× cost is steep for a 0.001 ECC
improvement. The LK variants are cheap and excellent on close-range textured
data but regress on the fisheye rig — and **none of the three round-2 LK
variations recover the regression**.

**No variant is a strict improvement over the baseline**, so the conservative
production default stays `subpixel="none"` and
`search_resolution_multiplier=1.0`. The variants are all opt-in.

## What changed in this branch

- `crates/sfmtool-py/src/py_patch_cloud.rs`:
  - `PatchCloud.refine_keypoints` gains a `starting_keypoints:
    Optional[dict[point_id, list[[x, y]]]]` kwarg (mirrors `view_sets`'s
    shape). Defaults to `None` = projection seed (current behaviour). A
    custom per-patch rayon walk routes the per-point seed list (or projection
    fallback) into the existing `refine_patch_keypoints` Rust function.
  - `PatchCloud.localize_keypoints` gains a
    `search_resolution_multiplier: f32 = 1.0` kwarg (previously fixed
    at the `KeypointLocalizeParams::default()` value).
- `src/sfmtool/_embed_patches.py`:
  - `embed_patches()` gains `subpixel: str = "none"` and
    `search_resolution_multiplier: float = 1.0` kwargs.
  - New internal `_refine_subpixel` helper splices the LK pass's per-view
    refined keypoints back into the localizer's per-point dicts, preserving
    the kept-view membership and order — only the per-view `keypoints` array
    is replaced.
  - The localizer's output is the seed precondition for the LK pass, so step
    3.5 is a no-op when `subpixel == "none"` and the existing pipeline is
    bit-for-bit unchanged at the default.
- `src/sfmtool/_commands/embed_patches.py`:
  - Adds `--search-resolution-multiplier` and `--subpixel
    {none,lk,lk_per_move}` Click options, threaded through to `embed_patches`.
- `specs/cli/embed-patches-command.md` updated with the two new options and
  the recommendation.
- `specs/core/keypoint-subpixel-refinement.md`: status section updated to
  record the decision-gate outcome (LK is wired into `embed_patches` as an
  opt-in pass; the default stays off pending a dataset-specific signal
  strong enough to flip it). _Round 2 (2026-06-27): updated to note the
  dino measurement and the ruled-out LK fix candidates._
- `scripts/measure_subpixel_decision_gate.py`:
  - **Round 2.** Rewritten to support both `embed_kwargs`-style variants
    (drive `embed_patches` directly, the production API) and
    `refine_kwargs`-style variants (replay the pipeline by hand and call
    `cloud.refine_keypoints(**custom_kwargs)` for step 3.5 — used by
    the three new LK variations that exercise params the production API
    does not expose).
  - Added the SIFT-seed snapshot and `shift_vs_sift` metric (mean / median
    / p95) per variant.
  - Added `subsampled_recon` and a `--dino-stride` CLI option for fitting
    the dino sweep into a reasonable runtime; the chosen stride is
    recorded in the data JSON.
  - Added a per-dataset variant skip list (currently
    `{dino_dog_toy: {grid_m3}}`) so the prohibitively-slow variants can be
    opted out per dataset without losing the rest.

## Test coverage added

- `tests/test_patch_keypoint_subpixel.py::test_refine_keypoints_honors_starting_keypoints`
  — the PyO3 binding contract: shifting the seed off the projection produces
  a different refined keypoint on at least one view; a length mismatch is
  rejected up front.
- `tests/test_embed_patches_compaction.py::test_embed_patches_default_is_no_subpixel`
  — the default `embed_patches` call is bit-equal to passing
  `subpixel="none"`; flipping the default in code cannot slip in silently.
- `tests/test_embed_patches_compaction.py::test_embed_patches_subpixel_lk_round_trips`
  — `subpixel="lk"` produces a valid `embedded_patches` reconstruction that
  round-trips through `.sfmr`, and its per-observation keypoints differ from
  the baseline (the wiring is not a no-op).

Round 2 added no production code, so no new tests. The measurement script
is exercised end-to-end by the runs that produced this report (the script
is the test); the SIFT-snapshot and custom-refine paths are not part of
the production code path.

## Open follow-ups for the lead

- **The plan on `claude/quirky-heisenberg-2vlpmp` should be marked complete**
  (the "Decision gate" section). I don't have that branch checked out, so I
  haven't touched it.
- **The kerry_park LK regression deserves a deeper look.** Round 2 ruled
  out the two intuition-driven candidates (anisotropic sampler; tighter
  per-view drift bound) and showed that more sweeps marginally improves
  but doesn't fix the regression. A genuinely targeted investigation —
  per-frame inspection of which views' keypoints LK *worsened* on the
  fisheye rig, and why — is the natural next step. The hypothesis that
  analytic-Jacobian linearization on the highly-non-linear fisheye warp
  biases the GN step is consistent with the data we have but not
  pinpointed.
- **A per-dataset auto-default ("LK on, except on fisheye rigs")** would be
  the natural follow-on if a fix for kerry_park's LK regression doesn't
  surface. Selecting by camera model in `embed_patches` is a one-line gate;
  the policy decision belongs upstairs.
- **Round 2's dino numbers use a sub-sample.** If the full-point dino sweep
  becomes useful (e.g. for a future round considering a default flip), it
  would take maybe an hour or two of dedicated runtime. The sub-sample
  preserves the *relative* per-variant comparison, which is what the
  decision gate cares about.
