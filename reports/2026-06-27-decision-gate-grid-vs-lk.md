# Decision Gate: Supersampled Grid vs LK Subpixel Refinement

_Date: 2026-06-27._ The final phase of the keypoint-localization implementation
plan ([`reports/2026-06-26-keypoint-localization-implementation-plan.md` on the
plan branch][plan]). Phases 3A (AVX2 grid kernel, #139) and 3B (analytic
Jacobian + per-move incremental consensus, #140 + #142) have landed. This report
takes the optimized variants head-to-head on the four checked-in datasets and
recommends a production default for the new `embed_patches(..., subpixel=...)`
knob.

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

Both `lk` and `lk_per_move` are now **opt-in** behind `subpixel="lk"` and
`subpixel="lk_per_move"` (CLI: `sfm embed-patches --subpixel lk`); the
supersampled grid is opt-in behind `--search-resolution-multiplier 2`. The data
is documented here so the next round (different datasets, or a per-dataset auto
mode) has the numbers it needs.

The "do nothing" verdict is the boring answer, but it's what the data supports —
two of four code paths landed (LK is the cheapest at ~1.5× and has a real
regression case, the grid is a positive uniform but a multiplicative cost), the
production knobs are wired and tested, and the production default doesn't move
until a clearer signal emerges.

## Variant matrix

| Variant | `search_resolution_multiplier` | `subpixel` | Notes |
|---|---|---|---|
| `baseline` | 1.0 | `"none"` | Current production behaviour |
| `grid_m2` | 2.0 | `"none"` | Supersampled grid only |
| `grid_m3` | 3.0 | `"none"` | Even more aggressive grid |
| `lk_per_sweep` | 1.0 | `"lk"` | LK + per-sweep consensus, `max_outer_sweeps=1` |
| `lk_per_move` | 1.0 | `"lk_per_move"` | LK + per-move (Gauss-Seidel), `max_outer_sweeps=5` |
| `grid_m2_then_lk` | 2.0 | `"lk"` | Supersampled grid + LK refinement on top |

## Methodology

For each (dataset × variant) the script
(`scripts/measure_subpixel_decision_gate.py`) runs the full
`embed_patches(...)` pipeline end-to-end and captures:

- **Wall time (`wall_secs`)** — time from before `embed_patches` until after it
  returns. Includes image load via `read_workspace_image` (~~~10s for
  dino_dog_toy~~, negligible for the others). Single sample per (dataset,
  variant); no warm-up pass.
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
- **Mean / median keypoint shift vs baseline (`mean_shift_px`,
  `median_shift_px`)** — per-observation distance in source-image px between
  the variant's stored keypoint and the baseline's, joined by **(world
  position, image index)** (the compaction renumbers survivor point indices,
  so joining by output index is wrong; joining by the source world position
  is the invariant key since compaction never re-bundles).

### What the metrics CAN and CAN'T tell you

- ECC is computed against a consensus rebuilt at the **stored** keypoints, so
  it rewards self-consistency — a variant that shifted all views by the same
  amount could score well even if it walked off the truth together. The
  cross-view shift bounds are what stop this in practice (the localizer's
  `max_shift_px` gate and the LK refiner's `max_offset_px`), but ECC alone
  doesn't tell you ground truth. Three-decimal-place ECC deltas between
  variants are small relative to the dataset-to-dataset variation (the
  textured close-range seoul_bull scores ~0.94, the fisheye kerry_park ~0.92).
- The keypoint shift vs baseline measures how much each variant *changed* the
  keypoint locations relative to the no-refinement path — it's a magnitude
  signal, not a "right vs wrong" signal. A high shift with high ECC means
  "the refiner moved things and they agree better with the rebuilt
  consensus"; a high shift with **lower** ECC (which happens on kerry_park
  for the LK variants) means the refiner moved things and they agree
  **worse** — a real regression.
- Wall time is single-sample, no statistics. Numbers within ~10% of each
  other should not be treated as significant. The big multiplicative gaps
  (`m=3` is ~17-26× baseline) are real signal.
- The localizer's `loo_zncc` field and the refiner's `scores` field are not
  directly comparable (different consensus rebuilds, different normalizations
  in different contexts) — that's why the ECC metric is recomputed
  end-of-pipeline against a single objective for every variant. A simpler
  measurement of "report whichever field is set" would have been an
  apples-to-oranges comparison and was rejected.
- `dino_dog_toy` (85 images at 2040×1536) was not run due to time —
  extrapolating from the slower variants on seattle_backyard
  (350 s for `grid_m3`), it would take ~30-60 min for the full sweep. The
  three other datasets cover the dataset-size and texture-difficulty
  spectrum (17 small close-range, 26 medium outdoor, 48 fisheye rig).

## Per-dataset results

### seoul_bull_sculpture (17 images @ 270×480, close-range, textured)

| Variant | Wall (s) | Points | Mean ECC | ΔECC | Mean shift (px) | Median shift (px) |
|---|---:|---:|---:|---:|---:|---:|
| baseline        |  2.15 | 1029 | 0.9418 |  —      | 0.000 | 0.000 |
| grid_m2         | 13.84 | 1028 | 0.9427 | +0.0009 | 0.326 | 0.206 |
| grid_m3         | 52.35 | 1028 | 0.9429 | +0.0011 | 0.301 | 0.201 |
| lk_per_sweep    |  3.20 | 1029 | 0.9445 | **+0.0027** | 0.100 | 0.063 |
| lk_per_move     |  6.37 | 1029 | 0.9458 | **+0.0040** | 0.174 | 0.101 |
| grid_m2_then_lk | 14.94 | 1028 | 0.9447 | +0.0029 | 0.335 | 0.212 |

LK is the clear winner here. The dataset has rich texture and short baselines,
the regime LK was designed for. Both LK variants slightly outperform the
supersampled grid at a fraction of the cost.

### seattle_backyard (26 images @ 360×640, outdoor mixed)

| Variant | Wall (s) | Points | Mean ECC | ΔECC | Mean shift (px) | Median shift (px) |
|---|---:|---:|---:|---:|---:|---:|
| baseline        |  13.44 | 3295 | 0.9266 |  —      | 0.000 | 0.000 |
| grid_m2         |  90.41 | 3293 | 0.9280 | +0.0014 | 0.202 | 0.140 |
| grid_m3         | 350.21 | 3293 | 0.9282 | +0.0016 | 0.195 | 0.131 |
| lk_per_sweep    |  19.46 | 3295 | 0.9273 | +0.0007 | 0.070 | 0.041 |
| lk_per_move     |  35.42 | 3295 | 0.9276 | +0.0010 | 0.106 | 0.061 |
| grid_m2_then_lk |  96.50 | 3293 | 0.9280 | +0.0014 | 0.208 | 0.144 |

Closer call — `grid_m2` edges out the LK variants on mean ECC, at ~5× the cost.
`grid_m2_then_lk` is no better than `grid_m2` alone, suggesting the grid found
the basin LK was going to find anyway. Pure LK is the cheapest improvement but
also the smallest.

### kerry_park (48 images @ 480×480, fisheye rig)

| Variant | Wall (s) | Points | Mean ECC | ΔECC | Mean shift (px) | Median shift (px) |
|---|---:|---:|---:|---:|---:|---:|
| baseline        |  4.97 | 786 | 0.9214 |  —       | 0.000 | 0.000 |
| grid_m2         | 25.11 | 786 | 0.9217 | +0.0003  | 0.243 | 0.176 |
| grid_m3         | 97.25 | 786 | 0.9221 | +0.0007  | 0.230 | 0.170 |
| lk_per_sweep    |  5.87 | 786 | 0.9183 | **−0.0031** | 0.133 | 0.053 |
| lk_per_move     | 12.18 | 786 | 0.9190 | **−0.0024** | 0.202 | 0.104 |
| grid_m2_then_lk | 28.29 | 786 | 0.9183 | −0.0031  | 0.290 | 0.186 |

The regression that determines the verdict. **Both LK variants score
*worse* than the baseline** on the fisheye rig. The grid variants stay
positive but small. The composition `grid_m2_then_lk` lands at the LK
score, not the grid score — once LK runs, it pulls the keypoints back
toward whatever local optimum it was going to find regardless of where
the grid put them, and on this dataset that optimum is worse than the
discrete localizer's pick.

A plausible cause: the fisheye projection is highly non-linear; the LK
refiner uses an analytic Jacobian (Phase 3B) that linearizes the warp
about the seed, and on a wide-angle camera the linearization error
relative to the offset budget (`max_offset_px = 2.0`) is large enough to
bias the GN step in the wrong direction. This isn't a bug to fix in this
phase — it's the kind of dataset-dependence the decision gate exists to
surface.

The point count is identical for every variant on this dataset (786) —
the LK regression isn't a "the refiner dropped more points" failure mode,
it's the refiner doing its job but converging to a slightly worse
keypoint than the discrete localizer found, on fisheye geometry where
that's possible.

## Aggregate verdict

| Variant | Mean ΔECC (3 datasets) | Cost ratio (geo-mean) | Always non-regressing? |
|---|---:|---:|---|
| `baseline`        |  0      | 1.0×  | (the baseline) |
| `grid_m2`         | +0.0009 | ~5.8× | yes |
| `grid_m3`         | +0.0011 | ~21×  | yes |
| `lk_per_sweep`    | +0.0001 | ~1.4× | **no — kerry_park −0.0031** |
| `lk_per_move`     | +0.0009 | ~2.5× | **no — kerry_park −0.0024** |
| `grid_m2_then_lk` | +0.0004 | ~6.4× | **no — kerry_park −0.0031** |

`grid_m3` has the best mean ΔECC but at 17-26× baseline cost — not worth
defaulting to. `grid_m2` is the only variant with a uniformly positive ΔECC
across all three datasets, but its ~5.8× cost is steep for a 0.001 ECC
improvement. The LK variants are cheap and excellent on close-range textured
data but regress on the fisheye rig.

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
  strong enough to flip it).

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

## Open follow-ups for the lead

- **The plan on `claude/quirky-heisenberg-2vlpmp` should be marked complete**
  (the "Decision gate" section). I don't have that branch checked out, so I
  haven't touched it.
- **`dino_dog_toy` was skipped.** A full sweep would take ~30-60 minutes;
  worth running before any production-default flip, since dino is the largest
  / sharpest of the four datasets and might shift the verdict if its LK
  behaviour is closer to seoul_bull's (textured close-range, big LK win) than
  to kerry_park's (fisheye, LK regression).
- **The kerry_park LK regression deserves investigation.** It's plausible
  enough (fisheye + analytic-Jacobian linearization), but if the cause is
  fixable (e.g. tighter `max_offset_px` for high-distortion cameras), the LK
  default flip becomes feasible. I have not investigated this — the gate's
  charter is to measure, not to fix.
- **A per-dataset auto-default ("LK on, except on fisheye rigs")** would be
  the natural follow-on if a fix for kerry_park's LK regression doesn't
  surface. Selecting by camera model in `embed_patches` is a one-line gate;
  the policy decision belongs upstairs.
