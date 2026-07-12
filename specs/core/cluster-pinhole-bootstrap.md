# Cluster pinhole bootstrap — experiment notes

**Status: idea-stage notes.** One experiment on one dataset
(`scripts/exp_pinhole_bootstrap.py`, seoul bull sculpture). Nothing here is
production design; numbers are from single runs.

## Problem

Starting from images and a cluster-patch `.matches` file (sift extraction →
cluster matching → cluster-patches), and using only a pinhole camera model:
how effectively can we bootstrap a coarse 3D reconstruction into a `.sfmr`
file — with no COLMAP solver, no pairwise two-view geometry, and no prior
intrinsics?

The input observations are the patch clusters' refined member positions
(the stored affine warp applied to the reference keypoint,
`x_m = A·x_ref + t`). Each cluster becomes one candidate track; members with
status `reference`/`kept` become its observations.

## Result (seoul bull, 17 images @ 270×480, ~150° orbit)

2,889 clusters (span ≥ 2), 7,193 observations. Reference: the incremental
COLMAP solve in the same workspace (944 points, 3,191 obs, SIMPLE_RADIAL
f = 332.6, k1 = +0.013).

A useful ceiling: triangulating every cluster against the *reference*
cameras explains 72% of the observations at < 2 px (median 0.29 px) — the
other ~28% of cluster members are junk that no cameras can explain. Per-image
ceilings range 54–83%.

The bootstrap (single run, ~3.3 min in scipy, of which ~2.7 min is a
5-candidate focal scan):

- **78.5% of all observations < 2 px** (kept-set rms 0.39 px, median
  0.17 px) — at the ceiling; the small excess over 72% is the pinhole
  absorbing part of the reference's k1.
- Cameras vs reference after similarity alignment: **rotation err mean 2.9°
  (max 5.1°), center err mean 1.6% of scene diameter (max 2.9%)**.
- Focal recovered **345.7 px vs 332.6 px reference** (+3.9%, plausibly the
  k1 trade-off) from a blind five-point grid over 0.55–1.6 × max(w, h).
- 2,319 points / 5,674 observations written as a valid
  `cluster_bootstrap` `.sfmr` (SIMPLE_PINHOLE, canonical convention,
  integrity OK).

## Method that survived

1. **Windowed affine factorization.** ALS (missing-data-tolerant, trimmed)
   weak-perspective factorization of overlapping 5-frame windows (stride 2),
   Tomasi–Kanade metric upgrade per window, both reflection hypotheses kept.
2. **Seed + incremental growth.** Seed a perspective solve on the
   densest window with a small fixed-focal BA (its inlier fraction also
   resolves the reflection: 93.7% vs 66.4% on the seed window). Then add one
   image at a time: constant-velocity pose extrapolation from the two
   nearest posed frames → trimmed-iteration pose-only resection → DLT for
   clusters that now have ≥ 2 posed views → a short global BA every 3
   images.
3. **Focal by outer scan, not by BA.** Run the whole growth at each
   candidate focal with f held fixed; the all-observation inlier fraction
   peaks near the true focal (41.8 / 67.4 / **78.2** / 66.6 / 54.6% across
   the grid). Release f only in the final BA of the winner — it then moves
   432 → 345.7 and stays put.
4. **Staged robust BA throughout**: trim gross outliers and behind-camera
   observations before each solve, re-triangulate every cluster from the
   refined cameras between rounds (re-admission), thresholds 50 → 12 → 4 px.

## What failed on the way (the actual findings)

- **A single global factorization over the orbit**: rotations 50–90° wrong
  before BA (the sequence spans ~150°; weak-perspective holds only over
  ~40° windows). Known from the earlier bootstrap experiments; confirmed.
- **Free-focal BA from a weak init escapes to the affine limit.** With
  median ~14 px init error, the reprojection residual decreases
  monotonically as f → ∞ at init, and a free-f BA slides 576 → 2435 px while
  *improving* its kept-set rms to 1.1 px — a self-consistent telephoto
  collapse fitting 48% of the data. The wrong reflection hypothesis
  survives this way too (reflection is unobservable in the affine limit).
  Hence the fixed-f outer scan.
- **Chaining independently-solved windows drifts.** Registering windows by
  similarity on shared cluster points accumulated 83° mean rotation error
  across the orbit — the sparse middle windows (~90–150 factorizable
  clusters) are 17–37° wrong even after a window BA, and everything after
  them inherits the error. Growth by resection against the *global*
  structure replaced it.
- **Plain robust resection fails from a one-frame-away init.** With ~100
  observations, an adjacent-frame init (~10° off, median residual ~50 px)
  and ~17% junk: an L2 warm-up is dragged by the junk's leverage, and
  soft_l1 has near-zero gradient when every residual starts as an
  "outlier" — both land ~21° wrong with ~1% inliers. Trimmed iterations
  (refit L2 on the best-fitting 60%, five times) reach 0.05° / 83% inliers
  on the same input. (COLMAP uses RANSAC P3P here for the same reason.)

## Open questions

- Other datasets: dino_dog_toy (photos, not a video orbit — the
  constant-velocity resection init and window structure assume sequence
  order), seattle (video), kerry (fisheye: pinhole-only should fail
  informatively — how gracefully?).
- The ~28% junk-observation floor is consistent with the kerry
  contamination floor from the grid-distortion experiments; per-member
  vetting signals (ZNCC, consistency residual) are stored in the `.matches`
  file and are not yet used to pre-filter here.
- Runtime is dominated by the focal scan (5 × full growth). A coarse-to-fine
  scan (2 candidates + golden-section refine), or reusing the seed across
  candidates, would cut most of it. All of this is throwaway scipy; a
  production version would be a Rust kernel.
- Whether the `cluster_bootstrap` `.sfmr` is good enough to seed `sfm solve`
  (as a triangulation/pose prior) or the planned cluster-level geometric
  verifier — not tried.
