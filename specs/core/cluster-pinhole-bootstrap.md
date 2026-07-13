# Cluster pinhole bootstrap — experiment notes

**Status: idea-stage notes.** Experiments with
`scripts/exp_pinhole_bootstrap.py` (first run: seoul bull sculpture; a
multi-dataset campaign log is at the bottom). Nothing here is production
design; numbers are from single runs.

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

(As first built for seoul this used the video's frame order — consecutive
windows and constant-velocity resection inits; the campaign version below
replaced every use of sequence order with cluster covisibility, with
identical results on seoul.)

1. **Seed-group affine factorization.** Candidate seed groups grow greedily
   from the strongest image-covisibility edge (shared-cluster counts),
   maximizing the minimum shared count against the group. ALS
   (missing-data-tolerant, trimmed) weak-perspective factorization +
   Tomasi–Kanade metric upgrade per group, both reflection hypotheses kept.
2. **Seed + incremental growth.** Seed a perspective solve on the best
   group with a small fixed-focal BA (its inlier fraction also resolves
   the reflection: 93.7% vs 66.4% on the seoul seed). Then repeatedly pose
   the next-best-view image (most observations of valid points): trimmed
   pose-only resection initialised from its most-covisible posed images →
   DLT for clusters that now have ≥ 2 posed views → a short global BA
   every few images.
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

## Campaign log (2026-07-12)

Running log of the multi-dataset campaign; entries appended as tried.

- **Campaign datasets** (picked from `test-data/` and `C:/DataSets` for
  contrast): seoul bull (17-frame 270×480 phone video, ~150° orbit —
  baseline), dino_dog_toy (85 unordered 2040×1536 photos, tabletop object),
  COLMAP south-building (128 photos 3072×2304, building + vegetation,
  reference = existing v4 glomap solve), DinoLedge subset (120 of 1196 4K
  video frames, stride 2, outdoor walk — reference = existing full-sequence
  solve, cross-workspace by image name), Swivel_Chair subset (106 of 632
  portrait-4K video frames — the existing April solve is .sfmr v1 with
  failed integrity, so a fresh subset solve is the reference), Kerry Park
  fisheye
  (24 frames, OPENCV_FISHEYE ground truth from the rig config — the
  deliberate failure-mode probe for pinhole-only). Also added
  `MAX_CLUSTERS = 10000` (keep highest-span clusters) — dino_dog_toy
  produces 122k clusters / 368k obs, 50× seoul, unusable in scipy BAs.
- **Covisibility grouping replaces sequence order.** All uses of the frame
  index (consecutive-frame windows, orbit growth order, constant-velocity
  resection init) replaced by shared-cluster counts: seed groups grow
  greedily from the strongest covisibility edge maximizing the *minimum*
  shared count against the group (mutual covisibility, not hub-and-spokes);
  growth picks the unposed image with the most observations of valid points
  (next-best-view); resection inits from the top-3 most-covisible posed
  images' poses (early-accept at >40% inliers). Focal scan now caps growth
  at ~20 images and only the winner grows fully. Seoul parity run: 78.4%
  inliers, 3.5° mean rotation, f 344.9 (vs 78.5% / 2.9° / 345.7 with the
  video-order version — same result, and the two seed groups it picked by
  covisibility alone were exactly the two ends of the orbit).

- **dino_dog_toy (85 unordered photos, 2040×1536), 18.7 min.** Focal scan
  peaked sharply and alone at 1428 (86.6% vs 59.5/61.8% neighbours on the
  20-image scan subset); released f = 1431.9 vs reference 1475.0 (−2.9%).
  Full growth posed all 85 images: 74.8% of 80,814 observations < 2 px
  (kept rms 1.34 px), 9,164 points. Camera errors vs the fresh incremental
  solve: rotation mean 17.9°, max 112.7°; centers mean 9.5%, max 41.2% —
  i.e. the bulk is right but a tail of cameras is badly wrong (unordered
  photos have weakly-connected views where trimmed resection can lock onto
  a locally-consistent wrong pose that the global 4 px trim then never
  revisits). Compare printout now includes medians + a >10° count to size
  that tail on later runs. Cluster cap engaged: kept 10,000 of 122,769 by
  span.

- **COLMAP south-building (128 photos, 3072×2304), 12.6 min.** The
  strongest run of the campaign: all 128 posed, **95.5% of 85,534
  observations < 2 px** (kept rms 0.53 px), rotation mean 2.33° / max 6.55°
  (zero cameras > 10°), centers mean 1.17%, f 2622.9 vs reference 2561.2
  (+2.4%). The focal scan again peaked decisively (95.9 vs 77.5/63.1%).
  Notable: seed groups picked by covisibility straddle non-consecutive
  file-order images ([24,25,26,125,126]) — the capture loops back, and the
  covisibility grouping finds it where a frame-index window could not.
  Cluster observations here are unusually clean (95.5% explainable —
  building texture localizes well).

- **DinoLedge subset (120 of 1196 4K video frames, stride 2, forward walk),
  15.1 min.** Reference = the existing full-sequence solve (cross-workspace
  by image name). All 120 posed, 84.6% of 88,890 observations < 2 px (kept
  rms 1.24 px), rotation mean 2.61° / max 3.20°, centers mean 0.36% / max
  0.60%. Focal essentially exact: 2751.3 vs 2746.8 (**+0.16%**) — this
  camera is nearly distortion-free, so pinhole-only has no k1 to absorb and
  the recovered f matches. Walking-forward motion (not an orbit) works the
  same as orbits under the covisibility machinery.

- **Swivel_Chair subset (106 of 632 portrait-4K (2160×3840) video frames,
  stride 6, indoor object orbit), 12.7 min.** The existing April solve is
  .sfmr v1 with failed integrity — reference is a fresh incremental solve
  on the same subset. Sharp focal peak again (84.1% at 2688). 78.9% of
  80,615 observations < 2 px, 9,609 points; f 2709.9 vs 2740.4 (−1.1%).
  The most accurate cameras of the campaign: rotation mean **0.26°** / max
  0.70°, centers mean 0.24% / max 0.60%. A dense indoor orbit of a
  texture-rich object at video frame rate is the easy case. (Logged as
  1080p during the campaign; the v4 `image_dims` array caught the error.)

- **Kerry Park fisheye (24 frames 480×480, OPENCV_FISHEYE fx = 129.1),
  6.3 min — the failure-mode probe.** Pinhole-only degrades gracefully
  rather than crashing: 43.4% of 13,217 observations < 2 px, 2,471 points,
  f settling at 340.7 (not comparable to the equidistant fx). Inliers by
  radial band: 44/46/43% out to 0.5·rmax (matching this capture's ~50%
  contamination floor even at the center), 34% (median 9 px) at 0.5–0.7,
  and 5% (median 164 px) past 0.7. So the bootstrap silently keeps the
  central ~half-field where pinhole ≈ equidistant and drops the rim —
  graceful and radially ordered, which is exactly the structure a
  camera-correction stage (or a center-out unlock) could pick up from.

  > _Status (2026-07-12): Falsified by visual inspection — the kerry
  > reconstruction is a complete geometric failure._ The GUI shows a
  > tangle; diagnostics confirm: consecutive-frame rotation deltas of a
  > steady walking capture swing 2.5°–157°, camera spread is 10× the
  > scene scale, and camera 0 sees the *median* 3D point at negative
  > depth (most structure behind the camera — a mirror/degenerate
  > collapse). The paragraph above stands as a lesson in metric
  > circularity, not as a result: kerry has no pose reference, so its
  > inlier and radial-band numbers were computed against the bootstrap's
  > own broken cameras — a self-referential score that locally-consistent
  > wrong geometry passes. Pinhole-only on this fisheye fails outright,
  > and the pipeline's internal metrics cannot detect it. Cheap
  > self-diagnostics that would have caught it: cheirality fraction over
  > all observations, and per-camera structure-depth sign stats.

- **Visual inspection of dino_dog_toy (2026-07-12).** The misregistration
  tail is visible as a partial duplicate ("echo") of the dino in the point
  cloud: the >10° cameras are not randomly wrong but form a coherent
  wrongly-registered subset that re-triangulates its own copy of the
  object — consistent with resection locking onto a locally-consistent
  wrong pose and the global 4 px trim then keeping the echo's
  self-consistent observations.

### Cross-dataset summary

| dataset | input | imgs | obs | inlier<2px | rot err mean/max (deg) | center err mean/max (% diam) | f vs ref | time |
|---|---|---|---|---|---|---|---|---|
| seoul bull | 270×480 video orbit | 17 | 7,193 | 78.4% | 3.46 / 6.00 | 1.8 / 3.2 | +3.7% | 4 min |
| dino_dog_toy | 2040×1536 photos, unordered | 85 | 80,814 | 74.8% | 17.9 / 112.7 | 9.5 / 41.2 | −2.9% | 19 min |
| south-building | 3072×2304 photos | 128 | 85,534 | 95.5% | 2.33 / 6.55 | 1.2 / 4.5 | +2.4% | 13 min |
| DinoLedge (subset) | 4K video, forward walk | 120 | 88,890 | 84.6% | 2.61 / 3.20 | 0.4 / 0.6 | +0.16% | 15 min |
| Swivel_Chair (subset) | portrait-4K video orbit | 106 | 80,615 | 78.9% | 0.26 / 0.70 | 0.2 / 0.6 | −1.1% | 13 min |
| Kerry fisheye | 480×480 fisheye | 24 | 13,217 | 43.4%† | (no reference) | — | n/a | 6 min |

† Self-referential (no pose reference); visual inspection shows a complete
geometric failure — see the kerry status note above.

Cross-dataset observations:

- The fixed-f focal scan peaked **decisively and uniquely on every
  dataset** — the inlier-fraction-at-fixed-f signal appears robust across
  scene types, resolutions, and motion patterns.
- The recovered focal lands within ±4% of the reference everywhere, and the
  deviation tracks the reference's k1 (DinoLedge, nearly distortion-free:
  +0.16%; seoul, largest k1 relative to f: +3.7%).
- The one weak spot is **unordered photo sets** (dino_dog_toy): a tail of
  weakly-connected cameras resects onto locally-consistent wrong poses that
  the global trim never revisits. Videos and the loop-closing
  south-building set have no such tail (0 cameras > 10° on all four).
- Wall-clock is dominated by the 5-candidate focal scan and the scipy BAs;
  the per-dataset ~13–19 min is throwaway-prototype speed, not a statement
  about the method.

## Open questions

- The dino_dog_toy misregistration tail: weakly-connected views in
  unordered photo sets need either a resection acceptance gate (reject and
  retry later when more structure exists), a re-resection pass after the
  final BA, or RANSAC P3P. Everything else in the campaign has no tail.
- The junk-observation floor (~5–25% by dataset, ~50% on kerry) matches the
  contamination floor seen in the grid-distortion experiments; per-member
  vetting signals (ZNCC, consistency residual) are stored in the `.matches`
  file and are not yet used to pre-filter here.
- Runtime is dominated by the focal scan (5 × capped growth) and the scipy
  BAs. A coarse-to-fine scan (2 candidates + golden-section refine) would
  cut most of it. All of this is throwaway scipy; a production version
  would be a Rust kernel.
- Whether the `cluster_bootstrap` `.sfmr` is good enough to seed `sfm solve`
  (as a triangulation/pose prior) or the planned cluster-level geometric
  verifier — not tried.
- Self-diagnostics: kerry produced a completely broken reconstruction
  while its internal inlier metric read 43% — with no external reference,
  the bootstrap cannot currently tell success from locally-consistent
  failure. Cheap candidates that would have caught it: cheirality fraction
  over all observations, per-camera structure-depth sign statistics, and
  pose-path coherence on sequential captures.
