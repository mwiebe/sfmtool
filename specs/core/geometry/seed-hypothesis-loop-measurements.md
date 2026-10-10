# Seed Hypothesis Loop Measurements

This file records measurements behind [seed-hypothesis-loop.md](seed-hypothesis-loop.md), the seed stage that develops a set of candidate reconstructions from a cluster-patches `.matches` file and records which candidate it ranks first. The measurements bear on how stable that set and that first candidate are when the input file changes slightly, which stage of the exploration turns a small change of the input into a different candidate, and which of those stages behave as defects rather than as a response to the evidence. The changes they led to are described in [seed-hypothesis-loop.md](seed-hypothesis-loop.md) § "Member order", § "Capture-level measurements", § "Choice among seed groups", § "Focal scan" and § "Rank", and in [rotation-locked-resection.md](rotation-locked-resection.md).

## Pick stability under small changes to the cluster file (2026-10-09)

**Question.** On `KerryPark480`, two cluster-patches files that differ in the shapes of a few dozen of about 16,000 kept members gave the seed a first candidate that failed against the ground truth and one that passed ([cluster-patch-refinement-measurements.md](../patch/cluster-patch-refinement-measurements.md#subset-with-the-loop-as-the-default-2026-10-09)), and the gates at the refined shape changed which candidates were committed again ([the next section there](../patch/cluster-patch-refinement-measurements.md#gates-at-the-refined-shape-2026-10-09)). Which stage turns such a change into a different set, and is that a defect (a tie broken by order, a hard threshold the baseline sits next to, a result that depends on the order of the members) or the seed responding to evidence that changed?

**Method.**

- **Code.** Branch `bootstrap-core-migration` at `c69ed0ac`, with no change to the seed code. Seed entry `scripts/exp_fast_seed.py`.
- **Files.** For `SeoulBull` and `KerryPark480`, the clusters files of [Seed stage on the ground-truth entries](../patch/cluster-patch-refinement-measurements.md#seed-stage-on-the-ground-truth-entries) were refined three ways with `sfm cluster-patches --patch-size 12`: `--piecewise` at the defaults (the *base* file, identical in every member column to the file the [gates section](../patch/cluster-patch-refinement-measurements.md#gates-at-the-refined-shape-2026-10-09) seeded from), `--no-piecewise` (the cascade's shapes), and `--piecewise --no-regate-at-refined-shape --max-capped-cells 9` (no gate at the refined shape). The *movers* are the kept members whose position or shape differs between the base and the cascade file: 5 of 2,820 kept members on `SeoulBull`, 56 of 15,609 on `KerryPark480`.
- **Perturbations.** Twelve files per entry, each derived from the base file:

  | class | files | change |
  |---|---|---|
  | base | 2 | the base file, seeded twice |
  | random drop | 3 | 0.3% of kept members, drawn with seeds 1001 to 1003 (8 on `SeoulBull`, 47 on `KerryPark480`), marked `rejected_unlocalizable_refined` |
  | drop movers | 1 | every mover marked `rejected_unlocalizable_refined` |
  | drop non-movers | 3 | as many random non-movers as there are movers, seeds 2001 to 2003 |
  | swap movers | 1 | every mover given its cascade position and shape |
  | ungated | 1 | the file with no gate at the refined shape: the 133 and 377 members the gate refuses are kept |
  | reorder | 2 | the members of every cluster permuted, seeds 3001 and 3002, with the reference index remapped; the content is unchanged |

- **Seed.** `pixi run -e dev python scripts/exp_fast_seed.py <workspace copy>` with `SFMTOOL_RELAX=0`, `SFMTOOL_SEED_RUNG1=3000`, `SFMTOOL_SEED_EVAL` and `SFMTOOL_PINHOLE_RAW_VOTE` unset, and the instrumentation switches `SFMTOOL_SEED_EVO_DUMP`, `SFMTOOL_STAGE_DUMPS` and `SFMTOOL_RESECT_TRACE` set. None of them changes a decision, and the two base runs of each entry are identical.
- **Scoring.** Every finite candidate is scored against the entry's ground truth as in [Seed stage on the ground-truth entries](../patch/cluster-patch-refinement-measurements.md#seed-stage-on-the-ground-truth-entries): it passes with a median rotation error under 1°, a median centre error under 5% of the capture's extent and an equivalent-focal error under 5%.
- **Attribution.** For each run, the committed candidates were compared with the base run's in commit order, and the first stage at which a candidate differs was recorded: the seed group (probe), the posed frames after the widen ladder, the attempt's kept count and the commit bar, the focal scan's peak and winner, the release (focal and spline verdict), and the flags and qualification. Where a stage differed, the decision was reproduced offline from the same module functions on the run's own workspace: the per-pair focal votes of `geometry.estimate_intrinsics`, the rotation core's gauge pair from `_best_ray_pair`, and the rotation-locked growth of `rotation_core_rays` with each failed resection's survivor count.

**Result: the pick is the first committed candidate in every run, and what changes is that candidate's own exploration and the passes after it.**

The rank never reorders the set: `ladder_first` is `0` in all 26 runs, either because `h00` qualifies or because nothing qualifies and the rank falls back to `h00`. The first candidate on both entries is the first pass's, so the pick changes when the first pass produces a different candidate or when that candidate's qualification changes. The set changes when the coverage claim of an earlier candidate changes, because each later pass explores the complement the earlier claims leave.

`SeoulBull`. The base `h00` is seeded from frames 0 to 4, poses 8 frames, qualifies and passes. The vote is the pinhole pairwise vote.

| class | run | vote f (change) | `h00` seed | `h00` qualified | `h00` pass | qualified in set | first stage that differs |
|---|---|---|---|---|---|---|---|
| base | repeat | 309.6 | same | yes | PASS | 2 | none |
| random drop | s1 | 317.5 (+2.5%) | same | yes | PASS | 2 | `h01` focal scan |
| random drop | s2 | 280.7 (−9.3%) | same | no (`edge_scan`) | PASS | 0 | `h00` flags |
| random drop | s3 | 258.3 (−16.6%) | frames 12 to 16 | no (`edge_scan`) | PASS | 0 | `h00` probe |
| drop movers | | 313.2 (+1.2%) | same | yes | PASS | 2 | `h00` release (spline accepted) |
| drop non-movers | s1 | 316.6 (+2.3%) | same | yes | PASS | 2 | `h01` focal scan |
| drop non-movers | s2 | 264.7 (−14.5%) | frames 12 to 16 | no (`edge_scan`) | **fail** (focal −9.4%) | 0 | `h00` probe |
| drop non-movers | s3 | 283.0 (−8.6%) | same | no (`edge_scan`) | PASS | 0 | `h00` flags |
| swap movers | | 309.6 | same | yes | PASS | 2 | none |
| ungated | | 298.6 (−3.6%) | same | yes | PASS | 1 | `h00` release (spline accepted) |
| reorder | s1, s2 | 309.6 | same | yes | PASS | 2 | `h02` release |

`KerryPark480`. The base `h00` is the rotation core of the first pass, poses 28 frames, qualifies and passes; none of the six qualified complement candidates passes. The capture routes to the equidistant model, whose verdict focal is 138.3 px in every run.

| class | run | `h00` posed | `h00` pass | complement candidates passing | seed frames of `h01` to `h06` changed | first stage that differs |
|---|---|---|---|---|---|---|
| base | repeat | 28 | PASS | 0 | none | none |
| random drop | s1 | 28 | PASS | 1 | 5 of 6 | `h01` widen |
| random drop | s2 | 28 | PASS | 0 | 6 of 6 | `h00` probe (gauge pair) |
| random drop | s3 | 29 | PASS | 1 | 6 of 6 | `h00` probe |
| drop movers | | 23 | PASS | 1 | 6 of 6 | `h00` probe (locked resection) |
| drop non-movers | s1 | 15 | PASS | 2 | 6 of 6 | `h00` probe (locked resection) |
| drop non-movers | s2 | 28 | PASS | 0 | 5 of 6 | `h01` widen |
| drop non-movers | s3 | 19 | **fail** (centre 10.4%) | 3 | 6 of 6 | `h00` probe (locked resection) |
| swap movers | | 28 | PASS | 0 | 6 of 6 | `h00` probe (gauge pair) |
| ungated | | 26 | PASS | 1 | 6 of 6 | `h00` probe |
| reorder | s1 | 28 | PASS | 0 | 3 of 6 | `h03` widen |
| reorder | s2 | 28 | PASS | 1 | 4 of 6 | `h02` widen |

**Where the change enters.**

1. **The pairwise focal vote (`SeoulBull`).** Each image pair votes one Bougnoux focal from one RANSAC fundamental-matrix fit with seed 0, and the vote is the median of 13 such votes whose log-focal interquartile range is 0.33 to 0.43. Dropping 8 of 2,820 members changes a pair's correspondence list, so the fixed-seed RANSAC draws different samples, and single pair votes move by up to a factor of two: pair (6, 7) from 499 to 240 px, pair (1, 2) from 572 to 418 px. On the unchanged base file, changing only the RANSAC seed from 0 to 9 moves the vote between 251.2 and 323.4 px, the same range the dropped members produced (258.3 to 317.5). The vote's movement is the spread of one RANSAC draw per pair, which seed 0 hides on any one file; it is not a change in the geometry. Reordering the members leaves every pair vote unchanged.
2. **The probe focal and the commit bar (`SeoulBull`).** The probe grows each seed group to the core cap of 8 images (`SCAN_CAP`) at the vote's focal, the widen ladder adds no frame on this capture, and the commit bar requires 8 kept frames. The bar therefore equals the cap, and one frame lost in growth fails it. On the base file with the probe focal forced to 264.7 px (`SFMTOOL_F0`), group 0 to 4 grows to 7 frames, fails the bar, and frames 12 to 16 seed `h00`; forced to 258.3 px it grows to 8. The growth is not monotonic in focal, so the vote's noise reaches the pick through this bar.
3. **The scan grid and the flags read off it (`SeoulBull`).** The five-point scan grid is centred on the vote at a ratio of 1.15. When the vote is low by 9% or more, the true focal (about 336 px) sits at the grid's top point and the inlier fraction rises to it monotonically, for example 67.5, 70.7, 76.5, 83.2 and 83.6% from 214 to 374 px with the vote at 283.0. `edge_scan` fires and `h00` is disqualified, although it releases the same 338.9 px and passes against the ground truth. The scan winner also follows the vote: the candidate nearest the vote takes a refit slot when it is within 5 points of the leader, and wins when the two refits are within 5 points. On `h01`, the vote point scored 79.4% against a leader of 85.5% with the vote at 309.6 (outside the 5-point band) and 82.6% with the vote at 317.5 (inside it), so the winner moved from 356.0 to 317.5 px and the release from 356.0 px (focal +5.9%, fail) to 329.7 px (pass).
4. **The rotation-locked growth of the rotation core (`KerryPark480`).** The first pass seeds from the rotation core: a gauge pair, a rotation spanning tree, then rotation-locked translation resections in order of support. `resect_translation` fails when fewer than 10 observations survive any of its trim rounds at 8 px, and a failed image is removed from the skeleton for the rest of the core. In the base run the resections keep 13 to 53 survivors, and image 6 keeps 13 of 52. With 56 non-movers dropped (s1) image 6 keeps 9, fails, and the 11 images resected after it lose the structure it would have added and fail in turn: the core poses 9 frames instead of 22, and `h00` 15 instead of 28. In s3 image 28 fails with 12 survivors at the last round, below 10 at an earlier round, and `h00` poses 19 frames and fails against the ground truth. With the movers dropped, image 31 fails with 31 survivors at the last round. The bar is read on an intermediate round rather than on the converged set, and one failure removes an image permanently.
5. **The rotation core's gauge pair (`KerryPark480`).** The gauge is the pair with the largest key (parallax at least 1.0°, cheiral count, parallax). Pair (0, 1) leads with 177 cheiral correspondences against 130 for the next pair. In the swapped-mover file the seed-0 essential-matrix RANSAC on pair (0, 1) converges to a solution with 0.71° of parallax; seeds 1 to 11 give 1.48°. Below the 1.0° floor the pair drops behind every pair that clears it, and (3, 4) becomes the gauge. In random drop s2, every seed gives 0.69° once one of 227 correspondences is gone. The pair's translation direction is poorly constrained, and its parallax reading falls on either side of a hard floor.
6. **The choice among starved outcomes (`KerryPark480`).** No complement pass reaches the commit bar's 60% reach on its own admission, so every pass returns the deferred outcome with the most kept frames, keeping the first one tried on a tie. Kept frames saturate at 14, the core cap of 8 plus the widen ladder's 6 rungs, and two or three seed groups reach 14 in passes 1 to 5 of the base run. The candidate each pass commits is the first group to reach 14 in covisibility-group order, and the group order is a greedy maximum over shared-cluster counts that one dropped member, or a different claim from an earlier pass, can change.
7. **Member order.** Reordering the members of each cluster without changing them changes results. `core_parallax` measures each point's widest ray angle from its first listed observation rather than over every pair of its views (12.66° against 11.02° on `SeoulBull`'s first probe). The ladder's far-field reading and the rotation-only layers take their reference ray from the first observation row of each cluster (460 against 491 points at infinity on `SeoulBull`). The bundle adjustments stop at a fixed evaluation count (25 to 60) and sum their rows in observation order, so on `KerryPark480` an off-peak scan point reads 37.8% or 39.1%, and in pass 2 the widen ladder admits 14 or 10 frames from the same probe. On `SeoulBull` the two reflection hypotheses of one probe read 47.5% and 41.2%, or 47.5% and 47.6% once reordered, and the other one is kept. Through rule 6, reorder s2 commits a different candidate from pass 2 on.
8. **The complement queue.** Every later pass explores what the earlier candidates leave unclaimed, so any change to `h00`'s posed set or structure changes every complement after it. This is the loop's design, and it turns every difference in rules 1 to 7 into a different set.

**Classification.**

| kind | stage | example |
|---|---|---|
| (a) hard threshold next to the baseline | commit bar of 8 kept frames equal to the core cap of 8 | `SeoulBull` group 0 to 4 at 7 frames when probed at 264.7 px |
| (a) | `edge_scan` and the 5-point refit band, read on a grid centred on the vote | `SeoulBull` `h00` disqualified in 3 of 12 files with an unchanged release |
| (a) | locked resection's 10 survivors on any trim round, with a failed image removed for good | `KerryPark480` image 6: 13 survivors in the base file, 9 in drop non-movers s1, then 11 more failures |
| (a) | gauge pair parallax floor of 1.0° | pair (0, 1) at 1.48° or 0.69 to 0.71° |
| (b) tie or order | first deferred outcome kept among equal kept counts | `KerryPark480` passes 1 to 5: 2 or 3 groups at the 14-frame ceiling |
| (b) | one RANSAC draw with a fixed seed per pair vote | `SeoulBull` vote 251 to 323 px over seeds 0 to 9 on one file |
| (b) | first observation as reference in `core_parallax`, the far-field reading and the rotation-only layers; observation-order sums in iteration-capped adjustments | reorder runs on both entries |
| (c) change in the evidence | none among the 16 reproduced | 4 of the runs not reproduced (`KerryPark480` random drop s1 and s3, drop non-movers s2, ungated) first differ at the widen gate or the rotation core's posed count, which read geometry the dropped members do change; they may hold a genuine change, which then reaches the set through rule 6 and the complement queue |

**What this decides.**

- **The amplification is upstream of the rank.** The rank's rule (first qualified, commit order) reorders nothing on these two entries; the first candidate's own exploration changes, and the complement queue carries the change through the rest of the set.
- **It is mostly a defect, not the evidence.** 21 of the 22 perturbed runs differ from the base in at least one committed candidate; only the swap of `SeoulBull`'s 5 moved shapes changes nothing. The first stage that differs was located in all 21. In 16 the mechanism was reproduced (the vote over RANSAC seeds, the forced probe focal, the scan readings, the gauge pair over seeds, the survivor counts of the locked resections, and the member-order readings), and each is a hard threshold the baseline sits within a few counts or percent of, a tie broken by try order, a single fixed-seed RANSAC draw, or member order. The other 5 (`KerryPark480` random drop s1 and s3, drop non-movers s2, ungated and reorder s1) first differ at the rotation core's posed count or at a complement pass's widen, stages governed by rules 4, 6 and 7, but their individual mechanism was not reproduced; reorder s1 changes no member and is order dependence by construction. None of the 16 needed a load-bearing member to explain it. The swap of the 56 moved shapes on `KerryPark480` changes the outcome only through rule 5, a RANSAC draw on a pair whose other seeds agree with the base.
- **What to change**, as this section concluded: a vote measured over several RANSAC draws, a scan grid and probe focal that do not move with every change of the vote, a commit bar below the cap, choice among outcomes by a continuous score with a tie-break on content, a locked resection that judges its converged set and retries a failed image, and member order fixed at load. The perturbation suite here is the acceptance test for those changes, rerun in the two sections below.

## Pick stability after the deterministic fixes (2026-10-09)

**Question.** Three of the changes the [section above](#pick-stability-under-small-changes-to-the-cluster-file-2026-10-09) led to remove causes that are not evidence at all: member order fixed at load, each pair vote and each ray-space pair init read as the median of several RANSAC draws, and the commit bar's kept count one below the core cap. With them built, does the [perturbation suite](#pick-stability-under-small-changes-to-the-cluster-file-2026-10-09) give bit-identical products under member reordering, which of its runs still change the first candidate, and does either capture's base `h00` stop passing its ground truth?

**Method.**

- **Code.** Branch `bootstrap-core-migration` with the three changes: [seed-hypothesis-loop.md](seed-hypothesis-loop.md) § "Member order", § "Capture-level measurements" (draws, `SFMTOOL_VOTE_DRAWS` at its default 5) and § "Rank" (commit bar `min(8, cap - 1)`, so 7 kept frames). Everything else in the seed is as in the section above.
- **Files, seed runs and scoring.** The same 26 files per the section above (12 per capture plus a repeat of each base file), the same seed command, environment and instrumentation, and the same ground-truth scorer and pass bar.
- **Identity.** Two products are called bit-identical when their manifests match with the run stamp and the timing fields removed, and every entry of every release file matches except `written.json`, `metadata.json` and `content_hash.json`, which carry the write time.
- **Attribution.** As in the section above, with the rotation core's gauge pair and locked resections reproduced offline from the module functions on each run's own workspace (`_best_ray_pair` now reads its median over draws).
- **Vote spread.** On the unchanged `SeoulBull` base file, the vote at `K` draws over ten disjoint seed sets (base seeds `0, K, …, 9K`), for `K` in 1, 5, 9, 17 and 33, measured through `geometry.estimate_intrinsics` on the full admission in the canonical member order.

**Result: reordering is solved and the vote no longer follows one draw; the first candidate changes on as many drop files as before, through causes the remaining proposals address.**

Vote spread on the `SeoulBull` base file:

| `K` | vote range over 10 seed sets (px) | log range | log s.d. | pool log-IQR (median) |
|---|---|---|---|---|
| 1 | 251.2 to 323.4 | 0.253 | 0.096 | 0.397 |
| 5 | 245.0 to 298.0 | 0.196 | 0.063 | 0.366 |
| 9 | 254.3 to 289.9 | 0.131 | 0.037 | 0.338 |
| 17 | 247.0 to 281.2 | 0.130 | 0.035 | 0.372 |
| 33 | 258.5 to 276.3 | 0.067 | 0.023 | 0.390 |

`K = 1` reproduces the 251 to 323 px of the section above. The spread over seed sets falls roughly as the square root of `K`, and at every `K` it is under the pool's own log-IQR, so that criterion does not choose `K`. The vote settles near 268 px. The ground truth's equivalent focal is 336 px, so the vote is 20% low once the draw is averaged out, and the single seed-0 draw (309.6 px) had been reading high by chance. At `K = 5` the base file's vote is 283.0 px.

`SeoulBull`, pre-fix (from the section above) against post-fix. The base `h00` is seeded from frames 0 to 4, poses 8 frames and releases 338.9 px in both.

| class | run | vote f post (pre) | `h00` seed changed (pre → post) | `h00` pass (pre → post) | qualified in set (pre → post) | set differs from base (pre → post) | first stage that differs (pre → post) |
|---|---|---|---|---|---|---|---|
| base | repeat | 283.0 (309.6) | no → no | PASS → PASS | 2 → 0 | no → no | none → none |
| random drop | s1 | 247.0 (317.5) | no → no | PASS → PASS | 2 → 0 | yes → yes | `h01` focal scan → `h00` flags (`vote_divergence`) |
| random drop | s2 | 263.5 (280.7) | no → no | PASS → PASS | 0 → 0 | yes → yes | `h00` flags → `h00` focal scan |
| random drop | s3 | 263.5 (258.3) | yes → yes | PASS → PASS | 0 → 0 | yes → yes | `h00` probe → `h00` probe (reach) |
| drop movers | | 283.0 (313.2) | no → no | PASS → PASS | 2 → 0 | yes → yes | `h00` release → `h00` release (spline accepted) |
| drop non-movers | s1 | 269.8 (316.6) | no → no | PASS → PASS | 2 → 0 | yes → yes | `h01` focal scan → `h00` focal scan |
| drop non-movers | s2 | 264.0 (264.7) | yes → yes | **fail** (focal −9.4%) → **fail** (focal −5.1%) | 0 → 0 | yes → yes | `h00` probe → `h00` probe (reach) |
| drop non-movers | s3 | 283.0 (283.0) | no → no | PASS → PASS | 0 → 0 | yes → no | `h00` flags → none |
| swap movers | | 283.0 (309.6) | no → no | PASS → PASS | 2 → 0 | no → no | none → none |
| ungated | | 279.6 (298.6) | no → no | PASS → PASS | 1 → 0 | yes → yes | `h00` release → `h00` release (spline accepted) |
| reorder | s1, s2 | 283.0 (309.6) | no → no | PASS → PASS | 2 → 0 | yes → no (bit-identical) | `h02` release → none |

Every finite `SeoulBull` candidate post-fix: `h00` passes in 12 of 13 runs (all but drop non-movers s2), `h01` in 12 of 13 (all but drop non-movers s1, released at 310.3 px), and `h02`, a 3-frame flat-scan window, in none. No candidate qualifies in any run.

`KerryPark480`, pre-fix against post-fix. Post-fix the vote is 144.0 px in every run but the ungated file (160.3 px); pre-fix it read 176.3, 215.7, 166.3 and 216.1 px on random drops s1 and s2, drop non-movers s2 and the ungated file. The equidistant verdict focal is 138.3 px in every run; the base `h00` is the rotation core of the first pass, poses 28 frames and passes, pre and post.

| class | run | `h00` posed (pre → post) | `h00` pass (pre → post) | complement candidates passing (pre → post) | seed frames of `h01` to `h06` changed (pre → post) | first stage that differs (pre → post) |
|---|---|---|---|---|---|---|
| base | repeat | 28 → 28 | PASS → PASS | 0 → 0 | 0 → 0 | none → none |
| random drop | s1 | 28 → 28 | PASS → PASS | 1 → 1 | 5 → 5 | `h01` widen → `h01` widen |
| random drop | s2 | 28 → 28 | PASS → PASS | 0 → 0 | 4 → 4 | `h00` probe (gauge pair) → `h00` probe (gauge pair) |
| random drop | s3 | 29 → 29 | PASS → PASS | 1 → 2 | 6 → 6 | `h00` probe → `h00` probe (core of 23) |
| drop movers | | 23 → 23 | PASS → **fail** (centre 11.9%) | 1 → 2 | 6 → 6 | `h00` probe (locked resection) → same |
| drop non-movers | s1 | 15 → 15 | PASS → PASS | 2 → 2 | 6 → 6 | `h00` probe (locked resection) → same |
| drop non-movers | s2 | 28 → 28 | PASS → PASS | 0 → 0 | 4 → 4 | `h01` widen → `h01` widen |
| drop non-movers | s3 | 19 → 19 | **fail** (centre 10.4%) → **fail** (centre 10.6%, rotation 1.81°) | 3 → 3 | 6 → 5 | `h00` probe (locked resection) → same |
| swap movers | | 28 → 28 | PASS → PASS | 0 → 1 | 6 → 6 | `h00` probe (gauge pair) → `h00` widen |
| ungated | | 26 → 26 | PASS → PASS | 1 → 1 | 6 → 6 | `h00` probe → `h00` probe (core of 20) |
| reorder | s1 | 28 → 28 | PASS → PASS | 0 → 0 | 2 → 0 | `h03` widen → none (bit-identical) |
| reorder | s2 | 28 → 28 | PASS → PASS | 1 → 0 | 4 → 0 | `h02` widen → none (bit-identical) |

The seed-frame counts here compare the `h01` to `h06` seed groups position by position with the base run's, the same reading for both columns. Qualified candidates in the base set go from 7 to 6.

**What changed and what did not.**

1. **Member order.** Both reorder files of both captures give products bit-identical to the base file's. Before the change, every one of the four differed from the base in at least one candidate.
2. **The vote.** The `SeoulBull` vote now moves between 247 and 283 px across the perturbed files, against 258 to 318 px before, and its base reading is 283.0 px rather than 309.6 px. Five draws narrow the spread over seed sets by about a third, not to nothing: on random drop s1 the vote reads 247.0 px, 12.7% under the base file's, which is inside the spread five draws still have on one file (245 to 298 px).
3. **The gauge pair.** With the parallax read over five draws, pair (0, 1) on `KerryPark480`'s swapped-mover file keeps its 1.48° reading and stays the gauge; its rotation core is the base run's, and the run first differs at `h00`'s widen, which admits frame 23 where the base admits frame 39. Random drop s2 still moves the gauge to (3, 4): there pair (0, 1) reads 0.69° under every seed, so that change follows the evidence.
4. **The commit bar.** On `SeoulBull`, random drop s3 and drop non-movers s2 still seed `h00` from frames 12 to 16. Group 0 to 4 grows to 7 frames as before and now clears the kept count, but at 7 frames its capture-level reach is 59%, under the bar's 60%, so it still fails the bar. The reach floor is the next hard threshold the baseline sits next to.
5. **`edge_scan`.** At 283 px the vote is 16% under the ground truth, so the scan grid centred on it ends at about 374 px and the inlier fraction rises to its top point on every `SeoulBull` file: `edge_scan` fires on every `h00` and no candidate qualifies, although every `h00` but drop non-movers s2's passes against the ground truth. The single draw had hidden this by reading 309.6 px on the base file; it is the reading the lattice scan and its extension ([seed-hypothesis-loop.md](seed-hypothesis-loop.md#focal-scan), measured in the [next section](#pick-stability-after-milestone-b-2026-10-09)) address.
6. **The locked resection.** On `KerryPark480`, drop movers and drop non-movers s1 and s3 lose the same resections as before (images 31, 6 and 28, with 31, 9 and 12 survivors at the last round). Drop movers' `h00` again poses 23 frames, with frame 39 where the pre-fix run posed frame 37, and now fails against the ground truth with an 11.9% centre error where the pre-fix run passed. Its core lost image 31 and the frames resected after it in both runs, and it is a run the converged-set resection ([rotation-locked-resection.md](rotation-locked-resection.md#mechanism), measured in the [next section](#pick-stability-after-milestone-b-2026-10-09)) addresses.

**What this decides.**

- **Reorder acceptance is met.** Member order no longer reaches any decision on either capture.
- **Drop acceptance is not met yet.** The first candidate still changes on 2 of 6 drop files on `SeoulBull` (as before) and 4 of 6 on `KerryPark480` (as before). On the 22 perturbed files it changes on 9, against 9 before, and the first candidate fails the ground truth on 3, against 2 before. Every remaining change of the first candidate enters at a stage the draft's unbuilt changes cover: the probe through the reach floor (`SeoulBull`), the scan grid and `edge_scan` read off a vote 16% low (`SeoulBull`), the rotation-locked resection's survivor floor (`KerryPark480`), and the gauge pair where every draw agrees (`KerryPark480` random drop s2, a change in the evidence).
- **Neither base `h00` regresses.** Both base files' `h00` pass against their ground truth with the same seed frames and posed frames as before; `SeoulBull`'s loses its qualification to `edge_scan`.
- **`K = 5` stands for now.** The vote's remaining log s.d. at five draws, 0.063, is under half of a scan step (ln 1.15 = 0.14); more draws narrow it slowly (a log s.d. of 0.023 at 33 draws), and the lattice grid is what keeps a residual movement from moving a grid point. The seed's run time did not measurably change: at one draw and at five draws, and with the pre-change scripts, a `KerryPark480` base run takes 53 to 67 s on the same machine on the same day.

## Pick stability after milestone B (2026-10-09)

**Question.** The [section above](#pick-stability-after-the-deterministic-fixes-2026-10-09) left three causes of first-candidate changes: a focal scan whose grid, `edge_scan` verdict and refit band all move with the vote; a seed group chosen as the first one tried among equal outcomes; and a rotation-locked resection that drops an image when any trim round dips under its survivor floor. With changes for those built, which runs of the [perturbation suite](#pick-stability-under-small-changes-to-the-cluster-file-2026-10-09) still change the first candidate, does either capture's base `h00` regress, and does `SeoulBull`'s base set qualify again?

**Method.**

- **Code.** Branch `bootstrap-core-migration` with three changes on top of the section above, described in [seed-hypothesis-loop.md](seed-hypothesis-loop.md) § "Focal scan" and § "Choice among seed groups" and in [rotation-locked-resection.md](rotation-locked-resection.md):
  - the focal scan reads rungs of one lattice per capture (`max(w, h) * 1.15^k`), the 5 nearest the structure-free focal; it extends past a peak at either end by up to 3 rungs, takes the best refit among the rungs within 5 points of the best rung, lets the vote break only a tie within half a point, and keeps a release step that ties the winner within half a point;
  - an attempt probes every seed group before gating any, against the best probe of all of them, finishes every measurable group, and chooses among the finished outcomes by reach, then kept frames, then median residual, then sorted image names;
  - the resection's survivor floor is read on its final kept set, and the fisheye rotation core retries a failed image once with a 16 px trim gate, judged at 8 px.
- **Files, seed runs and scoring.** The same 13 files per capture, seed command, environment, instrumentation, ground-truth scorer and pass bar as the two sections above.
- **Columns.** *pre-A* is [the first section](#pick-stability-under-small-changes-to-the-cluster-file-2026-10-09), *post-A* [the second](#pick-stability-after-the-deterministic-fixes-2026-10-09), *post-B* this one. The first candidate is *changed* when `h00`'s seed frames or posed frames differ from the base file's in the same column. The *first diverging stage* comes from comparing every committed finite candidate with the base run's in commit order and naming the first difference: the probe (seed frames), the posed set (the rotation core or the widen), the commit bar (kept frames, or reach or spread on the other side of its bar), the focal scan (winners more than 4% apart), the release (spline verdict, or focals more than 1% apart), or the flags and qualification. The same automatic comparison fills all three columns, so a pre-A or post-A entry can name an earlier stage than the hand attribution of the sections above. A run *differs from the base* when this comparison finds a difference.

**Result: every `h00` on both captures passes its ground truth, `SeoulBull`'s base set qualifies again, and the first candidate still changes on the same `SeoulBull` files and on one more `KerryPark480` file.**

`SeoulBull`. The vote reads 283.0 px on the base file in both post columns. Every file but random drop s1 scans the same five rungs, 207.5 to 362.9 px; its scan peaks at the top rung and is extended by one rung to 417.4 px, where the inlier fraction falls. Random drop s1, with a vote of 247.0 px, scans the window one rung lower, 180.4 to 315.6 px, and is extended two rungs up to the same 417.4 px. No `h00` carries `edge_scan`. On the base file the refits of the rungs at 315.6, 362.9 and 417.4 px read 81.3%, 83.7% and 83.4%, 362.9 px wins, and the release walks from it to 338.9 px at the same inlier fraction.

| class | run | first candidate changed (pre-A / post-A / post-B) | `h00` against the ground truth (pre-A / post-A / post-B) | qualified in set (pre-A / post-A / post-B) | first diverging stage (pre-A → post-A → post-B) |
|---|---|---|---|---|---|
| base | repeat | no / no / no | PASS / PASS / PASS | 2 / 0 / 2 | none → none → none |
| random drop | s1 | no / no / no | PASS / PASS / PASS | 2 / 0 / 0 | `h01` focal scan → `h00` flags → `h00` flags (`vote_divergence`) |
| random drop | s2 | no / no / no | PASS / PASS / PASS | 0 / 0 / 2 | `h00` focal scan → `h00` focal scan → none |
| random drop | s3 | yes / yes / yes | PASS / PASS / PASS | 0 / 0 / 1 | `h00` probe → `h00` probe → `h00` probe |
| drop movers | | no / no / no | PASS / PASS / PASS | 2 / 0 / 2 | `h00` release → `h00` release → `h00` release (spline accepted) |
| drop non-movers | s1 | no / no / no | PASS / PASS / PASS | 2 / 0 / 2 | `h01` focal scan → `h00` focal scan → none |
| drop non-movers | s2 | yes / yes / yes | **fail** (focal −9.4%) / **fail** (focal −5.1%) / PASS | 0 / 0 / 1 | `h00` probe → `h00` probe → `h00` probe |
| drop non-movers | s3 | no / no / no | PASS / PASS / PASS | 0 / 0 / 2 | `h00` focal scan → none → none |
| swap movers | | no / no / no | PASS / PASS / PASS | 2 / 0 / 2 | none → none → none |
| ungated | | no / no / no | PASS / PASS / PASS | 1 / 0 / 2 | `h00` release → `h00` release → `h00` release (spline accepted) |
| reorder | s1, s2 | no / no / no | PASS / PASS / PASS | 2 / 0 / 2 | `h02` release → none → none (bit-identical) |

Every finite `SeoulBull` candidate post-B: `h00` passes in 13 of 13 runs and `h01` in 13 of 13, released between 338 and 352 px; `h02`, a 3-frame flat-scan window, passes in none. Random drop s1's `h00` releases 338.9 px and passes, but its vote reads 247.0 px and the release is 37% above it, outside the 35% divergence band, so it does not qualify.

`KerryPark480`. The vote is 144.0 px and the equidistant verdict 138.3 px in every run but the ungated file (vote 160.3 px). Every run scans 103.2 to 180.4 px, every `h00` scan peaks at the middle rung (136.4 px), and every `h00` releases between 135.9 and 136.5 px.

| class | run | first candidate changed (pre-A / post-A / post-B) | `h00` posed (pre-A / post-A / post-B) | `h00` against the ground truth (pre-A / post-A / post-B) | `h00` median per-frame centre error, % (pre-A / post-A / post-B) | `h00` maximum per-frame centre error, % (pre-A / post-A / post-B) | complement candidates passing (pre-A / post-A / post-B) | first diverging stage (pre-A → post-A → post-B) |
|---|---|---|---|---|---|---|---|---|
| base | repeat | no / no / no | 28 / 28 / 28 | PASS / PASS / PASS | 0.67 / 0.67 / 0.60 | 1.5 / 1.5 / 2.1 | 0 / 0 / 1 | none → none → none |
| random drop | s1 | no / no / yes | 28 / 28 / 30 | PASS / PASS / PASS | 1.40 / 1.36 / 1.26 | 36.4 / 36.4 / 36.2 ¹ | 1 / 1 / 2 | `h01` posed set → `h01` posed set → `h00` posed set |
| random drop | s2 | yes / yes / yes | 28 / 28 / 29 | PASS / PASS / PASS | 4.22 / 4.22 / 1.65 | 32.2 / 32.2 / 35.8 ¹ | 0 / 0 / 0 | `h00` posed set in all three |
| random drop | s3 | yes / yes / yes | 29 / 29 / 29 | PASS / PASS / PASS | 1.65 / 1.69 / 1.95 | 35.7 / 35.7 / 36.4 ¹ | 1 / 2 / 0 | `h00` posed set in all three |
| drop movers |  | yes / yes / yes | 23 / 23 / 28 | PASS / **fail** (centre 11.9%) / PASS | 3.38 / 11.93 / 0.79 | 38.2 / 36.4 / 2.7 | 1 / 2 / 0 | `h00` posed set in all three |
| drop non-movers | s1 | yes / yes / yes | 15 / 15 / 30 | PASS / PASS / PASS | 4.54 / 4.59 / 0.40 | 27.7 / 27.6 / 1.2 | 2 / 2 / 1 | `h00` posed set in all three |
| drop non-movers | s2 | no / no / no | 28 / 28 / 28 | PASS / PASS / PASS | 1.89 / 1.58 / 1.68 | 35.8 / 35.8 / 35.2 ¹ | 0 / 0 / 0 | `h01` posed set in all three |
| drop non-movers | s3 | yes / yes / yes | 19 / 19 / 26 | **fail** (centre 10.4%) / **fail** (centre 10.6%) / PASS | 10.37 / 10.55 / 4.25 | 35.7 / 35.6 / 38.7 ¹ | 3 / 3 / 1 | `h00` posed set in all three |
| swap movers |  | yes / yes / yes | 28 / 28 / 28 | PASS / PASS / PASS | 2.37 / 0.88 / 0.57 | 34.9 / 10.9 / 1.5 | 0 / 1 / 0 | `h00` posed set in all three |
| ungated |  | yes / yes / yes | 26 / 26 / 26 | PASS / PASS / PASS | 0.50 / 0.56 / 0.39 | 1.6 / 2.1 / 12.8 | 1 / 1 / 2 | `h00` posed set in all three |
| reorder | s1 | no / no / no | 28 / 28 / 28 | PASS / PASS / PASS | 0.67 / 0.67 / 0.60 | 1.5 / 1.5 / 2.1 | 0 / 0 / 1 | `h03` posed set → none → none (bit-identical) |
| reorder | s2 | no / no / no | 28 / 28 / 28 | PASS / PASS / PASS | 0.67 / 0.67 / 0.60 | 1.5 / 1.5 / 2.1 | 1 / 0 / 1 | `h02` probe → none → none (bit-identical) |

The two centre-error columns score each of `h00`'s posed frames against the ground truth after a similarity alignment of its camera centres to the ground truth's, as a percentage of the capture's extent; the median column is the scorer's centre error, which the pass verdict reads. ¹ A pass on the median can hide single misplaced frames. Drop non-movers s3 passes post-B at a 4.25% median, 0.75 points under the 5% bar, with `fisheye_right/frame_16` at 38.7% and `fisheye_left/frame_07` at 32.6%. Random drop s1 carries `fisheye_right/frame_16` at 36% in all three columns, and random drops s2 and s3 and drop non-movers s2 carry the same frame at 32% to 36% in all three; the ungated file's post-B `h00` has one frame at 12.8%, and on the base file no frame is over 2.1%.

Qualified candidates in the base set go from 6 post-A to 8 post-B, and the base set now holds a passing complement candidate, `h03`.

**What changed and what did not.**

1. **The scan.** On `SeoulBull` the vote's movement over the perturbed files (247.0 to 283.0 px) moves the scanned window on one file by one rung and the scan winner on none; `edge_scan` no longer fires, and the base set's two candidates qualify again. Two rules carry this beside the lattice. A peak at the window's top rung is extended past rather than flagged. And the release keeps a step that ties the winner within half a point: the scan's top is flat from 315.6 to 417.4 px, the winner is a rung, and with the earlier rule (keep a step only when it raises the inlier fraction) the release stayed at the 362.9 px rung, 8% over the ground truth, on the three files where the spline release was refused, all three of which then failed. Choosing the rung nearest the vote among all refits within 5 points instead put `h01`'s release at 315.6 px, 3 points under its best rung and 6% under the ground truth, which is why the vote breaks only a half-point tie.
2. **The locked resection.** On `KerryPark480` the images the post-A runs lost now resect: image 31 with 31 survivors on drop movers, read on its final round; image 28 with 12 on drop non-movers s3; and image 6 on drop non-movers s1 with 17 survivors within 8 px after the 16 px retry. Drop non-movers s1 also resects image 8 through a second retry, with 21 survivors within 8 px. The rotation cores of those three files pose 22, 24 and 20 frames, against 17, 9 and 13 post-A and 22 on the base file; their `h00` poses 28, 30 and 26 frames, and all three pass. The base file's core still drops images 20 and 37, and the retry recovers neither.
3. **Choice among seed groups.** On `SeoulBull`, groups 0 to 4 and 12 to 16 both clear the commit bar in the first pass of every file whose `h00` is seeded from 0 to 4, and group 0 to 4 ranks first on reach (71% against 65% on the base file). On `KerryPark480` the first pass is the rotation core, and the choice acts inside the complement passes, where outcomes tied at the 14-frame ceiling are now ordered by reach, kept frames and median residual instead of by try order.
4. **What still changes the first candidate.** Two of the changes are threshold edges: a reading one unit from a hard bar decides which side of it the run falls on.
   - On `SeoulBull`, random drop s3 and drop non-movers s2 still seed `h00` from frames 12 to 16, because group 0 to 4's reach there is 59% at 7 kept frames, one point under the commit bar's 60%. Their `h00` now passes.
   - On `KerryPark480` random drop s1, the rotation core's retry accepts image 20 with exactly 10 of its 62 observations within 8 px, at the floor of 10. The base file's core drops image 20, and one survivor fewer would drop it here too. This run's `h00` poses 30 frames against the base file's 28.

   The other `KerryPark480` changes are not at an edge. The first candidate's posed set differs from the base file's on every drop file and on the swapped-mover and ungated files: the rotation core's resections succeed on a different set of skeleton images (drop non-movers s1 also resects image 20, with 14 survivors at the ordinary gate), and the widen admits a different frame or two after it. Every one of those `h00` passes against the ground truth on its median centre error; footnote ¹ under the `KerryPark480` table names the single frames that do not.

**What this decides.**

- **Neither base `h00` regresses.** `SeoulBull`'s base `h00` is seeded from frames 0 to 4, poses 8 frames, releases 338.9 px and passes, with no `edge_scan`; `KerryPark480`'s is the rotation core's, poses 28 frames, releases 136.1 px and passes.
- **No perturbed file's first candidate fails its ground truth.** Over the 22 perturbed files `h00` fails on none, against 2 pre-A and 3 post-A.
- **The first candidate still changes on 10 of the 22 files** (2 on `SeoulBull`, 8 on `KerryPark480`), against 9 pre-A and 9 post-A, and 14 of the 22 runs differ from their base, against 21 pre-A and 16 post-A. Dropping 0.3% of the members still changes which skeleton images resect and which frames the widen admits on `KerryPark480`, and still moves one `SeoulBull` group across the 60% reach floor. None of these changes follows a RANSAC draw, a member or try order, or a vote-centred grid. Two of them are threshold edges (item 4): `SeoulBull`'s reach at 59% against the 60% bar, and `KerryPark480` random drop s1's retry at exactly the floor of 10 survivors. No further change is made for them: the reach floor, the resection floor and the widen gate read the evidence, and in this suite the first candidate's ground-truth verdict, read on the median centre error, does not depend on them.
- **The cost.** Finishing every seed group and the scan's extension rungs add no measurable wall time on these captures: a `SeoulBull` run takes 8 to 11 s and a `KerryPark480` run 44 to 64 s, against 10 to 20 s and 51 to 110 s for the same files post-A, on the same machine with other work running beside both. A copy of `fleetws` seeds in 23.7 s and a copy of `MurdoSmallAntiqueCat` in 36.6 s, against 23.2 s and 35.9 s for the same workspaces in the 2026-10-07 fleet run.

**The base files on the resection's extra solve.** After the runs above, the resection kernel gained two rules ([rotation-locked-resection.md](rotation-locked-resection.md#mechanism)): a trim round that keeps fewer than two observations fails, and when three rounds end with the kept set still changing, the translation is solved once more over the final kept set. The two base files were re-run on that kernel; the 22 perturbed files were not.

- `SeoulBull`'s product is bit-identical to the post-B base run.
- On `KerryPark480`, image 31's resection keeps the same 31 of 71 observations but does not stabilise in three rounds, so its returned translation moves. Image 6 then fails at 8 px and resects through the retry with 15 survivors, and image 20, which the post-B base run drops, resects through the retry with exactly 10 of 62 observations within 8 px, at the floor of 10. The rotation core poses 23 frames against 22, and `h00` 29 against 28. `h00` releases 136.2 px against 136.1 px, with a median rotation error of 0.42° in both and a median centre error of 0.53% against 0.60%, and passes; its worst frame, `fisheye_right/frame_14`, is at 11.3% where the post-B run's worst was at 2.1%. Two complement candidates, `h05` and `h07`, now pose 2 frames and do not qualify, so the base set's qualified candidates go from 8 to 6.

## Pick against the ground-truth captures (2026-10-09)

**Question.** The three sections above measured the pick-stability changes on two captures. Across all eight captures with an approved ground truth, does the seed's pick (`ladder_first`) agree with the ground truth more often, less often or as often with those changes as without them? The draft [seed-photometric-candidate-score.md](../../drafts/seed-photometric-candidate-score.md) quotes 4 of 8 for 2026-10-07.

**Method.**

- **Captures.** The seven approved entries of `C:/DataSets/workspace-prep/approved-gts.tsv` (`MurdoSmallAntiqueCat`, `DnDTabletop`, `OmniHilltop`, `fleetws`, `KerryPark360`, `OmniTemple1`, `KerryPark480`) and `SeoulBull`. `KerryPark480` and `SeoulBull` are scored against the checked-in `kerry_park_ground_truth.sfmr` and `seoul_bull_sculpture_ground_truth.sfmr`, the others against their approved files.
- **Files.** For each capture, a workspace copy in scratch held `.sfm-workspace.json`, `rig_config.json` where the capture has one, the capture's `*-clusters.matches` file, and only the images that file names, with their `.sift` files; no fleet workspace was written. `sfm cluster-patches --piecewise` at its other defaults (patch size 12, the gate at the refined shape, at most 8 capped cells) wrote one file per capture, and both seed runs of a capture read that same file. Kept members: 620,044 (`MurdoSmallAntiqueCat`), 1,370,494 (`DnDTabletop`), 506,437 (`OmniHilltop`), 66,003 (`fleetws`), 39,567 (`KerryPark360`, 82 images), 155,936 (`OmniTemple1`), 15,609 (`KerryPark480`) and 2,820 (`SeoulBull`). The last two are the kept counts of the base files of the sections above, and the after runs reproduce those base runs: `SeoulBull`'s `h00` releases 338.9 px, and `KerryPark480`'s poses 29 frames at 0.42° and 0.53% with its worst frame at 11.3%, as in [the base files on the resection's extra solve](#pick-stability-after-milestone-b-2026-10-09).
- **Code.** *Before* is `scripts/exp_fast_seed.py` at `25aa647a`, the commit before the changes, run from a git worktree. *After* is the same script at `0de454aa`, which carries the changes of `2da0e7bb`, `2ea211d8` and `0de454aa`.
- **The extension.** Both script versions import the one installed `sfmtool._sfmtool`, built from `0de454aa`. Between the two commits the extension gains a `draws` argument on `focal_vote` and `estimate_intrinsics` whose default of 1 is the single-draw kernel bit for bit, and the before script passes none, so the vote of the before runs is the pre-change vote. The extension also changes `resect_translation`: its survivor floor is read on the final kept set, a round keeping fewer than two observations fails, and a trim that has not settled is solved once more over its final set ([rotation-locked-resection.md](rotation-locked-resection.md#mechanism)). The before runs get that kernel without the 16 px retry the after script adds. To measure this one confound, the extension was also built from `25aa647a` into a separate target directory and imported ahead of the installed one, and the before script was run a second time on it (*before, own extension*).
- **Seed.** `pixi run -e dev python <script> <workspace copy>`, with `SFMTOOL_RELAX=0` and `SFMTOOL_SEED_RUNG1=3000`, with `SFMTOOL_SEED_EVAL`, `SFMTOOL_PINHOLE_RAW_VOTE` and `SFMTOOL_VOTE_DRAWS` unset, and with `SFMTOOL_SEED_EVO_DUMP`, `SFMTOOL_STAGE_DUMPS` and `SFMTOOL_RESECT_TRACE` set. Three seed runs ran at a time.
- **Scoring.** Every finite candidate was scored as in [Seed stage on the ground-truth entries](../patch/cluster-patch-refinement-measurements.md#seed-stage-on-the-ground-truth-entries): rotation error after a best-fit rotation, centre error after a similarity alignment of the camera centres (`sfmtool._sfmtool.analysis.estimate_alignment`) as a percentage of the ground truth's camera extent, and the equivalent focal over the ground truth's observed radii. A candidate passes with a median rotation error under 1°, a median centre error under 5% and an equivalent-focal error under 5%. The pick is the release file of the `ladder_first` candidate. A pick that poses no image of the ground truth does not agree with it.
- **Machine.** Windows 11, Intel Core i9-14900HX (32 logical processors), 64 GB RAM.

**Result: the pick agrees with the ground truth on 8 of 8 captures after the changes, against 4 of 8 before them, and 5 of 8 with the before scripts on their own extension.**

| capture | pick (before → after) | pick against the ground truth (before → after) | qualified (before → after) | finite candidates passing, of 8 (before → after) | `h00` posed (before → after) | `h00` median rotation, ° (before → after) | `h00` median / max centre error, % (before → after) | `h00` focal error, % (before → after) | wall time, s (before → after) |
|---|---|---|---|---|---|---|---|---|---|
| `MurdoSmallAntiqueCat` | `h00` → `h00` | no ground-truth image → PASS | 0 → 0 | 4 → 3 | 14 → 14 | – → 0.23 | – → 0.11 / 0.23 | – → −1.2 | 62.6 → 59.6 |
| `DnDTabletop` | `h00` → `h00` | **fail** (focal +5.8%) → PASS | 0 → 0 | 4 → 7 | 14 → 14 | 0.42 → 0.24 | 0.12 / 0.23 → 0.17 / 0.47 | +5.8 → −1.3 | 108.1 → 100.0 |
| `OmniHilltop` | `h00` → `h00` | PASS → PASS | 0 → 0 | 6 → 5 | 14 → 14 | 0.32 → 0.33 | 1.16 / 6.67 → 1.16 / 6.30 | +0.9 → +0.8 | 83.8 → 112.2 |
| `fleetws` | `h03` → `h00` | PASS → PASS | 3 → 0 | 8 → 6 | 14 → 13 | 0.28 → 0.30 | 0.74 / 1.42 → 0.57 / 1.90 | +1.9 → +3.6 | 50.0 → 38.1 |
| `KerryPark360` | `h01` → `h01` | **fail** (rotation 20.6°, focal +247%) → PASS | 1 → 3 | 0 → 5 | 14 → 12 | 29.38 → 0.16 | 3.90 / 12.96 → 0.04 / 0.14 | +262.4 → +0.1 | 44.2 → 78.2 |
| `OmniTemple1` | `h00` → `h00` | PASS → PASS | 7 → 8 | 2 → 3 | 14 → 14 | 0.30 → 0.29 | 4.42 / 23.02 → 4.38 / 23.20 | −0.9 → −1.1 | 73.3 → 78.4 |
| `KerryPark480` | `h00` → `h00` | **fail** (rotation 1.15°) ¹ → PASS | 8 → 6 | 2 → 1 | 17 → 29 | 1.15 → 0.42 | 1.64 / 12.73 → 0.53 / 11.29 | −0.0 → +0.0 | 44.1 → 48.3 |
| `SeoulBull` | `h00` → `h00` | PASS → PASS | 2 → 2 | 1 → 2 | 8 → 8 | 0.96 → 0.96 | 0.19 / 0.37 → 0.19 / 0.37 | +0.8 → +0.8 | 8.7 → 10.7 |

*Qualified* counts the qualified candidates among all 17 committed (7 on `SeoulBull`), rotation-only ones included. *Finite candidates passing* counts release files; a `-spline-refused` variant is not counted. `fleetws`'s pick before is `h03`, which passes at 0.25°, 0.61% and +2.0%; `KerryPark360`'s pick `h01` passes after at 0.89°, 0.52% and −0.7%. Every other pick is `h00`.

¹ With the before scripts on their own extension, `KerryPark480`'s `h00` poses 28 frames and passes (0.43°, 0.67% median and 1.54% maximum centre error, −0.0% focal), with 7 qualified candidates and 1 passing; that run takes 46.8 s. On the other seven captures the two before runs commit the same 17 candidates with the same focals, posed counts, flags and qualification, and take 46.8, 101.0, 79.5, 47.5, 41.6, 63.2 and 8.3 s in table order. The confound therefore moves one capture: the before column's 4 of 8 is 5 of 8 on the extension the before scripts were written for.

**Which captures moved and why.**

1. **`MurdoSmallAntiqueCat`: choice among seed groups.** Every seed group of the first attempt starves the widen, so the attempt commits a deferred outcome. Before, the outcome kept is the first tried among those with the most kept frames (14): one posing frames 1051 to 1261, at a capture-level reach of 4%. Those frames lie in the 24-frame stretch from 961 to 1306 that the ground truth leaves out, so the pick cannot be scored and counts as not agreeing. After, the outcomes are ordered by reach first, and the attempt commits one at 9% reach whose `h00` passes. The vote reads 2659.4 and 2659.3 px, so it plays no part.
2. **`DnDTabletop`: choice among seed groups, then the scan.** Before, the single-draw vote reads 2996.4 px and the first-tried deferred outcome (reach 10%) is committed. Its refit band holds 3962.7 px at 90.9% and the vote's 2996.4 px at 89.4%, the vote's point wins, and the release leaves the scan basin and keeps 2996.4 px, 5.8% over the ground truth. After, the vote over five draws reads 2922.1 px and the outcome at 13% reach is committed, which poses other frames. Its lattice window from 2195.5 to 3840.0 px peaks at the 2903.6 px rung, the release walks to 2828.3 px, and the spline rung to 2773.9 px, 1.3% under the ground truth.
3. **`KerryPark360`: the vote's draws decide the camera model.** The camera-model arbitration runs when the pinhole vote's pool holds at most 9 votes, among other triggers ([estimate-intrinsics.md](estimate-intrinsics.md)). Before, one draw per pair puts 13 votes in the pool at a log-focal IQR of 0.015, the vote stands at 434.0 px, and the capture is seeded as a pinhole: every finite candidate has a median rotation error over 6° and a focal more than 200% over the ground truth's equivalent focal. After, a pair votes only when more than half of its five draws give a focal, the pool holds 4 votes, `thin_pool` escalates to the camera-model columns, and the arbitration routes the capture to the equidistant model at 276.6 px. On the cascade file of the 2026-10-07 fleet run the single draw gave 5 votes and the same fisheye routing, so the piecewise file moved the single-draw pool across the cut.
4. **`KerryPark480`: the confound, not the changes.** On the extension the before scripts were written for, the rotation core poses 22 frames and `h00` passes with 28 posed. On the current extension, image 31's trim does not settle and is solved once more over its final set, image 6, which the after script recovers through its 16 px retry, does not resect, the rotation core poses 11 frames, and `h00` poses 17 and fails on rotation. After, the core poses 23 frames and `h00` 29, as in the section above.
5. **The captures that did not move.** `OmniHilltop`, `OmniTemple1` and `SeoulBull` pick an `h00` that poses the same images before and after, and its median errors differ by at most 0.01° and 0.04 points. `fleetws` picks `h03` before, the first qualified candidate; after, no candidate qualifies, the rank falls back to `h00`, and both pass.

**What this decides.**

- **The changes help.** The pick agrees with the ground truth on 8 of 8 captures, against 4 of 8 before them on the same files and the same extension, and 5 of 8 with the before scripts on their own extension. No capture moves from a passing pick to a failing one. Three captures move to a passing pick through the changes themselves: `MurdoSmallAntiqueCat` through the choice among seed groups by reach, `DnDTabletop` through that choice and the lattice scan, and `KerryPark360` through the vote over five draws. The fourth, `KerryPark480`, also passes before once the before scripts run on their own extension.
- **The 2026-10-07 figure is not this baseline.** Scored the same way, the 2026-10-07 fleet run's picks pass on 5 of its 7 entries (`MurdoSmallAntiqueCat`'s pick poses images of the same stretch the ground truth leaves out, and `KerryPark480`'s fails at a 14.9% centre error), and `SeoulBull`'s pick passed on its cascade file that day, so that run scores 6 of 8. It read each workspace's cascade file with the code of that day. The before column here reads the piecewise files at the current defaults, and it differs from that run on `DnDTabletop`, `KerryPark360` and `KerryPark480`.
- **The margins.** `SeoulBull`'s pick passes at a median rotation error of 0.96°, 0.04° under the bar, before and after; `OmniTemple1`'s at a median centre error of 4.38%, with a worst frame at 23%; `KerryPark480`'s with a worst frame at 11.3%. More finite candidates pass after on four captures and fewer on four; over the eight captures, 27 pass before and 32 after.
- **The cost.** Summed over the eight captures, the seed takes 475 s before and 526 s after, with three runs at a time (435 s for the before scripts on their own extension). `KerryPark360` takes 78 s after against 44 s before because it now runs the fisheye seed, and `OmniHilltop` 112 s against 84 s. On the other six captures the after run is faster, or slower by at most 10%, except `SeoulBull`, at 10.7 s against 8.7 s.
