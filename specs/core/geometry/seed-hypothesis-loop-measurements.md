# Seed Hypothesis Loop Measurements

This file records measurements behind [seed-hypothesis-loop.md](seed-hypothesis-loop.md), the seed stage that develops a set of candidate reconstructions from a cluster-patches `.matches` file and records which candidate it ranks first. The measurements bear on how stable that set and that first candidate are when the input file changes slightly, which stage of the exploration turns a small change of the input into a different candidate, and which of those stages behave as defects rather than as a response to the evidence. The changes they lead to are proposed in [seed-pick-stability.md](../../drafts/seed-pick-stability.md).

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
- **What to change** is proposed in [seed-pick-stability.md](../../drafts/seed-pick-stability.md): a vote measured over several RANSAC draws, a scan grid and probe focal that do not move with every change of the vote, a commit bar below the cap, choice among outcomes by a continuous score with a tie-break on content, a locked resection that judges its converged set and retries a failed image, and member order fixed at load. The perturbation suite here is its acceptance test.
