# Seed Pick Stability

**Status:** Draft. Decided:
- the seed's candidate set and first candidate should change only when the evidence changes, not when a RANSAC seed, a member's position in its cluster, or a try order changes;
- the changes below are accepted when the perturbation suite of [seed-hypothesis-loop-measurements.md](../core/geometry/seed-hypothesis-loop-measurements.md#pick-stability-under-small-changes-to-the-cluster-file-2026-10-09) gives bit-identical results under member reordering, and unchanged first candidates (seed frames, posed frames, qualification, ground-truth verdict) under the random and non-mover drops on `SeoulBull` and `KerryPark480`.

Not decided: the number of RANSAC draws per pair vote; whether the probe focal snaps to the scan lattice or the probe runs at two focals; the score that orders deferred outcomes. See [Open questions](#open-questions).

Amends:
- [core/geometry/seed-hypothesis-loop.md](../core/geometry/seed-hypothesis-loop.md) § "Rank", and the exploration it records the outcome of (§ "Capture-level measurements", § "Ladder dedup")

## Purpose

The seed stage explores a capture's cluster evidence and commits a set of candidate reconstructions, ranked first-qualified in commit order. Perturbing the cluster file by under 0.4% of its kept members changes that set on both checked-in ground truths, and in 2 of 22 perturbed runs turns the first candidate from passing the ground truth to failing it. The rank itself reorders nothing. The change enters in the first candidate's own exploration, and the complement queue then carries it through every later pass.

The measurements trace most of it to four kinds of instability that have nothing to do with the evidence:

1. a single RANSAC draw with a fixed seed, read as if it were the estimate (the pairwise focal vote, the rotation core's gauge pair);
2. a hard threshold that the baseline sits within a few counts or percent of (the commit bar of 8 kept frames, equal to the core cap of 8; the locked resection's 10 survivors on any trim round; the 1.0° gauge parallax floor; `edge_scan` and the 5-point refit band on a grid centred on the vote);
3. a tie broken by the order things were tried (the deferred outcome with the most kept frames, first on a tie, where kept frames saturate at the cap plus the widen ladder);
4. member order inside a cluster (first-observation reference rays; observation-order sums in iteration-capped adjustments).

This draft removes each of them without changing what the exploration measures.

## Changes

### Member order fixed at load

The loader sorts each cluster's members by image index before any stage reads them, so the file's member order cannot reach a solver. Every reading that takes "the first observation" of a cluster instead takes the cluster's reference member, which the file names: `core_parallax` measures the widest angle over every pair of a point's posed views, and the ladder's far-field reading and the rotation-only layers take the reference member's ray. Measured by the reorder files: the outputs are bit-identical to the base file's.

### A vote measured over several draws

Each pair's Bougnoux and rotation votes are estimated over `K` RANSAC seeds and the pair contributes the median of its `K` estimates. On `SeoulBull` one draw per pair moves the pooled vote between 251 and 323 px across seeds 0 to 9 on one unchanged file, which is the whole range the dropped members produced; the median of `K` draws removes the seed from the reading without changing what is measured. The rotation core's gauge pair is chosen the same way: its parallax and cheiral count are the medians over `K` seeds, so pair (0, 1) on `KerryPark480`, at 1.48° under 11 of 12 seeds, keeps its place.

### Readings of the vote that do not move with it

The vote stays the referee for the divergence band. Three readings stop using its exact value:

- **Scan grid on a fixed lattice.** The grid points are `f_ref · 1.15^k` for integer `k`, with `f_ref` fixed by the camera's field-of-view band, and the scan takes the five lattice points nearest the vote. A vote that moves by less than half a step moves no grid point.
- **Grid extension instead of `edge_scan`.** When the inlier fraction peaks at the grid's edge, the scan adds lattice points past that edge until it falls or the field-of-view band ends. `edge_scan` is set only when the band ends first, which is the upward affine escape the flag exists for.
- **Scan winner by refit, not by distance to the vote.** The two refits are compared on inlier fraction. The vote breaks a tie only when one candidate lies outside the vote's own precision band (`f_band`) and the other inside it.

The probe focal is the lattice point nearest the vote, so it shares the grid's step. See [Open questions](#open-questions).

### Commit bar below the cap

The commit bar's kept count is read against the core cap: an outcome commits on at least `min(8, cap − 1)` kept frames, so one frame lost in growth at the cap does not fail the bar. On `SeoulBull` the bar and the cap are both 8, and a probe focal of 264.7 px instead of 258.3 px grows group 0 to 4 to 7 frames and moves the first candidate to another group.

### Choice among outcomes by a continuous score

An attempt finishes every seed group it tries (it already stops at 8 groups, and the probe and widen memos carry repeats) and returns the outcome with the best `score` (spread over the observability bar, then coverage), with the probe's inlier fraction after that, and a tie broken on the seed group's image names, sorted. The deferred-outcome slots (`best`, `flat`) use the same comparator. The probe's measurability gate is read against the maximum consensus over every group's probe, not over the groups probed so far. On `KerryPark480` two or three groups tie at the 14-frame ceiling in five of the base run's passes, and the first one tried is committed.

### Locked resection judged on its converged set

`resect_translation`'s survivor floor applies to the converged set, not to every trim round. An image whose resection fails is kept in the skeleton and tried again after the next accepted image, and is dropped only when a full round over the remaining skeleton accepts nothing. On `KerryPark480`, image 6 keeps 13 survivors in the base file and 9 with 56 non-movers dropped, and its removal costs the core 13 of its 22 frames.

## Open questions

- **`K`.** The vote's cost is linear in it. The question is the smallest `K` at which the `SeoulBull` vote's spread over seed sets is under its pool's own interquartile range.
- **Probe focal.** Snapping to the lattice removes small vote movements but still flips at a step's midpoint. Probing at the two lattice points that bracket the vote and keeping the better outcome removes the flip and doubles the probe's cost.
- **The score among outcomes.** Spread then coverage is the ladder's existing comparator; whether the photometric candidate score of [seed-photometric-candidate-score.md](seed-photometric-candidate-score.md) should order outcomes inside an attempt too is open.
- **Iteration caps.** With member order fixed, the capped adjustments are deterministic on one file but still amplify small changes of the input near a gate (the widen ladder's `0.35 · med_inl`). Whether the widen gate needs a margin is measured after the changes above, by the drop files.
