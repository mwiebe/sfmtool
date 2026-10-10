# Seed Photometric Candidate Score

**Status:** Draft. Decided:
- every candidate the seed's hypothesis loop releases is scored photometrically, at the coarse patch tier, by resampling each posed view into the candidate's own patches and reading how well the views agree;
- the score enters the qualification rule and the rank beside the geometric signals they read today, and a candidate whose photometric score is poor cannot qualify however good its reach and focal spread;
- a release carries patch frames, built from the cluster file's per-cell displacements and the candidate's poses, so the score reads the same frames a downstream consumer would;
- the score is validated against the eight ground-truth-backed captures by how well it orders the candidates that agree with the ground truth above those that do not, within each capture, and by whether it orders the candidates of the perturbed cluster files alike; the number of captures on which the seed's pick agrees with the ground truth, 8 of 8 without the score, is the check that no pick regresses.

Not decided: the aggregation from per-patch agreement to one candidate score; how the score is weighed against reach in the rank; whether rotation-only members are scored the same way. See [Open questions](#open-questions). [Measured](#measured) records what the score as first built separates: no aggregation of it orders the candidates better than the count of points it scores.

Amends:
- [core/geometry/seed-hypothesis-loop.md](../core/geometry/seed-hypothesis-loop.md) § "Rank" and § "Product"
- [core/geometry/seed-candidate-evaluation.md](../core/geometry/seed-candidate-evaluation.md): a new channel beside "Exact photometric witness", which today serves rotation-only members only
- [core/geometry/seed-drive.md](../core/geometry/seed-drive.md): the rank the drive reads

Related drafts: [cell-plane-normal-precision-gate.md](cell-plane-normal-precision-gate.md) builds the frames from the cell displacements the piecewise refinement of [cluster-patch-refinement.md](../core/patch/cluster-patch-refinement.md#piecewise-refinement) stores and the normals [cell-plane-normals.md](../core/patch/cell-plane-normals.md) derives from them; [piece-gated-grid-normal.md](piece-gated-grid-normal.md) the normals in them; [two-tier-patch-density.md](two-tier-patch-density.md) the tier the score runs at and the handoff of the pick to the fine tier; [sharper-patch-bitmap.md](sharper-patch-bitmap.md) how a view is scored against a bitmap.

## Purpose

The seed stage takes a set of images and their clustered feature matches and produces a handful of small candidate reconstructions, each a few camera poses and the points they see. Downstream, one of those candidates is grown into the full reconstruction, so choosing the right one matters more than producing many.

Today the choice reads geometry only. A candidate qualifies when it poses enough images, when its focal scan showed a clear optimum, and when none of four flags fired; qualification does not read reach, the fraction of the capture a candidate can explain. The rank is the order of qualification. None of those signals looks at a photograph.

Scored against the eight captures that have an approved ground truth, the seed's pick agrees with the truth on all eight ([seed-hypothesis-loop-measurements.md](../core/geometry/seed-hypothesis-loop-measurements.md#reach-and-retry-floors-without-a-fixed-count-2026-10-10)). The released sets still hold candidates that do not agree, 31 of the 65 finite release files that pose at least three ground-truth images, and the pick rests on the order of qualification and on thresholds that the [perturbation suite](../core/geometry/seed-hypothesis-loop-measurements.md#pick-stability-under-small-changes-to-the-cluster-file-2026-10-09) showed the candidates lie close to. Several of the failing candidates share one shape: rotations within a degree and focal within a percent, with the camera centres collapsed to a quarter or a half of their true spread, or one camera flipped. Reprojection residuals do not see this, because a collapsed layout re-triangulates its points to fit. Resampling the photographs into the candidate's patches does see it, because a wrong layout places the same patch at the wrong depth and the views disagree about what it looks like.

This draft adds that reading to the seed's choice.

## The score

For each candidate, for each finite point with a patch frame, every posed view that sees it is rendered into the patch's grid at the coarse tier, and the member coherence of the rendered tiles is read: the median pairwise ZNCC, or the median against the reference view's tile, per [sharper-patch-bitmap.md](sharper-patch-bitmap.md). A point's agreement is one number. A candidate's score aggregates those over its points, weighted so that a candidate is not rewarded for posing few images or keeping few points.

A candidate with a collapsed layout places each patch at a depth that fits the reprojections but not the parallax, so views from different positions render the surface at the wrong scale and offset relative to one another, and the agreement falls. A candidate with a flipped camera renders that camera's views mirrored or from behind, and their agreement with the rest is near zero. A correct candidate at the right focal renders every view consistently.

### Where it enters

**Qualification.** A candidate qualifies only if its photometric score passes a gate derived from the scores of all candidates in the release, so a release where every candidate is wrong does not qualify any. This is the same capture-relative gating principle the evaluation channels already use.

**Rank.** The rank orders qualified candidates by score, with reach as the tiebreak, instead of by order of qualification alone.

### What the frames are

A release today carries no patch frames; the writer's surfel block runs only for debug snapshots. For the score to read frames, the release must build them. The source is the cluster-patches file the seed already reads: each kept member's affine shape and per-cell displacements, back-projected through the candidate's poses, give the patch's placement and, via the gated plane fit, its normal. Nothing is rendered to build a frame; rendering happens only when the score reads the views.

## Theory

### Why photometry sees what geometry does not

A seed candidate is small, and small reconstructions are weakly constrained: a family of layouts fits the same reprojections to within the noise, differing in the twist of the camera path and the scale of its translations. Reprojection error and reach are blind within that family. The appearance of a surface is not, because the mapping from one view to another through a patch depends on the patch's depth and normal, and those depend on where the cameras are. A collapsed layout can only fit by placing surfaces at depths that make the inter-view mappings wrong. That is exactly what the project's photometric vetting elsewhere already relies on, and the seed is the one stage that has not used it.

### Why the coarse tier

Up to eight candidates per capture, each with hundreds of points and a dozen views, is thousands of renders. At `R = 12` that runs in seconds on the fleet. At `R = 24` it is four times that, for a decision that only needs to separate right from wrong, not to localise anything. The pick is then the one candidate handed to the fine tier.

### Validation

A candidate's ground-truth verdict is the existing pass rule: median rotation error under 1°, median camera-centre error under 5% of the ground truth's extent after similarity alignment, focal within 5%. The eight captures are the seven approved fleet entries plus `seoul_bull_sculpture`, with `kerry_park` scored against the checked-in ground truth; `C:/DataSets/workspace-prep/pick_suite.py` seeds them and gives every finite candidate its verdict.

The seed's pick agrees with the truth on 8 of 8 captures without the score, so the pick count cannot show the score helping; it is the check that the score makes no pick worse. The score itself is measured by:

- **Separation.** Within each capture, a pair of a passing and a failing candidate is ordered correctly when the passing one scores higher. The fraction of all such pairs over the eight captures that are ordered correctly is the score's separation. Beside it is read, per capture, whether the highest-scoring candidate with a verdict passes.
- **Stability.** The 13 cluster files per capture of the perturbation suite on `SeoulBull` and `KerryPark480` change the seed's first candidate on 2 and 9 of them. The score should order the candidates of each file alike, and its highest-scoring candidate should pass on the files where the geometric first candidate changes.

### Measured

**Method.** [`scripts/seed_photometric_score.py`](../../scripts/seed_photometric_score.py) scores every finite release file of a seeded workspace; nothing it computes feeds the seed. A release carries no patch frame, so the script builds one per finite point: the centre is the candidate's point; the normal is the [cell plane normal](../core/patch/cell-plane-normals.md) of the point's cluster in the workspace's cluster-patches file, read through the candidate's own cameras and poses (`sfmtool.analysis.cell_plane_normals`), where both of its axes are fixed, and the mean viewing direction elsewhere; the extent is the cluster refinement's `refine_radius` times each observation's member-shape scale, back-projected, median over the views (`PatchCloud.from_tracks` with `extent="feature_size"`). Only the normal reads the cell displacements; a frame placed from the triangulated cells is not built. Every tile is rendered at `R = 12` through that frame, anchored on the observation's stored keypoint (the candidate's inline `keypoints_xy`, the cluster members' refined positions), never at a reprojection. A point's readings are the median, over its observations other than the reference, of the blur-matched and of the plain ZNCC against the bitmap of the reference view the [reference-view rule](../core/patch/reference-view.md) picks (`PatchCloud.render_bitmaps`, `score_against_bitmap`), and the median off-diagonal entry of the plain and of the blur-matched member-coherence matrix (`PatchCloud.validate_member_coherence` with `keypoint_anchor=True`). A candidate's score is the median over its points, the mean over its points without the highest and lowest 10%, or the share of its points above the median of that reading pooled over every candidate of the capture. The candidates are those of the eight-capture pick-suite run at `de51f807` and those of the 26 perturbation-suite runs under the qualification rule of [2026-10-10](../core/geometry/seed-hypothesis-loop-measurements.md#reach-and-retry-floors-without-a-fixed-count-2026-10-10), each with the verdicts of its own run. A point lying on an observing camera's centre cannot be framed and is left out: 118 points on the eight captures, in candidates of `KerryPark360`, `OmniTemple1` and `KerryPark480`.

**Result: the score orders passing above failing candidates no better than the number of points it scores.** On the eight captures there are 109 pairs of a passing and a failing candidate:

| reading | median over points | trimmed mean | share above the capture's median |
|---|---|---|---|
| against the reference bitmap, blur-matched | 0.716 (6 of 8) | 0.716 (5 of 8) | 0.716 (6 of 8) |
| against the reference bitmap, plain | 0.716 (6 of 8) | 0.725 (6 of 8) | 0.716 (6 of 8) |
| median pairwise, plain | 0.697 (3 of 8) | 0.716 (6 of 8) | 0.734 (7 of 8) |
| median pairwise, blur-matched | 0.697 (3 of 8) | 0.725 (6 of 8) | 0.734 (6 of 8) |

Each cell is the fraction of pairs ordered correctly and, in brackets, the number of captures whose highest-scoring candidate with a verdict passes. Ordering the candidates by the number of points scored, with no photograph read, gives 0.734 (7 of 8); by posed count, 0.683, with 55 of the pairs tied.

- **The readings lie close together.** A candidate's median agreement against the reference bitmap lies between 0.90 and 0.98, and within a capture the candidates with a verdict span 0.011 (`MurdoSmallAntiqueCat`, `DnDTabletop`) to 0.033 (`SeoulBull`). The candidates with a focal error over 25% score lowest on `SeoulBull` (`h02`, focal −53%) and `DnDTabletop` (`h02`, +39%), but on `MurdoSmallAntiqueCat` `h06` (focal +41%, rotation 6.1°) and `h07` (focal −27%) score above the passing `h01` and `h02`.
- **Collapsed layouts are not seen.** The worst-ordered pair is on `KerryPark480`: the pick `h00` (29 posed, rotation 0.42°, centre 0.53%) reads 0.945 against the bitmap and `h04` (10 posed, rotation 1.16°, centre 12.1%) reads 0.959. The 14-frame candidates of that capture at 10.0 to 19.9% centre error read 0.929 to 0.941, below `h00`. On `OmniHilltop` the pick `h00` (centre 1.16%) has the lowest median of the capture, below `h02` (centre 8.3%) and `h04` (centre 13.9%, rotation 3.9°).
- **On the same surfaces.** Compared over the clusters both candidates of a pair score, where they share at least 30, the passing candidate reads higher in 23 of 33 pairs against the bitmap and 21 of 35 pairwise. The other pairs share fewer than 30 clusters, because the candidates of a capture pose different parts of it.
- **The cell plane normal matters.** With the mean viewing direction as every point's normal, the same readings order 0.541 to 0.624 of the pairs correctly, and the highest-scoring candidate passes on 3 or 4 of the 8 captures.
- **Stability.** Over the 26 perturbation-suite runs, 165 pairs, the readings order 0.830 to 0.921 of the pairs correctly, against 0.964 for the scored-point count. On `SeoulBull` every reading orders every pair correctly on all 13 files, and on the two files whose first candidate is seeded from frames 12 to 16 the bitmap median ranks the other passing candidate, `h01`, first. On `KerryPark480` the highest-scoring candidate with a verdict passes on 2 of the 13 files under the bitmap median, 7 under the pairwise share and 12 under the scored-point count (one file has no passing candidate). The bitmap median ranks the failing `h04` first on the base, repeat and reorder files, and under the pairwise share a two-frame candidate with no ground-truth verdict scores highest on 10 of the 13 files.
- **Cost.** Building the frames takes under 2 s per capture. Rendering and reading every candidate, with the posed photographs decoded once per capture, takes 0.8 s (`SeoulBull`, 5 candidates) to 10 s (`MurdoSmallAntiqueCat`, 10), and 15 to 62 s over three runs on `DnDTabletop` (Windows 11, Intel Core i9-14900HX).

**What this settles.** Blur matching changes the score by too little to matter: the blur-matched and plain readings differ by at most one correctly ordered pair in 109 under every aggregation. None of the twelve combinations of reading and aggregation orders the candidates better than the number of points scored, so the score as built here does not give qualification or the rank a signal they lack.

## Implementation notes

- The score must sample views at the patch's stored keypoints or at the frames derived from them, never at reprojections of the candidate's points, per the project's rule for photometric vetting. A reprojection carries the candidate's error into the sample and hides what the score is meant to find.
- The seed's existing rotation-only photometric channel reads windows sized by the stored affine shapes. The new score reads rendered tiles through frames. The two should share the renderer, not the window code.
- The gate and the rank weights are capture-relative and derived from the release's own distribution of scores. No constant is chosen here.

## Parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `score_tier` | `R = 12` | Patch resolution the score renders at. |
| `score_gate_quantile` | to be measured | Quantile of the release's scores below which a candidate does not qualify. |
| `min_scored_points` | to be measured | A candidate with fewer scorable points is not gated on the score. |

## Testing

On `seoul_bull_sculpture`: the two correct half-capture candidates score above the three-frame candidate whose focal is 11% low. On `kerry_park` against the checked-in ground truth: the candidate with 28% centre error scores below the one with 1%, and the flipped-camera candidate scores lowest. Across the eight captures: the pick agrees with the truth on all eight, as it does without the score.

## Non-goals

Scoring at the fine tier. The score separates candidates; it does not localise.

Replacing the geometric signals. Reach and focal spread stay as gates; the score is added, not substituted.

## Open questions

- **Aggregation.** Median over points, a trimmed mean, or the share of points above a capture-relative bar. [Measured](#measured) on the eight captures: the share of pairwise agreement and the trimmed means lead, at 0.716 to 0.734 of the pairs ordered correctly, and none does better than the count of points scored. Which aggregation to use stays open until a reading separates beyond that count.
- **Blur-matched or plain.** Settled by [Measured](#measured): the two order the candidates alike, so the choice does not decide the score's separation.
- **Weight against reach.** Whether a high-reach candidate with a mediocre score outranks a low-reach one with a good score. The ground-truth set should decide.
- **Robustness to small changes in member shapes.** The current pick turns on sub-percent changes in the cluster-patches file: on `KerryPark480`, moving under 0.4% of the kept members changes which candidates are committed and which pass, and two files that differ in a few dozen moved members give a first qualified candidate that fails and one that passes ([measurements](../core/patch/cluster-patch-refinement-measurements.md#subset-with-the-loop-as-the-default-2026-10-09), [seed-hypothesis-loop.md](../core/geometry/seed-hypothesis-loop.md#rank)). The photometric score must rank the candidates of two such files alike, or rank better on the file whose shapes are better; how to test that beyond the two `KerryPark480` files is open.
- **Rotation-only members.** They have no depth, so a patch's rendering through them tests rotation and focal only. Either score them the same way with that understood, or leave them to the existing rotation channel.
