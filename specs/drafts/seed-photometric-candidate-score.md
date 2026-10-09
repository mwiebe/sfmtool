# Seed Photometric Candidate Score

**Status:** Draft. Decided:
- every candidate the seed's hypothesis loop releases is scored photometrically, at the coarse patch tier, by resampling each posed view into the candidate's own patches and reading how well the views agree;
- the score enters the qualification rule and the rank beside the geometric signals they read today, and a candidate whose photometric score is poor cannot qualify however good its reach and focal spread;
- a release carries patch frames, built from the cluster file's per-cell displacements and the candidate's poses, so the score reads the same frames a downstream consumer would;
- the score is validated against the eight ground-truth-backed captures by one number: how often the seed's pick is a candidate that agrees with the ground truth.

Not decided: the aggregation from per-patch agreement to one candidate score; whether the score reads blur-matched or plain ZNCC; how the score is weighed against reach in the rank; whether rotation-only members are scored the same way. See [Open questions](#open-questions).

Amends:
- [core/geometry/seed-hypothesis-loop.md](../core/geometry/seed-hypothesis-loop.md) § "Rank" and § "Product"
- [core/geometry/seed-candidate-evaluation.md](../core/geometry/seed-candidate-evaluation.md): a new channel beside "Exact photometric witness", which today serves rotation-only members only
- [core/geometry/seed-drive.md](../core/geometry/seed-drive.md): the rank the drive reads

Related drafts: [cell-plane-normal-precision-gate.md](cell-plane-normal-precision-gate.md) builds the frames from the cell displacements the piecewise refinement of [cluster-patch-refinement.md](../core/patch/cluster-patch-refinement.md#piecewise-refinement) stores and the normals [cell-plane-normals.md](../core/patch/cell-plane-normals.md) derives from them; [piece-gated-grid-normal.md](piece-gated-grid-normal.md) the normals in them; [two-tier-patch-density.md](two-tier-patch-density.md) the tier the score runs at and the handoff of the pick to the fine tier; [sharper-patch-bitmap.md](sharper-patch-bitmap.md) how a view is scored against a bitmap.

## Purpose

The seed stage takes a set of images and their clustered feature matches and produces a handful of small candidate reconstructions, each a few camera poses and the points they see. Downstream, one of those candidates is grown into the full reconstruction, so choosing the right one matters more than producing many.

Today the choice reads geometry only. A candidate qualifies when it poses enough images, when its reach, the fraction of the capture it can explain, is high enough, when its focal scan showed a clear optimum, and when none of four flags fired. The rank is the order of qualification. None of those signals looks at a photograph.

Scored against the eight captures that have an approved ground truth, that rule picks a candidate agreeing with the truth in four. In the other four a correct candidate exists in the released set and is passed over. The failing candidates share one shape: rotations within a degree and focal within a percent, with the camera centres collapsed to a quarter or a half of their true spread, or one camera flipped. Reprojection residuals do not see this, because a collapsed layout re-triangulates its points to fit. Resampling the photographs into the candidate's patches does see it, because a wrong layout places the same patch at the wrong depth and the views disagree about what it looks like.

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

The measure is the number of ground-truth-backed captures on which the seed's pick agrees with the truth under the existing pass rule: median rotation error under 1°, median camera-centre error under 5% of the ground truth's extent after similarity alignment, focal within 5%. Today that is 4 of 8. The eight are the seven approved fleet entries plus `seoul_bull_sculpture`, which the fleet should enlist. `kerry_park` in the fleet should score against the checked-in ground truth rather than its older references. Each is a fast run and the whole set is under twenty minutes.

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

On `seoul_bull_sculpture`: the two correct half-capture candidates score above the three-frame candidate whose focal is 11% low. On `kerry_park` against the checked-in ground truth: the candidate with 28% centre error scores below the one with 1%, and the flipped-camera candidate scores lowest. Across the eight captures: the pick agrees with the truth on more than four, with no capture regressing from a correct pick to a wrong one.

## Non-goals

Scoring at the fine tier. The score separates candidates; it does not localise.

Replacing the geometric signals. Reach and focal spread stay as gates; the score is added, not substituted.

## Open questions

- **Aggregation.** Median over points, a trimmed mean, or a count of points above an agreement threshold. Measure all three against the eight captures.
- **Blur-matched or plain.** Blur matching removes a penalty on far or blurred views that a collapsed layout does not cause. Likely blur-matched, since the point is to isolate layout error.
- **Weight against reach.** Whether a high-reach candidate with a mediocre score outranks a low-reach one with a good score. The ground-truth set should decide.
- **Robustness to small changes in member shapes.** The current pick turns on sub-percent changes in the cluster-patches file: on `KerryPark480`, moving under 0.4% of the kept members changes which candidates are committed and which pass, and two files that differ in a few dozen moved members give a first qualified candidate that fails and one that passes ([measurements](../core/patch/cluster-patch-refinement-measurements.md#subset-with-the-loop-as-the-default-2026-10-09), [seed-hypothesis-loop.md](../core/geometry/seed-hypothesis-loop.md#rank)). The photometric score must rank the candidates of two such files alike, or rank better on the file whose shapes are better; how to test that beyond the two `KerryPark480` files is open.
- **Rotation-only members.** They have no depth, so a patch's rendering through them tests rotation and focal only. Either score them the same way with that understood, or leave them to the existing rotation channel.
