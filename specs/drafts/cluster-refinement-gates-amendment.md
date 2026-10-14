# Cluster Refinement Gates Amendment

**Status:** Draft. Amends [core/patch/cluster-patch-refinement.md](../core/patch/cluster-patch-refinement.md) § "Acceptance" and § "The reference". Decided:
- a member whose template support runs partly off the photograph is correlated over its in-frame samples and carries its **coverage**, the share of the support window's samples inside the frame, instead of being refused as `NotEvaluated`. A coverage floor replaces the refusal; the measurement below puts it at 0.75, with 0.9 the level at which the score is indistinguishable from full support;
- a `NotEvaluated` member records why: no reference in its cluster, or support off the frame. Today the one status covers both;
- the coverage floor is not made a default until it is measured again with the member gate at the refined shape on, since the measurement below was taken with that gate off.

Not decided: whether shrinkage of the refined shape joins the shift gate as a second reading of a wandering fit; whether coverage also weights the member downstream (member coherence, the normal's agreement weights); whether a member with coverage under 1 can be a reference candidate, and how the reference rule reads coverage and clipping, which the track reference-view rule already does ([core/patch/reference-view.md](../core/patch/reference-view.md)); whether the ZNCC bar reads the reference's self-similarity radius, in the measured shape of a blur-matched bar near 0.74 at a radius of 0.5 or less, 0.80 at 0.5 to 1, 0.84 at 1 to 2 and 0.86 above 2, lower again at scale ratios of 1/4 to 0.7 for radii under 1 (near 0.70 and 0.79); on held-out ground-truth pairs the bar with the ratio term passes 1.0 to 5.4 points fewer lookalike impostors than a single bar at the same share of right pairs, and the radius-only bar between 0.8 points more and 8.1 fewer ([True pairs across scale](../core/patch/cluster-patch-refinement-measurements.md#true-pairs-across-scale-2026-10-10)). Also not decided: whether the cascade seeds its rotation by scoring the template at 12 rotations of the SIFT seed and cascading from the best one, which on ground-truth pairs brings 80.1% of the right pairs whose SIFT seed is more than 20° off within 20° of the true map, against 50.3% from the SIFT seed, for 3% more time and 0.1 to 0.3 points more lookalike impostors at 0.85, while raising the right pairs reaching 0.85 by at most 1.4 points per scale bin ([Rotation seeding](../core/patch/cluster-patch-refinement-measurements.md#rotation-seeding-2026-10-11)). See [Open questions](#open-questions).

The shift gate in units of the member's scale, which this draft first proposed, is built and described in [core/patch/cluster-patch-refinement.md § "Acceptance"](../core/patch/cluster-patch-refinement.md#acceptance), under "The shift gate".

Related: [cluster-review-strips.md](cluster-review-strips.md) (the sheet and review method these findings came from), [image-photometric-model.md](image-photometric-model.md) (the colour reading beside the ZNCC).

## Purpose

Cluster-patch refinement vets every member of a cluster against its reference and stores a status beside it. A member is not evaluated when any sample of its 12× support window falls outside the photograph, however small the part outside is, so a member near the frame edge is never scored. The support window grows with the feature, so the refused band along each edge is wider for large features and is a larger part of a small image.

On the SeoulBull capture, 17 images of 270×480 px, browsing every cluster on a sheet showed this: a member whose support lost a sliver at the left edge was never scored although it reads 0.98 against its reference once refined. Measurements against the ground-truth poses confirm it and give the replacement rule its numbers.

## The rules

### Coverage in place of the out-of-frame refusal

A member's **coverage** is the share of the support window's samples, both bilinear taps at the selected pyramid level, that fall inside the photograph at the member's geometry. A member with coverage below the floor is `NotEvaluated` with the reason *support off the frame*. A member at or above it is correlated over the in-frame samples only, and an evaluation during the search that lowers coverage by more than 0.05 below the seed's scores worst, so the optimiser cannot gain by moving the window off the image. Coverage is stored per member.

### Reasons on `NotEvaluated`

`NotEvaluated` splits into *no reference* (the cluster is unrefinable) and *support off the frame* (coverage under the floor). The split is a status code, so the `.matches` status table gains two names where it had one; a reader of an older file reads the single code as either.

## Measured

The measurement uses the epipolar residual of a member's stored position against its reference's under the ground-truth poses, with *wrong* above `max(3, 5 · median over kept members)` px, which is 3.0 px on SeoulBull. Wrong shares are lower bounds, since a slip along the epipolar line is invisible to the test.

### Coverage (SeoulBull, 2026-10-10)

The cluster file is refined with the member gate at the refined shape off. The kernel cannot score partial support, so a port of the refinement to Python correlated over in-frame samples; on 2,000 full-coverage controls it reproduces the kernel's stored ZNCC to three decimals and agrees on every accept and reject.

The 386 `NotEvaluated` members in refinable clusters all lost support at an edge. Coverage quantiles (10/50/90%): 0.76 / 0.91 / 0.99; 222 are at 0.9 or above, 353 at 0.75 or above, one under 0.5. Refined, 194 pass the 0.85 ZNCC bar, 22 of them wrong (11.3%), against 6.6% among kept members. The wrong ones sit at the highest coverages, so coverage does not explain the excess; right border members already sit further from their epipolar lines (median 0.32 px against 0.19 px for kept), and most are at the top and side edges in two-member clusters.

Whether a partial-support score is biased was tested by cutting a virtual straight frame edge through each of the 2,000 controls' support at a random angle and re-refining against its own full run, re-weighted to the file's kept-to-rejected mix:

| coverage | right members still accepted | wrong members newly accepted | wrong share among accepted |
|---|---|---|---|
| 1.0 | 100% | 0% | 6.6% |
| 0.95 | 98.6% | 0.6% | 6.3% |
| 0.9 | 98.1% | 0.7% | 6.4% |
| 0.75 | 97.0% | 3.5% | 7.1% |
| 0.6 | 94.1% | 6.3% | 7.6% |
| 0.5 | 90.5% | 8.0% | 8.1% |

Down to 0.9 there is no measurable bias. Below about 0.75 the missing samples lift wrong members more than right ones. A floor of 0.9 admits 110 of the 386 on this file, 0.75 admits 182.

Weighting the correlation by coverage does not replace the floor. Multiplying the score by coverage, or by its square root or fourth root, acts as a soft floor that also refuses right members at 0.9 coverage (keeping 27%, 72% and 86% of them). A Fisher lower confidence bound on the correlation, calibrated to 0.85 at full coverage, raises the bar by at most 0.02 at 0.35 coverage, because the in-frame sample count stays in the hundreds; at 0.5 coverage the wrong members it would need to remove score in the same range as the right ones. Under every rule the wrong share among accepted real border members stays between 11% and 15%.

## Open questions

- **Shrinkage.** Whether the refined-over-detected scale ratio `√|det S| / √|det A|` (`S` the refined shape, `A` the SIFT detection's) is read beside the shift gate as a second sign of a wandering fit, or whether the gates at the refined shape already cover it. Among SeoulBull's members refused by the earlier 3 px shift gate, the ratio is 0.66 at the median for the wrong ones and 1.13 for the right ones ([Shift in scale units](../core/patch/cluster-patch-refinement-measurements.md#shift-in-scale-units-2026-10-10)).
- **Coverage downstream.** Whether a member's coverage weights it in member coherence and in the normal's agreement weights, as the self-similarity radius does.
- **Edge members as references.** Today a candidate whose template cannot be built is skipped silently, and a cluster whose only gate-passing members are at the border is unrefinable (SeoulBull cluster 33). Whether the reference rule reads coverage and clipping as the track rule does, and what floor a reference needs.
- **The border population.** Accepted border members are wrong about twice as often as kept members at every coverage. Whether that is the detector at the frame edge or the two-member clusters they sit in.
- **The bar across scale.** Right ground-truth pairs refined as members reach 0.85 less often below a scale ratio of 0.7, most at sharp references and partly through the SIFT seed's rotation, and a bar that reads the reference's self-similarity radius trades right pairs for impostors only slightly better than a single bar ([True pairs across scale](../core/patch/cluster-patch-refinement-measurements.md#true-pairs-across-scale-2026-10-10)).
- **The rotation seed.** Scoring 12 rotations of the SIFT seed and cascading from the best recovers most badly seeded right pairs at almost no cost, and an angular-signature correlation or more cascades add at most 2 points of recovery for 1.5 to 4.7 times the time and more impostors; whether the recovered pairs move the bar question above is open, since the right pairs' pass rate below a ratio of 0.7 rises by only 0.1 to 1.4 points under it ([Rotation seeding](../core/patch/cluster-patch-refinement-measurements.md#rotation-seeding-2026-10-11)).

## Non-goals

- Changing the ZNCC bar or the self-similarity gates in this amendment. The measurements read them as they are; a bar that reads the self-similarity radius is listed as not decided above.
- A format change beyond the two new status names and the stored coverage.
