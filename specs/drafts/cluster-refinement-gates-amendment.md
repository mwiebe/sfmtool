# Cluster Refinement Gates Amendment

**Status:** Draft. Amends [core/patch/cluster-patch-refinement.md](../core/patch/cluster-patch-refinement.md) § "Acceptance" and § "The reference". Decided:
- a member whose template support runs partly off the photograph is correlated over its in-frame samples and carries its **coverage**, the share of the support window's samples inside the frame, instead of being refused as `NotEvaluated`. A coverage floor replaces the refusal; the measurement below puts it at 0.75, with 0.9 the level at which the score is indistinguishable from full support;
- the shift gate reads the refined translation in units of the member's own scale, `|t| ≤ k · √|det A|`, instead of a fixed count of source pixels. The measurement below supports `k = 2`;
- a `NotEvaluated` member records why: no reference in its cluster, or support off the frame. Today the one status covers both;
- no default changes until each rule is measured again with the member gate at the refined shape on, because that gate catches most of the members the shift gate catches, and the two overlap.

Not decided: whether shrinkage of the refined shape joins the shift as a second reading of a wandering fit; whether coverage also weights the member downstream (member coherence, the normal's agreement weights); whether a member with coverage under 1 can be a reference candidate, and how the reference rule reads coverage and clipping, which the track reference-view rule already does ([core/patch/reference-view.md](../core/patch/reference-view.md)). See [Open questions](#open-questions).

Related: [cluster-review-strips.md](cluster-review-strips.md) (the sheet and review method these findings came from), [image-photometric-model.md](image-photometric-model.md) (the colour reading beside the ZNCC).

## Purpose

Cluster-patch refinement vets every member of a cluster against its reference and stores a status beside it. Two of its gates refuse members by a fixed count of source pixels: a member is not evaluated when any sample of its 12× support window falls outside the photograph, and a member is refused when its refined translation from the SIFT detection exceeds 3 px. Both bars sit on quantities that scale with the feature and the image, so each refuses correct members on small images and small features and admits wrong ones on large ones.

On the SeoulBull capture, 17 images of 270×480 px, browsing every cluster on a sheet showed both at once: a member whose support lost a sliver at the left edge was never scored although it reads 0.98 against its reference once refined, and most of the 81 members refused for shift were visibly the same surface as their reference. Measurements against the ground-truth poses confirm both and give the replacement rules their numbers.

## The rules

### Coverage in place of the out-of-frame refusal

A member's **coverage** is the share of the support window's samples, both bilinear taps at the selected pyramid level, that fall inside the photograph at the member's geometry. A member with coverage below the floor is `NotEvaluated` with the reason *support off the frame*. A member at or above it is correlated over the in-frame samples only, and an evaluation during the search that lowers coverage by more than 0.05 below the seed's scores worst, so the optimiser cannot gain by moving the window off the image. Coverage is stored per member.

### Shift in scale units

The shift gate reads `|t| / √|det A|`, the refined translation over the member's refined scale, against `max_shift_scales` (default to be set from the measurement, `2`). A member over it is `RejectedShift` as before. The fixed `max_shift_px` is retired.

### Reasons on `NotEvaluated`

`NotEvaluated` splits into *no reference* (the cluster is unrefinable) and *support off the frame* (coverage under the floor). The split is a status code, so the `.matches` status table gains two names where it had one; a reader of an older file reads the single code as either.

## Measured

Both measurements use the epipolar residual of a member's stored position against its reference's under the ground-truth poses, with *wrong* above `max(3, 5 · median over kept members)` px: 3.0 px on SeoulBull and KerryPark480, 4.68 px on DnDTabletop. Wrong shares are lower bounds, since a slip along the epipolar line is invisible to the test.

### Coverage (SeoulBull, 2026-10-10)

The kernel cannot score partial support, so a port of the refinement to Python correlated over in-frame samples; on 2,000 full-coverage controls it reproduces the kernel's stored ZNCC to three decimals and agrees on every accept and reject.

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

### Shift (SeoulBull, KerryPark480, DnDTabletop, 2026-10-10)

Wrong share by group:

| | kept | rejected low ZNCC | rejected shift |
|---|---|---|---|
| SeoulBull | 6.7% (2,953) | 53.7% (2,288) | 37.0% (81) |
| KerryPark480 | 16.1% (15,609) | 64.9% (5,503) | 38.7% (194) |
| DnDTabletop | 9.5% (1,370,494) | 47.9% (41,599) | 17.9% (112,662) |

In pixels the wrong share barely changes across the 3 px bar: kept members at 2–3 px are wrong 22.8% / 27.4% / 13.1% of the time, shift-rejected members at 3–4 px 33% / 34% / 15.2%. In units of the member's refined scale the share climbs the same way among rejected and kept members alike:

| shift / member scale | SeoulBull | KerryPark480 | DnDTabletop |
|---|---|---|---|
| under 1 | 0 of 22 | 12.4% | 16.2% |
| 1 – 2 | 0 of 8 | 40.5% | 33.4% |
| 2 – 4 | 35% | 62.2% | 70.3% |
| over 4 | 74% | 83.9% | 91.0% |

Replacing the gate, applied to kept plus shift-rejected members, admitted count and wrong share of the kept population:

| rule | SeoulBull | KerryPark480 | DnDTabletop |
|---|---|---|---|
| `|t| ≤ 3 px` (current) | 2,953, 6.67% | 15,609, 16.07% | 1,370,494, 9.51% |
| no shift gate | 3,034, 7.48% | 15,803, 16.35% | 1,483,156, 10.15% |
| `|t| ≤ 1 × scale` | 2,786, 4.38% | 15,316, 15.27% | 1,474,659, 9.98% |
| `|t| ≤ 2 × scale` | 2,938, 5.51% | 15,687, 15.94% | 1,482,082, 10.10% |
| `|t| ≤ 3 × scale` | 2,971, 6.06% | 15,751, 16.14% | 1,482,705, 10.13% |

On SeoulBull `k = 2` admits 30 refused members, none wrong, and refuses 45 kept, 35 wrong. On DnDTabletop, where scales are about 10 px, every scale rule admits nearly all of the 113k shift-rejected members and the kept wrong share rises half a point.

Two readings beside the shift. Among SeoulBull's shift-rejected members the wrong ones have a refined-over-detected scale ratio of 0.66 at the median against 1.13 for the right, so the fit that wanders also shrinks. And the self-similarity radius separates right from wrong among them better than shift or ZNCC: those that would pass the bar of 2.5 are 12.2% / 28.8% / 17.2% wrong, those that would fail 75% / 90.3% / 68.6%. The member gate at the refined shape, which reads that radius, is off in the measured files.

## Open questions

- **Shrinkage.** Whether `√|det A_refined| / √|det A_detected|` joins the shift as a reading of a wandering fit, or whether the refined-shape gate already covers it.
- **Coverage downstream.** Whether a member's coverage weights it in member coherence and in the normal's agreement weights, as the self-similarity radius does.
- **Edge members as references.** Today a candidate whose template cannot be built is skipped silently, and a cluster whose only gate-passing members are at the border is unrefinable (SeoulBull cluster 33). Whether the reference rule reads coverage and clipping as the track rule does, and what floor a reference needs.
- **The border population.** Accepted border members are wrong about twice as often as kept members at every coverage. Whether that is the detector at the frame edge or the two-member clusters they sit in.
- **The DnDTabletop cost.** A scale-relative shift admits 113k members there at an 18% wrong share. Whether the refined-shape gate removes most of them, which decides whether `k = 2` holds across captures.

## Non-goals

- Changing the ZNCC bar or the self-similarity gates. The measurements read them as they are.
- A format change beyond the two new status names and the stored coverage.
