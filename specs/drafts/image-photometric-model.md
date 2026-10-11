# Image Photometric Model

**Status:** Draft. Decided:
- every photograph in a reconstruction carries a photometric record: a per-image part (exposure, chromaticity) that changes shot to shot, and a per-camera part (response curve, vignetting) that a camera's frames share. The record maps the photograph's stored pixel values to a scene radiance that is comparable across images;
- the model follows PPISP's four stages and their order, applied to linear radiance: exposure as a base-2 scale, per-channel radial vignetting about a centre, a chromaticity homography built from offsets of the R, G, B and white primaries, and a per-channel piecewise power response curve with a gamma. The time-based controller that PPISP trains to predict per-frame values for novel views is not adopted: every image here is a photograph, so its values are fitted, never predicted;
- two additions PPISP does not model: clipped pixels are detected on the photograph and excluded from every fit and counted, and a one-parameter blend toward luminance covers partial desaturation (haze, veiling glare, a flat picture style);
- the first consumer is a colour agreement reading between two blur-matched patch tiles: fit the smallest per-image difference that maps the member's tile onto the reference's, and read the residual after the fit and the plausibility of the fitted parameters. Per-channel ZNCC stays the geometric match score; the colour reading is a second channel beside it;
- the record is per image, and effects that belong to one view of one surface are not in it: glare, specular highlights, lens flare and shadow edges change with the viewing direction and the surface, not the frame. The pairwise fit between two tiles absorbs them, so a true correspondence under glare still reads as a match. Where the per-tile fit lands far from the image's record, the difference is a reading on that observation, a **local photometric deviation**, with glare its first case: a brightening with desaturation, often with clipped samples. It is reported per observation and is a candidate input to the reference-view rule, which today has no way to refuse a glared view as the stored bitmap, and to the agreement weights of normal refinement and member coherence;
- the global ambiguities are fixed by recentring, not by pinning one image: the mean log exposure over an image set is zero, and the mean chromaticity offset is zero;
- no bar is chosen in this draft. Every threshold is measured against the ground-truth epipolar verdicts and a blind human review in the cluster-strips convention before it enters a gate.

Not decided: whether the per-camera response and vignetting are fitted from the photographs alone or seeded from EXIF and camera metadata; the file layout of the record in `.sfmr` and `.matches`; whether the per-pair colour reading is stored per observation like the ZNCC columns; how the record is used by the renderer when it samples a view; the linearising gamma, the rank rule's chromaticity-scatter tolerance and the priors' weight in a pairwise fit, which step 1 measured at γ 1 and 2.2, tolerance 0.01 and weights 1 and 0 without settling any of them ([Measured](#measured-step-1-2026-10-10)). See [Open questions](#open-questions).

Amends, once built:
- [formats/sfmr-file-format.md](../formats/sfmr-file-format.md): a per-image and per-camera photometric record
- [core/patch/blur-matched-zncc.md](../core/patch/blur-matched-zncc.md): a colour agreement reading beside the blur-matched score
- [core/patch/cluster-patch-refinement.md](../core/patch/cluster-patch-refinement.md): the colour reading as a member statistic
- [core/patch/reference-view.md](../core/patch/reference-view.md): colour typicality as a reading the reference-view rule can consult

Related: [cluster-review-strips.md](cluster-review-strips.md) (the review method every bar here is measured by), [sharper-patch-bitmap.md](sharper-patch-bitmap.md) (the blur-matched tiles the pairwise fit runs on), [seed-photometric-candidate-score.md](seed-photometric-candidate-score.md).

## Purpose

A reconstruction's photographs are taken under different exposures, white balances, and often by different cameras with different response curves and lens falloff. The project's patch scores correlate each colour channel separately and normalise it by its own mean and contrast, so an exposure or white-balance change between two shots costs nothing. That choice also discards information. In the blind review of cross-scale cluster members on 2026-10-10 ([measurements](../core/patch/cluster-patch-refinement-measurements.md#blur-matched-scores-across-scales-2026-10-10)), six of the nineteen members the reviewer called wrong were called wrong because their colour did not match the reference, and every one of those was a wrong correspondence by the ground truth. The per-channel scores cannot see this, because a wrong surface with a different colour reads as a white-balance change.

A model of what a camera does to colour separates the two cases. A photographic change between two shots of the same surface is an exposure shift, a white-balance shift, a response curve and lens falloff: a few parameters, the same for every pixel of the image, and inside a range cameras produce. Two different surfaces need a different transform in each channel, or one no camera produces, or leave a residual the model cannot absorb. The model is fitted per image and stored with it, so every consumer that compares two views can ask the same question: after the photographic difference is removed, do these two tiles show the same colours?

## The model

A photograph's stored value `p` in a channel, at pixel `x`, relates to scene radiance `L` by four stages applied in order:

1. **Exposure.** `I_exp = L · 2^Δt`, one offset `Δt` per image.
2. **Vignetting.** `I_vig = I_exp · v(r)`, with `v(r) = clip[0,1](1 + α₁ r² + α₂ r⁴ + α₃ r⁶)` about a centre `μ`, radius `r` measured from `μ` in normalised image coordinates. Per camera, per channel, with `α_j ≤ 0`.
3. **Chromaticity.** A 3×3 homography `H` on RG chromaticity plus intensity. `H` is built from four 2D offsets `Δc_k`, `k ∈ {R, G, B, W}`, of the fixed source chromaticities `R = (1, 0)`, `G = (0, 1)`, `B = (0, 0)`, `W = (⅓, ⅓)`: `H = T · diag(k) · S⁻¹`, normalised so `H₃₃ = 1`, with an intensity normalisation so that `H` changes chromaticity only and leaves exposure to stage 1. Eight parameters per image. The identity is `Δc_k = 0`.
4. **Response.** Per camera, per channel, `g(x) = f₀(x)^γ` with `f₀` the piecewise power S-curve with inflection `ξ` and exponents `τ`, `η` of PPISP, C¹ at the inflection and monotone by construction. Four parameters per channel.

Two readings sit outside the chain:

- **Clipping.** A pixel at or near the top or bottom of range in any channel of the photograph is clipped. It is excluded from every fit and its share is reported (the glossary's **clipped share**). Clipping is read on the photograph's own pixels, never on a resampled tile.
- **Desaturation.** `I_desat = (1 − s) · I + s · Y(I)`, with `Y` the luminance of `I` and one `s ∈ [0, 1)` per image, applied after stage 3. It is zero for a camera that renders colour faithfully and is the one parameter PPISP's chain lacks that a photograph commonly needs.

### Why these stages

PPISP's argument for separating a chromaticity homography from exposure, rather than fitting a per-channel white-point gain, is that a per-channel gain entangles brightness with colour and forces the response curve to absorb value-dependent colour shifts; it measures a lower correlation between fitted exposure and colour offsets with the homography. The same entanglement is what a per-channel ZNCC removes blindly. A homography also covers the gamut difference between two cameras' sensors, which a diagonal cannot.

The response curve and vignetting are per camera because they are properties of the sensor and lens, not the shot. Fitting them once over all of a camera's frames, and only then fitting exposure and chromaticity per image, keeps the per-image fit at nine parameters plus desaturation. For two images of the same camera the response and vignetting cancel in a pairwise comparison and need not be fitted at all.

### The pairwise colour reading

Given two blur-matched tiles of the same patch, a reference `A` and a member `B`, both sampled at the patch resolution and in register, the reading:

1. **Linearises** both tiles through their cameras' response curves and divides out vignetting at their sample positions, where those are known; else treats both as linear with a shared unknown gamma.
2. **Masks** clipped samples on either side and reports the two clipped shares.
3. **Fits** the per-image difference `B → A` over the unmasked samples: exposure difference `ΔΔt`, chromaticity offsets `Δc_k` of `B` relative to `A`, desaturation `s`. The fit is robust (Huber) with PPISP's pulls toward identity as priors: δ = 0.1 on `Δt`, δ = 0.005 on each `Δc_k`.
4. **Reads two numbers.** The **colour residual**: the robust RMS of `A − fit(B)` over the unmasked samples, in the units of `A`'s normalised radiance. The **plausibility** of the fit: whether `|ΔΔt|` is within a stop range, whether the `Δc_k` are within a white-balance range and consistent with a single illuminant shift (the four offsets agree with a 2-parameter white-point shift to within a tolerance, or they do not), whether `s` is small, and whether the clipped pixels of one tile are the bright pixels of the other.

A tile pair with little chromatic variety cannot constrain eight chromaticity offsets. The chromaticity scatter's rank decides: a full homography when the unmasked samples span the chromaticity plane, a white-point shift (two parameters) when they lie along a line, exposure only when they are a point. The reading reports which was fitted.

### Where the record lives

Per image: `Δt`, the four `Δc_k`, `s`, and the image's clipped share. Per camera: `μ`, the three `α_j` per channel, and `ξ`, `τ`, `η`, `γ` per channel. Both in the `.sfmr` file beside the camera and image tables, and in a `.matches` file's cluster-patches block when the record was fitted at that stage. A reader without a record treats every image as identity, which is what every reader does today.

The recentring rule fixes the gauge: over the image set the mean `Δt` is zero and the mean `Δc_k` is zero for each `k`, so no image's exposure is treated as truth. Pinning one image to the identity was tried in an earlier colour-correction experiment and was worse: when the pinned image is a poor representative, every other image drifts to match it.

## Plan

Each step is measured before the next is designed, and no step changes a default.

1. **Pairwise reading on the reviewed cases.** A script fits the pairwise model on the 48 reviewed cases of 2026-10-10 and on the full rescued and accepted sets of KerryPark480 and DnDTabletop, without per-camera stages (same camera within each entry; DnDTabletop's response unknown, so a shared gamma is fitted). Readings: colour residual and plausibility against the epipolar verdicts and the human verdicts, per scale-ratio bin. Question: does the colour residual separate the 19 human-wrong rescued members from the 18 human-right ones, and the epipolar-wrong accepted members from the core, where per-channel ZNCC does not?
2. **Second blind review.** Drawn from the cases where colour and ZNCC disagree, in the cluster-strips convention, so the reviewer judges the margin.
3. **Per-camera stages.** Fit response and vignetting per camera over an entry's frames, from the photographs alone, and compare with EXIF-seeded fits where EXIF exists. Measure on kerry_park, whose two fisheyes share a rig and should agree, and on a mixed-camera fleet entry.
4. **The record in the file.** Format change in `.sfmr` and the cluster-patches block, written by `cluster-patches` and `embed-patches`, read by the renderer and the scores. Written as an amendment draft to the format spec once step 1 shows the reading is worth storing.
5. **Consumers.** In order of least risk: a member statistic in cluster-patch refinement (reported, not gated); a reading the reference-view rule can consult; the renderer sampling a view through the record so every tile is compared in radiance; a gate, measured.

### Measured: step 1 (2026-10-10)

**Method.** [`scripts/colour_reading.py`](../../scripts/colour_reading.py) `measure`, on the two cluster-patches files and the `score` results of [Blur-matched scores across scales](../core/patch/cluster-patch-refinement-measurements.md#blur-matched-scores-across-scales-2026-10-10) (branch `bootstrap-core-migration` at `b7794624`). Each pair is the reference's tile and the member's tile rendered as `cluster_strips.py` renders them, the sharper blurred by the pair rule; the plain and blur-matched ZNCC recomputed on these tiles equal the `score` values to a largest difference of 0.0 on every pair.

- **Populations.** Every rescued member (217 and 4,939), a sample of accepted members of up to 2,000 per scale-ratio bin, topped up from the largest bin to 5,000 (5,000 and 6,088), and the 48 reviewed cases. The accepted sample is stratified, so its *all* column weights each bin by its share of the entry's accepted members. Verdicts are the epipolar rule's (bar 3.00 px and 4.68 px).
- **Linearisation.** No response curve is fitted: each tile is read at γ = 1.0 and at γ = 2.2, `(v/255)^γ`, both tiles at the same γ, and divided by `A`'s window-weighted mean intensity. No shared γ was fitted per entry.
- **Clipping.** A photograph pixel is clipped when any channel is 0 or 255. The clipped masks are reduced with the photograph's box pyramid and sampled through the tile's grid, so a sample is masked when any photograph pixel under its taps is clipped; each tile's clipped share is read at full resolution inside its outline.
- **Fit.** Over the unmasked samples of the scoring window's support, weighted by the window: `ΔΔt`, the four `Δc_k` building `H = T diag(k) S⁻¹` (`H₃₃ = 1`, each sample's `R + G + B` kept), and `s ∈ [0, 0.99]`. The data term is the window-weighted mean of a Huber loss whose δ is 1.345 times the MAD scale of a least-squares fit's residuals; the priors are Huber losses with δ = 0.1 on `ΔΔt` and δ = 0.005 on each `Δc_k`'s length, at PPISP's weight of 1, and, in a second run, at weight 0. L-BFGS-B, numerical gradient. `ΔΔt` is bounded to ±8 stops; the weight-1 runs were made before that bound was added, and their largest fitted `ΔΔt` is 5.30 stops, inside it.
- **Rank rule.** The chromaticity scatter's standard deviations along its principal axes are read for each tile and the smaller tile's kept; with tolerance 0.01, both at least 0.01 fits all ten parameters, the first only fits `ΔΔt` and `Δc_W` (R, G and B fixed), neither fits `ΔΔt` alone.
- **Readings.** *Colour residual*: the window-weighted RMS of `A − fit(B)` per sample, the worst 10% by weight trimmed, in `A`'s normalised radiance. *Gain-only residual*: the same after three per-channel gains fitted by Huber IRLS with no priors. *Chromaticity residual*: the same on RG chromaticities. *White-point disagreement*: the RMS over the four primaries of the distance of the fitted targets from the best per-channel-gain homography's, which moves only W. *Any flag*: `ΔΔt` beyond 1 stop, any `Δc_k` longer than 0.05, white-point disagreement above 0.02, or `s` above 0.1. *Clip coincidence*: of the samples clipped in exactly one tile, the share whose other tile is in its top quartile of intensity (for a 255) or bottom quartile (for a 0).
- **Ordering quality** is the share of (right, wrong) pairs that a reading orders with the wrong member on the wrong side (a higher residual, a lower ZNCC, a raised flag), ties counted half: 0.5 is no separation, 1.0 full.

**Ordering quality against the epipolar verdict, prior weight 1.**

`KerryPark480`, accepted, members right / wrong by bin: < 1/4 5 / 77, 1/4–1/2 151 / 235, 1/2–0.7 841 / 291, 0.7–1.4 2907 / 489.

| reading | < 1/4 | 1/4–1/2 | 1/2–0.7 | 0.7–1.4 | all |
|---|---|---|---|---|---|
| blur-matched ZNCC | 0.239 | 0.583 | 0.639 | 0.711 | 0.699 |
| plain ZNCC | 0.317 | 0.542 | 0.624 | 0.710 | 0.704 |
| colour residual, γ 1 | 0.592 | 0.789 | 0.788 | 0.445 | 0.540 |
| colour residual, γ 2.2 | 0.610 | 0.747 | 0.732 | 0.433 | 0.523 |
| gain-only residual, γ 1 | 0.545 | 0.768 | 0.754 | 0.408 | 0.509 |
| gain-only residual, γ 2.2 | 0.592 | 0.756 | 0.729 | 0.398 | 0.497 |
| chromaticity residual, γ 2.2 | 0.512 | 0.699 | 0.711 | 0.451 | 0.535 |
| size of ΔΔt, γ 2.2 | 0.860 | 0.642 | 0.603 | 0.606 | 0.618 |
| length of Δc_W, γ 2.2 | 0.496 | 0.699 | 0.708 | 0.635 | 0.668 |
| white-point disagreement, γ 2.2 | 0.294 | 0.393 | 0.391 | 0.302 | 0.340 |
| s, γ 2.2 | 0.197 | 0.490 | 0.500 | 0.394 | 0.426 |
| any flag, γ 2.2 | 0.457 | 0.620 | 0.629 | 0.475 | 0.528 |
| clip coincidence | 0.420 | 0.480 | 0.490 | 0.440 | 0.448 |

`KerryPark480`, rescued, members right / wrong by bin: < 1/4 0 / 13, 1/4–1/2 27 / 45, 1/2–0.7 52 / 22, 0.7–1.4 40 / 18.

| reading | < 1/4 | 1/4–1/2 | 1/2–0.7 | 0.7–1.4 | all |
|---|---|---|---|---|---|
| blur-matched ZNCC | – | 0.602 | 0.517 | 0.303 | 0.499 |
| plain ZNCC | – | 0.471 | 0.690 | 0.603 | 0.620 |
| colour residual, γ 1 | – | 0.802 | 0.862 | 0.929 | 0.850 |
| colour residual, γ 2.2 | – | 0.733 | 0.820 | 0.932 | 0.822 |
| gain-only residual, γ 1 | – | 0.758 | 0.762 | 0.860 | 0.806 |
| gain-only residual, γ 2.2 | – | 0.754 | 0.765 | 0.854 | 0.796 |
| chromaticity residual, γ 2.2 | – | 0.690 | 0.798 | 0.872 | 0.788 |
| size of ΔΔt, γ 2.2 | – | 0.646 | 0.656 | 0.767 | 0.681 |
| length of Δc_W, γ 2.2 | – | 0.779 | 0.751 | 0.593 | 0.691 |
| white-point disagreement, γ 2.2 | – | 0.508 | 0.427 | 0.332 | 0.421 |
| s, γ 2.2 | – | 0.516 | 0.535 | 0.351 | 0.465 |
| any flag, γ 2.2 | – | 0.700 | 0.669 | 0.664 | 0.649 |
| clip coincidence | – | 0.707 | 0.473 | 0.309 | 0.516 |

`DnDTabletop`, accepted, members right / wrong by bin: < 1/4 1 / 216, 1/4–1/2 694 / 1174, 1/2–0.7 1364 / 636, 0.7–1.4 1823 / 177.

| reading | < 1/4 | 1/4–1/2 | 1/2–0.7 | 0.7–1.4 | all |
|---|---|---|---|---|---|
| blur-matched ZNCC | 0.880 | 0.692 | 0.594 | 0.598 | 0.607 |
| plain ZNCC | 0.889 | 0.662 | 0.587 | 0.599 | 0.609 |
| colour residual, γ 1 | 0.894 | 0.867 | 0.632 | 0.647 | 0.657 |
| colour residual, γ 2.2 | 0.370 | 0.875 | 0.659 | 0.649 | 0.659 |
| gain-only residual, γ 1 | 0.843 | 0.885 | 0.650 | 0.606 | 0.616 |
| gain-only residual, γ 2.2 | 0.755 | 0.888 | 0.656 | 0.612 | 0.622 |
| chromaticity residual, γ 2.2 | 0.861 | 0.815 | 0.607 | 0.632 | 0.641 |
| size of ΔΔt, γ 2.2 | 0.551 | 0.519 | 0.502 | 0.625 | 0.628 |
| length of Δc_W, γ 2.2 | 0.847 | 0.680 | 0.551 | 0.589 | 0.597 |
| white-point disagreement, γ 2.2 | 0.588 | 0.602 | 0.535 | 0.498 | 0.500 |
| s, γ 2.2 | 0.671 | 0.622 | 0.532 | 0.479 | 0.481 |
| any flag, γ 2.2 | 0.924 | 0.667 | 0.564 | 0.554 | 0.560 |
| clip coincidence | 0.094 | 0.368 | 0.426 | 0.519 | 0.513 |

`DnDTabletop`, rescued, members right / wrong by bin: < 1/4 1 / 79, 1/4–1/2 41 / 147, 1/2–0.7 234 / 177, 0.7–1.4 3274 / 986.

| reading | < 1/4 | 1/4–1/2 | 1/2–0.7 | 0.7–1.4 | all |
|---|---|---|---|---|---|
| blur-matched ZNCC | 0.329 | 0.680 | 0.677 | 0.522 | 0.567 |
| plain ZNCC | 0.468 | 0.613 | 0.532 | 0.490 | 0.511 |
| colour residual, γ 1 | 0.949 | 0.929 | 0.737 | 0.595 | 0.667 |
| colour residual, γ 2.2 | 0.949 | 0.941 | 0.759 | 0.598 | 0.672 |
| gain-only residual, γ 1 | 0.975 | 0.955 | 0.750 | 0.540 | 0.628 |
| gain-only residual, γ 2.2 | 0.975 | 0.958 | 0.747 | 0.535 | 0.624 |
| chromaticity residual, γ 2.2 | 0.937 | 0.906 | 0.712 | 0.618 | 0.677 |
| size of ΔΔt, γ 2.2 | 0.620 | 0.645 | 0.576 | 0.581 | 0.600 |
| length of Δc_W, γ 2.2 | 0.956 | 0.731 | 0.610 | 0.509 | 0.568 |
| white-point disagreement, γ 2.2 | 0.690 | 0.651 | 0.546 | 0.482 | 0.504 |
| s, γ 2.2 | 0.595 | 0.633 | 0.545 | 0.487 | 0.505 |
| any flag, γ 2.2 | 0.956 | 0.684 | 0.623 | 0.531 | 0.585 |
| clip coincidence | – | 0.087 | 0.325 | 0.467 | 0.401 |

**Medians, prior weight 1.**

`KerryPark480`, median of right / wrong members:

| population | reading | < 1/4 | 1/4–1/2 | 1/2–0.7 | 0.7–1.4 | all |
|---|---|---|---|---|---|---|
| accepted | blur-matched ZNCC | 0.888 / 0.903 | 0.927 / 0.912 | 0.928 / 0.908 | 0.937 / 0.904 | 0.934 / 0.907 |
| accepted | colour residual, γ 2.2 | 0.251 / 0.317 | 0.218 / 0.301 | 0.201 / 0.281 | 0.113 / 0.065 | 0.136 / 0.222 |
| accepted | gain-only residual, γ 2.2 | 0.234 / 0.314 | 0.212 / 0.300 | 0.185 / 0.274 | 0.106 / 0.046 | 0.128 / 0.205 |
| accepted | size of ΔΔt, γ 2.2 | 0.006 / 0.032 | 0.017 / 0.053 | 0.019 / 0.050 | 0.012 / 0.021 | 0.014 / 0.027 |
| accepted | length of Δc_W, γ 2.2 | 0.078 / 0.110 | 0.019 / 0.118 | 0.017 / 0.098 | 0.007 / 0.014 | 0.009 / 0.033 |
| rescued | blur-matched ZNCC | – / 0.865 | 0.874 / 0.864 | 0.867 / 0.866 | 0.865 / 0.874 | 0.867 / 0.867 |
| rescued | colour residual, γ 2.2 | – / 0.269 | 0.216 / 0.346 | 0.191 / 0.324 | 0.182 / 0.392 | 0.189 / 0.345 |
| rescued | gain-only residual, γ 2.2 | – / 0.251 | 0.215 / 0.351 | 0.151 / 0.305 | 0.181 / 0.346 | 0.172 / 0.340 |
| rescued | size of ΔΔt, γ 2.2 | – / 0.184 | 0.019 / 0.084 | 0.018 / 0.073 | 0.023 / 0.582 | 0.019 / 0.153 |
| rescued | length of Δc_W, γ 2.2 | – / 0.074 | 0.013 / 0.116 | 0.018 / 0.132 | 0.012 / 0.091 | 0.015 / 0.098 |

`DnDTabletop`, median of right / wrong members:

| population | reading | < 1/4 | 1/4–1/2 | 1/2–0.7 | 0.7–1.4 | all |
|---|---|---|---|---|---|---|
| accepted | blur-matched ZNCC | 0.950 / 0.909 | 0.937 / 0.907 | 0.944 / 0.932 | 0.974 / 0.966 | 0.957 / 0.919 |
| accepted | colour residual, γ 2.2 | 0.346 / 0.315 | 0.157 / 0.285 | 0.155 / 0.191 | 0.070 / 0.101 | 0.111 / 0.247 |
| accepted | gain-only residual, γ 2.2 | 0.241 / 0.304 | 0.098 / 0.266 | 0.093 / 0.125 | 0.061 / 0.074 | 0.079 / 0.220 |
| accepted | size of ΔΔt, γ 2.2 | 0.073 / 0.122 | 0.061 / 0.067 | 0.071 / 0.082 | 0.011 / 0.019 | 0.026 / 0.066 |
| accepted | length of Δc_W, γ 2.2 | 0.022 / 0.161 | 0.042 / 0.110 | 0.025 / 0.030 | 0.004 / 0.009 | 0.011 / 0.065 |
| rescued | blur-matched ZNCC | 0.866 / 0.874 | 0.877 / 0.867 | 0.886 / 0.872 | 0.886 / 0.883 | 0.886 / 0.878 |
| rescued | colour residual, γ 2.2 | 0.148 / 0.325 | 0.128 / 0.320 | 0.137 / 0.213 | 0.126 / 0.149 | 0.127 / 0.177 |
| rescued | gain-only residual, γ 2.2 | 0.076 / 0.284 | 0.087 / 0.288 | 0.087 / 0.145 | 0.088 / 0.092 | 0.088 / 0.105 |
| rescued | size of ΔΔt, γ 2.2 | 0.222 / 0.714 | 0.061 / 0.185 | 0.077 / 0.129 | 0.038 / 0.066 | 0.040 / 0.079 |
| rescued | length of Δc_W, γ 2.2 | 0.000 / 0.144 | 0.042 / 0.112 | 0.024 / 0.043 | 0.012 / 0.014 | 0.013 / 0.021 |

**Against the human verdict.** Ordering quality on the 48 reviewed cases (23 right and 22 wrong answered by the reviewer) and on the rescued members among them (18 right, 19 wrong), against the reviewer and against the epipolar verdict (22 right and 26 wrong over the 48, 17 and 23 over the 40 rescued), prior weight 1:

| reading | 48: human | 48: epipolar | rescued: human | rescued: epipolar | rescued median, human right / wrong |
|---|---|---|---|---|---|
| blur-matched ZNCC | 0.711 | 0.722 | 0.734 | 0.752 | 0.877 / 0.863 |
| plain ZNCC | 0.670 | 0.661 | 0.626 | 0.609 | 0.836 / 0.827 |
| colour residual, γ 1 | 0.879 | 0.869 | 0.857 | 0.847 | 0.071 / 0.232 |
| colour residual, γ 2.2 | 0.879 | 0.876 | 0.851 | 0.852 | 0.121 / 0.329 |
| gain-only residual, γ 1 | 0.872 | 0.867 | 0.860 | 0.857 | 0.051 / 0.188 |
| gain-only residual, γ 2.2 | 0.868 | 0.865 | 0.857 | 0.854 | 0.100 / 0.332 |
| chromaticity residual, γ 2.2 | 0.872 | 0.876 | 0.827 | 0.839 | 0.020 / 0.049 |
| size of ΔΔt, γ 2.2 | 0.709 | 0.678 | 0.667 | 0.632 | 0.019 / 0.170 |
| length of Δc_W, γ 2.2 | 0.700 | 0.694 | 0.693 | 0.680 | 0.025 / 0.122 |
| white-point disagreement, γ 2.2 | 0.526 | 0.545 | 0.551 | 0.574 | 0.000 / 0.000 |
| s, γ 2.2 | 0.509 | 0.539 | 0.490 | 0.528 | 0.000 / 0.000 |
| any flag, γ 2.2 | 0.711 | 0.729 | 0.700 | 0.715 | 0.000 / 1.000 |
| clip coincidence | 0.374 | 0.420 | 0.296 | 0.358 | 0.935 / 0.412 |

**Prior weight.** At weight 1 the priors outweigh the data term on a pair: on a synthetic pair whose member is the reference at 1.5 times the exposure, the fit returns `ΔΔt` = −0.442 (−0.585 is exact) and a residual of 0.111; at weight 0 it returns −0.585 and 0. At weight 0 the full homography is unconstrained: the median largest `Δc_k` length of `KerryPark480`'s right accepted members is 180.3 below 1/4, 5.52 from 1/4 to 1/2 and 2.74 from 1/2 to 0.7. Ordering quality over all bins:

| entry | population | reading | weight 1 | weight 0 |
|---|---|---|---|---|
| `KerryPark480` | accepted | colour residual, γ 2.2 | 0.523 | 0.494 |
| `KerryPark480` | accepted | length of Δc_W, γ 2.2 | 0.668 | 0.616 |
| `KerryPark480` | accepted | size of ΔΔt, γ 2.2 | 0.618 | 0.611 |
| `KerryPark480` | rescued | colour residual, γ 2.2 | 0.822 | 0.816 |
| `KerryPark480` | rescued | length of Δc_W, γ 2.2 | 0.691 | 0.652 |
| `KerryPark480` | rescued | size of ΔΔt, γ 2.2 | 0.681 | 0.730 |
| `DnDTabletop` | accepted | colour residual, γ 2.2 | 0.659 | 0.625 |
| `DnDTabletop` | accepted | length of Δc_W, γ 2.2 | 0.597 | 0.549 |
| `DnDTabletop` | accepted | size of ΔΔt, γ 2.2 | 0.628 | 0.650 |
| `DnDTabletop` | rescued | colour residual, γ 2.2 | 0.672 | 0.639 |
| `DnDTabletop` | rescued | length of Δc_W, γ 2.2 | 0.568 | 0.552 |
| `DnDTabletop` | rescued | size of ΔΔt, γ 2.2 | 0.600 | 0.630 |

At weight 0 the colour residual orders the 48 cases at 0.887 (γ 1) and 0.881 (γ 2.2) against the reviewer, and the rescued 37 at 0.871 and 0.863.

**The model against the gain baseline.** The model's intensity normalisation keeps each sample's `R + G + B`, so a per-channel gain is outside it: on a synthetic pair whose member is the reference times (1.3, 1.0, 0.7), the three-gain fit leaves 0.000 and the full model 0.061 at weight 0. On the measured pairs at weight 1, the model's untrimmed RMS is above the gain fit's on 4,673 of 5,219 `KerryPark480` pairs at γ 1 (rank 2: 1,019 of 1,426) and 4,165 at γ 2.2 (rank 2: 1,853 of 2,763), and on 10,501 of 11,031 `DnDTabletop` pairs at γ 1 (rank 2: 770 of 1,088) and 9,774 at γ 2.2 (rank 2: 3,600 of 4,648).

**Rank fallbacks.** Pairs fitted at rank 0 / 1 / 2, at chromaticity-scatter tolerances 0.005, 0.01 (the one fitted) and 0.02; the counts at 0.005 and 0.02 are read from the stored singular values:

| entry | population | pairs | γ 1, 0.005 | γ 1, 0.01 | γ 1, 0.02 | γ 2.2, 0.005 | γ 2.2, 0.01 | γ 2.2, 0.02 |
|---|---|---|---|---|---|---|---|---|
| `KerryPark480` | accepted | 5000 | 190 / 2230 / 2578 | 639 / 2981 / 1378 | 1848 / 2953 / 197 | 87 / 1445 / 3466 | 179 / 2166 / 2653 | 604 / 2935 / 1459 |
| `KerryPark480` | rescued | 217 | 39 / 75 / 103 | 62 / 108 / 47 | 96 / 113 / 8 | 24 / 34 / 159 | 37 / 72 / 108 | 59 / 111 / 47 |
| `DnDTabletop` | accepted | 6088 | 313 / 2970 / 2805 | 1789 / 3432 / 867 | 4323 / 1656 / 109 | 87 / 671 / 5330 | 391 / 2409 / 3288 | 2407 / 2816 / 865 |
| `DnDTabletop` | rescued | 4939 | 732 / 3235 / 972 | 2419 / 2299 / 221 | 4288 / 604 / 47 | 331 / 1110 / 3498 | 947 / 2635 / 1357 | 3507 / 1208 / 224 |

Two accepted `KerryPark480` pairs had fewer than 20 unmasked samples and no reading.

**Clipping.** The median clipped share is 0.009 (reference) and 0.007 (member) over `KerryPark480`'s accepted sample and 0 in its rescued set and in both `DnDTabletop` populations. A member tile with a clipped share above 0: 3,028 of 5,000 and 72 of 217 on `KerryPark480`, 1,013 of 6,088 and 550 of 4,939 on `DnDTabletop`.

**Glare-like pairs.** One tile at least 0.5 stop brighter than the other (window-weighted mean intensity at γ 2.2, over every sample with data, clipped ones included), its saturation (the mean distance of its chromaticities from white) at most 0.8 times the other's, and blur-matched ZNCC at least 0.85; the second pair of columns adds a clipped share above 0 on the brighter tile. The accepted counts are of the stratified sample. *Side* counts the brighter tile, reference / member; *s > 0.1* counts the clipped glare-like pairs whose fitted `s` exceeds 0.1, which can only happen when the reference is the brighter tile, since the fit desaturates the member only.

| entry | population | bright, desaturated: right | wrong | and clipped: right | wrong | side, reference / member | s > 0.1 |
|---|---|---|---|---|---|---|---|
| `KerryPark480` | accepted | 69 of 3906 | 174 of 1094 | 32 of 3906 | 62 of 1094 | 50 / 44 | 42 |
| `KerryPark480` | rescued | 7 of 119 | 27 of 98 | 2 of 119 | 8 of 98 | 4 / 6 | 4 |
| `DnDTabletop` | accepted | 459 of 3883 | 553 of 2205 | 23 of 3883 | 74 of 2205 | 54 / 43 | 29 |
| `DnDTabletop` | rescued | 338 of 3550 | 275 of 1389 | 30 of 3550 | 27 of 1389 | 45 / 12 | 7 |

**Timing.** One pair, both γ, rendering excluded, on one thread: median 24.2 ms and 90th percentile 111 ms on `KerryPark480`, 14.3 ms and 88 ms on `DnDTabletop`, at weight 1. The fit's median at γ 2.2 by rank: 2.0, 5.3 and 39.5 ms on `KerryPark480`, 2.0, 6.4 and 32.3 ms on `DnDTabletop`. At weight 0 a rank-2 fit takes a median of 152 ms and 143 ms. The blur takes a median under 0.1 ms.

**Step-2 set.** `scripts/colour_reading.py review`, seed `20261011`, at γ 2.2 and weight 1, drew 20 members whose colour residual is above the median of the entry's epipolar-wrong accepted and rescued members (0.236 and 0.213) while their blur-matched ZNCC reaches 0.85, from pools of 1,397 and 2,465, and 20 whose colour residual is at most the median of the epipolar-right ones (0.137 and 0.120) while their blur-matched ZNCC misses 0.85, from pools of 263 and 168. The second kind comes from a further `measure --rejected` run over rejected low-ZNCC members that pass the later gates and whose blur-matched score also misses the bar (up to 500 per bin: 1,716 and 1,993 members). Ten of each kind per entry; none is one of the 48 reviewed cases. By the epipolar rule, 12 of the first kind are wrong and 8 right, and 8 of the second kind wrong and 12 right. The cases are rendered in the cluster-strips convention, with `cases.csv` the key; the answers are in [Measured: step 2](#measured-step-2-2026-10-10).

### Measured: step 2 (2026-10-10)

**Data.** The 40 cases of the step-2 set drawn in step 1, answered blind by the maintainer on the cluster-strips page, *right*, *wrong* or *unsure* with an optional note, before the key was opened. The joined record is [cluster-patch-refinement-human-review-2026-10-10b.csv](../core/patch/cluster-patch-refinement-human-review-2026-10-10b.csv): case, dataset, type, scale-ratio bin, member and cluster, scale ratio, plain and blur-matched ZNCC, colour and gain-only residual (γ 2.2, prior weight 1), epipolar distance, the epipolar verdict (bar 3.00 px and 4.68 px), the answer and the note. *Colour wrong, ZNCC right* is a colour residual above the entry's epipolar-wrong median with blur-matched ZNCC at least 0.85; *colour right, ZNCC wrong* is a colour residual at most the epipolar-right median with blur-matched ZNCC below 0.85. No case was answered *unsure*.

**Answers by type and entry**, right / wrong:

| type | `KerryPark480` | `DnDTabletop` | both | epipolar verdict, both |
|---|---|---|---|---|
| colour wrong, ZNCC right | 7 / 3 | 6 / 4 | 13 / 7 | 8 / 12 |
| colour right, ZNCC wrong | 8 / 2 | 7 / 3 | 15 / 5 | 12 / 8 |

The reviewer sides with the ZNCC on 13 of the 20 *colour wrong, ZNCC right* cases and with the colour reading on 7; on the 20 *colour right, ZNCC wrong* cases, with the colour reading on 15 and with the ZNCC on 5. Over the 40, the colour reading on 22 and the ZNCC on 18. The epipolar verdict sides with the colour reading on 24 of the 40 (12 and 12).

**Against the epipolar verdict**, cases, epipolar verdict by answer:

| | answered right | answered wrong |
|---|---|---|
| all 40: epipolar right | 18 | 2 |
| all 40: epipolar wrong | 10 | 10 |
| colour wrong, ZNCC right: epipolar right | 8 | 0 |
| colour wrong, ZNCC right: epipolar wrong | 5 | 7 |
| colour right, ZNCC wrong: epipolar right | 10 | 2 |
| colour right, ZNCC wrong: epipolar wrong | 5 | 3 |

Agreement is 28 of 40: 15 of 20 and 13 of 20 by type. Per entry, epipolar right / wrong among the answered-right and answered-wrong cases: `KerryPark480` *colour wrong, ZNCC right* 6 / 1 and 0 / 3, *colour right, ZNCC wrong* 5 / 3 and 0 / 2; `DnDTabletop` 2 / 4 and 0 / 4, 5 / 2 and 2 / 1.

**The gain-only residual at the same split.** The gain-only residual is read against its own medians over the same members (γ 2.2): above the epipolar-wrong median (0.2203 and 0.1652) it says wrong, at most the epipolar-right median (0.1292 and 0.0837) it says right, and between them it says neither. On the *colour wrong, ZNCC right* cases it says wrong on 19 (answered right 13, wrong 6) and right on 1 (d24, answered wrong). On the *colour right, ZNCC wrong* cases it says right on 17 (answered right 13, wrong 4) and neither on 3 (d01, d05 and d25, all `DnDTabletop`, answered wrong, right and right). The gain-only and colour verdicts differ on these 4 cases and agree on the other 36.

**Ordering within a type**, ordering quality as in step 1 (a ZNCC orders a wrong member lower):

| reading | colour wrong, ZNCC right: human (13 / 7) | epipolar (8 / 12) | colour right, ZNCC wrong: human (15 / 5) | epipolar (12 / 8) |
|---|---|---|---|---|
| colour residual | 0.604 | 0.656 | 0.613 | 0.625 |
| gain-only residual | 0.505 | 0.458 | 0.667 | 0.302 |
| blur-matched ZNCC | 0.297 | 0.490 | 0.773 | 0.667 |
| plain ZNCC | 0.418 | 0.573 | 0.773 | 0.667 |

**Cases where the answer and the epipolar verdict differ:**

| case | entry | type | bin | member | cluster | blur-matched ZNCC | colour residual | gain-only residual | epipolar px | epipolar | answer | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| d01 | `DnDTabletop` | colour right, ZNCC wrong | 0.7–1.4 | 1736671 | 302139 | 0.845 | 0.104 | 0.103 | 1.130 | right | wrong | this patch crosses two surfaces with parallax, so parts match and others don't. would split left half and right half as two patches |
| d19 | `DnDTabletop` | colour right, ZNCC wrong | 0.7–1.4 | 1284220 | 186432 | 0.784 | 0.070 | 0.070 | 0.253 | right | wrong | similar repeating pattern |
| d11 | `DnDTabletop` | colour right, ZNCC wrong | 0.7–1.4 | 603363 | 70361 | 0.767 | 0.120 | 0.062 | 5.752 | wrong | right | |
| d12 | `DnDTabletop` | colour wrong, ZNCC right | 0.7–1.4 | 922497 | 118404 | 0.866 | 0.254 | 0.215 | 664.210 | wrong | right | |
| d13 | `KerryPark480` | colour right, ZNCC wrong | 0.7–1.4 | 955 | 100 | 0.846 | 0.061 | 0.047 | 6.697 | wrong | right | |
| d15 | `DnDTabletop` | colour wrong, ZNCC right | 1/4–1/2 | 1564184 | 252711 | 0.881 | 0.363 | 0.352 | 20.018 | wrong | right | |
| d21 | `DnDTabletop` | colour right, ZNCC wrong | 0.7–1.4 | 1237922 | 176814 | 0.825 | 0.079 | 0.064 | 10.873 | wrong | right | |
| d22 | `DnDTabletop` | colour wrong, ZNCC right | 0.7–1.4 | 2166201 | 470751 | 0.853 | 0.340 | 0.207 | 2271.843 | wrong | right | |
| d23 | `KerryPark480` | colour wrong, ZNCC right | 1/2–0.7 | 30281 | 9099 | 0.905 | 0.431 | 0.406 | 3.482 | wrong | right | |
| d31 | `KerryPark480` | colour right, ZNCC wrong | 0.7–1.4 | 5964 | 818 | 0.806 | 0.044 | 0.038 | 13.200 | wrong | right | |
| d32 | `DnDTabletop` | colour wrong, ZNCC right | 0.7–1.4 | 1142571 | 157480 | 0.864 | 0.229 | 0.229 | 6.927 | wrong | right | |
| d38 | `KerryPark480` | colour right, ZNCC wrong | 0.7–1.4 | 11859 | 2121 | 0.840 | 0.133 | 0.054 | 12.756 | wrong | right | |

**Notes.** Nine cases carry a note, all answered wrong; the three other wrong answers (d04, d20, d35) and every right answer carry none. Two notes name colour (d06, d37), one a surface boundary with parallax (d01), four a similar or repeating texture (d07, d19, d29, d36), one two different objects (d34), and one says the patches match well (d24). None names glare or blur.

**What it shows.** When the colour reading says wrong and the ZNCC says right, the reviewer takes the ZNCC's side on 13 of 20, and 5 of those 13 are wrong by the epipolar rule; when the colour reading says right and the ZNCC says wrong, the reviewer takes the colour reading's side on 15 of 20, and 5 of those 15 are wrong by the epipolar rule. Within each type the colour residual still orders the reviewer's wrong answers above the right ones at 0.60 to 0.61, against 0.85 to 0.88 on the step-1 cases. The gain-only residual gives the same verdict as the colour residual on 36 of 40 cases, so these cases do not separate the homography from per-channel gains. Two of the twelve wrong answers are explained by colour in the note, both on *colour wrong, ZNCC right* cases; the others name texture, surface boundaries or different objects, and none names glare or blur. Of the 12 cases where the reviewer and the epipolar verdict differ, 10 are answered right against an epipolar distance above the bar, 4 of them within 2.5 times the bar and 6 at 10 px or more.

## Open questions

- **Response from photographs alone.** Fitting a per-camera curve needs overlapping surfaces seen at different exposures; whether an SfM image set has enough of them, or whether EXIF exposure and a camera database are needed as seeds.
- **Linear or encoded storage.** Whether the record maps to a linear radiance or stops at an encoded, response-corrected space; the pairwise reading needs only the latter.
- **Where the per-observation reading goes.** Whether the colour residual joins the per-observation score columns ([sharper-patch-bitmap.md](sharper-patch-bitmap.md) Part 7), or is read fresh wherever the tiles are rendered.
- **The renderer.** Whether the sampler applies the record when it renders a tile, so downstream scores see radiance, or whether only the colour reading does.
- **The local deviation's form.** Whether glare is read as the distance between the per-tile fit and the image's record in the model's own parameters (exposure, desaturation, clipped share), or as a separate glare model with a view-direction term; and how a wrong match that fits as "glare" is bounded, since brightening with desaturation is plausible and a hue change is not.
- **Rig cameras.** Whether the two fisheyes of a rig share exposure per frame, which halves the per-image parameters on kerry_park.

## Non-goals

- Predicting photometric parameters for a view that has no photograph. Every image here is a photograph; the record is fitted, not predicted.
- Tone mapping, lens flare and other spatially varying effects.
- Changing the per-channel ZNCC. It stays the match score; the colour reading is beside it.
