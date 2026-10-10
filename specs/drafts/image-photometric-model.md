# Image Photometric Model

**Status:** Draft. Decided:
- every photograph in a reconstruction carries a photometric record: a per-image part (exposure, chromaticity) that changes shot to shot, and a per-camera part (response curve, vignetting) that a camera's frames share. The record maps the photograph's stored pixel values to a scene radiance that is comparable across images;
- the model follows PPISP's four stages and their order, applied to linear radiance: exposure as a base-2 scale, per-channel radial vignetting about a centre, a chromaticity homography built from offsets of the R, G, B and white primaries, and a per-channel piecewise power response curve with a gamma. The time-based controller that PPISP trains to predict per-frame values for novel views is not adopted: every image here is a photograph, so its values are fitted, never predicted;
- two additions PPISP does not model: clipped pixels are detected on the photograph and excluded from every fit and counted, and a one-parameter blend toward luminance covers partial desaturation (haze, veiling glare, a flat picture style);
- the first consumer is a colour agreement reading between two blur-matched patch tiles: fit the smallest per-image difference that maps the member's tile onto the reference's, and read the residual after the fit and the plausibility of the fitted parameters. Per-channel ZNCC stays the geometric match score; the colour reading is a second channel beside it;
- the global ambiguities are fixed by recentring, not by pinning one image: the mean log exposure over an image set is zero, and the mean chromaticity offset is zero;
- no bar is chosen in this draft. Every threshold is measured against the ground-truth epipolar verdicts and a blind human review in the cluster-strips convention before it enters a gate.

Not decided: whether the per-camera response and vignetting are fitted from the photographs alone or seeded from EXIF and camera metadata; the file layout of the record in `.sfmr` and `.matches`; whether the per-pair colour reading is stored per observation like the ZNCC columns; how the record is used by the renderer when it samples a view. See [Open questions](#open-questions).

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

## Open questions

- **Response from photographs alone.** Fitting a per-camera curve needs overlapping surfaces seen at different exposures; whether an SfM image set has enough of them, or whether EXIF exposure and a camera database are needed as seeds.
- **Linear or encoded storage.** Whether the record maps to a linear radiance or stops at an encoded, response-corrected space; the pairwise reading needs only the latter.
- **Where the per-observation reading goes.** Whether the colour residual joins the per-observation score columns ([sharper-patch-bitmap.md](sharper-patch-bitmap.md) Part 7), or is read fresh wherever the tiles are rendered.
- **The renderer.** Whether the sampler applies the record when it renders a tile, so downstream scores see radiance, or whether only the colour reading does.
- **Rig cameras.** Whether the two fisheyes of a rig share exposure per frame, which halves the per-image parameters on kerry_park.

## Non-goals

- Predicting photometric parameters for a view that has no photograph. Every image here is a photograph; the record is fitted, not predicted.
- Tone mapping, lens flare and other spatially varying effects.
- Changing the per-channel ZNCC. It stays the match score; the colour reading is beside it.
