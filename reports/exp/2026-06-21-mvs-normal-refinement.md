# Experiment: MVS-style normal refinement (all visible views)

_Date: 2026-06-21 · Dataset: seoul_bull (17 images, 1054 points)_

## Question

Normal refinement (`sfm xform --refine-normals`) scores a point's surfel over
its **track** — the images that matched a feature there. What happens if we
instead refine over **every image that geometrically sees the point** (projects
in front of the camera, inside the frame), turning the sparse-feature consensus
into a denser, MVS-like one? Does the normal refine more broadly, and how do the
robust (IRLS) inclusion weights distribute over the extra views?

## What was run

* New binding hook: `PatchCloud.refine_normals(view_indices=...)` overrides the
  track-based per-patch view sets with an explicit list (here: all
  geometrically-visible images per point).
* `scripts/exp_mvs_normal_strips.py` selects the points whose visible set most
  exceeds their track, refines each twice (track-only normal vs. all-visible
  "MVS" normal), and renders padded patch strips. Each tile is annotated with its
  image index (red = visible-but-not-tracked) and its IRLS inclusion weight; the
  green box marks the validated extent inside the padding.
* Output: `reports/exp/seoul_bull_mvs_strips.jpg`.

## Result

On this tight orbit, **every** rendered point geometrically projects into all 17
frames, so each gets 3 track + 14 extra = 17 views.

| pt | normal moved | consensus Φ (track-n → MVS-n) | inclusion weights (uniform = 0.059) |
|----|-------------:|-------------------------------|-------------------------------------|
| 0  | 16.7° | +0.276 → +0.281 | 0.056–0.060 |
| 1  | 34.0° | +0.147 → +0.237 | 0.058–0.060 |
| 2  | 11.1° | +0.142 → +0.150 | 0.057–0.060 |
| 3  |  0.0° | +0.450 → +0.450 | 0.057–0.060 |
| 7  | 35.5° | +0.061 → +0.110 | 0.058–0.059 |
| 10 | 31.7° | +0.107 → +0.131 | 0.057–0.060 |
| 16 | 32.5° | +0.127 → +0.075 | 0.057–0.060 |
| 18 | 25.6° | +0.389 → +0.421 | 0.054–0.060 |

**Two findings:**

1. **The normal does refine more broadly.** Refining over all 17 views moves the
   normal substantially (11–36°) and raises the consensus Φ over the broader set
   in 6 of 8 points. So even a noisy MVS view set pulls the surfel toward a plane
   that more views agree on (one point, #16, got worse — the expansion there
   degraded rather than helped).

2. **The robust inclusion weights stay essentially uniform** (spread ≤ 0.006
   around the 0.059 uniform share) — the IRLS consensus does **not** reject the
   extra views, even though many are self-occluded (a front-facing point still
   "projects into" rear cameras, where the bull's own body is in the way). The
   ZNCC robust weighting needs a dominant inlier majority to localize outliers;
   once enough occluded/disagreeing views are mixed in, the residual MAD inflates
   and the Tukey cutoff stops discriminating, so weights collapse back to uniform.

## Takeaway

Naive "all images that see the point" geometric visibility is **not** safe as an
MVS expansion here: it admits self-occluded views the soft IRLS weighting can't
cull. This is exactly the case for spec item #2 (occlusion-aware **good-view
set**) and #7 (per-pixel robustness) in
`specs/core/patch-normal-refinement.md` — a real occlusion test (depth/ray vs.
the cloud) should gate the view set _before_ the photometric consensus, rather
than relying on the consensus to reject occluders after the fact. The
`view_indices` hook makes such an occlusion-filtered view set easy to feed in as
a follow-up.

Reproduce:

```bash
pixi run python scripts/exp_mvs_normal_strips.py \
  seoul_bull_ws/sfmr/<solve>.sfmr -o reports/exp/seoul_bull_mvs_strips.jpg
```

---

# Follow-up: staged good-view refinement

Rather than gate views with a geometric occlusion test, vet them
**photometrically** — the iterative good-view set the takeaway above (and spec
item #2) calls for:

1. **Refine on the track** → a trusted seed normal + appearance template.
2. **Score the rest against it** — render every other visible view under that
   normal and correlate (windowed ZNCC) to the track template.
3. **Admit the good candidates** — accept extra views whose correlation clears a
   threshold tied to the track's own self-agreement (adapts to surfel texture).
4. **Re-refine** over track + admitted, repeating a couple of rounds.

Run with `scripts/exp_goodview_normal_strips.py` (`--accept-frac 0.7`,
`--iters 2`). Points are chosen by how many extra views they *admit* (`--rank
admits`, a quick pre-vet over a candidate pool) so the montage surfaces the
expansion working rather than the heaviest-occlusion cases. The top strip is the
track normal (each view's correlation to the track template; green = admitted,
red = rejected/occluded); the bottom strip is the final good-view normal
(admitted views carry their inclusion weight).

### seoul_bull (17 images) → `reports/exp/seoul_bull_goodview_strips.jpg`

| pt  | track | admitted | rejected | thr | normal moved | Φ over good set (track-n → good-n) |
|-----|------:|---------:|---------:|----:|-------------:|------------------------------------|
| 489 | 3 | 7 | 7 | 0.67 | 29.6° | +0.788 → +0.789 |
| 529 | 3 | 7 | 7 | 0.69 | 12.9° | +0.853 → +0.863 |
| 570 | 3 | 6 | 8 | 0.69 | 16.3° | +0.876 → +0.880 |
| 617 | 3 | 6 | 8 | 0.68 | 23.6° | +0.806 → +0.811 |
| 624 | 3 | 6 | 8 | 0.68 |  3.9° | +0.838 → +0.846 |
| 758 | 3 | 6 | 8 | 0.65 | 36.4° | +0.819 → +0.886 |
| 772 | 3 | 6 | 8 | 0.69 | 76.5° | +0.887 → +0.850 |
| 780 | 3 | 6 | 8 | 0.69 |  6.5° | +0.884 → +0.879 |

The vetting admits 6–7 of the ~14 extra views and rejects the rest, roughly
**tripling** the 3-view track support at Φ **0.79–0.89** (vs. 0.08–0.28 for the
naive all-visible set — the occluders are gone, so the score reflects an actual
surface, not an averaged-out blur). Re-refining over the good set sharpens the
normal where it helps (pt 758: +0.819 → +0.886) and leaves converged points
alone; pt 772 is the cautionary case — a 76° swing that *lowered* Φ, i.e. the
seed was good and the bigger set pulled it into a worse basin.

### seattle_backyard (26 images) → `reports/exp/seattle_backyard_goodview_strips.jpg`

| pt   | track | admitted | rejected | thr | normal moved | Φ over good set (track-n → good-n) |
|------|------:|---------:|---------:|----:|-------------:|------------------------------------|
| 168  | 3 | 22 | 0 | 0.69 |  4.1° | +0.853 → +0.848 |
| 408  | 3 | 21 | 0 | 0.68 |  8.3° | +0.878 → +0.880 |
| 1581 | 3 | 21 | 1 | 0.69 | 11.8° | +0.801 → +0.799 |
| 435  | 3 | 21 | 0 | 0.67 | 14.4° | +0.814 → +0.839 |
| 529  | 3 | 21 | 0 | 0.69 | 10.1° | +0.817 → +0.829 |
| 646  | 3 | 20 | 1 | 0.64 | 41.2° | +0.740 → +0.740 |
| 1232 | 3 | 20 | 2 | 0.66 |  0.0° | +0.744 → +0.744 |
| 275  | 3 | 21 | 0 | 0.61 | 24.0° | +0.737 → +0.669 |

The same scheme generalizes, and the dataset shows **both regimes**:

* **Broadly-visible planar surfaces** (the points above) admit **20–22 of ~24**
  extra views — the surfel really is seen and consistent everywhere — turning a
  3-view track into a ~24-view consensus at Φ 0.74–0.88. This is the case naive
  MVS would *also* get right; here the vetting just confirms it.
* **Occluded / cluttered points** (what `--rank extra` surfaces instead — points
  that geometrically project into ~25 frames but are mostly blocked) get **20–22
  of those views rejected**, the vetting keeping only the few that agree. This is
  the case naive MVS gets *wrong*, and where the staged scheme earns its keep.

So the admit count is itself a useful signal: high = a confidently-visible
surfel, low = an occlusion-bound point where only the track can be trusted.

### kerry_park (24 frames × 2 fisheyes) → `reports/exp/kerry_park_goodview_strips.jpg`

A 360° two-fisheye rig (`OPENCV_FISHEYE`), global solve (GLOMAP), 48 images.

| pt  | track | admitted | rejected | thr | normal moved | Φ over good set (track-n → good-n) |
|-----|------:|---------:|---------:|----:|-------------:|------------------------------------|
| 362 | 3 | 17 | 4 | 0.70 | 26.9° | +0.964 → +0.964 |
| 365 | 3 | 18 | 3 | 0.69 |  0.0° | +0.947 → +0.947 |
| 94  | 3 | 17 | 5 | 0.69 |  0.0° | +0.845 → +0.845 |
| 135 | 3 | 17 | 4 | 0.68 | 58.1° | +0.836 → +0.868 |
| 152 | 3 | 17 | 4 | 0.69 | 14.3° | +0.868 → +0.865 |
| 164 | 3 | 17 | 4 | 0.65 | 47.7° | +0.789 → +0.808 |
| 245 | 3 | 17 | 4 | 0.68 | 33.8° | +0.829 → +0.853 |
| 293 | 3 | 17 | 4 | 0.68 |  0.0° | +0.884 → +0.884 |

The scheme is camera-model-agnostic — the patch warp goes through the fisheye
`ray_to_pixel`, so the same vetting runs unchanged on the rig. These textured
urban facades give the **highest Φ of the three datasets (0.79–0.96)**: a 3-view
track grows to ~20 views, with only 3–5 of the visible extras rejected. Several
points take large normal swings that *raise* Φ (pt 135: 58° → +0.032; pt 164:
48° → +0.019) — the wide rig baseline gives the re-refinement real leverage to
correct a track normal once it has ~20 views to agree against.

### dino_dog_toy (85 images) → `reports/exp/dino_dog_toy_goodview_strips.jpg`

A dense 85-image capture of a small non-convex toy, incremental solve, 19 024
points. Points project geometrically into ~82–85 frames, so the strips are
subsampled to 30 tiles for legibility (`--max-views 30`; decisions/Φ are still
over the full set).

| pt    | track | admitted | rejected | thr | normal moved | Φ over good set (track-n → good-n) |
|-------|------:|---------:|---------:|----:|-------------:|------------------------------------|
| 10348 | 3 | 47 | 35 | 0.54 |  8.2° | +0.590 → +0.599 |
| 13611 | 3 | 41 | 41 | 0.67 | 14.4° | +0.823 → +0.809 |
| 14201 | 3 | 50 | 32 | 0.61 | 17.1° | +0.593 → +0.636 |
| 16974 | 3 | 41 | 41 | 0.68 |  2.6° | +0.810 → +0.808 |
| 13871 | 3 | 38 | 44 | 0.68 |  7.4° | +0.789 → +0.795 |
| 14214 | 3 | 40 | 42 | 0.66 |  6.7° | +0.744 → +0.752 |
| 14259 | 3 | 39 | 43 | 0.69 | 11.3° | +0.827 → +0.836 |

The hardest case, and the most telling: a small **non-convex** object shot from
all around means a surface point is geometrically in-frame almost everywhere but
actually unoccluded/front-facing only about **half** the time — so the vetting
admits ~40–50 and rejects ~32–44, a near 50/50 split, exactly the discrimination
naive all-visible cannot make. It still grows the 3-view track to ~40+ vetted
views at Φ 0.59–0.83, and re-refinement nudges Φ up on most points.

## Takeaway (staged)

Photometric vetting against a track-seeded template is a robust,
dataset-and-camera-agnostic way to grow the refinement view set: across four
datasets it tripled support on the sculpture orbit, ~8×'d it on broadly-visible
backyard surfaces, ~7×'d it on the fisheye rig (Φ up to 0.96), and on the dense
non-convex toy split the ~82 in-frame views almost 50/50 into the ~40 that
actually see the surface vs. the occluded rest — correctly refusing to expand
where it shouldn't. The remaining gap is seed dependence — a bad track
normal vets poorly, and even a good one can be pulled into a worse basin by the
larger set (seoul_bull pt 772). A geometric occlusion pre-filter (depth vs. the
cloud) to bootstrap, keeping the re-refinement only when Φ improves, and the
seed-free per-pixel robustness of spec item #7 are the natural next steps; the
`view_indices` hook already carries whatever view set any of them produce.

## How much of the expansion is non-feature coverage?

If the admitted views were just views the matcher *dropped*, the expansion would
add little a better matcher couldn't. So for each admitted view we projected the
point and measured the distance to the nearest detected SIFT keypoint
(`--sift-overlap`). The track views — which *are* keypoints — sit 0.2–0.7px from
their keypoint, so **2px is the "a feature is really here" radius**; 10px mostly
just measures feature density (in textured frames there's almost always *some*
keypoint within 10px of any pixel), so it is not a correspondence test.

| dataset | admitted views | median dist | ≤2px (feature at point) | ≤5px | ≤10px |
|---------|---------------:|------------:|------------------------:|-----:|------:|
| seoul_bull       |  50 | 3.1px | **28%** | 80% | 100% |
| seattle_backyard | 167 | 3.0px | **34%** | 76% |  99% |
| kerry_park       | 137 | 2.8px | **47%** | 84% |  99% |
| dino_dog_toy     | 346 | 5.8px | **13%** | 42% |  75% |
| _track baseline_ |  —  | 0.2–0.7px | 92–100% | — | — |

**Only ~13–47% of the added views had a SIFT keypoint at the point** (the strict
2px radius). So the majority — roughly half on the orbit/rig, ~two-thirds on the
backyard, and **~⅞ on the dense toy** — is genuinely feature-less coverage: views
where SIFT either didn't fire at that spot or wouldn't match it, yet the patch is
still photoconsistent. The ~30–47% within 2px are the recoverable "matcher
misses" (a feature was there, just not triangulated into this track); the rest is
true MVS gain. dino makes the point sharpest — its admitted views are mostly
oblique/foreshortened looks that SIFT can't repeat, so the expansion there is
almost entirely non-feature. The staged scheme is doing real MVS densification,
not just re-collecting dropped features.

## Did the normal actually get more robust?

More agreeing views *should* mean a better-pinned normal — but does it? Assessed
three ways with `scripts/exp_goodview_render_patches.py`, refining the same
eligible points track-only vs. good-view:

* **Φ-peakedness confidence** — the curvature of the photoconsistency at the
  optimum (how sharply the normal is determined).
* **Leave-one-out (LOO) stability** — re-refine dropping each view; the angular
  spread of the results is the normal's sensitivity to any single view.
* **Image overlays** — `sfm render-patches --mode texture/normal` on a source
  frame; a correct surfel's texture continues the image.

| dataset | views track→good | confidence track→good | LOO spread (good) | normal moved |
|---------|-----------------:|----------------------:|------------------:|-------------:|
| seoul_bull       | 3 → 6.0  | 0.092 → 0.136 (+48%) |  8.4° | 7.3° med |
| seattle_backyard | 3 → 11.6 | 0.086 → 0.093 (+8%)  | 12.7° | 8.7° med |
| kerry_park       | 3 → 12.3 | 0.058 → 0.058 (~0%)  |  5.0° | 18.1° med |

**The honest answer: redundancy improved a lot; the normal's _determinacy_ did
not.** What clearly got better is support — every point goes from a 3-view track
with *zero* redundancy (drop one view and it falls below `min_views` — every view
is load-bearing) to 6–12 vetted views, and the consensus Φ *value* rises to
0.7–0.9 (vs. 0.08–0.28 for naive all-visible). But the normal *direction* stays
soft: the Φ peak is shallow everywhere (confidence 0.06–0.14, barely moved on
seattle/kerry), and dropping a single view still swings the good-view normal
5–13°. The normals also *move* 7–18° from the track estimate, which without
ground truth we can't call "more correct," only "differently posed, backed by
more agreeing views." Robustness is **surface-dependent, not view-count-driven**:
kerry's sharp urban facades are the most stable (LOO 5°) despite the largest
moves, while seattle's flatter foliage is the least (12.7°) — these low-relief,
textured surfels are intrinsically weakly observable in the 2-DOF normal, and
adding views raises Φ and kills single-view fragility without sharpening the peak.

The image overlays agree: the good-view textures align about as well as track and
the normal field is only modestly smoother — consistent with a real
support/consensus gain rather than a step change in normal accuracy. Pinning the
normal harder needs a stronger signal than photoconsistency alone (anisotropic
sampling for unbiased Φ on oblique views, a cloud-smoothness prior, or the
per-pixel robust template of spec item #7), not just more views.

### High-fidelity check (anisotropic + exact, no fronto cache)

The shallow confidence above was partly a *measurement* artifact: bilinear
sampling depresses Φ on oblique views, which flattens the apparent peak. Re-ran
seattle_backyard's same 3266 points with `--sampler anisotropic --cache off` (the
exact path, ~12 min vs. seconds):

| metric | bilinear + fronto | anisotropic + exact |
|--------|------------------:|--------------------:|
| confidence — track     | 0.086 | **0.242** |
| confidence — good-view | 0.093 | **0.203** |
| good vs. track         | +8%   | **−16%**  |
| LOO spread (good)      | 12.7° | 10.9°     |
| normal moved track→good| 8.7°  | 11.0°     |

Two things change the story, one doesn't:

* **The default sampler was under-reporting peak conditioning ~2.5×** (confidence
  0.086 → 0.242). With honest Φ the normals are far better-conditioned than the
  bilinear metric implied — so part of the earlier "shallow peak" was the
  measurement, not the geometry.
* **More views still doesn't sharpen the normal** — and now *reverses*: under
  unbiased Φ the 3-view track has a slightly sharper peak than the 11-view good
  set (0.242 vs 0.203). Adding diverse, oblique views broadens the optimum even
  as it raises the consensus value; peakedness is not monotone in view count.
* **LOO stability barely moves** (12.7° → 10.9°): the normal is still soft to a
  view dropout regardless of sampler.

So anisotropic+exact is worth it for an *honest* Φ/confidence read-out, but it
does not turn the good-view expansion into a normal-sharpening step. The
expansion's value remains redundancy and consensus support, not determinacy.

### Patch extent 1× vs 5× (anisotropic + exact)

A smaller surfel is more local — less surface curvature averaged into the patch —
so it *might* sharpen the normal. Re-ran seattle at `--extent 1` (vs. the default
`5`), same anisotropic+exact path:

| metric | extent 5 | extent 1 |
|--------|---------:|---------:|
| good-view views        | 11.6 | 9.8 |
| confidence track→good  | 0.242 → 0.203 | 0.240 → 0.214 |
| LOO spread (good)      | 10.9° | **16.3°** |
| normal moved (med / 90th) | 11.0° / 45.7° | 14.9° / **64.5°** |

Smaller is **worse**, not sharper: confidence is unchanged (≈0.21–0.24), but the
normal gets *noisier* — LOO stability degrades 10.9° → 16.3° and 10% of points
swing >64°. A 1×-feature patch carries too little texture to constrain the normal,
so each view's estimate is noisy and the consensus wanders; the 5× patch wins by
integrating more signal (the curvature it averages in costs less than the texture
it gains). The overlay shows it directly
(`reports/exp/seattle_backyard_extent_compare.jpg`): extent 1 covers far less of
the frame and its normal field is visibly fragmented vs. extent 5's coherent
surfels. So on these surfaces the determinacy is bounded by patch *signal*, not by
locality — shrinking the patch trades coverage and stability for nothing.

Reproduce:

```bash
pixi run python scripts/exp_goodview_normal_strips.py \
  seoul_bull_ws/sfmr/<solve>.sfmr -o reports/exp/seoul_bull_goodview_strips.jpg
pixi run python scripts/exp_goodview_normal_strips.py \
  seattle_backyard_ws/sfmr/<solve>.sfmr -o reports/exp/seattle_backyard_goodview_strips.jpg
pixi run python scripts/exp_goodview_normal_strips.py \
  kerry_park_ws/sfmr/<solve>.sfmr -o reports/exp/kerry_park_goodview_strips.jpg
pixi run python scripts/exp_goodview_normal_strips.py \
  dino_dog_toy_ws/sfmr/<solve>.sfmr -o reports/exp/dino_dog_toy_goodview_strips.jpg \
  --max-views 30

# add --sift-overlap to any of the above for the keypoint-overlap stats

# robustness assessment + reconstructions for render-patches overlays
pixi run python scripts/exp_goodview_render_patches.py \
  seattle_backyard_ws/sfmr/<solve>.sfmr -o /tmp/seattle
sfm render-patches /tmp/seattle.good.sfmr -o renders/ --mode texture --opaque --images _12
```

Overlay comparison (track vs. good-view, texture + normal field):
`reports/exp/seattle_backyard_overlay_compare.jpg`.

---

# Is the soft normal a *registration* problem? (sub-pixel view localization)

Determinacy never sharpened with more views, and a natural suspect is residual
**mis-registration**: refinement assumes each view's projection lands exactly on
the surfel (the warp comes straight from the point's 3D position + pose). If the
point center is a little off, or a pose/distortion residual nudges a view, that
view's rendered patch is *translated* a fraction of a pixel relative to the
consensus — which blurs the template and flattens the Φ peak even when the normal
is right. `scripts/exp_view_localization.py` measures that translation directly.

For each view of a well-supported (good-view) point it:

1. Renders the surfel onto a **larger context tile** than the scored 32px core,
   so a search window can slide without running off the patch edge (the boundary
   problem the obvious fixed-size render has).
2. Builds a **leave-one-out** robust consensus reference from the *other* views
   (so a view is never aligned to a template its own pixels polluted).
3. **Coarse-to-fine, sub-pixel translation search** — a half-res pyramid pass for
   range/robustness (the "downsampling"), then a full-res integer refine + a
   separable parabolic sub-pixel fit — for the shift that maximizes windowed ZNCC
   to that reference.

It then decomposes the per-view shifts into a **common** in-plane component (the
point *center* is off — one re-centering would fix every view) vs. the
**residual** incoherent spread (per-view pose/distortion/depth error no single
re-centering can remove), and re-renders at the found offset (translating the
patch center — a *single* resample) to verify the gain isn't an artifact of
resampling the context tile.

| dataset (8 pts) | median shift (src px) | ZNCC gain | ΔΦ centered→aligned | common vs. residual | exact-confirmed |
|-----------------|----------------------:|----------:|--------------------:|---------------------|----------------:|
| seoul_bull       | 0.55 | +0.010 | +0.046 | **0.18 vs 0.92 px** | 44% |
| seattle_backyard | 0.31 | +0.005 | +0.055 | **0.09 vs 0.58 px** | 90% |
| kerry_park (fisheye) | 0.32 | +0.009 | +0.053 | **0.05 vs 0.48 px** | 74% |
| dino_dog_toy (non-convex) | 1.96 | +0.012 | +0.073 | **1.25 vs 5.98 px** | 37% |

**Findings:**

1. **On well-behaved surfaces the projections are already sub-pixel localized.**
   Median best-fit shift is 0.3–0.55 source px on the orbit / backyard / fisheye —
   the warp puts each view within half a pixel of the consensus, so there is no
   gross mis-registration to recover. **dino_dog_toy breaks this** (median 1.96 px,
   ΔΦ +0.073): the dense non-convex toy is the one scene where the views are *not*
   well registered — see the dino note below.

2. **Aligning barely helps.** Snapping every view to its leave-one-out template
   raises ZNCC by only **+0.005–0.010** and the consensus Φ by **~+0.05** — real
   but small, and *not* the missing factor behind the soft normal. This
   corroborates the earlier conclusion: determinacy is **signal/surface-bound**
   (patch texture and observability), not registration-bound.

3. **What residual there is (on the good scenes), is incoherent, not a center
   offset.** The shifts decompose into a tiny **common** component (0.05–0.18 px —
   re-centering the 3D point would buy almost nothing) and a much larger
   **residual** incoherent spread (0.48–0.92 px). So the per-view jitter is
   pose/distortion/depth noise sprinkled across views, not a systematic "the point
   is in the wrong place" error a re-triangulation could fix globally.

4. **The probe doubles as point QC.** The few points with a large *common* shift
   are the genuinely mis-located ones — kerry_park pt 187 (common 2.83 / residual
   10.85 px) and seattle pt 275 (common 0.50 px) stand out from the sub-pixel
   crowd, exactly the points whose center triangulation is off.

5. **The context + downsample search is sound where it matters.** On the
   redundant well-registered datasets the exact single-resample re-render confirms
   the context-tile ZNCC **74–90%** of the time; the low seoul_bull figure (44%) is
   the small-gain regime where the ±0.03 confirm slack is comparable to the ~0.01
   signal itself, not a double-resampling failure. The larger context did remove
   the boundary issue — windows never clip — and the half-res coarse pass found the
   same optima as full-res at a quarter of the candidate count.

### dino_dog_toy — the one scene that is genuinely mis-registered

The dense 85-image non-convex toy behaves nothing like the other three: median
shift **1.96 px** (4–6× the others), with several points far worse (pt 13611:
5.8 px shift, common 3.6 / residual 8.1; pt 13871: residual **15.1 px**). Both
the **common** (1.25 px) and **residual** (5.98 px) channels blow up, and the
exact re-render confirms only **37%** of the found shifts. That low confirm rate
is the diagnosis, not a search bug: on a non-convex object ~half the geometrically
in-frame views are occluded or grazing (the good-view experiment already saw the
~50/50 admit/reject), so the leave-one-out reference is itself blurred by bad
views and the translation search latches onto spurious correlation peaks the exact
re-render won't reproduce. The non-trivial common component (1.25 px, vs ≤0.18 px
elsewhere) says some of these points are also genuinely mis-*centered* — small
foreshortened triangulations on a cluttered surface. So dino is the case where
re-localization *would* move things, but it can't be trusted view-by-view until
the view set is occlusion-clean first; it argues for the geometric occlusion
pre-filter (and re-centering as a per-point QC gate) before any photometric shift.

**Takeaway.** Sub-pixel re-localization is a *dead end for sharpening the normal*
on clean scenes — the views are already registered to ~0.4 px and the recoverable
Φ is ~0.05 — and on the one dirty scene (dino) the shifts are large but
**untrustworthy** (37% exact-confirm) because the view set is occlusion-contaminated.
Its dependable value is diagnostic: the common-shift channel cleanly flags
mis-triangulated points, and the exact-confirm rate flags occlusion-contaminated
view sets. Pinning the normal harder still needs a stronger signal than photometric
translation (the anisotropic-Φ readout, a cloud-smoothness prior, or spec item #7's
per-pixel robustness) plus an occlusion pre-filter, not better per-view alignment.
Montages (per-view context tile, white = scored core, cyan = found shift):
`reports/exp/{seoul_bull,seattle,kerry,dino}_localization.jpg`.

Reproduce:

```bash
pixi run python scripts/exp_view_localization.py \
  seattle_backyard_ws/sfmr/<solve>.sfmr -o reports/exp/seattle_localization.jpg
# --context N  larger search tile · --search R  max shift px · --no-coarse  full-res only
```

---

# Follow-up: feed the shifts back — congealed reference patches

The localization probe found a measurable per-view residual (sub-pixel on clean
scenes, several px on dino). The natural next step the user asked for: *use* the
shifts. Render each view at its best-fit in-plane offset and rebuild the
consensus from those **registered** renders, iterating to convergence — group-wise
translation registration, a.k.a. **congealing**. `scripts/exp_reference_refine.py`.

Two methodological guards:

* **Single-resample renders.** Each iteration adds the view's residual shift to an
  *accumulated* offset and re-renders the patch *from the source* with the patch
  center translated by the new total — never re-sampling an already-warped tile,
  so applying shifts cannot compound blur.
* **Leave-one-out scoring.** The "did it actually sharpen?" metric is mean
  per-view ZNCC against the robust consensus of the **other** views. A mean
  fitting its own noise would inflate self-agreement but cannot lift LOO. This
  is the honest signal; Φ and consensus-image gradient sharpness are reported
  alongside.

| dataset (8 pts, 5 iters) | LOO ZNCC before→after | ΔLOO | ΔΦ | consensus sharpness |
|---|---|--:|--:|--:|
| seoul_bull       | +0.867 → +0.905 | **+0.021** | +0.034 | **×1.07** |
| seattle_backyard | +0.837 → +0.869 | **+0.031** | +0.048 | **×1.12** |
| kerry_park       | +0.892 → +0.921 | **+0.028** | +0.048 | **×1.17** |
| dino_dog_toy     | +0.642 → +0.669 | **+0.022** | +0.034 | **×1.09** |

**This works.** The honest LOO ZNCC rises by **+0.02–0.03** on every dataset —
2–3× the single-step localization gain (which only aligned to a fixed reference)
because the iteration re-tunes the reference as views congeal in. ΔΦ is +0.03–0.05
and the consensus *image itself* is **7–17% sharper** by gradient energy. Visual
inspection (the montage panels) is the clearest evidence: on seattle pt 275 the
before is a brown/green mush, the after has a crisp bright lobe and a clean
boundary — sharpness ×1.63. seoul pt 88 and kerry pt 68 show the same kind of
edge crispening.

**dino is the rough case.** Median dLOO is still positive (+0.022) but with two
clear failure modes: large-shift duplicates (pts 14201/14202 drift 19 px and LOO
*drops* −0.007, sharpness ×0.80 — they congealed off-surface) and the very
weak-signal point (14259, baseline LOO 0.255) where the iteration is destabilizing
rather than refining. The clean diagnostic is **shift magnitude × baseline LOO**:
points with a multi-pixel shift *and* a low baseline are the ones to gate out.

**What this changes.** Congealing puts a sharper reference template under
everything that consumes the consensus patch — the photometric vetting in
`vet()`, the per-pixel robust template of spec item #7, and ultimately the normal
search itself. Even though re-localization alone could not sharpen the normal,
iterating the alignment *into* the reference produces a measurably better
reference for the next round — the natural input to a re-refinement of the normal
on the cleaned stack. Worth wiring into the normal-refine loop with a per-point
gate (skip when shift > ~3 px or baseline LOO < ~0.5) for the dino-like cases.

Montages (consensus patch before vs after, per point):
`reports/exp/{seoul_bull,seattle,kerry,dino}_refref.jpg`.

Reproduce:

```bash
pixi run python scripts/exp_reference_refine.py \
  seattle_backyard_ws/sfmr/<solve>.sfmr -o reports/exp/seattle_refref.jpg
# --iters N  congealing rounds (default 5) · --search R  max total shift px
```

---

# Validation: prototype vs. shipped production kernels

_Added 2026-06-25, after the experiments above landed in production._

The two mechanisms these experiments validated — staged photometric **view
selection** and **congealing** keypoint localization — were since productionized
as Rust kernels (`select_views` / `select_patch_views`,
`localize_keypoints` / `localize_patch_keypoints`) and wired into the
`sift_files → embedded_patches` pipeline (`sfm embed-patches`). The prototype
did the vetting and congealing in Python on top of the one Rust hook
(`refine_normals(view_indices=…)`); the shipped code reimplements both in the
core.

To confirm the port is faithful, two head-to-head harnesses run the Python
prototype and the Rust production path on the **same points, same (track-)refined
normals, and the same admitted view set**, and diff the decisions:

* `scripts/cmp_view_selection.py` — prototype `vet()` vs. `select_views`.
* `scripts/cmp_keypoint_localization.py` — prototype `congeal()` vs.
  `localize_keypoints`.

## View selection — a faithful 1:1 port

| dataset | mean admitted-extra (proto → prod) | mean Jaccard of admitted sets |
|---------|-----------------------------------:|------------------------------:|
| seoul_bull (12 pts)       | 2.4 → 2.4   | **1.00** |
| seattle_backyard (12 pts) | 3.1 → 3.0   | **0.83** |
| dino_dog_toy (8 pts)      | 19.5 → 19.9 | **0.87** |

The production kernel admits the same views as the Python prototype: identical on
the tight seoul orbit, and ~0.83–0.87 set-overlap with matching aggregate counts
elsewhere. The few disagreements are **single borderline views** near the admit
threshold — attributable to `robust_iters` 3 (prod) vs. 5 (proto) and the
production **trust gate** (`min_self_agreement = 0.3`) vs. the prototype's
absolute corr floor (0.1). The adaptive bar (`0.7 × self-agreement`) and the
self-agreement values themselves match to ±0.01. Even on dino's ~82-candidate
sets the counts agree (19.5 vs. 19.9); the lower Jaccard there is just more
borderline views to flip, not a behavioural difference.

## Keypoint localization — the port, plus the gate the prototype only theorized

The prototype `congeal()` registers *all* given views. Production
`localize_keypoints` ports that loop **and** adds the in-loop view dropping the
"Follow-up: congealed reference patches" takeaway above called for
(`max_shift_px = 3.0`, drop on low leave-one-out ZNCC) — "worth wiring in with a
per-point gate (skip when shift > ~3 px or baseline LOO < ~0.5) for the
dino-like cases."

| dataset | median shift, src px (proto → prod) | mean final LOO ZNCC (proto → prod) |
|---------|------------------------------------:|-----------------------------------:|
| seoul_bull (10 pts)       | 0.72 → 0.49 | 0.837 → 0.872 |
| seattle_backyard (10 pts) | 0.40 → 0.35 | 0.831 → 0.874 |
| dino_dog_toy (8 pts)      | **4.76 → 1.55** | **0.737 → 0.901** |

On the clean scenes the two agree (both find the sub-pixel shifts measured
earlier) and production's dropping nudges final LOO up ~+0.04. **dino is where it
matters.** The prototype, congealing every view including the ~half that are
occluded/grazing, runs away — median 4.76 px with individual points drifting
9–17 px (congealed off-surface), and one weak-signal point (5095) where
congealing *lowers* LOO (0.433 → 0.388), both failure modes flagged in the dino
note above. Production drops those views in-loop, collapsing the shifts to
~1–1.8 px and lifting final LOO to **0.901 vs. 0.737** — a +0.16 gain on dino
vs. +0.04 on the clean scenes:

| pid | views | proto shift | prod shift | proto LOO 0→N | prod LOO N | prod kept |
|-----|------:|------------:|-----------:|---------------|-----------:|----------:|
| 1403 | 34 | 17.62 px | 1.73 px | 0.707 → 0.737 | **0.882** | 17/34 |
| 4332 | 34 | 10.15 px | 1.03 px | 0.677 → 0.783 | **0.931** | 21/34 |
| 4056 |  7 |  9.30 px | 1.82 px | 0.775 → 0.888 | **0.942** |  4/7  |
| 5095 | 36 |  6.39 px | 1.28 px | 0.433 → 0.388 | **0.929** | 34/36 |

## Takeaway (validation)

The shipped kernels reproduce the prototype's findings on the clean scenes and
**improve on them where the prototype was weakest**: view selection is a faithful
port (the borderline-view differences are tuning, not behaviour), and keypoint
localization adds the occlusion-robust view dropping that turns dino from a
runaway (median ~5 px, LOO regressions) into a clean ~1.5 px / 0.90-LOO result.
The one open caveat from the production map still stands — the reference *bitmap*
is rendered during the track-only normal refine, before selection and congealing,
so it does not yet benefit from the congealed, view-pruned stack these
comparisons show is materially better on hard scenes.

Reproduce:

```bash
pixi run python scripts/cmp_view_selection.py \
  dino_dog_toy_ws/sfmr/<solve>.sfmr -n 8
pixi run python scripts/cmp_keypoint_localization.py \
  dino_dog_toy_ws/sfmr/<solve>.sfmr -n 8
```

---

# Does the normal refinement want better keypoint positions?

_Added 2026-06-25._

Normal refinement positions every view's patch at the **reprojection of the
shared 3D point center** (`WarpMap::from_patch` with one `OrientedPatch` per
candidate — `crates/sfmtool-core/src/patch/normal_refine.rs`). It never sees a
per-view keypoint: not the congealed sub-pixel offsets `localize_keypoints`
finds, nor the original SIFT keypoints the solve triangulated from. Two probes
ask whether either would help, each running the **same** coarse-to-fine normal
search twice and changing only where the views are positioned. (The search is a
Python proxy of the Rust kernel — it omits the frozen-support masking, so its
argmax carries a few degrees of noise; the robust signal throughout is the
**paired ΔLOO**, measured at each method's own optimum, not the per-point move.)

## A. Congealed keypoints fed back into the normal

Render each view at its keypoint-congealed in-plane offset (over the
production-kept stack) instead of the raw projection, then re-search the normal.

| scene | median kp shift | median normal move | mean ΔLOO (cong − raw) |
|-------|----------------:|-------------------:|-----------------------:|
| seoul_bull (clean)     | 0.66 px | 4.8° | **−0.012** |
| dino_dog_toy (hard)    | 2.10 px | 7.7° | **+0.026** |

Clean scenes: no benefit (ΔLOO ≈ 0, slightly negative) — the normal "moves" only
because its Φ-peak is shallow, sliding a weakly-determined optimum without
sharpening it. dino: a modest **net** win (+0.026) but **mixed** — points like
4111 (+0.099), 4056 (+0.057), 2240 (+0.052) gain, while 1403 (−0.038) and 1572
(−0.032) lose even over the pruned stack. So congealing-into-the-normal is not a
reliable, uniform improvement.

## B. SIFT-solve keypoints at the initial refine

The solve already carries the **actual detected keypoint** for every track
observation; the reprojection differs from it by the triangulation/pose
residual. Position each track view's patch on its SIFT keypoint (first-order
in-plane Jacobian solve) and re-search the normal. **The benefit scales directly
with the solve's reprojection error:**

| scene | median reproj residual | median normal move | mean ΔLOO (sift − proj) |
|-------|-----------------------:|-------------------:|------------------------:|
| seattle_backyard | 0.16 px | 13.9° (noise) | **−0.002** |
| seoul_bull       | 0.30 px |  2.9°          | **−0.000** |
| dino_dog_toy     | 0.87 px |  1.4°          | **+0.039** |

On the tight solves there is nothing to fix — the reprojection already lands
within ~0.2–0.3 px of the detected feature. On dino (3× the residual) the gain is
real and sizable: pt 15 LOO 0.830 → **0.948**, pt 17 0.854 → 0.934, pt 7
0.764 → 0.818, only one point flat. Crucially the normal **direction barely moves**
(1.4°) — the improvement is in **consensus/template quality**, not where the
normal points.

## Takeaway

Both probes converge on the same conclusion the determinacy experiments reached:
the normal **direction** is signal/observability-bound — better registration
barely moves it. What better keypoint positions *do* improve is the **registered
consensus** (LOO/Φ), i.e. the reference template and bitmap the patch carries
downstream.

Of the two, **SIFT-keypoint seeding of the initial refine is the cleaner win**:
it is free (the keypoints already exist from the solve — the sift→embedded
conversion already maps every track observation to its feature), safe (neutral on
good solves, monotonic in reprojection error, no occluder contamination since the
track features *are* the matched correspondences), and positive exactly where the
solve is loose. Congealed-keypoint feedback is fragile by comparison. Both want
the same small enabling change — a per-view keypoint/offset input to
`refine_patch_normal`'s warp — but B drives it from existing solve data rather
than a congealing pass. The highest-value consumer of the cleaner stack remains
the reference **bitmap** (rendered today over the raw track projection).

Reproduce:

```bash
pixi run python scripts/exp_normal_from_congealed_keypoints.py \
  dino_dog_toy_ws/sfmr/<solve>.sfmr 6
pixi run python scripts/exp_normal_from_sift_keypoints.py \
  dino_dog_toy_ws/sfmr/<solve>.sfmr 8
```
