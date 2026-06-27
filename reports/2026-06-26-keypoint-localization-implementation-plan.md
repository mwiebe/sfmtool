# Implementation Plan: Keypoint-Localization Search Cache + Subpixel Refinement

_Date: 2026-06-26 (revised 2026-06-27). An incremental plan for the two design
specs that landed on `main` in #135:_

- _**Grid search (Spec A):**
  `specs/core/keypoint-localization-search-cache.md` — the integer cross-view
  congealing search, accelerated by a render-once per-view cache, a hand-rolled
  AVX2 windowed-ZNCC kernel (centered-f32, then integer-`i16`), integer-tracked
  reads, and an `f32` search-resolution multiplier `m`._
- _**Subpixel refinement (Spec B):**
  `specs/core/keypoint-subpixel-refinement.md` — a standalone, continuous
  photometric (ECC / forward-additive Gauss–Newton) local refiner; the
  high-accuracy / ground-truth-approximating option._

_The prototype scripts on this branch (`scripts/exp_view_localization.py`,
`exp_normal_from_congealed_keypoints.py`, `exp_reference_refine.py`,
`_normal_strip_lib.py`) are the validated references these production
implementations should reproduce, and serve as oracles in the validation steps._

> _**Revision (2026-06-27): optimize both before measuring the production-path
> fork.** The original plan put the supersampled-grid-vs-LK decision gate
> **before** Phase 3 optimizations. That was wrong: measuring there would compare
> two unfinished implementations against each other — accuracy and speed are both
> still moving targets, and the gate's actual question, "where is the discrete
> grid good enough?", is a per-µs question that only optimized code can answer
> honestly. So **optimize first, decide last**: build A's hot path (AVX2 / `i16`)
> and B's hot path (analytic sampler Jacobian, optional incremental consensus)
> independently and to their natural completion, then run the head-to-head on
> production-shaped code and let the production-path decision ride on real data._

## Sequencing: A first, interleave a minimal B, optimize both, then measure

**Recommendation: implement Spec A's *structure* first, interleave a *minimal*
Spec B early so both algorithms exist end-to-end, then optimize each
independently, and finally run the head-to-head measurement on the optimized
implementations to decide the production path.** Not pure-sequential, not
free-for-all interleave.

Why **A first**:

- **A is the producer.** B consumes A's output — the seed keypoint, the kept view
  set, and the consensus. There is little to build or test in B end-to-end until
  A's localization pipeline produces those.
- **A is partly landed.** The scalar correlation-accumulation `search_shift`
  (branch `optimize-search-shift`, ~2.5× kernel / ~1.9× embed, validated by
  `search_shift_matches_reference*`) is the foundation the cache and AVX2 build on.
- **A carries the larger, already-measured win.** On dino, localization is ~83% of
  `embed-patches`; within it `search_shift` ~56% and `render_context` ~33%. A's
  cache + AVX2 attacks both.

Why **interleave a minimal B early** (rather than finishing all of A first):

- The specs flag a real **architectural fork**: a *supersampled* grid (`m > 1`,
  Spec A) resolves sub-pixel directly and may be good enough, while the continuous
  LK refiner (Spec B) is the high-accuracy reference. Both need to exist end-to-end
  to be optimized and then measured.
- Standing up the minimal B early also surfaces any reuse / interface frictions
  with A's machinery cheaply, before those choices get baked into the optimized
  code.

Why **optimize both before the gate** (the 2026-06-27 revision):

- The gate's question — "when is the discrete grid good enough?" — is decided in
  per-µs terms: accuracy per unit time. Measuring on unoptimized code compares
  two implementations' unfinished states, not the algorithms.
- A's `m²`-growing supersampled cost and B's per-step FD render cost are exactly
  the things the optimizations remove; the un-optimized speed numbers would
  systematically penalize whichever side has the more expensive MVP shortcut.
- Both optimization paths are **independently justified**: A's AVX2 kernel is the
  natural follow-on to the cache regardless of B, and B's analytic sampler
  Jacobian removes the per-step 5× FD render cost that any production B usage
  would carry. Neither is speculative.

Why **not** pure interleave from step one: B's MVP still needs A's seed/consensus
to run on real data, so a minimal A precedes even the minimal B.

## Phases

### Phase 0 — Land the scalar foundation (done, merged as #136)

The scalar correlation-accumulation `search_shift` (was on `optimize-search-shift`).
Validated by `search_shift_matches_reference*`; merged into `main` as #136.

### Phase 1 — Per-view cache + integer-tracked reads (done, merged as #137)

The render-once structural win, on the scalar kernel (no AVX2 yet). Per-view
expanded frame-oriented cache (`R + 4·search`); the round loop reads cores and
search candidates from it; offset tracking is split into the integer read
accumulator `iacc` (only thing indexing the cache) and a parabolic sub-pixel
residual (convergence + final keypoint only, never fed back into reads);
`search_resolution_multiplier: f32` knob (default `1.0`) added. Validated:
21/21 Rust tests including strengthened `search_shift_matches_reference*` and two
multiplier tests; dino re-profile confirmed `render_context` collapses from
per-(view, round) to per-view (renders 1.44M → 399k). Merged as #137.

### Phase 2 — Minimal subpixel reference (done, on `keypoint-subpixel-refine`)

The simplest correct continuous refiner — forward-additive ECC Gauss–Newton,
per-view against a single frozen IRLS consensus, bilinear source-pyramid
sampling, analytic z-norm derivative with the raw image Jacobian finite-
differenced on the warp/sample coords (no new sampler value+gradient interface
— that's deferred to Phase 3B). Never-worse-than-seed backtracking guard;
infinity (`w=0`) patches refined. New `crates/sfmtool-core/src/patch/keypoint_subpixel.rs`
+ `PatchCloud.refine_keypoints` PyO3 binding. Validated: synthetic recovery
< 0.02 px (finite / two-view / infinity), consensus sharpness rises, all guard
cases keep the seed; 14 Rust tests + 4 Python tests; `cargo test/clippy/fmt`
green. Branch ready to merge; not wired into `_embed_patches.py` yet (optional,
follows the gate).

### Phase 3 — Optimize each, independently

Both optimization paths are built independently and validated. **Order within
Phase 3 is flexible** — they touch disjoint code (A in `keypoint_localize.rs`,
B in `keypoint_subpixel.rs` + `camera/remap.rs`/`warp_map.rs`) and can be done
in parallel or in either order; pick by appetite. **Each lands fully validated
on its own merits**, not against the other.

- **Phase 3A — Grid optimizations (Spec A):**
  - The hand-rolled **AVX2** centered-`f32` register-blocked windowed-ZNCC
    kernel, with an `avx2_matches_scalar` equivalence test (mirroring the fronto
    cache's `resample_avx2_matches_scalar`).
  - Then the **integer-`i16`** stage 2 (box-window normalization) gated on
    argmax agreement *and* a real speedup vs. the f32 kernel.
  - Per-(view, round) planar deinterleave moved into the cache build (so the
    AVX2 kernel reads planar, contiguous f32 directly — the spec's intended
    layout).
  - Exit: `search_shift` collapses; the per-call planar deinterleave is gone;
    `localize` re-profile shows the expected drop on dino.

- **Phase 3B — Subpixel optimizations (Spec B):**
  - The sampler value+gradient variants (`sample_bilinear_with_grad_u8`,
    `remap_aniso_with_grad`) + `WarpMap::get_jacobian`, per the subpixel spec's
    "Design details" section. Removes the per-GN-step 5× FD render cost the MVP
    carries.
  - The incremental running-sum consensus with low-frequency IRLS weight
    refresh (the spec's "per-move (Gauss-Seidel) incremental" variant) — gated
    on a measured convergence/accuracy improvement vs. the single-pass-frozen
    MVP.
  - Inverse-compositional ECC only if the precompute clearly amortizes against
    the refreshing consensus.
  - Exit: B's per-view cost matches its asymptotic shape; the analytic Jacobian
    matches the MVP's FD Jacobian within sampling noise (an equivalence test).

### Decision gate — measure the production-path fork (after Phase 3A and 3B)

With **optimized** A (incl. the `m > 1` supersampled variant on AVX2/`i16`) and
**optimized** B both available:

- Compare supersampled-grid sub-pixel (`m = 2, 3, …`) keypoints against the LK
  keypoints — the **accuracy gap** and the **wall-clock cost** of each on real
  datasets (dino, seoul_bull, seattle_backyard, kerry_park).
- Decide the production path: grid-only (`m > 1`) where its accuracy suffices,
  grid + LK where quality matters, or a selectable choice (the most likely
  outcome). **The LK stays the high-accuracy / ground-truth option regardless**
  — the gate sizes *when* the cheaper grid is good enough to be the default,
  not whether LK exists.
- Wire the chosen path into `_embed_patches.py` (if not already) and document
  the default in the relevant CLI / spec.

## Interleave vs. sequential — verdict (updated)

**Mostly sequential through the MVPs, then parallel through the optimizations,
then a single measurement.** Spec A's structural layer (Phases 0–1) went first
because it's the producer, was partly done, and carries the bigger raw win.
The minimal Spec B (Phase 2) followed so both algorithms exist end-to-end.
Phase 3A and 3B touch disjoint files and can run in parallel — neither is
speculative, both are independently justified, and the production-path gate
runs only on the optimized output of both, so its decision rides on real per-µs
data, not on which MVP happened to have the worse shortcut.

## Risks / notes

- **Behaviour change, not equivalence (Phase 1).** The cache's integer-tracked
  reads build the consensus from integer-aligned cores; validated by
  kept-view/registration agreement, not bit-equivalence. (The search *kernel*
  itself stays equivalent — `search_shift_matches_reference*` holds.)
- **Supersampled-grid cost grows ~`m²`** — Phase 3A's AVX2/`i16` kernel must
  carry this; the gate weighs accuracy gained against the (now optimized) cost.
- **B-opt equivalence (Phase 3B).** The analytic sampler Jacobian must match the
  MVP's FD Jacobian within sampling noise — an explicit equivalence test, so the
  switch is a pure speedup with no behaviour drift.
- **Prototype parity.** The scripts on this branch are the validated references;
  reproduce their behaviour and reuse them as oracles rather than re-deriving.
- **Maturin rebuild.** Both specs touch `sfmtool-core` reached through PyO3;
  `pixi run maturin develop --release` before Python-side tests/profiling.
