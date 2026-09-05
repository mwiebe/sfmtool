# Draft: CSR-native observations for the intrinsics stack, and a from-matches entry

**Status:** Draft

Change proposal. The structure-free intrinsics stack (`focal_vote` and
`estimate_intrinsics`) takes its observations in the cluster-grouped (CSR)
form the `.matches` backbone stores, instead of the expanded member-parallel
form; a `from_matches` entry in core makes one call of the whole
file-to-estimate path in both languages. Outputs are **bit-identical** — the
change is representational, and every existing output pin must pass unchanged
(a moved pin is a bug, not a pin update).

## Why

- The kernel's contract is already CSR in disguise: `cluster_indexes` "must be
  nondecreasing (each distinct cluster is a contiguous run)", and the pair
  builder reconstructs the runs by scanning
  ([`focal_vote.rs`](../../crates/sfmtool-core/src/geometry/focal_vote.rs)
  `~:645`, `while cluster_indexes[run_end] == cid`). Expanded input can
  express invalid states the docs must forbid; `cluster_starts` cannot, and
  validation drops from O(observations) to O(clusters).
- Every real caller starts from a selection, whose native layout IS
  `cluster_starts` + member-parallel arrays; the expansion exists only to feed
  this API and is re-derived in Python at every call site.
- The convenience must live in Rust, not the binding, so the Rust and Python
  surfaces stay one API (the `resect_images` precedent: a core geometry entry
  taking `&MatchesData`, the binding a thin forwarder).

## Core (`sfmtool-core::geometry`)

- `focal_vote`, `focal_vote_with_min_disp`, `focal_vote_with_options`,
  `estimate_intrinsics` take
  `cluster_starts: &[u32]` (length `n_clusters + 1`, `starts[0] == 0`,
  nondecreasing, `starts[last] == n_observations` — validated up front),
  `member_images: &[u32]`, `member_positions: &[[f64; 2]]` (both length
  `n_observations`), `width`, `height`, options as today. The run scan becomes
  direct iteration over `starts` windows; **no arithmetic changes anywhere**
  — same observations in the same order produce the same bits.
- Positions stay `f64`. The file stores `f32` and widening is exact, so the
  arithmetic contract is untouched; whether the kernel's *internal* arithmetic
  wants `f32` is a separate, output-changing question under its own protocol.
- New `estimate_intrinsics_from_matches(&MatchesData, &IntrinsicsOptions)`
  beside the array form (the `resect_images` pattern): requires a
  clusters-bearing file (error otherwise, matching how `resect_images` words
  matches-borne errors), borrows `cluster_starts`/`member_images` directly,
  widens `member_positions` `f32 → f64` once, and derives `(width, height)`
  from the image table — **rejecting non-uniform dimensions**, which states
  and checks the shared-camera/centred-principal-point contract in one place
  instead of every caller silently picking `dims[0]`.

## Binding (`sfmtool-py`)

`estimate_intrinsics` accepts two forms and only these two:

- **Object form**: first argument a `MatchesFile` (a selection included),
  plus the existing `seed=` / `columns=` keywords — forwards to
  `estimate_intrinsics_from_matches`.
- **CSR form**: `cluster_starts`, `member_images`, `member_positions` (f64),
  `width`, `height`, keywords as today.

The expanded-rows form is removed; callers migrate. The `focal_vote` binding
migrates to the same two forms for consistency.

## Callers to migrate

- `src/sfmtool/_commands/estimate_intrinsics.py`: the command passes the
  opened `MatchesFile` (object form); `_load_observations` shrinks to
  whatever the report still reads off the handle, or disappears.
- `scripts/exp_fast_seed.py`: the wrapper converts its expanded working
  arrays with one exact `np.searchsorted(obs_c, np.arange(n_cl + 1))` at the
  call site (obs_c is nondecreasing by construction); `rung.vote_obs` and
  everything else stay as they are.
- `scripts/exp_seed_opening.py`: step 2 becomes the object form — the whole
  point of that script is measuring exactly this collapse.
- Tests (`focal_vote/tests.rs`, `estimate_intrinsics/tests.rs`, the binding
  and CLI pytest suites): synthetic fixtures build `starts` instead of
  expanded ids; **every pinned output value stays byte-for-byte**.

## Specs to update in the same change

- [`focal-vote.md`](../core/geometry/focal-vote.md) and
  [`estimate-intrinsics.md`](../core/geometry/estimate-intrinsics.md): the
  input-contract paragraphs (CSR shape, validation, the from-matches entry
  and its uniform-dimensions rule).
- `specs/cli/image-feature/estimate-intrinsics-command.md` if it names the
  array form.

## Verification

- `pixi run cargo fmt && pixi run cargo clippy --workspace && pixi run doc`;
  `pixi run cargo test --workspace` — all existing focal-vote and
  estimate-intrinsics pins green **unchanged**.
- `pixi run -e test maturin develop --release`, then `pixi run test`.
- End-to-end bit-identity against the standing fleet probe records is run by
  the operator after the change (the datasets are outside the implementation
  environment).
