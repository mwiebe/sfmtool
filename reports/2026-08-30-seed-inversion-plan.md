# Seed pipeline inversion of control: candidate source + controlling loop

Date: 2026-08-30. Branch: `bootstrap-core-migration`. Prerequisite landed:
the legacy pipeline deletion (`7e127e79` / `658c935a`, parity PASS), which
left `scripts/exp_fast_seed.py` with exactly one pipeline whose hypothesis
loop is already generator-shaped.

## Motivation

Today candidate production decides when to stop (the candidate budget, the
complement-exhaustion test) before any candidate has been judged, and the
selection pass (`scripts/exp_seed_rung2.py`) reads the stored evidence after
the door has closed: it can rank, refuse and trim, but it cannot ask for
more. A refused candidate's coverage claim also stands forever, suppressing
exploration exactly where a garbage candidate claimed. The inversion makes
production a resumable, uncapped source and puts the judged consumer in the
driver's seat: pull candidates, co-evolve quality evaluation with structure
growth, and dip back for more while the surviving set leaves material
evidence unexplained.

## Target architecture

```
capture_context(ws)                      once: load_clusters, focal vote,
                                         capture_covisibility, the Rung
candidate_source(ctx, budget=None)       a Python generator
    yields one PASS at a time: the committed candidates, releases on disk
    owns: selection handle, claims, pass index, memos;
    complement admission happens between pulls

controller (rung 2's drive loop):
    while True:
        batch = next(source)                       or StopIteration
        far-field sibling per new finite candidate
        relax the new layers (the 9-stage chain)
        evaluation battery on the new members only
        verdicts = rank/refuse/trim over ALL members (fleet gates)
        claims: only SURVIVORS' footprints stand for the next complement
        if sufficient(verdicts) or source exhausted: break
    peer corroboration over the final set
    capture-wide far-field layer + its relaxation
    write_candidate_solves (verdicts recorded; nothing discarded)
```

## Adopted decisions

- Claims from survivors only; a refusal withdraws the member's claim so the
  complement re-opens where it claimed. A no-progress stop preserves
  termination (see Phase 2).
- A trimmed member's claim restates on its surviving core's clusters, not
  its original footprint.
- `CANDIDATE_BUDGET` stays 8 as the controller's total cost governor for the
  A/B; revisit with fleet evidence.
- The controller lives in `exp_seed_rung2.py` as a `drive` mode; the offline
  `select` / `derive-gates` / `channels` modes stay for stored releases.
- No gates file means refusals off: the controller degenerates to
  coverage-driven pulls only, which is the safe default for a bare run.
- EVERYTHING IS REVIEWABLE. Every candidate the run ever committed keeps
  its `.sfmr` in the product, refused members included, so a rejected
  candidate can be opened in the Explorer exactly like a kept one. Its
  classification travels with it: the manifest entry carries the verdict
  (keep / trim / refuse), the named readings it was taken on, and the
  drive round the member arrived in; the full selection report
  (`rung2.json`) is written into the product directory beside the
  manifest. A trimmed member keeps BOTH artifacts, the original and the
  trimmed core, cross-referenced. The verdicts annotate the set; they
  never delete from it.

## Phase 1: extraction (content-neutral)

Pure restructuring of `scripts/exp_fast_seed.py`:

- `capture_context()`: everything `main()` does before the hypothesis loop
  (load, vote, fisheye routing, `capture_covisibility`), returned as one
  object.
- `candidate_source(ctx, budget)`: the loop becomes a generator yielding
  each pass's committed batch after the commit/claim block, then forming
  the complement and continuing. The budget is a parameter with today's
  exact check semantics (mid-batch cutoff, `budget_overflow` evolution
  reasons); the default driver passes `CANDIDATE_BUDGET`, a controller
  passes None for an uncapped source.
- `main()` becomes the default driver: drain the source, then the unchanged
  tail (far-field layers, `relax_far_layers`, rank, `compare_to_reference`,
  `attach_evaluation`, `write_candidate_solves`).

Acceptance: the standing parity bar. `run_arm` arm `postphase1` against the
`postlegacy` baseline on the two default workspaces, `parity.py` VERDICT
PASS (content-identical `candidate_solves/` products).

> _Status (2026-08-30): Done -- commit 88af6988 (`CaptureContext` /
> `capture_context()` / `candidate_source(ctx, budget)` generator, `main()`
> as the draining driver). AST neutrality audit: identical statement
> sequences except the parameterized budget test. Parity postlegacy vs
> postphase1 on KerryPark360 + 20250906_211742965: VERDICT PASS
> (PARITY-postlegacy-vs-postphase1.txt)._

## Phase 2: the controller (behavior change, human-validated)

`exp_seed_rung2.py drive <workspace>`:

1. Pull one pass from `candidate_source` (batch grain = pass, the
   complement's natural resumption point).
2. Complete the new members: rotation-only sibling per finite candidate
   (independent fits, safe per batch), relax those layers, run the battery
   on the new members only. Peer corroboration is computed once over the
   final set before the product write.
3. Judge the accumulated set with the fleet gates: rank, refuse, trim.
4. Form the next complement from the SURVIVING members' claims only.
5. `sufficient()`: at least one member the gates do not refuse, AND the
   surviving set's complement is not a new admission (the source's own
   existing test). Stop also on: source exhausted; the total budget; and
   NO PROGRESS (a pass whose members are all refused and whose
   surviving-claims complement did not shrink), which restores the
   termination guarantee that claim withdrawal alone would break.
6. Product: everything ships in `sfmr/candidate_solves/`, refused members
   included; the verdicts, per-member rounds and the coverage report ride
   in the manifest, `rung2.json` lands beside it, and trimmed members keep
   original and core as sibling files (see EVERYTHING IS REVIEWABLE
   above).

Gates: `derive-gates` over the fleet-chain run's 42 releases, written
beside the study harness; the population must match the 9-stage chain.

> _Status (2026-08-30): Commits landed -- d823e5a6 (per-member claim /
> measurement / far-layer seams in exp_fast_seed.py) + 1075fac2 (`drive`
> mode in exp_seed_rung2.py, specs/core/geometry/seed-drive.md, the pull
> contract in seed-hypothesis-loop.md). Branch rebased onto upstream
> ea110409 (8-point solver null-space fix) with the bindings rebuilt.
> Parity postphase1 vs postphase2 on the two default workspaces: VERDICT
> PASS (PARITY-postphase1-vs-postphase2.txt), certifying commit 1 neutral
> on the default driver path and the upstream solver fix inert on these
> captures. Gates derived: shipped/gates-chain.json (42 releases, 1064
> members, 3 model families). Open: the specimen drive runs and their
> human Explorer verdicts, then the 42-entry fleet A/B._

Acceptance, in order (per the house validation rules):
- Specimens one at a time with human Explorer verdicts between: the
  20250906 entry, KerryPark360, and one capture the old loop demonstrably
  under-explored (a garbage-claim member from the fleet tail census).
- Then the 42-entry fleet A/B against the Phase 1 baseline: coverage,
  direction/centre/map metrics, wall cost (the controller adds per-round
  battery and relaxation cost but can stop earlier; both directions need
  numbers). Movers get human review before adoption.

## Phase 3: controller moves (recorded, out of scope here)

Single insertions once the loop exists: both-orientations requests for
hard-bit mirror members (the source accepts a seed-group plus orientation
constraint), cross-member merges, the runaway-report gate (refuse and
replace instead of relaxing garbage).

## Spec updates (with Phase 2)

- `specs/core/geometry/seed-hypothesis-loop.md`: the Loop section gains the
  pull contract (the source yields passes; the judged survivors' claims
  form the complement). Mechanism only.
- A new spec for the drive loop (`sufficient()`, stopping rules, gates
  input, product semantics), placed with the seed geometry specs.
