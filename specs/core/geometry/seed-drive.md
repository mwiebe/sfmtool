# Seed Drive Loop

The drive loop is the consumer that judges the candidate set while it is
still being produced. It pulls one exploration pass at a time from the
seed hypothesis source ([seed-hypothesis-loop.md](seed-hypothesis-loop.md)),
completes and measures the members that pass produced, judges the whole
accumulated set with the selection pass's gates
([seed-candidate-evaluation.md](seed-candidate-evaluation.md)), and hands
the source a claim map built from the SURVIVING members only. A refused
member's claim is withdrawn, so the next complement re-opens exactly
where that member's structure claimed.

The set is still the product. Verdicts annotate members; they never
remove one.

## Round

One round is one pull and everything taken on it.

**Pull.** The source yields the pass's committed candidates, in commit
order. The pass is the batch grain because it is the complement's
natural resumption point: the source forms the next admission at the
head of the next pull, from the claim map standing then.

**Complete.** Each new finite candidate gets its rotation-only
far-field sibling, fit independently of the candidate it pairs with and
of every other candidate, and each new layer is relaxed into a finite
sibling ([seed-relaxation.md](seed-relaxation.md)). Both are readings
of one member, so a round completes its own batch without waiting for
the rest of the set. The capture's own far-field layer is a reading of
no single candidate and is taken once, on the final set.

**Measure.** The evaluation battery runs on the round's NEW members
only. Members already measured keep the readings they were measured
with; nothing is re-measured because the set grew. The capture's
conditioning floors (the resection inlier floors and the rotation
support spread bar) are drawn on the first population that resolves
them and passed into every later round, so a floor is never re-derived
on one round's smaller population, which would move the readings of
every member already measured.

**Judge.** The whole accumulated set is judged every round, on the
stored evidence, by the same rank / refuse / trim machinery the offline
selection pass applies. A member's capture-relative readings are taken
over its own capture's median, so a verdict is a statement about the set
the member stands in and cannot be settled once and carried forward.

## Survivor claims

The claim map handed back for the next complement is the union of the
surviving members' own claims.

* A member the gates KEEP contributes the claim it stamped.
* A member the gates REFUSE contributes nothing. Its claim is
  withdrawn whole, and the complement re-opens where it claimed.
* A member the gates TRIM restates its claim on the frames its core
  kept: its retained clusters are re-counted on those frames, and a
  cluster the core no longer sees at least twice stops being an
  explained point and stops claiming. This is the same point cull the
  trim itself applies, taken on the claim triangulation rather than on
  the released one, so the restatement lives in the cluster space the
  original claim was stamped from.

Only finite candidates ever stamp a claim, so only they are read here.

## Sufficiency

The loop stops on sufficiency when both hold:

* at least one member the gates do not refuse, and
* the surviving-claims complement is not a new admission, which is the
  source's own exhaustion test (the complement is empty, or it is the
  admission that produced this pass).

## Stopping

The loop also stops when:

* **The source is exhausted.** No further pass is available: the
  complement is not a new admission, a pass produced no reconstruction,
  or the source hit the budget it was given.
* **The budget is reached.** The accumulated candidate count reaches
  the total the run was given, which is the same resource bound the
  default driver passes and the same cap the source applies mid-pass.
* **No progress.** A pass whose members are ALL refused, whose claims
  are therefore all withdrawn, and whose surviving-claims complement is
  back to the admission that produced it. Pulling again would explore
  the same admission. This is what keeps claim withdrawal from costing
  the loop the termination the accumulating map gives it.

## Gates

The fleet gates are the same file the offline selection pass reads.
Without a gates file, refusals are off: every member keeps a `keep`
verdict, the claim map is the union of every member's claim, and the
pulls are ordered by coverage alone.

## Product

The product is the ordinary candidate-solves directory
([seed-hypothesis-loop.md](seed-hypothesis-loop.md)), with everything
the run ever committed in it, refused members included, so a rejected
candidate opens in the viewer exactly like a kept one. Its
classification travels with it.

Each manifest entry carries a `drive` block:

| Field | Meaning |
|-------|---------|
| `round` | the drive round the member arrived in |
| `verdict` | `keep`, `trim` or `refuse` |
| `verdict_reason` | the sentence the verdict was taken on |
| `readings` | the named channel readings behind it |
| `conditioning_limited` | the readings that were non-measurements |
| `rank` | the member's rank within its model family |
| `trimmed_release_file` | the core's artifact, on a taken trim |
| `trimmed_from` | the member's own artifact, on a taken trim |
| `frames_dropped` | the frames the cut removed |

A trimmed member keeps BOTH artifacts as cross-referenced siblings: the
member as it was committed, and the core the trim states, each entry
naming the other file.

The manifest gains a `drive` block of its own naming the round count,
the rule the loop stopped on, the budget, whether gates were supplied,
the verdict counts and the surviving set's coverage report. The full
selection report is written as `rung2.json` beside the manifest, inside
the product directory.
