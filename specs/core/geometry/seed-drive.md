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
stored evidence, by the same rank / refuse / trim / cull machinery the
offline selection pass applies. A member's capture-relative readings are
taken over its own capture's median, so a verdict is a statement about
the set the member stands in and cannot be settled once and carried
forward.

## Verdicts

A member carries one of four verdicts. The first three are survivors:
they stand in the ranking, they claim coverage, and they are counted as
members the set keeps.

| Verdict | What it says |
|---------|--------------|
| `keep` | no gated channel fired against the member |
| `trim` | the frames the walk moved carried the defect -- dropped, or re-seated where the structure supports another pose -- and the core that is left passes the member gates |
| `cull` | the points a cut removed carried the defect, and the core that is left passes the member gates |
| `refuse` | the defect is the member |

A refusal STANDS OUTRIGHT where a **pose-instability** channel fired and
did not fire alone. Those channels say the member's own poses do not
hold still: the settling family and its rotation-only twin, the hold-out
pose deltas of both model families, and the rotation cycles. Nothing
that can be taken out of a member repairs poses that will not hold, and
one channel on its own is one reading.

Every other member the gates would refuse reaches the **salvage** tier:
the evidence against it is about structure, or about a single reading,
and a member can be relieved of structure. Two remedies are tried in
order, and each stands only where the core it leaves passes the member
gates in its own right.

1. **Frame trim**, on the frames at the far end of the member's OWN frame
   population. Each naming per-frame channel is read against a quantile
   of that member's own readings of it, and a frame past the bar on any
   of them is a candidate; the fleet's own named frames are candidates
   too. This runs whether or not the fleet's bars localized the defect:
   a defect the fleet's frame bars do not localize is not thereby spread
   evenly over the member's frames.

   The cut is taken **a step at a time, smallest first**. Every
   per-frame reading of a member is taken against the structure the
   whole member holds, so a catastrophically mis-registered frame makes
   its sound neighbours read badly too: re-resecting a good frame
   against a cloud that still carries the bad frame's points measures
   the bad frame. A single naming pass cannot tell that contamination
   from a genuine spread of defects, and only a fresh measurement after
   a partial cut can.

   One step is: name the candidates on the readings standing now, rank
   them by **loudness** -- the loudest naming reading against the frame,
   as a multiple of its own bar -- take the frames past the **cliff** in
   that loudness population as the step's group, re-seat what that group
   supports and drop the rest, and re-measure the core the step leaves.
   The cliff is the widest successive RATIO in the loudness order: a
   loudness is a reading over its own bar, so ratios between loudnesses
   are scale-free and compare readings taken in degrees, in fractions
   and in correlations, and the widest of them is the one place the
   member's own population separates. Ties take the earlier place, which
   is the smaller cut; a population with no separation to find cuts one
   frame, the loudest. The criterion reads the population's own spacing
   and carries no magnitude of its own.

   A step **re-seats before it drops**. A stored pose was fit jointly
   with the frames beside it, so a subtractive move asks one question of
   a loud frame -- does the member read better without it -- and a frame
   that is only SEATED WRONG answers a different one: what pose it takes
   against structure held out from the whole group. The group goes to
   the resection primitive
   ([resect-image.md](../../gui/resect-image.md)) against the state the
   walk holds, which is the member minus what earlier steps dropped with
   the group STILL IN IT: the primitive holds the target set out itself,
   re-triangulating everything the set observes from the non-target
   frames alone before it estimates anything, so a state the group is
   part of is what makes the estimate independent of the group. The
   capture's own match graph carries the correspondences where the run
   names one, and the member's stored observations carry them where it
   does not. No bundle adjustment runs.

   The primitive's own acceptance is the criterion. A frame it accepts
   is **re-seated**: it keeps its place in the core at the pose the
   primitive estimated, with the points that pose re-triangulated, and
   the re-seated state is measured as it stands. A frame it refuses is
   **dropped** -- a frame the held-out structure gives no support at all
   is a stranger to it, and a frame whose estimate misses the gate is a
   frame the member cannot state a pose for. Each frame is tried once
   ever: a re-seated frame the fresh readings name again is dropped
   rather than re-seated a second time, so every step either drops a
   frame or tries one that was never tried. A member whose release ships
   no arrays to state a re-seated core from is not re-seated, because a
   core that cannot be measured is not a core this pass judges.

   The member's own frame bars are then **re-derived on the fresh
   readings** and the core is re-judged against the member gates, on the
   same terms whether the step dropped frames, re-seated them, or both.
   A core that passes ends the walk, and the trim is the union of the
   steps taken. Frames still named continue it. The share a trim may
   remove caps the frames DROPPED -- a re-seated frame is a kept frame
   and costs the cap nothing -- so a member whose fresh readings keep
   naming frames walks into that cap and its refusal stands. Each step
   pays one re-measurement and no step is measured twice. The verdict
   record carries the walk -- per step, the candidates and their
   loudness, where the cliff cut, what the resection said about each
   frame in the group, which frames it re-seated and which it dropped,
   and what the re-measured core said -- so a reviewer replays it
   without re-running the pass.
2. **Point cull**, on the points the fired channels indict. A channel
   with a per-point form of its own is read at that form; a channel with
   none is read on what the member's own arrays say about a point
   regardless of any channel -- whether its own rays re-triangulate to
   where it is stored, how much baseline priced that depth, how well it
   reprojects. Every bar is a quantile of the MEMBER'S own population of
   that signal, because there is no fleet population of points to rank a
   point in; a point past the bar on any indicting signal goes, and the
   cut is capped at the share of the cloud a trim may take of the
   frames. The cull may not break a link the member's frames had, and it
   writes the surviving cloud as a culled sibling of the member.

Where neither remedy leaves a core that passes, the refusal stands and
its reason records what was tried. Every refusal reason carries how many
distinct channels fired and the loudest of them, each reading measured
against its own bar so channels of different units compare.

## Survivor claims

The claim map handed back for the next complement is the union of the
surviving members' own claims.

* A member the gates KEEP contributes the claim it stamped.
* A member the gates REFUSE contributes nothing. Its claim is
  withdrawn whole, and the complement re-opens where it claimed.
* A member the gates TRIM restates its claim on the frames its core
  kept, which are the frames it never dropped and the frames it
  re-seated: its retained clusters are re-counted on those frames, and a
  cluster the core no longer sees at least twice stops being an
  explained point and stops claiming. This is the same point cull the
  trim itself applies, taken on the claim triangulation rather than on
  the released one, so the restatement lives in the cluster space the
  original claim was stamped from.
* A member the gates CULL restates its claim on the clusters its core
  kept, by the same mechanism: a culled point is not an explained point,
  so it stops claiming, and the frames the cull emptied restate with it.
  The cull is stated in the member's own cluster space, which is the
  claim's space unless the member solved on a group-local re-admission;
  where the two differ the restatement stands on the frames alone and
  the verdict records that it did.

Only finite candidates ever stamp a claim, so only they are read here.

## Sufficiency

The loop stops on sufficiency when all three hold:

* at least one member the gates do not refuse,
* at least one surviving claim, and
* the surviving-claims complement is not a new admission, which is the
  source's own exhaustion test (the complement is empty, or it is the
  admission that produced this pass).

The middle condition is what an empty claim map costs. Only finite
members stamp claims, so a set whose surviving members claim nothing
leaves the complement equal to the whole admission, which IS the
admission the last pass was formed from, and the exhaustion test reads
true for the reason that nothing was explained rather than for the
reason that everything was. A set that explains nothing has not answered
the capture, however many members it holds, so emptiness is never
sufficiency. The round's own line says plainly when the surviving claim
map is empty.

## Stopping

The loop also stops when:

* **The source is exhausted.** No further pass is available: the
  complement is not a new admission, a pass produced no reconstruction,
  or the source hit the budget it was given.
* **The budget is reached.** The accumulated candidate count reaches
  the total the run was given, which is the same resource bound the
  default driver passes and the same cap the source applies mid-pass.
* **No progress.** A complement that cannot move, together with either
  a pass whose members are ALL refused, or a surviving set with no claim
  at all. Both are the same situation: the claims are withdrawn or were
  never stamped, the complement is back to the admission that produced
  the pass, and pulling again would explore that same admission. This is
  what keeps claim withdrawal from costing the loop the termination the
  accumulating map gives it, and it is what an empty claim map stops on
  instead of stopping on sufficiency.

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
| `verdict` | `keep`, `trim`, `cull` or `refuse` |
| `verdict_reason` | the sentence the verdict was taken on, carrying the fired-channel count and the loudest reading |
| `readings` | the named channel readings behind it |
| `conditioning_limited` | the readings that were non-measurements |
| `rank` | the member's rank within its model family |
| `trimmed_release_file` | the core's artifact, on a taken trim |
| `trimmed_from` | the member's own artifact, on a taken trim |
| `culled_release_file` | the core's artifact, on a taken cull |
| `culled_from` | the member's own artifact, on a taken cull |
| `points_culled` | how many points the cull removed |
| `cull_signals` | the per-point signals the cull was taken on |
| `frames_dropped` | the frames the cut removed |
| `frames_reseated` | the frames the core keeps at a re-seated pose |

A trimmed or culled member keeps BOTH artifacts as cross-referenced
siblings: the member as it was committed, and the core the remedy
states, each entry naming the other file.

The manifest gains a `drive` block of its own naming the round count,
the rule the loop stopped on, the budget, whether gates were supplied,
the verdict counts and the surviving set's coverage report. The full
selection report is written as `rung2.json` beside the manifest, inside
the product directory.
