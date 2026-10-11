# Cluster Review Strips

**Status:** Draft. Decided: what a cluster strip shows, how a blind review set
of cluster members is drawn and recorded, and which statistics are reported
beside it. Built as the prototype
[`scripts/cluster_strips.py`](../../scripts/cluster_strips.py). Not decided:
whether cluster strips become a mode of `sfm inspect --strips`, and whether a
review's verdicts become a bar a gate change must pass (see
[Open questions](#open-questions)).

## Purpose

A feature matcher groups detections from several photographs that it takes to
be the same point on a surface, and cluster-patch refinement
([cluster-patch-refinement.md](../core/patch/cluster-patch-refinement.md))
then decides, member by member, whether each photograph really shows the same
small piece of surface as the cluster's reference member. Two kinds of evidence
judge those decisions. Statistics under ground-truth camera poses say whether a
member lies on the epipolar line of its reference, which catches most wrong
members but not one that is wrong along that line. A person looking at the two
pieces of image side by side catches that too, but sees only a sample. A
**cluster strip** is a picture of one cluster as the refinement sees it: the
reference's patch, then each member's patch sampled in the reference's frame,
ordered by how much smaller or larger the member's footprint is, each labelled
with its scores and its verdict. A **blind review set** is a sample of members
drawn from a population a statistic picks out, shown without their scores or
verdicts, and judged right or wrong by a person; the record of the judgement is
kept beside the measurements it is compared with.

## Interface

The prototype is one script with three subcommands. Each reads a
cluster-patches `.matches` file, the output of `sfm cluster-patches`, and
resolves the file's workspace for the photographs.

```bash
# Score every measured member against its reference, with the epipolar
# distance under a ground truth's poses when one is given.
pixi run python scripts/cluster_strips.py score ws/matches/x-clusters-patches.matches \
    --gt ws/sfmr/ground_truth.sfmr -o x-scores.npz

# Render chosen or random clusters as strips.
pixi run python scripts/cluster_strips.py strips ws/matches/x-clusters-patches.matches \
    --scores x-scores.npz --random 20 --seed 1 -o x-strips.png

# Draw a blind review set over one or more entries.
pixi run python scripts/cluster_strips.py review -o review/ \
    --entry KerryPark480 kp-clusters-patches.matches kp-scores.npz \
    --entry DnDTabletop dnd-clusters-patches.matches dnd-scores.npz \
    --cases 40 --controls 8 --seed 20261010
```

`score` writes one member-parallel `.npz`: whether the member was measured, its
plain and blur-matched scores, which tile the blur-matched score blurred and by
how much, the self-similarity radius of its tile, its scale ratio, and its
epipolar distance (NaN without `--gt`). The scores are computed once and read
by both other subcommands, so a strip and a review case show the same numbers.
`score` prints how closely its plain score reproduces the file's stored
`member_zncc`, which checks that its tiles are the kernel's.

## What a cluster strip shows

One row per cluster: a label panel (cluster index, member count, the span of
the kept members' scales), the reference member's tile, then one tile per
measured member (status `kept`, `rejected_low_zncc`, `rejected_shift`,
`rejected_unlocalizable_refined` or `rejected_unlocalizable_cells`; the others
have no refined geometry and are not drawn).

- **The tile** is the refinement's template grid: `resolution` samples per side
  over `[-patch_size/2, patch_size/2]²` keypoint-frame units, both read from the
  file's `refine_options`, mapped into the member's photograph by its refined
  position `p` and absolute shape `S` (`x = p + S·u`), sampled bilinearly from
  the pyramid level the kernel's rule picks, with pixel centres at `+0.5`.
  Because `S` is the shape the member was scored at, a right member's tile shows
  the same piece of surface as the reference's, in the same orientation and at
  the same size. Tiles are enlarged by pixel replication, so each grid sample
  is visible as the square the score reads.
- **The order** is by scale ratio, `√|det S|` of the member over the
  reference's, smallest first. The reference is the cluster's largest-scale
  member, so most ratios are below 1, and the row reads from the members that
  saw the surface from furthest away to those closest to the reference's view.
- **Each label** carries the scale ratio, the status, the plain score `p`, the
  blur-matched score `b` and the image index. The frame colour is green for a
  kept member, blue for a rejected member whose blur-matched score clears the
  file's `min_zncc`, and red for any other rejected member.

## The scores

- **Plain** is the windowed ZNCC of the member's tile against the reference's,
  under the refinement's window (a Gaussian of sigma 0.5 confined to the
  inscribed disk), per colour channel, averaged over the reference's textured
  channels. It is the number the refinement's `min_zncc` judges, recomputed.
- **Blur-matched** is the same correlation after the sharper tile of the pair
  is blurred to the other's sharpness, by the rule of
  [blur-matched-zncc.md](../core/patch/blur-matched-zncc.md): a tile is the
  sharper when its self-similarity semi-major axis is shorter than the other's
  semi-minor axis; it is blurred so that its semi-major axis reaches that
  semi-minor axis, capped at 2 grid px, and only where that target is at least
  1.25 times its own semi-major axis; otherwise the pair is read plain. Either
  tile may be the one blurred, since the reference is chosen by scale, not by
  sharpness.
- **Self-similarity radius** is the overlap reading of the member's tile, the
  number the refinement's gate at the refined shape judges against 2.5.
- **Epipolar distance** is the distance, in px of the member's photograph, of
  its refined position from the epipolar half-line of its reference's position
  under the ground truth's poses, as in
  [cluster-patch-refinement-measurements.md](../core/patch/cluster-patch-refinement-measurements.md#gates-at-the-refined-shape-2026-10-09).
  A member is *wrong by the ground truth* above `max(3, 5·m)` px, `m` the
  entry's median over kept members.

## Drawing a blind review set

A review set answers one question about one population, and the population is
named by the statistic that raised the question. For the blur-matched score,
the population is the **disagreement set**: members the plain score rejects
(`rejected_low_zncc`) whose blur-matched score, with one tile of the pair
blurred, clears `min_zncc`, and which would pass the refinement's later gates
(shift at most `max_shift_px`, radius at most 2.5). A pair the rule reads plain
scores the same both ways, so a member read plain that clears the bar differs
from the kernel's score only by the script's sampling and is left out. The set
is drawn with a fixed seed:

- **Cases** are stratified by scale-ratio bin (`< 1/4`, `1/4–1/2`, `1/2–0.7`,
  `0.7–1.4`, `≥ 1.4`), one case at a time from each (entry, bin) pool in turn,
  so every bin and entry is represented while its pool lasts.
- **Controls** are drawn from the two populations where the scores agree: kept
  members that also clear the bar blur-matched, and rejected members that the
  blur-matched score also rejects. They say how often the reviewer calls an
  accepted member wrong and a rejected one right.
- **The order** of all cases is shuffled, and case ids are assigned after the
  shuffle, so neither the id nor the position names the population.

Each case is one image: the reference's tile and the member's tile side by
side, over crops of their two photographs with each footprint outlined and a
tick on the edge the tile's top row samples. Nothing in the image names a
score, a status, the scale ratio or the dataset. The question for a single
member is whether it is a right correspondence of the reference, the same
piece of surface, or a wrong one; the answers are *right*, *wrong* and
*unsure*, with an optional note. A review that compares two shapes of one
member, as the review of moved shapes did, shows the two as sides A and B in a
random order per case, and the answers are A, B, *same* and *both wrong*.

## Recording the review

`review` writes, into one directory:

- `c01.png`, `c02.png`, …, one per case;
- `cases.csv`, the key: case id, dataset, population (`disagreement`,
  `control_accepted_both`, `control_rejected_both`), member and cluster
  indexes, reference member, both image names, scale ratio and its bin,
  status, the plain score, the file's stored score, the blur-matched score,
  which tile was blurred and by how much, shift, radius, epipolar distance, and
  the member's position and shape;
- `review.csv`, the answer sheet: case id, choice and note, empty;
- `index.html`, a contact sheet of every case with a choice and a note field
  per case, which keeps the answers in the browser while the review is under
  way and exports the `review.csv` rows.

The reviewer opens `index.html` and does not open `cases.csv` until every case
is answered; no answer is changed after the key is opened. The two files are
then joined on the case id into one archive CSV,
`<spec>-human-review-<date>.csv`, beside the measurements file it supports:
the case id, dataset, population, scale-ratio bin, member and cluster
indexes, scale ratio, plain and blur-matched scores, epipolar distance and the
ground truth's verdict from the key, then the reviewer's answer and note, as
[cluster-patch-refinement-human-review-2026-10-10.csv](../core/patch/cluster-patch-refinement-human-review-2026-10-10.csv)
is. The images and the page are not kept; the member index, with the
cluster-patches file the measurements section names, renders any case again.

## Folding a review into the measurements

The review becomes a dated section of the measurements file of the spec it
bears on, in the form of
[Human review of moved shapes](../core/patch/cluster-patch-refinement-measurements.md#human-review-of-moved-shapes-2026-10-09):

- **Question**, in one paragraph: which population, and why the statistics
  alone do not settle it.
- **Data**: the files and the commit that wrote them, the population's size per
  entry, the seed, the number of cases and controls per bin and entry.
- **Protocol** and **limits**: blind, randomised, who reviewed, and what could
  still be recognised (a control whose two tiles are identical, say).
- **Result**: per scale-ratio bin, the counts of *right*, *wrong* and *unsure*
  among the cases; the controls' answers; and the agreement of the reviewer's
  verdict with the ground-truth verdict, case by case, with the cases where the
  two differ listed, since a member wrong along the epipolar line reads right
  by the ground truth.
- **What this decides**, in plain terms, or what result would settle it.

The section links the archive CSV. The first review drawn with the prototype
is folded in this way as
[Blur-matched scores across scales](../core/patch/cluster-patch-refinement-measurements.md#blur-matched-scores-across-scales-2026-10-10),
with the statistics of the next section beside it. The second, drawn from
members where the colour reading of
[image-photometric-model.md](image-photometric-model.md) and the blur-matched
ZNCC disagree, is folded in as
[Measured: step 2](image-photometric-model.md#measured-step-2-2026-10-10), with
the archive
[cluster-patch-refinement-human-review-2026-10-10b.csv](../core/patch/cluster-patch-refinement-human-review-2026-10-10b.csv).

## The statistics beside a review

A review samples a population; these say how large it is and what the ground
truth says about all of it. They are reported per entry, over all measured
members, with the scale ratio binned as for the cases (finer bins where the
data reach them):

- per bin, the number of kept members and of `rejected_low_zncc` members, and
  for each the median epipolar distance and the share wrong by the ground
  truth;
- per bin, the disagreement set (rejected plain, accepted blur-matched), with
  its median epipolar distance and wrong share, before and after the later
  gates; and the members kept plain that the blur-matched score would drop;
- per bin, how often the reference's tile is the one blurred and how often the
  member's;
- per cluster, the scale span of the accepted members (largest `√|det S|` over
  smallest, the reference included), plain against blur-matched, as
  percentiles and a histogram; and the number of clusters that accept both a
  coarse member (feature size, the mean column norm of `S`, above 40 px) and a
  fine one (at most 10 or 20 px).

The blur-matched accepted set for these counts keeps, per cluster and image,
the one member with the highest blur-matched score, as the refinement keeps one
member per image.

## Open questions

- **Where the strips live.** `sfm inspect --strips` renders a reconstruction's
  points ([inspect-command.md](../cli/reconstruction/inspect-command.md#point-strips---strips)),
  and its layout (labels, reference tile, then a row of tiles) is the layout a
  cluster strip needs. The strips subpackage exposes only the two commands'
  entry points, so the prototype draws its own rows with the scripts' shared
  drawing helpers. A cluster mode of `sfm inspect --strips` reading a
  `.matches` file would share the layout.
- **The review as a bar.** Whether a change to the refinement's gates must
  leave the disagreement set's reviewed *right* share above some level, and how
  large a sample makes that level meaningful.
