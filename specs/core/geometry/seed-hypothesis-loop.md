# Seed Hypothesis Loop

The seed stage develops and commits a SET of candidate reconstructions.
A candidate is one full seed exploration (probe, widen, photometric
verify, focal scan, release) over an admitted cluster selection. The
first pass admits the whole selection; each later pass admits the
coverage complement of everything the candidates before it claimed.
Every distinct finalist of every pass commits, and the loop chooses
between none of them: the product is the set.

## Coarse admission

The loader admits the N coarsest clusters of the file's selection, N
being the cluster budget (`SFMTOOL_SEED_RUNG1`, default 3000; an
explicit `0` keeps every cluster). A cluster's coarseness is its widest
member's patch half-extent in image pixels, read off the stored member
affines with no `.sift` access: half the file's refine radius times the
mean of the member affine's two column norms. Any member over the bar
qualifies the cluster.

The bar is stated as a POPULATION, not a threshold: rank by radius
descending with cluster id ascending among ties, keep the first N. A
threshold lets the kept population span orders of magnitude across a
fleet, so runs at one bar are not one working set. When N is at least
the cluster count the cut is a no-op and the handle stays exactly as
loaded.

Coarse features are the alias-free evidence. On a repeating scene
texture (a tiled floor, a brick wall, a railing) fine features match
self-consistently at false lattice offsets, and the aliased basin is
internally clean and can outnumber the true one. A feature wider than
the repeat period cannot alias that way, so the coarse admission is the
admission on which basin structure is legible.

## Member order

The loader fixes the order of every cluster's members by content before
any stage reads them: the cluster's stored reference member first, then
the others by image index, then by keypoint position. Every selection the
loop solves on (the coarse cut, a coverage complement, a group-local
re-admission) is repackaged in that order, and the referee's observation
arrays carry it too. The order a file lists members in is therefore never
read by a solver:

- a stage that takes one observation per cluster as its reference ray
  (the rotation-only layers, the ladder's far-field reading, and
  `core_parallax`, which measures each point's widest angle from that ray)
  takes the stored reference member when its frame is posed, and otherwise
  the posed member of lowest image index;
- the bundle adjustments, which stop at an evaluation cap and sum their
  rows in observation order, see one row order for one file content.

Permuting the members of every cluster of a file leaves the product
bit-identical ([measurements](seed-hypothesis-loop-measurements.md#pick-stability-after-the-deterministic-fixes-2026-10-09)).

## Capture-level measurements

The pairwise focal vote, and where escalation confirms one the
camera-model verdict, is computed once over the FULL admission's pair
graph (the population as it stood before the coarse cut), before any
pass runs. Every candidate reads the same vote. The vote is a property
of the capture, not of a candidate: it is the independent referee each
release is measured against, so no pass re-derives it from its own
restricted pair graph.

Each pair's vote is read over `SFMTOOL_VOTE_DRAWS` RANSAC draws (default
5) and is the log-space median of its draws
([focal-vote.md](focal-vote.md#draws-per-pair)), so the referee does not
follow one draw's sample sequence. The ray-space pair inits that choose
the rotation core's gauge pair and each fisheye seed group's starting
pair read the same way: a pair's parallax and cheiral count are the
medians over that many draws, and its pose is the draw at the median
parallax. One draw per pair moved `SeoulBull`'s vote between 251 and 323
px over seeds 0 to 9 of one file; five draws move it between 245 and 298
px, and it settles slowly with more draws (254 to 290 px at nine)
([measurements](seed-hypothesis-loop-measurements.md#pick-stability-after-the-deterministic-fixes-2026-10-09)).

## Coverage claim

A committed candidate claims the image area its retained structure
samples, at the resolution the evidence itself samples it.

The retained structure is every cluster with a finite triangulated
position in the released geometry's full triangulation. The claim is
TRANSITIVE over cluster membership: a retained cluster is an explained
3D point, so its members stamp in **every** image they appear in, posed
or not. A candidate that poses a handful of a long capture's frames
still claims its structure's footprint capture-wide.

The claim is an occupancy grid per image, not a pixel bitmap. Each
image's cell size is the median nearest-neighbour distance among the
retained members' keypoints in that image: coverage measured at the
spacing the matcher actually sampled the scene, fine on dense texture
and coarse on sparse, with the pixel scale of the capture divided out.
A cell holding at least one retained member's keypoint is claimed.
Images with fewer than two retained members claim nothing (no spacing
exists to measure). Claims accumulate across committed candidates as a
per-image union of claimed cells (each candidate stamps into its own
grid geometry; a later test evaluates against every accumulated grid).

## Complement admission

The next pass admits the source selection minus the claimed clusters. A
cluster is claimed when more than half of its members fall in claimed
cells of their images. The complement is expressed as a cluster-id
restriction of the stage's selection handle (`select_clusters` with
`restrict_cluster_ids`, see
[cluster-selection.md](../features/cluster-selection.md)), so a
complement is itself an ordinary derived selection: it carries
provenance, and every stage downstream reads it exactly like the
unrestricted one. No stage applies a claim predicate of its own; the
selection file is the admission.

## Loop

The first pass explores the full admission. The loop then derives the
claims of every candidate that pass committed, forms the complement
selection, and explores it. The claim ORDERS the complement queue and
never gates it: it decides where the next admission starts, not whether
there is one, so no candidate's footprint can veto a group.

The loop ends when the complement is not a new admission (it is empty,
or nothing was claimed at all), when a pass produces no reconstruction
(no seed group, or every release posed no image), or at the candidate
budget. The budget is a resource bound, not a judgment: the generator
has no opinion about which candidates are worth keeping, and the cap
only bounds what a pathological capture can cost.

Termination is structural while the claim map only accumulates: a
committed candidate claims at least the clusters whose members it
retained, so the complement strictly shrinks with every committed pass.

### Pull contract

Candidate production is a SOURCE, and the loop above is what a consumer
that pulls from it without judging sees. The source yields one PASS at
a time: that pass's committed candidates, in commit order. The
accumulated set, the accumulated claims and the pass index are the
source's own state, which a consumer reads between pulls. The coverage
complement is formed at the HEAD of a pull rather than at the tail of
the pass that stamped the claims it stands on, so a consumer that stops
pulling never pays for an admission it does not explore.

Each committed candidate's claim stands on its own: the per-image grid
map is derivable from that candidate's retained clusters alone, and the
working claim map is the union of a chosen set of such maps. The
cluster test ORs over every grid an image carries, so the map answers
the same question whichever order the members' claims entered it. A
consumer therefore withdraws a member's claim by leaving that member
out of the union it hands back, and a member whose frames were cut
restates its claim by re-counting its retained clusters on the frames
the cut left: a cluster the surviving frames no longer see twice is no
longer an explained point, and stops claiming.

The budget is a parameter of the source with the semantics above. A
consumer that passes no budget gets an uncapped source and owns the
stop, including the termination guarantee claim withdrawal costs; the
driving consumer and its stopping rules are
[seed-drive.md](seed-drive.md).

## Group-local re-admission

Each attempt re-admits clusters for its own image set: the seed groups'
frames plus everything covisible with them in the attempt's own graph.
A group of five frames on its own would rank coarseness off five
viewpoints and leave the widen ladder nothing beyond them; its covisible
neighbourhood is the part of the capture the seed can grow into.

Over that image set the `n_local` coarsest clusters are kept, with
coarseness measured ON THOSE IMAGES (a cluster's widest member there,
which is the question the window's own solve depends on), eligibility
the loader's span bar counted on the same images, and the same stable
ordering the capture-wide cut uses. The selection is derived from the
PRE-cut handle, so a group can reach clusters the capture-wide cut
dropped. Admitted clusters keep their FULL member lists, so a frame
outside the neighbourhood still contributes wherever it sees one: the
neighbourhood bounds what the ranking is measured on, not what the solve
may reach. `n_local` is the capture budget unless
`SFMTOOL_SEED_LOCAL_ADMISSION` overrides it. An image set that carries
nothing eligible leaves the attempt on its own working set.

## Choice among seed groups

A pass explores its admission as a ladder of ATTEMPTS, one per working set
(the admission itself, then successively thinned ones). An attempt takes up
to 8 seed groups in covisibility order. On a parallax-poor capture it first
tries the rotation core, and a rotation core that clears the commit bar is
the attempt's outcome. Otherwise every seed group is probed at the probe
focal before any of them is judged. The probe focal is the raw pairwise
focal vote, with no bias correction; under an equidistant context it is the
verdict's equidistant focal, and with no vote it is `0.9 * max(w, h)`. It is
not a rung of the [focal scan's lattice](#focal-scan); snapping it to the
nearest rung is proposed in
[seed-probe-lattice-amendment.md](../../drafts/seed-probe-lattice-amendment.md).

- A probe is MEASURABLE when its inlier fraction (2 px) reaches
  `max(15%, 0.5 * best)`, where `best` is the highest inlier fraction over
  every group's probe in the attempt. The gate therefore does not depend on
  which groups were probed before it.
- Every measurable probe whose core parallax clears the near-static gate is
  finished: widened, photometrically verified, and checked against the
  commit bar ([Rank](#rank)).

The attempt's outcome is the finished group that the OUTCOME ORDER puts
first among those that clear the commit bar. When none clears it, the
outcome is the first by the same order among the reach-healthy but
focal-blind outcomes (kept frames and reach at the bar, scan spread under
it), then among the starved ones, then the near-static and the
unmeasurable fallbacks. The outcome order compares, in turn:

1. exploration reach, higher first: the share of the pass's images the
   posed set is connected to in the pass's own covisibility graph (8 shared
   clusters to an edge);
2. kept frames after the photometric verify, more first;
3. the median reprojection residual of the posed observations with a
   triangulated point, at the probe focal, lower first;
4. the seed group's image names, sorted, compared lexicographically, first
   first.

Every term is a property of the finished outcome, so the order the groups
were tried in decides nothing.

The rotation core is also the fallback of the ladder's first attempt (the
admission itself; thinned working sets never try it). When no seed group
clears the commit bar and none has left a starved outcome, that attempt
tries the rotation core after the seed groups if it has not tried it
already, on any capture rather than only a parallax-poor one. A rotation
core that clears the commit bar is the attempt's outcome; otherwise its
finished outcome, or its probe when its parallax is under the near-static
gate, joins the fallbacks above.

The cost of the outcome order is the widen and verify of groups a
first-come rule would not have reached; the probe and finish memos carry a
group that a later pass repeats at no further cost. The candidates a pass
commits are still the ladder's finalists (below), pulled from the source
one pass at a time ([Pull contract](#pull-contract)); the order only decides
which finished group stands for an attempt.

## Focal scan

The focal scan reads the posed geometry's inlier fraction at a series of
fixed focals, and its results set the commit bar's scan spread, the focal
the release starts from, and the `flat_scan` and `edge_scan` flags.

**Lattice.** The scan's focals are rungs of one lattice per capture,
`max(w, h) * 1.15^k` px for integer `k`, kept inside the scan band (the
field-of-view band of [focal-vote.md](focal-vote.md) under an equidistant
context; the pinhole plausibility floor `0.3 * max(w, h)` and no upper bound
otherwise). The structure-free focal (the pinhole vote, the equidistant
verdict, or the nominal probe focal when there is no vote) chooses only the
WINDOW: the 5 rungs nearest it, shifted to stay inside the band. A vote that
moves by less than half a step from its nearest rung moves no scanned focal,
and two votes either side of a midpoint scan windows that share 4 rungs.
The commit bar's spread is measured over this window.

**Extension.** The release's scan evaluates the window at a light
adjustment budget. When the best rung is at an end of the window, the scan
adds the next rung past that end and repeats, until the best rung is
interior, 3 rungs have been added on that side, or the band ends. A window
that only sat too low or too high is extended past its peak this way.

**Winner.** The REFIT BAND is every rung of the extended scan within 5
points of inlier fraction of the best rung and at most 2 rungs from it, so
at most 5 rungs. Each is refit at a heavier adjustment budget, and the best
refit wins. The structure-free focal breaks a tie only: among the refits
within half a point of the best one (`SFMTOOL_REFIT_TIE`, default 0.005), the
rung nearest it in log-focal wins. Because the tied rungs are lattice rungs,
the vote moves the winner only by crossing the log-midpoint of two of them.

**Release.** The release walks the focal from the winner with a free-focal
adjustment, at most three rounds, and stops when the focal moves by under
1% or leaves 15% of the winner. It keeps the latest round whose inlier
fraction is within half a point of the best seen, the winner's own
included. The winner is a lattice rung rather than an optimum in focal, so
on a flat-topped scan a walk that ties it is kept, and a walk that loses to
it by more than the tie is not. The release reads the same
`SFMTOOL_REFIT_TIE`.

The half-point tie and the release's tying-step rule were set on the two
ground-truth captures, `SeoulBull` and `KerryPark480`
([measurements](seed-hypothesis-loop-measurements.md#pick-stability-after-milestone-b-2026-10-09)),
and are validated against the fleet's references in the next fleet
measurement.

**`edge_scan`.** The flag is set when the extended scan still peaks at its
top rung and rises to it monotonically (no drop of more than half a point),
which is the upward affine escape: the inlier fraction keeps improving as
the focal grows, so the structure does not bound the focal from above. A
peak at the bottom rung sets no flag.

## Ladder dedup

A pass runs its exploration over successively thinned working sets and
ends with several finalists. They are ranked by score (scan spread,
floored to zero below the observability bar, then coverage: posed count
times capture reach, saturated at 60%) and de-duplicated: a finalist
whose relative rotations agree with an already-kept one within the
pose-noise scale (median relative-rotation disagreement at or below 5°
over the frames they both pose) collapses into it, and the better-scored
copy stands for both. Finalists with disjoint posed sets, or fewer than
two shared frames, never collapse.

That is dedup, not judgment: it removes one answer found twice and never
chooses between two different ones. Everything surviving it commits.

## Rank

Each committed candidate records its released focal, released inlier
fraction, capture-level coverage reach, scan spread, confidence flags,
and the log-focal distance between its release and the bias-corrected
capture-level vote. A candidate QUALIFIES when the structure-trust gates
all hold: the commit bar (posed count, reach, scan spread), the release
inside the corrected vote band, and no flat-scan, edge-scan or
near-static-seed verdict. The commit bar's posed count is
`min(8, cap - 1)` kept frames for the core cap `SFMTOOL_SCAN_CAP`
(default 8, so 7 frames): the probe grows a seed group to the cap, and a
bar equal to the cap fails on the first frame lost in growth. The same
bar decides whether an attempt's outcome commits and whether the ladder
stops early. Coverage reach is measured on the
CAPTURE-LEVEL covisibility graph, the full admission's, for every
candidate alike: reach asks how much of the capture a solve connects to,
and a complement's smaller admission must not deflate the answer for a
solve that genuinely spans it.

The rank is the recorded order of the set: the first qualified candidate
first, commit order otherwise. It is ADVISORY and decides nothing.
Ranking, refusal and trimming belong to the selection pass that reads
the stored evaluation evidence, see
[seed-candidate-evaluation.md](seed-candidate-evaluation.md).

The pick is sensitive to sub-percent changes in the member shapes of the
cluster-patches file it reads: on `KerryPark480`, moving under 0.4% of the
kept members to shapes of higher ZNCC changes which candidates are committed
and which pass against the ground truth, and two such files that differ in a
few dozen members give a first qualified candidate that fails and one that
passes
([measurements](../patch/cluster-patch-refinement-measurements.md#subset-with-the-loop-as-the-default-2026-10-09)).
The rank does not cause it: the first candidate is the first pass's, and the
change enters that candidate's own exploration and reaches every later pass
through the complement queue
([measurements](seed-hypothesis-loop-measurements.md#pick-stability-under-small-changes-to-the-cluster-file-2026-10-09)).
Member order, the vote's single draw per pair, the commit bar equal to the
core cap, a scan grid centred on the vote, the first-tried choice among seed
groups and a resection floor read on every trim round do not reach it: each is
replaced by the rule its section above states, the last by the trim rule of
[rotation-locked-resection.md](rotation-locked-resection.md#mechanism). On the
22 perturbed files of the two ground-truth captures the first candidate passes the ground truth on every file, but it still
changes on 10 of them
([measurements](seed-hypothesis-loop-measurements.md#pick-stability-after-milestone-b-2026-10-09)):
on `KerryPark480` the dropped members change which skeleton images the
rotation core resects and which frames the widen admits, and on `SeoulBull`
they move one seed group under the commit bar's 60% reach floor. Those
stages read the evidence the dropped members carried, and a first candidate
that never moves under a 0.3% change of the cluster file is not a goal of
this stage.

None of the signals qualification and the rank read compares a candidate's
photographs resampled into its patches. A photometric candidate score that
gates qualification and orders the rank is proposed in
[seed-photometric-candidate-score.md](../../drafts/seed-photometric-candidate-score.md).

## Product

`sfmr/candidate_solves/` is the product: one `h<NN>.sfmr` per committed
member and a `manifest.json` naming them, self-contained, with no other
file to read and no stamp in any path. Members are the finite
candidates, the rotation-only far-field layers, and the relaxed siblings
those layers commit (see
[seed-relaxation.md](seed-relaxation.md)), in commit order.

A member's release is release-grade: poses and points, no consensus
bitmaps and no patch frames. The artifacts are written under the
capture's own camera model, so a fisheye capture's members densify and
reproject through the equidistant context, never the pinhole default,
and a member carrying a released lens stamps that lens.

Because a release carries no patch frames, nothing reads its patches before
the choice or moves the chosen member's patches to a finer resolution after
it. A cluster file written with `--piecewise` stores per-cell displacements
([cluster-patch-refinement.md](../patch/cluster-patch-refinement.md#piecewise-refinement)),
and [cell-plane-normals.md](../patch/cell-plane-normals.md) turns them into
patch normals once poses exist, but the seed reads neither. Frames built from
them behind a precision gate are proposed in
[cell-plane-normal-precision-gate.md](../../drafts/cell-plane-normal-precision-gate.md), a photometric score that reads them in
[seed-photometric-candidate-score.md](../../drafts/seed-photometric-candidate-score.md), and the
coarse and fine resolutions the score and the chosen member use in
[two-tier-patch-density.md](../../drafts/two-tier-patch-density.md).

The manifest carries the run's stamp, the coarse admission's population
figures, the vote block with the admission the referee measured on, the
advisory rank's first entry, and one entry per member: its model, its
camera and focals, its metrics and flags, the admission its solve ran
on, the frames it was seeded from and the frames it posed, and the
evaluation block the battery attached (see
[seed-candidate-evaluation.md](seed-candidate-evaluation.md)).

The directory is replaced whole. Releases are written into a staging
sibling as they are committed, the manifest joins them there, and only
then is the destination removed and the staging directory renamed onto
it. A reader therefore sees the previous product, or nothing, or this
one, and never a partial set or a manifest naming a release that is not
there.
