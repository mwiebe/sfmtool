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

## Capture-level measurements

The pairwise focal vote, and where escalation confirms one the
camera-model verdict, is computed once over the FULL admission's pair
graph (the population as it stood before the coarse cut), before any
pass runs. Every candidate reads the same vote. The vote is a property
of the capture, not of a candidate: it is the independent referee each
release is measured against, so no pass re-derives it from its own
restricted pair graph.

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
[cluster-selection.md](../../formats/cluster-selection.md)), so a
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

Termination is structural: a committed candidate claims at least the
clusters whose members it retained, so the complement strictly shrinks
with every committed pass.

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
near-static-seed verdict. Coverage reach is measured on the
CAPTURE-LEVEL covisibility graph, the full admission's, for every
candidate alike: reach asks how much of the capture a solve connects to,
and a complement's smaller admission must not deflate the answer for a
solve that genuinely spans it.

The rank is the recorded order of the set: the first qualified candidate
first, commit order otherwise. It is ADVISORY and decides nothing.
Ranking, refusal and trimming belong to the selection pass that reads
the stored evaluation evidence, see
[seed-candidate-evaluation.md](seed-candidate-evaluation.md).

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
