# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Density-seeded clustering of descriptors with mean-shift refinement.

Materialises clusters without transitive merging: process descriptors
densest-first, and let each unclaimed seed gather its within-`T` cross-image
neighbours (one per image, nearest the centre), then claim them. Membership is
"within `T` of the centre" — a fixed point — so chains A–B–C with A–C > `T`
cannot merge, unlike connected components. The result is a hard partition of
descriptors into candidate tracks.

With ``refine_iters > 0`` (and a ``forest`` to query) each cluster is refined by
mean-shift: re-query the index at the cluster's *mean* descriptor — a synthetic
point that need not coincide with any indexed descriptor — and re-gather the
within-`T` neighbours, one per image. Centring on the mean can reach members
closer to the cluster centroid than to the original seed, so the cluster drifts
toward its descriptor-space density centre. One step usually suffices.
"""

from __future__ import annotations

import numpy as np


def seed_claim_clusters(
    idx: np.ndarray,
    dst: np.ndarray,
    image_label: np.ndarray,
    T: float,
    min_size: int = 2,
    descriptors: np.ndarray | None = None,
    refine_iters: int = 0,
    forest=None,
    stats: dict | None = None,
) -> np.ndarray:
    """Cluster descriptors by density-ordered seeding, then mean-shift refine.

    Args:
        idx: ``(N, k)`` neighbour indices (includes self), from the index.
        dst: ``(N, k)`` neighbour distances aligned with ``idx``.
        image_label: ``(N,)`` source image index per descriptor.
        T: match radius; a neighbour joins iff its distance ≤ T.
        min_size: drop clusters spanning fewer than this many images.
        descriptors: ``(N, D)`` descriptors, required when ``refine_iters > 0``.
        refine_iters: mean-shift iterations (0 = seed only). Each step re-queries
            ``forest`` at the cluster mean; needs both ``descriptors`` and
            ``forest``.
        forest: a ``KdForest`` (or any object with ``query(X, k)``) used to
            re-query at the synthetic cluster mean during refinement.

    Returns:
        ``(N,)`` int64 cluster id per descriptor, or -1 if unclustered.
    """
    n, k = idx.shape
    img = image_label
    rng = np.arange(n)
    valid = (dst <= T) & (idx != rng[:, None]) & (img[idx] != img[:, None])
    order = np.argsort(-valid.sum(1), kind="stable")

    claimed = np.zeros(n, dtype=bool)
    labels = np.full(n, -1, dtype=np.int64)

    def one_per_image(members: np.ndarray, mdist: np.ndarray) -> np.ndarray:
        """Resolve candidate members to one feature per image (nearest first)."""
        if members.size == 0:
            return members
        mimg = img[members]
        o = np.lexsort((mdist, mimg))
        mi = mimg[o]
        keep = np.empty(len(o), dtype=bool)
        keep[0] = True
        keep[1:] = mi[1:] != mi[:-1]
        return members[o][keep]

    def gather(c: int) -> np.ndarray:
        """Centre `c`'s within-T cross-image unclaimed neighbours + `c`,
        resolved to one feature per image (nearest the centre)."""
        m = valid[c] & ~claimed[idx[c]]
        members = np.concatenate(([c], idx[c][m]))
        mdist = np.concatenate(([0.0], dst[c][m]))
        return one_per_image(members, mdist)

    def regather_mean(members: np.ndarray) -> np.ndarray:
        """Re-query the index at the cluster mean (a synthetic point) and resolve
        the within-T, unclaimed neighbours to one feature per image."""
        mean = descriptors[members].astype(np.float32).mean(0)
        mu8 = np.clip(np.round(mean), 0, 255).astype(np.uint8)
        qi, qd = forest.query(mu8[None, :], k=k)
        qi = qi[0].astype(np.int64)
        qd = qd[0]
        m = (qd <= T) & ~claimed[qi]
        return one_per_image(qi[m], qd[m])

    can_refine = refine_iters > 0 and descriptors is not None and forest is not None
    n_seeds = n_moved = members_changed = cid = 0
    for s in order:
        if claimed[s]:
            continue
        members = gather(s)
        before = frozenset(members.tolist()) if stats is not None else None
        moved = False
        if can_refine:
            for _ in range(refine_iters):
                if members.size < 2:
                    break
                new = regather_mean(members)
                if frozenset(new.tolist()) == frozenset(members.tolist()):
                    break
                members = new
                moved = True
        if stats is not None and members.size >= min_size:
            n_seeds += 1
            n_moved += int(moved)
            if before != frozenset(members.tolist()):
                members_changed += 1
        if members.size >= min_size:
            claimed[members] = True
            labels[members] = cid
            cid += 1
        else:
            claimed[s] = True  # isolated seed; drop
    if stats is not None:
        stats.update(
            clusters=n_seeds, center_moved=n_moved, members_changed=members_changed
        )
    return labels
