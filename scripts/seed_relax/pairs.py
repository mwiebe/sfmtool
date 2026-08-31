# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Candidate row pairs inside one image, within a per-row reach.

Both the reconciliation and the hand-over ask one question of an image's rows
-- which other rows lie within this row's own disk -- and differ only in what
they then test, so the enumeration is stated once and the tests stay in the
stages.

The enumeration is the core's ``spatial::keypoint_reach`` kernel, bound as
:func:`sfmtool._sfmtool.analysis.keypoint_pairs_within_reach`; this module is
the shape the stages read it in, one image at a time.  The kernel orders an
image's rows by column and takes each disk's candidates off one contiguous run
of that order, found by binary search at ``x -/+ reach`` and then filtered by
true Euclidean distance, so a pair that comes back is already inside the asking
row's reach and a stage's own containment test restates it rather than
narrowing it.  Batching is the kernel's own and is a work bound rather than a
threshold: the pairs are the same pairs in the same order at any grain.

See ``specs/core/analysis/keypoint-reach.md`` for the contract.
"""

from __future__ import annotations

import numpy as np


def image_candidates(x, y, reach):
    """Yield ``(i, j, d)``: every row ``j`` inside row ``i``'s reach.

    ``x`` and ``y`` are the rows' pixel coordinates and ``reach`` each row's
    own query radius; a row whose reach is not finite asks nothing.  ``d`` is
    the pair's separation in pixels, so a caller testing against ``reach[i]``,
    against ``reach[j]``, or against both reads one distance.  Row ``i`` is
    never its own candidate.

    The rows handed in are one image's, so the indices that come back are
    positions in these arrays.  The stream is read as a stream because the
    relation is quadratic in a crowded image and the caller consumes it batch
    by batch; nothing about the pairs depends on how many batches arrive."""
    from sfmtool._sfmtool.analysis import keypoint_pairs_within_reach

    reach = np.ascontiguousarray(np.asarray(reach, float))
    xy = np.ascontiguousarray(
        np.stack([np.asarray(x, float), np.asarray(y, float)], axis=1)
    )
    one_image = np.zeros(len(reach), np.int64)
    big, small, d = keypoint_pairs_within_reach(one_image, xy, reach)
    if len(big):
        yield big, small, d


def image_slices(slot_i):
    """``(image, row positions)`` for every image, in image order."""
    slot_i = np.asarray(slot_i, np.int64)
    for img in np.unique(slot_i):
        yield int(img), np.nonzero(slot_i == img)[0]
