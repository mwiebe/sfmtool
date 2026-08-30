# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Candidate row pairs inside one image, within a per-row reach.

Provenance: the expansion the hand-over stage was written with
(`evict.covered_by_finer`'s x-sorted run and its chunked repeat), lifted out so
the reconciliation reads the same neighbourhood the hand-over does.  Both
stages ask one question of an image's rows -- which other rows lie within this
row's own disk -- and differ only in what they then test, so the enumeration is
stated once and the tests stay in the stages.

The rows of an image are ordered by column, so the candidates of a disk are a
contiguous run: a disk of radius ``reach`` cannot reach past a column further
than ``reach`` away.  The run is expanded in chunks because the pair count is
quadratic in a crowded image and the arrays are built whole; the chunk bound is
a memory bound and not a threshold, since the pairs are the same pairs and in
the same order at any value of it.
"""

from __future__ import annotations

import numpy as np

#: How many candidate row pairs one pass of the expansion holds at a time.
PAIR_CHUNK = 1 << 22


def chunk_ranges(cnt, chunk=PAIR_CHUNK):
    """``(start, stop)`` row ranges whose candidate pairs fit one pass."""
    cnt = np.asarray(cnt, np.int64)
    ends = np.cumsum(cnt)
    a = 0
    while a < len(cnt):
        base = ends[a - 1] if a else 0
        b = int(np.searchsorted(ends, base + int(chunk), side="right"))
        b = max(b, a + 1)
        yield a, min(b, len(cnt))
        a = b


def image_candidates(x, y, reach, chunk=PAIR_CHUNK):
    """Yield ``(i, j, d)``: every row ``j`` inside row ``i``'s reach, chunked.

    ``x`` and ``y`` are the rows' pixel coordinates and ``reach`` each row's
    own query radius; a row whose reach is not finite asks nothing.  ``d`` is
    the pair's separation in pixels, so a caller testing against ``reach[i]``,
    against ``reach[j]``, or against both reads one distance.  Row ``i`` is
    always its own candidate: a caller drops the self pair with the same test
    it uses for everything else it does not want."""
    x = np.ascontiguousarray(np.asarray(x, float))
    y = np.ascontiguousarray(np.asarray(y, float))
    reach = np.asarray(reach, float)
    xo = np.argsort(x, kind="stable")
    xs = x[xo]
    lo = np.searchsorted(xs, x - reach, side="left")
    hi = np.searchsorted(xs, x + reach, side="right")
    cnt = (hi - lo).astype(np.int64)
    cnt[~np.isfinite(reach)] = 0
    for a, b in chunk_ranges(cnt, chunk):
        take = cnt[a:b]
        total = int(take.sum())
        if not total:
            continue
        big = np.repeat(np.arange(a, b, dtype=np.int64), take)
        offs = np.arange(total, dtype=np.int64) - np.repeat(
            np.cumsum(take) - take, take
        )
        small = xo[np.repeat(lo[a:b].astype(np.int64), take) + offs]
        dx = x[small] - x[big]
        dy = y[small] - y[big]
        yield big, small, np.sqrt(dx * dx + dy * dy)


def image_slices(slot_i):
    """``(image, row positions)`` for every image, in image order."""
    slot_i = np.asarray(slot_i, np.int64)
    for img in np.unique(slot_i):
        yield int(img), np.nonzero(slot_i == img)[0]
