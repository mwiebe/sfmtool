# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Relative baseline lengths, from the depths each pair's own solve implies.

A pair's baseline direction fixes the pair's geometry only up to how long the
baseline is, so the two-view depths that come out of it are read in units of
THAT baseline.  A cluster two edges both see from the same frame therefore has
two depths for one world distance, and their ratio is the ratio of the two
baselines: the depths carry exactly the relative scale the directions cannot.

The whole graph is one fit of ``log z(edge, frame, cluster) = D(frame, cluster)
- x(edge)``, where ``x`` is the log baseline length and ``D`` the log world
depth.  That fit is the core's
``sfmtool_core::geometry::translation_averaging`` kernel, reached through
``sfmtool._sfmtool.geometry.relative_lengths``; what is here is the depth rows
the pair solve produces and the dict-keyed shape the chain speaks.  See
``specs/core/geometry/translation-averaging.md``.
"""

from __future__ import annotations

import numpy as np


def two_view_depths(u_i, u_j, d):
    """Depths along each ray at the closest approach, and the midpoint.

    ``c_i`` sits at the origin and ``c_j`` at ``+d``, so both depths are in
    units of the pair's own baseline."""
    a11 = np.einsum("ij,ij->i", u_i, u_i)
    a12 = -np.einsum("ij,ij->i", u_i, u_j)
    a22 = np.einsum("ij,ij->i", u_j, u_j)
    b1 = u_i @ d
    b2 = -(u_j @ d)
    det = a11 * a22 - a12 * a12
    with np.errstate(invalid="ignore", divide="ignore"):
        s = (a22 * b1 - a12 * b2) / det
        t = (a11 * b2 - a12 * b1) / det
    mid = 0.5 * (s[:, None] * u_i + (d + t[:, None] * u_j))
    return s, t, mid


def relative_lengths(keys, depths, rounds=None, min_tied=None):
    """``(lengths, scatter, n_tied)`` per edge of ``keys``, in that order.

    ``depths`` maps an edge key to ``(frames, clusters, z)``: one row per
    (frame, cluster) the edge's solve gave a positive depth for, ``z`` in units
    of that edge's own baseline.  Returns the relative length of every edge
    with the ``min_tied`` rows it needs, gauged to a median of one, the median
    absolute log residual each edge's own rows leave, and how many of its rows
    another edge also saw.  An edge without a length reads back as ``nan``.

    ``rounds`` and ``min_tied`` of ``None`` take the kernel's own constants."""
    from sfmtool._sfmtool.geometry import relative_lengths as kernel

    n_edge = len(keys)
    ee, ff, cc, zz = [], [], [], []
    for e, k in enumerate(keys):
        row = depths.get(k)
        if row is None:
            continue
        frames, clusters, z = row
        ee.append(np.full(len(z), e, np.int64))
        ff.append(np.asarray(frames, np.int64))
        cc.append(np.asarray(clusters, np.int64))
        zz.append(np.asarray(z, float))
    if not ee:
        none = np.full(n_edge, np.nan)
        return none, none.copy(), np.zeros(n_edge, np.int64)
    kw = {}
    if rounds is not None:
        kw["rounds"] = int(rounds)
    if min_tied is not None:
        kw["min_tied"] = int(min_tied)
    return kernel(
        np.concatenate(ee),
        np.concatenate(ff),
        np.concatenate(cc),
        np.concatenate(zz),
        n_edge,
        **kw,
    )
