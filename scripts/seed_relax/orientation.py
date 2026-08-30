# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""The global orientation bit, by parallax-weighted cheirality.

Provenance: the study's `sign/signlib.py` `cheirality_readings` (493-556),
reduced to the two readings the chain uses.  The near-versus-far rank statistic
beside them is not carried: it is a second reading of the same bit, and the
weighted observation vote is the one the census settled on.

Pairwise baseline directions determine the constellation only up to the point
reflection ``c -> -c``: negating every centre negates every triangulated depth
and changes nothing else, because the angular least squares is linear in the
centres and its matrix does not contain them.  Structure sitting in FRONT of
the cameras is the one physical statement that separates the two, and a
cluster's cheirality statement is worth exactly the parallax it was measured
with: a cluster inside the member's bound is a bearing whose depth sign is a
coin toss.

The vote itself is the core's
``sfmtool_core::geometry::translation_averaging`` kernel, reached through
``sfmtool._sfmtool.geometry.orientation_reading``; what is here is the world
rays the member's rotations produce.  See
``specs/core/geometry/translation-averaging.md``.
"""

from __future__ import annotations

import numpy as np


def angw_bit(m, per_frame, placed, tol):
    """The orientation reading of one constellation.

    Returns a dict: ``angw`` (the parallax-weighted front-minus-behind vote,
    in radians), ``obs_front`` / ``obs_total`` / ``obs_frac`` (the same vote
    unweighted), ``angw_per_obs`` and ``margin_frac`` (the two readings in the
    units they are read in), and the graduation census the pass sees on the
    way.  ``angw < 0`` says the constellation should be reflected.

    The reading is exactly antisymmetric under ``c -> -c``, so one pass
    describes both orientations and the second is arithmetic that is already
    known."""
    from sfmtool._sfmtool.geometry import orientation_reading

    frames = sorted(placed)
    slot = {f: k for k, f in enumerate(frames)}
    rays, point_of_ray, frame_of_ray = [], [], []
    for f in frames:
        cl, local, _rows = per_frame[f]
        # The rays are emitted frame by frame in this order, so a cluster is
        # first seen from the lowest-numbered frame that saw it and its own
        # rays stay in frame order -- the grouping the vote reads.
        rays.append(np.asarray(local, float) @ np.asarray(m.rot[f], float))
        point_of_ray.append(np.asarray(cl, np.int64))
        frame_of_ray.append(np.full(len(cl), slot[f], np.int64))
    empty = np.zeros((0, 3))
    return dict(
        orientation_reading(
            np.ascontiguousarray(np.stack([placed[f] for f in frames]), float)
            if frames
            else empty,
            np.ascontiguousarray(np.concatenate(rays)) if rays else empty,
            np.concatenate(point_of_ray) if point_of_ray else np.zeros(0, np.int64),
            np.concatenate(frame_of_ray) if frame_of_ray else np.zeros(0, np.int64),
            float(tol),
        )
    )
