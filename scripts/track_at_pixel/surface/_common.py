# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""What the surface diagnostics share: loading a dataset, and small patch moves."""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from context import DatasetContext  # noqa: E402
from dataset import prepare  # noqa: E402

UP = np.array([0.0, 0.0, 1.0])


def load(dataset: str):
    """The dataset's context and an edited reconstruction holding every point."""
    from sfmtool._sfmtool.reconstruction import EditedReconstruction

    ds = DatasetContext(prepare(dataset, None, quiet=True))
    return ds, EditedReconstruction(ds.recon)


def unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def ang(a, b):
    """The angle between two directions, ignoring sign, in degrees."""
    return float(np.degrees(np.arccos(np.clip(abs(unit(a) @ unit(b)), -1, 1))))


def tangent(n):
    a = np.array([1.0, 0, 0]) if abs(n[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = unit(np.cross(n, a))
    return e1, np.cross(n, e1)


def to_cam(ds, track):
    c = ds.cameras[int(track.observations[0]["image"])].center
    return unit(c - np.asarray(track.position, float))


def tilt(ds, edited, track, n):
    """The track tilted to ``n`` (faced toward its first camera)."""
    from sfmtool._sfmtool import bench as B

    n = n if n @ to_cam(ds, track) > 0 else -n
    track, _ = B.tilt_patch(track, edited, tuple(float(x) for x in n))
    return track


def move(edited, track, d):
    """The track translated by the world vector ``d``."""
    from sfmtool._sfmtool import bench as B

    pl = track.placement
    u, v = unit(pl["u_halfvec"]), unit(pl["v_halfvec"])
    m = unit(np.cross(u, v))
    track, _ = B.translate_patch(
        track, edited, (float(d @ u), float(d @ v), float(d @ m))
    )
    return track


def zncc_and_shift(ds, edited, track):
    """Each observation's ZNCC and correlation peak offset (px), evaluated in place."""
    from sfmtool._sfmtool import bench as B

    track, _ = B.evaluate(track, edited, ds.pyramids)
    z = [o.get("track", {}).get("zncc") for o in track.observations]
    s = [o.get("track", {}).get("seed_shift_px") for o in track.observations]
    return (
        np.asarray([np.nan if x is None else x for x in z], float),
        np.asarray([np.nan if x is None else x for x in s], float),
    )
