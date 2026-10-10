# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Photometric score of seed release candidates (prototype).

For each finite candidate under ``<workspace>/sfmr/candidate_solves/``, build
a patch frame at every finite point, render every posed view that observes the
point into that frame at its stored keypoint, and read how well the views
agree. Nothing is sampled at a reprojection of a candidate point: every tile is
anchored on the observation's stored keypoint (the candidate's inline
``keypoints_xy``, the cluster-patches member positions).

Frames. A release carries no patch frame, so one is built per point:

- centre: the candidate's point;
- normal: the cell-plane normal (``sfmtool.analysis.cell_plane_normals``) of
  the point's cluster in the workspace's ``*-clusters-patches.matches`` file,
  read through the candidate's own cameras and poses, where both of its axes
  are fixed; elsewhere the mean viewing direction (``--normal mean_viewing``
  uses that everywhere);
- extent: the cluster refinement's ``refine_radius`` times each observation's
  member-shape scale, back-projected (``extent="feature_size"``), median over
  the views.

Readings per point, at resolution ``R`` (default 12):

- ``ref_bm`` / ``ref_plain``: the median, over the point's non-reference
  observations, of the blur-matched / plain ZNCC of the observation's tile
  against the point's bitmap, the tile of the reference view the
  reference-view rule picks (``PatchCloud.render_bitmaps``,
  ``score_against_bitmap``);
- ``pair_plain`` / ``pair_bm``: the median off-diagonal entry of the member
  coherence matrix, plain and blur-matched
  (``PatchCloud.validate_member_coherence`` with ``keypoint_anchor=True``).

A candidate's score under each reading is aggregated three ways: the median
over its points, the 10%-trimmed mean, and the share of points above a
capture-relative bar (the median of that reading pooled over every candidate
of the capture, computed by the caller from the per-point values).

Usage::

    pixi run -e dev python scripts/seed_photometric_score.py <workspace> [<workspace> ...]
        [--out scores.json] [--resolution 12] [--normal cell|mean_viewing]
        [--files h00.sfmr ...]

Writes one JSON entry per workspace: per candidate the scored points, posed
views, per-reading aggregates, per-point values, and the wall time.
"""

import argparse
import json
import sys
import time
from pathlib import Path

import cv2
import numpy as np

from sfmtool.analysis import cell_plane_normals
from sfmtool.fileio import MatchesFile, read_sift_metadata
from sfmtool.geometry import RigidTransform
from sfmtool.patches import (
    CameraViews,
    ImagePyramidSet,
    PatchCloud,
    score_against_bitmap,
)
from sfmtool.reconstruction import SfmrReconstruction
from sfmtool.sift.file import get_sift_path_for_image

READINGS = ("ref_bm", "ref_plain", "pair_plain", "pair_bm")
NO_REF = 0xFFFFFFFF


def load_rgb(path):
    img = cv2.imread(str(path), cv2.IMREAD_COLOR | cv2.IMREAD_IGNORE_ORIENTATION)
    if img is None:
        raise FileNotFoundError(path)
    return cv2.cvtColor(img, cv2.COLOR_BGR2RGB)


class Workspace:
    """A seeded workspace: its cluster-patches file, decoded images, hashes."""

    def __init__(self, ws):
        self.ws = Path(ws)
        mfs = sorted(self.ws.glob("matches/*-clusters-patches.matches"))
        if len(mfs) != 1:
            raise SystemExit(
                f"{ws}: expected one *-clusters-patches.matches, found {mfs}"
            )
        m = MatchesFile(mfs[0])
        self.matches = m
        self.m_names = list(m.image_names)
        self.m_index = {n: i for i, n in enumerate(self.m_names)}
        self.cluster_starts = np.asarray(m.cluster_starts, dtype=np.uint32)
        self.reference_members = np.asarray(m.reference_members, dtype=np.uint32)
        self.member_images = np.asarray(m.member_images, dtype=np.uint32)
        self.member_features = np.asarray(m.member_features, dtype=np.uint32)
        self.member_status = np.asarray(m.member_status, dtype=np.uint8)
        self.member_positions = np.asarray(m.member_positions(), dtype=np.float32)
        self.member_shapes = np.asarray(m.member_affine_shapes(), dtype=np.float32)
        self.has_cells = bool(m.has_member_cells)
        if self.has_cells:
            self.cell_shift = np.asarray(m.member_cell_shift_px, dtype=np.float32)
            self.cell_status = np.asarray(m.member_cell_status, dtype=np.uint8)
        self.refine_radius = float(m.refine_radius)
        self.patch_size = float(m.refine_options["patch_size"])
        self.refine_resolution = int(m.refine_options["resolution"])
        self.member_cluster = (
            np.searchsorted(
                self.cluster_starts, np.arange(len(self.member_images)), side="right"
            )
            - 1
        )
        key = (
            self.member_images.astype(np.uint64) << np.uint64(32) | self.member_features
        )
        order = np.argsort(key, kind="stable")
        self.key_sorted = key[order]
        self.key_member = order
        self.images = {}
        self.hashes = {}

    def image(self, name):
        if name not in self.images:
            self.images[name] = load_rgb(self.ws / name)
        return self.images[name]

    def file_hash(self, name):
        if name not in self.hashes:
            try:
                meta = read_sift_metadata(get_sift_path_for_image(self.ws / name))[
                    "metadata"
                ]
                self.hashes[name] = bytes.fromhex(meta["image_file_xxh128"])
            except Exception:
                self.hashes[name] = bytes(16)
        return self.hashes[name]

    def members_of(self, m_img, feat):
        """Member index of each (matches image, feature) pair, -1 where none."""
        key = m_img.astype(np.uint64) << np.uint64(32) | feat.astype(np.uint64)
        pos = np.searchsorted(self.key_sorted, key)
        pos = np.minimum(pos, len(self.key_sorted) - 1)
        hit = self.key_sorted[pos] == key
        return np.where(hit, self.key_member[pos], -1)


def cluster_normals(W, rec, point_cluster):
    """Cell-plane normal per point (NaN rows where no normal is both-axes fixed)."""
    n_pts = len(point_cluster)
    out = np.full((n_pts, 3), np.nan)
    if not W.has_cells:
        return out, {}
    clusters = np.unique(point_cluster[point_cluster >= 0])
    if len(clusters) == 0:
        return out, {}
    lo = W.cluster_starts[clusters].astype(np.int64)
    hi = W.cluster_starts[clusters + 1].astype(np.int64)
    sizes = hi - lo
    idx = np.concatenate([np.arange(a, b) for a, b in zip(lo, hi)])
    starts = np.concatenate([[0], np.cumsum(sizes)]).astype(np.uint32)
    ref = W.reference_members[clusters]
    ref_local = np.where(
        ref == NO_REF, NO_REF, ref.astype(np.int64) - lo + starts[:-1]
    ).astype(np.uint32)
    # Poses over the matches file's images, NaN for an image the candidate does not pose.
    n_img = len(W.m_names)
    q = np.full((n_img, 4), np.nan)
    t = np.full((n_img, 3), np.nan)
    cam = np.zeros(n_img, dtype=np.uint32)
    rq = np.asarray(rec.quaternions_wxyz, dtype=np.float64)
    rt = np.asarray(rec.translations, dtype=np.float64)
    rc = np.asarray(rec.camera_indexes, dtype=np.uint32)
    for i, name in enumerate(rec.image_names):
        j = W.m_index.get(name)
        if j is not None:
            q[j], t[j], cam[j] = rq[i], rt[i], rc[i]
    res = cell_plane_normals(
        starts,
        ref_local,
        W.member_images[idx],
        W.member_status[idx],
        W.member_positions[idx],
        W.member_shapes[idx],
        W.cell_shift[idx],
        W.cell_status[idx],
        W.patch_size,
        W.refine_resolution,
        list(rec.cameras),
        cam,
        q,
        t,
    )
    names = list(res["determinacy_names"])
    both = names.index("both_axes")
    det = np.asarray(res["determinacy"])
    normal = np.asarray(res["normal"])
    local = {int(c): k for k, c in enumerate(clusters)}
    for p in range(n_pts):
        c = point_cluster[p]
        if c >= 0:
            k = local[int(c)]
            if det[k] == both:
                out[p] = normal[k]
    counts = {nm: int(np.sum(det == i)) for i, nm in enumerate(names)}
    return out, counts


def median_offdiag(mat):
    k = mat.shape[0]
    if k < 2:
        return np.nan
    v = mat[~np.eye(k, dtype=bool)]
    v = v[np.isfinite(v)]
    return float(np.median(v)) if len(v) else np.nan


def score_candidate(W, path, resolution, normal_mode):
    t0 = time.perf_counter()
    rec = SfmrReconstruction.load(str(path))
    names = list(rec.image_names)
    xyzw = np.asarray(rec.positions_xyzw, dtype=np.float64)
    tp = np.asarray(rec.track_point_indexes, dtype=np.int64)
    ti = np.asarray(rec.track_image_indexes, dtype=np.int64)
    tf = np.asarray(rec.track_feature_indexes, dtype=np.int64)
    kp = np.asarray(rec.keypoints_xy, dtype=np.float32)
    finite = xyzw[:, 3] != 0
    # Keep finite points only, so every patch is a planar surfel.
    keep_obs = finite[tp]
    pts = np.nonzero(finite)[0]
    remap = np.full(len(xyzw), -1, dtype=np.int64)
    remap[pts] = np.arange(len(pts))
    tp_f, ti_f, tf_f, kp_f = (
        remap[tp[keep_obs]],
        ti[keep_obs],
        tf[keep_obs],
        kp[keep_obs],
    )
    order = np.argsort(tp_f, kind="stable")
    tp_f, ti_f, tf_f, kp_f = tp_f[order], ti_f[order], tf_f[order], kp_f[order]

    # Observation -> cluster member, through the matches file's image index.
    m_img = np.array(
        [W.m_index.get(names[i], -1) for i in range(len(names))], dtype=np.int64
    )
    obs_m_img = m_img[ti_f]
    members = np.where(obs_m_img >= 0, W.members_of(np.maximum(obs_m_img, 0), tf_f), -1)
    mapped = members >= 0
    shapes = W.member_shapes[np.maximum(members, 0)].astype(np.float64)
    scale = 0.5 * (
        np.linalg.norm(shapes[:, :, 0], axis=1)
        + np.linalg.norm(shapes[:, :, 1], axis=1)
    )
    scale = np.where(mapped, scale, np.nan)
    pc = np.where(mapped, W.member_cluster[np.maximum(members, 0)], -1)
    point_cluster = np.full(len(pts), -1, dtype=np.int64)
    for p, c in zip(tp_f, pc):
        if c >= 0 and point_cluster[p] < 0:
            point_cluster[p] = c

    if normal_mode == "cell":
        nrm, det_counts = cluster_normals(W, rec, point_cluster)
    else:
        nrm, det_counts = np.full((len(pts), 3), np.nan), {}
    cell_normal = np.isfinite(nrm).all(axis=1)
    normals = np.where(cell_normal[:, None], nrm, 0.0)
    # Front-facing: towards the mean of the observing cameras.
    rq = np.asarray(rec.quaternions_wxyz, dtype=np.float64)
    rt = np.asarray(rec.translations, dtype=np.float64)
    views = CameraViews(
        list(rec.cameras), rq, rt, np.asarray(rec.camera_indexes, dtype=np.uint32)
    )
    w, x, y, z = rq.T
    rot = np.stack(
        [
            np.stack(
                [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)], -1
            ),
            np.stack(
                [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)], -1
            ),
            np.stack(
                [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)], -1
            ),
        ],
        1,
    )
    centres = -np.einsum("nji,nj->ni", rot, rt)
    to_cam = centres[ti_f] - xyzw[pts][tp_f, :3]
    to_cam /= np.maximum(np.linalg.norm(to_cam, axis=1, keepdims=True), 1e-300)
    mean_dir = np.zeros((len(pts), 3))
    np.add.at(mean_dir, tp_f, to_cam)
    flip = np.einsum("ij,ij->i", normals, mean_dir) < 0
    normals[flip] *= -1

    # Track points with no readable scale are dropped (a NaN refuses the cloud),
    # and so is a point that sits on an observing camera's centre, which a frame
    # cannot be sized at.
    ok_pt = np.ones(len(pts), dtype=bool)
    ok_pt[tp_f[~np.isfinite(scale)]] = False
    dist = np.linalg.norm(centres[ti_f] - xyzw[pts][tp_f, :3], axis=1)
    # The patch cloud refuses a camera-frame distance of 1e-6 or less.
    on_camera = dist <= max(2e-6, 1e-4 * float(np.median(dist)))
    ok_pt[tp_f[on_camera]] = False
    sel_obs = ok_pt[tp_f]
    pts2 = np.nonzero(ok_pt)[0]
    remap2 = np.full(len(pts), -1, dtype=np.int64)
    remap2[pts2] = np.arange(len(pts2))
    tp2 = remap2[tp_f[sel_obs]]
    ti2, tf2, kp2, sc2 = ti_f[sel_obs], tf_f[sel_obs], kp_f[sel_obs], scale[sel_obs]
    pos2 = xyzw[pts][pts2]
    nrm2 = normals[pts2]

    cloud = PatchCloud.from_tracks(
        views,
        np.ascontiguousarray(pos2),
        np.ascontiguousarray(tp2.astype(np.uint32)),
        np.ascontiguousarray(ti2.astype(np.uint32)),
        keypoint_scales=np.ascontiguousarray(sc2),
        normals=np.ascontiguousarray(nrm2),
        normal="stored",
        extent="feature_size",
        extent_value=W.refine_radius,
        exclude_points_at_infinity=True,
    )
    counts = np.bincount(tp2, minlength=len(pts2)).astype(np.uint32)
    emb = rec.clone_with_changes(
        positions=np.ascontiguousarray(pos2),
        colors=np.ascontiguousarray(np.asarray(rec.colors, dtype=np.uint8)[pts][pts2]),
        errors=np.ascontiguousarray(
            np.asarray(rec.errors, dtype=np.float32)[pts][pts2]
        ),
        track_point_indexes=np.ascontiguousarray(tp2.astype(np.uint32)),
        track_image_indexes=np.ascontiguousarray(ti2.astype(np.uint32)),
        track_feature_indexes=np.ascontiguousarray(tf2.astype(np.uint32)),
        observation_counts=counts,
        normals=np.ascontiguousarray(nrm2, dtype=np.float32),
        feature_source="embedded_patches",
        keypoints_xy=np.ascontiguousarray(kp2),
        image_file_hashes=[W.file_hash(n) for n in names],
        patches=cloud,
    )
    t_frames = time.perf_counter() - t0
    pyr = ImagePyramidSet(emb, [W.image(n) for n in names])
    bitmaps, refs = cloud.render_bitmaps(emb, pyr, resolution=resolution)
    bitmaps = np.asarray(bitmaps)
    refs = np.asarray(refs)

    cams = list(rec.cameras)
    cam_idx = np.asarray(rec.camera_indexes)

    poses = [
        RigidTransform.from_wxyz_translation(rq[i], rt[i]) for i in range(len(names))
    ]
    starts = np.concatenate([[0], np.cumsum(counts)]).astype(np.int64)
    cloud_pts = np.asarray(cloud.point_indexes)
    per = {r: np.full(len(pts2), np.nan) for r in READINGS}
    n_views = np.diff(starts)
    for k in range(len(cloud)):
        p = int(cloud_pts[k])
        lo, hi = int(starts[p]), int(starts[p + 1])
        if hi - lo < 2 or refs[p] < 0:
            continue
        patch = cloud[k]
        tiles, valid = [], []
        for o in range(lo, hi):
            i = int(ti2[o])
            tl = patch.render_view_tile(
                cams[int(cam_idx[i])],
                poses[i],
                pyr,
                image_index=i,
                keypoint=(float(kp2[o, 0]), float(kp2[o, 1])),
                resolution=resolution,
            )
            tiles.append(tl["samples"])
            valid.append(tl["valid"])
        out = score_against_bitmap(
            bitmaps[p], np.stack(tiles), valid=np.stack(valid), reference=int(refs[p])
        )
        others = np.arange(hi - lo) != int(refs[p])
        bm = np.asarray(out["blur_matched_zncc"])[others]
        pl = np.asarray(out["plain_zncc"])[others]
        if np.isfinite(bm).any():
            per["ref_bm"][p] = float(np.nanmedian(bm))
        if np.isfinite(pl).any():
            per["ref_plain"][p] = float(np.nanmedian(pl))
    t_ref = time.perf_counter() - t0 - t_frames

    coh = cloud.validate_member_coherence(
        emb,
        pyr,
        resolution=resolution,
        keypoint_anchor=True,
        matching="blur_matched",
        return_matrix=True,
    )
    for d in coh:
        p = int(d["point_index"])
        per["pair_plain"][p] = median_offdiag(np.asarray(d["zncc"]))
        if "blur_matched_zncc" in d:
            per["pair_bm"][p] = median_offdiag(np.asarray(d["blur_matched_zncc"]))
    t_all = time.perf_counter() - t0
    return dict(
        file=path.name,
        posed=len(names),
        points=int(len(pts)),
        framed=int(len(cloud)),
        scored=int(np.isfinite(per["ref_bm"]).sum()),
        cell_normals=int(cell_normal[pts2].sum()),
        determinacy=det_counts,
        median_views=float(np.median(n_views)) if len(n_views) else 0.0,
        unmapped_obs=int((~mapped).sum()),
        on_camera_points=int(len(np.unique(tp_f[on_camera]))),
        t_frames_s=round(t_frames, 3),
        t_ref_s=round(t_ref, 3),
        t_total_s=round(t_all, 3),
        point_cluster=[int(c) for c in point_cluster[pts2]],
        per_point={
            r: [None if not np.isfinite(v) else round(float(v), 5) for v in per[r]]
            for r in READINGS
        },
    )


def aggregate(values, bar):
    v = np.asarray([x for x in values if x is not None], dtype=np.float64)
    if len(v) == 0:
        return dict(median=None, trimmed=None, share=None, n=0)
    s = np.sort(v)
    cut = int(0.1 * len(s))
    tr = s[cut : len(s) - cut] if len(s) - 2 * cut > 0 else s
    return dict(
        median=float(np.median(v)),
        trimmed=float(np.mean(tr)),
        share=float(np.mean(v > bar)),
        n=int(len(v)),
    )


def add_aggregates(cands):
    """Capture-relative bar per reading: the median pooled over every candidate."""
    bars = {}
    for r in READINGS:
        pool = [x for c in cands for x in c["per_point"][r] if x is not None]
        bars[r] = float(np.median(pool)) if pool else float("nan")
        for c in cands:
            c.setdefault("agg", {})[r] = aggregate(c["per_point"][r], bars[r])
    return bars


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("workspace", nargs="+")
    ap.add_argument("--out")
    ap.add_argument("--resolution", type=int, default=12)
    ap.add_argument("--normal", choices=("cell", "mean_viewing"), default="cell")
    ap.add_argument(
        "--files",
        nargs="*",
        help="candidate files to score (default: every finite one)",
    )
    ap.add_argument(
        "--keep-per-point",
        action="store_true",
        help="keep per-point values in the output",
    )
    a = ap.parse_args()
    out = {}
    for ws in a.workspace:
        t0 = time.perf_counter()
        W = Workspace(ws)
        cdir = W.ws / "sfmr" / "candidate_solves"
        man = json.loads((cdir / "manifest.json").read_text(encoding="utf-8"))
        finite = {
            h.get("release_file")
            for h in man["hypotheses"]
            if h.get("model") == "finite"
        }
        files = sorted(
            p
            for p in cdir.glob("*.sfmr")
            if (p.name in finite or p.name.replace("-spline-refused", "") in finite)
            and (not a.files or p.name in a.files)
        )
        cands = []
        for p in files:
            try:
                c = score_candidate(W, p, a.resolution, a.normal)
            except ValueError as err:
                print(f"{W.ws.name:30s} {p.name:26s} not scored: {err}", flush=True)
                continue
            cands.append(c)
            print(
                f"{W.ws.name:30s} {c['file']:26s} posed {c['posed']:3d} pts {c['points']:5d} "
                f"scored {c['scored']:5d} cellN {c['cell_normals']:5d} {c['t_total_s']:6.1f} s",
                flush=True,
            )
        bars = add_aggregates(cands)
        for c in cands:
            print(
                f"   {c['file']:26s} "
                + "  ".join(
                    f"{r} {c['agg'][r]['median'] if c['agg'][r]['median'] is not None else float('nan'):.3f}"
                    f"/{c['agg'][r]['share'] if c['agg'][r]['share'] is not None else float('nan'):.2f}"
                    for r in READINGS
                )
            )
            if not a.keep_per_point:
                c.pop("per_point")
                c.pop("point_cluster")
        out[str(W.ws)] = dict(
            bars=bars,
            candidates=cands,
            wall_s=round(time.perf_counter() - t0, 1),
            resolution=a.resolution,
            normal=a.normal,
        )
    if a.out:
        Path(a.out).write_text(json.dumps(out, indent=1), encoding="utf-8")
        print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
