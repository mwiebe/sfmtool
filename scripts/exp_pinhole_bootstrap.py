# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment: pinhole-only coarse reconstruction from cluster patches.

Starting from a workspace holding images and a `*-clusters-patches.matches`
file (sift extraction -> cluster matching -> cluster-patches), and using only
a pinhole camera model, bootstrap a coarse 3D reconstruction and write it to
a `.sfmr` file — no COLMAP solver involved.

Pipeline:
  1. Load patch clusters; refined member positions are read directly from
     the stored affines' last column (`member_affines[k][:, 2]` holds the
     absolute keypoint position since .matches format version 4), and the
     image dimensions from the images section — no per-image .sift reads.
  2. Group images by cluster covisibility (shared-cluster counts) — no
     sequence order is assumed.  Affine (weak-perspective) ALS factorization
     of candidate seed groups (a single global factorization breaks on wide
     baselines) + Tomasi–Kanade metric upgrade, both reflection hypotheses.
  3. Seed a perspective solve on the best group (a small fixed-focal BA
     also resolves the reflection), then grow incrementally: the
     next-best-view image (most observations of valid points) is resected
     pose-only against the global structure (trimmed iterations, most-
     covisible posed poses as inits), new clusters are triangulated as
     they gain posed views, short global BAs run every few images.
  4. Steps 2–3 run per candidate focal on a small grid with f held FIXED —
     the focal is unobservable from a weak init (the residual decreases
     monotonically toward the affine limit), but with a converged geometry
     the inlier fraction peaks near the true focal.  The scan caps growth
     at ~20 images; the winner grows fully and its BA then releases f.
  5. Report reprojection stats and, when a reference solve exists in the
     workspace, camera errors after similarity alignment; save the result
     as `sfmr/bootstrap-pinhole.sfmr`.

Run: pixi run -e dev python scripts/exp_pinhole_bootstrap.py <workspace> [ref.sfmr]

The optional second argument names the reference solve to compare against
(it may live in another workspace, e.g. a full-sequence solve when
bootstrapping a frame subset — images are matched by workspace-relative
name).  Default: the first non-bootstrap .sfmr in the workspace.
"""

import itertools
import os
import sys
import time
from pathlib import Path

import numpy as np
import scipy.optimize
import scipy.sparse
from scipy.spatial.transform import Rotation

WS = Path(sys.argv[1] if len(sys.argv) > 1 else "e_seoul_ws")
REF = Path(sys.argv[2]) if len(sys.argv) > 2 else None
_T0 = time.perf_counter()
MIN_SPAN_BA = 2  # min distinct images for a cluster to become a point
MAX_CLUSTERS = 10000  # cap for the scipy BAs
F_GRID = [0.55, 0.7, 0.9, 1.2, 1.6]  # focal candidates, in units of max(w, h)
TRIM_PX = 4.0  # BA inter-round observation trim threshold
# Cluster ordering for the cap and the admission tiers: "cons" ranks by the
# stored warp-consistency residual (max over members, ascending — measured
# AUC 0.79-0.92 for junk prediction across the campaign datasets), "span" is
# the original highest-span-first ordering.
ORDER = os.environ.get("SFMTOOL_ORDER", "span")
# Coarse-to-fine admission: comma-separated cumulative fractions of the
# admitted (quality-ordered) clusters.  "1.0" = everything at once (the
# original behavior); "0.35,1.0" seeds, scans, and grows on the best 35%
# and then admits the rest for the fine-tune BAs.
TIERS = [float(x) for x in os.environ.get("SFMTOOL_TIERS", "1.0").split(",")]
# Resection init from warp-determinant depth ratios: each member warp's
# sqrt|det| predicts the point's depth in the new image from its depth in
# the (posed) reference image, giving camera-frame 3D points -> closed-form
# trimmed Kabsch pose init (no neighbor-pose inits needed when it works).
DEPTH_INIT = os.environ.get("SFMTOOL_DEPTH_INIT", "0") == "1"
# Diagnostics: trace per-resection inliers in growth; optionally disable the
# periodic growth BA to attribute damage between resection and BA.
TRACE = os.environ.get("SFMTOOL_TRACE", "0") == "1"
GROW_BA = os.environ.get("SFMTOOL_GROW_BA", "1") == "1"


# ── Data loading ─────────────────────────────────────────────────────────────


def load_clusters():
    """Patch clusters as flat observation arrays with refined positions.

    Everything geometric comes straight from the .matches file: image
    dimensions from the images section and member positions from the stored
    affines' last column (the absolute refined keypoint position).
    """
    from sfmtool._sfmtool.io import read_matches

    patches = sorted(WS.glob("matches/*-clusters-patches.matches"))
    data = read_matches(patches[0])
    names = list(data["image_names"])
    dims = [(int(w), int(h)) for w, h in np.asarray(data["image_dims"])]

    starts = np.asarray(data["cluster_starts"])
    mi = np.asarray(data["member_images"])
    mf = np.asarray(data["member_features"])
    st = np.asarray(data["member_status"])
    refs = np.asarray(data["reference_members"])
    aff = np.asarray(data["member_affines"])
    cons = np.asarray(data["member_consistency_residual"], dtype=np.float64)

    # First pass: member selections, spans, and quality of every usable
    # cluster.  Quality is the worst (max) finite warp-consistency residual
    # over the selected members — lower is better; clusters where no member
    # entered the consistency fit rank last (inf).
    usable = []
    for c in range(len(starts) - 1):
        lo, hi = int(starts[c]), int(starts[c + 1])
        if refs[c] == np.iinfo(np.uint32).max:
            continue
        sel = np.nonzero((st[lo:hi] == 0) | (st[lo:hi] == 1))[0] + lo
        span = len(np.unique(mi[sel]))
        if span >= MIN_SPAN_BA:
            cq = cons[sel]
            cq = cq[np.isfinite(cq)]
            quality = float(cq.max()) if len(cq) else np.inf
            usable.append((span, c, sel, quality))

    # Admission order (best first) — used for both the cap and the tiers.
    # "span": highest span first (the original).  "cons": best consistency
    # first.  "cons_strat": best consistency first WITHIN each span stratum,
    # strata interleaved proportionally — any admission prefix then keeps
    # the span distribution (wide-baseline rigidity) while dropping the
    # worst-consistency clusters of every stratum first.
    spans = np.array([t[0] for t in usable])
    quals = np.array([t[3] for t in usable])
    cids = np.array([t[1] for t in usable])
    if ORDER == "span":
        order = np.lexsort((cids, -spans))
    elif ORDER == "cons":
        order = np.lexsort((cids, -spans, quals))
    elif ORDER == "cons_strat":
        p = np.empty(len(usable))
        for s in np.unique(spans):
            idx = np.nonzero(spans == s)[0]
            r = np.argsort(np.argsort(quals[idx], kind="stable"))
            p[idx] = (r + 0.5) / len(idx)
        order = np.lexsort((cids, quals, p))
    elif ORDER == "cons_rr":
        # Per-image round-robin by quality: every image repeatedly claims
        # its best not-yet-claimed cluster, so any admission prefix gives
        # every image its locally-best clusters (balanced coverage — a
        # global quality or stratified cap can disconnect a chain-shaped
        # capture: south-building fragmented at 36/128 under cons_strat).
        by_img = {}
        for k, (_span, _c, sel, _q) in enumerate(usable):
            for im in np.unique(mi[sel]):
                by_img.setdefault(int(im), []).append(k)
        for im in by_img:
            by_img[im].sort(key=lambda k: (quals[k], cids[k]))
        ptr = dict.fromkeys(by_img, 0)
        claimed = np.zeros(len(usable), bool)
        order = []
        img_ids = sorted(by_img)
        while len(order) < len(usable):
            progress = False
            for im in img_ids:
                lst = by_img[im]
                p_i = ptr[im]
                while p_i < len(lst) and claimed[lst[p_i]]:
                    p_i += 1
                ptr[im] = p_i
                if p_i < len(lst):
                    claimed[lst[p_i]] = True
                    order.append(lst[p_i])
                    ptr[im] = p_i + 1
                    progress = True
            if not progress:
                break
        order = np.asarray(order, dtype=np.int64)
    else:
        raise SystemExit(f"unknown SFMTOOL_ORDER {ORDER!r}")

    # No admission cap: growth and triangulation see every usable cluster
    # (a capped set can disconnect a chain-shaped capture — south-building
    # fragmented at 36/128).  The ordering instead selects which clusters'
    # observations enter the BAs (the top MAX_CLUSTERS by adm_rank).
    pos = {int(k): i for i, k in enumerate(order)}
    keep_idx = sorted(range(len(usable)), key=lambda k: usable[k][1])

    obs_c, obs_i, obs_f, obs_uv, obs_warp, obs_ref = [], [], [], [], [], []
    adm_rank = []
    n_cl = 0
    for k in keep_idx:
        _span, c, sel, q = usable[k]
        adm_rank.append(pos[int(k)])
        for k in sel:
            obs_c.append(n_cl)
            obs_i.append(int(mi[k]))
            obs_f.append(int(mf[k]))
            # The affine's last column is the member's absolute refined
            # keypoint position (identity | x_ref for the reference row);
            # the 2x2 block is the member<-reference patch warp.
            obs_uv.append(aff[k, :, 2])
            obs_warp.append(aff[k, :, :2])
            obs_ref.append(st[k] == 0)
        n_cl += 1

    return {
        "names": names,
        "dims": dims,
        "obs_c": np.asarray(obs_c),
        "obs_i": np.asarray(obs_i),
        "obs_f": np.asarray(obs_f),
        "obs_uv": np.asarray(obs_uv, dtype=np.float64),
        "obs_warp": np.asarray(obs_warp, dtype=np.float64),
        "obs_ref": np.asarray(obs_ref, dtype=bool),
        "adm_rank": np.asarray(adm_rank, dtype=np.int64),
        "refine_radius": float(data["refine_options"]["radius"]),
        "n_img": len(names),
        "n_cl": n_cl,
    }


# ── Covisibility grouping ────────────────────────────────────────────────────
#
# No sequence order is assumed: the natural grouping is how many clusters a
# pair of images shares.  High mutual covisibility implies nearby viewpoints,
# which is exactly what the weak-perspective factorization needs from a seed
# group, and the same counts drive the growth order and the resection inits.
# The counting and grouping live in the ClusterCovisibility binding; it is
# built from the loaded (span-filtered, capped) observation arrays rather
# than from the file so it sees exactly the clusters the bootstrap uses.


def build_covisibility(obs_c, obs_i, n_img, n_cl):
    """ClusterCovisibility over the loaded observation arrays."""
    from sfmtool._sfmtool.matching import ClusterCovisibility

    # obs_c is grouped by cluster in ascending order — derive the CSR starts.
    starts = np.searchsorted(obs_c, np.arange(n_cl + 1)).astype(np.uint32)
    return ClusterCovisibility.from_arrays(starts, obs_i.astype(np.uint32), n_img)


def window_spans(obs_c, obs_i, u, imgs, min_span):
    """Observation selection for clusters seen in >= min_span window images."""
    inw = np.isin(obs_i, imgs)
    cl, il = obs_c[inw], np.searchsorted(imgs, obs_i[inw])
    span = np.zeros(cl.max() + 1 if len(cl) else 1, int)
    for c in np.unique(cl):
        span[c] = len(np.unique(il[cl == c]))
    good_cl = np.nonzero(span >= min_span)[0]
    sel = inw.copy()
    sel[inw] = np.isin(cl, good_cl)
    uniq, c2 = np.unique(obs_c[sel], return_inverse=True)
    return sel, np.searchsorted(imgs, obs_i[sel]), uniq, c2


def factorize_window(obs_c, obs_i, u, imgs, min_span=3):
    """Factorize the clusters seen in >= min_span images of the window.

    Uses the affine-factorization bindings (ALS + Tomasi–Kanade metric
    upgrade).  Returns (both metric hypotheses as (rot, scale, t_aff) in
    the window's local frame, used mask, span-2 selection for the window
    mini-BA) or None when the window is too sparse.
    """
    from sfmtool._sfmtool.geometry import factorize_affine

    sel, il, uniq, c2 = window_spans(obs_c, obs_i, u, imgs, min_span)
    if sel.sum() < 30:
        return None
    fac = factorize_affine(
        c2.astype(np.uint32),
        il.astype(np.uint32),
        np.ascontiguousarray(u[sel]),
        len(imgs),
        len(uniq),
    )
    used = np.asarray(fac.used_images)
    t_aff = np.asarray(fac.translations)
    upgraded = fac.metric_upgrade()
    hyps = []
    if upgraded is not None:
        for hyp in upgraded:
            rot = np.asarray(hyp.rotations)
            scale = np.asarray(hyp.scales)
            if (scale[used] > 0).all():
                hyps.append((rot, scale, t_aff))
    ba_sel = window_spans(obs_c, obs_i, u, imgs, 2)
    return hyps, used, ba_sel


def kabsch_trimmed(x_world, x_cam, rounds=3, keep_q=0.6):
    """Rigid R, t with x_cam ~ R·x_world + t, trimmed to the best-fitting
    fraction each round (the depth predictions include junk members)."""
    m = np.ones(len(x_world), bool)
    r_fit = np.eye(3)
    t_fit = np.zeros(3)
    for _ in range(rounds):
        muw, muc = x_world[m].mean(axis=0), x_cam[m].mean(axis=0)
        h = (x_cam[m] - muc).T @ (x_world[m] - muw)
        uu, _, vt = np.linalg.svd(h)
        d = np.sign(np.linalg.det(uu @ vt))
        r_fit = uu @ np.diag([1.0, 1.0, d]) @ vt
        t_fit = muc - r_fit @ muw
        res = np.linalg.norm(x_world @ r_fit.T + t_fit - x_cam, axis=1)
        m = res <= np.quantile(res, keep_q)
    return r_fit, t_fit


# Per-image warp-depth coherence measured at resection acceptance
# (image, median |log(z_pose / z_warp_predicted)|, resection inlier frac).
_DEPTH_COH = []


def depth_init(s, obs_c, u, pts, rvec, tvec, posed, f0, i, aux):
    """Closed-form pose init for image ``i`` from warp-predicted depths.

    Each observation's sqrt|det warp| is the reference->member magnification,
    so the point's depth in image i is its depth in the (posed) reference
    image divided by it; backprojecting at those depths gives camera-frame
    points and a trimmed Kabsch solve gives the pose.  Returns (rvec0,
    tvec0, obs index array, predicted depths) or None when too few
    observations have a posed reference view."""
    ds, ref_img = aux
    si = np.nonzero(s)[0]
    rc = ref_img[obs_c[si]]
    okd = (rc >= 0) & (rc != i) & posed[np.maximum(rc, 0)]
    if okd.sum() < 8:
        return None
    x_w = pts[obs_c[si[okd]]]
    r_ref = Rotation.from_rotvec(rvec[rc[okd]]).as_matrix()
    z_ref = np.einsum("nij,nj->ni", r_ref, x_w)[:, 2] + tvec[rc[okd], 2]
    z_pred = z_ref / ds[si[okd]]
    good = z_pred > 1e-6
    if good.sum() < 8:
        return None
    sel = si[okd][good]
    x_cam = np.column_stack([u[sel] / f0, np.ones(good.sum())]) * z_pred[good, None]
    r_fit, t_fit = kabsch_trimmed(x_w[good], x_cam)
    return Rotation.from_matrix(r_fit).as_rotvec(), t_fit, sel, z_pred[good]


def pose_refine(uv, x_pts, rv0, tv0, f):
    """Pose-only resection of one image against known 3D points.

    Trimmed iterations: repeatedly refit L2 on the best-fitting 60% of the
    observations.  A plain L2 warm-up is dragged by the junk observations'
    leverage, and a robust loss has near-zero gradient when every residual
    starts as a 100 px "outlier" — trimming from a decent init has neither
    problem.  A final refit uses the < 3 px inliers."""

    def rn_of(x):
        xc = Rotation.from_rotvec(x[:3]).apply(x_pts) + x[3:]
        z = np.maximum(xc[:, 2], 1e-6)
        return np.linalg.norm(f * xc[:, :2] / z[:, None] - uv, axis=1)

    def fit(x0, mask):
        pts_sel, uv_sel = x_pts[mask], uv[mask]

        def resid(x):
            xc = Rotation.from_rotvec(x[:3]).apply(pts_sel) + x[3:]
            z = np.maximum(xc[:, 2], 1e-6)
            return (f * xc[:, :2] / z[:, None] - uv_sel).ravel()

        return scipy.optimize.least_squares(resid, x0, max_nfev=60).x

    x0 = np.concatenate([rv0, tv0])
    rn = rn_of(x0)
    for _ in range(5):
        x0 = fit(x0, rn <= np.quantile(rn, 0.6))
        rn = rn_of(x0)
    inl = rn < 3.0
    if inl.sum() >= 6:
        x0 = fit(x0, inl)
        rn = rn_of(x0)
    return x0[:3], x0[3:], float((rn < 3.0).mean())


def fill_new_points(pts, obs_c, obs_i, u, rvec, tvec, posed, f):
    """DLT-triangulate clusters that lack a point but now have >= 2 posed
    observations.  Existing points are left untouched."""
    need = np.isnan(pts[:, 0])[obs_c] & posed[obs_i]
    if not need.any():
        return pts
    uniq, c2 = np.unique(obs_c[need], return_inverse=True)
    rot = Rotation.from_rotvec(rvec).as_matrix()
    newp = triangulate(c2, obs_i[need], u[need], rot, tvec, posed, len(uniq), f)
    out = pts.copy()
    out[uniq] = newp
    return out


def grow_reconstruction(
    grp_data, f0, obs_c, obs_i, u, n_img, n_cl, covis, max_images=None,
    aux=None, ba=None,
):
    """Incremental bootstrap, no sequence order assumed.

    Seed a perspective solve on a candidate covisibility group (all groups x
    both reflection hypotheses, best inlier fraction wins), then repeatedly
    pose the unposed image with the most observations of valid 3D points
    (next-best-view) by trimmed pose-only resection — initialised from its
    most-covisible posed images — triangulating new clusters as they gain
    posed views, with a short global BA every few images.

    ``max_images`` caps growth (used by the focal scan, which only needs
    enough geometry to rank candidates).  Returns (rvec, tvec, pts, posed)
    or None.
    """
    grow_schedule = [(30.0, 3.0), (8.0, 1.5)]
    ba_every = max(3, min(8, n_img // 10))

    best = None
    for imgs, wd in grp_data:
        if wd is None:
            continue
        hyps, used, (sel, il, cl_ids, c2) = wd
        for rot0, scale, t_aff in hyps:
            trans0 = perspective_init(rot0, scale, t_aff, used, f0)
            pts_w = triangulate(c2, il, u[sel], rot0, trans0, used, len(cl_ids), f0)
            ok = ~np.isnan(pts_w[:, 0])[c2] & used[il]
            _, rvw, tvw, p_w, _, _, inl = bundle_adjust(
                c2[ok],
                il[ok],
                u[sel][ok],
                rot0,
                trans0,
                pts_w,
                f0,
                len(imgs),
                len(cl_ids),
                opt_f=False,
                verbose=False,
            )
            if best is None or inl > best[0]:
                best = (inl, imgs, used, cl_ids, rvw, tvw, p_w)
    if best is None:
        return None
    _, imgs, used, cl_ids, rvw, tvw, p_w = best

    rvec = np.zeros((n_img, 3))
    tvec = np.tile([0.0, 0.0, f0], (n_img, 1))
    posed = np.zeros(n_img, bool)
    pts = np.full((n_cl, 3), np.nan)
    for k, i in enumerate(imgs):
        if used[k]:
            rvec[i], tvec[i], posed[i] = rvw[k], tvw[k], True
    pts[cl_ids] = p_w

    return grow_loop(
        rvec, tvec, pts, posed, f0, obs_c, obs_i, u, n_img, n_cl, covis,
        max_images, aux, ba,
    )


def grow_loop(
    rvec, tvec, pts, posed, f0, obs_c, obs_i, u, n_img, n_cl, covis,
    max_images=None, aux=None, ba=None,
):
    """Next-best-view growth from an existing state (resumable: tier
    admission re-enters here after activating more clusters)."""
    grow_schedule = [(30.0, 3.0), (8.0, 1.5)]
    ba_every = max(3, min(8, n_img // 10))

    def run_grow_ba(rvec, tvec, pts):
        live = posed[obs_i] & ~np.isnan(pts[obs_c, 0])
        if ba is not None:
            live &= ba
        rot = Rotation.from_rotvec(rvec).as_matrix()
        out = bundle_adjust(
            obs_c[live], obs_i[live], u[live], rot, tvec, pts, f0,
            n_img, n_cl, opt_f=False, verbose=False, schedule=grow_schedule,
        )
        return out[1], out[2], out[3]

    def image_inl(i, rvec, tvec, pts):
        s = (obs_i == i) & ~np.isnan(pts[obs_c, 0])
        if not s.any():
            return 0.0
        xc = Rotation.from_rotvec(rvec[i]).apply(pts[obs_c[s]]) + tvec[i]
        z = np.maximum(xc[:, 2], 1e-6)
        rn = np.linalg.norm(f0 * xc[:, :2] / z[:, None] - u[s], axis=1)
        return float((rn < 3.0).mean())

    since_ba = 0
    accepted_inl = []
    blocked = set()
    force_tried = set()
    ba_retry = True
    while max_images is None or posed.sum() < max_images:
        # Next-best-view: most observations of currently-valid points.
        cand = ~posed[obs_i] & ~np.isnan(pts[obs_c, 0])
        if not cand.any():
            break
        cnt = np.bincount(obs_i[cand], minlength=n_img)
        cnt_all = cnt.copy()
        for j in blocked:
            cnt[j] = 0
        i = int(np.argmax(cnt))
        if cnt[i] < 6:
            # Every eligible image is blocked or too weak.  One BA +
            # retriangulation pass may repair the frontier; afterwards the
            # blocked images get a second chance.  (Ranking-only scan
            # growth skips the retry like it skips force-accept: it does
            # not need completion and each retry costs a BA.)
            if blocked and ba_retry and max_images is None:
                ba_retry = False
                blocked.clear()
                rvec, tvec, pts = run_grow_ba(rvec, tvec, pts)
                pts = fill_new_points(pts, obs_c, obs_i, u, rvec, tvec, posed, f0)
                since_ba = 0
                continue
            # Verified force-accept: low-inlier resections are often
            # BA-recoverable (ungated seoul carried imgs 0-5 to <= 6°
            # final error this way).  Accept the strongest blocked
            # candidate WITHOUT building points from it, BA, then verify:
            # keep it only if its inliers rose into the accepted band,
            # else unpose it for good.  Damage is bounded to one BA whose
            # trims already suppress a single wrong camera.  Skipped in
            # capped (focal-scan) growth: the scan ranks candidates, it
            # does not need completion, and each trial costs a BA.
            if max_images is not None:
                break
            trial = [j for j in blocked if j not in force_tried and cnt_all[j] >= 6]
            if trial:
                j = max(trial, key=lambda k: cnt_all[k])
                force_tried.add(j)
                blocked.discard(j)
                sj = (obs_i == j) & ~np.isnan(pts[obs_c, 0])
                posed_idx = np.nonzero(posed)[0].astype(np.uint32)
                inits = covis.rank_by_covisibility(j, posed_idx)[:3]
                best_j = None
                for k in inits:
                    rv, tv, inl = pose_refine(
                        u[sj], pts[obs_c[sj]], rvec[k], tvec[k], f0
                    )
                    if best_j is None or inl > best_j[0]:
                        best_j = (inl, rv, tv)
                _, rvec[j], tvec[j] = best_j
                posed[j] = True
                rvec, tvec, pts = run_grow_ba(rvec, tvec, pts)
                since_ba = 0
                inl_after = image_inl(j, rvec, tvec, pts)
                bar = 0.35 * float(np.median(accepted_inl)) if accepted_inl else 0.0
                if inl_after >= bar:
                    accepted_inl.append(inl_after)
                    pts = fill_new_points(
                        pts, obs_c, obs_i, u, rvec, tvec, posed, f0
                    )
                    ba_retry = True
                    blocked.clear()
                    if TRACE:
                        print(f"    force-accept img {j}: {best_j[0]:.0%} -> "
                              f"{inl_after:.0%} after BA (kept)")
                else:
                    posed[j] = False
                    if TRACE:
                        print(f"    force-reject img {j}: {best_j[0]:.0%} -> "
                              f"{inl_after:.0%} after BA (unposed)")
                continue
            break
        s = (obs_i == i) & ~np.isnan(pts[obs_c, 0])
        # Init candidates: the warp-depth Kabsch pose (when enabled and
        # enough observations have posed reference views), then the
        # most-covisible posed images' poses.  First init clearing 40%
        # inliers wins.
        di = depth_init(s, obs_c, u, pts, rvec, tvec, posed, f0, i, aux) \
            if aux is not None and DEPTH_INIT else None
        init_poses = [] if di is None else [(di[0], di[1])]
        posed_idx = np.nonzero(posed)[0].astype(np.uint32)
        inits = covis.rank_by_covisibility(i, posed_idx)[:3]
        if len(inits) == 0:
            inits = posed_idx[:1]
        init_poses += [(rvec[j], tvec[j]) for j in inits]
        found = None
        for rv0, tv0 in init_poses:
            rv, tv, inl = pose_refine(u[s], pts[obs_c[s]], rv0, tv0, f0)
            if found is None or inl > found[0]:
                found = (inl, rv, tv)
            if inl > 0.4:
                break
        # Acceptance gate: a resection far below the accepted-so-far level
        # is a misregistration in the making (the no-gate trace showed 0-7%
        # resections cascading into an 80° wreck), but the marginal band is
        # recoverable by the periodic BAs and carries the growth chain, so
        # the bar sits well below the median (seoul full-data trace:
        # accepted 49-81%, recoverable boundary 22%, poison 0-10%).  Defer
        # the image; it gets another chance after the frontier improves.
        if accepted_inl and found[0] < 0.35 * float(np.median(accepted_inl)):
            blocked.add(i)
            if TRACE:
                print(f"    defer  img {i}: inl {found[0]:.0%} on "
                      f"{int(s.sum())} obs (median accepted "
                      f"{float(np.median(accepted_inl)):.0%})")
            continue
        accepted_inl.append(found[0])
        _, rvec[i], tvec[i] = found
        posed[i] = True
        ba_retry = True
        if TRACE:
            print(f"    resect img {i}: inl {found[0]:.0%} on {int(s.sum())} obs")
        if di is not None:
            # Warp-depth coherence of the accepted pose (echo diagnostics):
            # a misregistered camera can look reprojection-consistent while
            # its pose-implied depths disagree with the warp-predicted ones.
            _, _, sel, z_pred = di
            xc = Rotation.from_rotvec(rvec[i]).apply(pts[obs_c[sel]]) + tvec[i]
            ok_z = (xc[:, 2] > 1e-6) & (z_pred > 1e-6)
            if ok_z.sum() >= 6:
                coh = float(np.median(np.abs(np.log(xc[ok_z, 2] / z_pred[ok_z]))))
                _DEPTH_COH.append((i, coh, found[0]))
        pts = fill_new_points(pts, obs_c, obs_i, u, rvec, tvec, posed, f0)
        since_ba += 1
        if GROW_BA and since_ba >= ba_every:
            since_ba = 0
            rvec, tvec, pts = run_grow_ba(rvec, tvec, pts)
    return rvec, tvec, pts, posed


# ── Perspective conversion + triangulation ───────────────────────────────────


def perspective_init(rot, scale, t_cam, used, f0):
    """Weak-perspective -> pinhole poses: depth of the point centroid is
    f0/s, lateral offset t/s (COLMAP-style camera frame, +z forward)."""
    trans = np.zeros((len(rot), 3))
    for i in np.nonzero(used)[0]:
        trans[i] = [t_cam[i, 0] / scale[i], t_cam[i, 1] / scale[i], f0 / scale[i]]
    return trans


def triangulate(obs_c, obs_i, u, rot, trans, used, n_cl, f):
    """Ray-midpoint triangulation of every cluster from the posed images,
    via the batch triangulation binding (clusters with < 2 posed
    observations stay NaN)."""
    from sfmtool._sfmtool.analysis import triangulate_batch

    pts = np.full((n_cl, 3), np.nan)
    sel = used[obs_i]
    if not sel.any():
        return pts
    oc, oi, uv = obs_c[sel], obs_i[sel], u[sel]
    # World-space unit rays and camera centers: x_cam = R x + t, so the
    # ray is Rᵀ·normalize([u/f, 1]) and the center -Rᵀ t.
    d_loc = np.column_stack([uv / f, np.ones(len(uv))])
    d_loc /= np.linalg.norm(d_loc, axis=1, keepdims=True)
    dirs = np.einsum("nji,nj->ni", rot[oi], d_loc)
    centers = -np.einsum("nji,nj->ni", rot[oi], trans[oi])
    # obs_c is cluster-sorted, so the selection is CSR-ready.
    uniq, counts = np.unique(oc, return_counts=True)
    offsets = np.concatenate([[0], np.cumsum(counts)]).astype(np.int64)
    result = triangulate_batch(
        np.ascontiguousarray(dirs), np.ascontiguousarray(centers), offsets
    )
    good = counts >= 2
    pts[uniq[good]] = np.asarray(result["points"])[good]
    return pts


# ── Bundle adjustment ────────────────────────────────────────────────────────


def project(f, rvec, tvec, pts, obs_i, obs_c):
    xc = Rotation.from_rotvec(rvec[obs_i]).apply(pts[obs_c]) + tvec[obs_i]
    z = np.maximum(xc[:, 2], 1e-6)
    return f * xc[:, :2] / z[:, None], xc[:, 2]


def solve_round(obs_c, obs_i, u, f, rvec, tvec, pts, n_img, f_scale, opt_f):
    """One robust sparse least-squares pass over the given observations.

    Images and points are compact-indexed over what the observations touch;
    untouched poses/points pass through unchanged."""
    uniq, oc2 = np.unique(obs_c, return_inverse=True)
    uimg, oi2 = np.unique(obs_i, return_inverse=True)
    p0 = pts[uniq]
    n_pt, n_im = len(uniq), len(uimg)
    nf = 1 if opt_f else 0

    def unpack(x):
        fv = x[0] if opt_f else f
        rv = x[nf : nf + 3 * n_im].reshape(-1, 3)
        tv = x[nf + 3 * n_im : nf + 6 * n_im].reshape(-1, 3)
        pv = x[nf + 6 * n_im :].reshape(-1, 3)
        return fv, rv, tv, pv

    def resid(x):
        fv, rv, tv, pv = unpack(x)
        proj, _ = project(fv, rv, tv, pv, oi2, oc2)
        return (proj - u).ravel()

    n_obs = len(oc2)
    spar = scipy.sparse.lil_matrix(
        (2 * n_obs, nf + 6 * n_im + 3 * n_pt), dtype=np.uint8
    )
    rows = np.arange(2 * n_obs)
    if opt_f:
        spar[:, 0] = 1
    for k in range(3):
        spar[rows, nf + 3 * np.repeat(oi2, 2) + k] = 1
        spar[rows, nf + 3 * n_im + 3 * np.repeat(oi2, 2) + k] = 1
        spar[rows, nf + 6 * n_im + 3 * np.repeat(oc2, 2) + k] = 1

    x0 = np.concatenate(
        [[f] if opt_f else [], rvec[uimg].ravel(), tvec[uimg].ravel(), p0.ravel()]
    )
    sol = scipy.optimize.least_squares(
        resid,
        x0,
        jac_sparsity=spar,
        loss="soft_l1",
        f_scale=f_scale,
        max_nfev=60,
        x_scale="jac",
        verbose=0,
    )
    f, rv, tv, p_new = unpack(sol.x)
    rvec, tvec, out = rvec.copy(), tvec.copy(), pts.copy()
    rvec[uimg], tvec[uimg], out[uniq] = rv, tv, p_new
    return f, rvec, tvec, out


def bundle_adjust(
    obs_c,
    obs_i,
    u,
    rot,
    trans,
    pts,
    f0,
    n_img,
    n_cl,
    opt_f,
    verbose=True,
    schedule=None,
):
    """Staged robust BA: trim gross outliers and behind-camera observations
    before each solve, re-triangulate every cluster from the refined cameras
    between rounds (re-admitting points the bad init lost)."""
    rvec = Rotation.from_matrix(rot).as_rotvec()
    tvec = trans.copy()
    f = f0
    all_img = np.ones(n_img, bool)
    if schedule is None:
        schedule = [(50.0, 5.0), (12.0, 2.0), (TRIM_PX, 1.0)]

    for rnd, (thresh, f_scale) in enumerate(schedule):
        if rnd > 0:
            rot_now = Rotation.from_rotvec(rvec).as_matrix()
            pts = triangulate(obs_c, obs_i, u, rot_now, tvec, all_img, n_cl, f)
        proj, depth = project(f, rvec, tvec, pts, obs_i, obs_c)
        rn = np.linalg.norm(proj - u, axis=1)
        with np.errstate(invalid="ignore"):
            keep = (rn < thresh) & (depth > 1e-3 * f) & ~np.isnan(pts[obs_c, 0])
        surv = np.bincount(obs_c[keep], minlength=n_cl)
        keep &= surv[obs_c] >= MIN_SPAN_BA
        f, rvec, tvec, pts = solve_round(
            obs_c[keep],
            obs_i[keep],
            u[keep],
            f,
            rvec,
            tvec,
            pts,
            n_img,
            f_scale,
            opt_f,
        )
        proj, depth = project(f, rvec, tvec, pts, obs_i, obs_c)
        rn = np.linalg.norm(proj - u, axis=1)
        if verbose:
            print(
                f"  BA round {rnd} (trim {thresh:.0f}px): f {f:.1f}, "
                f"median reproj {np.nanmedian(rn):.2f} px on {keep.sum()} obs"
            )

    with np.errstate(invalid="ignore"):
        keep = (rn < TRIM_PX) & (depth > 1e-3 * f) & ~np.isnan(pts[obs_c, 0])
    surv = np.bincount(obs_c[keep], minlength=n_cl)
    keep &= surv[obs_c] >= MIN_SPAN_BA
    res = np.where(np.isnan(rn), np.inf, rn)
    inlier2 = float((res < 2.0).mean())
    return f, rvec, tvec, pts, keep, res, inlier2


# ── Evaluation against a reference solve ─────────────────────────────────────


def umeyama(src, dst):
    """Similarity (s, R, t) minimizing ||s·R·src + t − dst||, no reflection."""
    mu_s, mu_d = src.mean(axis=0), dst.mean(axis=0)
    sc, dc = src - mu_s, dst - mu_d
    cov = dc.T @ sc / len(src)
    uu, dd, vv = np.linalg.svd(cov)
    sgn = np.sign(np.linalg.det(uu @ vv))
    r = uu @ np.diag([1, 1, sgn]) @ vv
    var = (sc**2).sum() / len(src)
    s = (dd * [1, 1, sgn]).sum() / var
    t = mu_d - s * r @ mu_s
    return s, r, t


def compare_to_reference(names, rvec, tvec, f_est, mask=None):
    """Compare against the first non-bootstrap solve in the workspace.

    Our BA poses are COLMAP-convention; the reference ``.sfmr`` is canonical.
    Convert ours to canonical first so the per-camera rotation errors are
    meaningful (the world-frame difference is absorbed by the alignment).
    ``mask`` restricts to a subset of images (e.g. the posed ones).
    """
    if mask is not None:
        names = [n for j, n in enumerate(names) if mask[j]]
        rvec, tvec = rvec[np.asarray(mask)], tvec[np.asarray(mask)]
    if REF is not None:
        ref_files = [REF]
    else:
        ref_files = sorted(
            p for p in WS.glob("sfmr/*.sfmr") if p.name != "bootstrap-pinhole.sfmr"
        )
    if not ref_files:
        print("no reference solve found; skipping comparison")
        return
    from sfmtool._sfmtool import SfmrReconstruction
    from sfmtool.colmap.io import _colmap_poses_points_to_canonical

    q_colmap = Rotation.from_rotvec(rvec).as_quat()[:, [3, 0, 1, 2]]
    q_wxyz, t_xyz, _ = _colmap_poses_points_to_canonical(
        q_colmap, tvec, np.zeros((0, 3)), apply_world_rotation=True
    )

    ref = SfmrReconstruction.load(ref_files[0])
    ref_names = list(ref.image_names)
    common = [n for n in names if n in ref_names]
    if len(common) < 3:
        # Cross-workspace fallback: match by basename against the ref
        # directory with the most unique matches (e.g. the bootstrap's
        # frames/ against a rig reference's fisheye_left/).
        from collections import defaultdict
        from pathlib import PurePosixPath

        groups = defaultdict(dict)
        for rn in ref_names:
            pp = PurePosixPath(rn)
            groups[str(pp.parent)][pp.name] = rn
        best = {}
        for g in groups.values():
            mm = {
                n: g[PurePosixPath(n).name]
                for n in names
                if PurePosixPath(n).name in g
            }
            if len(mm) > len(best):
                best = mm
        if len(best) >= 3:
            print(f"matched {len(best)} images by basename fallback")
            names = [n if n not in best else best[n] for n in names]
            common = [best[n] for n in best]
    if len(common) < 3:
        print(f"only {len(common)} common images with {ref_files[0].name}; skipping")
        return

    def centers_rots(qs, ts, order):
        rs = Rotation.from_quat(np.asarray(qs)[order][:, [1, 2, 3, 0]]).as_matrix()
        cs = -np.einsum("nij,ni->nj", rs, np.asarray(ts)[order])
        return cs, rs

    ei = np.array([names.index(n) for n in common])
    ri = np.array([ref_names.index(n) for n in common])
    c_est, r_est = centers_rots(q_wxyz, t_xyz, ei)
    c_ref, r_ref = centers_rots(ref.quaternions_wxyz, ref.translations, ri)

    s, r_al, t_al = umeyama(c_est, c_ref)
    c_fit = (s * (r_al @ c_est.T)).T + t_al
    diam = np.max(np.linalg.norm(c_ref[:, None, :] - c_ref[None, :, :], axis=2))
    cen_err = np.linalg.norm(c_fit - c_ref, axis=1) / diam
    rot_err = Rotation.from_matrix(
        np.einsum("nij,nkj->nik", r_ref, np.einsum("nij,kj->nik", r_est, r_al))
    ).magnitude() * (180 / np.pi)

    cam0 = ref.cameras[0].to_dict()
    f_ref = cam0["parameters"].get("focal_length", cam0["parameters"].get("fx"))
    print(f"\nvs reference {ref_files[0].name} ({len(common)} common images):")
    print(
        f"  camera rotation err: mean {rot_err.mean():.2f}, "
        f"median {np.median(rot_err):.2f}, max {rot_err.max():.2f} deg; "
        f"{(rot_err > 10).sum()} cams > 10 deg"
    )
    print(
        f"  camera center err:   mean {100 * cen_err.mean():.2f}%, "
        f"median {100 * np.median(cen_err):.2f}%, "
        f"max {100 * cen_err.max():.2f}% of scene diameter"
    )
    print(
        f"  focal: bootstrap {f_est:.1f} px vs reference {f_ref:.1f} px "
        f"({cam0['model']})"
    )


# ── Save as .sfmr ────────────────────────────────────────────────────────────


def save_sfmr(data, f, rvec, tvec, pts, keep, res, out_path):
    """Write the bootstrap as an ``embedded_patches`` reconstruction.

    The bootstrap's observations are the cluster patches' *refined*
    positions, not the SIFT detections, so they are stored inline as
    ``keypoints_xy`` rather than as feature indexes into the ``.sift``
    files (which would silently resolve back to the unrefined seeds).
    """
    from sfmtool._sfmtool import SfmrReconstruction
    from sfmtool._sfmtool.geometry import CameraIntrinsics
    from sfmtool._workspace import load_workspace_config
    from sfmtool.colmap.io import (
        _build_sfmr_data_dict,
        _colmap_poses_points_to_canonical,
        _resolve_workspace_and_sift,
        build_metadata,
        finite_positions_xyzw,
    )

    out_path = Path(out_path).resolve()
    names, dims = data["names"], data["dims"]
    obs_c, obs_i, obs_f = data["obs_c"], data["obs_i"], data["obs_f"]
    w, h = dims[0]

    # Surviving points, renumbered densely; observations grouped by point.
    alive = np.nonzero(np.bincount(obs_c[keep], minlength=len(pts)) >= 2)[0]
    remap = {int(c): k for k, c in enumerate(alive)}
    order = np.argsort(obs_c[keep], kind="stable")
    ko = np.nonzero(keep)[0][order]
    ko = ko[np.isin(obs_c[ko], alive)]

    track_img = obs_i[ko]
    track_feat = obs_f[ko]
    keypoints_xy = data["obs_uv"][ko].astype(np.float32)
    point_idx = np.array([remap[int(c)] for c in obs_c[ko]])
    obs_counts = np.bincount(point_idx, minlength=len(alive))

    positions = pts[alive]
    per_point_err = np.zeros(len(alive), dtype=np.float32)
    np.add.at(per_point_err, point_idx, res[ko].astype(np.float32))
    per_point_err /= np.maximum(obs_counts, 1)

    q_colmap = Rotation.from_rotvec(rvec).as_quat()[:, [3, 0, 1, 2]]
    q_can, t_can, p_can = _colmap_poses_points_to_canonical(
        q_colmap, tvec, positions, apply_world_rotation=True
    )

    (
        workspace_dir,
        _contents,
        resolved_names,
        ft_hashes,
        sc_hashes,
        thumbnails,
    ) = _resolve_workspace_and_sift(names, WS.resolve())

    # Colors from the .sift thumbnails at the (scaled) observation position.
    colors = np.zeros((len(alive), 3), dtype=np.uint8)
    uv = data["obs_uv"][ko]
    for k in range(len(ko)):
        th = np.asarray(thumbnails[track_img[k]])
        ty = int(np.clip(uv[k, 1] * th.shape[0] / h, 0, th.shape[0] - 1))
        tx = int(np.clip(uv[k, 0] * th.shape[1] / w, 0, th.shape[1] - 1))
        colors[point_idx[k]] = th[ty, tx]

    camera = CameraIntrinsics.from_dict(
        {
            "model": "SIMPLE_PINHOLE",
            "width": w,
            "height": h,
            "parameters": {
                "focal_length": float(f),
                "principal_point_x": w / 2,
                "principal_point_y": h / 2,
            },
        }
    )

    metadata = build_metadata(
        workspace_dir=workspace_dir,
        output_path=out_path,
        workspace_config=load_workspace_config(workspace_dir),
        operation="cluster_bootstrap",
        tool_name="sfmtool",
        tool_options={"camera_model": "SIMPLE_PINHOLE", "focal_grid": F_GRID},
        image_count=len(names),
        point_count=len(alive),
        observation_count=int(obs_counts.sum()),
        camera_count=1,
    )

    sfmr_dict = _build_sfmr_data_dict(
        cameras=[camera],
        image_names=resolved_names,
        camera_indexes=np.zeros(len(names), dtype=np.uint32),
        quaternions_wxyz=q_can,
        translations_xyz=t_can,
        positions_xyzw=finite_positions_xyzw(p_can),
        colors_rgb=colors,
        reprojection_errors=per_point_err,
        track_image_indexes=track_img,
        track_feature_indexes=track_feat,
        point_indexes=point_idx,
        observation_counts=obs_counts,
        feature_tool_hashes=ft_hashes,
        sift_content_hashes=sc_hashes,
        thumbnails=thumbnails,
        metadata=metadata,
    )

    recon = SfmrReconstruction.from_data(workspace_dir, sfmr_dict)

    # ── Surfel frames copied from the cluster patches ────────────────────
    # Each member's stored 2x2 warp is the projection of the cluster's
    # common surfel into that image, so the 3D patch frame is recoverable:
    # solve J_k·B = A_k per point (J_k the projection Jacobian at the
    # point, B the 3x2 map from reference-image pixels to 3D on the surfel
    # plane; the reference row contributes J_ref·B = I), then
    # u = B·(r, 0), v = B·(0, r) with r the refinement radius in reference
    # pixels (keypoint-frame radius x the reference feature's scale).
    from sfmtool._sfmtool import PatchCloud
    from sfmtool._sfmtool.io import read_sift, read_sift_metadata
    from sfmtool.colmap.convention import world_rotate_w
    from sfmtool.sift.file import get_sift_path_for_image

    feature_scales = {}
    image_file_hashes = []
    for i, name in enumerate(names):
        sp = get_sift_path_for_image(workspace_dir / name)
        meta = read_sift_metadata(sp)["metadata"]
        image_file_hashes.append(bytes.fromhex(meta["image_file_xxh128"]))
        shapes = np.asarray(read_sift(sp)["affine_shapes"], dtype=np.float64)
        feature_scales[i] = 0.5 * (
            np.linalg.norm(shapes[:, :, 0], axis=1)
            + np.linalg.norm(shapes[:, :, 1], axis=1)
        )

    rot_all = Rotation.from_rotvec(rvec).as_matrix()
    warps = data["obs_warp"][ko]
    is_ref = data["obs_ref"][ko]
    radius_kf = data["refine_radius"]
    half_u = np.zeros((len(alive), 3), dtype=np.float64)
    half_v = np.zeros((len(alive), 3), dtype=np.float64)
    normals = np.zeros((len(alive), 3), dtype=np.float64)
    p_starts = np.searchsorted(point_idx, np.arange(len(alive) + 1))
    # The reference constraint J_ref·B = I determines B up to a 2-vector
    # b_z — the surfel's out-of-plane slope in the reference camera frame
    # (B = R_refᵀ·[(z_r/f)·I + p_r·b_z ; b_z] with p_r the normalized ref
    # coords).  Each other member contributes A_k − (z_r/f)·M2 = c_k·b_z
    # with M = J_k·R_refᵀ = [M2 | m3] and c_k = M2·p_r + m3.  The tilt is
    # exactly the depth-like weakly-observed direction, so the solve gets a
    # fronto-parallel Tikhonov prior (weight relative to the members'
    # leverage) and a hard obliquity cap; these are what the photometric
    # normal refinement later polishes.
    tan_cap = np.tan(np.radians(80.0))
    for p in range(len(alive)):
        lo, hi = int(p_starts[p]), int(p_starts[p + 1])
        refs_here = np.nonzero(is_ref[lo:hi])[0]
        if len(refs_here) == 0:
            continue  # reference member trimmed: leave the zero (no-patch) frame
        k_ref = lo + int(refs_here[0])
        i_ref = int(track_img[k_ref])
        x_pt = positions[p]
        r_ref = rot_all[i_ref]
        xc_ref = r_ref @ x_pt + tvec[i_ref]
        z_ref = max(xc_ref[2], 1e-6)
        p_r = xc_ref[:2] / z_ref
        rows, rhs = [], []
        for k in range(lo, hi):
            if k == k_ref:
                continue
            i = int(track_img[k])
            xc = rot_all[i] @ x_pt + tvec[i]
            z = max(xc[2], 1e-6)
            j_proj = (f / z) * np.array(
                [[1.0, 0.0, -xc[0] / z], [0.0, 1.0, -xc[1] / z]]
            )
            m = j_proj @ rot_all[i] @ r_ref.T
            c_k = m[:, :2] @ p_r + m[:, 2]
            resid = warps[k] - (z_ref / f) * m[:, :2]
            for j in range(2):
                rows.append([c_k[0] * (1 - j), c_k[0] * j])
                rows.append([c_k[1] * (1 - j), c_k[1] * j])
                rhs.append(resid[0, j])
                rhs.append(resid[1, j])
        if not rows:
            continue
        rows = np.asarray(rows)
        rhs = np.asarray(rhs)
        # Fronto prior: damping rows scaled to a fraction of member leverage.
        lam = 0.3 * np.sqrt((rows**2).sum() / max(len(rows), 1))
        rows = np.vstack([rows, [[lam, 0.0], [0.0, lam]]])
        rhs = np.concatenate([rhs, [0.0, 0.0]])
        b_z = np.linalg.lstsq(rows, rhs, rcond=None)[0]
        # Obliquity cap: tan(tilt) = |b_z| / (z_ref / f).
        b_norm = np.linalg.norm(b_z)
        max_bz = tan_cap * z_ref / f
        if b_norm > max_bz:
            b_z *= max_bz / b_norm
        b_map = r_ref.T @ np.vstack(
            [(z_ref / f) * np.eye(2) + np.outer(p_r, b_z), b_z[None, :]]
        )
        r_px = radius_kf * feature_scales[i_ref][int(track_feat[k_ref])]
        u3 = b_map @ np.array([r_px, 0.0])
        v3 = b_map @ np.array([0.0, r_px])
        n3 = np.cross(u3, v3)
        norm = np.linalg.norm(n3)
        if norm < 1e-12:
            continue
        n3 /= norm
        cam_c = -r_ref.T @ tvec[i_ref]
        if np.dot(n3, cam_c - x_pt) < 0:
            u3, v3, n3 = v3, u3, -n3  # keep normal = normalize(u x v), front-facing
        half_u[p], half_v[p], normals[p] = u3, v3, n3

    # COLMAP -> canonical for the direction quantities (same W as the points).
    half_u = np.asarray(world_rotate_w(half_u), dtype=np.float32)
    half_v = np.asarray(world_rotate_w(half_v), dtype=np.float32)
    normals = np.asarray(world_rotate_w(normals), dtype=np.float32)

    cloud = PatchCloud.from_halfvec_arrays(half_u, half_v, np.asarray(p_can))
    recon = recon.clone_with_changes(
        feature_source="embedded_patches",
        keypoints_xy=keypoints_xy,
        image_file_hashes=image_file_hashes,
        normals=normals,
        patches=cloud,
    )
    recon.save(out_path)
    n_patched = int(np.count_nonzero(np.linalg.norm(half_u, axis=1) > 0))
    print(
        f"\nwrote {out_path} ({len(alive)} points, {int(obs_counts.sum())} obs, "
        f"{recon.feature_source}, {n_patched} warp-derived patch frames)"
    )
    return recon


# ── Main ─────────────────────────────────────────────────────────────────────


def main():
    data = load_clusters()
    all_c, all_i, u_px = data["obs_c"], data["obs_i"], data["obs_uv"]
    n_img, n_cl = data["n_img"], data["n_cl"]
    dims = np.asarray(data["dims"], dtype=np.float64)
    half = dims[data["obs_i"]] / 2
    all_u = u_px - half
    print(
        f"{WS}: {n_img} images, {n_cl} clusters (span >= {MIN_SPAN_BA}), "
        f"{len(all_c)} observations"
    )

    # Coarse-to-fine admission: quality-ordered cumulative tiers.  The seed
    # search, focal scan, and first full growth run on tier 0 only; each
    # later tier is triangulated in from the current poses and fine-tuned
    # with a shortened BA.  If tier 0 fails to produce a seed, the next
    # tier is folded in and seeding retries (incremental "until it catches").
    qual_order = np.argsort(data["adm_rank"], kind="stable")
    bounds = sorted({max(1, int(round(fr * n_cl))) for fr in TIERS} | {n_cl})
    # BA working set: the best MAX_CLUSTERS clusters in admission order.
    # Growth, resection, and triangulation always see every cluster
    # (connectivity must not starve); only the BAs are restricted to the
    # representative subset.
    ba_cl = data["adm_rank"] < MAX_CLUSTERS
    if n_cl > MAX_CLUSTERS:
        print(f"BA set: best {MAX_CLUSTERS} of {n_cl} clusters by {ORDER}")
    # Warp-depth aux data: per-obs sqrt|det| magnification and each
    # cluster's reference image (for the depth-ratio resection init).
    ds_all = np.sqrt(np.maximum(np.abs(np.linalg.det(data["obs_warp"])), 1e-12))
    ref_img = np.full(n_cl, -1, np.int64)
    ref_img[all_c[data["obs_ref"]]] = all_i[data["obs_ref"]]
    active_cl = np.zeros(n_cl, bool)
    tier = 0
    covis = grp_data = aux = bam = None
    while tier < len(bounds):
        active_cl[qual_order[: bounds[tier]]] = True
        act = active_cl[all_c]
        obs_c, obs_i, u = all_c[act], all_i[act], all_u[act]
        aux = (ds_all[act], ref_img)
        bam = ba_cl[all_c][act]
        print(
            f"tier 0..{tier}: {int(active_cl.sum())} clusters, "
            f"{len(obs_c)} observations"
        )

        # Covisibility grouping — no sequence order assumed anywhere.
        covis = build_covisibility(obs_c, obs_i, n_img, n_cl)
        grp_data = []
        for group in itertools.islice(covis.seed_groups(), 2):
            imgs = np.asarray(group)
            wd = factorize_window(obs_c, obs_i, u, imgs)
            grp_data.append((imgs, wd))
            state = "sparse" if wd is None else f"{len(wd[2][2])} span-2 clusters"
            print(f"seed group {[int(k) for k in imgs]}: {state}")
        if any(wd is not None for _, wd in grp_data):
            break
        print("no factorizable seed group at this tier; admitting the next")
        tier += 1

    # The focal is unobservable at init quality (the residual decreases
    # monotonically toward the affine limit), so grow the reconstruction at
    # each candidate focal with f held FIXED and let the inlier fraction
    # pick.  The seed group's inlier fraction resolves the reflection
    # hypothesis.  The scan caps growth at ~20 images (enough to rank);
    # only the winner grows fully, then releases f.
    f_grid = np.array(F_GRID) * dims.max()
    scan_cap = min(n_img, 20)
    best = None
    for f_try in f_grid:
        grown = grow_reconstruction(
            grp_data,
            f_try,
            obs_c,
            obs_i,
            u,
            n_img,
            n_cl,
            covis,
            max_images=scan_cap,
            aux=aux,
            ba=bam,
        )
        if grown is None:
            continue
        g_rvec, g_tvec, pts, posed = grown
        rot = Rotation.from_rotvec(g_rvec).as_matrix()
        ok = posed[obs_i] & ~np.isnan(pts[:, 0])[obs_c] & bam
        ba = bundle_adjust(
            obs_c[ok],
            obs_i[ok],
            u[ok],
            rot,
            g_tvec,
            pts,
            f_try,
            n_img,
            n_cl,
            opt_f=False,
            verbose=False,
        )
        res = ba[5]
        # Rank on ALL observations of the posed subset: clusters that failed
        # to triangulate under this candidate count as misses.
        inl_scan = float((res < 2.0).sum() / max((posed[obs_i] & bam).sum(), 1))
        print(
            f"f={f_try:6.1f}: poses {posed.sum()}/{n_img}, "
            f"inlier<2px {100 * inl_scan:5.1f}%, "
            f"median {np.median(res[np.isfinite(res)]):6.2f} px"
        )
        if best is None or inl_scan > best[0]:
            best = (inl_scan, f_try)

    _, f_try = best
    elapsed = time.perf_counter() - _T0
    print(f"\nwinner: f = {f_try:.1f}; growing fully, then releasing f "
          f"[scan done at {elapsed:.0f}s]")
    grown = grow_reconstruction(
        grp_data, f_try, obs_c, obs_i, u, n_img, n_cl, covis, aux=aux, ba=bam
    )
    g_rvec, trans, pts, posed = grown
    rot = Rotation.from_rotvec(g_rvec).as_matrix()
    ok = posed[obs_i] & ~np.isnan(pts[:, 0])[obs_c] & bam
    f, rvec, tvec, pts, keep, res, inl = bundle_adjust(
        obs_c[ok],
        obs_i[ok],
        u[ok],
        rot,
        trans,
        pts,
        f_try,
        n_img,
        n_cl,
        opt_f=True,
    )
    print(f"[tier 0..{tier} solved at {time.perf_counter() - _T0:.0f}s: "
          f"f {f:.1f}, inlier<2px {100 * inl:.1f}% of its {int(ok.sum())} obs "
          f"({len(obs_c)} in tier), {int(posed.sum())}/{n_img} posed]")
    compare_to_reference(data["names"], rvec, tvec, f, mask=posed)

    # Fine tiers: activate the remaining quality-ordered clusters,
    # triangulate them from the current cameras, resume growth for images
    # the coarse tier could not pose, and fine-tune with a shortened BA.
    for b in [x for x in bounds if x > bounds[tier]]:
        active_cl[qual_order[:b]] = True
        act = active_cl[all_c]
        obs_c, obs_i, u = all_c[act], all_i[act], all_u[act]
        print(f"\nadmitting tier: {int(active_cl.sum())} clusters, "
              f"{len(obs_c)} observations")
        pts = fill_new_points(pts, obs_c, obs_i, u, rvec, tvec, posed, f)
        covis = build_covisibility(obs_c, obs_i, n_img, n_cl)
        aux = (ds_all[act], ref_img)
        bam = ba_cl[all_c][act]
        posed_before = int(posed.sum())
        rvec, tvec, pts, posed = grow_loop(
            rvec, tvec, pts, posed, f, obs_c, obs_i, u, n_img, n_cl, covis,
            aux=aux, ba=bam,
        )
        rot = Rotation.from_rotvec(rvec).as_matrix()
        ok = posed[obs_i] & ~np.isnan(pts[:, 0])[obs_c] & bam
        # When the tier only added observations to an already-posed set,
        # the cameras are converged and admission is gated tightly BY the
        # current solve (an admitted observation disagreeing at > 2×TRIM_PX
        # under known-good cameras is junk, not unconverged).  When the
        # tier grew NEW poses, those need the full staged schedule.
        grew = int(posed.sum()) > posed_before
        if grew:
            print(f"  [after regrowth, before fine BA: "
                  f"{int(posed.sum())}/{n_img} posed]")
            compare_to_reference(data["names"], rvec, tvec, f, mask=posed)
        f, rvec, tvec, pts, keep, res, inl = bundle_adjust(
            obs_c[ok],
            obs_i[ok],
            u[ok],
            rot,
            tvec,
            pts,
            f,
            n_img,
            n_cl,
            opt_f=True,
            verbose=False,
            schedule=None if grew else [(2 * TRIM_PX, 1.5), (TRIM_PX, 1.0)],
        )
        if grew:
            print(f"  (tier grew {int(posed.sum()) - posed_before} new poses; "
                  f"full BA schedule)")
        print(f"[tier fine-tuned at {time.perf_counter() - _T0:.0f}s: "
              f"f {f:.1f}, inlier<2px {100 * inl:.1f}% of its obs]")

    act = active_cl[all_c]
    act_idx = np.nonzero(act)[0]
    full_keep = np.zeros(len(all_c), bool)
    full_keep[act_idx[np.nonzero(ok)[0]]] = keep
    full_res = np.full(len(all_c), np.inf)
    full_res[act_idx[np.nonzero(ok)[0]]] = res
    keep, res = full_keep, full_res
    rk = res[keep]
    n_pts = len(np.unique(all_c[keep]))
    print(
        f"\nbootstrap result: f = {f:.1f} px, {n_pts} points, "
        f"{keep.sum()}/{len(all_c)} observations kept, "
        f"{int(posed.sum())}/{n_img} images posed"
    )
    print(
        f"reprojection (kept): rms {np.sqrt((rk**2).mean()):.2f} px, "
        f"median {np.median(rk):.2f} px; inlier<2px {100 * (res < 2).mean():.1f}% "
        f"of all obs"
    )

    if _DEPTH_COH:
        # Warp-depth coherence at resection time (final-growth resections
        # only appear once each; scan-phase entries repeat per focal).
        coh = np.array([c for _, c, _ in _DEPTH_COH])
        worst = sorted(_DEPTH_COH, key=lambda t: -t[1])[:5]
        print(
            f"\nwarp-depth coherence at resection ({len(coh)} resections): "
            f"median {np.median(coh):.3f}, p90 {np.percentile(coh, 90):.3f} "
            f"|log depth ratio|"
        )
        print(
            "  worst: "
            + ", ".join(f"img {i} {c:.2f} (inl {v:.0%})" for i, c, v in worst)
        )

    compare_to_reference(data["names"], rvec, tvec, f, mask=posed)

    out = WS / "sfmr" / os.environ.get("SFMTOOL_OUT", "bootstrap-pinhole.sfmr")
    save_sfmr(data, f, rvec, tvec, pts, keep, res, out)


if __name__ == "__main__":
    main()
