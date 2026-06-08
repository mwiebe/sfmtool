# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 08 — proper mean-shift refinement of clusters, measured carefully.

The medoid re-centring in exp07 barely moved (the medoid shares the seed's
neighbour list), so it was effectively a no-op. This implements the real
centroid step: re-query the index at the cluster's *mean* descriptor (a synthetic
point), which can reach members closer to the true centroid than to any single
descriptor. We compare seed-only vs mean-shift refinement on the actual cluster
differences — overall recovery, recovery on large tracks, coverage, purity, and
how many clusters/members actually changed.

Usage:
    pixi run -e experiments python experiments/exp08_refine_meanshift.py seoul_bull_ws
"""

from __future__ import annotations

import argparse
import glob

import numpy as np

from exp05_cluster_match import K, derive_threshold
from sfm_descriptors import load_descriptor_bank


def _resolve_one_per_image(cand, cand_d, img, claimed, T):
    """Within-T, unclaimed candidates → one per image (nearest first)."""
    m = (cand_d <= T) & (~claimed[cand])
    ci, cd = cand[m], cand_d[m]
    if ci.size == 0:
        return ci
    o = np.lexsort((cd, img[ci]))
    ii = img[ci][o]
    keep = np.empty(len(o), dtype=bool)
    keep[0] = True
    keep[1:] = ii[1:] != ii[:-1]
    return ci[o][keep]


def cluster(forest, desc, img, idx, dst, T, refine_iters):
    """Density-seeded clustering; refine_iters>0 = mean-shift re-query."""
    n = len(desc)
    rng = np.arange(n)
    valid = (dst <= T) & (idx != rng[:, None]) & (img[idx] != img[:, None])
    order = np.argsort(-valid.sum(1), kind="stable")
    claimed = np.zeros(n, dtype=bool)
    labels = np.full(n, -1, dtype=np.int64)
    moved = changed = nseed = 0
    cid = 0
    for s in order:
        if claimed[s]:
            continue
        # initial gather around the seed (cross-image, one per image)
        m0 = valid[s] & ~claimed[idx[s]]
        cand = np.concatenate(([s], idx[s][m0]))
        cdd = np.concatenate(([0.0], dst[s][m0]))
        members = _resolve_one_per_image(cand, cdd, img, claimed, T)
        before = frozenset(members.tolist())
        did_move = False
        for _ in range(refine_iters):
            if members.size < 2:
                break
            mean = desc[members].astype(np.float32).mean(0)
            mu8 = np.clip(np.round(mean), 0, 255).astype(np.uint8)
            qi, qd = forest.query(mu8[None, :], k=K)
            new = _resolve_one_per_image(
                qi[0].astype(np.int64), qd[0], img, claimed, T
            )
            if frozenset(new.tolist()) == frozenset(members.tolist()):
                break
            members = new
            did_move = True
        if members.size >= 2:
            nseed += 1
            moved += int(did_move)
            changed += int(before != frozenset(members.tolist()))
            claimed[members] = True
            labels[members] = cid
            cid += 1
        else:
            claimed[s] = True
    return labels, dict(clusters=nseed, moved=moved, changed=changed)


def score(bank, labels):
    """Overall + large-track recovery, coverage, purity, false-merge."""
    point = bank.point_label
    in_track = point >= 0
    sizes = np.bincount(labels[labels >= 0], minlength=1) if (labels >= 0).any() else np.array([0])

    # purity / false-merge over in-track members in real clusters
    m = in_track & (labels >= 0)
    P = int(point.max()) + 1
    key = labels[m].astype(np.int64) * P + point[m]
    _, cnt = np.unique(key, return_counts=True)
    comp_of = (np.unique(key) // P).astype(np.int64)
    cu, tot, dom = _group(comp_of, cnt)
    dom_cov = dom.sum() / tot.sum()
    false_merge = float(np.mean(dom < tot))

    # recovery per track: max members in one cluster / total active members
    mt = in_track & (labels >= 0)
    Cn = int(labels.max()) + 1
    keyp = point[mt].astype(np.int64) * Cn + labels[mt]
    _, cntp = np.unique(keyp, return_counts=True)
    pt_of = (np.unique(keyp) // Cn).astype(np.int64)
    pu, tmem, maxc = _group(pt_of, cntp)
    rec = maxc / tmem
    overall = float(rec[tmem >= 2].mean())
    # large tracks: total *solve* track size >= 4
    full = np.bincount(point[in_track])
    large_mask = full[pu] >= 4
    large = float(rec[large_mask & (tmem >= 2)].mean()) if (large_mask & (tmem >= 2)).any() else float("nan")
    coverage = int((in_track & (labels >= 0)).sum()) / int(in_track.sum())
    return dict(dom_cov=dom_cov, false_merge=false_merge, recovery=overall,
                recovery_large=large, coverage=coverage,
                clusters=int((sizes >= 2).sum()))


def _group(keys, vals):
    o = np.argsort(keys, kind="stable")
    k = keys[o]
    v = vals[o]
    b = np.r_[0, np.flatnonzero(np.diff(k)) + 1]
    return k[b], np.add.reduceat(v, b), np.maximum.reduceat(v, b)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("workspace")
    ap.add_argument("--iters", type=int, default=3)
    args = ap.parse_args()
    from sfmtool import KdForest

    path = sorted(glob.glob(f"../{args.workspace}/sfmr/*solve*.sfmr"))[0]
    bank = load_descriptor_bank(path)
    desc = np.ascontiguousarray(bank.descriptors)
    idx, dst = KdForest(desc, preset="accurate").query(desc, k=K)
    idx = idx.astype(np.int64)
    forest = KdForest(desc, preset="accurate")
    T = derive_threshold(dst[:, 1], "otsu")
    img = bank.image_label

    print(f"{args.workspace}: n={bank.n} T={T:.0f}")
    print(f"  {'variant':<16} {'clusters':>8} {'moved':>6} {'changed':>7} "
          f"{'purity':>7} {'falseM':>7} {'recov':>6} {'recovL':>7} {'cover':>6}")
    for name, it in [("seed-only", 0), (f"mean-shift×{args.iters}", args.iters)]:
        labels, st = cluster(forest, desc, img, idx, dst, T, it)
        sc = score(bank, labels)
        print(f"  {name:<16} {sc['clusters']:>8} {st['moved']:>6} "
              f"{st['changed']:>7} {sc['dom_cov']:>7.3f} {sc['false_merge']:>7.1%} "
              f"{sc['recovery']:>6.3f} {sc['recovery_large']:>7.3f} "
              f"{sc['coverage']:>6.1%}")


if __name__ == "__main__":
    main()
