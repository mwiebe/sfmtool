# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 09 — start from the solve's tracks and iterate the cluster step.

Initialise each cluster to a ground-truth track (the descriptors the solve tied
to one 3-D point), then repeatedly: take the mean descriptor, re-query the index
at that mean, keep the within-`T` neighbours (one per image, nearest the mean).
This shows whether the solve's tracks are fixed points of our cluster iteration
and how they evolve — retained vs shed original members, added background
(solve-discarded) members, and contamination from other tracks (merges).

With ``--radii`` it sweeps the search radius (multiples of the Otsu `T`) and
reports the converged evolution for each, tracing the precision/recall tradeoff:
a tight radius sheds true members; a loose one absorbs background then merges.

No claiming: each track evolves independently, so we observe dynamics rather than
a partition.

Usage:
    pixi run -e experiments python experiments/exp09_track_evolution.py seoul_bull_ws
    pixi run -e experiments python experiments/exp09_track_evolution.py seoul_bull_ws \
        --radii 0.5 0.75 1.0 1.25 1.5 2.0
"""

from __future__ import annotations

import argparse
import glob

import numpy as np

from exp05_cluster_match import derive_threshold
from sfm_descriptors import load_descriptor_bank

QK = 64  # query width: large enough to hold a track + nearby candidates


def one_per_image(cand, cand_d, img, T):
    """Within-T candidates → one per image, nearest the query (mean)."""
    m = cand_d <= T
    ci, cd = cand[m], cand_d[m]
    if ci.size == 0:
        return ci
    o = np.lexsort((cd, img[ci]))
    ii = img[ci][o]
    keep = np.empty(len(o), dtype=bool)
    keep[0] = True
    keep[1:] = ii[1:] != ii[:-1]
    return ci[o][keep]


def evolve(forest, desc, img, point, tracks, track_pid, orig, T, max_iters=10):
    """Iterate the mean-shift cluster step from the tracks until stable.

    Returns (converged metrics dict, iterations used)."""
    members = [m.copy() for m in tracks]
    prev = None
    iters = 0
    for it in range(1, max_iters + 1):
        means = np.zeros((len(members), desc.shape[1]), dtype=np.uint8)
        nonempty = [i for i, m in enumerate(members) if m.size > 0]
        for i in nonempty:
            means[i] = np.clip(
                np.round(desc[members[i]].astype(np.float32).mean(0)), 0, 255
            )
        qi, qd = forest.query(means, k=QK)
        for i in nonempty:
            members[i] = one_per_image(qi[i].astype(np.int64), qd[i], img, T)
        cur = [frozenset(m.tolist()) for m in members]
        iters = it
        if cur == prev:
            break
        prev = cur

    sizes, retain, addbg, addother, merged = [], [], [], [], []
    for i, m in enumerate(members):
        cur_s = frozenset(m.tolist())
        o = orig[i]
        sizes.append(len(m))
        retain.append(len(o & cur_s) / len(o))
        added = cur_s - o
        if added:
            ap_ = point[np.fromiter(added, dtype=np.int64)]
            addbg.append(int((ap_ < 0).sum()))
            addother.append(int(((ap_ >= 0) & (ap_ != track_pid[i])).sum()))
        else:
            addbg.append(0)
            addother.append(0)
        it_members = m[point[m] >= 0]
        if it_members.size:
            vals, cnts = np.unique(point[it_members], return_counts=True)
            merged.append(int(vals[cnts.argmax()] != track_pid[i]))
        else:
            merged.append(0)
    return dict(
        size=float(np.mean(sizes)), retain=float(np.mean(retain)),
        addbg=float(np.mean(addbg)), addother=float(np.mean(addother)),
        merged=float(np.mean(merged)),
    ), iters


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("workspace")
    ap.add_argument("--radii", nargs="*", type=float,
                    default=[0.5, 0.75, 1.0, 1.25, 1.5, 2.0],
                    help="search radii as multiples of the Otsu T")
    args = ap.parse_args()
    from sfmtool import KdForest

    path = sorted(glob.glob(f"../{args.workspace}/sfmr/*solve*.sfmr"))[0]
    bank = load_descriptor_bank(path)
    desc = np.ascontiguousarray(bank.descriptors)
    img = bank.image_label
    point = bank.point_label
    forest = KdForest(desc, preset="accurate")
    T_base = derive_threshold(
        forest.query(desc[: min(bank.n, 50000)], k=2)[1][:, 1], "otsu"
    )

    in_track = np.flatnonzero(point >= 0)
    order = in_track[np.argsort(point[in_track], kind="stable")]
    pid = point[order]
    bounds = np.r_[0, np.flatnonzero(np.diff(pid)) + 1, len(order)]
    tracks = [order[s:e] for s, e in zip(bounds[:-1], bounds[1:])]
    track_pid = [int(point[m[0]]) for m in tracks]
    orig = [frozenset(m.tolist()) for m in tracks]

    print(f"{args.workspace}: n={bank.n} Otsu_T={T_base:.0f} tracks={len(tracks)} "
          f"mean_size={np.mean([len(m) for m in tracks]):.2f}")
    print(f"  {'×Otsu':>6} {'T':>5} {'size':>6} {'retain':>7} {'+bg':>6} "
          f"{'+other':>7} {'merged':>7} {'iters':>5}")
    for sc in args.radii:
        T = T_base * sc
        m, iters = evolve(forest, desc, img, point, tracks, track_pid, orig, T)
        print(f"  {sc:>6.2f} {T:>5.0f} {m['size']:>6.2f} {m['retain']:>7.1%} "
              f"{m['addbg']:>6.2f} {m['addother']:>7.3f} {m['merged']:>7.1%} "
              f"{iters:>5}")


if __name__ == "__main__":
    main()
