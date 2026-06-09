# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 18 — reciprocity (mutual k-NN) as a membership criterion.

exp17 showed an absolute *radius* can't separate co-observations from background
for many in-track descriptors. A radius is the wrong test: a true co-observation
pair are each other's near neighbours, whereas a background feature that happens
to sit near point A usually has its *own* nearest neighbours elsewhere. So test a
*relative* criterion — reciprocity: edge i→j is mutual iff i is also among j's R
nearest neighbours.

Using the (independent) reference solve's tracks as ground truth, we label each
cross-image neighbour edge of an in-track descriptor as co-observation or
background and report, for the radius baseline (Otsu d1) and for mutual-within-R
(R = 1,2,4,8,16), the co-observation recall vs background admitted. A good
criterion keeps co-obs and drops background.

Usage:
    pixi run -e experiments python experiments/exp18_reciprocity.py
"""

from __future__ import annotations

import argparse
import glob

import numpy as np

from exp03_radius_clusters import K
from exp05_cluster_match import derive_threshold
from sfm_descriptors import load_descriptor_bank

DATASETS = [
    "seoul_bull_ws",
    "seattle_backyard_ws",
    "kerry_park_ws",
    "dino_dog_toy_ws",
]
RS = [1, 2, 4, 8, 16]


def run_one(ws: str, preset: str) -> None:
    from sfmtool import KdForest

    path = sorted(glob.glob(f"../{ws}/sfmr/*solve*.sfmr"))[0]
    bank = load_descriptor_bank(path)
    desc = np.ascontiguousarray(bank.descriptors)
    idx, dst = KdForest(desc, preset=preset).query(desc, k=K)
    idx = idx.astype(np.int64)
    n = bank.n
    img = bank.image_label.astype(np.int64)
    point = bank.point_label
    pos = point >= 0
    T = derive_threshold(dst[:, 1], "otsu")

    # All directed edges i -> j over the 16 neighbours (drop self col 0).
    i = np.repeat(np.arange(n), K - 1)
    j = idx[:, 1:].ravel()
    d = dst[:, 1:].ravel()
    rank = np.tile(np.arange(1, K), n)  # neighbour rank of j in i's list (1..16)

    cross = img[i] != img[j]
    src = pos[i] & cross  # in-track source, cross-image edge
    same = (point[j] == point[i]) & (point[i] >= 0)
    co = src & same
    bg = src & ~same

    # Reverse-rank: at what rank does i appear in j's neighbour list? (0=not there)
    key = i * n + j
    rev = j * n + i
    order = np.argsort(key, kind="stable")
    ks = key[order]
    # position of each reverse edge among forward edges (if present)
    pos_in = np.searchsorted(ks, rev)
    pos_in = np.clip(pos_in, 0, len(ks) - 1)
    found = ks[pos_in] == rev
    rev_rank = np.zeros(len(i), dtype=np.int64)
    rev_rank[found] = rank[order][pos_in[found]]  # rank of i in j's list

    print(
        f"\n{ws.replace('_ws', '')}: n={n}  Otsu T={T:.0f}  "
        f"co-obs edges={int(co.sum())} bg edges={int(bg.sum())}"
    )
    print(f"  {'criterion':<16} {'co-obs kept':>12} {'bg admitted':>12}")
    keep = d <= T
    print(f"  {'radius (Otsu)':<16} {keep[co].mean():>11.1%} {keep[bg].mean():>12.1%}")
    for r in RS:
        mutual = found & (rev_rank <= r)
        print(
            f"  {'mutual R≤' + str(r):<16} {mutual[co].mean():>11.1%} "
            f"{mutual[bg].mean():>12.1%}"
        )
    # combined: mutual (R≤16) AND within radius
    both = found & keep
    print(
        f"  {'mutual & radius':<16} {both[co].mean():>11.1%} {both[bg].mean():>12.1%}"
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--preset", default="accurate")
    ap.add_argument("--datasets", nargs="*", default=DATASETS)
    args = ap.parse_args()
    print(
        "Reciprocity vs radius for separating co-obs from background "
        "(reference-solve tracks)."
    )
    for ws in args.datasets:
        run_one(ws, args.preset)


if __name__ == "__main__":
    main()
