# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 20 — extract labelled local-neighbourhood test cases.

Saves, per dataset, a fixture of seed neighbourhoods with everything a local
cluster-determination strategy needs, so candidate rules can be benchmarked
(exp21) in milliseconds without re-querying the index or a solve.

For each sampled seed we store its `KB` nearest neighbours (a *wider*
neighbourhood than the matcher's 16) with, per neighbour:

  - dist        : descriptor L2 distance from the seed (sorted ascending)
  - image       : the neighbour's image index
  - point       : its track id under the solve (-1 if untracked)
  - is_coobs    : true co-observation = same track, different image
  - same_image  : neighbour shares the seed's image (never a valid member)
  - rev_rank    : rank of the seed within this neighbour's own KB list
                  (1 = neighbour's nearest; 0 = seed not in its list) — for
                  reciprocity/mutual-NN strategies.

Ground truth is the solve named in `sfm_descriptors.SOLVES` (dino_large included).
Fixtures are written to experiments/cluster_cases/<name>.npz.

Usage:
    pixi run -e experiments python experiments/exp20_make_cluster_cases.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from sfm_descriptors import SOLVES, load_descriptor_bank, resolve_solve

KB = 48  # neighbours kept per seed (a wider neighbourhood than the matcher's 16)


def extract(path: str, name: str, preset: str, n_in: int, n_bg: int, seed: int):
    from sfmtool import KdForest

    bank = load_descriptor_bank(path)
    desc = np.ascontiguousarray(bank.descriptors)
    idx, dst = KdForest(desc, preset=preset).query(desc, k=KB + 1)
    idx = idx.astype(np.int64)
    img = bank.image_label.astype(np.int32)
    point = bank.point_label.astype(np.int64)
    pos = point >= 0

    nb = idx[:, 1:]  # (n, KB) neighbour ids
    nb_img = img[nb]
    nb_pt = point[nb]
    cross = nb_img != img[:, None]
    coobs = cross & (nb_pt == point[:, None]) & pos[:, None]

    rng = np.random.default_rng(seed)
    in_pool = np.flatnonzero(pos & coobs.any(1))  # in-track seeds with ≥1 co-obs
    bg_pool = np.flatnonzero(~pos)  # untracked seeds (right answer: no cluster)
    sel_in = rng.choice(in_pool, min(n_in, in_pool.size), replace=False)
    sel_bg = rng.choice(bg_pool, min(n_bg, bg_pool.size), replace=False)
    sel = np.concatenate([sel_in, sel_bg])
    rng.shuffle(sel)

    # rev_rank: where does the seed appear in each neighbour's own KB list?
    nbn = idx[nb[sel], 1:]  # (S, KB, KB) each selected neighbour's neighbours
    match = nbn == sel[:, None, None]
    has = match.any(2)
    rev_rank = np.where(has, match.argmax(2) + 1, 0).astype(np.int16)

    return dict(
        dataset=name,
        kb=KB,
        solve=path,
        seed_global=sel.astype(np.int64),
        seed_image=img[sel],
        seed_point=point[sel],
        nb_global=nb[sel].astype(np.int64),
        dist=dst[sel, 1:].astype(np.float32),
        image=nb_img[sel],
        point=nb_pt[sel].astype(np.int64),
        is_coobs=coobs[sel],
        same_image=~cross[sel],
        rev_rank=rev_rank,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--preset", default="accurate")
    ap.add_argument("--datasets", nargs="*", default=[n for n, _ in SOLVES])
    ap.add_argument("--in-seeds", type=int, default=2500)
    ap.add_argument("--bg-seeds", type=int, default=500)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", default="cluster_cases")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    for name, pathglob in SOLVES:
        if name not in args.datasets:
            continue
        case = extract(
            resolve_solve(pathglob),
            name,
            args.preset,
            args.in_seeds,
            args.bg_seeds,
            args.seed,
        )
        fp = out / f"{name}.npz"
        np.savez_compressed(fp, **case)
        S = len(case["seed_global"])
        n_co = int(case["is_coobs"].sum())
        print(
            f"{name}: {S} seeds, KB={KB}, {n_co} co-obs labels "
            f"({n_co / S:.1f}/seed avg)  -> {fp} ({fp.stat().st_size // 1024} KB)"
        )


if __name__ == "__main__":
    main()
