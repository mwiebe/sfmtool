# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 03 — bounded-radius descriptor clustering vs the solve's tracks.

Pipeline (the idea from exp02):

  1. for every descriptor, take its <=16 nearest neighbours,
  2. optionally drop "isolated" descriptors (d1/d5 ratio or d1 too large),
  3. keep edges within radius T, cross-image only,
  4. connected components = candidate tracks,
  5. score components against the solve's ground-truth tracks.

The 17-NN (16 + self) with true distances is computed once (exact, full corpus)
and cached, so the T sweep is cheap.

Usage:
    pixi run -e experiments python experiments/exp03_radius_clusters.py \
        ../dino_dog_toy_ws/sfmr/*.sfmr
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from sfm_descriptors import load_descriptor_bank

K = 17  # 16 neighbours + self


def knn_all(X: np.ndarray, k: int, batch: int = 512) -> tuple[np.ndarray, np.ndarray]:
    """Exact k nearest (incl self) with true L2 distances, for every row."""
    X = X.astype(np.float32)
    sq = np.einsum("ij,ij->i", X, X)
    n = len(X)
    idx = np.empty((n, k), dtype=np.int32)
    dst = np.empty((n, k), dtype=np.float32)
    for s in range(0, n, batch):
        qb = X[s : s + batch]
        d2 = sq[None, :] - 2.0 * (qb @ X.T) + np.einsum("ij,ij->i", qb, qb)[:, None]
        part = np.argpartition(d2, k - 1, axis=1)[:, :k]
        rows = np.arange(len(qb))[:, None]
        order = np.argsort(d2[rows, part], axis=1)
        nn = part[rows, order]
        idx[s : s + batch] = nn
        dst[s : s + batch] = np.sqrt(np.maximum(d2[rows, nn], 0.0))
    return idx, dst


def group_reduce(keys: np.ndarray, vals: np.ndarray):
    """Return (unique_keys, sum, max) over vals grouped by sorted-unique keys."""
    order = np.argsort(keys, kind="stable")
    k = keys[order]
    v = vals[order]
    bounds = np.r_[0, np.flatnonzero(np.diff(k)) + 1]
    uk = k[bounds]
    s = np.add.reduceat(v, bounds)
    m = np.maximum.reduceat(v, bounds)
    return uk, s, m


def score(bank, comp, sizes, active, T):
    """Score connected components vs the solve's tracks."""
    n = bank.n
    point = bank.point_label
    image = bank.image_label
    in_multi = sizes[comp] >= 2

    # --- component purity over in-track members ---
    m = (point >= 0) & active & in_multi
    res = {"T": T}
    if m.sum() == 0:
        return {**res, "components": 0}
    P = int(point.max()) + 1
    # (comp, point) -> count, grouped to per-(comp,point) counts
    key = comp[m].astype(np.int64) * P + point[m]
    ckp, _, _ = group_reduce(key, np.ones(m.sum(), dtype=np.int64))
    # counts per (comp,point):
    _, cnt = np.unique(key, return_counts=True)
    comp_of_pair = (np.unique(key) // P).astype(np.int64)
    cu, tot_in, dom = group_reduce(comp_of_pair, cnt)
    # dominant-track coverage: fraction of in-track members in their comp's
    # dominant track (1.0 = perfectly pure components)
    dom_cov = dom.sum() / tot_in.sum()
    n_multi = len(cu)
    false_merge = int(np.sum(dom < tot_in))  # comps mixing >1 solve point

    # --- image uniqueness of multi components (real tracks: <=1 per image) ---
    ma = active & in_multi
    key_ci = comp[ma].astype(np.int64) * n + image[ma]
    uniq_ci = len(np.unique(key_ci))
    # distinct images per comp vs members per comp
    cc, members, _ = group_reduce(comp[ma], np.ones(ma.sum(), dtype=np.int64))
    # distinct images per comp:
    ci_comp = (np.unique(key_ci) // n).astype(np.int64)
    cu2, dimg, _ = group_reduce(ci_comp, np.ones(uniq_ci, dtype=np.int64))
    # align cu2 (has all multi comps with active members) with members
    img_unique_frac = float(np.mean(dimg == members))

    # --- track recovery / fragmentation ---
    mt = (point >= 0) & active
    Cn = comp.max() + 1
    key_pc = point[mt].astype(np.int64) * Cn + comp[mt]
    _, cntp = np.unique(key_pc, return_counts=True)
    pt_of_pair = (np.unique(key_pc) // Cn).astype(np.int64)
    pu, tot_mem, max_in_comp = group_reduce(pt_of_pair, cntp)
    multi_track = tot_mem >= 2  # tracks with >=2 active members
    recovery = float(np.mean((max_in_comp / tot_mem)[multi_track]))
    # fragments per track = distinct comps a track's members span
    frag_uk, frag_cnt = np.unique(pt_of_pair, return_counts=True)
    fragments = float(np.mean(frag_cnt[tot_mem >= 2]))

    in_track_total = int((point >= 0).sum())
    covered = int(((point >= 0) & active & in_multi).sum())

    return {
        **res,
        "components": n_multi,
        "dom_cov": dom_cov,
        "false_merge_comps": false_merge,
        "false_merge_frac": false_merge / max(n_multi, 1),
        "img_unique_frac": img_unique_frac,
        "recovery": recovery,
        "fragments_per_track": fragments,
        "in_track_coverage": covered / in_track_total,
    }


def build_and_score(bank, idx, dst, T, active):
    n = bank.n
    i_rep = np.repeat(np.arange(n), K)
    j = idx.ravel().astype(np.int64)
    dd = dst.ravel()
    img = bank.image_label
    keep = (
        (dd <= T)
        & (i_rep != j)
        & active[i_rep]
        & active[j]
        & (img[i_rep] != img[j])
    )
    rows = i_rep[keep]
    cols = j[keep]
    g = csr_matrix(
        (np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(n, n)
    )
    _, comp = connected_components(g, directed=False)
    sizes = np.bincount(comp, minlength=n)
    return score(bank, comp, sizes, active, T)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("sfmr")
    ap.add_argument("--outdir", default="out")
    ap.add_argument("--cache", default="out/dino_knn17.npz")
    args = ap.parse_args()

    path = sorted(glob.glob(args.sfmr))[0]
    bank = load_descriptor_bank(path)
    print(f"loaded {path}: n={bank.n} in_track={int(bank.in_track.sum())}")

    cache = Path(args.cache)
    if cache.exists():
        print(f"loading cached 17-NN from {cache}")
        z = np.load(cache)
        idx, dst = z["idx"], z["dst"]
    else:
        print(f"computing exact 17-NN over {bank.n} (one-time)...")
        idx, dst = knn_all(bank.descriptors, K)
        cache.parent.mkdir(parents=True, exist_ok=True)
        np.savez(cache, idx=idx, dst=dst)
        print(f"cached to {cache}")

    # isolated-point prefilter: d1 = nearest other, d5 = 5th other (cols 1 & 5,
    # since col 0 is self at distance 0)
    d1 = dst[:, 1]
    d5 = dst[:, 5]
    ratio = d1 / np.maximum(d5, 1e-6)
    keep_active = ~((ratio > 0.85) | (d1 > 90.0))
    all_active = np.ones(bank.n, dtype=bool)
    print(
        f"prefilter keeps {keep_active.sum()}/{bank.n} "
        f"({keep_active.mean():.1%}); drops isolated points"
    )

    thresholds = [50, 60, 70, 80, 100, 120]
    print("\n=== WITHOUT prefilter ===")
    rows_no = [build_and_score(bank, idx, dst, T, all_active) for T in thresholds]
    print("\n=== WITH isolated-point prefilter ===")
    rows_pf = [build_and_score(bank, idx, dst, T, keep_active) for T in thresholds]

    def show(rows, title):
        print(f"\n{title}")
        print(f"  {'T':>4} {'comps':>7} {'dom_cov':>8} {'falseM%':>8} "
              f"{'imgUniq':>8} {'recov':>7} {'frags':>6} {'cover':>7}")
        for r in rows:
            print(f"  {r['T']:>4} {r.get('components', 0):>7} "
                  f"{r.get('dom_cov', 0):>8.3f} {r.get('false_merge_frac', 0):>8.1%} "
                  f"{r.get('img_unique_frac', 0):>8.1%} {r.get('recovery', 0):>7.3f} "
                  f"{r.get('fragments_per_track', 0):>6.2f} "
                  f"{r.get('in_track_coverage', 0):>7.1%}")

    show(rows_no, "WITHOUT prefilter:")
    show(rows_pf, "WITH prefilter:")

    # plot recovery vs purity (dom_cov) tradeoff across T
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(8, 5))
    for rows, lbl, mk in [(rows_no, "no prefilter", "o-"), (rows_pf, "prefilter", "s-")]:
        rec = [r.get("recovery", 0) for r in rows]
        pur = [r.get("dom_cov", 0) for r in rows]
        plt.plot(rec, pur, mk, label=lbl)
        for r, x, y in zip(rows, rec, pur):
            plt.annotate(f"T={r['T']}", (x, y), fontsize=8)
    plt.xlabel("track recovery (completeness)")
    plt.ylabel("dominant-track coverage (purity)")
    plt.title("dino: radius-cluster purity vs recovery, swept over T")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    f = outdir / "exp03_purity_recovery.png"
    plt.savefig(f, dpi=110)
    plt.close()
    print(f"\nwrote {f}")


if __name__ == "__main__":
    main()
