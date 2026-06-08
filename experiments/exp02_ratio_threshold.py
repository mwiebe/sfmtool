# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Experiment 02 — neighbor-distance structure and a clustering threshold.

For a sample of dino descriptors, find the 6 nearest (exact, over the full
corpus, excluding self) and study:

  * the d1/d5 ratio (nearest vs 5th-nearest) — how quickly distance grows past
    the closest neighbours, a proxy for "am I inside a tight cluster?",
  * the raw nearest distance d1,

split by the solve's ground truth (in-track vs background, and whether the
nearest neighbour is a co-track member). The question: does either signal give a
distance threshold that isolates real-track descriptors, usable as the radius for
a distance-bounded (max-16) neighbour lookup to seed clustering?

Usage:
    pixi run -e experiments python experiments/exp02_ratio_threshold.py \
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

from sfm_descriptors import load_descriptor_bank

RNG = np.random.default_rng(42)


def knn_with_dist(
    corpus: np.ndarray, q_idx: np.ndarray, k: int, batch: int = 512
) -> tuple[np.ndarray, np.ndarray]:
    """Exact k nearest (true L2) for the rows `q_idx`, excluding self.

    Returns (indices (Q,k), distances (Q,k)) sorted nearest-first.
    """
    X = corpus.astype(np.float32)
    sq = np.einsum("ij,ij->i", X, X)
    fetch = k + 1
    out_idx = np.empty((len(q_idx), k), dtype=np.int64)
    out_dst = np.empty((len(q_idx), k), dtype=np.float32)
    for start in range(0, len(q_idx), batch):
        qsel = q_idx[start : start + batch]
        qb = X[qsel]
        # true squared distances (include the query norm for real magnitudes)
        d2 = sq[None, :] - 2.0 * (qb @ X.T) + np.einsum("ij,ij->i", qb, qb)[:, None]
        part = np.argpartition(d2, fetch - 1, axis=1)[:, :fetch]
        rows = np.arange(len(qb))[:, None]
        order = np.argsort(d2[rows, part], axis=1)
        nn = part[rows, order]  # (b, fetch) ascending
        # drop self (the query's own global index) per row, keep k
        for r in range(len(qb)):
            row = nn[r]
            row = row[row != qsel[r]][:k]
            out_idx[start + r] = row
            out_dst[start + r] = np.sqrt(np.maximum(d2[r, row], 0.0))
    return out_idx, out_dst


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("sfmr")
    ap.add_argument("--sample", type=int, default=40000)
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()

    path = sorted(glob.glob(args.sfmr))[0]
    bank = load_descriptor_bank(path)
    print(f"loaded {path}: n={bank.n} in_track={int(bank.in_track.sum())}")

    n_q = min(args.sample, bank.n)
    q_idx = np.sort(RNG.choice(bank.n, size=n_q, replace=False))
    print(f"computing exact 6-NN for {n_q} sampled queries over {bank.n} corpus...")
    idx, dst = knn_with_dist(bank.descriptors, q_idx, k=6)

    d1 = dst[:, 0]
    d5 = dst[:, 4]
    ratio = d1 / np.maximum(d5, 1e-6)

    pid_q = bank.point_label[q_idx]
    in_track = pid_q >= 0
    pid_nn1 = bank.point_label[idx[:, 0]]
    cotrack1 = in_track & (pid_nn1 == pid_q)

    def stats(name, v, m):
        s = v[m]
        print(f"  {name:<26} n={s.size:>6}  med={np.median(s):7.3f}  "
              f"p10={np.percentile(s, 10):7.3f}  p90={np.percentile(s, 90):7.3f}")

    print("\nnearest distance d1:")
    stats("in-track", d1, in_track)
    stats("background", d1, ~in_track)
    print("d1/d5 ratio:")
    stats("in-track", ratio, in_track)
    stats("background", ratio, ~in_track)
    print(f"nearest neighbour is co-track: {cotrack1.sum()}/{in_track.sum()} "
          f"({cotrack1.sum() / max(in_track.sum(), 1):.1%} of in-track)")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Fig 1: d1/d5 ratio histogram, in-track vs background.
    plt.figure(figsize=(8, 5))
    bins = np.linspace(0, 1, 80)
    plt.hist(ratio[in_track], bins=bins, density=True, alpha=0.6, label="in-track")
    plt.hist(ratio[~in_track], bins=bins, density=True, alpha=0.6, label="background")
    plt.xlabel("d1 / d5  (nearest / 5th-nearest)")
    plt.ylabel("density")
    plt.title(f"dino: nearest/5th distance ratio (n_q={n_q})")
    plt.legend()
    plt.tight_layout()
    f1 = outdir / "exp02_ratio_hist.png"
    plt.savefig(f1, dpi=110)
    plt.close()

    # Fig 2: nearest distance d1 histogram, in-track vs background.
    plt.figure(figsize=(8, 5))
    bins = np.linspace(0, np.percentile(d1, 99), 80)
    plt.hist(d1[in_track], bins=bins, density=True, alpha=0.6, label="in-track")
    plt.hist(d1[~in_track], bins=bins, density=True, alpha=0.6, label="background")
    plt.xlabel("nearest-neighbour distance d1 (L2)")
    plt.ylabel("density")
    plt.title(f"dino: nearest-neighbour distance (n_q={n_q})")
    plt.legend()
    plt.tight_layout()
    f2 = outdir / "exp02_d1_hist.png"
    plt.savefig(f2, dpi=110)
    plt.close()

    # Fig 3: 2D density of (d1, ratio), in-track vs background side by side.
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)
    xmax = np.percentile(d1, 99)
    for ax, m, title in [
        (axes[0], in_track, "in-track"),
        (axes[1], ~in_track, "background"),
    ]:
        ax.hexbin(d1[m], ratio[m], gridsize=60, extent=(0, xmax, 0, 1),
                  bins="log", cmap="viridis")
        ax.set_xlabel("d1 (L2)")
        ax.set_title(title)
    axes[0].set_ylabel("d1 / d5")
    fig.suptitle("dino: (nearest distance, d1/d5 ratio)")
    fig.tight_layout()
    f3 = outdir / "exp02_d1_ratio_hexbin.png"
    plt.savefig(f3, dpi=110)
    plt.close()

    # Fig 4: mean per-rank distance d1..d6, in-track vs background.
    plt.figure(figsize=(8, 5))
    ranks = np.arange(1, 7)
    plt.plot(ranks, dst[in_track].mean(0), "o-", label="in-track")
    plt.plot(ranks, dst[~in_track].mean(0), "s-", label="background")
    plt.xlabel("neighbour rank")
    plt.ylabel("mean distance (L2)")
    plt.title("dino: mean distance by neighbour rank")
    plt.legend()
    plt.tight_layout()
    f4 = outdir / "exp02_rank_curve.png"
    plt.savefig(f4, dpi=110)
    plt.close()

    print(f"\nwrote: {f1}\n       {f2}\n       {f3}\n       {f4}")


if __name__ == "__main__":
    main()
