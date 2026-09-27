# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""Score runs on the Kerry Park cases of `specs/drafts/surface-co-solve.md`.

Each case is a set of points of the Kerry Park candidate ground truth
`tk113`. For each run, case and pass the table gives the good tracks over the
queries, the normal credit per query (`N` in `goal_score.py`), the score
`S = (G + 2N) / 3` and the median normal error of the good tracks. The ground
sheet around point 112 is split into its flat points and its sloped ones
(normal 5 degrees or more from vertical), read from the run's ground truth.

    pixi run -e test python scripts/track_at_pixel/surface/cases.py RUN [RUN ...]

A quick run over the case points alone, for iterating on a candidate:

    pixi run -e test python scripts/track_at_pixel/harness.py --dataset <tk113.sfmr> --point-ids $(pixi run -e test python scripts/track_at_pixel/surface/cases.py --ids)
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from goal_score import score_rows  # noqa: E402

CASES = {
    "rock 322": [322],
    "rock 321": [321],
    "grass": [348, 349, 350, 351, 352],
    "hedge": [12, *range(133, 141)],
    "window": [15, 121, 122, 124],
    "chimney 218": [218],
    "house side 20": [20],
    "sign": list(range(326, 334)),
}
# The ground and the hill around point 112: a sheet of ground-truth patches.
GROUND = [
    0,
    34,
    55,
    69,
    75,
    76,
    77,
    78,
    79,
    80,
    83,
    84,
    85,
    86,
    87,
    88,
    89,
    90,
    91,
    92,
    93,
    94,
    95,
    96,
    97,
    98,
    99,
    100,
    103,
    104,
    105,
    106,
    107,
    108,
    109,
    110,
    111,
    112,
    113,
    114,
    115,
    116,
    117,
    118,
    119,
    160,
    161,
    162,
    163,
    164,
    165,
    166,
    167,
    171,
    172,
    173,
    175,
    177,
    180,
    181,
    182,
    183,
    184,
    185,
    186,
    187,
    188,
    189,
    190,
    191,
    192,
    193,
    194,
    196,
    197,
    198,
    199,
    200,
    201,
    202,
    203,
    204,
    210,
    211,
    212,
    213,
    214,
    215,
    216,
    217,
    219,
    222,
    223,
    224,
    225,
    226,
    227,
    228,
    229,
    230,
    231,
    232,
    233,
    234,
    235,
    236,
    237,
    238,
    239,
    240,
    241,
    242,
    243,
    244,
    245,
    246,
    263,
    266,
    267,
    268,
    271,
    281,
    282,
    283,
    284,
    285,
    286,
    287,
    288,
    289,
    290,
    291,
    292,
    293,
    294,
    309,
    310,
    311,
    312,
    313,
    314,
    322,
    323,
    324,
    325,
    359,
    360,
    361,
]


def ground_split(dataset: str):
    """The ground sheet's flat points and its sloped ones."""
    from sfmtool._sfmtool.reconstruction import EditedReconstruction, SfmrReconstruction

    if not Path(dataset).is_file():
        # A harness run records the dataset's name; its cache holds the file.
        from dataset import default_cache_dir

        dataset = default_cache_dir() / dataset / f"{dataset}.sfmr"
    edited = EditedReconstruction(SfmrReconstruction.load(str(dataset)))
    flat, sloped = [], []
    for p in GROUND:
        rec = edited.point(p)
        n = np.cross(rec["patch_u_halfvec"], rec["patch_v_halfvec"])
        tilt = np.degrees(np.arccos(min(1.0, abs(n[2]) / np.linalg.norm(n))))
        (sloped if tilt >= 5 else flat).append(p)
    return flat, sloped


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("runs", nargs="*", type=Path)
    ap.add_argument("--ids", action="store_true", help="print the case points and exit")
    args = ap.parse_args()
    if args.ids:
        print(",".join(str(p) for p in sorted({*GROUND, *sum(CASES.values(), [])})))
        return
    cases = dict(CASES)
    for run in args.runs:
        dataset = json.load(open(run / "config.json"))["dataset"]
        flat, sloped = ground_split(dataset)
        cases["ground flat"], cases["ground sloped"] = flat, sloped
        rows = [json.loads(line) for line in open(run / "rows.jsonl")]
        print(f"\n{run.name}")
        print(
            f"{'case':14s} {'pts':>4s} | {'full q':>6s} {'good':>5s} {'N':>5s} {'S':>5s} "
            f"{'nrm°':>5s} | {'empty q':>7s} {'good':>5s} {'N':>5s} {'S':>5s} {'nrm°':>5s}"
        )
        for name, points in cases.items():
            want = set(points)
            cells = []
            for pas in ("full", "empty"):
                sel = [r for r in rows if r["pass"] == pas and r["point"] in want]
                if not sel:
                    cells.append(f"{'--':>6s} {'':5s} {'':5s} {'':5s} {'':5s}")
                    continue
                m = score_rows(sel)
                cells.append(
                    f"{m['queries']:6d} {m['good']:5d} {m['N']:5.2f} {m['S']:5.2f} "
                    f"{m['normal_med']:5.1f}"
                )
            print(f"{name:14s} {len(points):4d} | {cells[0]} | {cells[1]:>}")


if __name__ == "__main__":
    main()
