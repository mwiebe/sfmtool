# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""`seed_camera.save_sfmr`'s release-grade file, written into a small workspace.

Needs the native extension and a workspace with `.sift` files, which this
module makes for three of the seoul_bull images.
"""

import shutil
from pathlib import Path

import numpy as np
import pytest

from sfmtool import init_workspace
from sfmtool.fileio import read_sfmr
from sfmtool.sift.extract import image_files_to_sift_files

SEOUL_BULL_DIR = (
    Path(__file__).resolve().parents[2]
    / "test-data"
    / "images"
    / "seoul_bull_sculpture"
)
NAMES = [f"seoul_bull_sculpture_{k:02d}.jpg" for k in (1, 2, 3)]
WIDTH, HEIGHT = 270, 480


@pytest.fixture(scope="module")
def workspace(tmp_path_factory):
    workspace_dir = tmp_path_factory.mktemp("seed_writer")
    config = init_workspace(workspace_dir, feature_tool="sfmtool", max_num_features=50)
    images = []
    for name in NAMES:
        shutil.copy(SEOUL_BULL_DIR / name, workspace_dir / name)
        images.append(workspace_dir / name)
    image_files_to_sift_files(
        images,
        feature_tool=config["feature_tool"],
        feature_options=config["feature_options"],
        feature_prefix_dir=config["feature_prefix_dir"],
    )
    return workspace_dir


def test_release_grade_file_has_inline_keypoints_and_posed_images_only(
    workspace, tmp_path, monkeypatch
):
    import seed_camera

    monkeypatch.setattr(seed_camera, "WS", workspace)
    monkeypatch.setattr(seed_camera, "_CAM_WH", (WIDTH, HEIGHT))

    rng = np.random.default_rng(3)
    n_points = 4
    # Images 0 and 2 see every point; image 1 sees none, so it is not posed.
    obs_c = np.repeat(np.arange(n_points), 2)
    obs_i = np.tile([0, 2], n_points)
    obs_uv = np.column_stack(
        [
            rng.uniform(10, WIDTH - 10, len(obs_c)),
            rng.uniform(10, HEIGHT - 10, len(obs_c)),
        ]
    )
    data = {
        "names": NAMES,
        "dims": [(WIDTH, HEIGHT)] * len(NAMES),
        "obs_c": obs_c,
        "obs_i": obs_i,
        "obs_f": np.arange(len(obs_c)),
        "obs_uv": obs_uv,
    }
    pts = np.column_stack(
        [
            rng.uniform(-1, 1, n_points),
            rng.uniform(-1, 1, n_points),
            -5.0 * np.ones(n_points),
        ]
    )
    out = tmp_path / "seed.sfmr"
    seed_camera.save_sfmr(
        data,
        300.0,
        np.zeros((len(NAMES), 3)),
        np.zeros((len(NAMES), 3)),
        pts,
        np.ones(len(obs_c), bool),
        np.full(len(obs_c), 0.5),
        out,
        quiet=True,
        release_grade=True,
        operation="seed_writer_test",
    )

    written = read_sfmr(out)
    assert written["image_names"] == [NAMES[0], NAMES[2]]
    assert written["metadata"]["operation"] == "seed_writer_test"
    # Observations are grouped by point, already the order they were given in.
    np.testing.assert_array_equal(written["keypoints_xy"], obs_uv.astype(np.float32))
    np.testing.assert_array_equal(written["image_indexes"], np.tile([0, 1], n_points))
