# Copyright The SfM Tool Authors
# SPDX-License-Identifier: Apache-2.0

"""`seed_camera.save_sfmr`'s two files, written into a small workspace.

Needs the native extension and a workspace with `.sift` files, which this
module makes for three of the seoul_bull images.
"""

import shutil
from pathlib import Path

import numpy as np
import pytest

from sfmtool import init_workspace, xxh128_of_file
from sfmtool.fileio import SiftFile, read_sfmr
from sfmtool.sift.extract import image_files_to_sift_files
from sfmtool.sift.file import get_sift_path_for_image

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


@pytest.mark.parametrize("release_grade", [True, False])
def test_seed_file_has_inline_keypoints_and_posed_images_only(
    workspace, tmp_path, monkeypatch, release_grade
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
        # What the patch-bearing file is built from: each observation's 2x2
        # warp, which member is its cluster's reference (image 0 here), and the
        # refinement radius in keypoint frames.
        "obs_warp": np.broadcast_to(np.eye(2), (len(obs_c), 2, 2)),
        "obs_ref": obs_i == 0,
        "refine_radius": 3.0,
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
        release_grade=release_grade,
        operation="seed_writer_test",
    )

    written = read_sfmr(out)
    assert written["image_names"] == [NAMES[0], NAMES[2]]
    assert written["metadata"]["operation"] == "seed_writer_test"
    # Observations are grouped by point, already the order they were given in.
    np.testing.assert_array_equal(written["keypoints_xy"], obs_uv.astype(np.float32))
    np.testing.assert_array_equal(written["image_indexes"], np.tile([0, 1], n_points))
    # Each image's thumbnail, and in a `sift_files` file its `.sift` hashes,
    # come from that image's own `.sift` file.
    for k, name in enumerate(written["image_names"]):
        with SiftFile(get_sift_path_for_image(workspace / name)) as sift:
            np.testing.assert_array_equal(
                written["thumbnails_y_x_rgb"][k], sift.read_thumbnail()
            )
            if release_grade:
                assert written["feature_tool_hashes"][k] == bytes.fromhex(
                    sift.content_hash["feature_tool_xxh128"]
                )
                assert written["sift_content_hashes"][k] == bytes.fromhex(
                    sift.content_hash["content_xxh128"]
                )
    if release_grade:
        assert written["metadata"]["feature_source"] == "sift_files"
        assert written["patch_u_halfvec_xyz"] is None
        return

    # The patch-bearing file: embedded patches, each image's file hash, and a
    # square patch frame on every point, perpendicular to its normal.
    assert written["metadata"]["feature_source"] == "embedded_patches"
    assert written["image_file_hashes"] == [
        bytes.fromhex(xxh128_of_file(workspace / name))
        for name in written["image_names"]
    ]
    u = written["patch_u_halfvec_xyz"].astype(np.float64)
    v = written["patch_v_halfvec_xyz"].astype(np.float64)
    n = written["normals_xyz"].astype(np.float64)
    assert u.shape == v.shape == n.shape == (n_points, 3)
    assert np.all(np.linalg.norm(u, axis=1) > 0)
    np.testing.assert_allclose(np.linalg.norm(u, axis=1), np.linalg.norm(v, axis=1))
    np.testing.assert_allclose(np.einsum("ij,ij->i", u, n), 0, atol=1e-5)
    np.testing.assert_allclose(np.einsum("ij,ij->i", v, n), 0, atol=1e-5)
