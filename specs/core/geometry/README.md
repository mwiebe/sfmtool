# Geometry Specifications

Pose estimation, epipolar geometry, and optimization. Implemented in
`crates/sfmtool-core/src/geometry/`.

| Document | Description |
|----------|-------------|
| [absolute-pose.md](absolute-pose.md) | Camera pose from 2D-3D correspondences: P3P plus RANSAC. |
| [relative-pose.md](relative-pose.md) | Relative pose between two cameras from ray correspondences. |
| [baseline-direction.md](baseline-direction.md) | The direction between two centres from ray coplanarity, with the rotations held, batched over a graph. |
| [translation-averaging.md](translation-averaging.md) | Camera centres from pairwise baseline directions and relative lengths, read off the null space of the angular form; the orientation bit from cheirality. |
| [seed-hypothesis-loop.md](seed-hypothesis-loop.md) | Developing and committing the whole set of seed candidates a capture's cluster evidence supports; the set is the product. |
| [seed-candidate-evaluation.md](seed-candidate-evaluation.md) | The gauge-free battery every committed seed candidate is measured by before release. |
| [seed-drive.md](seed-drive.md) | Judging the candidate set while it is still being produced: pull, complete, judge, and form the next complement from the survivors' claims. |
| [seed-relaxation.md](seed-relaxation.md) | Relaxing a rotation-only member into finite geometry: baselines from its refused rows, fill-in from its source clusters, a lens read on the result. |
| [rotation-locked-resection.md](rotation-locked-resection.md) | Translation-only resection with the rotation held fixed. |
| [epipolar-estimation.md](epipolar-estimation.md) | Fundamental matrix from 2D-2D correspondences: 7-point with RANSAC, and Bougnoux focal recovery. |
| [focal-vote.md](focal-vote.md) | Structure-free shared-focal estimation, where image pairs vote independently through per-pair estimators. |
| [estimate-intrinsics.md](estimate-intrinsics.md) | One typed camera out of the focal vote: the model verdict, whether it is corroborated, the focal, and the votes behind them. |
| [rotation-init.md](rotation-init.md) | Far-field, parallax-free correspondences fixing rotations before any translation is known. |
| [reconstruction-growth.md](reconstruction-growth.md) | Registering the un-posed images of a cluster-track set against a seeded reconstruction, in batches. |
| [bundle-adjustment.md](bundle-adjustment.md) | Staged bundle adjustment over a shared camera. |
| [reprojection-residuals.md](reprojection-residuals.md) | Batched reprojection residuals and inlier fractions. |
| [affine-factorization.md](affine-factorization.md) | Joint weak-perspective factorization over sparse, partly-junk cluster observations. |
| [pose-verification.md](pose-verification.md) | Displacement-neighbourhood check that flags — and repairs — poses the structure disagrees with. |
