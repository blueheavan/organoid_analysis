"""Does a voxel-spacing error survive into the segmentation study's own metric?

The downstream-bias metric in SEGMENTATION_VALIDATION_PROTOCOL.md is a *relative*
discrepancy between a production mask and a reference mask measured on the SAME
voxel grid. This script measures how much of a spacing error survives that ratio.
Two concentric spheres stand in for production (r = 12 vox) and reference
(r = 12.6 vox); volume is the voxel count times the spacing product, area is the
production marching-cubes estimator.

Result (see 2026-09-13-solution-report-review.md M5): the relative volume
discrepancy is invariant to ANY diagonal spacing error, exactly, because the
spacing product is a common factor of numerator and denominator. The relative
area discrepancy is invariant to a common scale factor to numerical precision,
and moves by 0.007 percentage points under a 20% axial-only error. The residual
dependency is therefore not the bias metric but the stratification: `rho_in` and
anisotropy are functions of the z:xy ratio.

Run: .pixi/envs/default/bin/python docs/reviews/spacing_invariance_check.py
"""
import numpy as np
from skimage.measure import marching_cubes

z, y, x = np.mgrid[-20:21, -20:21, -20:21]
PROD = (z**2 + y**2 + x**2) <= 12.0**2
REF = (z**2 + y**2 + x**2) <= 12.6**2


def area(mask, spacing):
    verts, faces, _, _ = marching_cubes(np.pad(mask.astype(float), 1), level=.5, spacing=spacing)
    tri = verts[faces]
    return .5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1).sum()


def relative_discrepancy(spacing):
    vp, vr = PROD.sum() * np.prod(spacing), REF.sum() * np.prod(spacing)
    return (vp - vr) / vr, (area(PROD, spacing) - area(REF, spacing)) / area(REF, spacing)


if __name__ == "__main__":
    for label, spacing in [("baseline", (1., 1., 1.)), ("common +3%", (1.03, 1.03, 1.03)),
                           ("common x2", (2., 2., 2.)), ("z-only +3%", (1.03, 1., 1.)),
                           ("z-only +20%", (1.2, 1., 1.))]:
        dv, da = relative_discrepancy(spacing)
        print(f"{label:12s} rel_vol_disc={dv:+.9f}  rel_area_disc={da:+.9f}")
