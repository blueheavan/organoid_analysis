from __future__ import annotations

import tempfile
from pathlib import Path
import unittest

import numpy as np
import tifffile

from organoid_analysis.segmentation.parameter_estimation import (
    _downsample_2d,
    auto_anisotropy,
    estimate_diameter_from_stack,
)


class _CircleModel:
    """Mock Cellpose model that returns filled-circle masks.

    The circle radius (in full-resolution pixels) is scaled to the actual input
    size so downsampling does not change the derived full-res diameter.
    """

    radius_px = 15.0

    def __init__(self, radius_px: float = 15.0):
        self.radius_px = radius_px

    def eval(self, work, **kwargs):
        work = work[..., 0]
        h, w = work.shape
        scale = h / 100.0
        r = self.radius_px * scale
        yy, xx = np.mgrid[0:h, 0:w]
        c1 = ((yy - 30 * scale) ** 2 + (xx - 30 * scale) ** 2) < r**2
        c2 = ((yy - 70 * scale) ** 2 + (xx - 30 * scale) ** 2) < r**2
        c3 = ((yy - 50 * scale) ** 2 + (xx - 72 * scale) ** 2) < r**2
        masks = np.zeros(work.shape, dtype=int)
        masks[c1] = 1
        masks[c2] = 2
        masks[c3] = 3
        return [masks]


class DownsampleTests(unittest.TestCase):
    def test_full_resolution_passthrough(self) -> None:
        img = np.ones((100, 120))
        out, factor = _downsample_2d(img, 1.0)
        self.assertEqual(factor, 1.0)
        self.assertEqual(out.shape, (100, 120))

    def test_half_resolution_scales(self) -> None:
        img = np.ones((100, 120))
        out, factor = _downsample_2d(img, 0.5)
        self.assertAlmostEqual(factor, 0.5)
        self.assertLessEqual(out.shape[0], 51)
        self.assertLessEqual(out.shape[1], 61)


class AnisotropyTests(unittest.TestCase):
    def test_anisotropy_from_ome_units(self) -> None:
        # Write a real OME-TIFF: tifffile turns a dict of physical sizes into
        # embedded OME-XML with proper units.
        path = Path(tempfile.mkdtemp()) / "stack.ome.tiff"
        img = np.zeros((4, 8, 8), dtype=np.uint16)
        tifffile.imwrite(
            path,
            img,
            metadata={
                "PhysicalSizeX": 0.4,
                "PhysicalSizeXUnit": "µm",
                "PhysicalSizeY": 0.4,
                "PhysicalSizeYUnit": "µm",
                "PhysicalSizeZ": 5.0,
                "PhysicalSizeZUnit": "µm",
            },
            ome=True,
        )
        anisotropy = auto_anisotropy(path)
        # z=5um, xy=0.4um -> anisotropy = 5/0.4 = 12.5
        self.assertIsNotNone(anisotropy)
        self.assertAlmostEqual(anisotropy, 12.5, places=3)

    def test_anisotropy_unknown_when_no_metadata(self) -> None:
        path = Path(tempfile.mkdtemp()) / "plain.tiff"
        img = np.zeros((4, 8, 8), dtype=np.uint16)
        tifffile.imwrite(path, img)  # no OME/ImageJ/resolution metadata
        self.assertIsNone(auto_anisotropy(path))


class DiameterEstimateTests(unittest.TestCase):
    def test_estimate_recovers_filled_circle_diameter(self) -> None:
        stack = np.zeros((10, 100, 100), dtype=np.float32)
        est = estimate_diameter_from_stack(
            _CircleModel(radius_px=15.0), stack, downsample=0.5, max_slices=10
        )
        # Equivalent disk diameter of a radius-15 circle ~= 30 px.
        self.assertIsNotNone(est)
        self.assertTrue(25 <= est <= 35, est)

    def test_estimate_consistent_across_resolutions(self) -> None:
        stack = np.zeros((10, 100, 100), dtype=np.float32)
        full = estimate_diameter_from_stack(
            _CircleModel(radius_px=15.0), stack, downsample=1.0, max_slices=10
        )
        half = estimate_diameter_from_stack(
            _CircleModel(radius_px=15.0), stack, downsample=0.5, max_slices=10
        )
        self.assertIsNotNone(full)
        self.assertIsNotNone(half)
        self.assertAlmostEqual(full, half, delta=3.0)

    def test_rejects_non_3d(self) -> None:
        with self.assertRaises(ValueError):
            estimate_diameter_from_stack(_CircleModel(), np.zeros((100, 100)))


if __name__ == "__main__":
    unittest.main()
