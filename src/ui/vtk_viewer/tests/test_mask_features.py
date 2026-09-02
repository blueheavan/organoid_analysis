"""Unit tests for mask feature extraction (``mask_features``)."""

from __future__ import annotations

import unittest

import numpy as np

from ui.vtk_viewer.mask_features import (
    count_mask_objects,
    extract_mask_features,
    mask_to_binary,
    summarize_features,
    summarize_mask,
    validate_mask,
)


def _two_spheres(shape=(24, 32, 32)) -> np.ndarray:
    """Two solid unit spheres of different radius, labelled 1 and 2."""
    mask = np.zeros(shape, dtype=np.int64)
    z, y, x = shape
    zz, yy, xx = np.mgrid[0:z, 0:y, 0:x]
    mask[(zz - z / 2) ** 2 + (yy - y / 2) ** 2 + (xx - x / 2) ** 2 < 5.0**2] = 1
    mask[(zz - z / 2) ** 2 + (yy - y / 4) ** 2 + (xx - x * 0.7) ** 2 < 3.0**2] = 2
    return mask


class ValidateMaskTests(unittest.TestCase):
    def test_rejects_non_3d(self) -> None:
        with self.assertRaises(ValueError):
            validate_mask(np.zeros((8, 8), dtype=np.uint8))

    def test_rejects_tiny(self) -> None:
        with self.assertRaises(ValueError):
            validate_mask(np.zeros((1, 8, 8), dtype=np.uint8))

    def test_rejects_fractional_float_labels(self) -> None:
        mask = np.zeros((4, 5, 6), dtype=np.float32)
        mask[0, 0, 0] = 3.7
        with self.assertRaisesRegex(ValueError, "integer"):
            validate_mask(mask)

    def test_rejects_negative_labels(self) -> None:
        mask = np.zeros((4, 5, 6), dtype=np.int32)
        mask[0, 0, 0] = -1
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            validate_mask(mask)

    def test_preserves_valid_integer_dtype_without_copy(self) -> None:
        mask = np.zeros((4, 5, 6), dtype=np.uint16)
        self.assertIs(validate_mask(mask), mask)

    def test_sparse_ids_are_counted_by_occurrence(self) -> None:
        mask = np.zeros((4, 5, 6), dtype=np.uint32)
        mask[0, 0, 0] = 7
        mask[1, 1, 1] = 1_000_000
        self.assertEqual(count_mask_objects(mask), 2)


class FeatureTests(unittest.TestCase):
    def test_counts_two_objects_with_correct_labels(self) -> None:
        mask = _two_spheres()
        features = extract_mask_features(mask, spacing_um=(1.0, 1.0, 1.0))
        self.assertEqual(set(features.index.tolist()), {1, 2})
        self.assertEqual(features.shape[0], 2)

    def test_bigger_object_has_greater_volume(self) -> None:
        mask = _two_spheres()
        features = extract_mask_features(mask, spacing_um=(1.0, 1.0, 1.0))
        self.assertGreater(
            features.loc[1, "volume_um3"], features.loc[2, "volume_um3"]
        )

    def test_sphericity_close_to_one_for_spheres(self) -> None:
        mask = _two_spheres()
        features = extract_mask_features(mask, spacing_um=(1.0, 1.0, 1.0))
        for label in (1, 2):
            sphericity = features.loc[label, "sphericity"]
            self.assertGreater(sphericity, 0.9)
            self.assertLessEqual(sphericity, 1.0)

    def test_spacing_scales_volume(self) -> None:
        mask = _two_spheres()
        unit = extract_mask_features(mask, spacing_um=(1.0, 1.0, 1.0))
        scaled = extract_mask_features(mask, spacing_um=(2.0, 1.0, 1.0))
        self.assertAlmostEqual(
            scaled.loc[1, "volume_um3"], unit.loc[1, "volume_um3"] * 2.0, places=2
        )

    def test_empty_mask_returns_empty_frame(self) -> None:
        features = extract_mask_features(np.zeros((8, 8, 8), dtype=np.uint8))
        self.assertTrue(features.empty)

    def test_invalid_spacing_rejected(self) -> None:
        with self.assertRaises(ValueError):
            extract_mask_features(_two_spheres(), spacing_um=(1.0, 1.0))


class SummaryTests(unittest.TestCase):
    def test_summary_counts_and_mean_volume(self) -> None:
        mask = _two_spheres()
        summary = summarize_mask(mask, spacing_um=(1.0, 1.0, 1.0))
        self.assertEqual(summary.n_objects, 2)
        self.assertGreater(summary.total_volume_um3, 0)
        self.assertAlmostEqual(
            summary.mean_volume_um3,
            summary.total_volume_um3 / 2,
            places=3,
        )
        self.assertLessEqual(summary.min_volume_um3, summary.max_volume_um3)

    def test_summary_empty(self) -> None:
        summary = summarize_mask(np.zeros((8, 8, 8), dtype=np.uint8))
        self.assertEqual(summary.n_objects, 0)
        self.assertEqual(summary.total_volume_um3, 0.0)

    def test_existing_features_can_be_summarized_without_reextracting(self) -> None:
        features = extract_mask_features(_two_spheres())
        summary = summarize_features(features)
        self.assertEqual(summary.n_objects, 2)
        self.assertAlmostEqual(summary.total_volume_um3, features.volume_um3.sum())


class BinaryTests(unittest.TestCase):
    def test_mask_to_binary(self) -> None:
        mask = _two_spheres()
        binary = mask_to_binary(mask)
        self.assertEqual(binary.dtype, bool)
        self.assertEqual(int(binary.sum()), int((mask > 0).sum()))


if __name__ == "__main__":
    unittest.main(verbosity=2)
