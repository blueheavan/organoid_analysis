from __future__ import annotations

import json
from pathlib import Path
import tempfile
import threading
import time
import unittest
from unittest.mock import patch

import numpy as np
import tifffile

from segmentation.cellpose import (
    SegmentationConfig,
    _wrap_run_3d,
    normalize_preview,
    list_saved_results,
    read_stack,
    read_stack_multichannel,
    restore_saved_result,
    save_result,
    segment_stacks,
    validate_stacks,
)


class SegmentationTests(unittest.TestCase):
    def test_read_stack_accepts_3d_tiff(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "stack.tif"
            expected = np.ones((2, 4, 4), dtype=np.uint16)
            tifffile.imwrite(path, expected)
            np.testing.assert_array_equal(read_stack(path), expected)

    def test_multichannel_reader_does_not_reinterpret_cyx_as_zyx(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "channels.ome.tif"
            tifffile.imwrite(path, np.ones((2, 4, 4), np.uint16), ome=True,
                             photometric="minisblack", metadata={"axes": "CYX"})
            with self.assertRaisesRegex(ValueError, "Z/Y/X"):
                read_stack_multichannel(path)

    def test_multichannel_reader_does_not_fallback_for_metadata_tiff(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "stack.ome.tif"
            tifffile.imwrite(path, np.ones((3, 4, 4), np.uint16), ome=True,
                             photometric="minisblack", metadata={"axes": "ZYX"})
            with (patch("segmentation.cellpose.load_zstack", side_effect=ValueError("rejected")),
                  patch("segmentation.cellpose.read_axes", return_value="QYX"),
                  self.assertRaisesRegex(ValueError, "rejected")):
                read_stack_multichannel(path)

    def test_validate_stacks_rejects_shape_mismatch(self) -> None:
        with self.assertRaisesRegex(ValueError, "identical shapes"):
            validate_stacks(np.zeros((2, 4, 4)), np.zeros((3, 4, 4)))

    def test_normalize_preview_handles_constant_data(self) -> None:
        normalized = normalize_preview(np.ones((4, 4), dtype=np.uint16))
        self.assertTrue(np.array_equal(normalized, np.zeros((4, 4), dtype=np.float32)))

    def test_save_result_writes_archive_and_metadata(self) -> None:
        nuclei = np.array([[[0, 1], [2, 0]]], dtype=np.uint32)
        cells = np.array([[[0, 1], [1, 0]]], dtype=np.uint32)
        with tempfile.TemporaryDirectory() as temporary_directory:
            result = save_result(
                nuclei,
                cells,
                SegmentationConfig(),
                output_directory=Path(temporary_directory),
            )
            self.assertTrue(result.archive_path.exists())
            summary = json.loads(result.summary_path.read_text(encoding="utf-8"))
            self.assertEqual(summary["nuclei_count"], 2)
            self.assertEqual(summary["cell_count"], 1)

    def test_save_result_counts_sparse_labels_not_max_id(self) -> None:
        nuclei = np.array([[[0, 7], [1_000_000, 0]]], dtype=np.uint32)
        with tempfile.TemporaryDirectory() as temporary_directory:
            result = save_result(nuclei, None, SegmentationConfig(), Path(temporary_directory))
            self.assertEqual(result.nuclei_count, 2)

    def test_saved_result_can_be_listed_and_restored_after_session_loss(self) -> None:
        nuclei = np.array([[[0, 7], [9, 0]]], dtype=np.uint32)
        cells = np.array([[[0, 2], [2, 0]]], dtype=np.uint32)
        config = SegmentationConfig(anisotropy=2.0, xy_spacing_um=0.65)
        with tempfile.TemporaryDirectory() as temporary_directory:
            output = Path(temporary_directory)
            saved = save_result(nuclei, cells, config, output)
            self.assertEqual([item.run_directory for item in list_saved_results(output)], [saved.run_directory])
            restored = restore_saved_result(saved.run_directory)
            self.assertEqual(restored.config, config)
            np.testing.assert_array_equal(restored.nuclei_masks, nuclei)
            np.testing.assert_array_equal(restored.cell_masks, cells)

    def test_concurrent_model_inference_is_serialized(self) -> None:
        class FakeModel:
            def __init__(self):
                self.active = 0
                self.maximum_active = 0
                self.guard = threading.Lock()

            def eval(self, image, **kwargs):
                with self.guard:
                    self.active += 1
                    self.maximum_active = max(self.maximum_active, self.active)
                time.sleep(0.03)
                with self.guard:
                    self.active -= 1
                return (np.zeros(image.shape[:3], dtype=np.uint32),)

        model = FakeModel()
        image = np.zeros((4, 8, 8), dtype=np.uint16)
        threads = [threading.Thread(target=segment_stacks,
                                    args=(model, image, None, SegmentationConfig()))
                   for _ in range(2)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        self.assertEqual(model.maximum_active, 1)

    def test_wrap_run_3d_maps_plane_reports_to_stage_fractions(self) -> None:
        from cellpose import models

        # Cellpose's run_3D reports 25 / 40 / 55 after each orthogonal plane.
        def recorder(net, imgs, **kwargs):
            kwargs["progress"].setValue(25)
            kwargs["progress"].setValue(40)
            kwargs["progress"].setValue(55)
            return ("yf", "styles")

        original = models.run_3D
        models.run_3D = recorder
        try:
            # Nuclei stage when a cell pass follows: 10% -> 55%.
            progress = []
            restore = _wrap_run_3d(progress.append, 0.10, 0.55)
            models.run_3D(None, None)
            restore()
            self.assertEqual([round(p, 3) for p in progress], [0.1, 0.325, 0.55])

            # Nuclei-only stage: 10% -> 95% (last 5% guards mask post-processing).
            progress = []
            restore = _wrap_run_3d(progress.append, 0.10, 0.95)
            models.run_3D(None, None)
            restore()
            self.assertEqual([round(p, 3) for p in progress], [0.1, 0.525, 0.95])

            # Cell stage: 55% -> 100%.
            progress = []
            restore = _wrap_run_3d(progress.append, 0.55, 1.0)
            models.run_3D(None, None)
            restore()
            self.assertEqual([round(p, 3) for p in progress], [0.55, 0.775, 1.0])

            # Restored back to the patched-in recorder after each pass.
            self.assertIs(models.run_3D, recorder)

            # on_progress=None is a no-op, not a crash.
            progress = []
            restore = _wrap_run_3d(None, 0.10, 0.55)
            models.run_3D(None, None)
            restore()
            self.assertEqual(progress, [])
        finally:
            models.run_3D = original


if __name__ == "__main__":
    unittest.main()
