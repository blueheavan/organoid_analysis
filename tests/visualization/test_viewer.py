"""Smoke tests for the vtk.js volume viewer data-packing / HTML pipeline.

These verify the Python side (payload building, dimension ordering, quota and
validation) without needing a browser or GPU. The browser-level render is the
responsibility of the Playwright probe in ``scripts/render_probe.py``.
"""

from __future__ import annotations

import base64
import json
import unittest
import zlib

import numpy as np

from organoid_analysis.visualization.volume_viewer import (
    MAX_CHANNELS,
    build_surface_payload,
    build_viewer_payload,
    render_viewer_html,
    spec_to_dict,
    validate_spacing,
)
from organoid_analysis.visualization.volume_viewer.rendering_presets import PRESETS
from organoid_analysis.visualization.volume_viewer.viewer_payload import _voxel_major_zyxc


def _sphere_stack(shape=(16, 20, 24), seed=0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    vol = rng.random(shape).astype(np.float32)
    # put a bright blob in the middle so volume/percentile logic is exercised
    z, y, x = shape
    zz, yy, xx = np.mgrid[0:z, 0:y, 0:x]
    vol[(zz - z / 2) ** 2 + (yy - y / 2) ** 2 + (xx - x / 2) ** 2 < 5.0] = 1.0
    return vol


class SpacingTests(unittest.TestCase):
    def test_valid_spacing_step_x_y_z(self) -> None:
        self.assertEqual(validate_spacing((0.5, 1.0, 2.0)), (0.5, 1.0, 2.0))

    def test_rejects_wrong_length(self) -> None:
        with self.assertRaises(ValueError):
            validate_spacing((1.0, 1.0))

    def test_rejects_nonpositive(self) -> None:
        with self.assertRaises(ValueError):
            validate_spacing((0.0, 1.0, 1.0))
        with self.assertRaises(ValueError):
            validate_spacing((-1.0, 1.0, 1.0))


class ChannelQuotaTests(unittest.TestCase):
    def test_max_channels(self) -> None:
        self.assertEqual(MAX_CHANNELS, 4)

    def test_too_many_channels_rejected(self) -> None:
        vols = [_sphere_stack() for _ in range(MAX_CHANNELS + 1)]
        with self.assertRaises(ValueError):
            build_viewer_payload(vols, (1.0, 1.0, 1.0))

    def test_float64_rejected(self) -> None:
        vols = [_sphere_stack().astype(np.float64)]
        with self.assertRaises(ValueError):
            build_viewer_payload(vols, (1.0, 1.0, 1.0))

    def test_default_fluorescence_lut_order(self) -> None:
        spec = build_viewer_payload(
            [_sphere_stack((4, 5, 6), seed=i) for i in range(4)],
            (1.0, 1.0, 1.0),
        )
        self.assertEqual(
            [channel["lut_key"] for channel in spec.channels],
            ["cyan", "magenta", "yellow", "green"],
        )

    def test_fluorescence_luts_are_linear_pure_colors(self) -> None:
        expected = {
            "cyan": (1.0, 0.0, 1.0, 1.0),
            "magenta": (1.0, 1.0, 0.0, 1.0),
            "yellow": (1.0, 1.0, 1.0, 0.0),
            "green": (1.0, 0.0, 1.0, 0.0),
        }
        for name, endpoint in expected.items():
            self.assertEqual(len(PRESETS[name]), 2)
            self.assertEqual(PRESETS[name][0], (0.0, 0.0, 0.0, 0.0))
            self.assertEqual(PRESETS[name][1], endpoint)


class DimensionOrderingTests(unittest.TestCase):
    def test_dims_reported_x_y_z_and_interleaved(self) -> None:
        """A (Z, Y, X) = (16, 20, 24) input must become dims (X, Y, Z)."""
        vols = [_sphere_stack((16, 20, 24))]
        spec = build_viewer_payload(vols, (1.0, 1.0, 1.0))
        self.assertEqual(tuple(spec.dims), (24, 20, 16))
        self.assertEqual(spec.n_channels, 1)


class PayloadTests(unittest.TestCase):
    def test_payload_preserves_raw_values_and_full_domain_metadata(self) -> None:
        arr = np.arange(3 * 4 * 5, dtype=np.float32).reshape(3, 4, 5)
        arr[0, 0, 0] = 1000.0
        spec = build_viewer_payload([arr], (1.0, 1.0, 1.0))
        decoded = np.frombuffer(base64.b64decode(spec.data_b64), dtype=np.float32)
        np.testing.assert_array_equal(decoded.reshape(3, 4, 5), arr)
        channel = spec.channels[0]
        self.assertEqual(channel["actual_min"], float(arr.min()))
        self.assertEqual(channel["actual_max"], float(arr.max()))
        self.assertGreaterEqual(channel["suggested_lo"], channel["actual_min"])
        self.assertLessEqual(channel["suggested_hi"], channel["actual_max"])

    def test_single_channel_wire_roundtrip_preserves_zyx_order(self) -> None:
        arr = np.zeros((32, 64, 64), dtype=np.float32)
        arr[5, 10, 20] = 1000.0
        wire = _voxel_major_zyxc(arr[None, ...])
        flat = np.frombuffer(wire.tobytes(), dtype=np.float32)
        expected_index = 5 * 64 * 64 + 10 * 64 + 20
        self.assertEqual(int(np.flatnonzero(flat == 1000.0)[0]), expected_index)
        np.testing.assert_array_equal(flat.reshape(arr.shape), arr)

    def test_multichannel_wire_is_voxel_major(self) -> None:
        z, y, x = np.indices((3, 4, 5))
        c0 = (100 * z + 10 * y + x).astype(np.float32)
        c1 = (1000 + c0).astype(np.float32)
        wire = _voxel_major_zyxc(np.stack([c0, c1]))
        flat = np.frombuffer(wire.tobytes(), dtype=np.float32)
        decoded = flat.reshape(3, 4, 5, 2)
        np.testing.assert_array_equal(decoded[..., 0], c0)
        np.testing.assert_array_equal(decoded[..., 1], c1)

    def test_base64_decodes_to_expected_float32_size(self) -> None:
        vols = [_sphere_stack((16, 20, 24))]
        spec = build_viewer_payload(vols, (1.0, 1.0, 2.0))
        raw = base64.b64decode(spec.data_b64)
        self.assertGreater(len(raw), 0)
        self.assertEqual(len(raw) % 4, 0)
        self.assertEqual(spec.compressed, False)

    def test_compressed_payload_is_zlib(self) -> None:
        vols = [_sphere_stack((16, 20, 24))]
        spec = build_viewer_payload(vols, (1.0, 1.0, 2.0), compress=True)
        self.assertTrue(spec.compressed)
        raw = zlib.decompress(base64.b64decode(spec.data_b64))
        self.assertGreater(len(raw), 0)

    def test_multichannel_interleave_has_last_dim_c(self) -> None:
        vols = [_sphere_stack((8, 10, 12)), _sphere_stack((8, 10, 12))]
        spec = build_viewer_payload(vols, (1.0, 1.0, 2.0))
        self.assertEqual(spec.n_channels, 2)
        self.assertEqual(len(spec.channels), 2)

    def test_spec_survives_json_roundtrip(self) -> None:
        vols = [_sphere_stack((8, 10, 12))]
        spec = build_viewer_payload(vols, (1.0, 1.0, 2.0))
        d = spec_to_dict(spec)
        loaded = json.loads(json.dumps(d))
        self.assertEqual(loaded["dims"], list(spec.dims))
        self.assertEqual(loaded["mode"], "volume")


class SurfaceTests(unittest.TestCase):
    def test_surface_payload_uses_legacy_cells(self) -> None:
        mask = (np.random.default_rng(3).random((24, 24, 24)) > 0.7).astype(np.uint8)
        payload = build_surface_payload(mask, spacing_zyx=(2.0, 1.0, 1.0))
        self.assertGreater(payload["n_faces"], 0)
        self.assertGreater(payload["n_points"], 0)
        faces = np.frombuffer(
            base64.b64decode(payload["faces_b64"]), dtype=np.uint32
        )
        # flat legacy arrays are [3, i, j, k, 3, ...] -> length is 4 * n_faces
        self.assertEqual(faces.size, 4 * payload["n_faces"])
        self.assertTrue((faces[0::4] == 3).all())

    def test_surface_mode_requires_mask(self) -> None:
        with self.assertRaises(ValueError):
            build_viewer_payload(
                [_sphere_stack()], (1.0, 1.0, 2.0), render_mode="surface"
            )

    def test_all_foreground_surface_is_closed_by_padding(self) -> None:
        payload = build_surface_payload(np.ones((4, 5, 6), np.uint8), (2.0, 1.0, 1.0))
        self.assertGreater(payload["n_faces"], 0)


class MaskOverlayTests(unittest.TestCase):
    def _labeled(self, shape=(16, 20, 24)) -> np.ndarray:
        labels = np.zeros(shape, dtype=np.int64)
        labels[2:6, 3:9, 4:12] = 1
        labels[8:12, 12:18, 14:20] = 2
        return labels

    def test_payload_roundtrips_labels_and_shape(self) -> None:
        labels = self._labeled()
        spec = build_viewer_payload(
            [_sphere_stack(labels.shape)],
            (1.0, 1.0, 2.0),
            mask_overlay=labels,
            mask_overlay_alpha=0.4,
        )
        self.assertIsNotNone(spec.mask)
        mask = spec.mask
        self.assertEqual(mask["shape"], list(labels.shape))
        decoded = np.frombuffer(
            base64.b64decode(mask["labels_b64"]), dtype=np.uint32
        ).reshape(labels.shape)
        np.testing.assert_array_equal(decoded, labels)
        self.assertEqual(len(mask["palette"]), 2)
        self.assertAlmostEqual(mask["alpha"], 0.4)
        self.assertGreater(mask["surface"]["n_faces"], 0)

    def test_label_surfaces_are_provided_for_small_label_count(self) -> None:
        labels = self._labeled()
        spec = build_viewer_payload(
            [_sphere_stack(labels.shape)],
            (1.0, 1.0, 1.0),
            mask_overlay=labels,
        )
        self.assertIsNotNone(spec.mask)
        self.assertEqual(len(spec.mask["label_surfaces"]), 2)
        for ls in spec.mask["label_surfaces"]:
            self.assertIn("color", ls)
            self.assertEqual(len(ls["color"]), 3)
            self.assertGreater(ls["n_faces"], 0)

    def test_empty_mask_yields_no_overlay(self) -> None:
        spec = build_viewer_payload(
            [_sphere_stack((16, 20, 24))],
            (1.0, 1.0, 1.0),
            mask_overlay=np.zeros((16, 20, 24), dtype=np.uint8),
        )
        self.assertIsNone(spec.mask)

    def test_mask_shape_must_match_volume(self) -> None:
        with self.assertRaises(ValueError):
            build_viewer_payload(
                [_sphere_stack((16, 20, 24))],
                (1.0, 1.0, 1.0),
                mask_overlay=np.zeros((8, 8, 8), dtype=np.uint8),
            )

    def test_overlay_survives_json_roundtrip(self) -> None:
        labels = self._labeled()
        spec = build_viewer_payload(
            [_sphere_stack(labels.shape)],
            (1.0, 1.0, 1.0),
            mask_overlay=labels,
        )
        d = spec_to_dict(spec)
        loaded = json.loads(json.dumps(d))
        self.assertEqual(loaded["mask"]["shape"], list(labels.shape))
        self.assertEqual(len(loaded["mask"]["palette"]), 2)

    def test_sparse_labels_use_local_surfaces(self) -> None:
        labels = np.zeros((12, 16, 18), dtype=np.uint32)
        labels[1:5, 2:7, 3:8] = 7
        labels[7:11, 9:14, 10:16] = 1_000_000
        spec = build_viewer_payload([_sphere_stack(labels.shape)], (1.0, 1.0, 1.0),
                                    mask_overlay=labels)
        self.assertEqual(len(spec.mask["label_surfaces"]), 2)

    def test_overlay_rejects_invalid_labels_and_alpha(self) -> None:
        volume = _sphere_stack((8, 10, 12))
        with self.assertRaisesRegex(ValueError, "integer"):
            build_viewer_payload([volume], (1.0, 1.0, 1.0),
                                 mask_overlay=np.ones(volume.shape, np.float32))
        with self.assertRaisesRegex(ValueError, "alpha"):
            build_viewer_payload([volume], (1.0, 1.0, 1.0),
                                 mask_overlay=np.ones(volume.shape, np.uint8),
                                 mask_overlay_alpha=1.1)

    def test_overlay_alpha_one_remains_opaque(self) -> None:
        labels = self._labeled()
        spec = build_viewer_payload([_sphere_stack(labels.shape)], (1.0, 1.0, 1.0),
                                    mask_overlay=labels, mask_overlay_alpha=1.0)
        self.assertEqual(spec.mask["alpha"], 1.0)


class HtmlTests(unittest.TestCase):
    def test_html_is_standalone_and_singleton_script(self) -> None:
        vols = [_sphere_stack((8, 10, 12))]
        spec = build_viewer_payload(vols, (1.0, 1.0, 2.0))
        html = render_viewer_html(spec)
        # self-contained (no network deps) and a single merged <script> block.
        # Count structural CLOSING tags: the vendored UMD bundle contains a JS
        # string literal "<script>" that would inflate an opening-tag count, but
        # real elements close with "</script>" (the bundle escapes it as <\\/script>).
        self.assertIn("</html>", html)
        self.assertEqual(html.count("</script>"), 1)
        # vendored UMD is present
        self.assertIn("window.vtk", html) or self.assertIn("globalThis", html)


if __name__ == "__main__":
    unittest.main(verbosity=2)
