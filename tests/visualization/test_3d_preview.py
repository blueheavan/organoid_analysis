"""Tests for the Z-stack 3D preview module (reader, metadata, rendering)."""

from __future__ import annotations

import os
import tempfile
import unittest

import numpy as np
import pytest
import tifffile

from organoid_analysis.microscopy_io import Spacing, load_zstack
from organoid_analysis.microscopy_io.metadata import (
    _resolution_to_um_per_px,
    classify_axes,
    parse_spacing_imagej,
    parse_spacing_ome,
    reject_rgb,
)
from organoid_analysis.microscopy_io.zstack_reader import _reorder_to_czyx
from organoid_analysis.visualization import (
    build_surface,
    normalize_to_gray,
    render_mip_image,
    render_surface_image,
    render_volume_image,
    surface_to_dict,
)
from organoid_analysis.visualization.capabilities import probe_gpu_mapper


class AxisReorderingTests(unittest.TestCase):
    """Verify canonical (T,C,Z,Y,X) reordering never guesses or swaps axes."""

    def test_zyx_stays_zyx_with_t_c_singletons(self) -> None:
        r = _reorder_to_czyx(np.arange(24).reshape(2, 3, 4), "ZYX")
        self.assertEqual(r.shape, (1, 1, 2, 3, 4))

    def test_cyx_maps_channels_to_c_slot(self) -> None:
        a = np.arange(3 * 4 * 2).reshape(3, 4, 2)
        r = _reorder_to_czyx(a, "CYX")
        self.assertEqual(r.shape, (1, 3, 1, 4, 2))
        # T and Z are size-1; r[T=0, C, Z=0, Y, X] must equal a[C, Y, X]
        self.assertTrue((r[0, :, 0, :, :] == a).all())

    def test_tzcyx_reorders_to_tczyx_values_preserved(self) -> None:
        a = np.zeros((2, 3, 4, 5, 6), dtype=int)  # T,Z,C,Y,X
        for t in range(2):
            for z in range(3):
                for c in range(4):
                    a[t, z, c, :, :] = 1000 * t + 100 * z + c
        r = _reorder_to_czyx(a, "TZCYX")
        self.assertEqual(r.shape, (2, 4, 3, 5, 6))  # T,C,Z,Y,X
        # r[T=t, C=c, Z=z] must equal a[t, z, c]
        self.assertEqual(int(r[1, 3, 2, 0, 0]), int(a[1, 2, 3, 0, 0]))
        # Confirm Z and C are not accidentally swapped
        self.assertEqual(int(r[0, 1, 2, 0, 0]), int(a[0, 2, 1, 0, 0]))

    def test_unknown_filler_axis_rejected(self) -> None:
        with self.assertRaises(ValueError):
            _reorder_to_czyx(np.zeros((2, 2, 3, 4)), "TYXQ")


class SpacingMetadataTests(unittest.TestCase):
    def test_anisotropy_rejects_unequal_xy_spacing(self):
        with self.assertRaisesRegex(ValueError, "equal X/Y pixel sizes"):
            Spacing(0.5, 0.6, 1.0).anisotropy

    def test_ome_physical_sizes_to_um(self) -> None:
        xml = (
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">'
            '<Image ID="Image:0"><Pixels ID="Pixels:0" DimensionOrder="XYZCT" '
            'Type="uint16" SizeX="4" SizeY="3" SizeZ="2" SizeC="1" SizeT="1" '
            'PhysicalSizeX="0.41432" PhysicalSizeY="0.41432" '
            'PhysicalSizeZ="1.5"/></Image></OME>'
        )
        s = parse_spacing_ome(xml)
        self.assertAlmostEqual(s.x, 0.41432)
        self.assertAlmostEqual(s.y, 0.41432)
        self.assertAlmostEqual(s.z, 1.5)

    def test_ome_nm_unit_converted_to_um(self) -> None:
        xml = (
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">'
            '<Image ID="Image:0"><Pixels ID="Pixels:0" DimensionOrder="XYZCT" '
            'Type="uint16" SizeX="4" SizeY="3" SizeZ="2" SizeC="1" SizeT="1" '
            'PhysicalSizeX="500" PhysicalSizeXUnit="nm" PhysicalSizeZ="1500" '
            'PhysicalSizeZUnit="nm"/></Image></OME>'
        )
        s = parse_spacing_ome(xml)
        self.assertAlmostEqual(s.x, 0.5)
        self.assertAlmostEqual(s.z, 1.5)

    def test_ome_absent_returns_unknown(self) -> None:
        s = parse_spacing_ome(None)
        self.assertIsNone(s.x)
        self.assertIsNone(s.z)

    def test_nonuniform_z_positions_are_rejected(self) -> None:
        """Regression: load_zstack()'s OME parsing previously accepted a
        nonuniform Z grid silently (unlike the classical manifest pipeline's
        tiff_contract.ome_spacing(), which already rejected this). A stack
        whose actual plane spacing does not match PhysicalSizeZ must be
        surfaced as an error here too, not silently averaged into one value.
        """
        xml = (
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">'
            '<Image ID="Image:0"><Pixels ID="Pixels:0" DimensionOrder="XYZCT" '
            'Type="uint16" SizeX="4" SizeY="3" SizeZ="3" SizeC="1" SizeT="1" '
            'PhysicalSizeX="1" PhysicalSizeY="1" PhysicalSizeZ="2">'
            '<Plane TheZ="0" TheC="0" TheT="0" PositionZUnit="µm" PositionZ="0"/>'
            '<Plane TheZ="1" TheC="0" TheT="0" PositionZUnit="µm" PositionZ="2"/>'
            '<Plane TheZ="2" TheC="0" TheT="0" PositionZUnit="µm" PositionZ="7"/>'
            '</Pixels></Image></OME>'
        )
        with self.assertRaisesRegex(ValueError, "uniformly"):
            parse_spacing_ome(xml)

    def test_uniform_z_positions_are_accepted(self) -> None:
        xml = (
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">'
            '<Image ID="Image:0"><Pixels ID="Pixels:0" DimensionOrder="XYZCT" '
            'Type="uint16" SizeX="4" SizeY="3" SizeZ="3" SizeC="1" SizeT="1" '
            'PhysicalSizeX="1" PhysicalSizeY="1" PhysicalSizeZ="2">'
            '<Plane TheZ="0" TheC="0" TheT="0" PositionZUnit="µm" PositionZ="0"/>'
            '<Plane TheZ="1" TheC="0" TheT="0" PositionZUnit="µm" PositionZ="2"/>'
            '<Plane TheZ="2" TheC="0" TheT="0" PositionZUnit="µm" PositionZ="4"/>'
            '</Pixels></Image></OME>'
        )
        s = parse_spacing_ome(xml)
        self.assertAlmostEqual(s.z, 2.0)

    def test_imagej_z_step(self) -> None:
        s = parse_spacing_imagej({"spacing": 2.5, "unit": "micron"})
        self.assertEqual(s.z, 2.5)

    def test_resolution_to_um_per_px(self) -> None:
        # unit 3 = centimetre: 100 px/cm -> 10000 um / 100 = 100 um/px
        self.assertAlmostEqual(_resolution_to_um_per_px((100, 1), 3), 100.0)
        # unit 2 = inch: 500 px/inch -> 25400 um / 500 = 50.8 um/px
        self.assertAlmostEqual(_resolution_to_um_per_px((500, 1), 2), 50.8)
        # unit 5 = micrometre: 2 px/um -> 0.5 um/px
        self.assertAlmostEqual(_resolution_to_um_per_px((2, 1), 5), 0.5)

    def test_spacing_anisotropy(self) -> None:
        s = Spacing(0.5, 0.5, 1.5)
        self.assertAlmostEqual(s.anisotropy, 3.0)

    def test_rgb_rejection(self) -> None:
        self.assertTrue(reject_rgb("QYX", 3))
        self.assertTrue(reject_rgb("SYX", 3))
        self.assertFalse(reject_rgb("ZYX", 1))

    def test_classify_axes(self) -> None:
        self.assertEqual(classify_axes("TZCYX")["T"], 0)
        self.assertEqual(classify_axes("TZCYX")["C"], 2)
        self.assertEqual(classify_axes("TZCYX")["Z"], 1)


class ZStackReaderTests(unittest.TestCase):
    def _tmp(self, name: str) -> str:
        d = tempfile.mkdtemp()
        return os.path.join(d, name)

    def test_load_preserves_uint16_and_zorder(self) -> None:
        p = self._tmp("z.ome.tif")
        vol = np.arange(2 * 3 * 4).reshape(2, 3, 4)
        tifffile.imwrite(p, vol.astype(np.uint16), photometric="minisblack",
                         metadata={"axes": "ZYX"}, ome=True)
        zs = load_zstack(p)
        self.assertEqual(zs.volume.shape, (2, 3, 4))
        self.assertEqual(zs.dtype, np.dtype("uint16"))
        self.assertFalse(zs.is_multichannel)
        # Z-order preserved: zs.volume[z] == original slice at z
        self.assertTrue((zs.volume[1] == vol[1]).all())
        self.assertEqual(int(zs.volume[1, 2, 3]), int(vol[1, 2, 3]))

    def test_multichannel_czyx(self) -> None:
        p = self._tmp("c.ome.tif")
        vol = np.arange(2 * 3 * 4 * 5).reshape(2, 3, 4, 5).astype(np.uint16)
        tifffile.imwrite(p, vol, photometric="minisblack", metadata={"axes": "CZYX"}, ome=True)
        zs = load_zstack(p)
        self.assertTrue(zs.is_multichannel)
        self.assertEqual(zs.channels, 2)
        self.assertEqual(zs.volume.shape, (2, 3, 4, 5))
        self.assertEqual(int(zs.volume[1, 2, 3, 4]), int(vol[1, 2, 3, 4]))

    def test_tzcyx_time_series_splits(self) -> None:
        p = self._tmp("t.ome.tif")
        vol = np.arange(2 * 3 * 4 * 5 * 6).reshape(2, 3, 4, 5, 6).astype(np.uint16)
        tifffile.imwrite(p, vol, metadata={"axes": "TZCYX"}, ome=True)
        with self.assertRaisesRegex(ValueError, "time_index"):
            load_zstack(p)
        zs = load_zstack(p, time_index=0)
        self.assertEqual(zs.frames, 2)
        self.assertTrue(zs.is_multichannel)
        self.assertEqual(zs.volume.shape, (4, 3, 5, 6))  # (C,Z,Y,X) of frame 0

    def test_missing_axes_refuses_to_guess(self) -> None:
        p = self._tmp("bad.tif")
        tifffile.imwrite(p, np.zeros((4, 5), dtype=np.uint16),
                         photometric="minisblack", ome=True)
        with self.assertRaises(ValueError):
            load_zstack(p)

    def test_rgb_stack_rejected(self) -> None:
        p = self._tmp("rgb.tif")
        rgb = np.zeros((8, 8, 8, 3), dtype=np.uint8)
        tifffile.imwrite(p, rgb, photometric="rgb")
        with self.assertRaises(ValueError):
            load_zstack(p)


class SpacingResolutionIntegrationTests(unittest.TestCase):
    """Reader + spacing metadata end-to-end (via proper OME-TIFF)."""

    def test_tiff_resolution_round_trip_in_um(self) -> None:
        # OME-TIFF with inch resolution -> reader returns x=y=50.8 um/px.
        p = os.path.join(tempfile.mkdtemp(), "res.ome.tif")
        vol = np.zeros((2, 3, 4), dtype=np.uint16)
        tifffile.imwrite(p, vol, photometric="minisblack", metadata={"axes": "ZYX"},
                         ome=True, resolution=(500, 500, "INCH"))
        zs = load_zstack(p)
        self.assertAlmostEqual(zs.spacing.x, 50.8)
        self.assertAlmostEqual(zs.spacing.y, 50.8)
        self.assertIsNone(zs.spacing.z)

    def test_metadata_missing_no_crash(self) -> None:
        # No OME/ImageJ metadata: tifffile still embeds default resolution tags,
        # so the reader must load without crashing and keep the volume intact.
        p = os.path.join(tempfile.mkdtemp(), "plain.tif")
        vol = np.arange(2 * 3 * 4).reshape(2, 3, 4).astype(np.uint16)
        tifffile.imwrite(p, vol, photometric="minisblack", metadata={"axes": "ZYX"}, ome=False)
        zs = load_zstack(p)
        # tifffile writes a default XResolution with ResolutionUnit=1
        # (undefined), which the reader now treats as unknown spacing:
        # load succeeds without crashing, Z-order intact.
        self.assertIsNone(zs.spacing.x)
        self.assertIsNone(zs.spacing.z)
        self.assertEqual(zs.volume.shape, (2, 3, 4))
        self.assertTrue((zs.volume[1] == vol[1]).all())


@pytest.mark.render
class RenderingTests(unittest.TestCase):
    """Functional rendering checks (offscreen VTK)."""

    def test_volume_render_returns_rgb_image(self) -> None:
        vol = np.zeros((8, 12, 16), dtype=np.uint8)
        vol[2:6, 4:8, 6:10] = 255
        res = render_volume_image(vol, Spacing(0.5, 0.5, 1.5), preset="iso")
        self.assertIn(res.image.shape[-1], (3, 4))
        self.assertGreater(int(res.image[..., 0].max()), 0)

    def test_mip_render_dimensions_match_preset(self) -> None:
        vol = np.zeros((8, 12, 16), dtype=np.uint8)
        vol[2:6, 4:8, 6:10] = 255
        res = render_mip_image(vol, Spacing(0.5, 0.5, 1.5), preset="xz")
        self.assertEqual(res.mode, "mip")
        self.assertGreater(int(res.image[..., 0].max()), 0)

    def test_mip_is_direct_grayscale_projection(self) -> None:
        vol = np.array(
            [
                [[0, 10], [20, 30]],
                [[40, 5], [15, 25]],
            ],
            dtype=np.uint8,
        )
        res = render_mip_image(vol, preset="xy", norm=(0, 40))
        expected = normalize_to_gray(vol.max(axis=0), low=0, high=40)
        self.assertTrue(np.array_equal(res.image[..., 0], expected))
        self.assertTrue(np.array_equal(res.image[..., 0], res.image[..., 1]))
        self.assertTrue(np.array_equal(res.image[..., 1], res.image[..., 2]))

    def test_density_projection_is_direct_grayscale_mean(self) -> None:
        vol = np.array(
            [
                [[0, 10], [20, 30]],
                [[40, 30], [20, 10]],
            ],
            dtype=np.uint8,
        )
        res = render_volume_image(vol, preset="xy", norm=(0, 40))
        expected = normalize_to_gray(vol.mean(axis=0), low=0, high=40)
        self.assertTrue(np.array_equal(res.image[..., 0], expected))
        self.assertTrue(np.array_equal(res.image[..., 0], res.image[..., 2]))

    def test_constant_volume_normalizes_without_dividing_by_zero(self) -> None:
        gray = normalize_to_gray(np.full((2, 3, 4), 7, dtype=np.uint16))
        self.assertTrue(np.array_equal(gray, np.zeros((2, 3, 4), dtype=np.uint8)))

    def test_surface_sphere_has_faces(self) -> None:
        vol = np.zeros((16, 16, 16), dtype=np.uint8)
        for z in range(16):
            for y in range(16):
                for x in range(16):
                    if (x - 8) ** 2 + (y - 8) ** 2 + (z - 8) ** 2 <= 9:
                        vol[z, y, x] = 255
        surfs = build_surface(vol, Spacing(1, 1, 1), threshold=127)
        self.assertEqual(len(surfs), 1)
        self.assertGreater(surfs[0].polydata.n_faces, 0)
        img = render_surface_image(surfs, preset="iso")
        self.assertGreater(int(img[..., 0].max()), 0)

    def test_label_map_two_tangent_labels(self) -> None:
        L = np.zeros((14, 14, 14), dtype=np.uint8)
        L[2:6, 2:6, 2:6] = 1
        L[5:9, 5:9, 5:9] = 2  # tangent along a face
        surfs = build_surface(L, Spacing(1, 1, 1), labels=[1, 2])
        self.assertEqual(len(surfs), 2)
        labels = sorted(s.label for s in surfs)
        self.assertEqual(labels, [1, 2])
        self.assertIsNotNone(surfs[0].color)
        self.assertIsNotNone(surfs[1].color)
        self.assertNotEqual(surfs[0].color, surfs[1].color)
        for s in surfs:
            self.assertGreater(s.polydata.n_faces, 0)

    def test_surface_serializes_to_dict(self) -> None:
        L = np.zeros((10, 10, 10), dtype=np.uint8)
        L[2:6, 2:6, 2:6] = 1
        surfs = build_surface(L, Spacing(1, 1, 1), labels=[1])
        d = surface_to_dict(surfs[0])
        poly = d["polydata"]
        self.assertEqual(len(poly["points"]), poly["n_points"] * 3)
        self.assertEqual(len(poly["faces"]), poly["n_faces"] * 3)
        self.assertEqual(d["label"], 1)


class CapabilityTests(unittest.TestCase):
    def test_gpu_probe_returns_bool(self) -> None:
        result = probe_gpu_mapper()
        self.assertIsInstance(result, bool)
