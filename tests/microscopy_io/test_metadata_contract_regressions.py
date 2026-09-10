"""Analytical unit and time identity oracles for the September local audit."""
from __future__ import annotations

import numpy as np
import pytest
import tifffile

from organoid_analysis.microscopy_io.metadata import parse_spacing_imagej, parse_spacing_ome
from organoid_analysis.microscopy_io.tiff_contract import ome_spacing, write_labels
from organoid_analysis.microscopy_io.zstack_reader import load_zstack
from organoid_analysis.segmentation.cellpose_inference import read_stack_multichannel


def plane_xml(positions, unit='µm', size_z='2'):
    planes = ''.join(
        f'<Plane TheZ="{i}" TheC="0" TheT="0" PositionZ="{z}"'
        + (f' PositionZUnit="{unit}"' if unit else '') + '/>'
        for i, z in enumerate(positions)
    )
    return (
        '<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">'
        '<Image ID="Image:0"><Pixels ID="Pixels:0" DimensionOrder="XYZCT" '
        'Type="uint16" SizeX="4" SizeY="4" SizeZ="3" SizeC="1" SizeT="1" '
        f'PhysicalSizeX="1" PhysicalSizeY="1" PhysicalSizeZ="{size_z}">'
        + planes + '</Pixels></Image></OME>'
    )


@pytest.mark.parametrize('unit,value', [('nm', 2000), ('mm', .002), ('um', 2), ('micron', 2)])
def test_imagej_z_units_are_converted(unit, value):
    assert parse_spacing_imagej({'spacing': value, 'unit': unit}).z == pytest.approx(2, rel=1e-12)


@pytest.mark.parametrize('unit', [None, '', 'pixel', 'pixels'])
def test_uncalibrated_imagej_spacing_is_unknown(unit):
    assert parse_spacing_imagej({'spacing': 2, 'unit': unit}).z is None


@pytest.mark.parametrize('value', [-1, 0, np.inf, np.nan])
def test_invalid_physical_imagej_spacing_rejected(value):
    with pytest.raises(ValueError):
        parse_spacing_imagej({'spacing': value, 'unit': 'um'})


def test_imagej_unit_calibrates_resolution_none_tiff(tmp_path):
    path = tmp_path / 'calibrated.tif'
    volume = np.arange(3 * 6 * 8, dtype=np.uint16).reshape(3, 6, 8)
    tifffile.imwrite(path, volume, imagej=True, resolution=(2, 4),
                     metadata={'axes': 'ZYX', 'unit': 'mm', 'spacing': .002})
    result = load_zstack(path)
    assert (result.spacing.x, result.spacing.y, result.spacing.z) == pytest.approx((500, 250, 2))
    np.testing.assert_array_equal(result.volume, volume)


@pytest.mark.parametrize('reader', [parse_spacing_ome, ome_spacing])
@pytest.mark.parametrize('positions', [[0, 2, 0], [0, 4], [0, 0, 0]])
def test_nonuniform_or_nonmonotonic_z_rejected(reader, positions):
    with pytest.raises(ValueError, match='uniformly'):
        reader(plane_xml(positions))


@pytest.mark.parametrize('reader', [parse_spacing_ome, ome_spacing])
@pytest.mark.parametrize('unit', [None, 'reference frame'])
def test_reference_frame_is_not_silently_micrometres(reader, unit):
    with pytest.raises(ValueError, match='unit'):
        reader(plane_xml([0, 2, 4], unit))


@pytest.mark.parametrize('reader', [parse_spacing_ome, ome_spacing])
def test_uniform_descending_physical_z_is_supported(reader):
    result = reader(plane_xml([4, 2, 0]))
    assert (result.z if hasattr(result, 'z') else result[0]) == 2


def test_missing_xy_spacing_does_not_bypass_nonuniform_z_check():
    xml = plane_xml([0, 2, 9]).replace(' PhysicalSizeX="1"', '')
    with pytest.raises(ValueError, match='uniformly'):
        ome_spacing(xml)


def test_multiframe_upload_requires_explicit_time_and_preserves_selected_identity(tmp_path):
    path = tmp_path / 'time.ome.tif'
    volume = np.zeros((2, 3, 6, 8), np.uint16)
    volume[1] = 17
    tifffile.imwrite(path, volume, ome=True, photometric='minisblack', metadata={'axes': 'TZYX'})
    for reader in (load_zstack, read_stack_multichannel):
        with pytest.raises(ValueError, match='time_index'):
            reader(path)
        selected = reader(path, time_index=1)
        np.testing.assert_array_equal(selected.volume, volume[1])
        assert selected.meta['time_index'] == 1
        with pytest.raises(ValueError, match='time_index'):
            reader(path, time_index=2)


@pytest.mark.parametrize('labels', [
    np.full((3, 4, 5), 2**32, np.uint64),
    np.full((3, 4, 5), -1, np.int32),
    np.full((3, 4, 5), 1.5),
])
def test_label_export_rejects_lossy_cast_before_writing(tmp_path, labels):
    path = tmp_path / 'labels.ome.tif'
    with pytest.raises(ValueError):
        write_labels(path, labels, (2, 1, 1))
    assert not path.exists()


def test_uint32_max_id_export_roundtrip(tmp_path):
    labels = np.zeros((3, 4, 5), np.uint64)
    labels[1, 1, 1] = 2**32 - 1
    path = tmp_path / 'labels.ome.tif'
    write_labels(path, labels, (2, 1, 1))
    np.testing.assert_array_equal(tifffile.imread(path), labels)


@pytest.mark.parametrize('index', [-1, True, 0.5])
def test_tiff_series_index_must_be_nonnegative_integer(tmp_path, index):
    from organoid_analysis.microscopy_io.tiff_contract import read_tiff
    path = tmp_path / 'series.tif'
    tifffile.imwrite(path, np.zeros((3, 5, 6), np.uint16), photometric='minisblack', metadata={'axes': 'ZYX'})
    with pytest.raises(ValueError, match='series_index'):
        read_tiff(path, series_index=index)


@pytest.mark.parametrize('index', [-1, True, 0.5])
def test_canonical_time_index_rejects_invalid_scalar(index):
    from organoid_analysis.microscopy_io.tiff_contract import canonical_czyx
    with pytest.raises(ValueError, match='time_index'):
        canonical_czyx(np.zeros((2, 3, 5, 6)), 'TZYX', index)


def test_complex_tiff_cannot_become_real_intensity():
    from organoid_analysis.microscopy_io.tiff_contract import canonical_czyx
    with pytest.raises(ValueError, match='real numeric'):
        canonical_czyx(np.ones((3, 5, 6), complex), 'ZYX')
