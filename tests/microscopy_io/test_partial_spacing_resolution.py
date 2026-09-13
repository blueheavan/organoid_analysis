"""Regression tests for per-axis physical-spacing resolution.

The defect these pin: the OME reader used to return ``None`` whenever *any*
axis was unstated, discarding the axes the file did calibrate. Downstream, the
manifest-versus-metadata comparison was guarded on complete metadata, so a
partially calibrated file took the manifest's numbers for *every* axis without
any comparison -- including axes where the file stated a conflicting value.
That silently rescales volume, surface area, axis lengths and sphericity.

Completion and comparison are separate rules and both are asserted here:
an axis the file omits may be completed from the manifest or CLI; an axis the
file states is compared against any external value for that axis and a
disagreement is a rejection. Neither rule may be skipped because the *other*
axes are incomplete.
"""

import numpy as np
import pytest
import tifffile

from organoid_analysis.config import load_config
from organoid_analysis.microscopy_io.metadata import (
    merge_axis_spacings,
    resolve_zyx_spacing,
    spacing_axis_conflicts,
)
from organoid_analysis.microscopy_io.tiff_contract import (
    load_sample,
    load_truth_labels,
    ome_axis_spacing,
    ome_spacing,
)

MATCHING = {"PhysicalSizeZ": 2.0, "PhysicalSizeY": 1.0, "PhysicalSizeX": 1.0}
PARTIAL_XY_ONLY = {"PhysicalSizeY": 1.0, "PhysicalSizeX": 1.0}


def write_stack(path, metadata, shape=(7, 12, 12)):
    tifffile.imwrite(path, np.zeros(shape, np.uint16), ome=True, photometric="minisblack",
                     metadata={"axes": "ZYX", **metadata})
    return path


def manifest_row(path, **spacing):
    row = {"image_path": str(path), "time_index": "", "series_index": ""}
    for axis in "zyx":
        row[f"spacing_{axis}_um"] = spacing.get(axis, "")
    return row


# ---------------------------------------------------------------- case 1 of 6
def test_complete_metadata_agreeing_with_the_manifest_is_accepted(tmp_path):
    """Both sources state all three axes and agree: measured values are used."""
    path = write_stack(tmp_path / "full.ome.tif", MATCHING)
    sample = load_sample(manifest_row(path, z="2", y="1", x="1"), load_config())
    assert sample.spacing == (2.0, 1.0, 1.0)
    assert sample.metadata["spacing_source_by_axis"] == {
        "z": "OME+manifest", "y": "OME+manifest", "x": "OME+manifest"}


# ---------------------------------------------------------------- case 2 of 6
def test_complete_metadata_conflicting_with_the_manifest_is_rejected(tmp_path):
    path = write_stack(tmp_path / "full.ome.tif", MATCHING)
    with pytest.raises(ValueError, match="conflicts with OME"):
        load_sample(manifest_row(path, z="3", y="1", x="1"), load_config())


# ---------------------------------------------------------------- case 3 of 6
def test_partial_metadata_is_completed_without_losing_the_stated_axes(tmp_path):
    """Z absent from OME, supplied by the manifest; OME's Y/X still win.

    The manifest here states a Y/X that is *within tolerance* of the file's, so
    the run is legitimate; what is asserted is that the resolved Y/X are the
    file's measured values and that the provenance shows Z came from elsewhere.
    """
    path = write_stack(tmp_path / "partial.ome.tif", PARTIAL_XY_ONLY)
    assert ome_axis_spacing(tifffile.TiffFile(path).ome_metadata) == (None, 1.0, 1.0)
    sample = load_sample(manifest_row(path, z="2"), load_config())
    assert sample.spacing == (2.0, 1.0, 1.0)
    assert sample.metadata["spacing_source_by_axis"] == {"z": "manifest", "y": "OME", "x": "OME"}
    assert sample.metadata["spacing_source"] == "mixed:OME+manifest"


# ---------------------------------------------------------------- case 4 of 6
def test_partial_metadata_conflicting_on_a_stated_axis_is_rejected(tmp_path):
    """The case that used to pass silently.

    OME states Y and X but not Z. The manifest states all three and disagrees
    with the file on Y. Previously the missing Z collapsed the whole OME
    spacing to unknown, the comparison was skipped, and the manifest's wrong Y
    was used. It must be rejected, naming the axis.
    """
    path = write_stack(tmp_path / "partial.ome.tif", PARTIAL_XY_ONLY)
    with pytest.raises(ValueError, match="conflicts with OME spacing on y"):
        load_sample(manifest_row(path, z="2", y="4", x="1"), load_config())


# ---------------------------------------------------------------- case 5 of 6
def test_insufficient_combined_metadata_is_rejected_naming_the_axis(tmp_path):
    """Neither source resolves Z: refuse, rather than assume or ignore."""
    path = write_stack(tmp_path / "partial.ome.tif", PARTIAL_XY_ONLY)
    with pytest.raises(ValueError, match="unresolved for z"):
        load_sample(manifest_row(path, y="1", x="1"), load_config())


# ---------------------------------------------------------------- case 6 of 6
@pytest.mark.parametrize("stated", [{"z": "0"}, {"z": "-2"}, {"z": "nan"}, {"z": "inf"}])
def test_invalid_manifest_spacing_is_rejected(tmp_path, stated):
    path = write_stack(tmp_path / "partial.ome.tif", PARTIAL_XY_ONLY)
    with pytest.raises(ValueError, match="not a positive finite length"):
        load_sample(manifest_row(path, **stated), load_config())


def test_invalid_metadata_spacing_is_rejected(tmp_path):
    """A non-positive physical size in the file itself is still rejected."""
    xml = ('<OME><Image><Pixels PhysicalSizeZ="0" PhysicalSizeY="1" '
           'PhysicalSizeX="1"/></Image></OME>')
    with pytest.raises(ValueError, match="Invalid OME physical spacing"):
        ome_axis_spacing(xml)


# -------------------------------------------------------- resolution contract
def test_absence_is_neither_agreement_nor_conflict():
    assert spacing_axis_conflicts((None, 1.0, 1.0), (2.0, 1.0, 1.0)) == []
    assert spacing_axis_conflicts((None, 4.0, 1.0), (2.0, 1.0, 1.0)) == ["y"]


def test_resolution_never_returns_an_unverified_or_incomplete_spacing():
    resolved = resolve_zyx_spacing((None, 1.0, 1.0), (2.0, 1.0, 1.0))
    assert resolved.zyx == (2.0, 1.0, 1.0)
    assert resolved.provenance == ("manifest", "OME+manifest", "OME+manifest")
    with pytest.raises(ValueError, match="unresolved"):
        resolve_zyx_spacing((None, 1.0, 1.0), (None, None, None))
    with pytest.raises(ValueError, match="conflicts"):
        resolve_zyx_spacing((2.0, 1.0, None), (4.0, 1.0, 1.0))


def test_the_collapsing_reader_is_only_a_view_of_the_per_axis_reader():
    """``ome_spacing()`` keeps its complete-or-None contract for old callers."""
    xml = '<OME><Image><Pixels PhysicalSizeY="1" PhysicalSizeX="1"/></Image></OME>'
    assert ome_spacing(xml) is None
    assert ome_axis_spacing(xml) == (None, 1.0, 1.0)


def test_merging_several_sources_compares_against_the_running_union():
    """A conflict between two sources is caught even if the first omits the axis."""
    with pytest.raises(ValueError, match="differ on z"):
        merge_axis_spacings([(None, 1.0, 1.0), (2.0, 1.0, 1.0), (4.0, 1.0, 1.0)])
    assert merge_axis_spacings([(None, 1.0, None), (2.0, None, 3.0)]) == (2.0, 1.0, 3.0)


# ------------------------------------------- companion files: same rule applies
def test_partially_calibrated_companion_file_is_still_checked(tmp_path):
    """A secondary TIFF stating only Y/X is compared on Y/X, not skipped."""
    primary = write_stack(tmp_path / "primary.ome.tif", MATCHING)
    companion = write_stack(tmp_path / "calcein.ome.tif", {"PhysicalSizeY": 4.0, "PhysicalSizeX": 1.0})
    row = manifest_row(primary)
    row["calcein_path"] = str(companion)
    with pytest.raises(ValueError, match="calcein TIFF spacing differs from the primary image on y"):
        load_sample(row, load_config())


def test_partially_calibrated_companion_file_records_partial_verification(tmp_path):
    primary = write_stack(tmp_path / "primary.ome.tif", MATCHING)
    companion = write_stack(tmp_path / "calcein.ome.tif", {"PhysicalSizeY": 1.0, "PhysicalSizeX": 1.0})
    row = manifest_row(primary)
    row["calcein_path"] = str(companion)
    sample = load_sample(row, load_config())
    assert sample.metadata["channel_grid_verification"]["calcein"] == \
        "spacing_metadata_matches_primary_on_yx_only"


def test_partially_calibrated_annotation_is_still_checked(tmp_path):
    """Validation annotations get the same per-axis comparison."""
    annotation = write_stack(tmp_path / "truth.ome.tif", {"PhysicalSizeY": 4.0, "PhysicalSizeX": 1.0},
                             shape=(7, 12, 12))
    row = {"truth_labels_path": str(annotation), "truth_labels_axes": "", "time_index": ""}
    with pytest.raises(ValueError, match="differs from the primary image on y"):
        load_truth_labels(row, (7, 12, 12), (2.0, 1.0, 1.0))
