"""Deterministic software-correctness tests for multilevel 3D analysis."""
from __future__ import annotations

import json
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
from scipy import ndimage as ndi

from organoid_analysis.command_line.organoid_commands import main
from organoid_analysis.microscopy_io.tiff_contract import write_labels
from organoid_analysis.quantification.multilevel_relationships import Multilevel3DConfig
from organoid_analysis.quantification.multilevel_relationships.qc import _volume_outliers
from organoid_analysis.quantification.multilevel_relationships.spatial import cell_spatial_features
from organoid_analysis.result_export.measurement_tables import export_results
from organoid_analysis.workflows.multilevel_measurement_workflow import analyze_multilevel_3d


def synthetic_labels() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Two organoids: direct contacts, anucleate/multinucleated and border QC."""
    shape = (16, 32, 32)
    organoids = np.zeros(shape, np.uint32)
    cells = np.zeros(shape, np.uint32)
    nuclei = np.zeros(shape, np.uint32)
    organoids[2:14, 2:16, 2:16] = 10
    organoids[2:14, 18:32, 18:32] = 20  # touches Y/X image borders
    cells[4:9, 4:9, 4:9] = 101
    cells[4:9, 4:9, 9:14] = 102  # direct X-face contact with 101
    cells[9:13, 10:14, 4:8] = 103  # deliberately anucleate
    cells[5:10, 21:27, 21:27] = 201
    cells[0:3, 16:18, 2:5] = 999  # no organoid parent / border touching
    nuclei[5:7, 5:7, 5:7] = 1001
    nuclei[5:7, 5:7, 10:12] = 1002
    nuclei[6:8, 6:8, 11:13] = 1003  # second nucleus in cell 102
    nuclei[6:8, 22:24, 22:24] = 2001
    nuclei[0:2, 28:30, 2:4] = 9001  # no cell parent
    return organoids, cells, nuclei


def test_multilevel_contract_hierarchy_and_physical_features():
    organoids, cells, nuclei = synthetic_labels()
    result = analyze_multilevel_3d(organoids, cells, nuclei, (2.0, 0.65, 0.65),
                                   metadata={"sample_id": "synthetic", "well_id": "A01", "field_id": "1"})
    assert len(result.organoid_features) == 2
    assert len(result.cell_features) == 5
    assert len(result.nucleus_features) == 5
    cell = result.cell_features.set_index("cell_id")
    assert cell.loc[101, "organoid_id"] == 10
    assert cell.loc[101, "organoid_parent_overlap_fraction"] == pytest.approx(1.0)
    assert cell.loc[101, "cell_parent_overlap_fraction"] == pytest.approx(1.0)
    assert result.nucleus_features.set_index("nucleus_id").loc[1001, "nucleus_parent_overlap_fraction"] == pytest.approx(1.0)
    assert cell.loc[999, "parent_assignment_failed"]
    assert cell.loc[103, "anucleate"]
    assert cell.loc[102, "multinucleated"]
    assert cell.loc[102, "nucleus_count"] == 2
    assert cell.loc[102, "nucleus_to_cell_volume_ratio"] > 0
    assert cell.loc[101, "volume_um3"] == pytest.approx(125 * 2.0 * .65 * .65)
    assert np.isfinite(cell.loc[101, "surface_area_um2"])
    assert np.isfinite(cell.loc[101, "sphericity"])
    assert np.isfinite(cell.loc[101, "normalized_radial_position_equivalent_radius"])
    assert cell.loc[101, "distance_to_organoid_surface_um"] > 0
    assert {"sample_id", "well_id", "field_id"}.issubset(result.cell_features.columns)
    assert {"sample_id", "well_id", "field_id"}.issubset(result.qc_flags.columns)


def test_face_contact_topology_uses_anisotropic_physical_area():
    organoids, cells, nuclei = synthetic_labels()
    result = analyze_multilevel_3d(organoids, cells, nuclei, (2.0, .65, .65))
    edges = result.cell_topology_edges.set_index(["cell_id_1", "cell_id_2"])
    assert edges.loc[(101, 102), "organoid_id"] == 10
    # A 5-by-5 X-normal contact plane has 25 faces, each Z*Y = 2*.65 um2.
    assert edges.loc[(101, 102), "contact_area_um2"] == pytest.approx(25 * 2.0 * .65)
    cell = result.cell_features.set_index("cell_id")
    assert cell.loc[101, "degree"] == 1
    assert cell.loc[101, "total_contact_area_um2"] == pytest.approx(25 * 2.0 * .65)


def test_qc_mad_border_and_parent_assignment_flags_are_formal():
    organoids, cells, nuclei = synthetic_labels()
    result = analyze_multilevel_3d(organoids, cells, nuclei, (2, .65, .65),
                                   config=Multilevel3DConfig(minimum_voxels=10))
    flags = set(map(tuple, result.qc_flags[["object_type", "object_id", "flag"]].to_numpy()))
    assert ("organoid", 20, "touches_image_border") in flags
    assert ("cell", 999, "parent_assignment_failed") in flags
    assert ("cell", 103, "anucleate") in flags
    assert ("cell", 102, "multinucleated") in flags
    assert ("nucleus", 9001, "parent_assignment_failed") in flags


def test_nucleus_inherits_organoid_from_assigned_cell_and_retains_direct_overlap_audit():
    organoids = np.zeros((5, 5, 6), dtype=np.uint32)
    cells = np.zeros_like(organoids)
    nuclei = np.zeros_like(organoids)
    organoids[1:4, 1:4, 1:3] = 1
    organoids[1:4, 1:4, 3:5] = 2
    cells[2, 2, 1:4] = 10  # two voxels in organoid 1, one in organoid 2
    nuclei[2, 2, 3] = 20  # direct overlap is organoid 2, but its cell is in 1

    result = analyze_multilevel_3d(organoids, cells, nuclei, (1.0, 1.0, 1.0))
    nucleus = result.nucleus_features.set_index("nucleus_id").loc[20]
    assert nucleus["cell_id"] == 10
    assert nucleus["organoid_id"] == 1
    assert nucleus["direct_organoid_id"] == 2
    assert nucleus["direct_organoid_parent_mismatch"]
    counts = result.organoid_features.set_index("organoid_id")["nucleus_count"]
    assert counts.loc[1] == 1 and counts.loc[2] == 0
    flags = set(map(tuple, result.qc_flags[["object_type", "object_id", "flag"]].to_numpy()))
    assert ("nucleus", 20, "direct_organoid_parent_mismatch") in flags


def test_raw_label_geometry_keeps_sphericity_consistent_with_exported_volume_and_surface():
    organoids = np.zeros((7, 7, 7), dtype=np.uint32)
    organoids[2:5, 2:5, 2:5] = 1
    organoids[3, 3, 3] = 0  # a cavity must not make volume and sphericity disagree
    empty = np.zeros_like(organoids)
    result = analyze_multilevel_3d(organoids, empty, empty, (1.0, 1.0, 1.0))
    feature = result.organoid_features.iloc[0]
    expected = np.cbrt(np.pi) * (6.0 * feature.volume_um3) ** (2.0 / 3.0) / feature.surface_area_um2
    assert feature.voxel_count == 26
    assert feature.sphericity == pytest.approx(expected)
    assert feature.sphericity <= 1.0


def test_cell_intensity_features_pool_all_assigned_nuclear_voxels():
    organoids, cells, nuclei = synthetic_labels()
    intensity = np.zeros_like(nuclei, dtype=float)
    intensity[nuclei == 1002] = 2.0
    intensity[nuclei == 1003] = 6.0
    result = analyze_multilevel_3d(organoids, cells, nuclei, (2.0, .65, .65), nucleus_intensity=intensity)
    cell = result.cell_features.set_index("cell_id").loc[102]
    n2 = int((nuclei == 1002).sum())
    n3 = int((nuclei == 1003).sum())
    expected_mean = (2.0 * n2 + 6.0 * n3) / (n2 + n3)
    expected_std = np.sqrt((n2 * (2.0 - expected_mean) ** 2 + n3 * (6.0 - expected_mean) ** 2) / (n2 + n3))
    assert cell.mean_nuclear_intensity == pytest.approx(expected_mean)
    assert cell.min_nuclear_intensity == pytest.approx(2.0)
    assert cell.max_nuclear_intensity == pytest.approx(6.0)
    assert cell.integrated_nuclear_intensity == pytest.approx(2.0 * n2 + 6.0 * n3)
    assert cell.std_nuclear_intensity == pytest.approx(expected_std)
    assert cell.CV_chromatin == pytest.approx(expected_std / expected_mean)


def test_zero_mad_still_flags_a_clear_volume_outlier():
    outliers = _volume_outliers(pd.Series([100.0, 100.0, 100.0, 10_000.0]), 3.5)
    assert outliers.tolist() == [False, False, False, True]


def test_volume_outliers_nonzero_mad_matches_hand_computed_robust_z():
    """Regression test for a P3 audit finding: the primary (nonzero-MAD)
    numeric path previously had no independent-oracle test, only the
    mad==0 special case above."""
    values = pd.Series([10.0, 12.0, 11.0, 13.0, 100.0])
    median, mad = 12.0, 1.0  # median(|10-12|,|12-12|,|11-12|,|13-12|,|100-12|) = median(2,0,1,1,88) = 1
    expected_robust_z = 0.67448975 * (values - median) / mad
    outliers = _volume_outliers(values, threshold=3.5)
    assert outliers.tolist() == (expected_robust_z.abs() > 3.5).tolist()
    assert outliers.tolist() == [False, False, False, False, True]
    assert expected_robust_z.iloc[-1] == pytest.approx(0.67448975 * 88.0)


def test_spatial_depth_map_is_cached_per_organoid():
    labels = np.zeros((7, 7, 7), dtype=np.uint32)
    labels[1:6, 1:6, 1:6] = 1
    organoids = pd.DataFrame([{
        "organoid_id": 1, "volume_um3": 125.0 * 2.0 * .65 * .65,
        "centroid_z_um": 6.0, "centroid_y_um": 1.95, "centroid_x_um": 1.95,
    }])
    cells = pd.DataFrame([
        {"cell_id": 1, "organoid_id": 1, "centroid_z_um": 4.0, "centroid_y_um": 1.3, "centroid_x_um": 1.3},
        {"cell_id": 2, "organoid_id": 1, "centroid_z_um": 8.0, "centroid_y_um": 2.6, "centroid_x_um": 2.6},
    ])
    with patch("organoid_analysis.quantification.multilevel_relationships.spatial.ndi.distance_transform_edt", wraps=ndi.distance_transform_edt) as edt:
        features = cell_spatial_features(cells, organoids, labels, (2.0, 0.65, 0.65))
    assert edt.call_count == 1
    assert features.distance_to_organoid_surface_um.gt(0).all()


def test_empty_single_tiny_and_invalid_labels_are_handled():
    empty = np.zeros((3, 3, 3), dtype=np.uint32)
    result = analyze_multilevel_3d(empty, empty, empty, (1, 1, 1))
    assert result.organoid_features.empty and result.cell_features.empty and result.qc_flags.empty
    single = empty.copy()
    single[1, 1, 1] = 1
    result = analyze_multilevel_3d(single, single, single, (1, 1, 1))
    assert result.cell_features.loc[0, "too_small"]
    assert result.cell_features.loc[0, "voxel_count"] == 1
    with pytest.raises(ValueError, match="does not match"):
        analyze_multilevel_3d(empty, empty[:, :, :2], empty, (1, 1, 1))
    with pytest.raises(ValueError, match="positive finite"):
        analyze_multilevel_3d(empty, empty, empty, (1, 0, 1))


def test_bool_labels_rejected():
    """VR-1: boolean label masks must be rejected (not silently coerced)."""
    bool_mask = np.zeros((3, 3, 3), dtype=bool)
    bool_mask[1, 1, 1] = True
    uint_mask = np.zeros((3, 3, 3), dtype=np.uint32)
    uint_mask[1, 1, 1] = 1
    with pytest.raises(ValueError, match="nonnegative integer instance IDs"):
        analyze_multilevel_3d(bool_mask, uint_mask, uint_mask, (1, 1, 1))
    with pytest.raises(ValueError, match="nonnegative integer instance IDs"):
        analyze_multilevel_3d(uint_mask, bool_mask, uint_mask, (1, 1, 1))
    with pytest.raises(ValueError, match="nonnegative integer instance IDs"):
        analyze_multilevel_3d(uint_mask, uint_mask, bool_mask, (1, 1, 1))


def test_negative_labels_rejected():
    """VR-1: signed-integer labels with negative values must be rejected."""
    neg = np.zeros((3, 3, 3), dtype=np.int32)
    neg[1, 1, 1] = -1
    pos = np.zeros((3, 3, 3), dtype=np.uint32)
    pos[1, 1, 1] = 1
    with pytest.raises(ValueError, match="nonnegative integer instance IDs"):
        analyze_multilevel_3d(neg, pos, pos, (1, 1, 1))
    with pytest.raises(ValueError, match="nonnegative integer instance IDs"):
        analyze_multilevel_3d(pos, neg, pos, (1, 1, 1))
    with pytest.raises(ValueError, match="nonnegative integer instance IDs"):
        analyze_multilevel_3d(pos, pos, neg, (1, 1, 1))


def test_nonfinite_nucleus_intensity_rejected():
    """VR-1: nucleus_intensity with NaN/Inf must be rejected."""
    labels = np.zeros((3, 3, 3), dtype=np.uint32)
    labels[1, 1, 1] = 1
    intensity_nan = labels.astype(float)
    intensity_nan[1, 1, 1] = np.nan
    intensity_inf = labels.astype(float)
    intensity_inf[1, 1, 1] = np.inf
    with pytest.raises(ValueError, match="finite numeric data"):
        analyze_multilevel_3d(labels, labels, labels, (1, 1, 1), nucleus_intensity=intensity_nan)
    with pytest.raises(ValueError, match="finite numeric data"):
        analyze_multilevel_3d(labels, labels, labels, (1, 1, 1), nucleus_intensity=intensity_inf)


def test_parquet_export_summary_and_cli_round_trip(tmp_path):
    organoids, cells, nuclei = synthetic_labels()
    result = analyze_multilevel_3d(organoids, cells, nuclei, (2, .65, .65), metadata={"sample_id": "synthetic"})
    output = tmp_path / "Result"
    paths = export_results(output, organoids=result.organoid_features, cells=result.cell_features,
                           nuclei=result.nucleus_features, edges=result.cell_topology_edges,
                           qc_flags=result.qc_flags, summary=result.summary)
    restored = pd.read_parquet(paths["cell_features"])
    assert len(restored) == 5
    assert restored.set_index("cell_id").loc[101, "organoid_id"] == 10
    numeric = restored.select_dtypes("number").to_numpy()
    assert not np.isinf(numeric).any()
    summary = json.loads((output / "summary" / "analysis_summary.json").read_text())
    assert summary["input_shape_zyx"] == [16, 32, 32]

    organoid_path, cell_path, nucleus_path = (tmp_path / name for name in ("organoids.ome.tif", "cells.ome.tif", "nuclei.ome.tif"))
    write_labels(organoid_path, organoids, (2, .65, .65))
    write_labels(cell_path, cells, (2, .65, .65))
    write_labels(nucleus_path, nuclei, (2, .65, .65))
    cli_out = tmp_path / "cli_result"
    assert main(["analyze-3d", "--organoid-labels", str(organoid_path), "--cell-labels", str(cell_path),
                 "--nucleus-labels", str(nucleus_path), "--sample-id", "synthetic", "--out", str(cli_out)]) == 0
    assert (cli_out / "features" / "organoid_features.parquet").exists()
    assert (cli_out / "masks" / "cell_labels.ome.tif").exists()

    # Regression check for the P2 audit finding: the analyze-3d CLI route
    # previously exported zero code/environment provenance, unlike the
    # classical analyze() pipeline.
    cli_summary = json.loads((cli_out / "summary" / "analysis_summary.json").read_text())
    provenance = cli_summary["provenance"]
    assert provenance["source_code_sha256"]  # nonempty: multilevel3d/*.py hashes present
    assert provenance["packages"]["numpy"]
    assert "python" in provenance and "platform" in provenance


def test_multilevel_determinism_two_run_identical():
    """VR-9: pipeline must be deterministic; two identical runs produce byte-identical outputs."""
    organoids, cells, nuclei = synthetic_labels()
    spacing = (2.0, 0.65, 0.65)
    meta = {"sample_id": "det", "well_id": "B01", "field_id": "2"}
    r1 = analyze_multilevel_3d(organoids, cells, nuclei, spacing, metadata=meta)
    r2 = analyze_multilevel_3d(organoids, cells, nuclei, spacing, metadata=meta)
    pd.testing.assert_frame_equal(r1.organoid_features, r2.organoid_features)
    pd.testing.assert_frame_equal(r1.cell_features, r2.cell_features)
    pd.testing.assert_frame_equal(r1.nucleus_features, r2.nucleus_features)
    pd.testing.assert_frame_equal(r1.cell_topology_edges, r2.cell_topology_edges)
    pd.testing.assert_frame_equal(r1.qc_flags, r2.qc_flags)
    # summary contains runtime_seconds which differs per run; compare without it
    s1 = {k: v for k, v in r1.summary.items() if k != "runtime_seconds"}
    s2 = {k: v for k, v in r2.summary.items() if k != "runtime_seconds"}
    assert s1 == s2
