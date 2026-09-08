import numpy as np
import pytest

from organoid_analysis.quantification.cellular_measurements import analyze_cells, cell_geometry, cell_neighborhood, pair_and_filter_cells
from organoid_analysis.command_line.organoid_commands import main
from organoid_analysis.microscopy_io.tiff_contract import write_labels


def spheres():
    z, y, x = np.indices((48, 48, 48))
    cells = np.zeros((48, 48, 48), np.uint32)
    nuclei = np.zeros_like(cells)
    for label, center in ((3, 14), (9, 34)):
        cells[(z-center)**2 + (y-center)**2 + (x-center)**2 <= 10**2] = label
        nuclei[(z-center)**2 + (y-center)**2 + (x-center)**2 <= 4**2] = label + 20
    return cells, nuclei


def test_pairing_uses_noncompact_ids_and_complete_nucleus_volume():
    cells, nuclei = spheres()
    filtered_cells, filtered_nuclei, pairing, summary = pair_and_filter_cells(cells, nuclei, (1, 1, 1))
    assert summary == {
        "cells_before": 2, "nuclei_before": 2, "cells_after": 2,
        "nuclei_after": 2, "paired_cells": 2, "unpaired_cells": 0, "removed": {},
    }
    assert set(np.unique(filtered_cells)) == {0, 1, 2}
    assert set(np.unique(filtered_nuclei)) == {0, 1, 2}
    assert pairing["original_cell_id"].tolist() == [3, 9]
    expected = np.count_nonzero(nuclei == 23) / np.count_nonzero(cells == 3)
    assert pairing.loc[0, "nc_ratio"] == pytest.approx(expected)


def test_pairing_applies_volume_ratio_containment_and_nucleusless_qc():
    cells = np.zeros((20, 30, 30), np.uint32)
    nuclei = np.zeros_like(cells)
    cells[2:12, 2:12, 2:12] = 1
    cells[2:12, 15:25, 2:12] = 2
    cells[2:4, 2:4, 20:22] = 3
    nuclei[4:8, 4:8, 4:8] = 11
    nuclei[4:8, 22:28, 4:8] = 12  # Only half is inside cell 2.
    _, _, pairing, summary = pair_and_filter_cells(
        cells, nuclei, (1, 1, 1), min_cell_volume_um3=20,
        min_nucleus_volume_um3=20, min_nucleus_containment=0.75,
    )
    assert summary["removed"] == {
        "nucleus_not_contained": 1, "too_small_cell": 1,
    }
    assert pairing.set_index("original_cell_id").loc[2, "reject_reason"] == "nucleus_not_contained"
    rejected = pairing.set_index("original_cell_id").loc[2]
    assert rejected["original_nucleus_id"] == 0
    assert rejected["candidate_nucleus_id"] == 12
    assert rejected["nucleus_volume_um3"] == np.count_nonzero(nuclei == 12)

    cells[5:15, 15:25, 15:25] = 4
    _, _, pairing, summary = pair_and_filter_cells(
        cells, nuclei, (1, 1, 1), min_cell_volume_um3=20,
        min_nucleus_volume_um3=20, require_nucleus=False,
        min_nucleus_containment=0.75,
    )
    assert pairing.set_index("original_cell_id").loc[4, "kept"]
    assert summary["unpaired_cells"] == 1


def test_same_nucleus_is_never_assigned_to_two_cells():
    cells = np.zeros((8, 12, 12), np.uint32)
    nuclei = np.zeros_like(cells)
    cells[:, :, :6] = 1
    cells[:, :, 6:] = 2
    nuclei[2:6, 3:9, 4:8] = 7
    _, filtered_nuclei, pairing, summary = pair_and_filter_cells(
        cells, nuclei, (1, 1, 1), min_cell_volume_um3=1,
        min_nucleus_volume_um3=1, min_nucleus_containment=0.4,
    )
    assert summary["paired_cells"] == 1
    assert summary["removed"] == {"nucleus_already_paired": 1}
    assert set(np.unique(filtered_nuclei)) == {0, 1}
    assert pairing["nucleus_id"].ne(0).sum() == 1
    assert pairing["candidate_nucleus_id"].tolist() == [7, 7]  # Audit candidate retained for rejection.


def test_matching_maximizes_number_of_pairs_in_crossed_case():
    cells = np.zeros((3, 3, 20), np.uint32)
    nuclei = np.zeros_like(cells)
    cells[:, :, :10] = 1
    cells[:, :, 10:] = 2
    nuclei[:, :, 4:15] = 11
    nuclei[:, :, :4] = 12
    _, _, pairing, summary = pair_and_filter_cells(
        cells, nuclei, (1, 1, 1), min_cell_volume_um3=1,
        min_nucleus_volume_um3=1, min_nucleus_containment=0, max_nc_ratio=2,
    )
    assert summary["paired_cells"] == 2
    assigned = pairing.set_index("original_cell_id")["original_nucleus_id"].to_dict()
    assert assigned == {1: 12, 2: 11}


def test_matching_maximizes_overlap_after_cardinality():
    cells = np.zeros((3, 3, 17), np.uint32)
    nuclei = np.zeros_like(cells)
    counts = [(1, 11, 3), (1, 12, 4), (2, 11, 1), (2, 12, 9)]
    offset = 0
    for cell_id, nucleus_id, count in counts:
        cells[1, 1, offset:offset+count] = cell_id
        nuclei[1, 1, offset:offset+count] = nucleus_id
        offset += count
    _, _, pairing, summary = pair_and_filter_cells(
        cells, nuclei, (1, 1, 1), min_cell_volume_um3=1,
        min_nucleus_volume_um3=1, min_nucleus_containment=0, max_nc_ratio=10,
    )
    assert summary["paired_cells"] == 2
    assigned = pairing.set_index("original_cell_id")["original_nucleus_id"].to_dict()
    assert assigned == {1: 11, 2: 12}


def test_below_threshold_cell_cannot_reserve_shared_nucleus():
    cells = np.zeros((3, 3, 20), np.uint32)
    nuclei = np.zeros_like(cells)
    cells[:, :, :7] = 1
    cells[:, :, 10:] = 2
    nuclei[:, :, :7] = 5
    nuclei[:, :, 10:14] = 5
    _, _, pairing, summary = pair_and_filter_cells(
        cells, nuclei, (1, 1, 1), min_cell_volume_um3=70,
        min_nucleus_volume_um3=1, min_nucleus_containment=0.3,
        max_nc_ratio=2,
    )
    assert summary["paired_cells"] == 1
    assert pairing.set_index("original_cell_id").loc[2, "kept"]


def test_geometry_has_correct_zyx_centroid_and_moment_equivalent_axes():
    spacing = (2.0, 1.0, 1.0)
    z, y, x = np.indices((25, 45, 55))
    cell = ((z-12)*2/6)**2 + ((y-20)/4)**2 + ((x-30)/3)**2 <= 1
    labels = cell.astype(np.uint32)
    measured = cell_geometry(labels, None, spacing).iloc[0]
    assert measured["centroid_z_um"] == pytest.approx(24, abs=0.1)
    assert measured["centroid_y_um"] == pytest.approx(20, abs=0.1)
    assert measured["centroid_x_um"] == pytest.approx(30, abs=0.1)
    assert measured["principal_axis_major_um"] == pytest.approx(12, rel=0.08)
    assert measured["principal_axis_intermediate_um"] == pytest.approx(8, rel=0.12)
    assert measured["principal_axis_minor_um"] == pytest.approx(6, rel=0.15)


def test_geometry_and_neighborhood_accept_large_sparse_label_ids():
    labels = np.zeros((8, 8, 8), np.uint32)
    labels[2:6, 2:6, 2:6] = 1_000_000_000
    measured = cell_geometry(labels, None, (1, 1, 1))
    neighbors = cell_neighborhood(labels, (1, 1, 1), radius_um=2)
    assert measured["cell_id"].tolist() == [1_000_000_000]
    assert neighbors["cell_id"].tolist() == [1_000_000_000]


def test_nucleus_centroid_distance_uses_voxel_face_not_background_center():
    cells = np.zeros((5, 5, 5), np.uint32)
    nuclei = np.zeros_like(cells)
    cells[2, 2, 2] = 1
    nuclei[2, 2, 2] = 1
    measured = cell_geometry(cells, nuclei, (4, 2, 1)).iloc[0]
    assert measured["nucleus_centroid_to_cell_border_um"] == pytest.approx(0.5)


def test_neighborhood_uses_physical_distance_and_reports_edge_qc():
    labels = np.zeros((20, 40, 40), np.uint32)
    labels[8:11, 8:11, 8:11] = 1
    labels[8:11, 8:11, 28:31] = 2
    near = cell_neighborhood(labels, (2, 1, 1), radius_um=25)
    assert near["n_neighbors"].tolist() == [1, 1]
    assert near["nearest_neighbor_distance_um"].tolist() == pytest.approx([20, 20])
    assert not near["neighborhood_complete"].any()
    assert near["local_density_per_mm3"].isna().all()
    far = cell_neighborhood(labels, (2, 1, 1), radius_um=15)
    assert far["n_neighbors"].tolist() == [0, 0]

    interior = np.zeros((20, 60, 60), np.uint32)
    interior[4:7, 29:32, 29:32] = 1
    interior[8:11, 29:32, 29:32] = 2
    anisotropic = cell_neighborhood(interior, (2, 1, 1), radius_um=9)
    assert anisotropic["nearest_neighbor_distance_um"].tolist() == pytest.approx([8, 8])
    assert anisotropic["neighborhood_complete"].all()
    assert anisotropic["local_density_per_mm3"].notna().all()


def test_analyze_cells_merges_pairing_geometry_and_neighborhood():
    cells, nuclei = spheres()
    result = analyze_cells(cells, nuclei, (1, 1, 1), neighbor_radius_um=40)
    assert len(result.features) == 2
    assert result.features["original_cell_id"].tolist() == [3, 9]
    assert result.features["n_neighbors"].tolist() == [1, 1]
    assert result.summary["paired_cells"] == 2
    assert result.summary["qc"]["min_nucleus_containment"] == 0.5


def test_analyze_cells_handles_empty_and_all_rejected_masks():
    empty = np.zeros((5, 5, 5), np.uint32)
    result = analyze_cells(empty, empty, (1, 1, 1))
    assert result.features.empty
    assert result.pairing.empty
    assert result.summary["cells_after"] == 0

    cells = empty.copy()
    cells[2, 2, 2] = 9
    rejected = analyze_cells(cells, empty, (1, 1, 1))
    assert rejected.features.empty
    assert rejected.pairing["reject_reason"].tolist() == ["too_small_cell"]


def test_invalid_inputs_are_rejected():
    labels = np.zeros((3, 3, 3), np.uint32)
    with pytest.raises(ValueError, match="matching 3D"):
        pair_and_filter_cells(labels, labels[0], (1, 1, 1))
    with pytest.raises(ValueError, match="positive finite"):
        pair_and_filter_cells(labels, labels, (1, 0, 1))
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        pair_and_filter_cells(labels, labels, (1, 1, 1), min_nucleus_containment=2)


def test_cells_cli_reads_ome_spacing_and_writes_outputs(tmp_path):
    cells, nuclei = spheres()
    cell_path = tmp_path / "cells.ome.tif"
    nucleus_path = tmp_path / "nuclei.ome.tif"
    write_labels(cell_path, cells, (2, 1, 1))
    write_labels(nucleus_path, nuclei, (2, 1, 1))
    out = tmp_path / "output"
    assert main([
        "cells", "--cell-labels", str(cell_path), "--nucleus-labels", str(nucleus_path),
        "--out", str(out),
    ]) == 0
    assert {path.name for path in out.iterdir()} == {
        "cell_labels.filtered.ome.tif", "nucleus_labels.filtered.ome.tif",
        "cell_features.csv", "pairing_qc.csv", "cell_analysis_summary.json",
    }
