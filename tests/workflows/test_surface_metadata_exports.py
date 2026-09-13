"""Domain metadata must survive workflow joins and exported artifacts."""
from __future__ import annotations

import io
import json
import zipfile

import numpy as np
import pandas as pd

from organoid_analysis.quantification.features import SURFACE_METADATA_COLUMNS
from organoid_analysis.quantification.mask_features import extract_mask_features
from organoid_analysis.result_export.mask_feature_bundle import build_mask_feature_bundle
from organoid_analysis.result_export.measurement_tables import export_results
from organoid_analysis.workflows.multilevel_measurement_workflow import analyze_multilevel_3d


def test_multilevel_exports_preserve_domain_metadata_for_every_level(tmp_path):
    labels = np.zeros((5, 5, 5), np.uint32)
    labels[2, 2, 2] = 71
    result = analyze_multilevel_3d(labels, labels, labels, (2., 1., 1.))
    paths = export_results(tmp_path / "result", organoids=result.organoid_features,
                           cells=result.cell_features, nuclei=result.nucleus_features,
                           edges=result.cell_topology_edges, qc_flags=result.qc_flags,
                           summary=result.summary)
    for level in ("organoid", "cell", "nucleus"):
        frame = getattr(result, f"{level}_features")
        restored = pd.read_parquet(paths[f"{level}_features"])
        pd.testing.assert_frame_equal(frame, restored)
        row = restored.iloc[0]
        assert set(SURFACE_METADATA_COLUMNS).issubset(restored)
        assert row.volume_um3 == 2.
        assert row.equivalent_diameter_um > 0
        assert row.minor_axis_um > 0
        assert not row.surface_in_qualified_domain
        assert "surface_outside_qualified_domain" in row.qc_flags
        assert result.summary["measurement_policies"][level]["segmentation_validation_status"] == "NOT ASSESSED"


def test_nuclear_only_workflow_and_csv_keep_small_objects():
    labels = np.zeros((5, 5, 5), np.uint32)
    labels[2, 2, 2] = 71
    empty = np.zeros_like(labels)
    result = analyze_multilevel_3d(empty, empty, labels, (2., 1., 1.))
    assert result.organoid_features.empty and result.cell_features.empty
    row = result.nucleus_features.iloc[0]
    assert row.nucleus_id == 71 and row.volume_um3 == 2.
    assert row.parent_assignment_failed and not row.surface_in_qualified_domain
    features = extract_mask_features(labels, (1., 1., 2.))
    content = build_mask_feature_bundle(labels, features, spacing_um=(1., 1., 2.),
                                        object_type="nucleus", fill_holes=False, segmentation_config={})
    with zipfile.ZipFile(io.BytesIO(content)) as archive:
        csv = pd.read_csv(io.BytesIO(archive.read("features.csv")), keep_default_na=False,
                          float_precision="round_trip")
        provenance = json.loads(archive.read("measurement_provenance.json"))
    for name in SURFACE_METADATA_COLUMNS:
        assert csv.iloc[0][name] == features.iloc[0][name], name
    assert csv.volume_um3.tolist() == [2.]
    assert csv.equivalent_diameter_um.tolist() == features.equivalent_diameter_um.tolist()
    policy = provenance["measurement_policy"]
    assert policy["object_type"] == "nucleus"
    assert policy["surface_domain_fields_apply_to"] == ["surface_area", "sphericity", "surface_to_volume_ratio"]
    assert policy["biological_measurement_validity_status"] == "NOT ASSESSED"
