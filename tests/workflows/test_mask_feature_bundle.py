import hashlib
import io
import json
import zipfile
from dataclasses import asdict

import numpy as np
import pandas as pd
import pytest

from organoid_analysis.quantification.mask_features import extract_mask_features
from organoid_analysis.result_export.mask_feature_bundle import build_mask_feature_bundle
from organoid_analysis.segmentation.cellpose_inference import SegmentationConfig


@pytest.mark.parametrize("object_type,fill_holes", [("nucleus", False), ("cell", True)])
def test_feature_bundle_preserves_values_ids_calibration_and_lineage(object_type, fill_holes):
    mask = np.zeros((7, 7, 7), np.uint32)
    mask[1:6, 1:6, 1:6] = 2**24 + 19
    mask[2, 2, 2] = 0
    features = extract_mask_features(mask, (.01, .01, .03), fill_holes=fill_holes)
    config = SegmentationConfig(xy_spacing_um=.01, anisotropy=3., xy_spacing_source="metadata", anisotropy_source="user_override")
    kwargs = dict(spacing_um=(.01, .01, .03), object_type=object_type, fill_holes=fill_holes,
                  segmentation_config=asdict(config), segmentation_provenance={"nuclei_input": {"sha256": "source-file-digest"}})
    content = build_mask_feature_bundle(mask, features, **kwargs)
    assert content == build_mask_feature_bundle(mask, features, **kwargs)
    with zipfile.ZipFile(io.BytesIO(content)) as archive:
        csv_bytes = archive.read("features.csv")
        table = pd.read_csv(io.BytesIO(csv_bytes), float_precision="round_trip")
        meta = json.loads(archive.read("measurement_provenance.json"))
    assert table.Label.tolist() == [2**24 + 19]
    assert table.volume_um3[0] == pytest.approx((125 if fill_holes else 124) * .000003, rel=1e-12)
    assert meta["files"]["features.csv"]["sha256"] == hashlib.sha256(csv_bytes).hexdigest()
    assert meta["object_type"] == object_type
    assert meta["measurement_basis"] == ("filled_envelope" if fill_holes else "raw_label")
    assert meta["calibration_status"] == "provided_not_independently_verified"
    assert meta["spacing_xyz_um"] == [.01, .01, .03]
    assert meta["segmentation_provenance"]["nuclei_input"]["sha256"] == "source-file-digest"
    # The exported surface values are tied to a versioned estimator definition.
    assert meta["surface_area_method"]["method_version"] == "crofton_minimax_sym_v3"
    # The exported identity carries the domain the accuracy claim is limited to.
    assert "rho_in >= 10.0" in meta["surface_area_method"]["domain"]
    assert "crease" in meta["surface_area_method"]["not_qualified"]
    header = {"shape_zyx": [7, 7, 7], "dtype": mask.dtype.str, "order": "C"}
    assert meta["selected_mask"]["array_sha256"] == hashlib.sha256(json.dumps(header, sort_keys=True).encode("ascii") + mask.tobytes()).hexdigest()


def test_empty_bundle_and_unknown_calibration_are_explicit():
    mask = np.zeros((3, 3, 3), np.uint8)
    content = build_mask_feature_bundle(mask, extract_mask_features(mask), spacing_um=(1, 1, 1),
                                        object_type="nucleus", fill_holes=False, segmentation_config={})
    with zipfile.ZipFile(io.BytesIO(content)) as archive:
        meta = json.loads(archive.read("measurement_provenance.json"))
        assert pd.read_csv(io.BytesIO(archive.read("features.csv"))).empty
    assert meta["n_objects"] == 0 and meta["n_review_objects"] == 0
    assert meta["calibration_status"] == "assumed_or_unknown"
    assert meta["segmentation_provenance"] is None


def test_mismatched_mask_or_geometry_cannot_be_exported():
    mask = np.ones((3, 3, 3), np.uint8)
    features = extract_mask_features(mask)
    kwargs = dict(spacing_um=(1, 1, 1), object_type="nucleus", fill_holes=False, segmentation_config={})
    with pytest.raises(ValueError, match="instance IDs"):
        build_mask_feature_bundle(mask * 2, features, **kwargs)
    with pytest.raises(ValueError, match="measurement basis"):
        build_mask_feature_bundle(mask, features, **{**kwargs, "fill_holes": True})
    with pytest.raises(ValueError, match="voxel spacing"):
        build_mask_feature_bundle(mask, features, **{**kwargs, "spacing_um": (2, 2, 2)})
