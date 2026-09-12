"""Audit counterexamples; does not modify production source or claim validation."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import tempfile
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from sklearn.metrics import silhouette_score

from organoid_analysis.config import load_config
from organoid_analysis.microscopy_io.metadata import Spacing, compare_registered_grid
from organoid_analysis.microscopy_io.tiff_contract import load_sample
from organoid_analysis.quantification.features import geometry
from organoid_analysis.quantification.mask_features import extract_mask_features
from organoid_analysis.statistics.exploration import (
    determine_optimal_clusters,
    perform_kmeans_clustering,
    silhouette_for,
)
from organoid_analysis.workflows.multilevel_measurement_workflow import analyze_multilevel_3d

HERE = Path(__file__).resolve().parent
RESULTS: dict = {}


def record(name, fn):
    try:
        RESULTS[name] = fn()
    except Exception:
        RESULTS[name] = {"unexpected_exception": traceback.format_exc()}
    print(name, json.dumps(RESULTS[name], allow_nan=False), flush=True)


def spheres():
    records = []
    for r in (10.0, 20.0, 40.0):
        for spacing in ((1.0, 1.0, 1.0), (2.0, 1.0, 1.0)):
            half = np.ceil(r / np.asarray(spacing)).astype(int) + 2
            z, y, x = np.ogrid[tuple(slice(-h, h + 1) for h in half)]
            mask = (z * spacing[0]) ** 2 + y ** 2 + x ** 2 <= r ** 2
            g, _ = geometry(mask, spacing)
            volume_truth, area_truth = 4 * np.pi * r ** 3 / 3, 4 * np.pi * r ** 2
            scaled, _ = geometry(mask, tuple(3 * v for v in spacing))
            records.append({
                "radius_um": r, "spacing_zyx_um": spacing,
                "volume_relative_error": g["volume_um3"] / volume_truth - 1,
                "surface_relative_error": g["surface_area_um2"] / area_truth - 1,
                "sphericity": g["sphericity"],
                "volume_criterion_pass": abs(g["volume_um3"] / volume_truth - 1) < .01,
                "surface_criterion_pass": abs(g["surface_area_um2"] / area_truth - 1) < .05,
                "scale_volume_ratio": scaled["volume_um3"] / g["volume_um3"],
                "scale_area_ratio": scaled["surface_area_um2"] / g["surface_area_um2"],
            })
    return records


def partial_spacing():
    with tempfile.TemporaryDirectory(prefix="organoid-science-audit-") as tmp:
        image = Path(tmp) / "partial.ome.tif"
        tifffile.imwrite(image, np.ones((7, 12, 12), np.uint16), ome=True,
                         photometric="minisblack", metadata={"axes": "ZYX",
                         "PhysicalSizeX": 1., "PhysicalSizeY": 1.,
                         "PhysicalSizeXUnit": "µm", "PhysicalSizeYUnit": "µm"})
        row = {"image_path": str(image), "spacing_z_um": "2",
               "spacing_y_um": "2", "spacing_x_um": "2"}
        digest = hashlib.sha256(image.read_bytes()).hexdigest()
        try:
            sample = load_sample(row, load_config())
            primary = {"rejected": False, "accepted_spacing_zyx_um": sample.spacing,
                       "metadata": sample.metadata}
        except ValueError as exc:
            primary = {"rejected": True, "reason": str(exc)}
        return {"known_ome_xy_um": [1., 1.], "manifest_xy_um": [2., 2.],
                "fixture_sha256": digest, "primary": primary,
                "partial_grid_comparison": compare_registered_grid(
                    Spacing(1., 1., None), Spacing(2., 2., 2.))}


def clustering_scales():
    rng = np.random.default_rng(20260911)
    centers = np.repeat([-3., 0., 3.], 60)
    df = pd.DataFrame({"Value of signal": centers + rng.normal(0, .25, len(centers)),
                       "Value of other": rng.normal(0, 1, len(centers))})
    columns = list(df)
    records = []
    for scale in (1., 1000.):
        data = df.copy()
        data[columns[1]] *= scale
        recommended, _, raw_scores = determine_optimal_clusters(data.to_numpy(), max_k=8)
        clustered, standardized, _, _ = perform_kmeans_clustering(data, columns, 3)
        std_k, _, std_scores = determine_optimal_clusters(standardized, max_k=8)
        records.append({"second_feature_unit_scale": scale,
                        "recommended_k_ui_raw": recommended,
                        "recommended_k_standardized": std_k,
                        "ui_raw_candidate_scores": raw_scores,
                        "standardized_candidate_scores": std_scores,
                        "displayed_silhouette": silhouette_for(clustered, columns),
                        "silhouette_on_fitted_features": float(silhouette_score(standardized, clustered.Cluster)),
                        "labels": clustered.Cluster.tolist()})
    return {"rows": records, "fitted_labels_identical": records[0]["labels"] == records[1]["labels"]}


def parent_identity():
    shape = (12, 12, 12)
    organoids = np.zeros(shape, np.int64)
    cells = np.zeros(shape, np.int64)
    nuclei = np.zeros(shape, np.int64)
    parent = 2 ** 53 + 1
    organoids[2:8, 2:8, 2:8] = parent
    cells[3:7, 3:7, 3:7] = 11
    nuclei[4:6, 4:6, 4:6] = 21
    nuclei[9:11, 9:11, 9:11] = 22
    result = analyze_multilevel_3d(organoids, cells, nuclei, (1., 1., 1.))
    return {"expected_parent_id": parent,
            "organoids": result.organoid_features[["organoid_id", "cell_count", "nucleus_count"]].to_dict("records"),
            "cells": result.cell_features[["cell_id", "organoid_id"]].to_dict("records"),
            "nuclei": result.nucleus_features[["nucleus_id", "cell_id", "organoid_id", "direct_organoid_id", "qc_flags"]].to_dict("records")}


def geometric_qc():
    labels = np.zeros((7, 7, 7), np.int64)
    labels[2:4, 2:4, 2:4] = 1
    result = analyze_multilevel_3d(labels, np.zeros_like(labels), np.zeros_like(labels), (1., 1., 1.))
    web = extract_mask_features(labels)
    return {"multilevel": result.organoid_features[["voxel_count", "sphericity", "qc_status", "qc_flags"]].to_dict("records"),
            "web_mask": web[["voxel_count", "sphericity", "qc_status", "qc_flags"]].to_dict("records")}


RESULTS["environment"] = {"python": platform.python_version(), "platform": platform.platform(),
    "packages": {name: importlib.metadata.version(name) for name in
                 ("organoid-analysis", "numpy", "scipy", "pandas", "scikit-image", "scikit-learn", "statsmodels", "tifffile", "cellpose")}}
record("C1_geometry", spheres)
record("C2_partial_metadata", partial_spacing)
record("C3_clustering", clustering_scales)
record("C4_parent_identity", parent_identity)
record("C5_qc", geometric_qc)
(HERE / "probe.json").write_text(json.dumps(RESULTS, indent=2, allow_nan=False) + "\n")
