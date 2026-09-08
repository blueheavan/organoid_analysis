"""Parquet feature database and JSON-summary export."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


def _json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Not JSON serializable: {type(value).__name__}")


def export_results(output_directory: str | Path, *, organoids: pd.DataFrame, cells: pd.DataFrame,
                   nuclei: pd.DataFrame, edges: pd.DataFrame, qc_flags: pd.DataFrame, summary: dict,
                   label_masks: dict[str, object] | None = None,
                   spacing_zyx_um: tuple[float, float, float] | None = None) -> dict[str, str]:
    """Write the stable Result-like hierarchy; output must be new or empty."""
    output = Path(output_directory).resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {output}")
    (output / "features").mkdir(parents=True, exist_ok=True)
    (output / "qc").mkdir(parents=True, exist_ok=True)
    (output / "summary").mkdir(parents=True, exist_ok=True)
    paths = {
        "organoid_features": output / "features" / "organoid_features.parquet",
        "cell_features": output / "features" / "cell_features.parquet",
        "nucleus_features": output / "features" / "nucleus_features.parquet",
        "cell_topology_edges": output / "features" / "cell_topology_edges.parquet",
        "qc_flags": output / "qc" / "qc_flags.parquet",
        "analysis_summary": output / "summary" / "analysis_summary.json",
    }
    organoids.to_parquet(paths["organoid_features"], index=False)
    cells.to_parquet(paths["cell_features"], index=False)
    nuclei.to_parquet(paths["nucleus_features"], index=False)
    edges.to_parquet(paths["cell_topology_edges"], index=False)
    qc_flags.to_parquet(paths["qc_flags"], index=False)
    paths["analysis_summary"].write_text(json.dumps(summary, indent=2, default=_json_default), encoding="utf-8")
    if label_masks is not None:
        if spacing_zyx_um is None:
            raise ValueError("spacing_zyx_um is required when exporting label masks")
        from organoid_analysis.microscopy_io.tiff_contract import write_labels
        for level, labels in label_masks.items():
            if level not in {"organoid", "cell", "nucleus"}:
                raise ValueError(f"Unknown label-mask level: {level}")
            path = output / "masks" / f"{level}_labels.ome.tif"
            path.parent.mkdir(parents=True, exist_ok=True)
            write_labels(path, labels, spacing_zyx_um)
            paths[f"{level}_labels"] = path
    return {key: str(value) for key, value in paths.items()}
