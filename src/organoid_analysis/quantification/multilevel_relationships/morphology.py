"""True physical 3D morphology for individual label instances."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

from organoid_analysis.quantification.features import (
    DERIVED_GEOMETRY_COLUMNS,
    SURFACE_METADATA_COLUMNS,
    geometry,
)
from organoid_analysis.quantification.labels import (
    bbox_touches_volume_boundary,
    compact_instance_labels,
)

MORPHOLOGY_COLUMNS = [
    "voxel_count", "volume_um3", "surface_area_um2", "sphericity", "centroid_z_um", "centroid_y_um",
    "centroid_x_um", "major_axis_um", "intermediate_axis_um", "minor_axis_um",
    "equivalent_diameter_um", "axis_ratio_minor_to_major", "measurement_basis",
    *DERIVED_GEOMETRY_COLUMNS, *SURFACE_METADATA_COLUMNS, "touches_image_border", "fragmented_object",
    "connected_component_count",
]


def measure_instances(labels: np.ndarray, spacing_zyx_um: tuple[float, float, float], *, object_type: str) -> pd.DataFrame:
    """Measure every present sparse label in physical units, never from a MIP."""
    if not np.any(labels):
        return pd.DataFrame(columns=[f"{object_type}_id", *MORPHOLOGY_COLUMNS])
    compact_labels, source_id_by_compact_id = compact_instance_labels(labels)

    # ``find_objects`` indexes boxes by ID; compaction preserves sparse source
    # IDs in the exported table without making that index depend on max(label).
    rows: list[dict] = []
    for compact_id, bbox in enumerate(ndi.find_objects(compact_labels), start=1):
        if bbox is None:
            continue
        object_id = source_id_by_compact_id[compact_id]
        mask = compact_labels[bbox] == compact_id
        # All reported measurements use the raw label, not a filled envelope.
        measured, _ = geometry(mask, spacing_zyx_um, tuple(part.start for part in bbox), fill_holes=False)
        axes = [measured["principal_axis_major_um"], measured["principal_axis_intermediate_um"], measured["principal_axis_minor_um"]]
        components = int(ndi.label(mask, structure=ndi.generate_binary_structure(3, 1))[1])
        voxel_count = int(measured["segmented_voxels"])
        volume = float(measured["volume_um3"])
        rows.append({
            f"{object_type}_id": int(object_id),
            "voxel_count": voxel_count,
            "volume_um3": volume,
            "surface_area_um2": measured["surface_area_um2"],
            "sphericity": measured["sphericity"],
            "centroid_z_um": measured["centroid_z_um"],
            "centroid_y_um": measured["centroid_y_um"],
            "centroid_x_um": measured["centroid_x_um"],
            "major_axis_um": axes[0],
            "intermediate_axis_um": axes[1],
            "minor_axis_um": axes[2],
            **{name: measured[name] for name in (
                "equivalent_diameter_um", "axis_ratio_minor_to_major", "measurement_basis",
                *DERIVED_GEOMETRY_COLUMNS, *SURFACE_METADATA_COLUMNS,
            )},
            "touches_image_border": bbox_touches_volume_boundary(bbox, labels.shape),
            "fragmented_object": components > 1,
            "connected_component_count": components,
        })
    return pd.DataFrame(rows, columns=[f"{object_type}_id", *MORPHOLOGY_COLUMNS])
