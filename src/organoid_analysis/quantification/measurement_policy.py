"""Object reporting policies; all numerical measurements belong to geometry()."""
from __future__ import annotations

MEASUREMENT_POLICY_VERSION = "1.0"
MEASUREMENT_INTERPRETATION = (
    "Volume is voxel count times voxel volume; equivalent diameter, moment axes, "
    "axis ratios and elongation are primary mask measurements. Crofton surface area, "
    "sphericity and surface-to-volume ratio are conditional measurements. "
    "surface_in_qualified_domain checks numerical resolution and anisotropy only; "
    "sphericity and surface-to-volume ratio inherit that surface gate. Smoothness "
    "is NOT ASSESSED from the mask. The gate does not validate segmentation boundaries, "
    "biological accuracy or volume accuracy. All objects remain reported outside the "
    "surface domain. Summaries include those objects and do not establish validation. "
    "Weights solved or cached outside the packaged tables are not the frozen "
    "evidence-bearing vectors; their provenance does not expand the qualified scope."
)


def measurement_policy(object_type: str, measurement_basis: str) -> dict[str, object]:
    """Describe a labelled object's estimand without inferring biological validity."""
    if object_type not in {"organoid", "cell", "nucleus"}:
        raise ValueError("object_type must be organoid, cell or nucleus")
    if measurement_basis not in {"raw_label", "filled_envelope"}:
        raise ValueError("measurement_basis must be raw_label or filled_envelope")
    return {
        "policy_version": MEASUREMENT_POLICY_VERSION,
        "object_type": object_type,
        "measurement_basis": measurement_basis,
        "volume_semantics": f"voxel-count volume of the {object_type} {measurement_basis}",
        "primary_metrics": ["volume", "equivalent_diameter", "principal_axes", "axis_ratios", "elongation"],
        "conditional_metrics": ["surface_area", "sphericity", "surface_to_volume_ratio"],
        "surface_domain_fields_apply_to": ["surface_area", "sphericity", "surface_to_volume_ratio"],
        "surface_scope_conditions_status": "NOT ASSESSED",
        "segmentation_validation_status": "NOT ASSESSED",
        "biological_measurement_validity_status": "NOT ASSESSED",
        "volume_accuracy_status": "NOT ASSESSED",
        "interpretation": MEASUREMENT_INTERPRETATION,
    }
