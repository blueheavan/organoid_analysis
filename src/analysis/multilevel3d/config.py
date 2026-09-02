"""Configuration for multilevel 3D measurement and QC."""
from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class Multilevel3DConfig:
    """Centralized physical-QC thresholds for registered label volumes.

    Fractions use the child's full voxel volume as their denominator.  Radial
    thresholds deliberately use *equivalent-radius* normalization, rather
    than claiming an object-specific local radius.
    """

    minimum_voxels: int = 5
    low_parent_overlap_fraction: float = 0.5
    mad_z_threshold: float = 3.5
    core_max_normalized_radial_position: float = 0.5
    peripheral_min_normalized_radial_position: float = 0.8

    def as_dict(self) -> dict[str, float | int]:
        return asdict(self)

    def validate(self) -> None:
        if self.minimum_voxels < 1:
            raise ValueError("minimum_voxels must be at least one")
        if not 0 <= self.low_parent_overlap_fraction <= 1:
            raise ValueError("low_parent_overlap_fraction must be in [0, 1]")
        if self.mad_z_threshold <= 0:
            raise ValueError("mad_z_threshold must be positive")
        if not 0 <= self.core_max_normalized_radial_position <= 1:
            raise ValueError("core_max_normalized_radial_position must be in [0, 1]")
        if not 0 <= self.peripheral_min_normalized_radial_position <= 1:
            raise ValueError("peripheral_min_normalized_radial_position must be in [0, 1]")
        if self.core_max_normalized_radial_position > self.peripheral_min_normalized_radial_position:
            raise ValueError("core threshold cannot exceed peripheral threshold")
