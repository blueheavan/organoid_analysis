"""Z-stack TIFF IO and metadata handling."""

from .metadata import (
    GRID_CONFLICT,
    GRID_CONSISTENT,
    GRID_PARTIAL,
    Spacing,
    compare_registered_grid,
    parse_axes,
    parse_spacing_imagej,
    parse_spacing_ome,
    reject_rgb,
    resolve_spacing,
)
from .voxel_spacing import (
    isotropic_xy_size_um,
    resolve_spacing_source,
    validate_voxel_spacing_xyz,
)
from .zstack_reader import (
    ZStack,
    load_zstack,
    read_axes,
)

__all__ = [
    "Spacing",
    "GRID_CONSISTENT",
    "GRID_CONFLICT",
    "GRID_PARTIAL",
    "compare_registered_grid",
    "parse_spacing_ome",
    "parse_spacing_imagej",
    "resolve_spacing",
    "parse_axes",
    "reject_rgb",
    "ZStack",
    "load_zstack",
    "read_axes",
    "isotropic_xy_size_um",
    "resolve_spacing_source",
    "validate_voxel_spacing_xyz",
]
