"""Z-stack TIFF IO and metadata handling."""

from .metadata import (
    Spacing,
    parse_spacing_ome,
    parse_spacing_imagej,
    resolve_spacing,
    parse_axes,
    reject_rgb,
)
from .zstack_reader import (
    ZStack,
    load_zstack,
    read_axes,
)

__all__ = [
    "Spacing",
    "parse_spacing_ome",
    "parse_spacing_imagej",
    "resolve_spacing",
    "parse_axes",
    "reject_rgb",
    "ZStack",
    "load_zstack",
    "read_axes",
]
