"""Headless 3D segmentation for the organoid pipeline.

Package layout:
* ``cellpose``    — Cellpose 3D segmentation service (nuclei + cells).
* ``auto_config`` — automatic inference of segmentation parameters.
* ``io``          — Z-stack TIFF IO and physical-metadata handling.
* ``paths``       — project path constants and torch acceleration detection.
"""

from . import io  # noqa: F401
from .cellpose import (  # noqa: F401
    MODEL_OPTIONS,
    SegmentationConfig,
    SegmentationResult,
    RestoredSegmentationRun,
    create_model,
    get_accelerator,
    normalize_preview,
    read_stack,
    read_stack_multichannel,
    list_saved_results,
    restore_saved_result,
    save_result,
    segment_stacks,
    validate_stacks,
)
from .auto_config import (  # noqa: F401
    auto_anisotropy,
    estimate_diameter_from_stack,
    suggest_config,
)

__all__ = [
    "MODEL_OPTIONS",
    "SegmentationConfig",
    "SegmentationResult",
    "RestoredSegmentationRun",
    "auto_anisotropy",
    "create_model",
    "estimate_diameter_from_stack",
    "get_accelerator",
    "normalize_preview",
    "read_stack",
    "read_stack_multichannel",
    "list_saved_results",
    "restore_saved_result",
    "save_result",
    "segment_stacks",
    "suggest_config",
    "validate_stacks",
]
