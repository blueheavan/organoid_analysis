"""Headless 3D segmentation for the organoid pipeline.

Package layout:
* ``cellpose_inference``  — Cellpose 3D segmentation service (nuclei + cells).
* ``watershed_instances`` — classical (non-deep-learning) Otsu + watershed segmentation.
* ``parameter_estimation`` — automatic inference of segmentation parameters.
* ``paths``               — project path constants and torch acceleration detection.

Z-stack TIFF IO and physical-metadata handling now lives in
:mod:`organoid_analysis.microscopy_io`, alongside it.
"""

from .cellpose_inference import (  # noqa: F401
    MODEL_OPTIONS,
    RestoredSegmentationRun,
    SegmentationConfig,
    SegmentationResult,
    create_model,
    get_accelerator,
    list_saved_results,
    normalize_preview,
    read_stack,
    read_stack_multichannel,
    restore_saved_result,
    save_result,
    segment_stacks,
    validate_stacks,
)
from .parameter_estimation import (  # noqa: F401
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
