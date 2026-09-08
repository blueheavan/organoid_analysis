from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime
import json
from pathlib import Path
import threading
from typing import Callable
import zipfile

import numpy as np
import tifffile

from organoid_analysis.microscopy_io import Spacing, ZStack, load_zstack, read_axes
from organoid_analysis.segmentation.paths import SEGMENTATION_OUTPUT_DIR, detect_torch_acceleration


@dataclass(frozen=True)
class SegmentationConfig:
    model_type: str = "cpdino-vitb"
    nuclei_diameter: float = 30.0
    cell_diameter: float = 50.0
    anisotropy: float = 2.9
    # Provenance for `anisotropy`: "metadata" (read from the uploaded TIFF via
    # Auto-detect), "default" (no metadata was available, using the fallback
    # above), or "user_override" (a human typed a different value). Frozen
    # into the saved run config so downstream volume/area/distance figures can
    # be traced back to whether they used a real or an assumed voxel size.
    anisotropy_source: str = "default"
    metadata_anisotropy: float | None = None
    xy_spacing_um: float = 0.414
    xy_spacing_source: str = "default"
    metadata_xy_spacing_um: float | None = None
    nuclei_flow_threshold: float = 0.4
    cell_flow_threshold: float = 0.6
    nuclei_cellprob_threshold: float = 0.0
    cell_cellprob_threshold: float = 0.0
    flow3d_smooth: float = 1.0
    batch_size: int = 8
    xy_downsample: float = 1.0


# Working models in this Cellpose v4 build, named and ordered exactly as in the
# official MODELS list (models.MODEL_NAMES). The classic fast U-Net models
# (cyto3/nuclei) no longer exist in v4, so only the foundation models are
# offered. All of them support 3D segmentation via the orthogonal-plane run_3D.
MODEL_OPTIONS = ("cpsam_v2", "cpdino", "cpdino-vitb", "cpsam")

# All model options support 3D; kept for clarity/membership tests.
_3D_CAPABLE = set(MODEL_OPTIONS)

# Cellpose's progress adapter temporarily replaces a process-global function,
# and Streamlit caches one model across sessions. Serialize model inference so
# concurrent sessions cannot use one model or overwrite one another's adapter.
_INFERENCE_LOCK = threading.Lock()


@dataclass(frozen=True)
class SegmentationResult:
    run_directory: Path
    nuclei_mask_path: Path
    cell_mask_path: Path | None
    summary_path: Path
    archive_path: Path
    nuclei_count: int
    cell_count: int


@dataclass(frozen=True)
class RestoredSegmentationRun:
    """A persisted Cellpose run reloaded independently of Streamlit state."""

    result: SegmentationResult
    config: SegmentationConfig
    nuclei_masks: np.ndarray
    cell_masks: np.ndarray | None


def read_stack(path: str | Path) -> np.ndarray:
    stack = tifffile.imread(path)
    if stack.ndim != 3:
        raise ValueError(f"Expected a 3D TIFF stack (Z, Y, X), received shape {stack.shape}.")
    if min(stack.shape) < 2:
        raise ValueError(f"TIFF stack is too small for 3D segmentation: {stack.shape}.")
    return stack


def read_stack_multichannel(path: str | Path) -> ZStack:
    """Axis-aware TIFF load: resolves OME/ImageJ metadata instead of assuming
    a plain 3D array. Returns a :class:`ZStack` whose ``volume`` is ``(Z,Y,X)``
    for single-channel input or ``(C,Z,Y,X)`` when ``is_multichannel`` is set,
    letting the caller pick a channel before segmentation.

    Falls back to the legacy plain-3D read when the TIFF's page sequence
    carries no Z/channel metadata at all (``load_zstack`` refuses to guess
    axis order for those rather than misreading them) -- most raw single-
    channel Z-stacks saved without ImageJ/OME tags fall in this case, and
    should keep working exactly as they did before this function existed.
    """
    try:
        zstack = load_zstack(path)
    except ValueError:
        # Legacy plain multipage grayscale TIFFs are commonly reported as QYX
        # or IYX. Only those explicitly metadata-free layouts may use the old
        # reader; never reinterpret rejected CYX/RGB/ambiguous metadata as ZYX.
        with tifffile.TiffFile(path) as handle:
            if (read_axes(path).upper() not in {"QYX", "IYX"}
                    or handle.ome_metadata is not None
                    or getattr(handle, "imagej_metadata", None) is not None):
                raise
        stack = read_stack(path)
        return ZStack(
            volume=stack, spacing=Spacing(None, None, None), dtype=stack.dtype,
            axes="ZYX", channels=1, frames=1, is_multichannel=False,
            meta={"source": str(path)},
        )
    spatial_shape = zstack.volume.shape[-3:]
    if min(spatial_shape) < 2:
        raise ValueError(f"TIFF stack is too small for 3D segmentation: {spatial_shape}.")
    return zstack


def validate_stacks(nuclei: np.ndarray, cells: np.ndarray) -> None:
    if nuclei.ndim != 3 or cells.ndim != 3:
        raise ValueError("Nuclei and cell inputs must both be 3D TIFF stacks in (Z, Y, X) order.")
    if nuclei.shape != cells.shape:
        raise ValueError(
            f"Nuclei and cell stacks must have identical shapes; got {nuclei.shape} and {cells.shape}."
        )


def normalize_preview(image: np.ndarray) -> np.ndarray:
    low, high = np.percentile(image, (1, 99))
    if high <= low:
        return np.zeros(image.shape, dtype=np.float32)
    return np.clip((image.astype(np.float32) - low) / (high - low), 0, 1)


def get_accelerator() -> str:
    import torch

    return detect_torch_acceleration(torch)[0]


def create_model(model_type: str = ""):
    import torch
    from cellpose import models

    if not model_type:
        model_type = SegmentationConfig().model_type
    if model_type not in MODEL_OPTIONS:
        raise ValueError(
            f"Unsupported model_type {model_type!r}. Choose one of {MODEL_OPTIONS}."
        )
    accelerator, _ = detect_torch_acceleration(torch)
    return models.CellposeModel(
        device=torch.device(accelerator),
        pretrained_model=model_type,
        use_bfloat16=accelerator != "mps",
    )


def segment_stacks(
    model,
    nuclei: np.ndarray,
    cells: np.ndarray | None,
    config: SegmentationConfig,
    on_progress: Callable[[float], None] | None = None,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Run 3D segmentation, optionally segmenting cell masks when ``cells`` is given.

    ``on_progress`` receives a fraction in [0, 1]. Progress is reported at the
    stage boundaries (nuclei then cell passes) and additionally *inside* each 3D
    network pass: Cellpose's ``run_3D`` processes the volume along its three
    orthogonal planes sequentially, and each completed plane is translated into a
    ``on_progress`` fraction, so the UI percentage climbs continuously instead of
    staying pinned until the pass finishes.
    """
    if cells is not None:
        validate_stacks(nuclei, cells)
    if on_progress:
        on_progress(0.0)

    downsample = config.xy_downsample
    if not 0 < downsample <= 1.0:
        raise ValueError("xy_downsample must be in (0, 1], with 1.0 meaning full resolution.")
    work_nuclei = nuclei
    work_cells = cells
    if downsample < 1.0:
        from scipy import ndimage

        factors = (1.0, downsample, downsample)
        work_nuclei = ndimage.zoom(nuclei, factors, order=1)
        if cells is not None:
            work_cells = ndimage.zoom(cells, factors, order=1)

    common = {
        "z_axis": 0,
        "do_3D": config.model_type in _3D_CAPABLE,
        "anisotropy": config.anisotropy * downsample,
        "flow3D_smooth": config.flow3d_smooth,
        "batch_size": config.batch_size,
    }
    if on_progress:
        on_progress(0.1)

    # The 3D network pass is monolithic, so funnel its per-plane events into the
    # outer progress callback. Nuclei-only runs map planes to 10% -> 95% (the
    # remaining 5% guards mask post-processing); with a cell pass, nuclei planes
    # fill 10% -> 55% and the cell planes fill 55% -> 100%.
    nuclei_end = 0.55 if cells is not None else 0.95
    cell_masks = None
    with _INFERENCE_LOCK:
        restore = _wrap_run_3d(on_progress, 0.10, nuclei_end)
        try:
            nuclei_masks = _validated_mask(
                model.eval(
                    work_nuclei,
                    diameter=config.nuclei_diameter,
                    flow_threshold=config.nuclei_flow_threshold,
                    cellprob_threshold=config.nuclei_cellprob_threshold,
                    **common,
                )[0],
                work_nuclei.shape,
                "nuclei",
            )
        finally:
            restore()
        if cells is not None:
            if on_progress:
                on_progress(0.55)
            cell_input = np.stack((work_cells, work_nuclei), axis=-1)
            restore = _wrap_run_3d(on_progress, 0.55, 1.0)
            try:
                cell_masks = _validated_mask(
                    model.eval(
                        cell_input,
                        channel_axis=-1,
                        diameter=config.cell_diameter,
                        flow_threshold=config.cell_flow_threshold,
                        cellprob_threshold=config.cell_cellprob_threshold,
                        **common,
                    )[0],
                    work_cells.shape,
                    "cell",
                )
            finally:
                restore()
    if cells is None:
        if downsample < 1.0:
            nuclei_masks = _upsample_masks(nuclei_masks, nuclei.shape)
        if on_progress:
            on_progress(1.0)
        return nuclei_masks, None
    if downsample < 1.0:
        nuclei_masks = _upsample_masks(nuclei_masks, nuclei.shape)
        cell_masks = _upsample_masks(cell_masks, cells.shape)
    if on_progress:
        on_progress(1.0)
    return nuclei_masks, cell_masks


def _validated_mask(mask: np.ndarray, expected_shape: tuple[int, ...], name: str) -> np.ndarray:
    """Validate model output before it reaches persistence or measurements."""
    result = np.asarray(mask)
    if result.shape != expected_shape:
        raise RuntimeError(f"Cellpose returned {name} mask shape {result.shape}; expected {expected_shape}")
    if not np.issubdtype(result.dtype, np.integer) or np.any(result < 0):
        raise RuntimeError(f"Cellpose returned invalid {name} instance labels")
    if result.size and int(result.max()) > np.iinfo(np.uint32).max:
        raise RuntimeError(f"Cellpose returned {name} labels outside the uint32 range")
    return result.astype(np.uint32, copy=False)


def _instance_count(mask: np.ndarray) -> int:
    """Count occurring positive labels instead of assuming IDs are compact."""
    values = np.unique(mask)
    return int(np.count_nonzero(values))


def _upsample_masks(masks: np.ndarray, full_shape: tuple[int, ...]) -> np.ndarray:
    from scipy import ndimage

    factors = tuple(full / orig for full, orig in zip(full_shape, masks.shape))
    return ndimage.zoom(masks, factors, order=0, mode="nearest").astype(np.uint32)


def _wrap_run_3d(
    on_progress: Callable[[float], None] | None, start: float, end: float
):
    """Route Cellpose's ``run_3D`` plane reports into ``on_progress``.

    Cellpose's 3D inference walks the volume along its three orthogonal planes
    (YX, ZY, ZX) and calls ``progress.setValue(25 + 15 * p)`` (i.e. 25/40/55)
    after each plane. ``model.eval`` itself exposes no progress hook for the
    network forward pass, so this wrapper temporarily rebinds ``run_3D`` (the
    name ``CellposeModel._run_net`` resolves) to translate those three events
    into fractions within ``[start, end]`` of the overall stage.

    Returns a no-op restore guard; call `.restore()` once the eval is done.
    """
    from cellpose import models as _models

    original = _models.run_3D

    class _PlaneProgress:
        __slots__ = ("start", "end")

        def __init__(self):
            self.start = start
            self.end = end

        def setValue(self, value: int) -> None:
            if on_progress is None:
                return
            local = max(0.0, min(1.0, (value - 25) / 30.0))
            on_progress(self.start + (self.end - self.start) * local)

    def wrapped(net, imgs, *args, **kwargs):
        kwargs["progress"] = _PlaneProgress()
        return original(net, imgs, *args, **kwargs)

    _models.run_3D = wrapped

    def restore() -> None:
        _models.run_3D = original

    return restore


def save_result(
    nuclei_masks: np.ndarray,
    cell_masks: np.ndarray | None,
    config: SegmentationConfig,
    output_directory: Path = SEGMENTATION_OUTPUT_DIR,
) -> SegmentationResult:
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S-%f")
    run_directory = output_directory / f"run-{timestamp}"
    run_directory.mkdir(parents=True, exist_ok=False)

    nuclei_mask_path = run_directory / "nuclei_masks_3d.tif"
    summary_path = run_directory / "summary.json"
    archive_path = run_directory / "segmentation_results.zip"
    nuclei_masks = _validated_mask(nuclei_masks, nuclei_masks.shape, "nuclei")
    paths = [nuclei_mask_path]
    tifffile.imwrite(nuclei_mask_path, nuclei_masks, compression="zlib")
    cell_count = 0
    if cell_masks is not None:
        cell_masks = _validated_mask(cell_masks, nuclei_masks.shape, "cell")
        cell_mask_path = run_directory / "cell_masks_3d.tif"
        tifffile.imwrite(cell_mask_path, cell_masks, compression="zlib")
        cell_count = _instance_count(cell_masks)
        paths.append(cell_mask_path)
    else:
        cell_mask_path = None

    summary = {
        "cellpose_model": config.model_type,
        "config": asdict(config),
        "nuclei_count": _instance_count(nuclei_masks),
        "cell_count": cell_count,
        "shape": list(nuclei_masks.shape),
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    paths.append(summary_path)
    with zipfile.ZipFile(archive_path, "w") as archive:
        for path in paths:
            compression = zipfile.ZIP_DEFLATED if path == summary_path else zipfile.ZIP_STORED
            archive.write(path, path.name, compress_type=compression)

    return SegmentationResult(
        run_directory=run_directory,
        nuclei_mask_path=nuclei_mask_path,
        cell_mask_path=cell_mask_path,
        summary_path=summary_path,
        archive_path=archive_path,
        nuclei_count=summary["nuclei_count"],
        cell_count=cell_count,
    )


def list_saved_results(output_directory: Path = SEGMENTATION_OUTPUT_DIR) -> list[SegmentationResult]:
    """List complete saved runs newest first without decoding their masks."""
    if not output_directory.is_dir():
        return []
    results = []
    for run_directory in output_directory.glob("run-*"):
        summary_path = run_directory / "summary.json"
        nuclei_path = run_directory / "nuclei_masks_3d.tif"
        if not summary_path.is_file() or not nuclei_path.is_file():
            continue
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            cell_path = run_directory / "cell_masks_3d.tif"
            results.append(SegmentationResult(
                run_directory=run_directory, nuclei_mask_path=nuclei_path,
                cell_mask_path=cell_path if cell_path.is_file() else None,
                summary_path=summary_path, archive_path=run_directory / "segmentation_results.zip",
                nuclei_count=int(summary["nuclei_count"]),
                cell_count=int(summary.get("cell_count", 0)) if cell_path.is_file() else 0,
            ))
        except (OSError, ValueError, KeyError, TypeError):
            continue
    return sorted(results, key=lambda item: item.run_directory.name, reverse=True)


def restore_saved_result(run_directory: str | Path) -> RestoredSegmentationRun:
    """Reload masks and saved physical configuration from a completed run."""
    directory = Path(run_directory).resolve()
    available = {result.run_directory.resolve(): result for result in list_saved_results(directory.parent)}
    if directory not in available:
        raise FileNotFoundError(f"Not a complete saved segmentation run: {directory}")
    result = available[directory]
    summary = json.loads(result.summary_path.read_text(encoding="utf-8"))
    try:
        config = SegmentationConfig(**summary["config"])
    except (KeyError, TypeError) as error:
        raise ValueError(f"Saved run has invalid segmentation configuration: {directory}") from error
    nuclei = _validated_mask(tifffile.imread(result.nuclei_mask_path), tuple(summary["shape"]), "nuclei")
    cells = None
    if result.cell_mask_path is not None:
        cells = _validated_mask(tifffile.imread(result.cell_mask_path), nuclei.shape, "cell")
    return RestoredSegmentationRun(result=result, config=config, nuclei_masks=nuclei, cell_masks=cells)
