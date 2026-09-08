"""TIFF axis, spacing, and experimental-design checks. Internal order: Z,Y,X."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import hashlib
import re
import subprocess
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import tifffile

# Tolerance for treating two physical-spacing/position values (µm) as the same
# acquisition metadata rather than a conflict to reject. Heuristic engineering
# judgment ("close enough to be metadata rounding, not a real mismatch"), not
# independently calibrated; see docs/PARAMETERS.md. Reused by cli.py for the
# same spacing-agreement checks in the `cells`/`analyze-3d` CLI routes.
SPACING_RTOL = 0.01
SPACING_ATOL_UM = 1e-5

REQUIRED = ["sample_id", "condition", "biological_replicate", "image_path"]
PATH_FIELDS = ["image_path", "structure_path", "calcein_path", "pi_path", "probability_path", "labels_path", "truth_labels_path"]
META_FIELDS = ["sample_id", "condition", "biological_replicate", "unit_id", "batch_id", "control"]


def read_manifest(path: str | Path) -> pd.DataFrame:
    path = Path(path).resolve()
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = sorted(set(REQUIRED) - set(df.columns))
    if missing or df.empty:
        raise ValueError(f"Manifest must have rows and columns {REQUIRED}; missing={missing}")
    df = df.map(lambda value: value.strip() if isinstance(value, str) else value)
    for field in REQUIRED:
        if (df[field] == "").any():
            raise ValueError(f"Manifest column {field} contains empty values")
    if df.sample_id.duplicated().any():
        raise ValueError("sample_id values must be unique")
    if not df.sample_id.map(lambda x: bool(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,80}", x))).all():
        raise ValueError("sample_id must contain only letters, digits, underscore, or hyphen (1-81 characters)")
    for field, default in [("unit_id", None), ("batch_id", "batch1"), ("control", "sample")]:
        if field not in df:
            df[field] = ""
        replacement = df.sample_id if default is None else default
        df[field] = df[field].mask(df[field] == "", replacement)
    if not set(df.control) <= {"sample", "live", "dead"}:
        raise ValueError("control must be sample, live, or dead")
    for _, unit in df.groupby(["batch_id", "unit_id"], sort=False):
        for field in ["condition", "biological_replicate", "control"]:
            if unit[field].nunique() != 1:
                raise ValueError(f"A unit_id must map to one {field} within a batch")
    for field in PATH_FIELDS:
        if field not in df:
            df[field] = ""
        df[field] = df[field].map(lambda x: str((path.parent / x).resolve()) if x else "")
        for value in df[field]:
            if value and not Path(value).is_file():
                raise FileNotFoundError(f"{field}: {value}")
    for field in ["time_index", "series_index"]:
        if field not in df:
            df[field] = ""
        for value in df[field]:
            if value and (not value.isdigit()):
                raise ValueError(f"{field} must be a nonnegative integer")
    identity = pd.DataFrame({"path": df.image_path,
                             "t": df.time_index.replace("", "0"),
                             "s": df.series_index.replace("", "0")})
    if identity.duplicated().any():
        raise ValueError("The same TIFF series/timepoint appears twice; do not duplicate fields of view")
    for row in df.to_dict("records"):
        values = [row.get(f"spacing_{axis}_um", "") for axis in "zyx"]
        if any(values):
            if not all(values):
                raise ValueError("Provide all three spacing_z_um, spacing_y_um, spacing_x_um values")
            spacing = np.asarray(values, dtype=float)
            if not np.isfinite(spacing).all() or (spacing <= 0).any():
                raise ValueError("Voxel spacing must be finite and positive, in micrometres")
    return df


def canonical_czyx(array: np.ndarray, axes: str, time_index: int | None = None) -> np.ndarray:
    axes = axes.upper()
    if len(axes) != array.ndim or len(set(axes)) != len(axes):
        raise ValueError(f"Axes '{axes}' do not uniquely describe shape {array.shape}")
    if "T" in axes:
        axis = axes.index("T")
        if array.shape[axis] > 1 and time_index is None:
            raise ValueError("Multiple timepoints: explicitly specify time_index in the manifest")
        index = 0 if time_index is None else time_index
        if index < 0 or index >= array.shape[axis]:
            raise ValueError("time_index is outside the TIFF time axis")
        array = np.take(array, index, axis=axis)
        axes = axes.replace("T", "")
    elif time_index not in {None, 0}:
        raise ValueError("time_index was given for a TIFF without a time axis")
    for axis in range(len(axes) - 1, -1, -1):
        if axes[axis] not in "CZYX":
            if array.shape[axis] != 1:
                raise ValueError(f"Ambiguous TIFF axis '{axes[axis]}' in {axes}. Set an explicit axes override; do not guess depth or RGB channels.")
            array = np.take(array, 0, axis=axis)
            axes = axes[:axis] + axes[axis + 1:]
    if not all(axis in axes for axis in "ZYX"):
        raise ValueError(f"A genuine 3D Z-stack is required; found axes={axes}")
    if array.shape[axes.index("Z")] < 3:
        raise ValueError("At least three Z planes are required; coarse stacks may still fail morphology QC")
    if "C" not in axes:
        array = array[np.newaxis]
        axes = "C" + axes
    array = np.transpose(array, [axes.index(axis) for axis in "CZYX"])
    if not np.issubdtype(array.dtype, np.number) and array.dtype != bool:
        raise ValueError("Only numeric microscopy images are supported")
    if np.issubdtype(array.dtype, np.inexact) and not np.isfinite(array).all():
        raise ValueError("Image contains NaN or infinite voxels")
    return array


def _unit_scale(unit: str) -> float:
    scales = {"µm": 1., "μm": 1., "um": 1., "micrometer": 1.,
              "nm": .001, "mm": 1000., "m": 1e6}
    if unit not in scales:
        raise ValueError(f"Unsupported OME physical unit: {unit}")
    return scales[unit]


def ome_spacing(xml: str | None, series_index: int = 0, time_index: int | None = None) -> tuple | None:
    if not xml:
        return None
    pixels = ET.fromstring(xml).findall(".//{*}Image/{*}Pixels")
    if series_index >= len(pixels):
        return None
    pixel = pixels[series_index]
    values = []
    for axis in "ZYX":
        value = pixel.get(f"PhysicalSize{axis}")
        if value is None:
            return None
        values.append(float(value) * _unit_scale(pixel.get(f"PhysicalSize{axis}Unit", "µm")))
    values = np.asarray(values, float)
    if not np.isfinite(values).all() or (values <= 0).any():
        raise ValueError("Invalid OME physical spacing")
    # When plane positions are available, reject a nonuniform Z grid.
    positions = {}
    for plane in pixel.findall("{*}Plane"):
        if int(plane.get("TheT", "0")) == (time_index or 0) and int(plane.get("TheC", "0")) == 0:
            if plane.get("PositionZ") is not None:
                positions[int(plane.get("TheZ", "0"))] = float(plane.get("PositionZ")) * _unit_scale(plane.get("PositionZUnit", "µm"))
    if len(positions) >= 3:
        indices = np.array(sorted(positions))
        zvalues = np.array([positions[k] for k in indices])
        steps = np.abs(np.diff(zvalues) / np.diff(indices))
        if not np.allclose(steps, values[0], rtol=SPACING_RTOL, atol=SPACING_ATOL_UM):
            raise ValueError("OME plane positions disagree with a uniformly spaced Z grid; resample before analysis")
    return tuple(values)


def read_tiff(path: str | Path, axes: str = "", time_index: int | None = None,
              series_index: int = 0) -> tuple[np.ndarray, tuple | None, str]:
    with tifffile.TiffFile(path) as handle:
        if series_index >= len(handle.series):
            raise ValueError("series_index is outside the TIFF series list")
        series = handle.series[series_index]
        source_axes = axes or series.axes
        array = canonical_czyx(series.asarray(), source_axes, time_index)
        spacing = ome_spacing(handle.ome_metadata, series_index, time_index)
    return array, spacing, source_axes


@dataclass
class Sample:
    structure: np.ndarray
    calcein: np.ndarray | None
    pi: np.ndarray | None
    probability: np.ndarray | None
    imported_labels: np.ndarray | None
    spacing: tuple[float, float, float]
    metadata: dict


def load_sample(row: dict, cfg: dict) -> Sample:
    time_index = int(row["time_index"]) if row.get("time_index") else None
    series_index = int(row.get("series_index") or 0)
    main_key = (row["image_path"], row.get("axes", ""), time_index, series_index)
    main, metadata_spacing, source_axes = read_tiff(*main_key)
    # Different roles can select channels from one separate multichannel TIFF.
    # Cache each distinct series once per sample instead of decoding it per role.
    loaded = {main_key: (main, metadata_spacing, source_axes)}
    explicit = [row.get(f"spacing_{axis}_um", "") for axis in "zyx"]
    spacing = tuple(float(x) for x in explicit) if all(explicit) else metadata_spacing
    if spacing is None:
        raise ValueError("Physical spacing is missing: provide OME metadata or all spacing_*_um columns")
    if metadata_spacing and all(explicit) and not np.allclose(spacing, metadata_spacing, rtol=SPACING_RTOL, atol=SPACING_ATOL_UM):
        raise ValueError(f"Manifest spacing {spacing} conflicts with OME spacing {metadata_spacing}; correct the source metadata or manifest")
    volumes = {}
    channel_sources = {}
    for role in ["structure", "calcein", "pi", "probability", "labels"]:
        separate = row.get(f"{role}_path", "")
        if separate:
            key = (separate, row.get(f"{role}_axes", ""), time_index, 0)
            if key not in loaded:
                loaded[key] = read_tiff(*key)
            volume, other_spacing, _ = loaded[key]
            index = int(row.get(f"{role}_channel") or 0)
            if other_spacing and not np.allclose(spacing, other_spacing, rtol=SPACING_RTOL, atol=SPACING_ATOL_UM):
                raise ValueError(f"{role} TIFF spacing differs from the primary image")
        else:
            index = cfg["channels"].get(role)
            volume = main
        if index is None:
            volumes[role] = None
            continue
        source = (str(Path(separate or row["image_path"]).resolve()), 0 if separate else series_index, time_index or 0, index)
        if role in {"structure", "calcein", "pi"}:
            for previous_role, previous_source in channel_sources.items():
                if source == previous_source:
                    raise ValueError(f"{role} and {previous_role} point to the same channel; use independent structural and viability inputs")
            channel_sources[role] = source
        if index < 0 or index >= volume.shape[0]:
            raise ValueError(f"{role} channel {index} does not exist; image has {volume.shape[0]} channel(s)")
        volumes[role] = volume[index]
        if volumes[role].shape != main.shape[1:]:
            raise ValueError(f"{role} stack shape differs from primary image; register/resample channels first")
    if volumes["structure"] is None:
        raise ValueError("Supply a structure channel or structure_path, including when importing labels")
    if cfg["segmentation"]["method"] == "probability" and volumes["probability"] is None:
        raise ValueError("probability mode requires a probability_path column")
    if cfg["segmentation"]["method"] == "labels" and volumes["labels"] is None:
        raise ValueError("labels mode requires a labels_path column")
    metadata = {"input_axes": source_axes, "shape_zyx": list(main.shape[1:]),
                "spacing_source": "manifest" if all(explicit) else "OME", "spacing_zyx_um": list(spacing)}
    return Sample(volumes["structure"], volumes["calcein"], volumes["pi"], volumes["probability"], volumes["labels"], spacing, metadata)


def load_truth_labels(row: dict, expected_shape: tuple[int, int, int], spacing: tuple[float, float, float]) -> np.ndarray | None:
    """Load an optional annotated instance mask for post-segmentation validation.

    The annotation is intentionally separate from ``labels_path``: the latter is
    a possible *input* segmentation, while this field is reserved for an
    independent reference used only to calculate validation metrics.
    """
    path = row.get("truth_labels_path", "")
    if not path:
        return None
    if path == row.get("labels_path", ""):
        raise ValueError("truth_labels_path must be independent from labels_path")
    labels, annotation_spacing, _ = read_tiff(path, row.get("truth_labels_axes", ""),
                                               int(row["time_index"]) if row.get("time_index") else None)
    if labels.shape[0] != 1:
        raise ValueError("truth_labels_path must contain exactly one instance-label channel")
    truth = labels[0]
    if truth.shape != expected_shape:
        raise ValueError("truth_labels_path must be registered on the primary image ZYX grid")
    if not np.issubdtype(truth.dtype, np.integer) or (truth < 0).any():
        raise ValueError("truth_labels_path must contain nonnegative integer instance labels")
    if annotation_spacing and not np.allclose(spacing, annotation_spacing, rtol=SPACING_RTOL, atol=SPACING_ATOL_UM):
        raise ValueError("truth_labels_path voxel spacing differs from the primary image")
    return truth


def write_labels(path: Path, labels: np.ndarray, spacing: tuple) -> None:
    tifffile.imwrite(path, labels.astype(np.uint32), ome=True, photometric="minisblack",
                     compression="zlib", metadata={"axes": "ZYX", "PhysicalSizeZ": spacing[0],
                     "PhysicalSizeY": spacing[1], "PhysicalSizeX": spacing[2],
                     "PhysicalSizeZUnit": "µm", "PhysicalSizeYUnit": "µm", "PhysicalSizeXUnit": "µm"})


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_commit_hash(repo_root: Path | None = None) -> str | None:
    """Best-effort current git commit hash, with a "-dirty" suffix when the
    working tree has uncommitted changes. Returns None (never raises) when
    not run from a git repository or git is unavailable -- provenance capture
    must never fail a run.
    """
    root = repo_root or Path(__file__).resolve().parents[3]
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, capture_output=True,
                                text=True, timeout=5, check=False)
        if commit.returncode != 0 or not commit.stdout.strip():
            return None
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=root, capture_output=True,
                               text=True, timeout=5, check=False)
        suffix = "-dirty" if dirty.returncode == 0 and dirty.stdout.strip() else ""
        return commit.stdout.strip() + suffix
    except (OSError, subprocess.SubprocessError):
        return None
