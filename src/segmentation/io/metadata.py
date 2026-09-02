"""Metadata parsing for Z-stack TIFF files (independent of tifffile loading).

Parses voxel spacing (X, Y, Z physical sizes, normalised to micrometres) and
axis semantics (which dimension is Z / C / T / Y / X).

Spacing source priority:
    1. OME-XML  (PixelPhysicalSizeX/Y/Z + units)          -- most reliable
    2. ImageJ   (spacing + unit; XY from TIFF resolution) -- second
    3. none found -> metadata returns no spacing (caller prompts the user)

We deliberately never invent a default spacing (e.g. never assume 1.0 um).

OME-XML is parsed with ``ome_types`` -- never hand-rolled XML regex.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

# μm conversion factors: value_in_unit * FACTOR = value_in_um
# OME's default unit when PhysicalSize*Unit is absent/empty is micrometre.
_UNIT_TO_UM: dict[str, float] = {
    "": 1.0,  # bare numbers / missing unit default to micrometres per OME
    "m": 1.0e6,
    "cm": 1.0e4,
    "mm": 1.0e3,
    "µm": 1.0,
    "um": 1.0,
    "nm": 1.0e-3,
    "Å": 1.0e-4,
}


def _to_um(value: float, unit: str | None) -> float:
    unit = (unit or "").strip() or "µm"
    factor = _UNIT_TO_UM.get(unit)
    if factor is None:
        # Unknown unit: be conservative and treat as metres only if explicitly SI.
        raise ValueError(f"Unsupported physical-size unit: {unit!r}")
    return value * factor


@dataclass(frozen=True)
class Spacing:
    """Physical voxel spacing in micrometres. ``None`` means unknown (not 1.0)."""

    x: float | None
    y: float | None
    z: float | None

    @property
    def complete(self) -> bool:
        return self.x is not None and self.y is not None and self.z is not None

    @property
    def anisotropy(self) -> float | None:
        """Z / XY ratio if both are known (XY = mean of X and Y)."""
        if self.z is None or self.x is None or self.y is None:
            return None
        xy = (self.x + self.y) / 2.0
        if xy <= 0:
            return None
        return self.z / xy


def parse_spacing_ome(ome_metadata: str | None) -> Spacing:
    """Parse voxel spacing from OME-XML via ``ome_types``.

    Returns a ``Spacing`` with per-axis values in micrometres. Any axis that
    cannot be resolved is ``None`` (unknown). Raises if the XML is malformed in
    a way ``ome_types`` cannot parse.
    """
    x = y = z = None
    if not ome_metadata:
        return Spacing(x, y, z)

    from ome_types import from_xml

    ome = from_xml(ome_metadata)
    if not ome.images:
        return Spacing(x, y, z)

    pixels = ome.images[0].pixels
    for attr, slot in (
        ("physical_size_x", "x"),
        ("physical_size_y", "y"),
        ("physical_size_z", "z"),
    ):
        value = getattr(pixels, attr)
        unit = getattr(pixels, f"{attr}_unit")
        if value is None:
            continue
        try:
            um = _to_um(float(value), unit.value if unit is not None else "µm")
        except ValueError:
            continue
        if slot == "x":
            x = um
        elif slot == "y":
            y = um
        elif slot == "z":
            z = um
    return Spacing(x, y, z)


def parse_spacing_imagej(imagej_metadata: dict[str, Any] | None) -> Spacing:
    """Parse spacing from ImageJ hyperstack metadata.

    ImageJ usually stores ``spacing`` (Z step) and ``unit``; XY pixel size must
    come from the TIFF XResolution/YResolution tags (handled in
    :func:`resolve_spacing`). Here we only extract what the ImageJ metadata
    yields directly (typically Z).
    """
    z = None
    if not imagej_metadata:
        return Spacing(None, None, z)
    spacing = imagej_metadata.get("spacing")
    if spacing is not None:
        try:
            z = float(spacing)
        except (TypeError, ValueError):
            z = None
    return Spacing(None, None, z)


def _resolution_to_um_per_px(resolution: float | tuple | None, unit: int | None) -> float | None:
    """TIFF stores resolution as pixels-per-unit; invert to µm/pixel.

    ``unit`` is the TIFF ResolutionUnit: 1=undefined, 2=inch, 3=centimetre,
    4=millimetre, 5=micrometre.
    """
    if resolution is None:
        return None
    if isinstance(resolution, (tuple, list)):
        if not resolution or not resolution[0]:
            return None
        px_per_unit = float(resolution[0]) / float(resolution[1])
    else:
        px_per_unit = float(resolution)
    if px_per_unit <= 0:
        return None
    if unit == 3:  # pixels per centimetre -> µm per pixel
        return (1.0e4) / px_per_unit
    if unit == 4:  # pixels per millimetre
        return (1.0e3) / px_per_unit
    if unit == 5:  # pixels per micrometre
        return 1.0 / px_per_unit
    if unit == 2:  # pixels per inch -> 25.4mm/inch = 25400 µm/inch
        return (25400.0) / px_per_unit
    # unit 1 = undefined / no absolute unit (tifffile's default when a file,
    # e.g. a stacked Z-series, carries no real physical resolution). Without an
    # absolute unit the stored density is meaningless, so treat spacing as
    # unknown rather than inventing a distorted value.
    return None


def resolve_spacing(
    ome_metadata: str | None,
    imagej_metadata: dict[str, Any] | None,
    tiff_tags: dict[str, Any] | None = None,
) -> Spacing:
    """Combine the best available spacing sources into a single ``Spacing``.

    Priority per axis: OME > ImageJ Z-step (for Z) / TIFF resolution (for X,Y)
    > unknown. ``tiff_tags`` is a dict with ``XResolution``/``YResolution`` and
    ``ResolutionUnit`` keys (as read from ``tifffile.TiffPage.tags``).
    """
    spacing = parse_spacing_ome(ome_metadata)
    if not spacing.complete:
        imagej = parse_spacing_imagej(imagej_metadata)
        tags = tiff_tags or {}
        xres = tags.get("XResolution")
        yres = tags.get("YResolution")
        resunit = tags.get("ResolutionUnit")
        # OME wins per-axis; fill gaps from ImageJ/TIFF.
        sx = spacing.x if spacing.x is not None else _resolution_to_um_per_px(xres, resunit)
        sy = spacing.y if spacing.y is not None else _resolution_to_um_per_px(yres, resunit)
        sz = spacing.z if spacing.z is not None else imagej.z
        spacing = Spacing(sx, sy, sz)
    return spacing


def parse_axes(axes: str | None) -> dict[str, int]:
    """Map dimension labels to their index in the TIFF series.

    ``axes`` is the tifffile series axes string (e.g. ``"ZYX"``, ``"TZCYX"``,
    ``"CYX"``, ``"QYX"``). Returns a dict ``{label: index}``.
    """
    axes = axes or ""
    return {ch: i for i, ch in enumerate(axes)}


def classify_axes(axes: str | None) -> dict[str, int]:
    """Return which axis indices are Z, C, T, Y, X given the axes string.

    This is the semantic resolution used by the reader to decide how to
    rearrange a series into our canonical layouts:
      * no channel: volume.shape == (Z, Y, X)
      * channel:    volume.shape == (C, Z, Y, X)
    Returns ``{"Z": idx, "Y": idx, "X": idx, "C": idx|None, "T": idx|None,
    "S": idx|None, "Q": idx|None}``.
    """
    labels = parse_axes(axes)
    out = {
        "Z": labels.get("Z"),
        "Y": labels.get("Y"),
        "X": labels.get("X"),
        "C": labels.get("C"),
        "T": labels.get("T"),
        "S": labels.get("S"),
        "Q": labels.get("Q"),
    }
    # 'Q' (quarto, RGB) and 'S' (sample) mean samples-per-pixel > 1 -> reject.
    return out


def reject_rgb(axes: str | None, samples_per_pixel: int | None = None) -> bool:
    """True if the image is RGB/RGBA (samples-per-pixel > 1) -- not quantitative.

    Z-stack intensity data must be a single grey channel per X/Y location.
    """
    if samples_per_pixel and samples_per_pixel > 1:
        return True
    labels = classify_axes(axes)
    return labels["S"] is not None or labels["Q"] is not None
