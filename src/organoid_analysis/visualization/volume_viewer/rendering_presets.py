"""Preset color transfer functions for vtk.js volume rendering.

Each preset is a list of ``(position, r, g, b)`` control points on the
normalized ``[0, 1]`` scalar axis. These are shipped to the browser and used to
build a ``vtkColorTransferFunction`` for a channel.

Cyan and magenta are the classical fluorescence channel colors (GFP / RFP
overlays); the rest are useful alternatives.
"""

from __future__ import annotations

from collections.abc import Sequence

LutPoint = tuple[float, float, float, float]
Lut = Sequence[LutPoint]

PRESETS: dict[str, Lut] = {
    # (position, r, g, b) on [0,1]
    "cyan": [
        (0.0, 0.0, 0.0, 0.0),
        (1.0, 0.0, 1.0, 1.0),
    ],
    "magenta": [
        (0.0, 0.0, 0.0, 0.0),
        (1.0, 1.0, 0.0, 1.0),
    ],
    "green": [
        (0.0, 0.0, 0.0, 0.0),
        (1.0, 0.0, 1.0, 0.0),
    ],
    "red": [
        (0.0, 0.0, 0.0, 0.0),
        (1.0, 1.0, 0.0, 0.0),
    ],
    "blue": [
        (0.0, 0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0, 1.0),
    ],
    "yellow": [
        (0.0, 0.0, 0.0, 0.0),
        (1.0, 1.0, 1.0, 0.0),
    ],
    "gray": [
        (0.0, 0.0, 0.0, 0.0),
        (1.0, 1.0, 1.0, 1.0),
    ],
}

DEFAULT_LUTS: tuple[str, ...] = ("cyan", "magenta", "yellow", "green")
