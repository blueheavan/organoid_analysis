"""Provenance of instance masks: were any objects drawn or corrected by hand?

Automatic segmentation corrected by hand where it failed is common practice
(e.g. Kim et al. 2026, Cell Rep. Methods 6:101601, Figure S7), but a manually
edited object is a different measurement process from an automatic one: it
carries operator dependence and is not covered by any automatic-segmentation
validation (SG-3). This module never changes a measurement. It attaches, per
object, whether its mask is known to be automatic, known to be edited, or of
unknown origin, together with the SHA-256 of the label array it came from, so
that results can be stratified or re-run without edited objects.

Declarations
------------
Masks produced by this pipeline in the same run are ``pipeline_segmentation``
and therefore not edited. Imported labels carry what the user declares:

* ``none``      -- no object was edited (the imported masks are automatic);
* ``edited``    -- some objects were edited; list their original label IDs in
                   ``edited_ids`` to flag only those, otherwise every object is
                   flagged;
* ``manual``    -- every object was drawn by hand;
* not declared  -- origin unknown; every object is flagged ``unknown``.
"""
from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import numpy as np

DECLARATIONS = ("none", "edited", "manual")
EDIT_STATES = ("no", "yes", "unknown")
PROVENANCE_COLUMNS = ["mask_origin", "mask_edit_declaration", "mask_manually_edited", "labels_sha256"]
LEVELS = ("organoid", "cell", "nucleus")


def labels_sha256(labels: np.ndarray) -> str:
    """Content hash of a label array, independent of the file it was stored in."""
    array = np.ascontiguousarray(labels)
    digest = hashlib.sha256()
    digest.update(f"{array.dtype.str}|{array.shape}".encode())
    digest.update(array.tobytes())
    return digest.hexdigest()


def parse_ids(value: object) -> frozenset[int] | None:
    """Parse ``"3;7;12"`` (or an iterable of ints) into label IDs; empty -> None."""
    if value is None:
        return None
    if isinstance(value, str):
        parts = [part.strip() for part in value.replace(",", ";").split(";") if part.strip()]
    elif isinstance(value, Iterable):
        parts = [str(part) for part in value]
    else:
        raise ValueError(f"edited_ids must be a list of label IDs, got {value!r}")
    if not parts:
        return None
    ids = []
    for part in parts:
        if not part.isdigit() or int(part) <= 0:
            raise ValueError(f"edited_ids must be positive integer label IDs, got {part!r}")
        ids.append(int(part))
    return frozenset(ids)


@dataclass(frozen=True)
class MaskProvenance:
    origin: str                              # pipeline_segmentation | imported_labels
    declaration: str                         # none | edited | manual | not_declared | not_applicable
    edited_ids: frozenset[int] | None = None
    sha256: str = ""

    @classmethod
    def pipeline(cls, labels: np.ndarray) -> MaskProvenance:
        return cls("pipeline_segmentation", "not_applicable", None, labels_sha256(labels))

    @classmethod
    def imported(cls, labels: np.ndarray, declaration: str = "", edited_ids: object = None) -> MaskProvenance:
        declaration = (declaration or "").strip().lower()
        if declaration and declaration not in DECLARATIONS:
            raise ValueError(f"manual-edit declaration must be one of {DECLARATIONS} or empty, got {declaration!r}")
        ids = parse_ids(edited_ids)
        if ids is not None and declaration != "edited":
            raise ValueError("edited_ids may only accompany the declaration 'edited'")
        present = set(np.unique(labels[labels > 0]).astype(int).tolist())
        if ids is not None and not ids <= present:
            missing = sorted(ids - present)
            raise ValueError(f"edited_ids not present in the label image: {missing[:10]}")
        return cls("imported_labels", declaration or "not_declared", ids, labels_sha256(labels))

    def state(self, label_id: int) -> str:
        """``no``, ``yes`` or ``unknown`` for one original label ID."""
        if self.origin == "pipeline_segmentation" or self.declaration == "none":
            return "no"
        if self.declaration == "manual":
            return "yes"
        if self.declaration == "edited":
            return "yes" if self.edited_ids is None or int(label_id) in self.edited_ids else "no"
        return "unknown"

    def columns(self, label_id: int) -> dict[str, str]:
        return {"mask_origin": self.origin, "mask_edit_declaration": self.declaration,
                "mask_manually_edited": self.state(label_id), "labels_sha256": self.sha256}

    def summary(self, label_ids: Iterable[int]) -> dict[str, object]:
        states = [self.state(i) for i in label_ids]
        return {"origin": self.origin, "declaration": self.declaration, "labels_sha256": self.sha256,
                "n_objects": len(states), **{f"n_edited_{s}": states.count(s) for s in EDIT_STATES},
                "edited_ids": sorted(self.edited_ids) if self.edited_ids is not None else None}


def load_edit_log(path: str | Path) -> dict[str, dict[str, object]]:
    """Read an ``analyze-3d`` edit log: ``{"organoid": {"declaration": "edited", "edited_ids": [3, 7]}, ...}``."""
    document = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(document, dict) or not set(document) <= set(LEVELS):
        raise ValueError(f"the label edit log must be a JSON object keyed by {LEVELS}")
    for level, entry in document.items():
        if not isinstance(entry, dict) or not set(entry) <= {"declaration", "edited_ids"}:
            raise ValueError(f"{level}: entries take only 'declaration' and 'edited_ids'")
    return document


def stratification_note(frame_states: Iterable[str]) -> str | None:
    """Warning text when any object is edited or of unknown origin, else None."""
    states = list(frame_states)
    edited, unknown = states.count("yes"), states.count("unknown")
    if not edited and not unknown:
        return None
    return (f"{edited} object(s) have manually edited masks and {unknown} have masks of undeclared origin. "
            "Automatic-segmentation validation does not cover them; report results with and without them "
            "(column mask_manually_edited).")
