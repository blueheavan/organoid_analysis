"""Shared fixtures."""
from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from organoid_analysis.validation.analytical_geometry_evidence import (
    current_manifest_relpath,
    referenced_paths,
)
from organoid_analysis.validation.evidence_manifest import load_manifest

PROJECT_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def evidence_copy(tmp_path: Path) -> tuple[Path, str]:
    """Byte-identical copy of the current analytical-geometry record and every file it depends on.

    Tests tamper with the copy, never with the repository's evidence.
    """
    relpath = current_manifest_relpath(PROJECT_ROOT)
    manifest = load_manifest(PROJECT_ROOT / relpath)
    for path in sorted(referenced_paths(manifest, relpath)):
        target = tmp_path / path
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(PROJECT_ROOT / path, target)
    return tmp_path, relpath
