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


def _copy_record(tmp_path: Path, module) -> tuple[Path, str]:
    from organoid_analysis.validation.record_contract import current_manifest_relpath as relpath_of
    from organoid_analysis.validation.record_contract import verify_referenced

    try:
        relpath = relpath_of(PROJECT_ROOT, module.SPEC)
    except ValueError:
        pytest.skip(f"no {module.CONTRACT_ID} record in this checkout")
    manifest = load_manifest(PROJECT_ROOT / relpath)
    for path in sorted(verify_referenced(manifest, relpath, module.SPEC)):
        target = tmp_path / path
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(PROJECT_ROOT / path, target)
    return tmp_path, relpath


@pytest.fixture
def qualification_copy(tmp_path: Path) -> tuple[Path, str]:
    """Byte-identical copy of the current SG-1a/SG-2a qualification record and its inputs."""
    from organoid_analysis.validation import analytical_qualification_evidence

    return _copy_record(tmp_path, analytical_qualification_evidence)


@pytest.fixture
def viability_copy(tmp_path: Path) -> tuple[Path, str]:
    """Byte-identical copy of the current SG-4A record and its inputs."""
    from organoid_analysis.validation import viability_rule_evidence

    return _copy_record(tmp_path, viability_rule_evidence)
