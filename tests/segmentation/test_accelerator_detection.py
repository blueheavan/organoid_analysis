"""Regression tests for accelerator detection and the OpenMP override removal.

Two independent defects are guarded here:

1. ``detect_torch_acceleration`` had no test, so its precedence (MPS before
   CUDA before CPU) and its tolerance of a torch build without ``backends.mps``
   were unverified.
2. The Web entry points set ``KMP_DUPLICATE_LIB_OK=TRUE`` unconditionally. That
   is documented by the OpenMP runtime as unsafe and can silently produce
   incorrect results. It was removed after the import order was verified; this
   test fails if it is reintroduced anywhere in the production package.
"""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from organoid_analysis.segmentation.paths import detect_torch_acceleration

PACKAGE_ROOT = Path(__file__).resolve().parents[2] / "src" / "organoid_analysis"


def _fake_torch(*, mps_available: bool, cuda_available: bool, has_mps: bool = True) -> object:
    backends = SimpleNamespace()
    if has_mps:
        backends.mps = SimpleNamespace(is_available=lambda: mps_available)
    return SimpleNamespace(
        backends=backends,
        cuda=SimpleNamespace(is_available=lambda: cuda_available),
    )


def test_mps_takes_precedence_over_cuda() -> None:
    assert detect_torch_acceleration(_fake_torch(mps_available=True, cuda_available=True)) == ("mps", True)


def test_cuda_is_used_when_mps_is_unavailable() -> None:
    assert detect_torch_acceleration(_fake_torch(mps_available=False, cuda_available=True)) == ("cuda", True)


def test_cpu_is_the_unaccelerated_fallback() -> None:
    assert detect_torch_acceleration(_fake_torch(mps_available=False, cuda_available=False)) == ("cpu", False)


def test_torch_build_without_mps_backend_does_not_crash() -> None:
    assert detect_torch_acceleration(
        _fake_torch(mps_available=False, cuda_available=False, has_mps=False)
    ) == ("cpu", False)


def _sets_openmp_override(path: Path) -> bool:
    """True when a non-comment source line references the unsafe override.

    Comments are ignored so the removal can be documented in place; an actual
    ``os.environ`` assignment or ``setdefault`` is what this guards against.
    """
    for line in path.read_text(encoding="utf-8").splitlines():
        code = line.split("#", 1)[0]
        if "KMP_DUPLICATE_LIB_OK" in code:
            return True
    return False


def test_production_package_does_not_set_the_unsafe_openmp_override() -> None:
    offenders = [
        str(path.relative_to(PACKAGE_ROOT.parents[1]))
        for path in PACKAGE_ROOT.rglob("*.py")
        if _sets_openmp_override(path)
    ]
    assert not offenders, (
        "KMP_DUPLICATE_LIB_OK is an unsafe OpenMP workaround and must not appear in "
        f"production code: {offenders}. Fix the runtime linkage instead."
    )
