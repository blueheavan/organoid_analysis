"""Report the execution environment that determines CPU/MPS/CUDA and OpenMP.

    pixi run python scripts/runtime_diagnostics.py

This is a diagnostic, not a validation step and not a gate. It exists because
the public-data trial (docs/REAL_DATA_TRIALS_2026-09-14.md) could not determine
whether segmentation ran on the accelerator: a torch probe aborted on a
duplicate OpenMP runtime before reporting device status. This script reports
the resolved accelerator, the OpenMP inputs that can cause that abort, and the
presence of the unsafe ``KMP_DUPLICATE_LIB_OK`` override -- so device provenance
is observable rather than inferred.
"""
from __future__ import annotations

import os
import platform
import sys
from pathlib import Path


def _openmp_providers(env_root: Path) -> list[str]:
    candidates = [
        env_root / "lib" / "libomp.dylib",
        env_root / "lib" / "libgomp.dylib",
        env_root / "lib" / "libiomp5.dylib",
    ]
    site_packages = env_root / "lib" / f"python{sys.version_info.major}.{sys.version_info.minor}" / "site-packages"
    if site_packages.is_dir():
        candidates += sorted(site_packages.rglob("libomp.dylib"))
        candidates += sorted(site_packages.rglob("libgomp*.dylib"))
    return [str(path) for path in candidates if path.is_file()]


def main() -> int:
    print("Runtime diagnostics (diagnostic only; not a validation gate)")
    print(f"python:   {sys.version.split()[0]}")
    print(f"platform: {platform.platform()}")
    print(f"machine:  {platform.machine()}")

    override = os.environ.get("KMP_DUPLICATE_LIB_OK")
    print(f"KMP_DUPLICATE_LIB_OK: {override!r}"
          + ("  <- UNSAFE override is set" if override else "  (unset)"))

    env_root = Path(sys.prefix)
    providers = _openmp_providers(env_root)
    print(f"OpenMP runtimes found in {env_root}: {len(providers)}")
    for path in providers:
        print(f"  - {path}")
    if len(providers) > 1:
        print("  NOTE: more than one OpenMP runtime is present. Import torch after the")
        print("        scientific stack (numpy/scipy) and do not set KMP_DUPLICATE_LIB_OK.")

    try:
        import numpy  # noqa: F401 - imported before torch by design; see paths.detect_torch_acceleration
        import torch

        from organoid_analysis.segmentation.paths import detect_torch_acceleration
    except Exception as error:  # noqa: BLE001 - diagnostics must report, not crash mysteriously
        print(f"torch import/device detection failed: {error!r}")
        return 1

    accelerator, accelerated = detect_torch_acceleration(torch)
    print(f"torch: {torch.__version__}")
    print(f"mps built: {getattr(torch.backends.mps, 'is_built', lambda: None)()}")
    print(f"detected accelerator: {accelerator} (accelerated={accelerated})")
    if accelerator == "cpu":
        print("  NOTE: CPU was selected. A CPU-only run is an environment-performance")
        print("        blocker for DINOv3 3D inference, not a scientific result.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
