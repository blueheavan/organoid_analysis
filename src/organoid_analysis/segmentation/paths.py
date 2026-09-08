from pathlib import Path
from typing import Any

# Analysis input now comes from pipelines' exported results (e.g. results/),
# not an external data/ folder. DATA_DIR is kept only as a legacy constant so
# existing importers do not break; it is no longer auto-created on import.
# paths.py lives under src/organoid_analysis/segmentation/, so the repo root
# is four levels up.
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
DATA_DIR = PROJECT_ROOT / "data"
DATA_DIR_EXISTS = DATA_DIR.is_dir()
OUTPUTS_DIR = PROJECT_ROOT / "results"
SEGMENTATION_OUTPUT_DIR = OUTPUTS_DIR / "segmentation_output"
FEATURES_DIR = OUTPUTS_DIR / "features"


def ensure_output_directories() -> None:
    """Create runtime output directories without mutating disk at import time."""
    for path in (OUTPUTS_DIR, SEGMENTATION_OUTPUT_DIR, FEATURES_DIR):
        path.mkdir(parents=True, exist_ok=True)


def detect_torch_acceleration(torch_module: Any) -> tuple[str, bool]:
    mps_backend = getattr(torch_module.backends, "mps", None)
    if mps_backend is not None and mps_backend.is_available():
        return "mps", True
    if torch_module.cuda.is_available():
        return "cuda", True
    return "cpu", False
