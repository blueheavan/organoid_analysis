from pathlib import Path
from typing import Any


# Analysis input now comes from pipelines' exported results (e.g. results/),
# not an external data/ folder. DATA_DIR is kept only as a legacy constant so
# existing importers do not break; it is no longer auto-created on import.
# paths.py lives under src/segmentation/, so the repo root is three levels up.
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
IMAGES_DIR = PROJECT_ROOT / "data" / "images"
DATA_DIR = PROJECT_ROOT / "data"
DATA_DIR_EXISTS = DATA_DIR.is_dir()
OUTPUTS_DIR = PROJECT_ROOT / "results"
SEGMENTATION_OUTPUT_DIR = OUTPUTS_DIR / "segmentation_output"
FEATURES_DIR = OUTPUTS_DIR / "features"


for path in (OUTPUTS_DIR, SEGMENTATION_OUTPUT_DIR, FEATURES_DIR):
    path.mkdir(parents=True, exist_ok=True)


def detect_torch_acceleration(torch_module: Any) -> tuple[str, bool]:
    mps_backend = getattr(torch_module.backends, "mps", None)
    if mps_backend is not None and mps_backend.is_available():
        return "mps", True
    if torch_module.cuda.is_available():
        return "cuda", True
    return "cpu", False


def require_file(path_like: str | Path) -> Path:
    path = Path(path_like)
    if not path.exists():
        raise FileNotFoundError(
            f"Missing required file: {path}\n"
            "Check the project's data/images/ or results/ directory."
        )
    return path
