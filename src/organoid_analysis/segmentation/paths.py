import os
from pathlib import Path
from typing import Any

# ORGANOID_ANALYSIS_PROJECT_ROOT lets a CLI wrapper, deployment config, or a
# user running outside the pixi-managed checkout point these paths at the
# right place explicitly, instead of only ever trusting this file's location
# in the package tree. That inference is still the default (nothing else to
# fall back on for a bare `pixi run web`/`pytest`), but a prior package move
# already broke it once silently (an extra directory level shifted "four
# parents up" from the repo root to src/, which created stray
# src/results/{segmentation_output,features}/ directories on import before
# anyone noticed) -- the env var is the escape hatch for the next time the
# package layout changes, and for anyone who wants explicit control instead
# of implicit inference.
_ROOT_ENV_VAR = "ORGANOID_ANALYSIS_PROJECT_ROOT"


def _resolve_project_root() -> Path:
    override = os.environ.get(_ROOT_ENV_VAR)
    if override:
        return Path(override).resolve()
    # paths.py lives under src/organoid_analysis/segmentation/, so the repo
    # root is four levels up.
    return Path(__file__).resolve().parent.parent.parent.parent


PROJECT_ROOT = _resolve_project_root()
# Analysis input now comes from pipelines' exported results (e.g. results/),
# not an external data/ folder. DATA_DIR is kept only as a legacy constant so
# existing importers do not break; it is no longer auto-created on import.
DATA_DIR = PROJECT_ROOT / "data"
DATA_DIR_EXISTS = DATA_DIR.is_dir()
OUTPUTS_DIR = PROJECT_ROOT / "results"
SEGMENTATION_OUTPUT_DIR = OUTPUTS_DIR / "segmentation_output"
# Actual directory creation happens lazily where these paths are used
# (save_result() creates its own run directory with parents=True;
# list_saved_results() treats a missing directory as "no runs yet") -- import
# never mutates the filesystem.


def detect_torch_acceleration(torch_module: Any) -> tuple[str, bool]:
    mps_backend = getattr(torch_module.backends, "mps", None)
    if mps_backend is not None and mps_backend.is_available():
        return "mps", True
    if torch_module.cuda.is_available():
        return "cuda", True
    return "cpu", False
