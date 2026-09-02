"""Streamlit UI layer (the optional final step of the pipeline).

Depends on the headless ``segmentation`` and ``analysis`` packages.
Package layout:
* ``upload``       — image upload + Cellpose 3D segmentation UI.
* ``preview_2d``   — 2D (and server-side VTK) preview.
* ``vtk_viewer``   — browser-GPU (vtk.js) 3D volume viewer + mask overlay.
* ``analysis_ui``  — statistical analysis UI (Tutorials 2-5).
* ``analysis``     — reusable statistical analysis engine.
* ``app``          — unified one-page entry point.
"""

from pathlib import Path
import sys


def _prioritize_project_packages() -> None:
    """Place ``src`` before ``src/ui`` to avoid the ``analysis`` name collision.

    Streamlit adds the directory containing the launched script to ``sys.path``.
    When launched as ``streamlit run src/ui/app.py``, that makes the UI helper
    ``ui/analysis.py`` shadow the sibling core package ``analysis/``. The mask
    feature layer imports ``analysis.features``, so the core package must always
    be searched first regardless of which UI entry point launched the process.
    """
    project_packages = str(Path(__file__).resolve().parent.parent)
    while project_packages in sys.path:
        sys.path.remove(project_packages)
    sys.path.insert(0, project_packages)


_prioritize_project_packages()

__version__ = "1.0.0"
