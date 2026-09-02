"""Unified 3D organoid pipeline — one page, full loop.

Wraps ``ui/vtk_viewer/integrated_app.py`` (image upload, Cellpose 3D
segmentation, browser-GPU vtk.js preview/mask overlay, and per-object mask
feature extraction) into a single entry point, and additionally exposes:

  1. Object features measured from the segmentation mask (``mask_features``).
  2. Full statistical workflows from exported Excel tables (``ui.analysis_ui``).

Run:
    streamlit run src/ui/app.py
(or ``pixi run app``)"""

from __future__ import annotations

import os
import sys
from pathlib import Path

# torch can abort at import if two OpenMP runtimes are linked; allow it so the
# app (and Cellpose segmentation) can start from any shell / launcher.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import streamlit as st

_PROJECT_ROOT = Path(__file__).resolve().parent.parent
# ``streamlit run src/ui/app.py`` puts ``src/ui`` first. Move ``src`` ahead of
# it so the core ``analysis`` package wins over the UI helper ``analysis.py``.
while str(_PROJECT_ROOT) in sys.path:
    sys.path.remove(str(_PROJECT_ROOT))
sys.path.insert(0, str(_PROJECT_ROOT))

from segmentation.cellpose import get_accelerator  # noqa: E402

try:  # noqa: E402  (importable from the ui.vtk_viewer package)
    from ui.vtk_viewer.integrated_app import (
        render_analysis_tab as _render_mask_stats,
    )
    from ui.vtk_viewer.integrated_app import (
        render_preview_tab as _render_preview,
    )
    from ui.vtk_viewer.integrated_app import (
        render_results_tab as _render_results,
    )
    from ui.vtk_viewer.integrated_app import render_sidebar as _render_sidebar
    HAS_VTK = True
except Exception:  # noqa: BLE001  (vtk viewer/GPU unavailable -> degraded mode)
    HAS_VTK = False

try:  # noqa: E402  (statistical analysis layer)
    from ui.analysis_ui import render_analysis as _render_tutorial_stats
    HAS_STATS = True
except Exception:  # noqa: BLE001
    HAS_STATS = False

from ui.multilevel_results import render_multilevel_results_section as _render_multilevel_results


def _summary_line() -> str:
    if HAS_VTK and HAS_STATS:
        return "GPU 3D viewer + segmentation + object features + Excel statistics."
    if HAS_VTK:
        return "GPU 3D viewer + segmentation + object features (statistics unavailable)."
    return "Degraded mode: only basic preview/statistics (vtk viewer unavailable)."


def main() -> None:
    st.set_page_config(page_title="Unified 3D Organoid Pipeline", layout="wide")
    st.title("3D Organoid Pipeline: Upload → Preview → Segment → Analyze")
    try:
        accelerator = get_accelerator().upper()
    except Exception:  # torch/Cellpose not importable -> preview-only mode
        accelerator = "unavailable"
    st.caption(f"Compute device: `{accelerator}`. {_summary_line()}")
    with st.expander("Analysis guide and interpretation limits", expanded=False):
        st.markdown(
            "**1. Upload and preview.** Confirm that the Z-stack has the correct orientation, "
            "the Z/XY anisotropy matches the microscope metadata, and the selected channel shows "
            "the structures you intend to segment. Display contrast is only for visualization.\n\n"
            "**2. Segmentation results.** Colored contours and translucent surfaces identify mask "
            "boundaries, not biological classes. Inspect missed objects, merges, splits and objects "
            "touching the image border before using any measurements. A Cellpose mask is not an "
            "independent accuracy validation; validate final settings against registered annotations "
            "with the headless analysis CLI when those annotations are available.\n\n"
            "**3. Object features.** Volumes, areas and axes use the entered physical voxel size. "
            "Sphericity and solidity describe the segmentation mask geometry, so changes in focus, "
            "thresholding or object clipping can change them without a biological change. One uploaded "
            "field produces descriptive measurements only, not treatment statistics.\n\n"
            "**4. Statistical analysis.** The Excel tools are exploratory tutorial workflows. Define "
            "the independent experimental unit (well, biological replicate, donor or batch) before "
            "interpreting p-values or classifier scores. Do not treat every cell or organoid in a field "
            "as an independent replicate."
        )

    # Shared segmentation settings (sidebar) only relevant when the vtk-backed
    # segmentation tabs are present.
    if HAS_VTK:
        config = _render_sidebar()
    else:
        config = None

    tab_preview, tab_results, tab_features, tab_multilevel, tab_stats = st.tabs(
        [
            "Upload & preview",
            "Segmentation results",
            "Object features",
            "3D Analysis results",
            "Statistical analysis (Excel)",
        ]
    )

    with tab_preview:
        if HAS_VTK:
            _render_preview(config)
        else:
            st.error("vtk.js viewer unavailable; cannot render the preview tab.")

    with tab_results:
        if HAS_VTK:
            _render_results(config)
        else:
            st.info("Run a segmentation (needs the vtk viewer) for results to appear here.")

    with tab_features:
        if HAS_VTK:
            _render_mask_stats()
        else:
            st.info("Per-object feature analysis requires the vtk-backed segmentation.")

    with tab_multilevel:
        _render_multilevel_results()

    with tab_stats:
        if HAS_STATS:
            _render_tutorial_stats()
        else:
            st.error("Statistical analysis layer unavailable.")


if __name__ == "__main__":
    main()
