"""Headless multilevel 3D organoid analysis.

This package accepts already segmented, registered instance masks.  It has no
dependency on a segmentation model or UI, which keeps measurement reproducible
and makes it usable by the CLI, batch workflows, and the Streamlit frontend.

The orchestration entry point (``Multilevel3DResult``, ``analyze_multilevel_3d``)
lives in :mod:`organoid_analysis.workflows.multilevel_measurement_workflow`,
which wires this package's hierarchy/morphology/spatial/topology/qc modules
together; import it from there rather than from here, to avoid a circular
import between the quantification and workflows layers.
"""

from .config import Multilevel3DConfig

__all__ = ["Multilevel3DConfig"]
