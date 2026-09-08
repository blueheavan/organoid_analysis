"""Streamlit UI layer (the optional final step of the pipeline).

Depends on the headless ``segmentation`` and ``quantification`` packages.
Package layout:
* ``organoid_workspace``      — unified one-page entry point.
* ``segmentation_workspace``  — upload + Cellpose 3D segmentation UI.
* ``preview_2d``              — 2D (and server-side VTK) preview.
* ``analysis_ui``             — statistical analysis UI (Tutorials 2-5).
* ``multilevel_results``      — read-only view of exported multilevel results.
* ``demo_app``                — standalone vtk.js viewer demo.
"""

__version__ = "1.0.0"
