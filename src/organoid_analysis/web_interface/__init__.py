"""Streamlit UI layer (the optional final step of the pipeline).

Depends on the headless ``segmentation`` and ``quantification`` packages.
Package layout:
* ``organoid_workspace``      — unified one-page entry point.
* ``segmentation_workspace``  — upload + preview + Cellpose 3D segmentation UI.
* ``analysis_ui``             — statistical analysis UI (Tutorials 2-5).
* ``multilevel_results``      — read-only view of exported multilevel results.
* ``demo_app``                — standalone vtk.js viewer demo.
"""

__version__ = "1.0.0"
