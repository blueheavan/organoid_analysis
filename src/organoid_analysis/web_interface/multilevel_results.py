"""Read-only Streamlit integration for exported multilevel 3D results."""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

RESULT_TABLES = {
    "Organoids": ("features/organoid_features.parquet", "organoid_features.parquet"),
    "Cells": ("features/cell_features.parquet", "cell_features.parquet"),
    "Nuclei": ("features/nucleus_features.parquet", "nucleus_features.parquet"),
    "Cell topology": ("features/cell_topology_edges.parquet", "cell_topology_edges.parquet"),
    "QC flags": ("qc/qc_flags.parquet", "qc_flags.parquet"),
}


def read_registered_organoid_labels(path: str | Path, expected_shape: tuple[int, ...],
                                   spacing_zyx_um: tuple[float, float, float]) -> np.ndarray:
    """Check uploaded mask axes and physical grid before joining hierarchies."""
    from organoid_analysis.microscopy_io.tiff_contract import (
        SPACING_ATOL_UM,
        SPACING_RTOL,
        read_tiff,
    )
    stack, spacing, _ = read_tiff(path)
    if stack.shape[0] != 1 or stack.shape[1:] != expected_shape:
        raise ValueError("Organoid mask must be a single channel on the cell/nucleus ZYX grid")
    if spacing is None:
        raise ValueError("Organoid mask spacing is missing; use a calibrated OME-TIFF or the analyze-3d CLI with explicit spacing")
    if not np.allclose(spacing, spacing_zyx_um, rtol=SPACING_RTOL, atol=SPACING_ATOL_UM):
        raise ValueError("Organoid mask spacing conflicts with the saved cell/nucleus grid")
    return np.asarray(stack[0])


def _session_spacing_zyx() -> tuple[float, float, float] | None:
    """Return saved Cellpose spacing when its current-session masks exist."""
    import streamlit as st

    config = st.session_state.get("segmentation_config")
    if config is None:
        return None
    xy = float(config.xy_spacing_um)
    return (xy * float(config.anisotropy), xy, xy)


def _render_current_session_runner() -> None:
    """Offer analysis from the masks already produced in the current UI session.

    Cellpose supplies cell/nucleus masks only.  This deliberately requires a
    separately segmented, registered organoid label mask rather than deriving
    one from cells or a projection and misrepresenting it as an organoid mask.
    """
    import streamlit as st

    from organoid_analysis.segmentation.cellpose_inference import (
        list_saved_results,
        restore_saved_result,
    )

    saved_runs = list_saved_results()
    if saved_runs:
        selected_run = st.selectbox(
            "Restore saved Cellpose segmentation",
            saved_runs,
            format_func=lambda item: f"{item.run_directory.name} — nuclei: {item.nuclei_count}, cells: {item.cell_count}",
            key="multilevel_saved_segmentation_run",
        )
        if st.button("Restore selected segmentation", key="multilevel_restore_segmentation"):
            try:
                restored = restore_saved_result(selected_run.run_directory)
                st.session_state["segmentation_result"] = restored.result
                st.session_state["segmentation_config"] = restored.config
                st.session_state["nuclei_masks"] = restored.nuclei_masks
                if restored.cell_masks is None:
                    st.session_state.pop("cell_masks", None)
                else:
                    st.session_state["cell_masks"] = restored.cell_masks
                st.success("Saved segmentation restored into this browser session.")
                st.rerun()
            except (OSError, RuntimeError, ValueError) as error:
                st.error(str(error))

    cell_labels = st.session_state.get("cell_masks")
    nucleus_labels = st.session_state.get("nuclei_masks")
    if cell_labels is None or nucleus_labels is None:
        st.info("Current session has no paired cell and nucleus masks. Run Cellpose with the optional cell/cytoplasm stack first, or use an exported Result directory below.")
        return
    st.subheader("Analyze current segmentation")
    st.caption("Cell and nucleus masks are reused from this session. Upload a separately segmented organoid label TIFF registered to the same ZYX grid.")
    organoid_upload = st.file_uploader("Organoid instance-label TIFF (required)", type=["tif", "tiff"], key="multilevel_organoid_upload")
    output = st.text_input("New/empty analysis output directory", key="multilevel_output_directory",
                           placeholder="/absolute/path/to/Result")
    metadata_columns = st.columns(3)
    sample_id = metadata_columns[0].text_input("sample_id (optional)", key="multilevel_sample_id")
    well_id = metadata_columns[1].text_input("well_id (optional)", key="multilevel_well_id")
    field_id = metadata_columns[2].text_input("field_id (optional)", key="multilevel_field_id")
    spacing = _session_spacing_zyx()
    if spacing:
        st.caption(f"Using saved Cellpose spacing: Z/Y/X = {spacing[0]:.4g}/{spacing[1]:.4g}/{spacing[2]:.4g} µm.")
    if st.button("Run multilevel 3D analysis", type="primary", key="multilevel_run_current"):
        if organoid_upload is None or not output:
            st.error("Provide both the registered organoid label TIFF and a new/empty output directory.")
            return
        if spacing is None:
            st.error("Saved voxel spacing is unavailable; use the analyze-3d CLI with explicit spacing instead.")
            return
        try:
            from organoid_analysis.microscopy_io.tiff_contract import sha256, source_code_hashes
            from organoid_analysis.result_export.measurement_tables import export_results
            from organoid_analysis.workflows.multilevel_measurement_workflow import (
                analyze_multilevel_3d,
            )

            with tempfile.NamedTemporaryFile(suffix=".tif") as handle:
                handle.write(organoid_upload.getbuffer())
                handle.flush()
                organoid_labels = read_registered_organoid_labels(handle.name, cell_labels.shape, spacing)
                organoid_sha256 = sha256(handle.name)
            metadata = {"sample_id": sample_id, "well_id": well_id, "field_id": field_id}
            with st.spinner("Measuring hierarchy, morphology, topology, spatial features, and QC…"):
                result = analyze_multilevel_3d(organoid_labels, cell_labels, nucleus_labels, spacing, metadata=metadata)
                import hashlib
                result.summary["provenance"] = {
                    "source_code_sha256": source_code_hashes(),
                    "organoid_tiff_sha256": organoid_sha256,
                    "session_masks": {name: {"sha256": hashlib.sha256(np.ascontiguousarray(mask).tobytes()).hexdigest(),
                                               "shape": list(mask.shape), "dtype": str(mask.dtype), "order": "C"}
                                      for name, mask in (("cell", cell_labels), ("nucleus", nucleus_labels))},
                }
                export_results(output, organoids=result.organoid_features, cells=result.cell_features,
                               nuclei=result.nucleus_features, edges=result.cell_topology_edges,
                               qc_flags=result.qc_flags, summary=result.summary,
                               label_masks={"organoid": organoid_labels, "cell": cell_labels, "nucleus": nucleus_labels},
                               spacing_zyx_um=spacing)
            st.session_state["multilevel_result_directory"] = str(Path(output).expanduser().resolve())
            st.success("Multilevel analysis completed. Result tables are loaded below.")
        except (OSError, RuntimeError, ValueError) as error:
            st.error(str(error))


def load_multilevel_result_tables(result_directory: str | Path) -> tuple[dict, dict[str, pd.DataFrame]]:
    """Load a completed `analyze-3d` result directory without recomputation."""
    root = Path(result_directory).expanduser().resolve()
    if (root / "RUN_INCOMPLETE.txt").exists():
        raise ValueError("This result export is incomplete; do not interpret partial outputs")
    summary_path = root / "summary" / "analysis_summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError("Select a Result directory containing summary/analysis_summary.json")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    tables = {}
    for title, (relative, _) in RESULT_TABLES.items():
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(f"Missing expected analysis output: {path.relative_to(root)}")
        tables[title] = pd.read_parquet(path)
    return summary, tables


def _load_exported_labels(root: Path) -> tuple[dict[str, np.ndarray], tuple[float, float, float]] | None:
    """Load optional CLI-copied label masks for QC visualization."""
    from organoid_analysis.microscopy_io.tiff_contract import read_tiff

    masks: dict[str, np.ndarray] = {}
    spacings = []
    for level in ("organoid", "cell", "nucleus"):
        path = root / "masks" / f"{level}_labels.ome.tif"
        if not path.is_file():
            return None
        stack, spacing, _ = read_tiff(path)
        if stack.shape[0] != 1 or spacing is None:
            return None
        masks[level] = stack[0]
        spacings.append(spacing)
    if not all(np.allclose(spacings[0], value) for value in spacings[1:]):
        raise ValueError("Exported label-mask spacing is inconsistent")
    return masks, tuple(float(value) for value in spacings[0])


def _render_label_slice(level: str, labels: np.ndarray, table: pd.DataFrame) -> None:
    """Render selectable Mask/Boundary/ID/QC overlays for a true Z slice."""
    import streamlit as st
    from skimage.color import label2rgb
    from skimage.segmentation import find_boundaries

    id_column = f"{level}_id"
    z = st.slider("Z slice", 0, labels.shape[0] - 1, labels.shape[0] // 2, key=f"multilevel_z_{level}")
    mode = st.selectbox("2D QC layer", ["Mask", "Boundary", "Object ID", "QC flags", "Raw"], key=f"multilevel_mode_{level}")
    visible_ids = np.unique(labels)
    visible_ids = visible_ids[visible_ids > 0]
    selected = st.selectbox("Selected object ID", visible_ids.tolist() or [0], key=f"multilevel_object_{level}")
    slice_labels = labels[z]
    if mode == "Raw":
        st.info("Raw intensity was not exported with this label-only analysis. Use the existing Upload & preview tab for raw Z-stack review.")
        image = np.zeros((*slice_labels.shape, 3), dtype=np.uint8)
    elif mode == "Boundary":
        image = np.zeros((*slice_labels.shape, 3), dtype=np.uint8)
        image[find_boundaries(slice_labels, mode="outer")] = (255, 255, 255)
    elif mode == "QC flags":
        flagged = set(table.loc[table.qc_status.eq("flagged"), id_column].astype(int)) if "qc_status" in table else set()
        image = np.zeros((*slice_labels.shape, 3), dtype=np.uint8)
        image[(slice_labels > 0) & np.isin(slice_labels, list(flagged))] = (230, 65, 65)
        image[(slice_labels > 0) & ~np.isin(slice_labels, list(flagged))] = (75, 180, 90)
    else:
        image = (label2rgb(slice_labels, bg_label=0) * 255).astype(np.uint8)
        if mode == "Object ID":
            image[slice_labels == selected] = (255, 255, 0)
    st.image(image, caption=f"{level.title()} Z={z}; selected object={selected}", width="stretch")
    row = table.loc[table[id_column].eq(selected)]
    if len(row):
        st.dataframe(row, use_container_width=True, hide_index=True)


def _render_multilevel_geometry(root: Path, tables: dict[str, pd.DataFrame]) -> None:
    """Reuse the vtk.js viewer, one independently toggleable real 3D level at a time."""
    import streamlit as st
    loaded = _load_exported_labels(root)
    if loaded is None:
        st.info("Masks were not copied into this Result directory; rerun analyze-3d to enable label QC geometry views.")
        return
    masks, spacing_zyx = loaded
    st.subheader("Multilevel mask QC")
    level = st.selectbox("2D object level", ["organoid", "cell", "nucleus"], key="multilevel_qc_level")
    _render_label_slice(level, masks[level], tables[{"organoid": "Organoids", "cell": "Cells", "nucleus": "Nuclei"}[level]])
    st.subheader("True 3D object geometry")
    st.caption("Use visibility toggles below. Every displayed surface comes from its 3D instance mask, never a MIP or density projection.")
    from organoid_analysis.visualization.volume_viewer import ChannelConfig, st_volume_viewer
    spacing_xyz = (spacing_zyx[2], spacing_zyx[1], spacing_zyx[0])
    for level, label in (("organoid", "Organoids"), ("cell", "Cells"), ("nucleus", "Nuclei")):
        if st.checkbox(f"Show {label}", value=level == "organoid", key=f"multilevel_show_{level}"):
            st_volume_viewer([masks[level].astype(np.float32)], spacing_um=spacing_xyz,
                             channels=[ChannelConfig(lut="gray", opacity=0.18, name=label)],
                             mask_overlay=masks[level], mask_overlay_alpha=0.65, height=440)


def render_multilevel_results_section() -> None:
    """Add a compact, filterable result reader to the existing unified app."""
    import streamlit as st

    st.subheader("3D Analysis results")
    st.caption("Read a completed `analysis analyze-3d` Result directory. This view never changes masks or measurements.")
    _render_current_session_runner()
    st.divider()
    st.subheader("Open exported Result")
    result_dir = st.text_input("Result directory", key="multilevel_result_directory",
                               placeholder="/absolute/path/to/Result")
    if not result_dir:
        st.info("Run the headless analyze-3d command, then enter its Result directory here.")
        return
    root = Path(result_dir).expanduser().resolve()
    try:
        summary, tables = load_multilevel_result_tables(result_dir)
    except (OSError, ValueError) as error:
        st.error(str(error))
        return
    columns = st.columns(4)
    columns[0].metric("Organoids", summary.get("organoid_count", 0))
    columns[1].metric("Cells", summary.get("cell_count", 0))
    columns[2].metric("Nuclei", summary.get("nucleus_count", 0))
    columns[3].metric("QC flags", len(tables["QC flags"]))
    st.caption(f"Input shape ZYX: {summary.get('input_shape_zyx')} | runtime: {summary.get('runtime_seconds', 0):.2f} s")
    selected = st.selectbox("Feature table", list(RESULT_TABLES), key="multilevel_table")
    table = tables[selected]
    filter_text = st.text_input("Filter rows (case-insensitive text)", key="multilevel_filter")
    if filter_text and not table.empty:
        mask = table.astype(str).apply(lambda column: column.str.contains(filter_text, case=False, na=False)).any(axis=1)
        table = table.loc[mask]
    st.dataframe(table, use_container_width=True, hide_index=True)
    _, filename = RESULT_TABLES[selected]
    st.download_button(f"Download {filename}", data=table.to_parquet(index=False), file_name=filename,
                       mime="application/vnd.apache.parquet")
    _render_multilevel_geometry(root, tables)
