"""Integrated 3D organoid pipeline: Upload -> Preview -> Segment -> Analyze.

Closes the loop inside the ``ui.vtk_viewer`` 3D viewer:
  1. Upload a nuclei (and optional cytoplasm) Z-stack TIFF + set voxel spacing.
  2. Preview the intensity volume with the browser GPU volume renderer.
  3. Run Cellpose 3D segmentation (reusing ``segmentation.cellpose``) in a background
     thread with a live progress bar.
  4. Show the segmentation mask overlaid on the intensity volume (2D per-label
     contours + 3D surface) and present nuclei/cell counts.
  5. Extract per-object features (``mask_features``) and run on the derived
     feature table streaming distributions / statistical summaries.

The heavy Cellpose model is cached per process via ``st.cache_resource`` exactly
like the original ``ui/upload.py``, so only the first run downloads weights.

Run:
    streamlit run src/ui/vtk_viewer/integrated_app.py
"""

from __future__ import annotations

import os

# torch can abort at import if two OpenMP runtimes are linked; allow it so the
# app (and Cellpose segmentation) can start from any shell / launcher.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import queue
import hashlib
import sys
import tempfile
import threading
import time
from pathlib import Path

import numpy as np
import streamlit as st
from streamlit.components.v1 import html as _st_html

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
# This module can also be launched directly. Keep sibling core packages before
# ``src/ui`` so ``analysis.features`` resolves to the package, not ui/analysis.py.
while str(_PROJECT_ROOT) in sys.path:
    sys.path.remove(str(_PROJECT_ROOT))
sys.path.insert(0, str(_PROJECT_ROOT))

from segmentation.cellpose import (  # noqa: E402
    SegmentationConfig,
    create_model,
    get_accelerator,
    read_stack_multichannel,
    save_result,
    segment_stacks,
    validate_stacks,
)
from segmentation.auto_config import estimate_diameter_from_stack  # noqa: E402
from ui.vtk_viewer import ChannelConfig, st_volume_viewer  # noqa: E402
from ui.vtk_viewer.mask_features import (  # noqa: E402
    count_mask_objects,
    extract_mask_features,
    summarize_features,
)


@st.cache_resource(show_spinner="Loading Cellpose model weights...")
def get_model(model_type: str = ""):
    return create_model(model_type)


def _set_tab_busy(busy: bool, slot=None) -> None:
    """Blink the browser tab's favicon + title while a long-running analysis
    (segmentation) is in progress, so it's visible from another tab/window.

    Injected JS reaches through the component iframe to ``window.parent`` --
    same-origin, so this is permitted -- and stores its interval/original
    state on ``window.parent`` itself so start/stop calls across separate
    Streamlit component renders can find each other. Pass the *same*
    ``st.empty()`` placeholder for both the start and stop call (``slot``) --
    reusing one DOM slot instead of appending a fresh iframe wherever the
    script happens to be when each call runs is what makes the stop signal
    reliable, since some browsers defer loading an iframe that ends up
    scrolled out of view before it's ever rendered.
    """
    if busy:
        script = """
        (function() {
          try {
            var doc = window.parent.document, w = window.parent;
            if (w.__organoidBusyTimer) { return; }
            var icon = doc.querySelector('link[rel~="icon"]');
            if (!icon) { icon = doc.createElement('link'); icon.rel = 'icon'; doc.head.appendChild(icon); }
            w.__organoidOrigIcon = icon.getAttribute('href');
            w.__organoidOrigTitle = doc.title;
            var bright = 'data:image/svg+xml,' + encodeURIComponent(
              '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 32 32">'
              + '<circle cx="16" cy="16" r="13" fill="red"/></svg>');
            var dim = 'data:image/svg+xml,' + encodeURIComponent(
              '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 32 32">'
              + '<circle cx="16" cy="16" r="6" fill="red"/></svg>');
            var on = false;
            w.__organoidBusyTimer = setInterval(function() {
              on = !on;
              icon.setAttribute('href', on ? bright : dim);
              doc.title = (on ? '\\uD83D\\uDD34 ' : '\\u26AA ') + 'Running... ' + w.__organoidOrigTitle;
            }, 500);
            // Safety net: clear on its own if the explicit "stop" signal is
            // ever missed (e.g. this tab lost focus before it could render),
            // so the indicator can never blink forever.
            if (w.__organoidBusyTimeout) { clearTimeout(w.__organoidBusyTimeout); }
            w.__organoidBusyTimeout = setTimeout(function() {
              if (w.__organoidBusyTimer) { clearInterval(w.__organoidBusyTimer); w.__organoidBusyTimer = null; }
              if (icon && w.__organoidOrigIcon != null) { icon.setAttribute('href', w.__organoidOrigIcon); }
              if (w.__organoidOrigTitle != null) { doc.title = w.__organoidOrigTitle; }
            }, 30 * 60 * 1000);
          } catch (e) {}
        })();
        """
    else:
        script = """
        (function() {
          try {
            var doc = window.parent.document, w = window.parent;
            if (w.__organoidBusyTimer) { clearInterval(w.__organoidBusyTimer); w.__organoidBusyTimer = null; }
            if (w.__organoidBusyTimeout) { clearTimeout(w.__organoidBusyTimeout); w.__organoidBusyTimeout = null; }
            var icon = doc.querySelector('link[rel~="icon"]');
            if (icon && w.__organoidOrigIcon != null) { icon.setAttribute('href', w.__organoidOrigIcon); }
            if (w.__organoidOrigTitle != null) { doc.title = w.__organoidOrigTitle; }
          } catch (e) {}
        })();
        """
    html = f"<script>{script}</script>"
    if slot is not None:
        with slot.container():
            _st_html(html, height=0, width=0)
    else:
        _st_html(html, height=0, width=0)


def config_spacing(config: SegmentationConfig) -> tuple[float, float, float]:
    """Physical spacing (x, y, z) in µm from XY pixel size and anisotropy."""
    sy = sx = config.xy_spacing_um
    return (sx, sy, sx * config.anisotropy)


# Streamlit's websocket message-size limit defaults to 200MB; the viewer HTML
# also carries mesh/surface data alongside the raw arrays, so budget well
# under that for the base64-encoded intensity + mask volumes alone.
_VIEWER_PAYLOAD_BUDGET_BYTES = 130_000_000


def _fit_viewer_payload_budget(
    volumes: list[np.ndarray],
    mask: np.ndarray | None,
    spacing: tuple[float, float, float],
) -> tuple[list[np.ndarray], np.ndarray | None, tuple[float, float, float]]:
    """Downsample XY (never Z) so the vtk.js payload fits Streamlit's message
    size limit on large uploads. Every array the viewer receives becomes
    4 bytes/voxel (float32 intensity or uint32 labels) plus ~1.34x for base64,
    so budget on that. Only the returned *preview* copies are downsampled --
    callers keep the full-resolution mask for feature extraction and
    downloads; label identity is preserved exactly via nearest-neighbor.
    """
    if not volumes:
        return volumes, mask, spacing
    z, y, x = volumes[0].shape
    n_arrays = len(volumes) + (1 if mask is not None else 0)
    estimated_bytes = z * y * x * 4 * n_arrays * 1.34
    if estimated_bytes <= _VIEWER_PAYLOAD_BUDGET_BYTES:
        return volumes, mask, spacing
    factor = max(0.1, min(1.0, (_VIEWER_PAYLOAD_BUDGET_BYTES / estimated_bytes) ** 0.5))
    from scipy import ndimage

    zoom_xy = (1.0, factor, factor)
    volumes = [
        ndimage.zoom(v, zoom_xy, order=1).astype(np.float32, copy=False) for v in volumes
    ]
    if mask is not None:
        mask = ndimage.zoom(mask, zoom_xy, order=0).astype(mask.dtype, copy=False)
    sx, sy, sz = spacing
    spacing = (sx / factor, sy / factor, sz)
    return volumes, mask, spacing


def _read_upload(upload, path: Path, label: str, role: str) -> tuple[np.ndarray, str, int, object]:
    """Decode an upload once per content digest and select its active channel."""
    digest = hashlib.sha256(upload.getbuffer()).hexdigest()
    digest_key = f"{role}_upload_digest"
    stack_key = f"{role}_decoded_zstack"
    if st.session_state.get(digest_key) != digest or stack_key not in st.session_state:
        path.write_bytes(upload.getbuffer())
        st.session_state[stack_key] = read_stack_multichannel(path)
        st.session_state[digest_key] = digest
    zstack = st.session_state[stack_key]
    if not zstack.is_multichannel:
        return zstack.volume, digest, 0, zstack
    channel = int(st.selectbox(
        f"{label} channel", options=list(range(zstack.channels)), key=f"{role}_channel"
    ))
    return zstack.volume[channel], digest, channel, zstack


def render_preview_tab(config: SegmentationConfig) -> None:
    st.subheader("Upload stacks")
    left, right = st.columns(2)
    with left:
        nuclei_upload = st.file_uploader(
            "Nuclei TIFF stack (required)", type=["tif", "tiff"]
        )
    with right:
        cells_upload = st.file_uploader(
            "Cell/cytoplasm TIFF stack (optional)", type=["tif", "tiff"]
        )

    if nuclei_upload is None:
        st.info("Upload a nuclei Z-stack TIFF to preview and segment.")
        return

    with tempfile.TemporaryDirectory() as temporary_directory:
        temporary_path = Path(temporary_directory)
        try:
            nuclei, nuclei_digest, nuclei_channel, nuclei_zstack = _read_upload(
                nuclei_upload, temporary_path / "nuclei_input.tif", "Nuclei", "nuclei"
            )
            cells = None
            cells_digest = None
            cells_channel = None
            if cells_upload is not None:
                cells, cells_digest, cells_channel, _ = _read_upload(
                    cells_upload, temporary_path / "cells_input.tif", "Cell/cytoplasm", "cells"
                )
                validate_stacks(nuclei, cells)
        except (OSError, ValueError) as error:
            st.error(str(error))
            return

        input_key = (
            nuclei_digest,
            nuclei_channel,
            cells_digest,
            cells_channel,
        )
        if st.session_state.get("input_key") != input_key:
            st.session_state["input_key"] = input_key
            for key in ("nuclei_masks", "cell_masks", "segmentation_result",
                        "segmentation_config", "orig_nuclei", "orig_cells",
                        "features", "summary", "feature_key"):
                st.session_state.pop(key, None)

        # Store originals so the results / analysis tabs can show the
        # intensity volume alongside the mask overlay.
        st.session_state["orig_nuclei"] = nuclei
        if cells is not None:
            st.session_state["orig_cells"] = cells
        else:
            st.session_state.pop("orig_cells", None)

        spacing = config_spacing(config)
        if cells is not None:
            vols = [nuclei, cells]
            channels = [
                ChannelConfig(lut="cyan", opacity=0.9, name="Nuclei"),
                ChannelConfig(lut="magenta", opacity=0.6, name="Cytoplasm"),
            ]
        else:
            vols = [nuclei]
            channels = [ChannelConfig(lut="cyan", opacity=0.9, name="Nuclei")]

        preview_vols, _, preview_spacing = _fit_viewer_payload_budget(vols, None, spacing)
        st_volume_viewer(preview_vols, spacing_um=preview_spacing, channels=channels, height=620)
        st.caption(
            "Intensity preview. Set segmentation parameters in the sidebar and run below."
        )

        auto_col = st.columns([1, 3])
        with auto_col[0]:
            auto_detect = st.button("Auto-detect parameters", type="secondary")
        with auto_col[1]:
            st.caption(
                "Estimates nuclei diameter from a fast 2D pass and reads Z/XY "
                "anisotropy from the TIFF metadata."
            )

        if auto_detect:
            auto_col0 = st.columns(2)
            with auto_col0[0]:
                with st.spinner("Estimating nuclei diameter (fast 2D pass)…"):
                    est = estimate_diameter_from_stack(
                        get_model(config.model_type), nuclei
                    )
            with auto_col0[1]:
                anisotropy = nuclei_zstack.spacing.anisotropy
            suggestion = {"nuclei_diameter": est, "cell_diameter": est}
            if anisotropy is not None:
                suggestion["anisotropy"] = anisotropy
            st.session_state["auto_suggest"] = suggestion
            msg = "Auto-detected: "
            if est:
                msg += f"nuclei diameter ≈ **{est:.1f} px**"
            else:
                msg += "could not estimate diameter (no objects found)"
            if anisotropy is not None:
                msg += f", anisotropy ≈ **{anisotropy:.2f}** (from metadata)"
            else:
                msg += ". No Z/XY spacing in metadata — please set anisotropy manually."
            st.success(msg)
            if suggestion:
                st.rerun()

        if st.button("Run 3D segmentation", type="primary"):
            _run_segmentation(nuclei, cells, config)


def _run_segmentation(nuclei, cells, config: SegmentationConfig) -> None:
    task = {"start": time.monotonic(), "frac": 0.0, "estimated": None}
    result_queue: queue.Queue = queue.Queue()

    def report(fraction: float) -> None:
        now = time.monotonic()
        task["frac"] = fraction
        if fraction >= 0.5 and task["estimated"] is None:
            task["estimated"] = now - task["start"]

    def worker() -> None:
        try:
            masks = segment_stacks(
                get_model(config.model_type), nuclei, cells, config, on_progress=report
            )
            saved = save_result(masks[0], masks[1], config)
            result_queue.put(("ok", saved, masks))
        except Exception as error:  # noqa: BLE001 - propagate worker failures to the UI
            result_queue.put(("error", error, None))

    progress_bar = st.progress(0.0, text="Starting Cellpose…")
    tab_indicator = st.empty()
    thread = threading.Thread(target=worker, daemon=True)
    thread.start()
    _set_tab_busy(True, tab_indicator)
    try:
        while thread.is_alive():
            elapsed = time.monotonic() - task["start"]
            fraction = min(task["frac"], 1.0)
            message = f"{int(fraction * 100)}% | Elapsed: {int(elapsed)}s"
            if task["estimated"] is not None:
                message += f" | ETA: ~{int(max(task['estimated'] - elapsed, 0))}s"
            progress_bar.progress(fraction, text=message)
            time.sleep(0.25)
    finally:
        # Leave the stop script's (now invisible, zero-size) iframe in place
        # rather than clearing the placeholder right after -- removing it
        # immediately risks the browser dropping the element before its
        # script has a chance to actually run.
        _set_tab_busy(False, tab_indicator)
    progress_bar.empty()

    try:
        status, payload, masks = result_queue.get_nowait()
    except queue.Empty:
        st.error("Segmentation worker stopped without returning a result. Restart the app before retrying.")
        return
    if status == "error":
        st.error(f"Segmentation failed: {payload}")
        return
    result, nuclei_masks, cell_masks = payload, masks[0], masks[1]
    st.session_state["segmentation_result"] = result
    st.session_state["segmentation_config"] = config
    st.session_state["nuclei_masks"] = nuclei_masks
    if cell_masks is not None:
        st.session_state["cell_masks"] = cell_masks
    else:
        st.session_state.pop("cell_masks", None)
    for key in ("features", "summary", "feature_key"):
        st.session_state.pop(key, None)
    msg = f"Finished: {result.nuclei_count} nuclei detected"
    if cell_masks is not None:
        msg += f" and {result.cell_count} cells"
    st.success(msg + ".")
    # A transient browser notification makes completion visible even when the
    # user is watching another part of the long-running segmentation view.
    st.toast("Segmentation Done", icon="✅")


def render_results_tab(config: SegmentationConfig) -> None:
    nuclei_masks = st.session_state.get("nuclei_masks")
    if nuclei_masks is None:
        st.info("Run a segmentation first; results appear here.")
        return
    result = st.session_state.get("segmentation_result")
    cell_masks = st.session_state.get("cell_masks")
    result_config = st.session_state.get("segmentation_config", config)
    spacing = config_spacing(result_config)
    if result_config != config:
        st.info("Sidebar settings changed after this run. Measurements below use the saved segmentation settings.")

    st.subheader("Mask overlay + per-object features")
    metric_cols = st.columns(2)
    with metric_cols[0]:
        st.metric("Nuclei detected", result.nuclei_count if result is not None else count_mask_objects(nuclei_masks))
    if cell_masks is not None:
        with metric_cols[1]:
            st.metric("Cells detected", result.cell_count if result is not None else count_mask_objects(cell_masks))

    # Use the original intensity image (if stored) so the overlay is meaningful;
    # fall back to the mask values when originals are unavailable.
    orig_nuclei = st.session_state.get("orig_nuclei")
    orig_cells = st.session_state.get("orig_cells")
    if orig_nuclei is not None:
        if orig_cells is not None:
            vols = [
                orig_nuclei,
                orig_cells,
            ]
            chs = [
                ChannelConfig(lut="cyan", opacity=0.9, name="Nuclei"),
                ChannelConfig(lut="magenta", opacity=0.6, name="Cytoplasm"),
            ]
        else:
            vols = [orig_nuclei]
            chs = [ChannelConfig(lut="cyan", opacity=0.9, name="Nuclei")]
    else:
        vols = [nuclei_masks.astype(np.float32, copy=False)]
        chs = [ChannelConfig(lut="gray", opacity=0.5, name="Mask")]

    preview_vols, preview_mask, preview_spacing = _fit_viewer_payload_budget(vols, nuclei_masks, spacing)
    st_volume_viewer(
        preview_vols,
        spacing_um=preview_spacing,
        channels=chs,
        render_mode="volume",
        mask_overlay=preview_mask,
        mask_overlay_alpha=0.6,
        height=620,
    )
    if preview_vols[0].shape != nuclei_masks.shape:
        st.caption(
            f"3D preview downsampled to {preview_vols[0].shape} for browser display "
            f"(full resolution {nuclei_masks.shape} used for all measurements below)."
        )
    st.caption("Colored contours = object boundaries; translucent shell = 3D surface.")
    with st.expander("How to review this segmentation", expanded=False):
        st.markdown(
            "Review the volume in all three directions, not only the brightest projection. Look for "
            "one contour spanning two nuclei, one nucleus split into several labels, labels that disappear "
            "between Z slices, and partial objects at image borders. The object count is a model output, "
            "not a cell count validated against manual annotation. Re-run with a representative parameter "
            "set at full XY resolution before exporting measurements. For whole-organoid morphology, use "
            "a mask that outlines the organoid boundary; a nuclei mask measures nuclei, not organoid size."
        )

    with st.expander("Per-object features (measured from the mask)", expanded=True):
        feature_key = (
            str(result.run_directory) if result is not None else id(nuclei_masks),
            tuple(spacing),
        )
        if st.session_state.get("feature_key") != feature_key:
            st.session_state["features"] = extract_mask_features(nuclei_masks, spacing_um=spacing)
            st.session_state["feature_key"] = feature_key
        features = st.session_state["features"]
        st.session_state["features"] = features
        if features.empty:
            st.warning("No objects found in the mask.")
        else:
            st.dataframe(features.reset_index(), use_container_width=True)

    summary = summarize_features(features)
    st.session_state["summary"] = summary
    mcols = st.columns(4)
    mcols[0].metric("Objects", summary.n_objects)
    mcols[1].metric("Total volume (µm³)", f"{summary.total_volume_um3:,.1f}")
    mcols[2].metric("Mean volume (µm³)", f"{summary.mean_volume_um3:,.1f}")
    mcols[3].metric("Mean sphericity", f"{summary.mean_sphericity:.3f}")

    if result is not None and result.archive_path.exists():
        st.download_button(
            "Download masks + summary",
            data=result.archive_path.read_bytes(),
            file_name=result.archive_path.name,
            mime="application/zip",
        )


def render_analysis_tab() -> None:
    st.subheader("Feature statistics")
    features = st.session_state.get("features")
    if features is None or features.empty:
        st.info("Run a segmentation first to populate the per-object feature table.")
        return
    summary = st.session_state.get("summary")

    with st.expander("What these feature statistics mean", expanded=True):
        st.markdown(
            "This tab summarizes labels from the current uploaded image only. It is descriptive and does "
            "not compare treatments or estimate viability. **Volume** is the hole-filled mask volume in "
            "µm³; **surface area** is a marching-cubes estimate in µm²; **sphericity** compares the mask "
            "with an equal-volume sphere; **solidity** compares the labeled object with its convex hull; "
            "and principal axes describe a moment-equivalent ellipsoid. Values depend on the entered voxel "
            "spacing and mask quality. Border-clipped, merged, or split labels should be excluded upstream "
            "rather than interpreted as unusual biology."
        )

    numeric = features.select_dtypes(include=[np.number])
    st.dataframe(numeric.describe(), use_container_width=True)

    st.markdown("**Volume distribution**")
    stat_col, chart_col = st.columns([1, 2])
    with stat_col:
        if summary is not None:
            st.markdown(
                f"- **Median volume**: {summary.median_volume_um3:,.2f} µm³\n"
                f"- **Min / max**: {summary.min_volume_um3:,.2f} / {summary.max_volume_um3:,.2f} µm³\n"
                f"- **Mean sphericity**: {summary.mean_sphericity:.3f}\n"
                f"- **Mean solidity**: {summary.mean_solidity:.3f}"
            )
    with chart_col:
        st.bar_chart(features["volume_um3"].sort_values().reset_index(drop=True))

    st.markdown("**Sphericity vs solidity (segmentation quality)**")
    st.scatter_chart(
        features[["sphericity", "solidity"]],
        x="sphericity",
        y="solidity",
        use_container_width=True,
    )

    st.caption(
        "For group statistics, export the feature table and use a design that preserves well, biological "
        "replicate and batch identifiers. The canonical analysis CLI performs this hierarchy-aware summary; "
        "the Excel tutorial page is exploratory and must not treat all objects as independent replicates."
    )


def render_sidebar() -> SegmentationConfig:
    model_labels = {
        "cpsam_v2": "cpsam_v2 — CellposeSAM (SAM-ViTL)",
        "cpdino": "cpdino — CellposeDINO (DINOv3-ViTL)",
        "cpdino-vitb": "cpdino-vitb — CellposeDINO (ViT-B, fastest)",
        "cpsam": "cpsam — original CellposeSAM",
    }
    # Auto-detected suggestions (filled by the "Auto-detect" button in the
    # preview tab) become the defaults for the number inputs.
    suggestions = st.session_state.get("auto_suggest") or {}

    def _default_for(key: str, fallback: float) -> float:
        value = suggestions.get(key)
        if isinstance(value, (int, float)) and value > 0:
            return float(value)
        return fallback

    # Single source of truth for every fallback below, so a future edit to
    # SegmentationConfig's dataclass defaults cannot silently desynchronize
    # from what this sidebar displays/uses (docs/PARAMETERS.md, "Cellpose 3D
    # segmentation parameters").
    backend_defaults = SegmentationConfig()

    with st.sidebar:
        st.header("Segmentation settings")
        model_type = st.selectbox(
            "Model",
            list(model_labels),
            index=list(model_labels).index(backend_defaults.model_type),
            format_func=lambda m: model_labels[m],
        )
        if suggestions:
            st.success(
                "Auto-detected: "
                f"diameter ≈ {suggestions['nuclei_diameter']:.1f} px"
                + (
                    f", anisotropy ≈ {suggestions['anisotropy']:.2f}"
                    if suggestions.get("anisotropy")
                    else ""
                )
            )
        return SegmentationConfig(
            model_type=model_type,
            nuclei_diameter=st.number_input(
                "Nuclei diameter (px)", 5.0, 200.0, _default_for("nuclei_diameter", backend_defaults.nuclei_diameter),
                help="Set manually, or press 'Auto-detect' after uploading to "
                "estimate it from the image.",
            ),
            cell_diameter=st.number_input(
                "Cell diameter (px)", 5.0, 300.0, _default_for("cell_diameter", backend_defaults.cell_diameter)
            ),
            anisotropy=st.number_input(
                "Z / XY anisotropy", 0.1, 20.0, _default_for("anisotropy", backend_defaults.anisotropy),
                help="Auto-filled from TIFF metadata when available; otherwise "
                "set from your microscope Z step / XY pixel size.",
            ),
            xy_spacing_um=st.number_input(
                "XY pixel size (µm)", 0.001, 100.0, backend_defaults.xy_spacing_um,
                help="Physical XY voxel spacing. Set to your microscope's pixel "
                "size so volumes/areas are in real units.",
            ),
            nuclei_flow_threshold=st.slider("Nuclei flow threshold", 0.0, 1.0, backend_defaults.nuclei_flow_threshold),
            cell_flow_threshold=st.slider("Cell flow threshold", 0.0, 1.0, backend_defaults.cell_flow_threshold),
            nuclei_cellprob_threshold=st.slider(
                "Nuclei cellprob threshold", -6.0, 6.0, backend_defaults.nuclei_cellprob_threshold
            ),
            cell_cellprob_threshold=st.slider("Cell cellprob threshold", -6.0, 6.0, backend_defaults.cell_cellprob_threshold),
            flow3d_smooth=st.slider("3D flow smoothing", 0.0, 5.0, backend_defaults.flow3d_smooth),
            batch_size=st.select_slider("Batch size", options=[1, 2, 4, 8, 16], value=backend_defaults.batch_size),
            xy_downsample={
                "Full resolution": 1.0,
                "1/2 resolution": 0.5,
                "1/4 resolution": 0.25,
            }[st.radio(
                "XY speed / resolution",
                ["Full resolution", "1/2 resolution", "1/4 resolution"],
                index=0,
            )],
        )
        with st.expander("Organoid tips", expanded=False):
            st.markdown(
                "**Organoid images** are large multicellular spheroids. "
                "Segmentation splits each 3D stack into individual nuclei (or "
                "cells). If the result does not look like organoids, press "
                "**Auto-detect** after uploading, then check:\n\n"
                "- **Nuclei diameter** — one nucleus width in pixels.\n"
                "- **Z / XY anisotropy** — Z step ÷ XY pixel size (auto-filled "
                "from metadata when present).\n"
                "- **Resolution** — down-sample for speed, then re-run at full "
                "resolution for final analysis.\n"
                "- **Model** — `cpdino-vitb` is fastest; `cpsam_v2` is most "
                "accurate for densely packed structures.\n"
                "- **QC** — inspect every representative condition for missed, merged and "
                "Z-fragmented objects. Changing a model or threshold changes the measured population.\n"
                "- **Scope** — nuclei/cell masks support nuclei/cell measurements. They do "
                "not automatically define a whole-organoid outer boundary or a viability state."
            )


def main() -> None:
    st.set_page_config(page_title="3D Organoid Pipeline", layout="wide")
    st.title("3D Organoid Pipeline: Upload → Preview → Segment → Analyze")
    try:
        accelerator = get_accelerator().upper()
    except Exception:  # torch/Cellpose not importable -> preview-only mode
        accelerator = "unavailable (preview + results only)"
    st.caption(
        f"Compute device: `{accelerator}`. Browser GPU rendering via vtk.js."
    )
    config = render_sidebar()

    tab_preview, tab_results, tab_analysis = st.tabs(
        ["Upload & preview", "Segmentation results", "Analysis"]
    )
    with tab_preview:
        render_preview_tab(config)
    with tab_results:
        render_results_tab(config)
    with tab_analysis:
        render_analysis_tab()


if __name__ == "__main__":
    main()
