"""Integrated 3D organoid pipeline: Upload -> Preview -> Segment -> Analyze.

Closes the loop inside the ``organoid_analysis.visualization.volume_viewer`` 3D viewer:
  1. Upload a nuclei (and optional cytoplasm) Z-stack TIFF + set voxel spacing.
  2. Preview the intensity volume with the browser GPU volume renderer.
  3. Run Cellpose 3D segmentation (reusing ``organoid_analysis.segmentation.cellpose_inference``)
     in a background thread with a live progress bar.
  4. Show the segmentation mask overlaid on the intensity volume (2D per-label
     contours + 3D surface) and present nuclei/cell counts.
  5. Extract per-object features (``mask_features``) and run on the derived
     feature table streaming distributions / statistical summaries.

The heavy Cellpose model is cached per process via ``st.cache_resource``, so
only the first run downloads weights.

Run:
    streamlit run src/organoid_analysis/web_interface/segmentation_workspace.py
"""

from __future__ import annotations

import os

# torch can abort at import if two OpenMP runtimes are linked; allow it so the
# app (and Cellpose segmentation) can start from any shell / launcher.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import hashlib
import json
import queue
import tempfile
import threading
import time
from dataclasses import asdict
from pathlib import Path
from typing import TypedDict

import numpy as np
import streamlit as st
from streamlit.components.v1 import html as _st_html
from streamlit.delta_generator import DeltaGenerator
from streamlit.runtime.uploaded_file_manager import UploadedFile

from organoid_analysis.microscopy_io import (  # noqa: E402
    GRID_CONFLICT,
    GRID_PARTIAL,
    ZStack,
    compare_registered_grid,
    isotropic_xy_size_um,
    resolve_spacing_source,
)
from organoid_analysis.microscopy_io.resampling import (
    isotropic_xy_scale,
    resampled_spacing_xyz,
    xy_downsample_shape,
)
from organoid_analysis.quantification.mask_features import (  # noqa: E402
    count_mask_objects,
    extract_mask_features,
    summarize_features,
)
from organoid_analysis.result_export.mask_feature_bundle import build_mask_feature_bundle
from organoid_analysis.segmentation.cellpose_inference import (  # noqa: E402
    SegmentationConfig,
    config_spacing,
    create_model,
    get_accelerator,
    read_stack_multichannel,
    save_result,
    segment_stacks,
    validate_stacks,
)
from organoid_analysis.segmentation.parameter_estimation import (
    estimate_diameter_from_stack,  # noqa: E402
)
from organoid_analysis.visualization.volume_viewer import (  # noqa: E402
    ChannelConfig,
    st_volume_viewer,
)


@st.cache_resource(show_spinner="Loading Cellpose model weights...")
def get_model(model_type: str = ""):
    return create_model(model_type)


def _set_tab_busy(busy: bool, slot: DeltaGenerator | None = None) -> None:
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
    source_shape = volumes[0].shape
    if len(source_shape) != 3 or any(v.shape != source_shape for v in volumes) or (mask is not None and mask.shape != source_shape):
        raise ValueError("All preview channels and labels must share one ZYX grid")
    z, y, x = source_shape
    n_arrays = len(volumes) + (1 if mask is not None else 0)
    estimated_bytes = z * y * x * 4 * n_arrays * 1.34
    if estimated_bytes <= _VIEWER_PAYLOAD_BUDGET_BYTES:
        return volumes, mask, spacing
    max_xy_pixels = int(_VIEWER_PAYLOAD_BUDGET_BYTES / (z * 4 * n_arrays * 1.34))
    if max_xy_pixels < 4 or min(y, x) < 2:
        raise ValueError("Preview exceeds the browser budget even at two XY pixels; select a smaller Z-stack")
    factor = (max_xy_pixels / (y * x)) ** 0.5
    target_y = max(2, min(y, int(y * factor)))
    target_x = max(2, min(x, int(x * factor)))
    if target_y * target_x > max_xy_pixels:
        if target_y == 2:
            target_x = max_xy_pixels // 2
        else:
            target_y = max_xy_pixels // target_x
    from scipy import ndimage

    zoom_xy = (1.0, target_y / y, target_x / x)
    volumes = [
        ndimage.zoom(v, zoom_xy, order=1, output=np.float32, grid_mode=False) for v in volumes
    ]
    if mask is not None:
        mask = ndimage.zoom(mask, zoom_xy, order=0, grid_mode=False)
    spacing = resampled_spacing_xyz(source_shape, volumes[0].shape, spacing)
    return volumes, mask, spacing


def resolve_auto_suggestion(
    raw_suggestion: dict | None, current_identity: tuple[str, int] | None
) -> dict:
    """Return ``raw_suggestion`` only if it was produced from the currently
    loaded nuclei stack + channel; otherwise return an empty suggestion.

    ``raw_suggestion`` must carry the ``source_digest``/``source_channel`` of
    the nuclei stack it was computed from (see the "Auto-detect parameters"
    button in ``render_preview_tab``). Without this check, an auto-detected
    diameter/spacing/anisotropy computed for one uploaded stack would keep
    being applied as if it were measured from a later, different upload or a
    different channel of the same TIFF -- silently mislabeling the new run's
    physical spacing as "from metadata" when it is actually stale data from
    an unrelated file. A pure function so this identity check is testable
    without spinning up Streamlit.
    """
    if not raw_suggestion or current_identity is None:
        return {}
    source_identity = (raw_suggestion.get("source_digest"), raw_suggestion.get("source_channel"))
    if source_identity != current_identity:
        return {}
    return raw_suggestion


def _read_upload(
    upload: UploadedFile, path: Path, label: str, role: str
) -> tuple[np.ndarray, str, int, ZStack]:
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
    with st.expander("Image suitability for 3D measurement", expanded=False):
        st.markdown(
            "Use the original quantitative Z-stack, not an RGB screenshot, MIP, or a "
            "pre-contrast-adjusted display export. The selected nuclei channel must show "
            "nuclei through the full depth of the object, with no widespread saturation and "
            "with background and object signal visibly separable. Preserve OME/ImageJ voxel "
            "spacing metadata, or supply measured Z/Y/X spacing before interpreting µm, µm², "
            "or µm³ outputs.\n\n"
            "For membrane/cytoplasm segmentation, acquire the companion channel on the same "
            "physical grid and register it before upload. Do not infer 3D morphology from a "
            "single slice or projection. Border-truncated, fragmented, low-SNR, or strongly "
            "depth-attenuated objects require QC review and may be outside the validated use "
            "domain. A successful segmentation run is not biological validation."
        )
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
            cells_zstack = None
            if cells_upload is not None:
                cells, cells_digest, cells_channel, cells_zstack = _read_upload(
                    cells_upload, temporary_path / "cells_input.tif", "Cell/cytoplasm", "cells"
                )
                validate_stacks(nuclei, cells)
        except (OSError, ValueError) as error:
            st.error(str(error))
            return

        nuclei_identity = (nuclei_digest, nuclei_channel)
        if st.session_state.get("nuclei_identity") != nuclei_identity:
            st.session_state["nuclei_identity"] = nuclei_identity
            # A previously auto-detected diameter/spacing/anisotropy was
            # measured from a *different* nuclei stack or channel; it must
            # never be silently carried over into this one and mislabeled as
            # "from metadata" (see resolve_auto_suggestion, P1-1 audit finding).
            st.session_state.pop("auto_suggest", None)

        if cells is not None and cells_zstack is not None:
            grid_status = compare_registered_grid(nuclei_zstack.spacing, cells_zstack.spacing)
            if grid_status == GRID_CONFLICT:
                st.error(
                    "Nuclei and cell/cytoplasm stacks report different physical "
                    f"voxel spacing in their metadata (nuclei={nuclei_zstack.spacing}, "
                    f"cell={cells_zstack.spacing}). They cannot be treated as "
                    "registered channels of the same acquisition."
                )
                return
            if grid_status == GRID_PARTIAL:
                st.warning(
                    "Voxel spacing metadata is incomplete on at least one of the "
                    "nuclei/cell stacks, so agreement on the physical acquisition "
                    "grid cannot be confirmed from metadata."
                )
                confirmed = st.checkbox(
                    "I confirm the nuclei and cell/cytoplasm stacks were acquired "
                    "on the same physical grid (same field of view, voxel size, "
                    "and Z step).",
                    key=f"grid_confirm_{nuclei_digest}_{nuclei_channel}_{cells_digest}_{cells_channel}",
                )
                if not confirmed:
                    st.info(
                        "Confirm the acquisition grid above to continue with a "
                        "cell/cytoplasm stack."
                    )
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
                        "features", "summary", "feature_key", "feature_bundle"):
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

        try:
            preview_vols, _, preview_spacing = _fit_viewer_payload_budget(vols, None, spacing)
        except ValueError as error:
            st.warning(f"3D preview unavailable: {error}")
        else:
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
                    nuclei_diameter_est = estimate_diameter_from_stack(
                        get_model(config.model_type), nuclei
                    )
                # A cell's actual diameter cannot be derived from the nuclei
                # stack alone (P2-5 audit finding); only estimate it when an
                # independent cell/cytoplasm stack was actually uploaded.
                cell_diameter_est = None
                if cells is not None:
                    with st.spinner("Estimating cell diameter (fast 2D pass)…"):
                        cell_diameter_est = estimate_diameter_from_stack(
                            get_model(config.model_type), cells
                        )
            with auto_col0[1]:
                try:
                    anisotropy = nuclei_zstack.spacing.anisotropy
                    xy_spacing_um = None
                    if nuclei_zstack.spacing.complete:
                        assert nuclei_zstack.spacing.x is not None
                        assert nuclei_zstack.spacing.y is not None
                        xy_spacing_um = isotropic_xy_size_um(
                            nuclei_zstack.spacing.x, nuclei_zstack.spacing.y
                        )
                except ValueError as error:
                    st.error(str(error))
                    return
            # Tagged with the exact nuclei stack + channel this suggestion was
            # computed from, so a later upload/channel switch can never reuse
            # it (see resolve_auto_suggestion, P1-1 audit finding).
            suggestion = {
                "nuclei_diameter": nuclei_diameter_est,
                "source_digest": nuclei_digest,
                "source_channel": nuclei_channel,
            }
            if cell_diameter_est is not None:
                suggestion["cell_diameter"] = cell_diameter_est
            if anisotropy is not None:
                suggestion["anisotropy"] = anisotropy
            if xy_spacing_um is not None:
                suggestion["xy_spacing_um"] = xy_spacing_um
            st.session_state["auto_suggest"] = suggestion
            msg = "Auto-detected: "
            if nuclei_diameter_est:
                msg += f"nuclei diameter ≈ **{nuclei_diameter_est:.1f} px**"
            else:
                msg += "could not estimate nuclei diameter (no objects found)"
            if cell_diameter_est:
                msg += f", cell diameter ≈ **{cell_diameter_est:.1f} px**"
            elif cells is not None:
                msg += ", could not estimate cell diameter (no objects found)"
            else:
                msg += ". No cell/cytoplasm stack uploaded — cell diameter left as manual/default."
            if anisotropy is not None:
                msg += (
                    f", XY spacing ≈ **{xy_spacing_um:.4g} µm**, "
                    f"anisotropy ≈ **{anisotropy:.2f}** (from metadata)"
                )
            else:
                msg += ". No Z/XY spacing in metadata — please set anisotropy manually."
            st.success(msg)
            st.rerun()

        if st.button("Run 3D segmentation", type="primary"):
            _run_segmentation(
                nuclei, cells, config,
                nuclei_digest=nuclei_digest, nuclei_channel=nuclei_channel,
                cells_digest=cells_digest, cells_channel=cells_channel,
            )


class _ProgressState(TypedDict):
    """``task``'s real per-key contract: ``start``/``frac`` are always set at
    construction and only ever reassigned to a float; ``estimated`` is the
    one key that is genuinely optional (unset until progress reaches 50%).
    A plain ``dict`` literal here would give every key the same widened
    ``float | None`` type, which is what previously made ``task["start"]``
    and ``task["frac"]`` look possibly-None everywhere they were read.
    """

    start: float
    frac: float
    estimated: float | None


def _run_segmentation(
    nuclei: np.ndarray, cells: np.ndarray | None, config: SegmentationConfig,
    *, nuclei_digest: str, nuclei_channel: int,
    cells_digest: str | None, cells_channel: int | None,
) -> None:
    try:
        working_shape = xy_downsample_shape(nuclei.shape, config.xy_downsample)
        isotropic_xy_scale(nuclei.shape, working_shape)
    except ValueError as error:
        st.error(str(error))
        return
    task: _ProgressState = {"start": time.monotonic(), "frac": 0.0, "estimated": None}
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
            saved = save_result(
                masks[0], masks[1], config,
                nuclei_input_sha256=nuclei_digest, nuclei_channel=nuclei_channel,
                nuclei_dtype=str(nuclei.dtype),
                cell_input_sha256=cells_digest,
                cell_channel=cells_channel,
                cell_dtype=str(cells.dtype) if cells is not None else None,
            )
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
    for key in ("features", "summary", "feature_key", "feature_bundle"):
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

    layer = st.selectbox(
        "Object layer", ["Nuclei", "Cells"] if cell_masks is not None else ["Nuclei"],
        key="feature_object_layer",
    )
    selected_masks = cell_masks if layer == "Cells" and cell_masks is not None else nuclei_masks
    object_type = "cell" if layer == "Cells" else "nucleus"
    fill_holes = st.checkbox(
        "Measure filled outer envelopes", value=False, key="feature_fill_holes",
        help="Default: measure the raw label voxels. Filling enclosed holes changes volume, surface, centroid, axes and solidity together.",
    )
    st.caption(
        f"Measuring {layer.lower()} · {'filled outer envelope' if fill_holes else 'raw label voxels'} · "
        f"voxel size XYZ = {spacing} µm. All objects are retained with review flags."
    )
    if "default" in (result_config.xy_spacing_source, result_config.anisotropy_source):
        st.warning("Physical measurements use assumed voxel spacing. Verify the acquisition calibration before using them in a publication.")

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
        vols = [selected_masks.astype(np.float32, copy=False)]
        chs = [ChannelConfig(lut="gray", opacity=0.5, name="Mask")]

    try:
        preview_vols, preview_mask, preview_spacing = _fit_viewer_payload_budget(vols, selected_masks, spacing)
    except ValueError as error:
        st.warning(f"3D preview unavailable: {error}")
    else:
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
            object_type,
            id(selected_masks),
            fill_holes,
            result_config,
        )
        if st.session_state.get("feature_key") != feature_key:
            features = extract_mask_features(selected_masks, spacing_um=spacing, fill_holes=fill_holes)
            segmentation_provenance = None
            if result is not None:
                provenance_path = result.run_directory / "provenance.json"
                if provenance_path.is_file():
                    try:
                        segmentation_provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
                        if not isinstance(segmentation_provenance, dict):
                            raise ValueError("Run provenance must be a JSON object")
                    except (OSError, ValueError):
                        segmentation_provenance = None
                        st.warning("Saved run provenance could not be read; the download will mark it unavailable.")
            st.session_state["feature_bundle"] = build_mask_feature_bundle(
                selected_masks, features, spacing_um=spacing, object_type=object_type,
                fill_holes=fill_holes, segmentation_config=asdict(result_config),
                segmentation_provenance=segmentation_provenance,
            )
            st.session_state["features"] = features
            st.session_state["feature_key"] = feature_key
        features = st.session_state["features"]
        st.session_state["features"] = features
        if features.empty:
            st.warning("No objects found in the mask.")
        else:
            st.dataframe(features.reset_index(), use_container_width=True)
            reviewed = int(features["qc_status"].eq("review").sum())
            if reviewed:
                st.warning(f"{reviewed} objects need review; consult qc_flags before interpretation. Summaries include these objects.")

        st.download_button(
            "Download feature CSV + provenance", data=st.session_state["feature_bundle"],
            file_name=f"{object_type}_{'envelope' if fill_holes else 'raw'}_features.zip",
            mime="application/zip",
        )

    summary = summarize_features(features)
    st.session_state["summary"] = summary
    mcols = st.columns(4)
    mcols[0].metric("Objects", summary.n_objects)
    mcols[1].metric("Total volume (µm³)", f"{summary.total_volume_um3:,.1f}")
    mcols[2].metric("Mean volume (µm³)", f"{summary.mean_volume_um3:,.1f}")
    mcols[3].metric("Mean sphericity (conditional)", f"{summary.mean_sphericity:.3f}")
    from organoid_analysis.quantification.measurement_policy import MEASUREMENT_INTERPRETATION

    st.caption(MEASUREMENT_INTERPRETATION)

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
            "not compare treatments or estimate viability. **Volume** uses the selected raw-label or filled-envelope geometry in "
            "µm³; **surface area** is a versioned Crofton estimate in µm²; **sphericity** compares the mask "
            "with an equal-volume sphere; **solidity** compares the labeled object with its convex hull; "
            "and principal axes describe a moment-equivalent ellipsoid. Values depend on the entered voxel "
            "spacing and mask quality. QC flags remain in the exported table and these summaries include "
            "flagged objects; define and document study-specific exclusions before group analysis."
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

    st.markdown("**Sphericity vs solidity (descriptive mask shape)**")
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
    # preview tab) become the defaults for the number inputs -- but only when
    # they were computed from the nuclei stack/channel currently loaded (see
    # resolve_auto_suggestion, P1-1 audit finding: a stale suggestion from a
    # previous upload must never silently carry over to a new one).
    suggestions = resolve_auto_suggestion(
        st.session_state.get("auto_suggest"), st.session_state.get("nuclei_identity")
    )

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
        if suggestions and suggestions.get("nuclei_diameter") is not None:
            st.success(
                "Auto-detected: "
                f"diameter ≈ {suggestions['nuclei_diameter']:.1f} px"
                + (
                    f", anisotropy ≈ {suggestions['anisotropy']:.2f}"
                    if suggestions.get("anisotropy")
                    else ""
                )
            )
        elif suggestions:
            st.info(
                "Auto-detect ran but could not estimate a nuclei diameter "
                "(no objects found); using manual/default values below."
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
                "- **Model** — `cpdino-vitb` is the fastest option; `cpsam_v2` is "
                "another available Cellpose foundation model. This pipeline has no "
                "independent target-domain benchmark comparing them, so pick based "
                "on your own QC on representative images rather than a claimed "
                "accuracy ranking.\n"
                "- **QC** — inspect every representative condition for missed, merged and "
                "Z-fragmented objects. Changing a model or threshold changes the measured population.\n"
                "- **Scope** — nuclei/cell masks support nuclei/cell measurements. They do "
                "not automatically define a whole-organoid outer boundary or a viability state."
            )
        nuclei_diameter = st.number_input(
            "Nuclei diameter (px)", 5.0, 200.0, _default_for("nuclei_diameter", backend_defaults.nuclei_diameter),
            help="Set manually, or press 'Auto-detect' after uploading to "
            "estimate it from the image.",
        )
        cell_diameter = st.number_input(
            "Cell diameter (px)", 5.0, 300.0, _default_for("cell_diameter", backend_defaults.cell_diameter)
        )
        anisotropy = st.number_input(
            "Z / XY anisotropy", 0.1, 20.0, _default_for("anisotropy", backend_defaults.anisotropy),
            help="Auto-filled from TIFF metadata when available; otherwise "
            "set from your microscope Z step / XY pixel size.",
        )
        xy_spacing_um = st.number_input(
            "XY pixel size (µm)", 0.001, 100.0,
            _default_for("xy_spacing_um", backend_defaults.xy_spacing_um),
            help="Physical XY voxel spacing. Set to your microscope's pixel "
            "size so volumes/areas are in real units.",
        )
        # Record where anisotropy/xy_spacing actually came from -- metadata
        # (accepted as auto-detected), the hardcoded fallback (no metadata was
        # available), or a manual override -- so the frozen run config never
        # lets an assumed voxel size look the same as a measured one.
        metadata_anisotropy = suggestions.get("anisotropy")
        metadata_xy_spacing_um = suggestions.get("xy_spacing_um")
        anisotropy_source = resolve_spacing_source(
            anisotropy, metadata_anisotropy, backend_defaults.anisotropy
        )
        xy_spacing_source = resolve_spacing_source(
            xy_spacing_um, metadata_xy_spacing_um, backend_defaults.xy_spacing_um
        )
        _SOURCE_LABEL = {
            "metadata": "from TIFF metadata",
            "default": "fallback default, no metadata read",
            "user_override": "manually overridden",
        }
        if xy_spacing_source == "default" or anisotropy_source == "default":
            st.caption(
                "⚠️ No usable Z/XY spacing metadata was read for this stack — "
                "volumes/areas will use the manually set values above, which "
                "default to an assumed voxel size, not a measured one."
            )
        elif "user_override" in (xy_spacing_source, anisotropy_source):
            st.caption(
                f"XY spacing: {_SOURCE_LABEL[xy_spacing_source]}"
                f" · anisotropy: {_SOURCE_LABEL[anisotropy_source]}"
            )
        return SegmentationConfig(
            model_type=model_type,
            nuclei_diameter=nuclei_diameter,
            cell_diameter=cell_diameter,
            anisotropy=anisotropy,
            anisotropy_source=anisotropy_source,
            metadata_anisotropy=metadata_anisotropy,
            xy_spacing_um=xy_spacing_um,
            xy_spacing_source=xy_spacing_source,
            metadata_xy_spacing_um=metadata_xy_spacing_um,
            nuclei_flow_threshold=st.slider("Nuclei flow threshold (unused in 3D)", 0.0, 1.0, backend_defaults.nuclei_flow_threshold, disabled=True),
            cell_flow_threshold=st.slider("Cell flow threshold (unused in 3D)", 0.0, 1.0, backend_defaults.cell_flow_threshold, disabled=True),
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
