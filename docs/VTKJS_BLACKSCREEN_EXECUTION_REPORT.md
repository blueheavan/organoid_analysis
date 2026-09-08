# EXECUTION REPORT — vtk.js Viewer Black-Screen Fix

## 1. Symptom
Uploading the real TIFF (`ZeroG-Breast-Cancer-Spheroid-C1.tif`, (56,512,512), uint16,
单通道) produced an **all-black** viewer. The **synthetic 3-channel demo was also black**,
which was the decisive observation: the synthetic data is NOT sparse and reaches
opacity 0.9, so the black screen could not be blamed on data sparsity or opacity.

Console always printed `[viewer] rendered mode=volume dims=...` with **no** vtk.js /
WebGL errors, `errGl: 0`.

## 2. Investigations performed (and results)
| Hypothesis | Test | Result |
|---|---|---|
| Opacity ~ 0 (percentile ramp hides signal) | CPU percentile analysis | **Refuted** — DAPI shell mapped to opacity 0.9 |
| Data domain [0,1] vs TF domain [0,255]/[0,65535] mismatch | Decoded real payload, compared to opacity_pts | **Refuted** — data `max=1.0`, `opacity_pts` x-domain `[0,1]`; they agree |
| Camera not framing volume | on-page diag `camPos/camFocal/bounds` | **Refuted** — camera at (41,41,530) looking at center (41,41,95), dist 435; bounds correct |
| vtk.js version regression (34+ multi-component) | A/B vendored **36.9.2 vs 33.3.2** | **Refuted** — both black in real browser |
| Multi-component volume, `setIndependentComponents(true)` | Interim per-channel actor workaround | **SUPERSEDED** — test occurred before wire-layout correction |

## 3. Corrected root cause
The earlier attribution to vtk.js independent-component compositing was incorrect:
it was made before the Z-fastest wire-layout defect was found. After correcting the
payload to x-fastest `(Z,Y,X,C)`, the required multi-component single-volume path
renders correctly. The objective root cause is the wire-layout defect documented in
Round 2, not `IndependentComponents`.

## 4. Current implementation
The interim per-channel actor workaround has been removed. Current implementation is
**Path A**: one interleaved multi-component `vtkImageData`, one mapper, one `vtkVolume`,
`IndependentComponents=true`, and per-component color/opacity transfer functions.

Additional robustness (Task 3): NaN/Inf guard in `pack_channel_stack`
(`nan_to_num` before packing). No hardcoded intensity domain exists. Raw float32 scalar
values are preserved; percentiles select only the initial display window.

## 5. Evidence (automated, no monitor)
Screenshots + Pillow analysis in `vtk_streamlit_volume/tests/artifacts/`:

| Scenario | non_black_fraction | cyan_fraction | max_px |
|---|---|---|---|
| Synthetic 3-ch (128x128x96) | **0.213** | **0.188** | 255 |
| Sparse TIFF-like (5.2% nz, uint16, 0-65535) | **0.270** | **0.039** | 255 |

All well above the `non_black > 0.005` threshold.

## 6. Regression tests
`tests/test_render_smoke.py` — headless Chromium + SwiftShader, screenshots `#viewer`,
asserts `non_black_fraction > 0.005` and `cyan_fraction > 0` for both the synthetic and
the sparse/TIFF-like volume. Added to the suite.

```
$ pytest vtk_streamlit_volume/tests/
16 passed
```
(14 unit + 2 render smoke)

## 7. Acceptance checklist
- [x] Sparse TIFF-like volume renders as a visible cyan 3D blob (non_black=0.27, cyan=0.039)
- [x] All-zero / degenerate / NaN data are guarded (no hardcoded domain; opacity guard,
      nan_to_num)
- [x] Multi-channel rendered as one multi-component volume with independent TFs
- [x] No hardcoded strength domain in logic
- [x] `pytest tests/` passes (16)
- [x] EXECUTION_REPORT recorded

## 8. Corrected Task 5 (automated visual confirm, no monitor)
Screenshots captured via Playwright headless + SwiftShader and analyzed with Pillow;
artifacts retained under `tests/artifacts/` for manual review by the user:
`smoke_synthetic.png`, `smoke_sparse.png`, `demo_synthetic.png`.

---

# Round 2 — Voxel Layout + Split 2D/3D Viewer

## 9. Confirmed layout root cause

The colored-but-scrambled output was caused by the Python wire layout, not by
base64, dtype decoding, or transfer functions.

The old path transposed each normalized channel from `(Z,Y,X)` to `(X,Y,Z)`,
made that array C-contiguous, and serialized it. Its interleaved linear index was:

```
q_old = c + C * (z + Z * (y + Y * x))
```

This makes **Z the fastest spatial dimension**. vtkImageData dimensions `(X,Y,Z)`
require X-fastest point ids:

```
q_vtk = c + C * (x + X * (y + Y * z))
```

For the real `(56,512,512)` stack, the wrong 56-voxel stride appeared as the
measured ~82-pixel periodic stripes.

## 10. Minimal layout fix and code review

- `core.py:141` `_voxel_major_zyxc()` now converts `(C,Z,Y,X)` to contiguous
  `(Z,Y,X,C)` voxel-major wire data.
- `core.py:203` serializes that array directly. No spatial transpose/reshape/
  Fortran ravel remains in the volume path.
- `core.py:295` alone converts metadata to `dims=(X,Y,Z)`.
- JS attaches the voxel-major interleaved array directly to one multi-component
  `vtkDataArray`; it performs no deinterleave or spatial reshape/transpose.
- JS now hard-fails if `data.length !== X*Y*Z*C` (`core.py:528`).
- Pure Python corner roundtrip verifies value 1000 at
  `5*64*64 + 10*64 + 20`; multi-channel coordinate-coded roundtrip verifies
  voxel-major interleaving.

## 11. Objective layout evidence

All browser tests use headless Chromium + SwiftShader screenshots and Pillow/
NumPy analysis. Artifacts are under `vtk_streamlit_volume/tests/artifacts/`.

| Test | Result | Artifact |
|---|---:|---|
| Corner voxel `(5,10,20)` | projection error **0.65 px** (latest full run) | `layout_corner_viewer.png` |
| Gradient `z*100+y+x/100` | NCC **0.7935**, y corr **0.8103**, x corr **0.0665** | `layout_gradient_viewer.png` |
| Orientation | one documented **global Y flip** (VTK view-up convention) | same |
| Real TIFF MIP vs `arr.max(axis=0)` | NCC **0.9911** | `spheroid_c1_mip.png`, `spheroid_c1_mip_viewer.png` |
| Real TIFF viewer MIP non-black | **0.1563** | `spheroid_c1_mip_viewer.png` |

The gradient NCC threshold is 0.75. The requested `x/100` term occupies less
than one 8-bit screenshot level after full-volume normalization, so direction is
also protected by the exact Python coordinate roundtrip and a positive column
correlation assertion.

## 12. Split 2D/3D preview

The generated viewer now has:

- left: ordinary 2D canvas displaying one Z slice and a client-side slice slider;
- right: existing vtk.js GPU volume/MIP viewport;
- shared current channel, LUT, raw-domain intensity window, opacity, and visibility;
- client-side updates only (no Streamlit rerun);
- manual `gray` LUT for brightfield-style display; no unreliable auto-detection;
- default fluorescence LUT order: C0 cyan, C1 magenta, C2 yellow, C3 green.

Shared-state screenshot test changes LUT to magenta, window low to 0.2, and slice
to z=10, then checks both sides:

| Pane | non-black fraction | magenta fraction | Artifact |
|---|---:|---:|---|
| 2D | **0.1663** | **0.1331** | `split_2d_magenta.png` |
| 3D | **0.0966** | **0.0650** | `split_3d_magenta.png` |

## 13. Round 2 acceptance

- [x] corner voxel at expected location (0.65 px error)
- [x] three-axis gradient direction correct; documented global Y flip only
- [x] real TIFF MIP structural match (NCC 0.9911)
- [x] real TIFF non-black fraction > 0.005 (0.1563)
- [x] 2D/3D split both non-black and shared state verified
- [x] screenshots retained in `tests/artifacts/`
- [x] full test suite passes: **23 tests** at the recorded Round 2 run

Final Streamlit iframe probe (`http://127.0.0.1:8501`) also passed after a fresh
server restart: 2D non-black **0.4907**, 3D non-black **0.3585**. Artifacts:
`demo_served_2d.png`, `demo_served_3d.png`.

## 14. Preview-only LUT correction

After user comparison with ImageJ, fluorescence LUTs were changed to linear,
pure-color preview mappings: black→cyan, black→magenta, black→yellow, and
black→green. The high end no longer drifts toward white. This changes only the
browser transfer functions; source TIFF values and packed scalar data are not
modified. The intensity Window remains intentionally shared between 2D and 3D.

## 15. Path A and rendering-quality correction

Runtime assertions now lock the architecture and quality invariants:

- implementation: `multi-component-single-volume`;
- one `vtkVolume`, three components in the 3-channel test;
- `IndependentComponents=true`, per-component TF and opacity-unit distance;
- scalar values remain raw float32; metadata carries `actual_min`, `actual_max`,
  `suggested_lo`, and `suggested_hi` from the same channel array;
- Window sliders span the full actual data domain and never rewrite the payload;
- trilinear interpolation and `shade=false` in every quality mode;
- `autoAdjustSampleDistances=false`;
- opacity unit distance is the minimum spacing (0.414 µm in the real dataset);
- quality differs only by sample distance: Smooth 0.5, Balanced 0.2 (default),
  Sharp 0.1 µm;
- vtk.js 33.3.2 WebGL has no public `setJittering` method: jitter is built into
  its OpenGL volume shader and always active. The compatibility call enables it
  explicitly on vtk.js variants that expose the method;
- `GenericRenderWindow.resize()` already allocates at DPR. At DPR=2 the automated
  assertion measured framebuffer/CSS ratios `[2,2]`; an explicit fallback covers
  backends that return CSS-pixel buffers.

Artifact: `multi_component_path_a.png`. Final suite at this stage: **26 passed**.

Fresh Streamlit iframe runtime probe confirmed the shipped page (not only standalone
test HTML): implementation `multi-component-single-volume`, volume count 1,
component count 3, independent components on, auto-adjust off, sample distance 0.2,
shade off. Served non-black fractions were 2D **0.5968** and 3D **0.3695**.

## 16. Mask overlay + closed-loop integration

Goal (user-confirmed Round 3): close the loop inside the vtk viewer —
**Upload → Preview → Segment → Analyze**. Two new deliverables:

### 16.1 `mask_features.py` — per-object features from a label mask
Pure functions; full test coverage (`tests/test_mask_features.py`, 12 tests):

- `extract_mask_features(mask, spacing_um)` → per-label `DataFrame` with
  `volume_um3`, `surface_area_um2`, `sphericity`, `solidity`,
  `major/minor/least_axis_um`, `elongation`, and physical `centroid_*_um`.
- `summarize_mask(mask, spacing_um)` → `MaskSummary` (count, total/mean/median/
  min/max volume, mean sphericity/solidity, mean axis lengths).
- **skimage 0.26.0 notes (verified by test):** `regionprops` in 3D has no
  `perimeter`/`surface_area`/`sphericity`/`axis_least_length`; passing `spacing`
  to `p.area` already yields **physical µm³** (so volume is NOT multiplied by
  voxel volume again). Surface area is computed via `marching_cubes` +
  `mesh_surface_area`; sphericity = `π^(1/3)(6V)^(2/3)/A`; principal-axis
  diameters derive from physical `inertia_tensor_eigvals`
  (`semi_axis² = 5·eigval/n`).

### 16.2 Mask overlay in the viewer
- `ViewerSpec.mask` + `build_mask_overlay_payload(mask, spacing_um, alpha)`:
  packs shape, **uncompressed** base64 uint32 labels (browser zlib unavailable),
  a marching-cubes surface, a per-label palette, and blend alpha.
- JS: `getMaskLabels()` decodes labels; `drawMaskContours()` paints per-label
  boundary pixels on the 2D slice (palette-coded, alpha-blended via double
  `putImageData`); `buildMaskOverlay()` adds a translucent 3D surface actor;
  `install()` draws overlays in the non-surface branch.
- `build_viewer_payload(..., mask_overlay=...)`, `st_volume_viewer(mask_overlay,
  mask_overlay_alpha)`, and `spec_to_dict()["mask"]` added.

**Browser poly-data mapper fix (root-caused this round):** `DD.vtkPolyDataMapper`
(OpenGL) exposes **no** `setInputData`, and `RC.vtkPolyDataMapper` does not exist —
so the previous surface/mask mapper selection would throw
`TypeError: mapper.setInputData is not a function`. The mapper is now resolved as
`RC.vtkMapper` (which provides `setInputData`) — verified: mask overlay renders
without console errors, 2D contour pixels appear (label 2 = blue, 2094 px at the
default mid slice), 3D overlay adds a translucent surface (non-black jumps from
0.04 → 0.41 on the standalone render smoke).

### 16.3 `segmentation_workspace.py` — closed-loop UI
`streamlit run vtk_streamlit_volume/segmentation_workspace.py`:

1. **Upload & preview** — nuclei (required) + cytoplasm (optional) TIFF, spacing
   from anisotropy (`z = 0.414·anisotropy`); preview with the browser volume
   renderer. Reuses `app.segmentation` (`read_stack`, `validate_stacks`) and
   `ChannelConfig`.
2. **Run 3D segmentation** — background thread via `segment_stacks`/
   `create_model`/`save_result` with a live progress bar + ETA; keeps the heavy
   Cellpose model in `st.cache_resource`. `KMP_DUPLICATE_LIB_OK` set at import so
   torch/Cellpose start reliably; `get_accelerator` guarded so the app still runs
   preview-only if Cellpose is unavailable.
3. **Segmentation results** — nuclei/cell counts + intensity mask volume rendered
   **with mask overlay** (2D per-label contours + 3D translucent surface) + a
   downloadable masks/summary zip.
4. **Analysis** — per-object feature table (`extract_mask_features`), 4 summary
   metrics, volume distribution bar chart, sphericity-vs-solidity scatter.

Acceptance this round:
- mask-feature unit tests: **12**; mask-overlay payload tests: **4**; render
  smoke: **1** → full suite **43 passed**.
- Integrated app boots at `http://localhost:8599` (HTTP 200, headless), renders
  title / uploader / sidebar / three tabs with **zero** `stException` and **zero**
  console errors. Segmentation button requires Cellpose weights (validated upstream
  in `app.segmentation`); preview + analysis are fully testable without weights.
