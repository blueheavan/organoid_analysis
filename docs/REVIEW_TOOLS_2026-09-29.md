# Review tools adopted from the MCAS paper — 2026-09-29

Source reviewed: Kim et al., *Parallelized bright-field and fluorescence imaging
platform for organoids*, Cell Rep. Methods 6:101601 (2026),
doi:10.1016/j.crmeth.2026.101601 (read from the full-text PDF; statements below
were checked against it). Three ideas were adopted. **None changes a
measurement, an eligibility rule or any gate item**; they add provenance, a
review flag and a viewer.

| # | Idea in the paper | Implementation here | Status |
|---|---|---|---|
| 5 | SAM masks were replaced by manual segmentation where automatic segmentation failed (Fig. S7) | `quantification/mask_provenance.py`: per-object `mask_manually_edited` (`no`/`yes`/`unknown`), `mask_edit_declaration`, `mask_origin`, `labels_sha256`; manifest columns `labels_manual_edit`, `labels_edited_ids`; `analyze-3d --label-edit-log` | provenance only |
| 6 | In-focus plane chosen with a 7×7 Laplacian contrast metric; later-stage organoids exceed a single plane's depth of field | `quantification/focus_qc.py`: LoG (σ 1.5 px) focus and background-corrected signal profiles per object; `z_coverage_suspect` when the focus peak is at the first/last plane, an edge plane keeps ≥ 50 % of peak signal, or the mask touches it | review flag; heuristic thresholds |
| 7 | Whole-plate images browsed with OpenSeadragon deep zoom | `visualization/deepzoom_viewer/` (OpenSeadragon 6.1.1 vendored, BSD-3, SHA-256 in `static/openseadragon.PROVENANCE.txt`); web results tab → *Full-resolution plane review* | display only |

## Why these are review aids, not qualified measurements

* **Manual edits.** An edited object is measured by the same code but produced by
  a different, operator-dependent process that no automatic-segmentation study
  (SG-3) covers. The flag lets results be stratified or re-run without edited
  objects; undeclared imported masks are `unknown`, never assumed automatic.
* **Z coverage.** Border truncation was already flagged when the *mask* reaches
  the first or last plane. The new flag adds the case of a dim object whose mask
  stops short of an edge that still cuts it, and stacks whose sharpest plane is at
  the edge. The 50 % signal ratio and the LoG scale are heuristics chosen for
  review; they have no validated sensitivity or specificity. Other objects are
  excluded plane by plane after XY dilation by the footprint margin, so organoids
  stacked at the same XY position do not contaminate one another (tested).
* **Deep-zoom viewer.** 8-bit, one intensity window per stack, lossless PNG tiles
  by default. Readout uses the voxel-centre convention (voxel `(z, y, x)` centred
  at `(z·sz, y·sy, x·sx)` µm). Non-square XY pixels are shown without aspect
  correction and labelled as such. It is not an annotation tool: SG-3 still
  requires blinded native-resolution 3D annotation under a written SOP.

## Verification

| Check | Result |
|---|---|
| Unit tests | `tests/quantification/test_focus_qc.py` (focus peak inside vs at edge, mask stopping short of a cut edge, stacked organoids, too few planes, flat stack), `tests/quantification/test_mask_provenance.py` |
| Workflow tests | `tests/workflows/test_organoid_measurement_workflow.py` (pipeline vs imported masks, per-ID edits, unknown origin, rejected declarations), `tests/workflows/test_multilevel_measurement_workflow.py` (measurements identical with and without declarations; CLI log hashed into provenance) |
| Viewer | `tests/visualization/test_deepzoom_viewer.py`: level geometry, lossless tile reassembly, contour tiles, budget, one-script self-contained HTML, bundle hash; headless Chromium render with pointer readout, ID markers and Z switching |
| Evidence scope | no file bound by a qualification record was changed; gate unchanged at 4/9 |
