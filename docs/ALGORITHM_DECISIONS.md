# Algorithm Decisions — Organoid Pipeline

Version: 1.0.0
Date: 2026-09-02
Applies to: multilevel 3D analysis (`analysis.analyze-3d`) and classical organoid morphology/viability analysis.

## Decision provenance

This document records the scientific rationale for consequential algorithmic choices in the pipeline. For each decision: the scientific task, credible alternatives considered, the selected method, the alternatives rejected and why, and the validation status. It is a living record and is updated when a method is added or changed.

---

## D1. Parent assignment by maximum voxel overlap

**Task:** Assign each cell to an organoid, and each nucleus to a cell, within a registered 3D instance hierarchy.

**Alternatives considered:**
- Nearest centroid (Euclidean) assignment.
- Boundary/contact overlap.
- Voxel-count majority.

**Selected:** Maximum voxel overlap (argmax of shared voxel count), with deterministic tie-break to the smallest parent ID.

**Rationale:** Voxel overlap measures physical co-occupancy of the child within each parent. For nested 3D structures with potentially non-convex shapes, overlap reflects where the object's mass actually lies. Centroid distance can mis-assign a large, elongated or fragmented child whose centroid falls outside its true parent (e.g., a nucleus that protrudes across a boundary). Overlap is robust to such geometry and is a well-established approach for label-map ancestry.

**Rejected:**
- Nearest centroid — fails for elongated/border-crossing objects; the centroid may lie in the neighbor's volume.
- Boundary contact — measures surface adjacency, not containment; two overlapping objects may share no boundary face.

**Failure modes:** A child straddling multiple parents gets the parent with plurality overlap. This is surfaced to the user via `crosses_multiple_parents`, `parent_overlap_fraction < low_parent_overlap_fraction`, and `parent_assignment_failed` QC flags rather than being silently dropped.

**Validation:** `PARTIAL — PASS for algorithmic correctness on controlled phantoms (test_multilevel3d.py); real-world assignment accuracy against expert annotation NOT ASSESSED.`

**Strength labeling:** established method (maximum overlap assignment).

---

## D2. Nucleus-organoid membership inherited through the cell

**Task:** Determine the organoid that contains each nucleus.

**Selected:** The nucleus's organoid is the organoid assigned to its parent cell. The nucleus-to-organoid direct overlap is retained as an audit field (`direct_organoid_id`) but never allowed to override the cell-mediated hierarchy.

**Rationale:** The two-level hierarchy Cell→Organoid, Nucleus→Cell defines a deterministic nesting. A nucleus should belong to the organoid of the cell that houses it, even if the nucleus's own voxels straddle an organoid boundary. This enforces a consistent interpretation of "this nucleus belongs to this cell, which belongs to this organoid."

**Rejected:** Direct nucleus→organoid overlap as the governing assignment — this would break the transitive cell-chain consistency and create contradictory memberships.

**Failure modes:** A mismatch between the cell-derived organoid and the direct nucleus overlap is recorded as the `direct_organoid_parent_mismatch` QC flag, not silently resolved.

**Validation:** `PASS — controlled phantom (test_multilevel3d.py::test_nucleus_inherits_organoid_from_assigned_cell_and_retains_direct_overlap_audit).`

**Strength labeling:** heuristic (design choice to enforce a consistent hierarchy).

---

## D3. True-3D morphology via marching cubes level-0.5 surface

**Task:** Compute physical surface area and sphericity from a 3D label mask.

**Selected:** `skimage.measure.marching_cubes` at `level=0.5`, native physical spacing, step size 1, one-voxel zero padding, no mesh smoothing. Volume = voxel count × voxel volume. Sphericity = `π^(1/3) (6V)^(2/3) / A`.

**Rationale:** Marching cubes at native spacing with no smoothing preserves the underlying voxel geometry and is the standard surface estimator for binary masks. Zero padding closes the crop surface so an object fully inside the image reports a closed surface. SPHERICITY uses the classic isoperimetric ratio; perfect sphere ⇒ 1.0.

**Rejected:**
- MIP / single-slice / density projections — explicitly excluded by project principle (never treat 2D projection as 3D measurement).
- Smoothed or downsampled mesh for measurements — used only for preview, never for reported area. Documented discretization bias is retained rather than hidden.

**Failure modes / known bias:** Voxel discretization inflates area and can push sphericity slightly above 1.0. Values >1.05 are flagged (`sphericity_above_geometric_range`), never silently clipped into [0,1]. Empirical check on a 15-voxel-radius sphere (spacing 1): volume ratio ≈ 1.013, area ratio ≈ 1.098, sphericity ≈ 0.919. This is the documented discretization bias, not a defect.

**Validation:** `PARTIAL — PASS for volume/centroid/axis analytical agreement on phantoms; surface-area absolute accuracy on real biological organoids NOT ASSESSED (requires independent surface measurement reference).`

**Strength labeling:** established method (skimage marching cubes surface).

---

## D4. Face-contact topology (cell adjacency), not centroid/kNN/Delaunay

**Task:** Define which cells are "neighbors" for a contact graph.

**Selected:** Voxel-face adjacency only. An edge exists only when two cell labels share at least one voxel face. Edge weight = physical contact area computed from the anisotropic face area (Z-face → Y×X, Y-face → Z×X, X-face → Z×Y).

**Rationale:** Face contact is a direct statement of physical adjacency/contact. Centroid distance, k-nearest-neighbor, and Delaunay triangulation infer proximity, not contact, and are deliberately not treated as contact.

**Rejected:** centroid-distance, kNN, Delaunay-based adjacency for contact claims.

**Validation:** `PASS — controlled phantom verifies X-face contact area = 25 × Z×Y for a 5×5×5 cell pair with 5×5 shared X-face (test_multilevel3d.py).`

**Strength labeling:** defined operational definition.

---

## D5. Robust MAD volume-outlier detection

**Task:** Flag unusually large/small organoids/cells/nuclei by volume.

**Selected:** Robust z-score: `0.67448975 × (V − median) / MAD` with a threshold of 3.5 (constant `mad_z_threshold`), i.e. equivalent to ~3.5 z-units on the robust scale (`MAD * 1.4826`). With MAD = 0 (common in near-identical synthetic populations), finite values differing from the median are flagged.

**Rationale:** MAD is robust to outliers, unlike mean/std which are inflated by the very outliers being detected. This prevents a single extreme object from masking the next-level outliers.

**Validation:** `PASS — controlled phantoms; real-data outlier calibration NOT ASSESSED.`

**Strength labeling:** established robust-statistics technique (median + MAD).

---

## D6. Radial position normalized to equivalent-sphere radius

**Task:** Locate a cell within its organoid (core vs periphery).

**Selected:** `normalized_radial_position_equivalent_radius = distance(cell_centroid, organoid_centroid) / r_eq`, where `r_eq = (3V/(4π))^(1/3)` is the equivalent-sphere radius of the organoid's volume.

**Rationale:** Normalizing to an equivalent sphere yields a dimensionless 0..~1+ index. Core = ≤0.5, periphery = ≥0.8 (configurable `core_max_normalized_radial_position`, `peripheral_min_normalized_radial_position`).

**Rejected:** Local-radius normalization — the codebase does not estimate a per-direction local organoid radius and does not claim to.

**Validation:** `PASS — computing path is deterministic/gradient-tested; biological meaning NOT ASSESSED.`

**Strength labeling:** heuristic; explicitly documented as equivalent-radius, not local-radius.

---

## D7. Classical morphology marker measurements & viability state gating

**Task:** Classical organoid-level Calcein/PI viability-like states.

**Selected:** Marker means measured on raw (unnormalized) arrays within the outer envelope; background = median in a physical shell; corrected mean = foreground mean − background median; saturation detection uses an explicit camera limit or integer dtype max; control-scaled endpoints learned per batch from live/dead controls; marker scaled as `(mean − low)/(high − low)`; states gated at 0.60 (high) and 0.30 (low) with indeterminate between.

**Rationale, rejected alternatives, and caveats are fully documented in docs/METHODS.md.** This is an intensity-based signal-pattern characterization, not a cell-count or a validated cell-viability fraction.

**Validation:** `PARTIAL — controlled phantom agreement on simplified marker patterns; biological calibration against matched controls and independent viability measurements NOT ASSESSED.`

**Strength labeling:** empirically calibrated / heuristic.

---

## Rejected overall approaches

- **Replace the existing segmentation model** — preserved by requirement. Segmentation is upstream; this pipeline consumes its output.
- **Use centroid-based hierarchy** — rejected as less robust than overlap (D1).
- **Use 2D projections for any quantitative 3D measurement** — rejected by project principle.

## Traceability

| Decision | Code location |
|---|---|
| D1 | `src/analysis/multilevel3d/hierarchy.py` |
| D2 | `src/analysis/multilevel3d/pipeline.py` (lines ~193-212) |
| D3 | `src/analysis/features.py::geometry`, `::surface_mesh` |
| D4 | `src/analysis/multilevel3d/topology.py` |
| D5 | `src/analysis/multilevel3d/qc.py::_volume_outliers` |
| D6 | `src/analysis/multilevel3d/spatial.py` |
| D7 | `src/analysis/features.py::marker_measurements`, `src/analysis/viability.py`, `src/analysis/report.py` |
