# Validation Plan — Organoid Pipeline (Multilevel 3D Analysis)

Version: 1.0.0
Date: 2026-09-02
Scope: Scientific V&V plan for the `analyze-3d` multilevel 3D organoid analysis (hierarchy, morphology, topology, spatial, QC). The classical morphology/viability validation is documented in `docs/VALIDATION.md`.

This plan predeclares, for each critical requirement, the reference standard, evaluation unit, dataset, independence level, minimum N, metric, and acceptance criteria. Status is filled in during execution; see `docs/VALIDATION_REPORT.md`.

---

## VR-1 Data-contract and identity integrity

- **Requirement:** Registered organoid/cell/nucleus label masks are accepted only when they share a common physical grid and valid integer IDs; identity is never silently relabeled or reshaped.
- **Reference standard:** Defined input contract (`validate_labels`, `validate_inputs`).
- **Evaluation unit:** per input volume.
- **Dataset:** synthetic + any real registered volumes.
- **Independence:** N/A (definitional).
- **Minimum N:** 1.
- **Metric:** error raised on shape mismatch, non-3D, negative or boolean labels, non-finite intensity.
- **Acceptance:** `PASS` if contract is enforced and tests cover each violation type.
- **Status:** `PASS` — enforced by `validation.py`; tests in `test_multilevel3d.py`.

## VR-2 Hierarchy correctness (maximum-overlap parent assignment)

- **Requirement:** Each cell is assigned to its max-overlap organoid; each nucleus to its max-overlap cell; nucleus organoid membership inherits through its cell; direct-overlap audit is retained.
- **Reference standard:** Hand-computed overlap counts on a controlled phantom with known geometry.
- **Evaluation unit:** per object.
- **Dataset:** `synthetic_labels()` phantom (test_multilevel3d.py).
- **Independence level:** independent author computes expected assignments by hand.
- **Minimum evaluable N:** 5 cells, 5 nuclei spanning multiple parents.
- **Metric:** fraction of objects with correct parent; QC flags `crosses_multiple_parents`, `parent_assignment_failed`, `low_parent_overlap`, `direct_organoid_parent_mismatch` set correctly.
- **Acceptance criterion:** 100% correct assignments; all flagged states correct.
- **Uncertainty/tolerance:** deterministic; no stochastic tolerance.
- **Status:** `PASS`.

## VR-3 Physical volume and centroid

- **Requirement:** `volume_um3` = voxel_count × voxel volume; centroid in physical µm.
- **Reference standard:** analytical voxel volume formula, independent centroid computation.
- **Dataset:** synthetic phantom + analytical sphere.
- **Metric:** volume exact (integer × spacing product); centroid agreement.
- **Acceptance:** volume must equal exact voxel volume; centroid within one voxel of independent computation.
- **Status:** `PASS`.

## VR-4 Surface area and sphericity (numerical correctness)

- **Requirement:** surface area from native-spacing level-0.5 marching cubes; sphericity = `π^(1/3) (6V)^(2/3)/A`; no silent clipping to [0,1].
- **Reference standard:** independent analytical sphere phantom for known bias; correctness of formula.
- **Dataset:** 15-voxel-radius sphere (spacing 1).
- **Metric:** measured vs analytical V and A ratios; sphericity value; behavior when sphericity >1.05.
- **Acceptance:** document discretization bias (measured V/A ratios); confirm sphericity >1.05 is flagged, not clipped. Do NOT require bias-free area (bias is expected and documented).
- **Status:** `PARTIAL — formula/numerical correctness PASS on phantom; absolute surface accuracy on real organoids INSUFFICIENT EVIDENCE (requires independent surface reference).`

## VR-5 Contact topology

- **Requirement:** Cell adjacency = voxel-face contact only; anisotropic physical contact area correct.
- **Reference standard:** hand-computed face counts and face areas on a controlled phantom.
- **Dataset:** synthetic_labels() phantom (101/102 X-contact).
- **Metric:** edge existence, contact_area_um2, degree, total_contact_area_um2.
- **Acceptance:** X-normal 5×5 shared plane → contact_area = 25 × (Z×Y); degree correct.
- **Status:** `PASS`.

## VR-6 Spatial features (radial position, surface depth)

- **Requirement:** centroid-to-organoid-centroid distance in µm; normalized radial position using equivalent-sphere radius; distance-to-surface from EDT depth.
- **Reference standard:** independent computation on a controlled phantom; defined operational semantics.
- **Dataset:** synthetic phantom.
- **Metric:** distance values; normalized radial positions in plausible range; depth ≥0.
- **Acceptance:** deterministic values match independent computation; depth semantics documented (0 at/outside surface).
- **Status:** `PARTIAL — computational correctness PASS on phantom; biological meaning of radial bins NOT ASSESSED.`
- **Known edge condition:** a cell centroid outside its organoid mask reports `distance_to_organoid_surface_um = 0` (on/outside surface). This is a documented operational choice, not a failure.

## VR-7 QC flags

- **Requirement:** Formal QC flags generated as a feature-database table (`qc_flags.parquet`), including MAD volume outliers, border touch, too-small, parent assignment, anucleate/multinucleated, direct-organoid mismatch.
- **Reference standard:** independent implementation of robust z-score and flag logic.
- **Dataset:** synthetic phantom with engineered border/outlier/anucleate/multinucleated/mismatch cases.
- **Metric:** exact flag set per object.
- **Acceptance:** all engineered flags produced; MAD==0 degenerate case handled.
- **Status:** `PASS`.

## VR-8 Representative real-data end-to-end behavior

- **Requirement:** The pipeline runs end-to-end on representative registered real organoid/cell/nucleus label volumes and produces physically plausible feature tables.
- **Reference standard:** expert review of feature distributions; no ground-truth claims.
- **Dataset:** real registered volumes (when available/permitted).
- **Independence:** domain expertise for plausibility review.
- **Minimum N:** TBD by available dataset.
- **Metric:** successful export; QC flag rates; physically plausible volume/area/ranges.
- **Acceptance:** runs without error; outputs nonempty features; plausible ranges.
- **Status:** `INSUFFICIENT EVIDENCE — depends on available real registered data; not fully assessed in this pass.`

## VR-9 Reproducibility / determinism

- **Requirement:** The measurement pipeline is deterministic given identical inputs and config.
- **Reference standard:** identical re-run reproducibility.
- **Dataset:** synthetic phantom.
- **Metric:** identical feature tables across two runs.
- **Acceptance:** byte-identical outputs.
- **Status:** `PASS — pipeline is deterministic (no RNG in measurement).`

## VR-10 Failure modes and out-of-domain input

- **Requirement:** Missing spacing, conflicting metadata, shape mismatch, non-finite intensity, empty labels are handled (error or explicit empty result), never a silently wrong number.
- **Reference standard:** defined failure behavior.
- **Dataset:** synthetic malformed inputs.
- **Metric:** correct error/exception or empty-DataFrame path.
- **Acceptance:** each case behaves as defined.
- **Status:** `PASS for tested cases (see VALIDATION_REPORT.md for full list).`

---

## Predeclared acceptance for this pass (multilevel 3D)

Based on the controlled phantoms and the deterministic, well-scoped nature of the measurement, we predeclare:

- VR-1, 2, 3, 5, 7, 9: **PASS** via controlled phantoms.
- VR-4, VR-6: **PARTIAL** — numerical/computational correctness PASS; absolute accuracy / biological meaning INSUFFICIENT EVIDENCE.
- VR-8: **INSUFFICIENT EVIDENCE** — requires representative real registered data.

## Independence level

Audit separation for this pass: Tier D (same agent and context as implementation) → final report records `LIMITED INDEPENDENCE`. Oracle independence: controlled phantoms use hand/independent expected values for hierarchy/topology/QC; morphology uses independent analytical sphere constants.

## Traceability

- Plan mapping to code: see `docs/ALGORITHM_DECISIONS.md` traceability table.
- Tests: `tests/workflows/test_multilevel_measurement_workflow.py`.
- Result outputs: `export_results` → Parquet/JSON hierarchy described in README.
