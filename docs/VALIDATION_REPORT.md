# Validation Report — Multilevel 3D Analysis

Version: 1.1.0
Date: 2026-09-09
Classification: S3 (Inferential) for the pipeline as a whole (see docs/SCIENTIFIC_SPEC.md for evidence); this report's scope remains the `analyze-3d` multilevel measurement step specifically, which is S2 on its own — the S3-triggering component (`stats.py`) belongs to the classical `analyze` pipeline and is validated separately (docs/ALGORITHM_DECISIONS.md D11).
Audit separation tier: **Tier D** (same agent and context as implementation) → `LIMITED INDEPENDENCE`.

> 2026-09-09 update: This report's reproducibility baseline is refreshed from `3c9793c` to HEAD `3b1fda9` (5 commits ahead), and an S1 data-handling verification entry (VR-11) covering the unified TIFF axis/spacing reader behavior (docs/ALGORITHM_DECISIONS.md D12) is added to `docs/VALIDATION_PLAN.md`. No measurement algorithm changed; D12 is a data-reader/unit-handling behavior unification verified by its own regression tests and re-loading of the real sample images.
>
> 2026-09-09 second update (P1/P2/P3 audit-remediation pass): baseline refreshed again to HEAD `e62e052` (2 commits ahead of `3b1fda9`). This pass fixed two P1 blocking scientific-integrity defects in the Web Cellpose upload path (a cross-sample auto-suggestion leak, and missing physical-acquisition-grid validation between a nuclei stack and a separately uploaded cell/cytoplasm stack), added Cellpose run provenance (provenance.json), stopped fabricating a nucleus-derived cell-diameter estimate, made the pairwise statistical-contrast CSV self-describing about its effect-size scale and FDR scope, fixed a CI defect that silently prevented `pixi run ci` from ever reaching `typecheck`, and downgraded the VR-11 real-sample-reload claim to PARTIAL with a data-identity manifest (`docs/validation_data_manifest.csv`) since `data/images/*` is gitignored and not repository-reproducible. See finding F-6 below and the two commit messages at `b61cc0f` and `e62e052` for full detail. No classical-pipeline or Cellpose numerical output changed in this pass.

This report records the scientific V&V of the `analyze-3d` multilevel 3D organoid analysis against `docs/VALIDATION_PLAN.md`. It separates **software correctness** (verification, this pass) from **biological validity** (which remains NOT ASSESSED / INSUFFICIENT EVIDENCE unless a fully independent real-data study is performed). The classical morphology/viability workflow is validated separately in `docs/VALIDATION.md` and is not re-validated here.

> **Bounded verdict caveat:** This pass evaluates algorithmic/numerical correctness on controlled phantoms. It does **not** establish biological or clinical accuracy on real organoid tissue. Real-data V&V is required before the multilevel outputs are used to draw biological conclusions.

---

## Intended-use scope (from SCIENTIFIC_SPEC.md §2)

Research-use quantitative morphology of registered 3D organoid→cell→nucleus label hierarchies. Not for clinical diagnosis, treatment selection, or go/no-go decisions.

---

## Evidence statuses

Per item: `PASS`, `PARTIAL`, `FAIL`, `NOT ASSESSED`, `INSUFFICIENT EVIDENCE`, `NOT APPLICABLE`.

## 1. Software correctness (verification)

| Item | Result | Evidence |
|---|---|---|
| Data-contract enforcement (VR-1) | `PASS` | `validation.py` rejects shape mismatch, non-3D, negative/bool labels, non-finite intensity; covered by `tests/workflows/test_multilevel_measurement_workflow.py` |
| Hierarchy assignment (VR-2) | `PASS` | independent phantom: organoids 10/20; cell101/102→organoid 10 (overlap 1.0); nucleus1001→cell101, nucleus1002→cell102; both cells fully within organoid 10 |
| Physical volume + centroid (VR-3) | `PASS` | unit-scaling oracle: `np.ones((9,13,17))` × (2,1,1) → volume exactly 8×, area exactly 4×, sphericity invariant (test_geometry.py); multilevel phantom: 125-voxel cell (spacing 2.0,0.65,0.65) → 125×2.0×0.65×0.65 = 105.625 µm³ (test_multilevel_measurement_workflow.py) |
| Surface area/sphericity formula (VR-4) | `PARTIAL` | analytical 15-voxel sphere: volume ratio 1.0007, area ratio 1.0895, sphericity 0.9183 — documented discretization bias; >1.05 flagged not clipped. Absolute real-tissue area accuracy: `INSUFFICIENT EVIDENCE` |
| Face-contact topology (VR-5) | `PASS` | independent phantom: 5×5 shared X-normal face → contact area exactly 25×(2.0×0.65)=32.5 µm²; degree=1 |
| Spatial features (VR-6) | `PARTIAL` | computational correctness PASS on phantom; biological meaning of core/periphery bins `NOT ASSESSED` |
| Formal QC flags (VR-7) | `PASS` | border, MAD-outlier, too-small, parent-failed, anucleate, multinucleated, direct-mismatch all produced; MAD==0 path handled |
| Reproducibility/determinism (VR-9) | `PASS` | identical feature tables across two runs (organoid/cell/nucleus/topology/qc) |
| Failure modes / out-of-domain (VR-10) | `PASS` | missing spacing, conflicting metadata, shape mismatch, non-finite intensity handled explicitly |
| Full automated suite | `PASS` | `pixi run test` at HEAD `e62e052`: 222 passed, 8 skipped (Playwright-gated headless render tests), 8 render-marked deselected by default; `pixi run test-render`: 8 passed. `pixi run check` (notebook validation): exits 0 (P3-1 fix -- previously exited 1 with no notebooks tracked, which silently blocked `pixi run ci`'s depends-on chain from ever reaching `typecheck`). Lint (`ruff`) clean. `pixi run python -m organoid_analysis --help` (CLI smoke): exits 0. |
| Type check backlog (P3-2, non-blocking by design) | `PARTIAL` | `pixi run typecheck` (mypy) now genuinely executes as part of `pixi run ci` (see the `check` row above) rather than being silently skipped: **194 errors** at HEAD `e62e052`, down from the **197-error** pre-existing baseline at `e0ed85e` (verified by a line-number-independent diff of both runs' error signatures). This pass's own new code introduced zero net new mypy errors; the -3 net change is an incidental fix (a return-typed-`object` ZStack parameter this pass had to touch anyway) unrelated to the audited defects. The remaining 194 errors are unrelated pre-existing debt, not re-triaged in this pass, and this row is intentionally reported `PARTIAL` rather than folded into the `PASS` row above -- mypy is not clean and this must not be read as "all checks pass." |

### Independent-oracle note (VR-2/3/5)

Hierarchy, volume, and contact-area acceptance were computed **independently by hand** from the phantom geometry, not from the implementation itself. This satisfies the "independent evidence where feasible" principle for those items.

## 2. Numerical correctness details (VR-4)

Analytical sphere, radius 15 voxels, spacing 1:
- `V_measured / V_true = 1.0007` (voxel discretization; measured 14147 voxels vs analytical 14137.2).
- `A_measured / A_true = 1.0895` (marching-cubes discretization bias, documented).
- `sphericity = 0.9183`.

The pipeline does not force sphericity into [0,1]; values >1.05 are flagged. This is the documented and intended behavior, not a defect.

## 3. Representative real-data validation (VR-8)

**Result: `INSUFFICIENT EVIDENCE — NOT ASSESSED` for biological accuracy.**

- **Reader path (S1, D12, PARTIAL — not repository-reproducible):** 5 real sample images (identity, SHA256, and size recorded in `docs/validation_data_manifest.csv`) were loaded through the reader path (`load_zstack`) at HEAD `3b1fda9` on the machine that ran this validation; the P1/P2/P3 pass that advanced HEAD to `e62e052` did not modify `load_zstack`, `parse_spacing_ome`, `canonical_czyx`, or `to_um`, so this result is unaffected and was not re-run at `e62e052`. All returned correct `ZYX` volumes (uint16), with physical spacing read from embedded metadata where present: `Human-Colon-Organoids-C1.tif` → 2.0 µm, `ZeroG-Breast-Cancer-Spheroid-C1.tif` → 1.5 µm; the remaining samples carry no spacing metadata and correctly report `None` rather than fabricating a value. **P2-7 audit finding:** `data/images/*` is gitignored and was never committed (`git ls-files data/images/` returns only `.gitkeep`), so a fresh clone of this repository cannot reproduce this specific re-load on its own -- only the manifest's recorded identity is repository-reproducible. This is downgraded from the earlier unqualified "PASS — real-data-backed" wording for exactly that reason.
- **Multilevel measurement (VR-8):** still `INSUFFICIENT EVIDENCE`. No real **registered organoid+cell+nucleus label volume set** exists in the repository (all prior Cellpose runs under `results/segmentation_output/` report `cell_count: 0`, so they cannot form the three-level hierarchy the pipeline requires). Lifting this requires user-provided real registered label volumes; the synthetic phantom over-estimates real-tissue conditions (high contrast, regular ellipsoids). **No real-data biological accuracy claim is made.**

## 4. Parameter / calibration validity

- All consequential parameters are catalogued with provenance in `docs/PARAMETERS.md` (status: `Partially validated` / `Not validated on real data`).
- Thresholds such as `mad_z_threshold=3.5`, core/peripheral radial bins (0.5/0.8), and `low_parent_overlap_fraction=0.5` are heuristic operational definitions, **not** biologically calibrated cutoffs. Sensitivity analysis on real data: `NOT ASSESSED`.

## 5. Statistical validity

The multilevel measurement step makes no inferential claims (no hypothesis tests, no predictive models). Statistical summaries for the classical workflow are addressed in `docs/VALIDATION.md`. Statistical unit and pseudoreplication caveats are documented in `docs/SCIENTIFIC_SPEC.md` §5–6 (`NOT APPLICABLE` to the measurement-only multilevel step in this pass).

## 6. Data integrity

Metadata supplied by the caller is preserved verbatim and never invented (`_add_metadata`). Units are explicit (µm, µm², µm³). Axes are canonicalized to ZYX.

## 7. Audit separation

**Tier D** — performed by the same agent and context as the implementation review. `LIMITED INDEPENDENCE`. Compensating evidence: controlled-phantom acceptance criteria were independently computed by hand for the main claims (hierarchy, volume, topology), reducing (but not eliminating) the self-referential risk.

## 8. Known failure modes / unsupported claims / limitations

- Surface-area absolute accuracy on real tissue is not established (known discretization and segmentation dependence).
- A cell centroid outside its organoid reports `distance_to_organoid_surface_um = 0` (on/outside surface) — operational definition, documented.
- Real-data segmentation/registration quality is an upstream dependency not controlled by this module; the module refuses mismatched grids rather than resampling silently.
- No clinical, diagnostic, or biological-endpoint claim is made.

## 9. Unresolved findings

| ID | Severity | Description | Status |
|---|---|---|---|
| F-1 | P4 | `cli.py` had a no-op `key.replace("_","_")` — misleading, fixed | Resolved (P3/P4, no scientific impact) |
| F-2 | P4 | In-progress uncommitted label-compaction refactor (`src/organoid_analysis/quantification/labels.py`) deduplicates logic across 4 files; behavior preserved (all tests green) | Open — refactor itself is beneficial; commit when ready |
| F-3 | P2 | `analyze-3d`'s exported `analysis_summary.json` had zero code/environment provenance (no git commit, source hash, or package versions), unlike the classical `analyze` pipeline | Resolved 2026-09-03 — `_run_multilevel` now writes a `provenance` block (git commit + dirty flag, `multilevel3d/*.py` source hashes, package versions); regression test in `test_multilevel3d.py` |
| F-4 | P2 | This report's own Gate 11 baseline cited git commit `280bd06`, which does not exist in this repository's actual history | Resolved 2026-09-08 — Gate 11 now cites `3c9793c` (an actual commit in this repository's history, `git log` verified) with package versions read directly from the installed pixi environment, not carried forward from the unreproducible historical entry |
| F-5 | P3 | This report's Gate 11 baseline and full-suite row went stale after HEAD advanced past `3c9793c` (D12 TIFF reader unification and subsequent commits), and the S1 data-handling change (D12) lacked its own verification entry | Resolved 2026-09-09 — Gate 11 refreshed to HEAD `3b1fda9`; full-suite counts and lint status updated; VR-11 added to `docs/VALIDATION_PLAN.md` covering the D12 axis/spacing/unit behavior (see the S1 verification entry) |
| F-6 | P1 (2 findings) | Independent audit found two blocking scientific-integrity defects in the Web Cellpose upload path: (1) an auto-detected diameter/spacing/anisotropy suggestion could silently leak across a new upload or a channel switch, including mislabeling stale spacing as "from metadata"; (2) a nuclei stack and a separately uploaded cell/cytoplasm stack were only checked for matching array shape, never for matching physical voxel spacing, so mismatched acquisitions could be silently treated as registered channels | Resolved 2026-09-09 at HEAD `e62e052` (commit `b61cc0f`) — suggestions are now tagged with and validated against the source stack's content digest + channel (`resolve_auto_suggestion`); nuclei/cell metadata is compared via `compare_registered_grid` using the existing SPACING_RTOL/SPACING_ATOL_UM contract, blocking on a real conflict and requiring explicit per-upload user confirmation when metadata is incomplete. Regression tests: `test_segmentation_workspace.py`, `test_voxel_spacing.py` (compare_registered_grid cases). Gate 11 baseline refreshed to `e62e052`; full-suite counts updated below. |

No P0/P1 findings **for this report's scope** (the `analyze-3d` multilevel module). An independent scientific-software-development-validation audit on 2026-09-03 found two P1 findings in the separately-validated *classical* `analyze` pipeline (`viability.calibrate()`'s control condition-scoping, and anti-conservative small-sample p-values in `stats.py`) — both fixed and regression-tested; see `docs/ALGORITHM_DECISIONS.md` D7/D11 and `docs/SCIENTIFIC_SPEC.md` (reclassified S2→S3) for that pipeline's own findings and current status. No scientific algorithm/threshold was changed merely to satisfy a test; changes were fixes for confirmed defects, each with its own regression test asserting the corrected (not the old, defective) behavior.

## 10. Overall verdict

For the **specified research use** (research-use quantitative morphology of registered 3D organoid hierarchies, scoped per SCIENTIFIC_SPEC.md):

> **READY FOR THE SPECIFIED RESEARCH USE WITH LIMITATIONS**

Rationale: Software/numerical correctness is established (PASS) on controlled phantoms for hierarchy, volume, topology, QC, determinism, and failure modes. Real-data biological accuracy and parameter sensitivity are **NOT ASSESSED / INSUFFICIENT EVIDENCE**, and audit separation is **Tier D / LIMITED INDEPENDENCE**. The verdict is bounded to the controlled phantom evidence and research-use scope; it must be lifted (or downgraded) as real-data V&V evidence is gathered.

This is a research-readiness assessment only. It does not establish conformity with IEC 62304, ISO 13485, ISO 14971, GAMP 5, CLIA LDT, IVDR, or any SaMD regulatory pathway. If the software is intended to support clinical decisions, a separate regulatory assessment is required and has not been performed here.

## Reproducibility / baseline (Gate 11)

| Component | Version |
|---|---|
| Git commit | `e62e052` (verified present via `git log`; supersedes the earlier `3c9793c` value, which superseded the historical `280bd06` value that was never a real commit in this repository -- see F-4) |
| Platform | macOS Apple Silicon (osx-arm64), pixi environment |
| Python | 3.12 (pixi `python = "3.12.*"`) |
| NumPy | 2.5.2 |
| SciPy | 1.17.1 |
| scikit-image | 0.26.0 |
| tifffile | 2026.3.3 |
| pandas | 2.3.3 |
| Matplotlib | 3.11.1 |
| PyYAML | 6.0.3 |

Environment is pinned via `pixi.lock`. The multilevel measurement step is deterministic (no RNG); determinism re-verified in this pass (VR-9). Package versions above were read directly from the installed pixi environment (`pixi run python -c "import numpy, ..."`) at the commit cited, not carried forward from a prior report.

## Traceability matrix

| Deliverable | Path |
|---|---|
| Scientific specification | `docs/SCIENTIFIC_SPEC.md` |
| Algorithm decisions | `docs/ALGORITHM_DECISIONS.md` |
| Parameters catalog | `docs/PARAMETERS.md` |
| Validation plan | `docs/VALIDATION_PLAN.md` |
| This report | `docs/VALIDATION_REPORT.md` |
| Classical validation record | `docs/VALIDATION.md` |
| Tests | `tests/workflows/test_multilevel_measurement_workflow.py` |
