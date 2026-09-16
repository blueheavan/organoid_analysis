# Validation update — 2026-09-16 (gap-closure pass)

This pass closes defects that can be resolved from existing evidence and adds
tested study *infrastructure*. It creates **no study result**, adopts **no
threshold**, fills **no `[OWNER]` placeholder**, and moves **no gate status**.
The overall verdict remains `NOT READY FOR THE SPECIFIED RESEARCH USE`; the
science gate still reports 0/9 qualification items PASS.

> **Superseded in part, same day.** A subsequent pass recorded all fourteen
> owner decisions (D-1 … D-14) in `docs/OWNER_DECISIONS.md` on the literature
> basis in `docs/LITERATURE_BASIS_FOR_DECISIONS.md`, so the `[OWNER]`
> placeholder for the indeterminate ceiling is now filled and the acceptance
> criteria for SG-3, SG-4B, SG-5 and SG-6 are predeclared. That pass still
> created no study result and moved no gate status — the gate remains `NOT
> PASSED (0/9)` — and the statements in this document remain accurate for the
> pass they describe.

## 1. Evidence classes kept separate

| Class | What this pass produced |
|---|---|
| Engineering verification | 26 new tests (agreement metrics, spacing-ratio workflow, accelerator detection, inference-scope guard); the unsafe OpenMP override removed and guarded; runtime diagnostics script |
| Analytical characterization | nothing new; the volume/surface characterization is unchanged |
| Qualification evidence | **none.** No predeclared, independent, one-time confirmation record was produced |
| Scientific conclusions | **none.** SG-1a…SG-6B retain their prior statuses |

## 2. Confirmed defects fixed (Phase 0)

| ID | Defect | Fix |
|---|---|---|
| F1 | `ROADMAP.md` §5 and `SCIENTIFIC_VALIDATION_MASTER_PLAN.md` §11 routed SG-6B (absolute calibration) into real-object morphology, overstating the SG-3 dependency; the SG-3 bias metric is a same-grid ratio, so absolute scale cancels | E-6 in `evidence/2026-09-13-volume-audit/ERRATA.md`; both dependency diagrams corrected to `SG-6A → SG-3` and `SG-6B → absolute units` |
| F2 | Open owner decisions were scattered across six documents | consolidated into `docs/OWNER_DECISIONS.md` (D-1…D-14); `ROADMAP.md` §6 now points there; no value decided |
| F3 | The SG-4B reference circularity (a "nuclear dead-stain" reference whose only implementation is the unvalidated production nuclear path) was implicit | `INTENDED_USE_AND_ESTIMANDS.md` §5.5 Q1 now names the independence requirement explicitly and routes the choice to D-8 |
| F4 | The bibliography had no access-depth statement and the external report's citations were unresolved (review B1) | `SCIENTIFIC_VALIDATION_MASTER_PLAN.md` §13 records the traceability rule; the external report is explicitly not filed as evidence |
| F5 | SG-3's downstream-bias threshold is inherited from the area estimator (review M6 / R-3) | recorded as open decision D-6 rather than silently reused |
| F6 | Docs said the viability denominator migration was pending, but the code already uses the frozen `N_classifiable` denominator and a contract test asserts it (stale docs) | `INTENDED_USE_AND_ESTIMANDS.md` §5.3 and master plan §9 corrected to match `statistics/aggregation.py:19, 28-31` |

## 3. Confirmed defects fixed (Phase 4, engineering)

| ID | Defect | Fix |
|---|---|---|
| F7 | `KMP_DUPLICATE_LIB_OK=TRUE` was set unconditionally in both Web entry points — an OpenMP-documented unsafe override that can silently produce incorrect results | removed from `segmentation_workspace.py` and `organoid_workspace.py`; the import-order requirement is documented at `cellpose_inference.get_accelerator`; a source-level regression test forbids reintroduction |
| F8 | `detect_torch_acceleration` had no test | `tests/segmentation/test_accelerator_detection.py` covers MPS > CUDA > CPU precedence and a torch build without `backends.mps` |
| F9 | The environment's OpenMP conflicts and device selection were not observable | `scripts/runtime_diagnostics.py` reports the resolved accelerator and every OpenMP runtime on disk |

Root cause of F7/F9, verified on this host: the environment contains **five**
OpenMP runtime files (conda `libomp`/`libgomp`/`libiomp5`, plus bundled copies
under `sklearn/.dylibs` and `torch/lib`). A bare `import torch` aborts with
OpenMP error #15 (exit 134); importing the scientific stack first (`numpy`,
`scipy`, `sklearn`) is safe and MPS is then selected. The Web entry points
import `numpy`/`streamlit` before any torch import, so the override was
unnecessary. The public-data trial's reported CPU selection could **not** be
reproduced here (current `detect_torch_acceleration()` returns `mps`); that
remainder is recorded as an unresolved environment question, not explained away.

**End-to-end engineering check after the fix:** `pixi run smoke` (which runs
the real Cellpose path on `data/images/PDAC-C1/-C2.tif`, 3 × 256 × 256) now
completes with the override unset and selects MPS: `Accelerator: mps`, 29
nuclei / 23 cell masks. This is an **engineering capability** result only. It
does not qualify SG-3 (no reference labels), SG-4 (no biological reference) or
SG-6 (metadata scale not independently verified), and it is not a scientific
claim.

The **durable** fix — one OpenMP provider by taking PyTorch/torchvision from
conda-forge — was **not** applied: it changes `pixi.lock`, which is inside the
analytical-geometry evidence's dependency scope and would invalidate the
canonical record. It requires an owner-authorised environment change plus a
fresh `pixi run validation-record`, so it is reported, not performed.

## 4. Study infrastructure added (Phases 1–2), no results

| Module | Purpose |
|---|---|
| `validation/agreement_metrics.py` | SG-4B per-object metrics: confusion matrix, per-class sensitivity/precision with Wilson intervals, macro-F1, quadratic-weighted Cohen's kappa, inter-rater pairwise kappa. No Spearman for a per-object claim; no thresholds |
| `validation/spacing_calibration.py` | SG-6A: per-session X:Y and Z:XY ratios, between-session mean/SD/CV and a normal-approximation interval. `stability_verdict` returns `NOT ASSESSED` unless a tolerance is supplied; no default tolerance exists |

Existing infrastructure confirmed present and unchanged: topology diagnostics
(`filled_void_voxels`, `filled_void_components`, `open_cavity_suspected`,
`topology_flags` in `quantification/features.py`), Hungarian instance matching
(`validation/segmentation_metrics.py`), and the partial-spacing provenance
(`spacing_source_by_axis`).

## 5. Statistical inference audit (Phase 3)

`docs/STATISTICAL_INFERENCE_AUDIT.md` reconstructs the shipped path and records
that **no primary method is frozen**, so the shipped branch is not the
qualification object. No method was chosen. A regression guard keeps
`exploration.py` out of the reported inference path. SG-5 stays
`INSUFFICIENT EVIDENCE`.

## 6. SG-1a record regeneration — performed 2026-09-16

The original text of this section claimed the regeneration was "blocked, not
performed" and that a commit was needed first. **That reasoning was wrong on the
blocker.** `scripts/record_analytical_geometry_validation.py:68-71` refuses only
when an *evidence-critical* file differs from HEAD; a dirty tree elsewhere is
permitted, and `VALIDATION_RECORDS.md` §1 grants such a run the `scope-clean`
class, which policy allows to support FAIL / INSUFFICIENT EVIDENCE. Crucially,
`VALIDATION_RECORDS.md` §5 names the remedy for `EVIDENCE REJECTED` as repeating
`pixi run validation-record`, which needs no commit. Only a `canonical`
(PASS-capable) record requires a whole-tree-clean run.

Performed: `pixi run validation-record --run-id 2026-09-16-analytical-geometry-record`.
No evidence-critical file was dirty, so the builder proceeded. The new record
verifies with zero problems at HEAD `87d4df3e` and is class `scope-clean`. Its
per-case measurement values are identical to the superseded record — the
`SG-1`/`SG-2` result objects are equal and the only differing raw column is the
wall-clock `runtime_s` — and the run reproduces its own raw CSV to a maximum
relative difference of 4.6e-10 with 0 mismatches (SG-1 max |error| 25.81%,
SG-2 max 18.26%, both `FAIL`). The pointer now names it; the 2026-09-13
record is retained unmodified.

Consequence for the analytical items: SG-1a/SG-2a now report `INSUFFICIENT
EVIDENCE` instead of `FAIL` (evidence rejected). The evidence chain is current
again; what is still missing is a **domain-restricted qualification record**
(D-14), not a valid record. No gate status changes: the gate remains `NOT PASSED
(0/9)`.

## 7. Verification

- `ruff check src tests scripts`: PASS.
- `pytest`: **537 passed, 8 skipped, 8 deselected** (baseline 511 passed; +26
  new tests, 0 regressions).
- `typecheck_ratchet.py`: baseline 110, current 110, new 0, fixed 0.
- `scripts/scientific_validation_gate.py`: `GATE: NOT PASSED (0/9)`. Same
  counted verdict as the baseline, with one status change: SG-1a/SG-2a report
  `INSUFFICIENT EVIDENCE` rather than `FAIL (EVIDENCE REJECTED)` once the record
  was regenerated (§6). No item became PASS.
- `scripts/runtime_diagnostics.py`: MPS selected; five OpenMP runtimes present;
  override unset.
- `pixi run smoke` (`scripts/smoke_segmentation.py`): PASS on MPS with the
  override unset (29 nuclei / 23 cell masks); engineering capability only.
- `scripts/validate_notebooks.py`: PASS (no notebooks).
- Packaging: the package imports from the editable `src/` tree; `pip` is not
  installed in the pixi environment, so `pip check`/wheel build is not
  applicable here.

## 8. What this pass explicitly did not do

- No confirmation run, no phantom domain frozen, no acceptance threshold set.
- No algorithm, estimator, threshold, intended-use or validation claim changed.
- No owner decision made.
- No dependency or environment change (the `pixi.lock`-bound conda-forge
  PyTorch migration remains open).
