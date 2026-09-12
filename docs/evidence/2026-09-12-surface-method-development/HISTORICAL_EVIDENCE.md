# Class: method-development evidence (not a validation record)

This directory holds the development and confirmation work that produced the
adopted surface-area estimator `crofton_minimax_sym_v3`. It is **not** a sealed
validation record and the scientific gate does not read it. The record that
validates the adopted estimator against this repository is
[`../2026-09-12-surface-crofton-v3/`](../2026-09-12-surface-crofton-v3); the
gate reads that one.

It is kept because the record alone does not show *how* the method was
obtained, and because the qualifying claim rests on evidence produced here: the
confirmation set, the freeze that preceded it, and the constants that were
fixed before any confirmation result was seen.

## Provenance and its limits

The work ran outside this repository, in a standalone workspace, against
independently generated phantoms — not through the repo's frozen V&V harness.
Nothing here is hashed by any record and none of it is re-executed by CI. Treat
it as a lab notebook: the numbers are as produced, the code is as run, and the
freeze records state what was fixed at which point. What binds the method to
this repository is the validation record, the packaged weight tables, and the
tests.

## Two rounds, and the first one failed

- **Round 1** (`FREEZE_RECORD.md`, `SURFACE_VV_REPORT.md`, `CANDIDATES.md`,
  `MECHANISM.md`) qualified nothing. No candidate met the predeclared rule. It
  is here in full, including the mechanism analysis that explained why, because
  the round-2 protocol was designed from it and because "no qualified
  estimator" was a legitimate outcome that had to be reportable.
- **Round 2** (`ROUND2_PROTOCOL.md`, `FREEZE_RECORD_V3.md`,
  `ROUND2_REPORT.md`, `confirm2_verdict.json`) qualified
  `crofton_minimax_sym_v3` within a declared domain. The report's addendum
  records what integration into this repository later revealed, including the
  fact that the a-priori budget is predictive rather than a proven bound.

## Reading order

1. `ROUND2_REPORT.md` — the verdict, the evidence, the domain, the addendum.
2. `FREEZE_RECORD_V3.md` — what was frozen, and when, relative to the
   confirmation run.
3. `ROUND2_PROTOCOL.md` — the predeclared rules, including the PASS/FAIL rule
   and the pre-registered predictions.
4. `confirm2_verdict.json`, `vv_confirm2.csv` — the confirmation outcome and
   its per-case table.
5. `vv_dev2_all.csv` — the development table. Tuning used this; the
   confirmation table was untouched until the single confirmation run.
6. `crease_indicator_rejected.md` — a machine-checkable crease indicator that
   was tried and rejected, which is why the smoothness scope is declared in
   prose and left to the caller instead of being enforced in code.
7. `lp_symmetry_verification.csv` — the orbit-symmetry reduction checked
   against unrestricted solutions.
8. `code/` — the modules and step scripts as run. `surface_area_v3.py` and
   `estimators_v3.py` are the development ancestors of
   `src/organoid_analysis/quantification/surface_crofton.py`; the port was
   cross-checked by replaying confirmation cases through both, which reproduced
   the recorded areas exactly.

## One thing this directory cannot tell you

The weight vectors here and those packaged in
`src/organoid_analysis/quantification/crofton_weights/` are the same vectors,
not merely the same program's output: the minimax optimum is not unique, so
re-solving gives a different optimal vertex. The packaged tables are the ones
that produced the confirmation evidence. See that directory's `PROVENANCE.md`.
