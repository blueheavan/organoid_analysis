# Validation records, evidence integrity and engineering gates

This document defines how scientific evidence is recorded, bound to the code it
validates, and checked by the scientific gate. It also describes the engineering
CI structure. It changes **no** scientific algorithm, threshold, criterion,
intended use, statistical model or measurement definition.

## 1. Source commit versus validation record

A validation result is only meaningful for the code that produced it. Two
objects are kept distinct:

| Object | What it is | Where it lives |
|---|---|---|
| **Source commit** (`validated_source_commit`, `validated_source_tree`) | An immutable commit containing the code under validation | Git history |
| **Validation record** | A sealed `evidence_manifest.json` plus the raw and derived artifacts, generated *after* the source commit | `docs/evidence/<run-id>/`, committed later (ideally in an evidence-only commit) |

A record names `validated_source_commit` but **never** the commit that stores
it, because a commit cannot contain its own hash. That avoids the circular claim
that "this commit validated itself". To find the storing commit, run
`git log --format=%H -1 -- docs/evidence/<run-id>/evidence_manifest.json`.

**Record classes**

- `canonical`: the whole working tree was clean at `validated_source_commit`
  when the V&V ran. Only a canonical record can make a gate item PASS.
- `scope-clean`: every evidence-critical file (§3) was identical to
  `validated_source_commit`, but other paths were dirty. They are listed in
  `source.dirty_paths_at_execution`, and those that affect the result
  (tooling) are hashed in `vv_implementation`. A scope-clean record can
  support FAIL / INSUFFICIENT EVIDENCE, never PASS.
- The builder refuses to write any record when an evidence-critical file
  differs from HEAD.

**Historical evidence.** Directories with legacy `snapshot.json` files
(schema v0) ran on dirty working trees and are labeled historical/pre-release
([index](evidence/README.md)). They are not validation records of any commit.

## 2. Workflow

```bash
# 1. Commit the code under validation (source commit). Branch first if on main.
git commit ...
# 2. Run the frozen V&V and write the record from the clean tree.
pixi run validation-record --run-id <YYYY-MM-DD>-analytical-geometry-record
# 3. Commit only the new record directory and the pointer (evidence-only commit).
git add docs/evidence/<run-id> docs/evidence/analytical_geometry_current_record.json
git commit -m "Validation record for <source commit>"
# 4. Verify.
pixi run science-gate
```

The builder:
1. refuses if the run directory exists (evidence is never overwritten) or if an
   evidence-critical file differs from HEAD;
2. checks the plan and harness against the freeze record written before any
   result existed;
3. copies the byte-identical frozen harness into the run directory and runs it;
4. re-serializes the harness summary as strict JSON;
5. confirms that no evidence-critical file changed during the run;
6. derives SG-1/SG-2 from the raw CSV;
7. re-executes the production estimator on every phantom;
8. writes and seals the manifest;
9. updates `docs/evidence/analytical_geometry_current_record.json`;
10. verifies the new record before reporting success.

The pointer file is deliberately not hashed by any record. Hashing it would be
circular: updating the pointer would invalidate the record it points to. It can
only choose among records that verify against the current repository.

## 3. Manifest schema (`organoid-analysis/evidence-manifest/1`)

| Field | Content |
|---|---|
| `schema_version`, `contract_id`, `validation_run_id`, `validation_scope` | Identity of the record and of its evidence contract (`analytical-geometry/1`) |
| `record_class` | `canonical` or `scope-clean` (§1) |
| `source` | `validated_source_commit`, `validated_source_tree`, `whole_tree_clean_at_execution`, `dirty_paths_at_execution`, `scope_clean_at_execution`, `record_storage` |
| `dependency_scope.files` | Path, role, rationale and SHA-256 of every evidence-critical source file |
| `dependency_scope.excluded` | What was deliberately left out, and why |
| `environment` | `pixi.lock` SHA-256, Python, platform, versions of numpy/scipy/scikit-image/pandas |
| `vv_plan` | SHA-256 of the frozen plan and of its freeze record |
| `vv_implementation` | SHA-256 of the executed harness (byte-identical to the frozen one) and of the derivation/verification/recording tooling |
| `raw_artifacts`, `derived_artifacts` | SHA-256 of the per-case CSV and run log, and of the summary JSON/MD (with the strict-JSON replacements made) |
| `estimator` | Production `SURFACE_AREA_METHOD` identity and the volume definition |
| `acceptance_criteria` | SG-1 `< 0.05`, SG-2 `< 0.01` (SCIENTIFIC_SPEC §9), per case |
| `commands`, `invocation` | Commands executed, with exit codes and times |
| `results`, `resulting_status` | Statuses and values derived from the raw CSV |
| `reproduction` | Re-execution summary (cases, tolerance, maximum difference, mismatches) |
| `manifest_content_sha256` | Seal over the canonical JSON of all other fields |

### Analytical-geometry dependency scope

These are the files whose content can change an SG-1/SG-2 value, or change what
that value means in an exported measurement
(`analytical_geometry_evidence.DEPENDENCY_SCOPE`):

| File | Role |
|---|---|
| `quantification/features.py` | Estimator core: `geometry()` volume, surface area and sphericity; `surface_mesh()` level, padding and spacing; `outer_envelope()`; `SURFACE_AREA_METHOD` |
| `quantification/mask_features.py` | Web route: raw vs filled support, and the (x,y,z)→(z,y,x) spacing conversion |
| `quantification/cellular_measurements.py` | Calls `geometry()` with the filled-envelope default for cells |
| `quantification/multilevel_relationships/morphology.py` | Calls `geometry()` on raw labels |
| `microscopy_io/voxel_spacing.py` | Spacing validation and axis order for the mask-feature route |

Excluded, with reasons:
- File-level acquisition-spacing parsing is gate item SG-6.
- Workflows only orchestrate.
- `watershed_instances`/`labels` provide types, crops and flags.
- Third-party numerics are bound by `pixi.lock`, package versions and re-execution.
- Package `__init__` modules contain no numerical code.
- The gate script only presents results.

## 4. What rejects a record

`analytical_geometry_evidence.verify_record` rejects a record on any of the
following:
- a missing or changed scope file, plan, freeze record, harness, tooling file,
  lockfile, raw artifact or derived artifact;
- a required artifact that is missing from the manifest, or a shrunken scope;
- a plan or harness that is not the frozen version;
- a different estimator identity, for example a version bump, which does not
  inherit old evidence;
- a criterion other than §9;
- different numerical package or Python versions;
- a failed validation command;
- a scope that was dirty at execution;
- a broken seal;
- results or status that do not re-derive from the raw CSV;
- a summary that disagrees with the raw derivation;
- an imported production module that is not the hashed file;
- a source commit that is not in the repository, a tree that does not match
  it, or recorded scope hashes that are not the content of that commit;
- raw values that are not reproduced by re-executing the production estimator
  and the closed-form oracles on all 188 phantoms (relative tolerance 1e-9,
  set by the CSV's 10 significant digits);
- non-strict JSON.

When the record is rejected, the gate reports SG-1 and SG-2 as FAIL with
`EVIDENCE REJECTED` until `pixi run validation-record` is repeated.

**PASS rule.** An item is PASS only when it is backed by a verified
**canonical** record whose results report PASS for that item. Items without a
record (SG-3 to SG-6) cannot be PASS.

**Limits.** The seal is not a signature and authenticates no one. Someone with
write access can always edit the gate itself. The protection against a fully
consistent forgery of hashes, seal and derivation is re-execution of the
production code, together with the binding of scope hashes to a real commit.
Git provenance needs full history: CI uses `fetch-depth: 0`. The record also
bounds numerical reproduction to this environment; a legitimate
floating-point difference on another platform would be reported, not ignored.

## 5. Current record

`docs/evidence/2026-09-11-analytical-geometry-record/`:
- class **scope-clean**;
- `validated_source_commit` `5143be4001a5f5430b1ef46cb8c887aca3ba45de`.

The uncommitted paths at execution were validation tooling, tests, CI/config
and the historical strict-JSON conversion. No evidence-critical file was among
them.

It is not canonical because the tooling was not yet committed. Repeat §2 after
committing to obtain a canonical record.

| Item | Status | Basis |
|---|---|---|
| SG-1 surface area | **FAIL** | 188 cases; max \|error\| 18.74% (cylinder r=24 µm, spacing 3×1×1); median signed +9.72%; 0 non-estimable |
| SG-2 voxel volume | **FAIL** | fails for ρ 6–12 (max 2.57%), ρ 3–6 (2.74%), ρ<3 (18.26%); passes ρ≥12 |
| Re-execution | 0 mismatches | 188 cases, max relative difference 4.8e-10 |

The raw results reproduce the historical dirty-tree run exactly. All 1,128
(case, estimator) rows match in every field except `runtime_s`, compared as
strings against `2026-09-11-measurement-vv/surface_vv_dev.csv`.

## 6. Adversarial regression tests

`tests/validation/test_evidence_integrity.py` tampers with a byte-identical
copy of the record and asserts that the tampered input is named.

| Required case | Test |
|---|---|
| 1. Surface change without a method-version bump | `test_surface_change_without_method_version_bump_is_rejected` (plus every scope file, parametrized) |
| 2. Volume change | `test_volume_change_is_rejected` |
| 3. Plan change | `test_plan_change_is_rejected` |
| 4. V&V script change | `test_vv_harness_change_is_rejected`, `test_derivation_tooling_change_is_rejected` |
| 5. Raw result change | `test_raw_result_change_is_rejected` |
| 6. Summary-only change | `test_summary_only_change_is_rejected` |
| 7. Missing artifact | `test_removed_run_artifact_is_rejected`, `test_removed_protocol_or_environment_file_is_rejected`, `test_artifact_dropped_from_the_manifest_is_rejected` |
| 8. Fabricated PASS | `test_hand_edited_pass_breaks_the_seal_and_the_raw_derivation`, `test_relaxed_criterion_in_the_record_is_not_honored`, `test_consistently_forged_pass_is_caught_by_reexecution`, gate tests `test_pass_*` |
| 9. Change invisible to the old sphere sentinel | `test_change_invisible_to_the_old_sphere_sentinel_is_still_rejected`: the mutant reproduces the old sentinel area exactly, but raw-label support has silently become the envelope |

The same file also covers:
- lockfile changes, package versions, failed commands, a dirty scope and
  non-strict JSON;
- a shrunken scope and a missing pointer;
- Git provenance: a matching commit is accepted; a nonexistent commit, a wrong
  tree, and re-hashed source that is not the validated commit are rejected.

## 7. Strict JSON

Evidence JSON must be RFC 8259 JSON:
- It is written with `allow_nan=False` and parsed with a strict parser.
- A non-finite value becomes `null` plus a sibling `<key>_non_finite`
  (`"NaN"`, `"Infinity"` or `"-Infinity"`), next to the existing
  `n_nonestimable` and `status` fields.
- Three historical files were converted. Their original hashes and the
  `git show` commands that recover them are listed in
  `evidence/STRICT_JSON_CONVERSION.json`.

`tests/validation/test_evidence_json_strict.py` strict-parses every
`docs/evidence/**/*.json` file and checks that each conversion is lossless and
recoverable.

## 8. Engineering CI

| | Before | After |
|---|---|---|
| Jobs | One job, `pixi run ci` = lint → test → check → typecheck; always red because of mypy debt | `engineering-regression` (`pixi run regression`: lint, tests, notebook check; full Git history) and a separate `typecheck-debt` job (`pixi run typecheck-ratchet`) |
| mypy | Pass/fail on the raw count | Full `mypy src/organoid_analysis`, unchanged config, compared as a multiset of `path \| code \| message` signatures against `scripts/typecheck_debt_baseline.json` |

The ratchet:
- fails on any new signature, including "fix one, add another";
- fails when debt was reduced but not locked in;
- has an `--update` mode that only ever shrinks the baseline.

Moving code that has baselined errors to another file counts as new errors
there.

The baseline is 110 errors, classified by deterministic rule:
- 34 internal typing defects;
- 44 missing third-party stubs (pandas, scipy, PyYAML, seaborn);
- 32 scientific-library typing boundaries (statsmodels, scikit-learn, cellpose
  and mpl_toolkits without types; tifffile, pyvista, vtk and matplotlib
  signature mismatches; Any returned from untyped calls).

Mark `engineering-regression` as a required status check in branch protection;
that setting lives outside the repository.
