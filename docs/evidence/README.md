# Evidence index

Each directory is either a **validation record** (a sealed
`evidence_manifest.json` that the scientific gate verifies) or **historical
evidence** (kept unchanged as produced; never used by the gate). See
[`docs/VALIDATION_RECORDS.md`](../VALIDATION_RECORDS.md) for the schema and workflow.

| Directory | Class | Source state | Used by the gate |
|---|---|---|---|
| `2026-09-16-analytical-geometry-record/` | Validation record, class `scope-clean` | V&V re-executed at HEAD `87d4df3e`; every evidence-critical file identical to that commit, and the 70 recorded dirty paths touch no dependency-scope file, the frozen plan/freeze/harness or `pixi.lock` (47 of them are agent scratch under `.wisp/`). Supersedes the 2026-09-13 record, whose scope hashes stopped matching the tree at commit `34b6f3c`. Per-case measurement values are identical to both earlier records — the `SG-1`/`SG-2` result objects are equal and the only differing raw column is the wall-clock `runtime_s` (max |error| 25.81% / 18.26%) — and the run reproduces its own raw CSV to a maximum relative difference of 4.6e-10 with 0 mismatches, so this is a re-execution that restores the binding between record and source under test, **not** new independent confirmation | Yes: SG-1, SG-2 (named by `analytical_geometry_current_record.json`) |
| `2026-09-13-analytical-geometry-record/` | Validation record, class `canonical` (superseded 2026-09-16) | V&V re-executed at HEAD `be223ba` with a clean working tree, after the per-axis spacing fix. Per-case measurement values are identical to `2026-09-12-surface-crofton-v3/` (the `SG-1`/`SG-2` result objects are equal; the only differing raw column is the wall-clock `runtime_s`): a re-execution that restored the binding between record and source under test, **not** new independent confirmation. Retained unmodified; its scope hashes ceased to match the tree when commit `34b6f3c` changed five `quantification/` files | No (superseded as the gate's record; still the provenance cited for its per-case numbers, which the successor reproduces exactly) |
| `2026-09-13-volume-audit/` | Analysis of existing evidence; no new measurement | Audits the volume estimator from the two records above and proposes studies. Round 4 added `ERRATA.md` (five corrections to this directory's own round-3 documents), the volume study design that follows from the estimand decision in `docs/INTENDED_USE_AND_ESTIMANDS.md`, derived calibration propagation with tiered acceptance criteria, the SG-3 addenda and the SG-5 design. Declares no domain and changes no algorithm | No (SG-2 stays non-PASS — `INSUFFICIENT EVIDENCE` for the qualification item, its 188-case unrestricted characterization being descriptive only; SG-3, SG-5 and SG-6 stay `INSUFFICIENT EVIDENCE` or `NOT ASSESSED`) |
| `2026-09-12-surface-crofton-v3/` | Validation record, class `canonical` (superseded 2026-09-13) | V&V run at HEAD `b659e2e89378` with a clean working tree; the run that re-validated the adopted surface estimator `crofton_minimax_sym_v3`. Retained unmodified; its per-case results equal those of the 2026-09-13 record | No (superseded as the gate's record; kept as the first record of the adopted estimator) |
| `2026-09-12-surface-method-development/` | Method-development evidence ([label](2026-09-12-surface-method-development/HISTORICAL_EVIDENCE.md)) | Produced outside this repository, against independently generated phantoms; not hashed by any record | No (it is how the adopted estimator was obtained and confirmed, not a validation of this repository) |
| `2026-09-11-analytical-geometry-canonical/` | Validation record, class `canonical` | Superseded 2026-09-12 by the record above: it validates the superseded surface estimator and no longer matches the repository | No (kept as the evidence that `marching_cubes_binary_lewiner_v1` fails the §9 area criterion) |
| `2026-09-11-analytical-geometry-record/` | Validation record, class `scope-clean` (superseded) | V&V run at HEAD `5143be4`; every evidence-critical file identical to that commit. The uncommitted changes (listed in the manifest) were limited to validation tooling, tests, CI/config and the historical strict-JSON conversion | No (superseded) |
| `2026-09-11-measurement-vv/` | Historical / pre-release ([label](2026-09-11-measurement-vv/HISTORICAL_EVIDENCE.md)) | Dirty working tree on `ed81d74`; code later committed as `5143be4` | No, except the frozen plan, freeze record and harness, which records verify against |
| `2026-09-11/` | Historical / pre-release | Dirty working tree (`snapshot.json`, legacy schema v0) | No |
| `2026-09-10/` | Historical / pre-release | Dirty working tree (`final-snapshot.json`, legacy schema v0) | No |
| `2026-09-09/` | Historical only | See `docs/VALIDATION_REPORT.md` | No |

Legacy `snapshot.json` files do not distinguish the source under validation
from the commit that stores the evidence. They must not be cited as validation
of a commit.

`STRICT_JSON_CONVERSION.json` lists the historical JSON files that were
re-serialized as strict JSON. It gives each original SHA-256 and the commit
from which the original bytes can be retrieved.
