# Evidence index

Each directory is either a **validation record** (a sealed
`evidence_manifest.json` that the scientific gate verifies) or **historical
evidence** (kept unchanged as produced; never used by the gate). See
[`docs/VALIDATION_RECORDS.md`](../VALIDATION_RECORDS.md) for the schema and workflow.

| Directory | Class | Source state | Used by the gate |
|---|---|---|---|
| `2026-09-12-surface-crofton-v3/` | Validation record, class `canonical` | V&V run at HEAD `b659e2e89378` with a clean working tree; the run that re-validated the adopted surface estimator `crofton_minimax_sym_v3` | Yes: SG-1, SG-2 (named by `analytical_geometry_current_record.json`) |
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
