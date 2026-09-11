# Historical / pre-release evidence — not a canonical validation record

This directory is preserved as produced on 2026-09-11. It is **historical,
pre-release evidence**, not a validation record of any commit:

- **Source state when it ran:** HEAD `ed81d74ac9d8c3b7b2c5efe6e9e9549d00bc7856`
  with an **uncommitted (dirty) working tree** (`snapshot.json`: `"dirty": true`).
  `snapshot.json` uses the legacy snapshot layout (schema v0), which does not
  separate the source under validation from the commit storing the evidence.
- **Relation to a commit (verified 2026-09-11, after the fact):** the tracked
  patch it recorded (`tracked.patch`, SHA-256 `23c4215f…ee2b`) is byte-identical
  to `git diff ed81d74 5143be4` over the same 23 tracked files, and every
  untracked-file hash in `snapshot.json` equals the corresponding blob of commit
  `5143be4001a5f5430b1ef46cb8c887aca3ba45de`. So the code it exercised is the
  code later committed in `5143be4`. That equivalence was established
  afterwards; the run itself did not happen at a clean commit.
- **Strict JSON:** `surface_vv_dev_summary.json` and `runtime_probe.json`
  contained NaN/Infinity and were re-serialized as strict JSON. Their original
  bytes and hashes are listed in `../STRICT_JSON_CONVERSION.json`, so the
  `snapshot.json` hashes remain checkable.
- **Still authoritative here:** the frozen protocol (`SURFACE_VV_PLAN.md`,
  `surface_vv_freeze.txt`, `surface_vv.py`). New records verify their plan and
  harness against this freeze record.

**Superseded for gate purposes by** the sealed record
[`../2026-09-11-analytical-geometry-record/evidence_manifest.json`](../2026-09-11-analytical-geometry-record/evidence_manifest.json),
which re-ran the same frozen harness from a scope-clean checkout of `5143be4`.
See [`docs/VALIDATION_RECORDS.md`](../../VALIDATION_RECORDS.md).
