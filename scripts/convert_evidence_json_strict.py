"""Re-serialize non-strict evidence JSON (NaN/Infinity) as strict JSON, recording each change.

    pixi run python scripts/convert_evidence_json_strict.py

Historical evidence written with Python's default ``json.dumps`` may contain the
non-standard tokens NaN and Infinity, which strict parsers reject. Each such
file under docs/evidence is rewritten losslessly (``null`` plus a sibling
``<key>_non_finite`` marker, see ``evidence_manifest.sanitize_non_finite``). Its
original SHA-256, the commit from which the original bytes can be retrieved,
and every replaced JSON pointer are appended to
docs/evidence/STRICT_JSON_CONVERSION.json, so historical snapshot hashes remain
checkable. Files that are already strict are left untouched.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from organoid_analysis.validation.evidence_manifest import (
    EvidenceError,
    dumps_strict,
    git_blob_sha256,
    loads_strict,
    sanitize_non_finite,
    sha256_file,
)

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "docs" / "evidence"
RECORD = EVIDENCE / "STRICT_JSON_CONVERSION.json"


def main() -> int:
    record = loads_strict(RECORD.read_text(encoding="utf-8")) if RECORD.is_file() else {
        "schema": "organoid-analysis/strict-json-conversion/1",
        "policy": "Evidence JSON is strict RFC 8259 JSON. Non-finite numbers are written as null with a sibling "
                  "'<key>_non_finite' marker ('NaN', 'Infinity' or '-Infinity'); list elements are recorded here.",
        "conversions": []}
    converted = 0
    for path in sorted(EVIDENCE.rglob("*.json")):
        text = path.read_text(encoding="utf-8")
        try:
            loads_strict(text)
            continue
        except EvidenceError:
            pass
        rel = path.relative_to(ROOT).as_posix()
        original = sha256_file(path)
        last = subprocess.run(["git", "log", "-1", "--format=%H", "--", rel], cwd=ROOT,
                              capture_output=True, text=True, check=True).stdout.strip()
        clean, replaced = sanitize_non_finite(json.loads(text))
        path.write_text(dumps_strict(clean), encoding="utf-8")
        record["conversions"].append({
            "path": rel, "original_sha256": original,
            "original_retrievable_with": f"git show {last}:{rel}",
            "original_matches_that_commit": git_blob_sha256(ROOT, last, rel) == original,
            "converted_sha256": sha256_file(path),
            "replacements": [{"pointer": pointer, "original": token} for pointer, token in replaced]})
        converted += 1
        print(f"converted {rel}: {len(replaced)} non-finite values")
    RECORD.write_text(dumps_strict(record), encoding="utf-8")
    print(f"{converted} files converted")
    return 0


if __name__ == "__main__":
    sys.exit(main())
