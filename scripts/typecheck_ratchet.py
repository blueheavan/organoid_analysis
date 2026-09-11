"""Ratchet around the existing mypy debt: no new errors, reductions are locked in.

    pixi run typecheck-ratchet            # CI: fail on any new error signature
    pixi run typecheck-ratchet --update   # after fixing errors: shrink the baseline

Full `mypy src/organoid_analysis` runs with the unchanged project configuration.
Each error is reduced to a signature ``path | code | message`` (line numbers
dropped, so unrelated edits do not churn the baseline) and compared with
``scripts/typecheck_debt_baseline.json`` as a MULTISET: removing one error and
adding a different one is a new error, not unchanged debt. The job fails when

* any signature occurs more often than in the baseline (new error), or
* the baseline lists errors that no longer occur (debt was paid down but not
  locked in; run ``--update`` so it cannot silently return).

``--update`` only ever removes signatures; it refuses while new errors exist.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT / "scripts" / "typecheck_debt_baseline.json"
TARGET = "src/organoid_analysis"
SCHEMA = "organoid-analysis/typecheck-debt-baseline/1"

ERROR_LINE = re.compile(r"^(?P<path>[^:\n]+):(?P<line>\d+): error: (?P<message>.*?)  \[(?P<code>[a-z0-9-]+)\]$")
SUMMARY_LINE = re.compile(r"^Found (?P<n>\d+) errors? in \d+ files?")

# Classification is rule-based and deterministic so that it can be recomputed
# and reviewed. Order matters: the first matching rule wins.
STUBS_MISSING = "missing_third_party_stubs"  # stub packages exist but are not installed
LIBRARY_BOUNDARY = "scientific_library_typing_boundary"  # library ships no types / its signatures disagree
INTERNAL = "internal_typing_defect"
_THIRD_PARTY_SIGNATURE = re.compile(r'TiffPage|TiffFrame|pyvista|"Plotter"|vtk|"Axes"|_Wrapped')


def normalize(path: str, code: str, message: str) -> str:
    message = re.sub(r"\bline \d+\b", "line N", message)
    return f"{Path(path).as_posix()} | {code} | {message}"


def parse_mypy_output(text: str) -> Counter[str]:
    """Return the multiset of error signatures; refuse output that does not add up."""
    signatures: Counter[str] = Counter()
    reported: int | None = None
    for line in text.splitlines():
        if match := ERROR_LINE.match(line):
            signatures[normalize(match["path"], match["code"], match["message"])] += 1
        elif match := SUMMARY_LINE.match(line):
            reported = int(match["n"])
    if reported is not None and reported != sum(signatures.values()):
        raise ValueError(f"mypy reported {reported} errors but {sum(signatures.values())} were parsed")
    return signatures


def classify(signature: str) -> str:
    _, code, message = signature.split(" | ", 2)
    if code == "import-untyped" and message.startswith("Library stubs not installed"):
        return STUBS_MISSING
    if code == "import-untyped" or code == "no-any-return" or _THIRD_PARTY_SIGNATURE.search(message):
        # no-any-return here is an Any leaking out of an untyped scientific
        # library call (scipy.ndimage, sklearn, ...) through a typed function.
        return LIBRARY_BOUNDARY
    return INTERNAL


def compare(baseline: Counter[str], current: Counter[str]) -> tuple[Counter[str], Counter[str]]:
    """Return (new, fixed): signatures above and below their baseline multiplicity."""
    return current - baseline, baseline - current


def baseline_document(signatures: Counter[str], mypy_version: str) -> dict[str, object]:
    classes = Counter({name: 0 for name in (INTERNAL, STUBS_MISSING, LIBRARY_BOUNDARY)})
    for signature, count in signatures.items():
        classes[classify(signature)] += count
    return {
        "schema": SCHEMA,
        "command": f"mypy {TARGET}",
        "mypy_version": mypy_version,
        "total_errors": sum(signatures.values()),
        "class_counts": dict(classes),
        "signatures": [{"signature": s, "count": signatures[s], "class": classify(s)} for s in sorted(signatures)],
    }


def load_baseline(path: Path) -> Counter[str]:
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema") != SCHEMA:
        raise ValueError(f"{path}: unexpected schema {document.get('schema')!r}")
    return Counter({entry["signature"]: int(entry["count"]) for entry in document["signatures"]})


def run_mypy() -> tuple[str, str]:
    done = subprocess.run([sys.executable, "-m", "mypy", TARGET], cwd=ROOT, capture_output=True, text=True, check=False)
    if done.returncode not in (0, 1):
        raise RuntimeError(f"mypy did not complete (exit {done.returncode}):\n{done.stdout}\n{done.stderr}")
    version = subprocess.run([sys.executable, "-m", "mypy", "--version"], capture_output=True, text=True, check=True)
    return done.stdout, version.stdout.split()[1]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--update", action="store_true", help="shrink the baseline to the current errors (never grows it)")
    args = parser.parse_args(argv)

    output, version = run_mypy()
    print(output, end="")
    current = parse_mypy_output(output)
    if not BASELINE.is_file():
        if not args.update:
            print(f"TYPECHECK RATCHET: no baseline at {BASELINE.relative_to(ROOT)}; create it with --update")
            return 1
        BASELINE.write_text(json.dumps(baseline_document(current, version), indent=2) + "\n", encoding="utf-8")
        print(f"TYPECHECK RATCHET: baseline created with {sum(current.values())} errors")
        return 0

    baseline = load_baseline(BASELINE)
    new, fixed = compare(baseline, current)
    counts = Counter()
    for signature, count in current.items():
        counts[classify(signature)] += count
    print(f"\nTYPECHECK RATCHET (mypy {version}): baseline {sum(baseline.values())}, current {sum(current.values())}, "
          f"new {sum(new.values())}, fixed {sum(fixed.values())}")
    print("current debt by class: " + ", ".join(f"{name} {counts[name]}" for name in (INTERNAL, STUBS_MISSING, LIBRARY_BOUNDARY)))
    if new:
        print("\nNEW mypy errors (not allowed; fix them, do not add them to the baseline):")
        for signature, count in sorted(new.items()):
            print(f"  +{count}  {signature}")
        return 1
    if fixed:
        if args.update:
            BASELINE.write_text(json.dumps(baseline_document(current, version), indent=2) + "\n", encoding="utf-8")
            print(f"baseline reduced by {sum(fixed.values())} errors; commit {BASELINE.relative_to(ROOT)}")
            return 0
        print("\nDebt was reduced but the baseline still allows it. Lock the reduction in with "
              "`pixi run typecheck-ratchet --update` and commit the baseline:")
        for signature, count in sorted(fixed.items()):
            print(f"  -{count}  {signature}")
        return 1
    print("no new errors; baseline unchanged")
    return 0


if __name__ == "__main__":
    sys.exit(main())
