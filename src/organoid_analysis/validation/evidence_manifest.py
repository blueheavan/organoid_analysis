"""Tamper-evident evidence manifests for scientific validation records.

A validation record binds a scientific result to every input that could have
changed it: the source files of its declared dependency scope, the frozen V&V
plan and implementation, the environment lockfile and package versions, and
the raw and derived result artifacts. Verification re-hashes each declared
input; any difference rejects the record.

The content seal detects edits that were not re-sealed. It is not a signature
and does not authenticate anyone; binding the recorded hashes to a real Git
commit (``verify_git_provenance``) and re-executing the production code are the
stronger checks. See docs/VALIDATION_RECORDS.md.

All evidence JSON is strict (RFC 8259): no NaN or Infinity.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import importlib.util
import json
import math
import subprocess
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any, NoReturn

SCHEMA_VERSION = "organoid-analysis/evidence-manifest/1"
SEAL_FIELD = "manifest_content_sha256"
NON_FINITE_SUFFIX = "_non_finite"


class EvidenceError(ValueError):
    """Evidence that cannot be read or used as recorded."""


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def _reject_constant(token: str) -> NoReturn:
    raise EvidenceError(f"non-strict JSON constant {token!r}: evidence JSON must not contain NaN or Infinity")


def loads_strict(text: str) -> Any:
    """Parse JSON, rejecting the NaN/Infinity extensions that ``json`` accepts by default."""
    return json.loads(text, parse_constant=_reject_constant)


def dumps_strict(value: object) -> str:
    """Serialize as strict JSON; a non-finite float raises instead of writing NaN."""
    return json.dumps(value, indent=2, allow_nan=False) + "\n"


def _non_finite_label(value: float) -> str:
    return "NaN" if math.isnan(value) else ("Infinity" if value > 0 else "-Infinity")


def _is_non_finite(value: object) -> bool:
    return isinstance(value, float) and not math.isfinite(value)


def sanitize_non_finite(value: Any, pointer: str = "") -> tuple[Any, list[tuple[str, str]]]:
    """Return a strict-JSON copy of ``value`` and the ``(JSON pointer, original)`` replacements.

    A non-finite float becomes ``null``. Inside an object the key also gains a
    sibling ``<key>_non_finite`` holding ``"NaN"``, ``"Infinity"`` or
    ``"-Infinity"``, so the conversion is lossless and a reader sees that the
    quantity was non-estimable rather than absent.
    """
    if _is_non_finite(value):
        return None, [(pointer, _non_finite_label(value))]
    if isinstance(value, dict):
        out: dict[str, Any] = {}
        replaced: list[tuple[str, str]] = []
        for key, item in value.items():
            child = f"{pointer}/{str(key).replace('~', '~0').replace('/', '~1')}"
            if _is_non_finite(item):
                marker = f"{key}{NON_FINITE_SUFFIX}"
                if marker in value:
                    raise EvidenceError(f"cannot mark {child}: key {marker!r} already exists")
                out[key] = None
                out[marker] = _non_finite_label(item)
                replaced.append((child, out[marker]))
            else:
                out[key], nested = sanitize_non_finite(item, child)
                replaced += nested
        return out, replaced
    if isinstance(value, list):
        items: list[Any] = []
        replaced = []
        for index, item in enumerate(value):
            clean, nested = sanitize_non_finite(item, f"{pointer}/{index}")
            items.append(clean)
            replaced += nested
        return items, replaced
    return value, []


def content_sha256(manifest: Mapping[str, Any]) -> str:
    """SHA-256 of the canonical strict-JSON form of every field except the seal."""
    body = {key: value for key, value in manifest.items() if key != SEAL_FIELD}
    canonical = json.dumps(body, sort_keys=True, separators=(",", ":"), allow_nan=False, ensure_ascii=True)
    return sha256_bytes(canonical.encode("ascii"))


def seal(manifest: Mapping[str, Any]) -> dict[str, Any]:
    sealed = {key: value for key, value in manifest.items() if key != SEAL_FIELD}
    sealed[SEAL_FIELD] = content_sha256(sealed)
    return sealed


def verify_seal(manifest: Mapping[str, Any]) -> list[str]:
    recorded = manifest.get(SEAL_FIELD)
    if not isinstance(recorded, str):
        return ["manifest: content seal missing"]
    if recorded != content_sha256(manifest):
        return ["manifest: content was edited after sealing (content seal mismatch)"]
    return []


def load_manifest(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise EvidenceError(f"evidence manifest missing: {path}")
    try:
        document = loads_strict(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise EvidenceError(f"evidence manifest is not valid JSON: {error}") from error
    if not isinstance(document, dict):
        raise EvidenceError("evidence manifest must be a JSON object")
    return document


def _safe_relpath(relpath: object) -> str:
    if not isinstance(relpath, str) or not relpath:
        raise EvidenceError(f"artifact path must be a non-empty string, got {relpath!r}")
    pure = PurePosixPath(relpath)
    if pure.is_absolute() or ".." in pure.parts:
        raise EvidenceError(f"artifact path must be repository-relative without '..': {relpath!r}")
    return relpath


def file_record(root: Path, relpath: str, **extra: Any) -> dict[str, Any]:
    return {"path": _safe_relpath(relpath), "sha256": sha256_file(root / relpath), **extra}


def verify_file_records(root: Path, records: object, label: str) -> list[str]:
    """Check that every recorded artifact exists with exactly its recorded content."""
    if not isinstance(records, list) or not records:
        return [f"{label}: no artifacts recorded"]
    problems = []
    for record in records:
        if not isinstance(record, dict) or not isinstance(record.get("sha256"), str):
            problems.append(f"{label}: malformed artifact record {record!r}")
            continue
        try:
            relpath = _safe_relpath(record.get("path"))
        except EvidenceError as error:
            problems.append(f"{label}: {error}")
            continue
        path = root / relpath
        if not path.is_file():
            problems.append(f"{label}: required artifact missing: {relpath}")
        elif sha256_file(path) != record["sha256"]:
            problems.append(f"{label}: {relpath} changed after the validation run (SHA-256 mismatch)")
    return problems


def recorded_paths(records: object) -> set[str]:
    if not isinstance(records, list):
        return set()
    return {str(record.get("path")) for record in records if isinstance(record, dict)}


def installed_versions(names: Iterable[str]) -> dict[str, str | None]:
    versions: dict[str, str | None] = {}
    for name in names:
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    return versions


def module_origin(relpath: str) -> Path | None:
    """Return the file Python actually imports for a ``src/`` module path, if importable."""
    pure = PurePosixPath(relpath)
    if pure.parts[:1] != ("src",) or pure.suffix != ".py":
        return None
    parts = list(pure.with_suffix("").parts[1:])
    if parts[-1] == "__init__":
        parts.pop()
    spec = importlib.util.find_spec(".".join(parts))
    return Path(spec.origin) if spec is not None and spec.origin else None


# ----------------------------------------------------------------------- git
def _git(root: Path, *args: str) -> bytes | None:
    try:
        done = subprocess.run(["git", *args], cwd=root, capture_output=True, check=False)
    except OSError:
        return None
    return done.stdout if done.returncode == 0 else None


def _git_text(root: Path, *args: str) -> str | None:
    output = _git(root, *args)
    return output.decode("utf-8").strip() if output is not None else None


@dataclass(frozen=True)
class GitSourceState:
    head: str
    tree: str
    dirty_paths: tuple[str, ...]


def git_source_state(root: Path) -> GitSourceState:
    head = _git_text(root, "rev-parse", "--verify", "HEAD")
    tree = _git_text(root, "rev-parse", "--verify", "HEAD^{tree}")
    status = _git(root, "status", "--porcelain=v1", "-z", "--untracked-files=all")
    if head is None or tree is None or status is None:
        raise EvidenceError(f"{root} is not a readable Git work tree")
    dirty: list[str] = []
    entries = iter(status.decode("utf-8").split("\0"))
    for entry in entries:
        if not entry:
            continue
        dirty.append(entry[3:])
        if entry[0] in "RC":  # renames/copies carry the original path as the next entry
            dirty.append(next(entries, ""))
    return GitSourceState(head=head, tree=tree, dirty_paths=tuple(sorted(set(dirty))))


def git_commit_tree(root: Path, commit: str) -> str | None:
    return _git_text(root, "rev-parse", "--verify", "--quiet", f"{commit}^{{tree}}")


def git_blob_sha256(root: Path, commit: str, relpath: str) -> str | None:
    blob = _git(root, "cat-file", "blob", f"{commit}:{relpath}")
    return sha256_bytes(blob) if blob is not None else None


def verify_git_provenance(root: Path, source: object, scope_records: object) -> list[str]:
    """Bind the recorded scope hashes to the claimed, existing source commit."""
    if not isinstance(source, dict):
        return ["source: provenance block missing"]
    commit, tree = source.get("validated_source_commit"), source.get("validated_source_tree")
    if not isinstance(commit, str) or not isinstance(tree, str):
        return ["source: validated_source_commit / validated_source_tree missing"]
    actual_tree = git_commit_tree(root, commit)
    if actual_tree is None:
        return [f"source: validated_source_commit {commit} is not in this repository "
                "(fabricated, or a shallow clone); provenance cannot be verified"]
    problems = []
    if actual_tree != tree:
        problems.append(f"source: recorded tree {tree} is not the tree of commit {commit} ({actual_tree})")
    for record in scope_records if isinstance(scope_records, list) else []:
        if not isinstance(record, dict):
            continue
        relpath = str(record.get("path"))
        blob = git_blob_sha256(root, commit, relpath)
        if blob is None:
            problems.append(f"source: {relpath} does not exist at validated_source_commit {commit[:12]}")
        elif blob != record.get("sha256"):
            problems.append(f"source: recorded hash of {relpath} is not its content at validated_source_commit {commit[:12]}")
    return problems
