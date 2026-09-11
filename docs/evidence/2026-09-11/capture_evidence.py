"""Reproduce scoped probes and identify the local tested source state.

Run from the repository root with Pixi. Real TIFFs and saved masks are optional,
untracked local runtime fixtures, never independent biological ground truth.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import inspect
import json
import platform
import subprocess
import sys
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
from scipy.ndimage import zoom

from organoid_analysis.microscopy_io.tiff_contract import sha256
from organoid_analysis.quantification.features import geometry
from organoid_analysis.quantification.mask_features import extract_mask_features
from organoid_analysis.result_export.mask_feature_bundle import build_mask_feature_bundle
from organoid_analysis.segmentation.cellpose_inference import (
    config_spacing,
    read_stack_multichannel,
    restore_saved_result,
)

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent


def write_json(name: str, record: dict) -> None:
    (OUT / name).write_text(json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def probes() -> None:
    from organoid_analysis.web_interface import segmentation_workspace as workspace

    z, y, x = np.indices((25, 49, 49))
    radius = 18.0
    sphere = ((z - 12) * 2) ** 2 + (y - 24) ** 2 + (x - 24) ** 2 <= radius**2
    values, _ = geometry(sphere, (2., 1., 1.))
    volume_error = abs(values["volume_um3"] / (4 * np.pi * radius**3 / 3) - 1)
    area_error = abs(values["surface_area_um2"] / (4 * np.pi * radius**2) - 1)
    write_json("analytical-geometry.json", {
        "oracle": "Continuous radius-18 um sphere rasterized on spacing ZYX=(2,1,1) um",
        "volume_relative_error": volume_error, "area_relative_error": area_error,
        "volume_limit": .01, "surface_limit": .05,
        "volume_status": "PASS" if volume_error < .01 else "FAIL",
        "surface_status": "PASS" if area_error < .05 else "FAIL",
        "measurements": values,
    })
    print(f"Sphere: volume error {volume_error:.8%}; surface error {area_error:.8%}", flush=True)

    real_path = ROOT / "data/images/Human-Colon-Organoids-C1.tif"
    if real_path.is_file():
        stack = read_stack_multichannel(real_path)
        assert stack.spacing.complete
        spacing = (stack.spacing.x, stack.spacing.y, stack.spacing.z)
        before = hashlib.sha256(memoryview(np.ascontiguousarray(stack.volume)).cast("B")).hexdigest()
        old_budget = workspace._VIEWER_PAYLOAD_BUDGET_BYTES
        try:
            workspace._VIEWER_PAYLOAD_BUDGET_BYTES = 1_000_000
            images, _, preview_spacing = workspace._fit_viewer_payload_budget([stack.volume], None, spacing)
        finally:
            workspace._VIEWER_PAYLOAD_BUDGET_BYTES = old_budget
        source_extent = (np.array(stack.volume.shape[::-1]) - 1) * spacing
        preview_extent = (np.array(images[0].shape[::-1]) - 1) * preview_spacing
        np.testing.assert_allclose(preview_extent, source_extent, rtol=1e-12)
        after = hashlib.sha256(memoryview(np.ascontiguousarray(stack.volume)).cast("B")).hexdigest()
        assert before == after
        write_json("real-preview.json", {
            "path": str(real_path.relative_to(ROOT)), "file_sha256": sha256(real_path),
            "shape_zyx": list(stack.volume.shape), "spacing_xyz_um": spacing,
            "preview_shape_zyx": list(images[0].shape), "preview_spacing_xyz_um": preview_spacing,
            "center_extent_xyz_um": source_extent.tolist(), "array_budget_bytes": 1_000_000,
            "estimated_array_payload_bytes": int(images[0].size * 4 * 1.34),
            "source_array_sha256": before, "source_unchanged": True,
            "runtime_status": "PASS", "segmentation_accuracy": "NOT ASSESSED",
        })
        print(f"Real TIFF preview: {stack.volume.shape} -> {images[0].shape}; source unchanged", flush=True)
    else:
        write_json("real-preview.json", {"status": "NOT ASSESSED", "reason": "Local TIFF unavailable"})

    run = ROOT / "results/segmentation_output/run-20260910-092826-672784"
    if run.is_dir():
        restored = restore_saved_result(run)
        spacing = config_spacing(restored.config)
        table = extract_mask_features(restored.nuclei_masks, spacing)
        raw_volume = np.count_nonzero(restored.nuclei_masks) * np.prod(spacing)
        assert np.isclose(table.volume_um3.sum(), raw_volume, rtol=1e-12)
        provenance = json.loads((run / "provenance.json").read_text(encoding="utf-8"))
        bundle = build_mask_feature_bundle(
            restored.nuclei_masks, table, spacing_um=spacing, object_type="nucleus",
            fill_holes=False, segmentation_config=asdict(restored.config),
            segmentation_provenance=provenance,
        )
        (OUT / "saved-run-nucleus-features.zip").write_bytes(bundle)
        write_json("saved-run.json", {
            "run": str(run.relative_to(ROOT)), "mask_sha256": sha256(restored.result.nuclei_mask_path),
            "shape_zyx": list(restored.nuclei_masks.shape), "n_objects": len(table),
            "n_review_objects": int(table.qc_status.eq("review").sum()),
            "summed_raw_volume_um3": float(table.volume_um3.sum()), "counted_raw_volume_um3": float(raw_volume),
            "bundle_sha256": hashlib.sha256(bundle).hexdigest(), "runtime_status": "PASS",
            "input_domain": "Existing saved segmentation; biological or synthetic origin not independently qualified",
            "biological_accuracy": "NOT ASSESSED",
        })
        print(f"Saved run: {len(table)} objects; raw volume sum {raw_volume:g} um^3; bundle exported", flush=True)
    else:
        write_json("saved-run.json", {"status": "NOT ASSESSED", "reason": "Local saved run unavailable"})
    (OUT / "installed-resampling.txt").write_text(
        f"scipy={importlib.metadata.version('scipy')}\n" + inspect.getsource(zoom), encoding="utf-8"
    )


def snapshot() -> None:
    patch = subprocess.check_output(["git", "diff", "--binary", "HEAD"], cwd=ROOT)
    (OUT / "tracked.patch").write_bytes(patch)
    status = subprocess.check_output(["git", "status", "--porcelain=v1"], cwd=ROOT, text=True)
    paths = [ROOT / "pyproject.toml", ROOT / "pixi.lock"]
    for folder in ("src", "tests", "scripts", "configs"):
        paths.extend(p for p in (ROOT / folder).rglob("*") if p.is_file() and "__pycache__" not in p.parts
                     and (p.suffix in {".py", ".toml", ".yaml", ".yml", ".sh"} or "static" in p.parts))
    source = {str(p.relative_to(ROOT)): sha256(p) for p in sorted(set(paths))}
    untracked = {}
    for name in subprocess.check_output(["git", "ls-files", "--others", "--exclude-standard"], cwd=ROOT, text=True).splitlines():
        path = ROOT / name
        if path.is_file() and path not in {OUT / "snapshot.json", OUT / "tracked.patch"}:
            untracked[name] = sha256(path)
    write_json("snapshot.json", {
        "captured_utc": datetime.now(UTC).isoformat(),
        "HEAD": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "dirty": bool(status.strip()), "git_status": status.splitlines(),
        "tracked_patch_sha256": hashlib.sha256(patch).hexdigest(),
        "source_files_sha256": source,
        "source_manifest_sha256": hashlib.sha256(json.dumps(source, sort_keys=True).encode()).hexdigest(),
        "relevant_untracked_files_sha256": untracked,
        "hash_exclusions": ["snapshot.json and tracked.patch excluded from untracked hashes to avoid self-reference"],
        "environment": {"python": sys.version, "platform": platform.platform(), "pixi_lock_sha256": sha256(ROOT / "pixi.lock")},
        "packages": {name: importlib.metadata.version(name) for name in ("numpy", "scipy", "scikit-image", "pandas", "tifffile", "streamlit")},
        "distribution_artifact": "NOT APPLICABLE: no software release built or published in this scoped change",
        "audit_separation": "Tier D; LIMITED INDEPENDENCE",
    })
    print("Source/evidence snapshot captured", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--snapshot-only", action="store_true")
    args = parser.parse_args()
    if not args.snapshot_only:
        probes()
    snapshot()
