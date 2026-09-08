from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from organoid_analysis import __version__


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="3D organoid morphology and control-scaled viability; research use only.")
    parser.add_argument("--version",action="version",version=__version__)
    sub = parser.add_subparsers(dest="command",required=True)
    run = sub.add_parser("analyze",help="Analyze TIFF stacks described by a sample manifest")
    run.add_argument("--manifest",required=True,help="Sample-design CSV; image paths are relative to this CSV")
    run.add_argument("--config",help="YAML configuration; unspecified keys use documented defaults")
    run.add_argument("--out",required=True,help="New/empty output directory")
    run.add_argument("--keep-going",action="store_true",help="Continue failed samples; return exit code 2 for partial runs")
    demo = sub.add_parser("demo",help="Generate and analyze labeled synthetic phantoms")
    demo.add_argument("--out",default="data/synthetic_demo",help="New/empty demo directory")
    demo.add_argument("--seed",type=int,default=1729)
    cells = sub.add_parser("cells", help="Pair cell/nucleus masks and measure cell-level features")
    cells.add_argument("--cell-labels", required=True, help="3D ZYX cell instance-label TIFF")
    cells.add_argument("--nucleus-labels", required=True, help="Registered 3D ZYX nucleus instance-label TIFF")
    cells.add_argument("--axes", default="", help="Axes override for both TIFFs, for example ZYX")
    cells.add_argument("--spacing-z-um", type=float)
    cells.add_argument("--spacing-y-um", type=float)
    cells.add_argument("--spacing-x-um", type=float)
    cells.add_argument("--min-cell-volume-um3", type=float, default=200.0)
    cells.add_argument("--min-nucleus-volume-um3", type=float, default=25.0)
    cells.add_argument("--max-nc-ratio", type=float, default=1.0)
    cells.add_argument("--min-nucleus-containment", type=float, default=0.5)
    cells.add_argument("--allow-nucleusless", action="store_true",
                       help="Keep cells with no nucleus overlap; explicit nucleus QC failures still reject")
    cells.add_argument("--neighbor-radius-um", type=float, default=25.0)
    cells.add_argument("--out", required=True, help="New/empty output directory")
    multilevel = sub.add_parser("analyze-3d", help="Analyze registered organoid, cell, and nucleus 3D labels")
    multilevel.add_argument("--organoid-labels", required=True, help="3D ZYX organoid instance-label TIFF")
    multilevel.add_argument("--cell-labels", required=True, help="3D ZYX cell instance-label TIFF")
    multilevel.add_argument("--nucleus-labels", required=True, help="3D ZYX nucleus instance-label TIFF")
    multilevel.add_argument("--nucleus-intensity", help="Optional registered raw nucleus-intensity TIFF")
    multilevel.add_argument("--axes", default="", help="Axes override for all supplied TIFFs, for example ZYX")
    multilevel.add_argument("--spacing-z-um", type=float)
    multilevel.add_argument("--spacing-y-um", type=float)
    multilevel.add_argument("--spacing-x-um", type=float)
    multilevel.add_argument("--sample-id")
    multilevel.add_argument("--well-id")
    multilevel.add_argument("--field-id")
    multilevel.add_argument("--batch-id")
    multilevel.add_argument("--condition")
    multilevel.add_argument("--treatment")
    multilevel.add_argument("--dose")
    multilevel.add_argument("--timepoint")
    multilevel.add_argument("--minimum-voxels", type=int, default=5)
    multilevel.add_argument("--low-parent-overlap-fraction", type=float, default=0.5)
    multilevel.add_argument("--mad-z-threshold", type=float, default=3.5)
    multilevel.add_argument("--core-max-normalized-radial-position", type=float, default=0.5)
    multilevel.add_argument("--peripheral-min-normalized-radial-position", type=float, default=0.8)
    multilevel.add_argument("--out", required=True, help="New/empty output directory")
    arguments = parser.parse_args(argv)
    try:
        if arguments.command == "demo":
            from organoid_analysis.workflows.demo_workflow import run_demo
            result = run_demo(arguments.out,arguments.seed)
        elif arguments.command == "cells":
            result = _run_cells(arguments)
        elif arguments.command == "analyze-3d":
            result = _run_multilevel(arguments)
        else:
            from organoid_analysis.workflows.organoid_measurement_workflow import analyze
            result = analyze(arguments.manifest,arguments.out,arguments.config,arguments.keep_going)
        print(json.dumps(result,indent=2))
        return 2 if result["status"] == "partial" else 0
    except Exception as error:
        print(f"ERROR: {type(error).__name__}: {error}",file=sys.stderr)
        return 1


def _run_cells(arguments) -> dict:
    from organoid_analysis.quantification.cellular_measurements import analyze_cells
    from organoid_analysis.microscopy_io.tiff_contract import SPACING_ATOL_UM, SPACING_RTOL, read_tiff, write_labels

    cell_stack, cell_spacing, _ = read_tiff(arguments.cell_labels, arguments.axes)
    nucleus_stack, nucleus_spacing, _ = read_tiff(arguments.nucleus_labels, arguments.axes)
    if cell_stack.shape[0] != 1 or nucleus_stack.shape[0] != 1:
        raise ValueError("Cell and nucleus label TIFFs must each contain exactly one channel")
    metadata_spacings = [value for value in (cell_spacing, nucleus_spacing) if value is not None]
    if len(metadata_spacings) == 2 and not np.allclose(metadata_spacings[0], metadata_spacings[1], rtol=SPACING_RTOL, atol=SPACING_ATOL_UM):
        raise ValueError("Cell and nucleus OME voxel spacings differ")

    explicit = (arguments.spacing_z_um, arguments.spacing_y_um, arguments.spacing_x_um)
    if any(value is not None for value in explicit) and not all(value is not None for value in explicit):
        raise ValueError("Provide all three spacing values or omit all three")
    if all(value is not None for value in explicit):
        spacing = explicit
        if metadata_spacings and not np.allclose(spacing, metadata_spacings[0], rtol=SPACING_RTOL, atol=SPACING_ATOL_UM):
            raise ValueError("Explicit spacing conflicts with OME metadata")
    elif metadata_spacings:
        spacing = metadata_spacings[0]
    else:
        raise ValueError("Voxel spacing is missing; provide OME metadata or all three --spacing-*-um values")

    out = Path(arguments.out).resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {out}")
    out.mkdir(parents=True, exist_ok=True)
    result = analyze_cells(
        cell_stack[0], nucleus_stack[0], spacing,
        min_cell_volume_um3=arguments.min_cell_volume_um3,
        min_nucleus_volume_um3=arguments.min_nucleus_volume_um3,
        require_nucleus=not arguments.allow_nucleusless,
        max_nc_ratio=arguments.max_nc_ratio,
        min_nucleus_containment=arguments.min_nucleus_containment,
        neighbor_radius_um=arguments.neighbor_radius_um,
    )
    write_labels(out / "cell_labels.filtered.ome.tif", result.cell_labels, spacing)
    write_labels(out / "nucleus_labels.filtered.ome.tif", result.nucleus_labels, spacing)
    result.features.to_csv(out / "cell_features.csv", index=False, float_format="%.10g")
    result.pairing.to_csv(out / "pairing_qc.csv", index=False, float_format="%.10g")
    summary = {
        **result.summary,
        "status": "complete",
        "spacing_zyx_um": list(spacing),
        "cell_labels": str(Path(arguments.cell_labels).resolve()),
        "nucleus_labels": str(Path(arguments.nucleus_labels).resolve()),
        "output": str(out),
        "neighborhood_incomplete_cells": int((~result.features["neighborhood_complete"]).sum()) if len(result.features) else 0,
    }
    (out / "cell_analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def _run_multilevel(arguments) -> dict:
    """CLI adapter: load registered labels/metadata then invoke pure analysis."""
    import importlib.metadata
    import platform

    from organoid_analysis.quantification.multilevel_relationships import Multilevel3DConfig
    from organoid_analysis.workflows.multilevel_measurement_workflow import analyze_multilevel_3d
    from organoid_analysis.result_export.measurement_tables import export_results
    from organoid_analysis.microscopy_io.tiff_contract import SPACING_ATOL_UM, SPACING_RTOL, git_commit_hash, read_tiff, sha256

    def read_single(path: str, role: str) -> tuple[np.ndarray, tuple | None]:
        stack, spacing, _ = read_tiff(path, arguments.axes)
        if stack.shape[0] != 1:
            raise ValueError(f"{role} TIFF must contain exactly one channel")
        return stack[0], spacing

    organoids, organoid_spacing = read_single(arguments.organoid_labels, "Organoid labels")
    cells, cell_spacing = read_single(arguments.cell_labels, "Cell labels")
    nuclei, nucleus_spacing = read_single(arguments.nucleus_labels, "Nucleus labels")
    intensity = None
    intensity_spacing = None
    if arguments.nucleus_intensity:
        intensity, intensity_spacing = read_single(arguments.nucleus_intensity, "Nucleus intensity")
    metadata_spacings = [item for item in (organoid_spacing, cell_spacing, nucleus_spacing, intensity_spacing) if item is not None]
    if len(metadata_spacings) > 1 and any(not np.allclose(metadata_spacings[0], item, rtol=SPACING_RTOL, atol=SPACING_ATOL_UM) for item in metadata_spacings[1:]):
        raise ValueError("Input TIFF voxel spacings differ; register/resample before multilevel analysis")
    explicit = (arguments.spacing_z_um, arguments.spacing_y_um, arguments.spacing_x_um)
    if any(value is not None for value in explicit) and not all(value is not None for value in explicit):
        raise ValueError("Provide all three spacing values or omit all three")
    if all(value is not None for value in explicit):
        spacing = explicit
        if metadata_spacings and not np.allclose(spacing, metadata_spacings[0], rtol=SPACING_RTOL, atol=SPACING_ATOL_UM):
            raise ValueError("Explicit spacing conflicts with OME metadata")
    elif metadata_spacings:
        spacing = metadata_spacings[0]
    else:
        raise ValueError("Voxel spacing is missing; provide OME metadata or all three --spacing-*-um values")
    config = Multilevel3DConfig(
        minimum_voxels=arguments.minimum_voxels,
        low_parent_overlap_fraction=arguments.low_parent_overlap_fraction,
        mad_z_threshold=arguments.mad_z_threshold,
        core_max_normalized_radial_position=arguments.core_max_normalized_radial_position,
        peripheral_min_normalized_radial_position=arguments.peripheral_min_normalized_radial_position,
    )
    metadata = {key: getattr(arguments, key) for key in
                ("sample_id", "well_id", "field_id", "batch_id", "condition", "treatment", "dose", "timepoint")}
    result = analyze_multilevel_3d(organoids, cells, nuclei, spacing, config=config, metadata=metadata,
                                   nucleus_intensity=intensity)
    # Provenance capture belongs at this orchestration layer, not inside
    # analyze_multilevel_3d (deliberately headless/pure per multilevel3d's
    # module docstring). Mirrors the classical analyze() pipeline's provenance.
    _package_root = Path(__file__).resolve().parent.parent
    multilevel3d_dir = _package_root / "quantification" / "multilevel_relationships"
    multilevel3d_files = sorted(multilevel3d_dir.glob("*.py")) + [
        _package_root / "workflows" / "multilevel_measurement_workflow.py",
        _package_root / "result_export" / "measurement_tables.py",
    ]
    result.summary["provenance"] = {
        "pipeline_version": __version__,
        "git_commit": git_commit_hash(),
        "python": sys.version,
        "platform": platform.platform(),
        "source_code_sha256": {p.name: sha256(p) for p in multilevel3d_files},
        "packages": {package: importlib.metadata.version(package) for package in
                    ["numpy", "scipy", "scikit-image", "tifffile", "pandas"]},
    }
    output_paths = export_results(arguments.out, organoids=result.organoid_features, cells=result.cell_features,
                                  nuclei=result.nucleus_features, edges=result.cell_topology_edges,
                                  qc_flags=result.qc_flags, summary=result.summary,
                                  label_masks={"organoid": organoids, "cell": cells, "nucleus": nuclei},
                                  spacing_zyx_um=spacing)
    return {**result.summary, "output": str(Path(arguments.out).resolve()), "output_files": output_paths}
