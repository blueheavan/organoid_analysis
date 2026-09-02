"""One-field-at-a-time orchestration with provenance and explicit failure status."""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import importlib.metadata
import json
import platform
import sys
import time
import numpy as np
import pandas as pd
import yaml
from . import __version__
from .config import load_config
from .io import META_FIELDS, PATH_FIELDS, read_manifest, load_sample, load_truth_labels, write_labels, sha256
from .segmentation import QC_COLUMNS, compare_instance_qc, instance_qc_summary, segment
from .evaluation import match_instances
from .features import GEOMETRY_COLUMNS, MARKER_COLUMNS, measure_instances
from .viability import calibrate, classify
from .summary import make_summaries
from .stats import condition_pairwise_tests
from .report import plot_qc, plot_size_comparison, plot_morphology_viability, write_html


def _is_substantive(path: Path) -> bool:
    """True for entries that represent real content (exclude macOS/OS junk)."""
    if path.is_dir():
        return True
    return path.name not in {".DS_Store", "Thumbs.db", ".gitkeep"}


def _segment_with_qc(sample, cfg: dict) -> tuple:
    """Run the configured segmentation and, optionally, an independent QC reference.

    A reference watershed pass is useful when the primary mask comes from an
    imported model or probability map. It is never a fallback that silently
    changes measurements: disagreement only adds review flags and a row to the
    exported ``segmentation_qc.csv`` table.
    """
    primary = segment(sample.structure, sample.spacing, cfg["segmentation"],
                      sample.probability, sample.imported_labels)
    primary_summary = instance_qc_summary(primary.labels, sample.spacing)
    qc = {
        "reference_method": "none",
        "primary_instances": primary_summary["instances"],
        "reference_instances": np.nan,
        "primary_median_z_extent_um": primary_summary["median_z_extent_um"],
        "reference_median_z_extent_um": np.nan,
        "object_count_difference_fraction": np.nan,
        "z_extent_ratio": np.nan,
        "requires_review": False,
        "qc_flags": "",
    }
    if cfg["segmentation"]["qc_reference_method"] == "watershed":
        # Use identical physical parameters but an independent intensity-driven
        # route. Clearing its own QC option prevents recursive comparisons.
        reference_cfg = dict(cfg["segmentation"])
        reference_cfg.update(method="watershed", qc_reference_method="none")
        reference = segment(sample.structure, sample.spacing, reference_cfg)
        reference_summary = instance_qc_summary(reference.labels, sample.spacing)
        qc = {
            "reference_method": "watershed",
            **compare_instance_qc(
                primary.labels, reference.labels, sample.spacing,
                cfg["segmentation"]["qc_count_difference_threshold"],
                cfg["segmentation"]["qc_min_z_extent_ratio"],
                primary_summary=primary_summary,
                reference_summary=reference_summary,
            ),
        }
        if qc["requires_review"]:
            primary.flags.append("segmentation_method_disagreement_review")
    return primary, qc


def analyze(manifest: str | Path, out: str | Path, config: str | Path | None = None,
            keep_going: bool = False, synthetic: bool = False) -> dict:
    start = time.perf_counter()
    cfg = load_config(config)
    manifest_path = Path(manifest).resolve()
    design = read_manifest(manifest_path)
    out = Path(out).resolve()
    if out.exists() and any(_is_substantive(p) for p in out.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {out}. Choose a new output path; existing results are never overwritten.")
    out.mkdir(parents=True, exist_ok=True)
    for subdir in ["labels", "qc"]:
        (out / subdir).mkdir()
    sentinel = out / "RUN_INCOMPLETE.txt"
    sentinel.write_text("This run has not completed. Do not interpret partial results.\n", encoding="utf-8")
    (out / "config.resolved.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
    design.to_csv(out / "manifest.resolved.csv", index=False)
    origin = "synthetic_phantom" if synthetic else "user_images_unvalidated"
    provenance = {"pipeline_version": __version__, "started_utc": datetime.now(timezone.utc).isoformat(),
                  "data_origin": origin, "python": sys.version, "platform": platform.platform(),
                  "manifest": str(manifest_path), "manifest_sha256": sha256(manifest_path),
                  "configuration": cfg, "status": "incomplete",
                  "source_code_sha256": {p.name: sha256(p) for p in sorted(Path(__file__).parent.glob("*.py"))},
                  "packages": {package: importlib.metadata.version(package) for package in ["numpy", "scipy", "scikit-image", "tifffile", "pandas", "matplotlib", "PyYAML"]},
                  "input_files": [], "sample_metadata": {}}
    paths = sorted({row[field] for row in design.to_dict("records") for field in PATH_FIELDS if row[field]})
    for path in paths:
        provenance["input_files"].append({"path": path, "sha256": sha256(path), "size_bytes": Path(path).stat().st_size})
    rows, sample_rows, failures, qc_rows = [], [], [], []
    validation_matches, validation_metrics = [], []
    for position, row in enumerate(design.to_dict("records"), start=1):
        meta = {key: row[key] for key in META_FIELDS}
        meta["data_origin"] = origin
        try:
            sample = load_sample(row, cfg)
            segmentation, segmentation_qc = _segment_with_qc(sample, cfg)
            mesh_dir = out / "meshes" / row["sample_id"] if cfg["report"]["save_meshes"] else None
            features, preview_meshes = measure_instances(sample, segmentation, cfg, mesh_dir)
            write_labels(out / "labels" / f'{row["sample_id"]}.labels.ome.tif', segmentation.labels, sample.spacing)
            truth = load_truth_labels(row, segmentation.labels.shape, sample.spacing)
            if truth is not None:
                matches, metrics = match_instances(truth, segmentation.labels)
                validation_matches.extend({**meta, **match} for match in matches.to_dict("records"))
                validation_metrics.append({**meta, **metrics})
            plot_qc(sample, segmentation, preview_meshes, row["sample_id"], out / "qc" / f'{row["sample_id"]}.png', synthetic)
            rows.extend({**meta, **feature} for feature in features)
            qc_rows.append({**meta, **segmentation_qc})
            sample_rows.append({**meta, "status": "ok", "sample_flags": ";".join(segmentation.flags),
                                "foreground_fraction": segmentation.foreground_fraction,
                                "segmentation_threshold": segmentation.threshold,
                                "removed_small_instances": segmentation.removed_small_instances,
                                "segmentation_qc_reference": segmentation_qc["reference_method"],
                                "segmentation_qc_flags": segmentation_qc["qc_flags"],
                                "spacing_z_um": sample.spacing[0], "spacing_y_um": sample.spacing[1], "spacing_x_um": sample.spacing[2]})
            provenance["sample_metadata"][row["sample_id"]] = sample.metadata
            print(f'[{position}/{len(design)}] {row["sample_id"]}: {len(features)} objects; {sum(f["morphology_eligible"] for f in features)} pass geometry QC',flush=True)
        except Exception as error:
            failure = {"sample_id": row["sample_id"], "error_type": type(error).__name__, "message": str(error)}
            failures.append(failure)
            (out / "errors.json").write_text(json.dumps(failures, indent=2), encoding="utf-8")
            sample_rows.append({**meta, "status": "failed", "sample_flags": str(error),
                                 "foreground_fraction": np.nan, "segmentation_threshold": np.nan,
                                 "removed_small_instances": np.nan,
                                 "segmentation_qc_reference": "", "segmentation_qc_flags": "",
                                 "spacing_z_um": np.nan, "spacing_y_um": np.nan, "spacing_x_um": np.nan})
            if not keep_going:
                provenance["failures"] = failures
                (out / "provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
                raise
            print(f'[{position}/{len(design)}] FAILED {row["sample_id"]}: {error}',file=sys.stderr,flush=True)
    columns = META_FIELDS + ["data_origin"] + GEOMETRY_COLUMNS + MARKER_COLUMNS
    objects = pd.DataFrame(rows, columns=columns)
    # Dtypes matter when all fields are empty: keep a consistent machine-readable schema.
    for column in ["morphology_eligible", "viability_measurement_eligible", "touches_border"]:
        objects[column] = objects[column].astype(bool)
    failed_units = {(r["batch_id"],r["unit_id"]) for r in sample_rows if r["status"] != "ok"}
    calibration_input = objects[np.array([(r.batch_id,r.unit_id) not in failed_units for r in objects.itertuples()],dtype=bool)]
    calibration = calibrate(calibration_input, list(dict.fromkeys(design.batch_id)), cfg["viability"])
    objects = classify(objects, calibration, cfg["viability"])
    sample_table, units, replicates, condition_table = make_summaries(objects, pd.DataFrame(sample_rows), cfg["report"])
    for table in [units, replicates, condition_table, calibration]:
        table["data_origin"] = origin
    tables = {"organoids": objects, "sample_summary": sample_table, "unit_summary": units,
              "replicate_summary": replicates, "condition_summary": condition_table, "calibration": calibration,
              "segmentation_qc": pd.DataFrame(qc_rows, columns=META_FIELDS + QC_COLUMNS)}
    if validation_metrics:
        tables["segmentation_validation_metrics"] = pd.DataFrame(validation_metrics)
        tables["segmentation_validation_matches"] = pd.DataFrame(validation_matches)
    stats_omnibus = {}
    if cfg["stats"]["enabled"]:
        pairwise_contrasts, stats_omnibus = condition_pairwise_tests(
            objects, tuple(cfg["stats"]["features"]), cfg["stats"]["min_replicates_per_condition"])
        tables["pairwise_contrasts"] = pairwise_contrasts
    for name, table in tables.items():
        table.to_csv(out / f"{name}.csv", index=False, float_format="%.10g")
    if cfg["stats"]["enabled"]:
        (out / "stats_results.json").write_text(json.dumps(stats_omnibus, indent=2), encoding="utf-8")
    plot_size_comparison(objects, units, replicates, condition_table, out, synthetic)
    plot_morphology_viability(objects, units, replicates, condition_table, out, synthetic)
    write_html(out, objects, sample_table, condition_table, calibration, synthetic, failures, cfg)
    provenance.update(status="partial" if failures else "complete", failures=failures,
                      completed_utc=datetime.now(timezone.utc).isoformat(), elapsed_seconds=time.perf_counter()-start,
                       n_detected_organoids=len(objects), n_geometry_eligible_organoids=int(objects.morphology_eligible.sum()),
                       n_samples=len(design), n_failed_samples=len(failures),
                       n_samples_with_segmentation_truth=len(validation_metrics))
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    if failures:
        sentinel.rename(out / "PARTIAL_RUN.txt")
        (out / "PARTIAL_RUN.txt").write_text("Some fields failed. See errors.json and the report. Incomplete wells are excluded from replicate summaries.\n",encoding="utf-8")
    else:
        sentinel.unlink()
    return {"status": provenance["status"], "out": str(out), "n_samples": len(design),
            "n_detected_organoids": len(objects), "n_geometry_eligible_organoids": int(objects.morphology_eligible.sum()),
            "n_failed_samples": len(failures), "elapsed_seconds": provenance["elapsed_seconds"]}
