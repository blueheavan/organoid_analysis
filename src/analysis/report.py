"""Publication-exportable plots, orthogonal QC, 3D surfaces, and an offline HTML report."""
from __future__ import annotations

from pathlib import Path
import base64
import html
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy import ndimage as ndi
from skimage.segmentation import find_boundaries
from .viability import STATES, STATE_COLORS

COLORS = ["#277DA8", "#D77632", "#665EA8", "#329E82", "#BA5578", "#9C823E", "#607C8E", "#6D9361"]


def _style() -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 12,
                         "axes.labelsize": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "svg.fonttype": "none", "pdf.fonttype": 42, "savefig.facecolor": "white"})


def _export(fig, stem: Path) -> None:
    for suffix in ["png", "svg", "pdf"]:
        fig.savefig(stem.with_suffix(f".{suffix}"), dpi=220, bbox_inches="tight")
    plt.close(fig)


def _norm(array: np.ndarray | None, shape=None) -> np.ndarray:
    if array is None:
        return np.zeros(shape, dtype=np.float32)
    low, high = np.percentile(array, [1, 99.8])
    if high <= low:
        return np.zeros(array.shape, dtype=np.float32)
    return np.clip((array.astype(np.float32) - low) / (high - low), 0, 1)


def _overlay(intensity, labels):
    base = np.repeat(_norm(intensity)[..., None], 3, axis=-1)
    base[find_boundaries(labels, mode="inner")] = [1., .79, .20]
    return base


def plot_qc(sample, segmentation, meshes: list, sample_id: str, path: Path, synthetic=False) -> None:
    _style()
    labels = segmentation.labels
    if labels.any():
        counts = np.bincount(labels.ravel())
        counts[0] = 0
        largest_id = counts.argmax()
        z, y, x = np.mean(np.argwhere(labels == largest_id), axis=0).astype(int)
    else:
        z, y, x = np.array(labels.shape) // 2
    sz, sy, sx = sample.spacing
    nz, ny, nx = labels.shape
    fig = plt.figure(figsize=(12, 7.8), layout="constrained")
    grid = fig.add_gridspec(2, 3)
    views = [(sample.structure[z], labels[z], [0, nx*sx, ny*sy, 0], f"XY slice · z = {z*sz:g} µm", "X (µm)", "Y (µm)"),
             (sample.structure[:, y, :], labels[:, y, :], [0, nx*sx, nz*sz, 0], f"XZ slice · y = {y*sy:g} µm", "X (µm)", "Z (µm)"),
             (sample.structure[:, :, x], labels[:, :, x], [0, ny*sy, nz*sz, 0], f"YZ slice · x = {x*sx:g} µm", "Y (µm)", "Z (µm)")]
    for col, (intensity, label, extent, title, xlabel, ylabel) in enumerate(views):
        ax = fig.add_subplot(grid[0, col])
        ax.imshow(_overlay(intensity, label), extent=extent, interpolation="nearest")
        ax.set(title=title, xlabel=xlabel, ylabel=ylabel)
    ax = fig.add_subplot(grid[1, 0])
    ax.imshow(_norm(sample.structure.max(axis=0)), cmap="gray", extent=[0,nx*sx,ny*sy,0])
    # Project only each object's bounding box. The previous full-volume
    # ``labels == id`` scan made report generation scale as objects × voxels.
    label_values = np.unique(labels)
    label_values = label_values[label_values != 0]
    if len(label_values) and int(label_values[-1]) > max(100_000, 4 * len(label_values)):
        local_labels = (np.searchsorted(label_values, labels) + 1).astype(np.uint32)
        local_labels[labels == 0] = 0
        entries = enumerate(label_values, start=1)
    else:
        local_labels = labels
        entries = ((int(value), value) for value in label_values)
    boxes = ndi.find_objects(local_labels)
    for local_id, object_id in entries:
        bbox = boxes[int(local_id) - 1]
        if bbox is None:
            continue
        footprint = (local_labels[bbox] == local_id).any(axis=0)
        if footprint.any() and not footprint.all():
            y0, y1 = bbox[1].start, bbox[1].stop
            x0, x1 = bbox[2].start, bbox[2].stop
            ax.contour(footprint, levels=[.5], colors=[COLORS[(int(object_id)-1) % len(COLORS)]],
                       linewidths=.9, extent=[x0*sx,x1*sx,y1*sy,y0*sy], origin="upper")
            yy, xx = np.argwhere(footprint).mean(axis=0)
            ax.text((xx+x0)*sx, (yy+y0)*sy, str(int(object_id)), fontsize=9, color="white", ha="center",
                    bbox={"facecolor":"#243242", "alpha":.7,"edgecolor":"none","pad":1.2})
    ax.set(title="Structure MIP + projected instance outlines", xlabel="X (µm)", ylabel="Y (µm)")
    ax = fig.add_subplot(grid[1, 1])
    merge = np.zeros((ny, nx, 3), float)
    if sample.calcein is not None:
        merge[:, :, 1] = _norm(sample.calcein.max(axis=0))
    if sample.pi is not None:
        merge[:, :, 0] = _norm(sample.pi.max(axis=0))
    ax.imshow(merge, extent=[0,nx*sx,ny*sy,0])
    ax.set(title="Calcein (green) / PI (red) MIP", xlabel="X (µm)", ylabel="Y (µm)")
    ax.text(.5, -.22, "Display contrast only; raw intensities quantify markers.", transform=ax.transAxes,
            ha="center", fontsize=8, color="#596577")
    ax = fig.add_subplot(grid[1, 2], projection="3d")
    for object_id, vertices, faces in meshes:
        xyz = vertices[:, ::-1]
        collection = Poly3DCollection(xyz[faces], alpha=.82, linewidths=0,
                                     facecolor=COLORS[(object_id-1) % len(COLORS)], edgecolor="none")
        ax.add_collection3d(collection)
    ax.set(xlim=(0,nx*sx), ylim=(0,ny*sy), zlim=(0,nz*sz), title="3D outer surfaces · categorical colors",
           xlabel="X (µm)", ylabel="Y (µm)", zlabel="Z (µm)")
    ax.set_box_aspect((nx*sx, ny*sy, nz*sz))
    ax.view_init(elev=25, azim=-60)
    ax.tick_params(labelsize=8)
    prefix = "SYNTHETIC DEMO · " if synthetic else ""
    fig.suptitle(f"{prefix}{sample_id} · segmentation quality control", fontsize=15, weight="bold")
    fig.savefig(path, dpi=155, bbox_inches="tight")
    plt.close(fig)


def _plot_inputs(objects, units, replicates):
    complete = set(map(tuple, units.loc[units.unit_complete, ["batch_id", "unit_id"]].to_numpy()))
    mask = np.array([(r.batch_id, r.unit_id) in complete for r in objects.itertuples()], bool)
    included = objects[objects.morphology_eligible.astype(bool) & (objects.control == "sample") & mask]
    reps = replicates[replicates.control == "sample"]
    conditions = list(dict.fromkeys(reps.condition))
    palette = {condition: COLORS[i % len(COLORS)] for i, condition in enumerate(conditions)}
    return included, reps, conditions, palette


def plot_size_comparison(objects, units, replicates, conditions_table, out: Path, synthetic=False) -> None:
    _style()
    included, reps, conditions, palette = _plot_inputs(objects, units, replicates)
    fig, axes = plt.subplots(1, 2, figsize=(max(11, len(conditions)*1.8), 4.9), layout="constrained")
    a, b = axes
    for condition in conditions:
        values = np.sort(included.loc[included.condition == condition, "volume_um3"].to_numpy(float))
        if values.size:
            a.step(values, np.arange(1,len(values)+1)/len(values), where="post", color=palette[condition],
                   lw=2.2, label=f"{condition} (n = {len(values)} organoids)")
    if len(included):
        a.set_xscale("log")
        a.legend(loc="upper left", frameon=False, fontsize=9)
    else:
        a.text(.5,.5,"No QC-eligible treatment organoids",ha="center",transform=a.transAxes)
    a.set(xlabel="Organoid outer-envelope volume (µm³; log scale)", ylabel="Cumulative fraction of organoids",
          ylim=(0,1.04), title="A  Size distributions · descriptive, pooled objects")
    rng = np.random.default_rng(0)
    has_replicate_values = False
    for index, condition in enumerate(conditions):
        selected = reps[reps.condition == condition]
        values = selected.median_of_unit_medians_volume_um3.dropna().to_numpy(float)
        if len(values):
            has_replicate_values = True
            jitter = rng.uniform(-.11,.11,len(values))
            b.scatter(index+jitter, values, s=65, color=palette[condition], edgecolor="white", linewidth=.8, zorder=3)
            center = values.mean()
            b.plot([index-.18,index+.18],[center,center],color="#202D3C",lw=2,zorder=4)
            summary = conditions_table[(conditions_table.condition == condition) & (conditions_table.control == "sample")].iloc[0]
            low, high = summary.ci95_low_volume_um3, summary.ci95_high_volume_um3
            if np.isfinite(low) and np.isfinite(high):
                b.vlines(index,low,high,color="#202D3C",lw=1.2,zorder=2)
                b.hlines([low,high],index-.07,index+.07,color="#202D3C",lw=1.2,zorder=2)
    if has_replicate_values:
        b.set_yscale("log")
    ticklabels = [f"{textwrap.fill(str(c), 18)}\nn = {int(reps.loc[reps.condition==c, 'median_of_unit_medians_volume_um3'].notna().sum())} replicates" for c in conditions]
    b.set_xticks(range(len(conditions)), ticklabels)
    b.set(ylabel="Replicate median of well medians (µm³; log scale)", title="B  Biological replicates · wells weighted equally")
    for axis in axes:
        axis.grid(axis="y", color="#E6EAF0",lw=.7)
        axis.set_axisbelow(True)
    label = "SYNTHETIC DEMO — not experimental results" if synthetic else "Treatment comparison — exploratory analysis"
    fig.suptitle(label, fontsize=15, weight="bold")
    fig.supxlabel("B: points = independent biological replicates; bar = mean; whiskers = replicate bootstrap 95% CI (n ≥ 3).",fontsize=9,color="#546171")
    _export(fig, out / "size_comparison")


def plot_morphology_viability(objects, units, replicates, conditions_table, out: Path, synthetic=False) -> None:
    _style()
    included, reps, conditions, palette = _plot_inputs(objects, units, replicates)
    fig, axes = plt.subplots(1,2,figsize=(max(11,len(conditions)*1.8),4.9),layout="constrained")
    a,b = axes
    for condition in conditions:
        values = included[included.condition == condition]
        a.scatter(values.volume_um3, values.sphericity, s=29, alpha=.7, color=palette[condition], label=str(condition))
    if len(included):
        a.set_xscale("log")
        a.legend(frameon=False,fontsize=9)
    a.set(xlabel="Outer-envelope volume (µm³; log scale)", ylabel="Sphericity (dimensionless)",
          title="A  Size and shape · pooled objects")
    bottoms = np.zeros(len(conditions))
    for state in STATES:
        values = []
        for condition in conditions:
            series = conditions_table[(conditions_table.condition == condition) & (conditions_table.control == "sample")][f"mean_replicate_fraction_{state}"]
            values.append(float(series.iloc[0]) if len(series) and pd.notna(series.iloc[0]) else 0.)
        b.bar(range(len(conditions)), values, bottom=bottoms, color=STATE_COLORS[state], width=.65,
              label=state.replace("_"," "), edgecolor="white", linewidth=.6)
        bottoms += values
    b.set_xticks(range(len(conditions)), [textwrap.fill(str(c),18) for c in conditions])
    b.set(ylabel="Mean fraction of organoids per biological replicate", ylim=(0,1.07),
          title="B  Control-scaled viability signal states")
    b.legend(loc="upper center",bbox_to_anchor=(.5,-.10),ncol=2,frameon=False,fontsize=9)
    for axis in axes:
        axis.grid(axis="y",color="#E6EAF0",lw=.7)
        axis.set_axisbelow(True)
    fig.suptitle("SYNTHETIC DEMO — morphology and marker states" if synthetic else "Morphology and marker-defined viability states",fontsize=15,weight="bold")
    fig.supxlabel("Viability states are organoid-level marker patterns, not cell survival percentages or apoptosis diagnoses.",fontsize=9,color="#546171")
    _export(fig,out / "morphology_viability")


def _image_data(path: Path) -> str:
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode("ascii")


def write_html(out: Path, objects: pd.DataFrame, samples: pd.DataFrame, conditions: pd.DataFrame,
               calibration: pd.DataFrame, synthetic: bool, failures: list, cfg: dict) -> None:
    title = "SYNTHETIC DEMO: 3D organoid analysis" if synthetic else "3D organoid analysis report"
    banner = "Synthetic phantoms only. These figures do not describe experimental treatment effects." if synthetic else "Research output. Inspect segmentation and validate marker gates before biological interpretation."
    banner += " PARTIAL RUN: some images failed; incomplete wells are excluded from replicate summaries." if failures else ""
    treatment = conditions[conditions.control == "sample"].copy()
    columns = ["condition","n_biological_replicates_total","n_biological_replicates_with_size_data",
               "n_included_organoids_complete_units","mean_replicate_median_volume_um3",
               "ci95_low_volume_um3","ci95_high_volume_um3"]
    summary_html = treatment[columns].to_html(index=False,escape=True,na_rep="—",float_format=lambda x:f"{x:,.1f}")
    calibration_html = calibration[["batch_id","status","reason","live_control_replicates","dead_control_replicates"]].to_html(index=False,escape=True,na_rep="—")
    details = []
    for row in samples.to_dict("records"):
        sid = row["sample_id"]
        qc = out / "qc" / f"{sid}.png"
        if qc.is_file():
            details.append(f'<details><summary>{html.escape(sid)} · {html.escape(str(row["condition"]))} · {row["n_included"]} included / {row["n_detected"]} detected</summary><img alt="Orthogonal slices, marker MIP and physical 3D surfaces" src="{_image_data(qc)}"></details>')
    object_columns = ["sample_id","condition","organoid_id","volume_um3","surface_area_um2","sphericity","viability_state","morphology_flags"]
    object_html = objects[object_columns].head(100).to_html(index=False,escape=True,na_rep="—",float_format=lambda x:f"{x:,.3f}")
    failure_html = "" if not failures else pd.DataFrame(failures).to_html(index=False,escape=True)
    links = ["organoids.csv","sample_summary.csv","unit_summary.csv","replicate_summary.csv","condition_summary.csv","calibration.csv","provenance.json"]
    qc_path = out / "segmentation_qc.csv"
    qc_html = ""
    if qc_path.is_file():
        qc = pd.read_csv(qc_path)
        reviewed = qc[qc.reference_method.ne("none")] if "reference_method" in qc else qc.iloc[0:0]
        if len(reviewed):
            qc_columns = ["sample_id", "reference_method", "primary_instances", "reference_instances",
                          "primary_median_z_extent_um", "reference_median_z_extent_um",
                          "object_count_difference_fraction", "z_extent_ratio", "requires_review", "qc_flags"]
            qc_table_html = reviewed[qc_columns].to_html(
                index=False, escape=True, na_rep="—", float_format=lambda x: f"{x:,.3f}"
            )
            qc_html = (
                "<section><h2>Independent segmentation-method QC</h2>"
                "<p>The primary mask remains the measurement mask. The watershed reference is an independent "
                "screen for large object-count or Z-extent disagreement, not ground truth and not an automatic "
                "method selector. Rows marked for review require inspection of the saved orthogonal QC panels.</p>"
                f'<div class="scroll">{qc_table_html}</div></section>'
            )
            links.append("segmentation_qc.csv")
    metrics_path = out / "segmentation_validation_metrics.csv"
    validation_html = ""
    if metrics_path.is_file():
        metrics = pd.read_csv(metrics_path)
        validation_columns = ["sample_id", "n_true", "n_predicted", "true_positives", "false_negatives",
                              "false_positives", "precision", "recall", "mean_dice_matched",
                              "mean_iou_matched", "panoptic_quality", "possible_splits", "possible_merges"]
        validation_table_html = metrics[validation_columns].to_html(
            index=False, escape=True, na_rep="—", float_format=lambda x: f"{x:,.3f}"
        )
        validation_html = (
            "<section><h2>Annotated segmentation validation</h2>"
            "<p>These metrics compare the exported mask with an independently registered instance annotation using "
            "one-to-one Hungarian matching at the recorded IoU threshold. Dice and IoU are reported only for matched "
            "objects; precision, recall, split/merge flags and panoptic quality expose detection errors that foreground "
            "overlap alone can hide. Validation is available only for fields with <code>truth_labels_path</code>.</p>"
            f'<div class="scroll">{validation_table_html}</div></section>'
        )
        links.extend(["segmentation_validation_metrics.csv", "segmentation_validation_matches.csv"])
    downloads = " · ".join(f'<a href="{name}">{name}</a>' for name in links)
    css = "body{font:15px/1.55 system-ui,sans-serif;color:#233143;max-width:1200px;margin:36px auto;padding:0 22px;background:#f6f8fb}h1,h2{line-height:1.2}h1{font-size:32px}h2{font-size:23px;margin-top:30px}p{max-width:1100px}.notice{padding:16px 20px;background:#fff1d6;border-left:5px solid #cb8a2a;border-radius:5px}section{background:white;padding:22px;border-radius:9px;margin:20px 0;border:1px solid #e0e5ee}img{width:100%;height:auto}table{border-collapse:collapse;white-space:nowrap;font-size:12px}th,td{padding:9px 12px;text-align:right;border-bottom:1px solid #e2e7ee}th{background:#eaf0f7;color:#233143}th:first-child,td:first-child{text-align:left}.scroll{overflow-x:auto}summary{cursor:pointer;font-weight:600;padding:12px;border-bottom:1px solid #e0e5ee}a{color:#1b6d99}.small{font-size:13px;color:#5c6674}code{background:#edf0f4;padding:2px 5px;border-radius:3px}"
    document = f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{html.escape(title)}</title><style>{css}</style></head><body>
<h1>{html.escape(title)}</h1><p class="notice">{html.escape(banner)}</p>
<p>Measurements use calibrated Z/Y/X spacing. Volume and surface area describe the outer envelope, including fully enclosed lumens. Border-truncated objects and other geometry QC failures remain in the object table but are excluded from size summaries.</p>
<section><h2>Treatment size comparisons</h2><img src="{_image_data(out/'size_comparison.png')}" alt="Volume ECDFs and biological replicate summaries"><p class="small">The pooled distribution describes measured objects and can overweight organoid-rich wells. The replicate plot first pools fields within wells, then uses the median of well medians within each biological replicate. Confidence intervals resample biological replicates only; no object-level significance tests are used. Bootstrap intervals with few replicates are unstable. Empty wells remain in count summaries but have no size median.</p><div class="scroll">{summary_html}</div></section>
<section><h2>Morphology and viability</h2><img src="{_image_data(out/'morphology_viability.png')}" alt="Size versus sphericity and organoid signal states"><p>Calcein and PI means are measured from raw data after local background subtraction. Each marker is scaled separately to matched live/dead controls within its acquisition batch. Gates ({cfg['viability']['high_gate']:g} high; {cfg['viability']['low_gate']:g} low) are configurable starting rules, not universal biological thresholds. Both-low signals, missing controls, saturation, or other measurement failures remain indeterminate. Mixed signal does not establish apoptosis or a cell survival fraction.</p><div class="scroll">{calibration_html}</div></section>
 <section><h2>Segmentation quality control</h2><p class="small">Expand a sample to inspect true XY/XZ/YZ slices, MIPs, and reconstructed surfaces. Colors of surfaces identify instances; they are not heatmaps. Preview meshes are coarser than the meshes used for measurements. Intensity contrast is adjusted for display only. Review border-truncated objects, implausible joins/splits, weak Z continuity, and method-disagreement flags before interpreting morphology or marker states.</p>{''.join(details)}{failure_html}</section>
 {qc_html}{validation_html}
<section><h2>Structured results</h2><p>{downloads}</p><p class="small">First 100 object rows below; the CSV contains every object, QC flags, marker readouts, and replicate identifiers. Keep this report beside its tables when using download links. All figures and sample QC images are embedded and work offline.</p><div class="scroll">{object_html}</div></section>
<p class="small">Methods: <a href="https://scikit-image.org/docs/stable/api/skimage.measure.html#skimage.measure.marching_cubes">spacing-aware marching cubes</a>; <a href="https://scikit-image.org/docs/stable/auto_examples/segmentation/plot_watershed.html">3D watershed</a>; <a href="https://www.cgohlke.com/docs/tifffile/">TIFF metadata</a>. See README.md and docs/METHODS.md in the package for assumptions and limitations.</p></body></html>'''
    (out / "report.html").write_text(document,encoding="utf-8")
