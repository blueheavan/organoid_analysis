"""Seeded synthetic phantoms. These are NOT real organoid microscopy observations."""
from __future__ import annotations

from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import ndimage as ndi
import tifffile
import yaml
from .config import load_config
from .pipeline import analyze
from .evaluation import match_instances


def generate_demo(out: str | Path, seed: int = 1729) -> Path:
    out = Path(out).resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"Demo directory is not empty: {out}; select a new directory")
    (out / "images").mkdir(parents=True, exist_ok=True)
    (out / "truth").mkdir()
    rng = np.random.default_rng(seed)
    spacing = (3.0, 1.0, 1.0)
    shape = (32, 128, 144)
    zz, yy, xx = np.indices(shape, dtype=float)
    zz, yy, xx = zz * spacing[0], yy * spacing[1], xx * spacing[2]
    design, truth_rows = [], []
    conditions = [("Vehicle", "sample", 1.0, [0.80, 0.15, 0.05]),
                  ("Treatment_A", "sample", .86, [0.20, 0.60, 0.20]),
                  ("Treatment_B", "sample", .72, [0.05, 0.15, 0.80]),
                  ("Live_control", "live", 1.0, [1., 0., 0.]),
                  ("Dead_control", "dead", 1.0, [0., 0., 1.])]
    states = ["viable_like", "mixed_signal", "compromised_like"]
    # Three independent synthetic biological replicates, paired across conditions.
    for replicate in range(1,4):
        replicate_size = [0.96, 1.05, 1.0][replicate-1]
        for condition, role, size_factor, state_probabilities in conditions:
            sample_id = f"SYN_{condition}_R{replicate}"
            truth = np.zeros(shape, np.uint16)
            channels = np.zeros((3, *shape), np.float32)
            for object_id, (cy, cx) in enumerate([(37,25),(37,72),(37,119),(91,25),(91,72),(91,119)], start=1):
                radius = rng.uniform(13,18) * size_factor * replicate_size
                radii = radius * rng.uniform(.86,1.14,size=3)
                cz = 46 + rng.uniform(-3,3)
                norm = ((zz-cz)/radii[0])**2 + ((yy-cy)/radii[1])**2 + ((xx-cx)/radii[2])**2
                mask = norm <= 1
                if (truth[mask] != 0).any():
                    raise RuntimeError("Synthetic objects unexpectedly overlap")
                truth[mask] = object_id
                state = rng.choice(states,p=state_probabilities)
                c_level,p_level = {"viable_like":(1250,50), "mixed_signal":(650,570), "compromised_like":(70,1150)}[state]
                scale = rng.uniform(.93,1.07)
                texture = 1 + .07*np.sin(xx/2.8)*np.cos(yy/3.4)
                channels[0][mask] = (1800*texture)[mask]
                channels[1][mask] = (c_level*scale*texture)[mask]
                channels[2][mask] = (p_level*scale*texture)[mask]
                truth_rows.append({"data_origin":"synthetic_phantom","sample_id":sample_id,"truth_id":object_id,
                                   "true_state":state,"voxelized_volume_um3":int(mask.sum())*np.prod(spacing),
                                   "analytic_ellipsoid_volume_um3":4*np.pi*np.prod(radii)/3,
                                   "center_z_um":cz,"center_y_um":cy,"center_x_um":cx})
            for channel in range(3):
                channels[channel] = ndi.gaussian_filter(channels[channel],sigma=(.35,.6,.6))
                channels[channel] += 90 + rng.normal(0,9,size=shape)
            channels = np.clip(np.rint(channels),0,65535).astype(np.uint16)
            metadata = {"axes":"CZYX","PhysicalSizeZ":spacing[0],"PhysicalSizeY":spacing[1],"PhysicalSizeX":spacing[2],
                        "PhysicalSizeZUnit":"µm","PhysicalSizeYUnit":"µm","PhysicalSizeXUnit":"µm",
                        "Description":"SYNTHETIC PHANTOM, NOT EXPERIMENTAL DATA",
                        "Channel":{"Name":["Synthetic structure","Synthetic Calcein","Synthetic PI"]}}
            tifffile.imwrite(out / "images" / f"{sample_id}.ome.tif", channels,ome=True,photometric="minisblack",compression="zlib",metadata=metadata)
            tifffile.imwrite(out / "truth" / f"{sample_id}.truth.ome.tif",truth,ome=True,photometric="minisblack",compression="zlib",
                             metadata={key:value for key,value in {**metadata,"axes":"ZYX"}.items() if key!="Channel"})
            design.append({"sample_id":sample_id,"condition":condition,"biological_replicate":f"Synthetic_R{replicate}",
                           "unit_id":sample_id,"batch_id":"synthetic_batch_1","control":role,
                           "image_path":f"images/{sample_id}.ome.tif","axes":"CZYX",
                           "spacing_z_um":spacing[0],"spacing_y_um":spacing[1],"spacing_x_um":spacing[2]})
    manifest = out / "manifest.csv"
    pd.DataFrame(design).to_csv(manifest,index=False)
    pd.DataFrame(truth_rows).to_csv(out / "ground_truth.csv",index=False)
    cfg = load_config()
    cfg["channels"] = {"structure":0,"calcein":1,"pi":2}
    cfg["segmentation"].update(gaussian_sigma_um=.6, closing_radius_um=0., min_volume_um3=1000., seed_min_distance_um=12.)
    cfg["viability"]["mode"] = "controls"
    (out / "config.yaml").write_text(yaml.safe_dump(cfg,sort_keys=False),encoding="utf-8")
    (out / "SYNTHETIC_DATA_NOTICE.txt").write_text("All images, treatment assignments, marker states, and differences in this folder are deliberately simulated. They are software test inputs, not biological evidence.\n",encoding="utf-8")
    return manifest


def validate_demo(out: Path) -> dict:
    measured = pd.read_csv(out / "results" / "organoids.csv")
    truth = pd.read_csv(out / "ground_truth.csv")
    comparisons, detection_rows = [], []
    for sample_id, expected in truth.groupby("sample_id"):
        predicted = measured[measured.sample_id == sample_id]
        labels = tifffile.imread(out / "results" / "labels" / f"{sample_id}.labels.ome.tif")
        truth_labels = tifffile.imread(out / "truth" / f"{sample_id}.truth.ome.tif")
        matching, detection = match_instances(truth_labels,labels,iou_threshold=.5)
        detection_rows.append({"sample_id":sample_id,**detection})
        matches = matching[matching.match_status=="TP"].set_index("true_id")
        # Bounding boxes once per image instead of a full-volume `labels == id`
        # scan per matched object.
        pred_boxes = ndi.find_objects(labels)
        truth_boxes = ndi.find_objects(truth_labels)
        for row in expected.itertuples():
            if row.truth_id not in matches.index:
                comparisons.append({"sample_id":sample_id,"truth_id":row.truth_id,"matched":False})
                continue
            match = int(matches.loc[row.truth_id,"predicted_id"])
            candidate = predicted[predicted.organoid_id == match].iloc[0]
            crop = tuple(
                slice(min(p.start, t.start), max(p.stop, t.stop))
                for p, t in zip(pred_boxes[match - 1], truth_boxes[row.truth_id - 1])
            )
            a,b = labels[crop] == match, truth_labels[crop] == row.truth_id
            dice = float(2*np.count_nonzero(a&b)/(a.sum()+b.sum()))
            comparisons.append({"sample_id":sample_id,"truth_id":row.truth_id,"matched":True,
                                "dice":dice,"relative_volume_error":abs(candidate.volume_um3-row.voxelized_volume_um3)/row.voxelized_volume_um3,
                                "true_state":row.true_state,"predicted_state":candidate.viability_state,
                                "state_agrees":candidate.viability_state==row.true_state})
    validation = pd.DataFrame(comparisons)
    validation.to_csv(out / "synthetic_validation_per_object.csv",index=False)
    detections = pd.DataFrame(detection_rows)
    detections.to_csv(out / "synthetic_detection_per_image.csv",index=False)
    # Disaggregate state agreement by whether the sample's own controls fed
    # its batch's calibration endpoints (live/dead) or were genuinely held
    # out (treatment conditions) -- the pooled figure alone mixes the two and
    # can look better than the held-out (treatment) agreement really is.
    sample_control = measured[["sample_id", "control"]].drop_duplicates().set_index("sample_id")["control"]
    matched = validation[validation.matched.astype(bool)]
    is_control_sample = matched.sample_id.map(sample_control).isin(["live", "dead"])
    treatment_agreement = matched.loc[~is_control_sample, "state_agrees"]
    control_agreement = matched.loc[is_control_sample, "state_agrees"]
    summary = {"data_origin":"synthetic_phantom", "limitation":"Software verification on simple phantoms; not real-image accuracy or independent classifier validation. Control scaling and simulation intentionally share simple marker assumptions.",
               "expected_objects":len(truth),"detected_objects":len(measured),"matched_objects":int(validation.matched.sum()),
               "instance_matching":"Hungarian assignment at IoU >= 0.5",
               "false_positives":int(detections.false_positives.sum()),"false_negatives":int(detections.false_negatives.sum()),
               "mean_image_detection_precision":float(detections.precision.mean()),"mean_image_detection_recall":float(detections.recall.mean()),
               "median_dice":float(validation.dice.median()),"minimum_dice":float(validation.dice.min()),
               "median_absolute_relative_volume_error":float(validation.relative_volume_error.median()),
               "maximum_absolute_relative_volume_error":float(validation.relative_volume_error.max()),
               "state_agreement_fraction_matched_objects":float(validation.state_agrees.mean()),
               "state_agreement_fraction_treatment_conditions":float(treatment_agreement.mean()) if len(treatment_agreement) else float("nan"),
               "state_agreement_fraction_live_dead_control_conditions":float(control_agreement.mean()) if len(control_agreement) else float("nan")}
    (out / "synthetic_validation.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    return summary


def run_demo(out: str | Path, seed: int = 1729) -> dict:
    out = Path(out).resolve()
    manifest = generate_demo(out,seed)
    result = analyze(manifest,out/"results",out/"config.yaml",synthetic=True)
    result["synthetic_validation"] = validate_demo(out)
    return result
