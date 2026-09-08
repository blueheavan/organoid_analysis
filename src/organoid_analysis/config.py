"""Explicit defaults, strict configuration validation, and units."""
from __future__ import annotations

import math
from copy import deepcopy
from pathlib import Path

import yaml

DEFAULTS = {
    "channels": {"structure": 0, "calcein": None, "pi": None},
    "segmentation": {
        "method": "watershed",  # watershed | probability | labels
        "polarity": "bright",  # dark is an exploratory brightfield baseline
        "gaussian_sigma_um": 1.0,
        "background_sigma_um": 0.0,
        "threshold": "otsu",  # number is in corrected raw intensity units
        "probability_threshold": 0.5,
        "closing_radius_um": 1.0,
        "fill_enclosed_holes": True,
        "min_volume_um3": 2000.0,
        "seed_h_um": 2.0,
        "seed_min_distance_um": 15.0,
        "split_touching": True,
        "max_foreground_fraction": 0.70,
        # A second, independent watershed pass can flag disagreement with an
        # imported-label or probability-map primary segmentation. It never
        # replaces the primary mask automatically.
        "qc_reference_method": "none",  # none | watershed
        "qc_count_difference_threshold": 0.40,
        "qc_min_z_extent_ratio": 0.65,
    },
    "quality": {
        "exclude_border": True,
        "min_z_slices": 5,
        "max_volume_um3": None,
        "background_inner_um": 2.0,
        "background_outer_um": 7.0,
        "min_background_voxels": 100,
        "max_saturated_fraction": 0.01,
        "calcein_saturation_value": None,
        "pi_saturation_value": None,
    },
    "viability": {
        "mode": "uncalibrated",  # uncalibrated | controls
        "min_control_replicates": 2,
        "min_control_separation_snr": 3.0,
        "high_gate": 0.60,
        "low_gate": 0.30,
    },
    "report": {"bootstrap_iterations": 2000, "seed": 20260831,
               "save_meshes": True, "max_meshes_in_preview": 20},
    "stats": {"enabled": True, "features": ["volume_um3", "sphericity"],
              "min_replicates_per_condition": 3},
}


def merge_strict(base: dict, update: dict, prefix: str = "") -> dict:
    for key, value in update.items():
        if key not in base:
            raise ValueError(f"Unknown configuration key: {prefix}{key}")
        if isinstance(base[key], dict):
            if not isinstance(value, dict):
                raise ValueError(f"{prefix}{key} must be a mapping")
            merge_strict(base[key], value, f"{prefix}{key}.")
        else:
            base[key] = value
    return base


def load_config(path: str | Path | None = None) -> dict:
    cfg = deepcopy(DEFAULTS)
    if path:
        with Path(path).open(encoding="utf-8") as handle:
            update = yaml.safe_load(handle) or {}
        if not isinstance(update, dict):
            raise ValueError("Configuration must be a YAML mapping")
        merge_strict(cfg, update)
    validate_config(cfg)
    return cfg


def validate_config(cfg: dict) -> None:
    def number(section, key, minimum=0, maximum=None, strict=False, optional=False):
        value = cfg[section][key]
        if optional and value is None:
            return
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError(f"{section}.{key} must be a finite number")
        if value < minimum or (strict and value == minimum) or (maximum is not None and value > maximum):
            raise ValueError(f"{section}.{key} is outside its allowed range")

    if cfg["segmentation"]["method"] not in {"watershed", "probability", "labels"}:
        raise ValueError("segmentation.method must be watershed, probability, or labels")
    if cfg["segmentation"]["qc_reference_method"] not in {"none", "watershed"}:
        raise ValueError("segmentation.qc_reference_method must be none or watershed")
    if (cfg["segmentation"]["qc_reference_method"] == "watershed"
            and cfg["segmentation"]["method"] == "watershed"):
        raise ValueError("A watershed primary segmentation cannot use watershed as its QC reference")
    if cfg["segmentation"]["polarity"] not in {"bright", "dark"}:
        raise ValueError("segmentation.polarity must be bright or dark")
    threshold = cfg["segmentation"]["threshold"]
    if threshold != "otsu":
        number("segmentation", "threshold", minimum=-math.inf)
    for key in ["gaussian_sigma_um", "background_sigma_um", "closing_radius_um", "seed_min_distance_um"]:
        number("segmentation", key)
    for key in ["min_volume_um3", "seed_h_um"]:
        number("segmentation", key, strict=True)
    for key in ["probability_threshold", "max_foreground_fraction"]:
        number("segmentation", key, maximum=1, strict=True)
    number("segmentation", "qc_count_difference_threshold", maximum=1)
    number("segmentation", "qc_min_z_extent_ratio", maximum=1, strict=True)
    for section, keys in {"segmentation": ["fill_enclosed_holes", "split_touching"],
                          "quality": ["exclude_border"], "report": ["save_meshes"],
                          "stats": ["enabled"]}.items():
        for key in keys:
            if not isinstance(cfg[section][key], bool):
                raise ValueError(f"{section}.{key} must be true or false")
    for key in ["background_inner_um", "background_outer_um"]:
        number("quality", key)
    if cfg["quality"]["background_outer_um"] <= cfg["quality"]["background_inner_um"]:
        raise ValueError("Background outer radius must exceed inner radius")
    number("quality", "max_volume_um3", strict=True, optional=True)
    number("quality", "max_saturated_fraction", maximum=1)
    for key in ["calcein_saturation_value", "pi_saturation_value"]:
        number("quality", key, strict=True, optional=True)
    if cfg["viability"]["mode"] not in {"uncalibrated", "controls"}:
        raise ValueError("viability.mode must be uncalibrated or controls")
    for key in ["high_gate", "low_gate"]:
        number("viability", key, maximum=1)
    if cfg["viability"]["low_gate"] >= cfg["viability"]["high_gate"]:
        raise ValueError("Viability low_gate must be below high_gate")
    number("viability", "min_control_separation_snr", strict=True)
    for section, key, minimum in [("quality", "min_z_slices", 3), ("quality", "min_background_voxels", 1),
                                  ("viability", "min_control_replicates", 2),
                                  ("report", "bootstrap_iterations", 100),
                                  ("report", "max_meshes_in_preview", 1), ("report", "seed", 0),
                                  ("stats", "min_replicates_per_condition", 2)]:
        value = cfg[section][key]
        if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
            raise ValueError(f"{section}.{key} must be an integer >= {minimum}")
    features = cfg["stats"]["features"]
    if not isinstance(features, list) or not features or not all(isinstance(f, str) for f in features):
        raise ValueError("stats.features must be a non-empty list of feature-name strings")
    for key, value in cfg["channels"].items():
        if value is not None and (isinstance(value, bool) or not isinstance(value, int) or value < 0):
            raise ValueError(f"channels.{key} must be a nonnegative index or null")
    used = [x for x in cfg["channels"].values() if x is not None]
    if len(used) != len(set(used)):
        raise ValueError("Structure, Calcein, and PI must not reuse the same channel")
