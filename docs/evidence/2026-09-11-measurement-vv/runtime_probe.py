"""Bounded runtime/robustness probes for the Web mask-feature path (P2 section 7).

    pixi run python docs/evidence/2026-09-11-measurement-vv/runtime_probe.py

These are engineering observations on synthetic labels. They are not accuracy
evidence and not a qualification of arbitrary production volumes.
"""
from __future__ import annotations

import json
import resource
import sys
import time
import zipfile
import io
from pathlib import Path

import numpy as np

from organoid_analysis.quantification.mask_features import extract_mask_features
from organoid_analysis.result_export.mask_feature_bundle import build_mask_feature_bundle

OUT = Path(__file__).resolve().parent
results: dict = {}


def peak_rss_mb() -> float:
    # macOS reports ru_maxrss in bytes, Linux in KiB.
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value / 2**20 if sys.platform == "darwin" else value / 2**10


def spheres(shape: tuple[int, int, int], n: int, radius: int, seed: int, dtype=np.uint32) -> np.ndarray:
    rng = np.random.default_rng(seed)
    labels = np.zeros(shape, dtype)
    ball = np.indices((2 * radius + 1,) * 3) - radius
    ball = (ball ** 2).sum(0) <= radius ** 2
    placed = 0
    for _ in range(n * 20):
        if placed == n:
            break
        c = [int(rng.integers(radius + 1, s - radius - 1)) for s in shape]
        sl = tuple(slice(ci - radius, ci + radius + 1) for ci in c)
        if labels[sl][ball].any():
            continue
        labels[sl][ball] = placed + 1
        placed += 1
    return labels


# 1. Large label volume and object count (Web mask-feature path).
rss0 = peak_rss_mb()
big = spheres((64, 1024, 1024), 3000, 4, 1)
t0 = time.perf_counter()
features = extract_mask_features(big, (0.5, 0.5, 2.0))
results["large_volume"] = {
    "shape_zyx": list(big.shape), "dtype": str(big.dtype), "array_mb": big.nbytes / 2**20,
    "n_objects_placed": int(big.max()), "n_rows": int(len(features)),
    "seconds": time.perf_counter() - t0, "peak_rss_mb_before": rss0, "peak_rss_mb_after": peak_rss_mb(),
    "all_rows_finite_volume": bool(np.isfinite(features.volume_um3).all()),
}
del big, features

# 2. Sparse uint32 identities near the dtype limit are preserved.
sparse = np.zeros((6, 8, 8), np.uint32)
ids = [1, 2**31 + 7, 2**32 - 2]
for k, i in enumerate(ids):
    sparse[1:4, 1:3, 1 + 2 * k:2 + 2 * k] = i
f = extract_mask_features(sparse)
results["sparse_uint32_ids"] = {"input_ids": ids, "output_ids": [int(x) for x in f.index],
                                "preserved": [int(x) for x in f.index] == ids}

# 3. uint16 and uint32 give identical measurements for the same labels.
small = spheres((20, 40, 40), 6, 3, 2, np.uint16)
f16 = extract_mask_features(small, (0.4, 0.4, 1.2))
f32 = extract_mask_features(small.astype(np.uint32), (0.4, 0.4, 1.2))
results["uint16_vs_uint32_identical"] = bool(f16.equals(f32))

# 4. Extreme anisotropy (Z/XY = 100): computed without warning; accuracy unqualified.
# 13 planes with the centre at plane 6 keep a background plane on both Z ends
# (an earlier 9-plane draft touched the Z border and was not a fair probe).
z, y, x = np.indices((13, 211, 211))
r = 10.0
ball = ((z - 6) * 2.5) ** 2 + ((y - 105) * .1) ** 2 + ((x - 105) * .1) ** 2 <= r ** 2
fa = extract_mask_features(ball.astype(np.uint8), (0.1, 0.1, 2.5))
results["extreme_anisotropy_25x"] = {
    "spacing_xyz_um": [0.1, 0.1, 2.5], "n_z_slices_with_object": int(ball.any(axis=(1, 2)).sum()),
    "area_rel_err": float(fa.surface_area_um2.iloc[0] / (4 * np.pi * r * r) - 1),
    "volume_rel_err": float(fa.volume_um3.iloc[0] / (4 / 3 * np.pi * r ** 3) - 1),
    "qc_flags": str(fa.qc_flags.iloc[0]), "warning_emitted": False,
}

# 5. Degenerate objects: single voxel, single plane, single line.
deg = np.zeros((5, 7, 7), np.uint8)
deg[2, 3, 3] = 1
deg[1, 1:4, 1:4] = 2
deg[3, 5, 1:6] = 3
fd = extract_mask_features(deg, (1, 1, 1))
results["degenerate_objects"] = fd[["volume_um3", "surface_area_um2", "sphericity", "solidity", "qc_flags"]] \
    .reset_index().to_dict("records")

# 6. Empty mask: explicit empty table.
results["empty_mask_rows"] = int(len(extract_mask_features(np.zeros((4, 4, 4), np.uint8))))

# 7. Export determinism: identical inputs produce byte-identical bundles.
cfg = {"xy_spacing_source": "metadata", "anisotropy_source": "metadata"}
b1 = build_mask_feature_bundle(small, f16, spacing_um=(0.4, 0.4, 1.2), object_type="nucleus",
                               fill_holes=False, segmentation_config=cfg)
b2 = build_mask_feature_bundle(small, extract_mask_features(small, (0.4, 0.4, 1.2)), spacing_um=(0.4, 0.4, 1.2),
                               object_type="nucleus", fill_holes=False, segmentation_config=cfg)
prov = json.loads(zipfile.ZipFile(io.BytesIO(b1)).read("measurement_provenance.json"))
results["export_bytes_identical"] = b1 == b2
results["export_surface_area_method"] = prov.get("surface_area_method")

(OUT / "runtime_probe.json").write_text(json.dumps(results, indent=2, default=str))
print(json.dumps(results, indent=2, default=str))
