"""Known 12-bit clipping in uint16 storage, for audit C8."""
import copy
import json
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

from organoid_analysis.config import load_config
from organoid_analysis.quantification.features import marker_measurements

labels = np.zeros((31, 31, 31), np.uint32)
labels[12:19, 12:19, 12:19] = 1
calcein = np.full(labels.shape, 100, np.uint16)
pi = np.full(labels.shape, 100, np.uint16)
calcein[labels > 0] = 4095
pi[labels > 0] = 200
bbox = ndi.find_objects(labels)[0]
cfg = load_config()["quality"]
explicit = copy.deepcopy(cfg)
explicit.update(calcein_saturation_value=4095, pi_saturation_value=4095)
result = {
    "detector_bit_depth": 12, "storage_dtype": "uint16", "known_detector_ceiling": 4095,
    "true_calcein_saturated_fraction": 1.0,
    "default": marker_measurements(labels, 1, bbox, {"calcein": calcein, "pi": pi}, (1., 1., 1.), cfg),
    "explicit_acquisition_ceiling": marker_measurements(labels, 1, bbox, {"calcein": calcein, "pi": pi}, (1., 1., 1.), explicit),
}
text = json.dumps(result, indent=2, allow_nan=False)
Path(__file__).with_suffix(".json").write_text(text + "\n")
print(text)
