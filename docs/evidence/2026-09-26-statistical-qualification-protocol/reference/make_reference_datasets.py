"""Deterministic object-level datasets for the SG-5 reference comparison.

Each dataset is a nested design (objects within biological replicates within
conditions). Values are written at full precision; R (lmerTest, clubSandwich)
and the Python implementation read the same file.
"""
import csv
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
# name: (replicates per condition, object-count rule, ICC, seed)
DESIGNS = {
    "k2_g3_bal_icc30": ([3, 3], ("fixed", 10), 0.30, 1),
    "k2_g4g5_unbal_icc50": ([4, 5], ("range", 3, 25), 0.50, 2),
    "k3_g345_unbal_icc20": ([3, 4, 5], ("range", 2, 30), 0.20, 3),
    "k2_g2_bal_icc60": ([2, 2], ("fixed", 5), 0.60, 4),
    "k3_g6_unbal_icc05": ([6, 6, 6], ("range", 20, 60), 0.05, 5),
    "k2_g10_unbal_icc10": ([10, 10], ("range", 4, 80), 0.10, 6),
    "k4_g3_unbal_icc40": ([3, 3, 4, 3], ("range", 5, 15), 0.40, 7),
    "k2_g3_bal_icc00": ([3, 3], ("fixed", 8), 0.00, 8),
    "k2_g6_unbal_icc00": ([6, 7], ("range", 3, 20), 0.00, 9),
    "k3_g4_unbal_icc02": ([4, 4, 4], ("range", 5, 12), 0.02, 10),
}


def main() -> None:
    rows = []
    for name, (reps, count, icc, seed) in DESIGNS.items():
        rng = np.random.default_rng(seed)
        shift = rng.normal(0, 0.5, len(reps))
        for k, g in enumerate(reps):
            for r in range(g):
                n = count[1] if count[0] == "fixed" else int(rng.integers(count[1], count[2] + 1))
                b = rng.normal(0, np.sqrt(icc))
                for value in shift[k] + b + rng.normal(0, np.sqrt(1 - icc), n):
                    rows.append({"dataset": name, "condition": f"C{k}", "replicate": f"C{k}R{r}", "y": repr(float(value))})
    with (HERE / "reference_datasets.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["dataset", "condition", "replicate", "y"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
