"""Development: re-estimate the development cells that deviated most, at n_sim = 20000.

Run after dev_coverage.csv was inspected; development only (different seed).
"""
import csv
import sys
from multiprocessing import Pool

import coverage_harness as h

CELLS = [
    {"K": 3, "G": 4, "m": 5, "balance": "unbalanced", "icc": 0.0},
    {"K": 2, "G": 4, "m": 5, "balance": "unbalanced", "icc": 0.0},
    {"K": 3, "G": 4, "m": 150, "balance": "unbalanced", "icc": 0.0},
    {"K": 3, "G": 6, "m": 150, "balance": "unbalanced", "icc": 0.0},
    {"K": 3, "G": 20, "m": 40, "balance": "balanced", "icc": 0.1},
    {"K": 2, "G": 10, "m": 5, "balance": "unbalanced", "icc": 0.1},
    {"K": 3, "G": 10, "m": 40, "balance": "unbalanced", "icc": 0.0},
    {"K": 2, "G": 3, "m": 40, "balance": "unbalanced", "icc": 0.0},
    {"K": 3, "G": 3, "m": 5, "balance": "unbalanced", "icc": 0.1},
]

if __name__ == "__main__":
    jobs = [({**c, "dist": "normal"}, 20000, 202) for c in CELLS]
    with Pool(4) as pool, open("dev_targeted.csv", "w", newline="") as handle:
        writer = None
        for row in pool.imap(h.run_cell, jobs):
            writer = writer or csv.DictWriter(handle, fieldnames=list(row), lineterminator="\n")
            if handle.tell() == 0:
                writer.writeheader()
            writer.writerow(row)
            print({k: row[k] for k in ("K", "G", "m", "balance", "icc", "type1", "coverage", "fallback_fraction")},
                  flush=True)
    sys.exit(0)
