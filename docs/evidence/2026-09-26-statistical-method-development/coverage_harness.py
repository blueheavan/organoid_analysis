"""Type-I error and 95% CI coverage of the frozen contrast method (SG-5; D-10).

Data are generated from the random-intercept model on the analysis scale with
total variance 1: replicate effects with variance ICC, object residuals with
variance 1 - ICC. Because the fit uses replicate sufficient statistics exactly,
each dataset is drawn as those statistics (replicate means and the pooled
within-replicate sum of squares, chi-square distributed), which is
distributionally identical to drawing every object.

Per cell: ``n_sim`` null datasets (all condition means equal) give the type-I
error of the contrast condition 1 vs condition 0 at alpha 0.05 and of the
omnibus F; ``n_sim`` further datasets with condition 1 shifted by EFFECT give
the coverage of the 95% interval for that contrast. Designs are redrawn for
every dataset. Branch use (LMM vs CR2 fallback) is counted and reported per
branch.

Usage: python coverage_harness.py --seed S --n-sim N --out results.csv [--processes P] [--robustness]
"""
from __future__ import annotations

import argparse
import csv
import itertools
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))

from organoid_analysis.statistics import small_sample as ss  # noqa: E402

GRID = {"K": (2, 3), "G": (2, 3, 4, 6, 10, 20), "m": (5, 40, 150),
        "balance": ("balanced", "unbalanced"), "icc": (0.0, 0.1, 0.3, 0.7)}
ROBUSTNESS = [{"K": 2, "G": g, "m": 40, "balance": "balanced", "icc": icc, "dist": "skewed_b"}
              for g in (3, 4, 6, 10, 20) for icc in (0.3, 0.7)]
EFFECT = 0.5
ALPHA = 0.05
UNBALANCED_SIGMA = 0.6  # log-normal spread of objects per replicate in unbalanced cells


def cells(robustness: bool) -> list[dict]:
    if robustness:
        return ROBUSTNESS
    return [dict(zip(GRID, values), dist="normal") for values in itertools.product(*GRID.values())]


def _design(cell: dict, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    k, g, m = cell["K"], cell["G"], cell["m"]
    unbalanced = cell["balance"] == "unbalanced"
    per_condition = [g + (c % 2 if unbalanced else 0) for c in range(k)]
    group = np.repeat(np.arange(k), per_condition)
    if unbalanced:
        n = np.maximum(2, np.rint(m * rng.lognormal(0.0, UNBALANCED_SIGMA, group.size)))
    else:
        n = np.full(group.size, float(m))
    return n, group


def _replicate_effects(cell: dict, size: int, rng: np.random.Generator) -> np.ndarray:
    if cell["dist"] == "skewed_b":
        raw = rng.lognormal(0.0, 1.0, size)
        raw = (raw - np.exp(0.5)) / np.sqrt((np.exp(1.0) - 1.0) * np.exp(1.0))
        return raw * np.sqrt(cell["icc"])
    return rng.normal(0.0, np.sqrt(cell["icc"]), size)


def _dataset(cell: dict, rng: np.random.Generator, shift: float) -> ss.ClusterData:
    n, group = _design(cell, rng)
    sigma_e2 = 1.0 - cell["icc"]
    mu = np.zeros(cell["K"])
    mu[1] = shift
    mean = mu[group] + _replicate_effects(cell, group.size, rng) + rng.normal(0.0, np.sqrt(sigma_e2 / n))
    ssw = sigma_e2 * rng.chisquare(n.sum() - group.size)
    return ss.ClusterData(n=n, mean=mean, group=group, ssw=float(ssw), n_conditions=cell["K"])


def run_cell(args: tuple[dict, int, int]) -> dict:
    cell, n_sim, seed = args
    key = [seed, cell["K"], cell["G"], cell["m"], GRID["balance"].index(cell["balance"]),
           int(round(cell["icc"] * 1000)), ["normal", "skewed_b"].index(cell["dist"])]
    rng = np.random.default_rng(np.random.SeedSequence(key))
    contrast = np.zeros(cell["K"])
    contrast[1], contrast[0] = 1.0, -1.0
    counts = {"null_n": 0, "null_reject": 0, "omnibus_reject": 0, "alt_n": 0, "alt_cover": 0, "failures": 0,
              "null_fallback": 0, "null_reject_fallback": 0, "alt_fallback": 0, "alt_cover_fallback": 0}
    for _ in range(n_sim):
        try:
            model = ss.fit(_dataset(cell, rng, 0.0))
            result = ss.contrast(model, contrast)
            omnibus_p = ss.omnibus(model)[3]
        except (ValueError, np.linalg.LinAlgError):
            counts["failures"] += 1
            continue
        fallback = model.branch == ss.FALLBACK_BRANCH
        reject = bool(result.p_value < ALPHA)
        counts["null_n"] += 1
        counts["null_reject"] += reject
        counts["omnibus_reject"] += bool(omnibus_p < ALPHA)
        counts["null_fallback"] += fallback
        counts["null_reject_fallback"] += reject and fallback
    for _ in range(n_sim):
        try:
            model = ss.fit(_dataset(cell, rng, EFFECT))
            result = ss.contrast(model, contrast)
        except (ValueError, np.linalg.LinAlgError):
            counts["failures"] += 1
            continue
        fallback = model.branch == ss.FALLBACK_BRANCH
        cover = bool(result.ci_low <= EFFECT <= result.ci_high)
        counts["alt_n"] += 1
        counts["alt_cover"] += cover
        counts["alt_fallback"] += fallback
        counts["alt_cover_fallback"] += cover and fallback
    return {**cell, "n_sim": n_sim, "seed": seed, **counts,
            "type1": counts["null_reject"] / max(counts["null_n"], 1),
            "coverage": counts["alt_cover"] / max(counts["alt_n"], 1),
            "omnibus_type1": counts["omnibus_reject"] / max(counts["null_n"], 1),
            "fallback_fraction": counts["null_fallback"] / max(counts["null_n"], 1)}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--n-sim", type=int, required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--processes", type=int, default=4)
    parser.add_argument("--robustness", action="store_true")
    args = parser.parse_args(argv)
    jobs = [(cell, args.n_sim, args.seed) for cell in cells(args.robustness)]
    print(f"{len(jobs)} cells x {args.n_sim} x 2 datasets", flush=True)
    t0 = time.time()
    from multiprocessing import Pool
    with Pool(args.processes) as pool, open(args.out, "w", newline="", encoding="utf-8") as handle:
        writer = None
        for k, row in enumerate(pool.imap(run_cell, jobs)):
            if writer is None:
                writer = csv.DictWriter(handle, fieldnames=list(row), lineterminator="\n")
                writer.writeheader()
            writer.writerow({key: (f"{value:.10g}" if isinstance(value, float) else value) for key, value in row.items()})
            handle.flush()
            if k % 20 == 0:
                print(f"{k:4d}/{len(jobs)} {time.time() - t0:7.0f}s", flush=True)
    print(f"done in {time.time() - t0:.0f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
