"""Regression guard: exploratory statistics must stay out of reported inference.

`statistics/exploration.py` is exploratory by name (normality pre-checks,
supervised classifiers). A reported p-value or interval may not be sourced from
it. This is a source-level contract so a future edit cannot quietly import it
into the inferential path.
"""
from __future__ import annotations

from pathlib import Path

STATISTICS_DIR = Path(__file__).resolve().parents[2] / "src" / "organoid_analysis" / "statistics"


def test_inferential_modules_do_not_import_exploration() -> None:
    offenders: list[str] = []
    for name in ("inference.py", "aggregation.py"):
        text = (STATISTICS_DIR / name).read_text(encoding="utf-8")
        for line in text.splitlines():
            code = line.split("#", 1)[0]
            if "import" in code and "exploration" in code:
                offenders.append(f"{name}: {line.strip()}")
    assert not offenders, (
        "exploratory statistics must not enter the reported inference path: "
        + "; ".join(offenders)
    )
