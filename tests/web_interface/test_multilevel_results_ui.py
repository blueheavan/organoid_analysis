from __future__ import annotations

import json

import pandas as pd

from organoid_analysis.web_interface.multilevel_results import (
    RESULT_TABLES,
    load_multilevel_result_tables,
)


def test_load_multilevel_result_tables(tmp_path):
    for relative, _ in RESULT_TABLES.values():
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"object_id": [1]}).to_parquet(path, index=False)
    summary = {"organoid_count": 1, "cell_count": 2, "nucleus_count": 3}
    summary_path = tmp_path / "summary" / "analysis_summary.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary), encoding="utf-8")
    loaded_summary, tables = load_multilevel_result_tables(tmp_path)
    assert loaded_summary == summary
    assert set(tables) == set(RESULT_TABLES)
