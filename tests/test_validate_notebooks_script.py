"""Regression test for scripts/validate_notebooks.py (P3-1 audit finding).

Notebooks are optional in this repository (``notebooks/`` currently holds
none). Before this fix, the script exited 1 when it found none, which made
`pixi run check` fail and -- because pixi's `ci` depends-on chain stops at
the first failure -- silently prevented `pixi run ci` from ever reaching
`typecheck`, without any actual test or validation failing.

The script derives its own project root from ``__file__`` rather than the
process cwd, so this test runs it against the real repository rather than an
isolated fixture tree.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_ROOT / "scripts" / "validate_notebooks.py"


def test_no_tracked_notebooks_still_exits_zero():
    notebooks = sorted((PROJECT_ROOT / "notebooks").glob("*.ipynb"))
    assert notebooks == [], (
        "This regression test assumes no notebooks are currently tracked; "
        "if that changes, this test (and the audited P3-1 behavior) should "
        "be revisited rather than silently skipped."
    )
    result = subprocess.run(
        [sys.executable, str(SCRIPT)],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stderr
