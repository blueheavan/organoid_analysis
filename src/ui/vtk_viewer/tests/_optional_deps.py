"""Shared optional-dependency probes for headless browser-render tests."""
from __future__ import annotations

try:
    import playwright  # noqa: F401
    HAS_PLAYWRIGHT = True
except ImportError:
    HAS_PLAYWRIGHT = False
