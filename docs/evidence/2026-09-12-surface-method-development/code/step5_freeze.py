"""Emit freeze_record.json: SHA-256 of every frozen file, plus the frozen
constants read back from the code itself (not retyped)."""
import hashlib, json, platform, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path

import numpy, scipy, skimage

import common, surface_area_v2 as v2, step6_vv_run as vv

HERE = Path(__file__).resolve().parent
FILES = ["common.py", "estimators_v2.py", "surface_area_v2.py", "step6_vv_run.py",
         "step1_baseline.py", "step2_mechanism.py", "step3_stencil_sweep.py",
         "step4_dev_explore.py", "out/FREEZE_RECORD.md", "out/CANDIDATES.md",
         "out/MECHANISM.md"]

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

rec = {
    "frozen_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    "confirmation_grid_executed_before_this_record": False,
    "estimator": {
        "method_name": v2.METHOD_NAME,
        "scheme": v2.SCHEME,
        "m_candidates": list(v2.M_CANDIDATES),
        "domain_rho_eff_min": v2.DOMAIN_RHO_EFF_MIN,
        "domain_aniso_max": v2.DOMAIN_ANISO_MAX,
    },
    "acceptance": {
        "area": vv.AREA_CRITERION,
        "volume": vv.VOLUME_CRITERION,
        "sphericity": vv.SPHERICITY_CRITERION,
        "rule": "every in-domain case must satisfy all three; one violation = FAIL",
    },
    "legacy_estimator_preserved": {
        "name": "marching_cubes_binary_lewiner_v1",
        "action": "kept as default; candidate added alongside, never in place of",
        "evidence": "out/E0_legacy_baseline.json",
    },
    "grids": {g: {k: (list(v) if isinstance(v, list) else v)
                  for k, v in common.GRIDS[g].items()} for g in ("dev", "confirm")},
    "sha256": {f: sha(HERE / f) for f in FILES},
    "environment": {
        "python": sys.version.split()[0], "platform": platform.platform(),
        "numpy": numpy.__version__, "scipy": scipy.__version__,
        "scikit_image": skimage.__version__,
        "repo_head": subprocess.run(["git", "rev-parse", "HEAD"], cwd=common.REPO,
                                    capture_output=True, text=True).stdout.strip(),
    },
}
Path(HERE / "out/freeze_record.json").write_text(json.dumps(rec, indent=2, default=str))
print(json.dumps({k: rec[k] for k in ("frozen_at_utc", "estimator", "acceptance")}, indent=2))
print("sha256 entries:", len(rec["sha256"]))
