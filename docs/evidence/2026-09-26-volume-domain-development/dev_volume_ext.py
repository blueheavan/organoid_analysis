"""Development EXTENSION (2026-09-26), run after dev_volume.py, before any confirmation.

dev_volume.py's pre-declared rule selected rho_vol_min = 36, the TOP of its
ladder: the threshold would sit exactly where coverage ends, the error the
2026-09-12 round-1 record documents (FREEZE_RECORD_V3.md section 8.2). This
extension samples the two binding classes (torus, capsule) above 36 with the
same spacings and placement scheme, so the declared threshold is interpolated
inside measured development evidence. The rule is unchanged; this only checks
that it does not rest on the last ladder value. Development data only.
"""
from __future__ import annotations

import json
import sys

import dev_volume as dv

dv.LADDER = (36, 40, 44, 48)
dv.SHAPES = ("torus", "capsule")
dv.SEED = 20260926 + 1

if __name__ == "__main__":
    rows = dv.run_grid(dv.HERE / "volume_dev_ext.csv")
    gated = [r for r in rows if r["rho_in"] >= 10 and r["anisotropy"] <= 4]
    out = {"n": len(gated), "by_level": {}}
    for level in dv.LADDER:
        sel = [r for r in gated if r["rho_vol"] >= level]
        w = max(sel, key=lambda r: abs(r["volume_rel_err"]))
        out["by_level"][level] = {"n": len(sel), "worst_abs": abs(w["volume_rel_err"]), "worst_case": w["case"],
                                  "worst_shape": w["shape"], "worst_spacing": w["spacing"], "worst_pose": w["pose"]}
    out["max_abs_all"] = float(max(abs(r["volume_rel_err"]) for r in gated))
    out["within_margin"] = bool(out["max_abs_all"] <= dv.MARGIN_LIMIT)
    (dv.HERE / "volume_dev_ext_summary.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    sys.exit(0)
