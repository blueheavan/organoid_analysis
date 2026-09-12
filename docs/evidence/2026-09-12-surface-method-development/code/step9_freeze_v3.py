"""Hash every frozen artefact immediately before the confirmation run."""
import hashlib, json, pathlib, sys

FILES = ["common.py", "common2.py", "estimators_v2.py", "estimators_v3.py",
         "surface_area_v3.py", "step8_dev2_run.py", "step10_confirm2_run.py",
         "step9_freeze_v3.py", "out/FREEZE_RECORD_V3.md"]


def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def snapshot():
    h = {f: sha(f) for f in FILES}
    w = sorted(pathlib.Path("out/crofton_weights_v3").glob("*.json"))
    h["out/crofton_weights_v3/"] = hashlib.sha256(
        b"".join(p.read_bytes() for p in w)).hexdigest()
    h["_n_weight_files"] = len(w)
    return h


if __name__ == "__main__":
    cur = snapshot()
    out = pathlib.Path("out/freeze_record_v3.json")
    if len(sys.argv) > 1 and sys.argv[1] == "--verify":
        old = json.load(open(out))
        diff = [k for k in old if k != "_n_weight_files" and old[k] != cur.get(k)]
        print(json.dumps({"V3_hashes_match": not diff, "changed": diff}, indent=1))
    else:
        json.dump(cur, open(out, "w"), indent=1)
        print(json.dumps(cur, indent=1))
