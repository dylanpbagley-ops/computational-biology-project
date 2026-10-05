from __future__ import annotations
import csv, hashlib, json
from datetime import datetime, timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STAGE=Path(__file__).resolve().parent
IGNORE_PARTS={"__pycache__", ".cache", "cache", "data", "dependencies"}

def fail_constant(value):
    raise ValueError(f"non-finite JSON constant: {value}")
def digest(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""): h.update(block)
    return h.hexdigest()
def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"),parse_constant=fail_constant)
def check_old_manifest(label: str, manifest_rel: str, base_rel: str):
    m=load_json(ROOT/manifest_rel); failures=[]
    assert m["file_count"]==len(m["files"]),(label,"count")
    base=ROOT/base_rel
    for row in m["files"]:
        p=base/row["path"]
        if not p.is_file() or p.stat().st_size!=row["bytes"] or digest(p)!=row["sha256"]:
            failures.append(row["path"])
    return {"label":label,"listed":len(m["files"]),"checked":len(m["files"]),"mismatches":failures,"pass":not failures}
def main():
    json_files=sorted(STAGE.rglob("*.json"))
    for p in json_files: load_json(p)
    widths={}
    for p in sorted(STAGE.glob("*.tsv")):
        with p.open(encoding="utf-8",newline="") as f: rows=list(csv.reader(f,delimiter="\t"))
        assert rows and len({len(row) for row in rows})==1,(p.name,sorted({len(row) for row in rows}))
        widths[p.name]=len(rows[0])
    manifest=load_json(STAGE/"output_hashes.json")
    listed={row["path"]:row for row in manifest["files"]}
    assert len(listed)==len(manifest["files"]),"duplicate stage manifest path"
    actual={p.relative_to(STAGE).as_posix():p for p in STAGE.rglob("*") if p.is_file() and p.name!="output_hashes.json" and not any(part in IGNORE_PARTS for part in p.relative_to(STAGE).parts)}
    assert set(actual)==set(listed),("stage coverage mismatch",sorted(set(actual)-set(listed)),sorted(set(listed)-set(actual)))
    for rel,row in listed.items():
        p=actual[rel]
        assert p.stat().st_size==row["bytes"] and digest(p)==row["sha256"],rel
    before=load_json(STAGE/"preservation/before_hashes.json"); changed=[]
    for row in before["files"]:
        p=ROOT/row["path"]
        if not p.is_file() or p.stat().st_size!=row["bytes"] or digest(p)!=row["sha256"]: changed.append(row["path"])
    assert not changed,("direct old-source changes",changed)
    preservation=load_json(STAGE/"preservation/verification.json")
    assert preservation["longitudinal_preservation_verifier"]["preservation_pass"]
    assert preservation["longitudinal_preservation_verifier"]["failures"]==[]
    prior=[
      check_old_manifest("strengthening_155","strengthening_execution_Oct2/next_stage/delivery/output_hashes.json","strengthening_execution_Oct2/next_stage"),
      check_old_manifest("longitudinal_201","longitudinal_extension_Oct3/delivery/output_hashes.json","longitudinal_extension_Oct3")]
    assert all(row["pass"] for row in prior),prior
    print(json.dumps({"checked_utc":datetime.now(timezone.utc).isoformat(),"json_finite_checked":len(json_files),"tsv_columns":widths,"stage_files_hash_checked":len(listed),"stage_files_complete":True,"direct_old_files_compared":len(before["files"]),"direct_old_files_unchanged":True,"prior_output_manifests":prior,"preservation_pass":True},indent=2,allow_nan=False))
if __name__=="__main__": main()
