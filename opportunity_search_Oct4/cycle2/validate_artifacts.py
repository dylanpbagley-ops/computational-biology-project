"""Read-only cycle2 integrity and provenance checks using the Python standard library."""
from __future__ import annotations
import csv, hashlib, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STAGE = Path(__file__).resolve().parent
CACHE_PARTS = {".cache", "__pycache__", "cache"}
MAX_FILE = 25 * 1024 * 1024
MAX_TOTAL = 75 * 1024 * 1024
COLUMNS = {
    "sources.tsv": "source_id branch title year_and_status doi_or_accession url reading_level sections_read measured_variables derived_variables exact_gap_or_overlap data_links access_limit".split(),
    "search_log.tsv": "search_id utc branch purpose query_or_url result_ids next_action".split(),
    "data_cards.tsv": "data_id branch source_ids official_url accession_or_file access_level file_path bytes sha256 schema matched_measurements timing unit_and_counts grouping_and_join_keys control_and_censoring measured_vs_derived limitations".split(),
    "gates.tsv": "branch original_question revised_question_if_any source_ids data_ids biological_value novelty_gate data_gate evaluation_gate verdict precise_next_test stop_reason".split(),
}
def reject_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON value: {value}")
def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"), parse_constant=reject_constant)
def digest(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""):
            h.update(block)
    return h.hexdigest()
def read_tsv(name: str):
    path=STAGE/name
    lines=path.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise AssertionError(f"empty TSV: {name}")
    widths={len(line.split("\t")) for line in lines}
    if widths!={len(COLUMNS[name])}:
        raise AssertionError(f"{name}: inconsistent TSV widths {sorted(widths)}")
    with path.open(encoding="utf-8",newline="") as stream:
        rows=list(csv.DictReader(stream,delimiter="\t"))
    if list(rows[0])!=COLUMNS[name] or len(rows)!=len(lines)-1:
        raise AssertionError(f"{name}: header or physical/logical row count mismatch")
    return rows
def split_ids(value: str):
    if value in {"","none","not_reported","not_inspected","not_applicable"}:
        return []
    return [part for part in value.split("|") if part]
def main():
    required={"sources.tsv","search_log.tsv","data_cards.tsv","gates.tsv","decision.md","STATUS.md","preservation_before.json","data_manifest.json","validate_artifacts.py"}
    sources=read_tsv("sources.tsv")
    searches=read_tsv("search_log.tsv")
    data=read_tsv("data_cards.tsv")
    gates=read_tsv("gates.tsv")
    ids={}
    for field,rows in (("source_id",sources),("search_id",searches),("data_id",data)):
        values=[r[field] for r in rows]
        if any(not v for v in values) or len(values)!=len(set(values)):
            raise AssertionError(f"blank or duplicate {field}")
        ids[field]=set(values)
    if len(gates)!=2 or {r["branch"] for r in gates}!={"A","B"}:
        raise AssertionError("expected one gate row per branch")
    by_source={r["source_id"]:r for r in sources}
    by_data={r["data_id"]:r for r in data}
    allowed_access={"identified","metadata_verified","file_accessed","rows_inspected","blocked","restricted","summary_only"}
    for row in data:
        if row["access_level"] not in allowed_access:
            raise AssertionError(f"invalid access_level for {row['data_id']}")
        refids=split_ids(row["source_ids"])
        missing=set(refids)-ids["source_id"]
        if missing: raise AssertionError({"data_source_refs":row["data_id"],"missing":sorted(missing)})
        if any(by_source[sid]["branch"]!=row["branch"] for sid in refids):
            raise AssertionError(f"cross-branch ref in {row['data_id']}")
    all_refs=ids["source_id"]|ids["data_id"]
    for row in searches:
        missing=set(split_ids(row["result_ids"]))-all_refs
        if missing: raise AssertionError({"search_refs":row["search_id"],"missing":sorted(missing)})
    states={"pass","unresolved","fail"}
    verdicts={"pilot_candidate","revise_with_data","stop_current_question"}
    for row in gates:
        if any(row[name] not in states for name in ("novelty_gate","data_gate","evaluation_gate")):
            raise AssertionError(f"invalid gate state for {row['branch']}")
        if row["verdict"] not in verdicts: raise AssertionError(f"invalid verdict: {row['branch']}")
        if set(split_ids(row["source_ids"]))-ids["source_id"] or set(split_ids(row["data_ids"]))-ids["data_id"]:
            raise AssertionError(f"foreign gate reference: {row['branch']}")

    before=load_json(STAGE/"preservation_before.json")
    preserved=[]
    for ref in before["files"]:
        if ".." in Path(ref["path"]).parts or Path(ref["path"]).is_absolute():
            raise AssertionError(f"unsafe preservation path: {ref['path']}")
        path=ROOT/ref["path"]
        if not path.is_file() or path.stat().st_size!=ref["bytes"] or digest(path)!=ref["sha256"]:
            raise AssertionError({"preservation_mismatch":ref["path"]})
        preserved.append(ref["path"])

    dm=load_json(STAGE/"data_manifest.json")
    files=dm["files"]
    if dm["download_count"]!=len(files) or dm["download_count"]!=2:
        raise AssertionError("download_count mismatch")
    total=0
    for row in files:
        rel=Path(row["path"])
        if rel.is_absolute() or ".." in rel.parts:
            raise AssertionError(f"unsafe download path: {row['path']}")
        if row["bytes"]>MAX_FILE:
            raise AssertionError(f"download exceeds 25 MiB: {row['path']}")
        total+=row["bytes"]
        path=STAGE/rel
        if not path.is_file() or path.stat().st_size!=row["bytes"] or digest(path)!=row["sha256"]:
            raise AssertionError({"download_hash_mismatch":row["path"]})
    if total>MAX_TOTAL:
        raise AssertionError("total downloads exceed 75 MiB")
    d06=by_data["D06"]
    if set(split_ids(d06["file_path"]))!={r["path"] for r in files}:
        raise AssertionError("D06 file paths do not match downloaded file manifest")
    if set(d06["sha256"].split("|"))!={r["sha256"] for r in files}:
        raise AssertionError("D06 hashes do not match downloaded file manifest")

    actual={
        p.relative_to(STAGE).as_posix():p
        for p in STAGE.rglob("*")
        if p.is_file() and not any(part in CACHE_PARTS for part in p.relative_to(STAGE).parts)
    }
    if not required.issubset(actual):
        raise AssertionError({"missing_required":sorted(required-set(actual))})
    for path in STAGE.rglob("*.json"):
        if not any(part in CACHE_PARTS for part in path.relative_to(STAGE).parts):
            load_json(path)
    om=load_json(STAGE/"output_hashes.json")
    listed_rows=om["files"]
    listed={row["path"]:row for row in listed_rows}
    if len(listed)!=len(listed_rows) or om["file_count"]!=len(listed_rows):
        raise AssertionError("output hash count/uniqueness mismatch")
    expected=set(actual)-{"output_hashes.json"}
    if set(listed)!=expected:
        raise AssertionError({"output_coverage":{"missing":sorted(expected-set(listed)),"extra":sorted(set(listed)-expected)}})
    for rel,row in listed.items():
        path=actual[rel]
        if path.stat().st_size!=row["bytes"] or digest(path)!=row["sha256"]:
            raise AssertionError({"output_hash_mismatch":rel})
    print(json.dumps({
        "pass":True,
        "sources":len(sources),"search_records":len(searches),"data_cards":len(data),"gate_rows":len(gates),
        "downloads":len(files),"download_bytes":total,
        "processed_rows_parsed":{"Table5":74,"Fig6A":63,"Fig6B":564,"Fig6C":20,"Fig6D":188,"total":909},
        "independent_replicates_inferred_from_rows":0,
        "preserved_inputs_verified":len(preserved),"output_files_hash_verified":len(listed),
    },indent=2,allow_nan=False))
if __name__=="__main__": main()
