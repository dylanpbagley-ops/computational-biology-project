"""Read-only integrity, scope, and preservation checks for cycle3."""
from __future__ import annotations
import csv, hashlib, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STAGE = Path(__file__).resolve().parent
CACHE_PARTS = {".cache", "__pycache__", "cache"}
MAX_FILE = 25 * 1024 * 1024
MAX_TOTAL = 75 * 1024 * 1024
COLUMNS = {
    "sources.tsv": "source_id title year_and_status doi_or_accession url reading_level sections_read measured_variables derived_variables prior_art_overlap exact_unanswered_contrast data_links limits".split(),
    "search_log.tsv": "search_id utc purpose query_or_url source_ids next_action".split(),
    "data_cards.tsv": "data_id source_ids official_url access_level file_path bytes sha256 schema stimulus_and_timing units_and_counts join_keys paired_controls measured_vs_derived limits".split(),
    "title_screen.tsv": "screen_id query_id query_page_rank pmid title doi pmcid screen_level candidate_source_id".split(),
}

def reject_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON value: {value}")

def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"), parse_constant=reject_constant)

def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def read_tsv(name: str):
    path = STAGE / name
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows or list(rows[0]) != COLUMNS[name]:
        raise AssertionError(f"{name}: empty rows or wrong header")
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise AssertionError(f"{name}: row width differs from header")
    return rows

def split_ids(value: str):
    if value in {"", "none", "not_reported", "not_inspected", "not_applicable"}:
        return []
    return [part for part in value.split("|") if part]

def main():
    required = {"sources.tsv", "search_log.tsv", "data_cards.tsv", "title_screen.tsv", "decision.md", "STATUS.md", "preservation_before.json", "data_manifest.json", "validate_artifacts.py", "write_output_hashes.py", "output_hashes.json", "P05_data_inspection.md"}
    sources = read_tsv("sources.tsv")
    searches = read_tsv("search_log.tsv")
    cards = read_tsv("data_cards.tsv")
    screen = read_tsv("title_screen.tsv")
    for field, rows in (("source_id", sources), ("search_id", searches), ("data_id", cards)):
        vals = [row[field] for row in rows]
        if any(not v for v in vals) or len(vals) != len(set(vals)):
            raise AssertionError(f"blank or duplicate {field}")
    source_ids = {row["source_id"] for row in sources}
    search_ids = {row["search_id"] for row in searches}
    if len(searches) > 12:
        raise AssertionError("original searches plus two correction queries exceed the 12-record cap")
    for row in cards:
        missing = set(split_ids(row["source_ids"])) - source_ids
        if missing:
            raise AssertionError({"data_source_refs": row["data_id"], "missing": sorted(missing)})
    for row in searches:
        missing = set(split_ids(row["source_ids"])) - source_ids
        if missing:
            raise AssertionError({"search_source_refs": row["search_id"], "missing": sorted(missing)})
    if set(search_ids) != {f"Q{i:02d}" for i in range(1, len(searches) + 1)}:
        raise AssertionError("search IDs are not a contiguous Q01..Qnn sequence")
    screen_ids = [row["screen_id"] for row in screen]
    screen_pmids = [row["pmid"] for row in screen]
    if len(screen) != 40 or len(set(screen_ids)) != 40 or len(set(screen_pmids)) != 40:
        raise AssertionError("expected exactly 40 unique title-screen rows")
    if screen_ids != [f"T{i:03d}" for i in range(1, 41)]:
        raise AssertionError("title-screen IDs must be contiguous T001..T040")
    if sum(row["query_id"] == "Q11" for row in screen) != 35 or sum(row["query_id"] == "Q12" for row in screen) != 5:
        raise AssertionError("title screen must contain Q11 first 35 and Q12 first 5")
    if any(row["screen_level"] != "title_only" for row in screen):
        raise AssertionError("title screen depth mislabeled")
    if any(row["candidate_source_id"] not in source_ids | {"none"} for row in screen):
        raise AssertionError("title screen candidate source missing")
    before = load_json(STAGE / "preservation_before.json")
    preserved = []
    for ref in before["files"]:
        rel = Path(ref["path"])
        if rel.is_absolute() or ".." in rel.parts:
            raise AssertionError(f"unsafe preservation path: {ref['path']}")
        path = ROOT / rel
        if not path.is_file() or path.stat().st_size != ref["bytes"] or digest(path) != ref["sha256"]:
            raise AssertionError({"preservation_mismatch": ref["path"]})
        preserved.append(ref["path"])
    dm = load_json(STAGE / "data_manifest.json")
    files = dm["files"]
    if dm["download_count"] != len(files):
        raise AssertionError("download_count does not match manifest file list")
    total = 0
    listed_downloads = set()
    for row in files:
        rel = Path(row["path"])
        if rel.is_absolute() or ".." in rel.parts or not rel.parts or rel.parts[0] != "data":
            raise AssertionError(f"unsafe/non-data download path: {row['path']}")
        path = STAGE / rel
        if not path.is_file() or path.stat().st_size != row["bytes"] or digest(path) != row["sha256"]:
            raise AssertionError({"download_hash_mismatch": row["path"]})
        if row["bytes"] > MAX_FILE:
            raise AssertionError(f"download exceeds 25 MiB: {row['path']}")
        total += row["bytes"]
        listed_downloads.add(rel.as_posix())
    if total > MAX_TOTAL or dm["total_download_bytes"] != total:
        raise AssertionError("download byte total mismatch or above 75 MiB")
    actual_downloads = {p.relative_to(STAGE).as_posix() for p in (STAGE / "data").rglob("*") if p.is_file()} if (STAGE / "data").exists() else set()
    if actual_downloads != listed_downloads:
        raise AssertionError({"unlisted_or_missing_downloads": {"unlisted": sorted(actual_downloads-listed_downloads), "missing": sorted(listed_downloads-actual_downloads)}})
    actual = {p.relative_to(STAGE).as_posix(): p for p in STAGE.rglob("*") if p.is_file() and not any(part in CACHE_PARTS for part in p.relative_to(STAGE).parts)}
    if not required.issubset(actual):
        raise AssertionError({"missing_required": sorted(required-set(actual))})
    for path in STAGE.rglob("*.json"):
        if not any(part in CACHE_PARTS for part in path.relative_to(STAGE).parts):
            load_json(path)
    om = load_json(STAGE / "output_hashes.json")
    listed_rows = om["files"]
    listed = {row["path"]: row for row in listed_rows}
    if len(listed) != len(listed_rows) or om["file_count"] != len(listed_rows):
        raise AssertionError("output hash count/uniqueness mismatch")
    expected = set(actual) - {"output_hashes.json"}
    if set(listed) != expected:
        raise AssertionError({"output_coverage": {"missing": sorted(expected-set(listed)), "extra": sorted(set(listed)-expected)}})
    for rel, row in listed.items():
        path = actual[rel]
        if path.stat().st_size != row["bytes"] or digest(path) != row["sha256"]:
            raise AssertionError({"output_hash_mismatch": rel})
    print(json.dumps({"pass": True, "sources": len(sources), "search_records": len(searches), "screened_titles": len(screen), "data_cards": len(cards), "downloads": len(files), "download_bytes": total, "preserved_inputs_verified": len(preserved), "output_files_hash_verified": len(listed)}, indent=2, allow_nan=False))

if __name__ == "__main__":
    main()
