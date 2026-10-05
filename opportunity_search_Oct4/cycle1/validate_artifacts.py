"""Read-only checks for the cycle 1 literature/data opportunity map."""
from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
STAGE = Path(__file__).resolve().parent
CACHE_PARTS = {".cache", "__pycache__", "cache"}

TSV_COLUMNS = {
    "search_log.tsv": "search_id utc purpose query_or_url results_used next_action_or_stop_reason".split(),
    "papers.tsv": "paper_id type title year doi_or_accession url problem_family reading_level sections_read measured_inputs measured_outputs independent_units key_limitation_or_gap candidate_links access_limit".split(),
    "data_cards.tsv": "data_id paper_ids accession official_url access_status file_examined bytes_if_known biological_system input_measurements outcome_measurements timing independent_unit independent_unit_count grouping_keys join_keys censoring_or_limits controls evaluation_resource measured_vs_derived candidate_links unresolved_access_or_schema".split(),
    "candidates.tsv": "candidate_id family question biological_value relationship_to_existing paper_ids closest_prior_art exact_increment later_work_search_ids novelty_status data_ids data_status independent_unit primary_endpoint comparator evaluation_design minimal_decisive_test success_gate revise_or_stop_gate main_risks rank_or_rejection_reason".split(),
}


def reject_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON value: {value}")


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"), parse_constant=reject_constant)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_manifest(label: str, manifest_rel: str, base_rel: str) -> dict:
    manifest = load_json(ROOT / manifest_rel)
    rows = manifest["files"]
    if manifest["file_count"] != len(rows):
        raise AssertionError(f"{label}: file_count mismatch")
    failures = []
    base = ROOT / base_rel
    for row in rows:
        rel = Path(row["path"])
        if rel.is_absolute() or ".." in rel.parts:
            failures.append({"path": row["path"], "error": "unsafe relative path"})
            continue
        path = base / rel
        if not path.is_file() or path.stat().st_size != row["bytes"] or sha256(path) != row["sha256"]:
            failures.append({"path": row["path"], "error": "missing or hash mismatch"})
    if failures:
        raise AssertionError({"manifest": label, "failures": failures})
    return {"manifest": label, "listed": len(rows), "verified": len(rows), "pass": True}


def split_ids(value: str) -> list[str]:
    if value in {"", "none", "not_reported", "not_inspected"}:
        return []
    return [part for part in value.split("|") if part]


def read_tsv(path: Path, expected_header: list[str]) -> list[dict[str, str]]:
    raw_lines = path.read_text(encoding="utf-8").splitlines()
    if not raw_lines:
        raise AssertionError(f"empty TSV: {path.name}")
    widths = [len(line.split("\t")) for line in raw_lines]
    if set(widths) != {len(expected_header)}:
        raise AssertionError(f"{path.name}: inconsistent/tab-containing fields; widths={sorted(set(widths))}")
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if list(rows[0].keys()) != expected_header:
        raise AssertionError(f"{path.name}: unexpected header")
    if len(rows) != len(raw_lines) - 1:
        raise AssertionError(f"{path.name}: physical/logical row count differs")
    return rows


def main() -> None:
    required = set(TSV_COLUMNS) | {
        "STATUS.md", "opportunity_map.md", "preservation_before.json",
        "data_manifest.json", "output_hashes.json", "validate_artifacts.py",
    }
    actual_stage = {
        path.relative_to(STAGE).as_posix(): path
        for path in STAGE.rglob("*")
        if path.is_file() and not any(part in CACHE_PARTS for part in path.relative_to(STAGE).parts)
    }
    if not required.issubset(actual_stage):
        raise AssertionError({"missing_required": sorted(required - set(actual_stage))})

    json_paths = [path for path in STAGE.rglob("*.json") if not any(part in CACHE_PARTS for part in path.relative_to(STAGE).parts)]
    for path in json_paths:
        load_json(path)

    tables = {name: read_tsv(STAGE / name, header) for name, header in TSV_COLUMNS.items()}
    id_fields = {
        "search_log.tsv": "search_id",
        "papers.tsv": "paper_id",
        "data_cards.tsv": "data_id",
        "candidates.tsv": "candidate_id",
    }
    ids: dict[str, set[str]] = {}
    for name, field in id_fields.items():
        values = [row[field] for row in tables[name]]
        if any(not value for value in values) or len(values) != len(set(values)):
            raise AssertionError(f"{name}: blank or duplicate {field}")
        ids[field] = set(values)
    valid_refs = {
        "papers": ids["paper_id"],
        "reviews_and_papers": ids["paper_id"],
        "candidates": ids["candidate_id"],
        "data": ids["data_id"],
        "searches": ids["search_id"],
    }
    for row in tables["papers.tsv"]:
        missing = set(split_ids(row["candidate_links"])) - valid_refs["candidates"]
        if missing:
            raise AssertionError({"paper_candidate_refs": row["paper_id"], "missing": sorted(missing)})
    for row in tables["data_cards.tsv"]:
        missing = set(split_ids(row["paper_ids"])) - valid_refs["papers"]
        missing |= set(split_ids(row["candidate_links"])) - valid_refs["candidates"]
        if missing:
            raise AssertionError({"data_card_refs": row["data_id"], "missing": sorted(missing)})
    for row in tables["candidates.tsv"]:
        missing = set(split_ids(row["paper_ids"])) - valid_refs["reviews_and_papers"]
        missing |= set(split_ids(row["data_ids"])) - valid_refs["data"]
        missing |= set(split_ids(row["later_work_search_ids"])) - valid_refs["searches"]
        if missing:
            raise AssertionError({"candidate_refs": row["candidate_id"], "missing": sorted(missing)})
    allowed_reading = {"metadata_only", "abstract", "targeted_sections", "full_main_text", "main_plus_relevant_supplement"}
    if any(row["reading_level"] not in allowed_reading for row in tables["papers.tsv"]):
        raise AssertionError("invalid paper reading_level")
    allowed_access = {"identified", "metadata_verified", "file_accessed", "rows_inspected", "restricted", "blocked", "summary_only"}
    if any(row["access_status"] not in allowed_access for row in tables["data_cards.tsv"]):
        raise AssertionError("invalid data access_status")
    allowed_novelty = {"overlap", "unresolved", "plausible_gap"}
    allowed_data = {"inadequate", "partial", "schema_verified"}
    allowed_biology = {"minor", "potentially_meaningful", "well_motivated"}
    allowed_relationship = {"natural_extension", "related_pivot", "broader_biology_pivot"}
    for row in tables["candidates.tsv"]:
        if row["novelty_status"] not in allowed_novelty or row["data_status"] not in allowed_data:
            raise AssertionError({"candidate_labels": row["candidate_id"]})
        if row["relationship_to_existing"] not in allowed_relationship:
            raise AssertionError({"candidate_relationship": row["candidate_id"]})
        biology_label = row["biological_value"].split(":", 1)[0]
        if biology_label not in allowed_biology:
            raise AssertionError({"candidate_biology_label": row["candidate_id"], "label": biology_label})
    # The log mixes prose/search terms with IDs; validate any explicit record-like token.
    all_known = set().union(*valid_refs.values())
    for row in tables["search_log.tsv"]:
        for token in re.findall(r"\b[PRDC]\d{2}\b", row["results_used"]):
            if token not in all_known:
                raise AssertionError({"search_log_ref": row["search_id"], "missing": token})

    before = load_json(STAGE / "preservation_before.json")
    direct_failures = []
    for row in before["files"]:
        path = ROOT / row["path"]
        if not path.is_file() or path.stat().st_size != row["bytes"] or sha256(path) != row["sha256"]:
            direct_failures.append(row["path"])
    if direct_failures:
        raise AssertionError({"preserved_inputs_changed": direct_failures})

    data_manifest = load_json(STAGE / "data_manifest.json")
    if data_manifest["download_count"] != len(data_manifest["files"]):
        raise AssertionError("download_count mismatch")
    if data_manifest["download_count"] != 0:
        for row in data_manifest["files"]:
            path = STAGE / row["path"]
            if not path.is_file() or path.stat().st_size != row["bytes"] or sha256(path) != row["sha256"]:
                raise AssertionError({"download_hash_mismatch": row["path"]})

    output_manifest = load_json(STAGE / "output_hashes.json")
    output_rows = output_manifest["files"]
    listed = {row["path"]: row for row in output_rows}
    if len(listed) != len(output_rows) or output_manifest["file_count"] != len(output_rows):
        raise AssertionError("output manifest count or uniqueness mismatch")
    if "output_hashes.json" in listed:
        raise AssertionError("output manifest must exclude itself")
    expected_files = {name for name in actual_stage if name != "output_hashes.json"}
    if set(listed) != expected_files:
        raise AssertionError({"stage_file_coverage": {
            "unlisted": sorted(expected_files - set(listed)),
            "missing_from_stage": sorted(set(listed) - expected_files),
        }})
    for rel, row in listed.items():
        path = STAGE / rel
        if path.stat().st_size != row["bytes"] or sha256(path) != row["sha256"]:
            raise AssertionError({"stage_hash_mismatch": rel})

    prior = [
        verify_manifest("functional_extension_Oct3", "functional_extension_Oct3/output_hashes.json", "functional_extension_Oct3"),
        verify_manifest("longitudinal_extension_Oct3", "longitudinal_extension_Oct3/delivery/output_hashes.json", "longitudinal_extension_Oct3"),
        verify_manifest("strengthening_execution_Oct2", "strengthening_execution_Oct2/next_stage/delivery/output_hashes.json", "strengthening_execution_Oct2/next_stage"),
    ]
    print(json.dumps({
        "json_files_finite": len(json_paths),
        "tsv_rows": {name: len(rows) for name, rows in tables.items()},
        "tsv_widths": {name: len(header) for name, header in TSV_COLUMNS.items()},
        "stage_files_hash_verified": len(listed),
        "direct_preservation_inputs_unchanged": len(before["files"]),
        "downloaded_files": data_manifest["download_count"],
        "prior_manifests": prior,
        "pass": True,
    }, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
