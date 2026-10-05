"""Seal the first completed selected-cell extraction before deterministic regeneration."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
LOCK = RESULTS / "external_full_selected_cells_first_parse_lock.json"
CELL_TABLE = RESULTS / "external_full_selected_cells.tsv"
ACTUAL_MANIFEST = RESULTS / "external_full_actual_download_manifest.tsv"
SUMMARY = RESULTS / "external_full_result_summary.json"
JOIN_AUDIT = RESULTS / "external_full_file_join_audit.tsv"
UNMATCHED = RESULTS / "external_full_unmatched_selected_annotations.tsv"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def count_rows(path: Path) -> int:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return sum(1 for _ in csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    if LOCK.exists() or LOCK.with_suffix(LOCK.suffix + ".sha256").exists():
        raise FileExistsError("First-parse selected-cell lock already exists; it is immutable")
    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    if summary.get("selected_processed_tables") != 456 or summary.get("exact_joined_unique_cells") != 62811:
        raise ValueError("Completed first-parse summary does not show the expected 456-table/62,811-cell run")
    if summary.get("actual_download_manifest_sha256") != "9275d41c22cc9359aa2c21b1b1b72e3fb2bad080534f229f4c5b146ae214f15f":
        raise ValueError("Completed parse does not point to the locked actual download manifest")
    if count_rows(ACTUAL_MANIFEST) != 456 or count_rows(JOIN_AUDIT) != 456:
        raise ValueError("Actual manifest and per-file join audit must each contain 456 rows")
    if count_rows(CELL_TABLE) != 62811 or count_rows(UNMATCHED) != 2697:
        raise ValueError("Selected-cell or unmatched annotation row count differs from completed parse")
    lock = {
        "state": "immutable_first_parse_cell_table_lock",
        "first_parse_policy_outputs_state": "provisional_until_sorted_training_regeneration",
        "selected_cell_table_sha256": sha256(CELL_TABLE),
        "selected_cell_rows": count_rows(CELL_TABLE),
        "annotation_cell_rows": summary["selected_annotation_cell_ids"],
        "unmatched_annotation_rows": count_rows(UNMATCHED),
        "actual_download_manifest_sha256": sha256(ACTUAL_MANIFEST),
        "file_join_audit_sha256": sha256(JOIN_AUDIT),
        "unmatched_annotation_file_sha256": sha256(UNMATCHED),
        "result_summary_sha256_at_seal": sha256(SUMMARY),
        "full_v3_protocol_sha256": summary["protocol_sha256"],
        "full_selection_lock_sha256": summary["full_selection_lock_sha256"],
        "note": "This lock freezes the first successfully completed parser extract before policy-table deterministic regeneration. It is not a pristine preregistration claim.",
    }
    LOCK.write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    LOCK.with_suffix(LOCK.suffix + ".sha256").write_text(f"{sha256(LOCK)}  {LOCK.name}\n", encoding="ascii")
    print(json.dumps(lock, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
