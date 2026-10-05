"""Freeze the marker-blind full MARS/BM/first-biopsy metadata and file manifest."""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
from collections import defaultdict
from pathlib import Path


HERE = Path(__file__).resolve().parent
NEXT = HERE.parent
RESULTS = HERE / "results"
DATA = NEXT / "data"
ANNOTATION = DATA / "GSE305971_cells_annotation_for_geo_ver_20250306.xlsx"
CANDIDATE = RESULTS / "external_pilot_manifest.tsv"
GEO_FILE_INDEX = RESULTS / "external_geo_count_file_index.tsv"
GEO_SAMPLE_INDEX = RESULTS / "external_geo_sample_index.json"
PILOT_MANIFEST = RESULTS / "external_pilot_manifest_locked.tsv"
PILOT_SUMMARY = RESULTS / "external_pilot_manifest_locked.json"
PILOT_DOWNLOAD_MANIFEST = RESULTS / "external_pilot_download_manifest.tsv"
ESTIMATE_TSV = RESULTS / "external_full_size_estimate.tsv"
ESTIMATE_JSON = RESULTS / "external_full_size_estimate.json"
FULL_MANIFEST = RESULTS / "external_full_manifest_locked.tsv"
FULL_DOWNLOADS = RESULTS / "external_full_download_manifest_locked.tsv"
LOCK_JSON = RESULTS / "external_full_selection_lock.json"
MAX_TOTAL = 500 * 1024 * 1024
PREEXPOSED_CODES = {"tlv53"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def save_sidecar(path: Path) -> str:
    digest = sha256(path)
    path.with_suffix(path.suffix + ".sha256").write_text(f"{digest}  {path.name}\n", encoding="ascii")
    return digest


def main() -> None:
    candidate_rows = read_tsv(CANDIDATE)
    pilot = json.loads(PILOT_SUMMARY.read_text(encoding="utf-8"))
    size = json.loads(ESTIMATE_JSON.read_text(encoding="utf-8"))
    size_rows = read_tsv(ESTIMATE_TSV)
    if not size["budget_pass"] or size["head_errors"] != 0:
        raise ValueError("Full processed-table set did not pass the bounded size check")
    if size["metadata_manifest_sha256"] != sha256(CANDIDATE):
        raise ValueError("Size estimate was not produced from the exact candidate manifest")
    if len(size_rows) != size["unique_processed_table_files"]:
        raise ValueError("Size estimate file-list row count does not match summary")

    codes_by_group: dict[str, set[str]] = defaultdict(set)
    unique_files: dict[str, dict[str, str]] = {}
    for row in candidate_rows:
        codes_by_group[row["analysis_group"]].add(row["Hospital.Code"])
        previous = unique_files.get(row["filename"])
        if previous is not None and (previous["download_url"], previous["GSM"], previous["Amp.Batch.ID"]) != (
            row["download_url"], row["GSM"], row["Amp.Batch.ID"]
        ):
            raise ValueError(f"Filename has inconsistent GEO keys: {row['filename']}")
        unique_files[row["filename"]] = row

    healthy = sorted(codes_by_group["healthy_normal_pc"])
    mm = sorted(codes_by_group["mm_malignant"])
    if len(healthy) != 17 or len(mm) != 106:
        raise ValueError(f"Unexpected metadata-only code frame: {len(healthy)} Healthy, {len(mm)} MM")
    pilot_healthy = sorted(pilot["selected_healthy_codes"])
    pilot_mm = sorted(pilot["selected_mm_codes"])
    if len(pilot_healthy) != 8 or len(pilot_mm) != 16:
        raise ValueError("Existing marker-blind pilot allocation no longer matches the expected 8/16 codes")
    if not set(pilot_healthy) <= set(healthy) or not set(pilot_mm) <= set(mm):
        raise ValueError("Pilot training/exposure codes are not contained in the locked metadata frame")
    if "tlv53" not in pilot_mm or not PREEXPOSED_CODES <= set(mm):
        raise ValueError("Previously inspected tlv53 code is not represented in the MM frame")

    new_healthy = sorted(set(healthy) - set(pilot_healthy))
    new_mm = sorted(set(mm) - set(pilot_mm))
    secondary_mm = sorted(set(mm) - PREEXPOSED_CODES)
    if (len(new_healthy), len(new_mm), len(secondary_mm)) != (9, 90, 105):
        raise ValueError("Unexpected new-code evaluation or exposure-excluded full-code sensitivity counts")

    shutil.copyfile(CANDIDATE, FULL_MANIFEST)
    manifest_sha = save_sidecar(FULL_MANIFEST)
    with FULL_DOWNLOADS.open("w", encoding="utf-8", newline="") as handle:
        fields = ("filename", "GSM", "Amp.Batch.ID", "download_url", "expected_bytes", "status", "error")
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        writer.writerows(size_rows)
    download_sha = save_sidecar(FULL_DOWNLOADS)

    lock = {
        "selection_is_marker_blind": True,
        "frozen_before_full_count_table_download_and_marker_parsing": True,
        "annotation_workbook_sha256": sha256(ANNOTATION),
        "GEO_count_file_index_sha256": sha256(GEO_FILE_INDEX),
        "GEO_sample_index_sha256": sha256(GEO_SAMPLE_INDEX),
        "candidate_manifest_sha256": sha256(CANDIDATE),
        "locked_metadata_manifest_sha256": manifest_sha,
        "locked_download_manifest_sha256": download_sha,
        "size_estimate_sha256": sha256(ESTIMATE_JSON),
        "pilot_manifest_sha256": sha256(PILOT_MANIFEST),
        "pilot_selection_summary_sha256": sha256(PILOT_SUMMARY),
        "pilot_download_manifest_sha256": sha256(PILOT_DOWNLOAD_MANIFEST),
        "metadata_candidate_rows": len(candidate_rows),
        "unique_processed_table_files": len(unique_files),
        "unique_gsm_accessions": len({row["GSM"] for row in unique_files.values()}),
        "existing_processed_tables": size["files_already_present"],
        "new_processed_tables": size["missing_files_head_requested"],
        "full_download_expected_bytes_all_tables": size["all_selected_file_bytes_estimated_including_present"],
        "existing_data_bytes": size["existing_next_stage_data_bytes"],
        "estimated_final_data_bytes": size["estimated_next_stage_data_bytes_after_missing_tables"],
        "data_budget_bytes": MAX_TOTAL,
        "budget_pass": size["budget_pass"],
        "code_universe": {
            "healthy_normal_pc_all": healthy,
            "mm_malignant_all": mm,
            "pilot_healthy_training_codes_seen_in_pilot": pilot_healthy,
            "new_healthy_evaluation_codes_not_seen_in_pilot": new_healthy,
            "pilot_mm_codes_excluded_from_confirmatory_evaluation": pilot_mm,
            "new_mm_evaluation_codes_not_seen_in_pilot": new_mm,
            "previously_marker_exposed_code": sorted(PREEXPOSED_CODES),
            "secondary_full_code_cv_healthy_codes": healthy,
            "secondary_full_code_cv_mm_codes_excluding_preexposed": secondary_mm,
        },
        "metadata_counts_by_analysis_group": {
            group: len(codes) for group, codes in sorted(codes_by_group.items())
        },
        "file_identity_rule": "Unique source tables are keyed by (GSM, Amp.Batch.ID, filename), not GSM alone; four GSMs legitimately cover two Amp.Batch.ID/file partitions each.",
        "selection_scope": "First Biopsy.Sequence, BM, MARS; source-annotated Healthy/Normal_PC training/evaluation, NDMM or RRMM/Malignant MM evaluation; retain all metadata-eligible cell counts without a minimum-cell adequacy gate.",
        "label_rule": "Healthy/Malignant annotation is an uncertain source-label group, separately counted and used only in the prespecified sensitivity; it is not relabeled as normal or cancer.",
        "putative_person_key": "Hospital.Code, retained as a conservative grouping key; global person uniqueness and non-overlap with prior cohorts are unverified.",
    }
    LOCK_JSON.write_text(json.dumps(lock, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    lock_sha = save_sidecar(LOCK_JSON)
    print(json.dumps({
        "full_manifest_sha256": manifest_sha,
        "download_manifest_sha256": download_sha,
        "selection_lock_sha256": lock_sha,
        "files": len(unique_files),
        "healthy_training_new_evaluation": [len(pilot_healthy), len(new_healthy)],
        "mm_excluded_new_evaluation_secondary": [len(pilot_mm), len(new_mm), len(secondary_mm)],
        "estimated_final_data_bytes": lock["estimated_final_data_bytes"],
        "budget_pass": lock["budget_pass"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
