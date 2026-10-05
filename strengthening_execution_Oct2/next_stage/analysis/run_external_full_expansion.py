"""Run the v3 exposure-aware GSE305971 count-table analysis after all locks verify."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
NEXT = HERE.parent
DATA = NEXT / "data"
RESULTS = HERE / "results"
sys.path.insert(0, str(HERE))
import run_external_pilot as pilot  # noqa: E402


ANNOTATION = DATA / "GSE305971_cells_annotation_for_geo_ver_20250306.xlsx"
GEO_FILE_INDEX = RESULTS / "external_geo_count_file_index.tsv"
GEO_SAMPLE_INDEX = RESULTS / "external_geo_sample_index.json"
PILOT_MANIFEST = RESULTS / "external_pilot_manifest_locked.tsv"
PILOT_SUMMARY = RESULTS / "external_pilot_manifest_locked.json"
PILOT_DOWNLOADS = RESULTS / "external_pilot_download_manifest.tsv"
FULL_MANIFEST = RESULTS / "external_full_manifest_locked.tsv"
SIZE_LIST = RESULTS / "external_full_download_manifest_locked.tsv"
LOCK_JSON = RESULTS / "external_full_selection_lock.json"
ACTUAL_DOWNLOADS = RESULTS / "external_full_actual_download_manifest.tsv"
DOWNLOAD_SUMMARY = RESULTS / "external_full_download.json"
PROTOCOL_V2 = HERE / "selected_protocol_v2.md"
V2_ADDENDUM = HERE / "selected_protocol_v2_postpilot_serialization_addendum.md"
PROTOCOL_V3 = HERE / "selected_protocol_v3_full_expansion.md"
MAX_TOTAL = 500 * 1024 * 1024
EXPECTED = {
    "annotation_sha256": "2fc3f9e19ab82097f6a50689de7bf4083aeb81aba442f263ea81afe8b150b94d",
    "geo_file_index_sha256": "2b92d049656bb4c8960c2e4811ea46347687f8b7fab8f24c076ef8a9d031e8ec",
    "geo_sample_index_sha256": "78b35804c7f4a743b969433999934aabafec6af99240ee70f5a355ba827ec943",
    "pilot_manifest_sha256": "43f17abd39a15136aeaaecb183887752c6a76d38d96e1626014abd60953508cb",
    "pilot_download_sha256": "04f24edf9473468838bd377477ea08590399c9cf1820272d5eaddca44fea5864",
    "full_manifest_sha256": "e3e21adee57ab689b6a69086b4e288d0e7176c7e733d3c01d55e718db29dcbe1",
    "size_list_sha256": "c2e8e14beecd2b7aca89caa9e94adfb814c455f92951a4e42ac4a5742fad40ad",
    "selection_lock_sha256": "c67c4e8e2fccf1ca0b3a36e61a5064b49cc446b6e0a411b191fe3f845ddb07c0",
    "v2_sha256": "aa7b16216786b308bd8fb49f375fb481793fc3169bd3e4719bc9966892a7029b",
    "v2_addendum_sha256": "0a1e667e83883cfd644be48ac24aee885fb2b4dd808b6d6da5d677711b20b833",
    "v3_sha256": "62eb79d20342c5fa220b75375c3b847108873830486664583ee2ed4fe31a4461",
}
PILOT_EXPOSED_CODES = set(json.loads(PILOT_SUMMARY.read_text(encoding="utf-8"))["selected_healthy_codes"])
PILOT_EXPOSED_MM_CODES = set(json.loads(PILOT_SUMMARY.read_text(encoding="utf-8"))["selected_mm_codes"])
PRIOR_SINGLE_FILE_EXPOSED = {"tlv53"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        writer.writerows(rows)


def verify_sidecar(path: Path, expected: str) -> None:
    actual = sha256(path)
    sidecar = path.with_suffix(path.suffix + ".sha256").read_text(encoding="ascii").split()[0]
    if actual != expected or sidecar != expected:
        raise ValueError(f"Protocol SHA/sidecar mismatch for {path.name}: {actual} / {sidecar}")


def evaluate_fixed_training(
    training_codes: set[str],
    heldout_healthy_codes: set[str],
    heldout_mm_codes: set[str],
    by_group_code: dict[tuple[str, str], list[dict[str, object]]],
    *,
    include_uncertain_reference: bool = False,
) -> dict[str, list[dict[str, object]]]:
    """Fit on fixed training codes, then summarize all new held-out codes once."""
    if not training_codes:
        raise ValueError("Fixed training set is empty")
    if training_codes & (heldout_healthy_codes | heldout_mm_codes):
        raise ValueError("Training and held-out code sets must be disjoint")
    if heldout_healthy_codes & heldout_mm_codes:
        raise ValueError("Healthy and MM held-out code sets must be disjoint")
    fit_reference = "healthy_normal_pc_plus_uncertain_label_cells" if include_uncertain_reference else "healthy_normal_pc"
    code_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    for marker in pilot.MARKERS:
        score_cache = {
            key: np.asarray([float(row[f"{marker}_per_10k"]) for row in values], dtype=np.float64)
            for key, values in by_group_code.items()
        }
        for budget in (0.10, 0.05, 0.20, 0.50):
            training_normal = {
                code: score_cache[("healthy_normal_pc", code)]
                for code in sorted(training_codes)
            }
            empty_training = [code for code, values in training_normal.items() if values.size == 0]
            if empty_training:
                raise ValueError(f"Fixed training codes lack exact-joined Healthy/Normal_PC cells: {empty_training}")
            if include_uncertain_reference:
                training = {
                    code: np.concatenate((
                        values,
                        score_cache.get(("healthy_ambiguous_malignant_label", code), np.empty(0, dtype=np.float64)),
                    ))
                    for code, values in training_normal.items()
                }
            else:
                training = training_normal
            training_max = max(float(np.max(values)) for values in training.values())
            training_cells = sum(int(values.size) for values in training.values())

            for policy in ("pooled", "equal_code", "donor_envelope"):
                threshold, training_constraint = pilot.fit_thresholds(policy, training, budget)
                threshold_value = pilot.format_threshold(threshold)
                threshold_state = pilot.classify_threshold(threshold, training_max)
                groups = (
                    ("healthy_normal_pc", heldout_healthy_codes, "heldout_new_healthy_normal_pc"),
                    ("mm_malignant", heldout_mm_codes, "heldout_new_mm_malignant_label"),
                    ("healthy_ambiguous_malignant_label", heldout_healthy_codes, "heldout_new_uncertain_healthy_malignant_label"),
                )
                for source_group, codes, heldout_group in groups:
                    per_code_rates: list[float] = []
                    positive_cells = 0
                    total_cells = 0
                    for code in sorted(codes):
                        values = score_cache.get((source_group, code), np.empty(0, dtype=np.float64))
                        if values.size == 0:
                            continue
                        positive = pilot.threshold_positive(values, threshold)
                        rate = positive / int(values.size)
                        per_code_rates.append(rate)
                        total_cells += int(values.size)
                        positive_cells += positive
                        code_rows.append({
                            "analysis_set": "fixed_pilot_reference_new_codes",
                            "fit_reference": fit_reference,
                            "marker": marker,
                            "budget": budget,
                            "policy": policy,
                            "heldout_group": heldout_group,
                            "Hospital.Code": code,
                            "n_cells": int(values.size),
                            "one_cell_fraction_resolution": 1.0 / int(values.size),
                            "positive_cells": positive,
                            "positive_fraction": rate,
                            "threshold_per_10k": threshold_value,
                            "threshold_state": threshold_state,
                            "previously_marker_exposed_code": False,
                            "prior_single_file_exposed_code": False,
                        })
                    summary_rows.append({
                        "analysis_set": "fixed_pilot_reference_new_codes",
                        "fit_reference": fit_reference,
                        "marker": marker,
                        "budget": budget,
                        "policy": policy,
                        "threshold_per_10k": threshold_value,
                        "threshold_state": threshold_state,
                        "training_healthy_codes": len(training),
                        "training_healthy_cells": training_cells,
                        "training_policy_constraint_fraction_at_or_above_threshold": training_constraint,
                        "heldout_group": heldout_group,
                        "heldout_codes": len(per_code_rates),
                        "heldout_cells": total_cells,
                        "heldout_positive_cells": positive_cells,
                        "heldout_cell_weighted_positive_fraction": positive_cells / total_cells if total_cells else "NR",
                        "heldout_equal_code_positive_fraction": float(np.mean(per_code_rates)) if per_code_rates else "NR",
                        "minimum_code_positive_fraction": min(per_code_rates) if per_code_rates else "NR",
                        "median_code_positive_fraction": float(np.median(per_code_rates)) if per_code_rates else "NR",
                        "maximum_code_positive_fraction": max(per_code_rates) if per_code_rates else "NR",
                        "minimum_one_cell_fraction_resolution": min((1.0 / int(score_cache[(source_group, code)].size)
                            for code in codes if score_cache.get((source_group, code), np.empty(0)).size), default="NR"),
                        "healthy_codes_above_normal_budget": sum(rate > budget for rate in per_code_rates) if source_group == "healthy_normal_pc" else "NR",
                        "mm_codes_with_coverage_above_illustrative_fraction": sum(rate > budget for rate in per_code_rates) if source_group == "mm_malignant" else "NR",
                    })
    return {"summary_rows": summary_rows, "code_rows": code_rows}


def main() -> None:
    verify_sidecar(PROTOCOL_V2, EXPECTED["v2_sha256"])
    verify_sidecar(V2_ADDENDUM, EXPECTED["v2_addendum_sha256"])
    verify_sidecar(PROTOCOL_V3, EXPECTED["v3_sha256"])
    expected_paths = {
        ANNOTATION: EXPECTED["annotation_sha256"],
        GEO_FILE_INDEX: EXPECTED["geo_file_index_sha256"],
        GEO_SAMPLE_INDEX: EXPECTED["geo_sample_index_sha256"],
        PILOT_MANIFEST: EXPECTED["pilot_manifest_sha256"],
        PILOT_DOWNLOADS: EXPECTED["pilot_download_sha256"],
        FULL_MANIFEST: EXPECTED["full_manifest_sha256"],
        SIZE_LIST: EXPECTED["size_list_sha256"],
        LOCK_JSON: EXPECTED["selection_lock_sha256"],
    }
    for path, expected in expected_paths.items():
        if sha256(path) != expected:
            raise ValueError(f"Locked input hash mismatch: {path.name}")

    lock = json.loads(LOCK_JSON.read_text(encoding="utf-8"))
    download_summary = json.loads(DOWNLOAD_SUMMARY.read_text(encoding="utf-8"))
    if download_summary["protocol_sha256"] != EXPECTED["v3_sha256"]:
        raise ValueError("Downloaded data summary does not point to the frozen v3 protocol")
    if not download_summary["cap_pass"] or download_summary["unique_processed_tables_verified"] != 456:
        raise ValueError("Downloaded data are incomplete or exceed the frozen data cap")
    actual_download_sha = sha256(ACTUAL_DOWNLOADS)
    if actual_download_sha != download_summary["actual_table_manifest_sha256"]:
        raise ValueError("Actual table manifest hash differs from retrieval summary")
    actual_rows = read_tsv(ACTUAL_DOWNLOADS)
    actual_by_filename = {row["filename"]: row for row in actual_rows}
    if len(actual_by_filename) != 456:
        raise ValueError("Actual download manifest does not contain 456 unique tables")
    source_by_filename = {row["filename"]: row for row in read_tsv(FULL_MANIFEST)}
    if set(actual_by_filename) != set(source_by_filename):
        raise ValueError("Retrieved and marker-blind metadata manifests list different filenames")
    data_bytes = sum(path.stat().st_size for path in DATA.iterdir() if path.is_file())
    if data_bytes != int(download_summary["actual_next_stage_data_bytes"]) or data_bytes > MAX_TOTAL:
        raise ValueError("Downloaded data folder byte count differs from locked run summary or cap")

    for row in actual_rows:
        path = DATA / row["filename"]
        if not path.is_file() or path.stat().st_size != int(row["expected_bytes"]) or sha256(path) != row["sha256"]:
            raise ValueError(f"Downloaded count table failed size/hash validation: {row['filename']}")
    pilot_download_rows = read_tsv(PILOT_DOWNLOADS)
    pilot_by_filename = {row["filename"]: row for row in pilot_download_rows}
    for filename, row in pilot_by_filename.items():
        if filename not in actual_by_filename or actual_by_filename[filename]["sha256"] != row["sha256"]:
            raise ValueError(f"Pilot table hash was not preserved in full run: {filename}")

    universe = lock["code_universe"]
    healthy_all = set(universe["healthy_normal_pc_all"])
    mm_all = set(universe["mm_malignant_all"])
    train_healthy = set(universe["pilot_healthy_training_codes_seen_in_pilot"])
    new_healthy = set(universe["new_healthy_evaluation_codes_not_seen_in_pilot"])
    pilot_mm = set(universe["pilot_mm_codes_excluded_from_confirmatory_evaluation"])
    new_mm = set(universe["new_mm_evaluation_codes_not_seen_in_pilot"])
    secondary_mm = set(universe["secondary_full_code_cv_mm_codes_excluding_preexposed"])
    if train_healthy | new_healthy != healthy_all or train_healthy & new_healthy:
        raise ValueError("Healthy training/held-out code split differs from the marker-blind lock")
    if pilot_mm | new_mm != mm_all or pilot_mm & new_mm:
        raise ValueError("MM exposure/held-out code split differs from the marker-blind lock")
    if secondary_mm != mm_all - pilot.PREEXPOSED_CODES:
        raise ValueError("Secondary MM analysis does not apply the locked pre-exposure exclusion")

    with PILOT_SUMMARY.open("r", encoding="utf-8") as handle:
        pilot_summary = json.load(handle)
    selected_codes = healthy_all | mm_all
    batch_records: dict[str, dict[str, str]] = {}
    for row in actual_rows:
        if row["Amp.Batch.ID"] in batch_records:
            raise ValueError(f"Duplicate processed table for batch {row['Amp.Batch.ID']}")
        batch_records[row["Amp.Batch.ID"]] = row
    batch_ids = set(batch_records)
    annotations, annotated_counts = pilot.read_selected_annotations(selected_codes, batch_ids)
    annotations_by_batch: dict[str, set[str]] = defaultdict(set)
    for cell_id, annotation in annotations.items():
        annotations_by_batch[annotation["Amp.Batch.ID"]].add(cell_id)

    cells: list[dict[str, object]] = []
    file_audit: list[dict[str, object]] = []
    unmatched_annotation_rows: list[dict[str, object]] = []
    batch_mismatch_rows: list[dict[str, object]] = []
    seen_cell_ids: set[str] = set()
    for file_index, batch_id in enumerate(sorted(batch_records), start=1):
        record = {
            **batch_records[batch_id],
            "GEO_series": source_by_filename[batch_records[batch_id]["filename"]]["GEO_series"],
        }
        path = DATA / record["filename"]
        cell_ids, _, totals, marker_values, gene_count = pilot.parse_count_file(path, set(annotations))
        selected_positions = [i for i, cell_id in enumerate(cell_ids) if cell_id in annotations]
        matches = pilot.exact_batch_matches(cell_ids, annotations, batch_id)
        same_batch_positions = [position for position, _, _ in matches]
        batch_mismatches = [i for i in selected_positions if annotations[cell_ids[i]]["Amp.Batch.ID"] != batch_id]
        if batch_mismatches:
            for position in batch_mismatches:
                cell_id = cell_ids[position]
                batch_mismatch_rows.append({
                    "GSM": record["GSM"],
                    "count_file_Amp.Batch.ID": batch_id,
                    "cell_id": cell_id,
                    "annotation_Amp.Batch.ID": annotations[cell_id]["Amp.Batch.ID"],
                    "Hospital.Code": annotations[cell_id]["Hospital.Code"],
                    "analysis_group": annotations[cell_id]["analysis_group"],
                })
            write_tsv(RESULTS / "external_full_batch_key_mismatch.tsv",
                      ("GSM", "count_file_Amp.Batch.ID", "cell_id", "annotation_Amp.Batch.ID", "Hospital.Code", "analysis_group"),
                      batch_mismatch_rows)
            raise ValueError(f"Selected annotation cell IDs occur in the wrong batch: {record['filename']} ({len(batch_mismatches)})")
        cell_ids_in_file = set(cell_ids)
        unmatched_for_batch = sorted(annotations_by_batch[batch_id] - cell_ids_in_file)
        for cell_id in unmatched_for_batch:
            unmatched_annotation_rows.append({
                "Amp.Batch.ID": batch_id,
                "cell_id": cell_id,
                "Hospital.Code": annotations[cell_id]["Hospital.Code"],
                "analysis_group": annotations[cell_id]["analysis_group"],
            })
        for pos in same_batch_positions:
            cell_id = cell_ids[pos]
            if cell_id in seen_cell_ids:
                raise ValueError(f"Selected annotation cell appears in multiple matrices: {cell_id}")
            seen_cell_ids.add(cell_id)
            annotation = annotations[cell_id]
            code = annotation["Hospital.Code"]
            total = int(totals[pos])
            if total <= 0:
                raise ValueError(f"Selected cell has no captured molecules: {cell_id}")
            cell = {
                "cell_id": cell_id,
                "GSM": record["GSM"],
                "GEO_series": record["GEO_series"],
                **annotation,
                "total_mRNA_molecule_counts": total,
                "prior_single_file_exposed_code": code in pilot.PREEXPOSED_CODES,
                "pilot_marker_exposed_code": code in (PILOT_EXPOSED_CODES | PILOT_EXPOSED_MM_CODES),
            }
            for marker in pilot.MARKERS:
                raw = int(marker_values[marker][pos])
                if raw > total:
                    raise ValueError(f"{marker} molecule count exceeds the all-gene total for {cell_id}")
                cell[f"{marker}_molecules"] = raw
                cell[f"{marker}_per_10k"] = (raw / total) * 10000.0
            cells.append(cell)
        file_audit.append({
            "GSM": record["GSM"],
            "Amp.Batch.ID": batch_id,
            "file": record["filename"],
            "matrix_cell_columns": len(cell_ids),
            "gene_rows": gene_count,
            "selected_annotation_cell_ids_for_batch": len(annotations_by_batch[batch_id]),
            "count_columns_with_selected_annotation": len(selected_positions),
            "exact_joined_selected_same_batch_cells": len(same_batch_positions),
            "cross_batch_annotation_matches_refused": 0,
            "unmatched_selected_annotation_cell_ids": json.dumps(unmatched_for_batch, separators=(",", ":")),
            "unmatched_selected_annotation_cell_count": len(unmatched_for_batch),
            "count_columns_without_selected_annotation": len(cell_ids) - len(selected_positions),
            "table_sha256": record["sha256"],
            "pilot_exposed_code_cells_in_batch": sum(
                1 for _, cell_id, _ in matches if annotations[cell_id]["Hospital.Code"] in (PILOT_EXPOSED_CODES | PILOT_EXPOSED_MM_CODES)
            ),
            "prior_single_file_exposed_code_cells_in_batch": sum(
                1 for _, cell_id, _ in matches if annotations[cell_id]["Hospital.Code"] in pilot.PREEXPOSED_CODES
            ),
        })
        if file_index % 40 == 0 or file_index == len(batch_records):
            print(f"parsed {file_index}/{len(batch_records)} tables; exact-joined selected cells={len(seen_cell_ids)}", flush=True)

    cells_fields = (
        "cell_id", "GSM", "GEO_series", "Amp.Batch.ID", "Hospital.Code", "Biopsy.Sequence", "Disease", "Tissue", "Method",
        "Populations", "Seq.Batch.ID", "analysis_group", "total_mRNA_molecule_counts", "TNFRSF17_molecules", "TNFRSF17_per_10k",
        "GPRC5D_molecules", "GPRC5D_per_10k", "SLAMF7_molecules", "SLAMF7_per_10k", "prior_single_file_exposed_code", "pilot_marker_exposed_code",
    )
    write_tsv(RESULTS / "external_full_selected_cells.tsv", cells_fields, cells)
    code_metrics, by_group_code = pilot.metric_rows(cells, selected_codes)

    code_flow: list[dict[str, object]] = []
    for code in sorted(selected_codes):
        for group in (
            "healthy_normal_pc", "healthy_ambiguous_malignant_label", "mm_malignant", "mm_normal_pc_prior_art_reference"
        ):
            matched = len(by_group_code.get((group, code), []))
            annotation_n = annotated_counts.get((group, code), 0)
            code_flow.append({
                "Hospital.Code": code,
                "analysis_group": group,
                "workbook_eligible_cells": annotation_n,
                "exact_joined_cells": matched,
                "one_cell_fraction_resolution": 1.0 / matched if matched else "NR",
                "meets_50_cell_sensitivity": matched >= pilot.MIN_CELLS,
                "pilot_marker_exposed_code": code in (PILOT_EXPOSED_CODES | PILOT_EXPOSED_MM_CODES),
                "prior_single_file_exposed_code": code in pilot.PREEXPOSED_CODES,
                "excluded_from_new_code_primary": code in (PILOT_EXPOSED_CODES | PILOT_EXPOSED_MM_CODES),
                "exclusion_or_limit": (
                    "prior single-file marker exposure; excluded from all policy analysis" if code in pilot.PREEXPOSED_CODES else
                    "pilot marker exposure; excluded from new-code evaluation" if code in (PILOT_EXPOSED_CODES | PILOT_EXPOSED_MM_CODES) else
                    "no exact-joined cells in this source group" if matched == 0 else
                    "below 50-cell screen; retained in primary" if 0 < matched < pilot.MIN_CELLS else "none"
                ),
            })
    for row in code_metrics:
        code = str(row["Hospital.Code"])
        row["pilot_marker_exposed_code"] = code in (PILOT_EXPOSED_CODES | PILOT_EXPOSED_MM_CODES)
        row["prior_single_file_exposed_code"] = code in pilot.PREEXPOSED_CODES

    write_tsv(RESULTS / "external_full_file_join_audit.tsv",
              ("GSM", "Amp.Batch.ID", "file", "matrix_cell_columns", "gene_rows", "selected_annotation_cell_ids_for_batch",
               "count_columns_with_selected_annotation", "exact_joined_selected_same_batch_cells", "cross_batch_annotation_matches_refused",
               "unmatched_selected_annotation_cell_ids", "unmatched_selected_annotation_cell_count", "count_columns_without_selected_annotation", "table_sha256", "pilot_exposed_code_cells_in_batch",
               "prior_single_file_exposed_code_cells_in_batch"), file_audit)
    write_tsv(RESULTS / "external_full_unmatched_selected_annotations.tsv",
              ("Amp.Batch.ID", "cell_id", "Hospital.Code", "analysis_group"), unmatched_annotation_rows)
    write_tsv(RESULTS / "external_full_code_flow.tsv",
              ("Hospital.Code", "analysis_group", "workbook_eligible_cells", "exact_joined_cells", "one_cell_fraction_resolution",
               "meets_50_cell_sensitivity", "pilot_marker_exposed_code", "prior_single_file_exposed_code",
               "excluded_from_new_code_primary", "exclusion_or_limit"), code_flow)
    write_tsv(RESULTS / "external_full_code_metrics.tsv",
              ("analysis_group", "Hospital.Code", "marker", "n_exact_joined_cells", "meets_50_cell_sensitivity", "preexposed_code",
               "median_total_mRNA_molecules", "mean_raw_molecules", "raw_detection_fraction", "mean_per_10k", "median_per_10k", "q90_per_10k",
               "pilot_marker_exposed_code", "prior_single_file_exposed_code"), code_metrics)

    primary_training = train_healthy
    primary_healthy_nonempty = {code for code in new_healthy if by_group_code.get(("healthy_normal_pc", code), [])}
    primary_mm_nonempty = {code for code in new_mm if by_group_code.get(("mm_malignant", code), [])}
    primary_uncertain = evaluate_fixed_training(primary_training, new_healthy, new_mm, by_group_code)
    uncertain_ref_sensitivity = evaluate_fixed_training(
        primary_training, new_healthy, new_mm, by_group_code, include_uncertain_reference=True
    )

    full_healthy_nonempty = {code for code in healthy_all if by_group_code.get(("healthy_normal_pc", code), [])}
    full_mm_nonempty = {code for code in secondary_mm if by_group_code.get(("mm_malignant", code), [])}
    secondary_cv = pilot.evaluate_policy_set(
        "exposure_affected_full_17H_105MM_secondary", full_healthy_nonempty, full_mm_nonempty, by_group_code
    )
    sensitivity_healthy = {code for code in full_healthy_nonempty if len(by_group_code[("healthy_normal_pc", code)]) >= pilot.MIN_CELLS}
    sensitivity_mm = {code for code in full_mm_nonempty if len(by_group_code[("mm_malignant", code)]) >= pilot.MIN_CELLS}
    secondary_50 = pilot.evaluate_policy_set(
        "exposure_affected_min50_sensitivity", sensitivity_healthy, sensitivity_mm, by_group_code
    )

    fixed_summary = primary_uncertain["summary_rows"]
    fixed_codes = primary_uncertain["code_rows"]
    augmented_summary = uncertain_ref_sensitivity["summary_rows"]
    augmented_codes = uncertain_ref_sensitivity["code_rows"]
    full_cv_fold_rows = list(secondary_cv["cv_fold_rows"]) + list(secondary_50["cv_fold_rows"])
    full_cv_code_rows = list(secondary_cv["cv_code_rows"]) + list(secondary_50["cv_code_rows"])
    full_cv_summary_rows = list(secondary_cv["cv_summary_rows"]) + list(secondary_50["cv_summary_rows"])
    for row in full_cv_code_rows:
        code = str(row["Hospital.Code"])
        row["pilot_marker_exposed_code"] = code in (PILOT_EXPOSED_CODES | PILOT_EXPOSED_MM_CODES)
        row["prior_single_file_exposed_code"] = code in pilot.PREEXPOSED_CODES
    for row in full_cv_fold_rows:
        row["pilot_training_reference_exposed"] = True

    write_tsv(RESULTS / "external_full_newcode_thresholds.tsv",
              ("analysis_set", "fit_reference", "marker", "budget", "policy", "threshold_per_10k", "threshold_state",
               "training_healthy_codes", "training_healthy_cells", "training_policy_constraint_fraction_at_or_above_threshold",
               "heldout_group", "heldout_codes", "heldout_cells", "heldout_positive_cells", "heldout_cell_weighted_positive_fraction",
               "heldout_equal_code_positive_fraction", "minimum_code_positive_fraction", "median_code_positive_fraction",
               "maximum_code_positive_fraction", "minimum_one_cell_fraction_resolution", "healthy_codes_above_normal_budget",
               "mm_codes_with_coverage_above_illustrative_fraction"), fixed_summary)
    write_tsv(RESULTS / "external_full_newcode_by_code.tsv",
              ("analysis_set", "fit_reference", "marker", "budget", "policy", "heldout_group", "Hospital.Code", "n_cells",
               "one_cell_fraction_resolution", "positive_cells", "positive_fraction", "threshold_per_10k", "threshold_state",
               "previously_marker_exposed_code", "prior_single_file_exposed_code"), fixed_codes)
    write_tsv(RESULTS / "external_full_uncertain_label_sensitivity.tsv",
              ("analysis_set", "fit_reference", "marker", "budget", "policy", "threshold_per_10k", "threshold_state",
               "training_healthy_codes", "training_healthy_cells", "training_policy_constraint_fraction_at_or_above_threshold",
               "heldout_group", "heldout_codes", "heldout_cells", "heldout_positive_cells", "heldout_cell_weighted_positive_fraction",
               "heldout_equal_code_positive_fraction", "minimum_code_positive_fraction", "median_code_positive_fraction",
               "maximum_code_positive_fraction", "minimum_one_cell_fraction_resolution", "healthy_codes_above_normal_budget",
               "mm_codes_with_coverage_above_illustrative_fraction"), augmented_summary)
    write_tsv(RESULTS / "external_full_uncertain_label_by_code.tsv",
              ("analysis_set", "fit_reference", "marker", "budget", "policy", "heldout_group", "Hospital.Code", "n_cells",
               "one_cell_fraction_resolution", "positive_cells", "positive_fraction", "threshold_per_10k", "threshold_state",
               "previously_marker_exposed_code", "prior_single_file_exposed_code"), augmented_codes)
    write_tsv(RESULTS / "external_full_cv_thresholds.tsv",
              ("analysis_set", "marker", "budget", "policy", "fold", "threshold_per_10k", "threshold_state", "training_healthy_codes",
               "training_healthy_cells", "training_policy_constraint_fraction_at_or_above_threshold", "heldout_group", "heldout_codes",
               "heldout_cells", "heldout_positive_cells", "heldout_cell_weighted_positive_fraction", "heldout_equal_code_positive_fraction",
               "heldout_healthy_codes_above_normal_budget", "heldout_mm_codes_with_coverage_above_illustrative_fraction",
               "pilot_training_reference_exposed"), full_cv_fold_rows)
    write_tsv(RESULTS / "external_full_cv_by_code.tsv",
              ("analysis_set", "marker", "budget", "policy", "fold", "holdout_group", "Hospital.Code", "n_cells",
               "one_cell_fraction_resolution", "positive_cells", "positive_fraction", "threshold_per_10k", "threshold_state",
               "preexposed_code", "pilot_marker_exposed_code", "prior_single_file_exposed_code"), full_cv_code_rows)
    write_tsv(RESULTS / "external_full_cv_summary.tsv",
              ("analysis_set", "marker", "budget", "policy", "heldout_group", "heldout_codes", "heldout_cells",
               "cell_weighted_positive_fraction", "equal_code_positive_fraction", "minimum_code_positive_fraction",
               "median_code_positive_fraction", "maximum_code_positive_fraction", "minimum_one_cell_fraction_resolution",
               "maximum_one_cell_fraction_resolution", "healthy_codes_above_normal_budget",
               "mm_codes_with_coverage_above_illustrative_fraction", "precision_note"), full_cv_summary_rows)
    write_tsv(RESULTS / "external_full_cv_uncertain_sensitivity.tsv",
              ("analysis_set", "marker", "budget", "policy", "fold", "sensitivity_group", "n_codes", "n_cells",
               "one_cell_fraction_resolution_min", "positive_cells", "cell_weighted_positive_fraction", "equal_code_positive_fraction",
               "augmented_training_threshold_per_10k", "augmented_training_threshold_state", "augmented_training_policy_constraint_fraction"),
              list(secondary_cv["uncertain_rows"]) + list(secondary_50["uncertain_rows"]))

    primary_pairs = pilot.paired_contrasts(by_group_code, new_mm)
    secondary_pairs = pilot.paired_contrasts(by_group_code, secondary_mm)
    for row in primary_pairs:
        row["analysis_set"] = "new_90_MM_codes_not_seen_in_pilot"
        row["pilot_marker_exposed_code"] = False
    for row in secondary_pairs:
        row["analysis_set"] = "exposure_affected_secondary_full_105_MM_codes"
        row["pilot_marker_exposed_code"] = row["Hospital.Code"] in PILOT_EXPOSED_MM_CODES
    pair_fields = ("analysis_set", "Hospital.Code", "marker", "n_malignant_cells", "n_within_code_Normal_PC_cells",
                   "paired_mean_difference_malignant_minus_Normal_PC", "normal_reference_sd", "standardized_mean_contrast",
                   "standardized_contrast_state", "paired_contrast_50_cell_sensitivity", "preexposed_code", "pilot_marker_exposed_code")
    write_tsv(RESULTS / "external_full_paired_contrasts.tsv", pair_fields, primary_pairs + secondary_pairs)

    counts_by_group = dict(sorted(Counter(str(row["analysis_group"]) for row in cells).items()))
    summary = {
        "protocol_sha256": EXPECTED["v3_sha256"],
        "v2_pilot_protocol_sha256": EXPECTED["v2_sha256"],
        "v2_pilot_serialization_addendum_sha256": EXPECTED["v2_addendum_sha256"],
        "full_selection_lock_sha256": EXPECTED["selection_lock_sha256"],
        "full_manifest_sha256": EXPECTED["full_manifest_sha256"],
        "download_manifest_sha256": EXPECTED["size_list_sha256"],
        "actual_download_manifest_sha256": actual_download_sha,
        "annotation_workbook_sha256": EXPECTED["annotation_sha256"],
        "GEO_file_index_sha256": EXPECTED["geo_file_index_sha256"],
        "selected_processed_tables": len(actual_rows),
        "selected_annotation_cell_ids": len(annotations),
        "exact_joined_unique_cells": len(cells),
        "selected_annotation_cells_unmatched_in_count_tables": len(unmatched_annotation_rows),
        "unmatched_selected_cell_ids_by_batch": [
            {"Amp.Batch.ID": row["Amp.Batch.ID"], "cell_id": row["cell_id"]}
            for row in unmatched_annotation_rows
        ],
        "exact_joined_cells_by_source_group": counts_by_group,
        "workbook_eligible_annotations_by_source_group": dict(sorted(Counter(
            key[0] for key, count in annotated_counts.items() for _ in range(count)
        ).items())),
        "prior_single_file_exposed_tlv53_cells_excluded_from_all_policy_analyses": sum(
            int(row["n_exact_joined_cells"]) for row in code_metrics if row["Hospital.Code"] in pilot.PREEXPOSED_CODES and row["marker"] == "TNFRSF17"
        ),
        "pilot_exposed_codes": {"healthy_training": sorted(PILOT_EXPOSED_CODES), "mm_excluded_from_primary": sorted(PILOT_EXPOSED_MM_CODES)},
        "primary_new_code_sets": {
            "training_healthy_codes": sorted(train_healthy),
            "new_healthy_evaluation_codes": sorted(new_healthy),
            "pilot_mm_codes_excluded_from_new_code_evaluation": sorted(pilot_mm),
            "new_mm_evaluation_codes": sorted(new_mm),
            "nonempty_joined_healthy_codes": len(primary_healthy_nonempty),
            "nonempty_joined_mm_codes": len(primary_mm_nonempty),
            "empty_or_unjoined_healthy_codes": sorted(new_healthy - primary_healthy_nonempty),
            "empty_or_unjoined_mm_codes": sorted(new_mm - primary_mm_nonempty),
        },
        "secondary_exposure_affected_code_cv": {
            "healthy_codes_nonempty": len(full_healthy_nonempty),
            "mm_codes_nonempty_excluding_tlv53": len(full_mm_nonempty),
            "healthy_codes_ge_50": len(sensitivity_healthy),
            "mm_codes_ge_50": len(sensitivity_mm),
            "four_fold_available": secondary_cv["available"],
            "unavailable_reason": secondary_cv["reason"],
            "pilot_exposed_codes_included": True,
        },
        "method": "integer molecule count normalization: 10,000 * marker mRNA molecule counts / all-gene mRNA molecule counts per cell; strict nonnegative integer parser and exact cell ID plus Amp.Batch.ID join",
        "interpretation": "same-study, assay-level held-out-code policy comparison; no clinical endpoint, external replication, independent patient count, safety claim, or novel target discovery",
        "outputs": [
            "external_full_selected_cells.tsv", "external_full_file_join_audit.tsv", "external_full_code_flow.tsv",
            "external_full_unmatched_selected_annotations.tsv",
            "external_full_code_metrics.tsv", "external_full_newcode_thresholds.tsv", "external_full_newcode_by_code.tsv",
            "external_full_uncertain_label_sensitivity.tsv", "external_full_uncertain_label_by_code.tsv",
            "external_full_cv_thresholds.tsv", "external_full_cv_by_code.tsv", "external_full_cv_summary.tsv",
            "external_full_cv_uncertain_sensitivity.tsv", "external_full_paired_contrasts.tsv",
        ],
    }
    (RESULTS / "external_full_result_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "selected_processed_tables": summary["selected_processed_tables"],
        "selected_annotation_cell_ids": summary["selected_annotation_cell_ids"],
        "exact_joined_unique_cells": summary["exact_joined_unique_cells"],
        "source_group_counts": counts_by_group,
        "new_code_sets": summary["primary_new_code_sets"],
        "secondary_cv_available": summary["secondary_exposure_affected_code_cv"]["four_fold_available"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
