"""Verify the completed full selected-cell extract and regenerate policy outputs deterministically."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
NEXT = HERE.parent
DATA = NEXT / "data"
RESULTS = HERE / "results"
sys.path.insert(0, str(HERE))
import run_external_full_expansion as full  # noqa: E402
import run_external_pilot as pilot  # noqa: E402


CELL_TABLE = RESULTS / "external_full_selected_cells.tsv"
LOCK = RESULTS / "external_full_selection_lock.json"
ACTUAL = RESULTS / "external_full_actual_download_manifest.tsv"
DOWNLOAD_SUMMARY = RESULTS / "external_full_download.json"
FILE_AUDIT = RESULTS / "external_full_file_join_audit.tsv"
UNMATCHED = RESULTS / "external_full_unmatched_selected_annotations.tsv"
SUMMARY = RESULTS / "external_full_result_summary.json"
PROVENANCE = RESULTS / "external_full_deterministic_reanalysis.json"
CELL_LOCK = RESULTS / "external_full_selected_cells_first_parse_lock.json"
POLICY_FILES = (
    "external_full_newcode_thresholds.tsv", "external_full_newcode_by_code.tsv",
    "external_full_uncertain_label_sensitivity.tsv", "external_full_uncertain_label_by_code.tsv",
    "external_full_cv_thresholds.tsv", "external_full_cv_by_code.tsv", "external_full_cv_summary.tsv",
    "external_full_cv_uncertain_sensitivity.tsv", "external_full_paired_contrasts.tsv",
)
ARCHIVE = RESULTS / "inflight_policy_outputs_before_deterministic_regeneration"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        writer.writerows(rows)


def archive_inflight_outputs() -> dict[str, str]:
    if ARCHIVE.exists():
        manifest_path = ARCHIVE / "archive_manifest.json"
        if not manifest_path.exists():
            raise ValueError("Inflight archive already exists without a manifest; inspect before writing")
        recorded_manifest_sha = manifest_path.with_suffix(".json.sha256")
        archived = json.loads(manifest_path.read_text(encoding="utf-8"))["files"]
        for name, expected_sha in archived.items():
            path = ARCHIVE / name
            if not path.exists() or sha256(path) != expected_sha:
                raise ValueError(f"Preserved inflight policy output was changed: {name}")
        if not recorded_manifest_sha.exists():
            # Seal the older first-write archive only after every stored file matches its embedded hash.
            recorded_manifest_sha.write_text(f"{sha256(manifest_path)}  {manifest_path.name}\n", encoding="ascii")
        if recorded_manifest_sha.read_text(encoding="ascii").split()[0] != sha256(manifest_path):
            raise ValueError("Inflight archive manifest sidecar is invalid")
        return archived
    ARCHIVE.mkdir(parents=True)
    files: dict[str, str] = {}
    for name in (*POLICY_FILES, "external_full_result_summary.json"):
        source = RESULTS / name
        if not source.exists():
            raise FileNotFoundError(f"Missing inflight policy result: {name}")
        target = ARCHIVE / name
        shutil.copy2(source, target)
        files[name] = sha256(target)
    manifest_path = ARCHIVE / "archive_manifest.json"
    manifest_path.write_text(
        json.dumps({
            "status": "provisional_inflight_outputs_superseded",
            "reason": "fixed-training policy tables were initially produced before code-keyed reference dictionaries were sorted; preserve those outputs but regenerate all policy result tables deterministically from the verified cell extract",
            "files": files,
        }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest_path.with_suffix(".json.sha256").write_text(f"{sha256(manifest_path)}  {manifest_path.name}\n", encoding="ascii")
    return files


def main() -> None:
    # Verify fixed source/protocol locks and downloaded raw sources before reanalysis.
    for path, key in (
        (full.ANNOTATION, "annotation_sha256"),
        (full.GEO_FILE_INDEX, "geo_file_index_sha256"),
        (full.GEO_SAMPLE_INDEX, "geo_sample_index_sha256"),
        (full.PILOT_MANIFEST, "pilot_manifest_sha256"),
        (full.PILOT_DOWNLOADS, "pilot_download_sha256"),
        (full.FULL_MANIFEST, "full_manifest_sha256"),
        (full.SIZE_LIST, "size_list_sha256"),
        (full.LOCK_JSON, "selection_lock_sha256"),
        (full.PROTOCOL_V2, "v2_sha256"),
        (full.V2_ADDENDUM, "v2_addendum_sha256"),
        (full.PROTOCOL_V3, "v3_sha256"),
    ):
        if sha256(path) != full.EXPECTED[key]:
            raise ValueError(f"Locked input changed: {path.name}")
    for protocol, key in ((full.PROTOCOL_V2, "v2_sha256"), (full.V2_ADDENDUM, "v2_addendum_sha256"), (full.PROTOCOL_V3, "v3_sha256")):
        full.verify_sidecar(protocol, full.EXPECTED[key])

    actual_manifest_sha = sha256(ACTUAL)
    if actual_manifest_sha != "9275d41c22cc9359aa2c21b1b1b72e3fb2bad080534f229f4c5b146ae214f15f":
        raise ValueError("Actual downloaded-file manifest hash differs from the completed full retrieval")
    actual_rows = read_tsv(ACTUAL)
    if len(actual_rows) != 456:
        raise ValueError(f"Expected 456 downloaded count tables, observed {len(actual_rows)}")
    raw_files_verified = 0
    raw_bytes = 0
    for row in actual_rows:
        path = DATA / row["filename"]
        size = path.stat().st_size
        if size != int(row["bytes"]) or sha256(path) != row["sha256"]:
            raise ValueError(f"Downloaded source does not match its manifest: {row['filename']}")
        raw_files_verified += 1
        raw_bytes += size
    download_summary = json.loads(DOWNLOAD_SUMMARY.read_text(encoding="utf-8"))
    data_bytes = sum(p.stat().st_size for p in DATA.rglob("*") if p.is_file())
    if data_bytes != int(download_summary["actual_next_stage_data_bytes"]) or data_bytes > full.MAX_TOTAL:
        raise ValueError("Data directory bytes differ from the locked download summary or 500 MiB cap")

    # Validate the original row-level cell extract against annotation keys and file audit.
    if not CELL_TABLE.exists():
        raise FileNotFoundError(CELL_TABLE)
    cell_table_sha = sha256(CELL_TABLE)
    if not CELL_LOCK.exists() or not CELL_LOCK.with_suffix(CELL_LOCK.suffix + ".sha256").exists():
        raise FileNotFoundError("First-parse selected-cell table lock and SHA sidecar must be sealed before reanalysis")
    cell_lock_sha = sha256(CELL_LOCK)
    sidecar_sha = CELL_LOCK.with_suffix(CELL_LOCK.suffix + ".sha256").read_text(encoding="ascii").split()[0]
    cell_lock = json.loads(CELL_LOCK.read_text(encoding="utf-8"))
    if sidecar_sha != cell_lock_sha or cell_table_sha != cell_lock["selected_cell_table_sha256"]:
        raise ValueError("Selected-cell table differs from its independently sealed first-parse lock")
    if (sha256(ACTUAL) != cell_lock["actual_download_manifest_sha256"]
            or sha256(FILE_AUDIT) != cell_lock["file_join_audit_sha256"]
            or sha256(UNMATCHED) != cell_lock["unmatched_annotation_file_sha256"]):
        raise ValueError("First-parse input/audit files differ from the immutable cell-extract lock")
    cells = read_tsv(CELL_TABLE)
    if len(cells) != 62811 or len({r["cell_id"] for r in cells}) != len(cells):
        raise ValueError("Selected-cell table has an unexpected row count or duplicate cell ID")
    lock = json.loads(LOCK.read_text(encoding="utf-8"))
    universe = lock["code_universe"]
    healthy_all = set(universe["healthy_normal_pc_all"])
    mm_all = set(universe["mm_malignant_all"])
    train_healthy = set(universe["pilot_healthy_training_codes_seen_in_pilot"])
    new_healthy = set(universe["new_healthy_evaluation_codes_not_seen_in_pilot"])
    pilot_mm = set(universe["pilot_mm_codes_excluded_from_confirmatory_evaluation"])
    new_mm = set(universe["new_mm_evaluation_codes_not_seen_in_pilot"])
    secondary_mm = set(universe["secondary_full_code_cv_mm_codes_excluding_preexposed"])
    selected_codes = healthy_all | mm_all
    actual_by_filename = {r["filename"]: r for r in actual_rows}
    # The frozen 717-row workbook manifest is one row per eligible annotation stratum;
    # the actual 456-row file manifest is one row per unique count-table batch.
    batch_records = actual_rows
    batch_ids = {r["Amp.Batch.ID"] for r in batch_records}
    if len(batch_records) != 456 or len(batch_ids) != 456:
        raise ValueError("Full batch manifest does not have 456 unique batch partitions")
    annotations, annotated_counts = pilot.read_selected_annotations(selected_codes, batch_ids)
    batch_to_source = {row["Amp.Batch.ID"]: row for row in batch_records}
    cell_by_id = {row["cell_id"]: row for row in cells}
    if set(cell_by_id) - set(annotations):
        raise ValueError("A selected cell ID is absent from the source annotation workbook")
    compare_fields = ("Amp.Batch.ID", "Hospital.Code", "Biopsy.Sequence", "Disease", "Tissue", "Method", "Populations", "Seq.Batch.ID", "analysis_group")
    for cell_id, cell in cell_by_id.items():
        annotation = annotations[cell_id]
        for field in compare_fields:
            if str(cell[field]) != str(annotation[field]):
                raise ValueError(f"Selected cell annotation mismatch for {cell_id}: {field}")
        source_record = batch_to_source[str(cell["Amp.Batch.ID"])]
        if str(cell["GSM"]) != str(source_record["GSM"]):
            raise ValueError(f"Selected cell GSM does not match the locked batch manifest: {cell_id}")
        total = int(cell["total_mRNA_molecule_counts"])
        if total <= 0:
            raise ValueError(f"Non-positive molecule total in selected cell {cell_id}")
        for marker in pilot.MARKERS:
            raw = int(cell[f"{marker}_molecules"])
            norm = float(cell[f"{marker}_per_10k"])
            if raw < 0 or raw > total or not math.isclose(norm, raw / total * 10000.0, rel_tol=0.0, abs_tol=1e-10):
                raise ValueError(f"Raw/normalized molecule count check failed in {cell_id}, {marker}")
    unmatched_rows = read_tsv(UNMATCHED)
    expected_unmatched_ids = set(annotations) - set(cell_by_id)
    recorded_unmatched_ids = {r["cell_id"] for r in unmatched_rows}
    if recorded_unmatched_ids != expected_unmatched_ids or len(unmatched_rows) != len(recorded_unmatched_ids):
        raise ValueError("Unmatched selected annotation IDs do not equal the annotation minus exact-join set")
    if len(expected_unmatched_ids) != 2697:
        raise ValueError(f"Unexpected count of unmatched selected annotations: {len(expected_unmatched_ids)}")
    audit = read_tsv(FILE_AUDIT)
    if len(audit) != 456 or sum(int(r["exact_joined_selected_same_batch_cells"]) for r in audit) != len(cells):
        raise ValueError("Per-file join audit does not reconcile with the selected-cell table")
    if sum(int(r["unmatched_selected_annotation_cell_count"]) for r in audit) != len(unmatched_rows):
        raise ValueError("Per-file unmatched counts do not reconcile with annotation IDs")
    audit_by_file = {r["file"]: r for r in audit}
    for filename, table in actual_by_filename.items():
        file_row = audit_by_file.get(filename)
        if not file_row or file_row["table_sha256"] != table["sha256"]:
            raise ValueError(f"Selected-file join audit lacks the actual source hash: {filename}")
        if int(file_row["cross_batch_annotation_matches_refused"]) != 0:
            raise ValueError("Cross-batch matches were not fully refused")

    cells_numeric: list[dict[str, object]] = []
    by_group_code: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for raw_cell in cells:
        cell: dict[str, object] = dict(raw_cell)
        cell["total_mRNA_molecule_counts"] = int(raw_cell["total_mRNA_molecule_counts"])
        for marker in pilot.MARKERS:
            cell[f"{marker}_molecules"] = int(raw_cell[f"{marker}_molecules"])
            cell[f"{marker}_per_10k"] = float(raw_cell[f"{marker}_per_10k"])
        cells_numeric.append(cell)
        by_group_code[(str(cell["analysis_group"]), str(cell["Hospital.Code"]))].append(cell)

    # Keep the first in-flight outputs intact and clearly labeled as provisional.
    archived_hashes = archive_inflight_outputs()

    primary_healthy_nonempty = {c for c in new_healthy if by_group_code.get(("healthy_normal_pc", c), [])}
    primary_mm_nonempty = {c for c in new_mm if by_group_code.get(("mm_malignant", c), [])}
    fixed = full.evaluate_fixed_training(train_healthy, new_healthy, new_mm, by_group_code)
    uncertain = full.evaluate_fixed_training(train_healthy, new_healthy, new_mm, by_group_code, include_uncertain_reference=True)
    full_healthy_nonempty = {c for c in healthy_all if by_group_code.get(("healthy_normal_pc", c), [])}
    full_mm_nonempty = {c for c in secondary_mm if by_group_code.get(("mm_malignant", c), [])}
    secondary_cv = pilot.evaluate_policy_set("exposure_affected_full_17H_105MM_secondary", full_healthy_nonempty, full_mm_nonempty, by_group_code)
    sensitivity_healthy = {c for c in full_healthy_nonempty if len(by_group_code[("healthy_normal_pc", c)]) >= pilot.MIN_CELLS}
    sensitivity_mm = {c for c in full_mm_nonempty if len(by_group_code[("mm_malignant", c)]) >= pilot.MIN_CELLS}
    secondary_50 = pilot.evaluate_policy_set("exposure_affected_min50_sensitivity", sensitivity_healthy, sensitivity_mm, by_group_code)
    full_cv_fold_rows = list(secondary_cv["cv_fold_rows"]) + list(secondary_50["cv_fold_rows"])
    full_cv_code_rows = list(secondary_cv["cv_code_rows"]) + list(secondary_50["cv_code_rows"])
    full_cv_summary_rows = list(secondary_cv["cv_summary_rows"]) + list(secondary_50["cv_summary_rows"])
    for row in full_cv_code_rows:
        code = str(row["Hospital.Code"])
        row["pilot_marker_exposed_code"] = code in (full.PILOT_EXPOSED_CODES | full.PILOT_EXPOSED_MM_CODES)
        row["prior_single_file_exposed_code"] = code in pilot.PREEXPOSED_CODES
    for row in full_cv_fold_rows:
        row["pilot_training_reference_exposed"] = True

    write_tsv(RESULTS / "external_full_newcode_thresholds.tsv",
        ("analysis_set", "fit_reference", "marker", "budget", "policy", "threshold_per_10k", "threshold_state",
         "training_healthy_codes", "training_healthy_cells", "training_policy_constraint_fraction_at_or_above_threshold",
         "heldout_group", "heldout_codes", "heldout_cells", "heldout_positive_cells", "heldout_cell_weighted_positive_fraction",
         "heldout_equal_code_positive_fraction", "minimum_code_positive_fraction", "median_code_positive_fraction",
         "maximum_code_positive_fraction", "minimum_one_cell_fraction_resolution", "healthy_codes_above_normal_budget",
         "mm_codes_with_coverage_above_illustrative_fraction"), fixed["summary_rows"])
    write_tsv(RESULTS / "external_full_newcode_by_code.tsv",
        ("analysis_set", "fit_reference", "marker", "budget", "policy", "heldout_group", "Hospital.Code", "n_cells",
         "one_cell_fraction_resolution", "positive_cells", "positive_fraction", "threshold_per_10k", "threshold_state",
         "previously_marker_exposed_code", "prior_single_file_exposed_code"), fixed["code_rows"])
    write_tsv(RESULTS / "external_full_uncertain_label_sensitivity.tsv",
        ("analysis_set", "fit_reference", "marker", "budget", "policy", "threshold_per_10k", "threshold_state",
         "training_healthy_codes", "training_healthy_cells", "training_policy_constraint_fraction_at_or_above_threshold",
         "heldout_group", "heldout_codes", "heldout_cells", "heldout_positive_cells", "heldout_cell_weighted_positive_fraction",
         "heldout_equal_code_positive_fraction", "minimum_code_positive_fraction", "median_code_positive_fraction",
         "maximum_code_positive_fraction", "minimum_one_cell_fraction_resolution", "healthy_codes_above_normal_budget",
         "mm_codes_with_coverage_above_illustrative_fraction"), uncertain["summary_rows"])
    write_tsv(RESULTS / "external_full_uncertain_label_by_code.tsv",
        ("analysis_set", "fit_reference", "marker", "budget", "policy", "heldout_group", "Hospital.Code", "n_cells",
         "one_cell_fraction_resolution", "positive_cells", "positive_fraction", "threshold_per_10k", "threshold_state",
         "previously_marker_exposed_code", "prior_single_file_exposed_code"), uncertain["code_rows"])
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
        row["pilot_marker_exposed_code"] = row["Hospital.Code"] in full.PILOT_EXPOSED_MM_CODES
    write_tsv(RESULTS / "external_full_paired_contrasts.tsv",
        ("analysis_set", "Hospital.Code", "marker", "n_malignant_cells", "n_within_code_Normal_PC_cells",
         "paired_mean_difference_malignant_minus_Normal_PC", "normal_reference_sd", "standardized_mean_contrast",
         "standardized_contrast_state", "paired_contrast_50_cell_sensitivity", "preexposed_code", "pilot_marker_exposed_code"),
        primary_pairs + secondary_pairs)

    # Add exact missing-cell worst-case bounds at the frozen threshold for every new code.
    flow = {(r["Hospital.Code"], r["analysis_group"]): r for r in read_tsv(RESULTS / "external_full_code_flow.tsv")}
    training_reference_annotated = sum(int(flow[(code, "healthy_normal_pc")]["workbook_eligible_cells"]) for code in train_healthy)
    training_reference_joined = sum(len(by_group_code.get(("healthy_normal_pc", code), [])) for code in train_healthy)
    training_reference_missing = training_reference_annotated - training_reference_joined
    if training_reference_missing < 0:
        raise ValueError("Training reference joined cell count exceeds annotation count")
    by_code_output = read_tsv(RESULTS / "external_full_newcode_by_code.tsv")
    missing_bound_rows = []
    for row in by_code_output:
        group = {
            "heldout_new_healthy_normal_pc": "healthy_normal_pc",
            "heldout_new_mm_malignant_label": "mm_malignant",
            "heldout_new_uncertain_healthy_malignant_label": "healthy_ambiguous_malignant_label",
        }.get(row["heldout_group"])
        if group is None:
            raise ValueError(f"Unknown fixed-training held-out group {row['heldout_group']}")
        flow_row = flow[(row["Hospital.Code"], group)]
        annotated_n = int(flow_row["workbook_eligible_cells"])
        matched_n = int(row["n_cells"])
        positive_n = int(row["positive_cells"])
        missing_n = annotated_n - matched_n
        if missing_n < 0:
            raise ValueError("Joined cell count exceeds source-eligible annotation count")
        lower = positive_n / annotated_n if annotated_n else "NR"
        upper = (positive_n + missing_n) / annotated_n if annotated_n else "NR"
        missing_bound_rows.append({
            "analysis_set": row["analysis_set"], "marker": row["marker"], "budget": row["budget"], "policy": row["policy"],
            "heldout_group": row["heldout_group"], "Hospital.Code": row["Hospital.Code"],
            "annotated_cells": annotated_n, "matched_cells": matched_n, "unmatched_annotation_cells": missing_n,
            "positive_matched_cells": positive_n, "fraction_lower_all_missing_negative": lower,
            "fraction_upper_all_missing_positive": upper, "threshold_per_10k": row["threshold_per_10k"],
            "threshold_state": row["threshold_state"],
            "interpretation": "worst-case bound among source-eligible annotated cells at the fixed threshold; not an interval estimate or imputed expression",
        })
    write_tsv(RESULTS / "external_full_newcode_missingness_bounds.tsv",
        ("analysis_set", "marker", "budget", "policy", "heldout_group", "Hospital.Code", "annotated_cells", "matched_cells",
         "unmatched_annotation_cells", "positive_matched_cells", "fraction_lower_all_missing_negative",
         "fraction_upper_all_missing_positive", "threshold_per_10k", "threshold_state", "interpretation"), missing_bound_rows)

    summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
    summary["exact_joined_unique_cells"] = len(cells)
    summary["selected_annotation_cell_ids"] = len(annotations)
    summary["selected_annotation_cells_unmatched_in_count_tables"] = len(expected_unmatched_ids)
    summary["source_group_counts"] = dict(sorted(Counter(row["analysis_group"] for row in cells).items()))
    summary["exact_joined_cells_by_source_group"] = summary["source_group_counts"]
    summary["external_full_selected_cells_sha256"] = cell_table_sha
    summary["deterministic_policy_regeneration"] = {
        "status": "complete_from_verified_selected_cell_table",
        "source_cell_table_sha256": cell_table_sha,
        "actual_count_table_manifest_sha256": actual_manifest_sha,
        "actual_count_tables_verified": raw_files_verified,
        "actual_count_table_bytes_verified": raw_bytes,
        "fixed_training_order": "sorted training Hospital.Code keys for reference-array and aggregation reproducibility",
        "first_inflight_policy_outputs": "preserved under inflight_policy_outputs_before_deterministic_regeneration and superseded",
        "cross_batch_match_refusal": "zero accepted; count table and annotation batch keys checked",
        "missing_selected_annotation_ids": len(expected_unmatched_ids),
    }
    summary["outputs"] = list(dict.fromkeys(summary.get("outputs", []) + ["external_full_newcode_missingness_bounds.tsv"]))
    summary["output_sha256"] = {name: sha256(RESULTS / name) for name in summary["outputs"] if (RESULTS / name).is_file()}
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")

    provenance = {
        "status": "deterministic_reanalysis_complete",
        "command": "python strengthening_execution_Oct2/next_stage/analysis/reanalyze_external_full_from_selected_cells.py",
        "source_manifest_sha256": actual_manifest_sha,
        "source_tables_verified": raw_files_verified,
        "source_table_bytes_verified": raw_bytes,
        "data_directory_bytes": data_bytes,
        "selection_lock_sha256": full.EXPECTED["selection_lock_sha256"],
        "v3_protocol_sha256": full.EXPECTED["v3_sha256"],
        "selected_cell_table_sha256": cell_table_sha,
        "selected_cell_table_lock_sha256": cell_lock_sha,
        "selected_cell_rows": len(cells),
        "selected_annotation_rows": len(annotations),
        "unmatched_annotation_rows": len(expected_unmatched_ids),
        "unmatched_by_group": dict(sorted(Counter(r["analysis_group"] for r in unmatched_rows).items())),
        "unmatched_by_code": dict(sorted(Counter(r["Hospital.Code"] for r in unmatched_rows).items())),
        "inflight_policy_files_archived_sha256": archived_hashes,
        "corrected_policy_output_sha256": {name: sha256(RESULTS / name) for name in POLICY_FILES},
        "selected_cell_and_annotation_key_validations": "passed: unique IDs, exact annotation field match, exact source batch, integer counts, marker <= total, normalized counts recomputed, unmatched ID set reconciled",
        "representative_cell_id_namespace_check": "In count table GSM9188838_AB11173.txt.gz the header begins W4313817, W4313818; unmatched selected workbook IDs for AB11173 begin WMC2972705, WMC2972707. No supplied workbook/index crosswalk was verified, so no positional or prefix-based identity mapping was inferred.",
        "missingness_bound_definition": "For each fixed threshold and held-out code, interval is [matched-positive/annotation-eligible, (matched-positive+unmatched-annotation)/annotation-eligible], assigning only worst-case bound endpoints; not a CI, imputation, or part of the frozen primary estimand.",
        "primary_fixed_training_reference_annotation_cells": training_reference_annotated,
        "primary_fixed_training_reference_joined_cells": training_reference_joined,
        "primary_fixed_training_reference_unmatched_cells": training_reference_missing,
    }
    PROVENANCE.write_text(json.dumps(provenance, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "selected_cell_rows": len(cells),
        "annotation_rows": len(annotations),
        "unmatched_annotation_rows": len(expected_unmatched_ids),
        "raw_count_tables_verified": raw_files_verified,
        "source_group_counts": summary["source_group_counts"],
        "policy_output_sha256": provenance["corrected_policy_output_sha256"],
        "missingness_bounds_rows": len(missing_bound_rows),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
