"""Create a marker-blind manifest for the GSE305971 MARS/BM pilot."""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from openpyxl import load_workbook

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
PAIR_FILE = RESULTS / "external_annotation_by_hospital_biopsy.tsv"
SAMPLE_INDEX = RESULTS / "external_geo_count_file_index.tsv"
ANNOTATION = HERE.parent / "data" / "GSE305971_cells_annotation_for_geo_ver_20250306.xlsx"
OUT = RESULTS / "external_pilot_manifest.tsv"
SUMMARY = RESULTS / "external_pilot_manifest.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    with SAMPLE_INDEX.open("r", encoding="utf-8", newline="") as f:
        index = list(csv.DictReader(f, delimiter="\t"))

    # Method labels identify which GEO subseries is the relevant processed table.
    # MARS (not MARS_v2) is the prespecified assay; this avoids duplicate batch IDs
    # that appear once in each subseries.
    eligible_batches: dict[str, Counter[tuple[str, str, str, str, str, str]]] = defaultdict(Counter)
    summary_by_group: Counter[str] = Counter()
    codes_by_group: dict[str, set[str]] = defaultdict(set)
    workbook = load_workbook(ANNOTATION, read_only=True, data_only=True)
    sheet = workbook.active
    iterator = sheet.iter_rows(values_only=True)
    header = next(iterator)
    columns = {name: i for i, name in enumerate(header)}
    required = ("Amp.Batch.ID", "Hospital.Code", "Biopsy.Sequence", "Disease", "Tissue", "Method", "Populations")
    if any(name not in columns for name in required):
        raise ValueError("Annotation workbook lacks a required selection field")
    for row in iterator:
        batch_id = str(row[columns["Amp.Batch.ID"]]).strip()
        hospital = str(row[columns["Hospital.Code"]]).strip()
        sequence = str(row[columns["Biopsy.Sequence"]]).strip()
        disease = str(row[columns["Disease"]]).strip()
        tissue = str(row[columns["Tissue"]]).strip()
        method = str(row[columns["Method"]]).strip()
        population = str(row[columns["Populations"]]).strip()
        if sequence != "1" or tissue != "BM" or method != "MARS":
            continue
        if disease == "Healthy" and population == "Normal_PC":
            group = "healthy_normal_pc"
        elif disease == "Healthy" and population == "Malignant":
            group = "healthy_ambiguous_malignant_label"
        elif disease in {"NDMM", "RRMM"} and population == "Malignant":
            group = "mm_malignant"
        elif disease in {"NDMM", "RRMM"} and population == "Normal_PC":
            group = "mm_normal_pc_prior_art_reference"
        else:
            # In particular, do not reinterpret Healthy/Malignant cells as
            # either normal controls or disease cases.
            continue
        eligible_batches[batch_id][(group, hospital, sequence, disease, tissue, population)] += 1
        summary_by_group[group] += 1
        codes_by_group[group].add(hospital)
    workbook.close()

    index_by_batch: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in index:
        batch = row["Amp.Batch.ID"]
        if row["GEO_series"] in {"GSE305967", "GSE305968"}:
            index_by_batch[batch].append(row)

    output_rows = []
    unresolved = {}
    for batch_id, labels in sorted(eligible_batches.items()):
        records = index_by_batch.get(batch_id, [])
        # Multiple same-modality GSM records are not resolved by expression values.
        if len(records) != 1:
            unresolved[batch_id] = {"matching_MARS_GEO_records": len(records)}
            continue
        record = records[0]
        for (group, hospital, sequence, disease, tissue, population), n_cells in sorted(labels.items()):
            output_rows.append({
                "analysis_group": group,
                "Hospital.Code": hospital,
                "Biopsy.Sequence": sequence,
                "Disease": disease,
                "Tissue": tissue,
                "Method": "MARS",
                "Populations": population,
                "eligible_annotation_cells_in_batch": str(n_cells),
                "Amp.Batch.ID": batch_id,
                "GEO_series": record["GEO_series"],
                "GSM": record["GSM"],
                "filename": record["filename"],
                "download_url": record["https_download_url"],
                "source_title": record["GSM_title"],
            })

    fields = (
        "analysis_group", "Hospital.Code", "Biopsy.Sequence", "Disease", "Tissue", "Method",
        "Populations", "eligible_annotation_cells_in_batch", "Amp.Batch.ID", "GEO_series", "GSM",
        "filename", "download_url", "source_title",
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        writer.writerows(output_rows)

    unique_files = {(r["GSM"], r["filename"], r["download_url"]) for r in output_rows}
    file_counts = Counter(r["analysis_group"] for r in output_rows)
    summary = {
        "selection_is_marker_blind": True,
        "selection_rule": "Biopsy.Sequence=1; Tissue=BM; Method=MARS; Healthy+Normal_PC controls; NDMM/RRMM+Malignant cases; NDMM/RRMM+Normal_PC cells retained only for same-study paired-expression comparator; Healthy+Malignant cells retained only for uncertain-label sensitivity; putative participant grouping by Hospital.Code.",
        "anomalous_label_rule": "Healthy+Malignant annotations are separately inventoried for a prespecified secondary sensitivity only; primary controls use Healthy+Normal_PC. The ambiguous cells are not called true malignant or normal cells.",
        "putative_codes_by_group": {key: len(value) for key, value in sorted(codes_by_group.items())},
        "annotated_candidate_cells_by_group_before_count_table_join": dict(sorted(summary_by_group.items())),
        "selected_GSM_file_count": len(unique_files),
        "manifest_rows": len(output_rows),
        "manifest_rows_by_group": dict(sorted(file_counts.items())),
        "eligible_batches_without_exactly_one_MARS_GEO_file": unresolved,
        "selected_manifest": str(OUT),
        "selected_manifest_sha256": sha256(OUT),
        "important_limit": "Hospital.Code is used conservatively to keep possible repeat specimens together; the source does not validate it as a unique person identifier. Workbook annotation counts are candidates only until exact cell-ID joins to count matrices are measured.",
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
