"""Validate one public GSE305967 molecule-count table against GEO annotations."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

from openpyxl import load_workbook

from metadata_diagnostics import parse_matrix_counts_row


HERE = Path(__file__).resolve().parent
NEXT = HERE.parent
DATA = NEXT / "data"
COUNT_FILE = DATA / "GSM9189168_AB14564.txt.gz"
ANNOTATION_FILE = DATA / "GSE305971_cells_annotation_for_geo_ver_20250306.xlsx"
GENES = ("TNFRSF17", "GPRC5D", "SLAMF7")
OUT_CELLS = HERE / "results" / "gse305967_GSM9189168_cells.tsv"
OUT_SUMMARY = HERE / "results" / "gse305967_GSM9189168_validation.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_count_file(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        cell_ids = next(reader, None)
        if not cell_ids or len(set(cell_ids)) != len(cell_ids):
            raise ValueError("Count file has missing or duplicate cell IDs")
        totals = [0] * len(cell_ids)
        markers = {}
        gene_count = 0
        for line_no, row in enumerate(reader, start=2):
            if len(row) != len(cell_ids) + 1:
                raise ValueError(f"Malformed count row {line_no}")
            gene = row[0]
            values = parse_matrix_counts_row(row[1:], line_no)
            gene_count += 1
            for index, count in enumerate(values):
                totals[index] += int(count)
            if gene in GENES:
                if gene in markers:
                    raise ValueError(f"Duplicate marker row: {gene}")
                markers[gene] = values
    missing = sorted(set(GENES) - set(markers))
    if missing:
        raise ValueError(f"Marker rows missing from count file: {missing}")
    return cell_ids, totals, markers, gene_count


def read_annotations(path: Path, wanted_cell_ids: set[str]) -> dict[str, dict[str, str]]:
    workbook = load_workbook(path, read_only=True, data_only=True)
    sheet = workbook.active
    rows = sheet.iter_rows(values_only=True)
    header = next(rows, None)
    if not header:
        raise ValueError("Annotation workbook has no header")
    columns = {name: index for index, name in enumerate(header)}
    required = ("index", "Amp.Batch.ID", "Hospital.Code", "Biopsy.Sequence", "Populations", "Disease", "Tissue", "Method", "Seq.Batch.ID")
    if any(name not in columns for name in required):
        raise ValueError(f"Required annotation columns missing: {required}")
    selected = {}
    for row_number, row in enumerate(rows, start=2):
        cell_id = row[columns["index"]]
        if cell_id is None or str(cell_id) not in wanted_cell_ids:
            continue
        cell_id = str(cell_id)
        if cell_id in selected:
            raise ValueError(f"Duplicate annotation cell ID: {cell_id}")
        selected[cell_id] = {
            name: "NR" if row[columns[name]] is None else str(row[columns[name]]).strip()
            for name in required[1:]
        }
    return selected


def main() -> None:
    cell_ids, totals, markers, gene_count = read_count_file(COUNT_FILE)
    annotations = read_annotations(ANNOTATION_FILE, set(cell_ids))
    shared = set(cell_ids) & set(annotations)
    missing = sorted(set(cell_ids) - shared)
    annotation_cell_counts = Counter(
        (annotations[cell_id]["Disease"], annotations[cell_id]["Tissue"], annotations[cell_id]["Populations"])
        for cell_id in sorted(shared)
    )

    OUT_CELLS.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CELLS.open("w", encoding="utf-8", newline="") as handle:
        fields = ("cell_id", "matrix_total_mrna_molecules", *(f"{gene}_molecules" for gene in GENES),
                  "Amp.Batch.ID", "Hospital.Code", "Biopsy.Sequence", "Disease", "Tissue", "Populations", "Method", "Seq.Batch.ID")
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        for index, cell_id in enumerate(cell_ids):
            if cell_id not in annotations:
                continue
            row = {"cell_id": cell_id, "matrix_total_mrna_molecules": totals[index]}
            row.update({f"{gene}_molecules": int(markers[gene][index]) for gene in GENES})
            row.update(annotations[cell_id])
            writer.writerow(row)

    marker_summary = {}
    for gene in GENES:
        values = markers[gene]
        marker_summary[gene] = {
            "cells_in_count_file": len(cell_ids),
            "nonzero_cells": int((values > 0).sum()),
            "molecule_sum": int(values.sum()),
        }
    result = {
        "count_file": str(COUNT_FILE),
        "count_file_bytes": COUNT_FILE.stat().st_size,
        "count_file_sha256": sha256(COUNT_FILE),
        "annotation_file": str(ANNOTATION_FILE),
        "annotation_file_sha256": sha256(ANNOTATION_FILE),
        "table_shape": {"gene_rows": gene_count, "cell_columns": len(cell_ids)},
        "verified_raw_units": "nonnegative integer mRNA molecule counts per cell (GEO sample-level description)",
        "annotation_join": {"matched_cells": len(shared), "unmatched_count_cells": len(missing), "phenotype_population_tissue_counts": {"|".join(key): count for key, count in sorted(annotation_cell_counts.items())}},
        "marker_slice": marker_summary,
        "allowed_inference": "single processed sample feasibility slice only; not donor-policy validation",
        "output_cells": str(OUT_CELLS),
    }
    OUT_SUMMARY.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
