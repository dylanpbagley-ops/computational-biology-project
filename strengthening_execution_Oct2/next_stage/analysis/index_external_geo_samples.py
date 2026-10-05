"""Index sample-level GEO count files and join them to the cell annotations."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from pathlib import Path


HERE = Path(__file__).resolve().parent
NEXT = HERE.parent
DATA = NEXT / "data"
RESULTS = HERE / "results"
SOFT_FILES = tuple(DATA / f"GSE305{n}_family.soft.gz" for n in (967, 968, 969, 970))
BATCH_INVENTORY = RESULTS / "external_annotation_batch_inventory.tsv"
OUT_INDEX = RESULTS / "external_geo_count_file_index.tsv"
OUT_CANDIDATES = RESULTS / "external_candidate_biopsies.tsv"
OUT_SUMMARY = RESULTS / "external_geo_sample_index.json"
BATCH_PATTERN = re.compile(r"_(AB\d+)\.txt\.gz$", re.IGNORECASE)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_soft(path: Path) -> list[dict[str, object]]:
    series = path.name.split("_", 1)[0]
    samples: list[dict[str, object]] = []
    current: dict[str, object] | None = None

    def save() -> None:
        if current is not None:
            samples.append(current.copy())

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith("^SAMPLE = "):
                save()
                current = {
                    "GEO_series": series,
                    "GSM": line.split(" = ", 1)[1].strip(),
                    "title": "NR",
                    "supplementary_files": [],
                }
            elif current is not None and line.startswith("!Sample_title = "):
                current["title"] = line.split(" = ", 1)[1].strip()
            elif current is not None and re.match(r"!Sample_supplementary_file(?:_\d+)? = ", line):
                current["supplementary_files"].append(line.split(" = ", 1)[1].strip())
    save()

    parsed = []
    for sample in samples:
        for url in sample["supplementary_files"]:
            if not url.lower().endswith(".txt.gz"):
                continue
            match = BATCH_PATTERN.search(url)
            if not match:
                continue
            parsed.append({
                "GEO_series": sample["GEO_series"],
                "GSM": sample["GSM"],
                "GSM_title": sample["title"],
                "Amp.Batch.ID": match.group(1).upper(),
                "supplementary_url": url,
                "https_download_url": url.replace("ftp://", "https://", 1),
                "filename": url.rsplit("/", 1)[-1],
            })
    return parsed


def main() -> None:
    for path in SOFT_FILES:
        if not path.is_file():
            raise FileNotFoundError(path)
    sample_files = [sample for path in SOFT_FILES for sample in parse_soft(path)]
    by_batch: dict[str, list[dict[str, object]]] = {}
    for sample in sample_files:
        by_batch.setdefault(str(sample["Amp.Batch.ID"]), []).append(sample)

    with BATCH_INVENTORY.open("r", encoding="utf-8", newline="") as handle:
        batch_rows = list(csv.DictReader(handle, delimiter="\t"))
    if len({row["Amp.Batch.ID"] for row in batch_rows}) != len(batch_rows):
        raise ValueError("Duplicate Amp.Batch.ID in annotation inventory")
    batch_by_id = {row["Amp.Batch.ID"]: row for row in batch_rows}

    fields = (
        "Amp.Batch.ID", "GEO_series", "GSM", "GSM_title", "filename", "supplementary_url",
        "https_download_url", "annotated_cells", "Hospital.Code_counts", "Biopsy.Sequence_counts",
        "Disease_counts", "Tissue_counts", "Method_counts", "Seq.Batch.ID_counts",
        "Populations_counts", "Hospital_Biopsy_Disease_counts",
        "Disease_Tissue_Method_Population_counts", "batch_join_status",
    )
    RESULTS.mkdir(parents=True, exist_ok=True)
    with OUT_INDEX.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        for sample in sorted(sample_files, key=lambda r: (r["GEO_series"], r["GSM"])):
            batch_id = str(sample["Amp.Batch.ID"])
            annotation = batch_by_id.get(batch_id)
            row = {**sample}
            if annotation is None:
                row.update({field: "NR" for field in fields if field not in row})
                row["batch_join_status"] = "no_annotation_batch_id_match"
            else:
                row.update(annotation)
                row["batch_join_status"] = "one_batch_match" if len(by_batch[batch_id]) == 1 else "multiple_GSM_records_for_batch"
            writer.writerow(row)

    candidate_rows = []
    unique_count = CounterLike()
    for batch_id, matches in sorted(by_batch.items()):
        if len(matches) != 1 or batch_id not in batch_by_id:
            continue
        row = batch_by_id[batch_id]
        joint = json.loads(row["Disease_Tissue_Method_Population_counts"])
        for combo, count in joint.items():
            disease, tissue, method, population = combo.split("|", 3)
            if tissue != "BM" or method != "MARS" or population not in ("Normal_PC", "Malignant"):
                continue
            if disease not in ("Healthy", "NDMM", "RRMM"):
                continue
            code_counts = json.loads(row["Hospital.Code_counts"])
            biopsy_counts = json.loads(row["Biopsy.Sequence_counts"])
            if len(code_counts) != 1 or len(biopsy_counts) != 1:
                continue
            code, code_cells = next(iter(code_counts.items()))
            biopsy, biopsy_cells = next(iter(biopsy_counts.items()))
            candidate_rows.append({
                "Hospital.Code": code,
                "Biopsy.Sequence": biopsy,
                "Disease": disease,
                "Tissue": tissue,
                "Method": method,
                "Populations": population,
                "annotated_cells_in_batch_class": count,
                "Amp.Batch.ID": batch_id,
                "GEO_series": matches[0]["GEO_series"],
                "GSM": matches[0]["GSM"],
                "filename": matches[0]["filename"],
                "supplementary_url": matches[0]["supplementary_url"],
            })
            unique_count.add((code, biopsy, disease, population), int(count))

    candidate_fields = (
        "Hospital.Code", "Biopsy.Sequence", "Disease", "Tissue", "Method", "Populations",
        "annotated_cells_in_batch_class", "Amp.Batch.ID", "GEO_series", "GSM", "filename", "supplementary_url",
    )
    with OUT_CANDIDATES.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=candidate_fields, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        writer.writerows(candidate_rows)

    unmatched_annotation = sorted(set(batch_by_id) - set(by_batch))
    unmatched_geo = sorted(set(by_batch) - set(batch_by_id))
    duplicates = {batch: len(records) for batch, records in by_batch.items() if len(records) != 1}
    summary = {
        "soft_inputs": [{"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)} for path in SOFT_FILES],
        "batch_inventory": {"path": str(BATCH_INVENTORY), "sha256": sha256(BATCH_INVENTORY), "batch_ids": len(batch_rows)},
        "GEO_count_file_records": len(sample_files),
        "unique_count_file_batch_ids": len(by_batch),
        "duplicate_GEO_records_by_batch": duplicates,
        "annotation_batch_ids_without_GEO_count_file": len(unmatched_annotation),
        "GEO_batch_ids_without_annotation": len(unmatched_geo),
        "candidate_batch_rows": len(candidate_rows),
        "candidate_counts_by_group": unique_count.compact(),
        "candidate_files": str(OUT_CANDIDATES),
        "full_index": str(OUT_INDEX),
    }
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, sort_keys=True))


class CounterLike:
    def __init__(self) -> None:
        self._data: dict[tuple[str, ...], int] = {}

    def add(self, key: tuple[str, ...], value: int) -> None:
        self._data[key] = self._data.get(key, 0) + value

    def compact(self) -> dict[str, int]:
        return {"|".join(key): value for key, value in sorted(self._data.items())}


if __name__ == "__main__":
    main()
