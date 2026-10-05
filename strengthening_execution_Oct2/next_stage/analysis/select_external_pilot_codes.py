"""Select a small whole-code external feasibility pilot without marker values."""

from __future__ import annotations

import csv
import hashlib
import json
import random
from collections import Counter
from pathlib import Path


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
FULL = RESULTS / "external_pilot_manifest.tsv"
OUT = RESULTS / "external_pilot_manifest_locked.tsv"
SUMMARY = RESULTS / "external_pilot_manifest_locked.json"
SEED = 20261003
HEALTHY_N = 8
MM_N = 16


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    with FULL.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    healthy_codes = sorted({r["Hospital.Code"] for r in rows if r["analysis_group"] == "healthy_normal_pc"})
    mm_codes = sorted({r["Hospital.Code"] for r in rows if r["analysis_group"] == "mm_malignant"})
    if len(healthy_codes) < HEALTHY_N or len(mm_codes) < MM_N:
        raise ValueError("Not enough eligible putative codes for the locked feasibility pilot")
    rng = random.Random(SEED)
    selected_healthy = sorted(rng.sample(healthy_codes, HEALTHY_N))
    selected_mm = sorted(rng.sample(mm_codes, MM_N))
    selected = set(selected_healthy) | set(selected_mm)
    pilot_rows = [r for r in rows if r["Hospital.Code"] in selected]

    # Retain only the appropriate comparator groups for each selected source key.
    pilot_rows = [
        r for r in pilot_rows
        if (r["Hospital.Code"] in selected_healthy and r["analysis_group"] in {
            "healthy_normal_pc", "healthy_ambiguous_malignant_label"
        }) or (r["Hospital.Code"] in selected_mm and r["analysis_group"] in {
            "mm_malignant", "mm_normal_pc_prior_art_reference"
        })
    ]
    fields = tuple(pilot_rows[0])
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        writer.writerows(pilot_rows)

    unique_gsms = {(r["GSM"], r["filename"], r["download_url"]) for r in pilot_rows}
    cells_by_group = Counter()
    codes_by_group: dict[str, set[str]] = {}
    for row in pilot_rows:
        cells_by_group[row["analysis_group"]] += int(row["eligible_annotation_cells_in_batch"])
        codes_by_group.setdefault(row["analysis_group"], set()).add(row["Hospital.Code"])
    summary = {
        "selection_is_marker_blind": True,
        "random_seed": SEED,
        "selection_algorithm": "Python random.Random(seed).sample(sorted eligible Hospital.Code values, n), after metadata-only eligibility filters; all corresponding MARS/BM/first-biopsy batches are retained.",
        "healthy_normal_pc_eligible_codes": len(healthy_codes),
        "mm_malignant_eligible_codes": len(mm_codes),
        "selected_healthy_codes": selected_healthy,
        "selected_mm_codes": selected_mm,
        "selected_source_row_count": len(pilot_rows),
        "selected_GSM_files": len(unique_gsms),
        "annotation_cells_by_group_before_count_table_join": dict(sorted(cells_by_group.items())),
        "putative_codes_by_group": {key: len(value) for key, value in sorted(codes_by_group.items())},
        "pilot_manifest": str(OUT),
        "pilot_manifest_sha256": sha256(OUT),
        "full_metadata_candidate_manifest": str(FULL),
        "full_metadata_candidate_manifest_sha256": sha256(FULL),
        "scope": "Small feasibility pilot only; not an atlas-wide estimate. Hospital.Code remains an unverified putative-person grouping key.",
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
