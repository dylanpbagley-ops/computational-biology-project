"""Estimate the full metadata-selected processed-table download size with bounded HEAD requests."""

from __future__ import annotations

import concurrent.futures
import csv
import hashlib
import json
import urllib.error
import urllib.request
from pathlib import Path


HERE = Path(__file__).resolve().parent
NEXT = HERE.parent
DATA = NEXT / "data"
RESULTS = HERE / "results"
MANIFEST = RESULTS / "external_pilot_manifest.tsv"
OUT = RESULTS / "external_full_size_estimate.tsv"
SUMMARY = RESULTS / "external_full_size_estimate.json"
MAX_TOTAL = 500 * 1024 * 1024


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def head_file(row: dict[str, str]) -> dict[str, object]:
    request = urllib.request.Request(
        row["download_url"],
        method="HEAD",
        headers={"User-Agent": "Codex public data research/1.0"},
    )
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            length = response.headers.get("Content-Length")
            if not length or not length.isdigit():
                raise ValueError("HEAD response lacks an integer Content-Length")
            return {
                "filename": row["filename"],
                "GSM": row["GSM"],
                "Amp.Batch.ID": row["Amp.Batch.ID"],
                "download_url": row["download_url"],
                "expected_bytes": int(length),
                "status": f"HTTP_{response.status}",
                "error": "",
            }
    except Exception as exc:
        return {
            "filename": row["filename"],
            "GSM": row["GSM"],
            "Amp.Batch.ID": row["Amp.Batch.ID"],
            "download_url": row["download_url"],
            "expected_bytes": "NR",
            "status": "ERROR",
            "error": f"{type(exc).__name__}: {exc}",
        }


def main() -> None:
    with MANIFEST.open("r", encoding="utf-8", newline="") as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    by_filename: dict[str, dict[str, str]] = {}
    for row in records:
        prior = by_filename.get(row["filename"])
        if prior and (prior["download_url"], prior["GSM"], prior["Amp.Batch.ID"]) != (
            row["download_url"], row["GSM"], row["Amp.Batch.ID"]
        ):
            raise ValueError(f"Filename maps to inconsistent source keys: {row['filename']}")
        by_filename[row["filename"]] = row

    existing_data_bytes = sum(path.stat().st_size for path in DATA.iterdir() if path.is_file())
    existing: list[dict[str, object]] = []
    missing: list[dict[str, str]] = []
    for filename, row in sorted(by_filename.items()):
        path = DATA / filename
        if path.is_file():
            existing.append({
                "filename": filename,
                "GSM": row["GSM"],
                "Amp.Batch.ID": row["Amp.Batch.ID"],
                "download_url": row["download_url"],
                "expected_bytes": path.stat().st_size,
                "status": "already_present_local_size",
                "error": "",
            })
        else:
            missing.append(row)

    responses: list[dict[str, object]] = list(existing)
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        futures = [executor.submit(head_file, row) for row in missing]
        for future in concurrent.futures.as_completed(futures):
            responses.append(future.result())
    responses.sort(key=lambda row: str(row["filename"]))
    fields = ("filename", "GSM", "Amp.Batch.ID", "download_url", "expected_bytes", "status", "error")
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        writer.writerows(responses)

    errors = [row for row in responses if row["status"] == "ERROR"]
    estimate = sum(int(row["expected_bytes"]) for row in responses if row["expected_bytes"] != "NR")
    data_after_full = existing_data_bytes + sum(
        int(row["expected_bytes"]) for row in responses if row["status"] == "HTTP_200"
    )
    summary = {
        "metadata_manifest_sha256": sha256(MANIFEST),
        "metadata_rows": len(records),
        "unique_processed_table_files": len(by_filename),
        "files_already_present": len(existing),
        "missing_files_head_requested": len(missing),
        "head_successes": sum(row["status"] == "HTTP_200" for row in responses),
        "head_errors": len(errors),
        "all_selected_file_bytes_estimated_including_present": estimate,
        "existing_next_stage_data_bytes": existing_data_bytes,
        "estimated_next_stage_data_bytes_after_missing_tables": data_after_full if not errors else "INCOMPLETE",
        "budget_bytes": MAX_TOTAL,
        "budget_pass": not errors and data_after_full <= MAX_TOTAL,
        "errors": errors,
        "output_tsv": str(OUT),
    }
    SUMMARY.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
