"""Download the locked public processed count tables within the 500 MB budget."""

from __future__ import annotations

import csv
import hashlib
import json
import time
import urllib.request
from pathlib import Path


HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
RESULTS = HERE / "results"
MANIFEST = RESULTS / "external_pilot_manifest_locked.tsv"
OUT = RESULTS / "external_pilot_download_manifest.tsv"
SUMMARY = RESULTS / "external_pilot_download.json"
MAX_TOTAL = 500 * 1024 * 1024


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    with MANIFEST.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    files = {row["filename"]: row for row in rows}
    used_bytes = sum(p.stat().st_size for p in DATA.iterdir() if p.is_file())
    if used_bytes > MAX_TOTAL:
        raise RuntimeError(f"Existing public-data downloads already exceed {MAX_TOTAL} bytes")

    results = []
    for index, (filename, row) in enumerate(sorted(files.items()), start=1):
        destination = DATA / filename
        if destination.exists():
            size = destination.stat().st_size
            status = "already_present"
        else:
            request = urllib.request.Request(
                row["download_url"], headers={"User-Agent": "Codex public data research/1.0"}
            )
            error = None
            for attempt in range(3):
                try:
                    with urllib.request.urlopen(request, timeout=60) as response:
                        length = response.headers.get("Content-Length")
                        expected = int(length) if length and length.isdigit() else None
                        if expected is not None and used_bytes + expected > MAX_TOTAL:
                            raise RuntimeError("500 MB cumulative public-data budget would be exceeded")
                        temp = destination.with_suffix(destination.suffix + ".part")
                        size = 0
                        with temp.open("wb") as out:
                            while True:
                                block = response.read(1 << 20)
                                if not block:
                                    break
                                size += len(block)
                                if used_bytes + size > MAX_TOTAL:
                                    raise RuntimeError("500 MB cumulative public-data budget would be exceeded")
                                out.write(block)
                        if expected is not None and size != expected:
                            raise IOError(f"Content-Length {expected} differs from retrieved bytes {size}")
                        temp.replace(destination)
                    used_bytes += size
                    status = "downloaded"
                    error = None
                    break
                except Exception as exc:
                    error = f"{type(exc).__name__}: {exc}"
                    time.sleep(1.0 + attempt)
            if error is not None:
                raise RuntimeError(f"Could not retrieve {filename}: {error}")
        results.append({
            "GSM": row["GSM"],
            "GEO_series": row["GEO_series"],
            "Amp.Batch.ID": row["Amp.Batch.ID"],
            "filename": filename,
            "bytes": size,
            "sha256": sha256(destination),
            "download_url": row["download_url"],
            "status": status,
        })
        if index % 25 == 0 or index == len(files):
            print(json.dumps({"processed": index, "total": len(files), "cumulative_data_bytes": used_bytes}), flush=True)

    fields = ("GSM", "GEO_series", "Amp.Batch.ID", "filename", "bytes", "sha256", "download_url", "status")
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        writer.writerows(results)
    summary = {
        "selected_files": len(files),
        "downloaded_files": sum(result["status"] == "downloaded" for result in results),
        "preexisting_selected_files": sum(result["status"] == "already_present" for result in results),
        "selected_processed_table_bytes": sum(result["bytes"] for result in results),
        "all_next_stage_data_bytes": used_bytes,
        "budget_bytes": MAX_TOTAL,
        "budget_pass": used_bytes <= MAX_TOTAL,
        "table_manifest": str(OUT),
    }
    SUMMARY.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
