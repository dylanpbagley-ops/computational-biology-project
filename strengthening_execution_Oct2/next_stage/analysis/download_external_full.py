"""Retrieve the locked metadata-selected GEO processed tables within the approved size cap."""

from __future__ import annotations

import concurrent.futures
import csv
import hashlib
import json
import urllib.request
from pathlib import Path


HERE = Path(__file__).resolve().parent
NEXT = HERE.parent
DATA = NEXT / "data"
RESULTS = HERE / "results"
LOCK_JSON = RESULTS / "external_full_selection_lock.json"
SIZE_LIST = RESULTS / "external_full_download_manifest_locked.tsv"
PILOT_DOWNLOADS = RESULTS / "external_pilot_download_manifest.tsv"
ACTUAL_TSV = RESULTS / "external_full_actual_download_manifest.tsv"
SUMMARY_JSON = RESULTS / "external_full_download.json"
PROTOCOL = HERE / "selected_protocol_v3_full_expansion.md"
PROTOCOL_SHA = "62eb79d20342c5fa220b75375c3b847108873830486664583ee2ed4fe31a4461"
LOCK_SHA = "c67c4e8e2fccf1ca0b3a36e61a5064b49cc446b6e0a411b191fe3f845ddb07c0"
MAX_TOTAL = 500 * 1024 * 1024


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def retrieve(row: dict[str, str]) -> dict[str, object]:
    destination = DATA / row["filename"]
    expected = int(row["expected_bytes"])
    if destination.exists():
        if destination.stat().st_size != expected:
            raise ValueError(f"Preexisting size differs from lock: {destination.name}")
        return {**row, "bytes": expected, "sha256": sha256(destination), "status": "already_present"}

    temp = destination.with_suffix(destination.suffix + ".part")
    if temp.exists():
        raise ValueError(f"Stale partial download exists; inspect before retry: {temp.name}")
    request = urllib.request.Request(
        row["download_url"], headers={"User-Agent": "Codex public data research/1.0"}
    )
    try:
        with urllib.request.urlopen(request, timeout=90) as response:
            if response.status != 200:
                raise IOError(f"Expected HTTP 200, received {response.status}")
            content_length = response.headers.get("Content-Length")
            if content_length is None or not content_length.isdigit() or int(content_length) != expected:
                raise IOError(f"GET Content-Length {content_length!r} differs from locked size {expected}")
            size = 0
            with temp.open("wb") as out:
                while True:
                    block = response.read(1 << 20)
                    if not block:
                        break
                    size += len(block)
                    if size > expected:
                        raise IOError(f"Downloaded more than locked size for {destination.name}")
                    out.write(block)
            if size != expected:
                raise IOError(f"Downloaded {size} bytes; locked size is {expected}")
        temp.replace(destination)
    except Exception:
        if temp.exists():
            temp.unlink()
        raise
    return {**row, "bytes": size, "sha256": sha256(destination), "status": "downloaded"}


def main() -> None:
    if sha256(PROTOCOL) != PROTOCOL_SHA:
        raise ValueError("Full expansion protocol hash mismatch")
    lock = json.loads(LOCK_JSON.read_text(encoding="utf-8"))
    expected_inputs = {
        "full_manifest_sha256": "e3e21adee57ab689b6a69086b4e288d0e7176c7e733d3c01d55e718db29dcbe1",
        "download_manifest_sha256": "c2e8e14beecd2b7aca89caa9e94adfb814c455f92951a4e42ac4a5742fad40ad",
        "selection_lock_sha256": LOCK_SHA,
    }
    if lock["locked_metadata_manifest_sha256"] != expected_inputs["full_manifest_sha256"]:
        raise ValueError("Locked metadata manifest no longer matches the protocol frame")
    if lock["locked_download_manifest_sha256"] != expected_inputs["download_manifest_sha256"]:
        raise ValueError("Locked download file list no longer matches the protocol frame")
    if sha256(LOCK_JSON) != expected_inputs["selection_lock_sha256"]:
        raise ValueError("Full selection lock hash mismatch")
    if not lock["budget_pass"] or int(lock["estimated_final_data_bytes"]) > MAX_TOTAL:
        raise ValueError("Locked full dataset no longer passes the 500 MiB data cap")

    pilot_rows = read_tsv(PILOT_DOWNLOADS)
    pilot_hashes = {row["filename"]: row for row in pilot_rows}
    if sha256(PILOT_DOWNLOADS) != lock["pilot_download_manifest_sha256"]:
        raise ValueError("Pilot download manifest hash mismatch")
    rows = read_tsv(SIZE_LIST)
    if sha256(SIZE_LIST) != expected_inputs["download_manifest_sha256"]:
        raise ValueError("Locked full-size file list hash mismatch")
    if len(rows) != lock["unique_processed_table_files"]:
        raise ValueError("Locked full-size file list count mismatch")
    expected_by_name = {row["filename"]: row for row in rows}
    if len(expected_by_name) != len(rows):
        raise ValueError("Duplicate filenames in the exact locked full-size list")
    current_data_bytes = sum(path.stat().st_size for path in DATA.iterdir() if path.is_file())
    missing_bytes = sum(
        int(row["expected_bytes"]) for row in rows if not (DATA / row["filename"]).exists()
    )
    if current_data_bytes + missing_bytes > MAX_TOTAL:
        raise ValueError("Current directory plus locked missing tables would exceed the cap")

    def fetch(row: dict[str, str]) -> dict[str, object]:
        result = retrieve(row)
        previous = pilot_hashes.get(row["filename"])
        if previous:
            if result["status"] != "already_present":
                raise ValueError(f"Pilot table unexpectedly absent before expansion: {row['filename']}")
            if result["sha256"] != previous["sha256"] or int(result["bytes"]) != int(previous["bytes"]):
                raise ValueError(f"Pilot-table preservation hash failed: {row['filename']}")
            result["status"] = "pilot_preexisting_verified"
        elif result["status"] == "already_present":
            result["status"] = "recovered_after_manifest_write_failure"
        return result

    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        future_by_name = {
            executor.submit(fetch, row): row["filename"] for row in rows
        }
        results: list[dict[str, object]] = []
        for index, future in enumerate(concurrent.futures.as_completed(future_by_name), start=1):
            filename = future_by_name[future]
            try:
                results.append(future.result())
            except Exception as exc:
                raise RuntimeError(f"Bounded retrieval failed for {filename}: {type(exc).__name__}: {exc}") from exc
            if index % 40 == 0 or index == len(rows):
                total = sum(int(row["bytes"]) for row in results)
                print(json.dumps({"verified": index, "total": len(rows), "selected_table_bytes_verified": total}), flush=True)

    results.sort(key=lambda row: str(row["filename"]))
    fields = ("filename", "GSM", "Amp.Batch.ID", "download_url", "expected_bytes", "bytes", "sha256", "status", "error")
    with ACTUAL_TSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        writer.writeheader()
        writer.writerows(results)
    actual_data_bytes = sum(path.stat().st_size for path in DATA.iterdir() if path.is_file())
    summary = {
        "protocol_sha256": PROTOCOL_SHA,
        "full_selection_lock_sha256": sha256(LOCK_JSON),
        "locked_manifest_sha256": lock["locked_metadata_manifest_sha256"],
        "locked_download_manifest_sha256": lock["locked_download_manifest_sha256"],
        "unique_processed_tables_verified": len(results),
        "new_processed_tables_downloaded": sum(row["status"] in {"downloaded", "recovered_after_manifest_write_failure"} for row in results),
        "preexisting_tables_verified": sum(row["status"] == "pilot_preexisting_verified" for row in results),
        "selected_table_bytes_verified": sum(int(row["bytes"]) for row in results),
        "actual_next_stage_data_bytes": actual_data_bytes,
        "cap_bytes": MAX_TOTAL,
        "cap_pass": actual_data_bytes <= MAX_TOTAL,
        "actual_table_manifest_sha256": sha256(ACTUAL_TSV),
        "output_manifest": str(ACTUAL_TSV),
    }
    ACTUAL_TSV.with_suffix(ACTUAL_TSV.suffix + ".sha256").write_text(
        f"{summary['actual_table_manifest_sha256']}  {ACTUAL_TSV.name}\n", encoding="ascii"
    )
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
