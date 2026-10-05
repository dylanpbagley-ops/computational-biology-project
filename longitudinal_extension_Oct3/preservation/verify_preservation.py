"""Verify all protected and prior files against the captured stage-start hashes."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RECORD = Path(__file__).resolve().parent / "input_hashes.json"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify_group(rows: list[dict], label: str, failures: list[dict]) -> int:
    for row in rows:
        path = ROOT / row["path"]
        expected = row.get("start_sha256", row.get("sha256"))
        observed = digest(path) if path.is_file() else None
        if observed != expected:
            failures.append({"group": label, "path": row["path"],
                             "expected_start_sha256": expected, "observed_sha256": observed})
    return len(rows)


def main() -> None:
    record = json.loads(RECORD.read_text(encoding="utf-8"))
    failures: list[dict] = []
    total = 0
    total += verify_group(record["protected_inputs"], "protected_inputs", failures)
    total += verify_group(record["prior_snapshot_65"], "prior_snapshot_65", failures)
    total += verify_group(record["prior_outputs_155"], "prior_outputs_155", failures)
    total += verify_group(record["directly_used_prior_files"], "directly_used_prior_files", failures)
    for row in record["status_suffix_reconciliation"]:
        current = ROOT / row["path"]
        backup = ROOT / "strengthening_execution_Oct2" / "next_stage" / row["backup"]
        current_bytes = current.read_bytes() if current.is_file() else b""
        backup_bytes = backup.read_bytes() if backup.is_file() else b""
        suffix = current_bytes[len(backup_bytes):] if current_bytes.startswith(backup_bytes) else b""
        valid = (
            digest(current) == row["observed_start_sha256"]
            and digest(backup) == row["backup_sha256"]
            and current_bytes.startswith(backup_bytes)
            and hashlib.sha256(suffix).hexdigest() == row["suffix_sha256"]
        )
        if not valid:
            failures.append({"group": "status_suffix_reconciliation", "path": row["path"]})
    result = {"files_checked_with_duplicates": total,
              "status_suffix_reconciliations": len(record["status_suffix_reconciliation"]),
              "preexisting_snapshot_mismatches": record["preexisting_mismatches"],
              "failures": failures,
              "preservation_pass": not failures}
    print(json.dumps(result, indent=2, allow_nan=False))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
