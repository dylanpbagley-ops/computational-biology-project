"""Capture immutable source hashes for the October 3 extension stage."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
STAGE = Path(__file__).resolve().parents[1]
OLD = ROOT / "strengthening_execution_Oct2"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def file_record(path: Path, relative: str | None = None) -> dict[str, object]:
    return {
        "path": relative or str(path.relative_to(ROOT)).replace("\\", "/"),
        "bytes": path.stat().st_size,
        "sha256": digest(path),
    }


def main() -> None:
    destination = STAGE / "preservation" / "input_hashes.json"
    if destination.exists():
        raise SystemExit(f"Refusing to replace starting-state record: {destination}")

    previous_lock_path = OLD / "next_stage" / "input_hashes.json"
    prior_output_path = OLD / "next_stage" / "delivery" / "output_hashes.json"
    protected_manifest_path = OLD / "input_hashes.json"
    lock = json.loads(previous_lock_path.read_text(encoding="utf-8"))
    output_manifest = json.loads(prior_output_path.read_text(encoding="utf-8"))
    protected_manifest = json.loads(protected_manifest_path.read_text(encoding="utf-8"))

    if len(lock["protected_inputs"]) != 63 or len(lock["current_execution_snapshot"]) != 65:
        raise SystemExit("Prior protection manifest counts differ from the approved starting values")
    if len(output_manifest["files"]) != 155:
        raise SystemExit("Prior output inventory count differs from the approved starting value")
    if digest(prior_output_path) != "b63865ce9249efbaf19e3f275b5e5d773f2ffdb9c51730c8f376b51b7d8c4305":
        raise SystemExit("Prior output inventory manifest SHA-256 differs from the handoff")

    protected = []
    for relative, expected in lock["protected_inputs"].items():
        path = ROOT / relative
        actual = digest(path)
        protected.append({"path": relative, "expected_sha256": expected,
                          "start_sha256": actual, "start_match": actual == expected,
                          "bytes": path.stat().st_size})

    snapshot = []
    for relative, expected in lock["current_execution_snapshot"].items():
        path = ROOT / relative
        actual = digest(path)
        snapshot.append({"path": relative, "prior_snapshot_sha256": expected,
                         "start_sha256": actual, "prior_snapshot_match": actual == expected,
                         "bytes": path.stat().st_size})

    outputs = []
    for row in output_manifest["files"]:
        relative = "strengthening_execution_Oct2/next_stage/" + row["path"]
        path = ROOT / relative
        actual = digest(path)
        outputs.append({"path": relative, "bytes": path.stat().st_size,
                        "prior_inventory_bytes": row["bytes"],
                        "start_sha256": actual, "prior_inventory_sha256": row["sha256"],
                        "start_match": actual == row["sha256"] and path.stat().st_size == row["bytes"]})

    direct_paths = [
        "docs/superpowers/plans/2026-10-03-longitudinal-target-extension.md",
        "EV-Capture-CURRENT.md",
        "strengthening_execution_Oct2/input_hashes.json",
        "strengthening_execution_Oct2/next_stage/input_hashes.json",
        "strengthening_execution_Oct2/next_stage/delivery/output_hashes.json",
        "strengthening_execution_Oct2/next_stage/feasibility/contribution_gate.md",
        "strengthening_execution_Oct2/next_stage/feasibility/path_decision_v2.md",
        "strengthening_execution_Oct2/next_stage/feasibility/fallback_cards.md",
        "strengthening_execution_Oct2/next_stage/feasibility/prior_art_delta.tsv",
        "strengthening_execution_Oct2/next_stage/feasibility/external_data_assessment.md",
        "strengthening_execution_Oct2/next_stage/delivery/reproduction.md",
        "strengthening_execution_Oct2/next_stage/analysis/results/external_annotation_identity_audit.json",
        "strengthening_execution_Oct2/next_stage/analysis/results/external_annotation_by_hospital_biopsy.tsv",
        "strengthening_execution_Oct2/next_stage/data/GSE305971_cells_annotation_for_geo_ver_20250306.xlsx",
        "strengthening_execution_Oct2/next_stage/delivery/root_status_append_record.json",
        "strengthening_execution_Oct2/next_stage/before/PROJECT_STATUS_Oct2.md",
        "strengthening_execution_Oct2/next_stage/before/execution_STATUS.md",
        "strengthening_execution_Oct2/next_stage/before/next_stage_STATUS_before_final_handoff_2026-10-03.md",
    ]
    direct = [file_record(ROOT / relative) for relative in direct_paths]

    append_record_path = OLD / "next_stage" / "delivery" / "root_status_append_record.json"
    append_record = json.loads(append_record_path.read_text(encoding="utf-8"))
    reconciliations = []
    for row in append_record["root_status_records"]:
        current_path = ROOT / row["path"]
        backup_path = OLD / "next_stage" / row["backup"]
        current_bytes = current_path.read_bytes()
        backup_bytes = backup_path.read_bytes()
        suffix = current_bytes[len(backup_bytes):]
        valid = (
            digest(backup_path) == row["original_sha256"] == row["backup_sha256"]
            and digest(current_path) == row["current_sha256"]
            and current_bytes.startswith(backup_bytes)
            and hashlib.sha256(suffix).hexdigest() == row["suffix_sha256"]
            and row["append_prefix_preserved"] is True
        )
        reconciliations.append({"path": row["path"], "backup": row["backup"],
                                "snapshot_sha256": row["original_sha256"],
                                "observed_start_sha256": digest(current_path),
                                "backup_sha256": digest(backup_path),
                                "suffix_sha256": hashlib.sha256(suffix).hexdigest(),
                                "prefix_preserved": current_bytes.startswith(backup_bytes),
                                "append_record_reconciled": valid})

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "root": str(ROOT),
        "starting_inventory": {
            "original_protected_count": len(protected),
            "prior_65_snapshot_count": len(snapshot),
            "prior_155_output_count": len(outputs),
            "original_protected_manifest": file_record(protected_manifest_path),
            "prior_snapshot_manifest": file_record(previous_lock_path),
            "prior_output_manifest": file_record(prior_output_path),
            "prior_output_manifest_handoff_sha256": "b63865ce9249efbaf19e3f275b5e5d773f2ffdb9c51730c8f376b51b7d8c4305",
        },
        "protected_inputs": protected,
        "prior_snapshot_65": snapshot,
        "prior_outputs_155": outputs,
        "directly_used_prior_files": direct,
        "status_suffix_reconciliation": reconciliations,
        "preexisting_mismatches": [row["path"] for row in snapshot if not row["prior_snapshot_match"]],
        "notes": [
            "The original 63 protected-input hashes are checked against the prior lock.",
            "The observed start hashes, rather than stale prior-snapshot hashes, are the preservation baseline for all 65 snapshot files.",
            "Two status files differ from the old snapshot; the existing append record and backups are checked for exact suffix-only appends.",
            "The 155 inventoried prior outputs are checked against both recorded sizes and hashes.",
        ],
    }
    if not all(row["start_match"] for row in protected):
        raise SystemExit("A protected input differs at stage start; see no output was written")
    if not all(row["start_match"] for row in outputs):
        raise SystemExit("A prior output differs at stage start; see no output was written")
    if not all(row["append_record_reconciled"] for row in reconciliations):
        raise SystemExit("A preexisting status append did not reconcile to its backup record")
    if len(record["preexisting_mismatches"]) != 2:
        raise SystemExit("Unexpected mismatch count in prior 65-file snapshot")

    destination.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "protected_inputs_match": len(protected),
        "prior_snapshot_files_captured": len(snapshot),
        "preexisting_status_suffix_mismatches": record["preexisting_mismatches"],
        "suffix_reconciliations_pass": all(row["append_record_reconciled"] for row in reconciliations),
        "prior_outputs_match": len(outputs),
        "direct_references_hashed": len(direct),
        "input_hashes_sha256": digest(destination),
    }, indent=2))


if __name__ == "__main__":
    main()
