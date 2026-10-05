"""Append final pointers only after verifying the preserved original prefixes."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
NEXT = Path(__file__).resolve().parents[1]
BACKUP_DIR = NEXT / "before"
SNAPSHOT = NEXT / "input_hashes.json"
NEXT_STATUS = NEXT / "STATUS.md"
NEXT_STATUS_BACKUP = BACKUP_DIR / "next_stage_STATUS_before_final_handoff_2026-10-03.md"
RECORD = NEXT / "delivery" / "root_status_append_record.json"


def digest(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def file_digest(path: Path) -> str:
    return digest(path.read_bytes())


def append_bytes(path: Path, suffix: bytes) -> bytes:
    original = path.read_bytes()
    path.write_bytes(original + suffix)
    if not path.read_bytes().startswith(original):
        raise ValueError(f"Append-only prefix check failed for {path}")
    return original


def main() -> None:
    if RECORD.exists() or NEXT_STATUS_BACKUP.exists():
        raise FileExistsError("Final status append already recorded; refusing to append a duplicate pointer")
    snapshot = json.loads(SNAPSHOT.read_text(encoding="utf-8"))["current_execution_snapshot"]
    targets = [
        ("PROJECT_STATUS_Oct2.md", "PROJECT_STATUS_Oct2.md", "\r\n## Next-stage completion pointer — October 3, 2026\r\n\r\nFinal gate: no demonstrated strong extension. Completed evidence, policy analysis, F1/F2 screens and clinical-source normalization (clinical review pending) are recorded at `strengthening_execution_Oct2/next_stage/feasibility/path_decision_v2.md` and `strengthening_execution_Oct2/next_stage/delivery/reproduction.md`. Source results and file hashes are indexed in `strengthening_execution_Oct2/next_stage/delivery/output_hashes.json`. Original labels and submitted prose were not changed.\r\n".encode("utf-8")),
        ("strengthening_execution_Oct2/STATUS.md", "execution_STATUS.md", "\r\n## Final next-stage pointer — October 3, 2026\r\n\r\nThe bounded extension execution is complete. Contribution/path decision: no demonstrated strong extension. See `next_stage/feasibility/path_decision_v2.md`, `next_stage/feasibility/contribution_gate.md`, `next_stage/delivery/reproduction.md`, and `next_stage/delivery/output_hashes.json`. Parent numerical/source review is recorded; independent clinical review and student authorship remain pending.\r\n".encode("utf-8")),
    ]
    root_records = []
    for rel_path, backup_name, _ in targets:
        expected = snapshot[rel_path]
        backup = BACKUP_DIR / backup_name
        live = ROOT / Path(rel_path)
        if file_digest(backup) != expected:
            raise ValueError(f"Verified original backup hash does not match snapshot: {backup}")
        if file_digest(live) != expected or live.read_bytes() != backup.read_bytes():
            raise ValueError(f"Refusing to append after an unexpected starting-prefix change: {live}")

    if not NEXT_STATUS.is_file():
        raise FileNotFoundError(NEXT_STATUS)
    original_next_status = NEXT_STATUS.read_bytes()
    NEXT_STATUS_BACKUP.parent.mkdir(parents=True, exist_ok=True)
    NEXT_STATUS_BACKUP.write_bytes(original_next_status)

    for rel_path, backup_name, suffix in targets:
        live = ROOT / Path(rel_path)
        original = append_bytes(live, suffix)
        backup = BACKUP_DIR / backup_name
        root_records.append({
            "path": rel_path.replace("\\", "/"),
            "backup": backup.relative_to(NEXT).as_posix(),
            "original_sha256": digest(original),
            "backup_sha256": file_digest(backup),
            "suffix_sha256": digest(suffix),
            "current_sha256": file_digest(live),
            "append_prefix_preserved": live.read_bytes().startswith(backup.read_bytes()),
        })

    next_suffix = (
        "\n## Final status — October 3, 2026\n\n"
        "Bounded execution is complete; the contribution gate selects no demonstrated strong extension. "
        "The GSE305971 same-study policy comparison is supporting evidence, and F1/F2 screens did not pass as primary branches. "
        "Clinical facts remain descriptive with all 34 endpoint observation states indeterminate pending independent review. "
        "See `feasibility/path_decision_v2.md`, `feasibility/contribution_gate.md`, `delivery/reproduction.md`, and `delivery/output_hashes.json`.\n"
    ).encode("utf-8")
    append_bytes(NEXT_STATUS, next_suffix)

    analysis_dir = NEXT / "analysis"
    sys.path.insert(0, str(analysis_dir))
    import metadata_diagnostics

    protected_count = metadata_diagnostics.verify_input_snapshot()
    if protected_count != 63:
        raise ValueError(f"Expected 63 protected input hashes, verified {protected_count}")

    record = {
        "status": "append_only_pointers_verified",
        "protected_inputs_verified": protected_count,
        "root_status_records": root_records,
        "next_stage_status_backup": NEXT_STATUS_BACKUP.relative_to(NEXT).as_posix(),
        "next_stage_status_backup_sha256": file_digest(NEXT_STATUS_BACKUP),
        "next_stage_status_current_sha256": file_digest(NEXT_STATUS),
        "current_snapshot_append_validation": "passed after append; expected original bytes remain a prefix",
    }
    RECORD.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
