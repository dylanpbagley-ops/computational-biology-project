"""Build a deterministic SHA-256 inventory of next-stage outputs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

NEXT = Path(__file__).resolve().parents[1]
OUTPUT = NEXT / "delivery" / "output_hashes.json"
EXCLUDED_PARTS = {"data", "dependencies", "__pycache__", ".cache", "cache"}


def digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            hasher.update(block)
    return hasher.hexdigest()


def main() -> None:
    files = []
    for path in sorted(NEXT.rglob("*")):
        if not path.is_file() or path == OUTPUT or path.is_symlink():
            continue
        relative = path.relative_to(NEXT)
        if any(part in EXCLUDED_PARTS for part in relative.parts) or path.suffix == ".pyc":
            continue
        files.append({"path": relative.as_posix(), "bytes": path.stat().st_size, "sha256": digest(path)})
    record = {
        "manifest_scope": "all regular next_stage files except this manifest, downloaded data, dependencies and caches",
        "excluded_parts": sorted(EXCLUDED_PARTS),
        "downloaded_data_integrity_source": "analysis/results/external_full_actual_download_manifest.tsv (SHA-256 recorded in reproduction.md)",
        "file_count": len(files),
        "files": files,
    }
    OUTPUT.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {OUTPUT} with {len(files)} files")


if __name__ == "__main__":
    main()
