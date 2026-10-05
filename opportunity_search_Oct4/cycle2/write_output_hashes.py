"""Create the final SHA-256 inventory for cycle2 outputs."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

STAGE = Path(__file__).resolve().parent
CACHE_PARTS = {".cache", "__pycache__", "cache"}

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def main() -> None:
    files = []
    for path in sorted(STAGE.rglob("*")):
        if not path.is_file() or path.name == "output_hashes.json":
            continue
        if any(part in CACHE_PARTS for part in path.relative_to(STAGE).parts):
            continue
        files.append({
            "path": path.relative_to(STAGE).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    target = STAGE / "output_hashes.json"
    target.write_text(json.dumps({"schema_version": 1, "file_count": len(files), "files": files}, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"file_count": len(files), "hash_file": target.relative_to(STAGE).as_posix()}, indent=2))

if __name__ == "__main__":
    main()
