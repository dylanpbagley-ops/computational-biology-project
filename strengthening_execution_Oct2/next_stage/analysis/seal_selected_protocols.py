"""Write sidecar hashes for the preserved post-pilot protocol artifacts."""

from __future__ import annotations

import hashlib
from pathlib import Path


HERE = Path(__file__).resolve().parent
FILES = (
    HERE / "selected_protocol_v2_postpilot_serialization_addendum.md",
    HERE / "selected_protocol_v3_full_expansion.md",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    for path in FILES:
        digest = sha256(path)
        path.with_suffix(path.suffix + ".sha256").write_text(
            f"{digest}  {path.name}\n", encoding="ascii"
        )
        print(f"{path.name}\t{digest}")


if __name__ == "__main__":
    main()
