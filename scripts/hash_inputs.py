"""Verify local source files against the stored SHA-256 digests."""

import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import anchormtm_paths as paths


def main():
    for line in (paths.ROOT / "data/source_sha256.txt").read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        path = paths.DATA / Path(name).relative_to("data")
        if not path.is_file():
            raise SystemExit(f"Missing input: {name}")
        sha256 = hashlib.sha256()
        with path.open("rb") as source:
            for block in iter(lambda: source.read(1024 * 1024), b""):
                sha256.update(block)
        if digest != sha256.hexdigest():
            raise SystemExit(f"SHA-256 mismatch: {name}")
        print(name)


if __name__ == "__main__":
    main()
