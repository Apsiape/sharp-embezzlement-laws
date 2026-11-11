"""Verify every distributed payload against SHA256SUMS.json."""
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath


def main():
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / "SHA256SUMS.json").read_text())
    failures = []
    for name, expected in manifest.items():
        relative = PurePosixPath(name)
        if relative.is_absolute() or ".." in relative.parts:
            raise SystemExit("Invalid manifest path")
        path = root / relative
        if not path.is_file():
            failures.append(f"Missing: {name}")
        elif sha256(path.read_bytes()).hexdigest() != expected:
            failures.append(f"Changed: {name}")
    if failures:
        raise SystemExit("\n".join(failures))
    print(f"Verified {len(manifest)} distributed file checksums")


if __name__ == "__main__":
    main()
