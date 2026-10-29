#!/usr/bin/env python3
"""Verify packaged file hashes. Does not load study data or execute analyses."""
from pathlib import Path
import hashlib
import json
import sys

root = Path(__file__).resolve().parent.parent
manifest = json.loads((root / "FILE_MANIFEST.json").read_text())
failures = []
for entry in manifest["files"]:
    path = root / entry["path"]
    if not path.is_file():
        failures.append(entry["path"] + ": missing")
        continue
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != entry["sha256"]:
        failures.append(entry["path"] + ": changed")
if failures:
    print("\n".join(failures))
    sys.exit(1)
print(f"Verified {len(manifest['files'])} packaged files. This is not a statistical validation.")
