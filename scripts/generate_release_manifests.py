#!/usr/bin/env python3
"""
Generate MANIFEST.tsv and MANIFEST.json for the public release.

Run this after all archives (user_input.tar.gz and per-analysis .tar.gz)
are finalized.

Output:
  - MANIFEST.tsv: Two-column manifest (path, size_bytes)
  - MANIFEST.json: Structured manifest for catalog browsing

Usage:
  python3 scripts/generate_release_manifests.py
"""

import os
import json
import subprocess
from pathlib import Path
from collections import OrderedDict

MASTER = Path("/scratch/ims86/NCLEdetector_Datastore")
MANIFEST_TSV = MASTER / "MANIFEST.tsv"
MANIFEST_JSON = MASTER / "MANIFEST.json"

EXCLUDE_PATTERNS = {
    "MANIFEST.tsv",
    "MANIFEST.json",
    "SHA256SUMS",
    "outputs/workflow2/OP.tar.gz",
}


def should_exclude(path_str):
    """Check if path matches any exclude pattern."""
    for pattern in EXCLUDE_PATTERNS:
        if path_str.endswith(pattern):
            return True
    return False


def get_archive_info(filename):
    """Determine archive relationship."""
    if filename == "user_input.tar.gz":
        return {"archive_type": "top-level", "is_archive": True}
    elif filename.endswith(".tar.gz"):
        base = filename[:-7]  # Remove .tar.gz
        if "_" in base:
            # Try to extract analysis name from workflow pattern
            parts = base.split("_", 1)  # Split on first underscore
            if len(parts) == 2 and parts[0].startswith("workflow"):
                return {
                    "archive_type": "analysis",
                    "is_archive": True,
                    "analysis_name": parts[1],
                }
        return {"archive_type": "other", "is_archive": True}
    else:
        return {"archive_type": None, "is_archive": False}


def format_size_human(size_bytes):
    """Format bytes to human-readable size."""
    for unit in ["B", "KiB", "MiB", "GiB", "TiB"]:
        if size_bytes < 1024:
            return f"{size_bytes:.1f}{unit}"
        size_bytes /= 1024
    return f"{size_bytes:.1f}PiB"


def main():
    print("Generating release manifests...")
    print()

    # ========================================================================
    # Generate MANIFEST.tsv
    # ========================================================================
    print("Creating MANIFEST.tsv...")

    entries = []
    for root, dirs, files in os.walk(MASTER):
        for filename in sorted(files):
            filepath = Path(root) / filename
            rel_path = filepath.relative_to(MASTER)
            rel_path_str = str(rel_path)

            if should_exclude(rel_path_str):
                continue

            try:
                size_bytes = filepath.stat().st_size
                entries.append((rel_path_str, size_bytes))
            except (OSError, FileNotFoundError):
                continue

    # Sort entries lexically
    entries.sort(key=lambda x: x[0])

    # Write TSV
    with open(MANIFEST_TSV, "w") as f:
        for path, size in entries:
            f.write(f"{path}\t{size}\n")

    print(f"✓ MANIFEST.tsv created with {len(entries)} entries")
    print()

    # ========================================================================
    # Generate MANIFEST.json
    # ========================================================================
    print("Creating MANIFEST.json...")

    manifest_data = []
    for path_str, size_bytes in entries:
        if path_str.startswith("MANIFEST") or path_str == "SHA256SUMS":
            continue

        filename = Path(path_str).name
        archive_info = get_archive_info(filename)

        entry = OrderedDict([
            ("path", path_str),
            ("filename", filename),
            ("size_bytes", size_bytes),
            ("size_human", format_size_human(size_bytes)),
        ])
        entry.update(archive_info)
        manifest_data.append(entry)

    # Write JSON
    with open(MANIFEST_JSON, "w") as f:
        json.dump(manifest_data, f, indent=2)

    with open(MANIFEST_JSON, "r") as f:
        json_lines = len(f.readlines())

    print(f"✓ MANIFEST.json created with {json_lines} lines")
    print()

    # ========================================================================
    # Summary
    # ========================================================================
    print("=" * 80)
    print("Manifest Generation Complete")
    print("=" * 80)
    print()
    print("Files created:")
    print(f"  - {MANIFEST_TSV} ({len(entries)} entries)")
    print(f"  - {MANIFEST_JSON}")
    print()
    print("Sample entries from MANIFEST.tsv:")
    with open(MANIFEST_TSV) as f:
        for i, line in enumerate(f):
            if i >= 5:
                break
            print(f"  {line.rstrip()}")
    print()
    print("Sample from MANIFEST.json (first 3 entries):")
    with open(MANIFEST_JSON) as f:
        data = json.load(f)
        for entry in data[:3]:
            print(f"  {entry['filename']:50s} {entry['size_human']:>10s}")


if __name__ == "__main__":
    main()
