#!/usr/bin/env python3
"""
Generate SHA256SUMS for all release files.

Run this after archives and manifests are finalized.

Output:
  - SHA256SUMS: SHA-256 checksums for all files except themselves

Usage:
  python3 scripts/generate_release_checksums.py
"""

import os
import hashlib
from pathlib import Path

MASTER = Path("/scratch/ims86/NCLEdetector_Datastore")
CHECKSUMS_FILE = MASTER / "SHA256SUMS"

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


def compute_sha256(filepath, chunk_size=8192*1024):
    """Compute SHA256 checksum of a file."""
    sha256_hash = hashlib.sha256()
    try:
        with open(filepath, "rb") as f:
            while True:
                chunk = f.read(chunk_size)
                if not chunk:
                    break
                sha256_hash.update(chunk)
        return sha256_hash.hexdigest()
    except (OSError, IOError):
        return None


def main():
    print("Generating SHA256SUMS...")
    print()

    # Collect all files to checksum
    files_to_checksum = []
    for root, dirs, files in os.walk(MASTER):
        for filename in sorted(files):
            filepath = Path(root) / filename
            rel_path = filepath.relative_to(MASTER)
            rel_path_str = str(rel_path)

            if should_exclude(rel_path_str):
                continue

            files_to_checksum.append((rel_path_str, filepath))

    # Sort lexically
    files_to_checksum.sort(key=lambda x: x[0])

    # Compute checksums and write
    print(f"Computing checksums for {len(files_to_checksum)} files...")
    checksums = []
    for i, (rel_path_str, filepath) in enumerate(files_to_checksum):
        if (i + 1) % 5000 == 0:
            print(f"  {i + 1}/{len(files_to_checksum)}...")

        checksum = compute_sha256(filepath)
        if checksum:
            checksums.append((checksum, rel_path_str))

    print(f"✓ Computed {len(checksums)} checksums")
    print()

    # Write SHA256SUMS file
    with open(CHECKSUMS_FILE, "w") as f:
        for checksum, rel_path_str in checksums:
            f.write(f"{checksum}  {rel_path_str}\n")

    print(f"✓ SHA256SUMS created with {len(checksums)} entries")
    print()
    print(f"File size: {CHECKSUMS_FILE.stat().st_size / (1024*1024):.1f} MiB")
    print()
    print("Sample checksums:")
    for checksum, rel_path_str in checksums[:5]:
        print(f"  {checksum}  {rel_path_str}")
    print()
    print("To verify checksums later, run:")
    print(f"  cd {MASTER}")
    print("  sha256sum -c SHA256SUMS")


if __name__ == "__main__":
    main()
