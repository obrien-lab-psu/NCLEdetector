#!/bin/bash
################################################################################
# Generate SHA256SUMS for all release files
#
# Run this after archives and manifests are finalized.
#
# Output:
#   - SHA256SUMS: SHA-256 checksums for all files except themselves
#
# Usage:
#   bash scripts/generate_release_checksums.sh
#
################################################################################

set -eo pipefail

MASTER="/scratch/ims86/NCLEdetector_Datastore"
CHECKSUMS="$MASTER/SHA256SUMS"

echo "Generating SHA256SUMS..."
echo ""

cd "$MASTER"

# Generate checksums for all files except the checksums file itself, 
# manifests, and other metadata
find . -type f \
    ! -path './MANIFEST.tsv' \
    ! -path './MANIFEST.json' \
    ! -path './SHA256SUMS' \
    -print0 \
    | LC_ALL=C sort -z \
    | xargs -0 sha256sum > "$CHECKSUMS"

CHECKSUM_COUNT=$(wc -l < "$CHECKSUMS")

echo "✓ SHA256SUMS created with $CHECKSUM_COUNT entries"
echo ""
echo "File size: $(du -h $CHECKSUMS | cut -f1)"
echo ""
echo "Sample checksums:"
head -5 "$CHECKSUMS"
echo ""
echo "To verify checksums later, run:"
echo "  cd $MASTER"
echo "  sha256sum -c SHA256SUMS"
