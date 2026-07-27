#!/bin/bash
################################################################################
# Generate MANIFEST.tsv and MANIFEST.json for the public release
#
# Run this after all archives (user_input.tar.gz and per-analysis .tar.gz)
# are finalized.
#
# Output:
#   - MANIFEST.tsv: Two-column manifest (path, size_bytes)
#   - MANIFEST.json: Structured manifest for catalog browsing
#
# Usage:
#   bash scripts/generate_release_manifests.sh
#
################################################################################

set -eo pipefail

MASTER="/scratch/ims86/NCLEdetector_Datastore"
MANIFEST_TSV="$MASTER/MANIFEST.tsv"
MANIFEST_JSON="$MASTER/MANIFEST.json"

echo "Generating release manifests..."
echo ""

# ============================================================================
# Generate MANIFEST.tsv
# ============================================================================

echo "Creating MANIFEST.tsv..."

find "$MASTER" -type f \
    ! -path './MANIFEST.tsv' \
    ! -path './MANIFEST.json' \
    ! -path './SHA256SUMS' \
    ! -path './outputs/workflow2/OP.tar.gz' \
    -printf '%P\t%s\n' \
    | LC_ALL=C sort > "$MANIFEST_TSV"

TSV_LINES=$(wc -l < "$MANIFEST_TSV")
echo "✓ MANIFEST.tsv created with $TSV_LINES entries"
echo ""

# ============================================================================
# Generate MANIFEST.json
# ============================================================================

echo "Creating MANIFEST.json..."

# Function to determine archive relationship
get_archive_info() {
    local path="$1"
    local filename=$(basename "$path")
    
    if [[ "$filename" == "user_input.tar.gz" ]]; then
        echo '{"archive_type": "top-level", "is_archive": true}'
    elif [[ "$filename" =~ \.tar\.gz$ ]]; then
        # Extract workflow and analysis name
        local base="${filename%.tar.gz}"
        if [[ "$base" =~ ^workflow[0-9]+_(.+)$ ]]; then
            local analysis="${BASH_REMATCH[1]}"
            echo "{\"archive_type\": \"analysis\", \"is_archive\": true, \"analysis_name\": \"$analysis\"}"
        else
            echo '{"archive_type": "other", "is_archive": true}'
        fi
    else
        echo '{"archive_type": null, "is_archive": false}'
    fi
}

cat > "$MANIFEST_JSON" << 'JSONSTART'
[
JSONSTART

first=true
while IFS=$'\t' read -r path size_bytes; do
    if [[ "$path" == "MANIFEST"* ]] || [[ "$path" == "SHA256SUMS" ]]; then
        continue
    fi
    
    # Add comma before each entry except the first
    if [[ "$first" == "false" ]]; then
        echo "," >> "$MANIFEST_JSON"
    fi
    first=false
    
    filename=$(basename "$path")
    archive_info=$(get_archive_info "$path")
    
    cat >> "$MANIFEST_JSON" << JSONENTRY
  {
    "path": "$path",
    "filename": "$filename",
    "size_bytes": $size_bytes,
    "size_human": "$(numfmt --to=iec-i --suffix=B $size_bytes 2>/dev/null || echo 'N/A')",
    $archive_info
  }
JSONENTRY
done < "$MANIFEST_TSV"

echo "" >> "$MANIFEST_JSON"
echo "]" >> "$MANIFEST_JSON"

JSON_LINES=$(wc -l < "$MANIFEST_JSON")
echo "✓ MANIFEST.json created with $JSON_LINES lines"
echo ""

# ============================================================================
# Summary
# ============================================================================

echo "================================================================================"
echo "Manifest Generation Complete"
echo "================================================================================"
echo ""
echo "Files created:"
echo "  - $MANIFEST_TSV ($TSV_LINES entries)"
echo "  - $MANIFEST_JSON"
echo ""
echo "Sample entries from MANIFEST.tsv:"
head -5 "$MANIFEST_TSV"
echo ""
echo "Sample from MANIFEST.json:"
head -20 "$MANIFEST_JSON"
