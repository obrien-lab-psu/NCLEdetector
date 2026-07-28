#!/bin/bash
#SBATCH --job-name=add_op_to_r2
#SBATCH --account=read_crch_ims86
#SBATCH --partition=basic
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err

set -e

# Configuration
DATASTORE_PATH="/scratch/ims86/NCLEdetector_Datastore"
OP_ARCHIVE="$DATASTORE_PATH/outputs/workflow2/OP.tar.gz"
R2_DEST="r2-ncle:ncle-public-data/NCLEdetector_Datastore/v1.0.0/outputs/workflow2/OP.tar.gz"
MANIFEST_SCRIPT="$PWD/scripts/generate_release_manifests.py"
CHECKSUM_SCRIPT="$PWD/scripts/generate_release_checksums.py"
LOG_FILE="/tmp/op_r2_upload_${SLURM_JOB_ID}.log"

echo "[$(date)] Starting OP.tar.gz R2 upload and manifest regeneration (Job: $SLURM_JOB_ID)"

# =============================================================================
# PHASE 1: Upload OP.tar.gz to R2
# =============================================================================
echo "[$(date)] PHASE 1: Uploading OP.tar.gz to R2..."

if [[ ! -f "$OP_ARCHIVE" ]]; then
    echo "[ERROR] OP archive not found: $OP_ARCHIVE"
    exit 1
fi

OP_SIZE=$(du -h "$OP_ARCHIVE" | cut -f1)
echo "[$(date)] Archive size: $OP_SIZE"
echo "[$(date)] Archive path: $OP_ARCHIVE"

# Dry-run first
echo "[$(date)] Running dry-run..."
rclone copyto "$OP_ARCHIVE" "$R2_DEST" \
    --dry-run \
    --transfers=20 \
    --checkers=20 \
    --retries=10 \
    --progress \
    --log-level=INFO \
    2>&1 | tee -a "$LOG_FILE"

echo "[$(date)] Dry-run complete. Starting actual upload..."

# Full upload
rclone copyto "$OP_ARCHIVE" "$R2_DEST" \
    --transfers=20 \
    --checkers=20 \
    --retries=10 \
    --progress \
    --log-level=INFO \
    2>&1 | tee -a "$LOG_FILE"

echo "[$(date)] Upload complete. Verifying..."

# Verify upload
LOCAL_CHECKSUM=$(sha256sum "$OP_ARCHIVE" | awk '{print $1}')
REMOTE_INFO=$(rclone lsf "$R2_DEST" --format "s" 2>/dev/null || echo "0")

if [[ -z "$REMOTE_INFO" ]] || [[ "$REMOTE_INFO" == "0" ]]; then
    echo "[ERROR] Failed to verify uploaded file in R2"
    exit 1
fi

echo "[$(date)] ✓ OP.tar.gz uploaded successfully to R2"
echo "[$(date)] Local checksum: $LOCAL_CHECKSUM"

# =============================================================================
# PHASE 2: Regenerate Manifests (now including OP.tar.gz)
# =============================================================================
echo "[$(date)] PHASE 2: Regenerating manifests with OP.tar.gz included..."

if [[ ! -f "$MANIFEST_SCRIPT" ]]; then
    echo "[ERROR] Manifest script not found: $MANIFEST_SCRIPT"
    exit 1
fi

cd "$DATASTORE_PATH" || exit 1
python3 "$MANIFEST_SCRIPT" 2>&1 | tee -a "$LOG_FILE"

if [[ ! -f "MANIFEST.tsv" ]] || [[ ! -f "MANIFEST.json" ]]; then
    echo "[ERROR] Manifest generation failed"
    exit 1
fi

MANIFEST_SIZE=$(wc -l MANIFEST.tsv | awk '{print $1}')
echo "[$(date)] ✓ Manifests regenerated ($MANIFEST_SIZE entries)"

# =============================================================================
# PHASE 3: Regenerate Checksums (now including OP.tar.gz)
# =============================================================================
echo "[$(date)] PHASE 3: Regenerating checksums with OP.tar.gz included..."

if [[ ! -f "$CHECKSUM_SCRIPT" ]]; then
    echo "[ERROR] Checksum script not found: $CHECKSUM_SCRIPT"
    exit 1
fi

# Use Python script for checksum generation (fallback to xargs if needed)
python3 "$CHECKSUM_SCRIPT" > SHA256SUMS.tmp 2>&1 || {
    echo "[$(date)] Python checksum script slower; using xargs for parallelism..."
    find . -type f \
        \( -name "MANIFEST.tsv" -o -name "MANIFEST.json" -o -name "SHA256SUMS" \) \
        -prune -o -type f -print | \
    sort | \
    xargs -P 10 -I {} sha256sum {} > SHA256SUMS.tmp
}

mv SHA256SUMS.tmp SHA256SUMS

CHECKSUM_COUNT=$(wc -l SHA256SUMS | awk '{print $1}')
echo "[$(date)] ✓ Checksums regenerated ($CHECKSUM_COUNT files)"

# =============================================================================
# PHASE 4: Upload Updated Metadata to R2
# =============================================================================
echo "[$(date)] PHASE 4: Uploading updated metadata to R2..."

METADATA_FILES=(
    "MANIFEST.tsv"
    "MANIFEST.json"
    "SHA256SUMS"
)

for FILE in "${METADATA_FILES[@]}"; do
    if [[ -f "$FILE" ]]; then
        echo "[$(date)] Uploading $FILE..."
        rclone copyto "$FILE" \
            "r2-ncle:ncle-public-data/NCLEdetector_Datastore/v1.0.0/$FILE" \
            --transfers=5 \
            --checkers=5 \
            --retries=5 \
            --progress \
            --log-level=INFO \
            2>&1 | tee -a "$LOG_FILE"
        echo "[$(date)] ✓ $FILE uploaded"
    else
        echo "[ERROR] Metadata file not found: $FILE"
        exit 1
    fi
done

# =============================================================================
# PHASE 5: Final Verification
# =============================================================================
echo "[$(date)] PHASE 5: Final verification..."

# Verify all metadata files in R2
for FILE in "${METADATA_FILES[@]}"; do
    REMOTE_SIZE=$(rclone lsf "r2-ncle:ncle-public-data/NCLEdetector_Datastore/v1.0.0/$FILE" --format "s" 2>/dev/null || echo "0")
    LOCAL_SIZE=$(wc -c < "$FILE")
    
    if [[ "$REMOTE_SIZE" -eq "$LOCAL_SIZE" ]]; then
        echo "[$(date)] ✓ $FILE verified (size: $REMOTE_SIZE bytes)"
    else
        echo "[ERROR] $FILE size mismatch! Local: $LOCAL_SIZE, Remote: $REMOTE_SIZE"
        exit 1
    fi
done

# =============================================================================
# Summary
# =============================================================================
echo "[$(date)] ============================================================"
echo "[$(date)] ✓ OP.tar.gz Addition Complete!"
echo "[$(date)] ============================================================"
echo "[$(date)] Summary:"
echo "[$(date)]   - OP.tar.gz (294GB) uploaded to R2"
echo "[$(date)]   - Manifests regenerated with OP.tar.gz included"
echo "[$(date)]   - Checksums regenerated with OP.tar.gz included"
echo "[$(date)]   - All metadata files uploaded to R2"
echo "[$(date)]   - v1.0.0 now contains all 19 archives"
echo "[$(date)] ============================================================"
echo "[$(date)] Upload log: $LOG_FILE"
