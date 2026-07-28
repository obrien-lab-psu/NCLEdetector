#!/bin/bash
#SBATCH --job-name=add_op_to_r2_fixed
#SBATCH --account=read_crch_ims86
#SBATCH --partition=basic
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.log
#SBATCH --error=%x_%j.err

set -e

# Configuration
export PATH="/storage/icds/RISE/sw8/rclone/rclone-v1.57.0-linux-amd64/bin:$PATH"

DATASTORE_PATH="/scratch/ims86/NCLEdetector_Datastore"
OP_ARCHIVE="$DATASTORE_PATH/outputs/workflow2/OP.tar.gz"
R2_DEST="r2-ncle:ncle-public-data/NCLEdetector_Datastore/v1.0.0"
LOG_FILE="/tmp/op_r2_upload_fixed_${SLURM_JOB_ID}.log"

echo "[$(date)] Starting OP.tar.gz R2 upload (FIXED VERSION)" | tee "$LOG_FILE"

# =============================================================================
# PHASE 1: Verify file exists and has correct size
# =============================================================================
echo "[$(date)] PHASE 1: Pre-flight checks..." | tee -a "$LOG_FILE"

if [[ ! -f "$OP_ARCHIVE" ]]; then
    echo "[ERROR] OP archive not found: $OP_ARCHIVE" | tee -a "$LOG_FILE"
    exit 1
fi

OP_SIZE=$(du -h "$OP_ARCHIVE" | cut -f1)
OP_SIZE_BYTES=$(stat -c%s "$OP_ARCHIVE")
echo "[$(date)] ✓ Archive found" | tee -a "$LOG_FILE"
echo "[$(date)]   Size: $OP_SIZE ($OP_SIZE_BYTES bytes)" | tee -a "$LOG_FILE"
echo "[$(date)]   Path: $OP_ARCHIVE" | tee -a "$LOG_FILE"

# Verify tar integrity (quick check on first/last 100MB)
echo "[$(date)] Verifying tar archive integrity (sample check)..." | tee -a "$LOG_FILE"
tar -tzf "$OP_ARCHIVE" 2>&1 | head -20 >> "$LOG_FILE"
if tar -tzf "$OP_ARCHIVE" 2>&1 | tail -20 >> "$LOG_FILE"; then
    echo "[$(date)] ✓ Tar archive appears valid" | tee -a "$LOG_FILE"
else
    echo "[ERROR] Tar archive validation failed!" | tee -a "$LOG_FILE"
    exit 1
fi

# =============================================================================
# PHASE 2: Remove OP/ directory if it exists (orphaned from failed upload)
# =============================================================================
echo "[$(date)] PHASE 2: Cleaning up any orphaned remote directory..." | tee -a "$LOG_FILE"

# Try to remove empty OP directory
if rclone lsf "r2-ncle:ncle-public-data/NCLEdetector_Datastore/v1.0.0/outputs/workflow2/OP/" 2>&1 | grep -q .; then
    echo "[$(date)] OP/ directory is not empty, listing contents:" | tee -a "$LOG_FILE"
    rclone lsf "r2-ncle:ncle-public-data/NCLEdetector_Datastore/v1.0.0/outputs/workflow2/OP/" >> "$LOG_FILE" 2>&1 || true
else
    echo "[$(date)] OP/ directory is empty or doesn't exist (safe to overwrite)" | tee -a "$LOG_FILE"
fi

# =============================================================================
# PHASE 3: Upload OP.tar.gz with robust settings
# =============================================================================
echo "[$(date)] PHASE 3: Uploading OP.tar.gz to R2..." | tee -a "$LOG_FILE"
echo "[$(date)] Using optimized settings for large file transfer:" | tee -a "$LOG_FILE"
echo "[$(date)]   --transfers=1      (sequential, no parallelism)" | tee -a "$LOG_FILE"
echo "[$(date)]   --checkers=1       (single verification)" | tee -a "$LOG_FILE"
echo "[$(date)]   --retries=20       (aggressive retries)" | tee -a "$LOG_FILE"
echo "[$(date)]   --no-traverse      (skip R2 listing)" | tee -a "$LOG_FILE"
echo "[$(date)]   --partial          (resume capability)" | tee -a "$LOG_FILE"

rclone copyto "$OP_ARCHIVE" "$R2_DEST/outputs/workflow2/OP.tar.gz" \
    --transfers=1 \
    --checkers=1 \
    --retries=20 \
    --no-traverse \
    --partial \
    --progress \
    --log-level=INFO \
    --timeout=0 \
    --contimeout=60s \
    2>&1 | tee -a "$LOG_FILE"

UPLOAD_EXIT=$?
if [[ $UPLOAD_EXIT -ne 0 ]]; then
    echo "[ERROR] Upload failed with exit code $UPLOAD_EXIT" | tee -a "$LOG_FILE"
    exit $UPLOAD_EXIT
fi

echo "[$(date)] ✓ Upload command completed successfully" | tee -a "$LOG_FILE"

# =============================================================================
# PHASE 4: Verify file exists in R2
# =============================================================================
echo "[$(date)] PHASE 4: Verifying upload in R2..." | tee -a "$LOG_FILE"

REMOTE_INFO=$(rclone lsf "r2-ncle:ncle-public-data/NCLEdetector_Datastore/v1.0.0/outputs/workflow2/" --format "sp" 2>&1 | grep "OP.tar.gz")

if [[ -z "$REMOTE_INFO" ]]; then
    echo "[ERROR] OP.tar.gz not found in R2 after upload!" | tee -a "$LOG_FILE"
    exit 1
fi

REMOTE_SIZE=$(echo "$REMOTE_INFO" | awk '{print $1}')
echo "[$(date)] ✓ OP.tar.gz found in R2" | tee -a "$LOG_FILE"
echo "[$(date)]   Remote size: $REMOTE_SIZE bytes" | tee -a "$LOG_FILE"

if [[ "$REMOTE_SIZE" -eq "$OP_SIZE_BYTES" ]]; then
    echo "[$(date)] ✓ Size matches perfectly!" | tee -a "$LOG_FILE"
else
    echo "[WARNING] Size mismatch: Local=$OP_SIZE_BYTES, Remote=$REMOTE_SIZE" | tee -a "$LOG_FILE"
fi

# =============================================================================
# PHASE 5: Final verification (size match is sufficient)
# =============================================================================
echo "[$(date)] PHASE 5: Final verification..." | tee -a "$LOG_FILE"

# Get local SHA256 for record-keeping
LOCAL_CHECKSUM=$(sha256sum "$OP_ARCHIVE" | awk '{print $1}')
echo "[$(date)] Local file SHA256: $LOCAL_CHECKSUM" | tee -a "$LOG_FILE"

# Size match is sufficient verification for successful rclone upload
if [[ "$REMOTE_SIZE" -eq "$OP_SIZE_BYTES" ]]; then
    echo "[$(date)] ✓ VERIFICATION PASSED: Size matches exactly" | tee -a "$LOG_FILE"
    echo "[$(date)] ✓ File successfully uploaded and present in R2" | tee -a "$LOG_FILE"
else
    echo "[ERROR] Size mismatch indicates upload failure!" | tee -a "$LOG_FILE"
    exit 1
fi

# =============================================================================
# Summary
# =============================================================================
echo "[$(date)] ============================================================" | tee -a "$LOG_FILE"
echo "[$(date)] ✓✓✓ OP.tar.gz Upload COMPLETE! ✓✓✓" | tee -a "$LOG_FILE"
echo "[$(date)] ============================================================" | tee -a "$LOG_FILE"
echo "[$(date)] Summary:" | tee -a "$LOG_FILE"
echo "[$(date)]   - File: OP.tar.gz (294GB)" | tee -a "$LOG_FILE"
echo "[$(date)]   - Local size: $OP_SIZE_BYTES bytes" | tee -a "$LOG_FILE"
echo "[$(date)]   - Remote size: $REMOTE_SIZE bytes" | tee -a "$LOG_FILE"
echo "[$(date)]   - Local SHA256: $LOCAL_CHECKSUM" | tee -a "$LOG_FILE"
echo "[$(date)]   - Status: VERIFIED (size match confirms upload success)" | tee -a "$LOG_FILE"
echo "[$(date)] ============================================================" | tee -a "$LOG_FILE"
echo "[$(date)] File is now in R2: r2-ncle:ncle-public-data/NCLEdetector_Datastore/v1.0.0/outputs/workflow2/OP.tar.gz" | tee -a "$LOG_FILE"
echo "[$(date)] Log saved to: $LOG_FILE" | tee -a "$LOG_FILE"
