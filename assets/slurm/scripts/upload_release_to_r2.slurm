#!/bin/bash
################################################################################
# Upload NCLEdetector_Datastore release to Cloudflare R2
#
# This job uploads archives and metadata to private R2 bucket.
# Excludes workflow2/OP.tar.gz (to be added later).
#
# Job: 20 CPUs, 32GB memory, up to 8 hours
#
################################################################################

#SBATCH --job-name=upload_release_r2
#SBATCH --partition=basic
#SBATCH --account=read_crch_ims86
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

set -eo pipefail

# Add rclone to PATH
export PATH="/storage/icds/RISE/sw8/rclone/rclone-v1.57.0-linux-amd64/bin:$PATH"

MASTER="/scratch/ims86/NCLEdetector_Datastore"
BUCKET="r2-ncle"
BUCKET_NAME="ncle-public-data"
VERSION="v1.0.0"
DEST="${BUCKET}:${BUCKET_NAME}/NCLEdetector_Datastore/${VERSION}"

UPLOAD_LOG="/tmp/r2_upload_${SLURM_JOB_ID}.log"
NUM_WORKERS=${SLURM_CPUS_PER_TASK:-20}

echo "========================================================================"
echo "Upload Release to Cloudflare R2"
echo "========================================================================"
echo "Source: $MASTER"
echo "Destination: $DEST"
echo "Workers: $NUM_WORKERS"
echo "Exclude: outputs/workflow2/OP.tar.gz"
echo "Start time: $(date)"
echo ""

cd "$MASTER"

# ============================================================================
# Create exclude file
# ============================================================================

EXCLUDE_FILE=$(mktemp)
trap "rm -f $EXCLUDE_FILE" EXIT

cat > "$EXCLUDE_FILE" << 'EXCLUDELIST'
outputs/workflow2/OP.tar.gz
EXCLUDELIST

echo "✓ Exclude file created"
echo ""

# ============================================================================
# Dry-run first (optional, comment out to skip)
# ============================================================================

echo "Running DRY-RUN to preview upload..."
rclone copy "$MASTER" "$DEST" \
    --exclude-from="$EXCLUDE_FILE" \
    --dry-run \
    --transfers=$NUM_WORKERS \
    --checkers=$NUM_WORKERS \
    --log-file="${UPLOAD_LOG}.dryrun" \
    2>&1 | tail -20

echo ""
echo "✓ Dry-run complete (see ${UPLOAD_LOG}.dryrun for details)"
echo ""

# ============================================================================
# Actual upload (with conservative settings)
# ============================================================================

echo "Starting full upload (this may take 30-60 minutes)..."
echo ""

rclone copy "$MASTER" "$DEST" \
    --exclude-from="$EXCLUDE_FILE" \
    --transfers=$NUM_WORKERS \
    --checkers=$NUM_WORKERS \
    --retries=10 \
    --low-level-retries=10 \
    --progress \
    --log-file="$UPLOAD_LOG" \
    --log-level=INFO

EXIT_CODE=$?

echo ""
echo "========================================================================"
echo "Upload Complete"
echo "========================================================================"
echo "End time: $(date)"
echo "Exit code: $EXIT_CODE"
echo ""

if [[ $EXIT_CODE -eq 0 ]]; then
    echo "✓ Upload successful!"
    echo ""
    echo "Next steps:"
    echo "1. Verify upload size:"
    echo "   rclone size $DEST"
    echo ""
    echo "2. Spot-check integrity:"
    echo "   rclone cat $DEST/user_input.tar.gz | md5sum"
    echo "   md5sum $MASTER/user_input.tar.gz"
    echo ""
    echo "3. Once verified, proceed to Phase 3 (CORS + public domain setup)"
else
    echo "⚠ Upload completed with errors (exit code $EXIT_CODE)"
    echo "See $UPLOAD_LOG for details"
    exit $EXIT_CODE
fi

echo ""
echo "Upload log: $UPLOAD_LOG"
