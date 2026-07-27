#!/bin/bash
################################################################################
# Generate SHA256SUMS in parallel using SLURM
#
# This job parallelizes checksum computation across multiple cores to speed up
# processing of 59k+ files across ~100GB of archives.
#
# Job: 4 CPUs, 16GB memory, up to 4 hours
#
################################################################################

#SBATCH --job-name=generate_checksums
#SBATCH --partition=basic
#SBATCH --account=read_crch_ims86
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

set -eo pipefail

MASTER="/scratch/ims86/NCLEdetector_Datastore"
CHECKSUMS_FILE="$MASTER/SHA256SUMS"
NUM_WORKERS=${SLURM_CPUS_PER_TASK:-10}

echo "========================================================================"
echo "Parallel SHA256SUMS Generation"
echo "========================================================================"
echo "Master: $MASTER"
echo "Workers: $NUM_WORKERS"
echo "Start time: $(date)"
echo ""

cd "$MASTER"

# ============================================================================
# Generate list of files to checksum (exclude metadata)
# ============================================================================

TEMP_FILELIST=$(mktemp)
trap "rm -f $TEMP_FILELIST" EXIT

echo "Collecting files to checksum..."
find "$MASTER" -type f \
    ! -path '*MANIFEST.tsv' \
    ! -path '*MANIFEST.json' \
    ! -path '*SHA256SUMS' \
    ! -path '*outputs/workflow2/OP.tar.gz' \
    -print0 \
    | LC_ALL=C sort -z \
    | tr '\0' '\n' > "$TEMP_FILELIST"

FILE_COUNT=$(wc -l < "$TEMP_FILELIST")
echo "✓ Found $FILE_COUNT files to process"
echo ""

# ============================================================================
# Compute checksums in parallel
# ============================================================================

echo "Computing checksums in parallel ($NUM_WORKERS workers)..."
echo "This may take 10-30 minutes for 59k+ files..."
echo ""

# Use xargs with parallel workers to compute checksums
cat "$TEMP_FILELIST" \
    | xargs -P "$NUM_WORKERS" -I {} sha256sum {} \
    | sort > "$CHECKSUMS_FILE"

CHECKSUM_COUNT=$(wc -l < "$CHECKSUMS_FILE")
CHECKSUM_SIZE=$(du -h "$CHECKSUMS_FILE" | cut -f1)

echo ""
echo "✓ Computed checksums for $CHECKSUM_COUNT files"
echo "✓ SHA256SUMS size: $CHECKSUM_SIZE"
echo ""

# ============================================================================
# Summary
# ============================================================================

echo "========================================================================"
echo "Checksum Generation Complete"
echo "========================================================================"
echo "End time: $(date)"
echo ""
echo "Sample checksums (first 5):"
head -5 "$CHECKSUMS_FILE"
echo ""
echo "To verify checksums later:"
echo "  cd $MASTER"
echo "  sha256sum -c SHA256SUMS"
echo ""
echo "Log file: assets/slurm/logs/generate_checksums-$SLURM_JOB_ID.out"
