#!/bin/bash
#SBATCH -J create_workflow2_FoldingPathway_random
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 00:30:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=16G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

set -eo pipefail

MASTER="/scratch/ims86/NCLEdetector_Datastore"
WORKFLOW="workflow2"
ANALYSIS="FoldingPathway_random"
ARCHIVE_DIR="$MASTER/outputs/$WORKFLOW"
SOURCE_DIR="$ARCHIVE_DIR/$ANALYSIS"
ARCHIVE="$ARCHIVE_DIR/$ANALYSIS.tar.gz"
PIGZ_THREADS=${SLURM_NPROCS:-4}

echo "Creating $WORKFLOW/$ANALYSIS.tar.gz"
echo "Start: $(date)"

if [[ -f "$ARCHIVE" ]]; then
    echo "Archive already exists: $ARCHIVE"
    exit 0
fi

tar -c --use-compress-program="pigz -p $PIGZ_THREADS" \
    -f "$ARCHIVE" \
    -C "$ARCHIVE_DIR" "$ANALYSIS"

# Validate
gzip -t "$ARCHIVE" && tar -tzf "$ARCHIVE" >/dev/null && echo "✓ Validated"

echo "Size: $(du -h "$ARCHIVE" | cut -f1)"
echo "End: $(date)"
