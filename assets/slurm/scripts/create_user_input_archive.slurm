#!/bin/bash
#SBATCH -J create_user_input_archive
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 01:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=16G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

################################################################################
# Create user_input.tar.gz for public release
#
# Compresses /scratch/ims86/NCLEdetector_Datastore/user_input into a single
# archive for distribution.
#
# Expected runtime: 10-15 minutes
# Memory: 16 GB
# CPUs: 4 (pigz parallel compression)
#
# Usage:
#   sbatch assets/slurm/scripts/create_user_input_archive.slurm
#
################################################################################

set -eo pipefail

echo "================================================================================"
echo "Creating user_input.tar.gz"
echo "================================================================================"
echo "Start time: $(date)"
echo ""

MASTER="/scratch/ims86/NCLEdetector_Datastore"
ARCHIVE="$MASTER/user_input.tar.gz"
PIGZ_THREADS=${SLURM_NPROCS:-4}

if [[ -f "$ARCHIVE" ]]; then
    echo "Archive already exists: $ARCHIVE"
    exit 0
fi

echo "Compressing $MASTER/user_input..."
tar -c --use-compress-program="pigz -p $PIGZ_THREADS" \
    -f "$ARCHIVE" \
    -C "$MASTER" "user_input"

echo ""
echo "Validating archive..."
gzip -t "$ARCHIVE" && tar -tzf "$ARCHIVE" >/dev/null && echo "✓ Validation passed"

echo ""
echo "Archive size: $(du -h "$ARCHIVE" | cut -f1)"
echo "End time: $(date)"
