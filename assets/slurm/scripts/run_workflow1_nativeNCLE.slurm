#!/bin/bash
#SBATCH -J workflow1_nativeNCLE
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=16G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

################################################################################
# Workflow 1: Native NCLE (Single Structure) SLURM Script
#
# Runs native non-covalent lasso entanglement (NCLE) analysis for one reference
# structure using scripts/run_nativeNCLE.py, matching the Workflow 1 tutorial.
#
# Expected runtime: ~2-20 minutes total
#   - Native entanglement detection: ~1-10 minutes
#   - High-quality filtering: < 1 minute
#   - Clustering: < 1 minute
#   - Feature generation: ~1-5 minutes
#
# Memory: 16 GB
# CPUs: 1
#
# Prerequisites:
#   - Structure file:
#       /scratch/ims86/EntDetect_Datastore/user_input/reference_structures/1zmr_model_clean.pdb
#   - Conda env: entdetect
#   - Repo path: /storage/group/epo2/default/ims86/git_repos/EntDetect
#
# Monitor progress:
#   tail -f assets/slurm/logs/workflow1_nativeNCLE-<JOBID>.out
#
# Usage:
#   sbatch assets/slurm/scripts/run_workflow1_nativeNCLE.slurm
#
################################################################################

# Set up error handling (exit on error, but allow unbound vars in external sourced scripts)
set -eo pipefail

# Print job info
echo "================================================================================"
echo "SLURM Job Information"
echo "================================================================================"
echo "Job ID: ${SLURM_JOB_ID:-unknown}"
echo "Job Name: ${SLURM_JOB_NAME:-unknown}"
echo "Partition: ${SLURM_PARTITION:-unknown}"
echo "Nodes: ${SLURM_NNODES:-1}"
echo "CPUs: ${SLURM_NPROCS:-1}"
echo "Memory: ${SLURM_MEM_PER_NODE:-unknown} MB"
echo "Time limit: ${SLURM_TIMELIMIT:-unknown} minutes"
echo "Start time: $(date)"
echo "================================================================================"
echo ""

# Navigate to project directory
cd /storage/group/epo2/default/ims86/git_repos/EntDetect

# Activate conda environment
echo "Activating conda environment..."
source ~/.bashrc
conda activate entdetect

# Verify required files exist
echo ""
echo "Verifying prerequisites..."

DATASTORE=/scratch/ims86/EntDetect_Datastore
OUTDIR=$DATASTORE/outputs/workflow1
STRUCT_FILE=$DATASTORE/user_input/reference_structures/1zmr_model_clean.pdb

if [[ ! -f "$STRUCT_FILE" ]]; then
    echo "ERROR: Structure file not found: $STRUCT_FILE"
    exit 1
fi
echo "✓ Found structure file: $STRUCT_FILE"

if [[ ! -f "scripts/run_nativeNCLE.py" ]]; then
    echo "ERROR: Missing script: scripts/run_nativeNCLE.py"
    exit 1
fi
echo "✓ Found script: scripts/run_nativeNCLE.py"

echo ""
echo "================================================================================"
echo "Running Workflow 1 Native NCLE Analysis"
echo "================================================================================"
echo "Start time: $(date)"
echo "Structure: $STRUCT_FILE"
echo "Output dir: $OUTDIR"
echo "ID: 1ZMR"
echo "Chain: A"
echo "Organism: Ecoli"
echo "Accession: P00558"
echo "Model: EXP"
echo "================================================================================"
echo ""

# Create output and log directories
mkdir -p "$OUTDIR" "$OUTDIR/logs"

python scripts/run_nativeNCLE.py \
    --struct "$STRUCT_FILE" \
    --outdir "$OUTDIR" \
    --ID 1ZMR \
    --chain A \
    --organism Ecoli \
    --Accession P00558 \
    --model EXP \
    2>&1 | tee "$OUTDIR/logs/workflow1_nativeNCLE.log"

# Check exit status
EXIT_CODE=${PIPESTATUS[0]}

echo ""
echo "================================================================================"
echo "Job Completion"
echo "================================================================================"
echo "Exit code: $EXIT_CODE"
echo "End time: $(date)"
echo "================================================================================"

if [[ $EXIT_CODE -eq 0 ]]; then
    echo ""
    echo "✓ Workflow 1 native NCLE analysis completed successfully!"
    echo ""
    echo "Output files are expected under: $OUTDIR"
    echo "  - NCLE/"
    echo "  - HQ_NCLE/"
    echo "  - clustered_HQ_NCLE/"
    echo "  - clustered_HQ_NCLE_features/"
    echo ""
    echo "Full stdout log: $OUTDIR/logs/workflow1_nativeNCLE.log"
    echo ""
    echo "View results:"
    echo "  tail -n 50 $OUTDIR/logs/workflow1_nativeNCLE.log"
else
    echo ""
    echo "✗ Workflow 1 native NCLE analysis failed with exit code $EXIT_CODE"
    echo "Check logs:"
    echo "  tail -n 100 $OUTDIR/logs/workflow1_nativeNCLE.log"
fi

exit $EXIT_CODE
