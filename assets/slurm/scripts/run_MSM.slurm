#!/bin/bash
#SBATCH -J workflow2_MSM
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 04:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=64G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

################################################################################
# Workflow 2: MSM Construction SLURM Script
#
# Builds a Markov state model (MSM) from pre-computed Q/G order parameters using
# scripts/run_MSM.py, matching the Workflow 2 tutorial single-script command.
#
# Expected runtime: minutes to ~1 hour (depends on trajectory/frame count)
#
# Memory: 64 GB
# CPUs: 8
#
# Prerequisites:
#   - OP directory with Q/ and G/ subdirectories:
#       /scratch/ims86/EntDetect_Datastore/outputs/workflow2/OP/
#   - Per-trajectory Q files in OP/Q/*.Q
#   - Per-trajectory G files in OP/G/*.G
#   - Conda env: entdetect
#   - Repo path: /storage/group/epo2/default/ims86/git_repos/EntDetect
#
# Monitor progress:
#   tail -f assets/slurm/logs/workflow2_MSM-<JOBID>.out
#
# Usage:
#   sbatch assets/slurm/scripts/run_MSM.slurm
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
OPPATH=$DATASTORE/outputs/workflow2/OP
OUTDIR=$DATASTORE/outputs/workflow2/MSM
LOGDIR=$OUTDIR/logs

if [[ ! -d "$OPPATH" ]]; then
    echo "ERROR: OP directory not found: $OPPATH"
    exit 1
fi

if [[ ! -d "$OPPATH/Q" ]]; then
    echo "ERROR: Missing Q directory: $OPPATH/Q"
    exit 1
fi

if [[ ! -d "$OPPATH/G" ]]; then
    echo "ERROR: Missing G directory: $OPPATH/G"
    exit 1
fi

if [[ -z $(ls -1 "$OPPATH"/Q/*.Q 2>/dev/null | head -1) ]]; then
    echo "ERROR: No Q files found in $OPPATH/Q"
    exit 1
fi
n_q=$(ls -1 "$OPPATH"/Q/*.Q 2>/dev/null | wc -l)
echo "✓ Found $n_q Q files"

if [[ -z $(ls -1 "$OPPATH"/G/*.G 2>/dev/null | head -1) ]]; then
    echo "ERROR: No G files found in $OPPATH/G"
    exit 1
fi
n_g=$(ls -1 "$OPPATH"/G/*.G 2>/dev/null | wc -l)
echo "✓ Found $n_g G files"

if [[ ! -f "scripts/run_MSM.py" ]]; then
    echo "ERROR: Missing script: scripts/run_MSM.py"
    exit 1
fi
echo "✓ Found script: scripts/run_MSM.py"

echo ""
echo "================================================================================"
echo "Running Workflow 2 MSM Construction"
echo "================================================================================"
echo "Start time: $(date)"
echo "OP path: $OPPATH"
echo "Output dir: $OUTDIR"
echo "ID: 1ZMR_prod"
echo "start: 0"
echo "n_large_states: 10"
echo "lagtime: 20"
echo "================================================================================"
echo ""

# Create output and log directories
mkdir -p "$OUTDIR" "$LOGDIR"

python scripts/run_MSM.py \
    --outdir "$OUTDIR" \
    --ID 1ZMR_prod \
    --OPpath "$OPPATH/" \
    --start 0 \
    --n_large_states 10 \
    --lagtime 20 \
    --logdir "$LOGDIR" \
    2>&1 | tee "$LOGDIR/workflow2_MSM.log"

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
    echo "✓ Workflow 2 MSM construction completed successfully!"
    echo ""
    echo "Output files expected in: $OUTDIR"
    echo "  - 1ZMR_prod_MSMmapping.csv"
    echo "  - 1ZMR_prod_meta_set.csv"
    echo "  - 1ZMR_prod_meta_dist.npy"
    echo "  - 1ZMR_prod_StateAndFEplot.png"
    echo ""
    echo "Full stdout log: $LOGDIR/workflow2_MSM.log"
    echo ""
    echo "View results:"
    echo "  tail -n 50 $LOGDIR/workflow2_MSM.log"
else
    echo ""
    echo "✗ Workflow 2 MSM construction failed with exit code $EXIT_CODE"
    echo "Check logs:"
    echo "  tail -n 100 $LOGDIR/workflow2_MSM.log"
fi

exit $EXIT_CODE
