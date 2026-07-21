#!/bin/bash
#SBATCH -J workflow3_consistency
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 36:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=100G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

################################################################################
# Workflow 3: Consistency Test (Collect + Run) SLURM Script
# 
# Collects per-trajectory SASA/XP files into arrays and runs the MassSpec
# consistency test comparing experimental LiP-MS/XL-MS data to the ensemble.
#
# Expected runtime: 4-5 hours total
#   - SASA collection: ~1 minute
#   - Jwalk collection: ~3-4 hours
#   - Consistency test: ~1-3 hours
#
# Memory: 100 GB
# CPUs: 10
#
# Prerequisites:
#   - Per-traj SASA files: /scratch/ims86/EntDetect_Datastore/outputs/workflow2/OP_AA_last335/SASA/
#   - Per-traj XP files: /scratch/ims86/EntDetect_Datastore/outputs/workflow2/OP_AA_last335/XP/
#   - Experimental Excel files: /scratch/ims86/EntDetect_Datastore/exp_data/
#   - MSM data: /scratch/ims86/EntDetect_Datastore/outputs/workflow2/MSM/
#
# Monitor progress:
#   tail -f assets/slurm/logs/workflow3_consistency-<JOBID>.out
#
# Usage:
#   sbatch assets/slurm/scripts/run_workflow3_consistency.slurm
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
OUTDIR=$DATASTORE/outputs/workflow3


echo ""
echo "================================================================================"
echo "Running Consistency Test (Collect + Run Mode)"
echo "================================================================================"
echo "Start time: $(date)"
echo ""

# Create logs directory
mkdir -p $OUTDIR/logs

python scripts/run_compare_sim2exp.py --config /storage/group/epo2/default/ims86/git_repos/EntDetect/scripts/configs/workflow3_consistency_config.json --restart

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
    echo "✓ Consistency test completed successfully!"
    echo ""
    echo "Output files saved to: $OUTDIR/MassSpec_ConsistencyTest_last335/"
    echo "  - SASA.npy"
    echo "  - native_M_data.npz"
    echo "  - CollectOP.log"
    echo "  - consist_signal_struct_data.npz"
    echo "  - Consistent_structures_v8.xlsx"
    echo "  - LiPMS_XLMS_consist_data_v9.npz"
    echo "  - LiPMS_XLMS_consist_pvalues_metastates_v11_down_sample_lag20.xlsx"
    echo "  - 1ZMR.log"
    echo ""
    echo "Full stdout log: $OUTDIR/logs/consistency_test.log"
    echo ""
    echo "View results:"
    echo "  tail -n 50 $OUTDIR/logs/consistency_test.log"
else
    echo ""
    echo "✗ Consistency test failed with exit code $EXIT_CODE"
    echo "Check logs:"
    echo "  tail -n 100 $OUTDIR/logs/consistency_test.log"
fi

exit $EXIT_CODE
