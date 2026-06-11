#!/bin/bash
#SBATCH -J workflow3_consistency
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 06:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=64G
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
# Memory: 64 GB
# CPUs: 8
#
# Prerequisites:
#   - Per-traj SASA files: /scratch/ims86/EntDetect_Datastore/outputs/workflow2/OP_AA/SASA/
#   - Per-traj XP files: /scratch/ims86/EntDetect_Datastore/outputs/workflow2/OP_AA/XP/
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
SASA_DIR=$DATASTORE/outputs/workflow2/OP_AA/SASA
XP_DIR=$DATASTORE/outputs/workflow2/OP_AA/XP

if [[ ! -d "$SASA_DIR" ]] || [[ -z $(ls -1 $SASA_DIR/*.SASA 2>/dev/null | head -1) ]]; then
    echo "ERROR: SASA files not found in $SASA_DIR"
    exit 1
fi
n_sasa=$(ls -1 $SASA_DIR/*.SASA 2>/dev/null | wc -l)
echo "✓ Found $n_sasa SASA files"

if [[ ! -d "$XP_DIR" ]] || [[ -z $(ls -1 $XP_DIR/*.XP 2>/dev/null | head -1) ]]; then
    echo "ERROR: XP files not found in $XP_DIR"
    exit 1
fi
n_xp=$(ls -1 $XP_DIR/*.XP 2>/dev/null | wc -l)
echo "✓ Found $n_xp XP files"

EXP_DATA_DIR="$DATASTORE/user_input/experimental_data"
if [[ ! -d "$EXP_DATA_DIR" ]] || [[ -z $(ls -1 $EXP_DATA_DIR/*.xlsx 2>/dev/null) ]]; then
    echo "ERROR: No experimental Excel files found in $EXP_DATA_DIR"
    exit 1
fi
echo "✓ Experimental data found"

if [[ ! -f "$DATASTORE/outputs/workflow2/MSM/1ZMR_prod_MSMmapping.csv" ]] || [[ ! -f "$DATASTORE/outputs/workflow2/MSM/1ZMR_prod_meta_dist.npy" ]]; then
    echo "ERROR: MSM data not found"
    exit 1
fi
echo "✓ MSM data found"

echo ""
echo "================================================================================"
echo "Running Consistency Test (Collect + Run Mode)"
echo "================================================================================"
echo "Start time: $(date)"
echo "SASA dir: $SASA_DIR"
echo "XP dir: $XP_DIR"
echo "Output dir: $OUTDIR"
echo ""

# Create logs directory
mkdir -p $OUTDIR/logs

python scripts/run_compare_sim2exp.py \
    --sasa_dir "$SASA_DIR" \
    --xp_dir "$XP_DIR" \
    --n_traj 1000 \
    --n_frames 335 \
    --msm_data_file "$DATASTORE/outputs/workflow2/MSM/1ZMR_prod_MSMmapping.csv" \
    --meta_dist_file "$DATASTORE/outputs/workflow2/MSM/1ZMR_prod_meta_dist.npy" \
    --LiPMS_exp_file "$EXP_DATA_DIR/ecPGK_significant_LiPMS_peptide_R1_merged.xlsx" \
    --XLMS_exp_file "$EXP_DATA_DIR/ecPGK_significant_XLMS_peptide_R1_merged.xlsx" \
    --cluster_data_file "$DATASTORE/outputs/workflow2/nonnative_clustering/cluster_data_topoly_linking_number.npz" \
    --OPpath "$DATASTORE/outputs/workflow2/OP_AA/" \
    --AAdcd_dir "$DATASTORE/user_input/aa_trajectories/" \
    --native_AA_pdb "$DATASTORE/user_input/reference_structures/1zmr_model_clean.pdb" \
    --state_idx_list 4 6 8 \
    --prot_len 387 \
    --last_num_frames 335 \
    --rm_traj_list 65 75 155 162 199 231 264 286 296 314 354 417 448 472 473 474 577 579 591 703 704 732 758 812 833 870 876 944 967 \
    --native_state_idx 9 \
    --outdir "$OUTDIR/MassSpec_ConsistencyTest" \
    --ID 1ZMR \
    --start 0 \
    --end -1 \
    --stride 1 \
    --num_perm 1000 \
    --n_boot 100 \
    --lag_frame 20 \
    --nproc 8 \
    2>&1 | tee "$OUTDIR/logs/consistency_test.log"

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
    echo "Output files saved to: $OUTDIR"
    echo "  - consistency_pvalues.npy"
    echo "  - representative_structures.pkl"
    echo "  - consistency_results.png"
    echo "  - logs/consistency_test.log"
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
