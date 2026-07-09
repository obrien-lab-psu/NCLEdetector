#!/bin/bash
#SBATCH -J workflow2_NNclust_last67
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 12:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=128G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

################################################################################
# Workflow 2: Non-native Entanglement Clustering (OP_last67) SLURM Script
#
# Runs scripts/run_nonnative_entanglement_clustering.py using the OP_last67
# order-parameter output folder and regenerates trajnum2file.txt from:
#   $DATASTORE/outputs/workflow2/OP_last67/G/Combined_NCLE/*.pkl
#
# Expected runtime: hours (depends on number of trajectories and nproc)
# Memory: 128 GB
# CPUs: 8
#
# Usage:
#   sbatch assets/slurm/scripts/run_nonnative_clustering_OP_last67.slurm
#
################################################################################

set -eo pipefail

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

cd /storage/group/epo2/default/ims86/git_repos/EntDetect

echo "Activating conda environment..."
source ~/.bashrc
conda activate entdetect

echo ""
echo "Verifying prerequisites..."

DATASTORE=/scratch/ims86/EntDetect_Datastore
OP_ROOT=$DATASTORE/outputs/workflow2/OP_last67
PKL_DIR=$OP_ROOT/G/Combined_NCLE
TRAJNUM2FILE=$OP_ROOT/trajnum2file.txt
OUTDIR=$DATASTORE/outputs/workflow2/nonnative_clustering_last67
LOGDIR=$OUTDIR/logs
TRAJ_DIR_PREFIX=$DATASTORE/user_input/cg_trajectories

# For OP_last67 we typically process all frames available in these pkl files.
START_FRAME=0
END_FRAME=9999999
NPROC=8

if [[ ! -d "$PKL_DIR" ]]; then
    echo "ERROR: Missing Combined_NCLE directory: $PKL_DIR"
    exit 1
fi

if [[ ! -f scripts/run_nonnative_entanglement_clustering.py ]]; then
    echo "ERROR: Missing script: scripts/run_nonnative_entanglement_clustering.py"
    exit 1
fi

if [[ ! -d "$TRAJ_DIR_PREFIX" ]]; then
    echo "ERROR: Missing trajectory directory: $TRAJ_DIR_PREFIX"
    exit 1
fi

mkdir -p "$OUTDIR" "$LOGDIR"

echo ""
echo "Rebuilding trajnum2file mapping from: $PKL_DIR"

if ! find "$PKL_DIR" -maxdepth 1 -type f -name '*.pkl' -print -quit | grep -q .; then
    echo "ERROR: No .pkl files found in $PKL_DIR"
    exit 1
fi

{
    echo "trajnum,pklfile"
    find "$PKL_DIR" -maxdepth 1 -type f -name '*.pkl' -print \
        | awk '{n=$0; sub(/^.*traj/, "", n); sub(/[^0-9].*$/, "", n); if (n ~ /^[0-9]+$/) print n "," $0}' \
        | sort -t, -k1,1n \
        | awk -F, '!seen[$1]++'
} > "$TRAJNUM2FILE"

if [[ $(wc -l < "$TRAJNUM2FILE") -le 1 ]]; then
    echo "ERROR: Failed to generate trajectory mapping from $PKL_DIR"
    exit 1
fi

N_MAP_ROWS=$(($(wc -l < "$TRAJNUM2FILE") - 1))
echo "✓ Wrote mapping file: $TRAJNUM2FILE"
echo "✓ Trajectories mapped: $N_MAP_ROWS"

echo ""
echo "================================================================================"
echo "Running Workflow 2 Non-native Clustering (OP_last67)"
echo "================================================================================"
echo "Start time: $(date)"
echo "Order-parameter root: $OP_ROOT"
echo "Combined_NCLE source: $PKL_DIR"
echo "trajnum2file: $TRAJNUM2FILE"
echo "Trajectory prefix: $TRAJ_DIR_PREFIX"
echo "Output dir: $OUTDIR"
echo "start_frame: $START_FRAME"
echo "end_frame: $END_FRAME"
echo "nproc: $NPROC"
echo "================================================================================"
echo ""

python scripts/run_nonnative_entanglement_clustering.py \
    --outdir "$OUTDIR" \
    --trajnum2pklfile_path "$TRAJNUM2FILE" \
    --traj_dir_prefix "$TRAJ_DIR_PREFIX" \
    --start_frame "$START_FRAME" \
    --end_frame "$END_FRAME" \
    --nproc "$NPROC" \
    --log_level INFO \
    --logdir "$LOGDIR" \
    2>&1 | tee "$LOGDIR/workflow2_nonnative_clustering_OP_last67.log"

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
    echo "✓ Workflow 2 non-native clustering (OP_last67) completed successfully"
    echo ""
    echo "Output directory: $OUTDIR"
    echo "Mapping file used: $TRAJNUM2FILE"
    echo "Main log: $LOGDIR/workflow2_nonnative_clustering_OP_last67.log"
else
    echo ""
    echo "✗ Workflow 2 non-native clustering (OP_last67) failed with exit code $EXIT_CODE"
    echo "Check log: $LOGDIR/workflow2_nonnative_clustering_OP_last67.log"
fi

exit $EXIT_CODE
