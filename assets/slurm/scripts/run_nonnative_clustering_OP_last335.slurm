#!/bin/bash
#SBATCH -J workflow2_NNclust_last335
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 48:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=128G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

################################################################################
# Workflow 2: Non-native Entanglement Clustering (OP_last335) SLURM Script
#
# Runs scripts/run_nonnative_entanglement_clustering.py using a JSON config
# file consumed by the runner's --config flag.
#
# Expected runtime: hours (depends on number of trajectories and nproc)
# Memory: 128 GB
# CPUs: 8
#
# Usage:
#   sbatch assets/slurm/scripts/run_nonnative_clustering_OP_last335.slurm
#
# Config used by default:
#   scripts/configs/workflow2_nonnative_clustering_config.json
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

CONFIG_PATH=scripts/configs/workflow2_nonnative_clustering_config.json

if [[ ! -f "$CONFIG_PATH" ]]; then
    echo "ERROR: Missing config file: $CONFIG_PATH"
    exit 1
fi

# Parse key paths from config for validation and logging.
readarray -t CFG < <(python - <<'PY' "$CONFIG_PATH"
import json, sys
cfg_path = sys.argv[1]
with open(cfg_path, 'r', encoding='utf-8') as fh:
    cfg = json.load(fh)

print(cfg.get('outdir', ''))
print(cfg.get('trajnum2pklfile_path', ''))
print(cfg.get('traj_dir_prefix', ''))
print(cfg.get('logdir', ''))
print(str(cfg.get('start_frame', 0)))
print(str(cfg.get('end_frame', 9999999)))
print(str(cfg.get('nproc', 1)))
PY
)

OUTDIR=${CFG[0]}
TRAJNUM2FILE=${CFG[1]}
TRAJ_DIR_PREFIX=${CFG[2]}
LOGDIR=${CFG[3]}
START_FRAME=${CFG[4]}
END_FRAME=${CFG[5]}
NPROC=${CFG[6]}

if [[ -z "$OUTDIR" || -z "$TRAJNUM2FILE" || -z "$TRAJ_DIR_PREFIX" ]]; then
    echo "ERROR: Config missing required keys (outdir, trajnum2pklfile_path, traj_dir_prefix): $CONFIG_PATH"
    exit 1
fi

if [[ -z "$LOGDIR" ]]; then
    LOGDIR=$OUTDIR
fi

OP_ROOT=$(dirname "$(dirname "$TRAJNUM2FILE")")

if [[ ! -f "$TRAJNUM2FILE" ]]; then
    echo "ERROR: Missing trajnum2pklfile file from config: $TRAJNUM2FILE"
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
echo "================================================================================"
echo "Running Workflow 2 Non-native Clustering (OP_last335)"
echo "================================================================================"
echo "Start time: $(date)"
echo "Order-parameter root: $OP_ROOT"
echo "Config file: $CONFIG_PATH"
echo "trajnum2file: $TRAJNUM2FILE"
echo "Trajectory prefix: $TRAJ_DIR_PREFIX"
echo "Output dir: $OUTDIR"
echo "start_frame: $START_FRAME"
echo "end_frame: $END_FRAME"
echo "nproc: $NPROC"
echo "================================================================================"
echo ""

python scripts/run_nonnative_entanglement_clustering.py \
    --config "$CONFIG_PATH" \
    2>&1 | tee "$LOGDIR/workflow2_nonnative_clustering_OP_last335.log"

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
    echo "✓ Workflow 2 non-native clustering (OP_last335) completed successfully"
    echo ""
    echo "Output directory: $OUTDIR"
    echo "Mapping file used: $TRAJNUM2FILE"
    echo "Main log: $LOGDIR/workflow2_nonnative_clustering_OP_last335.log"
else
    echo ""
    echo "✗ Workflow 2 non-native clustering (OP_last335) failed with exit code $EXIT_CODE"
    echo "Check log: $LOGDIR/workflow2_nonnative_clustering_OP_last335.log"
fi

exit $EXIT_CODE
