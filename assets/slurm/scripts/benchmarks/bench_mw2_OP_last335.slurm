#!/bin/bash
#SBATCH -J bench_mw2_OP_last335
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 01:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=16G
#SBATCH --array=0-249%50
#SBATCH -o assets/slurm/logs/%x-%A_%a.out
#SBATCH -e assets/slurm/logs/%x-%A_%a.err

################################################################################
# Benchmark - Minimal Workflow 2 (OP): last-10ns window, frames [6600, 6667).
# 5 proteins x 50 trajectories = 250 array tasks (throttled to 50 concurrent).
#   task_id -> protein = task_id / 50 ; traj = task_id % 50
# Config-driven: scripts/configs/benchmarks/<protein>/OP_last335.json
# Per-task overrides: --Traj <t> --DCD <source quench dcd>
################################################################################
set -eo pipefail

REPO=/storage/group/epo2/default/ims86/git_repos/EntDetect
SRC=/storage/group/epo2/default/ims86/git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics
cd "$REPO"
source ~/.bashrc
conda activate entdetect

PROTEINS=(P61175_6XZ7_S P0A6T9_3A7L_A P45748_1HRU_A P0A6L2_1DHP_A P0AD61_4YNG_C)
PIDX=$(( SLURM_ARRAY_TASK_ID / 50 ))
TRAJ=$(( SLURM_ARRAY_TASK_ID % 50 ))
P=${PROTEINS[$PIDX]}
DCD="$SRC/$P/Quenching/${P}_t${TRAJ}_quench.dcd"

echo "=== MW2 OP last335 | protein=$P | traj=$TRAJ | $(date) ==="
echo "DCD: $DCD"
python scripts/run_OP_on_simulation_traj.py \
    --config "scripts/configs/benchmarks/$P/OP_last335.json" \
    --Traj "$TRAJ" \
    --DCD "$DCD"
echo "=== done $P traj $TRAJ | $(date) ==="
