#!/bin/bash
#SBATCH -J bench_mw4_MSM
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 04:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=64G
#SBATCH --array=0-4
#SBATCH -o assets/slurm/logs/%x-%A_%a.out
#SBATCH -e assets/slurm/logs/%x-%A_%a.err

################################################################################
# Benchmark - Minimal Workflow 4: MSM construction.
# One array task per protein. Requires MW2 OP_full to have completed (uses the
# full-window Q/G outputs). No mirror trajectories are excluded (rm_traj_list=[]).
# Config-driven: scripts/configs/benchmarks/<protein>/MSM.json
################################################################################
set -eo pipefail

REPO=/storage/group/epo2/default/ims86/git_repos/EntDetect
cd "$REPO"
source ~/.bashrc
conda activate entdetect

PROTEINS=(P61175_6XZ7_S P0A6T9_3A7L_A P45748_1HRU_A P0A6L2_1DHP_A P0AD61_4YNG_C)
P=${PROTEINS[$SLURM_ARRAY_TASK_ID]}

echo "=== MW4 MSM | protein=$P | $(date) ==="
python scripts/run_MSM.py --config "scripts/configs/benchmarks/$P/MSM.json"
echo "=== done $P | $(date) ==="
