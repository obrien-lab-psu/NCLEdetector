#!/bin/bash
#SBATCH -J bench_mw1_nativeNCLE
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=16G
#SBATCH --array=0-4
#SBATCH -o assets/slurm/logs/%x-%A_%a.out
#SBATCH -e assets/slurm/logs/%x-%A_%a.err

################################################################################
# Benchmark - Minimal Workflow 1: native NCLE on CG reference structures.
# One array task per selected protein (5 total).
# Config-driven: scripts/configs/benchmarks/<protein>/nativeNCLE.json
################################################################################
set -eo pipefail

REPO=/storage/group/epo2/default/ims86/git_repos/EntDetect
cd "$REPO"
source ~/.bashrc
conda activate entdetect

PROTEINS=(P61175_6XZ7_S P0A6T9_3A7L_A P45748_1HRU_A P0A6L2_1DHP_A P0AD61_4YNG_C)
P=${PROTEINS[$SLURM_ARRAY_TASK_ID]}

echo "=== MW1 native NCLE | protein=$P | $(date) ==="
python scripts/run_nativeNCLE.py --config "scripts/configs/benchmarks/$P/nativeNCLE.json"
echo "=== done $P | $(date) ==="
