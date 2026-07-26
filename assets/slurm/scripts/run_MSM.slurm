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
#       /scratch/ims86/NCLEdetector_Datastore/outputs/workflow2/OP/
#   - Per-trajectory Q files in OP/Q/*.Q
#   - Per-trajectory G files in OP/G/*.G
#   - Conda env: entdetect
#   - Repo path: /storage/group/epo2/default/ims86/git_repos/NCLEdetector
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
cd /storage/group/epo2/default/ims86/git_repos/NCLEdetector

# Activate conda environment
echo "Activating conda environment..."
source ~/.bashrc
conda activate ncledetector

python scripts/run_MSM.py --config scripts/configs/workflow2_MSM_config.json 

exit $EXIT_CODE
