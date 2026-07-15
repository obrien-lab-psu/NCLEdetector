#!/bin/bash
#SBATCH -J ct_op_traj420
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 24:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=20G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

set -eo pipefail

REPO=/storage/group/epo2/default/ims86/git_repos/EntDetect
DATASTORE=/scratch/ims86/EntDetect_Datastore
SIF=${SIF:-$REPO/entdetect-v2.0.1.sif}
CFG=scripts/configs/workflow2_OP_last67_container_config.json

cd "$REPO"

if [[ ! -f "$SIF" ]]; then
  echo "ERROR: container image not found: $SIF" >&2
  exit 1
fi

mkdir -p "$DATASTORE/outputs/workflow2/OP_last67/container/logs"

apptainer exec \
  --bind "$DATASTORE:$DATASTORE" \
  --bind "$PWD:$PWD" \
  --pwd "$PWD" \
  "$SIF" \
  python scripts/run_OP_on_simulation_traj.py --config "$CFG"
