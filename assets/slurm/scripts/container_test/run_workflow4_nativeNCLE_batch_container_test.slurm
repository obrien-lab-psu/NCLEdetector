#!/bin/bash
#SBATCH -J ct_wf4_native_batch
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 36:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=100G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

set -eo pipefail

REPO=/storage/group/epo2/default/ims86/git_repos/EntDetect
DATASTORE=/scratch/ims86/EntDetect_Datastore
SIF=${SIF:-$REPO/entdetect-v2.0.1.sif}
CFG=scripts/configs/workflow4_nativeNCLE_batch_container_config.json

cd "$REPO"

if [[ ! -f "$SIF" ]]; then
  echo "ERROR: container image not found: $SIF" >&2
  exit 1
fi

mkdir -p "$DATASTORE/outputs/workflow4/nativeNCLE_all/container"
mkdir -p "$DATASTORE/outputs/workflow4/container/logs"

apptainer exec \
  --bind "$DATASTORE:$DATASTORE" \
  --bind "$PWD:$PWD" \
  --pwd "$PWD" \
  "$SIF" \
  python scripts/run_workflow4_nativeNCLE_batch.py --config "$CFG"
