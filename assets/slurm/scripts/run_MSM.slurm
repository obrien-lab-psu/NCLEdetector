#!/bin/bash
#SBATCH -J MSM
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=50G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err


cd /storage/group/epo2/default/ims86/git_repos/EntDetect
source ~/.bashrc
conda activate entdetect

set -euo pipefail

DATASTORE=/scratch/ims86/EntDetect_Datastore

mkdir -p $DATASTORE/outputs/workflow2/MSM
mkdir -p $DATASTORE/outputs/workflow2/MSM/logs

# run MSM
python scripts/run_MSM.py \
    --outdir         $DATASTORE/outputs/workflow2/MSM \
    --ID             1ZMR_prod \
    --OPpath         $DATASTORE/outputs/workflow2/OP_demo/ \
    --start          0 \
    --n_large_states 10 \
    --lagtime        20 \
    --logdir         $DATASTORE/outputs/workflow2/MSM/logs
