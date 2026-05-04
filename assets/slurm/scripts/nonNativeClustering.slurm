#!/bin/bash
#SBATCH -J NonNativeClustering
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 72:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=100G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err


cd /storage/group/epo2/default/ims86/git_repos/EntDetect
source ~/.bashrc
conda activate entdetect

set -euo pipefail

DATASTORE=/scratch/ims86/EntDetect_Datastore
REFSTRUCT=$DATASTORE/user_input/reference_structures

mkdir -p $DATASTORE/outputs/workflow2/nonnative_clustering
mkdir -p $DATASTORE/outputs/workflow2/nonnative_clustering/logs

# run nonnative clustering
python scripts/run_nonnative_entanglement_clustering.py \
    --outdir               $DATASTORE/outputs/workflow2/nonnative_clustering \
    --pkl_file_path        $DATASTORE/outputs/workflow2/OP_demo/G/Combined_GE/ \
    --trajnum2pklfile_path $DATASTORE/user_input/metadata/trajnum2file.txt \
    --traj_dir_prefix      $DATASTORE/user_input/cg_trajectories \
    --start_frame          0 \
    --logdir               $DATASTORE/outputs/workflow2/nonnative_clustering/logs \
    --nproc               10