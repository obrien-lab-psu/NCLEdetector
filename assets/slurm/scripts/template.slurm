#!/bin/bash
#SBATCH -J job_name
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 24:00:00
#SBATCH --account=epo2_cr_default
#SBATCH --mem=20G
#SBATCH -o slurm/logs/%x-%j.out
#SBATCH -e slurm/logs/%x-%j.err

# Compute Q, G, K (full CG trajectory) and SASA, XP (all-atom last-frames) for traj 420
# CG:  6667 frames, start=0  → full trajectory used for Q/G/K
# AA:  335 frames,  aa_start=0 → back-mapped last-frames used for SASA/XP
#
# Expected runtime: 12–20 h for G (Topoly, 6667 frames, nproc=10)
# Submit: sbatch assets/slurm/run_OP_traj420.slurm

cd /storage/group/epo2/default/ims86/git_repos/EntDetect
source ~/.bashrc
conda activate entdetect

set -euo pipefail

DATASTORE=/scratch/ims86/EntDetect_Datastore
REFSTRUCT=$DATASTORE/user_input/reference_structures

mkdir -p $DATASTORE/outputs/workflow2/OP_demo/{G,Q,K} \
         $DATASTORE/outputs/workflow2/OP_demo_AA/{SASA,XP}

python scripts/run_OP_on_simulation_traj.py \
  --Traj 420 \
  --PSF  $REFSTRUCT/1zmr_model_clean_ca.psf \
  --COR  $REFSTRUCT/1zmr_model_clean_ca.cor \
  --DCD  $DATASTORE/user_input/cg_trajectories/420_prod.dcd \
  --ID   1ZMR \
  --sec_elements $REFSTRUCT/secondary_struc_defs.txt \
  --domain       $REFSTRUCT/domain_def.dat \
  --outdir       $DATASTORE/outputs/workflow2/OP_demo \
  --start    0 \
  --ent_detection_method 2 \
  --nproc    10 \
  --aa_psf   $REFSTRUCT/1zmr_model_clean.pdb \
  --aa_dcd   $DATASTORE/user_input/aa_trajectories/420_prod_aa.dcd \
  --aa_start 0 \
  --xp_pdb   $REFSTRUCT/1zmr_model_clean.pdb
