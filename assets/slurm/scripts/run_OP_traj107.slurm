#!/bin/bash
#SBATCH -J OP_traj107
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n 10
#SBATCH -t 72:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=20G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

# Step 1: Compute Q, G, K on the full CG trajectory
# Step 2: Compute SASA and XP on the all-atom back-mapped frames
#
# Expected runtime: 12–20 h for G (Topoly, 6667 frames, nproc=10)
# Submit: sbatch assets/slurm/scripts/run_OP_traj107.slurm

cd /storage/group/epo2/default/ims86/git_repos/EntDetect
source ~/.bashrc
conda activate entdetect

set -euo pipefail

DATASTORE=/scratch/ims86/EntDetect_Datastore
REFSTRUCT=$DATASTORE/user_input/reference_structures

mkdir -p $DATASTORE/outputs/workflow2/OP/G
mkdir -p $DATASTORE/outputs/workflow2/OP/Q
mkdir -p $DATASTORE/outputs/workflow2/OP/K
mkdir -p $DATASTORE/outputs/workflow2/OP/logs
mkdir -p $DATASTORE/outputs/workflow2/OP_AA/SASA
mkdir -p $DATASTORE/outputs/workflow2/OP_AA/XP
mkdir -p $DATASTORE/outputs/workflow2/OP_AA/logs

# --- CG: Q, G, K ---
python scripts/run_OP_on_simulation_traj.py   --Traj 107   --PSF  $REFSTRUCT/1zmr_model_clean_ca.psf   --COR  $REFSTRUCT/1zmr_model_clean_ca.cor   --DCD  $DATASTORE/user_input/cg_trajectories/107_prod.dcd   --resolution cg   --contacts calpha   --ID   1ZMR   --sec_elements $REFSTRUCT/secondary_struc_defs.txt   --domain       $REFSTRUCT/domain_def.dat   --outdir       $DATASTORE/outputs/workflow2/OP   --logdir       $DATASTORE/outputs/workflow2/OP/logs   --start        0   --ent_detection_method 1   --nproc        10   --ops Q G K   --no_topoly   --chunk_frames 100 --chunk_suffix _chunk 

# --- AA: SASA, XP ---
python scripts/run_OP_on_simulation_traj.py   --Traj   107   --ID     1ZMR   --PSF    $REFSTRUCT/1zmr_model_clean.pdb   --COR    $REFSTRUCT/1zmr_model_clean.pdb   --DCD    $DATASTORE/user_input/aa_trajectories/107_prod_aa.dcd   --resolution aa   --contacts calpha   --sec_elements $REFSTRUCT/secondary_struc_defs.txt   --domain       $REFSTRUCT/domain_def.dat   --outdir $DATASTORE/outputs/workflow2/OP_AA   --logdir $DATASTORE/outputs/workflow2/OP_AA/logs   --start  0   --nproc        10   --xp_pdb $REFSTRUCT/1zmr_model_clean.pdb   --ops SASA XP
