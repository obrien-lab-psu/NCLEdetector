import os,sys
import glob

####################################################################################################################################
## TEMPLATES #######################################################################################################################
####################################################################################################################################
run_OP_on_simulation_traj_template_slurm = """#!/bin/bash
#SBATCH -J {job_name}
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n {nproc}
#SBATCH -t 72:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=20G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

# Step 1: Compute Q, G, K on the full CG trajectory
# Step 2: Compute SASA and XP on the all-atom back-mapped frames
#
# Expected runtime: 12–20 h for G (Topoly, 6667 frames, nproc=10)
# Submit: sbatch assets/slurm/scripts/run_OP_traj{traj_num}.slurm

cd /storage/group/epo2/default/ims86/git_repos/NCLEdetector
source ~/.bashrc
conda activate ncledetector

set -euo pipefail

DATASTORE=/scratch/ims86/NCLEdetector_Datastore
REFSTRUCT=$DATASTORE/user_input/reference_structures

mkdir -p $DATASTORE/outputs/workflow2/OP/G
mkdir -p $DATASTORE/outputs/workflow2/OP/Q
mkdir -p $DATASTORE/outputs/workflow2/OP/K
mkdir -p $DATASTORE/outputs/workflow2/OP/logs
mkdir -p $DATASTORE/outputs/workflow2/OP_AA/SASA
mkdir -p $DATASTORE/outputs/workflow2/OP_AA/XP
mkdir -p $DATASTORE/outputs/workflow2/OP_AA/logs

# --- CG: Q, G, K ---
python scripts/run_OP_on_simulation_traj.py \
  --Traj {traj_num} \
  --PSF  $REFSTRUCT/1zmr_model_clean_ca.psf \
  --COR  $REFSTRUCT/1zmr_model_clean_ca.cor \
  --DCD  $DATASTORE/user_input/cg_trajectories/{traj_num}_prod.dcd \
  --resolution cg \
  --contacts calpha \
  --ID   1ZMR \
  --sec_elements $REFSTRUCT/secondary_struc_defs.txt \
  --domain       $REFSTRUCT/domain_def.dat \
  --outdir       $DATASTORE/outputs/workflow2/OP \
  --logdir       $DATASTORE/outputs/workflow2/OP/logs \
  --start        0 \
  --ent_detection_method 1 \
  --nproc        {nproc} \
  --ops Q G K \
  --no_topoly \
  --chunk_frames 100 --chunk_suffix _chunk 

# --- AA: SASA, XP ---
python scripts/run_OP_on_simulation_traj.py \
  --Traj   {traj_num} \
  --ID     1ZMR \
  --PSF    $REFSTRUCT/1zmr_model_clean.pdb \
  --COR    $REFSTRUCT/1zmr_model_clean.pdb \
  --DCD    $DATASTORE/user_input/aa_trajectories/{traj_num}_prod_aa.dcd \
  --resolution aa \
  --contacts calpha \
  --sec_elements $REFSTRUCT/secondary_struc_defs.txt \
  --domain       $REFSTRUCT/domain_def.dat \
  --outdir $DATASTORE/outputs/workflow2/OP_AA \
  --logdir $DATASTORE/outputs/workflow2/OP_AA/logs \
  --start  0 \
  --nproc        {nproc} \
  --xp_pdb $REFSTRUCT/1zmr_model_clean.pdb \
  --ops SASA XP
"""
####################################################################################################################################

####################################################################################################################################
run_OP_on_simulation_traj_template_slurm_v2 = """#!/bin/bash
#SBATCH -J {job_name}
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n {nproc}
#SBATCH -t 72:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=20G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

# Step 2: Compute SASA and XP on the all-atom back-mapped frames
#
# Expected runtime: 12–20 h for G (Topoly, 6667 frames, nproc=10)
# Submit: sbatch assets/slurm/scripts/run_OP_traj{traj_num}.slurm

cd /storage/group/epo2/default/ims86/git_repos/NCLEdetector
source ~/.bashrc
conda activate ncledetector

set -euo pipefail

DATASTORE=/scratch/ims86/NCLEdetector_Datastore
REFSTRUCT=$DATASTORE/user_input/reference_structures

mkdir -p $DATASTORE/outputs/workflow2/OP_AA/SASA
mkdir -p $DATASTORE/outputs/workflow2/OP_AA/XP
mkdir -p $DATASTORE/outputs/workflow2/OP_AA/logs

# --- AA: SASA, XP ---
python scripts/run_OP_on_simulation_traj.py \
  --Traj   {traj_num} \
  --ID     1ZMR \
  --PSF    $REFSTRUCT/1zmr_model_clean.pdb \
  --COR    $REFSTRUCT/1zmr_model_clean.pdb \
  --DCD    $DATASTORE/user_input/aa_trajectories/{traj_num}_prod_aa.dcd \
  --resolution aa \
  --contacts calpha \
  --sec_elements $REFSTRUCT/secondary_struc_defs.txt \
  --domain       $REFSTRUCT/domain_def.dat \
  --outdir $DATASTORE/outputs/workflow2/OP_AA \
  --logdir $DATASTORE/outputs/workflow2/OP_AA/logs \
  --start  0 \
  --ent_detection_method 2 \
  --nproc        {nproc} \
  --xp_pdb $REFSTRUCT/1zmr_model_clean.pdb \
  --ops XP
"""
####################################################################################################################################

####################################################################################################################################
run_OP_on_simulation_traj_template_slurm_v3 = """#!/bin/bash
#SBATCH -J {job_name}
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n {nproc}
#SBATCH -t 72:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=20G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

# Step 1: Compute G on the full CG trajectory
# Expected runtime: 12–20 h for G (Topoly, 6667 frames, nproc=10)
# Submit: sbatch assets/slurm/scripts/run_OP_traj{traj_num}.slurm

cd /storage/group/epo2/default/ims86/git_repos/NCLEdetector
source ~/.bashrc
conda activate ncledetector

set -euo pipefail

DATASTORE=/scratch/ims86/NCLEdetector_Datastore
REFSTRUCT=$DATASTORE/user_input/reference_structures

mkdir -p $DATASTORE/outputs/workflow2/OP/G
mkdir -p $DATASTORE/outputs/workflow2/OP/logs

# --- CG: G ---
python scripts/run_OP_on_simulation_traj.py \
  --Traj {traj_num} \
  --PSF  $REFSTRUCT/1zmr_model_clean_ca.psf \
  --COR  $REFSTRUCT/1zmr_model_clean_ca.cor \
  --DCD  $DATASTORE/user_input/cg_trajectories/{traj_num}_prod.dcd \
  --resolution cg \
  --contacts calpha \
  --ID   1ZMR \
  --sec_elements $REFSTRUCT/secondary_struc_defs.txt \
  --domain       $REFSTRUCT/domain_def.dat \
  --outdir       $DATASTORE/outputs/workflow2/OP \
  --logdir       $DATASTORE/outputs/workflow2/OP/logs \
  --start        0 \
  --ent_detection_method 1 \
  --nproc        {nproc} \
  --ops G \
  --no_topoly \
  --chunk_frames 100 --chunk_suffix _chunk 
"""

#############################################################################################################
#############################################################################################################

nproc = 10
for i in range(1, 1001):
    job_name = f"OP_traj{i}"
    slurm_script = run_OP_on_simulation_traj_template_slurm_v2.format(job_name=job_name, traj_num=i, nproc=nproc)
    print(f"Generating SLURM script for trajectory {i}...")
    # print(slurm_script)
    outfile = f"assets/slurm/scripts/run_OP_traj{i}.slurm"
    with open(outfile, "w") as f:
        f.write(slurm_script)
    print(f"SLURM script for trajectory {i} written to {outfile}\n")