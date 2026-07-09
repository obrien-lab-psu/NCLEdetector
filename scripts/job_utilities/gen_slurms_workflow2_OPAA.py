import os,sys
import glob

####################################################################################################################################
## TEMPLATES #######################################################################################################################
####################################################################################################################################

outdir = 'OP_AA'
start = 0
####################################################################################################################################
run_OP_on_simulation_traj_template_slurm_XPonly = """#!/bin/bash
#SBATCH -J {job_name}
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n {nproc}
#SBATCH -t 10:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=20G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

# Step 2: Compute XP on the all-atom back-mapped frames
#
# Expected runtime: 12–20 h for G (Topoly, 6667 frames, nproc=10)
# Submit: sbatch assets/slurm/scripts/run_OP_traj{traj_num}.slurm

cd /storage/group/epo2/default/ims86/git_repos/EntDetect
source ~/.bashrc
conda activate entdetect

set -euo pipefail

DATASTORE=/scratch/ims86/EntDetect_Datastore
REFSTRUCT=$DATASTORE/user_input/reference_structures

mkdir -p $DATASTORE/outputs/workflow2/{outdir}/SASA
mkdir -p $DATASTORE/outputs/workflow2/{outdir}/XP
mkdir -p $DATASTORE/outputs/workflow2/{outdir}/logs

# --- AA: XP ---
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
  --outdir $DATASTORE/outputs/workflow2/{outdir} \
  --logdir $DATASTORE/outputs/workflow2/{outdir}/logs \
  --start  {start} \
  --ent_detection_method 2 \
  --nproc        {nproc} \
  --xp_pdb $REFSTRUCT/1zmr_model_clean.pdb \
  --ops XP
"""


####################################################################################################################################
run_OP_on_simulation_traj_template_slurm_SASAonly = """#!/bin/bash
#SBATCH -J {job_name}
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n {nproc}
#SBATCH -t 10:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=20G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

# Step 2: Compute SASA on the all-atom back-mapped frames
#
# Expected runtime: 12–20 h for G (Topoly, 6667 frames, nproc=10)
# Submit: sbatch assets/slurm/scripts/run_OP_traj{traj_num}.slurm

cd /storage/group/epo2/default/ims86/git_repos/EntDetect
source ~/.bashrc
conda activate entdetect

set -euo pipefail

DATASTORE=/scratch/ims86/EntDetect_Datastore
REFSTRUCT=$DATASTORE/user_input/reference_structures

mkdir -p $DATASTORE/outputs/workflow2/{outdir}/SASA
mkdir -p $DATASTORE/outputs/workflow2/{outdir}/XP
mkdir -p $DATASTORE/outputs/workflow2/{outdir}/logs

# --- AA: SASA ---
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
  --outdir $DATASTORE/outputs/workflow2/{outdir} \
  --logdir $DATASTORE/outputs/workflow2/{outdir}/logs \
  --start  {start} \
  --ent_detection_method 2 \
  --nproc        {nproc} \
  --xp_pdb $REFSTRUCT/1zmr_model_clean.pdb \
  --ops SASA
"""
####################################################################################################################################

####################################################################################################################################
run_OP_on_simulation_traj_template_slurm_XPSASA = """#!/bin/bash
#SBATCH -J {job_name}
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n {nproc}
#SBATCH -t 10:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=20G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

# Step 2: Compute SASA and XP on the all-atom back-mapped frames
#
# Expected runtime: 12–20 h for G (Topoly, 6667 frames, nproc=10)
# Submit: sbatch assets/slurm/scripts/run_OP_traj{traj_num}.slurm

cd /storage/group/epo2/default/ims86/git_repos/EntDetect
source ~/.bashrc
conda activate entdetect

set -euo pipefail

DATASTORE=/scratch/ims86/EntDetect_Datastore
REFSTRUCT=$DATASTORE/user_input/reference_structures

mkdir -p $DATASTORE/outputs/workflow2/{outdir}/SASA
mkdir -p $DATASTORE/outputs/workflow2/{outdir}/XP
mkdir -p $DATASTORE/outputs/workflow2/{outdir}/logs

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
  --outdir $DATASTORE/outputs/workflow2/{outdir} \
  --logdir $DATASTORE/outputs/workflow2/{outdir}/logs \
  --start  {start} \
  --ent_detection_method 2 \
  --nproc        {nproc} \
  --xp_pdb $REFSTRUCT/1zmr_model_clean.pdb \
  --ops XP SASA
"""

#############################################################################################################
#############################################################################################################

nproc = 1
for i in range(1, 1001):
    job_name = f"OP_traj{i}"

    ## check if .XP for this traj already exists /scratch/ims86/EntDetect_Datastore/outputs/workflow2/{outdir}/XP/1ZMR_Traj1.XP and if so use the SASA-only template
    ## else use the SASA+XP template
    xp_file = f"/scratch/ims86/EntDetect_Datastore/outputs/workflow2/{outdir}/XP/1ZMR_Traj{i}.XP"
    SASA_file = f"/scratch/ims86/EntDetect_Datastore/outputs/workflow2/{outdir}/SASA/1ZMR_Traj{i}.SASA"

    if os.path.exists(SASA_file) and os.path.exists(xp_file):
        print(f"Both SASA and XP files for trajectory {i} already exist at {SASA_file} and {xp_file}. Skipping SLURM script generation.")
        continue

    if os.path.exists(xp_file) and not os.path.exists(SASA_file):
        print(f"XP file for trajectory {i} already exists at {xp_file}. Using SASA-only template.")
        slurm_script = run_OP_on_simulation_traj_template_slurm_SASAonly.format(job_name=job_name, traj_num=i, nproc=nproc, outdir=outdir, start=start)

    elif not os.path.exists(xp_file) and os.path.exists(SASA_file):
        print(f"SASA file for trajectory {i} already exists at {SASA_file}. Using XP-only template.")
        slurm_script = run_OP_on_simulation_traj_template_slurm_XPonly.format(job_name=job_name, traj_num=i, nproc=nproc, outdir=outdir, start=start)

    else:
        print(f"XP & SASA file for trajectory {i} does NOT exist. Using SASA+XP template.")
        slurm_script = run_OP_on_simulation_traj_template_slurm_XPSASA.format(job_name=job_name, traj_num=i, nproc=nproc, outdir=outdir, start=start)

    print(f"Generating SLURM script for trajectory {i}...")
    # print(slurm_script)
    outfile = f"assets/slurm/scripts/run_{outdir}_traj{i}.slurm"
    with open(outfile, "w") as f:
        f.write(slurm_script)
    print(f"SLURM script for trajectory {i} written to {outfile}\n")