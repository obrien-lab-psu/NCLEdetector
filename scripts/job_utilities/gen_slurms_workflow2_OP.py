import os,sys
import glob

####################################################################################################################################
## TEMPLATES #######################################################################################################################
####################################################################################################################################


####################################################################################################################################
run_OP_on_simulation_traj_template_slurm = """#!/bin/bash
#SBATCH -J {job_name}
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n {nproc}
#SBATCH -t 24:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=20G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

# Submit: sbatch assets/slurm/scripts/run_OP_traj{traj_num}.slurm

cd /storage/group/epo2/default/ims86/git_repos/EntDetect
source ~/.bashrc
conda activate entdetect

set -euo pipefail

DATASTORE=/scratch/ims86/EntDetect_Datastore
REFSTRUCT=$DATASTORE/user_input/reference_structures

mkdir -p $DATASTORE/outputs/workflow2/OP_last67/G
mkdir -p $DATASTORE/outputs/workflow2/OP_last67/logs

# --- CG: G Q K ---
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
  --outdir       $DATASTORE/outputs/workflow2/OP_last67 \
  --logdir       $DATASTORE/outputs/workflow2/OP_last67/logs \
  --start        6600 \
  --ent_detection_method 2 \
  --nproc        {nproc} \
  --ops G Q K 
"""

#############################################################################################################
#############################################################################################################

nproc = 10
for i in range(1, 1001):
    job_name = f"OP_traj{i}"
    slurm_script = run_OP_on_simulation_traj_template_slurm.format(job_name=job_name, traj_num=i, nproc=nproc)
    print(f"Generating SLURM script for trajectory {i}...")
    # print(slurm_script)
    outfile = f"assets/slurm/scripts/run_OP_traj{i}.slurm"
    with open(outfile, "w") as f:
        f.write(slurm_script)
    print(f"SLURM script for trajectory {i} written to {outfile}\n")