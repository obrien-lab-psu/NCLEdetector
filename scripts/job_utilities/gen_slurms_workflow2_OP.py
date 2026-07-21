import os

####################################################################################################################################
## TEMPLATES #######################################################################################################################
####################################################################################################################################


####################################################################################################################################
run_OP_on_simulation_traj_cg_template_slurm = """#!/bin/bash
#SBATCH -J {job_name}
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n {nproc}
#SBATCH -t 24:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=10G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

# Submit: sbatch assets/slurm/scripts/run_OP_traj{traj_num}.slurm

cd /storage/group/epo2/default/ims86/git_repos/EntDetect
source ~/.bashrc
conda activate entdetect

set -euo pipefail

DATASTORE=/scratch/ims86/EntDetect_Datastore
REFSTRUCT=$DATASTORE/user_input/reference_structures
CG_TRAJ_DIR=$DATASTORE/user_input/cg_trajectories

CFG_OP_LAST335=scripts/configs/workflow2_OP_last335_config.json

mkdir -p $DATASTORE/outputs/workflow2/OP_last335/G
mkdir -p $DATASTORE/outputs/workflow2/OP_last335/logs

# 2) CG last-335-frame run (Q, G, K)
python scripts/run_OP_on_simulation_traj.py \
  --config $CFG_OP_LAST335 \
  --Traj {traj_num} \
  --DCD $CG_TRAJ_DIR/{traj_num}_prod.dcd
"""

run_OP_on_simulation_traj_aa_template_slurm = """#!/bin/bash
#SBATCH -J {job_name}
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n {nproc}
#SBATCH -t 24:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=10G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

# Submit: sbatch assets/slurm/scripts/run_OP_AA_traj{traj_num}.slurm

cd /storage/group/epo2/default/ims86/git_repos/EntDetect
source ~/.bashrc
conda activate entdetect

set -euo pipefail

DATASTORE=/scratch/ims86/EntDetect_Datastore
AA_TRAJ_DIR=$DATASTORE/user_input/aa_trajectories

CFG_OP_AA_LAST335=scripts/configs/workflow2_OP_AA_last335_config.json

mkdir -p $DATASTORE/outputs/workflow2/OP_AA_last335/SASA
mkdir -p $DATASTORE/outputs/workflow2/OP_AA_last335/XP
mkdir -p $DATASTORE/outputs/workflow2/OP_AA_last335/logs

# 3) AA run (SASA, XP)
python scripts/run_OP_on_simulation_traj.py \
  --config $CFG_OP_AA_LAST335 \
  --Traj {traj_num} \
  --DCD $AA_TRAJ_DIR/{traj_num}_prod_aa.dcd
"""

#############################################################################################################
#############################################################################################################

nproc = 10
for i in range(1, 1001):
  print(f"Generating SLURM scripts for trajectory {i}...")

  cg_job_name = f"OP_traj{i}"
  cg_script = run_OP_on_simulation_traj_cg_template_slurm.format(job_name=cg_job_name, traj_num=i, nproc=nproc)
  cg_outfile = f"assets/slurm/scripts/run_OP_traj{i}.slurm"
  with open(cg_outfile, "w") as f:
    f.write(cg_script)

  aa_job_name = f"OPAA_traj{i}"
  aa_script = run_OP_on_simulation_traj_aa_template_slurm.format(job_name=aa_job_name, traj_num=i, nproc=nproc)
  aa_outfile = f"assets/slurm/scripts/run_OP_AA_traj{i}.slurm"
  with open(aa_outfile, "w") as f:
    f.write(aa_script)

  print(f"  CG script written: {cg_outfile}")
  print(f"  AA script written: {aa_outfile}\n")