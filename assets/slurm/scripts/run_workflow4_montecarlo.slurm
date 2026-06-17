#!/bin/bash
#SBATCH -J workflow4_montecarlo
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 24:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=100G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

################################################################################
# Workflow 4: Step 6 Monte Carlo SLURM Script
#
# Runs scripts/run_montecarlo.py using a single combined residue-level design matrix.
#
# CPUs: 20
# Memory: 100 GB
#
# Prerequisites:
#   - Design matrix: /scratch/ims86/EntDetect_Datastore/outputs/workflow4/residue_dataframes_workflow4.csv
#   - Gene list:
#     /scratch/ims86/EntDetect_Datastore/user_input/experimental_data/Gene_lists/AF/
#     AF_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt
#
# Monitor progress:
#   tail -f assets/slurm/logs/workflow4_montecarlo-<JOBID>.out
#
# Usage:
#   sbatch assets/slurm/scripts/run_workflow4_montecarlo.slurm
#
################################################################################

set -eo pipefail

echo "================================================================================"
echo "SLURM Job Information"
echo "================================================================================"
echo "Job ID: ${SLURM_JOB_ID:-unknown}"
echo "Job Name: ${SLURM_JOB_NAME:-unknown}"
echo "Partition: ${SLURM_PARTITION:-unknown}"
echo "Nodes: ${SLURM_NNODES:-1}"
echo "CPUs: ${SLURM_NPROCS:-1}"
echo "Memory: ${SLURM_MEM_PER_NODE:-unknown} MB"
echo "Time limit: ${SLURM_TIMELIMIT:-unknown} minutes"
echo "Start time: $(date)"
echo "================================================================================"
echo ""

cd /storage/group/epo2/default/ims86/git_repos/EntDetect

echo "Activating conda environment..."
source ~/.bashrc
conda activate entdetect

echo ""
echo "Verifying prerequisites..."

DATASTORE=/scratch/ims86/EntDetect_Datastore
OUTROOT=$DATASTORE/outputs/workflow4
DATAFRAME_FILE=$OUTROOT/residue_dataframes_workflow4.csv
GENE_LIST=$DATASTORE/user_input/experimental_data/Gene_lists/AF/AF_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt
OUTDIR=$OUTROOT/monte_carlo
LOGDIR=$OUTROOT/logs
TAG=Ecoli_population_mc

# For faster test runs, lower STEPS (for example 10000) before submitting.
STEPS=100000
N_GROUPS=4
C1=1.0
C2=2.5
BETA=0.05

if [[ ! -f "$DATAFRAME_FILE" ]]; then
    echo "ERROR: Residue design matrix file not found: $DATAFRAME_FILE"
    echo "Build the combined matrix first (Workflow 4 Step 4)."
    exit 1
fi
echo "Found residue design matrix: $DATAFRAME_FILE"

if [[ ! -f "$GENE_LIST" ]]; then
    echo "ERROR: Gene list not found: $GENE_LIST"
    exit 1
fi
n_gene=$(grep -cv '^\s*$' "$GENE_LIST")
echo "Found gene list: $GENE_LIST ($n_gene non-empty lines)"

mkdir -p "$OUTDIR" "$LOGDIR"

echo ""
echo "================================================================================"
echo "Running Workflow 4 Step 6: Monte Carlo"
echo "================================================================================"
echo "Start time: $(date)"
echo "Design matrix: $DATAFRAME_FILE"
echo "Gene list: $GENE_LIST"
echo "Output dir: $OUTDIR"
echo "Log dir: $LOGDIR"
echo "Tag: $TAG"
echo "Steps: $STEPS"
echo "N groups: $N_GROUPS"
echo "================================================================================"
echo ""

python scripts/run_montecarlo.py \
    --dataframe_files "$DATAFRAME_FILE" \
    --outpath "$OUTDIR" \
    --gene_list "$GENE_LIST" \
    --tag "$TAG" \
    --reg_formula "cut_C_Rall ~ region + AA" \
    --response_var cut_C_Rall \
    --test_var region \
    --n_groups "$N_GROUPS" \
    --steps "$STEPS" \
    --C1 "$C1" \
    --C2 "$C2" \
    --beta "$BETA" \
    --log_level INFO \
    --logdir "$LOGDIR" \
    2>&1 | tee "$LOGDIR/workflow4_montecarlo.log"

EXIT_CODE=${PIPESTATUS[0]}

echo ""
echo "================================================================================"
echo "Job Completion"
echo "================================================================================"
echo "Exit code: $EXIT_CODE"
echo "End time: $(date)"
echo "================================================================================"

if [[ $EXIT_CODE -eq 0 ]]; then
    echo ""
    echo "Workflow 4 Step 6 completed successfully"
    echo "Outputs saved to: $OUTDIR"
    echo "Main log: $LOGDIR/workflow4_montecarlo.log"
else
    echo ""
    echo "Workflow 4 Step 6 failed with exit code $EXIT_CODE"
    echo "Check logs:"
    echo "  tail -n 100 $LOGDIR/workflow4_montecarlo.log"
fi

exit $EXIT_CODE
