#!/bin/bash
#SBATCH -J workflow4_nativeNCLE_batch
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 36:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=100G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

################################################################################
# Workflow 4: Native NCLE Batch (Steps 3+4) SLURM Script
#
# Runs scripts/run_workflow4_nativeNCLE_batch.py to process all PDBs in a
# directory that are also present in a supplied gene list, launching
# run_nativeNCLE.py jobs in parallel.
# Also builds Workflow 4 Step 4 regression design matrices by combining
# NCLE-derived region labels with per-residue feature tables.
#
# CPUs: 20
# Memory: 100 GB
#
# Prerequisites:
#   - PDB dir: /scratch/ims86/EntDetect_Datastore/user_input/proteome_structures/AF/
#   - Gene list:
#     /scratch/ims86/EntDetect_Datastore/user_input/experimental_data/Gene_lists/AF/
#     AF_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt
#
# Monitor progress:
#   tail -f assets/slurm/logs/workflow4_nativeNCLE_batch-<JOBID>.out
#
# Usage:
#   sbatch assets/slurm/scripts/run_workflow4_nativeNCLE_batch.slurm
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
PDB_DIR=$DATASTORE/user_input/proteome_structures/AF
GENE_LIST=$DATASTORE/user_input/experimental_data/Gene_lists/AF/AF_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt
RES_FEATURES_AF=$DATASTORE/user_input/experimental_data/PDB_residue_features/AF/residueFeatures.csv
OUTDIR=$DATASTORE/outputs/workflow4/nativeNCLE_all
LOGDIR=$DATASTORE/outputs/workflow4/logs

if [[ ! -d "$PDB_DIR" ]] || [[ -z $(ls -1 "$PDB_DIR"/*.pdb 2>/dev/null | head -1) ]]; then
    echo "ERROR: No .pdb files found in $PDB_DIR"
    exit 1
fi
n_pdb=$(ls -1 "$PDB_DIR"/*.pdb 2>/dev/null | wc -l)
echo "✓ Found $n_pdb PDB files"

if [[ ! -f "$GENE_LIST" ]]; then
    echo "ERROR: Gene list not found: $GENE_LIST"
    exit 1
fi
n_gene=$(grep -cv '^\s*$' "$GENE_LIST")
echo "✓ Found gene list: $GENE_LIST ($n_gene non-empty lines)"

if [[ ! -f "$RES_FEATURES_AF" ]]; then
    echo "ERROR: AF residue features file not found: $RES_FEATURES_AF"
    exit 1
fi
echo "✓ Found residue features: $RES_FEATURES_AF"

mkdir -p "$OUTDIR" "$LOGDIR"

echo ""
echo "================================================================================"
echo "Running Workflow 4 Native NCLE Batch"
echo "================================================================================"
echo "Start time: $(date)"
echo "PDB dir: $PDB_DIR"
echo "Gene list: $GENE_LIST"
echo "Residue features: $RES_FEATURES_AF"
echo "Output dir: $OUTDIR"
echo "Log dir: $LOGDIR"
echo "Parallel jobs (--nproc): 20"
echo "================================================================================"
echo ""

python scripts/run_workflow4_nativeNCLE_batch.py \
    --pdb_dir "$PDB_DIR" \
    --gene_list "$GENE_LIST" \
    --outdir "$OUTDIR" \
    --organism Ecoli \
    --model AF \
    --contacts heavy \
    --resolution aa \
    --ent_detection_method 3 \
    --residue_features_file "$RES_FEATURES_AF" \
    --reg_formula "cut_C_Rall ~ AA + region" \
    --design_matrix_file "$DATASTORE/outputs/workflow4/residue_dataframes_workflow4.csv" \
    --nproc 20 \
    --log_level INFO \
    --logdir "$LOGDIR" \
    2>&1 | tee "$LOGDIR/workflow4_nativeNCLE_batch.log"

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
    echo "✓ Workflow 4 nativeNCLE batch completed successfully!"
    echo ""
    echo "Outputs saved to: $OUTDIR"
    echo "Main log: $LOGDIR/workflow4_nativeNCLE_batch.log"
    echo ""
    echo "View results:"
    echo "  tail -n 50 $LOGDIR/workflow4_nativeNCLE_batch.log"
else
    echo ""
    echo "✗ Workflow 4 nativeNCLE batch failed with exit code $EXIT_CODE"
    echo "Check logs:"
    echo "  tail -n 100 $LOGDIR/workflow4_nativeNCLE_batch.log"
fi

exit $EXIT_CODE
