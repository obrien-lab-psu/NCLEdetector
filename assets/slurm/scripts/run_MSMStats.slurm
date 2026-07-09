#!/bin/bash
#SBATCH -J workflow2_MSMStats
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 01:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=32G
#SBATCH -o assets/slurm/logs/%x-%j.out
#SBATCH -e assets/slurm/logs/%x-%j.err

################################################################################
# Workflow 2: MSM State Probability Statistics SLURM Script
#
# Computes metastable-state probability evolution for two contrasting cases:
# Case 1: Biologically-informed split (Q >= 0.8 and G <= 0.05)
# Case 2: Random split (negative control)
#
# Uses scripts/run_MSMStats.py with config files for each case.
#
# Expected runtime: 5–10 minutes per case
#
# Memory: 32 GB
# CPUs: 4
#
# Prerequisites:
#   - Base MSM mapping file: 1ZMR_prod_MSMmapping.csv
#   - Labeled MSM files: 1ZMR_prod_MSMmapping_QG_native.csv
#                        1ZMR_prod_MSMmapping_random.csv
#   - Meta set file: 1ZMR_prod_meta_set.csv
#   - All in: /scratch/ims86/EntDetect_Datastore/outputs/workflow2/MSM/
#   - Conda env: entdetect
#   - Repo path: /storage/group/epo2/default/ims86/git_repos/EntDetect
#
# Monitor progress:
#   tail -f assets/slurm/logs/workflow2_MSMStats-<JOBID>.out
#
# Usage:
#   sbatch assets/slurm/scripts/run_MSMStats.slurm
#
################################################################################

# Set up error handling (exit on error)
set -eo pipefail

# Print job info
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

# Navigate to project directory
cd /storage/group/epo2/default/ims86/git_repos/EntDetect

# Activate conda environment
echo "Activating conda environment..."
source ~/.bashrc
conda activate entdetect

# Verify required files exist
echo ""
echo "Verifying prerequisites..."

DATASTORE=/scratch/ims86/EntDetect_Datastore
MSMDIR=$DATASTORE/outputs/workflow2/MSM
STATSDIR=$DATASTORE/outputs/workflow2/MSM_StateProbabilityStats
LOGDIR=$MSMDIR/logs

# Check base MSM files
if [[ ! -f "$MSMDIR/1ZMR_prod_MSMmapping.csv" ]]; then
    echo "ERROR: Base MSM mapping not found: $MSMDIR/1ZMR_prod_MSMmapping.csv"
    exit 1
fi
echo "✓ Found base MSM mapping"

if [[ ! -f "$MSMDIR/1ZMR_prod_meta_set.csv" ]]; then
    echo "ERROR: Meta set file not found: $MSMDIR/1ZMR_prod_meta_set.csv"
    exit 1
fi
echo "✓ Found meta set file"

# Check labeled MSM files
if [[ ! -f "$MSMDIR/1ZMR_prod_MSMmapping_QG_native.csv" ]]; then
    echo "ERROR: Case 1 labeled MSM file not found: $MSMDIR/1ZMR_prod_MSMmapping_QG_native.csv"
    echo "Tip: Run 'python scripts/generate_labeled_msm_cases.py' first"
    exit 1
fi
echo "✓ Found Case 1 labeled MSM file"

if [[ ! -f "$MSMDIR/1ZMR_prod_MSMmapping_random.csv" ]]; then
    echo "ERROR: Case 2 labeled MSM file not found: $MSMDIR/1ZMR_prod_MSMmapping_random.csv"
    echo "Tip: Run 'python scripts/generate_labeled_msm_cases.py' first"
    exit 1
fi
echo "✓ Found Case 2 labeled MSM file"

# Check config files
CFG_CASE1=scripts/configs/workflow2_MSMStats_case1_config.json
CFG_CASE2=scripts/configs/workflow2_MSMStats_case2_config.json

if [[ ! -f "$CFG_CASE1" ]]; then
    echo "ERROR: Config file not found: $CFG_CASE1"
    exit 1
fi
echo "✓ Found Case 1 config file"

if [[ ! -f "$CFG_CASE2" ]]; then
    echo "ERROR: Config file not found: $CFG_CASE2"
    exit 1
fi
echo "✓ Found Case 2 config file"

# Check script
if [[ ! -f "scripts/run_MSMStats.py" ]]; then
    echo "ERROR: Missing script: scripts/run_MSMStats.py"
    exit 1
fi
echo "✓ Found script: scripts/run_MSMStats.py"

echo ""
echo "================================================================================"
echo "Running Workflow 2 MSMStats Analysis"
echo "================================================================================"
echo "Start time: $(date)"
echo "MSM directory: $MSMDIR"
echo "Output base: $STATSDIR"
echo "================================================================================"
echo ""

# Create output and log directories
mkdir -p "$STATSDIR" "$LOGDIR"

# Case 1: Biologically-informed split
echo "================================================================================"
echo "Case 1: Biologically-informed split (Q >= 0.8 and G <= 0.05)"
echo "================================================================================"
echo "Start time: $(date)"
echo ""

python scripts/run_MSMStats.py \
    --config "$CFG_CASE1" \
    2>&1 | tee "$LOGDIR/workflow2_MSMStats_case1.log"

EXIT_CODE_CASE1=${PIPESTATUS[0]}
echo ""
echo "Case 1 exit code: $EXIT_CODE_CASE1"
echo "Completion time: $(date)"
echo ""

# Case 2: Random split (negative control)
echo "================================================================================"
echo "Case 2: Random split (negative control)"
echo "================================================================================"
echo "Start time: $(date)"
echo ""

python scripts/run_MSMStats.py \
    --config "$CFG_CASE2" \
    2>&1 | tee "$LOGDIR/workflow2_MSMStats_case2.log"

EXIT_CODE_CASE2=${PIPESTATUS[0]}
echo ""
echo "Case 2 exit code: $EXIT_CODE_CASE2"
echo "Completion time: $(date)"
echo ""

# Overall exit code (fail if either case failed)
EXIT_CODE=0
if [[ $EXIT_CODE_CASE1 -ne 0 ]] || [[ $EXIT_CODE_CASE2 -ne 0 ]]; then
    EXIT_CODE=1
fi

echo "================================================================================"
echo "Job Completion"
echo "================================================================================"
echo "Exit code: $EXIT_CODE"
echo "Case 1 exit code: $EXIT_CODE_CASE1"
echo "Case 2 exit code: $EXIT_CODE_CASE2"
echo "End time: $(date)"
echo "================================================================================"

if [[ $EXIT_CODE -eq 0 ]]; then
    echo ""
    echo "✓ Workflow 2 MSMStats analysis completed successfully!"
    echo ""
    echo "Output directories:"
    echo "  - Case 1 (QG_native): $STATSDIR"_QG_native
    echo "  - Case 2 (random):    $STATSDIR"_random
    echo ""
    echo "Expected output files per case:"
    echo "  - Plot_StateProbs_TimeSeries.pdf"
    echo "  - 1ZMR_prod_StateProbs_A.dat"
    echo "  - 1ZMR_prod_StateProbs_B.dat"
    echo ""
    echo "Log files:"
    echo "  - $LOGDIR/workflow2_MSMStats_case1.log"
    echo "  - $LOGDIR/workflow2_MSMStats_case2.log"
    echo ""
    echo "View results:"
    echo "  tail -n 50 $LOGDIR/workflow2_MSMStats_case1.log"
    echo "  tail -n 50 $LOGDIR/workflow2_MSMStats_case2.log"
else
    echo ""
    echo "✗ Workflow 2 MSMStats analysis failed!"
    echo ""
    echo "Debug:"
    echo "  tail -n 100 $LOGDIR/workflow2_MSMStats_case1.log"
    echo "  tail -n 100 $LOGDIR/workflow2_MSMStats_case2.log"
fi

exit $EXIT_CODE
