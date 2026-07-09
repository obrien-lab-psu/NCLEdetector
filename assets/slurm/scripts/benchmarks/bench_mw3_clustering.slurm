#!/bin/bash
#SBATCH -J bench_mw3_clustering
#SBATCH --partition=basic
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 06:00:00
#SBATCH --account=read_crch_ims86
#SBATCH --mem=64G
#SBATCH --array=0-4
#SBATCH -o assets/slurm/logs/%x-%A_%a.out
#SBATCH -e assets/slurm/logs/%x-%A_%a.err

################################################################################
# Benchmark - Minimal Workflow 3: nonnative entanglement clustering.
# One array task per protein. Requires MW2 OP_last67 to have completed.
# Builds the trajnum->pkl map from the last67 Combined_GE outputs, then runs the
# clustering with the config.
# Config-driven: scripts/configs/benchmarks/<protein>/nonnative_clustering.json
################################################################################
set -eo pipefail

REPO=/storage/group/epo2/default/ims86/git_repos/EntDetect
DATASTORE=/scratch/ims86/EntDetect_Datastore
cd "$REPO"
source ~/.bashrc
conda activate entdetect

PROTEINS=(P61175_6XZ7_S P0A6T9_3A7L_A P45748_1HRU_A P0A6L2_1DHP_A P0AD61_4YNG_C)
ACCS=(P61175 P0A6T9 P45748 P0A6L2 P0AD61)
P=${PROTEINS[$SLURM_ARRAY_TASK_ID]}
ACC=${ACCS[$SLURM_ARRAY_TASK_ID]}

MAP="$DATASTORE/user_input/benchmarks/$P/OP_last67_trajnum2file.txt"
CGE="$DATASTORE/outputs/workflow2/benchmarks/$P/OP_last67/G/Combined_GE"

echo "=== MW3 clustering | protein=$P | building trajnum map | $(date) ==="
echo "trajnum,pklfile" > "$MAP"
for t in $(seq 0 49); do
    pkl="$CGE/${ACC}_traj${t}_GE.pkl"
    if [[ -f "$pkl" ]]; then
        echo "${t},${pkl}" >> "$MAP"
    else
        echo "WARNING: missing pkl for traj $t: $pkl" >&2
    fi
done
echo "map rows: $(( $(wc -l < "$MAP") - 1 ))"

python scripts/run_nonnative_entanglement_clustering.py \
    --config "scripts/configs/benchmarks/$P/nonnative_clustering.json"
echo "=== done $P | $(date) ==="
