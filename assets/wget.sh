#!/usr/bin/env bash
set -eu

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <DATASTORE_PATH>"
    echo "  Example: $0 /scratch/myuser/NCLEdetector_Datastore"
    exit 1
fi

DATASTORE="$1"

# Create target directory
mkdir -p "$DATASTORE"
cd "$DATASTORE"
echo "Downloading NCLEdetector Datastore to $DATASTORE"

# Download all archives (version 1.0.0)
BASE_URL="https://downloads.ncleweb.org/NCLEdetector_Datastore/v1.0.0"

# User input data (36 GB) - Required for all tutorials 
wget "$BASE_URL/user_input.tar.gz" && tar -xzf user_input.tar.gz

# Workflow 1 analyses (16 KB total)
mkdir -p "${DATASTORE}/outputs/workflow1/"
wget "$BASE_URL/outputs/workflow1/NCLE.tar.gz" && tar -xzf NCLE.tar.gz -C "${DATASTORE}/outputs/workflow1/"
wget "$BASE_URL/outputs/workflow1/HQ_NCLE.tar.gz" && tar -xzf HQ_NCLE.tar.gz -C "${DATASTORE}/outputs/workflow1/"
wget "$BASE_URL/outputs/workflow1/clustered_HQ_NCLE.tar.gz" && tar -xzf clustered_HQ_NCLE.tar.gz -C "${DATASTORE}/outputs/workflow1/"
wget "$BASE_URL/outputs/workflow1/clustered_HQ_NCLE_features.tar.gz" && tar -xzf clustered_HQ_NCLE_features.tar.gz -C "${DATASTORE}/outputs/workflow1/"

# Workflow 2 analyses (48 GB total)
mkdir -p "${DATASTORE}/outputs/workflow2/"
# wget"$BASE_URL/outputs/workflow2/OP.tar.gz" && tar -xzf outputs/workflow2/OP.tar.gz ## (OPTIONAL) if you want to be able to recreate the full analysis in the paper
wget "$BASE_URL/outputs/workflow2/OP_AA_last335.tar.gz" && tar -xzf OP_AA_last335.tar.gz -C "${DATASTORE}/outputs/workflow2/"
wget "$BASE_URL/outputs/workflow2/OP_last335.tar.gz" && tar -xzf OP_last335.tar.gz -C "${DATASTORE}/outputs/workflow2/"
wget "$BASE_URL/outputs/workflow2/OP_last10.tar.gz" && tar -xzf OP_last10.tar.gz -C "${DATASTORE}/outputs/workflow2/"
wget "$BASE_URL/outputs/workflow2/MSM.tar.gz" && tar -xzf MSM.tar.gz -C "${DATASTORE}/outputs/workflow2/"
wget "$BASE_URL/outputs/workflow2/MSM_StateProbabilityStats_QG_native.tar.gz" && tar -xzf MSM_StateProbabilityStats_QG_native.tar.gz -C "${DATASTORE}/outputs/workflow2/"
wget "$BASE_URL/outputs/workflow2/MSM_StateProbabilityStats_random.tar.gz" && tar -xzf MSM_StateProbabilityStats_random.tar.gz -C "${DATASTORE}/outputs/workflow2/"
wget "$BASE_URL/outputs/workflow2/nonnative_clustering_last335.tar.gz" && tar -xzf nonnative_clustering_last335.tar.gz -C "${DATASTORE}/outputs/workflow2/"
wget "$BASE_URL/outputs/workflow2/FoldingPathway_QG_native.tar.gz" && tar -xzf FoldingPathway_QG_native.tar.gz -C "${DATASTORE}/outputs/workflow2/"
wget "$BASE_URL/outputs/workflow2/FoldingPathway_random.tar.gz" && tar -xzf FoldingPathway_random.tar.gz -C "${DATASTORE}/outputs/workflow2/"

# Workflow 3 analyses (1.2 GB)
mkdir -p "${DATASTORE}/outputs/workflow3/"
wget "$BASE_URL/outputs/workflow3/MassSpec_ConsistencyTest.tar.gz" && tar -xzf MassSpec_ConsistencyTest.tar.gz -C "${DATASTORE}/outputs/workflow3/"

# Workflow 4 analyses (8.5 MB)
mkdir -p "${DATASTORE}/outputs/workflow4/"
wget "$BASE_URL/outputs/workflow4/nativeNCLE_all.tar.gz" && tar -xzf nativeNCLE_all.tar.gz -C "${DATASTORE}/outputs/workflow4/"
wget "$BASE_URL/outputs/workflow4/monte_carlo.tar.gz" && tar -xzf monte_carlo.tar.gz -C "${DATASTORE}/outputs/workflow4/"
wget "$BASE_URL/outputs/workflow4/population_modeling.tar.gz" && tar -xzf population_modeling.tar.gz -C "${DATASTORE}/outputs/workflow4/"