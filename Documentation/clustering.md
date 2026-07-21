# Clustering Module

The `clustering` module provides comprehensive functionality for clustering both native and non-native entanglements, as well as building Markov State Models (MSMs) from simulation ensembles.

## Overview

This module enables:
- **Native entanglement clustering**: Remove degeneracies from identified entanglements
- **Non-native entanglement clustering**: Analyze dynamic entanglement changes across trajectories  
- **MSM construction**: Build Markov State Models for folding pathway analysis
- **Cross-organism analysis**: Support for different protein datasets (E. coli, Human, Yeast)

## Classes

### ClusterNativeEntanglements

Clusters native entanglements to remove degeneracies and identify unique topological features.

#### Initialization
```python
from NCLEdetector.clustering import ClusterNativeEntanglements

clustering = ClusterNativeEntanglements(organism='Ecoli', cut_off=57)
```

**Parameters:**
- `organism` (str): Target organism ('Ecoli', 'Human', 'Yeast') (default: 'Ecoli')
- `cut_off` (int): Distance cutoff for clustering (default: 57, organism-specific)

#### Key Methods

##### Cluster_NativeEntanglements()
Clusters high-quality native entanglements based on structural similarity.

```python
result = clustering.Cluster_NativeEntanglements(
    GE_filepath="hq_entanglements.csv",
    outdir="results/clustered/",
    outfile="clustered_entanglements.csv"
)
```

**Parameters:**
- `GE_filepath` (str): Path to high-quality entanglement data
- `outdir` (str): Output directory (default: './')
- `outfile` (str): Output filename (default: 'Cluster_NativeEntanglements.txt')

**Returns:**
- Dictionary with clustered entanglement data and file paths

### ClusterNonNativeEntanglements

Analyzes dynamic entanglement changes across multiple simulation trajectories.

#### Initialization  
```python
from NCLEdetector.clustering import ClusterNonNativeEntanglements

nn_clustering = ClusterNonNativeEntanglements(
    pkl_file_path="pkl_files/",
    trajnum2pklfile_path="mapping.txt", 
    traj_dir_prefix="trajectories/",
    outdir="results/nonnative/"
)
```

**Parameters:**
- `pkl_file_path` (str): Directory containing pickled entanglement data
- `trajnum2pklfile_path` (str): File mapping trajectory numbers to pkl files
- `traj_dir_prefix` (str): Prefix path for trajectory directories
- `outdir` (str): Output directory

#### Key Methods

##### cluster()
Performs clustering analysis of non-native entanglement changes.

```python
nn_clustering.cluster(start_frame=6332)
```

**Parameters:**
- `start_frame` (int): Starting frame for analysis (default: 0)

### MSMNonNativeEntanglementClustering

Builds Markov State Models from ensemble simulation data with entanglement analysis.

#### Initialization
```python
from NCLEdetector.clustering import MSMNonNativeEntanglementClustering

MSM = MSMNonNativeEntanglementClustering(
    outdir="results/MSM/",
    ID="protein_msm", 
    OPpath="OP_data/",
    start=0,
    n_large_states=10,
    rm_traj_list=[65, 75, 155],
    lagtime=20
)
```

**Parameters:**
- `outdir` (str): Output directory for MSM results
- `ID` (str): Identifier for MSM analysis
- `OPpath` (str): Path to order parameter data
- `start` (int): Starting frame for analysis
- `n_large_states` (int): Number of large states for MSM
- `rm_traj_list` (list): Trajectory indices to exclude
- `lagtime` (int): Lag time for MSM construction

#### Key Methods

##### run()
Executes complete MSM construction and analysis pipeline.

```python
MSM.run()
```

## Usage Examples

### Native Entanglement Clustering

```python
from NCLEdetector.clustering import ClusterNativeEntanglements

# Initialize for E. coli proteins
clustering = ClusterNativeEntanglements(organism='Ecoli')

# Cluster high-quality native entanglements
result = clustering.Cluster_NativeEntanglements(
    GE_filepath="results/hq_entanglements.csv",
    outdir="results/clustered/",
    outfile="protein_clustered.csv"
)

print(f"Clustered data saved to: {result['outfile']}")
```

### Non-Native Entanglement Analysis

```python
from NCLEdetector.clustering import ClusterNonNativeEntanglements

# Set up non-native clustering
nn_clustering = ClusterNonNativeEntanglements(
    pkl_file_path="simulation_data/pkl/",
    trajnum2pklfile_path="trajectory_mapping.txt",
    traj_dir_prefix="trajectories/",
    outdir="results/nonnative_analysis/"
)

# Perform clustering starting from equilibration
nn_clustering.cluster(start_frame=6332)
```

### MSM Construction

```python
from NCLEdetector.clustering import MSMNonNativeEntanglementClustering

# Build MSM with entanglement analysis
MSM = MSMNonNativeEntanglementClustering(
    outdir="results/MSM_analysis/",
    ID="1ZMR_folding",
    OPpath="results/OP_data/", 
    start=0,
    n_large_states=10,
    rm_traj_list=[65, 75, 155, 162],  # Remove problematic trajectories
    lagtime=20
)

# Run complete MSM pipeline
MSM.run()
```

## Integration with Scripts

### run_nativeNCLE.py
Complete native entanglement identification and clustering:

```bash
python scripts/run_nativeNCLE.py \
  --struct structure.pdb \
  --outdir results/native_analysis/ \
  --ID protein_name \
  --organism Ecoli
```

### run_nonnative_entanglement_clustering.py
Analyze dynamic entanglement changes across trajectories:

```bash
python scripts/run_nonnative_entanglement_clustering.py \
  --outdir results/nonnative_clustering/ \
  --pkl_file_path simulation_data/pkl/ \
  --trajnum2pklfile_path trajectory_mapping.txt \
  --traj_dir_prefix trajectories/
```

### run_MSM.py  
Build comprehensive MSMs with entanglement tracking:

```bash
python scripts/run_MSM.py \
  --outdir results/MSM/ \
  --ID protein_msm \
  --OPpath OP_analysis/ \
  --start 0 \
  --n_large_states 10 \
  --lagtime 20 \
  --rm_traj_list 65 75 155
```

## Output Files

### Native Clustering Output
- `*_clustered_entanglements.csv`: Clustered entanglement data
- `cluster_representatives.csv`: Representative structures for each cluster
- `clustering_metrics.txt`: Clustering statistics and quality measures

### Non-Native Clustering Output  
- `nonnative_clusters.pkl`: Pickled clustering results
- `entanglement_transitions.csv`: Dynamic entanglement change data
- `trajectory_analysis.txt`: Per-trajectory entanglement statistics

### MSM Output
- `msm_states.csv`: MSM state definitions and populations
- `transition_matrix.npy`: State transition probability matrix
- `representative_structures/`: PDB files for state representatives
- `pathway_analysis.txt`: Folding pathway statistics

## Algorithm Details

### Native Clustering Algorithm
1. **Feature Extraction**: Calculate structural descriptors for each entanglement
2. **Distance Matrix**: Compute pairwise similarities between entanglements  
3. **Hierarchical Clustering**: Group similar entanglements using linkage criteria
4. **Cluster Validation**: Assess cluster quality and select optimal number

### Non-Native Clustering
1. **Trajectory Alignment**: Synchronize entanglement data across trajectories
2. **Change Detection**: Identify formation/breaking of entanglements
3. **Pattern Recognition**: Cluster similar entanglement change patterns
4. **Temporal Analysis**: Analyze timing and frequency of changes

### MSM Construction
1. **State Definition**: Define microstates using order parameters
2. **Trajectory Discretization**: Assign trajectory frames to states  
3. **Transition Counting**: Count transitions between states
4. **Rate Matrix Estimation**: Calculate transition probabilities
5. **Coarse-Graining**: Identify metastable macrostates

## Performance Notes

- **Scalability**: Native clustering scales O(N²) with entanglement number
- **Memory Usage**: Non-native clustering requires substantial RAM for large ensembles
- **Parallelization**: MSM construction benefits from multiprocessing
- **Storage**: Output files can be large for extensive simulation data

## Dependencies

- Scikit-learn: Clustering algorithms
- NumPy/SciPy: Numerical computations  
- Pandas: Data manipulation
- PyEMMA: Markov model construction
- NetworkX: Graph analysis for pathways



