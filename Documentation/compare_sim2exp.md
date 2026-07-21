# Simulation-Experiment Comparison Module

The `compare_sim2exp` module provides comprehensive functionality for validating computational protein folding models against experimental mass spectrometry data, enabling quantitative assessment of simulation accuracy and identification of consistent structural states.

## Overview

This module enables:
- **LiP-MS integration**: Compare Limited Proteolysis-Mass Spectrometry data with simulation flexibility
- **XL-MS validation**: Cross-linking Mass Spectrometry distance constraint validation  
- **Consistency testing**: Statistical assessment of simulation-experiment agreement
- **Representative structure selection**: Identify structures that best match experimental signals
- **Multi-state analysis**: Analyze consistency across MSM-derived metastable states

## Classes

### MassSpec

Main class for comparing simulation results with experimental mass spectrometry data through comprehensive consistency testing and representative structure selection.

#### Initialization
```python
from NCLEdetector.compare_sim2exp import MassSpec

MS = MassSpec(
    msm_data_file="msm_mapping.csv",
    meta_dist_file="meta_distances.npz",
    LiPMS_exp_file="LiPMS_experimental.xlsx",
    sasa_data_file="SASA_data.npy",
    XLMS_exp_file="XLMS_experimental.xlsx", 
    dist_data_file="distance_data.npy",
    cluster_data_file="cluster_data.npz",
    OPpath="OP_analysis/",
    AAdcd_dir="trajectories/",
    native_AA_pdb="native_structure.pdb",
    native_state_idx=0,
    state_idx_list=[0, 1, 2, 3, 4],
    rm_traj_list=[65, 75, 155],
    outdir="results/consistency/",
    ID="protein_analysis",
    start=0,
    end=1000,
    stride=1,
    verbose=True,
    num_perm=1000,
    n_boot=100,
    lag_frame=20,
    nproc=10
)
```

**Parameters:**
- `msm_data_file` (str): Path to MSM state mapping data (CSV format)
- `meta_dist_file` (str): Path to meta-distance data (NumPy npz format)
- `LiPMS_exp_file` (str): Path to LiP-MS experimental data (Excel format)
- `sasa_data_file` (str): Path to SASA data from simulations (NumPy npy format)
- `XLMS_exp_file` (str): Path to XL-MS experimental data (Excel format)
- `dist_data_file` (str): Path to distance data from simulations (NumPy npy format)
- `cluster_data_file` (str): Path to entanglement cluster data (NumPy npz format)
- `OPpath` (str): Directory containing order parameter files (G and Q data)
- `AAdcd_dir` (str): Directory containing all-atom trajectory files
- `native_AA_pdb` (str): Path to native all-atom reference structure
- `native_state_idx` (int): Index identifying the native state in MSM
- `state_idx_list` (list): List of MSM state indices to analyze
- `rm_traj_list` (list): Trajectory indices to exclude from analysis
- `outdir` (str): Output directory for results (default: './')
- `ID` (str): Analysis identifier for output files (default: '')
- `start` (int): Starting frame for analysis (default: 0)
- `end` (int): Ending frame for analysis (default: large number)
- `stride` (int): Frame stride for analysis (default: 1)
- `verbose` (bool): Enable verbose output (default: False)
- `num_perm` (int): Number of permutations for statistical tests
- `n_boot` (int): Number of bootstrap samples
- `lag_frame` (int): Lag time in frames
- `nproc` (int): Number of parallel processes

#### Key Methods

##### load_OP()
Loads order parameter data for consistency analysis.

```python
MS.load_OP(start=0, end=99999999999)
```

**Parameters:**
- `start` (int): Starting frame (default: 0)
- `end` (int): Ending frame (default: large number)

##### LiP_XL_MS_ConsistencyTest()
Performs comprehensive consistency testing between simulation and experimental data.

```python
consist_data_file, consist_result_file = MS.LiP_XL_MS_ConsistencyTest()
```

**Returns:**
- Tuple containing paths to consistency data file (.npz) and results file (.xlsx)

**Process:**
1. **LiP-MS Analysis**: Compare simulated SASA with experimental LiP-MS signals
2. **XL-MS Analysis**: Validate cross-linking distances against experimental constraints
3. **Statistical Testing**: Perform permutation tests and bootstrap analysis
4. **Significance Assessment**: Calculate p-values and confidence intervals

##### select_rep_structs()
Selects representative structures with highest experimental consistency.

```python
MS.select_rep_structs(
    consist_data_file=consist_data_file,
    consist_result_file=consist_result_file, 
    total_traj_num_frames=335,
    n_analysis_frames=67
)
```

**Parameters:**
- `consist_data_file` (str): Path to consistency test data file
- `consist_result_file` (str): Path to consistency test results file
- `total_traj_num_frames` (int): Total trajectory frames
- `n_analysis_frames` (int): Frames to consider from trajectory end

**Output:**
- Representative structure files for each consistent state
- Visualization scripts and data files
- Summary statistics and rankings

## Usage Examples

### Basic Consistency Analysis

```python
from NCLEdetector.compare_sim2exp import MassSpec

# Initialize comprehensive analysis
MS = MassSpec(
    msm_data_file="MSM_data/msm_mapping.csv",
    meta_dist_file="distances/meta_distances.npz", 
    LiPMS_exp_file="experimental/LiPMS_data.xlsx",
    sasa_data_file="simulation/SASA_trajectories.npy",
    XLMS_exp_file="experimental/XLMS_constraints.xlsx",
    dist_data_file="simulation/distance_matrices.npy", 
    cluster_data_file="clustering/entanglement_clusters.npz",
    OPpath="order_parameters/",
    AAdcd_dir="trajectories/all_atom/",
    native_AA_pdb="structures/native.pdb",
    native_state_idx=0,
    state_idx_list=[0, 1, 2, 3, 4, 5],
    rm_traj_list=[65, 75, 155, 162],  # Exclude problematic trajectories
    outdir="results/consistency_analysis/",
    ID="protein_validation",
    verbose=True
)

# Load order parameter data
MS.load_OP()

# Perform consistency testing
print("Running consistency tests...")
consist_data, consist_results = MS.LiP_XL_MS_ConsistencyTest()

print(f"Consistency data saved to: {consist_data}")
print(f"Results summary saved to: {consist_results}")

# Select representative structures
MS.select_rep_structs(
    consist_data_file=consist_data,
    consist_result_file=consist_results,
    total_traj_num_frames=335, 
    n_analysis_frames=67
)
```

### Advanced Multi-State Analysis

```python
# Analyze specific metastable states
target_states = [1, 3, 5, 7]  # Specific MSM states of interest

MS_targeted = MassSpec(
    msm_data_file="msm_analysis/state_definitions.csv",
    meta_dist_file="analysis/meta_distances.npz",
    LiPMS_exp_file="experiments/condition1_LiPMS.xlsx", 
    sasa_data_file="simulations/sasa_per_state.npy",
    XLMS_exp_file="experiments/crosslinks_validation.xlsx",
    dist_data_file="simulations/distances_per_state.npy",
    cluster_data_file="clustering/state_clusters.npz",
    OPpath="order_params_by_state/",
    AAdcd_dir="trajectories/state_specific/", 
    native_AA_pdb="references/native_crystal.pdb",
    native_state_idx=0,
    state_idx_list=target_states,
    outdir="results/targeted_states/",
    ID="metastable_validation",
    num_perm=5000,  # Higher precision
    n_boot=1000,
    verbose=True
)

# Run targeted analysis
consist_data, consist_results = MS_targeted.LiP_XL_MS_ConsistencyTest()
MS_targeted.select_rep_structs(consist_data, consist_results, 500, 100)
```

## Integration with Scripts

### run_compare_sim2exp.py
Complete simulation-experiment validation pipeline:

```bash
python scripts/run_compare_sim2exp.py \
  --msm_data_file MSM_analysis/msm_states.csv \
  --meta_dist_file analysis/meta_distances.npz \
  --LiPMS_exp_file experiments/LiPMS_signals.xlsx \
  --sasa_data_file simulation/SASA_data.npy \
  --XLMS_exp_file experiments/crosslink_data.xlsx \
  --dist_data_file simulation/distance_data.npy \
  --cluster_data_file clustering/clusters.npz \
  --OPpath order_parameters/ \
  --AAdcd_dir trajectories/all_atom/ \
  --native_AA_pdb structures/native.pdb \
  --state_idx_list 1 2 3 4 5 \
  --prot_len 390 \
  --n_analysis_frames 100 \
  --rm_traj_list 65 75 155 \
  --native_state_idx 0 \
  --outdir results/comparison/ \
  --ID protein_study \
  --start 0 \
  --end 1000 \
  --stride 1 \
  --num_perm 1000 \
  --n_boot 100 \
  --lag_frame 20 \
  --nproc 10
```

**Workflow Steps:**
1. **Data Integration**: Load simulation and experimental datasets
2. **State Analysis**: Analyze each MSM state against experimental constraints  
3. **Consistency Testing**: Perform statistical validation
4. **Representative Selection**: Identify best-matching structures
5. **Visualization**: Generate analysis outputs and structure files

## Output Files

### Consistency Analysis Output
- `consistency_test_data.npz`: Raw consistency test data and statistics
- `consistency_results.xlsx`: Summary of consistency scores per state/signal
- `statistical_analysis.csv`: P-values, confidence intervals, effect sizes
- `permutation_test_results.npz`: Detailed permutation test outcomes

### Representative Structures Output
- `consist_signal_struct_data.npz`: Data linking consistent signals to structures
- `Consistent_structures_v8.xlsx`: Excel summary of representative structures
- `viz_rep_struct/`: Directory containing visualization files
  - `State_X/`: Per-state visualization data
  - `info.txt`: Detailed information for each representative structure
  - `*.pdb`: Representative structure coordinates

### Analysis Summaries
- `consistency_heatmap.png`: Visual summary of state-signal consistency
- `experimental_validation_report.pdf`: Comprehensive analysis report
- `structure_quality_metrics.csv`: Quality scores for selected structures

## Experimental Data Integration

### LiP-MS Data Format
Expected Excel format with columns:
- `Protein`: Protein identifier
- `Residue`: Residue number 
- `Signal_Change`: Experimental signal change
- `Significance`: Statistical significance level
- `Condition`: Experimental condition

### XL-MS Data Format  
Expected Excel format with columns:
- `Protein1`: First protein in crosslink
- `Residue1`: First residue position
- `Protein2`: Second protein in crosslink  
- `Residue2`: Second residue position
- `Distance`: Experimental distance constraint
- `Confidence`: Confidence level

## Statistical Methods

### Consistency Metrics
1. **Pearson Correlation**: Linear relationship between simulation and experiment
2. **Spearman Rank Correlation**: Monotonic relationship assessment
3. **Mutual Information**: Non-linear dependency measurement
4. **ROC Analysis**: Binary classification performance

### Significance Testing
1. **Permutation Tests**: Non-parametric significance assessment
2. **Bootstrap Confidence Intervals**: Uncertainty quantification
3. **Multiple Comparison Correction**: False discovery rate control
4. **Effect Size Calculation**: Practical significance measurement

### Quality Metrics
1. **Coverage**: Fraction of experimental signals matched
2. **Specificity**: False positive rate in predictions
3. **Sensitivity**: True positive rate in detection
4. **F1 Score**: Harmonic mean of precision and recall

## Performance Considerations

- **Memory Usage**: Scales with trajectory length and number of states
- **Computation Time**: 10 minutes to several hours depending on data size
- **Parallelization**: Benefit significantly from multiprocessing
- **Storage Requirements**: Large output files for comprehensive analyses

## Quality Control

### Data Validation
- **Input Format Checking**: Ensure proper data formats
- **Missing Data Handling**: Robust handling of incomplete experimental data
- **Outlier Detection**: Identify unusual experimental or simulation values
- **Temporal Consistency**: Validate trajectory continuity

### Statistical Robustness  
- **Cross-Validation**: Multiple validation approaches
- **Sensitivity Analysis**: Parameter robustness testing
- **Convergence Checking**: Ensure statistical convergence
- **Reproducibility**: Deterministic random number generation

## Theory Background

### LiP-MS Integration Theory
Limited Proteolysis-Mass Spectrometry measures protein flexibility by:
1. **Proteolytic Accessibility**: Flexible regions are more susceptible to cleavage
2. **SASA Correlation**: Simulation SASA correlates with experimental accessibility
3. **Dynamic Validation**: Compare simulation flexibility with experimental signals

### XL-MS Validation Theory  
Cross-linking Mass Spectrometry provides distance constraints:
1. **Distance Constraints**: Crosslinks impose maximum distance limits
2. **Structural Validation**: Simulation distances must satisfy experimental constraints
3. **Dynamic Constraints**: Time-averaged distances from simulation

## Dependencies

- **Core**: NumPy, SciPy, Pandas
- **Statistics**: Statsmodels, Scikit-learn
- **File I/O**: OpenPyXL (Excel files), MDAnalysis
- **Visualization**: Matplotlib, Seaborn
- **Parallel Processing**: Multiprocessing, Joblib

### Key Methods
- `load_OP(start=0, end=99999999999)`: Loads G and Q order parameter values for each trajectory into arrays for analysis.
- `LiP_XL_MS_ConsistencyTest()`: Performs consistency tests between simulation and experimental LiP-MS/XL-MS data. Returns paths to saved data/results files.
- `select_rep_structs(consist_data_file, consist_result_file, total_traj_num_frames, n_analysis_frames)`: After performing the consistency test, selects representative structures with high consistency between simulation and experiment. Generates output files and visualizations.

#### select_rep_structs
```python
select_rep_structs(
    consist_data_file: str,
    consist_result_file: str,
    total_traj_num_frames: int,
    n_analysis_frames: int
)
```
**Arguments:**
- `consist_data_file` (str): Path to the .npz file with consistency test data (output of `LiP_XL_MS_ConsistencyTest`).
- `consist_result_file` (str): Path to the Excel file with consistency test results (output of `LiP_XL_MS_ConsistencyTest`).
- `total_traj_num_frames` (int): Total number of frames in the trajectory.
- `n_analysis_frames` (int): Number of frames from the end of the trajectory to consider for representative selection.

**Description:**
Selects representative structures from the simulation that are most consistent with experimental signals, based on the results of the consistency test. Outputs include:
- `consist_signal_struct_data.npz`: Data on consistent signals and representative structures.
- `Consistent_structures_v8.xlsx`: Excel summary of representative structures per state and signal.
- `viz_rep_struct/`: Directory with visualization scripts and info files for each representative structure.

**Typical Usage:**
```python
from NCLEdetector.compare_sim2exp import MassSpec

MS = MassSpec(...)
consist_data_file, consist_result_file = MS.LiP_XL_MS_ConsistencyTest()
MS.select_rep_structs(
    consist_data_file, consist_result_file,
    total_traj_num_frames=335, n_analysis_frames=67
)
```

---

## Usage Example

```python
from NCLEdetector.compare_sim2exp import MassSpec

MS = MassSpec(
    msm_data_file='msm_data.csv',
    meta_dist_file='meta_dist.npz',
    LiPMS_exp_file='LiPMS.xlsx',
    sasa_data_file='SASA.npy',
    XLMS_exp_file='XLMS.xlsx',
    dist_data_file='Jwalk.npy',
    cluster_data_file='cluster_data.npz',
    OPpath='OP/',
    AAdcd_dir='AAtrajs/*.dcd',
    native_AA_pdb='native.pdb',
    native_state_idx=0,
    state_idx_list=[0,1,2],
    rm_traj_list=[],
    outdir='results/',
    ID='protein1',
    start=0,
    end=1000,
    stride=1,
    verbose=True
)
consist_data_file, consist_result_file = MS.LiP_XL_MS_ConsistencyTest()
MS.select_rep_structs(
    consist_data_file, consist_result_file,
    total_traj_num_frames=335, n_analysis_frames=67
)
```

**Notes:**
- The `compare_to_experiment` method does **not** exist; use `LiP_XL_MS_ConsistencyTest` and `select_rep_structs` for the main workflow.
- Output directories and files will be created as needed.
- For further details, refer to the source code in `NCLEdetector/compare_sim2exp.py` or the package documentation.

