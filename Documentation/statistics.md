# Statistics Module

The `statistics` module provides comprehensive statistical analysis tools for protein entanglement studies, enabling population-level modeling, Monte Carlo optimization, and folding pathway analysis across large proteomic datasets.

## Overview

This module enables:
- **Population modeling**: Logistic regression analysis across heterogeneous protein datasets
- **Monte Carlo optimization**: Identify protein groups with similar entanglement properties
- **Folding pathway analysis**: Statistical characterization of protein folding trajectories
- **Experimental correlation**: Statistical validation against experimental data
- **Proteome-wide studies**: Large-scale analysis across entire proteomes

## Classes

### ProteomeLogisticRegression

Performs logistic regression analysis to quantify associations between experimental signals and native entangled regions across heterogeneous protein datasets.

#### Initialization
```python
from NCLEdetector.statistics import ProteomeLogisticRegression

ProtRegression = ProteomeLogisticRegression(
    dataframe_files="proteome_data/",
    outdir="results/regression/",
    gene_list="gene_list.txt", 
    tag="proteome_study",
    reg_formula="cut_C_Rall ~ AA + region"
)
```

**Parameters:**
- `dataframe_files` (str): Directory containing protein feature dataframes (CSV format)
- `outdir` (str): Output directory for regression results
- `gene_list` (str): Path to file listing genes/proteins to include
- `tag` (str): Identifier for the analysis (used in output files)
- `reg_formula` (str): Regression formula specification

#### Key Methods

##### load_data()
Loads and preprocesses proteome-wide feature data.

```python
ProtRegression.load_data(
    sep='|',
    reg_var=['AA', 'region'],
    response_var='cut_C_Rall', 
    var2binarize=['cut_C_Rall', 'region'],
    mask_column='mapped_resid'
)
```

**Parameters:**
- `sep` (str): File separator for CSV files
- `reg_var` (list): Regression variables to include
- `response_var` (str): Response variable name
- `var2binarize` (list): Variables to convert to binary encoding
- `mask_column` (str): Column for filtering valid residues

##### run()
Executes complete regression analysis pipeline.

```python
regression_results = ProtRegression.run()
```

**Returns:**
- DataFrame containing regression coefficients, p-values, and statistics

### MonteCarlo

Performs Monte Carlo optimization to identify protein groups with poor per-protein statistics but significant population-level signals in experimental data.

#### Initialization
```python
from NCLEdetector.statistics import MonteCarlo

MC = MonteCarlo(
    dataframe_files="proteome_data/",
    outpath="results/monte_carlo/",
    gene_list="gene_list.txt",
    tag="mc_analysis", 
    reg_formula="cut_C_Rall ~ region + AA",
    response_var="cut_C_Rall",
    test_var="region",
    random=False,
    n_groups=4,
    steps=100000,
    C1=1.0,
    C2=2.5, 
    beta=0.05,
    linearT=False
)
```

**Parameters:**
- `dataframe_files` (str): Directory with protein feature files
- `outpath` (str): Output directory for Monte Carlo results
- `gene_list` (str): Path to gene/protein list file
- `tag` (str): Analysis identifier
- `reg_formula` (str): Regression formula for group analysis
- `response_var` (str): Response variable for regression
- `test_var` (str): Variable to test in regression
- `random` (bool): Whether to randomize response variable (control)
- `n_groups` (int): Number of protein groups for optimization
- `steps` (int): Number of Monte Carlo steps
- `C1` (float): Energy function coefficient (odds ratio term)
- `C2` (float): Energy function coefficient (size distribution term)
- `beta` (float): Initial temperature for simulated annealing
- `linearT` (bool): Use linear vs exponential temperature scaling

#### Key Methods

##### load_data()
Loads and preprocesses data for Monte Carlo analysis.

```python
MC.load_data(
    sep='|',
    reg_var=['AA', 'region'],
    response_var='cut_C_Rall',
    var2binarize=['cut_C_Rall', 'region'], 
    mask_column='mapped_resid',
    ID_column='gene',
    Length_column='uniprot_length'
)
```

##### run()
Executes Monte Carlo optimization with simulated annealing.

```python
MC.run(encoded_df=MC.data, ID_column='gene')
```

**Process:**
1. **Initialization**: Create random protein groups
2. **Energy calculation**: Compute regression metrics for each group
3. **Group swapping**: Propose moves between groups
4. **Acceptance criteria**: Accept/reject based on energy improvement
5. **Temperature annealing**: Gradually reduce temperature
6. **Convergence**: Monitor convergence and save optimal groupings

### FoldingPathwayStats

Analyzes statistical properties of protein folding pathways from temperature quenching simulations and MSM data.

#### Initialization
```python
from NCLEdetector.statistics import FoldingPathwayStats

FP = FoldingPathwayStats(
    msm_data=msm_dataframe,
    meta_set_file="meta_set.csv",
    traj_type_col="trajectory_type",
    outdir="results/folding_analysis/",
    rm_traj_list=[65, 75, 155, 162]
)
```

**Parameters:**
- `msm_data` (DataFrame): MSM trajectory data
- `meta_set_file` (str): Path to trajectory metadata file
- `traj_type_col` (str): Column specifying trajectory classification
- `outdir` (str): Output directory for pathway analysis
- `rm_traj_list` (list): Trajectory indices to exclude

#### Key Methods

##### analyze_pathways()
Performs comprehensive folding pathway analysis.

```python
pathway_stats = FP.analyze_pathways()
```

**Analysis includes:**
- Pathway probability distributions
- Transition state analysis
- Folding time statistics
- Entanglement evolution patterns

## Usage Examples

### Population-Level Logistic Regression

```python
from NCLEdetector.statistics import ProteomeLogisticRegression

# Initialize regression analysis for proteome study
ProtRegression = ProteomeLogisticRegression(
    dataframe_files="proteome_features/",
    outdir="results/population_modeling/", 
    gene_list="ecoli_proteins.txt",
    tag="ecoli_proteome",
    reg_formula="cut_C_Rall ~ AA + region"
)

# Load and preprocess data
ProtRegression.load_data(
    sep='|',
    reg_var=['AA', 'region'], 
    response_var='cut_C_Rall',
    var2binarize=['cut_C_Rall', 'region'],
    mask_column='mapped_resid'
)

# Run regression analysis
results = ProtRegression.run()

# Display key results
print(f"Regression coefficients:")
print(results[['Variable', 'Coefficient', 'P_value', 'Odds_Ratio']])

# Save results
results.to_csv("population_regression_results.csv", index=False)
```

### Monte Carlo Group Optimization

```python
from NCLEdetector.statistics import MonteCarlo

# Set up Monte Carlo optimization
MC = MonteCarlo(
    dataframe_files="proteome_features/",
    outpath="results/monte_carlo_analysis/",
    gene_list="candidate_proteins.txt", 
    tag="entanglement_groups",
    reg_formula="cut_C_Rall ~ region + AA",
    response_var="cut_C_Rall",
    test_var="region",
    n_groups=4,
    steps=50000,
    C1=1.0,
    C2=2.5,
    beta=0.05
)

# Load data
MC.load_data(
    sep='|',
    reg_var=['AA', 'region'],
    response_var='cut_C_Rall', 
    var2binarize=['cut_C_Rall', 'region'],
    mask_column='mapped_resid',
    ID_column='gene',
    Length_column='uniprot_length'
)

print("Starting Monte Carlo optimization...")
MC.run(encoded_df=MC.data, ID_column='gene')

print("Optimization complete! Check output directory for results.")
```

### Folding Pathway Analysis

```python
import pandas as pd
from NCLEdetector.statistics import FoldingPathwayStats

# Load MSM trajectory data
msm_data = pd.read_csv("MSM_analysis/trajectory_data.csv")

# Initialize folding pathway analysis
FP = FoldingPathwayStats(
    msm_data=msm_data,
    meta_set_file="MSM_analysis/meta_set.csv",
    traj_type_col="folding_type",
    outdir="results/pathway_analysis/",
    rm_traj_list=[65, 75, 155, 162, 199]  # Remove problematic trajectories
)

# Analyze folding pathways
pathway_results = FP.analyze_pathways()

print("Folding pathway analysis complete!")
print(f"Identified {len(pathway_results['pathways'])} distinct pathways")
```

## Integration with Scripts

### run_population_modeling.py
Population-level logistic regression analysis:

```bash
python scripts/run_population_modeling.py \
  --dataframe_files proteome_data/ \
  --outdir results/population_analysis/ \
  --gene_list protein_list.txt \
  --tag proteome_regression \
  --reg_formula "cut_C_Rall ~ AA + region"
```

### run_montecarlo.py
Monte Carlo optimization for protein grouping:

```bash
python scripts/run_montecarlo.py \
  --dataframe_files proteome_data/ \
  --outpath results/monte_carlo/ \
  --gene_list candidate_genes.txt \
  --tag mc_optimization \
  --steps 100000 \
  --n_groups 4 \
  --C1 1.0 \
  --C2 2.5 \
  --beta 0.05
```

### run_Foldingpathway.py
Statistical analysis of folding pathways:

```bash
python scripts/run_Foldingpathway.py \
  --msm_data_file MSM_data/trajectory_meta.csv \
  --meta_set_file MSM_data/meta_set.csv \
  --traj_type_col trajectory_classification \
  --outdir results/pathway_stats/ \
  --rm_traj_list 65 75 155
```

## Output Files

### Population Modeling Output
- `regression_results_<tag>.csv`: Complete regression analysis results
- `coefficient_summary.csv`: Summary of significant coefficients  
- `odds_ratios.csv`: Odds ratios with confidence intervals
- `model_diagnostics.txt`: Regression model validation metrics

### Monte Carlo Output
- `Final_step_reg_<tag>.csv`: Final regression results for optimal groups
- `State<X>_final_genelist_<tag>.txt`: Protein lists for each optimized group
- `<tag>.log`: Detailed optimization progress log
- `energy_trajectory.csv`: Energy function evolution during optimization
- `group_assignments.csv`: Final group assignments for all proteins

### Folding Pathway Output
- `pathway_statistics.csv`: Statistical summary of identified pathways
- `transition_probabilities.csv`: State transition probability matrices
- `folding_times.csv`: Distribution of folding times per pathway
- `entanglement_evolution.csv`: Entanglement changes along pathways

### Visualization Output
- `regression_coefficients.png`: Coefficient plots with confidence intervals
- `monte_carlo_convergence.png`: Optimization convergence plots  
- `pathway_network.png`: Network representation of folding pathways
- `correlation_heatmap.png`: Feature correlation matrices

## Statistical Methods

### Logistic Regression Framework
**Model specification:**
```
logit(P(Y=1)) = β₀ + β₁X₁ + β₂X₂ + ... + βₚXₚ
```

Where:
- Y: Binary response variable (e.g., experimental signal change)
- X: Predictor variables (amino acid type, entanglement region)
- β: Regression coefficients

**Key metrics:**
- **Odds Ratios**: exp(β) - relative odds of outcome
- **P-values**: Statistical significance of coefficients
- **Confidence Intervals**: Uncertainty quantification
- **AIC/BIC**: Model selection criteria

### Monte Carlo Optimization
**Energy function:**
```
E = -C₁ * log(OR) + C₂ * Σ|nᵢ - n̄|
```

Where:
- OR: Odds ratio from group regression
- nᵢ: Size of group i
- n̄: Mean group size
- C₁, C₂: Weighting parameters

**Simulated annealing:**
- **Temperature schedule**: T(t) = T₀ * exp(-t/τ) or linear
- **Acceptance probability**: P = exp(-ΔE/T)
- **Convergence criteria**: Energy stabilization

### Pathway Analysis Methods
- **Transition state theory**: First passage time analysis
- **Network analysis**: Graph-based pathway identification
- **Bootstrap sampling**: Uncertainty estimation
- **Clustering algorithms**: Pathway classification

## Quality Control

### Model Validation
- **Cross-validation**: K-fold validation for regression models
- **Residual analysis**: Model assumption checking
- **Outlier detection**: Identify influential observations
- **Multicollinearity**: Variance inflation factor analysis

### Monte Carlo Validation
- **Convergence diagnostics**: Energy trajectory analysis
- **Multiple runs**: Reproducibility assessment
- **Control experiments**: Randomization tests
- **Parameter sensitivity**: Robustness analysis

### Statistical Robustness
- **Multiple comparison correction**: False discovery rate control
- **Effect size calculation**: Practical significance assessment
- **Bootstrap confidence intervals**: Non-parametric uncertainty
- **Permutation tests**: Distribution-free significance testing

## Performance Considerations

- **Memory Usage**: Scales with dataset size and number of proteins
- **Computation Time**: 
  - Population modeling: ~10-60 minutes
  - Monte Carlo: ~1-24 hours depending on steps
  - Pathway analysis: ~10-120 minutes
- **Parallelization**: Limited parallelization in current implementation
- **Convergence**: Monitor for proper convergence in iterative methods

## Applications

### Proteome-Wide Studies
1. **Entanglement mapping**: Identify entangled proteins across proteomes
2. **Functional analysis**: Correlate entanglement with protein function
3. **Evolutionary studies**: Analyze entanglement conservation
4. **Disease association**: Link entanglement to pathological states

### Drug Discovery
1. **Target identification**: Find entangled proteins in disease pathways
2. **Allosteric sites**: Identify regulatory sites in entangled regions
3. **Drug mechanism**: Understand how drugs affect protein entanglement
4. **Side effect prediction**: Predict off-target effects

### Biotechnology Applications
1. **Protein design**: Design proteins with desired entanglement properties
2. **Enzyme engineering**: Optimize catalytic sites considering entanglement
3. **Stability engineering**: Enhance protein stability through entanglement
4. **Folding optimization**: Design proteins with efficient folding pathways

## Theory Background

### Population Modeling Theory
Population-level analysis enables detection of weak signals that may not be significant at the individual protein level but show consistent trends across many proteins.

**Statistical power**: N proteins increases power to detect small effect sizes
**Heterogeneity modeling**: Account for protein-specific variation
**Meta-analysis principles**: Combine evidence across multiple proteins

### Monte Carlo Theory
Simulated annealing enables escape from local optima in the optimization landscape, allowing discovery of globally optimal protein groupings.

**Metropolis-Hastings algorithm**: Balanced detailed balance
**Simulated annealing**: Temperature-controlled exploration
**Energy landscape**: Complex multi-modal optimization surface

## Dependencies

- **Core**: NumPy, SciPy, Pandas
- **Statistics**: Statsmodels, Scikit-learn
- **Optimization**: SciPy optimization, custom algorithms  
- **Visualization**: Matplotlib, Seaborn
- **Data handling**: Pickle, CSV readers
- **Parallel processing**: Multiprocessing (limited usage)

