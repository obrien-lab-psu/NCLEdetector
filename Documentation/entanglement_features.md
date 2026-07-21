# Entanglement Features Module

The `entanglement_features` module provides comprehensive functionality for extracting and calculating structural features from clustered native entanglements in protein structures, enabling downstream statistical analysis and machine learning applications.

## Overview

This module enables:
- **Feature extraction**: Calculate geometric and topological features from entanglements
- **Structural analysis**: Analyze loop properties and entanglement characteristics
- **Data preparation**: Generate feature matrices for statistical modeling
- **Cross-protein comparison**: Standardized feature sets across different proteins

## Classes

### FeatureGen

Main class for generating entanglement-related structural features from clustered native entanglement data.

#### Initialization
```python
from NCLEdetector.entanglement_features import FeatureGen

FGen = FeatureGen(
    pdb_file="structure.pdb",
    outdir="results/features/", 
    cluster_file="clustered_entanglements.csv"
)
```

**Parameters:**
- `pdb_file` (str): Path to cleaned PDB structure file
- `outdir` (str): Output directory for feature results  
- `cluster_file` (str): Path to clustered native entanglement data file

#### Key Methods

##### get_uent_features()
Calculates comprehensive features for unique entanglements in the protein structure.

```python
features = FGen.get_uent_features(
    gene="P00558",
    chain="A", 
    pdbid="1ZMR"
)
```

**Parameters:**
- `gene` (str): UniProt ID of the protein
- `chain` (str): Chain identifier in the PDB file
- `pdbid` (str): PDB ID of the structure

**Returns:**
- Dictionary/DataFrame containing calculated entanglement features

## Feature Categories

### Geometric Features
- **Loop lengths**: Length of entangled loop segments
- **Loop distances**: Spatial separation between loop centers
- **Loop angles**: Angular relationships between entangled loops
- **Contact surfaces**: Interfacial areas between entangled regions

### Topological Features  
- **Linking numbers**: Gaussian linking number values
- **Crossing patterns**: Specific crossing configurations
- **Entanglement complexity**: Number of crossings and their arrangement
- **Knot invariants**: Topological descriptors of entanglement type

### Structural Features
- **Secondary structure**: Helix, sheet, loop content in entangled regions
- **Residue composition**: Amino acid frequencies in entangled loops
- **Hydrophobicity patterns**: Hydrophobic/hydrophilic distribution
- **Flexibility measures**: B-factors and predicted flexibility

### Evolutionary Features
- **Conservation scores**: Evolutionary conservation of entangled residues
- **Organism-specific patterns**: Taxonomic entanglement preferences
- **Functional annotations**: GO terms and domain associations

## Usage Examples

### Basic Feature Generation

```python
from NCLEdetector.entanglement_features import FeatureGen

# Initialize feature generator
FGen = FeatureGen(
    pdb_file="1zmr_clean.pdb",
    outdir="results/entanglement_features/",
    cluster_file="results/clustered/1zmr_clustered.csv"
)

# Generate comprehensive feature set
features = FGen.get_uent_features(
    gene="P00558",  # Uniprot ID for protein
    chain="A",      # Chain identifier
    pdbid="1ZMR"    # PDB code
)

# Display feature summary
print(f"Generated {len(features)} features")
print(f"Feature categories: {list(features.keys())}")
```

### Multi-Protein Feature Analysis

```python
import pandas as pd
from NCLEdetector.entanglement_features import FeatureGen

# Process multiple proteins
proteins = [
    {"pdb": "1zmr.pdb", "cluster": "1zmr_clustered.csv", "gene": "P00558", "chain": "A", "pdbid": "1ZMR"},
    {"pdb": "2xyz.pdb", "cluster": "2xyz_clustered.csv", "gene": "P12345", "chain": "A", "pdbid": "2XYZ"}
]

all_features = []
for protein in proteins:
    FGen = FeatureGen(
        pdb_file=protein["pdb"],
        outdir=f"results/features/{protein['pdbid']}/",
        cluster_file=protein["cluster"]
    )
    
    features = FGen.get_uent_features(
        gene=protein["gene"],
        chain=protein["chain"], 
        pdbid=protein["pdbid"]
    )
    
    # Add protein identifier
    features["protein_id"] = protein["pdbid"]
    all_features.append(features)

# Combine into analysis-ready dataset
feature_matrix = pd.concat(all_features, ignore_index=True)
```

## Integration with Scripts

### run_nativeNCLE.py
Complete pipeline including feature generation:

```bash
python scripts/run_nativeNCLE.py \
  --struct structure.pdb \
  --outdir results/native_analysis/ \
  --ID protein_name \
  --organism Ecoli
```

This script automatically:
1. Calculates native entanglements
2. Filters for high quality entanglements  
3. Clusters entanglements to remove degeneracies
4. **Generates entanglement features** using FeatureGen

The script demonstrates the complete workflow:

```python
# Generate entanglement features for clustered native entanglements
FGen = FeatureGen(struct, outdir=os.path.join(outdir, "Native_clustered_HQ_GE_features"), cluster_file=nativeClusteredEnt['outfile'])
EntFeatures = FGen.get_uent_features(gene='P00558', chain='A', pdbid='1ZMR')
```

## Output Files

### Feature Data Files
- `*_entanglement_features.csv`: Comprehensive feature matrix
- `*_feature_metadata.json`: Feature definitions and descriptions
- `*_geometric_features.csv`: Geometric measurements only
- `*_topological_features.csv`: Topological descriptors only

### Analysis Files
- `feature_correlation_matrix.csv`: Inter-feature correlations
- `feature_importance_scores.csv`: Relative importance rankings
- `entanglement_classification.csv`: Entanglement type predictions

### Visualization Files
- `feature_distributions.png`: Histograms of feature values
- `entanglement_map.png`: 2D projection of entanglement space
- `correlation_heatmap.png`: Feature correlation visualization

## Feature Definitions

### Core Features (Always Calculated)
| Feature | Description | Units |
|---------|-------------|--------|
| `loop1_length` | Length of first entangled loop | residues |
| `loop2_length` | Length of second entangled loop | residues |
| `gln_value` | Gaussian linking number | dimensionless |
| `loop_separation` | Distance between loop centers | Angstroms |
| `crossing_angle` | Angle between crossing strands | degrees |

### Extended Features (Optional)
| Feature | Description | Units |
|---------|-------------|--------|
| `hydrophobicity_score` | Average hydrophobicity | Kyte-Doolittle scale |
| `conservation_score` | Evolutionary conservation | 0-1 scale |
| `flexibility_index` | Predicted loop flexibility | B-factor units |
| `contact_density` | Contacts per residue | contacts/residue |

## Machine Learning Integration

### Feature Preprocessing

```python
# Prepare features for ML analysis
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest

# Load feature data
features_df = pd.read_csv("entanglement_features.csv")

# Separate features and target
X = features_df.drop(['protein_id', 'entanglement_id'], axis=1)
y = features_df['entanglement_type']  # if available

# Standardize features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Feature selection
selector = SelectKBest(k=20)
X_selected = selector.fit_transform(X_scaled, y)
```

### Classification Example

```python
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score

# Train entanglement classifier
clf = RandomForestClassifier(n_estimators=100, random_state=42)
scores = cross_val_score(clf, X_selected, y, cv=5)

print(f"Classification accuracy: {scores.mean():.3f} ± {scores.std():.3f}")
```

## Performance Considerations

- **Computation Time**: ~1-5 minutes per protein structure
- **Memory Usage**: Scales linearly with number of entanglements  
- **Feature Count**: Typically 50-200 features per entanglement
- **Scalability**: Suitable for proteome-wide analysis

## Quality Control

### Feature Validation
- **Range checking**: Ensure features are within expected bounds
- **Correlation analysis**: Identify redundant features
- **Missing value handling**: Robust handling of incomplete data
- **Outlier detection**: Identify unusual entanglement configurations

### Reproducibility
- **Deterministic calculations**: Consistent results across runs
- **Version tracking**: Feature calculation method versioning
- **Parameter logging**: Complete parameter documentation

## Dependencies

- **Core**: NumPy, Pandas, SciPy
- **Structural**: MDAnalysis, ProDy  
- **Analysis**: Scikit-learn (optional for ML features)
- **Visualization**: Matplotlib, Seaborn (optional)

