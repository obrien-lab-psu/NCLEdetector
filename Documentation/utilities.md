# Utilities Module

The `utilities` module provides essential preprocessing and data manipulation tools for protein structure analysis, specializing in PDB file cleaning, trajectory processing, and data format conversions required for entanglement analysis workflows.

## Overview

This module enables:
- **Structure preprocessing**: Clean and validate PDB files for analysis
- **Trajectory processing**: Handle molecular dynamics trajectory data
- **Format conversions**: Convert between different structural data formats
- **Quality control**: Validate structural data integrity
- **Data standardization**: Ensure consistent data formatting across workflows

## Classes

### PDBcleaner

Comprehensive PDB file preprocessing and cleaning utilities for removing duplicate residues, handling missing atoms, and preparing structures for entanglement analysis.

#### Initialization
```python
from NCLEdetector.utilities import PDBcleaner

cleaner = PDBcleaner(
    pdb_file="protein_structure.pdb",
    outdir="cleaned_structures/",
    chain_selection=None,
    keep_hetero=False,
    remove_water=True,
    validate_structure=True
)
```

**Parameters:**
- `pdb_file` (str): Path to input PDB file requiring cleaning
- `outdir` (str, optional): Output directory for processed files. Default: current directory
- `chain_selection` (list, optional): Specific chains to retain (e.g., ['A', 'B']). Default: all chains
- `keep_hetero` (bool, optional): Retain heteroatoms and ligands. Default: False
- `remove_water` (bool, optional): Remove water molecules. Default: True
- `validate_structure` (bool, optional): Perform structure validation. Default: True

#### Key Methods

##### remove_duplicates()
Identifies and removes duplicate residues while preserving structural integrity.

```python
cleaned_pdb = cleaner.remove_duplicates(
    tolerance=0.1,
    preserve_occupancy=True,
    log_removed=True
)
```

**Parameters:**
- `tolerance` (float): Distance tolerance for duplicate detection (Å)
- `preserve_occupancy` (bool): Keep highest occupancy residue for duplicates
- `log_removed` (bool): Log details of removed residues

**Returns:**
- `clean_pdb` (str): Path to cleaned PDB file

##### validate_structure()
Performs comprehensive structure validation and quality assessment.

```python
validation_report = cleaner.validate_structure()
```

**Returns:**
- Dictionary containing validation metrics:
  - `missing_residues`: List of missing residues in sequence
  - `chain_breaks`: Detected chain discontinuities
  - `unusual_bonds`: Potentially problematic bond lengths
  - `quality_score`: Overall structure quality assessment

##### standardize_format()
Standardizes PDB format for consistent downstream processing.

```python
standardized_pdb = cleaner.standardize_format(
    renumber_residues=True,
    fix_nomenclature=True,
    add_missing_atoms=False
)
```

**Parameters:**
- `renumber_residues` (bool): Renumber residues sequentially
- `fix_nomenclature` (bool): Standardize atom and residue names
- `add_missing_atoms` (bool): Add missing backbone atoms

##### extract_chains()
Extract specific protein chains for focused analysis.

```python
chain_files = cleaner.extract_chains(
    chains=['A', 'B'],
    separate_files=True,
    suffix="_chain"
)
```

### TrajectoryProcessor

Handles molecular dynamics trajectory preprocessing and frame extraction for entanglement analysis.

#### Initialization
```python
from NCLEdetector.utilities import TrajectoryProcessor

processor = TrajectoryProcessor(
    topology_file="system.pdb",
    trajectory_file="trajectory.xtc",
    outdir="processed_trajectories/",
    selection="protein"
)
```

**Parameters:**
- `topology_file` (str): Reference topology/structure file
- `trajectory_file` (str): Molecular dynamics trajectory file
- `outdir` (str): Output directory for processed frames
- `selection` (str): Atom selection string for processing

#### Key Methods

##### extract_frames()
Extract specific frames from trajectory for analysis.

```python
frame_files = processor.extract_frames(
    frame_indices=[0, 100, 200, 500, 1000],
    format="pdb",
    align_to_reference=True
)
```

##### process_trajectory()
Complete trajectory processing pipeline.

```python
processed_data = processor.process_trajectory(
    stride=10,
    align=True,
    center=True,
    output_format="pdb"
)
```

### DataConverter

Utilities for converting between different structural and analysis data formats.

#### Initialization
```python
from NCLEdetector.utilities import DataConverter

converter = DataConverter(
    input_format="pdb",
    output_format="xyz",
    precision=3
)
```

#### Key Methods

##### convert_structure()
Convert between structural file formats.

```python
converted_file = converter.convert_structure(
    input_file="structure.pdb",
    output_file="structure.xyz"
)
```

##### convert_analysis_data()
Convert analysis results between formats.

```python
converted_data = converter.convert_analysis_data(
    data=analysis_results,
    from_format="pickle",
    to_format="json"
)
```

## Usage Examples

### Basic PDB Cleaning

```python
from NCLEdetector.utilities import PDBcleaner

# Initialize cleaner for problematic PDB file
cleaner = PDBcleaner(
    pdb_file="raw_structure.pdb",
    outdir="cleaned_pdbs/",
    remove_water=True,
    keep_hetero=False
)

# Remove duplicate residues
print("Cleaning PDB file...")
cleaned_pdb = cleaner.remove_duplicates(
    tolerance=0.1,
    preserve_occupancy=True,
    log_removed=True
)

# Validate cleaned structure
print("Validating structure...")
validation = cleaner.validate_structure()
print(f"Structure quality score: {validation['quality_score']}")

# Standardize format for analysis
print("Standardizing format...")
standardized_pdb = cleaner.standardize_format(
    renumber_residues=True,
    fix_nomenclature=True
)

print(f"Cleaned structure saved to: {standardized_pdb}")
```

### Comprehensive Preprocessing Pipeline

```python
from NCLEdetector.utilities import PDBcleaner, DataConverter
import os

def preprocess_pdb_collection(pdb_directory, output_directory):
    """
    Comprehensive preprocessing pipeline for multiple PDB files.
    """
    pdb_files = [f for f in os.listdir(pdb_directory) if f.endswith('.pdb')]
    
    processed_files = []
    
    for pdb_file in pdb_files:
        print(f"Processing {pdb_file}...")
        
        # Initialize cleaner
        cleaner = PDBcleaner(
            pdb_file=os.path.join(pdb_directory, pdb_file),
            outdir=output_directory,
            remove_water=True,
            validate_structure=True
        )
        
        # Clean and validate
        try:
            cleaned_pdb = cleaner.remove_duplicates()
            validation = cleaner.validate_structure()
            
            # Only proceed if structure quality is acceptable
            if validation['quality_score'] > 0.7:
                standardized = cleaner.standardize_format()
                processed_files.append({
                    'original': pdb_file,
                    'cleaned': standardized,
                    'quality': validation['quality_score']
                })
                print(f"  ✓ Successfully processed (quality: {validation['quality_score']:.2f})")
            else:
                print(f"  ✗ Structure quality too low: {validation['quality_score']:.2f}")
                
        except Exception as e:
            print(f"  ✗ Error processing {pdb_file}: {e}")
    
    return processed_files

# Run preprocessing pipeline
processed_structures = preprocess_pdb_collection(
    pdb_directory="raw_pdbs/",
    output_directory="processed_pdbs/"
)

print(f"\nProcessed {len(processed_structures)} structures successfully")
```

### Trajectory Processing for Entanglement Analysis

```python
from NCLEdetector.utilities import TrajectoryProcessor
from NCLEdetector.gaussian_entanglement import GaussianEntanglement

# Process MD trajectory for entanglement analysis
processor = TrajectoryProcessor(
    topology_file="system.pdb", 
    trajectory_file="production.xtc",
    outdir="trajectory_frames/",
    selection="protein and name CA"
)

# Extract representative frames
print("Extracting trajectory frames...")
frame_files = processor.extract_frames(
    frame_indices=list(range(0, 10000, 100)),  # Every 100 frames
    format="pdb",
    align_to_reference=True
)

# Analyze entanglement evolution
entanglement_evolution = []

for i, frame_file in enumerate(frame_files):
    print(f"Analyzing frame {i+1}/{len(frame_files)}")
    
    # Calculate entanglement for this frame
    ge = GaussianEntanglement(frame_file)
    ge.calculate_entanglement()
    
    entanglement_evolution.append({
        'frame': i * 100,
        'time_ns': i * 100 * 0.002,  # 2 ps timestep
        'entanglement_complexity': ge.entanglement_complexity,
        'linking_number': ge.linking_number
    })

print("Entanglement evolution analysis complete!")
```

### Multi-Format Data Conversion

```python
from NCLEdetector.utilities import DataConverter
import pickle
import json

# Convert analysis results between formats
converter = DataConverter()

# Load pickle data
with open('analysis_results.pkl', 'rb') as f:
    analysis_data = pickle.load(f)

# Convert to JSON for web visualization
json_data = converter.convert_analysis_data(
    data=analysis_data,
    from_format="pickle", 
    to_format="json"
)

# Save JSON output
with open('analysis_results.json', 'w') as f:
    json.dump(json_data, f, indent=2)

# Convert PDB collection to XYZ format for external analysis
pdb_files = ["structure1.pdb", "structure2.pdb", "structure3.pdb"]

for pdb_file in pdb_files:
    xyz_file = converter.convert_structure(
        input_file=pdb_file,
        output_file=pdb_file.replace('.pdb', '.xyz')
    )
    print(f"Converted {pdb_file} -> {xyz_file}")
```

## Integration with Scripts

The utilities module integrates seamlessly with all NCLEdetector scripts for preprocessing and data preparation:

### Script Integration Examples

```bash
# Clean PDB files before entanglement analysis
python scripts/run_entanglement_identification.py \
  --pdb_file cleaned_structures/protein_cleaned.pdb \
  --outdir results/entanglement_analysis/

# Process trajectory frames for order parameter analysis  
python scripts/run_OP_CGtrajAnal.py \
  --trajectory_dir trajectory_frames/ \
  --topology system.pdb \
  --outdir results/order_parameters/

# Use processed structures for clustering analysis
python scripts/run_nonnative_entanglement_clustering.py \
  --input_dir processed_pdbs/ \
  --outdir results/clustering_analysis/
```

## Output Files

### PDB Cleaning Output
- `<filename>_cleaned.pdb`: Main cleaned PDB file
- `<filename>_removed_residues.log`: Log of removed duplicate residues
- `<filename>_validation_report.txt`: Structure validation summary
- `<filename>_standardized.pdb`: Format-standardized structure

### Trajectory Processing Output
- `frame_<N>.pdb`: Individual trajectory frames
- `trajectory_summary.csv`: Frame extraction metadata
- `alignment_rmsd.csv`: Structural alignment quality metrics
- `processed_trajectory.xtc`: Cleaned and aligned trajectory

### Data Conversion Output
- Format-specific output files (`.xyz`, `.json`, `.csv`, etc.)
- `conversion_log.txt`: Conversion process summary
- `format_validation.txt`: Output format validation report

### Quality Control Output
- `structure_quality_report.csv`: Comprehensive quality metrics
- `preprocessing_summary.txt`: Complete preprocessing pipeline summary
- `error_log.txt`: Any errors or warnings encountered

## Quality Control Features

### Structure Validation
- **Geometric validation**: Bond lengths, angles, and chirality checks
- **Sequence validation**: Missing residues and chain break detection
- **Coordinate validation**: Reasonable coordinate ranges and precision
- **Format validation**: PDB format compliance checking

### Data Integrity Checks
- **File format verification**: Ensure correct format parsing
- **Coordinate consistency**: Validate coordinate transformations
- **Metadata preservation**: Maintain important structural annotations
- **Error logging**: Comprehensive error tracking and reporting

### Best Practices
- **Backup original files**: Always preserve original data
- **Validate output**: Check cleaned structures before analysis
- **Monitor quality metrics**: Track validation scores
- **Document processing**: Maintain preprocessing logs

## Performance Considerations

### Memory Usage
- **Large structures**: Efficient handling of structures >10,000 residues
- **Trajectory processing**: Chunked processing for large trajectories
- **Memory optimization**: Minimal memory footprint for batch processing

### Processing Speed
- **Parallel processing**: Multi-core support for batch operations
- **Efficient algorithms**: Optimized duplicate detection and cleaning
- **Caching**: Intermediate result caching for repeated operations

### Scalability
- **Batch processing**: Handle hundreds to thousands of structures
- **Pipeline integration**: Seamless workflow integration
- **Resource management**: Automatic cleanup and resource management

## Applications

### High-Throughput Studies
1. **Proteome-wide analysis**: Clean entire structural proteomes
2. **Database preparation**: Standardize PDB collections
3. **Quality filtering**: Remove low-quality structures
4. **Format standardization**: Ensure consistent data formats

### Molecular Dynamics Analysis
1. **Trajectory preprocessing**: Clean and align MD trajectories  
2. **Frame extraction**: Select representative conformations
3. **Structural evolution**: Track conformational changes
4. **Quality assessment**: Monitor simulation quality

### Structural Bioinformatics
1. **Database integration**: Prepare structures for analysis pipelines
2. **Comparative studies**: Standardize structures for comparison
3. **Method development**: Clean datasets for algorithm testing
4. **Validation studies**: Ensure data quality for benchmarking

## Theory Background

### Structure Preprocessing Theory
**Duplicate Detection**: Uses spatial proximity and chemical identity to identify redundant residues while preserving biological relevance.

**Geometric Validation**: Implements standard geometric criteria for protein structure validation including:
- Bond length validation (1.2-2.0 Å for common bonds)
- Bond angle validation (90-180° for common angles)  
- Chirality checking for amino acid stereochemistry

**Format Standardization**: Ensures compliance with PDB format specifications while maintaining structural information integrity.

### Quality Assessment Methods
**Structure Quality Metrics**:
- **Completeness**: Fraction of expected residues present
- **Geometric quality**: Bond length and angle deviations
- **Chemical quality**: Proper atom nomenclature and connectivity
- **Overall score**: Weighted combination of individual metrics

**Validation Algorithms**:
- **RMSD calculations**: Root mean square deviation for alignment quality
- **Ramachandran analysis**: Backbone conformation validation
- **Clash detection**: Steric conflict identification
- **B-factor analysis**: Temperature factor reasonableness

## Dependencies

- **Core**: NumPy, SciPy, Pandas
- **Structure handling**: BioPython, MDAnalysis (optional)
- **File I/O**: Standard library modules, format-specific parsers
- **Validation**: Custom geometric validation algorithms
- **Parallel processing**: Multiprocessing, concurrent.futures
- **Logging**: Standard logging framework

