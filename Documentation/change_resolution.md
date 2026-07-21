# Resolution Conversion Module

The `change_resolution` module provides comprehensive functionality for converting between all-atom and coarse-grained protein representations, enabling multi-scale modeling and simulation workflows with seamless resolution transitions.

## Overview

This module enables:
- **Coarse-graining**: Convert all-atom structures to coarse-grained representations
- **Back-mapping**: Reconstruct all-atom models from coarse-grained structures
- **Force field generation**: Create OpenMM-compatible force fields for coarse-grained models
- **Energy minimization**: Optimize structures at both resolutions
- **Multi-scale workflows**: Seamlessly transition between resolution scales

## Classes

### CoarseGrain

Converts all-atom protein structures to coarse-grained representations with customizable scaling factors and generates corresponding force field parameters.

#### Initialization
```python
from NCLEdetector.change_resolution import CoarseGrain

CoarseGrainer = CoarseGrain(
    pdbfile="protein.pdb",
    ID="protein_name", 
    nscal=2.0,
    outdir="results/coarse_graining/",
    fnn=1,
    potential_name='bt',
    casm=0,
    domain_file="domain_def.dat",
    ca_prefix='A',
    sc_prefix='B'
)
```

**Parameters:**
- `pdbfile` (str): Path to all-atom PDB structure file
- `ID` (str): Identifier for coarse-graining process
- `nscal` (float): Coarse-graining scale factor (default: 1.5)
- `outdir` (str): Output directory for results (default: './')
- `fnn` (int): Number of nearest neighbors for interactions (default: 1)
- `potential_name` (str): Interaction potential type (default: 'bt')
- `casm` (int): Use Cα-Cβ model if 1, Cα-only if 0 (default: 0)
- `domain_file` (str): Path to domain definition file (default: 'None')
- `ca_prefix` (str): Prefix for backbone atoms (default: 'A')
- `sc_prefix` (str): Prefix for sidechain atoms (default: 'B')

#### Key Methods

##### run()
Executes complete coarse-graining workflow.

```python
CG_files = CoarseGrainer.run()
```

**Returns:**
- Dictionary containing paths to generated files:
  - `'cor'`: Coarse-grained coordinates file
  - `'prm'`: Parameter file
  - `'psf'`: Protein structure file
  - `'top'`: Topology file

##### parse_cg_prm()
Generates OpenMM-compatible XML force field from parameter files.

```python
CoarseGrainer.parse_cg_prm(
    prmfile=CG_files['prm'], 
    topfile=CG_files['top']
)
```

**Parameters:**
- `prmfile` (str): Path to parameter file
- `topfile` (str): Path to topology file

### BackMapping

Reconstructs all-atom structures from coarse-grained representations using reference templates and energy minimization.

#### Initialization
```python
from NCLEdetector.change_resolution import BackMapping

backMapper = BackMapping(
    nproc=4,
    outdir="results/back_mapping/"
)
```

**Parameters:**
- `nproc` (int): Number of parallel processes (default: 1)
- `outdir` (str): Output directory for results (default: './')

#### Key Methods

##### backmap()
Reconstructs all-atom structure from coarse-grained coordinates.

```python
backMapper.backmap(
    cg_pdb="cg_structure.pdb",
    aa_pdb="reference_structure.pdb",
    TAG="protein_backmapped"
)
```

**Parameters:**
- `cg_pdb` (str): Path to coarse-grained structure
- `aa_pdb` (str): Path to reference all-atom structure
- `TAG` (str): Output file identifier

##### clean_pdb()
Cleans and standardizes PDB files for processing.

```python
clean_file = backMapper.clean_pdb(
    pdb="input.pdb",
    out_dir="cleaned/",
    name="protein_clean"
)
```

##### create_cg_model()
Creates CASM coarse-grained model from PDB file.

```python
cg_model = backMapper.create_cg_model(
    pdb="structure.pdb",
    ID="protein_cg"
)
```

##### OpenMM_vacuum_minimization()
Performs all-atom energy minimization using OpenMM.

```python
minimized = backMapper.OpenMM_vacuum_minimization(
    input_pdb="structure.pdb",
    maxcyc=1000
)
```

## Usage Examples

### Basic Coarse-Graining Workflow

```python
from NCLEdetector.change_resolution import CoarseGrain

# Initialize coarse-graining with 2x scaling
CoarseGrainer = CoarseGrain(
    pdbfile="1zmr_clean.pdb",
    ID="1ZMR",
    nscal=2.0,
    outdir="results/coarse_graining/",
    domain_file="domain_def.dat"
)

# Run coarse-graining process
print("Coarse-graining structure...")
CG_files = CoarseGrainer.run()

print(f"Generated files:")
for file_type, path in CG_files.items():
    print(f"  {file_type}: {path}")

# Generate OpenMM force field
CoarseGrainer.parse_cg_prm(
    prmfile=CG_files['prm'],
    topfile=CG_files['top']
)

print("Coarse-graining complete!")
```

### Complete Round-Trip Workflow

```python
from NCLEdetector.change_resolution import CoarseGrain, BackMapping

# Step 1: Coarse-grain all-atom structure
CoarseGrainer = CoarseGrain(
    pdbfile="original_structure.pdb",
    ID="protein_study",
    nscal=2.0,
    outdir="cg_output/",
    domain_file="domains.dat"
)

CG_files = CoarseGrainer.run()
CoarseGrainer.parse_cg_prm(CG_files['prm'], CG_files['top'])

# Step 2: Back-map to all-atom resolution
backMapper = BackMapping(
    nproc=4,
    outdir="backmapping_output/"
)

backMapper.backmap(
    cg_pdb=CG_files['cor'],
    aa_pdb="original_structure.pdb",  # Use original as reference
    TAG="reconstructed"
)

print("Round-trip conversion complete!")
```

### High-Throughput Processing

```python
import os
from NCLEdetector.change_resolution import CoarseGrain

# Process multiple structures
protein_list = [
    {"pdb": "1zmr.pdb", "id": "1ZMR", "domains": "1zmr_domains.dat"},
    {"pdb": "2xyz.pdb", "id": "2XYZ", "domains": "2xyz_domains.dat"},
    {"pdb": "3abc.pdb", "id": "3ABC", "domains": "3abc_domains.dat"}
]

for protein in protein_list:
    print(f"Processing {protein['id']}...")
    
    CoarseGrainer = CoarseGrain(
        pdbfile=protein["pdb"],
        ID=protein["id"],
        nscal=2.0,
        outdir=f"cg_results/{protein['id']}/",
        domain_file=protein["domains"]
    )
    
    CG_files = CoarseGrainer.run()
    CoarseGrainer.parse_cg_prm(CG_files['prm'], CG_files['top'])
    
    print(f"  Completed {protein['id']}")
```

## Integration with Scripts

### run_change_resolution.py
Complete resolution conversion pipeline:

```bash
python scripts/run_change_resolution.py \
  --outdir results/resolution_conversion/ \
  --pdbfile structure.pdb \
  --nscal 2 \
  --domain_file domain_def.dat \
  --ID protein_study
```

**Workflow Steps:**
1. **Coarse-graining**: Convert all-atom to coarse-grained representation
2. **Force field generation**: Create OpenMM-compatible parameters  
3. **Back-mapping**: Reconstruct all-atom structure
4. **Energy minimization**: Optimize reconstructed structure

The script demonstrates the complete workflow:

```python
# Coarse-grain the all-atom structure
CoarseGrainer = CoarseGrain(
    outdir=outdir,
    ID=ID,
    pdbfile=pdbfile,
    nscal=nscal,
    domain_file=domain_file
)

CGfiles = CoarseGrainer.run()

# Parse the prm and top file to make OpenMM compatible .xml file
CoarseGrainer.parse_cg_prm(prmfile=CGfiles['prm'], topfile=CGfiles['top'])

# Backmap the coarse-grained structure to all-atom  
backMapper = BackMapping(outdir=outdir)
backMapper.backmap(cg_pdb=CGfiles['cor'], aa_pdb=pdbfile, TAG=ID)
```

## Output Files

### Coarse-Graining Output
- `*_cg.cor`: Coarse-grained coordinates in CHARMM format
- `*_cg.psf`: Protein structure file with CG topology
- `*_cg.prm`: Parameter file with CG force field parameters
- `*_cg.top`: Topology file with atom type definitions
- `*_cg.xml`: OpenMM-compatible force field file
- `coarse_graining_log.txt`: Detailed process log

### Back-Mapping Output
- `*_backmapped.pdb`: Reconstructed all-atom structure
- `*_minimized.pdb`: Energy-minimized all-atom structure
- `*_clean.pdb`: Cleaned input structure
- `backmapping_log.txt`: Back-mapping process log
- `energy_minimization.log`: Minimization statistics

### Analysis Files
- `resolution_comparison.csv`: Structural comparison metrics
- `rmsd_analysis.txt`: RMSD between original and reconstructed
- `quality_metrics.csv`: Structure quality assessment

## Coarse-Graining Theory

### Scale Factor Selection
The coarse-graining scale factor (`nscal`) determines the level of detail:

- **nscal = 1.5**: Minimal coarse-graining, retains fine details
- **nscal = 2.0**: Standard coarse-graining for most applications  
- **nscal = 3.0**: Aggressive coarse-graining for large-scale simulations

### Interaction Potentials
Available potential types:
- **'bt'**: Betancourt-Thirumalai statistical potential
- **'kb'**: Knowledge-based statistical potential
- **'custom'**: User-defined interaction parameters

### Domain-Aware Coarse-Graining
Domain definitions enable:
- **Differential scaling**: Different domains with different scale factors
- **Interface preservation**: Maintain domain-domain interactions
- **Functional site protection**: Preserve critical functional regions

## Back-Mapping Methodology

### Template-Based Reconstruction
1. **Fragment library**: Use reference structure as template
2. **Coordinate alignment**: Align CG beads to reference positions
3. **Sidechain reconstruction**: Rebuild sidechain conformations
4. **Loop modeling**: Reconstruct flexible loop regions

### Energy Minimization Pipeline
1. **Steepest descent**: Initial clash removal
2. **Conjugate gradient**: Refined local optimization  
3. **Limited-memory BFGS**: Final structure optimization
4. **Restraint annealing**: Gradual restraint removal

## Quality Assessment

### Structural Metrics
- **RMSD**: Root mean square deviation from reference
- **GDT-TS**: Global distance test score
- **Local geometry**: Bond lengths, angles, dihedrals
- **Ramachandran analysis**: Backbone conformation validation

### Physical Validation
- **Energy analysis**: Potential energy components
- **Contact preservation**: Native contact recovery
- **Radius of gyration**: Overall compactness
- **Secondary structure**: DSSP assignment comparison

### Entanglement Preservation
- **GLN comparison**: Gaussian linking number conservation
- **Topological invariants**: Knot type preservation
- **Loop crossing patterns**: Entanglement configuration maintenance

## Performance Considerations

- **Memory Usage**: Scales linearly with protein size
- **Computation Time**: 
  - Coarse-graining: ~1-5 minutes per structure
  - Back-mapping: ~5-30 minutes depending on size
- **Parallelization**: Back-mapping benefits from multiple processors
- **Storage**: Intermediate files can be substantial

## Applications

### Multi-Scale Simulations
1. **Equilibration**: Use CG for rapid equilibration
2. **Conformational sampling**: Explore large-scale motions  
3. **Back-mapping**: Refine to atomic detail
4. **Production**: Run all-atom simulations

### Folding Studies
1. **CG folding**: Study folding pathways at coarse resolution
2. **State identification**: Find metastable states
3. **All-atom refinement**: Add atomic detail to key states
4. **Entanglement analysis**: Study topological changes

### Drug Design Applications
1. **Binding site mapping**: CG screening of binding sites
2. **Conformational states**: Sample receptor conformations
3. **All-atom docking**: Detailed drug-protein interactions
4. **Free energy calculations**: Precise binding affinity prediction

## Dependencies

- **Core**: NumPy, SciPy, MDAnalysis
- **Structure processing**: ProDy, PDBFixer
- **Energy minimization**: OpenMM
- **Force fields**: CHARMM parameter files
- **Parallel processing**: Multiprocessing
- **File I/O**: BioPython (optional)

