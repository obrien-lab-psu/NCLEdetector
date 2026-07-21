# Gaussian Entanglement Module

The `gaussian_entanglement` module provides the core functionality for calculating and analyzing non-covalent lasso entanglements in protein structures using Gaussian entanglement theory.

## Overview

This module implements methods to:
- Calculate Gaussian linking numbers between protein loops
- Identify native entanglements in protein structures  
- Filter high-quality entanglements based on structural criteria
- Support both all-atom and coarse-grained representations

## Classes

### GaussianEntanglement

Main class for entanglement calculations and analysis.

#### Initialization
```python
from NCLEdetector.gaussian_entanglement import GaussianEntanglement

ge = GaussianEntanglement(g_threshold=0.6, density=0.0, Calpha=False, CG=False)
```

**Parameters:**
- `g_threshold` (float): Threshold for Gaussian entanglement score (default: 0.6)
- `density` (float): Density parameter for contact filtering (default: 0.0)  
- `Calpha` (bool): Use Cα atoms for native contacts vs 4.5Å heavy atom cutoff (default: False)
- `CG` (bool): Enable coarse-grained mode (default: False)

#### Key Methods

##### calculate_native_entanglements()
Identifies all native entanglements in a protein structure.

```python
result = ge.calculate_native_entanglements(
    pdb_file="/path/to/structure.pdb",
    outdir="results/native_entanglements/", 
    ID="protein_name",
    chain="A"
)
```

**Parameters:**
- `pdb_file` (str): Path to cleaned PDB structure file
- `outdir` (str): Output directory for results
- `ID` (str): Identifier for the analysis (used in output filename)
- `chain` (str, optional): Chain identifier (default: processes all chains)

**Returns:**
- Dictionary with output file paths and entanglement data

**Output File:**
- Saved as: `{outdir}/{ID}_GE.csv`
- Format: Pipe-separated CSV with columns: `ID|chain|i|j|crossingsN|crossingsC|gn|gc|GLNn|GLNc|TLNn|TLNc|CCbond|ENT`

##### select_high_quality_entanglements()
Filters entanglements based on structural quality criteria (removes slipknots).

```python
hq_result = ge.select_high_quality_entanglements(
    entanglement_file="native_entanglements.csv",
    pdb_file="/path/to/structure.pdb",
    outdir="results/hq_entanglements/",
    ID="protein_name", 
    model="EXP",
    chain="A"
)
```

**Parameters:**
- `entanglement_file` (str): Path to entanglement data file (from calculate_native_entanglements)
- `pdb_file` (str): Path to structure file
- `outdir` (str): Output directory
- `ID` (str): Analysis identifier  
- `model` (str): Model type ("EXP" for experimental, "AF" for AlphaFold)
- `chain` (str, optional): Chain identifier

**Returns:**
- Dictionary with output file paths and filtered entanglement data

**Output File:**
- Saved as: `{outdir}/{ID}.csv`
- Format: Same as input with additional `Slipknot_N` and `Slipknot_C` columns

## Usage Examples

### Basic Native Entanglement Analysis

```python
from NCLEdetector.gaussian_entanglement import GaussianEntanglement

# Initialize for experimental structure analysis
ge = GaussianEntanglement(g_threshold=0.6, Calpha=False, CG=False)

# Calculate native entanglements for chain A
native_ents = ge.calculate_native_entanglements(
    pdb_file="structure.pdb",
    outdir="results/native/", 
    ID="1ZMR_A",
    chain="A"
)

# Filter for high quality entanglements
hq_ents = ge.select_high_quality_entanglements(
    entanglement_file=native_ents['outfile'],
    pdb_file="structure.pdb",
    outdir="results/hq/",
    ID="1ZMR_A",
    model="EXP",
    chain="A"
)
```

### Multi-Chain Analysis

```python
from NCLEdetector.gaussian_entanglement import GaussianEntanglement
import MDAnalysis as mda

# Initialize
ge = GaussianEntanglement(g_threshold=0.6, Calpha=False, CG=False)

# Get all chains from structure
u = mda.Universe("multi_chain.pdb")
chains = sorted(set([atom.chainID for atom in u.atoms]))

# Process each chain
for chain_id in chains:
    # Use chain suffix for file naming
    chain_id_suffix = f"1PKL_{chain_id}"
    
    native_ents = ge.calculate_native_entanglements(
        pdb_file="multi_chain.pdb",
        outdir="results/Native_GE/",
        ID=chain_id_suffix,
        chain=chain_id
    )
    # Results saved as: results/Native_GE/1PKL_A_GE.csv, 1PKL_B_GE.csv, etc.
```

### Coarse-Grained Analysis

```python
# Initialize for coarse-grained trajectory analysis
ge_cg = GaussianEntanglement(g_threshold=0.6, Calpha=True, CG=True)

# Process coarse-grained structure
cg_result = ge_cg.calculate_native_entanglements(
    pdb_file="cg_structure.pdb",
    outdir="results/cg_native/",
    ID="1ZMR_CG",
    chain="A"
)
```

### Converting CHARMM Files to PDB

For coarse-grained models in CHARMM format (.cor/.psf), use the conversion utility:

```bash
# Convert CHARMM COR/PSF to PDB
python scripts/convert_cor_psf_to_pdb.py \
  --cor structure.crd \
  --psf structure.psf \
  --output structure.pdb

# Then analyze with --cg flag
python scripts/run_nativeNCLE.py \
  --struct structure.pdb \
  --outdir results/ \
  --cg \
  --organism Ecoli
```

**Note:** The conversion utility automatically sets chain ID to 'A' for compatibility with NCLEdetector.

## Integration with Scripts

This module is used extensively in the standardized workflow scripts:

### run_nativeNCLE.py
Complete pipeline for native entanglement analysis:

```bash
# All-atom analysis (default)
python scripts/run_nativeNCLE.py \
  --struct structure.pdb \
  --outdir results/native_analysis/ \
  --ID protein_name \
  --organism Ecoli \
  --Accession P00558

# Coarse-grained (C-alpha only) analysis
python scripts/run_nativeNCLE.py \
  --struct cg_structure.pdb \
  --outdir results/cg_analysis/ \
  --cg \
  --organism Ecoli \
  --Accession P00558

# Single chain analysis
python scripts/run_nativeNCLE.py \
  --struct structure.pdb \
  --outdir results/chain_A/ \
  --chain A \
  --organism Ecoli

# Multi-chain analysis (processes all chains automatically)
python scripts/run_nativeNCLE.py \
  --struct multi_chain.pdb \
  --outdir results/multi_chain/ \
  --organism Ecoli
```

**Parameters:**
- `--struct`: Path to PDB structure file (required)
- `--outdir`: Output directory for results (required)
- `--ID`: Analysis identifier (optional, defaults to structure basename)
- `--chain`: Specific chain to analyze (optional, analyzes all chains if omitted)
- `--organism`: Organism name for clustering: Ecoli, Human, or Yeast (default: Ecoli)
- `--Accession`: UniProt Accession for the protein (default: P00558)
- `--cg`: Flag to indicate coarse-grained (C-alpha only) model

### run_OP_on_simulation_traj.py  
Uses entanglement calculations for G parameter computation:

```bash
python scripts/run_OP_on_simulation_traj.py \
  --Traj 1 \
  --PSF structure.psf \
  --DCD trajectory.dcd \
  --ID protein_name \
  --COR structure.cor \
  --sec_elements secondary.txt \
  --domain domain_def.dat \
  --outdir results/OP/ \
  --start 0
```

## Output Files

### Native Entanglement Files
Files are saved in pipe-separated CSV format (`.csv` extension with `|` delimiter):

**Directory Structure:**
- `Native_GE/`: Raw entanglement data
  - `{ID}_GE.csv`: Gaussian entanglement calculations
  - For multi-chain: `{ID}_A_GE.csv`, `{ID}_B_GE.csv`, etc.
- `Native_HQ_GE/`: High-quality filtered entanglements
  - `{ID}.csv`: Slipknot-removed entanglements
- `Native_clustered_HQ_GE/`: Clustered entanglements
  - `{ID}.csv`: Degeneracy-removed entanglements
- `Native_clustered_HQ_GE_features/`: Feature calculations
  - `{Accession}_{ID}_{chain}_uent_features.csv`: Unique entanglement features

### File Format
Pipe-separated CSV files (`|` delimiter) with columns:

**Gaussian Entanglement Output (`*_GE.csv`):**
- `ID`: Structure identifier
- `chain`: Chain identifier
- `i`: Loop start residue
- `j`: Loop end residue
- `crossingsN`: N-terminal crossing residues (e.g., "+108")
- `crossingsC`: C-terminal crossing residues (e.g., "+18")
- `gn`: N-terminal Gaussian linking number
- `gc`: C-terminal Gaussian linking number
- `GLNn`: N-terminal topological linking number (integer)
- `GLNc`: C-terminal topological linking number (integer)
- `TLNn`: N-terminal total linking number
- `TLNc`: C-terminal total linking number
- `CCbond`: Disulfide bond flag (True/False)
- `ENT`: Entanglement flag (True/False)

**High-Quality Output (`Native_HQ_GE/*.csv`):**
Same columns as GE output plus:
- `Slipknot_N`: N-terminal slipknot flag
- `Slipknot_C`: C-terminal slipknot flag

**Clustered Output (`Native_clustered_HQ_GE/*.csv`):**
Same columns as HQ output with degeneracies removed

**Feature Output (`Native_clustered_HQ_GE_features/*_uent_features.csv`):**
Extended feature set including:
- `gene`: UniProt Accession
- `PDB`: Structure identifier
- `chain`: Chain identifier
- `ENT-ID`: Entanglement cluster ID
- `loopsize`: Size of the loop
- `ent_region`: Entangled region
- `prot_size`: Total protein size
- `ACO`, `RCO`: Contact order metrics
- Additional geometric and topological features

### Multi-Chain Behavior
When processing multi-chain structures:
- All chains write to the same `Native_GE/` directory
- Files are named with chain suffix: `1PKL_A_GE.csv`, `1PKL_B_GE.csv`, etc.
- Each chain is processed independently through the full pipeline
- Feature files are generated per chain: `{Accession}_{ID}_A_uent_features.csv`

## Theory Background

The Gaussian entanglement method calculates linking numbers between protein loops using:

1. **Loop Identification**: Protein backbone segments between secondary structures
2. **GLN Calculation**: Gaussian linking number between loop pairs  
3. **Threshold Application**: Entanglements identified above g_threshold
4. **Quality Filtering**: Structural validation based on contact patterns

## Performance Notes

- **Memory Usage**: Scales quadratically with protein size
- **Computation Time**: ~1-10 minutes for typical protein structures
- **Parallelization**: Thread-safe for multiple structure analysis
- **File I/O**: Optimized CSV output for large datasets

## Dependencies

- NumPy: Numerical computations
- Pandas: Data manipulation  
- MDAnalysis: Structure file parsing
- SciPy: Statistical functions
select_high_quality_entanglements(entanglement_file, pdb_file, outdir, ID, model, mapping)
```
Filters and selects high-quality entanglements based on additional criteria.
- `entanglement_file` (str): Path to the file containing raw entanglement data.
- `pdb_file` (str): Path to the cleaned PDB file.
- `outdir` (str): Output directory for filtered results.
- `ID` (str): Identifier for the current analysis.
- `model` (str): Model type (e.g., 'EXP' or 'AF').
- `mapping` (str): Mapping file or 'None'.

**Returns:**
A dictionary containing output file paths and filtered entanglement data.

---

## Usage Example
```python
from NCLEdetector.gaussian_entanglement import GaussianEntanglement

ge = GaussianEntanglement(g_threshold=0.6, density=0.0, Calpha=False, CG=False)
native_ent = ge.calculate_native_entanglements("protein_clean.pdb", outdir="results/Native_GE", ID="protein1_A", chain="A")
hq_ent = ge.select_high_quality_entanglements(native_ent['outfile'], "protein_clean.pdb", outdir="results/Native_HQ_GE", ID="protein1_A", model="AF", mapping="None", chain="A")
```

**Notes:**
- The class is designed to be used in protein structure analysis pipelines, especially for large-scale or automated workflows.
- Output directories will be created if they do not exist.
- For best results, input PDB files should be pre-processed and cleaned.
- **Chain Tracking**: All output files include a `chain` column to track which chain each entanglement belongs to.
- **Crossing Format**: Crossings are split into `crossingsN` (N-terminal) and `crossingsC` (C-terminal) columns.
- **File Format**: All output files use pipe-separated CSV format (`.csv` extension with `|` delimiter).
- **Multi-Chain Support**: When processing structures with multiple chains, each chain is analyzed separately and saved with a chain suffix (e.g., `1PKL_A_GE.csv`, `1PKL_B_GE.csv`).
- **Coarse-Grained Models**: Use `--cg` flag with run_nativeNCLE.py and set `Calpha=True, CG=True` for GaussianEntanglement initialization.
- **CHARMM Format**: Convert .cor/.psf files to .pdb using the `convert_cor_psf_to_pdb.py` utility before analysis.

For further details, refer to the source code in `NCLEdetector/gaussian_entanglement.py` or the package documentation.
