# Order Parameters Module

The `order_params` module provides comprehensive functionality for calculating structural order parameters from protein simulation trajectories, enabling quantitative analysis of folding progress and conformational changes.

## Overview

This module calculates essential order parameters for protein folding analysis:
- **Q**: Fraction of native contacts
- **G**: Fraction of native contacts with entanglement changes  
- **K**: Mirror symmetry parameter for coarse-grained artifacts
- **SASA**: Solvent accessible surface area
- **Jwalk**: Surface distance between residues
- **XP**: Cross-linking propensity

## Classes

### CalculateOP

Main class for order parameter calculations across simulation trajectories.

#### Initialization
```python
from EntDetect.order_params import CalculateOP

CalcOP = CalculateOP(
    outdir="results/OP/",
    Traj="1", 
    ID="protein_name",
    psf="structure.psf",
    cor="structure.cor", 
    dcd="trajectory.dcd",
    sec_elements="secondary_structure.txt",
    domain="domain_def.dat",
    start=0,
    end=-1,
    stride=1
)
```

**Parameters:**
- `outdir` (str): Output directory for results
- `Traj` (str): Trajectory identifier  
- `ID` (str): Analysis identifier
- `psf` (str): Path to PSF topology file
- `cor` (str): Path to reference coordinates (native structure)
- `dcd` (str): Path to trajectory file
- `sec_elements` (str): Secondary structure definition file
- `domain` (str): Domain definition file  
- `start` (int): Starting frame (default: 0)
- `end` (int): Ending frame (default: -1 for all frames)
- `stride` (int): Frame stride (default: 1)

#### Key Methods

##### Q()
Calculates fraction of native contacts throughout the trajectory.

```python
Q_data = CalcOP.Q()
```

**Returns:**
- Dictionary containing Q values for each frame and metadata

##### G()
Calculates fraction of native contacts with entanglement status changes.

```python
G_data = CalcOP.G(topoly=True, Calpha=True, CG=True, nproc=10)
```

**Parameters:**
- `topoly` (bool): Use topological linking number (True) vs GLN (False)
- `Calpha` (bool): Use Cα atoms for contacts vs heavy atoms
- `CG` (bool): Coarse-grained trajectory mode
- `nproc` (int): Number of parallel processes

**Returns:**
- Dictionary with G values and entanglement change data

##### K()
Calculates mirror symmetry parameter for detecting CG artifacts.

```python
K_data = CalcOP.K()
```

**Returns:**
- Dictionary with K values indicating chirality preservation

## Usage Examples

### Basic Order Parameter Analysis

```python
from EntDetect.order_params import CalculateOP

# Initialize for coarse-grained trajectory
CalcOP = CalculateOP(
    outdir="results/OP_analysis/",
    Traj="1",
    ID="1ZMR", 
    psf="1zmr_ca.psf",
    cor="1zmr_ca.cor",
    dcd="trajectory.dcd",
    sec_elements="secondary_struct.txt",
    domain="domain_def.dat",
    start=0
)

# Calculate native contacts
Q_data = CalcOP.Q()
print(f"Q data keys: {Q_data.keys()}")

# Calculate entanglement changes  
G_data = CalcOP.G(topoly=True, Calpha=True, CG=True, nproc=10)
print(f"G data keys: {G_data.keys()}")

# Check for mirror artifacts
K_data = CalcOP.K()
print(f"K data keys: {K_data.keys()}")
```

### All-Atom Trajectory Analysis

```python
# Initialize for all-atom trajectory
CalcOP_AA = CalculateOP(
    outdir="results/AA_analysis/",
    Traj="1",
    ID="1ZMR_AA",
    psf="1zmr_allatom.psf", 
    cor="1zmr_allatom.pdb",
    dcd="aa_trajectory.dcd",
    sec_elements="secondary_struct.txt",
    domain="domain_def.dat"
)

# Use GLN instead of topological linking
G_data_AA = CalcOP_AA.G(topoly=False, Calpha=False, CG=False, nproc=8)
```

## Integration with Scripts

This module is central to the standardized workflow scripts:

### run_OP_on_simulation_traj.py
Complete order parameter analysis pipeline:

```bash
python scripts/run_OP_on_simulation_traj.py \
  --Traj 1 \
  --PSF structure.psf \
  --DCD trajectory.dcd \
  --ID protein_name \
  --COR structure.cor \
  --sec_elements secondary_struct.txt \
  --domain domain_def.dat \
  --outdir results/OP/ \
  --start 0
```

The script demonstrates different analysis modes:
- **Case 1**: Topological linking + Cα contacts + CG trajectory
- **Case 2**: GLN + Cα contacts + CG trajectory  
- **Case 3**: GLN + Heavy atoms + All-atom trajectory
- **Case 4**: Topological + Heavy atoms + All-atom trajectory

## Output Files

### Q Parameter Files
- `Q_data.pkl`: Pickled Q data dictionary
- `Q_timeseries.csv`: Q values vs simulation time
- `native_contacts.csv`: Native contact definitions

### G Parameter Files  
- `G_data.pkl`: Pickled G data with entanglement information
- `G_timeseries.csv`: G values throughout trajectory
- `entanglement_changes.csv`: Detailed entanglement transition data
- `cluster_data/`: Individual cluster entanglement evolution

### K Parameter Files
- `K_data.pkl`: Mirror symmetry data
- `K_timeseries.csv`: Chirality preservation over time
- `secondary_structure_chirality.csv`: Per-element chirality tracking

## Theory Background

### Q Parameter (Fraction of Native Contacts)
Measures folding progress by calculating the fraction of native contacts present:

```
Q(t) = N_formed(t) / N_native
```

Where N_formed(t) are contacts formed at time t and N_native are total native contacts.

### G Parameter (Entanglement-Sensitive Contacts)
Extension of Q that only counts contacts whose formation/breaking involves entanglement changes:

```
G(t) = N_entangled_contacts(t) / N_native
```

### K Parameter (Mirror Symmetry)
Detects coarse-grained model artifacts where secondary structures pack with opposite chirality:

```
K(t) = |χ_left(t) - χ_right(t)| / (χ_left(t) + χ_right(t))
```

## Performance Considerations

- **Memory**: Scales linearly with trajectory length
- **CPU**: G parameter calculation benefits significantly from parallelization
- **Storage**: Output files can be large for long trajectories
- **I/O**: Efficient trajectory reading with MDAnalysis

## Dependencies

- MDAnalysis: Trajectory analysis
- NumPy: Numerical computations  
- Pandas: Data handling
- Multiprocessing: Parallel G parameter calculation
- Pickle: Data serialization
- `G(topoly=True)`: Calculates the fraction of native contacts with a change in entanglement status (G).
- `K()`: Calculates the mirror symmetry parameter (K) for the trajectory.
- `SASA()`: Calculates the solvent accessible surface area (SASA) for each frame (all-atom trajectories).
- `XP(pdb)`: Calculates the cross-linking propensity (XP) for all pairs of cross-linkable residues in a PDB file.

## Usage

The module is typically used via the provided scripts:
- `scripts/run_OP_CGtrajAnal.py` for coarse-grained trajectories.
- `scripts/run_OP_AAtrajAnal.py` for all-atom trajectories.

These scripts call the relevant methods in `CalculateOP` to compute and output the desired order parameters for further analysis.

## Example Workflow

```python
from EntDetect.order_params import CalculateOP

CalcOP = CalculateOP(
    outdir='results',
    Traj=1,
    ID='protein1',
    psf='protein.psf',
    cor='protein.cor',
    dcd='traj.dcd',
    sec_elements='sec_elements.txt',
    domain='domain.txt',
    start=0,
    end=1000,
    stride=1
)

Qdata = CalcOP.Q()
Gdata = CalcOP.G(topoly=True)
Kdata = CalcOP.K()
SASAdata = CalcOP.SASA()
XPdata = CalcOP.XP(pdb_file='protein.pdb')
```

**Notes:**
- The class is designed to be used in automated pipelines for large-scale simulation analysis.
- Output files and formats are determined by the calling scripts.
- For further details, refer to the source code in `EntDetect/order_params.py` or the package documentation.

