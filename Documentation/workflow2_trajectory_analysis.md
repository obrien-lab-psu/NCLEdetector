# Workflow 2: Detect Changes in Entanglement Across Simulation Trajectories

← [Back to Master Index](index.md)

---

## Goal

Compute order parameters (Q, G, K, SASA, XP) across simulation trajectories, cluster non-native entanglement changes, build a Markov state model (MSM), and analyze metastable state behavior.

---

## Typical runtime

| Step | Runtime |
|------|---------|
| Q, G, K for one CG trajectory (nproc=10) | 12–20 hours |
| SASA and XP for one AA trajectory (nproc=10) | 2–6 hours |
| Q, G, K for all 1000 trajectories (cluster) | Hours to days |
| Non-native clustering | Hours (memory-intensive) |
| MSM construction | Minutes |
| MSM statistics / folding pathways | Minutes |

---

## Pre-computed outputs

All 1000 CG and all-atom trajectories have already been analyzed and their outputs are stored in the DATASTORE:

```
$DATASTORE/outputs/workflow2/
├── OP/                        # CG order parameters (all 1000 trajectories)
│   ├── Q/1ZMR_Traj{N}.Q
│   ├── G/1ZMR_Traj{N}.G
│   │   └── Combined_GE/            # Per-trajectory entanglement pkl files
│   │       └── 1ZMR_traj{N}_GE.pkl
│   └── K/K_{N}_prod.dat
└── OP_AA/                     # All-atom order parameters (all 1000 trajectories)
    ├── SASA/1ZMR.SASA
    └── XP/
        ├── Jwalk_results/
        ├── XLresidue_pairs_Full.csv
        └── XLresidue_pairs.txt
```

Where:
```bash
DATASTORE=/scratch/ims86/EntDetect_Datastore
OUTDIR=$DATASTORE/outputs/workflow2
```

---

## Required input files

| File | Path | Notes |
|------|------|-------|
| Cα PSF topology | `$REFSTRUCT/1zmr_model_clean_ca.psf` | CG Q/G/K |
| Cα COR reference | `$REFSTRUCT/1zmr_model_clean_ca.cor` | CG Q/G/K |
| CG trajectory (demo) | `$CG_TRAJ_DIR/420_prod.dcd` | Single trajectory for testing |
| CG trajectories (all) | `$CG_TRAJ_DIR/{N}_prod.dcd` (N=1–1000) | Full production run |
| All-atom PDB topology | `$REFSTRUCT/1zmr_model_clean.pdb` | SASA/XP |
| AA trajectory (demo) | `$AA_TRAJ_DIR/420_prod_aa.dcd` | AA demo |
| AA trajectories (all) | `$AA_TRAJ_DIR/{N}_prod_aa.dcd` (N=1–1000) | Full AA production run |
| Secondary structure defs | `$REFSTRUCT/secondary_struc_defs.txt` | **Required for Q/G/K** |
| Domain boundary file | `$REFSTRUCT/domain_def.dat` | **Required for Q/G/K** |
| Trajectory-to-pkl mapping | `$DATASTORE/user_input/metadata/trajnum2file.txt` | **Required for clustering** |

```bash
REFSTRUCT=$DATASTORE/user_input/reference_structures
CG_TRAJ_DIR=$DATASTORE/user_input/cg_trajectories
AA_TRAJ_DIR=$DATASTORE/user_input/aa_trajectories
```

### `trajnum2file.txt` — maps trajectory numbers to pkl files

This file maps trajectory numbers to their **G-order-parameter pkl files** (not DCD files). The format is comma-separated:

```
trajnum,pklfile
<trajectory_number>,<path_to_pkl_file>
```

The pre-populated copy at `$DATASTORE/user_input/metadata/trajnum2file.txt` maps all 1000 trajectories to the Combined_GE pkl files in the DATASTORE. It can be regenerated if needed:

```bash
echo "trajnum,pklfile" > $DATASTORE/user_input/metadata/trajnum2file.txt
for pkl in $DATASTORE/outputs/workflow2/OP/G/Combined_GE/*.pkl; do
    num=$(basename $pkl | sed 's/1ZMR_traj\([0-9]*\)_GE.pkl/\1/')
    echo "$num,$pkl"
done >> $DATASTORE/user_input/metadata/trajnum2file.txt
```

---

## Step 1. Activate your environment and set paths

```bash
source ~/.bashrc
conda activate entdetect

DATASTORE=/scratch/ims86/EntDetect_Datastore
CG_TRAJ_DIR=$DATASTORE/user_input/cg_trajectories
AA_TRAJ_DIR=$DATASTORE/user_input/aa_trajectories
REFSTRUCT=$DATASTORE/user_input/reference_structures
OUTDIR=$DATASTORE/outputs/workflow2

mkdir -p $OUTDIR/OP/G $OUTDIR/OP/Q $OUTDIR/OP/K \
         $OUTDIR/OP_AA/SASA $OUTDIR/OP_AA/XP \
         $OUTDIR/nonnative_clustering $OUTDIR/MSM \
         $OUTDIR/MSM_StateProbabilityStats
```

---

## Step 2. Compute order parameters Q, G, K on the CG trajectory

Compute the three canonical order parameters for the coarse-grained trajectory using the parameters from `assets/slurm/scripts/run_OP_traj420.slurm`.

### 2a. Initialize `CalculateOP` for the CG trajectory

```python
from EntDetect.order_params import CalculateOP

DATASTORE    = "/scratch/ims86/EntDetect_Datastore"
REFSTRUCT    = f"{DATASTORE}/user_input/reference_structures"
CG_TRAJ_DIR  = f"{DATASTORE}/user_input/cg_trajectories"
OUTDIR       = f"{DATASTORE}/outputs/workflow2"

Traj         = 420
PSF          = f"{REFSTRUCT}/1zmr_model_clean_ca.psf"
COR          = f"{REFSTRUCT}/1zmr_model_clean_ca.cor"
DCD          = f"{CG_TRAJ_DIR}/420_prod.dcd"
ID           = "1ZMR"
sec_elements = f"{REFSTRUCT}/secondary_struc_defs.txt"
domain       = f"{REFSTRUCT}/domain_def.dat"
outdir       = f"{OUTDIR}/OP"

# start: first frame to include (0-indexed).
# Adjust to skip early equilibration frames for production runs.
start = 0

CalcOP = CalculateOP(
    outdir=outdir,
    Traj=Traj,
    ID=ID,
    psf=PSF,
    cor=COR,
    sec_elements=sec_elements,
    dcd=DCD,
    domain=domain,
    start=start,
    ent_detection_method=1,   # 1 = GLN-only; matches production run settings
)
```

### 2b. Compute Q — fraction of native contacts

```python
Qdata_dict = CalcOP.Q()
```

`Q` measures how many of the native residue–residue contacts present in the reference structure are also present in each trajectory frame. A value near 1.0 indicates a native-like conformation.

**Output:** A `.Q` file in `$OUTDIR/OP/Q/`

### 2c. Compute G — entanglement order parameter

```python
Gdata_dict = CalcOP.G(topoly=False, Calpha=True, CG=True, nproc=10, chunk_frames=100, chunk_suffix='_chunk')
```

`G` captures the fraction of native contacts that exhibit a **change in entanglement state** relative to the native structure.

| Argument | Value | Meaning |
|----------|-------|----------|
| `topoly` | `False` | Use GLN-only workflow (no Topoly linking numbers) |
| `Calpha` | `True` | Use Cα-defined contacts (appropriate for CG trajectories) |
| `CG` | `True` | Input trajectory is coarse-grained |
| `nproc` | `10` | Number of CPU cores to use |
| `chunk_frames` | `100` | Write intermediate results every N frames (reduces memory usage on large trajectories) |
| `chunk_suffix` | `'_chunk'` | Suffix for chunked output files |

**Output:** A `.G` file and per-frame entanglement metadata `.pkl` in `$OUTDIR/OP/G/`

> **Runtime note:** `G` is the most expensive order parameter to compute. Expect 12–20 hours per CG trajectory at ~6700 frames on 10 cores. Submit via the cluster for a full 1000-trajectory run (see [Running the full script](#running-the-full-workflow-as-a-single-script)).

### 2d. Compute K — mirror-symmetry order parameter

```python
Kdata_dict = CalcOP.K()
```

`K` detects frames where the protein has adopted a **mirror-image conformation** relative to the native structure. These frames are artifacts that must be removed before clustering.

**Output:** A `K_*.dat` file in `$OUTDIR/OP/K/`

---

## Step 3. Compute order parameters SASA and XP on the all-atom trajectory

> **Why separate from Q, G, K?** The CG Cα-only representation does not carry enough atomic detail for accurate SASA (Shrake-Rupley requires explicit side-chain/backbone atoms) or cross-link SASD (Jwalk uses solvent-accessible surface geometry). The AA trajectories are produced by back-mapping the CG frames. Q, G, K are computed only on CG trajectories (Step 2).

### 3a. Initialize `CalculateOP` for the all-atom trajectory

```python
from EntDetect.order_params import CalculateOP

DATASTORE    = "/scratch/ims86/EntDetect_Datastore"
REFSTRUCT    = f"{DATASTORE}/user_input/reference_structures"
AA_TRAJ_DIR  = f"{DATASTORE}/user_input/aa_trajectories"
OUTDIR       = f"{DATASTORE}/outputs/workflow2"

Traj     = 420
AA_PDB   = f"{REFSTRUCT}/1zmr_model_clean.pdb"    # all-atom topology / reference
AA_DCD   = f"{AA_TRAJ_DIR}/420_prod_aa.dcd"
ID       = "1ZMR"

CalcOP = CalculateOP(
    outdir=f"{OUTDIR}/OP_AA",
    Traj=Traj,
    ID=ID,
    psf=AA_PDB,
    cor=AA_PDB,
    dcd=AA_DCD,
    sec_elements=f"{REFSTRUCT}/secondary_struc_defs.txt",
    domain=f"{REFSTRUCT}/domain_def.dat",
    start=0,
)
```

### 3b. Compute SASA — solvent-accessible surface area per residue

```python
SASAdata_dict = CalcOP.SASA()
```

Uses the Shrake-Rupley algorithm (mdtraj, `probe_radius=0.14 nm`, 1000 sphere points) to compute per-residue SASA for every frame. Used downstream (Workflow 3) to test LiP-MS signals: high SASA residues are more accessible to the protease.

| Key | Contents |
|-----|----------|
| `outfile` | Path to `{ID}.SASA` CSV written to `$OUTDIR/OP_AA/SASA/` |
| `result`  | DataFrame with columns `Time(ns)`, `Frame`, `resid`, `SASA(nm^2)` |

### 3c. Compute XP — Jwalk cross-link probability

```python
XPdata_dict = CalcOP.XP(pdb=AA_PDB)
```

`XP` computes the **solvent-accessible surface distance (SASD)** between all pairs of cross-linkable residue types (K, S, T, Y, M) and converts each SASD to a cross-link probability score using a Gaussian parameterised on K–K linker geometry. Used downstream to test XL-MS signals.

**Output:** Per-residue-pair XP scores in `$OUTDIR/OP_AA/XP/Jwalk_results/{stem}_crosslink_list.txt`

---

## Running the order parameter analysis as a single script

Order parameters are computed by `scripts/run_OP_on_simulation_traj.py`. Q, G, K are computed for CG trajectories; SASA and XP are computed for AA trajectories. This is the script used for the full 1000-trajectory production run.

### CG run — Q, G, K

```bash
source ~/.bashrc
conda activate entdetect

DATASTORE=/scratch/ims86/EntDetect_Datastore
REFSTRUCT=$DATASTORE/user_input/reference_structures

python scripts/run_OP_on_simulation_traj.py \
    --Traj 420 \
    --PSF  $REFSTRUCT/1zmr_model_clean_ca.psf \
    --COR  $REFSTRUCT/1zmr_model_clean_ca.cor \
    --DCD  $DATASTORE/user_input/cg_trajectories/420_prod.dcd \
    --resolution cg \
    --contacts calpha \
    --ID   1ZMR \
    --sec_elements $REFSTRUCT/secondary_struc_defs.txt \
    --domain       $REFSTRUCT/domain_def.dat \
    --outdir       $DATASTORE/outputs/workflow2/OP \
    --logdir       $DATASTORE/outputs/workflow2/OP/logs \
    --start        0 \
    --ent_detection_method 1 \
    --nproc        10 \
    --ops Q G K \
    --no_topoly \
    --chunk_frames 100 \
    --chunk_suffix _chunk
```

### AA run — SASA, XP

```bash
python scripts/run_OP_on_simulation_traj.py \
    --Traj   420 \
    --ID     1ZMR \
    --PSF    $REFSTRUCT/1zmr_model_clean.pdb \
    --COR    $REFSTRUCT/1zmr_model_clean.pdb \
    --DCD    $DATASTORE/user_input/aa_trajectories/420_prod_aa.dcd \
    --resolution aa \
    --contacts calpha \
    --sec_elements $REFSTRUCT/secondary_struc_defs.txt \
    --domain       $REFSTRUCT/domain_def.dat \
    --outdir $DATASTORE/outputs/workflow2/OP_AA \
    --logdir $DATASTORE/outputs/workflow2/OP_AA/logs \
    --start  0 \
    --nproc        10 \
    --xp_pdb $REFSTRUCT/1zmr_model_clean.pdb \
    --ops SASA XP
```

### Submitting all 1000 trajectories to the cluster

For a full production run, wrap each command in a SLURM script (see `assets/slurm/scripts/run_OP_traj{N}.slurm` for the template). Use `assets/slurm/scripts/gen_slurms.py` to generate one script per trajectory and submit with:

```bash
for i in $(seq 1 1000); do sbatch assets/slurm/scripts/run_OP_traj${i}.slurm; done
```

| Flag | CG value | AA value | Notes |
|------|---------|---------|-------|
| `--resolution` | `cg` | `aa` | Determines contact-type default and CG flag |
| `--ent_detection_method` | `1` (GLN) | `2` (TLN) | Detection strategy |
| `--no_topoly` | present | absent | GLN-only for CG; TLN available for AA |
| `--ops` | `Q G K` | `SASA XP` | OPs computed |

---

## Step 4. Identify and remove artificial mirror conformations

Before clustering entanglement changes, remove trajectories or frames that are mirror-image artifacts.

### Recommended procedure

1. Load Q and K time series for all trajectories.
2. Identify trajectories where K is persistently high and Q is elevated or anomalous.
3. Visually inspect flagged frames in VMD.
4. Record the trajectory numbers that are confirmed mirrors in a list (used as `rm_traj_list` in later steps).

> **Important:** The cutoff values for flagging mirror conformations must be tuned for your system by examining the Q and K distributions. 

---

## Step 5. Cluster non-native entanglement changes

This step identifies non-redundant changes in entanglement topology across all trajectories. The pkl file paths are specified in the trajectory-to-pkl mapping CSV (`trajnum2pklfile_path`), which serves as the single source of truth for which pkl files to analyze.

> **Memory warning:** This step can require tens of gigabytes of RAM. Run on a high-memory node.

### 5a. Build the clustering with Python

```python
from EntDetect.clustering import ClusterNonNativeEntanglements

DATASTORE            = "/scratch/ims86/EntDetect_Datastore"
OUTDIR               = f"{DATASTORE}/outputs/workflow2"
CG_TRAJ_DIR          = f"{DATASTORE}/user_input/cg_trajectories"

trajnum2pklfile_path = f"{DATASTORE}/user_input/metadata/trajnum2file.txt"
traj_dir_prefix      = CG_TRAJ_DIR
outdir               = f"{OUTDIR}/nonnative_clustering"

clustering_NNents = ClusterNonNativeEntanglements(
    trajnum2pklfile_path=trajnum2pklfile_path,
    traj_dir_prefix=traj_dir_prefix,
    outdir=outdir,
)

# start_frame, end_frame: frame range to analyze (0-indexed)
# Default is to use all frames; here we use the last 67 frames for a faster demo
clustering_NNents.cluster(start_frame=6600, end_frame=6667)
```

### 5b. Using the command-line interface

For convenience, use the `scripts/run_nonnative_entanglement_clustering.py` script directly (see "Running non-native clustering as a single script" below for full details).

### 5c. Inspecting clustering results

```python
import pandas as pd
import numpy as np

DATASTORE  = "/scratch/ims86/EntDetect_Datastore"
clust_dir  = f"{DATASTORE}/outputs/workflow2/nonnative_clustering"

# Representative entanglement changes
rep_df = pd.read_csv(f"{clust_dir}/rep_chg_ent_topoly_linking_number.csv")
print(f"Number of representative entanglement changes: {len(rep_df)}")
print(rep_df.head())

# Per-frame structural assignments
chg_df = pd.read_csv(f"{clust_dir}/chg_ent_struct_topoly_linking_number.csv")
print(chg_df.head())

# Cluster data array
cluster_data = np.load(f"{clust_dir}/cluster_data_topoly_linking_number.npz", allow_pickle=True)
print("Available arrays:", list(cluster_data.keys()))
```

### 5d. Expected outputs

| File | Contents |
|------|----------|
| `rep_chg_ent_topoly_linking_number.csv` | Representative entanglement changes (one per cluster) |
| `chg_ent_struct_topoly_linking_number.csv` | Per-frame cluster assignment |
| `cluster_data_topoly_linking_number.npz` | Compressed cluster data array |
| `cluster_tree_topoly_linking_number.dat` | Text representation of the clustering hierarchy |
| `rep_chg_ent_list_topoly_linking_number.pkl` | List of representative entanglement objects |
| `chg_ent_topoly_linking_number_distribution.pdf` | Distribution plot of loop/crossing residues |

---

## Running non-native clustering as a single script

The `scripts/run_nonnative_entanglement_clustering.py` script automates the clustering workflow. For production runs with all 1000 trajectories, submit via SLURM using `assets/slurm/scripts/nonNativeClustering.slurm`.

```bash
source ~/.bashrc
conda activate entdetect

DATASTORE=/scratch/ims86/EntDetect_Datastore

mkdir -p $DATASTORE/outputs/workflow2/nonnative_clustering
mkdir -p $DATASTORE/outputs/workflow2/nonnative_clustering/logs

python scripts/run_nonnative_entanglement_clustering.py \
    --outdir               $DATASTORE/outputs/workflow2/nonnative_clustering \
    --trajnum2pklfile_path $DATASTORE/user_input/metadata/trajnum2file.txt \
    --traj_dir_prefix      $DATASTORE/user_input/cg_trajectories \
    --start_frame          6600 \
    --end_frame            6667 \
    --logdir               $DATASTORE/outputs/workflow2/nonnative_clustering/logs \
    --nproc                4 \
    --log_level            DEBUG
```

For full production clustering (all frames, all trajectories):

```bash
sbatch assets/slurm/scripts/nonNativeClustering.slurm
```

---

## Step 6. Build a Markov state model (MSM)

Organize simulation frames into microstates and metastable states using the full-trajectory order-parameter data using parameters from `assets/slurm/scripts/run_MSM.slurm`.

### 6a. Build the MSM with Python

```python
from EntDetect.clustering import MSMNonNativeEntanglementClustering

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow2"

outdir         = f"{OUTDIR}/MSM"
ID             = "1ZMR_prod"

# OPpath must point to a directory containing Q/, G/, K/ subdirectories
# with per-trajectory OP files for ALL trajectories.
OPpath         = f"{OUTDIR}/OP/"
n_large_states = 10   # number of metastable macro-states requested
lagtime        = 20   # lag time in frames

MSM = MSMNonNativeEntanglementClustering(
    outdir=outdir,
    ID=ID,
    OPpath=OPpath,
    start=0,
    n_large_states=n_large_states,
    lagtime=lagtime
)
MSM.run()
```

### 6b. Using the command-line interface

For convenience, use the `scripts/run_MSM.py` script directly (see "Running MSM construction as a single script" below for full details).

### 6c. Critical notes

- Try **multiple values** of `n_large_states` (e.g., 5, 10, 15). The protocol notes ≤15 often works well.
- The final number of states may be lower than requested if empty states are discarded.
- The default lag time of 1 frame is appropriate for visualization and exploratory grouping. For kinetic interpretation, test for Markovian behavior explicitly and choose lag time accordingly.

### 6d. Expected outputs

| File | Contents |
|------|----------|
| `1ZMR_prod_MSMmapping.csv` | Per-frame microstate and metastable-state assignments |
| `1ZMR_prod_meta_set.csv` | Metastable-state summary |
| `1ZMR_prod_meta_dist.npy` | Metastable-state probability distribution |
| `1ZMR_prod_StateAndFEplot.png` | Order-parameter landscape and state assignments |

---

## Running MSM construction as a single script

The `scripts/run_MSM.py` script automates MSM construction. For production runs, submit via SLURM using `assets/slurm/scripts/run_MSM.slurm`.

```bash
source ~/.bashrc
conda activate entdetect

DATASTORE=/scratch/ims86/EntDetect_Datastore

mkdir -p $DATASTORE/outputs/workflow2/MSM
mkdir -p $DATASTORE/outputs/workflow2/MSM/logs

python scripts/run_MSM.py \
    --outdir         $DATASTORE/outputs/workflow2/MSM \
    --ID             1ZMR_prod \
    --OPpath         $DATASTORE/outputs/workflow2/OP/ \
    --start          0 \
    --n_large_states 10 \
    --lagtime        20 \
    --logdir         $DATASTORE/outputs/workflow2/MSM/logs
```

For full production MSM:

```bash
sbatch assets/slurm/scripts/run_MSM.slurm
```

---

## Step 7. Analyze metastable-state behavior

### 7a. Plot state probability evolution over simulation time

```python
from EntDetect.statistics import MSMStats

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow2"

# Point to your Step 6 (MSM) outputs
msm_meta_file  = f"{OUTDIR}/MSM/1ZMR_prod_MSMmapping_A80pctNative.csv"
meta_set_file  = f"{OUTDIR}/MSM/1ZMR_prod_meta_set.csv"

outdir         = f"{OUTDIR}/MSM_StateProbabilityStats"
traj_type_col  = "traj_type_A80pctNative"
traj_type_list = ['A', 'B']    # trajectory type labels present in the MSM file
rm_traj_list   = []             # trajectories excluded in Step 4 (none here)

MS = MSMStats(
    outdir=outdir,
    msm_data_file=msm_meta_file,
    meta_set_file=meta_set_file,
    traj_type_col=traj_type_col,
    rm_traj_list=rm_traj_list,
    traj_type_list=traj_type_list
)

df = MS.StateProbabilityStats()
MS.Plot_StateProbabilityStats(df=df)
```

**Output:** A time series describing the population (probability) of each metastable state over the simulation. Useful for identifying which states are transiently vs persistently populated.

### 7b. Visualize representative metastable structures

Identify the representative frame for each metastable state from `1ZMR_prod_MSMmapping.csv`, extract the corresponding structure from the trajectory, and inspect it in VMD.

---

## Step 8. Compute folding pathway statistics

Analyse how trajectories transition between metastable states and quantify the divergence between trajectory-type populations using Jensen-Shannon (JS) divergence.

### 8a. Annotate the MSM mapping with trajectory-type labels

The MSM mapping CSV produced in Step 6 does not contain trajectory-type labels. Add a column classifying each trajectory based on a Q threshold. Trajectories that ever reach ≥ 80 % of native contacts are labelled 'A' (native-like); the remainder are 'B':

```python
import pandas as pd

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow2"

msm_mapping = pd.read_csv(f"{OUTDIR}/MSM/1ZMR_prod_MSMmapping.csv")

# A = trajectories whose max Q >= 0.80; B = all others
max_q = msm_mapping.groupby('traj')['Q'].transform('max')
msm_mapping['traj_type_A80pctNative'] = max_q.ge(0.80).map({True: 'A', False: 'B'})

annotated_file = f"{OUTDIR}/MSM/1ZMR_prod_MSMmapping_A80pctNative.csv"
msm_mapping.to_csv(annotated_file, index=False)
print(msm_mapping.groupby('traj')['traj_type_A80pctNative'].first().value_counts())
```

> **Note:** The 80 % threshold is system-specific. Adjust to match the Q distribution of your system.

### 8b. Compute folding pathways and Jensen-Shannon divergence with Python

```python
from EntDetect.statistics import FoldingPathwayStats

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow2"

msm_data      = pd.read_csv(f"{OUTDIR}/MSM/1ZMR_prod_MSMmapping_A80pctNative.csv")
meta_set_file = f"{OUTDIR}/MSM/1ZMR_prod_meta_set.csv"
outdir        = f"{OUTDIR}/FoldingPathway_A80pctNative"

FP = FoldingPathwayStats(
    msm_data=msm_data,
    meta_set_file=meta_set_file,
    tarj_type_col='traj_type_A80pctNative',
    traj_type_list=['A', 'B'],
    outdir=outdir,
    rm_traj_list=[],   # trajectories excluded in Step 4 (none here)
)

folding_pathways = FP.post_trans()
JS_divergence    = FP.JS_divergence()
```

`post_trans()` traces state-to-state transitions for each trajectory, removing loops to yield the minimal directed pathway. `JS_divergence()` computes a windowed Jensen-Shannon divergence between the A and B trajectory-type populations over simulation time.

| JS divergence | Interpretation |
|---------------|----------------|
| Near 0 | A and B explore similar state distributions |
| Near 1 | A and B have divergent state usage |

### 8c. Using the command-line interface

For convenience, use the `scripts/run_Foldingpathway.py` script directly (see "Running folding pathway analysis as a single script" below for full details).

### 8d. Expected outputs

| File | Contents |
|------|----------|
| `FoldingPathways_metastablestate_A-B.csv` | Per-type folding pathway probabilities |
| `JS_div_metastablestate_A-B.dat` | Windowed JS divergence time series |

---

## Running folding pathway analysis as a single script

The `scripts/run_Foldingpathway.py` script computes both folding pathways and JS divergence in a single call. It requires the MSM mapping CSV to already contain the trajectory-type column (see Step 8a).

```bash
source ~/.bashrc
conda activate entdetect

DATASTORE=/scratch/ims86/EntDetect_Datastore

mkdir -p $DATASTORE/outputs/workflow2/FoldingPathway_A80pctNative
mkdir -p $DATASTORE/outputs/workflow2/FoldingPathway_A80pctNative/logs

python scripts/run_Foldingpathway.py \
    --msm_data_file $DATASTORE/outputs/workflow2/MSM/1ZMR_prod_MSMmapping_A80pctNative.csv \
    --meta_set_file $DATASTORE/outputs/workflow2/MSM/1ZMR_prod_meta_set.csv \
    --traj_type_col traj_type_A80pctNative \
    --traj_type_list A B \
    --outdir        $DATASTORE/outputs/workflow2/FoldingPathway_A80pctNative \
    --logdir        $DATASTORE/outputs/workflow2/FoldingPathway_A80pctNative/logs \
    --log_level     INFO
```

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `G()` hangs or is very slow | Normal for long trajectories; CG + Topoly is expensive | Reduce `nproc` if memory-limited; use a shorter demo trajectory |
| `SASA()` returns NaN for some frames | Bad backmapped AA structure or topology mismatch | Inspect suspect frames in VMD; these frames are filtered in Workflow 3 |
| Jwalk error: PDB not found | `xp_pdb` path incorrect or file missing | Verify `1zmr_model_clean.pdb` is in `$DATASTORE/user_input/reference_structures/` |
| Jwalk error: `freesasa` not found | Package not installed in env | `pip install freesasa` inside the `entdetect` conda env |
| `Combined_GE/` is empty | Full OP run not yet complete | Use pre-computed results from `$DATASTORE/outputs/workflow2/OP/G/Combined_GE/` |
| MSM produces fewer states than `n_large_states` | Empty states discarded | Normal; try a higher `n_large_states` |
| `secondary_struc_defs.txt` not found | File not in `$DATASTORE/user_input/reference_structures/` | Verify files were rsynced to DATASTORE |
| `domain_def.dat` not found | File not in `$DATASTORE/user_input/reference_structures/` | Verify files were rsynced to DATASTORE |

---

← [Workflow 1](workflow1_native_ncle.md) | [Back to Master Index](index.md) | Next → [Workflow 3: Sim-to-Experiment](workflow3_sim2exp.md)
