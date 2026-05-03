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
| Q, G, K, SASA, XP for one AA trajectory (nproc=10) | 2–6 hours |
| Q, G, K for all 1000 trajectories (cluster) | Hours to days |
| Non-native clustering | Hours (memory-intensive) |
| MSM construction | Minutes |
| MSM statistics / folding pathways | Minutes |

---

## Pre-computed outputs

All 1000 CG and all-atom trajectories have already been analyzed and their outputs are stored in the DATASTORE:

```
$DATASTORE/outputs/workflow2/
├── OP_demo/                        # CG order parameters (all 1000 trajectories)
│   ├── Q/1ZMR_Traj{N}.Q
│   ├── G/1ZMR_Traj{N}.G
│   │   └── Combined_GE/            # Per-trajectory entanglement pkl files
│   │       └── 1ZMR_traj{N}_GE.pkl
│   └── K/K_{N}_prod.dat
└── OP_demo_AA/                     # All-atom order parameters (all 1000 trajectories)
    ├── Q/1ZMR_Traj{N}.Q
    ├── G/1ZMR_Traj{N}.G
    ├── K/K_{N}_prod.dat
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
| All-atom PDB topology | `$REFSTRUCT/1zmr_model_clean.pdb` | Q/G/K/SASA/XP |
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
for pkl in $DATASTORE/outputs/workflow2/OP_demo/G/Combined_GE/*.pkl; do
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

mkdir -p $OUTDIR/OP_demo/G $OUTDIR/OP_demo/Q $OUTDIR/OP_demo/K \
         $OUTDIR/OP_demo_AA/G $OUTDIR/OP_demo_AA/Q $OUTDIR/OP_demo_AA/K \
         $OUTDIR/OP_demo_AA/SASA $OUTDIR/OP_demo_AA/XP \
         $OUTDIR/nonnative_clustering $OUTDIR/MSM \
         $OUTDIR/MSM_StateProbabilityStats
```

---

## Step 9. Compute order parameters Q, G, and K on the CG trajectory

This demo runs on trajectory 420 to verify the CG portion of the workflow end-to-end.

> **Production strategy:** All 1000 CG trajectories have already been analyzed on the cluster. Their outputs are in `$OUTDIR/OP_demo/`. Steps 11–13 use those pre-computed outputs directly.

### 9a. Initialize `CalculateOP` for the CG trajectory

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
outdir       = f"{OUTDIR}/OP_demo"

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

### 9b. Compute Q — fraction of native contacts

```python
Qdata_dict = CalcOP.Q()
```

`Q` measures how many of the native residue–residue contacts present in the reference structure are also present in each trajectory frame. A value near 1.0 indicates a native-like conformation.

**Output:** A `.Q` file in `$OUTDIR/OP_demo/Q/`

### 9c. Compute G — entanglement order parameter

```python
Gdata_dict = CalcOP.G(topoly=False, Calpha=True, CG=True, nproc=10)
```

`G` captures the fraction of native contacts that exhibit a **change in entanglement state** relative to the native structure.

| Argument | Value | Meaning |
|----------|-------|---------|
| `topoly` | `False` | Use GLN-only workflow (no Topoly linking numbers) |
| `Calpha` | `True` | Use Cα-defined contacts (appropriate for CG trajectories) |
| `CG` | `True` | Input trajectory is coarse-grained |
| `nproc` | `10` | Number of CPU cores to use |

**Output:** A `.G` file and per-frame entanglement metadata `.pkl` in `$OUTDIR/OP_demo/G/`

> **Runtime note:** `G` is the most expensive order parameter to compute. Expect 12–20 hours per CG trajectory at ~6700 frames on 10 cores. Submit via the cluster for a full 1000-trajectory run (see [Running the full script](#running-the-full-workflow-as-a-single-script)).

### 9d. Compute K — mirror-symmetry order parameter

```python
Kdata_dict = CalcOP.K()
```

`K` detects frames where the protein has adopted a **mirror-image conformation** relative to the native structure. These frames are artifacts that must be removed before clustering.

**Output:** A `K_*.dat` file in `$OUTDIR/OP_demo/K/`

---

## Steps 9e–9h. Q, G, K, SASA, and XP on the all-atom trajectory

Q, G, and K can also be computed on all-atom trajectories. SASA and XP additionally require full atomic detail so a second `CalculateOP` instance is initialized from the all-atom trajectory.

> **Why separate?** The CG Cα-only representation does not carry enough atomic detail for accurate SASA (Shrake-Rupley requires explicit side-chain/backbone atoms) or cross-link SASD (Jwalk uses solvent-accessible surface geometry). The AA trajectories are produced by back-mapping the CG frames (see Workflow 3, Step 15).

### 9e. Initialize `CalculateOP` for the all-atom trajectory

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
    outdir=f"{OUTDIR}/OP_demo_AA",
    Traj=Traj,
    ID=ID,
    psf=AA_PDB,
    cor=AA_PDB,
    dcd=AA_DCD,
    sec_elements=f"{REFSTRUCT}/secondary_struc_defs.txt",
    domain=f"{REFSTRUCT}/domain_def.dat",
    start=0,
    ent_detection_method=2,   # 2 = any TLN
)
```

### 9f. Compute Q, G, K on the all-atom trajectory

```python
Qdata_dict_AA = CalcOP.Q()
Gdata_dict_AA = CalcOP.G(topoly=False, Calpha=False, CG=False, nproc=10)
Kdata_dict_AA = CalcOP.K()
```

> **Note:** `Calpha=False` uses heavy-atom contacts appropriate for all-atom structures. `CG=False` indicates an all-atom input.

### 9g. Compute SASA — solvent-accessible surface area per residue

```python
SASAdata_dict = CalcOP.SASA()
```

Uses the Shrake-Rupley algorithm (mdtraj, `probe_radius=0.14 nm`, 1000 sphere points) to compute per-residue SASA for every frame. Used downstream (Workflow 3) to test LiP-MS signals: high SASA residues are more accessible to the protease.

| Key | Contents |
|-----|----------|
| `outfile` | Path to `{ID}.SASA` CSV written to `$OUTDIR/OP_demo_AA/SASA/` |
| `result`  | DataFrame with columns `Time(ns)`, `Frame`, `resid`, `SASA(nm^2)` |

### 9h. Compute XP — Jwalk cross-link probability

```python
XPdata_dict = CalcOP.XP(pdb=AA_PDB)
```

`XP` computes the **solvent-accessible surface distance (SASD)** between all pairs of cross-linkable residue types (K, S, T, Y, M) and converts each SASD to a cross-link probability score using a Gaussian parameterised on K–K linker geometry. Used downstream to test XL-MS signals.

**Output:** Per-residue-pair XP scores in `$OUTDIR/OP_demo_AA/XP/Jwalk_results/{stem}_crosslink_list.txt`

---

## Running the full workflow as a single script

All order parameters (Q, G, K for CG; Q, G, K, SASA, XP for AA) are computed by `scripts/run_OP_on_simulation_traj.py`. This is the script used for the full 1000-trajectory production run.

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
    --outdir       $DATASTORE/outputs/workflow2/OP_demo \
    --logdir       $DATASTORE/outputs/workflow2/OP_demo/logs \
    --start        0 \
    --ent_detection_method 1 \
    --nproc        10 \
    --ops Q G K \
    --no_topoly
```

### AA run — Q, G, K, SASA, XP

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
    --outdir $DATASTORE/outputs/workflow2/OP_demo_AA \
    --logdir $DATASTORE/outputs/workflow2/OP_demo_AA/logs \
    --start  0 \
    --ent_detection_method 2 \
    --nproc        10 \
    --xp_pdb $REFSTRUCT/1zmr_model_clean.pdb \
    --ops Q G K SASA XP
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
| `--ops` | `Q G K` | `Q G K SASA XP` | OPs computed |

---

## Step 10. Identify and remove artificial mirror conformations

Before clustering entanglement changes, remove trajectories or frames that are mirror-image artifacts.

### Recommended procedure

1. Load Q and K time series for all trajectories.
2. Identify trajectories where K is persistently high and Q is elevated or anomalous.
3. Visually inspect flagged frames in VMD.
4. Record the trajectory numbers that are confirmed mirrors in a list (used as `rm_traj_list` in later steps).

```python
import numpy as np
import pandas as pd
import glob

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow2"

# Load Q and K time series for all trajectories
Q_files = sorted(glob.glob(f"{OUTDIR}/OP_demo/Q/*.Q"))
K_files = sorted(glob.glob(f"{OUTDIR}/OP_demo/K/*.dat"))

# Example: inspect one trajectory
Q_data = pd.read_csv(Q_files[0])
K_data = pd.read_csv(K_files[0], header=None)

print(Q_data.describe())
print(K_data.describe())
```

> **Important:** The cutoff values for flagging mirror conformations must be tuned for your system by examining the Q and K distributions. Do **not** transfer cutoffs unchanged from a different protein.

### Result for the 1ZMR example

After inspection, no trajectories were excluded from the ecPGK production run:
```python
rm_traj_list = []
```

---

## Step 11. Cluster non-native entanglement changes

This step identifies non-redundant changes in entanglement topology across all trajectories. It requires the full set of per-trajectory G pkl files in `Combined_GE/` (present in the DATASTORE for all 1000 trajectories).

> **Memory warning:** This step can require tens of gigabytes of RAM. Run on a high-memory node.

### 11a. Build the clustering with Python

```python
from EntDetect.clustering import ClusterNonNativeEntanglements

DATASTORE            = "/scratch/ims86/EntDetect_Datastore"
OUTDIR               = f"{DATASTORE}/outputs/workflow2"
CG_TRAJ_DIR          = f"{DATASTORE}/user_input/cg_trajectories"

pkl_file_path        = f"{OUTDIR}/OP_demo/G/Combined_GE/"
trajnum2pklfile_path = f"{DATASTORE}/user_input/metadata/trajnum2file.txt"
traj_dir_prefix      = CG_TRAJ_DIR
outdir               = f"{OUTDIR}/nonnative_clustering"

clustering_NNents = ClusterNonNativeEntanglements(
    pkl_file_path=pkl_file_path,
    trajnum2pklfile_path=trajnum2pklfile_path,
    traj_dir_prefix=traj_dir_prefix,
    outdir=outdir,
)

# start_frame: matches the --start value used when computing G (0 here)
# end_frame: omit or set to a large number to include all frames
clustering_NNents.cluster(start_frame=0, end_frame=9999999)
```

### 11b. Run via the analysis script

```bash
DATASTORE=/scratch/ims86/EntDetect_Datastore

python scripts/run_nonnative_entanglement_clustering.py \
    --outdir               $DATASTORE/outputs/workflow2/nonnative_clustering \
    --pkl_file_path        $DATASTORE/outputs/workflow2/OP_demo/G/Combined_GE/ \
    --trajnum2pklfile_path $DATASTORE/user_input/metadata/trajnum2file.txt \
    --traj_dir_prefix      $DATASTORE/user_input/cg_trajectories \
    --start_frame          0 \
    --logdir               $DATASTORE/outputs/workflow2/nonnative_clustering/logs
```

### Inspecting clustering results

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

### Expected outputs

| File | Contents |
|------|----------|
| `rep_chg_ent_topoly_linking_number.csv` | Representative entanglement changes (one per cluster) |
| `chg_ent_struct_topoly_linking_number.csv` | Per-frame cluster assignment |
| `cluster_data_topoly_linking_number.npz` | Compressed cluster data array |
| `cluster_tree_topoly_linking_number.dat` | Text representation of the clustering hierarchy |
| `rep_chg_ent_list_topoly_linking_number.pkl` | List of representative entanglement objects |
| `chg_ent_topoly_linking_number_distribution.pdf` | Distribution plot of loop/crossing residues |

---

## Step 12. Build a Markov state model (MSM)

Organize simulation frames into microstates and metastable states using the full-trajectory order-parameter data.

```python
from EntDetect.clustering import MSMNonNativeEntanglementClustering

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow2"

outdir         = f"{OUTDIR}/MSM"
ID             = "1ZMR"

# OPpath must point to a directory containing Q/, G/, K/ subdirectories
# with per-trajectory OP files for ALL trajectories.
OPpath         = f"{OUTDIR}/OP_demo/"
n_large_states = 10   # number of metastable macro-states requested

MSM = MSMNonNativeEntanglementClustering(
    outdir=outdir,
    ID=ID,
    OPpath=OPpath,
    n_large_states=n_large_states
)
MSM.run()
```

### Run via the analysis script

```bash
DATASTORE=/scratch/ims86/EntDetect_Datastore

python scripts/run_MSM.py \
    --outdir  $DATASTORE/outputs/workflow2/MSM \
    --ID      1ZMR_prod \
    --OPpath  $DATASTORE/outputs/workflow2/OP_demo/ \
    --start   0 \
    --n_large_states 10 \
    --lagtime 20
```

### Critical notes

- Try **multiple values** of `n_large_states` (e.g., 5, 10, 15). The protocol notes ≤15 often works well.
- The final number of states may be lower than requested if empty states are discarded.
- The default lag time of 1 frame is appropriate for visualization and exploratory grouping. For kinetic interpretation, test for Markovian behavior explicitly and choose lag time accordingly.

### Expected outputs

| File | Contents |
|------|----------|
| `1ZMR_prod_MSMmapping.csv` | Per-frame microstate and metastable-state assignments |
| `1ZMR_prod_meta_set.csv` | Metastable-state summary |
| `1ZMR_prod_meta_dist.npy` | Metastable-state probability distribution |
| `1ZMR_prod_StateAndFEplot.png` | Order-parameter landscape and state assignments |

---

## Step 13. Analyze metastable-state behavior

### 13a. Plot state probability evolution over simulation time

```python
from EntDetect.statistics import MSMStats

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow2"

# Point to your Step 12 outputs
msm_meta_file  = f"{OUTDIR}/MSM/1ZMR_prod_meta_set_A80%Native.csv"
meta_set_file  = f"{OUTDIR}/MSM/1ZMR_prod_meta_set.csv"

outdir         = f"{OUTDIR}/MSM_StateProbabilityStats"
traj_type_col  = "traj_type_A80%Native"
traj_type_list = ['A', 'B']    # trajectory type labels present in the MSM file
rm_traj_list   = []             # trajectories excluded in Step 10 (none here)

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

### 13b. Visualize representative metastable structures

Identify the representative frame for each metastable state from `1ZMR_prod_MSMmapping.csv`, extract the corresponding structure from the trajectory, and inspect it in VMD.

### 13c. Compare folding pathways and Jensen-Shannon divergence

```python
import pandas as pd
from EntDetect.statistics import FoldingPathwayStats

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow2"

msm_meta_file = f"{OUTDIR}/MSM/1ZMR_prod_meta_set_A80%Native.csv"
meta_set_file = f"{OUTDIR}/MSM/1ZMR_prod_meta_set.csv"
traj_type_col = "traj_type_A80%Native"
rm_traj_list  = []
outdir        = f"{OUTDIR}/Foldingpathway_A80%Native"

msm_data = pd.read_csv(msm_meta_file)

FP = FoldingPathwayStats(
    msm_data=msm_data,
    meta_set_file=meta_set_file,
    traj_type_col=traj_type_col,
    outdir=outdir,
    traj_list=rm_traj_list
)

folding_pathways = FP.post_trans()
JS_divergence    = FP.JS_divergence()

print("Folding pathways:")
print(folding_pathways)
print("\nJensen-Shannon divergence time series:")
print(JS_divergence)
```

### Interpreting Jensen-Shannon divergence

| JSD value | Interpretation |
|-----------|----------------|
| Near 0 | Conditions A and B explore similar state distributions |
| Near 1 | Conditions A and B have divergent state usage |

Compare output to:
```
$OUTDIR/Foldingpathway_A80%Native/FoldingPathways_metastablestate_A-B.csv
$OUTDIR/Foldingpathway_A80%Native/JS_div_metastablestate_A-B.dat
```

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `G()` hangs or is very slow | Normal for long trajectories; CG + Topoly is expensive | Reduce `nproc` if memory-limited; use a shorter demo trajectory |
| `SASA()` returns NaN for some frames | Bad backmapped AA structure or topology mismatch | Inspect suspect frames in VMD; these frames are filtered in Workflow 3 |
| Jwalk error: PDB not found | `xp_pdb` path incorrect or file missing | Verify `1zmr_model_clean.pdb` is in `$DATASTORE/user_input/reference_structures/` |
| Jwalk error: `freesasa` not found | Package not installed in env | `pip install freesasa` inside the `entdetect` conda env |
| `Combined_GE/` is empty | Full OP run not yet complete | Use pre-computed results from `$DATASTORE/outputs/workflow2/OP_demo/G/Combined_GE/` |
| MSM produces fewer states than `n_large_states` | Empty states discarded | Normal; try a higher `n_large_states` |
| `secondary_struc_defs.txt` not found | File not in `$DATASTORE/user_input/reference_structures/` | Verify files were rsynced to DATASTORE |
| `domain_def.dat` not found | File not in `$DATASTORE/user_input/reference_structures/` | Verify files were rsynced to DATASTORE |

---

← [Workflow 1](workflow1_native_ncle.md) | [Back to Master Index](index.md) | Next → [Workflow 3: Sim-to-Experiment](workflow3_sim2exp.md)
