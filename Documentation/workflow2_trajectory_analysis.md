# Workflow 2: Detect Changes in Entanglement Across Simulation Trajectories

← [Back to Master Index](index.md)

---

## Goal

Compute order parameters (Q, G, K) across simulation trajectories, cluster non-native entanglement changes, build a Markov state model (MSM), and analyze metastable state behavior.

---

## Typical runtime

| Step | Runtime |
|------|---------|
| Q, G, K for one trajectory (single demo) | 10–60 minutes |
| Q, G, K for all 1000 trajectories | Hours to days |
| Non-native clustering | Hours (memory-intensive) |
| MSM construction | Minutes |
| MSM statistics / folding pathways | Minutes |

> **Strategy for this tutorial:** Steps 9–10 are demonstrated on a single trajectory (`420_prod.dcd`) to verify the workflow is functional. The full-run outputs for all 1000 trajectories are pre-computed and stored in `$TESTDIR`. Steps 11–13 are run directly from those pre-computed outputs.

---

## Pre-computed reference outputs

```
$TESTDIR/
├── OP_Full/                  # OPs computed over all frames (1000 trajectories)
│   ├── Q/1ZMR_prod_Traj1.Q
│   ├── G/1ZMR_prod_Traj1.G
│   └── K/K_1_prod.dat
├── OP_last67_density1/       # OPs computed over last 67 frames only
│   ├── Q/
│   ├── G/
│   └── K/
├── nonnative_entanglement_clustering/
│   ├── chg_ent_struct_topoly_linking_number.csv
│   ├── chg_ent_topoly_linking_number_distribution.pdf
│   ├── cluster_data_topoly_linking_number.npz
│   ├── cluster_tree_topoly_linking_number.dat
│   ├── rep_chg_ent_list_topoly_linking_number.pkl
│   └── rep_chg_ent_topoly_linking_number.csv
└── MSM/
    ├── 1ZMR_prod_MSMmapping.csv
    ├── 1ZMR_prod_meta_dist.npy
    ├── 1ZMR_prod_meta_set.csv
    ├── 1ZMR_prod_meta_set_A80%Native.csv
    ├── 1ZMR_prod_meta_set_random.csv
    └── 1ZMR_prod_StateAndFEplot.png
```

Where:
```bash
TESTDIR=/storage/group/epo2/default/ims86/git_repos/EntanGoPy/TestingGrounds
```

Tutorial outputs from re-running steps are written to:
```bash
DATASTORE=/scratch/ims86/EntDetect_Datastore
OUTDIR=$DATASTORE/outputs/workflow2
```

---

## Required input files

| File | Path | Notes |
|------|------|-------|
| Cα PSF topology | `$DATASTORE/user_input/reference_structures/1zmr_model_clean_ca.psf` | |
| Cα COR reference | `$DATASTORE/user_input/reference_structures/1zmr_model_clean_ca.cor` | |
| CG trajectory (demo) | `$CG_TRAJ_DIR/420_prod.dcd` | Single trajectory for testing |
| CG trajectories (all) | `$CG_TRAJ_DIR/{N}_prod.dcd` (N=1–1000) | Full production run |
| Secondary structure defs | See warning below | **Must be provided** |
| Domain boundary file | See warning below | **Must be provided** |
| Trajectory-to-file mapping | See warning below | **Must be provided** |

```bash
CG_TRAJ_DIR=$DATASTORE/user_input/cg_trajectories
```

> **Metadata files location:** The following files are required for the full workflow. `secondary_struc_defs.txt` and `domain_def.dat` are in `$DATASTORE/user_input/reference_structures/`. The `trajnum2file.txt` file is in `$DATASTORE/user_input/metadata/`.

### `trajnum2file.txt` — already available

This file maps trajectory numbers to their **G-order-parameter pkl files** (not DCD files). The format is comma-separated:

```
<trajectory_number>,<path_to_pkl_file>
```

A pre-populated copy has been placed at:

```
$DATASTORE/user_input/metadata/trajnum2file.txt
```

The paths inside point to the pre-computed pkl files in the backup:
```
/storage/group/epo2/default/ims86/git_repos/EntanGoPy_bak/TestingGrounds/FinalTesting_v2.1/OP_last67_density1/G/Combined_GE/
```

All 1000 pkl files are present there. No regeneration is needed.

---

## Step 1. Activate your environment and set paths

```bash
source ~/.bashrc
conda activate entdetect

DATASTORE=/scratch/ims86/EntDetect_Datastore
TESTDIR=/storage/group/epo2/default/ims86/git_repos/EntanGoPy/TestingGrounds
CG_TRAJ_DIR=$DATASTORE/user_input/cg_trajectories
REFSTRUCT=$DATASTORE/user_input/reference_structures
OUTDIR=$DATASTORE/outputs/workflow2

mkdir -p $OUTDIR/OP_demo/G $OUTDIR/OP_demo/Q $OUTDIR/OP_demo/K \
         $OUTDIR/nonnative_clustering $OUTDIR/MSM \
         $OUTDIR/MSM_StateProbabilityStats \
         "$OUTDIR/Foldingpathway_A80%Native"
```

---

## Step 9. Compute order parameters Q, G, and K (single-trajectory demo)

This demo runs on trajectory 420 to verify the workflow end-to-end without committing to a full 1000-trajectory run.

### 9a. Initialize `CalculateOP`

```python
from EntDetect.order_params import CalculateOP

DATASTORE   = "/scratch/ims86/EntDetect_Datastore"
TESTDIR     = "/storage/group/epo2/default/ims86/git_repos/EntanGoPy/TestingGrounds"
REFSTRUCT   = f"{DATASTORE}/user_input/reference_structures"
CG_TRAJ_DIR = f"{DATASTORE}/user_input/cg_trajectories"
OUTDIR      = f"{DATASTORE}/outputs/workflow2"

Traj        = 420
PSF         = f"{REFSTRUCT}/1zmr_model_clean_ca.psf"
COR         = f"{REFSTRUCT}/1zmr_model_clean_ca.cor"
DCD         = f"{CG_TRAJ_DIR}/420_prod.dcd"
ID          = "1ZMR"
sec_elements = f"{REFSTRUCT}/secondary_struc_defs.txt"
domain      = f"{REFSTRUCT}/domain_def.dat"
outdir      = f"{OUTDIR}/OP_demo"

# start: first frame to include (0-indexed; 0 means the very first frame).
# For a production run, the protocol uses start=6600 to skip the first
# portion of the trajectory. For this demo we start from frame 0.
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
    start=start
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
Gdata_dict = CalcOP.G(topoly=True, Calpha=True, CG=True, nproc=10)
```

`G` captures the fraction of native contacts that exhibit a **change in entanglement state** relative to the native structure.

| Argument | Value | Meaning |
|----------|-------|---------|
| `topoly` | `True` | Use Topoly to compute linking numbers |
| `Calpha` | `True` | Use Cα-defined contacts (appropriate for CG trajectories) |
| `CG` | `True` | Input trajectory is coarse-grained |
| `nproc` | `10` | Number of CPU cores to use |

**Output:** A `.G` file and per-frame entanglement metadata `.pkl` in `$OUTDIR/OP_demo/G/`

> **Memory and runtime note:** `G` is the most expensive order parameter to compute. For a full production run across 1000 trajectories, this step is parallelized across a cluster. The demo on a single trajectory should complete within 10–60 minutes depending on trajectory length and available cores.

### 9d. Compute K — mirror-symmetry order parameter

```python
Kdata_dict = CalcOP.K()
```

`K` detects frames where the protein has adopted a **mirror-image conformation** relative to the native structure. These frames are artifacts that must be removed before clustering.

**Output:** A `K_*.dat` file in `$OUTDIR/OP_demo/K/`

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

TESTDIR = "/storage/group/epo2/default/ims86/git_repos/EntanGoPy/TestingGrounds"

# Load pre-computed full Q and K outputs
q_file = f"{TESTDIR}/OP_Full/Q/1ZMR_prod_Traj1.Q"
k_file = f"{TESTDIR}/OP_Full/K/K_1_prod.dat"

Q_data = pd.read_csv(q_file)
K_data = pd.read_csv(k_file, header=None)

# Inspect the distributions
print(Q_data.describe())
print(K_data.describe())
```

> **Important:** The cutoff values for flagging mirror conformations in the published ecPGK protocol were tuned specifically for that system. Do **not** transfer them unchanged to a different protein. Derive appropriate cutoffs by examining the distributions for your own data.

### Result for the 1ZMR example

After mirror-artifact removal, the list of excluded trajectories is supplied as `rm_traj_list = []` in the downstream steps — meaning no trajectories were excluded from the ecPGK production run in this example.

---

## Step 11. Cluster non-native entanglement changes

This step identifies non-redundant changes in entanglement across all trajectories.

> **This step requires the full-run Combined_GE pkl files**, which are generated by running Step 9 on all 1000 trajectories and then combining the per-trajectory `.G` metadata into a single `Combined_GE/` directory. Because generating these files requires the full compute run, this step is demonstrated using the pre-computed clustering results already present in `$TESTDIR/nonnative_entanglement_clustering/`.

### Re-running from scratch (full production run only)

If `Combined_GE/` has been populated by running Step 9 across all trajectories, create the `trajnum2file.txt` mapping file (see above) and run:

```python
from EntDetect.clustering import ClusterNonNativeEntanglements

DATASTORE   = "/scratch/ims86/EntDetect_Datastore"
TESTDIR     = "/storage/group/epo2/default/ims86/git_repos/EntanGoPy/TestingGrounds"
OUTDIR      = f"{DATASTORE}/outputs/workflow2"
CG_TRAJ_DIR = f"{DATASTORE}/user_input/cg_trajectories"

# pkl_file_path: directory containing per-trajectory G pkl files
    pkl_file_path        = "/storage/group/epo2/default/ims86/git_repos/EntanGoPy_bak/TestingGrounds/FinalTesting_v2.1/OP_last67_density1/G/Combined_GE/"
    trajnum2pklfile_path = f"{DATASTORE}/user_input/metadata/trajnum2file.txt"
traj_dir_prefix      = CG_TRAJ_DIR
outdir               = f"{OUTDIR}/nonnative_clustering"

clustering_NNents = ClusterNonNativeEntanglements(
    pkl_file_path=pkl_file_path,
    trajnum2pklfile_path=trajnum2pklfile_path,
    traj_dir_prefix=traj_dir_prefix,
    outdir=outdir
)

# start_frame: matches the `start` value used when computing G
clustering_NNents.cluster(start_frame=6600)
```

### Inspecting pre-computed clustering results

The pre-computed outputs in `$TESTDIR/nonnative_entanglement_clustering/` can be examined immediately:

```python
import pandas as pd
import pickle
import numpy as np

TESTDIR = "/storage/group/epo2/default/ims86/git_repos/EntanGoPy/TestingGrounds"
clust_dir = f"{TESTDIR}/nonnative_entanglement_clustering"

# Representative entanglement changes
rep_df = pd.read_csv(f"{clust_dir}/rep_chg_ent_topoly_linking_number.csv")
print(f"Number of representative entanglement changes: {len(rep_df)}")
print(rep_df.head())

# Structural assignments for each frame
chg_df = pd.read_csv(f"{clust_dir}/chg_ent_struct_topoly_linking_number.csv")
print(chg_df.head())

# Cluster data array
cluster_data = np.load(f"{clust_dir}/cluster_data_topoly_linking_number.npz", allow_pickle=True)
print("Available arrays:", list(cluster_data.keys()))
```

### Expected outputs (from a full run)

| File | Contents |
|------|----------|
| `rep_chg_ent_topoly_linking_number.csv` | Representative entanglement changes (one per cluster) |
| `chg_ent_struct_topoly_linking_number.csv` | Per-frame cluster assignment |
| `cluster_data_topoly_linking_number.npz` | Compressed cluster data array |
| `cluster_tree_topoly_linking_number.dat` | Text representation of the clustering hierarchy |
| `rep_chg_ent_list_topoly_linking_number.pkl` | List of representative entanglement objects |
| `chg_ent_topoly_linking_number_distribution.pdf` | Distribution plot of loop/crossing residues |

> **Memory warning:** This step can require tens of gigabytes of RAM if the number of raw entanglement-change events is large. Monitor memory usage and consider running on a high-memory node.

---

## Step 12. Build a Markov state model (MSM)

Organize simulation frames into microstates and metastable states using the full-trajectory order-parameter data.

```python
from EntDetect.clustering import MSMNonNativeEntanglementClustering

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
TESTDIR   = "/storage/group/epo2/default/ims86/git_repos/EntanGoPy/TestingGrounds"
OUTDIR    = f"{DATASTORE}/outputs/workflow2"

outdir         = f"{OUTDIR}/MSM"
ID             = "1ZMR"

# OPpath must point to a directory containing Q/, G/, K/ subdirectories
# with per-trajectory OP files for ALL trajectories. Use pre-computed
# OP_Full data here.
OPpath         = f"{TESTDIR}/OP_Full/"
n_large_states = 10   # number of metastable macro-states requested

MSM = MSMNonNativeEntanglementClustering(
    outdir=outdir,
    ID=ID,
    OPpath=OPpath,
    n_large_states=n_large_states
)
MSM.run()
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

Compare to pre-computed references in `$TESTDIR/MSM/`.

---

## Step 13. Analyze metastable-state behavior

### 13a. Plot state probability evolution over simulation time

```python
from EntDetect.statistics import MSMStats

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
TESTDIR   = "/storage/group/epo2/default/ims86/git_repos/EntanGoPy/TestingGrounds"
OUTDIR    = f"{DATASTORE}/outputs/workflow2"

# Point to pre-computed MSM files (or your Step 12 outputs)
msm_meta_file  = f"{TESTDIR}/MSM/1ZMR_prod_meta_set_A80%Native.csv"
meta_set_file  = f"{TESTDIR}/MSM/1ZMR_prod_meta_set.csv"

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
TESTDIR   = "/storage/group/epo2/default/ims86/git_repos/EntanGoPy/TestingGrounds"
OUTDIR    = f"{DATASTORE}/outputs/workflow2"

msm_meta_file = f"{TESTDIR}/MSM/1ZMR_prod_meta_set_A80%Native.csv"
meta_set_file = f"{TESTDIR}/MSM/1ZMR_prod_meta_set.csv"
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

Compare output to pre-computed references:
```
$TESTDIR/Foldingpathway_A80%Native/FoldingPathways_metastablestate_A-B.csv
$TESTDIR/Foldingpathway_A80%Native/JS_div_metastablestate_A-B.dat
```

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `G()` hangs or is very slow | Normal for long trajectories; CG + Topoly is expensive | Reduce `nproc` if memory-limited; use a shorter demo trajectory |
| `Combined_GE/` is empty | Full OP run not yet complete | Use pre-computed clustering results from `$TESTDIR` |
| MSM produces fewer states than `n_large_states` | Empty states discarded | Normal; try a higher `n_large_states` |
| `secondary_struc_defs.txt` not found | File not in `$DATASTORE/user_input/reference_structures/` | Verify files were rsynced to DATASTORE |
| `domain_def.dat` not found | File not in `$DATASTORE/user_input/reference_structures/` | Verify files were rsynced to DATASTORE |

---

← [Workflow 1](workflow1_native_ncle.md) | [Back to Master Index](index.md) | Next → [Workflow 3: Sim-to-Experiment](workflow3_sim2exp.md)
