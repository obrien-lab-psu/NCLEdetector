# Workflow 3: Compare Structural Ensembles to High-Throughput Experimental Data

← [Back to Master Index](index.md)

---

## Goal

Test whether specific metastable states from simulation are statistically consistent with experimental conformational signals (LiP-MS and XL-MS). Identify and extract representative structures from the best-supported metastable states.

---

## Typical runtime

| Step | Runtime |
|------|---------|
| Back-mapping (per frame) | 5–30 minutes per frame |
| SASA calculation | Hours |
| Jwalk distances | Hours |
| XP (cross-link probability) | Hours |
| LiP-MS / XL-MS consistency test | Minutes to hours |

> **Strategy for this tutorial:** SASA, Jwalk, and XP outputs are pre-computed and available in `$TESTDIR/compare_sim2exp/`. The consistency test (Step 17) is fully runnable from those pre-computed inputs.

---

## Pre-computed reference outputs

```
$TESTDIR/compare_sim2exp/
├── SASA.npy                                               # per-frame SASA data
├── Jwalk.npy                                              # per-frame Jwalk distances
├── ecPGK_significant_LiPMS_peptide_R1_merged.xlsx         # processed LiP-MS experimental data
├── ecPGK_significant_XLMS_peptide_R1_merged.xlsx          # processed XL-MS experimental data
├── LiPMS_XLMS_consist_pvalues_metastates_v11_...xlsx      # consistency test p-values
├── Consistent_structures_v8.xlsx                          # selected representative structures
├── consist_signal_struct_data.npz                         # raw consistency test arrays
├── LiPMS_XLMS_consist_data_v9.npz                        # consistency data arrays
└── native_M_data.npz                                      # native state structural data
```

MSM files also required:
```
$TESTDIR/MSM/
├── 1ZMR_prod_MSMmapping.csv
└── 1ZMR_prod_meta_dist.npy
```

Where:
```bash
TESTDIR=/storage/group/epo2/default/ims86/git_repos/EntanGoPy/TestingGrounds
```

Tutorial outputs are written to:
```bash
DATASTORE=/scratch/ims86/EntDetect_Datastore
OUTDIR=$DATASTORE/outputs/workflow3
```

---

## Required input files

| File | Path | Notes |
|------|------|-------|
| Native all-atom PDB | `$DATASTORE/user_input/reference_structures/1zmr_model_clean.pdb` | All-atom reference |
| Cα COR coordinates | `$DATASTORE/user_input/reference_structures/1zmr_model_clean_ca.cor` | Cα reference for back-mapping |
| AA trajectory DCDs | `$AA_TRAJ_DIR/{N}_prod_aa.dcd` (N=1–1000) | All-atom MD trajectories |
| SASA data | `$TESTDIR/compare_sim2exp/SASA.npy` | Pre-computed |
| Jwalk distances | `$TESTDIR/compare_sim2exp/Jwalk.npy` | Pre-computed |
| LiP-MS experimental file | `$DATASTORE/user_input/experimental_data/ecPGK_significant_LiPMS_peptide_R1_merged.xlsx` | |
| XL-MS experimental file | `$DATASTORE/user_input/experimental_data/ecPGK_significant_XLMS_peptide_R1_merged.xlsx` | |
| MSM mapping | `$TESTDIR/MSM/1ZMR_prod_MSMmapping.csv` | From Workflow 2 |
| MSM meta distribution | `$TESTDIR/MSM/1ZMR_prod_meta_dist.npy` | From Workflow 2 |
| Non-native clustering data | `$TESTDIR/nonnative_entanglement_clustering/cluster_data_topoly_linking_number.npz` | From Workflow 2 |

```bash
AA_TRAJ_DIR=$DATASTORE/user_input/aa_trajectories
```

---

## Step 1. Activate your environment and set paths

```bash
source ~/.bashrc
conda activate entdetect

DATASTORE=/scratch/ims86/EntDetect_Datastore
TESTDIR=/storage/group/epo2/default/ims86/git_repos/EntanGoPy/TestingGrounds
AA_TRAJ_DIR=$DATASTORE/user_input/aa_trajectories
OUTDIR=$DATASTORE/outputs/workflow3

mkdir -p $OUTDIR/BackMapping $OUTDIR/MassSpec_ConsistencyTest
```

---

## Step 14. Prepare experimental inputs

The processed experimental files are already available in `$TESTDIR/compare_sim2exp/`. For your own system, you would:

1. Run your LiP-MS statistical analysis pipeline to identify significantly changed peptides.
2. Map significant peptides back to residue-level assignments.
3. Format the output as an Excel file matching the expected column structure.
4. Repeat for XL-MS cross-link data.

> **Critical:** There are many valid ways to define statistical significance in LiP-MS data. The choice of significance threshold and mapping strategy can strongly affect downstream conclusions. Be consistent and document your decisions.

To verify the experimental input files are accessible:

```python
import pandas as pd

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
TESTDIR   = "/storage/group/epo2/default/ims86/git_repos/EntanGoPy/TestingGrounds"

lipms = pd.read_excel(f"{DATASTORE}/user_input/experimental_data/ecPGK_significant_LiPMS_peptide_R1_merged.xlsx")
xlms  = pd.read_excel(f"{DATASTORE}/user_input/experimental_data/ecPGK_significant_XLMS_peptide_R1_merged.xlsx")

print("LiP-MS shape:", lipms.shape)
print(lipms.head())
print("\nXL-MS shape:", xlms.shape)
print(xlms.head())
```

---

## Step 15. Back-map coarse-grained trajectories to all-atom structures

> **This step is computationally expensive and has already been completed** for the 1ZMR example system. The all-atom DCD trajectories in `$AA_TRAJ_DIR` are the result. This section describes the procedure for reference and for running on your own system.

Back-mapping reconstructs full all-atom models from Cα-only (coarse-grained) trajectory frames using the native all-atom structure as a template.

### 15a. Save individual Cα PDB frames from a trajectory

The frame extraction step is system-specific and typically done with MDAnalysis or VMD. Example with MDAnalysis:

```python
import MDAnalysis as mda

DATASTORE   = "/scratch/ims86/EntDetect_Datastore"
TESTDIR     = "/storage/group/epo2/default/ims86/git_repos/EntanGoPy/TestingGrounds"
CG_TRAJ_DIR = f"{DATASTORE}/user_input/cg_trajectories"
OUTDIR      = f"{DATASTORE}/outputs/workflow3"

import os
os.makedirs(f"{OUTDIR}/cg_frames", exist_ok=True)

psf = f"{DATASTORE}/user_input/reference_structures/1zmr_model_clean_ca.psf"
dcd = f"{CG_TRAJ_DIR}/420_prod.dcd"

u = mda.Universe(psf, dcd)
protein = u.select_atoms("protein")

for i, ts in enumerate(u.trajectory[-67:]):   # last 67 frames
    out_pdb = f"{OUTDIR}/cg_frames/frame_{i:04d}.pdb"
    protein.write(out_pdb)
```

### 15b. Back-map each Cα frame to all-atom

```python
from EntDetect.change_resolution import BackMapping

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
TESTDIR   = "/storage/group/epo2/default/ims86/git_repos/EntanGoPy/TestingGrounds"
OUTDIR    = f"{DATASTORE}/outputs/workflow3"

# cg_pdb : a single CG (Cα-only) PDB frame extracted above
cg_pdb  = f"{OUTDIR}/cg_frames/frame_0000.pdb"
aa_pdb  = f"{DATASTORE}/user_input/reference_structures/1zmr_model_clean.pdb"
ID      = "1ZMR"

backMapper = BackMapping(outdir=f"{OUTDIR}/BackMapping")
backMapper.backmap(cg_pdb=cg_pdb, aa_pdb=aa_pdb, ID=ID)
```

### What back-mapping produces

Depending on configuration, outputs include:

- reconstructed all-atom structures (`.pdb`);
- intermediate PD2 / Pulchra outputs;
- OpenMM energy-minimization logs and energy-minimized final structures.

> **After back-mapping:** Inspect a representative subset of reconstructed structures in VMD before proceeding. Collate the validated per-frame PDBs into an all-atom DCD for downstream use.

---

## Step 16. Calculate SASA, Jwalk distances, and cross-link probability

> **These outputs are pre-computed** for the 1ZMR system. The files `SASA.npy` and `Jwalk.npy` in `$TESTDIR/compare_sim2exp/` were generated using the last 67 frames (downsampled every 20 frames) of all 1000 trajectories. This section describes how to reproduce them on a single trajectory for verification.

### 16a. Initialize `CalculateOP` for all-atom data

```python
from EntDetect.order_params import CalculateOP

DATASTORE   = "/scratch/ims86/EntDetect_Datastore"
TESTDIR     = "/storage/group/epo2/default/ims86/git_repos/EntanGoPy/TestingGrounds"
AA_TRAJ_DIR = f"{DATASTORE}/user_input/aa_trajectories"
OUTDIR      = f"{DATASTORE}/outputs/workflow3"

Traj   = 420
# For all-atom OP, PSF is replaced with the all-atom PDB as topology
PSF    = f"{DATASTORE}/user_input/reference_structures/1zmr_model_clean.pdb"
DCD    = f"{AA_TRAJ_DIR}/420_prod_aa.dcd"
ID     = "1ZMR"
outdir = f"{OUTDIR}/OP_AA_demo"
start  = 6600   # skip first portion; match the value used in production

CalcOP = CalculateOP(
    outdir=outdir,
    Traj=Traj,
    ID=ID,
    psf=PSF,
    dcd=DCD,
    start=start
)
```

### 16b. Compute SASA

```python
CalcOP.SASA()
```

Computes the **solvent-accessible surface area** for each residue in each frame. Used to test LiP-MS signals (protease accessibility is correlated with SASA).

**Output:** A compressed `.npy` array in `$OUTDIR/OP_AA_demo/SASA/`

Compare to pre-computed reference: `$TESTDIR/compare_sim2exp/SASA.npy`

### 16c. Compute Jwalk distances

```python
# Provide the path to one representative all-atom back-mapped PDB for Jwalk initialization
CalcOP.runJwalk(f"{DATASTORE}/user_input/reference_structures/1zmr_model_clean.pdb")
```

Computes **solvent-accessible surface distances** between lysine residues, which are used to assess XL-MS cross-link compatibility.

**Output:** A compressed `.npy` array in `$OUTDIR/OP_AA_demo/Jwalk/`

Compare to pre-computed reference: `$TESTDIR/compare_sim2exp/Jwalk.npy`

### 16d. Compute XP — cross-link probability

```python
XPdata_dict = CalcOP.XP(pdb=f"{DATASTORE}/user_input/reference_structures/1zmr_model_clean.pdb")
```

Computes a **cross-linking probability metric** for each pair of reactive residues based on spatial distance distributions across the trajectory ensemble.

> **Downsampling guidance:** The production protocol analyzed only the last 50 ns and downsampled every 20 frames. Select a downsampling interval that reduces autocorrelation while retaining the structural states of interest. This choice may differ for other proteins or simulation setups.

---

## Step 17. Run the LiP-MS / XL-MS consistency test

This is the key integrative step. It tests whether specific metastable states identified in Workflow 2 are statistically consistent with the experimental signals.

### 17a. Initialize `MassSpec`

```python
from EntDetect.compare_sim2exp import MassSpec

DATASTORE   = "/scratch/ims86/EntDetect_Datastore"
TESTDIR     = "/storage/group/epo2/default/ims86/git_repos/EntanGoPy/TestingGrounds"
AA_TRAJ_DIR = f"{DATASTORE}/user_input/aa_trajectories"
OUTDIR      = f"{DATASTORE}/outputs/workflow3"

msm_data_file    = f"{TESTDIR}/MSM/1ZMR_prod_MSMmapping.csv"
meta_dist_file   = f"{TESTDIR}/MSM/1ZMR_prod_meta_dist.npy"
LiPMS_exp_file   = f"{DATASTORE}/user_input/experimental_data/ecPGK_significant_LiPMS_peptide_R1_merged.xlsx"
sasa_data_file   = f"{TESTDIR}/compare_sim2exp/SASA.npy"
XLMS_exp_file    = f"{DATASTORE}/user_input/experimental_data/ecPGK_significant_XLMS_peptide_R1_merged.xlsx"
dist_data_file   = f"{TESTDIR}/compare_sim2exp/Jwalk.npy"
cluster_data_file = f"{TESTDIR}/nonnative_entanglement_clustering/cluster_data_topoly_linking_number.npz"
OPpath           = f"{TESTDIR}/OP_last67_density1/"
AAdcd_dir        = AA_TRAJ_DIR
native_AA_pdb    = f"{DATASTORE}/user_input/reference_structures/1zmr_model_clean.pdb"
outdir           = f"{OUTDIR}/MassSpec_ConsistencyTest"

# Protocol-specific parameters for the 1ZMR / ecPGK example system
state_idx_list   = [4, 6, 8]   # metastable state indices to test
prot_len         = 387          # protein length in residues
last_num_frames  = 335          # number of frames per trajectory analyzed
rm_traj_list     = []           # trajectories excluded in Workflow 2 Step 10
native_state_idx = 9            # index of the native/folded metastable state
start            = 6600         # first frame index used in OP calculations
ID               = "1ZMR"

MS = MassSpec(
    msm_data_file=msm_data_file,
    meta_dist_file=meta_dist_file,
    LiPMS_exp_file=LiPMS_exp_file,
    sasa_data_file=sasa_data_file,
    XLMS_exp_file=XLMS_exp_file,
    dist_data_file=dist_data_file,
    cluster_data_file=cluster_data_file,
    OPpath=OPpath,
    AAdcd_dir=AAdcd_dir,
    native_AA_pdb=native_AA_pdb,
    state_idx_list=state_idx_list,
    prot_len=prot_len,
    last_num_frames=last_num_frames,
    rm_traj_list=rm_traj_list,
    native_state_idx=native_state_idx,
    outdir=outdir,
    ID=ID,
    start=start
)
```

### 17b. Run the consistency test

```python
consist_data_file, consist_result_file = MS.LiP_XL_MS_ConsistencyTest()
print("Consistency data:", consist_data_file)
print("Consistency results:", consist_result_file)
```

### 17c. Select representative structures

```python
MS.select_rep_structs(
    consist_data_file,
    consist_result_file,
    total_traj_num_frames=335,   # total frames per trajectory in the full run
    last_num_frames=67            # frames from which representative structs are selected
)
```

### What the consistency test does

For each metastable state in `state_idx_list`, the test evaluates whether:

- LiP-MS signals (elevated protease cleavage) are concentrated in residues that are **more solvent-exposed** in that state than in the native state.
- XL-MS signals (cross-links formed / lost) are consistent with the **inter-residue distances** in that state.

A p-value is computed for each state/experiment combination.

### Expected outputs

| File | Contents |
|------|----------|
| `LiPMS_XLMS_consist_pvalues_metastates_*.xlsx` | per-state p-values for LiP-MS and XL-MS consistency |
| `consist_signal_struct_data.npz` | raw per-residue consistency arrays |
| `Consistent_structures_v8.xlsx` | selected representative structures from consistent states |

Compare to pre-computed reference files in `$TESTDIR/compare_sim2exp/`.

---

## Step 18. Visualize representative structures

Using the structure identifiers in `Consistent_structures_v8.xlsx`, load the corresponding all-atom trajectory frames into VMD:

```tcl
# In the VMD TkConsole — load the all-atom trajectory for the representative trajectory
set traj_num 420
mol new /scratch/ims86/EntDetect_Datastore/user_input/reference_structures/1zmr_model_clean.pdb
mol addfile /scratch/ims86/EntDetect_Datastore/user_input/aa_trajectories/420_prod_aa.dcd
# Navigate to the frame index from the results file
animate goto <frame_index>
```

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `FileNotFoundError` for SASA or Jwalk npy | Path to `$TESTDIR/compare_sim2exp/` incorrect | Verify `TESTDIR` and run `ls compare_sim2exp/` |
| `MassSpec` raises shape mismatch | `prot_len`, `last_num_frames`, or `start` inconsistent with pre-computed data | Confirm values match those used during OP calculation |
| Consistency test returns no significant states | `state_idx_list` out of range | Check valid state indices in `1ZMR_prod_MSMmapping.csv` |
| Back-mapping produces clashing structures | Pulchra or PD2 reconstruction failure | Inspect CG input frame; remove frames with very distorted Cα geometry |

---

← [Workflow 2](workflow2_trajectory_analysis.md) | [Back to Master Index](index.md) | Next → [Workflow 4: Population Analysis](workflow4_population.md)
