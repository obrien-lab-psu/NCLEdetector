# Workflow 3: Compare Structural Ensembles to High-Throughput Experimental Data

← [Back to Master Index](index.md)

---

## Goal

Test whether specific metastable states from simulation are statistically consistent with experimental conformational signals (LiP-MS and XL-MS). Identify and extract representative structures from the best-supported metastable states.

---

## Table of Contents

- [Step 0. Activate your environment and set paths](#step-0-activate-your-environment-and-set-paths)
- [Step 1. Test metastable state consistency with experimental signals](#step-1-test-metastable-state-consistency-with-experimental-signals)
- [Step 2. Select representative structures](#step-2-select-representative-structures)
- [About experimental input files](#about-experimental-input-files)
- [Minimal Workflow – 7: Consistency test between experimental signals and changes in NCLE status](#minimal-workflow--7-consistency-test-between-experimental-signals-and-changes-in-ncle-status)

---

## Typical runtime

| Step | Runtime |
|------|---------|
| LiP-MS / XL-MS consistency test | 1–3 hours |
| Structure selection | Minutes |

> **Note:** Per-trajectory SASA and XP (Jwalk) calculations are computed in **Workflow 2** and stored as individual files. When you initialize `MassSpec`, it will automatically collect those per-trajectory files into consolidated `.npy` arrays if they don't already exist (cached), then run the consistency test.

> **Simplified workflow:** Back-mapping to all-atom structures is no longer an explicit step. All-atom trajectories (`{N}_prod_aa.dcd`) are assumed to already exist in the datastore from Workflow 2. If you need to generate them, use a separate MDAnalysis or VMD workflow outside this pipeline.

---

## Pre-computed inputs from Workflow 2

All upstream inputs from Workflow 2 are available in the datastore:

```
$DATASTORE/outputs/workflow2/
├── MSM/
│   ├── 1ZMR_prod_MSMmapping.csv
│   └── 1ZMR_prod_meta_dist.npy
├── OP_AA/SASA/
│   ├── 1ZMR_Traj1.SASA
│   ├── 1ZMR_Traj2.SASA
│   └── ... (one per trajectory)
├── OP_AA/XP/
│   ├── 1ZMR_Traj1.XP
│   ├── 1ZMR_Traj2.XP
│   └── ... (one per trajectory)
└── nonnative_clustering/
    └── cluster_data_topoly_linking_number.npz
```

This workflow writes its outputs to:

```
$DATASTORE/outputs/workflow3/
├── SASA.npy                                          # auto-collected when MassSpec is initialized
├── Jwalk.npy                                         # optional; auto-collected if requested
└── MassSpec_ConsistencyTest/
    ├── LiPMS_XLMS_consist_pvalues_metastates_*.xlsx
    ├── LiPMS_XLMS_consist_data_*.npz
    ├── Consistent_structures_*.xlsx
    ├── viz_rep_struct/                               # extracted all-atom structure snapshots
    └── 1ZMR.log
```

Where:
```bash
DATASTORE=/scratch/ims86/EntDetect_Datastore
OUTDIR=$DATASTORE/outputs/workflow3
```

---

## Required input files (from Workflow 2)

| File | Path | Description |
|------|------|-------------|
| Per-traj SASA files | `$DATASTORE/outputs/workflow2/OP_AA/SASA/{ID}_Traj{N}.SASA` | Solvent-accessible surface area (one per trajectory) |
| Per-traj XP files | `$DATASTORE/outputs/workflow2/OP_AA/XP/{ID}_Traj{N}.XP` | Jwalk cross-link predictions (one per trajectory) |
| MSM mapping | `$DATASTORE/outputs/workflow2/MSM/1ZMR_prod_MSMmapping.csv` | Per-frame metastable state assignments |
| MSM meta distribution | `$DATASTORE/outputs/workflow2/MSM/1ZMR_prod_meta_dist.npy` | Distance matrix between metastable states |
| Non-native clustering data | `$DATASTORE/outputs/workflow2/nonnative_clustering/cluster_data_topoly_linking_number.npz` | Clustering results for validation |
| Native all-atom PDB | `$DATASTORE/user_input/reference_structures/1zmr_model_clean.pdb` | Native reference structure |
| AA trajectory DCDs | `$DATASTORE/user_input/aa_trajectories/{N}_prod_aa.dcd` | All-atom trajectories (N=1–1000) |

**Experimental inputs:**

| File | Path |
|------|------|
| LiP-MS data | `$DATASTORE/user_input/experimental_data/ecPGK_significant_LiPMS_peptide_R1_merged.xlsx` |
| XL-MS data | `$DATASTORE/user_input/experimental_data/ecPGK_significant_XLMS_peptide_R1_merged.xlsx` |

---

## Step 0. Activate your environment and set paths

```bash
source ~/.bashrc
conda activate entdetect

DATASTORE=/scratch/ims86/EntDetect_Datastore
OUTDIR=$DATASTORE/outputs/workflow3

mkdir -p $OUTDIR/MassSpec_ConsistencyTest $OUTDIR/logs
```

---

## Step 1. Test metastable state consistency with experimental signals 

This is the key integrative step. It tests whether specific metastable states identified in Workflow 2 are statistically consistent with the experimental signals.

When you initialize `MassSpec`, it will automatically:

1. Check if `SASA.npy` already exists in the output directory (cached from a previous run)
2. If **not** found, run `CollectOP` to consolidate per-trajectory SASA and XP files into arrays
3. Run the consistency test

### Initialize `MassSpec` and run the consistency test

```python
from EntDetect.compare_sim2exp import MassSpec

# ── Paths ──────────────────────────────────────────────────────────────────
DATASTORE   = "/scratch/ims86/EntDetect_Datastore"
OUTDIR      = f"{DATASTORE}/outputs/workflow3"
sasa_dir    = f"{DATASTORE}/outputs/workflow2/OP_AA/SASA"  # per-traj SASA files
xp_dir      = f"{DATASTORE}/outputs/workflow2/OP_AA/XP"    # per-traj XP files
OPpath      = f"{DATASTORE}/outputs/workflow2/OP_last67/"
AAdcd_dir   = f"{DATASTORE}/user_input/aa_trajectories"
TestOutdir  = f"{OUTDIR}/MassSpec_ConsistencyTest"

# ── Inputs ──────────────────────────────────────────────────────────────────
msm_data_file    = f"{DATASTORE}/outputs/workflow2/MSM/1ZMR_prod_MSMmapping.csv"
meta_dist_file   = f"{DATASTORE}/outputs/workflow2/MSM/1ZMR_prod_meta_dist.npy"
LiPMS_exp_file   = f"{DATASTORE}/user_input/experimental_data/ecPGK_significant_LiPMS_peptide_R1_merged.xlsx"
XLMS_exp_file    = f"{DATASTORE}/user_input/experimental_data/ecPGK_significant_XLMS_peptide_R1_merged.xlsx"
cluster_data_file = f"{DATASTORE}/outputs/workflow2/nonnative_clustering/cluster_data_topoly_linking_number.npz"
native_AA_pdb    = f"{DATASTORE}/user_input/reference_structures/1zmr_model_clean.pdb"

# Protocol-specific parameters for the 1ZMR / ecPGK example system
state_idx_list   = [4, 6, 8]   # metastable state indices to test
prot_len         = 387          # protein length in residues
last_num_frames  = 67          # number of frames per trajectory analyzed
n_traj           = 1000        # total number of trajectories
n_frames         = 67          # frames per trajectory in Workflow 2 output
rm_traj_list   = [65, 75, 155, 162, 199, 231, 264, 286, 296, 314, 354, 417, 448, 472, 473, 474, 577, 579, 591, 703, 704, 732, 758, 812, 833, 870, 876, 944, 967]
native_state_idx = 9            # index of the native/folded metastable state
start            = 6600         # first frame index used in OP calculations
ID               = "1ZMR"

# ── Run  ───────────────────────────────────────────────────────────────
# MassSpec will automatically collect SASA/XP from sasa_dir if not cached
MS = MassSpec(
    msm_data_file=msm_data_file,
    meta_dist_file=meta_dist_file,
    LiPMS_exp_file=LiPMS_exp_file,
    XLMS_exp_file=XLMS_exp_file,
    cluster_data_file=cluster_data_file,
    OPpath=OPpath,
    AAdcd_dir=AAdcd_dir,
    native_AA_pdb=native_AA_pdb,
    sasa_dir=sasa_dir,           # MassSpec will handle collection internally
    xp_dir=xp_dir,               # for streaming or collection
    n_traj=n_traj,
    n_frames=n_frames,
    state_idx_list=state_idx_list,
    prot_len=prot_len,
    last_num_frames=last_num_frames,
    rm_traj_list=rm_traj_list,
    native_state_idx=native_state_idx,
    outdir=TestOutdir,
    ID=ID,
    start=start,
    num_perm=1000,
    n_boot=100,
    lag_frame=20,
    nproc=8,
)

# Run the consistency test
consist_data_file, consist_result_file = MS.LiP_XL_MS_ConsistencyTest()
```

### What the consistency test does

For each metastable state in `state_idx_list`, the test evaluates whether:

- **LiP-MS signals** (elevated protease cleavage) are concentrated in residues that are **more solvent-exposed** in that state than in the native state.
- **XL-MS signals** (cross-links formed / lost) are consistent with the **inter-residue distances** in that state.

A p-value and effect size are computed for each state/experiment combination via permutation testing.

### Expected outputs

| File | Contents |
|------|----------|
| `LiPMS_XLMS_consist_pvalues_metastates_*.xlsx` | Per-state p-values, effect sizes, and significance flags for LiP-MS and XL-MS consistency |
| `LiPMS_XLMS_consist_data_*.npz` | Raw per-residue and per-pair consistency arrays for further analysis |
| `1ZMR.log` | Detailed execution log with state-by-state statistics |

---

## Step 2. Select representative structures

Using the consistency test results, select representative all-atom structures from the best-supported metastable states.

```python
# Select representative structures from states with significant signals
MS.select_rep_structs(
    consist_data_file,
    consist_result_file,
    total_traj_num_frames=335,   # total frames per trajectory in the full run
    last_num_frames=67            # frames from which representatives are selected
)
```

### What this produces

| File / Directory | Contents |
|---|---|
| `Consistent_structures_*.xlsx` | Table of (trajectory, frame) pairs selected as representative from each consistent state |
| `viz_rep_struct/` | Directory containing extracted all-atom structure snapshots for visualization (PDB format) |

### Visualize selected structures

Using the structure identifiers in the output Excel file, load the corresponding all-atom trajectory frames into VMD:

```tcl
# In VMD TkConsole — load the all-atom trajectory
mol new 1ZMR_Traj42_prod_aa.dcd type dcd
mol addfile 1zmr_model_clean.pdb type pdb

# Go to a specific frame (e.g., frame 50)
animate goto 50

# Color by B-factor to highlight consistency signals
color Name CA white
color Name {"LiP"} red    ;# highlight LiP-MS consistent residues
color Name {"XL"}  blue   ;# highlight XL-MS consistent residues
```

---

## About experimental input files

The processed experimental files are already available in `$DATASTORE/user_input/experimental_data/`. 

### To prepare your own experimental inputs:

1. **Run your LiP-MS statistical analysis** to identify significantly changed peptides.
   - Define significance (e.g., q-value < 0.05, fold-change > 2×).
   - Map peptides back to residue-level assignments.
   
2. **Prepare an Excel file** with columns:
   - `Peptide`: identifier
   - `Start_res`, `End_res`: residue range
   - `LFC` (log fold-change) or similar quantitative signal

3. **Repeat for XL-MS data** with cross-link pairs instead of peptides.

> **Critical:** The choice of significance threshold and mapping strategy can strongly affect downstream conclusions. Be consistent and document your decisions.

---

## Minimal Workflow – 7: Consistency test between experimental signals and changes in NCLE status 

All steps above are handled by `scripts/run_compare_sim2exp.py`. Run it from the repo root.

You can provide parameters either directly as CLI flags, via a `--config` file, or both.
When both are provided, **CLI flags override config values** for the same parameter.

### Example with direct CLI flags:

```bash
source ~/.bashrc
conda activate entdetect

DATASTORE=/scratch/ims86/EntDetect_Datastore

python scripts/run_compare_sim2exp.py \
    --msm_data_file   $DATASTORE/outputs/workflow2/MSM/1ZMR_prod_MSMmapping.csv \
    --meta_dist_file  $DATASTORE/outputs/workflow2/MSM/1ZMR_prod_meta_dist.npy \
    --LiPMS_exp_file  $DATASTORE/user_input/experimental_data/ecPGK_significant_LiPMS_peptide_R1_merged.xlsx \
    --XLMS_exp_file   $DATASTORE/user_input/experimental_data/ecPGK_significant_XLMS_peptide_R1_merged.xlsx \
    --cluster_data_file $DATASTORE/outputs/workflow2/nonnative_clustering/cluster_data_topoly_linking_number.npz \
    --OPpath          $DATASTORE/outputs/workflow2/OP_last67/ \
    --AAdcd_dir       $DATASTORE/user_input/aa_trajectories/ \
    --native_AA_pdb   $DATASTORE/user_input/reference_structures/1zmr_model_clean.pdb \
    --sasa_dir        $DATASTORE/outputs/workflow2/OP_AA/SASA \
    --xp_dir          $DATASTORE/outputs/workflow2/OP_AA/XP \
    --state_idx_list  4 6 8 \
    --prot_len        387 \
    --last_num_frames 67 \
    --rm_traj_list    65 75 155 162 199 231 264 286 296 314 354 417 448 472 473 474 577 579 591 703 704 732 758 812 833 870 876 944 967 \
    --native_state_idx 9 \
    --outdir          $DATASTORE/outputs/workflow3/MassSpec_ConsistencyTest \
    --ID              1ZMR \
    --start           6600 \
    --end             -1 \
    --stride          1 \
    --n_traj          1000 \
    --n_frames        67 \
    --num_perm        1000 \
    --n_boot          100 \
    --lag_frame       20 \
    --nproc           8
```

### Example with config plus CLI override:

```bash
CONFIG=scripts/configs/workflow3_consistency_config.json

# Here, --state_idx_list 5 7 overrides state_idx_list from the config file.
python scripts/run_compare_sim2exp.py \
    --config $CONFIG \
    --state_idx_list 5 7
```

### Config file example (matches `scripts/configs/workflow3_consistency_config.json`):

```json
{
    "msm_data_file": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/MSM/1ZMR_prod_MSMmapping.csv",
    "meta_dist_file": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/MSM/1ZMR_prod_meta_dist.npy",
    "LiPMS_exp_file": "/scratch/ims86/EntDetect_Datastore/user_input/experimental_data/ecPGK_significant_LiPMS_peptide_R1_merged.xlsx",
    "XLMS_exp_file": "/scratch/ims86/EntDetect_Datastore/user_input/experimental_data/ecPGK_significant_XLMS_peptide_R1_merged.xlsx",
    "cluster_data_file": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/nonnative_clustering/cluster_data_topoly_linking_number.npz",
    "OPpath": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/OP_last67/",
    "AAdcd_dir": "/scratch/ims86/EntDetect_Datastore/user_input/aa_trajectories/",
    "native_AA_pdb": "/scratch/ims86/EntDetect_Datastore/user_input/reference_structures/1zmr_model_clean.pdb",
    "sasa_dir": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/OP_AA/SASA",
    "xp_dir": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/OP_AA/XP",
    "outdir": "/scratch/ims86/EntDetect_Datastore/outputs/workflow3/MassSpec_ConsistencyTest",
    "ID": "1ZMR",
    "state_idx_list": [4, 6, 8],
    "prot_len": 387,
    "last_num_frames": 67,
    "native_state_idx": 9,
    "rm_traj_list": [65, 75, 155, 162, 199, 231, 264, 286, 296, 314, 354, 417, 448, 472, 473, 474, 577, 579, 591, 703, 704, 732, 758, 812, 833, 870, 876, 944, 967],
    "start": 6600,
    "end": -1,
    "stride": 1,
    "n_traj": 1000,
    "n_frames": 67,
    "num_perm": 1000,
    "n_boot": 100,
    "lag_frame": 20,
    "nproc": 8,
    "verbose": false
}
```

The JSON/YAML config keys and their matching CLI flags are listed below. For the command-line wrapper, use the same key name with a `--` prefix; the only extra wrapper-only flag is `--config`.

| Term | Definition |
|------|------------|
| `msm_data_file` (`--msm_data_file`) | Input MSM mapping CSV from Workflow 2 containing the per-frame microstate and metastable-state assignments used to align simulation frames with experimental comparisons. |
| `meta_dist_file` (`--meta_dist_file`) | NumPy array of MSM-derived distances between metastable states, used when ranking and selecting representative structures after the consistency test. |
| `LiPMS_exp_file` (`--LiPMS_exp_file`) | Processed LiP-MS workbook containing the experimentally significant proteolysis signals to compare against simulated solvent exposure. |
| `XLMS_exp_file` (`--XLMS_exp_file`) | Processed XL-MS workbook containing the experimentally significant cross-link signals to compare against simulated inter-residue distances. |
| `cluster_data_file` (`--cluster_data_file`) | Non-native entanglement clustering archive used to link experiment-consistent structures back to specific entanglement-change classes. |
| `OPpath` (`--OPpath`) | Directory of Workflow 2 order-parameter outputs used when annotating and ranking representative structures by Q and G. |
| `AAdcd_dir` (`--AAdcd_dir`) | Directory or glob for all-atom trajectory DCD files used to extract the representative structure snapshots. |
| `native_AA_pdb` (`--native_AA_pdb`) | Native all-atom PDB used as the structural reference and topology for representative-structure extraction. |
| `sasa_dir` (`--sasa_dir`) | Directory of per-trajectory `.SASA` files used in collection mode to supply the SASA data needed for LiP-MS scoring. |
| `xp_dir` (`--xp_dir`) | Directory of per-trajectory `.XP` files used to stream XL-MS-compatible distance data during the consistency test and representative-structure reporting. |
| `outdir` (`--outdir`) | Output directory for the consistency-test workbooks, NumPy archives, representative-structure summary workbook, and default `ID.log` file. |
| `ID` (`--ID`) | System identifier used in output filenames and as the default log-file stem. |
| `state_idx_list` (`--state_idx_list`) | List of metastable-state indices to test against the experimental data. |
| `prot_len` (`--prot_len`) | Protein length in residues, used when building residue-index mappings and interpreting LiP-MS or XL-MS signals. |
| `last_num_frames` (`--last_num_frames`) | Number of final frames retained from each trajectory for the experiment-versus-simulation comparison. |
| `native_state_idx` (`--native_state_idx`) | Metastable-state index treated as the native or near-native reference state in the consistency statistics. |
| `rm_traj_list` (`--rm_traj_list`) | Trajectory IDs to exclude before scoring, typically mirror-image or otherwise invalid trajectories that should not contribute to the statistics. |
| `start` (`--start`) | First frame index associated with the analyzed order-parameter window in the original trajectories. |
| `end` (`--end`) | Last frame index associated with the analyzed order-parameter window; `-1` means use the full remaining trajectory range. |
| `stride` (`--stride`) | Frame stride used when interpreting the selected trajectory window. |
| `n_traj` (`--n_traj`) | Total number of trajectories expected in collection-mode inputs. In the shipped example this matches the number of Workflow 2 trajectories. |
| `n_frames` (`--n_frames`) | Number of frames per trajectory expected in the per-trajectory collection-mode files. |
| `num_perm` (`--num_perm`) | Number of permutation samples used for the LiP-MS and XL-MS hypothesis tests. |
| `n_boot` (`--n_boot`) | Number of bootstrap samples used when estimating confidence intervals and resampled consistency statistics. |
| `lag_frame` (`--lag_frame`) | Down-sampling interval in frames used by the consistency test when building the workbook outputs. |
| `nproc` (`--nproc`) | Number of worker processes used for the parallel parts of the consistency analysis. |
| `verbose` (`--verbose`) | Enables verbose execution behavior for the analysis object. |

If you want to use pre-collected arrays instead of the per-trajectory directories, the wrapper also accepts `sasa_data_file` and `dist_data_file` as alternate config or CLI keys.

## I/O Reference for run_compare_sim2exp.py

Below, `$OUTDIR` refers to the `outdir` value passed to `run_compare_sim2exp.py`. In the shipped example config, that is `$DATASTORE/outputs/workflow3/MassSpec_ConsistencyTest`.

Each file below is listed once, followed by the column-level or archive-level description table for structured outputs.

### `$DATASTORE/outputs/workflow2/MSM/1ZMR_prod_MSMmapping.csv`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/outputs/workflow2/MSM/1ZMR_prod_MSMmapping.csv` | Per-frame MSM assignment table used to identify which trajectories and frames belong to each metastable state tested against experiment. |

### `$DATASTORE/outputs/workflow2/MSM/1ZMR_prod_meta_dist.npy`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/outputs/workflow2/MSM/1ZMR_prod_meta_dist.npy` | NumPy array containing the MSM-derived distance or transition metric used to score representative structures within each metastable state. |

### `$DATASTORE/user_input/experimental_data/ecPGK_significant_LiPMS_peptide_R1_merged.xlsx`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/user_input/experimental_data/ecPGK_significant_LiPMS_peptide_R1_merged.xlsx` | Processed LiP-MS signal table defining the peptide or residue-level experimental perturbations to test against simulated solvent exposure. |

### `$DATASTORE/user_input/experimental_data/ecPGK_significant_XLMS_peptide_R1_merged.xlsx`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/user_input/experimental_data/ecPGK_significant_XLMS_peptide_R1_merged.xlsx` | Processed XL-MS signal table defining the experimental cross-links or cross-link changes to test against simulated distances. |

### `$DATASTORE/outputs/workflow2/nonnative_clustering/cluster_data_topoly_linking_number.npz`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/outputs/workflow2/nonnative_clustering/cluster_data_topoly_linking_number.npz` | Compressed clustering data describing non-native entanglement-change classes used when linking experimental consistency back to structural mechanisms. |

### `$DATASTORE/outputs/workflow2/OP_last67/`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/outputs/workflow2/OP_last67/` | Directory containing the Q and G order-parameter files for the last analyzed frames of each trajectory, used when reporting representative state structures. |

### `$DATASTORE/user_input/aa_trajectories/<TRAJ>_prod_aa.dcd`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/user_input/aa_trajectories/<TRAJ>_prod_aa.dcd` | All-atom trajectory file used to recover the structural snapshots selected as representative of experiment-consistent metastable states. |

### `$DATASTORE/user_input/reference_structures/1zmr_model_clean.pdb`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/user_input/reference_structures/1zmr_model_clean.pdb` | Native all-atom structure used as the structural reference for SASA or XP comparisons and for extracting representative coordinates. |

### `$DATASTORE/outputs/workflow2/OP_AA/SASA/`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/outputs/workflow2/OP_AA/SASA/` | Directory of per-trajectory SASA files used in collection mode to assemble the cached SASA array consumed by the consistency analysis. |

### `$DATASTORE/outputs/workflow2/OP_AA/XP/`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/outputs/workflow2/OP_AA/XP/` | Directory of per-trajectory XP files used either to stream XL-MS distance data directly or to build the optional cached Jwalk array. |

### `$OUTDIR/SASA.npy`

| I/O | File | File Description |
|---|---|---|
| Input or Output | `$OUTDIR/SASA.npy` | Cached three-dimensional NumPy array of SASA values with axes `(trajectory, frame, residue)`. In collection-mode runs this cache is created or reused under `outdir`; in direct-array mode the same file can be supplied explicitly via `sasa_data_file`. |

### `$OUTDIR/Jwalk.npy`

| I/O | File | File Description |
|---|---|---|
| Optional Input | `$OUTDIR/Jwalk.npy` | Optional cached three-dimensional NumPy array of XL-MS-compatible distance data with axes `(trajectory, frame, signal)`. When provided via `dist_data_file`, it replaces the need to stream per-trajectory `.XP` files during XL-MS scoring and representative-structure annotation. |

### `$OUTDIR/LiPMS_XLMS_consist_pvalues_metastates_v11_down_sample_lag<LAG>.xlsx`

| I/O | File | File Description |
|---|---|---|
| Output | `$OUTDIR/LiPMS_XLMS_consist_pvalues_metastates_v11_down_sample_lag<LAG>.xlsx` | Multi-sheet Excel workbook summarizing LiP-MS and XL-MS consistency statistics for each tested metastable state. Sheets include `LiPMS`, `XLMS`, and `All states`, with repeated signal blocks separated by blank rows. |

| Workbook Element | Description |
|---|---|
| `LiPMS` sheet | Block-formatted LiP-MS significance results for each experimental cut site and each tested metastable state. |
| `XLMS` sheet | Block-formatted XL-MS significance results for each experimental cross-link signal and each tested metastable state. |
| `All states` sheet | Combined state-level summary table used for quick comparison across the tested metastable states. |

### Columns in `LiPMS` workbook blocks

| Column Name | Column Description |
|---|---|
| index (unlabeled) | LiP-MS signal identifier written as the Excel row index, typically a residue or peptide label such as `D8`. |
| Near-native state | Metastable-state index being compared to the native or near-native reference state. |
| Sample size | Number of sampled state frames compared against the native-state reference pool, reported as `state_count vs. reference_count`. |
| `<M>` | Summary statistic for the LiP-MS consistency metric in the tested state versus the native reference, reported together with its confidence bounds. |
| `p (!=)` | Two-sided p-value testing whether the state differs from the reference. |
| `Adjusted p (!=)` | Multiple-testing-adjusted version of the two-sided p-value. |
| `p (>)` | One-sided p-value testing whether the state shows increased LiP-MS-like signal relative to the reference. |
| `Adjusted p (>)` | Multiple-testing-adjusted version of the one-sided enrichment p-value. |

### Columns in `XLMS` workbook blocks

| Column Name | Column Description |
|---|---|
| index (unlabeled) | XL-MS signal identifier written as the Excel row index, typically a residue-pair label such as `S2-K179`. |
| Near-native state | Metastable-state index being compared to the native or near-native reference state. |
| Sample size | Number of sampled state frames compared against the native-state reference pool, reported as `state_count vs. reference_count`. |
| `<M>` | Summary statistic for the XL-MS consistency metric in the tested state versus the native reference, reported together with its confidence bounds. |
| `p (!=)` | Two-sided p-value testing whether the state differs from the reference. |
| `Adjusted p (!=)` | Multiple-testing-adjusted version of the two-sided p-value. |
| `p (<)` | One-sided p-value testing whether the simulated state becomes more cross-link-consistent than the reference. |
| `Adjusted p (<)` | Multiple-testing-adjusted version of the one-sided depletion or distance-shortening p-value. |

### `$OUTDIR/LiPMS_XLMS_consist_data_v9.npz`

| I/O | File | File Description |
|---|---|---|
| Output | `$OUTDIR/LiPMS_XLMS_consist_data_v9.npz` | Compressed NumPy archive storing the raw consistency tensors that back the Excel summary workbook. |

| Stored Array | Array Description |
|---|---|
| `M_LiPMS` | Three-dimensional array with axes trajectory, frame, and LiP-MS signal index, containing the state-aligned LiP-MS consistency values. |
| `M_XLMS` | Three-dimensional array with axes trajectory, frame, and XL-MS signal index, containing the state-aligned XL-MS consistency values. |

### `$OUTDIR/native_M_data.npz`

| I/O | File | File Description |
|---|---|---|
| Output | `$OUTDIR/native_M_data.npz` | Compressed NumPy archive containing the native-state LiP-MS and XL-MS consistency distributions used to define the reference confidence bounds in the downstream summaries. |

| Stored Array | Array Description |
|---|---|
| `M_LiPMS` | Native-reference LiP-MS consistency values with axes `(sample, LiP-MS signal index)`. |
| `M_XLMS` | Native-reference XL-MS consistency values with axes `(sample, XL-MS signal index)`. |

### `$OUTDIR/consist_signal_struct_data.npz`

| I/O | File | File Description |
|---|---|---|
| Output | `$OUTDIR/consist_signal_struct_data.npz` | Compressed NumPy archive storing the grouped signal-to-structure assignments used to build the representative-structure workbook and visualization output. |

| Stored Object | Object Description |
|---|---|
| `last_num_frames` | Number of final frames per trajectory retained for representative-structure selection. |
| `total_traj_num_frames` | Total frame count per trajectory in the source simulation used to map local frame indices back to absolute trajectory frames. |
| `LIPMS_consist_data` | LiP-MS signal-level consistency summary objects used during grouping and ranking. |
| `XLMS_consist_data` | XL-MS signal-level consistency summary objects used during grouping and ranking. |
| `LIPMS_struct_data` | Mapping from LiP-MS signals to the trajectory/frame subsets that satisfy the consistency criteria. |
| `XLMS_struct_data` | Mapping from XL-MS signals to the trajectory/frame subsets that satisfy the consistency criteria. |
| `dtrajs_MS` | Object array assigning each retained trajectory frame to the set of consistency-supported signals present there. |
| `sorted_consist_signal_dict` | Grouped and sorted dictionary summarizing recurring combinations of consistency-supported experimental signals. |
| `group_dict` | Grouped structure assignments for each state, consistency-signal set, and entanglement-change class. |
| `rep_group_dict` | Selected representative trajectory/frame indices for each state, signal group, and entanglement-change class. |

### `$OUTDIR/Consistent_structures_v8.xlsx`

| I/O | File | File Description |
|---|---|---|
| Output | `$OUTDIR/Consistent_structures_v8.xlsx` | Multi-sheet Excel workbook listing the most experiment-consistent structural representatives overall and for each tested metastable state. Sheets include `Total` plus one sheet per state such as `State 5`. |

| Workbook Element | Description |
|---|---|
| `Total` sheet | Cross-state summary of recurrent signal-supported structural subgroups. |
| `State <IDX>` sheets | Per-state representative-structure selections and their associated signal and entanglement annotations. |

### Columns in `Total` sheet

| Column Name | Column Description |
|---|---|
| Consistent signals | Comma-separated list of LiP-MS and XL-MS signals jointly supported by a structural subgroup. |
| IDs of Changes in Entanglements | Comma-separated identifiers of the clustered non-native entanglement changes associated with that subgroup. |
| Number of consistent signals | Count of distinct experimental signals represented in the subgroup. |
| Number of Structures | Number of trajectory frames assigned to that subgroup. |

### Columns in `State <IDX>` sheets

| Column Name | Column Description |
|---|---|
| Consistent signals | Comma-separated list of LiP-MS and XL-MS signals jointly supported by the representative structure. |
| IDs of Changes in Entanglements | Comma-separated identifiers of the clustered non-native entanglement changes present in that representative structure. |
| Number of consistent signals | Count of distinct experimental signals represented in the structure-level subgroup. |
| Number of Structures | Number of frames in the subgroup represented by the selected structure. |
| `Representative Structure (Traj #, Frame #)` | Selected all-atom trajectory and frame index used as the representative structure for the subgroup. |
| Prob | MSM-derived probability or weighting term used to prioritize the representative within the state. |
| Q | Native-contact order parameter of the representative frame. |
| G | Entanglement-change order parameter of the representative frame. |

### `$OUTDIR/1ZMR.log`

| I/O | File | File Description |
|---|---|---|
| Output | `$OUTDIR/1ZMR.log` | Default MassSpec execution log written with the configured `ID` as the log-file stem. It records loading, caching, consistency-test progress, and representative-structure extraction details. |

### `./viz_rep_struct/`

| I/O | File | File Description |
|---|---|---|
| Output | `./viz_rep_struct/` | Directory tree containing extracted representative structure files organized by state for downstream visualization in VMD or PyMOL. The current implementation creates this folder relative to the working directory from which `run_compare_sim2exp.py` is launched, rather than inside `outdir`. |

### Or via SLURM:

```bash
sbatch assets/slurm/scripts/run_workflow3_consistency.slurm
```

---

[**Next:** Workflow 4 — Population Modeling (optional)](workflow4_population.md)
