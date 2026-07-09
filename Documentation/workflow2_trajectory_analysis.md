# Workflow 2: Detect Changes in Entanglement Across Trajectories

← [Back to Master Index](index.md)

---

## Goal

Compute order parameters (Q, G, K, SASA, XP) across simulation trajectories, cluster non-native entanglement changes, build a Markov state model (MSM), and quantify differences in metastable-state behavior.

---

## Table of Contents

- [Step 0. Activate your environment and set paths](#step-0-activate-your-environment-and-set-paths)
- [Step 1. Compute order parameters Q, G, K on the CG trajectory](#step-1-compute-order-parameters-q-g-k-on-the-cg-trajectory)
- [Step 2. Compute order parameters SASA and XP on the all-atom trajectory](#step-2-compute-order-parameters-sasa-and-xp-on-the-all-atom-trajectory)
- [Minimal Workflow – 2: Running the OP calculations as single script](#minimal-workflow--2-running-the-op-calculations-as-single-script)
- [Step 3. Identify and remove artificial mirror conformations](#step-3-identify-and-remove-artificial-mirror-conformations)
- [Step 4. Cluster changes of NCLE status to remove redundancy](#step-4-cluster-changes-of-ncle-status-to-remove-redundancy)
- [Minimal Workflow – 3: Running the change in NCLE status clustering as single script](#minimal-workflow--3-running-the-change-in-ncle-status-clustering-as-single-script)
- [Step 5. Build a Markov state model of the GvQ probability surface](#step-5-build-a-markov-state-model-of-the-gvq-probability-surface)
- [Minimal Workflow – 4: Markov state modeling](#minimal-workflow--4-markov-state-modeling)
- [Step 5b. Label MSM Data and Define Analysis Cases](#step-5b-label-msm-data-and-define-analysis-cases)
- [Step 6. Visualize state distribution, state probability evolution, representative state structures and folding pathways (OPTIONAL)](#step-6-visualize-state-distribution-state-probability-evolution-representative-state-structures-and-folding-pathways-optional)
- [Minimal Workflow – 5: Analyzing MSM metadata state probability](#minimal-workflow--5-analyzing-msm-metadata-state-probability)
- [Step 7. Folding Pathways and Jensen-Shannon Divergence (OPTIONAL)](#step-7-folding-pathways-and-jensen-shannon-divergence-optional)
- [Minimal Workflow – 6: Folding pathway analysis and JS divergence](#minimal-workflow--6-folding-pathway-analysis-and-js-divergence)

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

Tutorial outputs from re-running the steps below are written to:

```bash
DATASTORE=/scratch/ims86/EntDetect_Datastore
OUTDIR=$DATASTORE/outputs/workflow2
```

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

Where these paths refer to the same `DATASTORE` and `OUTDIR` values shown above.

---

## Required input files

| File | Path | Notes |
|------|------|-------|
| Cα PSF topology | `$REFSTRUCT/1zmr_model_clean_ca.psf` | CG Q/G/K |
| Cα COR reference | `$REFSTRUCT/1zmr_model_clean_ca.cor` | CG Q/G/K |
| CG trajectories (all) | `$CG_TRAJ_DIR/{N}_prod.dcd` (N=1–1000) | Full production run |
| All-atom PDB topology | `$REFSTRUCT/1zmr_model_clean.pdb` | SASA/XP |
| AA trajectories (all) | `$AA_TRAJ_DIR/{N}_prod_aa.dcd` (N=1–1000) | Full AA production run |
| Secondary structure defs | `$REFSTRUCT/secondary_struc_defs.txt` | **Required for Q/G/K** |
| Domain boundary file | `$REFSTRUCT/domain_def.dat` | **Required for Q/G/K** |

```bash
REFSTRUCT=$DATASTORE/user_input/reference_structures
CG_TRAJ_DIR=$DATASTORE/user_input/cg_trajectories
AA_TRAJ_DIR=$DATASTORE/user_input/aa_trajectories
```

---

## Step 0. Activate your environment and set paths
```bash
source ~/.bashrc
conda activate entdetect
```

---

## Step 1. Compute order parameters Q, G, K on the CG trajectory

Compute the three canonical order parameters for the coarse-grained trajectory using the parameters from `assets/slurm/scripts/run_OP_traj420.slurm`.

```python
from EntDetect.order_params import CalculateOP

# ── Paths ──────────────────────────────────────────────────────────────────
DATASTORE    = "/scratch/ims86/EntDetect_Datastore"
REFSTRUCT    = f"{DATASTORE}/user_input/reference_structures"
CG_TRAJ_DIR  = f"{DATASTORE}/user_input/cg_trajectories"
AA_TRAJ_DIR  = f"{DATASTORE}/user_input/aa_trajectories"
OUTDIR       = f"{DATASTORE}/outputs/workflow2"
OP_outdir       = f"{OUTDIR}/OP"
OP_last67_outdir       = f"{OUTDIR}/OP_last67"
OP_AA_last67_outdir       = f"{OUTDIR}/OP_AA_last67"

# ── Inputs ──────────────────────────────────────────────────────────────────
Traj         = 420
PSF          = f"{REFSTRUCT}/1zmr_model_clean_ca.psf"
COR          = f"{REFSTRUCT}/1zmr_model_clean_ca.cor"
DCD          = f"{CG_TRAJ_DIR}/420_prod.dcd"
ID           = "1ZMR"
sec_elements = f"{REFSTRUCT}/secondary_struc_defs.txt"
domain       = f"{REFSTRUCT}/domain_def.dat"

pdb_file   = f"{REFSTRUCT}/1zmr_model_clean.pdb"    # all-atom topology / reference
dcd   = f"{AA_TRAJ_DIR}/420_prod_aa.dcd"

# start: first frame to include (0-indexed).
# Adjust to skip early equilibration frames for production runs.
OP_start = 0
OP_last67_start = 6600
OP_AA_last67_start = 268

# ── Initialize and Run across full traj ────────────────────────────────────────────────────
CalcOP = CalculateOP(outdir=OP_outdir, Traj=Traj, ID=ID, psf=PSF, cor=COR, sec_elements=sec_elements, dcd=DCD, domain=domain, start=OP_start, ent_detection_method=1)

Qdata_dict = CalcOP.Q()
Gdata_dict = CalcOP.G(topoly=False, Calpha=True, CG=True, nproc=10, chunk_frames=100, chunk_suffix='_chunk')
Kdata_dict = CalcOP.K()

# ── Initialize and Run across last 67 ────────────────────────────────────────────────────
CalcOP = CalculateOP(outdir=OP_last67_outdir, Traj=Traj, ID=ID, psf=PSF, cor=COR, sec_elements=sec_elements, dcd=DCD, domain=domain, start=OP_last67_start, ent_detection_method=2)

Qdata_dict = CalcOP.Q()
Gdata_dict = CalcOP.G(topoly=True, Calpha=True, CG=True, nproc=10)
Kdata_dict = CalcOP.K()
```

`Q` measures how many of the native residue–residue contacts present in the reference structure are also present in each trajectory frame. A value near 1.0 indicates a native-like conformation.

**Output:** A `.Q` file in `$OUTDIR/OP/Q/`

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

> **Runtime note:** `G` is the most expensive order parameter to compute. Expect 12–20 hours per CG trajectory at ~6700 frames on 10 cores. Submit via the cluster for a full 1000-trajectory run (see [Minimal Workflow – 2: Running the OP calculations as single script](#minimal-workflow--2-running-the-op-calculations-as-single-script)).

`K` detects frames where the protein has adopted a **mirror-image conformation** relative to the native structure. These frames are artifacts that must be removed before clustering.

**Output:** A `K_*.dat` file in `$OUTDIR/OP/K/`

> **NOTE** For each trajectory you should run the analysis twice.  
> 1. Across all frames of the trajectory using just the GLN for NCLE identification `ent_detection_method = 1` with `topoloy = False` (which is much faster than the other two methods).  
> 2. Last 10 ns (in this case 67 frames) of each trajectory using the same GLN identification method but with `topoloy = True`  
  
---

## Step 2. Compute order parameters SASA and XP on the all-atom trajectory

> **Why separate from Q, G, K?** The CG Cα-only representation does not carry enough atomic detail for accurate SASA (Shrake-Rupley requires explicit side-chain/backbone atoms) or cross-link SASD (Jwalk uses solvent-accessible surface geometry). The AA trajectories are produced by back-mapping the CG frames. Q, G, K are computed only on CG trajectories (Step 1).

```python
# ── Initialize and Run ───────────────────────────────────────────────────────────────
CalcOP = CalculateOP(
    outdir=OP_AA_last67_outdir,
    Traj=Traj,
    ID=ID,
    psf=pdb_file,
    dcd=AA_DCD,
    start=268,
)

SASAdata_dict = CalcOP.SASA()
XPdata_dict = CalcOP.XP(pdb_file=AA_PDB)
```

For the AA-only SASA/XP workflow, `CalculateOP` only needs the AA topology/PDB (`psf`), the AA trajectory (`dcd`), and the start frame. The CG-specific `cor`, `sec_elements`, and `domain` inputs are not required for Step 2.

`SASA` uses the Shrake-Rupley algorithm (mdtraj, `probe_radius=0.14 nm`, 1000 sphere points) to compute per-residue SASA for every frame. Used downstream (Workflow 3) to test LiP-MS signals: high SASA residues are more accessible to the protease.

| Key | Contents |
|-----|----------|
| `outfile` | Path to `{ID}.SASA` CSV written to `$OUTDIR/OP_AA/SASA/` |
| `result`  | DataFrame with columns `Time(ns)`, `Frame`, `resid`, `SASA(nm^2)` |


`XP` computes the **solvent-accessible surface distance (SASD)** between all pairs of cross-linkable residue types (K, S, T, Y, M) and converts each SASD to a cross-link probability score using a Gaussian parameterised on K–K linker geometry. Used downstream to test XL-MS signals.

**Output:** Per-residue-pair XP scores in `$OUTDIR/OP_AA/XP/Jwalk_results/{stem}_crosslink_list.txt`

---

## Minimal Workflow – 2: Running the OP calculations as single script

Order parameters are computed by `scripts/run_OP_on_simulation_traj.py`. Q, G, K are computed for CG trajectories; SASA and XP are computed for AA trajectories. This is the script used for the full 1000-trajectory production run.

You can provide parameters either directly as CLI flags, via a `--config` file, or both.
When both are provided, **CLI flags override config values** for the same parameter.

For direct CLI usage, pass `--CG --Calpha` for the CG `Q/G/K` runs. For the AA `SASA/XP` run, omit both flags.

```bash
# ── Activation ──────────────────────────────────────────────────────────────────
source ~/.bashrc
conda activate entdetect

# ── Config files ───────────────────────────────────────────────────────────────
CFG_OP=scripts/configs/workflow2_OP_config.json
CFG_OP_LAST67=scripts/configs/workflow2_OP_last67_config.json
CFG_OP_AA_LAST67=scripts/configs/workflow2_OP_AA_last67_config.json

# ── Run ───────────────────────────────────────────────────────────────
# 1) CG full trajectory run (Q, G, K)
python scripts/run_OP_on_simulation_traj.py \
        --config $CFG_OP

# 2) CG last-67-frame run (Q, G, K)
python scripts/run_OP_on_simulation_traj.py \
        --config $CFG_OP_LAST67

# 3) AA run (SASA, XP)
python scripts/run_OP_on_simulation_traj.py \
        --config $CFG_OP_AA_LAST67

# Example CLI override on top of config (override trajectory index only):
python scripts/run_OP_on_simulation_traj.py \
        --config $CFG_OP_LAST67 \
        --Traj 421 \
        --DCD /scratch/ims86/EntDetect_Datastore/user_input/cg_trajectories/421_prod.dcd
```

Config file example 1 (matches `scripts/configs/workflow2_OP_config.json`):

```json
{
    "Traj": "420",
    "ID": "1ZMR",
    "PSF": "/scratch/ims86/EntDetect_Datastore/user_input/reference_structures/1zmr_model_clean_ca.psf",
    "COR": "/scratch/ims86/EntDetect_Datastore/user_input/reference_structures/1zmr_model_clean_ca.cor",
    "DCD": "/scratch/ims86/EntDetect_Datastore/user_input/cg_trajectories/420_prod.dcd",
    "sec_elements": "/scratch/ims86/EntDetect_Datastore/user_input/reference_structures/secondary_struc_defs.txt",
    "domain": "/scratch/ims86/EntDetect_Datastore/user_input/reference_structures/domain_def.dat",
    "outdir": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/OP",
    "logdir": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/OP/logs",
    "start": 0,
    "ops": ["Q", "G", "K"],
    "CG": true,
    "Calpha": true,
    "ent_detection_method": 1,
    "no_topoly": true,
    "nproc": 10,
    "chunk_frames": 100,
    "chunk_suffix": "_chunk",
    "log_level": "INFO"
}
```

Config file example 2 (matches `scripts/configs/workflow2_OP_last67_config.json`):

```json
{
    "Traj": "420",
    "ID": "1ZMR",
    "PSF": "/scratch/ims86/EntDetect_Datastore/user_input/reference_structures/1zmr_model_clean_ca.psf",
    "COR": "/scratch/ims86/EntDetect_Datastore/user_input/reference_structures/1zmr_model_clean_ca.cor",
    "DCD": "/scratch/ims86/EntDetect_Datastore/user_input/cg_trajectories/420_prod.dcd",
    "sec_elements": "/scratch/ims86/EntDetect_Datastore/user_input/reference_structures/secondary_struc_defs.txt",
    "domain": "/scratch/ims86/EntDetect_Datastore/user_input/reference_structures/domain_def.dat",
    "outdir": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/OP_last67",
    "logdir": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/OP_last67/logs",
    "start": 6600,
    "ops": ["Q", "G", "K"],
    "CG": true,
    "Calpha": true,
    "ent_detection_method": 2,
    "nproc": 10,
    "log_level": "INFO"
}
```

Config file example 3 (matches `scripts/configs/workflow2_OP_AA_last67_config.json`):

```json
{
    "Traj": "420",
    "ID": "1ZMR",
    "PSF": "/scratch/ims86/EntDetect_Datastore/user_input/reference_structures/1zmr_model_clean.pdb",
    "DCD": "/scratch/ims86/EntDetect_Datastore/user_input/aa_trajectories/420_prod_aa.dcd",
    "outdir": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/OP_AA_last67",
    "logdir": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/OP_AA_last67/logs",
    "start": 268,
    "ops": ["SASA", "XP"],
    "CG": false,
    "Calpha": false,
    "pdb_file": "/scratch/ims86/EntDetect_Datastore/user_input/reference_structures/1zmr_model_clean.pdb",
    "nproc": 10,
    "log_level": "INFO"
}
```

For this AA-only config, `COR`, `sec_elements`, and `domain` are intentionally omitted because they are only required when computing `Q`, `G`, or `K`.

The JSON/YAML config keys and their matching CLI flags are listed below. For the command-line wrapper, use the same key name with a `--` prefix; the only extra wrapper-only flag is `--config`.

| Term | Definition |
|------|------------|
| `Traj` (`--Traj`) | Trajectory number used in output filenames and logging labels. |
| `ID` (`--ID`) | Base identifier prepended to generated order-parameter files. |
| `PSF` (`--PSF`) | Topology/reference file used to load the trajectory: a CG PSF for `Q`, `G`, `K`, or an AA PDB for `SASA` and `XP`. |
| `DCD` (`--DCD`) | Trajectory file to analyze. Use the CG DCD for `Q`, `G`, `K` and the AA DCD for `SASA`, `XP`. |
| `outdir` (`--outdir`) | Output directory where OP subfolders and result files will be written. |
| `start` (`--start`) | First frame index to analyze, using 0-based indexing. |
| `ops` (`--ops`) | Order parameters to compute. Allowed values are `Q`, `G`, `K`, `SASA`, `XP`. |
| `CG` (`--CG`) | Boolean flag indicating that the trajectory is coarse-grained. Use `true` / `--CG` for `Q`, `G`, `K` workflows; use `false` / omit the flag for `SASA`, `XP`. If omitted in the config, the wrapper defaults this from `ops`. |
| `Calpha` (`--Calpha`) | Boolean flag indicating that `G` should use C-alpha contacts rather than heavy-atom contacts. Typical Workflow 2 settings are `true` for CG runs and `false` for AA runs. If omitted in the config, the wrapper defaults this from `CG`. |
| `ent_detection_method` (`--ent_detection_method`) | Entanglement-change criterion for `G`: `1` for GLN-only, `2` for TLN-based detection, `3` for GLN and TLN both nonzero on the same terminus. |
| `no_topoly` (`--no_topoly`) | Boolean switch that disables Topoly crossing detection and runs the GLN-only workflow for `G`. |
| `nproc` (`--nproc`) | Number of CPU cores used by the heavier calculations, especially `G` and trajectory-mode `XP`. |
| `COR` (`--COR`) | CG reference coordinate file required when computing `Q`, `G`, or `K`. |
| `sec_elements` (`--sec_elements`) | Secondary-structure definition file required when computing `Q`, `G`, or `K`. |
| `domain` (`--domain`) | Domain-boundary definition file required when computing `Q`, `G`, or `K`. |
| `pdb_file` (`--pdb_file`) | All-atom reference PDB supplied to the `XP` calculation for residue-pair generation and Jwalk SASD calculations. Required when `XP` is included in `ops`. |
| `chunk_frames` (`--chunk_frames`) | Optional number of frames per chunk when writing `Combined_GE` outputs from `G`; use this to reduce memory pressure on long trajectories. |
| `chunk_suffix` (`--chunk_suffix`) | Filename suffix appended to chunked `Combined_GE` pickle outputs. |
| `log_level` (`--log_level`) | Logging verbosity for the wrapper: `DEBUG`, `INFO`, `WARNING`, or `ERROR`. |
| `logdir` (`--logdir`) | Directory where the log file is written. If omitted, logging defaults to the same directory as `outdir`. |

## I/O Reference for run_OP_on_simulation_traj.py

Each file below is listed once, followed by the column-level description table for tabular outputs.

### `$DATASTORE/user_input/reference_structures/1zmr_model_clean_ca.psf`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/user_input/reference_structures/1zmr_model_clean_ca.psf` | CG topology defining atoms/beads and bonds used to interpret CG trajectories. |

### `$DATASTORE/user_input/reference_structures/1zmr_model_clean_ca.cor`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/user_input/reference_structures/1zmr_model_clean_ca.cor` | Reference CG coordinates used as the baseline for Q/G/K calculations. |

### `$DATASTORE/user_input/cg_trajectories/<TRAJ>_prod.dcd`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/user_input/cg_trajectories/<TRAJ>_prod.dcd` | CG production trajectory used to compute Q, G, and K over frames. |

### `$DATASTORE/user_input/reference_structures/secondary_struc_defs.txt`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/user_input/reference_structures/secondary_struc_defs.txt` | Secondary-structure definition file required for Q/G/K setup. |

### `$DATASTORE/user_input/reference_structures/domain_def.dat`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/user_input/reference_structures/domain_def.dat` | Domain-boundary definition file required for Q/G/K setup. |

### `$DATASTORE/user_input/reference_structures/1zmr_model_clean.pdb`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/user_input/reference_structures/1zmr_model_clean.pdb` | All-atom topology used to load the AA trajectory for SASA calculations and as the input PDB for XP/Jwalk calculations. |

### `$DATASTORE/user_input/aa_trajectories/<TRAJ>_prod_aa.dcd`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/user_input/aa_trajectories/<TRAJ>_prod_aa.dcd` | All-atom production trajectory used for SASA and XP calculations. |

### `$OUTDIR/OP/Q/1ZMR_Traj<TRAJ>.Q`

| Column Name | Column Description |
|---|---|
| D_1 | Fraction of native contacts formed within domain/region 1 for the frame. Values range from 0 to 1. |
| D_2 | Fraction of native contacts formed within domain/region 2 for the frame. Values range from 0 to 1. |
| `1\|2` | Fraction of native contacts formed between domains/regions 1 and 2 for the frame. Values range from 0 to 1. |
| total | Overall fraction of native contacts formed in the frame across the full structure. Values range from 0 to 1. |
| Frame | Frame index in the trajectory. |

### `$OUTDIR/OP/G/1ZMR_Traj<TRAJ>.G`

| Column Name | Column Description |
|---|---|
| Frame | Frame index in the trajectory. |
| L-C~ | Fraction or count of contacts in class "loss of linking number and chirality switch" for the frame. |
| L-C# | Fraction or count of contacts in class "loss of linking number with chirality retained" for the frame. |
| L+C~ | Fraction or count of contacts in class "gain of linking number and chirality switch" for the frame. |
| L+C# | Fraction or count of contacts in class "gain of linking number with chirality retained" for the frame. |
| L#C~ | Fraction or count of contacts in class "linking number unchanged with chirality switch" for the frame. |
| L#C# | Fraction or count of contacts in class "no topology change" for the frame. |
| G | Aggregate G order parameter for the frame. |

### `$OUTDIR/OP/K/K_<TRAJ>_prod.dat`

This is a space-delimited text file with one row per analyzed frame, written in trajectory order starting from the requested `start` frame. Unlike the `Q` and `G` outputs, the file does not include an explicit `Frame` column.

| Column Name | Column Description |
|---|---|
| D_1 | Mirror-symmetry score for domain/region 1. Values range from 0 to 1, where lower values are more native-like and higher values indicate stronger mirror-image character. |
| D_2 | Mirror-symmetry score for domain/region 2, on the same 0 to 1 scale. |
| `1\|2` | Mirror-symmetry score for contacts spanning domains/regions 1 and 2. |
| total | Overall mirror-symmetry score for the frame. This is the main value used to flag possible mirror-image artifacts. |

### `$OUTDIR/OP_AA/SASA/1ZMR_Traj<TRAJ>.SASA`

| Column Name | Column Description |
|---|---|
| Time(ns) | Simulation time in nanoseconds for the frame. |
| Frame | Frame index in the trajectory. |
| resid | Residue index associated with the SASA value. |
| SASA(nm^2) | Solvent-accessible surface area for the residue in nm². |

### `$OUTDIR/OP_AA/XP/1ZMR_Traj<TRAJ>.XP`

| Column Name | Column Description |
|---|---|
| Frame | Frame index in the trajectory. |
| Index | Row index for the residue-pair measurement. |
| Model | Model identifier used by Jwalk for this row. |
| Atom1 | First residue/atom label in the cross-linkable pair. |
| Atom2 | Second residue/atom label in the cross-linkable pair. |
| SASD | Solvent-accessible surface distance between Atom1 and Atom2. |
| Euclidean Distance | Straight-line distance between Atom1 and Atom2. |
| XP | Cross-link propensity score derived from distance features. |

### `$OUTDIR/OP/G/Combined_GE/1ZMR_traj<TRAJ>_chunk_<CHUNK>.pkl`

This file is a binary Python pickle, not a flat table. Each chunk stores the native-reference entanglement fingerprints once under the `ref` key, plus a dictionary for each frame in the chunk. It is the main intermediate used by downstream non-native entanglement clustering.

#### Top-level dictionary layout

| Key | Value Type | Description |
|---|---|---|
| `ref` | `dict` | Reference/native entanglement data derived from the static reference structure. |
| `<frame_number>` | `dict` | Per-frame entanglement data for each analyzed trajectory frame in this chunk, keyed by the integer frame index (for example `6600`, `6601`, ...). |

#### Second-level dictionary for `ref` and each frame

| Key | Value Type | Description |
|---|---|---|
| `ent_fingerprint` | `dict` | Native-contact keyed entanglement fingerprint dictionary for the reference structure or frame. |
| `chg_ent_fingerprint` | `dict` or `None` | Per-contact change-of-entanglement fingerprints relative to the reference. This is `None` for `ref` and populated for trajectory frames. |
| `G_dict` | `dict` or `None` | Per-frame counts of the six entanglement-change classes (`L-C~`, `L-C#`, `L+C~`, `L+C#`, `L#C~`, `L#C#`). This is `None` for `ref`. |
| `G` | `float` or `None` | Scalar G order parameter for that frame. This is `None` for `ref`. |

#### `ent_fingerprint` dictionary

The `ent_fingerprint` dictionary is keyed by the native-contact residue pair `(i, j)`. Each value is a dictionary describing the entanglement state of that contact in the reference or in the current frame.

| Nested Key | Value Type | Description |
|---|---|---|
| `native_contact` | `list[int, int]` | The loop-defining native contact as the residue pair `[i, j]`. |
| `linking_value` | `list[float, float]` | Raw Gaussian entanglement values for the N-terminal and C-terminal sides of the loop. |
| `crossing_resid` | `list[list[int], list[int]]` | Crossing residue indices for the N-terminal and C-terminal sides of the loop. |
| `crossing_pattern` | `list[...]` | Crossing chirality/sign pattern aligned to `crossing_resid` for the N-terminal and C-terminal sides. |
| `gauss_linking_number` | `list[int, int]` | Rounded/discrete Gaussian linking numbers for the N-terminal and C-terminal sides. |
| `topoly_linking_number` | `list[int or NaN, int or NaN]` | Topological linking numbers for the N-terminal and C-terminal sides. Values may be `NaN` when not assigned. |

#### `chg_ent_fingerprint` dictionary

For each trajectory frame, `chg_ent_fingerprint` is keyed by the same native-contact tuple `(i, j)`. Each value stores the frame fingerprint side-by-side with the matching reference fingerprint and the resulting change classification.

| Nested Key | Value Type | Description |
|---|---|---|
| `type` | `list[str, str]` | Human-readable change labels for the N-terminal and C-terminal sides, such as `no change` or `loss of linking number & no change of linking chirality`. |
| `code` | `list[str, str]` | Compact change-class codes for the N-terminal and C-terminal sides, using the six-class scheme (`L-C~`, `L-C#`, `L+C~`, `L+C#`, `L#C~`, `L#C#`). |
| `native_contact` | `list[int, int]` | The frame-native contact residue pair `[i, j]`. |
| `linking_value` | `list[float, float]` | Frame Gaussian entanglement values for the N-terminal and C-terminal sides. |
| `crossing_resid` | `list[list[int], list[int]]` | Frame crossing residue indices for the N-terminal and C-terminal sides. |
| `crossing_pattern` | `list[...]` | Frame crossing chirality/sign pattern aligned to `crossing_resid`. |
| `gauss_linking_number` | `list[int, int]` | Frame Gaussian linking numbers for the two termini. |
| `topoly_linking_number` | `list[int or NaN, int or NaN]` | Frame topological linking numbers for the two termini. |
| `ref_native_contact` | `list[int, int]` | Reference/native contact residue pair `[i, j]`. |
| `ref_linking_value` | `list[float, float]` | Reference Gaussian entanglement values for the two termini. |
| `ref_crossing_resid` | `list[list[int], list[int]]` | Reference crossing residue indices for the two termini. |
| `ref_crossing_pattern` | `list[...]` | Reference crossing chirality/sign pattern aligned to `ref_crossing_resid`. |
| `ref_gauss_linking_number` | `list[int, int]` | Reference Gaussian linking numbers for the two termini. |
| `ref_topoly_linking_number` | `list[int or NaN, int or NaN]` | Reference topological linking numbers for the two termini. |
| `ent_detection_method` | `int` | Entanglement-change detection mode used when assigning the change code (`1`, `2`, or `3`). |

#### `G_dict` and `G`

| Key | Value Type | Description |
|---|---|---|
| `G_dict['L-C~']` | `int` | Count of termini showing loss of linking number with chirality switch in the frame. |
| `G_dict['L-C#']` | `int` | Count of termini showing loss of linking number with chirality retained. |
| `G_dict['L+C~']` | `int` | Count of termini showing gain of linking number with chirality switch. |
| `G_dict['L+C#']` | `int` | Count of termini showing gain of linking number with chirality retained. |
| `G_dict['L#C~']` | `int` | Count of termini showing unchanged linking number with chirality switch. |
| `G_dict['L#C#']` | `int` | Count of termini showing no entanglement-state change. |
| `G` | `float` | Aggregate G order parameter for the frame, computed from the non-`L#C#` change counts and normalized by twice the number of native contacts in the reference structure. |

---

## Step 3. Identify and remove artificial mirror conformations

Before clustering entanglement changes, remove trajectories or frames that are mirror-image artifacts.

### Recommended procedure

1. Load Q and K time series for all trajectories.
2. Identify trajectories where K is persistently low (<=0.6) and Q is >=0.2.
3. Visually inspect flagged frames in VMD for mirror-image conformations.
4. Record the trajectory numbers that are confirmed mirrors in a list (used as `rm_traj_list` in [Workflow 3: Sim-to-Experiment](workflow3_sim2exp.md)).

In this tutorial we identified the mirror artifact trajectories as: 65, 75, 155, 162, 199, 231, 264, 286, 296, 314, 354, 417, 448, 472, 473, 474, 577, 579, 591, 703, 704, 732, 758, 812, 833, 870, 876, 944, 967.

> **Important:** The cutoff values for flagging mirror conformations must be tuned for your system by examining the Q and K distributions.

---

## Step 4. Cluster changes of NCLE status to remove redundancy

This step identifies non-redundant changes in entanglement topology across all trajectories. The pkl file paths are specified in the trajectory-to-pkl mapping CSV (`trajnum2pklfile_path`), which serves as the single source of truth for which pkl files to analyze.

### `trajnum2file.txt` — maps trajectory numbers to pkl files

This file maps trajectory numbers to their **G-order-parameter pkl files** (not DCD files). The format is comma-separated:

```
trajnum,pklfile
<trajectory_number>,<path_to_pkl_file>
```

The pre-populated copy at `$DATASTORE/user_input/metadata/trajnum2file.txt` maps all 1000 trajectories to the Combined_GE pkl files in the DATASTORE. It can be regenerated if needed:

```bash
echo "trajnum,pklfile" > $DATASTORE/user_input/metadata/trajnum2file.txt
for pkl in $DATASTORE/outputs/workflow2/OP_last67/G/Combined_GE/*.pkl; do
    num=$(basename $pkl | sed 's/1ZMR_traj\([0-9]*\)_GE.pkl/\1/')
    echo "$num,$pkl"
done >> $DATASTORE/user_input/metadata/trajnum2file.txt
```
If you used the `--chunk_frames` argument when calculating G, the `.pkl` files will be chunked and you must specify which chunk to use for each trajectory. EntDetect currently does not support using multiple chunks for a single trajectory number.

> **Memory warning:** This step can require tens of gigabytes of RAM. Run on a high-memory node or reduce the number of frames per trajectory in the clustering pool.

```python
from EntDetect.clustering import ClusterNonNativeEntanglements

# ── Paths ──────────────────────────────────────────────────────────────────
DATASTORE   = "/scratch/ims86/EntDetect_Datastore"
OUTDIR      = f"{DATASTORE}/outputs/workflow2"
traj_dir_prefix = f"{DATASTORE}/user_input/cg_trajectories"

# ── Inputs ──────────────────────────────────────────────────────────────────
trajnum2pklfile_path = f"{DATASTORE}/user_input/metadata/OP_last67_trajnum2file.txt"
clust_outdir         = f"{OUTDIR}/nonnative_clustering_last67"

# ── Initialize and Run ──────────────────────────────────────────────────────
clustering_NNents = ClusterNonNativeEntanglements(
    trajnum2pklfile_path=trajnum2pklfile_path,
    traj_dir_prefix=traj_dir_prefix,
    outdir=clust_outdir,
)
clustering_NNents.cluster()
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

---

## Minimal Workflow – 3: Running the change in NCLE status clustering as single script

The `scripts/run_nonnative_entanglement_clustering.py` script automates the clustering workflow.

You can provide parameters either directly as CLI flags, via a `--config` file, or both.
When both are provided, **CLI flags override config values** for the same parameter.

```bash
# ── Activation ──────────────────────────────────────────────────────────────────
source ~/.bashrc
conda activate entdetect

# ── Config file ─────────────────────────────────────────────────────────────────
CFG=scripts/configs/workflow2_nonnative_clustering_config.json

# ── Run ───────────────────────────────────────────────────────────────
python scripts/run_nonnative_entanglement_clustering.py \
    --config $CFG

# Example CLI override (swap in a different trajnum2file):
python scripts/run_nonnative_entanglement_clustering.py \
    --config $CFG \
    --trajnum2pklfile_path /scratch/ims86/EntDetect_Datastore/outputs/workflow2/OP_last67/trajnum2file.txt
```

Config file (matches `scripts/configs/workflow2_nonnative_clustering_config.json`):

```json
{
  "outdir": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/nonnative_clustering_last67",
  "trajnum2pklfile_path": "/scratch/ims86/EntDetect_Datastore/user_input/metadata/OP_last67_trajnum2file.txt",
  "traj_dir_prefix": "/scratch/ims86/EntDetect_Datastore/user_input/cg_trajectories",
  "start_frame": 0,
  "end_frame": 9999999,
  "nproc": 4,
  "log_level": "INFO",
  "logdir": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/nonnative_clustering_last67/logs"
}
```

The JSON/YAML config keys and their matching CLI flags are listed below. For the command-line wrapper, use the same key name with a `--` prefix; the only extra wrapper-only flag is `--config`.

| Term | Definition |
|------|------------|
| `outdir` (`--outdir`) | Output directory where clustering results and derived files are written. |
| `trajnum2pklfile_path` (`--trajnum2pklfile_path`) | CSV manifest mapping each trajectory number to the `Combined_GE` pickle file that should be clustered for that trajectory. This file is the source of truth for which entanglement pickles are analyzed. |
| `traj_dir_prefix` (`--traj_dir_prefix`) | Directory containing the CG trajectory DCD files. These paths are used to resolve representative structures and populate downstream trajectory-to-file mappings. |
| `start_frame` (`--start_frame`) | First frame index to include when clustering, using 0-based indexing. |
| `end_frame` (`--end_frame`) | Last frame index to include when clustering. The example config uses a large sentinel value so all available frames are included. |
| `nproc` (`--nproc`) | Number of parallel worker threads used when loading pickles and clustering entanglement keywords. |
| `log_level` (`--log_level`) | Logging verbosity for the wrapper: `DEBUG`, `INFO`, `WARNING`, or `ERROR`. |
| `logdir` (`--logdir`) | Directory where the clustering log file is written. If omitted, logging defaults to the same directory as `outdir`. |

Submit through SLURM with the dedicated OP_last67 wrapper:

```bash
sbatch assets/slurm/scripts/run_nonnative_clustering_OP_last67.slurm
```

## I/O Reference for run_nonnative_entanglement_clustering.py

Each file below is listed once, followed by column-level details when applicable.

### `$DATASTORE/user_input/metadata/OP_last67_trajnum2file.txt`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/user_input/metadata/OP_last67_trajnum2file.txt` | Mapping of trajectory numbers to OP_last67 Combined_NCLE pickle files used as the source of truth for which pkl files to analyze. |

| Column Name | Column Description |
|---|---|
| trajnum | Integer trajectory identifier (for example 1..1000). |
| pklfile | Absolute path to the trajectory-specific Combined_NCLE pickle file. |

### `$DATASTORE/user_input/cg_trajectories/<TRAJ>_prod.dcd`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/user_input/cg_trajectories/<TRAJ>_prod.dcd` | CG trajectory used to resolve per-frame DCD coordinates; stored in `idx2trajfile` inside the output npz for downstream structure extraction. |

### `$OUTDIR/nonnative_clustering_last67/rep_chg_ent_topoly_linking_number.csv`

| Column Name | Column Description |
|---|---|
| State ID | 1-based index of the representative entanglement-change cluster. |
| Keywords | String encoding the entanglement type and topology classification for the cluster. |
| Trajectory | Path to the DCD file for the trajectory that supplied the representative frame. |
| Frame | Frame index within that trajectory for the representative structure. |
| Native Contact (Residue Index) | Loop-defining residue pair `(i, j)` for the representative entanglement. |
| Ref N-ter Crossing | N-terminal crossing residues in the native reference conformation. |
| Ref C-ter Crossing | C-terminal crossing residues in the native reference conformation. |
| N-ter Crossing | N-terminal crossing residues in the representative trajectory frame. |
| C-ter Crossing | C-terminal crossing residues in the representative trajectory frame. |
| Ref N-ter GLN | Gaussian linking number for the N-terminus in the native reference. |
| Ref C-ter GLN | Gaussian linking number for the C-terminus in the native reference. |
| N-ter GLN | Gaussian linking number for the N-terminus in the representative frame. |
| C-ter GLN | Gaussian linking number for the C-terminus in the representative frame. |
| Ref N-ter Linking Number | Topological linking number for the N-terminus in the native reference. |
| Ref C-ter Linking Number | Topological linking number for the C-terminus in the native reference. |
| N-ter Linking Number | Topological linking number for the N-terminus in the representative frame. |
| C-ter Linking Number | Topological linking number for the C-terminus in the representative frame. |

### `$OUTDIR/nonnative_clustering_last67/chg_ent_struct_topoly_linking_number.csv`

| Column Name | Column Description |
|---|---|
| State ID | 1-based index of the entanglement-change structural state. |
| Rep_chg_ents | List of representative entanglement-change cluster indices (1-based) that define this structural state. |
| Num of structures | Number of trajectory frames assigned to this structural state. |
| Probability | Fraction of total analyzed frames occupying this state. |
| Rep trajectory | Path to the DCD file containing the representative frame for this state. |
| Rep frame | Frame index of the representative structure within that trajectory. |
| Max Q | Maximum Q order-parameter value among all frames in this state. |
| Median Q | Median Q order-parameter value among all frames in this state. |

### `$OUTDIR/nonnative_clustering_last67/cluster_data_topoly_linking_number.npz`

| I/O | File | File Description |
|---|---|---|
| Output | `$OUTDIR/nonnative_clustering_last67/cluster_data_topoly_linking_number.npz` | Compressed archive containing all clustering artifacts: per-frame entanglement fingerprints, Q values, cluster assignments (`dtrajs`), representative structures, trajectory-to-file mapping (`idx2trajfile`), and structural state data. Used as input to Workflow 3. |

This is a NumPy `.npz` archive, not a flat table. It stores the full in-memory clustering state needed for downstream representative-structure extraction and Workflow 3 analyses.

| Key | Value Type | Description |
|---|---|---|
| `chg_ent_fingerprint_list` | `list[dict]` | Per-trajectory dictionaries of frame-indexed changed-entanglement fingerprints loaded from the `Combined_GE` pickle inputs. |
| `Q_list` | `list[dict]` | Per-trajectory mapping from frame index to Q value, used when selecting representative structural states. |
| `chg_ent_keyword_dict` | `dict[str, list]` | Mapping from entanglement-change keyword to all observed `(traj_idx, frame_idx, i, j)` observations carrying that keyword. |
| `chg_ent_keyword_list` | `list[str]` | Sorted list of unique entanglement-change keywords used as the top-level clustering categories. |
| `idx2trajfile` | `list[str]` | Mapping from internal trajectory index to the absolute DCD path for that trajectory. |
| `idx2frame` | `list[list[int]]` | Mapping from internal trajectory index to the analyzed frame numbers retained from that trajectory. |
| `ent_cluster_data` | `dict[str, list]` | Keyword-grouped representative entanglement clusters; each cluster contains the member observations assigned to that cluster. |
| `ent_cluster_tree` | `dict[str, list]` | Hierarchical clustering trace for each keyword, used to write `cluster_tree_topoly_linking_number.dat`. |
| `rep_chg_ent_list` | `list[list[int]]` | Representative entanglement observation for each entanglement-change cluster, stored as `[traj_idx, frame_idx, i, j]`. |
| `dtrajs` | `object array` | Discrete per-frame cluster assignments for each trajectory, where each frame stores the list of entanglement-cluster IDs present in that frame. |
| `rep_chg_ent_dtrajs` | `object array` | Per-frame representative entanglement fingerprints, keyed by cluster ID, for downstream structure visualization/extraction. |
| `sorted_chg_ent_structure_keyword_list` | `list[str]` | Structural-state labels, sorted by population, where each label encodes the set of representative entanglement clusters present in that state. |
| `chg_ent_structure_cluster_data` | `dict[str, list]` | Mapping from each structural-state label to the list of `(traj_idx, frame_idx)` observations assigned to that state. |
| `rep_struct_data` | `dict[str, list[int, int]]` | Representative `(traj_idx, frame_idx)` chosen for each structural state, selected by maximum Q within that state. |

### `$OUTDIR/nonnative_clustering_last67/cluster_tree_topoly_linking_number.dat`

| I/O | File | File Description |
|---|---|---|
| Output | `$OUTDIR/nonnative_clustering_last67/cluster_tree_topoly_linking_number.dat` | Plain-text representation of the hierarchical clustering tree showing how individual entanglement-change observations were merged into representative clusters. |

This is a human-readable text summary of the staged clustering procedure for each entanglement-change keyword.

| Section | Description |
|---|---|
| Keyword header line | The entanglement-change keyword being clustered. |
| `After clustering on N crossing:` | Cluster memberships after grouping by N-terminal crossing similarity. |
| `After clustering on C crossing:` | Cluster memberships after grouping by C-terminal crossing similarity. |
| `After clustering on loop:` | Final cluster memberships after loop-level clustering/merging. |
| Cluster entries like `[1, 4, 7]` | Representative entanglement-cluster IDs grouped together at that stage. These IDs correspond to the 1-based cluster numbering used in `rep_chg_ent_topoly_linking_number.csv`. |

### `$OUTDIR/nonnative_clustering_last67/rep_chg_ent_list_topoly_linking_number.pkl`

| I/O | File | File Description |
|---|---|---|
| Output | `$OUTDIR/nonnative_clustering_last67/rep_chg_ent_list_topoly_linking_number.pkl` | Serialized Python list of representative entanglement-change objects, one per cluster, containing the full fingerprint data for the chosen representative frame. |

This is a Python pickle containing the minimal representative observation for each entanglement-change cluster.

| Item Layout | Value Type | Description |
|---|---|---|
| List element | `list[int, int, int, int]` | One representative entanglement observation per cluster. |
| Element `[0]` | `int` | Internal trajectory index (`traj_idx`) used within the clustering run. Resolve to the DCD path with `idx2trajfile` from the `.npz` archive. |
| Element `[1]` | `int` | Representative frame index within that trajectory. |
| Element `[2]` | `int` | Native-contact residue index `i` for the representative changed entanglement. |
| Element `[3]` | `int` | Native-contact residue index `j` for the representative changed entanglement. |

The richer per-contact fingerprint metadata for these representatives is written in tabular form to `rep_chg_ent_topoly_linking_number.csv` and is also recoverable through `chg_ent_fingerprint_list` in the `.npz` archive.

### `$OUTDIR/nonnative_clustering_last67/chg_ent_topoly_linking_number_distribution.pdf`

| I/O | File | File Description |
|---|---|---|
| Output | `$OUTDIR/nonnative_clustering_last67/chg_ent_topoly_linking_number_distribution.pdf` | Plot showing the loop and crossing residue positions for each cluster, coloured by reference (green) vs. trajectory (blue) crossings. |

This PDF visualizes the residue-level geometry of each representative entanglement-change cluster.

| Plot Element | Description |
|---|---|
| X-axis | 1-based entanglement-change cluster ID. |
| Y-axis | Residue index. |
| Red vertical line | The loop span defined by the native-contact pair `(i, j)` for a sampled member of the cluster. |
| Green tick marks | Reference/native crossing residues for that sampled cluster member. |
| Blue tick marks | Crossing residues observed in the trajectory frame for that sampled cluster member. |

When a cluster contains many members, the plot shows a sampled subset rather than every structure, so it is best used as a qualitative overview of loop and crossing-position diversity across clusters.

---

## Step 5. Build a Markov state model of the GvQ probability surface

Organize simulation frames into microstates and metastable states using the full-trajectory order-parameter data using parameters from `assets/slurm/scripts/run_MSM.slurm`.

### Build the MSM with Python

```python
from EntDetect.clustering import MSMNonNativeEntanglementClustering

# ── Paths ────────────────────────────────────────────────
DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow2"
outdir         = f"{OUTDIR}/MSM"
OPpath         = f"{OUTDIR}/OP/"

# ── Inputs ────────────────────────────────────────────────
ID             = "1ZMR_prod"
n_large_states = 10   # number of metastable macro-states requested
lagtime        = 20   # lag time in frames

# ── Initialize and Run ─────────────────────────────────────────
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

### Using the command-line interface

For convenience, use the `scripts/run_MSM.py` script directly (see "Minimal Workflow – 4: Markov state modeling" below for full details).

### Critical notes

- Try **multiple values** of `n_large_states` (e.g., 5, 10, 15). The protocol notes ≤15 often works well.
- The final number of states may be lower than requested if empty states are discarded.
- The default lag time of 1 frame is appropriate for visualization and exploratory grouping. For kinetic interpretation, test for Markovian behavior explicitly and choose lag time accordingly.

### Expected outputs

This step produces a frame-level MSM mapping CSV with microstate and metastable-state assignments, a metastable-state composition CSV listing the microstates grouped into each metastable state, a `.npy` array containing metastable-state membership distributions over all MSM microstates for downstream weighting/aggregation analyses, and a PNG summary figure showing both the free-energy landscape in `(Q, G)` space and the metastable-state assignment map.

| File | Contents |
|------|----------|
| `1ZMR_prod_MSMmapping.csv` | Per-frame microstate and metastable-state assignments |
| `1ZMR_prod_meta_set.csv` | Metastable-state summary |
| `1ZMR_prod_meta_dist.npy` | Metastable-state probability distribution |
| `1ZMR_prod_StateAndFEplot.png` | Order-parameter landscape and state assignments |

---

## Minimal Workflow – 4: Markov state modeling

The `scripts/run_MSM.py` script automates MSM construction.

You can provide parameters either directly as CLI flags, via a `--config` file, or both.
When both are provided, **CLI flags override config values** for the same parameter.

```bash
# ── Activation ──────────────────────────────────────────────────────────────────
source ~/.bashrc
conda activate entdetect

# ── Config file ─────────────────────────────────────────────────────────────────
CFG=scripts/configs/workflow2_MSM_config.json

# ── Run ───────────────────────────────────────────────────────────────
python scripts/run_MSM.py --config $CFG

# Example CLI override (try a different number of metastable states):
python scripts/run_MSM.py --config $CFG --n_large_states 15
```

Config file (matches `scripts/configs/workflow2_MSM_config.json`):

```json
{
  "outdir": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/MSM",
  "ID": "1ZMR_prod",
  "OPpath": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/OP/",
  "start": 0,
  "n_large_states": 10,
  "lagtime": 20,
  "rm_traj_list": [65, 75, 155, 162, 199, 231, 264, 286, 296, 314, 354, 417, 448, 472, 473, 474, 577, 579, 591, 703, 704, 732, 758, 812, 833, 870, 876, 944, 967],
  "log_level": "INFO",
  "logdir": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/MSM/logs"
}
```

The JSON/YAML config keys and their matching CLI flags are listed below. For the command-line wrapper, use the same key name with a `--` prefix; the only extra wrapper-only flag is `--config`.

| Term | Definition |
|------|------------|
| `outdir` (`--outdir`) | Output directory where MSM results and summary plots are written. |
| `ID` (`--ID`) | Base identifier used as the prefix for all MSM output filenames. |
| `OPpath` (`--OPpath`) | Path to the order-parameter output directory containing the per-trajectory `Q` and `G` files used as MSM features. |
| `start` (`--start`) | First frame index to include when building the MSM, using 0-based indexing. |
| `n_large_states` (`--n_large_states`) | Requested number of metastable states for the largest connected MSM component. Empty states may be discarded, so the final number can be smaller. |
| `lagtime` (`--lagtime`) | MSM lag time in frames used to build the transition model. |
| `rm_traj_list` (`--rm_traj_list`) | Optional list of trajectory numbers to exclude before MSM construction, typically mirror-image artifacts identified from the `K` analysis. |
| `log_level` (`--log_level`) | Logging verbosity for the wrapper: `DEBUG`, `INFO`, `WARNING`, or `ERROR`. |
| `logdir` (`--logdir`) | Directory where the MSM log file is written. If omitted, logging defaults to the same directory as `outdir`. |

For full production MSM:

```bash
sbatch assets/slurm/scripts/run_MSM.slurm
```

## I/O Reference for run_MSM.py

Each file below is listed once, followed by column-level details for CSV outputs and structure details for non-tabular outputs.

### `$OUTDIR/OP/Q/1ZMR_Traj<TRAJ>.Q`

| I/O | File | File Description |
|---|---|---|
| Input | `$OUTDIR/OP/Q/1ZMR_Traj<TRAJ>.Q` | Per-trajectory Q order-parameter time series used as MSM features. |

### `$OUTDIR/OP/G/1ZMR_Traj<TRAJ>.G`

| I/O | File | File Description |
|---|---|---|
| Input | `$OUTDIR/OP/G/1ZMR_Traj<TRAJ>.G` | Per-trajectory G/topology-change time series used as MSM features. |

### `$OUTDIR/MSM/<ID>_MSMmapping.csv`

| Column Name | Column Description |
|---|---|
| traj | Trajectory identifier for the row. |
| frame | Frame index within the trajectory. |
| microstate | Assigned MSM microstate label for this frame. |
| metastablestate | Assigned metastable-state label for this frame. |
| Q | Q order-parameter value used in state construction/visualization. |
| G | G order-parameter value used in state construction/visualization. |

### `$OUTDIR/MSM/<ID>_meta_set.csv`

| Column Name | Column Description |
|---|---|
| metastable_state | Metastable-state identifier. |
| microstates | Collection of microstate IDs grouped into the metastable state. |

### `$OUTDIR/MSM/<ID>_meta_dist.npy`

| I/O | File | File Description |
|---|---|---|
| Output | `$OUTDIR/MSM/<ID>_meta_dist.npy` | NumPy array of metastable-state membership distributions over MSM microstates. Used downstream to map metastable populations onto microstate-indexed observables. |

This file is a NumPy `.npy` array (not a flat table). It stores one probability distribution per metastable state.

| Axis / Element | Value Type | Description |
|---|---|---|
| `meta_dist` | `np.ndarray` (`float`) | 2D array with shape `(n_metastable_states, n_microstates)`. |
| `meta_dist[s, m]` | `float` | Probability weight for microstate `m` within metastable state `s`. |
| Row index `s` | `int` | Metastable-state ID, consistent with `metastablestate` labels in `<ID>_MSMmapping.csv` and `metastable_state` in `<ID>_meta_set.csv`. |
| Column index `m` | `int` | Global microstate ID from MSM clustering, consistent with `microstate` in `<ID>_MSMmapping.csv`. |

| Property | Description |
|---|---|
| Row normalization | Each metastable-state row is normalized so that probabilities over its assigned microstates sum to 1 (rows contain zeros for microstates not in that metastable state). |
| Sparse pattern | Most entries are zero because each metastable state occupies only a subset of all microstates. |
| Fallback behavior | If a connected component collapses to a single metastable state, the row is estimated from observed microstate occupancies in that component. |

### `$OUTDIR/MSM/<ID>_StateAndFEplot.png`

| I/O | File | File Description |
|---|---|---|
| Output | `$OUTDIR/MSM/<ID>_StateAndFEplot.png` | Two-panel summary figure generated during MSM construction: free-energy landscape in (Q, G) space and metastable-state assignment map over the same coordinates. |

| Plot Element | Description |
|---|---|
| Left panel | 2D free-energy surface from the frame density in `(Q, G)` coordinates. |
| Right panel | Scatter of frames in `(Q, G)` colored by assigned metastable state. |
| Colorbar | Discrete metastable-state IDs used in the right panel. |
| Axes | `Q` and `G` order parameters used as MSM input features. |

---

## Step 5b. Label MSM Data and Define Analysis Cases

Before computing statistics or folding pathways, each trajectory needs a **type label** (e.g. `A` or `B`) representing the biological comparison of interest. Any column with consistent string labels per trajectory can serve as the type column downstream.

This tutorial demonstrates two contrasting cases:

| Case | Labeling rule | Expected signal |
|------|---------------|----------------|
| **Case 1 — Biologically-informed** | A = max Q ≥ 0.80 **and** max G ≤ 0.05 (native-like, non-entangled); B = all others | High JS divergence |
| **Case 2 — Random (negative control)** | A/B assigned randomly per trajectory (seed=42) | JS divergence ≈ 0 |

> **Define your own cases:** Replace the labeling rule below with whatever biological comparison makes sense for your data, such as different temperature conditions, mutation variants, or folding outcomes. The only requirement is a column with consistent string labels per trajectory.

### Generate annotated MSM mapping files

```python
import pandas as pd
import numpy as np

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow2"

msm_mapping = pd.read_csv(f"{OUTDIR}/MSM/1ZMR_prod_MSMmapping.csv")
```

**Case 1 — Biologically-informed split (Q ≥ 0.8 and G ≤ 0.05):**

```python
# A = trajectories whose max Q >= 0.80 AND max G <= 0.05 (native-like and non-entangled)
# B = all others (misfolded or with significant entanglement changes)
traj_stats = msm_mapping.groupby('traj').agg(max_Q=('Q', 'max'), max_G=('G', 'max'))
native_trajs = traj_stats[(traj_stats['max_Q'] >= 0.80) & (traj_stats['max_G'] <= 0.05)].index
msm_mapping['traj_type'] = msm_mapping['traj'].isin(native_trajs).map({True: 'A', False: 'B'})

annotated_file_case1 = f"{OUTDIR}/MSM/1ZMR_prod_MSMmapping_QG_native.csv"
msm_mapping[['traj', 'frame', 'microstate', 'metastablestate', 'Q', 'G', 'traj_type']].to_csv(
    annotated_file_case1, index=False)
print("Case 1 — trajectory type distribution:")
print(msm_mapping.groupby('traj')['traj_type'].first().value_counts())
```

> **Note:** The Q and G thresholds are system-specific. Adjust to match the folding and entanglement distributions of your system.

**Case 2 — Random split (negative control):**

```python
# Randomly assign A/B labels to trajectories (fixed seed for reproducibility)
rng = np.random.default_rng(seed=42)
all_trajs = msm_mapping['traj'].unique()
random_labels = dict(zip(all_trajs, rng.choice(['A', 'B'], size=len(all_trajs))))
msm_mapping['traj_type'] = msm_mapping['traj'].map(random_labels)

annotated_file_case2 = f"{OUTDIR}/MSM/1ZMR_prod_MSMmapping_random.csv"
msm_mapping[['traj', 'frame', 'microstate', 'metastablestate', 'Q', 'G', 'traj_type']].to_csv(
    annotated_file_case2, index=False)
print("Case 2 — random type distribution:")
print(msm_mapping.groupby('traj')['traj_type'].first().value_counts())
```

---

## Step 6. Visualize state distribution, state probability evolution, representative state structures and folding pathways (OPTIONAL)

Using `MSMStats`, compute how each trajectory population (A and B) distributes across metastable states over simulation time. Run separately for each labelled case to compare the biological signal against the random baseline.

**Input:** Annotated MSM mapping CSV files from Step 5b  
**Output:** `$OUTDIR/MSM_StateProbabilityStats_{case}/` — probability plots and summary tables

### Compute state probability statistics

```python
from EntDetect.statistics import MSMStats
import os

# ── Paths ──────────────────────────────────────────────────────────────────
DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow2"
outdir_stats_case1 = f"{OUTDIR}/MSM_StateProbabilityStats_QG_native"
outdir_stats_case2 = f"{OUTDIR}/MSM_StateProbabilityStats_random"

# ── Inputs ──────────────────────────────────────────────────────────────────
traj_type_list = ['A', 'B']
rm_traj_list   = [65, 75, 155, 162, 199, 231, 264, 286, 296, 314, 354, 417, 448, 472, 473, 474, 577, 579, 591, 703, 704, 732, 758, 812, 833, 870, 876, 944, 967]
msm_data_case1 = f"{OUTDIR}/MSM/1ZMR_prod_MSMmapping_QG_native.csv"
msm_data_case2 = f"{OUTDIR}/MSM/1ZMR_prod_MSMmapping_random.csv"

# ── Initialize once, run both cases ─────────────────────────────────────────────────────
MS = MSMStats(rm_traj_list=rm_traj_list)
# `outdir` in MSMStats(...) is the fallback default; per-case method `outdir` values below control where each case is written.

df1 = MS.StateProbabilityStats(
    msm_data_file=msm_data_case1,
    traj_type_col='traj_type',
    traj_type_list=traj_type_list,
    outdir=outdir_stats_case1,
)
MS.Plot_StateProbabilityStats(
    df=df1,
    traj_type_col='traj_type',
    traj_type_list=traj_type_list,
    outdir=outdir_stats_case1,
)

# ── Run Case 2 — Random split (negative control) ────────────────────────────────────────
df2 = MS.StateProbabilityStats(
    msm_data_file=msm_data_case2,
    traj_type_col='traj_type',
    traj_type_list=traj_type_list,
    outdir=outdir_stats_case2,
)
MS.Plot_StateProbabilityStats(
    df=df2,
    traj_type_col='traj_type',
    traj_type_list=traj_type_list,
    outdir=outdir_stats_case2,
)
```

**Output:** Time series plots showing the population (probability) of each metastable state over simulation time, one plot per case. Useful for identifying which states are transiently vs persistently populated by each trajectory subpopulation.

---

## Minimal Workflow – 5: Analyzing MSM metadata state probability

The `scripts/run_MSMStats.py` script computes state probability statistics in a single call. It requires the MSM mapping CSV to already contain the trajectory-type column (see Step 5b). Run it separately for each case.

You can provide parameters either directly as CLI flags, via a `--config` file, or both.
When both are provided, **CLI flags override config values** for the same parameter.

```bash
# ── Activation ──────────────────────────────────────────────────────────────────
source ~/.bashrc
conda activate entdetect

# ── Config files ─────────────────────────────────────────────────────────────────
CFG_CASE1=scripts/configs/workflow2_MSMStats_case1_config.json
CFG_CASE2=scripts/configs/workflow2_MSMStats_case2_config.json

# ── Run ───────────────────────────────────────────────────────────────
# Case 1 — biologically-informed split (Q >= 0.8 and G <= 0.05)
python scripts/run_MSMStats.py --config $CFG_CASE1

# Case 2 — random split (negative control)
python scripts/run_MSMStats.py --config $CFG_CASE2
```

Config file example 1 (matches `scripts/configs/workflow2_MSMStats_case1_config.json`):

```json
{
  "msm_data_file": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/MSM/1ZMR_prod_MSMmapping_QG_native.csv",
  "outdir": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/MSM_StateProbabilityStats_QG_native",
  "traj_type_col": "traj_type",
  "traj_type_list": ["A", "B"],
  "rm_traj_list": [65, 75, 155, 162, 199, 231, 264, 286, 296, 314, 354, 417, 448, 472, 473, 474, 577, 579, 591, 703, 704, 732, 758, 812, 833, 870, 876, 944, 967],
  "log_level": "INFO",
  "logdir": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/MSM_StateProbabilityStats_QG_native/logs"
}
```

Config file example 2 (matches `scripts/configs/workflow2_MSMStats_case2_config.json`):

```json
{
  "msm_data_file": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/MSM/1ZMR_prod_MSMmapping_random.csv",
  "outdir": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/MSM_StateProbabilityStats_random",
  "traj_type_col": "traj_type",
  "traj_type_list": ["A", "B"],
  "rm_traj_list": [65, 75, 155, 162, 199, 231, 264, 286, 296, 314, 354, 417, 448, 472, 473, 474, 577, 579, 591, 703, 704, 732, 758, 812, 833, 870, 876, 944, 967],
  "log_level": "INFO",
  "logdir": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/MSM_StateProbabilityStats_random/logs"
}
```

The JSON/YAML config keys and their matching CLI flags are listed below. For the command-line wrapper, use the same key name with a `--` prefix; the only extra wrapper-only flag is `--config`.

| Term | Definition |
|------|------------|
| `msm_data_file` (`--msm_data_file`) | Input MSM mapping CSV produced by Step 5b, containing frame-level state assignments plus the trajectory-type label column used for group comparison. |
| `outdir` (`--outdir`) | Output directory where state-probability statistics tables and plots are written. |
| `traj_type_col` (`--traj_type_col`) | Column name in `msm_data_file` that stores trajectory population labels (for example `traj_type` with values `A` and `B`). |
| `traj_type_list` (`--traj_type_list`) | Ordered list of trajectory-type labels to include in the probability analysis and plotting. |
| `rm_traj_list` (`--rm_traj_list`) | Optional list of trajectory IDs to exclude before computing state-probability statistics (for example known mirror artifacts). |
| `log_level` (`--log_level`) | Logging verbosity for the wrapper: `DEBUG`, `INFO`, `WARNING`, or `ERROR`. |
| `logdir` (`--logdir`) | Directory where the `run_MSMStats` log file is written. If omitted, logging defaults to `outdir`. |

## I/O Reference for run_MSMStats.py

Each file below is listed once, followed by column-level details for tabular outputs.

### `$OUTDIR/MSM/1ZMR_prod_MSMmapping_<case>.csv`

| I/O | File | File Description |
|---|---|---|
| Input | `$OUTDIR/MSM/1ZMR_prod_MSMmapping_<case>.csv` | Labeled frame-level MSM mapping with trajectory-type column used to compute population-specific state probabilities. |

### `$OUTDIR/MSM_StateProbabilityStats_<case>/MSTS.csv`

| I/O | File | File Description |
|---|---|---|
| Output | `$OUTDIR/MSM_StateProbabilityStats_<case>/MSTS.csv` | Long-format state-probability summary table containing one row per `(traj_type, time, state)` with bootstrap confidence intervals. |

| Column Name | Column Description |
|---|---|
| traj_type | Trajectory-population label (for example A or B). |
| Time(s) | Simulation time in seconds for the frame index. |
| State | Metastable-state identifier. |
| Probability | Population probability for this state at this time, averaged across all trajectories of the given type. |
| Lower CI | Lower bound of the bootstrap 95% confidence interval for state probability. |
| Upper CI | Upper bound of the bootstrap 95% confidence interval for state probability. |

### `$OUTDIR/MSM_StateProbabilityStats_<case>/<traj_type>_MSTS_plot.png`

| I/O | File | File Description |
|---|---|---|
| Output | `$OUTDIR/MSM_StateProbabilityStats_<case>/<traj_type>_MSTS_plot.png` | One PNG per trajectory type, showing metastable-state probability trajectories over time with 95% bootstrap confidence intervals. |

| Plot Element | Description |
|---|---|
| X-axis | Simulation time in seconds (`Time(s)`). |
| Y-axis | State population probability. |
| Colored lines | Mean probability trajectory for each metastable state. |
| Shaded bands | Bootstrap 95% confidence interval (`Lower CI`, `Upper CI`) around each state trajectory. |
| Per-file scope | Each file contains only one trajectory type (for example `A` or `B`), with all metastable states overlaid. |

---

### Visualize representative metastable structures

Identify the representative frame for each metastable state from `1ZMR_prod_MSMmapping.csv`, extract the corresponding structure from the trajectory, and inspect it in VMD.

## Step 7. Folding Pathways and Jensen-Shannon Divergence (OPTIONAL)

Analyse how trajectories transition between metastable states and quantify the divergence between the two populations using Jensen-Shannon (JS) divergence.

We demonstrate **two contrasting cases** to illustrate what the JS signal looks like under biologically meaningful vs. meaningless partitioning.

| Case | Labeling rule | Expected JS signal |
|------|---------------|-------------------|
| **Case 1 — Biologically-informed split** | A = correctly folded (max Q ≥ 0.8 **and** max G ≤ 0.05); B = misfolded/entangled | High divergence — A and B follow distinct metastable-state progressions |
| **Case 2 — Random split (negative control)** | A/B assigned randomly per trajectory (fixed seed) | Near-zero divergence — populations are statistically identical by construction |

Running both cases back-to-back makes the biological signal immediately recognizable against baseline noise.

Run `FoldingPathwayStats` for each case. Both methods take the annotated MSM mapping CSV path directly and load it internally. `post_trans()` traces state-to-state transitions for each trajectory, removing loops to yield the minimal directed pathway, while `JS_divergence()` computes a windowed Jensen-Shannon divergence between the two populations over simulation time.

**Case 1 — Biologically-informed split:**

```python
from EntDetect.statistics import FoldingPathwayStats

# ── Paths ──────────────────────────────────────────────────────────────────
DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow2"
outdir_case1   = f"{OUTDIR}/FoldingPathway_QG_native"
outdir_case2   = f"{OUTDIR}/FoldingPathway_random"

# ── Inputs ──────────────────────────────────────────────────────────────────
msm_data_file_case1 = f"{OUTDIR}/MSM/1ZMR_prod_MSMmapping_QG_native.csv"
msm_data_file_case2 = f"{OUTDIR}/MSM/1ZMR_prod_MSMmapping_random.csv"
meta_set_file  = f"{OUTDIR}/MSM/1ZMR_prod_meta_set.csv"
rm_traj_list = [65, 75, 155, 162, 199, 231, 264, 286, 296, 314, 354, 417, 448, 472, 473, 474, 577, 579, 591, 703, 704, 732, 758, 812, 833, 870, 876, 944, 967],

# ── Initialize once, run both cases ─────────────────────────────────────────────────────
FP = FoldingPathwayStats(outdir=f"{OUTDIR}/FoldingPathway", rm_traj_list=rm_traj_list)

folding_pathways_case1 = FP.post_trans(
    msm_data_file=msm_data_file_case1,
    traj_type_col='traj_type',
    traj_type_list=['A', 'B'],
    outdir=outdir_case1,
)
FP.JS_divergence(
    msm_data_file=msm_data_file_case1,
    traj_type_col='traj_type',
    traj_type_list=['A', 'B'],
    meta_set_file=meta_set_file,
    outdir=outdir_case1,
)

# ── Run Case 2 — Random split (negative control) ────────────────────────────────────────
folding_pathways_case2 = FP.post_trans(
    msm_data_file=msm_data_file_case2,
    traj_type_col='traj_type',
    traj_type_list=['A', 'B'],
    outdir=outdir_case2,
)
FP.JS_divergence(
    msm_data_file=msm_data_file_case2,
    traj_type_col='traj_type',
    traj_type_list=['A', 'B'],
    meta_set_file=meta_set_file,
    outdir=outdir_case2,
)
```

| JS divergence | Interpretation |
|---------------|----------------|
| Near 0 | A and B explore similar state distributions |
| Near 1 | A and B have divergent state usage |

**Expected outcome:**
- Case 1 should show **elevated JS divergence** — the native-like (A) and misfolded (B) populations traverse metastable states differently.
- Case 2 should show **JS divergence near 0** throughout — random labels produce no systematic separation.

### Expected outputs

Each case produces the same file set in its respective output directory:

| File | Contents |
|------|----------|
| `FoldingPathways_metastablestate_A-B.csv` | Per-type folding pathway probabilities |
| `JS_div_metastablestate_A-B.dat` | Windowed JS divergence time series |

---

## Minimal Workflow – 6: Folding pathway analysis and JS divergence

The `scripts/run_Foldingpathway.py` script computes both folding pathways and JS divergence in a single call. It requires the MSM mapping CSV to already contain the trajectory-type column (see Step 6). Run it separately for each case.

You can provide parameters either directly as CLI flags, via a `--config` file, or both.
When both are provided, **CLI flags override config values** for the same parameter.

```bash
# ── Activation ──────────────────────────────────────────────────────────────────
source ~/.bashrc
conda activate entdetect

# ── Config files ─────────────────────────────────────────────────────────────────
CFG_CASE1=scripts/configs/workflow2_FoldingPathway_case1_config.json
CFG_CASE2=scripts/configs/workflow2_FoldingPathway_case2_config.json

# ── Run ───────────────────────────────────────────────────────────────
# Case 1 — biologically-informed split (Q >= 0.8 and G <= 0.05)
python scripts/run_Foldingpathway.py --config $CFG_CASE1

# Case 2 — random split (negative control)
python scripts/run_Foldingpathway.py --config $CFG_CASE2
```

Config file example 1 (matches `scripts/configs/workflow2_FoldingPathway_case1_config.json`):

```json
{
  "msm_data_file": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/MSM/1ZMR_prod_MSMmapping_QG_native.csv",
  "meta_set_file": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/MSM/1ZMR_prod_meta_set.csv",
  "outdir": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/FoldingPathway_QG_native",
  "traj_type_col": "traj_type",
  "traj_type_list": ["A", "B"],
  "rm_traj_list": [65, 75, 155, 162, 199, 231, 264, 286, 296, 314, 354, 417, 448, 472, 473, 474, 577, 579, 591, 703, 704, 732, 758, 812, 833, 870, 876, 944, 967],
  "log_level": "INFO",
  "logdir": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/FoldingPathway_QG_native/logs"
}
```

Config file example 2 (matches `scripts/configs/workflow2_FoldingPathway_case2_config.json`):

```json
{
  "msm_data_file": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/MSM/1ZMR_prod_MSMmapping_random.csv",
  "meta_set_file": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/MSM/1ZMR_prod_meta_set.csv",
  "outdir": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/FoldingPathway_random",
  "traj_type_col": "traj_type",
  "traj_type_list": ["A", "B"],
  "rm_traj_list": [65, 75, 155, 162, 199, 231, 264, 286, 296, 314, 354, 417, 448, 472, 473, 474, 577, 579, 591, 703, 704, 732, 758, 812, 833, 870, 876, 944, 967],
  "log_level": "INFO",
  "logdir": "/scratch/ims86/EntDetect_Datastore/outputs/workflow2/FoldingPathway_random/logs"
}
```

The JSON/YAML config keys and their matching CLI flags are listed below. For the command-line wrapper, use the same key name with a `--` prefix; the only extra wrapper-only flag is `--config`.

| Term | Definition |
|------|------------|
| `msm_data_file` (`--msm_data_file`) | Input MSM mapping CSV produced by Step 5b, containing frame-level state assignments plus the trajectory-type label column used for pathway and divergence analysis. |
| `meta_set_file` (`--meta_set_file`) | Metastable-state composition CSV produced by the MSM workflow. Required when `state_type` is `microstate`, and retained as a standard input in the wrapper configuration. |
| `outdir` (`--outdir`) | Output directory where folding pathway summaries and JS-divergence results are written. |
| `traj_type_col` (`--traj_type_col`) | Column name in `msm_data_file` that stores trajectory population labels (for example `traj_type` with values `A` and `B`). |
| `traj_type_list` (`--traj_type_list`) | Ordered list of trajectory-type labels to compare in the folding-pathway and JS-divergence analyses. |
| `rm_traj_list` (`--rm_traj_list`) | Optional list of trajectory IDs to exclude before analysis (for example known mirror artifacts). |
| `n_window` (`--n_window`) | Window length used by the JS-divergence/state-occupancy calculation for handling trajectory-end state persistence and smoothing-related behavior. |
| `n_traj` (`--n_traj`) | Total number of trajectories assumed for the analysis context. This is stored in the analysis object even though the current metastable-state workflow derives active counts from the input data. |
| `state_type` (`--state_type`) | State representation to analyze: `metastablestate` or `microstate`. |
| `log_level` (`--log_level`) | Logging verbosity for the wrapper: `DEBUG`, `INFO`, `WARNING`, or `ERROR`. |
| `logdir` (`--logdir`) | Directory where the `run_Foldingpathway` log file is written. If omitted, logging defaults to `outdir`. |

## I/O Reference for run_Foldingpathway.py

Each file below is listed once, followed by column-level details for tabular outputs.

### `$OUTDIR/MSM/<ID>_MSMmapping_<case>.csv`

| I/O | File | File Description |
|---|---|---|
| Input | `$OUTDIR/MSM/<ID>_MSMmapping_<case>.csv` | Labeled frame-level MSM mapping used to compute population-specific folding pathways. |

### `$OUTDIR/MSM/<ID>_meta_set.csv`

| I/O | File | File Description |
|---|---|---|
| Input | `$OUTDIR/MSM/<ID>_meta_set.csv` | Metastable-state composition file used to aggregate transitions and pathway probabilities. |

### `$OUTDIR/FoldingPathways_metastablestate_A-B.csv`

| Column Name | Column Description |
|---|---|
| Time | Time coordinate (or window center) for the pathway estimate. |
| traj_type | Trajectory-population label (for example A or B). |
| pathway | Encoded metastable-state transition pathway representation. |
| probability | Estimated probability for the pathway at the given time/window. |

### `$OUTDIR/JS_div_metastablestate_A-B.dat`

| Column Name | Column Description |
|---|---|
| Time | Time coordinate (or window center) for the divergence calculation. |
| JS_divergence | Jensen-Shannon divergence between trajectory populations at that time/window. |

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `G()` hangs or is very slow | Normal for long trajectories; CG + Topoly is expensive | Reduce `nproc` if memory-limited; use a shorter demo trajectory |
| `SASA()` returns NaN for some frames | Bad backmapped AA structure or topology mismatch | Inspect suspect frames in VMD; these frames are filtered in Workflow 3 |
| Jwalk error: PDB not found | `pdb_file` path incorrect or file missing | Verify `1zmr_model_clean.pdb` is in `$DATASTORE/user_input/reference_structures/` |
| Jwalk error: `freesasa` not found | Package not installed in env | `pip install freesasa` inside the `entdetect` conda env |
| `Combined_GE/` is empty | Full OP run not yet complete | Use pre-computed results from `$DATASTORE/outputs/workflow2/OP/G/Combined_GE/` |
| MSM produces fewer states than `n_large_states` | Empty states discarded | Normal; try a higher `n_large_states` |
| `secondary_struc_defs.txt` not found | File not in `$DATASTORE/user_input/reference_structures/` | Verify files were rsynced to DATASTORE |
| `domain_def.dat` not found | File not in `$DATASTORE/user_input/reference_structures/` | Verify files were rsynced to DATASTORE |

---

← [Workflow 1](workflow1_native_ncle.md) | [Back to Master Index](index.md) | Next → [Workflow 3: Sim-to-Experiment](workflow3_sim2exp.md)

