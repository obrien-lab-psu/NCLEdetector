# Workflow 1: Identify Native Non-Covalent Lasso Entanglements (NCLEs)

← [Back to Master Index](index.md)

---

## Goal

Identify all native NCLEs present in a reference protein structure, filter for high-confidence entanglements, cluster redundant variants, and compute structural features for the representative set.

---

## Table of Contents

- [Step 0. Activate your environment and set paths](#step-0-activate-your-environment-and-set-paths)
- [Step 1. Prepare a cleaned structure](#step-1-prepare-a-cleaned-structure)
- [Step 2. Detect native non-covalent lasso entanglements (NCLEs)](#step-2-detect-native-non-covalent-lasso-entanglements-ncles)
- [Step 3. Filter for high-quality entanglements](#step-3-filter-for-high-quality-entanglements)
- [Step 4. Cluster redundant NCLEs into representative entanglements](#step-4-cluster-redundant-ncles-into-representative-entanglements)
- [Step 5. Compute structural features for representative NCLEs](#step-5-compute-structural-features-for-representative-ncles)
- [Step 6. Visualize the representative NCLEs](#step-6-visualize-the-representative-ncles)
- [Minimal Workflow – 1: Running the NCLE identification pipeline as single script](#minimal-workflow--1-running-the-ncle-identification-pipeline-as-single-script)

---

## Typical runtime

| Step | Runtime |
|------|---------|
| Native entanglement detection | 1–10 minutes |
| High-quality filtering | < 1 minute |
| Clustering | < 1 minute |
| Feature generation | 1–5 minutes |

---

Tutorial outputs from re-running the steps below are written to:
```bash
DATASTORE=/scratch/ims86/EntDetect_Datastore
OUTDIR=$DATASTORE/outputs/workflow1
```

---

## Required input files
| Cleaned all-atom PDB | `$DATASTORE/user_input/reference_structures/1zmr_model_clean.pdb` | Input structure |
| Cα PSF topology | `$DATASTORE/user_input/reference_structures/1zmr_model_clean_ca.psf` | Topology reference |
| Cα COR coordinates | `$DATASTORE/user_input/reference_structures/1zmr_model_clean_ca.cor` | Coordinate reference |
native_HQ_outdir = f"{OUTDIR}/HQ_NCLE"

> *Note*: This tutorial shows you how to handel both all-atom structures but it can handle alpha carbon coarse-grained structures as well.  
  
---

## Step 0. Activate your environment and set paths

```bash
source ~/.bashrc
conda activate entdetect
```

---

## Step 1. Prepare a cleaned structure

The `1zmr_model_clean.pdb` file in `$DATASTORE/user_input/reference_structures/` is already cleaned and ready to use. For your own protein you would:

- rebuild missing residues and atoms (e.g., with **Modeller** or **CHARMM**);
- remove duplicate residue records;
- verify chain labels and residue numbering are sensible;
- ensure the PDB contains only the protein chain(s) of interest.

---

## Step 2. Detect native non-covalent lasso entanglements (NCLEs)

Create a Python script `run_nativeNCLE.py` or run interactively:

> *Note*: There is a complete version of this python script ready to use located at the end of the tutorial: [Jump to full workflow script section](#minimal-workflow--1-running-the-ncle-identification-pipeline-as-single-script)  

```python
from EntDetect.gaussian_entanglement import GaussianEntanglement
from EntDetect.clustering import ClusterNativeEntanglements
from EntDetect.entanglement_features import FeatureGen

# ── Paths ──────────────────────────────────────────────────────────────────
DATASTORE = "/Path/to/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow1"
NCLE_outdir  = f"{OUTDIR}/NCLE"
HQ_NCLE_outdir = f"{OUTDIR}/HQ_NCLE"
clustered_HQ_NCLE_outdir = f"{OUTDIR}/clustered_HQ_NCLE"
clustered_HQ_NCLE_feature_outdir = f"{OUTDIR}/clustered_HQ_NCLE_features"

# ── Inputs ──────────────────────────────────────────────────────────────────
pdb_file            = f"{DATASTORE}/user_input/reference_structures/1zmr_model_clean.pdb"
chain          = "A"
ID             = "1ZMR"
gene           = "P00558"

# ── Initialize and Run ───────────────────────────────────────────────────────────────
# g_threshold         : Gaussian entanglement score cutoff (0.6 is the standard value)
# Calpha              : set True only if user wants to use the Calpha distance to define contacts. False will use heavy atoms between residues. 
# CG                  : set True only if input is a coarse-grained model
# ent_detection_method: 2 = any nonzero TLN for either termini  ← tutorial default (matches class default)
ge = GaussianEntanglement(g_threshold=0.6, density=1.0, Calpha=False, CG=False,
                          ent_detection_method=2)
ge.calculate_native_entanglements(pdb_file=pdb_file, outdir=NCLE_outdir, ID=ID, chain=chain)
```

### What this step does

Scans the input structure and identifies all **native non-covalent lasso entanglements** present in the reference conformation. Each NCLE is characterized by a loop closed by a disulfide or non-covalent contact and one or more chain termini that thread through it.

### Expected output

A CSV file in `NCLE_outdir` named `{ID}_GE.csv`:

```
$OUTDIR/NCLE/1ZMR_GE.csv
```

### Practical note

- If your structure is all-atom (conventional PDB), keep `Calpha=False` and `CG=False`.
- If your structure is all-atom and you wish to use the alpha carbons to define contacts rather than the residue heavy atoms use `Calpha=True` and `CG=False`.
- If the structure contains only Cα atoms (e.g., a alpha carbon coarse-grained model), set `Calpha=True` and `CG=True`.
- The `g_threshold=0.6` is the standard cutoff used in the published protocol. Lower values will detect more (but potentially noisier) entanglements.

---

## Step 3. Filter for high-quality entanglements

This step is especially important when using **AlphaFold** or other predicted structures, or when you want to remove slipknots and low-confidence entanglements from the set.

```python
# ── Inputs from Step 3 ──────────────────────────────────────────────────────────────────
rawNCLE_file = f"{OUTDIR}/NCLE/1ZMR_GE.csv"

# ── Initialize and Run ───────────────────────────────────────────────────────────────
ge = GaussianEntanglement(g_threshold=0.6, density=1.0, Calpha=False, CG=False,
                          ent_detection_method=2)
ge.select_high_quality_entanglements(
    rawNCLE_file=rawNCLE_file,
    pdb_file=pdb_file,
    outdir=HQ_NCLE_outdir,
    ID=ID,
    model="EXP",  # "EXP" for experimental structures; "AF" for AlphaFold models
    chain=chain)
```

### `model` argument

| Value | Use when |
|-------|----------|
| `"EXP"` | Input structure is from X-ray crystallography, cryo-EM, or NMR |
| `"AF"` | Input structure is an AlphaFold prediction |

### Expected output

A CSV file of filtered high-quality entanglements:

```
$OUTDIR/HQ_NCLE/1ZMR.csv
```

> *Note*: (non-canonical regions) If the structure contains residues that are not part of the canonical protein sequence (e.g., expression tags, fusion segments, or engineered regions), a **mapping file** can be supplied to restrict the filtered set to the canonical region only. See the module-level documentation for `select_high_quality_entanglements` for the expected format.

---

## Step 4. Cluster redundant NCLEs into representative entanglements

Many entanglements are structurally degenerate variants of the same topological motif. This step collapses them into a non-redundant representative set.

```python
# ── Inputs from Step 4 ──────────────────────────────────────────────────────────────────
HQ_NCLE_file = f"{OUTDIR}/HQ_NCLE/1ZMR.csv"

# ── Initialize and Run ───────────────────────────────────────────────────────────────
clustering = ClusterNativeEntanglements(organism="Ecoli", cut_off=None)
clustering.Cluster_NativeEntanglements(
    HQ_NCLE_file=HQ_NCLE_file,
    outdir=clustered_HQ_NCLE_outdir,
    ID="1ZMR",
    chain=chain)
```

> *Note*: The `organism` parameter sets the distance threshold used in clustering to that optimized on one of three model orgnaisms. Use `"Ecoli"` for *E. coli* proteins, `"Human"` for *H. sapiens* proteins, and `"Yeast"` for *S. cerevessa* proteins. If you are analyzing a protein outside these model organisms or wish to adjust the cut-off you can specify the `cut_off` parameter which will overide the organism defaults.  

### Expected output

```
$OUTDIR/clustered_HQ_NCLE/1ZMR.csv
```

This CSV contains:

- **representative entanglements** (one row per unique NCLE topology);
- associated **degenerate loop** variants;
- a reduced set suitable for downstream feature generation and interpretation.

---

## Step 5. Compute structural features for representative NCLEs

Once you have the representative set, calculate structural features that characterize each entanglement's complexity.

```python
# ── Inputs from Step 5 ──────────────────────────────────────────────────────────────────
cluster_file = f"{OUTDIR}/clustered_HQ_NCLE/1ZMR.csv"

# ── Initialize and Run ───────────────────────────────────────────────────────────────
FGen = FeatureGen()
# gene       : UniProt accession for the protein
# chain      : PDB chain identifier
# ID         : structure identifier used throughout the workflow
EntFeatures = FGen.get_uent_features(
    pdb_file=pdb_file,
    outdir=clustered_HQ_NCLE_feature_outdir,
    cluster_file=cluster_file,
    gene=gene,
    chain=chain,
    pdbid=ID)
```

### What this step computes

$OUTDIR/clustered_HQ_NCLE_features/P00558_1ZMR_A_uent_features.csv

- supercoiling characteristics of the entangled terminus;
- contact patterns between the loop and the rest of the structure;
- other structural descriptors capturing topological complexity.

### Expected output

```
$OUTDIR/clustered_HQ_NCLE_features/P00558_1ZMR_A_uent_features.csv
```

---

## Step 6. Visualize the representative NCLEs

Use **VMD** or **PyMOL** to inspect the representative NCLEs visually.

### Why visualization matters

- Validates that identified entanglements make structural sense.
- Reveals the loop/crossing geometry in 3D.
- Essential for preparing publication-quality figures.
- Can catch preprocessing errors (e.g., wrong chain, renumbered residues) that tables will not expose.

---

## Minimal Workflow – 1: Running the NCLE identification pipeline as single script

All four steps above are handled by `scripts/run_nativeNCLE.py`. Run it from the repo root.

You can now provide parameters either directly as CLI flags, via a `--config` file, or both.
When both are provided, **CLI flags override config values** for the same parameter.

Example with direct CLI flags:

```bash
source ~/.bashrc
conda activate entdetect

DATASTORE=/Path/to/EntDetect_Datastore

python scripts/run_nativeNCLE.py \
    --pdb_file  $DATASTORE/user_input/reference_structures/1zmr_model_clean.pdb \
    --outdir  $DATASTORE/outputs/workflow1 \
    --ID      1ZMR \
    --chain   A \
    --organism Ecoli \
    --gene    P00558 \
    --model   EXP
```

Example with config plus CLI override:

```bash
CONFIG=scripts/configs/workflow1_nativeNCLE_config.json

# Here, --ent_detection_method 2 overrides ent_detection_method=3 from the config file.
python scripts/run_nativeNCLE.py \
        --config $CONFIG \
    --ent_detection_method 2
```

Container equivalent (same config and override):

```bash
CONFIG=scripts/configs/workflow1_nativeNCLE_config.json
DATASTORE=/Path/to/EntDetect_Datastore

# Here, --ent_detection_method 2 overrides ent_detection_method=3 from the config file.
apptainer exec \
    --bind "$DATASTORE:$DATASTORE" \
    --bind "$PWD:$PWD" \
    --pwd "$PWD" \
    entdetect-latest.sif \
    python scripts/run_nativeNCLE.py \
        --config "$CONFIG" \
        --ent_detection_method 2
```

Config file example (matches `scripts/configs/workflow1_nativeNCLE_config.json`):

```json
{
  "pdb_file": "/scratch/ims86/EntDetect_Datastore/user_input/reference_structures/1zmr_model_clean.pdb",
  "outdir": "/scratch/ims86/EntDetect_Datastore/outputs/tmp/workflow1",
  "ID": "1ZMR",
  "chain": "A",
  "organism": "Ecoli",
  "gene": "P00558",
  "model": "EXP",
  "CG": false,
  "Calpha": false,
  "ent_detection_method": 2,
  "log_level": "INFO",
  "g_threshold": 0.6,
  "density": 1.0,
  "cut_off": null
}
```

The JSON/YAML config keys and their matching CLI flags are listed below. For the command-line wrapper, use the same key name with a `--` prefix; the only extra wrapper-only flag is `--config`.

| Term | Definition |
|------|------------|
| `pdb_file` (`--pdb_file`) | Absolute path to the PDB file that was cleaned in Step 1 and is ready to be analyzed for NCLEs. |
| `outdir` (`--outdir`) | Absolute path to the output directory where workflow files will be written. |
| `ID` (`--ID`) | String identifier for the structure that is appended to output filenames. |
| `chain` (`--chain`) | Molecule chain in the input structure; omit it to process all chains. |
| `organism` (`--organism`) | Source organism used to choose the default clustering threshold: `Human`, `Yeast`, or `Ecoli`. |
| `gene` (`--gene`) | UniProt accession for the protein's source gene, used in feature-file naming. |
| `model` (`--model`) | Whether the structure is experimental (`EXP`) or predicted (`AF`). |
| `CG` (`--CG`) | Whether to treat the input as a coarse-grained C-alpha model with one bead per residue. |
| `Calpha` (`--Calpha`) | Whether to use alpha carbons, rather than heavy sidechain atoms, to define native non-covalent contacts. |
| `ent_detection_method` (`--ent_detection_method`) | NCLE detection criterion: `1` for any nonzero GLN, `2` for any nonzero TLN, `3` for both GLN and TLN nonzero for the same terminus. |
| `log_level` (`--log_level`) | Logging verbosity for the single-run workflow: `DEBUG`, `INFO`, `WARNING`, or `ERROR`. |
| `g_threshold` (`--g_threshold`) | Threshold used when rounding the discrete Gaussian linking integral into GLN values. |
| `density` (`--density`) | Topoly triangulation density used when building loop surfaces for piercing detection. |
| `cut_off` (`--cut_off`) | Optional custom clustering threshold in Angstroms; use this to override the organism-specific default. |

---

## I/O Reference for run_nativeNCLE.py

Each file below is listed once, followed by the column-level description table for that file.

### `$DATASTORE/user_input/reference_structures/1zmr_model_clean.pdb`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/user_input/reference_structures/1zmr_model_clean.pdb` | Cleaned all-atom structure used as the geometric source for NCLE discovery. |

### `$OUTDIR/NCLE/1ZMR_GE.csv`

| Column Name | Column Description |
|---|---|
| ID | Protein or structure identifier for the entanglement record. |
| chain | Chain identifier in the input structure. |
| i | First residue index of the native contact that defines the loop start. |
| j | Second residue index of the native contact that defines the loop end. |
| crossingsN | Comma-separated N-terminal crossing residues, stored with chirality signs such as +109 or -258. |
| crossingsC | Comma-separated C-terminal crossing residues, stored with chirality signs such as +109 or -258. |
| gn | Raw Gaussian entanglement score for the N-terminal side of the loop. |
| gc | Raw Gaussian entanglement score for the C-terminal side of the loop. |
| GLNn | Rounded Gaussian linking number for the N-terminal side, derived from gn. |
| GLNc | Rounded Gaussian linking number for the C-terminal side, derived from gc. |
| TLNn | Topological linking number for the N-terminal side after Topoly-based crossing assignment and buffer filtering. |
| TLNc | Topological linking number for the C-terminal side after Topoly-based crossing assignment and buffer filtering. |
| CCbond | Boolean flag indicating that the loop contact is a covalent C-C bond, typically a disulfide-linked loop. |
| ENT | Boolean flag showing whether the native contact passed the chosen entanglement detection criterion. |

### `$OUTDIR/HQ_NCLE/1ZMR.csv`

| Column Name | Column Description |
|---|---|
| ID | Protein or structure identifier for the filtered entanglement record. |
| chain | Chain identifier in the input structure. |
| i | First residue index of the loop-defining contact. |
| j | Second residue index of the loop-defining contact. |
| crossingsN | N-terminal crossing residues retained after slipknot filtering. |
| crossingsC | C-terminal crossing residues retained after slipknot filtering. |
| gn | Raw Gaussian entanglement score for the N-terminal side. |
| gc | Raw Gaussian entanglement score for the C-terminal side. |
| GLNn | Rounded Gaussian linking number for the N-terminal side. |
| GLNc | Rounded Gaussian linking number for the C-terminal side. |
| TLNn | Topological linking number for the N-terminal side. |
| TLNc | Topological linking number for the C-terminal side. |
| CCbond | Boolean flag indicating a covalent C-C loop closure. |
| ENT | Boolean flag showing whether the record is still classified as an entanglement after quality filtering. |
| Slipknot_N | Boolean flag indicating that the N-terminal crossing pattern is consistent with a slipknot and should be treated cautiously. |
| Slipknot_C | Boolean flag indicating that the C-terminal crossing pattern is consistent with a slipknot and should be treated cautiously. |

### `$OUTDIR/clustered_HQ_NCLE/1ZMR.csv`

| Column Name | Column Description |
|---|---|
| ID | Protein or structure identifier for the clustered representative. |
| chain | Chain identifier in the input structure. |
| i | First residue index of the representative loop contact. |
| j | Second residue index of the representative loop contact. |
| crossingsN | Representative N-terminal crossing residues for the cluster. |
| crossingsC | Representative C-terminal crossing residues for the cluster. |
| gn | Representative raw Gaussian entanglement score for the N-terminal side. |
| gc | Representative raw Gaussian entanglement score for the C-terminal side. |
| GLNn | Representative rounded Gaussian linking number for the N-terminal side. |
| GLNc | Representative rounded Gaussian linking number for the C-terminal side. |
| TLNn | Representative topological linking number for the N-terminal side. |
| TLNc | Representative topological linking number for the C-terminal side. |
| num_contacts | Number of raw entanglement records merged into the representative cluster. |
| contacts | Semicolon-separated list of the loop residue pairs that were merged into the representative cluster. |
| CCBond | Boolean flag indicating whether the representative loop closure is a covalent C-C bond. |

### `$OUTDIR/clustered_HQ_NCLE_features/P00558_1ZMR_A_uent_features.csv`

| Column Name | Column Description |
|---|---|
| gene | UniProt accession used for protein-level annotation and output naming. |
| PDB | Four-letter PDB identifier for the structure. |
| chain | Chain identifier used for feature extraction. |
| ENT-ID | Row index of the representative entanglement within the clustered input file. |
| gn | Representative raw Gaussian entanglement score for the N-terminal side. |
| N_term_thread | Number of N-terminal crossing residues threading through the loop. |
| gc | Representative raw Gaussian entanglement score for the C-terminal side. |
| C_term_thread | Number of C-terminal crossing residues threading through the loop. |
| i | First residue index of the representative loop contact. |
| j | Second residue index of the representative loop contact. |
| NC | Core loop contact residues, stored as the pair i,j. |
| NC_wbuff | Loop residues expanded by the fixed residue buffer used for neighborhood calculations. |
| NC_region | Residues whose side-chain heavy atoms lie within the neighborhood cutoff of the buffered loop residues. |
| crossingsN | Core N-terminal crossing residues, without the buffer expansion. |
| crossingsC | Core C-terminal crossing residues, without the buffer expansion. |
| crossingsN_wbuff | N-terminal crossing residues expanded by the fixed residue buffer. |
| crossingsC_wbuff | C-terminal crossing residues expanded by the fixed residue buffer. |
| crossingsN_region | Residues whose side-chain heavy atoms lie within the neighborhood cutoff of the buffered N-terminal crossings. |
| crossingsC_region | Residues whose side-chain heavy atoms lie within the neighborhood cutoff of the buffered C-terminal crossings. |
| ent_region | Union of the loop neighborhood, crossing neighborhoods, core loop residues, and core crossing residues. |
| loopsize | Loop length in residues, computed as j minus i. |
| num_zipper_nc | Number of raw entanglement contacts that were merged into the representative cluster. |
| perc_bb_loop | Loop length normalized by the full protein length. |
| num_loop_contacting_res | Number of residues whose side-chain heavy atoms contact the loop neighborhood. |
| num_cross_nearest_neighbors | Total number of residues near the N- and C-terminal crossing neighborhoods. |
| ent_coverage | Fraction of the protein covered by the full entangled region. |
| min_N_prot_depth_left | Fractional distance of the closest N-terminal crossing from the protein N-terminus. |
| min_N_thread_depth_left | Fractional distance of the closest N-terminal crossing from the loop start. |
| min_N_thread_slippage_left | Residue offset of the closest N-terminal crossing from the loop start. |
| min_C_prot_depth_right | Fractional distance of the closest C-terminal crossing from the protein C-terminus. |
| min_C_thread_depth_right | Fractional distance of the closest C-terminal crossing from the loop end. |
| min_C_thread_slippage_right | Residue offset of the closest C-terminal crossing from the protein C-terminus. |
| prot_size | Full protein length from the UniProt canonical sequence. |
| ACO | Average contact order of the merged loop contacts in the representative cluster. |
| RCO | Relative contact order, computed as ACO divided by prot_size. |
| CCBond | Boolean flag indicating whether the representative loop closure is a covalent C-C bond. |

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `FileNotFoundError` on PDB | Wrong path | Verify the file exists |
| Empty output file from Step 3 | Structure has no detectable NCLEs, or `g_threshold` too high | Lower `g_threshold` to 0.4 and rerun |
| Step 6 raises `KeyError` on gene or chain | `gene`, `chain`, or `pdbid` arguments do not match the structure | Check the PDB header for the correct chain ID |
| `ClusterNativeEntanglements` returns only 1 cluster | Normal for small proteins or proteins with a single entanglement topology | Proceed; the representative is the single detected NCLE |

---

← [Back to Master Index](index.md) | Next → [Workflow 2: Trajectory Analysis](workflow2_trajectory_analysis.md)
