# Workflow 1: Identify Native Non-Covalent Lasso Entanglements (NCLEs)

← [Back to Master Index](index.md)

---

## Goal

Identify all native NCLEs present in a reference protein structure, filter for high-confidence entanglements, cluster redundant variants, and compute structural features for the representative set.

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

| File | Path | Role |
|------|------|------|
| Cleaned all-atom PDB | `$DATASTORE/user_input/reference_structures/1zmr_model_clean.pdb` | Input structure |
| Cα PSF topology | `$DATASTORE/user_input/reference_structures/1zmr_model_clean_ca.psf` | Topology reference |
| Cα COR coordinates | `$DATASTORE/user_input/reference_structures/1zmr_model_clean_ca.cor` | Coordinate reference |

> *Note*: This tutorial shows you how to handel both all-atom structures but it can handle alpha carbon coarse-grained structures as well.  
  
---

## Step 1. Activate your environment and set paths

```bash
source ~/.bashrc
conda activate entdetect

DATASTORE=/scratch/ims86/EntDetect_Datastore
OUTDIR=$DATASTORE/outputs/workflow1

mkdir -p $OUTDIR/Native_GE $OUTDIR/Native_HQ_GE $OUTDIR/Native_clustered_HQ_GE $OUTDIR/Native_clustered_HQ_GE_features
```

---

## Step 2. Prepare a cleaned structure

The `1zmr_model_clean.pdb` file in `$DATASTORE/user_input/reference_structures/` is already cleaned and ready to use. For your own protein you would:

- rebuild missing residues and atoms (e.g., with **Modeller** or **CHARMM**);
- remove duplicate residue records;
- verify chain labels and residue numbering are sensible;
- ensure the PDB contains only the protein chain(s) of interest.

---

## Step 3. Detect native entanglements

Create a Python script `run_nativeNCLE.py` or run interactively:

> *Note*: There is a complete version of this python script ready to use located at the end of the tutorial: [Jump to full workflow script section](#running-the-full-workflow-as-a-single-script)  

```python
from EntDetect.gaussian_entanglement import GaussianEntanglement

# ── Inputs ──────────────────────────────────────────────────────────────────
DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow1"

pdb            = f"{DATASTORE}/user_input/reference_structures/1zmr_model_clean.pdb"
native_outdir  = f"{OUTDIR}/Native_GE"
ID             = "1ZMR"

# ── Initialize ───────────────────────────────────────────────────────────────
# g_threshold         : Gaussian entanglement score cutoff (0.6 is the standard value)
# density             : contact density filter (0.0 = no filter)
# Calpha              : set True only if user wants to use the Calpha distance to define contacts. False will use heavy atoms between residues. 
# CG                  : set True only if input is a coarse-grained model
# ent_detection_method: 2 = any nonzero TLN for either termini  ← tutorial default (matches class default)
#                       3 = both GLN and TLN nonzero for same termini  ← recommended for production
ge = GaussianEntanglement(g_threshold=0.6, density=0.0, Calpha=False, CG=False,
                          ent_detection_method=2)

# ── Run ──────────────────────────────────────────────────────────────────────
ge.calculate_native_entanglements(pdb_file=pdb, outdir=native_outdir, ID=ID, chain='A')
```

### What this step does

Scans the input structure and identifies all **native non-covalent lasso entanglements** present in the reference conformation. Each NCLE is characterized by a loop closed by a disulfide or non-covalent contact and one or more chain termini that thread through it.

### Expected output

A CSV file in `native_outdir` named `{ID}_GE.csv`:

```
$OUTDIR/Native_GE/1ZMR_GE.csv
```

### Practical note

- If your structure is all-atom (conventional PDB), keep `Calpha=False` and `CG=False`.
- If your structure is all-atom and you wish to use the alpha carbons to define contacts rather than the residue heavy atoms use `Calpha=True` and `CG=False`.
- If the structure contains only Cα atoms (e.g., a alpha carbon coarse-grained model), set `Calpha=True` and `CG=True`.
- The `g_threshold=0.6` is the standard cutoff used in the published protocol. Lower values will detect more (but potentially noisier) entanglements.

---

## Step 4. Filter for high-quality entanglements

This step is especially important when using **AlphaFold** or other predicted structures, or when you want to remove slipknots and low-confidence entanglements from the set.

```python
from EntDetect.gaussian_entanglement import GaussianEntanglement

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow1"

pdb              = f"{DATASTORE}/user_input/reference_structures/1zmr_model_clean.pdb"
native_HQ_outdir = f"{OUTDIR}/Native_HQ_GE"
ID               = "1ZMR"

# Point to the output from Step 3 (or the pre-computed reference)
NCLE_file = f"{OUTDIR}/Native_GE/1ZMR_GE.csv"

ge = GaussianEntanglement(g_threshold=0.6, density=0.0, Calpha=False, CG=False,
                          ent_detection_method=2)

ge.select_high_quality_entanglements(
    NCLE_file,
    pdb,
    outdir=native_HQ_outdir,
    ID=ID,
    model="EXP",  # "EXP" for experimental structures; "AF" for AlphaFold models
    chain='A'
)
```

### `model` argument

| Value | Use when |
|-------|----------|
| `"EXP"` | Input structure is from X-ray crystallography, cryo-EM, or NMR |
| `"AF"` | Input structure is an AlphaFold prediction |

### Expected output

A CSV file of filtered high-quality entanglements:

```
$OUTDIR/Native_HQ_GE/1ZMR.csv
```

> *Note*: (non-canonical regions) If the structure contains residues that are not part of the canonical protein sequence (e.g., expression tags, fusion segments, or engineered regions), a **mapping file** can be supplied to restrict the filtered set to the canonical region only. See the module-level documentation for `select_high_quality_entanglements` for the expected format.

---

## Step 5. Cluster redundant NCLEs into representative entanglements

Many entanglements are structurally degenerate variants of the same topological motif. This step collapses them into a non-redundant representative set.

```python
from EntDetect.clustering import ClusterNativeEntanglements

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow1"

native_clustered_HQ_outdir = f"{OUTDIR}/Native_clustered_HQ_GE"
outfile = "1ZMR.csv"

# Point to Step 4 output (or pre-computed reference)
NCLE_file = f"{OUTDIR}/Native_HQ_GE/1ZMR.csv"

clustering = ClusterNativeEntanglements(organism="Ecoli", cut_off=None)
clustering.Cluster_NativeEntanglements(
    NCLE_file,
    outdir=native_clustered_HQ_outdir,
    outfile=outfile,
    chain='A'
)
```

> *Note*: The `organism` parameter sets the distance threshold used in clustering to that optimized on one of three model orgnaisms. Use `"Ecoli"` for *E. coli* proteins, `"Human"` for *H. sapiens* proteins, and `"Yeast"` for *S. cerevessa* proteins. If you are analyzing a protein outside these model organisms or wish to adjust the cut-off you can specify the `cut_off` parameter which will overide the organism defaults.  

### Expected output

```
$OUTDIR/Native_clustered_HQ_GE/1ZMR.csv
```

This CSV contains:

- **representative entanglements** (one row per unique NCLE topology);
- associated **degenerate loop** variants;
- a reduced set suitable for downstream feature generation and interpretation.

---

## Step 6. Compute structural features for representative NCLEs

Once you have the representative set, calculate structural features that characterize each entanglement's complexity.

```python
from EntDetect.entanglement_features import FeatureGen

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR    = f"{DATASTORE}/outputs/workflow1"

pdb                      = f"{DATASTORE}/user_input/reference_structures/1zmr_model_clean.pdb"
native_GQ_feature_outdir = f"{OUTDIR}/Native_clustered_HQ_GE_features"

# Point to Step 5 output (or pre-computed reference)
cluster_file = f"{OUTDIR}/Native_clustered_HQ_GE/1ZMR.csv"

FGen = FeatureGen(pdb, outdir=native_GQ_feature_outdir, cluster_file=cluster_file)

# gene    : UniProt accession for ecPGK
# chain   : PDB chain identifier
# pdbid   : four-letter PDB code (uppercase)
EntFeatures = FGen.get_uent_features(gene='P00558', chain='A', pdbid='1ZMR')
print(EntFeatures)
```

### What this step computes

For each representative NCLE:

- supercoiling characteristics of the entangled terminus;
- contact patterns between the loop and the rest of the structure;
- other structural descriptors capturing topological complexity.

### Expected output

```
$OUTDIR/Native_clustered_HQ_GE_features/P00558_1ZMR_A_uent_features.csv
```

---

## Step 7. Visualize the representative entanglements

Use **VMD** or **PyMOL** to inspect the representative NCLEs visually.

### Why visualization matters

- Validates that identified entanglements make structural sense.
- Reveals the loop/crossing geometry in 3D.
- Essential for preparing publication-quality figures.
- Can catch preprocessing errors (e.g., wrong chain, renumbered residues) that tables will not expose.

---

## Running the full workflow as a single script

All four steps above are handled by `scripts/run_nativeNCLE.py`. Run it from the repo root:

```bash
source ~/.bashrc
conda activate entdetect

DATASTORE=/scratch/ims86/EntDetect_Datastore

python scripts/run_nativeNCLE.py \
    --struct  $DATASTORE/user_input/reference_structures/1zmr_model_clean.pdb \
    --outdir  $DATASTORE/outputs/workflow1 \
    --ID      1ZMR \
    --chain   A \
    --organism Ecoli \
    --Accession P00558 \
    --model   EXP
```

| Flag | Value | Notes |
|------|-------|-------|
| `--struct` | path to cleaned PDB | All-atom experimental structure |
| `--outdir` | output root directory | Sub-dirs are created automatically |
| `--ID` | four-letter PDB code (uppercase) | Used for output filenames |
| `--chain` | `A` | Omit to process all chains |
| `--organism` | `Ecoli` | Reference proteome for clustering |
| `--Accession` | `P00558` | UniProt accession for feature generation |
| `--model` | `EXP` | `EXP` for crystal/cryo-EM; `AF` for AlphaFold |

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
