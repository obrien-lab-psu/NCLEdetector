# Workflow 4: Population-Level Detection of Misfolding Involving Native Entanglements

← [Back to Master Index](index.md)

---

## Goal

Use proteome-scale structure and LiP-MS-derived residue labels to test whether residues in native entangled regions are enriched for experimental conformational-change signals, then identify high-risk protein subpopulations by Monte Carlo optimization.

---

## Table of Contents

- [Step 1. Activate your environment and set paths](#step-1-activate-your-environment-and-set-paths)
- [Step 2. Verify proteome structure inputs](#step-2-verify-proteome-structure-inputs)
- [Step 3. Run Workflow 1 in batch to generate native NCLEs per protein](#step-3-run-workflow-1-in-batch-to-generate-native-ncles-per-protein)
- [Step 3c. Single-run batch script for Steps 3 and 4 (recommended)](#step-3c-single-run-batch-script-for-steps-3-and-4-recommended)
- [Step 4. Build residue-level modeling tables](#step-4-build-residue-level-modeling-tables)
- [Step 5. Run proteome-level logistic regression](#step-5-run-proteome-level-logistic-regression)
- [Running proteome regression as a single script](#running-proteome-regression-as-a-single-script)
- [Step 6. Run Monte Carlo subpopulation selection](#step-6-run-monte-carlo-subpopulation-selection)
- [Running Monte Carlo as a single script](#running-monte-carlo-as-a-single-script)
- [Step 7. Rank and refine candidate protein groups](#step-7-rank-and-refine-candidate-protein-groups)

---

## Typical runtime

| Step | Runtime |
|------|---------|
| Batch native NCLE generation (proteome-wide) | Hours to days |
| Residue-table construction | Minutes to hours |
| Logistic regression | Minutes |
| Monte Carlo (single run) | Hours |
| Monte Carlo (multiple independent runs) | Days |

---

## Scope and assumptions

This workflow is population-level and is independent of the single-protein trajectory analysis in Workflows 2 and 3.

It assumes you can provide, for each protein in your study:

1. native entanglement region annotations (from Workflow 1 outputs), and
2. residue-level experimental labels (for example, LiP-MS significant vs non-significant residues).

---

## Available structure libraries

Pre-downloaded proteome structure libraries are available at:

```bash
DATASTORE=/scratch/ims86/EntDetect_Datastore
STRUCTDIR=$DATASTORE/user_input/proteome_structures

$STRUCTDIR/AF    # AlphaFold structures (example filenames: A0A385XJ53.pdb)
$STRUCTDIR/EXP   # Experimental structures (example filenames: A5A618-6RKO_H.pdb)
                 # Optional mapping files: *_resid_mapping.txt
```

Tutorial outputs are written to:

```bash
OUTDIR=$DATASTORE/outputs/workflow4
```

---

## Required input files

| File / directory | Path | Role |
|------------------|------|------|
| AlphaFold proteome structures | `$DATASTORE/user_input/proteome_structures/AF/` | Native structure library (AF mode) |
| Experimental proteome structures | `$DATASTORE/user_input/proteome_structures/EXP/` | Native structure library (EXP mode) |
| Optional residue mapping files | `$DATASTORE/user_input/proteome_structures/EXP/*_resid_mapping.txt` | Map structure residue numbering to canonical indices |
| Gene list file (tutorial test case) | `$DATASTORE/user_input/experimental_data/Gene_lists/AF/AF_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt` | Controls cohort used in modeling |
| AF per-residue features table | `$DATASTORE/user_input/experimental_data/PDB_residue_features/AF/residueFeatures.csv` | Provides response/predictor columns (e.g., `cut_C_Rall`, `AA`) |
| Combined residue design matrix | Built in Step 4 (`$OUTDIR/residue_dataframes_workflow4.csv`) | Input for Step 5 regression |

---

## Step 1. Activate your environment and set paths

```bash
source ~/.bashrc
conda activate entdetect

DATASTORE=/scratch/ims86/EntDetect_Datastore
STRUCTDIR=$DATASTORE/user_input/proteome_structures
OUTDIR=$DATASTORE/outputs/workflow4
GENELIST_AF=$DATASTORE/user_input/experimental_data/Gene_lists/AF/AF_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt
RES_FEATURES_AF=$DATASTORE/user_input/experimental_data/PDB_residue_features/AF/residueFeatures.csv

mkdir -p $OUTDIR/nativeNCLE_all \
         $OUTDIR/population_modeling \
         $OUTDIR/monte_carlo \
         $OUTDIR/logs
```

---

## Step 2. Verify proteome structure inputs

Check both structure libraries and confirm you can enumerate structures before launching batch jobs.

```bash
ls $STRUCTDIR/AF  | head -n 5
ls $STRUCTDIR/EXP | head -n 5
```

Expected pattern examples:

- AF: `A0A385XJ53.pdb`
- EXP: `A5A618-6RKO_H.pdb`
- EXP mapping: `A5A618-6RKO_H_resid_mapping.txt`

---

## Step 3. Run Workflow 1 in batch to generate native NCLEs per protein

For each protein structure, run `scripts/run_nativeNCLE.py` and collect outputs under a per-protein directory.

### 3a. Single-protein example

```bash
python scripts/run_nativeNCLE.py \
    --struct $STRUCTDIR/AF/A0A385XJ53.pdb \
    --outdir $OUTDIR/nativeNCLE_all/A0A385XJ53 \
    --ID A0A385XJ53 \
    --organism Ecoli \
    --model AF \
    --contacts heavy \
    --resolution aa \
    --ent_detection_method 3 \
    --logdir $OUTDIR/logs
```

### 3b. Batch example (AF structures)

```bash
for pdb in $STRUCTDIR/AF/*.pdb; do
    id=$(basename "$pdb" .pdb)
    python scripts/run_nativeNCLE.py \
        --struct "$pdb" \
        --outdir "$OUTDIR/nativeNCLE_all/$id" \
        --ID "$id" \
        --organism Ecoli \
        --model AF \
        --contacts heavy \
        --resolution aa \
        --ent_detection_method 3 \
        --logdir $OUTDIR/logs
done
```

### 3c. Single-run batch script for Steps 3 and 4 (recommended)

Use the dedicated Workflow 4 batch runner to:

1. run native NCLE analysis in parallel for proteins in the gene list, and
2. build the final Step 4 design matrix from per-residue feature inputs plus NCLE-derived `region` labels.

```bash
python scripts/run_workflow4_nativeNCLE_batch.py \
    --pdb_dir $STRUCTDIR/AF \
    --gene_list $GENELIST_AF \
    --outdir $OUTDIR/nativeNCLE_all \
    --organism Ecoli \
    --model AF \
    --contacts heavy \
    --resolution aa \
    --ent_detection_method 3 \
    --residue_features_file $RES_FEATURES_AF \
    --reg_formula "cut_C_Rall ~ AA + region" \
    --design_matrix_file $OUTDIR/residue_dataframes_workflow4.csv \
    --nproc 16 \
    --logdir $OUTDIR/logs
```

Notes:

- This wrapper runs `scripts/run_nativeNCLE.py` in parallel.
- It accepts all core `run_nativeNCLE.py` options except `--struct`, which is set per PDB automatically.
- Selection rule: the PDB stem (filename without `.pdb`) must be present in `--gene_list`.
- Per-protein NCLE outputs are created under `$OUTDIR/nativeNCLE_all/<ID>/`.
- By default, accessions in generated NCLE feature filenames now use each protein ID (no fixed `P00558` prefix).
- Step 4 output is generated automatically as a single combined matrix:
    - `$OUTDIR/residue_dataframes_workflow4.csv`

### Expected per-protein outputs

```
$OUTDIR/nativeNCLE_all/<ID>/
├── Native_GE/
├── Native_HQ_GE/
├── Native_clustered_HQ_GE/
└── Native_clustered_HQ_GE_features/
```

---

## Step 4. Build residue-level modeling table

Create one combined table across proteins. This matrix is the direct input for Step 5 regression.

Recommended path: use Step 3c command above to build these automatically.

### Required columns

| Column | Meaning |
|--------|---------|
| `gene` | Protein identifier matching entries in `gene_list` |
| `mapped_resid` | Canonical residue index used as inclusion mask (non-null rows are used) |
| `uniprot_length` | Full canonical protein length |
| `AA` | Amino-acid single-letter code |
| `region` | Binary: 1 if residue is in native entangled region, else 0 |
| `cut_C_Rall` | Binary response: 1 if residue has significant experimental signal, else 0 |

In this tutorial, `cut_C_Rall` and `AA` come from:

- `$DATASTORE/user_input/experimental_data/PDB_residue_features/AF/residueFeatures.csv`

`region` is computed on-the-fly from NCLE output `ent_region` values in:

- `$OUTDIR/nativeNCLE_all/<ID>/Native_clustered_HQ_GE_features/*_uent_features.csv`

### File format expectations

- File separator: `|`
- Single combined matrix file: `$OUTDIR/residue_dataframes_workflow4.csv`
- Must contain `gene` values matching entries in `gene_list`

### Minimal example row

```text
gene|mapped_resid|uniprot_length|AA|region|cut_C_Rall
A0A385XJ53|145|312|K|1|1
```

> Important: The quality of this workflow is dominated by how you build `cut_C_Rall` and `region` labels. Keep your peptide-to-residue mapping, significance thresholds, and masking rules consistent and documented.

---

## Step 5. Run proteome-level logistic regression

Use `ProteomeLogisticRegression` to test enrichment while controlling for confounders.

### Python API example

```python
from EntDetect.statistics import ProteomeLogisticRegression

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR = f"{DATASTORE}/outputs/workflow4"
GENELIST_AF = f"{DATASTORE}/user_input/experimental_data/Gene_lists/AF/AF_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt"

dataframe_files = f"{OUTDIR}/residue_dataframes_workflow4.csv"
gene_list = GENELIST_AF
reg_outdir = f"{OUTDIR}/population_modeling"
ID = "Ecoli_population"
reg_formula = "cut_C_Rall ~ AA + region"

ProtRegression = ProteomeLogisticRegression(
    dataframe_files=dataframe_files,
    outdir=reg_outdir,
    gene_list=gene_list,
    ID=ID,
    reg_formula=reg_formula,
)

ProtRegression.load_data(
    sep='|',
    reg_var=['AA', 'region'],
    response_var='cut_C_Rall',
    var2binarize=['cut_C_Rall', 'region'],
    mask_column='mapped_resid',
)

reg_df = ProtRegression.run()
print(reg_df)
```

### Interpreting the main term

For `region`, compute odds ratio as:

$$
\text{OR}_{region} = e^{\beta_{region}}
$$

- `OR > 1` with significant p-value: enrichment in entangled regions
- `OR < 1` with significant p-value: depletion in entangled regions
- non-significant p-value: no detectable association in this model

---

## Running proteome regression as a single script

```bash
python scripts/run_population_modeling.py \
    --dataframe_files $OUTDIR/residue_dataframes_workflow4.csv \
    --outdir $OUTDIR/population_modeling \
    --gene_list $GENELIST_AF \
    --tag Ecoli_population \
    --reg_formula "cut_C_Rall ~ AA + region" \
    --logdir $OUTDIR/logs
```

Expected output:

```bash
$OUTDIR/population_modeling/regression_results_Ecoli_population.csv
```

---

## Step 6. Run Monte Carlo subpopulation selection

Use `MonteCarlo` to partition proteins into groups and optimize an objective combining enrichment signal and group-size distribution constraints.

### Python API example

```python
from EntDetect.statistics import MonteCarlo

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR = f"{DATASTORE}/outputs/workflow4"
GENELIST_AF = f"{DATASTORE}/user_input/experimental_data/Gene_lists/AF/AF_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt"

mc_outdir = f"{OUTDIR}/monte_carlo"

MC = MonteCarlo(
    dataframe_files=f"{OUTDIR}/residue_dataframes_workflow4.csv",
    outdir=mc_outdir,
    gene_list=GENELIST_AF,
    ID="Ecoli_population_mc",
    reg_formula="cut_C_Rall ~ region + AA",
    response_var="cut_C_Rall",
    test_var="region",
    random=False,
    n_groups=4,
    steps=100000,
    C1=1.0,
    C2=2.5,
    beta=0.05,
    linearT=False,
)

MC.load_data(
    sep='|',
    reg_var=['AA', 'region'],
    response_var='cut_C_Rall',
    var2binarize=['cut_C_Rall', 'region'],
    mask_column='mapped_resid',
    ID_column='gene',
    Length_column='uniprot_length',
)

MC.run(encoded_df=MC.data, ID_column='gene')
```

---

## Running Monte Carlo as a single script

```bash
python scripts/run_montecarlo.py \
    --dataframe_files $OUTDIR/residue_dataframes_workflow4.csv \
    --outpath $OUTDIR/monte_carlo \
    --gene_list $GENELIST_AF \
    --tag Ecoli_population_mc \
    --reg_formula "cut_C_Rall ~ region + AA" \
    --response_var cut_C_Rall \
    --test_var region \
    --n_groups 4 \
    --steps 100000 \
    --C1 1.0 \
    --C2 2.5 \
    --beta 0.05 \
    --logdir $OUTDIR/logs
```

Typical Monte Carlo outputs include:

- `Final_step_reg_<tag>.csv`
- `State0_final_genelist_<tag>.txt` ... `StateN_final_genelist_<tag>.txt`
- `State0_final_traj_<tag>.csv` ... `StateN_final_traj_<tag>.csv`

---

## Step 7. Rank and refine candidate protein groups

Use the final state files from repeated Monte Carlo runs to build robust candidate sets:

1. Run Monte Carlo multiple times with different random seeds or independent starts.
2. For each run, identify the top-enrichment state by final odds ratio.
3. Count per-protein inclusion frequency in top states across runs.
4. Keep proteins with high reproducibility (for example, >70% top-state inclusion).

This produces a stable, high-confidence subpopulation for downstream biological interpretation.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| No regression rows produced | `gene_list` does not match residue dataframe filenames | Ensure each gene in `gene_list.txt` appears in exactly one dataframe filename |
| Many warnings: "No residue feature file found" | Running older code path expecting per-gene files | Update to latest scripts/SLURM to use the combined matrix workflow |
| Regression fails due to missing columns | Required columns absent in residue tables | Add `gene`, `mapped_resid`, `uniprot_length`, `AA`, `region`, `cut_C_Rall` |
| Monte Carlo crashes on `log(OR)` | OR becomes non-positive from degenerate contingency tables | Increase cohort size, adjust label sparsity, or filter extremely small proteins |
| Monte Carlo outputs not found where expected | `--outpath` missing trailing slash behavior in script outputs | Use a dedicated empty output directory and inspect generated filenames directly |

---

← [Workflow 3](workflow3_sim2exp.md) | [Back to Master Index](index.md) | Next → [Tutorial Examples](tutorial_examples.md)
