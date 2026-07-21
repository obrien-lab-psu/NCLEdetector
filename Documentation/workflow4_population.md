# Workflow 4: Population-Level Detection of Misfolding Involving Native Entanglements

← [Back to Master Index](index.md)

---

## Goal

Use proteome-scale structure and LiP-MS-derived residue labels to test whether residues in native entangled regions are enriched for experimental conformational-change signals, then identify high-risk protein subpopulations by Monte Carlo optimization.

---

## Table of Contents

- [Step 0. Activate your environment and set paths](#step-0-activate-your-environment-and-set-paths)
- [Step 1. Run Workflow 1 in batch to generate native NCLEs per protein](#step-1-run-workflow-1-in-batch-to-generate-native-ncles-per-protein)
- [Step 2. Run proteome-level logistic regression](#step-2-run-proteome-level-logistic-regression)
- [Minimal Workflow – 8: Proteome-wide regression for predicting experimental signals](#minimal-workflow--8-proteome-wide-regression-for-predicting-experimental-signals)
- [Step 3. Run Monte Carlo subpopulation selection](#step-3-run-monte-carlo-subpopulation-selection)
- [Minimal Workflow – 9: Candidate selection by Monte Carlo](#minimal-workflow--9-candidate-selection-by-monte-carlo)
- [Step 4. Rank and refine candidate protein groups](#step-4-rank-and-refine-candidate-protein-groups)

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
DATASTORE=/scratch/ims86/NCLEdetector_Datastore
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
| Combined residue design matrix | Built in Step 1 (`$OUTDIR/residue_dataframes_workflow4.csv`) | Input for Step 2 regression and Step 3 Monte Carlo |

---

## Step 0. Activate your environment and set paths

```bash
source ~/.bashrc
conda activate ncledetector

DATASTORE=/scratch/ims86/NCLEdetector_Datastore
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

## Step 1. Run Workflow 1 in batch to generate native NCLEs per protein

Use the dedicated Workflow 4 batch runner to:

1. run native NCLE analysis in parallel for proteins in the gene list, and
2. build the final residue-level design matrix from per-residue feature inputs plus NCLE-derived `region` labels.

### Example with direct CLI flags:

```bash
# ── Activation ─────────────────────────────────────────────────────────────
source ~/.bashrc
conda activate ncledetector

# ── Paths ──────────────────────────────────────────────────────────────────
DATASTORE=/scratch/ims86/NCLEdetector_Datastore
STRUCTDIR=$DATASTORE/user_input/proteome_structures
OUTDIR=$DATASTORE/outputs/workflow4
GENELIST_AF=$DATASTORE/user_input/experimental_data/Gene_lists/AF/AF_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt
RES_FEATURES_AF=$DATASTORE/user_input/experimental_data/PDB_residue_features/AF/residueFeatures.csv

# ── Run ────────────────────────────────────────────────────────────────────
python scripts/run_workflow4_nativeNCLE_batch.py \
    --pdb_dir $STRUCTDIR/AF \
    --gene_list $GENELIST_AF \
    --outdir $OUTDIR/nativeNCLE_all \
    --organism Ecoli \
    --model AF \
    --ent_detection_method 3 \
    --g_threshold 0.6 \
    --density 1.0 \
    --residue_features_file $RES_FEATURES_AF \
    --reg_formula "cut_C_Rall ~ AA + region" \
    --design_matrix_file $OUTDIR/residue_dataframes_workflow4.csv \
    --nproc 16 \
    --logdir $OUTDIR/logs
```

Notes:

- This wrapper runs `scripts/run_nativeNCLE.py` in parallel.
- It accepts the current native-NCLE options used in Workflow 1, including `organism`, `model`, `ent_detection_method`, `g_threshold`, `density`, `cut_off`, `CG`, and `Calpha`. The per-structure `pdb_file` is set automatically for each matched PDB.
- Selection rule: the PDB stem (filename without `.pdb`) must be present in `--gene_list`.
- Per-protein NCLE outputs are created under `$OUTDIR/nativeNCLE_all/<ID>/`.
- By default, the wrapper passes each matched protein ID as the `gene` value to `run_nativeNCLE.py`, so generated NCLE feature filenames are keyed to the PDB stem unless you override `--gene`.
- The combined design matrix is generated automatically as a single file:
    - `$OUTDIR/residue_dataframes_workflow4.csv`

---

## Step 2. Run proteome-level logistic regression

Use `ProteomeLogisticRegression` to test enrichment while controlling for confounders.

### Python API example

```python
from NCLEdetector.statistics import ProteomeLogisticRegression

# ── Paths ──────────────────────────────────────────────────────────────────
DATASTORE    = "/scratch/ims86/NCLEdetector_Datastore"
OUTDIR       = f"{DATASTORE}/outputs/workflow4"
GENELIST_AF  = f"{DATASTORE}/user_input/experimental_data/Gene_lists/AF/AF_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt"

# ── Inputs ──────────────────────────────────────────────────────────────────
dataframe_file = f"{OUTDIR}/residue_dataframes_workflow4.csv"
reg_outdir     = f"{OUTDIR}/population_modeling"
ID             = "Ecoli_population"
reg_formula    = "cut_C_Rall ~ AA + region"

# ── Initialize and Run ───────────────────────────────────────────────────────────────
PReg = ProteomeLogisticRegression(
    dataframe_files=dataframe_file,
    outdir=reg_outdir,
    gene_list=GENELIST_AF,
    ID=ID,
    reg_formula=reg_formula,
)

PReg.load_data(
    sep='|',
    reg_var=['AA', 'region'],
    response_var='cut_C_Rall',
    var2binarize=['cut_C_Rall', 'region'],
    mask_column='mapped_resid',
)

reg_df = PReg.run()
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

## Minimal Workflow – 8: Proteome-wide regression for predicting experimental signals

You can provide parameters either directly as CLI flags, via a `--config` file, or both.
When both are provided, **CLI flags override config values** for the same parameter.

If you provide the optional batch-preprocessing arguments below, `scripts/run_population_modeling.py` will first call `scripts/run_workflow4_nativeNCLE_batch.py` to execute Workflow 4 Step 1, generate the per-protein native-NCLE outputs, and build the combined residue-level design matrix before running the proteome-level regression.

The config example below targets `scripts/run_population_modeling.py`, and the `batch_*` keys are the mechanism by which that wrapper forwards Step 1 parameters into `scripts/run_workflow4_nativeNCLE_batch.py`.

### Wrapper command

Use `scripts/run_population_modeling.py` to test whether the NCLE-derived `region` label is enriched for the experimental residue-level signal while controlling for amino-acid composition.

### Example with direct CLI flags:

```bash
# ── Activation ─────────────────────────────────────────────────────────────
source ~/.bashrc
conda activate ncledetector

# ── Paths ──────────────────────────────────────────────────────────────────
DATASTORE=/scratch/ims86/NCLEdetector_Datastore
STRUCTDIR=$DATASTORE/user_input/proteome_structures
OUTDIR=$DATASTORE/outputs/workflow4
GENELIST_AF=$DATASTORE/user_input/experimental_data/Gene_lists/AF/AF_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt
RES_FEATURES_AF=$DATASTORE/user_input/experimental_data/PDB_residue_features/AF/residueFeatures.csv

# ── Run ────────────────────────────────────────────────────────────────────
python scripts/run_population_modeling.py \
    --run_batch_native_ncle \
    --batch_pdb_dir $STRUCTDIR/AF \
    --batch_outdir $OUTDIR/nativeNCLE_all \
    --batch_organism Ecoli \
    --batch_model AF \
    --batch_ent_detection_method 3 \
    --batch_g_threshold 0.6 \
    --batch_density 1.0 \
    --batch_residue_features_file $RES_FEATURES_AF \
    --batch_design_matrix_file $OUTDIR/residue_dataframes_workflow4.csv \
    --batch_nproc 16 \
    --dataframe_files $OUTDIR/residue_dataframes_workflow4.csv \
    --outdir $OUTDIR/population_modeling \
    --gene_list $GENELIST_AF \
    --ID Ecoli_population \
    --reg_formula "cut_C_Rall ~ AA + region" \
    --logdir $OUTDIR/logs
```

### Example with config plus CLI override:

```bash
CONFIG=scripts/configs/workflow4_population_modeling_config.json

# Here, --batch_nproc 32 overrides batch_nproc from the config file.
python scripts/run_population_modeling.py \
    --config $CONFIG \
    --batch_nproc 32
```

Container equivalent (same config and optional CLI override):

```bash
CONFIG=scripts/configs/workflow4_population_modeling_config.json
DATASTORE=/scratch/ims86/NCLEdetector_Datastore

apptainer exec \
    --bind "$DATASTORE:$DATASTORE" \
    --bind "$PWD:$PWD" \
    --pwd "$PWD" \
    ncledetector-latest.sif \
    python scripts/run_population_modeling.py \
        --config "$CONFIG" \
        --batch_nproc 32
```

Config file example (matches `scripts/configs/workflow4_population_modeling_config.json`):

```json
{
    "run_batch_native_ncle": true,
    "batch_pdb_dir": "/scratch/ims86/NCLEdetector_Datastore/user_input/proteome_structures/AF",
    "batch_outdir": "/scratch/ims86/NCLEdetector_Datastore/outputs/workflow4/nativeNCLE_all",
    "batch_nproc": 16,
    "batch_organism": "Ecoli",
    "batch_model": "AF",
    "batch_ent_detection_method": 3,
    "batch_g_threshold": 0.6,
    "batch_density": 1.0,
    "batch_residue_features_file": "/scratch/ims86/NCLEdetector_Datastore/user_input/experimental_data/PDB_residue_features/AF/residueFeatures.csv",
    "batch_design_matrix_file": "/scratch/ims86/NCLEdetector_Datastore/outputs/workflow4/residue_dataframes_workflow4.csv",
    "dataframe_files": "/scratch/ims86/NCLEdetector_Datastore/outputs/workflow4/residue_dataframes_workflow4.csv",
    "outdir": "/scratch/ims86/NCLEdetector_Datastore/outputs/workflow4/population_modeling",
    "gene_list": "/scratch/ims86/NCLEdetector_Datastore/user_input/experimental_data/Gene_lists/AF/AF_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt",
    "ID": "Ecoli_population",
    "reg_formula": "cut_C_Rall ~ AA + region",
    "sep": "|",
    "reg_var": ["AA", "region"],
    "response_var": "cut_C_Rall",
    "var2binarize": ["cut_C_Rall", "region"],
    "mask_column": "mapped_resid",
    "logdir": "/scratch/ims86/NCLEdetector_Datastore/outputs/workflow4/logs",
    "log_level": "INFO"
}
```

Because the schema is unified, you can also run Workflow 4 Step 1 directly with the same config file:

```bash
python scripts/run_workflow4_nativeNCLE_batch.py \
    --config scripts/configs/workflow4_population_modeling_config.json \
```

The JSON/YAML config keys and their matching CLI flags are listed below. For the command-line wrapper, use the same key name with a `--` prefix; the only extra wrapper-only flag is `--config`.

| Term | Definition |
|------|------------|
| `run_batch_native_ncle` (`--run_batch_native_ncle`) | Optional switch telling `run_population_modeling.py` to launch Workflow 4 Step 1 via `run_workflow4_nativeNCLE_batch.py` before fitting the regression. |
| `batch_pdb_dir` (`--batch_pdb_dir`) | PDB directory forwarded to `run_workflow4_nativeNCLE_batch.py --pdb_dir` when Step 1 preprocessing is enabled. |
| `batch_outdir` (`--batch_outdir`) | Native-NCLE batch output directory forwarded to `run_workflow4_nativeNCLE_batch.py --outdir`. |
| `batch_nproc` (`--batch_nproc`) | Parallel native-NCLE job count forwarded to `run_workflow4_nativeNCLE_batch.py --nproc`. |
| `batch_allow_prefix_match` (`--batch_allow_prefix_match`) | Optional switch forwarded to `run_workflow4_nativeNCLE_batch.py --allow_prefix_match`. |
| `batch_dry_run` (`--batch_dry_run`) | Optional switch forwarded to `run_workflow4_nativeNCLE_batch.py --dry_run`. |
| `batch_chain` (`--batch_chain`) | Optional chain identifier forwarded to `run_workflow4_nativeNCLE_batch.py --chain`. |
| `batch_gene` (`--batch_gene`) | Optional gene/accession override forwarded to `run_workflow4_nativeNCLE_batch.py --gene`. |
| `batch_CG` (`--batch_CG`) | Optional switch forwarded to `run_workflow4_nativeNCLE_batch.py --CG`. |
| `batch_Calpha` (`--batch_Calpha`) | Optional switch forwarded to `run_workflow4_nativeNCLE_batch.py --Calpha`. |
| `batch_g_threshold` (`--batch_g_threshold`) | Gaussian entanglement score cutoff forwarded to the batch native-NCLE run. |
| `batch_density` (`--batch_density`) | Topoly triangulation density forwarded to the batch native-NCLE run. |
| `batch_cut_off` (`--batch_cut_off`) | Optional clustering cutoff forwarded to `run_workflow4_nativeNCLE_batch.py --cut_off`. |
| `batch_model` (`--batch_model`) | Structure source label forwarded to the batch native-NCLE run, typically `AF` or `EXP`. |
| `batch_ent_detection_method` (`--batch_ent_detection_method`) | Entanglement detection mode forwarded to the batch native-NCLE run. |
| `batch_residue_features_file` (`--batch_residue_features_file`) | Pipe-delimited residue feature table used by the batch wrapper to build the combined design matrix. |
| `batch_design_matrix_file` (`--batch_design_matrix_file`) | Output path for the combined residue-level design matrix produced by the batch wrapper. If `dataframe_files` is omitted, `run_population_modeling.py` uses this path as the regression input. |
| `batch_organism` (`--batch_organism`) | Organism setting forwarded to `run_workflow4_nativeNCLE_batch.py --organism`. |
| `batch_logdir` (`--batch_logdir`) | Optional log directory forwarded only to the batch wrapper. If omitted, the batch step reuses `logdir` when available. |
| `dataframe_files` (`--dataframe_files`) | Input residue-level design matrix path. This can be the single combined Workflow 4 matrix file or, in older workflows, a directory of per-protein residue tables. |
| `outdir` (`--outdir`) | Output directory where the regression summary CSV is written. |
| `gene_list` (`--gene_list`) | Cohort-definition file listing the protein identifiers to retain in the proteome-wide regression. The same file is also forwarded to the batch preprocessing step when `run_batch_native_ncle` is enabled. |
| `ID` (`--ID`) | Analysis identifier appended to output filenames and used as the regression run label. |
| `reg_formula` (`--reg_formula`) | Statsmodels-style logistic-regression formula used for the enrichment test. In the tutorial example, this models `cut_C_Rall` as a function of amino-acid identity and the NCLE-region label. The same formula is forwarded to the batch step to decide which columns must be retained in the design matrix. |
| `sep` (`--sep`) | Delimiter used by `run_population_modeling.py` when loading the design matrix file before regression preprocessing. |
| `reg_var` (`--reg_var`) | Predictor columns passed to `ProteomeLogisticRegression.load_data`, used when encoding regression covariates. |
| `response_var` (`--response_var`) | Response column passed to `ProteomeLogisticRegression.load_data` before model fitting. |
| `var2binarize` (`--var2binarize`) | Columns binarized during `load_data` preprocessing prior to fitting the regression model. |
| `mask_column` (`--mask_column`) | Residue index column used as the structural mask key during regression data loading. |
| `log_level` (`--log_level`) | Logging verbosity for the wrapper: `DEBUG`, `INFO`, `WARNING`, or `ERROR`. When batch preprocessing is enabled, the same level is forwarded to the batch wrapper unless `batch_logdir` is overridden separately. |
| `logdir` (`--logdir`) | Directory where the `run_population_modeling` log file is written. If omitted, logging defaults to `outdir`. When batch preprocessing is enabled, this directory is also reused for the batch wrapper unless `batch_logdir` is supplied. |

### Expected outputs

```bash
$OUTDIR/nativeNCLE_all/<ID>/
$OUTDIR/residue_dataframes_workflow4.csv
$OUTDIR/population_modeling/regression_results_Ecoli_population.csv
```

## I/O Reference for run_workflow4_nativeNCLE_batch.py

Each file below is listed once, followed by the column-level description table for tabular outputs.

### `$DATASTORE/user_input/proteome_structures/AF/*.pdb`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/user_input/proteome_structures/AF/*.pdb` | Proteome-scale AlphaFold structure library scanned by the batch runner. Each file provides the native coordinates for one protein processed by `run_nativeNCLE.py`. |

### `$DATASTORE/user_input/experimental_data/Gene_lists/AF/AF_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/user_input/experimental_data/Gene_lists/AF/AF_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt` | Cohort-definition file used to select which PDB stems are processed and retained in the final Workflow 4 design matrix. |

### `$DATASTORE/user_input/experimental_data/PDB_residue_features/AF/residueFeatures.csv`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/user_input/experimental_data/PDB_residue_features/AF/residueFeatures.csv` | Residue-level feature table providing canonical residue indices, amino-acid identities, and experimental-response columns that are merged with NCLE-derived region labels. |

### `$OUTDIR/nativeNCLE_all/<ID>/Native_clustered_HQ_GE_features/<GENE>_<ID>_<CHAIN>_uent_features.csv`

| I/O | File | File Description |
|---|---|---|
| Output | `$OUTDIR/nativeNCLE_all/<ID>/Native_clustered_HQ_GE_features/<GENE>_<ID>_<CHAIN>_uent_features.csv` | Per-protein representative-entanglement feature table produced by the Workflow 1 pipeline and used here to derive the binary `region` label. Here `<GENE>` is the value passed to `run_nativeNCLE.py --gene`; in this batch wrapper it defaults to the matched PDB stem unless you override it. This file uses the same schema documented in Workflow 1. |

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

### `$OUTDIR/residue_dataframes_workflow4.csv`

| I/O | File | File Description |
|---|---|---|
| Output | `$OUTDIR/residue_dataframes_workflow4.csv` | Combined residue-level design matrix written after all batch native-NCLE jobs complete successfully. This matrix is the direct input for the regression and Monte Carlo scripts below. |

| Column Name | Column Description |
|---|---|
| gene | Protein identifier matching the selected entries in `gene_list`. |
| mapped_resid | Canonical residue index used as the residue-level inclusion mask and lookup key for the `region` assignment. |
| uniprot_length | Full canonical protein length used downstream for population-size normalization. |
| AA | Single-letter amino-acid code for the residue. |
| region | Binary indicator equal to 1 when `mapped_resid` falls inside the NCLE-derived entangled region, otherwise 0. |
| cut_C_Rall | Binary experimental-response label indicating whether the residue carries a significant conformational-change signal. |

### Design matrix details

The batch runner executes the design-matrix build internally after all per-protein native-NCLE jobs finish successfully.

The matrix written to `$OUTDIR/residue_dataframes_workflow4.csv` is the direct input for the regression and Monte Carlo analyses below.

### Expected per-protein outputs

```
$OUTDIR/nativeNCLE_all/<ID>/
├── Native_GE/
├── Native_HQ_GE/
├── Native_clustered_HQ_GE/
└── Native_clustered_HQ_GE_features/
```

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

## I/O Reference for run_population_modeling.py

Each file below is listed once, followed by the column-level description table for tabular outputs.

### `$OUTDIR/residue_dataframes_workflow4.csv`

| I/O | File | File Description |
|---|---|---|
| Input | `$OUTDIR/residue_dataframes_workflow4.csv` | Residue-level design matrix loaded by `ProteomeLogisticRegression` and filtered to the genes in the requested cohort. |

### `$DATASTORE/user_input/experimental_data/Gene_lists/AF/AF_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/user_input/experimental_data/Gene_lists/AF/AF_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt` | Inclusion list constraining which proteins contribute residue rows to the proteome-wide regression fit. |

### `$OUTDIR/population_modeling/regression_results_<ID>.csv`

| I/O | File | File Description |
|---|---|---|
| Output | `$OUTDIR/population_modeling/regression_results_<ID>.csv` | Pipe-delimited regression summary written from the fitted statsmodels table after Workflow 4 proteome-level logistic regression completes. |

| Column Name | Column Description |
|---|---|
| var | Regression term name, including the intercept, amino-acid dummy variables, and the test variable `region`. |
| coef | Estimated regression coefficient for the term on the log-odds scale. |
| std err | Standard error of the coefficient estimate. |
| z | Wald z-statistic for the coefficient. |
| `P>\|z\|` | Two-sided p-value for the coefficient, recomputed at full precision from the z-statistic. |
| `[0.025` | Lower bound of the 95% confidence interval, stored on the odds-ratio scale in the saved table. |
| `0.975]` | Upper bound of the 95% confidence interval, stored on the odds-ratio scale in the saved table. |
| OR | Odds ratio for the fitted term, computed as `exp(coef)`. |
| ID | Analysis identifier copied from `--ID`. |
| n | Number of proteins retained in the fitted cohort after filtering. |

---

## Step 3. Run Monte Carlo subpopulation selection

Use `MonteCarlo` to partition proteins into groups and optimize an objective combining enrichment signal and group-size distribution constraints.

### Python API example

```python
from NCLEdetector.statistics import MonteCarlo

# ── Paths ──────────────────────────────────────────────────────────────────
DATASTORE    = "/scratch/ims86/NCLEdetector_Datastore"
OUTDIR       = f"{DATASTORE}/outputs/workflow4"
GENELIST_AF  = f"{DATASTORE}/user_input/experimental_data/Gene_lists/AF/AF_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt"

# ── Inputs ──────────────────────────────────────────────────────────────────
dataframe_file = f"{OUTDIR}/residue_dataframes_workflow4.csv"
mc_outdir      = f"{OUTDIR}/monte_carlo"
ID             = "Ecoli_population_mc"
reg_formula    = "cut_C_Rall ~ region + AA"

# ── Initialize and Run ───────────────────────────────────────────────────────────────
MC = MonteCarlo(
    dataframe_files=dataframe_file,
    outdir=mc_outdir,
    gene_list=GENELIST_AF,
    ID=ID,
    reg_formula=reg_formula,
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
    var2binarize=['cut_C_Rall', 'region'],
    mask_column='mapped_resid',
    ID_column='gene',
    Length_column='uniprot_length',
)

MC.run(encoded_df=MC.data, ID_column='gene')
```

---

## Minimal Workflow – 9: Candidate selection by Monte Carlo

You can provide parameters either directly as CLI flags, via a `--config` file, or both.
When both are provided, **CLI flags override config values** for the same parameter.

### Example with direct CLI flags:

```bash
# ── Activation ─────────────────────────────────────────────────────────────
source ~/.bashrc
conda activate ncledetector

# ── Paths ──────────────────────────────────────────────────────────────────
DATASTORE=/scratch/ims86/NCLEdetector_Datastore
OUTDIR=$DATASTORE/outputs/workflow4
GENELIST_AF=$DATASTORE/user_input/experimental_data/Gene_lists/AF/AF_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt

# ── Run ────────────────────────────────────────────────────────────────────
python scripts/run_montecarlo.py \
    --dataframe_files $OUTDIR/residue_dataframes_workflow4.csv \
    --outdir $OUTDIR/monte_carlo \
    --gene_list $GENELIST_AF \
    --ID Ecoli_population_mc \
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

### Example with config plus CLI override:

```bash
CONFIG=scripts/configs/workflow4_monte_carlo_config.json

# Here, --n_groups 6 overrides n_groups from the config file.
python scripts/run_montecarlo.py \
    --config $CONFIG \
    --n_groups 6
```

Container equivalent (same config and optional CLI override):

```bash
CONFIG=scripts/configs/workflow4_monte_carlo_config.json
DATASTORE=/scratch/ims86/NCLEdetector_Datastore

apptainer exec \
    --bind "$DATASTORE:$DATASTORE" \
    --bind "$PWD:$PWD" \
    --pwd "$PWD" \
    ncledetector-latest.sif \
    python scripts/run_montecarlo.py \
        --config "$CONFIG" \
        --n_groups 6
```

Config file example (matches `scripts/configs/workflow4_monte_carlo_config.json`):

```json
{
    "dataframe_files": "/scratch/ims86/NCLEdetector_Datastore/outputs/workflow4/residue_dataframes_workflow4.csv",
    "outdir": "/scratch/ims86/NCLEdetector_Datastore/outputs/workflow4/monte_carlo",
    "gene_list": "/scratch/ims86/NCLEdetector_Datastore/user_input/experimental_data/Gene_lists/AF/AF_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt",
    "ID": "Ecoli_population_mc",
    "reg_formula": "cut_C_Rall ~ region + AA",
    "response_var": "cut_C_Rall",
    "sep": "|",
    "reg_var": ["AA", "region"],
    "var2binarize": ["cut_C_Rall", "region"],
    "mask_column": "mapped_resid",
    "ID_column": "gene",
    "Length_column": "uniprot_length",
    "test_var": "region",
    "n_groups": 4,
    "steps": 100000,
    "C1": 1.0,
    "C2": 2.5,
    "beta": 0.05,
    "logdir": "/scratch/ims86/NCLEdetector_Datastore/outputs/workflow4/logs",
    "log_level": "INFO"
}
```

The JSON/YAML config keys and their matching CLI flags are listed below. For the command-line wrapper, use the same key name with a `--` prefix; the only extra wrapper-only flag is `--config`.

| Term | Definition |
|------|------------|
| `dataframe_files` (`--dataframe_files`) | Input residue-level design matrix path used to score candidate protein subpopulations. This can be the single combined Workflow 4 matrix file or a legacy directory of per-protein residue tables. |
| `outdir` (`--outdir`) | Output directory where the Monte Carlo state summaries, final gene lists, and optimization trajectories are written. |
| `gene_list` (`--gene_list`) | Universe of proteins eligible for assignment to Monte Carlo states or groups. |
| `ID` (`--ID`) | Analysis identifier appended to all Monte Carlo output filenames and log files and forwarded to `MonteCarlo(..., ID=...)`. |
| `reg_formula` (`--reg_formula`) | Regression formula used when scoring each candidate Monte Carlo state during optimization. |
| `sep` (`--sep`) | Delimiter used by `run_montecarlo.py` when loading the residue-level design matrix before preprocessing. |
| `reg_var` (`--reg_var`) | Predictor columns passed to `MonteCarlo.load_data`, used when encoding regression covariates for each candidate state. |
| `response_var` (`--response_var`) | Response variable forwarded to `MonteCarlo(..., response_var=...)` and then reused internally during `load_data` preprocessing and downstream Monte Carlo scoring. |
| `var2binarize` (`--var2binarize`) | Columns binarized during `load_data` preprocessing before state-specific model fitting. |
| `mask_column` (`--mask_column`) | Residue index column used as the structural mask key during Monte Carlo data loading. |
| `ID_column` (`--ID_column`) | Protein identifier column used during `load_data` preprocessing and when grouping rows during `MC.run(...)`. |
| `Length_column` (`--Length_column`) | Protein length column used during `load_data` preprocessing for the size-aware Monte Carlo workflow. |
| `test_var` (`--test_var`) | Primary predictor of interest whose enrichment effect is tracked during optimization. This value is used later by the Monte Carlo scoring methods, not just at initialization, and the class ensures the column is loaded even if it is omitted from `reg_var`. |
| `random` (`--random`) | Optional switch enabling random sampling mode instead of the structured Monte Carlo optimization path. |
| `n_groups` (`--n_groups`) | Number of protein groups or states to infer during the Monte Carlo partitioning. |
| `steps` (`--steps`) | Total number of Monte Carlo updates to perform. |
| `C1` (`--C1`) | Weight on the enrichment-driven component of the Monte Carlo objective function. |
| `C2` (`--C2`) | Weight on the population-size or distribution-penalty component of the Monte Carlo objective function. |
| `beta` (`--beta`) | Inverse-temperature or annealing parameter controlling acceptance behavior during optimization. |
| `linearT` (`--linearT`) | Optional switch using a linear temperature schedule during Monte Carlo optimization. |
| `log_level` (`--log_level`) | Logging verbosity for the wrapper: `DEBUG`, `INFO`, `WARNING`, or `ERROR`. |
| `logdir` (`--logdir`) | Directory where the `run_montecarlo` log file is written. If omitted, logging defaults to `outdir`. |

Typical Monte Carlo outputs include:

- `Final_step_reg_<ID>.csv`
- `State0_final_genelist_<ID>.txt` ... `StateN_final_genelist_<ID>.txt`
- `State0_final_traj_<ID>.csv` ... `StateN_final_traj_<ID>.csv`

## I/O Reference for run_montecarlo.py

Each file below is listed once, followed by the column-level description table for tabular outputs.

### `$OUTDIR/residue_dataframes_workflow4.csv`

| I/O | File | File Description |
|---|---|---|
| Input | `$OUTDIR/residue_dataframes_workflow4.csv` | Residue-level design matrix used to score candidate subpopulations during Monte Carlo optimization. |

### `$DATASTORE/user_input/experimental_data/Gene_lists/AF/AF_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt`

| I/O | File | File Description |
|---|---|---|
| Input | `$DATASTORE/user_input/experimental_data/Gene_lists/AF/AF_0.6g_C_Rall_spa50_LiPMScov50_ent_genes.txt` | Universe of proteins that can be assigned to Monte Carlo states or groups. |

### `$OUTDIR/monte_carlo/Final_step_reg_<ID>.csv`

| I/O | File | File Description |
|---|---|---|
| Output | `$OUTDIR/monte_carlo/Final_step_reg_<ID>.csv` | Final per-state regression summary at the last Monte Carlo step, combining the fitted enrichment term with optimization diagnostics and state-size statistics. |

| Column Name | Column Description |
|---|---|
| var | Name of the fitted regression term reported for the state summary, typically the primary test variable `region`. |
| coef | Final-step regression coefficient for the reported term. |
| std err | Standard error of the coefficient estimate. |
| z | Wald z-statistic for the final-step coefficient. |
| `P>|z|` | Two-sided p-value associated with the z-statistic. |
| `[0.025` | Lower bound of the 95% confidence interval on the log-odds scale. |
| `0.975]` | Upper bound of the 95% confidence interval on the log-odds scale. |
| state | Integer state or group identifier in the Monte Carlo partition. |
| beta | Final annealing or inverse-temperature value used when the state summary was written. |
| ks_stat_size | Kolmogorov-Smirnov statistic quantifying deviation of the state-size distribution from the target distribution. |
| E | Final Monte Carlo objective value for the state after applying enrichment and size-distribution penalties. |
| psize_mean | Mean sampled population size for the state across retained configurations. |
| psize_lb | Lower bound of the state-size interval tracked during optimization. |
| psize_ub | Upper bound of the state-size interval tracked during optimization. |
| step | Monte Carlo step index at which the final state summary was recorded. |
| OR | Odds ratio for the reported term, computed as `exp(coef)`. |
| OR_lb | Lower confidence bound for the odds ratio. |
| OR_ub | Upper confidence bound for the odds ratio. |
| ID | Analysis identifier copied from `--ID`. |
| n | Number of proteins assigned to the final state. |

### `$OUTDIR/monte_carlo/State*_final_genelist_<ID>.txt`

| I/O | File | File Description |
|---|---|---|
| Output | `$OUTDIR/monte_carlo/State*_final_genelist_<ID>.txt` | Plain-text file listing the final protein members of one inferred Monte Carlo state, one protein identifier per line. |

### `$OUTDIR/monte_carlo/State*_final_traj_<ID>.csv`

| I/O | File | File Description |
|---|---|---|
| Output | `$OUTDIR/monte_carlo/State*_final_traj_<ID>.csv` | Per-step trajectory of optimization statistics for one inferred Monte Carlo state, showing how enrichment and size penalties evolved over the run. |

| Column Name | Column Description |
|---|---|
| state | Integer state or group identifier tracked by this trajectory file. |
| step | Monte Carlo step index for the recorded summary row. |
| OR | Odds ratio of the primary test variable for that step. |
| pvalue | P-value associated with the per-step enrichment model. |
| psize_mean | Mean population size of the state at that step. |
| psize_lb | Lower bound of the state-size interval tracked at that step. |
| psize_ub | Upper bound of the state-size interval tracked at that step. |
| ks_stat_size | Kolmogorov-Smirnov statistic describing how well the state-size distribution matches the target at that step. |
| E | Monte Carlo objective value for the state at that step. |
| beta | Annealing or inverse-temperature value used at that step. |

---

## Step 4. Rank and refine candidate protein groups

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
| Monte Carlo outputs not found where expected | `--outdir` points to an unexpected location or mixes runs | Use a dedicated empty output directory and inspect generated filenames directly |

---

← [Workflow 3](workflow3_sim2exp.md) | [Back to Master Index](index.md) | Next → [Tutorial Examples](tutorial_examples.md)
