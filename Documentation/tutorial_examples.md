# NCLEdetector User Tutorial

A detailed, beginner-friendly guide to the three main analysis workflows in the **NCLEdetector** Python package.

---

## Table of Contents

1. [What NCLEdetector is for](#what-ncledetector-is-for)
2. [Who this tutorial is for](#who-this-tutorial-is-for)
3. [Before you begin](#before-you-begin)
4. [Recommended project layout](#recommended-project-layout)
5. [Workflow overview](#workflow-overview)
6. [Workflow 1: Identify native non-covalent lasso entanglements (NCLEs)](#workflow-1-identify-native-non-covalent-lasso-entanglements-ncles)
7. [Workflow 2: Detect changes in entanglement across simulation trajectories](#workflow-2-detect-changes-in-entanglement-across-simulation-trajectories)
8. [Workflow 3: Compare structural ensembles to high-throughput experimental data](#workflow-3-compare-structural-ensembles-to-high-throughput-experimental-data)
9. [Workflow 4: Population-level detection of misfolding involving native entanglements from high-throughput experiments](#workflow-4-population-level-detection-of-misfolding-involving-native-entanglements-from-high-throughput-experiments)
10. [Common pitfalls and practical advice](#common-pitfalls-and-practical-advice)
11. [What to verify against the package documentation](#what-to-verify-against-the-package-documentation)

---

## What NCLEdetector is for

**NCLEdetector** is a Python package for studying **non-covalent lasso entanglements (NCLEs)** in proteins and protein ensembles. It supports three major use cases:

1. **Identifying native entanglements** in a single protein structure.
2. **Tracking changes in entanglement** across simulation trajectories and organizing those changes into structural states.
3. **Connecting simulated structural ensembles to experimental data**, especially high-throughput conformational readouts such as **LiP-MS** and **XL-MS**.

A fourth related use case is also covered here:

4. **Testing whether experimentally observed conformational changes are enriched in proteins or regions containing native entanglements** at the proteome or population level.

This tutorial is written as a practical guide for a **new user** who is learning how to run these analyses from the command line and from Python.

---

## Who this tutorial is for

This tutorial assumes you:

- are working in a **Linux environment**;
- have already installed the software dependencies listed in the NCLEdetector materials/setup documentation;
- have created a **Miniconda or Conda environment** with NCLEdetector installed;
- are comfortable running terminal commands and basic Python scripts.

Even if you are new to NCLEdetector, the goal of this guide is to walk you through the logic of each workflow, the required inputs, the core API calls, and the expected outputs.

---

## Before you begin

### 1. Activate your environment

Open a terminal and move to your working directory:

```bash
cd /path/to/base/directory
conda activate NCLEdetector_env
```

This tutorial assumes all relative paths are defined from your current working directory.

### 2. Confirm your inputs are prepared correctly

Depending on the workflow, you may need some combination of the following:

- a cleaned **PDB** or AlphaFold structure;
- coarse-grained or all-atom trajectory files such as **DCD**;
- topology/reference files such as **PSF** and **COR**;
- secondary structure definitions from **STRIDE**;
- domain definitions from **CATH**;
- processed **LiP-MS** and/or **XL-MS** experimental files;
- mapping files connecting trajectory numbers to files;
- metadata tables for population-level modeling.

### 3. Clean your structural inputs first

For structure-based analyses, make sure your input structure has been preprocessed to:

- rebuild missing residues and atoms if needed;
- remove duplicate residues;
- ensure chain and residue numbering are sensible;
- match the coordinate representation expected by the analysis step.

The original protocol specifically recommends rebuilding missing atoms/residues using tools such as **Modeller** or **CHARMM** before running the native entanglement workflow.

---

## Recommended project layout

A consistent directory structure makes NCLEdetector much easier to use and debug. One reasonable layout is:

```text
project/
├── structures/
│   ├── 1zmr_model_clean.pdb
│   ├── 1zmr_model_clean_ca.psf
│   ├── 1zmr_model_clean_ca.cor
│   └── secondary_struc_defs.txt
├── trajectories/
│   ├── 1_prod.dcd
│   └── 1_prod_aa.dcd
├── metadata/
│   ├── domain_def.dat
│   ├── trajnum2file.txt
│   └── gene_list.txt
├── experimental/
│   ├── ecPGK_significant_LiPMS_peptide_R1_merged.xlsx
│   └── ecPGK_significant_XLMS_peptide_R1_merged.xlsx
├── nativeNCLE/
├── OP/
├── MSM/
├── BackMapping/
├── MassSpec_ConsistencyTest/
└── population_modeling/
```

You do not have to use this exact structure, but keeping structures, trajectories, metadata, and results separated will save a lot of confusion later.

---

## Workflow overview

The full NCLEdetector analysis logic can be thought of as four connected layers:

### A. Native structure layer
Start with one experimental or predicted structure and identify its **native NCLEs**.

### B. Trajectory layer
For simulation data, compute order parameters such as **Q**, **G**, and **K**, then cluster entanglement changes and identify **metastable states**.

### C. Simulation-to-experiment layer
If you have LiP-MS or XL-MS data, test whether specific metastable states are structurally consistent with the experimental signals.

### D. Population/statistical layer
Across many proteins, ask whether residues or proteins associated with native entanglements are more likely to show experimentally observed conformational changes.

The sections below walk through each of these in detail.

---

# Workflow 1: Identify native non-covalent lasso entanglements (NCLEs)

## Goal

This workflow identifies native NCLEs in a single reference structure, optionally filters them, clusters redundant entanglements, and computes structural features for the representative entanglements.

## Typical runtime

Usually **less than 1 minute to 10 minutes**, depending on protein size.

## Step 1. Prepare and activate your environment

```bash
cd /path/to/base/directory
conda activate NCLEdetector_env
```

## Step 2. Prepare a cleaned structure

Before running NCLEdetector, prepare a cleaned structure file:

- input can be a **PDB** or AlphaFold-derived structure;
- rebuild missing residues and atoms if needed;
- remove duplicate residues;
- verify chain labels and residue numbering;
- save the cleaned structure to a known path.

Example:

```text
./1zmr_model_clean.pdb
```

## Step 3. Start a Python session or create a script

You can run the native workflow interactively:

```bash
python
```

or place the commands below into a Python script such as:

```text
run_native_ncle.py
```

Using a script is usually better for reproducibility.

## Step 4. Calculate native entanglements

Import the required class and define your inputs:

```python
from NCLEdetector.gaussian_entanglement import GaussianEntanglement
from NCLEdetector.clustering import ClusterNativeEntanglements

# Input paths and identifiers
pdb = "./1zmr_model_clean.pdb"
native_outdir = "./nativeNCLE/Native_GE"
ID = "1zmr"

# Initialize objects
ge = GaussianEntanglement(g_threshold=0.6, density=0.0, Calpha=False, CG=False)
clustering = ClusterNativeEntanglements(organism="Ecoli")

# Run the native entanglement calculation
ge.calculate_native_entanglements(clean_pdb=pdb, outdir=native_outdir, ID=ID)
```

### What this step does

This step scans the input structure and identifies the set of **native non-covalent lasso entanglements** present in the reference conformation.

### Main output

At minimum, this step produces a file describing the detected NCLEs. In the original protocol this is described as a CSV-like output; filenames may vary depending on your package version.

### Practical note

If your structure is all-atom, keep `Calpha=False` and `CG=False` unless your specific use case requires otherwise.

## Step 5. Optionally filter for high-quality entanglements

This step is particularly important when using **predicted structures** such as AlphaFold models, or when you want to remove slipknots or low-confidence entanglements.

```python
# Define additional inputs
native_HQ_outdir = "./nativeNCLE/Native_HQ_GE"
NCLE_file = "./nativeNCLE/Native_GE/1zmr_model_clean_ca_GE.txt"

# Filter the NCLEs
ge.select_high_quality_entanglements(
    NCLE_file,
    pdb,
    outdir=native_HQ_outdir,
    ID=ID,
    model="EXP"
)
```

### What this step does

This filtering stage removes lower-quality or potentially artifactual entanglements.

### Main output

A filtered entanglement file containing the higher-confidence NCLEs.

### Important note about non-canonical regions

If the structure includes residues that are **not part of the canonical protein sequence**—for example:

- chimeric fusion segments,
- inserted tags,
- engineered regions,
- unresolved sequence-to-structure mismatches,

then the protocol notes that a **mapping file** can be supplied so that only entanglements involving the canonical region are retained.

## Step 6. Cluster NCLEs into non-redundant representative entanglements

Many entanglements can be structurally redundant or degenerate variants of the same topological motif. This step reduces the set to a more interpretable non-redundant representation.

```python
from NCLEdetector.clustering import ClusterNativeEntanglements

native_clustered_HQ_outdir = "./nativeNCLE/Native_clustered_HQ_GE"
NCLE_file = "./nativeNCLE/Native_HQ_GE/1zmr.csv"
outfile = "1zmr.csv"

clustering = ClusterNativeEntanglements(organism="Ecoli")
clustering.Cluster_NativeEntanglements(
    NCLE_file,
    outdir=native_clustered_HQ_outdir,
    outfile=outfile
)
```

### What this step does

This clustering step groups redundant native entanglements and selects **representative NCLEs**.

### Main output

A file containing:

- representative entanglements;
- associated degenerate loops;
- a reduced set suitable for downstream feature generation and interpretation.

## Step 7. Compute structural features for representative NCLEs

Once you have representative entanglements, you can calculate features that describe their structural complexity.

```python
from NCLEdetector.entanglement_features import FeatureGen

pdb = "./1zmr_model_clean.pdb"
native_GQ_feature_outdir = "./nativeNCLE/Native_clustered_HQ_GE_features"
cluster_file = "./nativeNCLE/Native_clustered_HQ_GE/1zmr.csv"

FGen = FeatureGen(pdb, outdir=native_GQ_feature_outdir, cluster_file=cluster_file)
EntFeatures = FGen.get_uent_features(gene='P00558', chain='A', pdbid='1ZMR')
```

### What this step does

This computes descriptive features for each representative NCLE, such as:

- supercoiling characteristics of the entangled terminus;
- contact patterns between the loop and the rest of the structure;
- other structural descriptors capturing entanglement complexity.

### Main output

A feature table for the representative entanglements.

## Step 8. Visualize the entanglements

Use molecular graphics software such as:

- **VMD**
- **PyMOL**

to inspect the representative NCLEs visually.

### Why this matters

Visual inspection is important for:

- validating that the identified entanglements make structural sense;
- understanding loop/crossing geometry;
- preparing figures for presentations or publications.

---

# Workflow 2: Detect changes in entanglement across simulation trajectories

## Goal

This workflow analyzes simulation trajectories to detect changes in entanglements over time, remove redundancy, and build a coarse structural-state model using order parameters and Markov state modeling.

## Typical runtime

This can take **hours to several days**, depending on:

- number of trajectories,
- trajectory length,
- protein size,
- number of frames analyzed,
- whether coarse-grained or all-atom data are used.

## Step 9. Compute order parameters Q, G, and K

This is the starting point for trajectory-based analysis.

### Required inputs

You will generally need:

- a **PSF** file;
- a **DCD** trajectory;
- a **COR** reference structure;
- a file defining **secondary structure elements**;
- a file defining **domain boundaries**.

### Example setup

```python
from NCLEdetector.order_params import CalculateOP

Traj = 1
PSF = "./1zmr_model_clean_ca.psf"
DCD = "./1_prod.dcd"
ID = "1ZMR"
COR = "./1zmr_model_clean_ca.cor"
sec_elements = "./secondary_struc_defs.txt"
domain = "./domain_def.dat"
outdir = "./run_OP_on_simulation_traj_last335frames/"
start = 6332

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

### Compute Q

```python
Qdata_dict = CalcOP.Q()
```

`Q` is the **fraction of native contacts**.

### Compute G

```python
Gdata_dict = CalcOP.G(topoly=True, Calpha=True, CG=True, nproc=10)
```

`G` captures the **fraction of native contacts with a change in entanglement**.

In the original protocol, this step uses:

- **Topoly**,
- Cα-defined native contacts,
- a coarse-grained trajectory.

### Compute K

```python
Kdata_dict = CalcOP.K()
```

`K` is the **mirror-symmetry order parameter**.

### Main outputs

This stage creates per-order-parameter output directories containing at minimum:

- a CSV time series for each metric;
- additional metadata, especially for `G`, often stored in binary `pkl` files.

### Pause point

For large trajectories, this step can take a long time. If needed, **downsample the trajectory** to reduce runtime.

## Step 10. Identify and remove artificial mirror conformations

Before clustering entanglement changes, remove frames or trajectories that are mirror artifacts rather than physically relevant structural states.

### Recommended approach

1. Examine the steady-state region of your trajectories.
2. Apply cutoffs on **Q** and **K** to flag candidate mirror states.
3. Visually inspect flagged structures.
4. Remove trajectories or frames that are clearly mirror images.

### Important note

The protocol states that the cutoff values used in the example were tuned for **ecPGK**. Do not assume they transfer perfectly to a different protein system.

## Step 11. Cluster non-native entanglement changes

Once you have computed `G` and filtered trajectory artifacts, cluster the changes in entanglement across all trajectories.

### Example setup

```python
from NCLEdetector.clustering import ClusterNonNativeEntanglements

pkl_file_path = "./OP/G/Combined_GE/"
trajnum2pklfile_path = "./trajnum2file.txt"
traj_dir_prefix = "/path/to/dir/containing/dcds/"
outdir = "./nonnative_entanglement_clustering"

clustering_NNents = ClusterNonNativeEntanglements(
    pkl_file_path=pkl_file_path,
    trajnum2pklfile_path=trajnum2pklfile_path,
    traj_dir_prefix=traj_dir_prefix,
    outdir=outdir
)

clustering_NNents.cluster(start_frame=6332)
```

### What this step does

This runs a multistep clustering procedure across all selected trajectories to identify **non-redundant changes in entanglement**.

### Main outputs

The protocol indicates this stage produces:

- summary tables of representative entanglement changes;
- structural fingerprints;
- cluster memberships;
- cluster probabilities;
- a compressed archive of clustering inputs and mappings;
- a distribution plot of loop/crossing residues;
- a text representation of the clustering tree.

### Visualization

After clustering, inspect non-redundant entanglement changes in VMD or another structure viewer.

### Pause point

This step can be both **time-intensive** and **memory-intensive**. Make sure your workstation has enough RAM if the number of raw entanglement-change events is large.

## Step 12. Build a Markov state model (MSM)

This step uses the order-parameter outputs to organize structures into microstates and metastable states.

### Example setup

```python
from NCLEdetector.clustering import MSMNonNativeEntanglementClustering

outdir = "./run_MSM"
ID = "1ZMR"
OPpath = "./run_OP_on_simulation_traj_Allframes/"
n_large_states = 10

MSM = MSMNonNativeEntanglementClustering(
    outdir=outdir,
    ID=ID,
    OPpath=OPpath,
    n_large_states=n_large_states
)
MSM.run()
```

### What this step does

This builds an MSM-based decomposition of the trajectory data into:

- microstates;
- metastable states;
- a coarse map of the accessible conformational landscape.

### Main outputs

The protocol describes outputs including:

- a CSV mapping each structure to a microstate and metastable state;
- plots of the order-parameter landscape and metastable-state assignments.

### Critical step

Try **multiple values** of `n_large_states`.

Too few states may oversimplify the landscape. Too many may reduce interpretability. The protocol notes that **15 or fewer** often works well in practice.

Also note:

- the final number of states can be lower than the requested value if empty states are discarded;
- the default lag time is often set to **1 frame** for coarse visualization or exploratory grouping;
- if you want to interpret actual kinetics, you should explicitly test for **Markovian behavior** and choose lag time carefully.

## Step 13. Analyze metastable-state behavior

Once the MSM has been built, you can summarize how states behave over time and compare pathways between simulation conditions.

### 13a. State distribution visualization

The state-distribution surface is automatically generated by the MSM step.

### 13b. Plot state probability evolution

```python
from NCLEdetector.statistics import MSMStats

outdir = "./MSM_StateProbabilityStats"
msm_meta_file = "MSM/1ZMR_prod_meta_set_A80%Native.csv"
meta_set_file = "MSM/1ZMR_prod_meta_set.csv"
traj_type_col = "traj_type_A80%Native"
traj_type_list = ['A', 'B']
rm_traj_list = []

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

### Main output

A time series describing the probabilities of the metastable states over simulation time.

### 13c. Visualize representative metastable structures

Inspect representative structures for each metastable state in VMD.

### 13d. Compare folding pathways and Jensen-Shannon divergence

You can compare structural-ensemble evolution between two conditions by analyzing folding pathways and their divergence.

```python
import pandas as pd
from NCLEdetector.statistics import FoldingPathwayStats

outdir = "./Foldingpathway_A80%Native"
msm_meta_file = "MSM/1ZMR_prod_meta_set_A80%Native.csv"
meta_set_file = "MSM/1ZMR_prod_meta_set.csv"
traj_type_col = "traj_type_A80%Native"
rm_traj_list = []

msm_data = pd.read_csv(msm_meta_file)
FP = FoldingPathwayStats(
    msm_data=msm_data,
    meta_set_file=meta_set_file,
    traj_type_col=traj_type_col,
    outdir=outdir,
    traj_list=rm_traj_list
)

folding_pathways = FP.post_trans()
JS_divergence = FP.JS_divergence()
```

### Main outputs

This step produces:

- folding pathways and their probabilities;
- a time series of Jensen-Shannon divergence (JSD).

### Interpretation

Low JSD suggests the two conditions explore similar state distributions or pathway behavior, while higher JSD suggests divergence between the ensembles.

---

# Workflow 3: Compare structural ensembles to high-throughput experimental data

## Goal

This workflow asks whether structural ensembles derived from simulation are **consistent with experimental conformational signals**, especially from **LiP-MS** and **XL-MS**.

## Typical runtime

Roughly **5 to 24 hours**, depending on:

- ensemble size,
- number of states,
- number of frames,
- number of experimental signals,
- whether back-mapping is required.

## Step 14. Prepare experimental inputs

Before running the simulation-to-experiment consistency test, prepare processed experimental files.

Examples include:

- LiP-MS peptide or residue-level files;
- XL-MS peptide or cross-link files.

Make sure these are already cleaned, statistically filtered, and formatted the way the package expects.

## Step 15. Back-map coarse-grained trajectories to all-atom structures (if needed)

If your simulations are coarse-grained Cα trajectories, you must reconstruct all-atom models before computing some structure-to-experiment metrics.

### 15a. Save framewise PDBs

For each frame of interest, save a structure file.

### 15b. Back-map the structures

```python
from NCLEdetector.change_resolution import BackMapping

Outdir = "./BackMapping/"
cg_pdb = "./1zmr_model_clean_ca.cor"
aa_pdb = "./1zmr_model_clean.pdb"
ID = "1ZMR"

backMapper = BackMapping(outdir=Outdir)
backMapper.backmap(cg_pdb=cg_pdb, aa_pdb=aa_pdb, ID=ID)
```

### What this step does

This reconstructs full all-atom structures consistent with:

- the coarse-grained Cα input;
- the native reference structure.

### Main outputs

Depending on configuration, this can produce:

- reconstructed all-atom structures;
- intermediate PD2/Pulchra reconstruction outputs;
- OpenMM minimization logs;
- energy-minimized final structures.

### Pause point

Inspect reconstructed structures before using them downstream. After validation, collate them into an all-atom DCD if needed.

## Step 16. Calculate SASA, Jwalk distances, and cross-link probability

These metrics are needed for comparing structural ensembles to LiP-MS and XL-MS data.

### Example setup

```python
from NCLEdetector.order_params import CalculateOP

Traj = 1
PSF = "./1zmr_model_clean.pdb"
DCD = "./1_prod_aa.dcd"
ID = "1ZMR"
outdir = "./run_OP_on_simulation_traj_last335frames/"
start = 6332

CalcOP = CalculateOP(
    outdir=outdir,
    Traj=Traj,
    ID=ID,
    psf=PSF,
    dcd=DCD,
    start=start
)

CalcOP.SASA()
CalcOP.runJwalk('/path/to/backmapped/pdb')
XPdata_dict = CalcOP.XP(pdb_file='/path/to/AA/ref/PDBfile')
```

### What this step does

This computes:

- **SASA** for solvent exposure;
- **Jwalk distances** relevant to cross-link accessibility;
- **XP**, a cross-linking probability metric.

### Main output

The protocol describes these outputs as compressed binary files or other structured metric outputs required for the consistency test.

### Pause point

These analyses can be expensive. The original protocol recommends, for the example system:

- analyzing only the **last 50 ns**;
- downsampling every **20 frames**.

That exact rate may not be optimal for your system. Choose a downsampling interval that reasonably reduces autocorrelation while retaining the long-lived structural states you care about.

## Step 17. Run the consistency test between metastable states and experiments

This is the key step that integrates simulation and experiment.

### Example setup

```python
from NCLEdetector.compare_sim2exp import MassSpec

outdir = "./MassSpec_ConsistencyTest/"
msm_data_file = "./MSM/1ZMR_prod_MSMmapping.csv"
meta_dist_file = "./MSM/1ZMR_prod_meta_dist.npy"
LiPMS_exp_file = "./ecPGK_significant_LiPMS_peptide_R1_merged.xlsx"
sasa_data_file = "./run_OP_on_simulation_traj_last335frames/SASA/SASA.npy"
XLMS_exp_file = "./ecPGK_significant_XLMS_peptide_R1_merged.xlsx"
dist_data_file = "./run_OP_on_simulation_traj_last335frames/Jwalk/Jwalk.npy"
cluster_data_file = "./nonnative_entanglement_clustering/cluster_data_topoly_linking_number.npz"
OPpath = "./run_OP_on_simulation_traj_last335frames/"
AAdcd_dir = "/path/to/backmapped/dcds/"
native_AA_pdb = "./1zmr_model_clean.pdb"
state_idx_list = [4, 6, 8]
prot_len = 387
n_analysis_frames = 335
rm_traj_list = []
native_state_idx = 9
start = 6332
ID = "1ZMR"

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
    n_analysis_frames=n_analysis_frames,
    rm_traj_list=rm_traj_list,
    native_state_idx=native_state_idx,
    outdir=outdir,
    ID=ID,
    start=start
)

consist_data_file, consist_result_file = MS.LiP_XL_MS_ConsistencyTest()
MS.select_rep_structs(consist_data_file, consist_result_file, total_traj_num_frames=335, n_analysis_frames=n_analysis_frames)
```

### What this step does

This tests whether specific metastable states are statistically consistent with the supplied experimental LiP-MS and XL-MS signals.

### Main outputs

The protocol indicates this stage produces:

- Excel workbooks summarizing the statistical tests;
- consistency-test outputs for each metastable state;
- selected representative structures from the best-supported ensembles.

### Final step

Visualize the representative structures chosen from the experimentally consistent ensembles.

---

# Workflow 4: Population-level detection of misfolding involving native entanglements from high-throughput experiments

## Goal

This workflow is designed for **large-scale experimental datasets**, such as proteome-wide differential LiP-MS studies, where you want to test whether proteins or regions involving native entanglements are statistically associated with conformational change.

## Typical runtime

Anywhere from **~5 hours to multiple days**, depending on dataset size and the number of proteins analyzed.

## Conceptual overview

This section is independent of the trajectory-based simulation-to-experiment comparison above.

The general logic is:

1. identify native entanglements for all proteins of interest;
2. process experimental data to identify significant conformational signals;
3. define confounders;
4. build residue-level data tables;
5. run regression models;
6. use Monte Carlo methods to identify subpopulations with especially strong enrichment.

## Step 1. Compute representative NCLEs and features for all proteins

For every protein observed in the high-throughput experiment, run the native-structure workflow from **Workflow 1**.

In practice, this means repeating the following for each protein:

- structure preparation;
- native NCLE detection;
- optional filtering;
- native NCLE clustering;
- feature generation.

## Step 2. Process the LiP-MS data

Convert raw or processed LiP-MS results into a set of residues or sites showing **statistically significant conformational change**.

### Critical point

The protocol emphasizes that there are **many valid ways** to process LiP-MS data, and that this choice can strongly affect downstream statistical conclusions.

Be consistent in how you define significance, mapping, and residue-level assignments.

## Step 3. Define confounding variables

Decide which covariates you need to control for.

Examples given in the protocol include:

- amino acid identity or composition;
- solvent accessibility.

Other confounders may be relevant depending on the experiment.

## Step 4. Build residue-level data frames

For each protein, construct a table where each row corresponds to a residue and the columns include:

- a **response variable** indicating whether the residue shows a significant conformational change;
- a binary variable indicating whether the residue is in an **entangled region**;
- one or more confounding variables.

### Example conceptual columns

```text
residue_index | cut_C_Rall | region | AA | SASA | mapped_resid | gene | uniprot_length
```

where, for example:

- `cut_C_Rall = 1` means a significant conformational signal was observed;
- `region = 1` means the residue belongs to the entangled region.

## Step 5. Run proteome-level logistic regression

Use the `ProteomeLogisticRegression` class to model whether the odds of conformational change depend on entanglement-associated regions while adjusting for covariates.

### Example setup

```python
from NCLEdetector.statistics import ProteomeLogisticRegression

dataframe_files = "/path/to/dataframe/files"
outdir = "./population_modeling/"
gene_list = "/path/to/gene/list.txt"
ID = "Ecoli_noChaperones"
reg_formula = "cut_C_Rall ~ AA + region"

ProtRegession = ProteomeLogisticRegression(
    dataframe_files=dataframe_files,
    outdir=outdir,
    gene_list=gene_list,
    ID=ID,
    reg_formula=reg_formula
)

ProtRegession.load_data(
    sep='|',
    reg_var=['AA', 'region'],
    response_var='cut_C_Rall',
    var2binarize=['cut_C_Rall', 'region'],
    mask_column='mapped_resid'
)

reg_df = ProtRegession.run()
```

### What this step does

This fits a regression model using a formula of the form:

```text
y ~ x1 + x2 + ... + xn
```

to estimate the association between entangled regions and experimentally observed conformational change.

### How to interpret the results

For a coefficient in the fitted model:

- `exp(coefficient)` gives an **odds ratio**;
- the coefficient p-value reflects the significance of that association.

If the coefficient for `region` is positive and significant, that suggests residues in entangled regions are more likely to show conformational-change signals, after adjusting for the included covariates.

### Critical step

If you perform these tests across multiple conditions, subsets, or experiments, apply an **FDR correction** to control false positives.

## Step 6. Use Monte Carlo selection to identify high-risk protein subpopulations

The Monte Carlo workflow is meant to identify subsets of proteins whose entanglements are especially associated with misfolding-like experimental signals.

### Example setup

```python
from NCLEdetector.statistics import MonteCarlo

dataframe_files = "/path/to/dataframe/files"
outdir = "./monte_carlo/"
gene_list = "/path/to/gene/list.txt"
ID = "Ecoli_noChaperones"
reg_formula = "cut_C_Rall ~ AA + region"

MC = MonteCarlo(
    dataframe_files=dataframe_files,
    outdir=outdir,
    gene_list=gene_list,
    ID=ID,
    reg_formula=reg_formula
)

MC.load_data(
    sep='|',
    reg_var=['AA', 'region'],
    response_var='cut_C_Rall',
    var2binarize=['cut_C_Rall', 'region'],
    mask_column='mapped_resid',
    ID_column='gene',
    Length_column='uniprot_length'
)

MC.run(encoded_df=MC.data, ID_column='gene')
```

### Main outputs

The protocol states this step outputs:

- objective-function statistics across Monte Carlo steps;
- odds-ratio statistics for each bin or group over the course of the simulation.

## Step 7. Rank the groups after convergence

After the Monte Carlo simulation reaches a steady state:

1. take the **last 100 Monte Carlo steps**;
2. compute the average odds ratio for each group;
3. rank the groups by this average.

The top-ranking group should contain proteins most strongly associated with experimentally observed conformational changes in regions involving native entanglements.

## Step 8. Refine candidate selection across repeated simulations

To make the selection more robust:

- run **multiple independent Monte Carlo simulations**;
- keep only proteins that appear in the top odds-ratio group in **more than 70%** of runs.

This helps reduce sensitivity to stochastic effects in the Monte Carlo procedure.

---

# Common pitfalls and practical advice

## 1. Path mismatches are the most common source of errors

Double-check every input path before running a long analysis.

Especially verify:

- structure files;
- trajectory files;
- output directories;
- mapping files;
- experimental spreadsheet paths.

## 2. Keep naming consistent

Use the same identifier consistently across files when possible:

- `ID = "1ZMR"`
- `Traj = 1`
- filenames containing the same protein or trajectory number.

This makes downstream bookkeeping much easier.

## 3. Save scripts instead of relying only on interactive sessions

Interactive use is helpful for testing, but serious analyses should be placed in version-controlled Python scripts.

## 4. Downsample when appropriate

Some steps, especially order-parameter calculation and structural ensemble analysis, can become very slow on large trajectories.

If the science allows it, downsample to reduce:

- runtime,
- disk usage,
- memory usage.

## 5. Always visually inspect key structural results

Do not rely only on tables and scores. At several stages, visualization is essential:

- native NCLE detection;
- mirror-artifact removal;
- non-native entanglement clustering;
- metastable-state interpretation;
- representative-structure validation.

## 6. Choose model parameters deliberately

Important tunable parameters include:

- `g_threshold` for native entanglement detection;
- `start` frame for trajectory analysis;
- `n_large_states` for MSM construction;
- downsampling intervals;
- significance thresholds for LiP-MS/XL-MS preprocessing.

These are analysis choices, not just technical settings.

## 7. Monitor memory usage for clustering steps

The clustering of entanglement changes can require substantial memory, especially when the number of candidate changes is large.

## 8. Keep intermediate outputs

Do not delete intermediate files too early. Many downstream steps depend on:

- combined `G` metadata;
- clustering archives;
- MSM mapping files;
- SASA and Jwalk outputs;
- processed experimental tables.

---

# What to verify against the package documentation

This tutorial is based on the written protocol and preserves its intended workflow, but you should still verify several details against the current NCLEdetector codebase and GitHub documentation before running a production analysis.

In particular, confirm:

1. **Exact class and method names** in your installed package version.
2. **Required argument names and capitalization**.
3. **Expected file extensions and output filenames**.
4. **Whether certain steps require additional imports**, such as `pandas`.
5. **Any package-version-specific changes** to clustering, MSM, or mass-spec comparison APIs.

This is especially important because protocols often describe a stable conceptual workflow while implementation details may change between package versions.

---

# Final checklist for a new user

Before starting a full run, make sure you can answer **yes** to all of the following:

- Is my Conda environment active?
- Are my structure files cleaned and consistent?
- Do I know whether I am working with all-atom or coarse-grained data?
- Are my trajectory and metadata files in the expected locations?
- Have I tested the workflow on a small example first?
- Do I know which outputs from one step are required by the next?
- Have I planned enough storage, memory, and runtime for the larger steps?

If yes, you are ready to begin using NCLEdetector productively.

---

## Suggested first run for beginners

If you are completely new to NCLEdetector, the best starting point is:

1. run **Workflow 1** on a single cleaned structure;
2. inspect the representative NCLEs visually;
3. compute features;
4. only then move on to the trajectory and experimental workflows.

That approach will help you learn the package in manageable stages.

---

## Source note

This tutorial was derived from the uploaded NCLEdetector execution protocol and reorganized into a beginner-facing markdown guide.
