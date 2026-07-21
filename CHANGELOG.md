# Changelog

All notable changes to EntDetect will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- **BREAKING** `EntDetect.compare_sim2exp.MassSpec` — renamed `last_num_frames` → `n_analysis_frames` (trailing-frame analysis window shared by the MSM/G/Q/SASA/XP slicing) and `EntDetect.order_params.CollectOP` / `MassSpec`'s `n_frames` → `sasa_xp_frames_per_traj` (frames per trajectory stored in each SASA/XP file), to disambiguate the two previously same-named-sounding but distinct concepts. Updated `run_compare_sim2exp.py` CLI flags (`--last_num_frames`→`--n_analysis_frames`, `--n_frames`→`--sasa_xp_frames_per_traj`), config JSON keys, the workflow3 SLURM script, and documentation accordingly.

### Fixed
- `scripts/run_compare_sim2exp.py` — `select_rep_structs(...)` was called with hardcoded `last_num_frames=67` instead of the config-derived value, causing representative-structure frame indices (and `info.txt` frame numbers) to be offset from the correct trailing-frame window whenever the configured analysis window differed from 67 frames.

## [2.0.0] - 2026-07-08

### Added
- **Containerization** (`container/`): reproducible Docker (`Dockerfile`) and Singularity/Apptainer (`apptainer.def`) images that expose all minimal-workflow console-scripts directly, an exact `environment.lock.yml`, a `smoke_test.sh`, and a `README.md`. GitHub Actions workflow (`.github/workflows/container.yml`) builds and pushes both to GHCR on version tags.
- `scripts/run_OP_on_simulation_traj.py` — new `--end` flag (backward-compatible; default = end of trajectory) so a frame window can be selected from the config together with `--start`.
- **`--config` (JSON/YAML) support** with CLI-over-config merge for `run_compare_sim2exp.py`, `run_population_modeling.py`, and `run_montecarlo.py`, via the new shared helper `scripts/_cli_config.py`.
- `run_population_modeling.py` — optional integrated batch preprocessing (`--run_batch_native_ncle` and `batch_*` keys) that runs Workflow 4 Step 1 before the regression, plus configurable `load_data(...)` parameters (`sep`, `reg_var`, `response_var`, `var2binarize`, `mask_column`).
- `run_montecarlo.py` — configurable `load_data(...)` parameters (`sep`, `reg_var`, `response_var`, `var2binarize`, `mask_column`, `ID_column`, `Length_column`).
- `MassSpec` — automatic collection of per-trajectory SASA/XP arrays via `CollectOP` when a `sasa_dir`/`xp_dir` is provided (SASA always; Jwalk only with `--collect_jwalk_npy`, otherwise XP streaming).
- **CG benchmark harness**: `scripts/benchmarks/setup_cg_benchmarks.py`, per-protein configs under `scripts/configs/benchmarks/`, config-driven SLURM arrays under `assets/slurm/scripts/benchmarks/` (MW1–4), and `Documentation/Benchmarks.md`.

### Changed
- **BREAKING** `run_montecarlo.py` — renamed `--outpath`→`--outdir` and `--tag`→`--ID` to match the underlying `MonteCarlo` API and the rest of the toolkit.
- **BREAKING** `run_population_modeling.py` — renamed `--tag`→`--ID`.
- **BREAKING** `run_compare_sim2exp.py` — removed the dead `--start`/`--end`/`--stride` flags (they were never used in the consistency test).
- **BREAKING** `EntDetect.statistics.MonteCarlo.load_data(...)` — no longer takes `response_var` (it is reused from the constructor); `test_var` is auto-appended to the loaded columns when omitted from `reg_var`.
- **BREAKING** `EntDetect.statistics` — `MSMStats` and `FoldingPathwayStats` moved to a one-init / path-based method API.
- `EntDetect.compare_sim2exp.MassSpec.__init__` — added `sasa_dir`, `n_traj`, `n_frames`, `collect_jwalk_npy` parameters for in-constructor collection.
- Documentation: substantial Workflow 2/3/4 tutorial and I/O-reference updates; fixed all broken table-of-contents anchors across the workflow tutorials.

### Fixed
- `EntDetect/order_params.py` — G computation crashed with `Cluster_NativeEntanglements() got an unexpected keyword argument 'outfile'`; corrected the call (`ID=`), unblocking all G/OP runs on CG structures.
- `EntDetect/compare_sim2exp.py` — Workflow 3 consistency test aligned SASA/XP OP arrays to the MSM by **trajectory number** instead of MSM position, fixing a silent misalignment (and dropped trajectories) when the MSM mapping already excludes mirror trajectories.
- `EntDetect/statistics.py` — `MonteCarlo` now handles single-entry gene lists (`np.atleast_1d`).

## [1.2.0] - 2026-06-11

### Added
- `Documentation/workflow3_sim2exp.ipynb` — comprehensive Jupyter notebook for Workflow 3 (Sim-to-Experiment consistency testing). Covers experimental data loading, SASA/XP array collection, LiP-MS/XL-MS consistency test execution on full production ensemble (1000 trajectories), and result visualization. Includes prominent runtime warnings (1–3 hours) to set user expectations.
- `collect_jwalk_optimized.py` — optimized Jwalk collection script for standalone cluster execution; designed for submission as long-running SLURM job.
- `assets/slurm/scripts/run_workflow3_jwalk_collect.slurm` — SLURM job script for collecting Jwalk data overnight with 6-hour walltime.
- `scripts/test_workflow3_collection.py` — comprehensive test script for verifying SASA/XP files and testing CollectOP on full 1000-trajectory ensemble.
- Enhanced documentation with table of contents in Workflow 1, 2, and 3 markdown files for improved navigation.
- Documentation: `workflow3_collection_test_results.md` — detailed test results and performance analysis for OP collection (superseded by notebook).

### Fixed
- `EntDetect/order_params.py` — **NaN coordinate handling in SASA and XP methods**: Added robust detection and fallback mechanism for trajectory frames with NaN coordinates (corrupted frames). When NaN detected:
  - **SASA method**: Uses previous frame's per-residue SASA values with current frame's timestamp; preserves frame indexing
  - **XP method**: Copies previous frame's cross-link DataFrame and updates Frame column to current index; supports both sequential and parallel execution modes
  - Frames with NaN are no longer skipped; frame counter increments normally to maintain alignment with original trajectory indexing
  - Comprehensive logging reports which frames/atoms contain NaN and whether fallback was successful
  - All 1000 trajectories now complete SASA/XP collection without errors or gaps
- `EntDetect/order_params.py` — `CollectOP.collect_SASA()` now correctly converts units from nm² to Ų (×100) and properly handles trajectories with missing files (fills with NaN).
- `EntDetect/order_params.py` — `CollectOP.collect_Jwalk()` now correctly parses XP file column names and maps `SASD` → `'Jwalk'` in output dictionary.

### Changed
- **Workflow 2 documentation** (`workflow2_trajectory_analysis.md`):
  - Removed duplicate "Using the command-line interface" section (9b) that duplicated the comprehensive "Running folding pathway analysis as a single script" section
  - Renumbered former section 9c to 9b (Expected outputs)
  - Added "Running folding pathway analysis as a single script" to table of contents for complete TOC coverage
- **Workflow 1, 2, 3 documentation**: Enhanced table of contents with all "Running X as a single script" sections for improved discoverability and user navigation.
- **Code organization**: All test/debug scripts for NaN handling consolidated into single diagnostic test suite.

### Performance
- SASA collection: ~60 seconds for 1000 trajectories (✓ verified)
- Jwalk collection: ~14 seconds per trajectory, estimated 3–4 hours for full ensemble
- Consistency test: 1–3 hours on full production data (1000 trajectories × 335 frames each)

## [1.1.7] - 2026-05-15

### Added
- `Documentation/workflow2_trajectory_analysis.ipynb` — new interactive Jupyter notebook covering the full Workflow 2 pipeline (Q/G/K/SASA/XP order parameters, non-native entanglement clustering, MSM construction, folding pathway statistics). Includes inline figures: G vs Q free-energy/metastable-state map and JS divergence time series plot.
- `CalculateOP.G()` now accepts `chunk_frames` and `chunk_suffix` keyword arguments for chunked GE pickle output, matching production SLURM settings.

### Fixed
- `CalculateOP.XP()`: re-reading cached Jwalk output files failed with `AttributeError` when the saved file was 7-column tab-separated (written by a prior XP call) rather than the original 6-column space-padded format. Parser now uses `sep=r'\s+'` with `engine='python'` and `index_col=False` so both formats are read correctly.
- `CalculateOP.SASA()` and `CalculateOP.XP()` output paths are now namespaced per trajectory (`1ZMR_Traj{N}.SASA`, `XLresidue_pairs` prefixed with ID+Traj, Jwalk results subdirectory named `Jwalk_results_{ID}_Traj{N}/`) to avoid cross-trajectory collisions when running on multiple trajectories in the same output directory.

### Changed
- SLURM scripts (`assets/slurm/scripts/run_OP_traj*.slurm`): corrected output directory paths from `OP_demo/` → `OP/` and `OP_demo_AA/` → `OP_AA/` to match production layout.
- `EntDetect/clustering.py`: parallelised non-native entanglement clustering; fixed MSM logger output ordering.
- `EntDetect/gaussian_entanglement.py`: refactored GE combination logic to support chunked pickle merging.
- `Documentation/workflow2_trajectory_analysis.md`: corrected output-file schemas, column names, and parameter names throughout (e.g. `tarj_type_col` typo in `MSMStats`/`FoldingPathwayStats`).

## [1.1.4] - 2026-03-03

### Added
- macOS-friendly setup path via `environment-mac.yml` and README installation notes.
- New tutorial document with runnable examples using bundled `assets/`: `Documentation/tutorial_examples.md`.

### Changed
- `run_OP_on_simulation_traj.py` now accepts CLI flags to select CG vs all-atom mode and C-alpha vs heavy-atom contact definitions (defaults preserve prior behavior).

## [1.1.5] - 2026-03-03

### Fixed
- Removed tracked `__pycache__/*.pyc` artifacts from version control (bytecode is now ignored and no longer shipped).

## [1.1.6] - 2026-03-03

### Changed
- `run_nativeNCLE.py` now supports `--resolution {aa,cg}` and `--contacts {heavy,calpha}` for consistency with `run_OP_on_simulation_traj.py` (legacy `--cg`/`--Calpha` flags still work).
- Tutorial examples updated to use `--resolution/--contacts` consistently.

## [1.1.3] - 2026-03-03

### Added
- `run_OP_on_simulation_traj.py` now supports `--ent_detection_method` (mirrors `run_nativeNCLE.py`) and threads it through the OP/GE pipeline.

### Changed
- Trajectory entanglement CSV outputs (`Traj_GE/*_GE.csv`) now only report per-frame contacts that are also present in the same-run reference contact set (`Native_GE/*_GE.csv`).
- Packaging/install workflow improvements to ensure bundled resources are available via `importlib.resources` and to reduce friction for fresh installs.

## [1.1.2] - 2026-01-26

### Added
- **Configurable ENT Detection Methods** - Three distinct methods for determining entanglement status:
  - Method 1: GLN-based (any nonzero Gaussian Linking Number)
  - Method 2: TLN-based (any nonzero Topological Linking Number) - default
  - Method 3: Consensus (both GLN and TLN must agree on same terminus)
  - New `--ent_detection_method` command-line argument in `run_nativeNCLE.py`
  - New `determine_ent_status()` method in `GaussianEntanglement` class

- **Memory Optimization for AlphaFold Models**
  - Auto-enable C-alpha mode for AF structures to reduce memory usage
  - Prevents out-of-memory errors on large structures

- **Version Tracking**
  - Added `__version__` to `EntDetect/__init__.py`
  - Added CHANGELOG.md for release documentation

### Fixed
- **Crossing Residue Validation Bug** - Critical fix in `find_crossing()` method:
  - Original code used alternating pair comparison which missed consecutive crossings
  - Replaced with greedy algorithm ensuring all retained crossings are ≥10 residues apart
  - Applied fix to both N-terminal and C-terminal crossing validation

### Removed
- **Placeholder '?' Logic** - Removed no-longer-needed sentinel value system:
  - Removed `mark_absent_crossings()` method from `GaussianEntanglement`
  - Removed `replace('?', '-100000')` from clustering pipeline
  - Crossing strings now only contain actual identified residues or remain empty

### Changed
- Default `ent_detection_method` in `run_nativeNCLE.py` changed to 3 (consensus method)
- Crossing validation now uses greedy filtering instead of alternating pair comparison

### Technical Details
- **Files Modified**:
  - `EntDetect/gaussian_entanglement.py` - Added ENT detection methods, fixed crossing validation
  - `EntDetect/clustering.py` - Removed placeholder logic
  - `scripts/run_nativeNCLE.py` - Added CLI argument, AF model auto-optimization
  - `EntDetect/__init__.py` - Added version and package metadata
  - `setup.py` - Updated version to 0.2.0

## [0.1.0] - Initial Release

### Features
- Native Gaussian entanglement detection in protein structures
- Support for both experimental (PDB) and AlphaFold (AF) structures
- Coarse-grained and all-atom representation options
- Entanglement clustering to remove degeneracies
- High-quality entanglement filtering
- Entanglement feature generation
- Disulfide bond detection
- Support for C-alpha (CA) and coarse-grained (CG) models
