# Changelog

All notable changes to EntDetect will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
