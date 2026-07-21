# NCLEdetector

A comprehensive Python package for studying non-covalent lasso entanglements in protein folding through molecular dynamics simulations and experimental data analysis.

## Overview

NCLEdetector provides a complete toolkit for analyzing protein entanglements across multiple scales - from individual structures to large-scale proteomic datasets. The package enables researchers to:

- **Identify and characterize** native entanglements in protein structures
- **Calculate order parameters** for simulation trajectories (Q, G, K, SASA, Jwalk, XP)
- **Build Markov State Models** from coarse-grained simulation ensembles
- **Compare simulations to experiments** using LiP-MS and XL-MS data
- **Perform population-level analysis** across heterogeneous protein datasets
- **Coarse-grain and back-map** between atomic resolutions

## Key Features

- **Multi-scale Analysis**: From single proteins to proteome-wide studies
- **Experimental Integration**: Direct comparison with mass spectrometry data
- **Advanced Clustering**: Non-native entanglement clustering and MSM construction
- **Statistical Methods**: Monte Carlo simulations and logistic regression modeling
- **Flexible Resolution**: Seamless conversion between all-atom and coarse-grained representations

## Installation

NCLEdetector supports two installation methods. Choose the one that best matches your workflow.

### Method 1: Conda environment install (Linux and macOS)

Use this when you want a local Python environment and direct script/package development.

From the repo root:

```bash
conda env create -f environment.yml
conda activate ncledetector
```

Notes:
- The provided conda environment targets Python 3.11 for best compatibility with the scientific stack.
- Run `conda env create` from the NCLEdetector repo root (the environment file uses `pip -e .`).
- If you prefer installing into an existing env, use `pip install -e .` from the repo root.

macOS notes:
- The default `environment.yml` is intended to work on both Linux and macOS via conda-forge.
- If you hit solver/build issues on macOS, try:

```bash
conda env create -f environment-mac.yml
conda activate ncledetector
```

- On Apple Silicon, if a dependency is missing for `osx-arm64`, a common workaround is:

```bash
CONDA_SUBDIR=osx-64 conda env create -f environment-mac.yml
conda activate ncledetector
```

### Method 2: Container install (Docker and Apptainer/Singularity)

Use this for reproducible runs, HPC workflows, or when you do not want to manage local dependencies.

Pull prebuilt images:

```bash
# Docker
docker pull ghcr.io/obrien-lab-psu/ncledetector:latest

# Apptainer/Singularity
apptainer pull ncledetector.sif docker://ghcr.io/obrien-lab-psu/ncledetector:latest
```

For full run examples and bind-mount guidance, see [Documentation/container_usage.md](Documentation/container_usage.md).

## Tutorials

Step-by-step, runnable tutorials covering all four analysis workflows are in `Documentation/`. Start here:

- [Documentation/index.md](Documentation/index.md) — master index with environment setup, path variables, and links to all workflows

For quick CLI reference, every script supports `--help`:

```bash
python scripts/run_nativeNCLE.py --help
python scripts/run_OP_on_simulation_traj.py --help
```

## Package Structure

```
NCLEdetector/
├── NCLEdetector/                    # Main package
│   ├── __init__.py
│   ├── gaussian_entanglement.py  # Core entanglement calculations
│   ├── clustering.py             # Entanglement clustering methods
│   ├── order_params.py          # Order parameter calculations
│   ├── compare_sim2exp.py       # Simulation-experiment comparison
│   ├── statistics.py            # Statistical analysis methods
│   ├── entanglement_features.py # Feature generation
│   ├── change_resolution.py     # Resolution conversion
│   └── utilities.py             # Helper functions
├── scripts/                     # Example workflow scripts
├── Documentation/               # Detailed module documentation
└── TestingGrounds/             # Test data and examples
```

## Core Modules

- **`gaussian_entanglement`**: Calculate Gaussian linking numbers and identify entanglements
- **`clustering`**: Cluster native and non-native entanglements, build MSMs
- **`order_params`**: Compute Q, G, K, SASA, Jwalk, and cross-linking propensity
- **`compare_sim2exp`**: Integrate LiP-MS and XL-MS experimental data
- **`statistics`**: Population modeling and Monte Carlo analysis
- **`entanglement_features`**: Generate structural features for entanglements
- **`change_resolution`**: Convert between all-atom and coarse-grained representations

## Documentation

Detailed documentation for each module:

- [Tutorial index](Documentation/index.md)
- [Workflow 1: Native NCLE detection](Documentation/workflow1_native_ncle.md)
- [Workflow 2: Trajectory analysis](Documentation/workflow2_trajectory_analysis.md)
- [Workflow 3: Sim-to-experiment comparison](Documentation/workflow3_sim2exp.md)
- [Workflow 4: Population-level analysis](Documentation/workflow4_population.md)
- [Container usage (Docker and Apptainer/Singularity)](Documentation/container_usage.md)
- [Gaussian Entanglement](Documentation/gaussian_entanglement.md)
- [Clustering](Documentation/clustering.md)
- [Order Parameters](Documentation/order_params.md)
- [Simulation-Experiment Comparison](Documentation/compare_sim2exp.md)
- [Statistical Analysis](Documentation/statistics.md)
- [Entanglement Features](Documentation/entanglement_features.md)
- [Resolution Conversion](Documentation/change_resolution.md)
- [Utilities](Documentation/utilities.md)

## Requirements

- Python 3.8+
- NumPy
- Pandas
- SciPy
- MDAnalysis
- OpenMM (for force field operations)
- Matplotlib (for visualization)
- See `environment.yml` for complete dependencies

## Citation

If you use NCLEdetector in your research, please cite:

```bibtex
@software{ncledetector2024,
  title={NCLEdetector: A Python Package for Protein Entanglement Analysis},
  author={Ian Sitarik},
  author={Yang Jiang},
  author={Hyebin Song},
  author={Edward O'Brien},
  year={2026},
  url={https://github.com/obrien-lab-psu/NCLEdetector}
}
```

## Contributing

Contributions are welcome! Please see our contributing guidelines and submit pull requests for any improvements.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Support

For questions and support, please open an issue on GitHub or contact the developers. 
  

