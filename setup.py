from setuptools import setup, find_packages

with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='NCLEdetector',
    version='2.0.0',
    description='Non-covalent Lasso-like Entanglement (NCLE) Detection in Protein Structures and trajectories',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/obrien-lab-psu/NCLEdetector',
    project_urls={
        'Source': 'https://github.com/obrien-lab-psu/NCLEdetector',
        'Bug Tracker': 'https://github.com/obrien-lab-psu/NCLEdetector/issues',
    },
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'NCLEdetector': [
            'resources/calc_Q.pl',
            'resources/calc_K.pl',
            'resources/pulchra',
            'resources/stride',
            'resources/shared_files/*.dat',
            'Jwalk/naccess.config.txt',
        ]
    },
    install_requires=[
        'biopython',
        'numpy',
        'scipy',
        'pandas',
        'MDAnalysis',
        'mdtraj',
        'parmed',
        'numba',
        'topoly',
        'geom_median',
        'matplotlib',
        'seaborn',
        'scikit-learn',
        'networkx',
        'pyyaml',
        'tqdm',
        'requests',
        'statsmodels',
    ],
    entry_points={
        'console_scripts': [
            'convert_cor_psf_to_pdb=scripts.convert_cor_psf_to_pdb:main',
            'run_nativeNCLE=scripts.run_nativeNCLE:main',
            'run_OP_on_simulation_traj=scripts.run_OP_on_simulation_traj:main',
            'run_change_resolution=scripts.run_change_resolution:main',
            'run_nonnative_entanglement_clustering=scripts.run_nonnative_entanglement_clustering:main',
            'run_MSM=scripts.run_MSM:main',
            'run_compare_sim2exp=scripts.run_compare_sim2exp:main',
            'run_population_modeling=scripts.run_population_modeling:main',
            'run_montecarlo=scripts.run_montecarlo:main',
            'run_Foldingpathway=scripts.run_Foldingpathway:main',
        ],
    },
)