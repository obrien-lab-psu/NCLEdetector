#!/usr/bin/env python
"""
Generate labeled MSM mapping files for Case 1 (biologically-informed) and Case 2 (random).

This script implements Step 6a from workflow2_trajectory_analysis.md:
- Case 1: A = max Q >= 0.80 AND max G <= 0.05 (native-like, non-entangled); B = all others
- Case 2: A/B assigned randomly per trajectory (seed=42)
"""

import pandas as pd
import numpy as np
import sys

DATASTORE = "/scratch/ims86/EntDetect_Datastore"
OUTDIR = f"{DATASTORE}/outputs/workflow2"

# Load the base MSM mapping
print("Loading base MSMmapping...")
msm_mapping = pd.read_csv(f"{OUTDIR}/MSM/1ZMR_prod_MSMmapping.csv")
print(f"Loaded {len(msm_mapping)} rows")

# ── Case 1: Biologically-informed split (Q >= 0.8 and G <= 0.05) ──────────────────────
print("\n--- Generating Case 1: Biologically-informed split ---")
traj_stats = msm_mapping.groupby('traj').agg(max_Q=('Q', 'max'), max_G=('G', 'max'))
native_trajs = traj_stats[(traj_stats['max_Q'] >= 0.80) & (traj_stats['max_G'] <= 0.05)].index
msm_mapping['traj_type_QG_native'] = msm_mapping['traj'].isin(native_trajs).map({True: 'A', False: 'B'})

annotated_file_case1 = f"{OUTDIR}/MSM/1ZMR_prod_MSMmapping_QG_native.csv"
msm_mapping[['traj', 'frame', 'microstate', 'metastablestate', 'Q', 'G', 'traj_type_QG_native']].to_csv(
    annotated_file_case1, index=False)
print(f"Case 1 — trajectory type distribution:")
print(msm_mapping.groupby('traj')['traj_type_QG_native'].first().value_counts())
print(f"Saved to: {annotated_file_case1}")

# ── Case 2: Random split (negative control) ────────────────────────────────────────────
print("\n--- Generating Case 2: Random split (negative control) ---")
rng = np.random.default_rng(seed=42)
all_trajs = msm_mapping['traj'].unique()
random_labels = dict(zip(all_trajs, rng.choice(['A', 'B'], size=len(all_trajs))))
msm_mapping['traj_type_random'] = msm_mapping['traj'].map(random_labels)

annotated_file_case2 = f"{OUTDIR}/MSM/1ZMR_prod_MSMmapping_random.csv"
msm_mapping[['traj', 'frame', 'microstate', 'metastablestate', 'Q', 'G', 'traj_type_random']].to_csv(
    annotated_file_case2, index=False)
print(f"Case 2 — random type distribution:")
print(msm_mapping.groupby('traj')['traj_type_random'].first().value_counts())
print(f"Saved to: {annotated_file_case2}")

print("\n✓ Both labeled MSM files generated successfully")
