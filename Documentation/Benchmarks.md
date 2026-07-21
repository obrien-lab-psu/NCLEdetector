## Benchmark Wall Times

Wall times measured from `NORMAL TERMINATION` markers in the SLURM logs (`assets/slurm/logs/`) and the per-analysis application logs at the datastore output destinations (for example `/scratch/ims86/EntDetect_Datastore/outputs/workflow2/*/logs/`) on the PSU Roar Collab cluster (basic partition). Times reflect processing of a 387-residue protein (ecPGK, 1ZMR) with the parameters used in the tutorials, except Workflow 4 which is proteome-scale. The `minimal_workflow` column matches the numbered Minimal Workflow sections in the tutorial markdowns. Averages and ranges are based on successful runs (logs containing NORMAL TERMINATION).

| workflow | minimal_workflow | script | n_cpus | wall_time_avg | wall_time_range | notes |
|---|---|---|---|---|---|---|
| workflow1 | Minimal Workflow 1 | run_nativeNCLE.py | 1 | <1 min | <1 min | Datastore log shows 0.7 s for a resumed run on ecPGK (1ZMR); full de novo NCLE runtime depends on structure size and entanglement complexity. |
| workflow2 | Minimal Workflow 2 (CG full traj) | run_OP_on_simulation_traj.py | 10 | ~26 min | 14–37 min | Per-trajectory across 1000 CG trajectories (full DCD, ~6667 frames), ops=[Q,G,K]; config `workflow2_OP_config.json`. Stats use the per-file final run to exclude partial reruns. |
| workflow2 | Minimal Workflow 2 (CG last-335) | run_OP_on_simulation_traj.py | 10 | ~17 min | 3.5–30 min | Per-trajectory across 1000 CG trajectories (last 335 of 6667 frames), ops=[Q,G,K]; config `workflow2_OP_last335_config.json`. Variance driven by entanglement complexity per traj. |
| workflow2 | Minimal Workflow 2 (AA SASA/XP last-335) | run_OP_on_simulation_traj.py | 10 | ~43 min | 22–168 min | Per-trajectory across 1000 all-atom trajectories (last 335 frames), ops=[SASA,XP]; config `workflow2_OP_AA_last335_config.json`. Larger cost and spread than CG; stats use the per-file final run. |
| workflow2 | Minimal Workflow 3 | run_nonnative_entanglement_clustering.py | 8 | ~88 min | 62–114 min | 2 datastore runs: full trajectory set (61.8 min) and last-335-frame set (114.1 min); `nproc=8`. Runtime scales with number of detected entanglements, not frame count. |
| workflow2 | Minimal Workflow 4 | run_MSM.py | 8 | ~5.9 min | 5.3–6.8 min | 3 datastore runs (317–410 s); Q/G MSM construction. |
| workflow2 | Minimal Workflow 5 | run_MSMStats.py | 4 | ~2.8 min | 2.7–2.9 min | Full runs across QG_native and random conditions (163–175 s); excludes 1.5 s partial re-runs. |
| workflow2 | Minimal Workflow 6 | run_Foldingpathway.py | 4 | <1 min | <1 min | 3 datastore runs (14.6–16.7 s) across QG_native, random, and A80pctNative conditions. |
| workflow3 | Minimal Workflow 7 | run_compare_sim2exp.py | 8 | ~56 min | 56 min | 1 successful run of 15 attempts; other jobs failed precheck or were killed. Run includes collect + consistency test path. |
| workflow4 | Step 1 (MW8 preprocessing) | run_workflow4_nativeNCLE_batch.py | 20 | ~17 min | 17–18 min | 3 successful runs of 4 attempts; processed 544 proteins per run in AF model mode. Runs Workflow 1 in batch and builds the combined design matrix. |
| workflow4 | Minimal Workflow 8 | run_population_modeling.py | 8 | <1 min | <1 min | 2 successful runs of 3 attempts (one failed precheck due to missing input directory). Regression only; batch preprocessing timed separately above. |
| workflow4 | Minimal Workflow 9 | run_montecarlo.py | 20 | ~27.5 min | 27.4–27.6 min | Monte Carlo sampling for population modeling (`steps=10000`, `n_groups=4`). |
