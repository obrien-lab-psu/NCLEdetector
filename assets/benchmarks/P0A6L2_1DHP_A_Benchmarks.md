## Benchmark Wall Times - P0A6L2_1DHP_A

Template benchmark table for the CG benchmark test protein P0A6L2 (1DHP, chain A). Fields are aligned with Documentation/Benchmarks.md. Fill wall times from successful NORMAL TERMINATION runs and application logs as needed.

| workflow | minimal_workflow | script | n_cpus | wall_time_avg | wall_time_range | notes |
|---|---|---|---|---|---|---|
| workflow1 | Minimal Workflow 1 | run_nativeNCLE.py | 1 | 1.8 min | 1.8 min-1.8 min | Completed in rerun job 54205932 (task 3). |
| workflow2 | Minimal Workflow 2 (CG full traj) | run_OP_on_simulation_traj.py | 8 | 8.4 min | 7.1 min-11 min | Completed for 50 trajectories (Q and chunked G outputs present). |
| workflow2 | Minimal Workflow 2 (CG last-67) | run_OP_on_simulation_traj.py | 8 | 1.1 min | 0.7 min-2.9 min | First-pass de novo estimate (base job 54185781 only); 50 completed trajectories in first pass. |
| workflow2 | Minimal Workflow 2 (AA SASA/XP last-67) | run_OP_on_simulation_traj.py |  |  |  |  |
| workflow2 | Minimal Workflow 3 | run_nonnative_entanglement_clustering.py | 4 |  |  | Not included in this benchmark summary (MW3 dropped). |
| workflow2 | Minimal Workflow 4 | run_MSM.py | 8 | 0.5 min | 0.5 min-0.5 min | Completed; MSM mapping and metastable outputs present. |
| workflow2 | Minimal Workflow 5 | run_MSMStats.py |  |  |  |  |
| workflow2 | Minimal Workflow 6 | run_Foldingpathway.py |  |  |  |  |
| workflow3 | Minimal Workflow 7 | run_compare_sim2exp.py |  |  |  |  |
| workflow4 | Step 1 (MW8 preprocessing) | run_workflow4_nativeNCLE_batch.py |  |  |  |  |
| workflow4 | Minimal Workflow 8 | run_population_modeling.py |  |  |  |  |
| workflow4 | Minimal Workflow 9 | run_montecarlo.py |  |  |  |  |
