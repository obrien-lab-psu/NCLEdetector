# Default Parameters Across Tutorials

This table tracks default parameters used in EntDetect tutorials and their implementation touchpoints.

| workflow | analysis_name | workflow_steps_covered | parameter | description | default_value | class/attribute/methods/functions that use this parameter across EntDetect |
|---|---|---|---|---|---|---|
| workflow1 | run_nativeNCLE.py | workflow1: steps 3-6 | g_threshold | Gaussian entanglement threshold for NCLE detection/filtering | 0.6 | EntDetect.gaussian_entanglement.GaussianEntanglement.__init__; EntDetect.order_params.CalculateOP.Gpy (local g_threshold=0.6); scripts/run_nativeNCLE.py tutorial usage |
| workflow1 | run_nativeNCLE.py | workflow1: steps 3-6 | density | Triangulation density for minimal loop surface | 1.0 | scripts/run_nativeNCLE.py; EntDetect.gaussian_entanglement.GaussianEntanglement.__init__ |
| workflow1 | run_nativeNCLE.py | workflow1: steps 3-6 | ent_detection_method | ENT criterion (1=GLN, 2=TLN, 3=GLN+TLN same termini) | 2 | EntDetect.gaussian_entanglement.GaussianEntanglement.__init__/determine_ent_status; scripts/run_nativeNCLE.py; scripts/run_workflow4_nativeNCLE_batch.py |
| workflow1 | run_nativeNCLE.py | workflow1: steps 3-6 | Calpha | Use C-alpha contact mode | False | EntDetect.gaussian_entanglement.GaussianEntanglement.__init__; scripts/run_nativeNCLE.py |
| workflow1 | run_nativeNCLE.py | workflow1: steps 3-6 | CG | Input is coarse-grained C-alpha model | False | EntDetect.gaussian_entanglement.GaussianEntanglement.__init__; scripts/run_nativeNCLE.py |
| workflow1 | run_nativeNCLE.py | workflow1: steps 3-6 | nproc | Worker count for frame-wise helper jobs | 10 | EntDetect.gaussian_entanglement.GaussianEntanglement.__init__ |
| workflow1 | run_nativeNCLE.py | workflow1: steps 3-6 | organism | Organism-specific clustering cutoff preset | Ecoli | EntDetect.clustering.ClusterNativeEntanglements.__init__; scripts/run_nativeNCLE.py; scripts/run_workflow4_nativeNCLE_batch.py |
| workflow1 | run_nativeNCLE.py | workflow1: steps 3-6 | cut_off | Native-clustering distance cutoff override | None (uses organism preset: Ecoli=57, Human=52, Yeast=49) | EntDetect.clustering.ClusterNativeEntanglements.__init__; scripts/run_nativeNCLE.py; scripts/run_workflow4_nativeNCLE_batch.py |
| workflow1 | run_nativeNCLE.py | workflow1: steps 3-6 | model | Structure model type for HQ filtering | EXP | scripts/run_nativeNCLE.py (default EXP; tutorial often overrides AF for AF structures) |
| workflow1 | run_nativeNCLE.py | workflow1: steps 3-6 | log_level | Logging verbosity | INFO | scripts/run_nativeNCLE.py; EntDetect._logging.setup_logger consumers across modules |
| workflow2 | run_OP_on_simulation_traj.py | workflow2: steps 2-3 | outdir | Base output directory for OP products | ./ | EntDetect.order_params.CalculateOP.__init__; scripts/run_OP_on_simulation_traj.py |
| workflow2 | run_OP_on_simulation_traj.py | workflow2: steps 2-3 | Traj | Trajectory number tag | required in run_OP_on_simulation_traj.py (fallback class default: 1) | scripts/run_OP_on_simulation_traj.py; EntDetect.order_params.CalculateOP.__init__ |
| workflow2 | run_OP_on_simulation_traj.py | workflow2: steps 2-3 | start | Start frame index (0-based) | 0 | EntDetect.order_params.CalculateOP.__init__; scripts/run_OP_on_simulation_traj.py |
| workflow2 | run_OP_on_simulation_traj.py | workflow2: steps 2-3 | end | End frame index (large sentinel = all frames) | 99999999999999 (class default; run_OP_on_simulation_traj.py does not expose `--end`) | EntDetect.order_params.CalculateOP.__init__; scripts/run_OP_on_simulation_traj.py |
| workflow2 | run_OP_on_simulation_traj.py | workflow2: steps 2-3 | stride | Frame stride for OP calculations | 1 | EntDetect.order_params.CalculateOP.__init__; scripts/run_OP_on_simulation_traj.py |
| workflow2 | run_OP_on_simulation_traj.py | workflow2: steps 2-3 | ent_detection_method | ENT criterion for G calculations | 2 (CalculateOP class), 1 (run_OP_on_simulation_traj.py default) | EntDetect.order_params.CalculateOP.__init__/Gpy; scripts/run_OP_on_simulation_traj.py |
| workflow2 | run_OP_on_simulation_traj.py | workflow2: steps 2-3 | ops | Selected OP outputs to compute | [Q, G, K] | scripts/run_OP_on_simulation_traj.py |
| workflow2 | run_OP_on_simulation_traj.py | workflow2: steps 2-3 | resolution | Trajectory resolution mode | cg | scripts/run_OP_on_simulation_traj.py |
| workflow2 | run_OP_on_simulation_traj.py | workflow2: steps 2-3 | contacts | Contact mode for topology operations | None (derived: calpha for cg, heavy for aa) | scripts/run_OP_on_simulation_traj.py |
| workflow2 | run_OP_on_simulation_traj.py | workflow2: steps 2-3 | nproc | CPU cores used for G calculation path | 10 | scripts/run_OP_on_simulation_traj.py; EntDetect.gaussian_entanglement.GaussianEntanglement |
| workflow2 | run_OP_on_simulation_traj.py | workflow2: steps 2-3 | chunk_frames | Chunk size for Combined_GE output splitting | None | scripts/run_OP_on_simulation_traj.py |
| workflow2 | run_OP_on_simulation_traj.py | workflow2: steps 2-3 | chunk_suffix | Naming suffix for chunked Combined_GE outputs | _chunk | scripts/run_OP_on_simulation_traj.py |
| workflow2 | run_MSM.py | workflow2: step 6 | start | First frame for MSM feature loading | 0 | scripts/run_MSM.py; EntDetect.clustering.MSMNonNativeEntanglementClustering.__init__ |
| workflow2 | run_MSM.py | workflow2: step 6 | end | Last frame for MSM feature loading | 99999999999 | scripts/run_MSM.py; EntDetect.clustering.MSMNonNativeEntanglementClustering.__init__ |
| workflow2 | run_MSM.py | workflow2: step 6 | stride | Frame stride for MSM feature loading | 1 | scripts/run_MSM.py; EntDetect.clustering.MSMNonNativeEntanglementClustering.__init__ |
| workflow2 | run_MSM.py | workflow2: step 6 | ITS | Toggle implied-timescale analysis mode | False (script string), class default 'False' | scripts/run_MSM.py; EntDetect.clustering.MSMNonNativeEntanglementClustering.__init__/run |
| workflow2 | run_MSM.py | workflow2: step 6 | lagtime | MSM lag time (frames) | 20 (script), 1 (class default) | scripts/run_MSM.py; EntDetect.clustering.MSMNonNativeEntanglementClustering.__init__/build_msm |
| workflow2 | run_MSM.py | workflow2: step 6 | n_cluster | Number of k-means microstates | 400 | scripts/run_MSM.py; EntDetect.clustering.MSMNonNativeEntanglementClustering.__init__/kmeans_clustering |
| workflow2 | run_MSM.py | workflow2: step 6 | kmean_stride | Frame stride during k-means clustering | 2 | scripts/run_MSM.py; EntDetect.clustering.MSMNonNativeEntanglementClustering.__init__ |
| workflow2 | run_MSM.py | workflow2: step 6 | n_small_states | Number of inactive-state macro clusters | 1 | scripts/run_MSM.py; EntDetect.clustering.MSMNonNativeEntanglementClustering.__init__ |
| workflow2 | run_MSM.py | workflow2: step 6 | n_large_states | Number of active metastable macro states | 10 | scripts/run_MSM.py; EntDetect.clustering.MSMNonNativeEntanglementClustering.__init__ |
| workflow2 | run_MSM.py | workflow2: step 6 | dt | MD timestep in ns for axis scaling | 0.015/1000 | scripts/run_MSM.py; EntDetect.clustering.MSMNonNativeEntanglementClustering.__init__ |
| workflow2 | run_MSM.py | workflow2: step 6 | rm_traj_list | Trajectories to exclude (e.g., mirror images) | [] | scripts/run_MSM.py; EntDetect.clustering.MSMNonNativeEntanglementClustering.__init__/load_OP |
| workflow2 | run_nonnative_entanglement_clustering.py | workflow2: step 5 | start_frame | First frame in nonnative entanglement clustering | 0 | scripts/run_nonnative_entanglement_clustering.py |
| workflow2 | run_nonnative_entanglement_clustering.py | workflow2: step 5 | end_frame | Last frame in nonnative entanglement clustering | 9999999 | scripts/run_nonnative_entanglement_clustering.py |
| workflow2 | run_nonnative_entanglement_clustering.py | workflow2: step 5 | nproc | Worker threads for nonnative clustering | 1 | scripts/run_nonnative_entanglement_clustering.py; EntDetect.clustering.ClusterNonNativeEntanglements.__init__ |
| workflow2 | run_Foldingpathway.py | workflow2: step 9 | traj_type_list | Trajectory labels compared in pathway statistics | [A, B] | scripts/run_Foldingpathway.py; EntDetect.statistics.FoldingPathwayStats |
| workflow2 | run_Foldingpathway.py | workflow2: step 9 | rm_traj_list | Trajectory exclusions for pathway analyses | [] | scripts/run_Foldingpathway.py; EntDetect.statistics.FoldingPathwayStats |
| workflow2 | run_Foldingpathway.py | workflow2: step 9 | n_window | Rolling window size for state probability smoothing | 200 | scripts/run_Foldingpathway.py; EntDetect.statistics.FoldingPathwayStats |
| workflow2 | run_Foldingpathway.py | workflow2: step 9 | n_traj | Total number of trajectories for normalization | 1000 | scripts/run_Foldingpathway.py; EntDetect.statistics.FoldingPathwayStats |
| workflow2 | run_Foldingpathway.py | workflow2: step 9 | state_type | State column used for pathway stats | metastablestate | scripts/run_Foldingpathway.py; EntDetect.statistics.FoldingPathwayStats |
| workflow3 | run_compare_sim2exp.py | workflow3: steps 4-5 | start | First frame index used when loading trajectory-level arrays | required in run_compare_sim2exp.py (fallback class default: 0) | scripts/run_compare_sim2exp.py; EntDetect.compare_sim2exp.MassSpec.__init__ |
| workflow3 | run_compare_sim2exp.py | workflow3: steps 4-5 | end | Last frame index used for loading | required in run_compare_sim2exp.py (fallback class default: 999999999999) | scripts/run_compare_sim2exp.py; EntDetect.compare_sim2exp.MassSpec.__init__ |
| workflow3 | run_compare_sim2exp.py | workflow3: steps 4-5 | stride | Frame stride during loading | required in run_compare_sim2exp.py (fallback class default: 1) | scripts/run_compare_sim2exp.py; EntDetect.compare_sim2exp.MassSpec.__init__ |
| workflow3 | run_compare_sim2exp.py | workflow3: steps 4-5 | verbose | Verbose mode flag | False | EntDetect.compare_sim2exp.MassSpec.__init__; scripts/run_compare_sim2exp.py |
| workflow3 | run_compare_sim2exp.py | workflow3: steps 4-5 | num_perm | Number of permutations in significance testing | required in run_compare_sim2exp.py (fallback class default: 10000) | EntDetect.compare_sim2exp.MassSpec.__init__/permutation_test/LiP_XL_MS_ConsistencyTest; scripts/run_compare_sim2exp.py |
| workflow3 | run_compare_sim2exp.py | workflow3: steps 4-5 | n_boot | Number of bootstrap samples | required in run_compare_sim2exp.py (fallback class default: 10000) | EntDetect.compare_sim2exp.MassSpec.__init__/bootstrap usage in LiP_XL_MS_ConsistencyTest; scripts/run_compare_sim2exp.py |
| workflow3 | run_compare_sim2exp.py | workflow3: steps 4-5 | lag_frame | Downsampling lag (frames) for state sampling | required in run_compare_sim2exp.py (fallback class default: 1) | EntDetect.compare_sim2exp.MassSpec.__init__/LiP_XL_MS_ConsistencyTest; scripts/run_compare_sim2exp.py |
| workflow3 | run_compare_sim2exp.py | workflow3: steps 4-5 | nproc | Parallel worker count for downstream representative-structure extraction | required in run_compare_sim2exp.py (fallback class default: 1) | EntDetect.compare_sim2exp.MassSpec.__init__/select_rep_structs; scripts/run_compare_sim2exp.py |
| workflow3 | run_compare_sim2exp.py | workflow3: steps 4-5 | rm_traj_list | Trajectory exclusions (mirror images) | [] (class) | EntDetect.compare_sim2exp.MassSpec.__init__/all selection/filtering paths; scripts/run_compare_sim2exp.py |
| workflow3 | run_compare_sim2exp.py | workflow3: steps 4-5 | resid2residueidx_map | Residue index remapping dictionary | {} (identity map inferred when empty) | EntDetect.compare_sim2exp.MassSpec.__init__; select_rep_structs mapping logic |
| workflow3 | run_compare_sim2exp.py | workflow3: steps 4-5 | collect_jwalk_npy | Build legacy Jwalk.npy from XP files | False (flag off by default) | scripts/run_compare_sim2exp.py; EntDetect.order_params.CollectOP.collect_Jwalk |
| workflow4 | run_workflow4_nativeNCLE_batch.py | workflow4: steps 3-4 | nproc | Parallel run_nativeNCLE jobs for batch wrapper | 8 | scripts/run_workflow4_nativeNCLE_batch.py |
| workflow4 | run_workflow4_nativeNCLE_batch.py | workflow4: steps 3-4 | allow_prefix_match | Allow gene IDs to match PDB stem by prefix token | False | scripts/run_workflow4_nativeNCLE_batch.py |
| workflow4 | run_workflow4_nativeNCLE_batch.py | workflow4: steps 3-4 | dry_run | Print selected proteins only, no execution | False | scripts/run_workflow4_nativeNCLE_batch.py |
| workflow4 | run_workflow4_nativeNCLE_batch.py | workflow4: steps 3-4 | organism | Organism mode forwarded to run_nativeNCLE | Ecoli | scripts/run_workflow4_nativeNCLE_batch.py; scripts/run_nativeNCLE.py; EntDetect.clustering.ClusterNativeEntanglements |
| workflow4 | run_workflow4_nativeNCLE_batch.py | workflow4: steps 3-4 | Accession | Accession identifier forwarded to feature generation | None in run_workflow4_nativeNCLE_batch.py (resolved per-protein to `ID`); P00558 in run_nativeNCLE.py | scripts/run_workflow4_nativeNCLE_batch.py; scripts/run_nativeNCLE.py; EntDetect.entanglement_features.FeatureGen.get_uent_features |
| workflow4 | run_workflow4_nativeNCLE_batch.py | workflow4: steps 3-4 | resolution | Resolution forwarded to run_nativeNCLE | None | scripts/run_workflow4_nativeNCLE_batch.py; scripts/run_nativeNCLE.py |
| workflow4 | run_workflow4_nativeNCLE_batch.py | workflow4: steps 3-4 | contacts | Contact mode forwarded to run_nativeNCLE | None | scripts/run_workflow4_nativeNCLE_batch.py; scripts/run_nativeNCLE.py |
| workflow4 | run_workflow4_nativeNCLE_batch.py | workflow4: steps 3-4 | cluster_cutoff | Native clustering cutoff override in batch path | None | scripts/run_workflow4_nativeNCLE_batch.py; scripts/run_nativeNCLE.py; EntDetect.clustering.ClusterNativeEntanglements |
| workflow4 | run_workflow4_nativeNCLE_batch.py | workflow4: steps 3-4 | model | Model type forwarded to HQ filter stage | AF (batch wrapper), EXP (run_nativeNCLE default) | scripts/run_workflow4_nativeNCLE_batch.py; scripts/run_nativeNCLE.py; EntDetect.gaussian_entanglement.select_high_quality_entanglements |
| workflow4 | run_workflow4_nativeNCLE_batch.py | workflow4: steps 3-4 | ent_detection_method | ENT criterion forwarded to run_nativeNCLE | 3 (batch wrapper), 3 (run_nativeNCLE) | scripts/run_workflow4_nativeNCLE_batch.py; scripts/run_nativeNCLE.py; EntDetect.gaussian_entanglement.GaussianEntanglement |
| workflow4 | run_population_modeling.py | workflow4: step 5 | reg_formula | Logistic model formula | cut_C_Rall ~ AA + region | scripts/run_population_modeling.py; EntDetect.statistics.ProteomeLogisticRegression.__init__/run |
| workflow4 | run_population_modeling.py | workflow4: step 5 | log_level | Logging verbosity | INFO | scripts/run_population_modeling.py; EntDetect.statistics.ProteomeLogisticRegression |
| workflow4 | run_montecarlo.py | workflow4: step 6 | reg_formula | Regression formula used during state evaluation | cut_C_Rall ~ region + AA | scripts/run_montecarlo.py; EntDetect.statistics.MonteCarlo.__init__/run |
| workflow4 | run_montecarlo.py | workflow4: step 6 | response_var | Response variable for model fitting | cut_C_Rall | scripts/run_montecarlo.py; EntDetect.statistics.MonteCarlo.__init__/load_data |
| workflow4 | run_montecarlo.py | workflow4: step 6 | test_var | Coefficient tested for enrichment objective | region | scripts/run_montecarlo.py; EntDetect.statistics.MonteCarlo.__init__/run |
| workflow4 | run_montecarlo.py | workflow4: step 6 | random | Enable random baseline mode | False | scripts/run_montecarlo.py; EntDetect.statistics.MonteCarlo.__init__/run |
| workflow4 | run_montecarlo.py | workflow4: step 6 | n_groups | Number of subpopulations | 4 | scripts/run_montecarlo.py; EntDetect.statistics.MonteCarlo.__init__/create_unique_subgroups |
| workflow4 | run_montecarlo.py | workflow4: step 6 | steps | Number of Monte Carlo optimization steps | 100000 | scripts/run_montecarlo.py; EntDetect.statistics.MonteCarlo.__init__/run |
| workflow4 | run_montecarlo.py | workflow4: step 6 | C1 | Objective coefficient on enrichment term | 1.0 | scripts/run_montecarlo.py; EntDetect.statistics.MonteCarlo.__init__/run |
| workflow4 | run_montecarlo.py | workflow4: step 6 | C2 | Objective coefficient on size-distribution penalty | 2.5 | scripts/run_montecarlo.py; EntDetect.statistics.MonteCarlo.__init__/run |
| workflow4 | run_montecarlo.py | workflow4: step 6 | beta | Initial inverse temperature for annealing | 0.05 | scripts/run_montecarlo.py; EntDetect.statistics.MonteCarlo.__init__/run |
| workflow4 | run_montecarlo.py | workflow4: step 6 | linearT | Use linear temperature schedule | False | scripts/run_montecarlo.py; EntDetect.statistics.MonteCarlo.__init__/run |
| workflow4 | run_montecarlo.py | workflow4: step 6 | log_level | Logging verbosity | INFO | scripts/run_montecarlo.py; EntDetect.statistics.MonteCarlo |

## Notes

- `analysis_name` is now the tutorial-facing run script filename (for example, `run_nativeNCLE.py`) rather than a conceptual analysis label.
- Default precedence policy used here: `scripts/run_*.py` defaults first, class/function constructor defaults second when scripts do not provide defaults.
- If a run script marks a parameter as required (no script default), the table lists it as `required in <script>` and then reports constructor fallback defaults when available.
- Where script defaults differ from class defaults, both are noted in `default_value`.
