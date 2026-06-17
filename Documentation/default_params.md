# Default Parameters Across Tutorials

This table tracks default parameters used in EntDetect tutorials and their implementation touchpoints.

| workflow | analysis_name | parameter | description | default_value | class/attribute/methods/functions that use this parameter across EntDetect |
|---|---|---|---|---|---|
| workflow1 | native_entanglement_detection | g_threshold | Gaussian entanglement threshold for NCLE detection/filtering | 0.6 | EntDetect.gaussian_entanglement.GaussianEntanglement.__init__; EntDetect.order_params.CalculateOP.Gpy (local g_threshold=0.6); scripts/run_nativeNCLE.py tutorial usage |
| workflow1 | native_entanglement_detection | density | Triangulation density for minimal loop surface | 0.0 | EntDetect.gaussian_entanglement.GaussianEntanglement.__init__ |
| workflow1 | native_entanglement_detection | ent_detection_method | ENT criterion (1=GLN, 2=TLN, 3=GLN+TLN same termini) | 2 (class default), 3 (run_nativeNCLE script default) | EntDetect.gaussian_entanglement.GaussianEntanglement.__init__/determine_ent_status; scripts/run_nativeNCLE.py; scripts/run_workflow4_nativeNCLE_batch.py |
| workflow1 | native_entanglement_detection | Calpha | Use C-alpha contact mode | False | EntDetect.gaussian_entanglement.GaussianEntanglement.__init__; scripts/run_nativeNCLE.py |
| workflow1 | native_entanglement_detection | CG | Input is coarse-grained C-alpha model | False | EntDetect.gaussian_entanglement.GaussianEntanglement.__init__; scripts/run_nativeNCLE.py |
| workflow1 | native_entanglement_detection | nproc | Worker count for frame-wise helper jobs | 10 | EntDetect.gaussian_entanglement.GaussianEntanglement.__init__ |
| workflow1 | native_clustering | organism | Organism-specific clustering cutoff preset | Ecoli | EntDetect.clustering.ClusterNativeEntanglements.__init__; scripts/run_nativeNCLE.py; scripts/run_workflow4_nativeNCLE_batch.py |
| workflow1 | native_clustering | cut_off | Native-clustering distance cutoff override | None (uses organism preset: Ecoli=57, Human=52, Yeast=49) | EntDetect.clustering.ClusterNativeEntanglements.__init__; scripts/run_nativeNCLE.py; scripts/run_workflow4_nativeNCLE_batch.py |
| workflow1 | native_hq_selection | model | Structure model type for HQ filtering | EXP | scripts/run_nativeNCLE.py (default EXP; tutorial often overrides AF for AF structures) |
| workflow1 | logging | log_level | Logging verbosity | INFO | scripts/run_nativeNCLE.py; EntDetect._logging.setup_logger consumers across modules |
| workflow2 | order_parameters | outdir | Base output directory for OP products | ./ | EntDetect.order_params.CalculateOP.__init__; scripts/run_OP_on_simulation_traj.py |
| workflow2 | order_parameters | Traj | Trajectory number tag | 1 | EntDetect.order_params.CalculateOP.__init__; scripts/run_OP_on_simulation_traj.py |
| workflow2 | order_parameters | start | Start frame index (0-based) | 0 | EntDetect.order_params.CalculateOP.__init__; scripts/run_OP_on_simulation_traj.py |
| workflow2 | order_parameters | end | End frame index (large sentinel = all frames) | 99999999999999 (class); script may pass None/-1 semantics externally | EntDetect.order_params.CalculateOP.__init__; scripts/run_OP_on_simulation_traj.py |
| workflow2 | order_parameters | stride | Frame stride for OP calculations | 1 | EntDetect.order_params.CalculateOP.__init__; scripts/run_OP_on_simulation_traj.py |
| workflow2 | order_parameters | ent_detection_method | ENT criterion for G calculations | 2 (CalculateOP class), 1 (run_OP_on_simulation_traj.py default) | EntDetect.order_params.CalculateOP.__init__/Gpy; scripts/run_OP_on_simulation_traj.py |
| workflow2 | order_parameters | ops | Selected OP outputs to compute | [Q, G, K] | scripts/run_OP_on_simulation_traj.py |
| workflow2 | order_parameters | resolution | Trajectory resolution mode | cg | scripts/run_OP_on_simulation_traj.py |
| workflow2 | order_parameters | contacts | Contact mode for topology operations | None (derived: calpha for cg, heavy for aa) | scripts/run_OP_on_simulation_traj.py |
| workflow2 | order_parameters | nproc | CPU cores used for G calculation path | 10 | scripts/run_OP_on_simulation_traj.py; EntDetect.gaussian_entanglement.GaussianEntanglement |
| workflow2 | order_parameters | chunk_frames | Chunk size for Combined_GE output splitting | None | scripts/run_OP_on_simulation_traj.py |
| workflow2 | order_parameters | chunk_suffix | Naming suffix for chunked Combined_GE outputs | _chunk | scripts/run_OP_on_simulation_traj.py |
| workflow2 | msm_building | start | First frame for MSM feature loading | 0 | scripts/run_MSM.py; EntDetect.clustering.MSMNonNativeEntanglementClustering.__init__ |
| workflow2 | msm_building | end | Last frame for MSM feature loading | 99999999999 | scripts/run_MSM.py; EntDetect.clustering.MSMNonNativeEntanglementClustering.__init__ |
| workflow2 | msm_building | stride | Frame stride for MSM feature loading | 1 | scripts/run_MSM.py; EntDetect.clustering.MSMNonNativeEntanglementClustering.__init__ |
| workflow2 | msm_building | ITS | Toggle implied-timescale analysis mode | False (script string), class default 'False' | scripts/run_MSM.py; EntDetect.clustering.MSMNonNativeEntanglementClustering.__init__/run |
| workflow2 | msm_building | lagtime | MSM lag time (frames) | 20 (script), 1 (class default) | scripts/run_MSM.py; EntDetect.clustering.MSMNonNativeEntanglementClustering.__init__/build_msm |
| workflow2 | msm_building | n_cluster | Number of k-means microstates | 400 | scripts/run_MSM.py; EntDetect.clustering.MSMNonNativeEntanglementClustering.__init__/kmeans_clustering |
| workflow2 | msm_building | kmean_stride | Frame stride during k-means clustering | 2 | scripts/run_MSM.py; EntDetect.clustering.MSMNonNativeEntanglementClustering.__init__ |
| workflow2 | msm_building | n_small_states | Number of inactive-state macro clusters | 1 | scripts/run_MSM.py; EntDetect.clustering.MSMNonNativeEntanglementClustering.__init__ |
| workflow2 | msm_building | n_large_states | Number of active metastable macro states | 10 | scripts/run_MSM.py; EntDetect.clustering.MSMNonNativeEntanglementClustering.__init__ |
| workflow2 | msm_building | dt | MD timestep in ns for axis scaling | 0.015/1000 | scripts/run_MSM.py; EntDetect.clustering.MSMNonNativeEntanglementClustering.__init__ |
| workflow2 | msm_building | rm_traj_list | Trajectories to exclude (e.g., mirror images) | [] | scripts/run_MSM.py; EntDetect.clustering.MSMNonNativeEntanglementClustering.__init__/load_OP |
| workflow2 | nonnative_clustering | start_frame | First frame in nonnative entanglement clustering | 0 | scripts/run_nonnative_entanglement_clustering.py |
| workflow2 | nonnative_clustering | end_frame | Last frame in nonnative entanglement clustering | 9999999 | scripts/run_nonnative_entanglement_clustering.py |
| workflow2 | nonnative_clustering | nproc | Worker threads for nonnative clustering | 1 | scripts/run_nonnative_entanglement_clustering.py; EntDetect.clustering.ClusterNonNativeEntanglements.__init__ |
| workflow2 | pathway_stats | traj_type_list | Trajectory labels compared in pathway statistics | [A, B] | scripts/run_Foldingpathway.py; EntDetect.statistics.FoldingPathwayStats |
| workflow2 | pathway_stats | rm_traj_list | Trajectory exclusions for pathway analyses | [] | scripts/run_Foldingpathway.py; EntDetect.statistics.FoldingPathwayStats |
| workflow2 | pathway_stats | n_window | Rolling window size for state probability smoothing | 200 | scripts/run_Foldingpathway.py; EntDetect.statistics.FoldingPathwayStats |
| workflow2 | pathway_stats | n_traj | Total number of trajectories for normalization | 1000 | scripts/run_Foldingpathway.py; EntDetect.statistics.FoldingPathwayStats |
| workflow2 | pathway_stats | state_type | State column used for pathway stats | metastablestate | scripts/run_Foldingpathway.py; EntDetect.statistics.FoldingPathwayStats |
| workflow3 | consistency_test | start | First frame index used when loading trajectory-level arrays | 0 (script usage), 0 (class default) | scripts/run_compare_sim2exp.py; EntDetect.compare_sim2exp.MassSpec.__init__ |
| workflow3 | consistency_test | end | Last frame index used for loading | -1 (script usage), 999999999999 (class default) | scripts/run_compare_sim2exp.py; EntDetect.compare_sim2exp.MassSpec.__init__ |
| workflow3 | consistency_test | stride | Frame stride during loading | 1 | scripts/run_compare_sim2exp.py; EntDetect.compare_sim2exp.MassSpec.__init__ |
| workflow3 | consistency_test | verbose | Verbose mode flag | False | EntDetect.compare_sim2exp.MassSpec.__init__; scripts/run_compare_sim2exp.py |
| workflow3 | consistency_test | num_perm | Number of permutations in significance testing | 10000 (class), script usually supplies 1000 | EntDetect.compare_sim2exp.MassSpec.__init__/permutation_test/LiP_XL_MS_ConsistencyTest; scripts/run_compare_sim2exp.py |
| workflow3 | consistency_test | n_boot | Number of bootstrap samples | 10000 (class), script usually supplies 100 | EntDetect.compare_sim2exp.MassSpec.__init__/bootstrap usage in LiP_XL_MS_ConsistencyTest; scripts/run_compare_sim2exp.py |
| workflow3 | consistency_test | lag_frame | Downsampling lag (frames) for state sampling | 1 (class), script usually supplies 20 | EntDetect.compare_sim2exp.MassSpec.__init__/LiP_XL_MS_ConsistencyTest; scripts/run_compare_sim2exp.py |
| workflow3 | consistency_test | nproc | Parallel worker count for downstream representative-structure extraction | 1 (class), script passes user value | EntDetect.compare_sim2exp.MassSpec.__init__/select_rep_structs; scripts/run_compare_sim2exp.py |
| workflow3 | consistency_test | rm_traj_list | Trajectory exclusions (mirror images) | [] (class) | EntDetect.compare_sim2exp.MassSpec.__init__/all selection/filtering paths; scripts/run_compare_sim2exp.py |
| workflow3 | consistency_test | resid2residueidx_map | Residue index remapping dictionary | {} (identity map inferred when empty) | EntDetect.compare_sim2exp.MassSpec.__init__; select_rep_structs mapping logic |
| workflow3 | consistency_test_collection | collect_jwalk_npy | Build legacy Jwalk.npy from XP files | False (flag off by default) | scripts/run_compare_sim2exp.py; EntDetect.order_params.CollectOP.collect_Jwalk |
| workflow4 | native_ncle_batch | nproc | Parallel run_nativeNCLE jobs for batch wrapper | 8 | scripts/run_workflow4_nativeNCLE_batch.py |
| workflow4 | native_ncle_batch | allow_prefix_match | Allow gene IDs to match PDB stem by prefix token | False | scripts/run_workflow4_nativeNCLE_batch.py |
| workflow4 | native_ncle_batch | dry_run | Print selected proteins only, no execution | False | scripts/run_workflow4_nativeNCLE_batch.py |
| workflow4 | native_ncle_batch | organism | Organism mode forwarded to run_nativeNCLE | Ecoli | scripts/run_workflow4_nativeNCLE_batch.py; scripts/run_nativeNCLE.py; EntDetect.clustering.ClusterNativeEntanglements |
| workflow4 | native_ncle_batch | Accession | Accession identifier forwarded to feature generation | P00558 | scripts/run_workflow4_nativeNCLE_batch.py; scripts/run_nativeNCLE.py; EntDetect.entanglement_features.FeatureGen.get_uent_features |
| workflow4 | native_ncle_batch | resolution | Resolution forwarded to run_nativeNCLE | None | scripts/run_workflow4_nativeNCLE_batch.py; scripts/run_nativeNCLE.py |
| workflow4 | native_ncle_batch | contacts | Contact mode forwarded to run_nativeNCLE | None | scripts/run_workflow4_nativeNCLE_batch.py; scripts/run_nativeNCLE.py |
| workflow4 | native_ncle_batch | cluster_cutoff | Native clustering cutoff override in batch path | None | scripts/run_workflow4_nativeNCLE_batch.py; scripts/run_nativeNCLE.py; EntDetect.clustering.ClusterNativeEntanglements |
| workflow4 | native_ncle_batch | model | Model type forwarded to HQ filter stage | AF (batch wrapper), EXP (run_nativeNCLE default) | scripts/run_workflow4_nativeNCLE_batch.py; scripts/run_nativeNCLE.py; EntDetect.gaussian_entanglement.select_high_quality_entanglements |
| workflow4 | native_ncle_batch | ent_detection_method | ENT criterion forwarded to run_nativeNCLE | 3 (batch wrapper), 3 (run_nativeNCLE) | scripts/run_workflow4_nativeNCLE_batch.py; scripts/run_nativeNCLE.py; EntDetect.gaussian_entanglement.GaussianEntanglement |
| workflow4 | proteome_regression | reg_formula | Logistic model formula | cut_C_Rall ~ AA + region | scripts/run_population_modeling.py; EntDetect.statistics.ProteomeLogisticRegression.__init__/run |
| workflow4 | proteome_regression | log_level | Logging verbosity | INFO | scripts/run_population_modeling.py; EntDetect.statistics.ProteomeLogisticRegression |
| workflow4 | monte_carlo | reg_formula | Regression formula used during state evaluation | cut_C_Rall ~ region + AA | scripts/run_montecarlo.py; EntDetect.statistics.MonteCarlo.__init__/run |
| workflow4 | monte_carlo | response_var | Response variable for model fitting | cut_C_Rall | scripts/run_montecarlo.py; EntDetect.statistics.MonteCarlo.__init__/load_data |
| workflow4 | monte_carlo | test_var | Coefficient tested for enrichment objective | region | scripts/run_montecarlo.py; EntDetect.statistics.MonteCarlo.__init__/run |
| workflow4 | monte_carlo | random | Enable random baseline mode | False | scripts/run_montecarlo.py; EntDetect.statistics.MonteCarlo.__init__/run |
| workflow4 | monte_carlo | n_groups | Number of subpopulations | 4 | scripts/run_montecarlo.py; EntDetect.statistics.MonteCarlo.__init__/create_unique_subgroups |
| workflow4 | monte_carlo | steps | Number of Monte Carlo optimization steps | 100000 | scripts/run_montecarlo.py; EntDetect.statistics.MonteCarlo.__init__/run |
| workflow4 | monte_carlo | C1 | Objective coefficient on enrichment term | 1.0 | scripts/run_montecarlo.py; EntDetect.statistics.MonteCarlo.__init__/run |
| workflow4 | monte_carlo | C2 | Objective coefficient on size-distribution penalty | 2.5 | scripts/run_montecarlo.py; EntDetect.statistics.MonteCarlo.__init__/run |
| workflow4 | monte_carlo | beta | Initial inverse temperature for annealing | 0.05 | scripts/run_montecarlo.py; EntDetect.statistics.MonteCarlo.__init__/run |
| workflow4 | monte_carlo | linearT | Use linear temperature schedule | False | scripts/run_montecarlo.py; EntDetect.statistics.MonteCarlo.__init__/run |
| workflow4 | monte_carlo | log_level | Logging verbosity | INFO | scripts/run_montecarlo.py; EntDetect.statistics.MonteCarlo |

## Notes

- This table focuses on defaults that are explicit in tutorial-facing scripts/classes used by Workflows 1–4.
- Some classes require caller-supplied values with no intrinsic defaults (for example, `MonteCarlo` constructor arguments). In those cases, the tutorial script defaults are listed.
- Where script defaults differ from class defaults, both are noted in `default_value`.
