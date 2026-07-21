from NCLEdetector.gaussian_entanglement import GaussianEntanglement
from NCLEdetector.clustering import ClusterNativeEntanglements, MSMNonNativeEntanglementClustering
from NCLEdetector.compare_sim2exp import MassSpec

if __package__ in {None, ""}:
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from scripts._cli_config import parse_args_with_config
else:
    from ._cli_config import parse_args_with_config

"""
Run the LiP-MS / XL-MS consistency test with automatic collection of SASA/XP files.

Provide parameters directly as CLI flags, through `--config` JSON/YAML, or both.
When both are provided, CLI flags override config values.

Two usage modes:

  1. Collection mode (provide --sasa_dir):
     MassSpec will check for cached SASA.npy and collect if missing.
     python scripts/run_compare_sim2exp.py
         --sasa_dir    /path/to/OP_AA/SASA
         --xp_dir      /path/to/OP_AA/XP
         --n_traj      1000
         --sasa_xp_frames_per_traj    335
         --msm_data_file ...  (remaining args as below)

  2. Pre-built mode (provide --sasa_data_file):
     Use pre-built SASA.npy directly (skip internal collection).
     python scripts/run_compare_sim2exp.py
         --sasa_data_file /path/to/SASA.npy
         --dist_data_file /path/to/Jwalk.npy (optional)
         --msm_data_file ...

Full example (mode 1):
python scripts/run_compare_sim2exp.py
--sasa_dir        /path/to/OP_AA/SASA
--xp_dir          /path/to/OP_AA/XP
--n_traj          1000
--sasa_xp_frames_per_traj 335
--msm_data_file   /path/to/msm_data.csv
--meta_dist_file  /path/to/meta_dist.npy
--LiPMS_exp_file  /path/to/LiPMS_exp.xlsx
--XLMS_exp_file   /path/to/XLMS_exp.xlsx
--cluster_data_file /path/to/cluster_data.npz
--OPpath          /path/to/OP_AA/
--AAdcd_dir       /path/to/aa_trajectories/
--native_AA_pdb   /path/to/native.pdb
--state_idx_list  4 6 8
--prot_len        387
--n_analysis_frames 335
--rm_traj_list    65 75 155
--native_state_idx 9
--outdir          /path/to/outdir/
--ID              1ZMR
--num_perm        1000
--n_boot          100
--lag_frame       20
--nproc           10
"""

def main(argv=None):

    import multiprocessing as mp
    import argparse
    import time

    start_time = time.time()
    
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--config", type=str, required=False, default=argparse.SUPPRESS,
                        help="Optional path to JSON/YAML config file. CLI flags override config values.")
    parser.add_argument("--msm_data_file", type=str, required=False, default=argparse.SUPPRESS, help="Path to MSM mapping file")
    parser.add_argument("--meta_dist_file", type=str, required=False, default=argparse.SUPPRESS, help="Path to meta-distance file")
    parser.add_argument("--LiPMS_exp_file", type=str, required=False, default=argparse.SUPPRESS, help="Path to LiP-MS experimental data file")
    parser.add_argument("--XLMS_exp_file", type=str, required=False, default=argparse.SUPPRESS, help="Path to XL-MS experimental data file")
    parser.add_argument("--cluster_data_file", type=str, required=False, default=argparse.SUPPRESS, help="Path to clustering data file")
    parser.add_argument("--OPpath", type=str, required=False, default=argparse.SUPPRESS, help="Path to order parameters directory")
    parser.add_argument("--AAdcd_dir", type=str, required=False, default=argparse.SUPPRESS, help="Path to all-atom DCD files directory")
    parser.add_argument("--native_AA_pdb", type=str, required=False, default=argparse.SUPPRESS, help="Path to native all-atom PDB file")
    parser.add_argument("--state_idx_list", type=int, nargs='+', required=False, default=argparse.SUPPRESS, help="List of state indices to analyze")
    parser.add_argument("--prot_len", type=int, required=False, default=argparse.SUPPRESS, help="Length of the protein")
    parser.add_argument("--n_analysis_frames", type=int, required=False, default=argparse.SUPPRESS, help="Number of trailing frames per trajectory actually analyzed (must match the MSM/G/Q/SASA/XP last-frame window)")
    parser.add_argument("--rm_traj_list", type=int, nargs='+', required=False, default=argparse.SUPPRESS, help="List of trajectory indices to remove")
    parser.add_argument("--native_state_idx", type=int, required=False, default=argparse.SUPPRESS, help="Index of the native state")
    parser.add_argument("--outdir", type=str, required=False, default=argparse.SUPPRESS, help="Output directory for results")
    parser.add_argument("--ID", type=str, required=False, default=argparse.SUPPRESS, help="An ID for the analysis")
    parser.add_argument("--verbose", action='store_true', default=argparse.SUPPRESS, help="Enable verbose output")
    parser.add_argument("--num_perm", type=int, required=False, default=argparse.SUPPRESS, help="Number of permutations for statistical tests")
    parser.add_argument("--n_boot", type=int, required=False, default=argparse.SUPPRESS, help="Number of bootstrap samples")
    parser.add_argument("--lag_frame", type=int, required=False, default=argparse.SUPPRESS, help="Lag time in frames")
    parser.add_argument("--nproc", type=int, required=False, default=argparse.SUPPRESS, help="Number of processes for parallel computation")
    # --- CollectOP arguments (optional: collect from per-traj files) ---
    parser.add_argument("--sasa_dir", type=str, default=argparse.SUPPRESS,
                        help="Directory of per-traj {ID}_Traj{N}.SASA files. "
                             "If provided together with --xp_dir, CollectOP is run "
                             "before MassSpec and --sasa_data_file / --dist_data_file "
                             "are set automatically.")
    parser.add_argument("--xp_dir", type=str, default=argparse.SUPPRESS,
                        help="Directory of per-traj {ID}_Traj{N}.XP files (used with --sasa_dir).")
    parser.add_argument("--n_traj", type=int, default=argparse.SUPPRESS,
                        help="Total number of trajectories for CollectOP (required with --sasa_dir).")
    parser.add_argument("--sasa_xp_frames_per_traj", type=int, default=argparse.SUPPRESS,
                        help="Frames per trajectory stored in each SASA/XP file (required with --sasa_dir).")
    parser.add_argument("--collect_jwalk_npy", action='store_true', default=argparse.SUPPRESS,
                        help="Also build Jwalk.npy with CollectOP (legacy path). "
                             "By default, XL-MS scoring streams directly from XP files to reduce memory use.")
    # --- Direct array paths (used when skipping collection) ---
    parser.add_argument("--sasa_data_file", type=str, default=argparse.SUPPRESS,
                        help="Path to pre-built SASA.npy. Required if --sasa_dir is not provided.")
    parser.add_argument("--dist_data_file", type=str, default=argparse.SUPPRESS,
                        help="Path to pre-built Jwalk.npy. Required if --xp_dir is not provided.")
    parser.add_argument("--restart", action='store_true', default=argparse.SUPPRESS,
                        help="Resume an interrupted viz_rep_struct run. Keeps the existing directory "
                             "and skips any group that already has a .done sentinel file.")
    args = parse_args_with_config(
        parser,
        argv,
        defaults={
            "sasa_dir": None,
            "xp_dir": None,
            "n_traj": None,
            "sasa_xp_frames_per_traj": None,
            "collect_jwalk_npy": False,
            "sasa_data_file": None,
            "dist_data_file": None,
            "verbose": False,
            "restart": False,
        },
        aliases={"id": "ID"},
    )

    required_fields = [
        "msm_data_file", "meta_dist_file", "LiPMS_exp_file", "XLMS_exp_file",
        "cluster_data_file", "OPpath", "AAdcd_dir", "native_AA_pdb",
        "state_idx_list", "prot_len", "n_analysis_frames", "rm_traj_list",
        "native_state_idx", "outdir", "ID",
        "num_perm", "n_boot", "lag_frame", "nproc",
    ]
    for field in required_fields:
        if not hasattr(args, field) or getattr(args, field) is None:
            parser.error(f"Missing required argument: --{field} (or provide it in --config)")

    if isinstance(args.state_idx_list, int):
        args.state_idx_list = [args.state_idx_list]
    if isinstance(args.rm_traj_list, int):
        args.rm_traj_list = [args.rm_traj_list]

    print(args)
    msm_data_file = args.msm_data_file
    meta_dist_file = args.meta_dist_file
    LiPMS_exp_file = args.LiPMS_exp_file
    sasa_data_file = args.sasa_data_file
    XLMS_exp_file = args.XLMS_exp_file
    dist_data_file = args.dist_data_file
    cluster_data_file = args.cluster_data_file
    OPpath = args.OPpath
    AAdcd_dir = args.AAdcd_dir
    native_AA_pdb = args.native_AA_pdb
    state_idx_list = args.state_idx_list
    prot_len = args.prot_len
    n_analysis_frames = args.n_analysis_frames
    rm_traj_list = args.rm_traj_list
    native_state_idx = args.native_state_idx
    outdir = args.outdir
    ID = args.ID
    verbose = args.verbose
    num_perm = args.num_perm
    n_boot = args.n_boot
    lag_frame = args.lag_frame
    nproc = args.nproc

    # ── validate input mode ────────────────────────────────────────────────
    collect_mode = args.sasa_dir is not None
    direct_mode  = args.sasa_data_file is not None

    if not collect_mode and not direct_mode:
        parser.error(
            "Provide either --sasa_dir (for automatic collection) or "
            "--sasa_data_file (for pre-built array)."
        )
    if collect_mode and args.n_traj is None:
        parser.error("--n_traj is required when using --sasa_dir (or provide it in --config)")
    if collect_mode and args.sasa_xp_frames_per_traj is None:
        parser.error("--sasa_xp_frames_per_traj is required when using --sasa_dir (or provide it in --config)")

    # ── Step 1: Initialize MassSpec (it handles collection internally if needed)
    MS = MassSpec(msm_data_file=msm_data_file,
                    meta_dist_file=meta_dist_file,
                    LiPMS_exp_file=LiPMS_exp_file,
                    XLMS_exp_file=XLMS_exp_file,
                    cluster_data_file=cluster_data_file,
                    OPpath=OPpath,
                    AAdcd_dir=AAdcd_dir,
                    native_AA_pdb=native_AA_pdb,
                    sasa_data_file=sasa_data_file,
                    dist_data_file=dist_data_file,
                    sasa_dir=args.sasa_dir,
                    n_traj=args.n_traj,
                    sasa_xp_frames_per_traj=args.sasa_xp_frames_per_traj,
                    collect_jwalk_npy=args.collect_jwalk_npy,
                    xp_dir=args.xp_dir,
                    state_idx_list=state_idx_list,
                    prot_len=prot_len,
                    n_analysis_frames=n_analysis_frames,
                    rm_traj_list=rm_traj_list,
                    native_state_idx=native_state_idx,
                    outdir=outdir,
                    ID=ID,
                    verbose=verbose,
                    num_perm=num_perm,
                    n_boot=n_boot,
                    lag_frame=lag_frame,
                    nproc=nproc)

    # run the consistency test
    consist_data_file, consist_result_file = MS.LiP_XL_MS_ConsistencyTest()
    print(f'consist_data_file: {consist_data_file}')
    print(f'consist_result_file: {consist_result_file}')

    # select the representative structures from the consistency test
    MS.select_rep_structs(consist_data_file, consist_result_file, total_traj_num_frames=335, n_analysis_frames=n_analysis_frames, restart=args.restart)
    
    print(f'NORMAL TERMINATION - {time.time() - start_time} seconds')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
