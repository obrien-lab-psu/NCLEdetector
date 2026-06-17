from EntDetect.gaussian_entanglement import GaussianEntanglement
from EntDetect.clustering import ClusterNativeEntanglements, MSMNonNativeEntanglementClustering
from EntDetect.order_params import CalculateOP, CollectOP
from EntDetect.compare_sim2exp import MassSpec

"""
Collect per-trajectory SASA/XP outputs (optional) and run the LiP-MS / XL-MS
consistency test.

Two usage modes:

  1. Collect + run  (provide --sasa_dir and --xp_dir):
     python scripts/run_compare_sim2exp.py
         --sasa_dir   /path/to/OP_AA/SASA
         --xp_dir     /path/to/OP_AA/XP
         --n_traj     1000
         --n_frames   335
         --msm_data_file ...  (remaining args as below)

  2. Skip collection (provide --sasa_data_file and --dist_data_file directly):
     python scripts/run_compare_sim2exp.py
         --sasa_data_file /path/to/SASA.npy
         --dist_data_file /path/to/Jwalk.npy
         --msm_data_file ...

Full example (mode 1):
python scripts/run_compare_sim2exp.py
--sasa_dir        /path/to/OP_AA/SASA
--xp_dir          /path/to/OP_AA/XP
--n_traj          1000
--n_frames        335
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
--last_num_frames 335
--rm_traj_list    65 75 155
--native_state_idx 9
--outdir          /path/to/outdir/
--ID              1ZMR
--start           6600
--end             -1
--stride          1
--num_perm        1000
--n_boot          100
--lag_frame       20
--nproc           10
"""

def main(argv=None):

    import multiprocessing as mp
    import sys, os
    import argparse
    import time

    start_time = time.time()
    
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--msm_data_file", type=str, required=True, help="Path to MSM mapping file")
    parser.add_argument("--meta_dist_file", type=str, required=True, help="Path to meta-distance file")
    parser.add_argument("--LiPMS_exp_file", type=str, required=True, help="Path to LiP-MS experimental data file")
    parser.add_argument("--XLMS_exp_file", type=str, required=True, help="Path to XL-MS experimental data file")
    parser.add_argument("--cluster_data_file", type=str, required=True, help="Path to clustering data file")
    parser.add_argument("--OPpath", type=str, required=True, help="Path to order parameters directory")
    parser.add_argument("--AAdcd_dir", type=str, required=True, help="Path to all-atom DCD files directory")
    parser.add_argument("--native_AA_pdb", type=str, required=True, help="Path to native all-atom PDB file")
    parser.add_argument("--state_idx_list", type=int, nargs='+', required=True, help="List of state indices to analyze")
    parser.add_argument("--prot_len", type=int, required=True, help="Length of the protein")
    parser.add_argument("--last_num_frames", type=int, required=True, help="Number of last frames to consider")
    parser.add_argument("--rm_traj_list", type=int, nargs='+', required=True, help="List of trajectory indices to remove")
    parser.add_argument("--native_state_idx", type=int, required=True, help="Index of the native state")
    parser.add_argument("--outdir", type=str, required=True, help="Output directory for results")
    parser.add_argument("--ID", type=str, required=True, help="An ID for the analysis")
    parser.add_argument("--start", type=int, required=True, help="Start frame index")
    parser.add_argument("--end", type=int, required=True, help="End frame index")
    parser.add_argument("--stride", type=int, required=True, help="Stride for frame selection")
    parser.add_argument("--verbose", action='store_true', help="Enable verbose output")
    parser.add_argument("--num_perm", type=int, required=True, help="Number of permutations for statistical tests")
    parser.add_argument("--n_boot", type=int, required=True, help="Number of bootstrap samples")
    parser.add_argument("--lag_frame", type=int, required=True, help="Lag time in frames")
    parser.add_argument("--nproc", type=int, required=True, help="Number of processes for parallel computation")
    # --- CollectOP arguments (optional: collect from per-traj files) ---
    parser.add_argument("--sasa_dir", type=str, default=None,
                        help="Directory of per-traj {ID}_Traj{N}.SASA files. "
                             "If provided together with --xp_dir, CollectOP is run "
                             "before MassSpec and --sasa_data_file / --dist_data_file "
                             "are set automatically.")
    parser.add_argument("--xp_dir", type=str, default=None,
                        help="Directory of per-traj {ID}_Traj{N}.XP files (used with --sasa_dir).")
    parser.add_argument("--n_traj", type=int, default=None,
                        help="Total number of trajectories for CollectOP (required with --sasa_dir).")
    parser.add_argument("--n_frames", type=int, default=None,
                        help="Frames per trajectory stored in each file (required with --sasa_dir).")
    parser.add_argument("--collect_jwalk_npy", action='store_true',
                        help="Also build Jwalk.npy with CollectOP (legacy path). "
                             "By default, XL-MS scoring streams directly from XP files to reduce memory use.")
    # --- Direct array paths (used when skipping collection) ---
    parser.add_argument("--sasa_data_file", type=str, default=None,
                        help="Path to pre-built SASA.npy. Required if --sasa_dir is not provided.")
    parser.add_argument("--dist_data_file", type=str, default=None,
                        help="Path to pre-built Jwalk.npy. Required if --xp_dir is not provided.")
    args = parser.parse_args(argv)
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
    last_num_frames = args.last_num_frames
    rm_traj_list = args.rm_traj_list
    native_state_idx = args.native_state_idx
    outdir = args.outdir
    ID = args.ID
    start = args.start
    end = args.end
    stride = args.stride
    verbose = args.verbose
    num_perm = args.num_perm
    n_boot = args.n_boot
    lag_frame = args.lag_frame
    nproc = args.nproc

    # ── validate input mode ────────────────────────────────────────────────
    collect_mode = args.sasa_dir is not None and args.xp_dir is not None
    direct_mode  = args.sasa_data_file is not None and args.dist_data_file is not None

    if not collect_mode and not direct_mode:
        parser.error(
            "Provide either (--sasa_dir + --xp_dir + --n_traj + --n_frames) "
            "to collect from per-trajectory files, or "
            "(--sasa_data_file + --dist_data_file) to use pre-built arrays."
        )

    # ── Step 1: collect per-trajectory outputs if requested ───────────────
    if collect_mode:
        if args.n_traj is None or args.n_frames is None:
            parser.error("--n_traj and --n_frames are required when using --sasa_dir / --xp_dir")

        os.makedirs(outdir, exist_ok=True)
        collector = CollectOP(
            sasa_dir  = args.sasa_dir,
            xp_dir    = args.xp_dir,
            outdir    = outdir,
            ID        = ID,
            n_traj    = args.n_traj,
            n_frames  = args.n_frames,
            prot_len  = prot_len,
        )
        sasa_data_file = collector.collect_SASA()
        if args.collect_jwalk_npy:
            dist_data_file = collector.collect_Jwalk()
        else:
            dist_data_file = None
        print(f'CollectOP SASA:  {sasa_data_file}')
        if dist_data_file is not None:
            print(f'CollectOP Jwalk: {dist_data_file}')
        else:
            print('CollectOP Jwalk: skipped (streaming XP mode enabled)')

    # ── Step 2: run the consistency test ──────────────────────────────────
    MS = MassSpec(msm_data_file=msm_data_file,
                    meta_dist_file=meta_dist_file,
                    LiPMS_exp_file=LiPMS_exp_file,
                    sasa_data_file=sasa_data_file,
                    XLMS_exp_file=XLMS_exp_file,
                    dist_data_file=dist_data_file,
                    cluster_data_file=cluster_data_file,
                    OPpath=OPpath,
                    AAdcd_dir=AAdcd_dir,
                    native_AA_pdb=native_AA_pdb,
                    xp_dir=args.xp_dir,
                    state_idx_list=state_idx_list,
                    prot_len=prot_len,
                    last_num_frames=last_num_frames,
                    rm_traj_list=rm_traj_list,
                    native_state_idx=native_state_idx,
                    outdir=outdir,
                    ID=ID,
                    start=start,
                    end=end,
                    stride=stride,
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
    MS.select_rep_structs(consist_data_file, consist_result_file, total_traj_num_frames=335, last_num_frames=67)
    
    print(f'NORMAL TERMINATION - {time.time() - start_time} seconds')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
