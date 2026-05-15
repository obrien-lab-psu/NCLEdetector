from EntDetect.clustering import ClusterNonNativeEntanglements
from EntDetect._logging import setup_logger

"""
Cluster non-native entanglement changes across an ensemble of simulation trajectories.

Reads per-trajectory entanglement pkl files produced by run_OP_on_simulation_traj.py
(located in the Combined_GE/ subdirectory of the G/ output folder), groups them into
non-redundant entanglement-change clusters, and writes representative structures and
per-frame cluster assignments to --outdir.

Examples
--------
Basic run:
  python scripts/run_nonnative_entanglement_clustering.py \\
    --outdir              $DATASTORE/outputs/workflow2/nonnative_clustering \\
    --trajnum2pklfile_path $DATASTORE/user_input/metadata/trajnum2file.txt \\
    --traj_dir_prefix     $DATASTORE/user_input/cg_trajectories

Flags
-----
  --outdir               Output directory for clustering results
  --trajnum2pklfile_path CSV file (source of truth) with columns: trajnum, pklfile
                         Users control exactly which pkl files to analyze via this file
  --traj_dir_prefix      Path prefix to the directory containing trajectory DCD files
  --start_frame          First frame index to include, 0-based (default: 0)
  --end_frame            Last frame index to include, 0-based (default: all frames)
  --nproc                Number of parallel worker threads (default: 1)
                         Parallelises both pkl loading (per trajectory) and
                         entanglement-keyword clustering (per unique keyword).
                         Use the number of available CPU cores for best speed.
  --log_level            Logging verbosity: DEBUG, INFO, WARNING, ERROR (default: INFO)
  --logdir               Directory for log file (default: same as --outdir)
"""


def main(argv=None):

    ###---------------------------------------------------------------------------------------------------------
    import sys, os
    import argparse
    import time
    import logging
    start_time = time.time()
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(
        description="Cluster non-native entanglement changes across simulation trajectories.")

    # --- identity / IO ---
    parser.add_argument("--outdir",               type=str, required=True,  help="Output directory for clustering results")
    parser.add_argument("--trajnum2pklfile_path",  type=str, required=True,  help="CSV file mapping trajectory numbers to pkl file paths (source of truth for which files to analyze)")
    parser.add_argument("--traj_dir_prefix",       type=str, required=True,  help="Path prefix to the directory containing trajectory DCD files")

    # --- frame selection ---
    parser.add_argument("--start_frame", type=int, default=0,           help="First frame index to include, 0-based (default: 0)")
    parser.add_argument("--end_frame",   type=int, default=9999999,     help="Last frame index to include, 0-based (default: all frames)")

    # --- parallelism ---
    parser.add_argument("--nproc",       type=int, default=1,           help="Number of parallel worker threads (default: 1)")

    # --- logging ---
    parser.add_argument("--log_level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging verbosity (default: INFO)")
    parser.add_argument("--logdir",    type=str, default=None, help="Directory for log file (default: same as --outdir)")

    args = parser.parse_args(argv)

    outdir = args.outdir
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    log_level = getattr(logging, args.log_level.upper(), logging.INFO)
    logdir    = args.logdir if args.logdir is not None else outdir

    logger = setup_logger('run_nonnative_clustering', outdir=logdir, ID='ClusterNonNativeEntanglements', log_level=log_level)
    setup_logger('ClusterNonNativeEntanglements', outdir=logdir, ID='ClusterNonNativeEntanglements', log_level=log_level)
    logger.info(f'args: {args}')
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    # --- input validation ---
    if not os.path.isfile(args.trajnum2pklfile_path):
        parser.error(f"--trajnum2pklfile_path does not exist: {args.trajnum2pklfile_path}")

    if not os.path.isdir(args.traj_dir_prefix):
        parser.error(f"--traj_dir_prefix does not exist or is not a directory: {args.traj_dir_prefix}")

    os.makedirs(outdir, exist_ok=True)
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    clustering_NNents = ClusterNonNativeEntanglements(
        trajnum2pklfile_path=args.trajnum2pklfile_path,
        traj_dir_prefix=args.traj_dir_prefix,
        outdir=outdir,
        log_level=log_level,
        logdir=logdir,
        nproc=args.nproc,
    )
    logger.info(f'ClusterNonNativeEntanglements: {clustering_NNents}')
    clustering_NNents.cluster(start_frame=args.start_frame, end_frame=args.end_frame)
    ###---------------------------------------------------------------------------------------------------------

    logger.info(f'NORMAL TERMINATION - {time.time() - start_time:.1f} seconds')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
