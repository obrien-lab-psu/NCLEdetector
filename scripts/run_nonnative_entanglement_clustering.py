from NCLEdetector.clustering import ClusterNonNativeEntanglements
from NCLEdetector._logging import setup_logger

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
    import json
    import time
    import logging
    start_time = time.time()
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(
        description="Cluster non-native entanglement changes across simulation trajectories.")

    parser.add_argument("--config", type=str, required=False, default=argparse.SUPPRESS,
                        help="Optional path to JSON or YAML config file. CLI flags override config values.")

    # --- identity / IO ---
    parser.add_argument("--outdir",               type=str, default=argparse.SUPPRESS, help="Output directory for clustering results")
    parser.add_argument("--trajnum2pklfile_path",  type=str, default=argparse.SUPPRESS, help="CSV file mapping trajectory numbers to pkl file paths (source of truth for which files to analyze)")
    parser.add_argument("--traj_dir_prefix",       type=str, default=argparse.SUPPRESS, help="Path prefix to the directory containing trajectory DCD files")

    # --- frame selection ---
    parser.add_argument("--start_frame", type=int, default=argparse.SUPPRESS, help="First frame index to include, 0-based (default: 0)")
    parser.add_argument("--end_frame",   type=int, default=argparse.SUPPRESS, help="Last frame index to include, 0-based (default: all frames)")

    # --- parallelism ---
    parser.add_argument("--nproc",       type=int, default=argparse.SUPPRESS, help="Number of parallel worker threads (default: 1)")

    # --- logging ---
    parser.add_argument("--log_level", default=argparse.SUPPRESS, choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging verbosity (default: INFO)")
    parser.add_argument("--logdir",    type=str, default=argparse.SUPPRESS, help="Directory for log file (default: same as --outdir)")

    def _load_config_file(cfg_path):
        if not os.path.isfile(cfg_path):
            parser.error(f"--config file does not exist: {cfg_path}")
        ext = os.path.splitext(cfg_path)[1].lower()
        try:
            with open(cfg_path, "r", encoding="utf-8") as fh:
                if ext == ".json":
                    cfg = json.load(fh)
                elif ext in {".yml", ".yaml"}:
                    try:
                        import yaml
                    except ImportError:
                        parser.error("YAML config requested but PyYAML is not installed. Use JSON or install PyYAML.")
                    cfg = yaml.safe_load(fh)
                else:
                    parser.error(f"Unsupported config extension '{ext}'. Use .json, .yml, or .yaml.")
        except Exception as exc:
            parser.error(f"Failed to load config file {cfg_path}: {exc}")
        if cfg is None:
            cfg = {}
        if not isinstance(cfg, dict):
            parser.error("Config file must define a top-level object/dictionary of key-value pairs.")
        return cfg

    cli_args = vars(parser.parse_args(argv))
    config_path = cli_args.pop("config", None)
    config_args = _load_config_file(config_path) if config_path else {}

    merged = {
        "outdir": None,
        "trajnum2pklfile_path": None,
        "traj_dir_prefix": None,
        "start_frame": 0,
        "end_frame": 9999999,
        "nproc": 1,
        "log_level": "INFO",
        "logdir": None,
    }
    merged.update(config_args)
    merged.update(cli_args)

    if merged["outdir"] is None:
        parser.error("Missing required argument: --outdir (or provide 'outdir' in --config)")
    if merged["trajnum2pklfile_path"] is None:
        parser.error("Missing required argument: --trajnum2pklfile_path (or provide it in --config)")
    if merged["traj_dir_prefix"] is None:
        parser.error("Missing required argument: --traj_dir_prefix (or provide it in --config)")

    args = argparse.Namespace(**merged)

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
