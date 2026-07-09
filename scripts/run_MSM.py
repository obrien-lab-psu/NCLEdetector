from EntDetect.clustering import MSMNonNativeEntanglementClustering
from EntDetect._logging import setup_logger

"""
Build a Markov state model (MSM) from pre-computed order-parameter data across
an ensemble of simulation trajectories.

Reads Q and G order-parameter files from --OPpath (which must contain Q/ and G/
subdirectories produced by run_OP_on_simulation_traj.py), groups frames into
microstates via k-means clustering, and then coarse-grains microstates into
metastable macro-states using PCCA+.

Examples
--------
Basic run — 10 metastable states, lag time 20:
  python scripts/run_MSM.py \\
    --outdir  $DATASTORE/outputs/workflow2/MSM \\
    --ID      1ZMR_prod \\
    --OPpath  $DATASTORE/outputs/workflow2/OP_demo/ \\
    --start   0 \\
    --n_large_states 10 \\
    --lagtime 20

Excluding mirror-image trajectories (identified in Step 10):
  python scripts/run_MSM.py \\
    --outdir  $DATASTORE/outputs/workflow2/MSM \\
    --ID      1ZMR_prod \\
    --OPpath  $DATASTORE/outputs/workflow2/OP_demo/ \\
    --start   0 \\
    --n_large_states 10 \\
    --lagtime 20 \\
    --rm_traj_list 65 75 155 162

Flags
-----
  --outdir           Output directory for MSM results
  --OPpath           Directory containing Q/ and G/ subdirectories of per-trajectory OP files
  --ID               Base name for output files
  --start            First frame index to include, 0-based (default: 0)
  --end              Last frame index to include, 0-based (default: all frames)
  --stride           Frame stride for loading OP data (default: 1)
  --n_large_states   Number of metastable macro-states requested from PCCA+ (default: 10)
  --n_small_states   Number of inactive micro-state clusters (default: 1)
  --n_cluster        Number of k-means microstates (default: 400)
  --kmean_stride     Frame stride used during k-means clustering (default: 2)
  --lagtime          MSM lag time in frames (default: 20)
  --dt               MD timestep in ns, used for time-axis labelling (default: 1.5e-5)
  --ITS              Run implied timescale analysis to validate lag time: True/False (default: False)
  --rm_traj_list     Trajectory numbers to exclude (e.g. confirmed mirror conformations)
  --log_level        Logging verbosity: DEBUG, INFO, WARNING, ERROR (default: INFO)
  --logdir           Directory for log file (default: same as --outdir)
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
        description="Build a Markov state model from pre-computed order-parameter trajectories.")

    parser.add_argument("--config", type=str, required=False, default=argparse.SUPPRESS,
                        help="Optional path to JSON or YAML config file. CLI flags override config values.")

    # --- identity / IO ---
    parser.add_argument("--outdir",  type=str, default=argparse.SUPPRESS, help="Output directory for MSM results")
    parser.add_argument("--OPpath",  type=str, default=argparse.SUPPRESS, help="Directory containing Q/ and G/ subdirectories of per-trajectory OP files")
    parser.add_argument("--ID",      type=str, default=argparse.SUPPRESS, help="Base name for output files")

    # --- frame selection ---
    parser.add_argument("--start",   type=int, default=argparse.SUPPRESS, help="First frame index to include, 0-based (default: 0)")
    parser.add_argument("--end",     type=int, default=argparse.SUPPRESS, help="Last frame index to include, 0-based (default: all frames)")
    parser.add_argument("--stride",  type=int, default=argparse.SUPPRESS, help="Frame stride for loading OP data (default: 1)")

    # --- MSM settings ---
    parser.add_argument("--n_large_states", type=int,   default=argparse.SUPPRESS, help="Number of metastable macro-states requested from PCCA+ (default: 10)")
    parser.add_argument("--n_small_states", type=int,   default=argparse.SUPPRESS, help="Number of inactive micro-state clusters (default: 1)")
    parser.add_argument("--n_cluster",      type=int,   default=argparse.SUPPRESS, help="Number of k-means microstates (default: 400)")
    parser.add_argument("--kmean_stride",   type=int,   default=argparse.SUPPRESS, help="Frame stride used during k-means clustering (default: 2)")
    parser.add_argument("--lagtime",        type=int,   default=argparse.SUPPRESS, help="MSM lag time in frames (default: 20)")
    parser.add_argument("--dt",             type=float, default=argparse.SUPPRESS, help="MD timestep in ns (default: 1.5e-5)")
    parser.add_argument("--ITS",            type=str,   default=argparse.SUPPRESS, help="Run implied timescale analysis: True/False (default: False)")

    # --- trajectory filtering ---
    parser.add_argument("--rm_traj_list", type=int, nargs='+', default=argparse.SUPPRESS, help="Trajectory numbers to exclude (e.g. confirmed mirror conformations)")

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
                        parser.error("YAML config requested but PyYAML is not installed.")
                    cfg = yaml.safe_load(fh)
                else:
                    parser.error(f"Unsupported config extension '{ext}'. Use .json, .yml, or .yaml.")
        except Exception as exc:
            parser.error(f"Failed to load config file {cfg_path}: {exc}")
        if cfg is None:
            cfg = {}
        if not isinstance(cfg, dict):
            parser.error("Config file must define a top-level object/dictionary.")
        return cfg

    cli_args = vars(parser.parse_args(argv))
    config_path = cli_args.pop("config", None)
    config_args = _load_config_file(config_path) if config_path else {}

    merged = {
        "outdir": None,
        "OPpath": None,
        "ID": None,
        "start": 0,
        "end": 99999999999,
        "stride": 1,
        "n_large_states": 10,
        "n_small_states": 1,
        "n_cluster": 400,
        "kmean_stride": 2,
        "lagtime": 20,
        "dt": 0.015/1000,
        "ITS": "False",
        "rm_traj_list": [],
        "log_level": "INFO",
        "logdir": None,
    }
    merged.update(config_args)
    merged.update(cli_args)

    if merged["outdir"] is None:
        parser.error("Missing required argument: --outdir (or provide it in --config)")
    if merged["OPpath"] is None:
        parser.error("Missing required argument: --OPpath (or provide it in --config)")
    if merged["ID"] is None:
        parser.error("Missing required argument: --ID (or provide it in --config)")

    args = argparse.Namespace(**merged)

    outdir = args.outdir
    OPpath = args.OPpath
    ID     = args.ID
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    log_level = getattr(logging, args.log_level.upper(), logging.INFO)
    logdir    = args.logdir if args.logdir is not None else outdir

    logger = setup_logger('run_MSM', outdir=logdir, ID=ID, log_level=log_level)
    setup_logger('MSMNonNativeEntanglementClustering', outdir=logdir, ID=ID, log_level=log_level)
    logger.info(f'args: {args}')
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    # --- input validation ---
    if not os.path.isdir(OPpath):
        parser.error(f"--OPpath does not exist or is not a directory: {OPpath}")

    for subdir in ('Q', 'G'):
        expected = os.path.join(OPpath, subdir)
        if not os.path.isdir(expected):
            parser.error(f"Expected subdirectory not found in --OPpath: {expected}")

    os.makedirs(outdir, exist_ok=True)
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    MSM = MSMNonNativeEntanglementClustering(
        outdir=outdir,
        ID=ID,
        OPpath=OPpath,
        start=args.start,
        end=args.end,
        stride=args.stride,
        n_large_states=args.n_large_states,
        n_small_states=args.n_small_states,
        n_cluster=args.n_cluster,
        kmean_stride=args.kmean_stride,
        lagtime=args.lagtime,
        dt=args.dt,
        ITS=args.ITS,
        rm_traj_list=args.rm_traj_list,
        log_level=log_level,
        logdir=logdir,
    )
    logger.info(f'MSMNonNativeEntanglementClustering: {MSM}')
    MSM.run()
    ###---------------------------------------------------------------------------------------------------------

    logger.info(f'NORMAL TERMINATION - {time.time() - start_time:.1f} seconds')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
