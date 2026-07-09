from EntDetect.statistics import MSMStats
from EntDetect._logging import setup_logger

"""
Compute metastable-state probability evolution (MSMStats) from an MSM
trajectory-type-annotated mapping CSV produced by run_MSM.py.

The input --msm_data_file must be a CSV with the columns produced by run_MSM.py
(traj, frame, microstate, metastablestate, Q, G) plus a user-added
trajectory-type column (--traj_type_col) that labels each trajectory as belonging
to one of the types in --traj_type_list (e.g. 'A' for folded, 'B' for unfolded).

This classification is typically added by the user based on a Q threshold, e.g.:
  df['traj_type_A80%Native'] = df.groupby('traj')['Q'].transform('max').ge(0.80).map({True:'A', False:'B'})

Examples
--------
Basic run — compute state probabilities for two trajectory types:
  python scripts/run_MSMStats.py \\
    --msm_data_file $DATASTORE/outputs/workflow2/MSM/1ZMR_prod_MSMmapping_QG_native.csv \\
    --traj_type_col traj_type_QG_native \\
    --traj_type_list A B \\
    --outdir        $DATASTORE/outputs/workflow2/MSM_StateProbabilityStats_QG_native

Excluding mirror-image trajectories:
  python scripts/run_MSMStats.py \\
    --msm_data_file $DATASTORE/outputs/workflow2/MSM/1ZMR_prod_MSMmapping_QG_native.csv \\
    --traj_type_col traj_type_QG_native \\
    --traj_type_list A B \\
    --outdir        $DATASTORE/outputs/workflow2/MSM_StateProbabilityStats_QG_native \\
    --rm_traj_list  65 75 155 162

Flags
-----
  --msm_data_file   CSV produced by run_MSM.py, annotated with a trajectory-type column
  --traj_type_col   Column name in msm_data_file that contains trajectory-type labels
  --traj_type_list  Space-separated list of trajectory-type labels to compare (default: A B)
  --outdir          Output directory for state probability results
  --rm_traj_list    Trajectory numbers to exclude (e.g. confirmed mirror conformations)
  --log_level       Logging verbosity: DEBUG, INFO, WARNING, ERROR (default: INFO)
  --logdir          Directory for log file (default: same as --outdir)
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
        description="Compute metastable-state probability evolution from MSM output.")

    parser.add_argument("--config", type=str, required=False, default=argparse.SUPPRESS,
                        help="Optional path to JSON or YAML config file. CLI flags override config values.")

    # --- IO ---
    parser.add_argument("--msm_data_file", type=str, default=argparse.SUPPRESS, help="CSV produced by run_MSM.py, annotated with a trajectory-type column")
    parser.add_argument("--outdir",        type=str, default=argparse.SUPPRESS, help="Output directory for state probability results")

    # --- trajectory classification ---
    parser.add_argument("--traj_type_col",  type=str,           default=argparse.SUPPRESS, help="Column name in msm_data_file containing trajectory-type labels")
    parser.add_argument("--traj_type_list", type=str, nargs='+', default=argparse.SUPPRESS, help="Trajectory-type labels to compare (default: A B)")

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
        "msm_data_file": None,
        "outdir": None,
        "traj_type_col": None,
        "traj_type_list": ["A", "B"],
        "rm_traj_list": [],
        "log_level": "INFO",
        "logdir": None,
    }
    merged.update(config_args)
    merged.update(cli_args)

    if merged["msm_data_file"] is None:
        parser.error("Missing required argument: --msm_data_file (or provide it in --config)")
    if merged["outdir"] is None:
        parser.error("Missing required argument: --outdir (or provide it in --config)")
    if merged["traj_type_col"] is None:
        parser.error("Missing required argument: --traj_type_col (or provide it in --config)")

    args = argparse.Namespace(**merged)

    outdir = args.outdir
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    log_level = getattr(logging, args.log_level.upper(), logging.INFO)
    logdir    = args.logdir if args.logdir is not None else outdir

    os.makedirs(outdir, exist_ok=True)

    logger = setup_logger('run_MSMStats', outdir=logdir, ID='MSMStats', log_level=log_level)
    setup_logger('MSMStats', outdir=logdir, ID='MSMStats', log_level=log_level)
    logger.info(f'args: {args}')
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    # --- input validation ---
    if not os.path.isfile(args.msm_data_file):
        parser.error(f"--msm_data_file does not exist: {args.msm_data_file}")
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    MS = MSMStats(
        outdir=outdir,
        rm_traj_list=args.rm_traj_list,
        log_level=log_level,
        logdir=logdir,
    )
    logger.info(f'MSMStats: {MS}')
    df = MS.StateProbabilityStats(
        msm_data_file=args.msm_data_file,
        traj_type_col=args.traj_type_col,
        traj_type_list=args.traj_type_list,
    )
    MS.Plot_StateProbabilityStats(
        df=df,
        traj_type_col=args.traj_type_col,
        traj_type_list=args.traj_type_list,
    )
    ###---------------------------------------------------------------------------------------------------------

    logger.info(f'NORMAL TERMINATION - {time.time() - start_time:.1f} seconds')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
