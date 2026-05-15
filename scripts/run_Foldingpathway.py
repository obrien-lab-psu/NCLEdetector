from EntDetect.statistics import FoldingPathwayStats
from EntDetect._logging import setup_logger
import pandas as pd

"""
Compute folding pathway statistics and Jensen-Shannon divergence from an MSM
trajectory-type-annotated mapping CSV produced by run_MSM.py.

The input --msm_data_file must be a CSV with the columns produced by run_MSM.py
(traj, frame, microstate, metastablestate, Q, G, StateSample) plus a user-added
trajectory-type column (--traj_type_col) that labels each trajectory as belonging
to one of the types in --traj_type_list (e.g. 'A' for folded, 'B' for unfolded).

This classification is typically added by the user based on a Q threshold, e.g.:
  df['traj_type_A80%Native'] = df.groupby('traj')['Q'].transform('max').ge(0.80).map({True:'A', False:'B'})

Examples
--------
Basic run — two trajectory types, no exclusions:
  python scripts/run_Foldingpathway.py \\
    --msm_data_file $DATASTORE/outputs/workflow2/MSM/1ZMR_prod_MSMmapping_A80pctNative.csv \\
    --meta_set_file $DATASTORE/outputs/workflow2/MSM/1ZMR_prod_meta_set.csv \\
    --traj_type_col traj_type_A80pctNative \\
    --traj_type_list A B \\
    --outdir        $DATASTORE/outputs/workflow2/FoldingPathway_A80pctNative

Excluding mirror-image trajectories identified in Step 4:
  python scripts/run_Foldingpathway.py \\
    --msm_data_file $DATASTORE/outputs/workflow2/MSM/1ZMR_prod_MSMmapping_A80pctNative.csv \\
    --meta_set_file $DATASTORE/outputs/workflow2/MSM/1ZMR_prod_meta_set.csv \\
    --traj_type_col traj_type_A80pctNative \\
    --traj_type_list A B \\
    --outdir        $DATASTORE/outputs/workflow2/FoldingPathway_A80pctNative \\
    --rm_traj_list  65 75 155 162

Flags
-----
  --msm_data_file   CSV produced by run_MSM.py, annotated with a trajectory-type column
  --meta_set_file   meta_set CSV produced by run_MSM.py (microstates per metastable state)
  --traj_type_col   Column name in msm_data_file that contains trajectory-type labels
  --traj_type_list  Space-separated list of trajectory-type labels to compare (default: A B)
  --outdir          Output directory for folding pathway and JS-divergence results
  --rm_traj_list    Trajectory numbers to exclude (e.g. confirmed mirror conformations)
  --n_window        Rolling window size for state probability smoothing (default: 200)
  --n_traj          Total number of trajectories in the ensemble (default: 1000)
  --state_type      State level to analyse: metastablestate or microstate (default: metastablestate)
  --log_level       Logging verbosity: DEBUG, INFO, WARNING, ERROR (default: INFO)
  --logdir          Directory for log file (default: same as --outdir)
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
        description="Compute folding pathway statistics and Jensen-Shannon divergence from MSM output.")

    # --- IO ---
    parser.add_argument("--msm_data_file", type=str, required=True,  help="CSV produced by run_MSM.py, annotated with a trajectory-type column")
    parser.add_argument("--meta_set_file", type=str, required=True,  help="meta_set CSV produced by run_MSM.py")
    parser.add_argument("--outdir",        type=str, required=True,  help="Output directory for folding pathway and JS-divergence results")

    # --- trajectory classification ---
    parser.add_argument("--traj_type_col",  type=str,          required=True,        help="Column name in msm_data_file containing trajectory-type labels")
    parser.add_argument("--traj_type_list", type=str, nargs='+', default=['A', 'B'], help="Trajectory-type labels to compare (default: A B)")

    # --- trajectory filtering ---
    parser.add_argument("--rm_traj_list", type=int, nargs='+', default=[], help="Trajectory numbers to exclude (e.g. confirmed mirror conformations)")

    # --- analysis settings ---
    parser.add_argument("--n_window",   type=int, default=200,              help="Rolling window size for state probability smoothing (default: 200)")
    parser.add_argument("--n_traj",     type=int, default=1000,             help="Total number of trajectories in the ensemble (default: 1000)")
    parser.add_argument("--state_type", type=str, default='metastablestate',
                        choices=['metastablestate', 'microstate'],          help="State level to analyse (default: metastablestate)")

    # --- logging ---
    parser.add_argument("--log_level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging verbosity (default: INFO)")
    parser.add_argument("--logdir",    type=str, default=None, help="Directory for log file (default: same as --outdir)")

    args = parser.parse_args(argv)

    outdir = args.outdir
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    log_level = getattr(logging, args.log_level.upper(), logging.INFO)
    logdir    = args.logdir if args.logdir is not None else outdir

    os.makedirs(outdir, exist_ok=True)

    logger = setup_logger('run_Foldingpathway', outdir=logdir, ID='FoldingPathwayStats', log_level=log_level)
    setup_logger('FoldingPathwayStats', outdir=logdir, ID='FoldingPathwayStats', log_level=log_level)
    logger.info(f'args: {args}')
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    # --- input validation ---
    if not os.path.isfile(args.msm_data_file):
        parser.error(f"--msm_data_file does not exist: {args.msm_data_file}")

    if not os.path.isfile(args.meta_set_file):
        parser.error(f"--meta_set_file does not exist: {args.meta_set_file}")
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    # Load MSM data and validate the trajectory-type column
    logger.info(f'Loading MSM data from {args.msm_data_file}')
    msm_data = pd.read_csv(args.msm_data_file)
    logger.info(f'msm_data shape: {msm_data.shape}, columns: {msm_data.columns.tolist()}')

    if args.traj_type_col not in msm_data.columns:
        parser.error(
            f"--traj_type_col '{args.traj_type_col}' not found in {args.msm_data_file}. "
            f"Available columns: {msm_data.columns.tolist()}"
        )

    present_types = set(msm_data[args.traj_type_col].unique())
    missing_types = [t for t in args.traj_type_list if t not in present_types]
    if missing_types:
        parser.error(
            f"--traj_type_list values {missing_types} not found in column '{args.traj_type_col}'. "
            f"Values present: {sorted(present_types)}"
        )
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    FP = FoldingPathwayStats(
        msm_data=msm_data,
        meta_set_file=args.meta_set_file,
        tarj_type_col=args.traj_type_col,
        traj_type_list=args.traj_type_list,
        outdir=outdir,
        n_window=args.n_window,
        n_traj=args.n_traj,
        state_type=args.state_type,
        rm_traj_list=args.rm_traj_list,
        log_level=log_level,
        logdir=logdir,
    )
    logger.info(f'FoldingPathwayStats: {FP}')

    folding_pathways = FP.post_trans()
    logger.info(f'folding_pathways:\n{folding_pathways}')

    JS_divergence = FP.JS_divergence()
    logger.info(f'JS_divergence:\n{JS_divergence}')
    ###---------------------------------------------------------------------------------------------------------

    logger.info(f'NORMAL TERMINATION - {time.time() - start_time:.1f} seconds')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
