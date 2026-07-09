#!/usr/bin/env python3
from EntDetect.statistics import MonteCarlo
from EntDetect._logging import setup_logger

if __package__ in {None, ""}:
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from scripts._cli_config import parse_args_with_config
else:
    from ._cli_config import parse_args_with_config

"""
Run Workflow 4 Monte Carlo subpopulation selection.

Provide parameters directly as CLI flags, through `--config` JSON/YAML, or both.
When both are provided, CLI flags override config values.

This script wraps EntDetect.statistics.MonteCarlo and optimizes population
partitions using a logistic-regression objective and penalty terms.

Example
-------
python scripts/run_montecarlo.py \
    --dataframe_files /path/to/residue_dataframes_workflow4.csv \
    --outdir /path/to/workflow4/monte_carlo/ \
    --gene_list /path/to/gene_list.txt \
    --ID Ecoli_population_mc \
    --steps 100000 \
    --n_groups 4 \
    --C1 1.0 \
    --C2 2.5 \
    --beta 0.05

Expected input schema
---------------------
- Either a single combined design matrix file OR one residue table per protein
- Files are pipe-delimited ("|")
- Required columns: gene, mapped_resid, uniprot_length, AA, region, cut_C_Rall
"""


def _normalize_list_arg(value):
    if isinstance(value, str):
        return [item.strip() for item in value.split(",") if item.strip()]
    return value

def main(argv=None):

    import os
    import argparse
    import time
    import logging

    start_time = time.time()

    parser = argparse.ArgumentParser(
        description="Run Workflow 4 Monte Carlo subpopulation selection."
    )
    parser.add_argument("--config", type=str, required=False, default=argparse.SUPPRESS,
                        help="Optional path to JSON/YAML config file. CLI flags override config values.")

    # --- required IO ---
    parser.add_argument("--dataframe_files", type=str, required=False, default=argparse.SUPPRESS,
                        help="Input design matrix path: either a directory of per-protein files or a single combined CSV")
    parser.add_argument("--outdir", type=str, required=False, default=argparse.SUPPRESS,
                        help="Output directory for Monte Carlo results")
    parser.add_argument("--gene_list", type=str, required=False, default=argparse.SUPPRESS,
                        help="Path to gene list file (one ID per line)")
    parser.add_argument("--ID", type=str, required=False, default=argparse.SUPPRESS,
                        help="Identifier used for output naming")

    # --- model options ---
    parser.add_argument("--reg_formula", type=str, default=argparse.SUPPRESS,
                        help="Regression formula used by state scoring")
    parser.add_argument("--sep", type=str, default=argparse.SUPPRESS,
                        help="Column separator used when loading residue data (default: '|')")
    parser.add_argument("--reg_var", nargs='+', default=argparse.SUPPRESS,
                        help="Predictor columns used for regression model encoding (default: AA region)")
    parser.add_argument("--response_var", type=str, default=argparse.SUPPRESS,
                        help="Response variable in regression")
    parser.add_argument("--var2binarize", nargs='+', default=argparse.SUPPRESS,
                        help="Columns binarized during preprocessing (default: cut_C_Rall region)")
    parser.add_argument("--mask_column", type=str, default=argparse.SUPPRESS,
                        help="Residue index column used as mask/filter key (default: mapped_resid)")
    parser.add_argument("--ID_column", type=str, default=argparse.SUPPRESS,
                        help="Protein identifier column used during Monte Carlo loading and grouping (default: gene)")
    parser.add_argument("--Length_column", type=str, default=argparse.SUPPRESS,
                        help="Protein length column used during Monte Carlo loading (default: uniprot_length)")
    parser.add_argument("--test_var", type=str, default=argparse.SUPPRESS,
                        help="Primary test variable")
    parser.add_argument("--random", action='store_true', default=argparse.SUPPRESS,
                        help="Use random sampling mode")
    parser.add_argument("--n_groups", type=int, default=argparse.SUPPRESS,
                        help="Number of groups (default: 4)")
    parser.add_argument("--steps", type=int, default=argparse.SUPPRESS,
                        help="Number of Monte Carlo steps (default: 100000)")
    parser.add_argument("--C1", type=float, default=argparse.SUPPRESS,
                        help="Monte Carlo objective weight C1")
    parser.add_argument("--C2", type=float, default=argparse.SUPPRESS,
                        help="Monte Carlo objective weight C2")
    parser.add_argument("--beta", type=float, default=argparse.SUPPRESS,
                        help="Inverse temperature/annealing parameter beta")
    parser.add_argument("--linearT", action='store_true', default=argparse.SUPPRESS,
                        help="Use linear temperature schedule")

    # --- logging ---
    parser.add_argument("--log_level", default=argparse.SUPPRESS, choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        help="Logging verbosity (default: INFO)")
    parser.add_argument("--logdir", type=str, default=argparse.SUPPRESS,
                        help="Directory for log files (default: same as --outdir)")

    args = parse_args_with_config(
        parser,
        argv,
        defaults={
            "reg_formula": 'cut_C_Rall ~ region + AA',
            "sep": '|',
            "reg_var": ['AA', 'region'],
            "response_var": 'cut_C_Rall',
            "var2binarize": ['cut_C_Rall', 'region'],
            "mask_column": 'mapped_resid',
            "ID_column": 'gene',
            "Length_column": 'uniprot_length',
            "test_var": 'region',
            "random": False,
            "n_groups": 4,
            "steps": 100000,
            "C1": 1.0,
            "C2": 2.5,
            "beta": 0.05,
            "linearT": False,
            "log_level": "INFO",
            "logdir": None,
        },
    )

    args.reg_var = _normalize_list_arg(args.reg_var)
    args.var2binarize = _normalize_list_arg(args.var2binarize)

    for field in ["dataframe_files", "outdir", "gene_list", "ID"]:
        if not hasattr(args, field) or getattr(args, field) is None:
            parser.error(f"Missing required argument: --{field} (or provide it in --config)")

    dataframe_files = args.dataframe_files
    outdir = args.outdir
    gene_list = args.gene_list
    ID = args.ID
    reg_formula = args.reg_formula
    response_var = args.response_var
    test_var = args.test_var
    random = args.random
    n_groups = args.n_groups
    steps = args.steps
    C1 = args.C1
    C2 = args.C2
    beta = args.beta
    linearT = args.linearT

    log_level = getattr(logging, args.log_level.upper(), logging.INFO)
    logdir = args.logdir if args.logdir is not None else outdir
    os.makedirs(logdir, exist_ok=True)

    logger = setup_logger('run_montecarlo', outdir=logdir, ID=ID, log_level=log_level)
    setup_logger('MonteCarlo', outdir=logdir, ID=ID, log_level=log_level)
    logger.info(f'args: {args}')

    # --- input validation ---
    if not os.path.exists(dataframe_files):
        parser.error(f"--dataframe_files does not exist: {dataframe_files}")
    if not os.path.isfile(gene_list):
        parser.error(f"--gene_list does not exist or is not a file: {gene_list}")
    if n_groups < 1:
        parser.error("--n_groups must be >= 1")
    if steps < 1:
        parser.error("--steps must be >= 1")
    os.makedirs(outdir, exist_ok=True)

    # --- step 1: initialize Monte Carlo object ---
    MC = MonteCarlo(
        dataframe_files=dataframe_files,
        outdir=outdir,
        gene_list=gene_list,
        ID=ID,
        reg_formula=reg_formula,
        response_var=response_var,
        test_var=test_var,
        random=random,
        n_groups=n_groups,
        steps=steps,
        C1=C1,
        C2=C2,
        beta=beta,
        linearT=linearT,
        log_level=log_level,
        logdir=logdir,
    )
    logger.info(f'MonteCarlo: {MC}')

    # --- step 2: load residue-level data ---
    MC.load_data(
        sep=args.sep,
        reg_var=args.reg_var,
        var2binarize=args.var2binarize,
        mask_column=args.mask_column,
        ID_column=args.ID_column,
        Length_column=args.Length_column,
    )

    # --- step 3: run simulation ---
    MC.run(encoded_df=MC.data, ID_column=args.ID_column)
    
    logger.info(f'NORMAL TERMINATION - {time.time() - start_time:.1f} seconds')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
