#!/usr/bin/env python3
from EntDetect.statistics import MonteCarlo
from EntDetect._logging import setup_logger

"""
Run Workflow 4 Monte Carlo subpopulation selection.

This script wraps EntDetect.statistics.MonteCarlo and optimizes population
partitions using a logistic-regression objective and penalty terms.

Example
-------
python scripts/run_montecarlo.py \
    --dataframe_files /path/to/residue_dataframes_workflow4.csv \
    --outpath /path/to/workflow4/monte_carlo/ \
    --gene_list /path/to/gene_list.txt \
    --tag Ecoli_population_mc \
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

def main(argv=None):

    import os
    import argparse
    import time
    import logging

    start_time = time.time()

    parser = argparse.ArgumentParser(
        description="Run Workflow 4 Monte Carlo subpopulation selection."
    )

    # --- required IO ---
    parser.add_argument("--dataframe_files", type=str, required=True,
                        help="Input design matrix path: either a directory of per-protein files or a single combined CSV")
    parser.add_argument("--outpath", type=str, required=True,
                        help="Output directory for Monte Carlo results")
    parser.add_argument("--gene_list", type=str, required=True,
                        help="Path to gene list file (one ID per line)")
    parser.add_argument("--tag", type=str, required=True,
                        help="Identifier tag for output naming")

    # --- model options ---
    parser.add_argument("--reg_formula", type=str, default='cut_C_Rall ~ region + AA',
                        help="Regression formula used by state scoring")
    parser.add_argument("--response_var", type=str, default='cut_C_Rall',
                        help="Response variable in regression")
    parser.add_argument("--test_var", type=str, default='region',
                        help="Primary test variable")
    parser.add_argument("--random", action='store_true',
                        help="Use random sampling mode")
    parser.add_argument("--n_groups", type=int, default=4,
                        help="Number of groups (default: 4)")
    parser.add_argument("--steps", type=int, default=100000,
                        help="Number of Monte Carlo steps (default: 100000)")
    parser.add_argument("--C1", type=float, default=1.0,
                        help="Monte Carlo objective weight C1")
    parser.add_argument("--C2", type=float, default=2.5,
                        help="Monte Carlo objective weight C2")
    parser.add_argument("--beta", type=float, default=0.05,
                        help="Inverse temperature/annealing parameter beta")
    parser.add_argument("--linearT", action='store_true',
                        help="Use linear temperature schedule")

    # --- logging ---
    parser.add_argument("--log_level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        help="Logging verbosity (default: INFO)")
    parser.add_argument("--logdir", type=str, default=None,
                        help="Directory for log files (default: same as --outpath)")

    args = parser.parse_args(argv)

    dataframe_files = args.dataframe_files
    outdir = args.outpath
    gene_list = args.gene_list
    tag = args.tag
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

    logger = setup_logger('run_montecarlo', outdir=logdir, ID=tag, log_level=log_level)
    setup_logger('MonteCarlo', outdir=logdir, ID=tag, log_level=log_level)
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
        ID=tag,
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
        sep='|',
        reg_var=['AA', 'region'],
        response_var='cut_C_Rall',
        var2binarize=['cut_C_Rall', 'region'],
        mask_column='mapped_resid',
        ID_column='gene',
        Length_column='uniprot_length',
    )

    # --- step 3: run simulation ---
    MC.run(encoded_df=MC.data, ID_column='gene')
    
    logger.info(f'NORMAL TERMINATION - {time.time() - start_time:.1f} seconds')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
