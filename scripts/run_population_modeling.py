#!/usr/bin/env python3
from EntDetect.statistics import ProteomeLogisticRegression
from EntDetect._logging import setup_logger

"""
Run Workflow 4 proteome-level logistic regression from residue feature tables.

This script wraps EntDetect.statistics.ProteomeLogisticRegression and writes
the final regression table as a pipe-delimited CSV.

Example
-------
python scripts/run_population_modeling.py \
    --dataframe_files /path/to/residue_dataframes_workflow4.csv \
    --outdir /path/to/workflow4/population_modeling/ \
    --gene_list /path/to/gene_list.txt \
    --tag Ecoli_population \
    --reg_formula "cut_C_Rall ~ AA + region"

Expected input schema
---------------------
- Either a single combined design matrix file OR one residue table per protein
- Files are pipe-delimited ("|")
- Required columns: gene, mapped_resid, AA, region, cut_C_Rall
"""

def main(argv=None):

    import os
    import argparse
    import time
    import logging

    start_time = time.time()

    parser = argparse.ArgumentParser(
        description="Run Workflow 4 proteome-level logistic regression from residue feature tables."
    )

    # --- required IO ---
    parser.add_argument("--dataframe_files", type=str, required=True,
                        help="Input design matrix path: either a directory of per-protein files or a single combined CSV")
    parser.add_argument("--outdir", type=str, required=True,
                        help="Output directory for regression results")
    parser.add_argument("--gene_list", type=str, required=True,
                        help="Path to gene list file (one ID per line)")
    parser.add_argument("--tag", type=str, required=True,
                        help="Identifier tag for output naming")

    # --- model options ---
    parser.add_argument("--reg_formula", type=str, default='cut_C_Rall ~ AA + region',
                        help="Regression formula (default: 'cut_C_Rall ~ AA + region')")

    # --- logging ---
    parser.add_argument("--log_level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        help="Logging verbosity (default: INFO)")
    parser.add_argument("--logdir", type=str, default=None,
                        help="Directory for log files (default: same as --outdir)")

    args = parser.parse_args(argv)

    dataframe_files = args.dataframe_files
    outdir = args.outdir
    gene_list = args.gene_list
    tag = args.tag
    reg_formula = args.reg_formula

    log_level = getattr(logging, args.log_level.upper(), logging.INFO)
    logdir = args.logdir if args.logdir is not None else outdir
    os.makedirs(logdir, exist_ok=True)

    logger = setup_logger('run_population_modeling', outdir=logdir, ID=tag, log_level=log_level)
    setup_logger('ProteomeLogisticRegression', outdir=logdir, ID=tag, log_level=log_level)
    logger.info(f'args: {args}')

    # --- input validation ---
    if not os.path.exists(dataframe_files):
        parser.error(f"--dataframe_files does not exist: {dataframe_files}")
    if not os.path.isfile(gene_list):
        parser.error(f"--gene_list does not exist or is not a file: {gene_list}")
    os.makedirs(outdir, exist_ok=True)

    ## initialize the regression object
    ProtRegression = ProteomeLogisticRegression(
        dataframe_files=dataframe_files,
        outdir=outdir,
        gene_list=gene_list,
        ID=tag,
        reg_formula=reg_formula,
        log_level=log_level,
        logdir=logdir,
    )
    logger.info(f'ProteomeLogisticRegression: {ProtRegression}')

    # --- step 1: load residue-level data ---
    ProtRegression.load_data(
        sep='|',
        reg_var=['AA', 'region'],
        response_var='cut_C_Rall',
        var2binarize=['cut_C_Rall', 'region'],
        mask_column='mapped_resid',
    )

    # --- step 2: run regression ---
    reg_df = ProtRegression.run()

    # --- step 3: persist results ---
    reg_outfile = os.path.join(outdir, f"regression_results_{tag}.csv")
    reg_df.to_csv(reg_outfile, index=False, sep='|')
    logger.info(f"SAVED: {reg_outfile}")
    
    logger.info(f'NORMAL TERMINATION - {time.time() - start_time:.1f} seconds')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
