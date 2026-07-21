#!/usr/bin/env python3
from NCLEdetector.statistics import ProteomeLogisticRegression
from NCLEdetector._logging import setup_logger

if __package__ in {None, ""}:
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from scripts._cli_config import parse_args_with_config
    import scripts.run_workflow4_nativeNCLE_batch as workflow4_batch
else:
    from ._cli_config import parse_args_with_config
    from . import run_workflow4_nativeNCLE_batch as workflow4_batch

"""
Run Workflow 4 proteome-level logistic regression from residue feature tables.

Provide parameters directly as CLI flags, through `--config` JSON/YAML, or both.
When both are provided, CLI flags override config values.

This script wraps NCLEdetector.statistics.ProteomeLogisticRegression and writes
the final regression table as a pipe-delimited CSV.

If batch-native-NCLE arguments are provided, the script first calls
scripts/run_workflow4_nativeNCLE_batch.py to generate the per-protein native
NCLE outputs and combined residue-level design matrix before running the
regression.

Example
-------
python scripts/run_population_modeling.py \
    --dataframe_files /path/to/residue_dataframes_workflow4.csv \
    --outdir /path/to/workflow4/population_modeling/ \
    --gene_list /path/to/gene_list.txt \
    --ID Ecoli_population \
    --reg_formula "cut_C_Rall ~ AA + region"

Expected input schema
---------------------
- Either a single combined design matrix file OR one residue table per protein
- Files are pipe-delimited ("|")
- Required columns: gene, mapped_resid, AA, region, cut_C_Rall
"""


def _default_batch_dataframe_path(batch_outdir: str) -> str:
    import os

    workflow4_root = os.path.dirname(os.path.abspath(batch_outdir.rstrip(os.sep)))
    return os.path.join(workflow4_root, "residue_dataframes_workflow4.csv")


def _build_batch_argv(args):
    batch_argv = [
        "--pdb_dir", args.batch_pdb_dir,
        "--gene_list", args.gene_list,
        "--outdir", args.batch_outdir,
        "--organism", args.batch_organism,
        "--model", args.batch_model,
        "--ent_detection_method", str(args.batch_ent_detection_method),
        "--g_threshold", str(args.batch_g_threshold),
        "--density", str(args.batch_density),
        "--nproc", str(args.batch_nproc),
        "--reg_formula", args.reg_formula,
        "--log_level", args.log_level,
    ]

    design_matrix_file = args.batch_design_matrix_file or args.dataframe_files
    if design_matrix_file is not None:
        batch_argv.extend(["--design_matrix_file", design_matrix_file])
    if args.batch_residue_features_file is not None:
        batch_argv.extend(["--residue_features_file", args.batch_residue_features_file])
    if args.batch_chain is not None:
        batch_argv.extend(["--chain", args.batch_chain])
    if args.batch_gene is not None:
        batch_argv.extend(["--gene", args.batch_gene])
    if args.batch_cut_off is not None:
        batch_argv.extend(["--cut_off", str(args.batch_cut_off)])
    if args.batch_logdir is not None:
        batch_argv.extend(["--logdir", args.batch_logdir])
    elif args.logdir is not None:
        batch_argv.extend(["--logdir", args.logdir])

    if args.batch_allow_prefix_match:
        batch_argv.append("--allow_prefix_match")
    if args.batch_dry_run:
        batch_argv.append("--dry_run")
    if args.batch_CG:
        batch_argv.append("--CG")
    if args.batch_Calpha:
        batch_argv.append("--Calpha")

    return batch_argv


def _run_batch_preprocessing(args, logger, parser):
    if not args.run_batch_native_ncle:
        return

    required_batch_fields = ["batch_pdb_dir", "batch_outdir"]
    for field in required_batch_fields:
        if getattr(args, field) is None:
            parser.error(f"Missing required argument: --{field} (or provide it in --config) when --run_batch_native_ncle is enabled")

    batch_argv = _build_batch_argv(args)
    logger.info("Launching workflow4 batch preprocessing before regression")
    logger.info(f"workflow4 batch args: {batch_argv}")
    batch_rc = workflow4_batch.main(batch_argv)
    if batch_rc != 0:
        raise RuntimeError(f"run_workflow4_nativeNCLE_batch.py failed with exit code {batch_rc}")


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
        description="Run Workflow 4 proteome-level logistic regression from residue feature tables."
    )
    parser.add_argument("--config", type=str, required=False, default=argparse.SUPPRESS,
                        help="Optional path to JSON/YAML config file. CLI flags override config values.")

    # --- required IO ---
    parser.add_argument("--dataframe_files", type=str, required=False, default=argparse.SUPPRESS,
                        help="Input design matrix path: either a directory of per-protein files or a single combined CSV")
    parser.add_argument("--outdir", type=str, required=False, default=argparse.SUPPRESS,
                        help="Output directory for regression results")
    parser.add_argument("--gene_list", type=str, required=False, default=argparse.SUPPRESS,
                        help="Path to gene list file (one ID per line)")
    parser.add_argument("--ID", type=str, required=False, default=argparse.SUPPRESS,
                        help="Identifier used for output naming")

    # --- model options ---
    parser.add_argument("--reg_formula", type=str, default=argparse.SUPPRESS,
                        help="Regression formula (default: 'cut_C_Rall ~ AA + region')")
    parser.add_argument("--sep", type=str, default=argparse.SUPPRESS,
                        help="Column separator used when loading residue data (default: '|')")
    parser.add_argument("--reg_var", nargs='+', default=argparse.SUPPRESS,
                        help="Predictor columns used for regression model encoding (default: AA region)")
    parser.add_argument("--response_var", type=str, default=argparse.SUPPRESS,
                        help="Response column name used in regression loading (default: cut_C_Rall)")
    parser.add_argument("--var2binarize", nargs='+', default=argparse.SUPPRESS,
                        help="Columns binarized during preprocessing (default: cut_C_Rall region)")
    parser.add_argument("--mask_column", type=str, default=argparse.SUPPRESS,
                        help="Residue index column used as mask/filter key (default: mapped_resid)")

    # --- optional batch native-NCLE preprocessing ---
    parser.add_argument("--run_batch_native_ncle", action='store_true', default=argparse.SUPPRESS,
                        help="Run scripts/run_workflow4_nativeNCLE_batch.py before regression to generate native NCLE outputs and the combined residue design matrix")
    parser.add_argument("--batch_pdb_dir", type=str, default=argparse.SUPPRESS,
                        help="Input PDB directory forwarded to run_workflow4_nativeNCLE_batch.py --pdb_dir")
    parser.add_argument("--batch_outdir", type=str, default=argparse.SUPPRESS,
                        help="Native-NCLE batch output directory forwarded to run_workflow4_nativeNCLE_batch.py --outdir")
    parser.add_argument("--batch_nproc", type=int, default=argparse.SUPPRESS,
                        help="Parallel job count forwarded to run_workflow4_nativeNCLE_batch.py --nproc")
    parser.add_argument("--batch_allow_prefix_match", action='store_true', default=argparse.SUPPRESS,
                        help="Forward --allow_prefix_match to run_workflow4_nativeNCLE_batch.py")
    parser.add_argument("--batch_dry_run", action='store_true', default=argparse.SUPPRESS,
                        help="Forward --dry_run to run_workflow4_nativeNCLE_batch.py")
    parser.add_argument("--batch_chain", type=str, default=argparse.SUPPRESS,
                        help="Chain identifier forwarded to run_workflow4_nativeNCLE_batch.py --chain")
    parser.add_argument("--batch_gene", type=str, default=argparse.SUPPRESS,
                        help="Gene override forwarded to run_workflow4_nativeNCLE_batch.py --gene")
    parser.add_argument("--batch_CG", action='store_true', default=argparse.SUPPRESS,
                        help="Forward --CG to run_workflow4_nativeNCLE_batch.py")
    parser.add_argument("--batch_Calpha", action='store_true', default=argparse.SUPPRESS,
                        help="Forward --Calpha to run_workflow4_nativeNCLE_batch.py")
    parser.add_argument("--batch_g_threshold", type=float, default=argparse.SUPPRESS,
                        help="Forward --g_threshold to run_workflow4_nativeNCLE_batch.py")
    parser.add_argument("--batch_density", type=float, default=argparse.SUPPRESS,
                        help="Forward --density to run_workflow4_nativeNCLE_batch.py")
    parser.add_argument("--batch_cut_off", type=float, default=argparse.SUPPRESS,
                        help="Forward --cut_off to run_workflow4_nativeNCLE_batch.py")
    parser.add_argument("--batch_model", type=str, default=argparse.SUPPRESS,
                        help="Forward --model to run_workflow4_nativeNCLE_batch.py")
    parser.add_argument("--batch_ent_detection_method", type=int, default=argparse.SUPPRESS,
                        help="Forward --ent_detection_method to run_workflow4_nativeNCLE_batch.py")
    parser.add_argument("--batch_residue_features_file", type=str, default=argparse.SUPPRESS,
                        help="Forward --residue_features_file to run_workflow4_nativeNCLE_batch.py")
    parser.add_argument("--batch_design_matrix_file", type=str, default=argparse.SUPPRESS,
                        help="Forward --design_matrix_file to run_workflow4_nativeNCLE_batch.py")
    parser.add_argument("--batch_organism", type=str, default=argparse.SUPPRESS,
                        help="Forward --organism to run_workflow4_nativeNCLE_batch.py")
    parser.add_argument("--batch_logdir", type=str, default=argparse.SUPPRESS,
                        help="Forward --logdir to run_workflow4_nativeNCLE_batch.py")

    # --- logging ---
    parser.add_argument("--log_level", default=argparse.SUPPRESS, choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        help="Logging verbosity (default: INFO)")
    parser.add_argument("--logdir", type=str, default=argparse.SUPPRESS,
                        help="Directory for log files (default: same as --outdir)")

    args = parse_args_with_config(
        parser,
        argv,
        defaults={
            "reg_formula": 'cut_C_Rall ~ AA + region',
            "sep": '|',
            "reg_var": ['AA', 'region'],
            "response_var": 'cut_C_Rall',
            "var2binarize": ['cut_C_Rall', 'region'],
            "mask_column": 'mapped_resid',
            "run_batch_native_ncle": False,
            "batch_pdb_dir": None,
            "batch_outdir": None,
            "batch_nproc": 8,
            "batch_allow_prefix_match": False,
            "batch_dry_run": False,
            "batch_chain": None,
            "batch_gene": None,
            "batch_CG": False,
            "batch_Calpha": False,
            "batch_g_threshold": 0.6,
            "batch_density": 1.0,
            "batch_cut_off": None,
            "batch_model": "AF",
            "batch_ent_detection_method": 3,
            "batch_residue_features_file": None,
            "batch_design_matrix_file": None,
            "batch_organism": "Ecoli",
            "batch_logdir": None,
            "log_level": "INFO",
            "logdir": None,
        },
    )

    args.reg_var = _normalize_list_arg(args.reg_var)
    args.var2binarize = _normalize_list_arg(args.var2binarize)

    if args.run_batch_native_ncle and getattr(args, "dataframe_files", None) is None:
        if args.batch_design_matrix_file is not None:
            args.dataframe_files = args.batch_design_matrix_file
        elif args.batch_outdir is not None:
            args.dataframe_files = _default_batch_dataframe_path(args.batch_outdir)

    for field in ["dataframe_files", "outdir", "gene_list", "ID"]:
        if not hasattr(args, field) or getattr(args, field) is None:
            parser.error(f"Missing required argument: --{field} (or provide it in --config)")

    dataframe_files = args.dataframe_files
    outdir = args.outdir
    gene_list = args.gene_list
    ID = args.ID
    reg_formula = args.reg_formula

    log_level = getattr(logging, args.log_level.upper(), logging.INFO)
    logdir = args.logdir if args.logdir is not None else outdir
    os.makedirs(logdir, exist_ok=True)

    logger = setup_logger('run_population_modeling', outdir=logdir, ID=ID, log_level=log_level)
    setup_logger('ProteomeLogisticRegression', outdir=logdir, ID=ID, log_level=log_level)
    logger.info(f'args: {args}')

    _run_batch_preprocessing(args, logger, parser)

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
        ID=ID,
        reg_formula=reg_formula,
        log_level=log_level,
        logdir=logdir,
    )
    logger.info(f'ProteomeLogisticRegression: {ProtRegression}')

    # --- step 1: load residue-level data ---
    ProtRegression.load_data(
        sep=args.sep,
        reg_var=args.reg_var,
        response_var=args.response_var,
        var2binarize=args.var2binarize,
        mask_column=args.mask_column,
    )

    # --- step 2: run regression ---
    reg_df = ProtRegression.run()

    # --- step 3: persist results ---
    reg_outfile = os.path.join(outdir, f"regression_results_{ID}.csv")
    reg_df.to_csv(reg_outfile, index=False, sep='|')
    logger.info(f"SAVED: {reg_outfile}")
    
    logger.info(f'NORMAL TERMINATION - {time.time() - start_time:.1f} seconds')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
