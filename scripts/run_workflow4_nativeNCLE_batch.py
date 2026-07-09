#!/usr/bin/env python3
from EntDetect._logging import setup_logger

if __package__ in {None, ""}:
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from scripts._cli_config import parse_args_with_config
else:
    from ._cli_config import parse_args_with_config

"""
Batch Workflow 4 helper for Step 2.

This script scans a directory of PDB files, filters to structures whose IDs are
present in a gene list, and runs scripts/run_nativeNCLE.py in parallel for each
selected structure.

It also builds the residue-level design matrix needed for downstream regression and Monte Carlo analysis by combining:
1) per-residue experimental/structural features from a residueFeatures CSV, and
2) region labels inferred from NCLE `ent_region` output columns.

Compared to run_nativeNCLE.py, this wrapper accepts --pdb_dir and --gene_list
instead of --pdb_file, then forwards the remaining supported nativeNCLE options.

Provide parameters directly as CLI flags, through `--config` JSON/YAML, or both.
When both are provided, CLI flags override config values.

Example
-------
python scripts/run_workflow4_nativeNCLE_batch.py \
  --pdb_dir /scratch/ims86/EntDetect_Datastore/user_input/proteome_structures/AF \
  --gene_list /scratch/ims86/EntDetect_Datastore/user_input/experimental_data/Gene_lists/AF/AF_0.6g_C_Rall_spa50_LiPMScov50_all_genes.txt \
  --outdir /scratch/ims86/EntDetect_Datastore/outputs/workflow4/nativeNCLE_all \
  --organism Ecoli \
  --model AF \
  --ent_detection_method 3 \
    --g_threshold 0.6 \
    --density 1.0 \
    --nproc 16 \
    --residue_features_file /scratch/ims86/EntDetect_Datastore/user_input/experimental_data/PDB_residue_features/AF/residueFeatures.csv \
    --reg_formula "cut_C_Rall ~ AA + region"
"""


def _read_gene_list(path):
    genes = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            item = line.strip()
            if not item:
                continue
            if item.startswith("#"):
                continue
            genes.append(item)
    return set(genes)


def _parse_formula(reg_formula: str):
    left, right = reg_formula.split("~", 1)
    response_var = left.strip()
    reg_vars = [v.strip() for v in right.split("+") if v.strip()]
    return response_var, reg_vars


def _parse_ent_region_to_set(ent_region_value):
    """Convert ent_region field to a set of integer residue indices."""
    if ent_region_value is None:
        return set()

    text = str(ent_region_value).strip()
    if text == "" or text.lower() == "nan":
        return set()

    region = set()
    for tok in text.split(","):
        tok = tok.strip()
        if tok == "":
            continue
        try:
            region.add(int(tok))
        except ValueError:
            continue
    return region


def _collect_ent_region_map(root_outdir, selected_gene_ids, logger):
    """Build gene -> set(mapped_resid in entangled regions) from NCLE feature files."""
    import glob
    import os
    import pandas as pd

    ent_region_map = {}
    missing_genes = []

    for gene in sorted(selected_gene_ids):
        feat_glob = os.path.join(root_outdir, gene, "Native_clustered_HQ_GE_features", "*_uent_features.csv")
        feat_files = sorted(glob.glob(feat_glob))
        if not feat_files:
            missing_genes.append(gene)
            ent_region_map[gene] = set()
            continue

        region_set = set()
        for fp in feat_files:
            try:
                df = pd.read_csv(fp, sep="|", usecols=["ent_region"])
            except Exception as exc:
                logger.warning(f"Could not read ent_region from {fp}: {exc}")
                continue

            for v in df["ent_region"].values:
                region_set.update(_parse_ent_region_to_set(v))

        ent_region_map[gene] = region_set

    if missing_genes:
        logger.warning(
            f"No NCLE feature files found for {len(missing_genes)} gene(s). "
            f"Their region labels will default to 0."
        )

    return ent_region_map


def _build_design_matrices(args, selected_gene_ids, logger):
    import os
    import pandas as pd

    if args.residue_features_file is None or args.reg_formula is None:
        logger.info("Design-matrix build skipped (provide both --residue_features_file and --reg_formula to enable).")
        return

    if not os.path.isfile(args.residue_features_file):
        raise FileNotFoundError(f"residue_features_file not found: {args.residue_features_file}")

    response_var, reg_vars = _parse_formula(args.reg_formula)
    logger.info(f"Building design matrices for formula: {response_var} ~ {' + '.join(reg_vars)}")

    if "region" not in reg_vars:
        logger.warning("Formula does not include 'region'; NCLE-derived region labels will not be used.")

    ent_region_map = _collect_ent_region_map(args.outdir, selected_gene_ids, logger)

    workflow4_root = os.path.dirname(os.path.abspath(args.outdir.rstrip(os.sep)))
    combined_outfile = args.design_matrix_file or os.path.join(workflow4_root, "residue_dataframes_workflow4.csv")

    req_cols = ["gene", "mapped_resid", "uniprot_length", response_var, *reg_vars]
    # Keep AA for downstream filters even if omitted in formula
    req_cols.extend(["AA"]) if "AA" not in req_cols else None
    # remove region from file-read requirements (it is built from ent_region)
    req_cols = [c for c in req_cols if c != "region"]

    logger.info(f"Reading residue features from: {args.residue_features_file}")
    df = pd.read_csv(args.residue_features_file, sep="|", low_memory=False)

    missing_cols = [c for c in req_cols if c not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns in residue_features_file: {missing_cols}")

    data = df[df["gene"].isin(selected_gene_ids)][req_cols].copy()
    logger.info(f"Rows after gene-list filter: {len(data)}")

    mapped_resid_num = pd.to_numeric(data["mapped_resid"], errors="coerce")
    data["mapped_resid"] = mapped_resid_num

    # Build region from NCLE ent_region sets
    region_values = []
    for gene, resid in zip(data["gene"].values, data["mapped_resid"].values):
        if pd.isna(resid):
            region_values.append(0)
            continue
        region_values.append(1 if int(resid) in ent_region_map.get(str(gene), set()) else 0)
    data["region"] = region_values

    final_cols = ["gene", "mapped_resid", "uniprot_length", *reg_vars, response_var]
    # preserve order and uniqueness
    seen = set()
    final_cols = [c for c in final_cols if not (c in seen or seen.add(c))]
    data = data[final_cols]

    os.makedirs(os.path.dirname(os.path.abspath(combined_outfile)), exist_ok=True)
    data.to_csv(combined_outfile, sep="|", index=False)

    logger.info(f"Design matrix build complete: single matrix file written to {combined_outfile}")


def _build_native_command(args, native_script, pdb_file, root_outdir, logdir):
    import os
    import sys

    pdb_name = os.path.basename(pdb_file)
    protein_id = os.path.splitext(pdb_name)[0]
    protein_outdir = os.path.join(root_outdir, protein_id)
    gene = args.gene if args.gene is not None else protein_id

    cmd = [
        sys.executable,
        native_script,
        "--pdb_file", pdb_file,
        "--outdir", protein_outdir,
        "--ID", protein_id,
        "--organism", args.organism,
        "--gene", gene,
        "--model", args.model,
        "--ent_detection_method", str(args.ent_detection_method),
        "--g_threshold", str(args.g_threshold),
        "--density", str(args.density),
        "--log_level", args.log_level,
        "--logdir", logdir,
    ]

    if args.chain is not None:
        cmd.extend(["--chain", args.chain])

    if args.cut_off is not None:
        cmd.extend(["--cut_off", str(args.cut_off)])

    if args.CG:
        cmd.append("--CG")

    if args.Calpha:
        cmd.append("--Calpha")

    return protein_id, cmd


def _run_one(job):
    import subprocess

    protein_id, cmd = job
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return protein_id, proc.returncode, proc.stdout, proc.stderr


def main(argv=None):
    import argparse
    import glob
    import logging
    import os
    import time
    from concurrent.futures import ThreadPoolExecutor, as_completed

    start_time = time.time()

    parser = argparse.ArgumentParser(
        description=(
            "Batch-run scripts/run_nativeNCLE.py over a PDB directory filtered by a gene list."
        )
    )
    parser.add_argument("--config", type=str, required=False, default=argparse.SUPPRESS,
                        help="Optional path to JSON/YAML config file. CLI flags override config values.")

    # --- required batch inputs ---
    parser.add_argument("--pdb_dir", type=str, required=False, default=argparse.SUPPRESS,
                        help="Directory containing input .pdb files")
    parser.add_argument("--gene_list", type=str, required=False, default=argparse.SUPPRESS,
                        help="Gene/accession list file (one ID per line)")
    parser.add_argument("--outdir", type=str, required=False, default=argparse.SUPPRESS,
                        help="Root output directory; each protein writes to outdir/<ID>")

    # --- parallelism / matching behavior ---
    parser.add_argument("--nproc", type=int, default=argparse.SUPPRESS,
                        help="Number of parallel nativeNCLE jobs (default: 8)")
    parser.add_argument("--allow_prefix_match", action="store_true", default=argparse.SUPPRESS,
                        help=(
                            "Allow gene IDs to match as prefix of PDB stem (useful for "
                            "filenames containing structure suffixes)."
                        ))
    parser.add_argument("--dry_run", action="store_true", default=argparse.SUPPRESS,
                        help="Print selected proteins and exit without running jobs")

    # --- forwarded run_nativeNCLE options (minus --pdb_file, which is set per PDB) ---
    parser.add_argument("--chain", type=str, default=argparse.SUPPRESS,
                        help="Chain identifier (optional)")
    parser.add_argument("--organism", type=str, default=argparse.SUPPRESS,
                        help="Organism for clustering: Ecoli | Human | Yeast")
    parser.add_argument("--gene", type=str, default=argparse.SUPPRESS,
                        help="Gene or accession value passed to run_nativeNCLE. If omitted, each protein uses its PDB stem.")
    parser.add_argument("--CG", action="store_true", dest="CG", default=argparse.SUPPRESS,
                        help="Pass --CG to run_nativeNCLE for coarse-grained C-alpha inputs")
    parser.add_argument("--Calpha", "--calpha", action="store_true", dest="Calpha", default=argparse.SUPPRESS,
                        help="Pass --Calpha to run_nativeNCLE to use C-alpha contacts")
    parser.add_argument("--g_threshold", type=float, default=argparse.SUPPRESS,
                        help="Gaussian entanglement score cutoff forwarded to run_nativeNCLE")
    parser.add_argument("--density", type=float, default=argparse.SUPPRESS,
                        help="Topoly triangulation density forwarded to run_nativeNCLE")
    parser.add_argument("--cut_off", type=float, default=argparse.SUPPRESS,
                        help="Cluster cutoff forwarded to run_nativeNCLE")
    parser.add_argument("--model", type=str, default=argparse.SUPPRESS,
                        help="Model type for HQ selection: EXP | AF")
    parser.add_argument("--ent_detection_method", type=int, default=argparse.SUPPRESS,
                        help="Entanglement detection method passed to run_nativeNCLE")

    # --- optional design-matrix build integrated into the batch workflow ---
    parser.add_argument("--residue_features_file", type=str, default=argparse.SUPPRESS,
                        help="Path to residue features CSV (e.g., .../PDB_residue_features/AF/residueFeatures.csv)")
    parser.add_argument("--reg_formula", type=str, default=argparse.SUPPRESS,
                        help="Regression formula for design matrix (e.g., 'cut_C_Rall ~ AA + region')")
    parser.add_argument("--design_matrix_file", type=str, default=argparse.SUPPRESS,
                        help="Output path for combined design matrix CSV (default: sibling of --outdir/residue_dataframes_workflow4.csv)")

    # --- logging ---
    parser.add_argument("--log_level", default=argparse.SUPPRESS, choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        help="Logging verbosity (default: INFO)")
    parser.add_argument("--logdir", type=str, default=argparse.SUPPRESS,
                        help="Directory for run log files (default: <outdir>/logs)")

    # --- unified-schema aliases (used by run_population_modeling integrated config) ---
    parser.add_argument("--batch_pdb_dir", type=str, required=False, default=argparse.SUPPRESS,
                        help="Alias for --pdb_dir when sharing one config with run_population_modeling.py")
    parser.add_argument("--batch_outdir", type=str, required=False, default=argparse.SUPPRESS,
                        help="Alias for --outdir when sharing one config with run_population_modeling.py")
    parser.add_argument("--batch_nproc", type=int, default=argparse.SUPPRESS,
                        help="Alias for --nproc")
    parser.add_argument("--batch_allow_prefix_match", action="store_true", default=argparse.SUPPRESS,
                        help="Alias for --allow_prefix_match")
    parser.add_argument("--batch_dry_run", action="store_true", default=argparse.SUPPRESS,
                        help="Alias for --dry_run")
    parser.add_argument("--batch_chain", type=str, default=argparse.SUPPRESS,
                        help="Alias for --chain")
    parser.add_argument("--batch_organism", type=str, default=argparse.SUPPRESS,
                        help="Alias for --organism")
    parser.add_argument("--batch_gene", type=str, default=argparse.SUPPRESS,
                        help="Alias for --gene")
    parser.add_argument("--batch_CG", action="store_true", default=argparse.SUPPRESS,
                        help="Alias for --CG")
    parser.add_argument("--batch_Calpha", action="store_true", default=argparse.SUPPRESS,
                        help="Alias for --Calpha")
    parser.add_argument("--batch_g_threshold", type=float, default=argparse.SUPPRESS,
                        help="Alias for --g_threshold")
    parser.add_argument("--batch_density", type=float, default=argparse.SUPPRESS,
                        help="Alias for --density")
    parser.add_argument("--batch_cut_off", type=float, default=argparse.SUPPRESS,
                        help="Alias for --cut_off")
    parser.add_argument("--batch_model", type=str, default=argparse.SUPPRESS,
                        help="Alias for --model")
    parser.add_argument("--batch_ent_detection_method", type=int, default=argparse.SUPPRESS,
                        help="Alias for --ent_detection_method")
    parser.add_argument("--batch_residue_features_file", type=str, default=argparse.SUPPRESS,
                        help="Alias for --residue_features_file")
    parser.add_argument("--batch_design_matrix_file", type=str, default=argparse.SUPPRESS,
                        help="Alias for --design_matrix_file")
    parser.add_argument("--batch_logdir", type=str, default=argparse.SUPPRESS,
                        help="Alias for --logdir")

    args = parse_args_with_config(
        parser,
        argv,
        defaults={
            "nproc": 8,
            "allow_prefix_match": False,
            "dry_run": False,
            "chain": None,
            "organism": "Ecoli",
            "gene": None,
            "CG": False,
            "Calpha": False,
            "g_threshold": 0.6,
            "density": 1.0,
            "cut_off": None,
            "model": "AF",
            "ent_detection_method": 3,
            "residue_features_file": None,
            "reg_formula": None,
            "design_matrix_file": None,
            "log_level": "INFO",
            "logdir": None,
        },
        aliases={"cg": "CG", "calpha": "Calpha"},
    )

    # Prefer unified batch_* keys when both prefixed and unprefixed entries exist.
    preferred_aliases = {
        "batch_pdb_dir": "pdb_dir",
        "batch_outdir": "outdir",
        "batch_nproc": "nproc",
        "batch_allow_prefix_match": "allow_prefix_match",
        "batch_dry_run": "dry_run",
        "batch_chain": "chain",
        "batch_organism": "organism",
        "batch_gene": "gene",
        "batch_CG": "CG",
        "batch_Calpha": "Calpha",
        "batch_g_threshold": "g_threshold",
        "batch_density": "density",
        "batch_cut_off": "cut_off",
        "batch_model": "model",
        "batch_ent_detection_method": "ent_detection_method",
        "batch_residue_features_file": "residue_features_file",
        "batch_design_matrix_file": "design_matrix_file",
        "batch_logdir": "logdir",
    }
    for batch_key, native_key in preferred_aliases.items():
        if hasattr(args, batch_key) and getattr(args, batch_key) is not None:
            setattr(args, native_key, getattr(args, batch_key))

    for field in ["pdb_dir", "gene_list", "outdir"]:
        if not hasattr(args, field) or getattr(args, field) is None:
            parser.error(f"Missing required argument: --{field} (or provide it in --config)")

    log_level = getattr(logging, args.log_level.upper(), logging.INFO)
    logdir = args.logdir if args.logdir is not None else os.path.join(args.outdir, "logs")
    os.makedirs(logdir, exist_ok=True)

    logger = setup_logger("run_workflow4_nativeNCLE_batch", outdir=logdir, ID="workflow4_batch", log_level=log_level)
    logger.info(f"args: {args}")

    # --- validation ---
    if not os.path.isdir(args.pdb_dir):
        parser.error(f"--pdb_dir does not exist or is not a directory: {args.pdb_dir}")
    if not os.path.isfile(args.gene_list):
        parser.error(f"--gene_list does not exist or is not a file: {args.gene_list}")
    if args.nproc < 1:
        parser.error("--nproc must be >= 1")

    os.makedirs(args.outdir, exist_ok=True)

    # Resolve native runner path relative to this script.
    native_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "run_nativeNCLE.py")
    if not os.path.isfile(native_script):
        parser.error(f"Could not locate run_nativeNCLE.py at expected path: {native_script}")

    gene_set = _read_gene_list(args.gene_list)
    if not gene_set:
        parser.error(f"--gene_list appears empty: {args.gene_list}")

    pdb_files = sorted(glob.glob(os.path.join(args.pdb_dir, "*.pdb")))
    if not pdb_files:
        parser.error(f"No .pdb files found in --pdb_dir: {args.pdb_dir}")

    selected = []
    skipped = 0
    for pdb_file in pdb_files:
        stem = os.path.splitext(os.path.basename(pdb_file))[0]
        if stem in gene_set:
            selected.append(pdb_file)
            continue

        if args.allow_prefix_match:
            # Accept if any gene is an exact prefix token of the filename stem.
            if any(stem.startswith(gene + "_") or stem.startswith(gene + "-") for gene in gene_set):
                selected.append(pdb_file)
                continue

        skipped += 1

    logger.info(f"Found {len(pdb_files)} pdb files")
    logger.info(f"Matched {len(selected)} structures against gene list; skipped {skipped}")

    if not selected:
        parser.error("No PDBs matched the provided gene list. Check naming conventions or use --allow_prefix_match.")

    jobs = [_build_native_command(args, native_script, pdb_file, args.outdir, logdir) for pdb_file in selected]

    if args.dry_run:
        logger.info("Dry-run mode enabled. Selected proteins:")
        for protein_id, _ in jobs:
            logger.info(f"  {protein_id}")
        logger.info(f"NORMAL TERMINATION - {time.time() - start_time:.1f} seconds")
        return 0

    logger.info(f"Launching {len(jobs)} nativeNCLE jobs with nproc={args.nproc}")

    failures = []
    completed = 0
    with ThreadPoolExecutor(max_workers=args.nproc) as executor:
        future_map = {executor.submit(_run_one, job): job[0] for job in jobs}
        for future in as_completed(future_map):
            protein_id = future_map[future]
            try:
                _, code, stdout, stderr = future.result()
            except Exception as exc:
                failures.append((protein_id, -1, "", f"internal runner error: {exc}"))
                logger.error(f"FAILED: {protein_id} (internal error)")
                continue

            completed += 1
            if code == 0:
                logger.info(f"DONE: {protein_id} ({completed}/{len(jobs)})")
            else:
                failures.append((protein_id, code, stdout, stderr))
                logger.error(f"FAILED: {protein_id} exit_code={code} ({completed}/{len(jobs)})")

    if failures:
        fail_log = os.path.join(args.outdir, "workflow4_nativeNCLE_batch_failures.log")
        with open(fail_log, "w", encoding="utf-8") as handle:
            for protein_id, code, stdout, stderr in failures:
                handle.write(f"#{'='*78}\n")
                handle.write(f"protein: {protein_id}\n")
                handle.write(f"exit_code: {code}\n")
                handle.write("--- stdout ---\n")
                handle.write(stdout or "")
                handle.write("\n--- stderr ---\n")
                handle.write(stderr or "")
                handle.write("\n")
        logger.error(f"{len(failures)} jobs failed. See: {fail_log}")
        logger.info(f"NORMAL TERMINATION WITH FAILURES - {time.time() - start_time:.1f} seconds")
        return 1

    logger.info(f"All {len(jobs)} jobs completed successfully")

    # Build the downstream design matrix from the completed per-protein outputs.
    selected_gene_ids = [job[0] for job in jobs]
    _build_design_matrices(args, selected_gene_ids, logger)

    logger.info(f"NORMAL TERMINATION - {time.time() - start_time:.1f} seconds")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
