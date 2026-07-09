#!/usr/bin/env python3
from EntDetect.gaussian_entanglement import GaussianEntanglement
from EntDetect.clustering import ClusterNativeEntanglements
from EntDetect.entanglement_features import FeatureGen
from EntDetect._logging import setup_logger

"""
Script to calculate native Gaussian entanglements in a given structure (PDB or COR file),
filter for high-quality entanglements, cluster them, and generate entanglement features.

Usage example (1ZMR / ecPGK):
    python scripts/run_nativeNCLE.py \\
        --pdb_file  /scratch/ims86/EntDetect_Datastore/user_input/reference_structures/1zmr_model_clean.pdb \\
        --outdir  /scratch/ims86/EntDetect_Datastore/outputs/workflow1 \\
        --ID      1ZMR \\
        --chain   A \\
        --organism Ecoli \\
        --gene P00558 \\
        --model   EXP

Arguments:
    --config              Optional path to JSON/YAML config file containing any
                          script arguments. CLI flags override config values.
    --pdb_file              Path to input PDB (or COR) structure file                    [required]
    --outdir              Root output directory; sub-dirs are created automatically    [required]
    --ID                  Identifier for the analysis (default: structure basename)
    --chain               Chain ID to process; omit to process all chains
    --organism            Reference proteome for clustering: Ecoli | Human | Yeast     (default: Ecoli)
    --gene           UniProt accession used in feature-file naming                (default: P00558)
    --model               Structure type for HQ filtering: EXP | AF                   (default: EXP)
    --CG                  Flag: input is a coarse-grained C-alpha model
    --Calpha              Flag: use C-alpha atoms for contact definition
    --cut_off      Clustering distance cutoff in Å; if omitted, uses the
                          organism-specific default (Ecoli: 57, Human: 52, Yeast: 49)
    --ent_detection_method
                          Entanglement detection criterion:
                            1 = any nonzero GLN for either termini
                            2 = any nonzero TLN for either termini  (class default)
                            3 = both GLN and TLN nonzero for same termini  (recommended; script default)
"""


def main(argv=None):

    import multiprocessing as mp
    import sys, os
    import argparse
    import json
    import time

    start_time = time.time()

    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--config", type=str, required=False, default=argparse.SUPPRESS, help="Optional path to JSON or YAML config file. CLI flags override config values.")
    parser.add_argument("--pdb_file", type=str, required=False, default=argparse.SUPPRESS, help="Path to PDB structure file")
    parser.add_argument("--outdir", type=str, required=False, default=argparse.SUPPRESS, help="output directory for results")
    parser.add_argument("--ID", type=str, required=False, default=argparse.SUPPRESS, help="An id for the analysis (defaults to structure basename)")
    parser.add_argument("--chain", type=str, required=False, default=argparse.SUPPRESS, help="Chain identifier (optional, processes all chains if not specified)")
    parser.add_argument("--organism", type=str, required=False, default=argparse.SUPPRESS, help="Organism name for clustering: {Ecoli, Human, Yeast}")
    parser.add_argument("--gene", type=str, required=False, default=argparse.SUPPRESS, help="UniProt Accession for the protein")
    parser.add_argument("--CG", "--cg", action='store_true', dest="CG", default=argparse.SUPPRESS, help="Indicate structure is coarse-grained (C-alpha only) model")
    parser.add_argument("--Calpha", "--calpha", action='store_true', dest="Calpha", default=argparse.SUPPRESS, help="Use C-alpha atoms for contact definition")
    parser.add_argument("--g_threshold", type=float, required=False, default=argparse.SUPPRESS, help="Gaussian entanglement score cutoff used when rounding GLN values")
    parser.add_argument("--density", type=float, required=False, default=argparse.SUPPRESS, help="Triangulation density used for minimal loop surface generation")
    parser.add_argument("--cut_off", type=float, required=False, default=argparse.SUPPRESS, help="Clustering distance cutoff in Å. If omitted, uses the organism-specific default (Ecoli: 57, Human: 52, Yeast: 49).",)
    parser.add_argument("--model", type=str, required=False, default=argparse.SUPPRESS, help="Model type for high-quality selection: {EXP, AF}")
    parser.add_argument("--ent_detection_method", type=int, required=False, default=argparse.SUPPRESS, help="ENT detection method: 1=any GLN, 2=any TLN (default), 3=both GLN and TLN same termini")
    parser.add_argument("--log_level", type=str, default=argparse.SUPPRESS, choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging verbosity level (default: INFO)")
    parser.add_argument("--logdir", type=str, default=argparse.SUPPRESS, help="Directory for log file. Defaults to --outdir if not specified.")

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

    # Merge precedence: defaults < config file < explicit CLI flags.
    merged = {
        "pdb_file": None,
        "outdir": None,
        "ID": None,
        "chain": None,
        "organism": "Ecoli",
        "gene": "P00558",
        "CG": False,
        "Calpha": False,
        "g_threshold": 0.6,
        "density": 1.0,
        "cut_off": None,
        "model": "EXP",
        "ent_detection_method": 2,
        "log_level": "INFO",
        "logdir": None,
    }
    merged.update(config_args)
    merged.update(cli_args)

    if merged["pdb_file"] is None:
        parser.error("Missing required argument: --pdb_file (or provide 'pdb_file' in --config)")
    if merged["outdir"] is None:
        parser.error("Missing required argument: --outdir (or provide 'outdir' in --config)")

    if str(merged["log_level"]).upper() not in {"DEBUG", "INFO", "WARNING", "ERROR"}:
        parser.error("log_level must be one of: DEBUG, INFO, WARNING, ERROR")

    args = argparse.Namespace(**merged)
    import logging
    log_level = getattr(logging, args.log_level.upper(), logging.INFO)
    
    pdb_file = args.pdb_file
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    ID = args.ID if args.ID is not None else os.path.splitext(os.path.basename(pdb_file))[0]
    logdir = args.logdir if args.logdir is not None else outdir
    # Pre-configure all EntDetect loggers for this run so they share one log file
    logger = setup_logger('run_nativeNCLE', outdir=logdir, ID=ID, log_level=log_level)
    for _cls in ['GaussianEntanglement', 'ClusterNativeEntanglements', 'FeatureGen']:
        setup_logger(_cls, outdir=logdir, ID=ID, log_level=log_level)
    logger.info(f'args: {args}')
    chain = args.chain
    organism = args.organism
    cut_off = args.cut_off
    model = args.model
    g_threshold = args.g_threshold
    density = args.density

    CG = bool(args.CG)
    Calpha = bool(args.Calpha)

    # Set up Gaussian Entanglement and Clustering objects
    ge = GaussianEntanglement(g_threshold=g_threshold, density=density, Calpha=Calpha, CG=CG, ent_detection_method=args.ent_detection_method, log_level=log_level, logdir=logdir)
    clustering = ClusterNativeEntanglements(organism=organism, cut_off=cut_off, log_level=log_level, logdir=logdir)

    # Determine which chains to process
    if chain is not None:
        chains_to_process = [chain]
    else:
        # Get all chains from the structure
        import MDAnalysis as mda
        u = mda.Universe(pdb_file)
        chains_to_process = sorted(set([atom.segid if atom.segid else 'A' for atom in u.atoms if atom.segid or atom.chainID]))
        if not chains_to_process or chains_to_process == ['']:
            # Fallback: use mdtraj to get chains
            import mdtraj as md
            traj = md.load(pdb_file)
            chains_to_process = sorted(set([c.chain_id for c in traj.topology.chains]))
        logger.info(f'Processing chains: {chains_to_process}')

    # Process each chain separately for all steps
    for chain_id in chains_to_process:
        logger.info(f"{'='*80}\nProcessing chain {chain_id}\n{'='*80}")
        
        # Use chain suffix for file naming when processing multiple chains
        if len(chains_to_process) > 1:
            hq_id = f"{ID}_{chain_id}"
        else:
            hq_id = ID
        
        # All chains use the same NCLE directory
        ge_outdir = os.path.join(outdir, 'NCLE')
        os.makedirs(ge_outdir, exist_ok=True)
        
        # Calculate native entanglements for this chain
        NativeEnt = ge.calculate_native_entanglements(pdb_file=pdb_file, outdir=ge_outdir, ID=hq_id, chain=chain_id)
        logger.info(f'Native entanglements saved to {NativeEnt["outfile"]}')
        
        # Optional steps: select high-quality entanglements 
        HQNativeEnt = ge.select_high_quality_entanglements(rawNCLE_file=NativeEnt['outfile'], pdb_file=pdb_file, outdir=os.path.join(outdir, "HQ_NCLE"), ID=hq_id, model=model, chain=chain_id)
        logger.info(f'High-quality native entanglements saved to {HQNativeEnt["outfile"]}')

        # Cluster the native entanglements to remove degeneracies
        nativeClusteredEnt = clustering.Cluster_NativeEntanglements(HQ_NCLE_file=HQNativeEnt['outfile'], outdir=os.path.join(outdir, "clustered_HQ_NCLE"), ID=hq_id, chain=chain_id)
        logger.info(f'Clustered native entanglements saved to {nativeClusteredEnt["outfile"]}')

        # Generate entanglement features for clustered native entanglements
        FGen = FeatureGen(log_level=log_level, logdir=logdir)
        EntFeatures = FGen.get_uent_features(pdb_file=pdb_file,
                             outdir=os.path.join(outdir, "clustered_HQ_NCLE_features"),
                             cluster_file=nativeClusteredEnt['outfile'],
                             gene=args.gene,
                             chain=chain_id,
                             pdbid=ID)
        logger.info(f'Entanglement features saved to {EntFeatures["outfile"]}')


    logger.info(f'NORMAL TERMINATION - {time.time() - start_time:.1f} seconds')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
