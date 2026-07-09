from EntDetect.order_params import CalculateOP
from EntDetect._logging import setup_logger

"""
Calculate any combination of order parameters on CG and/or all-atom trajectories.

Available OPs:  Q  G  K  SASA  XP
  Q    — fraction of native contacts
  G    — fraction of native contacts with a change of entanglement (+ entanglement features)
  K    — mirror symmetry order parameter
  SASA — solvent accessible surface area  (requires all-atom trajectory)
    XP   — Jwalk cross-link probability     (requires all-atom trajectory + --pdb_file)

For SASA/XP the all-atom trajectory is used; supply the AA topology and DCD as
--PSF and --DCD, and leave --CG / --Calpha unset.

Examples
--------
CG — Q, G, K only:
  python scripts/run_OP_on_simulation_traj.py \\
    --Traj 420 --ID 1ZMR \\
    --PSF  $REFSTRUCT/1zmr_model_clean_ca.psf \\
    --COR  $REFSTRUCT/1zmr_model_clean_ca.cor \\
    --DCD  $DATASTORE/cg_trajectories/420_prod.dcd \\
    --sec_elements $REFSTRUCT/secondary_struc_defs.txt \\
    --domain       $REFSTRUCT/domain_def.dat \\
    --outdir $DATASTORE/outputs/OP_demo \\
        --CG --Calpha \\
    --ops Q G K

AA trajectory — SASA and XP only:
  python scripts/run_OP_on_simulation_traj.py \\
    --Traj 420 --ID 1ZMR \\
    --PSF    $REFSTRUCT/1zmr_model_clean.pdb \\
    --DCD    $DATASTORE/aa_trajectories/420_prod_aa.dcd \\
    --outdir $DATASTORE/outputs/OP_demo_AA \\
    --ops SASA XP \\
    --pdb_file $REFSTRUCT/1zmr_model_clean.pdb

Flags
-----
    --config               Optional path to JSON/YAML config file. CLI flags override config values.
  --Traj                 Trajectory number (used in output filenames)
  --ID                   Base name for output files
  --PSF                  Topology file (CG PSF or AA PDB)
  --DCD                  DCD trajectory
  --outdir               Output directory (default: ./)
  --start                First frame index, 0-based (default: 0)
  --ops                  OPs to compute: Q G K SASA XP (default: Q G K)
    --CG                   Flag: input trajectory is coarse-grained
    --Calpha               Flag: use C-alpha contacts when computing G
  --ent_detection_method ENT detection: 1=GLN, 2=TLN (default), 3=GLN+TLN same termini
  --no_topoly            Disable topoly; use GLN-only workflow
  --nproc                CPU cores for G calculation (default: 10)
  --COR                  CG COR reference coordinates (required for Q, G, K)
  --sec_elements         STRIDE secondary structure definitions (required for Q, G, K)
  --domain               Domain boundary definitions (required for Q, G, K)
    --pdb_file             All-atom PDB for XP cross-link probability (required for XP)
  --chunk_frames         Frames per chunk for Combined_GE (default: None = single file)
  --chunk_suffix         Naming suffix for chunk files (default: _chunk)
  --log_level            Logging verbosity: DEBUG, INFO, WARNING, ERROR (default: INFO)
  --logdir               Directory for log file (default: same as --outdir)
"""

_CG_OPS  = {'Q', 'G', 'K'}
_AA_OPS  = {'SASA', 'XP'}
_ALL_OPS = _CG_OPS | _AA_OPS


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
        description="Calculate order parameters on CG and/or all-atom trajectories.")
    parser.add_argument("--config", type=str, required=False, default=argparse.SUPPRESS,
                        help="Optional path to JSON/YAML config file. CLI flags override config values.")
    # --- identity / IO ---
    parser.add_argument("--Traj",    type=str, required=False, default=argparse.SUPPRESS, help="Trajectory number (used in output filenames)")
    parser.add_argument("--ID",      type=str, required=False, default=argparse.SUPPRESS, help="Base name for output files")
    parser.add_argument("--PSF",     type=str, required=False, default=argparse.SUPPRESS, help="Topology file (CG PSF or AA PDB)")
    parser.add_argument("--DCD",     type=str, required=False, default=argparse.SUPPRESS, help="DCD trajectory")
    parser.add_argument("--outdir",  type=str, default=argparse.SUPPRESS, help="Output directory (default: ./)")
    parser.add_argument("--start",   type=int, default=argparse.SUPPRESS, help="First frame index, 0-based (default: 0)")
    parser.add_argument("--end",     type=int, default=argparse.SUPPRESS, help="Last frame index, exclusive (default: end of trajectory)")

    # --- which OPs ---
    parser.add_argument("--ops", nargs='+', default=argparse.SUPPRESS, choices=['Q', 'G', 'K', 'SASA', 'XP'], help="Order parameters to compute (default: Q G K)")
    
    # --- trajectory settings ---
    parser.add_argument("--CG", "--cg", action='store_true', dest="CG", default=argparse.SUPPRESS, help="Indicate trajectory is coarse-grained")
    parser.add_argument("--Calpha", "--calpha", action='store_true', dest="Calpha", default=argparse.SUPPRESS, help="Use C-alpha contacts when computing G")
    parser.add_argument("--ent_detection_method", type=int, default=argparse.SUPPRESS, help="ENT detection: 1=GLN, 2=TLN (default), 3=GLN+TLN same termini")
    parser.add_argument("--no_topoly", action="store_true", default=argparse.SUPPRESS, help="Disable topoly crossing detection (uses GLN-only workflow)")
    parser.add_argument("--nproc",  type=int, default=argparse.SUPPRESS, help="CPU cores for G (default: 10)")

    # --- CG-specific inputs (required for Q/G/K) ---
    parser.add_argument("--COR",          type=str, default=argparse.SUPPRESS, help="CG COR reference coordinates")
    parser.add_argument("--sec_elements", type=str, default=argparse.SUPPRESS, help="STRIDE secondary structure definitions file")
    parser.add_argument("--domain",       type=str, default=argparse.SUPPRESS, help="Domain boundary definitions file")

    parser.add_argument("--pdb_file",   type=str, default=argparse.SUPPRESS, help="All-atom PDB for XP (required for XP)")

    # --- G chunking (for large trajectories) ---
    parser.add_argument("--chunk_frames", type=int, default=argparse.SUPPRESS, help="Frames per chunk for Combined_GE output (default: None = single file)")
    parser.add_argument("--chunk_suffix", type=str, default=argparse.SUPPRESS, help="Naming suffix for chunked files (default: _chunk)")

    # --- logging ---
    parser.add_argument("--log_level", default=argparse.SUPPRESS, choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging verbosity (default: INFO)")
    parser.add_argument("--logdir", type=str, default=argparse.SUPPRESS, help="Directory for log file (default: same as --outdir)")

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

    # Support lowercase aliases in config files while preserving current CLI names.
    aliases = {
        "traj": "Traj",
        "id": "ID",
        "psf": "PSF",
        "dcd": "DCD",
        "cor": "COR",
        "cg": "CG",
        "calpha": "Calpha",
        "xp_pdb": "pdb_file",
    }
    for old_key, new_key in aliases.items():
        if old_key in config_args and new_key not in config_args:
            config_args[new_key] = config_args.pop(old_key)

    # Backward compatibility for older workflow2 configs.
    if "CG" not in config_args and "resolution" in config_args:
        config_args["CG"] = str(config_args.pop("resolution")).lower() == "cg"
    if "Calpha" not in config_args and "contacts" in config_args:
        config_args["Calpha"] = str(config_args.pop("contacts")).lower() == "calpha"

    merged = {
        "Traj": None,
        "ID": None,
        "PSF": None,
        "DCD": None,
        "outdir": "./",
        "start": 0,
        "end": 99999999999999,
        "ops": ["Q", "G", "K"],
        "CG": None,
        "Calpha": None,
        "ent_detection_method": 1,
        "no_topoly": False,
        "nproc": 10,
        "COR": None,
        "sec_elements": None,
        "domain": None,
        "pdb_file": None,
        "chunk_frames": None,
        "chunk_suffix": "_chunk",
        "log_level": "INFO",
        "logdir": None,
    }
    merged.update(config_args)
    merged.update(cli_args)

    if merged["Traj"] is None:
        parser.error("Missing required argument: --Traj (or provide 'Traj'/'traj' in --config)")
    if merged["ID"] is None:
        parser.error("Missing required argument: --ID (or provide 'ID'/'id' in --config)")
    if merged["PSF"] is None:
        parser.error("Missing required argument: --PSF (or provide 'PSF'/'psf' in --config)")
    if merged["DCD"] is None:
        parser.error("Missing required argument: --DCD (or provide 'DCD'/'dcd' in --config)")

    if isinstance(merged["ops"], str):
        merged["ops"] = [merged["ops"]]
    invalid_ops = set(merged["ops"]) - _ALL_OPS
    if invalid_ops:
        parser.error(f"Invalid --ops entries: {sorted(invalid_ops)}. Allowed values: {sorted(_ALL_OPS)}")

    if str(merged["log_level"]).upper() not in {"DEBUG", "INFO", "WARNING", "ERROR"}:
        parser.error("log_level must be one of: DEBUG, INFO, WARNING, ERROR")

    args = argparse.Namespace(**merged)

    ops = set(args.ops)
    traj   = args.Traj
    ID     = args.ID
    outdir = args.outdir
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    # --- resolve trajectory flags ---
    CG       = args.CG if args.CG is not None else not bool(ops & _AA_OPS)
    Calpha   = args.Calpha if args.Calpha is not None else CG
    topoly   = not args.no_topoly

    log_level = getattr(logging, args.log_level.upper(), logging.INFO)
    log_id    = f"{ID}_Traj{traj}"
    logdir    = args.logdir if args.logdir is not None else outdir

    # Pre-configure all EntDetect loggers so they share one log file
    logger = setup_logger('run_OP', outdir=logdir, ID=log_id, log_level=log_level)
    for _cls in ['CalculateOP', 'GaussianEntanglement']:
        setup_logger(_cls, outdir=logdir, ID=log_id, log_level=log_level)
    logger.info(f'args: {args}')
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    # --- input validation ---
    if ops & _CG_OPS and not all([args.COR, args.sec_elements, args.domain]):
        parser.error("--COR, --sec_elements, and --domain are required when computing Q, G, or K.")

    if ops & _AA_OPS:
        if CG:
            parser.error("SASA and XP require an all-atom trajectory: do not set --CG.")
        if 'XP' in ops and args.pdb_file is None:
            parser.error("--pdb_file is required when XP is in --ops.")
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    # --- instantiate CalculateOP for primary (CG or AA) trajectory ---
    CalcOP = CalculateOP(outdir=outdir,
                            Traj=traj,
                            ID=ID,
                            psf=args.PSF,
                            cor=args.COR,
                            sec_elements=args.sec_elements,
                            dcd=args.DCD,
                            domain=args.domain,
                            start=args.start,
                            end=args.end,
                            ent_detection_method=args.ent_detection_method,
                            log_level=log_level,
                            logdir=logdir)
    logger.info(f'CalculateOP (primary): {CalcOP}')

    if 'Q' in ops:
        Qdata_dict = CalcOP.Q()
        logger.info(f'Q keys: {list(Qdata_dict.keys())}')

    if 'G' in ops:
        Gdata_dict = CalcOP.G(topoly=topoly, Calpha=Calpha, CG=CG, nproc=args.nproc, chunk_frames=args.chunk_frames, chunk_suffix=args.chunk_suffix)
        logger.info(f'G keys: {list(Gdata_dict.keys())}')

    if 'K' in ops:
        Kdata_dict = CalcOP.K()
        logger.info(f'K keys: {list(Kdata_dict.keys())}')

    if 'SASA' in ops:
        SASAdata_dict = CalcOP.SASA()
        logger.info(f'SASA keys: {list(SASAdata_dict.keys())}')

    if 'XP' in ops:
        XPdata_dict = CalcOP.XP(pdb_file=args.pdb_file, use_traj=True, nproc=args.nproc)
        logger.info(f'XP keys: {list(XPdata_dict.keys())}')
    ###---------------------------------------------------------------------------------------------------------

    logger.info(f'NORMAL TERMINATION - {time.time() - start_time:.1f} seconds')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
