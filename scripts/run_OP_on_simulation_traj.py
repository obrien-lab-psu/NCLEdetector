from EntDetect.order_params import CalculateOP
from EntDetect._logging import setup_logger

"""
Calculate any combination of order parameters on CG and/or all-atom trajectories.

Available OPs:  Q  G  K  SASA  XP
  Q    — fraction of native contacts
  G    — fraction of native contacts with a change of entanglement (+ entanglement features)
  K    — mirror symmetry order parameter
  SASA — solvent accessible surface area  (requires all-atom trajectory)
  XP   — Jwalk cross-link probability     (requires all-atom trajectory + --xp_pdb)

For SASA/XP the all-atom trajectory is used; set --resolution aa and supply the
AA topology and DCD as --PSF and --DCD.

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
    --ops Q G K

AA trajectory — SASA and XP only:
  python scripts/run_OP_on_simulation_traj.py \\
    --Traj 420 --ID 1ZMR \\
    --PSF    $REFSTRUCT/1zmr_model_clean.pdb \\
    --DCD    $DATASTORE/aa_trajectories/420_prod_aa.dcd \\
    --resolution aa \\
    --outdir $DATASTORE/outputs/OP_demo_AA \\
    --ops SASA XP \\
    --xp_pdb $REFSTRUCT/1zmr_model_clean.pdb

Flags
-----
  --Traj                 Trajectory number (used in output filenames)
  --ID                   Base name for output files
  --PSF                  Topology file (CG PSF or AA PDB)
  --DCD                  DCD trajectory
  --outdir               Output directory (default: ./)
  --start                First frame index, 0-based (default: 0)
  --ops                  OPs to compute: Q G K SASA XP (default: Q G K)
  --resolution           Trajectory resolution: cg (default) or aa
  --contacts             Contact type: calpha or heavy (default: calpha for cg, heavy for aa)
  --ent_detection_method ENT detection: 1=GLN, 2=TLN (default), 3=GLN+TLN same termini
  --no_topoly            Disable topoly; use GLN-only workflow
  --nproc                CPU cores for G calculation (default: 10)
  --COR                  CG COR reference coordinates (required for Q, G, K)
  --sec_elements         STRIDE secondary structure definitions (required for Q, G, K)
  --domain               Domain boundary definitions (required for Q, G, K)
  --xp_pdb               All-atom PDB for XP cross-link probability (required for XP)
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
    import time
    import logging
    start_time = time.time()
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(
        description="Calculate order parameters on CG and/or all-atom trajectories.")
    # --- identity / IO ---
    parser.add_argument("--Traj",    type=str, required=True,  help="Trajectory number (used in output filenames)")
    parser.add_argument("--ID",      type=str, required=True,  help="Base name for output files")
    parser.add_argument("--PSF",     type=str, required=True,  help="Topology file (CG PSF or AA PDB)")
    parser.add_argument("--DCD",     type=str, required=True,  help="DCD trajectory")
    parser.add_argument("--outdir",  type=str, default='./',   help="Output directory (default: ./)")
    parser.add_argument("--start",   type=int, default=0,      help="First frame index, 0-based (default: 0)")

    # --- which OPs ---
    parser.add_argument("--ops", nargs='+', default=['Q', 'G', 'K'], choices=['Q', 'G', 'K', 'SASA', 'XP'], help="Order parameters to compute (default: Q G K)")
    
    # --- trajectory settings ---
    parser.add_argument("--resolution", choices=["cg", "aa"], default="cg", help="Trajectory resolution: cg (default) or aa")
    parser.add_argument("--contacts",   choices=["calpha", "heavy"], default=None, help="Contact type: calpha or heavy (default: calpha for cg, heavy for aa)")
    parser.add_argument("--ent_detection_method", type=int, default=1, help="ENT detection: 1=GLN, 2=TLN (default), 3=GLN+TLN same termini")
    parser.add_argument("--no_topoly", action="store_true", help="Disable topoly crossing detection (uses GLN-only workflow)")
    parser.add_argument("--nproc",  type=int, default=10, help="CPU cores for G (default: 10)")

    # --- CG-specific inputs (required for Q/G/K) ---
    parser.add_argument("--COR",          type=str, default=None, help="CG COR reference coordinates")
    parser.add_argument("--sec_elements", type=str, default=None, help="STRIDE secondary structure definitions file")
    parser.add_argument("--domain",       type=str, default=None, help="Domain boundary definitions file")

    parser.add_argument("--xp_pdb",   type=str, default=None, help="All-atom PDB for XP (required for XP)")

    # --- logging ---
    parser.add_argument("--log_level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging verbosity (default: INFO)")
    parser.add_argument("--logdir", type=str, default=None, help="Directory for log file (default: same as --outdir)")
    args = parser.parse_args(argv)

    ops = set(args.ops)
    traj   = args.Traj
    ID     = args.ID
    outdir = args.outdir
    ###---------------------------------------------------------------------------------------------------------

    ###---------------------------------------------------------------------------------------------------------
    # --- resolve derived settings ---
    contacts = args.contacts if args.contacts is not None else ("calpha" if args.resolution == "cg" else "heavy")
    Calpha   = contacts == "calpha"
    CG       = args.resolution == "cg"
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
            parser.error("SASA and XP require an all-atom trajectory: set --resolution aa.")
        if 'XP' in ops and args.xp_pdb is None:
            parser.error("--xp_pdb is required when XP is in --ops.")
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
                            ent_detection_method=args.ent_detection_method,
                            log_level=log_level,
                            logdir=logdir)
    logger.info(f'CalculateOP (primary): {CalcOP}')

    if 'Q' in ops:
        Qdata_dict = CalcOP.Q()
        logger.info(f'Q keys: {list(Qdata_dict.keys())}')

    if 'G' in ops:
        Gdata_dict = CalcOP.G(topoly=topoly, Calpha=Calpha, CG=CG, nproc=args.nproc)
        logger.info(f'G keys: {list(Gdata_dict.keys())}')

    if 'K' in ops:
        Kdata_dict = CalcOP.K()
        logger.info(f'K keys: {list(Kdata_dict.keys())}')

    if 'SASA' in ops:
        SASAdata_dict = CalcOP.SASA()
        logger.info(f'SASA keys: {list(SASAdata_dict.keys())}')

    if 'XP' in ops:
        XPdata_dict = CalcOP.XP(pdb=args.xp_pdb)
        logger.info(f'XP keys: {list(XPdata_dict.keys())}')
    ###---------------------------------------------------------------------------------------------------------

    logger.info(f'NORMAL TERMINATION - {time.time() - start_time:.1f} seconds')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
