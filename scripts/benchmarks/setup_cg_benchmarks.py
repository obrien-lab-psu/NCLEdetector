#!/usr/bin/env python3
"""
Set up "quick and dirty" CG benchmark inputs and config files for Minimal
Workflows 1-4 using additional Temp-Quench CG trajectories from the
Failure-to-Form Native Entanglements project.

This script ONLY prepares inputs and writes config JSON files. It does not run
any analysis and does not modify any pipeline scripts. The generated configs are
consumed by the benchmark SLURM scripts under
``assets/slurm/scripts/benchmarks/``.

For each selected protein it:
  1. Converts the CG CA structure (``*_ca.cor`` + ``*_ca.psf``) to a PDB for MW1.
  2. Symlinks the real domain-definition file and secondary-structure file.
  3. Creates the datastore output directories.
  4. Writes config JSONs for MW1 (native NCLE), MW2 (OP full + last335),
     MW3 (nonnative clustering), and MW4 (MSM).

Frame windows (per the benchmark design, matching PGK):
    - Full traj analysis   : frames [0, 6667)
    - Last-window analysis : frames [6332, 6667)  # 335 frames
The source trajectories have 26667 frames; the OP runner now supports ``--end``
so no trajectory files are duplicated.

Usage
-----
    python scripts/benchmarks/setup_cg_benchmarks.py            # write everything
    python scripts/benchmarks/setup_cg_benchmarks.py --dry_run  # print actions only
"""

import argparse
import json
import os
import subprocess
import sys
import glob

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
REPO = "/storage/group/epo2/default/ims86/git_repos/NCLEdetector"
SRC_BASE = ("/storage/group/epo2/default/ims86/git_slugs/"
            "Failure-to-Form_Native_Entanglements_slug/"
            "Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics")
DOMAINS_DIR = ("/storage/group/epo2/default/ims86/git_slugs/"
               "Failure-to-Form_Native_Entanglements_slug/"
               "Simulations_of_Native_Entanglement_Misfolding/"
               "Rebuild_AllAtom_structures/DOMAINS")
DATASTORE = "/scratch/ims86/NCLEdetector_Datastore"

BENCH_IN_ROOT = f"{DATASTORE}/user_input/benchmarks"
OUT1_ROOT = f"{DATASTORE}/outputs/workflow1/benchmarks"
OUT2_ROOT = f"{DATASTORE}/outputs/workflow2/benchmarks"
CONFIG_ROOT = f"{REPO}/scripts/configs/benchmarks"

# Frame windows
FULL_START, FULL_END = 0, 6667
LAST_START, LAST_END = 6332, 6667

N_TRAJ = 50            # t0 .. t49 per protein
NPROC = 8              # CPUs used for G (also SLURM -n in the slurm scripts)

# --------------------------------------------------------------------------- #
# Selected proteins: stratified-random across the size distribution (55-619).
# dirname = <accession>_<pdb>_<chain>
# --------------------------------------------------------------------------- #
PROTEINS = [
    # dirname,            accession, pdb,   chain, size
    ("P61175_6XZ7_S",     "P61175",  "6XZ7", "S",  110),
    ("P0A6T9_3A7L_A",     "P0A6T9",  "3A7L", "A",  129),
    ("P45748_1HRU_A",     "P45748",  "1HRU", "A",  190),
    ("P0A6L2_1DHP_A",     "P0A6L2",  "1DHP", "A",  292),
    ("P0AD61_4YNG_C",     "P0AD61",  "4YNG", "C",  470),
]


def _write_json(path, obj, dry_run):
    print(f"  config: {path}")
    if dry_run:
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2)


def _symlink(src, dst, dry_run):
    print(f"  link: {dst} -> {src}")
    if dry_run:
        return
    if not os.path.exists(src):
        raise FileNotFoundError(f"Source for symlink missing: {src}")
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.islink(dst) or os.path.exists(dst):
        os.remove(dst)
    os.symlink(src, dst)


def _convert_cor_psf_to_pdb(cor, psf, out_pdb, dry_run):
    print(f"  pdb:  {out_pdb}  (from {os.path.basename(cor)} + {os.path.basename(psf)})")
    if dry_run:
        return
    if os.path.exists(out_pdb):
        print("        (already exists, skipping conversion)")
        return
    # These are CHARMM CARD (.cor) CG Go-model coordinates; MDAnalysis needs the
    # format specified explicitly ('CRD'). Done inline to avoid touching the
    # shared convert utility.
    import warnings
    import MDAnalysis as mda
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u = mda.Universe(psf, cor, format="CRD")
        # Ensure the PDB chain column is 'A' so both MDAnalysis (segid) and
        # mdtraj (chainID) parsers agree on the chain during native NCLE.
        u.add_TopologyAttr("chainIDs")
        u.atoms.chainIDs = "A"
        u.atoms.write(out_pdb)


def setup_protein(dirname, accession, pdb, chain, size, dry_run):
    print(f"\n=== {dirname} (accession={accession}, size={size}) ===")

    src_setup = f"{SRC_BASE}/{dirname}/setup"
    src_quench = f"{SRC_BASE}/{dirname}/Quenching"
    psf = glob.glob(f"{src_setup}/*_ca.psf")
    cor = glob.glob(f"{src_setup}/*_ca.cor")
    sec = f"{src_setup}/secondary_struc_defs.txt"
    dom = f"{DOMAINS_DIR}/{dirname}.txt"
    if not psf or not cor:
        raise FileNotFoundError(f"Missing CA psf/cor in {src_setup}")
    psf, cor = psf[0], cor[0]
    for f in (sec, dom):
        if not os.path.exists(f):
            raise FileNotFoundError(f"Missing required input: {f}")

    # ---- benchmark input dir ---------------------------------------------- #
    bench_in = f"{BENCH_IN_ROOT}/{dirname}"
    pdb_out = f"{bench_in}/{accession}_ca.pdb"
    dom_link = f"{bench_in}/domain_def.dat"
    sec_link = f"{bench_in}/secondary_struc_defs.txt"
    if not dry_run:
        os.makedirs(bench_in, exist_ok=True)
    _convert_cor_psf_to_pdb(cor, psf, pdb_out, dry_run)
    _symlink(dom, dom_link, dry_run)
    _symlink(sec, sec_link, dry_run)

    # ---- output dirs ------------------------------------------------------ #
    out1 = f"{OUT1_ROOT}/{dirname}"
    op_full = f"{OUT2_ROOT}/{dirname}/OP"
    op_last = f"{OUT2_ROOT}/{dirname}/OP_last335"
    clust = f"{OUT2_ROOT}/{dirname}/nonnative_clustering"
    msm = f"{OUT2_ROOT}/{dirname}/MSM"
    if not dry_run:
        for d in (out1, op_full, op_last, clust, msm):
            os.makedirs(f"{d}/logs", exist_ok=True)

    # ---- configs ---------------------------------------------------------- #
    cfg_dir = f"{CONFIG_ROOT}/{dirname}"

    # MW1: native NCLE (CG structure)
    # The CG Go-model structures use segid 'A' regardless of the original
    # crystallographic chain, so the structure chain to process is 'A'.
    _write_json(f"{cfg_dir}/nativeNCLE.json", {
        "pdb_file": pdb_out,
        "outdir": out1,
        "ID": accession,
        "chain": "A",
        "organism": "Ecoli",
        "gene": accession,
        "model": "EXP",
        "CG": True,
        "Calpha": True,
        "ent_detection_method": 2,
        "g_threshold": 0.6,
        "density": 1.0,
        "cut_off": None,
        "log_level": "INFO",
        "logdir": f"{out1}/logs",
    }, dry_run)

    # Shared OP settings (DCD + Traj supplied per array task on the CLI)
    op_common = {
        "ID": accession,
        "PSF": psf,
        "COR": cor,
        "sec_elements": sec_link,
        "domain": dom_link,
        "ops": ["Q", "G", "K"],
        "CG": True,
        "Calpha": True,
        "nproc": NPROC,
        "log_level": "INFO",
    }

    # MW2 full: frames [0, 6667)
    _write_json(f"{cfg_dir}/OP_full.json", {
        **op_common,
        "outdir": op_full,
        "logdir": f"{op_full}/logs",
        "start": FULL_START,
        "end": FULL_END,
        "ent_detection_method": 1,
        "no_topoly": True,
        "chunk_frames": 100,
        "chunk_suffix": "_chunk",
    }, dry_run)

    # MW2 last335: frames [6332, 6667)
    _write_json(f"{cfg_dir}/OP_last335.json", {
        **op_common,
        "outdir": op_last,
        "logdir": f"{op_last}/logs",
        "start": LAST_START,
        "end": LAST_END,
        "ent_detection_method": 2,
    }, dry_run)

    # MW3: nonnative clustering (trajnum2file built by the slurm before running)
    _write_json(f"{cfg_dir}/nonnative_clustering.json", {
        "outdir": clust,
        "trajnum2pklfile_path": f"{bench_in}/OP_last335_trajnum2file.txt",
        "traj_dir_prefix": f"{SRC_BASE}/{dirname}/Quenching",
        "start_frame": 0,
        "end_frame": 9999999,
        "nproc": 4,
        "log_level": "INFO",
        "logdir": f"{clust}/logs",
    }, dry_run)

    # MW4: MSM from the full-traj OP output; no known mirror trajectories
    _write_json(f"{cfg_dir}/MSM.json", {
        "outdir": msm,
        "ID": f"{accession}_prod",
        "OPpath": f"{op_full}/",
        "start": 0,
        "n_large_states": 10,
        "lagtime": 20,
        "rm_traj_list": [],
        "log_level": "INFO",
        "logdir": f"{msm}/logs",
    }, dry_run)


def main():
    ap = argparse.ArgumentParser(description="Set up CG benchmark inputs/configs for MW1-4.")
    ap.add_argument("--dry_run", action="store_true", help="Print actions without writing files.")
    args = ap.parse_args()

    print(f"Selected {len(PROTEINS)} proteins x {N_TRAJ} trajectories.")
    print(f"Frame windows: full=[{FULL_START},{FULL_END})  last335=[{LAST_START},{LAST_END})")
    for row in PROTEINS:
        setup_protein(*row, dry_run=args.dry_run)
    print("\nDONE. Configs under:", CONFIG_ROOT)


if __name__ == "__main__":
    main()
