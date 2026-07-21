#!/usr/bin/env python3
"""
Convert CHARMM COR and PSF files to PDB format.

Usage:
    python scripts/convert_cor_psf_to_pdb.py --cor structure.cor --psf structure.psf --output structure.pdb

This script uses MDAnalysis to read CHARMM coordinate (.cor/.crd) and topology (.psf) files
and writes them as a single PDB file that can be used with NCLEdetector analysis tools.
"""

import argparse
import os
import sys

def convert_cor_psf_to_pdb(cor_file, psf_file, output_pdb):
    """
    Convert CHARMM COR and PSF files to PDB format.
    
    Parameters
    ----------
    cor_file : str
        Path to CHARMM coordinate file (.cor or .crd)
    psf_file : str
        Path to CHARMM PSF topology file (.psf)
    output_pdb : str
        Path to output PDB file
    """
    try:
        import MDAnalysis as mda
    except ImportError:
        print("Error: MDAnalysis is required for this conversion.")
        print("Install with: pip install MDAnalysis")
        sys.exit(1)
    
    # Validate input files exist
    if not os.path.exists(cor_file):
        print(f"Error: COR file not found: {cor_file}")
        sys.exit(1)
    if not os.path.exists(psf_file):
        print(f"Error: PSF file not found: {psf_file}")
        sys.exit(1)
    
    print(f"Reading topology from: {psf_file}")
    print(f"Reading coordinates from: {cor_file}")
    
    # Load the structure with PSF topology and COR coordinates
    u = mda.Universe(psf_file, cor_file)
    
    print(f"Loaded structure with {len(u.atoms)} atoms")
    
    # Set chainID to 'A' for all atoms if not already set
    if not hasattr(u.atoms, 'chainIDs') or all(c == '' or c == 'X' for c in u.atoms.chainIDs):
        print("Setting chainID to 'A' for all atoms")
        u.add_TopologyAttr('chainIDs', ['A'] * len(u.atoms))
    
    # Write to PDB format
    print(f"Writing PDB to: {output_pdb}")
    u.atoms.write(output_pdb)
    
    print("Conversion complete!")
    return output_pdb


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Convert CHARMM COR/PSF files to PDB format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  convert_cor_psf_to_pdb --cor model.cor --psf model.psf --output model.pdb
  convert_cor_psf_to_pdb --cor 1zmr_ca.crd --psf 1zmr_ca.psf --output 1zmr_ca.pdb
        """,
    )

    parser.add_argument(
        "--cor",
        "--crd",
        type=str,
        required=True,
        help="Input CHARMM coordinate file (.cor or .crd)",
    )
    parser.add_argument(
        "--psf",
        type=str,
        required=True,
        help="Input CHARMM PSF topology file (.psf)",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Output PDB file",
    )

    args = parser.parse_args(argv)
    convert_cor_psf_to_pdb(args.cor, args.psf, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
