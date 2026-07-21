import logging
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser, PDBIO
import pathlib
import os
from NCLEdetector._logging import setup_logger

class PDBcleaner:
    """
    Class to clean a PDB before entanglement analysis by: 
        1. Checking for duplicate residues
    """

    #############################################################
    def __init__(self, pdb:str, outdir:str='./', log_level:int=logging.INFO, logdir:str=None) -> None:
        """
        Load the PDB file and initate the PDBcleaner class
        """
        self.outdir = outdir
        self.logger = setup_logger('PDBcleaner', outdir=logdir if logdir is not None else outdir, log_level=log_level)

        parser = PDBParser()
        structure = parser.get_structure('protein', pdb)
        self.logger.debug(f'structure: {structure}')
        self.structure = structure
        self.pdb_filename = pathlib.Path(pdb).stem
        self.logger.debug(f'pdb_filename: {self.pdb_filename}')

        ## make a tmp directory to populate with cleaned PDBs if it doesnt already exists
        if not os.path.exists(outdir):
            os.mkdir(outdir)
    #############################################################

    #############################################################
    def remove_duplicates(self, pdb:str='None'):
        """
        Remove duplicate residues that are present in the PDB
        """
        ## load the PDB file if provided
        if pdb != 'None':
            self.logger.debug(f'Reading PDB file: {pdb}')
            parser = PDBParser()
            structure = parser.get_structure('protein', pdb)
            self.logger.debug(f'structure: {structure}')
            self.structure = structure
            self.pdb_filename = pathlib.Path(pdb).stem
            self.logger.debug(f'pdb_filename: {self.pdb_filename}')

        ## define the output pdb and directories
        clean_outdir = os.path.join(self.outdir, 'cleanPDB_tmp/')
        if not os.path.exists(clean_outdir):
            os.mkdir(clean_outdir)
            self.logger.debug(f'Made directory: {clean_outdir}')

        output_pdb = os.path.join(clean_outdir, f'{self.pdb_filename}_removed_duplicates.pdb')

        # Iterate over residues and identify disulfide bonds
        for model in self.structure:
            self.logger.debug(f'Model: {model}')

            for chain in model:
                self.logger.debug(f'    Chain: {chain}')
                residues_to_keep = []

                for residue in chain:
                    resname = residue.get_resname()
                    self.logger.debug(f'        Residue: {residue} {residue.get_id()} {resname} {residue.__repr__}')

                    ## check if there are any residues with alternate locs
                    residue_alts = False
                    for atom in residue:
                        self.logger.debug(f'            Atom: {atom} {atom.get_altloc()}')
                        if atom.get_altloc() != ' ':
                            residue_alts = True
                    self.logger.debug(f'            residue_alts: {residue_alts}')

                    ## if there are no alternate locs and this residues is not an insertion then keep it
                    hetflag, resseq, icode = residue.id
                    if icode == ' ' and residue_alts == False and hetflag == ' ':  # Keep only residues with no insertion code
                        residues_to_keep.append(residue)


                # Remove all residues from the chain, then re-add only the ones without insertion code
                for residue in list(chain):
                    chain.detach_child(residue.id)
                for residue in residues_to_keep:
                    chain.add(residue)

        # Save filtered structure
        io = PDBIO()
        io.set_structure(self.structure)
        io.save(output_pdb)
        self.logger.info(f'SAVED: {output_pdb}')

        return output_pdb
        #############################################################

    #############################################################
    def remove_incomplete(self, pdb:str='None'):
        """
        Remove incomplete residues that are present in the PDB
        """
        expected_heavy_atoms = {
            "ALA": 5,
            "ARG": 11,
            "ASN": 8,
            "ASP": 8,
            "CYS": 6,
            "GLU": 9,
            "GLN": 9,
            "GLY": 4,
            "HIS": 10,
            "ILE": 8,
            "LEU": 8,
            "LYS": 9,
            "MET": 8,
            "PHE": 11,
            "PRO": 7,
            "SER": 6,
            "THR": 7,
            "TRP": 14,
            "TYR": 12,
            "VAL": 7,
        }

        ## load the PDB file if provided
        if pdb != 'None':
            self.logger.debug(f'Reading PDB file: {pdb}')
            parser = PDBParser()
            structure = parser.get_structure('protein', pdb)
            self.logger.debug(f'structure: {structure}')
            self.structure = structure
            self.pdb_filename = pathlib.Path(pdb).stem
            self.logger.debug(f'pdb_filename: {self.pdb_filename}')

        ## define the output pdb and directories
        clean_outdir = os.path.join(self.outdir, 'cleanPDB_tmp/')
        if not os.path.exists(clean_outdir):
            os.mkdir(clean_outdir)
            self.logger.debug(f'Made directory: {clean_outdir}')

        output_pdb = os.path.join(clean_outdir, f'{self.pdb_filename}_removed_incomplete.pdb')

        # Iterate over residues and identify disulfide bonds
        for model in self.structure:
            self.logger.debug(f'Model: {model}')

            for chain in model:
                self.logger.debug(f'    Chain: {chain}')
                residues_to_keep = []

                for residue in chain:
                    resname = residue.get_resname()
                    self.logger.debug(f'        Residue: {residue} {residue.get_id()} {resname} {residue.__repr__}')

                    ## check if the residue is complete
                    residue_complete = False
                    if resname in expected_heavy_atoms:
                        num_heavy_atoms = sum(1 for atom in residue if atom.element != 'H')
                        self.logger.debug(f'            Number of heavy atoms: {num_heavy_atoms}')
                        if num_heavy_atoms == expected_heavy_atoms[resname]:
                            residue_complete = True
                            self.logger.debug(f'            Complete residue: {residue} with {num_heavy_atoms} heavy atoms, expected {expected_heavy_atoms[resname]}')
                        else:
                            self.logger.debug(f'            Incomplete residue: {residue} with {num_heavy_atoms} heavy atoms, expected {expected_heavy_atoms[resname]}')
                    else:
                        self.logger.debug(f'            Unknown residue type: {resname}, keeping it by default')
                        residue_complete = False

                    ## check if there are any residues with alternate locs
                    residue_alts = False
                    for atom in residue:
                        self.logger.debug(f'            Atom: {atom} {atom.get_altloc()}')
                        if atom.get_altloc() != ' ':
                            residue_alts = True
                    self.logger.debug(f'            residue_alts: {residue_alts}')

                    ## if there are no alternate locs and this residues is not an insertion then keep it
                    hetflag, resseq, icode = residue.id
                    if icode == ' ' and residue_alts == False and hetflag == ' ' and residue_complete:  # Keep only residues with no insertion code
                        residues_to_keep.append(residue)


                # Remove all residues from the chain, then re-add only the ones without insertion code
                for residue in list(chain):
                    chain.detach_child(residue.id)
                for residue in residues_to_keep:
                    chain.add(residue)

        # Save filtered structure
        io = PDBIO()
        io.set_structure(self.structure)
        io.save(output_pdb)
        self.logger.info(f'SAVED: {output_pdb}')

        return output_pdb
        #############################################################

#####################################################################
#####################################################################
