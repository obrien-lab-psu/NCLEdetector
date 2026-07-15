import os
import math
import re
import logging
import argparse
import numpy as np
import pandas as pd
from glob import glob
import mdtraj as md
from Bio.PDB import PDBParser, is_aa
from Bio import PDB
from scipy.spatial.distance import pdist, squareform
import MDAnalysis as mda
import requests, sys
from EntDetect._logging import setup_logger
np.set_printoptions(linewidth=np.inf, precision=4)
pd.set_option('display.max_rows', None)

class FeatureGen:
    """
    Processes biological data including PDB files, sequence data, and interaction potentials.
    """
    #############################################################################################################
    def __init__(self, log_level:int=logging.INFO, logdir:str=None):
        self.logger = setup_logger('FeatureGen', outdir=logdir if logdir is not None else './', log_level=log_level)
    #############################################################################################################

    #############################################################################################################
    def get_AA(self, pdb_file, gene):

        """
        Get the PDB resid to AA mapping for the provided PDB
        """
        three_to_one_letter = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'MSE': 'M', 'PHE': 'F', 
        'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 
        'VAL': 'V'}

        resid2AA = {}
        # Define the path to your PDB file

        # Create a PDB parser
        parser = PDBParser(QUIET=True)

        # Parse the PDB file
        structure = parser.get_structure("protein", pdb_file)

        # Initialize an empty list to store amino acid codes
        amino_acid_codes = []

        # Iterate through the structure and extract amino acid codes
        for model in structure:
            for chain in model:
                for residue in chain:
                    if is_aa(residue):
                        resname = residue.get_resname()
                        resid = residue.get_id()[1]
                        if resname in three_to_one_letter:
                            AA = three_to_one_letter[resname]
                        else:
                            AA = 'NC'
                        #print(resname, resid, AA)
                        resid2AA[resid] = AA
        self.resid2AA = resid2AA

        ## get the canonical uniprot sequence length
        # Define the URL for the UniProt API
        url = f"https://rest.uniprot.org/uniprotkb/{gene}.fasta"
        
        # Make a GET request to the UniProt API
        response = requests.get(url)
        
        # Check if the response is OK
        if response.status_code == 200:
            # Extract the sequence from the FASTA format
            fasta_data = response.text.splitlines()
            sequence = ''.join(fasta_data[1:])  # Skip the first line (header)
            
            # Return the length of the sequence
            self.prot_size = len(sequence)
            if self.prot_size == 0:
                self.logger.error(f'The size of the protein in Uniprot is {self.prot_size} == 0. This likely means this uniprot ID no longer exists. No entanglement features will be calculated')
                quit()
        else:
            raise ValueError(f"Error: Could not retrieve data for UniProt ID {uniprot_id}.")
    #############################################################################################################

    #############################################################################################################
    def split_on_nth_char(self, s, char, n):
        # Find the index of the nth occurrence of the char
        occurrence = 0
        index = -1
        for i, c in enumerate(s):
            if c == char:
                occurrence += 1
            if occurrence == n:
                index = i
                break
        
        # If the nth occurrence is found, split the string
        if index != -1:
            return s[:index], s[index+1:]
        else:
            return s, ""
    #############################################################################################################

    #############################################################################################################
    def get_uent_features(self, pdb_file:str, outdir:str='./', cluster_file:str='None', gene:str='None', pdbid:str='None', chain:str='A'):
        """
        Get the features for each unique entanglement provided in the clustered_unampped_GE file
        """

        self.pdb_file = pdb_file
        self.outdir = outdir
        self.cluster_file = cluster_file

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
            self.logger.debug(f'Made directory: {self.outdir}')

        self.traj = md.load(self.pdb_file)

        if os.path.exists(self.cluster_file):
            self.GE_data = pd.read_csv(self.cluster_file, sep='|', dtype={'c': str, 'crossingsN': str, 'crossingsC': str})
        else:
            raise ValueError(f"{self.cluster_file} does not exist")

        uent_df = {'gene':[],
                    'PDB':[],
                    'chain':[], 
                    'ENT-ID':[],
                    'gn':[],
                    'N_term_thread':[],
                    'gc':[],
                    'C_term_thread':[],
                    'i':[],
                    'j':[],
                    'NC':[],
                    'NC_wbuff':[],
                    'NC_region':[],
                    'crossingsN':[], 
                    'crossingsC':[],
                    'crossingsN_wbuff':[], 
                    'crossingsC_wbuff':[],
                    'crossingsN_region':[],
                    'crossingsC_region':[],
                    'ent_region':[],
                    'loopsize': [], 
                    'num_zipper_nc':[], 
                    'perc_bb_loop':[],
                    'num_loop_contacting_res':[],
                    'num_cross_nearest_neighbors':[],
                    'ent_coverage':[],
                    'min_N_prot_depth_left':[],
                    'min_N_thread_depth_left':[],
                    'min_N_thread_slippage_left':[],
                    'min_C_prot_depth_right':[],
                    'min_C_thread_depth_right':[],
                    'min_C_thread_slippage_right':[], 
                    'prot_size':[], 
                    'ACO':[],
                    'RCO':[],
                    'CCBond':[]}

        #############################################################################################################################################################################
        ### Load entanglement information if present
        topology = self.traj.topology

        # get mapping of chain letters to chain index 
        chain_ids = {chain.chain_id: chain.index for chain in topology.chains}
        self.logger.debug(f'chain_ids: {chain_ids}')
        if chain not in chain_ids:
            raise ValueError(f'chain {chain} not in PDB file')
        

        # Get the protein size from uniprot and a dictionary that maps resid to amino acid (one letter)
        self.get_AA(self.pdb_file, gene)
        self.logger.debug(f'gene: {gene}, chain: {chain}, pdbid: {pdbid}, prot_size: {self.prot_size}')

        ## parse lines to get native contacts, crossings,
        rbuffer = 3
        pdb_NC_list = [] # list of PDB native contact residues +/- rbuffer
        pdb_NC_core_list = [] # list of PDB natvie contact residues
        pdb_crossing_list = [] # list of PDB crossing residues +/- rbuffer
        pdb_crossing_core_list = [] # list of PDB crossing residues

        for rowi, row in self.GE_data.iterrows():
            #print(row)
            #print(f'#######: ENT-ID: {rowi}')
            ent_core = []

            ## check that the entanglement isnt in a non-mapped area. if so skip it
            #line = line[1].split(',')
            pdb_NCi_core = row['i']
            pdb_NCj_core = row['j']

            # Parse crossings from crossingsN and crossingsC columns  
            # Each column contains comma-separated crossing residues like "+109" or "+92,+93,+94"
            pdb_crossing_res_core_N = []
            pdb_crossing_res_core_C = []
            for col in ['crossingsN', 'crossingsC']:
                if col in row.index and pd.notna(row[col]) and row[col] != '':
                    crossings_str = str(row[col])
                    for cross in crossings_str.split(','):
                        if cross:  # Skip empty strings
                            # Remove +/- sign and convert to int, handling potential .0 float artifacts
                            cross_num = cross[1:].split('.')[0]  # Remove sign and any decimal part
                            cross_int = int(cross_num)
                            if col == 'crossingsN':
                                pdb_crossing_res_core_N.append(cross_int)
                            else:
                                pdb_crossing_res_core_C.append(cross_int)
            
            # Combined list for backward compatibility in calculations
            pdb_crossing_res_core = pdb_crossing_res_core_N + pdb_crossing_res_core_C
            #print(f'pdb_crossing_res_core_N: {pdb_crossing_res_core_N}, pdb_crossing_res_core_C: {pdb_crossing_res_core_C}')

            uent_df['gene'] += [gene]
            uent_df['PDB'] += [pdbid]
            uent_df['chain'] += [chain]
            uent_df['ENT-ID'] += [rowi]
            uent_df['i'] += [pdb_NCi_core]
            uent_df['j'] += [pdb_NCj_core]


            #########################################################################
            ## get Gn and Gc and if it is present the cluster size
            num_zipper_nc = row['num_contacts']
            CCBond = row['CCBond']
            gn = row['gn']
            gc = row['gc']

            
            # Calcualte the absolute and relative contact orders
            range_strings = row['contacts'].split(';')
            loops = []
            for l in range_strings:
                # if no negative residue was found
                if l.count('-') == 1:
                    x = l.split('-', 1)
                elif l.count('--') == 1 and l.count('-') == 2:
                    x = l.split('-', 1)
                elif l.count('-') == 2 and l.count('--') == 0:
                    x = self.split_on_nth_char(l, '-', 2)
                elif l.count('--') == 1 and l.count('-') == 3:
                    x = self.split_on_nth_char(l, '-', 2)
                loops += [(int(x[0]), int(x[1]))]   

            loop_sizes = [j-i for i,j in loops]
            #print(f'loop_sizes: {loop_sizes}')
            ACO = np.sum(loop_sizes)/len(loop_sizes)
            RCO = ACO/self.prot_size
            #print(f'gn: {gn} | gc: {gc} | num_zipper_nc: {num_zipper_nc} | ACO: {ACO} | RCO: {RCO} | CCBond: {CCBond}')

            uent_df['gn'] += [gn]
            uent_df['gc'] += [gc]
            uent_df['num_zipper_nc'] += [num_zipper_nc]
            uent_df['ACO'] += [ACO]
            uent_df['RCO'] += [RCO]
            uent_df['CCBond'] += [CCBond]


            #########################################################################
            #get PDB native contact and those +/- rbuffer along the primary structure
            pdb_NC_core = [pdb_NCi_core, pdb_NCj_core]
            pdb_NC_core_list += pdb_NC_core

            pdb_NCi = np.arange(pdb_NCi_core - rbuffer, pdb_NCi_core + rbuffer + 1)
            pdb_NCj = np.arange(pdb_NCj_core - rbuffer, pdb_NCj_core + rbuffer + 1)
            pdb_NC = np.hstack([pdb_NCi, pdb_NCj]).tolist()
            pdb_NC_list += pdb_NC

            #print(f'pdb_NC: {pdb_NC}')
            #print(f'pdb_NC_core: {pdb_NC_core}')
            uent_df['NC'] += [",".join([str(r) for r in pdb_NC_core])]
            uent_df['NC_wbuff'] += [",".join([str(r) for r in pdb_NC])]

            ## Calculate the NC_region using heavy atom distances
            NC_region = self.find_neighboring_residues(self.traj, pdb_NC)
            #print(f'NC_region: {NC_region}')
            uent_df['NC_region'] += [",".join([str(r) for r in NC_region])]
            

            loopsize = pdb_NCj_core - pdb_NCi_core
            loop_resids = np.arange(pdb_NCi_core, pdb_NCj_core + 1)
            loop_contacting_res = self.find_neighboring_residues(self.traj, loop_resids)
            num_loop_contacting_res = len(loop_contacting_res)
            #print(f'loop_contacting_res: {loop_contacting_res}')
            #print(f'num_loop_contacting_res: {num_loop_contacting_res}')

            uent_df['loopsize'] += [loopsize]
            uent_df['perc_bb_loop'] += [loopsize/self.prot_size]
            uent_df['num_loop_contacting_res'] += [num_loop_contacting_res]
            #########################################################################


            #########################################################################
            #get PDB crossings and those +/- rbuffer along the primary structure
            if pdb_crossing_res_core_N:
                pdb_crossing_res_N = np.hstack([np.arange(int(x) - rbuffer, int(x) + rbuffer + 1) for x in pdb_crossing_res_core_N]).tolist()
            else:
                pdb_crossing_res_N = []
                
            if pdb_crossing_res_core_C:
                pdb_crossing_res_C = np.hstack([np.arange(int(x) - rbuffer, int(x) + rbuffer + 1) for x in pdb_crossing_res_core_C]).tolist()
            else:
                pdb_crossing_res_C = []
                
            # Combined for overall calculations
            pdb_crossing_res = pdb_crossing_res_N + pdb_crossing_res_C
            #print(f'pdb_crossing_res_N: {pdb_crossing_res_N}')
            #print(f'pdb_crossing_res_C: {pdb_crossing_res_C}')
            #print(f'pdb_crossing_res_core_N: {pdb_crossing_res_core_N}, pdb_crossing_res_core_C: {pdb_crossing_res_core_C}')

            pdb_crossing_list += pdb_crossing_res
            pdb_crossing_core_list += pdb_crossing_res_core
            
            # Store separated crossings
            uent_df['crossingsN'] += [",".join([str(c) for c in pdb_crossing_res_core_N])]
            uent_df['crossingsC'] += [",".join([str(c) for c in pdb_crossing_res_core_C])]
            uent_df['crossingsN_wbuff'] += [",".join([str(c) for c in pdb_crossing_res_N])]
            uent_df['crossingsC_wbuff'] += [",".join([str(c) for c in pdb_crossing_res_C])]

            ### Get the crossing region using heavy atom distances
            crossing_region_N = self.find_neighboring_residues(self.traj, pdb_crossing_res_N)
            crossing_region_C = self.find_neighboring_residues(self.traj, pdb_crossing_res_C)
            #print(f'crossing_region_N: {crossing_region_N}')
            #print(f'crossing_region_C: {crossing_region_C}')
            uent_df['crossingsN_region'] += [",".join([str(r) for r in crossing_region_N])]
            uent_df['crossingsC_region'] += [",".join([str(r) for r in crossing_region_C])]

            num_cross_nearest_neighbors = len(crossing_region_N) + len(crossing_region_C)
            #print(f'num_cross_nearest_neighbors: {num_cross_nearest_neighbors}')
            uent_df['num_cross_nearest_neighbors'] += [num_cross_nearest_neighbors]
            #########################################################################


            #########################################################################
            ## Get number of threads in each termini and depth
            #print(f'prot_size: {self.prot_size}')
            N_term_thread = [c for c in pdb_crossing_res_core if c < pdb_NCi_core]            
            num_N_term_thread = len(N_term_thread)
            #print(f'num_N_term_thread: {num_N_term_thread}')
            
            C_term_thread = [c for c in pdb_crossing_res_core if c > pdb_NCj_core]            
            num_C_term_thread = len(C_term_thread)
            #print(f'num_C_term_thread: {num_C_term_thread}')

            #print(f'N_term_thread: {N_term_thread}')
            #print(f'C_term_thread: {C_term_thread}')
            uent_df['N_term_thread'] += [num_N_term_thread]
            uent_df['C_term_thread'] += [num_C_term_thread]

            if num_N_term_thread != 0:
                min_N_thread_slippage_left = min(N_term_thread)
                min_N_thread_depth_left = min_N_thread_slippage_left / pdb_NCi_core
                min_N_prot_depth_left = min_N_thread_slippage_left / self.prot_size
            else:
                min_N_thread_slippage_left = np.nan
                min_N_thread_depth_left = np.nan
                min_N_prot_depth_left = np.nan
            uent_df['min_N_thread_slippage_left'] += [min_N_thread_slippage_left]
            uent_df['min_N_thread_depth_left'] += [min_N_thread_depth_left]
            uent_df['min_N_prot_depth_left'] += [min_N_prot_depth_left]
            
            if num_C_term_thread != 0:
                min_C_thread_slippage_right = self.prot_size - max(C_term_thread)
                denom = self.prot_size - pdb_NCj_core
                min_C_thread_depth_right = min_C_thread_slippage_right / denom if denom != 0 else np.nan
                min_C_prot_depth_right = min_C_thread_slippage_right / self.prot_size
            else:
                min_C_thread_slippage_right = np.nan
                min_C_thread_depth_right = np.nan
                min_C_prot_depth_right = np.nan
            uent_df['min_C_thread_slippage_right'] += [min_C_thread_slippage_right]
            uent_df['min_C_thread_depth_right'] += [min_C_thread_depth_right]
            uent_df['min_C_prot_depth_right'] += [min_C_prot_depth_right]
            #########################################################################
            

            #########################################################################
            ### Get entangled residues = NC_region U crossing_region
            #print('Get total entangled region residues')
            #print(f'NC_region: {NC_region}')
            #print(f'crossing_region_N: {crossing_region_N}')
            #print(f'crossing_region_C: {crossing_region_C}')
            ent_region = set(NC_region).union(set(crossing_region_N)).union(set(crossing_region_C))
            ent_region = ent_region.union(set(pdb_NC))
            ent_region = ent_region.union(set(pdb_crossing_res))

            #print(f'ent_region: {ent_region}')
            uent_df['ent_region'] += [",".join([str(r) for r in ent_region])]
     
            uent_df['ent_coverage'] += [len(ent_region)/self.prot_size]
            uent_df['prot_size'] += [self.prot_size]
            #########################################################################

        ### save file for unique entanglement features
        uent_df = pd.DataFrame(uent_df)
        #print(f'uent_df:\n{uent_df}')
        outfile = os.path.join(self.outdir, f'{gene}_{pdbid}_{chain}_uent_features.csv')
        uent_df.to_csv(outfile, index=False, sep='|')
        self.logger.info(f'Unique entanglement features saved to {outfile}')
        
        return {'outfile':outfile, 'results': uent_df}
    ########################################################################################################################
       
    ########################################################################################################################
    def find_neighboring_residues(self, traj, target_resids, cutoff=0.45):
        """
        Find all residues whose side-chain heavy atoms are within `cutoff` Å
        of any side-chain heavy atom of the residues in `target_resids`.

        Parameters
        ----------
        traj : md.Trajectory
            The trajectory (or single-frame PDB) to search.
        target_resids : list of int
            List of topology residue indices (residue.index) to probe around.
        cutoff : float, optional
            Distance cutoff in Å (default 4.5).

        Returns
        -------
        neighbors : list of int
            Sorted residue indices (residue.index) whose side-chain heavy atoms
            lie within `cutoff` of the target side-chain atoms.
            (Target residues themselves are excluded from the result.)
        """
        #print(f'Finding neighboring residues for target residues: {target_resids}')
        # --- 1. Identify side-chain heavy atoms in the target residues  ---
        query_atoms = [
            atom.index
            for atom in traj.topology.atoms
            if atom.residue.index in target_resids
            and atom.element.symbol != 'H'
            and atom.name not in ('N', 'CA', 'C', 'O')
        ]

        # --- 2. Identify side‑chain heavy atoms in all residues  ---
        haystack_atoms = [
            atom.index
            for atom in traj.topology.atoms
            if atom.element.symbol != 'H'
            and atom.name not in ('N', 'CA', 'C', 'O')
        ]
        #print(f'haystack_atoms: {haystack_atoms} {len(haystack_atoms)}')

        # --- 3. Use MDTraj’s neighbor search  ---
        #    For each frame, this returns indices *into* `haystack_atoms`
        neighbors_per_frame = md.compute_neighbors(
            traj, cutoff=cutoff,
            query_indices=query_atoms,
            haystack_indices=haystack_atoms
        )
        #print(f'neighbors_per_frame: {neighbors_per_frame}')
        

        # --- 4. Map back to global atom indices and then to residues  ---
        neighbor_atom_indices = set(neighbors_per_frame[0])
        #for frame_indices in neighbors_per_frame:
        #    print(f'frame_indices: {frame_indices}')
        #    for idx_in_haystack in frame_indices:
        #        neighbor_atom_indices.add(haystack_atoms[idx_in_haystack])

        neighbor_residue_indices = {
            traj.topology.atom(atom_idx).residue.resSeq
            for atom_idx in neighbor_atom_indices
        }
        #print(f'neighbor_residue_indices: {neighbor_residue_indices} {len(neighbor_residue_indices)}')

        # --- 5. Exclude the original target residues  ---
        neighbor_residue_indices -= set(target_resids)

        return sorted(neighbor_residue_indices)
    ########################################################################################################################
