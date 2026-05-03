###################################################################################################
## Multiprocesing helper function(s)
def process_frame(args):
    frame_idx, dcd, PSF, chain, chain_res, resids, atom_names, topoly, ID, Calpha, CG, g_threshold, density, ent_detection_method = args
    import MDAnalysis as mda
    from EntDetect.gaussian_entanglement import GaussianEntanglement
    univ = mda.Universe(PSF, dcd)
    chain_atoms = univ.select_atoms(f"segid {chain}")
    univ.trajectory[frame_idx]
    coor = chain_atoms.positions
    ge = GaussianEntanglement(g_threshold=g_threshold, density=density, Calpha=Calpha, CG=CG, ent_detection_method=ent_detection_method)
    ent_result = ge.get_traj_entanglements(coor, chain_res, resids, atom_names, topoly=topoly)
    result_rows = []
    if not ent_result:
        result_rows.append((ID, chain, frame_idx, '', '', '', '', '', '', '', '', '', '', ''))
    else:
        for ij_gN_gC, crossings in ent_result.items():
            # print(ij_gN_gC, crossings)
        
            i, j = ij_gN_gC[0], ij_gN_gC[1]
            crossingsN, crossingsC = crossings
            crossingsN = ','.join(crossingsN)
            crossingsC = ','.join(crossingsC)

            result_rows.append((ID, chain, frame_idx, i, j, crossingsN, crossingsC, f'{ij_gN_gC[2]: .5f}', f'{ij_gN_gC[3]: .5f}', ij_gN_gC[4], ij_gN_gC[5], ij_gN_gC[6], ij_gN_gC[7]))
    logging.getLogger('EntDetect.GaussianEntanglement').info(f'Frame idx: {frame_idx} processed with {len(result_rows)} contact(s) found.')
    return result_rows
###################################################################################################

###################################################################################################
## load in necessary packages
import logging
import itertools
import math
import sys
from Bio.PDB import PDBParser
import os
from operator import itemgetter
from warnings import filterwarnings
import numpy as np
import MDAnalysis as mda
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from topoly import lasso_type  # used pip
import re
import pickle
from collections import defaultdict
import time
import multiprocessing as mp
from EntDetect._logging import setup_logger
from typing import Optional
filterwarnings("ignore")

class GaussianEntanglement:
    """
    Gaussian Entanglement Class for calculating entanglements in protein structures derived from both experiments and AlphaFold. 
    """
    def __init__(self, g_threshold: float = 0.6, density: float = 0.0, Calpha: bool = False, CG: bool = False, nproc: int = 10, ent_detection_method: int = 2, log_level: int = logging.INFO, logdir: str = None) -> None:
        """
        Constructor for GaussianEntanglement class.

        Parameters
        ----------
        g_threshold : float, optional
            Threshold for Gaussian entanglement, by default 0.6
        density : float, optional
            Density for triangulation of minimal loop surface, by default 0.0
        Calpha : bool, optional
            Whether to use C-alpha atoms or heavy-atoms, by default False
        CG : bool, optional
            Whether the CG model was used to generate the simulations or structures
        ent_detection_method : int, optional
            Method to define ENT status of a raw NCLE:
            1 = any nonzero GLN for either termini
            2 = any nonzero TLN for either termini (default)
            3 = both GLN and TLN must have nonzero for same termini
        """

        self.logger = setup_logger('GaussianEntanglement', outdir=logdir, log_level=log_level)
        self.g_threshold = g_threshold
        self.density = density
        self.Calpha = Calpha
        self.CG = CG
        self.nproc = nproc
        self.ent_detection_method = ent_detection_method
        self.logger.info(f'GaussianEntanglement initialized with g_threshold: {g_threshold}, density: {density}, Calpha: {Calpha}, CG: {CG}, nproc: {nproc}, ent_detection_method: {ent_detection_method}')

        self.change_codes = {'L-C~': 'loss of linking number & switched linking chirality', 
                        'L-C#': 'loss of linking number & no change of linking chirality', 
                        'L+C~': 'gain of linking number & switched linking chirality', 
                        'L+C#': 'gain of linking number & no change of linking chirality', 
                        'L#C~': 'no change of linking number & switched linking chirality', 
                        'L#C#': 'no change'}
        # print class initialization parameters
        print(f'GaussianEntanglement initialized with g_threshold: {g_threshold}, density: {density}, Calpha: {Calpha}, CG: {CG}, nproc: {nproc}, ent_detection_method: {ent_detection_method}')
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def determine_ent_status(self, gln_n: float, gln_c: float, tln_n: int, tln_c: int) -> bool:
        """
        Determine if a native contact is entangled based on the selected detection method.
        
        Parameters
        ----------
        gln_n : float
            Gaussian linking number for N-terminus
        gln_c : float
            Gaussian linking number for C-terminus
        tln_n : int
            Topological linking number for N-terminus
        tln_c : int
            Topological linking number for C-terminus
        
        Returns
        -------
        bool
            True if entangled according to the selected method, False otherwise
        """
        if self.ent_detection_method == 1:
            # Any nonzero GLN for either termini
            return (gln_n != 0) or (gln_c != 0)
        elif self.ent_detection_method == 2:
            # Any nonzero TLN for either termini (default)
            return (tln_n != 0) or (tln_c != 0)
        elif self.ent_detection_method == 3:
            # Both GLN and TLN must have nonzero for same termini
            # Both N, both C, or both N and C
            n_both = (gln_n != 0) and (tln_n != 0)
            c_both = (gln_c != 0) and (tln_c != 0)
            return n_both or c_both
        else:
            raise ValueError(f"Invalid ent_detection_method: {self.ent_detection_method}. Must be 1, 2, or 3.")
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def helper_dot(self, Runit: np.ndarray, dR_cross: np.ndarray) -> list:

        """
        Numba function to speed up dot product calculation. Ability 
        to use current GPU (if available) and CPU
        
        """

        return [np.dot(x,y) for x,y in zip(Runit,dR_cross)]
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def point_rounding(self, num: float) -> float:

        """
        Rounds n based on the specified threshold:
        - If n is NaN, returns NaN.
        - For positive n: if fractional part >= threshold ? ceil; else ? floor.
        - For negative n: if |fractional part| >= threshold ? more negative; else ? toward zero.
        """
        # Handle NaN
        if isinstance(num, float) and math.isnan(num):
            return num

        int_part = math.trunc(num)      # truncate toward zero
        frac     = num - int_part
        if num >= 0:
            return int_part + (1 if frac >= self.g_threshold else 0)
        else:
            return int_part - (1 if abs(frac) >= self.g_threshold else 0)
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def get_entanglements(self, coor: np.ndarray, l: int, pdb_file: str, resids: np.ndarray, 
            resnames: np.ndarray,resid_index_to_ref_allatoms_idx: dict, ca_coor: np.ndarray, resid_index_to_resid: dict,
            termini_threshold: list=[5,5], loop_thread_thresh: list=[4,4], topoly: bool = True) -> dict:

        """
        Find proteins containing non-covalent lasso entanglements.

        Entanglements are composed of loops (defined by native contacts) and crossing residue(s).

        """
        Nterm_thresh = termini_threshold[0]
        Cterm_thresh = termini_threshold[1]
        loop_Nthread_thresh = loop_thread_thresh[0]
        loop_Cthread_thresh = loop_thread_thresh[1]
        # print(f'Finding entanglements with Nterm_thresh: {Nterm_thresh}, Cterm_thresh: {Cterm_thresh}, loop_Nthread_thresh: {loop_Nthread_thresh}, loop_Cthread_thresh: {loop_Cthread_thresh}')
        # print(f'Coordinates shape: {coor.shape}\n{coor[:10]}\n...\n{coor[-10:]}')

        # make native contact contact map
        dist_matrix = squareform(pdist(coor))
  
        if self.Calpha == False:
            native_cmap = np.where(dist_matrix <= 4.5, 1, 0) # if true then 1 will appear otherwise zero
        elif self.Calpha == True:
            native_cmap = np.where(dist_matrix <= 8, 1, 0) # if true then 1 will appear otherwise zero
        native_cmap = np.triu(native_cmap, k=4) # element below the 4th diagonal starting from middle are all zeros; # protein contact map

        num_res = len(resid_index_to_ref_allatoms_idx.keys())

        assert num_res == len(resids), f"Something's wrong with {pdb_file} residues {num_res} != {len(resids)}"

        res_ncmap = np.zeros((num_res, num_res))
        resid_pairs = list(itertools.product(np.arange(num_res), np.arange(num_res)))

        for pair in resid_pairs:
            pair0_resid = resid_index_to_resid[pair[0]]
            pair1_resid = resid_index_to_resid[pair[1]]
            # check that the resid are greater than 4 apart
            if abs(pair1_resid - pair0_resid) > 4:
                if pair[0] in resid_index_to_ref_allatoms_idx and pair[1] in resid_index_to_ref_allatoms_idx:
                    res1_atoms = resid_index_to_ref_allatoms_idx[pair[0]]
                    res2_atoms = resid_index_to_ref_allatoms_idx[pair[1]]
                    res1_atoms_start = min(res1_atoms)
                    res1_atoms_end = max(res1_atoms)
                    res2_atoms_start = min(res2_atoms)
                    res2_atoms_end = max(res2_atoms)
                    sub_array = native_cmap[res1_atoms_start:res1_atoms_end + 1, res2_atoms_start:res2_atoms_end + 1]
                    contact = np.sum(sub_array)
                    dist_sub_array = dist_matrix[res1_atoms_start:res1_atoms_end + 1, res2_atoms_start:res2_atoms_end + 1]
                    min_dist = np.min(dist_sub_array)

                    if contact > 0:
                        res_ncmap[pair[0], pair[1]] = 1
                        #print(f'Found contact: {pair0_resid} & {pair1_resid} with min dist: {min_dist}')        
                    # else:
                        # print(f'No contact: {pair0_resid} & {pair1_resid} with min dist: {min_dist}')
        del native_cmap
        native_cmap = res_ncmap 

        nc_indexs = np.stack(np.nonzero(native_cmap)).T # stack indices based on rows

        # make R coordinate and gradient of R length N-1
        range_l = np.arange(0, l-1)
        range_next_l = np.arange(1,l)

        ca_coor = ca_coor.astype(np.float32)
        R = 0.5*(ca_coor[range_l] + ca_coor[range_next_l])
        dR = ca_coor[range_next_l] - ca_coor[range_l]

        #make dRcross matrix
        pair_array = np.asarray(list(itertools.product(dR,dR))) # combination of elements within array

        x = pair_array[:,0,:]
        y = pair_array[:,1,:]

        dR_cross = np.cross(x, y)

        #make Rnorm matrix
        pair_array = np.asarray(list(itertools.product(R,R)))
        diff = pair_array[:,0,:] - pair_array[:,1,:]
        diff = diff.astype(np.float32)

        Runit = diff / np.linalg.norm(diff, axis=1)[:,None]**3 
        Runit = Runit.astype(np.float32)

        #make final dot matrix
        dot_matrix = self.helper_dot(Runit, dR_cross)
        dot_matrix = np.asarray(dot_matrix)
        dot_matrix = dot_matrix.reshape((l-1,l-1))

        nc_gdict = {} 

        for i,j in nc_indexs:

            # loop_range = np.arange(i,j)
            # nterm_range = np.arange(Nterm_thresh,i-5)
            # cterm_range = np.arange(j+6,l-(Cterm_thresh + 1))
            loop_range = np.arange(i, j)
            nterm_range = np.arange(Nterm_thresh, i-loop_Nthread_thresh-1)
            cterm_range = np.arange(j+loop_Cthread_thresh+1, l-(Cterm_thresh + 1))

            gn_pairs_array = np.fromiter(itertools.chain(*itertools.product(nterm_range, loop_range)), int).reshape(-1, 2)
            gc_pairs_array = np.fromiter(itertools.chain(*itertools.product(loop_range, cterm_range)), int).reshape(-1, 2)

            if gn_pairs_array.size != 0:
                
                gn_vals = dot_matrix[gn_pairs_array[:,0],gn_pairs_array[:,1]]
                gn_vals = gn_vals[~np.isnan(gn_vals)] 
                gn_val = np.sum(gn_vals) / (4.0 * np.pi)
            
            else:
                gn_val = 0

            if gc_pairs_array.size != 0:
                
                gc_vals = dot_matrix[gc_pairs_array[:,0],gc_pairs_array[:,1]]
                gc_vals = gc_vals[~np.isnan(gc_vals)] 
                gc_val = np.sum(gc_vals) / (4.0 * np.pi)
            
            else:
                gc_val = 0
            
            rounded_gc_val = self.point_rounding(np.float64(gc_val))
            rounded_gn_val = self.point_rounding(np.float64(gn_val))            

            #if np.abs(rounded_gn_val) >= 1 or np.abs(rounded_gc_val) >= 1:
            #    #print(f'({i}, {j}) with gn: {gn_val} and gc: {gc_val}')
            #    nc_gdict[ (int(i), int(j)) ] = (gn_val, gc_val, rounded_gn_val, rounded_gc_val)
            nc_gdict[ (int(i), int(j)) ] = (gn_val, gc_val, rounded_gn_val, rounded_gc_val)

        missing_residues = self.find_missing_residues(resids)
        #print(f'missing_residues: {missing_residues}')

        filtered_nc_gdict = self.loop_filter(nc_gdict, resids, missing_residues)
        #print(f'size filtered_nc_gdict after accounting for missing residues: {len(filtered_nc_gdict)}\n{filtered_nc_gdict}')

        if topoly:
            entangled_res = self.find_crossing(ca_coor.tolist(), filtered_nc_gdict, resids)
            #print(f'entangled_res: {entangled_res}')

            filtered_entangled_res = self.crossing_filter(entangled_res, missing_residues)
            #print(f'filtered_entangled_res: {filtered_entangled_res}')
        else:
            filtered_entangled_res = {}
            for ij, values in filtered_nc_gdict.items():
                i, j = ij[0], ij[1]
                gn = values[0]
                gc = values[1]
                GLNn = values[2]
                GLNc = values[3]
                TLNn = np.nan
                TLNc = np.nan
                filtered_entangled_res[(resids[i], resids[j], gn, gc, GLNn, GLNc, TLNn, TLNc)] = [np.array([]), np.array([])]

        return filtered_entangled_res, missing_residues
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def find_missing_residues(self, resids:np.ndarray) -> np.ndarray:

        """
        Find missing residues in pdb file

        """

        check_all_resids = np.arange(resids[0], resids[-1] + 1)

        missing_residues = np.setdiff1d(check_all_resids, resids)

        return missing_residues
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def loop_filter(self, native_contacts: dict, resids: np.ndarray, missing_res: np.ndarray) -> dict:

        """
        Remove loops if there are three or more consecutive missing residues
        or the amount of any missing residues exceed 5% of the loop length 

        """

        for ij, values in native_contacts.items():

            native_i = resids[ij[0]]

            native_j = resids[ij[1]]

            rounded_gn = values[-2]

            rounded_gc = values[-1]

            check_loop = np.arange(native_i , native_j + 1) 

            loop_length = check_loop.size

            missing_res_loop = np.intersect1d(check_loop, missing_res)

            for index, diff_resid_index in itertools.groupby(enumerate(missing_res_loop), lambda ix : ix[0] - ix[1]):

                conseuctive_missing_residues = list(map(itemgetter(1), diff_resid_index))

                if len(conseuctive_missing_residues) >= 3 or len(missing_res_loop) > 0.05 * loop_length:

                    native_contacts[ij] = None

        return native_contacts
    ##########################################################################################################################################################


    ##########################################################################################################################################################
    def find_crossing(self, coor: np.ndarray, nc_data: dict, resids: np.ndarray, terminal_buff:int=5, loop_buff:int=4) -> dict:

        """
        Use Topoly to find crossing(s) based on partial linking number

        """

        entangled_res = {}

        native_contacts = [[ij[0], ij[1]] for ij, values in nc_data.items() if values is not None]

        data = lasso_type(coor, loop_indices=native_contacts, more_info=True, density=self.density, min_dist=[10, loop_buff, terminal_buff])
        # data = lasso_type(coor, loop_indices=native_contacts, more_info=True, density=self.density, min_dist=[0, 4, 5])
        for native_contact in native_contacts:
            # print(f'native_contact:\n{native_contact}')

            crossings = []

            native_contact = tuple(native_contact)
            i, j = resids[native_contact[0]], resids[native_contact[1]]
            # print('\n', (i,j), native_contact, data[native_contact])
            ## Parse the N terminal crossings
            # crossingN = [f"{cr[0]}{resids[int(cr[1:])]}" for cr in data[native_contact]["crossingsN"]]
            crossingN = []
            beforeN = [f"{cr[0]}{resids[int(cr[1:])]}" for cr in data[native_contact]["beforeN"] if '*' not in cr]
            if beforeN:
                # remove residues that violate the terminial and loop buffers
                Ncrossing_resids = sorted([int(cr) for cr in beforeN], key=lambda x: abs(x), reverse=True)
                # print(f'Ncrossing_resids before filtering:', Ncrossing_resids)
                filtered_crossingN_resids = []
                for c in Ncrossing_resids:
                    if abs(c) < (i - terminal_buff) and abs(c) > terminal_buff:
                        filtered_crossingN_resids.append(c)
                # print(f'Filtered crossings after checking termini and loop buffers: {filtered_crossingN_resids}')

                if len(filtered_crossingN_resids) > 1:
                    # check the distance between each pair of residues. If a pair is less than 10 residues apart, remove both crossings
                    # Use a greedy approach: iterate and keep crossings that are >=10 apart from the previous kept crossing
                    results = [filtered_crossingN_resids[0]]
                    for k in range(1, len(filtered_crossingN_resids)):
                        # Check if current crossing is >= 10 away from the last kept crossing
                        if abs(abs(filtered_crossingN_resids[k]) - abs(results[-1])) >= 10:
                            results.append(filtered_crossingN_resids[k])
                    filtered_crossingN_resids = results

                crossingN = []
                for c in filtered_crossingN_resids:
                    if c > 0:
                        crossingN.append(f'+{c}')
                    else:
                        crossingN.append(f'{c}') 
            # print(f'Final crossingN: {crossingN}\n')
            crossings += crossingN


            ## Parse the C terminal crossings
            # crossingC = [f"{cr[0]}{resids[int(cr[1:])]}" for cr in data[native_contact]["crossingsC"]]
            crossingC = []
            # print("Before C:", data[native_contact]["beforeC"])
            beforeC = [f"{cr[0]}{resids[int(cr[1:])]}" for cr in data[native_contact]["beforeC"] if '*' not in cr]
            # print("Before C (formatted):", beforeC)
            if beforeC:
                # remove residues that violate the terminial and loop buffers
                Ccrossing_resids = sorted([int(cr) for cr in beforeC], key=lambda x: abs(x))
                # print(f'Ccrossing_resids before filtering:', Ccrossing_resids)
                filtered_crossingC_resids = []
                for c in Ccrossing_resids:
                    if abs(c) > (j + terminal_buff) and abs(c) < (resids[-1] - terminal_buff):
                        filtered_crossingC_resids.append(c)

                # print(f'Filtered crossings after checking termini and loop buffers: {filtered_crossingC_resids}')
                if len(filtered_crossingC_resids) > 1:
                    # check the distance between each pair of residues. If a pair is less than 10 residues apart, remove both crossings
                    # Use a greedy approach: iterate and keep crossings that are >=10 apart from the previous kept crossing
                    results = [filtered_crossingC_resids[0]]
                    for k in range(1, len(filtered_crossingC_resids)):
                        # Check if current crossing is >= 10 away from the last kept crossing
                        if abs(abs(filtered_crossingC_resids[k]) - abs(results[-1])) >= 10:
                            results.append(filtered_crossingC_resids[k])
                    filtered_crossingC_resids = results

                crossingC = []
                for c in filtered_crossingC_resids:
                    if c > 0:
                        crossingC.append(f'+{c}')
                    else:
                        crossingC.append(f'{c}')

            # print(f'Final crossingC: {crossingC}\n')
            crossings += crossingC

            gn = nc_data[native_contact][0]
            GLNn = nc_data[native_contact][2]
            
            gc = nc_data[native_contact][1]
            GLNc = nc_data[native_contact][3]

            if crossingN: 
                TLNn_signs = [c[0] for c in crossingN if int(c[1:]) < i]
                TLNn = [1 for s in TLNn_signs if s == '+'] + [-1 for s in TLNn_signs if s == '-']
                TLNn = sum(TLNn)
            else: 
                TLNn = 0

            if crossingC: 
                TLNc_signs = [c[0] for c in crossingC if int(c[1:]) > j]
                TLNc = [1 for s in TLNc_signs if s == '+'] + [-1 for s in TLNc_signs if s == '-']
                TLNc = sum(TLNc)
            else:
                TLNc = 0

            ij_gN_gC = (resids[native_contact[0]], resids[native_contact[1]]) + (gn, gc) + (GLNn, GLNc) + (TLNn, TLNc)

            entangled_res[ij_gN_gC] = [np.unique(crossingN), np.unique(crossingC)]
        
        return entangled_res
    ##########################################################################################################################################################


    ##########################################################################################################################################################
    def crossing_filter(self, entanglements: dict, missing_res: np.ndarray) -> dict:

        """
        Remove entanglements if there are any missing residues plus-and-minus 10 of the crossing(s)

        """
        filtered_entanglements = {}
        for ij_gN_gC, crossings in entanglements.items():
            crossingsN, crossingsC = crossings
            
            # Convert to list for easier manipulation
            crossingsN_filtered = list(crossingsN)
            crossingsC_filtered = list(crossingsC)

            if len(crossingsN_filtered) > 0:
                crossings_to_remove = []
                for crossing in crossingsN_filtered:
                    reg_exp = re.split("\\+|-", crossing, maxsplit=1) # split the chirality
                    check_crossing = np.arange(int(reg_exp[1]) - 10 , int(reg_exp[1]) + 11)
                    missing_res_cr = np.intersect1d(check_crossing, missing_res)

                    if missing_res_cr.size:
                        crossings_to_remove.append(crossing)
                
                for crossing in crossings_to_remove:
                    crossingsN_filtered.remove(crossing)

            if len(crossingsC_filtered) > 0:
                crossings_to_remove = []
                for crossing in crossingsC_filtered:
                    reg_exp = re.split("\\+|-", crossing, maxsplit=1) # split the chirality
                    check_crossing = np.arange(int(reg_exp[1]) - 10 , int(reg_exp[1]) + 11)
                    missing_res_cr = np.intersect1d(check_crossing, missing_res)

                    if missing_res_cr.size:
                        crossings_to_remove.append(crossing)
                
                for crossing in crossings_to_remove:
                    crossingsC_filtered.remove(crossing)

            filtered_entanglements[ij_gN_gC] = [np.array(crossingsN_filtered), np.array(crossingsC_filtered)]

        # filtered_entanglements = {nc: re_cr for nc, re_cr in entanglements.items() if re_cr is not None and len(re_cr) > 0} 
        # filtered_entanglements = {nc: re_cr for nc, re_cr in entanglements.items() if re_cr is not None} 

        return filtered_entanglements
    ##########################################################################################################################################################


    ##########################################################################################################################################################
    def check_disulfideBonds(self, pdb_file): 

        # Parse the PDB structure
        parser = PDBParser()
        structure = parser.get_structure('protein', pdb_file)

        disulfide_bonds = []
        # Iterate over residues and identify disulfide bonds
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_resname() == 'CYS':
                        if 'SG' in residue:
                            sg_atom = residue['SG']
                        else:
                            continue
                        # Check for disulfide bonds with distance threshold (e.g., <2.2 Å)
                        for model2 in structure:
                            for chain2 in model2:
                                for residue2 in chain2:
                                    if residue2.get_resname() == 'CYS' and residue != residue2:
                                        if 'SG' in residue2:
                                            sg_atom2 = residue2['SG']
                                        else:
                                            continue
                                        distance = sg_atom - sg_atom2
                                        if distance < 2.2:
                                            self.logger.info(f"Disulfide bond between {residue} and {residue2} at distance {distance:.2f} Å")
                                            i,j = residue.get_id()[1], residue2.get_id()[1] 
                                            
                                            if (i,j) not in disulfide_bonds and (j,i) not in disulfide_bonds:
                                                disulfide_bonds += [(i,j)]
        return disulfide_bonds 
    ##########################################################################################################################################################


    ##########################################################################################################################################################
    def calculate_native_entanglements(self, pdb_file: str, outdir: str, ID: str='', chain: str=None, topoly: bool = True) -> dict:

        """
        Driver function that outputs native lasso-like self entanglements and missing residues for pdb and all of its chains if any
        """

        ## set up the outdir for this calculation
        #outdir = f"{os.getcwd()}/{outdir}"
        if not os.path.isdir(outdir):
            os.mkdir(f"{outdir}") 
            self.logger.info(f"Creating directory: {outdir}")

        if not os.path.isdir(f"{outdir}/unmapped_missing_residues"):
            os.mkdir(f"{outdir}/unmapped_missing_residues")
            self.logger.info(f"Creating directory: {outdir}/unmapped_missing_residues")


        ## get the PDB file name
        pdb = pdb_file.split('/')[-1].split(".")[0]
        if ID == '':
            ID = pdb
        self.logger.info(f"\n{'#'*100}\nCOMPUTING ENTANGLEMENTS FOR \033[4m{pdb}\033[0m with ID {ID}")

        ## Define the outfile and check if it exists. If so load it else create it
        outfile = os.path.join(f'{outdir}', f'{ID}_GE.csv')
        if os.path.exists(outfile):
            self.logger.info(f'{outfile} ALREADY EXISTS AND WILL BE LOADED')
            outdf = pd.read_csv(outfile, sep='|', dtype={'c': str})
            return {'outfile':outfile, 'ent_result':outdf}

        ## Load the reference universe and use the MDA parser for it with the special excapetion for our CG files that end in .cor
        if pdb_file.endswith('.cor'):
            ref_univ = mda.Universe(f'{pdb_file}', format='CRD')
        else:
            ref_univ = mda.Universe(f"{pdb_file}")  
        self.logger.debug(f'ref_univ: {ref_univ}')


        ### Find any disulfide bonds
        self.logger.info(f'Checking for disulfide bonds')
        disulfide_bonds = self.check_disulfideBonds(pdb_file)
        self.disulfide_bonds = disulfide_bonds
        self.logger.debug(f'disulfide_bonds: {disulfide_bonds}')
        

        ## Get only the heavy atoms or CA atoms depending on what type of contact we are looking for
        if self.Calpha == False and self.CG == False:
            self.logger.info('All-atom model and contacts: Selecting all heavy atoms (no hydrogens) for entanglement calculations')
            ref_allatoms_dups = ref_univ.select_atoms("not name H* and protein")
        elif self.Calpha == True and self.CG == False:
            self.logger.info('All-atom model and Calpha contacts: Selecting only CA atoms for entanglement calculations')
            ref_allatoms_dups = ref_univ.select_atoms("name CA and protein")
        elif self.CG == True:
            self.logger.info('Coarse-grained model: Selecting all atoms for entanglement calculations')
            ref_allatoms_dups = ref_univ.select_atoms("all")
        #print(f'ref_allatoms_dups: {ref_allatoms_dups} {len(ref_allatoms_dups)}')

        chains_to_analyze = set(ref_univ.segments.segids)
        if chain is not None:
            chains_to_analyze = {chain} if chain in chains_to_analyze else set()
            if not chains_to_analyze:
                raise ValueError(f"Chain {chain} not found in structure. Available chains: {set(ref_univ.segments.segids)}")

        for chain in chains_to_analyze:

            ## Check for duplicate residues
            atom_data = []
            check = set()

            for atom in ref_allatoms_dups.select_atoms(f"segid {chain}").atoms:

                atom_data.append((atom.resid, atom.name))
                check.add((atom.resid, atom.name))

            temp_df = pd.DataFrame(atom_data, columns=["resid", "name"])

            unique_rows = temp_df.drop_duplicates()
            unique_indices = unique_rows.index.tolist()

            assert len(check) == len(unique_indices), "You did not remove dup atoms!"

            ref_allatoms_unique = ref_allatoms_dups.select_atoms(f"segid {chain}")[unique_indices]
            #print(f'ref_allatoms_unique: {ref_allatoms_unique} {len(ref_allatoms_unique)}')

            ## select only those unique residue alpha carbons
            if self.CG == False:
                ref_ca_unique = ref_allatoms_unique.select_atoms("name CA")
            else:
                ref_ca_unique = ref_allatoms_unique.select_atoms("all")
            #print(f'ref_ca_unique: {ref_ca_unique} {len(ref_ca_unique)}')
  

            resid_index_to_ref_allatoms_idx = {}
            resid_index_to_resid = {}
            ref_allatoms_idx_to_resid = {}
            atom_ix = 0
            res_ix = 0
            PDB_resids = ref_ca_unique.resids
            #print(f'PDB_resids: {PDB_resids}')
            new_atm_idx = []

            ## QC if the chain has only one alpha carbon or none
            if len(PDB_resids) == 0 or len(PDB_resids) == 1:
                raise ValueError(f"Skipping over chain {chain} for \033[4m{pdb}\033[0m since chain has only one alpha carbon or none")


            for atom in ref_allatoms_unique.atoms:
                new_atm_idx.append(atom_ix)
                ref_allatoms_idx_to_resid[atom_ix] = [atom.resid]

                if atom_ix == new_atm_idx[0]:
                    resid = atom.resid
                    resid_index_to_ref_allatoms_idx[res_ix] = [atom_ix]
                    resid_index_to_resid[res_ix] = resid
                    atom_ix += 1
                else:
                    if atom.resid == resid:
                        resid_index_to_ref_allatoms_idx[res_ix] += [atom_ix]
                        resid_index_to_resid[res_ix] = resid
                        resid = atom.resid
                        atom_ix += 1
                    else:
                        res_ix += 1
                        resid_index_to_ref_allatoms_idx[res_ix] = [atom_ix]
                        resid_index_to_resid[res_ix] = resid
                        atom_ix += 1
                        resid = atom.resid
            
            ## Quality check that if Calpha is True there is 1-to-1 mapping of resid index to allatom indexs
            if self.Calpha == True:
                for k,v in resid_index_to_ref_allatoms_idx.items():
                    if len(v) != 1:
                        raise ValueError(f'When Calpha is specified there should only be one resid index for each all atom index: resid index {k} has {v}')

            assert len(new_atm_idx) == np.concatenate(list(resid_index_to_ref_allatoms_idx.values())).size, f"Not enough atom indicies! {pdb_file}"
            
            # x y z cooridnates of chain
            coor = ref_allatoms_unique.atoms.positions[new_atm_idx]

            for resid_idx, all_atom_idx in resid_index_to_ref_allatoms_idx.items():

                resid = PDB_resids[resid_idx]
                #print(resid_idx, resid, all_atom_idx)

                check_coor = coor[all_atom_idx]
                #print(f'check_coor:\n{check_coor}')

                structure_coor = ref_allatoms_unique.select_atoms(f"resid {resid}").positions
                #print(f'structure_coor:\n{structure_coor}')
                
                try:
                    np.all(check_coor == structure_coor)
                except:
                    raise ValueError(f'Error in checking residue coordinates: most likely caused by resides with letters after them or specifying CG=True when it is infact an allatom model. Check resid: {resid} in PDB')

                if not np.all(check_coor == structure_coor):
                    raise ValueError(f"Coordinates do not match up! Resid {resid} {pdb_file}")


            ca_coor = ref_ca_unique.positions

            resnames = ref_ca_unique.resnames

            chain_res = PDB_resids.size

            if PDB_resids.size:

                ent_result, missing_residues = self.get_entanglements(coor, chain_res, pdb_file, PDB_resids, resnames, resid_index_to_ref_allatoms_idx, ca_coor, resid_index_to_resid, topoly=topoly)
                # print(f'Number of ENTs found: {len(ent_result)}\n{ent_result}')
      
                
                ## If there is Native entanglement then save a file
                if ent_result: 
                    outfile = os.path.join(f'{outdir}', f'{ID}_GE.csv')
                    #print(f'WRITING: {outfile}')
                    outdf = {'ID':[], 'chain':[], 'i':[], 'j': [], 'crossingsN': [], 'crossingsC': [], 'gn':[], 'gc':[], 'GLNn':[], 'GLNc':[], 'TLNn':[], 'TLNc':[], 'CCbond':[]}

                    for ij_gN_gC, crossings in ent_result.items():
                        # print(f'ij_gN_gC: {ij_gN_gC} with crossings: {crossings}')
                        crossingsN, crossingsC = crossings
                        crossingsN = ','.join(crossingsN)
                        crossingsC = ','.join(crossingsC)
                        i, j, gn, gc, GLNn, GLNc, TLNn, TLNc = ij_gN_gC

                        if (i,j) in disulfide_bonds or (j,i) in disulfide_bonds:
                            CCbond = True
                        else:
                            CCbond = False

                        # print(f'Contact: ({i}, {j}) with GLNn: {GLNn} and GLNc: {GLNc} has crossings: {crossingsN} {crossingsC} and CCbond: {CCbond}')

                        outdf['ID'] += [ID]
                        outdf['chain'] += [chain]
                        outdf['i'] += [i]
                        outdf['j'] += [j]
                        outdf['crossingsN'] += [crossingsN]
                        outdf['crossingsC'] += [crossingsC]
                        outdf['gn'] += [f'{gn: .5f}']
                        outdf['gc'] += [f'{gc: .5f}']
                        outdf['GLNn'] += [GLNn]
                        outdf['GLNc'] += [GLNc]
                        outdf['TLNn'] += [TLNn]
                        outdf['TLNc'] += [TLNc]
                        outdf['CCbond'] += [CCbond]
                    
                    outdf = pd.DataFrame(outdf)
                    outdf['ENT'] = outdf.apply(lambda row: self.determine_ent_status(row['GLNn'], row['GLNc'], row['TLNn'], row['TLNc']), axis=1)
                    outdf.to_csv(outfile, sep='|', index=False)
                    self.logger.info(f'SAVED: {outfile}')
                    ent_result = pd.read_csv(outfile, sep='|', dtype={'crossingsN': str, 'crossingsC': str})
                else:
                    self.logger.info(f'NO CONTACTS DETECTED for {pdb}')
                    ent_result = pd.DataFrame({'ID':[], 'chain':[], 'i':[], 'j': [], 'crossingsN': [], 'crossingsC': [], 'gn':[], 'gc':[], 'GLNn':[], 'GLNc':[], 'TLNn':[], 'TLNc':[], 'CCbond':[], 'ENT':[]})
                    ent_result.to_csv(outfile, sep='|', index=False)
                    self.logger.info(f'SAVED: {outfile}')
                    ent_result = pd.read_csv(outfile, sep='|', dtype={'c': str})
                                
                if len(missing_residues):
                    self.logger.info(f'WRITING: {pdb}_M.txt')
                    with open(f"{outdir}/unmapped_missing_residues/{pdb}_M.txt", "a") as f:
                        f.write(f"Chain {chain}: ")
                        for m_res in missing_residues:
                            f.write(f"{m_res} ")
                        f.write("\n")
    
        ## Return a dictionary with the outfile and the results
        return {'outfile':outfile, 'ent_result':ent_result}
            

    ##########################################################################################################################################################
    def calculate_traj_entanglements(
        self,
        dcd: str,
        PSF: str,
        outdir: str = './',
        ID: str = '',
        topoly: bool = True,
        start: int = 0,
        stop: int = 999999999,
        stride: int = 1,
        ref_contact_file: Optional[str] = None,
    ) -> dict:

        """
        Driver function that takes a CA coarse grained MD trajectory and looks for entanglements. 
        """

        ## set up the outdir for this calculation
        #outdir = f"{os.getcwd()}/{outdir}"
        if not os.path.isdir(outdir):
            os.mkdir(f"{outdir}") 
            self.logger.info(f"Creating directory: {outdir}")

        ## get the DCD file name
        dcd_name = dcd.split('/')[-1].split(".")[0]
        if ID == '':
            ID = dcd_name
        self.logger.info(f"\n{'#'*100}\nCOMPUTING ENTANGLEMENTS FOR \033[4m{dcd_name}\033[0m with ID {ID}")
        self.logger.debug(f'Topoly: {topoly} {type(topoly)}')
    

        ## Define the outfile and check if it exists. If so load it else create it
        outfile = os.path.join(f'{outdir}', f'{ID}_GE.csv')
        if os.path.exists(outfile):
            self.logger.info(f'{outfile} ALREADY EXISTS AND WILL BE LOADED')
            outdf = pd.read_csv(outfile, sep='|', dtype={'c': str})
            if ref_contact_file is not None:
                try:
                    ref_df = pd.read_csv(ref_contact_file, sep='|', usecols=['chain', 'i', 'j'])
                    ref_min = np.minimum(ref_df['i'].astype(int).to_numpy(), ref_df['j'].astype(int).to_numpy())
                    ref_max = np.maximum(ref_df['i'].astype(int).to_numpy(), ref_df['j'].astype(int).to_numpy())
                    ref_keys = (
                        ref_df['chain'].astype(str).to_numpy().astype(str)
                        + ':'
                        + ref_min.astype(str)
                        + '-'
                        + ref_max.astype(str)
                    )
                    ref_keys = set(ref_keys.tolist())

                    out_min = np.minimum(outdf['i'].astype(int).to_numpy(), outdf['j'].astype(int).to_numpy())
                    out_max = np.maximum(outdf['i'].astype(int).to_numpy(), outdf['j'].astype(int).to_numpy())
                    out_keys = (
                        outdf['chain'].astype(str).to_numpy().astype(str)
                        + ':'
                        + out_min.astype(str)
                        + '-'
                        + out_max.astype(str)
                    )
                    mask = pd.Series(out_keys).isin(ref_keys).to_numpy()
                    filtered = outdf.loc[mask].reset_index(drop=True)
                    if len(filtered) != len(outdf):
                        self.logger.info(
                            f'Filtered existing Traj_GE output to reference contacts: '
                            f'{len(outdf)} -> {len(filtered)} rows (ref: {ref_contact_file})'
                        )
                        filtered.to_csv(outfile, sep='|', index=False)
                        outdf = filtered
                except Exception as e:
                    self.logger.warning(f'WARNING: Failed to filter Traj_GE output using ref_contact_file={ref_contact_file}: {e}')

            return {'outfile': outfile, 'ent_result': outdf}
    
        ## Else analyze the traj and create the outfile
        univ = mda.Universe(PSF, dcd)    
        self.logger.debug(f'univ: {univ}')

        chains_to_analyze = set(univ.segments.segids)

        ## define the output dataframe
        # outdf = {'ID': [], 'chain':[], 'frame':[], 'i':[], 'j': [], 'c': [], 'gn':[], 'gc':[], 'Gn':[], 'Gc':[]}
        rows = []
        
        for chain in chains_to_analyze:
            self.logger.info(f'Analyzing chain {chain}')

            # Get the coordinates of the chain
            chain_atoms = univ.select_atoms(f"segid {chain}")

            resids = chain_atoms.resids
            self.logger.debug(f'resids: {resids}')

            resnames = chain_atoms.resnames
            self.logger.debug(f'resnames: {resnames}')

            chain_res = resids.size
            self.logger.debug(f'chain_res: {chain_res}')

            atom_names = chain_atoms.names
            self.logger.debug(f'atom_names: {atom_names}')


            frame_indices = []
            for ts in univ.trajectory[start:stop:stride]:
                frame_indices.append(ts.frame)
            self.logger.info(f'Analyzing frames from {start} to {stop} with stride {stride}')
            self.logger.info(f'Total frames to analyze: {len(frame_indices)}')
            self.logger.debug(f'frame_indices: {frame_indices[:10]} ... {frame_indices[-10:]}')

            pool_args = [
                (frame_idx, dcd, PSF, chain, chain_res, resids, atom_names, topoly, ID, self.Calpha, self.CG, self.g_threshold, self.density, self.ent_detection_method)
                for frame_idx in frame_indices
            ]
            import multiprocessing as mp
            with mp.get_context("spawn").Pool(processes=self.nproc) as pool:
                all_rows = pool.map(process_frame, pool_args)

            for frame_rows in all_rows:
                rows.extend(frame_rows)

        # outdf = pd.DataFrame(outdf)
        outdf = pd.DataFrame(rows, columns=['ID','chain','frame','i','j','crossingsN','crossingsC','gn','gc','GLNn','GLNc','TLNn','TLNc'])
        outdf['frame'] = pd.to_numeric(outdf['frame'], errors='coerce')
        outdf = outdf.sort_values(by='frame', ascending=True).reset_index(drop=True)

        # ENT = []
        # for idx, row in outdf.iterrows():
        #     if (row['TLNn'] != 0) or (row['TLNc'] != 0):
        #         ENT.append(True)
        #     else:
        #         ENT.append(False)
        # outdf['ENT'] = ENT
        # print(f'outdf:\n{outdf.to_string()}')
        outdf['ENT'] = outdf.apply(lambda row: self.determine_ent_status(row['GLNn'], row['GLNc'], row['TLNn'], row['TLNc']), axis=1)

        if ref_contact_file is not None:
            try:
                ref_df = pd.read_csv(ref_contact_file, sep='|', usecols=['chain', 'i', 'j'])
                ref_min = np.minimum(ref_df['i'].astype(int).to_numpy(), ref_df['j'].astype(int).to_numpy())
                ref_max = np.maximum(ref_df['i'].astype(int).to_numpy(), ref_df['j'].astype(int).to_numpy())
                ref_keys = (
                    ref_df['chain'].astype(str).to_numpy().astype(str)
                    + ':'
                    + ref_min.astype(str)
                    + '-'
                    + ref_max.astype(str)
                )
                ref_keys = set(ref_keys.tolist())

                out_min = np.minimum(outdf['i'].astype(int).to_numpy(), outdf['j'].astype(int).to_numpy())
                out_max = np.maximum(outdf['i'].astype(int).to_numpy(), outdf['j'].astype(int).to_numpy())
                out_keys = (
                    outdf['chain'].astype(str).to_numpy().astype(str)
                    + ':'
                    + out_min.astype(str)
                    + '-'
                    + out_max.astype(str)
                )
                mask = pd.Series(out_keys).isin(ref_keys).to_numpy()
                before = len(outdf)
                outdf = outdf.loc[mask].reset_index(drop=True)
                after = len(outdf)
                self.logger.info(f'Filtered Traj_GE to reference contacts: {before} -> {after} rows (ref: {ref_contact_file})')
            except Exception as e:
                self.logger.warning(f'WARNING: Failed to filter Traj_GE output using ref_contact_file={ref_contact_file}: {e}')

        self.logger.info(f'outdf:\n{outdf.to_string()}')
        outdf.to_csv(outfile, sep='|', index=False)
        self.logger.info(f'SAVED: {outfile}')

        # Multiprocessing disables per-frame timing, so skip frame_times and mean_time reporting
        # Return a dictionary with the outfile and the results
        return {'outfile':outfile, 'ent_result':outdf}
    ##########################################################################################################################################################

    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def get_traj_entanglements(self, coor: np.ndarray, l: int, resids: np.ndarray, atom_names: np.ndarray, topoly:bool=True, dist_cutoff:float=8.0, termini_threshold: list=[5,5], loop_thread_thresh: list=[4,4]) -> dict:

        """
        Find proteins containing non-covalent lasso entanglements.

        Entanglements are composed of loops (defined by native contacts) and crossing residue(s).

        """

        ## Check that if topoly is False then ent_detection_method is not 2 or 3 since those require TLN which requires topoly
        if not topoly and self.ent_detection_method in (2, 3):
            self.logger.warning(f'topoly=False but ent_detection_method={self.ent_detection_method} requires TLN — no NCLEs will be detected. Use ent_detection_method=1 (GLN) when topoly is disabled.')
            raise ValueError(f'topoly=False but ent_detection_method={self.ent_detection_method} requires TLN — no NCLEs will be detected. Use ent_detection_method=1 (GLN) when topoly is disabled.')

        Nterm_thresh = termini_threshold[0]
        Cterm_thresh = termini_threshold[1]
        loop_Nthread_thresh = loop_thread_thresh[0]
        loop_Cthread_thresh = loop_thread_thresh[1]
        self.logger.info(f'Finding entanglements with Nterm_thresh: {Nterm_thresh}, Cterm_thresh: {Cterm_thresh}, loop_Nthread_thresh: {loop_Nthread_thresh}, loop_Cthread_thresh: {loop_Cthread_thresh}')

        # make native contact contact map
        native_cmap, bb_coor, dist_matrix = self.processes_coor(coor, resids, atom_names, CG=self.CG, Calpha=self.Calpha)
        l = len(bb_coor)
        # print(f'native_cmap: {native_cmap.shape} {native_cmap}')
        # print(f'bb_coor: {bb_coor.shape} {bb_coor}')
        # print(f'dist_matrix: {dist_matrix.shape} {dist_matrix}')

        nc_indexs = np.stack(np.nonzero(native_cmap)).T # stack indices based on rows
        # print(f'nc_indexs: {nc_indexs.shape} {nc_indexs}')

        # make R coordinate and gradient of R length N-1
        range_l = np.arange(0, l-1)
        range_next_l = np.arange(1,l)

        bb_coor = bb_coor.astype(np.float32)
        R = 0.5*(bb_coor[range_l] + bb_coor[range_next_l])
        dR = bb_coor[range_next_l] - bb_coor[range_l]

        #make dRcross matrix
        pair_array = np.asarray(list(itertools.product(dR,dR))) # combination of elements within array

        x = pair_array[:,0,:]
        y = pair_array[:,1,:]

        dR_cross = np.cross(x, y)

        #make Rnorm matrix
        pair_array = np.asarray(list(itertools.product(R,R)))
        diff = pair_array[:,0,:] - pair_array[:,1,:]
        diff = diff.astype(np.float32)

        Runit = diff / np.linalg.norm(diff, axis=1)[:,None]**3 
        Runit = Runit.astype(np.float32)

        #make final dot matrix
        dot_matrix = self.helper_dot(Runit, dR_cross)
        dot_matrix = np.asarray(dot_matrix)
        dot_matrix = dot_matrix.reshape((l-1,l-1))

        nc_gdict = {} 

        for i,j in nc_indexs:

            # loop_range = np.arange(i,j)
            # nterm_range = np.arange(Nterm_thresh,i-5)
            # cterm_range = np.arange(j+6,l-(Cterm_thresh + 1))
            loop_range = np.arange(i, j)
            nterm_range = np.arange(Nterm_thresh, i-loop_Nthread_thresh-1)
            cterm_range = np.arange(j+loop_Cthread_thresh+1, l-(Cterm_thresh + 1))

            gn_pairs_array = np.fromiter(itertools.chain(*itertools.product(nterm_range, loop_range)), int).reshape(-1, 2)
            gc_pairs_array = np.fromiter(itertools.chain(*itertools.product(loop_range, cterm_range)), int).reshape(-1, 2)

            if gn_pairs_array.size != 0:
                
                gn_vals = dot_matrix[gn_pairs_array[:,0],gn_pairs_array[:,1]]
                gn_vals = gn_vals[~np.isnan(gn_vals)] 
                gn_val = np.sum(gn_vals) / (4.0 * np.pi)
            
            else:
                gn_val = 0

            if gc_pairs_array.size != 0:
                
                gc_vals = dot_matrix[gc_pairs_array[:,0],gc_pairs_array[:,1]]
                gc_vals = gc_vals[~np.isnan(gc_vals)] 
                gc_val = np.sum(gc_vals) / (4.0 * np.pi)
            
            else:
                gc_val = 0
            
            rounded_gc_val = self.point_rounding(np.float64(gc_val))
            rounded_gn_val = self.point_rounding(np.float64(gn_val))    

            #if np.abs(rounded_gn_val) >= 1 or np.abs(rounded_gc_val) >= 1:
            #    #print(f'({i}, {j}) with gn: {gn_val} and gc: {gc_val}')
            #    nc_gdict[ (int(i), int(j)) ] = (gn_val, gc_val, rounded_gn_val, rounded_gc_val)
            nc_gdict[ (int(i), int(j)) ] = (gn_val, gc_val, rounded_gn_val, rounded_gc_val)

        ## check for crossings if there are entanglements and topoly==True
        if len(nc_gdict) == 0:
            #print(f'No entanglements found')
            return {}
        
        else:
            if topoly == True:
                entangled_res = self.find_crossing(bb_coor.tolist(), nc_gdict, resids)
                # print(f'entangled_res: {entangled_res}')
                # for k, v in entangled_res.items():
                #     print(f'{k}: {v}')
                # quit()
                return entangled_res
            
            if topoly == False:
                entangled_res = {}
                for ij, values in nc_gdict.items():
                    i,j = ij[0], ij[1]
                    gn = values[0]
                    gc = values[1]
                    GLNn = values[2]
                    GLNc = values[3]
                    TLNn = np.nan
                    TLNc = np.nan
                    entangled_res[(resids[i], resids[j], gn, gc, GLNn, GLNc, TLNn, TLNc)] = [[], []]
                return entangled_res
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def processes_coor(self, coor: np.ndarray, resids: np.ndarray, atom_names: np.ndarray, CG: bool=False, Calpha: bool=False) -> tuple:
        """
        Processes the coordinates of the protein to create a residue level contact map and backbone coordinates.
        If CG is True, it uses the full coor array as the backbone coordinates and the value of Calpha does not matter.
        If CG is False and Calpha is True, it uses the coordinates of the alpha carbons to determine if two residues are in contact.
        If CG is False and Calpha is False, it uses the coordinates of the heavy atoms to determine if two residues are in contact.
        In all cases the backbone coordinates are the alphacarbon coordinates of the residues.
        Returns a tuple of the native contact map and backbone coordinates.
        """
        heavy_atom_names = ['C', 'O', 'N', 'CA', 'CB', 'CG', 'CD', 'CE', 'CZ', 'SD', 'SG'] # heavy atoms for all atom models
        # for idx, c in enumerate(coor):
        #     print(f'Atom {idx}: {atom_names[idx]} at position {c}')
        if CG == True:
            bb_coor = coor
            # make native contact contact map
            dist_matrix = squareform(pdist(bb_coor))
            # print the index pairs and distances where distance is not 0
            # for i in range(dist_matrix.shape[0]):
            #     for j in range(dist_matrix.shape[1]):
            #         print(f'Contact between residues at indexs {i} and {j} at distance {dist_matrix[i, j]} Å')
            native_cmap = np.where(dist_matrix <= 8.0, 1, 0) # if true then 1 will appear otherwise zero
            native_cmap = np.triu(native_cmap, k=4) # element below the 4th diagonal starting from middle are all zeros; # protein contact map
            return native_cmap, bb_coor, dist_matrix
        
        if CG == False:
            CA_idx = [i for i, resname in enumerate(atom_names) if resname == 'CA']
            # print(f'CA_idx: {CA_idx}')
            bb_coor = coor[CA_idx]
            # print(f'bb_coor: {bb_coor.shape} {bb_coor}')

            if Calpha == True:
                # make native contact contact map
                dist_matrix = squareform(pdist(bb_coor))
                # print the index pairs and distances where distance is not 0
                # for i in range(dist_matrix.shape[0]):
                #     for j in range(dist_matrix.shape[1]):
                #         print(f'Contact between residues at indexs {i} and {j} at distance {dist_matrix[i, j]} Å')
                native_cmap = np.where(dist_matrix <= 8.0, 1, 0) # if true then 1 will appear otherwise zero
                native_cmap = np.triu(native_cmap, k=4) # element below the 4th diagonal starting from middle are all zeros; # protein contact map
                return native_cmap, bb_coor, dist_matrix

            if Calpha == False:

                # Select heavy atoms
                heavy_mask = np.isin(atom_names, heavy_atom_names)
                heavy_coor = coor[heavy_mask]
                heavy_resids = resids[heavy_mask]

                # Distance matrix between heavy atoms
                heavy_dist_matrix = squareform(pdist(heavy_coor))
                heavy_native_cmap = heavy_dist_matrix <= 4.5

                # Precompute which atoms belong to which residue
                residue_to_atom_indices = defaultdict(list)
                for idx, resid in enumerate(heavy_resids):
                    residue_to_atom_indices[resid].append(idx)

                # Unique residues and index mapping
                unique_resids = np.unique(resids)
                resid_to_index = {resid: i for i, resid in enumerate(unique_resids)}
                n = len(unique_resids)
                native_cmap = np.zeros((n, n), dtype=int)

                # Only compute upper triangle (excluding diagonal)
                for i, resid_i in enumerate(unique_resids):
                    atoms_i = residue_to_atom_indices.get(resid_i, [])
                    for j in range(i + 4, n):  # skip close-in-sequence residues
                        resid_j = unique_resids[j]
                        atoms_j = residue_to_atom_indices.get(resid_j, [])

                        # Skip if either residue has no heavy atoms
                        if not atoms_i or not atoms_j:
                            continue

                        # Use any contact between atoms of residues i and j
                        contact_exists = np.any(heavy_native_cmap[np.ix_(atoms_i, atoms_j)])
                        if contact_exists:
                            native_cmap[i, j] = 1
                            native_cmap[j, i] = 1
               
                return native_cmap, bb_coor, heavy_dist_matrix
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def combine_ref_traj_GE(self, RefFile: dict, TrajFile: dict, outdir: str='./', ID: str=''):
        """
        This function collects the reference and trajectory entanglements and combines them into a single dataframe and pickle file and saves both. 
        The dataframe is saved as a csv and is easy to read while the .pkl saves as a binary format and saves space and is required for downstream clustering
        of non-native entanglements. 
        The format of this pickle file is as follows

        LEVEL 1: ref, 6600, 6601, ...

        LEVEL 2 for ref: ent_fingerprint, chg_ent_fingerprint, G_dict, G
            ent_fingerprint is a dictionary with the following structure
                KEY: (368, 372) 
                VALUE: {'GLN': [0.04118376194383041, 0.020784770839564926], 
                        'crossing_resid': [[], []], 
                        'crossing_pattern': ['', ''], 
                        'linking_number': [0, 0], 
                        'native_contact': [368, 372]}

            chg_ent_fingerprint = None (because it is a refernece state)

            G_dict = None (because it is a refernece state)

            G = None (because it is a refernece state)

        LEVEL 2 for 6600 or any frame: 
            ent_fingerprint is a dictionary with the following structure
                KEY: (368, 372) 
                VALUE: {'GLN': [0.044431527151880604, 0.010669483331061843], 'crossing_resid': [[], []], 'crossing_pattern': ['', ''], 'linking_number': [0, 0], 'native_contact': [368, 372], 'surrounding_resid': [[], []]}

            chg_ent_fingerprint is a dictionary with the following structure
                KEY: (368, 372) 
                VALUE: {'type': ['no change', 'no change'], 'code': ['L#C#', 'L#C#'], 'GLN': [0.044431527151880604, 0.010669483331061843], 'crossing_resid': [[], []], 'crossing_pattern': ['', ''], 'linking_number': [0, 0], 'native_contact': [368, 372], 'surrounding_resid': [[], []], 'ref_GLN': [0.04118376194383041, 0.020784770839564926], 'ref_crossing_resid': [[], []], 'ref_crossing_pattern': ['', ''], 'ref_linking_number': [0, 0], 'ref_native_contact': [368, 372], 'ref_surrounding_resid': [[], []]}

            G_dict is a dictionary with the following structure
                KEY: change code
                VALUE: count of change code
                EXP: {'L-C~': 0, 'L-C#': 2, 'L+C~': 0, 'L+C#': 0, 'L#C~': 0, 'L#C#': 1630}

            G is the float value of the fraction of native contacts that have a change in entanglement in the frame
        """
        ## set up the outdir for this calculation
        #outdir = f"{os.getcwd()}/{outdir}"
        if not os.path.isdir(outdir):
            os.mkdir(f"{outdir}") 
            self.logger.info(f"Creating directory: {outdir}")

        ## Define the outfile and check if it exists. If so load it else create it
        outfile = os.path.join(f'{outdir}', f'{ID}_GE.pkl')
        #Ref = Ref['ent_result']
        Ref = pd.read_csv(RefFile, sep='|', dtype={'crossingsN': str, 'crossingsC': str})
        Traj = pd.read_csv(TrajFile, sep='|', dtype={'crossingsN': str, 'crossingsC': str})
        #print(Ref.keys())
        #Traj = Traj['ent_result']
        self.logger.info(f'Ref {RefFile}')
        # print(Ref.to_string())
        self.logger.info(f'Traj {TrajFile}')
        # print(Traj.to_string())

        Master = {}
        ##########################################################################################
        ## Get the native state info into the dictionary
        ref = {'ent_fingerprint':{}, 'chg_ent_fingerprint':None, 'G_dict':None, 'G':None}
        Num_native_contacts = len(Ref)
        self.logger.debug(f'Num_native_contacts: {Num_native_contacts}')
        for rowi, row in Ref.iterrows():
            # print(row)
            # Make the key: (i, j)
            i = int(row['i'])
            j = int(row['j'])
            key = (i, j)
            #print(f'\nREF key: {key}')

            # Make the VALUE: {'GLN': [0.041183, 0.02078], 'crossing_resid': [[], []], 'crossing_pattern': ['', ''], 'linking_number': [0, 0], 'native_contact': [368, 372]}
            crossing_resid, crossing_pattern = self.processes_crossings(row)

            gn = float(row['gn'])
            gc = float(row['gc'])
            GLNn = int(row['GLNn'])
            GLNc = int(row['GLNc'])
            TLNn = np.nan if pd.isna(row['TLNn']) else int(row['TLNn'])
            TLNc = np.nan if pd.isna(row['TLNc']) else int(row['TLNc'])
            value = {'linking_value': [gn, gc], 'crossing_resid': crossing_resid, 'crossing_pattern': crossing_pattern, 'gauss_linking_number': [GLNn, GLNc], 'topoly_linking_number': [TLNn, TLNc], 'native_contact': [i, j]}
            # print(f'Ref value: {value}')

            # Update the Master dict
            ref['ent_fingerprint'][key] = value
      
        Master['ref'] = ref

        ##########################################################################################
        ## Get the traj info into the dictionary
        Gdf = {'Frame':[], 'L-C~':[], 'L-C#':[], 'L+C~':[], 'L+C#':[], 'L#C~':[], 'L#C#':[], 'G':[]}
        for frame, frame_df in Traj.groupby('frame'):
            frame_dict = {'ent_fingerprint': {}, 
                     'chg_ent_fingerprint': {}, 
                     'G_dict': {'L-C~': 0, 'L-C#': 0, 'L+C~': 0, 'L+C#': 0, 'L#C~': 0, 'L#C#': 0}, 
                     'G': None}
            
            ## Get the ent_fingerprint data for the frame
            for rowi, row in frame_df.iterrows():
                # print(row)
                # Make the key: (i, j)
                i = int(row['i'])
                j = int(row['j'])
                key = (i, j)
                gn = float(row['gn'])
                gc = float(row['gc'])
                GLNn = int(row['GLNn'])
                GLNc = int(row['GLNc'])
                TLNn = np.nan if pd.isna(row['TLNn']) else int(row['TLNn'])
                TLNc = np.nan if pd.isna(row['TLNc']) else int(row['TLNc'])
                #print(f'\nFRAME {frame} key: {key}')

                # Make the VALUE: {'linking_value': [0.041183, 0.02078], 'crossing_resid': [[], []], 'crossing_pattern': ['', ''], 'gauss_linking_number': [0, 0], 'topoly_linking_number': [0, 0], 'native_contact': [368, 372]}
                crossing_resid, crossing_pattern = self.processes_crossings(row)

                value = {'linking_value': [gn, gc], 'crossing_resid': crossing_resid, 'crossing_pattern': crossing_pattern, 'gauss_linking_number': [GLNn, GLNc], 'topoly_linking_number': [TLNn, TLNc], 'native_contact': [i, j]}
                # print(f'Traj value: {value}')

                # Update the Master dict
                frame_dict['ent_fingerprint'][key] = value

                ## Get the chg_ent_fingerprint data for the frame for those native contacts present
                if key in Master['ref']['ent_fingerprint']:
                    chg_ent_fingerprint = self.get_chg_ent_fingerprint(ref = Master['ref']['ent_fingerprint'][key], frame = value)
                    # print(chg_ent_fingerprint)
                    frame_dict['chg_ent_fingerprint'][key] = chg_ent_fingerprint

                    ## update the G_dict
                    frame_dict['G_dict'][chg_ent_fingerprint['code'][0]] += 1
                    frame_dict['G_dict'][chg_ent_fingerprint['code'][1]] += 1

            ## Calculate G and then update the master dictionary
            # print(f'FRAME {frame} CHANGE SUMMARY: {frame_dict["G_dict"]}')
            Gdf['Frame'] += [frame]
            G = 0
            for code in ['L-C~', 'L-C#', 'L+C~', 'L+C#', 'L#C~']:
                G += frame_dict['G_dict'][code]
                Gdf[code] += [frame_dict['G_dict'][code]]
            Gdf['L#C#'] += [frame_dict['G_dict']['L#C#']]
            G /= (Num_native_contacts*2)
            # print(f'G: {G}')
            Gdf['G'] += [G]
            frame_dict['G'] = G
            Master[frame] = frame_dict
        #print('Master:', Master)

        ## save the master as a pickle file
        with open(outfile, 'wb') as fw:
            pickle.dump(Master, fw)
        self.logger.info(f'SAVED: {outfile}')

        return {'outfile':outfile, 'Combined_ref_traj_dict':Master, 'G':pd.DataFrame(Gdf)}
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def processes_crossings(self, row: pd.Series) -> tuple:
        """
        This function takes the crossing string and processes it into a list of crossing residues and a list of crossing patterns
        1. Split the crossings on , and determine which are N terminal and C terminal and then extracts the crossing number and sign

        Example: 
            i = 15 and j = 101
            crossing_str = +-15,-10,+108,-150 
            
            processes_crossings -> crossing_resid = [[-15, 10], [108, 150]] and crossing_pattern = ['+-', '+-']
        """
        crossing_resid = [[], []]
        crossing_pattern = [[], []]

        crossingsN = row['crossingsN']
        if isinstance(crossingsN, str):
            crossingsN = crossingsN.split(',')
            for crossing in crossingsN:
                if crossing == '?':
                    crossing_resid[0] += []
                    crossing_pattern[0] += ['?']
                else:
                    sign = crossing[0]
                    num = int(crossing[1:])
                    crossing_resid[0] += [num]
                    crossing_pattern[0] += [sign]
            crossing_pattern[0] = ''.join(crossing_pattern[0])
        else:
            crossing_pattern[0] = ''

        crossingsC = row['crossingsC']
        if isinstance(crossingsC, str):
            crossingsC = crossingsC.split(',')
            for crossing in crossingsC:
                if crossing == '?':
                    crossing_resid[1] += []
                    crossing_pattern[1] += ['?']
                else:
                    sign = crossing[0]
                    num = int(crossing[1:])
                    crossing_resid[1] += [num]
                    crossing_pattern[1] += [sign]
            crossing_pattern[1] = ''.join(crossing_pattern[1])
        else:
            crossing_pattern[1] = ''

        return crossing_resid, crossing_pattern
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def get_chg_ent_fingerprint(self, ref: dict, frame: dict) -> dict:
        """
        This function takes the Gn and Gc numbers for a given native contact and determines if ther eis a change in the frame. 
        The types of changes are represented by the following codes: 'L' = 'Linking number'; 'C' = 'Chirality'; '+' = 'Gain'; '-' = 'Loss'; '~' = 'Switch'; '#' = 'No change'. 
        For example, a change with a code "L+C~" refers to a "Gain of entanglement with a switched chirality".
        {'L-C~': 0, 'L-C#': 2, 'L+C~': 0, 'L+C#': 0, 'L#C~': 0, 'L#C#': 1630}

        # G1: L-C~, loss of linking number & switched linking chirality
        # G2: L-C#, loss of linking number & no change of linking chirality
        # G3: L+C~, gain of linking number & switched linking chirality
        # G4: L+C#, gain of linking number & no change of linking chirality
        # G5: L#C~, no change of linking number & switched linking chirality
        # G6: L#C#, loop formed & no change
        # G: Number of change of entanglement (G1+...+G5) / (2 x Number of native contacts in reference structure)
        # Number of native contact in the reference structure: 1066       

        Finally the chg_ent_fingerprint looks like this
        {'type': ['loss of linking number & no change of linking chirality', 'no change'],
        'code': ['L-C#', 'L#C#'],
        'GLN': [-0.1950360202409954, 0.40138207852044694],
        'crossing_resid': [[], []],
        'crossing_pattern': ['', ''],
        'linking_number': [0, 0],
        'native_contact': [267, 316],
        'surrounding_resid': [[], []],
        'ref_GLN': [-0.8566778684626731, 0.38085875463868263],
        'ref_crossing_resid': [[256], []],
        'ref_crossing_pattern': ['-', ''],
        'ref_linking_number': [-1, 0],
        'ref_native_contact': [267, 316]}
        """
        #print(f'ref: {ref}')
        #print(f'frame {frame}')

        ## Check for N/C terminal change — compute only what ent_detection_method requires
        ###-------------------------------------------------------------------------------------
        if self.ent_detection_method in (1, 3):
            GLN_ref_N_G = ref['gauss_linking_number'][0]
            GLN_ref_C_G = ref['gauss_linking_number'][1]
            GLN_frame_N_G = frame['gauss_linking_number'][0]
            GLN_frame_C_G = frame['gauss_linking_number'][1]
            print(f'\nGLN_ref_N_G: {GLN_ref_N_G}, GLN_ref_C_G: {GLN_ref_C_G}, GLN_frame_N_G: {GLN_frame_N_G}, GLN_frame_C_G: {GLN_frame_C_G}')

            ## get the change code for GLN N terminal
            if abs(GLN_frame_N_G) < abs(GLN_ref_N_G):
                GLN_N_link = '-'
            elif abs(GLN_frame_N_G) > abs(GLN_ref_N_G):
                GLN_N_link = '+'
            elif abs(GLN_frame_N_G) == abs(GLN_ref_N_G):
                GLN_N_link = '#'

            ## get the change code for GLN C terminal
            if abs(GLN_frame_C_G) < abs(GLN_ref_C_G):
                GLN_C_link = '-'
            elif abs(GLN_frame_C_G) > abs(GLN_ref_C_G):
                GLN_C_link = '+'
            elif abs(GLN_frame_C_G) == abs(GLN_ref_C_G):
                GLN_C_link = '#'

            # check if the signs of frame_G and ref_G are the same for N terminus
            if GLN_frame_N_G * GLN_ref_N_G >= 0:
                GLN_N_chiral = '#'
            elif GLN_frame_N_G * GLN_ref_N_G < 0:
                GLN_N_chiral = '~'

            # check if the signs of frame_G and ref_G are the same for C terminus
            if GLN_frame_C_G * GLN_ref_C_G >= 0:
                GLN_C_chiral = '#'
            elif GLN_frame_C_G * GLN_ref_C_G < 0:
                GLN_C_chiral = '~'

            print(f'GLN_N_link: {GLN_N_link}, GLN_N_chiral: {GLN_N_chiral}, GLN_C_link: {GLN_C_link}, GLN_C_chiral: {GLN_C_chiral}')
        ###-------------------------------------------------------------------------------------

        ###-------------------------------------------------------------------------------------
        if self.ent_detection_method in (2, 3):
            TLN_ref_N_G = ref['topoly_linking_number'][0]
            TLN_ref_C_G = ref['topoly_linking_number'][1]
            TLN_frame_N_G = frame['topoly_linking_number'][0]
            TLN_frame_C_G = frame['topoly_linking_number'][1]
            print(f'TLN_ref_N_G: {TLN_ref_N_G}, TLN_ref_C_G: {TLN_ref_C_G}, TLN_frame_N_G: {TLN_frame_N_G}, TLN_frame_C_G: {TLN_frame_C_G}')

            ## get the change code for TLN N terminal
            if abs(TLN_frame_N_G) < abs(TLN_ref_N_G):
                TLN_N_link = '-'
            elif abs(TLN_frame_N_G) > abs(TLN_ref_N_G):
                TLN_N_link = '+'
            elif abs(TLN_frame_N_G) == abs(TLN_ref_N_G):
                TLN_N_link = '#'

            ## get the change code for TLN C terminal
            if abs(TLN_frame_C_G) < abs(TLN_ref_C_G):
                TLN_C_link = '-'
            elif abs(TLN_frame_C_G) > abs(TLN_ref_C_G):
                TLN_C_link = '+'
            elif abs(TLN_frame_C_G) == abs(TLN_ref_C_G):
                TLN_C_link = '#'

            # check if the signs of frame_G and ref_G are the same for N terminus
            if TLN_frame_N_G * TLN_ref_N_G >= 0:
                TLN_N_chiral = '#'
            elif TLN_frame_N_G * TLN_ref_N_G < 0:
                TLN_N_chiral = '~'

            # check if the signs of frame_G and ref_G are the same for C terminus
            if TLN_frame_C_G * TLN_ref_C_G >= 0:
                TLN_C_chiral = '#'
            elif TLN_frame_C_G * TLN_ref_C_G < 0:
                TLN_C_chiral = '~'

            print(f'TLN_N_link: {TLN_N_link}, TLN_N_chiral: {TLN_N_chiral}, TLN_C_link: {TLN_C_link}, TLN_C_chiral: {TLN_C_chiral}')
        ###-------------------------------------------------------------------------------------

        ###-------------------------------------------------------------------------------------
        ## Determine the overall change code for the frame based on the ent_detection_method
        Ncode = ''
        Ccode = ''
        if self.ent_detection_method == 1:
            # Any nonzero GLN for either termini
            Ncode = f'L{GLN_N_link}C{GLN_N_chiral}'
            Ccode = f'L{GLN_C_link}C{GLN_C_chiral}'
            codes = [Ncode, Ccode]

            Ntype = self.change_codes[Ncode]
            Ctype = self.change_codes[Ccode]
            types = [Ntype, Ctype]


        elif self.ent_detection_method == 2:
            # Any nonzero TLN for either termini (default)
            Ncode = f'L{TLN_N_link}C{TLN_N_chiral}'
            Ccode = f'L{TLN_C_link}C{TLN_C_chiral}'
            codes = [Ncode, Ccode]

            Ntype = self.change_codes[Ncode]
            Ctype = self.change_codes[Ccode]
            types = [Ntype, Ctype]

        elif self.ent_detection_method == 3:
            # if both GLN and TLN changes are the same then use that code, if they are different then use the GLN code (because it is more sensitive)
            GLN_N_code = f'L{GLN_N_link}C{GLN_N_chiral}'
            GLN_C_code = f'L{GLN_C_link}C{GLN_C_chiral}'
            TLN_N_code = f'L{TLN_N_link}C{TLN_N_chiral}'
            TLN_C_code = f'L{TLN_C_link}C{TLN_C_chiral}'
            if GLN_N_code == TLN_N_code:
                Ncode = GLN_N_code
            else:
                Ncode = 'L#C#'
            if GLN_C_code == TLN_C_code:
                Ccode = GLN_C_code
            else:                
                Ccode = 'L#C#'
            codes = [Ncode, Ccode]

            Ntype = self.change_codes[Ncode]
            Ctype = self.change_codes[Ccode]
            types = [Ntype, Ctype]
        ###-------------------------------------------------------------------------------------

        #print(f'codes: {codes}')
        #print(f'types: {types}')

        chg_ent_fingerprint = {'type': types,
        'code': codes,
        'native_contact': frame['native_contact'],
        'linking_value': frame['linking_value'],
        'crossing_resid': frame['crossing_resid'],
        'crossing_pattern': frame['crossing_pattern'],
        'gauss_linking_number': frame['gauss_linking_number'],
        'topoly_linking_number': frame['topoly_linking_number'],
        'ref_native_contact': ref['native_contact'],
        'ref_linking_value': ref['linking_value'],
        'ref_crossing_resid': ref['crossing_resid'],
        'ref_crossing_pattern': ref['crossing_pattern'],
        'ref_gauss_linking_number': ref['gauss_linking_number'],
        'ref_topoly_linking_number': ref['topoly_linking_number'], 
        'ent_detection_method': self.ent_detection_method}

        return chg_ent_fingerprint
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def select_high_quality_entanglements(self, GE_filepath: str, pdb: str, outdir: str='./', ID: str='', model: str='EXP', mapping: str='None', chain: str=None) -> dict:
        """
        This function takes the GE file and selects the high quality entanglements based on the following criteria:
        1. Remove any native NCLE's that are predicted to be pure slipknots (crossings with a net sign that cancels out)
        2. if the model is EXP it will try and only grab those ENT that are mapped to a uniprot sequence if the user specifies a mapping file
        3. if the model is AF then also check that the i, j, and k meet our criteria
        """
        ## set up the outdir for this calculation
        #outdir = f"{os.getcwd()}/{outdir}"
        if not os.path.isdir(outdir):
            os.mkdir(f"{outdir}") 
            self.logger.info(f"Creating directory: {outdir}")

        ## load the dataframe
        GE_data = pd.read_csv(GE_filepath, sep='|', dtype={'crossingsN': str, 'crossingsC': str})
        GE_data = GE_data[GE_data['ENT'] == True].reset_index(drop=True)
        self.logger.info(f'GE FILE: {GE_filepath}')
        # print(f'RAW GE_data:\n{GE_data}')

        ## select only those entanglements that are mapped for the EXP model
        if model == 'EXP' and mapping != 'None':
            GE_data = self.remove_slipknots(GE_data)
            #print(f'No Slipknot GE_data:\n{GE_data}')

            GE_data = self.mappingPDB2Uniprot(GE_data, mapping)
            #print(f'No Slipknot mapped GE_data:\n{GE_data}')            

        ## select only those entanglements that are mapped for the EXP model
        if model == 'EXP' and mapping == 'None':
            GE_data = self.remove_slipknots(GE_data)
            # print(f'No Slipknot GE_data:\n{GE_data}')
          
        ## select only those entanglements that meet our pLDDT thresholds for AF model
        if model == 'AF':
            #GE_data = self.remove_slipknots(GE_data)
            #print(f'No Slipknot GE_data:\n{GE_data}')
            GE_data = self.remove_low_quality_AF_entanglements(GE_data, pdb)
            #print(f'No Slipknot GE_data:\n{GE_data}')


        outfile = os.path.join(outdir, f'{ID}.csv')
        GE_data.to_csv(outfile, index=False, sep='|')
        self.logger.info(f'SAVED: {outfile}')
        GE_data = pd.read_csv(outfile, sep='|', dtype={'c': str})
        # print(f'HQ GE_data:\n{GE_data}')
        return {'outfile':outfile, 'GE_data':GE_data}
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def remove_slipknots(self, df):
        """
        Checks each raw entanglement for crossings where the signs sum to a net of 0. 
        """        
        self.logger.info(f'\n{"#"*50}\nRemoving slipknots...')
        new_df = {'ID':[], 'chain':[], 'i':[], 'j':[], 'crossingsN':[], 'crossingsC':[], 'gn':[], 'gc':[], 'GLNn':[], 'GLNc':[], 'TLNn':[], 'TLNc':[], 'CCbond':[], 'ENT':[], 'Slipknot_N':[], 'Slipknot_C':[]}
        for rowi, row in df.iterrows():
            # print(row)
            ID = row['ID']
            chain = row['chain']
            i = row['i']
            j = row['j']
            if isinstance(row['crossingsN'], float):
                rN = ['']
            else:
                rN = row['crossingsN'].split(',')
                ENT = row['ENT']
            if isinstance(row['crossingsC'], float):
                rC = ['']
            else:
                rC = row['crossingsC'].split(',')
                ENT = row['ENT']
            gn = row['gn']
            gc = row['gc']
            GLNn = row['GLNn']
            GLNc = row['GLNc']
            TLNn = row['TLNn']
            TLNc = row['TLNc']
            CCbond = row['CCbond']
            ENT = row['ENT']
            # print(ID, i, j, rN, rC, gn, gc, GLNn, GLNc, TLNn, TLNc, CCbond, ENT)
        

            # Check each termini for slipknots
            slipknot_dict = {'N': False, 'C': False}
            for termini, crossing_list in {'N': rN, 'C': rC}.items():
                if len(crossing_list) > 1:
                    crossings_signs = []
                    crossings = []
                    crossings_resids = []
                    for cross in crossing_list:
                        #print(cross)
                        cross_sign = cross[0]
                        cross_int = int(cross[1:])
                        #print(i, j, cross)

                        #check if the crossing in N terminal
                        if cross_int < i:
                            crossings += [cross]
                            crossings_resids += [cross_int]
                            if cross_sign == '+':
                                crossings_signs += [1]
                            elif cross_sign == '-':
                                crossings_signs += [-1]
                            else:
                                raise ValueError(f'The crossing sign was not + or - {cross}')
                    
                    # check if either termini has duplicate crossings and empty out those crossings lists as we will not have any confidence
                    if len(crossings_resids) != len(set(crossings_resids)):
                        raise ValueError(f'Duplicate crossings found in N terminus: {crossings_resids}')

                    # get the sum of the N terminal crossings
                    if len(crossings_signs) != 0:
                        slipknot_dict[termini] = True
            # print(f'slipknot_dict: {slipknot_dict}')

            # update the new df 
            new_df['i'] += [i]
            new_df['j'] += [j]
            new_df['ID'] += [ID]
            new_df['chain'] += [chain]
            new_df['crossingsN'] += [','.join(rN)]
            new_df['crossingsC'] += [','.join(rC)]
            new_df['gn'] += [gn]
            new_df['gc'] += [gc]
            new_df['GLNn'] += [GLNn]
            new_df['GLNc'] += [GLNc]
            new_df['TLNn'] += [TLNn]
            new_df['TLNc'] += [TLNc]
            new_df['CCbond'] += [CCbond]
            new_df['ENT'] += [ENT]
            new_df['Slipknot_N'] += [slipknot_dict['N']]
            new_df['Slipknot_C'] += [slipknot_dict['C']]
        
        new_df = pd.DataFrame(new_df)
        return new_df
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def mappingPDB2Uniprot(self, df, mapping):
        """
        Maps the PDB level ENT to the uniprot resid if desired
        """
        ## Check if the mapping file exists
        if os.path.exists(mapping):
            mapping = np.loadtxt(mapping, dtype='O')
            mapping = np.vstack([x[1:] for x in mapping if ('Mapped' in x[0] or 'Modifed_Residue' in x[0] or 'Missense' in x[0])]).astype(int)
            mapping_pdb2uniprot = {pdb:uni for pdb, uni in mapping}
        else:
            raise ValueError(f'Mapping file {mapping} could not be found!')

        new_df = {'ID':[], 'chain':[], 'i':[], 'j':[], 'crossingsN':[], 'crossingsC':[], 'gn':[], 'gc':[], 'GLNn':[], 'GLNc':[], 'TLNn':[], 'TLNc':[], 'CCbond':[], 'ENT':[]}
        for rowi, row in df.iterrows():
            #print(row)
            ID = row['ID']
            chain = row['chain']
            i = row['i']
            j = row['j']
            
            # Parse crossings from crossingsN and crossingsC
            crossingsN = row['crossingsN'] if pd.notna(row['crossingsN']) and row['crossingsN'] != '' and row['crossingsN'] != '?' else ''
            crossingsC = row['crossingsC'] if pd.notna(row['crossingsC']) and row['crossingsC'] != '' and row['crossingsC'] != '?' else ''
            
            rN = crossingsN.split(',') if crossingsN else []
            rC = crossingsC.split(',') if crossingsC else []
            r = rN + rC
            
            ENT = row['ENT']
            crossings = [int(c[1:]) for c in r if c != '']
            
            gn = row['gn']
            gc = row['gc']
            GLNn = row['GLNn']
            GLNc = row['GLNc']
            TLNn = row['TLNn']
            TLNc = row['TLNc']
            CCbond = row['CCbond']
            key_res = [i, j] + crossings
            #print(ID, i, j, r, crossings, gn, gc, GLNn, GLNc, TLNn, TLNc, CCbond, ENT, key_res)
            
            mapped = True
            for res in key_res:
                if res not in mapping_pdb2uniprot:
                    self.logger.debug(f'Res: {res} not mapped! this entanglement will be discarded from {ID}')
                    mapped = False
                                    
            if mapped:
                new_df['ID'] += [ID]
                new_df['chain'] += [chain]
                new_df['i'] += [i]
                new_df['j'] += [j]
                new_df['crossingsN'] += [crossingsN]
                new_df['crossingsC'] += [crossingsC]
                new_df['gn'] += [gn]
                new_df['gc'] += [gc]
                new_df['GLNn'] += [GLNn]
                new_df['GLNc'] += [GLNc]
                new_df['TLNn'] += [TLNn]
                new_df['TLNc'] += [TLNc]
                new_df['CCbond'] += [CCbond]
                new_df['ENT'] += [ENT]
            else:
                self.logger.info(f'Entanglement {rowi} was not mapped and will be discarded')

        new_df = pd.DataFrame(new_df)
        return new_df
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def remove_low_quality_AF_entanglements(self, df, pdb):
        """
        # (1) check if both i and j have pLDDt >= 70. if so continue else completely ignore the ent
        # (2) starting from the loop base get the set of ordered crossings that have pLDDT > 70 and discard any after the first crossings that fails this. 
        """
        self.logger.info(f'\n{"#"*50}\nRemoving low quality AF entanglements...')
        avg_pLDDT, pLDDT_df = self.average_pLDDT(pdb)
        #print(f'avg_pLDDT: {avg_pLDDT}\n{pLDDT_df}')

        new_df = {'ID':[], 'chain':[], 'i':[], 'j':[], 'crossingsN':[], 'crossingsC':[], 'gn':[], 'gc':[], 'GLNn':[], 'GLNc':[], 'TLNn':[], 'TLNc':[], 'CCbond':[], 'ENT':[], 'Quality':[], 'Reason':[]}
        for rowi, row in df.iterrows():
            # print(row)
            ID = row['ID']
            chain = row['chain']
            i = row['i']
            j = row['j']
            
            # Parse crossings from crossingsN and crossingsC
            crossingsN = row['crossingsN'] if pd.notna(row['crossingsN']) and row['crossingsN'] != '' and row['crossingsN'] != '?' else ''
            crossingsC = row['crossingsC'] if pd.notna(row['crossingsC']) and row['crossingsC'] != '' and row['crossingsC'] != '?' else ''
            
            rN = crossingsN.split(',') if crossingsN else []
            rC = crossingsC.split(',') if crossingsC else []
            # Filter out empty strings from combined list
            r = [x for x in (rN + rC) if x != '']
            
            ENT = row['ENT']
            gn = row['gn']
            gc = row['gc']
            GLNn = row['GLNn']
            GLNc = row['GLNc']
            TLNn = row['TLNn']
            TLNc = row['TLNc']
            CCbond = row['CCbond']
            #print(ID, i, j, r, gn, gc, GLNn, GLNc, TLNn, TLNc, CCbond, ENT)


            # (1) check if both i and j have pLDDt >= 70. if so continue else completely ignore the ent
            NC_pLDDT = pLDDT_df[pLDDT_df['resid'].isin([i,j])]['pLDDT'].values
            if all(NC_pLDDT >= 70):
                #print(f'Native contact pLDDT are >= 70 {NC_pLDDT}')
                if ENT == False: ## if no entanglement but the native contact has a pLDDT >= 70 still return the contact as HQ even though there is no ent
                    new_df['ID'] += [ID]
                    new_df['chain'] += [chain]
                    new_df['i'] += [i]
                    new_df['j'] += [j]
                    new_df['crossingsN'] += [crossingsN]
                    new_df['crossingsC'] += [crossingsC]
                    new_df['gn'] += [gn]
                    new_df['gc'] += [gc]
                    new_df['GLNn'] += [GLNn]
                    new_df['GLNc'] += [GLNc]
                    new_df['TLNn'] += [TLNn]
                    new_df['TLNc'] += [TLNc]
                    new_df['CCbond'] += [CCbond]
                    new_df['ENT'] += [ENT]
                    new_df['Quality'] += ['High']
                    new_df['Reason'] += ['NC pLDDT >= 70']  
                    continue
            else:
                #print(f'Native contact pLDDT are < 70 {NC_pLDDT}')
                new_df['ID'] += [ID]
                new_df['chain'] += [chain]
                new_df['i'] += [i]
                new_df['j'] += [j]
                new_df['crossingsN'] += [crossingsN]
                new_df['crossingsC'] += [crossingsC]
                new_df['gn'] += [gn]
                new_df['gc'] += [gc]
                new_df['GLNn'] += [GLNn]
                new_df['GLNc'] += [GLNc]
                new_df['TLNn'] += [TLNn]
                new_df['TLNc'] += [TLNc]
                new_df['CCbond'] += [CCbond]
                new_df['ENT'] += [ENT]
                new_df['Quality'] += ['Low']
                new_df['Reason'] += ['NC pLDDT < 70']
                continue
                
            # (2) starting from the loop base get the set of ordered crossings that have pLDDT > 70 and discard any after the first crossings that fails this. 
            #print(f'Getting HQ N-terminal entanglements')
            Ncrossings_resids, Ncrossings, Ncrossings_signs = self.parse_crossings(r, i=i, j=j, term='N')
            #print(Ncrossings_resids, Ncrossings, Ncrossings_signs)
            Ncrossings_resids, Ncrossings, Ncrossings_signs, Ndup_flag = self.remove_duplicates(Ncrossings_resids, Ncrossings, Ncrossings_signs)
            HQ_Ncrossings = []
            HQ_Ncrossings_resids = []
            HQ_Ncrossings_signs = []
            if len(Ncrossings_resids) != 0:
                sorted_indices = np.argsort(Ncrossings_resids)[::-1]  # [::-1] reverses the order
                Ncrossings_resids = [Ncrossings_resids[i] for i in sorted_indices]
                Ncrossings = [Ncrossings[i] for i in sorted_indices]
                Ncrossings_signs = [Ncrossings_signs[i] for i in sorted_indices]
                Ncrossings_pLDDTs = pLDDT_df[pLDDT_df['resid'].isin(Ncrossings_resids)]['pLDDT'].values
                #print(Ncrossings_resids, Ncrossings, Ncrossings_signs, Ncrossings_pLDDTs)
                for cross_i, cross in enumerate(Ncrossings_resids):
                    if Ncrossings_pLDDTs[cross_i] >=70:
                        HQ_Ncrossings += [Ncrossings[cross_i]]
                        HQ_Ncrossings_resids += [cross]
                        HQ_Ncrossings_signs += [Ncrossings_signs[cross_i]]
                    else:
                        break
                
                # check that the remaining HQ Nterminal crossigns are not a slipknot
                if sum(HQ_Ncrossings_signs) == 0:
                    #print(f'SlipKnot foundin N terminus after HQ search')
                    HQ_Ncrossings = []
                    HQ_Ncrossings_resids = []
                    HQ_Ncrossings_signs = []
            
            
            #print(f'Getting HQ C-terminal entanglements')
            Ccrossings_resids, Ccrossings, Ccrossings_signs = self.parse_crossings(r, i=i, j=j, term='C')
            #print(Ccrossings_resids, Ccrossings, Ccrossings_signs)
            Ccrossings_resids, Ccrossings, Ccrossings_signs, Cdup_flag = self.remove_duplicates(Ccrossings_resids, Ccrossings, Ccrossings_signs)
            HQ_Ccrossings = []
            HQ_Ccrossings_resids = []
            HQ_Ccrossings_signs = []
            if len(Ccrossings_resids) != 0:
                sorted_indices = np.argsort(Ccrossings_resids) # [::-1] reverses the order
                Ccrossings_resids = [Ccrossings_resids[i] for i in sorted_indices]
                Ccrossings = [Ccrossings[i] for i in sorted_indices]
                Ccrossings_signs = [Ccrossings_signs[i] for i in sorted_indices]
                Ccrossings_pLDDTs = pLDDT_df[pLDDT_df['resid'].isin(Ccrossings_resids)]['pLDDT'].values
                #print(Ccrossings_resids, Ccrossings, Ccrossings_signs, Ccrossings_pLDDTs)
                for cross_i, cross in enumerate(Ccrossings_resids):
                    if Ccrossings_pLDDTs[cross_i] >=70:
                        HQ_Ccrossings += [Ccrossings[cross_i]]
                        HQ_Ccrossings_resids += [cross]
                        HQ_Ccrossings_signs += [Ccrossings_signs[cross_i]]
                    else:
                        break

                # check that the remaining HQ Cterminal crossigns are not a slipknot
                if sum(HQ_Ccrossings_signs) == 0:
                    #print(f'SlipKnot foundin N terminus after HQ search')
                    HQ_Ccrossings = []
                    HQ_Ccrossings_resids = []
                    HQ_Ccrossings_signs = []

            HQ_crossings = HQ_Ccrossings + HQ_Ncrossings
            HQ_crossings = np.asarray(HQ_crossings, dtype=str)
            #print(f'HQ_crossings: {HQ_crossings}')
            
            # Separate HQ crossings into N and C terminal
            HQ_crossingsN = ','.join(HQ_Ncrossings) if len(HQ_Ncrossings) > 0 else ''
            HQ_crossingsC = ','.join(HQ_Ccrossings) if len(HQ_Ccrossings) > 0 else ''

            if len(HQ_crossings) != 0:
                    new_df['ID'] += [ID]
                    new_df['chain'] += [chain]
                    new_df['i'] += [i]
                    new_df['j'] += [j]
                    new_df['crossingsN'] += [HQ_crossingsN]
                    new_df['crossingsC'] += [HQ_crossingsC]
                    new_df['gn'] += [gn]
                    new_df['gc'] += [gc]
                    new_df['GLNn'] += [GLNn]
                    new_df['GLNc'] += [GLNc]
                    new_df['TLNn'] += [TLNn]
                    new_df['TLNc'] += [TLNc]
                    new_df['CCbond'] += [CCbond]
                    new_df['ENT'] += [ENT]
                    new_df['Quality'] += ['High']
                    new_df['Reason'] += ['NC and Crossings pLDDT >= 70']  
            else:
                #print(f'No crossing remaining in either the N or C terminus after removals. ENT will be ignored')
                new_df['ID'] += [ID]
                new_df['chain'] += [chain]
                new_df['i'] += [i]
                new_df['j'] += [j]
                new_df['crossingsN'] += [crossingsN]
                new_df['crossingsC'] += [crossingsC]
                new_df['gn'] += [gn]
                new_df['gc'] += [gc]
                new_df['GLNn'] += [GLNn]
                new_df['GLNc'] += [GLNc]
                new_df['TLNn'] += [TLNn]
                new_df['TLNc'] += [TLNc]
                new_df['CCbond'] += [CCbond]
                new_df['ENT'] += [ENT]
                new_df['Quality'] += ['Low']
                new_df['Reason'] += ['No HQ crossings']  

        new_df = pd.DataFrame(new_df)
        return new_df
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def parse_crossings(self, r:list, i=int, j=int, term:str='N'):
        crossings_resids, crossings, crossings_signs = [], [], []
        #print(i, j, r, term)
        for cross in r:
            if cross[0] == '+':
                cross_sign = 1
            elif cross[0] == '-':
                cross_sign = -1
            else:
                cross_sign = 0

            cross_resid = int(cross[1:])

            if term == 'N':
                if cross_resid < i:
                    crossings_resids += [cross_resid]
                    crossings += [cross]
                    crossings_signs += [cross_sign]
            if term == 'C':
                if cross_resid > j:
                    crossings_resids += [cross_resid]
                    crossings += [cross]
                    crossings_signs += [cross_sign]      

        return crossings_resids, crossings, crossings_signs
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def remove_duplicates(self, crossings_resids, crossings, crossings_signs):
        # Create a list to keep track of unique elements in A and their corresponding elements in B
        unique_crossings_resids = []
        unique_crossings = []
        unique_crossings_signs = []

        # Dictionary to count occurrences of elements in A
        counts_A = {item: crossings_resids.count(item) for item in crossings_resids}

        # Remove elements from both lists where duplicates are found in A
        dup_flag = False
        for i in range(len(crossings_resids)):
            if counts_A[crossings_resids[i]] == 1:  # If the element in A is unique
                unique_crossings_resids.append(crossings_resids[i])
                unique_crossings.append(crossings[i])
                unique_crossings_signs.append(crossings_signs[i])
            else:
                unique_crossings_resids = []
                unique_crossings = []
                unique_crossings_signs = []
                dup_flag = True
                break

        return unique_crossings_resids, unique_crossings, unique_crossings_signs, dup_flag
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def average_pLDDT(self, pdb_filename):
        
        # Create a PDB parser
        parser = PDBParser(QUIET=True)
        
        # Parse the PDB structure
        structure = parser.get_structure('PDB_structure', pdb_filename)
        
        # List to hold pLDDT
        pLDDTs = []
        pLDDT_df = {'resid':[], 'pLDDT':[]}
        # Iterate over all atoms in the structure to extract pLDDT
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        if atom.get_name() == 'CA':
                            pLDDTs.append(atom.get_bfactor())
                            pLDDT_df['resid'] += [residue.get_id()[1]]
                            pLDDT_df['pLDDT'] += [atom.get_bfactor()]
        
        # Calculate the average pLDDT
        if len(pLDDTs) > 0:
            avg_pLDDT = sum(pLDDTs) / len(pLDDTs)
            pLDDT_df = pd.DataFrame(pLDDT_df)
            return avg_pLDDT, pLDDT_df
        else:
            return None   
    ##########################################################################################################################################################
