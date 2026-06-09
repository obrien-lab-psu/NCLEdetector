#!/usr/bin/env python3
from collections import defaultdict
import numpy as np
import itertools
from geom_median.numpy import compute_geometric_median
from scipy.spatial.distance import cdist, squareform
from functools import cache
import re
import random
import pandas as pd
import pickle
import logging
import sys, getopt, math, os, time, traceback, glob, copy
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from EntDetect._logging import setup_logger
from scipy.cluster.hierarchy import fcluster, linkage, cophenet
try:
    import parmed as pmd
    import mdtraj as mdt
except ImportError:
    pmd = None
    mdt = None
import matplotlib
import matplotlib.pyplot as plt
try:
    import pyemma as pem
    import deeptime
except ImportError:
    pem = None
    deeptime = None
from matplotlib.cm import get_cmap
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.colors as mcolors
import seaborn as sns
from scipy.stats import mode
import pathlib
from dataclasses import dataclass
import json

matplotlib.use('Agg')
pd.set_option('display.max_rows', 500)

class ClusterNativeEntanglements:
    """
    Class to calculate native entanglements given either a file path to an entanglement file or an entanglement object
    """

    ##########################################################################################################################################################
    def __init__(self, organism: str = 'Ecoli', cut_off: int = None, outdir: str = None, log_level: int = logging.INFO, logdir: str = None) -> None:
        """
        Constructor for GaussianEntanglement class.

        Parameters
        ----------
        """

        if organism == 'Human':
            self.cut_off = 52
        elif organism == 'Ecoli':
            self.cut_off = 57
        elif organism == 'Yeast':
            self.cut_off = 49
        
        if cut_off is not None:
            self.cut_off = cut_off
        self.organism = organism
        self.logger = setup_logger('ClusterNativeEntanglements', outdir=logdir if logdir is not None else outdir, log_level=log_level)
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def loop_distance(self, entangled_A: tuple, entangled_B: tuple):

        # remove chiralites then perform euclidean distance
        new_cr_A = [int(cr_A[1:]) for cr_A in entangled_A[3:-1]]
        new_entangled_A = (entangled_A[0], entangled_A[1], entangled_A[2], *new_cr_A)

        new_cr_B = [int(cr_B[1:]) for cr_B in entangled_B[3:-1]]
        new_entangled_B = (entangled_B[0], entangled_B[1], entangled_B[2], *new_cr_B)

        return math.dist(new_entangled_A[1:], new_entangled_B[1:])
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def check_step_ij_kl_range(self, ent1: tuple, ent2: tuple):

        # check if i or j of (i,j) reside within the range (inclusive) of (k,l), and vice versa

        i1, j1 = ent1[1], ent1[2]
        i2, j2 = ent2[1], ent2[2]
        return (i2 <= i1 <= j2) or (i2 <= j1 <= j2) or (i1 <= i2 <= j1) or (i1 <= j2 <= j1)
    ##########################################################################################################################################################

    ##########################################################################################################################################################

    ##########################################################################################################################################################
    #@cache
    def Cluster_NativeEntanglements(self, GE_filepath: str, outdir: str='./', outfile: str='Cluster_NativeEntanglements.txt', chain: str=None):
        
        """
        PARAMS:
            GE_file: str 
            cut_off: int

        1. Identify all unique "residue crossing set and chiralites"

            1b. sort the residues along with the chiralities
        
        2. Find the minimal loop encompassing a given "residue crossing set and chiralites"
            
            2b
            
            i. Identify entanglements that have any crossing residues between them that are 
                less than or equal to 3 residues apart and have the same chirality.
        
            ii. Then check if i or j of (i,j) reside within the range (inclusive) of (k,l), or vice versa;
            
            iii. If yes, then check if any crossing residues are in the range of min(i,j,k,l) to max(i,j,k,l); 
            if yes, skip rest of 2
            
            iv. If no, then check if the number of crossing residues, in each residue set, are different;
            
            v. All crossing residue(s) in the entanglement with the fewer crossings need to have a distance <= 20 
            with the crossings in the other entanglement. Do this by the "brute force" approach and 
            the true distance formula. This means, calculate the distances and take the minimal distance 
            as the distance you check that is less than or equal to 20.
            
            If yes, then keep the {i,j} {r} that have the greatest number of crossing residues;
            If not, then keep the two entanglements separate. 
        
        3. For at least two entanglements each with 1 or more crossings. 
            Loop over the entanglments two at time (avoid double counting)
            Check if i or j of (i,j) reside within the range (inclusive) of (k,l), or  vice versa;
            If yes, check if number of crossing residues is the same (and it is 1 or more)
            If yes, calculate the distances between all crossing residues
                and if both have the same chiralities. 
                (Do NOT use brute force, just compare 1-to-1 index of crossing residues).
            If all the distances are less than or equal to 20, then determine which 
                entanglement has the smaller loop, remove the entanglement with the larger loop
                
        4. Spatially cluster those outputs that have (i) the same number of crossings and (ii) the same chiralities

        """
        self.logger.info(f'Clustering {self.organism} Native Entanglements with dist_cutoff: {self.cut_off}')
        GE_file = GE_filepath.split('/')[-1]
        self.logger.debug(f'{GE_file} cut_off={self.cut_off} outdir={outdir}')

        full_entanglement_data = defaultdict(list)

        ent_data = defaultdict(list)

        rep_ID_ent = defaultdict(list) 

        grouped_entanglement_data = defaultdict(list)

        Before_cr_dist = defaultdict(list)
        
        After_cr_dist = defaultdict(list)

        entanglement_partial_g_data = {}
        
        ## Check if the clustering file is already made and if so use it
        outfilepath = os.path.join(f'{outdir}', f'{outfile}')
        if os.path.exists(outfilepath):
            self.logger.info(f'{outfilepath} ALREADY EXISTS AND WILL BE LOADED')
            outdf = pd.read_csv(outfilepath, sep='|')
            return {'outfile':outfilepath, 'ent_result':outdf}

        self.logger.info(f'Loading {GE_filepath}')
        GE_data = pd.read_csv(GE_filepath, sep='|', dtype={'crossingsN': str, 'crossingsC': str})
        GE_data = GE_data[GE_data['ENT'] == True].reset_index(drop=True)
        # if Quality is a column name then only get the High Quality raw entanglements
        if 'Quality' in GE_data.keys():
            GE_data = GE_data[GE_data['Quality'] == 'High'].reset_index(drop=True)
        
        GE_data = GE_data.replace(np.nan, '', regex=True)
        self.logger.debug(GE_data)
      
        ### STEP 1 INITAL LOADING AND MERGING ################################################################################################################
        ############################################################################################
        ## parse the entanglement file and
        ## get those native contacts that are disulfide bonds
        self.logger.info(f'# Step 1')
        CCBonds = []
        num_raw_ents = {}
        chain_info = {}  # Track chain for each ID
        for rowi, row in GE_data.iterrows():
            # print(row)
   
            ID = row['ID']
            chain = row['chain'] if 'chain' in row else 'A'  # Default to 'A' if chain not present
            chain_info[ID] = chain
            native_contact_i, native_contact_j = row['i'], row['j']
            
            # Store separate N and C crossings
            crossingsN = row['crossingsN'] if pd.notna(row['crossingsN']) and row['crossingsN'] != '' else ''
            crossingsC = row['crossingsC'] if pd.notna(row['crossingsC']) and row['crossingsC'] != '' else ''
            
            crossing_res = [cr for cr in [crossingsN, crossingsC] if cr != '']
            crossing_res = ','.join(crossing_res)
            self.logger.debug(f'{native_contact_i} {native_contact_j} {crossing_res}')

            gn, gc = row['gn'], row['gc']
            GLNn, GLNc = row['GLNn'], row['GLNc']
            TLNn, TLNc = row['TLNn'], row['TLNc']
            # Handle empty strings from NaN replacement: convert to np.nan
            TLNn = np.nan if (isinstance(TLNn, str) and TLNn == '') else TLNn
            TLNc = np.nan if (isinstance(TLNc, str) and TLNc == '') else TLNc
            CCbond = row['CCbond']

            # keep track of number of raw ents for QC purposes
            if ID not in num_raw_ents:
                num_raw_ents[ID] = 1
            else:
                num_raw_ents[ID] += 1

            #native_contact_i, native_contact_j, crossing_res = line[1], line[2], line[3]
            #native_contact_i = int(native_contact_i)
            #native_contact_j = int(native_contact_j)

            reformat_cr = crossing_res.split(',') if crossing_res else []
            
            # Filter out any empty strings from the split
            reformat_cr = [cr for cr in reformat_cr if cr]

            if reformat_cr:
                reformat_cr = sorted(reformat_cr, key = lambda x: int(re.split("\\+|-", x, maxsplit= 1)[1]))
            #print(native_contact_i, native_contact_j, reformat_cr)


            # Step 1 and 1b
            grouped_entanglement_data[(ID, *reformat_cr)].append((native_contact_i, native_contact_j))

            entanglement_partial_g_data[(native_contact_i, native_contact_j, *reformat_cr)] = (gn, gc, GLNn, GLNc, TLNn, TLNc, crossingsN, crossingsC)

            #print(f'CCbond: {CCbond}')
            if CCbond == True:
                CCBonds += [(native_contact_i, native_contact_j)]

        #print(f'num_raw_ents: {num_raw_ents}')        

        #print(f'Step 1 results')
        Step1_QC_counter = 0
        for k,v in grouped_entanglement_data.items():
            #print(k,v)
            Step1_QC_counter += len(v)
        
        # STEP 1 SUMMARY
        self.logger.info(f'\n{"="*100}')
        self.logger.info(f'STEP 1 SUMMARY: Data Loading and Grouping by Unique Crossing Sets')
        self.logger.info(f'{"="*100}')
        self.logger.info(f'Total raw entanglements loaded: {Step1_QC_counter}')
        self.logger.info(f'Number of protein IDs: {len(num_raw_ents)}')
        for prot_id, count in sorted(num_raw_ents.items()):
            self.logger.info(f'  - {prot_id}: {count} raw entanglements')
        self.logger.info(f'Unique crossing patterns identified: {len(grouped_entanglement_data)}')
        self.logger.info(f'Disulfide bonds found: {len(CCBonds)}')
        if CCBonds:
            self.logger.info(f'  - Disulfide bond pairs: {CCBonds}')

        ### STEP 2 ################################################################################################################
        ############################################################################################
        # Step 2 Get the minimal loop encompassing each set of unique crossings
        self.logger.info(f'\n# Step 2a')
        for ID_cr, loops in grouped_entanglement_data.items():

            ID = ID_cr[0]

            crossings = np.asarray(list(ID_cr[1:]))

            loop_lengths = [nc[1] - nc[0] for nc in loops]

            minimum_loop_length = min(loop_lengths)

            minimum_loop_length_index = loop_lengths.index(minimum_loop_length)

            minimum_loop_nc_i, minimum_loop_nc_j = loops[minimum_loop_length_index]

            ent_data[ID].append((len(loops), minimum_loop_nc_i, minimum_loop_nc_j, *crossings, ';'.join(['-'.join([str(loop[0]), str(loop[1])]) for loop in loops])))

        # STEP 2a SUMMARY
        self.logger.info(f'{"="*100}')
        self.logger.info(f'STEP 2A SUMMARY: Minimal Loop Identification for Unique Crossing Sets')
        self.logger.info(f'{"="*100}')
        for ID, ents in ent_data.items():
            Step2a_QC_counter = 0
            for ent_i, ent in enumerate(ents):
                #print(ID, ent_i, ent)
                Step2a_QC_counter += ent[0]
            
            self.logger.info(f'{ID}: {Step2a_QC_counter} raw entanglements grouped into {len(ents)} representative loops')
            for ent_i, ent in enumerate(ents):
                num_loops = ent[0]
                loop_i, loop_j = ent[1], ent[2]
                crossings = [str(c) for c in ent[3:-1]]
                self.logger.info(f'  Representative {ent_i+1}: Loop ({loop_i}, {loop_j}), ' +
                      f'Crossings={crossings if crossings else "none"}, ' +
                      f'Represents {num_loops} raw entanglement(s)')

            ## QC that the number of tracked entanglements after step 2a is still valid
            #print(f'Step2a_QC_counter: {Step2a_QC_counter} should = {num_raw_ents[ID]}')
            if Step2a_QC_counter != num_raw_ents[ID]:
                raise ValueError(f'The number of tracked entaglements after Step 2a {Step2a_QC_counter} != {num_raw_ents[ID]}')

        ############################################################################################
        # Step 2b: 
        self.logger.info(f'{"="*100}')
        self.logger.info(f'STEP 2B: Merging of Entanglements Based on Crossing Proximity')
        self.logger.info(f'{"="*100}')
        merged_ents = []
        for ID, ents in ent_data.items():
            orig_ents = ents.copy()
            comb_ents = itertools.combinations(ents, 2)

            # for each pair of ents
            for each_ent_pair in comb_ents:
                #print(f'\nAnalyzing pair: {each_ent_pair}')
                if each_ent_pair[0] == each_ent_pair[1]:
                    self.logger.info(f'Ents are the same: {each_ent_pair}')
                    continue

                distance_thresholds = []

                ent1 = each_ent_pair[0]
                ent2 = each_ent_pair[1]

                # get crossings from ent pair without chiralities
                cr1 = set([int(ent_cr_1[1:]) for ent_cr_1 in list(ent1[3:-1])])
                cr2 = set([int(ent_cr_2[1:]) for ent_cr_2 in list(ent2[3:-1])])
                #print(cr1, cr2)

                # get all possible pairs of the crossings
                all_cr_pairs = itertools.product(ent1[3:-1], ent2[3:-1])

                # get the distances between all pairs of crossings
                cr_dist_same_chiral = np.abs([int(pr[0][1:]) - int(pr[1][1:]) for pr in all_cr_pairs if pr[0][0] == pr[1][0]])
                #print(cr_dist_same_chiral)
                
                
                # if any of those distances is less than 3 and the number of crossings is not the same and ij in range of kl
                dist_check = np.any(cr_dist_same_chiral <= 3)
                loop_check = self.check_step_ij_kl_range(ent1, ent2)
                diff_cross_size_check = len(cr1) != len(cr2)
                #print(dist_check, loop_check, diff_cross_size_check)

                if np.any(cr_dist_same_chiral <= 3) and len(cr1) != len(cr2) and self.check_step_ij_kl_range(ent1, ent2):
                    #print(f'\nAnalyzing pair: {each_ent_pair}')
                    #print(f'step 2b conditions met')

                    minumum_loop_base = min(ent1[1], ent1[2], ent2[1], ent2[2])

                    maximum_loop_base = max(ent1[1], ent1[2], ent2[1], ent2[2])

                    all_crossings = cr1.union(cr2)

                    min_max_loop_base_range = set(range(minumum_loop_base, maximum_loop_base + 1))

                    # if the crossings are not within the min max loop range covering both entanglements
                    if not min_max_loop_base_range.intersection(all_crossings):

                        fewer_cr = min(cr1, cr2, key = len)
                        more_cr = max(cr1, cr2, key = len)


                        # Direct permutations - eliminates combinatorial explosion
                        # Equivalent to: generate all valid injective mappings from fewer_cr to more_cr
                        # This replaces itertools.product(*groupings) + column-uniqueness filter
                        all_pair_groupings = set(
                            tuple(zip(sorted(fewer_cr), perm))
                            for perm in itertools.permutations(more_cr, len(fewer_cr))
                        )


                        for condensed_pair in all_pair_groupings:

                            if isinstance(condensed_pair[0], int):

                                # when dealing with ent with one crossing 

                                condensed_pair = [condensed_pair]
                            
                            dist = np.sqrt(sum([(each_ele[0] - each_ele[1]) ** 2 for each_ele in condensed_pair]))

                            distance_thresholds.append(dist)

                        # all_pair_groupings and distance thresholds have the same size
                        if min(distance_thresholds) <= 20:

                            min_ent = min(ent1, ent2, key = len)
                            max_ent = max(ent1, ent2, key = len)
                            if min_ent == max_ent:
                                self.logger.warning(f'WARNING: Ents are the same. Setting min_ent = ent1 and max_ent = ent2')
                                min_ent = ent1
                                max_ent = ent2
                            #print(f'ent1: {ent1} | ent2: {ent2}')
                            #print(f'min_ent: {min_ent}')
                            #print(f'max_ent: {max_ent}')

                            if min_ent in ents and len(ents) > 1:
                            #if min_ent in ents and max_ent in ents and len(ents) > 1:
                                #min_ent_num_ncs = min_ent[0]
                    
                                #print(f'Removing: min_ent {min_ent} at index {ents.index(min_ent)}')
                                del ents[ents.index(min_ent)]
                                #del ents[ents.index(max_ent)]

                                if max_ent == min_ent:
                                    raise ValueError(f'WARNING: Ents are the same\n{min_ent} == {max_ent}')
                                else:
                                    merged_ents += [(max_ent, min_ent)]

        
        #print(f'\nStep 2b results')
        # results foor the end of step 2
        for ID, ents in ent_data.items():
            ent_dict = {ent_idx:[ent] for ent_idx, ent in enumerate(ents)}
            #for ent_idx, ent in ent_dict.items():
            #    print(ent_idx, ent)

            ### Update entanglement list with those that got merged
            self.logger.info(f'\n  Processing {ID}: Analyzing {len(ents)} representative entanglements for merging...')
            self.logger.info(f'  Before merge: {len(ents)} representatives')
            merge_count = 0
            while len(merged_ents) != 0:
                pre_num_merged = len(merged_ents)
                #print(f'# merged_ents: {pre_num_merged}')
                for m_ent in merged_ents:
                    #print(f'\n{m_ent}')
                    for ent_idx, ent in ent_dict.copy().items():
                        #print(ent_idx, ent)
                        if m_ent[0] in ent:
                            #print(f'FOUND MATCH for kept ent {ent_idx}')
                            ent_dict[ent_idx] += [m_ent[1]]
                            merged_ents.remove(m_ent)
                            merge_count += 1

                #print(f'# merged_ents: {len(merged_ents)}')

                # QC to ensure you dont enter an infinite loop
                if len(merged_ents) == pre_num_merged:
                    raise ValueError('Failed to find a match this cycle and entering infi loop')
                    
            self.logger.info(f'  After merge: {merge_count} merges completed')
            self.logger.info(f'  Reformatting entanglement data...')
            updated_ents = []
            for ent_idx, ent in ent_dict.items():
                #print(ent_idx, ent)
                if len(ent) > 1:
                    num_loops = np.sum([e[0] for e in ent])
                    NCs = ';'.join([e[-1] for e in ent])
                    #print(ent, num_loops, NCs)
                    ent = (num_loops, *ent[0][1:-1], NCs)
                    updated_ents += [ent]
                else:
                    updated_ents += ent

            #print(f'Results after adding those that got merged to each representative entanglement')
            Step2b_QC_counter = 0
            for uent in updated_ents:
                #print(uent)
                Step2b_QC_counter += uent[0]

            self.logger.info(f'  Result: {len(updated_ents)} representative entanglements after merging (tracking {Step2b_QC_counter} total raw)')
            self.logger.debug(updated_ents)
            ## QC that the number of tracked entanglements after step 2a is still valid
            if Step2b_QC_counter != num_raw_ents[ID]:
                raise ValueError(f'The number of tracked entaglements after Step 2b {Step2b_QC_counter} != {num_raw_ents[ID]}')        

            ent_data[ID] = updated_ents

        ### STEP 3 ################################################################################################################
        # Step 3
        self.logger.info(f'{"="*100}')
        self.logger.info(f'STEP 3: Removing Duplicate Entanglements (Same Crossings, Different Loop Sizes)')
        self.logger.info(f'{"="*100}')
        for ID, processed_ents in ent_data.items():
            self.logger.info(f'  {ID}: Checking {len(processed_ents)} entanglements for duplicates...')

            comb_processed_ents = itertools.combinations(processed_ents, 2)

            keep_track_of_larger_proc_ent = []
            keep_track_of_shorter_proc_ent = []
            removal_count = 0

            for each_processed_ent_pair in comb_processed_ents:
                #print(f'\npair of ents: {each_processed_ent_pair}')

                proc_ent1 = each_processed_ent_pair[0]
                proc_ent2 = each_processed_ent_pair[1]

                proc_ent1_ijr = proc_ent1[1:-1]
                proc_ent2_ijr = proc_ent2[1:-1]
                #print(proc_ent1_ijr, proc_ent2_ijr)

                if proc_ent1_ijr not in keep_track_of_larger_proc_ent and proc_ent2_ijr not in keep_track_of_larger_proc_ent:

                    # without chiralites
                    proc_cr1 = np.asarray([int(ent_cr_1[1:]) for ent_cr_1 in list(proc_ent1[3:-1])])
                    proc_cr2 = np.asarray([int(ent_cr_2[1:]) for ent_cr_2 in list(proc_ent2[3:-1])])

                    if len(proc_ent1[3:-1]) == len(proc_ent2[3:-1]):

                        chirality1 = [chir1[0] for chir1 in proc_ent1[3:-1]]
                        chirality2 = [chir2[0] for chir2 in proc_ent2[3:-1]]

                        if chirality1 == chirality2 and self.check_step_ij_kl_range(proc_ent1, proc_ent2) and np.all(np.abs(proc_cr1 - proc_cr2) <= 20):
                            #print(proc_ent1, proc_ent2)

                            loop1_length = proc_ent1[2] - proc_ent1[1]
                            loop2_length = proc_ent2[2] - proc_ent2[1]
                            #print(loop1_length, loop2_length)

                            check = [loop1_length, loop2_length]

                            maximum_loop_length = max(loop1_length, loop2_length)
                            minimum_loop_length = min(loop1_length, loop2_length)

                            if maximum_loop_length == minimum_loop_length:
                                longer_loop_ent = proc_ent1
                                shorter_loop_ent = proc_ent2
                            else:
                                longer_loop_ent = each_processed_ent_pair[check.index(maximum_loop_length)]
                                shorter_loop_ent = each_processed_ent_pair[check.index(minimum_loop_length)]
                            longer_loop_ent_ijr = longer_loop_ent[1:-1]
                            shorter_loop_ent_ijr = shorter_loop_ent[1:-1]

                            if len(processed_ents) > 1:

                                for long_proc_ent_index, long_proc_ent in enumerate(processed_ents):
                                    long_proc_ent_ijr = long_proc_ent[1:-1]
                                    #print(long_proc_ent_index, long_proc_ent, long_proc_ent_ijr)
                                    if long_proc_ent_ijr == longer_loop_ent_ijr:
                                        break
                                del processed_ents[long_proc_ent_index]
                                # remove the one with larger loop

                                # find the shorter loop and remove it
                                for short_proc_ent_index, short_proc_ent in enumerate(processed_ents):
                                    short_proc_ent_ijr = short_proc_ent[1:-1]
                                    if short_proc_ent_ijr == shorter_loop_ent_ijr:
                                        break
                                del processed_ents[short_proc_ent_index]

                                updated_ent = (short_proc_ent[0] + long_proc_ent[0],  *short_proc_ent[1:-1], short_proc_ent[-1] + ';' + long_proc_ent[-1])
                                processed_ents += [updated_ent]

                            keep_track_of_larger_proc_ent.append(longer_loop_ent_ijr)
                            keep_track_of_shorter_proc_ent.append(shorter_loop_ent_ijr)
        
        # STEP 3 FINAL SUMMARY
        self.logger.info(f'\nSTEP 3 RESULTS:')
        for ID, ents in ent_data.items():
            Step3_QC_counter = 0
            for ent in ents:
                #print(ent)
                Step3_QC_counter += ent[0]
            
            self.logger.info(f'  {ID}: {len(ents)} representative entanglements remaining (tracking {Step3_QC_counter} raw)')
            # QC to ensure number of raw ents was preserved after step 3
            if Step3_QC_counter != num_raw_ents[ID]:
                raise ValueError(f'The number of tracked entaglements after Step 3 {Step3_QC_counter} != {num_raw_ents[ID]}')

        ### STEP 4 SPATIAL CLUSTERING ################################################################################################################
        # Step 4 prep
        self.logger.info(f'{"="*100}')
        self.logger.info(f'STEP 4 PREP: Grouping Entanglements by Number and Chirality of Crossings')
        self.logger.info(f'{"="*100}')
        for ID, new_ents in ent_data.items():

            for ent in new_ents:

                number_of_crossings = len(ent[3:-1])

                chiralites = [each_cr[0] for each_cr in ent[3:-1]]

                ID_num_chirality_key = f"{ID}_{number_of_crossings}_{chiralites}"

                full_entanglement_data[ID_num_chirality_key].append(ent)
        
        self.logger.info(f'\nGrouping Summary:')
        for group_key in sorted(full_entanglement_data.keys()):
            ents = full_entanglement_data[group_key]
            self.logger.info(f'  {group_key}: {len(ents)} entanglements')

        reset_counter = []

        # Step 4
        self.logger.info(f'{"="*100}')
        self.logger.info(f'STEP 4: Primary Structure Clustering Within Each Group')
        self.logger.info(f'{"="*100}')
        for ID_num_chiral in full_entanglement_data.keys():
            #print(ID_num_chiral)
            #ID = ID_num_chiral.split("_")[0]

            if ID not in reset_counter:

                reset_counter.append(ID)

                split_cluster_counter = 0

            length_key = defaultdict(list)
            loop_dist = defaultdict(list)
            dups = []
            clusters = {} 
            cluster_count = 0

            pairwise_entanglements = list(itertools.combinations(full_entanglement_data[ID_num_chiral], 2))
            
            self.logger.info(f'\n  Group: {ID_num_chiral}')
            self.logger.info(f'    Total entanglements: {len(full_entanglement_data[ID_num_chiral])}')
            self.logger.info(f'    Pairwise comparisons: {len(pairwise_entanglements)}')

            if pairwise_entanglements:

                for i, pairwise_ent in enumerate(pairwise_entanglements):

                    dist = self.loop_distance(pairwise_ent[0], pairwise_ent[1])
                    
                    if dist <= self.cut_off and pairwise_ent[0] not in dups and pairwise_ent[1] not in dups:
                        # 1. pair must be <= self.cut_off
                        # 2. the neighbor cannot be the next key and it cannot be captured by another key

                        loop_dist[pairwise_ent[0]].append(pairwise_ent[1])
                        dups.append(pairwise_ent[1])
                
                key_list = list(loop_dist.keys())

                for key in key_list:

                    length_key[len(loop_dist[key])].append(key)

                # create clusters

                while len(length_key.values()) > 0:

                    max_neighbor = max(length_key.keys())

                    selected_ent = random.choice(length_key[max_neighbor])

                    cluster = copy.deepcopy(loop_dist[selected_ent])
                    cluster.append(selected_ent)

                    clusters[cluster_count] = cluster
                    cluster_count += 1

                    length_key[max_neighbor].remove(selected_ent)

                    if len(length_key[max_neighbor]) == 0:
                        length_key.pop(max_neighbor)
            
            # create single clusters
            if clusters:
                clusters_ijr_values = list(itertools.chain.from_iterable(list(clusters.values())))
            else:
                clusters_ijr_values = []

            full_ent_values = np.asarray(full_entanglement_data[ID_num_chiral], dtype=object)

            difference_ent = np.zeros(len(full_ent_values), dtype=bool)

            for k, ijr in enumerate(full_ent_values):

                if tuple(ijr) in clusters_ijr_values:
                    difference_ent[k] = True
                else:
                    difference_ent[k] = False

            i = np.unique(np.where(difference_ent == False)[0])

            next_cluster_count = cluster_count

            for single_cluster in full_ent_values[i]:
                
                single_cluster_list = []
                single_cluster_list.append(tuple(single_cluster))

                clusters[next_cluster_count] = single_cluster_list

                next_cluster_count += 1

            # pick representative entanglement per cluster
            self.logger.info(f'    Primary structure clusters formed: {len(clusters)}')
            for counter, ijr_values in clusters.items():
                #print(f'\nCluster {counter} {ijr_values}')

                # clusters contain many entanglements
                if len(ijr_values) > 1:

                    ijr = np.asarray(ijr_values)
                    #print(f'cluster ijr:\n{ijr}')

                    cr_values = np.asarray([[int(r_value[0][1:])] for r_value in ijr[:, 3:-1]])
                    #print(f'cr_values: {cr_values}')

                    median_cr = compute_geometric_median(cr_values).median
                    #print(f'median_cr: {median_cr}')

                    distances = cdist(cr_values, [median_cr])

                    minimum_distances_i = np.where(distances == min(distances))[0]
                    #print(f'minimum_distances_i: {minimum_distances_i}')

                    possible_cand = ijr[minimum_distances_i]
                    #print(f'possible_cand:\n{possible_cand}')

                    loop_lengths = np.abs(possible_cand[:, 1].astype(int) - possible_cand[:, 2].astype(int))
                    #print(f'loop_lengths: {loop_lengths}')

                    smallest_loop_length = min(loop_lengths)
                    #print(f'smallest_loop_length: {smallest_loop_length}')

                    num_nc_in_cluster = np.sum([int(n[0]) for n in ijr_values])
                    #print(f'num_nc_in_cluster: {num_nc_in_cluster}')

                    loops = ';'.join([n[-1] for n in ijr_values])
                    #print(f'loops: {loops}')
                    
                    #rep_entanglement = possible_cand[random.choice(np.where(smallest_loop_length == loop_lengths)[0])]
                    rep_entanglement = possible_cand[random.choice(np.where(smallest_loop_length == loop_lengths))[0]]
                    rep_entanglement = [str(num_nc_in_cluster), *rep_entanglement[1:-1], loops]
                    #rep_ID_ent[f"{ID}_{split_cluster_counter}"].append(rep_entanglement)
                    rep_ID_ent[(ID, split_cluster_counter)].append(rep_entanglement)

                # clusters with a single entnalgement
                else:
                    #rep_ID_ent[f"{ID}_{split_cluster_counter}"].append(ijr_values[0])
                    rep_ID_ent[(ID, split_cluster_counter)].append(ijr_values[0])
                    if counter == list(clusters.keys())[-1]:  # Print only for last cluster to avoid clutter
                        num_single = sum(1 for c_vals in clusters.values() if len(c_vals) == 1)
                        self.logger.info(f'    Single-entanglement clusters: {num_single}')
                
                split_cluster_counter += 1
        
        ## QC Step 4 results
        self.logger.info(f'\n{"="*100}')
        self.logger.info(f'STEP 4 FINAL RESULTS: Primary Structure Clustering Summary')
        self.logger.info(f'{"="*100}')
        num_raw_ents_FINAL = {}
        for ID_counter, ijrs in rep_ID_ent.items():
            #print(ID_counter, ijrs)
            ID, counter = ID_counter
            #print(ID_counter, ID, counter, ijrs)

            if ID not in num_raw_ents_FINAL:
                num_raw_ents_FINAL[ID] = 0

            for ijr in ijrs:
                num_nc = int(ijr[0])
                num_raw_ents_FINAL[ID] += num_nc

        ## check the final tracking of raw ents
        for ID, count in num_raw_ents.items():
            final_count = num_raw_ents_FINAL[ID]
            num_clusters = len([ijrs for (c_id, _), ijrs in rep_ID_ent.items() if c_id == ID])
            self.logger.info(f'{ID}: {count} raw → {final_count} raw in {num_clusters} final clusters')
            if count != final_count:
                raise ValueError(f'The FINAL # of raw ents {final_count} != the starting {count} for ID {ID}')

        ### STEP 5 OUTPUT FILE ################################################################################################################
        # Step 5
        self.logger.info(f'{"="*100}')
        self.logger.info(f'STEP 5: Writing Output File')
        self.logger.info(f'{"="*100}')

        ## set up the outdir for this calculation
        #outdir = f"{os.getcwd()}/{outdir}"
        if not os.path.isdir(outdir):
            os.mkdir(f"{outdir}") 
            self.logger.info(f"Creating directory: {outdir}")

        outfilepath = os.path.join(f'{outdir}', f'{outfile}')

        with open(outfilepath, "w") as f:

            f.write(f'ID|chain|i|j|crossingsN|crossingsC|gn|gc|GLNn|GLNc|TLNn|TLNc|num_contacts|contacts|CCBond\n')
            for ID_counter, ijrs in rep_ID_ent.items():

                ID, counter = ID_counter
                chain = chain_info.get(ID, 'A')  # Get chain for this ID, default to 'A'

                for ijr in ijrs:

                    new_ijr = (int(ijr[1]), int(ijr[2]), *list(ijr[3:-1]))

                    num_nc = int(ijr[0])

                    gn, gc, GLNn, GLNc, TLNn, TLNc, crossingsN_stored, crossingsC_stored = entanglement_partial_g_data[new_ijr]
                    gn = float(gn)
                    gc = float(gc)
                    GLNn = int(GLNn)
                    GLNc = int(GLNc)
                    # Handle NaN/empty TLN values: convert to int only if not NaN
                    TLNn = np.nan if pd.isna(TLNn) else int(TLNn)
                    TLNc = np.nan if pd.isna(TLNc) else int(TLNc)
                    
                    # Separate crossings into N and C terminal
                    all_crossings = list(ijr[3:-1])
                    crossingsN = []
                    crossingsC = []
                    i_val = int(ijr[1])
                    j_val = int(ijr[2])
                    for cross in all_crossings:
                        cross_resid = int(cross[1:])
                        if cross_resid < i_val:
                            crossingsN.append(cross)
                        elif cross_resid > j_val:
                            crossingsC.append(cross)
                    
                    crossingsN_str = ','.join(crossingsN) if crossingsN else ''
                    crossingsC_str = ','.join(crossingsC) if crossingsC else ''

                    ## check for disulfide bonds
                    CCBond_flag = False
                    for CCBond in CCBonds:
                        check1 = f'{CCBond[0]}-{CCBond[1]}' 
                        check2 = f'{CCBond[1]}-{CCBond[0]}'
                        if check1 in ijr[-1] or check2 in ijr[-1]:
                            CCBond_flag = True

                    line = f"{ID}|{chain}|{int(ijr[1])}|{int(ijr[2])}|{crossingsN_str}|{crossingsC_str}|{gn:.5f}|{gc:.5f}|{GLNn}|{GLNc}|{TLNn}|{TLNc}|{num_nc}|{ijr[-1]}|{CCBond_flag}"
                    #print(line)
                    f.write(f"{line}\n")
        self.logger.info(f'SAVED: {outfilepath}')
        outdf = pd.read_csv(outfilepath, sep='|')
        
        # FINAL CLUSTERING SUMMARY
        self.logger.info(f'\n{"="*100}')
        self.logger.info(f'CLUSTERING COMPLETE: Final Summary')
        self.logger.info(f'{"="*100}')
        self.logger.info(f'Total raw entanglements processed: {sum(num_raw_ents.values())}')
        self.logger.info(f'Total final representative entanglements: {len(outdf)}')
        # print(f'Compression ratio: {sum(num_raw_ents.values())/len(outdf):.2f}x (raw → final)')
        self.logger.info(f'Clustering by organism: {self.organism}')
        self.logger.info(f'Spatial distance cutoff: {self.cut_off}')
        self.logger.info(f'Output file: {outfilepath}')
        self.logger.info(f'{"="*100}\n')
        
        return {'outfile':outfilepath, 'ent_result':outdf}
    ##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################


##########################################################################################################################################################
##########################################################################################################################################################
class ClusterNonNativeEntanglements:
    """
    Class to calculate Non-native entanglements given either a file path to an entanglement file or an entanglement object
    """

    ##########################################################################################################################################################
    def __init__(self, trajnum2pklfile_path:str, traj_dir_prefix:str='./', outdir:str='./ClusterNonNativeEntanglements/', log_level:int=logging.INFO, logdir:str=None, nproc:int=1) -> None:
        """
        Constructor for GaussianEntanglement class.

        Parameters
        ----------
        """
        self.classify_key = ['topoly_linking_number']
        self.cluster_method = ['average', 'average', 'average']
        # cluster_dist_cutoff = [20, 1.0, 0.6] # Allow contamination
        self.cluster_dist_cutoff = [20, 1.0, 0.1] # No contamination
        self.memory_cutoff = 6.4e10 # 64 Gb
        self.max_plot_samples = 1000

        # matplotlib.rcParams['mathtext.fontset'] = 'stix'
        # matplotlib.rcParams['font.sans-serif'] = ['Arial']
        matplotlib.rcParams['axes.labelsize'] = 'small'
        matplotlib.rcParams['axes.linewidth'] = 1
        matplotlib.rcParams['lines.markersize'] = 4
        matplotlib.rcParams['xtick.major.width'] = 1
        matplotlib.rcParams['ytick.major.width'] = 1
        matplotlib.rcParams['xtick.labelsize'] = 'x-small'
        matplotlib.rcParams['ytick.labelsize'] = 'x-small'
        matplotlib.rcParams['legend.fontsize'] = 'x-small'
        matplotlib.rcParams['figure.dpi'] = 600

        self.nproc = max(1, int(nproc))
        self.logger = setup_logger('ClusterNonNativeEntanglements', outdir=logdir if logdir is not None else outdir, log_level=log_level)

        self.traj_dir_prefix = traj_dir_prefix
        
        ## Load the dataframe that maps trajectory numbers to pkl file paths
        self.trajnum2pklfile_path = trajnum2pklfile_path
        self.trajnum2pklfile = pd.read_csv(self.trajnum2pklfile_path)
        
        ## Extract pkl file paths from the manifest (source of truth)
        if 'pklfile' not in self.trajnum2pklfile.columns:
            self.logger.error('Error: trajnum2pklfile CSV must contain a "pklfile" column')
            sys.exit()
        
        self.ent_data_file_list = self.trajnum2pklfile['pklfile'].tolist()
        
        # Verify all pkl files exist
        missing_files = [f for f in self.ent_data_file_list if not os.path.isfile(f)]
        if missing_files:
            self.logger.error(f'Error: {len(missing_files)} pkl files not found:')
            for f in missing_files:
                self.logger.error(f'  {f}')
            sys.exit()
        
        self.logger.info(f'FOUND {len(self.ent_data_file_list)} .pkl files to cluster from manifest')

        ## Set up the outdir for this calculation
        self.outdir = outdir
        if not os.path.isdir(self.outdir):
            os.mkdir(f"{self.outdir}") 
            self.logger.info(f"Creating directory: {self.outdir}")
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def save_pickle(self, filename, mode, data, protocol=4):
        with open(filename, mode) as fh:
            pickle.dump(data, fh, protocol=protocol)
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def load_pickle(self, filename, start_frame=None, end_frame=None):
        """Load pickle file and optionally filter frames.
        
        Parameters
        ----------
        filename : str
            Path to pickle file
        start_frame : int, optional
            Minimum frame index to keep (inclusive)
        end_frame : int, optional
            Maximum frame index to keep (inclusive)
        
        Returns
        -------
        dict
            Dictionary with frame keys and 'ref' key. Filtered to frame range if specified.
        
        Note: Large unfiltered dictionaries are explicitly deleted to ensure timely 
        garbage collection, especially important when multiple threads load in parallel.
        """
        import gc
        data_dict = {}
        self.logger.debug(f'Loading pickle file: {filename}')
        with open(filename, 'rb') as fr:
            try:
                while True:
                    chunk = pickle.load(fr)
                    self.logger.debug(f"Loaded chunk with size: {len(chunk)}")
                    if start_frame is not None or end_frame is not None:
                        self.logger.debug(f"Filtering chunk for frames between {start_frame} and {end_frame}")
                        filtered_chunk = {}

                        for k, v in chunk.items():
                            if k == "ref":
                                filtered_chunk[k] = v

                            elif ((start_frame is None or k >= start_frame) and (end_frame is None or k <= end_frame)):
                                filtered_chunk[k] = v

                        del chunk  # CRITICAL: Explicitly free large unfiltered dict before update
                        chunk = filtered_chunk

                    self.logger.debug(f"Loaded filtered chunk with size: {len(chunk)}")
                    data_dict.update(chunk)
                    del chunk  # Free memory after update

            except EOFError:
                pass
        
        self.logger.debug(f'Total frames loaded: {len(data_dict) - 1}')  # Exclude 'ref' key from frame count
        
        # Force garbage collection to prevent memory accumulation in parallel loads
        gc.collect()
        return data_dict
    ##########################################################################################################################################################
    
    ##########################################################################################################################################################
    def extract_traj_number(self, f):
        f_match = self.trajnum2pklfile[self.trajnum2pklfile['pklfile'] == f]
        if f_match.empty:
            self.logger.error(f'Error: {f} not found in {self.trajnum2pklfile_path}')
            raise ValueError(f'Error: {f} not found in {self.trajnum2pklfile_path}')
        match = f_match['trajnum'].values[0]
        return match
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def pdist_loop_overlap(self, data_array_1, data_array_2):
        M_1 = np.repeat(data_array_1.reshape((data_array_1.shape[0], 1, data_array_1.shape[1])), data_array_2.shape[0], axis=1)
        M_2 = np.repeat(data_array_2.reshape((1, data_array_2.shape[0], data_array_2.shape[1])), data_array_1.shape[0], axis=0)
        M = np.concatenate((M_1, M_2), axis=-1)
        del M_1, M_2
        dist_M = (np.max(M, axis=-1) - np.min(M, axis=-1) + 1) / (M[:,:,1]-M[:,:,0]+M[:,:,3]-M[:,:,2])
        # Make distance between inclusive loops to be minimum
        dist_M[(M[:,:,2]-M[:,:,0])*(M[:,:,1]-M[:,:,3]) >= 0] = 0.5 
        return dist_M
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def pdist_thread_deviation(self, data_array_1, data_array_2):
        M_1 = np.repeat(data_array_1.reshape((data_array_1.shape[0], 1, data_array_1.shape[1])), data_array_2.shape[0], axis=1)
        M_2 = np.repeat(data_array_2.reshape((1, data_array_2.shape[0], data_array_2.shape[1])), data_array_1.shape[0], axis=0)
        M = np.concatenate((M_1, M_2), axis=-1)
        del M_1, M_2
        dist_M = np.abs((M[:,:,2:4]-M[:,:,0:2]))
        dist_M[M[:,:,2:4]*M[:,:,0:2] < 0] = 10 # Make distance between no crossing and crossings to be small
        del M
        dist = np.max(dist_M, axis=-1)
        return dist
    ##########################################################################################################################################################
        
    ##########################################################################################################################################################
    def pdist_cross_contamination(self, data_array_1, data_array_2):
        # data looks like [nc_1, nc_2, cross_N1, cross_N2, ..., cross_C1, cross_C2, ...]
        M_1 = np.repeat(data_array_1.reshape((data_array_1.shape[0], 1, data_array_1.shape[1])), data_array_2.shape[0], axis=1)
        M_2 = np.repeat(data_array_2.reshape((1, data_array_2.shape[0], data_array_2.shape[1])), data_array_1.shape[0], axis=0)
        M = np.concatenate((M_1, M_2), axis=-1)
        del M_1, M_2

        # Distance for cross_2 contaminate loop_1
        idx_array_1 = np.zeros((data_array_2.shape[1]-2,2), dtype=int)
        idx_array_1[:,0] = 1
        idx_array_1[:,1] = np.arange(data_array_1.shape[1]+2, data_array_1.shape[1]+data_array_2.shape[1], dtype=int)
        idx_array_2 = np.zeros((data_array_2.shape[1]-2,2), dtype=int)
        idx_array_2[:,0] = np.arange(data_array_1.shape[1]+2, data_array_1.shape[1]+data_array_2.shape[1], dtype=int)
        idx_array_2[:,1] = 0
        L = (M[:,:,1]-M[:,:,0]).reshape((M.shape[0], M.shape[1], 1))

        dist_M_1 = np.min(M[:,:,idx_array_1]-M[:,:,idx_array_2], axis=-1) / L
        dist_M_1[dist_M_1 <= 0] = 0
        dist_M_1[dist_M_1 >= 1] = 0
        
        # Distance for cross_1 contaminate loop_2
        idx_array_1 = np.zeros((data_array_1.shape[1]-2,2), dtype=int)
        idx_array_1[:,0] = data_array_1.shape[1]+1
        idx_array_1[:,1] = np.arange(2, data_array_1.shape[1], dtype=int)
        idx_array_2 = np.zeros((data_array_1.shape[1]-2,2), dtype=int)
        idx_array_2[:,0] = np.arange(2, data_array_1.shape[1], dtype=int)
        idx_array_2[:,1] = data_array_1.shape[1]
        L = (M[:,:,data_array_1.shape[1]+1]-M[:,:,data_array_1.shape[1]]).reshape((M.shape[0], M.shape[1], 1))

        dist_M_2 = np.min(M[:,:,idx_array_1]-M[:,:,idx_array_2], axis=-1) / L
        dist_M_2[dist_M_2 <= 0] = 0
        dist_M_2[dist_M_2 >= 1] = 0
        
        del M
        dist_M = np.max(np.concatenate((dist_M_1, dist_M_2), axis=-1), axis=-1)
        return dist_M
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def agglomerative_clustering(self, dist, cluster_method, cluster_dist_cutoff, num_perm):
        min_SSDIFN = np.inf
        best_Z = None
        best_perm_idx_list = None
        pdist = squareform(dist, checks=False)
        if np.sum(pdist**2) == 0:
            best_Z = linkage(pdist, method=cluster_method)
            best_perm_idx_list = np.arange(dist.shape[0])
        else:
            # permuCLUSTER
            for idx_perm in range(np.max([1, num_perm])):
                perm_idx_list = np.random.permutation(np.arange(dist.shape[0]))
                pdist = squareform(dist[perm_idx_list,:][:,perm_idx_list], checks=False)
                Z = linkage(pdist, method=cluster_method)
                cdist = cophenet(Z)
                SSDIFN = np.sum((pdist - cdist)**2)/np.sum(pdist**2)
                if SSDIFN < min_SSDIFN:
                    min_SSDIFN = SSDIFN
                    best_Z = Z
                    best_perm_idx_list = perm_idx_list
        cluster_id_list = fcluster(best_Z, cluster_dist_cutoff, criterion='distance')
        backmap_list = np.zeros(len(best_perm_idx_list), dtype=int)
        for i, j in enumerate(best_perm_idx_list):
            backmap_list[j] = i
        cluster_id_list = cluster_id_list[backmap_list]
        return cluster_id_list
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def do_clustering(self, map_list, chg_ent_fingerprint_list, key, pdist_fun, cluster_method, cluster_dist_cutoff, num_perm=100):
        data = []
        # Get max number of crossings
        max_n_cross = 1
        if 'cross_contamination' in key:
            for map_idx in map_list:
                fingerprint = chg_ent_fingerprint_list[map_idx[0]][map_idx[1]][tuple(map_idx[2:4])]
                cr = fingerprint['crossing_resid']
                for ci, c in enumerate(cr):
                    if len(c) > max_n_cross:
                        max_n_cross = len(c)
        # Prepare clustering data for distance calculation
        for map_idx in map_list:
            fingerprint = chg_ent_fingerprint_list[map_idx[0]][map_idx[1]][tuple(map_idx[2:4])]
            if 'crossing_resid' in key:
                ter_idx = int(key.split('_')[-1])
                mc = []
                cr = fingerprint['crossing_resid']
                ref_cr = fingerprint['ref_crossing_resid']
                for c in [ref_cr[ter_idx], cr[ter_idx]]:
                    if len(c) == 0:
                        mc.append(-1)
                    else:
                        mc.append(np.median(c))
                data.append(mc)
            elif 'native_contact' in key:
                nc = fingerprint['native_contact']
                data.append(nc)
            elif 'cross_contamination' in key:
                nc = fingerprint['native_contact']
                mc = []
                cr = fingerprint['crossing_resid']
                for ci, c in enumerate(cr):
                    for cii in range(max_n_cross):
                        if cii >= len(c):
                            mc.append(-1)
                        else:
                            mc.append(c[cii])
                data.append(nc + mc)
            else:
                self.logger.error('Error: Unknown key specified for do_clustering(), %s'%(key))
                sys.exit()
        data = np.array(data)
        # Reduce data size (saving memory usage) by combining duplicated data points
        reduced_data = np.unique(data, axis=0)
        reduced_data_map = np.array([np.all(data == d, axis=1).nonzero()[0].tolist() for d in reduced_data], dtype=object)
        
        # If all data are the same, group them into a single cluster and return
        if len(reduced_data) == 1:
            cluster_data = [map_list]
            return cluster_data
        
        # Chunk data to reduce memory usage if the expanded matrix occupy >= memory_cutoff
        if reduced_data.nbytes ** 2 >= self.memory_cutoff:
            n_chunk = int(np.ceil(reduced_data.nbytes / np.sqrt(self.memory_cutoff/2)))
            len_chunk = int(np.ceil(len(reduced_data) / n_chunk))
            dist = np.zeros((len(reduced_data), len(reduced_data)))
            for i in range(n_chunk):
                i_1 = i*len_chunk
                i_2 = np.min([(i+1)*len_chunk,len(reduced_data)])
                for j in range(n_chunk):
                    j_1 = j*len_chunk
                    j_2 = np.min([(j+1)*len_chunk,len(reduced_data)])
                    dist[i_1:i_2, j_1:j_2] = pdist_fun(reduced_data[i_1:i_2], reduced_data[j_1:j_2])
        else:
            dist = pdist_fun(reduced_data, reduced_data)
        
        if 'cross_contamination' in key:
            # Do divisive clustering
            cluster_idx_mapping = [list(np.arange(dist.shape[0]))]
            while True:
                cluster_1 = cluster_idx_mapping[-1]
                cluster_2 = []
                cluster_0 = copy.deepcopy(cluster_1)
                for i in range(len(cluster_0)-1):
                    rm_idx_list = np.where(dist[cluster_0[i],cluster_0[i+1:]] >= cluster_dist_cutoff)[0]
                    if len(rm_idx_list) > 0:
                        cluster_1.remove(cluster_0[i])
                        cluster_2.append(cluster_0[i])
                if len(cluster_2) > 0:
                    cluster_idx_mapping.append(cluster_2)
                else:
                    break
            # Do agglomerative clustering
            if cluster_method == 'single':
                dist_fun = np.min
            elif cluster_method == 'complete':
                dist_fun = np.max
            elif cluster_method == 'average':
                dist_fun = np.mean
            else:
                dist_fun = np.mean
            if len(cluster_idx_mapping) > 1:
                dist_0 = np.zeros((len(cluster_idx_mapping),len(cluster_idx_mapping)))
                for i in range(len(dist_0)-1):
                    for j in range(i+1, len(dist_0)):
                        dist_0[i,j] = dist_fun(dist[cluster_idx_mapping[i],:][:,cluster_idx_mapping[j]])
                cluster_id_list_0 = self.agglomerative_clustering(dist_0, cluster_method, cluster_dist_cutoff, num_perm)
            else:
                cluster_id_list_0 = np.array([1], dtype=int)
            cluster_id_list = np.zeros(dist.shape[0], dtype=int)
            for cluster_idx, mapping in enumerate(cluster_idx_mapping):
                cluster_id_list[mapping] = cluster_id_list_0[cluster_idx]
        else:
            # Do agglomerative clustering
            cluster_id_list = self.agglomerative_clustering(dist, cluster_method, cluster_dist_cutoff, num_perm)
        n_cluster = np.max(cluster_id_list)
        
        # Back-Mapping indices 
        cluster_data = []
        for cluster_id in range(n_cluster):
            idx = np.where(cluster_id_list == cluster_id+1)[0]
            idx_list = reduced_data_map[idx].tolist()
            idx_list_1 = []
            for i in idx_list:
                idx_list_1 += i
            idx_list_1 = sorted(idx_list_1)
            cluster_data.append(np.array(map_list)[idx_list_1].tolist())
        return cluster_data
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def cluster_chg_ent(self, chg_ent_keyword_dict, chg_ent_fingerprint_list, cluster_method=['average', 'average', 'average'], cluster_dist_cutoff=[20, 1.0, 0.1]):
        cluster_data_keys = sorted(list(chg_ent_keyword_dict.keys()))
        cluster_data = {key: [] for key in cluster_data_keys}
        cluster_tree = {key: [] for key in cluster_data_keys}

        def _process_key(key):
            """Run the 4-level hierarchical clustering pipeline for one keyword.

            All four ``do_clustering`` calls are fully independent across keywords,
            so this function is safe to run in a ThreadPoolExecutor.  numpy/scipy
            release the GIL during heavy array operations, giving real parallelism.

            Note: when ``self.nproc > 1``, ``agglomerative_clustering`` uses
            ``np.random.permutation`` which draws from the shared global numpy
            random state.  Results are therefore non-deterministic across runs with
            multiple workers, but clustering quality is unaffected.
            """
            map_list = chg_ent_keyword_dict[key]
            local_cluster_data = []
            backtrace_idx_list = []
            idx_1 = 0
            idx_2 = 0
            idx_3 = 0

            # First clustering based on N-ter crossing residues
            N_cr_cluster_data = self.do_clustering(
                map_list, chg_ent_fingerprint_list, 'crossing_resid_0',
                self.pdist_thread_deviation, cluster_method[0], cluster_dist_cutoff[0])

            # Second clustering based on C-ter crossing residues
            for map_list_N_cr in N_cr_cluster_data:
                C_cr_cluster_data = self.do_clustering(
                    map_list_N_cr, chg_ent_fingerprint_list, 'crossing_resid_1',
                    self.pdist_thread_deviation, cluster_method[0], cluster_dist_cutoff[0])

                # Third clustering based on loop
                for map_list_C_cr in C_cr_cluster_data:
                    nc_cluster_data = self.do_clustering(
                        map_list_C_cr, chg_ent_fingerprint_list, 'native_contact',
                        self.pdist_loop_overlap, cluster_method[1], cluster_dist_cutoff[1])

                    # Fourth clustering based on cross contamination
                    for map_list_nc in nc_cluster_data:
                        final_cluster_data = self.do_clustering(
                            map_list_nc, chg_ent_fingerprint_list, 'cross_contamination',
                            self.pdist_cross_contamination, cluster_method[2], cluster_dist_cutoff[2])
                        for final_cluster in final_cluster_data:
                            local_cluster_data.append(final_cluster)
                            backtrace_idx_list.append([idx_1, idx_2, idx_3])
                        idx_3 += 1
                    idx_2 += 1
                idx_1 += 1

            backtrace_idx_list = np.array(backtrace_idx_list, dtype=int)
            local_tree = [
                [np.where(backtrace_idx_list[:, i] == j)[0].tolist()
                 for j in range(backtrace_idx_list[:, i].max() + 1)]
                for i in range(backtrace_idx_list.shape[1])
            ]
            n_cluster = len(local_cluster_data)
            self.logger.info('Found %d cluster(s) for %s' % (n_cluster, key))
            return key, local_cluster_data, local_tree

        n_workers = min(self.nproc, len(cluster_data_keys)) if cluster_data_keys else 1
        self.logger.info(
            f'Clustering {len(cluster_data_keys)} keyword(s) '
            f'(nproc={n_workers})...')
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            for key, cd, ct in executor.map(_process_key, cluster_data_keys):
                cluster_data[key] = cd
                cluster_tree[key] = ct

        return (cluster_data, cluster_tree)
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def find_representative_entanglement(self, cluster_data, ent_cluster_idx_map):
        # most probable loop midpoint
        rep_ent_list = []
        for [key, idx] in ent_cluster_idx_map:
            cluster = np.array(cluster_data[key][idx])
            loop_midpoint_list = np.mean(cluster[:,2:4], axis=1)
            loop_max = np.max(cluster[:,3])
            loop_min = np.min(cluster[:,2])
            # mode
            bins = np.arange(loop_min, loop_max, 5)
            if len(bins) == 1:
                bins = np.arange(loop_min, loop_max, 1)
            hist, edges = np.histogram(loop_midpoint_list, bins=bins)
            idx = np.argmax(hist)
            min_loop_len = 1e6
            idx0 = np.where(loop_midpoint_list >= edges[idx])[0]
            idx1 = np.where(loop_midpoint_list[idx0] < edges[idx+1])[0]
            for iidx in idx0[idx1]:
                if cluster[iidx,3]-cluster[iidx,2] < min_loop_len:
                    rep_ent = cluster[iidx]
                    min_loop_len = rep_ent[3] - rep_ent[2]
            rep_ent_list.append(rep_ent)
        return rep_ent_list
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def _process_traj_file(self, traj_idx, ent_data_file, start_frame, end_frame):
        """Load and pre-process one trajectory pkl file for clustering.

        Called in parallel (one call per trajectory) by ``cluster()``.

        Returns
        -------
        tuple : (traj_idx, fingerprint_dict, Q_dict, frame_list, traj_file, keyword_entries)
            * fingerprint_dict  – {frame: {nc: fingerprint}}
            * Q_dict            – {frame: Q_value}
            * frame_list        – sorted list of in-range frame indices
            * traj_file         – resolved path to the matching .dcd file
            * keyword_entries   – list of (keyword_str, entry) pairs for building
                                  chg_ent_keyword_dict after all trajectories are loaded
        """
        traj = self.extract_traj_number(ent_data_file)
        self.logger.debug(
            f'Processing {ent_data_file} {traj} ({traj_idx + 1} / {len(self.ent_data_file_list)})...')

        # Load pickle with frame filtering applied
        ent_data = self.load_pickle(ent_data_file, start_frame, end_frame)
        # Frame list is now already filtered, just exclude 'ref'
        frame_list = sorted([frame for frame in ent_data.keys() if frame != 'ref'])
        self.logger.debug(f'frame_list: {frame_list} {len(frame_list)}')

        # Locate the matching trajectory DCD file
        traj_file = os.path.join(self.traj_dir_prefix, f'{traj}_*.dcd')
        traj_file = glob.glob(traj_file)
        if len(traj_file) == 0:
            raise ValueError(f'No trajectory file found for {ent_data_file}.')
        elif len(traj_file) > 1:
            raise ValueError(f'More than 1 trajectory file found for {ent_data_file}.')
        else:
            traj_file = traj_file[0]

        fingerprint_dict = {}
        Q_dict = {}
        keyword_entries = []  # collected as (keyword_str, entry) to merge after parallel load

        for frame in frame_list:
            fingerprint_dict[frame] = {}
            Q_dict[frame] = (
                np.sum(list(ent_data[frame]['G_dict'].values()))
                / len(list(ent_data['ref']['ent_fingerprint'].keys())) / 2
            )
            for nc, fingerprint in ent_data[frame]['chg_ent_fingerprint'].items():
                # Skip if no change of entanglement
                if fingerprint['type'] == ['no change', 'no change']:
                    continue
                fingerprint_dict[frame][nc] = fingerprint
                chg_ent_keyword = fingerprint['code'].copy()
                for ck in self.classify_key:
                    if type(fingerprint[ck]) == list:
                        chg_ent_keyword += fingerprint[ck]
                    else:
                        chg_ent_keyword += [fingerprint[ck]]
                chg_ent_keyword = str(chg_ent_keyword)
                keyword_entries.append((chg_ent_keyword, [traj_idx, frame] + list(nc)))

        # Explicitly free the full deserialized pkl dict — it can be ~10 GB and the
        # semaphore in cluster() is held until this function returns, so freeing it
        # here lets CPython recycle the memory before the next load starts.
        del ent_data

        return traj_idx, fingerprint_dict, Q_dict, frame_list, traj_file, keyword_entries
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def cluster(self, start_frame:int=0, end_frame:int=9999999):

        ## Define the .npz file name 
        npz_data_file = f'cluster_data_{"_".join(self.classify_key)}.npz'
        npz_data_file = os.path.join(self.outdir, npz_data_file)
        self.logger.info(f'Checking for {npz_data_file}')

        if not os.path.exists(npz_data_file):
            # Classify changes of entanglement based on the keyword 
            # "[change_code, classify_key_1_N, classify_key_1_C, classify_key_2_N, classify_key_2_C, ...]"
            self.logger.debug('Reading pkl data and classify changes of entanglement...')
            chg_ent_fingerprint_list = [None] * len(self.ent_data_file_list)
            Q_list = [None] * len(self.ent_data_file_list)
            idx2frame = [None] * len(self.ent_data_file_list)
            idx2trajfile = [None] * len(self.ent_data_file_list)
            dtrajs = [None] * len(self.ent_data_file_list)
            chg_ent_keyword_dict = {}
            chg_ent_keyword_list = []
            # combined_traj = None

            n_workers = min(self.nproc, len(self.ent_data_file_list))
            # Each ~1 GB pkl file inflates to ~10 GB of live Python objects during
            # deserialization.  Use memory_cutoff as a proxy for available RAM to
            # derive a safe concurrency limit: memory_cutoff / 1e10 (10 GB/file).
            n_load_workers = min(n_workers, max(1, int(self.memory_cutoff // 1e10)))
            self.logger.info(
                f'Loading {len(self.ent_data_file_list)} pkl files '
                f'(nproc={n_workers}, concurrent_loads={n_load_workers})...')
            # Semaphore caps how many threads may simultaneously hold a deserialized
            # pkl in memory.  n_workers threads are still available to start new
            # loads as soon as a slot frees, so throughput is not reduced.
            _load_sem = threading.Semaphore(n_load_workers)

            def _throttled_process(ti, f):
                with _load_sem:
                    return self._process_traj_file(ti, f, start_frame, end_frame)

            with ThreadPoolExecutor(max_workers=n_workers) as executor:
                futures = {
                    executor.submit(_throttled_process, ti, f): ti
                    for ti, f in enumerate(self.ent_data_file_list)
                }
                for fut in as_completed(futures):
                    ti, fingerprint_dict, Q_dict, frame_list, traj_file, keyword_entries = fut.result()
                    chg_ent_fingerprint_list[ti] = fingerprint_dict
                    Q_list[ti] = Q_dict
                    idx2frame[ti] = frame_list
                    idx2trajfile[ti] = traj_file
                    dtrajs[ti] = [[] for _ in frame_list]
                    # Merge keyword entries into the shared dicts (sequential, no lock needed)
                    for kw, entry in keyword_entries:
                        if kw not in chg_ent_keyword_list:
                            chg_ent_keyword_dict[kw] = []
                            chg_ent_keyword_list.append(kw)
                        chg_ent_keyword_dict[kw].append(entry)

            self.logger.info('%d data files have been read.' % len(self.ent_data_file_list))
            
            # cluster changes of entanglements found in the trajectories
            self.logger.info('Clustering changes of entanglement for %d keywords...'%(len(chg_ent_keyword_list)))
            ent_cluster_data, ent_cluster_tree = self.cluster_chg_ent(chg_ent_keyword_dict, chg_ent_fingerprint_list, cluster_method=self.cluster_method, cluster_dist_cutoff=self.cluster_dist_cutoff)
            chg_ent_keyword_list = sorted(chg_ent_keyword_list)
            # Save calculted data in case job is unexpectedly terminated
            np.savez(npz_data_file,
                    chg_ent_fingerprint_list=chg_ent_fingerprint_list,
                    Q_list=Q_list,
                    chg_ent_keyword_dict=chg_ent_keyword_dict,
                    chg_ent_keyword_list=chg_ent_keyword_list,
                    idx2trajfile=idx2trajfile,
                    idx2frame=idx2frame,
                    ent_cluster_data=ent_cluster_data,
                    ent_cluster_tree=ent_cluster_tree)
            self.logger.info(f'SAVED: {npz_data_file}')
            
        else:
            self.logger.info(f'Reading clustering data from {npz_data_file}...')
            npz_data = np.load(npz_data_file, allow_pickle=True)
            chg_ent_fingerprint_list = npz_data['chg_ent_fingerprint_list'].tolist()
            Q_list = npz_data['Q_list'].tolist()        
            chg_ent_keyword_dict = npz_data['chg_ent_keyword_dict'].item()
            chg_ent_keyword_list = npz_data['chg_ent_keyword_list'].tolist()
            idx2frame = npz_data['idx2frame'].tolist()
            idx2trajfile = npz_data['idx2trajfile'].tolist()
            dtrajs = [[[] for frame in chg_ent_fingerprint.keys()] for chg_ent_fingerprint in chg_ent_fingerprint_list]
            ent_cluster_data = npz_data['ent_cluster_data'].item()
            ent_cluster_tree = npz_data['ent_cluster_tree'].item()


        ent_cluster_idx_map = []
        for ent_keyword in chg_ent_keyword_list:
            for i  in range(len(ent_cluster_data[ent_keyword])):
                ent_cluster_idx_map.append([ent_keyword, i])

        ## Print and save cluster tree
        cluster_headers = ['After clustering on N crossing', 'After clustering on C crossing', 'After clustering on loop']
        cluster_tree_file = f'cluster_tree_{"_".join(self.classify_key)}.dat'
        cluster_tree_file = os.path.join(self.outdir, cluster_tree_file)
        self.logger.info(f'Making {cluster_tree_file}')
        with open(cluster_tree_file, 'w') as f:
            for ent_keyword in chg_ent_keyword_list:
                f.write(ent_keyword+'\n')
                clusters = ent_cluster_tree[ent_keyword]
                for i in range(len(clusters)):
                    f.write(' '*4 + cluster_headers[i] + ':\n')
                    for cluster in clusters[i]:
                        f.write(' '*8 + '[')
                        for ci, c in enumerate(cluster):
                            cluster_id = ent_cluster_idx_map.index([ent_keyword, c])+1
                            if ci == 0:
                                f.write('%d'%cluster_id)
                            else:
                                f.write(', %d'%cluster_id)
                        f.write(']\n')
                f.write('\n')
        self.logger.info(f'SAVED: {cluster_tree_file}')

        # Find representative changes of entanglement in each cluster
        rep_chg_ent_list_file = f'rep_chg_ent_list_{"_".join(self.classify_key)}.pkl'
        rep_chg_ent_list_file = os.path.join(self.outdir, rep_chg_ent_list_file)
        rep_chg_ent_data_file = f'rep_chg_ent_{"_".join(self.classify_key)}.csv'
        rep_chg_ent_data_file = os.path.join(self.outdir, rep_chg_ent_data_file)
        if os.path.exists(rep_chg_ent_list_file) and os.path.exists(rep_chg_ent_data_file):
            self.logger.debug('Reading representative changes of entanglement...')
            with open(rep_chg_ent_list_file, 'rb') as f:
                rep_chg_ent_list = pickle.load(f)
            self.logger.debug(f'Loaded: {rep_chg_ent_data_file} into rep_chg_ent_list')

        else:
            self.logger.debug('Finding representative changes of entanglement...')
            rep_chg_ent_list = self.find_representative_entanglement(ent_cluster_data, ent_cluster_idx_map)
            with open(rep_chg_ent_list_file, 'wb') as f: # save the list as a pickle file
                pickle.dump(rep_chg_ent_list, f)
            self.logger.info(f'SAVED: {rep_chg_ent_list_file}')

            # Create dataframe and save data
            data = []
            column_list = ['Keywords', 'Trajectory', 'Frame', 'Native Contact (Residue Index)',
                        'Ref N-ter Crossing', 'Ref C-ter Crossing', 'N-ter Crossing', 'C-ter Crossing',
                        'Ref N-ter GLN', 'Ref C-ter GLN', 'N-ter GLN', 'C-ter GLN',
                        'Ref N-ter Linking Number', 'Ref C-ter Linking Number', 'N-ter Linking Number', 'C-ter Linking Number']
            index_list = []
            for state_id, rep_chg_ent in enumerate(rep_chg_ent_list):
                index_list.append(state_id+1)
                [traj_idx, frame_idx] = rep_chg_ent[:2]
                nc = tuple(rep_chg_ent[2:])
                keyword = ent_cluster_idx_map[state_id][0]
                chg_ent_fingerprint = chg_ent_fingerprint_list[traj_idx][frame_idx][nc]
                cross = []
                for i in range(len(chg_ent_fingerprint['crossing_resid'])):
                    cross.append([])
                    for j in range(len(chg_ent_fingerprint['crossing_resid'][i])):
                        cross[-1].append(chg_ent_fingerprint['crossing_pattern'][i][j]+str(chg_ent_fingerprint['crossing_resid'][i][j]))
                ref_cross = []
                for i in range(len(chg_ent_fingerprint['ref_crossing_resid'])):
                    ref_cross.append([])
                    for j in range(len(chg_ent_fingerprint['ref_crossing_resid'][i])):
                        ref_cross[-1].append(chg_ent_fingerprint['ref_crossing_pattern'][i][j]+str(chg_ent_fingerprint['ref_crossing_resid'][i][j]))
                GLN = chg_ent_fingerprint['linking_value']
                ref_GLN = chg_ent_fingerprint['ref_linking_value']
                LN = chg_ent_fingerprint['topoly_linking_number']
                ref_LN = chg_ent_fingerprint['ref_topoly_linking_number']

                data_0 = [keyword, idx2trajfile[traj_idx], frame_idx, nc,
                        ref_cross[0], ref_cross[1], cross[0], cross[1],
                        ref_GLN[0], ref_GLN[1], GLN[0], GLN[1],
                        ref_LN[0], ref_LN[1], LN[0], LN[1]]
                data.append(data_0)

            df = pd.DataFrame(data, columns=column_list, index=index_list)
            rep_chg_ent_data_file = f'rep_chg_ent_{"_".join(self.classify_key)}.csv'
            rep_chg_ent_data_file = os.path.join(self.outdir, rep_chg_ent_data_file)
            df.to_csv(rep_chg_ent_data_file, index_label='State ID')
            self.logger.info(f'SAVED: {rep_chg_ent_data_file}')

        # plot entanglement distribution
        n_cluster = len(ent_cluster_idx_map)
        fig = plt.figure(figsize=(np.max([6, 0.3*n_cluster]),5))
        ax = fig.add_subplot(1,1,1)
        window_width = 0.8
        for state_id, [key, cluster_idx] in enumerate(ent_cluster_idx_map):
            cluster = ent_cluster_data[key][cluster_idx]
            nc_list = [c[2:4] for c in cluster]
            sort_index = [i for i, x in sorted(enumerate(nc_list), key=lambda x: (x[1][1]-x[1][0], x[1][0]))]
            sort_index = np.array(sort_index)
            if len(sort_index) <= self.max_plot_samples:
                plot_idx = np.arange(0, len(sort_index), 1, dtype=int)
            else:
                plot_idx = np.linspace(0, len(sort_index)-1, self.max_plot_samples, dtype=int)
            for idx, ci in enumerate(sort_index[plot_idx]):
                c = cluster[ci]
                nc = c[2:4]
                traj_idx = c[0]
                frame_idx = c[1]
                fingerprint = chg_ent_fingerprint_list[traj_idx][frame_idx][tuple(nc)]
                crossings = fingerprint['crossing_resid']
                ref_crossings = fingerprint['ref_crossing_resid']
                # plot loop
                x = state_id+1-window_width/2 + (idx+1)*window_width/(len(plot_idx)+1)
                ax.plot([x,x], nc, '-', color='tomato', linewidth=0.5, alpha=0.4)
                # plot crossings
                x = state_id+1-window_width/2 + (idx+1)*window_width/(len(plot_idx)+1)
                for ccr in ref_crossings:
                    for cc in ccr:
                        ax.plot([x, x], [cc-0.5, cc+0.5], '-', color='green', linewidth=0.5, alpha=0.4)
                for ccr in crossings:
                    for cc in ccr:
                        ax.plot([x, x], [cc-0.5, cc+0.5], '-', color='blue', linewidth=0.5, alpha=0.4)
        ax.set_xticks(np.arange(1,n_cluster+1,1), np.arange(1,n_cluster+1,1))
        ax.set_xlim([0, n_cluster+1])
        ax.set_xlabel('Cluster')
        ax.set_ylabel('Residue index')
        chg_dist_data_file = f'chg_ent_{"_".join(self.classify_key)}_distribution.pdf'
        chg_dist_data_file = os.path.join(self.outdir, chg_dist_data_file)
        fig.savefig(chg_dist_data_file, bbox_inches='tight')
        self.logger.info(f'SAVED: {chg_dist_data_file}')
        del fig

        self.logger.debug('Clustering structures with unique combinations of changes of entanglements...')
        # Assign entanglement clusters (list of ent_cluster_idx) in discrete trajectories
        for key, ent_clusters in ent_cluster_data.items():
            for i, ent_cluster in enumerate(ent_clusters):
                cluster_id = ent_cluster_idx_map.index([key, i])
                for chg_ent_keyword in ent_cluster:
                    traj_idx, frame = chg_ent_keyword[:2]
                    dtrajs[traj_idx][idx2frame[traj_idx].index(frame)].append(cluster_id)
        self.logger.info(f'Assigned entanglement clusters in discrete trajectories')

        # Strip same cluster ids in each frame
        chg_ent_structure_keyword_list = []
        for dtraj in dtrajs:
            for i, cluster_id_list in enumerate(dtraj):
                dtraj[i] = sorted(list(set(cluster_id_list)))
                if str(dtraj[i]) not in chg_ent_structure_keyword_list:
                    chg_ent_structure_keyword_list.append(str(dtraj[i]))
        chg_ent_structure_keyword_list = sorted(chg_ent_structure_keyword_list)
        self.logger.info(f'Stripped same cluster ids in each frame')

        # cluster trajectory frames with different combinations of changes in entanglement
        chg_ent_structure_cluster_data = {chg_ent_structure_keyword: [] for chg_ent_structure_keyword in chg_ent_structure_keyword_list}
        for traj_idx, dtraj in enumerate(dtrajs):
            for frame_idx, cluster_id_list in enumerate(dtraj):
                chg_ent_structure_cluster_data[str(cluster_id_list)].append([traj_idx, idx2frame[traj_idx][frame_idx]])
        Num_struct_list = [len(chg_ent_structure_cluster_data[keyword]) for keyword in chg_ent_structure_keyword_list]
        sort_idx = np.argsort(-np.array(Num_struct_list, dtype=int))
        sorted_chg_ent_structure_keyword_list = [chg_ent_structure_keyword_list[idx] for idx in sort_idx]
        sorted_Num_struct_list = [Num_struct_list[idx] for idx in sort_idx]
        self.logger.info(f'Cluster trajectory frames with different combinations of changes in entanglement')

        self.logger.debug('Find representative combinations of changes of entanglements in structures...')
        # Find representative changes of entanglement (minimal loop) in each frame
        rep_chg_ent_dtrajs = []
        for traj_idx, dtraj in enumerate(dtrajs):
            #print(f'\nTRAJ IDX: {traj_idx} with {len(dtraj)} frames')
            rep_chg_ent_dtrajs.append([])
            for frame_idx, cluster_id_list in enumerate(dtraj):
                #print(f'FRAME IDX: {frame_idx} {cluster_id_list}')
                frame_idx_0 = idx2frame[traj_idx][frame_idx]
                rep_chg_ent_dtrajs[-1].append({})
                for cluster_id in cluster_id_list:
                    [ent_keyword, idx] = ent_cluster_idx_map[cluster_id]
                    nc_list = []
                    for element in ent_cluster_data[ent_keyword][idx]:
                        if traj_idx == element[0] and frame_idx_0 == element[1]:
                            nc_list.append(tuple(element[2:]))
                    rep_nc = nc_list[0]
                    for nc in nc_list:
                        if nc[1]-nc[0] < rep_nc[1]-rep_nc[0]:
                            rep_nc = nc
                    rep_chg_ent_dtrajs[-1][-1][cluster_id] = chg_ent_fingerprint_list[traj_idx][frame_idx_0][rep_nc]
            
        # Find representative structures (max Q) for each combination
        self.logger.info(f'Find representative structures (max Q) for each combination')
        rep_struct_data = {}
        for keyword in sorted_chg_ent_structure_keyword_list:
            Q_0 = 0
            rep_struct_data[keyword] = chg_ent_structure_cluster_data[keyword][0]
            for [traj_idx, frame_idx] in chg_ent_structure_cluster_data[keyword]:
                Q = Q_list[traj_idx][frame_idx]
                if Q > Q_0:
                    Q_0 = Q
                    rep_struct_data[keyword] = [traj_idx, frame_idx]
        self.logger.debug('Found representative structures (max Q) for each combination')

        # Create dataframe and save data
        #chg_ent_data_file = f'chg_ent_struct_{"_".join(self.classify_key)}.csv'
        #chg_ent_data_file = os.path.join(self.outdir, chg_ent_data_file)
        #if os.path.exitst(chg_ent_data_file):
        #    print(f'Reading {chg_ent_data_file}')
        #    df = pd.read_csv(chg_ent_data_file, index_col='State ID')
        data = []
        column_list = ['Rep_chg_ents', 'Num of structures', 'Probability',
                    'Rep trajectory', 'Rep frame', 'Max Q', 'Median Q']
        index_list = []
        tot_num_frames = 0
        for traj_idx, dtraj in enumerate(dtrajs):
            tot_num_frames += len(dtraj)
        for state_id, keyword in enumerate(sorted_chg_ent_structure_keyword_list):
            index_list.append(state_id+1)
            Rep_chg_ents = str(list(np.array(eval(keyword))+1))
            Num = len(chg_ent_structure_cluster_data[keyword])
            Prob = Num / tot_num_frames
            Q_0_list = [Q_list[cd[0]][cd[1]] for cd in chg_ent_structure_cluster_data[keyword]]
            max_Q = np.max(Q_0_list)
            median_Q = np.median(Q_0_list)
            data_0 = [Rep_chg_ents, Num, Prob, idx2trajfile[rep_struct_data[keyword][0]], rep_struct_data[keyword][1], max_Q, median_Q]
            data.append(data_0)
        df = pd.DataFrame(data, columns=column_list, index=index_list)
        chg_ent_data_file = f'chg_ent_struct_{"_".join(self.classify_key)}.csv'
        chg_ent_data_file = os.path.join(self.outdir, chg_ent_data_file)
        df.to_csv(chg_ent_data_file, index_label='State ID')
        self.logger.info(f'SAVED: {chg_ent_data_file}')

        ## determine if there is any issue with the item shapes before saving
        save_items = {
        "chg_ent_fingerprint_list": chg_ent_fingerprint_list,
        "Q_list": Q_list,
        "chg_ent_keyword_dict": chg_ent_keyword_dict,
        "chg_ent_keyword_list": chg_ent_keyword_list,
        "idx2trajfile": idx2trajfile,
        "idx2frame": idx2frame,
        # "RMSD_array": RMSD_array,
        "ent_cluster_data": ent_cluster_data,
        "ent_cluster_tree": ent_cluster_tree,
        "rep_chg_ent_list": rep_chg_ent_list,
        "dtrajs": dtrajs,
        "rep_chg_ent_dtrajs": rep_chg_ent_dtrajs,
        "sorted_chg_ent_structure_keyword_list": sorted_chg_ent_structure_keyword_list,
        "chg_ent_structure_cluster_data": chg_ent_structure_cluster_data,
        "rep_struct_data": rep_struct_data}

        for key, value in save_items.items():
            try:
                shape_info = f", shape={np.shape(value)}"
            except Exception as e:
                shape_info = f", shape=unavailable ({e})"
            self.logger.info(f"{key}: type={type(value)}{shape_info}")


        # Save data
        np.savez(npz_data_file,
                chg_ent_fingerprint_list=chg_ent_fingerprint_list,
                Q_list=Q_list,
                chg_ent_keyword_dict=chg_ent_keyword_dict,
                chg_ent_keyword_list=chg_ent_keyword_list,
                idx2trajfile=idx2trajfile,
                idx2frame=idx2frame,
                # RMSD_array=RMSD_array,
                ent_cluster_data=ent_cluster_data,
                ent_cluster_tree=ent_cluster_tree,
                rep_chg_ent_list=rep_chg_ent_list,
                dtrajs=np.array(dtrajs, dtype=object),
                rep_chg_ent_dtrajs=np.array(rep_chg_ent_dtrajs, dtype=object),
                sorted_chg_ent_structure_keyword_list=sorted_chg_ent_structure_keyword_list,
                chg_ent_structure_cluster_data=chg_ent_structure_cluster_data,
                rep_struct_data=rep_struct_data)
        self.logger.info(f'SAVED: {npz_data_file}')
        self.logger.info(f'Clustering Complete')
    ##########################################################################################################################################################
    
    ##########################################################################################################################################################
    def viz_rep_changes(self, ):
        psf_file = None
        if_viz = 1
        if_backmap = 0
        pulchra_only = False
        native_AA_pdb = None
        top_struct = 0.01

        if if_viz:
            ## Generate visualiztion for representative changes of entanglement in each cluster
            self.logger.debug('Generate visualization for representative changes of entanglement...')
            os.system('mkdir viz_rep_chg_ent_%s'%('_'.join(self.classify_key)))
            for state_id, rep_chg_ent in enumerate(rep_chg_ent_list):
                state_cor = mdt.load(idx2trajfile[rep_chg_ent[0]], top=psf_file)[rep_chg_ent[1]].center_coordinates().xyz*10
                nc = (rep_chg_ent[2], rep_chg_ent[3])
                chg_ent_fingerprint = chg_ent_fingerprint_list[rep_chg_ent[0]][rep_chg_ent[1]][nc]
                rep_ent_dict = {tuple(chg_ent_fingerprint['code']): [chg_ent_fingerprint]}
                os.chdir('viz_rep_chg_ent_%s'%('_'.join(self.classify_key)))
                gen_state_visualizion(state_id+1, psf_file, state_cor, native_AA_pdb, rep_ent_dict, if_backmap=if_backmap, pulchra_only=pulchra_only)
                os.chdir('../')
            
            # Generate visualiztion for representative changes of entanglements in each structural cluster
            self.logger.debug('Generate visualization for unique entangled structures...')
            if top_struct >= 1:
                viz_dir = 'viz_chg_ent_struct_%s_%d'%('_'.join(self.classify_key), top_struct)
            else:
                viz_dir = 'viz_chg_ent_struct_%s_%.4f'%('_'.join(self.classify_key), top_struct)
            os.system('mkdir %s'%viz_dir)
            for state_id, keyword in enumerate(sorted_chg_ent_structure_keyword_list):
                if top_struct >= 1 and state_id >= top_struct:
                    break
                elif top_struct < 1 and sorted_Num_struct_list[state_id]/tot_num_frames < top_struct:
                    break
                [traj_idx, frame_idx] = rep_struct_data[keyword]
                state_cor = mdt.load(idx2trajfile[traj_idx], top=psf_file)[frame_idx].center_coordinates().xyz*10
                frame_idx_0 = idx2frame[traj_idx].index(frame_idx)
                rep_chg_ent_dict = rep_chg_ent_dtrajs[traj_idx][frame_idx_0]
                rep_ent_dict = {tuple(v['code']): [] for k, v in rep_chg_ent_dict.items()}
                for k, v in rep_chg_ent_dict.items():
                    rep_ent_dict[tuple(v['code'])].append(v)
                os.chdir(viz_dir)
                gen_state_visualizion(state_id+1, psf_file, state_cor, native_AA_pdb, rep_ent_dict, if_backmap=if_backmap, pulchra_only=pulchra_only)
                os.chdir('../')
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def gen_state_visualizion(self, state_id, psf, state_cor, native_AA_pdb, rep_ent_dict, if_backmap=True, pulchra_only=False):
        def idx2sel(idx_list):
            if len(idx_list) == 0:
                return ''
            else:
                sel = 'index'
                idx_0 = idx_list[0]
                idx_1 = idx_list[0]
                sel_0 = ' %d'%idx_0
                for i in range(1, len(idx_list)):
                    if idx_list[i] == idx_list[i-1] + 1:
                        idx_1 = idx_list[i]
                    else:
                        if idx_1 > idx_0:
                            sel_0 += ' to %d'%idx_1
                        sel += sel_0
                        idx_0 = idx_list[i]
                        idx_1 = idx_list[i]
                        sel_0 = ' %d'%idx_0
                if idx_1 > idx_0:
                    sel_0 += ' to %d'%idx_1
                sel += sel_0
                return sel

        AA_name_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
                        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
                        'HIE', 'HID', 'HIP']
        protein_colorid_list = [6, 6]
        loop_colorid_list = [1, 1]
        thread_colorid_list = [0, 0]
        nc_colorid_list = [3, 3]
        crossing_colorid_list = [8, 8]
        thread_cutoff=3
        terminal_cutoff=3
        
        self.logger.info('Generate visualization of state %d'%(state_id))

        struct = pmd.load_file(psf)
        struct.coordinates = state_cor

        # backmap
        if if_backmap:
            if pulchra_only:
                pulchra_only = '1'
            else:
                pulchra_only = '0'
            struct.save('tmp.pdb', overwrite=True)
            os.system('backmap.py -i '+native_AA_pdb+' -c tmp.pdb -p '+pulchra_only)
            os.system('mv tmp_rebuilt.pdb state_%d.pdb'%state_id)
            os.system('rm -f tmp.pdb')
            os.system('rm -rf ./rebuild_tmp/')
        else:
            struct.save('state_%d.pdb'%state_id, overwrite=True)
        
        ref_struct = pmd.load_file(native_AA_pdb)
        current_struct = pmd.load_file('state_%d.pdb'%state_id)

        if len(list(rep_ent_dict.keys())) == 0:
            # no change of entaglement
            f = open('vmd_s%d_none.tcl'%(state_id), 'w')
            f.write('# Entanglement type: no change\n')
            f.write('''display rendermode GLSL
                        axes location off

                        color Display {Background} white

                        mol new ./'''+('state_%d.pdb'%state_id)+''' type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
                        mol delrep 0 top
                        mol representation NewCartoon 0.300000 10.000000 4.100000 0
                        mol color ColorID '''+str(protein_colorid_list[1])+'''
                        mol selection {all}
                        mol material AOChalky
                        mol addrep top
                        ''')
            f.close()

        # Create vmd script for each type of change
        for ent_code, rep_ent_list in rep_ent_dict.items():
            pmd_struct_list = [ref_struct, current_struct]
            struct_dir_list = [native_AA_pdb, './state_%d.pdb'%state_id]
            key_prefix_list = ['ref_', '']
            repres_list = ['', '']
            align_sel_list = ['', '']
            
            vmd_script = '''# Entanglement type: '''+str(rep_ent_list[0]['type'])+'''
                            package require topotools
                            display rendermode GLSL
                            axes location off

                            color Display {Background} white

                            '''
            for struct_idx, pmd_struct in enumerate(pmd_struct_list):
                struct_dir = struct_dir_list[struct_idx]
                protein_colorid = protein_colorid_list[struct_idx]
                loop_colorid = loop_colorid_list[struct_idx]
                thread_colorid = thread_colorid_list[struct_idx]
                nc_colorid = nc_colorid_list[struct_idx]
                crossing_colorid = crossing_colorid_list[struct_idx]
                key_prefix = key_prefix_list[struct_idx]

                # Clean ligands
                clean_sel_idx = np.zeros(len(pmd_struct.atoms))
                for res in pmd_struct.residues:
                    if res.name in AA_name_list:
                        for atm in res.atoms:
                            clean_sel_idx[atm.idx] = 1
                pmd_clean_struct = pmd_struct[clean_sel_idx]
                clean_idx_to_idx = np.where(clean_sel_idx == 1)[0]

                # vmd selection string for protein
                idx_list = []
                for res in pmd_struct.residues:
                    if res.name in AA_name_list:
                        idx_list += [atm.idx for atm in res.atoms]
                vmd_sel = idx2sel(idx_list)

                repres = '''mol new '''+struct_dir+''' type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
                            mol delrep 0 top
                            mol representation NewCartoon 0.300000 10.000000 4.100000 0
                            mol color ColorID '''+str(protein_colorid)+'''
                            mol selection {'''+vmd_sel+'''}
                            mol material AOChalky
                            mol addrep top
                            '''
                align_sel = vmd_sel
                for chg_ent_fingerprint in rep_ent_list:
                    nc = chg_ent_fingerprint[key_prefix+'native_contact']

                    idx_list = []
                    for res in pmd_clean_struct.residues:
                        if res.idx in nc:
                            idx_list += [atm.idx for atm in res.atoms if atm.name == 'CA']
                    nc_sel = idx2sel(clean_idx_to_idx[idx_list])

                    idx_list = []
                    for res in pmd_clean_struct.residues:
                        if res.idx >= nc[0] and res.idx <= nc[1]:
                            idx_list += [atm.idx for atm in res.atoms]
                    loop_sel = idx2sel(clean_idx_to_idx[idx_list])

                    align_sel += ' and not (%s)'%loop_sel
                    ref_coss_resid = chg_ent_fingerprint['ref_crossing_resid']
                    cross_resid = chg_ent_fingerprint['crossing_resid']
                    thread = []
                    thread_sel_list = []
                    for ter_idx in range(len(ref_coss_resid)):
                        thread_0 = []
                        resid_list = ref_coss_resid[ter_idx] + cross_resid[ter_idx]
                        if len(resid_list) > 0:
                            thread_0 = [np.min(resid_list)-5, np.max(resid_list)+5]
                            if ter_idx == 0:
                                thread_0[0] = np.max([thread_0[0], terminal_cutoff])
                                thread_0[1] = np.min([thread_0[1], nc[0]-thread_cutoff])
                            else:
                                thread_0[0] = np.max([thread_0[0], nc[1]+thread_cutoff])
                                thread_0[1] = np.min([thread_0[1], len(struct.atoms)-1-terminal_cutoff])
                            idx_list = []
                            for res in pmd_clean_struct.residues:
                                if res.idx >= thread_0[0] and res.idx <= thread_0[1]:
                                    idx_list += [atm.idx for atm in res.atoms]
                            thread_0_sel = idx2sel(clean_idx_to_idx[idx_list])
                            thread_sel_list.append(thread_0_sel)
                            align_sel += ' and not (%s)'%thread_0_sel
                        else:
                            thread_sel_list.append('')
                        thread.append(thread_0)

                    ln = chg_ent_fingerprint[key_prefix+'topoly_linking_number']
                    cross = []
                    for i in range(len(chg_ent_fingerprint[key_prefix+'crossing_resid'])):
                        cross.append([])
                        for j in range(len(chg_ent_fingerprint[key_prefix+'crossing_resid'][i])):
                            cross[-1].append(chg_ent_fingerprint[key_prefix+'crossing_pattern'][i][j]+str(chg_ent_fingerprint[key_prefix+'crossing_resid'][i][j]))
                    repres += '# idx: native contact %s, linking number %s, crossings %s.\n'%(str(nc), str(ln), str(cross))
                    repres +='''mol representation NewCartoon 0.350000 10.000000 4.100000 0
                                mol color ColorID '''+str(loop_colorid)+'''
                                mol selection {'''+loop_sel+'''}
                                mol material Opaque
                                mol addrep top
                                mol representation VDW 1.000000 12.000000
                                mol color ColorID '''+str(nc_colorid)+'''
                                mol selection {'''+nc_sel+'''}
                                mol material Opaque
                                mol addrep top
                                set sel [atomselect top "'''+nc_sel+'''"]
                                set idx [$sel get index]
                                topo addbond [lindex $idx 0] [lindex $idx 1]
                                mol representation Bonds 0.300000 12.000000
                                mol color ColorID '''+str(nc_colorid)+'''
                                mol selection {'''+nc_sel+'''}
                                mol material Opaque
                                mol addrep top
                                '''
                    for ter_idx, thread_resid in enumerate(thread):
                        if len(thread_resid) == 0:
                            continue
                        repres += '''mol representation NewCartoon 0.350000 10.000000 4.100000 0
                                    mol color ColorID '''+str(thread_colorid)+'''
                                    mol selection {'''+thread_sel_list[ter_idx]+'''}
                                    mol material Opaque
                                    mol addrep top
                                    '''
                        if len(chg_ent_fingerprint[key_prefix+'crossing_resid'][ter_idx]) > 0:
                            idx_list = []
                            for res in pmd_clean_struct.residues:
                                if res.idx in chg_ent_fingerprint[key_prefix+'crossing_resid'][ter_idx]:
                                    idx_list += [atm.idx for atm in res.atoms if atm.name == 'CA']
                            crossing_sel = idx2sel(clean_idx_to_idx[idx_list])
                            repres += '''mol representation VDW 1.000000 12.000000
                                        mol color ColorID '''+str(crossing_colorid)+'''
                                        mol selection {'''+crossing_sel+'''}
                                        mol material Opaque
                                        mol addrep top
                                        '''
                    
                if struct_idx == 0:
                    repres += '''mol representation VDW 1.000000 12.000000
                                mol color Name
                                mol selection {not ('''+vmd_sel+''') and not water}
                                mol material Opaque
                                mol addrep top
                                '''
                repres_list[struct_idx] = repres
                align_sel_list[struct_idx] = align_sel

            vmd_script += '\n'.join(repres_list)
            vmd_script += '''
                            set sel1 [atomselect 0 "'''+align_sel_list[0]+''' and name CA"]
                            set sel2 [atomselect 1 "'''+align_sel_list[1]+''' and name CA"]
                            set trans_mat [measure fit $sel1 $sel2]
                            set move_sel [atomselect 0 "all"]
                            $move_sel move $trans_mat
                            '''
            f = open('vmd_s%d_n%s_c%s.tcl'%(state_id, ent_code[0], ent_code[1]), 'w')
            f.write(vmd_script)
            f.close()
    ##########################################################################################################################################################

##########################################################################################################################################################
##########################################################################################################################################################

##########################################################################################################################################################
class MSMNonNativeEntanglementClustering:
    """
    Build a markov state model across an ensemble of protein structures with non-native entanglements to identify metastable states
    """

    #######################################################################################
    def __init__(self, OPpath:str = './', outdir:str = './', ID:str = '', 
        start:int= 0 , end:int = 99999999999, stride:int = 1, ITS:str = 'False', lagtime:int = 1,
        n_cluster:int = 400, kmean_stride:int = 2, n_small_states:int = 1, n_large_states:int = 10, dt:float = 0.015/1000, rm_traj_list:list = [], log_level:int = logging.INFO, logdir:str = None):
        
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        ("--outdir", type=str, required=True, help="Path to output directory")
        ("--OPpath", type=str, required=True, help="Path to directory containing G and Q directories created by GQ.py")
        ("--ID", type=str, required=True, help="base name for output files")
        ("--start", type=int, required=False, help="First frame to analyze 0 indexed", default=0)
        ("--end", type=int, required=False, help="Last frame to analyze 0 indexed", default=-1)
        ("--stride", type=int, required=False, help="Frame stride", default=1)
        ("--ITS", type=str, required=False, help="Find optimal lag time with ITS", default='False')
        ("--lagtime", type=int, required=False, help="lagtime to build the model", default=1)
        ("--n_cluster", type=int, required=False, help="Number of k-means clusters to group. Default is 400.")
        ("--kmean_stride", type=int, required=False, help="Stride of reading trajectory frame when clustring by k-means. Default is 2.")
        ("--n_small_states", type=int, required=False, help="Number of clusters for the inactive microstates after MSM modeling to be clustered into", default=1)
        ("--n_large_states", type=int, required=False, help="Number of clusters for the active microstates after MSM modeling to be clustered into", default=10)
        ("--dt", type=float, required=False, help="timestep used in MD simulations in ns", default=0.015/1000)
        ("--rm_traj_list", type=str, nargs='+', required=False, help="List of trajectory numbers to remove from analysis", default=[])
        """

        # parse the parameters 
        self.outdir = outdir
        self.ID = ID
        self.logger = setup_logger('MSMNonNativeEntanglementClustering', outdir=logdir if logdir is not None else outdir, ID=ID, log_level=log_level)

        self.OPpath = OPpath
        self.logger.debug(f'OPpath: {self.OPpath}')
        self.logger.debug(f'outdir: {self.outdir}')
        self.logger.debug(f'ID: {self.ID}')

        self.ITS = ITS
        self.logger.debug(f'ITS: {ITS}')

        self.lagtime = lagtime
        self.logger.debug(f'lagtime: {lagtime}')

        #self.dcds = args.dcds 
        #print(f'dcds: {self.dcds}')

        self.start = start
        self.end = end
        self.stride = stride
        self.logger.debug(f'START: {self.start} | END: {self.end} | STRIDE: {self.stride}')

        self.n_cluster = n_cluster # Number of k-means clusters to group. Default is 400.
        self.kmean_stride = kmean_stride # Stride of reading trajectory frame when clustring by k-means. 
        self.n_small_states = n_small_states # Number of clusters for the inactive microstates after MSM modeling to be clustered into
        self.n_large_states = n_large_states  # Adjust based on your system
        self.dt = dt # timestep used in MD simulations in ns

        self.rm_traj_list = rm_traj_list
        self.logger.info(f'Trajectories to ignore: {self.rm_traj_list}')

    #######################################################################################

    #######################################################################################
    def load_OP(self,):
        """
        Loads the GQ values of each trajectory into a 2D array and then appends it to a list
        The list should have Nt = number of trajectories and each array should be n x 2 where n is the number of frames
        """
        self.logger.info(f'Loading G and Q order parameters...')
        cor_list = []
        cor_list_idx_2_traj = {}
        Qfiles = glob.glob(os.path.join(self.OPpath, 'Q/*.Q'))
        QTrajs = [int(pathlib.Path(Qf).stem.split('Traj')[-1]) for Qf in Qfiles]

        Gfiles = glob.glob(os.path.join(self.OPpath, 'G/*.G'))
        GTrajs = [int(pathlib.Path(Gf).stem.split('Traj')[-1]) for Gf in Gfiles]

        shared_Trajs = set(QTrajs).intersection(GTrajs)
        #print(f'Shared Traj between Q and G: {shared_Trajs} {len(shared_Trajs)}')
        self.logger.info(f'Number of Q files found: {len(Qfiles)} | Number of G files found: {len(Gfiles)}')
        self.logger.info(f'Number of shared Traj between Q and G: {len(shared_Trajs)}')


        ## remove trajectories that are in the rm_traj_list
        if len(self.rm_traj_list) > 0:
            self.logger.info(f'Removing trajectories: {self.rm_traj_list}')
            shared_Trajs = [traj for traj in shared_Trajs if traj not in self.rm_traj_list]
            self.logger.info(f'Number of shared Traj after removing: {len(shared_Trajs)}')


        # loop through the Qfiles and find matching Gfile
        # then load the Q and G time series into a 2D array
        idx = 0
        QFrames = {}
        GFrames = {}
        for traj in shared_Trajs:
            #print(f'Traj: {traj}')

            # get the cooresponding G and Q file
            Qf = [f for f in Qfiles if f.endswith(f'Traj{traj}.Q')]
            Gf = [f for f in Gfiles if f.endswith(f'Traj{traj}.G')]
            self.logger.debug(f'Qf: {Qf}')
            self.logger.debug(f'Gf: {Gf}')

            ## Quality check to assert that only a single G and Q file were found
            assert len(Qf) == 1, f"the number of Q files {len(Qf)} should equal 1 for Traj {traj}"
            assert len(Gf) == 1, f"the number of G files {len(Gf)} should equal 1 for Traj {traj}"

            # load the G Q data and extract only the time series column
            Qdata = pd.read_csv(Qf[0], sep=',')
            if self.start < 0:
                Qdata = Qdata.iloc[self.start:self.end + 1]
            else:
                Qdata = Qdata[(Qdata['Frame'] >= self.start) & (Qdata['Frame'] <= self.end)]
            #print(Qdata)
            QFrames[traj] = Qdata['Frame'].values
            Qdata = Qdata['total'].values.astype(float)
            


            Gdata = pd.read_csv(Gf[0])
            if self.start < 0:
                Gdata = Gdata.iloc[self.start:self.end + 1]
            else:
                Gdata = Gdata[(Gdata['Frame'] >= self.start) & (Gdata['Frame'] <= self.end)]
            #print(Gdata)
            GFrames[traj] = Gdata['Frame'].values
            Gdata = Gdata['G'].values.astype(float)
            #print(f'Shape of OP: Q {Qdata.shape} G {Gdata.shape}')

            ## Quality check that QFrames == GFrames.  
            if set(QFrames[traj]) != set(GFrames[traj]):
                raise ValueError(f'The frames in Q {QFrames[traj]} do not match the frames in G {GFrames[traj]} for Traj {traj}. Please check your data files.')

            ## Quality check that the G and Q data has the same number of frames
            if Qdata.shape != Gdata.shape:
                self.logger.warning(f"WARNING: The number of frames in Q {Qdata.shape} should equal the number of frames in G {Gdata.shape} in Traj {traj}")
                continue
            
            ## Check and ensure that Qdata or Gdata has no nan values
            if np.isnan(Qdata).any():
                raise ValueError(f'There is a NaN value in this Qdata')
             
            if np.isnan(Gdata).any():
                raise ValueError(f'There is a NaN value in this Gdata')
         
            data = np.stack((Qdata, Gdata)).T
            data = data.astype(float)
            cor_list.append(data)

            cor_list_idx_2_traj[idx] = int(traj)
            idx += 1

        self.logger.info(f'Number of trajecotry OP coordinate loaded: {len(cor_list)}')
        self.cor_list = cor_list
        self.QFrames = QFrames
        self.GFrames = GFrames

        ## Quality check that the number of trajectories loaded is equal to the number of Q and G files
        assert len(cor_list) == len(shared_Trajs), f"The # of coordinates loaded {len(cor_list)} does not equal the number of Q and G files with shared traj after removal of mirror images {len(shared_Trajs)}"

        self.logger.info(f'Mapping of cor_list index to trajID in file names: {cor_list_idx_2_traj}')
        self.cor_list_idx_2_traj = cor_list_idx_2_traj
    #######################################################################################  

    #######################################################################################
    def standardize(self,):
        """
        Standardizes your OP by taking the mean and std Q and G across all traj data and rescaling each trajectorys data by
        Z = (d - mean)/std 
        """
        data_con = self.cor_list[0]
        for i in range(1, len(self.cor_list)):
            data_con = np.vstack((data_con, self.cor_list[i]))
        self.data_mean = np.mean(data_con, axis=0)
        self.data_std = np.std(data_con, axis=0)
        self.standard_cor_list = [(d - self.data_mean) / self.data_std for d in self.cor_list]
    #######################################################################################

    #######################################################################################
    def unstandardize(self, data):
        """
        Unstandardizes your OP by taking the mean and std Q and G across all traj data and rescaling each trajectorys data by
        Z*std + mean = d
        """
        return data*self.data_std + self.data_mean
    #######################################################################################

    #######################################################################################
    def cluster(self,):
        """
        Cluster the GQ data across all trajectories using kmeans. 
        dtrajs contains the resulting kmeans cluster labels for each trajectory time series 
        centers contains the standardized GQ coordinates of the cluster centers

        if the number of unique centers found is not equal to self.n_cluster then adjust it to reflect the number found. This can happen if you have data that has a narrow distribution. 
        """
        self.clusters = pem.coordinates.cluster_kmeans(self.standard_cor_list, k=self.n_cluster, max_iter=5000, stride=self.kmean_stride)
        
        # Get the microstate tagged trajectories and their state counts
        self.dtrajs = self.clusters.dtrajs
        self.logger.debug(f'dtrajs: {len(self.dtrajs)} {self.dtrajs[0].shape}\n{self.dtrajs[0][:10]}')
        clusterIDs, counts = np.unique(self.dtrajs, return_counts=True)
        self.logger.info(f'Number of unique microstate IDs: {len(clusterIDs)} {clusterIDs}')
        
        state_counts = {}
        for i,c in zip(clusterIDs, counts):
            state_counts[i] = c
        self.logger.debug(f'state_counts: {state_counts}')
        
        # Quality check that all microstate ids are assigned
        # If not renumber from 0
        if len(clusterIDs) != self.n_cluster:
            self.logger.info(f'The number of microstate IDs assigned does not match the number specified: {len(clusterIDs)} != {self.n_cluster}')

            mapping_dict = {}
            for new,old in enumerate(clusterIDs):
                mapping_dict[old] = new
            self.logger.debug(f'mapping_dict: {mapping_dict}')

            # Convert the dictionary to a numpy array for efficient mapping
            max_key = max(mapping_dict.keys())
            mapping_array = np.zeros(max_key + 1, dtype=int)
            for key, value in mapping_dict.items():
                mapping_array[key] = value

            # Map the arrays using the mapping array
            self.dtrajs = [mapping_array[arr] for arr in self.dtrajs]
            
            clusterIDs, counts = np.unique(self.dtrajs, return_counts=True)
            self.logger.info(f'Number of unique microstate IDs after mapping: {len(clusterIDs)} {clusterIDs}')
            state_counts = {}
            for i,c in zip(clusterIDs, counts):
                state_counts[i] = c
            self.logger.debug(f'state_counts: {state_counts}')

            self.n_cluster = len(clusterIDs)
        

        standard_centers = self.clusters.clustercenters
        unstandard_centers = self.unstandardize(standard_centers)
        self.logger.info(f'unstandard_centers:\n{unstandard_centers} {unstandard_centers.shape}')
        self.logger.info(f'self.n_cluster: {self.n_cluster}')
        
    #######################################################################################

    #######################################################################################
    def build_msm(self, lagtime=1):
        self.logger.info(f'Building MSM model with a lag time of {lagtime}')

        # Get count matrix and connective groups of microstates
        c_matrix = deeptime.markov.tools.estimation.count_matrix(self.dtrajs, lagtime).toarray()
        self.logger.info(f'c_matrix:\n{c_matrix} {c_matrix.shape}')
        
        sub_groups = deeptime.markov.tools.estimation.connected_sets(c_matrix)
        self.logger.info(f'Total number of sub_groups: {len(sub_groups)}\n{sub_groups}')
        
        # Build the MSM models for any connected sets that have more than 1 microstate
        msm_list = []        
        for sg in sub_groups:
            cm = deeptime.markov.tools.estimation.largest_connected_submatrix(c_matrix, lcc=sg)
            self.logger.info(f'For sub_group: {sg}')
            if len(cm) == 1:
                msm = None
            else:
                self.logger.info(f'Building Transition matrix and MSM model')
                T = deeptime.markov.tools.estimation.transition_matrix(cm, reversible=True)
                msm = pem.msm.markov_model(T, dt_model=str(self.dt)+' ns')
            msm_list.append(msm)
        self.logger.info(f'Number of models: {len(msm_list)}')

        # Coarse grain out the metastable macrostates in the models
        self.logger.info(f'Coarse grain out the metastable macrostates in the models')
        meta_dist = []
        meta_set = []
        eigenvalues_list = []
        for idx_msm, msm in enumerate(msm_list):

            # the first model should contain the largest connected state so use the largest number of metastable states
            # for every other subgroup use the smallest
            if idx_msm == 0:
                n_states = self.n_large_states
            else:
                n_states = self.n_small_states

            if msm == None:
                eigenvalues_list.append(None)
                dist = np.zeros(self.n_cluster)
                iidx = sub_groups[idx_msm][0]
                dist[iidx] = 1.0
                meta_dist.append(dist)
                meta_set.append(sub_groups[idx_msm])

            else:
                eigenvalues_list.append(msm.eigenvalues())
                # coarse-graining 
                while n_states > 1:
                    tag_empty = False
                    pcca = msm.pcca(n_states)
                    for ms in msm.metastable_sets:
                        if ms.size == 0:
                            tag_empty = True
                            break
                    if not tag_empty:
                        break
                    else:
                        n_states -= 1
                        self.logger.info('Reduced number of states to %d for active group %d'%(n_states, idx_msm+1))
                if n_states == 1:
                    # use observation prob distribution for non-active set
                    dist = np.zeros(self.n_cluster)
                    for nas in sub_groups[idx_msm]:
                        for dtraj in dtrajs:
                            dist[nas] += np.count_nonzero(dtraj == nas)
                    dist /= np.sum(dist)
                    meta_dist.append(dist)
                    meta_set.append(sub_groups[idx_msm])
                else:
                    for i, md in enumerate(msm.metastable_distributions):
                        dist = np.zeros(self.n_cluster)
                        s = np.sum(md[msm.metastable_sets[i]])
                        set_0 = []
                        for idx in msm.metastable_sets[i]:
                            iidx = sub_groups[idx_msm][idx]
                            dist[iidx] = md[idx]
                            set_0.append(iidx)
                        dist = dist / s
                        meta_dist.append(dist)
                        meta_set.append(set_0)
        meta_dist = np.array(meta_dist)
        self.logger.debug(f'meta_dist: {len(meta_dist)} {meta_dist.shape}')
        meta_dist_outfile = os.path.join(self.outdir, f'{self.ID}_meta_dist.npy')
        np.save(meta_dist_outfile, meta_dist, allow_pickle=True)
        self.logger.info(f'SAVED: {meta_dist_outfile}')

        # print(f'meta_set: {meta_set}')
        meta_set_df = {'metastable_state':[], 'microstates':[]}
        for i, ms in enumerate(meta_set):
            # print(f'Metastable state {i}: {ms} with {len(ms)} microstates')
            for m in ms:
                meta_set_df['metastable_state'].append(i)
                meta_set_df['microstates'].append(m)

        meta_set_df = pd.DataFrame(meta_set_df)
        self.logger.info(f'Meta set DataFrame:\n{meta_set_df}')
        meta_set_outfile = os.path.join(self.outdir, f'{self.ID}_meta_set.csv')
        meta_set_df.to_csv(meta_set_outfile, index=False)
        self.logger.info(f'SAVED: {meta_set_outfile}')

   
        ## make microstate to metastable state mapping object
        self.logger.info(f'\nMetastable state assignment')
        meta_mapping = {}
        for metaID, microstates in enumerate(meta_set):
            #print(metaID, microstates)
            for m in microstates:
                if m not in meta_mapping:
                    meta_mapping[m] = metaID
                else:
                    raise ValueError(f'Microstate {m} already in a metastable state!')
        self.logger.debug(f'meta_mapping: {meta_mapping} {len(meta_mapping)}')

        # map those microstate states to the metastable state
        metastable_dtraj = []
        for dtraj_idx, dtraj in enumerate(self.dtrajs):
            mapped_dtraj = []
            for d in dtraj:
                mapped_dtraj.append(meta_mapping[d])

            #rint(mapped_dtraj)
            metastable_dtraj += [np.asarray(mapped_dtraj)]

        self.logger.info(f'Metastable state mapping:')
        for dtraj_idx, dtraj in enumerate(metastable_dtraj):
            self.logger.debug(f'dtraj_idx={dtraj_idx} dtrajs[:10]={self.dtrajs[dtraj_idx][:10]} dtraj[:10]={dtraj[:10]} shape={dtraj.shape}')
      

        ## get samples of metastable states by most populated microstates
        self.logger.debug(f'num_dtrajs={len(self.dtrajs)} dtrajs[0].shape={self.dtrajs[0].shape}')
        cluster_indexes = deeptime.markov.sample.compute_index_states(self.dtrajs)
        self.logger.debug(f'cluster_indexes: {len(cluster_indexes)}')


        samples = deeptime.markov.sample.indices_by_distribution(cluster_indexes, meta_dist, 5)
        self.logger.debug(f'samples: {samples} {len(samples)}')
        
        ## Make the output dataframe that has assignments for each frame of each traj
        df = {'traj':[], 'frame':[], 'microstate':[], 'metastablestate':[], 'Q':[], 'G':[], 'StateSample':[]}
        self.logger.info(f'Active & inactive metastable state mapping')
        for k,v in enumerate(metastable_dtraj):
            traj = self.cor_list_idx_2_traj[k]
            #print(k, traj, v[:10])
            for frame, macrostate in enumerate(v):
                microstate = self.dtrajs[k][frame]
                Q = self.cor_list[k][frame, 0]
                G = self.cor_list[k][frame, 1]

                if [k, frame] in samples[macrostate].tolist():
                    StateSample = True
                else:
                    StateSample = False
                #print(k, frame, microstate, macrostate, StateSample)
                df['traj'] += [traj]
                df['frame'] += [self.QFrames[traj][frame]]
                df['microstate'] += [microstate]
                df['metastablestate'] += [macrostate]
                df['Q'] += [Q]
                df['G'] += [G]
                df['StateSample'] += [StateSample]

        df = pd.DataFrame(df)
        #df['frame'] = self.start # correct the frame index to start from the start specified by the user as this frame index starts from 0
        self.logger.info(f'Final MSM mapping DF:\n{df}')
        df_outfile = os.path.join(self.outdir, f'{self.ID}_MSMmapping.csv')
        df.to_csv(df_outfile, index=False)
        self.logger.info(f'SAVED: {df_outfile}')

        # Plot the metastable state membership and free energy surface
        xall = np.hstack([dtraj[:, 0] for dtraj in self.cor_list])
        yall = np.hstack([dtraj[:, 1] for dtraj in self.cor_list])
        states = np.hstack(metastable_dtraj)
        self.logger.debug(f'xall: {xall} {xall.shape}')
        self.logger.debug(f'yall: {yall} {yall.shape}')
        self.logger.debug(f'states: {states} {states.shape}')

        stateplot_outfile = os.path.join(self.outdir, f'{self.ID}_StateAndFEplot.png')
        self.plot_state_map_and_FE(xall, yall, states, stateplot_outfile)
    #######################################################################################

    #######################################################################################
    def plot_state_map_and_FE(self, x, y, states, outfile, cmap='viridis', point_size=50, alpha=0.85, title='State Map'):
        """
        Plots a state map using x and y values colored by state assignments with labeled colorbar.
        Parameters:
            x (array-like): The x-coordinates of the points.
            y (array-like): The y-coordinates of the points.
            states (array-like): The state assignment for each point.
            cmap (str or Colormap): Colormap for state coloring (default is 'viridis').
            point_size (int): Size of the scatter plot points (default is 50).
            alpha (float): Transparency of the points (default is 0.7).
            title (str): Title of the plot (default is 'State Map with Labels').
        """
        #############################################################################################
        # Create a figure and subplots with 1 row and 2 columns
        fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

        ### plot FE surface on left plot
        # Define the number of bins for the 2D histogram
        num_bins = 20

        # Calculate the 2D histogram
        hist, xedges, yedges = np.histogram2d(x, y, bins=num_bins, density=True)

        # Calculate the probability as the histogram values
        probability = hist / np.sum(hist)
        #print(f'probability: {probability} {np.unique(probability)}')

        # Compute the free energy as -log10(probability)
        with np.errstate(divide='ignore'):  # Ignore divide-by-zero warnings
            free_energy = -np.log10(probability)
            #free_energy[np.isinf(free_energy)] = np.nan  # Set infinities to NaN for better plotting
            free_energy[np.isinf(free_energy) | np.isnan(free_energy)] = np.nanmax(free_energy[np.isfinite(free_energy)]) #+ 1  # Replace NaN/Inf with a large value
        self.logger.debug(f'free_energy: {free_energy} {free_energy.shape} {np.unique(free_energy)}')

        # Create the meshgrid for the contour plot
        X, Y = np.meshgrid(xedges[:-1], yedges[:-1])
        self.logger.debug(f'X: {X.shape}\nY: {Y.shape}')

        # Create a custom colormap
        #cmap = plt.cm.viridis
        #cmap = plt.cm.magma
        cmap = plt.cm.gist_ncar

        # Plotting the contour plot
        contour = axes[0].contourf(X, Y, free_energy.T, levels=100, cmap=cmap)  # Transpose to align axes
        fig.colorbar(contour, ax=axes[0], label='Free Energy (-log10 Probability)')
        axes[0].set_xlabel('Q')
        axes[0].set_ylabel('G')
        axes[0].set_title('2D Free Energy Contour Plot')
        axes[0].set_xlim(0,1)


        #############################################################################################
        ## Plot state map
        #_, axes[1], _ = pem.plots.plot_state_map(x, y, states)

        # Create a 2D histogram to determine the bin index for each (x, y) pair
        #############################################################################################
        # Step 1: Identify unique states
        unique_states = np.unique(states)
        n_states = len(unique_states)
        self.logger.debug(f'unique_states: {unique_states} {n_states}')

        # Step 2: Create a colormap with one color per unique state
        # You can use any colormap, or define specific colors if desired
        colors = plt.cm.get_cmap('tab20', n_states)  # 'tab10' has up to 10 colors; change if needed
        cmap = ListedColormap([colors(i) for i in range(n_states)])

        # Step 3: Map states to color indices
        state_to_index = {state: i for i, state in enumerate(unique_states)}
        color_indices = np.vectorize(state_to_index.get)(states)

        # # Step 4: Create scatter plot
        scatter = axes[1].scatter(x, y, c=color_indices, cmap=cmap, s=50, edgecolor='k')  # Customize marker size, etc.

        # Step 5: Add a colorbar with labels
        cbar = plt.colorbar(scatter, ax=axes[1], ticks=np.linspace(0.5, n_states - 1.5, num=n_states), label=f'Metastable States')
        cbar.ax.set_yticklabels(unique_states)  # Label colorbar with the unique state values
        #############################################################################################

        axes[1].set_xlabel('Q')
        axes[1].set_ylabel('G')
        axes[1].set_title('2D state map')
        axes[1].set_xlim(0,1)

        #plt.tight_layout()
        plt.savefig(outfile)
        self.logger.info(f'SAVED: {outfile}')
        plt.clf()

    #######################################################################################

    #######################################################################################
    def plot_implied_timescales(self,):
        """
        Should be done first before building the model to find a proper lagtime for which the timescales (eignenvalues of the transition matrix) of the model are no longer dependant.
        Look for the point or range where the implied timescales stop changing significantly with increasing lag times. 
        This lag time is generally a good choice for building your MSM, as it suggests the dynamics are being captured without undue dependence on the initial conditions.
        """
        #nits = -1
        lag_times = np.arange(1, 100, 10)  # adjust the range based on your system
        n_states = len(np.unique(self.dtrajs))  # or a predefined number of states
        its = pem.msm.its(self.dtrajs, lags=lag_times, errors='bayes')
        pem.plots.plot_implied_timescales(its)
        ITS_outfile = os.path.join(self.outdir, f'{self.ID}_ITS.png')
        plt.savefig(ITS_outfile)
        self.logger.info(f'SAVED: {ITS_outfile}')
    #######################################################################################

    #######################################################################################
    def run(self, ):

        ## make output folder
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
            self.logger.info(f'Made directory: {self.outdir}')
 
        # load the G and Q data
        self.load_OP()

        # apply the standard scalar transformation to the data 
        self.standardize()

        # cluster the standardized data using kmeans clustering with a stride of 10, change this if necessary
        self.cluster()

        # genereate the implied timescales plot to check for a suitable lag time
        # Should be done first before building the model to choose a suitable lag time
        if self.ITS == 'True':
            anal.plot_implied_timescales()
            self.logger.info(f'Analysis terminated since ITS was selected. Check the figure and choose an approrate lagtime')
            quit()
        
        # Build the MSM model with the choosen lagtime
        self.build_msm(lagtime=self.lagtime)
    #######################################################################################

