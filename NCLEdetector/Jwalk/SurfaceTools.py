# ===============================================================================
#     This file is part of Jwalk (Python 3).
#     
#     Jwalk - A tool to calculate the solvent accessible surface distance (SASD) 
#     between crosslinked residues.
#     
#     Copyright 2016 Josh Bullock and Birkbeck College University of London.
# 
#     Jwalk is available under Public Licence.
#     This software is made available under GPL V3
#
#     Please cite your use of Jwalk in published work:
#     
#     J.Bullock, J. Schwab, K. Thalassinos, M. Topf (2016)
#     The importance of non-accessible crosslinks and solvent accessible surface distance
#     in modelling proteins with restraints from crosslinking mass spectrometry. 
#     Molecular and Cellular Proteomics (15) pp.2491-2500
#
# ===============================================================================

import os
import freesasa

def update_crosslink_pairs(crosslink_pairs, aa1_CA, aa2_CA, remove_aa1, remove_aa2):
    
    '''Removes buried residues from crosslink_pairs'''
    
    buried_residues = []
    index_to_delete = []
    
    for i in range(len(crosslink_pairs)): # for each residue pair, check both are solvent accessible
        
        xl_pair_1, xl_pair_2 = crosslink_pairs[i]
        
        if xl_pair_1 not in aa1_CA:
            index_to_delete.append(i)
            if xl_pair_1 not in buried_residues:
                buried_residues.append(xl_pair_1)
            if xl_pair_2 not in aa2_CA and xl_pair_2 not in buried_residues:
                buried_residues.append(xl_pair_2)
        elif xl_pair_2 not in aa2_CA:
            index_to_delete.append(i)
            if xl_pair_2 not in buried_residues:  
                buried_residues.append(xl_pair_2)
                
        if [xl_pair_1[0],xl_pair_1[1]] in remove_aa1:
            index_to_delete.append(i)
            if xl_pair_1 not in buried_residues:
                buried_residues.append(xl_pair_1)
            if xl_pair_2 in remove_aa2 and not xl_pair_2 in buried_residues:
                buried_residues.append(xl_pair_2)
        
        elif [xl_pair_2[0],xl_pair_2[1]] in remove_aa2:        
            index_to_delete.append(i)
            if xl_pair_2 not in buried_residues:
                buried_residues.append(xl_pair_2)
    
    no_sasd_possible = []
    crosslink_pairs_final = []
    for i in range(len(crosslink_pairs)):
        if i not in index_to_delete:
            crosslink_pairs_final.append(crosslink_pairs[i])
        else:
            no_sasd_possible.append(crosslink_pairs[i])
    
    if len(no_sasd_possible) > 0:
        print("the following crosslinks cannot be calculated:")
        for s in no_sasd_possible:
            print("{}-{}-{} - {}-{}-{}".format(s[0][2],s[0][0],s[0][1],s[1][2],s[1][0],s[1][1]))
                    
    return crosslink_pairs_final
    
def check_solvent_accessibility_freesasa(prot, aa_CA, xl_list, aa_dict, ncpus):
    
    freesasa.Parameters().setNSlices(50)
    freesasa.Parameters().setNThreads(ncpus)

    pt = os.path.dirname(os.path.realpath(__file__))
    classifier = freesasa.Classifier(os.path.join(pt,"naccess.config.txt"))
    structure = freesasa.Structure(os.path.normpath(prot), classifier)
    result = freesasa.calc(structure)

    solv_access_residue = {}
    for chain, residue in result.residueAreas().items():
        for res_sasa_info in residue.values():
            if res_sasa_info.total > 7.0: # if total residue SASA is greater than 7.0 ... 
                solv_access_residue[(int(res_sasa_info.residueNumber), chain, res_sasa_info.residueType)] = True

    surface_solv_access_residue = {}

    for res_num, chain, res_name in aa_CA:
        if (res_num, chain, res_name) in solv_access_residue:
            surface_solv_access_residue[(res_num, chain, res_name)] = aa_CA[(res_num, chain, res_name)]
            sd_res = res_name
        else:
            print("Residue {}-{}-{} is buried".format(res_num, chain, res_name))
            sd_res = res_name
    
    # inform user on buried resiudes
    if xl_list != "NULL":
        pass
    elif sd_res == "LYS":
        print("{} {} and 1 N-terminus of which {} are on the surface".format(len(aa_CA)-1, aa_dict[sd_res], len(surface_solv_access_residue)))            
    else:            
        print("{} {} of which {} are on the surface".format(len(aa_CA), aa_dict[sd_res], len(surface_solv_access_residue)))
            
    return surface_solv_access_residue


def check_solvent_accessibility_freesasa_both(prot, aa1_CA, aa2_CA, xl_list, aa_dict, ncpus):
    """
    Run freesasa ONCE for the PDB and return filtered surface-accessible dicts
    for both aa1_CA and aa2_CA.  Avoids the duplicate freesasa call that occurs
    when the function is called separately for each residue set.

    Returns: (surface_aa1_CA, surface_aa2_CA)
    """
    freesasa.Parameters().setNSlices(50)
    freesasa.Parameters().setNThreads(ncpus)

    pt = os.path.dirname(os.path.realpath(__file__))
    classifier = freesasa.Classifier(os.path.join(pt, "naccess.config.txt"))
    structure = freesasa.Structure(os.path.normpath(prot), classifier)
    result = freesasa.calc(structure)

    solv_access_residue = {}
    for chain, residue in result.residueAreas().items():
        for res_sasa_info in residue.values():
            if res_sasa_info.total > 7.0:
                solv_access_residue[(int(res_sasa_info.residueNumber), chain, res_sasa_info.residueType)] = True

    def _filter(aa_CA):
        out = {}
        for res_num, chain, res_name in aa_CA:
            if (res_num, chain, res_name) in solv_access_residue:
                out[(res_num, chain, res_name)] = aa_CA[(res_num, chain, res_name)]
            else:
                print("Residue {}-{}-{} is buried".format(res_num, chain, res_name))
        if xl_list == "NULL" and out:
            sd_res = next(iter(out))[2]
            if sd_res == "LYS":
                print("{} {} and 1 N-terminus of which {} are on the surface".format(
                    len(aa_CA) - 1, aa_dict[sd_res], len(out)))
            else:
                print("{} {} of which {} are on the surface".format(
                    len(aa_CA), aa_dict[sd_res], len(out)))
        return out

    return _filter(aa1_CA), _filter(aa2_CA)
    