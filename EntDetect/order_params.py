#!/usr/bin/env python3
import requests, logging, os, sys
from EntDetect._logging import setup_logger
import time
import argparse
import pandas as pd
import numpy as np
import glob
import MDAnalysis as mda
import mdtraj as md
from scipy.spatial.distance import pdist, squareform
from topoly import lasso_type  # used pip
import itertools
import concurrent.futures
from EntDetect.gaussian_entanglement import GaussianEntanglement
from EntDetect.clustering import ClusterNativeEntanglements
from EntDetect.Jwalk import PDBTools, GridTools, SurfaceTools, SASDTools
from importlib.resources import files
import subprocess
import pathlib
from multiprocessing import cpu_count
from scipy.stats import norm

pd.set_option('display.max_rows', 5000)

class CalculateOP:
    """
    A class to handel the analyssis of a C-alpha CG trajectory. 
    Current analysis available:
    (1) - Fraction of native contacts (Q)
    (2) - Fraction of native contacts with a change in entanglement (G)
    (3) - Solvant Accessible Surface Area (SASA)
    (4) - Mirror symmetry order parameter (K)
    (5) - Cross linking probability score (XP)
    (6) - Jwalk SASD
    """
    #######################################################################################
    def __init__(self, outdir:str='./', ID:str='', Traj:int=1, psf:str='', cor:str='', dcd:str='', sec_elements:str='', domain:str='', start:int=0, end:int=99999999999999, stride:int=1, ent_detection_method:int=2, log_level:int=logging.INFO, logdir:str=None):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        ("--outdir", type=str, required=True, help="Path to output directory")
        ("--ID", type=str, required=True, help="base name for output files")
        ("--Traj", type=int, required=True, help="trajectory index")
        ("--psf", type=str, required=True, help="Path to CA protein structure file")
        ("--cor", type=str, required=True, help="Path to CA native coordinates file")
        ("--dcd", type=str, required=True, help="Path to trajectory to analyze")
        ("--sec_elements", type=str, required=True, help="Path to STRIDE secondary structure elements file")
        ("--domain", type=str, required=True, help="Path to domain definition file")
        ("--start", type=int, required=False, help="First frame to analyze 0 indexed", default=0)
        ("--end", type=int, required=False, help="Last frame to analyze 0 indexed", default=-1)
        ("--stride", type=int, required=False, help="Frame stride", default=1)
        """

        # parse the parameters
        self.logger = setup_logger('CalculateOP', outdir=logdir if logdir is not None else outdir, ID=ID, log_level=log_level)
        self.outdir = outdir
        self.logger.debug(f'outdir: {self.outdir}')

        self.ID = ID
        self.logger.debug(f'ID: {self.ID}')

        self.Traj = Traj
        self.logger.debug(f'Traj: {Traj}')

        self.psf = psf
        self.logger.debug(f'psf: {self. psf}')

        self.sec_elements = sec_elements
        self.logger.debug(f'sec_elements: {self.sec_elements}')

        self.domain = domain
        self.logger.debug(f'domain: {self.domain}')

        self.cor = cor
        self.logger.debug(f'cor: {self.cor}')

        self.dcd = dcd
        self.logger.debug(f'dcd: {self.dcd}')

        if self.cor != '' and self.cor.endswith('.cor'):
            self.ref_universe = mda.Universe(self.psf, self.cor, format='CRD')
            self.logger.debug(f'ref_universe: {self.ref_universe}')

        self.traj_universe = mda.Universe(self.psf, self.dcd, format='DCD')
        self.logger.debug(f'traj_universe: {self.traj_universe}')

        self.start = start
        self.end = end
        self.stride = stride
        self.logger.debug(f'START: {self.start} | END: {self.end} | STRIDE: {self.stride}')

        self.ent_detection_method = ent_detection_method
        self.logger.debug(f'ent_detection_method: {self.ent_detection_method}')

        self.three_to_one = {
                        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
                        "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
                        "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
                        "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
                        "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
                        "SEC": "U",  # Selenocysteine
                        "PYL": "O",  # Pyrrolysine
                        "ASX": "B",  # Asp or Asn (ambiguous)
                        "GLX": "Z",  # Glu or Gln (ambiguous)
                        "XAA": "X",  # Any/unknown amino acid
                        "TER": "*"} # Stop codon
    #######################################################################################

    #######################################################################################
    def Qpy(self, ):
        self.logger.info(f'Calculating the fraction of native contacts (Q)')
        """
        Calculate the fraction of native contacts in each frame of the DCD where a native contact is defined between secondary structures 
        and for residues atleast that are atleast 3 residues apart. So if i = 1 then j at a minimum can be 5. 
        For a contact to be present the distance between i and j must be less than 8A in the native structure and in a trajectory frame be less than 1.2*native distance.
        """
        # make directory for Q data if it doesnt exist
        self.Qpath = os.path.join(self.outdir, 'Q')
        if not os.path.exists(self.Qpath):
            os.makedirs(self.Qpath)
            self.logger.info(f'Made directory: {self.Qpath}')

        # Step 0: load the reference structure and topology
        ref_coor = self.ref_universe.atoms.positions
        #print(f'ref_coor:\n{ref_coor} {ref_coor.shape}')        


        # Step 1: Get the secondary structure information
        # get both those resides in the secondary structures and those not
        self.logger.info(f'Step 1: Get the secondary structure information')
        resid_in_sec_elements = np.loadtxt(self.sec_elements, dtype=int)
        resid_in_sec_elements = [np.arange(x[1], x[2] + 1) for x in resid_in_sec_elements]
        resid_in_sec_elements = np.hstack(resid_in_sec_elements)
        #print(f'resid_in_sec_elements: {resid_in_sec_elements}')

        resid_not_in_sec_elements = np.asarray([r for r in range(1, len(ref_coor) + 1) if r not in resid_in_sec_elements]) # residue ID not in secondary structures
        #print(f'resid_not_in_sec_elements: {resid_not_in_sec_elements}')


        # Step 2: Get the native distance map for the native state cordinates
        self.logger.info(f'Step 2: Get the native distance map for the native state cordinates')
        # Zero the resulting distance map up to the 4th diagonal so only those residues with more than 3 residues between them can be in contact
        # Zero out any secondary structure element residues
        # Zero out any distance not less than 8A
        ref_distances = np.triu(squareform(pdist(ref_coor)), k=4)
        ref_distances[resid_not_in_sec_elements - 1, :] = 0
        ref_distances[:, resid_not_in_sec_elements - 1] = 0
        ref_distances[ref_distances > 8] = 0
        NumNativeContacts = np.count_nonzero(ref_distances)
        self.logger.debug(f'NumNativeContacts: {NumNativeContacts}')
        self.logger.debug(f'NumNativeContacts: {NumNativeContacts}')

        # Step 3: Analyze each frame of the traj_universe and get the distance map
        self.logger.info(f'Step 3: Analyze each frame of the traj_universe and calc Q')
        # then determine the fraction of native contacts by those distances less than 1.2*native distance
        Qoutput = {'Time(ns)':[], 'Frame':[], 'FrameNumNativeContacts':[], 'Q':[]}
        for ts in self.traj_universe.trajectory[self.start:self.end:self.stride]:
            frame_coor = self.traj_universe.atoms.positions
            frame_distances = np.triu(squareform(pdist(frame_coor)), k=4)
            frame_distances[resid_not_in_sec_elements - 1, :] = 0
            frame_distances[:, resid_not_in_sec_elements - 1] = 0

            cond = (frame_distances <= 1.2*ref_distances) & (ref_distances != 0)

            FrameNumNativeContacts = np.sum(cond)
            #print(f'FrameNumNativeContacts: {FrameNumNativeContacts} for frame {ts.frame}')

            Q = FrameNumNativeContacts/NumNativeContacts
            #print(f'Q: {Q} for frame {ts.frame}')

            frame_time = ts.time/1000
            Qoutput['Frame'] += [ts.frame]
            Qoutput['FrameNumNativeContacts'] += [FrameNumNativeContacts]
            Qoutput['Q'] += [Q]
            Qoutput['Time(ns)'] += [frame_time]
        
        # Step 4: save Q output 
        self.logger.info(f'Step 4: save Q output')
        Qoutput = pd.DataFrame(Qoutput)
        Qoutfile = os.path.join(self.Qpath, f'{self.ID}.Q')
        Qoutput.to_csv(Qoutfile, index=False)
        self.logger.info(f'SAVED: {Qoutfile}')
        self.logger.info(f'SAVED: {Qoutfile}')
        return {'outfile':Qoutfile, 'result':Qoutput}
    #######################################################################################  

    #######################################################################################
    def Q(self,):
        """
        Calculate the fraction of native contacts (Q) using Yang's perl code which goes further and uses the domain definitions as well as the secondary structure elements
        it will return the fraction of native contacts overall (same as what Qpy) will give you as well as the Q within each domain and between them
        """
        # make directory for Q data if it doesnt exist
        self.Qpath = os.path.join(self.outdir, 'Q')
        if not os.path.exists(self.Qpath):
            os.makedirs(self.Qpath)
            self.logger.info(f'Made directory: {self.Qpath}')

        # Check if the Q output file exists. else make it
        dcdname = self.dcd.split('/')[-1].split('.')[0]
        self.logger.debug(f'dcdname: {dcdname}')
        outfilename = os.path.join(self.Qpath, f'Q_{dcdname}.dat')
        self.logger.debug(f'outfilename: {outfilename}')
        renamed_outfile = os.path.join(self.Qpath, f'{self.ID}_Traj{self.Traj}.Q') ## This is the new name of the calc_Q.pl output script after it has had the Frames added

        u = mda.Universe(self.psf, self.dcd)
        self.logger.debug(u)
        frames = [ts.frame for ts in u.trajectory]
        #print(f'frames: {frames}')
        if self.start < 0:
            self.start = frames[self.start]
        self.logger.debug(f'START: {self.start}')
        
        if os.path.exists(renamed_outfile):
            self.logger.info(f'Q outfile exists: {renamed_outfile}')
            Qoutput = pd.read_csv(renamed_outfile, sep = ',')
            #print(f'Qoutput:\n{Qoutput}')

        else:        
            script_path = files('EntDetect.resources').joinpath('calc_Q.pl')
            self.logger.debug(f'script_path: {script_path}')

            cmd = f'perl {script_path} -i {self.cor} -t {self.dcd} -d {self.domain} -s {self.sec_elements} -b {self.start + 1} -e {self.end} -o {self.Qpath}'
            self.logger.debug(f'cmd: {cmd}')
        
            result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
            if result.returncode != 0:
                raise RuntimeError(f"Perl script failed:\n{result.stderr}")
            self.logger.debug(result.stdout)

            ## rename the file to match the standard OP file format {ID}_Traj{traj}.Q
            os.rename(outfilename, renamed_outfile)
            self.logger.debug(f'Renamed: {outfilename} -> {renamed_outfile}')

            ## read the Q file back in and add the Frame column
            Qoutput = pd.read_csv(renamed_outfile, delim_whitespace=True)
            #print(f'Qoutput:\n{Qoutput}')
            
            #print(frames[self.start:self.end])
            sel_frames = frames[self.start:self.end]
            #print(f'sel_frames: {sel_frames}')
            Qoutput['Frame'] = sel_frames
            #print(f'Qoutput:\n{Qoutput}')

            Qoutput.to_csv(renamed_outfile, index=False, sep = ',')
            self.logger.info(f'SAVED: {renamed_outfile}')

        return {'outfile':renamed_outfile, 'result':Qoutput}
    #######################################################################################

    #######################################################################################
    def G(self, topoly:bool=True, Calpha:bool=True, CG:bool=True, nproc: int = 10, chunk_frames:int=None, chunk_suffix:str='_chunk') -> dict:
        self.logger.info(f'Calculating the G entanglement order parameter')
        """
        Calculate the G entanglement order parameter for each frame of the DCD 
        """
        # make directory for G data if it doesnt exist
        self.Gpath = os.path.join(self.outdir, 'G')
        if not os.path.exists(self.Gpath):
            os.makedirs(self.Gpath)
            self.logger.info(f'Made directory: {self.Gpath}')

        # parse some of the default parameters
        g_threshold = 0.6
        density = 1.0
        
        self.logger.debug(f'g_threshold: {g_threshold}')
        self.logger.debug(f'density: {density}')
        self.logger.debug(f'Calpha: {Calpha}')
        self.logger.debug(f'CG: {CG}')
        self.logger.debug(f'nproc: {nproc}')
        self.logger.debug(f'chunk_frames: {chunk_frames}')
        self.logger.debug(f'chunk_suffix: {chunk_suffix}')

        ## initialize the entanglement object
        ge = GaussianEntanglement(
            g_threshold=g_threshold,
            density=density,
            Calpha=Calpha,
            CG=CG,
            nproc=nproc,
            ent_detection_method=self.ent_detection_method,
        ) # for CG structures and trajectories
        #ge = GaussianEntanglement(g_threshold=g_threshold, density=density, Calpha=False, CG=False) # for all-atom structures
        self.logger.debug(ge)
    
        ## initialize the clustering object
        clustering = ClusterNativeEntanglements(organism='Ecoli')
        self.logger.debug(clustering)

        ## Get the native entanglements from a CG model
        self.logger.info(f'Calculating the native entanglements...')
        NativeEnt = ge.calculate_native_entanglements(self.cor, outdir=os.path.join(self.Gpath,'Native_GE/'), ID=f'{self.ID}_native', topoly=topoly)
        #print(NativeEnt)
        
        ## Cluster the native entanglements
        self.logger.info(f'Clustering the native entanglements...')
        nativeClusteredEnt = clustering.Cluster_NativeEntanglements(NativeEnt['outfile'], outdir=os.path.join(self.Gpath,'Native_clustered_GE/'), outfile=f'{self.ID}_NativeEntClusters.txt')
        #print(nativeClusteredEnt)
        
        ## Get the trajectory entanglements
        self.logger.info(f'Calculating the trajectory entanglements...')
        TrajEnt = ge.calculate_traj_entanglements(
            self.dcd,
            self.psf,
            outdir=os.path.join(self.Gpath, 'Traj_GE/'),
            ID=f'{self.ID}_traj{self.Traj}',
            start=self.start,
            stop=self.end,
            topoly=topoly,
            ref_contact_file=NativeEnt['outfile'],
        )
        #print(TrajEnt)
        
        ## Create the combined .pkl file required for clustering non-native entanglements
        ## Will also calculate G at the same time
        self.logger.info(f'Creating the combined .pkl file required for clustering non-native entanglements...')
        Combined_data = ge.combine_ref_traj_GE(NativeEnt['outfile'], TrajEnt['outfile'], outdir=os.path.join(self.Gpath,'Combined_GE/'), ID=f'{self.ID}_traj{self.Traj}', chunk_frames=chunk_frames, chunk_suffix=chunk_suffix)
        G = Combined_data['G']
        Goutfile = os.path.join(self.Gpath, f'{self.ID}_Traj{self.Traj}.G')
        G.to_csv(Goutfile, index=False)
        self.logger.info(f'SAVED: {Goutfile}')
        return {'outfile':Goutfile, 'result':G}
    #######################################################################################

    #######################################################################################
    def SASA(self,) -> dict:
        """
        Calculate the solvent accessible surface area (SASA) for each frame of the DCD
        """
        # make directory for SASA data if it doesnt exist
        self.SASAPATH = os.path.join(self.outdir, 'SASA')
        if not os.path.exists(self.SASAPATH):
            os.makedirs(self.SASAPATH)
            self.logger.info(f'Made directory: {self.SASAPATH}')

        # Step -1: get the resid list from the MDAnalysis universe self.traj_universe
        # this is the list of residues in the trajectory
        resids = self.traj_universe.atoms.residues.resids
        #print(f'resids: {resids}')
  
        # Step 0: load the dcd and psf into a mdtraj trajectory
        traj = md.load(self.dcd, top=self.psf)
        #print(f'traj: {traj}')

        # Get the total frames and then adjust the frame_number to start from there
        total_frames = len(traj)
        self.logger.debug(f'total_frames: {total_frames}')
        if self.start >= 0:
            frame_number = self.start
        else:
            frame_number = total_frames + self.start
        self.logger.debug(f'frame_number: {frame_number}')
    
        # Step 1: loop through the trajectory and calculate the SASA for each frame
        self.logger.info(f'Step 1: loop through the trajectory and calculate the SASA for each frame')

        SASAoutput = {'Time(ns)':[], 'Frame':[], 'resid':[], 'SASA(nm^2)':[]}
        for ts in traj[self.start:self.end:self.stride]:
            # calculate the SASA for the current frame
            sasa = md.shrake_rupley(ts, mode='residue', probe_radius=0.14, n_sphere_points=1000)
            #print(f'sasa: {sasa} {sasa.shape}')

            # get the time and frame number
            frame_time = ts.time[0]/1000
            #print(frame_time, frame_number)

            # add the results to the output dictionary
            for resididx, res_sasa in enumerate(sasa[0]):
                #print(res_sasa)
                SASAoutput['Time(ns)'] += [frame_time]
                SASAoutput['Frame'] += [frame_number]
                SASAoutput['resid'] += [resids[resididx]]
                SASAoutput['SASA(nm^2)'] += [res_sasa]
                #print(f'SASAoutput: {SASAoutput}')

            frame_number += 1
        
        # Step 2: save the SASA output
        self.logger.info(f'Step 2: save the SASA output')
        SASAoutput = pd.DataFrame(SASAoutput)
        self.logger.info(f'SASAoutput:\n{SASAoutput}')
        SASAoutfile = os.path.join(self.SASAPATH, f'{self.ID}_Traj{self.Traj}.SASA')
        SASAoutput.to_csv(SASAoutfile, index=False)
        self.logger.info(f'SAVED: {SASAoutfile}')

        return {'outfile':SASAoutfile, 'result':SASAoutput}
    #######################################################################################

    #######################################################################################
    def K(self,) -> dict:
        """
        Calculate the mirror symmetry order parameter K for each frame of the DCD
        """
        # make directory for SASA data if it doesnt exist
        self.KPATH = os.path.join(self.outdir, 'K')
        if not os.path.exists(self.KPATH):
            os.makedirs(self.KPATH)
            self.logger.info(f'Made directory: {self.KPATH}')

        script_path = files('EntDetect.resources').joinpath('calc_K.pl')
        #print(f'script_path: {script_path}')

        cmd = f'perl {script_path} -i {self.cor} -t {self.dcd} -d {self.domain} -s {self.sec_elements} -b {self.start} -e {self.end} -o {self.KPATH}'
        #print(f'cmd: {cmd}')
    
        result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
        if result.returncode != 0:
            raise RuntimeError(f"Perl script failed:\n{result.stderr}")
        #print(result)

        ## outfile will follow the following format K_{name}.dat where name is the name of the DCD read in
        dcdname = self.dcd.split('/')[-1].split('.')[0]
        #print(f'dcdname: {dcdname}')
        outfilename = os.path.join(self.KPATH, f'K_{dcdname}.dat')
        #print(f'outfilename: {outfilename}')

        if os.path.exists(outfilename):
            self.logger.info(f'K outfile exists: {outfilename}')
            Koutput = pd.read_csv(outfilename, delim_whitespace=True)
            self.logger.info(f'Koutput:\n{Koutput}')
            return {'outfile':outfilename, 'result':Koutput}
        else:
            self.logger.info(f'K outfile does not exist: {outfilename}')
            raise FileNotFoundError(f'K outfile does not exist: {outfilename}')
    #######################################################################################

    #######################################################################################
    def XP(self, pdb:str='None') -> dict:
        """
        Calculates the cross-linking probability (XP) for an all atom pdb for all pairs of the following amino acid types [K S, T, Y, M]
        """
        
        # make directory for Q data if it doesnt exist
        self.XPpath = os.path.join(self.outdir, 'XP')
        if not os.path.exists(self.XPpath):
            os.makedirs(self.XPpath)
            self.logger.info(f'Made directory: {self.XPpath}')

        # Make the xl_list file for all possible cross-linkable residue combinations
        xl_list = os.path.join(self.XPpath, f'{self.ID}_Traj{self.Traj}_XLresidue_pairs.txt')
        self.find_residue_pairs(pdb, output_file=xl_list)

        ## Check that the Jwalk results directory exists and load them
        #  else calculate the Jwalk distance and then load the files
        pdbObj = pathlib.Path(pdb)
        if not pdbObj.exists():
            self.logger.error(f'ERROR: The input file supplied cannot be found. Please enter a .pdb file type')
            sys.exit(2)
        jwalk_results_dir = os.path.join(self.XPpath, f'Jwalk_results_{self.ID}_Traj{self.Traj}')
        Jwalk_outfile = os.path.join(jwalk_results_dir, 'Jwalk_results', f'{pdbObj.stem}_crosslink_list.txt')
        if os.path.exists(Jwalk_outfile):
            self.logger.info(f'Jwalk outfile exists: {Jwalk_outfile}')

        else:       
            ## then calculate the Jwalk distance
            # pdb: Input path to .pdb file
            # xl_list: OPTIONAL - Input path to crosslink list (default: Finds all Lys-to-Lys crosslinks)
            # aa1: OPTIONAL - Specify inital crosslink amino acid three letter code (default: LYS)
            # aa2: OPTIONAL - Specify ending crosslink amino acid three letter code (default: LYS)
            # max_dist: OPTIONAL - Specify maximum crosslink distance cutoff in Angstroms (default: Keeps all distances)
            # jwalk_results_dir: OPTIONAL - Output path for Jwalk results (default: Out to "./Jwalk_results" in the current working directory)
            # vox: OPTIONAL - Specify voxel resolution to use in Angstrom (default: 1 Angstrom)
            # ncpus: OPTIONAL - Specify number of cpus to use (default: 1)

            #self.runJwalk(pdb, xl_list='NULL', aa1='LYS', aa2='LYS', max_dist=sys.float_info.max, jwalk_results_dir=self.XPpath, vox=1, ncpus=1) # No precompiled list of residues to check
            self.runJwalk(pdb, xl_list=xl_list, max_dist=sys.float_info.max, jwalk_results_dir=jwalk_results_dir, vox=1, ncpus=1)
            self.logger.debug('Jwalk calculated')

        ## Check that the Jwalk results directory exists and load them
        col_names = ["Index", "Model", "Atom1", "Atom2", "SASD", "Euclidean Distance"]
        Jwalk_df = pd.read_csv(Jwalk_outfile, sep=r'\s+', names=col_names, skiprows=1, engine='python', index_col=False)

        ## Calculate the XP score
        XP_scores = []
        for rowi, row in Jwalk_df.iterrows():
            AA1 = self.three_to_one[row['Atom1'].split('-')[0][0:3]]
            AA2 = self.three_to_one[row['Atom2'].split('-')[0][0:3]]
            pair_AA = (AA1, AA2)
            JWalk_dist = row['SASD']
            score = self.score_XL(pair_AA, JWalk_dist)
            XP_scores.append(score)

        Jwalk_df['XP'] = XP_scores
        # save the updated Jwalk file
        Jwalk_df.to_csv(Jwalk_outfile, index=False, sep='\t')
        self.logger.info(f'SAVED: {Jwalk_outfile}')
        return {'outfile':Jwalk_outfile, 'result':Jwalk_df}
    #######################################################################################

    #######################################################################################
    def find_residue_pairs(self, pdb_path, output_file="XLresidue_pairs.txt"):
        """
        Finds all unique residue pairs from amino acids [K, S, T, Y, M] in a PDB file.
        Writes output as: resnum1|chain1|resnum2|chain2|
        """
        u = mda.Universe(pdb_path)
        
        # Define one-letter code set and their three-letter equivalents
        aa_of_interest = {'LYS', 'SER', 'THR', 'TYR', 'MET'}
        
        # Select relevant residues
        selection = u.select_atoms("protein and (" + " or ".join(f"resname {aa}" for aa in aa_of_interest) + ")")
        residues = selection.residues
        
        # Create all unique, unordered pairs (no double-counting)
        pairs = list(itertools.combinations(residues, 2))
        
        full_pairs = {'resid1': [], 'resname1':[], 'chain1': [], 'resid2': [], 'resname2':[], 'chain2': []}
        with open(output_file, "w") as f:
            for res1, res2 in pairs:
                line = f"{res1.resid}|{res1.segid or res1.chain}|{res2.resid}|{res2.segid or res2.chain}|\n"
                f.write(line)

                full_pairs['resid1'] += [res1.resid]
                full_pairs['resname1'] += [res1.resname]
                full_pairs['chain1'] += [res1.segid or res1.chain]
                full_pairs['resid2'] += [res2.resid]
                full_pairs['resname2'] += [res2.resname]
                full_pairs['chain2'] += [res2.segid or res2.chain]

        # Convert to DataFrame
        full_pairs_df = pd.DataFrame(full_pairs)
        # Save to CSV   
        full_pairs_df.to_csv(output_file.replace('.txt', '_Full.csv'), index=False)
        self.logger.info(f"Residue pairs saved to '{output_file.replace('.txt', '_Full.csv')}'")
        
        self.logger.info(f"Found {len(pairs)} residue pairs and wrote to '{output_file}'")
    #######################################################################################

    #######################################################################################
    def score_XL(self, pair_AA, JWalk_dist, XL_offset:float=1.1):
        """
        Calculates the cross-linking probability score using the Jwalk distance and the amino acid types
        """
        sc_length = {'K': 6.3,
                    'S': 2.5,
                    'T': 2.5,
                    'Y': 6.5,
                    'M': 1.5,}
        
        KK_mu = 18.6
        KK_sigma = 6.0
        KK_threshold = 33

        KK_mu += XL_offset
        KK_sigma = (XL_offset + 3*KK_sigma) / 3
        KK_threshold += XL_offset

        mu = KK_mu + (sc_length[pair_AA[0]] + sc_length[pair_AA[1]]) - 2*sc_length['K']
        sigma = (mu - (KK_mu - 3*KK_sigma)) / 3
        threshold = KK_threshold + mu - KK_mu

        N = norm(mu, sigma)

        if JWalk_dist == -1:
            score = 0
        elif JWalk_dist <= threshold:
            score = N.pdf(JWalk_dist)
        else:
            score = 0
        return score
    #######################################################################################

    #######################################################################################
    def runJwalk(self, pdb, xl_list:str='NULL', aa1:str='LYS', aa2:str='LYS', max_dist:float=sys.float_info.max, jwalk_results_dir:str='./', vox:int=1, ncpus:int=1):
        """
            Execute Jwalk with processed command line options
                
            pdb: Input path to .pdb file
            xl_list: OPTIONAL - Input path to crosslink list (default: Finds all Lys-to-Lys crosslinks)
            aa1: OPTIONAL - Specify inital crosslink amino acid three letter code (default: LYS)
            aa2: OPTIONAL - Specify ending crosslink amino acid three letter code (default: LYS)
            max_dist: OPTIONAL - Specify maximum crosslink distance cutoff in Angstroms (default: Keeps all distances)
            jwalk_results_dir: OPTIONAL - Output path for Jwalk results (default: Out to "./Jwalk_results" in the current working directory)
            vox: OPTIONAL - Specify voxel resolution to use in Angstrom (default: 1 Angstrom)
            ncpus: OPTIONAL - Specify number of cpus to use (default: 1)            

            J.Bullock, J. Schwab, K. Thalassinos, M. Topf (2016)
                The importance of non-accessible crosslinks and solvent accessible surface distance
                in modelling proteins with restraints from crosslinking mass spectrometry. 
                Molecular and Cellular Proteomics (15) pp.2491–2500
        """
        self.logger.info("Running Jwalk with the following parameters:")
        self.logger.debug(f"pdb: {pdb}")
        self.logger.debug(f"xl_list: {xl_list}")
        self.logger.debug(f"aa1: {aa1}")
        self.logger.debug(f"aa2: {aa2}")
        self.logger.debug(f"max_dist: {max_dist}")
        self.logger.debug(f"jwalk_results_dir: {jwalk_results_dir}")
        self.logger.debug(f"vox: {vox}")
        self.logger.debug(f"ncpus: {ncpus}")

        # check if the number of cpus is greater than the number of available cpus
        max_cpus = cpu_count()
        amino_acids = {"LYS":"lysines",         "CYS":"cysteines",      "ASP":"aspartates",  "GLU":"glutamates",
                    "VAL":"valines",         "ILE":"isoleucines",    "LEU":"leucines",    "ARG":"arginines",
                    "PRO":"prolines",        "GLY":"glycines",       "ALA":"alanines",    "TRP":"tryptophans",
                    "PHE":"phenylalanines",  "SER":"serines",        "GLN":"glutamines",  "HIS":"histidines",
                    "MET":"methionines",     "THR":"threonines",     "ASN":"asparagines", "TYR":"tyrosines"}
        
        # checking if pdb file supplied exists and is of type .pdb
        if os.path.exists(pdb) and pdb.endswith(".pdb"):
            self.logger.info("PDB file supplied is valid")
            pass
        elif not os.path.exists(pdb):
            self.logger.error("ERROR: The input file supplied cannot be found. Please enter a .pdb file type")
            sys.exit(2)
        elif not pdb.endswith(".pdb"):
            self.logger.error("ERROR: The input file supplied is not supported. Please enter a .pdb file type")
            sys.exit(2)
        else:
            self.logger.error("ERROR: The input file supplied is not supported. Please enter a .pdb file type")
            sys.exit(2)

        # creating result output directory (defaulting to creating it in the working directory)
        if os.path.exists(jwalk_results_dir) and os.path.isdir(jwalk_results_dir):
            self.logger.warning(f"WARNING: {jwalk_results_dir} already exists. Overwriting directory")
            pass
        else:
            self.logger.warning(f"WARNING: {jwalk_results_dir} not found. Creating directory {jwalk_results_dir}")
            os.mkdir(jwalk_results_dir)
            pass
        
        # checking if an xl_list was provided
        # if none is provided use the aa1 and aa2 inputs (default is LYS-LYS crosslinks)
        if os.path.normpath(xl_list) == "NULL" or xl_list == "NULL":
            aa1 = aa1.upper()
            aa2 = aa2.upper()
            xl_list = "NULL"

            if aa1 not in amino_acids or aa2 not in amino_acids:
                self.logger.error("ERROR: Please type amino acid in three letter code format")
                self.logger.debug(amino_acids.keys())
                sys.exit(2)
            else:
                self.logger.info("Calculating all {}-to-{} crosslinks".format(aa1,aa2))
                pass
        # accepting xl_list
        elif os.path.exists(xl_list) and os.path.isfile(xl_list):
            self.logger.info(f"Calculating all crosslinks found in {xl_list}")
            aa1 = "NULL" 
            aa2 = "NULL"
            pass

        # load pdb into Jwalk
        structure_instance = PDBTools.read_PDB_file(pdb)
        # generate grid of voxel size (vox) that encapsulates pdb
        grid = GridTools.makeGrid(structure_instance, vox)

        # mark C-alpha positions on grid
        if xl_list != "NULL": # if specific crosslinks need to be calculated
            crosslink_pairs, aa1_CA, aa2_CA = GridTools.mark_CAlphas_pairs(grid, structure_instance, xl_list)
        else:
            crosslink_pairs = [] # na if searching every combination between residue types
            aa1_CA, aa2_CA = GridTools.mark_CAlphas(grid, structure_instance, aa1, aa2)
        
        # check more rigorously if residues are solvent accessible or not
        aa1_CA = SurfaceTools.check_solvent_accessibility_freesasa(pdb, aa1_CA, xl_list, amino_acids, ncpus)
        aa2_CA = SurfaceTools.check_solvent_accessibility_freesasa(pdb, aa2_CA, xl_list, amino_acids, ncpus)

        dens_map = GridTools.generate_solvent_accessible_surface(grid, structure_instance, aa1_CA, aa2_CA)    
        # identify which residues are on the surface
        aa1_voxels, remove_aa1 = GridTools.find_surface_voxels(aa1_CA, dens_map, xl_list)
        aa2_voxels, remove_aa2 = GridTools.find_surface_voxels(aa2_CA, dens_map, xl_list)
        
        crosslink_pairs = SurfaceTools.update_crosslink_pairs(crosslink_pairs, aa1_CA, aa2_CA, remove_aa1, remove_aa2)
        
        # calculate sasds
        sasds = SASDTools.parallel_BFS(aa1_voxels, aa2_voxels, dens_map, aa1_CA, aa2_CA, crosslink_pairs, max_dist, vox, ncpus, xl_list)

        # remove duplicates
        sasds = GridTools.remove_duplicates(sasds)
        sasds = SASDTools.get_euclidean_distances(sasds, pdb, aa1, aa2)
            
        # output sasds to .txt file (the .pdb visualisation file is skipped — the
        # chain-counter in write_sasd_to_pdb overflows for large residue-pair sets)
        PDBTools.write_sasd_to_txt(sasds, pdb,jwalk_results_dir)
        self.logger.info(f"{len(sasds)} SASDs calculated")
    #######################################################################################
     

## Round GaussLink values
def custom_round(number):
    if number >= 0:
        # For positive numbers, round up if fractional part >= 0.6
        return np.ceil(number) if number % 1 >= 0.6 else np.floor(number)
    else:
        # For negative numbers, round down if the absolute fractional part >= 0.6
        # need to take the abs of the number first else the modulus does work right for negative numbers?
        return np.floor(number) if abs(abs(number) % 1) >= 0.6 else np.ceil(number)

def process_frame(frame_data):
    frame_coor, nc_list, ref_nc_gdict, frame_time, frame, GaussLink, GetLinkChanges, Nnative = frame_data
    # Call GaussLink function
    t1 = time.time()
    frame_nc_gdict = GaussLink(frame_coor, contact_mask=nc_list)
    #print(f'FRAME: {frame} GaussLink time: {time.time() - t1}')

    # Call GetLinkChanges function
    t1 = time.time()
    change_info, count_info = GetLinkChanges(ref_nc_gdict, frame_nc_gdict, frame_time, frame, Nnative)
    #print(f'FRAME: {frame} GetLinkChanges time: {time.time() - t1}')
    return change_info, count_info
 

