#!/usr/bin/env python3
import requests, logging, os, sys
from NCLEdetector._logging import setup_logger
import time
import argparse
import pandas as pd
import numpy as np
import glob
import MDAnalysis as mda
import mdtraj as md
import freesasa
from scipy.spatial.distance import pdist, squareform
from topoly import lasso_type  # used pip
import itertools
import concurrent.futures
from NCLEdetector.gaussian_entanglement import GaussianEntanglement
from NCLEdetector.clustering import ClusterNativeEntanglements
from NCLEdetector.Jwalk import PDBTools, GridTools, SurfaceTools, SASDTools
from importlib.resources import files
import subprocess
import pathlib
import shutil
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
            use_traj=False (default): single static PDB
            use_traj=True:  per-frame from DCD, respects self.start/end/stride
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

        if self.cor is not None and self.cor != '' and self.cor.endswith('.cor'):
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
    def Qpy(self, chunk_frames:int=None, chunk_suffix:str='_chunk'):
        self.logger.info(f'Calculating the fraction of native contacts (Q)')
        """
        Calculate the fraction of native contacts in each frame of the DCD where a native contact is defined between secondary structures 
        and for residues atleast that are atleast 3 residues apart. So if i = 1 then j at a minimum can be 5. 
        For a contact to be present the distance between i and j must be less than 8A in the native structure and in a trajectory frame be less than 1.2*native distance.

        If chunk_frames is None (default): accumulates every frame's result in memory before a
        single write (backward compatible).
        If chunk_frames > 0: flushes accumulated rows to a restart-friendly part file every
        chunk_frames frames instead of holding the whole trajectory's results in memory. Parts
        are concatenated into the final {ID}.Q once complete, then removed. If interrupted and
        restarted, any part files already on disk are reused instead of recomputed.
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

        Qoutfile = os.path.join(self.Qpath, f'{self.ID}.Q')

        # Step 3: Analyze each frame of the traj_universe and get the distance map
        self.logger.info(f'Step 3: Analyze each frame of the traj_universe and calc Q')
        # then determine the fraction of native contacts by those distances less than 1.2*native distance
        def _calc_frame_Q(ts):
            frame_coor = self.traj_universe.atoms.positions
            frame_distances = np.triu(squareform(pdist(frame_coor)), k=4)
            frame_distances[resid_not_in_sec_elements - 1, :] = 0
            frame_distances[:, resid_not_in_sec_elements - 1] = 0

            cond = (frame_distances <= 1.2*ref_distances) & (ref_distances != 0)
            FrameNumNativeContacts = np.sum(cond)
            Q = FrameNumNativeContacts/NumNativeContacts
            frame_time = ts.time/1000
            return {'Time(ns)': frame_time, 'Frame': ts.frame, 'FrameNumNativeContacts': FrameNumNativeContacts, 'Q': Q}

        if chunk_frames is None:
            Qoutput = {'Time(ns)':[], 'Frame':[], 'FrameNumNativeContacts':[], 'Q':[]}
            for ts in self.traj_universe.trajectory[self.start:self.end:self.stride]:
                row = _calc_frame_Q(ts)
                for k, v in row.items():
                    Qoutput[k] += [v]

            # Step 4: save Q output
            self.logger.info(f'Step 4: save Q output')
            Qoutput = pd.DataFrame(Qoutput)
            Qoutput.to_csv(Qoutfile, index=False)
            self.logger.info(f'SAVED: {Qoutfile}')
            return {'outfile':Qoutfile, 'result':Qoutput}

        # Chunked mode: flush every chunk_frames frames to a restart-friendly part file
        part_dir = os.path.join(self.Qpath, f'.{self.ID}{chunk_suffix}_parts')
        os.makedirs(part_dir, exist_ok=True)

        part_files = []
        chunk_idx = 0
        rows_in_chunk = []

        def _flush(idx, rows):
            if not rows:
                return
            part_file = os.path.join(part_dir, f'{self.ID}{chunk_suffix}_{idx:04d}.Q')
            tmp_part_file = part_file + '.tmp'
            pd.DataFrame(rows).to_csv(tmp_part_file, index=False)
            os.replace(tmp_part_file, part_file)
            self.logger.info(f'SAVED chunk {idx}: {part_file}')
            part_files.append(part_file)

        frame_pos = 0
        for ts in self.traj_universe.trajectory[self.start:self.end:self.stride]:
            expected_part_file = os.path.join(part_dir, f'{self.ID}{chunk_suffix}_{chunk_idx:04d}.Q')
            if os.path.exists(expected_part_file):
                # Chunk already computed by a previous (interrupted) run; skip recompute for it
                if expected_part_file not in part_files:
                    part_files.append(expected_part_file)
                    self.logger.info(f'Chunk {chunk_idx} already computed, skipping: {expected_part_file}')
                frame_pos += 1
                if frame_pos % chunk_frames == 0:
                    chunk_idx += 1
                continue

            rows_in_chunk.append(_calc_frame_Q(ts))
            frame_pos += 1
            if frame_pos % chunk_frames == 0:
                _flush(chunk_idx, rows_in_chunk)
                rows_in_chunk = []
                chunk_idx += 1

        _flush(chunk_idx, rows_in_chunk)

        # Merge all chunk parts into the single final output file, one chunk in memory at a time
        for i, part_file in enumerate(part_files):
            part_df = pd.read_csv(part_file)
            part_df.to_csv(Qoutfile, mode='a' if i > 0 else 'w', header=(i == 0), index=False)
        self.logger.info(f'SAVED: {Qoutfile}')

        shutil.rmtree(part_dir, ignore_errors=True)
        Qoutput = pd.read_csv(Qoutfile)
        return {'outfile':Qoutfile, 'result':Qoutput}
    #######################################################################################  

    #######################################################################################
    def Q(self, chunk_frames:int=None, chunk_suffix:str='_chunk'):
        """
        Calculate the fraction of native contacts (Q) using Yang's perl code which goes further and uses the domain definitions as well as the secondary structure elements
        it will return the fraction of native contacts overall (same as what Qpy) will give you as well as the Q within each domain and between them

        If chunk_frames is None (default): a single perl invocation covers the whole
        [start, end) frame range (backward compatible).
        If chunk_frames > 0: the perl script is invoked once per sub-range of chunk_frames
        frames, bounding both perl's and Python's peak memory. Each chunk's output is written
        to a restart-friendly part file; parts are merged into the single final
        {ID}_Traj{N}.Q once complete, then removed. If interrupted and restarted, any part
        files already on disk are reused instead of recomputed.
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
            return {'outfile':renamed_outfile, 'result':Qoutput}

        script_path = files('NCLEdetector.resources').joinpath('calc_Q.pl')
        self.logger.debug(f'script_path: {script_path}')

        if chunk_frames is None:
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
            sel_frames = frames[self.start:self.end]
            Qoutput['Frame'] = sel_frames
            Qoutput.to_csv(renamed_outfile, index=False, sep = ',')
            self.logger.info(f'SAVED: {renamed_outfile}')
            return {'outfile':renamed_outfile, 'result':Qoutput}

        # Chunked mode: invoke the perl script over frame sub-ranges
        sel_frames = frames[self.start:self.end]
        total_frames = len(sel_frames)
        num_chunks = (total_frames + chunk_frames - 1) // chunk_frames
        part_dir = os.path.join(self.Qpath, f'.{self.ID}_Traj{self.Traj}{chunk_suffix}_parts')
        os.makedirs(part_dir, exist_ok=True)

        part_files = []
        for chunk_idx in range(num_chunks):
            start_idx = chunk_idx * chunk_frames
            end_idx = min(start_idx + chunk_frames, total_frames)
            chunk_frame_nums = sel_frames[start_idx:end_idx]
            part_file = os.path.join(part_dir, f'{self.ID}_Traj{self.Traj}{chunk_suffix}_{chunk_idx:04d}.Q')
            part_files.append(part_file)

            if os.path.exists(part_file):
                self.logger.info(f'Chunk {chunk_idx} already computed, skipping: {part_file}')
                continue

            chunk_tmp_dir = os.path.join(part_dir, f'tmp_{chunk_idx:04d}')
            os.makedirs(chunk_tmp_dir, exist_ok=True)
            chunk_b = self.start + start_idx + 1
            chunk_e = self.start + end_idx
            cmd = f'perl {script_path} -i {self.cor} -t {self.dcd} -d {self.domain} -s {self.sec_elements} -b {chunk_b} -e {chunk_e} -o {chunk_tmp_dir}'
            self.logger.debug(f'cmd (chunk {chunk_idx}): {cmd}')
            result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
            if result.returncode != 0:
                raise RuntimeError(f"Perl script failed on chunk {chunk_idx}:\n{result.stderr}")

            chunk_raw_file = os.path.join(chunk_tmp_dir, f'Q_{dcdname}.dat')
            chunk_df = pd.read_csv(chunk_raw_file, delim_whitespace=True)
            chunk_df['Frame'] = chunk_frame_nums
            tmp_part_file = part_file + '.tmp'
            chunk_df.to_csv(tmp_part_file, index=False, sep=',')
            os.replace(tmp_part_file, part_file)
            self.logger.info(f'SAVED chunk {chunk_idx}: {part_file}')
            shutil.rmtree(chunk_tmp_dir, ignore_errors=True)

        # Merge all chunk parts into the single final output file, one chunk in memory at a time
        for i, part_file in enumerate(part_files):
            chunk_df = pd.read_csv(part_file)
            chunk_df.to_csv(renamed_outfile, mode='a' if i > 0 else 'w', header=(i == 0), index=False)
        self.logger.info(f'SAVED: {renamed_outfile}')

        shutil.rmtree(part_dir, ignore_errors=True)
        Qoutput = pd.read_csv(renamed_outfile, sep=',')
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
        NativeEnt = ge.calculate_native_entanglements(self.cor, outdir=os.path.join(self.Gpath,'Native_NCLE/'), ID=f'{self.ID}_native', topoly=topoly)
        #print(NativeEnt)
        
        ## Cluster the native entanglements
        self.logger.info(f'Clustering the native entanglements...')
        nativeClusteredEnt = clustering.Cluster_NativeEntanglements(NativeEnt['outfile'], outdir=os.path.join(self.Gpath,'Native_clustered_NCLE/'), ID=f'{self.ID}_native')
        #print(nativeClusteredEnt)
        
        ## Get the trajectory entanglements
        self.logger.info(f'Calculating the trajectory entanglements...')
        TrajEnt = ge.calculate_traj_entanglements(
            self.dcd,
            self.psf,
            outdir=os.path.join(self.Gpath, 'Traj_NCLE/'),
            ID=f'{self.ID}_traj{self.Traj}',
            start=self.start,
            stop=self.end,
            topoly=topoly,
            ref_contact_file=NativeEnt['outfile'],
            chunk_frames=chunk_frames,
            chunk_suffix=chunk_suffix,
        )
        #print(TrajEnt)
        
        ## Create the combined .pkl file required for clustering non-native entanglements
        ## Will also calculate G at the same time
        self.logger.info(f'Creating the combined .pkl file required for clustering non-native entanglements...')
        Combined_data = ge.combine_ref_traj_NCLE(NativeEnt['outfile'], TrajEnt['outfile'], outdir=os.path.join(self.Gpath,'Combined_NCLE/'), ID=f'{self.ID}_traj{self.Traj}', chunk_frames=chunk_frames, chunk_suffix=chunk_suffix)
        G = Combined_data['G']
        Goutfile = os.path.join(self.Gpath, f'{self.ID}_Traj{self.Traj}.G')
        G.to_csv(Goutfile, index=False)
        self.logger.info(f'SAVED: {Goutfile}')
        return {'outfile':Goutfile, 'result':G}
    #######################################################################################

    #######################################################################################
    def SASA(self,) -> dict:
        """
        Calculate the solvent accessible surface area (SASA) for each frame of the DCD using freesasa.
        Uses freesasa library which is robust to coordinate artifacts (e.g., overlapping atoms).
        """
        import tempfile
        
        # make directory for SASA data if it doesnt exist
        self.SASAPATH = os.path.join(self.outdir, 'SASA')
        if not os.path.exists(self.SASAPATH):
            os.makedirs(self.SASAPATH)
            self.logger.info(f'Made directory: {self.SASAPATH}')

        # Step -1: get the resid list from the MDAnalysis universe self.traj_universe
        # this is the list of residues in the trajectory
        resids = self.traj_universe.atoms.residues.resids
  
        # Step 0: load the dcd and psf into a mdtraj trajectory for frame iteration
        traj = md.load(self.dcd, top=self.psf)

        # Get the total frames and then adjust the frame_number to start from there
        total_frames = len(traj)
        self.logger.debug(f'total_frames: {total_frames}')
        if self.start >= 0:
            frame_number = self.start
        else:
            frame_number = total_frames + self.start
        self.logger.debug(f'frame_number: {frame_number}')
    
        # Step 1: loop through the trajectory and calculate the SASA for each frame using freesasa
        self.logger.info(f'Step 1: loop through the trajectory and calculate the SASA for each frame using freesasa')

        SASAoutput = {'Time(ns)':[], 'Frame':[], 'resid':[], 'SASA(nm^2)':[]}
        last_valid_sasa = None  # Store last valid SASA results for fallback
        
        for ts in traj[self.start:self.end:self.stride]:
            # Save frame to temporary PDB
            with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
                tmp_pdb = tmp.name
            
            try:
                # Check for NaN coordinates before attempting to save PDB
                positions = ts.xyz
                if np.any(np.isnan(positions)):
                    nan_atoms = np.where(np.any(np.isnan(positions), axis=1))[0]
                    self.logger.warning(f'Frame {frame_number} has NaN coordinates in {len(nan_atoms)} atoms. Using SASA from previous frame.')
                    
                    # Use last valid SASA results if available
                    if last_valid_sasa is not None:
                        frame_time = ts.time[0]/1000
                        for resididx, res_sasa in enumerate(last_valid_sasa):
                            SASAoutput['Time(ns)'] += [frame_time]
                            SASAoutput['Frame'] += [frame_number]
                            SASAoutput['resid'] += [resids[resididx]]
                            SASAoutput['SASA(nm^2)'] += [res_sasa]
                    else:
                        self.logger.error(f'Frame {frame_number} has NaN coordinates but no previous valid SASA to fall back to. Skipping this frame.')
                    
                    frame_number += 1
                    continue
                
                ts.save_pdb(tmp_pdb)
                
                # Load with freesasa and calculate SASA
                structure = freesasa.Structure(tmp_pdb)
                result = freesasa.calc(structure)
                
                # Get per-residue SASA from freesasa result
                # residueAreas() returns nested dict: {chain: {res_num: ResidueArea_object}}
                res_areas = result.residueAreas()
                sasa_per_residue = []
                
                # Extract residues in order across all chains
                for chain in sorted(res_areas.keys()):
                    for res_num in sorted(res_areas[chain].keys(), key=lambda x: int(x)):
                        # Get total SASA for this residue (in Angstroms^2)
                        res_sasa = res_areas[chain][res_num].total
                        sasa_per_residue.append(res_sasa)
                
                # Convert from Angstroms^2 to nm^2 (1 nm^2 = 100 Angstroms^2)
                sasa_per_residue = np.array(sasa_per_residue) / 100.0
                
                # Store this as last valid SASA for fallback
                last_valid_sasa = sasa_per_residue
                
                # get the time
                frame_time = ts.time[0]/1000
                
                # add the results to the output dictionary
                for resididx, res_sasa in enumerate(sasa_per_residue):
                    SASAoutput['Time(ns)'] += [frame_time]
                    SASAoutput['Frame'] += [frame_number]
                    SASAoutput['resid'] += [resids[resididx]]
                    SASAoutput['SASA(nm^2)'] += [res_sasa]
                
                frame_number += 1
                
            finally:
                # Clean up temporary file
                if os.path.exists(tmp_pdb):
                    os.remove(tmp_pdb)
        
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
    def K(self, chunk_frames:int=None, chunk_suffix:str='_chunk') -> dict:
        """
        Calculate the mirror symmetry order parameter K for each frame of the DCD.

        Output is normalized to match Q()'s convention: a Frame column is added and the
        final file is renamed to {ID}_Traj{N}.K (the perl script's raw K_{dcdname}.dat is
        treated as an intermediate file, same as Q()'s handling of Q_{dcdname}.dat).

        If chunk_frames is None (default): a single perl invocation covers the whole
        [start, end) frame range (backward compatible aside from the renamed/Frame-tagged
        output described above).
        If chunk_frames > 0: the perl script is invoked once per sub-range of chunk_frames
        frames, bounding both perl's and Python's peak memory. Each chunk's output is written
        to a restart-friendly part file; parts are merged into the single final
        {ID}_Traj{N}.K once complete, then removed. If interrupted and restarted, any part
        files already on disk are reused instead of recomputed.
        """
        # make directory for K data if it doesnt exist
        self.KPATH = os.path.join(self.outdir, 'K')
        if not os.path.exists(self.KPATH):
            os.makedirs(self.KPATH)
            self.logger.info(f'Made directory: {self.KPATH}')

        dcdname = self.dcd.split('/')[-1].split('.')[0]
        outfilename = os.path.join(self.KPATH, f'K_{dcdname}.dat')
        renamed_outfile = os.path.join(self.KPATH, f'{self.ID}_Traj{self.Traj}.K')

        u = mda.Universe(self.psf, self.dcd)
        frames = [ts.frame for ts in u.trajectory]
        if self.start < 0:
            self.start = frames[self.start]

        if os.path.exists(renamed_outfile):
            self.logger.info(f'K outfile exists: {renamed_outfile}')
            Koutput = pd.read_csv(renamed_outfile, sep=',')
            return {'outfile':renamed_outfile, 'result':Koutput}

        script_path = files('NCLEdetector.resources').joinpath('calc_K.pl')

        if chunk_frames is None:
            cmd = f'perl {script_path} -i {self.cor} -t {self.dcd} -d {self.domain} -s {self.sec_elements} -b {self.start + 1} -e {self.end} -o {self.KPATH}'
            result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
            if result.returncode != 0:
                raise RuntimeError(f"Perl script failed:\n{result.stderr}")

            if not os.path.exists(outfilename):
                raise FileNotFoundError(f'K outfile does not exist: {outfilename}')

            Koutput = pd.read_csv(outfilename, delim_whitespace=True)
            sel_frames = frames[self.start:self.end]
            Koutput['Frame'] = sel_frames
            Koutput.to_csv(renamed_outfile, index=False, sep=',')
            self.logger.info(f'SAVED: {renamed_outfile}')
            return {'outfile':renamed_outfile, 'result':Koutput}

        # Chunked mode: invoke the perl script over frame sub-ranges
        sel_frames = frames[self.start:self.end]
        total_frames = len(sel_frames)
        num_chunks = (total_frames + chunk_frames - 1) // chunk_frames
        part_dir = os.path.join(self.KPATH, f'.{self.ID}_Traj{self.Traj}{chunk_suffix}_parts')
        os.makedirs(part_dir, exist_ok=True)

        part_files = []
        for chunk_idx in range(num_chunks):
            start_idx = chunk_idx * chunk_frames
            end_idx = min(start_idx + chunk_frames, total_frames)
            chunk_frame_nums = sel_frames[start_idx:end_idx]
            part_file = os.path.join(part_dir, f'{self.ID}_Traj{self.Traj}{chunk_suffix}_{chunk_idx:04d}.K')
            part_files.append(part_file)

            if os.path.exists(part_file):
                self.logger.info(f'Chunk {chunk_idx} already computed, skipping: {part_file}')
                continue

            chunk_tmp_dir = os.path.join(part_dir, f'tmp_{chunk_idx:04d}')
            os.makedirs(chunk_tmp_dir, exist_ok=True)
            chunk_b = self.start + start_idx + 1
            chunk_e = self.start + end_idx
            cmd = f'perl {script_path} -i {self.cor} -t {self.dcd} -d {self.domain} -s {self.sec_elements} -b {chunk_b} -e {chunk_e} -o {chunk_tmp_dir}'
            result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
            if result.returncode != 0:
                raise RuntimeError(f"Perl script failed on chunk {chunk_idx}:\n{result.stderr}")

            chunk_raw_file = os.path.join(chunk_tmp_dir, f'K_{dcdname}.dat')
            if not os.path.exists(chunk_raw_file):
                raise FileNotFoundError(f'K outfile does not exist: {chunk_raw_file}')
            chunk_df = pd.read_csv(chunk_raw_file, delim_whitespace=True)
            chunk_df['Frame'] = chunk_frame_nums
            tmp_part_file = part_file + '.tmp'
            chunk_df.to_csv(tmp_part_file, index=False, sep=',')
            os.replace(tmp_part_file, part_file)
            self.logger.info(f'SAVED chunk {chunk_idx}: {part_file}')
            shutil.rmtree(chunk_tmp_dir, ignore_errors=True)

        # Merge all chunk parts into the single final output file, one chunk in memory at a time
        for i, part_file in enumerate(part_files):
            chunk_df = pd.read_csv(part_file)
            chunk_df.to_csv(renamed_outfile, mode='a' if i > 0 else 'w', header=(i == 0), index=False)
        self.logger.info(f'SAVED: {renamed_outfile}')

        shutil.rmtree(part_dir, ignore_errors=True)
        Koutput = pd.read_csv(renamed_outfile, sep=',')
        return {'outfile':renamed_outfile, 'result':Koutput}
    #######################################################################################

    #######################################################################################
    def XP(self, pdb_file:str='None', use_traj:bool=False, nproc:int=1, **kwargs) -> dict:
        """
        Calculates the cross-linking probability (XP) for all pairs of amino acid types [K, S, T, Y, M].

        use_traj=False (default):
            Runs on the single static PDB supplied as `pdb_file`.
            Output: XP/Jwalk_results_{ID}_Traj{N}/Jwalk_results/{stem}_crosslink_list.txt

        use_traj=True:
            Iterates over DCD frames [self.start : self.end : self.stride]. For each frame a
            temporary per-frame PDB is written, Jwalk is run, XP is scored, and the per-frame PDB
            is deleted immediately. Produces a single combined tab-separated file:
            XP/{ID}_Traj{N}.XP  with columns: Frame | Index | Model | Atom1 | Atom2 | SASD | Euclidean Distance | XP
            nproc > 1 parallelises frame-level Jwalk runs via ThreadPoolExecutor.
        """
        legacy_pdb = kwargs.pop('pdb', None)
        if kwargs:
            unexpected = ', '.join(sorted(kwargs.keys()))
            raise TypeError(f"Unexpected keyword argument(s): {unexpected}")
        if legacy_pdb is not None and pdb_file == 'None':
            pdb_file = legacy_pdb

        traj_nproc = 1 # number of processors in the pool for the trajectory mode — Jwalk is not thread safe so must be run with nproc=1, but we can parallelise across frames with ThreadPoolExecutor

        # make output directory
        self.XPpath = os.path.join(self.outdir, 'XP')
        if not os.path.exists(self.XPpath):
            os.makedirs(self.XPpath)
            self.logger.info(f'Made directory: {self.XPpath}')

        col_names = ["Index", "Model", "Atom1", "Atom2", "SASD", "Euclidean Distance"]

        if not use_traj:
            # ── single-PDB path (original behaviour, unchanged) ───────────────
            xl_list = os.path.join(self.XPpath, f'{self.ID}_Traj{self.Traj}_XLresidue_pairs.txt')
            self.find_residue_pairs(pdb_file, output_file=xl_list)

            pdbObj = pathlib.Path(pdb_file)
            if not pdbObj.exists():
                self.logger.error(f'ERROR: The input file supplied cannot be found. Please enter a .pdb file type')
                sys.exit(2)
            jwalk_results_dir = os.path.join(self.XPpath, f'Jwalk_results_{self.ID}_Traj{self.Traj}')
            Jwalk_outfile = os.path.join(jwalk_results_dir, 'Jwalk_results', f'{pdbObj.stem}_crosslink_list.txt')
            if os.path.exists(Jwalk_outfile):
                self.logger.info(f'Jwalk outfile exists: {Jwalk_outfile}')
            else:
                self.runJwalk(pdb_file, xl_list=xl_list, max_dist=50.0,
                              jwalk_results_dir=jwalk_results_dir, vox=1, ncpus=nproc)
                self.logger.debug('Jwalk calculated')

            Jwalk_df = pd.read_csv(Jwalk_outfile, sep=r'\s+', names=col_names,
                                   skiprows=1, engine='python', index_col=False)
            XP_scores = []
            for rowi, row in Jwalk_df.iterrows():
                AA1 = self.three_to_one[row['Atom1'].split('-')[0][0:3]]
                AA2 = self.three_to_one[row['Atom2'].split('-')[0][0:3]]
                XP_scores.append(self.score_XL((AA1, AA2), row['SASD']))
            Jwalk_df['XP'] = XP_scores
            Jwalk_df.to_csv(Jwalk_outfile, index=False, sep='\t')
            self.logger.info(f'SAVED: {Jwalk_outfile}')
            return {'outfile': Jwalk_outfile, 'result': Jwalk_df}

        else:
            # ── trajectory mode ───────────────────────────────────────────────

            # skip-if-exists guard
            combined_outfile = os.path.join(self.XPpath, f'{self.ID}_Traj{self.Traj}.XP')
            if os.path.exists(combined_outfile):
                self.logger.info(f'XP outfile exists, loading: {combined_outfile}')
                return {'outfile': combined_outfile,
                        'result': pd.read_csv(combined_outfile, sep='\t')}

            # compute residue pairs once from the reference PDB (topology is frame-invariant)
            xl_list = os.path.join(self.XPpath, f'{self.ID}_Traj{self.Traj}_XLresidue_pairs.txt')
            self.find_residue_pairs(pdb_file, output_file=xl_list)

            # temporary directory for per-frame PDB files
            frames_dir = os.path.join(self.XPpath, f'frames_Traj{self.Traj}')
            os.makedirs(frames_dir, exist_ok=True)

            # parent directory for per-frame Jwalk outputs
            # runJwalk uses os.mkdir so the parent must already exist
            jwalk_base_dir = os.path.join(self.XPpath, f'Jwalk_results_{self.ID}_Traj{self.Traj}')
            os.makedirs(jwalk_base_dir, exist_ok=True)

            # per-frame worker — closure over self, safe for ThreadPoolExecutor
            def _run_frame(task):
                frame_idx, frame_pdb, frame_jwalk_dir = task
                pdb_stem = pathlib.Path(frame_pdb).stem
                jwalk_outfile = os.path.join(frame_jwalk_dir, 'Jwalk_results',
                                             f'{pdb_stem}_crosslink_list.txt')
                print(f'\nProcessing frame {frame_idx} | PDB: {frame_pdb} | Jwalk out: {jwalk_outfile}')

                if not os.path.exists(jwalk_outfile):
                    self.runJwalk(frame_pdb, xl_list=xl_list, max_dist=50.0,
                                  jwalk_results_dir=frame_jwalk_dir, vox=1, ncpus=nproc)
                    print(f'Jwalk completed for frame {frame_idx}')

                frame_df = pd.read_csv(jwalk_outfile, sep=r'\s+', names=col_names,
                                       skiprows=1, engine='python', index_col=False)
                print(f'Jwalk results loaded for frame {frame_idx}, calculating XP...')

                xp_scores = [
                    self.score_XL(
                        (self.three_to_one[row['Atom1'].split('-')[0][0:3]],
                         self.three_to_one[row['Atom2'].split('-')[0][0:3]]),
                        row['SASD']
                    )
                    for _, row in frame_df.iterrows()
                ]
                frame_df['XP'] = xp_scores
                frame_df['Frame'] = frame_idx
                # delete per-frame PDB immediately after use
                if os.path.exists(frame_pdb):
                    os.remove(frame_pdb)
                self.logger.debug(f'Frame {frame_idx}: XP computed, per-frame PDB removed')
                return frame_df

            frame_tasks = []
            results = []
            last_valid_frame_result = None  # Store last valid frame results for fallback
            
            for ts in self.traj_universe.trajectory[self.start:self.end:self.stride]:
                frame_idx = ts.frame
                
                # Check for NaN coordinates before attempting to write PDB
                positions = self.traj_universe.atoms.positions
                if np.any(np.isnan(positions)):
                    nan_atoms = np.where(np.any(np.isnan(positions), axis=1))[0]
                    self.logger.warning(f'Frame {frame_idx} has NaN coordinates in {len(nan_atoms)} atoms. Using XP from previous frame.')
                    
                    # Use last valid frame results if available
                    if last_valid_frame_result is not None:
                        # Copy previous frame's results but update frame number
                        fallback_frame_result = last_valid_frame_result.copy()
                        fallback_frame_result['Frame'] = frame_idx
                        results.append(fallback_frame_result)
                        self.logger.debug(f'Frame {frame_idx}: Used fallback XP from previous frame')
                    else:
                        self.logger.error(f'Frame {frame_idx} has NaN coordinates but no previous valid XP to fall back to. Skipping this frame.')
                    
                    continue
                
                frame_pdb = os.path.join(frames_dir, f'frame_{frame_idx}.pdb')
                with mda.Writer(frame_pdb, self.traj_universe.atoms.n_atoms) as W:
                    W.write(self.traj_universe.atoms)


                if traj_nproc > 1:
                    # parallel: write all frame PDBs first, then process with ThreadPoolExecutor
                    # (trajectory iteration must be sequential; Jwalk runs are independent)
                    frame_tasks.append((frame_idx, frame_pdb, os.path.join(jwalk_base_dir, f'frame_{frame_idx}')))

                else:
                    # sequential: write PDB → run Jwalk → score → delete, one frame at a time
                    frame_df = _run_frame((frame_idx, frame_pdb, os.path.join(jwalk_base_dir, f'frame_{frame_idx}')))
                    results.append(frame_df)
                    last_valid_frame_result = frame_df  # Update fallback with this valid result
                        
            if traj_nproc > 1:
                print(f'\nRunning frame-level Jwalk in parallel with {traj_nproc} workers...')
                with concurrent.futures.ThreadPoolExecutor(max_workers=traj_nproc) as executor:
                    results = list(executor.map(_run_frame, frame_tasks))
                # Update fallback with last result from parallel execution
                if results:
                    last_valid_frame_result = results[-1]
            else:
                print(f'Jwalk run completed for all frames in sequential mode.')


            combined_df = pd.concat(results, ignore_index=True)
            combined_df = combined_df[['Frame'] + [c for c in combined_df.columns if c != 'Frame']]
            combined_df.to_csv(combined_outfile, index=False, sep='\t')
            self.logger.info(f'SAVED: {combined_outfile}')
            try:
                if os.path.isdir(jwalk_base_dir):
                    shutil.rmtree(jwalk_base_dir)
                    self.logger.info(f'Removed temporary Jwalk directory: {jwalk_base_dir}')
            except Exception as exc:
                self.logger.warning(f'Could not remove temporary Jwalk directory {jwalk_base_dir}: {exc}')
            return {'outfile': combined_outfile, 'result': combined_df}
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
    def runJwalk(self, pdb, xl_list:str='NULL', aa1:str='LYS', aa2:str='LYS', max_dist:float=50.0, jwalk_results_dir:str='./', vox:int=1, ncpus:int=1):
        """
            Execute Jwalk with processed command line options
                
            pdb: Input path to .pdb file
            xl_list: OPTIONAL - Input path to crosslink list (default: Finds all Lys-to-Lys crosslinks)
            aa1: OPTIONAL - Specify inital crosslink amino acid three letter code (default: LYS)
            aa2: OPTIONAL - Specify ending crosslink amino acid three letter code (default: LYS)
            max_dist: OPTIONAL - Specify maximum crosslink distance cutoff in Angstroms (default: 50.0 Angstroms)
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
        print(f"Generating grid for PDB: {pdb}")
        grid = GridTools.makeGrid(structure_instance, vox)

        # mark C-alpha positions on grid
        print(f"Marking C-alpha positions on grid for PDB: {pdb}")
        if xl_list != "NULL": # if specific crosslinks need to be calculated
            crosslink_pairs, aa1_CA, aa2_CA = GridTools.mark_CAlphas_pairs(grid, structure_instance, xl_list)
        else:
            crosslink_pairs = [] # na if searching every combination between residue types
            aa1_CA, aa2_CA = GridTools.mark_CAlphas(grid, structure_instance, aa1, aa2)
        
        # check more rigorously if residues are solvent accessible or not
        print(f"Checking solvent accessibility for C-alpha positions for PDB: {pdb}")
        aa1_CA, aa2_CA = SurfaceTools.check_solvent_accessibility_freesasa_both(
            pdb, aa1_CA, aa2_CA, xl_list, amino_acids, ncpus
        )

        dens_map = GridTools.generate_solvent_accessible_surface(grid, structure_instance, aa1_CA, aa2_CA)    
        # identify which residues are on the surface
        aa1_voxels, remove_aa1 = GridTools.find_surface_voxels(aa1_CA, dens_map, xl_list)
        aa2_voxels, remove_aa2 = GridTools.find_surface_voxels(aa2_CA, dens_map, xl_list)
        
        crosslink_pairs = SurfaceTools.update_crosslink_pairs(crosslink_pairs, aa1_CA, aa2_CA, remove_aa1, remove_aa2)
        
        # calculate sasds
        print(f"Calculating SASDs for PDB: {pdb} with len(crosslink_pairs): {len(crosslink_pairs)}")
        sasds = SASDTools.parallel_BFS(aa1_voxels, aa2_voxels, dens_map, aa1_CA, aa2_CA, crosslink_pairs, max_dist, vox, ncpus, xl_list)

        # remove duplicates
        print(f"Removing duplicate SASDs for PDB: {pdb}")
        sasds = GridTools.remove_duplicates(sasds)
        sasds = SASDTools.get_euclidean_distances(sasds, pdb, aa1, aa2)
            
        # output sasds to .txt file (the .pdb visualisation file is skipped — the
        # chain-counter in write_sasd_to_pdb overflows for large residue-pair sets)
        PDBTools.write_sasd_to_txt(sasds, pdb,jwalk_results_dir)
        self.logger.info(f"{len(sasds)} SASDs calculated")
        print(f"Jwalk completed for PDB: {pdb} | Results saved to: {jwalk_results_dir}")
    #######################################################################################
     

#########################################################################################
class CollectOP:
    """
    Aggregate per-trajectory CalculateOP outputs into the single .npy arrays
    expected by MassSpec (compare_sim2exp).

    Reads
    -----
    {sasa_dir}/{ID}_Traj{N}.SASA   – CSV written by CalculateOP.SASA()
    {xp_dir}/{ID}_Traj{N}.XP       – TSV written by CalculateOP.XP()

    Writes
    ------
    SASA.npy  : float64 array  (n_traj, sasa_xp_frames_per_traj, prot_len)   units Å²
    Jwalk.npy : object  array  (n_traj, sasa_xp_frames_per_traj)
                each element is a dict
                { 'RESNUM|CHAIN-RESNUM|CHAIN' : {'Euclidean': float, 'Jwalk': float} }

    Trajectories whose output file is missing are filled with NaN (SASA) or
    left as None (Jwalk) so that MassSpec can skip them via its existing NaN
    filtering logic.
    """

    def __init__(self, sasa_dir: str, xp_dir: str, outdir: str, ID: str,
                 n_traj: int, sasa_xp_frames_per_traj: int, prot_len: int):
        """
        Parameters
        ----------
        sasa_dir  : directory containing {ID}_Traj{N}.SASA files
        xp_dir    : directory containing {ID}_Traj{N}.XP   files
        outdir    : directory where SASA.npy and Jwalk.npy are written
        ID        : protein ID used in file naming (e.g. '1ZMR')
        n_traj    : total number of trajectories (files named 1 … n_traj)
        sasa_xp_frames_per_traj : number of frames per trajectory stored in each SASA/XP file
        prot_len  : number of residues in the protein
        """
        self.sasa_dir = sasa_dir
        self.xp_dir   = xp_dir
        self.outdir   = outdir
        self.ID       = ID
        self.n_traj   = n_traj
        self.sasa_xp_frames_per_traj = sasa_xp_frames_per_traj
        self.prot_len = prot_len
        os.makedirs(outdir, exist_ok=True)
        self.logger   = setup_logger('CollectOP', outdir)

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _atom_to_key_part(atom_str: str) -> str:
        """Parse Jwalk atom string to key fragment.

        'MET-1-A-CA'  →  '1|A'
        """
        parts = atom_str.split('-')
        return f'{parts[1]}|{parts[2]}'

    # ------------------------------------------------------------------
    def collect_SASA(self, outfile: str = 'SASA.npy') -> str:
        """Read all {ID}_Traj{N}.SASA CSVs, convert nm² → Å² (×100), pivot
        each to (sasa_xp_frames_per_traj, prot_len), stack into
        (n_traj, sasa_xp_frames_per_traj, prot_len) and save.  Missing
        trajectory files are filled with NaN.

        Returns the absolute path to the saved .npy file.
        """
        out_path = os.path.join(self.outdir, outfile)
        sasa_arr = np.full(
            (self.n_traj, self.sasa_xp_frames_per_traj, self.prot_len),
            np.nan,
            dtype=np.float64,
        )

        for traj_num in range(1, self.n_traj + 1):
            fpath = os.path.join(self.sasa_dir, f'{self.ID}_Traj{traj_num}.SASA')
            if not os.path.isfile(fpath):
                self.logger.warning(f'Missing SASA file: {fpath}')
                continue

            df = pd.read_csv(fpath)

            # pivot to (sasa_xp_frames_per_traj, prot_len): rows = frames (sorted), cols = resids (sorted)
            pivot = (
                df.pivot_table(index='Frame', columns='resid',
                               values='SASA(nm^2)', aggfunc='first')
                  .sort_index()        # ascending frame order
                  .sort_index(axis=1)  # ascending resid order
            )
            arr = pivot.values  # shape (sasa_xp_frames_per_traj, prot_len)

            if arr.shape != (self.sasa_xp_frames_per_traj, self.prot_len):
                self.logger.warning(
                    f'Traj {traj_num}: unexpected shape {arr.shape}, '
                    f'expected ({self.sasa_xp_frames_per_traj}, {self.prot_len}) – skipping'
                )
                continue

            sasa_arr[traj_num - 1] = arr * 100.0  # nm² → Å²
            self.logger.info(f'Collected SASA: Traj {traj_num}')

        np.save(out_path, sasa_arr)
        self.logger.info(f'SAVED: {out_path}  shape={sasa_arr.shape}')
        return out_path

    # ------------------------------------------------------------------
    def collect_Jwalk(self, outfile: str = 'Jwalk.npy') -> str:
        """Read all {ID}_Traj{N}.XP TSVs and reconstruct the per-frame dict
        structure used by MassSpec.  Save an object array of shape
        (n_traj, sasa_xp_frames_per_traj).

        Each array element is a dict::

            { 'RESNUM|CHAIN-RESNUM|CHAIN' : {'Euclidean': float, 'Jwalk': float} }

        The 'SASD' column from the XP file maps to the 'Jwalk' key; the
        'Euclidean Distance' column maps to 'Euclidean'.
        Missing trajectory files leave the corresponding row as None entries.

        Returns the absolute path to the saved .npy file.
        """
        out_path  = os.path.join(self.outdir, outfile)
        jwalk_arr = np.empty((self.n_traj, self.sasa_xp_frames_per_traj), dtype=object)

        for traj_num in range(1, self.n_traj + 1):
            fpath = os.path.join(self.xp_dir, f'{self.ID}_Traj{traj_num}.XP')
            if not os.path.isfile(fpath):
                self.logger.warning(f'Missing XP file: {fpath}')
                continue

            df = pd.read_csv(
                fpath,
                sep='\t',
                usecols=['Frame', 'Atom1', 'Atom2', 'Euclidean Distance', 'SASD'],
                dtype={
                    'Frame': np.int32,
                    'Euclidean Distance': np.float32,
                    'SASD': np.float32,
                    'Atom1': 'string',
                    'Atom2': 'string',
                },
            )
            frames = sorted(df['Frame'].unique())

            if len(frames) != self.sasa_xp_frames_per_traj:
                self.logger.warning(
                    f'Traj {traj_num}: found {len(frames)} frames, '
                    f'expected {self.sasa_xp_frames_per_traj} – skipping'
                )
                continue

            for frame_idx, (_, fdf) in enumerate(df.groupby('Frame', sort=True)):
                fdict = {}
                for _, row in fdf.iterrows():
                    k1  = self._atom_to_key_part(row['Atom1'])
                    k2  = self._atom_to_key_part(row['Atom2'])
                    fdict[f'{k1}-{k2}'] = {
                        'Euclidean': float(row['Euclidean Distance']),
                        'Jwalk':     float(row['SASD']),
                    }
                jwalk_arr[traj_num - 1, frame_idx] = fdict

            del df

            self.logger.info(f'Collected Jwalk/XP: Traj {traj_num}')

        np.save(out_path, jwalk_arr, allow_pickle=True)
        self.logger.info(f'SAVED: {out_path}  shape={jwalk_arr.shape}')
        return out_path
#########################################################################################


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
 

