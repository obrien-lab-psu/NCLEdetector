#!/usr/bin/env python3
try:
    from openmm.app import *
    from openmm import *
    from openmm.unit import *
except:
    from simtk.openmm.app import *
    from simtk.openmm import *
    from simtk.unit import *
from sys import stdout, exit, stderr
import getopt, os, time, random, math, traceback, io, sys, string
import parmed as pmd
import numpy as np
from importlib.resources import files
import xml.etree.cElementTree as ET
import xml.dom.minidom as MD
import numpy
import pathlib
import subprocess
import logging
from NCLEdetector._logging import setup_logger

sys.setrecursionlimit(int(1e6))

class CoarseGrain:
    """
    Processes biological data including PDB files, sequence data, and interaction potentials.
    """
    #############################################################################################################
    def __init__(self, pdbfile:str, ID:str='ID', nscal:int = 1.5, outdir:str = './', fnn:int = 1,
                 potential_name:str = 'bt', casm:int = 0, domain_file:str = 'None', ca_prefix:str = 'A', sc_prefix:str = 'B', log_level:int = logging.INFO, logdir:str = None):
        
        self.pdbfile = pdbfile
        self.ID = ID
        self.nscal = nscal
        self.outdir = outdir
        self.logger = setup_logger('CoarseGrain', outdir=logdir if logdir is not None else outdir, ID=ID, log_level=log_level)
        self.fnn = fnn
        self.potential_name = potential_name
        self.casm = casm
        self.domain_file = domain_file
        self.ca_prefix = ca_prefix
        self.sc_prefix = sc_prefix
        self.heav_cut = 4.5

        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
            self.logger.info(f'Made directory: {self.outdir}')

        ######################## Data #########################
        ## Loop-up table for uniquely indentifying residues #
        self.aa = ["GLY","ALA","VAL","LEU","ILE","MET","PHE","PRO","SER","THR","CYS","ASN","GLN","TYR","TRP","ASP","GLU","HIS","LYS","ARG"]
        if len(self.aa) != 20:
            self.logger.error('ERROR')
            sys.exit()
        res2n = {}
        n2res = {}
        for i, a in enumerate(self.aa):
            res2n[a] = i
            n2res[i] = a
        self.res2n = res2n
        self.n2res = n2res

        self.Mass = {"N": 14.0067,
                "H": 1.00794,
                "C": 12.011,
                "O": 15.9994,
                "S": 32.06,}

        # number of heavy atoms in sidechains
        self.refNscat = {"ALA": 1,
                    "CYS": 2,
                    "ASP": 4,
                    "GLU": 5,
                    "PHE": 7,
                    "GLY": 0,
                    "HIS": 6,
                    "HSD": 6,
                    "HSE": 6,
                    "HSP": 6,
                    "ILE": 4,
                    "LYS": 5,
                    "LEU": 4,
                    "MET": 4,
                    "ASN": 4,
                    "PRO": 3,
                    "GLN": 5,
                    "ARG": 7,
                    "SER": 2,
                    "THR": 3,
                    "VAL": 3,
                    "TRP": 10,
                    "TYR": 8}

        # charges on side chains at pH 7
        self.refcharge = {"ALA": 0.0,
                    "CYS": 0.0,
                    "ASP": -1.0,
                    "GLU": -1.0,
                    "PHE": 0.0,
                    "GLY": 0.0,
                    "HIS": 0.0,
                    "HSD": 0.0,
                    "HSE": 0.0,
                    "HSP": 0.0,
                    "ILE": 0.0,
                    "LYS": 1.0,
                    "LEU": 0.0,
                    "MET": 0.0,
                    "ASN": 0.0,
                    "PRO": 0.0,
                    "GLN": 0.0,
                    "ARG": 1.0,
                    "SER": 0.0,
                    "THR": 0.0,
                    "VAL": 0.0,
                    "TRP": 0.0,
                    "TYR": 0.0}

        # Generic C_alpha side-chain center of mass distance
        self.lbs_nongo = {"ASP": 2.46916481058687,
                    "PRO": 1.87381801537346,
                    "LYS": 3.49738414814426,
                    "ILE": 2.25260184847053,
                    "TRP": 3.58251993741888,
                    "CYS": 2.06666004558289,
                    "HSD": 3.15209719417679,
                    "PHE": 3.38385541816659,
                    "HSP": 3.15209719417679,
                    "GLN": 3.08654121335,
                    "SER": 1.89840600762153,
                    "ASN": 2.46916481058687,
                    "VAL": 1.93953811063784,
                    "LEU": 2.56580983973678,
                    "TYR": 3.38981664391425,
                    "GLU": 3.07971386504681,
                    "ARG": 3.39687572938579,
                    "THR": 1.931721703272,
                    "ALA": 1.51146031725997,
                    "MET": 2.95389402456081,
                    "HIS": 3.15209719417679,
                    "HSE": 3.15209719417679}

        self.improper_nongo = {"ASP": 14.655341300544,
                        "PRO": 26.763068425539,
                        "LYS": 12.765248692601,
                        "ILE": 13.5446902008313,
                        "TRP": 11.4483488626106,
                        "CYS": 20.0484470024042,
                        "HSD": 14.9962640689562,
                        "PHE": 10.9217771918902,
                        "HSP": 14.9962640689562,
                        "GLN": 17.3050853491068,
                        "SER": 20.1390130256255,
                        "ASN": 14.655341300544,
                        "VAL": 13.3216022614598,
                        "LEU": 11.8137180266206,
                        "TYR": 12.2715081962165,
                        "GLU": 15.4130821146834,
                        "ARG": 15.5451613009777,
                        "THR": 16.2956083930276,
                        "ALA": 16.8418866013662,
                        "MET": 12.7046284165739,
                        "HIS": 14.9962640689562,
                        "HSE": 14.9962640689562}

        self.ang_sb_nongo = {"ASP": 120.380153696218,
                        "PRO": 125.127927161651,
                        "LYS": 119.523270610009,
                        "ILE": 118.791108398805,
                        "TRP": 130.018548241749,
                        "CYS": 110.512719347428,
                        "HSD": 116.815900172681,
                        "PHE": 122.937540996701,
                        "HSP": 116.815900172681,
                        "GLN": 116.182123224059,
                        "SER": 107.971234136647,
                        "ASN": 120.380153696218,
                        "VAL": 112.877421898116,
                        "LEU": 123.32179171436,
                        "TYR": 116.783314494739,
                        "GLU": 116.659068554985,
                        "ARG": 119.709740783191,
                        "THR": 111.719883260793,
                        "ALA": 108.623605160075,
                        "MET": 116.636559053295,
                        "HIS": 116.815900172681,
                        "HSE": 116.815900172681}

        self.ang_bs_nongo = {"ASP": 116.629356207687,
                        "PRO": 79.4932105625367,
                        "LYS": 119.779735484239,
                        "ILE": 116.923861483529,
                        "TRP": 100.858690902849,
                        "CYS": 114.816253227757,
                        "HSD": 115.848569293979,
                        "PHE": 112.804608190743,
                        "HSP": 115.848569293979,
                        "GLN": 119.106753006548,
                        "SER": 116.361829754186,
                        "ASN": 116.629356207687,
                        "VAL": 121.299281732077,
                        "LEU": 117.587011217416,
                        "TYR": 116.72484692836,
                        "GLU": 119.507585037498,
                        "ARG": 117.532816176021,
                        "THR": 117.044133956143,
                        "ALA": 120.747734648009,
                        "MET": 123.234171432545,
                        "HIS": 115.848569293979,
                        "HSE": 115.848569293979}

        # segment id relationships
        self.alphabet = list(map(chr, range(ord('A'), ord('Z')+1)))
        segid2num = {}
        for nseg, letter in enumerate(self.alphabet):
            segid2num[letter] = nseg
        self.segid2num = segid2num
            
        # mass of amino acids
        # UNSURE! about pro, arg, his and cys weights
        self.aaSCmass = {"ALA": 71.000000,
                    "CYS": 114.000000,
                    "ASP": 114.000000,
                    "GLU": 128.000000,
                    "PHE": 147.000000,
                    "GLY": 57.000000,
                    "HIS": 114.000000,
                    "HSD": 114.000000,
                    "HSE": 114.000000,
                    "HSP": 114.000000,
                    "ILE": 113.000000,
                    "LYS": 128.000000,
                    "LEU": 113.000000,
                    "MET": 131.000000,
                    "ASN": 114.000000,
                    "PRO": 114.000000,
                    "GLN": 128.000000,
                    "ARG": 114.000000,
                    "SER": 87.000000,
                    "THR": 101.000000,
                    "VAL": 99.000000,
                    "TRP": 186.000000,
                    "TYR": 163.000000}

        # vdw radius of sidechains
        self.rvdw = {"ALA": 2.51958406732374,
                "CYS": 2.73823091624513,
                "ASP": 2.79030096923572,
                "GLU": 2.96332591119925,
                "PHE": 3.18235414984794,
                "GLY": 2.25450393833984,
                "HIS": 3.04273820988499,
                "HSD": 3.04273820988499,
                "HSE": 3.04273820988499,
                "HSP": 3.04273820988499,
                "ILE": 3.09345983013354,
                "LYS": 3.18235414984794,
                "LEU": 3.09345983013354,
                "MET": 3.09345983013354,
                "ASN": 2.84049696898525,
                "PRO": 2.78004241717965,
                "GLN": 3.00796101305807,
                "ARG": 3.28138980397453,
                "SER": 2.59265585208464,
                "THR": 2.81059478021734,
                "VAL": 2.92662460060742,
                "TRP": 3.38869998431408,
                "TYR": 3.22881842919248}
        

        ## Check dependency installation ##
        # find stride resource
        self.stride_path = files('NCLEdetector.resources').joinpath('stride')
        self.logger.debug(f'stride_path: {self.stride_path}')

        #if os.popen('stride 2>&1').readlines()[0].strip().endswith('command not found'):
        if os.popen(f'{self.stride_path} 2>&1').readlines()[0].strip().endswith('command not found'):
            self.logger.error('Error: Essential software "stride" is not installed.\nPlease install stride before coarse-graining.')
            sys.exit()
        else:
            self.logger.info(f'STRIDE found')

        Header = f"""

        # Build CG Protein Model: Python version #
        #   Yang Jiang & Edward P. O'Brien Jr.   #
        #            Dept. of Chemistry          #
        #          Penn State University         #
                                                
        Configuration:
        pdbfile = {self.pdbfile}
        casm = {self.casm}
        nscal = {self.nscal}
        fnn = {self.fnn}
        potential_name = {self.potential_name}
        domain_file = {self.domain_file}
        sc_prefix = {self.sc_prefix}
        ca_prefix = {self.ca_prefix}    
        """
        self.logger.debug(Header)

        if self.domain_file != "None":
            self.nscal_0 = '1'
            self.nscal = 1
            self.logger.info('domain_file is defined, nscal will be ignored.\n')

        if self.casm != 0 and self.casm != 1:
            self.logger.error('ERROR: casm can only be either 0 (ca model) or 1 (ca-sidechain model).')

        if self.potential_name.upper().startswith('GENERIC'):
            words = self.potential_name.split('-')
            if len(words) == 1:
                self.logger.error("ERROR: Generic potential keyword must be invoked as 'generic-bt'")
                sys.exit()
            else:
                if words[-1].upper() != 'BT' and words[-1].upper() != 'MJ' and words[-1].upper() != 'KGS':
                    self.logger.error("ERROR: You can only invoke Generic potential keyword as 'generic-bt' or 'generic-mj' or 'generic-kgs'")
                    sys.exit()
                else:
                    self.potential_name = self.potential_name.upper()
                    self.logger.error("ERROR: The generic potential is not supported in this version.\nCoarse-graining terminated.")
                    sys.exit()
        else:
            self.potential_name = self.potential_name.upper()
        ### END: get info from control file ###

        ## BEGIND: Conditional Defaults ##
        if self.casm == 1:  
            self.ene_bsc = 0.37 # energy of a backbone-sidechain native contact (0.03 in old version)
            self.single_hbond_ene = 0.75 # energy of a hydrogen bond for everthing but helices (0.50 in old version)
            self.single_hbond_ene_helix = 0.75 # energy of a hydrogen bond in a helix (0.50 in old version)
            self.bondlength_go = 0 # non-Go bond length
            self.angle_dw = 0 # Go angle potential
            self.dihedral_go = 1 # Go dihedral potential
            self.improperdihed_go = 1 # Go improper dihedral potential
            
        else: 
            self.ene_bsc = 0.37;  
            self.single_hbond_ene = 0.75; # energy of a hydrogen bond for everthing but helices
            self.single_hbond_ene_helix = 0.75; # energy of a hydrogen bond in a helix
            self.bondlength_go = 0 # non-Go bond length
            self.angle_dw = 1 # double-well angle potential
            self.dihedral_go = 0 # non-Go dihedral potential
            self.improperdihed_go = 0 # non-Go improper dihedral potential

        # read domain nscal values if domain is defined
        dom_nscal = []
        ndomain = 0
        dom = []
        if self.domain_file != "None":
            if not os.path.exists(self.domain_file):
                self.logger.error("ERROR: File %s does not exist"%self.domain_file)
                sys.exit()
            f = open(self.domain_file)
            lines = f.readlines()
            f.close()
            for line in lines:
                line = line.strip()
                if line.startswith('scale factor'):
                    words = line.split('=')
                    dom_nscal.append(float(words[-1]))
                if line.startswith('domain'):
                    ndomain += 1
                    words = line.split('=')[-1].split('-')
                    words = [int(w) for w in words]
                    dom.append(words)
                    if words[0] > words[1]:
                        self.logger.error("ERROR: When defining the domains in the interface file, index %d is Greater than %d!"%(words[0], words[1]))
                        sys.exit()
            self.logger.info('%d domain(s) defined in the Domain file %s'%(ndomain, self.domain_file))
            if ndomain == 0:
                self.logger.error("ERROR: No domain definitions were read. Check the domain definition file!")
                sys.exit()
            self.logger.info("Domain information:")
            for i, d in enumerate(dom):
                self.logger.info("Domain %d: %d to %d"%(i+1, d[0], d[1]))
            self.logger.info("")
            if len(dom_nscal) != (1+ndomain)*ndomain/2:
                self.logger.error("ERROR: Incorrect number of interfaces assigned. (%d, should be %d)"%(len(dom_nscal)-ndomain, (ndomain-1)*ndomain/2))
                sys.exit()
        self.dom_nscal = dom_nscal
        self.ndomain = ndomain
        self.dom = dom
        # END read domain nscal values if domain is defined

        # initialize nonbonding potential
        root_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        if self.potential_name.startswith('MJ'):
            miya = files('NCLEdetector.resources.shared_files').joinpath('mj_contact_potential.dat')
        elif self.potential_name.startswith('KGS'):
            miya = files('NCLEdetector.resources.shared_files').joinpath('kgs_contact_potential.dat')
        elif self.potential_name.startswith('BT'):
            miya = files('NCLEdetector.resources.shared_files').joinpath('bt_contact_potential.dat')
        else:
            self.logger.error("ERROR: Unrecognized force-field %s"%self.potential_name)
            sys.exit()
        self.logger.debug(miya)

        eps = np.zeros((20,20))

        f = open(miya)
        lines = f.readlines()
        f.close()
        nrows = 0
        avg_mj = 0
        nmj = 0
        for line in lines:
            line = line.strip()
            if line.startswith('#'):
                continue
            if line.startswith('AA'):
                words = line.split()
                vec = []
                for w in words[1:]:
                    vec.append(self.res2n[w.upper()])
                if len(vec) != 20:
                    self.logger.error("ERROR: missing residues in file %s"%miya)
                    sys.exit()
            else:
                words = line.split()
                for tc, w in enumerate(words):
                    w = float(w)
                    if self.potential_name.startswith('MJ'):
                        eps[vec[nrows]][vec[tc]] = nscal * abs(w-1.2)
                        eps[vec[tc]][vec[nrows]] = nscal * abs(w-1.2)
                        avg_mj += nscal * abs(w-1.2)
                    elif self.potential_name.startswith('BT'):
                        eps[vec[nrows]][vec[tc]] = nscal * abs(w-0.6)
                        eps[vec[tc]][vec[nrows]] = nscal * abs(w-0.6)
                        avg_mj += nscal * abs(w-0.6)
                    elif self.potential_name.startswith('KGS'):
                        eps[vec[nrows]][vec[tc]] = nscal * abs(w-1.8)
                        eps[vec[tc]][vec[nrows]] = nscal * abs(w-1.8)
                        avg_mj += nscal * abs(w-1.8)
                    nmj += 1
                nrows += 1
                if nrows > 20:
                    self.logger.error("ERROR 2: missing residues in file %s: %d"%(miya, nrows))
                    sys.exit()
                if len(words) != nrows:
                    self.logger.error("ERROR 3: missing residues in file %s, %d != %d"%(miya, len(words), nrows))
                    sys.exit()
        self.eps = eps
        avg_mj = avg_mj/nmj
        self.avg_mj = avg_mj
        self.logger.info("The average %s interaction energy is %.4f\n"%(self.potential_name, self.avg_mj))
        # END initialize nonbonding potential

        # Read in the generic backbone dihedral potential of CL Brooks if NON-GO dihedrals
        # requested by user.
        if self.dihedral_go == 0:
            dihedb_nongo = [[[] for j in range(20)] for i in range(20)]
            kpot_f = files('NCLEdetector.resources.shared_files').joinpath('karanicolas_dihe_parm.dat')
            f = open(kpot_f)
            lines = f.readlines()
            f.close()
            nphi = 0
            r1_old = None
            r2_old = None
            for line in lines:
                line = line.strip()
                dat = line.split()
                r1 = dat[0].upper()
                r2 = dat[1].upper()
                if r1 != r1_old or r2 != r2_old:
                    nphi = 0
                dihedb_nongo[self.res2n[r1]][self.res2n[r2]].append([0.756*float(dat[2]), int(dat[3]), float(dat[4])])
                nphi += 1
                r1_old = r1
                r2_old = r2
                if nphi > 4:
                    self.logger.error("ERROR: nphi = %d upon reading in generic dihedral file"%nphi)
                    self.logger.debug(line)
                    sys.exit()
            self.dihedb_nongo = dihedb_nongo
        # END Read in the generic backbone dihedral potential
    #############################################################################################################

    ###################################################################################################
    # generate charmm .psf
    def create_psf(self, struct, ca_list, name):
        # creat backbone bonds
        for i in range(len(ca_list)-1):
            segid_list = [ca_list[i+j].residue.segid for j in range(2)]
            segid_list = list(set(segid_list))
            if len(segid_list) == 1:
                struct.bonds.append(pmd.topologyobjects.Bond(ca_list[i], ca_list[i+1]))
        # creat backbone-sidechain bonds if exist
        for ca_atom in ca_list:
            if len(ca_atom.residue.atoms) > 1:
                b_bead = ca_atom.residue.atoms[1]
                struct.bonds.append(pmd.topologyobjects.Bond(ca_atom, b_bead))
        # create Angles
        for atm in struct.atoms:
            bond_list = atm.bond_partners
            if len(bond_list) > 1:
                for i in range(len(bond_list)-1):
                    for j in range(i+1, len(bond_list)):
                        struct.angles.append(pmd.topologyobjects.Angle(bond_list[i], atm, bond_list[j]))
        # create Dihedrals
        for i in range(len(ca_list)-3):
            segid_list = [ca_list[i+j].residue.segid for j in range(4)]
            segid_list = list(set(segid_list))
            if len(segid_list) == 1:
                struct.dihedrals.append(pmd.topologyobjects.Dihedral(ca_list[i], ca_list[i+1], ca_list[i+2], ca_list[i+3]))
        # create Impropers
        for i in range(1, len(ca_list)-1):
            segid_list = [ca_list[i+j-1].residue.segid for j in range(3)]
            segid_list = list(set(segid_list))
            if len(segid_list) == 1 and len(ca_list[i].residue.atoms) > 1:
                b_bead = ca_list[i].residue.atoms[1]
                struct.impropers.append(pmd.topologyobjects.Improper(ca_list[i], ca_list[i-1], ca_list[i+1], b_bead))
        psffile = os.path.join(self.outdir, name+'.psf')
        self.logger.info(f'Writing {psffile}')
        struct.save(psffile, overwrite=True, vmd=False)
        return psffile
    # END generate charmm .psf
    ###################################################################################################

    ###################################################################################################
    # generate charmm .top
    def Create_rtf(self, struct, out_name):
        #global self.pdbfile, self.casm
        topfile = os.path.join(self.outdir, out_name+'.top')
        self.logger.info(f'Writing {topfile}')
        fo = open(topfile, 'w')
        if self.casm == 1:
            fo.write('* This CHARMM .top file describes a Ca-Cb Go model of %s\n*\n20 1\n'%self.pdbfile)
        else:
            fo.write('* This CHARMM .top file describes a Ca Go model of %s\n*\n20 1\n'%self.pdbfile)
        # MASS section
        fo.write('! backbone masses\n')
        for idx, atm in enumerate(struct.atoms):
            fo.write('MASS %-4s %-8s %.6f\n'%(str(idx+1), atm.type, atm.mass))
        fo.write('\n')
        fo.write('DECL +%s\n'%struct[0].name)
        fo.write('DECL -%s\n'%struct[0].name)
        fo.write('DECL #%s\n'%struct[0].name)
        # residue section
        for res in struct.residues:
            res_charge = 0
            for atm in res.atoms:
                res_charge += atm.charge
            fo.write('RESI %-6s %.1f\n'%(res.name, res_charge))
            fo.write('GROUP\n')
            for atm in res.atoms:
                fo.write('ATOM %s %-6s %.1f\n'%(atm.name, atm.type, atm.charge))
            if self.casm == 1 and len(res.atoms) != 1:
                fo.write("Bond %s %s  %s +%s\n"%(res.atoms[0].name, res.atoms[1].name, 
                                                res.atoms[0].name, res.atoms[0].name))
                fo.write("Angle -%s %s %s  %s %s +%s  -%s %s +%s\n"%(res.atoms[0].name, res.atoms[0].name, res.atoms[1].name, 
                                                                    res.atoms[1].name, res.atoms[0].name, res.atoms[0].name, 
                                                                    res.atoms[0].name, res.atoms[0].name, res.atoms[0].name))
                fo.write("DIHE -%s %s +%s #%s\n"%(res.atoms[0].name, res.atoms[0].name, res.atoms[0].name, res.atoms[0].name))
                fo.write("IMPH %s -%s +%s %s\n\n"%(res.atoms[0].name, res.atoms[0].name, res.atoms[0].name, res.atoms[1].name))
            else:
                fo.write('Bond %s +%s\n'%(res.atoms[0].name, res.atoms[0].name))
                fo.write('Angle -%s %s +%s\n'%(res.atoms[0].name, res.atoms[0].name, res.atoms[0].name))
                fo.write('DIHE -%s %s +%s #%s\n\n'%(res.atoms[0].name, res.atoms[0].name, res.atoms[0].name, res.atoms[0].name))
        # end section
        fo.write('END\n')
        fo.close()
        return topfile
    # END generate charmm .top
    ###################################################################################################

    ###################################################################################################
    def calc_distance(self, atom_1, atom_2):
        dist = ((atom_1.xx - atom_2.xx)**2 + (atom_1.xy - atom_2.xy)**2 + (atom_1.xz - atom_2.xz)**2)**0.5
        return dist
    ###################################################################################################

    ###################################################################################################
    def cg_energy_minimization(self, cor, prefix, prm_file):
        temp = 310
        np = '1'
        timestep = 0.015*picoseconds
        fbsolu = 0.05/picosecond
        temp = temp*kelvin

        psf_pmd = pmd.charmm.CharmmPsfFile(prefix+'.psf')
        psf = CharmmPsfFile(prefix+'.psf')
        top = psf.topology
        os.system('parse_cg_cacb_prm.py -p '+prm_file+' -t '+prefix+'.top')
        name = prm_file.split('.prm')[0]
        forcefield = ForceField(name+'.xml')
        
        template_map = {}
        for chain in top.chains():
            for res in chain.residues():
                template_map[res] = res.name
                    
        
        system = forcefield.createSystem(top, nonbondedCutoff=2.0*nanometer, 
                                        constraints=None, removeCMMotion=False, ignoreExternalBonds=True,
                                        residueTemplates=template_map)
        custom_nb_force = system.getForce(4)
        custom_nb_force.setUseSwitchingFunction(True)
        custom_nb_force.setSwitchingDistance(1.8*nanometer)
        custom_nb_force.setNonbondedMethod(custom_nb_force.CutoffNonPeriodic)
        
        # add position restraints
        force = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        force.addPerParticleParameter("k")
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")
        system.addForce(force)
        # END add position restraints
        
        # add position restraints for CA
        force = system.getForces()[-1]
        k = 100*kilocalorie/mole/angstrom**2
        for atm in top.atoms():
            if atm.name == 'A':
                force.addParticle(atm.index, (k, cor[atm.index][0], cor[atm.index][1], cor[atm.index][2]))
        
        integrator = LangevinIntegrator(temp, fbsolu, timestep)
        integrator.setConstraintTolerance(0.00001)
        # prepare simulation
        platform = Platform.getPlatformByName('CPU')
        properties = {'Threads': np}
        simulation = Simulation(top, system, integrator, platform, properties)
        simulation.context.setPositions(cor)
        simulation.context.setVelocitiesToTemperature(temp)
        energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalorie/mole)
        getEnergyDecomposition(stdout, simulation.context, system)
        self.logger.info('   Potential energy before minimization: %.4f kcal/mol'%energy)
        simulation.minimizeEnergy(tolerance=0.1*kilocalories_per_mole)
        energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalorie/mole)
        getEnergyDecomposition(stdout, simulation.context, system)
        self.logger.info('   Potential energy after minimization: %.4f kcal/mol'%energy)
        current_cor = simulation.context.getState(getPositions=True).getPositions()
        return current_cor
    ###################################################################################################

    ###################################################################################################
    # remove bond constraints of 0 mass atoms
    def rm_cons_0_mass(self, system):
        tag = 0
        while tag == 0 and system.getNumConstraints() != 0:
            for i in range(system.getNumConstraints()):
                con_i = system.getConstraintParameters(i)[0]
                con_j = system.getConstraintParameters(i)[1]
                mass_i = system.getParticleMass(con_i).value_in_unit(dalton)
                mass_j = system.getParticleMass(con_j).value_in_unit(dalton)
                if mass_i == 0 and mass_j == 0:
                    system.removeConstraint(i)
                    #print('Constraint %d is removed, range is %d'%(i, system.getNumConstraints()))
                    tag = 0
                    break
                elif mass_i == 0 or mass_j == 0:
                    system.removeConstraint(i)
                    #print('Constraint %d is removed, range is %d'%(i, system.getNumConstraints()))
                    system.getForce(0).addBond(con_i, con_j, 3.81*angstroms, 50*kilocalories/mole/angstroms**2)
                    tag = 0
                    break
                else:
                    tag = 1
    # END remove bond constraints of 0 mass atoms
    ###################################################################################################

    ###################################################################################################
    # energy decomposition 
    def forcegroupify(self, system):
        forcegroups = {}
        for i in range(system.getNumForces()):
            force = system.getForce(i)
            force.setForceGroup(i)
            f = str(type(force))
            s = f.split('\'')
            f = s[1]
            s = f.split('.')
            f = s[-1]
            forcegroups[i] = f
        return forcegroups
    ###################################################################################################

    ###################################################################################################
    def getEnergyDecomposition(self, handle, context, system):
        forcegroups = forcegroupify(system)
        energies = {}
        for i, f in forcegroups.items():
            try:
                states = context.getState(getEnergy=True, groups={i})
            except ValueError as e:
                self.logger.debug(str(e))
                energies[i] = Quantity(np.nan, kilocalories/mole)
            else:
                energies[i] = states.getPotentialEnergy()
        results = energies
        handle.write('    Potential Energy:\n')
        for idd in energies.keys():
            handle.write('      %s: %.4f kcal/mol\n'%(forcegroups[idd], energies[idd].value_in_unit(kilocalories/mole))) 
        return results
    ###################################################################################################

    ###################################################################################################
    def parse_cg_prm(self, prmfile:str, topfile:str):
        """
        Parse CHARMM parameter file and generate OpenMM XML file.
        """

        top_file_list = topfile.strip().split()

        command = 'pmd.charmm.CharmmParameterSet('
        for tf in top_file_list:
            command += '"'+tf + '", '
        command += 'prmfile)'
        self.logger.debug(command)
  
        param=eval(command)

        name = prmfile.split('.prm')
        file_name = name[0]

        openmm_param=pmd.openmm.parameters.OpenMMParameterSet.from_parameterset(param)
        openmm_param.write(file_name+'_tmp.xml', skip_duplicates=False)
        self.logger.info(f'Writing {file_name}_tmp.xml')

        dom = MD.parse(file_name+'_tmp.xml')
        root = dom.documentElement
        atom_type = root.getElementsByTagName('AtomTypes')
        residue = root.getElementsByTagName('Residues')
        os.remove(file_name+'_tmp.xml')

        root = ET.Element("ForceField")

        pf = open(prmfile, 'r')
        section = None
        node = None
        nbxmod = None
        ep = None
        kc = 138.935485
        ld = 1 # 10 Angstrom
        atom_type_list = [];
        num_atom = 0
        dihedral_array = []
        acoef_array = None
        bcoef_array = None
        ccoef_array = None
        nb_table = []
        nbfix_table = []
        try:
            for line in pf:
                line = line.strip()
                if not line:
                    # This is a blank line
                    continue
                if line.startswith('!'):
                    # This is a comment line
                    continue
                if line.startswith('ATOM'):
                    section = 'ATOM'
                    #node = ET.SubElement(root, "AtomTypes")
                    continue
                if line.startswith('BOND'):
                    section = 'BOND'
                    node = ET.SubElement(root, 'HarmonicBondForce')
                    continue
                if line.startswith('ANGLE'):
                    section = 'ANGLE'
                    node = ET.SubElement(root, 'CustomAngleForce', 
                        energy='-1/gamma*log(e); e=exp(-gamma*(k_alpha*(theta-theta_alpha)^2+epsilon_alpha))+exp(-gamma*k_betta*(theta-theta_betta)^2)')
                    ET.SubElement(node, 'PerAngleParameter', name='k_alpha')
                    ET.SubElement(node, 'PerAngleParameter', name='theta_alpha')
                    ET.SubElement(node, 'PerAngleParameter', name='k_betta')
                    ET.SubElement(node, 'PerAngleParameter', name='theta_betta')
                    ET.SubElement(node, 'PerAngleParameter', name='gamma')
                    ET.SubElement(node, 'PerAngleParameter', name='epsilon_alpha')
                    continue
                if line.startswith('DIHEDRAL'):
                    section = 'DIHEDRAL'
                    node = ET.SubElement(root, 'PeriodicTorsionForce')
                    continue
                if line.startswith('IMPHI'):
                    section = 'IMPROPER'
                    if len(dihedral_array) != 0:
                        proper_node = ET.SubElement(node, 'Proper', type1=dihedral_array[0], type2=dihedral_array[1], 
                                type3=dihedral_array[2], type4=dihedral_array[3])
                    n0 = 1
                    for index in range(4, len(dihedral_array), 3):
                        proper_node.set('k'+str(n0), dihedral_array[index])
                        proper_node.set('periodicity'+str(n0), dihedral_array[index+1])
                        proper_node.set('phase'+str(n0), dihedral_array[index+2])
                        n0 += 1
                    node = ET.SubElement(root, 'CustomTorsionForce', 
                        energy='k*min(dtheta, 2*pi-dtheta)^2; dtheta = abs(theta-theta0); pi = 3.1415926535')
                    ET.SubElement(node, 'PerTorsionParameter', name='k')
                    ET.SubElement(node, 'PerTorsionParameter', name='theta0')
                    continue
                if line.startswith('NONBONDED'):
                    section = 'NONBONDED'
                    words = line.split()
                    nbxmod = int(words[2])
                    continue
                if line.startswith('CUTNB'):
                    words = line.split()
                    ep = float(words[7])
                    node = ET.SubElement(root, 'CustomNonbondedForce', 
                        energy='ke*charge1*charge2/ep/r*exp(-r/ld)+kv*(a/r^12+b/r^10+c/r^6); '+
                        'ke=ke1*ke2; ep=ep1*ep2; ld=ld1*ld2; kv=kv1*kv2; '+
                        'a=acoef(index1, index2); b=bcoef(index1, index2); c=ccoef(index1, index2)',
                        bondCutoff=str(nbxmod-1))
                    ET.SubElement(node, 'PerParticleParameter', name='ke')
                    ET.SubElement(node, 'PerParticleParameter', name='kv')
                    ET.SubElement(node, 'PerParticleParameter', name='ep')
                    ET.SubElement(node, 'PerParticleParameter', name='ld')
                    ET.SubElement(node, 'PerParticleParameter', name='charge')
                    ET.SubElement(node, 'PerParticleParameter', name='index')
                    acoef_array = numpy.zeros((num_atom, num_atom))
                    bcoef_array = numpy.zeros((num_atom, num_atom))
                    ccoef_array = numpy.zeros((num_atom, num_atom))
                    nb_table = [[] for i in atom_type_list]
                    continue
                if line.startswith('NBFIX'):
                    section = 'NBFIX'
                    continue
                # It seems like files? sections? can be terminated with 'END'
                if line.startswith('END'): # should this be case-insensitive?
                    section = None
                    continue
                # If we have no section, skip
                if section is None: continue
                # Now handle each section specifically
                if section == 'ATOM':
                    words = line.split()
                    idx = int(words[1])
                    name = words[2]
                    mass = float(words[3])
                    #atom_node = ET.SubElement(node, 'Type', name=name, element='C', mass=str(mass))
                    #atom_node.set('class', name)
                    num_atom += 1
                    atom_type_list.append(name)
                if section == 'BOND':
                    words = line.split()
                    ET.SubElement(node, 'Bond', type1=words[0], type2=words[1], length=str(float(words[3])/10), k=str(float(words[2])*4.184*100*2))
                if section == 'ANGLE':
                    words = line.split()
                    ET.SubElement(node, 'Angle', type1=words[0], type2=words[1], type3= words[2],
                        k_alpha=str(float(words[3])*4.184), theta_alpha=str(float(words[4])/180*math.pi),
                        k_betta=str(float(words[5])*4.184), theta_betta=str(float(words[6])/180*math.pi),
                        gamma=str(float(words[7])/4.184), epsilon_alpha=str(float(words[8])*4.184))
                if section == 'DIHEDRAL':
                    words = line.split()
                    type1 = words[0]
                    type2 = words[1]
                    type3 = words[2]
                    type4 = words[3]
                    k = str(float(words[4])*4.184)
                    n = words[5]
                    phase = str(float(words[6])/180*math.pi)
                    if len(dihedral_array) == 0:
                        dihedral_array.append(type1)
                        dihedral_array.append(type2)
                        dihedral_array.append(type3)
                        dihedral_array.append(type4)
                        dihedral_array.append(k)
                        dihedral_array.append(n)
                        dihedral_array.append(phase)
                    elif (type1 == dihedral_array[0] and type2 == dihedral_array[1] and 
                        type3 == dihedral_array[2] and type4 == dihedral_array[3]):
                        dihedral_array.append(k)
                        dihedral_array.append(n)
                        dihedral_array.append(phase)
                    else:
                        proper_node = ET.SubElement(node, 'Proper', type1=dihedral_array[0], type2=dihedral_array[1], 
                            type3=dihedral_array[2], type4=dihedral_array[3])
                        n0 = 1
                        for index in range(4, len(dihedral_array), 3):
                            proper_node.set('k'+str(n0), dihedral_array[index])
                            proper_node.set('periodicity'+str(n0), dihedral_array[index+1])
                            proper_node.set('phase'+str(n0), dihedral_array[index+2])
                            n0 += 1
                        dihedral_array = []
                        dihedral_array.append(type1)
                        dihedral_array.append(type2)
                        dihedral_array.append(type3)
                        dihedral_array.append(type4)
                        dihedral_array.append(k)
                        dihedral_array.append(n)
                        dihedral_array.append(phase)
                if section == 'IMPROPER':
                    # No improper torsion energy term for Ca model
                    continue
                if section == 'NONBONDED':
                    words = line.split()
                    name = words[0]
                    epsilon = -float(words[2])*4.184
                    R_min_half = float(words[3])/10
                    index = atom_type_list.index(name)
                    nb_table[index] = [epsilon, R_min_half]
                if section == 'NBFIX':
                    words = line.split()
                    type1 = words[0]
                    type2 = words[1]
                    index1 = atom_type_list.index(type1)
                    index2 = atom_type_list.index(type2)
                    epsilon = -float(words[2])*4.184
                    R_min_half = float(words[3])/10
                    nbfix_table.append([index1, index2, epsilon, R_min_half])
        finally:
            pf.close()

        #Build acoef, bcoef, ccoef tables
        for index1 in range(num_atom):
            epsilon1 = nb_table[index1][0]
            R_min1 = nb_table[index1][1]
            for index2 in range(num_atom):
                epsilon2 = nb_table[index2][0]
                R_min2 = nb_table[index2][1]
                epsilon = numpy.sqrt(epsilon1 * epsilon2)
                R_min = R_min1 + R_min2
                a = 13 * epsilon * pow(R_min, 12)
                b = -18 * epsilon * pow(R_min, 10)
                c = 4 * epsilon * pow(R_min, 6)
                acoef_array[index1, index2] = a
                bcoef_array[index1, index2] = b
                ccoef_array[index1, index2] = c
        for nbfix_list in nbfix_table:
            index1 = nbfix_list[0]
            index2 = nbfix_list[1]
            epsilon = nbfix_list[2]
            R_min = nbfix_list[3]
            a = 13 * epsilon * pow(R_min, 12)
            b = -18 * epsilon * pow(R_min, 10)
            c = 4 * epsilon * pow(R_min, 6)
            acoef_array[index1, index2] = a
            bcoef_array[index1, index2] = b
            ccoef_array[index1, index2] = c
            acoef_array[index2, index1] = a
            bcoef_array[index2, index1] = b
            ccoef_array[index2, index1] = c

        #build tabulated function for acoef, bcoef, ccoef
        acoef_node = ET.SubElement(node, "Function", name='acoef', type='Discrete2D',
            xsize=str(num_atom), ysize=str(num_atom))
        text = ''
        for index1 in range(num_atom):
            for index2 in range(num_atom):
                text += str(acoef_array[index1, index2]) + " "
        acoef_node.text = text

        bcoef_node = ET.SubElement(node, "Function", name='bcoef', type='Discrete2D',
            xsize=str(num_atom), ysize=str(num_atom))
        text = ''
        for index1 in range(num_atom):
            for index2 in range(num_atom):
                text += str(bcoef_array[index1, index2]) + " "
        bcoef_node.text = text

        ccoef_node = ET.SubElement(node, "Function", name='ccoef', type='Discrete2D',
            xsize=str(num_atom), ysize=str(num_atom))
        text = ''
        for index1 in range(num_atom):
            for index2 in range(num_atom):
                text += str(ccoef_array[index1, index2]) + " "
        ccoef_node.text = text

        #add custom nonbond parameters
        ET.SubElement(node, 'UseAttributeFromResidue', name='charge')
        for index in range(num_atom):
            name = atom_type_list[index]
            ET.SubElement(node, 'Atom', type=name, index=str(index), ke=str(kc**0.5), ep=str(ep**0.5), ld=str(ld**0.5), kv='1')

        dom = MD.parseString(ET.tostring(root))
        root = dom.documentElement
        root = root.toprettyxml(indent=' ', newl='\n')
        dom = MD.parseString(root)
        root = dom.documentElement
        bond = root.getElementsByTagName('HarmonicBondForce')
        root.insertBefore(atom_type[0], bond[0])
        if len(residue) > 0:
            root.insertBefore(residue[0], bond[0])

        xf = open(file_name+'.xml', 'w')
        #dom.writexml(xf, indent='', addindent=' ', newl='\n')
        dom.writexml(xf, indent='')
    ###################################################################################################

    ###################################################################################################
    def parse_cg_cacb_prm(self, prmfile:str, topfile:str):


        top_file_list = topfile.strip().split()

        command = 'pmd.charmm.CharmmParameterSet('
        for tf in top_file_list:
            command += '"'+tf + '", '
        command += 'prmfile)'

        param=eval(command)

        name = prmfile.split('.prm')
        file_name = name[0]

        openmm_param=pmd.openmm.parameters.OpenMMParameterSet.from_parameterset(param)
        openmm_param.write(file_name+'_tmp.xml', skip_duplicates=False)
        dom = MD.parse(file_name+'_tmp.xml')
        root = dom.documentElement
        atom_type = root.getElementsByTagName('AtomTypes')
        residue = root.getElementsByTagName('Residues')
        os.remove(file_name+'_tmp.xml')

        root = ET.Element("ForceField")

        pf = open(prmfile, 'r')
        section = None
        node = None
        nbxmod = None
        ep = None
        kc = 138.935485
        ld = 1 # 10 Angstrom
        atom_type_list = [];
        num_atom = 0
        dihedral_array = []
        acoef_array = None
        bcoef_array = None
        ccoef_array = None
        nb_table = []
        nbfix_table = []
        try:
            for line in pf:
                line = line.strip()
                if not line:
                    # This is a blank line
                    continue
                if line.startswith('!'):
                    # This is a comment line
                    continue
                if line.startswith('ATOM'):
                    section = 'ATOM'
                    #node = ET.SubElement(root, "AtomTypes")
                    continue
                if line.startswith('BOND'):
                    section = 'BOND'
                    node = ET.SubElement(root, 'HarmonicBondForce')
                    continue
                if line.startswith('ANGLE'):
                    section = 'ANGLE'
                    node = ET.SubElement(root, 'HarmonicAngleForce')
                    continue
                if line.startswith('DIHEDRAL'):
                    section = 'DIHEDRAL'
                    node = ET.SubElement(root, 'PeriodicTorsionForce')
                    continue
                if line.startswith('IMPHI'):
                    section = 'IMPROPER'
                    if len(dihedral_array) != 0:
                        proper_node = ET.SubElement(node, 'Proper', type1=dihedral_array[0], type2=dihedral_array[1], 
                                type3=dihedral_array[2], type4=dihedral_array[3])
                    n0 = 1
                    for index in range(4, len(dihedral_array), 3):
                        proper_node.set('k'+str(n0), dihedral_array[index])
                        proper_node.set('periodicity'+str(n0), dihedral_array[index+1])
                        proper_node.set('phase'+str(n0), dihedral_array[index+2])
                        n0 += 1
                    node = ET.SubElement(root, 'CustomTorsionForce', 
                        energy='k*min(dtheta, 2*pi-dtheta)^2; dtheta = abs(theta-theta0); pi = 3.1415926535')
                    ET.SubElement(node, 'PerTorsionParameter', name='k')
                    ET.SubElement(node, 'PerTorsionParameter', name='theta0')
                    continue
                if line.startswith('NONBONDED'):
                    section = 'NONBONDED'
                    words = line.split()
                    nbxmod = int(words[2])
                    continue
                if line.startswith('CUTNB'):
                    words = line.split()
                    ep = float(words[7])
                    node = ET.SubElement(root, 'CustomNonbondedForce', 
                        energy='ke*charge1*charge2/ep/r*exp(-r/ld)+kv*(a/r^12+b/r^6); '+
                        'ke=ke1*ke2; ep=ep1*ep2; ld=ld1*ld2; kv=kv1*kv2; '+
                        'a=acoef(index1, index2); b=bcoef(index1, index2)',
                        bondCutoff=str(nbxmod-1))
                    ET.SubElement(node, 'PerParticleParameter', name='ke')
                    ET.SubElement(node, 'PerParticleParameter', name='kv')
                    ET.SubElement(node, 'PerParticleParameter', name='ep')
                    ET.SubElement(node, 'PerParticleParameter', name='ld')
                    ET.SubElement(node, 'PerParticleParameter', name='charge')
                    ET.SubElement(node, 'PerParticleParameter', name='index')
                    acoef_array = numpy.zeros((num_atom, num_atom))
                    bcoef_array = numpy.zeros((num_atom, num_atom))
                    nb_table = [[] for i in atom_type_list]
                    continue
                if line.startswith('NBFIX'):
                    section = 'NBFIX'
                    continue
                # It seems like files? sections? can be terminated with 'END'
                if line.startswith('END'): # should this be case-insensitive?
                    section = None
                    continue
                # If we have no section, skip
                if section is None: continue
                # Now handle each section specifically
                if section == 'ATOM':
                    words = line.split()
                    idx = int(words[1])
                    name = words[2]
                    mass = float(words[3])
                    #atom_node = ET.SubElement(node, 'Type', name=name, element='C', mass=str(mass))
                    #atom_node.set('class', name)
                    num_atom += 1
                    atom_type_list.append(name)
                if section == 'BOND':
                    words = line.split()
                    ET.SubElement(node, 'Bond', type1=words[0], type2=words[1], length=str(float(words[3])/10), k=str(float(words[2])*4.184*100*2))
                if section == 'ANGLE':
                    words = line.split()
                    ET.SubElement(node, 'Angle', type1=words[0], type2=words[1], type3= words[2],
                        k=str(float(words[3])*4.184*2), angle=str(float(words[4])/180*math.pi))
                if section == 'DIHEDRAL':
                    words = line.split()
                    type1 = words[0]
                    type2 = words[1]
                    type3 = words[2]
                    type4 = words[3]
                    k = str(float(words[4])*4.184)
                    n = words[5]
                    phase = str(float(words[6])/180*math.pi)
                    if len(dihedral_array) == 0:
                        dihedral_array.append(type1)
                        dihedral_array.append(type2)
                        dihedral_array.append(type3)
                        dihedral_array.append(type4)
                        dihedral_array.append(k)
                        dihedral_array.append(n)
                        dihedral_array.append(phase)
                    elif (type1 == dihedral_array[0] and type2 == dihedral_array[1] and 
                        type3 == dihedral_array[2] and type4 == dihedral_array[3]):
                        dihedral_array.append(k)
                        dihedral_array.append(n)
                        dihedral_array.append(phase)
                    else:
                        proper_node = ET.SubElement(node, 'Proper', type1=dihedral_array[0], type2=dihedral_array[1], 
                            type3=dihedral_array[2], type4=dihedral_array[3])
                        n0 = 1
                        for index in range(4, len(dihedral_array), 3):
                            proper_node.set('k'+str(n0), dihedral_array[index])
                            proper_node.set('periodicity'+str(n0), dihedral_array[index+1])
                            proper_node.set('phase'+str(n0), dihedral_array[index+2])
                            n0 += 1
                        dihedral_array = []
                        dihedral_array.append(type1)
                        dihedral_array.append(type2)
                        dihedral_array.append(type3)
                        dihedral_array.append(type4)
                        dihedral_array.append(k)
                        dihedral_array.append(n)
                        dihedral_array.append(phase)
                if section == 'IMPROPER':
                    words = line.split()
                    type1 = words[0]
                    type2 = words[1]
                    type3 = words[2]
                    type4 = words[3]
                    k = str(float(words[4])*4.184)
                    phase = str((float(words[6])-180)/180*math.pi)
                    improper_node = ET.SubElement(node, 'Improper', type1=type1, type2=type2, 
                        type3=type3, type4=type4, k=k,theta0=phase)
                    continue
                if section == 'NONBONDED':
                    words = line.split()
                    name = words[0]
                    epsilon = -float(words[2])*4.184
                    R_min_half = float(words[3])/10
                    index = atom_type_list.index(name)
                    nb_table[index] = [epsilon, R_min_half]
                if section == 'NBFIX':
                    words = line.split()
                    type1 = words[0]
                    type2 = words[1]
                    index1 = atom_type_list.index(type1)
                    index2 = atom_type_list.index(type2)
                    epsilon = -float(words[2])*4.184
                    R_min_half = float(words[3])/10
                    nbfix_table.append([index1, index2, epsilon, R_min_half])
        finally:
            pf.close()

        #Build acoef, bcoef, ccoef tables
        for index1 in range(num_atom):
            epsilon1 = nb_table[index1][0]
            R_min1 = nb_table[index1][1]
            for index2 in range(num_atom):
                epsilon2 = nb_table[index2][0]
                R_min2 = nb_table[index2][1]
                epsilon = numpy.sqrt(epsilon1 * epsilon2)
                R_min = R_min1 + R_min2
                a = epsilon * pow(R_min, 12)
                b = -2 * epsilon * pow(R_min, 6)
                acoef_array[index1, index2] = a
                bcoef_array[index1, index2] = b
        for nbfix_list in nbfix_table:
            index1 = nbfix_list[0]
            index2 = nbfix_list[1]
            epsilon = nbfix_list[2]
            R_min = nbfix_list[3]
            a = epsilon * pow(R_min, 12)
            b = -2 * epsilon * pow(R_min, 6)
            acoef_array[index1, index2] = a
            bcoef_array[index1, index2] = b
            acoef_array[index2, index1] = a
            bcoef_array[index2, index1] = b

        #build tabulated function for acoef, bcoef, ccoef
        acoef_node = ET.SubElement(node, "Function", name='acoef', type='Discrete2D',
            xsize=str(num_atom), ysize=str(num_atom))
        text = ''
        for index1 in range(num_atom):
            for index2 in range(num_atom):
                text += str(acoef_array[index1, index2]) + " "
        acoef_node.text = text

        bcoef_node = ET.SubElement(node, "Function", name='bcoef', type='Discrete2D',
            xsize=str(num_atom), ysize=str(num_atom))
        text = ''
        for index1 in range(num_atom):
            for index2 in range(num_atom):
                text += str(bcoef_array[index1, index2]) + " "
        bcoef_node.text = text

        #add custom nonbond parameters
        ET.SubElement(node, 'UseAttributeFromResidue', name='charge')
        for index in range(num_atom):
            name = atom_type_list[index]
            ET.SubElement(node, 'Atom', type=name, index=str(index), ke=str(kc**0.5), ep=str(ep**0.5), ld=str(ld**0.5), kv='1')

        dom = MD.parseString(ET.tostring(root))
        root = dom.documentElement
        root = root.toprettyxml(indent=' ', newl='\n')
        dom = MD.parseString(root)
        root = dom.documentElement
        bond = root.getElementsByTagName('HarmonicBondForce')
        root.insertBefore(atom_type[0], bond[0])
        if len(residue) > 0:
            root.insertBefore(residue[0], bond[0])

        xf = open(file_name+'.xml', 'w')
        #dom.writexml(xf, indent='', addindent=' ', newl='\n')
        dom.writexml(xf, indent='')
    ###################################################################################################

    ###################################################################################################
    def run(self,):

        resname_prefix = 'G'
        atomname_prefix = ''
                
        # Read PDB file
        cg_structure = pmd.Structure()
        self.logger.info("Reading in PDB file %s"%self.pdbfile)

        struct = pmd.load_file(self.pdbfile)
        sel_idx = np.zeros(len(struct.atoms))
        for idx, res in enumerate(struct.residues):
            res.number = idx+1
            if res.name in self.aa:
                for atm in res.atoms:
                    if atm.element != 1:
                        sel_idx[atm.idx] = 1
        heavy_protein = struct[sel_idx]

        for idx, res in enumerate(heavy_protein.residues):
            num_backbone = 0
            num_sidechain = 0
            for atm in res.atoms:
                if atm.name in ['C', 'N', 'O', 'CA']:
                    num_backbone += 1
                elif atm.name != 'OXT':
                    num_sidechain += 1
            if num_backbone != 4:
                self.logger.error("ERROR: In pdb the number of backbone atoms in residue %d is incorrect: %d != 4"%(idx+1, num_backbone))
                sys.exit()
            if num_sidechain != self.refNscat[res.name]:
                self.logger.error("ERROR: In pdb the number of sidechain atoms in residue %d is incorrect: %d != %d"%(idx+1, num_sidechain, self.refNscat[res.name]))
                sys.exit()

        idx_atm = 0
        ca_list = []
        chain_id_list = []
        for res in heavy_protein.residues:
            if not res.chain in chain_id_list:
                chain_id_list.append(res.chain)
        if len(chain_id_list) > len(self.alphabet):
            self.logger.error('ERROR: The number of chains in pdb file (%d) exceeds the maximum (%d)'%(len(chain_id_list), len(self.alphabet)))
            sys.exit()
        resid = 0
        chainid = chain_id_list[0]
        for idx, res in enumerate(heavy_protein.residues):
            if res.segid == '':
                segid = self.alphabet[chain_id_list.index(res.chain)]
            else:
                segid = res.segid
            
            if res.chain != chainid:
                chainid = res.chain
                resid = 1
            else:
                resid += 1
            
            SC_Mass = self.aaSCmass[res.name] - self.aaSCmass['GLY']
            CA_Mass = self.aaSCmass['GLY']
            SC_COM = np.zeros(3)
            CA_COM = np.zeros(3)
            sum_SC_Mass = 0
            
            for atm in res.atoms:
                if atm.name not in ['C', 'N', 'O', 'CA', 'OXT']:
                    sum_SC_Mass += atm.mass
                    SC_COM += atm.mass * np.array([atm.xx, atm.xy, atm.xz])
                elif atm.name == 'CA':
                    CA_COM[0] = atm.xx
                    CA_COM[1] = atm.xy
                    CA_COM[2] = atm.xz
            if sum_SC_Mass == 0:
                is_gly = True
            else:
                is_gly = False
                SC_COM /= sum_SC_Mass
            
            if self.casm == 0:
                cg_atm = pmd.topologyobjects.Atom(name=atomname_prefix+self.ca_prefix, 
                                                type=self.ca_prefix+str(idx+1), charge=self.refcharge[res.name], 
                                                mass=self.aaSCmass[res.name], number=idx_atm+1)
                cg_atm.xx = CA_COM[0]
                cg_atm.xy = CA_COM[1]
                cg_atm.xz = CA_COM[2]
                cg_structure.add_atom(cg_atm, resname_prefix+str(idx+1), resid, segid=segid, chain=res.chain)
                idx_atm += 1
                ca_list.append(cg_atm)
            else:
                ca_atm = pmd.topologyobjects.Atom(name=atomname_prefix+self.ca_prefix, 
                                                type=self.ca_prefix+str(idx+1), charge=0.0, 
                                                mass=CA_Mass, number=idx_atm+1)
                ca_atm.xx = CA_COM[0]
                ca_atm.xy = CA_COM[1]
                ca_atm.xz = CA_COM[2]
                cg_structure.add_atom(ca_atm, resname_prefix+str(idx+1), resid, segid=segid, chain=res.chain)
                idx_atm += 1
                ca_list.append(ca_atm)
                
                if not is_gly:
                    sc_atm = pmd.topologyobjects.Atom(name=atomname_prefix+self.sc_prefix, 
                                                    type=self.sc_prefix+str(idx+1), charge=self.refcharge[res.name], 
                                                    mass=SC_Mass, number=idx_atm+1)
                    sc_atm.xx = SC_COM[0]
                    sc_atm.xy = SC_COM[1]
                    sc_atm.xz = SC_COM[2]
                    cg_structure.add_atom(sc_atm, resname_prefix+str(idx+1), resid, segid=segid, chain=res.chain)
                    idx_atm += 1

        # Assign domain id to atom
        if self.ndomain != 0:
            self.logger.debug('Assign domain id to each atom')
            id_domain = []
            for atm in cg_structure.atoms:
                res_id = atm.residue.idx+1
                found = False
                for i, di in enumerate(self.dom):
                    if res_id >= di[0] and res_id <= di[1]:
                        id_domain.append(i)
                        found = True
                        break
                if not found:
                    self.logger.error('ERROR: %s is not located in any domain.'%atm)
                    sys.exit()
            self.logger.debug('')

        # Write psf, cor and top
        output_prefix = self.pdbfile.strip().split('/')[-1].split('.pdb')[0]
        if self.casm == 1:
            output_prefix += '_ca-cb'
        else:
            output_prefix += '_ca'
        self.logger.info('Create psf')
        psffile = self.create_psf(cg_structure, ca_list, output_prefix)

        self.logger.debug('Create cor')
        corfile = os.path.join(self.outdir, output_prefix+'.cor')
        self.logger.info(f'Writing {corfile}')
        cg_structure.save(corfile, overwrite=True, format='charmmcrd')

        self.logger.debug('Create top')
        topfile = self.Create_rtf(cg_structure, output_prefix)
            
        # Prepare FF parameters
        self.logger.info("Determining native contacts")
        dist_map = np.zeros((len(cg_structure.atoms), len(cg_structure.atoms)))
        for idx_1, atm_1 in enumerate(cg_structure.atoms):
            for idx_2, atm_2 in enumerate(cg_structure.atoms):
                dist_map[idx_1, idx_2] = self.calc_distance(atm_1, atm_2)
        self.logger.info("Finished calculating distance matrix")

        ## Compute native contacts between side-chains
        self.logger.info("Determining side-chains - side-chains contacts")
        native_ss_map = np.zeros((len(cg_structure.residues), len(cg_structure.residues)))
        for i in range(len(cg_structure.residues)-3):
            res_1 = heavy_protein.residues[i]
            for j in range(i+3, len(cg_structure.residues)): # separate by 2 residues
                res_2 = heavy_protein.residues[j]
                found = False
                for atm_1 in res_1.atoms:
                    for atm_2 in res_2.atoms:
                        if not atm_1.name in ['C', 'N', 'O', 'CA', 'OXT'] and not atm_2.name in ['C', 'N', 'O', 'CA', 'OXT']:
                            dij = self.calc_distance(atm_1, atm_2)
                            if dij <= self.heav_cut:
                                native_ss_map[i,j] = 1
                                native_ss_map[j,i] = 1
                                found = True
                                break
                    if found:
                        break
        ## Compute native contacts between backbone and side-chains
        self.logger.info("Determining backbone - side-chains contacts")
        native_bsc_map = np.zeros((len(cg_structure.residues), len(cg_structure.residues)))
        for i in range(len(cg_structure.residues)):
            res_1 = heavy_protein.residues[i]
            for j in range(len(cg_structure.residues)):
                res_2 = heavy_protein.residues[j]
                if i < j-2 or i > j+2: # separate by 2 residues
                    found = False
                    for atm_1 in res_1.atoms:
                        for atm_2 in res_2.atoms:
                            if atm_1.name in ['C', 'N', 'O', 'CA', 'OXT'] and atm_2.name not in ['C', 'N', 'O', 'CA', 'OXT']:
                                dij = self.calc_distance(atm_1, atm_2)
                                if dij <= self.heav_cut:
                                    native_bsc_map[i,j] = 1
                                    found = True
                                    break
                        if found:
                            break
        self.logger.info('# nat sc-sc contacts %d, # nat bb-sc contacts %d, and  # non-nat sc-sc %d' % (np.sum(native_ss_map)/2,
            np.sum(native_bsc_map), (len(cg_structure.residues)-3)*(len(cg_structure.residues)-2)/2 - np.sum(native_ss_map)/2))
            
        ## Determine hydrogen bonds that are present using STRIDE,
        ## and assign to Calpha-Calpha pairs. Also secondary structural elements
        ## within the native structure.
        self.logger.info("Determining the presence of hydrogen bonds using STRIDE")
        native_hb_map = np.zeros((len(cg_structure.residues), len(cg_structure.residues)))
        helical_list = np.zeros(len(cg_structure.residues))
        hb_ene_map = np.zeros((len(cg_structure.residues), len(cg_structure.residues)))
        #screen_out = os.popen(f'stride -h %s'%self.pdbfile).readlines()
        stride_cmd = f'{self.stride_path} -h {self.pdbfile}'
        #print(stride_cmd)
        screen_out = subprocess.run(stride_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        screen_out = screen_out.stdout.decode('utf-8', errors='replace').splitlines()
        #screen_out = os.popen(stride_cmd).read().splitlines()
        #print(screen_out)
    
        for line in screen_out:
            line = line.strip()
            resid = 0
            if line.startswith('ASD '):
                if 'Helix' in line.split()[6]:
                    helical_list[resid] = 1
                resid += 1
            if line.startswith('ACC ') or line.startswith('DNR '):
                # Get H-bonding info
                resid_1 = int(line[16:20])+1
                resid_2 = int(line[36:40])+1
                chainid_1 = line[8:10].strip()
                if chainid_1 == '-':
                    chainid_1 = ''
                chainid_2 = line[28:30].strip()
                if chainid_2 == '-':
                    chainid_2 = ''
                found = [0, 0]
                for idx, res in enumerate(cg_structure.residues):
                    if res.number == resid_1 and res.chain == chainid_1:
                        idx_1 = idx
                        found[0] = 1
                    elif res.number == resid_2 and res.chain == chainid_2:
                        idx_2 = idx
                        found[1] = 1
                    if sum(found) == 2:
                        break
                if sum(found) != 2:
                    self.logger.error("ERROR: Cannot find residue in parmed structure according to the Hbond info.\n  %s"%line)
                    sys.exit()
                if chainid_1 == chainid_2:
                    if idx_1 < idx_2:
                        if native_hb_map[idx_1, idx_2] == 1:
                            if helical_list[idx_1] == 1 and helical_list[idx_2] == 1:
                                hb_ene_map[idx_1, idx_2] = 2*self.single_hbond_ene_helix
                                hb_ene_map[idx_2, idx_1] = 2*self.single_hbond_ene_helix
                            else:
                                hb_ene_map[idx_1, idx_2] = 2*self.single_hbond_ene
                                hb_ene_map[idx_2, idx_1] = 2*self.single_hbond_ene
                        else:
                            native_hb_map[idx_1, idx_2] = 1
                            native_hb_map[idx_2, idx_1] = 1
                            if helical_list[idx_1] == 1 and helical_list[idx_2] == 1:
                                hb_ene_map[idx_1, idx_2] = self.single_hbond_ene_helix
                                hb_ene_map[idx_2, idx_1] = self.single_hbond_ene_helix
                            else:
                                hb_ene_map[idx_1, idx_2] = self.single_hbond_ene
                                hb_ene_map[idx_2, idx_1] = self.single_hbond_ene
                else:
                    native_hb_map[idx_1, idx_2] = 1
                    native_hb_map[idx_2, idx_1] = 1
                    hb_ene_map[idx_1, idx_2] = self.single_hbond_ene
                    hb_ene_map[idx_2, idx_1] = self.single_hbond_ene
        num_hb = 0
        for i in range(len(cg_structure.residues)-1):
            for j in range(i+1, len(cg_structure.residues)):
                if native_hb_map[i,j] == 1:
                    num_hb += 1
                    #print('%d %.4f, %d %d'%(num_hb, hb_ene_map[i,j], i+1, j+1))
        self.logger.info('# of unique Hbonds %d'%num_hb)
        native_contact_map = np.zeros((len(cg_structure.residues), len(cg_structure.residues)))
        for i in range(len(cg_structure.residues)):
            for j in range(len(cg_structure.residues)):
                if native_ss_map[i,j] == 1 or native_bsc_map[i,j] == 1 or native_hb_map[i,j] == 1:
                    native_contact_map[i,j] == 1

        ## Write prm file ##
        self.logger.debug('Create prm')
        prmfile = self.pdbfile.strip().split('/')[-1].split('.pdb')[0] + '_nscal' + str(self.nscal) + '_fnn' + str(self.fnn) + '_go_' + self.potential_name.lower() + '.prm'
        prmfile = os.path.join(self.outdir, prmfile)
        self.logger.info(f'Writing {prmfile}')

        f = open(prmfile, 'w')
        f.write('* This CHARMM .param file describes a Go model of %s\n'%(self.pdbfile.split('/')[-1]))
        f.write('*\n\n')
        # Atomic mass
        f.write('ATOM\n')
        for idx, atm in enumerate(cg_structure.atoms):
            f.write('MASS %-5s %-8s %-10.6f\n'%(str(idx+1), atm.type, atm.mass))
        f.write('\n')
        # Bond section (non-go bondlength for both models)
        f.write('BOND\n')
        kb = 50.0
        for idx, bond in enumerate(cg_structure.bonds):
            if self.bondlength_go == 0:
                if bond.atom2.name == (atomname_prefix+self.sc_prefix):
                    res_idx = bond.atom2.residue.idx
                    bond_length = self.lbs_nongo[heavy_protein.residues[res_idx].name]
                    f.write('%-8s%-10s%-12.6f%-9.6f\n'%(bond.atom1.type, bond.atom2.type, kb, bond_length))
                else:
                    f.write('%-8s%-10s%-12.6f%-9.6f\n'%(bond.atom1.type, bond.atom2.type, kb, 3.81))
            else:
                f.write('%-8s%-10s%-12.6f%-9.6f\n'%(bond.atom1.type, bond.atom2.type, kb, bond.measure()))
        f.write('\n')
        # Angle section
        f.write('ANGLE\n')
        ka = 30.0
        for idx, angle in enumerate(cg_structure.angles):
            if self.angle_dw == 0:
                f.write('%-8s%-8s%-10s%11.6f%11.6f\n'%(angle.atom1.type, angle.atom2.type, angle.atom3.type, 
                                                    ka, angle.measure()))
            else:
                if angle.atom1 == (atomname_prefix+self.sc_prefix):
                    res_idx = angle.atom1.residue.idx
                    angle_value = self.ang_sb_nongo[heavy_protein.residues[res_idx].name]
                    f.write('%-8s%-8s%-10s%11.6f%11.6f\n'%(angle.atom1.type, angle.atom2.type, angle.atom3.type, 
                                                        ka, angle_value))
                elif angle.atom3 == (atomname_prefix+self.sc_prefix):
                    res_idx = angle.atom3.residue.idx
                    angle_value = self.ang_bs_nongo[heavy_protein.residues[res_idx].name]
                    f.write('%-8s%-8s%-10s%11.6f%11.6f\n'%(angle.atom1.type, angle.atom2.type, angle.atom3.type, 
                                                        ka, angle_value))
                else:
                    f.write('%-8s%-8s%-10s  106.4 91.7 26.3 130.0 0.1 4.3\n'%(angle.atom1.type, angle.atom2.type, angle.atom3.type))
        f.write('\n')
        # Dihedral section
        f.write('DIHEDRAL\n')
        f.write('! backbone dihedrals\n')
        for idx, dihedral in enumerate(cg_structure.dihedrals):
            if self.dihedral_go == 1: # Use Go backbone dihedral angles
                delta = 1*dihedral.measure()-180
                if self.casm == 1:
                    if helical_list[dihedral.atom2.residue.idx] == 1 and helical_list[dihedral.atom3.residue.idx] == 1: # helical
                        kd = 0.30
                    else: # not helical
                        kd = 0.55
                    f.write('%-5s %-5s %-5s %-7s%-10.6f%-3d%-10.5f\n'%(dihedral.atom1.type, dihedral.atom2.type, dihedral.atom3.type, 
                                                                    dihedral.atom4.type, kd, 1, delta))
                    delta = 3*dihedral.measure()-180
                    if helical_list[dihedral.atom2.residue.idx] == 1 and helical_list[dihedral.atom3.residue.idx] == 1: # helical
                        kd = 0.15
                    else: # not helical
                        kd = 0.275
                    f.write('%-5s %-5s %-5s %-7s%-10.6f%-3d%-10.5f\n'%(dihedral.atom1.type, dihedral.atom2.type, dihedral.atom3.type, 
                                                                    dihedral.atom4.type, kd, 3, delta))
                else:
                    if helical_list[dihedral.atom2.residue.idx] == 1 and helical_list[dihedral.atom3.residue.idx] == 1: # helical
                        kd = 0.75
                    else: # not helical
                        kd = 0.75
                    f.write('%-5s %-5s %-5s %-7s%-10.6f%-3d%-10.5f\n'%(dihedral.atom1.type, dihedral.atom2.type, dihedral.atom3.type, 
                                                                    dihedral.atom4.type, kd, 1, delta))
                    delta = 3*dihedral.measure()-180
                    if helical_list[dihedral.atom2.residue.idx] == 1 and helical_list[dihedral.atom3.residue.idx] == 1: # helical
                        kd = 0.275
                    else: # not helical
                        kd = 0.275
                    f.write('%-5s %-5s %-5s %-7s%-10.6f%-3d%-10.5f\n'%(dihedral.atom1.type, dihedral.atom2.type, dihedral.atom3.type, 
                                                                    dihedral.atom4.type, kd, 3, delta))
            else: # Use Non-go dihedrals
                for i in range(4):
                    res_idx_1 = dihedral.atom2.residue.idx
                    res_idx_2 = dihedral.atom3.residue.idx
                    [kd, period, delta] = self.dihedb_nongo[self.res2n[heavy_protein.residues[res_idx_1].name]][self.res2n[heavy_protein.residues[res_idx_2].name]][i]
                    f.write('%-5s %-5s %-5s %-7s%-10.6f%-3d%-10.5f\n'%(dihedral.atom1.type, dihedral.atom2.type, dihedral.atom3.type, 
                                                                    dihedral.atom4.type, kd, period, delta))
        f.write('\n')
        # Improper dihedral section
        f.write('IMPHI\n')
        f.write('! sidechain improper dihedrals to maintain chirality\n')
        if self.casm == 1:
            for idx, improper in enumerate(cg_structure.impropers):
                if self.improperdihed_go == 1:
                    angle = improper.measure()
                else:
                    res_idx = improper.atom1.residue.idx
                    angle = self.improper_nongo[heavy_protein.residues[res_idx].name] # use transferable improper dihedral
                delta = angle + 180
                kd = 20*abs(self.avg_mj)
                f.write('%-5s %-5s %-5s %-7s%.6f %-3d%-10.5f\n'%(improper.atom1.type, improper.atom2.type, improper.atom3.type, 
                                                                improper.atom4.type, kd, 1, delta))
        f.write('\n')

        ## nonbonded section
        f.write('NONBONDED NBXMOD 3 ATOM CDIEL SWITCH VATOM VDISTANCE VSWITCH -\n')
        f.write('CUTNB 32 CTOFNB 20 CTONNB 18 EPS 78.5 WMIN 1.5 E14FAC 1.0\n')
        f.write('!atom           e_min   r_min/2\n')
        # if using the C-alpha only model do some preprocessing to determine the collision
        # diameter of non-native interactions according to the Karanacolis-Brooks
        # algorithm
        if self.casm != 1:
            sigmin = 1000000*np.ones(len(cg_structure.residues))
            if self.potential_name.startswith('GENERIC'):
                for idx, res in enumerate(cg_structure.residues):
                    sigmin[idx] = 2*self.rvdw[heavy_protein.residues[idx].name]
            else:
                # determine the collision diameter
                for i in range(len(cg_structure.residues)):
                    for j in range(len(cg_structure.residues)):
                        if native_contact_map[i,j] != 1 and (j < i-2 or j > i+2):
                            if dist_map[i,j] < sigmin[i]:
                                sigmin[i] = dist_map[i,j]
            for idx, atm in enumerate(cg_structure.atoms):
                eps2 = -0.000132
                rmin2 = sigmin[idx]*2**(1/6)/2
                temp = self.fnn*rmin2
                f.write("%-9s%-5.1f%-9.6f    %-10.6f\n"%(atm.type, 0.0, eps2, temp))
        else:
            eps2 = '-1e-12' #!!!! SYSTem dependent !!!!!!!!
            rmin2 = 20.0
            for idx, atm in enumerate(cg_structure.atoms):
                if atm.name == (atomname_prefix+self.ca_prefix):
                    f.write("%-9s%-5.1f%-s    %-10.6f\n"%(atm.type, 0.0, eps2, rmin2))
                else:
                    t1 = 1
                    t2 = (t1*(2*self.rvdw[heavy_protein.residues[atm.residue.idx].name]*2**(1/6))**12/(1e-12))**(1/12)
                    temp = self.fnn*t2/2
                    f.write("%-9s%-5.1f%-s    %-10.6f\n"%(atm.type, 0.0, eps2, temp))
        f.write('\n')
        ## NBFIX section
        f.write('NBFIX\n')
        ### native side-chain pairs and backbone Hbonding
        if self.casm == 1:
            f.write('! b-b due to Hbonding\n')
            totene_bb = 0
            for i in range(len(cg_structure.residues)-1):
                for j in range(i+1, len(cg_structure.residues)):
                    if native_hb_map[i,j] == 1:
                        atm_i = cg_structure.residues[i].atoms[0]
                        atm_j = cg_structure.residues[j].atoms[0]
                        comment = ''
                        if self.ndomain == 0: # No domain defined
                            ene = hb_ene_map[i,j]
                        else: # Domain defined
                            if id_domain[atm_i.idx] == id_domain[atm_j.idx]: # in the same domain
                                di = id_domain[atm_i.idx]
                                comment = '! in Domain %d'%(di+1)
                                ene = hb_ene_map[i,j]
                            else: # in the interface
                                di = id_domain[atm_i.idx]
                                dj = id_domain[atm_j.idx]
                                comment = '! in Interface %d | %d'%(di+1, dj+1)
                                ii =  int((2*self.ndomain - min(di, dj)) * (min(di, dj) + 1) / 2 + abs(di - dj) - 1)
                                #ene = self.dom_nscal[ii] # ??? Use nscal at interface
                                ene = hb_ene_map[i,j] # ??? Use the same energy
                        f.write('%-8s%-11s%-13.6f%-11.6f%s\n'%(atm_i.type, atm_j.type, -ene, dist_map[atm_i.idx, atm_j.idx], comment))
                        totene_bb += ene
            
            totene_sc = 0
            totene_bsc = 0
            if self.potential_name.startswith('GENERIC'): # C-alpha - side chain model Generic non-bond interactions
                f.write('!Generic interactions between unstructured portions of this protein\n')
                # Print out NBFIX energy values
                for i in range(len(cg_structure.residues)-3):
                    resname_1 = heavy_protein.residues[cg_structure.residues[i].idx].name
                    for j in range(i+3, len(cg_structure.residues)):
                        resname_2 = heavy_protein.residues[cg_structure.residues[j].idx].name
                        atm_i = cg_structure.residues[i].atoms[1] # ??? should be side-chain
                        atm_j = cg_structure.residues[j].atoms[1] # ??? should be side-chain
                        temp = self.rvdw[resname_1] + self.rvdw[resname_2]
                        ene=(0.3/10)*self.eps[self.res2n[resname_1]][self.res2n[resname_2]]
                        f.write('%-8s%-11s%-13.6f%-11.6f\n'%(atm_i.type, atm_j.type, -ene, temp))
            else: # Go non-bond interactions 
                f.write('! native side-chain interactions\n')
                for i in range(len(cg_structure.residues)-1):
                    resname_1 = heavy_protein.residues[cg_structure.residues[i].idx].name
                    for j in range(i+1, len(cg_structure.residues)):
                        resname_2 = heavy_protein.residues[cg_structure.residues[j].idx].name
                        if native_ss_map[i,j] == 1:
                            atm_i = cg_structure.residues[i].atoms[1]
                            atm_j = cg_structure.residues[j].atoms[1]
                            if self.eps[self.res2n[resname_1]][self.res2n[resname_2]] == 0:
                                self.logger.error('ERROR 1: Well depth equal to zero!!! %s - %s'%(resname_1, resname_2))
                                sys.exit()
                            comment = ''
                            if self.ndomain == 0: # No domain defined
                                ene = self.eps[self.res2n[resname_1]][self.res2n[resname_2]]
                            else: # If domain is defined
                                if id_domain[atm_i.idx] == id_domain[atm_j.idx]: # in the same domain
                                    di = id_domain[atm_i.idx]
                                    comment = '! in Domain %d'%(di+1)
                                    ene = self.eps[self.res2n[resname_1]][self.res2n[resname_2]] * self.dom_nscal[di] 
                                else: # in the interface
                                    di = id_domain[atm_i.idx]
                                    dj = id_domain[atm_j.idx]
                                    comment = '! in Interface %d | %d'%(di+1, dj+1)
                                    ii =  int((2*self.ndomain - min(di, dj)) * (min(di, dj) + 1) / 2 + abs(di - dj) - 1)
                                    ene = self.eps[self.res2n[resname_1]][self.res2n[resname_2]] * self.dom_nscal[ii] 
                            f.write('%-8s%-11s%-13.6f%-11.6f%s\n'%(atm_i.type, atm_j.type, -ene, dist_map[atm_i.idx, atm_j.idx], comment))
                            totene_sc += ene
                            
                f.write('! backbone-sidechain interactions\n')
                for i in range(len(cg_structure.residues)): 
                    resname_1 = heavy_protein.residues[cg_structure.residues[i].idx].name
                    for j in range(len(cg_structure.residues)):
                        resname_2 = heavy_protein.residues[cg_structure.residues[j].idx].name
                        if native_bsc_map[i,j] == 1:
                            atm_i = cg_structure.residues[i].atoms[0] # backbone
                            atm_j = cg_structure.residues[j].atoms[1] # sidechain
                            comment = ''
                            if self.ndomain == 0: # No domain defined
                                ene = self.ene_bsc
                            else: # If domain is defined
                                if id_domain[atm_i.idx] == id_domain[atm_j.idx]: # in the same domain
                                    di = id_domain[atm_i.idx]
                                    comment = '! in Domain %d'%(di+1)
                                    ene = self.ene_bsc
                                else: # in the interface
                                    di = id_domain[atm_i.idx]
                                    dj = id_domain[atm_j.idx]
                                    comment = '! in Interface %d | %d'%(di+1, dj+1)
                                    ii =  int((2*self.ndomain - min(di, dj)) * (min(di, dj) + 1) / 2 + abs(di - dj) - 1)
                                    ene = self.ene_bsc * self.dom_nscal[ii] # Rescaled energy
                            f.write('%-8s%-11s%-13.6f%-11.6f%s\n'%(atm_i.type, atm_j.type, -ene, dist_map[atm_i.idx, atm_j.idx], comment))
                            totene_bsc += ene
                            
            f.write('\n')
            f.write('! %.4f, %.4f, %.4f\n'%(totene_bb, totene_sc, totene_bsc))
        else:
            if not self.potential_name.startswith('GENERIC'): # C-alpha model
                f.write('! b-b due to Hbonding plus native side-chain interactions plus backbone-sidechain interactions\n')
                # Add up non-bonded energies
                for i in range(len(cg_structure.residues)-1): 
                    resname_1 = heavy_protein.residues[cg_structure.residues[i].idx].name
                    for j in range(i+1, len(cg_structure.residues)):
                        resname_2 = heavy_protein.residues[cg_structure.residues[j].idx].name
                        atm_i = cg_structure.residues[i].atoms[0]
                        atm_j = cg_structure.residues[j].atoms[0]
                        ene = 0
                        # hydrogen bonds
                        if native_hb_map[i,j] == 1:
                            if self.ndomain == 0: # No domain defined
                                ene += hb_ene_map[i,j]
                            else: # Domain defined
                                if id_domain[atm_i.idx] == id_domain[atm_j.idx]: # in the same domain
                                    di = id_domain[atm_i.idx]
                                    ene += hb_ene_map[i,j]
                                else: # in the interface
                                    di = id_domain[atm_i.idx]
                                    dj = id_domain[atm_j.idx]
                                    ii =  int((2*self.ndomain - min(di, dj)) * (min(di, dj) + 1) / 2 + abs(di - dj) - 1)
                                    ene += hb_ene_map[i,j] # Use the same energy
                        # sc-sc interactions
                        if native_ss_map[i,j] == 1:
                            if self.eps[self.res2n[resname_1]][self.res2n[resname_2]] == 0:
                                self.logger.error('ERROR 1: Well depth equal to zero!!! %s - %s'%(resname_1, resname_2))
                                sys.exit()
                            if self.ndomain == 0: # No domain defined
                                ene += self.eps[self.res2n[resname_1]][self.res2n[resname_2]]
                            else: # Domain defined
                                if id_domain[atm_i.idx] == id_domain[atm_j.idx]: # in the same domain
                                    di = id_domain[atm_i.idx]
                                    ene += self.eps[self.res2n[resname_1]][self.res2n[resname_2]] * self.dom_nscal[di]
                                else: # in the interface
                                    di = id_domain[atm_i.idx]
                                    dj = id_domain[atm_j.idx]
                                    ii =  int((2*self.ndomain - min(di, dj)) * (min(di, dj) + 1) / 2 + abs(di - dj) - 1)
                                    ene += self.eps[self.res2n[resname_1]][self.res2n[resname_2]] * self.dom_nscal[ii] 
                        # b-sc interactions
                        if native_bsc_map[i,j] == 1:
                            if self.ndomain == 0: # No domain defined
                                ene += self.ene_bsc
                            else: # Domain defined
                                if id_domain[atm_i.idx] == id_domain[atm_j.idx]: # in the same domain
                                    di = id_domain[atm_i.idx]
                                    ene += self.ene_bsc
                                else: # in the interface
                                    di = id_domain[atm_i.idx]
                                    dj = id_domain[atm_j.idx]
                                    ii =  int((2*self.ndomain - min(di, dj)) * (min(di, dj) + 1) / 2 + abs(di - dj) - 1)
                                    ene += self.ene_bsc # Use the same energy
                        if native_bsc_map[j,i] == 1:
                            if self.ndomain == 0: # No domain defined
                                ene += self.ene_bsc
                            else: # Domain defined
                                if id_domain[atm_i.idx] == id_domain[atm_j.idx]: # in the same domain
                                    di = id_domain[atm_i.idx]
                                    ene += self.ene_bsc
                                else: # in the interface
                                    di = id_domain[atm_i.idx]
                                    dj = id_domain[atm_j.idx]
                                    ii =  int((2*self.ndomain - min(di, dj)) * (min(di, dj) + 1) / 2 + abs(di - dj) - 1)
                                    ene += self.ene_bsc # Use the same energy
                        
                        # Write NBFIX
                        if ene != 0:
                            comment = ''
                            if self.ndomain != 0:
                                if id_domain[atm_i.idx] == id_domain[atm_j.idx]: # in the same domain
                                    di = id_domain[atm_i.idx]
                                    comment = '! in Domain %d'%(di+1)
                                else: # in the interface
                                    di = id_domain[atm_i.idx]
                                    dj = id_domain[atm_j.idx]
                                    comment = '! in Interface %d | %d'%(di+1, dj+1)
                            f.write('%-8s%-11s%-13.6f%-11.6f%s\n'%(atm_i.type, atm_j.type, -ene, dist_map[atm_i.idx, atm_j.idx], comment))
            else:
                f.write('!Generic interactions between unstructured portions of this protein\n')
                # Print out NBFIX energy values
                for i in range(len(cg_structure.residues)-3):
                    resname_1 = heavy_protein.residues[cg_structure.residues[i].idx].name
                    for j in range(i+3, len(cg_structure.residues)):
                        resname_2 = heavy_protein.residues[cg_structure.residues[j].idx].name
                        atm_i = cg_structure.residues[i].atoms[0]
                        atm_j = cg_structure.residues[j].atoms[0]
                        temp = self.rvdw[resname_1] + self.rvdw[resname_2]
                        ene=(0.3/10)*self.eps[self.res2n[resname_1]][self.res2n[resname_2]]
                        f.write('%-8s%-11s%-13.6f%-11.6f\n'%(atm_i.type, atm_j.type, -ene, temp))
        f.write('\nEND\n')
        f.close()

        self.logger.debug('All done.')
        return {'cor': corfile, 'prm': prmfile, 'psf': psffile, 'top': topfile,}
    ##########################################################################################


class BackMapping:
    """
    Take a C-alpha coarse grained structure and backmap it to the all-atom resolution
    """
    #############################################################################################################
    def __init__(self, nproc:int=1, outdir:str='./'):
        
        self.nproc = str(nproc)
        self.outdir = outdir
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
            self.logger.info(f'Made directory: {self.outdir}')

    #######################################################################################################

    #######################################################################################################
    def backmap(self, cg_pdb:str, aa_pdb:str, ID:str):
        """
        Backmap a C-alpha coarse-grained structure to an all-atom structure.
        """

        ##########################################################################
        self.logger.info(f"-> Cleaning PDB file {aa_pdb}")
        name = pathlib.Path(cg_pdb).stem + f'_{ID}'
        work_dir = os.path.join(self.outdir, 'rebuild_'+name)
        self.logger.info(name, work_dir)
        
        if not os.path.exists(work_dir):
            os.makedirs(work_dir)

        aa_clean_pdb, aa_clean_pdb_outfile = self.clean_pdb(aa_pdb, work_dir, name)
        os.chdir(work_dir)
        self.logger.debug('   Done')
        ##########################################################################

        ##########################################################################
        # buld ca-cb model
        self.logger.info(f"-> Building ca-cb model for {aa_clean_pdb_outfile}")
        (prefix, prm_file) = self.create_cg_model(aa_clean_pdb_outfile, ID)
        self.logger.debug('   Done')

        cacb_struct = pmd.load_file(prefix+'.psf')
        cacb_cor = pmd.load_file(prefix+'.cor')
        cacb_struct.coordinates = cacb_cor.coordinates
        #print(f'cacb_struct.coordinates: {cacb_struct.coordinates[:10]}')
        ##########################################################################
    
        ##########################################################################
        # add SC beads to cg pdb
        self.logger.info("-> Adding side chain beads")
        target_name = name
        cg_sc_struct = self.add_sc_beads(cg_pdb, cacb_struct)
        self.logger.debug('   Done')
        ##########################################################################

        ##########################################################################
        # run energy minimization for cacb model
        self.logger.info("-> Running energy minimization for ca-cb model")
        cg_sc_min_cor = self.cacb_energy_minimization(cg_sc_struct.positions, prefix, prm_file)
        aa_pdb_struct = pmd.load_file(aa_clean_pdb)
        for index in range(len(aa_pdb_struct.residues)):
            res_name = aa_pdb_struct.residues[index].name
            cg_sc_struct.residues[index].name = res_name
        for atm in cg_sc_struct.atoms:
            if atm.name == 'A':
                atm.name = ' CA'
            elif atm.name == 'B':
                atm.name = ' SC'

        # Remove units to perform arithmetic
        coordinates_values = cg_sc_min_cor.value_in_unit(nanometer)

        # Calculate the geometric center (center of mass could be calculated similarly if masses are available)
        geometric_center = np.mean(coordinates_values, axis=0)

        # Shift coordinates to center at the origin
        centered_coordinates = coordinates_values - geometric_center

        # Reapply the original unit
        centered_coordinates = centered_coordinates * nanometer

        #cg_sc_struct.positions = cg_sc_min_cor
        cg_sc_struct.positions = centered_coordinates
        self.logger.info(f'cg_sc_struct.positions: {cg_sc_struct.positions[:10]}')
        target_name_mini_pdb = target_name+'_mini.pdb'
        cg_sc_struct.save(target_name_mini_pdb, overwrite=True)
        self.logger.debug(f'SAVED: {target_name_mini_pdb}')
        self.logger.debug('   Done')

        #output_from_PD2 = target_name+'_mini.pdb'
        self.logger.info(f"-> Running Pulchra for {target_name_mini_pdb}")
        output_from_Pultra = self.Call_Pulchra(target_name_mini_pdb)
 

        ## remove the left over SC atoms that cause a template issue
        output_from_Pultra_cleaned = output_from_Pultra.replace('.pdb', '_cleaned.pdb')
        self.remove_sc_beads(output_from_Pultra, output_from_Pultra_cleaned)
 
        try:
            rec_pdb = self.OpenMM_vacuum_minimization(output_from_Pultra_cleaned, 500000)
            os.system('cp '+rec_pdb+' ../'+rec_pdb)
        except Exception as e:
            self.logger.info(traceback.print_exc(), e)
            self.logger.debug('Failed to run OpenMM minimization. Use Pulchra result instead.')
            rec_pdb = target_name+'_rebuilt.pdb'
            os.system('cp '+output_from_Pultra_cleaned+' ../'+rec_pdb)

        os.chdir('../')
        self.logger.info(f'Backmapping from {cg_pdb} -> {rec_pdb}')
    #######################################################################################################

    #######################################################################################################
    def clean_pdb(self, pdb, out_dir, name):
        AA_name_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
                        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
                        'HIE', 'HID', 'HIP']
        #name = pdb.split('/')[-1].split('.pdb')[0]
        struct = pmd.load_file(pdb)
        sel_idx = np.zeros(len(struct.atoms))
        for idx, res in enumerate(struct.residues):
            res.number = idx+1
            if res.name in AA_name_list:
                for atm in res.atoms:
                    sel_idx[atm.idx] = 1

        clean_pdb_outfile = os.path.join(out_dir, f'{name}_clean.pdb')
        self.logger.info(f'Writing {clean_pdb_outfile}')
        struct[sel_idx].save(clean_pdb_outfile, overwrite=True)
        return f'{name}_clean.pdb', clean_pdb_outfile
    #######################################################################################################

    #######################################################################################################    
    def create_psf(self, name):
        segid = 'A'
        parm = pmd.charmm.CharmmParameterSet(name+'.top')
        f = open(name+'.seq','r')
        seq = f.readlines()[0].strip().split()
        f.close()
        struct = pmd.Structure()
        for resname in seq:
            struct += parm.residues[resname].to_structure()
        ca_list = []
        for atm in struct.atoms:
            atm.mass = parm.atom_types[atm.type].mass
            if atm.name == 'A':
                ca_list.append(atm)
        # creat backbond bonds
        for i in range(len(ca_list)-1):
            struct.bonds.append(pmd.topologyobjects.Bond(ca_list[i], ca_list[i+1]))
        # create Angles
        for atm in struct.atoms:
            bond_list = atm.bond_partners
            if len(bond_list) > 1:
                for i in range(len(bond_list)-1):
                    for j in range(i+1, len(bond_list)):
                        struct.angles.append(pmd.topologyobjects.Angle(bond_list[i], atm, bond_list[j]))
        # create Dihedrals
        for i in range(len(ca_list)-3):
            struct.dihedrals.append(pmd.topologyobjects.Dihedral(ca_list[i], ca_list[i+1], ca_list[i+2], ca_list[i+3]))
        # create Impropers
        for i in range(1, len(ca_list)-1):
            if len(ca_list[i].residue.atoms) > 1:
                b_bead = ca_list[i].residue.atoms[1]
                struct.impropers.append(pmd.topologyobjects.Improper(ca_list[i], ca_list[i-1], ca_list[i+1], b_bead))
        for res in struct.residues:
            res.segid = segid
        struct.save(name+'.psf', overwrite=True)
    #######################################################################################################

    #######################################################################################################
    def create_cg_model(self, pdb, ID):
        self.logger.info(f'Creating CASM CG model from {pdb}')
        if not os.path.exists('./create_model'):
            os.makedirs('./create_model')
        os.chdir("./create_model")
        self.logger.debug(os.getcwd())

        CoarseGrainer = CoarseGrain(outdir='./',
                                        ID=ID,
                                        pdbfile=pdb,
                                        nscal=10.0,
                                        potential_name='mj',
                                        casm=1)
        self.logger.debug(CoarseGrainer)
    
        CGfiles = CoarseGrainer.run()
        self.logger.debug(os.getcwd())
 
        
        name = pathlib.Path(pdb).stem
        prefix = name+'_ca-cb'
        prm_name = name + '_nscal10.0_fnn1_go_mj.prm'
        self.logger.debug(f'name: {name}')
        self.logger.debug(f'prefix: {prefix}')
        self.logger.debug(f'prm_name: {prm_name}')
        
        if os.path.exists(prefix+'.psf'):
            os.system('cp *.psf ../')
            os.system('cp *.cor ../')
            os.system('cp *.top ../')
            os.system('cp *.prm ../')
            os.chdir('../')
        else:
            self.logger.error("Error: failed to create CG model from %s\n\n"%pdb)
            sys.exit()
        return (prefix, prm_name)
    #######################################################################################################

    #######################################################################################################
    def add_sc_beads(self, cg_pdb, cacb_struct):
        self.logger.info(f'Adding SC beads to {cg_pdb}')
        cor = pmd.load_file(cg_pdb)
        cor = cor.coordinates
        cor = cor[0]
        self.logger.info(cor, cor.shape)
        new_cacb_struct = cacb_struct.copy(pmd.Structure)
        idx = 0
        for res in new_cacb_struct.residues:
            res.atoms[0].xx = cor[idx,0]
            res.atoms[0].xy = cor[idx,1]
            res.atoms[0].xz = cor[idx,2]
            if len(res.atoms) > 1:
                cor1 = cacb_struct.coordinates[res.atoms[0].idx,:]
                cor2 = cacb_struct.coordinates[res.atoms[1].idx,:]
                bond_length = np.sum((cor1-cor2)**2)**0.5
                res.atoms[1].xx = cor[idx,0] + bond_length
                res.atoms[1].xy = cor[idx,1]
                res.atoms[1].xz = cor[idx,2]
            idx += 1
        return new_cacb_struct
    #######################################################################################################

    #######################################################################################################
    def cacb_energy_minimization(self, cor, prefix, prm_file):
        global nproc
        temp = 310
        timestep = 0.015*picoseconds
        fbsolu = 0.05/picosecond
        temp = temp*kelvin

        psf_pmd = pmd.charmm.CharmmPsfFile(prefix+'.psf')
        psf = CharmmPsfFile(prefix+'.psf')
        top = psf.topology

        # parse the cg cacb prm file
        topfile = f'{prefix}.top'
        self.logger.debug(os.getcwd())
        self.logger.debug(f'prm_file: {prm_file}')
        self.logger.debug(f'topfile: {topfile}')
        CoarseGrain.parse_cg_cacb_prm(self, prmfile=prm_file, topfile=topfile)
        xml_file = prm_file.split('.prm')[0]+'.xml'
        self.logger.debug(f'xml_file: {xml_file}')
        if not os.path.exists(xml_file):
            raise ValueError(f"Error: {xml_file} not found. Please check the file path.")
        
        #os.system('parse_cg_cacb_prm.py -p '+prm_file+' -t '+prefix+'.top')
        #name = prm_file.split('.prm')[0]
        forcefield = ForceField(xml_file)
        self.logger.debug(f'forcefield: {forcefield}')
        
        # re-name residues that are changed by openmm
        for resid, res in enumerate(top.residues()):
            if res.name != psf_pmd.residues[resid].name:
                res.name = psf_pmd.residues[resid].name
        
        template_map = {}
        for chain in top.chains():
            for res in chain.residues():
                template_map[res] = res.name
                    
        
        system = forcefield.createSystem(top, nonbondedCutoff=2.0*nanometer, constraints=None, 
                                        removeCMMotion=False, ignoreExternalBonds=True,
                                        residueTemplates=template_map)
        for force in system.getForces():
            if force.getName() == 'CustomNonbondedForce':
                custom_nb_force = force
                break
        # custom_nb_force = system.getForce(4)
        custom_nb_force.setUseSwitchingFunction(True)
        custom_nb_force.setSwitchingDistance(1.8*nanometer)
        custom_nb_force.setNonbondedMethod(custom_nb_force.CutoffNonPeriodic)
        
        # add position restraints
        force = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        force.addPerParticleParameter("k")
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")
        system.addForce(force)
        # END add position restraints
        
        # add position restraints for CA
        force = system.getForces()[-1]
        k = 100*kilocalorie/mole/angstrom**2
        for atm in top.atoms():
            if atm.name == 'A':
                force.addParticle(atm.index, (k, cor[atm.index][0], cor[atm.index][1], cor[atm.index][2]))
        
        integrator = LangevinIntegrator(temp, fbsolu, timestep)
        integrator.setConstraintTolerance(0.00001)
        # prepare simulation
        platform = Platform.getPlatformByName('CPU')
        properties = {'Threads': self.nproc}
        simulation = Simulation(top, system, integrator, platform, properties)
        simulation.context.setPositions(cor)
        simulation.context.setVelocitiesToTemperature(temp)
        energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalorie/mole)
        self.getEnergyDecomposition(stdout, simulation.context, system)
        self.logger.info('   Potential energy before minimization: %.4f kcal/mol'%energy)
        simulation.minimizeEnergy(tolerance=0.1*kilocalories_per_mole)
        energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalorie/mole)
        self.getEnergyDecomposition(stdout, simulation.context, system)
        self.logger.info('   Potential energy after minimization: %.4f kcal/mol'%energy)
        current_cor = simulation.context.getState(getPositions=True).getPositions()
        #print(f'current_cor:\n{current_cor[:10]}')
        return current_cor
    #######################################################################################################

    #######################################################################################################
    # energy decomposition 
    def forcegroupify(self, system):
        forcegroups = {}
        for i in range(system.getNumForces()):
            force = system.getForce(i)
            force.setForceGroup(i)
            f = str(type(force))
            s = f.split('\'')
            f = s[1]
            s = f.split('.')
            f = s[-1]
            forcegroups[i] = f
        return forcegroups
    #######################################################################################################

    #######################################################################################################
    def getEnergyDecomposition(self, handle, context, system):
        forcegroups = self.forcegroupify(system)
        energies = {}
        for i, f in forcegroups.items():
            try:
                states = context.getState(getEnergy=True, groups={i})
            except ValueError as e:
                self.logger.debug(str(e))
                energies[i] = Quantity(np.nan, kilocalories/mole)
            else:
                energies[i] = states.getPotentialEnergy()
        results = energies
        handle.write('    Potential Energy:\n')
        for idd in energies.keys():
            handle.write('      %s: %.4f kcal/mol\n'%(forcegroups[idd], energies[idd].value_in_unit(kilocalories/mole))) 
        return results
    #######################################################################################################

    #######################################################################################################
    def Call_Pulchra(self, rebult_pdb):
        self.logger.info("-> Calling pulchra to reconstruct all-atom PDB")
        self.pulchra = files('NCLEdetector.resources').joinpath('pulchra')
        self.logger.debug(f'pulchra: {self.pulchra}')
        pulchra_cmd = f'{self.pulchra} -v -g -q {rebult_pdb} > pulchra.log'
        self.logger.debug(f'CALL: {pulchra_cmd}')
        os.system(pulchra_cmd)

        pdb_code = rebult_pdb.split('.pdb')[0]
        old_name = pdb_code + ".rebuilt.pdb"
        new_name = pdb_code + "_pulchra.pdb"
        os.system("mv "+old_name+" "+new_name)
        self.logger.info("   Reconstructed all-atom PDB "+new_name)

        return new_name
    #######################################################################################################

    #######################################################################################################
    def OpenMM_vacuum_minimization(self, input_pdb, maxcyc):
        global nproc
        pdb_code = input_pdb.split('.pdb')[0]

        self.logger.info("-> Running all-atom energy minimization for %d steps in vacuum via OpenMM"%maxcyc)

        #platform = Platform.getPlatformByName('CUDA')
        #properties = {'CudaPrecision': 'mixed'}
        platform = Platform.getPlatformByName('CPU')
        properties = {'Threads': self.nproc}

        forcefield = ForceField('amber14-all.xml')
        self.logger.debug(f'input_pdb: {input_pdb}')
        pdb = pdbfile.PDBFile(input_pdb)
        self.logger.debug('FF made and PDB file loaded')

        # Check if the end residue has missing OXT atom and add if needed
        for chain in pdb.topology.chains():
            end_res = list(chain.residues())[-1]
            found = False
            for atom in end_res.atoms():
                if atom.name == 'OXT':
                    found = True
                elif atom.name == 'C':
                    C_atom = atom
                elif atom.name == 'CA':
                    CA_atom = atom
                elif atom.name == 'O':
                    O_atom = atom
            C_position = np.array(pdb.positions[C_atom.index].value_in_unit(nanometer))
            CA_position = np.array(pdb.positions[CA_atom.index].value_in_unit(nanometer))
            O_position = np.array(pdb.positions[O_atom.index].value_in_unit(nanometer))
            if not found:
                new_atom = pdb.topology.addAtom('OXT', element.oxygen, end_res)
                pdb.topology.addBond(C_atom, new_atom)
                new_position = np.dot(self.rotation_matrix(C_position-CA_position, np.pi), O_position-C_position) + C_position
                new_position = Quantity(value=Vec3(x=new_position[0], y=new_position[1], z=new_position[2]), unit=nanometer)
                pdb.positions.insert(O_atom.index+1, new_position)
        self.logger.debug('QC for OXT complete')

        model = modeller.Modeller(pdb.topology, pdb.positions)
        self.logger.debug(f'model: {model}')
        model.addHydrogens(forcefield=forcefield, pH=7.0)
        #model.addHydrogens(forcefield=forcefield, pH=7.0, variants=None, platform=platform)
        self.logger.debug('Hydrogens added')

        top = model.topology
        structure = pmd.openmm.load_topology(top)
        cor = model.positions
        #structure.positions = cor
        #structure.save('111.pdb', overwrite=True)
        
        system = forcefield.createSystem(top, nonbondedMethod=NoCutoff, constraints=None)
        self.logger.debug('System created')

        # add position restraints
        force = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        force.addPerParticleParameter("k")
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")
        system.addForce(force)
        self.logger.debug('Position restraints added')
        # END add position restraints
        
        # add position restraints for CA
        force = system.getForces()[-1]
        k = 500*kilocalorie/mole/angstrom**2
        for atm in top.atoms():
            if atm.name == 'CA':
                force.addParticle(atm.index, (k, cor[atm.index][0], cor[atm.index][1], cor[atm.index][2]))
        
        integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
        integrator.setConstraintTolerance(0.00001)
        self.logger.debug('Integrator set')

        simulation = Simulation(top, system, integrator, platform, properties)
        simulation.context.setPositions(cor)
        energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalorie/mole)
        self.getEnergyDecomposition(stdout, simulation.context, system)
        self.logger.info('   Potential energy before minimization: %.4f kcal/mol'%energy)

        simulation.minimizeEnergy(maxIterations=maxcyc)
        energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalorie/mole)
        self.getEnergyDecomposition(stdout, simulation.context, system)
        self.logger.info('   Potential energy after minimization: %.4f kcal/mol'%energy)
        current_cor = simulation.context.getState(getPositions=True).getPositions()
        
        structure.positions = current_cor
        outfile = pdb_code+'_OpenMM_min.pdb'
        structure['!@/H'].save(outfile, overwrite=True)
        self.logger.debug(f'SAVED: {outfile}')
        return outfile
    #######################################################################################################

    #######################################################################################################
    def remove_sc_beads(self, input_pdb, output_pdb):
        """
        Removes any atoms named 'SC' from a PDB file and writes the cleaned file.

        Parameters:
            input_pdb (str): Path to the input PDB file.
            output_pdb (str): Path to save the cleaned PDB file.
        """
        with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
            for line in infile:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    atom_name = line[12:16].strip()  # Extract the atom name
                    if atom_name == "SC":
                        continue  # Skip this line if the atom name is 'SC'
                outfile.write(line)  # Write all other lines

        self.logger.info(f"Cleaned PDB file saved to {output_pdb}")
    #######################################################################################################

    #######################################################################################################
    def rotation_matrix(self, axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        """
        axis = np.asarray(axis)
        axis = axis / math.sqrt(np.dot(axis, axis))
        a = math.cos(theta / 2.0)
        b, c, d = -axis * math.sin(theta / 2.0)
        aa, bb, cc, dd = a * a, b * b, c * c, d * d
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
        return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                        [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                        [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])
    #######################################################################################################
