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

import os, pathlib
from string import ascii_uppercase
from numpy import array, append
from Bio.PDB import PDBParser as PDBParserBiopy

class Vector:
    """A class representing Cartesian 3-dimensonal vectors."""
    
    def __init__(self, x,y,z):
        """x, y, z = Cartesian co-ordinates of vector."""
        self.x = x
        self.y = y
        self.z = z
        
    def copy(self):
        """
        Return:
            A copy of Vector instance
        """
        return Vector(self.x, self.y, self.z)
        
    def to_atom(self):
        """
        Create an Atom instance based on Vector instance.
        
        Return:
            Atom instance
        """
        atom = BioPyAtom([])
        atom.x = self.x
        atom.y = self.y
        atom.z = self.z
        return atom

class BioPy_Structure:
    
    
    """
    
    A class representing a bjectStructure o, as read from a PDB file using Bio.PDB in Biopython.
    
    
    """

    def __init__(self, atomList, filename='Unknown', header='', footer =''):
        """
        
        Initialise using a string of the relevant pdb file name or a numpy array of Atom objects.
        
        Arguments:
            *pdbFileOrList*
                String of pdb file name or array of Atom objects
            
        """
        self.header = header
        self.footer = footer
        self.filename = filename
        self.atomList = array(atomList)
        #Centre of mass calculations
        self.CoM = self.calculate_centre_of_mass()
        self.initCoM = self.CoM.copy()

    
    def __getitem__(self, index):
        return self.atomList[index]

    def __len__(self):
        return len(self.atomList)

    def __repr__(self):
        if not self.filename == 'Unknown':
            repr_str =  'Filename: ' + self.filename + '\n'
        else: 
            repr_str = ''
        repr_str += 'No Of Atoms: ' + str(len(self))  + '\n'
        repr_str += 'First Atom: ' + str(self.atomList[0]) + '\n'
        repr_str += 'Last Atom: ' + str(self.atomList[-1]) + '\n'
        return repr_str

    def copy(self):
        """
        
        Return:
            Copy of Structure instance.
        
        """
        newAtomList = []
        for atom in self.atomList:
            newAtomList.append(atom.copy())
        return BioPy_Structure(newAtomList)

    def calculate_centre_of_mass(self):
        """
        
        Return:    
            Center of mass of structure as a Vector instance.
        
        """
        x_momentTotal = 0.0
        y_momentTotal = 0.0
        z_momentTotal = 0.0
        massTotal = 0.0
        for atom in self.atomList:
            x = atom.get_x()
            y = atom.get_y()
            z = atom.get_z()
            m = atom.get_mass()
            x_momentTotal += x*m
            y_momentTotal += y*m
            z_momentTotal += z*m
            massTotal += m
            x_CoM = x_momentTotal/massTotal
            y_CoM = y_momentTotal/massTotal
            z_CoM = z_momentTotal/massTotal
        return Vector(x_CoM, y_CoM, z_CoM)

    def get_extreme_values(self):
        """
        
        Return:
            A 6-tuple containing the minimum and maximum of x, y and z co-ordinates of the structure.
            Given in order (min_x, max_x, min_y, max_y, min_z, max_z).
        
        """
        min_x = self.atomList[0].get_x()
        max_x = self.atomList[0].get_x()
        min_y = self.atomList[0].get_y()
        max_y = self.atomList[0].get_y()
        min_z = self.atomList[0].get_z()
        max_z = self.atomList[0].get_z()
        for atom in self.atomList[1:]:
            if atom.get_x() < min_x:
                min_x = atom.get_x()
            if atom.get_x() > max_x:
                max_x = atom.get_x()
            if atom.get_y() < min_y:
                min_y = atom.get_y()
            if atom.get_y() > max_y:
                max_y = atom.get_y()
            if atom.get_z() < min_z:
                min_z = atom.get_z()
            if atom.get_z() > max_z:
                max_z = atom.get_z()
        return (min_x, max_x, min_y, max_y, min_z, max_z)
    
class BioPyAtom:
    """
    
    A class representing an atom, as read from a PDB file using Biopython.
    
    """

    def __init__(self, atom):
        """Atom from BioPython"""
        if atom == []:
            return

        #http://deposit.rcsb.org/adit/docs/pdb_atom_format.html
        #print "bioatom",atom#'bioatom <Atom O>'
        if atom.get_parent().get_id()[0][0] == "W" or atom.get_parent().id[0][0]=="H":
            self.record_name = "HETATM"
        else:
            self.record_name = "ATOM" # was pdbString[:6].strip() as "ATOM"
#             res.id[0] == "W" or res.id[0][0]=="H": #skip water and hetero residues
        self.serial = atom.get_serial_number()
        self.atom_name = atom.get_name()
        self.alt_loc = atom.get_altloc() #Return alternative location specifier.
        self.fullid=atom.get_full_id()
        #('3ukr_test', 0, 'G', (' ', 113, ' '), ('CA', ' '))
        self.res = atom.get_parent().get_resname()
        self.chain = atom.get_full_id()[2]
        self.res_no = int(self.fullid[3][1])
        self.icode = ""
        if atom.is_disordered()==1:
            self.icode = "D"
             # 1 if the residue has disordered atoms
#            self.icode = pdbString[26].strip()#code for insertion residues
#             # Starting co-ordinates of atom.
        self.init_x = atom.get_coord()[0]
        self.init_y = atom.get_coord()[1]
        self.init_z = atom.get_coord()[2]
#             # Current co-ordinates of atom.
        self.x = float(atom.get_coord()[0])
        self.y = float(atom.get_coord()[1])
        self.z = float(atom.get_coord()[2])
#             
        self.occ = atom.get_occupancy()
        self.temp_fac = atom.get_bfactor()
        try:
            self.elem = atom.get_element()
        except:
                self.elem=""
        self.charge=""  
            #Mass of atom as given by atomicMasses global constant. Defaults to 1.
        self.mass = 1.0
        
#             # True if atom is the terminal of a chain. Automatically false until modified.
        self.isTerm = False
    
    def __repr__(self):
        return '('+ self.get_res() +' '+ str(self.res_no) + ' '+self.chain + ': ' + str(self.x) + ', ' + str(self.y) + ', ' + str(self.z) + ')'


    def copy(self):
        """
        
        Return:
            Copy of the Atom instance.
        """
        atom = BioPyAtom([])
        atom.record_name = self.record_name
        atom.serial = self.serial
        atom.atom_name = self.atom_name
        atom.alt_loc = self.alt_loc
        atom.res = self.res
        atom.chain = self.chain
        atom.res_no = self.res_no
        atom.icode = self.icode
        atom.init_x = self.init_x
        atom.init_y = self.init_y
        atom.init_z = self.init_z
        atom.x = self.x
        atom.y = self.y
        atom.z = self.z
        atom.occ =self.occ
        atom.temp_fac = self.temp_fac
        atom.elem = self.elem
        atom.charge = self.charge
        atom.mass = self.mass
        atom.isTerm = self.isTerm
        return atom

    def get_mass(self):
        """
        
        Return:
            Atom mass.
        """
        return self.mass

    def map_grid_position(self, densMap):
        """
                          
        Arguments:   
            *densMap*
                EM map object consisting the 3D grid of density values.
                
        Return:
              The co-ordinates and density value of the grid point in a density map closest to this atom.
              Return 0 if atom is outside of map.
        """
        x_origin = densMap.x_origin
        y_origin = densMap.y_origin
        z_origin = densMap.z_origin
        apix = densMap.apix
        x_size = densMap.x_size
        y_size = densMap.y_size
        z_size = densMap.z_size
        x_pos = int((self.getX()-x_origin)/apix)
        y_pos = int((self.getY()-y_origin)/apix)
        z_pos = int((self.getZ()-z_origin)/apix)
        if((x_size > x_pos >= 0) and (y_size > y_pos >= 0) and (z_size > z_pos >= 0)):
            return (x_pos, y_pos, z_pos, self.mass)
        else:
            return 0

    def get_x(self):
        """
        
        Return:
            x co-ordinate of atom.
        """
        return float(self.x)
    
    def get_y(self):
        """
        
        Return: 
            y co-ordinate of atom.
        """
        return float(self.y)
    
    def get_z(self):
        """
        
        Return:
            z co-ordinate of atom.
        """
        return float(self.z)
    
    
    
        
    def get_name(self):
        """
        atom name (ie. 'CA' or 'O')
        
        Return: 
            atom name.
        """
        return self.atom_name
    
    def get_res(self):
        """
        
        Return:
            three letter residue code corresponding to the atom (i.e 'ARG').
        """
        return self.res

    def get_res_no(self):
        """
        
        Return:
            residue number corresponding to the atom.
        """
        return self.res_no

    def get_id_no(self):
        """
        
        Return: 
            string of atom serial number.
        """
        return self.serial

    def write_to_PDB(self):
        """
        
        Writes a PDB ATOM record based in the atom attributes to a file.
        """
        line = ''
        line += self.record_name.ljust(6) 
        line += str(self.serial).rjust(5)+' '
        line += self.atom_name.center(4)
        line += self.alt_loc.ljust(1)
        line += self.res.ljust(3)+' '
        line += self.chain.ljust(1)
        line += str(self.res_no).rjust(4)
        line += str(self.icode).ljust(1)+'   '
        x = '%.3f' % self.x
        y = '%.3f' % self.y
        z = '%.3f' % self.z
        line += x.rjust(8)
        line += y.rjust(8)
        line += z.rjust(8)
        occ = '%.2f'% float(self.occ)
        temp_fac = '%.2f'% float(self.temp_fac)
        line += occ.rjust(6)
        line += temp_fac.rjust(6)+'          '
        line += self.elem.strip().rjust(2)
        line += self.charge.strip().ljust(2)
        return line + '\n'

def read_PDB_file(filename,hetatm=False,water=False):
        struct_file = open(filename, "r")
        # hydrogens are omitted.
        p=PDBParserBiopy(QUIET=True)#permissive default True
        structure=p.get_structure("id", struct_file)
        
        atomList = []
        hetatomList=[]
        wateratomList=[]
        footer = ''
        header = ''
        
        residues = structure[0].get_residues()
        for res in residues:
            hetfield=res.get_id()[0]
            if hetfield[0]=="H":
                for atom in res:
                    BioPyAtom(atom)
                    hetatomList.append(BioPyAtom(atom))
            elif hetfield[0]=="W":
                for atom in res:
                    BioPyAtom(atom)
                    wateratomList.append(BioPyAtom(atom))
            else:
                for atom in res:
                    if atom.id[0] != "H":
                        BioPyAtom(atom)
                        atomList.append(BioPyAtom(atom))
        if hetatm:
            atomList = append(atomList, hetatomList)
        if water:
            atomList = append(atomList, wateratomList)
        
        return BioPy_Structure(atomList, filename, header, footer)
        
def write_sasd_to_txt(sasds,pdb,result_dir):
    
    """
    
    Outputs sasds to .txt file
    
    Arguments:
    
       *sasds*
           dictionary of sasds
       *pdb*
           .pdb file sasds were calculated on
    """ 

    jwalk_pure_path = pathlib.PurePath(result_dir, 'Jwalk_results')
    jwalk_path = pathlib.Path(jwalk_pure_path)
    if not jwalk_path.exists():
        os.mkdir(jwalk_path)
    
    pdb = pathlib.Path(pdb)
    write_pure_path = pathlib.PurePath(jwalk_pure_path,'{}_crosslink_list.txt'.format(pdb.stem))
    write_path = pathlib.Path(write_pure_path)
    with open(write_path,'w') as outf:
        
        outf.write(' '.join('{0:<13}'.format(col) for col in ['Index','Model','Atom1','Atom2','SASD','Euclidean Distance']))
        outf.write('\n')
        index = 1
        
        for xl in sasds:
            (aa1,chain1,res1)=xl[0]
            (aa2,chain2,res2)=xl[1]
            atom1 = ('%s-%d-%s-CA' % (res1,aa1,chain1) )
            atom2 = ('%s-%d-%s-CA' % (res2,aa2,chain2) )
            sasd=xl[2]
            ed=xl[3]
            outf.write(' '.join('{0:<13}'.format(col) for col in [index,pdb.stem,atom1,atom2,sasd,ed]))
            outf.write('\n')
            index +=1
        
def write_sasd_to_pdb(dens_map,sasds,pdb,result_dir):
    
    """
    
    Outputs sasds to .pdb file
    
    Arguments:
       
       *dens_map*
           Solvent accessible surface on masked array
       *sasds*
           dictionary of sasds
       *pdb*
           .pdb file sasds were calculated on
    """ 
    jwalk_pure_path = pathlib.PurePath(result_dir, 'Jwalk_results')
    jwalk_path = pathlib.Path(jwalk_pure_path)
    if not jwalk_path.exists():
        os.mkdir(jwalk_path)
    
    apix = dens_map.apix
    origin = dens_map.origin
    path_coord = {}
    
    for xl in sasds:
        a = []
        for (x,y,z) in sasds[xl]:
            a.append([(x*apix)+origin[0], (y*apix)+origin[1], (z*apix)+origin[2]])
        
        path_coord[xl] = a
        
    pdb = pathlib.Path(pdb)
    write_pure_path = pathlib.PurePath(jwalk_path,'{}_crosslinks.pdb'.format(pdb.stem))
    write_path = pathlib.Path(write_pure_path)
    with open(write_path,'w') as pdb:
        # little trick to uniquely id all crosslinks with unique flase ATOM (X[A-Z]) / CHAIN ([A-Z]) name pairs 
        atom_cnt = 0
        chain_cnt = 0
        model_cnt = 1
        for xl in path_coord:
            (aa1,chain1,res1)=xl[0]
            (aa2,chain2,res2)=xl[1]

            atom_count_per_model = 1
            pdb.write('# MODEL {:d} {:s}{:d}{:s}-{:s}{:d}{:s}\n'.format(model_cnt,res1,aa1,chain1,res2,aa2,chain2))
            model_cnt += 1

            for (x,y,z) in path_coord[xl]:

                if atom_cnt > 25:
                    chain_cnt += 1
                    atom_cnt = 0
                
                atom_tmp = ascii_uppercase[atom_cnt]+'X'
                chain_tmp = ascii_uppercase[chain_cnt]

                p=Vector(x,y,z)
                a=p.to_atom()
                a.record_name = 'ATOM'
                a.serial = atom_count_per_model
                a.atom_name = atom_tmp
                a.alt_loc = ''
                a.res = chain_tmp+atom_tmp
                a.chain = chain_tmp
                a.res_no = atom_count_per_model
                a.icode = ''
                a.occ = 1
                a.temp_fac = 0
                a.elem = 'X'
                a.charge = ''
                #print a.__dict__
                #atom = BioPyAtom(a)
                pdb.write(a.write_to_PDB())

                atom_count_per_model += 1
            atom_cnt += 1

            # accounting for the extra count+=1 in the previous for loop
            atom_count_per_model -= 1
            # added for better visualization in pymol
            for i in range(1,atom_count_per_model):
                pdb.write('CONECT {:4d} {:4d}\n'.format(i,i+1))
            pdb.write('END\n')
            