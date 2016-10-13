import numpy as np
import mdtraj as md
import warnings
from os import path
from six.moves.urllib import request

def fetch(object, skip_atoms, skip_residues):
    '''
    Input:
    ------
    object: str
        The object can pass the path to a topology file. If the given object is
        of length 4 and does not contain a standard topology extension, fetch
        will search the RCSB Protein Data Bank and return a pdb file if found.
    '''
    if object.endswith('.pdb') or object.endswith('.gro'):
        file = object
        ID = object[:-4]
    elif len(object) == 4:
        ID = object
        if path.isfile(ID.upper()+'.pdb') or path.isfile(ID.lower()+'.pdb'):
            file = ID+'.pdb'
            warnings.warn('''File exists in current working directory. This will
            be used in place of fetching new file.''')
        else:
            try:
                # Fetch from RCSB
                print('Searching for structure from RCSB:PDB')
                url = 'http://www.rcsb.org/pdb/files/%s.pdb' % ID
                request.urlretrieve(url, ID.upper() + '.pdb')
                file = ID.upper()+'.pdb'
            except:
                raise ValueError('No PDB file retrieved.')
                pass

    molecule = solvent(md.load(file))
    molecule = check_hydrogens(molecule, ID)
    molecule = atoms(molecule, skip_atoms)
    molecule = residues(molecule, skip_residues)

    if not molecule.xyz.shape[0] == 1:
        molecule = molecule.superpose(molecule[0])
    return molecule.center_coordinates()

def trajectory(trajfile, topfile, skip_atoms, skip_residues):
    if topfile is None:
        molecule = solvent(md.load_frame(trajfile,0))
    else:
        molecule = solvent(md.load(topfile))
    molecule = atoms(molecule, skip_atoms)
    molecule = residues(molecule, skip_residues)
    molecule = molecule.superpose(molecule[0])
    return molecule.center_coordinates()

def topology(self):
    self.topology = self._mdtraj.top
    self.atoms = np.squeeze([str(atom)[str(atom).index('-') + 1:]
                             for atom in self.topology.atoms])
    self.n_atoms = len(self.atoms)
    self.residues = np.squeeze([str(residue)[:3]
                                for residue in self.topology.residues])
    self.n_residues = len(self.residues)
    if self._method == 'ensemble':
        self.n_conformers = self._mdtraj.n_frames

def solvent(molecule):
    return molecule.remove_solvent()

def check_hydrogens(molecule, ID):
    # Check that Hydrogens are in structure
    if len(molecule.top.select('name == H')) == 0:
        # If absent, then add Hydrogens using the Amber99sb force-field
        try:
            from simtk.openmm.app import PDBFile, Modeller, ForceField
            pdb = PDBFile(ID + '.pdb')
            modeller = Modeller(pdb.topology, pdb.positions)
            forcefield = ForceField('amber99sb.xml','tip3p.xml')
            modeller.addHydrogens(forcefield)
            PDBFile.writeFile(modeller.topology, modeller.positions, open(ID + '.pdb', 'w'))
            molecule = md.load(ID + '.pdb').remove_solvent()
        except:
            raise ValueError('''PDB topology missing Hydrogens. Either manually add
            or install OpenMM through SIMTK to automatically correct.''')
            pass
    return molecule

def atoms(molecule, skip_atoms):
    # Remove any selected atoms
    if skip_atoms is not None:
        for atom in skip_atoms:
            if atom == 'CA' or atom == 'N' or atom == 'H':
                warnings.warn('''Cannot skip atoms from the peptide backbone as
                these are coarse-graining sites. These atom types will be ignored.''')
                pass
            else:
                selection_criteria = 'name != %s' % atom
                molecule = molecule.atom_slice(molecule.top.select(selection_criteria))
    return molecule

def residues(molecule, skip_residues):
    # Remove any select residues
    if skip_residues is not None:
        for residue in skip_residues:
            selection_criteria = 'resname != %s' % residue
            molecule = molecule.atom_slice(molecule.top.select(selection_criteria))
    return molecule
