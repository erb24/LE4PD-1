import os.path
from warnings import warn
from urllib.request import urlretrieve
import numpy as np
import mdtraj as md
from LE4PD import properties

"""Module to manage topology of varying molecule types. Currently only proteins
are implemented. Future versions will contain nucleic acids such as RNA and DNA.
"""

class protein(object):
	"""
	Input:
	------
	fetch: str
		Fetch will read .pdb files and load the topology into an MDtraj object.
		If only an ID is given (if there is no .pdb extension), fetch will
		search [RCSB Protein Data Bank](http://www.rcsb.org/) for a protein
		structure.
	traj: dict
		Traj is a dictionary that stores the local directories of multiple
		MD trajectory files of any extension supported by MDtraj. Keys in
		dictionary should be indexed from 0.

		traj = {0: '/path/to/trajectory/0.xtc',
				1: '/path/to/trajectory/1.xtc',
							 ...
				N: '/path/to/trajectory/N.xtc'}

	 top: dict
		Top is a dictionary in the same formatting as traj, which stores the
		local directories of multiple topology files. Note that top is only
		necessary if the trajectory file format does not store topology
		information.

		top =  {0: '/path/to/topology/0.pdb',
				1: '/path/to/topology/1.pdb',
							 ...
				N: '/path/to/topology/N.pdb'}
	skip_atoms: list
		Skip_atoms is a list of strings that identifies any atom types that
		should be excluded from the analysis. This is ultimately to ignore
		non-standard atom types such as metals, halogens, etc. that are not
		recognized by MDtraj. Note: "CA", "N", nor "H" can be skipped as these
		are coarse-graining sites.
	skip_residues: list
		Skip_residues is a list of strings that identifies any residues types
		that should be excluded from the analysis. This is ultimately to ignore
		non-standard residue types such as ligands, amino acid derivatives, etc.
		that are not recognized by MDtraj. Note: If a residue is an amino acid
		without non-standard atom types, relabel the residue to their
		corresponding amino acid and use skip_atoms to remove any atoms not
		recognized by MDtraj.

	Attributes:
	-----------
	_MD: mdtraj
		MDtraj generator which processes full topology of given input coordinate
		files, or full MD simulation files. See MDtraj documentation for full
		details.
	xyz: mdtraj
		Cartesian coordinates of supplied topology or trajectory data. See
		MDtraj documentation for full details.
	top: mdtraj
		Topology information of supplied topology data. Trajectory selection can
		be used to slice portions of the desired trajectory. See MDtraj
		documentation for full details.
	atoms: str
		List of ordered atoms from topology.
	n_atoms: int
		Total number of atoms from topology.
	residues: str
		List of ordered residues from topology.
	n_residues: int
		Total number of residues from topology.
	n_conformers: int
		Total number of conformations from topology.

	"""

	def __init__(self, fetch=None, traj=None, top=None, skip_atoms=None, skip_residues=None):
		if fetch != None and traj != None:
			warn("""Please select between building an ensemble from multiple coordinate files, or generating from MD simulations. LE4PD does not support building ensembles from both methods.
			""")

		# Load structure from PDB file or RCSB Protein Data Bank
		if fetch is not None:
			if fetch.endswith(".pdb"):
				fetch = fetch[:len(fetch)-4]
			# Check if topology file already present in working directory
			else:
				try:
					# Fetch from RCSB
					print("Fetching structure from RCSB")
					url = 'http://www.rcsb.org/pdb/files/%s.pdb' % fetch
					urlretrieve(url, fetch + ".pdb")
				except:
					pass

			# Load structure and remove any solvent
			self._MD = md.load(fetch + ".pdb").remove_solvent()

			# Check that Hydrogens are in structure
			if len(self._MD.top.select("name == H")) == 0:
				# If absent, then add Hydrogens using the Amber99sb force-field
				warn("""Hydrogen atoms are not located within the topology file. Protein structure will be corected using Amber99sb.xml force-field""")
				from simtk.openmm.app import PDBFile, Modeller, ForceField
				pdb = PDBFile(fetch + ".pdb")
				modeller = Modeller(pdb.topology, pdb.positions)
				forcefield = ForceField('amber99sb.xml','tip3p.xml')
				modeller.addHydrogens(forcefield)
				PDBFile.writeFile(modeller.topology, modeller.positions, open(fetch + ".pdb", 'w'))
				self._MD = md.load(fetch + ".pdb").remove_solvent()

			if self._MD.xyz.shape[0] == 1:
				warn("""PDB file contains only one conformation. For accurate results, it is recommended to supply conformation data from NMR or CryoEM structures. Proceed with caution.
				""")

		# Load structure from trajectory files
		if traj is not None:
			if top is None:
				self._MD = md.load(traj).remove_solvent()
			else:
				self._MD = md.load(traj, top=top).remove_solvent()

		# Remove any selected atoms
		if skip_atoms is not None:
			for atom in skip_atoms:
				if atom == "CA" or atom == "N" or atom == "H":
					warn("""Cannot skip atoms from the peptide backbone as these are coarse-graining sites. These atom types will be ignored.
					""")
					pass
				else:
					selection_criteria = "name != %s" % atom
					self._MD = self._MD.atom_slice(self._MD.top.select(selection_criteria))

		# Remove any select residues
		if skip_residues is not None:
			for residue in skip_residues:
				selection_criteria = "resname != %s" % residue
				self._MD = self._MD.atom_slice(self._MD.top.select(selection_criteria))

		self.top = self._MD.top
		self.xyz = self._MD.xyz
		self.atoms = np.squeeze([str(atom)[str(atom).index('-') + 1:]
								 for atom in self.top.atoms])
		self.n_atoms = len(self.atoms)
		self.residues = np.squeeze([str(residue)[:3]
									for residue in self.top.residues])
		self.n_residues = len(self.residues)
		self.n_conformers = self._MD.xyz.shape[0]

	def predict(self, temp=298, fD2O=0.0, int_visc=2.71828):
		from LE4PD.dynamics import dynamics
		self.dynamics = dynamics(self)
		self.dynamics.predict

	def calculate_rmsd(self, reference=0, atom_indices=None, precentered=False):
		properties.calculate_rmsd(self, reference=reference,atom_indices=atom_indices, precentered=precentered)

	def calculate_COM(self, method="CA"):
		if method == "CA":
			properties.calculate_CA_COM(self)
		elif method == "All":
			print("Developer is lazy! This feature has not been added.")
		else:
			print("Developer is lazy! This feature has not been added.")
