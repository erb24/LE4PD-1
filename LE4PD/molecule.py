from urllib.request import urlretrieve
import numpy as np
import mdtraj as md
import os.path
from warnings import warn


class protein(object):
	"""
	Some type of description of what the protein object is/does.
	"""

	def __init__(self, fetch=None, traj=None, top=None):
		if fetch != None and traj != None:
			warn("""\n
			Please select between building an ensemble from multiple coordinate files, or generating from MD simulations. LE4PD does not support building ensembles from both methods.
			""")

		if fetch is not None:
			if fetch.endswith(".pdb"):
				# Load structure and remove any solvent
				self._MD = md.load(fetch).remove_solvent()
			else:
				if not os.path.isfile(fetch + ".pdb"):
					print("Fetching structure from RCSB")
					url = 'http://www.rcsb.org/pdb/files/%s.pdb' % fetch
					urlretrieve(url, fetch + ".pdb")

				# Load structure and remove any solvent
				self._MD = md.load(fetch + ".pdb").remove_solvent()

			# Check that Hydrogens are in structure
			if len(self._MD.top.select("name == H")) == 0:
				warn("""Hydrogen atoms are not located within the topology file. Protein structure will be corected using Amber99sb.xml force-field""")
				from simtk.openmm.app import PDBFile, Modeller, ForceField
				pdb = PDBFile(fetch + ".pdb")
				modeller = Modeller(pdb.topology, pdb.positions)
				forcefield = ForceField('amber99sb.xml','tip3p.xml')
				modeller.addHydrogens(forcefield)
				PDBFile.writeFile(modeller.topology, modeller.positions, open(fetch + ".pdb", 'w'))
				self._MD = md.load(fetch + ".pdb").remove_solvent()

			self.n_conformers = self._MD.xyz.shape[0]
			if self.n_conformers == 1:
				warn(""" \n
				PDB file contains only one conformation. For accurate results, it is recommended to supply conformation data from NMR or CryoEM structures. Proceed with caution.
				""")
			self._method = "ensemble"

		if traj is not None:
			if top is None:
				self._MD = md.load(traj)
			else:
				self._MD = md.load(traj, top=top)
			self._method = "trajectory"

		self.top = self._MD.top
		self.xyz = self._MD.xyz
		self.atoms = np.squeeze([str(atom)[str(atom).index('-') + 1:]
								 for atom in self.top.atoms])
		self.n_atoms = len(self.atoms)
		self.residues = np.squeeze([str(residue)[:3]
									for residue in self.top.residues])
		self.n_residues = len(self.residues)

	def predict(self, temp=298, fD2O=0.0, int_visc=2.71828):
		if self._method == 'ensemble':
			from LE4PD.ensembles.dynamics import dynamics
		elif self._method == 'trajectory':
			from LE4PD.trajectories.dynamics import dynamics
		self.dynamics = dynamics(self)
		self.dynamics.predict

	def calculate_rmsd(self, reference=0, atom_indices=None, precentered=False):
		rmsd = md.rmsd(self._MD, self._MD, reference,atom_indices=atom_indices, precentered=precentered)
		self.rmsd = rmsd
