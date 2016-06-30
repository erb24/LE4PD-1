import numpy as np
from LE4PD.util import properties


class dynamics(object):
	"""
	Attributes:
	-----------
	molecule: object
			This imports the molecule class which contains either the topology
			of protein or nucleic classes.
		temp: float
			Desired temperature for thermal white noise.
		_fD2O: float
			Fraction of heavy water within solvent.
		_fH2O: float
			Fraction of standard water within solvent.
		_internal_viscosity: float
			Internal viscosity term.
	"""

	def __init__(self, molecule, temp=298, fD2O=0.0, int_visc=2.71828,
				 mass_factor=1.0, NHfactor=1.0, n_iter=200):

		self.temp = temp
		self._fD2O = fD2O
		self._fH2O = 1.0 - self._fD2O
		self._internal_viscosity = int_visc
		self._NHfactor = NHfactor
		self._n_iter = n_iter

		# Import molecule class attributes
		self._MD = molecule._MD
		self.top = molecule.top
		self.xyz = molecule.xyz
		self.n_conformers = molecule.n_conformers
		self.atoms = molecule.atoms
		self.n_atoms = molecule.n_atoms
		self.residues = molecule.residues
		self.n_residues = molecule.n_residues
		if hasattr(molecule, 'rmsd'):
			self.rmsd = molecule.rmsd

	@property
	def predict(self):
		properties.calculate_MSA(self)
		properties.calculate_SASA(self)
		properties.calculate_bond_vectors(self)
		properties.calculate_R_matrix(self)
		properties.calculate_gamma_contacts(self)
		properties.calculate_MSF(self)
		properties.calculate_U_matrix(self)
		properties.calculate_friction_coefficients(self)
		properties.calculate_sigma(self)
		properties.calculate_H_matrix(self)
		properties.calculate_M_matrix(self)
		properties.calculate_a_matrix(self)
		properties.calculate_L_matrix(self)
		properties.calculate_Q_matrix(self)
		properties.calculate_eigenvalues(self)
		properties.calculate_P2(self)

	def save_modes_pdb(self, max_modes=10):
		properties.save_modes_pdb(self, max_modes=max_modes)

	@property
	def calculate_NMR_observables(self):
		properties.calculate_NMR_observables(self)

	def calculate_rmsd(self, reference=0, atom_indices=None, precentered=False):
		properties.calculate_rmsd(self, reference=reference,atom_indices=atom_indices, precentered=precentered)
