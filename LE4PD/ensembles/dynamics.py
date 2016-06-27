import numpy as np
from LE4PD.ensembles.util import properties


class dynamics(object):
	"""
	"""

	def __init__(self, protein, temp=298, fD2O=0.0, int_visc=2.71828, n_iter=200):
		self.temp = temp
		self._fD2O = fD2O
		self._fH2O = 1.0 - self._fD2O
		self._internal_viscosity = int_visc
		self._n_iter = n_iter

		self._MD = protein._MD
		self.top = protein.top
		self.xyz = protein.xyz
		self.n_conformers = protein.n_conformers
		self.atoms = protein.atoms
		self.n_atoms = protein.n_atoms
		self.residues = protein.residues
		self.n_residues = protein.n_residues
		if hasattr(protein, 'rmsd'):
			self.rmsd = protein.rmsd

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

	@property
	def save_modes_pdb(self):
		properties.save_modes_pdb(self)
