import numpy as np
from LE4PD import properties


class dynamics(object):
	"""
	Attributes:
	-----------
	molecule: object
			This imports the molecule class which contains either the topology
			of protein or nucleic classes.
		temp: float (Default: 298 K)
			Desired temperature for thermal white noise.
		_fD2O: float
			Fraction of heavy water within solvent.
		_fH2O: float
			Fraction of standard water within solvent.
		_internal_viscosity: float
			Internal viscosity term.
		dt: float (Default: 1 ps)
			Time step in picoseconds used during simulation.
	Arguments:
	----------
		t0: float (Default: 0 ps)
			Start time in picoseconds. This defines time to start slicing the
			trajectory.
		tf: float (Default: -1)
			Stop time in picoseconds. This defines time to stop slice the
			trajectory.
	"""

	def __init__(self, molecule, temp=298, fD2O=0.0, internal_viscosity=2.71828,
				 mass_factor=1.0, NHfactor=1.0, n_iter=200, t0=0, tf=-1, dt=1):

		self.temp = temp
		self._fD2O = fD2O
		self._fH2O = 1.0 - self._fD2O
		self._internal_viscosity = internal_viscosity
		self._NHfactor = NHfactor
		self._n_iter = n_iter
		self.dt = dt

		# Import molecule class attributes
		self._MD = molecule._MD[t0:tf]

		# Remove any selected atoms
		if molecule._skip_atoms is not None:
			for atom in molecule._skip_atoms:
				if atom == "CA" or atom == "N" or atom == "H":
					warn("""Cannot skip atoms from the peptide backbone as these are coarse-graining sites. These atom types will be ignored.
					""")
					pass
				else:
					selection_criteria = "name != %s" % atom
					self._MD = self._MD.atom_slice(self._MD.top.select(selection_criteria))

		# Remove any select residues
		if molecule._skip_residues is not None:
			for residue in molecule._skip_residues:
				selection_criteria = "resname != %s" % residue
				self._MD = self._MD.atom_slice(self._MD.top.select(selection_criteria))

		self.top = self._MD.top
		self.xyz = self._MD.xyz
		self.n_conformers = self._MD.xyz.shape[0]
		self.atoms = molecule.atoms
		self.n_atoms = molecule.n_atoms
		self.residues = molecule.residues
		self.n_residues = molecule.n_residues

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
		properties.calculate_aspherocity(self)
		properties.calculate_P2(self)
		properties.calculate_mode_trajectory(self)
		properties.calculate_NMR_observables(self)

	def calculate_FES(self, bins=100):
		properties.calculate_FES(self, bins=bins)

	def save_modes_pdb(self, max_mode=10, max_conf=20):
		properties.save_modes_pdb(self, max_mode=max_mode, max_conf=max_conf)
