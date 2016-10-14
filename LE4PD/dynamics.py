import numpy as np
from LE4PD.util import matrix, prepare

class ensemble(object):
	'''
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
	'''

	def __init__(self, molecule, temp=298, fD2O=0.0, internal_viscosity=2.71828,
				 mass_factor=1.0, NHfactor=1.0, n_iter=200, t0=0, tf=-1):

		self.temp = temp
		self._fD2O = fD2O
		self._fH2O = 1.0 - self._fD2O
		self._internal_viscosity = internal_viscosity
		self._NH_factor = NHfactor
		self._n_iter = n_iter

		# Import molecule class attributes
		self._mdtraj = molecule._mdtraj
		self._topfile = molecule._topfile
		self._skip_atoms = molecule._skip_atoms
		self._skip_residues = molecule._skip_residues

		# Import LE4PD topology
		self.topology = molecule.topology
		self.atoms = molecule.atoms
		self.n_atoms = molecule.n_atoms
		self.residues = molecule.residues
		self.n_residues = molecule.n_residues
		self.n_conformers = molecule.n_conformers

	def predict(self, timescale=4, probe_radius=0.14, n_sphere_points=250):
		from LE4PD.util import analyze_ensemble as analyze
		analyze.system(self, timescale=timescale, probe_radius=probe_radius,
						n_sphere_points=n_sphere_points)
		self.matrix = matrix.matrix(self)

	def calculate_FES(self, bins=100):
		from LE4PD.util import analyze_ensemble as analyze
		analyze.calculate_FES(self, bins=bins)

	def save_modes_pdb(self, max_mode=10, max_conf=20):
		from LE4PD.util import analyze_ensemble as analyze
		analyze.save_modes_pdb(self, max_mode=max_mode, max_conf=max_conf)

class simulation(object):
	'''
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
	'''

	def __init__(self, molecule, temp=298, fD2O=0.0, internal_viscosity=2.71828,
				 mass_factor=1.0, NHfactor=1.0, n_iter=200, t0=0, tf=-1):

		self.temp = temp
		self._fD2O = fD2O
		self._fH2O = 1.0 - self._fD2O
		self._internal_viscosity = internal_viscosity
		self._NH_factor = NHfactor
		self._n_iter = n_iter

		# Import molecule class attributes
		self._chunk_size = molecule._chunk_size
		self._mdtraj = molecule._mdtraj
		self._trajfile = molecule._trajfile
		self._topfile = molecule._topfile
		self._skip_atoms = molecule._skip_atoms
		self._skip_residues = molecule._skip_residues

		# Import LE4PD topology
		self.topology = molecule.topology
		self.atoms = molecule.atoms
		self.n_atoms = molecule.n_atoms
		self.residues = molecule.residues
		self.n_residues = molecule.n_residues

	def predict(self, timescale=4, probe_radius=0.14, n_sphere_points=250,stride=None):
		from LE4PD.util import analyze_simulation as analyze
		analyze.system(self, timescale=timescale, probe_radius=probe_radius,
						n_sphere_points=n_sphere_points, stride=stride)
		self.matrix = matrix.matrix(self)

	def calculate_FES(self, bins=100):
		from LE4PD.util import analyze_simulation as analyze
		analyze.calculate_FES(self, bins=bins)

	def save_modes_pdb(self, max_mode=10, max_conf=20):
		from LE4PD.util import analyze_simulation as analyze
		analyze.save_modes_pdb(self, max_mode=max_mode, max_conf=max_conf)
