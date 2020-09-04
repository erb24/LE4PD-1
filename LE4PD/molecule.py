import warnings
import numpy as np
import subprocess
import os 

import mdtraj as md
from LE4PD.util import prepare, _m1int, _m1rot
from LE4PD.codes import LE4PD

'''Module to manage topology of varying molecule types. Currently only proteins
are implemented. Future versions will contain nucleic acids such as RNA and DNA.
'''

class protein(object):
	'''
	Input:
	------
	method: str
		The method determines what type of data will be supplied to predict the
		dynamics under the LE4PD formalism. This can be supplied in two forms:
		'simulation': Trajectory file supported by MDTraj and if necessary a
					  topology file.
		'ensemble':	  Any topology format that is supported by MDTraj,
					  additionally, structures can be fetched from the RCSB:PDB
					  [RCSB Protein Data Bank](http://www.rcsb.org/).
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
		recognized by MDtraj. Note: 'CA', 'N', nor 'H' can be skipped as these
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
	_mdtraj: mdtraj
		MDtraj generator which processes full topology of given input coordinate
		files, or full MD simulation files. See MDtraj documentation for full
		details.
	xyz: mdtraj
		Cartesian coordinates of supplied topology or trajectory data. See
		MDtraj documentation for full details. This is stored in as an H5PY
		datasets, to avoid memory constraints. (DEPRECATED)
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
	'''
	def __init__(self, method='simulation', T = 298):
		self._method = method
		self.temp = T

	def load(self, traj, top = None, skip_atoms = None, skip_residues = None, chunk_size = 1000):
		if self._method == 'ensemble':
			self._topfile = traj
			self._skip_atoms = skip_atoms
			self._skip_residues = skip_residues
			self._mdtraj = prepare.fetch(traj, skip_atoms, skip_residues)
			prepare.topology(self)

		elif self._method == 'simulation':
			#It isn't worth it right now to allow the trajectory conversion here. Probably better to give
			#the user a .sh file to do the conversion themselves. The modularized Python yahoos can get over it.
			#if traj.split('.')[-1] == "xtc":
				#Convert to g96 format by calling GROMACS using the subprocess module. 
				#Note that this action requires GROMACS be installed locally.
				#protname = traj.split('.')[0]
				#status = subprocess.call("echo '3' | `which gmx` -f " + str(traj) + " -s " + str(top) + " -o " + str(protname) + ".g96", shell = True)
			if traj.split('.')[-1] != "g96":
				print("Only accepting .g96 file formats at this time. If you don't have a .g96 trajectory file for the alpha-carbons of \n")
				print("the protein, please run the 'process.sh' file in the main 'LE4PD' directory.")
			else:
				self._trajfile = traj
				self._topfile = top
				self.protname = traj.split('.')[0]
				N, NFRS, NATOMS = LE4PD.gen_protinfo(self.protname, self._trajfile, self._topfile)
				self.nres = N
				self.nframes = NFRS
				self.natoms = NATOMS
				#self._skip_atoms = skip_atoms
				#self._skip_residues = skip_residues
				#self._chunk_size = chunk_size
				#self._mdtraj = prepare.trajectory(traj, top, skip_atoms, skip_residues)
				#self._n_frames = md.load(traj, top = top).n_frames
				#prepare.topology(self)

	'''def predict(self, temp=298, fD2O=0.0, internal_viscosity=2.71828,
				t0=0, tf=-1, dt=1, timescale=4,
				probe_radius=0.14, n_sphere_points=250,
				max_conf=100,stride=None):

		if self._method == 'ensemble':
			from LE4PD.dynamics import ensemble as dynamics
		elif self._method == 'simulation':
			from LE4PD.dynamics import simulation as dynamics

		# Prepare dynamics object
		self.dynamics = dynamics(self, temp=temp, fD2O=fD2O,
					internal_viscosity=internal_viscosity,
					t0=int(t0/dt), tf=int(tf/dt))

		# Predict dynamics with LE4PD
		if self._method == 'ensemble':
			if stride is not None:
				warnings('Stride feature is not supported for ensemble method. Stride set to None.')
			self.dynamics.predict(timescale=timescale, probe_radius=probe_radius,
								n_sphere_points=n_sphere_points)
		elif self._method == 'simulation':
			self.dynamics.predict(timescale=timescale, probe_radius=probe_radius,
							n_sphere_points=n_sphere_points,stride=stride)'''

	def prepare_traj(self):
		self.unformatted_traj = LE4PD.convert_traj(self._trajfile)
		self.traj = LE4PD.format_traj(self.unformatted_traj, self.nres, self.nframes)

	def calc_umatrix(self, protname, N, nfrs, natoms):
		self.umatrix, self.Rinv, self.avbl, self.avblsq = LE4PD.Umatrix(self.traj, self.protname, self.nres, self.nframes, self.natoms)

	def fric_calc(self, protname, N, nfrs, natoms, intv = 2.71828, viscosity = 1e-3, fd20 = 0.0):
		self.fratio, self.sigma, self.fric, self.avfr = LE4PD.fric_calc(self._topfile, self.protname, self.nres, self.nframes, 
																		self.natoms, self.avblsq, self.temp, intv = intv, viscosity = viscosity, fd20 = fd20)
	def lui_calc(self):
		self.UILI, self.hmatrix, self.Q, self.QINV, self.lambda_eig, self.mu_eig = LE4PD.LUI_calc(self.protname, self.nres, 
																				self.nframes, self.umatrix, self.fratio, 
																				self.avblsq, self.sigma, self.fric, self.Rinv, 
																				self.temp)
	def mode_mad(self, nmodes = 10):
		self.barriers, self.xi_traj, self.theta_phi_traj = LE4PD.mode_mad(self.traj, self.protname, self.nres, 
															self.nframes, self.Q, self.QINV, self.temp, nmodes = nmodes)

	def tau_convert(self):
		self.tau, self.tau_scaled = LE4PD.tau_convert(self.lambda_eig, self.sigma, self.barriers, self.temp)

	def LML(self):
		self.LML = LE4PD.LML(self.Q, self.avbl, self.mu_eig)

	'''def calculate_rmsd(self, frame = 0, atom_indices=None, precentered=False):
		if self._trajfile is None:
			self.rmsd = md.rmsd(self._mdtraj, reference=self._mdtraj,
								frame = frame,
								atom_indices=atom_indices,
								precentered=precentered)
		else:
			rmsd = []
			for chunk in md.iterload(self._trajfile, top = self.topology, chunk=self._chunk_size):
				chunk = prepare.atoms(chunk, self._skip_atoms)
				chunk = prepare.atoms(chunk, self._skip_residues)
				dummy = md.rmsd(chunk, reference=self._mdtraj, 
								frame = frame,
								atom_indices=atom_indices,
								precentered=precentered)
				for k in range(len(dummy)):
					rmsd.append(dummy[k])
		self.rmsd = np.array(rmsd)'''
