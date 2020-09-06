import warnings
import numpy as np
import subprocess
import os 

import mdtraj as md
from LE4PD.util import prepare
from LE4PD.codes import LE4PD, LE4PD_cmp

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
		'simulation': Trajectory file in .g96 format and an accompanying topology 
					  file in PDB format.
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
	def __init__(self, method='simulation', comp = False, isotropic = True, T = 298):
		self._method = method
		self.comp = comp 
		self.isotropic = isotropic
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

			if not self.comp:
				print("Assuming a single protein chain in the topology...")
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
			else:
				print("Multiple chains in the trajectory. Analyzing a protein complex...")
				if traj.split('.')[-1] != "g96":
					print("Only accepting .g96 file formats at this time. If you don't have a .g96 trajectory file for the alpha-carbons of \n")
					print("the protein, please run the 'process.sh' file in the main 'LE4PD' directory.")
				else:
					self._trajfile = traj
					self._topfile = top
					self.protname = traj.split('.')[0]
					N, NFRS, NATOMS = LE4PD_cmp.gen_protinfo(self.protname, self._trajfile, self._topfile)
					self.nmol, self.nres_list, self.natoms_list = LE4PD_cmp.get_chain_info(self.protname, self._topfile)
					self.nres = N
					self.nframes = NFRS
					self.natoms = NATOMS



	def prepare_traj(self):
		self.unformatted_traj = LE4PD.convert_traj(self._trajfile)
		self.traj = LE4PD.format_traj(self.unformatted_traj, self.nres, self.nframes)

	def calc_umatrix(self):
		if not self.comp:
			self.umatrix, self.Rinv, self.avbl, self.avblsq = LE4PD.Umatrix(self.traj, self.protname, self.nres, self.nframes, self.natoms)
		else:
			self.umatrix, self.Rinv, self.avbl, self.avblsq = LE4PD_cmp.Umatrix(self.traj, self.protname, self.nres, self.nframes, self.nmol, self.nres_list)

	def fric_calc(self, intv = 2.71828, viscosity = 1e-3, fd20 = 0.0, path_to_resarea = './'):
		if not self.comp:
			self.fratio, self.sigma, self.fric, self.avfr = LE4PD.fric_calc(self._topfile, self.protname, self.nres, self.nframes, 
																		self.natoms, self.avblsq, self.temp, intv = intv, viscosity = viscosity, fd20 = fd20,
																		path_to_resarea = path_to_resarea)
		else:
			self.fratio, self.sigma, self.fric, self.avfr = LE4PD_cmp.fric_calc(self._topfile, self.protname, self.nres, self.nframes, 
																		self.natoms, self.avblsq, self.temp, intv = intv, viscosity = viscosity, fd20 = fd20,
																		path_to_resarea = path_to_resarea)

	def lui_calc(self):
		if not self.comp:
			self.UILI, self.hmatrix, self.Q, self.QINV, self.lambda_eig, self.mu_eig = LE4PD.LUI_calc(self.protname, self.nres, 
																				self.nframes, self.umatrix, self.fratio, 
																				self.avblsq, self.sigma, self.fric, self.Rinv, 
																				self.temp)
		else:
			self.UILI, self.hmatrix, self.Q, self.QINV, self.lambda_eig, self.mu_eig = LE4PD_cmp.LUI_calc(self.protname, self.nres, 
																				self.nframes, self.nmol, self.nres_list, self.umatrix,
																				self.fratio, self.avblsq, self.sigma, self.fric, 
																				self.Rinv, self.temp)

	def mode_mad(self, nmodes = 10):
		if not self.comp:
			self.barriers, self.xi_traj, self.theta_phi_traj = LE4PD.mode_mad(self.traj, self.protname, self.nres, 
															self.nframes, self.Q, self.QINV, self.temp, nmodes = nmodes)
		else:
			self.barriers, self.xi_traj, self.theta_phi_traj = LE4PD_cmp.mode_mad(self.traj, self.protname, self.nres, 
															self.nframes, self.nmol, self.nres_list, self.Q, self.QINV, 
															self.temp, nmodes = nmodes)


	def tau_convert(self):
		if not self.comp:
			self.tau, self.tau_scaled = LE4PD.tau_convert(self.lambda_eig, self.sigma, self.barriers, self.temp)
		else:
			self.tau, self.tau_scaled = LE4PD_cmp.tau_convert(self.lambda_eig, self.sigma, self.barriers, self.temp)

	def LML(self):
		if not self.comp:
			self.LML = LE4PD.LML(self.Q, self.avbl, self.mu_eig)
		else:
			self.LML = LE4PD_cmp.LML(self.Q, self.avbl, self.mu_eig)