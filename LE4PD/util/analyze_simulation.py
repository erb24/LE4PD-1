import warnings
import numpy as np
import mdtraj as md
from LE4PD.util import prepare
from LE4PD.util import _simulation_dynamics as _dynamics

def system(self, timescale=4, probe_radius=0.14, n_sphere_points=250, stride=None):
    calculate_bond_vectors(self, stride=stride)
    calculate_MSA(self)
    calculate_SASA(self, probe_radius=probe_radius, n_sphere_points=n_sphere_points, stride=stride)
    calculate_U_matrix(self, stride=stride)
    calculate_friction_coefficients(self)
    calculate_sigma(self)
    calculate_H_matrix(self)
    calculate_M_matrix(self)
    calculate_a_matrix(self)
    calculate_L_matrix(self)
    calculate_Q_matrix(self)
    calculate_eigenvalues(self)
    calculate_P2(self, timescale=timescale)
    calculate_mode_trajectory(self)
    calculate_NMR_observables(self)

def calculate_a_matrix(self):
    a = np.zeros((self.n_residues - 1, self.n_residues),
                 dtype=float, order='F')
    M = np.array(self._M, dtype=float, order='F')
    _dynamics.calculate_a_matrix(a, M, self.n_residues)
    self._a = a

def calculate_bond_vectors(self, stride=None):
    self.n_conformers = 0
    if self._trajfile is None:
        traj = self._topfile
        top = None
    else:
        traj = self._trajfile
        top = self._topfile

    for chunk in md.iterload(traj, top=top, chunk=self._chunk_size, stride=stride):
        chunk = prepare.solvent(chunk)
        chunk = prepare.atoms(chunk, self._skip_atoms)
        chunk = prepare.residues(chunk, self._skip_residues)
        self.n_conformers += chunk.n_frames
        # Set atomic slice range for alpha Carbons, Hydrogens, and Nitrogens from trajectory
        H_index = [atom.index for atom in self.topology.atoms if(
                    (str(atom.residue)[:3] == 'PRO' and
                    atom.name == 'CD') or
                    (atom.name == 'H'))]
        N_index = self.topology.select('name == N')
        C_index = self.topology.select('name == CA')

        # Convert slices into ndarray
        H = np.array(chunk.xyz[:, H_index, :], dtype=float, order='F')
        N = np.array(chunk.xyz[:, N_index, :], dtype=float, order='F')
        C = np.array(chunk.xyz[:, C_index, :], dtype=float, order='F')

        # Project bond dynamics on NH bond vectors
        mat = np.zeros((self.n_residues - 1, self.n_residues - 1), dtype=float, order='F')
        sbl = np.zeros((self.n_residues-1), dtype=float, order='F')
        _dynamics.calculate_bonds(mat, sbl, H, N, C, chunk.n_frames, self.n_residues)

    self._NH_matrix = mat/self.n_conformers
    self._average_bond_lengths = sbl/self.n_conformers


def calculate_eigenvalues(self):
    # Calculate lambda eigenvalues and Q-matrix
    w_lambda, Q = np.linalg.eig(self._LU)

    # Sort by eigenvalues
    idx = np.argsort(abs(w_lambda))[::-1]
    w_lambda = np.array(1 / (w_lambda[idx].real), dtype=float, order='F')
    Q = np.array(Q[:, idx], dtype=float, order='F')
    try:
        QI = np.array(np.linalg.inv(Q), dtype=float, order='F')
    except:
        QI = np.array(np.linalg.pinv(Q), dtype=float, order='F')

    if w_lambda.min() < 0:
        io = 0.0
    else:
        io = 1.0

    # Calculate mu eigenvalues
    w_mu = np.zeros(self.n_residues - 1, dtype=float, order='F')
    U = np.array(self._UI, dtype=float, order='F')
    _dynamics.calculate_nu_eigenvalues(w_mu, QI, U, self.n_residues)

    # Calculate nu eigenvalues
    w_nu = np.array(np.diag(np.dot(QI, np.dot(self._L, Q))),
                    dtype=float, order='F')

    self.lambda_eigenvalues = w_lambda
    self.mu_eigenvalues = w_mu
    self.nu_eigenvalues = w_nu

'''MODIFY THIS SO THAT FES ARE CALCULATED FROM H5PY DATASET'''
def calculate_FES(self,bins=50):
    Rb = 0.00198
    FES = {}
    FES_extent = {}
    phi = np.genfromtxt('phi')
    theta = np.genfromtxt('theta')
    for n in range(self.n_residues-1):
        z,x,y = np.histogram2d(phi[:,n], theta[:,n], bins=bins)
        FES[n] = -Rb*self.temp*np.ma.log(z.T)
        FES_extent[n] = [x.min(),x.max(),y.min(),y.max()]
    self.FES = FES
    self.FES_extent = FES_extent


def calculate_friction_coefficients(self):
    aKb = 1.38066e-23
    # Define solvent viscosity from fit to NIST at 1.0 atm
    vw = 0.2131590 - 1.96290e-3 * self.temp + \
        (0.00246411 * self.temp)**2 + (-0.0018462 * self.temp)**3

    # Define protein viscosity
    vp = vw * (1.23 * self._fD2O + self._fH2O) * self._internal_viscosity

    '''Check if Solvent Accessible Surface Area (SASA) is greater than or less than
	expected atomic surface area from VdW radii as defined by Miller Values. This
	can be used to determine the friction for proteins (fp) and solvent (fw).
	aw: Solvent Accessible Surface Area (SASA)
	av: VdW Surface Area from Miller values
	ap: Protein Accessible Surface Area
	'''

    if not hasattr(self, 'sasa'):
        SASA(self)
    if not hasattr(self, 'msa'):
        MSA(self)

    # Hydrophobic and solvent accessible surface area
    aw = self.sasa
    av = self.msa
    ap = np.zeros(self.n_residues)

    # Friction coefficients
    fp = np.zeros(self.n_residues)
    fw = np.zeros(self.n_residues)
    for n in range(self.n_residues):
        if aw[n] < av[n]:
            ap[n] = (av[n]**2 - aw[n]**2)**0.5
        fw[n] = 0.6 * np.pi * aw[n] * vw * (1.23 * self._fD2O + self._fH2O)
        fp[n] = 0.6 * np.pi * ap[n] * vp + fw[n]
    self.fw = fw
    self.fp = fp


def calculate_H_matrix(self):
    H = np.zeros((self.n_residues, self.n_residues), dtype=float, order='F')
    R = np.array(self._R, dtype=float, order='F')
    fratio = self.fw.mean() / self.fp.mean()
    _dynamics.calculate_h_matrix(H, R, self.fp, fratio, self.n_residues)
    self._H = H


def calculate_Hh_matrix(self, fratio):
    H = np.zeros((self.n_residues, self.n_residues), dtype=float, order='F')
    R = np.array(self._R, dtype=float, order='F')
    _dynamics.calculate_h_matrix(H, R, self.fp, fratio, self.n_residues)
    return H


def calculate_L_matrix(self):
    self._L = np.dot(self._a, np.dot(self._H, self._a.T))
    try:
        self._LI = np.linalg.inv(self._L)
    except:
        self._LI = np.linalg.pinv(self._L)
    self._LU = np.dot(self._UI, self._LI)


def calculate_M_matrix(self):
    M = np.zeros((self.n_residues, self.n_residues))
    for i in range(self.n_residues):
        M[i, i] = 1.0
    for i in range(self.n_residues - 1):
        M[i + 1, i] = -1.0
    M[0, :] = 1.0 / self.n_residues
    self._M = M


'''
Add to _dynamics.calculate_mode_traj to pass a histogram matrix.
'''
def calculate_mode_trajectory(self):
    QI = np.array(self._QI, dtype=float, order='F')
    if self._trajfile is None:
        traj = self._topfile
        top = None
    else:
        traj = self._trajfile
        top = self._topfile

    open('phi','w').close()
    open('theta','w').close()
    open('CA.xyz','w').close()
    for chunk in md.iterload(traj, top=top, chunk=self._chunk_size):
        chunk = prepare.solvent(chunk)
        chunk = prepare.atoms(chunk, self._skip_atoms)
        chunk = prepare.residues(chunk, self._skip_residues)

        # Set atomic slice range for alpha Carbons
        C_index = self.topology.select('name == CA')

        # Convert slices into ndarray
        C = np.array(chunk.xyz[:, C_index, :], dtype=float, order='F')

        _dynamics.calculate_mode_traj(C, QI, chunk.n_frames, self.n_residues)


def calculate_MSA(self):
    miller = {'ALA': (113.0 / (4 * np.pi))**0.5,
              'ARG': (241.0 / (4 * np.pi))**0.5,
              'ASN': (158.0 / (4 * np.pi))**0.5,
              'ASP': (151.0 / (4 * np.pi))**0.5,
              'CYS': (140.0 / (4 * np.pi))**0.5,
              'GLN': (189.0 / (4 * np.pi))**0.5,
              'GLU': (183.0 / (4 * np.pi))**0.5,
              'GLY': (85.0 / (4 * np.pi))**0.5,
              'HIS': (194.0 / (4 * np.pi))**0.5,
              'ILE': (182.0 / (4 * np.pi))**0.5,
              'LEU': (180.0 / (4 * np.pi))**0.5,
              'LYS': (211.0 / (4 * np.pi))**0.5,
              'MET': (204.0 / (4 * np.pi))**0.5,
              'PHE': (218.0 / (4 * np.pi))**0.5,
              'PRO': (143.0 / (4 * np.pi))**0.5,
              'SER': (122.0 / (4 * np.pi))**0.5,
              'THR': (146.0 / (4 * np.pi))**0.5,
              'TRP': (259.0 / (4 * np.pi))**0.5,
              'TYR': (229.0 / (4 * np.pi))**0.5,
              'VAL': (160.0 / (4 * np.pi))**0.5}

    self.msa = np.squeeze([miller[residue] for residue in self.residues])


def calculate_NMR_observables(self):
    timescale = (len(self.P2[:,0]) - 1) / 1000

    T1 = np.zeros(self.n_residues - 1, dtype=float, order='F')
    T2 = np.zeros(self.n_residues - 1, dtype=float, order='F')
    NOE = np.zeros(self.n_residues - 1, dtype=float, order='F')

    P2 = np.array(self.P2[:,1:], dtype=float, order='F')
    time = np.array(self.P2[:,0], dtype=float, order='F')
    _dynamics.calculate_nmr_observables(
        T1, T2, NOE, P2, time, self._NH_factor, timescale, self.n_residues)

    self.T1 = T1
    self.T2 = T2
    self.NOE = NOE


def calculate_Q_matrix(self):
    ''' Tells if the minimum eigenvalue of LU**-1 is positive or negative.
    r$\lambda_{min} < 0 \rightarrow io = 0
    r$\lambda_{min} > 0 \rightarrow io = 1
    '''

    def eigen_decomposition(self, fratio):
        U = np.array(self._UI, dtype=float, order='F')
        H = np.array(calculate_Hh_matrix(self, fratio), dtype=float, order='F')
        L = np.array(np.dot(self._a, np.dot(H, self._a.T)),
                     dtype=float, order='F')
        try:
            LI = np.array(np.linalg.inv(L), dtype=float, order='F')
        except:
            LI = np.array(np.linalg.pinv(L), dtype=float, order='F')
        LU = np.array(np.dot(U, LI), dtype=float, order='F')

        # Calculate lambda eigenvalues and Q-matrix
        w_lambda, Q = np.linalg.eig(LU)

        # Sort by eigenvalues
        idx = np.argsort(abs(w_lambda))[::-1]
        w_lambda = np.array(1 / (w_lambda[idx].real), dtype=float, order='F')
        Q = np.array(Q[:, idx], dtype=float, order='F')
        try:
            QI = np.array(np.linalg.inv(Q), dtype=float, order='F')
        except:
            QI = np.array(np.linalg.pinv(Q), dtype=float, order='F')

        if w_lambda.min() < 0:
            io = 0.0
        else:
            io = 1.0

        # Calculate mu eigenvalues
        w_mu = np.zeros(self.n_residues - 1, dtype=float, order='F')
        _dynamics.calculate_nu_eigenvalues(w_mu, QI, U, self.n_residues)

        # Calculate nu eigenvalues
        w_nu = np.array(np.diag(np.dot(QI, np.dot(L, Q))),
                        dtype=float, order='F')

        return Q, QI, w_lambda, w_mu, w_nu, io

    fratio = 0.27
    for i in range(self._n_iter):
        # Prepare changing H, L, LI, LU matrices
        Q, QI, w_lambda, w_mu, w_nu, io = eigen_decomposition(self, fratio)
        if io == 1.0:
            fratio -= 0.003
        fratio -= 0.0015
    fratio = self.fw.mean() / self.fp.mean()
    Q, QI, w_lambda, w_mu, w_nu, io = eigen_decomposition(self, fratio)
    self._Q = Q
    self._QI = QI

'''TODO: The time units prepared by the calculate_P2 method, can be
unnecessarily small, and unevenly spaced. This could potentially be a bug and
should be looked into.
'''
def calculate_P2(self, timescale=4):
    # Averaged Squared Bond Length
    blsq = 1.002 * (self._average_bond_lengths.mean())**2

    # Prepare bond matrix in NH basis as f2py input
    NH = np.array(self._NH_matrix, dtype=float, order='F')

    # Prepare Q matrices as f2py input
    Q = np.array(self._Q, dtype=float, order='F')
    QI = np.array(self._QI, dtype=float, order='F')

    # Prepare eigenvalues as f2py input
    w_lambda = np.array(self.lambda_eigenvalues, dtype=float, order='F')
    w_mu = np.array(self.mu_eigenvalues, dtype=float, order='F')

    # Prepare barriers as f2py input
    barriers = np.zeros(self.n_residues - 1, dtype=float, order='F')

    # Prepare mode lengthscale as f2py input
    modelength = np.zeros(
        (self.n_residues - 1, self.n_residues - 1), dtype=float, order='F')

    # Prepare time parameters as f2py input
    tau = np.zeros(self.n_residues - 1, dtype=float, order='F')
    taum = np.zeros(self.n_residues - 1, dtype=float, order='F')
    time = np.zeros(1000 * timescale + 1, dtype=float, order='F')

    # Prepare P2 correlation function as f2py input
    P2 = np.zeros((1000 * timescale + 1, self.n_residues - 1),
                  dtype=float, order='F')

    _dynamics.calculate_p2(P2, time, tau, taum, barriers, modelength, w_lambda, Q, QI,
                      NH, blsq, w_mu, self.sigma, self.temp, self.n_residues, timescale)

    self.lambda_eigenvalues = w_lambda
    self.mu_eigenvalues = w_mu
    self.barriers = barriers
    self.mode_length = modelength.T
    self.tau = tau
    self.tau_m1 = taum
    self.P2 = np.insert(P2, 0, time, axis=1)


def calculate_rmsd(self, reference=0, atom_indices=None, precentered=False):
    rmsd = md.rmsd(self._mdtraj, self._mdtraj, reference,
                   atom_indices=atom_indices, precentered=precentered)
    self.rmsd = rmsd


def calculate_SASA(self, probe_radius=0.14, n_sphere_points=250, stride=None):
    self.sasa = np.zeros((self.n_residues))
    if self._trajfile is None:
        traj = self._topfile
        top = None
    else:
        traj = self._trajfile
        top = self._topfile
    for chunk in md.iterload(traj, top=top, chunk=self._chunk_size, stride=stride):
        chunk = prepare.solvent(chunk)
        chunk = prepare.atoms(chunk, self._skip_atoms)
        chunk = prepare.residues(chunk, self._skip_residues)
        sasa = md.shrake_rupley(chunk, mode='residue',
                                probe_radius=probe_radius,
                                n_sphere_points=n_sphere_points)
        self.sasa += ((sasa*100/(4*np.pi))**0.5).sum(axis=0)
    self.sasa = self.sasa/self.n_conformers


def calculate_sigma(self):
    aKb = 1.38066e-23
    blsq = 1.002 * (self._average_bond_lengths.mean())**2
    self.sigma = (3 * aKb * self.temp * 1e-12 /
                  blsq / 1e-18) / (self.fp.mean() * 1e-9)


def calculate_U_matrix(self,stride=None):
    contacts = np.zeros((self.n_residues, self.n_residues), dtype=float, order='F')
    R = np.zeros((self.n_residues, self.n_residues), dtype=float, order='F')
    MSF = np.zeros((self.n_residues, self.n_residues), dtype=float, order='F')
    UI = np.zeros((self.n_residues-1, self.n_residues-1), dtype=float, order='F')

    if self._trajfile is None:
        traj = self._topfile
        top = None
    else:
        traj = self._trajfile
        top = self._topfile

    for chunk in md.iterload(traj, top=top, chunk=self._chunk_size, stride=stride):
        chunk = prepare.solvent(chunk)
        chunk = prepare.atoms(chunk, self._skip_atoms)
        chunk = prepare.residues(chunk, self._skip_residues)

        # Set atomic slice range for alpha Carbons
        C_index = self.topology.select('name == CA')

        # Convert slices into ndarray
        C = np.array(chunk.xyz[:, C_index, :], dtype=float, order='F')

        _dynamics.calculate_u_matrix(contacts, MSF, R, UI, C,
                                self.temp, chunk.n_frames, self.n_residues)

    # Average Gamma contacts
    self._contacts = contacts

    # Average B-factors and Mean Squared Fluctuations
    Rb = 1.98720412e-3  # Boltzmann's constant in (kcal/mol*K)
    cc = (3 * Rb * self.temp) / 6
    self.MSF = MSF/self.n_conformers
    self.bfactors = 100*(8/3)*(np.pi**2)*((3*Rb*self.temp)/6)*np.diag(MSF)

    # Average R matrix
    self._R = R/self.n_conformers

    # Average UI (inversed U) matrix
    UI = UI/self.n_conformers

    '''NOTE: The factor of 100 is a conversion between angstroms in the pdb
	file to nm. Bond lengths are in nm'''
    # Average inverse U matrix across all conformers
    for i in range(self.n_residues - 1):
        for j in range(self.n_residues - 1):
            UI[i, j] /= (self._average_bond_lengths[i] *
                        self._average_bond_lengths[j] * 100)
    self._UI = UI

    # Average U matrix
    try:
        self._U = np.linalg.inv(UI)
    except:
        self._U = np.linalg.pinv(UI)


def save_modes_pdb(self, max_mode=10, max_conf=20, file=None):
    '''
	Attributes:
	-----------
	self: object
			Imports the protein.dynamics object.
		max_mode: int (Default: 10)
			Maximum number of modes to be saved to pdb file.
		max_conf: int
			Maximum number of conformers to be saved to pdb files. For cases
            where there are more available conformations than max_conf
            (typically the case for simulation data), this method will randomly
            select conformations.
	'''
    # Create atomic count array, with the number of Atoms Per Residue (n_atom_per_res)
    n_atom_per_res = np.squeeze([len(self.topology.select('resid ' + str(i)))
                      for i in range(self.n_residues)])
    n_atom_per_res = np.array(n_atom_per_res, dtype=int, order='F')

    # Prepare input mode length as fortran contiguous array
    mode_length_in = np.array(self.mode_length, dtype=float, order='F')

    # Create output array as mode length per atomic site
    mode_length_out = np.zeros((self.n_residues - 1, self.n_atoms),
                        dtype=float, order='F')

    # Calculate output mode length
    _dynamics.calculate_mode_length_array(
        mode_length_out, mode_length_in, n_atom_per_res, self.n_atoms, self.n_residues)

    # Save atomic coordinates to pdb files with mode length as bfactor
    if max_mode > self.n_residues - 1:
        warnings.warn('''Requested number of copies exceeds the number of available modes
        ''')
        max_mode = self.n_residues - 1

    if max_conf < self.n_conformers:
        idx = np.random.randint(0,self.n_conformers, size=max_conf)
    for n in range(3, max_mode):
        if file is None:
            filename = 'mode_' + str(n + 1) + '.pdb'
        else:
            filename = file+'_'+str(n+1)+'.pdb'

        self._mdtraj.save_pdb(filename, bfactors=mode_length_out[n, :])
