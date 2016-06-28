import numpy as np
import mdtraj as md
import warnings
from LE4PD.ensembles.util import properties_util as util
from tqdm import trange


def calculate_a_matrix(self):
    a = np.zeros((self.n_residues - 1, self.n_residues),
                 dtype=float, order='F')
    M = np.array(self.M, dtype=float, order='F')
    util.calculate_a_matrix(a, M, self.n_residues)
    self.a = a


def calculate_bfactors(self):
    if not hasattr(self, 'bfactors'):
        Rb = 1.98720412e-3  # Boltzmann's constant in (kcal/mol*K)
        cc = (3 * Rb * self.temp) / 6
        self.bfactors = np.squeeze(
            [(8 / 3) * (np.pi**2) * cc * np.diag(self.MSF[n]) * 100 for n in range(self.n_conformers)])


def calculate_bond_vectors(self):
    """
    DESCRIBE INDEXING FOR THE THE ALPHA CARBONS, HYDROGENS, AND NITROGENS
    """

    # Slice alpha Carbons, Hydrogens, and Nitrogens from trajectory
    H_index = [atom.index for atom in self.top.atoms if(
        (str(atom.residue)[:3] == 'PRO' and atom.name == 'CD') or (atom.name == 'H'))]
    H = self._MD.atom_slice(H_index)
    C = self._MD.atom_slice(self._MD.top.select('name == CA'))
    N = self._MD.atom_slice(self._MD.top.select('name == N'))

    # Convert mdtraj slices into ndarray
    H = np.array(H.xyz, dtype=float, order='F')
    C = np.array(C.xyz, dtype=float, order='F')
    N = np.array(N.xyz, dtype=float, order='F')

    # Create alpha carbon bond vectors and magnitudes
    C_bond = np.zeros((self.n_conformers, self.n_residues - 1, 3),
                      dtype=float, order='F')
    C_length = np.zeros((self.n_conformers, self.n_residues - 1),
                        dtype=float, order='F')
    util.calculate_bond_vectors(
        C_bond, C_length, C, self.n_conformers, self.n_residues)

    # Create NH bond vectors and magnitudes
    NH_bond = np.zeros((self.n_conformers, self.n_residues - 1, 3),
                       dtype=float, order='F')
    NH_length = np.zeros((self.n_conformers, self.n_residues - 1),
                         dtype=float, order='F')
    util.calculate_nh_vectors(
        NH_bond, NH_length, N, H, self.n_conformers, self.n_residues)

    # Project bond dynamics on NH bond vectors
    NH_matrix = np.empty((self.n_conformers, self.n_residues - 1,
                          self.n_residues - 1), dtype=float, order='F')
    util.calculate_nh_matrix(
        NH_matrix, C_bond, NH_bond, C_length, NH_length, self.n_conformers, self.n_residues)

    self._NH_bonds = NH_bond
    self._NH_bond_length = NH_length
    self._NH_matrix = NH_matrix
    self._bonds = C_bond
    self._bond_length = C_length
    self._average_bond_length = C_length.mean(axis=0)


def calculate_eigenvalues(self):
    # Calculate lambda eigenvalues and Q-matrix
    w_lambda, Q = np.linalg.eig(self.LU)

    # Sort by eigenvalues
    idx = np.argsort(abs(w_lambda))[::-1]
    w_lambda = np.array(1/(w_lambda[idx].real), dtype=float, order='F')
    Q = np.array(Q[:, idx], dtype=float, order='F')
    QI = np.array(np.linalg.inv(Q), dtype=float, order='F')

    if w_lambda.min() < 0:
        io = 0.0
    else:
        io = 1.0

    # Calculate mu eigenvalues
    w_mu = np.zeros(self.n_residues-1, dtype=float, order='F')
    U = np.array(self.U, dtype=float, order='F')
    util.calculate_nu_eigenvalues(w_mu, QI, U, self.n_residues)

    # Calculate nu eigenvalues
    w_nu = np.array(np.diag(np.dot(QI, np.dot(self.L, Q))), dtype=float, order='F')

    self.lambda_eigenvalues = w_lambda
    self.mu_eigenvalues = w_mu
    self.nu_eigenvalues = w_nu


def calculate_friction_coefficients(self):
    aKb = 1.38066e-23
    # Define solvent viscosity from fit to NIST at 1.0 atm
    vw = 0.2131590 - 1.96290e-3 * self.temp + \
        (0.00246411 * self.temp)**2 + (-0.0018462 * self.temp)**3

    # Define protein viscosity
    vp = vw * (1.23 * self._fD2O + self._fH2O) * self._internal_viscosity

    """Check if Solvent Accessible Surface Area (SASA) is greater than or less than
	expected atomic surface area from VdW radii as defined by Miller Values. This
	can be used to determine the friction for proteins (fp) and solvent (fw).

	aw: Solvent Accessible Surface Area (SASA)
	av: VdW Surface Area from Miller values
	ap: Protein Accessible Surface Area
	"""

    if not hasattr(self, 'sasa'):
        calculate_SASA(self)
    if not hasattr(self, 'msa'):
        calculate_MSA(self)

    # Hydrophobic and solvent accessible surface area
    aw = self.sasa.mean(axis=0)
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


def calculate_gamma_contacts(self):
    if not hasattr(self, '_R'):
        calculate_R_matrix(self)
    # The gamma contacts for the Gaussian Network Model
    contacts = np.zeros((self.n_conformers, self.n_residues,
                         self.n_residues), dtype=float, order='F')
    eigenvalues = np.zeros(
        (self.n_conformers, self.n_residues), dtype=float, order='F')
    eigenvectors = np.zeros(
        (self.n_conformers, self.n_residues, self.n_residues), dtype=float, order='F')
    util.calculate_contacts(contacts, eigenvalues, eigenvectors, np.array(
        self.R, order='F'), self.n_conformers, self.n_residues)

    self._gamma_contacts = contacts
    self._gamma_eigenvalues = eigenvalues
    self._gamma_eigenvectors = eigenvectors


def calculate_H_matrix(self):
    H = np.zeros((self.n_residues, self.n_residues), dtype=float, order='F')
    R = np.array(self.R.mean(axis=0), dtype=float, order='F')
    fratio = self.fw.mean() / self.fp.mean()
    util.calculate_h_matrix(H, R, self.fp, fratio, self.n_residues)
    self.H = H


def calculate_H(self, fratio):
    H = np.zeros((self.n_residues, self.n_residues), dtype=float, order='F')
    R = np.array(self.R.mean(axis=0), dtype=float, order='F')
    util.calculate_h_matrix(H, R, self.fp, fratio, self.n_residues)
    return H


def calculate_L_matrix(self):
    self.L = np.dot(self.a, np.dot(self.H, self.a.T))
    self.LI = np.linalg.inv(self.L)
    self.LU = np.dot(self.U,self.LI)

def calculate_M_matrix(self):
    M = np.zeros((self.n_residues, self.n_residues))
    for i in range(self.n_residues):
        M[i, i] = 1.0
    for i in range(self.n_residues - 1):
        M[i + 1, i] = -1.0
    M[0, :] = 1.0 / self.n_residues
    self.M = M


def calculate_MSA(self):
    miller = {"ALA": (113.0 / (4 * np.pi))**0.5,
              "ARG": (241.0 / (4 * np.pi))**0.5,
              "ASN": (158.0 / (4 * np.pi))**0.5,
              "ASP": (151.0 / (4 * np.pi))**0.5,
              "CYS": (140.0 / (4 * np.pi))**0.5,
              "GLN": (189.0 / (4 * np.pi))**0.5,
              "GLU": (183.0 / (4 * np.pi))**0.5,
              "GLY": ( 85.0 / (4 * np.pi))**0.5,
              "HIS": (194.0 / (4 * np.pi))**0.5,
              "ILE": (182.0 / (4 * np.pi))**0.5,
              "LEU": (180.0 / (4 * np.pi))**0.5,
              "LYS": (211.0 / (4 * np.pi))**0.5,
              "MET": (204.0 / (4 * np.pi))**0.5,
              "PHE": (218.0 / (4 * np.pi))**0.5,
              "PRO": (143.0 / (4 * np.pi))**0.5,
              "SER": (122.0 / (4 * np.pi))**0.5,
              "THR": (146.0 / (4 * np.pi))**0.5,
              "TRP": (259.0 / (4 * np.pi))**0.5,
              "TYR": (229.0 / (4 * np.pi))**0.5,
              "VAL": (160.0 / (4 * np.pi))**0.5}

    self.msa = np.squeeze([miller[residue] for residue in self.residues])


def calculate_MSF(self):
    Rb = 1.98720412e-3  # Boltzmann's constant in (kcal/mol*K)
    cc = (3 * Rb * self.temp) / 6
    if not hasattr(self, '_gamma_contacts'):
        calculate_gamma_contacts(self)
    eigenvalues = np.array(self._gamma_eigenvalues, dtype=float, order='F')
    eigenvectors = np.array(self._gamma_eigenvectors, dtype=float, order='F')
    MSF = np.zeros((self.n_conformers, self.n_residues,
                    self.n_residues), dtype=float, order='F')
    util.calculate_msf(MSF, eigenvectors, eigenvalues, self.temp,
                                  self.n_conformers, self.n_residues)
    self.bfactors = np.squeeze(
        [(8 / 3) * (np.pi**2) * cc * np.diag(MSF[n]) * 100 for n in range(self.n_conformers)])
    self.MSF = MSF


def calculate_Q_matrix(self):
    """ Tells if the minimum eigenvalue of LU**-1 is positive or negative.
    r$\lambda_{min} < 0 \rightarrow io = 0
    r$\lambda_{min} > 0 \rightarrow io = 1
    """

    def eigen_decomposition(self, fratio):
        U = np.array(self.U, dtype=float, order='F')
        H = np.array(calculate_H(self, fratio), dtype=float, order='F')
        L = np.array(np.dot(self.a, np.dot(H, self.a.T)), dtype=float, order='F')
        LI = np.array(np.linalg.inv(L), dtype=float, order='F')
        LU = np.array(np.dot(U,LI), dtype=float, order='F')

        # Calculate lambda eigenvalues and Q-matrix
        w_lambda, Q = np.linalg.eig(LU)

        # Sort by eigenvalues
        idx = np.argsort(abs(w_lambda))[::-1]
        w_lambda = np.array(1/(w_lambda[idx].real), dtype=float, order='F')
        Q = np.array(Q[:,idx], dtype=float, order='F')
        QI = np.array(np.linalg.inv(Q), dtype=float, order='F')

        if w_lambda.min() < 0:
            io = 0.0
        else:
            io = 1.0

        # Calculate mu eigenvalues
        w_mu = np.zeros(self.n_residues-1, dtype=float, order='F')
        util.calculate_nu_eigenvalues(w_mu, QI, U, self.n_residues)

        # Calculate nu eigenvalues
        w_nu = np.array(np.diag(np.dot(QI, np.dot(L, Q))), dtype=float, order='F')

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
    self.Q = Q
    self.QI = QI

def calculate_P2(self, order=4):
    # Prepare bond matrix in NH basis as f2py input
    NH = np.array(self._NH_matrix.mean(axis=0), dtype=float, order='F')

    # Prepare Q matrices as f2py input
    Q = np.array(self.Q, dtype=float, order='F')
    QI = np.array(self.QI, dtype=float, order='F')

    # Prepare eigenvalues as f2py input
    w_lambda = np.array(self.lambda_eigenvalues, dtype=float, order='F')
    w_mu = np.array(self.mu_eigenvalues, dtype=float, order='F')

    # Prepare barriers as f2py input
    barriers = np.zeros(self.n_residues-1, dtype=float, order='F')

    # Prepare mode lengthscale as f2py input
    modelength = np.zeros((self.n_residues-1, self.n_residues-1), dtype=float, order='F')

    # Prepare time parameters as f2py input
    tau = np.zeros(self.n_residues-1, dtype=float, order='F')
    taum = np.zeros(self.n_residues-1, dtype=float, order='F')
    time = np.zeros(1000*order+1, dtype=float, order='F')

    # Prepare P2 correlation function as f2py input
    P2 = np.zeros((1000*order+1, self.n_residues-1), dtype=float, order='F')

    util.calculate_p2(P2, time, tau, taum, barriers, modelength, w_lambda, Q, QI,
                        NH, self._blsq, w_mu, self.sigma, self.temp, self.n_residues, order)

    self.lambda_eigenvalues = w_lambda
    self.mu_eigenvalues = w_mu
    self.barriers = barriers
    self.mode_length = modelength.T
    self.tau = tau
    self.tau_m1 = taum
    self.P2_time = time
    self.P2 = P2

def calculate_R_matrix(self):
    # Slice alpha carbons from trajectory
    C = self._MD.atom_slice(self._MD.top.select('name == CA'))
    # Convert mdtraj slices to ndarray
    C = np.array(C.xyz, dtype=float, order='F')
    R = np.zeros((self.n_conformers, self.n_residues,
                  self.n_residues), dtype=float, order='F')
    util.calculate_r_matrix(
        R, C, self.n_conformers, self.n_residues)
    self.R = R


def calculate_SASA(self, probe_radius=0.14, n_sphere_points=960):
    sasa = md.shrake_rupley(self._MD, mode='residue',
                            probe_radius=probe_radius, n_sphere_points=n_sphere_points)
    self.sasa = (sasa*100 / (4*np.pi))**0.5


def calculate_sigma(self):
    aKb = 1.38066e-23
    blsq = 1.002 * (self._bond_length.mean())**2
    self._blsq = blsq
    self.sigma = (3 * aKb * self.temp * 1e-12 /
                  blsq / 1e-18) / (self.fp.mean() * 1e-9)


def calculate_U_matrix(self):
    U = np.zeros((self.n_conformers, self.n_residues - 1,
                  self.n_residues - 1), dtype=float, order='F')
    if not hasattr(self, 'MSF'):
        calculate_MSF(self)
    util.calculate_u_matrix(U, np.array(self._bonds * 10, dtype=float, order='F'),
                                       np.array(
        self.MSF, dtype=float, order='F'),
        self.temp, self.n_conformers, self.n_residues)
    """NOTE: The fortran code actually returns the inverse U matrix. This is
	necessary for the LUI codes. For this reason, self.U is the TRUE matrix and
	self._U is the inverse matrix."""
    self._U = U         # INVERSE U MATRIX

    """NOTE: The factor of 100 is a conversion between angstroms in the pdb
	file to nm. Bond lengths are in nm"""
    # Average inverse U matrix across all conformers
    U = self._U.mean(axis=0)
    for i in range(self.n_residues - 1):
        for j in range(self.n_residues - 1):
            U[i, j] /= (self._average_bond_length[i] *
                        self._average_bond_length[j] * 100)
    self.U = U

def save_modes_pdb(self, max_modes=10):
	# Create atomic count array, with the number of Atoms Per Residue (apr)
    apr = np.squeeze([len(self.top.select('resid '+str(i))) for i in range(self.n_residues)])
    apr = np.array(apr, dtype=int, order='F')

    # Prepare input mode length as fortran contiguous array
    mlen_in = np.array(self.mode_length, dtype=float, order='F')

    # Create output array as mode length per atomic site
    mlen_out = np.zeros((self.n_residues-1,self.n_atoms), dtype=float, order='F')

    # Calculate output mode length
    util.calculate_mode_length_array(mlen_out, mlen_in, apr, self.n_atoms, self.n_residues)

    # Save atomic coordinates to pdb files with mode length as bfactor
    if max_modes > self.n_residues-1:
        warnings.warn("""Requested number of copies exceeds the number of available modes""")
        max_modes = self.n_residues-1
    try:
        for n in trange(3, max_modes):
            filename = "mode_"+str(n+1)+".pdb"
            self._MD.save_pdb(filename, bfactors=mlen_out[n,:])
    except:
        for n in range(3, max_modes):
            filename = "mode_"+str(n+1)+".pdb"
            self._MD.save_pdb(filename, bfactors=mlen_out[n,:])
