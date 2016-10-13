import numpy as np
import mdtraj as md
import warnings
from LE4PD.util import _ensemble_dynamics as util


def system(self, timescale=4, probe_radius=0.14, n_sphere_points=250):
    calculate_MSA(self)
    calculate_SASA(self, probe_radius=probe_radius, n_sphere_points=n_sphere_points)
    calculate_bond_vectors(self)
    calculate_R_matrix(self)
    calculate_MSF(self)
    calculate_U_matrix(self)
    calculate_friction_coefficients(self)
    calculate_sigma(self)
    calculate_H_matrix(self)
    calculate_M_matrix(self)
    calculate_a_matrix(self)
    calculate_L_matrix(self)
    calculate_Q_matrix(self)
    calculate_eigenvalues(self)
    calculate_aspherocity(self)
    calculate_P2(self, timescale=timescale)
    calculate_mode_trajectory(self)
    calculate_NMR_observables(self)


def calculate_a_matrix(self):
    a = np.zeros((self.n_residues - 1, self.n_residues),
                 dtype=float, order='F')
    M = self._M.astype(float, order='F')
    util.calculate_a_matrix(a, M, self.n_residues)
    self._a = a


def calculate_aspherocity(self):
    aspherocity = np.zeros(self.n_conformers)
    for n in range(self.n_conformers):
        aspherocity[n] = self._gamma_eigenvalues[n, 0] - 0.5 * \
            (self._gamma_eigenvalues[n, 1] + self._gamma_eigenvalues[n, 2])
    self._aspherocity = aspherocity


def calculate_bfactors(self):
    if not hasattr(self, 'bfactors'):
        Rb = 1.98720412e-3  # Boltzmann's constant in (kcal/mol*K)
        cc = (3 * Rb * self.temp) / 6
        self.bfactors = np.squeeze(
            [(8 / 3) * (np.pi**2) * cc * np.diag(self.MSF[n]) * 100
            for n in range(self.n_conformers)])


def calculate_bond_vectors(self):
    """
    DESCRIBE INDEXING FOR THE THE ALPHA CARBONS, HYDROGENS, AND NITROGENS
    """

    # Slice alpha Carbons, Hydrogens, and Nitrogens from trajectory
    H_index = [atom.index for atom in self.topology.atoms if(
        (str(atom.residue)[:3] == 'PRO' and atom.name == 'CD') or (atom.name == 'H'))]
    H = self._mdtraj.atom_slice(H_index)
    C = self._mdtraj.atom_slice(self._mdtraj.top.select('name == CA'))
    N = self._mdtraj.atom_slice(self._mdtraj.top.select('name == N'))

    # Convert mdtraj slices into ndarray
    H = H.xyz.astype(float, order='F')
    C = C.xyz.astype(float, order='F')
    N = N.xyz.astype(float, order='F')

    # Create alpha carbon bond vectors and magnitudes
    C_bond = np.zeros((self.n_conformers, self.n_residues - 1, 3))
    C_length = np.zeros((self.n_conformers, self.n_residues - 1))

    # Calculate bond dynamics
    util.calculate_bond_vectors(
        C_bond.astype(float, order='F'),
        C_length.astype(float, order='F'),
        C, self.n_conformers, self.n_residues)

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

    self._NH_bonds = NH_bond.mean(axis=0)
    self._NH_bond_length = NH_length
    self._NH_matrix = NH_matrix
    self._bonds = C_bond
    self._bond_length = C_length
    self._average_bond_length = C_length.mean(axis=0)


def calculate_CA_COM(self):
    C = self._mdtraj.atom_slice(self._mdtraj.top.select('name == CA'))
    C = np.array(C.xyz, dtype=float, order='F')
    com = np.zeros((self.n_conformers, 3), dtype=float, order='F')
    util.calculate_ca_com(com, C, self.n_conformers, self.n_residues)
    self.com = com


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
    UI = np.array(self._UI, dtype=float, order='F')
    util.calculate_nu_eigenvalues(w_mu, QI, UI, self.n_residues)

    # Calculate nu eigenvalues
    w_nu = np.array(np.diag(np.dot(QI, np.dot(self._L, Q))),
                    dtype=float, order='F')

    self.lambda_eigenvalues = w_lambda
    self.mu_eigenvalues = w_mu
    self.nu_eigenvalues = w_nu


def calculate_FES(self, bins=50):
    Rb = 0.00198
    FES = {}
    FES_extent = {}
    for n in range(self.n_residues - 1):
        z, x, y = np.histogram2d(
            self.modes[:, n, 0], self.modes[:, n, 1], bins=bins)
        FES[n] = -Rb * self.temp * np.ma.log(z.T)
        FES_extent[n] = [x.min(), x.max(), y.min(), y.max()]
    self.FES = FES
    self.FES_extent = FES_extent


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


def calculate_H_matrix(self):
    H = np.zeros((self.n_residues, self.n_residues), dtype=float, order='F')
    R = np.array(self._R, dtype=float, order='F')
    fratio = self.fw.mean() / self.fp.mean()
    util.calculate_h_matrix(H, R, self.fp, fratio, self.n_residues)
    self._H = H


def calculate_H(self, fratio):
    H = np.zeros((self.n_residues, self.n_residues), dtype=float, order='F')
    R = np.array(self._R, dtype=float, order='F')
    util.calculate_h_matrix(H, R, self.fp, fratio, self.n_residues)
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


def calculate_mode_trajectory(self, bins=100):
    bonds = np.array(self._bonds, dtype=float, order='F')
    QI = np.array(self._QI, dtype=float, order='F')
    modes = np.zeros((self.n_conformers, self.n_residues -
                      1, 2), dtype=float, order='F')
    util.calculate_mode_traj(
        modes, bonds, QI, self.n_conformers, self.n_residues)
    self.modes = modes


def calculate_MSA(self):
    miller = {"ALA": (113.0 / (4 * np.pi))**0.5,
              "ARG": (241.0 / (4 * np.pi))**0.5,
              "ASN": (158.0 / (4 * np.pi))**0.5,
              "ASP": (151.0 / (4 * np.pi))**0.5,
              "CYS": (140.0 / (4 * np.pi))**0.5,
              "GLN": (189.0 / (4 * np.pi))**0.5,
              "GLU": (183.0 / (4 * np.pi))**0.5,
              "GLY": (85.0 / (4 * np.pi))**0.5,
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


def calculate_NMR_observables(self):
    timescale = (len(self.P2[:, 0]) - 1) / 1000

    T1 = np.zeros(self.n_residues - 1, dtype=float, order='F')
    T2 = np.zeros(self.n_residues - 1, dtype=float, order='F')
    NOE = np.zeros(self.n_residues - 1, dtype=float, order='F')

    P2 = np.array(self.P2[:, 1:], dtype=float, order='F')
    time = np.array(self.P2[:, 0], dtype=float, order='F')
    util.calculate_nmr_observables(
        T1, T2, NOE, P2, time, self._NH_factor, timescale, self.n_residues)

    self.T1 = T1
    self.T2 = T2
    self.NOE = NOE


def calculate_Q_matrix(self):
    """ Tells if the minimum eigenvalue of LU**-1 is positive or negative.
    r$\lambda_{min} < 0 \rightarrow io = 0
    r$\lambda_{min} > 0 \rightarrow io = 1
    """

    def eigen_decomposition(self, fratio):
        UI = np.array(self._UI, dtype=float, order='F')
        H = np.array(calculate_H(self, fratio), dtype=float, order='F')
        L = np.array(np.dot(self._a, np.dot(H, self._a.T)),
                     dtype=float, order='F')
        try:
            LI = np.array(np.linalg.inv(L), dtype=float, order='F')
        except:
            LI = np.array(np.linalg.pinv(L), dtype=float, order='F')
        LU = np.array(np.dot(UI, LI), dtype=float, order='F')

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
        util.calculate_nu_eigenvalues(w_mu, QI, UI, self.n_residues)

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

"""TODO: The time units prepared by the calculate_P2 method, can be
unnecessarily small, and unevenly spaced. This could potentially be a bug and
should be looked into. Additionally, it makes no sense to have a separate
attribute (P2_time and P2). This should be concatenated such that P2_time is the
first column, and the remainder all correspond to a bond P2. This would give
P2.shape = (1000 * order + 1, self.n_residues)
"""


def calculate_P2(self, timescale=4):
    # Prepare bond matrix in NH basis as f2py input
    NH = np.array(self._NH_matrix.mean(axis=0), dtype=float, order='F')

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

    util.calculate_p2(P2, time, tau, taum, barriers, modelength, w_lambda, Q, QI,
                      NH, self._blsq, w_mu, self.sigma, self.temp, self.n_residues, timescale)

    self.lambda_eigenvalues = w_lambda
    self.mu_eigenvalues = w_mu
    self.barriers = barriers
    self.mode_length = modelength.T
    self.tau = tau
    self.tau_m1 = taum
    self.P2 = np.insert(P2, 0, time, axis=1)


def calculate_R_matrix(self):
    # Slice alpha carbons from trajectory
    C = self._mdtraj.atom_slice(self._mdtraj.top.select('name == CA'))
    # Convert mdtraj slices to ndarray
    C = np.array(C.xyz, dtype=float, order='F')
    R = np.zeros((self.n_conformers, self.n_residues,
                  self.n_residues), dtype=float, order='F')
    util.calculate_r_matrix(
        R, C, self.n_conformers, self.n_residues)

    # The gamma contacts for the Gaussian Network Model
    contacts = np.zeros((self.n_conformers, self.n_residues,
                         self.n_residues), dtype=float, order='F')
    eigenvalues = np.zeros(
        (self.n_conformers, self.n_residues), dtype=float, order='F')
    eigenvectors = np.zeros(
        (self.n_conformers, self.n_residues, self.n_residues), dtype=float, order='F')
    util.calculate_contacts(contacts, eigenvalues, eigenvectors, np.array(
        R, order='F'), self.n_conformers, self.n_residues)

    self._gamma_contacts = contacts
    self._gamma_eigenvalues = eigenvalues
    self._gamma_eigenvectors = eigenvectors
    self._R = R.mean(axis=0)


def calculate_rmsd(self, reference=0, atom_indices=None, precentered=False):
    rmsd = md.rmsd(self._mdtraj, self._mdtraj, reference,
                   atom_indices=atom_indices, precentered=precentered)
    self.rmsd = rmsd


def calculate_SASA(self, probe_radius=0.14, n_sphere_points=250):
    sasa = md.shrake_rupley(self._mdtraj, mode='residue',
                            probe_radius=probe_radius, n_sphere_points=n_sphere_points)
    self.sasa = (sasa * 100 / (4 * np.pi))**0.5


def calculate_sigma(self):
    aKb = 1.38066e-23
    blsq = 1.002 * (self._bond_length.mean())**2
    self._blsq = blsq
    self.sigma = (3 * aKb * self.temp * 1e-12 /
                  blsq / 1e-18) / (self.fp.mean() * 1e-9)


def calculate_U_matrix(self):
    UI = np.zeros((self.n_conformers, self.n_residues - 1,
                   self.n_residues - 1), dtype=float, order='F')
    if not hasattr(self, 'MSF'):
        calculate_MSF(self)
    util.calculate_u_matrix(UI, np.array(self._bonds * 10, dtype=float, order='F'),
                            np.array(
        self.MSF, dtype=float, order='F'),
        self.temp, self.n_conformers, self.n_residues)
    """NOTE: The fortran code actually returns the inverse U matrix. This is
	necessary for the LUI codes. For this reason, self._UI is the TRUE matrix and
	self._UI is the inverse matrix."""
    self._UI = UI         # INVERSE U MATRIX

    """NOTE: The factor of 100 is a conversion between angstroms in the pdb
	file to nm. Bond lengths are in nm"""
    # Average inverse U matrix across all conformers
    UI = self._UI.mean(axis=0)
    for i in range(self.n_residues - 1):
        for j in range(self.n_residues - 1):
            UI[i, j] /= (self._average_bond_length[i] *
                         self._average_bond_length[j] * 100)
    self._UI = UI

    # U matrix
    try:
        self._U = np.linalg.inv(UI)
    except:
        self._U = np.linalg.pinv(UI)


def save_modes_pdb(self, max_mode=10, max_conf=20):
    """
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
        """
    # Create atomic count array, with the number of Atoms Per Residue (apr)
    apr = np.squeeze([len(self.topology.select('resid ' + str(i)))
                      for i in range(self.n_residues)])
    apr = np.array(apr, dtype=int, order='F')

    # Prepare input mode length as fortran contiguous array
    mode_length_in = np.array(self.mode_length, dtype=float, order='F')

    # Create output array as mode length per atomic site
    mode_length_out = np.zeros((self.n_residues - 1, self.n_atoms),
                               dtype=float, order='F')

    # Calculate output mode length
    util.calculate_mode_length_array(
        mode_length_out, mode_length_in, apr, self.n_atoms, self.n_residues)

    # Save atomic coordinates to pdb files with mode length as bfactor
    if max_mode > self.n_residues - 1:
        warnings.warn("""Requested number of copies exceeds the number of available modes
        """)
        max_mode = self.n_residues - 1

    if max_conf < self.n_conformers:
        idx = np.random.randint(0, self.n_conformers, size=max_conf)
    for n in range(3, max_mode):
        filename = "mode_" + str(n + 1) + ".pdb"
        if max_conf >= self.n_conformers:
            self._mdtraj.save_pdb(filename, bfactors=mode_length_out[n, :])
        elif max_conf < self.n_conformers:
            self._mdtraj[idx].save_pdb(filename, bfactors=mode_length_out[n, :])
