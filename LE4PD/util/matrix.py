import numpy as np

class matrix(object):
    """
    dynamics: object
            This imports the dynamics class, which contains the predicted
            dynamics for the molecule object using the LE4PD analysis.
    """
    def __init__(self, dynamics):
        self.a = dynamics._a
        self.H = dynamics._H
        self.L = dynamics._L
        self.LI = dynamics._LI
        self.LU = dynamics._LU
        self.M = dynamics._M
        self.NH = dynamics._NH_matrix
        self.Q = dynamics._Q
        self.QI = dynamics._QI
        self.R = dynamics._R
        self.U = dynamics._U
        self.UI = dynamics._UI
