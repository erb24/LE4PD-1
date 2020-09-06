import numpy as np
import subprocess
import os
from LE4PD.util import _m1int, _m1rot, _m1int_cmp, _m1rot_cmp

def m1int(protname, nres, sigma, temp, barriers, lambda_eig, mu_eig, Q, QINV, comp = False, nmol = 1):

	if not comp:
		_m1int.p2m1(protname, nres, sigma, temp, barriers, lambda_eig, mu_eig, Q, QINV) 
	else:
		_m1int_cmp.p2m1(protname, nres, nmol, sigma, temp, barriers, lambda_eig, mu_eig, Q, QINV)

	if os.path.exists('mode_analysis'):
		pass
	else:
		os.mkdir('mode_analysis/')

	path = 'mode_analysis/'
	for i in range(nres - nmol):
		subprocess.call("mv -v m1int_" + str(i+1) + " " + path, shell = True)
	subprocess.call("mv -v avm1int "  + path, shell = True)

def m1rot(protname, nres, sigma, temp, barriers, lambda_eig, mu_eig, Q, QINV, comp = False, nmol = 1):
	if not comp:
		_m1rot.p2m1(protname, nres, sigma, temp, barriers, lambda_eig, mu_eig, Q, QINV) 
	else:
		_m1rot_cmp.p2m1(protname, nres, nmol, sigma, temp, barriers, lambda_eig, mu_eig, Q, QINV)

	if os.path.exists('mode_analysis'):
		pass
	else:
		os.mkdir('mode_analysis/')

	path = 'mode_analysis/'
	for i in range(nres - nmol):
		subprocess.call("mv -v m1CArot_" + str(i+1) + " " + path, shell = True)
	subprocess.call("mv -v avm1CArot "  + path, shell = True)