#!/bin/bash

BD=$PWD

pip install -e .
cd LE4PD/util/
rm -rfv __pycache__
rm -rfv *.so

f2py -m _ensemble_dynamics --fcompiler=gfortran --link-llapack -c _ensemble_dynamics.f95
f2py -m _simulation_dynamics --fcompiler=gfortran --link-llapack -c _simulation_dynamics.f95

f2py --quiet -m _m1int --fcompiler=gfortran -c p2_writeCAint.f
f2py --quiet -m _m1rot --fcompiler=gfortran -c p2_writeCArot.f

f2py --quiet -m _m1int_cmp --fcompiler=gfortran -c p2_writeCAint_cmp.f
f2py --quiet -m _m1rot_cmp --fcompiler=gfortran -c p2_writeCArot_cmp.f

cd $BD
