#!/bin/bash

BD=$PWD

pip install -e .
cd LE4PD/util/
rm -rfv __pycache__
rm -rfv *.so

f2py -m _ensemble_dynamics --fcompiler=gfortran --link-llapack -c _ensemble_dynamics.f95
f2py -m _simulation_dynamics --fcompiler=gfortran --link-llapack -c _simulation_dynamics.f95

f2py -m _m1int --fcompiler=gfortran -c p2_writeCAint.f
f2py -m _m1rot --fcompiler=gfortran -c p2_writeCArot.f

cd $BD
