#!/bin/bash

BD=$PWD

pip install -e .
cd LE4PD/util/
f2py -m _ensemble_dynamics --fcompiler=gfortran --link-llapack -c _ensemble_dynamics.f95
f2py -m _simulation_dynamics --fcompiler=gfortran --link-llapack -c _simulation_dynamics.f95

cd $BD
