# Langevin Equation for Protein Dynamics

This API provides a Python interface to use the updated LE4PD fortran codes developed by __Jeremy Copperman__, and __Marina Guenza__. The codes use a coase-grained Langevin formalism to obtain site-specific dynamics of proteins. This formalism takes the structural ensemble available to the protein as input, and generates the relaxation dynamics analytically.

# Installation
The installation can be done simply with pip...

    pip install -e .

The source files are written in Fortran, and while pre-compiled files are in LE4PD/util/, it is recommended the user compile these manually.

    f2py -m properties_util --fcompiler=gfortran --link-llapack -c properties_util.f95


# Examples
Examples for using LE4PD are located in LE4PD/examples. Two examples are present: ensemble and trajectory. LE4PD_ensemble_example shows how to set up the analysis for an ensemble of NMR experimental data from local files or fetched from the [Protein Data Bank](http://www.rcsb.org). LE4PD_trajectory_example performs the analysis on a 1ns simulation of Ubiquitin with modified residues.

# References
* __J. Copperman__ and __M. G. Guenza__ _"Mode Localization in the Cooperative Dynamics of Protein Recognition"_" Biophysical Journal (2015) (submitted) arXiv:1509.08913.
* __J. Copperman__ and __M. G. Guenza__ _"Predicting protein dynamics from structural ensembles” Journal of Chemical Physics, Invited Contribution for Special Topics Issue on Coarse Graining of Macromolecules, Biopolymers, and Membranes"_" 143, 243131-12 (2015) doi: 10.1063/1.4935575.
* __J. Copperman__, __M. G. Guenza__ _"A Coarse-Grained Langevin Equation for Protein Dynamics: Global anisotropy and a mode approach to local complexity"_" J. Phys. Chem. B, Festschrift issue honoring Branka Ladanyi, 119,  9195 (2015). DOI:10.1021/jp509473z.
* __E. Caballero-Manrique__, __M. G. Guenza__ _"A Langevin Equation Approach to Bridge Different Timescales of Relaxation in Protein Dynamics"_" 51st Annual Meeting of the Biophysical Society, MAR 03-07, 2007. Biophysical Journal 377A (2007).
