# Langevin Equation for Protein Dynamics

This API provides a Python interface to use the updated LE4PD fortran codes developed by __Jeremy Copperman__, and __Marina Guenza__. The codes use a coase-grained Langevin formalism to obtain site-specific dynamics of proteins. This formalism takes the structural ensemble available to the protein as input, and generates the relaxation dynamics analytically.

# Dependencies
The following "non-standard" Python libraries MUST be installed for this code to work correctly:
* physt (https://github.com/janpipek/physt)

# Installation
The installation can be done simply with pip...

    pip install -e .

One can also use the 'install.sh' file, which will also complile the FORTRAN files in the LE4PD/util directory.

# Requirements
To run the analysis, all that is required is a trajectory file in .g96 format and a .pdb file with the structure of the protein of interest. The .g96 file should have been processed to remove rigid body rotation and translational motions. These codes are set up to analyze MD trajectories of proteins generated using GROMACS (http://www.gromacs.org/). If you have a trajectory generated from another MD simulation engine, you can use a piece of software such as Open Babel (http://openbabel.org) to convert to .g96 format or 2) write your own code to convert to .g96 format. Future implementations might have to ability to accept other types of MD trajectory files

# Examples
Examples for using LE4PD are located in LE4PD/examples. Two examples are present: ensemble and trajectory. LE4PD_ensemble_example shows how to set up the analysis for an ensemble of NMR experimental data from local files or fetched from the [Protein Data Bank](http://www.rcsb.org). LE4PD_trajectory_example performs the analysis on a 1ns simulation of Ubiquitin with modified residues.

# Computational Performance

Running the full LE4PD analysis (including save the model to file) takes about 1 minute per 50 000 frames of trajectory data and scales approximately linearly up to 1 500 000 frames of trajectory data. These benchmarking data were generated running the code on the Comet supercomputer at teh San Diego Supercomputing Center (https://www.sdsc.edu/support/user_guides/comet.html). 

# Caveats
This implementation of the LE4PD code is still a work in progress. Right now, the code is only set up to perform an analysis of a trajectory of a single protein. That is, we have not yet added the ability to analyze protein complexes or an ensemble of structures from a PDB file. Coming soon! 

It is also highly probable that the f2py settings in the install.sh file will have to be tweaked based on individual machine architectures and installations. 

# References
* __J. Copperman__ and __M. G. Guenza__ _"Mode Localization in the Cooperative Dynamics of Protein Recognition"_" Biophysical Journal (2015) (submitted) arXiv:1509.08913.
* __J. Copperman__ and __M. G. Guenza__ _"Predicting protein dynamics from structural ensembles” Journal of Chemical Physics, Invited Contribution for Special Topics Issue on Coarse Graining of Macromolecules, Biopolymers, and Membranes"_" 143, 243131-12 (2015) doi: 10.1063/1.4935575.
* __J. Copperman__, __M. G. Guenza__ _"A Coarse-Grained Langevin Equation for Protein Dynamics: Global anisotropy and a mode approach to local complexity"_" J. Phys. Chem. B, Festschrift issue honoring Branka Ladanyi, 119,  9195 (2015). DOI:10.1021/jp509473z.
* __E. Caballero-Manrique__, __M. G. Guenza__ _"A Langevin Equation Approach to Bridge Different Timescales of Relaxation in Protein Dynamics"_" 51st Annual Meeting of the Biophysical Society, MAR 03-07, 2007. Biophysical Journal 377A (2007).
