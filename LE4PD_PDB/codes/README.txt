Codes with a _cmp in them are set up to run for a protein complex, which can be one protein or many proteins.

average_cmp.f03
average.f03
Calculates averages of quantities over the conformers.

entropy_cmp.f03
entropy.f03
Quick code which calculates the estimated configurational entropy from the quasiharmonic approximation.

error.f
calculates error to experimental NMR relaxation data, e.g. PROTNAME_T1_exp.dat

gather_cmp.f
gather.f
LUI_calc_cmp.f03
LUI_calc.f03
T1T2_NH_cmp.f
T1T2_NH.f
Pretty much the same as for the MD codes.

p2_lb_cmp.f
p2_lb.f
Calculates p2 and m1 functions, uses energy barriers calculated from mode length as eps*Lp.

p2_lbscale.f
Calculates p2 and m1 functions, uses energy barriers calculated from scaling of energy barriers with mode number.

pdb_calc_cmp.f03
pdb_calc.f03
Calculates local flexibility around each conformer using GNM, and also necessary inputs from PDB file.

pdb_mode_cmp.f03
pdb_modesum.f03
Takes LML and outputs into B-factor slot of PDB so you can visualize it. Only set to write for first 6 modes. In pymol, visualize with the following commands to make pictures like in "Mode localization in the cooperative dynamics of protein recognition" --
pymol PROTNAME_m*.pdb
hide lines, all
show cartoon, all
cartoon putty
spectrum b, rainbow
set cartoon_putty_transform, 7
set cartoon_putty_scale_min, .06
set all_states, on

pdb_split_cmp.f03
pdb_split.f03
Splits pdb into separate files for each molecule and conformer

read_rsa_cmp.f03
read_rsa.f03
reads surface area output from NACCESS output

