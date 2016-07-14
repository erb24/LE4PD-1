#!/bin/bash
# executes bridget codes
rm protname.txt
rm nres*.dat
rm natoms*.dat
prot="2m6i"
echo $prot > prot
echo 298 > temp #temperature in kelvin
echo 0.0 > fd20 #fraction D20 content
echo 2.71828 > internalv #internal viscosity rescaling from solvent viscosity
echo $prot > protname.txt
gfortran ../codes/pdb_split_cmp.f03 -o pdb_split_cmp.exe
./pdb_split_cmp.exe #splits pdb file into individual conformers
gfortran ../codes/pdb_calc_cmp.f03 -o pdb_calc_cmp.exe -llapack -fcheck-bounds
gfortran ../codes/read_rsa_cmp.f03 -o read_rsa.exe
while read line
do
 nc=$line
 echo "$nc conformers"
done < "ncopies.dat"
while read line
do
 nmol=$line
done < "nmol.dat"
echo "$nmol molecules"
i=1
while [ $i -le $nc ]
do
 protname=${prot}${i}
 rm protname.txt
 echo $protname > protname.txt
 echo $protname
 ./pdb_calc_cmp.exe #calculates for each conformer (GNM,Umatrix,Rij,etc.)
 echo "done pdb_calc"
 j=1
 echo "COMPLEX" > complex.pdb
 while [ $j -le $nmol ]
 do
  cat ${prot}${i}_${j}.pdb >> complex.pdb
  echo "TER" >> complex.pdb
  j=$(( $j+1 ))
 done
 /home/jeremy/naccess2.1.1/naccess complex.pdb
 ./read_rsa.exe
 mv resrad resrad${i}
 mv Umatrix Umatrix${i}
 mv Rij Rij${i}
 mv bfactors.dat bfactors${i}.dat
 mv ldot.dat ldot${i}.dat
 mv avbl avbl${i}
 mv length length${i}
 i=$(( $i+1 ))
done

gfortran ../codes/average_cmp.f03 -o average_cmp.exe
./average_cmp.exe #calculate averages over conformers
gfortran ../codes/LUI_calc_cmp.f03 -o LUI_calc_cmp.exe -llapack

./LUI_calc_cmp.exe #friction coefficients,builds and diagonalizes matrices
gfortran ../codes/p2_lb_cmp.f -o p2_lb_cmp.exe
./p2_lb_cmp.exe
echo 1.0 > NHfactor.dat
gfortran ../codes/T1T2_NH_cmp.f -o T1T2_NH_cmp.exe
gfortran ../codes/gather_cmp.f -o gather_cmp.exe
./T1T2_NH_cmp.exe
./gather_cmp.exe
gfortran ../codes/pdb_mode_cmp.f03 -o pdb_mode_cmp.exe
./pdb_mode_cmp.exe
