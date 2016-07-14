#!/bin/bash
echo 'load "fes.plt"' > lp
for a in 1 2 3 4 5 6 7 8 9 10
do
 echo $a
 cp fes.plt fes${a}.plt
 echo 'load "fes'"${a}.plt"'"' > lp
 str="s/NMODE/${a}/g"
 sed -i ${str} fes${a}.plt
 gnuplot < lp
 a=$(( a+1 ))
done
