#!/bin/bash
echo "Mode number:"
read mode
echo $mode > mode
echo "#FE along path" > pathFE_${mode}.dat
filefe=../fe_${mode}.dat
echo $filefe
theta="84."
phi="160."
filest=../av_m${mode}_theta${theta}_phi${phi}
echo $filefe > filefe
echo $filest > filest
#filest1=""
points=""
rm points_m${mode}.dat
gfortran conv_rect.f -o c.exe
/usr/bin/gnuplot <<-EOF
set terminal wxt size 1300,1000
unset key
set pm3d map
splot "${filefe}" with pm3d
#splot "${filefe}" with pm3d; set mouse; pause mouse
set mouse
print "Click on origin of FE map"
pause mouse
print "Keystroke ", MOUSE_KEY, " at ", MOUSE_X, " ", MOUSE_Y
set print "xmin"
print MOUSE_X
set print "ymin"
print MOUSE_Y
EOF
/usr/bin/gnuplot <<-EOF
set terminal wxt size 1300,1000
unset key
#splot "${filest}" w l
set pm3d map
splot "${filefe}" with pm3d
#splot "${filefe}" with pm3d; pause mouse
set mouse
print "Click on upper right of FE map"
pause mouse
print "Keystroke ", MOUSE_KEY, " at ", MOUSE_X, " ", MOUSE_Y
set print "xmax"
print MOUSE_X
set print "ymax"
print MOUSE_Y
EOF

while [ $mode -le 75 ]
do
echo $mode > mode
filefe=../fe_${mode}.dat
echo "CLICK ON EQUILIBRIUM"
/usr/bin/gnuplot <<-EOF
set terminal wxt size 1300,1000
unset key
set pm3d map
set contour both
unset clabel
set cntrparam levels incremental 0.0, .15, 1.8
set cbrange[0:3.6]
splot "${filefe}" with pm3d linecolor rgb 'gray'
#splot "${filefe}" with pm3d; pause mouse
set mouse
pause mouse
print "Keystroke ", MOUSE_KEY, " at ", MOUSE_X, " ", MOUSE_Y
set print "x"
print MOUSE_X
set print "y"
print MOUSE_Y
EOF

./c.exe
read theta < theta
read phi < phi
Result=`echo "$phi > 365.0" | bc`
if [ $Result -eq 0 ]
then
read feline < feline
head -n ${feline} ${filefe} | tail -n 1 -f > FEmin_${mode}.dat

echo "click on barrier"
/usr/bin/gnuplot <<-EOF
set terminal wxt size 1300,1000
unset key
set pm3d map
set contour both
unset clabel
set cntrparam levels incremental 0.0, .15, 1.8
set cbrange[0:3.6]
splot "${filefe}" with pm3d linecolor rgb 'gray'${points}
#splot "${filefe}" with pm3d${points}; pause mouse
set mouse
pause mouse
print "Keystroke ", MOUSE_KEY, " at ", MOUSE_X, " ", MOUSE_Y
set print "x"
print MOUSE_X
set print "y"
print MOUSE_Y
EOF
else
mode=1000
fi

./c.exe
read theta < theta
read phi < phi
Result=`echo "$phi > 365.0" | bc`
if [ $Result -eq 0 ]
then
read feline < feline
head -n ${feline} ${filefe} | tail -n 1 -f > FEbarrier_${mode}.dat

mode=$(( $mode+1 ))
fi

done

