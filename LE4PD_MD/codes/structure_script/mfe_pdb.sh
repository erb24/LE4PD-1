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
gfortran convert.f -o c.exe
gfortran pdbmaker.f -o pm.exe
/usr/bin/gnuplot <<-EOF
set terminal wxt size 1300,1000
set size 1.5,1.0
set origin 0.1,0.1
set multiplot
set size .5,.5
set origin 0.1,0.5
unset key
#splot "${filest}" w l;
set size .5,.5
set origin 0.1,0.1
set pm3d map
splot "${filefe}" with pm3d
#splot "${filefe}" with pm3d; set mouse; pause mouse
unset multiplot
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
set size 1.5,1.0
set origin 0.1,0.1
set multiplot
set size .5,.5
set origin 0.1,0.5
unset key
#splot "${filest}" w l
set size .5,.5
set origin 0.1,0.1
set pm3d map
splot "${filefe}" with pm3d
#splot "${filefe}" with pm3d; pause mouse
unset multiplot
set mouse
print "Click on upper right of FE map"
pause mouse
print "Keystroke ", MOUSE_KEY, " at ", MOUSE_X, " ", MOUSE_Y
set print "xmax"
print MOUSE_X
set print "ymax"
print MOUSE_Y
EOF

#use min and max coordinates to specify theta, phi values from clicks
/usr/bin/gnuplot <<-EOF
set terminal wxt size 1300,1000
set size 1.5,1.0
set origin 0.1,0.1
set multiplot
set size .5,.5
set origin 0.1,0.5
unset key
#splot "${filest}" w l
set size .5,.5
set origin 0.1,0.1
set pm3d map
splot "${filefe}" with pm3d
#splot "${filefe}" with pm3d; pause mouse
unset multiplot
set mouse
pause mouse
print "Keystroke ", MOUSE_KEY, " at ", MOUSE_X, " ", MOUSE_Y
set print "x"
print MOUSE_X
set print "y"
print MOUSE_Y
EOF
q=""
i=0
while [ "$q" != "q" ]
do
./c.exe
read theta < theta
read phi < phi
if [ -e ../av_m${mode}_theta${theta}_phi${phi} ]
then
echo "${theta}     ${phi}     0.0" > points_m${mode}_${i}.dat
echo "${theta}     ${phi}     0.0" >> points_m${mode}.dat
animate[$i]='set size 1.5,1.0; set origin 0.1,0.1; set multiplot; set size .5,.5; set origin 0.1,0.1; set xrange[0.0:180.0]; set yrange[0:360]; set zrange[-1.0:8]; set pm3d map; unset key; splot "'${filefe}'" with pm3d, "points_m'${mode}'_'${i}'.dat" with points linecolor rgb "white"; unset pm3d; set size .5,.5; set origin 0.1,0.5; set xrange[1.0:4.5]; set yrange[.3:6.0]; set zrange[1.2:5.0]; splot "../'av_m${mode}_theta${theta}_phi${phi}'" w l; unset multiplot'
filest1=$filest1', "../'av_m${mode}_theta${theta}_phi${phi}'" w l'
points=${points}', "points_m'${mode}'_'${i}'.dat" with points linecolor rgb "white"'
i=$(( $i+1 ))
echo 'conformation '$i
echo $i > frame
echo "../av_m${mode}_theta${theta}_phi${phi}" > framefile
./pm.exe
echo $filest1
read feline < feline
head -n ${feline} ${filefe} | tail -n 1 -f >> pathFE_${mode}.dat
/usr/bin/gnuplot <<-EOF
set terminal wxt size 1300,1000
set size 1.5,1.0
set origin 0.1,0.1
set multiplot
set size .5,.5
set origin 0.1,0.5
unset key
splot${filest1:1}
set size .5,.5
set origin 0.1,0.1
unset key
set pm3d map
set cbrange[0:3.6]
splot "${filefe}" with pm3d${points}
#splot "${filefe}" with pm3d${points}; pause mouse
unset multiplot
set mouse
pause mouse
print "Keystroke ", MOUSE_KEY, " at ", MOUSE_X, " ", MOUSE_Y
set print "x"
print MOUSE_X
set print "y"
print MOUSE_Y
EOF
else
Result=`echo "$phi > 500.0" | bc`
if [ $Result -eq 1 ]
then
echo "plotting..."
/usr/bin/gnuplot -persist <<-EOF
set terminal wxt size 1300,1000
unset key
unset xtics
unset ytics
unset ztics
unset border
splot${filest1:1}
save "conformations_m${mode}.plt"
EOF

/usr/bin/gnuplot -persist <<-EOF
set terminal wxt size 1500,1000
set origin 0.1,0.1
set size 1.5,1.0
set multiplot
set size .5,.5
set origin 0.0,0.15
unset key
set pm3d map interpolate 3,3
set xlabel "{/Symbol q} (deg)"
set ylabel "{/Symbol f} (deg)" 
set cblabel "kcal/mol" 
set cbrange[0:3.6]
splot "${filefe}" with pm3d, "points_m${mode}.dat" with points linecolor rgb "white", "points_m${mode}.dat" w l linecolor rgb "white"
set size .5,.5
set origin 0.5,0.15
set xrange[1.5:5.0]
set yrange[1.5:5.0]
set zrange[1.5:5.5]
unset key
unset pm3d
set xlabel ""
set ylabel ""
set zlabel ""
unset xtics
unset ytics
unset ztics
splot${filest1}
unset multiplot
save "mode_m${mode}.plt"
EOF

#write animation gnuplot script
j=0
echo "set terminal gif animate delay 10 size 1280,960" > animate_m${mode}.plt
echo 'set output "animate_m'${mode}'.gif"' >> animate_m${mode}.plt
while [ $j -le $i ]
do
echo ${animate[$j]} >> animate_m${mode}.plt
j=$(( $j+1 ))
done
/usr/bin/gnuplot <<-EOF
load "animate_m${mode}.plt"
EOF
animate animate_m${mode}.gif

echo ${filefe} > filefe
echo ${mode} > mode
echo ${filest1:1} > filest1

echo "-100.0" > y
echo "Enter q to quit:"
read q
else
echo "NO FILE EXISTS!"
/usr/bin/gnuplot <<-EOF
set terminal wxt size 1300,1000
set size 1.5,1.0
set origin 0.1,0.1
set multiplot
set size .5,.5
set origin 0.1,0.5
unset key
splot${filest1:1}
set size .5,.5
set origin 0.1,0.1
set pm3d map
splot "${filefe}" with pm3d${points}
#splot "${filefe}" with pm3d${points}; pause mouse
unset multiplot
set mouse
pause mouse
print "Keystroke ", MOUSE_KEY, " at ", MOUSE_X, " ", MOUSE_Y
set print "x"
print MOUSE_X
set print "y"
print MOUSE_Y
EOF
echo "q to quit"
read q
fi
fi
done

