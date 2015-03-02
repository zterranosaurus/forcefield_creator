#! /bin/bash


Npoints=`grep Npoints FFbuild.in | awk '{print $2}'`
stepsize=`grep stepsize FFbuild.in | awk '{print $2}'`
natoms=`grep natoms FFbuild.in | awk '{print $2}'`
pots=$((121))

cd ../fitting/dlpoly/field/
bond=`grep Nbond_fit input| awk '{print $3}'`
angle=`grep Nangle_fit input| awk '{print $3}'`
dihed=`grep Ndihe_fit input| awk '{print $3}'`
#tailer=$(($bond*3+$angle*2+$dihed))
#echo $tailer $bond $angle $dihed > ../../nparm
tailer=$(($bond*3+$angle*2+$dihed))
echo $tailer $bond $angle $dihed > ../..//nparm
tailer=$(($tailer+11*$pots+5))
cp input ../../../setup/
cd ../../

tail -$tailer best_fit.txt > best_data
ifort check.f -o check
grep SSQ best_fit.txt |tail -1
slots=`grep SSQ best_fit.txt | wc -l`
echo "only $slots iterations have passed"
./check $pots $Npoints
cp pot.dat	../setup
cp ga.inp ../setup/ga.tmp
cp best_params.dat ../setup

cd ../setup
xmgrace -nxy pot.dat -viewport .15 .15 .85 .85 -world 0.0 0.0 100.0 20.0 -autoscale x -batch supp/bfile -hardcopy -printfile "fig.png" -hdevice PNG -saveall "Compare10.agr"

grep BONDS -A$bond force_field.dat | tail -$bond > paramsfitting
grep ANGLES -A$angle force_field.dat | tail -$angle >> paramsfitting
grep DIHEDRALS -A$dihed force_field.dat | tail -$dihed >> paramsfitting


python stitch.py
