#! /bin/bash

if [ -z "$1" ]
  then
	echo ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
    echo "Tell me how many bonds, angles, and dihedrals you want to fit!"
	nb=`cat unstats | awk '{print $2}'`
	na=`cat unstats | awk '{print $3}'`
	nd=`cat unstats | awk '{print $4}'`
	echo ""
	echo "Recall there are a total of $nb bonds, $na angles, and $nd dihedrals that involve the atom of interest."
	echo ""
	nb=`cat stats | awk '{print $2}'`
	na=`cat stats | awk '{print $3}'`
	nd=`cat stats | awk '{print $4}'`
	echo "And there are a total of $nb bonds, $na angles, and $nd dihedrals overall."
	echo ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
	exit 0
fi
cp cutoffs.dat cutoffs.bak
nb=$1
na=$2
nd=$3
# Create a nparm file for gamaker.f90
tailer=$(($nb*3+$na*2+$nd))
echo $tailer $nb $na $nd > nparm

#Run code that generates the parameter file for the GA and move the files

./specsort 

mv prepfield.f ../fitting/ga_code/
mv check.f ../fitting/

Npoints=`grep Npoints FFbuild.in | awk '{print $2}'`
stepsize=`grep stepsize FFbuild.in | awk '{print $2}'`
natoms=`grep natoms FFbuild.in | awk '{print $2}'`
pots=$(($Npoints*$Npoints))

let tmp=$Npoints-1
let tmp=$tmp/2
qmin=`echo "-1*$tmp*$stepsize" | bc`

cd ../fitting/ga_code/
make
cp fit_ff.exe ../
cd ../

cd dlpoly/md
echo $pots > ord.in
echo $Npoints >> ord.in
echo $qmin >> ord.in
echo $stepsize >> ord.in

cd ../../../setup/
cp input.DAT input
echo 'CRYST1  100.000  100.000  100.000  90.00  90.00 120.00' > liner
cat liner mimic.pdb > reference.pdb
rm liner
#This creates a FIELD_TMP with the dummy input
./build_ff_fit > stdout

echo TOTAL ATOMS  
grep atom stdout| awk '{print $3}'  
echo TOTAL BONDS  
grep bonds stdout| awk '{print $3}'
bond=`grep bonds stdout| awk '{print $3}'`
echo TOTAL ANGLES  
grep angles stdout| awk '{print $3}'
angle=`grep angles stdout| awk '{print $3}'`
echo TOTAL DIHEDRALS and IMPROPERS
grep impropers stdout| awk '{print $3}'
dihed=`grep impropers stdout| awk '{print $3}'`
echo TOTAL VDW  
grep vdw stdout| awk '{print $3}'
vdw=`grep vdw stdout| awk '{print $3}'`

sed -i "s/Nbond =  0/Nbond = $bond/" input.DAT
sed -i "s/Nangle = 0/Nangle = $angle/" input.DAT
sed -i "s/Ndihe =  0/Ndihe = $dihed/" input.DAT
sed -i "s/Nvdw =   0/Nvdw = $vdw/" input.DAT

cp input.DAT input
#This creates a FIELD_TMP with the real input
./build_ff_fit > stdout


cp force_field.dat ../fitting/dlpoly/field
cp input ../fitting/dlpoly/field
cp new_params.dat ../fitting
cp new_params.dat ../fitting/dlpoly/field/
cp reference.pdb ../fitting/dlpoly/field/
cp FIELD_TMP ../fitting/dlpoly/field/FIELD

./gamaker
