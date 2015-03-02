#! /bin/bash

Npoints=`grep Npoints ../setup/FFbuild.in | awk '{print $2}'`
stepsize=`grep stepsize ../setup/FFbuild.in | awk '{print $2}'`
natoms=`grep natoms ../setup/FFbuild.in | awk '{print $2}'`
pots=$(($Npoints*$Npoints))

let tmp=$Npoints-1
let tmp=$tmp/2
qmin=`echo "-1*$tmp*$stepsize" | bc`

cd ../fitting/
cp new_params.dat dlpoly/field

cd dlpoly/field
./build_ff_fit
cp FIELD_TMP FIELD
cp FIELD ../md

cd ../ref_config
./make_config
cp CONFIG ../md

cd ../md

#echo $pots > ord.in
#echo $Npoints >> ord.in
#echo $qmin >> ord.in
#echo $stepsize >> ord.in

./run_dlpoly
./ord_pot
mv pot_dlpoly.dat ../..

cd ../..
echo "Updated FIELD with new parameters"
