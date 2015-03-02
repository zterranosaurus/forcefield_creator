#! /bin/bash

./build_ff_fit
cp input ../fitting/dlpoly/field/
cp new_params.dat ../fitting
cp new_params.dat ../fitting/dlpoly/field/
cp reference.pdb ../fitting/dlpoly/field/
cp FIELD_TMP ../fitting/dlpoly/field/FIELD
mv ga.tmp ../fitting/ga.inp


cd ../fitting/
touch ga.out ga.restart best_fit.txt	
rm ga.out ga.restart best_fit.txt
./fit_ff.exe
