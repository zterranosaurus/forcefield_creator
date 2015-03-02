#! /bin/bash

#Check and see if your system can run forcefield_creator
PYV=`python -c "import sys;t='{v[0]}.{v[1]}'.format(v=list(sys.version_info[:2]));sys.stdout.write(t)";`
echo "Python version "$PYV "installed"

tar -xvzf DLCLASS.tar.gz > install_notes.txt
rm install_notes.txt

mkdir -p fitting gaussian setup
mkdir -p setup/supp
mkdir -p fitting/dlpoly
mkdir -p fitting/dlpoly/field fitting/dlpoly/md fitting/dlpoly/ref_config
mkdir -p fitting/ga_code

mkdir -p gaussian
mkdir -p gaussian/optimization
mkdir -p gaussian/scans

pwd=`echo $PWD`

cd setup/
ifort gamaker.f90 -o gamaker
ifort specsort.f90 -o specsort
ifort findnongaffbonds.f90 -o findnongaffbonds
ifort build_ff_fit.f90 -o build_ff_fit
ifort genconfig.f90 -o genconfig
cp build_ff_fit ../fitting/dlpoly/field
cd ../

cd $pwd/gaussian/scans/
ifort compute_pot.f90 -o compute_pot

cd
cd $pwd/fitting/dlpoly/ref_config/
ifort make_config.f90 -o make_config
cd
cd $pwd/fitting/dlpoly/md
ifort ord_pot.f90 -o ord_pot
touch run_dlpoly
rm run_dlpoly

echo '#!/bin/bash' >>run_dlpoly
echo ' '>>run_dlpoly
echo 'export DO_ENG_FRC=1'>>run_dlpoly
echo "mpirun -np 1 $pwd/dl_class/execute/dl_class_mpi.exe">>run_dlpoly
echo ' '>>run_dlpoly
chmod +x run_dlpoly

#mkdir -p MOFsystems

echo "Abbreviated instructions:"
echo ""
echo "Move into the setup folder, execute:"
echo "FFinit.sh, createPES.sh, gaffprep.sh, createFIELD.sh"
