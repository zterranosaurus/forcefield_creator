#! /bin/bash

num=`grep 'natoms' FFbuild.in | awk '{print $2 + 5}'`
npoints=`grep Npoints FFbuild.in |awk '{print $2}'`
stepsize=`grep stepsize FFbuild.in |awk '{print $2}'`
unatom=`grep unatom FFbuild.in |awk '{print $2}'`
natoms=`grep natoms FFbuild.in |awk '{print $2}'`
scan=`grep scan FFbuild.in |awk '{print $2}'`
aoipos=`grep aoi_pos FFbuild.in |awk '{print $2}'`
pdbnam=`grep pdbnam FFbuild.in |awk '{print $2}'`

i="1"
N=$((npoints*npoints*npoints))

grep 'Atom  No    Charge' *.log -A$num -m1 > charge.tmp
antechamber -fi gout -fo mol2 -i *.log -o mimic.mol2 -j 4 -s 2 -pf Y 
antechamber -fi gout -fo gcrt -i *.log -o tmp -j 4 -s 2 -pf Y
#This assumes the tmp file has 8 lines before the coordinates start
./findnongaffbonds $aoipos $natoms

./genconfig gaussian $npoints $stepsize tmpxyz $scan $natoms $aoipos


mv *.pbs set$scan
cd set$scan
 while [ $i -le $N ]; do
	cat ../gaussianINPUTS/ffHEADER config$i.com ../gaussianINPUTS/basis > config$i.tmp
	mv config$i.tmp config$i.com
    let i=$i+1
 done

i="1"
N=$((npoints*npoints*npoints))

while [ $i -le $N ]; do
	cat ../gaussianINPUTS/carverHEAD jobPES$i.pbs ../gaussianINPUTS/carverTAIL > jobPES$i.tmp
	mv jobPES$i.tmp jobPES$i.pbs
   let i=$i+100
done

cd ..
mv set$scan ../gaussian/scans/
cd ../gaussian/scans/
CWD=$(pwd)
echo 'The PES scans inputs are located :: '$CWD
echo ""

echo "Ready to scp set$scan over to carver.nersc.gov to run the gaussian jobs"
echo ""
cd ../../setup
cp POSITIONS ../fitting/dlpoly/md/

#cleanup
rm atommoved.tmp posmoved.tmp tmpxyz tmp


