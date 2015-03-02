#! /bin/bash

npoints=`grep Npoints FFbuild.in |awk '{print $2}'`
stepsize=`grep stepsize FFbuild.in |awk '{print $2}'`
unatom=`grep unatom FFbuild.in |awk '{print $2}'`
natoms=`grep natoms FFbuild.in |awk '{print $2}'`
scan=`grep scan FFbuild.in |awk '{print $2}'`

i="1"
N=$((npoints*npoints*npoints))

cd ../gaussian/scans/
touch target_pot.dat
rm target_pot.dat
cp compute_pot set$scan
cd set$scan
touch enrg.dat
rm enrg*.dat

 while [ $i -le $N ]; do
     enrg=`grep 'SCF Done' config$i.out | sed 's/ SCF Done:  E(UM062X) = //g' | awk '{print $1;}'`
     echo "$i  $enrg" >> enrg.dat
     sed -n "11,62 p"  config$i.com > config$i.xyz
     echo "config$i done..."
     let i=$i+1
  done
  mv enrg.dat ..
  cd ..
  echo $N > parm.dat
  echo $npoints >> parm.dat
  echo $stepsize >> parm.dat
  ./compute_pot
  mv pot.dat target_pot.tmp
  rm parm.dat
  N=$((npoints*npoints))
  echo $N $npoints >tmp
 cat tmp target_pot.tmp > target_pot.dat
 rm target_pot.tmp tmp
 mv target_pot.dat ../../fitting/
