#! /bin/bash

#Identify the unique atom of interest's position and name:
unatom=`grep unatom FFbuild.in |awk '{print $2}'`
pdbnam=`grep pdbnam FFbuild.in |awk '{print $2}'`
atype=$unatom
ligname=mimic
num=`grep 'natoms' FFbuild.in | awk '{print $2 + 1}'`
grep 'Atom  No    Charge' *.log -A$num -m1 > charge.tmp 
echo ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
echo 'Unique atom of interest is:' $atype

# Run antechamber/pamchk on the mol2 to print out all bonds, angles, and dihedral terms 
# necessary.
parmchk -f mol2 -i mimic.mol2 -o mimic.frcmod -a Y -pf Y
natoms=`grep natoms FFbuild.in |awk '{print $2}'`

# Extract the charges that were previously being stored in a .tmp file
tail -$natoms charge.tmp > charge.clear

# Adjust map.py have the correct number of atoms
cp map.py map.run.py
sed  "s/natoms=XX/natoms=$natoms"/ map.py > map.run.py
sed -i "s/benz.pdb/$pdbnam/" map.run.py

python map.run.py

nun=`cat stats | awk '{print $1}'`
nb=`cat stats | awk '{print $2}'`
na=`cat stats | awk '{print $3}'`
nd=`cat stats | awk '{print $4}'`
nimp=`cat stats | awk '{print $5}'`


#echo $nun $nb $na $nd $nimp
# Asjust sort.py to have the correct number of atoms
cp sort.py sort.run.py
sed  -i s/"AOI='XX'"/"AOI='$atype'"/g sort.run.py             
sed  -i s/"natoms=XX"/"natoms=$natoms"/ sort.run.py                
sed  -i s/"totalbond=XX"/"totalbond=$nb"/ sort.run.py               
sed  -i s/"totalangle=XX"/"totalangle=$na"/ sort.run.py            
sed  -i s/"totaldihedral=XX"/"totaldihedral=$nd"/ sort.run.py      
sed  -i s/"totalimproper=XX"/"totalimproper=$nimp"/ sort.run.py   
sed  -i s/"totalvdw=XX"/"totalvdw=$nun"/ sort.run.py               

#Remove the atom of interest from the ffappend.tmp to enable sorting
grep -w $atype ffappend.tmp > specificID.txt
grep -v -w $atype ffappend.tmp > force_field.tmp

# Creates a force_field.dat
python ./sort.run.py

sed -i "s/BONDS/BONDS     $nb/" force_field.dat
sed -i "s/ANGLES/ANGLES     $na/" force_field.dat
sed -i "s/DIHEDRALS/DIHEDRALS     $nd/" force_field.dat
sed -i "s/IMPROPERS/IMPROPERS     $nimp/" force_field.dat
sed -i "s/VDW/VDW     $nun/" force_field.dat

#Moves the appropriate DLPOLY CONFIG files in place
cp ATOM_LIST GAUSSIAN refconfiginput ../fitting/dlpoly/ref_config

#Clean up
rm *run.py

##  Will output a force_field.dat and cutoff.dat with the appropriate gaff terms in 

echo ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
echo "Now check over the force_field.dat file and verify the atom names and parameters are all appropriate."
echo ""
echo "Also, a file called cutoffs.dat was created. "
echo ""
echo "This file has a default value of 1.6 Angstroms for each atom-atom pair"
echo "and will search for atoms that are less than 1.6 Angstroms apart and assume they are bonded."
echo ""
echo "Update the values as necessary for each bonding type if you expect to have longer bonds (e.g. metal-oxygen bonds)."
echo ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"












