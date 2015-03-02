#! /usr/bin/env python

# This program obtains the mass, charges, and atom types from the mol2 file
# and the initial pdb.  It also will compile all of the gaff information from 
# the .frcmod file and output it into a temporary file called ffappend.tmp.  
# This file contains all the force field parameters in a DL_POLY format.

import periodic as per
import math
#print per.table
natoms=XX

OUTPUT=open('ffappend.tmp', 'w+')
OUTPUT.write("%-8s\t %d\n"  % ('ATOMS', natoms))
OUTPUT2=open('cutoffs.dat','w+')


f=open('mimic.mol2')
f2=open('benz.pdb')
f3=open('charge.clear')
line = f2.readline()
for lines in range (0,8):
	line = f.readline()
for lines in range (0,natoms):
	line = f.readline()
	charge=f3.readline().split()
	q=float(charge[2])
	list1 = line.split()
	line2 = f2.readline()
	list2 = line2.split()
	atom_type=(list1[5])
	element = per.element(list2[10])
	mass=element.mass
	OUTPUT.write("%-2s \t%2s\t%12.4f\t%12.4f\n" % (atom_type,atom_type,q,mass))
	
f.close()
f2.close()
f3.close()

#
nun=0; nbond=0; nangle=0; ndihedral=0; nimproper=0; nvdw=0
f=open('mimic.frcmod')
line=f.readline()
line=f.readline()
line=f.readline().split()
#Read the .frcmod until it finds a blank line
while len(line) > 0:
	atom=(line[0])
	nun=nun+1
	line = f.readline().split()
	if line in ['\n', '\r\n']:
		break
line=f.readline()
line=f.readline().replace('-', ' ').split()
OUTPUT.write('\n%s\n' % 'BONDS')
while len(line) > 0:
	bond_1=(line[0])
	bond_2=(line[1])
	bond_k=float(line[2])
	bond_r=float(line[3])
	nbond=nbond+1
	OUTPUT.write("%-2s %2s \t%4s %12.4f%12.4f\n" % (bond_1,bond_2,'harm',bond_k*2.0,bond_r))
	OUTPUT2.write("%-2s %2s \t%f\n" % (bond_1,bond_2,1.6))
	line = f.readline().replace('-', ' ').split()
	if line in ['\n', '\r\n']:
		break
line=f.readline()
line=f.readline().replace('-', ' ').split()
OUTPUT.write('\n%s\n' % 'ANGLES')
while len(line) > 0:
	nangle=nangle+1
	angle_1=(line[0])
	angle_2=(line[1])
	angle_3=(line[2])
	angle_k=float(line[3])
	angle_r=float(line[4])
	OUTPUT.write("%-2s %2s %2s \t %4s %12.4f%12.4f\n" % (angle_1,angle_2,angle_3,'harm',angle_k*2.0,angle_r))
	line = f.readline().replace('-', ' ').split()
	if line in ['\n', '\r\n']:
		break
#
line=f.readline()
line=f.readline().replace('-', ' ').split()
OUTPUT.write('\n%s\n' % 'DIHEDRALS')
while len(line) > 0:
	ndihedral=ndihedral+1
	dihedral_1=(line[0])
	dihedral_2=(line[1])
	dihedral_3=(line[2])
	dihedral_4=(line[3])
	dihedral_k=float(line[5])
	dihedral_a=float(line[6])
	dihedral_p=float(line[7])
	OUTPUT.write("%-2s %2s %2s %2s \t %3s %12.4f %12.4f %12.4f %d\n" % (dihedral_1,dihedral_2,dihedral_3,dihedral_4,'cos',
	dihedral_k*2.0,dihedral_a,dihedral_p,0))
	line = f.readline().replace('-', ' ').split()
	if line in ['\n', '\r\n']:
		break
line=f.readline()
line=f.readline().replace('-', ' ').split()
OUTPUT.write('\n%s\n' % 'IMPROPERS')
while len(line) > 0:
	nimproper=nimproper+1
	dihedral_1=(line[0])
	dihedral_2=(line[1])
	dihedral_3=(line[2])
	dihedral_4=(line[3])
	dihedral_k=float(line[4])
	dihedral_a=float(line[5])
	dihedral_p=float(line[6])
	OUTPUT.write("%-2s %2s %2s %2s \t %3s %12.4f %12.4f %12.4f %d\n" % (dihedral_1,dihedral_2,dihedral_3,dihedral_4,'cos',
	dihedral_k*2.0,dihedral_a,dihedral_p,0))
	line = f.readline().replace('-', ' ').split()
	if line in ['\n', '\r\n']:
		break
#
Rtosig=math.pow(2.0,0.166666666)
line=f.readline()
line=f.readline().replace('-', ' ').split()
OUTPUT.write('\n%s\n' % 'VDW')
while len(line) > 0:
	nvdw=nvdw+1
	atom_1=(line[0])
	sigma=float(line[1])
	eps=float(line[2])
	OUTPUT.write("%-2s \t %2s%12.8f%12.8f\n" % (atom_1,'lj',sigma/Rtosig,eps))
	line=f.readline().replace('-', ' ').split()
	if line in ['\n', '\r\n']:
		break

print 'unique atoms:', nun
print 'bonds:', nbond
print 'angles:', nangle
print 'dihedrals:', ndihedral
print 'impropers:', nimproper
print 'VDW:', nvdw
OUTPUT.close()	
OUTPUT2.close()
f.close()


## Write out the output stats to be used by sort.py
OUTPUT=open('stats', 'w+')
OUTPUT.write('%d\t%d\t%d\t%d\t%d\t%d' % (nun, nbond, nangle, ndihedral, nimproper, nvdw))
OUTPUT.close()

#Convert the pdb to a format for DLPOLY
F=open("mimic.pdb")
OUTPUT=open("ATOM_LIST", "w+")
OUTPUT2=open("GAUSSIAN", "w+")
OUTPUT3=open("refconfiginput", "w+")

OUTPUT3.write("\n")
OUTPUT3.write("%s\n\n" % "&input")
OUTPUT3.write("%s %d\n" % ("Natom =  ", natoms))
OUTPUT3.write("%s\n" % "box_x= 500.00")
OUTPUT3.write("%s\n" % "box_y = 500.00")
OUTPUT3.write("%s\n\n" % "box_z = 500.00")
OUTPUT3.write("%s\n" % "/")

for line in F:
	liner=line.split()
	OUTPUT.write("%s\n" % liner[2])
	OUTPUT2.write("%s\t%s\t%s\t%s\n" % (liner[2],liner[4],liner[5],liner[6]))

OUTPUT3.close()
OUTPUT2.close()
OUTPUT.close()
F.close()

