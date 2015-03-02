#! /usr/bin/env python

AOI='XX'
natoms=XX
totalbond=XX
totalangle=XX
totaldihedral=XX
totalimproper=XX
totalvdw=XX


F=open('ffappend.tmp')
atom=0; bond=0; angle=0; dihedral=0; vdw=0
a=[]; b=[]; c=[]; d=[]; e=[]
for line in F:
	if AOI in line:
#		print line, len(line)
		if len(line) == 33:
			atom=atom+1
			a.append(line)
		if len(line) == 37:
			bond=bond+1
			b.append(line)
		if len(line) == 41:
			angle=angle+1
			c.append(line)
		if len(line) == 59:
			dihedral=dihedral+1
			d.append(line)
		if len(line) == 32:
			vdw=vdw+1
			e.append(line)

F.close()

#print atom, bond, angle, dihedral, vdw	
OUTPUT=open('unstats', 'w+')
OUTPUT.write('%d\t%d\t%d\t%d\t%d' % (atom, bond, angle, dihedral, vdw))
OUTPUT.close()

F=open('ffappend.tmp')
F2=open('force_field.tmp')

OUTPUT=open('force_field.dat','w+')
for i in range (0,natoms+3):
	line=F.readline()
	OUTPUT.write(line)
for i in range (0,natoms-atom+3):
	line2=F2.readline()
for i in range (0,bond):
	OUTPUT.write('%s\t %12.6f\t %12.6f\n' % (b[i].strip(), 0.0, 0.0))
for i in range (0,totalbond-bond):
	line2=F2.readline()
	OUTPUT.write('%s\t %12.6f\t %12.6f\n' % (line2.strip(), 0.0, 0.0))
line2=F2.readline()
OUTPUT.write(line2)
line2=F2.readline()
OUTPUT.write(line2)
for i in range (0,angle):
	OUTPUT.write('%s\t %12.6f\t %12.6f\n' % (c[i].strip(), 0.0, 0.0))
for i in range (0,totalangle-angle):
	line2=F2.readline()
	OUTPUT.write('%s\t %12.6f\t %12.6f\n' % (line2.strip(), 0.0, 0.0))
line2=F2.readline()
OUTPUT.write(line2)
line2=F2.readline()
OUTPUT.write(line2)
for i in range (0,dihedral):
	OUTPUT.write('%s\t\n' % (d[i].strip()))
for line in F2:
	OUTPUT.write(line)
for i in range (0,vdw):
	OUTPUT.write('%s' % e[i])


F.close()
F2.close()
	






