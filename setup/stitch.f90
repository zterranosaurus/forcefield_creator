program combineandspecandbest
implicit none

integer :: i, j, k,bonds,angles,dihedrals,nparm
real :: c1,c2,c3,c4
character(2) :: el1,el2,el3,el4
character(4) :: harm
character(3) :: mors

open(11,file='nparm')
read(11,*) nparm,bonds,angles,dihedrals
close(11)

open(11,file='best_params.dat')
open(12,file='force_field.dat')
open(13,file='BESTPARAMS')
	do i=1,21
	read(12,*)
	enddo
	do i=1,bonds
	read(11,*) c1,c2,c3,c4
	read(12,*) el1,el2,harm
	write(13,'(a2,5x,a2,5x,a4,5x,4f12.5)') el1,el2,'mors',c1,c2,c3,c4
	enddo
	read(12,*)
	read(12,*)
	do i=1,angles
	  read(11,*) c1,c2,c3,c4
        read(12,*) el1,el2,el3,harm
        write(13,'(a2,5x,a2,5x,a2,5x,a4,5x,4f12.5)') el1,el2,el3,'harm',c1,c2,c3,c4
        enddo
	read(12,*)
	read(12,*)
	do i=1,dihedrals
        read(11,*) c1,c2,c3,c4
        read(12,*) el1,el2,el3,el4,harm
        write(13,'(a2,5x,a2,5x,a2,5x,a2,5x,a3,5x,4f12.5)') el1,el2,el3,el4,'cos',c1,c2,c3,c4
        enddo
close(11)
close(12)
close(13)

end program combineandspecandbest
