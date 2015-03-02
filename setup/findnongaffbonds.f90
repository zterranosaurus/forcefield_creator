program findthem
  implicit none

  integer :: i, j, k,coloc,c,aoipos,natom
  integer, parameter:: poss=100
  real :: x(100),y(100),z(100),rij,rx(100),ry(100),rz(100),q(poss),newx,newy,newz,tmpx,tmpy,tmpz
  character(10) :: arg1,arg2
  character(3) :: mol(100),dum3
  character(2) :: el(100),el2(100)

  open(11,file='mimic.mol2')
  open(12,file='tmp')
  open(13,file='mimic.pdb')
  open(14,file='tmpxyz')

  do i=1,8
     read(11,*)
     read(12,*)
  enddo

  call getarg(1,arg1)
  read(arg1,*) aoipos
  call getarg(2,arg2)
  read(arg2,*) natom

  do i=1,natom
     read(12,*) el2(i),x(i),y(i),z(i)
  enddo
  close(12)
  
  do i=1,natom
     newx=x(i)-x(aoipos); newy=y(i)-y(aoipos); newz=z(i)-z(aoipos)
     write(14,*) el2(i),newx,newy,newz
  enddo
  do i=1,natom
     read(11,*) k,mol(i),x(i),y(i),z(i),el(i),k,dum3,q(i)
     write(13,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,10X,2A2)') 'HETATM',i,el(i),' ','   ',' ',i,'',x(i),y(i),z(i)
  enddo
  close(13)
  close(11)


!!$  coloc=1
!!$  do i=1,natom
!!$     if(i.ne.coloc) then
!!$        call dis(x(i),y(i),z(i),x(coloc),y(coloc),z(coloc),rij)
!!$        if(rij.lt.2.4) print *,1,coloc,i,1
!!$     endif
!!$  enddo
!!$
!!$  coloc=3
!!$  do i=1,52
!!$     if(i.ne.coloc) then
!!$        call dis(x(i),y(i),z(i),x(coloc),y(coloc),z(coloc),rij)
!!$        if(rij.lt.2.4) print *,1,coloc,i,1
!!$     endif
!!$  enddo
!!$
!!$  coloc=14
!!$  do i=1,52
!!$     if(i.ne.coloc) then
!!$        call dis(x(i),y(i),z(i),x(coloc),y(coloc),z(coloc),rij)
!!$        if(rij.lt.2.4) print *,1,coloc,i,1
!!$     endif
!!$  enddo




end program findthem
subroutine dis(xs,ys,zs,ref1,ref2,ref3,r_ij)
  implicit none
  integer i, j, k
  real:: xs,ys,zs,ref1,ref2,ref3,r_ij
  real :: delx,dely,delz

  delx=ref1-xs
  dely=ref2-ys
  delz=ref3-zs

  r_ij=sqrt(delx*delx+dely*dely+delz*delz)

end subroutine dis
  
