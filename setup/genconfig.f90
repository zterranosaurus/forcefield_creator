program genconfig
  use IFPORT
  implicit none

  integer:: i,Npoints,filesfin,set,natoms,ret,tmp,c,j,k,icount
  integer,parameter :: poss=100
  real :: x(poss),y(poss),z(poss),stepsize,xtmp,ytmp,ztmp,scale
  character(2) :: atom(poss)
  character(22):: infile
  CHARACTER(20) :: arg1, arg2, arg3,arg5,arg6,arg7
  character(10):: int2char
  CHARACTER(100) :: command, line, atommoved(5000),posline(5000)

  call getarg(1,arg1)

  if(arg1.ne.'gaussian') then
     print *, './configure.o gaussian Npoints, stepsize, xyzfile, set#, natoms atomofinterestpos'
     STOP
  endif
  call getarg(2,arg2)
  read(arg2,*) Npoints
  call getarg(3,arg3)
  read(arg3,*) stepsize
  call getarg(4,infile)
  open(10, file=infile, action="read")
  call getarg(5,arg5)
  read(arg5,*) set
  call getarg(6,arg6)
  read(arg6,*) natoms
  call getarg(7,arg7)
  read(arg7,*) c
  

  do i=1,natoms
     read(10,*) atom(i),x(i),y(i),z(i)
  enddo
  close(10)

  scale=(Npoints-1)/2*stepsize
  print *, 'max deviation from atom''s original place =', scale


  open(12,file='atommoved.tmp')
  open(13,file='posmoved.tmp')
  icount=0
  xtmp=x(c)-scale
  do k=1,Npoints
     ytmp=y(c)-scale
     do j=1,Npoints
        ztmp=z(c)-scale
        do i=1,Npoints
           write(12,'(A2,3F12.6)') atom(c),xtmp,ytmp,ztmp
           write(13,'(3e20.9,a8)') xtmp,ytmp,ztmp,atom(c)
           ztmp=ztmp+stepsize
           icount=icount+1
        enddo
        ytmp=ytmp+stepsize
     enddo
     xtmp=xtmp+stepsize
  enddo

  print *, 'Total configurations: ',icount
  close(12)
  close(13)

  open(13,file='atommoved.tmp')
  open(14,file='posmoved.tmp')
  do i=1,icount
     read(13,'(A)') atommoved(i)
     read(14,'(A)') posline(i)
  enddo
  close(13)
  close(14)

  open(31,file='POSITIONS')
  open(33,file='PESscan.xyz')

  do j=1,icount
  write(31,'(3e20.9)') 500.0, 0.0, 0.0
  write(31,'(3e20.9)') 0.0, 500.0, 0.0
  write(31,'(3e20.9)') 0.0, 0.0, 500.0
     call openfile(23,j,'config','.com')
     call openfile(25,j,'config','.xyz')
     write(33,*) 52
     write(33,*) 
     do i=1,c-1
        write(23,'(A2,3F12.6)') atom(i),x(i),y(i),z(i)
        write(25,'(A2,3F12.6)') atom(i),x(i),y(i),z(i)
        write(31,'(3e20.9,a8)') x(i),y(i),z(i),atom(i)
        write(33,'(A2,3F12.6)') atom(i),x(i),y(i),z(i)
     enddo
     write(23,'(A)') TRIM(atommoved(j))
     write(25,'(A)') TRIM(atommoved(j))
     write(31,'(A)') TRIM(posline(j))
     write(33,'(A)') TRIM(atommoved(j))
     do i=c+1,natoms
        write(23,'(A2,3F12.6)') atom(i),x(i),y(i),z(i)
        write(25,'(A2,3F12.6)') atom(i),x(i),y(i),z(i)
        write(31,'(3e20.9,a8)') x(i),y(i),z(i),atom(i)
        write(33,'(A2,3F12.6)') atom(i),x(i),y(i),z(i)
     enddo
     write(23,'(A)') '                                                          '
     close(23)
     close(25)
  enddo
  close(31)
  close(33)

  command='mkdir -p set'//arg5
  ret=system(command)

  command='mv *.com *.xyz set'//arg5
  print *, TRIM(command)
  ret=system(command)

  
  do i=1,Npoints*Npoints*Npoints,100
     call openfile(23,i,'jobPES','.pbs')
     write(23,*)
     write(23,*)
     do j=i,i+99
        if(j.le.Npoints*Npoints*Npoints) then
           write(23,'(a)') 'g09l < config'//TRIM(int2char(j))//'.com > config'//TRIM(int2char(j))//'.out'
        endif
     enddo
     write(23,*)
     close(23)
  enddo


  close(10)

end program genconfig

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Converts an integer to a character... (written by JR Schmidt)

character*10 function int2char(n)
  implicit none     
  integer n,pow
  int2char = ""
  do pow = 0, 9
     if (10**pow .le. n) then
        int2char = char(48 + mod(n/10**pow, 10))//int2char
     endif
  enddo
end function int2char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Opens a file
subroutine openfile(iunit, n, pre, suf)
  character*10 int2char, cfile
  character*1 cfile1
  character*2 cfile2
  character*3 cfile3
  character*4 cfile4
  character*6 pre
  character*4 suf      
  cfile=int2char(n)
  if (n.le.9) then
     cfile1=cfile
     open(iunit, file=pre//cfile1//suf)
  endif
  if (n.ge.10.and.n.le.99) then
     cfile2=cfile
     open(iunit, file=pre//cfile2//suf)
  endif
  if (n.ge.100.and.n.le.999) then
     cfile3=cfile
     open(iunit, file=pre//cfile3//suf)
  endif
  if (n.ge.1000.and.n.le.9999) then
     cfile4=cfile
     open(iunit, file=pre//cfile4//suf)
  endif
end subroutine openfile



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

