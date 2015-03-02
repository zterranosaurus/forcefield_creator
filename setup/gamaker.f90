program makegainp
  use IFPORT
  implicit none

  integer :: i, j, k,nparm,bond,angle,dihed
  real :: x
  character(10):: int2char


  open(11,file='nparm')
  read(11,*) nparm,bond,angle,dihed
  close(11)

  open(19,file='ga.tmp')
  write(19,'(a)')' $ga            '
  write(19,'(a)')' icreep = 1,    '
  write(19,'(a)')' idum = -10000, '
  write(19,'(a)')' ielite = 1,    '
  write(19,'(a)')' iniche = 1,    '
  write(19,'(a)')' irestrt = 0,   '
  write(19,'(a)')' iskip = 0,     '
  write(19,'(a)')' itourny=1,     '
  write(19,'(a)')' iunifrm = 1,   '
  write(19,'(a)')' kountmx = 1,   '
  write(19,'(a)')' maxgen =  100, '
  write(19,'(a)')' microga = 1,   '
  write(19,'(a)')' nchild = 2,    '
  write(19,'(a)')' nichflg = 2,   '
  write(19,'(a)')' nowrite = 1,   '
  write(19,'(a10,I5,a1)') 'nparam = ',nparm,','
  write(19,'(a)')' npopsiz = 100,'
  write(19,'(a)')' nposibl ='//TRIM(int2char(nparm))//'*1048576,'
  do i=1,bond
     if(i.eq.1) then
        write(19,'(a)') 'parmin = 0.0,    0.0,    0.0,'
     else 
        write(19,'(a)') '         0.0,    0.0,    0.0,'
     endif
  enddo
  do i=1,angle
     write(19,'(a)') '         0.0,    0.0,'
  enddo
  do i=1,dihed
     write(19,'(a)') '         0.0,'
  enddo
  do i=1,bond
     if(i.eq.1) then
        write(19,'(a)') 'parmax = 500.0,    5.0,    5.0,'
     else 
        write(19,'(a)') '         500.0,    5.0,    5.0,'
     endif
  enddo
  do i=1,angle
     write(19,'(a)') '         300.0,    180.0,'
  enddo
  do i=1,dihed
     write(19,'(a)') '         8.0,'
  enddo
   write(19,'(a)')'pcreep = 0.04,'
   write(19,'(a)')'pcross = 0.50,'
  write(19,'(a)')' pmutate = 0.02,'
  write(19,'(a)')' $end'
  
  
  close(19)

end program makegainp
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
