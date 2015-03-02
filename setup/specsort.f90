program specifictxt
  use IFPORT
  implicit none

  integer :: i,num,ioerr,a,b,c,d,e,natoms,dumi,count(10),tmpa,j,start,ret
  integer :: totbonds,totangle,totdihed,totvdw,totimps,counter,nparm
  integer :: fitbond,fitangle,fitdihedral
  real :: q, m
  character(3)::c3
  character(2):: c2
  character(256) :: line(200),all
  character(256) :: line2,clip
  character(1024) :: SORTNMAKEUNIQUE

  SORTNMAKEUNIQUE ='grep ''natoms'' FFbuild.in | awk ''{print $2}'' > natoms.tmp'
  ret=system(SORTNMAKEUNIQUE)
  open(11,file='natoms.tmp')
  read(11,*) natoms
  close(11)
 ! natoms=52

  open(11,file='unstats')
  read(11,*) a,b,c,d,e
  close(11)

  open(19,file='nparm')
  read(19,*) nparm,fitbond,fitangle,fitdihedral
  close(19)

  b=fitbond
  c=fitangle
  d=fitdihedral


  !Write a dummy input for the buildff_fit
  open(19,file='input.DAT')
  write(19,*)
  write(19,'(A)')'&input'
  write(19,*)
  write(19,'(A7,I2)') 'Nuni = ',natoms
  write(19,'(A8,I2)') 'Natom = ',natoms
  write(19,'(A10)') 'Nbond =  0'
  write(19,'(A10)') 'Nangle = 0'
  write(19,'(A10)') 'Ndihe =  0'
  write(19,'(A,I2)') 'Nbond_fit = ',fitbond
  write(19,'(A,I2)') 'Nangle_fit = ',fitangle
  write(19,'(A,I2)') 'Ndihe_fit = ',fitdihedral
  write(19,'(A,I2)') 'Nimp_fit = ',0
  write(19,'(A7,1X,I2)') 'Nvdw = ',0
  write(19,'(A,I2)') 'Nfac = ',4
  write(19,'(A,I2)') 'wifacok = ',1
  write(19,*)
  write(19,'(a)') '/'

  close(19)


  open(20,file='new_params.dat')

  do i=1,fitbond
     write(20,'(4F20.10)') 0.0000, 0.0000, 0.0000, 0.0000
  enddo
  do i=1,fitangle
     write(20,'(4F20.10)') 0.0000, 0.0000, 0.0000, 0.0000
  enddo
  do i=1,fitdihedral
     write(20,'(4F20.10)') 0.0000, 180.0000, 2.0000, 2.0000
  enddo

  close(20)
!!$
!!$
  open(19,file='prepfield.f')
  open(21,file='check.f')
  write(21,'(a)')'       program prepfield                                '
  write(21,'(a)')'       implicit none                                    '
  write(21,'(a)')'       real(8), allocatable :: parent(:)                '
  write(21,'(a)')'       real :: Nq,pot1,pot2,zero,one,two,three,pi       '
  write(21,'(a)')'       integer :: c,i,j,k,Npot,nparm,i0,i1,i2,i3        '
  write(21,'(a)')'       character(20):: arg1,arg2,arg3,arg4              '
  write(21,'(a)')'                                                        '
  write(21,'(a)')'       call getarg(1,arg1)                              '
  write(21,'(a)')'       read(arg1,*) Npot                                '
  write(21,'(a)')'       call getarg(2,arg2)                              '
  write(21,'(a)')'       read(arg2,*) Nq                                  '
  write(21,'(a)')'                                                        '
  write(21,'(a)')'                                                        '
  write(21,'(a)')'       open(11,file=''nparm'')                            '
  write(21,'(a)')'       read(11,*) nparm                                 '
  write(21,'(a)')'                                                        '            
  write(21,'(a)')'       open(unit=2,file=''best_data'')                    '
  write(21,'(a)')'       open(unit=101,file=''pot.dat'')                    '
  write(21,'(a)')'       open(unit=102,file=''best_params.dat'')            '
  write(21,'(a)')'                                                        '                      
  write(21,'(a)')'                                                        '
  write(21,'(a)')'       allocate(parent(Nparm))                          '
  write(21,'(a)')'        c=1                                             '
  write(21,'(a)')'       do i = 1, Npot                                   '
  write(21,'(a)')'          do j = 1, Nq                                  '
  write(21,'(a)')'             read(2,*) pot1, pot2                       '
  write(21,'(a)')'             write(101,''(I,2f15.6)'') c, pot1, pot2      '
  write(21,'(a)')'                 c=c+1                                  '
  write(21,'(a)')'          end do                                        '
  write(21,'(a)')'          write(101,*)                                  '
  write(21,'(a)')'       end do                                           '
  write(21,'(a)')'       read(2,*)                                        '
  write(21,'(a)')'       read(2,*)                                        '
  write(21,'(a)')'       read(2,*)                                        '
  write(21,'(a)')'       read(2,*)                                        '
  write(21,'(a)')'       read(2,*)                                        '
  write(21,'(a)')'       i0 = 0                      '
  write(21,'(a)')'       i1 = 1                      '
  write(21,'(a)')'       i2 = 2                      '
  write(21,'(a)')'       i3 = 3                      '
  write(21,'(a)')'       zero = 0.d0                 '
  write(21,'(a)')'       one = 1.d0                  '
  write(21,'(a)')'       two = 2.d0                  '
  write(21,'(a)')'       three = 3.d0                '
  write(21,'(a)')'       pi = 180.d0                 '
  write(21,'(a)')'                                   '
  write(21,'(a)')'       do i = 1, Nparm             '
  write(21,'(a)')'          read(2,*) parent(i)      '
  write(21,'(a)')'       end do                      '

  write(19,'(a)')'      subroutine prepfield(ind)                                      '
  write(19,'(a)')'      implicit double precision (a-h,o-z)                            '
  write(19,'(a)')'                                                                     '
  write(19,'(a)')'      include ''ga_params.h''                                          '
  write(19,'(a)')'                                                                     '
  write(19,'(a)')'      common /evbpcon/ nv,ntmol,natms,numoen                         '
  write(19,'(a)')'      common /ecvpran/ oexmax, oexmin, pexmaxr, pexminr              '
  write(19,'(a)')'      common /ga3/ parent(nparmax,indmax),iparent(nchrmax,indmax)    '
  write(19,'(a)')'                                                                     '
  write(19,'(a)')'      open(80,file=''new_params.dat'')                                 '
  write(19,'(a)')'                                                                     '
  write(19,'(a)')'      i0 = 0                                                         '
  write(19,'(a)')'      i1 = 1                                                         '
  write(19,'(a)')'      i2 = 2                                                         '
  write(19,'(a)')'      i3 = 3                                                         '
  write(19,'(a)')'      zero = 0.d0                                                    '
  write(19,'(a)')'      one = 1.d0                                                     '
  write(19,'(a)')'      two = 2.d0                                                     '
  write(19,'(a)')'      three = 3.d0                                                   '
  write(19,'(a)')'      pi = 180.d0                                                    '
  write(19,'(a)')'                                                                     '
  write(19,'(a)')'c     Parameters.                                                    '

  j=1
  do i=1,b
     write(19,'(a36,I4,a1,I4,a11)')'      write(80,''(4f15.6)'') parent(',j,':',j+2,',ind), zero'
     write(21,'(a38,I4,a1,I4,a8)')'      write(102,''(4f15.6)'') parent(',j,':',j+2,TRIM('), zero')
     j=j+3
  enddo
  do i=b+1,b+c
     write(19,'(a36,I4,a1,I4,a17)')'      write(80,''(4f15.6)'') parent(',j,':',j+1,',ind), zero, zero'
     write(21,'(a38,I4,a1,I4,a13)')'      write(102,''(4f15.6)'') parent(',j,':',j+1,'), zero, zero'
     j=j+2
  enddo
  do i=b+c+1,b+c+d
     write(19,'(a41,I4,a22)')'      write(80,''(3f15.6,i10)'') parent(',j,',ind), 180.0,   two, 0'
     write(21,'(a42,I4,a22)')'      write(102,''(3f15.6,i10)'') parent(',j,'), 180.0,   two, 0'
     j=j+1
  enddo

  write(19,'(a)')'      close(80)'
  write(19,'(a)')
  write(19,'(a)')'      return'
  write(19,'(a)')'      end'
  close(19)
  write(21,'(a)')'      close(102)'
  write(21,'(a)')
  write(21,'(a)')'      end'
  close(121)

end program specifictxt




