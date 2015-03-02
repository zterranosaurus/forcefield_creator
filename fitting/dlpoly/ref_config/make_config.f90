program make_config

  integer, parameter :: Ndim=3, levcfg=0, imcon=1, izero=0
  real(8), parameter :: zero = 0.d0
  
  character(60) :: title
  character(10) :: atname
  integer :: iatom, Natom
  real(8) :: box_x, box_y, box_z, pos(Ndim)

  namelist/input/Natom, box_x, box_y, box_z

  open(unit=1, file='input', status="old", iostat=ierr)
  read(1,nml=input)
  close(1)

  open(unit=2,file='ATOM_LIST')
  open(unit=11,file='GAUSSIAN')
  open(unit=101,file='CONFIG')
  open(unit=102,file='CONFIG.xyz')




  !! CONFIG file.
  write(102,*) Natom
  write(102,*)
  title = 'mof fragment'
  write(101,'(a60)') title
  !write(101,'(3i10,f15.6)') izero, izero, Natom, zero
  write(101,'(3i10,f20.10)') levcfg, imcon, Natom, zero
  write(101,'(3f20.10)') 500.0, zero, zero
  write(101,'(3f20.10)') zero, 500.0, zero
  write(101,'(3f20.10)') zero, zero, 500.0

  do iatom = 1, Natom
     read(11,*) atname, pos(:)
     write(102,'(A2,5x,3f12.6)') atname, pos(:)
     read(2,*) atname
     write(101,'(a8,i10)') atname, iatom
     write(101,'(3f20.10)') pos(:)
  end do
  close(101)
  close(102)

end program make_config
