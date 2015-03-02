program compute_pot

  implicit none

  real(8), parameter :: energ_au2kcal = 627.509d0
  real(8), parameter :: qscale = 0.01d0

  integer :: iq, Nq
  
  real(8) :: enrgeq, weight
  real(8), allocatable :: enrg(:), q(:)

  open(unit=1,file='parm.dat')
  open(unit=2,file='enrg.dat')
  open(unit=101,file='pot.dat')

  read(1,*) Nq
  read(1,*) weight
  allocate(q(Nq), enrg(Nq))

  enrgeq = 0.d0
  do iq = 1, Nq
     read(2,*) q(iq), enrg(iq)
     if (enrg(iq).le.enrgeq) enrgeq = enrg(iq)
  end do

  do iq = 1, Nq
     write(101,'(2f15.6,f10.4)') q(iq)*qscale, (enrg(iq)-enrgeq)*energ_au2kcal, weight
  end do

end 

