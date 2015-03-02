program compute_pot

  implicit none

  real(8), parameter :: energ_au2kcal = 627.509d0
  real(8), parameter :: qscale = 0.01d0

  integer :: iq, Nq, iscan, Nscan, N

  real(8) :: enrgeq, weight,icount,tmp,dq
  real(8), allocatable :: enrg(:), q(:)

  open(unit=1,file='parm.dat')
  open(unit=2,file='enrg.dat')
  open(unit=101,file='pot.dat')

  read(1,*) N 
  read(1,*) Nq
  read(1,*) dq
  allocate(q(N), enrg(N))

  Nscan = N / Nq
  print*, 'Nscan = ', Nscan

  do iscan = 1, Nscan
     icount = -((Nq-1)/2)*dq	
     enrgeq = 0.d0
     do iq = 1, Nq
        read(2,*) q(iq), enrg(iq)
        if (enrg(iq).le.enrgeq) enrgeq = enrg(iq)
     end do
     do iq = 1, Nq
        tmp=(enrg(iq)-enrgeq)*energ_au2kcal
        if(tmp.gt.100) weight = 1
        if(tmp.gt.10.and.tmp.lt.100) weight = 10 
        if(tmp.gt.1.and.tmp.lt.10) weight = 100
        if(tmp.lt.1) weight = 1000
        write(101,'(f10.4,f15.6,f10.4)') icount, (enrg(iq)-enrgeq)*energ_au2kcal, weight
        icount = icount + dq
     end do
     icount=0
     
  end do


end program compute_pot

