program ord_pot

  implicit none

  integer :: iq, Nq, imode, Nmode, itmp
  real(8) :: pot0, q_min, dq, q
  real(8), allocatable :: pot(:)

  open(unit=1,file='ord.in')
  open(unit=101,file='ENERGIES')
  open(unit=102,file='pot_dlpoly.dat')

  read(1,*) Nmode
  read(1,*) Nq
  read(1,*) q_min
  read(1,*) dq
  allocate(pot(Nq))

  do imode = 1, Nmode
     pot0 = 1.d20
     do iq = 1, Nq
        read(101,*) pot(iq)
        if (pot(iq).le.pot0) pot0 = pot(iq)
     end do
     do iq = 1, Nq
        q = q_min + dble(iq-1)*dq
        write(102,*) q, pot(iq)-pot0
     end do
  end do
  close(101)
  close(102)

end
