c        Calls David Carroll's genetic algorithm program
c        to optimize PES parameters

      program fit

      implicit double precision (a-h,o-z)

      include 'ga_params.h'

      common /fitobj/ valabinitio(maxnv),weight(maxnv)
      common /fitcon/ nv

c
c     read in ab initio data
      open(11,file='target_pot.dat')
      read(11,*) Nmode, Nq
      nv = Nmode*Nq
      i = 0
      do imode=1,Nmode
         do iq = 1, Nq
            i = i + 1
            read(11,*) tmp, valabinitio(i), weight(i)
         end do
      end do
      close(11)

      open(24,file='ga.out')
      rewind(24)

      call ga_carroll

      close(24)

      stop
      end

