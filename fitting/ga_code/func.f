      subroutine func(ind,ssq,ssqq)
      implicit double precision (a-h,o-z)
c
c     Energies are minimized.
c
c     individual    index for the current parameter set 
c     ssq           sum of squares of deviations from ab initio data
c
      include 'ga_params.h'

      common /ga2/ nparam,nchrome
      common /ga3/ parent(nparmax,indmax),iparent(nchrmax,indmax)
      common /bestfit/ bestfit
      common /fitobj/ valabinitio(maxnv), weight(maxnv)
      common /fitcon/ nv

      dimension valdum(nv)

      character*150 :: call_buffer
      LOGICAL FIRST
      DATA FIRST /.TRUE./
      save first

      call prepfield(ind)
c
c     compute normal modes
      call system('../setup/computeMODES.sh')
c
      open(unit=10,file='pot_dlpoly.dat')
      do i=1,nv
         read(10,*) itmp, valdum(i)
      enddo
      close(10)

c     assign weighting factors and calc SSQ      
      ssq=0.d0
      do i=1,nv
c        if (valabinitio(i).eq.0.d0) then
            ssd = (valdum(i)-valabinitio(i))
c        else
c           ssd = (valdum(i)-valabinitio(i)) / valabinitio(i)
c        end if
         ssq = ssq + weight(i)*ssd*ssd
      enddo

c     Change to negative because GA maximizes instead of minimizes
      ssq=-ssq
c
c       write the value of the currently best fit
c
      if (first) then
         bestfit=abs(ssq)
         first=.false.
      else
         if (bestfit.gt.abs(ssq)) then
             bestfit=abs(ssq)
         else
            goto 999
         endif
      endif

      open(77,file='best_fit.txt',position='append')

      write(77,*)
      write(77,*)
      write(77,*)' calculated, target, difference'
      do i=1,nv
         write(77,'(3f15.3)')valdum(i),valabinitio(i),
     *                       valdum(i)-valabinitio(i)
      enddo
      write(77,*)
      write(77,*) 'SSQ =',ssq

      write(77,*)
      write(77,*)'   ========== Parameters =========='
      write(77,*)

      do i=1,nparam
         write(77,'(f20.10)') parent(i,ind)
      enddo
      close(77)

  999 continue
c
      return
      end

