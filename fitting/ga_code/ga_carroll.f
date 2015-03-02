c#######################################################################
c
      subroutine ga_carroll
c
c  This is version 1.7a, last updated on 4/2/2001.
c  Last minor bug found 6/14/00.
c
c  Copyright David L. Carroll; this code may not be reproduced for sale
c  or for use in part of another code for sale without the express 
c  written permission of David L. Carroll.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  David L. Carroll
c  CU Aerospace
c  2004 South Wright Street Extended
c  Urbana, IL  61802
c
c  e-mail: carroll@cuaerospace.com
c  Phone:  217-333-8274
c  fax:    217-244-7757
c
c  This genetic algorithm (GA) driver is free for public use.  My only
c  request is that the user reference and/or acknowledge the use of this
c  driver in any papers/reports/articles which have results obtained
c  from the use of this driver.  I would also appreciate a copy of such
c  papers/articles/reports, or at least an e-mail message with the
c  reference so I can get a copy.  Thanks.
c
c  This program is a FORTRAN version of a genetic algorithm driver.
c  This code initializes a random sample of individuals with different
c  parameters to be optimized using the genetic algorithm approach, i.e.
c  evolution via survival of the fittest.  The selection scheme used is
c  tournament selection with a shuffling technique for choosing random
c  pairs for mating.  The routine includes binary coding for the
c  individuals, jump mutation, creep mutation, and the option for
c  single-point or uniform crossover.  Niching (sharing) and an option
c  for the number of children per pair of parents has been added.
c  An option to use a micro-GA is also included.
c
c  For companies wishing to link this GA driver with an existing code,
c  I am available for some consulting work.  Regardless, I suggest
c  altering this code as little as possible to make future updates
c  easier to incorporate.
c
c  Any users new to the GA world are encouraged to read David Goldberg's
c  "Genetic Algorithms in Search, Optimization and Machine Learning,"
c  Addison-Wesley, 1989.
c
c  Other associated files are:  ga.inp
c                               ga.out
c                               ga.restart
c                               ga_params.h
c                               ReadMe
c                               ga2.inp (w/ different namelist identifier)
c
c  I have provided a sample subroutine "func", but ultimately
c  the user must supply this subroutine "func" which should be your
c  cost function.  You should be able to run the code with the
c  sample subroutine "func" and the provided ga.inp file and obtain
c  the optimal function value of 1.0000 at generation 187 with the
c  uniform crossover micro-GA enabled (this is 935 function evaluations).
c
c  The code is presently set for a maximum population size of 200,
c  30 chromosomes (binary bits) and 8 parameters.  These values can be
c  changed in ga_params.h as appropriate for your problem.  Correspondingly
c  you will have to change a few 'write' and 'format' statements if you
c  change nchrome and/or nparam.  In particular, if you change nchrome
c  and/or nparam, then you should change the 'format' statement numbers
c  1050, 1075, 1275, and 1500 (see ReadMe file).
c
c  Please feel free to contact me with questions, comments, or errors
c  (hopefully none of latter).
c
c  Disclaimer:  this program is not guaranteed to be free of error
c  (although it is believed to be free of error), therefore it should
c  not be relied on for solving problems where an error could result in
c  injury or loss.  If this code is used for such solutions, it is
c  entirely at the user's risk and the author disclaims all liability.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit real*8 (a-h,o-z)
      save
c
      include 'ga_params.h'
      dimension parent(nparmax,indmax),child(nparmax,indmax)
      dimension fitness(indmax),nposibl(nparmax),nichflg(nparmax)
      dimension iparent(nchrmax,indmax),ichild(nchrmax,indmax)
      dimension g0(nparmax),g1(nparmax),ig2(nparmax)
      dimension ibest(nchrmax)
      dimension parmax(nparmax),parmin(nparmax),pardel(nparmax)
      dimension geni(1000000),genavg(1000000),genmax(1000000)
CX    
      dimension ssqq(1000000)
CX
c      real*4 cpu,cpu0,cpu1,tarray(2)
c
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga4   / fitness
      common / ga5   / g0,g1,ig2
      common / ga6   / parmax,parmin,pardel,nposibl
      common / ga7   / child,ichild
      common / ga8   / nichflg
      common /inputga/ pcross,pmutate,pcreep,maxgen,idum,irestrt,
     +                 itourny,ielite,icreep,iunifrm,iniche,
     +                 iskip,iend,nchild,microga,kountmx
CX
      COMMON / SEED   / NSEED
CX
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Input variable definitions:
c
c  icreep   = 0 for no creep mutations
c           = 1 for creep mutations; creep mutations are recommended.
c  idum     The initial random number seed for the GA run.  Must equal
c           a negative integer, e.g. idum=-1000.
c  ielite   = 0 for no elitism (best individual not necessarily
c               replicated from one generation to the next).
c           = 1 for elitism to be invoked (best individual replicated
c               into next generation); elitism is recommended.
c  iend         = 0 for normal GA run (this is standard).
c           = number of last population member to be looked at in a set
c             of individuals.  Setting iend-0 is only used for debugging
c             purposes and is commonly used in conjunction with iskip.
c  iniche   = 0 for no niching
c           = 1 for niching; niching is recommended.
c  irestrt  = 0 for a new GA run, or for a single function evaluation
c           = 1 for a restart continuation of a GA run.
c  iskip    = 0 for normal GA run (this is standard).
c           = number in population to look at a specific individual or
c             set of individuals.  Setting iskip-0 is only used for
c             debugging purposes.
c  itourny  No longer used.  The GA is presently set up for only
c           tournament selection.
c  iunifrm  = 0 for single-point crossover
c           = 1 for uniform crossover; uniform crossover is recommended.
c  kountmx  = the maximum value of kount before a new restart file is
c             written; presently set to write every fifth generation.
c             Increasing this value will reduce I/O time requirements
c             and reduce wear and tear on your storage device
c  maxgen   The maximum number of generations to run by the GA.
c           For a single function evaluation, set equal to 1.
c  microga  = 0 for normal conventional GA operation
c           = 1 for micro-GA operation (this will automatically reset
c             some of the other input flags).  I recommend using
c             npopsiz=5 when microga=1.
c  nchild   = 1 for one child per pair of parents (this is what I
c               typically use).
c           = 2 for two children per pair of parents (2 is more common
c               in GA work).
c  nichflg  = array of 1/0 flags for whether or not niching occurs on
c             a particular parameter.  Set to 0 for no niching on
c             a parameter, set to 1 for niching to operate on parameter.
c             The default value is 1, but the implementation of niching
c             is still controlled by the flag iniche.
c  nowrite  = 0 to write detailed mutation and parameter adjustments
c           = 1 to not write detailed mutation and parameter adjustments
c  nparam   Number of parameters (groups of bits) of each individual.
c           Make sure that nparam matches the number of values in the
c           parmin, parmax and nposibl input arrays.
c  npopsiz  The population size of a GA run (typically 100 works well).
c           For a single calculation, set equal to 1.
c  nposibl  = array of integer number of possibilities per parameter.
c             For optimal code efficiency set nposibl=2**n, i.e. 2, 4,
c             8, 16, 32, 64, etc.
c  parmax   = array of the maximum allowed values of the parameters
c  parmin   = array of the minimum allowed values of the parameters
c  pcreep   The creep mutation probability.  Typically set this
c           = (nchrome/nparam)/npopsiz.
c  pcross   The crossover probability.  For single-point crossover, a
c           value of 0.6 or 0.7 is recommended.  For uniform crossover,
c           a value of 0.5 is suggested.
c  pmutate  The jump mutation probability.  Typically set = 1/npopsiz.
c
c
c  For single function evaluations, set npopsiz=1, maxgen=1, & irestrt=0.
c
c  My favorite initial choices of GA parameters are:
c     microga=1, npopsiz=5, iunifrm=1, maxgen=200
c     microga=1, npopsiz=5, iunifrm=0, maxgen=200
c  I generally get good performance with both the uniform and single-
c  point crossover micro-GA.
c
c  For those wishing to use the more conventional GA techniques,
c  my old favorite choice of GA parameters was:
c     iunifrm=1, iniche=1, ielite=1, itourny=1, nchild=1
c  For most problems I have dealt with, I get good performance using
c     npopsiz=100, pcross=0.5, pmutate=0.01, pcreep=0.02, maxgen=26
c  or
c     npopsiz= 50, pcross=0.5, pmutate=0.02, pcreep=0.04, maxgen=51
c
c  Any negative integer for idum should work.  I typically arbitrarily
c  choose idum=-10000 or -20000.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Code variable definitions (those not defined above):
c
c  best     = the best fitness of the generation
c  child    = the floating point parameter array of the children
c  cpu      = cpu time of the calculation
c  cpu0,cpu1= cpu times associated with 'etime' timing function
c  creep    = +1 or -1, indicates which direction parameter creeps
c  delta    = del/nparam
c  diffrac  = fraction of total number of bits which are different
c             between the best and the rest of the micro-GA population.
c             Population convergence arbitrarily set as diffrac<0.05.
c  evals    = number of function evaluations
c  fbar     = average fitness of population
c  fitness  = array of fitnesses of the parents
c  fitsum   = sum of the fitnesses of the parents
c  genavg   = array of average fitness values for each generation
c  geni     = generation array
c  genmax   = array of maximum fitness values for each generation
c  g0       = lower bound values of the parameter array to be optimized.
c             The number of parameters in the array should match the
c             dimension set in the above parameter statement.
c  g1       = the increment by which the parameter array is increased
c             from the lower bound values in the g0 array.  The minimum
c             parameter value is g0 and the maximum parameter value
c             equals g0+g1*(2**g2-1), i.e. g1 is the incremental value
c             between min and max.
c  ig2      = array of the number of bits per parameter, i.e. the number
c             of possible values per parameter.  For example, ig2=2 is
c             equivalent to 4 (=2**2) possibilities, ig2=4 is equivalent
c             to 16 (=2**4) possibilities.
c  ig2sum   = sum of the number of possibilities of ig2 array
c  ibest    = binary array of chromosomes of the best individual
c  ichild   = binary array of chromosomes of the children
c  icount   = counter of number of different bits between best
c             individual and other members of micro-GA population
c  icross   = the crossover point in single-point crossover
c  indmax   = maximum # of individuals allowed, i.e. max population size
c  iparent  = binary array of chromosomes of the parents
c  istart   = the generation to be started from
c  jbest    = the member in the population with the best fitness
c  jelite   = a counter which tracks the number of bits of an individual
c             which match those of the best individual
c  jend     = used in conjunction with iend for debugging
c  jstart   = used in conjunction with iskip for debugging
c  kount    = a counter which controls how frequently the restart
c             file is written
c  kelite   = kelite set to unity when jelite=nchrome, indicates that
c             the best parent was replicated amongst the children
c  mate1    = the number of the population member chosen as mate1
c  mate2    = the number of the population member chosen as mate2
c  nchrmax  = maximum # of chromosomes (binary bits) per individual
c  nchrome  = number of chromosomes (binary bits) of each individual
c  ncreep   = # of creep mutations which occurred during reproduction
c  nmutate  = # of jump mutations which occurred during reproduction
c  nparmax  = maximum # of parameters which the chromosomes make up
c  paramav  = the average of each parameter in the population
c  paramsm  = the sum of each parameter in the population
c  parent   = the floating point parameter array of the parents
c  pardel   = array of the difference between parmax and parmin
c  rand     = the value of the current random number
c  npossum  = sum of the number of possible values of all parameters
c  tarray   = time array used with 'etime' timing function
c  time0    = clock time at start of run
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Subroutines:
c  ____________
c
c  code     = Codes floating point value to binary string.
c  crosovr  = Performs crossover (single-point or uniform).
c  decode   = Decodes binary string to floating point value.
c  evalout  = Evaluates the fitness of each individual and outputs
c             generational information to the 'ga.out' file.
c  func     = The function which is being evaluated.
c  gamicro  = Implements the micro-GA technique.
c  input    = Inputs information from the 'ga.inp' file.
c  initial  = Program initialization and inputs information from the
c             'ga.restart' file.
c  mutate   = Performs mutation (jump and/or creep).
c  newgen   = Writes child array back into parent array for new
c             generation; also checks to see if best individual was
c             replicated (elitism).
c  niche    = Performs niching (sharing) on population.
c  possibl  = Checks to see if decoded binary string falls within
c             specified range of parmin and parmax.
c  ran3     = The random number generator.
c  restart  = Writes the 'ga.restart' file.
c  select   = A subroutine of 'selectn'.
c  selectn  = Performs selection; tournament selection is the only
c             option in this version of the code.
c  shuffle  = Shuffles the population randomly for selection.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      call etime(tarray)
c      write(6,*) tarray(1),tarray(2)
c      cpu0=tarray(1)
c
CX
      DO NSEED=28408,28408
CX
c  Call the input subroutine.
c      TIME0=SECNDS(0.0)
      call input
c
c  Perform necessary initialization and read the ga.restart file.
      call initial(istart,npossum,ig2sum)
c
c  $$$$$ Main generational processing loop. $$$$$
      kount=0
      do 20 i=istart,maxgen+istart-1
         write (6,1111) i
         write (24,1111) i
         write(24,1050)
c
c  Evaluate the population, assign fitness, establish the best
c  individual, and write output information.
         call evalout(iskip,iend,ibest,fbar,best,bestssqq)
         geni(i)=float(i)
         genavg(i)=fbar
         genmax(i)=best
CX
         ssqq(i)=bestssqq
CX
         if(npopsiz.eq.1 .or. iskip.ne.0) then
            close(24)
            stop
         endif
c
c  Implement "niching".
         if (iniche.ne.0) call niche
c
c  Enter selection, crossover and mutation loop.
         ncross=0
         ipick=npopsiz
         do 45 j=1,npopsiz,nchild
c
c  Perform selection.
            call selectn(ipick,j,mate1,mate2)
c
c  Now perform crossover between the randomly selected pair.
            call crosovr(ncross,j,mate1,mate2)
 45      continue
         write(6,1225) ncross
         write(24,1225) ncross
c
c  Now perform random mutations.  If running micro-GA, skip mutation.
         if (microga.eq.0) call mutate
c
c  Write child array back into parent array for new generation.  Check
c  to see if the best parent was replicated.
         call newgen(ielite,npossum,ig2sum,ibest)
c
c  Implement micro-GA if enabled.
         if (microga.ne.0) call gamicro(i,npossum,ig2sum,ibest)
c
c  Write to restart file.
         call restart(i,istart,kount)
 20   continue
c  $$$$$ End of main generational processing loop. $$$$$
c 999  continue
      write(24,3000)
      do 100 i=istart,maxgen+istart-1
         evals=float(npopsiz)*geni(i)
         write(24,3100) geni(i),evals,genavg(i),genmax(i)
         write(24,*) NSEED,genmax(maxgen),ssqq(maxgen)
 100  continue
c      call etime(tarray)
c      write(6,*) tarray(1),tarray(2)
c      cpu1=tarray(1)
c      cpu=(cpu1-cpu0)
c      write(6,1400) cpu,cpu/60.0
c      write(24,1400) cpu,cpu/60.0
c      CLOSE (24)
c
 1050 format(1x,' #      Binary Code',16x,'Param1  Param2  Fitness')
 1111 format(//'#################  Generation',i5,'  #################')
 1225 format(/'  Number of Crossovers      =',i5)
c 1400 format(2x,'CPU time for all generations=',e12.6,' sec'/
c     +       2x,'                             ',e12.6,' min')
 3000 format(2x//'Summary of Output'/
     +       2x,'Generation   Evaluations   Avg.Fitness   Best Fitness')
 3100 format(2x,3(e10.4,4x),e11.5)
c
CX 
      END DO
CX
      stop
      end
c
c#######################################################################
      subroutine input
c
c  This subroutine inputs information from the ga.inp (gafort.in) file.
c
      implicit real*8 (a-h,o-z)
      save
c
      include 'ga_params.h'
      dimension nposibl(nparmax),nichflg(nparmax)
      dimension parmax(nparmax),parmin(nparmax),pardel(nparmax)
c
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga6   / parmax,parmin,pardel,nposibl
      common / ga8   / nichflg
      common /inputga/ pcross,pmutate,pcreep,maxgen,idum,irestrt,
     +                 itourny,ielite,icreep,iunifrm,iniche,
     +                 iskip,iend,nchild,microga,kountmx
c
      namelist / ga   / irestrt,npopsiz,pmutate,maxgen,idum,pcross,
     +                  itourny,ielite,icreep,pcreep,iunifrm,iniche,
     +                  iskip,iend,nchild,nparam,parmin,parmax,nposibl,
     +                  nowrite,nichflg,microga,kountmx
ccd   namelist / ga   / irestrt,npopsiz,pmutate,maxgen,idum,pcross,
ccd  +                  itourny,ielite,icreep,pcreep,iunifrm,iniche,
ccd  +                  iskip,iend,nchild,nparam,              nposibl,
ccd  +                  nowrite,nichflg,microga,kountmx
c
      kountmx=5
      irestrt=0
      itourny=0
      ielite=0
      iunifrm=0
      iniche=0
      iskip=0
      iend=0
      nchild=1
      do 2 i=1,nparmax
         nichflg(i)=1
 2    continue
      microga=0
c
      OPEN (UNIT=23, FILE='ga.inp', STATUS='OLD')
      READ (23, NML = ga)
      CLOSE (23)
      itourny=1
c      if (itourny.eq.0) nchild=2
c
c  Check for array sizing errors.
      if (npopsiz.gt.indmax) then
         write(6,1600) npopsiz
         write(24,1600) npopsiz
         close(24)
         stop
      endif
      if (nparam.gt.nparmax) then
         write(6,1700) nparam
         write(24,1700) nparam
         close(24)
         stop
      endif
c
c  If using the microga option, reset some input variables
      if (microga.ne.0) then
         pmutate=0.0d0
         pcreep=0.0d0
         itourny=1
         ielite=1
         iniche=0
         nchild=1
         if (iunifrm.eq.0) then
            pcross=1.0d0
         else
            pcross=0.5d0
         endif
      endif
c
 1600 format(1x,'ERROR: npopsiz > indmax.  Set indmax = ',i6)
 1700 format(1x,'ERROR: nparam > nparmax.  Set nparmax = ',i6)
c
      return
      end
c
c#######################################################################
      subroutine initial(istart,npossum,ig2sum)
c
c  This subroutine sets up the program by generating the g0, g1 and
c  ig2 arrays, and counting the number of chromosomes required for the
c  specified input.  The subroutine also initializes the random number
c  generator, parent and iparent arrays (reads the ga.restart file).
      implicit real*8 (a-h,o-z)
      save
c
      include 'ga_params.h'
      dimension parent(nparmax,indmax),iparent(nchrmax,indmax)
      dimension nposibl(nparmax)
      dimension g0(nparmax),g1(nparmax),ig2(nparmax)
      dimension parmax(nparmax),parmin(nparmax),pardel(nparmax)
c
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga5   / g0,g1,ig2
      common / ga6   / parmax,parmin,pardel,nposibl
      common /inputga/ pcross,pmutate,pcreep,maxgen,idum,irestrt,
     +                 itourny,ielite,icreep,iunifrm,iniche,
     +                 iskip,iend,nchild,microga,kountmx
CX
      COMMON / SEED   / NSEED
CX
c
c
      open(24,file='ga.out',access='append')      ! ccd ---------
      do 3 i=1,nparam
         g0(i)=parmin(i)
         pardel(i)=parmax(i)-parmin(i)
         g1(i)=pardel(i)/dble(nposibl(i)-1)
 3    continue
      do 6 i=1,nparam
         do 7 j=1,30
            n2j=2**j
            if (n2j.ge.nposibl(i)) then
               ig2(i)=j
               goto 8
            endif
            if (j.ge.30) then
               write(6,2000)
               write(24,2000)
               close(24)
               stop
            endif
 7       continue
 8       continue
 6    continue
c
c  Count the total number of chromosomes (bits) required
      nchrome=0
      npossum=0
      ig2sum=0
      do 9 i=1,nparam
         nchrome=nchrome+ig2(i)
         npossum=npossum+nposibl(i)
         ig2sum=ig2sum+(2**ig2(i))
 9    continue
      if (nchrome.gt.nchrmax) then
         write(6,1800) nchrome
         write(24,1800) nchrome
         close(24)
         stop
      endif
c
      if (npossum.lt.ig2sum .and. microga.ne.0) then
         write(6,2100)
         write(24,2100)
      endif
c
c  Initialize random number generator
CX
      KSEED=-NSEED
      call ran3(KSEED,rand)
CX
c
      if(irestrt.eq.0) then
c  Initialize the random distribution of parameters in the individual
c  parents when irestrt=0.
         istart=1
         OPEN (UNIT=25, FILE='ga.restart',status='new')
         close(25)
         do 10 i=1,npopsiz
            do 15 j=1,nchrome
               call ran3(1,rand)
               iparent(j,i)=1
               if(rand.lt.0.5d0) iparent(j,i)=0
 15         continue
 10      continue
         if (npossum.lt.ig2sum) call possibl(parent,iparent)
      else
c  If irestrt.ne.0, read from restart file.
         OPEN (UNIT=25, FILE='ga.restart', STATUS='OLD')
         rewind 25
         read(25,*) istart,npopsiz
         do 1 j=1,npopsiz
            read(25,*) k,(iparent(l,j),l=1,nchrome)
 1       continue
         CLOSE (25)
      endif
c
      if(irestrt.ne.0) call ran3(idum-istart,rand)
c
 1800 format(1x,'ERROR: nchrome > nchrmax.  Set nchrmax = ',i6)
 2000 format(1x,'ERROR: You have a parameter with a number of '/
     +       1x,'   possibilities > 2**30!  If you really desire this,'/
     +       1x,'   change the DO loop 7 statement and recompile.'//
     +       1x,'   You may also need to alter the code to work with'/
     +       1x,'   REAL numbers rather than INTEGER numbers; Fortran'/
     +       1x,'   does not like to compute 2**j when j>30.')
 2100 format(1x,'WARNING: for some cases, a considerable performance'/
     +       1x,'   reduction has been observed when running a non-'/
     +       1x,'   optimal number of bits with the micro-GA.'/
     +       1x,'   If possible, use values for nposibl of 2**n,'/
     +       1x,'   e.g. 2, 4, 8, 16, 32, 64, etc.  See ReadMe file.')
c
      return
      end
c
c#######################################################################
      subroutine evalout(iskip,iend,ibest,fbar,best,bestssqq)
c
c  This subroutine evaluates the population, assigns fitness,
c  establishes the best individual, and outputs information.
      implicit real*8 (a-h,o-z)
      save
c
      include 'ga_params.h'
      dimension parent(nparmax,indmax),iparent(nchrmax,indmax)
      dimension fitness(indmax)
      dimension paramsm(nparmax),paramav(nparmax),ibest(nchrmax)
c
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga4   / fitness
 
      fitsum=0.0d0
      best=-1.0d10
      do 29 n=1,nparam
         paramsm(n)=0.0d0
 29   continue
      jstart=1
      jend=npopsiz
      if(iskip.ne.0) jstart=iskip
      if(iend.ne.0) jend=iend
      do 30 j=jstart,jend
         call decode(j,parent,iparent)
         if(iskip.ne.0 .and. iend.ne.0 .and. iskip.eq.iend)
     +   write(6,1075) j,(iparent(k,j),k=1,nchrome),
     +                   (parent(kk,j),kk=1,nparam),0.0
 
c  Call function evaluator, write out individual and fitness, and add
c  to the summation for later averaging.
         call func(j,funcval,ssqqval)
         fitness(j)=funcval
         write(24,1075) j,(iparent(k,j),k=1,nchrome),
     +                  (parent(kk,j),kk=1,nparam),fitness(j)
         write(24,1075) j,fitness(j),(parent(kk,j),kk=1,nparam)
         fitsum=fitsum+fitness(j)
         do 22 n=1,nparam
            paramsm(n)=paramsm(n)+parent(n,j)
 22      continue
 
c  Check to see if fitness of individual j is the best fitness.
         if (fitness(j).gt.best) then
            best=fitness(j)
  
            bestssqq=ssqqval
  
            jbest=j
            do 24 k=1,nchrome
               ibest(k)=iparent(k,j)
 24         continue
         endif
 30   continue
 
c  Compute parameter and fitness averages.
      fbar=fitsum/dble(npopsiz)
      do 23 n=1,nparam
         paramav(n)=paramsm(n)/dble(npopsiz)
 23   continue
 
c  Write output information
      if (npopsiz.eq.1) then
         write(24,1075) 1,(iparent(k,1),k=1,nchrome),
     +                  (parent(k,1),k=1,nparam),fitness(1)
c         write(24,1075) 1,fitness(1),(parent(kk,1),kk=1,nparam)
         write(24,*) ' Average Values:'
         write(24,1275) (parent(k,1),k=1,nparam),fbar
         write(24,'(6f12.6)') (parent(k,1),k=1,nparam)
      else
         write(24,*)' Average parameter values:'
         write(24,1275) (paramav(k),k=1,nparam),fbar
         write(24,'(6f12.6)') (paramav(k),k=1,nparam)
      endif
      write(6,1100) fbar
      write(24,1100) fbar
      write(6,1200) best
      write(24,1200) best
c
 1075 format(i3,1x,30i1,2(1x,f7.4),1x,f8.5)
c 1075 format(i5,'   fitness =',e16.5,/,(6f12.6))
 1100 format(1x,'Average fitness of generation =',e16.5)
c 1200 format(1x,'Maximum Function Value              =',f8.5/)
 1200 format(1x,'Best fitness so far           =',e16.5)
 1275 format(/' Average Values:',18x,2(1x,f7.4),1x,f8.5/)
      return
      end
c
c#######################################################################
      subroutine niche
c
c  Implement "niching" through Goldberg's multidimensional phenotypic
c  sharing scheme with a triangular sharing function.  To find the
c  multidimensional distance from the best individual, normalize all
c  parameter differences.
c
      implicit real*8 (a-h,o-z)
      save
c
      include 'ga_params.h'
      dimension parent(nparmax,indmax),iparent(nchrmax,indmax)
      dimension fitness(indmax),nposibl(nparmax),nichflg(nparmax)
      dimension parmax(nparmax),parmin(nparmax),pardel(nparmax)
c
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga4   / fitness
      common / ga6   / parmax,parmin,pardel,nposibl
      common / ga8   / nichflg
c
c   Variable definitions:
c
c  alpha   = power law exponent for sharing function; typically = 1.0
c  del     = normalized multidimensional distance between ii and all
c            other members of the population
c            (equals the square root of del2)
c  del2    = sum of the squares of the normalized multidimensional
c            distance between member ii and all other members of
c            the population
c  nniche  = number of niched parameters
c  sigshar = normalized distance to be compared with del; in some sense,
c            1/sigshar can be viewed as the number of regions over which
c            the sharing function should focus, e.g. with sigshar=0.1,
c            the sharing function will try to clump in ten distinct
c            regions of the phase space.  A value of sigshar on the
c            order of 0.1 seems to work best.
c  share   = sharing function between individual ii and j
c  sumshar = sum of the sharing functions for individual ii
c
c      alpha=1.0
      sigshar=0.1d0
      nniche=0
      do 33 jj=1,nparam
         nniche=nniche+nichflg(jj)
 33   continue
      if (nniche.eq.0) then
         write(6,1900)
         write(24,1900)
         close(24)
         stop
      endif
      do 34 ii=1,npopsiz
         sumshar=0.0d0
         do 35 j=1,npopsiz
            del2=0.0d0
            do 36 k=1,nparam
               if (nichflg(k).ne.0) then
                  del2=del2+((parent(k,j)-parent(k,ii))/pardel(k))**2
               endif
 36         continue
            del=(dsqrt(del2))/dble(nniche)
            if (del.lt.sigshar) then
c               share=1.0-((del/sigshar)**alpha)
               share=1.0d0-(del/sigshar)
            else
               share=0.0d0
            endif
            sumshar=sumshar+share/dble(npopsiz)
 35      continue
         if (sumshar.ne.0.0d0) fitness(ii)=fitness(ii)/sumshar
 34   continue
c
 1900 format(1x,'ERROR: iniche=1 and all values in nichflg array = 0'/
     +       1x,'       Do you want to niche or not?')
c
      return
      end
c
c#######################################################################
      subroutine selectn(ipick,j,mate1,mate2)
c
c  Subroutine for selection operator.  Presently, tournament selection
c  is the only option available.
c
      implicit real*8 (a-h,o-z)
      save
c
      include 'ga_params.h'
      dimension parent(nparmax,indmax),child(nparmax,indmax)
      dimension fitness(indmax)
      dimension iparent(nchrmax,indmax),ichild(nchrmax,indmax)
c
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga4   / fitness
      common / ga7   / child,ichild
      common /inputga/ pcross,pmutate,pcreep,maxgen,idum,irestrt,
     +                 itourny,ielite,icreep,iunifrm,iniche,
     +                 iskip,iend,nchild,microga,kountmx
c
c  If tournament selection is chosen (i.e. itourny=1), then
c  implement "tournament" selection for selection of new population.
      if(itourny.eq.1) then
         call select(mate1,ipick)
         call select(mate2,ipick)
c        write(3,*) mate1,mate2,fitness(mate1),fitness(mate2)
         do 46 n=1,nchrome
            ichild(n,j)=iparent(n,mate1)
            if(nchild.eq.2) ichild(n,j+1)=iparent(n,mate2)
 46      continue
      endif
c
      return
      end
c
c#######################################################################
      subroutine crosovr(ncross,j,mate1,mate2)
c
c  Subroutine for crossover between the randomly selected pair.
      implicit real*8 (a-h,o-z)
      save
c
      include 'ga_params.h'
      dimension parent(nparmax,indmax),child(nparmax,indmax)
      dimension iparent(nchrmax,indmax),ichild(nchrmax,indmax)
c
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga7   / child,ichild
      common /inputga/ pcross,pmutate,pcreep,maxgen,idum,irestrt,
     +                 itourny,ielite,icreep,iunifrm,iniche,
     +                 iskip,iend,nchild,microga,kountmx
c
      if (iunifrm.eq.0) then
c  Single-point crossover at a random chromosome point.
         call ran3(1,rand)
         if(rand.gt.pcross) goto 69
         ncross=ncross+1
         call ran3(1,rand)
         icross=2+dint(dble(nchrome-1)*rand)
         do 50 n=icross,nchrome
            ichild(n,j)=iparent(n,mate2)
            if(nchild.eq.2) ichild(n,j+1)=iparent(n,mate1)
 50      continue
      else
c  Perform uniform crossover between the randomly selected pair.
         do 60 n=1,nchrome
            call ran3(1,rand)
            if(rand.le.pcross) then
               ncross=ncross+1
               ichild(n,j)=iparent(n,mate2)
               if(nchild.eq.2) ichild(n,j+1)=iparent(n,mate1)
            endif
 60      continue
      endif
 69   continue
c
      return
      end
c
c#######################################################################
      subroutine mutate
c
      implicit real*8 (a-h,o-z)
      save
c
      include 'ga_params.h'
      dimension nposibl(nparmax)
      dimension child(nparmax,indmax),ichild(nchrmax,indmax)
      dimension g0(nparmax),g1(nparmax),ig2(nparmax)
      dimension parmax(nparmax),parmin(nparmax),pardel(nparmax)
c
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga5   / g0,g1,ig2
      common / ga6   / parmax,parmin,pardel,nposibl
      common / ga7   / child,ichild
      common /inputga/ pcross,pmutate,pcreep,maxgen,idum,irestrt,
     +                 itourny,ielite,icreep,iunifrm,iniche,
     +                 iskip,iend,nchild,microga,kountmx
c
c  This subroutine performs mutations on the children generation.
c  Perform random jump mutation if a random number is less than pmutate.
c  Perform random creep mutation if a different random number is less
c  than pcreep.

c     open(24,file='ga.out',status='old')      ! ccd ---------
      nmutate=0
      ncreep=0
      do 70 j=1,npopsiz
         do 75 k=1,nchrome
c  Jump mutation
            call ran3(1,rand)
            if (rand.le.pmutate) then
               nmutate=nmutate+1
               if(ichild(k,j).eq.0) then
                  ichild(k,j)=1
               else
                  ichild(k,j)=0
               endif
               if (nowrite.eq.0) write(6,1300) j,k
               if (nowrite.eq.0) write(24,1300) j,k
            endif
 75      continue
c  Creep mutation (one discrete position away).
         if (icreep.ne.0) then
            do 76 k=1,nparam
               call ran3(1,rand)
               if(rand.le.pcreep) then
                  call decode(j,child,ichild)
                  ncreep=ncreep+1
                  creep=1.0d0
                  call ran3(1,rand)
                  if (rand.lt.0.5d0) creep=-1.0d0
                  child(k,j)=child(k,j)+g1(k)*creep
                  if (child(k,j).gt.parmax(k)) then
                     child(k,j)=parmax(k)-1.0d0*g1(k)
                  elseif (child(k,j).lt.parmin(k)) then
                     child(k,j)=parmin(k)+1.0d0*g1(k)
                  endif
                  call code(j,k,child,ichild)
                  if (nowrite.eq.0) write(6,1350) j,k
                  if (nowrite.eq.0) write(24,1350) j,k
               endif
 76         continue
         endif
 70   continue
      write(6,1250) nmutate,ncreep
      write(24,1250) nmutate,ncreep
c
 1250 format(/'  Number of Jump Mutations  =',i5/
     +        '  Number of Creep Mutations =',i5)
 1300 format('*** Jump mutation performed on individual  ',i4,
     +       ', chromosome ',i3,' ***')
 1350 format('*** Creep mutation performed on individual ',i4,
     +       ', parameter  ',i3,' ***')
c
      return
      end
c
c#######################################################################
      subroutine newgen(ielite,npossum,ig2sum,ibest)
c
c  Write child array back into parent array for new generation.  Check
c  to see if the best parent was replicated; if not, and if ielite=1,
c  then reproduce the best parent into a random slot.
c
      implicit real*8 (a-h,o-z)
      save
c
      include 'ga_params.h'
      dimension parent(nparmax,indmax),child(nparmax,indmax)
      dimension iparent(nchrmax,indmax),ichild(nchrmax,indmax)
      dimension ibest(nchrmax)
c
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga7   / child,ichild
c
      if (npossum.lt.ig2sum) call possibl(child,ichild)
      kelite=0
      do 94 j=1,npopsiz
         jelite=0
         do 95 n=1,nchrome
            iparent(n,j)=ichild(n,j)
            if (iparent(n,j).eq.ibest(n)) jelite=jelite+1
            if (jelite.eq.nchrome) kelite=1
 95      continue
 94   continue
      if (ielite.ne.0 .and. kelite.eq.0) then
         call ran3(1,rand)
         irand=1d0+dint(dble(npopsiz)*rand)
         do 96 n=1,nchrome
            iparent(n,irand)=ibest(n)
 96      continue
         write(24,1260) irand
      endif
c
 1260 format('  Elitist Reproduction on Individual ',i4)
c
      return
      end
c
c#######################################################################
      subroutine gamicro(i,npossum,ig2sum,ibest)
c
c  Micro-GA implementation subroutine
c
      implicit real*8 (a-h,o-z)
      save
c
      include 'ga_params.h'
      dimension parent(nparmax,indmax),iparent(nchrmax,indmax)
      dimension ibest(nchrmax)
c
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
c
c  First, check for convergence of micro population.
c  If converged, start a new generation with best individual and fill
c  the remainder of the population with new randomly generated parents.
c
c  Count number of different bits from best member in micro-population
      icount=0
      do 81 j=1,npopsiz
         do 82 n=1,nchrome
            if(iparent(n,j).ne.ibest(n)) icount=icount+1
 82      continue
 81   continue
c
c  If icount less than 5% of number of bits, then consider population
c  to be converged.  Restart with best individual and random others.
      diffrac=dble(icount)/dble((npopsiz-1)*nchrome)
      if (diffrac.lt.0.05d0) then
      do 87 n=1,nchrome
         iparent(n,1)=ibest(n)
 87   continue
      do 88 j=2,npopsiz
         do 89 n=1,nchrome
            call ran3(1,rand)
            iparent(n,j)=1
            if(rand.lt.0.5d0) iparent(n,j)=0
 89      continue
 88   continue
      if (npossum.lt.ig2sum) call possibl(parent,iparent)
      write(6,1375) i
      write(24,1375) i
      endif
c
 1375 format(//'%%%%%%%  Restart micro-population at generation',
     +       i5,'  %%%%%%%')
c
      return
      end
c
c#######################################################################
      subroutine select(mate,ipick)
c
c  This routine selects the better of two possible parents for mating.
c
      implicit real*8 (a-h,o-z)
      save
c
      include 'ga_params.h'
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga4   / fitness
      dimension parent(nparmax,indmax),iparent(nchrmax,indmax)
      dimension fitness(indmax)
c
      if(ipick+1.gt.npopsiz) call shuffle(ipick)
      ifirst=ipick
      isecond=ipick+1
      ipick=ipick+2
      if(fitness(ifirst).gt.fitness(isecond)) then
         mate=ifirst
      else
         mate=isecond
      endif
c     write(3,*)'select',ifirst,isecond,fitness(ifirst),fitness(isecond)
c
      return
      end
c
c#######################################################################
      subroutine shuffle(ipick)
c
c  This routine shuffles the parent array and its corresponding fitness
c
      implicit real*8 (a-h,o-z)
      save
c
      include 'ga_params.h'
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga4   / fitness
      dimension parent(nparmax,indmax),iparent(nchrmax,indmax)
      dimension fitness(indmax)
c
      ipick=1
      do 10 j=1,npopsiz-1
         call ran3(1,rand)
         iother=j+1+dint(dble(npopsiz-j)*rand)
         do 20 n=1,nchrome
            itemp=iparent(n,iother)
            iparent(n,iother)=iparent(n,j)
            iparent(n,j)=itemp
 20      continue
         temp=fitness(iother)
         fitness(iother)=fitness(j)
         fitness(j)=temp
 10   continue
c
      return
      end
c
c#######################################################################
      subroutine decode(i,array,iarray)
c
c  This routine decodes a binary string to a real number.
c
      implicit real*8 (a-h,o-z)
      save
c
      include 'ga_params.h'
      common / ga2   / nparam,nchrome
      common / ga5   / g0,g1,ig2
      dimension array(nparmax,indmax),iarray(nchrmax,indmax)
      dimension g0(nparmax),g1(nparmax),ig2(nparmax)
c
      l=1
      do 10 k=1,nparam
         iparam=0
         m=l
         do 20 j=m,m+ig2(k)-1
            l=l+1
            iparam=iparam+iarray(j,i)*(2**(m+ig2(k)-1-j))
 20      continue
         array(k,i)=g0(k)+g1(k)*dble(iparam)
 10   continue
c
      return
      end
c
c#######################################################################
      subroutine code(j,k,array,iarray)
c
c  This routine codes a parameter into a binary string.
c
      implicit real*8 (a-h,o-z)
      save
c
      include 'ga_params.h'
      common / ga2   / nparam,nchrome
      common / ga5   / g0,g1,ig2
      dimension array(nparmax,indmax),iarray(nchrmax,indmax)
      dimension g0(nparmax),g1(nparmax),ig2(nparmax)
c
c  First, establish the beginning location of the parameter string of
c  interest.
      istart=1
      do 10 i=1,k-1
         istart=istart+ig2(i)
 10   continue
c
c  Find the equivalent coded parameter value, and back out the binary
c  string by factors of two.
      m=ig2(k)-1
      if (g1(k).eq.0.0d0) return
      iparam=nint((array(k,j)-g0(k))/g1(k))
      do 20 i=istart,istart+ig2(k)-1
         iarray(i,j)=0
         if ((iparam+1).gt.(2**m)) then
            iarray(i,j)=1
            iparam=iparam-2**m
         endif
         m=m-1
 20   continue
c     write(3,*)array(k,j),iparam,(iarray(i,j),i=istart,istart+ig2(k)-1)
c
      return
      end
c
c#######################################################################
c
      subroutine possibl(array,iarray)
c
c  This subroutine determines whether or not all parameters are within
c  the specified range of possibility.  If not, the parameter is
c  randomly reassigned within the range.  This subroutine is only
c  necessary when the number of possibilities per parameter is not
c  optimized to be 2**n, i.e. if npossum < ig2sum.
c
      implicit real*8 (a-h,o-z)
      save
c
      include 'ga_params.h'
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga5   / g0,g1,ig2
      common / ga6   / parmax,parmin,pardel,nposibl
      dimension array(nparmax,indmax),iarray(nchrmax,indmax)
      dimension g0(nparmax),g1(nparmax),ig2(nparmax),nposibl(nparmax)
      dimension parmax(nparmax),parmin(nparmax),pardel(nparmax)
c
      do 10 i=1,npopsiz
         call decode(i,array,iarray)
         do 20 j=1,nparam
            n2ig2j=2**ig2(j)
            if(nposibl(j).ne.n2ig2j .and. array(j,i).gt.parmax(j)) then
               call ran3(1,rand)
               irand=dint(dble(nposibl(j))*rand)
               array(j,i)=g0(j)+dble(irand)*g1(j)
               call code(i,j,array,iarray)
               if (nowrite.eq.0) write(6,1000) i,j
                       if (nowrite.eq.0) write(24,1000) i,j
            endif
 20      continue
 10   continue
c
 1000 format('*** Parameter adjustment to individual     ',i4,
     +       ', parameter  ',i3,' ***')
c
      return
      end
c
c#######################################################################
      subroutine restart(i,istart,kount)
c
c  This subroutine writes restart information to the ga.restart file.
c
      implicit real*8 (a-h,o-z)
      save
c
      include 'ga_params.h'
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      dimension parent(nparmax,indmax),iparent(nchrmax,indmax)
      common /inputga/ pcross,pmutate,pcreep,maxgen,idum,irestrt,
     +                 itourny,ielite,icreep,iunifrm,iniche,
     +                 iskip,iend,nchild,microga,kountmx

      kount=kount+1
      if(i.eq.maxgen+istart-1 .or. kount.eq.kountmx) then
         OPEN (UNIT=25, FILE='ga.restart', STATUS='OLD')
         rewind 25
         write(25,*) i+1,npopsiz
         do 80 j=1,npopsiz
            write(25,1500) j,(iparent(l,j),l=1,nchrome)
 80      continue
         CLOSE (25)
         kount=0
      endif
c
 1500 format(i5,3x,30i2)
c
      return
      end
c
c#######################################################################
      subroutine ran3(idum,rand)
c
c  Returns a uniform random deviate between 0.0 and 1.0.  Set idum to
c  any negative value to initialize or reinitialize the sequence.
c  This function is taken from W.H. Press', "Numerical Recipes" p. 199.
c
      implicit real*8 (a-h,m,o-z)
             save
c      implicit real*4(m)
      parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=1./mbig)
c     parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)
c
c  According to Knuth, any large mbig, and any smaller (but still large)
c  mseed can be substituted for the above values.
      dimension ma(55)
      data iff /0/
      if (idum.lt.0 .or. iff.eq.0) then
         iff=1
         mj=mseed-dble(iabs(idum))
         mj=dmod(mj,mbig)
         ma(55)=mj
         mk=1
         do 11 i=1,54
            ii=mod(21*i,55)
            ma(ii)=mk
            mk=mj-mk
            if(mk.lt.mz) mk=mk+mbig
            mj=ma(ii)
 11      continue
         do 13 k=1,4
            do 12 i=1,55
               ma(i)=ma(i)-ma(1+mod(i+30,55))
               if(ma(i).lt.mz) ma(i)=ma(i)+mbig
 12         continue
 13      continue
         inext=0
         inextp=31
         idum=1
      endif
      inext=inext+1
      if(inext.eq.56) inext=1
      inextp=inextp+1
      if(inextp.eq.56) inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz) mj=mj+mbig
      ma(inext)=mj
      rand=mj*fac
      return
      end
c
c#######################################################################
