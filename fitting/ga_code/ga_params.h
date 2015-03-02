c        Header file for GA

      parameter (indmax=2000,nparmax=1000)
      parameter (maxbits_per_param=30)  ! must be <= 30
      parameter (nchrmax=maxbits_per_param*nparmax)
      parameter (maxnv=1000)

c  indmax  = maximum # of individuals, i.e. max population size
c            (individual means parameter set)
c  maxbits_per_param = maximum # of bits per parameter
c  nchrmax = maximum # of chromosomes (binary bits) per individual
c  nparmax = maximum # of parameters which the chromosomes make up
