      subroutine prepfield(ind)                                      
      implicit double precision (a-h,o-z)                            
                                                                     
      include 'ga_params.h'                                          
                                                                     
      common /evbpcon/ nv,ntmol,natms,numoen                         
      common /ecvpran/ oexmax, oexmin, pexmaxr, pexminr              
      common /ga3/ parent(nparmax,indmax),iparent(nchrmax,indmax)    
                                                                     
      open(80,file='new_params.dat')                                 
                                                                     
      i0 = 0                                                         
      i1 = 1                                                         
      i2 = 2                                                         
      i3 = 3                                                         
      zero = 0.d0                                                    
      one = 1.d0                                                     
      two = 2.d0                                                     
      three = 3.d0                                                   
      pi = 180.d0                                                    
                                                                     
c     Parameters.                                                    
        write(80,'(4f15.6)') parent(   1:   3,ind), zero
        write(80,'(4f15.6)') parent(   4:   6,ind), zero
        write(80,'(4f15.6)') parent(   7:   9,ind), zero
        write(80,'(4f15.6)') parent(  10:  12,ind), zero
        write(80,'(4f15.6)') parent(  13:  15,ind), zero
        write(80,'(4f15.6)') parent(  16:  17,ind), zero, zero
        write(80,'(4f15.6)') parent(  18:  19,ind), zero, zero
        write(80,'(4f15.6)') parent(  20:  21,ind), zero, zero
        write(80,'(4f15.6)') parent(  22:  23,ind), zero, zero
        write(80,'(4f15.6)') parent(  24:  25,ind), zero, zero
        write(80,'(4f15.6)') parent(  26:  27,ind), zero, zero
        write(80,'(4f15.6)') parent(  28:  29,ind), zero, zero
        write(80,'(4f15.6)') parent(  30:  31,ind), zero, zero
        write(80,'(4f15.6)') parent(  32:  33,ind), zero, zero
        write(80,'(4f15.6)') parent(  34:  35,ind), zero, zero
        write(80,'(4f15.6)') parent(  36:  37,ind), zero, zero
        write(80,'(4f15.6)') parent(  38:  39,ind), zero, zero
        write(80,'(4f15.6)') parent(  40:  41,ind), zero, zero
        write(80,'(4f15.6)') parent(  42:  43,ind), zero, zero
        write(80,'(4f15.6)') parent(  44:  45,ind), zero, zero
        write(80,'(4f15.6)') parent(  46:  47,ind), zero, zero
         write(80,'(3f15.6,i10)') parent(  48,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  49,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  50,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  51,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  52,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  53,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  54,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  55,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  56,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  57,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  58,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  59,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  60,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  61,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  62,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  63,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  64,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  65,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  66,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  67,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  68,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  69,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  70,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  71,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  72,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  73,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  74,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  75,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  76,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  77,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  78,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  79,ind), 180.0,   two, 0
         write(80,'(3f15.6,i10)') parent(  80,ind), 180.0,   two, 0
      close(80)

      return
      end
