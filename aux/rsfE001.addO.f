C**** Combine final restart file of PRESCRIBED OCEAN DATA model run
C**** with mean & bottom temperature of 2nd ocean layer to create
C**** Initial Conditions for a PREDICTED OCEAN DATA model run.
C**** Input:  unit  1 - restart file of a run with prescribed ODATA
C****         unit 27 - mixed lyr depth,TO2avg,TO2bot
C**** Output: unit  9 - augmented restart file
C****

      USE MODEL_COM, only: im,jm,lm,iowrite_mon,irerun  
      USE TIMINGS, only : ntimeacc,timing,timestr 
      USE OCEAN, only : tocean
      IMPLICIT NONE
      INTEGER IARGC,iu_AIC,I,J,L,N,ioerr
      INTEGER ItimeX
      REAL*8 O1(IM,JM,5)

      iu_AIC = 1
      call io_rsf(iu_AIC,ItimeX,irerun,ioerr)
      close (iu_AIC)

      READ(27) O1
C* 
C***  Set the ocean temperature below the mixed layer 
C* 
      do j = 1,jm 
        do i = 1,im 
          tocean(2,i,j) = o1(i,j,4) 
          tocean(3,i,j) = o1(i,j,5) 
        end do 
      end do
C* 
      ntimeacc = 1
      timing = 0 
      timestr = " " 
C* 
      iu_AIC = 9      
      call io_rsf(iu_AIC,ItimeX,iowrite_mon,ioerr)
      close (iu_AIC)

      WRITE(0,*) ' NORMAL END'
      END










