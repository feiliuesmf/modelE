C**** Combine final restart file of PRESCRIBED OCEAN DATA model run
C**** with mean & bottom temperature of 2nd ocean layer to create
C**** Initial Conditions for a PREDICTED OCEAN DATA model run.
C**** Input:  unit  1 - restart file of a run with prescribed ODATA
C****         unit 27 - mixed lyr depth,TO2avg,TO2bot
C**** Output: unit  9 - augmented restart file
C****

c      PARAMETER (IM=72,JM=46,LM=12, KP1=4*LM+1,
c     *    KP2=16 + 4*6+5 + 12 + 174 + 6*LM + 18*LM + 3 ,
c     *    KP3=2*(LM+1)+ 4 + 2)
CB05  PARAMETER (IM=12,JM=24,LM=9, KP1=4*LM+1+3,KP2=14+8+3)
c      CHARACTER*4 C(39)
c      DIMENSION JC(100),KY(2100)
c      REAL*8 TAU,RC(161),U(IM,JM,KP1),O(IM,JM,5),G(IM,JM,KP2),
c     *       S(IM,JM,KP3),T50(IM,JM) 
c      REAL*8 TAUX

c      USE E001M12_COM, only : im,jm,iowrite_mon,ioread_mon
      use model_com, only: im,jm 
      USE E001M12_COM, only : lm,wm,u,v,t,p,q,jc,rc,clabel
     *     ,iowrite_mon,irerun  
      use timings, only : ntimeacc,timing,timestr 
      USE SOMTQ_COM
      USE GHYCOM, only : ghdata,snowe,tearth,wearth,aiearth,snoage
      USE RADNCB, only : rqt,lm_req
      USE CLD01_COM_E001, only : ttold,qtold,svlhx,rhsav,cldsav
      USE DAGCOM, only : keynr
      USE PBLCOM, only : uabl,vabl,tabl,qabl,eabl,cm=>cmgs,ch=>chgs,cq
     *     =>cqgs,ipbl,wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg,ustar
c      USE OCEAN, only : tocean 
      USE OCEAN, only : tocean,z1o
      USE SEAICE_COM, only : rsi,msi,hsi,snowi
      USE LANDICE_COM, only : tlandi,snowli
c      USE LAKES_COM, only : t50
      IMPLICIT NONE
c      CHARACTER infile*60, outfile*60
      INTEGER IARGC,iu_AIC,I,J,L,N,ioerr
      INTEGER ItimeX
      REAL*8 O1(IM,JM,5)

c      CALL GETARG(1,infile)
c      CALL GETARG(2,outfile)

      iu_AIC = 1
      call io_rsf(iu_AIC,ItimeX,irerun,ioerr)
      close (iu_AIC)

c      READ(1)  TAU,JC,C,RC,KY,U,O,G,S,TAUX
c      IF(TAU.NE.TAUX) THEN
c         WRITE (0,*) TAU,TAUX
c         STOP
c      ENDIF
C* 
      READ(27) O1
c      WRITE(9) TAU,JC,C,RC,KY,U,O,G,T50,S 
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










