      program convert_rsf_files
C**** convert rsf files from model II' (B399) format to modelE format 
C**** compile with: gmake conv_rsf.o
C**** f90 -o conv_rsf conv_rsf.o *.o -O2 -64 -mips4 -static \
C****                      -OPT:reorg_comm=off -w2 -listing
C**** Note that since it uses modules and routines from the model, it
C**** must be compiled after the model
      USE CONSTANT, only : lhm,shi
      USE E001M12_COM, only : im,jm,lm,wm,u,v,t,p,q,jc,rc,clabel
     *     ,iowrite_mon 
      USE SOMTQ_COM
      USE GHYCOM, only : ghdata,snowe,tearth,wearth,aiearth,snoage
      USE RADNCB, only : rqt,lm_req
      USE CLD01_COM_E001, only : ttold,qtold,svlhx,rhsav,cldsav
      USE DAGCOM, only : keynr,tsfrez
      USE PBLCOM, only : uabl,vabl,tabl,qabl,eabl,cm=>cmgs,ch=>chgs,cq
     *     =>cqgs,ipbl,wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg,ustar
     *     ,egcm
      USE OCEAN, only : tocean,z1o
      USE SEAICE_COM, only : rsi,msi,hsi,snowi
      USE SEAICE, only : ace1i,xsi,ac2oim
      USE LANDICE_COM, only : tlandi,snowli
      USE LAKES_COM, only : flake
      IMPLICIT NONE
      CHARACTER infile*60, outfile*60
      INTEGER IARGC,iu_AIC,I,J,L,N,ioerr,iu_TOPO
      REAL*8 TAUX,X   ! ? temporary for compatibility only
      REAL*8 MSI1
      INTEGER ItimeX
!@ egcm_init_max maximum initial vaule of egcm
      real*8, parameter :: egcm_init_max=0.5
      real*8 :: tmp   ! temperary variable

      
      IF (IARGC().lt.2) THEN
        PRINT*,"Convert rsf files from old format to new"
        PRINT*,"conv_rsf filename output_file"
        STOP
      END IF

      CALL GETARG(1,infile)
      CALL GETARG(2,outfile)

      iu_AIC=9
      OPEN(iu_AIC,FILE=trim(infile),FORM="UNFORMATTED",STATUS="OLD")

      READ (iu_AIC,ERR=800,END=810) TAUX,JC,CLABEL,RC,KEYNR,
     *     U,V,T,P,Q,
     2     ((TOCEAN(1,I,J),I=1,IM),J=1,JM),RSI,MSI,
     *     (((TOCEAN(L,I,J),I=1,IM),J=1,JM),L=2,3),
     *     SNOWI,SNOWE,
     *     ((HSI(1,I,J),I=1,IM),J=1,JM),TEARTH,WEARTH,AIEARTH,
     *     ((HSI(2,I,J),I=1,IM),J=1,JM),((X,I=1,IM),J=1,JM),
     *     (((SNOAGE(L,I,J),I=1,IM),J=1,JM),L=1,3),SNOWLI,
     *     (((TLANDI(L,I,J),I=1,IM),J=1,JM),L=1,2),
     *     (((HSI(L,I,J),I=1,IM),J=1,JM),L=3,4),GHDATA,
     *     wsavg,tsavg,qsavg,dclev,Z1O,usavg,vsavg,tauavg,ustar,
     *     uabl,vabl,tabl,qabl,eabl,cm,ch,cq,ipbl,
     A     (((TTOLD(L,I,J),I=1,IM),J=1,JM),L=1,LM),
     B     (((QTOLD(L,I,J),I=1,IM),J=1,JM),L=1,LM),
     C     (((SVLHX(L,I,J),I=1,IM),J=1,JM),L=1,LM),
     D     (((RHSAV(L,I,J),I=1,IM),J=1,JM),L=1,LM),WM,
     E     (((CLDSAV(L,I,J),I=1,IM),J=1,JM),L=1,LM),
     4     ((((TMOM(N,I,J,L),I=1,IM),J=1,JM),L=1,LM),N=1,9),
     4     ((((QMOM(N,I,J,L),I=1,IM),J=1,JM),L=1,LM),N=1,9),
     6     (((RQT(L,I,J),I=1,IM),J=1,JM),L=1,LM_REQ)
      CLOSE (iu_AIC)

      ItimeX=NINT(TAUX)
      print*,ItimeX

C**** convert sea ice temperatures into enthalpy
      DO J=1,JM
        DO I=1,IM
          IF (RSI(I,J).gt.0) THEN
            MSI1=SNOWI(I,J)+ACE1I
            HSI(1,I,J) = (SHI*MIN(HSI(1,I,J),0d0)-LHM)*XSI(1)*MSI1
            HSI(2,I,J) = (SHI*MIN(HSI(2,I,J),0d0)-LHM)*XSI(2)*MSI1
            HSI(3,I,J) = (SHI*MIN(HSI(3,I,J),0d0)-LHM)*XSI(3)*MSI(I,J)
            HSI(4,I,J) = (SHI*MIN(HSI(4,I,J),0d0)-LHM)*XSI(4)*MSI(I,J)
          ELSE
            MSI(I,J)=AC2OIM
            SNOWI(I,J)=0.
            HSI(1,I,J) = -LHM*XSI(1)*ACE1I
            HSI(2,I,J) = -LHM*XSI(2)*ACE1I
            HSI(3,I,J) = -LHM*XSI(3)*AC2OIM
            HSI(4,I,J) = -LHM*XSI(4)*AC2OIM
          END IF
        END DO
      END DO

c     initialize the 3-d turbulent kinetic enery to be used in
c     the subroutine diffus.
      do l=1,lm
        tmp=egcm_init_max/(float(l)**2)
        do j=1,jm
          do i=1,im
            egcm(i,j,l)=tmp
          end do
        end do
      end do

C**** initialize TSFREZ to defaults
      TSFREZ(:,:,1:2)=365.
      TSFREZ(:,:,3:4)=-999.

C**** read in FLAKE data
      iu_TOPO=10
      OPEN(iu_TOPO,FILE="/u/cmrun/Z72X46N.cor4",FORM="UNFORMATTED"
     *     ,STATUS="OLD") 
      READ (iu_TOPO)
      CALL READT (iu_TOPO,0,FLAKE,IM*JM,FLAKE,1) ! Lake fraction
      close (iu_TOPO)

      OPEN(iu_AIC,FILE=trim(outfile),
     *     FORM="UNFORMATTED",STATUS="UNKNOWN")

      call io_rsf(iu_AIC,ItimeX,iowrite_mon,ioerr)
      close (iu_AIC)
      print*,ioerr
      print*,"New rsf file written out to ",trim(outfile)
      stop
 800  print*,"Error reading in file"
 810  print*,"Error reading in file"
      stop
      end 
