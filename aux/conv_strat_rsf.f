      program convert_rsf_files
C**** convert rsf files from model II' (B399) format to modelE format
C**** compile with: gmake conv_rsf.o
C**** f90 -o conv_rsf conv_rsf.o *.o -O2 -64 -mips4 -static \
C****                      -OPT:reorg_comm=off -w2 -listing
C**** Note that since it uses modules and routines from the model, it
C**** must be compiled after the model
      USE CONSTANT, only : lhm,shi
      USE MODEL_COM, only : im,jm,lm,wm,u,v,t,p,q,xlabel
     *     ,iowrite_mon,focean,airx,lmc
     *     ,nday,itime,itimei,itimee,itime0,iyear1
      USE SOMTQ_COM
      USE GHYCOM, only : snowe,tearth,wearth,aiearth,snoage,wbare,wvege
     *     ,htbare,htvege,snowbv,ngm
      USE RADNCB, only : rqt,lm_req
      USE CLD01_COM_E001, only : ttold,qtold,svlhx,rhsav,cldsav
      USE DAGCOM, only : keynr,tsfrez
      USE PBLCOM, only : uabl,vabl,tabl,qabl,eabl,cm=>cmgs,ch=>chgs,cq
     *     =>cqgs,ipbl,wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg,ustar
     *     ,egcm
      USE OCEAN, only : tocean,z1o,tfo
      USE SEAICE_COM, only : rsi,msi,hsi,snowi,ssi
      USE SEAICE, only : ace1i,xsi,ac2oim,ssi0
      USE LANDICE_COM, only : tlandi,snowli
      USE LAKES_COM, only : flake
      IMPLICIT NONE
      CHARACTER infile*60, outfile*60               ,clabel*156
      INTEGER IARGC,iu_AIC,I,J,L,N,ioerr,iu_TOPO    ,jc(100)
      REAL*8 TAUX,X                                 ,rc(169)
      REAL*8 MSI1
      INTEGER ItimeX
!@ egcm_init_max maximum initial vaule of egcm
      real*8, parameter :: egcm_init_max=0.5

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
     *     (((HSI(L,I,J),I=1,IM),J=1,JM),L=3,4),
     *     (((wbare(L,I,J),I=1,IM),J=1,JM),L=1,NGM),
     *     (((wvege(L,I,J),I=1,IM),J=1,JM),L=0,NGM),
     *     (((htbare(L,I,J),I=1,IM),J=1,JM),L=0,NGM),
     *     (((htvege(L,I,J),I=1,IM),J=1,JM),L=0,NGM),
     *     (((snowbv(L,I,J),I=1,IM),J=1,JM),L=1,2),
     *     wsavg,tsavg,qsavg,dclev,Z1O,usavg,vsavg,tauavg
     *     ,ustar,uabl,vabl,tabl,qabl,eabl,cm,ch,cq,ipbl,(((TTOLD(L,I,J)
     *     ,I=1,IM),J=1,JM),L=1,LM),(((QTOLD(L,I,J),I=1,IM),J=1,JM),L=1
     *     ,LM),(((SVLHX(L,I,J),I=1,IM),J=1,JM),L=1,LM),(((RHSAV(L,I,J)
     *     ,I=1,IM),J=1,JM),L=1,LM),WM,(((CLDSAV(L,I,J),I=1,IM),J=1,JM)
     *     ,L=1,LM),((((TMOM(N,I,J,L),I=1,IM),J=1,JM),L=1,LM),N=1,9)
     *     ,((((QMOM(N,I,J,L),I=1,IM),J=1,JM),L=1,LM),N=1,9),(((RQT(L,I
     *     ,J),I=1,IM),J=1,JM),L=1,LM_REQ)
     *     ,(x,I=1,IM*JM*(2*lm+8)),AIRX,
     *    (((LMC(l,i,j),i=1,IM),J=1,JM),L=1,2)
      CLOSE (iu_AIC)

      XLABEL=CLABEL(1:132)
      ItimeX=NINT(TAUX)
      NDAY = 24
      IYEAR1 = JC(41)
      itimei = itimex
      itimee = itimex
      itime0 = itimex
      print*,iyear1,ItimeX,xlabel

C**** read in FLAKE/FOCEAN data
      iu_TOPO=10
      OPEN(iu_TOPO,FILE="/u/cmrun/Z72X46N.cor4",FORM="UNFORMATTED"
     *     ,STATUS="OLD")
      CALL READT (iu_TOPO,0,FOCEAN,IM*JM,FOCEAN,1) ! ocean fraction
      CALL READT (iu_TOPO,0,FLAKE,IM*JM,FLAKE,1) ! Lake fraction
      close (iu_TOPO)

C**** convert sea ice temperatures into enthalpy
C**** and initialize sea ice salinity to 3.2 ppt (0 in snow & lake ice).
      DO J=1,JM
        DO I=1,IM
          IF (RSI(I,J).gt.0) THEN
            MSI1=SNOWI(I,J)+ACE1I
            IF (FOCEAN(I,J).gt.0) THEN
              HSI(1:2,I,J)=(SHI*MIN(HSI(1:2,I,J),TFO)-LHM)*XSI(1:2)*MSI1
              HSI(3:4,I,J)=(SHI*MIN(HSI(3:4,I,J),TFO)-LHM)*XSI(3:4)
     *             *MSI(I,J)
              IF (ACE1I*XSI(1).gt.SNOWI(I,J)*XSI(2)) THEN
                SSI(1,I,J)=SSI0 * (ACE1I-(ACE1I+SNOWI(I,J))* XSI(2))
                SSI(2,I,J)=SSI0 * (ACE1I+SNOWI(I,J))* XSI(2)
              ELSE
                SSI(1,I,J)=0.
                SSI(2,I,J)=SSI0 * ACE1I
              END IF
            ELSE
              HSI(1:2,I,J)=(SHI*MIN(HSI(1:2,I,J),0d0)-LHM)*XSI(1:2)*MSI1
              HSI(3:4,I,J)=(SHI*MIN(HSI(3:4,I,J),0d0)-LHM)*XSI(3:4)
     *             *MSI(I,J)
              SSI(1:4,I,J) = 0
            END IF
          ELSE
            MSI(I,J)=AC2OIM
            SNOWI(I,J)=0.
            IF (FOCEAN(I,J).gt.0) THEN
              HSI(1:2,I,J) = (SHI*TFO-LHM)*XSI(1:2)*ACE1I
              HSI(3:4,I,J) = (SHI*TFO-LHM)*XSI(3:4)*AC2OIM
              SSI(1:2,I,J)=SSI0 * ACE1I  * XSI(1:2)
              SSI(3:4,I,J)=SSI0 * AC2OIM * XSI(3:4)
            ELSE
              HSI(1:2,I,J) = -LHM*XSI(1:2)*ACE1I
              HSI(3:4,I,J) = -LHM*XSI(3:4)*AC2OIM
              SSI(1:4,I,J) = 0.
            END IF
          END IF
        END DO
      END DO


c     initialize the 3-d turbulent kinetic enery to be used in
c     the subroutine diffus.
      do j=1,jm
        do i=1,im
          do l=1,lm
            egcm(l,i,j)=egcm_init_max/(float(l)**2)
          end do
        end do
      end do

C**** initialize TSFREZ to defaults
      TSFREZ(:,:,1:2)=365.
      TSFREZ(:,:,3:4)=-999.


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
