      MODULE LAKES
!@sum  LAKES subroutines for Lakes and Rivers
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0 (based on LB265)
      USE CONSTANT, only : grav,bygrav,shw,rhow,lhm,shi
      USE MODEL_COM, only : IM,JM
      IMPLICIT NONE
      SAVE
C****
C**** Changes from Model III: MO -> MWL (kg), G0M -> GML (J),
C****                         GZM -> TLAKE (deg C)
!@var KDIREC directions for river flow
C**** (0 no flow, 1-8 anti-clockwise from top RH corner
      INTEGER, DIMENSION(IM,JM) :: KDIREC
!@var IDPOLE,JDPOLE special directions for south pole flow
      INTEGER :: IDPOLE,JDPOLE
!@var RATE rate of river flow downslope (fraction)
!@var DHORZ horizontal distance to downstream box (m)
      REAL*8, DIMENSION(IM,JM) :: RATE,DHORZ
!@var IFLOW,JFLOW grid box indexes for downstream direction
      INTEGER IFLOW(IM,JM),JFLOW(IM,JM)
      INTEGER, PARAMETER :: NRVR = 42 !@param Number of specified rivers
!@var IRVRMTH,JRVRMTH indexes for specified river mouths
      INTEGER, DIMENSION(NRVR) :: IRVRMTH,JRVRMTH

!@param MINMLD minimum mixed layer depth in lake (m)
      REAL*8, PARAMETER :: MINMLD = 1.
!@param TMAXRHO temperature of maximum density (pure water) (C)
      REAL*8, PARAMETER :: TMAXRHO = 4.
!@param KVLAKE lake diffusion constant at mixed layer depth (m^2/s)
      REAL*8, PARAMETER :: KVLAKE = 1d-5
!@param TFL freezing temperature for lakes (=0 C)
      REAL*8, PARAMETER :: TFL = 0.
!@param AC1LMIN, AC2LMIN minimum ice thickness for lake ice (kg/m^2)
      REAL*8, PARAMETER :: AC1LMIN = 0.1, AC2LMIN=0.1  ! (not used)
!@param FLEADLK lead fraction for lakes
      REAL*8, PARAMETER :: FLEADLK = 0.
!@param BYZETA reciprocal of solar rad. extinction depth for lake (1/m)
      REAL*8, PARAMETER :: BYZETA = 1./0.35d0

      CONTAINS

      SUBROUTINE LKSOURC (I0,J0,ROICE,MLAKE,ELAKE,RUN0,F0DT,F2DT,SROX
     *     ,FSR2,EVAPO,ENRGFO,ACEFO,ACEFI,ENRGFI)
!@sum  LKSOURC applies fluxes to lake in ice-covered and ice-free areas
!@auth Gary Russell
!@ver  1.0
      IMPLICIT NONE
!@var MLAKE,ELAKE mass and energy in lake layers (kg,J /m^2)
      REAL*8, INTENT(INOUT), DIMENSION(2) :: MLAKE,ELAKE
      INTEGER, INTENT(IN) :: I0,J0
      REAL*8, INTENT(IN) :: ROICE, EVAPO, RUN0
      REAL*8, INTENT(IN) :: F0DT, F2DT, SROX(2)
      REAL*8, INTENT(OUT) :: ENRGFO, ACEFO, ENRGFI, ACEFI
      REAL*8 ENRGF1, ACEF1, ENRGF2, ACEF2, FHO, FHI, FH0, FH1, FH2, FSR2
      REAL*8 ENRGI, ENRGI2, ENRGO, ENRGO2, RUNO, RUNI, TLK2, DM2, DH2
      REAL*8 FRATO,FRATI,E2O,E2I
C**** initiallize output
      ENRGFO=0. ; ACEFO=0. ; ACEFI=0. ; ENRGFI=0.

C**** Calculate heat and mass fluxes to lake
      ENRGO = F0DT-SROX(1)*FSR2 ! in open water
      ENRGO2=     +SROX(1)*FSR2 ! in open water, second layer
      ENRGI = F2DT-SROX(2)*FSR2 ! under ice
      ENRGI2=     +SROX(2)*FSR2 ! under ice, second layer
      RUNO  =-EVAPO
      RUNI  = RUN0
C**** Bring up mass from second layer if required/allowed
      IF (MLAKE(1)+RUNO.lt.MINMLD*RHOW.and.MLAKE(2).gt.0) THEN
        DM2 = MIN(MLAKE(2),MINMLD*RHOW-(MLAKE(1)+RUNO))
        DH2 = DM2*(ELAKE(2)+(1.-ROICE)*ENRGO2+ROICE*ENRGI2)/MLAKE(2)
      ELSE
        DM2 = 0.
        DH2 = 0.
      END IF

C**** Apply fluxes to 2nd layer
      IF (DM2.lt.MLAKE(2)) THEN
        MLAKE(2)=MLAKE(2) - DM2
        ELAKE(2)=ELAKE(2) - DH2 + (1.-ROICE)*ENRGO2 + ROICE*ENRGI2
      ELSE
        MLAKE(2)=0.
        ELAKE(2)=0.
      END IF

      E2O = 0. ; E2I = 0.

C**** Calculate energy in mixed layer (open ocean)
      IF (ROICE.LT.1d0) THEN
        FHO=ELAKE(1)+ENRGO+DH2-(MLAKE(1)+DM2+RUNO)*TFL*SHW
        IF (FHO.LT.0) THEN ! FLUXES COOL WATER TO FREEZING, FORM ICE
          ACEFO =FHO/(TFL*(SHI-SHW)-LHM)
          ACEFO =MIN(ACEFO,MAX(MLAKE(1)+DM2+RUNO-MINMLD*RHOW,0d0))
          ENRGFO=ACEFO*(TFL*SHI-LHM)
          E2O=FHO-ENRGFO
        END IF
      END IF

      IF (ROICE.GT.0) THEN
C**** Calculate energy in mixed layer (under ice)
        FHI=ELAKE(1)+DH2+ENRGI-(MLAKE(1)+DM2+RUNI)*TFL*SHW
        IF (FHI.LT.0) THEN ! FLUXES COOL WATER TO FREEZING, FORM ICE
          ACEFI =FHI/(TFL*(SHI-SHW)-LHM)
          ACEFI =MIN(ACEFI,MAX(MLAKE(1)+DM2+RUNI-MINMLD*RHOW,0d0))
          ENRGFI=ACEFI*(TFL*SHI-LHM)
          E2I=FHI-ENRGFI
        END IF
      END IF

C**** Update first layer variables
      MLAKE(1)=MLAKE(1)+DM2+(1.-ROICE)*(RUNO -ACEFO)+ROICE*(RUNI-ACEFI)
      ELAKE(1)=ELAKE(1)+DH2+(1.-ROICE)*(ENRGO-ENRGFO)+
     *                                             ROICE*(ENRGI-ENRGFI)

      ACEF1=0. ; ACEF2=0. ; ENRGF1=0. ; ENRGF2=0.
C**** Take remaining energy and start to freeze second layer
      FH2= ELAKE(1)-MLAKE(1)*TFL*SHW
      IF (FH2.LT.0) THEN
        IF (MLAKE(2).gt.0) THEN
C**** FH2=-ACEF2*(TLK2-TFL)*SHW+ACEF2*LHM
          TLK2    =ELAKE(2)/(MLAKE(2)*SHW)
          ACEF2   =-FH2/(TLK2*SHW-TFL*SHI+LHM)
          ACEF2   =MIN(ACEF2,MLAKE(2))
          ENRGF2  =ACEF2*(TFL*SHI-LHM)
          ELAKE(1)=ELAKE(1)+ACEF2*TLK2*SHW-ENRGF2
          ELAKE(2)=ELAKE(2)-ACEF2*TLK2*SHW
          MLAKE(2)=MLAKE(2)-ACEF2
        END IF
        FH1= ELAKE(1)-MLAKE(1)*TFL*SHW
        IF (FH1.lt.0) THEN      ! all layer 2 froze, freeze layer 1
          ACEF1   =FH1/(TFL*(SHI-SHW)-LHM)
C**** limit freezing if lake is between 50 and 20cm depth
          IF (MLAKE(1).lt.0.5d0*RHOW)
     *         ACEF1=MIN(ACEF1,MAX(0.5*(MLAKE(1)-0.2d0*RHOW),0d0))
          ENRGF1  =ACEF1*(TFL*SHI-LHM)
          ELAKE(1)=ELAKE(1)-ENRGF1
          MLAKE(1)=MLAKE(1)-ACEF1
          FH0     =ELAKE(1)-MLAKE(1)*TFL*SHW
          IF (FH0.lt.-1d-8) THEN ! max. amount of lake frozen, cool ice
            PRINT*,"Minimum lake level reached: rsi,mlake,elake",i0,j0
     *           ,roice,mlake(1)/rhow,elake(1)
            ENRGF1  =ENRGF1+FH0
            ELAKE(1)=MLAKE(1)*TFL*SHW
          END IF
        END IF
      END IF

C**** combine mass and energy fluxes for output
C**** Note that output fluxes are over open water/ice covered fractions
C**** distribute ice fluxes according to flux amounts
      FRATO = 1d0
      FRATI = 1d0
      IF (E2I+E2O.lt.0) THEN
        FRATO = E2O/(E2I*ROICE+E2O*(1.-ROICE))
        FRATI = E2I/(E2I*ROICE+E2O*(1.-ROICE))
      END IF
      ACEFO =ACEFO + (ACEF1 +ACEF2 )*FRATO
      ACEFI =ACEFI + (ACEF1 +ACEF2 )*FRATI
      ENRGFO=ENRGFO+ (ENRGF1+ENRGF2)*FRATO
      ENRGFI=ENRGFI+ (ENRGF1+ENRGF2)*FRATI

      RETURN
      END SUBROUTINE LKSOURC

      SUBROUTINE LKMIX(MLAKE,ELAKE,HLAKE,TKE,ROICE,DTSRC)
!@sum  LKMIX calculates mixing and entrainment in lakes
!@auth Gavin Schmidt
!@ver  1.0
      IMPLICIT NONE
!@var MLAKE,ELAKE mass and energy in lake layers (kg,J /m^2)
      REAL*8, INTENT(INOUT), DIMENSION(2) :: MLAKE,ELAKE
!@var TKE turbulent kinetic energy input at surface of lake (J/m^2)
!@var ROICE ice fraction
      REAL*8, INTENT(IN) :: TKE,ROICE
!@var HLAKE sill depth for lake (m)
      REAL*8, INTENT(IN) :: HLAKE
!@var DTSRC source time step (s)
      REAL*8, INTENT(IN) :: DTSRC
!@param MAXRHO,RHO0,BFAC freshwater density function approximation
      REAL*8, PARAMETER :: MAXRHO=1d3, RHO0=999.842594d0,
     *     BFAC=(MAXRHO-RHO0)/16d0

      REAL*8 TLK1, TLK2, HLT, MLT, DTK, E1N, E2N, ATKE, H1, H2,
     *      DRHO, DML, DHML

C**** Only mix if there is a second layer!
      IF (MLAKE(2).gt.0) THEN
        TLK1=ELAKE(1)/(MLAKE(1)*SHW)
        TLK2=ELAKE(2)/(MLAKE(2)*SHW)
        HLT=ELAKE(1)+ELAKE(2)
        MLT=MLAKE(1)+MLAKE(2)
C**** Test for static stability
        IF ((TMAXRHO-TLK1)*(TLK2-TLK1).lt.0) THEN
C**** mix uniformly and set MLD to minimum
          MLAKE(1)=MAX(MINMLD*RHOW,(MLAKE(1)+MLAKE(2))-HLAKE*RHOW)
          MLAKE(2)=MLT-MLAKE(1)
          ELAKE(1)=HLT*MLAKE(1)/MLT
          ELAKE(2)=HLT*MLAKE(2)/MLT
        ELSE ! not unstable, implicitly diffuse heat + entrain
C**** reduce mixing if there is ice cover, no mixing if temperature
C**** gradient is negative and deep water is at max density.
          IF (TLK2.lt.TLK1 .and. TLK2.gt.TMAXRHO) THEN
            DTK=2.*KVLAKE*(1.-ROICE)*DTSRC*RHOW**2
            E1N=(ELAKE(1)+DTK*HLT/(MLT*MLAKE(2)))/
     *           (1.+DTK/(MLAKE(1)*MLAKE(2)))
            E2N=(ELAKE(2)+DTK*HLT/(MLT*MLAKE(1)))/
     *           (1.+DTK/(MLAKE(1)*MLAKE(2)))
            ELAKE(1)=E1N
            ELAKE(2)=E2N
          END IF
C**** entrain deep water if there is available TKE
C**** take a factor of TKE and calculate change in PE
          IF (TKE.gt.0) THEN
            ATKE=0.2d0*TKE      ! 20% of TKE is available for mixing
            H1=MLAKE(1)/RHOW
            H2=MLAKE(2)/RHOW
C**** DRHO=RHO(TLK2)-RHO(TLK1)~=(TLK2-TLK1)*dRHOdT(TLK1)
C**** Assumes a parabolic density function going through MAXRHO at
C**** TMAXRHO, and RHO0 at T=0. (reasonable up to about 12 C)
            DRHO=(TLK2-TLK1)*2d0*BFAC*(TMAXRHO-TLK1)
            DML=ATKE*BYGRAV/(DRHO*0.5*H1)
            IF (DML*RHOW.lt.MLAKE(2)) THEN
              DHML=DML*ELAKE(2)/H2
              ELAKE(1)=ELAKE(1)+DHML
              ELAKE(2)=ELAKE(2)-DHML
              MLAKE(1)=MLAKE(1)+DML*RHOW
              MLAKE(2)=MLAKE(2)-DML*RHOW
            ELSE                ! entire second layer is entrained
              MLAKE(1)=MLT
              MLAKE(2)=0
              ELAKE(1)=HLT
              ELAKE(2)=0
            END IF
          END IF
        END IF
      END IF
      RETURN
C****
      END SUBROUTINE LKMIX

      END MODULE LAKES

      SUBROUTINE init_LAKES(inilake)
!@sum  init_LAKES initiallises lake variables
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : rhow,shw,tf,pi,grav
      USE MODEL_COM, only : im,jm,flake0,zatmo,dtsrc,flice,hlake,ftype
     *     ,itlake,itlkice
      USE GEOM, only : dxyp,dxv,dyv,dxp,imaxj
      USE SEAICE_COM, only : rsi
      USE PBLCOM, only : tsavg
      USE LAKES
      USE LAKES_COM
      USE DAGCOM, only : npts,icon_LKM,icon_LKE,tsfrez,tf_lkon,tf_lkoff
     *     ,title_con
      USE FLUXES, only : gtemp
      USE FILEMANAGER

      IMPLICIT NONE
      LOGICAL inilake
!@var I,J,I72,IU,JU,ID,JD,IMAX loop variables
      INTEGER I,J,I72,IU,JU,ID,JD,IMAX,INM
      INTEGER iu_RVR  !@var iu_RVR unit number for river direction file
      CHARACTER TITLEI*80, CDIREC(IM,JM)*1
      REAL*8 SPMIN,SPMAX,SPEED0,SPEED,DZDH,DZDH1,MLK1
      LOGICAL :: QCON(NPTS), T=.TRUE. , F=.FALSE.
C****
C**** LAKECB  MWL      Mass of water in lake (kg)
C****         GML      Liquid lake enthalpy (J)
C****         TLAKE    Temperature of lake surface (C)
C****         HLAKE    Lake sill depth (m)
C****         ALAKE    Lake surface area (m^2)
C****         TANLK    Lake slope (tan(alpha)) (1)
C****
C**** FIXDCB  FLAKE0    Original Lake fraction (1)
C****
      FLAKE = FLAKE0

      IF (INILAKE) THEN
C**** Set lake variables from surface temperature
C**** This is just an estimate for the initiallisation
        DO J=1,JM
          DO I=1,IM
c            FLAKE(I,J) = FLAKE0(I,J)
            IF (FLAKE(I,J).gt.0) THEN
              TLAKE(I,J) = MAX(0d0,TSAVG(I,J)-TF)
              MWL(I,J)   = RHOW*HLAKE(I,J)*FLAKE(I,J)*DXYP(J)
              MLK1       = MINMLD*RHOW*FLAKE(I,J)*DXYP(J)
              GML(I,J)   = SHW*(MLK1*TLAKE(I,J)
     *             +(MWL(I,J)-MLK1)*MAX(TLAKE(I,J),4d0))
              MLDLK(I,J) = MINMLD
            ELSE
              TLAKE(I,J) = 0.
              MWL(I,J)   = 0.
              GML(I,J)   = 0.
              MLDLK(I,J) = MINMLD
            END IF
            TSFREZ(I,J,3:4) = -999.
          END DO
        END DO
      END IF

C**** Set geometric variables
C**** TANLK=TAN(ALPHA) = R/H for a conical lake of equivalent volume
      DO J=1,JM
        DO I=1,IM
          IF (FLAKE0(I,J).gt.0) THEN
            TANLK(I,J) = SQRT(FLAKE0(I,J)*DXYP(J)/PI)/(3d0*HLAKE(I,J))
          ELSE
            TANLK(I,J) = 2d3   ! reasonable average value
          END IF
        END DO
      END DO

      CALL PRINTLK("IN")
C**** Set FTYPE array for lakes
      DO J=1,JM
        DO I=1,IM
          IF (FLAKE(I,J).gt.0) THEN
            FTYPE(ITLKICE,I,J)=FLAKE(I,J)*RSI(I,J)
            FTYPE(ITLAKE ,I,J)=FLAKE(I,J)-FTYPE(ITLKICE,I,J)
            GTEMP(1,1,I,J)=TLAKE(I,J)
            GTEMP(2,1,I,J)=(GML(I,J)-TLAKE(I,J)*SHW*FLAKE(I,J)*DXYP(J))
     *           /((MWL(I,J)-MLDLK(I,J)*RHOW*FLAKE(I,J)*DXYP(J))*SHW)
          END IF
        END DO
      END DO

C****
C**** Always initiallise River direction and Rate
C**** Read in CDIREC: Number = octant direction, Letter = river mouth
      call getunit("RVR",iu_RVR,.false.,.true.)
      READ  (iu_RVR,910) TITLEI
      WRITE (6,*) 'River Direction file read: ',TITLEI
      READ  (iu_RVR,910)
      DO I72=72,IM,72
        DO J=JM,1,-1
          READ  (iu_RVR,911) (CDIREC(I,J),I=I72-71,I72)
        END DO
      END DO
      CLOSE (iu_RVR)
C**** Create integral direction array KDIREC from CDIREC
      INM=0
      DO J=1,JM
      DO I=1,IM
        KDIREC(I,J) = ICHAR(CDIREC(I,J)) - 48
        IF(KDIREC(I,J).lt.0 .or. KDIREC(I,J).gt.8)  KDIREC(I,J) = 0
C**** Check for specified river mouths
        IF (ICHAR(CDIREC(I,J)).GE.65 .AND. ICHAR(CDIREC(I,J)).LE.90)
     *       THEN
          INM=INM+1
          IRVRMTH(INM)=I
          JRVRMTH(INM)=J
        END IF
      END DO
      END DO
C****
C**** From each box calculate the downstream river box
C****
      DO J=2,JM-1
        DO I=1,IM
          SELECT CASE (KDIREC(I,J))
          CASE (0)
            IFLOW(I,J) = I
            JFLOW(I,J) = J
            DHORZ(I,J) = 0.
          CASE (1)
            IFLOW(I,J) = I+1
            JFLOW(I,J) = J+1
            DHORZ(I,J) = SQRT(DXV(J)*DXV(J)+DYV(J)*DYV(J))
            IF(I.eq.IM)  IFLOW(I,J) = 1
          CASE (2)
            IFLOW(I,J) = I
            JFLOW(I,J) = J+1
            DHORZ(I,J) = DYV(J)
          CASE (3)
            IFLOW(I,J) = I-1
            JFLOW(I,J) = J+1
            DHORZ(I,J) = SQRT(DXV(J)*DXV(J)+DYV(J)*DYV(J))
            IF(I.eq.1)  IFLOW(I,J) = IM
          CASE (4)
            IFLOW(I,J) = I-1
            JFLOW(I,J) = J
            DHORZ(I,J) = DXP(J)
            IF(I.eq.1)  IFLOW(I,J) = IM
          CASE (5)
            IFLOW(I,J) = I-1
            JFLOW(I,J) = J-1
            DHORZ(I,J) = SQRT(DXV(J-1)*DXV(J-1)+DYV(J-1)*DYV(J-1))
            IF(I.eq.1)  IFLOW(I,J) = IM
          CASE (6)
            IFLOW(I,J) = I
            JFLOW(I,J) = J-1
            DHORZ(I,J) = DYV(J-1)
          CASE (7)
            IFLOW(I,J) = I+1
            JFLOW(I,J) = J-1
            DHORZ(I,J) = SQRT(DXV(J-1)*DXV(J-1)+DYV(J-1)*DYV(J-1))
            IF(I.eq.IM)  IFLOW(I,J) = 1
          CASE (8)
            IFLOW(I,J) = I+1
            JFLOW(I,J) = J
            DHORZ(I,J) = DXP(J)
            IF(I.eq.IM)  IFLOW(I,J) = 1
          END SELECT
        END DO
      END DO
C**** South Pole is a special case
      DO I=1,IM
        IF(KDIREC(I,1).eq.2)  THEN
          IDPOLE     = I
          JDPOLE     = 2
          IFLOW(1,1) = I
          JFLOW(1,1) = 2
          DHORZ(1,1) = DYV(1)
        END IF
        IF(KDIREC(I,2).eq.6)  THEN
          IFLOW(I,2) = 1
          JFLOW(I,2) = 1
        END IF
      END DO
C****
C**** Calculate river flow RATE (per source time step)
C****
      SPEED0= .35d0
      SPMIN = .15d0
      SPMAX = 5.
      DZDH1 = .00005
      DO JU=1,JM-1
        IMAX=IMAXJ(JU)
        DO IU=1,IMAX
          IF(KDIREC(IU,JU).gt.0) THEN
            JD=JFLOW(IU,JU)
            ID=IFLOW(IU,JU)
            DZDH  = (ZATMO(IU,JU)-ZATMO(ID,JD)) / (GRAV*DHORZ(IU,JU))
            SPEED = SPEED0*DZDH/DZDH1
            IF(SPEED.lt.SPMIN)  SPEED = SPMIN
            IF(SPEED.gt.SPMAX)  SPEED = SPMAX
            RATE(IU,JU) = DTsrc*SPEED/DHORZ(IU,JU)
          END IF
        END DO
      END DO

C**** Set conservation diagnostics for Lake mass and energy
      QCON=(/ F, F, F, T, T, T, F, F, T, F, F/)
      CALL SET_CON(QCON,"LAK MASS","(10**10 KG)     ",
     *     "(10**3 KG/S)    ",1d-10,1d-3,icon_LKM)
      QCON=(/ F, F, F, T, T, T, F, F, T, F, F/)
      CALL SET_CON(QCON,"LAK ENRG","(10**14 J)      ",
     *     "(10**8 J/S)     ",1d-14,1d-8,icon_LKE)

      RETURN
C****
 910  FORMAT (A72)
 911  FORMAT (72A1)
      END SUBROUTINE init_LAKES

      SUBROUTINE RIVERF
!@sum  RIVERF transports lake water from each grid box downstream
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0 (based on LB265)
      USE CONSTANT, only : grav,shw,rhow
      USE MODEL_COM, only : im,jm,focean,zatmo,hlake,itlake,itlkice
     *     ,ftype
      USE GEOM, only : dxyp
      USE LAKES, only : kdirec,idpole,jdpole,rate,iflow,jflow
      USE LAKES_COM, only : tlake,gml,mwl,mldlk,flake
      USE SEAICE_COM, only : rsi
      USE FLUXES, only : flowo,eflowo,gtemp
      USE DAGCOM, only : aij,ij_ervr,ij_mrvr

      IMPLICIT NONE
!@var I,J,IU,JU,ID,JD loop variables
      INTEGER I,J,IU,JU,ID,JD
      REAL*8 MWLSILL,DMM,DGM,HLK1
      REAL*8, DIMENSION(IM,JM) :: FLOW,EFLOW
C****
C**** LAKECB  MWL  Liquid lake mass  (kg)
C****         GML  Liquid lake enthalpy  (J)
C****         TLAKE  Lake surface temperature (C)
C****
C**** Calculate net mass and energy changes due to river flow
C****
      FLOW = 0. ; EFLOW = 0.
      FLOWO = 0. ; EFLOWO = 0.
      DO JU=2,JM-1
        DO IU=1,IM
          IF(KDIREC(IU,JU).gt.0)  THEN
            JD=JFLOW(IU,JU)
            ID=IFLOW(IU,JU)
C**** Only overflow if lake mass is above sill height (HLAKE (m))
            MWLSILL = RHOW*HLAKE(IU,JU)*FLAKE(IU,JU)*DXYP(JU)
            IF(MWL(IU,JU).gt.MWLSILL) THEN
              DMM = (MWL(IU,JU)-MWLSILL)*RATE(IU,JU)
c              IF (FLAKE(IU,JU).gt.0) THEN
c                MLM=RHOW*MLDLK(IU,JU)*FLAKE(IU,JU)*DXYP(JU)
c                DMM=MIN(DMM,MLM)   ! not necessary since MLM>TOTD-HLAKE
c              END IF
              DGM=TLAKE(IU,JU)*DMM*SHW  ! TLAKE always defined
              FLOW(IU,JU) =  FLOW(IU,JU) - DMM
              EFLOW(IU,JU) = EFLOW(IU,JU) - DGM
              AIJ(ID,JD,IJ_MRVR)=AIJ(ID,JD,IJ_MRVR) +  DMM
              AIJ(ID,JD,IJ_ERVR)=AIJ(ID,JD,IJ_ERVR)+
     *             DGM+DMM*(ZATMO(IU,JU)-ZATMO(ID,JD))
              IF(FOCEAN(ID,JD).le.0.) THEN
                FLOW(ID,JD) =  FLOW(ID,JD) + DMM
                EFLOW(ID,JD) = EFLOW(ID,JD) + DGM
     *               +DMM*(ZATMO(IU,JU)-ZATMO(ID,JD))
              ELSE ! Save river mouth flow to for output to oceans
                FLOWO(ID,JD)=FLOWO(ID,JD)+DMM
                EFLOWO(ID,JD)=EFLOWO(ID,JD)+DGM
     *               +DMM*(ZATMO(IU,JU)-ZATMO(ID,JD))
              END IF
            END IF
          END IF
        END DO
      END DO
C****
C**** Calculate river flow at the South Pole
C****
       FLOW(1,1) =  FLOW(1,1)/IM
      EFLOW(1,1) = EFLOW(1,1)/IM
C**** Only overflow if lake mass is above sill height (HLAKE (m))
      MWLSILL = RHOW*HLAKE(1,1)*FLAKE(1,1)*DXYP(1)
      IF(MWL(1,1).gt.MWLSILL) THEN
        DMM = (MWL(1,1)-MWLSILL)*RATE(1,1)
c        IF (FLAKE(1,1).gt.0) THEN
c         MLM=RHOW*MLDLK(1,1)*FLAKE(1,1)*DXYP(1)
c         DMM=MIN(DMM,MLM)   ! not necessary since MLM>TOTD-HLAKE
c        END IF
        DGM=TLAKE(1,1)*DMM*SHW ! TLAKE always defined
        FLOW(1,1) =  FLOW(1,1) - DMM
        EFLOW(1,1) = EFLOW(1,1) - DGM
        AIJ(IDPOLE,JDPOLE,IJ_MRVR)=AIJ(IDPOLE,JDPOLE,IJ_MRVR) + DMM
        AIJ(IDPOLE,JDPOLE,IJ_ERVR)=AIJ(IDPOLE,JDPOLE,IJ_ERVR) + DGM
     *       +DMM*(ZATMO(1,1)-ZATMO(IDPOLE,JDPOLE))
         FLOW(IDPOLE,JDPOLE) =  FLOW(IDPOLE,JDPOLE) + IM*DMM
        EFLOW(IDPOLE,JDPOLE) = EFLOW(IDPOLE,JDPOLE) +
     +       IM*(DGM+DMM*(ZATMO(1,1)-ZATMO(IDPOLE,JDPOLE)))
      END IF
C****
C**** Apply net river flow to continental reservoirs
C****
      DO J=2,JM-1
        DO I=1,IM
          IF(FOCEAN(I,J).le.0.) THEN
            MWL(I,J) = MWL(I,J) +  FLOW(I,J)
            GML(I,J) = GML(I,J) + EFLOW(I,J)
            IF (FLAKE(I,J).gt.0) THEN
              HLK1=(MLDLK(I,J)*RHOW)*TLAKE(I,J)*SHW
              MLDLK(I,J)=MLDLK(I,J)+FLOW(I,J)/(RHOW*FLAKE(I,J)*DXYP(J))
              TLAKE(I,J)=(HLK1*FLAKE(I,J)*DXYP(J)+EFLOW(I,J))
     *             /(MLDLK(I,J)*RHOW*FLAKE(I,J)*DXYP(J)*SHW)
            ELSE
              TLAKE(I,J)=GML(I,J)/(SHW*MWL(I,J)+1d-20)
            END IF
          END IF
        END DO
      END DO
      MWL(1,1) = MWL(1,1) +  FLOW(1,1)
      GML(1,1) = GML(1,1) + EFLOW(1,1)
      IF (FLAKE(1,1).gt.0) THEN
        HLK1=(MLDLK(1,1)*RHOW)*TLAKE(1,1)*SHW
        MLDLK(1,1)=MLDLK(1,1)+FLOW(1,1)/(RHOW*FLAKE(1,1)*DXYP(1))
        TLAKE(1,1)=(HLK1*FLAKE(1,1)*DXYP(1)+EFLOW(1,1))
     *       /(MLDLK(1,1)*RHOW*FLAKE(1,1)*DXYP(1)*SHW)
      ELSE
        TLAKE(1,1)=GML(1,1)/(SHW*MWL(1,1)+1d-20)
      END IF

      CALL PRINTLK("RV")
C**** Set FTYPE array for lakes
      DO J=1,JM
        DO I=1,IM
          IF (FLAKE(I,J).gt.0) THEN
            FTYPE(ITLKICE,I,J)=FLAKE(I,J)*RSI(I,J)
            FTYPE(ITLAKE ,I,J)=FLAKE(I,J)-FTYPE(ITLKICE,I,J)
            GTEMP(1,1,I,J)=TLAKE(I,J)
c            GTEMP(2,1,I,J)=(GML(I,J)-TLAKE(I,J)*SHW*FLAKE(I,J)*DXYP(J))
c     *           /((MWL(I,J)-MLDLK(I,J)*RHOW*FLAKE(I,J)*DXYP(J))*SHW)
          END IF
        END DO
      END DO

      RETURN
C****
      END SUBROUTINE RIVERF

      SUBROUTINE diag_RIVER
!@sum  diag_RIVER prints out the river outflow for various rivers
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : rhow,sday
      USE MODEL_COM, only : jyear0,amon0,jdate0,jhour0,jyear,amon
     *     ,jdate,jhour,itime,dtsrc,idacc,itime0,nday,jdpery,jmpery
      USE DAGCOM, only : aij,ij_mrvr
      USE LAKES, only : irvrmth,jrvrmth,nrvr
      IMPLICIT NONE
!@var NAMERVR Names of specified rivers
      CHARACTER*8, PARAMETER :: NAMERVR(NRVR) = (/
     *     "Murray  ","Parana  ","Orange  ","Limpopo ","Zambezi ",
     *     "SaoFranc","Congo   ","Amazon  ","Niger   ","Orinoco ",
     *     "Mekon   ","Magdalen","Godavari","Irrawadi","Brahmapu",
     *     "Indus   ","HsiChian","RioGrand","Tigris  ","Colorado",
     *     "Mississi","Yangtze ","Nile    ","Yellow  ","Volga   ",
     *     "Columbia","Loire   ","Danube  ","Fraser  ","StLawren",
     *     "Rhine   ","Amur    ","Nelson  ","Elbe    ","Yukon   ",
     *     "Severnay","Mackenzi","Kolym   ","Ob      ","Yenesei ",
     *     "Lena    ","Indigirk"/)
      REAL*8 RVROUT(6), SCALERVR, DAYS
      INTEGER INM,I

      DAYS=(Itime-Itime0)/DFLOAT(nday)
      WRITE(6,900) JYEAR0,AMON0,JDATE0,JHOUR0,JYEAR,AMON,JDATE,JHOUR
     *     ,ITIME,DAYS
C**** convert kg/(source time step) to km^3/mon
      SCALERVR = 1d-9*SDAY*JDPERY/(JMPERY*RHOW*DTSRC)
      DO INM=1,NRVR,6
        DO I=1,6
          RVROUT(I) = SCALERVR*AIJ(IRVRMTH(I-1+INM),JRVRMTH(I-1+INM)
     *         ,IJ_MRVR)/IDACC(1)
        END DO
        WRITE(6,901) (NAMERVR(I-1+INM),RVROUT(I),I=1,6)
      END DO

      RETURN
C****
 900  FORMAT ('0* River Outflow (km^3/mon) **  From:',I6,A6,I2,',  Hr'
     *     ,I3,6X,'To:',I6,A6,I2,', Hr',I3,'  Model-Time:',I9,5X
     *     ,'Dif:',F7.2,' Days')
 901  FORMAT (' ',A8,':',F8.3,5X,A8,':',F8.3,5X,A8,':',F8.3,5X,
     *            A8,':',F8.3,5X,A8,':',F8.3,5X,A8,':',F8.3)
      END SUBROUTINE diag_RIVER

      SUBROUTINE CHECKL (SUBR)
!@sum  CHECKL checks whether the lake variables are reasonable.
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0 (based on LB265)
      USE CONSTANT, only : rhow
      USE MODEL_COM, only : im,jm,hlake,fearth
      USE GEOM, only : dxyp
      USE DAGCOM, only : qcheck
      USE LAKES
      USE LAKES_COM

      IMPLICIT NONE

      INTEGER I,J !@var I,J loop variables
      CHARACTER*6, INTENT(IN) :: SUBR
      LOGICAL QCHECKL

C**** Check for NaN/INF in lake data
      CALL CHECK3(MWL,IM,JM,1,SUBR,'mwl')
      CALL CHECK3(GML,IM,JM,1,SUBR,'gml')
      CALL CHECK3(MLDLK,IM,JM,1,SUBR,'mld')
      CALL CHECK3(TLAKE,IM,JM,1,SUBR,'tlk')

      QCHECKL = .FALSE.
      DO J=2,JM-1
      DO I=1,IM
        IF(FEARTH(I,J).gt.0.) THEN
C**** check for negative mass
          IF (MWL(I,J).lt.0 .or. MLDLK(I,J).lt.0) THEN
            WRITE(6,*) 'After ',SUBR,': I,J,TSL,MWL,GML,MLD=',I,J
     *           ,TLAKE(I,J),MWL(I,J),GML(I,J),MLDLK(I,J)
            QCHECKL = .TRUE.
          END IF
C**** check for reasonable lake surface temps
          IF (TLAKE(I,J).ge.40 .or. TLAKE(I,J).lt.-0.5) THEN
            WRITE(6,*) 'After ',SUBR,': I,J,TSL=',I,J,TLAKE(I,J)
c            QCHECKL = .TRUE.
          END IF
        END IF
C**** Check total lake mass (<10%, >10x orig depth)
        IF(FLAKE(I,J).gt.0.) THEN
          IF(MWL(I,J).lt.0.1d0*RHOW*HLAKE(I,J)*DXYP(J)*FLAKE(I,J)) THEN
            WRITE (6,*) 'After ',SUBR,
     *           ': I,J,FLAKE,HLAKE,lake level low=',I,J,FLAKE(I,J),
     *           HLAKE(I,J),MWL(I,J)/(RHOW*DXYP(J)*FLAKE(I,J))
          END IF
          IF(MWL(I,J).gt.RHOW*MAX(10.*HLAKE(I,J),3d1)*DXYP(J)*FLAKE(I,J)
     *         )THEN
            WRITE (6,*) 'After ',SUBR,
     *           ': I,J,FLAKE,HLAKE,lake level high=',I,J,FLAKE(I,J),
     *           HLAKE(I,J),MWL(I,J)/(RHOW*DXYP(J)*FLAKE(I,J))
          END IF
        END IF
      END DO
      END DO
      IF (QCHECKL) STOP 'CHECKL: Lake variables out of bounds'
      RETURN
C****
      END SUBROUTINE CHECKL

      SUBROUTINE daily_LAKE(IEND)
!@sum  daily_LAKE does lake things at the end of every day
!@auth G. Schmidt
!@ver  1.0
      USE CONSTANT, only : shw,rhow,pi,by3
      USE MODEL_COM, only : im,jm,ftype,itlake,itlkice,jday
      USE LAKES_COM, only : tlake,mwl,gml,mldlk,flake,tanlk
      USE SEAICE_COM, only : rsi,msi,hsi,snowi,ssi
      USE SEAICE, only : simelt,lmi
      USE GEOM, only : imaxj,dxyp
      USE FLUXES, only : gtemp
      USE DAGCOM, only : tsfrez,tf_lkon,tf_lkoff,aij,ij_lkon,ij_lkoff
      IMPLICIT NONE
      INTEGER IEND,IMAX,I,J,L
!@var FDAILY fraction of energy available to be used for melting
      REAL*8 :: FDAILY = BY3
      REAL*8, DIMENSION(LMI) :: HSIL,TSIL,SSIL
      REAL*8 MSI2,ROICE,SNOW,ENRGW,ENRGUSED,ANGLE,RUN0,SALT,POCEAN,TFO

C**** set and initiallise freezing diagnostics
C**** Note that TSFREZ saves the last day of no-ice and some-ice.
C**** The AIJ diagnostics are set once a year (zero otherwise)
C**** To ensure correct averaging they are multiplied by the number of
C**** months in a year.
      DO J=1,JM
        DO I=1,IMAXJ(J)
          IF (J.le.JM/2) THEN
C**** initiallise/save South. Hemi. on Jan 1
            IF (JDAY.eq.1 .and. TSFREZ(I,J,TF_LKOFF).ne.-999) THEN
              AIJ(I,J,IJ_LKON)  = 12.*TSFREZ(I,J,TF_LKON)
              AIJ(I,J,IJ_LKOFF) = 12.*TSFREZ(I,J,TF_LKOFF)
              TSFREZ(I,J,TF_LKOFF) = -999.
            END IF
          ELSE
C**** initiallise/save North. Hemi. on Jul 1
            IF (JDAY.eq.182 .and. TSFREZ(I,J,TF_LKOFF).ne.-999) THEN
              AIJ(I,J,IJ_LKON)  = 12.*TSFREZ(I,J,TF_LKON)
              AIJ(I,J,IJ_LKOFF) = 12.*TSFREZ(I,J,TF_LKOFF)
              TSFREZ(I,J,TF_LKOFF) = -999.
            END IF
          END IF
C**** set ice on/off days
          IF (FLAKE(I,J).gt.0) THEN
            IF (RSI(I,J).eq.0.and.TSFREZ(I,J,TF_LKOFF).eq.-999.)
     *           TSFREZ(I,J,TF_LKON)=JDAY
            IF (RSI(I,J).gt.0) TSFREZ(I,J,TF_LKOFF)=JDAY
          END IF
        END DO
      END DO

C**** Melt lake ice if energy is available in mixed layer
      DO J=1,JM
        IMAX=IMAXJ(J)
        DO I=1,IMAX
          IF (FLAKE(I,J)*RSI(I,J) .GT. 0.) THEN
C**** REDUCE ICE EXTENT IF LAKE TEMPERATURE IS GREATER THAN ZERO
C**** (MELTING POINT OF ICE)
C**** Also remove ice fractions less than 0.0001
            IF (TLAKE(I,J).GT.0. .or. RSI(I,J).lt.1d-4) THEN
              ROICE=RSI(I,J)
              MSI2 =MSI(I,J)
              SNOW =SNOWI(I,J)   ! snow mass
              HSIL =HSI(:,I,J) ! sea ice enthalpy
              SSIL =0.         ! sea ice salt (always 0)
              POCEAN=0.        ! ocean fraction (always 0)
              TFO = 0.         ! freezing point (always 0) 
C**** energy of water available for melting
              ENRGW=TLAKE(I,J)*MLDLK(I,J)*SHW*RHOW*FDAILY
              CALL SIMELT(ROICE,SNOW,MSI2,HSIL,SSIL,POCEAN,TFO,TSIL
     *             ,ENRGW,ENRGUSED,RUN0,SALT)
C**** SALT always 0 for lakes
C**** RESAVE PROGNOSTIC QUANTITIES
              GML(I,J)=GML(I,J)-FLAKE(I,J)*DXYP(J)*ENRGUSED
              MWL(I,J)=MWL(I,J)+FLAKE(I,J)*DXYP(J)*RUN0
              MLDLK(I,J)=MLDLK(I,J)+RUN0/RHOW
              TLAKE(I,J)=TLAKE(I,J)-ENRGUSED/(MLDLK(I,J)*SHW*RHOW)
              RSI(I,J)=ROICE
              MSI(I,J)=MSI2
              SNOWI(I,J)=SNOW
              HSI(:,I,J)=HSIL(:)
              SSI(:,I,J)=0.
C**** set ftype/gtemp arrays
              FTYPE(ITLKICE,I,J)=FLAKE(I,J)*RSI(I,J)
              FTYPE(ITLAKE ,I,J)=FLAKE(I,J)-FTYPE(ITLKICE,I,J)
              GTEMP(1:2,2,I,J) = TSIL(1:2)
              GTEMP(1  ,1,I,J) = TLAKE(I,J)
            END IF
          END IF
        END DO
      END DO

C**** Experimental code: not yet functional
C**** Update lake fraction as a function of lake mass at end of day
C**** Assume lake is conical
C****   => A = pi*(h*tanlk)^2, M=(1/3)*pi*rho*h*(h*tanlk)^2
C****
c      IF (IEND.gt.0) THEN
c        PRINT*,"TEST FLAKE CHANGE"
c        DO J=1,JM
c        DO I=1,IMAXJ(J)
c          IF (FLAKE(I,J).gt.0) THEN
c            PRINT*,"tan(a)",I,J,TANLK(I,J)
c            PRINT*,FLAKE(I,J),(9d0*PI*(TANLK(I,J)*MWL(I,J)/RHOW)**2)
c     *           **BY3/DXYP(J)
C**** remove lakes that are too small ??
c            IF (FLAKE(I,J).lt.0.005) THEN
c          PRINT*,"Lake->0 (too small)",I,J,MWL(I,J),GML(I,J),MLDLK(I,J)
c          MLDLK(I,J)=MIN(MLDLK(I,J),MWL(I,J)/(RHOW*FLAKE(I,J)*DXYP(J)))
c              FLAKE(I,J) = 0.
c            FTYPE(ITLKICE,I,J)=FLAKE(I,J)*RSI(I,J)
c            FTYPE(ITLAKE ,I,J)=FLAKE(I,J)-FTYPE(ITLKICE,I,J)
c          END IF
c        END DO
c        END DO
c      END IF
C****
      CALL PRINTLK("DY")
C****
      RETURN
      END SUBROUTINE daily_LAKE

      SUBROUTINE PRECIP_LK
!@sum  PRECIP_LK driver for applying precipitation/melt to lake fraction
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : rhow,shw
      USE MODEL_COM, only : im,jm,fland,flice,itlake,itlkice
      USE GEOM, only : imaxj,dxyp
      USE SEAICE_COM, only : rsi
      USE LAKES_COM, only : mwl,gml,tlake,mldlk,flake
      USE FLUXES, only : runpsi,runoli,prec,eprec,gtemp
      USE DAGCOM, only : aj,aij,ij_f0oc
      IMPLICIT NONE

      REAL*8 PRCP,ENRGP,DXYPJ,PLICE,PLKICE,RUN0,ERUN0,POLAKE,HLK1
      INTEGER I,J,IMAX

      CALL PRINTLK("PR")

      DO J=1,JM
      IMAX=IMAXJ(J)
      DXYPJ=DXYP(J)
      DO I=1,IMAX
      IF (FLAKE(I,J)+FLICE(I,J).gt.0) THEN
        POLAKE=(1.-RSI(I,J))*FLAKE(I,J)
        PLKICE=RSI(I,J)*FLAKE(I,J)
        PLICE=FLICE(I,J)
        PRCP=PREC(I,J)
        ENRGP=EPREC(I,J)        ! energy of precipitation

C**** save diagnostics
        IF (FLAKE(I,J).gt.0) THEN
          AIJ(I,J,IJ_F0OC)=AIJ(I,J,IJ_F0OC)+ENRGP*POLAKE
        END IF

C**** calculate fluxes over whole box
        RUN0 =POLAKE*PRCP  + PLKICE* RUNPSI(I,J) + PLICE* RUNOLI(I,J)
        ERUN0=POLAKE*ENRGP ! PLKICE*ERUNPSI(I,J) + PLICE*ERUNOLI(I,J) = 0

        MWL(I,J) = MWL(I,J) +  RUN0*DXYPJ
        GML(I,J) = GML(I,J) + ERUN0*DXYPJ

        IF (FLAKE(I,J).gt.0) THEN
          HLK1=TLAKE(I,J)*MLDLK(I,J)*RHOW*SHW
          MLDLK(I,J)=MLDLK(I,J) + RUN0/(FLAKE(I,J)*RHOW)
          TLAKE(I,J)=(HLK1*FLAKE(I,J)+ERUN0)/(MLDLK(I,J)*FLAKE(I,J)
     *         *RHOW*SHW)
          GTEMP(1,1,I,J)=TLAKE(I,J)
        ELSE
          TLAKE(I,J)=GML(I,J)/(MWL(I,J)*SHW+1d-20)
        END IF

      END IF
      END DO
      END DO
      RETURN
C****
      END SUBROUTINE PRECIP_LK

      SUBROUTINE GROUND_LK
!@sum  GROUND_LK driver for applying surface fluxes to lake fraction
!@auth Gavin Schmidt
!@ver  1.0
!@calls
      USE CONSTANT, only : rhow,shw
      USE MODEL_COM, only : im,jm,flice,fland,hlake
     *     ,fearth,dtsrc,itlake,itlkice
      USE GEOM, only : imaxj,dxyp
      USE FLUXES, only : runosi, erunosi, e0,evapor, dmsi,dhsi,dssi,
     *     runoli, runoe, erunoe, solar, dmua, dmva,gtemp
      USE SEAICE_COM, only : rsi
      USE PBLCOM, only : ustar
      USE DAGCOM, only : aj,aij,areg,jreg,ij_f0oc,j_run2
     *     ,j_dwtr2,j_tg1,j_tg2,j_evap,j_oht,j_omlt,j_erun2,j_imelt
     *     ,ij_tgo,ij_tg1,ij_evap,ij_evapo,j_type
      USE LAKES_COM, only : mwl,gml,tlake,mldlk,flake
      USE LAKES, only : lkmix,lksourc,byzeta
      IMPLICIT NONE
C**** grid box variables
      REAL*8 ROICE, POLAKE, PLKICE, PEARTH, PLICE, DXYPJ
!@var MLAKE,ELAKE mass and energy /m^2 for lake model layers
      REAL*8, DIMENSION(2) :: MLAKE,ELAKE
C**** fluxes
      REAL*8 EVAPO, F2DT, F0DT, RUN0, ERUN0, RUNLI, RUNE, ERUNE,
     *     HLK1,TLK1,TLK2,TKE,SROX(2),FSR2, U2RHO
C**** output from LKSOURC
      REAL*8 ENRGFO, ACEFO, ACEFI, ENRGFI

      INTEGER I,J,IMAX,JR

      CALL PRINTLK("GR")

      DO J=1,JM
      IMAX=IMAXJ(J)
      DXYPJ=DXYP(J)
      DO I=1,IMAX
      JR=JREG(I,J)
      ROICE=RSI(I,J)
      PLKICE=FLAKE(I,J)*ROICE
      POLAKE=FLAKE(I,J)*(1.-ROICE)

C**** Add land ice and surface runoff to lake variables
      IF (FLAND(I,J).gt.0) THEN
        PLICE =FLICE(I,J)
        PEARTH=FEARTH(I,J)
        RUNLI=RUNOLI(I,J)
        RUNE =RUNOE(I,J)
        ERUNE=ERUNOE(I,J)
C**** calculate flux over whole box
        RUN0 =RUNLI*PLICE + RUNE*PEARTH
        ERUN0=             ERUNE*PEARTH
        MWL(I,J) = MWL(I,J) + RUN0*DXYPJ
        GML(I,J) = GML(I,J) +ERUN0*DXYPJ
        IF (FLAKE(I,J).gt.0) THEN
          HLK1=TLAKE(I,J)*MLDLK(I,J)*RHOW*SHW
          MLDLK(I,J)=MLDLK(I,J) + RUN0/(FLAKE(I,J)*RHOW)
          TLAKE(I,J)=(HLK1*FLAKE(I,J)+ERUN0)/(MLDLK(I,J)*FLAKE(I,J)
     *         *RHOW*SHW)
        ELSE
          TLAKE(I,J)=GML(I,J)/(MWL(I,J)*SHW+1d-20)
        END IF
      END IF

      IF (FLAKE(I,J).gt.0) THEN
        TLK1 =TLAKE(I,J)
        EVAPO=EVAPOR(I,J,1)     ! evap/dew over open lake (kg/m^2)
        F0DT =E0(I,J,1)         ! net heat over open lake (J/m^2)
        SROX(1)=SOLAR(1,I,J)      ! solar radiation open lake (J/m^2)
        SROX(2)=SOLAR(3,I,J)      ! solar radiation through ice (J/m^2)
        FSR2 =EXP(-MLDLK(I,J)*BYZETA)
C**** get ice-ocean fluxes from sea ice routine (over ice fraction)
        RUN0 =RUNOSI(I,J)        ! includes ACE2M term
        F2DT =ERUNOSI(I,J)
C**** calculate kg/m^2, J/m^2 from saved variables
        MLAKE(1)=MLDLK(I,J)*RHOW
        MLAKE(2)=MAX(MWL(I,J)/(FLAKE(I,J)*DXYPJ)-MLAKE(1),0d0)
        ELAKE(1)=TLK1*SHW*MLAKE(1)
        ELAKE(2)=GML(I,J)/(FLAKE(I,J)*DXYPJ)-ELAKE(1)
        IF (MLAKE(2).lt.1d-10) THEN
          MLAKE(2)=0.
          ELAKE(2)=0.
        END IF
C**** Limit FSR2 in the case of thin second layer
        FSR2=MIN(FSR2,MLAKE(2)/(MLAKE(1)+MLAKE(2)))

        AJ(J,J_TG1, ITLAKE) =AJ(J,J_TG1, ITLAKE) +TLK1 *POLAKE
        AJ(J,J_EVAP,ITLAKE) =AJ(J,J_EVAP,ITLAKE) +EVAPO*POLAKE
        AJ(J,J_TYPE,ITLAKE) =AJ(J,J_TYPE,ITLAKE) +      POLAKE
        AJ(J,J_TYPE,ITLKICE)=AJ(J,J_TYPE,ITLKICE)+      PLKICE
        IF (JR.ne.24) AREG(JR,J_TG1)=AREG(JR,J_TG1)+TLK1*POLAKE*DXYPJ
        AIJ(I,J,IJ_TGO)  =AIJ(I,J,IJ_TGO)  +TLK1
        AIJ(I,J,IJ_TG1)  =AIJ(I,J,IJ_TG1)  +TLK1 *POLAKE
        AIJ(I,J,IJ_EVAP) =AIJ(I,J,IJ_EVAP) +EVAPO*POLAKE
        AIJ(I,J,IJ_EVAPO)=AIJ(I,J,IJ_EVAPO)+EVAPO*POLAKE

C**** Apply fluxes and calculate the amount of ice formation
        CALL LKSOURC (I,J,ROICE,MLAKE,ELAKE,RUN0,F0DT,F2DT,SROX,FSR2
     *       ,EVAPO,ENRGFO,ACEFO,ACEFI,ENRGFI)

C**** Mixing and entrainment
C**** Calculate turbulent kinetic energy for lake
        U2rho=(1.-ROICE)*SQRT(DMUA(I,J,1)**2+DMVA(I,J,1)**2)/DTSRC
c       TKE=0.5 * (19.3)^(2/3) * U2rho /rhoair ! (m/s)^2
        TKE=0.  ! 3.6d0*U2rho/rhoair*MLAKE(1)  ! (J/m^2)

        CALL LKMIX (MLAKE,ELAKE,HLAKE(I,J),TKE,ROICE,DTSRC)

C**** Resave prognostic variables
        MWL(I,J)  =(MLAKE(1)+MLAKE(2))*(FLAKE(I,J)*DXYPJ)
        GML(I,J)  =(ELAKE(1)+ELAKE(2))*(FLAKE(I,J)*DXYPJ)
        MLDLK(I,J)= MLAKE(1)/RHOW
        TLAKE(I,J)= ELAKE(1)/(SHW*MLAKE(1))
        IF (MLAKE(2).gt.0) THEN
          TLK2    = ELAKE(2)/(SHW*MLAKE(2))
        ELSE
          TLK2    = TLAKE(I,J)
        END IF
        GTEMP(1,1,I,J)=TLAKE(I,J)
        GTEMP(2,1,I,J)=TLK2       ! not used

C**** Open lake diagnostics
        AJ(J,J_TG2, ITLAKE)=AJ(J,J_TG2, ITLAKE)+TLK2*POLAKE
        IF (JR.ne.24) AREG(JR,J_TG2)=AREG(JR,J_TG2)+TLK2
     *       *POLAKE*DXYPJ
        AIJ(I,J,IJ_F0OC)=AIJ(I,J,IJ_F0OC)+F0DT*POLAKE
C**** Ice-covered ocean diagnostics
        AJ(J,J_ERUN2,ITLAKE) =AJ(J,J_ERUN2,ITLAKE) -ENRGFO*POLAKE
        AJ(J,J_IMELT,ITLAKE) =AJ(J,J_IMELT,ITLAKE) -ACEFO *POLAKE
        AJ(J,J_ERUN2,ITLKICE)=AJ(J,J_ERUN2,ITLKICE)-ENRGFI*PLKICE
        AJ(J,J_IMELT,ITLKICE)=AJ(J,J_IMELT,ITLKICE)-ACEFI *PLKICE

C**** Store mass and energy fluxes for formation of sea ice
        DMSI(1,I,J)=ACEFO
        DMSI(2,I,J)=ACEFI
        DHSI(1,I,J)=ENRGFO
        DHSI(2,I,J)=ENRGFI
        DSSI(:,I,J)=0.     ! assume zero salinity
      END IF
      END DO
      END DO

      CALL PRINTLK("G2")

      RETURN
C****
      END SUBROUTINE GROUND_LK


      SUBROUTINE PRINTLK(STR)
!@sum  PRINTLK print out selected diagnostics from specified lakes
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : lhm,byshi,rhow,shw
      USE GEOM, only : dxyp
      USE LAKES_COM, only : tlake,mwl,mldlk,gml,flake
      USE SEAICE_COM, only : rsi,hsi,msi,snowi
      USE SEAICE, only : xsi,ace1i,rhoi
      USE DAGCOM, only : qcheck
      IMPLICIT NONE
      CHARACTER*2, INTENT(IN) :: STR
      INTEGER, PARAMETER :: NDIAG=6
      INTEGER I,J,N,IDIAG(NDIAG),JDIAG(NDIAG)
      DATA IDIAG/57,53,15,18,22,23/,  !10/,
     *     JDIAG/31,32,37,33,25,22/  !,40/
      REAL*8 HLK2,TLK2

      IF (.NOT.QCHECK) RETURN

      DO N=1,NDIAG
        I=IDIAG(N)
        J=JDIAG(N)
        HLK2 = MWL(I,J)/(RHOW*FLAKE(I,J)*DXYP(J)) - MLDLK(I,J)
        IF (HLK2.gt.0) THEN
        TLK2 = (GML(I,J)/(SHW*RHOW*FLAKE(I,J)*DXYP(J)) -
     *       TLAKE(I,J)*MLDLK(I,J))/HLK2
        ELSE
          TLK2=0.
        END IF
        WRITE(99,*) STR,I,J,TLAKE(I,J),TLK2,MLDLK(I,J),HLK2
     *       ,RSI(I,J),MSI(I,J)/RHOI,SNOWI(I,J)/RHOW
      END DO

      RETURN
      END  SUBROUTINE PRINTLK

      SUBROUTINE conserv_LKM(LKM)
!@sum  conserv_LKM calculates total lake mass
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm
      USE GEOM, only : imaxj
      USE LAKES_COM, only : mwl
      IMPLICIT NONE
      REAL*8, DIMENSION(JM) :: LKM
      INTEGER :: I,J
C****
C**** LAKE MASS (kg)
C****
      DO J=1,JM
        LKM(J)=0.
        DO I=1,IMAXJ(J)
          LKM(J)=LKM(J)+MWL(I,J)
        END DO
      END DO
      RETURN
      END SUBROUTINE conserv_LKM

      SUBROUTINE conserv_LKE(LKE)
!@sum  conserv_LKE calculates total lake energy
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : grav
      USE MODEL_COM, only : im,jm,zatmo
      USE GEOM, only : imaxj
      USE LAKES_COM, only : gml,mwl
      IMPLICIT NONE
      REAL*8, DIMENSION(JM) :: LKE
      INTEGER :: I,J
C****
C**** LAKE ENERGY (J) (includes potential energy)
C****
      DO J=1,JM
        LKE(J)=0.
        DO I=1,IMAXJ(J)
          LKE(J)=LKE(J)+GML(I,J)+ZATMO(I,J)*MWL(I,J)
        END DO
      END DO
      RETURN
      END SUBROUTINE conserv_LKE
