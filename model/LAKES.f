      MODULE LAKES
!@sum  LAKES subroutines for Lakes and Rivers
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0 (based on LB265)
      USE E001M12_COM, only : IM,JM
      IMPLICIT NONE
      
C**** 
C**** Changes from Model III: MO -> MWL (kg), G0M -> GML (J), 
C****    GZM -> TLAKE (deg C),  RSI -> RLI (fract),
C****    MSI -> MLI (kg/m^2), HSI -> HLI (J/m^2)
C****
!@var LMLI Number of layers of lake ice model
      INTEGER, PARAMETER :: LMLI=4
!@var XLI1,XLI2,XLI3,XLI4 layering for lake ice model
      REAL*8, PARAMETER :: XLI1=.50,XLI2=.50,XLI3=.50,XLI4=.50
!@var BYXLI1,BYXLI2,BYXLI3,BYXLI4 reciprocals of XLIx
      REAL*8, PARAMETER :: BYXLI1=1./XLI1, BYXLI2=1./XLI2,
     *                     BYXLI3=1./XLI3, BYXLI4=1./XLI4
!@var KDIREC directions for river flow 
!**** (0 no flow, 1-8 anti-clockwise from top RH corner
!@var IDPOLE,JDPOLE special directions for south pole flow
      INTEGER KDIREC(IM,JM),IDPOLE,JDPOLE
!@var RATE rate of river flow downslope (m/s)
!@var DHORZ horizontal distance to downstream box (m)
      REAL*8, DIMENSION(IM,JM) :: RATE,DHORZ
!@var FLOWO,EFLOWO runoff and energy of runoff into ocean 
      REAL*8, DIMENSION(IM,JM) :: FLOWO,EFLOWO

      END MODULE LAKES

      SUBROUTINE init_LAKES(inilake)
!@sum  init_LAKES initiallises lake variables
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : shi,lhm
      USE E001M12_COM, only : im,jm,gdata,flake,zatmo,dt,flice,hlake
     *     ,kocean
      USE OCEAN, only : odata
      USE SEAICE, only : ace1i
      USE GEOM, only : dxyp,dxv,dyv,dxp,imaxj
      USE LAKES
      USE LAKES_COM
      USE FILEMANAGER
      
      IMPLICIT NONE
      LOGICAL inilake
      REAL*8 HLI1,HLI2,HLI3,HLI4
!@var I,J,I72,IU,JU,ID,JD,IMAX loop variables
      INTEGER I,J,I72,IU,JU,ID,JD,IMAX 
      INTEGER iu_RVR
      INTEGER*4 IFLOW(IM,JM),JFLOW(IM,JM)
      CHARACTER TITLEI*80, CDIREC(IM,JM)*1
      REAL*8 SPMIN,SPMAX,SPEED0,SPEED,DZDH,DZDH1
C**** 
C**** LAKECB  MWL      Mass of water in lake (kg)
C****         GML      Liquid lake enthalpy (J)
C****         TLAKE      Temperature of lake surface (C)
C****         RLI      Horizontal ratio of lake ice to lake (1)
C****         MLI      Mass of sea ice (kg/m**2)
C****         HLI      Heat including latent of sea ice (J/m**2)
C****         HLAKE    Lake sill depth (m)
C****
C**** FIXDCB  FLAKE    Lake fraction (1)  
C****

C**** Always set lake ice extent consistent with T50
c      IF (KOCEAN.eq.1) CALL set_LAKEICE

      IF (INILAKE) THEN
C**** Set lake variables from model blocks ODATA and GDATA
C**** GDATA(3) and (7) are first and second thermal layers
C**** This is just an estimate for the initiallisation
      DO J=2,JM-1
       DO I=1,IM
         IF (FLAKE(I,J).gt.0) THEN
           TLAKE(I,J)   = ODATA(I,J,1)
c           RLI(I,J)   = ODATA(I,J,2)
c           MLI(I,J,2) = ODATA(I,J,3)
c           MLI(I,J,1) = GDATA(I,J,1) + ACE1I
c           HLI1=(SHI*GDATA(I,J,3)-LHM)*MLI(I,J,1)
c           HLI(I,J,1) = XLI1*HLI1
c           HLI(I,J,2) = XLI2*HLI1
c           HLI2=(SHI*GDATA(I,J,7)-LHM)*MLI(I,J,2)
c           HLI(I,J,3) = XLI3*HLI2
c           HLI(I,J,4) = XLI4*HLI2
           MWL(I,J) = 1d3*HLAKE(I,J)*FLAKE(I,J)*DXYP(J)
           GML(I,J) = MWL(I,J)*MAX(TLAKE(I,J),4d0)
C**** reset GDATA correctly
c           GDATA(I,J,3) = (HLI(I,J,1)/(XLI1*MLI(I,J,1))+LHM)/SHI
c           GDATA(I,J,7) = (HLI(I,J,2)/(XLI2*MLI(I,J,1))+LHM)/SHI
         ELSE
           TLAKE(I,J)   = 0.
c           RLI(I,J)   = 0.
c           MLI(I,J,1:2) = 0.
c           HLI(I,J,1:4) = 0.
           MWL(I,J) = 0.
           GML(I,J) = 0.
         END IF
       END DO
      END DO
      END IF
C****
C**** Always initiallise River direction and Rate
C****
C**** Read in CDIREC: Number = octant direction, Letter = river mouth
      call getunit("RVR",iu_RVR,.FALSE.)
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
      DO I=1,IM*JM
      KDIREC(I,1) = ICHAR(CDIREC(I,1)) - 48
      IF(KDIREC(I,1).lt.0 .or. KDIREC(I,1).gt.8)  KDIREC(I,1) = 0
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
            DZDH  = (ZATMO(IU,JU)-ZATMO(ID,JD)) / DHORZ(IU,JU)
            SPEED = SPEED0*DZDH/DZDH1
            IF(SPEED.lt.SPMIN)  SPEED = SPMIN
            IF(SPEED.gt.SPMAX)  SPEED = SPMAX
            RATE(IU,JU) = DT*SPEED/DHORZ(IU,JU)
          END IF
        END DO
      END DO
C**** 
C**** Set runoff temperature of glacial ice to be 0 (C)
C****
      DO J=1,JM
        DO I=1,IM
          IF(FLICE(I,J).GT.0.)  TLAKE(I,J) = 0
        END DO
      END DO

      RETURN
C****
 910  FORMAT (A72)
 911  FORMAT (72A1)
      END SUBROUTINE init_LAKES


      SUBROUTINE RIVERF
!@sum  LAKES subroutines for Lakes and Rivers
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0 (based on LB265)
C****
C**** RIVERF transports lake water from each GCM grid box to its
C**** downstream neighbor according to the river direction file.
C****
      USE CONSTANT, only : shi,lhm,grav,shw
      USE E001M12_COM, only : im,jm,flake,dt,focean,zatmo,hlake
      USE GEOM, only : dxyp
      USE LAKES
      USE LAKES_COM
      USE DAGCOM, only : aij,ij_ervr,ij_mrvr
     
      IMPLICIT NONE
!@var I,J,IU,JU,ID,JD loop variables
      INTEGER I,J,IU,JU,ID,JD
      INTEGER iu_RVR  !@var iu_RVR unit number for river direction file
      INTEGER IFLOW(IM,JM),JFLOW(IM,JM)
c      REAL*8 XK(8),YK(8)
c      DATA XK/.707107d0,0.,-.707107d0,-1.,-.707107d0, 0., .707107d0,1./,
c     *     YK/.707107d0,1., .707107d0, 0.,-.707107d0,-1.,-.707107d0,0./
      REAL*8 MWLSILL,DMM,DGM
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
            MWLSILL = 1d3*HLAKE(IU,JU)*FLAKE(IU,JU)*DXYP(JU)
            IF(MWL(IU,JU).gt.MWLSILL) THEN
              DMM = (MWL(IU,JU)-MWLSILL)*RATE(IU,JU)
C             DGM = GML(IU,JU)*DMM/MWL(IU,JU)
              DGM = TLAKE(IU,JU)*DMM*SHW
              FLOW(IU,JU) =  FLOW(IU,JU) - DMM
              EFLOW(IU,JU) = EFLOW(IU,JU) - DGM
              AIJ(IU,JU,IJ_MRVR)=AIJ(IU,JU,IJ_MRVR) +  DMM
              AIJ(IU,JU,IJ_ERVR)=AIJ(IU,JU,IJ_ERVR)+
     *             (DGM+GRAV*DMM*ZATMO(IU,JU))
              IF(FOCEAN(ID,JD).le.0.) THEN
                FLOW(ID,JD) =  FLOW(ID,JD) + DMM
                EFLOW(ID,JD) = EFLOW(ID,JD) +
     +               DGM+GRAV*DMM*(ZATMO(IU,JU)-ZATMO(ID,JD))
              ELSE ! Save river mouth flow to for output to oceans
c                K=KDIREC(IU,JU)
                FLOWO(ID,JD)=FLOWO(ID,JD)+DMM
                EFLOWO(ID,JD)=EFLOWO(ID,JD)+DGM
     *               +GRAV*DMM*(ZATMO(IU,JU)-ZATMO(ID,JD))
C**** move application to ocean to ocean code
c      GXM(ID,JD,1) = GXM(ID,JD,1) - XK(K)*(DGM -
c     * (G0M(ID,JD,1)+XK(K)*GXM(ID,JD,1))*(DMM/(MO(ID,JD,1)*DXYPO(JD))))
c      GYM(ID,JD,1) = GYM(ID,JD,1) - YK(K)*(DGM -
c     * (G0M(ID,JD,1)+YK(K)*GYM(ID,JD,1))*(DMM/(MO(ID,JD,1)*DXYPO(JD))))
c      SXM(ID,JD,1) = SXM(ID,JD,1) + XK(K)*
c     * (S0M(ID,JD,1)+XK(K)*SXM(ID,JD,1))*(DMM/(MO(ID,JD,1)*DXYPO(JD)))
c      SYM(ID,JD,1) = SYM(ID,JD,1) + YK(K)*
c     * (S0M(ID,JD,1)+YK(K)*SYM(ID,JD,1))*(DMM/(MO(ID,JD,1)*DXYPO(JD)))
c       MO(ID,JD,1) =  MO(ID,JD,1) + BYDXYPO(JD)* DMM
c      G0M(ID,JD,1) = G0M(ID,JD,1) +
c     *               DGM+GRAV*DMM*(ZATMO(IU,JU)-ZATMO(ID,JD))
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
      MWLSILL = 1d3*HLAKE(1,1)*FLAKE(1,1)*DXYP(1)
      IF(MWL(1,1).gt.MWLSILL) THEN
        DMM = (MWL(1,1)-MWLSILL)*RATE(1,1)
C       DGM = GML(1,1)*DMM/MWL(1,1)
        DGM = TLAKE(1,1)*DMM*SHW
        FLOW(1,1) =  FLOW(1,1) - DMM
        EFLOW(1,1) = EFLOW(1,1) - DGM
          AIJ(1,1,IJ_MRVR)=AIJ(1,1,IJ_MRVR) +  DMM
          AIJ(1,1,IJ_ERVR)=AIJ(1,1,IJ_ERVR) + (DGM+GRAV*DMM*ZATMO(1,1))
         FLOW(IDPOLE,JDPOLE) =  FLOW(IDPOLE,JDPOLE) + IM*DMM
        EFLOW(IDPOLE,JDPOLE) = EFLOW(IDPOLE,JDPOLE) +
     +       IM*(DGM+GRAV*DMM*(ZATMO(1,1)-ZATMO(IDPOLE,JDPOLE)))
      END IF
C****
C**** Apply net river flow to continental reservoirs
C****
      DO J=2,JM-1
        DO I=1,IM
          IF(FOCEAN(I,J).le.0.)  THEN
            MWL(I,J) = MWL(I,J) +  FLOW(I,J)
            GML(I,J) = GML(I,J) + EFLOW(I,J)
          END IF
        END DO
      END DO
      MWL(1,1) = MWL(1,1) +  FLOW(1,1)
      GML(1,1) = GML(1,1) + EFLOW(1,1)
      RETURN
C**** 
      END SUBROUTINE RIVERF

      SUBROUTINE CHECKL (SUBRN)
!@sum  CHECKL checks whether the lake variables are reasonable.
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0 (based on LB265)
      USE CONSTANT, only : shi,lhm
      USE E001M12_COM, only : IM,JM,FLAKE,HLAKE
      USE GEOM, only : DXYP
      USE DAGCOM, only : QCHECK
      USE LAKES
      USE LAKES_COM

      IMPLICIT NONE
      
      INTEGER I,J !@var I,J loop variables
      CHARACTER*6 SUBRN
      LOGICAL QCHECKL

      QCHECKL = .FALSE.
      DO J=2,JM-1
        DO I=1,IM
          IF(FLAKE(I,J).gt.0.) THEN
C**** check for reasonable lake surface temps
            IF (TLAKE(I,J).ge.40 .or. TLAKE(I,J).le.-80) THEN
              WRITE(6,*) 'After ',SUBRN,': I,J,TSL=',I,J,TLAKE(I,J)
              QCHECKL = .TRUE.
            END IF
C**** Check potential specific enthalpy of first lake layer
c     GO1 = GML(I,J)/MWL(I,J)
c     IF(GO1.lt.0. .or. GO1.gt.200000.)  THEN
c     WRITE (6,*) 'After ',SUBRN,': I,J,GO1=',I,J,GO1,MWL(i,j),GML(i,j)
c     QCHECKL = .TRUE.
c     END IF
C**** Check first layer lake mass (<10%, >10x orig depth)
            IF(MWL(I,J).lt.1d2*HLAKE(I,J)*DXYP(J)*FLAKE(I,J)) THEN
              WRITE (6,*) 'After ',SUBRN,
     *             ': I,J,FLAKE,HLAKE,lake level low=',I,J,FLAKE(I,J),
     *             HLAKE(I,J),1d-3*MWL(I,J)/(DXYP(J)*FLAKE(I,J))
c     QCHECKL = .TRUE.
            END IF
            IF(MWL(I,J).gt.1d4*MAX(HLAKE(I,J),1d0)*DXYP(J)*FLAKE(I,J))
     *           THEN
              WRITE (6,*) 'After ',SUBRN,
     *             ': I,J,FLAKE,HLAKE,lake level high=',I,J,FLAKE(I,J),
     *             HLAKE(I,J),1d-3*MWL(I,J)/(DXYP(J)*FLAKE(I,J))
c     QCHECKL = .TRUE.
            END IF
          END IF
        END DO
      END DO
      IF (QCHECKL) STOP 'CHECKL: Lake variables out of bounds'
      RETURN
      END SUBROUTINE CHECKL

      SUBROUTINE daily_LAKE(IEND)
!@sum  daily_LAKE sets Lake Ice extent consistent with 50-day SAT
!@auth L. Nazarenko
!@ver  1.0
      USE CONSTANT, only : rhoi
      USE E001M12_COM, only : IM,JM,FLAKE,GDATA,KOCEAN,TAU,TAUI,DT
      USE GEOM, only : IMAXJ
!      USE LAKES, only : ZIMIN,ZIMAX,T_ICE,T_NOICE,DRSIDT
      USE LAKES_COM, only : T50
      USE OCEAN, only : ODATA,DM
      USE SEAICE, only : Z1I
      IMPLICIT NONE
      REAL*8  ZIMIN,ZIMAX,T_ICE,T_NOICE,DRSIDT
      INTEGER I,J,K,IMAX,IEND
      REAL*8 RSINEW

      IF (KOCEAN.eq.0) RETURN
      IF (IEND.eq.0 .and. TAU.GT.TAUI+DT/7200.) RETURN
      ZIMIN=.5     ! minimum lake ice thickness
      ZIMAX=2.     ! maximum lake ice thickness
      T_ICE = -8.  ! surface air temperature for 100% ice cover
      T_NOICE = 0. ! surface air temperature for no ice cover
      DRSIDT = 1./(T_ICE-T_NOICE)
      DO J=1,JM
        IMAX=IMAXJ(J)
        DO I=1,IMAX
          IF (FLAKE(I,J) .GT. 0.) THEN ! linear fit for -8< T50 <0
            RSINEW = MIN(1d0,MAX(0d0,(T50(I,J)-T_NOICE)*DRSIDT))
            ODATA(I,J,2)=RSINEW
            ODATA(I,J,3)=RHOI*(ZIMIN-Z1I+(ZIMAX-ZIMIN)*RSINEW*DM(I,J))
            IF (RSINEW.LE.0.) THEN
              GDATA(I,J,1)=0.
              GDATA(I,J,3)=0.
              GDATA(I,J,7)=0.
              GDATA(I,J,15)=0.
              GDATA(I,J,16)=0.
            END IF
          END IF
        END DO
      END DO

      RETURN
      END SUBROUTINE daily_LAKE

