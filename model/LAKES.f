      MODULE LAKES
!@sum  LAKES subroutines for Lakes and Rivers
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0 (based on LB265)
      USE E001M12_COM, only : IM,JM
      IMPLICIT NONE
      
C**** 
C**** Changes from Model III: MO -> MWL (kg), G0M -> GML (J), 
C****    GZM -> TSL (deg C),  RSI -> RLI (fract),
C****    MSI -> MLI (kg/m^2), HSI -> HLI (J/m^2)
C****
!@var LMLI Number of layers of lake ice model
      INTEGER, PARAMETER :: LMLI=4
!@var XLI1,XLI2,XLI3,XLI4 layering for lake ice model
      REAL*8, PARAMETER :: XLI1=.50,XLI2=.50,XLI3=.50,XLI4=.50
!@var BYXLI1,BYXLI2,BYXLI3,BYXLI4 reciprocals of XLIx
      REAL*8, PARAMETER :: BYXLI1=1./XLI1, BYXLI2=1./XLI2,
     *                     BYXLI3=1./XLI3, BYXLI4=1./XLI4
!@var TFL freezing temperature for lakes (=0 C)
      REAL*8, PARAMETER :: TFL = 0
!@var FLEAD lead fraction for lakes
      REAL*8, PARAMETER :: FLEAD=.05d0   ! = 0?
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
      USE E001M12_COM, only : IM,JM,GDATA,FLAKE,ZATMO,DT,FLICE,HLAKE
      USE OCEAN, only : ODATA,ACE1I
      USE GEOM
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
C****         TSL      Temperature of lake surface (C)
C****         RLI      Horizontal ratio of lake ice to lake (1)
C****         MLI      Mass of sea ice (kg/m**2)
C****         HLI      Heat including latent of sea ice (J/m**2)
C****         HLAKE    Lake sill depth (m)
C****
C**** FIXDCB  FLAKE    Lake fraction (1)  
C****

      IF (INILAKE) THEN
C**** Set lake variables from model blocks ODATA and GDATA
C**** GDATA(3) and (7) are first and second thermal layers
C**** This is just an estimate for the initiallisation
      DO J=2,JM-1
       DO I=1,IM
         IF (FLAKE(I,J).gt.0) THEN
           TSL(I,J)   = ODATA(I,J,1)
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
           GML(I,J) = MWL(I,J)*MAX(TSL(I,J),4d0)
C**** reset GDATA correctly
c           GDATA(I,J,3) = (HLI(I,J,1)/(XLI1*MLI(I,J,1))+LHM)/SHI
c           GDATA(I,J,7) = (HLI(I,J,2)/(XLI2*MLI(I,J,1))+LHM)/SHI
         ELSE
           TSL(I,J)   = 0.
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
          IF(FLICE(I,J).GT.0.)  TSL(I,J) = 0
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
      USE E001M12_COM, only : IM,JM,FLAKE,DT,FOCEAN,ZATMO,HLAKE
      USE GEOM
      USE LAKES
      USE LAKES_COM
      USE DAGCOM, only : AIJ,IJ_ERVR,IJ_MRVR
     
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
C****         TSL  Lake surface temperature (C)
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
              DGM = TSL(IU,JU)*DMM*SHW
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
        DGM = TSL(1,1)*DMM*SHW
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
      REAL*8 TLI1,TLI2,TLI3,TLI4 !@var TLIn lake ice temperatures (C)
      CHARACTER*6 SUBRN
      LOGICAL QCHECKL
C****
C**** Check lake reservoirs for reasonable values
C****
      QCHECKL = .FALSE.
      DO J=2,JM-1
      DO I=1,IM
      IF(FLAKE(I,J).gt.0.) THEN
C**** check for reasonable lake surface temps
      IF (TSL(I,J).ge.40 .or. TSL(I,J).le.-80) THEN
      WRITE(6,*) 'After ',SUBRN,': I,J,TSL=',I,J,TSL(I,J)
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
      WRITE (6,*) 'After ',SUBRN,': I,J,FLAKE,HLAKE,lake level low=',
     *     I,J,FLAKE(I,J),HLAKE(I,J),1d-3*MWL(I,J)/(DXYP(J)*FLAKE(I,J))
c      QCHECKL = .TRUE.
      END IF
      IF(MWL(I,J).gt.1d4*MAX(HLAKE(I,J),1d0)*DXYP(J)*FLAKE(I,J)) THEN
      WRITE (6,*) 'After ',SUBRN,': I,J,FLAKE,HLAKE,lake level high=',
     *     I,J,FLAKE(I,J),HLAKE(I,J),1d-3*MWL(I,J)/(DXYP(J)*FLAKE(I,J))
c      QCHECKL = .TRUE.
      END IF
c      IF(FLAKE(I,J)*RLI(I,J).gt.0.)  THEN
C**** Make sure mLI1 and mLI2 are  +ve
c      IF (MLI(I,J,1).GT.500.OR.MLI(I,J,1).LT.-.001) THEN
c      WRITE(6,*) 'After', SUBRN,': I,J,MLI1=',I,J,MLI(I,J,1),RLI(I,J)
c      QCHECKL = .TRUE.
c      ENDIF
c      IF (FLAKE(I,J)*RLI(I,J)*MLI(I,J,2).LT.0.) THEN
c      WRITE(6,*) 'After', SUBRN,': I,J,MLI2=',I,J,MLI(I,J,2),RLI(I,J)
c      QCHECKL = .TRUE.
c      ENDIF
C**** Check temperature of first layer over lake ice
c      TLI1 = (HLI(I,J,1)/(XLI1*MLI(I,J,1)) + LHM) / SHI
c      TLI2 = (HLI(I,J,2)/(XLI2*MLI(I,J,1)) + LHM) / SHI
c      TLI3 = (HLI(I,J,3)/(XLI3*MLI(I,J,2)) + LHM) / SHI
c      TLI4 = (HLI(I,J,4)/(XLI4*MLI(I,J,2)) + LHM) / SHI
c      IF(TLI1.gt.1.D-12 .or. TLI2.gt.1.D-12 .or. TLI1.lt.-80. .or.
c     *   TLI3.gt.1.D-12 .or. TLI4.gt.1.D-12)  THEN
c      WRITE (6,*) 'After ',SUBRN,': I,J,MLI1,TICE,FLAKE,RLI=',
c     *       I,J,MLI(I,J,1),TLI1,TLI2,TLI3,TLI4,FLAKE(I,J),RLI(I,J)
c      QCHECKL = .TRUE.
c      END IF
C**** Check horizontal lake ice ratio to ocean
c      IF(RLI(I,J).lt.0. .or. RLI(I,J).gt.1d0)  THEN
c      WRITE (6,*) 'After ',SUBRN,': I,J,RLI=',I,J,RLI(I,J)
c      QCHECKL = .TRUE.
c      END IF
c      END IF
      END IF
      END DO
      END DO
      IF (QCHECKL) STOP 'CHECKL: Lake variables out of bounds'
      RETURN
      END SUBROUTINE CHECKL

c$$$      SUBROUTINE PRECLI
c$$$!@sum  PRECLI applies precipitation to lake ice 
c$$$!@auth Gary Russell
c$$$!@ver  1.0 (based on LB265)
c$$$      USE CONSTANT, only : shi,lhm
c$$$      USE E001M12_COM, only : IM,JM,GDATA,FLAKE
c$$$      USE OCEAN, only : ODATA
c$$$      USE LAKES
c$$$      USE LAKES_COM
c$$$      USE CLD01_COM_E001, only : PRCP=>PREC,TPRCP=>TPREC
c$$$
c$$$      IMPLICIT NONE
c$$$
c$$$      INTEGER I,J !@var I,J loop variables
c$$$      REAL*8 MELT1 !@var MELT1 melt from first layer
c$$$C****
c$$$C**** LAKECB  MWL    Liquid lake mass  (kg)
c$$$C****         GML    Liquid lake enthalpy  (J)
c$$$C****         FLAKE  Lake fraction (1)
c$$$C****         RLI    Horizontal ratio of lake ice to lake (1)
c$$$C****         MLI    Mass of sea ice (kg/m**2)
c$$$C****         HLI    Heat including latent of sea ice (J/m**2)
c$$$C****
c$$$C**** WORK1   PREC   Precipitation from atmosphere (kg/m**2)
c$$$C****         EPREC  Energy of precipitation (J/m**2)
c$$$C****
c$$$C**** Outside loop over J and I
c$$$C****
c$$$      DO 600 J=2,JM-1
c$$$      IMAX=IM
c$$$      IF(J.EQ.1 .OR. J.EQ.JM)  IMAX=1
c$$$      DO 600 I=1,IM
c$$$      FLICE  = FLAKE(I,J)*RLI(I,J)
c$$$      PREC = PRCP(I,J)
c$$$      IF(PREC.LE.0.OR.FLICE.LE.0.)  GO TO 600
c$$$      IF (TPRCP(I,J).LT.0.) THEN
c$$$C       EPRCP=PREC*TPRCP(I,J)*SHI
c$$$        EPRCP=0.
c$$$        EPREC=EPRCP-PREC*LHM
c$$$      ELSE
c$$$C       EPRCP=PREC*TPRCP(I,J)*SHW
c$$$        EPRCP=0.
c$$$        EPREC=EPRCP
c$$$      END IF
c$$$         JR=JREG(I,J)
c$$$C****
c$$$C**** Apply precipitation heat flux to HLI1
c$$$C****
c$$$c     MLI(I,J,1) = MLI(I,J,1) + ACE1I
c$$$      HLI(I,J,1) = HLI(I,J,1) + EPREC
c$$$      IF(EPREC.LE.-PREC*LHM)  GO TO 300
c$$$      IF (EPREC.LE.0.) GO TO 200
c$$$C****
c$$$C**** All precipitation is rain above 0 deg C
c$$$C**** Rain compresses snow amount into ice
c$$$C**** Calculate snow and ice melted or rain frozen in layer 1
c$$$C****
c$$$      RAIN = PREC
c$$$      IF(HLI(I,J,1)/LHM+XLI1*MLI(I,J,1) .le. 0.)  GO TO 140
c$$$C**** Warm rain cools to 0 C and melts some snow or ice
c$$$      MELT1 = HLI(I,J,1)/LHM + XLI1*MLI(I,J,1)
c$$$      MWL(I,J) = MWL(I,J) + PLKICE*(RAIN+MELT1)*DXYP(J)
c$$$      CJ(J,54) = CJ(J,54) + PLKICE*(RAIN+MELT1)
c$$$      DJ(JR,54) = DJ(JR,54) + PLKICE*(RAIN+MELT1)*DXYP(J)
c$$$      IF(MLI(I,J,1)-ACE1I .le. MELT1)  GO TO 120
c$$$C**** Rain and melting compresses snow into ice
c$$$      FMLI2 = MIN (dSNdRN*(RAIN+MELT1), MLI(I,J,1)-ACE1I-MELT1)
c$$$      FHLI1 = - LHM*(XLI1*FMLI2 - XLI2*MELT1)
c$$$      FHLI2 = HLI(I,J,2)*FMLI2*RXLI2 / MLI(I,J,1)
c$$$      MLI(I,J,1) = MLI(I,J,1) - MELT1 - FMLI2
c$$$      GO TO 400
c$$$C**** All snow and some ice has melted
c$$$C     FMLI1 = XLI1*(MLI1-ACE1I) - MELT1     !  < 0
c$$$  120 FMLI2 = MLI(I,J,1) - ACE1I - MELT1  !  < 0
c$$$      FHLI1 = HLI(I,J,2)*((XLI1/XLI2)*FMLI2-MELT1) / MLI(I,J,1)
c$$$      FHLI2 = HLI(I,J,3)*FMLI2*RXLI3 / MLI(I,J,2)
c$$$      FHLI3 = HLI(I,J,4)*FMLI2    / MLI(I,J,2)
c$$$      MLI(I,J,1) = ACE1I
c$$$      GO TO 420
c$$$C**** Rain compresses snow into ice, some rain will freeze
c$$$  140 CMPRS = MIN (dSNdRN*RAIN, MLI(I,J,1)-ACE1I)
c$$$      IF(-HLI(I,J,1)/LHM-XLI1*MLI(I,J,1) .lt. RAIN)  GO TO 150
c$$$C**** All rain freezes in layer 1
c$$$C     FREZ1 = RAIN
c$$$C     FMLI1 = XLI1*CMPRS + FREZ1
c$$$      FMLI2 = CMPRS + RAIN
c$$$      FHLI1 = HLI(I,J,1)*(XLI1*CMPRS+RAIN) / (XLI1*MLI(I,J,1)+RAIN)
c$$$      FHLI2 = HLI(I,J,2)*FMLI2*RXLI2 / MLI(I,J,1)
c$$$      MLI(I,J,1) = MLI(I,J,1) - CMPRS
c$$$      GO TO 400
c$$$C**** Not all rain freezes in layer 1
c$$$  150 FREZ1 = - HLI(I,J,1)/LHM - XLI1*MLI(I,J,1)
c$$$      MWL(I,J) = MWL(I,J) + PLKICE*(RAIN-FREZ1)*DXYP(J)
c$$$      CJ(J,54) = CJ(J,54) + PLKICE*(RAIN-FREZ1)
c$$$      DJ(JR,54) = DJ(JR,54) + PLKICE*(RAIN-FREZ1)*DXYP(J)
c$$$C     FMLI1 = XLI1*CMPRS + FREZ1
c$$$      FMLI2 = CMPRS + FREZ1
c$$$      FHLI1 = - LHM*(XLI1*CMPRS + FREZ1)
c$$$      FHLI2 = HLI(I,J,2)*FMLI2*RXLI2 / MLI(I,J,1)
c$$$      MLI(I,J,1) = MLI(I,J,1) - CMPRS
c$$$      GO TO 400
c$$$C****
c$$$C**** Precipitation is a mixture of rain and snow at 0C
c$$$C**** Rain compresses snow into ice, some rain will freeze
c$$$C**** Note that snow amount may exceed SNOMAX
c$$$C****
c$$$  200 SNOW = -EPREC/LHM
c$$$      RAIN  = PREC - SNOW
c$$$      CMPRS = MIN (dSNdRN*RAIN, MLI(I,J,1)+SNOW-ACE1I)
c$$$      IF(-HLI(I,J,1)/LHM-XLI1*MLI(I,J,1)-SNOW .lt. RAIN)  GO TO 210
c$$$C**** All rain freezes in layer 1
c$$$C     FREZ1 = RAIN
c$$$      FMLI1 = XLI1*CMPRS + XLI2*SNOW + RAIN
c$$$      FMLI2 = CMPRS + RAIN
c$$$      FHLI1 = HLI(I,J,1)*FMLI1 / (XLI1*MLI(I,J,1)+PREC)
c$$$      FHLI2 = HLI(I,J,2)*FMLI2*RXLI2 / MLI(I,J,1)
c$$$      MLI(I,J,1) = MLI(I,J,1) + (SNOW-CMPRS)
c$$$      GO TO 400
c$$$C**** Not all rain freezes in layer 1
c$$$  210 FREZ1 = - HLI(I,J,1)/LHM - XLI1*MLI(I,J,1) - SNOW
c$$$      MWL(I,J) = MWL(I,J) + PLKICE*(RAIN-FREZ1)*DXYP(J)
c$$$      CJ(J,54) = CJ(J,54) + PLKICE*(RAIN-FREZ1)
c$$$      DJ(JR,54) = DJ(JR,54) + PLKICE*(RAIN-FREZ1)*DXYP(J)
c$$$C     FMLI1 = XLI1*CMPRS + XLI2*SNOW + FREZ1
c$$$      FMLI2 = CMPRS + FREZ1
c$$$      FHLI1 = - LHM*(XLI1*CMPRS + XLI2*SNOW + FREZ1)
c$$$      FHLI2 = HLI(I,J,2)*FMLI2*RXLI2 / MLI(I,J,1)
c$$$      MLI(I,J,1) = MLI(I,J,1) + (SNOW-CMPRS)
c$$$      GO TO 400
c$$$C****
c$$$C**** All precipitation is snow which increases snow amount
c$$$C****
c$$$  300 IF(MLI(I,J,1)+PREC .gt. SNOMAX+ACE1I)  GO TO 310
c$$$C     FMLI1 = XLI2*PREC  !  > 0
c$$$      FHLI1 = HLI(I,J,1)*XLI2*PREC /(XLI1*MLI(I,J,1)+PREC)
c$$$      MLI(I,J,1) = MLI(I,J,1) + PREC
c$$$      HLI(I,J,1) = HLI(I,J,1) - FHLI1
c$$$      HLI(I,J,2) = HLI(I,J,2) + FHLI1
c$$$      GO TO 500
c$$$C**** Too much snow has accumulated, compress some into ice
c$$$  310 FMLI2 = MLI(I,J,1) + PREC - (.9d0*SNOMAX+ACE1I)  !  > 0
c$$$C     FMLI1 = XLI2*PREC + XLI1*FMLI2                        !  > 0
c$$$      FHLI1 = HLI(I,J,1)*(XLI2*PREC+XLI1*FMLI2) /
c$$$     /        (XLI1*MLI(I,J,1)+PREC)
c$$$      FHLI2 = HLI(I,J,2)*FMLI2*RXLI2 / MLI(I,J,1)
c$$$      MLI(I,J,1) = .9d0*SNOMAX+ACE1I
c$$$C****
c$$$C**** Advect ice (usually downward)
c$$$C****
c$$$  400 FHLI3 = HLI(I,J,3)*FMLI2*(XLI4/XLI3) / MLI(I,J,2)
c$$$  420 MLI(I,J,2) = MLI(I,J,2) +  FMLI2
c$$$      HLI(I,J,1) = HLI(I,J,1) -  FHLI1
c$$$      HLI(I,J,2) = HLI(I,J,2) + (FHLI1 - FHLI2)
c$$$      HLI(I,J,3) = HLI(I,J,3) + (FHLI2 - FHLI3)
c$$$      HLI(I,J,4) = HLI(I,J,4) +  FHLI3
c$$$  500 CONTINUE
c$$$C**** Update model blocks ODATA and GDATA
c$$$      ODATA(I,J,3) = MLI(I,J,2)
c$$$      GDATA(I,J,1) = MLI(I,J,1)-ACE1I
c$$$      GDATA(I,J,3) = (HLI(I,J,1)/(XLI1*MLI(I,J,1))+LHM)/SHI
c$$$      GDATA(I,J,7) = (HLI(I,J,2)/(XLI2*MLI(I,J,1))+LHM)/SHI
c$$$  600 CONTINUE
c$$$      RETURN
c$$$      END SUBROUTINE PRECLI
c$$$
c$$$      SUBROUTINE LAKICE
c$$$!@sum  LAKES subroutines for Lakes and Rivers
c$$$!@auth Gavin Schmidt/Gary Russell
c$$$!@ver  1.0 (based on LB265)
c$$$C****
c$$$C**** LAKICE receives the atmospheric fluxes of heat and water and
c$$$C**** applies them to the lake ice, updating the heat, mass and snow
c$$$C****
c$$$      USE CONSTANT, only : shi,lhm
c$$$      USE E001M12_COM, only : IM,JM,FLAKE
c$$$      USE LAKES
c$$$      USE LAKES_COM
c$$$
c$$$      IMPLICIT NONE
c$$$
c$$$      INTEGER I,J !@var I,J loop variables
c$$$      REAL*8 MELT1,MELT4
c$$$      PARAMETER (ALPHA=1.)
c$$$      COMMON/WORK3/F0DT(IM,JM,4),F1DT(IM,JM,4),EVAPOR(IM,JM,4)
c$$$C****
c$$$C**** LAKECB  MWL    Liquid lake mass  (kg)
c$$$C****         GML    Liquid lake enthalpy  (J)
c$$$C****         FLAKE  Lake fraction (1)
c$$$C****         RLI    Horizontal ratio of lake ice to lake (1)
c$$$C****         MLI    Mass of sea ice (kg/m**2)
c$$$C****         HLI    Enthalpy minus latent of sea ice (J/m**2)
c$$$C****
c$$$C**** WORK3   E0(2)  Energy received from atmosphere (J/m**2)
c$$$C****         E1(2)  Energy from layer 1 to layer 2 (J/m**2)
c$$$C****         EVAP(2)  Evaporation from ground (kg/m**2)
c$$$C****
c$$$C**** Outside loop over J and I
c$$$C****
c$$$      DO 500 J=2,JM-1
c$$$      IMAX=IM
c$$$      IF(J.EQ.1 .OR. J.EQ.JM)  IMAX=1
c$$$      DO 500 I=1,IM
c$$$      FLICE  = FLAKE(I,J)*RLI(I,J)
c$$$      IF(FLICE.LE.0.)  GO TO 500
c$$$      E0=F0DT(I,J,2)
c$$$      E1=F1DT(I,J,2)
c$$$      EVAP=EVAPOR(I,J,2)
c$$$      JR=JREG(I,J)
c$$$c     MLI(I,J,1) = MLI(I,J,1) + ACE1I
c$$$      TLI1 = (HLI(I,J,1)*RXLI1/MLI(I,J,1) + LHM)/SHI
c$$$      TLI2 = (HLI(I,J,2)*RXLI2/MLI(I,J,1) + LHM)/SHI
c$$$      TLI3 = (HLI(I,J,3)*RXLI3/MLI(I,J,2) + LHM)/SHI
c$$$      TLI4 = (HLI(I,J,4)*RXLI4/MLI(I,J,2) + LHM)/SHI
c$$$C**** Accumulate diagnostics
c$$$      CJ(J,30) = CJ(J,30) + PLKICE
c$$$      CJ(J,50) = CJ(J,50) + PLKICE*ACE1I
c$$$      CJ(J,52) = CJ(J,52) + PLKICE*MLI(I,J,2)
c$$$      CJ(J,18) = CJ(J,18) + PLKICE*TLI1   
c$$$      CJ(J,17) = CJ(J,17) + PLKICE*TLI2
c$$$      CJ(J,42) = CJ(J,42) + PLKICE*E1
c$$$      CJ(J,19) = CJ(J,19) + PLKICE*EVAP
c$$$c      CJ(J,39) = CJ(J,39) + PLKICE*EPREC
c$$$      CJ(J,53)   = CJ(J,53)   + PLKICE*(MLI(I,J,1)-ACE1I)
c$$$      IF (JR.NE.24) THEN
c$$$       DJ(JR,30) = DJ(JR,30) + PLKICE*DXYP(J)
c$$$       DJ(JR,53) = DJ(JR,30) + PLKICE*DXYP(J)*(MLI(I,J,1)-ACE1I)
c$$$      DJ(JR,50)= DJ(JR,50) + PLKICE*DXYP(J)*ACE1I
c$$$       DJ(JR,52) = DJ(JR,52) + PLKICE*DXYP(J)*MLI(I,J,2)
c$$$      ENDIF
c$$$c     AIJ(I,J,1) = AIJ(I,J, 1) + PLKICE
c$$$c     AIJ(I,J,62) = AIJ(I,J,62) + PLKICE*EVAP
c$$$c     AIJ(I,J,66) = AIJ(I,J,66) + PLKICE*E0
c$$$c     AIJ(I,J,58) = AIJ(I,J,58) + PLKICE*MLI(I,J,2)
c$$$C****
c$$$C**** Calculate and apply diffusive and surface energy fluxes
c$$$C****
c$$$   10 HCI2 = SHI*XLI2*MLI(I,J,1)
c$$$      HCI3 = SHI*XLI3*MLI(I,J,2)
c$$$      HCI4 = SHI*XLI4*MLI(I,J,2)
c$$$      dF2dTI = ALAMI*RHOI*DT / (.5*XLI2*MLI(I,J,1)+.5*XLI3*MLI(I,J,2))
c$$$      dF3dTI = ALAMI*RHOI*DT*2.    / MLI(I,J,2)
c$$$      dF4dTI = ALAMI*RHOI*DT*2.*RXLI4 / MLI(I,J,2)
c$$$CEXP  F2 = dF2dTI*(TLI2-TLI3)
c$$$CEXP  F3 = dF3dTI*(TLI3-TLI4)
c$$$CEXP  F4 = dF4dTI*(TLI4-TFL)
c$$$      F2 = dF2dTI*(HCI2*(TLI2-TLI3) + ALPHA*E1) /
c$$$     /     (HCI2 + ALPHA*dF2dTI)
c$$$      F3 = dF3dTI*(HCI3*(TLI3-TLI4) + ALPHA*F2) / (HCI3 + ALPHA*dF3dTI)
c$$$      F4 = dF4dTI*(HCI4*(TLI4-TFL) + ALPHA*F3) / (HCI4 + ALPHA*dF4dTI)
c$$$      HLI(I,J,1) = HLI(I,J,1) + (E0 - E1)
c$$$      HLI(I,J,2) = HLI(I,J,2) + (E1 - F2)
c$$$      HLI(I,J,3) = HLI(I,J,3) + (F2 - F3)
c$$$      HLI(I,J,4) = HLI(I,J,4) + (F3 - F4)
c$$$      GML(I,J) = GML(I,J) + PLKICE*F4*DXYP(J)
c$$$         CJ(J,15) = CJ(J,15)  + PLKICE*F4
c$$$      IF(HLI(I,J,1)/LHM+XLI1*MLI(I,J,1)-EVAP .le. 0.)  GO TO 100
c$$$C****
c$$$C**** Fluxes heat layer 1 to freezing point and melt some snow or ice
c$$$C****
c$$$      MELT1 = HLI(I,J,1)/LHM + XLI1*MLI(I,J,1) - EVAP
c$$$      MWL(I,J) = MWL(I,J) + PLKICE*MELT1*DXYP(J)
c$$$        DJ(JR,54) = DJ(JR,54) + PLKICE*MELT1*DXYP(J)
c$$$        CJ(J,54) = CJ(J,54) + PLKICE*MELT1
c$$$      IF(EVAP .lt. 0.)  GO TO 30
c$$$C**** Evaporation reduces snow or ice amount, MELT1 > 0
c$$$      IF(MLI(I,J,1)-ACE1I-EVAP .gt. MELT1)  GO TO 20
c$$$C**** All snow and some ice evaporates and melts
c$$$C**** Ice advection is upward into layer 2 from layer 3
c$$$      FMLI1 = XLI1*(MLI(I,J,1)-ACE1I) - EVAP - MELT1  !  < 0
c$$$      FMLI2 =       MLI(I,J,1)-ACE1I  - EVAP - MELT1  !  < 0
c$$$      GO TO 110
c$$$C**** Evaporation and melting reduce snow amount
c$$$C**** Ice advection is downward from layer 2 into layer 3
c$$$   20 FMLI2 = MIN (dSNdML*MELT1, MLI(I,J,1)-ACE1I-EVAP-MELT1) !>0
c$$$      FMLI1 = XLI1*FMLI2 - XLI2*(EVAP + MELT1)
c$$$      IF(FMLI1.lt.0.)  FHLI1 = HLI(I,J,2)*FMLI1*RXLI2 / MLI(I,J,1)
c$$$      IF(FMLI1.ge.0.)  FHLI1 = - LHM*FMLI1
c$$$      FHLI2 = HLI(I,J,2)*FMLI2*RXLI2 / MLI(I,J,1)
c$$$      MLI(I,J,1) = MLI(I,J,1) - EVAP - MELT1 - FMLI2
c$$$      GO TO 310
c$$$C**** Dew increases ice amount, MELT1 > 0
c$$$   30 IF(MLI(I,J,1)-ACE1I .gt. MELT1)  GO TO 40
c$$$C**** All snow and some ice melts
c$$$      FMLI1 = XLI1*(MLI(I,J,1)-ACE1I) - EVAP - MELT1
c$$$      FMLI2 = MLI(I,J,1) - ACE1I - EVAP - MELT1
c$$$      IF(FMLI2.le.0.)  GO TO 110
c$$$C**** Ice advection is downward from layer 2 into layer 3
c$$$      IF(FMLI1.lt.0.)  FHLI1 = HLI(I,J,2)*FMLI1*RXLI2 / MLI(I,J,1)
c$$$      IF(FMLI1.ge.0.)  FHLI1 = - LHM*FMLI1
c$$$      FHLI2 = HLI(I,J,2)*FMLI2*RXLI2 / MLI(I,J,1)
c$$$      MLI(I,J,1) = ACE1I
c$$$      GO TO 310
c$$$C**** Melting reduces snow amount
c$$$C**** Ice advection is downward from layer 2 into layer 3
c$$$   40 CMPRS = MIN (dSNdML*MELT1, MLI(I,J,1)-ACE1I-MELT1)
c$$$      FMLI1 = -EVAP + XLI1*CMPRS - XLI2*MELT1
c$$$      FMLI2 = -EVAP + CMPRS  !  > 0
c$$$      IF(FMLI1.lt.0.)  FHLI1 = HLI(I,J,2)*FMLI1*RXLI2 / MLI(I,J,1)
c$$$      IF(FMLI1.ge.0.)  FHLI1 = - LHM*FMLI1
c$$$      FHLI2 = HLI(I,J,2)*FMLI2*RXLI2 / MLI(I,J,1)
c$$$      MLI(I,J,1) = MLI(I,J,1) - MELT1 - CMPRS
c$$$      GO TO 310
c$$$C****
c$$$C**** No snow or ice melts in layer 1
c$$$C****
c$$$  100 IF(EVAP .lt. 0.)  GO TO 300
c$$$      IF(MLI(I,J,1)-ACE1I-EVAP .ge. 0.)  GO TO 200
c$$$C****
c$$$C**** All snow and some ice evaporates
c$$$C**** Ice advection is upward into layer 2 from layer 3
c$$$C****
c$$$      FMLI1 = XLI1*(MLI(I,J,1)-ACE1I) - EVAP  !  < 0
c$$$      FMLI2 = MLI(I,J,1) - ACE1I - EVAP       !  < 0
c$$$  110 FHLI1 = HLI(I,J,2)*FMLI1*RXLI2 / MLI(I,J,1)
c$$$      FHLI2 = HLI(I,J,3)*FMLI2*RXLI3 / MLI(I,J,2)
c$$$      MLI(I,J,1) = ACE1I
c$$$      IF(HLI(I,J,4)/LHM+XLI4*MLI(I,J,2) .le. 0.)  GO TO 120
c$$$C**** Fluxes heat layer 4 to freezing point and melt some ice
c$$$      MELT4 = HLI(I,J,4)/LHM + XLI4*MLI(I,J,2)
c$$$      MWL(I,J) = MWL(I,J) + PLKICE*MELT4*DXYP(J)
c$$$        CJ(J,46) = CJ(J,46) + PLKICE*MELT4
c$$$      FMLI3 = XLI4*FMLI2 + XLI3*MELT4
c$$$      IF(FMLI3.le.0.)  FHLI3 = - LHM*FMLI3
c$$$      IF(FMLI3.gt.0.)  FHLI3 = HLI(I,J,3)*FMLI3*RXLI3 / MLI(I,J,2)
c$$$      MLI(I,J,2) = MLI(I,J,2) + FMLI2 - MELT4
c$$$      GO TO 330
c$$$C**** No ice melts in layer 4
c$$$C     FMLI3 = XLI4*FMLI2
c$$$  120 FHLI3 = HLI(I,J,4)*FMLI2 / MLI(I,J,2)
c$$$      MLI(I,J,2) = MLI(I,J,2) +  FMLI2
c$$$      GO TO 330
c$$$C****
c$$$C**** Some snow evaporates
c$$$C**** No ice advection between layers 2 and 3
c$$$C****
c$$$C     FMLI1 = -XLI2*EVAP  !  < 0
c$$$C     FMLI2 = 0
c$$$  200 FHLI1 = -HLI(I,J,2)*EVAP / MLI(I,J,1)
c$$$      MLI(I,J,1) = MLI(I,J,1) - EVAP
c$$$      HLI(I,J,1) = HLI(I,J,1) - FHLI1
c$$$      HLI(I,J,2) = HLI(I,J,2) + FHLI1
c$$$      IF(HLI(I,J,4)/LHM+XLI4*MLI(I,J,2) .le. 0.)  GO TO 400
c$$$C**** Fluxes heat layer 4 to freezing point and melt some ice
c$$$      MELT4 = HLI(I,J,4)/LHM + XLI4*MLI(I,J,2)
c$$$      MWL(I,J) = MWL(I,J) + PLKICE*MELT4*DXYP(J)
c$$$c       CJ(J,22) = CJ(J,22) + PLKICE*MELT4*FDATA(I,J,1)
c$$$c       DJ(J,22) = DJ(J,22) - PLKICE*MELT4*FDATA(I,J,1)
c$$$C       CJ(J,38) = CJ(J,38) + PLKICE*MELT4  !  in DIAGJ
c$$$        CJ(J,46) = CJ(J,46) + PLKICE*MELT4
c$$$C     FMLI3 = XLI3*MELT4  !  > 0
c$$$      FHLI3 = HLI(I,J,3)*MELT4 / MLI(I,J,2)
c$$$      MLI(I,J,2) = MLI(I,J,2) - MELT4
c$$$      HLI(I,J,3) = HLI(I,J,3) - FHLI3
c$$$      HLI(I,J,4) = HLI(I,J,4) + FHLI3
c$$$      GO TO 400
c$$$C****
c$$$C**** Dew increases ice amount, advect ice downward
c$$$C**** Ice advection is downward from layer 2 into layer 3
c$$$C****
c$$$C     FMLI1 = -EVAP  !  > 0
c$$$  300 FMLI2 = -EVAP  !  > 0
c$$$      FHLI1 = HLI(I,J,1)*FMLI2 / (XLI1*MLI(I,J,1)-EVAP)
c$$$      FHLI2 = HLI(I,J,2)*FMLI2*RXLI2 / MLI(I,J,1)
c$$$  310 IF(HLI(I,J,4)/LHM+XLI4*MLI(I,J,2) .le. 0.)  GO TO 320
c$$$C**** Fluxes heat layer 4 to freezing point and melt some ice
c$$$      MELT4 = HLI(I,J,4)/LHM + XLI4*MLI(I,J,2)
c$$$      MWL(I,J) = MWL(I,J) + PLKICE*MELT4*DXYP(J)
c$$$        CJ(J,46) = CJ(J,46) + PLKICE*MELT4
c$$$C     FMLI3 = XLI4*FMLI2 + XLI3*MELT4  !  > 0
c$$$      FHLI3 = HLI(I,J,3)*((XLI4/XLI3)*FMLI2+MELT4) / MLI(I,J,2)
c$$$      MLI(I,J,2) = MLI(I,J,2) + FMLI2 - MELT4
c$$$      GO TO 330
c$$$C**** No ice melts in layer 4
c$$$C     FMLI3 = XLI4*FMLI2  !  > 0
c$$$  320 FHLI3 = HLI(I,J,3)*(XLI4/XLI3)*FMLI2 / MLI(I,J,2)
c$$$      MLI(I,J,2) = MLI(I,J,2) +  FMLI2
c$$$  330 HLI(I,J,1) = HLI(I,J,1) -  FHLI1
c$$$      HLI(I,J,2) = HLI(I,J,2) + (FHLI1 - FHLI2)
c$$$      HLI(I,J,3) = HLI(I,J,3) + (FHLI2 - FHLI3)
c$$$      HLI(I,J,4) = HLI(I,J,4) +  FHLI3
c$$$c 395 MLI(I,J,1) = MLI(I,J,1) -ACE1I
c$$$  400 continue
c$$$C**** Update model blocks ODATA and GDATA
c$$$      ODATA(I,J,3) = MLI(I,J,2)
c$$$      GDATA(I,J,1) = MLI(I,J,1)-ACE1I
c$$$      GDATA(I,J,3) = (HLI(I,J,1)/(XLI1*MLI(I,J,1))+LHM)/SHI
c$$$      GDATA(I,J,7) = (HLI(I,J,2)/(XLI2*MLI(I,J,1))+LHM)/SHI
c$$$  500 CONTINUE
c$$$      RETURN
c$$$      END SUBROUTINE LAKICE
c$$$
c$$$      SUBROUTINE LKCLIM(I,J)
c$$$!@sum  LAKES subroutines for Lakes and Rivers
c$$$!@auth Gavin Schmidt/Gary Russell
c$$$!@ver  1.0 (based on LB265)
c$$$C****
c$$$C**** LKCLIM calculates new lake arrays from interpolated monthly
c$$$C**** mean climatological values.
c$$$C**** Output: TSL (C)      : lake surface temperature
c$$$C****         RLI (%)      : horizontal ratio of lake ice to lake
c$$$C****         MLI (kg/m**2): mass of sea ice
c$$$C****         HLI (J/m**2) : lake ice enthalpy + latent energy
c$$$C****         MWL (kg)     : liquid lake mass 
c$$$C****         GML (J)      : liquid lake enthalpy
c$$$C****
c$$$      USE CONSTANT, only : shi,lhm
c$$$      USE E001M12_COM, only : IM,JM,FLAKE
c$$$      USE LAKES
c$$$      USE LAKES_COM
c$$$
c$$$      IMPLICIT NONE
c$$$
c$$$      INTEGER I,J !@var I,J loop variables
c$$$      REAL*8 MLINEW
c$$$C****
c$$$      TSL(I,J) = ODATA(I,J,1)
c$$$      RLINEW   = ODATA(I,J,2)
c$$$      MLINEW   = ODATA(I,J,3)
c$$$      IF(RLINEW.le.0. .and. RLI(I,J).le.0.)  GO TO 320
c$$$c     IF(RLINEW.gt.0.) RLINEW = MIN(RLINEW,1.-FLEAD*RHOI/(ACE1I+MLINEW))
c$$$      IF(RLINEW.gt.RLI(I,J))  GO TO 310
c$$$C**** Some lake ice has melted
c$$$      FHLI3 = HLI(I,J,3)*(MLI(I,J,2)-MLINEW)/(MLI(I,J,2)+1d-20)
c$$$      IF(MLINEW.gt.MLI(I,J,2))
c$$$     *FHLI3 = HLI(I,J,4)*(MLI(I,J,2)-MLINEW)*(XLI4/XLI3) /(MLI(I,J,2)
c$$$     + +1d-20)
c$$$      FHLI4 = HLI(I,J,4)*(MLI(I,J,2)-MLINEW)/XLI4/(MLI(I,J,2)+1d-20)
c$$$      MWL(I,J) = MWL(I,J) + FLAKE(I,J)*(RLINEW*(MLI(I,J,2)-MLINEW) +
c$$$     +         (RLI(I,J)-RLINEW)*(MLI(I,J,1)+MLI(I,J,2)))*DXYP(J)
c$$$      GML(I,J) = GML(I,J) + FLAKE(I,J)*(RLINEW*FHLI4 +
c$$$     +  (RLI(I,J)-RLINEW)*(HLI(I,J,1)+HLI(I,J,2)+HLI(I,J,3)+HLI(I,J,4)))
c$$$     *  *DXYP(J)
c$$$c      DJ(J,36) = DJ(J,36) - FLAKE(I,J)*(RLINEW*(MLI(I,J,2)-MLINEW) +
c$$$c    +             (RLI(I,J)-RLINEW)*(MLI(I,J,1)+MLI(I,J,2)))
c$$$c     DJ(J,45) = DJ(J,45) - FLAKE(I,J)*(RLINEW*FHLI4 +
c$$$c    +  (RLI(I,J)-RLINEW)*(HLI(I,J,1)+HLI(I,J,2)+HLI(I,J,3)+HLI(I,J,4))
c$$$      RLI(I,J)   = RLINEW
c$$$      MLI(I,J,2) = MLINEW
c$$$      HLI(I,J,3) = HLI(I,J,3) -  FHLI3
c$$$      HLI(I,J,4) = HLI(I,J,4) + (FHLI3 - FHLI4)
c$$$      GO TO 320
c$$$C**** Some lake ice has been frozen
c$$$  310 FREZ12 = (RLINEW-RLI(I,J))*ACE1I
c$$$      FREZ34 = (RLINEW-RLI(I,J))*MLINEW
c$$$      FREZ4  = RLI(I,J)*(MLINEW-MLI(I,J,2))
c$$$      FHLI3  = HLI(I,J,4)*(MLINEW-MLI(I,J,2))*(XLI3/XLI4) / MLI(I,J,2)
c$$$      IF(MLINEW.lt.MLI(I,J,2))
c$$$     *FHLI3  = HLI(I,J,3)*(MLINEW-MLI(I,J,2)) / MLI(I,J,2)
c$$$      MWL(I,J) = MWL(I,J) - FLAKE(I,J)*(FREZ12+FREZ34+FREZ4)*DXYP(J)
c$$$      GML(I,J) = GML(I,J) +
c$$$     +        FLAKE(I,J)*(FREZ12+FREZ34+FREZ4)*LHM*DXYP(J)
c$$$      IF (MWL(I,J).lt.0) THEN
c$$$        WRITE(99,*) "Lake frozen up:",I,J,MWL(I,J)
c$$$        MWL(I,J) = 0.
c$$$        GML(I,J) = 0.
c$$$      END IF
c$$$c     DJ(J,37) = DJ(J,37) + FLAKE(I,J)*(FREZ12+FREZ34+FREZ4)
c$$$c     DJ(J,46) = DJ(J,46) - FLAKE(I,J)*(FREZ12+FREZ34+FREZ4)*LHM
c$$$      MLI(I,J,1) = ((RLI(I,J)*MLI(I,J,1)+FREZ12) / RLINEW)
c$$$      MLI(I,J,2) =  MLINEW
c$$$      HLI(I,J,1) = (RLI(I,J)*HLI(I,J,1)-.5*FREZ12*LHM) / RLINEW
c$$$      HLI(I,J,2) = (RLI(I,J)*HLI(I,J,2)-.5*FREZ12*LHM) / RLINEW
c$$$      HLI(I,J,3) = (RLI(I,J)*(HLI(I,J,3)+FHLI3)- .5*FREZ34*LHM)/RLINEW
c$$$      HLI(I,J,4) = (RLI(I,J)*(HLI(I,J,4)-FHLI3)-(.5*FREZ34+FREZ4)*LHM)/
c$$$     /              RLINEW
c$$$      RLI(I,J)   =  RLINEW
c$$$C**** load GDATA with updated quantities
c$$$ 320  GDATA(I,J,1) = MLI(I,J,1)-ACE1I
c$$$      GDATA(I,J,3)=(HLI(I,J,1)/(XLI1*MLI(I,J,1))+LHM)/SHI
c$$$      GDATA(I,J,7)=(HLI(I,J,2)/(XLI2*MLI(I,J,1))+LHM)/SHI
c$$$C****
c$$$      RETURN
c$$$      END SUBROUTINE LKCLIM
