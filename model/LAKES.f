      MODULE LAKES
!@sum  LAKES subroutines for Lakes and Rivers
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0 (based on LB265)
      USE E001M12_COM, only : IM,JM
      IMPLICIT NONE
C****
C**** Changes from Model III: MO -> MWL (kg), G0M -> GML (J),
C****    GZM -> TLAKE (deg C)
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

      END MODULE LAKES

      SUBROUTINE init_LAKES(inilake)
!@sum  init_LAKES initiallises lake variables
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : rhow
      USE E001M12_COM, only : im,jm,flake,zatmo,dtsrc,flice,hlake
      USE OCEAN, only : tocean
      USE GEOM, only : dxyp,dxv,dyv,dxp,imaxj
      USE LAKES
      USE LAKES_COM
      USE FILEMANAGER

      IMPLICIT NONE
      LOGICAL inilake
!@var I,J,I72,IU,JU,ID,JD,IMAX loop variables
      INTEGER I,J,I72,IU,JU,ID,JD,IMAX,INM
      INTEGER iu_RVR  !@var iu_RVR unit number for river direction file
      CHARACTER TITLEI*80, CDIREC(IM,JM)*1
      REAL*8 SPMIN,SPMAX,SPEED0,SPEED,DZDH,DZDH1
C****
C**** LAKECB  MWL      Mass of water in lake (kg)
C****         GML      Liquid lake enthalpy (J)
C****         TLAKE    Temperature of lake surface (C)
C****         RLI      Horizontal ratio of lake ice to lake (1)
C****         MLI      Mass of sea ice (kg/m**2)
C****         HLI      Heat including latent of sea ice (J/m**2)
C****         HLAKE    Lake sill depth (m)
C****
C**** FIXDCB  FLAKE    Lake fraction (1)
C****

      IF (INILAKE) THEN
C**** Set lake variables from model block TOCEAN
C**** This is just an estimate for the initiallisation
        DO J=2,JM-1
          DO I=1,IM
            IF (FLAKE(I,J).gt.0) THEN
              TLAKE(I,J) = TOCEAN(1,I,J)
              MWL(I,J) = RHOW*HLAKE(I,J)*FLAKE(I,J)*DXYP(J)
              GML(I,J) = MWL(I,J)*MAX(TLAKE(I,J),4d0)
            ELSE
              TLAKE(I,J) = 0.
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
            DZDH  = (ZATMO(IU,JU)-ZATMO(ID,JD)) / DHORZ(IU,JU)
            SPEED = SPEED0*DZDH/DZDH1
            IF(SPEED.lt.SPMIN)  SPEED = SPMIN
            IF(SPEED.gt.SPMAX)  SPEED = SPMAX
            RATE(IU,JU) = DTsrc*SPEED/DHORZ(IU,JU)
          END IF
        END DO
      END DO
C****
C**** Set runoff temperature of glacial ice to be 0 (C)
C**** Is this necessary? ERUN0==0 already from FLICE
c      DO J=1,JM
c        DO I=1,IM
c          IF(FLICE(I,J).GT.0.)  TLAKE(I,J) = 0
c        END DO
c      END DO

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
      USE CONSTANT, only : grav,shw,rhow
      USE E001M12_COM, only : im,jm,flake,focean,zatmo,hlake
      USE GEOM, only : dxyp
      USE LAKES, only : kdirec,idpole,jdpole,rate,iflow,jflow
      USE LAKES_COM, only : tlake,gml,mwl
      USE FLUXES, only : flowo,eflowo
      USE DAGCOM, only : aij,ij_ervr,ij_mrvr

      IMPLICIT NONE
!@var I,J,IU,JU,ID,JD loop variables
      INTEGER I,J,IU,JU,ID,JD
c      REAL*8 XK(8),YK(8)
c      DATA XK/.707107d0,0.,-.707107d0,-1.,-.707107d0, 0., .707107d0,1./
c     *    ,YK/.707107d0,1., .707107d0, 0.,-.707107d0,-1.,-.707107d0,0./
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
            MWLSILL = RHOW*HLAKE(IU,JU)*FLAKE(IU,JU)*DXYP(JU)
            IF(MWL(IU,JU).gt.MWLSILL) THEN
              DMM = (MWL(IU,JU)-MWLSILL)*RATE(IU,JU)
C             DGM = GML(IU,JU)*DMM/MWL(IU,JU)   ! if well mixed 
              DGM = TLAKE(IU,JU)*DMM*SHW        ! surface only flow
              FLOW(IU,JU) =  FLOW(IU,JU) - DMM
              EFLOW(IU,JU) = EFLOW(IU,JU) - DGM
              AIJ(ID,JD,IJ_MRVR)=AIJ(ID,JD,IJ_MRVR) +  DMM
              AIJ(ID,JD,IJ_ERVR)=AIJ(ID,JD,IJ_ERVR)+
     *             (DGM+GRAV*DMM*(ZATMO(IU,JU)-ZATMO(ID,JD)))
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
      MWLSILL = RHOW*HLAKE(1,1)*FLAKE(1,1)*DXYP(1)
      IF(MWL(1,1).gt.MWLSILL) THEN
        DMM = (MWL(1,1)-MWLSILL)*RATE(1,1)
C       DGM = GML(1,1)*DMM/MWL(1,1)
        DGM = TLAKE(1,1)*DMM*SHW
        FLOW(1,1) =  FLOW(1,1) - DMM
        EFLOW(1,1) = EFLOW(1,1) - DGM
        AIJ(IDPOLE,JDPOLE,IJ_MRVR)=AIJ(IDPOLE,JDPOLE,IJ_MRVR) +  DMM
        AIJ(IDPOLE,JDPOLE,IJ_ERVR)=AIJ(IDPOLE,JDPOLE,IJ_ERVR) + (DGM
     *       +GRAV*DMM*(ZATMO(1,1)-ZATMO(IDPOLE,JDPOLE)))
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

      SUBROUTINE diag_RIVER
!@sum  diag_RIVER prints out the river outflow for various rivers
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : rhow
      USE E001M12_COM, only : jyear0,amon0,jdate0,jhour0,jyear,amon
     *     ,jdate,jhour,itime,dtsrc,idacc,itime0,nday 
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
      WRITE(6,902) JYEAR0,AMON0,JDATE0,JHOUR0,JYEAR,AMON,JDATE,JHOUR
     *     ,ITIME,DAYS
C**** convert kg to km^3/mon
      SCALERVR = 2.628d-2/(RHOW*DTSRC) 
      DO INM=1,NRVR,6
        DO I=1,6
          RVROUT(I) = SCALERVR*AIJ(IRVRMTH(I-1+INM),JRVRMTH(I-1+INM)
     *         ,IJ_MRVR)/IDACC(1)
        END DO
        WRITE(6,903) (NAMERVR(I-1+INM),RVROUT(I),I=1,6)
      END DO

 902  FORMAT ('0* River Outflow (km^3/mon) **  From:',I6,A6,I2,',  Hr'
     *     ,I3,6X,'To:',I6,A6,I2,', Hr',I3,'  Model-Time:',I9,5X
     *     ,'Dif:',F7.2,' Days')
 903  FORMAT (' ',A8,':',F8.3,5X,A8,':',F8.3,5X,A8,':',F8.3,5X,
     *            A8,':',F8.3,5X,A8,':',F8.3,5X,A8,':',F8.3)

      END SUBROUTINE diag_RIVER

      SUBROUTINE CHECKL (SUBRN)
!@sum  CHECKL checks whether the lake variables are reasonable.
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0 (based on LB265)
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
      USE E001M12_COM, only : im,jm,flake,kocean,itime,itimei
      USE GEOM, only : imaxj
!      USE LAKES, only : zimin,zimax,t_ice,t_noice,bydtmp
      USE LAKES_COM, only : t50
      USE OCEAN, only : dm
      USE SEAICE_COM, only : rsi,msi,tsi,snowi
      USE SEAICE, only : z1i
      IMPLICIT NONE
      REAL*8  ZIMIN,ZIMAX,T_ICE,T_NOICE,byDTMP
      INTEGER I,J,K,IMAX,IEND
      REAL*8 RSINEW

      IF (KOCEAN.eq.0) RETURN

      IF (IEND.eq.0 .and. ITime.gt.ITimeI) RETURN
      ZIMIN=.5d0                ! minimum lake ice thickness
      ZIMAX=2d0                 ! maximum lake ice thickness
      T_ICE = -8.  ! surface air temperature for 100% ice cover
      T_NOICE = 0. ! surface air temperature for no ice cover
      byDTMP = 1./(T_ICE-T_NOICE)
      DO J=1,JM
        IMAX=IMAXJ(J)
        DO I=1,IMAX
          IF (FLAKE(I,J) .GT. 0.) THEN ! linear fit for -8< T50 <0
            RSINEW = MIN(1d0,MAX(0d0,(T50(I,J)-T_NOICE)*byDTMP))
            RSI(I,J)=RSINEW
            MSI(I,J)=RHOI*(ZIMIN-Z1I+(ZIMAX-ZIMIN)*RSINEW*DM(I,J))
            IF (RSINEW.LE.0.) THEN
              SNOWI(I,J)=0.
              TSI(:,I,J)=0.
            END IF
          END IF
        END DO
      END DO

      RETURN
      END SUBROUTINE daily_LAKE

      SUBROUTINE PRECIP_LK
!@sum  PRECIP_LK driver for applying precipitation/melt to lake fraction
!@auth Original Development team
!@ver  1.0
      USE CONSTANT, only : rhow,shw
      USE E001M12_COM, only : im,jm,fland,flice,flake,kocean,hlake
      USE GEOM, only : imaxj,dxyp
      USE FLUXES, only : runosi,runoli,prec,eprec
      USE SEAICE, only : ace1i
      USE SEAICE_COM, only : rsi,msi,snowi
      USE LAKES_COM, only : mwl,gml,tlake
      USE DAGCOM, only : aj,cj,aij,j_eprcp,ij_f0oc,j_run2,j_dwtr2
      IMPLICIT NONE

      REAL*8 PRCP,ENRGP,DXYPJ,PLICE,PLKICE,RUN0,ERUN0,POLAKE,ROICE
     *     ,TGW,WTRO,ERUN4,ENRGO,SNOW,SMSI,ENGRW,WTRW0,WTRW,RUN4,ENRGW
      INTEGER I,J,IMAX

      DO J=1,JM
      IMAX=IMAXJ(J)
      DXYPJ=DXYP(J)
      DO I=1,IMAX
      POLAKE=(1.-RSI(I,J))*FLAKE(I,J)
      PLKICE=RSI(I,J)*FLAKE(I,J)
      PLICE=FLICE(I,J)
      ROICE=RSI(I,J)
      PRCP=PREC(I,J)
      ENRGP=EPREC(I,J)        ! energy of precipitation
      IF (FLAKE(I,J)+FLICE(I,J).gt.0) THEN

        IF (FLAKE(I,J).gt.0) THEN
          AJ(J,J_EPRCP)=AJ(J,J_EPRCP)+ENRGP*POLAKE
          AIJ(I,J,IJ_F0OC)=AIJ(I,J,IJ_F0OC)+ENRGP*POLAKE
        END IF

        RUN0 =POLAKE*PRCP  + PLKICE* RUNOSI(I,J) + PLICE* RUNOLI(I,J)
        ERUN0=POLAKE*ENRGP ! PLKICE*ERUNOSI(I,J) + PLICE*ERUNOLI(I,J) = 0
        
        MWL(I,J) = MWL(I,J) +  RUN0*DXYP(J)
        GML(I,J) = GML(I,J) + ERUN0*DXYP(J)

C**** This is here for continuity only
        IF (FLAKE(I,J).gt.0 .and. KOCEAN .EQ. 1) THEN
          TGW=TLAKE(I,J)
          WTRO=HLAKE(I,J)*RHOW   != MLD
          RUN0=RUNOSI(I,J)
          SNOW=SNOWI(I,J)
          SMSI=MSI(I,J)+ACE1I+SNOW
          RUN4=PRCP
          ERUN4=RUN4*TGW*SHW

          IF (PLKICE.LE.0.) THEN
            ENRGW=TGW*WTRO*SHW + ENRGP - ERUN4
            WTRW =WTRO
          ELSE
            WTRW0=WTRO -ROICE*SMSI
            ENRGW=WTRW0*TGW*SHW + (1.-ROICE)*ENRGP - ERUN4   
            WTRW =WTRW0+ROICE*(RUN0-RUN4)
          END IF
          TGW=ENRGW/(WTRW*SHW)
          TLAKE(I,J)=TGW
          AJ(J,J_RUN2) =AJ(J,J_RUN2) +RUN4 *POLAKE
          AJ(J,J_DWTR2)=AJ(J,J_DWTR2)+ERUN4*POLAKE
          CJ(J,J_RUN2) =CJ(J,J_RUN2) +RUN4 *PLKICE
          CJ(J,J_DWTR2)=CJ(J,J_DWTR2)+ERUN4*PLKICE
        END IF

      END IF
      END DO
      END DO
C****
      END SUBROUTINE PRECIP_LK

      SUBROUTINE GROUND_LK
!@sum  GROUND_LK driver for applying surface fluxes to lake fraction
!@auth Original Development Team
!@ver  1.0
!@calls 
      USE CONSTANT, only : twopi,rhow,shw,edpery,tf
      USE E001M12_COM, only : im,jm,flake,flice,kocean,fland,hlake
     *     ,fearth
      USE GEOM, only : imaxj,dxyp
      USE FLUXES, only : runosi, erunosi, e0,e1,evapor, dmsi,dhsi,
     *     runoli, runoe, erunoe
      USE OCEAN, only : osourc
      USE SEAICE_COM, only : rsi,msi,snowi
      USE SEAICE, only : ace1i
      USE PBLCOM, only : tsavg
      USE DAGCOM, only : aj,cj,aij,areg,jreg,j_eprcp,ij_f0oc,j_run2
     *     ,j_dwtr2,j_tg1,j_tg2,j_evap,j_oht,j_omlt,j_erun2,j_imelt
     *     ,ij_tgo,ij_tg1,ij_evap,ij_evapo
      USE LAKES_COM, only : mwl,gml,t50,tlake,tfl
      IMPLICIT NONE
C**** grid box variables
      REAL*8 POLAKE, PLKICE, PEARTH, PLICE, DXYPJ
C**** prognostic varaibles 
      REAL*8 TGW, WTRO, SMSI, ROICE
C**** fluxes 
      REAL*8 EVAPO, EVAPI, F2DT, F0DT, OTDT, RVRRUN, RVRERUN, RUN0,
     *     RUNLI, RUNE, ERUNE
C**** output from OSOURC
      REAL*8 ERUN4I, ERUN4O, RUN4I, RUN4O, ENRGFO, ACEFO, ACE2F, ENRGFI 

      REAL*8, PARAMETER :: FACT_T50  =1.-1./(24.*50.),
     *                     FACT_TSAVG=1./(24.*50.)
      INTEGER I,J,IMAX,JR

      DO J = 1,JM
      DO I = 1,IM
        T50(I,J)=T50(I,J)*FACT_T50+(TSAVG(I,J)-TF)*FACT_TSAVG
      END DO
      END DO

      DO J=1,JM
      IMAX=IMAXJ(J)
      DXYPJ=DXYP(J)
      DO I=1,IMAX
      JR=JREG(I,J)
      ROICE=RSI(I,J)
      PLKICE=FLAKE(I,J)*RSI(I,J)
      POLAKE=FLAKE(I,J)*(1.-RSI(I,J))

C**** Add land ice and surface runoff to lake variables
      IF (FLAND(I,J).gt.0) THEN
        PLICE=FLICE(I,J)
        PEARTH=FEARTH(I,J)
        RUNLI=RUNOLI(I,J)
        RUNE=RUNOE(I,J)
        ERUNE=ERUNOE(I,J)
        MWL(I,J) = MWL(I,J) + (RUNLI*PLICE + RUNE*PEARTH)*DXYPJ
        GML(I,J) = GML(I,J) + (             ERUNE*PEARTH)*DXYPJ
      END IF

      IF (FLAKE(I,J).gt.0) THEN

        TGW  =TLAKE(I,J)
        EVAPO=EVAPOR(I,J,1)
        EVAPI=EVAPOR(I,J,2) ! evaporation/dew at the ice surface (kg/m^2)
        F0DT =E0(I,J,1)
        SMSI =MSI(I,J)+ACE1I+SNOWI(I,J)
C**** get ice-ocean fluxes from sea ice routine
        RUN0=RUNOSI(I,J)        ! includes ACE2M term 
        F2DT=ERUNOSI(I,J)

        AJ(J,J_TG1) =AJ(J,J_TG1) +TGW  *POLAKE
        AJ(J,J_EVAP)=AJ(J,J_EVAP)+EVAPO*POLAKE
        IF (JR.ne.24) AREG(JR,J_TG1)=AREG(JR,J_TG1)+TGW*POLAKE*DXYPJ
        AIJ(I,J,IJ_TGO)  =AIJ(I,J,IJ_TGO)  +TGW
        AIJ(I,J,IJ_TG1)  =AIJ(I,J,IJ_TG1)  +TGW  *POLAKE
        AIJ(I,J,IJ_EVAP) =AIJ(I,J,IJ_EVAP) +EVAPO*POLAKE
        AIJ(I,J,IJ_EVAPO)=AIJ(I,J,IJ_EVAPO)+EVAPO*POLAKE
        
        MWL(I,J) = MWL(I,J) - EVAPO*POLAKE*DXYPJ
        GML(I,J) = GML(I,J) + F0DT *POLAKE*DXYPJ

        IF (KOCEAN .EQ. 1) THEN
          WTRO=HLAKE(I,J)*RHOW  

          OTDT=0.
          RVRRUN=0.
          RVRERUN=0.
C**** Calculate the amount of ice formation            
            CALL OSOURC (ROICE,SMSI,TGW,WTRO,OTDT,RUN0,F0DT,F2DT,RVRRUN
     *           ,RVRERUN,EVAPO,EVAPI,TFL,RUN4O,ERUN4O,RUN4I,ERUN4I
     *           ,ENRGFO,ACEFO,ACE2F,ENRGFI)

C**** Resave prognostic variables
          TLAKE(I,J)=TGW
C**** Open lake diagnostics 
          AJ(J,J_TG2)  =AJ(J,J_TG2)  +TLAKE(I,J)*POLAKE
          AJ(J,J_RUN2) =AJ(J,J_RUN2) +RUN4O     *POLAKE
          AJ(J,J_OMLT) =AJ(J,J_OMLT) +TLAKE(I,J)*POLAKE
          AJ(J,J_DWTR2)=AJ(J,J_DWTR2)+ERUN4O    *POLAKE
          IF (JR.ne.24) AREG(JR,J_TG2)=AREG(JR,J_TG2)+TLAKE(I,J)
     *         *POLAKE*DXYPJ
          AIJ(I,J,IJ_F0OC)=AIJ(I,J,IJ_F0OC)+F0DT*POLAKE
C**** Ice-covered ocean diagnostics
          AJ(J,J_ERUN2)=AJ(J,J_ERUN2)-ENRGFO*POLAKE
          AJ(J,J_IMELT)=AJ(J,J_IMELT)-ACEFO *POLAKE
          CJ(J,J_RUN2) =CJ(J,J_RUN2) +RUN4I *PLKICE
          CJ(J,J_DWTR2)=CJ(J,J_DWTR2)+ERUN4I*PLKICE
          CJ(J,J_ERUN2)=CJ(J,J_ERUN2)-ENRGFI*PLKICE
          CJ(J,J_IMELT)=CJ(J,J_IMELT)-ACE2F *PLKICE
        ELSE
          ACEFO=0 ; ACE2F=0. ; ENRGFO=0. ; ENRGFI=0.
          AJ(J,J_TG2)  =AJ(J,J_TG2)  +TGW   *POLAKE
          IF (JR.ne.24) AREG(JR,J_TG2)=AREG(JR,J_TG2)+TGW*POLAKE*DXYPJ
        END IF

C**** Add ice mass/energy fluxes to lakes
        MWL(I,J) = MWL(I,J) + ((RUN0-ACE2F) *PLKICE - ACEFO*POLAKE)
     *       *DXYPJ
        GML(I,J) = GML(I,J) + ((F2DT-ENRGFI)*PLKICE -ENRGFO*POLAKE)
     *       *DXYPJ

C**** Store mass and energy fluxes for formation of sea ice
        DMSI(1,I,J)=ACEFO
        DMSI(2,I,J)=ACE2F
        DHSI(1,I,J)=ENRGFO
        DHSI(2,I,J)=ENRGFI
      END IF
      END DO
      END DO
      END SUBROUTINE GROUND_LK

