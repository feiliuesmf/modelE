C**** PE001M12 E001M12 SOMTQ PB395M12
C**** OPT(3)
C**** semi-random cloud overlap (+snow age updates+computed opt.d+diagn)
C**** to be used with R99E or later radiation  routines.  carbon/2
C**** Constant pressure at L=LS1 and above (SIGE(LS1)=0., PLE(LS1)=PTOP)
C**** Using 5 harmonics for horizontal ocean heat transport, thinner ice
C**** Routines included:   COSZ0, RADIA, DRYCNV, SDRAG
*****
      SUBROUTINE COSZ0
C****
C**** calculates the Earth's zenith angle, weighted either
C**** by time or by sun light.
C****
      USE CONSTANT, only : twopi
      USE E001M12_COM
      USE GEOM, only : DLAT,DLON
      USE RADNCB, only : COSD,SIND
      IMPLICIT NONE

      REAL*8, DIMENSION(IM) :: LT1,LT2,SLT1,SLT2,S2LT1,S2LT2
      REAL*8, DIMENSION(IM,JM) :: COSZ,COSZA
      REAL*8, DIMENSION(IM) :: RI,SINI,COSI
      REAL*8, DIMENSION(JM) :: SINJ,COSJ
      COMMON/WORK5/LT1,LT2,SLT1,SLT2,S2LT1,S2LT2
C**** ZERO1 HAS TO EQUAL THE CUT-OFF VALUE FOR COSZ USED IN SOLAR
C**** COSZS WORKS CORRECTLY ONLY IF ZERO1 >> 1.D-3
      REAL*8 ZERO1
      DATA ZERO1/1.D-2/
      INTEGER I,J,L
      REAL*8 PHIS,PSHIS,CPHIS,SPHIS,PHIN,SPHIN,CPHIN,PHIM,S2DAWN,S2DUSK
     *     ,ECOSZ,ECOSQZ,CLT1,CLT2,ZERO2,CDUSK,DUSK,DAWN,SDUSK,SDAWN
     *     ,CJCD,SJSD,SR1,CR1,SR2,CR2,ROT1,ROT2,DROT

C**** COMPUTE THE AREA WEIGHTED LATITUDES AND THEIR SINES AND COSINES
      PHIS=-.25*TWOPI
      SPHIS=-1.
      CPHIS=0.
      DO 20 J=1,JM-1
      PHIN=DLAT*(J-.5*JM)
      SPHIN=SIN(PHIN)
      CPHIN=COS(PHIN)
      PHIM=(PHIN*SPHIN+CPHIN-PHIS*SPHIS-CPHIS)/(SPHIN-SPHIS)
      SINJ(J)=SIN(PHIM)
      COSJ(J)=COS(PHIM)
      PHIS=PHIN
      SPHIS=SPHIN
   20 CPHIS=CPHIN
      PHIN=.25*TWOPI
      SPHIN=1.
      CPHIN=0.
      PHIM=(PHIN*SPHIN+CPHIN-PHIS*SPHIS-CPHIS)/(SPHIN-SPHIS)
      SINJ(JM)=SIN(PHIM)
      COSJ(JM)=COS(PHIM)
C**** COMPUTE THE SINES AND COSINES OF LONGITUDE
      DO 40 I=1,IM
      RI(I)=DLON*(I-.5)
      SINI(I)=SIN(RI(I))
   40 COSI(I)=COS(RI(I))
      RETURN
C****
C****
      ENTRY COSZT (ROT1,ROT2,COSZ)
C****
C**** THIS ENTRY COMPUTES THE ZENITH ANGLE WEIGHTED BY DAYTIME
C**** HOURS FROM ROT1 TO ROT2, GREENWICH MEAN TIME IN RADIANS.  ROT1
C**** MUST BE BETWEEN 0 AND 2*PI.  ROT2 MUST BE BETWEEN ROT1 AND
C**** ROT1+2*PI.  I=1 MUST LIE ON THE INTERNATIONAL DATE LINE.
C****
      DROT=ROT2-ROT1
C**** COMPUTE THE SINES AND COSINES OF THE INITIAL AND FINAL GMT'S
  100 SR1=SIN(ROT1)
      CR1=COS(ROT1)
      SR2=SIN(ROT2)
      CR2=COS(ROT2)
C**** COMPUTE THE INITIAL AND FINAL LOCAL TIMES (MEASURED FROM NOON TO
C****   NOON) AND THEIR SINES AND COSINES
      DO 120 I=1,IM
      LT1(I)=ROT1+RI(I)
      SLT1(I)=SR1*COSI(I)+CR1*SINI(I)
      LT2(I)=ROT2+RI(I)
  120 SLT2(I)=SR2*COSI(I)+CR2*SINI(I)
C****
C**** CALCULATION FOR POLAR GRID BOXES
C****
      DO 200 J=1,JM,JM-1
      SJSD=SINJ(J)*SIND
      CJCD=COSJ(J)*COSD
      IF (SJSD+CJCD.LE.ZERO1) GO TO 180
      IF (SJSD-CJCD.GE.0.) GO TO 160
C**** AVERAGE COSZ FROM DAWN TO DUSK NEAR THE POLES
      DUSK=ACOS(-SJSD/CJCD)
      SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
      DAWN=-DUSK
      SDAWN=-SDUSK
      COSZ(1,J)=(SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN))/TWOPI
      GO TO 200
C**** CONSTANT DAYLIGHT NEAR THE POLES
  160 COSZ(1,J)=SJSD
      GO TO 200
C**** CONSTANT NIGHTIME NEAR THE POLES
  180 COSZ(1,J)=0.
  200 CONTINUE
C****
C**** LOOP OVER NON-POLAR LATITUDES
C****
      DO 500 J=2,JM-1
      SJSD=SINJ(J)*SIND
      CJCD=COSJ(J)*COSD
      IF (SJSD+CJCD.LE.ZERO1) GO TO 460
      IF (SJSD-CJCD.GE.0.) GO TO 420
C**** COMPUTE DAWN AND DUSK (AT LOCAL TIME) AND THEIR SINES
      DUSK=ACOS(-SJSD/CJCD)
      SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
      DAWN=-DUSK
      SDAWN=-SDUSK
C**** NEITHER CONSTANT DAYTIME NOR CONSTANT NIGHTIME AT THIS LATITUDE,
C**** LOOP OVER LONGITUDES
      ZERO2=ZERO1/CJCD
      DO 400 I=1,IM
C**** FORCE DUSK TO LIE BETWEEN LT1 AND LT1+2*PI
      IF (DUSK.GT.LT1(I)+ZERO2) GO TO 220
      DUSK=DUSK+TWOPI
      DAWN=DAWN+TWOPI
  220 IF (DAWN.LT.LT2(I)-ZERO2) GO TO 240
C**** CONTINUOUS NIGHTIME FROM INITIAL TO FINAL TIME
      COSZ(I,J)=0.
      GO TO 400
  240 IF (DAWN.GE.LT1(I)) GO TO 300
      IF (DUSK.LT.LT2(I)) GO TO 260
C**** CONTINUOUS DAYLIGHT FROM INITIAL TIME TO FINAL TIME
      COSZ(I,J)=SJSD+CJCD*(SLT2(I)-SLT1(I))/DROT
      GO TO 400
  260 IF (DAWN+TWOPI.LT.LT2(I)-ZERO2) GO TO 280
C**** DAYLIGHT AT INITIAL TIME AND NIGHT AT FINAL TIME
      COSZ(I,J)=(SJSD*(DUSK-LT1(I))+CJCD*(SDUSK-SLT1(I)))/DROT
      GO TO 400
C**** DAYLIGHT AT INITIAL AND FINAL TIMES WITH NIGHTIME IN BETWEEN
  280 COSZ(I,J)=(SJSD*(LT2(I)-DAWN-TWOPI+DUSK-LT1(I))+CJCD*
     *  (SLT2(I)-SDAWN+SDUSK-SLT1(I)))/DROT
      GO TO 400
  300 IF (DUSK.LT.LT2(I)) GO TO 320
C**** NIGHT AT INITIAL TIME AND DAYLIGHT AT FINAL TIME
      COSZ(I,J)=(SJSD*(LT2(I)-DAWN)+CJCD*(SLT2(I)-SDAWN))/DROT
      GO TO 400
C**** NIGHTIME AT INITIAL AND FINAL TIMES WITH DAYLIGHT IN BETWEEN
  320 COSZ(I,J)=(SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN))/DROT
  400 CONTINUE
      GO TO 500
C**** CONSTANT DAYLIGHT AT THIS LATITUDE
  420 DO 440 I=1,IM
  440 COSZ(I,J)=SJSD+CJCD*(SLT2(I)-SLT1(I))/DROT
      GO TO 500
C**** CONSTANT NIGHTIME AT THIS LATITUDE
  460 DO 480 I=1,IM
  480 COSZ(I,J)=0.
  500 CONTINUE
      RETURN
C****
C****
      ENTRY COSZS (ROT1,ROT2,COSZ,COSZA)
C****
C**** THIS ENTRY COMPUTES THE ZENITH ANGLE TWICE, FIRST WEIGHTED BY THE
C**** DAYTIME HOURS FROM ROT1 TO ROT2 AND SECONDLY WEIGHTED BY THE
C**** INCIDENT SUN LIGHT FROM ROT1 TO ROT2.  COSZT MUST HAVE BEEN
C**** CALLED JUST PREVIOUSLY.
C****
      DROT=ROT2-ROT1
C**** COMPUTE THE SINES AND COSINES OF THE INITIAL AND FINAL GMT'S
      SR1=SIN(ROT1)
      CR1=COS(ROT1)
      SR2=SIN(ROT2)
      CR2=COS(ROT2)
C**** COMPUTE THE INITIAL AND FINAL LOCAL TIMES (MEASURED FROM NOON TO
C****   NOON) AND THEIR SINES AND COSINES
      DO 520 I=1,IM
      LT1(I)=ROT1+RI(I)
      SLT1(I)=SR1*COSI(I)+CR1*SINI(I)
      CLT1=CR1*COSI(I)-SR1*SINI(I)
      S2LT1(I)=2.*SLT1(I)*CLT1
      LT2(I)=ROT2+RI(I)
      SLT2(I)=SR2*COSI(I)+CR2*SINI(I)
      CLT2=CR2*COSI(I)-SR2*SINI(I)
  520 S2LT2(I)=2.*SLT2(I)*CLT2
C****
C**** CALCULATION FOR POLAR GRID BOXES
C****
      DO 600 J=1,JM,JM-1
      SJSD=SINJ(J)*SIND
      CJCD=COSJ(J)*COSD
      IF (SJSD+CJCD.LE.ZERO1) GO TO 580
      IF (SJSD-CJCD.GE.0.) GO TO 560
C**** AVERAGE COSZ FROM DAWN TO DUSK NEAR THE POLES
      CDUSK=-SJSD/CJCD
      DUSK=ACOS(CDUSK)
      SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
      S2DUSK=2.*SDUSK*CDUSK
      DAWN=-DUSK
      SDAWN=-SDUSK
      S2DAWN=-S2DUSK
      ECOSZ=SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN)
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SDAWN)+
     *  .5*CJCD*(DUSK-DAWN+.5*(S2DUSK-S2DAWN)))
      COSZ(1,J)=ECOSZ/TWOPI
      COSZA(1,J)=ECOSQZ/ECOSZ
      GO TO 600
C**** CONSTANT DAYLIGHT NEAR THE POLES
  560 ECOSZ=SJSD*TWOPI
      ECOSQZ=SJSD*ECOSZ+.5*CJCD*CJCD*TWOPI
      COSZ(1,J)=ECOSZ/TWOPI
      COSZA(1,J)=ECOSQZ/ECOSZ
      GO TO 600
C**** CONSTANT NIGHTIME NEAR THE POLES
  580 COSZ(1,J)=0.
      COSZA(1,J)=0.
  600 CONTINUE
C****
C**** LOOP OVER NON-POLAR LATITUDES
C****
      DO 900 J=2,JM-1
      SJSD=SINJ(J)*SIND
      CJCD=COSJ(J)*COSD
      IF (SJSD+CJCD.LE.ZERO1) GO TO 860
      IF (SJSD-CJCD.GE.0.) GO TO 820
C**** COMPUTE DAWN AND DUSK (AT LOCAL TIME) AND THEIR SINES
      CDUSK=-SJSD/CJCD
      DUSK=ACOS(CDUSK)
      SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
      S2DUSK=2.*SDUSK*CDUSK
      DAWN=-DUSK
      SDAWN=-SDUSK
      S2DAWN=-S2DUSK
C**** NEITHER CONSTANT DAYTIME NOR CONSTANT NIGHTIME AT THIS LATITUDE,
C**** LOOP OVER LONGITUDES
      ZERO2=ZERO1/CJCD
      DO 800 I=1,IM
C**** FORCE DUSK TO LIE BETWEEN LT1 AND LT1+2*PI
      IF (DUSK.GT.LT1(I)+ZERO2) GO TO 620
      DUSK=DUSK+TWOPI
      DAWN=DAWN+TWOPI
  620 IF (DAWN.LT.LT2(I)-ZERO2) GO TO 640
C**** CONTINUOUS NIGHTIME FROM INITIAL TO FINAL TIME
      COSZ(I,J)=0.
      COSZA(I,J)=0.
      GO TO 800
  640 IF (DAWN.GE.LT1(I)) GO TO 700
      IF (DUSK.LT.LT2(I)) GO TO 660
C**** CONTINUOUS DAYLIGHT FROM INITIAL TIME TO FINAL TIME
      ECOSZ=SJSD*DROT+CJCD*(SLT2(I)-SLT1(I))
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SLT2(I)-SLT1(I))+
     *  .5*CJCD*(DROT+.5*(S2LT2(I)-S2LT1(I))))
      COSZ(I,J)=ECOSZ/DROT
      COSZA(I,J)=ECOSQZ/ECOSZ
      GO TO 800
  660 IF (DAWN+TWOPI.LT.LT2(I)-ZERO2) GO TO 680
C**** DAYLIGHT AT INITIAL TIME AND NIGHT AT FINAL TIME
      ECOSZ=SJSD*(DUSK-LT1(I))+CJCD*(SDUSK-SLT1(I))
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SLT1(I))+
     *  .5*CJCD*(DUSK-LT1(I)+.5*(S2DUSK-S2LT1(I))))
      COSZ(I,J)=ECOSZ/DROT
      COSZA(I,J)=ECOSQZ/ECOSZ
      GO TO 800
C**** DAYLIGHT AT INITIAL AND FINAL TIMES WITH NIGHTIME IN BETWEEN
  680 ECOSZ=SJSD*(DROT-DAWN-TWOPI+DUSK)+
     *  CJCD*(SLT2(I)-SDAWN+SDUSK-SLT1(I))
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SLT1(I)+SLT2(I)-SDAWN)+
     *  .5*CJCD*(DUSK+DROT-DAWN-TWOPI+
     *  .5*(S2DUSK-S2LT1(I)+S2LT2(I)-S2DAWN)))
      COSZ(I,J)=ECOSZ/DROT
      COSZA(I,J)=ECOSQZ/ECOSZ
      GO TO 800
  700 IF (DUSK.LT.LT2(I)) GO TO 720
C**** NIGHT AT INITIAL TIME AND DAYLIGHT AT FINAL TIME
      ECOSZ=SJSD*(LT2(I)-DAWN)+CJCD*(SLT2(I)-SDAWN)
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SLT2(I)-SDAWN)+
     *  .5*CJCD*(LT2(I)-DAWN+.5*(S2LT2(I)-S2DAWN)))
      COSZ(I,J)=ECOSZ/DROT
      COSZA(I,J)=ECOSQZ/ECOSZ
      GO TO 800
C**** NIGHTIME AT INITIAL AND FINAL TIMES WITH DAYLIGHT IN BETWEEN
  720 ECOSZ=SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN)
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SDAWN)+
     *  .5*CJCD*(DUSK-DAWN+.5*(S2DUSK-S2DAWN)))
      COSZ(I,J)=ECOSZ/DROT
      COSZA(I,J)=ECOSQZ/ECOSZ
  800 CONTINUE
      GO TO 900
C**** CONSTANT DAYLIGHT AT THIS LATITUDE
  820 DO 840 I=1,IM
      ECOSZ=SJSD*DROT+CJCD*(SLT2(I)-SLT1(I))
      ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SLT2(I)-SLT1(I))+
     *  .5*CJCD*(DROT+.5*(S2LT2(I)-S2LT1(I))))
      COSZ(I,J)=ECOSZ/DROT
  840 COSZA(I,J)=ECOSQZ/ECOSZ
      GO TO 900
C**** CONSTANT NIGHTIME AT THIS LATITUDE
  860 DO 880 I=1,IM
      COSZ(I,J)=0.
  880 COSZA(I,J)=0.
  900 CONTINUE
      RETURN
      END

      SUBROUTINE RADIA
!@sum  RADIA adds the radiation heating to the temperatures
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : grav,rgas,kapa,sday,lhe,lhs,twopi,tf,stbo
      USE E001M12_COM
      USE GEOM, Jdummy=>JLAT
      USE RADNCB, only : RQT,SRHR,TRHR,FSF,COSZ1,S0X,CO2,RSDIST
      USE RE001
     &  , only : setnew,rcomp1,writer,rcompx,rcompt ! routines
     &             ,lx
     &             ,FULGAS ,PTLISO ,KTREND ,LMR=>NL ,LMRP=>NLP
C     INPUT DATA
     &             ,PLE=>PLB ,TL=>TLM ,QL=>SHL
     &             ,TAUWC ,TAUIC ,SIZEWC ,SIZEIC
     &             ,POCEAN,PEARTH,POICE,PLICE,AGESN,SNOWE,SNOWOI,SNOWLI
     &             ,TGO,TGE,TGOI,TGLI,TS=>TSL,WS=>WMAG,WEARTH,PTOPTR
     &             ,S00WM2,RATLS0,S0,COSZ,PVT
     &             ,JYEARR=>JYEAR,JDAYR=>JDAY,JLAT,ILON,KCLDEM,KVEGA6
C     OUTPUT DATA
     &             ,TRDFLB ,TRNFLB ,TRFCRL
     &             ,SRDFLB ,SRNFLB ,SRFHRL
     &             ,PLAVIS ,PLANIR ,ALBVIS ,ALBNIR ,FSRNFG
     &             ,SRRVIS ,SRAVIS ,SRRNIR ,SRANIR
     &             ,BTEMPW
C     Some Aerosol Control Parameters not currently referenced
c    &             ,QXAERO ,QSAERO ,QCAERO ,ATAERO
c    &             ,QAER55 ,REAERO ,VEAERO ,ROAERO ,PI0MAX
c    &             ,FSAERO ,FTAERO ,VDGAER ,SSBTAU ,PIAERO
      USE RANDOM
      USE CLD01_COM_E001, only : tauss,taumc,svlhx,rhsav,svlat,cldsav,
     *     cldmc,cldss,csizmc,csizss
      USE PBLCOM, only : wsavg,tsavg
      USE DAGCOM, only : aj,bj,cj,areg,jreg,aij,ail,ajl,asjl,adaily,
     *     iwrite,jwrite,itwrite,ndlypt,j_pcldss,j_pcldmc,ij_pmccld,
     *     j_cdldep,j_pcld,ij_cldcv,ij_pcldl,ij_pcldm,ij_pcldh,
     *     ij_cldtppr,lm_req,j_srincp0,j_srnfp0,j_srnfp1,j_srincg,
     *     j_srnfg,j_brtemp,j_trincg,j_hsurf,j_hatm,j_plavis,ij_trnfp0,
     *     ij_srnfp0,ij_srincp0,ij_srnfg,ij_srincg,ij_btmpw,ij_srref
      USE DYNAMICS, only : pk,pedn,plij,pmid,pdsig
      USE OCEAN, only : tocean
      USE SEAICE_COM, only : rsi,snowi,tsi
      USE GHYCOM, only : snowe_com=>snowe,snoage,tearth,
     *     wearth_com=>wearth,aiearth
      USE LANDICE_COM, only : snowli_com=>snowli,tlandi
      USE FILEMANAGER

      IMPLICIT NONE

      REAL*8, DIMENSION(IM,JM) :: COSZ2,COSZA,TRINCG,BTMPW
      REAL*8, DIMENSION(4,IM,JM) :: SNFS,TNFS
      REAL*8, DIMENSION(LM_REQ,IM,JM) :: TRHRS,SRHRS
      REAL*8, DIMENSION(IM,JM,9) :: ALB

      REAL*8, DIMENSION(LM) :: TOTCLD

      REAL*8, DIMENSION(LM+3) :: COE

      INTEGER :: IFIRST = 1, JDLAST = -9
      INTEGER I,J,L,K,KR,LR,JYFIX,JDFIX,MADVEL,J50N,J70N,LLOW,LMID
     *     ,LHI,IHOUR,IMAX,IM1,JR,IH,INCH,LMID1,LHI1,JK
      REAL*8 COEX,CO2REF,ROT1,ROT2,PLAND,PIJ,RANDSS,RANDMC,CSS
     *     ,CMC,DEPTH,QSS,TAUSSL,TAUMCL,ELHX,CLDCV,DXYPJ,ASRHR,ATRHR
     *     ,ASNFS1,BSNFS1,CSNFS1,ATNFS1,BTNFS1,CTNFS1,SRNFLG,X
!@var NRFUN indices of unit numbers for radiation routines
      INTEGER NRFUN(14),IU
!@var RUNSTR names of files for radiation routines
      CHARACTER*5 :: RUNSTR(14) = (/"RADN1","RADN2","RADN3",
     *     "RADN4","RADN5","RADN6","RADN7","RADN8",
     *     "RADN9","RADNA","RADNB","RADNC","RADND",
     *     "RADNE"/)
!@var QBIN true if files for radiation input files are binary
      LOGICAL :: QBIN(14)=(/.TRUE.,.TRUE.,.FALSE.,.FALSE.,.TRUE.,.TRUE.
     *     ,.TRUE.,.TRUE.,.FALSE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE./)
C****
C**** FLAND     LAND COVERAGE (1)
C**** FLICE     LAND ICE COVERAGE (1)
C****
C**** TOCEAN(1)  OCEAN TEMPERATURE (C)
C****   RSI  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
C****
C**** GDATA  1  OCEAN ICE SNOW AMOUNT (KG/M**2)
C****        2  EARTH SNOW AMOUNT (KG/M**2)
C****        3  OCEAN ICE TEMPERATURE OF FIRST LAYER (C)
C****        4  EARTH TEMPERATURE OF FIRST LAYER (C)
C****        5  EARTH WATER OF FIRST LAYER (KG/M**2)
C****        6  EARTH ICE OF FIRST LAYER (KG/M**2)
C****        9  AGE OF SNOW OVER OCEAN ICE (DAYS)
C****       10  AGE OF SNOW OVER LAND ICE (DAYS)
C****       11  AGE OF SNOW OVER EARTH (DAYS)
C****       12  LAND ICE SNOW AMOUNT (KG/M**2)
C****       13  LAND ICE TEMPERATURE OF FIRST LAYER (C)
C****
C****
C**** VDATA  1-11 RATIOS FOR THE 11 VEGETATION TYPES (1)
C****
      REAL*8 QSAT
      IF (MODRD.EQ.0) IDACC(2)=IDACC(2)+1
      IF (IFIRST.NE.1) GO TO 50
      IFIRST=0
      CALL COSZ0
C**** SET THE CONTROL PARAMETERS FOR THE RADIATION (need mean pressures)
      LMR=LM+3
      LMRP=LMR+1
      COEX=.01*GRAV*KAPA/RGAS
      DO L=1,LM
         COE(L)=DTsrc*COEX/DSIG(L)
         PLE(L)=SIGE(L)*PSFMPT+PTOP
      END DO
      PLE(LM+1)=SIGE(LM+1)*PSFMPT+PTOP
      PLE(LM+2)=.5*PLE(LM+1)
      PLE(LMR)=.2*PLE(LM+1)
      PLE(LMR+1)=1.D-5
      PTOPTR=PTOP ! top of sigma-coord.system
      DO 40 LR=LM+1,LMR
   40 COE(LR)=DTsrc*NRAD*COEX/(PLE(LR)-PLE(LR+1))
      PTLISO=15.
      IF(CO2.LT.0.) KTREND=-NINT(CO2)
C**** Default: time-dependent So/GHG/O3/Trop-Aeros/Dust/Volc-Aeros
C****     For control runs e.g. with Jul 1,1951 atmosphere use
          JYFIX=1951
          JDFIX=182   !  Julian day (if JDFIX=0, annual cycle is used)
        CALL SETNEW(11,JYFIX,JDFIX, 1,0,0.D0) ! fix sol.const. - KSOLAR
        CALL SETNEW( 2,JYFIX,JDFIX, 0,0,0.D0) ! fix GHG (trend KTREND)
        CALL SETNEW(13,0    ,0    , 0,0,0.D0) ! no GHG-resets - MADGAS
        CALL SETNEW( 3,JYFIX,0    , 0,0,0.D0) ! quasi-fix O3 (ann.cycle)
        CALL SETNEW( 4,JYFIX,0    , 0,0,0.D0) ! quasi-fix Trop. Aerosols
        CALL SETNEW( 5,JYFIX,0    , 0,0,0.D0) ! quasi-fix Desert Dust
        CALL SETNEW( 6,JYFIX,JDFIX, 0,0,0.D0) ! fix Volc. Aerosols
C     An annual cycle is used for the data below, to prevent this, use
CnoAC   CALL SETNEW(7, 0,JDFIX, 0,0,0.D0) ! cloud heterogeneity - KCLDEP
CnoAC   CALL SETNEW(8, 0,JDFIX, 0,0,0.D0) !  surface albedo
C**** New options (currently not used)
      KCLDEM=0  ! 0:old 1:new LW cloud scattering scheme  -  KCLDEM
      KVEGA6=0  ! 0:2-band 1:6-band veg.albedo            -  KVEGA6
      MADVEL=123456         ! suppress reading i-th time series by i->0
C**** set up unit numbers for 14 radiation input files
      DO IU=1,14
         IF (IU.EQ.12.OR.IU.EQ.13) CYCLE    ! not used in GCM
         call getunit(RUNSTR(IU),NRFUN(IU),QBIN(IU),.true.)
      END DO
      CALL RCOMP1 (MADVEL,NRFUN) ! MAD 1-6: O3 TrAer Dust VAer Clds SoUV
      CO2REF=FULGAS(2)
      IF(CO2.GE.0.) FULGAS(2)=CO2REF*CO2
         CALL WRITER (6,0)
c         JEQ=1+JM/2
         J50N=(50.+90.)*(JM-1)/180.+1.5
         J70N=(70.+90.)*(JM-1)/180.+1.5
C**** CLOUD LAYER INDICES USED FOR DIAGNOSTICS
         DO 43 L=1,LM
         LLOW=L
         IF (.5*(PLE(L+1)+PLE(L+2)).LT.750.) GO TO 44 ! was 786. 4/16/97
   43    CONTINUE
   44    LMID1=LLOW+1
         DO 45 L=LMID1,LM
         LMID=L
         IF (.5*(PLE(L+1)+PLE(L+2)).LT.430.) GO TO 46
   45    CONTINUE
   46    LHI1=LMID+1
         LHI=LM
         IF (LHI1.GT.LHI) LHI=LHI1
         WRITE (6,47) LLOW,LMID1,LMID,LHI1,LHI
   47    FORMAT (' LOW CLOUDS IN LAYERS 1-',I2,'   MID LEVEL CLOUDS IN',
     *     ' LAYERS',I3,'-',I2,'   HIGH CLOUDS IN LAYERS',I3,'-',I2)
C**** Calculate mean cosine of zenith angle for the current physics step
   50 CONTINUE
      ROT1=(TWOPI*MOD(ITIME,NDAY))/NDAY  ! MOD(ITIME,NDAY)*TWOPI/NDAY ??
      ROT2=ROT1+TWOPI*DTsrc/SDAY
      CALL COSZT (ROT1,ROT2,COSZ1)
         IHOUR=1+JHOUR
      IF (MODRD.NE.0) GO TO 900
C****
C**** Interface with radiation routines, done only every NRAD time steps
C****
C**** Calculate mean cosine of zenith angle for the full radiation step
      ROT2=ROT1+TWOPI*NRAD*DTsrc/SDAY
      CALL COSZS (ROT1,ROT2,COSZ2,COSZA)
C****
C**** COMPUTE EARTH ALBEDOS AND OTHER PARAMETERS FOR BEGINNING OF DAY
C****
      JDAYR=JDAY
      JYEARR=JYEAR
      IF(JDAY.NE.JDLAST) CALL RCOMPT
      S0=S0X*S00WM2*RATLS0/RSDIST
      JDLAST=JDAY
C****
C**** MAIN J LOOP
C****
      DO 600 J=1,JM
      IMAX=IMAXJ(J)
      JLAT=NINT(1.+(J-1.)*45./(JM-1.))  !  j w.r.to 72x46 grid
C****
C**** MAIN I LOOP
C****
      IM1=IM
      DO 500 I=1,IMAX
         JR=JREG(I,J)
C**** DETERMINE FRACTIONS FOR SURFACE TYPES AND COLUMN PRESSURE
      PLAND=FLAND(I,J)
      POICE=RSI(I,J)*(1.-PLAND)
      POCEAN=(1.-PLAND)-POICE
      PLICE=FLICE(I,J)
      PEARTH=FEARTH(I,J)
C****
C**** DETERMINE CLOUDS (AND THEIR OPTICAL DEPTHS) SEEN BY RADIATION
C****
      RANDSS=RANDU(X)
      RANDMC=RANDU(X)
         CSS=0.
         CMC=0.
         DEPTH=0.
      DO 240 L=1,LM
      PIJ=PLIJ(L,I,J)
      QSS=Q(I,J,L)/(RHSAV(L,I,J)+1.D-20)
      QL(L)=QSS
      IF(CLDSAV(L,I,J).LT.1.)
     *  QL(L)=(Q(I,J,L)-QSS*CLDSAV(L,I,J))/(1.-CLDSAV(L,I,J))
      TL(L)=T(I,J,L)*PK(L,I,J)
      IF(CLDSS(L,I,J).EQ.0.) RANDSS=RANDU(X)
      TAUSSL=0.
      TAUMCL=0.
      TAUWC(L)=0.
      TAUIC(L)=0.
      SIZEWC(L)=CSIZMC(L,I,J)
      SIZEIC(L)=CSIZMC(L,I,J)
         TOTCLD(L)=0.
      IF (CLDSS(L,I,J).LT.RANDSS.OR.TAUSS(L,I,J).LE.0.) GO TO 220
      TAUSSL=TAUSS(L,I,J)
      QL(L)=QSS
         CSS=1.
         AJL(J,L,28)=AJL(J,L,28)+CSS
         TOTCLD(L)=1.
  220 IF (CLDMC(L,I,J).LT.RANDMC.OR.TAUMC(L,I,J).LE.0.) GO TO 230
         CMC=1.
         AJL(J,L,29)=AJL(J,L,29)+CMC
         TOTCLD(L)=1.
         DEPTH=DEPTH+PDSIG(L,I,J)
      IF(TAUMC(L,I,J).LE.TAUSSL) GO TO 230
      TAUMCL=TAUMC(L,I,J)
      ELHX=LHE
      IF(TL(L).LE.TF) ELHX=LHS
      QL(L)=QSAT(TL(L),ELHX,PMID(L,I,J))
  230    AJL(J,L,19)=AJL(J,L,19)+TOTCLD(L)
      IF(TAUSSL+TAUMCL.GT.0.) THEN
        IF(TAUMCL.GT.TAUSSL) THEN
          IF(SVLAT(L,I,J).EQ.LHE) THEN
            TAUWC(L)=TAUMCL
          ELSE
            TAUIC(L)=TAUMCL
          END IF
        ELSE
          IF(SVLHX(L,I,J).EQ.LHE) THEN
            TAUWC(L)=TAUSSL
            SIZEWC(L)=CSIZSS(L,I,J)
          ELSE
            TAUIC(L)=TAUSSL
            SIZEIC(L)=CSIZSS(L,I,J)
          END IF
        END IF
      END IF
  240 CONTINUE
         AJ(J,J_PCLDSS)=AJ(J,J_PCLDSS)+CSS*POCEAN
         BJ(J,J_PCLDSS)=BJ(J,J_PCLDSS)+CSS*PLAND
         CJ(J,J_PCLDSS)=CJ(J,J_PCLDSS)+CSS*POICE
         AREG(JR,J_PCLDSS)=AREG(JR,J_PCLDSS)+CSS*DXYP(J)
         AJ(J,J_PCLDMC)=AJ(J,J_PCLDMC)+CMC*POCEAN
         BJ(J,J_PCLDMC)=BJ(J,J_PCLDMC)+CMC*PLAND
         CJ(J,J_PCLDMC)=CJ(J,J_PCLDMC)+CMC*POICE
         AREG(JR,J_PCLDMC)=AREG(JR,J_PCLDMC)+CMC*DXYP(J)
         AIJ(I,J,IJ_PMCCLD)=AIJ(I,J,IJ_PMCCLD)+CMC
         AJ(J,J_CDLDEP)=AJ(J,J_CDLDEP)+DEPTH*POCEAN
         BJ(J,J_CDLDEP)=BJ(J,J_CDLDEP)+DEPTH*PLAND
         CJ(J,J_CDLDEP)=CJ(J,J_CDLDEP)+DEPTH*POICE
         AREG(JR,J_CDLDEP)=AREG(JR,J_CDLDEP)+DEPTH*DXYP(J)
         CLDCV=CMC+CSS-CMC*CSS
         AJ(J,J_PCLD)=AJ(J,J_PCLD)+CLDCV*POCEAN
         BJ(J,J_PCLD)=BJ(J,J_PCLD)+CLDCV*PLAND
         CJ(J,J_PCLD)=CJ(J,J_PCLD)+CLDCV*POICE
         AREG(JR,J_PCLD)=AREG(JR,J_PCLD)+CLDCV*DXYP(J)
         AIJ(I,J,IJ_CLDCV)=AIJ(I,J,IJ_CLDCV)+CLDCV
         DO 250 L=1,LLOW
         IF (TOTCLD(L).NE.1.) GO TO 250
         AIJ(I,J,IJ_PCLDL)=AIJ(I,J,IJ_PCLDL)+1.
         GO TO 255
  250    CONTINUE
  255    DO 260 L=LMID1,LMID
         IF (TOTCLD(L).NE.1.) GO TO 260
         AIJ(I,J,IJ_PCLDM)=AIJ(I,J,IJ_PCLDM)+1.
         GO TO 265
  260    CONTINUE
  265    DO 270 L=LHI1,LHI
         IF (TOTCLD(L).NE.1.) GO TO 270
         AIJ(I,J,IJ_PCLDH)=AIJ(I,J,IJ_PCLDH)+1.
         GO TO 275
  270    CONTINUE
  275    CONTINUE
         DO 280 L=LM,1,-1
         PIJ=PLIJ(L,I,J)
         IF (TOTCLD(L).NE.1.) GO TO 280
         AIJ(I,J,IJ_CLDTPPR)=AIJ(I,J,IJ_CLDTPPR)+SIGE(L+1)*PIJ+PTOP
         GO TO 285
  280    CONTINUE
  285    DO KR=1,NDLYPT
           IF (I.EQ.IJD6(1,KR).AND.J.EQ.IJD6(2,KR)) THEN
             IH=IHOUR
             DO INCH=1,NRAD
               IF (IH.GT.24) IH=IH-24
               ADAILY(IH,21,KR)=ADAILY(IH,21,KR)+TOTCLD(6)
               ADAILY(IH,22,KR)=ADAILY(IH,22,KR)+TOTCLD(5)
               ADAILY(IH,23,KR)=ADAILY(IH,23,KR)+TOTCLD(4)
               ADAILY(IH,24,KR)=ADAILY(IH,24,KR)+TOTCLD(3)
               ADAILY(IH,25,KR)=ADAILY(IH,25,KR)+TOTCLD(2)
               ADAILY(IH,26,KR)=ADAILY(IH,26,KR)+TOTCLD(1)
               ADAILY(IH,27,KR)=ADAILY(IH,27,KR)+CLDCV
               ADAILY(IH,53,KR)=ADAILY(IH,53,KR)+TOTCLD(7)
               ADAILY(IH,54,KR)=ADAILY(IH,54,KR)+TOTCLD(6)
               ADAILY(IH,55,KR)=ADAILY(IH,55,KR)+TOTCLD(5)
               ADAILY(IH,56,KR)=ADAILY(IH,56,KR)+TOTCLD(4)
               ADAILY(IH,57,KR)=ADAILY(IH,57,KR)+TOTCLD(3)
               ADAILY(IH,58,KR)=ADAILY(IH,58,KR)+TOTCLD(2)
               ADAILY(IH,59,KR)=ADAILY(IH,59,KR)+TOTCLD(1)
               ADAILY(IH,61,KR)=ADAILY(IH,61,KR)+CLDCV
               IH=IH+1
             END DO
           END IF
         END DO
C****
  300 CONTINUE
C****
C**** SET UP VERTICAL ARRAYS OMITTING THE I AND J INDICES
C****
C**** EVEN PRESSURES
      DO 340 L=1,LM
      PLE(L)=PEDN(L,I,J)
C**** TEMPERATURES
C---- TL(L)=T(I,J,L)*PK(L,I,J)     ! already defined
      IF(TL(L).LT.130..OR.TL(L).GT.370.) THEN
         WRITE(99,*) 'In Radia: Time,I,J,L,TL',ITime,I,J,L,TL(L)
         STOP 'In Radia: Temperature out of range'
      END IF
C**** MOISTURE VARIABLES
C---- QL(L)=Q(I,J,L)        ! already defined
  340 CONTINUE
C****
C**** RADIATION, SOLAR AND THERMAL
C****
      DO K=1,LM_REQ
        IF(RQT(K,I,J).LT.130..OR.RQT(K,I,J).GT.370.) THEN
        WRITE(99,*) 'In RADIA: Time,I,J,L,TL',ITime,I,J,LM+K,RQT(K,I,J)
        STOP 'In Radia: RQT out of range'
        END IF
        TL(LM+K)=RQT(K,I,J)
      END DO
      COSZ=COSZA(I,J)
      TGO =TOCEAN(1,I,J)+TF
      TGOI=TSI   (1,I,J)+TF
      TGLI=TLANDI(1,I,J)+TF
      TGE =TEARTH(  I,J)+TF
      TS=TSAVG(I,J)
      SNOWOI=SNOWI(I,J)
      SNOWLI=SNOWLI_COM(I,J)
      SNOWE=SNOWE_COM(I,J)
      AGESN(1)=SNOAGE(3,I,J)    ! land    ! why are these numbers so confusing?
      AGESN(2)=SNOAGE(1,I,J)    ! ocean ice
      AGESN(3)=SNOAGE(2,I,J)    ! land ice
      WEARTH=(WEARTH_COM(I,J)+AIEARTH(I,J))/(WFCS(I,J)+1.D-20)
      DO K=1,11
        PVT(K)=VDATA(I,J,K)
      END DO
      WS=WSAVG(I,J)
C-OLD FGOLDU(2)=XFRADJ*(1.-PEARTH)
C-OLD FGOLDU(3)=XFRADJ*PEARTH
      ILON=NINT(.5+(I-.5)*72./IM)
      CALL RCOMPX
      FSF(1,I,J)=FSRNFG(1)   !  ocean
      FSF(2,I,J)=FSRNFG(3)   !  ocean ice
      FSF(3,I,J)=FSRNFG(4)   !  land ice
      FSF(4,I,J)=FSRNFG(2)   !  soil
      IF(I.EQ.IWRITE.AND.J.EQ.JWRITE) CALL WRITER(6,ITWRITE)
      SRHR(1,I,J)=SRNFLB(1)
      TRHR(1,I,J)=STBO*(POCEAN*TGO**4+POICE*TGOI**4+PLICE*TGLI**4
     *  +PEARTH*TGE**4)-TRNFLB(1)
      DO 440 L=1,LM
      SRHR(L+1,I,J)=SRFHRL(L)
  440 TRHR(L+1,I,J)=-TRFCRL(L)
      DO 450 LR=1,LM_REQ
      SRHRS(LR,I,J)= SRFHRL(LM+LR)
  450 TRHRS(LR,I,J)=-TRFCRL(LM+LR)
      DO 460 K=1,4
      SNFS(K,I,J)=SRNFLB(K+LM)
  460 TNFS(K,I,J)=TRNFLB(K+LM)-TRNFLB(1)
         TRINCG(I,J)=TRDFLB(1)
         BTMPW(I,J)=BTEMPW-TF
         ALB(I,J,1)=SRNFLB(1)/(SRDFLB(1)+1.D-20)
         ALB(I,J,2)=PLAVIS
         ALB(I,J,3)=PLANIR
         ALB(I,J,4)=ALBVIS
         ALB(I,J,5)=ALBNIR
         ALB(I,J,6)=SRRVIS
         ALB(I,J,7)=SRRNIR
         ALB(I,J,8)=SRAVIS
         ALB(I,J,9)=SRANIR
  500 IM1=I
C****
C**** END OF MAIN LOOP FOR I INDEX
C****
  600 CONTINUE
C****
C**** END OF MAIN LOOP FOR J INDEX
C****
C****
C**** ACCUMULATE THE RADIATION DIAGNOSTICS
C****
         DO 780 J=1,JM
         DXYPJ=DXYP(J)
         IMAX=IMAXJ(J)
         DO L=1,LM
           ASRHR=0.
           ATRHR=0.
           DO I=1,IMAX
             ASRHR=ASRHR+SRHR(L+1,I,J)*COSZ2(I,J)
             ATRHR=ATRHR+TRHR(L+1,I,J)
           END DO
           AJL(J,L,9)=AJL(J,L,9)+ASRHR
           AJL(J,L,10)=AJL(J,L,10)+ATRHR
         END DO
         ASNFS1=0.
         BSNFS1=0.
         CSNFS1=0.
         ATNFS1=0.
         BTNFS1=0.
         CTNFS1=0.
         DO 770 I=1,IMAX
         COSZ=COSZ2(I,J)
         PLAND=FLAND(I,J)
         POICE=RSI(I,J)*(1.-PLAND)
         POCEAN=(1.-PLAND)-POICE
         JR=JREG(I,J)
         DO LR=1,LM_REQ
           ASJL(J,LR,3)=ASJL(J,LR,3)+SRHRS(LR,I,J)*COSZ
           ASJL(J,LR,4)=ASJL(J,LR,4)+TRHRS(LR,I,J)
         END DO
         DO KR=1,NDLYPT
           IF (I.EQ.IJD6(1,KR).AND.J.EQ.IJD6(2,KR)) THEN
             IH=IHOUR
             DO INCH=1,NRAD
               IF (IH.GT.24) IH=IH-24
               ADAILY(IH,2,KR)=ADAILY(IH,2,KR)+(1.-SNFS(4,I,J)/S0)
               ADAILY(IH,3,KR)=ADAILY(IH,3,KR)+(1.-ALB(I,J,1))
               ADAILY(IH,4,KR)=ADAILY(IH,4,KR)
     *              +((SNFS(4,I,J)-SNFS(1,I,J))*COSZ-TNFS(4,I,J)
     *              +TNFS(1,I,J))
               IH=IH+1
             END DO
           END IF
         END DO
  750    CONTINUE
         AJ(J,J_SRINCP0)=AJ(J,J_SRINCP0)+(S0*COSZ)*POCEAN
         BJ(J,J_SRINCP0)=BJ(J,J_SRINCP0)+(S0*COSZ)*PLAND
         CJ(J,J_SRINCP0)=CJ(J,J_SRINCP0)+(S0*COSZ)*POICE
         AREG(JR,J_SRINCP0)=AREG(JR,J_SRINCP0)+(S0*COSZ)*DXYPJ
         AJ(J,J_SRNFP0)=AJ(J,J_SRNFP0)+(SNFS(4,I,J)*COSZ)*POCEAN
         BJ(J,J_SRNFP0)=BJ(J,J_SRNFP0)+(SNFS(4,I,J)*COSZ)*PLAND
         CJ(J,J_SRNFP0)=CJ(J,J_SRNFP0)+(SNFS(4,I,J)*COSZ)*POICE
         AREG(JR,J_SRNFP0)=AREG(JR,J_SRNFP0)+(SNFS(4,I,J)*COSZ)*DXYPJ
         ASNFS1=ASNFS1+(SNFS(1,I,J)*COSZ)*POCEAN
         BSNFS1=BSNFS1+(SNFS(1,I,J)*COSZ)*PLAND
         CSNFS1=CSNFS1+(SNFS(1,I,J)*COSZ)*POICE
         AREG(JR,J_SRNFP1)=AREG(JR,J_SRNFP1)+(SNFS(1,I,J)*COSZ)*DXYPJ
         AJ(J,J_SRINCG)=AJ(J,J_SRINCG)+(SRHR(1,I,J)*COSZ/
     *        (ALB(I,J,1)+1.D-20))*POCEAN
         BJ(J,J_SRINCG)=BJ(J,J_SRINCG)+(SRHR(1,I,J)*COSZ/
     *        (ALB(I,J,1)+1.D-20))*PLAND
         CJ(J,J_SRINCG)=CJ(J,J_SRINCG)+(SRHR(1,I,J)*COSZ/
     *        (ALB(I,J,1)+1.D-20))*POICE
         AREG(JR,J_SRINCG)=AREG(JR,J_SRINCG)+
     *     (SRHR(1,I,J)*COSZ/(ALB(I,J,1)+1.D-20))*DXYPJ
         AJ(J,J_SRNFG)=AJ(J,J_SRNFG)+(FSF(1,I,J)*COSZ)*POCEAN
         SRNFLG=FSF(3,I,J)*FLICE(I,J)+FSF(4,I,J)*(PLAND-FLICE(I,J))
         BJ(J,J_SRNFG)=BJ(J,J_SRNFG)+(SRNFLG*COSZ)
         CJ(J,J_SRNFG)=CJ(J,J_SRNFG)+(FSF(2,I,J)*COSZ)*POICE
         AREG(JR,J_SRNFG)=AREG(JR,J_SRNFG)+(SRHR(1,I,J)*COSZ)*DXYPJ
         AJ(J,J_BRTEMP)=AJ(J,J_BRTEMP)+BTMPW(I,J)*POCEAN
         BJ(J,J_BRTEMP)=BJ(J,J_BRTEMP)+BTMPW(I,J)*PLAND
         CJ(J,J_BRTEMP)=CJ(J,J_BRTEMP)+BTMPW(I,J)*POICE
         AREG(JR,J_BRTEMP)=AREG(JR,J_BRTEMP)+BTMPW(I,J)*DXYPJ
         AJ(J,J_TRINCG)=AJ(J,J_TRINCG)+TRINCG(I,J)*POCEAN
         BJ(J,J_TRINCG)=BJ(J,J_TRINCG)+TRINCG(I,J)*PLAND
         CJ(J,J_TRINCG)=CJ(J,J_TRINCG)+TRINCG(I,J)*POICE
         AREG(JR,J_TRINCG)=AREG(JR,J_TRINCG)+TRINCG(I,J)*DXYPJ
         AJ(J,J_HSURF)=AJ(J,J_HSURF)-TNFS(4,I,J)*POCEAN
         BJ(J,J_HSURF)=BJ(J,J_HSURF)-TNFS(4,I,J)*PLAND
         CJ(J,J_HSURF)=CJ(J,J_HSURF)-TNFS(4,I,J)*POICE
         AREG(JR,J_HSURF)=AREG(JR,J_HSURF)-TNFS(4,I,J)*DXYPJ
         ATNFS1=ATNFS1-TNFS(1,I,J)*POCEAN
         BTNFS1=BTNFS1-TNFS(1,I,J)*PLAND
         CTNFS1=CTNFS1-TNFS(1,I,J)*POICE
         AREG(JR,J_HATM)=AREG(JR,J_HATM)-TNFS(1,I,J)*DXYPJ
         DO 760 K=2,9
           JK=K+J_PLAVIS-2     ! accumulate 8 radiation diags.
         AJ(J,JK)=AJ(J,JK)+(S0*COSZ)*ALB(I,J,K)*POCEAN
         BJ(J,JK)=BJ(J,JK)+(S0*COSZ)*ALB(I,J,K)*PLAND
         CJ(J,JK)=CJ(J,JK)+(S0*COSZ)*ALB(I,J,K)*POICE
  760    AREG(JR,JK)=AREG(JR,JK)+(S0*COSZ)*ALB(I,J,K)*DXYPJ
         AIJ(I,J,IJ_TRNFP0)=AIJ(I,J,IJ_TRNFP0)-TNFS(4,I,J)
         AIJ(I,J,IJ_SRNFP0)=AIJ(I,J,IJ_SRNFP0)+(SNFS(4,I,J)*COSZ)
         AIJ(I,J,IJ_SRINCP0)=AIJ(I,J,IJ_SRINCP0)+(S0*COSZ)
         AIJ(I,J,IJ_SRNFG)=AIJ(I,J,IJ_SRNFG)+(SRHR(1,I,J)*COSZ)
         AIJ(I,J,IJ_SRINCG)=AIJ(I,J,IJ_SRINCG)+(SRHR(1,I,J)*COSZ/
     *        (ALB(I,J,1)+1.D-20))
         AIJ(I,J,IJ_BTMPW)=AIJ(I,J,IJ_BTMPW)+BTMPW(I,J)
         AIJ(I,J,IJ_SRREF)=AIJ(I,J,IJ_SRREF)+S0*COSZ*ALB(I,J,2)
  770    CONTINUE
         AJ(J,J_SRNFP1)=AJ(J,J_SRNFP1)+ASNFS1
         BJ(J,J_SRNFP1)=BJ(J,J_SRNFP1)+BSNFS1
         CJ(J,J_SRNFP1)=CJ(J,J_SRNFP1)+CSNFS1
         AJ(J,J_HATM)=AJ(J,J_HATM)+ATNFS1
         BJ(J,J_HATM)=BJ(J,J_HATM)+BTNFS1
         CJ(J,J_HATM)=CJ(J,J_HATM)+CTNFS1
  780    CONTINUE
         DO 790 L=1,LM
         DO 790 I=1,IM
         AIL(I,L,7)=AIL(I,L,7)+((SRHR(L+1,I,JEQ-2)*COSZ2(I,JEQ-2)+
     *     TRHR(L+1,I,JEQ-2))*DXYP(JEQ-2)+(SRHR(L+1,I,JEQ-1)*
     *     COSZ2(I,JEQ-1)+TRHR(L+1,I,JEQ-1))*DXYP(JEQ-1)+
     *     (SRHR(L+1,I,JEQ)*COSZ2(I,JEQ)+TRHR(L+1,I,JEQ))*DXYP(JEQ))
         AIL(I,L,11)=AIL(I,L,11)+(SRHR(L+1,I,J50N)*COSZ2(I,J50N)+
     *     TRHR(L+1,I,J50N))*DXYP(J50N)
  790    AIL(I,L,15)=AIL(I,L,15)+(SRHR(L+1,I,J70N)*COSZ2(I,J70N)+
     *     TRHR(L+1,I,J70N))*DXYP(J70N)
C****
C**** Update radiative equilibrium temperatures
C****
      DO J=1,JM
        IMAX=IMAXJ(J)
        DO I=1,IMAX
          DO LR=1,LM_REQ
            RQT(LR,I,J)=RQT(LR,I,J)+(SRHRS(LR,I,J)*COSZ2(I,J)
     *           +TRHRS(LR,I,J))*COE(LR+LM)
          END DO
        END DO
      END DO
C****
C**** Update other temperatures every physics time step
C****
  900 DO J=1,JM
        IMAX=IMAXJ(J)
        DO I=1,IMAX
          DO L=1,LM
            T(I,J,L)=T(I,J,L)+(SRHR(L+1,I,J)*COSZ1(I,J)+TRHR(L+1,I,J))*
     *           COE(L)/(PLIJ(L,I,J)*PK(L,I,J))
          END DO
        END DO
      END DO
C**** daily diagnostics
      DO KR=1,NDLYPT
         ADAILY(IHOUR,1,KR)=ADAILY(IHOUR,1,KR)+S0*COSZ1(IJD6(1,KR)
     *        ,IJD6(2,KR))
      END DO

      RETURN
      END

      SUBROUTINE DRYCNV(LBASE_MIN,LBASE_MAX)
C****
C**** THIS SUBROUTINE MIXES AIR CAUSED BY DRY CONVECTION.
C**** THIS VERSION CHECKS BASE LAYERS LBASE_MIN TO LBASE_MAX.
C****
      USE CONSTANT, only : lhe,sha
      USE E001M12_COM
      USE GEOM
      USE QUSDEF, only : nmom,zmoms,xymoms
      USE SOMTQ_COM, only : tmom,qmom
      USE DAGCOM, only : ajl
      USE DYNAMICS, only : pk,pdsig,plij
      USE PBLCOM, only : dclev
      IMPLICIT NONE

      integer, intent(in) :: LBASE_MIN,LBASE_MAX
      REAL*8, DIMENSION(IM,JM,LM) :: UT,VT
      REAL*8, DIMENSION(LM) :: DP
      COMMON/WORK2/UT,VT,DP
      INTEGER, DIMENSION(IM) :: IDI,IDJ    !@var ID
      REAL*8, DIMENSION(IM) :: RA !@var
      REAL*8, DIMENSION(IM) :: UMS,VMS !@var
      LOGICAL POLE
      INTEGER I,J,L,K,IMAX,KMAX,IM1,LMAX,LMIN

      DOUBLE PRECISION, DIMENSION(NMOM) :: TMOMS,QMOMS
      REAL*8 RVX,DOK,PIJBOT,PIJ,PKMS,THPKMS,QMS
     *     ,TVMS,THETA,RDP,THM

      if(LBASE_MAX.GE.LM) stop 'DRYCNV: LBASE_MAX.GE.LM'
      RVX=0.
C**** LOAD U,V INTO UT,VT.  UT,VT WILL BE FIXED DURING DRY CONVECTION
C****   WHILE U,V WILL BE UPDATED.

      DO 50 L=1,LM
      DO 50 J=2,JM
      DO 50 I=1,IM
      UT(I,J,L)=U(I,J,L)
   50 VT(I,J,L)=V(I,J,L)
C**** OUTSIDE LOOPS OVER J AND I
      JLOOP: DO J=1,JM
      POLE=.FALSE.
      IF (J.EQ.1.OR.J.EQ.JM) POLE=.TRUE.

      IMAX=IMAXJ(J)
      KMAX=KMAXJ(J)
C****
C**** MAIN LOOP
C****
      IM1=IM
      ILOOP: DO I=1,IMAX
         DO K=1,KMAX
            RA(K)=RAJ(K,J)
            IDI(K)=IDIJ(K,I,J)
            IDJ(K)=IDJJ(K,J)
         END DO
      LMAX=LBASE_MIN-1
      lbase_loop: do while(lmax.lt.lbase_max)
      LMIN=LMAX+1
      LMAX=LMIN
      IF (T(I,J,LMIN)*(1.+Q(I,J,LMIN)*RVX).LE.
     *   T(I,J,LMIN+1)*(1.+Q(I,J,LMIN+1)*RVX)) cycle lbase_loop
C**** MIX HEAT AND MOISTURE THROUGHOUT THE UNSTABLE LAYERS
C**** MIX THROUGH TWO LOWER LAYERS
      PIJBOT=PLIJ(LMIN,I,J)
      DP(LMIN)=PDSIG(LMIN,I,J)
      PIJ=PLIJ(LMIN+1,I,J)
      DP(LMIN+1)=PDSIG(LMIN+1,I,J)
      PKMS=PK(LMIN,I,J)*DP(LMIN)+PK(LMIN+1,I,J)*DP(LMIN+1)
      THPKMS=T(I,J,LMIN)*(PK(LMIN,I,J)*DP(LMIN))
     *  +T(I,J,LMIN+1)*(PK(LMIN+1,I,J)*DP(LMIN+1))
      QMS=Q(I,J,LMIN)*DP(LMIN)+Q(I,J,LMIN+1)*DP(LMIN+1)
C**** sum moments to mix over unstable layers
      TMOMS(XYMOMS) =
     &     TMOM(XYMOMS,I,J,LMIN  )*(PK(LMIN  ,I,J)*DP(LMIN  ))  +
     &     TMOM(XYMOMS,I,J,LMIN+1)*(PK(LMIN+1,I,J)*DP(LMIN+1))
      QMOMS(XYMOMS) =
     &     QMOM(XYMOMS,I,J,LMIN  )*(DP(LMIN  ))  +
     &     QMOM(XYMOMS,I,J,LMIN+1)*(DP(LMIN+1))
      IF (LMIN+1.GE.LM) GO TO 150
      TVMS=T(I,J,LMIN)*(1.+Q(I,J,LMIN)*RVX)*(PK(LMIN,I,J)*DP(LMIN))
     *    +T(I,J,LMIN+1)*(1.+Q(I,J,LMIN+1)*RVX)
     *                                  *(PK(LMIN+1,I,J)*DP(LMIN+1))
      THETA=TVMS/PKMS
C**** MIX THROUGH SUBSEQUENT UNSTABLE LAYERS
      DO 140 L=LMIN+2,LM
      IF (THETA.LT.T(I,J,L)*(1.+Q(I,J,L)*RVX)) GO TO 160
      PIJ=PLIJ(L,I,J)
      DP(L)=PDSIG(L,I,J)
      PKMS=PKMS+(PK(L,I,J)*DP(L))
      THPKMS=THPKMS+T(I,J,L)*(PK(L,I,J)*DP(L))
      QMS=QMS+Q(I,J,L)*DP(L)
      TVMS=TVMS+T(I,J,L)*(1.+Q(I,J,L)*RVX)*(PK(L,I,J)*DP(L))
      TMOMS(XYMOMS) = TMOMS(XYMOMS) +
     &     TMOM(XYMOMS,I,J,L)*(PK(L,I,J)*DP(L))
      QMOMS(XYMOMS) = QMOMS(XYMOMS) +
     &     QMOM(XYMOMS,I,J,L)*DP(L)
  140 THETA=TVMS/PKMS
  150 L=LM+1
  160 LMAX=L-1
      RDP=1./(PIJBOT*SIGE(LMIN)-PIJ*SIGE(LMAX+1))
      THM=THPKMS/PKMS
      QMS=QMS*RDP
      DO 180 L=LMIN,LMAX
         AJL(J,L,12)=AJL(J,L,12)+(THM-T(I,J,L))*PK(L,I,J)*PLIJ(L,I,J)
         AJL(J,L,55)=AJL(J,L,55)+(QMS-Q(I,J,L))*PDSIG(L,I,J)*LHE/SHA
      T(I,J,L)=THM
      TMOM(XYMOMS,I,J,L)=TMOMS(XYMOMS)/PKMS
      TMOM(ZMOMS,I,J,L)=0.
      Q(I,J,L)=QMS
      QMOM(XYMOMS,I,J,L)=QMOMS(XYMOMS)*RDP
      QMOM(ZMOMS,I,J,L)=0.
  180 CONTINUE
C**** MIX MOMENTUM THROUGHOUT UNSTABLE LAYERS
      UMS(1:KMAX)=0.
      VMS(1:KMAX)=0.
      DO L=LMIN,LMAX
         DO K=1,KMAX
            UMS(K)=UMS(K)+UT(IDI(K),IDJ(K),L)*DP(L)
            VMS(K)=VMS(K)+VT(IDI(K),IDJ(K),L)*DP(L)
         ENDDO
      ENDDO
      UMS(1:KMAX)=UMS(1:KMAX)*RDP
      VMS(1:KMAX)=VMS(1:KMAX)*RDP
      DO L=LMIN,LMAX
         DO K=1,KMAX
            U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)
     &           +(UMS(K)-UT(IDI(K),IDJ(K),L))*RA(K)
            V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)
     &           +(VMS(K)-VT(IDI(K),IDJ(K),L))*RA(K)
c the following line gives bytewise different ajl
            AJL(IDJ(K),L,38)=AJL(IDJ(K),L,38)
     &           +(UMS(K)-UT(IDI(K),IDJ(K),L))*PLIJ(L,I,J)*RA(K)
         ENDDO
      ENDDO
      enddo lbase_loop
C**** ACCUMULATE BOUNDARY LAYER DIAGNOSTICS
      if(lbase_min.eq.1) then ! was called from surfce
         DCLEV(I,J)=LMAX
      endif
      IM1=I
      ENDDO ILOOP
      ENDDO JLOOP
      RETURN
      END

      SUBROUTINE SDRAG(DT1) 
C****
C**** THIS SUBROUTINE PUTS A DRAG ON THE WINDS ON THE TOP LAYER OF
C**** THE ATMOSPHERE
C****
      USE CONSTANT, only : grav,rgas
      USE E001M12_COM, only : im,jm,lm,psfmpt,u,v,sige,ptop,t,xcdlm
     *     ,bydsig,itime
c      USE GEOM
      USE DAGCOM, only : aij, ij_wlm,ajl,ij_sdrag
      USE DYNAMICS, only : pk
      IMPLICIT NONE

      INTEGER I,J,IP1,L
      REAL*8 PIJU,WL,RHO,CDN,X,BYPIJU
!@var DT1 time step (s)
      REAL*8, INTENT(IN) :: DT1

      PIJU = PSFMPT
      BYPIJU=1./PSFMPT
      DO L=LM,LM
      DO 100 J=2,JM
      I=IM
      DO 100 IP1=1,IM
      WL=SQRT(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))
      RHO=(PIJU*SIGE(L+1)+PTOP)/(RGAS*T(I,J,L)*PK(L,I,J))
      CDN=XCDLM(1)+XCDLM(2)*WL
         AIJ(I,J,IJ_WLM)=AIJ(I,J,IJ_WLM)+WL
      X=DT1*RHO*CDN*WL*GRAV*BYDSIG(L)*BYPIJU
      IF(X.GT.1) THEN
        write(99,*) 'SDRAG: ITime,I,J,PIJU,X,RHO,CDN,U(I,J,L),V(I,J,L)',
     *   ITime,I,J,PIJU,X,RHO,CDN,U(I,J,L),V(I,J,L),
     *   ' If problem persists, winds are too high! ',
     *   'Try setting XCDLM smaller.'
        X=1.
      END IF
         AJL(J,L,52) = AJL(J,L,52)-U(I,J,L)*X
c        IF(L.EQ.LM) AIJ(I,J,IJ_SDRAG)=AIJ(I,J,IJ_SDRAG)-U(I,J,L)*X
         AIJ(I,J,IJ_SDRAG)=AIJ(I,J,IJ_SDRAG)-U(I,J,L)*X
      U(I,J,L)=U(I,J,L)*(1.-X)
      V(I,J,L)=V(I,J,L)*(1.-X)
  100 I=IP1
      END DO
      RETURN
      END
