!@sum RAD_DRV contains drivers for the radiation related routines
!@ver  1.0
!@cont COSZ0, init_RAD, RADIA
C**** semi-random cloud overlap (computed opt.d+diagn)
C**** to be used with R99E or later radiation  routines.  carbon/2
C****
      SUBROUTINE COSZ0
!@sum  COSZ0 calculates Earth's zenith angle, weighted by time/sunlight
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : twopi
      USE MODEL_COM
      USE GEOM, only : dlat,dlon,lon,sinip,cosip
      USE RADNCB, only : cosd,sind,sinj,cosj
      IMPLICIT NONE
      SAVE
      REAL*8, DIMENSION(IM) :: LT1,LT2,SLT1,SLT2,S2LT1,S2LT2
      REAL*8, DIMENSION(IM,JM) :: COSZ,COSZA
      COMMON/WORK5/LT1,LT2,SLT1,SLT2,S2LT1,S2LT2
C**** ZERO1 HAS TO EQUAL THE CUT-OFF VALUE FOR COSZ USED IN SOLAR
C**** COSZS WORKS CORRECTLY ONLY IF ZERO1 >> 1.D-3
      REAL*8, PARAMETER :: ZERO1=1.D-2
      INTEGER I,J,L
      REAL*8 S2DAWN,S2DUSK,ECOSZ,ECOSQZ,CLT1,CLT2,ZERO2,CDUSK,DUSK,DAWN
     *     ,SDUSK,SDAWN,CJCD,SJSD,SR1,CR1,SR2,CR2,ROT1,ROT2,DROT
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
      SR1=SIN(ROT1)
      CR1=COS(ROT1)
      SR2=SIN(ROT2)
      CR2=COS(ROT2)
C**** COMPUTE THE INITIAL AND FINAL LOCAL TIMES (MEASURED FROM NOON TO
C****   NOON) AND THEIR SINES AND COSINES
      DO I=1,IM
        LT1(I)=ROT1+LON(I)
        SLT1(I)=SR1*COSIP(I)+CR1*SINIP(I)
        LT2(I)=ROT2+LON(I)
        SLT2(I)=SR2*COSIP(I)+CR2*SINIP(I)
      END DO
C****
C**** CALCULATION FOR POLAR GRID BOXES
C****
      DO J=1,JM,JM-1
        SJSD=SINJ(J)*SIND
        CJCD=COSJ(J)*COSD
        IF (SJSD+CJCD.GT.ZERO1) THEN
          IF (SJSD-CJCD.LT.0.) THEN
C**** AVERAGE COSZ FROM DAWN TO DUSK NEAR THE POLES
            DUSK=ACOS(-SJSD/CJCD)
            SDUSK=SQRT(CJCD*CJCD-SJSD*SJSD)/CJCD
            DAWN=-DUSK
            SDAWN=-SDUSK
            COSZ(1,J)=(SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN))/TWOPI
          ELSE
C**** CONSTANT DAYLIGHT NEAR THE POLES
            COSZ(1,J)=SJSD
          END IF
        ELSE
C**** CONSTANT NIGHTIME NEAR THE POLES
          COSZ(1,J)=0.
        END IF
      END DO
C****
C**** LOOP OVER NON-POLAR LATITUDES
C****
      DO 500 J=2,JM-1
      SJSD=SINJ(J)*SIND
      CJCD=COSJ(J)*COSD
      IF (SJSD+CJCD.GT.ZERO1) THEN
      IF (SJSD-CJCD.LT.0.) THEN
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
      ELSE
C**** CONSTANT DAYLIGHT AT THIS LATITUDE
        DO I=1,IM
          COSZ(I,J)=SJSD+CJCD*(SLT2(I)-SLT1(I))/DROT
        END DO
      END IF
      ELSE
C**** CONSTANT NIGHTIME AT THIS LATITUDE
        COSZ(1:IM,J)=0.
      END IF
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
      DO I=1,IM
        LT1(I)=ROT1+LON(I)
        SLT1(I)=SR1*COSIP(I)+CR1*SINIP(I)
        CLT1=CR1*COSIP(I)-SR1*SINIP(I)
        S2LT1(I)=2.*SLT1(I)*CLT1
        LT2(I)=ROT2+LON(I)
        SLT2(I)=SR2*COSIP(I)+CR2*SINIP(I)
        CLT2=CR2*COSIP(I)-SR2*SINIP(I)
        S2LT2(I)=2.*SLT2(I)*CLT2
      END DO
C****
C**** CALCULATION FOR POLAR GRID BOXES
C****
      DO J=1,JM,JM-1
        SJSD=SINJ(J)*SIND
        CJCD=COSJ(J)*COSD
        IF (SJSD+CJCD.GT.ZERO1) THEN
          IF (SJSD-CJCD.LT.0.) THEN
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
     *           .5*CJCD*(DUSK-DAWN+.5*(S2DUSK-S2DAWN)))
            COSZ(1,J)=ECOSZ/TWOPI
            COSZA(1,J)=ECOSQZ/ECOSZ
          ELSE
C**** CONSTANT DAYLIGHT NEAR THE POLES
            ECOSZ=SJSD*TWOPI
            ECOSQZ=SJSD*ECOSZ+.5*CJCD*CJCD*TWOPI
            COSZ(1,J)=ECOSZ/TWOPI
            COSZA(1,J)=ECOSQZ/ECOSZ
          END IF
        ELSE
C**** CONSTANT NIGHTIME NEAR THE POLES
          COSZ(1,J)=0.
          COSZA(1,J)=0.
        END IF
      END DO
C****
C**** LOOP OVER NON-POLAR LATITUDES
C****
      DO 900 J=2,JM-1
      SJSD=SINJ(J)*SIND
      CJCD=COSJ(J)*COSD
      IF (SJSD+CJCD.GT.ZERO1) THEN
      IF (SJSD-CJCD.LT.0.) THEN
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
  700 IF (DUSK.GE.LT2(I)) THEN
C**** NIGHT AT INITIAL TIME AND DAYLIGHT AT FINAL TIME
        ECOSZ=SJSD*(LT2(I)-DAWN)+CJCD*(SLT2(I)-SDAWN)
        ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SLT2(I)-SDAWN)+
     *       .5*CJCD*(LT2(I)-DAWN+.5*(S2LT2(I)-S2DAWN)))
        COSZ(I,J)=ECOSZ/DROT
        COSZA(I,J)=ECOSQZ/ECOSZ
      ELSE
C**** NIGHTIME AT INITIAL AND FINAL TIMES WITH DAYLIGHT IN BETWEEN
        ECOSZ=SJSD*(DUSK-DAWN)+CJCD*(SDUSK-SDAWN)
        ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SDUSK-SDAWN)+
     *       .5*CJCD*(DUSK-DAWN+.5*(S2DUSK-S2DAWN)))
        COSZ(I,J)=ECOSZ/DROT
        COSZA(I,J)=ECOSQZ/ECOSZ
      END IF
  800 CONTINUE
      ELSE
C**** CONSTANT DAYLIGHT AT THIS LATITUDE
        DO I=1,IM
          ECOSZ=SJSD*DROT+CJCD*(SLT2(I)-SLT1(I))
          ECOSQZ=SJSD*ECOSZ+CJCD*(SJSD*(SLT2(I)-SLT1(I))+
     *         .5*CJCD*(DROT+.5*(S2LT2(I)-S2LT1(I))))
          COSZ(I,J)=ECOSZ/DROT
          COSZA(I,J)=ECOSQZ/ECOSZ
        END DO
      END IF
C**** CONSTANT NIGHTIME AT THIS LATITUDE
      ELSE
        COSZ(1:IM,J)=0.
        COSZA(1:IM,J)=0.
      END IF
  900 CONTINUE
      RETURN
      END

      SUBROUTINE init_RAD
!@sum  init_RAD initialises radiation code
!@auth Original Development Team
!@ver  1.0
!@calls RE001:RCOMP1
      USE CONSTANT, only : grav,bysha,twopi
      USE MODEL_COM, only : jm,lm,dsig,sige,psfmpt,ptop,dtsrc,nrad
      USE GEOM, only : dlat
      USE RADNCB, only : s0x,co2,lm_req,llow,lmid,lhi,coe,sinj,cosj
      USE RE001, only : setnew,rcomp1,writer             ! routines
     &     ,FULGAS ,PTLISO ,KTREND ,LMR=>NL ,LMRP=>NLP, PLE=>PLB, PTOPTR
     *     ,KCLDEM,KVEGA6
      USE FILEMANAGER
      USE PARAM
      IMPLICIT NONE

      INTEGER J,L,LR,JYFIX,JDFIX,MADVEL
      REAL*8 COEX,CO2REF,SPHIS,CPHIS,PHIN,SPHIN,CPHIN,PHIM,PHIS
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

C**** sync radiation parameters from input
      call sync_param( "S0X", S0X )
      call sync_param( "CO2", CO2 )

C**** COMPUTE THE AREA WEIGHTED LATITUDES AND THEIR SINES AND COSINES
      PHIS=-.25*TWOPI
      SPHIS=-1.
      CPHIS=0.
      DO J=1,JM-1
        PHIN=DLAT*(J-.5*JM)
        SPHIN=SIN(PHIN)
        CPHIN=COS(PHIN)
        PHIM=(PHIN*SPHIN+CPHIN-PHIS*SPHIS-CPHIS)/(SPHIN-SPHIS)
        SINJ(J)=SIN(PHIM)
        COSJ(J)=COS(PHIM)
        PHIS=PHIN
        SPHIS=SPHIN
        CPHIS=CPHIN
      END DO
      PHIN=.25*TWOPI
      SPHIN=1.
      CPHIN=0.
      PHIM=(PHIN*SPHIN+CPHIN-PHIS*SPHIS-CPHIS)/(SPHIN-SPHIS)
      SINJ(JM)=SIN(PHIM)
      COSJ(JM)=COS(PHIM)
C****
C**** SET THE CONTROL PARAMETERS FOR THE RADIATION (need mean pressures)
C****
      LMR=LM+LM_REQ
      LMRP=LMR+1
      COEX=1d-2*GRAV*BYSHA
      DO L=1,LM
        COE(L)=DTsrc*COEX/DSIG(L)
        PLE(L)=SIGE(L)*PSFMPT+PTOP
      END DO
      PLE(LM+1)=SIGE(LM+1)*PSFMPT+PTOP
      PLE(LM+2)=.5*PLE(LM+1)
      PLE(LMR)=.2d0*PLE(LM+1)
      PLE(LMR+1)=1d-5
      PTOPTR=PTOP ! top of sigma-coord.system
      DO LR=LM+1,LMR
        COE(LR)=DTsrc*NRAD*COEX/(PLE(LR)-PLE(LR+1))
      END DO
      PTLISO=15.
      IF(CO2.LT.0.) KTREND=-NINT(CO2)
C**** Default: time-dependent So/GHG/O3/Trop-Aeros/Dust/Volc-Aeros
C****     For control runs e.g. with Jul 1,1951 atmosphere use
      JYFIX=1951
      JDFIX=182       !  Julian day (if JDFIX=0, annual cycle is used)
      CALL SETNEW(11,JYFIX,JDFIX, 1,0,0.D0) ! fix sol.const. - KSOLAR
      CALL SETNEW( 2,JYFIX,JDFIX, 0,0,0.D0) ! fix GHG (trend KTREND)
      CALL SETNEW(13,0    ,0    , 0,0,0.D0) ! no GHG-resets - MADGAS
      CALL SETNEW( 3,JYFIX,0    , 0,0,0.D0) ! quasi-fix O3 (ann.cycle)
      CALL SETNEW( 4,JYFIX,0    , 0,0,0.D0) ! quasi-fix Trop. Aerosols
      CALL SETNEW( 5,JYFIX,0    , 0,0,0.D0) ! quasi-fix Desert Dust
      CALL SETNEW( 6,JYFIX,JDFIX, 0,0,0.D0) ! fix Volc. Aerosols
C     An annual cycle is used for the data below, to prevent this, use
CnoAC CALL SETNEW(7, 0,JDFIX, 0,0,0.D0) ! cloud heterogeneity - KCLDEP
CnoAC CALL SETNEW(8, 0,JDFIX, 0,0,0.D0) !  surface albedo
C**** New options (currently not used)
      KCLDEM=0  ! 0:old 1:new LW cloud scattering scheme  -  KCLDEM
      KVEGA6=0  ! 0:2-band 1:6-band veg.albedo            -  KVEGA6
      MADVEL=123456         ! suppress reading i-th time series by i->0
C**** set up unit numbers for 14 radiation input files
      DO IU=1,14
        IF (IU.EQ.12.OR.IU.EQ.13) CYCLE ! not used in GCM
        call openunit(RUNSTR(IU),NRFUN(IU),QBIN(IU),.true.)
      END DO
      CALL RCOMP1 (MADVEL,NRFUN) ! MAD 1-6: O3 TrAer Dust VAer Clds SoUV
      DO IU=1,14
        IF (IU.EQ.12.OR.IU.EQ.13) CYCLE ! not used in GCM
        call closeunit(NRFUN(IU))
      END DO
      CO2REF=FULGAS(2)
      IF(CO2.GE.0.) FULGAS(2)=CO2REF*CO2
C**** write out now occurs after first pass to ensure correct ppm values
c      CALL WRITER (6,0)
C**** CLOUD LAYER INDICES USED FOR DIAGNOSTICS
      DO L=1,LM
        LLOW=L
        IF (.5*(PLE(L+1)+PLE(L+2)).LT.750.) GO TO 44 ! was 786. 4/16/97
      END DO
 44   DO L=LLOW+1,LM
        LMID=L
        IF (.5*(PLE(L+1)+PLE(L+2)).LT.430.) GO TO 46
      END DO
 46   LHI=LM
      IF (LMID+1.GT.LHI) LHI=LMID+1
      WRITE (6,47) LLOW,LLOW+1,LMID,LMID+1,LHI
 47   FORMAT (' LOW CLOUDS IN LAYERS 1-',I2,'   MID LEVEL CLOUDS IN',
     *     ' LAYERS',I3,'-',I2,'   HIGH CLOUDS IN LAYERS',I3,'-',I2)
C****
      END SUBROUTINE init_RAD

      SUBROUTINE RADIA
!@sum  RADIA adds the radiation heating to the temperatures
!@auth Original Development Team
!@ver  1.0
!@calls tropwmo, RE001:rcompt, RE001:rcompx, RE001:writer, coszs, coszt
      USE CONSTANT, only : sday,lhe,lhs,twopi,tf,stbo
      USE MODEL_COM
      USE GEOM
      USE RADNCB, only : rqt,srhr,trhr,fsf,cosz1,s0x,co2,rsdist,lm_req
     *     ,llow,lmid,lhi,coe
      USE RE001
     &  , only : writer,rcompx,rcompt ! routines
C     INPUT DATA
     &             ,PLE=>PLB ,TL=>TLM ,QL=>SHL
     &             ,TAUWC ,TAUIC ,SIZEWC ,SIZEIC
     &             ,POCEAN,PEARTH,POICE,PLICE,AGESN,SNOWE,SNOWOI,SNOWLI
     &             ,TGO,TGE,TGOI,TGLI,TS=>TSL,WS=>WMAG,WEARTH
     &             ,S00WM2,RATLS0,S0,COSZ,PVT
     &             ,JYEARR=>JYEAR,JDAYR=>JDAY,JLAT,ILON
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
      USE CLOUDS_COM, only : tauss,taumc,svlhx,rhsav,svlat,cldsav,
     *     cldmc,cldss,csizmc,csizss
      USE PBLCOM, only : wsavg,tsavg
      USE DAGCOM, only : aj,areg,jreg,aij,ail,ajl,asjl,adiurn,
     *     iwrite,jwrite,itwrite,ndiupt,j_pcldss,j_pcldmc,ij_pmccld,
     *     j_cdldep,j_pcld,ij_cldcv,ij_pcldl,ij_pcldm,ij_pcldh,
     *     ij_cldtppr,lm_req,j_srincp0,j_srnfp0,j_srnfp1,j_srincg,
     *     j_srnfg,j_brtemp,j_trincg,j_hsurf,j_hatm,j_plavis,ij_trnfp0,
     *     ij_srnfp0,ij_srincp0,ij_srnfg,ij_srincg,ij_btmpw,ij_srref,
     *     j50n,j70n,j_clrtoa,j_clrtrp,j_tottrp,il_req,il_r50n,il_r70n,
     *     ijdd,idd_cl7,idd_cl6,idd_cl5,idd_cl4,idd_cl3,idd_cl2,idd_cl1,
     *     idd_ccv,idd_isw,idd_palb,idd_galb,idd_absa,j5s,j5n,
     &     jl_srhr,jl_trcr,jl_totcld,jl_sscld,jl_mccld
      USE DYNAMICS, only : pk,pedn,plij,pmid,pdsig,ltropo
      USE SEAICE_COM, only : rsi,snowi
      USE GHYCOM, only : snowe_com=>snowe,snoage,wearth_com=>wearth
     *     ,aiearth
      USE LANDICE_COM, only : snowli_com=>snowli
      USE LAKES_COM, only : flake
      USE FLUXES, only : gtemp
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,JM) :: COSZ2,COSZA,TRINCG,BTMPW
      REAL*8, DIMENSION(4,IM,JM) :: SNFS,TNFS
      REAL*8, DIMENSION(LM_REQ,IM,JM) :: TRHRS,SRHRS
      REAL*8, DIMENSION(IM,JM,9) :: ALB
      REAL*8, DIMENSION(LM) :: TOTCLD

      INTEGER, SAVE :: JDLAST = -9, IFIRST = 1
      INTEGER I,J,L,K,KR,LR,IMAX,IM1,JR,IH,INCH,JK,IT
      REAL*8 ROT1,ROT2,PLAND,PIJ,RANDSS,RANDMC,CSS,CMC,DEPTH,QSS,TAUSSL
     *     ,TAUMCL,ELHX,CLDCV,DXYPJ,SRNFLG,X,OPNSKY
      REAL*8 QSAT
C****
C**** FLAND     LAND COVERAGE (1)
C**** FLICE     LAND ICE COVERAGE (1)
C****
C**** GTEMP(1)  GROUND TEMPERATURE ARRAY OVER ALL SURFACE TYPES (C)
C****   RSI  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
C****
C**** VDATA  1-11 RATIOS FOR THE 11 VEGETATION TYPES (1)
C****

C**** Calculate mean cosine of zenith angle for the current physics step
      ROT1=(TWOPI*MOD(ITIME,NDAY))/NDAY  ! MOD(ITIME,NDAY)*TWOPI/NDAY ??
      ROT2=ROT1+TWOPI*DTsrc/SDAY
      CALL COSZT (ROT1,ROT2,COSZ1)
      IF (MODRD.NE.0) GO TO 900
      IDACC(2)=IDACC(2)+1
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
      DO I=1,IMAX
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
         AJL(J,L,JL_SSCLD)=AJL(J,L,JL_SSCLD)+CSS
         TOTCLD(L)=1.
  220 IF (CLDMC(L,I,J).LT.RANDMC.OR.TAUMC(L,I,J).LE.0.) GO TO 230
         CMC=1.
         AJL(J,L,JL_MCCLD)=AJL(J,L,JL_MCCLD)+CMC
         TOTCLD(L)=1.
         DEPTH=DEPTH+PDSIG(L,I,J)
      IF(TAUMC(L,I,J).LE.TAUSSL) GO TO 230
      TAUMCL=TAUMC(L,I,J)
      ELHX=LHE
      IF(TL(L).LE.TF) ELHX=LHS
      QL(L)=QSAT(TL(L),ELHX,PMID(L,I,J))
  230    AJL(J,L,JL_TOTCLD)=AJL(J,L,JL_TOTCLD)+TOTCLD(L)
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
         CLDCV=CMC+CSS-CMC*CSS
         OPNSKY=1.-CLDCV
         DO IT=1,NTYPE
           AJ(J,J_PCLDSS,IT)=AJ(J,J_PCLDSS,IT)+CSS  *FTYPE(IT,I,J)
           AJ(J,J_PCLDMC,IT)=AJ(J,J_PCLDMC,IT)+CMC  *FTYPE(IT,I,J)
           AJ(J,J_CDLDEP,IT)=AJ(J,J_CDLDEP,IT)+DEPTH*FTYPE(IT,I,J)
           AJ(J,J_PCLD  ,IT)=AJ(J,J_PCLD  ,IT)+CLDCV*FTYPE(IT,I,J)
         END DO
         AREG(JR,J_PCLDSS)=AREG(JR,J_PCLDSS)+CSS  *DXYP(J)
         AREG(JR,J_PCLDMC)=AREG(JR,J_PCLDMC)+CMC  *DXYP(J)
         AREG(JR,J_CDLDEP)=AREG(JR,J_CDLDEP)+DEPTH*DXYP(J)
         AREG(JR,J_PCLD)  =AREG(JR,J_PCLD)  +CLDCV*DXYP(J)
         AIJ(I,J,IJ_PMCCLD)=AIJ(I,J,IJ_PMCCLD)+CMC
         AIJ(I,J,IJ_CLDCV) =AIJ(I,J,IJ_CLDCV) +CLDCV
         DO 250 L=1,LLOW
         IF (TOTCLD(L).NE.1.) GO TO 250
         AIJ(I,J,IJ_PCLDL)=AIJ(I,J,IJ_PCLDL)+1.
         GO TO 255
  250    CONTINUE
  255    DO 260 L=LLOW+1,LMID
         IF (TOTCLD(L).NE.1.) GO TO 260
         AIJ(I,J,IJ_PCLDM)=AIJ(I,J,IJ_PCLDM)+1.
         GO TO 265
  260    CONTINUE
  265    DO 270 L=LMID+1,LHI
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
  285    DO KR=1,NDIUPT
           IF (I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
             DO INCH=1,NRAD
               IH=1+MOD(JHOUR+INCH-1,24)
               ADIURN(IH,IDD_CL7,KR)=ADIURN(IH,IDD_CL7,KR)+TOTCLD(7)
               ADIURN(IH,IDD_CL6,KR)=ADIURN(IH,IDD_CL6,KR)+TOTCLD(6)
               ADIURN(IH,IDD_CL5,KR)=ADIURN(IH,IDD_CL5,KR)+TOTCLD(5)
               ADIURN(IH,IDD_CL4,KR)=ADIURN(IH,IDD_CL4,KR)+TOTCLD(4)
               ADIURN(IH,IDD_CL3,KR)=ADIURN(IH,IDD_CL3,KR)+TOTCLD(3)
               ADIURN(IH,IDD_CL2,KR)=ADIURN(IH,IDD_CL2,KR)+TOTCLD(2)
               ADIURN(IH,IDD_CL1,KR)=ADIURN(IH,IDD_CL1,KR)+TOTCLD(1)
               ADIURN(IH,IDD_CCV,KR)=ADIURN(IH,IDD_CCV,KR)+CLDCV
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
         WRITE(99,*) 'GTEMP:',GTEMP(:,:,I,J)
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
      TGO =GTEMP(1,1,I,J)+TF
      TGOI=GTEMP(1,2,I,J)+TF
      TGLI=GTEMP(1,3,I,J)+TF
      TGE =GTEMP(1,4,I,J)+TF
      TS=TSAVG(I,J)
      SNOWOI=SNOWI(I,J)
      SNOWLI=SNOWLI_COM(I,J)
      SNOWE=SNOWE_COM(I,J)
      AGESN(1)=SNOAGE(3,I,J)    ! land         ! ? why are these numbers
      AGESN(2)=SNOAGE(1,I,J)    ! ocean ice        so confusing ?
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
      DO L=1,LM
        SRHR(L+1,I,J)=SRFHRL(L)
        TRHR(L+1,I,J)=-TRFCRL(L)
      END DO
      DO LR=1,LM_REQ
        SRHRS(LR,I,J)= SRFHRL(LM+LR)
        TRHRS(LR,I,J)=-TRFCRL(LM+LR)
      END DO
      DO K=1,4
        SNFS(K,I,J)=SRNFLB(K+LM)
        TNFS(K,I,J)=TRNFLB(K+LM)-TRNFLB(1)
      END DO
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
C**** Save clear sky/tropopause diagnostics here
      DO IT=1,NTYPE
        AJ(J,J_CLRTOA,IT)=AJ(J,J_CLRTOA,IT)+OPNSKY*(SRNFLB(LM+LM_REQ+1)
     *     *COSZ2(I,J)-TRNFLB(LM+LM_REQ+1))*FTYPE(IT,I,J)
        AJ(J,J_CLRTRP,IT)=AJ(J,J_CLRTRP,IT)+OPNSKY*(SRNFLB(LTROPO(I,J))
     *     *COSZ2(I,J)-TRNFLB(LTROPO(I,J)))*FTYPE(IT,I,J)
        AJ(J,J_TOTTRP,IT)=AJ(J,J_TOTTRP,IT)+(SRNFLB(LTROPO(I,J))
     *     *COSZ2(I,J)-TRNFLB(LTROPO(I,J)))*FTYPE(IT,I,J)
      END DO
      AREG(JR,J_CLRTOA)=AREG(JR,J_CLRTOA)+OPNSKY*(SRNFLB(LM+LM_REQ+1)
     *     *COSZ2(I,J)-TRNFLB(LM+LM_REQ+1))*DXYP(J)
      AREG(JR,J_CLRTRP)=AREG(JR,J_CLRTRP)+OPNSKY*
     *     (SRNFLB(LTROPO(I,J))*COSZ2(I,J)-TRNFLB(LTROPO(I,J)))*DXYP(J)
      AREG(JR,J_TOTTRP)=AREG(JR,J_TOTTRP)+
     *     (SRNFLB(LTROPO(I,J))*COSZ2(I,J)-TRNFLB(LTROPO(I,J)))*DXYP(J)
C****
      IM1=I
      END DO
C****
C**** END OF MAIN LOOP FOR I INDEX
C****
  600 CONTINUE
C****
C**** END OF MAIN LOOP FOR J INDEX
C****
C**** Put this here so that print out values are correct
      if (ifirst.eq.1) then
        call writer(6,0)
        ifirst=0
      end if
C****
C**** ACCUMULATE THE RADIATION DIAGNOSTICS
C****
         DO 780 J=1,JM
         DXYPJ=DXYP(J)
         IMAX=IMAXJ(J)
         DO L=1,LM
           DO I=1,IMAX
             AJL(J,L,JL_SRHR)=AJL(J,L,JL_SRHR)+SRHR(L+1,I,J)*COSZ2(I,J)
             AJL(J,L,JL_TRCR)=AJL(J,L,JL_TRCR)+TRHR(L+1,I,J)
           END DO
         END DO
         DO 770 I=1,IMAX
         COSZ=COSZ2(I,J)
         JR=JREG(I,J)
         DO LR=1,LM_REQ
           ASJL(J,LR,3)=ASJL(J,LR,3)+SRHRS(LR,I,J)*COSZ
           ASJL(J,LR,4)=ASJL(J,LR,4)+TRHRS(LR,I,J)
         END DO
         DO KR=1,NDIUPT
           IF (I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
             DO INCH=1,NRAD
               IH=1+MOD(JHOUR+INCH-1,24) 
               ADIURN(IH,IDD_PALB,KR)=ADIURN(IH,IDD_PALB,KR)+
     *              (1.-SNFS(4,I,J)/S0)
               ADIURN(IH,IDD_GALB,KR)=ADIURN(IH,IDD_GALB,KR)+
     *              (1.-ALB(I,J,1))
               ADIURN(IH,IDD_ABSA,KR)=ADIURN(IH,IDD_ABSA,KR)+
     *              ((SNFS(4,I,J)-SNFS(1,I,J))*COSZ-TNFS(4,I,J)
     *              +TNFS(1,I,J))
             END DO
           END IF
         END DO

         DO IT=1,NTYPE
         AJ(J,J_SRINCP0,IT)=AJ(J,J_SRINCP0,IT)+(S0*COSZ)*FTYPE(IT,I,J)
         AJ(J,J_SRNFP0 ,IT)=AJ(J,J_SRNFP0 ,IT)+(SNFS(4,I,J)*COSZ)*
     *          FTYPE(IT,I,J)
         AJ(J,J_SRINCG ,IT)=AJ(J,J_SRINCG ,IT)+(SRHR(1,I,J)*COSZ/
     *          (ALB(I,J,1)+1.D-20))*FTYPE(IT,I,J)
         AJ(J,J_BRTEMP ,IT)=AJ(J,J_BRTEMP ,IT)+BTMPW(I,J) *FTYPE(IT,I,J)
         AJ(J,J_TRINCG ,IT)=AJ(J,J_TRINCG ,IT)+TRINCG(I,J)*FTYPE(IT,I,J)
         AJ(J,J_HSURF  ,IT)=AJ(J,J_HSURF  ,IT)-TNFS(4,I,J)*FTYPE(IT,I,J)
         AJ(J,J_SRNFP1 ,IT)=AJ(J,J_SRNFP1 ,IT)+SNFS(1,I,J)*COSZ
     *          *FTYPE(IT,I,J)
         AJ(J,J_HATM   ,IT)=AJ(J,J_HATM   ,IT)-TNFS(1,I,J)*FTYPE(IT,I,J)
         END DO
         AREG(JR,J_SRINCP0)=AREG(JR,J_SRINCP0)+(S0*COSZ)*DXYPJ
         AREG(JR,J_SRNFP0)=AREG(JR,J_SRNFP0)+(SNFS(4,I,J)*COSZ)*DXYPJ
         AREG(JR,J_SRNFP1)=AREG(JR,J_SRNFP1)+(SNFS(1,I,J)*COSZ)*DXYPJ
         AREG(JR,J_SRINCG)=AREG(JR,J_SRINCG)+
     *     (SRHR(1,I,J)*COSZ/(ALB(I,J,1)+1.D-20))*DXYPJ
C**** Note: confusing because the types for radiation are a subset
         AJ(J,J_SRNFG,ITOCEAN)=AJ(J,J_SRNFG,ITOCEAN)+(FSF(1,I,J)*COSZ)
     *        *FOCEAN(I,J)*(1.-RSI(I,J))
         AJ(J,J_SRNFG,ITLAKE) =AJ(J,J_SRNFG,ITLAKE) +(FSF(1,I,J)*COSZ)
     *        * FLAKE(I,J)*(1.-RSI(I,J))
         AJ(J,J_SRNFG,ITEARTH)=AJ(J,J_SRNFG,ITEARTH)+(FSF(4,I,J)*COSZ)
     *        *FEARTH(I,J)
         AJ(J,J_SRNFG,ITLANDI)=AJ(J,J_SRNFG,ITLANDI)+(FSF(3,I,J)*COSZ)
     *        * FLICE(I,J)
         AJ(J,J_SRNFG,ITOICE )=AJ(J,J_SRNFG,ITOICE )+(FSF(2,I,J)*COSZ)
     *        *FOCEAN(I,J)*RSI(I,J)
         AJ(J,J_SRNFG,ITLKICE)=AJ(J,J_SRNFG,ITLKICE)+(FSF(2,I,J)*COSZ)
     *        * FLAKE(I,J)*RSI(I,J)
C****
         AREG(JR,J_HATM)  =AREG(JR,J_HATM)  - TNFS(1,I,J)      *DXYPJ
         AREG(JR,J_SRNFG) =AREG(JR,J_SRNFG) +(SRHR(1,I,J)*COSZ)*DXYPJ
         AREG(JR,J_HSURF) =AREG(JR,J_HSURF) - TNFS(4,I,J)      *DXYPJ
         AREG(JR,J_BRTEMP)=AREG(JR,J_BRTEMP)+  BTMPW(I,J)      *DXYPJ
         AREG(JR,J_TRINCG)=AREG(JR,J_TRINCG)+ TRINCG(I,J)      *DXYPJ
         DO K=2,9
           JK=K+J_PLAVIS-2     ! accumulate 8 radiation diags.
           DO IT=1,NTYPE
             AJ(J,JK,IT)=AJ(J,JK,IT)+(S0*COSZ)*ALB(I,J,K)*FTYPE(IT,I,J)
           END DO
           AREG(JR,JK)=AREG(JR,JK)+(S0*COSZ)*ALB(I,J,K)*DXYPJ
         END DO
         AIJ(I,J,IJ_SRNFG)  =AIJ(I,J,IJ_SRNFG)  +(SRHR(1,I,J)*COSZ)
         AIJ(I,J,IJ_BTMPW)  =AIJ(I,J,IJ_BTMPW)  +BTMPW(I,J)
         AIJ(I,J,IJ_SRREF)  =AIJ(I,J,IJ_SRREF)  +S0*COSZ*ALB(I,J,2)
         AIJ(I,J,IJ_TRNFP0) =AIJ(I,J,IJ_TRNFP0) - TNFS(4,I,J)
         AIJ(I,J,IJ_SRNFP0) =AIJ(I,J,IJ_SRNFP0) +(SNFS(4,I,J)*COSZ)
         AIJ(I,J,IJ_SRINCG) =AIJ(I,J,IJ_SRINCG) +(SRHR(1,I,J)*COSZ/
     *        (ALB(I,J,1)+1.D-20))
         AIJ(I,J,IJ_SRINCP0)=AIJ(I,J,IJ_SRINCP0)+(S0*COSZ)
  770    CONTINUE
  780    CONTINUE
         DO L=1,LM
           DO I=1,IM
             DO J=J5S,J5N
               AIL(I,L,IL_REQ)=AIL(I,L,IL_REQ)+
     *              (SRHR(L+1,I,J)*COSZ2(I,J)+TRHR(L+1,I,J))*DXYP(J)
             END DO
             AIL(I,L,IL_R50N)=AIL(I,L,IL_R50N)+(SRHR(L+1,I,J50N)*COSZ2(I
     *            ,J50N)+TRHR(L+1,I,J50N))*DXYP(J50N)
             AIL(I,L,IL_R70N)=AIL(I,L,IL_R70N)+(SRHR(L+1,I,J70N)*COSZ2(I
     *            ,J70N)+TRHR(L+1,I,J70N))*DXYP(J70N)
           END DO
         END DO
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
      IH=1+JHOUR
      DO KR=1,NDIUPT
        ADIURN(IH,IDD_ISW,KR)=ADIURN(IH,IDD_ISW,KR)+
     *       S0*COSZ1(IJDD(1,KR),IJDD(2,KR))
      END DO

      RETURN
      END
