C**** DE001M12 E001M12 SOMTQ D_f90 DB285M12
C****
C**** Modified Model II diagnostics in double precision for work station
C**** subroutines in DB192SM15: All Diagnostics subroutines
C**** Additional diagnostics for moist energy fluxes
C**** changes for f90
C****
      MODULE DAGPCOM
!@sum  DAGCOMP Diagnostic model variables used in the printouts
!@auth Jean Lerner
!@var  1.0
      USE E001M12_COM, only : lm,sige,sig,psf,ptop
      USE RADNCB, only : LM_REQ
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION, DIMENSION(LM) :: PLE
      DOUBLE PRECISION, DIMENSION(LM+LM_REQ) :: PLM
      DOUBLE PRECISION PMTOP
      END MODULE DAGPCOM

C****                                                             IDACC
C**** CONTENTS OF AJ(J,N)  (SUM OVER LONGITUDE AND TIME OF)
C****   1  SRINCP0 (W/M**2)                                        2 RD
C****   2  SRNFP0 (W/M**2)                                         2 RD
C****   3  SRNFP1 (W/M**2)                                         2 RD
C****   4  SRABSATM=AJ(2)-AJ(6) (W/M**2)                           2 D1
C****   5  SRINCG (W/M**2)                                         2 RD
C****   6  SRNFG (W/M**2)                                          2 RD
C****   7  TRNFP0=AJ(74)+A2BYA1*AJ(9)/DTSRCE (W/M**2)              2 D1
C****   8  TRNFP1=AJ(75)+A2BYA1*AJ(9)/DTSRCE (W/M**2)              2 D1
C****   9  TRHDT (J/M**2)                                          1 SF
C****  10  RNFP0=AJ(2)+AJ(7) (W/M**2)                              2 D1
C****  11  RNFP1=AJ(3)+AJ(8) (W/M**2)                              2 D1
C****  12  RHDT=A1BYA2*AJ(6)*DTSRCE+AJ(9) (J/M**2)                 1 D1
C****  13  SHEATDT (J/M**2)                                        1 SF
C****  14  EVHDT (J/M**2)                                          1 SF
C****  15  F2DT (J/M**2)                                           1 GD
C****  16  HEATZ1=AJ(41)+AJ(42)                                    1 D1
C****  17  TG2 (K-273.16)                                          1 GD
C****  18  TG1 (K-273.16)                                          1 GD
C****  19  EVAP (KG/M**2)                                          1 GD
C****  20  PRCP=AJ(61)+AJ(62) (100 PA)                             1 D1
C****  21  TX (K-273.16)  (INTEGRAL OVER ATMOSPHERE OF)            4 DA
C****  22  TX1 (K-273.16)                                          4 DA
C****  23  TS (K-273.16)                                           3 SF
C****  24  DTH/DPHI  (STRATOSPHERE)                                4 DA
C****  25  DTH/DPHI  (TROPOSPHERE)                                 4 DA
C****  26  .0625*DTH*DLNP/(DU*DU+DV*DV)  (STRATOSPHERE)            4 DA
C****  27  .0625*DTH*DLNP/(DU*DU+DV*DV)  (TROPOSPHERE)             4 DA
C****  28  4*UMAX/(DX*SINJ)  (STRATOSPHERE)                        4 DA
C****  29  4*UMAX/(DX*SINJ)  (TROPOSPHERE)                         4 DA
C****  30  POICE (1)                                               1 GD
C****  31  PSNOW (1)                                               4 DA
C****  32  SW CORRECTION                                           2 RD
C****  33  OCEAN TRANSPORT                                         1 GD
C****  34  OCEAN TEMPERATURE AT MAX. MIXED LAYER DEPTH             1 GD
C****  35  T(J+1)-T(J-1)  (SUM OVER STRATOSPHERE OF)               4 DA
C****  36  T(J+1)-T(J-1)  (SUM OVER TROPOSPHERE OF)                4 DA
C****  37  SQRT(DTH/DLNP)/SINJ  (STRATOSPHERE)                     4 DA
C****  38  SQRT(DTH/DLNP)/SINJ  (TROPOSPHERE)                      4 DA
C****  39  ENERGP (J/M**2)                                         1 PR
C****  40  ERUN1 (J/M**2)                                          1 GP
C****  41  EDIFS (J/M**2)                                          1 GP
C****  42  F1DT (J/M**2)                                           1 GD
C****  43  ERUN2 (J/M**2)                                          1 GP
C****  44  HEATZ0=AJ(12)+AJ(13)+AJ(14)+AJ(39)-AJ(40) (J/M**2)      1 D1
C****  45  DIFS (KG/M**2)                                          1 GP
C****  46  AIFO ; BRUN2 ; CRUN2+CIFI                               1 GP
C****  47  RUN2 (KG/M**2)                                          1 GP
C****  48  DWTR2=AJ(45)-AJ(47) (KG/M**2)                           1 D1
C****  49  WTR1 (KG/M**2)                                          1 GD
C****  50  ACE1 (KG/M**2)                                          1 GD
C****  51  WTR2 (KG/M**2)                                          1 GD
C****  52  ACE2 (KG/M**2)                                          1 GD
C****  53  SNOW (KG/M**2)                                          1 GD
C****  54  RUN1 (KG/M**2)                                          1 GP
C****  55  BTEMPW-TF                                               2 RD
C****  56  HEATZ2=AJ(15)+AJ(43) (J/M**2)                           1 D1
C****  57  PCLDSS (1)  (COMPOSITE OVER ATMOSPHERE)                 2 RD
C****  58  PCLDMC (1)  (COMPOSITE OVER ATMOSPHERE)                 2 RD
C****  59  PCLD (1)  (COMPOSITE OVER ATMOSPHERE)                   2 RD
C****  60  CLDTOPMC=AJ(80)/AJ(58) (100 PA)                         0 D1
C****  61  PRCPSS (100 PA)                                         1 CN
C****  62  PRCPMC (100 PA)                                         1 CN
C****  63  Q*P (100 PA)  (INTEGRAL OVER ATMOSPHERE OF)             4 DA
C****  64  GAM  (K/M)  (*SIG(TROPOSPHERE)/GRAV)                    4 DA
C****  65  GAMM  (K-S**2/M**2)  (SIG(TROPOSPHERE)/GAMD)            4 DA
C****  66  GAMC  (K/M)                                             4 DA
C****  67  TRINCG (W/M**2)                                         2 RD
C****  68  ENERGY DIFFUSION INTO THERMOCLINE (W/M**2)           .5*9 MN
C****  69  FREE
C****  70  TRNFP0-TRNFG (W/M**2)                                   2 RD
C****  71  TRNFP1-TRNFG (W/M**2)                                   2 RD
C****  72  PLAVIS*S0*COSZ (W/M**2)                                 2 RD
C****  73  PLANIR*S0*COSZ (W/M**2)                                 2 RD
C****  74  ALBVIS*S0*COSZ (W/M**2)                                 2 RD
C****  75  ALBNIR*S0*COSZ (W/M**2)                                 2 RD
C****  76  SRRVIS*S0*COSZ (W/M**2)                                 2 RD
C****  77  SRRNIR*S0*COSZ (W/M**2)                                 2 RD
C****  78  SRAVIS*S0*COSZ (W/M**2)                                 2 RD
C****  79  SRANIR*S0*COSZ (W/M**2)                                 2 RD
C****  80  PBOTMC-PTOPMC (100 PA)                                  2 RD
C****
C**** CONTENTS OF APJ(J,N)  (SUM OVER LONGITUDE AND TIME OF)
C****   1  P (100 PA)                                              4 DA
C****   2  4*P4I (100 PA)  (UV GRID)                               4 DA
C****
C**** CONTENTS OF AJL(J,L,N)  (SUM OVER LONGITUDE AND TIME OF)
C****   1  FREE                                                    4 DA
C****   2  FREE                                                    4 DA
C****   3  FREE                                                    4 DA
C****   4  FREE                                                    4 DA
C****   5  FREE                                                    4 DA
C****   6  FREE                                                    4 DA
C****   7  FREE                                                    4 DA
C****   8  FMX(MC)*P (100 PA)                                      1 CN
C****   9  SRHR (W/M**2)                                           2 RD
C****  10  TRHR (W/M**2)                                           2 RD
C****  11  DTX(SS)*P (100 K*PA)                                    1 CN
C****  12  DT(DC)*P                                                1 CN
C****  13  DT(MC)*P (100 PA*K)  DRY HEATING                        1 CN
C****  14  FREE                                                    4 DA
C****  15  FREE                                                    4 DA
C****  16  (TH*SQRT(P)-THGM)**2/GMEAN(PR**(1-KAPA)*DTH/DPR)        4 DA
C****  17  FREE                                                    4 DA
C****  18  FREE                                                    4 DA
C****  19  PCLD*P (TOTAL)                                          1 CN
C****  20  DU/DT BY STRAT DRAG (M S-2)                             1 SD
C****  21  FREE                                                    4 DA
C****  22  FREE                                                    4 DA
C****  23  FREE                                                    4 DA
C****  24  FREE                                                    4 DA
C****  25  FREE                                                    4 DA
C****  26  FREE                                                    4 DA
C****  27  FREE                                                    4 DA
C****  28  PCLD*P (SS)                                             1 CN
C****  29  PCLD*P (MC)                                             1 CN
C****  30  FREE                                                    4 DA
C****  31  FREE                                                    4 DA
C****  32  FREE                                                    4 DA
C****  33  FREE                                                    4 DA
C****  34  FREE                                                    4 DA
C****  35  FREE                                                    4 DA
C****  36  (2F-2D(UDX))*16PV(TH-THMEAN)/(DTH/DSIG)+(SD-SDMEAN)*8U  4 DA
C****  37  P*V*((TH-THMEAN) * (DU/DP) / (DTH/DP) - U+UMEAN )       4 DA
C****  38  DU(DC)*P  (UV GRID)                                       GD
C****  39  DU(MC)*P (100 N/M/S)  (UV GRID)                         1 CN
C****  40  DU(ED)*P*(DTSURF*DSIG*ED/DZ**2)  (UV GRID)                SF
C****  41  U  (SUM OVER I, 135W to 110W)(PV GRID)                  4 DA
C****  42  V  (SUM OVER I, 135W to 110W)(PV GRID)                  4 DA
C****  43  SD  (SUM OVER I, 135W to 110W                           4 DA
C****  44  U  (SUM OVER I, 150E TO IM)   (PV GRID)                 4 DA
C****  45  V  (SUM OVER I, 150E TO IM)   (PV GRID)                 4 DA
C****  46  SD  (SUM OVER  I, 150E TO IM))                          4 DA
C****  47  V-V*  =D((V-VI)*(T-TI)/DTHDP)/DP                        4 DA
C****  48  4*PU4I*PV4I/P4I (100 N/S**2)  (UV GRID)                 4 DA
C****  49  4*PUV4I (100 N/S**2)  (UV GRID)                         4 DA
C****  50  DT(MC)*P (100 PA*K)  CHANGE OF PHASE                    1 CN
C****  51  CLHE*DQ(MC BEFORE COND)*P (100 PA*K)                    1 CN
C****  52  FREE                                                    4 DA
C****  53  FREE                                                    4 DA
C****  54  SIGMA  (VARIANCE FOR MOIST CONVECTION)                  1 CN
C****
C**** CONTENTS OF ASJL(J,L,N)  (SUM OVER LONGITUDE AND TIME OF)
C****   1  TX (C)                                                  4 DA
C****   2  PHI (M**2/S**2)                                         4 DA
C****   3  SRHR (W/M**2)                                           2 RD
C****   4  TRHR (W/M**2)                                           2 RD
C****
C**** CONTENTS OF AIJ(I,J,N)  (SUM OVER TIME OF)
C****   1  POICE (1)                                               1 GD
C****   2  PSNOW (1)                                               4 DA
C****   3  SNOW (KG/M**2)                                          4 DA
C****   4  SHDT (J/M**2)                                           1 SF
C****   5  PREC (KG/M**2)                                          1 PR
C****   6  EVAP (KG/M**2)                                          1 SF
C****   7  BETA (1)                                                1 GD
C****   8  SLP (100 PA-1000) (USING T1) (COMMENTED OUT)            4 DA
C****   8  4*P4 (100 PA)  (UV GRID)     (COMMENTED) (NO PRINTOUT)  4 DA
C****   8  PIJ (100 PA)                             (NO PRINTOUT)  4 DA
C****   9  PHI1000 (M**2/S**2)                                     4 DA
C****  10  PHI850 (M**2/S**2-1500*GRAV)                            4 DA
C****  11  PHI700-3000*GRAV                                        4 DA
C****  12  PHI500-5600*GRAV                                        4 DA
C****  13  PHI300-9500*GRAV                                        4 DA
C****  14  PHI100-16400*GRAV                                       4 DA
C****  15  PHI30-24000*GRAV                                        4 DA
C****  16  T850-273.16 (K-273.16)*GRAV)             (NO PRINTOUT)  4 DA
C****  17  PCLDMC (1)  (COMPOSITE OVER ATMOSPHERE)                 2 RD
C****  18  P-CLOUD TOP   (100 PA)                                  2 RD
C****  19  PCLD (1)  (COMPOSITE OVER ATMOSPHERE)                   2 RD
C****  20  16*P4*(SHA*T4+Z4)*V1*DSIG*DXV (100 W*M/S**2)  (UV GRID) 4 DA
C****  21  TRNFP0 (W/M**2)                                         2 RS
C****  22  SRHDT+TRHDT (J/M**2)                                    1 SF
C****  23  SRHDT+TRHDT+SHDT+EVHDT+ENRGP (J/M**2)                   1 SP
C****  24  SRNFP0 (W/M**2)                                         2 RD
C****  25  SRINCP0 (W/M**2)                                        2 RD
C****  26  SRNFG (W/M**2)                                          2 RD
C****  27  SRINCG (W/M**2)                                         2 RD
C****  28  TG1 (K-273.16)                                          1 GD
C****  29  POICE+PLICE+(IF SNOW)PEARTH                             4 DA
C****  30  DIURNAL DELTA TS (K) OVER SOIL       (NO PRINTOUT)   .5*9 MN
C****  31  DTHETA/DPHI (K S**2/M**2) IN TROPOSPHERE                4 DA
C****  32  RUN1 OVER EARTH  (KG/M**2)                              1 PG
C****  33  TS (K-273.16)  (USING LAPSE RATE FROM TX1)  (COMM'D OUT)4 DA
C****  33  RUN1 OVER LAND ICE  (KG/M**2)            (NO PRINTOUT)  1 PG
C****  34  SURFACE WIND SPEED (M/S)                                3 SF
C****  35  TS (K-273.16)                                           3 SF
C****  36  US (M/S)                                                3 SF
C****  37  VS (M/S)                                                3 SF
C****  38  PSL (100 PA-1000)  (USING TS)                           4 DA
C****  39  UJET (M/S)                                              4 DA
C****  40  VJET (M/S)                                              4 DA
C****  41  PCLD(LOW) (1)                                           2 RD
C****  42  PCLD(MID) (1)                                           2 RD
C****  43  PCLD(HIGH) (1)                                          2 RD
C****  44  BTEMPW-TF (K-273.16)                                    2 RD
C****  45  PLAVIS*S0*COSZ (W/M**2)                                 2 RD
C****  46  TGO2=ODATA(4)  (C)                                   .5*9 MN
C****  47  TAUS  (MOMENTUM SURFACE DRAG) (kg/m**2)  (NO PRINTOUT)  3 SF
C****  48  TAUUS (MOMENTUM SURFACE DRAG) (kg/m**2)  (NO PRINTOUT)  3 SF
C****  49  TAUVS (MOMENTUM SURFACE DRAG) (kg/m**2)  (NO PRINTOUT)  3 SF
C****  50  WATER1+WATER2+ICE1+ICE2  (FOR EARTH POINTS ONLY)        1 GD
C****  51  QS                                       (NO PRINTOUT)  3 SF
C****  52  MAX(0,33-1.8*DAILY MEAN ON TS IN C)                  .5*9 MN
C****  53  40.6+.72*(2TS(C)-(QSATS-QS)*LHA/SHA)                    3 SF
C****  54  18*(DEL(TG)/DEL(TS)-1), DEL=DIURNAL MAX-MIN          .5*9 MN
C****  55  8*P*U*Q (VERTICALLY INTEGRATED)  (12.5 PA*M/S)          4 DA
C****  56  8*P*V*Q (VERTICALLY INTEGRATED)  (12.5 PA*M/S)          4 DA
C****  57  TGO=ODATA(1)  (C)                                       1 GD
C****  58  ACE2OI=ODATA(3)*POICE  (KG/M**2)                        1 GD
C****  59  WIND SPEED IN TOP LAYER (M/S)                           1 SD
C****  60  TGO12=ODATA(5)  (C)                                  .5*9 MN
C****  61  EVAP*POCEAN  (KG/M**2)                                  1 GD
C****  62  EVAP*POICE  (KG/M**2)                                   1 GD
C****  63  EVAP OVER LAND ICE  (KG/M**2)                           1 GD
C****  64  EVAP OVER EARTH  (KG/M**2)                              1 GD
C****  65  F0DT*POCEAN, NET HEAT AT Z0  (J/M**2)                   1 GD
C****  66  F0DT*POICE, NET HEAT AT Z0  (J/M**2)                    1 GD
C****  67  F0DT, NET HEAT AT Z0 OVER LAND ICE  (J/M**2)            1 GD
C****  68  F0DT, NET HEAT AT Z0 OVER EARTH  (J/M**2)               1 GD
C****  69  F1DT OVER LAND ICE  (J/M**2)                            1 PG
C****  70  SNOW FALL  (KG/M**2)                                    1 PR
C****  71  SURF AIR TEMP OVER LAND ICE  (C)                  NSURF*1 SF
C****  72  F2DT OVER LAND ICE  (J/M**2)                            1 PG
C****  73  SHDT OVER LAND ICE  (J/M**2)                            3 SF
C****  74  EVHDT OVER LAND ICE  (J/M**2)                           3 SF
C****  75  TRHDT OVER LAND ICE  (J/M**2)                           3 SF
C****  76  MAX(COMPOSITE TS)                                      12 SF
C****  77  MIN(COMPOSITE TS)                                      12 SF
C****  78  MIN(DIURNAL MAX OF COMPOSITE TS)                       12 MN
C****  79  POTENTIAL EVAPORATION (KG/M**2)                         1 EA
C****  80  MAX TS OVER EARTH FOR CURRENT DAY (K)                .5*9 MN
C****  81  LIQUID WATER PATH (kg/M**2)                             1 CL
C****  82  SHALLOW CONVECTIVE CLOUD COVER  (1)                     1 CL
C****  83  DEEP CONVECTIVE CLOUD COVER     (1)                     1 CL
C****  84  DEEP CONVECTIVE CLOUD FREQUENCY (1)                     1 CL
C****  85  SHALLOW CONVECTIVE CLOUD FREQUENCY (1)                  1 CL
C****  86  INCIDENT MTN EAST MOMENTUM FLUX   (MB-M/S**2)           1 SD
C****  87  INCIDENT MTN SOUTH MOMENTUM FLUX  (MB-M/S**2)           1 SD
C****  88  EAST-WEST POTENTIAL ENTHALPY FLUX (W)            /3600.*1 DY
C****  89  NORTH-SOUTH POTENTIAL ENTHALPY FLUX (W)          /3600.*1 DY
C****  90  EAST-WEST MASS FLUX (KG/S)              100./GRAV/3600.*1 DY
C****  91  NORTH-SOUTH MASS FLUX (KG/S)            100./GRAV/3600.*1 DY
C****  92  EAST-WEST WATER VAPOR FLUX (KG/S)                /3600.*1 DY
C****  93  NORTH-SOUTH WATER VAPOR FLUX (KG/S)              /3600.*1 DY
C****  94  EAST-WEST GEOPOTENTIAL FLUX (W)                  /3600.*1 DY
C****  95  NORTH-SOUTH GEOPOTENTIAL FLUX (W)                /3600.*1 DY
C****  96  Energy Outflow by Rivers (10**10 W)     E-10/DTS*IDACC(1) RV
C****  97  Mass Outflow by Rivers (10**5 kg/s)      E-5/DTS*IDACC(1) RV
C****  98  FREE
C****  99  FREE
C**** 100  FREE
C****
C**** CONTENTS OF AIL(I,L,N)  (SUM OVER TIME OF)
C**** WE ARE NOT TAKING INTO ACCOUNT THE VARIATION OF MASS
C****   1  U (M/S) (SUM FOR J=JEQ+1,JEQ,JEQ-1,JEQ-2)  (PU GRID)    4 DA
C****   2  V (M/S) (SUM FOR J=JEQ+1,JEQ,JEQ-1,JEQ-2)  (PU GRID)    4 DA
C****   3  SD (100 N/S) (SUM FOR J=JEQ,JEQ-1,JEQ-2)                4 DA
C****   4  TX (K-273.16) (SUM FOR J=JEQ,JEQ-1,JEQ-2)               4 DA
C****   5  RH (1) (SUM FOR J=JEQ,JEQ-1,JEQ-2)                      4 DA
C****   6  DTX(MC)*P*DA (100 K*N) (SUM FOR J=JEQ,JEQ-1,JEQ-2)      1 CN
C****   7  (SRHR+TRHR)*DA (W) (SUM FOR J=JEQ,JEQ-1,JEQ-2)          2 RD
C****   9  SD (100 N/S) (AT LAT 50 N) (COMMENTED OUT)              4 DA
C****  10  TX-273.16  (AT LAT 50 N)                                4 DA
C****  11  SR+TR  (AT LAT 50 N)                                    2 RD
C****  12  2*U  (AT LAT 50 N)                                      4 DA
C****  13  SD  (AT LAT 70 N)         (COMMENTED OUT)               4 DA
C****  14  TX-273.16  (AT LAT 70 N)  (COMMENTED OUT)               4 DA
C****  15  SR+TR  (AT LAT 70 N)                                    2 RD
C****  16  2*U  (AT LAT 70 N)        (COMMENTED OUT)               4 DA
C****
C**** CONTENTS OF AIJG(I,J,K)  (SUM OVER TIME OF)
C****   1  LAYER 1 REL SAT OF BARE AMD VEG. SOIL (%)            .5*9 MN
C****   2  LAYER 2 REL SATURATION OF BARE SOIL (%)              .5*9 MN
C****   3  LAYER 3 REL SATURATION OF BARE SOIL (%)              .5*9 MN
C****   4  LAYER 4 REL SATURATION OF BARE SOIL (%)              .5*9 MN
C****   5  BETA (BARE SOIL) (%)                                    1 EA
C****   6  BETA (PENMAN) (%)                                       1 EA
C****   7  CANOPY  REL SATURATION (%)                           .5*9 MN
C****   8  LAYER 1 REL SATURATION OF VEG. SOIL (%)              .5*9 MN
C****   9  LAYER 2 REL SATURATION OF VEG. SOIL (%)               5*9 MN
C****  10  LAYER 3 REL SATURATION OF VEG. SOIL (%)               5*9 MN
C****  11  BETA (BARE SOIL & VEGETATION) (%)                       1 EA
C****  12  CONDUCTANCE OF ATMOSPHERE (.01M/S)                      1 EA
C****  13  CONDUCTANCE OF CANOPY (.01M/S)                          1 EA
C****  14  PENMAN POTENTIAL EVAPORATION (KG/M**2)                  1 EA
C****  15  TEMP OF LAYER 1 BARE SOIL AND SNOW (C)                  1 EA
C****  16  TEMP OF SOIL LAYER 2 - BARE SOIL (C)                    1 EA
C****  17  TEMP OF SOIL LAYER 3 - BARE SOIL (C)                    1 EA
C****  18  BARE SOIL EVAPORATION (KG/M**2)                         1 EA
C****  19  DRY CANOPY EVAPORATION (KG/M**2)                        1 EA
C****  20  WET CANOPY EVAPORATION (KG/M**2)                        1 EA
C****  21  TEMP OF CANOPY AND SNOW (C)                             1 EA
C****  23  TEMP OF SOIL LAYER 2 1 VEGETATED SOIL (C)               1 EA
C****  23  TEMP OF SOIL LAYER 2 - VEGETATED SOIL (C)               1 EA
C****  24  TEMP OF SOIL LAYER 3 - VEGETATED SOIL (C)               1 EA
C****  25  AVERAGE WATER TABLE (M)                                 1 EA
C****  26  BETAV, OVER VEGETATION (PERCENT)                        1 EA
C****  27  BETAT, TRANSPIRATION (PERCENT)                          1 EA
C****
C**** CONTENTS OF IDACC(N), NUMBER OF ACCUMULATION TIMES OF
C****   1  SOURCE TERMS  (DETERMINED BY NDYN)
C****   2  RADIATION SOURCE TERMS  (DETERMINED BY NRAD)
C****   3  SURFACE INTERACTION SOURCE TERMS  (DETERMINED BY NDASF)
C****   4  QUANTITIES IN DIAGA  (DETERMINED BY NDAA)
C****   5  ENERGY NUMBERS IN DIAG4  (DEYERMINED BY NDA4)
C****   6  KINETIC ENERGY IN DIAG5 FROM DYNAMICS  (DETERMINED BY NDA5K)
C****   7  ENERGY IN DIAG5 FROM DYNAMICS  (DETERMINED BY NDA5D)
C****   8  ENERGY IN DIAG5 FROM SOURCES  (DETERMINED BY NDA5S)
C****   9  WAVE ENERGY IN DIAG7  (EVERY 12 HOURS)
C****  10  ENERGY IN DIAG5 FROM FILTER  (DETERMINED BY NFILTR)
C****  11  USED FOR T-DIAGNOSTICS CHECKT:  0=OFF,1=ON
C****  12  ALWAYS =1 (UNLESS SEVERAL RESTART FILES WERE ACCUMULATED)
C****
C**** CONTENTS OF AUXILIARY ARRAYS (TSFREZ(I,J,1-2),TDIURN(I,J,N))
C****   1  FIRST DAY OF GROWING SEASON (JULIAN DAY)
C****   2  LAST DAY OF GROWING SEASON (JULIAN DAY)
C****
C****   1  MIN TG1 OVER EARTH FOR CURRENT DAY (C)
C****   2  MAX TG1 OVER EARTH FOR CURRENT DAY (C)
C****   3  MIN TS OVER EARTH FOR CURRENT DAY (K)
C****   4  MAX TS OVER EARTH FOR CURRENT DAY (K)
C****   5  SUM OF COMPOSITE TS OVER TIME FOR CURRENT DAY (C)
C****   6  MAX COMPOSITE TS FOR CURRENT DAY (K)
C****   7  MAX TG1 OVER OCEAN ICE FOR CURRENT DAY (C)
C****   8  MAX TG1 OVER LAND ICE FOR CURRENT DAY (C)
C****
      SUBROUTINE DIAGA (U,V,T,P,Q,PIT,SD)
!@sum  DIAGA accumulate various diagnostics during dynamics
!@auth Original Development Team
!@ve   1.0
      USE CONSTANT, only : grav,rgas,kapa,lhe,sha
      USE E001M12_COM, only :
     &     im,imh,fim,byim,jm,jeq,lm,ls1,idacc,psf,ptop,psfmpt,
     &     mdyn,mdiag,
     &     sig,sige,dsig,dsigo,
     &     zatmo,fland,flice,fearth, gdata,WM
      USE DYNAMICS, only : pk
      USE GEOM, only :
     &     AREAG,COSP,DLAT,DXV,DXYN,DXYP,DXYS,
     &     DXYV,DYP,FCOR,IMAXJ,RAVPN,RAVPS,SINP

      USE RADNCB, only : rqt,lm_req
      USE PBLCOM, only : tsavg
      USE DAGPCOM
      USE DAGCOM, only : aj,bj,cj,areg,jreg,apj,ajl,asjl,ail,
     &     aij,ij_dtdp,ij_pev,ij_phi1k,ij_pres,ij_puq,ij_pvq,
     &     ij_rsit,ij_rsnw,ij_slp,ij_snow,ij_t850,ij_ujet,ij_vjet,j_tx1,
     *     j_rsnow,j_tx,j_qp,j_dtdjt,j_dtdjs,j_dtdgtr,j_dtsgst,j_rictr,
     *     j_rostr,j_ltro,j_ricst,j_rosst,j_lstr,j_gamm,j_gam,j_gamc

      USE OCEAN, only : odata

      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: U,V,T,Q
      DOUBLE PRECISION, DIMENSION(IM,JM) :: P
      LOGICAL POLE
      DOUBLE PRECISION, DIMENSION(IM,JM) :: PIT
      DOUBLE PRECISION, DIMENSION(IM,JM,LM-1) :: SD
      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: W
      COMMON/WORK2/W
      DOUBLE PRECISION, DIMENSION(LM) :: GMEAN
      DOUBLE PRECISION, DIMENSION(JM) :: TIL,UI,UMAX,PI,EL,RI,
     &     SOCEAN,SLAND,SOICE,DUDVSQ
      DOUBLE PRECISION, DIMENSION(JM,LM) ::
     &     THJL,THSQJL,SPI,PHIPI,TPI
      DOUBLE PRECISION, DIMENSION(JM,LM-1) :: SDMEAN
      DOUBLE PRECISION, DIMENSION(IM,JM) :: PUV
      DOUBLE PRECISION, DIMENSION(LM_REQ) :: TRI
      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: PHI,TX
      DOUBLE PRECISION, DIMENSION(IM) :: THSEC,PSEC,SQRTP,PDA
      COMMON/WORK3/PHI,TX
      INTEGER, DIMENSION(LM) :: LUPA,LDNA
      DOUBLE PRECISION, DIMENSION(LM) :: D2SIG,PKS
      CHARACTER*16 TITLE
      DOUBLE PRECISION, DIMENSION(7), PARAMETER ::
     &     PMB=(/1000.,850.,700.,500.,300.,100.,30./),
     &     GHT=(/0.,1500.,3000.,5600.,9500.,16400.,24000./)
      INTEGER :: IFIRST
      DATA IFIRST/1/
      DOUBLE PRECISION, PARAMETER :: ONE=1.,ZERO20=1.E-20,P1000=1000.
      INTEGER :: I,IM1,J,K,L,J50N,J70N,JET,JR,KM,LDN,LUP,
     &     IP1,IMAX,LM1,LP1,LR,KBEGIN,KEND,
     &     I150E,I110W,I135W,J5NUV,J5SUV,J5N,J5S
      INTEGER MCLOCK ! external
      DOUBLE PRECISION THBAR ! external
      DOUBLE PRECISION ::
     &     ADTDL,AQ,ASS,AT1,ATX,BBYG,BBYGV,BDTDL,BETA,
     &     BQ,BSCOV,BSS,BT1,BTX,BYDXYV,BYSDSG,CDTDL,
     &     CQ,CSCOV,CSS,CT1,CTX,DLNP,DLNP12,DLNP23,
     &     DLNS,DP,DS,DT2,DTHDP,DU,DUDP,DUDX,DV,
     &     DXYPJ,ELX,EPSLON,ESEPS,FPHI,GAMC,GAMD,
     &     GAMM,GAMX,GBYRB,GMEANL,P1,P4,P4I,
     &     PDN,PE,PEARTH,PEQ,PEQM1,PEQM2,PHIRI,
     &     PIBYIM,PIJ,PITIJ,PITMN,PKE,PL,PLAND,PLICE,
     &     POCEAN,POICE,PRQ1,PRT,PSMPT4,PSNOW,PU4I,PUP,
     &     PUV4I,PV4I,PVTHP,
     &     QLH,RBBYG,ROSSX,SCOVL,SCOVOI,
     &     SDMN,SDPU,SMALL,SP,SP2,
     &     SS,T4,THETA,THGM,THMN,
     &     TPIL,TZL,UAMAX,UMN,
     &     UPE,VPE,X,Z4

C**** QSAT=(RAIR/RVAPOR)*6.1071*EXP((L/RVAPOR)*(1/TF-1/T))/P
      DOUBLE PRECISION, PARAMETER ::
     &     AQSAT=3.797915,BQSAT=7.93252E-6,CQSAT=2.166847E-3
      DOUBLE PRECISION :: QSAT,TM,PR,QL
      QSAT(TM,PR,QL)=AQSAT*EXP(QL*(BQSAT-CQSAT/TM))/PR
      KBEGIN = MCLOCK ()
      IDACC(4)=IDACC(4)+1
      IF (IFIRST.NE.1) GO TO 50
      IFIRST=0
C**** INITIALIZE CERTAIN QUANTITIES
      L=LM+1
    3 L=L-1
      IF (L.EQ.1) GO TO 4
      IF (.25*(SIGE(L-1)+2*SIGE(L)+SIGE(L+1))*PSFMPT+PTOP.LT.250.)
     *   GO TO 3
    4 JET=L
      WRITE (6,888) JET
  888 FORMAT (' JET WIND LEVEL FOR DIAG',I3)
c      BYIM=1./FIM
c      SHA=RGAS/KAPA
      BYSDSG=1./(1.-SIGE(LM+1))
      BETA=.0065
      BBYG=BETA/GRAV
      RBBYG=RGAS*BETA/GRAV
      GBYRB=GRAV/(RGAS*BETA)
      EPSLON=1.
      KM=0
      DO 5 K=1,7
      IF (PMTOP.GT.PMB(K)) GO TO 6
    5 KM=KM+1
    6 CONTINUE    ! JEQ=2.+.5*(JM-1.)
      J50N=(50.+90.)*(JM-1.)/180.+1.5
      J70N=(70.+90.)*(JM-1.)/180.+1.5
      I150E = IM*(180+150)/360+1   ! WEST EDGE OF 150 EAST
      I110W = IM*(180-110)/360+1   ! WEST EDGE OF 110 WEST
      I135W = IM*(180-135)/360+1   ! WEST EDGE OF 135 WEST
      J5NUV = (90.+5.)*(JM-1.)/180.+2.
      J5SUV = (90.-5.)*(JM-1.)/180.+2.
      J5N   = (90.+5.)*(JM-1.)/180.+1.5
      J5S   = (90.-5.)*(JM-1.)/180.+1.5
      PRQ1=.75*PMTOP
      DLNP12=LOG(.75/.35)
      DLNP23=LOG(.35/.1)
      DO 10 L=1,LM
      LUPA(L)=L+1
   10 LDNA(L)=L-1
      LDNA(1)=1
      LUPA(LM)=LM
      DO 20 L=1,LM
   20 D2SIG(L)=SIG(LUPA(L))-SIG(LDNA(L))
      DO 30 L=LS1,LM
   30 PKS(L)=(SIG(L)*PSFMPT+PTOP)**KAPA
      PSMPT4=4.*PSFMPT
   50 CONTINUE
C****
C**** FILL IN HUMIDITY AND SIGMA DOT ARRAYS AT THE POLES
C****
      DO 65 L=1,LM
      DO 65 I=2,IM
      Q(I,1,L)=Q(1,1,L)
   65 Q(I,JM,L)=Q(1,JM,L)
C****
C**** CALCULATE PK AND TX, THE REAL TEMPERATURE
C****
      DO 80 L=1,LS1-1
      PK(L,1,1)=(SIG(L)*P(1,1)+PTOP)**KAPA
      TX(1,1,L)=T(1,1,L)*PK(L,1,1)
      PK(L,1,JM)=(SIG(L)*P(1,JM)+PTOP)**KAPA
      TX(1,JM,L)=T(1,JM,L)*PK(L,1,JM)
      DO 70 I=2,IM
      T(I,1,L)=T(1,1,L)
      T(I,JM,L)=T(1,JM,L)
      PK(L,I,1)=PK(L,1,1)
      TX(I,1,L)=TX(1,1,L)
      PK(L,I,JM)=PK(L,1,JM)
   70 TX(I,JM,L)=TX(1,JM,L)
      DO 80 J=2,JM-1
      DO 80 I=1,IM
      PK(L,I,J)=(SIG(L)*P(I,J)+PTOP)**KAPA
   80 TX(I,J,L)=T(I,J,L)*PK(L,I,J)
      DO 83 L=LS1,LM
      DO 82 I=2,IM
      T(I,1,L)=T(1,1,L)
   82 T(I,JM,L)=T(1,JM,L)
      DO 83 J=1,JM
      DO 83 I=1,IM
      PK(L,I,J)=PKS(L)
   83 TX(I,J,L)=T(I,J,L)*PK(L,I,J)  !      PKS(L)
C****
C**** CALCULATE PUV, THE MASS WEIGHTED PRESSURE
C****
      DO 90 J=2,JM
      I=IM
      DO 85 IP1=1,IM
      PUV(I,J)=RAVPN(J-1)*(P(I,J-1)+P(IP1,J-1))+
     *         RAVPS(  J)*(P(I,  J)+P(IP1,  J))
   85 I=IP1
   90 CONTINUE
C****
C**** J LOOPS FOR ALL PRIMARY GRID ROWS
C****
      DO 190 J=1,JM
      POLE=.FALSE.
      IF (J.EQ.1.OR.J.EQ.JM) POLE=.TRUE.
      IMAX=IMAXJ(J)
      DXYPJ=DXYP(J)
C**** NUMBERS ACCUMULATED FOR A SINGLE LEVEL
      AT1=0.
      BT1=0.
      CT1=0.
      BSCOV=0.
      CSCOV=0.
      PI(J)=0.
      SLAND(J)=0.
      SOICE(J)=0.
      SOCEAN(J)=0.
      DO 120 I=1,IMAX
      JR=JREG(I,J)
      PLAND=FLAND(I,J)
      POICE=ODATA(I,J,2)*(1.-PLAND)
      PLICE=FLICE(I,J)
      POCEAN=(1.-PLAND)-POICE
      PEARTH=FEARTH(I,J)
      SLAND(J)=SLAND(J)+PLAND
      SOICE(J)=SOICE(J)+POICE
      SOCEAN(J)=SOCEAN(J)+POCEAN
      AT1=AT1+(TX(I,J,1)-273.16)*POCEAN
      BT1=BT1+(TX(I,J,1)-273.16)*PLAND
      CT1=CT1+(TX(I,J,1)-273.16)*POICE
      AREG(JR,J_TX1)=AREG(JR,J_TX1)+(TX(I,J,1)-273.16)*DXYPJ
      SCOVL=0.
      IF (GDATA(I,J,2).GT.0.) SCOVL=PEARTH
      IF (GDATA(I,J,12).GT.0.) SCOVL=SCOVL+PLICE
      BSCOV=BSCOV+SCOVL
      SCOVOI=0.
      IF (GDATA(I,J,1).GT.0.) SCOVOI=POICE
      CSCOV=CSCOV+SCOVOI
      AREG(JR,J_RSNOW)=AREG(JR,J_RSNOW)+(SCOVL+SCOVOI)*DXYPJ
      PI(J)=PI(J)+P(I,J)
      AIJ(I,J,IJ_RSNW)=AIJ(I,J,IJ_RSNW)+(SCOVOI+SCOVL)
      AIJ(I,J,IJ_SNOW)=AIJ(I,J,IJ_SNOW)+(GDATA(I,J,1)*POICE+GDATA(I,J,2)
     *     *PEARTH+GDATA(I,J,12)*PLICE)
C     TS=TX(I,J,1)*((P(I,J)+PTOP)/(SIG(1)*P(I,J)+PTOP))**RBBYG
C     AIJ(I,J,IJ_SLP1)=AIJ(I,J,IJ_SLP1)+((P(I,J)+PTOP)*(1.+BBYG*ZATMO(I
C    *  ,J)/TS)**GBYRB-P1000)
      AIJ(I,J,IJ_PRES)=AIJ(I,J,IJ_PRES)+ P(I,J)
      PSNOW=0.
      IF (GDATA(I,J,2).GT.0.) PSNOW=PEARTH
      AIJ(I,J,IJ_RSIT)=AIJ(I,J,IJ_RSIT)+POICE+PLICE+PSNOW
C     AIJ(I,J,IJ_TS1)=AIJ(I,J,IJ_TS1)+(TS-273.16)
      AIJ(I,J,IJ_SLP)=AIJ(I,J,IJ_SLP)+((P(I,J)+PTOP)*(1.+BBYG*ZATMO(I,J)
     *     /TSAVG(I,J))**GBYRB-P1000)
  120 CONTINUE
      AJ(J,J_TX1)=AJ(J,J_TX1)+AT1
      BJ(J,J_TX1)=BJ(J,J_TX1)+BT1
      CJ(J,J_TX1)=CJ(J,J_TX1)+CT1
      BJ(J,J_RSNOW)=BJ(J,J_RSNOW)+BSCOV
      CJ(J,J_RSNOW)=CJ(J,J_RSNOW)+CSCOV
      APJ(J,1)=APJ(J,1)+PI(J)
C**** GEOPOTENTIALS CALCULATED FOR EACH LAYER
      DO 160 I=1,IMAX
      PIJ=P(I,J)
      P1=SIG(1)*PIJ+PTOP
      PUP=SIG(2)*PIJ+PTOP
      IF (ABS(TX(I,J,2)-TX(I,J,1)).LT.EPSLON) GO TO 152
      BBYGV=LOG(TX(I,J,1)/TX(I,J,2))/(RGAS*LOG(P1/PUP))
      PHI(I,J,1)=ZATMO(I,J)+TX(I,J,1)
     *  *(((PIJ+PTOP)/P1)**(RGAS*BBYGV)-1.)/BBYGV
      PHI(I,J,2)=PHI(I,J,1)+(TX(I,J,1)-TX(I,J,2))/BBYGV
      GO TO 154
  152 PHI(I,J,1)=ZATMO(I,J)+RGAS*TX(I,J,1)*LOG((PIJ+PTOP)/P1)
      PHI(I,J,2)=PHI(I,J,1)+RGAS*.5*(TX(I,J,1)+TX(I,J,2))*LOG(P1/PUP)
  154 DO 160 L=3,LM
      IF(L.GE.LS1) PIJ=PSFMPT
      PDN=PUP
      PUP=SIG(L)*PIJ+PTOP
      IF (ABS(TX(I,J,L)-TX(I,J,L-1)).LT.EPSLON) GO TO 156
      BBYGV=LOG(TX(I,J,L-1)/TX(I,J,L))/(RGAS*LOG(PDN/PUP))
      PHI(I,J,L)=PHI(I,J,L-1)+(TX(I,J,L-1)-TX(I,J,L))/BBYGV
      GO TO 160
  156 PHI(I,J,L)=PHI(I,J,L-1)+RGAS*.5*(TX(I,J,L-1)+TX(I,J,L))
     *  *LOG(PDN/PUP)
  160 CONTINUE
      IF (.NOT.POLE) GO TO 170
      DO 162 L=1,LM
      DO 162 I=2,IM
  162 PHI(I,J,L)=PHI(1,J,L)
C**** CALCULATE GEOPOTENTIAL HEIGHTS AT SPECIFIC MILLIBAR LEVELS
  170 DO 180 I=1,IMAX
      PIJ=P(I,J)
      PL=SIG(1)*PIJ+PTOP
      K=1
      L=1
  172 L=L+1
      IF(L.GE.LS1) PIJ=PSFMPT
      PDN=PL
      PL=SIG(L)*PIJ+PTOP
      IF (PMB(K).LT.PL.AND.L.LT.LM) GO TO 172
      IF (ABS(TX(I,J,L)-TX(I,J,L-1)).LT.EPSLON) GO TO 176
      BBYGV=(TX(I,J,L-1)-TX(I,J,L))/(PHI(I,J,L)-PHI(I,J,L-1))
  174 AIJ(I,J,IJ_PHI1K-1+K)=AIJ(I,J,IJ_PHI1K-1+K)+(PHI(I,J,L)
     *     -TX(I,J,L)*((PMB(K)/PL)**(RGAS*BBYGV)-1.)/BBYGV-GHT(K)*GRAV)
      IF (K.EQ.2) AIJ(I,J,IJ_T850)=AIJ(I,J,IJ_T850)+(TX(I,J,L)-273.16
     *     +(TX(I,J,L-1)-TX(I,J,L))*LOG(PMB(K)/PL)/LOG(PDN/PL))
      IF (K.GE.KM) GO TO 180
      K=K+1
      IF (PMB(K).LT.PL.AND.L.LT.LM) GO TO 172
      GO TO 174
 176  AIJ(I,J,IJ_PHI1K-1+K)=AIJ(I,J,IJ_PHI1K-1+K)+(PHI(I,J,L)
     *     -RGAS*TX(I,J,L)*LOG(PMB(K)/PL)-GHT(K)*GRAV)
      IF (K.EQ.2) AIJ(I,J,IJ_T850)=AIJ(I,J,IJ_T850)+(TX(I,J,L)-273.16)
      IF (K.GE.KM) GO TO 180
      K=K+1
      IF (PMB(K).LT.PL.AND.L.LT.LM) GO TO 172
      GO TO 176
  180 CONTINUE
  190 CONTINUE
C**** ACCUMULATION OF TEMP., POTENTIAL TEMP., Q, AND RH
      DO 250 J=1,JM
      IMAX=IMAXJ(J)
      DXYPJ=DXYP(J)
      DO 230 L=1,LM
      ATX=0.
      BTX=0.
      CTX=0.
      TPI(J,L)=0.
      AQ=0.
      BQ=0.
      CQ=0.
      PHIPI(J,L)=0.
C     QPI=0.
      SPI(J,L)=0.
C     RHPI=0.
      DO 220 I=1,IMAX
      JR=JREG(I,J)
      PLAND=FLAND(I,J)
      POICE=ODATA(I,J,2)*(1.-PLAND)
      POCEAN=(1.-PLAND)-POICE
      PIJ=P(I,J)
      IF(L.GE.LS1) PIJ=PSFMPT
      ATX=ATX+(TX(I,J,L)-273.16)*POCEAN
      BTX=BTX+(TX(I,J,L)-273.16)*PLAND
      CTX=CTX+(TX(I,J,L)-273.16)*POICE
      AQ=AQ+(Q(I,J,L)+WM(I,J,L))*PIJ*POCEAN
      BQ=BQ+(Q(I,J,L)+WM(I,J,L))*PIJ*PLAND
      CQ=CQ+(Q(I,J,L)+WM(I,J,L))*PIJ*POICE
      AREG(JR,J_QP)=AREG(JR,J_QP)+(Q(I,J,L)+WM(I,J,L))*PIJ*DSIG(L)*DXYPJ
      AREG(JR,J_TX)=AREG(JR,J_TX)+(TX(I,J,L)-273.16)*DSIG(L)*BYSDSG
     *     *DXYPJ
      TPI(J,L)=TPI(J,L)+(TX(I,J,L)-273.16)*PIJ
      PHIPI(J,L)=PHIPI(J,L)+PHI(I,J,L)*PIJ
C     QPI=QPI+Q(I,J,L)*PIJ
      SPI(J,L)=SPI(J,L)+T(I,J,L)*PIJ
C     QLH=LHE
C     QSATL=QSAT(TX(I,J,L),SIG(L)*PIJ+PTOP,QLH)
C     IF (QSATL.GT.1.) QSATL=1.
C     RHPI=RHPI+Q(I,J,L)*PIJ/QSATL
  220 CONTINUE
      AJ(J,J_TX)=AJ(J,J_TX)+ATX*DSIG(L)*BYSDSG
      BJ(J,J_TX)=BJ(J,J_TX)+BTX*DSIG(L)*BYSDSG
      CJ(J,J_TX)=CJ(J,J_TX)+CTX*DSIG(L)*BYSDSG
      AJ(J,J_QP)=AJ(J,J_QP)+AQ*DSIG(L)
      BJ(J,J_QP)=BJ(J,J_QP)+BQ*DSIG(L)
      CJ(J,J_QP)=CJ(J,J_QP)+CQ*DSIG(L)
C     AJL(J,L,1)=AJL(J,L,1)+TPI(J,L)
C     AJL(J,L,2)=AJL(J,L,2)+PHIPI(J,L)
C     AJL(J,L,3)=AJL(J,L,3)+QPI
C     AJL(J,L,17)=AJL(J,L,17)+SPI(J,L)
C     AJL(J,L,18)=AJL(J,L,18)+RHPI
  230 CONTINUE
  250 CONTINUE
C****
C**** NORTHWARD GRADIENT OF TEMPERATURE: TROPOSPHERIC AND STRATOSPHERIC
C****
      DO 385 J=2,JM-1
C**** MEAN TROPOSPHERIC NORTHWARD TEMPERATURE GRADIENT
      DO 340 L=1,LS1-1
      ADTDL=0.
      BDTDL=0.
      CDTDL=0.
      DO 335 I=1,IM
      PLAND=FLAND(I,J)
      POICE=ODATA(I,J,2)*(1.-PLAND)
      POCEAN=(1.-PLAND)-POICE
      ADTDL=ADTDL+(TX(I,J+1,L)-TX(I,J-1,L))*POCEAN
      BDTDL=BDTDL+(TX(I,J+1,L)-TX(I,J-1,L))*PLAND
      CDTDL=CDTDL+(TX(I,J+1,L)-TX(I,J-1,L))*POICE
  335 CONTINUE
  338 AJ(J,J_DTDJT)=AJ(J,J_DTDJT)+ADTDL*DSIG(L)
      BJ(J,J_DTDJT)=BJ(J,J_DTDJT)+BDTDL*DSIG(L)
  340 CJ(J,J_DTDJT)=CJ(J,J_DTDJT)+CDTDL*DSIG(L)
C**** MEAN STRATOSPHERIC NORTHWARD TEMPERATURE GRADIENT
      IF (LS1.GT.LM) GO TO 380
      DO 370 L=LS1,LM
      ADTDL=0.
      BDTDL=0.
      CDTDL=0.
      DO 350 I=1,IM
      PLAND=FLAND(I,J)
      POICE=ODATA(I,J,2)*(1.-PLAND)
      POCEAN=(1.-PLAND)-POICE
      ADTDL=ADTDL+(TX(I,J+1,L)-TX(I,J-1,L))*POCEAN
      BDTDL=BDTDL+(TX(I,J+1,L)-TX(I,J-1,L))*PLAND
      CDTDL=CDTDL+(TX(I,J+1,L)-TX(I,J-1,L))*POICE
  350 CONTINUE
  360 AJ(J,J_DTDJS)=AJ(J,J_DTDJS)+ADTDL*DSIG(L)
      BJ(J,J_DTDJS)=BJ(J,J_DTDJS)+BDTDL*DSIG(L)
  370 CJ(J,J_DTDJS)=CJ(J,J_DTDJS)+CDTDL*DSIG(L)
  380 CONTINUE
  385 CONTINUE
C****
C**** STATIC STABILITIES: TROPOSPHERIC AND STRATOSPHERIC
C****
      DO 490 J=1,JM
      IMAX=IMAXJ(J)
      DXYPJ=DXYP(J)
C**** OLD TROPOSPHERIC STATIC STABILITY
      ASS=0.
      BSS=0.
      CSS=0.
      DO 390 I=1,IMAX
      JR=JREG(I,J)
      PLAND=FLAND(I,J)
      POICE=ODATA(I,J,2)*(1.-PLAND)
      POCEAN=(1.-PLAND)-POICE
      SS=(T(I,J,LS1-1)-T(I,J,1))/(PHI(I,J,LS1-1)-PHI(I,J,1)+ZERO20)
      ASS=ASS+SS*POCEAN
      BSS=BSS+SS*PLAND
      CSS=CSS+SS*POICE
      AREG(JR,J_DTDGTR)=AREG(JR,J_DTDGTR)+SS*DXYPJ
  390 AIJ(I,J,IJ_DTDP)=AIJ(I,J,IJ_DTDP)+SS
      AJ(J,J_DTDGTR)=AJ(J,J_DTDGTR)+ASS
      BJ(J,J_DTDGTR)=BJ(J,J_DTDGTR)+BSS
      CJ(J,J_DTDGTR)=CJ(J,J_DTDGTR)+CSS
C**** OLD STRATOSPHERIC STATIC STABILITY
      ASS=0.
      BSS=0.
      CSS=0.
      DO 440 I=1,IMAX
      JR=JREG(I,J)
      PLAND=FLAND(I,J)
      POICE=ODATA(I,J,2)*(1.-PLAND)
      POCEAN=(1.-PLAND)-POICE
      SS=(T(I,J,LM)-T(I,J,LS1-1))/((PHI(I,J,LM)-PHI(I,J,LS1-1))+ZERO20)
      ASS=ASS+SS*POCEAN
      BSS=BSS+SS*PLAND
      CSS=CSS+SS*POICE
      AREG(JR,J_DTSGST)=AREG(JR,J_DTSGST)+SS*DXYPJ
  440 CONTINUE
      AJ(J,J_DTSGST)=AJ(J,J_DTSGST)+ASS
      BJ(J,J_DTSGST)=BJ(J,J_DTSGST)+BSS
      CJ(J,J_DTSGST)=CJ(J,J_DTSGST)+CSS
C****
C**** NUMBERS ACCUMULATED FOR THE RADIATION EQUILIBRIUM LAYERS
C****
      DO 470 LR=1,LM_REQ
      TRI(LR)=0.
      DO 460 I=1,IMAX
  460 TRI(LR)=TRI(LR)+RQT(I,J,LR)
  470 ASJL(J,LR,1)=ASJL(J,LR,1)+(TRI(LR)-273.16*IMAX)
      PHIRI=0.
      DO 480 I=1,IMAX
  480 PHIRI=PHIRI+(PHI(I,J,LM)+RGAS*.5*(TX(I,J,LM)+RQT(I,J,1))
     *  *LOG((SIG(LM)*PSFMPT+PTOP)/PRQ1))
      ASJL(J,1,2)=ASJL(J,1,2)+PHIRI
      PHIRI=PHIRI+RGAS*.5*(TRI(1)+TRI(2))*DLNP12
      ASJL(J,2,2)=ASJL(J,2,2)+PHIRI
      PHIRI=PHIRI+RGAS*.5*(TRI(2)+TRI(3))*DLNP23
      ASJL(J,3,2)=ASJL(J,3,2)+PHIRI
  490 CONTINUE
C****
C**** RICHARDSON NUMBER , ROSSBY NUMBER , RADIUS OF DEFORMATION
C****
C**** NUMBERS ACCUMULATED OVER THE TROPOSPHERE
      DO 506 J=2,JM
      DUDVSQ(J)=0.
      UMAX(J)=0.
      DO 504 I=1,IM
      DU=U(I,J,LS1-1)-U(I,J,1)
      DV=V(I,J,LS1-1)-V(I,J,1)
      DUDVSQ(J)=DUDVSQ(J)+(DU*DU+DV*DV)*PUV(I,J)
  504 CONTINUE
  506 CONTINUE
      DO 510 J=2,JM-1
      PIBYIM=PI(J)*BYIM
      DLNP=LOG((SIG(1)*PIBYIM+PTOP)/(SIG(LS1-1)*PIBYIM+PTOP))
      DLNS=LOG(SPI(J,LS1-1)/SPI(J,1))
      DS=SPI(J,LS1-1)-SPI(J,1)
      EL(J)=SQRT(DLNS/DLNP)
      RI(J)=DS*DLNP/(.5*(DUDVSQ(J)+DUDVSQ(J+1)))
  510 CONTINUE
      DO 515 L=1,LS1-1
      DO 514 J=2,JM
      UI(J)=0.
      DO 512 I=1,IM
  512 UI(J)=UI(J)+U(I,J,L)
  514 CONTINUE
      DO 515 J=2,JM-1
      UAMAX=ABS(UI(J)+UI(J+1))
      IF (UAMAX.GT.UMAX(J)) UMAX(J)=UAMAX
  515 CONTINUE
      DO 520 J=2,JM-1
      ROSSX=DYP(J)/(DXYP(J)*SINP(J))
      ELX=1./SINP(J)
      AJ(J,J_RICTR)=AJ(J,J_RICTR)+RI(J)*SOCEAN(J)
      BJ(J,J_RICTR)=BJ(J,J_RICTR)+RI(J)*SLAND(J)
      CJ(J,J_RICTR)=CJ(J,J_RICTR)+RI(J)*SOICE(J)
      AJ(J,J_ROSTR)=AJ(J,J_ROSTR)+UMAX(J)*SOCEAN(J)*ROSSX
      BJ(J,J_ROSTR)=BJ(J,J_ROSTR)+UMAX(J)*SLAND(J)*ROSSX
      CJ(J,J_ROSTR)=CJ(J,J_ROSTR)+UMAX(J)*SOICE(J)*ROSSX
      AJ(J,J_LTRO)=AJ(J,J_LTRO)+EL(J)*SOCEAN(J)*ELX
      BJ(J,J_LTRO)=BJ(J,J_LTRO)+EL(J)*SLAND(J)*ELX
      CJ(J,J_LTRO)=CJ(J,J_LTRO)+EL(J)*SOICE(J)*ELX
  520 CONTINUE
C**** NUMBERS ACCUMULATED OVER THE STRATOSPHERE
CNOST IF (LS1.GT.LM) GO TO 551    NEEDED FOR RUNS WITHOUT A STRATOSPHERE
      DO 532 J=2,JM
      DUDVSQ(J)=0.
      UMAX(J)=0.
  532 CONTINUE
      DO 536 J=2,JM
      DO 534 I=1,IM
      DU=U(I,J,LM)-U(I,J,LS1-1)
      DV=V(I,J,LM)-V(I,J,LS1-1)
      DUDVSQ(J)=DUDVSQ(J)+(DU*DU+DV*DV)*PUV(I,J)
  534 CONTINUE
  536 CONTINUE
      DO 540 J=2,JM-1
      PIBYIM=PI(J)*BYIM
      DLNP=LOG((SIG(LS1-1)*PIBYIM+PTOP)/(SIG(LM)*PSFMPT+PTOP))
      DLNS=LOG(SPI(J,LM)/SPI(J,LS1-1))
      DS=SPI(J,LM)-SPI(J,LS1-1)
      EL(J)=SQRT(DLNS/DLNP)
      RI(J)=DS*DLNP/(.5*(DUDVSQ(J)+DUDVSQ(J+1)))
  540 CONTINUE
      DO 545 L=LS1,LM
      DO 544 J=2,JM
      UI(J)=0.
      DO 542 I=1,IM
  542 UI(J)=UI(J)+U(I,J,L)
  544 CONTINUE
      DO 545 J=2,JM-1
      UAMAX=ABS(UI(J)+UI(J+1))
      IF (UAMAX.GT.UMAX(J)) UMAX(J)=UAMAX
  545 CONTINUE
      DO 550 J=2,JM-1
      ROSSX=DYP(J)/(DXYP(J)*SINP(J))
      ELX=1./SINP(J)
      AJ(J,J_RICST)=AJ(J,J_RICST)+RI(J)*SOCEAN(J)
      BJ(J,J_RICST)=BJ(J,J_RICST)+RI(J)*SLAND(J)
      CJ(J,J_RICST)=CJ(J,J_RICST)+RI(J)*SOICE(J)
      AJ(J,J_ROSST)=AJ(J,J_ROSST)+UMAX(J)*SOCEAN(J)*ROSSX
      BJ(J,J_ROSST)=BJ(J,J_ROSST)+UMAX(J)*SLAND(J)*ROSSX
      CJ(J,J_ROSST)=CJ(J,J_ROSST)+UMAX(J)*SOICE(J)*ROSSX
      AJ(J,J_LSTR)=AJ(J,J_LSTR)+EL(J)*SOCEAN(J)*ELX
      BJ(J,J_LSTR)=BJ(J,J_LSTR)+EL(J)*SLAND(J)*ELX
      CJ(J,J_LSTR)=CJ(J,J_LSTR)+EL(J)*SOICE(J)*ELX
  550 CONTINUE
CN551 CONTINUE
C****
C**** MEAN TROPOSPHERIC LAPSE RATES:  MOIST CONVECTIVE, ACTUAL,
C****    DRY ADIABATIC
C****
      X=RGAS*LHE*LHE/(SHA*461.5)
      DO 570 J=1,JM
      GAMM=0.
      DO 560 L=1,LS1-1
      TZL=TPI(J,L)/PI(J)+273.16
      PRT=(SIG(L)*PI(J)*BYIM+PTOP)*RGAS*TZL
      ESEPS=QSAT(TZL,ONE,LHE)
      GAMM=GAMM+(PRT+LHE*ESEPS)/(PRT+X*ESEPS/TZL)*DSIG(L)
  560 CONTINUE
      AJ(J,J_GAMM)=AJ(J,J_GAMM)+GAMM*SOCEAN(J)
      BJ(J,J_GAMM)=BJ(J,J_GAMM)+GAMM*SLAND(J)
      CJ(J,J_GAMM)=CJ(J,J_GAMM)+GAMM*SOICE(J)
      GAMX=(TPI(J,1)-TPI(J,LS1-1))/(PHIPI(J,LS1-1)-PHIPI(J,1))
      AJ(J,J_GAM)=AJ(J,J_GAM)+GAMX*SOCEAN(J)
      BJ(J,J_GAM)=BJ(J,J_GAM)+GAMX*SLAND(J)
      CJ(J,J_GAM)=CJ(J,J_GAM)+GAMX*SOICE(J)
  570 CONTINUE
C**** DRY ADIABATIC LAPSE RATE
      GAMD=.0098
      DO 580 J=1,JM
      TPIL=0.
      DO 575 L=1,LS1-1
  575 TPIL=TPIL+TPI(J,L)*DSIG(L)
      TIL(J)=TPIL/(PI(J)*(SIGE(1)-SIGE(LS1)))
  580 CONTINUE
      DO 590 J=2,JM-1
      X=SINP(J)*GRAV/(COSP(J)*RGAS*2.*DLAT)
      DT2=TIL(J+1)-TIL(J-1)
      GAMC=GAMD+X*DT2/(TIL(J)+273.16)
      AJ(J,J_GAMC)=AJ(J,J_GAMC)+GAMC*SOCEAN(J)
      BJ(J,J_GAMC)=BJ(J,J_GAMC)+GAMC*SLAND(J)
      CJ(J,J_GAMC)=CJ(J,J_GAMC)+GAMC*SOICE(J)
  590 CONTINUE
C****
C**** EASTWARD TRANSPORTS
C****
      I=IM
      DO 600 L=1,LS1-1
      DO 600 J=2,JM-1
      DO 600 IP1=1,IM
      AIJ(I,J,IJ_PUQ)=AIJ(I,J,IJ_PUQ)+(P(I,J)+P(IP1,J))*(U(I,J,L)+U(I
     *        ,J+1,L))*(Q(I,J,L)+Q(IP1,J,L))*DSIG(L)
  600 I=IP1
      I=IM
      DO 601 L=LS1,LM
      DO 601 J=2,JM-1
      DO 601 IP1=1,IM
         AIJ(I,J,IJ_PUQ)=AIJ(I,J,IJ_PUQ)+2.*PSFMPT*(U(I,J,L)+U(I,J+1
     *        ,L))*(Q(I,J,L)+Q(IP1,J,L))*DSIG(L)
  601 I=IP1
C****
C**** MOMENTUM, KINETIC ENERGY, NORTHWARD TRANSPORTS, ANGULAR MOMENTUM
C****
      DO 640 J=2,JM
      P4I=0.
      I=IM
      DO 610 IP1=1,IM
      P4=P(I,J-1)+P(IP1,J-1)+P(I,J)+P(IP1,J)
      P4I=P4I+P4
C     AIJ(I,J,IJ_P4UV)=AIJ(I,J,IJ_P4UV)+P4
      AIJ(I,J,IJ_UJET)=AIJ(I,J,IJ_UJET)+U(I,J,JET)
      AIJ(I,J,IJ_VJET)=AIJ(I,J,IJ_VJET)+V(I,J,JET)
  610 I=IP1
      APJ(J,2)=APJ(J,2)+P4I
      DO 640 L=1,LM
      PU4I=0.
      PV4I=0.
C     PWW4I=0.
C     PT16I=0.
C     PTV16I=0.
C     PZ16I=0.
C     PZV16I=0.
C     PQ16I=0.
C     PQV16I=0.
C     PWWV4I=0.
      PUV4I=0.
      I=IM
      DO 620 IP1=1,IM
      P4=P(I,J-1)+P(IP1,J-1)+P(I,J)+P(IP1,J)
      IF(L.GE.LS1) P4=PSMPT4
      IF(L.EQ.LS1) P4I=FIM*P4
      PU4I=PU4I+P4*U(I,J,L)
      PV4I=PV4I+P4*V(I,J,L)
C     PWW4I=PWW4I+P4*(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))
C     PWWV4I=PWWV4I+P4*(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))*V(I,J,L)
      PUV4I=PUV4I+P4*U(I,J,L)*V(I,J,L)
      T4=TX(I,J-1,L)+TX(IP1,J-1,L)+TX(I,J,L)+TX(IP1,J,L)
C     PT16I=PT16I+P4*T4
C     PTV16I=PTV16I+P4*T4*V(I,J,L)
      Z4=PHI(I,J-1,L)+PHI(IP1,J-1,L)+PHI(I,J,L)+PHI(IP1,J,L)
C     PZ16I=PZ16I+P4*Z4
C     PZV16I=PZV16I+P4*Z4*V(I,J,L)
C     Q4=Q(I,J-1,L)+Q(IP1,J-1,L)+Q(I,J,L)+Q(IP1,J,L)
C     PQ16I=PQ16I+P4*Q4
C     PQV16I=PQV16I+P4*Q4*V(I,J,L)
      AIJ(I,J,IJ_PEV)=AIJ(I,J,IJ_PEV)+P4*(SHA*T4+Z4)*V(I,J,L)*DSIG(L)
     *     *DXV(J)
      SP2=P(IP1,J-1)+P(IP1,J)
      IF(L.GE.LS1) SP2=2.*PSFMPT
      AIJ(IP1,J,IJ_PVQ)=AIJ(IP1,J,IJ_PVQ)+SP2
     *  *(V(I,J,L)+V(IP1,J,L))*(Q(IP1,J-1,L)+Q(IP1,J,L))*DSIG(L)
  620 I=IP1
      IF(L.GE.LS1) P4I=FIM*PSMPT4
C     AJL(J,L,4)=AJL(J,L,4)+PU4I
C     AJL(J,L,5)=AJL(J,L,5)+PV4I
C     AJL(J,L,14)=AJL(J,L,14)+(PU4I*PU4I+PV4I*PV4I)/P4I
C     AJL(J,L,15)=AJL(J,L,15)+PWW4I
C     AJL(J,L,20)=AJL(J,L,20)+PT16I*PV4I/P4I
C     AJL(J,L,21)=AJL(J,L,21)+PTV16I
C     AJL(J,L,22)=AJL(J,L,22)+PZ16I*PV4I/P4I
C     AJL(J,L,23)=AJL(J,L,23)+PZV16I
C     AJL(J,L,24)=AJL(J,L,24)+PQ16I*PV4I/P4I
C     AJL(J,L,25)=AJL(J,L,25)+PQV16I
C     AJL(J,L,26)=AJL(J,L,26)+PWW4I*PV4I/P4I
C     AJL(J,L,27)=AJL(J,L,27)+PWWV4I
      AJL(J,L,48)=AJL(J,L,48)+PU4I*PV4I/P4I
      AJL(J,L,49)=AJL(J,L,49)+PUV4I
  640 CONTINUE
C****
C**** EVEN LEVEL GEOPOTENTIALS, VERTICAL WINDS AND VERTICAL TRANSPORTS
C****
      DO 655 J=1,JM
      IMAX=IMAXJ(J)
C     PITI=0.
C     DO 648 I=1,IMAX
C 648 PITI=PITI+PIT(I,J)
      DO 655 L=1,LM-1
C     SDI=0.
C     PZI=0.
C     SDZI=0.
C     PDSE2I=0.
C     SDDS2I=0.
C     PQ2I=0.
C     SDQ2I=0.
      DO 650 I=1,IMAX
C     SDI=SDI+SD(I,J,L)
      PIJ=P(I,J)
      IF(L.GE.LS1-1) PIJ=PSFMPT
c      PE=PEDN(L+1,I,J)  ! SIGE(L+1)*PIJ+PTOP
c      PKE=PEK(L+1,I,J)  ! PE**KAPA
      PE=SIGE(L+1)*PIJ+PTOP
      PKE=PE**KAPA
      THETA=THBAR(T(I,J,L+1),T(I,J,L))
      W(I,J,L)=SD(I,J,L)*THETA*PKE/PE
C     PHIE(I,J,L)=PHI(I,J,L)+SHA*THETA*(PK(L,I,J)-PKE)
C     PZI=PZI+PHIE(I,J,L)*P(I,J)
C     SDZI=SDZI+PHIE(I,J,L)*SD(I,J,L)
C     PDSE2I=PDSE2I+(SHA*(TX(I,J,L)+TX(I,J,L+1))+2.*PHIE(I,J,L))*P(I,J)
C     SDDS2I=SDDS2I+(SHA*(TX(I,J,L)+TX(I,J,L+1))+2.*PHIE(I,J,L))*
C    *  SD(I,J,L)
C     PQ2I=PQ2I+(Q(I,J,L)*Q(I,J,L+1)/(Q(I,J,L)+Q(I,J,L+1)+ZERO20))*
C    *   P(I,J)
C     SDQ2I=SDQ2I+(Q(I,J,L)*Q(I,J,L+1)/(Q(I,J,L)+Q(I,J,L+1)+
C    *  ZERO20))*SD(I,J,L)
  650 CONTINUE
C     SDMEAN(J,L)=SDI*BYIM
C     AJL(J,L,6)=AJL(J,L,6)+SDI+DSIG(L+1)*PITI
C     AJL(J,L,34)=AJL(J,L,34)+(SDZI-PZI*SDI/PI(J))
C     AJL(J,L,30)=AJL(J,L,30)+PDSE2I*SDI/PI(J)
C     AJL(J,L,31)=AJL(J,L,31)+SDDS2I
C     AJL(J,L,32)=AJL(J,L,32)+PQ2I*SDI/PI(J)
C     AJL(J,L,33)=AJL(J,L,33)+SDQ2I
  655 CONTINUE
C****
C**** VERTICAL TRANSPORT OF KINETIC ENERGY AND ANGULAR MOMENTUM
C****
C**** FILL IN AND/OR DOUBLE SD AND SDMEAN AT THE POLES
C     DO 657 L=1,LM-1
C     SDMEAN(1,L)=2.*FIM*SDMEAN(1,L)
C     SDMEAN(JM,L)=2.*FIM*SDMEAN(JM,L)
C     SDSP=2.*SD(1,1,L)
C     SDNP=2.*SD(1,JM,L)
C     DO 657 I=1,IM
C     SD(I,1,L)=SDSP
C 657 SD(I,JM,L)=SDNP
C     DO 670 J=2,JM
C     AMA=RADIUS*OMEGA*COSV(J)
C     DO 670 L=1,LM-1
C     TKEM=0.
C     TKET=0.
C     UM=0.
C     UT=0.
C     I=IM
C     DO 660 IP1=1,IM
C     SDU=SD(I,J,L)+SD(IP1,J,L)+SD(I,J-1,L)+SD(IP1,J-1,L)
C     UE=U(I,J,L)+U(I,J,L+1)
C     TKE=UE*UE+(V(I,J,L)+V(I,J,L+1))*(V(I,J,L)+V(I,J,L+1))
C     TKEM=TKEM+TKE
C     TKET=TKET+TKE*SDU
C     UM=UM+UE
C     UT=UT+UE*SDU
C 660 I=IP1
C     AJL(J,L,35)=AJL(J,L,35)+TKET
C     AJL(J,L,36)=AJL(J,L,36)+(UT-2.*UM*(SDMEAN(J,L)+SDMEAN(J-1,L)))
C     AJL(J,L,37)=AJL(J,L,37)+(UT+4*AMA*FIM*(SDMEAN(J,L)+SDMEAN(J-1,L)))
C 670 CONTINUE
C****
C**** AVAILABLE POTENTIAL ENERGY
C****
C**** SET UP FOR CALCULATION
      DO 710 L=1,LM
  710 GMEAN(L)=0.
      DO 740 J=1,JM
      IMAX=IMAXJ(J)
      DO 720 I=1,IMAX
  720 SQRTP(I)=SQRT(P(I,J))
C**** GMEAN CALCULATED FOR EACH LAYER, THJL, THSQJL ARRAYS FILLED
      DO 730 L=1,LM
      LDN=LDNA(L)
      LUP=LUPA(L)
      THJL(J,L)=0.
      THSQJL(J,L)=0.
      DO 730 I=1,IMAX
      THJL(J,L)=THJL(J,L)+T(I,J,L)*SQRTP(I)
      THSQJL(J,L)=THSQJL(J,L)+T(I,J,L)*T(I,J,L)*P(I,J)
  730 GMEAN(L)=GMEAN(L)+(SIG(L)*P(I,J)+PTOP)*(T(I,J,LUP)-T(I,J,LDN))*
     *  DXYP(J)/(P(I,J)*PK(L,I,J))
  740 CONTINUE
C**** CALCULATE APE
      DO 760 L=1,LM
      LP1=LUPA(L)
      LM1=LDNA(L)
      THJL(1,L)=THJL(1,L)*FIM
      THJL(JM,L)=THJL(JM,L)*FIM
      THSQJL(1,L)=THSQJL(1,L)*FIM
      THSQJL(JM,L)=THSQJL(JM,L)*FIM
      THGM=0.
      DO 750 J=1,JM
  750 THGM=THGM+THJL(J,L)*DXYP(J)
      THGM=THGM/AREAG
      GMEANL=GMEAN(L)/((SIG(LM1)-SIG(LP1))*AREAG)
      DO 760 J=1,JM
  760 AJL(J,L,16)=AJL(J,L,16)+(THSQJL(J,L)-2.*THJL(J,L)*THGM+THGM*THGM*
     *  FIM)/GMEANL
C****
C**** OMEGA'*ALPHA' ;  BAROCLINIC EKE GENERATION
C****
C     DO 770 L=1,LM
C     DO 770 J=2,JM-1
C     PWAI=0.
C     SPAI=0.
C     PWI=0.
C     IM1=IM
C     DO 766 I=1,IM
C     PL=SIG(L)*P(I,J)+PTOP
C     SPA=SIG(L)*P(I,J)*RGAS*TX(I,J,L)/PL
C     SVDX=(V(IM1,J+1,L)+V(I,J+1,L))*DXV(J+1)-(V(IM1,J,L)+V(I,J,L))*
C    *   DXV(J)
C     SUDY=(U(I,J+1,L)+U(I,J,L)-U(IM1,J+1,L)-U(IM1,J,L))*DYP(J)
C     PWA=-.5*SPA*(SUDY+SVDX)*DSIG(L)*P(I,J)
C     IF (L.NE.LM) PWA=PWA+SD(I,J,L)*((PHIE(I,J,L)-PHI(I,J,L))
C    *    +SPA)
C     IF (L.NE.1) PWA=PWA+SD(I,J,L-1)*((PHI(I,J,L)-PHIE(I,J,L-1))
C    *    -SPA)
C     PWAI=PWAI+PWA
C     SPAI=SPAI+SPA
C     PWI=PWI+PWA*P(I,J)/SPA
C 766 IM1=I
C 770 AJL(J,L,7)=AJL(J,L,7)-(PWAI-PWI*SPAI/PI(J))
C****
C**** P-K BY PRESSURE GRADIENT FORCE
C****
C     DO 774 L=1,LM
C     DO 774 J=2,JM
C     PDA4I=0.
C     VPDA4I=0.
C     DVTI=0.
C     VDVTI=0.
C     DUTI=0.
C     UDUTI=0.
C     UPDA4I=0.
C     I=IM
C     DO 772 IP1=1,IM
C     PDA4=(P(I,J)+P(IP1,J))*DXYP(J) + (P(I,J-1)+P(IP1,J-1))*DXYP(J-1)
C     PDA4I=PDA4I+PDA4
C     VPDA4I=VPDA4I+V(I,J,L)*PDA4
C     DVTI=DVTI+DVT(I,J,L)
C     VDVTI=VDVTI+V(I,J,L)*DVT(I,J,L)
C     DUTI=DUTI+DUT(I,J,L)
C     UDUTI=UDUTI+U(I,J,L)*DUT(I,J,L)
C     UPDA4I=UPDA4I+U(I,J,L)*PDA4
C 772 I=IP1
C 774 AJL(J,L,53)=AJL(J,L,53)+VDVTI+UDUTI
C    *   -(UPDA4I*DUTI+VPDA4I*DVTI)/PDA4I
C        IF (JM.NE.24) GO TO 850
C****
C**** CERTAIN HORIZONTAL WIND AVERAGES
C****
      DO 790 L=1,LM
      DO 780 J=2,JM
      DO 776 I=I135W,I110W                ! EAST PACIFIC
      AJL(J,L,41)=AJL(J,L,41)+U(I,J,L)
  776 AJL(J,L,42)=AJL(J,L,42)+V(I,J,L)
      DO 778 I=I150E,IM                   ! WEST PACIFIC
      AJL(J,L,44)=AJL(J,L,44)+U(I,J,L)
  778 AJL(J,L,45)=AJL(J,L,45)+V(I,J,L)
  780 CONTINUE
      DO 785 J=J5S,J5N
      DO 785 I=1,IM
      PIJ=P(I,J)
      IF(L.GE.LS1) PIJ=PSFMPT
      AIL(I,L,1)=AIL(I,L,1)+U(I,J,L)
      AIL(I,L,2)=AIL(I,L,2)+V(I,J,L)
      AIL(I,L,4)=AIL(I,L,4)+(TX(I,J,L)-273.16)
  785 AIL(I,L,5)=AIL(I,L,5)+Q(I,J,L)/QSAT(TX(I,J,L),SIG(L)*PIJ+PTOP,LHE)
      DO 788 I=1,IM
      AIL(I,L,10)=AIL(I,L,10)+(TX(I,J50N,L)-273.16)
      AIL(I,L,12)=AIL(I,L,12)+(U(I,J50N,L)+U(I,J50N+1,L))
      AIL(I,L,14)=AIL(I,L,14)+(TX(I,J70N,L)-273.16)
      AIL(I,L,16)=AIL(I,L,16)+(U(I,J70N,L)+U(I,J70N+1,L))
  788 CONTINUE
  790 CONTINUE
C****
C**** CERTAIN VERTICAL WIND AVERAGES
C****
      DO 820 L=1,LM-1
      DO 820 J=2,JM-1
      DO 810 I=I135W,I110W                ! EAST PACIFIC
  810 AJL(J,L,43)=AJL(J,L,43)+W(I,J,L)
      DO 815 I=I150E,IM                   ! WEST PACIFIC
  815 AJL(J,L,46)=AJL(J,L,46)+W(I,J,L)
  820 CONTINUE
      DO 840 L=1,LM-1
      DO 840 I=1,IM
      DO 835 J=J5SUV,J5NUV                ! +/1 5 DEG (APPROX.)
  835 AIL(I,L, 3)=AIL(I,L, 3)+W(I,J,L)
      AIL(I,L, 9)=AIL(I,L, 9)+W(I,J50N,L)
      AIL(I,L,13)=AIL(I,L,13)+W(I,J70N,L)
  840 CONTINUE
C****
C**** ELIASSEN PALM FLUX
C****
C**** NORTHWARD TRANSPORT
      DO 868 J=2,JM
      BYDXYV=1./DXYV(J)
      I=IM
      DO 862 IP1=1,IM
      PDA(I)=.5*((P(I,J)+P(IP1,J))*DXYS(J)+(P(I,J-1)+P(IP1,J-1))*
     *  DXYN(J-1))
      PSEC(I)=PDA(I)*BYDXYV
  862 I=IP1
      DO 868 L=1,LM
      DUDP=0.
      DTHDP=0.
      UMN=0.
      THMN=0.
      LDN=LDNA(L)
      LUP=LUPA(L)
      I=IM
      DO 864 IP1=1,IM
      DUDP=DUDP+U(I,J,LUP)-U(I,J,LDN)
      DTHDP=DTHDP+T(I,J,LUP)+T(I,J-1,LUP)-T(I,J,LDN)-T(I,J-1,LDN)
      UMN=UMN+U(I,J,L)
      THMN=THMN+T(I,J,L)+T(I,J-1,L)
      THSEC(I)=T(I,J,L)+T(IP1,J,L)+T(I,J-1,L)+T(IP1,J-1,L)
  864 I=IP1
      UMN=UMN*BYIM
      THMN=2.*THMN/FIM
      FPHI=0.
      SMALL=.0002*FIM*T(1,J,L)
      IF (DTHDP.LT.SMALL) WRITE (6,999) J,L,DTHDP,SMALL
      IF (DTHDP.LT.SMALL) DTHDP=SMALL
      DO 866 I=1,IM
      SP=PSEC(I)
      IF(L.GE.LS1) SP=PSFMPT
  866 FPHI=FPHI+SP*V(I,J,L)*(.5*(THSEC(I)-THMN)*DUDP/DTHDP
     *   -U(I,J,L)+UMN)
  868 AJL(J,L,37)=AJL(J,L,37)+FPHI
C**** VERTICAL TRANSPORT
      DO 878 J=2,JM-1
      PITMN=0.
      DO 870 I=1,IM
  870 PITMN=PITMN+PIT(I,J)
      PITMN=PITMN/FIM
      DO 878 L=1,LM-1
      IF(L.GE.LS1-1) PITMN=0.
      THMN=0.
      SDMN=0.
      DTHDP=0.
      DO 872 I=1,IM
      DTHDP=DTHDP+T(I,J,L+1)-T(I,J,L)
      THMN=THMN+T(I,J,L+1)+T(I,J,L)
  872 SDMN=SDMN+SD(I,J,L)
      SMALL=.0001*FIM*T(1,J,L+1)
      IF (DTHDP.LT.SMALL) WRITE (6,999) J,L,DTHDP,SMALL
      IF (DTHDP.LT.SMALL) DTHDP=SMALL
      THMN=THMN/FIM
      SDMN=SDMN/FIM
      DUDX=0.
      PVTHP=0.
      SDPU=0.
      IM1=IM
      DO 874 I=1,IM
      DUDX=DUDX+DXV(J+1)*(U(I,J+1,L)+U(I,J+1,L+1))-DXV(J)*
     *   (U(I,J,L)+U(I,J,L+1))
      UPE=U(IM1,J,L)+U(IM1,J+1,L)+U(I,J,L)+U(I,J+1,L)+
     *    U(IM1,J,L+1)+U(IM1,J+1,L+1)+U(I,J,L+1)+U(I,J+1,L+1)
      VPE=V(IM1,J,L)+V(IM1,J+1,L)+V(I,J,L)+V(I,J+1,L)+
     *    V(IM1,J,L+1)+V(IM1,J+1,L+1)+V(I,J,L+1)+V(I,J+1,L+1)
      DP=DSIGO(L)*P(I,J)
      IF(L.GE.LS1) DP=DSIGO(L)*PSFMPT
      IF(L.EQ.LS1-1) DP=P(I,J)*SIG(LS1-1)-PSFMPT*SIG(LS1)
      PVTHP=PVTHP+DP*VPE*(T(I,J,L)+T(I,J,L+1)-THMN)
      PITIJ=PIT(I,J)
      IF(L.GE.LS1-1) PITIJ=0.
      SDPU=SDPU+(SD(I,J,L)-SDMN+(PITIJ-PITMN)*SIGE(L+1))*UPE
  874 IM1=I
      AJL(J,L,36)=AJL(J,L,36)+(.5*FIM*FCOR(J)-.25*DUDX)*PVTHP
     *   /DTHDP + SDPU
  878 CONTINUE
C****
C**** POTENTIAL VORTICITY  (10/30/81)
C****
C     DO 895 L=1,LM
C     LUP=LUPA(L)
C     LDN=LDNA(L)
C     DO 895 J=2,JM-1
C     PV2I=0.
C     I=IM
C     DO 890 IP1=1,IM
C     DS2=T(I,J,LUP)+T(IP1,J,LUP)-T(I,J,LDN)-T(IP1,J,LDN)
C     PV2I=PV2I+DS2*(U(I,J+1,L)*DXV(J+1)-U(I,J,L)*DXV(J) - FCOR(J))
C 890 I=IP1
C 895 AJL(J,L,52)=AJL(J,L,52)+PV2I
C****
C**** TRANSFORMED STREAM FUNCTION
C****
C---- skip for now lines 1132-1163
C**** ACCUMULATE TIME USED IN DIAGA
      KEND = MCLOCK ()
      MDIAG = MDIAG + (KEND-KBEGIN)
      MDYN  = MDYN  - (KEND-KBEGIN)
C     WRITE (6,997) IDACC,MINC,MDIAG
      RETURN
  997 FORMAT (' DIAGNOSTICS ACCUMULATED ',12I4,15X,2I7)
  999 FORMAT (' DTHETA/DP IS TOO SMALL AT J=',I4,' L=',I4,2F15.6)
      END SUBROUTINE DIAGA

      SUBROUTINE DIAGB (U,V,T,P,Q,WM,DUT,DVT)
C****
C**** CONTENTS OF AJK(J,K,N)  (SUM OVER LONGITUDE AND TIME OF)
C**CP   1  DP  (PDN-PM(K+1);  PDN=MAX(PM(K+1),PS)
C**CP   2  DP4  (PDN-PM(K+1))   (UV GRID)
C***1   3  (TX-273.16)*DP
C***2   4  PHI*DP
C***3   5  Q*DP
C**17   6  TH*DP
C**18   7  RH*DP
C***4   8  U*DP4 (100 PA*M/S)  (UV GRID)
C***5   9  V*DP4 (100 PA*M/S)  (UV GRID)
C**14  10  (PUI*PUI+PVI*PVI)/DPI (100 N/S**2)  (UV GRID) ..I=SUM OVER I
C**15  11  (U*DP4*U*DP4+V*DP4*V*DP4)/DP4 (100 N/S**2)  (UV GRID)
C**20  12  4*PT4I*PVI/DPI (100 PA*K*M/S)  (UV GRID)
C**21  13  4*PTV4I (100 PA*K*M/S)  (UV GRID)
C**22  14  4*PZ4I*PVI/DPI (100 W/S**2)  (UV GRID)
C**23  15  4*PZV4I (100 W/S**2)  (UV GRID)
C**24  16  4*PQ4I*PVI/DPI (100 PA*M/S)  (UV GRID)
C**25  17  4*PQV4I (100 PA*M/S)  (UV GRID)
C**26  18  PWWI*PVI/DPI (100 W/S**2)  (UV GRID)
C**27  19  PWWVI (100 W/S**2)  (UV GRID)
C**48  20  PUI*PVI/DPI (100 N/S**2)  (UV GRID)
C**49  21  PUVI (100 N/S**2)  (UV GRID)
C**53  22  VDVT - VDPI*DVT/DPI
C****  23  DP**2  (UV GRID)
C****  24  IM  (UV GRID)
C***6  25  W*DA  (100 NT/S)
C**30  26  WI*(SHA*TI+ZI)/FIM
C**31  27  SHA*WTI+WZI
C**32  28  WI*QI/FIM
C**33  29  WQI
C**34  30  WZI-WI*ZI/FIM
C***7  31  2*(WPA2I-W2I*PAI/FIM)  (ODD LAYERS)
C**52  32  PV = STB*(D(UDX)-F*DA)
C****  33  WPV4I
C****  33  WPV4I - W2I*PV2I/FIM
C****  35  IM  (PT GRID)
C**35  36  WKE4I  (UV GRID)
C**36  37  WU4I - W4I*UKI/FIM  (UV GRID)
C**37  38  WU4I+W4I*UEARTH  (UV GRID)
C****  39     4*(PSV4I-PS4I*PVI/DPI)/STB (MB*MB M/S)
C****  40     ADVU=-(DUDY-F)*V-DUDP*W  (M/S/S)
C****  41     ADVT=-DTDY*V-DTDP*W   (DEGK/S)
C****  42     LADVU=-(DUDY-F)*V* - DUDP*W*  (M/S/S)
C****  43     LADVT=-DTDY*V* - DTDP*W*  (DEGK/S)
C****  44     FY= V'T'/DTDP*DUDP - U'V'
C****  45     FP= -U'W' - V'T'/DTDP*(DUDY-F)
C****  46     U1 = U AT IDACC(4)=1  (M/S) (NOT SUM OVER TIME)
C****  47     U-U1   (M/S)                (NOT SUM OVER TIME)
C****  48     T1 = T AT IDACC(4)=1  (DEG K  THETA) (NOT SUM OVER TIME)
C****  49     T-T1    (DEG K   THETA)              (NOT SUM OVER TIME)
C****  50     W'TH'      (DEG K-MB/S)
C****
C**** CONTENTS OF AIJK(I,J,K,N)   (SUM OVER TIME OF)
C****   1  DP4*U (100 PA*M/S)  (UV GRID)
C****   2  DP4*V (100 PA*M/S)  (UV GRID)
C****   3  4*DP4*(SHA*T4+Z4) (100 N/S**2)  (UV GRID)
C****   4  DP4 (100 PA)  (UV GRID)
C****   5  4*DP4*T4 (100 K*PA)  (UV GRID)
C****   6  4*DP4*Q4 (100 PA)  (UV GRID)
C****
      USE CONSTANT, only : lhe,omega,sha
      USE E001M12_COM, only :
     &     im,imh,fim,byim,jm,jeq,lm,ls1,idacc,psf,ptop,psfmpt,
     &     mdyn,mdiag, ndyn,ndaa, skipse,
     &     sig,sige,dsig, tofday, ijd6
      USE GEOM, only :
     &     COSV,DXV,DXYN,DXYP,DXYS,DXYV,DYP,DYV,FCOR,IMAXJ,RADIUS
      USE DAGCOM, only : ajk,aijk,aijl,ajlsp,speca,adaily,nspher,
     &     nwav_dag,ndlypt,hr_in_day
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: U,V,T,Q,WM
      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: DUT,DVT
      DOUBLE PRECISION, DIMENSION(IM,JM) :: P
      DOUBLE PRECISION, DIMENSION(IMH+1,NSPHER) :: KE
      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: W,ZX,STB
      COMMON/WORK2/W,ZX,STB
      DOUBLE PRECISION, DIMENSION(JM,LM) ::
     &     STJK,DPJK,UJK,VJK,WJK,TJK,
     &     PSIJK,UP,TY,PSIP,WTJK,UVJK,WUJK
      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: PHI,TX ! from diaga
      COMMON/WORK3/PHI,TX
      DOUBLE PRECISION, DIMENSION(IM) :: PSEC,X1
      DOUBLE PRECISION, DIMENSION(LM) :: SHETH,PMO,PLO,DPM,DTH
      DOUBLE PRECISION, DIMENSION(LM+1) :: PM,PL
      DOUBLE PRECISION, DIMENSION(0:IMH,JM,LM):: FCJKA,FCJKB,FCPVA,FCPVB
      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: AIJL1,AIJL2,AIJL3,AIJL4

      INTEGER :: MCLOCK ! external
      INTEGER ::
     &     I,IH,IHOUR,IM1,IMAX,INCH,INCHM,IP1,IZERO,J,J45N,
     &     JET,JHEMI,K,KDN,KM,KMM1,KR,KS,KS1,KSPHER,KUP,KX,L,LMP1,
     &     LUP,MBEGIN,MEND,N,NM

      DOUBLE PRECISION ::
     &     BYDP,BYFIM,DELP,DP,DPDN,DPDX,DPDY,
     &     DPE,DPI,DPK,DPSQI,DPUP,DPUV,DUTI,DUTK,DVTI,DVTK,FIMI,
     &     PAI,PAK,PDN,PHIPI,PMK,PQ4I,PQ4K,PQ4L,PQV4I,PS,PS4I,
     &     PS4K,PSIY,PSV4I,PT4I,PT4K,PT4L,PTK,PTV4I,PUI,PUK,PUP,
     &     PUVI,PV2,PV2I,PVI,PVK,PWWI,PWWVI,PY,PZ4I,PZ4K,PZ4L,PZM,
     &     PZV4I,QK,QKI,QLH,QPI,QSATL,RHPI,
     &     SMALL,SP,SP2,SQRTDP,THK,THKI,THPI,TK,TKI,TPI,
     &     UDUTI,UDX,UEARTH,UK,UKI,UY,VDVTI,VK,VSTAR,W2,W2I,W4,
     &     W4I,WI,WKE4I,WMPI,WNP,WPA2I,WPV4I,WQI,WSP,WSTAR,WTHI,
     &     WTI,WU4I,WUP,WZI,ZK,ZKI

      INTEGER :: IFIRST
      DATA IFIRST/1/
      DOUBLE PRECISION, PARAMETER :: ZERO20=1.E-20,BIG=1.E20
C**** QSAT=(RAIR/RVAPOR)*6.1071*EXP((L/RVAPOR)*(1/TF-1/T))/P
      DOUBLE PRECISION, PARAMETER ::
     &     AQSAT=3.797915,BQSAT=7.93252E-6,CQSAT=2.166847E-3
      DOUBLE PRECISION :: QSAT,TM,PR,QL
      QSAT(TM,PR,QL)=AQSAT*EXP(QL*(BQSAT-CQSAT/TM))/PR

      MBEGIN = MCLOCK ()
      IF (IFIRST.NE.1) GO TO 50
      IFIRST=0
C**** INITIALIZE CERTAIN QUANTITIES
      LMP1=LM+1
      JET=LS1-1
c      BYIM=1./FIM
c      SHA=RGAS/KAPA
      KM=LM
      KMM1=KM-1
      PM(1)=1200.
      DO 20 L=2,LM+1
      PL(L)=PSFMPT*SIGE(L)+PTOP
   20 PM(L)=PSFMPT*SIGE(L)+PTOP
      DO 30 L=1,LM
      PLO(L)=PSFMPT*SIG(L)+PTOP
   30 PMO(L)=.5*(PM(L)+PM(L+1))
   50 CONTINUE
C****
C**** INTERNAL QUANTITIES T,TH,Q,RH
C****
      QLH=LHE
      DO 170 J=1,JM
      IMAX=IMAXJ(J)
      DO 170 K=1,KM
      DPI=0.
      TPI=0.
      PHIPI=0.
      QPI=0.
      WMPI=0.
      THPI=0.
      RHPI=0.
      FIMI=0.
      DO 160 I=1,IMAX
C**** FIND L=L(K) AND LUP=L(K+1) S.T. P(LUP).GT.P(K+1)
      SP=P(I,J)
      IF(K.GE.LS1) SP=PSFMPT
      PS=SP+PTOP
      IF (PM(K+1).GE.PS) GO TO 160
      L=1
      PDN=PS
      IF (PM(K).GE.PS) GO TO 120
      PDN=PM(K)
  110 IF (PM(K).GT.SP*SIGE(L+1)+PTOP) GO TO 120
      L=L+1
      GO TO 110
  120 LUP=L
  130 IF (PM(K+1).GE.SP*SIGE(LUP+1)+PTOP) GO TO 140
      LUP=LUP+1
      GO TO 130
  140 CONTINUE
C**** ACCUMULATE HERE
      DPI=DPI+PDN-PM(K+1)
      FIMI=FIMI+1.
  150 PUP=SP*SIGE(L+1)+PTOP
      IF (LUP.EQ.L) PUP=PM(K+1)
      DP=PDN-PUP
      TPI=TPI+(TX(I,J,L)-273.16)*DP
      PHIPI=PHIPI+PHI(I,J,L)*DP
      QPI=QPI+Q(I,J,L)*DP
      WMPI=WMPI+WM(I,J,L)*DP
CW       IF(WMPI.GT.1.E-3) WRITE(6,169) I,J,L,DP,WM(I,J,L),WMPI
CW169 FORMAT(1X,'1616--',3I5,3E15.2)
      THPI=THPI+T(I,J,L)*DP
      QSATL=QSAT(TX(I,J,L),SIG(L)*SP+PTOP,QLH)
      IF (QSATL.GT.1.) QSATL=1.
      RHPI=RHPI+Q(I,J,L)*DP/QSATL
      IF (L.EQ.LUP) GO TO 160
      L=L+1
      PDN=SP*SIGE(L)+PTOP
      GO TO 150
  160 CONTINUE
      AJK(J,K,35)=AJK(J,K,35)+FIMI
      AJK(J,K,1)=AJK(J,K,1)+DPI
      AJK(J,K,3)=AJK(J,K,3)+TPI
      AJK(J,K,4)=AJK(J,K,4)+PHIPI
      AJK(J,K,5)=AJK(J,K,5)+QPI
      AJK(J,K,6)=AJK(J,K,6)+THPI
      AJK(J,K,7)=AJK(J,K,7)+RHPI
      AJK(J,K,51)=AJK(J,K,51)+WMPI
         TJK(J,K)=THPI/(DPI+ZERO20)
         IF (IDACC(4).EQ.1) AJK(J,K,48)=TJK(J,K)
         AJK(J,K,49)=TJK(J,K)-AJK(J,K,48)
  170 CONTINUE
C****
C**** CALCULATE STABILITY AT ODD LEVELS ON PU GRID
C****
      DO 230 J=1,JM
      IMAX=IMAXJ(J)
      I=IMAX
      DO 230 IP1=1,IMAX
      SP2=P(I,J)+P(IP1,J)
      SP=.5*SP2
      DO 175 L=1,LS1-1
      PLO(L)=SP*SIG(L)+PTOP
  175 PL(L)=SP*SIGE(L)+PTOP
      DO 180 L=1,LM-1
      DTH(L)=(T(I,J,L)+T(IP1,J,L)-T(I,J,L+1)-T(IP1,J,L+1))/
     *  (2.*(PLO(L)-PLO(L+1)))
  180 CONTINUE
      DO 220 K=1,KM
      STB(I,J,K)=0.
      IF (PM(K+1).GE.PL(1)) GO TO 220
      PMK=PMO(K)
      IF (PM(K).GT.PL(1)) PMK=.5*(SP+PTOP+PM(K+1))
      L=2
      IF (PMK.GE.PL(2)) GO TO 210
  190 LUP=L+1
      IF (L.EQ.LM) GO TO 210
      IF (PMK.GE.PL(LUP)) GO TO 200
      L=LUP
      GO TO 190
  200 DPUP=PMK-PL(LUP)
      DPDN=PL(L)-PMK
      STB(I,J,K)=(DTH(L-1)*DPUP+DTH(L)*DPDN)/(DPUP+DPDN+ZERO20)
      GO TO 220
C**** SPECIAL CASES,  L=2, L=LM
  210 STB(I,J,K)=DTH(L-1)
  220 CONTINUE
  230 I=IP1
C**** CALCULATE STJK; THE MEAN STATIC STABILITY
      DO 260 J=1,JM
      IMAX=IMAXJ(J)
      DO 260 K=1,KM
      STJK(J,K)=0.
      DPJK(J,K)=0.
      I=IMAX
      DO 250 IP1=1,IMAX
      PS=.5*(P(I,J)+P(IP1,J))+PTOP
      IF (PM(K+1).GT.PS) GO TO 250
      STJK(J,K)=STJK(J,K)+STB(I,J,K)
      DPJK(J,K)=DPJK(J,K)+1.
  250 I=IP1
      STJK(J,K)=STJK(J,K)/(DPJK(J,K)+ZERO20)
      SMALL=.0001
      IF (ABS(STJK(J,K)).LT.SMALL) STJK(J,K)=-SMALL
  260 CONTINUE
C****
C**** CONSTANT PRESSURE DIAGNOSTICS:  FLUX, ENERGY, ANGULAR MOMENTUM
C****
      DO 390 J=2,JM
      PZM=0.
      I=IM
      DO 280 IP1=1,IM
      PSEC(I)=.25*(P(I,J-1)+P(IP1,J-1)+P(I,J)+P(IP1,J))
      PZM=PZM+PSEC(I)
      DO 270 K=1,KM
  270 ZX(I,J,K)=0.
      DO 275 L=1,LS1-1
      DUT(I,J,L)=DUT(I,J,L)/(PSEC(I)*DXYV(J)*DSIG(L))
  275 DVT(I,J,L)=DVT(I,J,L)/(PSEC(I)*DXYV(J)*DSIG(L))
      DO 276 L=LS1,LM
      DUT(I,J,L)=DUT(I,J,L)/(PSFMPT*DXYV(J)*DSIG(L))
  276 DVT(I,J,L)=DVT(I,J,L)/(PSFMPT*DXYV(J)*DSIG(L))
  280 I=IP1
C**** ACCUMULATE AIJL ARRAYS
      PZM=PZM/FIM
      DO 285 L=1,LM
      I=IM
      IF(L.EQ.LS1) PZM=PSFMPT
      DELP=PZM*DSIG(L)
      DO 284 IP1=1,IM
      PT4L=DELP*(TX(I,J-1,L)+TX(IP1,J-1,L)+TX(I,J,L)+TX(IP1,J,L))
      PZ4L=DELP*(PHI(I,J-1,L)+PHI(IP1,J-1,L)+PHI(I,J,L)+PHI(IP1,J,L))
      PQ4L=DELP*(Q(I,J-1,L)+Q(IP1,J-1,L)+Q(I,J,L)+Q(IP1,J,L))
      AIJL(I,J,L,1)=AIJL(I,J,L,1)+DELP*U(I,J,L)
      AIJL(I,J,L,2)=AIJL(I,J,L,2)+DELP*V(I,J,L)
      AIJL(I,J,L,3)=AIJL(I,J,L,3)+SHA*PT4L+PZ4L
      AIJL(I,J,L,4)=AIJL(I,J,L,4)+PQ4L
      AIJL(I,J,L,5)=AIJL(I,J,L,5)+DELP
      AIJL1(I,J,L)=DELP*U(I,J,L)
      AIJL2(I,J,L)=V(I,J,L)
      AIJL3(I,J,L)=SHA*PT4L+PZ4L
      AIJL4(I,J,L)=PQ4L
  284 I=IP1
  285 CONTINUE
      DO 350 K=1,KM
      DPI=0.
      DPSQI=0.
      FIMI=0.
      PUI=0.
      PVI=0.
      PWWI=0.
      PT4I=0.
      PTV4I=0.
      PZ4I=0.
      PZV4I=0.
      PQ4I=0.
      PQV4I=0.
      PWWVI=0.
      PUVI=0.
      DVTI=0.
      VDVTI=0.
      DUTI=0.
      UDUTI=0.
      PS4I=0.
      PSV4I=0.
      I=IM
      DO 340 IP1=1,IM
      SP=PSEC(I)
      PS=SP+PTOP
      DO 286 L=1,LS1-1
  286 PL(L)=SP*SIGE(L)+PTOP
      IF (PM(K+1).GE.PS) GO TO 336
      L=1
      PDN=PS
      IF (PM(K).GE.PS) GO TO 300
      PDN=PM(K)
  290 IF (PM(K).GT.PL(L+1)) GO TO 300
      L=L+1
      GO TO 290
  300 LUP=L
  310 IF (PM(K+1).GE.PL(LUP+1)) GO TO 320
      LUP=LUP+1
      GO TO 310
  320 CONTINUE
      DPK=PDN-PM(K+1)
      PUK=0.
      PVK=0.
      PT4K=0.
      PZ4K=0.
      PQ4K=0.
      DUTK=0.
      DVTK=0.
      PS4K=0.
C**** INTERPOLATE HERE
  330 PUP=PL(L+1)
      IF (LUP.EQ.L) PUP=PM(K+1)
      DP=PDN-PUP
      PUK=PUK+DP*U(I,J,L)
      PVK=PVK+DP*V(I,J,L)
      PT4K=PT4K+DP*(TX(I,J-1,L)+TX(IP1,J-1,L)+TX(I,J,L)+TX(IP1,J,L))
      PZ4K=PZ4K+DP*(PHI(I,J-1,L)+PHI(IP1,J-1,L)+PHI(I,J,L)+PHI(IP1,J,L))
      PQ4K=PQ4K+DP*(Q(I,J-1,L)+Q(IP1,J-1,L)+Q(I,J,L)+Q(IP1,J,L))
      DUTK=DUTK+DP*DUT(I,J,L)
      DVTK=DVTK+DP*DVT(I,J,L)
      PS4K=PS4K+DP*(T(I,J-1,L)+T(IP1,J-1,L)+T(I,J,L)+T(IP1,J,L))
      IF (LUP.EQ.L) GO TO 332
      L=L+1
      PDN=PL(L)
      GO TO 330
C**** ACCUMULATE HERE
  332 FIMI=FIMI+1.
      DPI=DPI+DPK
      DPSQI=DPSQI+DPK*DPK
      IF (DPK.LT.ZERO20) DPK=ZERO20
      BYDP=1./DPK
      PUI=PUI+PUK
      PVI=PVI+PVK
      PWWI=PWWI+BYDP*(PUK*PUK+PVK*PVK)
      PWWVI=PWWVI+BYDP*BYDP*(PUK*PUK+PVK*PVK)*PVK
      PUVI=PUVI+BYDP*PUK*PVK
      PT4I=PT4I+PT4K
      PTV4I=PTV4I+BYDP*PT4K*PVK
      PZ4I=PZ4I+PZ4K
      PZV4I=PZV4I+BYDP*PZ4K*PVK
      PQ4I=PQ4I+PQ4K
      PQV4I=PQV4I+BYDP*PQ4K*PVK
      DVTI=DVTI+DVTK
      VDVTI=VDVTI+BYDP*PVK*DVTK
      DUTI=DUTI+DUTK
      UDUTI=UDUTI+BYDP*PUK*DUTK
      IF(SKIPSE.EQ.1.) GO TO 334
      AIJK(I,J,K,1)=AIJK(I,J,K,1)+PUK
      AIJK(I,J,K,2)=AIJK(I,J,K,2)+PVK
      AIJK(I,J,K,3)=AIJK(I,J,K,3)+SHA*PT4K+PZ4K
      AIJK(I,J,K,4)=AIJK(I,J,K,4)+DPK
      AIJK(I,J,K,5)=AIJK(I,J,K,5)+PT4K
      AIJK(I,J,K,6)=AIJK(I,J,K,6)+PQ4K
C**** EDDY TRANSPORT OF THETA;  VORTICITY
  334 PS4I=PS4I+PS4K
      PSV4I=PSV4I+BYDP*PVK*PS4K
      UDX=BYDP*PUK*DXV(J)
      ZX(I,J,K)=-UDX
      IF (ZX(I,J-1,K).LT.BIG) ZX(I,J-1,K)=ZX(I,J-1,K)+UDX
      IF (ZX(I,J-1,K).GE.BIG) ZX(I,J-1,K)=0.
      GO TO 340
  336 ZX(I,J,K)=BIG
      ZX(I,J-1,K)=0.
  340 I=IP1
      DPM(K)=DPI/(FIMI+ZERO20)
      DPJK(J,K)=DPI
      AJK(J,K,2)=AJK(J,K,2)+DPI
      AJK(J,K,23)=AJK(J,K,23)+DPSQI
      AJK(J,K,24)=AJK(J,K,24)+FIMI
      IF (DPI.LT.ZERO20) DPI=ZERO20
      AJK(J,K,8)=AJK(J,K,8)+PUI
      AJK(J,K,9)=AJK(J,K,9)+PVI
      AJK(J,K,10)=AJK(J,K,10)+(PUI*PUI+PVI*PVI)/DPI
      AJK(J,K,11)=AJK(J,K,11)+PWWI
      AJK(J,K,12)=AJK(J,K,12)+PT4I*PVI/DPI
      AJK(J,K,13)=AJK(J,K,13)+PTV4I
      AJK(J,K,14)=AJK(J,K,14)+PZ4I*PVI/DPI
      AJK(J,K,15)=AJK(J,K,15)+PZV4I
      AJK(J,K,16)=AJK(J,K,16)+PQ4I*PVI/DPI
      AJK(J,K,17)=AJK(J,K,17)+PQV4I
      AJK(J,K,18)=AJK(J,K,18)+PWWI*PVI/DPI
      AJK(J,K,19)=AJK(J,K,19)+PWWVI
      AJK(J,K,20)=AJK(J,K,20)+PUI*PVI/DPI
      AJK(J,K,21)=AJK(J,K,21)+PUVI
      AJK(J,K,22)=AJK(J,K,22)+VDVTI+UDUTI-
     *   (PUI*DUTI+PVI*DVTI)/DPI
      SHETH(K)=(PSV4I-PS4I*PVI/DPI)*DXYV(J)/(STJK(J-1,K)*DXYN(J-1)+
     *   STJK(J,K)*DXYS(J))
         UJK(J,K)=PUI/DPI
         VJK(J,K)=PVI/DPI
         PSIJK(J,K)=.25*SHETH(K)/DPI
         UVJK(J,K)=(PUVI-PUI*PVI/DPI)/DPI
         IF (IDACC(4).EQ.1) AJK(J,K,46)=UJK(J,K)
         AJK(J,K,47)=UJK(J,K)-AJK(J,K,46)
  350 AJK(J,K,39)=AJK(J,K,39)+SHETH(K)
      DO 345 L=1,LM
C**** SPECTRAL ANALYSIS OF DRY STATIC ENERGY FLUX, LATENT HEAT FLUX,
C**** AND ANGULAR MOMENTUM FLUX
      CALL FFT(AIJL2(1,J,L),FCPVA(0,J,L),FCPVB(0,J,L))
      CALL FFT(AIJL3(1,J,L),FCJKA(0,J,L),FCJKB(0,J,L))
      DO 342 N=1,NWAV_DAG
  342 AJLSP(J,L,N,1)=AJLSP(J,L,N,1)+.5*FIM*(FCPVA(N-1,J,L)*
     *  FCJKA(N-1,J,L)+FCPVB(N-1,J,L)*FCJKB(N-1,J,L))
      CALL FFT(AIJL4(1,J,L),FCJKA(0,J,L),FCJKB(0,J,L))
      DO 343 N=1,NWAV_DAG
  343 AJLSP(J,L,N,2)=AJLSP(J,L,N,2)+.5*FIM*(FCPVA(N-1,J,L)*
     *  FCJKA(N-1,J,L)+FCPVB(N-1,J,L)*FCJKB(N-1,J,L))
      CALL FFT(AIJL1(1,J,L),FCJKA(0,J,L),FCJKB(0,J,L))
      DO 344 N=1,NWAV_DAG
  344 AJLSP(J,L,N,3)=AJLSP(J,L,N,3)+.5*FIM*(FCPVA(N-1,J,L)*
     *  FCJKA(N-1,J,L)+FCPVB(N-1,J,L)*FCJKB(N-1,J,L))
  345 CONTINUE
  390 CONTINUE
C****
C**** VERTICAL MASS FLUXES  W(I,J,K)
C****
      DO 400 I=1,IM
      DO 400 J=1,JM
      DO 400 K=1,KM
      DUT(I,J,K)=0.
      DVT(I,J,K)=0.
  400 W(I,J,K)=0.
C**** EASTWARD MASS FLUX DUT (PU POINTS)
      DO 460 J=2,JM-1
      DO 460 K=1,KM
      I=IM
      DO 460 IP1=1,IM
      SP=.5*(P(I,J)+P(IP1,J))
      DO 405 L=1,LS1-1
  405 PL(L)=SP*SIGE(L)+PTOP
      IF (PM(K+1).GE.SP+PTOP) GO TO 460
      L=1
      PDN=SP+PTOP
      IF (PM(K).GE.SP+PTOP) GO TO 420
      PDN=PM(K)
  410 IF (PM(K).GT.PL(L+1)) GO TO 420
      L=L+1
      GO TO 410
  420 LUP=L
  430 IF (PM(K+1).GE.PL(LUP+1)) GO TO 440
      LUP=LUP+1
      GO TO 430
  440 CONTINUE
C**** CALCULATE HERE
  450 PUP=PL(L+1)
      IF (LUP.EQ.L) PUP=PM(K+1)
      DPDY=(PDN-PUP)*DYP(3)
      DUT(I,J,K)=DUT(I,J,K)+DPDY*(U(I,J,L)+U(I,J+1,L))
      IF (LUP.EQ.L) GO TO 460
      L=L+1
      PDN=PL(L)
      GO TO 450
  460 I=IP1
C**** NORTHWARD MASS FLUX DVT (PV POINTS)
      DO 520 J=2,JM
      DO 520 K=1,KM
      IM1=IM
      DO 520 I=1,IM
      SP=.5*(P(I,J-1)+P(I,J))
      DO 465 L=1,LS1-1
  465 PL(L)=SP*SIGE(L)+PTOP
      IF (PM(K+1).GE.SP+PTOP) GO TO 520
      L=1
      PDN=SP+PTOP
      IF (PM(K).GE.SP+PTOP) GO TO 480
      PDN=PM(K)
  470 IF (PM(K).GT.PL(L+1)) GO TO 480
      L=L+1
      GO TO 470
  480 LUP=L
  490 IF (PM(K+1).GE.PL(LUP+1)) GO TO 500
      LUP=LUP+1
      GO TO 490
  500 CONTINUE
C**** CALCULATE HERE
  510 PUP=PL(L+1)
      IF (LUP.EQ.L) PUP=PM(K+1)
      DPDX=(PDN-PUP)*DXV(J)
      DVT(I,J,K)=DVT(I,J,K)+DPDX*(V(IM1,J,L)+V(I,J,L))
      IF (LUP.EQ.L) GO TO 520
      L=L+1
      PDN=PL(L)
      GO TO 510
  520 IM1=I
C**** POLAR VERTICAL MASS FLUX
      DO 560 K=KM,1,-1
      W(1,1,K)=0.
      IF (K.LT.KM) W(1,1,K)=W(1,1,K+1)
      W(1,JM,K)=0.
      IF (K.LT.KM) W(1,JM,K)=W(1,JM,K+1)
  530 DO 540 I=1,IM
      W(1,1,K)=W(1,1,K)-.5*DVT(I,2,K)
  540 W(1,JM,K)=W(1,JM,K)+.5*DVT(I,JM,K)
C**** NON-POLAR VERTICAL MASS FLUX
      WUP=0.
      DO 560 J=2,JM-1
      IM1=IM
      DO 560 I=1,IM
      IF (K.LT.KM) WUP=W(I,J,K+1)
      W(I,J,K)=WUP+.5*(DUT(IM1,J,K)-DUT(I,J,K)+
     *  DVT(I,J,K)-DVT(I,J+1,K))
C**** ACCUMULATE ALL VERTICAL WINDS
  560 IM1=I
      INCHM=NDAA/NDYN
      IHOUR=1.5+TOFDAY
      DO 558 J=1,JM
      DO 558 I=1,IM
      DO KR=1,NDLYPT
         IF(I.EQ.IJD6(1,KR).AND.J.EQ.IJD6(2,KR)) THEN
            IH=IHOUR
            DO INCH=1,INCHM
               IF(IH.GT.HR_IN_DAY) IH=IH-HR_IN_DAY
               ADAILY(IH,60,KR)=ADAILY(IH,60,KR)+1.E5*W(I,J,3)/DXYP(J)
               IH=IH+1
            END DO
         END IF
      END DO
  558 CONTINUE
      DO 565 J=1,JM
      IMAX=IMAXJ(J)
      DO 565 K=1,KM
      WI=0.
      DO 562 I=1,IMAX
  562 WI=WI+W(I,J,K)
  565 AJK(J,K,25)=AJK(J,K,25)+WI
C**** ZERO OUT SUBSURFACE VERTICAL WINDS
      DO 568 J=1,JM
      DO 568 I=1,IM
      PS=P(I,J)+PTOP
      K=2
  566 IF (PM(K+1).LT.PS) GO TO 568
      W(I,J,K)=0.
      K=K+1
      GO TO 566
  568 CONTINUE
C****
C**** ACCUMULATE T,Z,Q VERTICAL TRANSPORTS
C****
      DO 610 J=1,JM
      IMAX=IMAXJ(J)
      DO 610 K=2,KM
      WI=0.
      TKI=0.
      QKI=0.
      ZKI=0.
      WTI=0.
      WQI=0.
      WZI=0.
         THKI=0.
         WTHI=0.
      FIMI=0.
      DO 600 I=1,IMAX
      SP=P(I,J)
      DO 569 L=1,LS1-1
  569 PLO(L)=SP*SIG(L)+PTOP
      IF (PM(K).GE.SP+PTOP) GO TO 600
      L=1
      IF (PM(K).GE.PLO(1)) GO TO 580
  570 LUP=L+1
      IF (L.EQ.LM) GO TO 580
      IF (PM(K).GE.PLO(LUP)) GO TO 575
      L=LUP
      GO TO 570
  575 DPUP=PM(K)-PLO(LUP)
      DPDN=PLO(L)-PM(K)
      BYDP=1./(DPDN+DPUP)
      TK=BYDP*(TX(I,J,L)*DPUP+TX(I,J,LUP)*DPDN)
      QK=Q(I,J,L)*Q(I,J,LUP)/(BYDP*(Q(I,J,L)*DPDN+Q(I,J,LUP)*DPUP)+
     *  ZERO20)
      ZK=BYDP*(PHI(I,J,L)*DPUP+PHI(I,J,LUP)*DPDN)
         THK=BYDP*(T(I,J,L)*DPUP+T(I,J,LUP)*DPDN)
      GO TO 590
C**** SPECIAL CASES;  L=1, L=LM
  580 TK=TX(I,J,L)
      QK=Q(I,J,L)
      ZK=PHI(I,J,L)
         THK=T(I,J,L)
C**** MERIDIONAL AVERAGING
  590 WI=WI+W(I,J,K)
      TKI=TKI+TK
      QKI=QKI+QK
      ZKI=ZKI+ZK
      WTI=WTI+W(I,J,K)*TK
      WQI=WQI+W(I,J,K)*QK
      WZI=WZI+W(I,J,K)*ZK
         THKI=THKI+THK
         WTHI=WTHI+W(I,J,K)*THK
      FIMI=FIMI+1.
  600 CONTINUE
      BYFIM=ZERO20
      IF (FIMI.GT.ZERO20) BYFIM=1./FIMI
      AJK(J,K-1,26)=AJK(J,K-1,26)+BYFIM*(SHA*TKI+ZKI)*WI
      AJK(J,K-1,27)=AJK(J,K-1,27)+SHA*WTI+WZI
      AJK(J,K-1,28)=AJK(J,K-1,28)+BYFIM*QKI*WI
      AJK(J,K-1,29)=AJK(J,K-1,29)+WQI
      AJK(J,K-1,30)=AJK(J,K-1,30)+WZI-BYFIM*WI*ZKI
C     AJK(J,K-1,31)=AJK(J,K-1,31)+WTI-BYFIM*WI*TKI
         WJK(J,K)=BYFIM*WI/DXYP(J)
         WTJK(J,K)=BYFIM*(WTHI-BYFIM*WI*THKI)/DXYP(J)
         AJK(J,K-1,50)=AJK(J,K-1,50)+WTJK(J,K)
  610 CONTINUE
C****
C**** BAROCLINIC EDDY KINETIC ENERGY GENERATION
C****
      DO 630 J=1,JM
      IMAX=IMAXJ(J)
      DO 630 K=1,KM
      FIMI=0.
      W2I=0.
      PAI=0.
      WPA2I=0.
      DO 626 I=1,IMAX
      SP=P(I,J)
      DO 611 L=1,LS1-1
  611 PL(L)=SP*SIGE(L)+PTOP
      PS=SP+PTOP
      IF (PM(K+1).GE.PS) GO TO 626
      L=1
      PDN=PS
      IF (PM(K).GE.PS) GO TO 614
      PDN=PM(K)
  612 IF (PM(K).GT.PL(L+1)) GO TO 614
      L=L+1
      GO TO 612
  614 LUP=L
  616 IF (PM(K+1).GE.PL(LUP+1)) GO TO 618
      LUP=LUP+1
      GO TO 616
  618 CONTINUE
      PTK=0.
C**** INTERPOLATE HERE
  620 PUP=PL(L+1)
      IF (LUP.EQ.L) PUP=PM(K+1)
      DP=PDN-PUP
      PTK=PTK+DP*TX(I,J,L)
      IF (LUP.EQ.L) GO TO 622
      L=L+1
      PDN=PL(L)
      GO TO 620
C**** ACCUMULATE HERE
  622 FIMI=FIMI+1.
      WUP=0.
      IF (K.LT.KM) WUP=W(I,J,K+1)
      W2I=W2I+W(I,J,K)+WUP
      PY=PMO(K)
      IF (PM(K).GE.PS) PY=.5*(PS+PM(K+1))
      PAK=PTK/PY
      PAI=PAI+PAK
      WPA2I=WPA2I+(W(I,J,K)+WUP)*PAK
  626 CONTINUE
  630 AJK(J,K,31)=AJK(J,K,31)-(WPA2I-W2I*PAI/(FIMI+ZERO20))
C****
C**** ACCUMULATE UV VERTICAL TRANSPORTS
C****
C**** DOUBLE POLAR WINDS
      DO 640 K=1,KM
      WSP=2.*W(1,1,K)/FIM
      WNP=2.*W(1,JM,K)/FIM
      DO 640 I=1,IM
      W(I,1,K)=WSP
  640 W(I,JM,K)=WNP
      DO 710 J=2,JM
      UEARTH=RADIUS*OMEGA*COSV(J)
      I=IM
      DO 650 IP1=1,IM
      PSEC(I)=.25*(P(I,J-1)+P(IP1,J-1)+P(I,J)+P(IP1,J))
  650 I=IP1
      DO 710 K=2,KM
      W4I=0.
      UKI=0.
      WU4I=0.
      WKE4I=0.
      FIMI=0.
      I=IM
      DO 700 IP1=1,IM
      SP=PSEC(I)
      DO 660 L=1,LS1-1
  660 PLO(L)=SP*SIG(L)+PTOP
      IF (PM(K).GE.SP+PTOP) GO TO 700
      L=1
      IF (PM(K).GE.PLO(1)) GO TO 680
  670 LUP=L+1
      IF (L.EQ.LM) GO TO 680
      IF (PM(K).GE.PLO(LUP)) GO TO 675
      L=LUP
      GO TO 670
  675 DPUP=PM(K)-PLO(LUP)
      DPDN=PLO(L)-PM(K)
      BYDP=1./(DPDN+DPUP)
      UK=BYDP*(U(I,J,L)*DPUP+U(I,J,LUP)*DPDN)
      VK=BYDP*(V(I,J,L)*DPUP+V(I,J,LUP)*DPDN)
      GO TO 690
C**** SPECIAL CASES;  L=1,L=LM
  680 UK=U(I,J,L)
      VK=V(I,J,L)
C**** MERIDIONAL AVERAGING
  690 W4=W(I,J-1,K)+W(IP1,J-1,K)+W(I,J,K)+W(IP1,J,K)
      W4I=W4I+W4
      UKI=UKI+UK
      WU4I=WU4I+W4*UK
      WKE4I=WKE4I+W4*(UK*UK+VK*VK)
      FIMI=FIMI+1.
  700 I=IP1
      BYFIM=1./(FIMI+ZERO20)
         WUJK(J,K)=.25*(WU4I-W4I*UKI*BYFIM)*BYFIM/DXYV(J)
      AJK(J,K-1,36)=AJK(J,K-1,36)+WKE4I
      AJK(J,K-1,37)=AJK(J,K-1,37)+WU4I-BYFIM*W4I*UKI
  710 AJK(J,K-1,38)=AJK(J,K-1,38)+WU4I+W4I*UEARTH
C****
C**** POTENTIAL VORTICITY AND VERTICAL TRANSPORT OF POT. VORT.
C****
      DO 760 J=2,JM-1
      JHEMI=1
      IF (J.LT.1+JM/2) JHEMI=-1
      DO 730 K=1,KM
      PVI=0.
      DO 720 I=1,IM
      DUT(I,J,K)=JHEMI*STB(I,J,K)*(ZX(I,J,K)-FCOR(J))
  720 PVI=PVI+DUT(I,J,K)
  730 AJK(J,K,32)=AJK(J,K,32)+PVI
      DO 760 K=2,KM
      W2I=0.
      PV2I=0.
      WPV4I=0.
      FIMI=0.
      I=IM
      DO 740 IP1=1,IM
      PS=.5*(P(I,J)+P(IP1,J))+PTOP
      IF (PM(K).GE.PS) GO TO 740
      W2=W(I,J,K)+W(IP1,J,K)
      W2I=W2I+W2
      PV2=DUT(I,J,K-1)+DUT(I,J,K)
      PV2I=PV2I+PV2
      WPV4I=WPV4I+W2*PV2
      FIMI=FIMI+1.
  740 I=IP1
      AJK(J,K-1,33)=AJK(J,K-1,33)+WPV4I
  760 AJK(J,K-1,34)=AJK(J,K-1,34)+WPV4I-W2I*PV2I/(FIMI+ZERO20)
C****
C**** SPECIAL MEAN/EDDY DIAGNOSTICS ARE CALCULATED
C****
      DO 770 J=2,JM
      DO 765 K=2,KM
      DPE=PMO(K)-PMO(K-1)
      UP(J,K)=(UJK(J,K)-UJK(J,K-1))/DPE
  765 PSIP(J,K)=(PSIJK(J,K)-PSIJK(J,K-1))/DPE
      UP(J,1)=UP(J,2)
      PSIP(J,1)=PSIP(J,2)
  770 CONTINUE
      DO 780 K=1,KM
      KUP=K+1
      IF (K.EQ.KM) KUP=KM
      KDN=K-1
      IF (K.EQ.1) KDN=1
      DO 780 J=2,JM
      TY(J,K)=(TJK(J,K)-TJK(J-1,K))/DYV(J)
C**** E-P FLUX NORTHWARD COMPONENT
      AJK(J,K,44)=AJK(J,K,44)+PSIJK(J,K)*(UJK(J,KUP)-UJK(J,KDN))/
     *  (PMO(KUP)-PMO(KDN))-UVJK(J,K)
  780 CONTINUE
      DO 800 J=2,JM-1
      DO 800 K=2,KMM1
      UY=(UJK(J+1,K)*DXV(J+1)-UJK(J,K)*DXV(J)-FCOR(J))/DXYP(J)
      PSIY=(PSIJK(J+1,K)*DXV(J+1)-PSIJK(J,K)*DXV(J))/DXYP(J)
C**** ZONAL MEAN MOMENTUM EQUATION   (MEAN ADVECTION)
      AJK(J,K,40)=AJK(J,K,40)-.5*UY*(VJK(J,K)+VJK(J+1,K))-
     *  .25*((UP(J+1,K+1)+UP(J,K+1))*WJK(J,K+1)+(UP(J+1,K)+UP(J,K))*
     *   WJK(J,K))
C**** ZONAL MEAN HEAT EQUATION   (MEAN ADVECTION)
      AJK(J,K,41)=AJK(J,K,41)-.5*(TY(J,K)*VJK(J,K)+TY(J+1,K)*VJK(J+1,K))
     *  -.5*STJK(J,K)*(WJK(J,K+1)+WJK(J,K))
C**** LAGRANGIAN MEAN MOMENTUM EQUATION  (MEAN ADVECTION)
      VSTAR=.5*(VJK(J,K)+VJK(J+1,K)-.5*(PSIP(J,K)+PSIP(J,K+1)
     *  +PSIP(J+1,K)+PSIP(J+1,K+1)))
      WSTAR=.5*(WJK(J,K)+WJK(J,K+1))+PSIY
      AJK(J,K,42)=AJK(J,K,42)-UY*VSTAR-.25*(UP(J,K)+UP(J+1,K)+
     *  UP(J,K+1)+UP(J+1,K+1))*WSTAR
      AJK(J,K,43)=AJK(J,K,43)-.5*(TY(J+1,K)+TY(J,K))*VSTAR-
     *  STJK(J,K)*WSTAR
C**** VERTICAL E-P FLUX
      AJK(J,K-1,45)=AJK(J,K-1,45)-WUJK(J,K)-.5*PSIJK(J,K)*UY
      AJK(J,K,45)=AJK(J,K,45)-.5*PSIJK(J,K)*UY
  800 CONTINUE
C****
C**** SPECTRAL ANALYSIS OF KINETIC ENERGIES AT CONSTANT PRESSURE
C****
      IZERO=0
      NM=1+IM/2
      J45N=2+.75*(JM-1.)
      KS1=LS1
C**** TOTAL THE KINETIC ENERGIES
      KE(:,:)=0.
      DO 2140 J=2,JM
      I=IM
      DO 2020 IP1=1,IM
      PSEC(I)=.25*(P(I,J-1)+P(IP1,J-1)+P(I,J)+P(IP1,J))
 2020 I=IP1
      DO 2140 K=1,KM
      KSPHER=2
      IF (K.GE.KS1) KSPHER=1
      IF (J.GT.JEQ) KSPHER=KSPHER+2
      DO 2140 KX=IZERO,LM,LM
      DO 2090 I=1,IM
      DPUV=0.
      SP=PSEC(I)
      DO 2025 L=1,LS1-1
      PLO(L)=SP*SIG(L)+PTOP
 2025 PL(L)=SP*SIGE(L)+PTOP
      PS=SP+PTOP
      IF (PM(K+1).GE.PLO(1)) GO TO 2090           ! really ?? not PL?
      L=1
      PDN=PS
      IF (PM(K).GE.PLO(1)) GO TO 2040             ! really ?? not PL?
      PDN=PM(K)
 2030 IF (PM(K).GT.PL(L+1)) GO TO 2040
      L=L+1
      GO TO 2030
 2040 LUP=L
 2050 IF (PM(K+1).GE.PL(LUP+1)) GO TO 2060
      LUP=LUP+1
      GO TO 2050
 2060 CONTINUE
C**** ACCUMULATE HERE
      SQRTDP=SQRT(PDN-PM(K+1))
 2070 PUP=PL(L+1)
      IF (LUP.EQ.L) PUP=PM(K+1)
      DP=PDN-PUP
      IF(KX.EQ.IZERO) DPUV=DPUV+DP*U(I,J,L)
      IF(KX.EQ.LM)    DPUV=DPUV+DP*V(I,J,L)
      IF (LUP.EQ.L) GO TO 2080
      L=L+1
      PDN=PL(L)
      GO TO 2070
 2080 IF (SQRTDP.EQ.0.) SQRTDP=ZERO20
      DPUV=DPUV/SQRTDP
 2090 X1(I)=DPUV
      CALL FFTE (X1,X1)
      IF (J.EQ.JEQ) GO TO 2120
      DO 2100 N=1,NM
 2100 KE(N,KSPHER)=KE(N,KSPHER)+X1(N)*DXYV(J)
      IF (J.NE.J45N) GO TO 2140
      DO 2110 N=1,NM
 2110 KE(N,KSPHER+4)=KE(N,KSPHER+4)+X1(N)*DXYV(J)
      GO TO 2140
 2120 DO 2130 N=1,NM
      KE(N,KSPHER+4)=KE(N,KSPHER+4)+X1(N)*DXYV(J)
      KE(N,KSPHER)=KE(N,KSPHER)+.5D0*X1(N)*DXYV(J)
 2130 KE(N,KSPHER+2)=KE(N,KSPHER+2)+.5D0*X1(N)*DXYV(J)
 2140 CONTINUE
      DO 2150 KS=1,NSPHER
      DO 2150 N=1,NM
 2150 SPECA(N,18,KS)=SPECA(N,18,KS)+KE(N,KS)
C**** ACCUMULATE TIME USED IN DIAGA
      MEND = MCLOCK ()
      MDIAG = MDIAG + (MEND-MBEGIN)
      MDYN  = MDYN  - (MEND-MBEGIN)
C     WRITE (6,997) IDACC,MINC,MDIAG
      RETURN
  997 FORMAT (' DIAGNOSTICS ACCUMULATED ',12I4,15X,2I7)
      END SUBROUTINE DIAGB

      SUBROUTINE DIAGJ
C****
C**** THIS SUBROUTINE PRODUCES AREA WEIGHTED STATISTICS OF
C****
C   K   N
C****
C***1   1  SOLAR RADIATION INCIDENT ON PLANET (W/M**2)
C****
C**1A   2/1  PLANETARY ALBEDO (10**-2)
C**1B  72/1  PLANETARY ALBEDO VISUAL (10**-2)
C**1C  73/1  PLANETARY ALBEDO NEAR IR (10**-2)
C**1D   6/5  GROUND ALBEDO (10**-2)
C**1E  74/1  GROUND ALBEDO VISUAL (10**-2)
C**1F  75/1  GROUND ALBEDO NEAR IR (10**-2)
C**1G  76/1  ATMOSPHERIC ALBEDO VISUAL (10**-2)
C**1H  77/1  ATMOSPHERIC ALBEDO NEAR IR (10**-2)
C**1I  78/1  ATMOSPHERIC ABSORPTION VISUAL (10**-2)
C**1J  79/1  ATMOSPHERIC ABSORPTION NEAR IR (10**-2)
C****
C***2   2  SOLAR RADIATION ABSORBED BY PLANET (W/M**2)
C***3   3  SOLAR RADIATION ABSORBED BELOW PTOP (W/M**2)
C***4   4  SOLAR RADIATION ABSORBED BY ATMOSPHERE (W/M**2)
C***5   5  SOLAR RADIATION INCIDENT ON GROUND (W/M**2)
C***6   6  SOLAR RADIATION ABSORBED BY GROUND (W/M**2)
C***7  32  SOLAR RADIATION WATER CORRECTION
C***8   7  THERMAL RADIATION EMITTED BY PLANET (W/M**2)
C***9   8  THERMAL RADIATION AT PTOP (W/M**2)
C**10   9  THERMAL RADIATION EMITTED BY GROUND (W/M**2)
C****
C**11  67  THERMAL RADIATION INCIDENT ON GROUND (W/M**2)
C****  55  BRIGHTNESS TEMPERATURE THROUGH WINDOW REGION (K-273.16)
C****  10  NET RADIATION ABSORBED BY PLANET (W/M**2)
C****  11  NET RADIATION ABSORBED BELOW PTOP (W/M**2)
C****  12  NET RADIATION ABSORBED BY GROUND (W/M**2)
C****  13  SENSIBLE HEAT FLUX INTO THE GROUND (W/M**2)
C****  14  EVAPORATION HEAT FLUX INTO THE GROUND (W/M**2)
C****  39  PRECIPITATION HEAT FLUX INTO THE GROUND (W/M**2)
C****  40  HEAT RUNOFF FROM FIRST GROUND LAYER (W/M**2)
C****  44  NET HEATING AT Z0 (W/M**2)
C****
C**21  42  CONDUCTION AT -Z1 (W/M**2)
C****  41  HEAT OF WATER OR ICE DUFFUSION AT -Z1 (W/M**2)
C****  16  NET HEATING AT -Z1 (W/M**2)
C****  15  CONDUCTION AT -Z1-Z2 (W/M**2)
C****  43  ENERGY OF ICE MELTING (OR TRANSPORTING) AT -Z1-Z2 (W/M**2)
C****  56  NET HEATING AT -Z1-Z2 (W/M**2)
C****  33  OCEAN TRANSPORT (W/M**2)
C****  48  HEAT RUNOFF THROUGH THE MIXED LAYER DEPTH (W/M**2)
C****  68  ENERGY DIFFUSION INTO THE THERMOCLINE (W/M**2)
C****  18  MEAN TEMPERATURE OF FIRST GROUND LAYER (.1 K-273.16)
C****
C**31  17  MEAN TEMPERATURE OF SECOND GROUND LAYER (.1 K-273.16)
C****  34  OCEAN TEMPERATURE AT THE MAXIMUM MIXED LAYER DEPTH
C****  23  SURFACE AIR TEMPERATURE (.1 K-273.16)
C****  22  FIRST LAYER AIR TEMPERATURE (.1 K-273.16)
C****  21  COMPOSITE AIR TEMPERATURE (.1 K-273.16)
C****  35  STRATO TEMPERATURE CHANGE PER DEGREE LATITUDE (10**-2 K)
C****  36  TROPO TEMPERATURE CHANGE PER DEGREE LATITUDE (10**-2 K)
C****  24  STRATOSPHERIC STATIC STABILITY (10**-3 K/M)
C****  25  TROPOSPHERIC STATIC STABILITY (10**-3 K/M)
C****  26  STRATOSPHERIC RICHARDSON NUMBER (1)
C****
C**41  27  TROPOSPHERIC RICHARDSON NUMBER (1)
C****  28  STRATOSPHERIC ROSSBY NUMBER (1)
C****  29  TROPOSPHERIC ROSSBY NUMBER (1)
C****  37  L IN THE STRATOSPHERE (10**5 M)
C****  38  L IN THE TROPOSPHERE (10**5 M)
C****  64  GAM  (10**-3 K/M)
C****  65  GAMM  (10**-3 K/M)
C****  66  GAMC  (10**-3 K/M)
C****  57  INTEGRATED SUPER-SATURATION CLOUD COVER (10**-2)
C****  58  INTEGRATED MOIST CONVECTIVE CLOUD COVER (10**-2)
C****
C**51  59  INTEGRATED TOTAL CLOUD COVER (10**-2)
C****  60  MOIST CONVECTIVE CLOUD DEPTH (100 N)
C****  61  SUPER SATURATION PRECIPITATION (KG/M**2/86400 S)
C****  62  MOIST CONVECTIVE PRECIPITATION (KG/M**2/86400 S)
C****  20  PRECIPITATION (KG/M**2/86400 S)
C****  19  EVAPORATION (KG/M**2/86400 S)
C****  63  WATER CONTENT OF ATMOSPHERE (KG/M**2)
C****  54  WATER RUNOFF AT Z0 (KG/M**2/86400 S)
C****  45  WATER OR ICE DIFFUSION AT -Z1 (KG/M**2/86400 S)
C****  46  ICE MELTING (OR TRANSPORTING) AT -Z1-Z2 (KG/M**2/86400 S)
C****
C**61  47  WATER RUNOFF THROUGH MIXED LAYER DEPTH (KG/M**2/86400 S)
C****  49  WATER CONTAINED IN FIRST GROUND LAYER (KG/M**2)
C****  50  ICE CONTAINED IN FIRST GROUND LAYER (KG/M**2)
C****  51  WATER CONTAINED IN SECOND GROUND LAYER (KG/M**2)
C****  52  ICE CONTAINED IN SECOND GROUND LAYER (KG/M**2)
C****  53  SNOW DEPTH (KG/M**2)
C****  31  SNOW COVER (10**-2)
C**68  30  OCEAN ICE COVER (10**-2)
C****
      USE CONSTANT, only : grav,rgas,sday,twopi,omega,kapa
      USE E001M12_COM,
     &     only : im,jm,lm,fim,
     &     DT,FLAND,IDACC,IDAY,IDAY0,JDATE,JDATE0,JMNTH0,JMONTH,
     &     JYEAR,JYEAR0,LS1,NCNDS,NDYN,SIGE,TAU,TAU0,TOFDAY,
     &     TOFDY0,XLABEL
      USE GEOM, only : DXYP,LAT,JLAT
      USE DAGCOM, only : aj,bj,cj,areg,jreg,kdiag,namreg,nreg,kaj,
     *     j_ctopp,j_cdldep,j_pcldmc,j_srabs,j_srnfp0,j_srnfg,j_trnfp0,
     *     j_hsurf,j_trhdt,j_trnfp1,j_hatm,j_rnfp0,j_rnfp1,j_srnfp1,
     *     j_rhdt,j_hz1,j_edifs,j_f1dt,j_prcp,j_prcpss,j_prcpmc,j_hz0,
     *     j_shdt,j_evhdt,j_eprcp,j_erun1,j_hz2,j_f2dt,j_erun2,j_rsi

      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(JM) ::
     &     CONTJ,CONTO,CONTL,CONTOI,S1,SPOCEN,SPOICE,SPLAND,FLAT
      INTEGER, DIMENSION(JM) :: MLAT
      DOUBLE PRECISION, DIMENSION(NREG) :: SAREA
      INTEGER, DIMENSION(2) :: MHEM
      DOUBLE PRECISION, DIMENSION(2) :: FHEM
      INTEGER III
C****
      CHARACTER*16 :: TERAIN(5) = (
     *     /'    (GLOBAL)','      (LAND)','     (OCEAN)',
     *      ' (OCEAN ICE)','   (REGIONS)'/)
      CHARACTER*16 TITLEA(10)
      CHARACTER(LEN=16), DIMENSION(KAJ) :: TITLE = (/
     1  ' INC SW(WT/M**2)', '0SW ABS BELOW P0', ' SW ABS BELOW P1',
     4  ' SW ABS BY ATMOS', ' SW INC ON Z0   ', ' SW ABS AT Z0   ',
     7  '0NET LW AT P0   ', ' NET LW AT P1   ', ' NET LW AT Z0   ',
     O  '0NET RAD AT P0  ', ' NET RAD AT P1  ', ' NET RAD AT Z0  ',
     3  '0SENSBL HEAT FLX', ' EVAPOR HEAT FLX', '0CONDC AT -Z1-Z2',
     6  ' NET HEAT AT -Z1', ' TG2 (.1 C)     ', '1TG1 (.1 C)     ',
     9  ' EVAPOR (MM/DAY)', ' PRECIP (MM/DAY)', ' T AIR (.1 C)   ',
     2  ' T1 (.1 C)      ', '0T SURF (.1 C)  ', '0STAT STB(STRAT)',
     5  ' STAT STB(TROPO)', '0RICH NUM(STRAT)', ' RICH NUM(TROPO)',
     8  ' ROSS NUM(STRAT)', ' ROSS NUM(TROPO)', ' OCEAN ICE COVER',
     1  '0SNOW COVER     ', ' SW CORRECTION  ', '0OCEAN TRANSPORT',
     4  ' TG3 (.1 C)     ', '0DT/DLAT(STRAT) ', ' DT/DLAT(TROPO) ',
     7  ' L(STRAT)(10**5)', ' L(TROP) (10**5)', ' PRECIP HEAT FLX',
     O  ' HEAT RUNOFF Z0 ', ' HT WTR DIFS -Z1', '0CONDUCTN AT -Z1',
     3  ' ICE ENRG -Z1-Z2', ' NET HEAT AT Z0 ', ' H2O DIFS AT -Z1',
     6  ' ICE THRU -Z1-Z2', ' WATR RUNOFF MLD', ' HEAT RUNOFF MLD',
     9  '0WATER IN G1    ', ' ICE IN G1      ', ' WATER IN G2    ',
     2  ' ICE IN G2      ', ' SNOW DEPTH     ', '0WATER RUNOFF Z0',
     5  ' LW WINDOW BTEMP', ' NET HEAT -Z1-Z2', '0TOT SUP SAT CLD',
     8  ' TOT MST CNV CLD', ' TOTAL CLD COVER', ' MC CLD DPTH(MB)',
     *  '0SS PRECIP(MM/D)', ' MC PRECIP(MM/D)', ' H2O OF ATM (MM)',
     4  '0GAM(K/KM)      ', ' GAMM(K/KM)     ', ' GAMC(K/KM)     ',
     *  ' LW INC ON Z0   ', ' HT INTO THRMOCL',
     & ('                ',iii=1,26)/)
      DATA TITLEA/' PLANETARY ALBDO',' PLAN ALB VISUAL',
     *  ' PLAN ALB NEARIR', ' SURFACE G ALBDO', ' SURF ALB VISUAL',
     *  ' SURF ALB NEARIR', '0ATMO ALB VISUAL', ' ATMO ALB NEARIR',
     *  ' ATMO ABS VISUAL', ' ATMO ABS NEARIR'/
C****
      DOUBLE PRECISION, DIMENSION(4) ::
     &     WTA=(/1.,0.,1.,0./),WTB=(/1.,1.,0.,0./),WTC=(/1.,0.,0.,1./)
      INTEGER, DIMENSION(70) :: NDEX
      DATA NDEX/1,2,3,4,5,6,32,7,8,9,  67,55,10,11,12,13,14,39,40,44,
     *  42,41,16,15,43,56,33,48,68,18,  17,34,23,22,21,35,36,24,25,26,
     *  27,28,29,37,38,64,65,66,57,58,  59,60,61,62,20,19,63,54,45,46,
     *  47,49,50,51,52,53,31,30,2*0/
      INTEGER, DIMENSION(10) :: INNUM,INDEN
      DATA INNUM/2,72,73,6,74,75,76,77,78,79/, INDEN/3*1,5,6*1/
C**** IA: 1 CONDENSATION, 2 RADIATION, 3 SURFACE, 4 DIAGA, 0 UNUSED
      INTEGER, DIMENSION(KAJ) :: IA
      DATA IA/6*2,  2,2,1,2,2,1,  6*1,  1,1,4,4,3,4,  5*4,1,
     *  4,2,1,1,4,4, 4,4,2*1,1,1, 1,1,3*1,1, 6*1, 2,1,4*2,  1,1,4*4,
     *  2,9,4*0,8*0,14*1/
      DOUBLE PRECISION, DIMENSION(KAJ) :: SCALE
      DATA SCALE/6*1.,  6*1.,  4*1.,2*10.,  2*1.,4*10.,  6*100.,
     *  100.,2*1.,10.,2*100.,  6*1.,  6*1.,  6*1.,  2*1.,3*100.,1.,
     *  6*1., 6*1., 22*1./

      DOUBLE PRECISION ::
     &     A1BYA2,A2BYA1,AMULT,BYA1,BYIACC,DTCNDS,DTSRCE,
     &     FGLOB,GSUM,GSUM2,GWT,HSUM,HSUM2,HWT,QDEN,QJ,QNUM,
     &     TAUDIF,WTX

      INTEGER ::
     &     I,IACC,IHOUR,IHOUR0,INC,J,JH,JHEMI,JMHALF,JR,
     &     K,KA,KD1M,M,M1,N,ND,NN,NSTPDJ

      INTEGER :: IFIRST
      DATA IFIRST/1/
      DOUBLE PRECISION, PARAMETER :: P1000=1000.
      IF (IFIRST.NE.1) GO TO 100
      IFIRST=0
C**** INITIALIZE CERTAIN QUANTITIES  (KD1M LE 69)
      KD1M=68
      INC=1+(JM-1)/24
      DTSRCE=DT*NDYN
      DTCNDS=DT*NCNDS
      JMHALF=JM/2
      DO 10 JR=1,NREG
   10 SAREA(JR)=0.
      DO 30 J=1,JM
      S1(J)=IM
      SPLAND(J)=0.
      DO 20 I=1,IM
      SPLAND(J)=SPLAND(J)+FLAND(I,J)
      JR=JREG(I,J)
   20 SAREA(JR)=SAREA(JR)+DXYP(J)
   30 CONTINUE
      S1(1)=1.
      S1(JM)=1.
      SPLAND(1)=FLAND(1,1)
      SPLAND(JM)=FLAND(1,JM)
      SCALE(9)=1./DTSRCE
      SCALE(12)=1./DTSRCE
      SCALE(13)=1./DTSRCE
      SCALE(14)=1./DTSRCE
      SCALE(15)=1./DTSRCE
      SCALE(16)=1./DTSRCE
      SCALE(19)=SDAY/DTSRCE
      SCALE(20)=100.*SDAY/(DTCNDS*GRAV)
      SCALE(24)=1.D3*GRAV*(P1000**KAPA)
      SCALE(25)=SCALE(24)
      SCALE(26)=16.*RGAS
      SCALE(27)=16.*RGAS
      SCALE(28)=.5/(2.*OMEGA*FIM)
      SCALE(29)=.5/(2.*OMEGA*FIM)
      SCALE(33)=1./DTSRCE
      SCALE(35)=.5D2*(JM-1.)/((SIGE(LS1)-SIGE(LM+1)+1.D-12)*180.)
      SCALE(36)=.5E2*(JM-1.)/((SIGE(1)-SIGE(LS1))*180.)
      SCALE(37)=1.D-5*SQRT(RGAS)/(2.*OMEGA)
      SCALE(38)=SCALE(37)
      SCALE(39)=1./DTSRCE
      SCALE(40)=1./DTSRCE
      SCALE(41)=1./DTSRCE
      SCALE(42)=1./DTSRCE
      SCALE(43)=1./DTSRCE
      SCALE(44)=1./DTSRCE
      SCALE(45)=SDAY/DTSRCE
      SCALE(46)=SDAY/DTSRCE
      SCALE(47)=SDAY/DTSRCE
      SCALE(48)=1./DTSRCE
      SCALE(54)=SDAY/DTSRCE
      SCALE(56)=1./DTSRCE
      SCALE(61)=SCALE(20)
      SCALE(62)=SCALE(20)
      SCALE(63)=100./GRAV
      SCALE(64)=1.D3*GRAV
      SCALE(65)=1.D3*.0098/(SIGE(1)-SIGE(LS1))
      SCALE(66)=1.D3
      SCALE(68)=2.E3*4185./SDAY
C**** CALCULATE THE DERIVED QUANTITIES
  100 BYA1=1./(IDACC(1)+1.D-20)
      A2BYA1=DFLOAT(IDACC(2))/DFLOAT(IDACC(1))
      A1BYA2=IDACC(1)/(IDACC(2)+1.D-20)
      DO 200 JR=1,23 ! only 23 will fit on a green sheet
      AREG(JR,J_SRABS)=AREG(JR,J_SRNFP0)-AREG(JR,J_SRNFG)
      AREG(JR,J_TRNFP0)=AREG(JR,J_HSURF)+A2BYA1*AREG(JR,J_TRHDT)/DTSRCE
      AREG(JR,J_TRNFP1)=AREG(JR,J_HATM )+A2BYA1*AREG(JR,J_TRHDT)/DTSRCE
      AREG(JR,J_RNFP0)=AREG(JR,J_SRNFP0)+AREG(JR,J_TRNFP0)
      AREG(JR,J_RNFP1)=AREG(JR,J_SRNFP1)+AREG(JR,J_TRNFP1)
      AREG(JR,J_RHDT)=A1BYA2*AREG(JR,J_SRNFG)*DTSRCE+AREG(JR,J_TRHDT)
      AREG(JR,J_HZ1)=AREG(JR,J_EDIFS)+AREG(JR,J_F1DT)
      AREG(JR,J_PRCP)=AREG(JR,J_PRCPSS)+AREG(JR,J_PRCPMC)
      AREG(JR,J_HZ0)=AREG(JR,J_RHDT)+AREG(JR,J_SHDT)+AREG(JR,J_EVHDT)+
     *  AREG(JR,J_EPRCP)-AREG(JR,J_ERUN1)
      AREG(JR,J_HZ2)=AREG(JR,J_F2DT)+AREG(JR,J_ERUN2)
  200 AREG(JR,J_CTOPP)=IDACC(2)*SAREA(JR)*AREG(JR,J_CDLDEP)/
     *     (AREG(JR,J_PCLDMC)+1.D-20)
      DO 210 J=1,JM
      SPOICE(J)=CJ(J,J_RSI)*BYA1
      SPOCEN(J)=S1(J)-SPLAND(J)-SPOICE(J)
c      AJ(J,J_SWCOR)=(1.-SRCOR)*AJ(J,J_SRNFG)   !obsolete
c      CJ(J,J_SWCOR)=(1.-SRCOR)*CJ(J,J_SRNFG)   !obsolete
      AJ(J,J_CTOPP)=IDACC(2)*SPOCEN(J)*AJ(J,J_CDLDEP)/(AJ(J,J_PCLDMC)
     *     +1.D-20)
      BJ(J,J_CTOPP)=IDACC(2)*SPLAND(J)*BJ(J,J_CDLDEP)/(BJ(J,J_PCLDMC)
     *     +1.D-20)
      CJ(J,J_CTOPP)=IDACC(2)*SPOICE(J)*CJ(J,J_CDLDEP)/(CJ(J,J_PCLDMC)
     *     +1.D-20)
      AJ(J,J_SRABS)=AJ(J,J_SRNFP0)-AJ(J,J_SRNFG)
      BJ(J,J_SRABS)=BJ(J,J_SRNFP0)-BJ(J,J_SRNFG)
      CJ(J,J_SRABS)=CJ(J,J_SRNFP0)-CJ(J,J_SRNFG)
      AJ(J,J_TRNFP0)=AJ(J,J_HSURF)+A2BYA1*AJ(J,J_TRHDT)/DTSRCE
      BJ(J,J_TRNFP0)=BJ(J,J_HSURF)+A2BYA1*BJ(J,J_TRHDT)/DTSRCE
      CJ(J,J_TRNFP0)=CJ(J,J_HSURF)+A2BYA1*CJ(J,J_TRHDT)/DTSRCE
      AJ(J,J_TRNFP1)=AJ(J,J_HATM)+A2BYA1*AJ(J,J_TRHDT)/DTSRCE
      BJ(J,J_TRNFP1)=BJ(J,J_HATM)+A2BYA1*BJ(J,J_TRHDT)/DTSRCE
      CJ(J,J_TRNFP1)=CJ(J,J_HATM)+A2BYA1*CJ(J,J_TRHDT)/DTSRCE
      AJ(J,J_RNFP0)=AJ(J,J_SRNFP0)+AJ(J,J_TRNFP0)
      BJ(J,J_RNFP0)=BJ(J,J_SRNFP0)+BJ(J,J_TRNFP0)
      CJ(J,J_RNFP0)=CJ(J,J_SRNFP0)+CJ(J,J_TRNFP0)
      AJ(J,J_RNFP1)=AJ(J,J_SRNFP1)+AJ(J,J_TRNFP1)
      BJ(J,J_RNFP1)=BJ(J,J_SRNFP1)+BJ(J,J_TRNFP1)
      CJ(J,J_RNFP1)=CJ(J,J_SRNFP1)+CJ(J,J_TRNFP1)
      AJ(J,J_RHDT)=A1BYA2*AJ(J,J_SRNFG)*DTSRCE+AJ(J,J_TRHDT)
      BJ(J,J_RHDT)=A1BYA2*BJ(J,J_SRNFG)*DTSRCE+BJ(J,J_TRHDT)
      CJ(J,J_RHDT)=A1BYA2*CJ(J,J_SRNFG)*DTSRCE+CJ(J,J_TRHDT)
      AJ(J,J_HZ1)=AJ(J,J_EDIFS)+AJ(J,J_F1DT)
      BJ(J,J_HZ1)=BJ(J,J_EDIFS)+BJ(J,J_F1DT)
      CJ(J,J_HZ1)=CJ(J,J_EDIFS)+CJ(J,J_F1DT)
      AJ(J,J_PRCP)=AJ(J,J_PRCPSS)+AJ(J,J_PRCPMC)
      BJ(J,J_PRCP)=BJ(J,J_PRCPSS)+BJ(J,J_PRCPMC)
      CJ(J,J_PRCP)=CJ(J,J_PRCPSS)+CJ(J,J_PRCPMC)
      AJ(J,J_HZ0)=AJ(J,J_RHDT)+AJ(J,J_SHDT)+AJ(J,J_EVHDT)+AJ(J,J_EPRCP)
     *           -AJ(J,J_ERUN1)
      BJ(J,J_HZ0)=BJ(J,J_RHDT)+BJ(J,J_SHDT)+BJ(J,J_EVHDT)+BJ(J,J_EPRCP)
     *           -BJ(J,J_ERUN1)
      CJ(J,J_HZ0)=CJ(J,J_RHDT)+CJ(J,J_SHDT)+CJ(J,J_EVHDT)+CJ(J,J_EPRCP)
     *           -CJ(J,J_ERUN1)
      AJ(J,J_HZ2)=AJ(J,J_F2DT)+AJ(J,J_ERUN2)
      BJ(J,J_HZ2)=BJ(J,J_F2DT)+BJ(J,J_ERUN2)
      CJ(J,J_HZ2)=CJ(J,J_F2DT)+CJ(J,J_ERUN2)
  210 CONTINUE
      IHOUR0=TOFDY0+.5
      IHOUR=TOFDAY+.5
      TAUDIF=TAU-TAU0
C****
C**** LOOP OVER SURFACE TYPES: GLOBAL, LAND, OCEAN, AND OCEAN ICE
C****
      IF (KDIAG(1).GT.7) GO TO 510
      M1=1
      NSTPDJ=1
      IF (KDIAG(1).GT.4) M1=KDIAG(1)-3
      IF (KDIAG(1).GT.3) NSTPDJ=4
      IF (KDIAG(1).EQ.3) NSTPDJ=2
      DO 500 M=M1,4,NSTPDJ
      IF (KDIAG(1).EQ.2.AND.M.GT.2) RETURN
      WRITE (6,901) XLABEL
      WRITE (6,902) TERAIN(M),IDAY0,IHOUR0,JDATE0,JMNTH0,JYEAR0,TAU0,
     *  IDAY,IHOUR,JDATE,JMONTH,JYEAR,TAU,TAUDIF
      WRITE (6,903) (JLAT(J,1),J=JM,INC,-INC)
      WRITE (6,905)
      DO 490 K=1,KD1M
      N=NDEX(K)
      IACC=IDACC(IA(N))
      GSUM=0.
      GWT=0.
      DO 320 JHEMI=1,2
      HSUM=0.
      HWT=0.
      DO 310 JH=1,JMHALF
      J=(JHEMI-1)*JMHALF+JH
      QJ=(AJ(J,N)*WTA(M)+BJ(J,N)*WTB(M)+CJ(J,N)*WTC(M))*SCALE(N)
      WTX=(SPOCEN (J)*WTA(M)+SPLAND(J)*WTB(M)+SPOICE(J)*WTC(M))*IACC
      FLAT(J)=QJ/(WTX+1.D-20)
      MLAT(J)=NINT(FLAT(J))
      HSUM=HSUM+QJ*DXYP(J)*(FIM+1.-S1(J))
  310 HWT=HWT+WTX*DXYP(J)*(FIM+1.-S1(J))
      FHEM(JHEMI)=HSUM/(HWT+1.D-20)
      GSUM=GSUM+HSUM
  320 GWT=GWT+HWT
      FGLOB=GSUM/(GWT+1.D-20)
         IF (M.EQ.1) CALL KEYDJ (N,FGLOB,FHEM(2))
CF       DO 323 J=1,JM
CF323    GBUDG(J,K,M)=FLAT(J)
CF       GBUDG(JM+1,K,M)=FHEM(1)
CF       GBUDG(JM+2,K,M)=FHEM(2)
CF       GBUDG(JM+3,K,M)=FGLOB
      GO TO (350,350,350,350,350,350,  350,350,350,350,350,350,
     *       350,350,348,348,350,350,  345,345,350,350,350,350,
     *       340,350,350,345,345,350,  350,340,348,350,350,350,
     *       345,345,348,345,348,348,  348,348,345,345,345,345,
     *       350,350,350,350,350,345,  350,348,350,350,350,350,
     *       345,345,350,345,345,345,  350,345,350,350,350,350),N
  340 WRITE (6,906) TITLE(N),FGLOB,FHEM(2),FHEM(1),
     *  (FLAT(J),J=JM,INC,-INC)
      GO TO 490
  345 WRITE (6,911) TITLE(N),FGLOB,FHEM(2),FHEM(1),
     *  (FLAT(J),J=JM,INC,-INC)
      GO TO 490
  348 WRITE (6,912) TITLE(N),FGLOB,FHEM(2),FHEM(1),
     *  (MLAT(J),J=JM,INC,-INC)
      GO TO 490
  350 WRITE (6,907) TITLE(N),FGLOB,FHEM(2),FHEM(1),
     *  (MLAT(J),J=JM,INC,-INC)
      IF (N.NE.1) GO TO 490
C**** CALCULATE AND PRINT ALBEDOS
  400 DO 430 KA=1,10
      NN=INNUM(KA)
      ND=INDEN(KA)
      AMULT=1.
      IF (KA.LE.1.OR.KA.EQ.4) AMULT=-1.
      GSUM=0.
      GSUM2=0.
      DO 420 JHEMI=1,2
      HSUM=0.
      HSUM2=0.
      DO 410 JH=1,JMHALF
      J=(JHEMI-1)*JMHALF+JH
      QNUM=AJ(J,NN)*WTA(M)+BJ(J,NN)*WTB(M)+CJ(J,NN)*WTC(M)
      QDEN=AJ(J,ND)*WTA(M)+BJ(J,ND)*WTB(M)+CJ(J,ND)*WTC(M)
      FLAT(J)=AMULT*(100.*     QNUM/(QDEN    +1.D-20)-50.)+50.
      MLAT(J)=FLAT(J)+.5
      HSUM=HSUM+QNUM*DXYP(J)*(FIM+1.-S1(J))
  410 HSUM2=HSUM2+QDEN*DXYP(J)*(FIM+1.-S1(J))
      FHEM(JHEMI)=50.+AMULT*(100.*HSUM/(HSUM2+1.D-20)-50.)
      GSUM=GSUM+HSUM
  420 GSUM2=GSUM2+HSUM2
      FGLOB=50.+AMULT*(100.*GSUM/(GSUM2+1.D-20)-50.)
         IF (M.EQ.1.AND.KA.EQ.1) CALL KEYDJA (FGLOB)
CF       DO 423 J=1,JM
CF423    GBUDG(J,KA+KD1M,M)=FLAT(J)
CF       GBUDG(JM+1,KA+KD1M,M)=FHEM(1)
CF       GBUDG(JM+2,KA+KD1M,M)=FHEM(2)
CF       GBUDG(JM+3,KA+KD1M,M)=FGLOB
      WRITE (6,912) TITLEA(KA),FGLOB,FHEM(2),FHEM(1),
     *  (MLAT(J),J=JM,INC,-INC)
  430 CONTINUE
  490 CONTINUE
      WRITE (6,903) (JLAT(J,1),J=JM,INC,-INC)
      WRITE (6,905)
      IF (KDIAG(1).GT.3) RETURN
  500 CONTINUE
  510 IF (KDIAG(1).GT.0.AND.KDIAG(1).NE.8) RETURN
C****
C**** PRODUCE REGIONAL STATISTICS
C****
      WRITE (6,901) XLABEL
      WRITE (6,902) TERAIN(5),IDAY0,IHOUR0,JDATE0,JMNTH0,JYEAR0,TAU0,
     *  IDAY,IHOUR,JDATE,JMONTH,JYEAR,TAU,TAUDIF
      WRITE(6,918)(NAMREG(1,K),K=1,23),(NAMREG(2,K),K=1,23)
      DO 700 K=1,KD1M
      N=NDEX(K)
      BYIACC=1./(IDACC(IA(N))+1.D-20)
      DO 520 JR=1,23
      FLAT(JR)=AREG(JR,N)*SCALE(N)*BYIACC/SAREA(JR)
  520 MLAT(JR)=NINT(FLAT(JR))
CF       DO 523 J=1,23
CF523    RBUDG(J,K)=FLAT(J)
      GO TO (550,550,550,550,550,550,  550,550,550,550,550,550,
     *       550,550,550,550,550,550,  540,540,550,550,550,550,
     *       540,550,550,540,540,540,  540,540,550,550,550,550,
     *       540,540,550,540,550,550,  550,550,540,540,540,540,
     *       540,540,550,550,540,540,  550,550,550,550,550,550,
     *       540,540,540,540,540,540,  550,540,550,550,550,550),N
  540 WRITE (6,910) TITLE(N),(FLAT(JR),JR=1,23)
      GO TO 700
  550 WRITE (6,909) TITLE(N),(MLAT(JR),JR=1,23)
      IF (N.NE.1) GO TO 700
      GO TO 600
C**** CALCULATE AND PRINT ALBEDOS FOR REGIONAL STATISTICS
  600 DO 620 KA=1,10
      NN=INNUM(KA)
      ND=INDEN(KA)
      AMULT=1.
      IF (KA.LE.1.OR.KA.EQ.4) AMULT=-1.
      DO 610 JR=1,23
      FLAT(JR)=AMULT*(100.*AREG(JR,NN)/(AREG(JR,ND)+1.D-20)-50.)+50.
  610 MLAT(JR)=FLAT(JR)+.5
CF       DO 613 J=1,23
CF613    RBUDG(J,KA+KD1M)=FLAT(J)
  620 WRITE (6,909) TITLEA(KA),(MLAT(JR),JR=1,23)
  700 CONTINUE
      WRITE (6,905)
      WRITE(6,918)(NAMREG(1,K),K=1,23),(NAMREG(2,K),K=1,23)
      RETURN
C****
  901 FORMAT ('1',A)
  902 FORMAT ('0** BUDGETS',A16,' **   DAY',I6,', HR',I2,' (',I2,A5,
     *  I4,')',F9.0,'   TO   DAY',I6,', HR',I2,' (',I2,A5,I4,')',
     *  F9.0,'   DIF',F5.0,' HR')
  903 FORMAT ('0',131('-')/20X,'G     NH    SH  ',24I4)
  904 FORMAT (A16,3I6,2X,24I4)
  905 FORMAT (1X,131('-'))
  906 FORMAT (A16,3F6.1,2X,24F4.1)
  907 FORMAT (A16,3F6.1,2X,24I4)
  909 FORMAT (A16,1X,23I5)
  910 FORMAT (A16,1X,23F5.1)
  911 FORMAT (A16,3F6.2,2X,24F4.1)
  912 FORMAT (A16,3F6.2,2X,24I4)
  918 FORMAT ('0',16X,23(1X,A4)/17X,23(1X,A4)/1X,131('-'))
      END SUBROUTINE DIAGJ

      BLOCK DATA BDJK
C****
C**** TITLES FOR SUBROUTINE DIAGJK
C****
      COMMON/DJKTTL/TITLE1,TITLE2,TITLE3,TITLE4,TITLE5
      CHARACTER*64 TITLE1(12)
      DATA TITLE1/
C****                                                              1-12
     1'TEMPERATURE (DEGREES CENTIGRADE)                               ',
     *'HEIGHT (HUNDREDS OF METERS)                                    ',
     3'SPECIFIC HUMIDITY (10**-5 KG H2O/KG AIR)                       ',
     *'RELATIVE HUMIDITY (PERCENT)                                    ',
     *'ZONAL WIND (U COMPONENT) (TENTHS OF METERS/SECOND)             ',
     6'MERIDIONAL WIND (V COMPONENT) (HUNDREDTHS OF METERS/SECOND)    ',
     *'                                                               ',
     *'                                                               ',
     9'BAROCLINIC EDDY KINETIC ENERGY GEN. (10**-1 WATTS/M**2/SIGMA)  ',
     A'NORTH. TRANS. OF EDDY Q-G POT. VORTICITY  (10**-6 M/S**2)      ',
     1'P-K BY EDDY PRESSURE GRADIENT FORCE  (10**-1 W/M**2/UNIT SIGMA ',
     2'DYNAMIC CONVERGENCE OF EDDY GEOPOTENTIAL (.1 WATTS/M**2/DSIGMA)'/
      CHARACTER*64 TITLE2(8)
      DATA TITLE2/
C****                                                             13-20
     3'NORTH. TRANS. OF SENSIBLE HEAT BY EDDIES (10**14 WATTS/DSIGMA) ',
     4'DYNAMIC CONVERGENCE OF DRY STATIC ENERGY (10 WATTS/M**2/DSIGMA)',
     *'                                                               ',
     *'                                                               ',
     7'STANDING EDDY KINETIC ENERGY (10**4 JOULES/M**2/UNIT SIGMA)    ',
     8'EDDY KINETIC ENERGY (10**4 JOULES/M**2/UNIT SIGMA)             ',
     9'TOTAL KINETIC ENERGY (10**4 JOULES/M**2/UNIT SIGMA)            ',
     *'                                                               '/
      CHARACTER*64 TITLE3(8)
      DATA TITLE3/
C****                                                             21-28
     1'POTENTIAL TEMPERATURE (DEGREES KELVIN)                          '
     *,
     2'NOR. TRANS. OF DRY STAT. ENERGY BY STAND. EDDIES (10**14 W/DSIG)'
     *,
     *'NORTH. TRANS. OF DRY STATIC ENERGY BY EDDIES (10**14 WATTS/DSIG)'
     *,
     4'TOTAL NORTH. TRANSPORT OF DRY STATIC ENERGY (10**15 WATTS/DSIG) '
     *,
     5'NORTHWARD TRANSPORT OF LATENT HEAT BY EDDIES (10**13 WATTS/DSIG)'
     *,
     6'TOTAL NORTHWARD TRANSPORT OF LATENT HEAT (10**14 WATTS/UNIT SIG)'
     *,
     7'NORTH.TRANSPORT OF STATIC ENERGY BY EDDIES (10**14 WATTS/DSIGMA)'
     *,
     8'TOTAL NORTHWARD TRANSPORT OF STATIC ENERGY (10**15 WATTS/DSIGMA)'
     */
      CHARACTER*64 TITLE4(8)
      DATA TITLE4/
C****                                                             29-36
     9'                                                               ',
     A'TOTAL NORTHWARD TRANSPORT OF KINETIC ENERGY (10**12 WATTS/DSIG)',
     1'NORTH. TRANS. OF ANG. MOMENTUM BY STAND. EDDIES (10**18 J/DSIG)',
     2'NORTH. TRANS. OF ANG. MOMENTUM BY EDDIES (10**18 JOULES/DSIGMA)',
     3'TOTAL NORTHWARD TRANSPORT OF ANG. MOMENTUM (10**19 JOULES/DSIG)',
     4'NORTHWARD WAVE ENERGY FLUX  (10**11 JOULES/METER/UNIT SIGMA)   ',
     5'VERTICAL WAVE ENERGY FLUX  (10**11 JOULES/METER)               ',
     6'DIVERGENCE OF THE WAVE ENERGY FLUX  (10**-5 M/S**2)            '/
      CHARACTER*64 TITLE5(8)
      DATA TITLE5/
C****                                                             37-44
     7'REFRACTION INDEX FOR WAVE NUMBER 1  (10**-8 PER METER**2)      ',
     8'REFRACTION INDEX FOR WAVE NUMBER 2  (10**-8 PER METER**2)      ',
     9'REFRACTION INDEX FOR WAVE NUMBER 3  (10**-8 PER METER**2)      ',
     O'REFRACTION INDEX FOR WAVE NUMBER 6  (10**-8 PER METER**2)      ',
     1'REFRACTION INDEX FOR WAVE NUMBER 9  (10**-8 PER METER**2)      ',
     2'Q-G POT. VORTICITY CHANGE OVER LATITUDES (10**-12 1/(SEC-M))   ',
     *'TOTAL CLOUD WATER CONTENT (10**-6 KG/KG)                       ',
     4'N. TRANSPORT OF LATENT HEAT BY STAND. EDDIES(10**13 WATTS/DSIG)'/
      END BLOCK DATA BDJK

      SUBROUTINE DIAGJK
      USE CONSTANT, only : grav,rgas,kapa,sday,lhe,twopi,omega,sha
      USE E001M12_COM, only : im,jm,lm,fim,
     &     BYIM,DSIG,DT,IDACC,IMH,LS1,NCNDS,NDAA,
     &     PSF,PTOP,PSFMPT,SIG,SIGE,SKIPSE,TOFDAY,TOFDY0
      USE GEOM, only :  
     &     AREAG,BYDXYP,COSP,COSV,DLON,DXV,DXYP,DXYV,DYP,FCOR,RADIUS,WTJ
      USE DAGPCOM
      USE RADNCB, only : LM_REQ
      USE DAGCOM, only : ajk,ajl,asjl,ajlsp,kdiag,aijl,aijk,
     &     nwav_dag,kajlsp
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(IM,JM) :: SENDEG
      common/ntdse/ sendeg
      DOUBLE PRECISION, DIMENSION(JM) ::
     &     BYP,BYPV,BYDAPO,BYPDA,DXCOSV,ONES,ONESPO
      DOUBLE PRECISION, DIMENSION(JM,LM) :: AX,BX,CX,DX,VX,EX
      DOUBLE PRECISION, DIMENSION(JM,LM_REQ) :: ARQX
      DOUBLE PRECISION, DIMENSION(LM_REQ) :: BYDPS,BYPKS
      DOUBLE PRECISION, DIMENSION(2,IMH+1) :: CN

      DOUBLE PRECISION, DIMENSION(JM,LM,2) :: DSJK
      DOUBLE PRECISION, DIMENSION(2,LM,2) :: DSHEM
      DOUBLE PRECISION, DIMENSION(LM,2) :: DSGLOB
      COMMON/WORK5/DSJK,DSHEM,DSGLOB

      INTEGER :: LINECT,JMHALF,INC,IHOUR0,IHOUR
      COMMON/DJLCOM/LINECT,JMHALF,INC,IHOUR0,IHOUR

      DOUBLE PRECISION, DIMENSION(JM,LM,NWAV_DAG,KAJLSP):: SPSTAD,SPTRAN
      DOUBLE PRECISION, DIMENSION(0:IMH,JM,LM) ::FCJKA,FCJKB,FCPVA,FCPVB
      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: AIJL2
      DOUBLE PRECISION, DIMENSION(LM) :: PM,PKM,PME
      DOUBLE PRECISION, DIMENSION(JM,2) :: PJ
      INTEGER, PARAMETER, DIMENSION(5) :: MW=(/1,2,3,6,9/)
      DOUBLE PRECISION, PARAMETER :: ONE=1.,P1000=1000.

      INTEGER ::
     &     I,J,J0,J1,JH,JHEMI,K,K1,KDN,KM,KMM1,KUP,L,LMP1,LR,M,N

      DOUBLE PRECISION ::
     &     BY100G,BYDP2,BYDPK,
     &     BYFSQ,BYIACN,BYIADA,BYIARD,BYIMDA,BYN,
     &     BYRCOS,DALPHA,DE4TI,DP,DPG,DPH,DPTI,
     &     DSSIG,DTHDP,DTHETA,EL4QI,
     &     FIMDA,GBYRSQ,GSQ,
     &     P1000K,PDN,PIG,PIH,PMK,
     &     PUP,PUTI,PVTI,SCALE,SCALES,SDDP,
     &     SKEI,SMALL,SN,SNAMI,SNDEGI,SNELGI,SQM,SQN,SZNDEG,SZNELG,
     &     THDN,THETA,THUP,TX,UDXN,UDXS,UX,WTKP1,WTN,WTNK,
     &     WTS,WTSK,XWON

C**** INITIALIZE CERTAIN QUANTITIES
      IF (KDIAG(2).GE.8) GO TO 120
      XWON=TWOPI/(DLON*FIM)
      LMP1=LM+1
      KM=LM
      KMM1=KM-1
      INC=1+(JM-1)/24
      JMHALF=JM/2
c      BYIM=1./FIM
      BY100G=.01/GRAV
c      SHA=RGAS/KAPA
      P1000K=P1000**KAPA
      DO 30 L=1,LM
      PKM(L)=PLM(L)**KAPA
      PME(L)=PSFMPT*SIGE(L)+PTOP
   30 PM(L)=PSFMPT*SIGE(L+1)+PTOP
      BYDPS(1)=1./(.5*PMTOP)
      BYDPS(2)=1./(.3*PMTOP)
      BYDPS(3)=1./(.2*PMTOP)
      BYPKS(1)=1./(.75*PMTOP)**KAPA
      BYPKS(2)=1./(.35*PMTOP)**KAPA
      BYPKS(3)=1./(.1*PMTOP)**KAPA
      DO 40 J=1,JM
      ONES(J)=1.
      ONESPO(J)=1.
c      BYDXYP(J)=1./DXYP(J)
      BYDAPO(J)=BYDXYP(J)
   40 CONTINUE
      ONESPO(1)=FIM
      ONESPO(JM)=FIM
      BYDAPO(1)=BYDAPO(1)*FIM
      BYDAPO(JM)=BYDAPO(JM)*FIM
      DO 50 J=2,JM
      DXCOSV(J)=DXV(J)*COSV(J)
   50 CONTINUE
      IHOUR0=TOFDY0+.5
      IHOUR=TOFDAY+.5
      LINECT=65
      WRITE (6,901)
      BYIACN=1./(IDACC(1)+1.D-20)
      BYIARD=1./(IDACC(2)+1.D-20)
      BYIADA=1./(IDACC(4)+1.D-20)
      BYIMDA=BYIADA*BYIM
      FIMDA=IDACC(4)*FIM
      DO 52 J=1,JM
      PJ(J,1)=0
      PJ(J,2)=0
      DO 51 K=1,KM
      PJ(J,1)=PJ(J,1)+AJK(J,K,1)
   51 PJ(J,2)=PJ(J,2)+AJK(J,K,2)
      BYP(J)=1./(PJ(J,1)+1.D-20)
      BYPDA(J)=BYP(J)*BYDXYP(J)
   52 BYPV(J)=1./(PJ(J,2)+1.D-20)
C****
C**** INITIALIZE DELTA SIGMA IN PRESSURE COORDINATES
C****
      DO 60 J1=1,2
      J0=J1-1
      DO 60 K=1,KM
      DPG=0.
      PIG=0.
      DO 58 JHEMI=1,2
      DPH=0.
      PIH=0.
      DO 55 JH=1,JMHALF
      J=(JHEMI-1)*(JMHALF-J0)+JH+J0
      DSJK(J,K,J1)=AJK(J,K,J1)/(PJ(J,J1)+1.D-20)
      DPH=DPH+AJK(J,K,J1)*WTJ(J,2,J1)
   55 PIH=PIH+PJ(J,J1)*WTJ(J,2,J1)
      DSHEM(JHEMI,K,J1)=DPH/(PIH+1.D-20)
      DPG=DPG+DPH
   58 PIG=PIG+PIH
   60 DSGLOB(K,J1)=DPG/(PIG+1.D-20)
C****
C**** V-V* IS D/DP(V'TH'/DTH/DP) , DX=4*V'TH'/DTH/DP AT INTERFACES
C****
      DO 70 J=2,JM
      KDN=1
      DO 70 K=1,KM
      CX(J,K)=FIM*IDACC(4)
      AX(J,K)=AJK(J,K,39)/(AJK(J,K,2)+1.D-20)
      KUP=K+1
      IF (K.EQ.KM) KUP=KM
      VX(J,K)=0.
      IF (AJK(J,K,2).EQ.0.) GO TO 70
      IF (AJK(J,KDN,2).EQ.0.) KDN=KDN+1
      VX(J,K)=AJK(J,K,2)*(AJK(J,KUP,39)/AJK(J,KUP,2)-
     *   AJK(J,KDN,39)/AJK(J,KDN,2))/(PM(KUP)-PM(KDN)+.5*
     *   (AJK(J,KUP,2)/AJK(J,KUP,24)-AJK(J,KDN,2)/AJK(J,KDN,24)))
   70 KDN=K
C     DX(J,K)=AX*CX=4*(TRANSFORMED STREAM FUNCTION-STREAM FUNCTION)
   90 DO 95 J=2,JM
      DX(J,KM)=AX(J,KM)*CX(J,KM)
      DO 95 K=1,KMM1
      WTKP1=AJK(J,K,2)/(AJK(J,K+1,2)+AJK(J,K,2)+1.D-20)
   95 DX(J,K)=(AX(J,K)*(1.-WTKP1)+AX(J,K+1)*WTKP1)*CX(J,K)
C****
C**** PROGNOSTIC QUANTITIES AT CONSTANT PRESSURE
C****
C**** # OF GRIDPOINTS, DELTA P, S.D. OF DELTA P
      CALL JLMAP (89,PLM,AJK(1,1,24),XWON*BYIADA,ONES,ONES,LS1-1,1,2)
      CALL JLMAP (91,PLM,AJK(1,1,2),BYIMDA,ONES,ONES,LS1-1,2,2)
      DO 98 J=2,JM
      DO 98 K=1,LS1-1
      BYN=1./(AJK(J,K,24)+1.D-10)
      AX(J,K)=0.
      SDDP=(AJK(J,K,23)-AJK(J,K,2)*AJK(J,K,2)*BYN)*BYN
   98 IF (SDDP.GT.0.) AX(J,K)=SQRT(SDDP)
      CALL JLMAP (68,PLM,AX,ONE,ONES,ONES,LS1-1,2,2)
C**** TEMPERATURE, HEIGHT, SPECIFIC AND RELATIVE HUMIDITY
      CALL JKMAPS (1,PLM,AJK(1,1,3),ONES,BYP,ONES,KM,2,1,
     *  ASJL,BYIMDA,ONESPO,ONES)
      SCALES=BYIMDA*BY100G
      CALL JKMAPS (2,PLM,AJK(1,1,4),BY100G,BYP,ONES,KM,2,1,
     *  ASJL(1,1,2),SCALES,ONESPO,ONES)
      SCALE=1.D5
      CALL JKMAP (3,PLM,AJK(1,1,5),SCALE,BYP,ONES,KM,2,1)
      SCALE=100.
      CALL JKMAP (4,PLM,AJK(1,1,7),SCALE,BYP,ONES,KM,2,1)
C**** LIQUID WATER CONTENT
      SCALE=1.D6
      CALL JKMAP (43,PLM,AJK(1,1,51),SCALE,BYP,ONES,KM,2,1)
C**** U AND V WINDS, STREAM FUNCTION
      SCALE=10.
      CALL JKMAP (5,PLM,AJK(1,1,8),SCALE,BYPV,ONES,KM,2,2)
      SCALE=100.
      CALL JKMAP (6,PLM,AJK(1,1,9),SCALE,BYPV,ONES,KM,2,2)
      DO 100 J=2,JM
      AX(J,1)=AJK(J,1,9)
      DO 100 K=2,KM
  100 AX(J,K)=AX(J,K-1)+AJK(J,K,9)
      SCALE=100.D-9*XWON*BYIADA/GRAV
      CALL JLMAP (61,PM,AX,SCALE,DXV,ONES,KM,2,2)
      DO 110 K=1,KM
      DO 110 J=2,JM
  110 BX(J,K)=AX(J,K)+.25*DX(J,K)
      CALL JLMAP (106,PM,BX,SCALE,DXV,ONES,KM,2,2)
C**** VERTICAL WINDS
      SCALE=-1.D5*BYIMDA
      CALL JLMAP (62,PME,AJK(1,1,25),SCALE,BYDXYP,ONES,KM,2,1)
C****
C**** CALCULATIONS FOR STANDING EDDIES
C****
      IF (SKIPSE.EQ.1.) GO TO 180
  120 DO 150 J=2,JM
      DO 150 K=1,KM
      DO 151 I=1,IM
      IF (AIJL(I,J,K,5).LE.1.D-20) THEN
        AIJL(I,J,K,1)=0.
        AIJL(I,J,K,2)=0.
        AIJL(I,J,K,3)=0.
        AIJL(I,J,K,4)=0.
      ENDIF
  151 AIJL2(I,J,K)=AIJL(I,J,K,2)/(AIJL(I,J,K,5)+1.D-20)
      EX(J,K)=0.
      AX(J,K)=0.
      BX(J,K)=0.
  150 CX(J,K)=0.
      DO 170 J=2,JM
      DO 155 I=1,IM
c**** COMPARE VERTICAL SUM OF AIJL AND AIJK
C     DEK(I,J)=0.
C     DEL(I,J)=0.
C     ELK(I,J)=0.
C     ELL(I,J)=0.
C     AMK(I,J)=0.
C     AML(I,J)=0.
C     DO 125 K=1,KM
C     DEL(I,J)=DEL(I,J)+AIJL(I,J,K,3)
C     ELL(I,J)=ELL(I,J)+AIJL(I,J,K,4)
C     AML(I,J)=AML(I,J)+AIJL(I,J,K,1)
C     DEK(I,J)=DEK(I,J)+AIJK(I,J,K,3)
C     ELK(I,J)=ELK(I,J)+AIJK(I,J,K,6)
C 125 AMK(I,J)=AMK(I,J)+AIJK(I,J,K,1)
C     IF(J.GT.40.OR.J.LT.38) GO TO 129
C     IF(I.GT.70) THEN
C       WRITE(6,126) I,J,DEK(I,J),DEL(I,J),ELK(I,J),ELL(I,J),
C    *               AMK(I,J),AML(I,J)
C 126   FORMAT(1X,'I J DEK DEL ELK ELL AMK AML= ',/1X,2I3,6E12.3)
C     ENDIF
C 129 CONTINUE
  155 SENDEG(I,J)=0.
      DO 170 K=1,KM
C**** SPECTRAL ANALYSIS OF STAND. AND TRANSIENT EDDY FLUXES
      DSSIG=DSIG(K)
      CALL FFT(AIJL2(1,J,K),FCPVA(0,J,K),FCPVB(0,J,K))
      CALL FFT(AIJL(1,J,K,3),FCJKA(0,J,K),FCJKB(0,J,K))
C     CALL FFTI(FCJKA(0,J,K),FCJKB(0,J,K),ACHK(1,J,K))
      DO 156 N=1,NWAV_DAG
  156 SPSTAD(J,K,N,1)=.5*FIM*(FCPVA(N-1,J,K)*FCJKA(N-1,J,K)+
     *   FCPVB(N-1,J,K)*FCJKB(N-1,J,K))/DSSIG
      CALL FFT(AIJL(1,J,K,4),FCJKA(0,J,K),FCJKB(0,J,K))
      DO 157 N=1,NWAV_DAG
  157 SPSTAD(J,K,N,2)=.5*FIM*(FCPVA(N-1,J,K)*FCJKA(N-1,J,K)+
     *   FCPVB(N-1,J,K)*FCJKB(N-1,J,K))/DSSIG
      CALL FFT(AIJL(1,J,K,1),FCJKA(0,J,K),FCJKB(0,J,K))
      DO 158 N=1,NWAV_DAG
  158 SPSTAD(J,K,N,3)=.5*FIM*(FCPVA(N-1,J,K)*FCJKA(N-1,J,K)+
     *   FCPVB(N-1,J,K)*FCJKB(N-1,J,K))/DSSIG
      DO 159 N=1,NWAV_DAG
      SPTRAN(J,K,N,1)=AJLSP(J,K,N,1)/DSSIG-SPSTAD(J,K,N,1)
      SPTRAN(J,K,N,2)=AJLSP(J,K,N,2)/DSSIG-SPSTAD(J,K,N,2)
  159 SPTRAN(J,K,N,3)=AJLSP(J,K,N,3)/DSSIG-SPSTAD(J,K,N,3)
      DPTI=0.
      PUTI=0.
      PVTI=0.
      DE4TI=0.
      EL4QI=0.
      SKEI=0.
      SNDEGI=0.
      SNELGI=0.
      SNAMI=0.
      DO 160 I=1,IM
      IF (AIJK(I,J,K,4).EQ.0.) GO TO 160
      DPTI=DPTI+AIJK(I,J,K,4)
      BYDPK=1./(AIJK(I,J,K,4)+1.D-20)
      PUTI=PUTI+AIJK(I,J,K,1)
      PVTI=PVTI+AIJK(I,J,K,2)
      DE4TI=DE4TI+AIJK(I,J,K,3)
      EL4QI=EL4QI+AIJK(I,J,K,6)
      SKEI=SKEI+(AIJK(I,J,K,1)*AIJK(I,J,K,1)
     *            +AIJK(I,J,K,2)*AIJK(I,J,K,2))*BYDPK
      SNDEGI=SNDEGI+(AIJK(I,J,K,3)*AIJK(I,J,K,2)*BYDPK)
      SENDEG(I,J)=SENDEG(I,J)
     *  +(AIJK(I,J,K,3)*AIJK(I,J,K,2)*BYDPK)
      SNELGI=SNELGI+(AIJK(I,J,K,6)*AIJK(I,J,K,2)*BYDPK)
      SNAMI=SNAMI+AIJK(I,J,K,1)*AIJK(I,J,K,2)*BYDPK
  160 CONTINUE
      AX(J,K)=SKEI-(PUTI*PUTI+PVTI*PVTI)/(DPTI+1.D-20)
      SZNDEG=DE4TI*PVTI/(DPTI+1.D-20)
      SZNELG=EL4QI*PVTI/(DPTI+1.D-20)
      DO 165 I=1,IM
  165 SENDEG(I,J)=SENDEG(I,J)-SZNDEG/FIM
      EX(J,K)=SNELGI-SZNELG
      BX(J,K)=SNDEGI-SZNDEG
  170 CX(J,K)=SNAMI-PUTI*PVTI/(DPTI+1.D-20)
      IF (KDIAG(2).GE.8) RETURN
C**** STANDING EDDY, EDDY AND TOTAL KINETIC ENERGY
  180 SCALE=50.D-4*BYIMDA/GRAV
      IF (SKIPSE.EQ.1.) GO TO 190
      CALL JKMAP (17,PLM,AX,SCALE,ONES,ONES,KM,2,2)
  190 DO 200 K=1,KM
      DO 200 J=2,JM
  200 AX(J,K)=AJK(J,K,11)-AJK(J,K,10)
      CALL JKMAP (18,PLM,AX,SCALE,ONES,ONES,KM,2,2)
      CALL JKMAP (19,PLM,AJK(1,1,11),SCALE,ONES,ONES,KM,2,2)
C**** POTENTIAL TEMPERATURE, POTENTIAL VORTICITY
      DO 205 LR=1,LM_REQ
      DO 205 J=1,JM
  205 ARQX(J,LR)=ASJL(J,LR,1)*BYIMDA*ONESPO(J)+273.16
      CALL JKMAPS (21,PLM,AJK(1,1,6),P1000K,BYP,ONES,KM,2,1,
     *  ARQX,P1000K,ONES,BYPKS)
      SCALE=1.D6*BYIMDA*P1000K
      CALL JLMAP (63,PLM,AJK(1,1,32),SCALE,BYDXYP,ONES,KM,2,1)
C****
C**** NORTHWARD TRANSPORTS AT CONSTANT PRESSURE
C****
C**** NORTHWARD TRANSPORT OF SENSIBLE HEAT BY EDDIES
      SCALE=25.D-14*XWON*SHA*BYIADA/GRAV
      DO 210 K=1,KM
      DO 210 J=2,JM
  210 AX(J,K)=AJK(J,K,13)-AJK(J,K,12)
      CALL JKMAP (13,PLM,AX,SCALE,DXV,ONES,KM,2,2)
C**** NORTHWARD TRANSPORT OF DRY STATIC ENERGY BY STANDING EDDIES,
C****   EDDIES, AND TOTAL
      SCALE=25.D-14*XWON*BYIADA/GRAV
      IF (SKIPSE.EQ.1.) GO TO 220
      CALL JKMAP (22,PLM,BX,SCALE,DXV,ONES,KM,2,2)
      DO 219 N=1,5
  219 CALL JLMAP(114+N,PLM,SPSTAD(1,1,N+1,1),SCALE,DXV,ONES,LM,2,2)
  220 DO 230 K=1,KM
      DO 230 J=2,JM
      AX(J,K)=SHA*(AJK(J,K,13)-AJK(J,K,12))+(AJK(J,K,15)-AJK(J,K,14))
  230 BX(J,K)=SHA*AJK(J,K,13)+AJK(J,K,15)
      CALL JKMAP (23,PLM,AX,SCALE,DXV,ONES,KM,2,2)
      DO 231 N=1,9
  231 CALL JLMAP(119+N,PLM,SPTRAN(1,1,N+1,1),SCALE,DXV,ONES,LM,2,2)
      SCALE=SCALE*.1
      CALL JKMAP (24,PLM,BX,SCALE,DXV,ONES,KM,2,2)
C**** NORTHWARD TRANSPORT OF LATENT HEAT BY STAND. EDDY, EDDIES AND TOTA
      DO 240 K=1,KM
      DO 240 J=2,JM
      DX(J,K)=AJK(J,K,17)-AJK(J,K,16)
  240 AX(J,K)=AX(J,K)+LHE*DX(J,K)
      SCALE=25.D-13*XWON*LHE*BYIADA/GRAV
      IF(SKIPSE.EQ.1.) GO TO 242
      CALL JKMAP(44,PLM,EX,SCALE,DXV,ONES,KM,2,2)
      DO 241 N=1,5
  241 CALL JLMAP(128+N,PLM,SPSTAD(1,1,N+1,2),SCALE,DXV,ONES,LM,2,2)
  242 CONTINUE
      CALL JKMAP (25,PLM,DX,SCALE,DXV,ONES,KM,2,2)
      DO 243 N=1,9
  243 CALL JLMAP(133+N,PLM,SPTRAN(1,1,N+1,2),SCALE,DXV,ONES,LM,2,2)
      SCALE=SCALE*.1
      CALL JKMAP (26,PLM,AJK(1,1,17),SCALE,DXV,ONES,KM,2,2)
C**** NORTHWARD TRANSPORT OF STATIC ENERGY BY EDDIES AND TOTAL
      DO 245 K=1,KM
      DO 245 J=2,JM
  245 DX(J,K)=BX(J,K)+LHE*AJK(J,K,17)
      SCALE=25.D-14*XWON*BYIADA/GRAV
      CALL JKMAP (27,PLM,AX,SCALE,DXV,ONES,KM,2,2)
      SCALE=SCALE*.1
      CALL JKMAP (28,PLM,DX,SCALE,DXV,ONES,KM,2,2)
C**** NORTHWARD TRANSPORT OF KINETIC ENERGY
      SCALE=50.D-12*XWON*BYIADA/GRAV
      CALL JKMAP (30,PLM,AJK(1,1,19),SCALE,DXV,ONES,KM,2,2)
C**** NOR. TRANS. OF ANG. MOMENTUM BY STANDING EDDIES, EDDIES AND TOTAL
      SCALE=100.D-18*XWON*RADIUS*BYIADA/GRAV
      IF (SKIPSE.EQ.1.) GO TO 250
      CALL JKMAP (31,PLM,CX,SCALE,DXCOSV,ONES,KM,2,2)
      DO 249 N=1,5
  249 CALL JLMAP(142+N,PLM,SPSTAD(1,1,N+1,3),SCALE,DXCOSV,ONES,LM,2,2)
  250 DO 260 K=1,KM
      DO 260 J=2,JM
      CX(J,K)=AJK(J,K,21)-AJK(J,K,20)
  260 DX(J,K)=AJK(J,K,21)+RADIUS*OMEGA*COSV(J)*AJK(J,K,9)
      CALL JKMAP (32,PLM,CX,SCALE,DXCOSV,ONES,KM,2,2)
      DO 261 N=1,9
  261 CALL JLMAP(147+N,PLM,SPTRAN(1,1,N+1,3),SCALE,DXCOSV,ONES,LM,2,2)
      SCALE=.1*SCALE
      CALL JKMAP (33,PLM,DX,SCALE,DXCOSV,ONES,KM,2,2)
C**** COMPARE SPSTAD AND SPTRAN OF D.E.
C     DO 267 J=31,34
C     DO 267 L=1,LM
C     WRITE(6,268) J,L,(SPSTAD(J,L,N,1),N=4,9),(SPTRAN(J,L,N,1),N=4,9)
C 267 CONTINUE
C 268 FORMAT(1X,'J L STAD DE= ',2I3,6E12.3/,5X,'TRAN DE= ',6X,6E12.3)
C****
C**** DYNAMIC CONVERGENCE OF ENERGY
C****
      SCALE=25.D-1*BYIMDA/GRAV
      DO 370 K=1,KM
      CX(1,K)=-BX(2,K)*SCALE*DXV(2)
      CX(JM,K)=BX(JM,K)*SCALE*DXV(JM)
      DX(1,K)=-(AJK(2,K,15)-AJK(2,K,14))*SCALE*DXV(2)
      DX(JM,K)=(AJK(JM,K,15)-AJK(JM,K,14))*SCALE*DXV(JM)
      DO 370 J=2,JM-1
      CX(J,K)=(BX(J,K)*DXV(J)-BX(J+1,K)*DXV(J+1))*SCALE
      DX(J,K)=((AJK(J,K,15)-AJK(J,K,14))*DXV(J) -
     *  (AJK(J+1,K,15)-AJK(J+1,K,14))*DXV(J+1))*SCALE
  370 CONTINUE
      SCALE=-100.D-1*BYIMDA/GRAV
      DO 380 K=1,KMM1
      DO 380 J=1,JM
      CX(J,K)=CX(J,K)-AJK(J,K,27)*SCALE
      CX(J,K+1)=CX(J,K+1)+AJK(J,K,27)*SCALE
      DX(J,K)=DX(J,K)-AJK(J,K,30)*SCALE
  380 DX(J,K+1)=DX(J,K+1)+AJK(J,K,30)*SCALE
      CALL JKMAP (14,PLM,CX,ONE,BYDXYP,ONES,KM,2,1)
      SCALE=100.
      CALL JKMAP (12,PLM,DX,SCALE,BYDXYP,ONES,KM,2,1)
C**** BAROCLINIC EKE GENERATION, P-K BY EDDY PRESSURE GRADIENT FORCE
      SCALE=50.E1*BYIMDA*RGAS/GRAV
      CALL JKMAP (9,PLM,AJK(1,1,31),SCALE,BYDXYP,ONES,KM,2,1)
      SCALE=100.E1*BYIMDA/(GRAV*2.*DT)
      CALL JKMAP (11,PLM,AJK(1,1,22),SCALE,ONES,ONES,KM,2,2)
C****
C**** VERTICAL TRANSPORTS
C****
C**** VERTICAL TRANSPORT OF GEOPOTENTIAL ENERGY BY EDDIES
      SCALE=-100.D-12*XWON*BYIADA/GRAV
      CALL JLMAP (99,PM,AJK(1,1,30),SCALE,ONES,ONES,KMM1,1,1)
C**** VERTICAL TRANSPORT OF DRY STATIC ENERGY BY EDDIES AND TOTAL
      DO 390 K=1,KMM1
      DO 390 J=1,JM
      AX(J,K)=AJK(J,K,27)-AJK(J,K,26)
  390 BX(J,K)=AJK(J,K,29)-AJK(J,K,28)
      SCALE=-100.D-12*XWON*BYIADA/GRAV
      CALL JLMAP (100,PM,AX,SCALE,ONES,ONES,KMM1,1,1)
      SCALE=SCALE*.01
      CALL JLMAP (101,PM,AJK(1,1,27),SCALE,ONES,ONES,KMM1,1,1)
C**** VERTICAL TRANSPORT OF LATENT HEAT BY EDDIES AND TOTAL
      SCALE=-100.D-12*XWON*LHE*BYIADA/GRAV
      CALL JLMAP (102,PM,BX,SCALE,ONES,ONES,KMM1,1,1)
      SCALE=SCALE*.1
      CALL JLMAP (103,PM,AJK(1,1,29),SCALE,ONES,ONES,KMM1,1,1)
C**** VERTICAL TRANSPORT OF STATIC ENERGY BY EDDIES AND TOTAL
      DO 420 K=1,KMM1
      DO 420 J=1,JM
      AX(J,K)=AX(J,K)+LHE*BX(J,K)
  420 BX(J,K)=AJK(J,K,27)+LHE*AJK(J,K,29)
      SCALE=-100.D-13*XWON*BYIADA/GRAV
      CALL JLMAP (104,PM,AX,SCALE,ONES,ONES,KMM1,1,1)
      SCALE=SCALE*.1
      CALL JLMAP (105,PM,BX,SCALE,ONES,ONES,KMM1,1,1)
C**** VERTICAL TRANSPORT OF KINETIC ENERGY
      SCALE=-12.5E-11*XWON*BYIADA/GRAV
      CALL JLMAP (110,PM,AJK(1,1,36),SCALE,ONES,ONES,KMM1,1,2)
C**** VERTICAL TRANSPORT OF ANGULAR MOMENTUM BY LARE SCALE MOTIONS
      SCALE=-25.D-16*XWON*RADIUS*BYIADA/GRAV
      CALL JLMAP (111,PM,AJK(1,1,37),SCALE,COSV,ONES,KMM1,1,2)
      SCALE=1.D-2*SCALE
      CALL JLMAP (112,PM,AJK(1,1,38),SCALE,COSV,ONES,KMM1,1,2)
C**** VERTICAL TRANSPORT OF POTENTIAL VORTICITY TOTAL AND BY EDDIES
      SCALE=-25.D-4*XWON*P1000K*BYIADA/GRAV
      CALL JLMAP (107,PM,AJK(1,1,33),SCALE,BYDXYP,ONES,KMM1,1,1)
      CALL JLMAP (108,PM,AJK(1,1,34),SCALE,BYDXYP,ONES,KMM1,1,1)
C**** NOR. TRANSPORT OF QUASI-GEOSTROPHIC POT. VORTICITY BY EDDIES
      DO 490 K=1,KM
      AX(1,K)=0.
      AX(JM,K)=0.
      DX(1,K)=0.
      DX(JM,K)=0.
      DO 490 J=2,JM-1
      AX(J,K)=((AJK(J,K,21)-AJK(J,K,20))*DXCOSV(J)/(AJK(J,K,2)+1.D-20)-
     *  (AJK(J+1,K,21)-AJK(J+1,K,20))*DXCOSV(J+1)/
     *  (AJK(J+1,K,2)+1.D-20))/COSP(J)
      DX(J,K)=FCOR(J)*(VX(J,K)+VX(J+1,K))/(AJK(J,K,2)+AJK(J+1,K,2)+
     *  1.D-20)
  490 CONTINUE
      DO 500 K=1,KM
      DO 500 J=2,JM-1
  500 AX(J,K)=AJK(J,K,1)*(AX(J,K)+.25*DX(J,K))
      SCALE=1.D6
      CALL JKMAP (10,PLM,AX,SCALE,BYPDA,ONES,KM,2,1)
C****
C**** ELIASSEN PALM FLUX:  NORTHWARD, VERTICAL, DIVERGENCE
C****
      SCALE=25.D-11*XWON*BYIADA/GRAV
      SMALL=1.D-20
      DO 510 K=1,KM
      AX(1,K)=0.
      DO 510 J=2,JM
      UX=AJK(J,K,8)/(AJK(J,K,2)+1.D-20)
      IF (ABS(UX).GE.SMALL) GO TO 510
      SN=+1.
      IF (UX.LT.0.) SN=-1.
      UX=SN*SMALL
  510 AX(J,K)=(AJK(J,K,15)-AJK(J,K,14))/UX*DXV(J)
      CALL JKMAP (34,PLM,AX,SCALE,ONES,ONES,KM,2,2)
      SCALE=-100.D-11*XWON*BYIADA/GRAV
      DO 520 K=1,KMM1
      BX(1,K)=0.
      BX(JM,K)=0.
      DO 520 J=1,JM
      IF (J.NE.1.AND.J.NE.JM) GO TO 516  ! corrected 4-25-2000
      IF (J.EQ.1) UX=.5*(AJK(J+1,K,8)+AJK(J+1,K+1,8))/
     *     (AJK(J+1,K,2)+AJK(J+1,K+1,2)+1.D-20)
      IF (J.EQ.JM) UX=.5*(AJK(J,K,8)+AJK(J,K+1,8))/
     *     (AJK(J,K,2)+AJK(J,K+1,2)+1.D-20)
      GO TO 518
  516 UX=(AJK(J,K,8)+AJK(J+1,K,8)+AJK(J,K+1,8)+AJK(J+1,K+1,8))/
     *  (AJK(J,K,2)+AJK(J+1,K,2)+AJK(J,K+1,2)+AJK(J+1,K+1,2)+1.D-20)
  518 IF (ABS(UX).GE.SMALL) GO TO 520
      SN=+1.
      IF (UX.LT.0.) SN=-1.
      UX=SN*SMALL
  520 BX(J,K)=AJK(J,K,30)/UX
      CALL JLMAP (113,PM,BX,SCALE,ONES,ONES,KMM1,1,1)
      DO 530 K=1,KM
      CX(1,K)=0.
      CX(JM,K)=0.
      DO 530 J=2,JM-1
  530 CX(J,K)=.25*(AX(J+1,K)-AX(J,K))
      DO 540 K=1,KMM1
      DO 540 J=1,JM
      CX(J,K)=CX(J,K)-BX(J,K)
  540 CX(J,K+1)=CX(J,K+1)+BX(J,K)
      SCALE=1.D5
      CALL JKMAP (36,PLM,CX,SCALE,BYPDA,ONES,KM,2,1)
C****
C**** D/DY OF Q-G POTENTIAL VORTICITY AND REFRACTION INDICES
C****
C**** PRELIMINARIES:  VERTICAL DERIVATIVES AND N**2
      GSQ=GRAV*GRAV
      GBYRSQ=GRAV*GRAV/(RGAS*RGAS)
      IF (AJK(2,KM,1).LT.1.D-20) GO TO 670  ! ISTART=4,...
      DO 600 J=1,JM
      K1=1
  580 IF (AJK(J,K1,1).GT.1.D-20) GO TO 590
      AX(J,K1)=0.
      BX(J,K1)=0.
      DX(J,K1)=0.
      K1=K1+1
      GO TO 580
  590 KDN=K1
      PDN=PM(KDN)+.5*AJK(J,KDN,1)/(AJK(J,KDN,35)+1.D-20)
      DO 600 K=K1,KM
      DP=AJK(J,K,1)
      PMK=PM(K)+.5*AJK(J,K,1)/(AJK(J,K,35)+1.D-20)
      KUP=K+1
      IF (K.EQ.KM) KUP=KM
      PUP=PM(KUP)+.5*AJK(J,KUP,1)/(AJK(J,KUP,35)+1.D-20)
      DALPHA=(AJK(J,KUP,3)/(AJK(J,KUP,1)+1.D-20)+273.16)/PUP-
     *  (AJK(J,KDN,3)/(AJK(J,KDN,1)+1.D-20)+273.16)/PDN
      DTHETA=AJK(J,KUP,6)/(AJK(J,KUP,1)+1.D-20)-
     *  AJK(J,KDN,6)/(AJK(J,KDN,1)+1.D-20)
      THETA=AJK(J,K,6)/(AJK(J,K,1)+1.D-20)
      TX=AJK(J,K,3)/(AJK(J,K,1)+1.D-20)+273.16
      IF (ABS(DTHETA).GE.1.D-20) GO TO 595
      SN=+1.
      IF (DTHETA.LT.0.) SN=-1.
      DTHETA=SN*1.D-20
  595 DX(J,K)=DP*FCOR(J)*PMK*THETA*(PUP-PDN)/(TX*DTHETA*DXYP(J))
      AX(J,K)=DALPHA/(PUP-PDN-1.D-20)
C**** CALCULATE N**2 AT PRESSURE LATITUDES
      BX(J,K)=-DP*GSQ*PMK*DTHETA/(RGAS*TX*THETA*(PUP-PDN-1.D-20))
      KDN=K
  600 PDN=PMK
C**** CALCULATE  Q12 = (D(UDX) + F*DA)/DA
      DO 620 K=1,KM
      UDXS=0.
      DO 610 J=1,JM-1
      UDXN=AJK(J+1,K,8)/(AJK(J+1,K,2)+1.D-20)*DXV(J+1)
      CX(J,K)=(UDXS-UDXN+FCOR(J))/DXYP(J)
  610 UDXS=UDXN
      CX(JM,K)=(UDXS+FCOR(JM))/DXYP(JM)
C**** FIND DQ/DY = (Q12(J)-Q12(J-1)+Q3(J)-Q3(J-1))/DY
      DO 620 J=JM,2,-1
      DP=AJK(J,K,2)
      AX(J,K)=DP*(CX(J,K)-CX(J-1,K) + (AX(J,K)-AX(J-1,K))*
     *  (DX(J,K)+DX(J-1,K))/(AJK(J,K,1)+AJK(J-1,K,1)+1.D-20))/DYP(3)
  620 CONTINUE
      SCALE=1.D12
      CALL JKMAP (42,PLM,AX,SCALE,BYPV,ONES,KM,2,2)
C**** TERMS FOR THE REFRACTION INDEX EXPRESSION
      DO 640 J=2,JM
      BYFSQ=2.*DXYV(J)*DXYV(J)/(FCOR(J-1)*FCOR(J-1)+FCOR(J)*FCOR(J))
      DO 640 K=1,KM
      BYDP2=1./(AJK(J-1,K,1)+AJK(J,K,1)+1.D-20)
      TX=BYDP2*(AJK(J-1,K,3)+AJK(J,K,3))+273.16
      DX(J,K)=GBYRSQ/(TX*TX)
      SQN=BYDP2*(BX(J-1,K)+BX(J,K))
      CX(J,K)=SQN*BYFSQ
      UX=AJK(J,K,8)
      IF (ABS(UX).GE.1.D-20) GO TO 635
      SN=+1.
      IF (UX.LT.0.) SN=-1.
      UX=SN*1.D-20
  635 AX(J,K)=AX(J,K)/UX
  640 CONTINUE
C**** COMBINE TERMS, PRINT OUT REFRACTION INDICES
      SCALE=1.D8
      DO 660 M=1,5
      SQM=MW(M)*MW(M)
      DO 650 J=2,JM
      BYRCOS=1./(RADIUS*RADIUS*COSV(J)*COSV(J))
      DO 650 K=1,KM
      DP=AJK(J,K,2)
  650 BX(J,K)=DP*(CX(J,K)*(AX(J,K)-SQM*BYRCOS)-.25*DX(J,K))
  660 CALL JKMAP (M+36,PLM,BX,SCALE,BYPV,ONES,KM,2,2)
  670 CONTINUE
C**** SKIP REMAINING MAPS IF DATA NOT AVAILABLE
      IF (AJK(1,1,44).NE.0.) GO TO 799
C****
C**** CHANGE OF THE MEAN FIELDS OF WIND AND TEMPERATURE
C****
      DO 710 K=1,KM
      DO 710 J=1,JM
      AX(J,K)=0.
      BX(J,K)=0.
      CX(J,K)=0.
  710 DX(J,K)=0.
      DO 720 K=2,KMM1
      DO 720 J=2,JM-1
      AX(J,K)=((AJK(J,K,21)-AJK(J,K,20))*DXV(J)-(AJK(J+1,K,21)-
     *  AJK(J+1,K,20))*DXV(J+1))/(AJK(J,K,1)*DXYP(J)+1.D-20)+
     *  .125*((AJK(J,K,37)-AJK(J,K-1,37))/(AJK(J,K,2)*DXYV(J)+1.D-20)+
     * (AJK(J+1,K,37)-AJK(J+1,K-1,37))/(AJK(J+1,K,2)*DXYV(J+1)+1.D-20))
      BX(J,K)=(AJK(J+1,K,44)*DXCOSV(J+1)-AJK(J,K,44)*DXCOSV(J))/
     *  (DXYP(J)*COSP(J))+.5*(AJK(J,K-1,45)-AJK(J,K,45)+
     *  AJK(J+1,K-1,45)-AJK(J+1,K,45))/(PM(K-1)-PM(K))
      CX(J,K)=.25*((AJK(J,K,13)-AJK(J,K,12))*DXV(J)-(AJK(J+1,K,13)-
     *  AJK(J+1,K,12))*DXV(J+1))/(AJK(J,K,1)*DXYP(J)+1.D-20)+
     *  (AJK(J,K,50)-AJK(J,K-1,50))/(PM(K-1)-PM(K))*BYIADA*PKM(K)
  720 CONTINUE
C**** WIND: RATE OF CHANGE, ADVECTION, EDDY CONVERGENCE
      IF (IDACC(4).LE.1) GO TO 730
      SCALE=1.D6/((IDACC(4)-1)*DT*NDAA)
      CALL JLMAP (40,PLM,AJK(1,1,47),SCALE,ONES,ONES,KM,2,2)
  730 SCALE=1.D6*BYIADA
      CALL JLMAP (59,PLM,AJK(1,1,40),SCALE,ONES,ONES,KMM1,2,1)
      SCALE=1.D6
      CALL JLMAP (60,PLM,AX,SCALE,ONES,ONES,KMM1,2,1)
C**** WIND: TRANSFORMED ADVECTION, LAGRANGIAN CONVERGENCE (DEL.F)
      SCALE=1.D6*BYIADA
      CALL JLMAP (64,PLM,AJK(1,1,42),SCALE,ONES,ONES,KMM1,2,1)
      CALL JLMAP (65,PLM,BX,SCALE,ONES,ONES,KMM1,2,1)
C**** WIND: DU/DT BY STRAT. DRAG
      SCALE=1.D6/(FIM*IDACC(1)*DT*NCNDS+1.E-20)
      CALL JLMAP (1,PLM,AJL(1,1,20),SCALE,ONES,ONES,LM,2,2)
C**** TEMPERATURE: RATE OF CHANGE, ADVECTION, EDDY CONVERGENCE
      IF (IDACC(4).LE.1) GO TO 750
      SCALE=1.D1*SDAY/((IDACC(4)-1)*DT*NDAA)
      CALL JLMAP (66,PLM,AJK(1,1,49),SCALE,ONES,PKM,KM,2,1)
  750 SCALE=1.D1*SDAY*BYIADA
      CALL JLMAP (67,PLM,AJK(1,1,41),SCALE,ONES,PKM,KMM1,2,1)
      SCALE=1.D1*SDAY
      CALL JLMAP (69,PLM,CX,SCALE,ONES,ONES,KMM1,2,1)
C**** TEMPERATURE: TRANSFORMED ADVECTION
      SCALE=1.D1*SDAY*BYIADA
      CALL JLMAP (70,PLM,AJK(1,1,43),SCALE,ONES,PKM,KMM1,2,1)
  799 CONTINUE
      RETURN
  901 FORMAT (
     *  '010**14 WATTS = .2067 * 10**19 CALORIES/DAY'/
     *  ' 10**18 JOULES = .864 * 10**30 GM*CM**2/SEC/DAY')
      END SUBROUTINE DIAGJK

      SUBROUTINE JKMAP (NT,PM,AX,SCALE,SCALEJ,SCALEK,KMAX,JWT,J1)
      USE E001M12_COM, only : jm,lm,
     &     JDATE,JDATE0,JMNTH0,JMONTH,JYEAR,JYEAR0,XLABEL
      USE GEOM, only : JLAT,WTJ
      IMPLICIT NONE

      INTEGER, DIMENSION(JM) :: MLAT
      DOUBLE PRECISION, DIMENSION(JM) :: FLAT,ASUM
      DOUBLE PRECISION, DIMENSION(2) :: AHEM
      DOUBLE PRECISION, DIMENSION(JM,LM) :: CX

      DOUBLE PRECISION, DIMENSION(JM,LM,2) :: DSJK
      DOUBLE PRECISION, DIMENSION(2,LM,2) :: DSHEM
      DOUBLE PRECISION, DIMENSION(LM,2) :: DSGLOB
      COMMON/WORK5/DSJK,DSHEM,DSGLOB

      INTEGER :: LINECT,JMHALF,INC,IHOUR0,IHOUR
      COMMON/DJLCOM/LINECT,JMHALF,INC,IHOUR0,IHOUR

      INTEGER :: J1,JWT,KMAX,NT
      DOUBLE PRECISION :: SCALE,SCALER
      DOUBLE PRECISION, DIMENSION(JM,LM) :: AX
      DOUBLE PRECISION, DIMENSION(JM,3) :: ARQX ! pass 3 as arg ??
      DOUBLE PRECISION, DIMENSION(JM) :: SCALEJ,SCALJR
      DOUBLE PRECISION, DIMENSION(LM) :: SCALEK
      DOUBLE PRECISION, DIMENSION(3) :: SCALLR ! pass 3 as arg ??
      DOUBLE PRECISION, DIMENSION(LM+3) :: PM ! pass 3 as arg ??

      COMMON/DJKTTL/TITLE(1)
      CHARACTER*4 DASH,WORD(4),TITLE*64
      DATA DASH/'----'/,WORD/'SUM','MEAN',' ','.1*'/

      INTEGER :: IWORD,J,J0,JH,JHEMI,JY,JY0,K,L
      DOUBLE PRECISION :: AGLOB,FGLOB,FLATJ,G1,H1,H2,SUMFAC

C****
C**** PRODUCE A LATITUDE BY LAYER TABLE OF THE ARRAY A
C****
   10 LINECT=LINECT+KMAX+7
      IF (LINECT.LE.60) GO TO 20
      JY0=JYEAR0-1900
      JY=JYEAR-1900
      WRITE (6,907) XLABEL(1:108),JDATE0,JMNTH0,JY0,JDATE,JMONTH,JY
      LINECT=KMAX+8
   20 WRITE (6,901) TITLE(NT),(DASH,J=J1,JM,INC)
      WRITE (6,904) WORD(JWT),(JLAT(J,J1),J=JM,J1,-INC)
      WRITE (6,905) (DASH,J=J1,JM,INC)
      J0=J1-1
  100 DO 110 J=J1,JM
      DO 110 K=1,KMAX
  110 CX(J,K)=AX(J,K)*SCALE*SCALEJ(J)*SCALEK(K)
C**** HORIZONTAL SUMS AND TABLE ENTRIES
      DO 140 K=KMAX,1,-1
      AGLOB=0.
      DO 130 JHEMI=1,2
      AHEM(JHEMI)=0.
      DO 120 JH=1,JMHALF
      J=(JHEMI-1)*(JMHALF-J0)+JH+J0
      FLAT(J)=CX(J,K)/(DSJK(J,K,J1)+1.D-20)
      MLAT(J)=NINT(FLAT(J))
  120 AHEM(JHEMI)=AHEM(JHEMI)+CX(J,K)*WTJ(J,JWT,J1)
  130 AGLOB=AGLOB+AHEM(JHEMI)/JWT
      H1=AHEM(1)/(DSHEM(1,K,J1)+1.D-20)
      H2=AHEM(2)/(DSHEM(2,K,J1)+1.D-20)
      G1=AGLOB/(DSGLOB(K,J1)+1.D-20)
      WRITE (6,902) PM(K),G1,H2,H1,(MLAT(J),J=JM,J1,-INC)
         IF (NT.EQ.5) CALL KEYJKJ (K,FLAT)
  140 CONTINUE
C**** VERTICAL SUMS
      WRITE (6,905) (DASH,J=J1,JM,INC)
C     IF (NT.GE.80.AND.NT.LE.87) RETURN
      SUMFAC=1.
      IWORD=3
      IF (NT.NE.1.AND.NT.NE.6.AND.NT.NE.24.AND.NT.NE.26.AND.NT.NE.28
     *  .AND.NT.NE.33) GO TO 160
      SUMFAC=10.
      IWORD=4
  160 CONTINUE
      DO 180 J=J1,JM
      ASUM(J)=0.
      DO 170 K=1,KMAX
  170 ASUM(J)=ASUM(J)+CX(J,K)
  180 MLAT(J)=NINT(ASUM(J)*SUMFAC)
      AGLOB=0.
      DO 200 JHEMI=1,2
      AHEM(JHEMI)=0.
      DO 190 JH=1,JMHALF
      J=(JHEMI-1)*(JMHALF-J0)+JH+J0
      AHEM(JHEMI)=AHEM(JHEMI)+ASUM(J)*WTJ(J,JWT,J1)*SUMFAC
  190 CONTINUE
  200 AGLOB=AGLOB+AHEM(JHEMI)/JWT
      WRITE (6,903) WORD(IWORD),AGLOB,AHEM(2),AHEM(1),
     *  (MLAT(J),J=JM,J1,-INC)
         IF (NT.EQ.1) CALL KEYJKT (AGLOB,ASUM)
         IF (NT.EQ.18.OR.NT.EQ.19) CALL KEYJKE (NT,AHEM,ASUM)
         IF (NT.GE.22.AND.NT.LE.33) CALL KEYJKN (NT,ASUM,SUMFAC)
      RETURN
C****
      ENTRY JKMAPS (NT,PM,AX,SCALE,SCALEJ,SCALEK,KMAX,JWT,J1,
     *  ARQX,SCALER,SCALJR,SCALLR)
      LINECT=LINECT+KMAX+10
      IF (LINECT.LE.60) GO TO 230
      JY0=JYEAR0-1900
      JY=JYEAR-1900
      WRITE (6,907) XLABEL(1:108),JDATE0,JMNTH0,JY0,JDATE,JMONTH,JY
      LINECT=KMAX+11
  230 J0=J1-1
C**** PRODUCE UPPER STRATOSPHERE NUMBERS FIRST
      WRITE (6,901) TITLE(NT),(DASH,J=J1,JM,INC)
      WRITE (6,904) WORD(JWT),(JLAT(J,J1),J=JM,J1,-INC)
      WRITE (6,905) (DASH,J=J1,JM,INC)
      DO 260 L=3,1,-1
      FGLOB=0.
      DO 250 JHEMI=1,2
      AHEM(JHEMI)=0.
      DO 240 JH=1,JMHALF
      J=(JHEMI-1)*(JMHALF-J0)+JH-J0
      FLATJ=ARQX(J,L)*SCALER*SCALJR(J)*SCALLR(L)
      MLAT(J)=NINT(FLATJ)
  240 AHEM(JHEMI)=AHEM(JHEMI)+FLATJ*WTJ(J,JWT,J1)
  250 FGLOB=FGLOB+AHEM(JHEMI)/JWT
  260 WRITE (6,902) PM(L+LM),FGLOB,AHEM(2),AHEM(1),
     *  (MLAT(J),J=JM,J1,-INC)
      GO TO 100
  901 FORMAT ('0',30X,A64,'  CP'/1X,30('-'),24A4)
  902 FORMAT (1X,F5.1,3F8.1,1X,24I4)
  903 FORMAT (A6,3F8.1,1X,24I4)
  904 FORMAT (' P(MB) ',A4,' G      NH      SH  ',24I4)
  905 FORMAT (1X,30('-'),24A4)
  907 FORMAT ('1',A,I4,1X,A3,I3,' TO',I3,1X,A3,I3)
      END SUBROUTINE JKMAP

      BLOCK DATA BDJL
C****
C**** TITLES FOR SUBROUTINE DIAGJL
C****
      COMMON/DJLTTL/
     *  TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6,TITLE7,
     *  TITLE8,TITLE9,TITLEA,TITLEB,TITLEC,TITLED,TITLEE,
     *  TITLEF,TITLEG,TITLEH
      CHARACTER*64 TITLE1(13)
      DATA TITLE1/
C****                                                              1-13
     1'ZONAL WIND CHANGE BY STRATOSPHERIC DRAG  (10**-6 M S-2)        ',
     2'HEIGHT (HUNDREDS OF METERS)                                    ',
     3'SPECIFIC HUMIDITY (10**-5 KG H2O/KG AIR)                       ',
     4'RELATIVE HUMIDITY (PERCENT)                                    ',
     5'ZONAL WIND (U COMPONENT) (TENTHS OF METERS/SECOND)             ',
     6'MERIDIONAL WIND (V COMPONENT) (HUNDREDTHS OF METERS/SECOND)    ',
     7'STREAM FUNCTION (10**9 KILOGRAMS/SECOND)                       ',
     8'VERTICAL VELOCITY (10**-5 MILLIBARS/SECOND)                    ',
     9'BAROCLINIC EDDY KINETIC ENERGY GEN. (10**-1 WATTS/M**2/SIGMA)  ',
     A'VERTICAL MASS EXCHANGE FROM MOIST CONVECTION (10**9 KG/SECOND) ',
     B'SOLAR RADIATION HEATING RATE (HUNDREDTHS OF DEGREES KELVIN/DAY)',
C    *   'EES KELVIN/DAY)',
     C'THERMAL RADIATION COOLING RATE (HUNDREDTHS OF DEGREES K/DAY)   ',
     D'TOTAL RADIATION COOLING RATE (10**13 WATTS/UNIT SIGMA)         '/
      CHARACTER*64 TITLE2(8)
      DATA TITLE2/
C****                                                             14-21
     1'HEATING BY LARGE SCALE CONDENSATION (10**13 WATTS/UNIT SIGMA)',
     2'HEATING BY DRY CONVECTION (10**13 WATTS/UNIT SIGMA)          ',
     6'CHANGE OF LATENT HEAT BY DRY CONV.  (10**14 W/UNIT SIGMA)    ',
     4'STANDING EDDY KINETIC ENERGY (10**4 JOULES/M**2/UNIT SIGMA)  ',
     5'EDDY KINETIC ENERGY (10**4 JOULES/M**2/UNIT SIGMA)           ',
     6'TOTAL KINETIC ENERGY (10**4 JOULES/M**2/UNIT SIGMA)          ',
     7'AVAILABLE POTENTIAL ENERGY (10**5 JOULES/M**2/UNIT SIGMA)    ',
     8'POTENTIAL TEMPERATURE (DEGREES KELVIN)                       '/
      CHARACTER*64 TITLE3(7)
      DATA TITLE3/
C****                                                             22-28
     2'NOR. TRANS. OF DRY STAT. ENERGY BY STAND. EDDIES (10**14 W/DSIG)'
     *,
     3'NORTH. TRANS. OF DRY STATIC ENERGY BY EDDIES (10**14 WATTS/DSIG)'
     *,
     4'TOTAL NORTH. TRANSPORT OF DRY STATIC ENERGY (10**15 WATTS/DSIG) '
     *,
     5'NORTHWARD TRANSPORT OF LATENT HEAT BY EDDIES (10**13 WATTS/DSIG)'
     *,
     5'TOTAL NORTHWARD TRANSPORT OF LATENT HEAT (10**14 WATTS/UNIT SIG)'
     *,
     6'NORTH.TRANSPORT OF STATIC ENERGY BY EDDIES (10**14 WATTS/DSIGMA)'
     *,
     9'TOTAL NORTHWARD TRANSPORT OF STATIC ENERGY (10**15 WATTS/DSIGMA)'
     *     /
      CHARACTER*64 TITLE4(5)
      DATA TITLE4/
C****                                                             29-33
     9'NORTH. TRANSPORT OF KINETIC ENERGY BY EDDIES (10**12 WATTS/DSIG)'
     *,
     *'TOTAL NORTHWARD TRANSPORT OF KINETIC ENERGY (10**12 WATTS/DSIG) '
     *,
     1'NORTH. TRANS. OF ANG. MOMENTUM BY STAND. EDDIES (10**18 J/DSIG) '
     *,
     2'NORTH. TRANS. OF ANG. MOMENTUM BY EDDIES (10**18 JOULES/DSIGMA) '
     *,
     *'TOTAL NORTHWARD TRANSPORT OF ANG. MOMENTUM (10**19 JOULES/DSIG) '
     *     /
      CHARACTER*64 TITLE5(6)
      DATA TITLE5/
C****                                                             34-39
     4'VERT. TRANS. OFDRY STATIC ENERGY BY EDDIES (10**12 WATTS)      ',
     5'TOT. LARGE SCALE VERT. TRANS. OF DRY STAT. ENER. (10**14 WATTS)',
     6'VERTICAL TRANSPORT OF LATENT HEAT BY EDDIES (10**12 WATTS)     ',
     7'TOTAL LARGE SCALE VERT. TRANS. OF LATENT HEAT (10**13 WATTS)   ',
     8'VERTICAL TRANSPORT OF STATIC ENERGY BY EDDIES (10**13 WATTS)   ',
     9'TOTAL LARGE SCALE VERT. TRANS. OF STATIC ENERGY (10**14 WATTS) '/
      CHARACTER*64 TITLE6(4)
      DATA TITLE6/
C****                                                             40-43
     *'DU/DT   TOTAL CHANGE  (10**-6 M/S/S)                 CP        ',
     1'TOTAL LARGE SCALE VERT. TRANS. OF KINETIC ENERGY (10**11 WATTS)',
     2'VERT. TRANS. OFANG. MOMENTUM BY EDDIES (10**16JOULES)          ',
     3'TOTAL LARGE SCALE VERT. TRANS. OF ANG. MOMENTUM (10**18 JOULES)'/
      CHARACTER*64 TITLE7(9)
      DATA TITLE7/
C****                                                             44-52
     4'CHANGE OF ANG. MOMENTUM BY DRY CONVEC (10**18 JOULE/UNIT SIGMA)',
     5'CHANGE OF ANG. MOMENTUM BY MOIST CONV (10**18 JOULE/UNIT SIGMA)',
     6'CHANGE OF ANG. MOMENTUM BY DIFFUSION (10**18 JOULES/UNIT SIGMA)',
     7'U WIND AVERAGED OVER I,EAST PACIFIC (TENTHS OF METERS/SECOND)',
     8'V WIND AVERAGED OVER EAST PACIFIC (TENTHS OF METERS/SECOND)',
     9'VERTICAL VELOCITY FOR EAST PACIFIC (10**-5 METERS/SECOND)',
     A'U WIND AVERAGED OVER I=WEST PACIFIC (TENTHS OF METERS/SECOND)',
     B'V WIND AVERAGED OVER I=WEST PACIFIC (TENTHS OF METERS/SECOND)',
     C'VERTICAL VELOCITY FOR WEST PACIFIC (10**-5 METERS/SECOND)'/
      CHARACTER*64 TITLE8(8)
      DATA TITLE8/
C****                                                             53-60
     3'POTENTIAL VORTICITY (10**-6 K/(MB-S))                          ',
     4'NORTHWARD TRANSPORT OF Q-G POT. VORTICITY  (10**18 JOULES/DSIG)',
     1'P-K BY EDDY PRESSURE GRADIENT FORCE  (10**-1 W/M**2/UNIT SIGMA)',
     6'Q-G POT. VORTICITY CHANGE OVER LATITUDES (10**-12 1/(SEC-M))   ',
     7'TRANSFORMED STREAM FUNCTION  (10**9 KILOGRAMS/SECOND)          ',
     8'DYNAMIC CONVERGENCE OF EDDY GEOPOTENTIAL (.1 WATTS/M**2/DSIGMA)',
     *'DU/DT BY MEAN ADVECTION  (10**-6 M/S/S)              CP        ',
     *'DU/DT BY EDDY CONVERGENCE  (10**-6 M/S/S)            CP        '/
      CHARACTER*64 TITLE9(12)
      DATA TITLE9/
C****                                                             61-72
     1'STREAM FUNCTION (10**9 KILOGRAMS/SECOND)             CP  ',
     2'VERTICAL VELOCITY (10**-5 MILLIBARS/SECOND)              ',
     3'POTENTIAL VORTICITY (10**-6 K/(MB-S))                CP  ',
     4'DU/DT BY TRANSFORMED ADVECTION  (10**-6 M/S/S)       CP  ',
     5'DU/DT BY ELIASSEN-PALM DIVERGENCE  (10**-6 M/S/S)    CP  ',
     6'DTEMP/DT   TOTAL CHANGE  (10**-1 DEG-K/DAY)          CP  ',
     7'DTEMP/DT BY MEAN ADVECTION  (10**-1 DEG-K/DAY)       CP  ',
     8'STANDARD DEVIATION OF PRESSURE DIFFERENCES  (MB)         ',
     9'DTEMP/DT BY EDDY CONVERGENCE  (10**-1 DEG-K/DAY)     CP  ',
     O'DTEMP/DT BY TRANSFORMED ADVECTION (10**-1 DEG-K/DAY) CP  ',
     1'REFRACTION INDEX FOR WAVE NUMBER 1  (10**-8 PER METER**2)',
     2'REFRACTION INDEX FOR WAVE NUMBER 2  (10**-8 PER METER**2)'/
      CHARACTER*64 TITLEA(11)
      DATA TITLEA/
C****                                                             73-83
     3'REFRACTION INDEX FOR WAVE NUMBER 3  (10**-8 PER METER**2)  ',
     4'REFRACTION INDEX FOR WAVE NUMBER 6  (10**-8 PER METER**2)  ',
     5'REFRACTION INDEX FOR WAVE NUMBER 9  (10**-8 PER METER**2)  ',
     6'                                                           ',
     7'TOTAL CLOUD COVER (PERCENT)                                ',
     8'SUPER SATURATION CLOUD COVER (PERCENT)                     ',
     9'MOIST CONVECTIVE CLOUD COVER (PERCENT)                     ',
     O'AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 1 (METERS)',
     1'AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 2 (METERS)',
     2'AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 3 (METERS)',
     3'AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 4 (METERS)'/
      CHARACTER*64 TITLEB(9)
      DATA TITLEB/
C****                                                             84-92
     4'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 1 (DEG WEST LONG) ',
     5'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 2 (DEG WEST LONG) ',
     6'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 3 (DEG WEST LONG) ',
     7'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 4 (DEG WEST LONG) ',
     8'NORTH. TRANS. OF SENSIBLE HEAT BY EDDIES (10**14 WATTS/DSIGMA) ',
C    9'TOTAL NORTHWARD TRANSPORT OF SENSIBLE HEAT (10**15 WATTS/DSIGMA)'
C    *,
     9'NUMBER OF GRIDPOINTS INCLUDED IN AVERAGE             CP        ',
     O'VERT. TRANS. OF GEOPOTENTIAL ENERGY BY EDDIES (10**12 WATTS)   ',
C    1'TOTAL LARGE SCALE VERT. TRANS. OF GEOPOTEN. ENER. (10**14 WATTS)'
C    *,
     1'PRESSURE DIFFERENCES  (MB)                           CP        ',
     2'SUBGRID SCALE TEMPERATURE DEVIATION (HUNDREDTHS OF DEGREES KEL)'
     */
      CHARACTER*64 TITLEC(6)
      DATA TITLEC/
C****                                                             93-98
     3'DYNAMIC CONVERGENCE OF DRY STATIC ENERGY (10 WATTS/M**2/DSIGMA)',
     4'DIVERGENCE OF THE ELIASSEN-PALM FLUX (10**17 JOULES/UNIT SIGMA)',
     5'NORTHWARD ELIASSEN-PALM FLUX (10**17 JOULES/UNIT SIGMA)',
     6'VERTICAL ELIASSEN-PALM FLUX (10**17 JOULES) ',
     7'TOTAL DRYING BY MOIST CONVECTION (Q2)        10*14 WATTS/DSIG) ',
     8'TOTAL HEATING BY MOIST CONVECTION (Q1)      (10*14 WATTS/DSIG) '/
      CHARACTER*64 TITLED(8)
      DATA TITLED/
C****                                                            99-106
     1'VERT. TRANS. OF GEOPOTENTIAL ENERGY BY EDDIES (10**12 WATTS)  CP'
     *,
     2'VERT. TRANS. OF DRY STATIC ENERGY BY EDDIES (10**12 WATTS)    CP'
     *,
     3'TOTAL LGE SCALE VERT. TRANS. OF DRY STAT. ENER. (10**14 WATTS)CP'
     *,
     4'VERTICAL TRANSPORT OF LATENT HEAT BY EDDIES (10**12 WATTS)    CP'
     *,
     5'TOTAL LARGE SCALE VERT. TRANS. OF LATENT HEAT (10**13 WATTS)  CP'
     *,
     6'VERTICAL TRANSPORT OF STATIC ENERGY BY EDDIES (10**13 WATTS)  CP'
     *,
     7'TOTAL LARGE SCALE VERT. TRANS. OF STATIC ENERGY (10**14 WATTS)CP'
     *,
     8'TRANSFORMED STREAM FUNCTION  (10**9 KG/SEC)                   CP'
     *     /
      CHARACTER*64 TITLEE(8)
      DATA TITLEE/
C****                                                           107-114
     7'VERT. TRANSPORT OF POTENTIAL VORTICITY (10**4 KG-DEG K/MB/S/S)CP'
     *,
     8'VERT. TRANS. OF POT. VORT. BY EDDIES (10**4 KG-DEG K/MB/S/S)  CP'
     *,
     9'                                                                '
     *,
     A'TOTAL LGE SCALE VERT. TRANS. OF KINETIC ENERGY (10**11 WATTS) CP'
     *,
     5'VERT. TRANS. OF ANG. MOMENTUM BY EDDIES (10**16 JOULES)       CP'
     *,
     6'TOTAL LGE SCALE VERT. TRANS. OF ANG. MOMENTUM (10**18 JOULES) CP'
     *,
     7'VERTICAL ELIASSEN-PALM FLUX (10**11 JOULES/METER)             CP'
     *,
     E'                                                                '
     *     /
      CHARACTER*64 TITLEF(14)
      DATA TITLEF/
C****                                                            115-128
     5'STAND. N. TRANS.OF D. STATIC ENERGY, WAVE #1(10**14 WATTS/DSIG)',
     6'STAND. N. TRANS.OF D. STATIC ENERGY, WAVE #2(10**14 WATTS/DSIG)',
     7'STAND. N. TRANS.OF D. STATIC ENERGY, WAVE #3(10**14 WATTS/DSIG)',
     8'STAND. N. TRANS.OF D. STATIC ENERGY, WAVE #4(10**14 WATTS/DSIG)',
     9'STAND. N. TRANS.OF D. STATIC ENERGY, WAVE #5(10**14 WATTS/DSIG)',
     O'TRNSNT N. TRANS.OF D. STATIC ENERGY, WAVE #1(10**14 WATTS/DSIG)',
     1'TRNSNT N. TRANS.OF D. STATIC ENERGY, WAVE #2(10**14 WATTS/DSIG)',
     2'TRNSNT N. TRANS.OF D. STATIC ENERGY, WAVE #3(10**14 WATTS/DSIG)',
     3'TRNSNT N. TRANS.OF D. STATIC ENERGY, WAVE #4(10**14 WATTS/DSIG)',
     4'TRNSNT N. TRANS.OF D. STATIC ENERGY, WAVE #5(10**14 WATTS/DSIG)',
     5'TRNSNT N. TRANS.OF D. STATIC ENERGY, WAVE #6(10**14 WATTS/DSIG)',
     6'TRNSNT N. TRANS.OF D. STATIC ENERGY, WAVE #7(10**14 WATTS/DSIG)',
     7'TRNSNT N. TRANS.OF D. STATIC ENERGY, WAVE #8(10**14 WATTS/DSIG)',
     8'TRNSNT N. TRANS.OF D. STATIC ENERGY, WAVE #9(10**14 WATTS/DSIG)'/
      CHARACTER*64 TITLEG(14)
      DATA TITLEG/
     9'STAND. N. TRANS. OF LATENT HEAT, WAVE #1 (10**13 WATTS/DSIG)   ',
     O'STAND. N. TRANS. OF LATENT HEAT, WAVE #2 (10**13 WATTS/DSIG)   ',
     1'STAND. N. TRANS. OF LATENT HEAT, WAVE #3 (10**13 WATTS/DSIG)   ',
     2'STAND. N. TRANS. OF LATENT HEAT, WAVE #4 (10**13 WATTS/DSIG)   ',
     3'STAND. N. TRANS. OF LATENT HEAT, WAVE #5 (10**13 WATTS/DSIG)   ',
     4'TRANSIENT N. TRANS. OF LATENT HEAT, WAVE #1 (10**13 WATTS/DSIG)',
     5'TRANSIENT N. TRANS. OF LATENT HEAT, WAVE #2 (10**13 WATTS/DSIG)',
     6'TRANSIENT N. TRANS. OF LATENT HEAT, WAVE #3 (10**13 WATTS/DSIG)',
     7'TRANSIENT N. TRANS. OF LATENT HEAT, WAVE #4 (10**13 WATTS/DSIG)',
     8'TRANSIENT N. TRANS. OF LATENT HEAT, WAVE #5 (10**13 WATTS/DSIG)',
     9'TRANSIENT N. TRANS. OF LATENT HEAT, WAVE #6 (10**13 WATTS/DSIG)',
     O'TRANSIENT N. TRANS. OF LATENT HEAT, WAVE #7 (10**13 WATTS/DSIG)',
     1'TRANSIENT N. TRANS. OF LATENT HEAT, WAVE #8 (10**13 WATTS/DSIG)',
     2'TRANSIENT N. TRANS. OF LATENT HEAT, WAVE #9 (10**13 WATTS/DSIG)'/
      CHARACTER*64 TITLEH(14)
      DATA TITLEH/
C****                                                            143-156
     3'STAND. N. TRANS. OF ANG. MOMENT., WAVE #1 (10**18 JOULES/DSIG)',
     4'STAND. N. TRANS. OF ANG. MOMENT., WAVE #2 (10**18 JOULES/DSIG)',
     5'STAND. N. TRANS. OF ANG. MOMENT., WAVE #3 (10**18 JOULES/DSIG)',
     6'STAND. N. TRANS. OF ANG. MOMENT., WAVE #4 (10**18 JOULES/DSIG)',
     7'STAND. N. TRANS. OF ANG. MOMENT., WAVE #5 (10**18 JOULES/DSIG)',
     8'TRNSNT N. TRANS. OF ANG. MOMENT., WAVE #1 (10**18 JOULES/DS)  ',
     9'TRNSNT N. TRANS. OF ANG. MOMENT., WAVE #2 (10**18 JOULES/DS)  ',
     O'TRNSNT N. TRANS. OF ANG. MOMENT., WAVE #3 (10**18 JOULES/DS)  ',
     1'TRNSNT N. TRANS. OF ANG. MOMENT., WAVE #4 (10**18 JOULES/DS)  ',
     2'TRNSNT N. TRANS. OF ANG. MOMENT., WAVE #5 (10**18 JOULES/DS)  ',
     3'TRNSNT N. TRANS. OF ANG. MOMENT., WAVE #6 (10**18 JOULES/DS)  ',
     4'TRNSNT N. TRANS. OF ANG. MOMENT., WAVE #7 (10**18 JOULES/DS)  ',
     5'TRNSNT N. TRANS. OF ANG. MOMENT., WAVE #8 (10**18 JOULES/DS)  ',
     6'TRNSNT N. TRANS. OF ANG. MOMENT., WAVE #9 (10**18 JOULES/DS)  '/
      END BLOCK DATA BDJL

      SUBROUTINE DIAGJL
      USE CONSTANT, only : grav,rgas,kapa,sday,twopi,sha
      USE E001M12_COM, only : im,jm,lm,FIM,
     &     BYDSIG,BYIM,DSIG,DT,IDACC,IMH,NCNDS,
     &     PSF,PTOP,SIG,SIGE,TOFDAY,TOFDY0
      USE GEOM, only :
     &     AREAG,BYDXYP,COSP,COSV,DLON,DXV,DXYP,DXYV,FCOR,RADIUS,WTJ
      USE DAGPCOM
      USE RADNCB, only : LM_REQ
      USE DAGCOM, only : ajl,apj,asjl,kdiag,aij
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(JM) ::
     &     BYP,BYPV,BYDAPO,BYPDA,DXCOSV,DACOSV,DXYPPO,ONES,ONESPO
      DOUBLE PRECISION, DIMENSION(JM,LM) :: AX,BX,CX,DX
      DOUBLE PRECISION, DIMENSION(JM,LM_REQ) :: ARQX
      DOUBLE PRECISION, DIMENSION(LM_REQ) :: BYDPS,BYPKS
      DOUBLE PRECISION, DIMENSION(0:IMH) :: AN,BN
      DOUBLE PRECISION, DIMENSION(JM,8,4) :: AMPLTD,PHASE

      INTEGER :: LINECT,JMHALF,INC,IHOUR0,IHOUR
      COMMON/DJLCOM/LINECT,JMHALF,INC,IHOUR0,IHOUR

      DOUBLE PRECISION, DIMENSION(LM) :: PKM,BYD2SG
      DOUBLE PRECISION, DIMENSION(LM+LM_REQ) :: PL

      DOUBLE PRECISION, DIMENSION(7), PARAMETER ::
     &     PMB=(/999.9,850.,700.,500.,300.,100.,30./)
      DOUBLE PRECISION, PARAMETER :: P1000=1000.

      INTEGER :: J,K,KM,L,LDN,LS,LUP,N

      DOUBLE PRECISION ::
     &     BDN,BUP,BY100G,BYIACN,BYIADA,BYIARD,BYIMDA,
     &     DAM4,DTCNDS,DXCVN,DXCVS,ELOFIM,FIMDA,P1000K,
     &     SCALE,SCALE2,SCALES,SCALEV,XWON

C**** INITIALIZE CERTAIN QUANTITIES
      XWON=TWOPI/(DLON*FIM)
      INC=1+(JM-1)/24
      JMHALF=JM/2
c      BYIM=1./FIM
      BY100G=.01/GRAV
c      SHA=RGAS/KAPA
      DTCNDS=NCNDS*DT
      P1000K=P1000**KAPA
      KM=0
      DO 5 K=1,7
      IF (PMTOP.GT.PMB(K)) GO TO 6
    5 KM=KM+1
    6 ELOFIM=.5*TWOPI-TWOPI/FIM
      DO 20 L=1,LM
      LUP=L+1
      LDN=L-1
      IF (L.EQ.LM) LUP=LM
      IF (L.EQ.1) LDN=1
      BYD2SG(L)=1./(SIG(LUP)-SIG(LDN))
   20 CONTINUE
      BYDPS(1)=1./(.5*PMTOP)
      BYDPS(2)=1./(.3*PMTOP)
      BYDPS(3)=1./(.2*PMTOP)
      BYPKS(1)=1./(.75*PMTOP)**KAPA
      BYPKS(2)=1./(.35*PMTOP)**KAPA
      BYPKS(3)=1./(.1*PMTOP)**KAPA
      DO 40 J=1,JM
      DXYPPO(J)=DXYP(J)
c      BYDXYP(J)=1./DXYP(J)
      BYDAPO(J)=BYDXYP(J)
      ONES(J)=1.
      ONESPO(J)=1.
   40 CONTINUE
      DXYPPO(JM)=DXYP(JM)*FIM
      DXYPPO(1)=DXYP(1)*FIM
      BYDAPO(1)=BYDAPO(1)*FIM
      BYDAPO(JM)=BYDAPO(JM)*FIM
      ONESPO(1)=FIM
      ONESPO(JM)=FIM
      DO 50 J=2,JM
      DXCOSV(J)=DXV(J)*COSV(J)
      DACOSV(J)=DXYV(J)*COSV(J)
   50 CONTINUE
      IHOUR0=TOFDY0+.5
      IHOUR=TOFDAY+.5
      LINECT=65
      BYIACN=1./(IDACC(1)+1.D-20)
      BYIARD=1./(IDACC(2)+1.D-20)
      BYIADA=1./(IDACC(4)+1.D-20)
      BYIMDA=BYIADA*BYIM
      FIMDA=IDACC(4)*FIM
      DO 120 J=1,JM
      BYPDA(J)=1./(APJ(J,1)*DXYP(J)+1.D-20)
      BYP(J)=1./(APJ(J,1)+1.D-20)
  120 BYPV(J)=1./(APJ(J,2)+1.D-20)
C****
C**** PROGNOSTIC QUANTITIES
C****
C**** TEMPERATURE, HEIGHT, SPECIFIC HUMIDITY, AND RELATIVE HUMIDITY
C     CALL JLMAPS (1,PLM,AJL,ONES,BYP,ONES,LM,2,1,
C    *  ASJL,BYIMDA,ONESPO,ONES)
C     SCALES=BYIMDA*BY100G
C     CALL JLMAPS (2,PLM,AJL(1,1,2),BY100G,BYP,ONES,LM,2,1,
C    *  ASJL(1,1,2),SCALES,ONESPO,ONES)
C     SCALE=1.D5
C     CALL JLMAP (3,PLM,AJL(1,1,3),SCALE,BYP,ONES,LM,2,1)
C     SCALE=100.
C     CALL JLMAP (4,PLM,AJL(1,1,18),SCALE,BYP,ONES,LM,2,1)
C**** U WIND, V WIND, AND STREAM FUNCTION
C     SCALE=10.
C     CALL JLMAP (5,PLM,AJL(1,1,4),SCALE,BYPV,ONES,LM,2,2)
C     SCALE=100.
C     CALL JLMAP (6,PLM,AJL(1,1,5),SCALE,BYPV,ONES,LM,2,2)
C     DO 220 J=2,JM
C     AX(J,1)=AJL(J,1,5)*DSIG(1)
C     BX(J,1)=(AJL(J,1,5)-.5*AJL(J,1,47)*FIM)*DSIG(1)
C     DO 220 L=2,LM
C     BX(J,L)=BX(J,L-1)+(AJL(J,L,5)-.5*AJL(J,L,47)*FIM)*DSIG(L)
C 220 AX(J,L)=AX(J,L-1)+AJL(J,L,5)*DSIG(L)
C     SCALE=25.D-9*BYIADA/GRAV
C     CALL JLMAP (7,PLE,AX,SCALE,DXV,ONES,LM,2,2)
C     CALL JLMAP (57,PLE,BX,SCALE,DXV,ONES,LM,2,2)
C**** VERTICAL VELOCITY AND MASS FLUX MOIST CONVECTION
C     SCALE=-1.D5*BYIMDA
C     CALL JLMAP (8,PLE,AJL(1,1,6),SCALE,BYDAPO,ONES,LM-1,2,1)
      SCALE=100.D-9*XWON*BYIACN/(GRAV*DTCNDS)
      CALL JLMAP (10,PLE,AJL(1,1,8),SCALE,DXYPPO,ONES,LM-1,1,1)
C****
C**** RADIATION, CONDENSATION AND CONVECTION
C****
C**** SOLAR AND THERMAL RADIATION HEATING
      SCALE=100.D-2*GRAV*SDAY*IDACC(4)*BYIARD/SHA
      SCALES=100.D-2*GRAV*SDAY*BYIM*BYIARD/SHA
      CALL JLMAPS (11,PLM,AJL(1,1,9),SCALE,BYP,BYDSIG,LM,2,1,
     *  ASJL(1,1,3),SCALES,ONESPO,BYDPS)
      SCALES=-SCALES
      SCALE=-SCALE
      CALL JLMAPS (12,PLM,AJL(1,1,10),SCALE,BYP,BYDSIG,LM,2,1,
     *  ASJL(1,1,4),SCALES,ONESPO,BYDPS)
      DO 250 J=1,JM
      DO 240 LS=1,LM_REQ
  240 ARQX(J,LS)=ASJL(J,LS,3)+ASJL(J,LS,4)
      DO 250 L=1,LM
  250 AX(J,L)=AJL(J,L,9)+AJL(J,L,10)
      SCALE=-1.D-13*XWON*BYIARD
      SCALES=SCALE*(PSF-PMTOP)
      CALL JLMAPS (13,PLM,AX,SCALE,DXYPPO,BYDSIG,LM,1,1,
     *  ARQX,SCALES,DXYPPO,BYDPS)
C**** TOTAL, SUPER SATURATION, AND CONVECTIVE CLOUD COVER
      SCALE=100.*BYIARD*BYIM
      CALL JLMAP (77,PLM,AJL(1,1,19),SCALE,ONESPO,ONES,LM,2,1)
      CALL JLMAP (78,PLM,AJL(1,1,28),SCALE,ONESPO,ONES,LM,2,1)
      CALL JLMAP (79,PLM,AJL(1,1,29),SCALE,ONESPO,ONES,LM,2,1)
C**** SUBGRID SCALE TEMPERATURE DEVIATION
      SCALE=1.D2*SQRT(2.)*BYIACN
      CALL JLMAP (92,PLM,AJL(1,1,54),SCALE,ONES,ONES,LM,2,1)
C**** HEATING BY LARGE SCALE CONDENSATION AND DRY AND MOIST CONVECTION
      SCALE=100.D-13*XWON*SHA*BYIACN/(GRAV*DTCNDS)
      CALL JLMAP (14,PLM,AJL(1,1,11),SCALE,DXYPPO,BYDSIG,LM,1,1)
      CALL JLMAP (15,PLM,AJL(1,1,12),SCALE,DXYPPO,ONES,LM,1,1)
      SCALE=0.1*SCALE
      CALL JLMAP (16,PLM,AJL(1,1,55),SCALE,DXYPPO,BYDSIG,LM,1,1)
      CALL JLMAP (98,PLM,AJL(1,1,56),SCALE,DXYPPO,ONES,LM,1,1)
      CALL JLMAP (97,PLM,AJL(1,1,57),SCALE,DXYPPO,ONES,LM,1,1)
C****
C**** ENERGY
C****
C**** AVAILABLE POTENTIAL ENERGY
      SCALE=50.D-5*RGAS*BYIMDA/GRAV
      CALL JLMAP (20,PLM,AJL(1,1,16),SCALE,ONES,ONES,LM,2,1)
C****
C**** NORTHWARD TRANSPORTS
C****
C**** NOR. TRANSPORT OF QUASI-GEOSTROPHIC POT. VORTICITY BY EDDIES
      DO 366 L=1,LM
      CX(1,L)=0.
      CX(2,L)=DXCOSV(2)*(AJL(2,L,49)-AJL(2,L,48))+.25*FIM*FCOR(2)*
     *  COSP(2)*(AJL(2,L,47)+AJL(3,L,47))
      DO 364 J=3,JM-1
      DAM4=DXCOSV(J)*(AJL(J,L,49)-AJL(J,L,48))
      CX(J,L)=DAM4+.25*FIM*FCOR(J)*COSP(J)*(AJL(J,L,47)+AJL(J-1,L,47))
      CX(J-1,L)=CX(J-1,L)-DAM4
  364 CONTINUE
      CX(JM-1,L)=CX(JM-1,L)-DXCOSV(JM)*(AJL(JM,L,49)-AJL(JM,L,48))
      CX(JM,L)=0.
  366 CONTINUE
      SCALE=25.D-18*XWON*BYIADA*RADIUS/GRAV * 0. ! AJL47 not done
      CALL JLMAP (54,PLM,CX,SCALE,ONES,ONES,LM,1,1)
C****
C**** VERTICAL TRANSPORTS
C****
C**** VERTICAL TRANSPORT OF ANGULAR MOMENTUM BY SMALL SCALE MOTIONS
      SCALE=100.D-18*XWON*RADIUS*BYIACN/(GRAV*DTCNDS)
      CALL JLMAP (44,PLM,AJL(1,1,38),SCALE,DACOSV,ONES,LM,1,2)
      CALL JLMAP (45,PLM,AJL(1,1,39),SCALE,DACOSV,BYDSIG,LM,1,2)
C     CALL JLMAP (46,PLM,AJL(1,1,40),SCALE,DACOSV,BYDSIG,LM,1,2)
C        IF (JM.NE.24) GO TO 500
C****
C**** MERIDIONAL LUNES
C****
C**** U, V AND W VELOCITY FOR EAST PACIFIC
      SCALE=.2E+1*BYIADA
      CALL JLMAP (47,PLM,AJL(1,1,41),SCALE,ONES,ONES,LM,2,2)
      CALL JLMAP (48,PLM,AJL(1,1,42),SCALE,ONES,ONES,LM,2,2)
      SCALE2=-1.D5*BYIADA*RGAS/(5.*GRAV)
      CALL JLMAP (49,PLE,AJL(1,1,43),SCALE2,BYDXYP,ONES,LM-1,2,1)
C**** U, V AND W VELOCITY FOR WEST PACIFIC
      CALL JLMAP (50,PLM,AJL(1,1,44),SCALE,ONES,ONES,LM,2,2)
      CALL JLMAP (51,PLM,AJL(1,1,45),SCALE,ONES,ONES,LM,2,2)
      CALL JLMAP (52,PLE,AJL(1,1,46),SCALE2,BYDXYP,ONES,LM-1,2,1)
  500 CONTINUE
C****
C**** ELIASSEN-PALM FLUX : NORTHWARD, VERTICAL, DIVERGENCE
C****
      SCALE=100.D-17*XWON*BYIADA*RADIUS/GRAV
      CALL JLMAP (95,PLM,AJL(1,1,37),SCALE,DXCOSV,ONES,LM,1,2)
      SCALEV=.125*SCALE
      CALL JLMAP (96,PLE,AJL(1,1,36),SCALEV,COSP,ONES,LM-1,1,1)
      DXCVS=DXCOSV(2)
      DO 540 J=2,JM-1
      BDN=0.
      DXCVN=DXCOSV(J+1)
      DO 530 L=1,LM
      BUP=AJL(J,L,36)*COSP(J)
      AX(J,L)=AJL(J+1,L,37)*DXCVN-AJL(J,L,37)*DXCVS+
     *   .125*(BUP-BDN)/DSIG(L)
  530 BDN=BUP
  540 DXCVS=DXCVN
      DO 550 L=1,LM
      AX(1,L)=0.
  550 AX(JM,L)=0.
      CALL JLMAP (94,PLM,AX,SCALE,ONES,ONES,LM,1,1)
C****
C**** FOURIER ANALYSIS OF GEOPOTENTIAL HEIGHTS FOR WAVE NUMBERS 1 TO 4,
C****   AMPLITUDE AND PHASE
C****
            LINECT=63
      DO 620 K=1,KM
      DO 610 N=1,4
      AMPLTD(1,K,N)=0.
      AMPLTD(JM,K,N)=0.
      PHASE(1,K,N)=0.
  610 PHASE(JM,K,N)=0.
      DO 620 J=2,JM-1
      CALL FFT (AIJ(1,J,8+K),AN,BN)
      DO 620 N=1,4
      AMPLTD(J,K,N)=SQRT(AN(N)*AN(N)+BN(N)*BN(N))
      PHASE(J,K,N)=(ATAN2(BN(N),AN(N))-TWOPI)/N+ELOFIM
      IF (PHASE(J,K,N).LE.-.5*TWOPI) PHASE(J,K,N)=PHASE(J,K,N)+TWOPI
      PHASE(J,K,N)=-PHASE(J,K,N)
  620 CONTINUE
      SCALE=BYIADA/GRAV
      DO 630 N=1,4
  630 CALL JLMAP (N+79,PMB,AMPLTD(1,1,N),SCALE,ONES,ONES,KM,2,1)
      SCALE=360./TWOPI
      DO 640 N=1,4
  640 CALL JLMAP (N+83,PMB,PHASE(1,1,N),SCALE,ONES,ONES,KM,2,1)
      RETURN
      END SUBROUTINE DIAGJL

      SUBROUTINE JLMAP (NT,PL,AX,SCALE,SCALEJ,SCALEL,LMAX,JWT,J1)
C****
C**** THIS SUBROUTINE PRODUCES LAYER BY LATITUDE TABLES ON THE LINE
C**** PRINTER.  THE INTERIOR NUMBERS OF THE TABLE ARE CALCULATED AS
C****               AX * SCALE * SCALEJ * SCALEL.
C**** WHEN JWT=1, THE INSIDE NUMBERS ARE NOT AREA WEIGHTED AND THE
C****    HEMISPHERIC AND GLOBAL NUMBERS ARE SUMMATIONS.
C**** WHEN JWT=2, ALL NUMBERS ARE PER UNIT AREA.
C**** J1 INDICATES PRIMARY OR SECONDARY GRID.
C**** THE BOTTOM LINE IS CALCULATED AS THE SUMMATION OF DSIG TIMES THE
C**** NUMBERS ABOVE (POSSIBLY MULTIPLIED BY A FACTOR OF 10)
C****
      USE E001M12_COM, only : jm,lm,
     &     DSIG,JDATE,JDATE0,JMNTH0,JMONTH,JYEAR,JYEAR0,SIGE,XLABEL
      USE GEOM, only : JLAT,WTJ
      IMPLICIT NONE

      INTEGER, DIMENSION(JM) :: MLAT
      DOUBLE PRECISION, DIMENSION(JM) :: FLAT,ASUM
      DOUBLE PRECISION, DIMENSION(2) :: FHEM,HSUM

      INTEGER :: LINECT,JMHALF,INC,IHOUR0,IHOUR
      COMMON/DJLCOM/LINECT,JMHALF,INC,IHOUR0,IHOUR

      INTEGER :: J1,JWT,LMAX,NT
      DOUBLE PRECISION :: SCALE,SCALER
      DOUBLE PRECISION, DIMENSION(JM,LM) :: AX
      DOUBLE PRECISION, DIMENSION(JM,3) :: ARQX ! pass 3 as arg ??
      DOUBLE PRECISION, DIMENSION(JM) :: SCALEJ,SCALJR
      DOUBLE PRECISION, DIMENSION(LM) :: SCALEL
      DOUBLE PRECISION, DIMENSION(3) :: SCALLR
      DOUBLE PRECISION, DIMENSION(LM+3) :: PL

      COMMON/DJLTTL/TITLE(1)
      CHARACTER*4 DASH,WORD(4),TITLE*64
      DATA DASH/'----'/,WORD/'SUM','MEAN',' ','.1*'/

      INTEGER :: IWORD,J,J0,JH,JHEMI,JY,JY0,K,L
      DOUBLE PRECISION :: FGLOB,FLATJ,GSUM,SDSIG,SUMFAC

C****
C**** PRODUCE A LATITUDE BY LAYER TABLE OF THE ARRAY A
C****
   10 LINECT=LINECT+LMAX+7
      IF (LINECT.LE.60) GO TO 20
      JY0=JYEAR0-1900
      JY=JYEAR-1900
      WRITE (6,907) XLABEL(1:108),JDATE0,JMNTH0,JY0,JDATE,JMONTH,JY
      LINECT=LMAX+8
   20 WRITE (6,901) TITLE(NT),(DASH,J=J1,JM,INC)
      WRITE (6,904) WORD(JWT),(JLAT(J,J1),J=JM,J1,-INC)
      WRITE (6,905) (DASH,J=J1,JM,INC)
      J0=J1-1
  100 SDSIG=1.-SIGE(LMAX+1)
      DO 110 J=1,JM
  110 ASUM(J)=0.
      HSUM(1)=0.
      HSUM(2)=0.
      GSUM=0.
      SUMFAC=1.
      IWORD=3
      IF (NT.NE.1.AND.NT.NE.6.AND.NT.NE.24.AND.NT.NE.26.AND.NT.NE.28
     *  .AND.NT.NE.33) GO TO 112
      SUMFAC=10.
      IWORD=4
  112 DO 140 L=LMAX,1,-1
      FGLOB=0.
      DO 130 JHEMI=1,2
      FHEM(JHEMI)=0.
      DO 120 JH=1,JMHALF
      J=(JHEMI-1)*(JMHALF-J0)+JH+J0
      FLAT(J)=AX(J,L)*SCALE*SCALEJ(J)*SCALEL(L)
      MLAT(J)=NINT(FLAT(J))
  115 ASUM(J)=ASUM(J)+FLAT(J)*DSIG(L)/SDSIG
  120 FHEM(JHEMI)=FHEM(JHEMI)+FLAT(J)*WTJ(J,JWT,J1)
  130 FGLOB=FGLOB+FHEM(JHEMI)/JWT
      WRITE (6,902) PL(L),FGLOB,FHEM(2),FHEM(1),(MLAT(J),J=JM,J1,-INC)
C        IF (NT.EQ.5) CALL KEYJLJ (L,FLAT)
         IF (NT.EQ.61) CALL KEYJLS (L,FLAT)
  136 HSUM(1)=HSUM(1)+FHEM(1)*SUMFAC*DSIG(L)/SDSIG
      HSUM(2)=HSUM(2)+FHEM(2)*SUMFAC*DSIG(L)/SDSIG
  140 GSUM=GSUM+FGLOB*SUMFAC*DSIG(L)/SDSIG
      WRITE (6,905) (DASH,J=J1,JM,INC)
      IF (NT.GE.80.AND.NT.LE.87) RETURN
      ASUM(JMHALF+1)=ASUM(JMHALF+1)/J1
      DO 150 J=J1,JM
  150 MLAT(J)=NINT(ASUM(J)*SUMFAC)
      WRITE (6,903) WORD(IWORD),GSUM,HSUM(2),HSUM(1),
     *  (MLAT(J),J=JM,J1,-INC)
C        IF (NT.EQ.1) CALL KEYJLT (GSUM,ASUM)
C        IF (NT.EQ.18.OR.NT.EQ.19) CALL KEYJLE (NT,HSUM,ASUM)
C        IF (NT.GE.22.AND.NT.LE.33) CALL KEYJLN (NT,ASUM,SUMFAC)
      RETURN
C****
      ENTRY JLMAPS (NT,PL,AX,SCALE,SCALEJ,SCALEL,LMAX,JWT,J1,
     *  ARQX,SCALER,SCALJR,SCALLR)
      LINECT=LINECT+LMAX+10
      IF (LINECT.LE.60) GO TO 200
      JY0=JYEAR0-1900
      JY=JYEAR-1900
      WRITE (6,907) XLABEL(1:108),JDATE0,JMNTH0,JY0,JDATE,JMONTH,JY
      LINECT=LMAX+11
  200 J0=J1-1
C**** PRODUCE UPPER STRATOSPHERE NUMBERS FIRST
      WRITE (6,901) TITLE(NT),(DASH,J=J1,JM,INC)
      WRITE (6,904) WORD(JWT),(JLAT(J,J1),J=JM,J1,-INC)
      WRITE (6,905) (DASH,J=J1,JM,INC)
      DO 230 L=3,1,-1
      FGLOB=0.
      DO 220 JHEMI=1,2
      FHEM(JHEMI)=0.
      DO 210 JH=1,JMHALF
      J=(JHEMI-1)*(JMHALF-J0)+JH-J0
      FLATJ=ARQX(J,L)*SCALER*SCALJR(J)*SCALLR(L)
      MLAT(J)=NINT(FLATJ)
  210 FHEM(JHEMI)=FHEM(JHEMI)+FLATJ*WTJ(J,JWT,J1)
  220 FGLOB=FGLOB+FHEM(JHEMI)/JWT
  230 WRITE (6,902) PL(L+LM),FGLOB,FHEM(2),FHEM(1),
     *  (MLAT(J),J=JM,J1,-INC)
      GO TO 100
  901 FORMAT ('0',30X,A64/1X,30('-'),24A4)
  902 FORMAT (F6.1,3F8.1,1X,24I4)
  903 FORMAT (A6,3F8.1,1X,24I4)
  904 FORMAT (' P(MB) ',A4,' G      NH      SH  ',24I4)
  905 FORMAT (1X,30('-'),24A4)
  907 FORMAT ('1',A,I4,1X,A3,I3,' TO',I3,1X,A3,I3)
      END SUBROUTINE JLMAP

      BLOCK DATA BDIL
C****
C**** TITLES FOR SUBROUTINE DIAGIL
C****
      COMMON/DILTTL/TITLE1,TITLE2
      CHARACTER*64 TITLE1(8)
      DATA TITLE1/
C****                                                              1-8
     1'ZONAL WIND (U COMPONENT) AROUND +/- 5 DEG (METERS/SECOND)',
     2'MERIDIONAL WIND (V COMPONENT) AROUND +/- 5 DEG (METERS/SECOND)',
     3'VERTICAL VELOCITY AROUND +/- 5 DEG (10**-4 METERS/SECOND)',
     4'TEMPERATURE AROUND +/- 5 DEG (DEGREES CENTIGRADE)',
     5'RELATIVE HUMIDITY AROUND +/- 5 DEG (PERCENT) ',
     6'MOIST CONVECTIVE HEATING AROUND +/- 5 DEG (10**13 WATTS/DSIG)',
     7'TOTAL RADIATIVE COOLING AROUND +/- 5 DEG (10**13 WATTS/DSIG)',
     8'                                                              '/
      CHARACTER*64 TITLE2(8)
      DATA TITLE2/
C****                                                              9-16
     9'VERTICAL VELOCITY AT 50 N (10**-4 METERS/SECOND)         ',
     O'TEMPERATURE AT 50 N (DEGREES CENTIGRADE)                 ',
     1'TOTAL RADIATIVE COOLING AT 50 N (10**13 WATTS/UNIT SIGMA)',
     2'ZONAL WIND AT 50 N (METERS/SECOND)                       ',
     3'VERTICAL VELOCITY AT 70 N (10**-4 METERS/SECOND)         ',
     4'TEMPERATURE AT 70 N (DEGREES CENTIGRADE)                 ',
     5'TOTAL RADIATIVE COOLING AT 70 N (10**13 WATTS/UNIT SIGMA)',
     6'ZONAL WIND AT 70 N (METERS/SECOND)                       '/
      END BLOCK DATA BDIL

      SUBROUTINE DIAGIL
      USE CONSTANT, only : grav,rgas,sha
      USE E001M12_COM, only : im,jm,lm,fim,
     &     BYDSIG,DSIG,DT,IDACC,JEQ,NCNDS,
     &     PSF,PTOP,SIG,SIGE,TOFDAY,TOFDY0
      USE GEOM, only : AREAG,DXYP,DXYV
      USE RADNCB, only : LM_REQ
      USE DAGPCOM
      USE DAGCOM, only : ail
      IMPLICIT NONE

      INTEGER :: LINECT,JMHALF,INC,IHOUR0,IHOUR
      COMMON/DILCOM/LINECT,JMHALF,INC,IHOUR0,IHOUR

      CHARACTER TITLE*64
      COMMON/DILTTL/TITLE(1)

      DOUBLE PRECISION, DIMENSION(LM) :: ONES
      DOUBLE PRECISION, DIMENSION(LM+LM_REQ) :: PL

      INTEGER :: J,J50N,J70N,L,K
      DOUBLE PRECISION :: BYIACN,BYIADA,BYIARD,DTCNDS,SCALE

      INTEGER, PARAMETER :: K_PIL=16
      INTEGER, DIMENSION(K_PIL) ::
     *  KNDEX=(/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/),
     *  CASEIL=(/1,1,2,3,3, 4,5,0,2,3, 0,1,2,3,0, 1/)
c    *  CASEIL=(/1,1,2,3,3, 4,5,0,2,3, 5,1,2,3,5, 1/)
      DOUBLE PRECISION SIL(K_PIL)
C**** INITIALIZE CERTAIN QUANTITIES
      INC=1+(JM-1)/24
      JMHALF=JM/2
      J50N=(50.+90.)*(JM-1.)/180.+1.5
      J70N=(70.+90.)*(JM-1.)/180.+1.5
      DTCNDS=NCNDS*DT
      DO 20 L=1,LM
      ONES(L)=1.
   20 CONTINUE
      IHOUR0=TOFDY0+.5
      IHOUR=TOFDAY+.5
      LINECT=65
      BYIACN=1./(IDACC(1)+1.D-20)
      BYIARD=1./(IDACC(2)+1.D-20)
      BYIADA=1./(IDACC(4)+1.D-20)
      SIL( 1)=BYIADA/3.
      SIL( 2)=BYIADA/3.
      SIL( 3)=
     *  -1.D4*BYIADA*RGAS/GRAV/(DXYP(JEQ)+DXYP(JEQ-1)+DXYP(JEQ-2))
      SIL( 4)=BYIADA/3.
      SIL( 5)=1.D2*BYIADA/3.
      SIL( 6)=100.D-13*SHA*BYIACN/(GRAV*DTCNDS)
      SIL( 7)=-1.D-13*BYIARD
      SIL( 9) =-1.D4*BYIADA*RGAS/(GRAV* DXYP(J50N))
      SIL(10)=BYIADA
      SIL(11)=-1.D-13*BYIARD
      SIL(12)=BYIADA/2.
      SIL(13) =-1.D4*BYIADA*RGAS/(GRAV* DXYP(J70N))
      SIL(14)=BYIADA
      SIL(15)=-1.D-13*BYIARD
      SIL(16)=BYIADA/2.

      DO K=1,K_PIL
      SELECT CASE (CASEIL(K))
      CASE DEFAULT
      CASE (1)  ! Centered in L; secondary grid; hor. mean; vert. sum
        CALL ILMAP(TITLE(K),PLM,AIL(1,1,KNDEX(K)),SIL(K),ONES,LM,2,2)
      CASE (2)  ! Vertical edges; primary grid; hor. mean; vert. sum
        CALL ILMAP(TITLE(K),PLE,AIL(1,1,KNDEX(K)),SIL(K),ONES,LM-1,2,1)
      CASE (3)  ! Centered in L; primary grid; hor. mean; vert. sum
        CALL ILMAP(TITLE(K),PLM,AIL(1,1,KNDEX(K)),SIL(K),ONES,LM,2,1)
      CASE (4)  ! Centered in L; primary grid; hor. sum; vert. sum
        CALL ILMAP(TITLE(K),PLM,AIL(1,1,KNDEX(K)),SIL(K),ONES,LM,1,1)
      CASE (5)  ! Centered in L; primary grid; hor. sum; vert. mean
        CALL ILMAP(TITLE(K),PLM,AIL(1,1,KNDEX(K)),SIL(K),BYDSIG,LM,1,1)
      END SELECT
      END DO

      RETURN
      END SUBROUTINE DIAGIL

      SUBROUTINE ILMAP (TITLE,PL,AX,SCALE,SCALEL,LMAX,JWT,ISHIFT)
      USE CONSTANT, only : twopi
      USE E001M12_COM, only : im,jm,lm,
     &     DSIG,JDATE,JDATE0,JMNTH0,JMONTH,JYEAR,JYEAR0,SIGE,XLABEL
     &    ,Q_GISS,Q_HDF,Q_PRT,Q_NETCDF
      USE GEOM, only : DLON,LON,JLAT
      IMPLICIT NONE

      DOUBLE PRECISION :: XIL(IM,LM),ZONAL(LM) ! used for post-proc
      CHARACTER XLB*80,CWORD*8
      CHARACTER*16, PARAMETER :: CBLANK=' ',
     &  CLAT='LONGITUDE',CPRES='PRESSURE (MB)'

      DOUBLE PRECISION :: FGLOB,FLON,GSUM,SDSIG, SCALE
      INTEGER :: I,JY,JY0,K,L, LMAX,JWT,ISHIFT

      INTEGER, DIMENSION(IM) :: MLON
      DOUBLE PRECISION, DIMENSION(IM) :: ASUM

      INTEGER :: LINECT,JMHALF,INC,IHOUR0,IHOUR
      COMMON/DILCOM/LINECT,JMHALF,INC,IHOUR0,IHOUR

      DOUBLE PRECISION, DIMENSION(IM,LM) :: AX
      DOUBLE PRECISION, DIMENSION(LM) :: PL,SCALEL

      CHARACTER*64 TITLE
      CHARACTER*4, PARAMETER :: DASH='----'
      CHARACTER*4, DIMENSION(2), PARAMETER :: WORD=(/'SUM ','MEAN'/)
C****
C**** PRODUCE A LONGITUDE BY LAYER TABLE OF THE ARRAY A
C****
C**** ISHIFT: When=2, print longitude indices off center (U-grid)
      LINECT=LINECT+LMAX+7
      IF (LINECT.GT.60) THEN
        JY0=JYEAR0-1900
        JY=JYEAR-1900
        WRITE (6,907) XLABEL(1:108),JDATE0,JMNTH0,JY0,JDATE,JMONTH,JY
        LINECT=LMAX+8
      END IF
      SDSIG=1.-SIGE(LMAX+1)
      WRITE (6,901) TITLE,(DASH,I=1,IM,INC)
      IF (ISHIFT.EQ.1) WRITE (6,904) WORD(JWT),(I,I=1,IM,INC)
      IF (ISHIFT.EQ.2) WRITE (6,906) WORD(JWT),(I,I=1,IM,INC)
      WRITE (6,905) (DASH,I=1,IM,INC)
      ASUM(:)=0.
      GSUM=0.
      DO 130 L=LMAX,1,-1
      FGLOB=0.
      DO 120 I=1,IM
      FLON=AX(I,L)*SCALE*SCALEL(L)
         XIL(I,L)=FLON
      MLON(I)=NINT(FLON)
      ASUM(I)=ASUM(I)+FLON*DSIG(L)/SDSIG
      FGLOB=FGLOB+FLON
  120 CONTINUE
      FGLOB=FGLOB/IM
      IF (JWT.EQ.1) FGLOB=FGLOB*TWOPI/DLON
         ZONAL(L)=FGLOB
      WRITE (6,902) PL(L),FGLOB,(MLON(I),I=1,IM,INC)
      GSUM=GSUM+FGLOB*DSIG(L)/SDSIG
  130 CONTINUE
      DO 140 I=1,IM
      MLON(I)=NINT(ASUM(I))
  140 CONTINUE
      WRITE (6,905) (DASH,I=1,IM,INC)
      WRITE (6,903) GSUM,(MLON(I),I=1,IM,INC)
C**** Output for post-processing
         CWORD=WORD(JWT)  ! pads out to 8 characters
         XLB=TITLE
         WRITE(XLB(65:),'(1X,A3,I5)') JMNTH0,JYEAR0
c     IF (Q_GISS) THEN
c        WRITE (85) XLB,IM,LMAX,1,1,
c    *     ((XIL(I,L),I=1,IM),L=1,LMAX),
c    *     (LON(I,ISHIFT),I=1,IM),(PL(L),L=1,LMAX),0.,0.,
c    *     CLAT,CPRES,CBLANK,CBLANK,CWORD,
c    *     (ASUM(I),I=1,IM),GSUM,(ZONAL(L),L=1,LMAX)
c     END IF
      IF (Q_GISS) CALL POUT_GISS
     *  ('IL',XLB,IM,LMAX,1,1,XIL,LON(1,ISHIFT),PL,0.,0.,
     *     CLAT,CPRES,CBLANK,CBLANK,CWORD,ASUM,GSUM,ZONAL)
      IF (Q_HDF) CALL POUT_HDF
     *  ('IL',XLB,IM,LMAX,1,1,XIL,LON(1,ISHIFT),PL,0.,0.,
     *     CLAT,CPRES,CBLANK,CBLANK,CWORD,ASUM,GSUM,ZONAL)
      IF (Q_NETCDF) CALL POUT_NETCDF
     *  ('IL',XLB,IM,LMAX,1,1,XIL,LON(1,ISHIFT),PL,0.,0.,
     *     CLAT,CPRES,CBLANK,CBLANK,CWORD,ASUM,GSUM,ZONAL)
      RETURN
  901 FORMAT ('0',30X,A64/1X,14('-'),36A3)
  902 FORMAT (F6.1,F8.1,1X,36I3)
  903 FORMAT (F14.1,1X,36I3)
  904 FORMAT (' P(MB)',4X,A4,1X,36I3)  ! U-grid (i.e., centers)
  905 FORMAT (1X,14('-'),36A3)
  906 FORMAT (' P(MB)',4X,A4,36I3)     ! V-grid (i.e., edges)
  907 FORMAT ('1',A,I4,1X,A3,I3,' TO',I3,1X,A3,I3)
      END SUBROUTINE ILMAP

      BLOCK DATA BDWP
C****
C**** TITLES FOR SUBROUTINE DIAG7
C****
      COMMON/D7COM/TITLE
      CHARACTER*64 TITLE(12)
      DATA TITLE/
     1'WAVE POWER FOR U NEAR 850 MB AND EQUATOR (DAY*(M/S)**2)       ',
     2'WAVE POWER FOR V NEAR 850 MB AND EQUATOR (DAY*(M/S)**2)       ',
     3'WAVE POWER FOR U NEAR 300 MB AND EQUATOR (10 DAY*(M/S)**2)    ',
     4'WAVE POWER FOR V NEAR 300 MB AND EQUATOR (DAY*(M/S)**2)       ',
     5'WAVE POWER FOR U NEAR 50 MB AND EQUATOR (10 DAY*(M/S)**2)     ',
     6'WAVE POWER FOR V NEAR 50 MB AND EQUATOR (DAY*(M/S)**2)        ',
     7'WAVE POWER FOR PHI AT 922 MB AND 50 DEG NORTH (10**3 DAY*M**2)',
     8'WAVE POWER FOR PHI AT 700 MB AND 50 DEG NORTH (10**3 DAY*M**2)',
     9'WAVE POWER FOR PHI AT 500 MB AND 50 DEG NORTH (10**3 DAY*M**2)',
     A'WAVE POWER FOR PHI AT 300 MB AND 50 DEG NORTH (10**3 DAY*M**2)',
     B'WAVE POWER FOR PHI AT 100 MB AND 50 DEG NORTH (10**4 DAY*M**2)',
     C'WAVE POWER FOR PHI AT 10 MB AND 50 DEG NORTH (10**4 DAY*M**2) '/
      END BLOCK DATA BDWP

      SUBROUTINE DIAG7A
C****
C**** THIS SUBROUTINE ACCUMULATES A TIME SEQUENCE FOR SELECTED
C**** QUANTITIES AND FROM THAT PRINTS A TABLE OF WAVE FREQUENCIES.
C****
      USE CONSTANT, only : grav
      USE E001M12_COM, only : im,imh,jm,lm,
     &     IDACC,JEQ,LS1,MDIAG,P,PSF,PTOP,PSFMPT,SIG,SIGE,U,V
      USE DAGCOM, only : wave,n12hrs_in_31day
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: PHI
      COMMON/WORK3/PHI

      DOUBLE PRECISION, DIMENSION(0:IMH) :: AN,BN

      INTEGER, PARAMETER :: KM=6,KQMAX=12
      INTEGER :: NMAX=9

      DOUBLE PRECISION, DIMENSION(IM,KM) :: HTRD
      INTEGER, DIMENSION(KM) :: JLKDEX

      DOUBLE PRECISION, DIMENSION(KM), PARAMETER ::
     &     PMB=(/922.,700.,500.,300.,100.,10./),
     &     GHT=(/500.,2600.,5100.,8500.,15400.,30000./)

      DOUBLE PRECISION ::
     &     PIJ50N,PL,PLE,PLM1,SLOPE

      INTEGER ::
     &     I,IDACC9,J50N,JLK,K,KQ,L,L300,L50,L850,LX,MINC,MNOW,N

      INTEGER :: IFIRST
      DATA IFIRST/1/

      IDACC9=IDACC(9)+1
      IDACC(9)=IDACC9
      IF (IDACC9.GT.N12HRS_IN_31DAY) RETURN
      NMAX=MIN0(9,IM/2)
      IF (IFIRST.NE.1) GO TO 100
      IFIRST=0
C      JEQ=1+JM/2
      J50N=(50.+90.)*(JM-1.)/180.+1.5
      L850=LM
      L300=LM
      L50=LM
      DO 10 L=2,LM
      LX=LM+1-L
      PLE=.25*(SIGE(LX)+2.*SIGE(LX+1)+SIGE(LX+2))*PSFMPT+PTOP
      IF (PLE.LT.850.) L850=LX
      IF (PLE.LT.300.) L300=LX
   10 IF (PLE.LT.50.) L50=LX
      WRITE (6,889) L850,L300,L50
  889 FORMAT (' LEVELS FOR WIND WAVE POWER DIAG  L850=',I3,
     *  ' L300=',I3,' L50=',I3)
      JLKDEX(1)=JEQ+JM*(L850-1)
      JLKDEX(2)=JEQ+JM*(L850-1+LM)
      JLKDEX(3)=JEQ+JM*(L300-1)
      JLKDEX(4)=JEQ+JM*(L300-1+LM)
      JLKDEX(5)=JEQ+JM*(L50-1)
      JLKDEX(6)=JEQ+JM*(L50-1+LM)
  100 DO KQ=1,5,2
         JLK=JLKDEX(KQ)
         CALL FFT (U(1,JLK,1),AN,BN)
         DO N=1,NMAX
            WAVE(1,IDACC9,N,KQ)=AN(N)
            WAVE(2,IDACC9,N,KQ)=BN(N)
         ENDDO
         CALL FFT (V(1,JLK,1),AN,BN)
         DO N=1,NMAX
            WAVE(1,IDACC9,N,KQ+1)=AN(N)
            WAVE(2,IDACC9,N,KQ+1)=BN(N)
         ENDDO
      ENDDO
      DO 150 I=1,IM
      PIJ50N=P(I,J50N)
      K=1
      L=1
      PL=SIG(1)*P(I,J50N)+PTOP
  130 L=L+1
      IF(L.GE.LS1) PIJ50N=PSFMPT
      PLM1=PL
      PL=SIG(L)*PIJ50N+PTOP
      IF (PMB(K).LT.PL.AND.L.LT.LM) GO TO 130
C**** ASSUME THAT PHI IS LINEAR IN LOG P
      SLOPE=(PHI(I,J50N,L-1)-PHI(I,J50N,L))/LOG(PLM1/PL)
  140 HTRD(I,K)=(PHI(I,J50N,L)+SLOPE*LOG(PMB(K)/PL))/GRAV-GHT(K)
      IF (K.GE.KM) GO TO 150
      K=K+1
      IF (PMB(K).LT.PL.AND.L.LT.LM) GO TO 130
      GO TO 140
  150 CONTINUE
      DO 160 KQ=7,KQMAX
      CALL FFT(HTRD(1,KQ-6),AN,BN)
      DO 160 N=1,NMAX
      WAVE(1,IDACC9,N,KQ)=AN(N)
  160 WAVE(2,IDACC9,N,KQ)=BN(N)
      CALL TIMER (MNOW,MINC,MDIAG)
      RETURN
      END SUBROUTINE DIAG7A

      SUBROUTINE DIAG7P
C****
C**** THIS ENTRY PRINTS THE TABLES
C****
      USE DAGCOM, only : wave,n12hrs_in_31day
      USE E001M12_COM, only : im,
     &     IDACC,JDATE,JDATE0,JMNTH0,JMONTH,JYEAR,JYEAR0,XLABEL
      IMPLICIT NONE

      INTEGER, DIMENSION(44) :: IPOWER
      DOUBLE PRECISION, DIMENSION(120) :: POWER
      DOUBLE PRECISION, DIMENSION(43) :: XPOWER
      DOUBLE PRECISION, DIMENSION(13) :: FPE
      DOUBLE PRECISION, DIMENSION(13,3,4) :: FFPE
      DOUBLE PRECISION, DIMENSION(43,9,3,4) :: FPOWER

      INTEGER, PARAMETER :: MMAX=12,NUAMAX=120,NUBMAX=15

      COMMON/D7COM/TITLE
      CHARACTER*64 TITLE(12)

      DOUBLE PRECISION, DIMENSION(12) :: SCALE
      DATA SCALE/1.,1., .1,1., .1,1., 4*1.D-3,1.D-4,1.D-5/

      DOUBLE PRECISION ::
     &     BYIA12,PNU,POWX,VAR

      INTEGER ::
     &     IA9X,IDACC9,JY,JY0,K,KPAGE,KQ,KTABLE,
     &     M,MMAXP1,N,NMAX,NS,NUA,NX

      NMAX=MIN0(9,IM/2)
      IDACC9=IDACC(9)
      IF (IDACC9.LE.MMAX) RETURN
C**** PATCH NEEDED IF SEVERAL RESTART FILES WERE ACCUMULATED
      IF (IDACC(12).LE.1) GO TO 320
      IDACC9=56
      IA9X=240*N12HRS_IN_31DAY
      BYIA12=1./IDACC(12)
      DO 310 N=1,IA9X
  310 WAVE(N,1,1,1)=WAVE(N,1,1,1)*BYIA12
  320 CONTINUE
      IF (IDACC9.GT.N12HRS_IN_31DAY) IDACC9=N12HRS_IN_31DAY
C****
C**** OUTPUT WAVE POWER AT THE EQUATOR
C****
      MMAXP1=MMAX+1
      DO 400 KPAGE=1,2
      JY0=JYEAR0-1900
      JY=JYEAR-1900
      WRITE (6,907) XLABEL(1:108),JDATE0,JMNTH0,JY0,JDATE,JMONTH,JY
      DO 390 KTABLE=1,3
      KQ=3*(KPAGE-1)+KTABLE
      WRITE (6,901) TITLE(KQ)
      DO 380 NX=1,NMAX
      N=NMAX+1-NX
      CALL MEM (WAVE(1,1,N,KQ),IDACC9,MMAX,NUAMAX,NUBMAX,POWER,FPE,
     *  VAR,PNU)
      POWX=.5*POWER(1)
      DO 330 NUA=2,27
  330 POWX=POWX+POWER(NUA)
      XPOWER(1)=SCALE(KQ)*POWX/26.5
      POWX=0.
      DO 340 NUA=28,34
  340 POWX=POWX+POWER(NUA)
      XPOWER(2)=SCALE(KQ)*POWX/7.
      XPOWER(3)=SCALE(KQ)*(POWER(35)+POWER(36)+POWER(37)+POWER(38))/4.
      XPOWER(4)=SCALE(KQ)*(POWER(39)+POWER(40))/2.
      DO 350 NUA=41,76
  350 XPOWER(NUA-36)=SCALE(KQ)*POWER(NUA)
      POWX=.5*POWER(1)
      DO 360 NUA=77,120
  360 POWX=POWX+POWER(NUA)
      XPOWER(41)=SCALE(KQ)*POWX/44.5
      XPOWER(42)=10.*SCALE(KQ)*VAR
      XPOWER(43)=1000.*SCALE(KQ)*(VAR-PNU)
      DO 370 NS=1,43
         FPOWER(NS,NX,KTABLE,KPAGE)=XPOWER(NS)
      IPOWER(NS)=XPOWER(NS)+.5
  370 CONTINUE
  380 WRITE (6,902) N,(IPOWER(NS),NS=1,43)
         DO 385 M=1,MMAXP1
  385    FFPE(M,KTABLE,KPAGE)=FPE(M)
  390 WRITE (6,903) (FPE(M),M=1,MMAXP1)
  400 CONTINUE
C****
C**** OUTPUT WAVE POWER AT 50 DEG NORTH
C****
      DO 500 KPAGE=3,4
      JY0=JYEAR0-1900
      JY=JYEAR-1900
      WRITE (6,907) XLABEL(1:108),JDATE0,JMNTH0,JY0,JDATE,JMONTH,JY
      DO 490 KTABLE=1,3
      KQ=3*(KPAGE-1)+KTABLE
  410 WRITE (6,911) TITLE(KQ)
      DO 480 NX=1,NMAX
      N=NMAX+1-NX
      CALL MEM (WAVE(1,1,N,KQ),IDACC9,MMAX,NUAMAX,NUBMAX,POWER,FPE,
     *  VAR,PNU)
      DO 420 M=1,MMAXP1
  420 FPE(M)=1000.*SCALE(KQ)*FPE(M)
      POWX=.5*POWER(1)
      DO 430 NUA=2,45
  430 POWX=POWX+POWER(NUA)
      XPOWER(1)=SCALE(KQ)*POWX/44.5
      DO 440 NUA=46,81
  440 XPOWER(NUA-44)=SCALE(KQ)*POWER(NUA)
      XPOWER(38)=SCALE(KQ)*(POWER(82)+POWER(83))/2.
      XPOWER(39)=SCALE(KQ)*(POWER(84)+POWER(85)+POWER(86)+POWER(87))/4.
      POWX=0.
      DO 450 NUA=88,94
  450 POWX=POWX+POWER(NUA)
      XPOWER(40)=SCALE(KQ)*POWX/7.
      POWX=.5*POWER(1)
      DO 460 NUA=95,120
  460 POWX=POWX+POWER(NUA)
      XPOWER(41)=SCALE(KQ)*POWX/26.5
      XPOWER(42)=10.*SCALE(KQ)*VAR
      XPOWER(43)=1000.*SCALE(KQ)*(VAR-PNU)
      DO 470 NS=1,43
         FPOWER(NS,NX,KTABLE,KPAGE)=XPOWER(NS)
      IPOWER(NS)=XPOWER(NS)+.5
  470 CONTINUE
  480 WRITE (6,902) N,(IPOWER(NS),NS=1,43)
         DO 485 M=1,MMAXP1
  485    FFPE(M,KTABLE,KPAGE)=FPE(M)
  490 WRITE (6,903) (FPE(M),M=1,MMAXP1)
  500 CONTINUE
      RETURN
C****
  901 FORMAT ('0',30X,A64,8X,'*1/60 (1/DAY)'/'   PERIOD EASTWARD--',
     *   35('---')/' N    -2      *-3   -3.3      -4       -5    -6   -7
     *.5  -10-12-15-20-30-60    60 30 20 15 12 10    7.5    6     5
     *   4*   VAR ERR'/'   --',40('---'))
  902 FORMAT (I2,41I3,I4,I4)
  903 FORMAT ('   --',40('---')/(1X,13F10.4))
  907 FORMAT ('1',A,I4,1X,A3,I3,' T0',I3,1X,A3,I3)
  911 FORMAT ('0',30X,A64,8X,'*1/60 (1/DAY)'/'   PERIOD EASTWARD--',
     *  35('---')/               ' N   *-4       -5    -6   -7.5  -10-12
     *-15-20-30-60    60 30 20 15 12 10    7.5    6     5        4
     * 3.3    3*       2    VAR ERR'/'   --',40('---'))
      END SUBROUTINE DIAG7P

      SUBROUTINE MEM (SERIES,ITM,MMAX,NUAMAX,NUBMAX,POWER,FPE,VAR,PNU)
      IMPLICIT NONE
      DIMENSION C(1800),S(1800),B1(62),B2(62),A(12),AA(11),P(13)
      DIMENSION SERIES(*),POWER(*),FPE(*)
      DOUBLE PRECISION PI,ARG,PP,POWERX,P,C,S,POWER,FPE
      COMPLEX*16 CI,CSUM,SS,A,AA,B1,B2,ANOM,ADEN
      COMPLEX*16 SERIES
      DOUBLE PRECISION :: PNU,VAR
      INTEGER ::
     &     I,ITM,L,M,MM1,MMAX,MMAXP1,NU,NUA,
     &     NUAMAX,NUB,NUBMAX,NUMAX,NUTM
      PI=3.141592653589793D0
      CI=DCMPLX(0.D0,1.D0)
      MMAXP1=MMAX+1
C**COSINE AND SINE FUNCTION
      NUMAX=NUAMAX*NUBMAX
      DO 20 NU=1,NUMAX
      ARG=2.0*PI*DFLOAT(NU)/DFLOAT(NUMAX)
      C(NU)=DCOS(ARG)
   20 S(NU)=DSIN(ARG)
   50 PP=0.0
      DO 60 I=1,ITM
   60 PP=PP+SERIES(I)*CONJG(SERIES(I))
      P(1)=PP/DFLOAT(ITM)
      VAR=P(1)
      M=1
      B1(1)=SERIES(1)
      B2(ITM-1)=SERIES(ITM)
      DO 70 I=2,ITM-1
      B1(I)=SERIES(I)
   70 B2(I-1)=SERIES(I)
      GO TO 80
  100 DO 110 I=1,M
  110 AA(I)=A(I)
      M=M+1
      DO 120 I=1,ITM-M
      B1(I)=B1(I)-CONJG(AA(M-1))*B2(I)
  120 B2(I)=B2(I+1)-AA(M-1)*B1(I+1)
   80 ANOM=DCMPLX(0.D0,0.D0)
      ADEN=DCMPLX(0.D0,0.D0)
      DO 90 I=1,ITM-M
      ANOM=ANOM+CONJG(B1(I))*B2(I)
   90 ADEN=ADEN+B1(I)*CONJG(B1(I))+B2(I)*CONJG(B2(I))
      A(M)=(ANOM+ANOM)/ADEN
      P(M+1)=P(M)*(1.0-CONJG(A(M))*A(M))
      IF (M.EQ.1) GO TO 100
  130 CONTINUE
      DO 140 I=1,M-1
  140 A(I)=AA(I)-A(M)*CONJG(AA(M-I))
      IF (M.LT.MMAX) GO TO 100
C**FINAL PREDICTION ERROR
      DO 150 M=1,MMAXP1
  150 FPE(M)=P(M)*DFLOAT(ITM+M-1)/DFLOAT(ITM-M+1)
      DO 180 NUA=1,NUAMAX
      POWERX=0.
C**FREQUENCY BAND AVERAGE
      DO 170 NUB=1,NUBMAX
      NU=NUB+NUA*NUBMAX+(NUMAX-3*NUBMAX-1)/2
      CSUM=1.
      DO 160 M=1,MMAX
      NUTM=MOD(NU*M-1,NUMAX)+1
  160 CSUM=CSUM-A(M)*(C(NUTM)-CI*S(NUTM))
  170 POWERX=POWERX+P(MMAXP1)/(CSUM*CONJG(CSUM))
      POWER(NUA)=.5*POWERX/DFLOAT(NUBMAX)
  180 CONTINUE
      PNU=0.0
      DO 210 L=1,NUAMAX
  210 PNU=PNU+POWER(L)
      PNU=PNU/(.5*NUAMAX)
      RETURN
      END SUBROUTINE MEM

      BLOCK DATA BDIJ
C****
C**** TITLES, LEGENDS AND CHARACTERS FOR DIAGIJ
C****
      COMMON/DIJCOM/TITLE1(18),TITLE2(18),TITLE3(18),
     *  LEGND1(10),LEGND2(14),ACHAR,BCHAR,CCHAR,DCHAR,ECHAR
C****
      CHARACTER*32 TITLE1
      DATA TITLE1/
     1   'TOPOGRAPHY (METERS)            ',
     *   'LAND COVERAGE                  ',
     *   'OCEAN ICE COVERAGE             ',
     *   'SNOW COVERAGE                  ',
     *   'SNOW DEPTH (MM H2O)            ',
     *   'SNOW AND ICE COVERAGE          ',
C
     7   'PRECIPITATION (MM/DAY)         ',
     *   'EVAPORATION (MM/DAY)           ',
     *   'SENSIBLE HEAT FLUX (WATTS/M**2)',
     *   'GROUND WETNESS                 ',
     *   'GROUND RUNOFF (MM/DAY)         ',
     *   'GROUND TEMPERATURE (DEGREES C) ',
C
     3   'SURFACE WIND SPEED (METERS/SEC)',
     *   'JET SPEED (METERS/SEC          ',
     *   'SURF WIND SPEED FROM U,V (M/S) ',
     *   'MTN WAVE MOM. FLUX (D/CM**2)   ',
     *   'JET DIRECTION (CW NOR)         ',
     *   'SURFACE WIND DIRECTION (CW NOR)'/
      CHARACTER*32 TITLE2
      DATA TITLE2/
     9   'TOTAL CLOUD COVER               ',
     *   'CONVECTIVE CLOUD COVER          ',
     *   'CLOUD TOP PRESSURE (MB)         ',
     *   'LOW LEVEL CLOUDINESS            ',
     *   'MIDDLE LEVEL CLOUDINESS         ',
     *   'HIGH LEVEL CLOUDINESS           ',
C
     5   'NET RAD. OF PLANET (WATTS/M**2) ',
     *   'NET RADIATION AT Z0 (WATTS/M**2)',
     *   'BRIGHTNESS TEMP THRU WNDW(DEG C)',
     *   'PLANETARY ALBEDO                ',
     *   'GROUND ALBEDO                   ',
     *   'VISUAL ALBEDO                   ',
C
     1   'NET THRML RADIATION (WATTS/M**2)',
     *   'NET HEAT AT Z0 (WATTS/M**2)     ',
     *   'TROP STATIC STABILITY (DEG K/KM)',
     *   'TOTAL NT DRY STAT ENR(10**14 WT)',
     *   'NT DRY STAT ENR BY ST ED(E14 WT)',
     *   'NT DRY STAT ENR BY TR ED(E14 WT)'/
      CHARACTER*32 TITLE3
      DATA TITLE3/
     7   '850 MB HEIGHT (METERS-1500)    ',
     *   '700 MB HEIGHT (METERS-3000)    ',
     *   '500 MB HEIGHT (METERS-5600)    ',
     *   '300 MB HEIGHT (METERS-9500)    ',
     *   '100 MB HEIGHT (METERS-16400)   ',
     *   ' 30 MB HEIGHT (METERS-24000)   ',
C
     3   'THICKNESS TEMPERATURE 1000-850 ',
     *   'THICKNESS TEMPERATURE 850-700  ',
     *   'THICKNESS TEMPERATURE 700-500  ',
     *   'THICKNESS TEMPERATURE 500-300  ',
     *   'THICKNESS TEMPERATURE 300-100  ',
     *   'THICKNESS TEMPERATURE 100-30   ',
C
     *   'TOTAL EARTH WATER (KG/M**2)    ',
     *   'LIQUID WATER PATH  (.1 KG/M**2)',
     *   'DEEP CONV CLOUD FREQUENCY      ',
     *   'SHALLOW CONV CLOUD FREQUENCY   ',
     *   'DEEP CONVECTIVE CLOUD COVER    ',
     *   'SHALLOW CONVECTIVE CLOUD COVER '/
C****
      CHARACTER*40 LEGND1
      DATA LEGND1/
     *  '0=0,1=5...9=45,A=50...K=100            ',
     *  '0=0...9=90,A=100...I=180...R=270       ',
     *  '1=.5...9=4.5,A=5...Z=17.5,+=MORE       ',
     *  '1=.1...9=.9,A=1...Z=3.5,+=MORE         ',
     *  '1=2...9=18,A=20...Z=70,+=MORE          ',
C
     *  '1=50...9=450,A=500...Z=1750,+=MORE     ',
     *  '1=100...9=900,A=1000...Z=3500,+=MORE   ',
     *  '1=20...9=180,A=200...Z=700,+=MORE      ',
     *  'A=1...Z=26,3=30...9=90,+=100-150,*=MORE',
     *  '0=0,A=.1...Z=2.6,3=3...9=9,+=10-15     '/
      CHARACTER*40 LEGND2
      DATA LEGND2/
     *  '-=LESS,Z=-78...0=0...9=27,+=MORE       ',
     *  '-=LESS,Z=-260...0=0...9=90,+=MORE      ',
     *  '-=LESS,Z=-520...0=0...9=180,+=MORE     ',
     *  '-=LESS,Z=-1300...0=0...9=450,+=MORE    ',
     *  '-=LESS,Z=-2600...0=0...9=900,+=MORE    ',
C
     *  '-=LESS,Z=-3900...0=0...9=1350,+=MORE   ',
     *  '-=LESS,Z=-5200...0=0...9=1800,+=MORE   ',
     *  '-=LESS,9=-.9...0=0,A=.1...Z=2.6,+=MORE ',
     *  '-=LESS,9=-45...0=0,A=5...K=45...+=MORE ',
     *  '-=LESS,9=-90...0=0,A=10...Z=260,+=MORE ',
C
     *  '-=LESS,9=-180...A=20...Z=520,+=MORE    ',
     *  '-=LESS,9=-9...0=0,A=1...Z=26,+=MORE    ',
     *  '-=LESS,9=-36...0=0,A=4...Z=104,+=MOR   ',
     *  '1=5...9=45,A=50...Z=175,+=MORE         '/
C****
      CHARACTER ACHAR*38,BCHAR*23,CCHAR*38,DCHAR*37,ECHAR*38
      DATA ACHAR/' 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ+'/
      DATA BCHAR/' 0123456789ABCDEFGHIJKX'/
      DATA CCHAR/'-9876543210ABCDEFGHIJKLMNOPQRSTUVWXYZ+'/
      DATA DCHAR/' 0ABCDEFGHIJKLMNOPQRSTUVWXYZ3456789+*'/
      DATA ECHAR/'-ZYXWVUTSRQPONMLKJIHGFEDCBA0123456789+'/
      END BLOCK DATA BDIJ

      SUBROUTINE DIAGIJ
C****
C**** THIS SUBROUTINE PRODUCES LATITUDE BY LONGITUDE MAPS OF
C****
C   K  IND                                                        IDACC
C****
C***1      TOPOGRAPHY (M)
C***2      LAND COVERAGE (10**-2)
C***3   1  OCEAN ICE COVERAGE (10**-2)                                4
C****   2  SNOW COVERAGE (10**-2)                                     4
C****   3  SNOW DEPTH (KG H2O/M**2)
C***6  29  SNOW AND ICE COVERAGE (PERCENT)
C****
C***7   5  PRECIPITATION (KG/M**2/86400 S)                            1
C****   6  EVAPORATION (KG/M**2/86400 S)                              1
C***9   4  SENSIBLE HEAT FLUX (WATTS/METER**2)
C**10   7  BETA, GROUND WETNESS (10**-2)                              3
C**11  32  GROUND RUNOFF FROM SURFACE (KG/M**2/86400 S)               1
C**12  28  FIRST LAYER GROUND TEMPERATURE (K-273.16)                  1
C****
C**13  34  SURFACE WIND SPEED (m/s)                                   3
C**14  39,40  JET SPEED (M/S)                                         4
C**15  36,37  SURFACE WIND SPEED FROM U,V (M/S)                       3
C**16  86  MOMENTUM FLUX BY MTN WAVES (.1 NT/M**2)
C**17  39,40  JET DIRECTION (CLOCKWISE FROM NORTH)                    0
C**18  36,37  SURFACE WIND DIRECTION (CLOCKWISE FROM NORTH)           0
C****
C**19  19  TOTAL CLOUD COVERAGE (PERCENT)
C**20  17  CLOUD COVERAGE FROM MOIST CONVECTION (PERCENT)
C**21  18/19   CLOUD TOP PRESSURE (MILLIBARS)
C**22  41  LOW LEVEL CLOUDINESS (PERCENT)
C**23  42  RMIDDLE LEVEL CLOUDINESS (PERCENT)
C**24  43  HIGH LEVEL CLOUDINESS (PERCENT)
C****
C**25  21+24  RADIATION BALANCE OF PLANET (WATTS/METER**2)
C**26  22  RADIATION BALANCE OF GROUND (WATTS/METER**2)
C**27  44  BRIGHTNESS TEMPERATURE THROUGH WINDOW REGION (K-273.16)
C**28  24/25  PLANETARY ALBEDO (PERCENT)
C**29  26/27  GROUND ALBEDO (PERCENT)
C**30  45/25  VISUAL ALBEDO (PERCENT)
C****
C**31  21  NET THERMAL RADIATION  (WATTA/METER**2)
C**32  23  NET HEAT AT GROUND (WATTS/METER**2)
C**33  31  TROPOSPHERIC STATIC STABILITY
C**34  20  TOTAL NORTH. TRANS. OF DRY STATIC ENERGY  (10**14 WATTS)
C**35      STAND. EDDY NORTH. TRANS. OF DRY STATIC ENERGY (10**14 WATTS)
C**36      TRANS. EDDY NORTH. TRANS. OF DRY STATIC ENERGY (10**14 WATTS)
C****
C**37  10  850 MB GEOPOTENTIAL HEIGHT (METERS-1500)
C****  11  700 MB GEOPOTENTIAL HEIGHT (METERS-3000)
C****  12  500 MB GEOPOTENTIAL HEIGHT (METERS-5600)
C****  13  300 MB GEOPOTENTIAL HEIGHT (METERS-9500)
C****  14  100 MB GEOPOTENTIAL HEIGHT (METERS-16400)
C****  15  30 MB GEOPOTENTIAL HEIGHT (METERS-24000)
C****
C**43   9,10  THICKNESS TEMPERATURE FROM 1000 TO 850 MB (DEGREES CENT.)
C****  10,11  THICKNESS TEMPERATURE FROM 850 TO 700 MB (DEGREES CENT.)
C****  11,12  THICKNESS TEMPERATURE FROM 700 TO 500 MB (DEGREES CENT.)
C****  12,13  THICKNESS TEMPERATURE FROM 500 TO 300 MB (DEGREES CENT.)
C****  13,14  THICKNESS TEMPERATURE FROM 300 TO 100 MB (DEGREES CENT.)
C****  14,15  THICKNESS TEMPERATURE FROM 100 TO 30 MB (DEGREES CENT.)
C****
C****  50  TOTAL EARTH WATER (KG H2O/M**2)
C****  81  LIQUID WATER PATH (.1 kg/M**2)
C****  84  DEEP CONVECTIVE CLOUD FREQUENCY (%)
C****  85  SHALLOW CONVECTIVE CLOUD FREQUENCY (%)
C****  83  DEEP CONVECTIVE CLOUD COVER (%)
C****  82  SHALLOW CONVECTIVE CLOUD COVER (%)
C****
      USE CONSTANT, only : grav,rgas,sday,twopi,sha,kapa
      USE E001M12_COM, only : im,jm,lm,fim,
     &     BYIM,DT,FLAND,IDACC,IDAY,IDAY0,JDATE,JDATE0,JEQ,
     &     JMNTH0,JMONTH,JYEAR,JYEAR0,NCNDS,NDYN,SKIPSE,TAU,
     &     TAU0,TOFDAY,TOFDY0,XLABEL,ZATMO
      USE GEOM, only : AREAG,DXV,DXYP,DXYV
      USE DAGCOM, only : ajk,kdiag,aij,kaij,aijk,tsfrez,
     &     IJ_PEV,IJ_TRNFP0,IJ_SRNFP0,IJ_SLP,IJ_TS !not a generic subr.
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(IM,JM,2) :: ENDE16
      common/ntdse/ende16 ! aka sendeg
      INTEGER, DIMENSION(3) :: MLAT,MGLOBE
      DOUBLE PRECISION, DIMENSION(3) :: FLAT,FNH,FGLOBE,GNUM,GDEN

      COMMON/DIJCOM/TITLE(3,18),LEGEND(10,24),ACHAR(38),BCHAR(23),
     *  CCHAR(38),DCHAR(37),ECHAR(38)
      CHARACTER*4 LEGEND,TITLE*32
      CHARACTER*1 ACHAR,BCHAR,CCHAR,DCHAR,ECHAR
C**** ACHAR/ ,0,1,...,8,9,A,B,...,Y,Z,+/
C**** BCHAR/ ,0,1,...,8,9,A,B,...,K,X/
C**** CCHAR/-,9,8,...,1,0,A,B,...,Y,Z,+/
C**** DCHAR/ ,0,A,B,...,Y,Z,3,4,...,8,9,+,*/
C**** ECHAR/-,Z,Y,...,B,A,0,1,...,8,9,+/
      CHARACTER*1 LINE(IM,3),LONGTD(36)
      DATA LONGTD/'+',35*' '/

      INTEGER, DIMENSION(54), PARAMETER :: IND=(/
     &   1, 1, 1, 2, 3,29,    5, 6, 4, 7,32,28,   34,39,36,86,39,36,
     &  19,17,18,41,42,43,   21,22,44,24,26,45,   21,23,31,20, 1, 2,
     &  10,11,12,13,14,15,    9,10,11,12,13,14,   50,81,84,85,83,82/)
      INTEGER, DIMENSION(54), PARAMETER :: IA=(/
     &   0, 0, 1, 4, 4, 4,    1, 1, 1, 1, 1, 1,    3, 4, 3, 1, 0, 0,
     &   2, 2, 0, 2, 2, 2,    2, 1, 2, 0, 0, 0,    2, 1, 4, 4, 4, 4,
     &   4, 4, 4, 4, 4, 4,    4, 4, 4, 4, 4, 4,    1, 1, 1, 1, 1, 1/)
      INTEGER, DIMENSION(3,18) :: ILEG
      DATA ILEG/7,3*1,9,1,   10,10,12, 1,18,11,    3, 5, 3, 4, 2, 2,
     *   1, 1, 6, 1, 1, 1,   13,20,11, 1, 1, 1,   13,13, 3,20,20,18,
     *  12,13,14,15,15,16,   11,11,11,11,11,11,    8, 3, 1, 1, 1, 1/
      DOUBLE PRECISION, DIMENSION(54) :: SCALE
      DATA SCALE/1.,3*100.,1.,100.,  3*1.,100.,2*1.,  6*1.,
     *  2*100.,1.,3*100.,  3*1.,3*100.,  2*1.,2.,21*1./
      DOUBLE PRECISION, DIMENSION(54) :: FAC
      DATA FAC/.01,3*.2,1.,.2,  2*10.,.1,.2,10.,.3333333,
     *  2.,.5,2.,10.,2*.1,  2*.2,.02,3*.2,  .05,.1,.3333333,3*.2,
     *  2*.05,2.,2*.1,10.,  .1,.05,.02,.01,.01,.006666667,  6*.3333333,
     *  .05,2.0,.2,.2,.2,.2/
      INTEGER, DIMENSION(KAIJ) :: JGRID
      DATA JGRID/19*1,2, 18*1,2*2, 40*1, 5*1,2*2,13*1/
      DOUBLE PRECISION, DIMENSION(7), PARAMETER :: PMB=(/
     &     1000.,850.,700.,500.,300.,100.,30. /)
      DOUBLE PRECISION, DIMENSION(7), PARAMETER :: GHT=(/
     &     0.,1500.,3000.,5600.,9500.,16400.,24000. /)
      DOUBLE PRECISION, PARAMETER :: P1000=1000.

      INTEGER ::
     &     I,IFRSTP,IHOUR,IHOUR0,ILINE,INC,IQ1,IQ2,IQ3,J,K,
     &     KC,KCOLMN,KPAGE,KR,KROW,KT,L,LASTP,M,N,NDEX,NDEX2,NDIAG,KM

      DOUBLE PRECISION ::
     &     A,BYIACC,BYIACN,BYIADA,DAREA,DTCNDS,DTSRCE,
     &     FDEN,FLATK,FNUM,PLAND,TAUDIF,ZNDE16,ZS,DPTI,PVTI,
     &     DE4TI,BYDPK,SZNDEG

C**** INITIALIZE CERTAIN QUANTITIES
c      SHA=RGAS/KAPA
      INC=1+(JM-1)/24
      ILINE=36*INC
      IQ1=1+IM/(4*INC)
      LONGTD(IQ1)=LONGTD(1)
      IQ2=1+IM/(2*INC)
      LONGTD(IQ2)=LONGTD(1)
      IQ3=1+3*IM/(4*INC)
      LONGTD(IQ3)=LONGTD(1)
C      JEQ=2.+.5*(JM-1.)
c      BYIM=1./FIM
      DTSRCE=NDYN*DT
      DTCNDS=NCNDS*DT
      SCALE(7)=SDAY/DTCNDS
      SCALE(8)=SDAY/DTSRCE
      SCALE(9)=1./DTSRCE
      SCALE(11)=SDAY/DTSRCE
      SCALE(13)=1.
      SCALE(16)=1E3/GRAV
      SCALE(26)=1./DTSRCE
      SCALE(32)=1./DTSRCE
      SCALE(33)=1.D3*GRAV*(P1000**KAPA)
      SCALE(34)=6.25E-14/GRAV
      SCALE(35)=SCALE(34)
      SCALE(36)=SCALE(34)
      DO 70 M=37,42
   70 SCALE(M)=1./GRAV
      DO 80 M=43,48
   80 SCALE(M)=1./(RGAS*LOG(PMB(M-42)/PMB(M-41)))
      SCALE(50)=10.
      SCALE(51)=100.
      SCALE(52)=100.
      SCALE(53)=100.
      SCALE(54)=100.
      IF (IDACC(12).LT.1) IDACC(12)=1
C****
      IHOUR0=TOFDY0+.5
      IHOUR = TOFDAY + .5
      TAUDIF=TAU-TAU0
CF*** NO PALMER INDEX FOR FINE GRID RUNS
      BYIADA=1./(IDACC(4)+1.D-20)
      IF (SKIPSE.EQ.1.) GO TO 160
C**** CACULATE STANDING AND TRANSIENT EDDY NORTHWARD TRANSPORT OF DSE
      KM=LM
      DO I=1,IM*JM*2
        ENDE16(I,1,1) = 0.
      END DO
      DO 115 J=2,JM
      DO 115 K=1,KM
      DPTI=0.
      PVTI=0.
      DE4TI=0.
      DO 105 I=1,IM
      IF (AIJK(I,J,K,4).EQ.0.) GO TO 105
      DPTI=DPTI+AIJK(I,J,K,4)
      BYDPK=1./(AIJK(I,J,K,4)+1.D-20)
      PVTI=PVTI+AIJK(I,J,K,2)
      DE4TI=DE4TI+AIJK(I,J,K,3)
      ENDE16(I,J,1)=ENDE16(I,J,1)
     *  +(AIJK(I,J,K,3)*AIJK(I,J,K,2)*BYDPK)
  105 CONTINUE
      SZNDEG=DE4TI*PVTI/(DPTI+1.D-20)
      DO 110 I=1,IM
  110 ENDE16(I,J,1)=ENDE16(I,J,1)-SZNDEG/FIM
  115 CONTINUE
      DO 130 J=2,JM
      ZNDE16=0.
      DO 120 L=1,LM
  120 ZNDE16=ZNDE16+(SHA*AJK(J,L,12)+AJK(J,L,14))
      ZNDE16=4.*ZNDE16*DXV(J)/FIM
      DO 130 I=1,IM
      ENDE16(I,J,1)=4.*ENDE16(I,J,1)*DXV(J)
  130 ENDE16(I,J,2)=AIJ(I,J,IJ_PEV)-ZNDE16-ENDE16(I,J,1)
C****
  160 NDIAG=54
      DO 180 N=1,KAIJ
      IF (JGRID(N).EQ.2) GO TO 180
      DO 170 I=1,IM
      AIJ(I,1,N)=AIJ(1,1,N)
  170 AIJ(I,JM,N)=AIJ(1,JM,N)
  180 CONTINUE
      IFRSTP=1
      LASTP=9
      IF (KDIAG(3).GT.0) LASTP=9-KDIAG(3)
      IF (KDIAG(3).LT.0) IFRSTP=-KDIAG(3)
      IF (KDIAG(3).LT.0) LASTP=IFRSTP
      DO 610 KPAGE=IFRSTP,LASTP
      WRITE (6,901) XLABEL
      WRITE (6,902) IDAY0,IHOUR0,JDATE0,JMNTH0,JYEAR0,TAU0,IDAY,IHOUR,
     *  JDATE,JMONTH,JYEAR,TAU,TAUDIF
      DO 610 KROW=1,2
      KR=2*(KPAGE-1)+KROW
      WRITE (6,903) (TITLE(K,KR),K=1,3)
      DO 200 KCOLMN=1,3
      FNH(KCOLMN)=0.
      FGLOBE(KCOLMN)=0.
      GNUM(KCOLMN)=0.
  200 GDEN(KCOLMN)=0.
      DO 550 J=JM,1,-1
      DO 210 I=1,IM*3
  210 LINE(I,1)=' '
      DO 510 KCOLMN=1,3
      FLATK=0.
      K=3*KR+KCOLMN-3
      NDEX=IND(K)
      BYIACC=1./(IDACC(IA(K))+1.D-20)
      GO TO (320,340,400,400,440,400, 440,440,460,400,420,460,
     *       380,300,300,475,240,240, 400,400,260,400,400,400,
     *       220,420,460,260,260,260, 460,460,380,420,280,280,
     *       460,460,460,460,460,460, 360,360,360,360,360,360,
     *       380,380,400,400,400,400),K
C**** SUM OF TWO ARRAYS
  220 DO 230 I=1,IM
      A=(AIJ(I,J,IJ_TRNFP0)+AIJ(I,J,IJ_SRNFP0))*SCALE(K)*BYIACC
CB       DIJMAP(I,J,K)=A
      FLATK=FLATK+A
      N=28.5+A*FAC(K)
      IF (N.LT.1 ) N=1
      IF (N.GT.38) N=38
  230 LINE(I,KCOLMN)=ECHAR(N)
      GO TO 500
C**** WIND DIRECTION
  240 IF (J.EQ.1) GO TO 500
      DO 250 I=1,IM
      A=360.*ATAN2(AIJ(I,J,NDEX)+1.D-20,AIJ(I,J,NDEX+1)+1.D-20)/TWOPI
CB       DIJMAP(I,J,K)=A
      FLATK=FLATK+A
      N=2.5+A*FAC(K)
      IF (N.LT.2) N=N+36
  250 LINE(I,KCOLMN)=ACHAR(N)
      GO TO 500
C**** RATIO OF 2 ARRAYS (MAINLY FOR ALBEDO)
  260 FNUM=0.
      FDEN=0.
      NDEX2=NDEX+1
      IF (NDEX.EQ.45) NDEX2=25
      DO 270 I=1,IM
      A=SCALE(K)*AIJ(I,J,NDEX)/(AIJ(I,J,NDEX2)+1.D-20)
      IF (NDEX.EQ.24 .OR. NDEX.EQ.26) A=100.-A
      FNUM=FNUM+AIJ(I,J,NDEX)
      FDEN=FDEN+AIJ(I,J,NDEX2)
      N=2.5+A*FAC(K)
      IF (A*FAC(K).GE.20.) N=23
      IF (AIJ(I,J,NDEX2).LE.0.) N=1
         IF (AIJ(I,J,NDEX2).LE.0.) A=0.
CB       DIJMAP(I,J,K)=A
  270 LINE(I,KCOLMN)=ACHAR(N)
      FLAT(KCOLMN)=SCALE(K)*FNUM/(FDEN+1.D-20)
      IF (NDEX.EQ.24 .OR. NDEX.EQ.26) FLAT(KCOLMN)=100.-FLAT(KCOLMN)
CB       DIJMAP(IM+1,J,K)=FLAT(KCOLMN)
      MLAT(KCOLMN)=NINT(FLAT(KCOLMN))
      GNUM(KCOLMN)=GNUM(KCOLMN)+FNUM*DXYP(J)
      GDEN(KCOLMN)=GDEN(KCOLMN)+FDEN*DXYP(J)
      IF (J.GT.INC) GO TO 510
      FGLOBE(KCOLMN)=SCALE(K)*GNUM(KCOLMN)/(GDEN(KCOLMN)+1.D-20)
      IF (NDEX.EQ.24.OR.NDEX.EQ.26) FGLOBE(KCOLMN)=100.-FGLOBE(KCOLMN)
      FGLOBE(KCOLMN)=FGLOBE(KCOLMN)*AREAG/ FIM
      GO TO 510
C**** STANDING AND TRANSIENT EDDY NORTHWARD TRANSPORTS OF DSE
  280 IF (SKIPSE.EQ.1.) GO TO 510
      DO 290 I=1,IM
      A=ENDE16(I,J,NDEX)*SCALE(K)*BYIACC
CB       DIJMAP(I,J,K)=A
      FLATK=FLATK+A
      N=11.5+A*FAC(K)
      IF (N.LT.1) N=1
      IF (N.GT.38) N=38
  290 LINE(I,KCOLMN)=CCHAR(N)
      FLAT(KCOLMN)=FLATK
      DAREA=DXYV(J)
      GO TO 505
C**** MAGNITUDE OF TWO PERPENDICULAR COMPONENTS
  300 IF (K.NE.15.AND.J.EQ.1) GO TO 500
      DO 310 I=1,IM
      A=SQRT(AIJ(I,J,NDEX)**2+AIJ(I,J,NDEX+1)**2)*SCALE(K)*BYIACC
CB       DIJMAP(I,J,K)=A
      FLATK=FLATK+A
      N=2.5+A*FAC(K)
      IF (N.GT.38) N=38
  310 LINE(I,KCOLMN)=ACHAR(N)
      GO TO 500
C**** SURFACE TOPOGRAPHY
  320 DO 330 I=1,IM
      ZS=ZATMO(I,J)/GRAV
      FLATK=FLATK+ZS
      N=2.5+.01*ZS
      IF (ZS.LE.0.) N=1
      IF (N.GT.38) N=38
  330 LINE(I,KCOLMN)=ACHAR(N)
      GO TO 500
C**** LAND COVERAGE
  340 DO 350 I=1,IM
      PLAND=FLAND(I,J)*100.
      FLATK=FLATK+PLAND
      N=2.5+PLAND*.2
      IF (PLAND.LE.0.) N=1
      IF (PLAND.GE.100.) N=23
  350 LINE(I,KCOLMN)=BCHAR(N)
      GO TO 500
C**** THICKNESS TEMPERATURES
  360 DO 370 I=1,IM
      A=((AIJ(I,J,NDEX+1)-AIJ(I,J,NDEX))*BYIACC
     *  +(GHT(NDEX-7)-GHT(NDEX-8))*GRAV)*SCALE(K)-273.16
CB       DIJMAP(I,J,K)=A
      FLATK=FLATK+A
      N=28.5+A*FAC(K)
      IF (N.LT.1) N=1
      IF (N.GT.38) N=38
  370 LINE(I,KCOLMN)=ECHAR(N)
      GO TO 500
C**** POSITIVE QUANTITIES UNIFORMLY SCALED
  380 DO 390 I=1,IM
      A=AIJ(I,J,NDEX)*SCALE(K)*BYIACC
CB       DIJMAP(I,J,K)=A
      FLATK=FLATK+A
      N=2.5+A*FAC(K)
      IF (A.EQ.0.) N=1
      IF (N.GT.38) N=38
  390 LINE(I,KCOLMN)=ACHAR(N)
      GO TO 500
C**** PERCENTAGES
  400 DO 410 I=1,IM
      A=AIJ(I,J,NDEX)*SCALE(K)*BYIACC
CB       DIJMAP(I,J,K)=A
      FLATK=FLATK+A
      N=2.5+A*FAC(K)
      IF (A.LE.0.) N=1
      IF (A*FAC(K).GE.20.) N=23
  410 LINE(I,KCOLMN)=BCHAR(N)
      GO TO 500
C**** SIGNED QUANTITIES UNIFORMLY SCALED (LETTERS +, NUMBERS -)
  420 DO 430 I=1,IM
      A=AIJ(I,J,NDEX)*SCALE(K)*BYIACC
CB       DIJMAP(I,J,K)=A
      FLATK=FLATK+A
      N=11.5+A*FAC(K)
      IF (N.LT.1) N=1
      IF (N.GT.38) N=38
  430 LINE(I,KCOLMN)=CCHAR(N)
      IF (K.EQ.34) FLATK=FLATK*FIM
      GO TO 500
C**** PRECIPITATION AND EVAPORATION
  440 DO 450 I=1,IM
      A=AIJ(I,J,NDEX)*SCALE(K)*BYIACC
CB       DIJMAP(I,J,K)=A
      FLATK=FLATK+A
      N=1
      IF (A.LE.0.) GO TO 450
      N=2.5+A*FAC(K)
      IF (N.GT.28) N=(N+263)/10
      IF (N.GT.35) N=(N+180)/6
      IF (N.GT.37) N=37
  450 LINE(I,KCOLMN)=DCHAR(N)
      GO TO 500
C**** SIGNED QUANTITIES UNIFORMLY SCALED (NUMBERS +, LETTERS -)
  460 DO 470 I=1,IM
      A=AIJ(I,J,NDEX)*SCALE(K)*BYIACC
CB       DIJMAP(I,J,K)=A
      FLATK=FLATK+A
      N=28.5+A*FAC(K)
      IF (N.LT.1 ) N=1
      IF (N.GT.38) N=38
  470 LINE(I,KCOLMN)=ECHAR(N)
      GO TO 500
C**** ABSOLUTE VALUE OF QUANTITIES UNIFORMLY SCALED
  475 DO 477 I=1,IM
      A=ABS(AIJ(I,J,NDEX)*SCALE(K)*BYIACC)
CB       DIJMAP(I,J,K)=A
      FLATK=FLATK+A
      N=2.5+A*FAC(K)
      IF(A.EQ.0.) N=1
      IF (N.GT.38) N=38
  477 LINE(I,KCOLMN)=ACHAR(N)
      GO TO 500
C**** POSITIVE QUANTITIES NON-UNIFORMLY SCALED
  480 DO 490 I=1,IM
      A=AIJ(I,J,NDEX)*SCALE(K)*BYIACC
CB       DIJMAP(I,J,K)=A
      FLATK=FLATK+A
      N=2.5+A*FAC(K)
      IF (N.GE.13) N=(N+123)/10
      IF (N.GT.38) N=38
  490 LINE(I,KCOLMN)=ACHAR(N)
      GO TO 500
C**** LENGTH OF GROWING SEASON
  491 DO 492 I=1,IM
      A=TSFREZ(I,J,2)-TSFREZ(I,J,1)
      IF (A.LT.0.) A=A+365.
      FLATK=FLATK+A
      N=2.5+A*FAC(K)
      IF (A.LE.0.) N=1
      IF (N.GT.38) N=38
  492 LINE(I,KCOLMN)=ACHAR(N)
      GO TO 500
C**** PALMER DROUGHT INDEX
  493 DO 494 I=1,IM
      A=0.
      FLATK=FLATK+A
      N=11.5+A*FAC(K)
      IF (N.LT.1 ) N=1
      IF (N.GT.38) N=38
  494 LINE(I,KCOLMN)=CCHAR(N)
  500 FLAT(KCOLMN)=FLATK*BYIM
      DAREA=DXYP(J)
      IF (JGRID(NDEX).EQ.2) DAREA=DXYV(J)
      IF (J.GE.JEQ) FNH(KCOLMN)=FNH(KCOLMN)+FLAT(KCOLMN)*DAREA
  505 FGLOBE(KCOLMN)=FGLOBE(KCOLMN)+FLAT(KCOLMN)*DAREA
CB       DIJMAP(IM+1,J,K)=FLAT(KCOLMN)
      MLAT(KCOLMN)=NINT(FLAT(KCOLMN))
  510 CONTINUE
      IF (MOD(J,INC).NE.0) GO TO 550
      GO TO (524,520, 520,520, 520,520, 521,520, 526,520, 526,524,
     *       527,527, 520,520, 524,520, 527,527),KR
  520 WRITE (6,910) (FLAT(KC),(LINE(I,KC),I=1,ILINE,INC),KC=1,3)
      GO TO 530
  521 WRITE (6,911) (FLAT(KC),(LINE(I,KC),I=1,ILINE,INC),KC=1,2),
     *  MLAT(3),(LINE(I,3),I=1,ILINE,INC)
      GO TO 530
  524 WRITE (6,914) MLAT(1),(LINE(I,1),I=1,ILINE,INC),
     *  (FLAT(KC),(LINE(I,KC),I=1,ILINE,INC),KC=2,3)
      GO TO 530
  526 WRITE (6,916) (MLAT(KC),(LINE(I,KC),I=1,ILINE,INC),KC=1,2),
     *  FLAT(3),(LINE(I,3),I=1,ILINE,INC)
      GO TO 530
  527 WRITE (6,917) (MLAT(KC),(LINE(I,KC),I=1,ILINE,INC),KC=1,3)
  530 DO 540 I=1,IM
      IF (FLAND(I,J).GE..5) GO TO 540
      LINE(I,1)=' '
      LINE(I,2)=' '
      LINE(I,3)=' '
  540 CONTINUE
C     WRITE (6,906) ((LINE(I,KC),I=1,ILINE,INC),KC=1,3)
C     WRITE (6,906) ((LINE(I,KC),I=1,ILINE,INC),KC=1,3)
      WRITE (6,906) ((LINE(I,KC),I=1,ILINE,INC),KC=1,3)
  550 CONTINUE
      DO 555 KCOLMN=1,3
      FNH(KCOLMN)=2.*FNH(KCOLMN)*FIM/AREAG
      FGLOBE(KCOLMN)=FGLOBE(KCOLMN)*FIM/AREAG
CB       DIJMPG(KCOLMN,KROW,KPAGE)=FGLOBE(KCOLMN)
  555 MGLOBE(KCOLMN)=NINT(FGLOBE(KCOLMN))
      IF (KR.EQ.2) CALL KEYIJ (FGLOBE(3),FNH(3))
      GO TO (574,570, 570,570, 570,570, 571,570, 577,570, 576,570,
     *       577,577, 570,570, 574,570, 577,577),KR
  570 WRITE (6,910) (FGLOBE(KC),LONGTD,KC=1,3)
      GO TO 600
  571 WRITE (6,911) FGLOBE(1),LONGTD,FGLOBE(2),LONGTD,MGLOBE(3),LONGTD
      GO TO 600
  574 WRITE (6,914) MGLOBE(1),LONGTD,FGLOBE(2),LONGTD,FGLOBE(3),LONGTD
      GO TO 600
  576 WRITE (6,916) MGLOBE(1),LONGTD,MGLOBE(2),LONGTD,FGLOBE(3),LONGTD
      GO TO 600
  577 WRITE (6,917) (MGLOBE(KC),LONGTD,KC=1,3)
  600 WRITE (6,909) ((LEGEND(KT,ILEG(KCOLMN,KR)),KT=1,10),KCOLMN=1,2),
     *  (LEGEND(KT,ILEG(3,KR)),KT=1,9)
  610 CONTINUE
  690 CONTINUE
C****
C**** PRODUCE FULL PAGE I,J MAPS
C****
      CALL IJMAP (1,AIJ(1,1,IJ_SLP),BYIADA)
      BYIACN=1./(IDACC(3)+1.D-20)
      CALL IJMAP (2,AIJ(1,1,IJ_TS),BYIACN)
C     CALL IJMAP (4,AIJ(1,1,IJ_SLP1),BYIADA)
C     CALL IJMAP (5,AIJ(1,1,IJ_TS1),BYIADA)
      RETURN
C****
  901 FORMAT ('1',A)
  902 FORMAT ('0',15X,'DAY',I6,', HR',I2,' (',I2,A5,I4,')',F9.0,
     *      '   TO   DAY',I6,', HR',I2,' (',I2,A5,I4,')',F9.0,
     *  '    DIF',F5.0,' HR')
  903 FORMAT ('0',6X,A32,13X,A32,13X,A32)
  906 FORMAT ('+',6X,36A1,9X,36A1,9X,36A1)
  909 FORMAT (7X,10A4,5X,10A4,5X,9A4)
  910 FORMAT (1X,F5.1,1X,36A1,F8.1,1X,36A1,F8.1,1X,36A1)
  911 FORMAT (1X,F5.1,1X,36A1,F8.1,1X,36A1,I8,1X,36A1)
  914 FORMAT (1X,I5,1X,36A1,F8.1,1X,36A1,F8.1,1X,36A1)
  916 FORMAT (1X,I5,1X,36A1,I8,1X,36A1,F8.1,1X,36A1)
  917 FORMAT (1X,I5,1X,36A1,I8,1X,36A1,I8,1X,36A1)
      END SUBROUTINE DIAGIJ

      SUBROUTINE IJMAP (NT,ARRAY,BYIACC)
      USE E001M12_COM, only : im,jm,
     &     BYIM,FLAND,IDAY,IDAY0,JDATE,JDATE0,JMNTH0,JMONTH,
     &     JYEAR,JYEAR0,TAU,TAU0,TOFDAY,TOFDY0,XLABEL
      USE GEOM, only : JLAT,LON
      IMPLICIT NONE

      INTEGER :: NT
      DOUBLE PRECISION :: BYIACC
      DOUBLE PRECISION, DIMENSION(IM,JM) :: ARRAY


      CHARACTER*1 IDX(12),BLANK,TITLE(5)*48
      CHARACTER(LEN=3), DIMENSION(IM) :: LINE
      CHARACTER(LEN=9) :: AVG
      DATA IDX/'0','1','2','3','4','5','6','7','8','9','-','*'/
      DATA BLANK/' '/
      DATA TITLE/
C****
C**** THIS SUBROUTINE PRODUCES NUMERICAL LATITUDE BY LONGITUDE MAPS OF
C****
     *  'SEA LEVEL PRESSURE (MB-1000) ',
     *  'SURFACE TEMPERATURE (DEGREES C) ',
     *  'INSTANTANEOUS 850 MB HEIGHTS (DEKAMETERS-100)',
     *  'SEA LEVEL PRESSURE (MB-1000)  (USING T1)',
     *  'SURFACE TEMPERATURE (DEG C)  (LAPSE RATE FROM T1'/

      DOUBLE PRECISION :: A,FLAT,TAUDIF
      INTEGER :: I,IA,IHOUR,IHOUR0,INC,J,JA,JX,K,LD

C****
C**** INITIALIZE CERTAIN QUANTITIES
C****
c      BYIM=1./IM
      INC=(IM+35)/36
      LD=360/IM
      DO 40 I=1,IM
      WRITE(LINE(I),'(I3)') I
   40 CONTINUE
      AVG='     MEAN'
      IHOUR0=TOFDY0+.5
      IHOUR = TOFDAY + .5
      TAUDIF=TAU-TAU0
      WRITE(6,901)XLABEL
      WRITE(6,902)IDAY0,IHOUR0,JDATE0,JMNTH0,JYEAR0,TAU0,IDAY,IHOUR,
     *  JDATE,JMONTH,JYEAR,TAU,TAUDIF
      WRITE(6,900) TITLE(NT)
      WRITE (6,910) (LINE(I),I=1,IM,INC),AVG
      WRITE(6,940)
      WRITE(6,940)
C**** OUTSIDE J LOOP
      DO 300 J=JM,1,-1
      FLAT=0.
      DO 250 I=1,IM
      A=ARRAY(I,J)*BYIACC
CB       QIJMAP(I,J,NT)=A
      FLAT=FLAT+A
      IF (A.LT.999.5.AND.A.GE.-99.5) GO TO 140
      DO 100 K=1,3
  100 LINE(I)(K:K)=IDX(12)
      GO TO 250
  140 DO 150 K=1,3
  150 LINE(I)(K:K)=BLANK
      JA=NINT(A)
      IA=IABS(JA)
      IF (IA.GT.99) GO TO 210
      IF (IA-9) 230,230,220
  210 LINE(I)(1:1)=IDX(IA/100+1)
      IA=MOD(IA,100)
  220 LINE(I)(2:2)=IDX(IA/10+1)
      IA=MOD(IA,10)
  230 LINE(I)(3:3)=IDX(IA+1)
      IF (JA.GE.0) GO TO 250
      IF (JA+9) 240,245,245
  240 LINE(I)(1:1)=IDX(11)
      GO TO 250
  245 LINE(I)(2:2)=IDX(11)
  250 CONTINUE
      FLAT=FLAT*BYIM
CB       QIJMAP(IM+1,J,NT)=FLAT
      WRITE(AVG,'(F9.2)') FLAT
      WRITE (6,920) JLAT(J,1),J,(LINE(I),I=1,IM,INC),AVG
      DO 260 I=1,IM
      IF (FLAND(I,J).GE..5) GO TO 260
      DO 255 K=1,3
  255 LINE(I)(K:K)=BLANK
  260 CONTINUE
C     WRITE (6,925) (LINE(I),I=1,IM,INC)
      WRITE (6,925) (LINE(I),I=1,IM,INC)
      WRITE (6,925) (LINE(I),I=1,IM,INC)
  300 IF (JM.LE.24) WRITE (6,940)
      WRITE (6,930) (LON(I,1),I=1,IM,INC*2)
      RETURN
C****
  900 FORMAT('0',45X,A48)
  901 FORMAT ('1',A)
  902 FORMAT ('0',15X,'DAY',I6,', HR',I2,' (',I2,A5,I4,')',F9.0,
     *      '   TO   DAY',I6,', HR',I2,' (',I2,A5,I4,')',F9.0,
     *  '    DIF',F5.0,' HR')
  910 FORMAT('0LAT  J/I  ',36A3,A9)
  920 FORMAT(2I4,3X,36A3,A9)
  925 FORMAT('+',10X,36A3)
  930 FORMAT('0  LONG',2X,19F6.1)
  940 FORMAT(' ')
      END SUBROUTINE IJMAP

      BLOCK DATA BDCNS
C****
C**** TITLES FOR SUBROUTINE DIAG9
C****
      COMMON/D9COM/TITLE1,TITLE2,TITLE3,TITLE4,TITLE5
      CHARACTER*32 TITLE1(11)
      DATA TITLE1/
     *  ' INSTANTANE AM (10**9 J*S/M**2) ',
     *  ' CHANGE OF AM BY ADVECTION      ',
     *  ' CHANGE OF AM BY CORIOLIS FORCE ',
     *  ' CHANGE OF AM BY ADVEC + COR    ',
     *  ' CHANGE OF AM BY PRESSURE GRAD  ',
     *  ' CHANGE OF AM BY DYNAMICS       ',
     *  ' CHANGE OF AM BY SURFACE FRIC   ',
     *  ' CHANGE OF AM BY STRATOS DRAG   ',
     *  ' CHANGE OF AM BY FILTER         ',
     *  ' CHANGE OF AM BY DAILY RESTOR   ',
     *  ' SUM OF CHANGES (10**2 J/M**2)  '/
      CHARACTER*32 TITLE2(12)
      DATA TITLE2/
     *  '0INSTANTANEOUS KE (10**3 J/M**2)',
     *  ' CHANGE OF KE BY ADVECTION      ',
     *  ' CHANGE OF KE BY CORIOLIS FORCE ',
     *  ' CHANGE OF KE BY ADVEC + COR    ',
     *  ' CHANGE OF KE BY PRESSURE GRAD  ',
     *  ' CHANGE OF KE BY DYNAMICS       ',
     *  ' CHANGE OF KE BY MOIST CONVEC   ',
     *  ' CHANGE OF KE BY SURF + DRY CONV',
     *  ' CHANGE OF KE BY STRATOS DRAG   ',
     *  ' CHANGE OF KE BY FILTER         ',
     *  ' CHANGE OF KE BY DAILY RESTOR   ',
     *  ' SUM OF CHANGES (10**-3 W/M**2) '/
      CHARACTER*32 TITLE3(5)
      DATA TITLE3/
     *  ' INSTANTANEOUS MASS (KG/M**2)   ',
     *  ' CHANGE OF MASS BY DYNAMICS     ',
     *  ' CHANGE OF MASS BY FILTER       ',
     *  ' CHANGE OF MASS BY DAILY RESTOR ',
     *  ' SUM CHANGES (10**-8 KG/S/M**2) '/
      CHARACTER*32 TITLE4(8)
      DATA TITLE4/
     *  '0INSTANTANE TPE (10**5 J/M**2)  ',
     *  ' CHANGE OF TPE BY DYNAMICS      ',
     *  ' CHANGE OF TPE BY CONDENSATION  ',
     *  ' CHANGE OF TPE BY RADIATION     ',
     *  ' CHANGE OF TPE BY SURFACE INTER ',
     *  ' CHANGE OF TPE BY FILTER        ',
     *  ' CHANGE OF TPE BY DAILY RESTOR  ',
     *  ' SUM OF CHANGES (10**-2 W/M**2) '/
      CHARACTER*32 TITLE5(6)
      DATA TITLE5/
     *  '0INSTANT WATER (10**-2 KG/M**2) ',
     *  ' CHANGE OF WATER BY DYNAMICS    ',
     *  ' CHANGE OF WATER BY CLOUDS      ',
     *  ' CHANGE OF WATER BY SURFACE EVAP',
     *  ' CHANGE OF WATER BY DAILY RESTOR',
     *  ' SUM CHANGES (10**-8 KG/S/M**2) '/
      END BLOCK DATA BDCNS

      SUBROUTINE DIAG9A (M)
!@sum  DIAG9A Keeps track of the conservation properties of angular
!@sum  momentum, kinetic energy, mass, total potential energy and water
!@auth Gary Russell
!@ver  1.0

      USE CONSTANT, only : kapa,omega,sha
      USE E001M12_COM, only : im,jm,lm,fim,
     &     DSIG,LS1,MDIAG,P,PTOP,T,U,V,ZATMO,Q,WM,PSFMPT,PSTRAT
      USE GEOM, only :
     &     COSV,DXYN,DXYS,DXYV,IMAXJ,RADIUS
      USE DYNAMICS, only : pk
      USE DAGCOM, only : consrv
      IMPLICIT NONE

!@var M index denoting from where DIAG9A is called
      INTEGER, INTENT(IN) :: M
!@var NAMOFM,NKEOFM,NMSOFM,NPEOFM,NWVOFM indices for CONSRV array
      INTEGER :: NAMOFM(8) = (/1,6,1,1,7,8,9,10/)
      INTEGER :: NKEOFM(8) = (/1,17,18,1,19,20,21,22/)
      INTEGER :: NMSOFM(8) = (/1,25,1,1,1,1,26,27/)
      INTEGER :: NPEOFM(8) = (/1,30,31,32,33,1,34,35/)
      INTEGER :: NWVOFM(8) = (/1,38,39,1,40,1,1,41/)

      DOUBLE PRECISION, DIMENSION(JM) :: PI,AM,RKE,RMASS,TPE,VAPOR
      INTEGER :: I,IMAX,IP1,J,L,MINC,MNOW,N
      DOUBLE PRECISION :: RKEI,RKEIL,SGEOI,TPEI,TPEIL,UMI,UMIL
C****
C**** THE PARAMETER M INDICATES WHEN DIAG9A IS BEING CALLED
C**** M=1  INITIALIZE CURRENT A.M., K.E., MASS, T.P.E. AND W.V.
C****   2  AFTER DYNAMICS
C****   3  AFTER CONDENSATION
C****   4  AFTER RADIATION
C****   5  AFTER SURFACE INTERACTION AND DRY CONVECTION
C****   6  AFTER STRATOSPHERIC DRAG
C****   7  AFTER FILTER
C****   8  AFTER DAILY RESTORATION
C****
      GO TO (100,100,200,400,100,100,100,100),M
C****
C**** ANGULAR MOMENTUM
C****
  100 PI(1)=FIM*P(1,1)
      PI(JM)=FIM*P(1,JM)
      DO J=2,JM-1
         PI(J)=0.
         DO I=1,IM
            PI(J)=PI(J)+P(I,J)
         END DO
      END DO
      DO J=2,JM
         UMIL=0.
         DO L=1,LM
            UMI=0.
            I=IM
            DO IP1=1,IM
               IF(L.LT.LS1) THEN
                  UMI=UMI+U(I,J,L)*((P(I,J-1)+P(IP1,J-1))*DXYN(J-1)
     *                 +(P(I,J)+P(IP1,J))*DXYS(J))
               ELSE
                  UMI=UMI+U(I,J,L)*(2.*PSFMPT*DXYV(J))
               END IF
               I=IP1
            END DO
            UMIL=UMIL+UMI*DSIG(L)
         END DO
         AM(J)=RADIUS*OMEGA*COSV(J)*((PI(J-1)*DXYN(J-1)+PI(J)*DXYS(J))
     *        +FIM*PSTRAT*DXYV(J))+.5*UMIL
      END DO
      IF (M.GT.1) THEN
         N=NAMOFM(M)
         DO J=2,JM
            CONSRV(J,N)=CONSRV(J,N)+(AM(J)-CONSRV(J,1))
         END DO
      END IF
      DO J=2,JM
         CONSRV(J,1)=AM(J)
      END DO
C****
C**** KINETIC ENERGY
C****
 200  DO J=2,JM
         RKEIL=0.
         DO L=1,LM
            RKEI=0.
            I=IM
            DO IP1=1,IM
               IF(L.LT.LS1) THEN
                  RKEI=RKEI+(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))
     *                 *((P(I,J-1)+P(IP1,J-1))*DXYN(J-1)+(P(I,J)+P(IP1,J
     *                 ))*DXYS(J))
               ELSE
                  RKEI=RKEI+(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))*
     *                 (2.*PSFMPT*DXYV(J))
               END IF
               I=IP1
            END DO
            RKEIL=RKEIL+RKEI*DSIG(L)
            RKE(J)=RKEIL
         END DO
      END DO
      IF (M.GT.1) THEN
         N=NKEOFM(M)
         DO J=2,JM
            CONSRV(J,N)=CONSRV(J,N)+(RKE(J)-CONSRV(J,12))
         END DO
      END IF
      DO J=2,JM
         CONSRV(J,12)=RKE(J)
      END DO
      IF (M.EQ.6) GO TO 495
      IF (M.EQ.3.OR.M.EQ.5) GO TO 400
C****
C**** MASS
C****
 300  RMASS(1)=FIM*(P(1,1)+PSTRAT)
      RMASS(JM)=FIM*(P(1,JM)+PSTRAT)
      DO J=2,JM-1
         RMASS(J)=FIM*PSTRAT
         DO I=1,IM
            RMASS(J)=RMASS(J)+P(I,J)
         END DO
      END DO
      IF (M.GT.1) THEN
         N=NMSOFM(M)
         DO J=1,JM
            CONSRV(J,N)=CONSRV(J,N)+(RMASS(J)-CONSRV(J,24))
         END DO
      END IF
      DO J=1,JM
         CONSRV(J,24)=RMASS(J)
      END DO
C****
C**** TOTAL POTENTIAL ENERGY
C****
 400  DO J=1,JM
         IMAX=IMAXJ(J)
         TPEIL=0.
         DO L=1,LM
            TPEI=0.
            DO I=1,IMAX
               IF(L.LT.LS1) THEN
                  TPEI=TPEI+T(I,J,L)*PK(L,I,J)*P(I,J)
               ELSE
                  TPEI=TPEI+T(I,J,L)*PK(L,I,J)*PSFMPT
               END IF
            END DO
            TPEIL=TPEIL+TPEI*DSIG(L)
         END DO
         SGEOI=0.
         DO I=1,IMAX
            SGEOI=SGEOI+ZATMO(I,J)*(P(I,J)+PTOP)
         END DO
         TPE(J)=SGEOI+TPEIL*SHA
      END DO
      TPE(1)=FIM*TPE(1)
      TPE(JM)=FIM*TPE(JM)
      IF (M.GT.1) THEN
         N=NPEOFM(M)
         DO J=1,JM
            CONSRV(J,N)=CONSRV(J,N)+(TPE(J)-CONSRV(J,29))
         END DO
      END IF
      DO J=1,JM
         CONSRV(J,29)=TPE(J)
      END DO
      IF (M.EQ.4.OR.M.EQ.7) GO TO 495
C****
C**** TOTAL WATER MASS
C****
      DO J=1,JM
         IMAX=IMAXJ(J)
         VAPOR(J) = 0.
         DO L=1,LS1-1
            DO I=1,IMAX
               VAPOR(J) = VAPOR(J) + (Q(I,J,L)+WM(I,J,L))*P(I,J)*DSIG(L)
            END DO
         END DO
         DO L=LS1,LM
            DO I=1,IMAX
               VAPOR(J) = VAPOR(J) + (Q(I,J,L)+WM(I,J,L))*PSFMPT*DSIG(L)
            END DO
         END DO
      END DO
      VAPOR(1) = VAPOR( 1)*IM
      VAPOR(JM)= VAPOR(JM)*IM
      IF (M.GT.1) THEN
         N=NWVOFM(M)
         DO J=1,JM
            CONSRV(J,N) = CONSRV(J,N) + (VAPOR(J)-CONSRV(J,37))
         END DO
      END IF
      DO J=1,JM
         CONSRV(J,37) = VAPOR(J)
      END DO
C****
 495  CALL TIMER (MNOW,MINC,MDIAG)
      RETURN
      END SUBROUTINE DIAG9A

      SUBROUTINE DIAG9D (M,DT1,UX,VX,DUT,DVT,PIT)
!@sum  DIAG9D Keeps track of the conservation properties of angular
!@sum  momentum and kinetic energy inside dynamics routines
!@auth Gary Russell
!@ver  1.0

      USE CONSTANT, only : omega
      USE E001M12_COM, only : im,jm,lm,fim,mdiag,mdyn
      USE GEOM, only : COSV,RADIUS,RAVPN,RAVPS
      USE DAGCOM, only : consrv
      IMPLICIT NONE
C****
C**** THE PARAMETER M INDICATES WHEN DIAG9D IS BEING CALLED
C**** M=1  AFTER ADVECTION IN DYNAMICS
C****   2  AFTER CORIOLIS FORCE IN DYNAMICS
C****   3  AFTER PRESSURE GRADIENT FORCE IN DYNAMICS
C****
!@var M index denoting from where DIAG9D is called
      INTEGER, INTENT(IN) :: M
!@var DT1 current time step
      DOUBLE PRECISION, INTENT(IN) :: DT1
!@var UX,VX current velocities
      DOUBLE PRECISION, INTENT(IN), DIMENSION(IM,JM,LM) :: UX,VX
!@var DUT,DVT current momentum changes
      DOUBLE PRECISION, INTENT(IN), DIMENSION(IM,JM,LM) :: DUT,DVT
!@var PIT current pressure tendency
      DOUBLE PRECISION, INTENT(IN), DIMENSION(IM,JM) :: PIT

      DOUBLE PRECISION, DIMENSION(JM) :: PI

      INTEGER :: MCLOCK ! external
      INTEGER :: I,J,L,MBEGIN,MEND

      DOUBLE PRECISION :: DUTI,DUTIL,RKEI,RKEIL

      MBEGIN = MCLOCK ()
C****
C**** CHANGE OF ANGULAR MOMENTUM AND KINETIC ENERGY BY ADVECTION
C****
      IF (M.eq.1) THEN
         PI(1)=FIM*PIT(1,1)
         PI(JM)=FIM*PIT(1,JM)
         DO J=2,JM-1
            PI(J)=0.
            DO I=1,IM
               PI(J)=PI(J)+PIT(I,J)
            END DO
         END DO
         DO J=2,JM
            DUTIL=0.
            RKEIL=0.
            DO L=1,LM
               DUTI=0.
               RKEI=0.
               DO I=1,IM
                  DUTI=DUTI+DUT(I,J,L)
                  RKEI=RKEI+(UX(I,J,L)*DUT(I,J,L)+VX(I,J,L)*DVT(I,J,L))
               END DO
               DUTIL=DUTIL+DUTI
               RKEIL=RKEIL+RKEI
            END DO
            CONSRV(J,2)=CONSRV(J,2)+(DUTIL+DT1*RADIUS*OMEGA*COSV(J)*
     *           (PI(J-1)*RAVPN(J-1)+PI(J)*RAVPS(J)))
            CONSRV(J,13)=CONSRV(J,13)+RKEIL
         END DO
      ELSE
C****
C**** CHANGE OF ANGULAR MOMENTUM AND KINETIC ENERGY BY CORIOLIS AND
C**** PRESSURE GRADIENT FORCES
C****
         DO J=2,JM
            DUTIL=0.
            RKEIL=0.
            DO L=1,LM
               DUTI=0.
               RKEI=0.
               DO I=1,IM
                  DUTI=DUTI+DUT(I,J,L)
                  RKEI=RKEI+(UX(I,J,L)*DUT(I,J,L)+VX(I,J,L)*DVT(I,J,L))
               END DO
               DUTIL=DUTIL+DUTI
               RKEIL=RKEIL+RKEI
            END DO
            CONSRV(J,2*M-1)=CONSRV(J,2*M-1)+DUTIL
            CONSRV(J,2*M+10)=CONSRV(J,2*M+10)+RKEIL
         END DO
      END IF
C****
      MEND = MCLOCK ()
      MDIAG = MDIAG + (MEND-MBEGIN)
      MDYN  = MDYN  - (MEND-MBEGIN)
      RETURN
      END SUBROUTINE DIAG9D
C****
C****
      SUBROUTINE DIAG9P
C****
C**** THIS ENTRY PRODUCES TABLES OF CONSERVATION QUANTITIES
C****
      USE CONSTANT, only : grav,sday,twopi

      USE E001M12_COM, only : im,jm,lm,fim,
     &     DT,IDACC,IDAY,IDAY0,JDATE,JDATE0,JEQ,JMNTH0,
     &     JMONTH,JYEAR,JYEAR0,NDYN,NFILTR,
     &     TAU,TAU0,TOFDAY,TOFDY0,XLABEL

      USE GEOM, only :
     &     AREAG,COSV,DLON,DXYP,DXYV,RADIUS,JLAT
      USE DAGCOM, only : consrv,kcon
      IMPLICIT NONE

      INTEGER, DIMENSION(JM) :: MAREA
      DOUBLE PRECISION, DIMENSION(KCON) :: SCALE,FGLOB
      DOUBLE PRECISION, DIMENSION(2,KCON) :: FHEM
      INTEGER, DIMENSION(JM,KCON) :: MLAT
      DOUBLE PRECISION, DIMENSION(JM+3,KCON) :: CNSLAT

      CHARACTER*4 :: HEMIS(2) = (/' SH ',' NH '/),DASH = ('----')

      CHARACTER*32 TITLE
      COMMON/D9COM/TITLE(KCON)

      INTEGER ::
     &     IHOUR,IHOUR0,INC,J,JHEMI,JNH,
     &     JP1,JPM,JSH,JV1,JVM,JX,N,NDAYS

      DOUBLE PRECISION ::
     &     AGLOB,AHEM,DTSRCE,FEQ,FNH,FSH,TAUDIF,XWON

C**** CALCULATE SCALING FACTORS
      XWON=TWOPI/(DLON*FIM)
      IF (IDACC(12).LT.1) IDACC(12)=1
      NDAYS=IDACC(9)/2
      DTSRCE=DT*NDYN
      SCALE(1)=100.D-9*RADIUS/(GRAV*IDACC(12))
      SCALE(2)=100.D-2*RADIUS/(GRAV*IDACC(6)*DTSRCE+1.D-20)
      SCALE(3)=SCALE(2)
      SCALE(4)=SCALE(2)
      SCALE(5)=SCALE(2)
      SCALE(6)=100.D-2*RADIUS/(DTSRCE*GRAV*IDACC(7)+1.D-20)
      SCALE(7)=100.D-2*RADIUS/(DTSRCE*GRAV*IDACC(8)+1.D-20)
      SCALE(8)=SCALE(7)
      SCALE(9)=100.D-2*RADIUS/(NFILTR*DT*GRAV*IDACC(10)+1.D-20)
      SCALE(10)=100.D-2*RADIUS/(SDAY*GRAV*NDAYS+1.D-20)
      SCALE(11)=1.
      SCALE(12)=25.D-3/(GRAV*IDACC(12))
      SCALE(13)=100.E3/(DTSRCE*GRAV*IDACC(6)+1.D-20)
      SCALE(14)=SCALE(13)
      SCALE(15)=SCALE(13)
      SCALE(16)=SCALE(13)
      SCALE(17)=25.E3/(DTSRCE*GRAV*IDACC(7)+1.D-20)
      SCALE(18)=25.E3/(DTSRCE*GRAV*IDACC(8)+1.D-20)
      SCALE(19)=SCALE(18)
      SCALE(20)=SCALE(18)
      SCALE(21)=25.E3/(NFILTR*DT*GRAV*IDACC(10)+1.D-20)
      SCALE(22)=25.E3/(SDAY*GRAV*NDAYS+1.D-20)
      SCALE(23)=1.
      SCALE(24)=100.E0/(GRAV*IDACC(12))
      SCALE(25)=100.E8/(DTSRCE*GRAV*IDACC(7)+1.D-20)
      SCALE(26)=100.E8/(NFILTR*DT*GRAV*IDACC(10)+1.D-20)
      SCALE(27)=100.E8/(SDAY*GRAV*NDAYS+1.D-20)
      SCALE(28)=1.
      SCALE(29)=100.D-5/(GRAV*IDACC(12))
      SCALE(30)=100.E2/(DTSRCE*GRAV*IDACC(7)+1.D-20)
      SCALE(31)=100.E2/(DTSRCE*GRAV*IDACC(8)+1.D-20)
      SCALE(32)=SCALE(31)
      SCALE(33)=SCALE(31)
      SCALE(34)=100.E2/(NFILTR*DT*GRAV*IDACC(10)+1.D-20)
      SCALE(35)=100.E2/(SDAY*GRAV*NDAYS+1.D-20)
      SCALE(36)=1.
      SCALE(37)=100.E2/(GRAV*IDACC(12))
      SCALE(38)=100.E8/(DTSRCE*GRAV*IDACC(7)+1.D-20)
      SCALE(39)=SCALE(38)
      SCALE(40)=SCALE(38)
      SCALE(41)=100.E8/(SDAY*GRAV*NDAYS+1.D-20)
      SCALE(42)=1.
C**** CALCULATE SUMMED QUANTITIES
      DO J=1,JM
         CONSRV(J,4)=CONSRV(J,2)+CONSRV(J,3)
         CONSRV(J,11)=CONSRV(J,6)*SCALE(6)+CONSRV(J,7)*SCALE(7)
     *        +CONSRV(J,8)*SCALE(8)+CONSRV(J,9)*SCALE(9)
     *        +CONSRV(J,10)*SCALE(10)
         CONSRV(J,15)=CONSRV(J,13)+CONSRV(J,14)
         CONSRV(J,23)=CONSRV(J,17)*SCALE(17)+CONSRV(J,18)*SCALE(18)
     *        +CONSRV(J,19)*SCALE(19)+CONSRV(J,20)*SCALE(20)
     *        +CONSRV(J,21)*SCALE(21)+CONSRV(J,22)*SCALE(22)
         CONSRV(J,28)=CONSRV(J,25)*SCALE(25)+CONSRV(J,26)*SCALE(26)
     *        +CONSRV(J,27)*SCALE(27)
         CONSRV(J,36)=CONSRV(J,30)*SCALE(30)+CONSRV(J,31)*SCALE(31)
     *        +CONSRV(J,32)*SCALE(32)+CONSRV(J,33)*SCALE(33)
     *        +CONSRV(J,34)*SCALE(34)+CONSRV(J,35)*SCALE(35)
         CONSRV(J,42)=CONSRV(J,38)*SCALE(38)+CONSRV(J,39)*SCALE(39)
     *        +CONSRV(J,40)*SCALE(40)+CONSRV(J,41)*SCALE(41)
      END DO
C**** CALCULATE FINAL ANGULAR MOMENTUM
      DO N=1,11
         FEQ=CONSRV(JEQ,N)*SCALE(N)*COSV(JEQ)
         FGLOB(N)=FEQ
         FHEM(1,N)=.5*FEQ
         FHEM(2,N)=.5*FEQ
         CNSLAT(JEQ,N)=FEQ/(FIM*DXYV(JEQ))
         DO JSH=2,JEQ-1
            JNH=2+JM-JSH
            FSH=CONSRV(JSH,N)*SCALE(N)*COSV(JSH)
            FNH=CONSRV(JNH,N)*SCALE(N)*COSV(JNH)
            FGLOB(N)=FGLOB(N)+(FSH+FNH)
            FHEM(1,N)=FHEM(1,N)+FSH
            FHEM(2,N)=FHEM(2,N)+FNH
            CNSLAT(JSH,N)=FSH/(FIM*DXYV(JSH))
            CNSLAT(JNH,N)=FNH/(FIM*DXYV(JNH))
         END DO
         FGLOB(N)=FGLOB(N)/AREAG
         FHEM(1,N)=FHEM(1,N)/(.5*AREAG)
         FHEM(2,N)=FHEM(2,N)/(.5*AREAG)
      END DO
C**** CALCULATE FINAL KINETIC ENERGY
      DO N=12,23
         FEQ=CONSRV(JEQ,N)*SCALE(N)
         FGLOB(N)=FEQ
         FHEM(1,N)=.5*FEQ
         FHEM(2,N)=.5*FEQ
         CNSLAT(JEQ,N)=FEQ/(FIM*DXYV(JEQ))
         DO JSH=2,JEQ-1
            JNH=2+JM-JSH
            FSH=CONSRV(JSH,N)*SCALE(N)
            FNH=CONSRV(JNH,N)*SCALE(N)
            FGLOB(N)=FGLOB(N)+(FSH+FNH)
            FHEM(1,N)=FHEM(1,N)+FSH
            FHEM(2,N)=FHEM(2,N)+FNH
            CNSLAT(JSH,N)=FSH/(FIM*DXYV(JSH))
            CNSLAT(JNH,N)=FNH/(FIM*DXYV(JNH))
         END DO
         FGLOB(N)=FGLOB(N)/AREAG
         FHEM(1,N)=FHEM(1,N)/(.5*AREAG)
         FHEM(2,N)=FHEM(2,N)/(.5*AREAG)
      END DO
C**** CALCUALTE FINAL MASS, TOTAL POTENTIAL ENERGY AND WATER
      DO N=24,KCON
         FGLOB(N)=0.
         FHEM(1,N)=0.
         FHEM(2,N)=0.
         DO JSH=1,JEQ-1
            JNH=1+JM-JSH
            FSH=CONSRV(JSH,N)*SCALE(N)
            FNH=CONSRV(JNH,N)*SCALE(N)
            FGLOB(N)=FGLOB(N)+(FSH+FNH)*DXYP(JSH)
            FHEM(1,N)=FHEM(1,N)+FSH*DXYP(JSH)
            FHEM(2,N)=FHEM(2,N)+FNH*DXYP(JNH)
            CNSLAT(JSH,N)=FSH/FIM
            CNSLAT(JNH,N)=FNH/FIM
         END DO
         FGLOB(N)=FGLOB(N)/AREAG
         FHEM(1,N)=FHEM(1,N)/(.5*AREAG)
         FHEM(2,N)=FHEM(2,N)/(.5*AREAG)
      END DO
      AGLOB=1.D-10*AREAG*XWON
      AHEM=1.D-10*(.5*AREAG)*XWON
C**** LOOP OVER HEMISPHERES
      INC=1+(JM-1)/24
      IHOUR0=TOFDY0+.5
      IHOUR=TOFDAY+.5
      TAUDIF=TAU-TAU0
      DO N=1,KCON
         DO J=1,JM
            MLAT(J,N)=NINT(CNSLAT(J,N))
         END DO
         CNSLAT(JM+1,N)=FHEM(1,N)
         CNSLAT(JM+2,N)=FHEM(2,N)
         CNSLAT(JM+3,N)=FGLOB(N)
      END DO
      DO JHEMI=2,1,-1
         WRITE (6,901) XLABEL
         WRITE (6,902) IDAY0,IHOUR0,JDATE0,JMNTH0,JYEAR0,TAU0,IDAY,IHOUR
     *        ,JDATE,JMONTH,JYEAR,TAU,TAUDIF
         JP1=1+(JHEMI-1)*(JEQ-1)
         JPM=JHEMI*(JEQ-1)
         JV1=2+(JHEMI-1)*(JEQ-2)
         JVM=JEQ+(JHEMI-1)*(JEQ-2)
C**** PRODUCE TABLES FOR ANGULAR MOMENTUM AND KINETIC ENERGY
         WRITE (6,903) (DASH,J=JV1,JVM,INC)
         WRITE (6,904) HEMIS(JHEMI),(JLAT(JX,2),JX=JVM,JV1,-INC)
         WRITE (6,903) (DASH,J=JV1,JVM,INC)
         DO N=1,23
            WRITE (6,905) TITLE(N),FGLOB(N),FHEM(JHEMI,N),
     *           (MLAT(JX,N),JX=JVM,JV1,-INC)
         END DO
         DO J=JV1,JVM
            MAREA(J)=1.D-10*XWON*FIM*DXYV(J)+.5
         END DO
         WRITE (6,906) AGLOB,AHEM,(MAREA(JX),JX=JVM,JV1,-INC)
C**** PRODUCE TABLES FOR MASS, TOTAL POTENTIAL ENERGY AND WATER
         WRITE (6,907)
         WRITE (6,903) (DASH,J=JP1,JPM,INC)
         WRITE (6,904) HEMIS(JHEMI),(JLAT(JX,1),JX=JPM,JP1,-INC)
         WRITE (6,903) (DASH,J=JP1,JPM,INC)
         DO N=24,KCON
            WRITE (6,905) TITLE(N),FGLOB(N),FHEM(JHEMI,N),
     *           (MLAT(JX,N),JX=JPM,JP1,-INC)
         END DO
         DO J=JP1,JPM
            MAREA(J)=1.D-10*XWON*FIM*DXYP(J)+.5
         END DO
         WRITE (6,906) AGLOB,AHEM,(MAREA(JX),JX=JPM,JP1,-INC)
      END DO
      RETURN
C****
  901 FORMAT ('1',A)
  902 FORMAT ('0CONSERVATION QUANTITIES   DAY',I6,', HR',I2,' (',I2,
     *  A5,I4,')',F9.0,'   TO   DAY',I6,', HR',I2,' (',I2,A5,I4,')',
     *  F9.0,'   DIF',F5.0,' HR'/)
  903 FORMAT (1X,25('--'),13(A4,'--'))
  904 FORMAT (35X,'GLOBAL',A7,2X,13I6)
  905 FORMAT (A32,2F9.2,1X,13I6)
  906 FORMAT ('0AREA (10**10 M**2)',F22.1,F9.1,1X,13I6)
  907 FORMAT ('0')
      END SUBROUTINE DIAG9P

      SUBROUTINE DIAG5D (M5,NDT,DUT,DVT)
      USE E001M12_COM, only : im,imh,jm,lm,fim,
     &     DSIG,JEQ,LS1,MDIAG,MDYN
      USE DAGCOM, only : speca,nspher
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: DUT,DVT

      DOUBLE PRECISION, DIMENSION(0:IMH,JM,LM,2) :: FCUVA,FCUVB
      COMMON/WORK7/FCUVA,FCUVB

      INTEGER :: M5,NDT

      DOUBLE PRECISION, DIMENSION(IMH+1) :: X
      DOUBLE PRECISION, DIMENSION(0:IMH) :: FA,FB
      DOUBLE PRECISION, DIMENSION(IMH+1,NSPHER) :: KE

      INTEGER :: MCLOCK ! external

      INTEGER ::
     &     J,J45N,KUV,KSPHER,L,
     &     MBEGIN,MEND,MKE,N,NM

      NM=1+IM/2
      J45N=2.+.75*(JM-1.)
      MKE=M5
C****
C**** KSPHER=1 SOUTHERN STRATOSPHERE       3 NORTHERN STRATOSPHERE
C****        2 SOUTHERN TROPOSPHERE        4 NORTHERN TROPOSPHERE
C****
C****        5 EQUATORIAL STRATOSPHERE     7 45 DEG NORTH STRATOSPHERE
C****        6 EQUATORIAL TROPOSPHERE      8 45 DEG NORTH TROPOSPHERE
C****
      GO TO (810,810,810,100,100,  100,810),M5
  810 WRITE (6,910) M5
  910 FORMAT ('0INCORRECT VALUE OF M5 WHEN CALLING DIAG5D.  M5=',I5)
      STOP 29
C****
C**** KINETIC ENERGY
C****
C**** TRANSFER RATES FOR KINETIC ENERGY IN THE DYNAMICS
  100 MBEGIN = MCLOCK ()
      KE(:,:)=0.

      DO 170 L=1,LM
      KSPHER=2
      IF (L.GE.LS1) KSPHER=1
      DO 170 J=2,JM
      DO 170 KUV=1,2 ! loop over u,v
      IF(KUV.EQ.1) CALL FFT(DUT(1,J,L),FA,FB)
      IF(KUV.EQ.2) CALL FFT(DVT(1,J,L),FA,FB)
      DO N=1,NM
         X(N)=.5*FIM*
     &        (FA(N-1)*FCUVA(N-1,J,L,KUV)+FB(N-1)*FCUVB(N-1,J,L,KUV))
      ENDDO
      X(1)=X(1)+X(1)
      X(NM)=X(NM)+X(NM)
      IF (J.NE.JEQ) KE(:,KSPHER)=KE(:,KSPHER)+X(:)*DSIG(L)
      IF (J.EQ.J45N) THEN     ! 45 N
         KE(:,KSPHER+4)=KE(:,KSPHER+4)+X(:)*DSIG(L)
      ELSE IF (J.EQ.JEQ) THEN ! EQUATOR
        DO N=1,NM
        KE(N,KSPHER+4)=KE(N,KSPHER+4)+     X(N)*DSIG(L)
        KE(N,KSPHER)  =KE(N,KSPHER)  +.5D0*X(N)*DSIG(L) ! CONTRIB TO SH
        KE(N,KSPHER+2)=KE(N,KSPHER+2)+.5D0*X(N)*DSIG(L) ! CONTRIB TO NH
        ENDDO
        IF (KUV.EQ.2) KSPHER=KSPHER+2
      ENDIF
  170 CONTINUE

      DO 180 KSPHER=1,NSPHER
      DO 180 N=1,NM
  180 SPECA(N,MKE,KSPHER)=SPECA(N,MKE,KSPHER)+KE(N,KSPHER)/NDT
      MEND = MCLOCK ()
      MDIAG = MDIAG + (MEND-MBEGIN)
      MDYN  = MDYN  - (MEND-MBEGIN)
      RETURN
      END SUBROUTINE DIAG5D

      SUBROUTINE DIAG5A (M5,NDT)
C****
C**** THIS DIAGNOSTICS ROUTINE PRODUCES A SPECTRAL ANALYSIS OF KINETIC
C**** AND AVAILABLE POTENTIAL ENERGIES AND THEIR TRANSFER RATES BY
C**** VARIOUS ATMOSPHERIC PROCESSES.
C****
C**** THE PARAMETER M INDICATES WHAT IS STORED IN SPECA(N,M,KSPHER),
C**** IT ALSO INDICATES WHEN DIAG5A IS BEING CALLED.
C**** M=1  MEAN STANDING KINETIC ENERGY            BEFORE SOURCES
C****   2  MEAN KINETIC ENERGY                     BEFORE DYNAMICS
C****   3  MEAN POTENTIAL ENERGY
C****   4  CONVERSION OF K.E. BY ADVECTION         AFTER ADVECTION
C****   5  CONVERSION OF K.E. BY CORIOLIS FORCE    AFTER CORIOLIS TERM
C****   6  CONVERSION FROM P.E. INTO K.E.          AFTER PRESS GRAD FORC
C****   7  CHANGE OF K.E. BY DYNAMICS              AFTER DYNAMICS
C****   8  CHANGE OF P.E. BY DYNAMICS
C****   9  CHANGE OF K.E. BY CONDENSATION          AFTER CONDENSATION
C****  10  CHANGE OF P.E. BY CONDENSATION
C****  11  CHANGE OF P.E. BY RADIATION             AFTER RADIATION
C****  12  CHANGE OF K.E. BY SURFACE               AFTER SURFACE
C****  13  CHANGE OF P.E. BY SURFACE
C****  14  CHANGE OF K.E. BY FILTER                AFTER FILTER
C****  15  CHANGE OF P.E. BY FILTER
C****  16  CHANGE OF K.E. BY DAILY                 AFTER DAILY
C****  17  CHANGE OF P.E. BY DAILY
C****  18  UNUSED
C****  19  LAST KINETIC ENERGY
C****  20  LAST POTENTIAL ENERGY
C****
      USE CONSTANT, only : kapa,sha
      USE E001M12_COM, only : im,imh,jm,lm,fim,
     &     DSIG,IDACC,JEQ,LS1,MDIAG,MDYN,
     &     P,PSF,PTOP,PSFMPT,SIG,T,U,V,ZATMO
      USE GEOM, only : AREAG,DXYN,DXYP,DXYS
      USE DAGCOM, only : speca,atpe,nspher,kspeca
      USE DYNAMICS, only : sqrtp,pk
      IMPLICIT NONE

      INTEGER :: M5,NDT

      DOUBLE PRECISION, DIMENSION(IM) :: X
      DOUBLE PRECISION, DIMENSION(IMH+1,NSPHER) :: KE,APE
      DOUBLE PRECISION, DIMENSION(IMH+1,4) :: VAR
      DOUBLE PRECISION, DIMENSION(2) :: TPE
      DOUBLE PRECISION, DIMENSION(IM,JM) :: SQRTM
      DOUBLE PRECISION, DIMENSION(LM) :: THJSP,THJNP,THGM

      INTEGER :: MCLOCK ! external

      INTEGER, PARAMETER :: IZERO=0

      INTEGER, DIMENSION(KSPECA), PARAMETER ::
     &     MTPEOF=(/0,0,1,0,0,0,0,2,0,3,  4,0,5,0,6,0,7,0,0,8/)

      INTEGER ::
     &     I,IJL2,IP1,J,J45N,
     &     JH,JHEMI,JP,K,KS,KSPHER,L,LDN,
     &     LUP,MAPE,MBEGIN,MEND,MINC,MKE,MNOW,MTPE,N,
     &     NM

      DOUBLE PRECISION ::
     &     GMEAN,GMSUM,SQRTPG,SUMI,SUMT,THGSUM,THJSUM

      INTEGER :: IFIRST
      DATA IFIRST/1/

      IF(IFIRST.NE.1) GO TO 50
      IFIRST=0
      SQRTPG = SQRT(PSFMPT)
      NM=1+IM/2
      J45N=2.+.75*(JM-1.)
      IJL2=IM*JM*LM*2
   50 CONTINUE

      MKE=M5
      MAPE=M5
C****
C**** KSPHER=1 SOUTHERN STRATOSPHERE       3 NORTHERN STRATOSPHERE
C****        2 SOUTHERN TROPOSPHERE        4 NORTHERN TROPOSPHERE
C****
C****        5 EQUATORIAL STRATOSPHERE     7 45 DEG NORTH STRATOSPHERE
C****        6 EQUATORIAL TROPOSPHERE      8 45 DEG NORTH TROPOSPHERE
C****
      GO TO (200,200,810,810,810,  810,200,810,205,810,
     *       296,205,810,205,810,  205,810,810,810,810),M5

  810 WRITE (6,910) M5
  910 FORMAT ('0INCORRECT VALUE OF M5 WHEN CALLING DIAG5A.  M5=',I5)
      STOP 29
C**** MASS FOR KINETIC ENERGY
  200 I=IM
      DO 202 J=2,JM
      DO 202 IP1=1,IM
      SQRTM(I,J)=SQRT(.5*((P(I,J)+P(IP1,J))*DXYS(J)+(P(I,J-1)+
     *  P(IP1,J-1))*DXYN(J-1)))
  202 I=IP1
C****
  205 MAPE=MKE+1
      KE(:,:)=0.
C**** CURRENT KINETIC ENERGY
      DO 240 L=1,LM
      KSPHER=2
      IF (L.GE.LS1) KSPHER=1
      DO 240 J=2,JM
      DO 240 K=IZERO,LM,LM
      IF(K.EQ.IZERO) X(1:IM)=U(1:IM,J,L)*SQRTM(1:IM,J)
      IF(K.EQ.LM)    X(1:IM)=V(1:IM,J,L)*SQRTM(1:IM,J)
      CALL FFTE (X,X)
      IF (J.EQ.JEQ) GO TO 225
      DO 220 N=1,NM
  220 KE(N,KSPHER)=KE(N,KSPHER)+X(N)*DSIG(L)
      IF (J.NE.J45N) GO TO 240
      DO 222 N=1,NM
  222 KE(N,KSPHER+4)=KE(N,KSPHER+4)+X(N)*DSIG(L)
      GO TO 240
  225 DO 230 N=1,NM
      KE(N,KSPHER+4)=KE(N,KSPHER+4)+X(N)*DSIG(L)
      KE(N,KSPHER)=KE(N,KSPHER)+.5D0*X(N)*DSIG(L)
  230 KE(N,KSPHER+2)=KE(N,KSPHER+2)+.5D0*X(N)*DSIG(L)
      IF (K.EQ.LM) KSPHER=KSPHER+2
  240 CONTINUE
      IF (NDT.EQ.0) GO TO 260
C**** TRANSFER RATES AS DIFFERENCES OF KINETIC ENERGY
      DO 250 KS=1,NSPHER
      DO 250 N=1,NM
  250 SPECA(N,MKE,KS)=SPECA(N,MKE,KS)+(KE(N,KS)-SPECA(N,19,KS))/NDT
  260 DO 270 KS=1,NSPHER
      DO 270 N=1,NM
  270 SPECA(N,19,KS)=KE(N,KS)
C****
C**** POTENTIAL ENERGY
C****
  296 CONTINUE
      APE(:,:)=0.
C**** CURRENT AVAILABLE POTENTIAL ENERGY
      LUP=0
  300 LUP=LUP+1
      THJSP(LUP)=T(1,1,LUP)*SQRTP(1,1)
      THJNP(LUP)=T(1,JM,LUP)*SQRTP(1,JM)
      IF(LUP.GE.LS1) THEN
      THJSP(LUP) = T(1,1,LUP)*SQRTPG
      THJNP(LUP) = T(1,JM,LUP)*SQRTPG
      ENDIF
      THGSUM=FIM*(THJSP(LUP)*DXYP(1)+THJNP(LUP)*DXYP(JM))
      DO 320 J=2,JM-1
      THJSUM=0.
      DO 310 I=1,IM
  310 THJSUM=THJSUM+T(I,J,LUP)*SQRTP(I,J)
  320 THGSUM=THGSUM+THJSUM*DXYP(J)
      THGM(LUP)=THGSUM/AREAG
      IF (LUP.GE.2) GO TO 350
      LDN=LUP
      L=LUP
      GO TO 300
  350 DO 360 JHEMI=1,2
      DO 360 N=2,NM
  360 VAR(N,JHEMI)=0.
      VAR(1,1)=.5*(THJSP(L)-THGM(L))**2*DXYP(1)*FIM
      VAR(1,2)=.5*(THJNP(L)-THGM(L))**2*DXYP(JM)*FIM
      GMEAN=((THJSP(LUP)-THJSP(LDN))*DXYP(1)*(SIG(L)*P(1,1)+PTOP)/
     *  (SQRTP(1,1)*P(1,1)*PK(L,1,1))+(THJNP(LUP)-THJNP(LDN))*DXYP(JM)*
     *  (SIG(L)*P(1,JM)+PTOP)/(SQRTP(1,JM)*P(1,JM)*PK(L,1,JM)))*FIM
      JHEMI=1
      DO 388 J=2,JM-1
      GMSUM=0.
      DO 370 I=1,IM
      X(I)=T(I,J,L)*SQRTP(I,J)-THGM(L)
  370 GMSUM=GMSUM+(T(I,J,LUP)-T(I,J,LDN))*(SIG(L)*P(I,J)+PTOP)/
     *  (P(I,J)*PK(L,I,J))
      GMEAN=GMEAN+GMSUM*DXYP(J)
      CALL FFTE (X,X)
      DO 380 N=1,NM
  380 VAR(N,JHEMI)=VAR(N,JHEMI)+X(N)*DXYP(J)
      IF (J.NE.JEQ-1) GO TO 384
      DO 382 N=1,NM
  382 VAR(N,3)=X(N)*DXYP(J)
      JHEMI=2
  384 IF (J.NE.J45N-1) GO TO 388
      DO 386 N=1,NM
  386 VAR(N,4)=X(N)*DXYP(J)
  388 CONTINUE
      GMEAN=DSIG(L)*AREAG*(SIG(LDN)-SIG(LUP))/GMEAN
      KS=2
      IF (L.GE.LS1) KS=1
      DO 400 JHEMI=1,4
      DO 390 N=1,NM
  390 APE(N,KS)=APE(N,KS)+VAR(N,JHEMI)*GMEAN
  400 KS=KS+2
      IF (L.EQ.LM) GO TO 450
      LDN=L
      L=LUP
      IF (LUP.LT.LM) GO TO 300
      GO TO 350
C**** CURRENT TOTAL POTENTIAL ENERGY
  450 DO 480 JHEMI=1,2
      JP=1+(JM-1.)*(JHEMI-1)
      SUMT=0.
      DO 455 L=1,LM
  455 SUMT=SUMT+T(1,JP,L)*PK(L,1,JP)*DSIG(L)
      TPE(JHEMI)=FIM*DXYP(JP)*(ZATMO(1,JP)*(P(1,JP)+PTOP)+
     *  SUMT*SHA*P(1,JP))
      DO 480 JH=2,JEQ-1
      J=JH+(JEQ-2)*(JHEMI-1)
      SUMI=0.
      DO 470 I=1,IM
      SUMT=0.
      DO 460 L=1,LM
  460 SUMT=SUMT+T(I,J,L)*PK(L,I,J)*DSIG(L)
  470 SUMI=SUMI+ZATMO(I,J)*(P(I,J)+PTOP)+SUMT*SHA*P(I,J)
  480 TPE(JHEMI)=TPE(JHEMI)+SUMI*DXYP(J)
      IF (NDT.EQ.0) GO TO 520
      MTPE=MTPEOF(MAPE)
C**** TRANSFER RATES AS DIFFERENCES FOR POTENTIAL ENERGY
      DO 510 KS=1,NSPHER
      DO 510 N=1,NM
  510 SPECA(N,MAPE,KS)=SPECA(N,MAPE,KS)+(APE(N,KS)-SPECA(N,20,KS))/NDT
      ATPE(MTPE,1)=ATPE(MTPE,1)+(TPE(1)-ATPE(8,1))/NDT
      ATPE(MTPE,2)=ATPE(MTPE,2)+(TPE(2)-ATPE(8,2))/NDT
  520 DO 530 KS=1,NSPHER
      DO 530 N=1,NM
  530 SPECA(N,20,KS)=APE(N,KS)
      ATPE(8,1)=TPE(1)
      ATPE(8,2)=TPE(2)
      CALL TIMER (MNOW,MINC,MDIAG)
      IF (M5.NE.2) RETURN
C**** ACCUMULATE MEAN KINETIC ENERGY AND MEAN POTENTIAL ENERGY
      IDACC(7)=IDACC(7)+1
      DO 550 KS=1,NSPHER
      DO 550 N=1,NM
      SPECA(N,2,KS)=SPECA(N,2,KS)+KE(N,KS)
  550 SPECA(N,3,KS)=SPECA(N,3,KS)+APE(N,KS)
      ATPE(1,1)=ATPE(1,1)+TPE(1)
      ATPE(1,2)=ATPE(1,2)+TPE(2)
      RETURN
      END SUBROUTINE DIAG5A
C****
      SUBROUTINE DIAG5F(UX,VX)
C**** FOURIER COEFFICIENTS FOR CURRENT WIND FIELD
C****
      USE E001M12_COM, only : im,imh,jm,lm,
     &     IDACC,MDIAG,MDYN
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: UX,VX

      DOUBLE PRECISION, DIMENSION(0:IMH,JM,LM,2) :: FCUVA,FCUVB
      COMMON/WORK7/FCUVA,FCUVB

      INTEGER :: MCLOCK ! external

      INTEGER ::
     &     J,L,MBEGIN,MEND,MINC

      MBEGIN = MCLOCK ()
      DO L=1,LM
         DO J=2,JM
            CALL FFT(UX(1,J,L),FCUVA(0,J,L,1),FCUVB(0,J,L,1))
            CALL FFT(VX(1,J,L),FCUVA(0,J,L,2),FCUVB(0,J,L,2))
         ENDDO
      ENDDO
      IDACC(6)=IDACC(6)+1
      MEND = MCLOCK ()
      MINC=MBEGIN-MEND
      MDIAG = MDIAG + (MEND-MBEGIN)
      MDYN  = MDYN  - (MEND-MBEGIN)
      RETURN
      END SUBROUTINE DIAG5F
C****
      SUBROUTINE DIAG5P
C**** THIS ENTRY PRINTS THE SPECTRAL ANALYSIS TABLES
C****
      USE CONSTANT, only : grav,rgas,twopi
      USE E001M12_COM, only : im,jm,lm,fim,
     &     DT,IDACC,IDAY,IDAY0,JDATE,JDATE0,JEQ,
     &     JMNTH0,JMONTH,JYEAR,JYEAR0,LS1,
     &     SKIPSE,TOFDAY,TOFDY0,XLABEL
      USE GEOM, only : DLON,DXYV
      USE DAGCOM, only : speca,atpe,ajk,aijk,kspeca,ktpe,nhemi,nspher
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(IM) :: X
      DOUBLE PRECISION, DIMENSION(KSPECA) :: SCALE,F0,FNSUM
      INTEGER, DIMENSION(KSPECA) :: MN

      DOUBLE PRECISION, DIMENSION(KTPE,NHEMI) :: FATPE

      INTEGER, PARAMETER :: IZERO=0

      INTEGER, DIMENSION(KTPE), PARAMETER ::
     &     MAPEOF=(/3,8,10,11,13,15,17,20/)

      CHARACTER*8 :: LATITD(4) = (/
     *     'SOUTHERN','NORTHERN',' EQUATOR','45 NORTH'/)
      CHARACTER*16 :: SPHERE(2) = (/
     *     'STRATOSPHERE    ','TROPOSPHERE     '/)

      INTEGER ::
     &     I,IHOUR,IHOUR0,IUNITJ,IUNITW,J,J45N,
     &     K,KPAGE,KROW,KSPHER,L,
     &     M,MAPE,MTPE,N,NM,NM1

      DOUBLE PRECISION ::
     &     FACTOR,FNM

      NM=1+IM/2
      IF (IDACC(12).LT.1) IDACC(12)=1
      IF (SKIPSE.GE.1.) GO TO 600
C      JEQ=1+JM/2
      J45N=2.+.75*(JM-1.)
C****
C**** STANDING KINETIC ENERGY
C****
      DO 710 K=1,NSPHER
      DO 710 N=1,NM
  710 SPECA(N,1,K)=0.
      DO 770 L=1,LM
      KSPHER=2
      IF (L.GE.LS1) KSPHER=1
      DO 770 J=2,JM
      IF (AJK(J,L,2).LE.1.D-20) GO TO 770
      FACTOR=FIM*DXYV(J)/AJK(J,L,2)
      DO 769 K=IZERO,LM,LM
      DO 720 I=1,IM
  720 X(I)=AIJK(I,J,L+K,1)
      CALL FFTE (X,X)
      IF (J.EQ.JEQ) GO TO 750
      DO 730 N=1,NM
  730 SPECA(N,1,KSPHER)=SPECA(N,1,KSPHER)+X(N)*FACTOR
      IF (J.NE.J45N) GO TO 769
      DO 740 N=1,NM
  740 SPECA(N,1,KSPHER+4)=SPECA(N,1,KSPHER+4)+X(N)*FACTOR
      GO TO 769
  750 DO 760 N=1,NM
      SPECA(N,1,KSPHER+4)=SPECA(N,1,KSPHER+4)+X(N)*FACTOR
      SPECA(N,1,KSPHER)=SPECA(N,1,KSPHER)+.5*X(N)*FACTOR
  760 SPECA(N,1,KSPHER+2)=SPECA(N,1,KSPHER+2)+.5*X(N)*FACTOR
      IF (K.EQ.LM) KSPHER=KSPHER+2
  769 CONTINUE
  770 CONTINUE
C****
  600 SCALE(1)=100.D-17/(GRAV*IDACC(4)+1.D-20)
      SCALE(19)=100.D-17/(GRAV*IDACC(12))
      SCALE(20)=SCALE(19)*RGAS
      SCALE(2)=SCALE(19)*IDACC(12)/(IDACC(7)+1.D-20)
      SCALE(3)=SCALE(2)*RGAS
      SCALE(4)=100.D-12/(GRAV*DT*IDACC(6)+1.D-20)
      SCALE(5)=SCALE(4)
      SCALE(6)=SCALE(4)
      SCALE(7)=100.D-12/(GRAV*DT*(IDACC(7)+1.D-20))
      SCALE(8)=SCALE(7)*RGAS
      SCALE(9)=100.D-12/(GRAV*DT*(IDACC(8)+1.D-20))
      SCALE(10)=SCALE(9)*RGAS
      SCALE(11)=SCALE(10)
      SCALE(12)=SCALE(9)
      SCALE(13)=SCALE(10)
      SCALE(14)=100.D-12/(GRAV*DT*(IDACC(10)+1.D-20))
      SCALE(15)=SCALE(14)*RGAS
      SCALE(16)=100.D-12/(GRAV*DT*(IDAY-IDAY0+1.D-20))
      SCALE(17)=SCALE(16)*RGAS
      SCALE(18)=100.D-17/(GRAV*IDACC(4)+1.D-20)
      DO 605 K=1,KSPECA
  605 SCALE(K)=(TWOPI/(DLON*FIM))*SCALE(K)
      IUNITJ=17
      IUNITW=12
      IHOUR0=TOFDY0+.5
      IHOUR=TOFDAY+.5
      DO 690 KPAGE=1,4
C**** WRITE HEADINGS
      WRITE (6,901) XLABEL
      WRITE (6,902) IDAY0,IHOUR0,JDATE0,JMNTH0,JYEAR0,IDAY,IHOUR,JDATE,
     *  JMONTH,JYEAR,IUNITJ,IUNITW
      DO 670 KROW=1,2
      IF (JM.GE.25.AND.KROW.EQ.2) WRITE (6,901)
      WRITE (6,903) LATITD(KPAGE),SPHERE(KROW)
      KSPHER=2*(KPAGE-1)+KROW
C**** WRITE KINETIC AND AVAILABLE POTENTIAL ENERGY BY WAVE NUMBER
      DO 610 M=1,KSPECA
      F0(M)=SPECA(1,M,KSPHER)*SCALE(M)
      MN(M)=NINT(F0(M))
  610 FNSUM(M)=0.
      WRITE (6,904) MN
      DO 630 N=2,NM
      KSPHER=2*(KPAGE-1)+KROW
      DO 620 M=1,KSPECA
      FNM=SPECA(N,M,KSPHER)*SCALE(M)
      MN(M)=NINT(FNM)
  620 FNSUM(M)=FNSUM(M)+FNM
      NM1=N-1
  630 WRITE (6,905) NM1,MN
      DO 640 M=1,KSPECA
  640 MN(M)=NINT(FNSUM(M))
      WRITE (6,906) MN
      DO 650 M=1,KSPECA
  650 MN(M)=NINT(FNSUM(M)+F0(M))
      WRITE (6,907) MN
  670 CONTINUE
      IF (KPAGE.GE.3) GO TO 690
C**** WRITE TOTAL POTENTIAL ENERGY
      DO 680 MTPE=1,KTPE
      MAPE=MAPEOF(MTPE)
         FATPE(MTPE,KPAGE)=ATPE(MTPE,KPAGE)*SCALE(MAPE)/RGAS
  680 MN(MTPE)=NINT(FATPE(MTPE,KPAGE))
      WRITE (6,909) (MN(MTPE),MTPE=1,8)
      IF (KPAGE.NE.2) GO TO 690
      DO 685 M=1,KSPECA
  685 SCALE(M)=SCALE(M)*10.
      IUNITJ=16
      IUNITW=11
  690 CONTINUE
      RETURN
C****
  901 FORMAT ('1',A)
  902 FORMAT ('0**  SPECTRAL ANALYSIS **  DAY',I6,', HR',I2,' (',I2,
     *  A5,I4,')   TO   DAY',I6,', HR',I2,' (',I2,A5,I4,
     *  ')   UNITS 10**',I2,' JOULES AND 10**',I2,' WATTS')
  903 FORMAT ('0',50X,A8,1X,A16/
     *  13X,'MEAN',19X,'DYNAMICS',25X,'SOURCES',16X,'FILTER',8X,
     *     'DAILY',4X,'PR SURF',5X,'LAST'/
     *'   N    SKE   KE   APE    KADV  KCOR   P-K  KDYN  PDYN   ',
     *     'KCNDS PCNDS   PRAD KSURF PSURF   KFIL  PFIL   KGMP  PGMP',
     *     '    KE',6X,'KE   APE')
  904 FORMAT ( '0  0',I7,I5,I6,I8,4I6,I8,I6,I7,2I6,I7,I6,I7,2I6,I8,I6/)
  905 FORMAT (     I4,I7,I5,I6,I8,4I6,I8,I6,I7,2I6,I7,I6,I7,2I6,I8,I6)
  906 FORMAT (' EDDY',I6,I5,I6,I8,4I6,I8,I6,I7,2I6,I7,I6,I7,2I6,I8,I6)
  907 FORMAT ('0TOTL',I6,I5,I6,I8,4I6,I8,I6,I7,2I6,I7,I6,I7,2I6,I8,I6)
  909 FORMAT (/'0TPE',I18,I32,I14,I7,I12,2I13,I20)
      END SUBROUTINE DIAG5P

      BLOCK DATA BDDLY
C****
C**** TITLES FOR SUBROUTINE DIAG6
C****
      COMMON/D6COM/TITLE
      CHARACTER*8 TITLE(63)
      DATA TITLE/
     *  '0INC SW ',' P ALBD ',' G ALBD ',' ABS ATM',' E CNDS ',
     *  '0SRF PRS',' PT 5   ',' PT 4   ',' PT 3   ',' PT 2   ',
     *  ' PT 1   ',' TS     ',' TG1    ','0Q 5    ',' Q 4    ',
     *  ' Q 3    ',' Q 2    ',' Q 1    ',' QS     ',' QG     ',
     *  '0CLD 6  ',' CLD 5  ',' CLD 4  ',' CLD 3  ',' CLD 2  ',
     *  ' CLD 1  ',' COVER  ','0SW ON G',' LW AT G',' SNSB HT',
     *  ' LAT HT ',' HEAT Z0','0UG*10  ',' VG*10  ',' WG*10  ',
     *  ' US*10  ',' VS*10  ',' WS*10  ',' ALPHA0 ','0RIS1*E2',
     *  ' RIGS*E2',' CDM*E4 ',' CDH*E4 ',' DGS*10 ',' EDS1*10',
     *  '0DBL    ',' DC FREQ',' LDC*10 ','0PRC*10 ',' EVP*10 ',
     *  ' DEEP MC',' SHLW MC',' CLD 7  ',' CLD 6  ',' CLD 5  ',
     *  ' CLD 4  ',' CLD 3  ',' CLD 2  ',' CLD 1  ',' W TO-5 ',
     *  ' C COVER',' SS P*10',' MC P*10'/
      END BLOCK DATA BDDLY

      SUBROUTINE DIAG6
C****
C**** THIS SUBROUTINE PRINTS THE DIURNAL CYCLE OF SOME QUANTITIES
C****
      USE CONSTANT, only : grav,rgas,kapa,sday,twopi
      USE E001M12_COM, only : dt,idacc,
     &     IJD6,JDATE,JDATE0,JMNTH0,JMONTH,JYEAR,JYEAR0,NAMD6,
     &     NCNDS,NDYN,NSURF,XLABEL
      USE DAGCOM, only : adaily,kdiag,ndlyvar,hr_in_day
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(NDLYVAR) :: SCALE
      DATA SCALE/1.,2*100.,2*1.,  5*1.,  3*1.,2*1.D5,  5*1.D5,
     *  5*100.,  2*100.,3*1.,  2*1.,3*10.,  3*10.,1.,100.,
     *  100.,2*1.D4,2*10.,  1.,100.,10.,2*1.,
     *  2*100., 7*100., 1.,100., 2*10./
      DOUBLE PRECISION, DIMENSION(HR_IN_DAY+1) :: XHOUR
      INTEGER, DIMENSION(HR_IN_DAY+1) :: MHOUR

      CHARACTER*8 TITLE
      COMMON/D6COM/TITLE(NDLYVAR)

      DOUBLE PRECISION :: AVE,AVED,AVEN,BYIDAC,DTCNDS,DTSURF
      INTEGER :: I,IH,IREGF,IREGL,IS,JY,JY0,K,KQ,KR,NDAYS

C****
      NDAYS=IDACC(9)/2
      IF (NDAYS.LE.0) RETURN
      DTCNDS=NCNDS*DT
      DTSURF=NDYN*DT/NSURF
      BYIDAC=1./NDAYS
      SCALE(5)=100.*RGAS/(KAPA*GRAV*DTCNDS)
      SCALE(28)=1./DTSURF
      SCALE(29)=1./DTSURF
      SCALE(30)=1./DTSURF
      SCALE(31)=1./DTSURF
      SCALE(32)=1./DTSURF
      SCALE(39)=360./TWOPI
      SCALE(49)=100.*100.*SDAY/(DTCNDS*GRAV)
      SCALE(62)=100.*100.*SDAY/(DTCNDS*GRAV)
      SCALE(63)=SCALE(62)
      SCALE(50)=100.*SDAY/DTSURF
C****
      IREGF=1
      IREGL=4
      IF (KDIAG(6).GT.0) IREGL=4-KDIAG(6)
      IF (KDIAG(6).LT.0.AND.KDIAG(6).GT.-5) IREGF=-KDIAG(6)
      IF (KDIAG(6).LT.0) IREGL=IREGF
      DO 500 KR=IREGF,IREGL
      JY0=JYEAR0-1900
      JY=JYEAR-1900
      WRITE (6,901) XLABEL(1:108),JDATE0,JMNTH0,JY0,JDATE,JMONTH,JY
      WRITE (6,903) NAMD6(KR),IJD6(1,KR),IJD6(2,KR),(I,I=1,HR_IN_DAY)
      DO 500 KQ=1,NDLYVAR
      IF (KQ.EQ.48) GO TO 200
      IF(KQ.GE.21.AND.KQ.LE.27) GO TO 500
C**** NORMAL QUANTITIES
      AVE=0.
      DO 120 IH=1,HR_IN_DAY
      AVE=AVE+ADAILY(IH,KQ,KR)
  120 XHOUR(IH)=ADAILY(IH,KQ,KR)*SCALE(KQ)*BYIDAC
      XHOUR(25)=AVE/FLOAT(HR_IN_DAY)*SCALE(KQ)*BYIDAC
      GO TO 480
C**** RATIO OF TWO QUANTITIES
  200 AVEN=0.
      AVED=0.
      DO 220 IH=1,HR_IN_DAY
      AVEN=AVEN+ADAILY(IH,KQ,KR)
      AVED=AVED+ADAILY(IH,KQ-1,KR)
  220 XHOUR(IH)=ADAILY(IH,KQ,KR)*SCALE(KQ)/(ADAILY(IH,KQ-1,KR)+1.D-20)
      XHOUR(25)=AVEN*SCALE(KQ)/(AVED+1.D-20)
  480 CONTINUE
      DO 490 IS=1,HR_IN_DAY+1
CB       FHOUR(IS,KQ,KR)=XHOUR(IS)
      MHOUR(IS)=NINT(XHOUR(IS))
  490 CONTINUE
      WRITE (6,904) TITLE(KQ),MHOUR
  500 CONTINUE
      RETURN
C****
  901 FORMAT ('1',A,I4,1X,A3,I3,' TO',I3,1X,A3,I3)
  903 FORMAT ('0',A4,I2,',',I2,' ',I2,23I5,'  AVE')
  904 FORMAT (A8,25I5)
      END SUBROUTINE DIAG6

      SUBROUTINE DIAG4A
C****
C**** THIS SUBROUTINE PRODUCES A TIME HISTORY OF ENERGIES
C****
      USE E001M12_COM, only : im,jm,lm,
     &     IDACC,JEQ,LS1,SKIPSE
      USE GEOM, only : DXYV
      USE DAGCOM, only : energy,speca,ajk,aijk
      IMPLICIT NONE

      INTEGER ::
     &     I,IDACC5,J,KS,L,N,NM
      DOUBLE PRECISION ::
     &     BYIADA,PU4TI,PV4TI,SEKE,SKE4I

      IF (IDACC(4).LE.0.OR.IDACC(7).LE.0) RETURN
C      JEQ=2.+.5*(JM-1.)
      NM=1+IM/2
C****
C**** LOAD ENERGIES INTO TIME HISTORY ARRAY
C****
      IDACC5=IDACC(5)+1
      IF (IDACC5.GT.100) RETURN
      IF (SKIPSE.EQ.1.) GO TO 540
C**** CALCULATE CURRENT SEKE
      BYIADA=1./IDACC(4)
      DO 530 L=1,LM
      KS=5
      IF (L.GE.LS1) KS=15
      DO 530 J=2,JM
      IF (AJK(J,L,2).LE.1.D-20) GO TO 530
      PU4TI=0.
      PV4TI=0.
      SKE4I=0.
      DO 510 I=1,IM
      PU4TI=PU4TI+AIJK(I,J,L,1)
      PV4TI=PV4TI+AIJK(I,J,L,2)
  510 SKE4I=SKE4I+(AIJK(I,J,L,1)*AIJK(I,J,L,1)
     *            +AIJK(I,J,L,2)*AIJK(I,J,L,2))/(AIJK(I,J,L,4)+1.D-20)
      SEKE=(SKE4I-(PU4TI*PU4TI+PV4TI*PV4TI)/AJK(J,L,2))*DXYV(J)*BYIADA
      IF (J.EQ.JEQ) GO TO 520
      ENERGY(KS,IDACC5)=ENERGY(KS,IDACC5)+SEKE
      GO TO 530
  520 ENERGY(KS,IDACC5)=ENERGY(KS,IDACC5)+.5*SEKE
      ENERGY(KS+1,IDACC5)=ENERGY(KS+1,IDACC5)+.5*SEKE
      KS=KS+1
  530 CONTINUE
C**** OTHER ENERGIES COME FROM LATEST SPECTRAL ANALYSIS
  540 ENERGY(1,IDACC5)=SPECA(1,19,2)
      ENERGY(2,IDACC5)=SPECA(1,19,4)
      ENERGY(7,IDACC5)=SPECA(1,20,2)
      ENERGY(8,IDACC5)=SPECA(1,20,4)
      ENERGY(11,IDACC5)=SPECA(1,19,1)
      ENERGY(12,IDACC5)=SPECA(1,19,3)
      ENERGY(17,IDACC5)=SPECA(1,20,1)
      ENERGY(18,IDACC5)=SPECA(1,20,3)
      DO 550 N=2,NM
      ENERGY(3,IDACC5)=ENERGY(3,IDACC5)+SPECA(N,19,2)
      ENERGY(4,IDACC5)=ENERGY(4,IDACC5)+SPECA(N,19,4)
      ENERGY(9,IDACC5)=ENERGY(9,IDACC5)+SPECA(N,20,2)
      ENERGY(10,IDACC5)=ENERGY(10,IDACC5)+SPECA(N,20,4)
      ENERGY(13,IDACC5)=ENERGY(13,IDACC5)+SPECA(N,19,1)
      ENERGY(14,IDACC5)=ENERGY(14,IDACC5)+SPECA(N,19,3)
      ENERGY(19,IDACC5)=ENERGY(19,IDACC5)+SPECA(N,20,1)
  550 ENERGY(20,IDACC5)=ENERGY(20,IDACC5)+SPECA(N,20,3)
      IDACC(5)=IDACC5
      RETURN
      END SUBROUTINE DIAG4A
C****
      SUBROUTINE DIAG4
C**** THIS ENTRY PRODUCES A TIME HISTORY TABLE OF ENERGIES
C****
      USE CONSTANT, only : grav,rgas,twopi
      USE E001M12_COM, only : im,jm,lm,fim,
     &     DT,IDACC,IDAY,IDAY0,JDATE,JDATE0,JMNTH0,JMONTH,
     &     JYEAR,JYEAR0,NDA4,NDYN,TAU0,TOFDAY,TOFDY0,XLABEL
      USE GEOM, only : DLON
      USE DAGCOM, only : energy,nehist,hist_days
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(NEHIST) :: SUM,SCALE
      INTEGER, DIMENSION(NEHIST) :: IK
      DOUBLE PRECISION, DIMENSION(NEHIST,HIST_DAYS+1) :: EHIST

      INTEGER ::
     &     I,IDACC5,IDAYX,IDAYXM,IHOUR,IHOUR0,K
      DOUBLE PRECISION ::
     &     TAUX,TOFDYX

      IDACC5=IDACC(5)
      IF (IDACC5.LE.0) RETURN
      IF (IDACC(12).LT.1) IDACC(12)=1
      SCALE(1)=100.D-18/GRAV
      SCALE(2)=SCALE(1)
      SCALE(3)=SCALE(1)
      SCALE(4)=SCALE(1)
      SCALE(5)=.5*SCALE(1)
      SCALE(6)=SCALE(5)
      SCALE(7)=SCALE(1)*RGAS
      SCALE(8)=SCALE(7)
      SCALE(9)=SCALE(7)
      SCALE(10)=SCALE(7)
      SCALE(11)=SCALE(1)
      SCALE(12)=SCALE(1)
      SCALE(13)=SCALE(1)
      SCALE(14)=SCALE(1)
      SCALE(15)=SCALE(5)
      SCALE(16)=SCALE(5)
      SCALE(17)=SCALE(7)
      SCALE(18)=SCALE(7)
      SCALE(19)=SCALE(7)
      SCALE(20)=SCALE(7)
      DO 60 K=1,NEHIST
   60 SCALE(K)=(TWOPI/(DLON*FIM))*SCALE(K)/IDACC(12)
C****
      IHOUR0=TOFDY0+.5
      IHOUR=TOFDAY+.5
      WRITE (6,901) XLABEL
      WRITE (6,902) IDAY0,IHOUR0,JDATE0,JMNTH0,JYEAR0,IDAY,IHOUR,JDATE,
     *  JMONTH,JYEAR
      DO 110 K=1,NEHIST
  110 SUM(K)=0.
      WRITE (6,903)
      DO 200 I=1,IDACC5
      TAUX=TAU0+(I*NDA4-NDYN)*DT/3600.
      IDAYX=1.+(TAUX+.001)/24.
      IDAYXM=MOD(IDAYX,100000)
      TOFDYX=TAUX-24.*(IDAYX-1)
      DO 150 K=1,NEHIST
         EHIST(K,I)=ENERGY(K,I)*SCALE(K)
      IK(K)=EHIST(K,I)+.5
  150 SUM(K)=SUM(K)+ENERGY(K,I)
      WRITE (6,904) IDAYXM,TOFDYX,IK
  200 CONTINUE
      DO 250 K=1,NEHIST
         EHIST(K,101)=SUM(K)*SCALE(K)/IDACC5
  250 IK(K)=EHIST(K,101)+.5
      WRITE (6,905) IK
         CALL KEYD4 (IK)
      RETURN
C****
  901 FORMAT ('1',A)
  902 FORMAT ('0** ENERGY HISTORY **   DAY',I6,', HR',I3,' (',I2,A5,I5,
     *  ')    TO    DAY',I6,', HR',I3,' (',I2,A5,I5,
     *  ')    UNITS OF 10**18 JOULES')
  903 FORMAT ('0',15X,21('-'),' TROPOSPHERE ',22('-'),5X,21('-'),
     *  ' STRATOSPHERE ',21('-')/8X,2(11X,'ZKE',8X,'EKE',7X,'SEKE',9X,
     * 'ZPE',10X,'EPE')/3X,'DAY  HOUR     SH   NH    SH   NH    SH   NH
     *    SH    NH     SH    NH      SH   NH    SH   NH    SH   NH     S
     *H    NH     SH    NH'/1X,132('='))
  904 FORMAT (I6,F6.1,1X,3(I6,I5),2(I7,I6),2X,3(I6,I5),2(I7,I6))
  905 FORMAT (1X,132('=')/8X,'MEAN ',3(I6,I5),2(I7,I6),2X,3(I6,I5),
     *  2(I7,I6))
      END SUBROUTINE DIAG4

      SUBROUTINE DIAGKS
C****
C**** THIS SUBROUTINE PRODUCES A SUMMARY OF KEY NUMBERS CALCULATED IN
C**** OTHER DIAGNOSTIC SUBROUTINES
C****
C**** CONTENTS OF KEYNR
C****
C   K   N
C****
C***1  1 MONTH
C***2  2 TOTAL CLOUD COVER (PERCENT)
C****  3 SNOW COVER--NORTHERN HEMSIPHERE (PERCENT)
C****  4 ICE COVER--NORTHERN HEMISPHERE (PERCENT)
C****  5 PLANETARY ALBEDO (PERCENT)
C****  6 SOLAR RADIATION ABSORBED BY ATMOSPHERE (WT/M**2)
C****  7 SOLAR RADIATION ABSORBED BY PLANET (WT/M**2)
C****  8 NET HEAT AT GROUND (WT/M**2)
C****  8 ANGULAR MOMENTUM PER UNIT AREA (10**10 J*SEC/M**2)
C****  9 EVAPORATION (.1 MM/DAY)
C****  9 PRECIPITATION (.1 MM/DAY)
C**** 10 SENSIBLE HEAT FLUX INTO GROUND (ABS.VALUE)
C**** 11 LATENT HEAT FLUX INTO GROUND (ABS.VALUE)
C**** 12 MEAN GROUND TEMPERATURE (DEGREES K)
C**** 13 MEAN GLOBAL ATMOSPHERIC TEMPERATURE (DEGREES K)
C**** 14 MERID. TEMPERATURE GRADIENT (N.HEMISPHERE)
C**** 15 MERID. TEMPERATURE GRADIENT (S.HEMISPHERE)
C**** 16 MEAN TROPOSPHERIC EKE-NORTHERN HEMISPHERE
C**** 17 MEAN TROPOSPHERIC EKE-SOUTHERNN HEMISPHERE
C**** 18 MEAN TROPOSPHERIC ZKE-NORTHERN HEMISPHERE
C**** 19 MEAN TROPOSPHERIC ZKE-SOUTHERN HEMISPHERE
C**** 20 MEAN TROPOSPHERIC EPE-NORTHERN HEMISPHERE
C**** 21 MEAN TROPOSPHERIC ZPE-NORTHERN HEMISPHERE
C**** 22 MEAN EDDY KINETIC ENERGY AT EQUATOR
C**** 23 MAX. MEAN EDDY KINETIC ENERGY IN MID NORTH LATITUDES
C**** 24 MAX. ZONAL WIND (U COMPONENT) IN TROPOSPHERE (NH), M/SEC
C**** 25 LATITUDE CORRESPONDING TO 24
C**** 26 MAX. ZONAL WIND (U COMPONENT) IN TROPOSPHERE (SH), M/SEC
C**** 27 LATITUDE CORRESPONDING TO 26
C**** 28-30: 29 IS LARGEST VALUE OF STREAM FUNCTION, POSITIVE OR
C****    NEGATIVE; 28 AND 30 ARE THE MAGNITUDES OF THE LARGEST VALUES OF
C****    OPPOSITE SIGN TO THE NORTH AND SOUTH RESPECTIVELY
C***3 31 SNOW AND ICE COVERAGE OF GLOBE (PERCENT)
C***4 32 SNOW AND ICE COVERAGE OF NORTHERN HEMISPHERE (PERCENT)
C****   33-39 REFER TO NORTHERN HEMISPHERE ONLY
C**** 33 MAX.NORTHWARD TRANS. OF DRY STATIC ENERGY BY STANDING EDDIES
C**** 34 MAX.NORTHWARD TRANS. OF DRY STATIC ENERGY BY EDDIES
C**** 35 MAX. TOTAL NORTH. TRANS. OF DRY STATIC ENERGY
C**** 36 MAX.NORTHWARD TRANS. OF STATIC ENERGY BY EDDIES
C**** 37 MAX.TOTAL NORTH. TRANS. OF STATIC ENERGY
C**** 38 LATITUDE CORRESPONDING TO 37
C**** 39 MAX. NORTH. TRANS. OF ANGULAR MOMENTUM BY STANDING EDDIES
C**** 40 MAX. NORTH. TRANS. OF ANGULAR MOMENTUM BY EDDIES
C**** 41 MAX. TOTAL NORTH. TRANS. OF ANGULAR MOMENTUM
C**** 42 LATITUDE CORRESPONDING TO 41
C****
      USE CONSTANT, only : twopi
      USE E001M12_COM, only : jm,lm,
     &     DT,IDAY,IDAY0,JDATE,JDATE0,JEQ,JMNTH0,JMONTH,JYEAR,JYEAR0,
     &     KEYCT,NDYN,PSF,PTOP,SIG,TAU,TAU0,TAUI,TOFDAY,TOFDY0,XLABEL,
     &     PSFMPT
      USE GEOM, only : DLAT,DXYP,JLAT
      USE DAGCOM, only : keynr,nehist,nkeynr
      IMPLICIT NONE

      integer, parameter :: int_ic=1229135904
c      = 256*(256*(256*iachar('I')+iachar('C'))+iachar(' '))+iachar(' ')
      integer, parameter :: int_jan=1245793824
c      = 256*(256*(256*iachar('J')+iachar('A'))+iachar('N'))+iachar(' ')

      INTEGER*4 :: NDEX(42) = (/
     *     1,2,31,32,3,   4,5,6,7,8,   9,10,11,12,13,
     *             14,15,16,17,18,  19,20,21,22,23,  24,25,26,27,28,
     *             29,30,33,34,35,  36,37,38,39,40,  41,42/)

      DOUBLE PRECISION, DIMENSION(JM) :: FLAT
      DOUBLE PRECISION, DIMENSION(JM,LM) :: FKEY
      DOUBLE PRECISION, DIMENSION(JM) :: ASUM
      DOUBLE PRECISION, DIMENSION(2) :: HSUM
      INTEGER, DIMENSION(NEHIST) :: IK

      INTEGER ::
     &     I,I35,I70,IEND,ISIGN,IHOUR,IHOUR0,INT_JMNTH0,
     &     J,J60,JMAX,JEQP1,JNDEX,JSTART,
     &     K,KEYMAX,KNDEX,
     &     L,LL,LMAX,LNLM,LNM,LSLM,LSM,
     &     N,NT,NTDIF

      DOUBLE PRECISION ::
     &     A,BIG,CEPT,CHECK,
     &     FGLOB,FNH,GSUM,HN,HS,PISG,PISN,
     &     SAVE,SUMFAC,TAUDIF,TEQ,TNOR,TSOU,
     &     UNLM,UNM,USLM,USM,X60

C****
C**** ENTRIES CALLED FROM DIAGJ
C****
      ENTRY KEYDJ (N,FGLOB,FNH)
      GO TO (                100,100,100,110,100, 100,100,100,100,115,
     *  100,100,120,125,100, 100,100,130,100,135, 100,100,100,100,100,
     *  100,100,100,100,140, 145,100,100,100,100, 100,100,100,100,100,
     *  100,100,100,150,100, 100,100,100,100,100, 100,100,100,100,100,
     *  100,100,100,155),N
  100 RETURN
  110 KEYNR(6,KEYCT)=NINT(FGLOB)
      RETURN
  115 KEYNR(7,KEYCT)=NINT(FGLOB)
      RETURN
  120 KEYNR(10,KEYCT)=NINT(-FGLOB)
      RETURN
  125 KEYNR(11,KEYCT)=NINT(-FGLOB)
      RETURN
  130 KEYNR(12,KEYCT)=NINT(.1*FGLOB)
      RETURN
  135 KEYNR(9,KEYCT)=NINT(10.*FGLOB)
      RETURN
  140 KEYNR(4,KEYCT)=NINT(FNH)
      RETURN
  145 KEYNR(3,KEYCT)=NINT(FNH)
      RETURN
  150 KEYNR(8,KEYCT)=NINT(FGLOB)
      RETURN
  155 KEYNR(2,KEYCT)=NINT(FGLOB)
      RETURN
C****
      ENTRY KEYDJA (FGLOB)
      KEYNR(5,KEYCT)=NINT(10.*FGLOB)
      RETURN
C****
C**** ENTRIES CALLED FROM DIAGJL VIA JLMAP OR FROM DIAGJK VIA JKMAP
C****
      ENTRY KEYJKT (GSUM,ASUM)
C**** TEMPERATURES
C      JEQ=2.+.5*(JM-1.)
      TEQ=.5*(ASUM(JEQ-1)+ASUM(JEQ))
      X60=TWOPI/(12.*DLAT)
      J60=.5+X60
      A=DXYP(J60+1)*(X60+.5-J60)
      TSOU=ASUM(J60+1)*A
      TNOR=ASUM(JM-J60)*A
      DO 210 J=1,J60
      A=A+DXYP(J)
      TSOU=TSOU+ASUM(J)*DXYP(J)
  210 TNOR=TNOR+ASUM(JM+1-J)*DXYP(J)
      KEYNR(14,KEYCT)=NINT(TEQ-TNOR/A)
      KEYNR(15,KEYCT)=NINT(TEQ-TSOU/A)
      KEYNR(13,KEYCT)=NINT(.1*GSUM)
      RETURN
C****
      ENTRY KEYJKJ (L,FLAT)
C**** JET STREAMS
      IF (L.LT.LM) GO TO 220
      DO 216 LL=1,LM
      IF (PSFMPT*SIG(LL)+PTOP.LT.200.) GO TO 218
  216 CONTINUE
  218 LMAX=LL-1
  220 IF (L.GT.LMAX) RETURN
      USLM=-999999.
      DO 222 J=3,JEQ
      IF (FLAT(J).LT.USLM) GO TO 222
      USLM=FLAT(J)
      JMAX=J
  222 CONTINUE
      CEPT=.5*(FLAT(JMAX-1)-FLAT(JMAX+1))/
     *  (FLAT(JMAX-1)-2.*FLAT(JMAX)+FLAT(JMAX+1))
      LSLM=INT((JMAX-1.5+CEPT)*DLAT*360/TWOPI+.5)-90
      UNLM=-999999.
      DO 224 J=JEQ,JM-1
      IF (FLAT(J).LT.UNLM) GO TO 224
      UNLM=FLAT(J)
      JMAX=J
  224 CONTINUE
      CEPT=.5*(FLAT(JMAX-1)-FLAT(JMAX+1))/
     *  (FLAT(JMAX-1)-2.*FLAT(JMAX)+FLAT(JMAX+1))
      LNLM=INT((JMAX-1.5+CEPT)*DLAT*360/TWOPI+.5)-90
      IF (L.LT.LMAX) GO TO 226
      USM=USLM
      LSM=LSLM
      UNM=UNLM
      LNM=LNLM
      RETURN
  226 IF (USLM.LT.USM) GO TO 228
      USM=USLM
      LSM=LSLM
  228 IF (UNLM.LT.UNM) GO TO 230
      UNM=UNLM
      LNM=LNLM
  230 IF (L.NE.1) RETURN
      KEYNR(24,KEYCT)=.1*UNM+.5
      KEYNR(25,KEYCT)=LNM
      KEYNR(26,KEYCT)=.1*USM+.5
      KEYNR(27,KEYCT)=-LSM
      RETURN
C****
      ENTRY KEYJLS (L,FLAT)
C**** STREAM FUNCTION
      DO 290 J=2,JM
  290 FKEY(J,L)=FLAT(J)
      IF (L.NE.1) RETURN
  300 SAVE=0.
      HS=0.
      HN=0.
      DO 310 K=1,LM
      DO 310 I=2,JM
      CHECK=ABS(FKEY(I,K))
      IF (CHECK.LT.SAVE) GO TO 310
      SAVE=CHECK
      JNDEX=I
      KNDEX=K
  310 CONTINUE
      SAVE=FKEY(JNDEX,KNDEX)
      ISIGN=1
      IF (SAVE.GT.0.0) ISIGN=-1
      IF (JNDEX.LT.4) GO TO 325
      IEND=JNDEX-1
      DO 320 K=1,LM
      DO 320 I=2,IEND
      CHECK=FKEY(I,K)*ISIGN
  320 IF (CHECK.GT.HS)HS=CHECK
  325 CONTINUE
      IF (JNDEX.GT.(JM-2))GO TO 335
      JSTART=JNDEX+1
      DO 330 K=1,LM
      DO 330 I=JSTART,JM
      CHECK=FKEY(I,K)*ISIGN
  330 IF (CHECK.GT.HN)HN=CHECK
  335 CONTINUE
      KEYNR(28,KEYCT)=ABS(HN)+0.5
      KEYNR(29,KEYCT)=NINT(SAVE)
      KEYNR(30,KEYCT)=ABS(HS)+0.5
      RETURN
C****
      ENTRY KEYJKE (NT,HSUM,ASUM)
C**** EDDY AND ZONAL KINETIC ENERGY
      IF (NT.EQ.19) GO TO 450
      KEYNR(16,KEYCT)=NINT(HSUM(2))
      KEYNR(17,KEYCT)=NINT(HSUM(1))
      KEYNR(18,KEYCT)=KEYNR(18,KEYCT)-NINT(HSUM(2))
      KEYNR(19,KEYCT)=KEYNR(19,KEYCT)-NINT(HSUM(1))
      KEYNR(22,KEYCT)=NINT(ASUM(JEQ))
      BIG=-99999.
      I35=2.+(JM-1.)*125./180.
      I70=2.+(JM-1.)*160./180.
      DO 440 I=I35,I70
      IF (ASUM(I).LT.BIG) GO TO 440
      BIG=ASUM(I)
  440 CONTINUE
      KEYNR(23,KEYCT)=NINT(BIG)
      RETURN
  450 KEYNR(18,KEYCT)=KEYNR(18,KEYCT)+NINT(HSUM(2))
      KEYNR(19,KEYCT)=KEYNR(19,KEYCT)+NINT(HSUM(1))
      RETURN
C****
      ENTRY KEYJKN (NT,ASUM,SUMFAC)
C**** NORTHWARD TRANSPORTS
  500 BIG=-99999.
      JEQP1=JEQ+1
      DO 510 I=JEQP1,JM
      IF (ASUM(I).LT.BIG) GO TO 510
      BIG=ASUM(I)
      JNDEX=I
  510 CONTINUE
      BIG=BIG*SUMFAC
      NTDIF=NT-21
      GO TO (392,392,392,390,390,396,394,390,390,400,400,398),NTDIF
  390 CONTINUE
  392 KEYNR(NT+11,KEYCT)=NINT(BIG)
      RETURN
  394 KEYNR(38,KEYCT)=JLAT(JNDEX,2)
  396 KEYNR(NT+9,KEYCT)=NINT(BIG)
      RETURN
  398 KEYNR(42,KEYCT)=JLAT(JNDEX,2)
  400 KEYNR(NT+8,KEYCT)=NINT(BIG)
      RETURN
C****
C**** ENTRY CALLED FROM DIAGIJ
C****
      ENTRY KEYIJ(PISG,PISN)
      KEYNR(31,KEYCT)=NINT(PISG)
      KEYNR(32,KEYCT)=NINT(PISN)
      RETURN
C****
C**** ENTRY CALLED FROM DIAG4
C****
      ENTRY KEYD4 (IK)
      KEYNR(20,KEYCT)=(IK(10)+IK(20)+5)/10
      KEYNR(21,KEYCT)=(IK(8)+IK(18)+5)/10
      RETURN
C****
      ENTRY DIAGKN
C**** PRINTS THE TABLE OF KEY NUMBERS
C****
      IHOUR0=TOFDY0+.5
      IHOUR=TOFDAY+.5
      TAUDIF=TAU-TAU0
      int_jmnth0 = 256*(256*(256*
     &     iachar(jmnth0(1:1))  +iachar(jmnth0(2:2)))+
     &     iachar(jmnth0(3:3))) +iachar(jmnth0(4:4))
      KEYNR(1,KEYCT)=int_jmnth0
      IF (TAU.LE.TAUI+(DT/3600.)*(NDYN+.5)) KEYNR(1,KEYCT)=INT_IC
      IF (KEYCT.GE.2.AND.KEYNR(1,KEYCT-1).EQ.INT_JMNTH0) KEYCT=KEYCT-1
      WRITE(6,901) XLABEL
      WRITE(6,910) IDAY0,IHOUR0,JDATE0,JMNTH0,JYEAR0,TAU0,IDAY,IHOUR,
     *  JDATE,JMONTH,JYEAR,TAU,TAUDIF
      WRITE(6,902)
      DO 810 I=1,KEYCT
      IF (KEYNR(1,I).EQ.INT_JAN) WRITE (6,905)
  810 WRITE(6,905) (KEYNR(NDEX(K),I),K=1,42)
      WRITE (6,915)
CB       DO 815 K=1,42
CB815    FKEYDS(K)=KEYNR(K,KEYCT)
      KEYCT=KEYCT+1
      KEYMAX=49
      IF (KEYNR(1,1).NE.INT_IC) KEYMAX=48
      IF (KEYCT.LE.KEYMAX) RETURN
C**** ROLL UP KEY NUMBERS 1 YEAR AT A TIME
      DO 820 K=1,36
      DO 820 I=1,NKEYNR
  820 KEYNR(I,K)=KEYNR(I,K+KEYMAX-36)
      DO 880 K=37,50
      DO 880 I=1,NKEYNR
  880 KEYNR(I,K)=0
      KEYCT=37
      RETURN
  901 FORMAT('1',A)
  902 FORMAT ('0',7X,'SN+IC NH NH AL AB NT NT PR        T   T-OF-ATM  EK
     *E   ZKE           EKE   JET-STREAMS STREAM-FN NOR-TRAN NOR-TRAN NO
     *RTH-TRANS'/
     *         5X,'CL GL    SN OI BE BY RD HT EC SN LAT OF  GL  GRAD ---
     *-- ----- EPE ZPE ------ NORTH SOUTH --------- DRY-STAT STAT-ENR AN
     *G MOMENTM'/
     *         5X,'CV OB NH CV CV DO AT P0 Z0 IP HT  HT GD  OB NH SH NH
     *SH NH SH  NH  NH EQ  ML VL LT VL LT NH MAX SH SE ED TL ED TL LT SE
     * ED TL LT'/)
  905 FORMAT (1X,A3,4I3,I2,I4,5I3,I4,I3,I4,6I3,2I4,I3,I4,5I3,I4,11I3)
  910 FORMAT ('0',13X,'DAY',I6,', HR',I3,' (',I2,A5,I5,')',F9.1,
     *  '   TO   DAY',I6,', HR',I3,' (',I2,A5,I5,')',F9.1,'   DIF',
     *  F6.1,' HR',7X,I5,I5)
  915 FORMAT('0')
      END

      SUBROUTINE DIAG10(IPFLAG)
C****
C**** THIS SUBROUTINE SAVES THE INSTANTANEOUS SEA LEVEL PRESSURES
C**** EVERY ABS(USESLP) HOURS. IF USESLP.LT.0 THE FIRST RECORD IS
C**** WRITTEN TO THE BEGINNING OF UNIT 16.
C****
      USE CONSTANT, only : grav,rgas
      USE E001M12_COM, only : im,jm,p,ptop,tau,zatmo
      USE PBLCOM, only : tsavg
      IMPLICIT NONE
      INTEGER :: IPFLAG
      DOUBLE PRECISION, DIMENSION(IM,JM) :: SLP
      INTEGER :: I,J
      DOUBLE PRECISION :: BBYG,GBYRB
      BBYG=.0065/GRAV
      GBYRB=GRAV/(.0065*RGAS)
      DO 10 J=1,JM
      DO 10 I=1,IM
   10 SLP(I,J)=(P(I,J)+PTOP)*(1.+BBYG*ZATMO(I,J)/TSAVG(I,J))**GBYRB
      WRITE (16) SNGL(TAU),((SNGL(SLP(I,J)),I=1,IM),J=1,JM),SNGL(TAU)
      ENDFILE 16
      BACKSPACE 16
      RETURN
      ENTRY DIAG8(IPFLAG)
      RETURN
      END SUBROUTINE DIAG10

      SUBROUTINE init_DIAG
!@sum  init_DIAG initiallises the diagnostics
!@auth Gavin Schmidt
!@ver  1.0
      USE E001M12_COM
      USE DAGCOM
      USE DAGPCOM
      IMPLICIT NONE
      INTEGER L

C**** AIJ diagnostic names:
C**** NAME     NO.    DESCRIPTION   (SCALE)*IDACC  LOCATION
C**** ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IJ_RSOI = 1  ! POICE (1)            1 GD
      IJ_RSNW = 2  ! PSNOW (1)            4 DA
      IJ_SNOW = 3  ! SNOW (KG/M**2)       4 DA
      IJ_SHDT = 4  ! SHDT (J/M**2)        1 SF
      IJ_PREC = 5  ! PREC (KG/M**2)       1 PR
      IJ_EVAP = 6  ! EVAP (KG/M**2)       1 SF
C**** IJ_SSAT = 7  ! SOIL SATURATION RATIO 1 GD ! OBS
      IJ_BETA = 7  ! BETA (1)             1 GD
C**** IJ_SLP1 = 8  ! SLP (100 PA-1000) (USING T1) 4 DA ! OBS
C**** IJ_P4UV = 8  ! 4*P4 (100 PA)  (UV GRID)     4 DA ! OBS
      IJ_PRES = 8  ! PIJ (100 PA)  (NO PRINTOUT)  4 DA
      IJ_PHI1K = 9 ! PHI1000 (M**2/S**2) 4 DA
      IJ_PHI850 = 10  ! PHI850 (M**2/S**2-1500*GRAV) 4 DA
      IJ_PHI700 = 11  ! PHI700-3000*GRAV  4 DA
      IJ_PHI500 = 12  ! PHI500-5600*GRAV  4 DA
      IJ_PHI300 = 13  ! PHI300-9500*GRAV  4 DA
      IJ_PHI100 = 14  ! PHI100-16400*GRAV 4 DA
      IJ_PHI30 = 15   ! PHI30-24000*GRAV   4 DA
      IJ_T850 = 16 ! T850-273.16 (K-273.16)*GRAV) (NO PRT) 4 DA
      IJ_PMCCLD = 17  ! PCLDMC (1)  (COMPOSITE OVER ATMOSPHERE) 2 RD
      IJ_CLDTPPR = 18 ! P-CLOUD TOP   (100 PA)                  2 RD
      IJ_CLDCV = 19   ! PCLD (1)  (COMPOSITE OVER ATMOSPHERE)   2 RD
      IJ_PEV  = 20
       ! 16*P4*(SHA*T4+Z4)*V1*DSIG*DXV (100 W*M/S**2) (UV GRID) 4 DA
      IJ_TRNFP0 = 21  ! TRNFP0 (W/M**2)                         2 RS
      IJ_SRTR = 22    ! SRHDT+TRHDT (J/M**2)                    1 SF
      IJ_NETH = 23    ! SRHDT+TRHDT+SHDT+EVHDT+ENRGP (J/M**2)   1 SP
      IJ_SRNFP0 = 24  ! SRNFP0 (W/M**2)                         2 RD
      IJ_SRINCP0 = 25 ! SRINCP0 (W/M**2)                        2 RD
      IJ_SRNFG = 26   ! SRNFG (W/M**2)                          2 RD
      IJ_SRINCG = 27  ! SRINCG (W/M**2)                         2 RD
      IJ_TG1  = 28  ! TG1 (K-273.16)                            1 GD
      IJ_RSIT = 29  ! POICE+PLICE+(IF SNOW)PEARTH               4 DA
      IJ_TDSL = 30 ! DIURNAL DELTA TS (K) OVER SOIL (NO PRT) .5*9 MN
      IJ_DTDP = 31  ! DTHETA/DPHI (K S**2/M**2) IN TROPOSPHERE  4 DA
      IJ_RUNE = 32  ! RUN1 OVER EARTH  (KG/M**2)                1 PG
C**** IJ_TS1  = 33 ! TS (K-273.16) (W/ LAPSE RATE FROM TX1) 4 DA ! OBS
      IJ_RUNLI = 33 ! RUN1 OVER LAND ICE  (KG/M**2) (NO PRT)    1 PG
      IJ_WS   = 34  ! SURFACE WIND SPEED (M/S)                  3 SF
      IJ_TS   = 35  ! TS (K-273.16)                             3 SF
      IJ_US   = 36  ! US (M/S)                                  3 SF
      IJ_VS   = 37  ! VS (M/S)                                  3 SF
      IJ_SLP  = 38  ! PSL (100 PA-1000)  (USING TS)             4 DA
      IJ_UJET = 39  ! UJET (M/S)                                4 DA
      IJ_VJET = 40  ! VJET (M/S)                                4 DA
      IJ_PCLDL = 41  ! PCLD(LOW) (1)                            2 RD
      IJ_PCLDM = 42  ! PCLD(MID) (1)                            2 RD
      IJ_PCLDH = 43  ! PCLD(HIGH) (1)                           2 RD
      IJ_BTMPW = 44  ! BTEMPW-TF (K-273.16)                     2 RD
      IJ_SRREF = 45  ! PLAVIS*S0*COSZ (W/M**2)                  2 RD
      IJ_ODATA4 = 46  ! TGO2= ODATA(4)  (C)                  .5*9 MN
      IJ_TAUS = 47  ! TAUS  (MOM. SURF. DRAG) (kg/m**2) (NO PRT)  3 SF
      IJ_TAUUS = 48 ! TAUUS (MOM. SURF. DRAG) (kg/m**2) (NO PRT)  3 SF
      IJ_TAUVS = 49 ! TAUVS (MOM. SURF. DRAG) (kg/m**2) (NO PRT)  3 SF
      IJ_GWTR = 50  ! WATER1+WATER2+ICE1+ICE2 (EARTH POINTS ONLY) 1 GD
      IJ_QS   = 51  ! QS                                (NO PRT)  3 SF
      IJ_STRNGTS = 52 ! MAX(0,33-1.8*DAILY MEAN ON TS IN C)    .5*9 MN
      IJ_ARUNU = 53  ! 40.6+.72*(2TS(C)-(QSATS-QS)*LHA/SHA)       3 SF
      IJ_DTGDTS = 54 ! 18*(DEL(TG)/DEL(TS)-1),DEL= DIURN MX-MN .5*9 MN
      IJ_PUQ  = 55  ! 8*P*U*Q (VERT. INTEGRATED) (12.5 PA*M/S) 4 DA
      IJ_PVQ  = 56  ! 8*P*V*Q (VERT. INTEGRATED) (12.5 PA*M/S) 4 DA
      IJ_TGO  = 57  ! TGO= ODATA(1)  (C)                       1 GD
      IJ_MSI2 = 58  ! ACE2OI= ODATA(3)*POICE  (KG/M**2)        1 GD
      IJ_WLM  = 59  ! WIND SPEED IN TOP LAYER (M/S)            1 SD
      IJ_TGO2 = 60  ! TGO12= ODATA(5)  (C)                  .5*9 MN
      IJ_EVAPO = 61  ! EVAP*POCEAN  (KG/M**2)                  1 GD
      IJ_EVAPI = 62  ! EVAP*POICE  (KG/M**2)                   1 GD
      IJ_EVAPLI = 63  ! EVAP OVER LAND ICE  (KG/M**2)          1 GD
      IJ_EVAPE = 64  ! EVAP OVER EARTH  (KG/M**2)              1 GD
      IJ_F0OC = 65  ! F0DT*POCEAN, NET HEAT AT Z0  (J/M**2)    1 GD
      IJ_F0OI = 66  ! F0DT*POICE, NET HEAT AT Z0  (J/M**2)     1 GD
      IJ_F0LI = 67  ! F0DT, NET HEAT AT Z0 OVER LAND ICE  (J/M**2) 1 GD
      IJ_F0E  = 68  ! F0DT, NET HEAT AT Z0 OVER EARTH  (J/M**2) 1 GD
      IJ_F1LI = 69  ! F1DT OVER LAND ICE  (J/M**2)             1 PG
      IJ_SNWF = 70  ! SNOW FALL  (KG/M**2)                     1 PR
      IJ_TSLI = 71  ! SURF AIR TEMP OVER LAND ICE  (C)   NSURF*1 SF
      IJ_ERUN2 = 72  ! F2DT OVER LAND ICE  (J/M**2)            1 PG
      IJ_SHDTLI = 73  ! SHDT OVER LAND ICE  (J/M**2)           3 SF
      IJ_EVHDT = 74  ! EVHDT OVER LAND ICE  (J/M**2)           3 SF
      IJ_TRHDT = 75  ! TRHDT OVER LAND ICE  (J/M**2)           3 SF
      IJ_TMAX  = 76  ! MAX(COMPOSITE TS)                      12 SF
      IJ_TMIN  = 77  ! MIN(COMPOSITE TS)                      12 SF
      IJ_TMNMX  = 78  ! MIN(DIURNAL MAX OF COMPOSITE TS)      12 MN
      IJ_PEVAP = 79  ! POTENTIAL EVAPORATION (KG/M**2)         1 EA
      IJ_TMAXE = 80  ! MAX TS OVER EARTH FOR CURRENT DAY (K).5*9 MN
      IJ_WMSUM = 81  ! LIQUID WATER PATH (kg/M**2)             1 CL
      IJ_PSCLD = 82  ! SHALLOW CONVECTIVE CLOUD COVER  (1)     1 CL
      IJ_PDCLD = 83  ! DEEP CONVECTIVE CLOUD COVER     (1)     1 CL
      IJ_DCNVFRQ = 84  ! DEEP CONVECTIVE CLOUD FREQUENCY (1)   1 CL
      IJ_SCNVFRQ = 85  ! SHALLOW CONVECTIVE CLOUD FREQUENCY (1) 1 CL
      IJ_EMTMOM = 86 ! INCIDENT MTN EAST MOM. FLUX (MB-M/S**2)  1 SD
      IJ_SMTMOM = 87 ! INCIDENT MTN SOUTH MOM. FLUX (MB-M/S**2) 1 SD
      IJ_FPEU = 88  ! EAST-WEST POT. ENTHALPY FLUX (W)   /3600.*1 DY
      IJ_FPEV = 89  ! NORTH-SOUTH POT. ENTHALPY FLUX (W) /3600.*1 DY
      IJ_FMU  = 90 ! EAST-WEST MASS FLUX (KG/S) 100./GRAV/3600.*1 DY
      IJ_FMV  = 91 ! NORTH-SOUTH MASS FLUX (KG/S) 100./GRAV/3600.*1 DY
      IJ_FQU = 92  ! EAST-WEST WATER VAPOR FLUX (KG/S)   /3600.*1 DY
      IJ_FQV = 93  ! NORTH-SOUTH WATER VAPOR FLUX (KG/S) /3600.*1 DY
      IJ_FGZU = 94 ! EAST-WEST GEOPOTENTIAL FLUX (W)    /3600.*1 DY
      IJ_FGZV = 95 ! NORTH-SOUTH GEOPOTENTIAL FLUX (W)  /3600.*1 DY
      IJ_ERVR=96 ! Energy Outflow by Rivers (10**10 W) E-10/DTS*1 RV
      IJ_MRVR=97 ! Mass Outflow by Rivers (10**5 kg/s)  E-5/DTS*1 RV
C**** Definiton of IA_IJ (IDACC index) and NAME_IJ
      NAME_IJ(IJ_RSOI)   = "RSOI" ; IA_IJ(IJ_RSOI) = 1
      NAME_IJ(IJ_RSNW)   = "RSNW" ; IA_IJ(IJ_RSNW) = 4
      NAME_IJ(IJ_SNOW)   = "SNOW" ; IA_IJ(IJ_SNOW) = 4
      NAME_IJ(IJ_SHDT)   = "SHDT" ; IA_IJ(IJ_SHDT) = 1
      NAME_IJ(IJ_PREC)   = "PREC" ; IA_IJ(IJ_PREC) = 1
      NAME_IJ(IJ_EVAP)   = "EVAP" ; IA_IJ(IJ_EVAP) = 1
C**** NAME_IJ(IJ_SSAT)   = "SSAT" ; IA_IJ(IJ_SSAT) = 1
      NAME_IJ(IJ_BETA)   = "BETA" ; IA_IJ(IJ_BETA) = 1
C**** NAME_IJ(IJ_SLP1)   = "SLP1" ; IA_IJ(IJ_SLP1) = 4
C**** NAME_IJ(IJ_P4UV)   = "P4UV" ; IA_IJ(IJ_P4UV) = 4
      NAME_IJ(IJ_PRES)   = "PRES" ; IA_IJ(IJ_PRES) = 4
      NAME_IJ(IJ_PHI1K)  = "PHI1K"  ; IA_IJ(IJ_PHI1K ) = 4
      NAME_IJ(IJ_PHI850) = "PHI850" ; IA_IJ(IJ_PHI850) = 4
      NAME_IJ(IJ_PHI700) = "PHI700" ; IA_IJ(IJ_PHI700) = 4
      NAME_IJ(IJ_PHI500) = "PHI500" ; IA_IJ(IJ_PHI500) = 4
      NAME_IJ(IJ_PHI300) = "PHI300" ; IA_IJ(IJ_PHI300) = 4
      NAME_IJ(IJ_PHI100) = "PHI100" ; IA_IJ(IJ_PHI100) = 4
      NAME_IJ(IJ_PHI30)  = "PHI30"  ; IA_IJ(IJ_PHI30 ) = 4
      NAME_IJ(IJ_T850)   = "T850"   ; IA_IJ(IJ_T850) = 4
      NAME_IJ(IJ_PMCCLD) = "PMCCLD" ; IA_IJ(IJ_PMCCLD) = 2
      NAME_IJ(IJ_CLDTPPR)= "CLDTPPR"; IA_IJ(IJ_CLDTPPR) = 2
      NAME_IJ(IJ_CLDCV)  = "CLDCV"  ; IA_IJ(IJ_CLDCV) = 2
      NAME_IJ(IJ_PEV)    = "PEV"    ; IA_IJ(IJ_PEV) = 4
      NAME_IJ(IJ_TRNFP0) = "TRNFP0" ; IA_IJ(IJ_TRNFP0) = 2
      NAME_IJ(IJ_SRTR)   = "SRTR"   ; IA_IJ(IJ_SRTR) = 1
      NAME_IJ(IJ_NETH)   = "NETH"   ; IA_IJ(IJ_NETH) = 1
      NAME_IJ(IJ_SRNFP0) = "SRNFP0" ; IA_IJ(IJ_SRNFP0) = 2
      NAME_IJ(IJ_SRINCP0)= "SRINCP0"; IA_IJ(IJ_SRINCP0) = 2
      NAME_IJ(IJ_SRNFG)  = "SRNFG"  ; IA_IJ(IJ_SRNFG) = 2
      NAME_IJ(IJ_SRINCG) = "SRINCG" ; IA_IJ(IJ_SRINCG) = 2
      NAME_IJ(IJ_TG1)    = "TG1"    ; IA_IJ(IJ_TG1) = 1
      NAME_IJ(IJ_RSIT)   = "RSIT"   ; IA_IJ(IJ_RSIT) = 4
      NAME_IJ(IJ_TDSL)   = "TDSL"   ; IA_IJ(IJ_TDSL) = 9
      NAME_IJ(IJ_DTDP)   = "DTDP"   ; IA_IJ(IJ_DTDP) = 4
      NAME_IJ(IJ_RUNE)   = "RUNE"   ; IA_IJ(IJ_RUNE) = 1
C**** NAME_IJ(IJ_TS1)    = "TS1"    ; IA_IJ(IJ_TS1) = 4
      NAME_IJ(IJ_RUNLI)  = "RUNLI"  ; IA_IJ(IJ_RUNLI) = 1
      NAME_IJ(IJ_WS)     = "WS"     ; IA_IJ(IJ_WS) = 3
      NAME_IJ(IJ_TS)     = "TS"     ; IA_IJ(IJ_TS) = 3
      NAME_IJ(IJ_US)     = "US"     ; IA_IJ(IJ_US) = 3
      NAME_IJ(IJ_VS)     = "VS"     ; IA_IJ(IJ_VS) = 3
      NAME_IJ(IJ_SLP)    = "SLP"    ; IA_IJ(IJ_SLP) = 4
      NAME_IJ(IJ_UJET)   = "UJET"   ; IA_IJ(IJ_UJET) = 4
      NAME_IJ(IJ_VJET)   = "VJET"   ; IA_IJ(IJ_VJET) = 4
      NAME_IJ(IJ_PCLDL)  = "PCLDL"  ; IA_IJ(IJ_PCLDL) = 2
      NAME_IJ(IJ_PCLDM)  = "PCLDM"  ; IA_IJ(IJ_PCLDM) = 2
      NAME_IJ(IJ_PCLDH)  = "PCLDH"  ; IA_IJ(IJ_PCLDH) = 2
      NAME_IJ(IJ_BTMPW)  = "BTMPW"  ; IA_IJ(IJ_BTMPW) = 2
      NAME_IJ(IJ_SRREF)  = "SRREF"  ; IA_IJ(IJ_SRREF) = 2
      NAME_IJ(IJ_ODATA4) = "ODATA4" ; IA_IJ(IJ_ODATA4) = 9
      NAME_IJ(IJ_TAUS)   = "TAUS"   ; IA_IJ(IJ_TAUS) = 3
      NAME_IJ(IJ_TAUUS)  = "TAUUS"  ; IA_IJ(IJ_TAUUS) = 3
      NAME_IJ(IJ_TAUVS)  = "TAUVS"  ; IA_IJ(IJ_TAUVS) = 3
      NAME_IJ(IJ_GWTR)   = "GWTR"   ; IA_IJ(IJ_GWTR) = 1
      NAME_IJ(IJ_QS)     = "QS"     ; IA_IJ(IJ_QS) = 3
      NAME_IJ(IJ_STRNGTS)= "STRNGTS"; IA_IJ(IJ_STRNGTS) = 9
      NAME_IJ(IJ_ARUNU)  = "ARUNU"  ; IA_IJ(IJ_ARUNU) = 3
      NAME_IJ(IJ_DTGDTS) = "DTGDTS" ; IA_IJ(IJ_DTGDTS) = 9
      NAME_IJ(IJ_PUQ)    = "PUQ"    ; IA_IJ(IJ_PUQ) = 4
      NAME_IJ(IJ_PVQ)    = "PVQ"    ; IA_IJ(IJ_PVQ) = 4
      NAME_IJ(IJ_TGO)    = "TGO"    ; IA_IJ(IJ_TGO) = 1
      NAME_IJ(IJ_MSI2)   = "MSI2"   ; IA_IJ(IJ_MSI2) = 1
      NAME_IJ(IJ_WLM)    = "WLM"    ; IA_IJ(IJ_WLM) = 1
      NAME_IJ(IJ_TGO2)   = "TGO2"   ; IA_IJ(IJ_TGO2) = 9
      NAME_IJ(IJ_EVAPO)  = "EVAPO"  ; IA_IJ(IJ_EVAPO) = 1
      NAME_IJ(IJ_EVAPI)  = "EVAPI"  ; IA_IJ(IJ_EVAPI) = 1
      NAME_IJ(IJ_EVAPLI) = "EVAPLI" ; IA_IJ(IJ_EVAPLI) = 1
      NAME_IJ(IJ_EVAPE)  = "EVAPE"  ; IA_IJ(IJ_EVAPE) = 1
      NAME_IJ(IJ_F0OC)   = "F0OC"   ; IA_IJ(IJ_F0OC) = 1
      NAME_IJ(IJ_F0OI)   = "F0OI"   ; IA_IJ(IJ_F0OI) = 1
      NAME_IJ(IJ_F0LI)   = "F0LI"   ; IA_IJ(IJ_F0LI) = 1
      NAME_IJ(IJ_F0E)    = "F0E"    ; IA_IJ(IJ_F0E) = 1
      NAME_IJ(IJ_F1LI)   = "F1LI"   ; IA_IJ(IJ_F1LI) = 1
      NAME_IJ(IJ_SNWF)   = "SNWF"   ; IA_IJ(IJ_SNWF) = 1
      NAME_IJ(IJ_TSLI)   = "TSLI"   ; IA_IJ(IJ_TSLI) = 1
      NAME_IJ(IJ_ERUN2)  = "ERUN2"  ; IA_IJ(IJ_ERUN2) = 1
      NAME_IJ(IJ_SHDTLI) = "SHDTLI" ; IA_IJ(IJ_SHDTLI) = 3
      NAME_IJ(IJ_EVHDT)  = "EVHDT"  ; IA_IJ(IJ_EVHDT) = 3
      NAME_IJ(IJ_TRHDT)  = "TRHDT"  ; IA_IJ(IJ_TRHDT) = 3
      NAME_IJ(IJ_TMAX)   = "TMAX"   ; IA_IJ(IJ_TMAX) = 12
      NAME_IJ(IJ_TMIN)   = "TMIN"   ; IA_IJ(IJ_TMIN) = 12
      NAME_IJ(IJ_TMNMX)  = "TMNMX"  ; IA_IJ(IJ_TMNMX) = 12
      NAME_IJ(IJ_PEVAP)  = "PEVAP"  ; IA_IJ(IJ_PEVAP) = 1
      NAME_IJ(IJ_TMAXE)  = "TMAXE"  ; IA_IJ(IJ_TMAXE) = 9
      NAME_IJ(IJ_WMSUM)  = "WMSUM"  ; IA_IJ(IJ_WMSUM) = 1
      NAME_IJ(IJ_PSCLD)  = "PSCLD"  ; IA_IJ(IJ_PSCLD) = 1
      NAME_IJ(IJ_PDCLD)  = "PDCLD"  ; IA_IJ(IJ_PDCLD) = 1
      NAME_IJ(IJ_DCNVFRQ)= "DCNVFRQ"; IA_IJ(IJ_DCNVFRQ) = 1
      NAME_IJ(IJ_SCNVFRQ)= "SCNVFRQ"; IA_IJ(IJ_SCNVFRQ) = 1
      NAME_IJ(IJ_EMTMOM) = "EMTMOM" ; IA_IJ(IJ_EMTMOM) = 1
      NAME_IJ(IJ_SMTMOM) = "SMTMOM" ; IA_IJ(IJ_SMTMOM) = 1
      NAME_IJ(IJ_FPEU)   = "FPEU"   ; IA_IJ(IJ_FPEU) = 1
      NAME_IJ(IJ_FPEV)   = "FPEV"   ; IA_IJ(IJ_FPEV) = 1
      NAME_IJ(IJ_FMU)    = "FMU"    ; IA_IJ(IJ_FMU) = 1
      NAME_IJ(IJ_FMV)    = "FMV"    ; IA_IJ(IJ_FMV) = 1
      NAME_IJ(IJ_FQU)    = "FQU"    ; IA_IJ(IJ_FQU) = 1
      NAME_IJ(IJ_FQV)    = "FQV"    ; IA_IJ(IJ_FQV) = 1
      NAME_IJ(IJ_FGZU)   = "FGZU"   ; IA_IJ(IJ_FGZU) = 1
      NAME_IJ(IJ_FGZV)   = "FGZV"   ; IA_IJ(IJ_FGZV) = 1
      NAME_IJ(IJ_ERVR)   = "ERVR"   ; IA_IJ(IJ_ERVR) = 1
      NAME_IJ(IJ_MRVR)   = "MRVR"   ; IA_IJ(IJ_MRVR) = 1

C**** AJ diagnostic names: (SUM OVER LONGITUDE AND TIME)
C**** NAME     NO.    DESCRIPTION   (SCALE)*IDACC  LOCATION 
C**** ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      J_SRINCP0= 1 ! SRINCP0 (W/M**2)                              2 RD
      J_SRNFP0=  2 ! SRNFP0 (W/M**2)                               2 RD
      J_SRNFP1=  3 ! SRNFP1 (W/M**2)                               2 RD
      J_SRABS =  4 ! SRABSATM=AJ(SRNFP0)-AJ(SRNFG) (W/M**2)        2 D1
      J_SRINCG=  5 ! SRINCG (W/M**2)                               2 RD
      J_SRNFG =  6 ! SRNFG (W/M**2)                                2 RD
      J_TRNFP0=  7 ! TRNFP0=AJ(ALBVIS)+A2BYA1*AJ(TRHDT)/DTS(W/M**2)2 D1
      J_TRNFP1=  8 ! TRNFP1=AJ(ALBNIR)+A2BYA1*AJ(TRHDT)/DTS(W/M**2)2 D1
      J_TRHDT =  9 ! TRHDT (J/M**2)                                1 SF
      J_RNFP0 = 10 ! RNFP0=AJ(SRNFP0)+AJ(TRNFP0) (W/M**2)          2 D1
      J_RNFP1 = 11 ! RNFP1=AJ(SRNFP1)+AJ(TRNFP1) (W/M**2)          2 D1
      J_RHDT  = 12 ! RHDT=A1BYA2*AJ(SRNFG)*DTS+AJ(TRHDT)(J/M^2)    1 D1
      J_SHDT  = 13 ! SHEATDT (J/M**2)                              1 SF
      J_EVHDT = 14 ! EVHDT (J/M**2)                                1 SF
      J_F2DT  = 15 ! F2DT (J/M**2)                                 1 GD
      J_HZ1   = 16 ! HEATZ1=AJ(EDIFS)+AJ(F1DT)                     1 D1
      J_TG2   = 17 ! TG2 (K-273.16)                                1 GD
      J_TG1   = 18 ! TG1 (K-273.16)                                1 GD
      J_EVAP  = 19 ! EVAP (KG/M**2)                                1 GD
      J_PRCP  = 20 ! PRCP=AJ(PRCPSS)+AJ(PRCPMC) (100 PA)           1 D1
      J_TX    = 21 ! TX (K-273.16)  (INTEGRAL OVER ATMOSPHERE OF)  4 DA
      J_TX1   = 22 ! TX1 (K-273.16)                                4 DA
      J_TSRF  = 23 ! TS (K-273.16)                                 3 SF
      J_DTSGST= 24 ! DTH/DPHI  (STRATOSPHERE)                      4 DA
      J_DTDGTR= 25 ! DTH/DPHI  (TROPOSPHERE)                       4 DA
      J_RICST = 26 ! .0625*DTH*DLNP/(DU*DU+DV*DV)  (STRATOSPHERE)  4 DA
      J_RICTR = 27 ! .0625*DTH*DLNP/(DU*DU+DV*DV)  (TROPOSPHERE)   4 DA
      J_ROSST = 28 ! 4*UMAX/(DX*SINJ)  (STRATOSPHERE)              4 DA
      J_ROSTR = 29 ! 4*UMAX/(DX*SINJ)  (TROPOSPHERE)               4 DA
      J_RSI   = 30 ! POICE (1)                                     1 GD
      J_RSNOW = 31 ! PSNOW (1)                                     4 DA
      J_SWCOR = 32 ! SW CORRECTION   (obsolete)                    2 RD
      J_OHT   = 33 ! OCEAN TRANSPORT                               1 GD
      J_OMLT  = 34 ! OCEAN TEMPERATURE AT MAX. MIXED LAYER DEPTH   1 GD
      J_DTDJS = 35 ! T(J+1)-T(J-1)  (SUM OVER STRATOSPHERE OF)     4 DA
      J_DTDJT = 36 ! T(J+1)-T(J-1)  (SUM OVER TROPOSPHERE OF)      4 DA
      J_LSTR  = 37 ! SQRT(DTH/DLNP)/SINJ  (STRATOSPHERE)           4 DA
      J_LTRO  = 38 ! SQRT(DTH/DLNP)/SINJ  (TROPOSPHERE)            4 DA
      J_EPRCP = 39 ! ENERGP (J/M**2)                               1 PR
      J_ERUN1 = 40 ! ERUN1 (J/M**2)                                1 GP
      J_EDIFS = 41 ! EDIFS (J/M**2)                                1 GP
      J_F1DT  = 42 ! F1DT (J/M**2)                                 1 GD
      J_ERUN2 = 43 ! ERUN2 (J/M**2)                                1 GP
      J_HZ0   = 44 ! HEATZ0=RHDT+SHDT+EVHDT+EPRCP-EDIFS (J/M**2)   1 D1
      J_DIFS  = 45 ! DIFS (KG/M**2)                                1 GP
      J_IMELT = 46 ! AIFO ; BRUN2 ; CRUN2+CIFI                     1 GP
      J_RUN2  = 47 ! RUN2 (KG/M**2)                                1 GP
      J_DWTR2 = 48 ! DWTR2=AJ(DIFS)-AJ(RUN2) (KG/M**2)             1 D1
      J_WTR1  = 49 ! WTR1 (KG/M**2)                                1 GD
      J_ACE1  = 50 ! ACE1 (KG/M**2)                                1 GD
      J_WTR2  = 51 ! WTR2 (KG/M**2)                                1 GD
      J_ACE2  = 52 ! ACE2 (KG/M**2)                                1 GD
      J_SNOW  = 53 ! SNOW (KG/M**2)                                1 GD
      J_RUN1  = 54 ! RUN1 (KG/M**2)                                1 GP
      J_BRTEMP= 55 ! BTEMPW-TF                                     2 RD
      J_HZ2   = 56 ! HEATZ2=AJ(F2DT)+AJ(ERUN2) (J/M**2)            1 D1
      J_PCLDSS= 57 ! PCLDSS (1)  (COMPOSITE OVER ATMOSPHERE)       2 RD
      J_PCLDMC= 58 ! PCLDMC (1)  (COMPOSITE OVER ATMOSPHERE)       2 RD
      J_PCLD  = 59 ! PCLD (1)  (COMPOSITE OVER ATMOSPHERE)         2 RD
      J_CTOPP = 60 ! CLDTOPMC=AJ(CDLDEP)/AJ(PCLDMC) (100 PA)       0 D1
      J_PRCPSS= 61 ! PRCPSS (100 PA)                               1 CN
      J_PRCPMC= 62 ! PRCPMC (100 PA)                               1 CN
      J_QP    = 63 ! Q*P (100 PA)  (INTEGRAL OVER ATMOSPHERE OF)   4 DA
      J_GAM   = 64 ! GAM  (K/M)  (*SIG(TROPOSPHERE)/GRAV)          4 DA
      J_GAMM  = 65 ! GAMM  (K-S**2/M**2)  (SIG(TROPOSPHERE)/GAMD)  4 DA
      J_GAMC  = 66 ! GAMC  (K/M)                                   4 DA
      J_TRINCG= 67 ! TRINCG (W/M**2)                               2 RD
      J_FTHERM= 68 ! ENERGY DIFFUSION INTO THERMOCLINE (W/M**2) .5*9 MN
      J_HSURF = 70 ! TRNFP0-TRNFG (W/M**2)                         2 RD
      J_HATM  = 71 ! TRNFP1-TRNFG (W/M**2)                         2 RD
C**** Note: next eight diagnostics must remain in order
      J_PLAVIS= 72 ! PLAVIS*S0*COSZ (W/M**2)                       2 RD
      J_PLANIR= 73 ! PLANIR*S0*COSZ (W/M**2)                       2 RD
      J_ALBVIS= 74 ! ALBVIS*S0*COSZ (W/M**2)                       2 RD
      J_ALBNIR= 75 ! ALBNIR*S0*COSZ (W/M**2)                       2 RD
      J_SRRVIS= 76 ! SRRVIS*S0*COSZ (W/M**2)                       2 RD
      J_SRRNIR= 77 ! SRRNIR*S0*COSZ (W/M**2)                       2 RD
      J_SRAVIS= 78 ! SRAVIS*S0*COSZ (W/M**2)                       2 RD
      J_SRANIR= 79 ! SRANIR*S0*COSZ (W/M**2)                       2 RD
      J_CDLDEP= 80 ! PBOTMC-PTOPMC (100 PA)                        2 RD

C**** Definiton of IA_J (IDACC index) and NAME_J
      NAME_J(J_SRINCP0)= "SRINCP0" ; IA_J(J_SRINCP0)= 2
      NAME_J(J_SRNFP0) = "SRNFP0"  ; IA_J(J_SRNFP0) = 2
      NAME_J(J_SRNFP1) = "SRNFP1"  ; IA_J(J_SRNFP1) = 2
      NAME_J(J_SRABS)  = "SRABS"   ; IA_J(J_SRABS)  = 2
      NAME_J(J_SRINCG) = "SRINCG"  ; IA_J(J_SRINCG) = 2
      NAME_J(J_SRNFG)  = "SRNFG"   ; IA_J(J_SRNFG)  = 2
      NAME_J(J_TRNFP0) = "TRNFP0"  ; IA_J(J_TRNFP0) = 2
      NAME_J(J_TRNFP1) = "TRNFP1"  ; IA_J(J_TRNFP1) = 2
      NAME_J(J_TRHDT)  = "TRHDT"   ; IA_J(J_TRHDT) = 1
      NAME_J(J_RNFP0)  = "RNFP0"   ; IA_J(J_RNFP0) = 2
      NAME_J(J_RNFP1)  = "RNFP1"   ; IA_J(J_RNFP1) = 2
      NAME_J(J_RHDT)   = "RHDT"    ; IA_J(J_RHDT)  = 1
      NAME_J(J_SHDT)   = "SHDT"    ; IA_J(J_SHDT)  = 1
      NAME_J(J_EVHDT)  = "EVHDT"   ; IA_J(J_EVHDT) = 1
      NAME_J(J_F2DT)   = "F2DT"    ; IA_J(J_F2DT)  = 1
      NAME_J(J_HZ1)    = "HZ1"     ; IA_J(J_HZ1) = 1
      NAME_J(J_TG2)    = "TG2"     ; IA_J(J_TG2) = 1
      NAME_J(J_TG1)    = "TG1"     ; IA_J(J_TG1) = 1
      NAME_J(J_EVAP)   = "EVAP"    ; IA_J(J_EVAP) = 1
      NAME_J(J_PRCP)   = "PRCP"    ; IA_J(J_PRCP) = 1
      NAME_J(J_TX)     = "TX"      ; IA_J(J_TX)   = 4
      NAME_J(J_TX1)    = "TX1"     ; IA_J(J_TX1) = 4
      NAME_J(J_TSRF)   = "TSRF"    ; IA_J(J_TSRF) = 3
      NAME_J(J_DTSGST) = "DTSGST"  ; IA_J(J_DTSGST) = 4
      NAME_J(J_DTDGTR) = "DTDGTR"  ; IA_J(J_DTDGTR) = 4
      NAME_J(J_RICST)  = "RICST"   ; IA_J(J_RICST) = 4
      NAME_J(J_RICTR)  = "RICTR"   ; IA_J(J_RICTR) = 4
      NAME_J(J_ROSST)  = "ROSST"   ; IA_J(J_ROSST) = 4
      NAME_J(J_ROSTR)  = "ROSTR"   ; IA_J(J_ROSTR) = 4
      NAME_J(J_RSI)    = "RSI"     ; IA_J(J_RSI) = 1
      NAME_J(J_RSNOW)  = "RSNOW"   ; IA_J(J_RSNOW) = 4
      NAME_J(J_SWCOR)  = "SWCOR"   ; IA_J(J_SWCOR) = 2
      NAME_J(J_OHT)    = "OHT"     ; IA_J(J_OHT) = 1
      NAME_J(J_OMLT)   = "OMLT"    ; IA_J(J_OMLT) = 1
      NAME_J(J_DTDJS)  = "DTDJS"   ; IA_J(J_DTDJS) = 4
      NAME_J(J_DTDJT)  = "DTDJT"   ; IA_J(J_DTDJT) = 4
      NAME_J(J_LSTR)   = "LSTR"    ; IA_J(J_LSTR) = 4
      NAME_J(J_LTRO)   = "LTRO"    ; IA_J(J_LTRO) = 4
      NAME_J(J_EPRCP)  = "EPRCP"   ; IA_J(J_EPRCP) = 1
      NAME_J(J_ERUN1)  = "ERUN1"   ; IA_J(J_ERUN1) = 1
      NAME_J(J_EDIFS)  = "EDIFS"   ; IA_J(J_EDIFS) = 1
      NAME_J(J_F1DT)   = "F1DT"    ; IA_J(J_F1DT)  = 1
      NAME_J(J_ERUN2)  = "ERUN2"   ; IA_J(J_ERUN2) = 1
      NAME_J(J_HZ0)    = "HZ0"     ; IA_J(J_HZ0) = 1
      NAME_J(J_DIFS)   = "DIFS"    ; IA_J(J_DIFS) = 1
      NAME_J(J_IMELT)  = "IMELT"   ; IA_J(J_IMELT) = 1
      NAME_J(J_RUN2)   = "RUN2"    ; IA_J(J_RUN2) = 1
      NAME_J(J_DWTR2)  = "DWTR2"   ; IA_J(J_DWTR2) = 1
      NAME_J(J_WTR1)   = "WTR1"    ; IA_J(J_WTR1) = 1
      NAME_J(J_ACE1)   = "ACE1"    ; IA_J(J_ACE1) = 1
      NAME_J(J_WTR2)   = "WTR2"    ; IA_J(J_WTR2) = 1
      NAME_J(J_ACE2)   = "ACE2"    ; IA_J(J_ACE2) = 1
      NAME_J(J_SNOW)   = "SNOW"    ; IA_J(J_SNOW) = 1
      NAME_J(J_RUN1)   = "RUN1"    ; IA_J(J_RUN1) = 1
      NAME_J(J_BRTEMP) = "BRTEMP"  ; IA_J(J_BRTEMP) = 2
      NAME_J(J_HZ2)    = "HZ2"     ; IA_J(J_HZ2) = 1
      NAME_J(J_PCLDSS) = "PCLDSS"  ; IA_J(J_PCLDSS) = 2
      NAME_J(J_PCLDMC) = "PCLDMC"  ; IA_J(J_PCLDMC) = 2
      NAME_J(J_PCLD)   = "PCLD"    ; IA_J(J_PCLD) = 2
      NAME_J(J_CTOPP)  = "CTOPP"   ; IA_J(J_CTOPP) = 0
      NAME_J(J_PRCPSS) = "PRCPSS"  ; IA_J(J_PRCPSS) = 1
      NAME_J(J_PRCPMC) = "PRCPMC"  ; IA_J(J_PRCPMC) = 1
      NAME_J(J_QP)     = "QP"      ; IA_J(J_QP) = 4
      NAME_J(J_GAM)    = "GAM"     ; IA_J(J_GAM) = 4
      NAME_J(J_GAMM)   = "GAMM"    ; IA_J(J_GAMM) = 4
      NAME_J(J_GAMC)   = "GAMC"    ; IA_J(J_GAMC) = 4
      NAME_J(J_TRINCG) = "TRINCG"  ; IA_J(J_TRINCG) = 2
      NAME_J(J_FTHERM) = "FTHERM"  ; IA_J(J_FTHERM) = 9
      NAME_J(J_HSURF)  = "HSURF"   ; IA_J(J_HSURF) = 2
      NAME_J(J_HATM)   = "HATM"    ; IA_J(J_HATM) = 2
      NAME_J(J_PLAVIS) = "PLAVIS"  ; IA_J(J_PLAVIS) = 2
      NAME_J(J_PLANIR) = "PLANIR"  ; IA_J(J_PLANIR) = 2
      NAME_J(J_ALBVIS) = "ALBVIS"  ; IA_J(J_ALBVIS) = 2
      NAME_J(J_ALBNIR) = "ALBNIR"  ; IA_J(J_ALBNIR) = 2
      NAME_J(J_SRRVIS) = "SRRVIS"  ; IA_J(J_SRRVIS) = 2
      NAME_J(J_SRRNIR) = "SRRNIR"  ; IA_J(J_SRRNIR) = 2
      NAME_J(J_SRAVIS) = "SRAVIS"  ; IA_J(J_SRAVIS) = 2
      NAME_J(J_SRANIR) = "SRANIR"  ; IA_J(J_SRANIR) = 2
      NAME_J(J_CDLDEP) = "CDLDEP"  ; IA_J(J_CDLDEP) = 2

C**** Ensure that diagnostics are reset at the beginning of the run
      IF (TAU.le.TAUI) THEN
         CALL reset_DIAG
         CALL daily_DIAG
      END IF

C**** Initialize certain arrays used by more than one print routine
      PMTOP=(PSF-PTOP)*SIGE(LM+1)+PTOP
      DO L=1,LM
        PLE(L)=SIGE(L+1)*(PSF-PTOP)+PTOP
        PLM(L)=SIG(L)*(PSF-PTOP)+PTOP
      END DO
      PLM(LM+1)=.75*PMTOP
      PLM(LM+2)=.35*PMTOP
      PLM(LM+3)=.1*PMTOP

      RETURN
      END SUBROUTINE init_DIAG

      SUBROUTINE reset_DIAG
!@sum  reset_DIAG resets/initiallises diagnostics
!@auth Original Development Team
!@ver  1.0
      USE E001M12_COM
      USE DAGCOM
      IMPLICIT NONE

      TAU0=TAU
      IDAY0=IDAY
      TOFDY0=TOFDAY
      JDATE0=JDATE
      JMNTH0=JMONTH
      JYEAR0=JYEAR

      IDACC(1:10)=0
      AJ=0    ; BJ=0   ; CJ=0     ; AREG=0
      APJ=0   ; AJL=0  ; ASJL=0   ; AIJ=0
      AIL=0   ; AIJG=0 ; ENERGY=0 ; CONSRV=0
      SPECA=0 ; ATPE=0 ; ADAILY=0 ; WAVE=0
      AJK=0   ; AIJK=0 ; AIJL=0   ; AJLSP=0

      AIJ(:,:,IJ_TMNMX)=1000.

      RETURN
      END SUBROUTINE reset_DIAG

      SUBROUTINE daily_DIAG
!@sum  daily_DIAG resets diagnostics at beginning of each day
!@auth Original Development Team
!@ver  1.0
      USE E001M12_COM
      USE DAGCOM
      IMPLICIT NONE

      INTEGER I,J
C**** INITIALIZE SOME ARRAYS AT THE BEGINNING OF SPECIFIED DAYS
      IF (JDAY.EQ.32) THEN
         DO J=1+JM/2,JM
            DO I=1,IM
               TSFREZ(I,J,1)=JDAY
            END DO
         END DO
         DO J=1,JM/2
            DO I=1,IM
               TSFREZ(I,J,2)=JDAY
            END DO
         END DO
      ELSEIF (JDAY.EQ.213) THEN
         DO J=1,JM/2
            DO I=1,IM
              TSFREZ(I,J,1)=JDAY
            END DO
         END DO
      END IF
C**** INITIALIZE SOME ARRAYS AT THE BEGINNING OF EACH DAY
      DO J=1,JM
         DO I=1,IM
            TDIURN(I,J,1)= 1000.
            TDIURN(I,J,2)=-1000.
            TDIURN(I,J,3)= 1000.
            TDIURN(I,J,4)=-1000.
            TDIURN(I,J,5)=    0.
            TDIURN(I,J,6)=-1000.
            TDIURN(I,J,7)=-1000.
            TDIURN(I,J,8)=-1000.
            IF (FEARTH(I,J).LE.0.) THEN
               TSFREZ(I,J,1)=365.
               TSFREZ(I,J,2)=365.
            END IF
         END DO
      END DO

       END SUBROUTINE daily_DIAG
