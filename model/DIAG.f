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
!@ver  1.0
      USE MODEL_COM, only : im,jm,lm
      USE RADNCB, only : LM_REQ
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION, DIMENSION(LM) :: PLE
      DOUBLE PRECISION, DIMENSION(LM) :: PLE_DN
      DOUBLE PRECISION, DIMENSION(LM+LM_REQ) :: PLM
!@var P1000K scaling to change reference pressure from 1mb to 1000mb
      DOUBLE PRECISION :: P1000K
!@var inci,incj print increments for i and j, so maps/tables fit on page
      integer, parameter :: inci=(im+35)/36,incj=(JM+23)/24, jmby2=jm/2
!@var linect = current line on page of print out
      integer linect

      END MODULE DAGPCOM

C****                                                             IDACC
C**** CONTENTS OF AJ(J,N)  (SUM OVER LONGITUDE AND TIME OF)
C****   1  SRINCP0 (W/M**2)                                        2 RD
C****   2  SRNFP0 (W/M**2)                                         2 RD
C****   3  SRNFP1 (W/M**2)                                         2 RD
C****   4  SRABSATM=AJ(2)-AJ(6) (W/M**2)                           2 D1
C****   5  SRINCG (W/M**2)                                         2 RD
C****   6  SRNFG (W/M**2)                                          2 RD
C****   7  TRNFP0=AJ(74)+A2BYA1*AJ(9)/DTSRC  (W/M**2)              2 D1
C****   8  TRNFP1=AJ(75)+A2BYA1*AJ(9)/DTSRC  (W/M**2)              2 D1
C****   9  TRHDT (J/M**2)                                          1 SF
C****  10  RNFP0=AJ(2)+AJ(7) (W/M**2)                              2 D1
C****  11  RNFP1=AJ(3)+AJ(8) (W/M**2)                              2 D1
C****  12  RHDT=A1BYA2*AJ(6)*DTSRC+AJ(9) (J/M**2)                  1 D1
C****  13  SHEATDT (J/M**2)                                        1 SF
C****  14  EVHDT (J/M**2)                                          1 SF
C****  15  F2DT (J/M**2)                                           1 GD
C****  16  HEATZ1=AJ(41)+AJ(42)                                    1 D1
C****  17  TG2 (K-TF)                                              1 GD
C****  18  TG1 (K-TF)                                              1 GD
C****  19  EVAP (KG/M**2)                                          1 GD
C****  20  PRCP=AJ(61)+AJ(62) (100 PA)                             1 D1
C****  21  TX (K-TF)  (INTEGRAL OVER ATMOSPHERE OF)                4 DA
C****  22  TX1 (K-TF)                                              4 DA
C****  23  TS (K-TF)                                               3 SF
C****  24  DTH/DPHI  (STRATOSPHERE)                                4 DA
C****  25  DTH/DPHI  (TROPOSPHERE)                                 4 DA
C****  26  .0625*DTH*DLNP/(DU*DU+DV*DV)  (STRATOSPHERE)            4 DA
C****  27  .0625*DTH*DLNP/(DU*DU+DV*DV)  (TROPOSPHERE)             4 DA
C****  28  4*UMAX/(DX*SINJ)  (STRATOSPHERE)                        4 DA
C****  29  4*UMAX/(DX*SINJ)  (TROPOSPHERE)                         4 DA
C****  30  POICE (1)                                               1 GD
C****  31  PSNOW (1)                                               1 GD
C****  32  SW CORRECTION                                           2 RD
C****  33  OCEAN TRANSPORT                                         1 GD
C****  34  OCEAN TEMPERATURE AT MAX. MIXED LAYER DEPTH             1 GD
C****  35  T(J+1)-T(J-1)  (SUM OVER STRATOSPHERE OF)               4 DA
C****  36  T(J+1)-T(J-1)  (SUM OVER TROPOSPHERE OF)                4 DA
C****  37  SQRT(DTH/DLNP)/SINJ  (STRATOSPHERE)                     4 DA
C****  38  SQRT(DTH/DLNP)/SINJ  (TROPOSPHERE)                      4 DA
C****  39  ENERGP (J/M**2)                                         1 CN
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
C****  69  PTYPE                                                   1 GD
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
C****  17  T-T0 (change of Th by dynamics)                         4 DA
C****  18  DU/DT BY STRAT DEFORM DRAG (M/S)                        1 SD
C****  19  PCLD*P (TOTAL)                                          1 CN
C****  20  DU/DT BY STRAT MTN DRAG (M S)                           1 SD
C****  21  DU/DT BY STRAT SHR DRAG (M S)                           1 SD
C****  22  DU/DT BY STRAT MC DRAG C=-10 (M/S)                      1 SD
C****  23  DU/DT BY STRAT MC DRAG C=+10 (M/S)                      1 SD
C****  24  DU/DT BY STRAT MC DRAG C=-40 (M/S)                      1 SD
C****  25  DU/DT BY STRAT MC DRAG C=+40 (M/S)                      1 SD
C****  26  DU/DT BY STRAT MC DRAG C=-20 (M/S)                      1 SD
C****  27  DU/DT BY STRAT MC DRAG C=+20 (M/S)                      1 SD
C****  28  PCLD*P (SS)                                             1 CN
C****  29  PCLD*P (MC)                                             1 CN
C****  30  FREE                                                    4 DA
C****  31  D  (M*M/S *       ) STRAT. DIFFUSION COEFF.             1 SD
C****  32  DU  (M/S)  STRATOSPHERIC DIFFUSION                      1 SD
C****  33  DT  (DEG-K/TIMESTEP)  HEATING BY STRATOSPHERIC DRAG     1 SD
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
C****  46  SD  (SUM OVER I, 150E TO IM)                            4 DA
C****  47  V-V*  =D((V-VI)*(T-TI)/DTHDP)/DP                        4 DA
C****  48  4*PU4I*PV4I/P4I (100 N/S**2)  (UV GRID)                 4 DA
C****  49  4*PUV4I (100 N/S**2)  (UV GRID)                         4 DA
C****  50  DT(MC)*P (100 PA*K)  CHANGE OF PHASE                    1 CN
C****  51  CLHE*DQ(MC BEFORE COND)*P (100 PA*K)                    1 CN
C****  52  DU/DT BY SDRAG (M S-2)                                  1 SD
C****  53  CHANGE OF LATENT HEAT BY MOSIT CONV.                    1 CN
C****  54  TURBULENT KINETIC ENERGY (W/M^2)                        1 AT
C****  55  CHANGE OF LATENT HEAT BY TURBULENCE                     1 AT
C****  56  TOTAL HEATING BY MOIST CONVECTION (Q1)                  1 CN
C****  57  TOTAL DRYING BY MOIST CONVECTION (Q2)                   1 CN
C****
C**** CONTENTS OF ASJL(J,L,N)  (SUM OVER LONGITUDE AND TIME OF)
C****   1  TX (C)                                                  4 DA
C****   2  PHI (M**2/S**2)                                         4 DA
C****   3  SRHR (W/M**2)                                           2 RD
C****   4  TRHR (W/M**2)                                           2 RD
C****
C**** CONTENTS OF AIJ(I,J,N)  (SUM OVER TIME OF)
C****   1  POICE (1)                                               1 GD
C****   2  PSNOW (1)                                               1 GD
C****   3  SNOW (KG/M**2)                                          1 GD
C****   4  SHDT (J/M**2)                                           1 SF
C****   5  PREC (KG/M**2)                                          1 CN
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
C****  16  T850-TF (K-TF)*GRAV)             (NO PRINTOUT)          4 DA
C****  17  PCLDMC (1)  (COMPOSITE OVER ATMOSPHERE)                 2 RD
C****  18  P-CLOUD TOP   (100 PA)                                  2 RD
C****  19  PCLD (1)  (COMPOSITE OVER ATMOSPHERE)                   2 RD
C****  20  16*P4*(SHA*T4+Z4)*V1*DSIG*DXV (100 W*M/S**2)  (UV GRID) 4 DA
C****  21  TRNFP0 (W/M**2)                                         2 RS
C****  22  SRHDT+TRHDT (J/M**2)                                    1 SF
C****  23  SRHDT+TRHDT+SHDT+EVHDT+ENRGP (J/M**2)                   1 SC
C****  24  SRNFP0 (W/M**2)                                         2 RD
C****  25  SRINCP0 (W/M**2)                                        2 RD
C****  26  SRNFG (W/M**2)                                          2 RD
C****  27  SRINCG (W/M**2)                                         2 RD
C****  28  TG1 (K-TF)                                              1 GD
C****  29  POICE+PLICE+(IF SNOW)PEARTH                             1 GD
C****  30  DIURNAL DELTA TS (K) OVER SOIL       (NO PRINTOUT)   .5*9 MN
C****  31  DTHETA/DPHI (K S**2/M**2) IN TROPOSPHERE                4 DA
C****  32  RUN1 OVER EARTH  (KG/M**2)                              1 PG
C****  33  TS (K-TF)  (USING LAPSE RATE FROM TX1)  (COMM'D OUT)    4 DA
C****  33  RUN1 OVER LAND ICE  (KG/M**2)            (NO PRINTOUT)  1 PG
C****  34  SURFACE WIND SPEED (M/S)                                3 SF
C****  35  TS (K-TF)                                               3 SF
C****  36  US (M/S)                                                3 SF
C****  37  VS (M/S)                                                3 SF
C****  38  PSL (100 PA-1000)  (USING TS)                           4 DA
C****  39  UJET (M/S)                                              4 DA
C****  40  VJET (M/S)                                              4 DA
C****  41  PCLD(LOW) (1)                                           2 RD
C****  42  PCLD(MID) (1)                                           2 RD
C****  43  PCLD(HIGH) (1)                                          2 RD
C****  44  BTEMPW-TF (K-TF)                                        2 RD
C****  45  PLAVIS*S0*COSZ (W/M**2)                                 2 RD
C****  46  TGO2=TOCEAN(2) (C)                                   .5*9 MN
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
C****  57  TGO=TOCEAN(1)  (C)                                      1 GD
C****  58  ACE2OI=MSI2*POICE  (KG/M**2)                            1 GD
C****  59  WIND SPEED IN TOP LAYER (M/S)                           1 SD
C****  60  TGO12=TOCEAN(3)  (C)                                 .5*9 MN
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
C****  71  SURF AIR TEMP OVER LAND ICE  (C)                 NISURF*1 SF
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
C****  98  DU/DT BY SDRAG (M S-2)                                  1 SD
C****  99  LAST DAY OF ICE-FREE LAKE (DAYS)                       12 DA
C**** 100  LAST DAY OF ICED-UP LAKE (DAYS)                        12 DA
C**** 101  LAYER 1 REL SAT OF BARE AMD VEG. SOIL (%)               1 EA
C**** 102  LAYER 2 REL SATURATION OF BARE SOIL (%)                 1 EA
C**** 103  LAYER 3 REL SATURATION OF BARE SOIL (%)                 1 EA
C**** 104  LAYER 4 REL SATURATION OF BARE SOIL (%)                 1 EA
C**** 105  BETA (BARE SOIL) (%)                                    1 EA
C**** 106  BETA (PENMAN) (%)                                       1 EA
C**** 107  CANOPY  REL SATURATION (%)                              1 EA
C**** 108  LAYER 1 REL SATURATION OF VEG. SOIL (%)                 1 EA
C**** 109  LAYER 2 REL SATURATION OF VEG. SOIL (%)                 1 EA
C**** 110  LAYER 3 REL SATURATION OF VEG. SOIL (%)                 1 EA
C**** 111  BETA (BARE SOIL & VEGETATION) (%)                       1 EA
C**** 112  CONDUCTANCE OF ATMOSPHERE (.01M/S)                      1 EA
C**** 113  CONDUCTANCE OF CANOPY (.01M/S)                          1 EA
C**** 114  PENMAN POTENTIAL EVAPORATION (KG/M**2)                  1 EA
C**** 115  TEMP OF LAYER 1 BARE SOIL AND SNOW (C)                  1 EA
C**** 116  TEMP OF SOIL LAYER 2 - BARE SOIL (C)                    1 EA
C**** 117  TEMP OF SOIL LAYER 3 - BARE SOIL (C)                    1 EA
C**** 118  BARE SOIL EVAPORATION (KG/M**2)                         1 EA
C**** 119  DRY CANOPY EVAPORATION (KG/M**2)                        1 EA
C**** 120  WET CANOPY EVAPORATION (KG/M**2)                        1 EA
C**** 121  TEMP OF CANOPY AND SNOW (C)                             1 EA
C**** 123  TEMP OF SOIL LAYER 2 1 VEGETATED SOIL (C)               1 EA
C**** 123  TEMP OF SOIL LAYER 2 - VEGETATED SOIL (C)               1 EA
C**** 124  TEMP OF SOIL LAYER 3 - VEGETATED SOIL (C)               1 EA
C**** 125  AVERAGE WATER TABLE (M)                                 1 EA
C**** 126  BETAV, OVER VEGETATION (PERCENT)                        1 EA
C**** 127  BETAT, TRANSPIRATION (PERCENT)                          1 EA
C**** 128  SNOW OVER BARE SOIL                                     1 EA
C**** 129  SNOW OVER VEGETATED SOIL                                1 EA
C****
C**** CONTENTS OF AIL(I,L,N)  (SUM OVER TIME OF)
C**** WE ARE NOT TAKING INTO ACCOUNT THE VARIATION OF MASS
C****   1  U (M/S) (SUM FOR J=J5SUV,J5NUV)            (PU GRID)    4 DA
C****   2  V (M/S) (SUM FOR J=J5SUV,J5NUV)            (PU GRID)    4 DA
C****   3  SD (100 N/S) (SUM FOR J=J5S,J5N)                        4 DA
C****   4  TX (K-TF) (SUM FOR J=J5S,J5N)                           4 DA
C****   5  RH (1) (SUM FOR J=J5S,J5N)                              4 DA
C****   6  DTX(MC)*P*DA (100 K*N) (SUM FOR J=J5S,J5N)              1 CN
C****   7  (SRHR+TRHR)*DA (W) (SUM FOR J=J5S,J5N)                  2 RD
C****   9  SD (100 N/S) (AT LAT 50 N)                              4 DA
C****  10  TX-TF  (AT LAT 50 N)                                    4 DA
C****  11  SR+TR  (AT LAT 50 N)  (COMMENTED OUT)                   2 RD
C****  12  2*U  (AT LAT 50 N)                                      4 DA
C****  13  SD  (AT LAT 70 N)                                       4 DA
C****  14  TX-TF  (AT LAT 70 N)                                    4 DA
C****  15  SR+TR  (AT LAT 70 N)  (COMMENTED OUT)                   2 RD
C****  16  2*U  (AT LAT 70 N)                                      4 DA
C****
C****
C**** CONTENTS OF IDACC(N), NUMBER OF ACCUMULATION TIMES OF
C****   1  SOURCE TERMS  (dt: DTSRC)
C****   2  RADIATION SOURCE TERMS  (dt: NRAD*DTsrc)
C****   3  SURFACE INTERACTION SOURCE TERMS  (dt: NDASF*DTsrc+DTsurf)
C****   4  QUANTITIES IN DIAGA  (dt: NDAA*DTsrc+2*DTdyn)
C****   5  ENERGY NUMBERS IN DIAG4  (DEYERMINED BY NDA4)
C****   6  KINETIC ENERGY IN DIAG5 FROM DYN'CS (dt: NDA5K*DTsrc+2*DTdyn)
C****   7  ENERGY IN DIAG5 FROM DYNAMICS  (dt: NDA5D*DTsrc)
C****   8  ENERGY IN DIAG5 FROM SOURCES  (DETERMINED BY NDA5S)
C****   9  WAVE ENERGY IN DIAG7  (dt: 12 HOURS)
C****  10  ENERGY IN DIAG5 FROM FILTER  (DT: NFILTR*DTsrc)
C****  11  NOT USED
C****  12  ALWAYS =1 (UNLESS SEVERAL RESTART FILES WERE ACCUMULATED)
C****
C**** CONTENTS OF AUXILIARY ARRAYS (TSFREZ(I,J,1-4),TDIURN(I,J,N))
C****   1  FIRST DAY OF GROWING SEASON (JULIAN DAY)
C****   2  LAST DAY OF GROWING SEASON (JULIAN DAY)
C****   3  LAST DAY OF ICE-FREE LAKE (JULIAN DAY)
C****   4  LAST DAY OF ICED-UP LAKE  (JULIAN DAY)
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

      SUBROUTINE DIAGA
!@sum  DIAGA accumulate various diagnostics during dynamics
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : grav,rgas,kapa,lhe,sha,bygrav,bbyg,gbyrb,tf
     *     ,rvap,gamd
      USE MODEL_COM, only : im,imh,fim,byim,jm,jeq,lm,ls1,idacc,ptop
     *     ,pmtop,psfmpt,mdyn,mdiag,sig,sige,dsig,zatmo,WM,ntype,ftype
     *     ,u,v,t,p,q
      USE GEOM, only : areag,cosp,dlat,dxv,dxyn,dxyp,dxys,dxyv,dyp,fcor
     *     ,imaxj,ravpn,ravps,sinp,bydxyv
      USE DAGCOM, only : aj,areg,jreg,apj,ajl,asjl,ail,j50n,j70n,J5NUV
     *     ,J5SUV,J5S,J5N,aij,ij_dtdp,ij_dsev,ij_phi1k,ij_pres,ij_puq
     *     ,ij_pvq,ij_slp,ij_t850,ij_ujet,ij_vjet,j_tx1,j_tx,j_qp
     *     ,j_dtdjt,j_dtdjs,j_dtdgtr,j_dtsgst,j_rictr,j_rostr,j_ltro
     *     ,j_ricst,j_rosst,j_lstr,j_gamm,j_gam,j_gamc,lstr,il_ueq
     *     ,il_veq,il_weq,il_teq,il_qeq,il_w50n,il_t50n,il_u50n,il_w70n
     *     ,il_t70n,il_u70n,   kgz,pmb,ght,
     &     JL_DTDYN,JL_ZMFNTMOM,JL_TOTNTMOM,JL_APE,JL_UEPAC,
     &     JL_VEPAC,JL_UWPAC,JL_VWPAC,JL_WEPAC,JL_WWPAC,
     &     JL_EPFLXN,JL_EPFLXV
      USE DYNAMICS, only : pk,phi,pmid,plij, pit,sd
      USE RADNCB, only : rqt,lm_req
      USE PBLCOM, only : tsavg

      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: W
      COMMON/WORK2/W
      DOUBLE PRECISION, DIMENSION(LM) :: GMEAN
      DOUBLE PRECISION, DIMENSION(JM) :: TIL,UI,UMAX,PI,EL,RI,DUDVSQ
      DOUBLE PRECISION, DIMENSION(NTYPE,JM) :: SPTYPE
      DOUBLE PRECISION, DIMENSION(JM,LM) ::
     &     THJL,THSQJL,SPI,PHIPI,TPI,TJL0
      DOUBLE PRECISION, DIMENSION(JM,LM), SAVE :: TJL0
      DOUBLE PRECISION, DIMENSION(JM,LM-1) :: SDMEAN
      DOUBLE PRECISION, DIMENSION(IM,JM) :: PUV
      DOUBLE PRECISION, DIMENSION(LM_REQ) :: TRI
      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: TX
      DOUBLE PRECISION, DIMENSION(IM) :: THSEC,PSEC,SQRTP,PDA
      COMMON/WORK3/TX
      INTEGER, DIMENSION(LM) :: LUPA,LDNA
      CHARACTER*16 TITLE
      INTEGER :: IFIRST = 1
      DOUBLE PRECISION, PARAMETER :: ONE=1.,ZERO20=1.E-20,P1000=1000.
      INTEGER :: I,IM1,J,K,L,JET,JR,KM,LDN,LUP,
     &     IP1,LM1,LP1,LR,MBEGIN,
     &     I150E,I110W,I135W,IT
      DOUBLE PRECISION THBAR ! external
      DOUBLE PRECISION ::
     &     BBYGV,BDTDL,BYSDSG,CDTDL,DLNP,DLNP12,DLNP23,DBYSD,
     &     DLNS,DP,DS,DT2,DTHDP,DU,DUDP,DUDX,DV,DXYPJ,ELX,EPSLON,
     *     ESEPS,FPHI,GAMC,GAMM,GAMX,GMEANL,P1,P4,P4I,
     &     PDN,PE,PEQ,PEQM1,PEQM2,PHIRI,PIBYIM,PIJ,PITIJ,PITMN,
     *     PKE,PL,PRQ1,PRT,PU4I,PUP,PUV4I,PV4I,PVTHP,
     *     QLH,ROSSX,SDMN,SDPU,SMALL,SP,SP2,SS,T4,THETA,THGM,THMN,TPIL,
     *     TZL,UAMAX,UMN,UPE,VPE,X,Z4,THI

      REAL*8 QSAT

      CALL GETTIME(MBEGIN)
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
      BYSDSG=1./(1.-SIGE(LM+1))
      EPSLON=1.
      KM=0
      DO 5 K=1,KGZ
      IF (PMTOP.GT.PMB(K)) GO TO 6
    5 KM=KM+1
    6 CONTINUE
      I150E = IM*(180+150)/360+1   ! WEST EDGE OF 150 EAST
      I110W = IM*(180-110)/360+1   ! WEST EDGE OF 110 WEST
      I135W = IM*(180-135)/360+1   ! WEST EDGE OF 135 WEST
      PRQ1=.75*PMTOP
      DLNP12=LOG(.75/.35)
      DLNP23=LOG(.35/.1)
      DO 10 L=1,LM
      LUPA(L)=L+1
   10 LDNA(L)=L-1
      LDNA(1)=1
      LUPA(LM)=LM
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
      DO 84 L=1,LM
      TX(1,1,L)=T(1,1,L)*PK(L,1,1)
      TX(1,JM,L)=T(1,JM,L)*PK(L,1,JM)
      DO 70 I=2,IM
      T(I,1,L)=T(1,1,L)
      T(I,JM,L)=T(1,JM,L)
      TX(I,1,L)=TX(1,1,L)
   70 TX(I,JM,L)=TX(1,JM,L)
      DO 80 J=2,JM-1
      DO 80 I=1,IM
   80 TX(I,J,L)=T(I,J,L)*PK(L,I,J)
   84 CONTINUE
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
      DXYPJ=DXYP(J)
C**** NUMBERS ACCUMULATED FOR A SINGLE LEVEL
      PI(J)=0.
      SPTYPE(:,J)=0
      DO 120 I=1,IMAXJ(J)
      JR=JREG(I,J)
      DO IT=1,NTYPE
        SPTYPE(IT,J)=SPTYPE(IT,J)+FTYPE(IT,I,J)
        AJ(J,J_TX1,IT)=AJ(J,J_TX1,IT)+(TX(I,J,1)-TF)*FTYPE(IT,I,J)
      END DO
      AREG(JR,J_TX1)=AREG(JR,J_TX1)+(TX(I,J,1)-TF)*DXYPJ
      PI(J)=PI(J)+P(I,J)
      AIJ(I,J,IJ_PRES)=AIJ(I,J,IJ_PRES)+ P(I,J)
      AIJ(I,J,IJ_SLP)=AIJ(I,J,IJ_SLP)+((P(I,J)+PTOP)*(1.+BBYG*ZATMO(I,J)
     *     /TSAVG(I,J))**GBYRB-P1000)
  120 CONTINUE
      APJ(J,1)=APJ(J,1)+PI(J)
C**** CALCULATE GEOPOTENTIAL HEIGHTS AT SPECIFIC MILLIBAR LEVELS
      DO 180 I=1,IMAXJ(J)
      K=1
      L=1
  172 L=L+1
      PDN=PMID(L-1,I,J)
      PL=PMID(L,I,J)
      IF (PMB(K).LT.PL.AND.L.LT.LM) GO TO 172
      IF (ABS(TX(I,J,L)-TX(I,J,L-1)).LT.EPSLON) GO TO 176
      BBYGV=(TX(I,J,L-1)-TX(I,J,L))/(PHI(I,J,L)-PHI(I,J,L-1))
  174 AIJ(I,J,IJ_PHI1K-1+K)=AIJ(I,J,IJ_PHI1K-1+K)+(PHI(I,J,L)
     *     -TX(I,J,L)*((PMB(K)/PL)**(RGAS*BBYGV)-1.)/BBYGV-GHT(K)*GRAV)
      IF (K.EQ.2) AIJ(I,J,IJ_T850)=AIJ(I,J,IJ_T850)+(TX(I,J,L)-TF
     *     +(TX(I,J,L-1)-TX(I,J,L))*LOG(PMB(K)/PL)/LOG(PDN/PL))
      IF (K.GE.KM) GO TO 180
      K=K+1
      IF (PMB(K).LT.PL.AND.L.LT.LM) GO TO 172
      GO TO 174
 176  AIJ(I,J,IJ_PHI1K-1+K)=AIJ(I,J,IJ_PHI1K-1+K)+(PHI(I,J,L)
     *     -RGAS*TX(I,J,L)*LOG(PMB(K)/PL)-GHT(K)*GRAV)
      IF (K.EQ.2) AIJ(I,J,IJ_T850)=AIJ(I,J,IJ_T850)+(TX(I,J,L)-TF)
      IF (K.GE.KM) GO TO 180
      K=K+1
      IF (PMB(K).LT.PL.AND.L.LT.LM) GO TO 172
      GO TO 176
  180 CONTINUE
  190 CONTINUE
C**** ACCUMULATION OF TEMP., POTENTIAL TEMP., Q, AND RH
      DO 250 J=1,JM
      DXYPJ=DXYP(J)
      DO 230 L=1,LM
      TPI(J,L)=0.
      PHIPI(J,L)=0.
C     QPI=0.
      SPI(J,L)=0.
      THI=0.
C     RHPI=0.
      DBYSD=DSIG(L)*BYSDSG
      DO 220 I=1,IMAXJ(J)
      JR=JREG(I,J)
      PIJ=PLIJ(L,I,J)
      DO IT=1,NTYPE
        AJ(J,J_TX,IT)=AJ(J,J_TX,IT)+(TX(I,J,L)-TF)*FTYPE(IT,I,J)*DBYSD
        AJ(J,J_QP,IT)=AJ(J,J_QP,IT)+(Q(I,J,L)+WM(I,J,L))*PIJ*DSIG(L)*
     *       FTYPE(IT,I,J)
      END DO
      AREG(JR,J_QP)=AREG(JR,J_QP)+(Q(I,J,L)+WM(I,J,L))*PIJ*DSIG(L)*DXYPJ
      AREG(JR,J_TX)=AREG(JR,J_TX)+(TX(I,J,L)-TF)*DBYSD*DXYPJ
      TPI(J,L)=TPI(J,L)+(TX(I,J,L)-TF)*PIJ
      PHIPI(J,L)=PHIPI(J,L)+PHI(I,J,L)*PIJ
C     QPI=QPI+Q(I,J,L)*PIJ
      SPI(J,L)=SPI(J,L)+T(I,J,L)*PIJ
      THI=THI+T(I,J,L)
C     QLH=LHE
C     QSATL=QSAT(TX(I,J,L),QLH,SIG(L)*PIJ+PTOP)
C     IF (QSATL.GT.1.) QSATL=1.
C     RHPI=RHPI+Q(I,J,L)*PIJ/QSATL
  220 CONTINUE
C     AJL(J,L,JL_01)=AJL(J,L,JL_01)+TPI(J,L)
C     AJL(J,L,JL_02)=AJL(J,L,JL_02)+PHIPI(J,L)
C     AJL(J,L,JL_03)=AJL(J,L,JL_03)+QPI
      AJL(J,L,JL_DTDYN)=AJL(J,L,JL_DTDYN)+THI-TJL0(J,L)
C     AJL(J,L,JL_18)=AJL(J,L,JL_18)+RHPI
  230 CONTINUE
  250 CONTINUE
C****
C**** NORTHWARD GRADIENT OF TEMPERATURE: TROPOSPHERIC AND STRATOSPHERIC
C****
      DO 385 J=2,JM-1
C**** MEAN TROPOSPHERIC NORTHWARD TEMPERATURE GRADIENT
      DO 340 L=1,LS1-1
      DO 335 I=1,IM
        DO IT=1,NTYPE
          AJ(J,J_DTDJT,IT)=AJ(J,J_DTDJT,IT)+(TX(I,J+1,L)-TX(I,J-1,L))
     *         *FTYPE(IT,I,J)*DSIG(L)
        END DO
  335 CONTINUE
 340  CONTINUE
C**** MEAN STRATOSPHERIC NORTHWARD TEMPERATURE GRADIENT
      IF (LS1.GT.LM) GO TO 380
      DO 370 L=LS1,LM
      DO 350 I=1,IM
        DO IT=1,NTYPE
          AJ(J,J_DTDJS,IT)=AJ(J,J_DTDJS,IT)+(TX(I,J+1,L)-TX(I,J-1,L))
     *         *FTYPE(IT,I,J)*DSIG(L)
        END DO
  350 CONTINUE
 370  CONTINUE
  380 CONTINUE
  385 CONTINUE
C****
C**** STATIC STABILITIES: TROPOSPHERIC AND STRATOSPHERIC
C****
      DO 490 J=1,JM
      DXYPJ=DXYP(J)
C**** OLD TROPOSPHERIC STATIC STABILITY
      DO 390 I=1,IMAXJ(J)
      JR=JREG(I,J)
      SS=(T(I,J,LS1-1)-T(I,J,1))/(PHI(I,J,LS1-1)-PHI(I,J,1)+ZERO20)
      DO IT=1,NTYPE
        AJ(J,J_DTDGTR,IT)=AJ(J,J_DTDGTR,IT)+SS*FTYPE(IT,I,J)
      END DO
      AREG(JR,J_DTDGTR)=AREG(JR,J_DTDGTR)+SS*DXYPJ
  390 AIJ(I,J,IJ_DTDP)=AIJ(I,J,IJ_DTDP)+SS
C**** OLD STRATOSPHERIC STATIC STABILITY
      DO 440 I=1,IMAXJ(J)
      JR=JREG(I,J)
      SS=(T(I,J,LM)-T(I,J,LS1-1))/((PHI(I,J,LM)-PHI(I,J,LS1-1))+ZERO20)
      DO IT=1,NTYPE
        AJ(J,J_DTSGST,IT)=AJ(J,J_DTSGST,IT)+SS*FTYPE(IT,I,J)
      END DO
      AREG(JR,J_DTSGST)=AREG(JR,J_DTSGST)+SS*DXYPJ
  440 CONTINUE
C****
C**** NUMBERS ACCUMULATED FOR THE RADIATION EQUILIBRIUM LAYERS
C****
      DO 470 LR=1,LM_REQ
      TRI(LR)=0.
      DO 460 I=1,IMAXJ(J)
  460 TRI(LR)=TRI(LR)+RQT(LR,I,J)
  470 ASJL(J,LR,1)=ASJL(J,LR,1)+(TRI(LR)-TF*IMAXJ(J))
      PHIRI=0.
      DO 480 I=1,IMAXJ(J)
  480 PHIRI=PHIRI+(PHI(I,J,LM)+RGAS*.5*(TX(I,J,LM)+RQT(1,I,J))
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
      DO J=2,JM
        DUDVSQ(J)=0.
        UMAX(J)=0.
        DO I=1,IM
          DU=U(I,J,LS1-1)-U(I,J,1)
          DV=V(I,J,LS1-1)-V(I,J,1)
          DUDVSQ(J)=DUDVSQ(J)+(DU*DU+DV*DV)*PUV(I,J)
        END DO
      END DO
      DO J=2,JM-1
        PIBYIM=PI(J)*BYIM
        DLNP=LOG((SIG(1)*PIBYIM+PTOP)/(SIG(LS1-1)*PIBYIM+PTOP))
        DLNS=LOG(SPI(J,LS1-1)/SPI(J,1))
        DS=SPI(J,LS1-1)-SPI(J,1)
        EL(J)=SQRT(DLNS/DLNP)
        RI(J)=DS*DLNP/(.5*(DUDVSQ(J)+DUDVSQ(J+1)))
      END DO
      DO L=1,LS1-1
        DO J=2,JM
          UI(J)=0.
          DO I=1,IM
            UI(J)=UI(J)+U(I,J,L)
          END DO
        END DO
        DO J=2,JM-1
          UAMAX=ABS(UI(J)+UI(J+1))
          IF (UAMAX.GT.UMAX(J)) UMAX(J)=UAMAX
        END DO
      END DO
      DO J=2,JM-1
        ROSSX=DYP(J)/(DXYP(J)*SINP(J))
        ELX=1./SINP(J)
        DO IT=1,NTYPE
          AJ(J,J_RICTR,IT)=AJ(J,J_RICTR,IT)+RI(J)  *SPTYPE(IT,J)
          AJ(J,J_ROSTR,IT)=AJ(J,J_ROSTR,IT)+UMAX(J)*SPTYPE(IT,J)*ROSSX
          AJ(J,J_LTRO ,IT)=AJ(J,J_LTRO ,IT)+EL(J)  *SPTYPE(IT,J)*ELX
        END DO
      END DO
C**** NUMBERS ACCUMULATED OVER THE LOWER STRATOSPHERE
C**** LSTR is approx. 10mb level. This maintains consistency over
C**** the different model tops
CNOST IF (LS1.GT.LSTR) SKIP THIS PART !NEEDED FOR NO-STRATOSPHERE MODELS
      DO J=2,JM
        DUDVSQ(J)=0.
        UMAX(J)=0.
        DO I=1,IM
          DU=U(I,J,LSTR)-U(I,J,LS1-1)
          DV=V(I,J,LSTR)-V(I,J,LS1-1)
          DUDVSQ(J)=DUDVSQ(J)+(DU*DU+DV*DV)*PUV(I,J)
        END DO
      END DO
      DO J=2,JM-1
        PIBYIM=PI(J)*BYIM
        DLNP=LOG((SIG(LS1-1)*PIBYIM+PTOP)/(SIG(LSTR)*PSFMPT+PTOP))
        DLNS=LOG(SPI(J,LSTR)/SPI(J,LS1-1))
        DS=SPI(J,LSTR)-SPI(J,LS1-1)
        EL(J)=SQRT(DLNS/DLNP)
        RI(J)=DS*DLNP/(.5*(DUDVSQ(J)+DUDVSQ(J+1)))
      END DO
      DO L=LS1,LSTR
        DO J=2,JM
          UI(J)=0.
          DO I=1,IM
            UI(J)=UI(J)+U(I,J,L)
          END DO
        END DO
        DO J=2,JM-1
          UAMAX=ABS(UI(J)+UI(J+1))
          IF (UAMAX.GT.UMAX(J)) UMAX(J)=UAMAX
        END DO
      END DO
      DO J=2,JM-1
        ROSSX=DYP(J)/(DXYP(J)*SINP(J))
        ELX=1./SINP(J)
        DO IT=1,NTYPE
          AJ(J,J_RICST,IT)=AJ(J,J_RICST,IT)+RI(J)  *SPTYPE(IT,J)
          AJ(J,J_ROSST,IT)=AJ(J,J_ROSST,IT)+UMAX(J)*SPTYPE(IT,J)*ROSSX
          AJ(J,J_LSTR ,IT)=AJ(J,J_LSTR ,IT)+EL(J)  *SPTYPE(IT,J)*ELX
        END DO
      END DO
C****
C**** MEAN TROPOSPHERIC LAPSE RATES:  MOIST CONVECTIVE, ACTUAL,
C****    DRY ADIABATIC
C****
      X=RGAS*LHE*LHE/(SHA*RVAP)
      DO 570 J=1,JM
      GAMM=0.
      DO 560 L=1,LS1-1
      TZL=TPI(J,L)/PI(J)+TF
      PRT=(SIG(L)*PI(J)*BYIM+PTOP)*RGAS*TZL
      ESEPS=QSAT(TZL,LHE,ONE)
      GAMM=GAMM+(PRT+LHE*ESEPS)/(PRT+X*ESEPS/TZL)*DSIG(L)
  560 CONTINUE
      GAMX=(TPI(J,1)-TPI(J,LS1-1))/(PHIPI(J,LS1-1)-PHIPI(J,1))
      DO IT=1,NTYPE
        AJ(J,J_GAMM,IT)=AJ(J,J_GAMM,IT)+GAMM*SPTYPE(IT,J)
        AJ(J,J_GAM ,IT)=AJ(J,J_GAM ,IT)+GAMX*SPTYPE(IT,J)
      END DO
  570 CONTINUE
C**** DRY ADIABATIC LAPSE RATE
      DO 580 J=1,JM
      TPIL=0.
      DO 575 L=1,LS1-1
  575 TPIL=TPIL+TPI(J,L)*DSIG(L)
      TIL(J)=TPIL/(PI(J)*(SIGE(1)-SIGE(LS1)))
  580 CONTINUE
      DO 590 J=2,JM-1
      X=SINP(J)*GRAV/(COSP(J)*RGAS*2.*DLAT)
      DT2=TIL(J+1)-TIL(J-1)
      GAMC=GAMD+X*DT2/(TIL(J)+TF)
      DO IT=1,NTYPE
        AJ(J,J_GAMC,IT)=AJ(J,J_GAMC,IT)+GAMC*SPTYPE(IT,J)
      END DO
  590 CONTINUE
C****
C**** EASTWARD TRANSPORTS
C****
      I=IM
      DO 600 L=1,LM
      DO 600 J=2,JM-1
      DO 600 IP1=1,IM
      AIJ(I,J,IJ_PUQ)=AIJ(I,J,IJ_PUQ)+(PLIJ(L,I,J)+PLIJ(L,IP1,J))*
     *   (U(I,J,L)+U(I,J+1,L))*(Q(I,J,L)+Q(IP1,J,L))*DSIG(L)
  600 I=IP1
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
      P4=PLIJ(L,I,J-1)+PLIJ(L,IP1,J-1)+PLIJ(L,I,J)+PLIJ(L,IP1,J)
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
      AIJ(I,J,IJ_DSEV)=AIJ(I,J,IJ_DSEV)+P4*(SHA*T4+Z4)*V(I,J,L)*DSIG(L)
     *     *DXV(J)
      SP2=PLIJ(L,IP1,J-1)+PLIJ(L,IP1,J)
      AIJ(IP1,J,IJ_PVQ)=AIJ(IP1,J,IJ_PVQ)+SP2
     *  *(V(I,J,L)+V(IP1,J,L))*(Q(IP1,J-1,L)+Q(IP1,J,L))*DSIG(L)
  620 I=IP1
C     AJL(J,L,JL_04)=AJL(J,L,JL_04)+PU4I
C     AJL(J,L,JL_05)=AJL(J,L,JL_05)+PV4I
C     AJL(J,L,JL_14)=AJL(J,L,JL_14)+(PU4I*PU4I+PV4I*PV4I)/P4I
C     AJL(J,L,JL_15)=AJL(J,L,JL_15)+PWW4I
C     AJL(J,L,JL_20)=AJL(J,L,JL_20)+PT16I*PV4I/P4I
C     AJL(J,L,JL_21)=AJL(J,L,JL_21)+PTV16I
C     AJL(J,L,JL_22)=AJL(J,L,JL_22)+PZ16I*PV4I/P4I
C     AJL(J,L,JL_23)=AJL(J,L,JL_23)+PZV16I
C     AJL(J,L,JL_24)=AJL(J,L,JL_24)+PQ16I*PV4I/P4I
C     AJL(J,L,JL_25)=AJL(J,L,JL_25)+PQV16I
C     AJL(J,L,JL_26)=AJL(J,L,JL_26)+PWW4I*PV4I/P4I
C     AJL(J,L,JL_27)=AJL(J,L,JL_27)+PWWV4I
      AJL(J,L,JL_ZMFNTMOM)=AJL(J,L,JL_ZMFNTMOM)+PU4I*PV4I/P4I
      AJL(J,L,JL_TOTNTMOM)=AJL(J,L,JL_TOTNTMOM)+PUV4I
  640 CONTINUE
C****
C**** EVEN LEVEL GEOPOTENTIALS, VERTICAL WINDS AND VERTICAL TRANSPORTS
C****
      DO 655 J=1,JM
C     PITI=0.
C     DO 648 I=1,IMAXJ(J)
C 648 PITI=PITI+PIT(I,J)
      DO 655 L=1,LM-1
C     SDI=0.
C     PZI=0.
C     SDZI=0.
C     PDSE2I=0.
C     SDDS2I=0.
C     PQ2I=0.
C     SDQ2I=0.
      DO 650 I=1,IMAXJ(J)
C     SDI=SDI+SD(I,J,L)
      PIJ=PLIJ(L,I,J)
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
C     AJL(J,L,JL_06)=AJL(J,L,JL_06)+SDI+DSIG(L+1)*PITI
C     AJL(J,L,JL_34)=AJL(J,L,JL_34)+(SDZI-PZI*SDI/PI(J))
C     AJL(J,L,JL_30)=AJL(J,L,JL_30)+PDSE2I*SDI/PI(J)
C     AJL(J,L,JL_31)=AJL(J,L,JL_31)+SDDS2I
C     AJL(J,L,JL_32)=AJL(J,L,JL_32)+PQ2I*SDI/PI(J)
C     AJL(J,L,JL_33)=AJL(J,L,JL_33)+SDQ2I
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
C     AJL(J,L,JL_35)=AJL(J,L,JL_35)+TKET
C     AJL(J,L,JL_EPFLXV)=
C    *   AJL(J,L,JL_EPFLXV)+(UT-2.*UM*(SDMEAN(J,L)+SDMEAN(J-1,L)))
C     AJL(J,L,JL_EPFLXN)=
C    *  AJL(J,L,JL_EPFLXN)+(UT+4*AMA*FIM*(SDMEAN(J,L)+SDMEAN(J-1,L)))
C 670 CONTINUE
C****
C**** AVAILABLE POTENTIAL ENERGY
C****
C**** SET UP FOR CALCULATION
      DO 710 L=1,LM
  710 GMEAN(L)=0.
      DO 740 J=1,JM
      DO 720 I=1,IMAXJ(J)
  720 SQRTP(I)=SQRT(P(I,J))
C**** GMEAN CALCULATED FOR EACH LAYER, THJL, THSQJL ARRAYS FILLED
      DO 730 L=1,LM
      LDN=LDNA(L)
      LUP=LUPA(L)
      THJL(J,L)=0.
      THSQJL(J,L)=0.
      DO 730 I=1,IMAXJ(J)
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
  760 AJL(J,L,JL_APE)=AJL(J,L,JL_APE)+
     &        (THSQJL(J,L)-2.*THJL(J,L)*THGM+THGM*THGM*
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
C 770 AJL(J,L,JL_07)=AJL(J,L,JL_07)-(PWAI-PWI*SPAI/PI(J))
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
C 774 AJL(J,L,JL_53)=AJL(J,L,JL_53)+VDVTI+UDUTI
C    *   -(UPDA4I*DUTI+VPDA4I*DVTI)/PDA4I
C        IF (JM.NE.24) GO TO 850
C****
C**** CERTAIN HORIZONTAL WIND AVERAGES
C****
      DO L=1,LM
      DO J=2,JM
      DO I=I135W,I110W      ! EAST PACIFIC
        AJL(J,L,JL_UEPAC)=AJL(J,L,JL_UEPAC)+U(I,J,L)
        AJL(J,L,JL_VEPAC)=AJL(J,L,JL_VEPAC)+V(I,J,L)
      END DO
      DO I=I150E,IM             ! WEST PACIFIC
        AJL(J,L,JL_UWPAC)=AJL(J,L,JL_UWPAC)+U(I,J,L)
        AJL(J,L,JL_VWPAC)=AJL(J,L,JL_VWPAC)+V(I,J,L)
      END DO
      END DO
      DO J=J5SUV,J5NUV
      DO I=1,IM
        AIL(I,L,IL_UEQ)=AIL(I,L,IL_UEQ)+U(I,J,L)
        AIL(I,L,IL_VEQ)=AIL(I,L,IL_VEQ)+V(I,J,L)
      END DO
      END DO
      DO J=J5S,J5N
      DO I=1,IM
        AIL(I,L,IL_TEQ)=AIL(I,L,IL_TEQ)+(TX(I,J,L)-TF)
        AIL(I,L,IL_QEQ)=AIL(I,L,IL_QEQ)+Q(I,J,L)/QSAT(TX(I,J,L),LHE
     *       ,PMID(L,I,J))
      END DO
      END DO
      DO I=1,IM
        AIL(I,L,IL_T50N)=AIL(I,L,IL_T50N)+(TX(I,J50N,L)-TF)
        AIL(I,L,IL_U50N)=AIL(I,L,IL_U50N)+(U(I,J50N,L)+U(I,J50N+1,L))
        AIL(I,L,IL_T70N)=AIL(I,L,IL_T70N)+(TX(I,J70N,L)-TF)
        AIL(I,L,IL_U70N)=AIL(I,L,IL_U70N)+(U(I,J70N,L)+U(I,J70N+1,L))
      END DO
      END DO
C****
C**** CERTAIN VERTICAL WIND AVERAGES
C****
      DO L=1,LM-1
      DO J=2,JM-1
        DO I=I135W,I110W        ! EAST PACIFIC
          AJL(J,L,JL_WEPAC)=AJL(J,L,JL_WEPAC)+W(I,J,L)
        END DO
        DO I=I150E,IM           ! WEST PACIFIC
          AJL(J,L,JL_WWPAC)=AJL(J,L,JL_WWPAC)+W(I,J,L)
        END DO
      END DO
      DO I=1,IM
        DO J=J5S,J5N        ! +/- 5 DEG (APPROX.)
          AIL(I,L,IL_WEQ) =AIL(I,L,IL_WEQ)+W(I,J,L)
        END DO
        AIL(I,L,IL_W50N)=AIL(I,L,IL_W50N)+W(I,J50N,L)
        AIL(I,L,IL_W70N)=AIL(I,L,IL_W70N)+W(I,J70N,L)
      END DO
      END DO
C****
C**** ELIASSEN PALM FLUX
C****
C**** NORTHWARD TRANSPORT
      DO 868 J=2,JM
      I=IM
      DO 862 IP1=1,IM
      PDA(I)=.5*((P(I,J)+P(IP1,J))*DXYS(J)+(P(I,J-1)+P(IP1,J-1))*
     *  DXYN(J-1))
      PSEC(I)=PDA(I)*BYDXYV(J)
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
  868 AJL(J,L,JL_EPFLXN)=AJL(J,L,JL_EPFLXN)+FPHI
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
c      IF (DTHDP.LT.SMALL) WRITE (6,999) J,L,DTHDP,SMALL
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
      DP=(SIG(L)-SIG(L+1))*P(I,J)
      IF(L.GE.LS1) DP=(SIG(L)-SIG(L+1))*PSFMPT
      IF(L.EQ.LS1-1) DP=P(I,J)*SIG(L)-PSFMPT*SIG(LS1)
      PVTHP=PVTHP+DP*VPE*(T(I,J,L)+T(I,J,L+1)-THMN)
      PITIJ=PIT(I,J)
      IF(L.GE.LS1-1) PITIJ=0.
      SDPU=SDPU+(SD(I,J,L)-SDMN+(PITIJ-PITMN)*SIGE(L+1))*UPE
  874 IM1=I
      AJL(J,L,JL_EPFLXV)=AJL(J,L,JL_EPFLXV)+
     &     (.5*FIM*FCOR(J)-.25*DUDX)*PVTHP/DTHDP + SDPU
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
C 895 AJL(J,L,JL_52)=AJL(J,L,JL_52)+PV2I
C****
C**** TRANSFORMED STREAM FUNCTION
C****
C---- skip for now lines 1132-1163
C**** ACCUMULATE TIME USED IN DIAGA

      CALL TIMEOUT(MBEGIN,MDIAG,MDYN)
      RETURN
  999 FORMAT (' DTHETA/DP IS TOO SMALL AT J=',I4,' L=',I4,2F15.6)

      ENTRY DIAGA0
C****
C**** INITIALIZE TJL0 ARRAY (FROM PRIOR TO ADVECTION)
C****
      DO L=1,LM
      DO J=1,JM
        THI=0.
        DO I=1,IMAXJ(J)
          THI=THI+T(I,J,L)
        END DO
        TJL0(J,L)=THI
      END DO
      END DO
      RETURN
C****
      END SUBROUTINE DIAGA

      SUBROUTINE DIAGB
C****
C**** CONTENTS OF AJK(J,K,N)  (SUM OVER LONGITUDE AND TIME OF)
C**CP   1  DP  (PDN-PM(K+1);  PDN=MAX(PM(K+1),PS)
C**CP   2  DP4  (PDN-PM(K+1))   (UV GRID)
C***1   3  (TX-TF)*DP
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
      USE CONSTANT, only : lhe,omega,sha,tf
      USE MODEL_COM, only :
     &     im,imh,fim,byim,jm,jeq,lm,ls1,idacc,ptop,psfmpt,
     &     mdyn,mdiag, ndaa,      !! skipse,
     &     sig,sige,dsig, Jhour
     *     ,u,v,t,p,q,wm
      USE GEOM, only :
     &     COSV,DXV,DXYN,DXYP,DXYS,DXYV,DYP,DYV,FCOR,IMAXJ,RADIUS
      USE DAGCOM, only : ajk,aijk,aijl,ajlsp,speca,adiurn,nspher,
     &     nwav_dag,ndiupt,hr_in_day,ijk_u,ijk_v,ijk_t,ijk_q,ijk_dp
     *     ,ijk_dse,klayer,idd_w,ijdd,ijl_u, ijl_v, ijl_dse, ijl_dp,
     *     ijl_q,
     &      JK_DPA,JK_DPB,JK_TEMP,JK_HGHT,JK_Q,JK_THETA,
     &      JK_RH,JK_U,JK_V,JK_ZMFKE,JK_TOTKE,JK_ZMFNTSH,
     &      JK_TOTNTSH,JK_ZMFNTGEO,JK_TOTNTGEO,JK_ZMFNTLH,
     &      JK_TOTNTLH,JK_ZMFNTKE,JK_TOTNTKE,JK_ZMFNTMOM,JK_TOTNTMOM,
     &      JK_P2KEDPGF,JK_DPSQR,JK_NPTSAVG,
     &      JK_VVEL,JK_ZMFVTDSE,JK_TOTVTDSE,JK_ZMFVTLH,JK_TOTVTLH,
     &      JK_VTGEOEDDY,JK_BAREKEGEN,JK_POTVORT,JK_VTPV,
     &      JK_VTPVEDDY,JK_NPTSAVG1,JK_TOTVTKE,
     &      JK_VTAMEDDY,JK_TOTVTAM,JK_SHETH,JK_DUDTMADV,JK_DTDTMADV,
     &      JK_DUDTTEM,JK_DTDTTEM,JK_EPFLXNCP,JK_EPFLXVCP,
     &      JK_UINST,JK_TOTDUDT,JK_TINST,
     &      JK_TOTDTDT,JK_EDDVTPT,JK_CLDH2O
      USE DYNAMICS, only : phi,dut,dvt,plij
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION, DIMENSION(IMH+1,NSPHER) :: KE
      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: W,ZX,STB
      COMMON/WORK2/W,ZX,STB
      DOUBLE PRECISION, DIMENSION(JM,LM) ::
     &     STJK,DPJK,UJK,VJK,WJK,TJK,
     &     PSIJK,UP,TY,PSIP,WTJK,UVJK,WUJK
      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: TX ! from diaga
      COMMON/WORK3/TX
      DOUBLE PRECISION, DIMENSION(IM) :: PSEC,X1
      DOUBLE PRECISION, DIMENSION(LM) :: SHETH,PMO,PLO,DPM,DTH
      DOUBLE PRECISION, DIMENSION(LM+1) :: PM,PL
      DOUBLE PRECISION, DIMENSION(0:IMH,JM,LM):: FCJKA,FCJKB,FCPVA,FCPVB
      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: AIJL1,AIJL2,AIJL3,AIJL4

      INTEGER ::
     &     I,IH,IM1,INCH,INCHM,IP1,IZERO,J,J45N,
     &     JET,JHEMI,K,KDN,KM,KMM1,KR,KS,KS1,KSPHER,KUP,KX,L,LMP1,
     &     LUP,MBEGIN,N,NM

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

      INTEGER :: IFIRST = 1
      DOUBLE PRECISION, PARAMETER :: ZERO20=1.E-20,BIG=1.E20
      DOUBLE PRECISION :: QSAT

      CALL GETTIME(MBEGIN)
      IF (IFIRST.NE.1) GO TO 50
      IFIRST=0
C**** INITIALIZE CERTAIN QUANTITIES
      LMP1=LM+1
      JET=LS1-1
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
      DO 170 K=1,KM
      DPI=0.
      TPI=0.
      PHIPI=0.
      QPI=0.
      WMPI=0.
      THPI=0.
      RHPI=0.
      FIMI=0.
      DO 160 I=1,IMAXJ(J)
C**** FIND L=L(K) AND LUP=L(K+1) S.T. P(LUP).GT.P(K+1)
      SP=PLIJ(K,I,J)
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
      TPI=TPI+(TX(I,J,L)-TF)*DP
      PHIPI=PHIPI+PHI(I,J,L)*DP
      QPI=QPI+Q(I,J,L)*DP
      WMPI=WMPI+WM(I,J,L)*DP
CW       IF(WMPI.GT.1.E-3) WRITE(6,169) I,J,L,DP,WM(I,J,L),WMPI
CW169 FORMAT(1X,'1616--',3I5,3E15.2)
      THPI=THPI+T(I,J,L)*DP
      QSATL=QSAT(TX(I,J,L),QLH,SIG(L)*SP+PTOP)
      IF (QSATL.GT.1.) QSATL=1.
      RHPI=RHPI+Q(I,J,L)*DP/QSATL
      IF (L.EQ.LUP) GO TO 160
      L=L+1
      PDN=SP*SIGE(L)+PTOP
      GO TO 150
  160 CONTINUE
      AJK(J,K,JK_NPTSAVG1)=AJK(J,K,JK_NPTSAVG1)+FIMI
      AJK(J,K,JK_DPA)=AJK(J,K,JK_DPA)+DPI
      AJK(J,K,JK_TEMP)=AJK(J,K,JK_TEMP)+TPI
      AJK(J,K,JK_HGHT)=AJK(J,K,JK_HGHT)+PHIPI
      AJK(J,K,JK_Q)=AJK(J,K,JK_Q)+QPI
      AJK(J,K,JK_THETA)=AJK(J,K,JK_THETA)+THPI
      AJK(J,K,JK_RH)=AJK(J,K,JK_RH)+RHPI
      AJK(J,K,JK_CLDH2O)=AJK(J,K,JK_CLDH2O)+WMPI
         TJK(J,K)=THPI/(DPI+ZERO20)
         IF (IDACC(4).EQ.1) AJK(J,K,JK_TINST)=TJK(J,K)
         AJK(J,K,JK_TOTDTDT)=TJK(J,K)-AJK(J,K,JK_TINST)
  170 CONTINUE
C****
C**** CALCULATE STABILITY AT ODD LEVELS ON PU GRID
C****
      DO 230 J=1,JM
      I=IMAXJ(J)
      DO 230 IP1=1,IMAXJ(J)
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
      DO 260 K=1,KM
      STJK(J,K)=0.
      DPJK(J,K)=0.
      I=IMAXJ(J)
      DO 250 IP1=1,IMAXJ(J)
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
      AIJL(I,J,L,IJL_U)=AIJL(I,J,L,IJL_U)+DELP*U(I,J,L)
      AIJL(I,J,L,IJL_V)=AIJL(I,J,L,IJL_V)+DELP*V(I,J,L)
      AIJL(I,J,L,IJL_DSE)=AIJL(I,J,L,IJL_DSE)+SHA*PT4L+PZ4L
      AIJL(I,J,L,IJL_Q)=AIJL(I,J,L,IJL_Q)+PQ4L
      AIJL(I,J,L,IJL_DP)=AIJL(I,J,L,IJL_DP)+DELP
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
!!    IF(SKIPSE.EQ.1.) GO TO 334
      AIJK(I,J,K,IJK_U)  =AIJK(I,J,K,IJK_U)  +PUK
      AIJK(I,J,K,IJK_V)  =AIJK(I,J,K,IJK_V)  +PVK
      AIJK(I,J,K,IJK_DSE)=AIJK(I,J,K,IJK_DSE)+SHA*PT4K+PZ4K
      AIJK(I,J,K,IJK_DP) =AIJK(I,J,K,IJK_DP) +DPK
      AIJK(I,J,K,IJK_T)  =AIJK(I,J,K,IJK_T)  +PT4K
      AIJK(I,J,K,IJK_Q)  =AIJK(I,J,K,IJK_Q)  +PQ4K
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
      AJK(J,K,JK_DPB)=AJK(J,K,JK_DPB)+DPI
      AJK(J,K,JK_DPSQR)=AJK(J,K,JK_DPSQR)+DPSQI
      AJK(J,K,JK_NPTSAVG)=AJK(J,K,JK_NPTSAVG)+FIMI
      IF (DPI.LT.ZERO20) DPI=ZERO20
      AJK(J,K,JK_U)=AJK(J,K,JK_U)+PUI
      AJK(J,K,JK_V)=AJK(J,K,JK_V)+PVI
      AJK(J,K,JK_ZMFKE)=AJK(J,K,JK_ZMFKE)+(PUI*PUI+PVI*PVI)/DPI
      AJK(J,K,JK_TOTKE)=AJK(J,K,JK_TOTKE)+PWWI
      AJK(J,K,JK_ZMFNTSH)=AJK(J,K,JK_ZMFNTSH)+PT4I*PVI/DPI
      AJK(J,K,JK_TOTNTSH)=AJK(J,K,JK_TOTNTSH)+PTV4I
      AJK(J,K,JK_ZMFNTGEO)=AJK(J,K,JK_ZMFNTGEO)+PZ4I*PVI/DPI
      AJK(J,K,JK_TOTNTGEO)=AJK(J,K,JK_TOTNTGEO)+PZV4I
      AJK(J,K,JK_ZMFNTLH)=AJK(J,K,JK_ZMFNTLH)+PQ4I*PVI/DPI
      AJK(J,K,JK_TOTNTLH)=AJK(J,K,JK_TOTNTLH)+PQV4I
      AJK(J,K,JK_ZMFNTKE)=AJK(J,K,JK_ZMFNTKE)+PWWI*PVI/DPI
      AJK(J,K,JK_TOTNTKE)=AJK(J,K,JK_TOTNTKE)+PWWVI
      AJK(J,K,JK_ZMFNTMOM)=AJK(J,K,JK_ZMFNTMOM)+PUI*PVI/DPI
      AJK(J,K,JK_TOTNTMOM)=AJK(J,K,JK_TOTNTMOM)+PUVI
      AJK(J,K,JK_P2KEDPGF)=AJK(J,K,JK_P2KEDPGF)+VDVTI+UDUTI-
     *   (PUI*DUTI+PVI*DVTI)/DPI
      SHETH(K)=(PSV4I-PS4I*PVI/DPI)*DXYV(J)/(STJK(J-1,K)*DXYN(J-1)+
     *   STJK(J,K)*DXYS(J))
         UJK(J,K)=PUI/DPI
         VJK(J,K)=PVI/DPI
         PSIJK(J,K)=.25*SHETH(K)/DPI
         UVJK(J,K)=(PUVI-PUI*PVI/DPI)/DPI
         IF (IDACC(4).EQ.1) AJK(J,K,JK_UINST)=UJK(J,K)
         AJK(J,K,JK_TOTDUDT)=UJK(J,K)-AJK(J,K,JK_UINST)
  350 AJK(J,K,JK_SHETH)=AJK(J,K,JK_SHETH)+SHETH(K)
      DO 345 L=1,LM
C**** SPECTRAL ANALYSIS OF DRY STATIC ENERGY FLUX, LATENT HEAT FLUX,
C**** AND ANGULAR MOMENTUM FLUX
      CALL FFT(AIJL2(1,J,L),FCPVA(0,J,L),FCPVB(0,J,L))
      CALL FFT(AIJL3(1,J,L),FCJKA(0,J,L),FCJKB(0,J,L))
      DO 342 N=0,NWAV_DAG
  342 AJLSP(J,L,N,1)=AJLSP(J,L,N,1)+.5*FIM*(FCPVA(N,J,L)*
     *  FCJKA(N,J,L)+FCPVB(N,J,L)*FCJKB(N,J,L))
      CALL FFT(AIJL4(1,J,L),FCJKA(0,J,L),FCJKB(0,J,L))
      DO 343 N=0,NWAV_DAG
  343 AJLSP(J,L,N,2)=AJLSP(J,L,N,2)+.5*FIM*(FCPVA(N,J,L)*
     *  FCJKA(N,J,L)+FCPVB(N,J,L)*FCJKB(N,J,L))
      CALL FFT(AIJL1(1,J,L),FCJKA(0,J,L),FCJKB(0,J,L))
      DO 344 N=0,NWAV_DAG
  344 AJLSP(J,L,N,3)=AJLSP(J,L,N,3)+.5*FIM*(FCPVA(N,J,L)*
     *  FCJKA(N,J,L)+FCPVB(N,J,L)*FCJKB(N,J,L))
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
      DO 558 J=1,JM
      DO 558 I=1,IM
      DO KR=1,NDIUPT
         IF(I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
C**** Warning:     This diagnostic has 3 flaws   (?)
C****          1 - It assumes that DTsrc=1hr, (DTsrc=3600.)
C****          2 - since DTdaa-Ndaa*DTsrc=2*DTdyn rather than 0,
C****              some hours are skipped once in a while
C****          3 - Some of the first Ndaa hours are skipped at the
C****              beginning of a month and overcounted at the end;
C****              this happens to balance out, if and only if
C****              mod(days_in_month,ndaa)=0  (i.e. February if Ndaa=7)
            IH=JHOUR+1
            DO INCH=1,NDAA
              IF(IH.GT.HR_IN_DAY) IH=IH-HR_IN_DAY
              ADIURN(IH,IDD_W,KR)=ADIURN(IH,IDD_W,KR)+1.E5*W(I,J,3)
     *             /DXYP(J)
              IH=IH+1
            END DO
         END IF
      END DO
  558 CONTINUE
      DO 565 J=1,JM
      DO 565 K=1,KM
      WI=0.
      DO 562 I=1,IMAXJ(J)
  562 WI=WI+W(I,J,K)
  565 AJK(J,K,JK_VVEL)=AJK(J,K,JK_VVEL)+WI
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
      DO 600 I=1,IMAXJ(J)
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
      AJK(J,K-1,JK_ZMFVTDSE)=AJK(J,K-1,JK_ZMFVTDSE)+
     &     BYFIM*(SHA*TKI+ZKI)*WI
      AJK(J,K-1,JK_TOTVTDSE)=AJK(J,K-1,JK_TOTVTDSE)+SHA*WTI+WZI
      AJK(J,K-1,JK_ZMFVTLH)=AJK(J,K-1,JK_ZMFVTLH)+BYFIM*QKI*WI
      AJK(J,K-1,JK_TOTVTLH)=AJK(J,K-1,JK_TOTVTLH)+WQI
      AJK(J,K-1,JK_VTGEOEDDY)=AJK(J,K-1,JK_VTGEOEDDY)+WZI-BYFIM*WI*ZKI
C     AJK(J,K-1,JK_BAREKEGEN)=AJK(J,K-1,JK_BAREKEGEN)+WTI-BYFIM*WI*TKI
         WJK(J,K)=BYFIM*WI/DXYP(J)
         WTJK(J,K)=BYFIM*(WTHI-BYFIM*WI*THKI)/DXYP(J)
         AJK(J,K-1,JK_EDDVTPT)=AJK(J,K-1,JK_EDDVTPT)+WTJK(J,K)
  610 CONTINUE
C****
C**** BAROCLINIC EDDY KINETIC ENERGY GENERATION
C****
      DO 630 J=1,JM
      DO 630 K=1,KM
      FIMI=0.
      W2I=0.
      PAI=0.
      WPA2I=0.
      DO 626 I=1,IMAXJ(J)
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
  630 AJK(J,K,JK_BAREKEGEN)=AJK(J,K,JK_BAREKEGEN)-
     &     (WPA2I-W2I*PAI/(FIMI+ZERO20))
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
      AJK(J,K-1,JK_TOTVTKE)=AJK(J,K-1,JK_TOTVTKE)+WKE4I
      AJK(J,K-1,JK_VTAMEDDY)=AJK(J,K-1,JK_VTAMEDDY)+WU4I-BYFIM*W4I*UKI
  710 AJK(J,K-1,JK_TOTVTAM)=AJK(J,K-1,JK_TOTVTAM)+WU4I+W4I*UEARTH
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
  730 AJK(J,K,JK_POTVORT)=AJK(J,K,JK_POTVORT)+PVI
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
      AJK(J,K-1,JK_VTPV)=AJK(J,K-1,JK_VTPV)+WPV4I
  760 AJK(J,K-1,JK_VTPVEDDY)=AJK(J,K-1,JK_VTPVEDDY)+
     &     WPV4I-W2I*PV2I/(FIMI+ZERO20)
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
      AJK(J,K,JK_EPFLXNCP)=AJK(J,K,JK_EPFLXNCP)+
     &     PSIJK(J,K)*(UJK(J,KUP)-UJK(J,KDN))/
     *  (PMO(KUP)-PMO(KDN))-UVJK(J,K)
  780 CONTINUE
      DO 800 J=2,JM-1
      DO 800 K=2,KMM1
      UY=(UJK(J+1,K)*DXV(J+1)-UJK(J,K)*DXV(J)-FCOR(J))/DXYP(J)
      PSIY=(PSIJK(J+1,K)*DXV(J+1)-PSIJK(J,K)*DXV(J))/DXYP(J)
C**** ZONAL MEAN MOMENTUM EQUATION   (MEAN ADVECTION)
      AJK(J,K,JK_DUDTMADV)=AJK(J,K,JK_DUDTMADV)-
     &     .5*UY*(VJK(J,K)+VJK(J+1,K))-
     *  .25*((UP(J+1,K+1)+UP(J,K+1))*WJK(J,K+1)+(UP(J+1,K)+UP(J,K))*
     *   WJK(J,K))
C**** ZONAL MEAN HEAT EQUATION   (MEAN ADVECTION)
      AJK(J,K,JK_DTDTMADV)=AJK(J,K,JK_DTDTMADV)-
     &     .5*(TY(J,K)*VJK(J,K)+TY(J+1,K)*VJK(J+1,K))
     *  -.5*STJK(J,K)*(WJK(J,K+1)+WJK(J,K))
C**** LAGRANGIAN MEAN MOMENTUM EQUATION  (MEAN ADVECTION)
      VSTAR=.5*(VJK(J,K)+VJK(J+1,K)-.5*(PSIP(J,K)+PSIP(J,K+1)
     *  +PSIP(J+1,K)+PSIP(J+1,K+1)))
      WSTAR=.5*(WJK(J,K)+WJK(J,K+1))+PSIY
      AJK(J,K,JK_DUDTTEM)=AJK(J,K,JK_DUDTTEM)-
     &     UY*VSTAR-.25*(UP(J,K)+UP(J+1,K)+
     *  UP(J,K+1)+UP(J+1,K+1))*WSTAR
      AJK(J,K,JK_DTDTTEM)=AJK(J,K,JK_DTDTTEM)-
     &     .5*(TY(J+1,K)+TY(J,K))*VSTAR-
     *  STJK(J,K)*WSTAR
C**** VERTICAL E-P FLUX
      AJK(J,K-1,JK_EPFLXVCP)=AJK(J,K-1,JK_EPFLXVCP)-
     &     WUJK(J,K)-.5*PSIJK(J,K)*UY
      AJK(J,K,JK_EPFLXVCP)=AJK(J,K,JK_EPFLXVCP)-.5*PSIJK(J,K)*UY
  800 CONTINUE
C****
C**** SPECTRAL ANALYSIS OF KINETIC ENERGIES AT CONSTANT PRESSURE
C****
      IZERO=0
      NM=1+IM/2
      J45N=2+.75*(JM-1.)
c      KS1=LS1
C**** TOTAL THE KINETIC ENERGIES
      KE(:,:)=0.
      DO 2140 J=2,JM
      I=IM
      DO 2020 IP1=1,IM
      PSEC(I)=.25*(P(I,J-1)+P(IP1,J-1)+P(I,J)+P(IP1,J))
 2020 I=IP1
      DO 2140 K=1,KM
        KSPHER=KLAYER(K)
c      KSPHER=2
c      IF (K.GE.KS1) KSPHER=1
      IF (J.GT.JEQ) KSPHER=KSPHER+1
      DO 2140 KX=IZERO,LM,LM
      DO 2090 I=1,IM
      DPUV=0.
      SP=PSEC(I)
      DO 2025 L=1,LS1-1
      PLO(L)=SP*SIG(L)+PTOP                       ! PL or PLO ??
 2025 PL(L)=SP*SIGE(L)+PTOP                       ! PLE or PL ??
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
 2110 KE(N,KSPHER+2)=KE(N,KSPHER+2)+X1(N)*DXYV(J)
      GO TO 2140
 2120 DO 2130 N=1,NM
      KE(N,KSPHER+2)=KE(N,KSPHER+2)+X1(N)*DXYV(J)
      KE(N,KSPHER)=KE(N,KSPHER)+.5D0*X1(N)*DXYV(J)
 2130 KE(N,KSPHER+1)=KE(N,KSPHER+1)+.5D0*X1(N)*DXYV(J)
 2140 CONTINUE
      DO 2150 KS=1,NSPHER
      DO 2150 N=1,NM
 2150 SPECA(N,18,KS)=SPECA(N,18,KS)+KE(N,KS)
C**** ACCUMULATE TIME USED IN DIAGA
      CALL TIMEOUT(MBEGIN,MDIAG,MDYN)
      RETURN
      END SUBROUTINE DIAGB

      SUBROUTINE DIAG7A
C****
C**** THIS SUBROUTINE ACCUMULATES A TIME SEQUENCE FOR SELECTED
C**** QUANTITIES AND FROM THAT PRINTS A TABLE OF WAVE FREQUENCIES.
C****
      USE CONSTANT, only : grav,bygrav
      USE MODEL_COM, only : im,imh,jm,lm,
     &     IDACC,JEQ,LS1,MDIAG,P,PTOP,PSFMPT,SIG,SIGE,U,V
      USE DYNAMICS, only : PHI
      USE DAGCOM, only : nwav_dag,wave,max12hr_sequ,j50n
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(0:IMH) :: AN,BN
      INTEGER, PARAMETER :: KM=6,KQMAX=12
      INTEGER :: NMAX=nwav_dag
      DOUBLE PRECISION, DIMENSION(IM,KM) :: HTRD
      DOUBLE PRECISION, DIMENSION(KM), PARAMETER ::
     &     PMB=(/922.,700.,500.,300.,100.,10./),
     &     GHT=(/500.,2600.,5100.,8500.,15400.,30000./)
      DOUBLE PRECISION :: PIJ50N,PL,PLE,PLM1,SLOPE
      INTEGER I,IDACC9,JLK,K,KQ,L,LX,MNOW,N
      INTEGER, SAVE, DIMENSION(KM) :: JLKDEX
      INTEGER, SAVE :: L300,L50,L850
      INTEGER, SAVE :: IFIRST = 1

      IDACC9=IDACC(9)+1
      IDACC(9)=IDACC9
      IF (IDACC9.GT.Max12HR_sequ) RETURN

      IF (IFIRST.EQ.1) THEN
        IFIRST=0
        L850=LM
        L300=LM
        L50=LM
        DO L=2,LM
          LX=LM+1-L
          PLE=.25*(SIGE(LX)+2.*SIGE(LX+1)+SIGE(LX+2))*PSFMPT+PTOP
          IF (PLE.LT.850.) L850=LX
          IF (PLE.LT.300.) L300=LX
          IF (PLE.LT.50.) L50=LX
        END DO
        WRITE (6,889) L850,L300,L50
 889    FORMAT (' LEVELS FOR WIND WAVE POWER DIAG  L850=',I3,
     *       ' L300=',I3,' L50=',I3)
        JLKDEX(1)=JEQ+JM*(L850-1)
        JLKDEX(2)=JEQ+JM*(L850-1+LM)
        JLKDEX(3)=JEQ+JM*(L300-1)
        JLKDEX(4)=JEQ+JM*(L300-1+LM)
        JLKDEX(5)=JEQ+JM*(L50-1)
        JLKDEX(6)=JEQ+JM*(L50-1+LM)
      END IF

      DO KQ=1,5,2
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
 130    L=L+1
        IF(L.GE.LS1) PIJ50N=PSFMPT
        PLM1=PL
        PL=SIG(L)*PIJ50N+PTOP
        IF (PMB(K).LT.PL.AND.L.LT.LM) GO TO 130
C**** ASSUME THAT PHI IS LINEAR IN LOG P
        SLOPE=(PHI(I,J50N,L-1)-PHI(I,J50N,L))/LOG(PLM1/PL)
 140    HTRD(I,K)=(PHI(I,J50N,L)+SLOPE*LOG(PMB(K)/PL))*BYGRAV-GHT(K)
        IF (K.GE.KM) GO TO 150
        K=K+1
        IF (PMB(K).LT.PL.AND.L.LT.LM) GO TO 130
        GO TO 140
 150  CONTINUE
      DO KQ=7,KQMAX
        CALL FFT(HTRD(1,KQ-6),AN,BN)
        DO N=1,NMAX
          WAVE(1,IDACC9,N,KQ)=AN(N)
          WAVE(2,IDACC9,N,KQ)=BN(N)
        END DO
      END DO
      CALL TIMER (MNOW,MDIAG)
      RETURN
      END SUBROUTINE DIAG7A

      SUBROUTINE DIAGCA (M)
!@sum  DIAGCA Keeps track of the conservation properties of angular
!@+    momentum, kinetic energy, mass, total potential energy and water
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : mdiag
      USE DAGCOM, only : icon_AM,icon_KE,icon_MS,icon_TPE,icon_WM
     *     ,icon_LKM,icon_LKE,icon_EWM,icon_WTG,icon_HTG,icon_MSI
     *     ,icon_HSI,icon_SSI,title_con
      USE SOIL_DRV, only: conserv_WTG,conserv_HTG
      IMPLICIT NONE
!@var M index denoting from where DIAGCA is called
      INTEGER, INTENT(IN) :: M
C****
C**** THE PARAMETER M INDICATES WHEN DIAGCA IS BEING CALLED
C**** M=1  INITIALIZE CURRENT QUANTITY
C****   2  AFTER DYNAMICS
C****   3  AFTER CONDENSATION
C****   4  AFTER RADIATION
C****   5  AFTER PRECIPITATION
C****   6  AFTER LAND SURFACE (INCL. RIVER RUNOFF)
C****   7  AFTER FULL SURFACE INTERACTION
C****   8  AFTER FILTER
C****   9  AFTER STRATOSPHERIC DRAG
C****  10  AFTER OCEAN DYNAMICS
C****  11  AFTER OCEAN SUB-GRIDSCALE PHYS
C****  12  AFTER DAILY
C****
      REAL*8, EXTERNAL :: conserv_AM,conserv_KE,conserv_MS,conserv_PE
     *     ,conserv_WM,conserv_EWM,conserv_LKM,conserv_LKE
     *     ,conserv_MSI,conserv_HSI,conserv_SSI
      REAL*8 MNOW

C**** ATMOSPHERIC ANGULAR MOMENTUM
      CALL conserv_DIAG(M,conserv_AM,icon_AM)

C**** ATMOSPHERIC KINETIC ENERGY
      CALL conserv_DIAG(M,conserv_KE,icon_KE)

C**** ATMOSPHERIC MASS
      CALL conserv_DIAG(M,conserv_MS,icon_MS)

C**** ATMOSPHERIC TOTAL POTENTIAL ENERGY
      CALL conserv_DIAG(M,conserv_PE,icon_TPE)

C**** ATMOSPHERIC TOTAL WATER MASS
      CALL conserv_DIAG(M,conserv_WM,icon_WM)

C**** ATMOSPHERIC TOTAL WATER ENERGY
      CALL conserv_DIAG(M,conserv_EWM,icon_EWM)

C**** LAKE MASS AND ENERGY
      CALL conserv_DIAG(M,conserv_LKM,icon_LKM)
      CALL conserv_DIAG(M,conserv_LKE,icon_LKE)

C**** SEAICE MASS, ENERGY, SALT
      CALL conserv_DIAG(M,conserv_MSI,icon_MSI)
      CALL conserv_DIAG(M,conserv_HSI,icon_HSI)
      CALL conserv_DIAG(M,conserv_SSI,icon_SSI)

C**** GROUND WATER AND ENERGY
      CALL conserv_DIAG(M,conserv_WTG,icon_WTG)
      CALL conserv_DIAG(M,conserv_HTG,icon_HTG)

C**** OCEAN CALLS ARE DEALT WITH SEPERATELY
      CALL DIAGCO (M)
C****
      CALL TIMER (MNOW,MDIAG)
      RETURN
      END SUBROUTINE DIAGCA

      SUBROUTINE DIAGCD (M,DT1,UX,VX,DUT,DVT,PIT)
!@sum  DIAGCD Keeps track of the conservation properties of angular
!@+    momentum and kinetic energy inside dynamics routines
!@auth Gary Russell
!@ver  1.0
      USE CONSTANT, only : omega,mb2kg
      USE MODEL_COM, only : im,jm,lm,fim,mdiag,mdyn
      USE GEOM, only : cosv,radius,ravpn,ravps
      USE DAGCOM, only : consrv
      IMPLICIT NONE
C****
C**** THE PARAMETER M INDICATES WHEN DIAGCD IS BEING CALLED
C**** M=1  AFTER ADVECTION IN DYNAMICS
C****   2  AFTER CORIOLIS FORCE IN DYNAMICS
C****   3  AFTER PRESSURE GRADIENT FORCE IN DYNAMICS
C****
!@var M index denoting from where DIAGCD is called
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
      INTEGER :: I,J,L,MBEGIN
      DOUBLE PRECISION :: DUTI,DUTIL,RKEI,RKEIL

      CALL GETTIME(MBEGIN)
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
     *         (PI(J-1)*RAVPN(J-1)+PI(J)*RAVPS(J)))*COSV(J)*RADIUS*mb2kg
          CONSRV(J,13)=CONSRV(J,13)+RKEIL*mb2kg
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
          CONSRV(J,2*M-1)=CONSRV(J,2*M-1)+DUTIL*COSV(J)*RADIUS*mb2kg
          CONSRV(J,2*M+10)=CONSRV(J,2*M+10)+RKEIL*mb2kg
        END DO
      END IF
C****
      CALL TIMEOUT(MBEGIN,MDIAG,MDYN)
      RETURN
      END SUBROUTINE DIAGCD

      SUBROUTINE conserv_DIAG (M,CONSFN,ICON)
!@sum  conserv_DIAG generic routine keeps track of conserved properties
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : jm
      USE DAGCOM, only : consrv,nofm
      IMPLICIT NONE
!@var M index denoting from where routine is called
      INTEGER, INTENT(IN) :: M
!@var ICON index for the quantity concerned
      INTEGER, INTENT(IN) :: ICON
!@var CONSFN external routine that calculates total conserved quantity
      EXTERNAL CONSFN
!@var TOTAL amount of conserved quantity at this time
      REAL*8, DIMENSION(JM) :: TOTAL
      INTEGER :: I,J,NM,NI

C**** NOFM contains the indexes of the CONSRV array where each
C**** change is to be stored for each quantity. If NOFM(M,ICON)=0,
C**** no calculation is done.
C**** NOFM(1,ICON) is the index for the instantaneous value.
      IF (NOFM(M,ICON).gt.0) THEN
C**** Calculate current value TOTAL
        CALL CONSFN(TOTAL)
        NM=NOFM(M,ICON)
        NI=NOFM(1,ICON)
C**** Accumulate difference from last time in CONSRV(NM)
        IF (M.GT.1) THEN
          DO J=1,JM
            CONSRV(J,NM)=CONSRV(J,NM)+(TOTAL(J)-CONSRV(J,NI))
          END DO
        END IF
C**** Save current value in CONSRV(NI)
        DO J=1,JM
          CONSRV(J,NI)=TOTAL(J)
        END DO
      END IF
      RETURN
C****
      END SUBROUTINE conserv_DIAG

      SUBROUTINE conserv_AM(AM)
!@sum  conserv_AM calculates total atmospheric angular momentum
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : omega,radius,mb2kg
      USE MODEL_COM, only : im,jm,lm,fim,ls1,dsig,p,u,psfmpt,pstrat
      USE GEOM, only : cosv,dxyn,dxys,dxyv
      IMPLICIT NONE

      REAL*8, DIMENSION(JM) :: AM,PI
      INTEGER :: I,IP1,J,L
      DOUBLE PRECISION :: UMI,UMIL
C****
C**** ANGULAR MOMENTUM
C****
      AM(1)=0.
      PI(1)=FIM*P(1,1)
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
     *             +(P(I,J)+P(IP1,J))*DXYS(J))
            ELSE
              UMI=UMI+U(I,J,L)*(2.*PSFMPT*DXYV(J))
            END IF
            I=IP1
          END DO
          UMIL=UMIL+UMI*DSIG(L)
        END DO
        AM(J)=(RADIUS*OMEGA*COSV(J)*((PI(J-1)*DXYN(J-1)+PI(J)*DXYS(J))
     *       +FIM*PSTRAT*DXYV(J))+.5*UMIL)*COSV(J)*RADIUS*mb2kg
      END DO
      RETURN
C****
      END SUBROUTINE conserv_AM

      SUBROUTINE conserv_KE(RKE)
!@sum  conserv_KE calculates total atmospheric kinetic energy
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : mb2kg
      USE MODEL_COM, only : im,jm,lm,fim,dsig,ls1,p,u,v,psfmpt
      USE GEOM, only : dxyn,dxys,dxyv
      IMPLICIT NONE

      REAL*8, DIMENSION(JM) :: RKE
      INTEGER :: I,IP1,J,L
      DOUBLE PRECISION :: RKEI,RKEIL
C****
C**** KINETIC ENERGY
C****
      RKE(1)=0.
      DO J=2,JM
        RKEIL=0.
        DO L=1,LM
          RKEI=0.
          I=IM
          DO IP1=1,IM
            IF(L.LT.LS1) THEN
              RKEI=RKEI+(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))
     *             *((P(I,J-1)+P(IP1,J-1))*DXYN(J-1)+(P(I,J)+P(IP1,J
     *             ))*DXYS(J))
            ELSE
              RKEI=RKEI+(U(I,J,L)*U(I,J,L)+V(I,J,L)*V(I,J,L))*
     *             (2.*PSFMPT*DXYV(J))
            END IF
            I=IP1
          END DO
          RKEIL=RKEIL+RKEI*DSIG(L)
        END DO
        RKE(J)=RKEIL*mb2kg
      END DO
      RETURN
C****
      END SUBROUTINE conserv_KE

      SUBROUTINE conserv_MS(RMASS)
!@sum  conserv_MA calculates total atmospheric mass
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : mb2kg
      USE MODEL_COM, only : im,jm,fim,p,pstrat
      IMPLICIT NONE
      REAL*8, DIMENSION(JM) :: RMASS
      INTEGER :: I,J
C****
C**** MASS
C****
      RMASS(1) =FIM*(P(1,1) +PSTRAT)*mb2kg
      RMASS(JM)=FIM*(P(1,JM)+PSTRAT)*mb2kg
      DO J=2,JM-1
        RMASS(J)=FIM*PSTRAT
        DO I=1,IM
          RMASS(J)=RMASS(J)+P(I,J)
        END DO
        RMASS(J)=RMASS(J)*mb2kg
      END DO
      RETURN
C****
      END SUBROUTINE conserv_MS

      SUBROUTINE conserv_PE(TPE)
!@sum  conserv_TPE calculates total atmospheric potential energy
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : sha,mb2kg
      USE MODEL_COM, only : im,jm,lm,fim,dsig,ls1,t,p,ptop,psfmpt
     *     ,zatmo
      USE GEOM, only : imaxj
      USE DYNAMICS, only : pk
      IMPLICIT NONE
      REAL*8, DIMENSION(JM) :: TPE
      INTEGER :: I,J,L
      DOUBLE PRECISION :: TPEI,TPEIL,SGEOI
C****
C**** TOTAL POTENTIAL ENERGY (J/m^2)
C****
 400  DO J=1,JM
        TPEIL=0.
        DO L=1,LM
          TPEI=0.
          DO I=1,IMAXJ(J)
            IF(L.LT.LS1) THEN
              TPEI=TPEI+T(I,J,L)*PK(L,I,J)*P(I,J)
            ELSE
              TPEI=TPEI+T(I,J,L)*PK(L,I,J)*PSFMPT
            END IF
          END DO
          TPEIL=TPEIL+TPEI*DSIG(L)
        END DO
        SGEOI=0.
        DO I=1,IMAXJ(J)
          SGEOI=SGEOI+ZATMO(I,J)*(P(I,J)+PTOP)
        END DO
        TPE(J)=(SGEOI+TPEIL*SHA)*mb2kg
      END DO
      TPE(1)=FIM*TPE(1)
      TPE(JM)=FIM*TPE(JM)
      RETURN
C****
      END SUBROUTINE conserv_PE

      SUBROUTINE conserv_WM(WATER)
!@sum  conserv_WM calculates total atmospheric water mass
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : mb2kg
      USE MODEL_COM, only : im,jm,lm,fim,wm,q
      USE GEOM, only : imaxj
      USE DYNAMICS, only : pdsig
      IMPLICIT NONE

      REAL*8, DIMENSION(JM) :: WATER
      INTEGER :: I,J,L
C****
C**** TOTAL WATER MASS (kg/m^2)
C****
      DO J=1,JM
        WATER(J) = 0.
        DO L=1,LM
          DO I=1,IMAXJ(J)
            WATER(J)=WATER(J)+(Q(I,J,L)+WM(I,J,L))*PDSIG(L,I,J)*mb2kg
          END DO
        END DO
      END DO
      WATER(1) = FIM*WATER(1)
      WATER(JM)= FIM*WATER(JM)
      RETURN
C****
      END SUBROUTINE conserv_WM

      SUBROUTINE conserv_EWM(EWATER)
!@sum  conserv_EWM calculates total atmospheric water energy
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : mb2kg,shv,grav
      USE MODEL_COM, only : im,jm,lm,fim,wm,t,q,p
      USE GEOM, only : imaxj
      USE DYNAMICS, only : pdsig, pmid, pk
      IMPLICIT NONE
      REAL*8, PARAMETER :: HSCALE = 7.8 ! km ?????
      REAL*8, DIMENSION(JM) :: EWATER
      INTEGER :: I,J,L
      REAL*8 W
C****
C**** TOTAL WATER ENERGY (J/m^2)
C****
      DO J=1,JM
        EWATER(J) = 0.
        DO L=1,LM
          DO I=1,IMAXJ(J)
            W = (Q(I,J,L)+WM(I,J,L))*PDSIG(L,I,J)*mb2kg
c this calculation needs to be checked!
            EWATER(J)=EWATER(J)+SHV*W*T(I,J,L)*PK(L,I,J)+W*GRAV*HSCALE
     *           *LOG(P(I,J)/PMID(L,I,J))
          END DO
        END DO
      END DO
      EWATER(1) = FIM*EWATER(1)
      EWATER(JM)= FIM*EWATER(JM)
      RETURN
C****
      END SUBROUTINE conserv_EWM

      SUBROUTINE DIAG5D (M5,NDT,DUT,DVT)
      USE MODEL_COM, only : im,imh,jm,lm,fim,
     &     DSIG,JEQ,LS1,MDIAG,MDYN
      USE DAGCOM, only : speca,nspher,klayer
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: DUT,DVT

      DOUBLE PRECISION, DIMENSION(0:IMH,JM,LM,2) :: FCUVA,FCUVB
      COMMON/WORK7/FCUVA,FCUVB

      INTEGER :: M5,NDT

      DOUBLE PRECISION, DIMENSION(IMH+1) :: X
      DOUBLE PRECISION, DIMENSION(0:IMH) :: FA,FB
      DOUBLE PRECISION, DIMENSION(IMH+1,NSPHER) :: KE

      INTEGER :: J,J45N,KUV,KSPHER,L,MBEGIN,MKE,N,NM

      NM=1+IM/2
      J45N=2.+.75*(JM-1.)
      MKE=M5

      GO TO (810,810,810,100,100,  100,810),M5
  810 WRITE (6,910) M5
  910 FORMAT ('0INCORRECT VALUE OF M5 WHEN CALLING DIAG5D.  M5=',I5)
      STOP 29
C****
C**** KINETIC ENERGY
C****
C**** TRANSFER RATES FOR KINETIC ENERGY IN THE DYNAMICS
  100 CALL GETTIME(MBEGIN)
      KE(:,:)=0.

      DO 170 L=1,LM
        KSPHER=KLAYER(L)
c      KSPHER=2
c      IF (L.GE.LS1) KSPHER=1
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
         KE(:,KSPHER+2)=KE(:,KSPHER+2)+X(:)*DSIG(L)
      ELSE IF (J.EQ.JEQ) THEN ! EQUATOR
        DO N=1,NM
        KE(N,KSPHER+2)=KE(N,KSPHER+2)+     X(N)*DSIG(L)
        KE(N,KSPHER)  =KE(N,KSPHER)  +.5D0*X(N)*DSIG(L) ! CONTRIB TO SH
        KE(N,KSPHER+1)=KE(N,KSPHER+1)+.5D0*X(N)*DSIG(L) ! CONTRIB TO NH
        ENDDO
        IF (KUV.EQ.2) KSPHER=KSPHER+1
      ENDIF
  170 CONTINUE

      DO 180 KSPHER=1,NSPHER
      DO 180 N=1,NM
  180 SPECA(N,MKE,KSPHER)=SPECA(N,MKE,KSPHER)+KE(N,KSPHER)/NDT
C****
      CALL TIMEOUT(MBEGIN,MDIAG,MDYN)
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
      USE CONSTANT, only : sha
      USE MODEL_COM, only : im,imh,jm,lm,fim,
     &     DSIG,IDACC,JEQ,LS1,MDIAG,
     &     P,PTOP,PSFMPT,SIG,T,U,V,ZATMO
      USE GEOM, only : AREAG,DXYN,DXYP,DXYS
      USE DAGCOM, only : speca,atpe,nspher,kspeca,klayer
      USE DYNAMICS, only : sqrtp,pk
      IMPLICIT NONE
      SAVE
      INTEGER :: M5,NDT

      DOUBLE PRECISION, DIMENSION(IM) :: X
      DOUBLE PRECISION, DIMENSION(IMH+1,NSPHER) :: KE,APE
      DOUBLE PRECISION, DIMENSION(IMH+1,4) :: VAR
      DOUBLE PRECISION, DIMENSION(2) :: TPE
      DOUBLE PRECISION, DIMENSION(IM,JM) :: SQRTM
      DOUBLE PRECISION, DIMENSION(LM) :: THJSP,THJNP,THGM

      INTEGER, PARAMETER :: IZERO=0

      INTEGER, DIMENSION(KSPECA), PARAMETER ::
     &     MTPEOF=(/0,0,1,0,0,0,0,2,0,3,  4,0,5,0,6,0,7,0,0,8/)

      INTEGER ::
     &     I,IJL2,IP1,J,J45N,
     &     JH,JHEMI,JP,K,KS,KSPHER,L,LDN,
     &     LUP,MAPE,MKE,MNOW,MTPE,N,
     &     NM

      DOUBLE PRECISION ::
     &     GMEAN,GMSUM,SQRTPG,SUMI,SUMT,THGSUM,THJSUM

      INTEGER :: IFIRST = 1

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
C**** Note: KSPHER has been re-arranged from previous models to better
C****       deal with optionally resolved stratosphere. The higher
C****       values are only used if the model top is high enough.
C****
C**** KSPHER=1 SOUTHERN TROPOSPHERE         2 NORTHERN TROPOSPHERE
C****        3 EQUATORIAL TROPOSPHERE       4 45 DEG NORTH TROPOSPHERE
C****
C****        5 SOUTHERN LOW STRATOSPHERE    6 NORTHERN LOW STRATOSPHERE
C****        7 EQUATORIAL LOW STRATOSPHERE  8 45 DEG NORTH LOW STRATOSPH
C****
C****        9 SOUTHERN MID STRATOSPHERE   10 NORTHERN MID STRATOSPHERE
C****       11 EQUATORIAL MID STRATOSPHERE 12 45 DEG NORTH MID STRATOSPH
C****
C****       13 SOUTHERN UPP STRATOSPHERE   14 NORTHERN UPP STRATOSPHERE
C****       15 EQUATORIAL UPP STRATOSPHERE 16 45 DEG NORTH UPP STRATOSPH
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
        KSPHER=KLAYER(L)
c      KSPHER=2
c      IF (L.GE.LS1) KSPHER=1
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
  222 KE(N,KSPHER+2)=KE(N,KSPHER+2)+X(N)*DSIG(L)
      GO TO 240
  225 DO 230 N=1,NM
      KE(N,KSPHER+2)=KE(N,KSPHER+2)+X(N)*DSIG(L)
      KE(N,KSPHER)=KE(N,KSPHER)+.5D0*X(N)*DSIG(L)
  230 KE(N,KSPHER+1)=KE(N,KSPHER+1)+.5D0*X(N)*DSIG(L)
      IF (K.EQ.LM) KSPHER=KSPHER+1
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
  350 CONTINUE

      VAR(2:NM,1:2)=0.
      VAR(1,1)=.5*(THJSP(L)-THGM(L))**2*DXYP(1)*FIM
      VAR(1,2)=.5*(THJNP(L)-THGM(L))**2*DXYP(JM)*FIM
      GMEAN=((THJSP(LUP)-THJSP(LDN))*DXYP(1)*(SIG(L)*P(1,1)+PTOP)/
     *  (SQRTP(1,1)*P(1,1)*PK(L,1,1))+(THJNP(LUP)-THJNP(LDN))*DXYP(JM)*
     *  (SIG(L)*P(1,JM)+PTOP)/(SQRTP(1,JM)*P(1,JM)*PK(L,1,JM)))*FIM
      JHEMI=1
      DO 388 J=2,JM-1
        GMSUM=0.
        DO I=1,IM
          X(I)=T(I,J,L)*SQRTP(I,J)-THGM(L)
          GMSUM=GMSUM+(T(I,J,LUP)-T(I,J,LDN))*(SIG(L)*P(I,J)+PTOP)/
     *         (P(I,J)*PK(L,I,J))
        END DO
        GMEAN=GMEAN+GMSUM*DXYP(J)
        CALL FFTE (X,X)
        DO N=1,NM
          VAR(N,JHEMI)=VAR(N,JHEMI)+X(N)*DXYP(J)
        END DO
        IF (J.NE.JEQ-1) GO TO 384
        DO N=1,NM
          VAR(N,3)=X(N)*DXYP(J)
        END DO
        JHEMI=2
 384    IF (J.NE.J45N-1) GO TO 388
        DO N=1,NM
          VAR(N,4)=X(N)*DXYP(J)
        END DO
 388  CONTINUE
      GMEAN=DSIG(L)*AREAG*(SIG(LDN)-SIG(LUP))/GMEAN
      KS=KLAYER(L)
c      KS=2
c      IF (L.GE.LS1) KS=1
      DO JHEMI=1,4
        DO N=1,NM
          APE(N,KS)=APE(N,KS)+VAR(N,JHEMI)*GMEAN
        END DO
        KS=KS+1
      END DO
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
      IF (M5.EQ.2) THEN
C**** ACCUMULATE MEAN KINETIC ENERGY AND MEAN POTENTIAL ENERGY
        DO KS=1,NSPHER
        DO N=1,NM
          SPECA(N,2,KS)=SPECA(N,2,KS)+KE(N,KS)
          SPECA(N,3,KS)=SPECA(N,3,KS)+APE(N,KS)
        END DO
        END DO
        ATPE(1,1)=ATPE(1,1)+TPE(1)
        ATPE(1,2)=ATPE(1,2)+TPE(2)
      END IF
      CALL TIMER (MNOW,MDIAG)
      RETURN
      END SUBROUTINE DIAG5A
C****
      SUBROUTINE DIAG5F(UX,VX)
C**** FOURIER COEFFICIENTS FOR CURRENT WIND FIELD
C****
      USE MODEL_COM, only : im,imh,jm,lm,
     &     IDACC,MDIAG,MDYN
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: UX,VX
      DOUBLE PRECISION, DIMENSION(0:IMH,JM,LM,2) :: FCUVA,FCUVB
      COMMON/WORK7/FCUVA,FCUVB
      INTEGER :: J,L,MBEGIN

      CALL GETTIME(MBEGIN)
      IDACC(6)=IDACC(6)+1
      DO L=1,LM
         DO J=2,JM
            CALL FFT(UX(1,J,L),FCUVA(0,J,L,1),FCUVB(0,J,L,1))
            CALL FFT(VX(1,J,L),FCUVA(0,J,L,2),FCUVB(0,J,L,2))
         ENDDO
      ENDDO
      CALL TIMEOUT(MBEGIN,MDIAG,MDYN)

      RETURN
      END SUBROUTINE DIAG5F

      SUBROUTINE DIAG4A
C****
C**** THIS SUBROUTINE PRODUCES A TIME HISTORY OF ENERGIES
C****
      USE MODEL_COM, only : im,jm,lm,
     &     IDACC,JEQ,LS1,ISTRAT          !! ,SKIPSE
      USE GEOM, only : DXYV
      USE DAGCOM, only : energy,speca,ajk,aijk,ijk_u,ijk_v,ijk_dp,ned
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
!!    IF (SKIPSE.EQ.1.) GO TO 540
C**** CALCULATE CURRENT SEKE   ! this diagnostic makes no sense
c      BYIADA=1./IDACC(4)
c      DO 530 L=1,LM
c      KS=5
c      IF (L.GE.LS1) KS=15
c      DO 530 J=2,JM
c      IF (AJK(J,L,JK_DPB).LE.1.D-20) GO TO 530
c      PU4TI=0.
c      PV4TI=0.
c      SKE4I=0.
c      DO 510 I=1,IM
c      PU4TI=PU4TI+AIJK(I,J,L,IJK_U)
c      PV4TI=PV4TI+AIJK(I,J,L,IJK_V)
c  510 SKE4I=SKE4I+(AIJK(I,J,L,IJK_U)*AIJK(I,J,L,IJK_U)
c     *            +AIJK(I,J,L,IJK_V)*AIJK(I,J,L,IJK_V))
c     *            /(AIJK(I,J,L,IJK_DP)+1.D-20)
c      SEKE=(SKE4I-(PU4TI*PU4TI+PV4TI*PV4TI)/AJK(J,L,JK_DPB))*
c     *   DXYV(J)*BYIADA
c      IF (J.EQ.JEQ) GO TO 520
c      ENERGY(KS,IDACC5)=ENERGY(KS,IDACC5)+SEKE
c      GO TO 530
c  520 ENERGY(KS,IDACC5)=ENERGY(KS,IDACC5)+.5*SEKE
c      ENERGY(KS+1,IDACC5)=ENERGY(KS+1,IDACC5)+.5*SEKE
c      KS=KS+1
c  530 CONTINUE
C**** OTHER ENERGIES COME FROM LATEST SPECTRAL ANALYSIS
c  540 CONTINUE
      DO I=0,1+ISTRAT  ! loop over number of 'spheres'
        ENERGY(1+NED*I,IDACC5)=SPECA(1,19,1+4*I)   ! SH
        ENERGY(2+NED*I,IDACC5)=SPECA(1,19,2+4*I)   ! NH
        ENERGY(5+NED*I,IDACC5)=SPECA(2,19,2+4*I)   ! NH wave 1
        ENERGY(6+NED*I,IDACC5)=SPECA(3,19,2+4*I)   ! NH wave 2
        ENERGY(7+NED*I,IDACC5)=SPECA(1,20,1+4*I)
        ENERGY(8+NED*I,IDACC5)=SPECA(1,20,2+4*I)
        DO N=2,NM
        ENERGY( 3+NED*I,IDACC5)=ENERGY( 3+10*I,IDACC5)+SPECA(N,19,1+4*I)
        ENERGY( 4+NED*I,IDACC5)=ENERGY( 4+10*I,IDACC5)+SPECA(N,19,2+4*I)
        ENERGY( 9+NED*I,IDACC5)=ENERGY( 9+10*I,IDACC5)+SPECA(N,20,1+4*I)
        ENERGY(10+NED*I,IDACC5)=ENERGY(10+10*I,IDACC5)+SPECA(N,20,2+4*I)
        END DO
      END DO
      IDACC(5)=IDACC5
      RETURN
C****
      END SUBROUTINE DIAG4A

      SUBROUTINE get_SLP(iu_SLP)
C****
C**** THIS SUBROUTINE SAVES THE INSTANTANEOUS SEA LEVEL PRESSURES
C**** EVERY ABS(NSLP) HOURS. IF NSLP.LT.0 THE FIRST RECORD IS
C**** WRITTEN TO THE BEGINNING OF UNIT 16.
C****
      USE CONSTANT, only : grav,rgas,bygrav,bbyg,gbyrb
      USE MODEL_COM, only : im,jm,p,ptop,Itime,zatmo
      USE PBLCOM, only : tsavg
      IMPLICIT NONE
      INTEGER :: iu_SLP
      REAL*4, DIMENSION(IM,JM) :: SLP
      INTEGER :: I,J
      DO 10 J=1,JM
      DO 10 I=1,IM
   10 SLP(I,J)=(P(I,J)+PTOP)*(1.+BBYG*ZATMO(I,J)/TSAVG(I,J))**GBYRB
      CALL WRITEI(iu_SLP,Itime,SLP,IM*JM)
      RETURN
      END SUBROUTINE get_SLP

      SUBROUTINE init_DIAG(ISTART)
!@sum  init_DIAG initiallises the diagnostics
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : sday,kapa
      USE MODEL_COM, only : lm,Itime,ItimeI,Itime0,sige,sig,ptop
     *     ,pmtop,psfmpt,nfiltr,dtsrc,jhour,jdate,jmon,amon,jyear
     *     ,jhour0,jdate0,jmon0,amon0,jyear0,idacc,ioread_single
     *     ,xlabel,iowrite_single,iyear1,nday,dtsrc,dt,nmonav
     *     ,ItimeE,lrunid
      USE DAGCOM
      USE DAGPCOM, only : ple,ple_dn,plm,p1000k
      USE PARAM
      USE FILEMANAGER
      IMPLICIT NONE
      INTEGER L,K,KL,iargc,ioerr,months,years,mswitch,ldate,iu_AIC
     *     ,ISTART,jday0,jday,moff,kb
      CHARACTER FILENM*100
      LOGICAL :: QCON(NPTS), T=.TRUE. , F=.FALSE.

      call sync_param( "NAMDD", NAMDD, 4 )
      call sync_param( "IJDD", IJDD(1:8,1), 8)

      IF(ISTART.LT.0) THEN
        call getdte(Itime0,Nday,Iyear1,Jyear0,Jmon0,Jday0,Jdate0,Jhour0
     *       ,amon0)
        call getdte(Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour
     *       ,amon)
        months=1 ; years=monacc(jmon0) ; mswitch=0 ; moff=0 ; kb=jmon0
        do kl=jmon0+1,jmon0+11
          k = kl
          if (k.gt.12) k=k-12
          if (monacc(k).eq.years) then
            months=months+1
          else if (monacc(k).ne.0) then
            write(6,*) 'uneven period:',monacc
            stop 'uneven period'
          end if
          if(monacc(k).ne.monacc(kb)) mswitch = mswitch+1
          if(mswitch.eq.2) moff = moff+1
          kb = k
        end do
        if (mswitch.gt.2) then
          write(6,*) 'non-consecutive period:',monacc
          stop 'non-consecutive period'
        end if
        call aPERIOD (JMON0,JYEAR0,months,moff, years,acc_period,Ldate)
        if (iargc().gt.1) then  ! save the summed acc-file
          write(6,*) iargc(),' files are summed up'
          keyct=1 ; KEYNR=0
          XLABEL(128:132)='     '
          XLABEL(120:132)=acc_period(1:3)//' '//acc_period(4:Ldate)
          OPEN (30,FILE=acc_period(1:Ldate)//'.acc'//XLABEL(1:LRUNID),
     *         FORM='UNFORMATTED')
          call io_rsf (30,Itime,iowrite_single,ioerr)
          CLOSE (30)
        end if
        ItimeE = -1
        close (6)
        open(6,file=acc_period(1:Ldate)//'.'//XLABEL(1:LRUNID)//'.PRT',
     *       FORM='FORMATTED')
      END IF

C**** Initialize certain arrays used by more than one print routine
      DO L=1,LM
        PLE(L)=SIGE(L+1)*PSFMPT+PTOP
        PLE_DN(L)=SIGE(L)*PSFMPT+PTOP
        PLM(L)=SIG(L)*PSFMPT+PTOP
      END DO
      PLM(LM+1)=.75d0*PMTOP
      PLM(LM+2)=.35d0*PMTOP
      PLM(LM+3)=.1d0*PMTOP
      p1000k=1000.0**kapa

C**** Initialize conservation diagnostics
C**** NCON=1:23 are special cases: Angular momentum and kinetic energy
      icon_AM=1
      NOFM(:,icon_AM) = (/  1, 6, 0, 0, 0, 0, 7, 8, 9,10, 0, 0/)
      icon_KE=2
      NOFM(:,icon_KE) = (/ 12,17,18, 0, 0, 0,19,20,21,22, 0, 0/)
      NSUM_CON(1:23) = (/-1, 4, 4, 0,-1,11,11,11,11,11, 0,
     *                   -1,15,15, 0,-1,23,23,23,23,23,23, 0/)
      IA_CON(1:23) =   (/12, 6, 6,12, 6, 7, 8,10, 8, 9,12,
     *                   12, 6, 6,12, 6, 7, 8, 8,10, 8, 9,12/)
      SCALE_CON(1)              = 1d-9
      SCALE_CON((/2,3,5,6,7,9/))= 1d-2/DTSRC
      SCALE_CON((/4,11,15,23/)) = 1.
      SCALE_CON(8)              = 1d-2/(NFILTR*DTSRC)
      SCALE_CON(10)             = 2d-2/SDAY
      SCALE_CON(12)             = 25d-5
      SCALE_CON((/13,14,16/))   = 1d3 /DTSRC
      SCALE_CON((/17,18,19,21/))= 25d1/DTSRC
      SCALE_CON(20)             = 25d1/(NFILTR*DTSRC)
      SCALE_CON(22)             = 50d1/SDAY
      TITLE_CON(1:23) = (/
     *  ' INSTANTANE AM (10**9 J*S/M**2) ',
     *  ' CHANGE OF AM BY ADVECTION      ',
     *  ' CHANGE OF AM BY CORIOLIS FORCE ',
     *  ' CHANGE OF AM BY ADVEC + COR    ',
     *  ' CHANGE OF AM BY PRESSURE GRAD  ',
     *  ' CHANGE OF AM BY DYNAMICS       ',
     *  ' CHANGE OF AM BY SURFACE FRIC   ',
     *  ' CHANGE OF AM BY FILTER         ',
     *  ' CHANGE OF AM BY STRATOS DRAG   ',
     *  ' CHANGE OF AM BY DAILY RESTOR   ',
     *  ' SUM OF CHANGES (10**2 J/M**2)  ',
     *  '0INSTANTANEOUS KE (10**3 J/M**2)',
     *  ' CHANGE OF KE BY ADVECTION      ',
     *  ' CHANGE OF KE BY CORIOLIS FORCE ',
     *  ' CHANGE OF KE BY ADVEC + COR    ',
     *  ' CHANGE OF KE BY PRESSURE GRAD  ',
     *  ' CHANGE OF KE BY DYNAMICS       ',
     *  ' CHANGE OF KE BY MOIST CONVEC   ',
     *  ' CHANGE OF KE BY SURF + DRY CONV',
     *  ' CHANGE OF KE BY FILTER         ',
     *  ' CHANGE OF KE BY STRATOS DRAG   ',
     *  ' CHANGE OF KE BY DAILY RESTOR   ',
     *  ' SUM OF CHANGES (10**-3 W/M**2) '/)
      name_consrv(1:23) = (/
     *     'inst_AM   ','chg_AM_ADV','chg_AM_COR','chg_AM_ADV',
     *     'chg_AM_PRE','chg_AM_DYN','chg_AM_SUR','chg_AM_FIL',
     *     'chg_AM_STR','chg_AM_DAI','sum_chg_AM',
     *     'inst_KE   ','chg_KE_ADV','chg_KE_COR','chg_KE_ADV',
     *     'chg_KE_PRE','chg_KE_DYN','chg_KE_MOI','chg_KE_SUR',
     *     'chg_KE_FIL','chg_KE_STR','chg_KE_DAI','sum_chg_KE'/)
      units_consrv(1)    ="10**9 J*S/M**2"
      units_consrv(2:11) ="10**2 J/M**2"
      units_consrv(12)   ="10**3 J/M**2"
      units_consrv(13:23)="10**-3 W/M**2"
      lname_consrv(1:23)=TITLE_CON(1:23)
C**** To add a new conservation diagnostic:
C****    i) Add 1 to NQUANT, and increase KCON in DAGCOM.f
C****   ii) Set up a QCON, and call SET_CON to allocate array numbers,
C****       set up scales, titles, etc. The icon_XX index must be
C****       declared in DAGCOM.f for the time being
C**** QCON denotes when the conservation diags should be done
C**** 1:NPTS ==> DYN,   COND,   RAD,   PREC,   LAND,  SURF,
C****            FILTER,STRDG/OCEAN, DAILY, OCEAN1, OCEAN2,
C****  iii) Write a conserv_XYZ routine that returns the zonal average
C****       of your quantity
C****   iv) Add a line to DIAGCA that calls conserv_DIAG (declared
C****       as external)
C****    v) Note that the conserv_XYZ routine, and call to SET_CON
C****       should be in the driver module for the relevant physics

C**** Set up atmospheric component conservation diagnostics
C**** Atmospheric mass
      QCON=(/ T, F, F, F, F, F, T, F, T, F, F/)
      CALL SET_CON(QCON,"MASS    ","(KG/M**2)       ",
     *     "(10^-8 KG/S/M^2)",1d0,1d8,icon_MS)
C**** Atmospheric total potential energy
      QCON=(/ T, T, T, F, F, T, T, F, F, F, F/)
      CALL SET_CON(QCON,"TPE     ","(10**5 J/M**2)  ",
     *     "(10**-2 W/M**2) ",1d-5,1d2,icon_TPE)
C**** Atmospheric water mass
      QCON=(/ T, T, F, F, F, T, F, F, F, F, F/)
      CALL SET_CON(QCON,"WATER   ","(10**-2 KG/M**2)",
     *     "(10^-8 KG/S/M^2)",1d2,1d8,icon_WM)
C**** Atmospheric water energy
C**** This is not currently a conserved quantity, but it should be.
C**** Hence this diagnostic gives the error
      QCON=(/ T, T, F, F, F, T, F, F, F, F, F/)
      CALL SET_CON(QCON,"ENRG WAT","(J/M**2)        ",
     *     "(10^-6 J/S/M^2) ",1d0,1d6,icon_EWM)

C**** initialise longitudinal diagnostic special latitudes
      J50N=(50.+90.)*(JM-1)/180.+1.5
      J70N=(70.+90.)*(JM-1)/180.+1.5
      J5NUV = (90.+5.)*(JM-1.)/180.+2.
      J5SUV = (90.-5.)*(JM-1.)/180.+2.
      J5N   = (90.+5.)*(JM-1.)/180.+1.5
      J5S   = (90.-5.)*(JM-1.)/180.+1.5

C**** initiallise layering for spectral diagnostics
C**** add in epsilon=1d-5 to avoid roundoff mistakes
      KL=1
      DO L=1,LM
        IF (PTOP+PSFMPT*SIGE(L+1)+1d-5.lt.PSPEC(KL) .and.
     *      PTOP+PSFMPT*SIGE(L)+1d-5.gt.PSPEC(KL)) THEN
          IF (KL.eq.2) LSTR = L  ! approx. 10mb height
          KL=KL+1
        END IF
        KLAYER(L)=4*(KL-1)+1
      END DO
      IF (KL*4 .gt. NSPHER) THEN
        WRITE(6,*) "Inconsistent definitions of stratosphere:"
        WRITE(6,*) "Adjust PSPEC, ISTRAT so that KL*4 = NSPHER"
        WRITE(6,*) "ISTRAT,PSPEC,NSPHER,KL=",ISTRAT,PSPEC,NSPHER,KL
        STOP "Stratospheric definition problem for spectral diags."
      END IF

C**** Calculate the max number of geopotential heights
      do k=1,kgz
        kgz_max = k
        if (pmb(k).le.pmtop) exit
      end do

c**** Initialize acc-array names, units, idacc-indices
      call def_acc

C**** Ensure that diagnostics are reset at the beginning of the run
      IF (Itime.le.ItimeI .and. ISTART.gt.0) THEN
        CALL reset_DIAG(0)
        CALL daily_DIAG
      END IF

      RETURN
      END SUBROUTINE init_DIAG

      SUBROUTINE reset_DIAG(isum)
!@sum  reset_DIAG resets/initiallises diagnostics
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : Itime,iyear1,nday,
     *     Itime0,jhour0,jdate0,jmon0,amon0,jyear0,idacc,u
      USE DAGCOM
      USE PARAM
      IMPLICIT NONE
      INTEGER :: isum !@var isum if =1 preparation to add up acc-files
      INTEGER jd0

      IDACC(1:12)=0
      AJ=0    ; AREG=0
      APJ=0   ; AJL=0  ; ASJL=0   ; AIJ=0
      AIL=0   ; ENERGY=0 ; CONSRV=0
      SPECA=0 ; ATPE=0 ; ADIURN=0 ; WAVE=0
      AJK=0   ; AIJK=0 ; AIJL=0   ; AJLSP=0
      call reset_ODIAG(isum)

      if (isum.eq.1) return

      AIJ(:,:,IJ_TMNMX)=1000. ; IDACC(12)=1

      Itime0=Itime
      call getdte(Itime0,Nday,Iyear1,Jyear0,Jmon0,Jd0,
     *     Jdate0,Jhour0,amon0)

      CALL EPFLXI (U)  ! strat

      RETURN
      END SUBROUTINE reset_DIAG

      SUBROUTINE daily_DIAG
!@sum  daily_DIAG resets diagnostics at beginning of each day
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : im,jm,jday,fearth
      USE DAGCOM, only : tsfrez,tdiurn
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

      SUBROUTINE SET_CON(QCON,NAME_CON,INST_UNIT,SUM_UNIT,INST_SC
     *     ,CHNG_SC,ICON)
!@sum  SET_CON assigns conservation diagnostic array indices
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : sday
      USE MODEL_COM, only : dtsrc,nfiltr
      USE DAGCOM, only : kcon,nquant,npts,title_con,scale_con,nsum_con
     *     ,nofm,ia_con,kcmx,ia_d5d,ia_d5s,ia_filt,ia_12hr,name_consrv
     *     ,lname_consrv,units_consrv
      IMPLICIT NONE
!@var QCON logical variable sets where conservation diags are saved
      LOGICAL, INTENT(IN),DIMENSION(NPTS) :: QCON
!@var INST_SC scale for instantaneous value
      REAL*8, INTENT(IN) :: INST_SC
!@var CHNG_SC scale for changes
      REAL*8, INTENT(IN) :: CHNG_SC
!@var NAME_CON name of conservation quantity
      CHARACTER*8, INTENT(IN) :: NAME_CON
!@var sname name of conservation quantity (no spaces)
      CHARACTER*8 :: sname
!@var INST_UNIT string for unit for instant. values
      CHARACTER*16, INTENT(IN) :: INST_UNIT
!@var SUM_UNIT string for unit for summed changes
      CHARACTER*16, INTENT(IN) :: SUM_UNIT
!@var ICON index for the conserved quantity
      INTEGER, INTENT(OUT) :: ICON
!@var CONPT titles for each point at which conservation diags. are done
      CHARACTER*10, DIMENSION(NPTS) :: CONPT = (/
     *     "DYNAMICS  ","CONDENSATN","RADIATION ","PRECIPITAT",
     *     "LAND SURFC","SURFACE   ","FILTER    ","OCEAN     ",
     *     "DAILY     ","OCN DYNAM ","OCEAN PHYS"/)
      INTEGER NI,NM,NS,N,k
      INTEGER, SAVE :: NQ = 2   ! first 2 special cases AM + KE

      NQ=NQ+1
      IF (NQ.gt.NQUANT) THEN
        WRITE(6,*) "Number of conserved quantities larger than NQUANT"
     *       ,NQUANT,NQ
        STOP "Change NQUANT in diagnostic common block"
      END IF
C**** remove spaces in NAME_CON for netcdf names
      sname=TRIM(NAME_CON)
      do k=1,len_trim(NAME_CON)
        if (sname(k:k).eq." ") sname(k:k)="_"
      end do
C****
      NI=KCMX+1
      NOFM(1,NQ) = NI
      TITLE_CON(NI) = "0INSTANT "//TRIM(NAME_CON)//" "//TRIM(INST_UNIT)
      SCALE_CON(NI) = INST_SC
      name_consrv(NI) ="inst_"//sname
      lname_consrv(NI) = "INSTANT "//TRIM(NAME_CON)
      units_consrv(NI) = INST_UNIT
      IA_CON(NI) = 12
      NM=NI
      DO N=1,NPTS
        IF (QCON(N)) THEN
          NM = NM + 1
          NOFM(N+1,NQ) = NM
          TITLE_CON(NM) = " CHANGE OF "//TRIM(NAME_CON)//" BY "//
     *         CONPT(N)
          name_consrv(NM) ="chg_"//trim(sname)//"_"//TRIM(CONPT(N)(1:3))
          lname_consrv(NM) = TITLE_CON(NM)
          units_consrv(NM) = SUM_UNIT
          SELECT CASE (N)
          CASE (1)
            SCALE_CON(NM) = CHNG_SC/DTSRC
            IA_CON(NM) = ia_d5d
          CASE (2,3,4,5,6,8,10,11)
            SCALE_CON(NM) = CHNG_SC/DTSRC
            IA_CON(NM) = ia_d5s
          CASE (7)
            SCALE_CON(NM) = CHNG_SC/(NFILTR*DTSRC)
            IA_CON(NM) = ia_filt
          CASE (9)
            SCALE_CON(NM) = CHNG_SC*2./SDAY
            IA_CON(NM) = ia_12hr
          END SELECT
        ELSE
          NOFM(N+1,NQ) = 0
        END IF
      END DO
      NS=NM+1
      TITLE_CON(NS) = " SUM OF CHANGES "//TRIM(SUM_UNIT)
      name_consrv(NS) ="sum_chg_"//trim(sname)
      lname_consrv(NS) = " SUM OF CHANGES OF "//TRIM(NAME_CON)
      units_consrv(NS) = SUM_UNIT
      SCALE_CON(NS) = 1.
      IA_CON(NS) = 12
      NSUM_CON(NI) = -1
      NSUM_CON(NI+1:NS-1) = NS
      NSUM_CON(NS) = 0
      KCMX=NS
      ICON=NQ
      IF (KCMX.gt.KCON) THEN
        WRITE(6,*) "KCON not large enough for extra conserv diags",
     *       KCON,NI,NM,NQ,NS,NAME_CON
        STOP "Change KCON in diagnostic common block"
      END IF
      RETURN
C****
      END SUBROUTINE set_con
