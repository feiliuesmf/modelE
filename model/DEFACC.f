      subroutine def_acc
c-----------------------------------------------------------------------
c define acc names, units, etc
c-----------------------------------------------------------------------
      use DAGCOM
      implicit none
      integer :: k
      call iparm_defs
      call dparm_defs
      call tsf_defs
      call j_defs
      name_reg=name_j
      do k=1,kaj ! to avoid naming conflicts, put a prefix
         name_j(k) = 'J_'//trim(name_j(k))
         name_reg(k) = 'reg_'//trim(name_reg(k))
      enddo
      call pj_defs
      call jl_defs
      call sjl_defs
      call ij_defs
      call il_defs
      call ijg_defs
      call consrv_defs
      call wave_defs
      call jk_defs
      call ijk_defs
      call ijl_defs
      call jlsp_defs
      return
      end subroutine def_acc

      subroutine iparm_defs
      use MODEL_COM
      use DAGCOM
      implicit none
      integer :: k
c
      niparm=0
      do k=1,niparm_max
         iparm_name(k) = 'nowrite'
      enddo
c contents of iparmb common block
      call set_iparm('IM0',IM0)
      call set_iparm('JM0',JM0)
      call set_iparm('LM0',LM0)
      call set_iparm('LS1',LS1)
      call set_iparm('KACC0',KACC0)
      call set_iparm('KTACC0',KTACC0)
      call set_iparm('Itime',Itime)
      call set_iparm('ItimeI',ItimeI)
      call set_iparm('ItimeE',ItimeE)
      call set_iparm('Itime0',Itime0)
      call set_iparm('KOCEAN',KOCEAN)
      call set_iparm('KDISK',KDISK)
      call set_iparm('KEYCT',KEYCT)
      call set_iparm('KCOPY',KCOPY)
      call set_iparm('IRAND',IRAND)
      call set_iparm('MFILTR',MFILTR)
      call set_iparm('Ndisk',Ndisk)
      call set_iparm('Kvflxo',Kvflxo)
      call set_iparm('Nslp',Nslp)
      call set_iparm('NIdyn',NIdyn)
      call set_iparm('NRAD',NRAD)
      call set_iparm('NIsurf',NIsurf)
      call set_iparm('NFILTR',NFILTR)
      call set_iparm('NDAY',NDAY)
      call set_iparm('NDAA',NDAA)
      call set_iparm('NDA5D',NDA5D)
      call set_iparm('NDA5K',NDA5K)
      call set_iparm('NDA5S',NDA5S)
      call set_iparm('NDA4',NDA4)
      call set_iparm('NDASF',NDASF)
      call set_iparm('MLAST',MLAST)
      call set_iparm('MDYN',MDYN)
      call set_iparm('MCNDS',MCNDS)
      call set_iparm('MRAD',MRAD)
      call set_iparm('MSURF',MSURF)
      call set_iparm('MDIAG',MDIAG)
      call set_iparm('MELSE',MELSE)
      call set_iparm('MODRD',MODRD)
      call set_iparm('MODD5K',MODD5K)
      call set_iparm('MODD5S',MODD5S)
      call set_iparm('IYEAR0',IYEAR0)
      call set_iparm('JYEAR',JYEAR)
      call set_iparm('JYEAR0',JYEAR0)
      call set_iparm('JMON',JMON)
      call set_iparm('JMON0',JMON0)
      call set_iparm('JDATE',JDATE)
      call set_iparm('JDATE0',JDATE0)
      call set_iparm('JHOUR',JHOUR)
      call set_iparm('JHOUR0',JHOUR0)
      call set_iparm('JDAY',JDAY)
      call set_iparm('NSSW',NSSW)
      call set_iparm('NSTEP',NSTEP)
      call set_iparm('MRCH',MRCH)
      call set_iparm('NIPRNT',NIPRNT)
      call set_iparm('NMONAV',NMONAV)
c
      return
      end subroutine iparm_defs

      subroutine set_iparm(parm_str,ival)
      use DAGCOM
      implicit none
      character(len=20) :: parm_str
      integer :: ival
      if(niparm.eq.niparm_max) stop 'set_iparm: too many iparms'
      niparm = niparm + 1
      iparm_name(niparm) = parm_str
      iparm(niparm) = ival
      return
      end subroutine set_iparm

      subroutine dparm_defs
      use CONSTANT
      use MODEL_COM
      use DAGCOM
      use RADNCB, only : co2
      use CLD01, only : u00wtr,u00ice
      implicit none
      integer :: k
c
      ndparm=0
      do k=1,ndparm_max
         dparm_name(k) = 'nowrite'
      enddo
c contents of rparmb common block
      call set_dparm('DTsrc',DTsrc)
      call set_dparm('DT',DT)
      call set_dparm('PTOP',PTOP)
      call set_dparm('PSF',PSF)
      call set_dparm('PSFMPT',PSFMPT)
      call set_dparm('PSTRAT',PSTRAT)
      call set_dparm('PSDRAG',PSDRAG)
      call set_dparm('SKIPSE',SKIPSE)
c from constants module
      call set_dparm('SDAY',SDAY)
      call set_dparm('GRAV',GRAV)
      call set_dparm('OMEGA',OMEGA)
      call set_dparm('RGAS',RGAS)
      call set_dparm('KAPA',KAPA)
      call set_dparm('TF',TF)
c cloud parameters
      call set_dparm('U00wtr',U00wtr)
      call set_dparm('U00ice',U00ice)
c radiation parameters
      call set_dparm('CO2',CO2)
c
      return
      end subroutine dparm_defs

      subroutine set_dparm(parm_str,dval)
      use DAGCOM
      implicit none
      character(len=20) :: parm_str
      double precision :: dval
      if(ndparm.eq.ndparm_max) stop 'set_dparm: too many dparms'
      ndparm = ndparm + 1
      dparm_name(ndparm) = parm_str
      dparm(ndparm) = dval
      return
      end subroutine set_dparm

      subroutine j_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c  AJ01
      k=k+1
      J_SRINCP0= k ! SRINCP0 (W/M**2)                              2 RD
      name_j(k) = 'SRINCP0'
      lname_j(k) = 'INC SW'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_rad
c  AJ02
      k=k+1
      J_SRNFP0=  k ! SRNFP0 (W/M**2)                               2 RD
      name_j(k) = 'SRNFP0'
      lname_j(k) = 'PLANETARY ALBEDO'
      units_j(k) = '1'
      ia_j(k) = ia_rad
c  AJ03
      k=k+1
      J_SRNFP1=  k ! SRNFP1 (W/M**2)                               2 RD
      name_j(k) = 'SRNFP1'
      lname_j(k) = 'SRNFP1'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_rad
c  AJ04
      k=k+1
      J_SRABS =  k ! SRABSATM=AJ(SRNFP0)-AJ(SRNFG) (W/M**2)        2 D1
      name_j(k) = 'SRABS'
      lname_j(k) = 'SRABSATM'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_rad
c  AJ05
      k=k+1
      J_SRINCG=  k ! SRINCG (W/M**2)                               2 RD
      name_j(k) = 'SRINCG'
      lname_j(k) = 'SW INC ON ZO'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_rad
c  AJ06
      k=k+1
      J_SRNFG =  k ! SRNFG (W/M**2)                                2 RD
      name_j(k) = 'SRNFG'
      lname_j(k) = 'SW ABS AT ZO'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_rad
c  AJ07
      k=k+1
      J_TRNFP0=  k ! TRNFP0=AJ(ALBVIS)+A2BYA1*AJ(TRHDT)/DTS(W/M**2)2 D1
      name_j(k) = 'TRNFP0'
      lname_j(k) = 'TRNFP0=AJ(74)+A2BYA1*AJ(9)/DTSRCE'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_rad
c  AJ08
      k=k+1
      J_TRNFP1=  k ! TRNFP1=AJ(ALBNIR)+A2BYA1*AJ(TRHDT)/DTS(W/M**2)2 D1
      name_j(k) = 'TRNFP1'
      lname_j(k) = 'TRNFP1=AJ(75)+A2BYA1*AJ(9)/DTSRCE'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_rad
c  AJ09
      k=k+1
      J_TRHDT =  k ! TRHDT (J/M**2)                                1 SF
      name_j(k) = 'TRHDT'
      lname_j(k) = 'NET LW AT ZO (positive downward)'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_src
c  AJ10
      k=k+1
      J_RNFP0 = k ! RNFP0=AJ(SRNFP0)+AJ(TRNFP0) (W/M**2)          2 D1
      name_j(k) = 'RNFP0'
      lname_j(k) = 'RNFP0=AJ(2)+AJ(7)'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_rad
c  AJ11
      k=k+1
      J_RNFP1 = k ! RNFP1=AJ(SRNFP1)+AJ(TRNFP1) (W/M**2)          2 D1
      name_j(k) = 'RNFP1'
      lname_j(k) = 'RNFP1=AJ(3)+AJ(8)'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_rad
c  AJ12
      k=k+1
      J_RHDT  = k ! RHDT=A1BYA2*AJ(SRNFG)*DTS+AJ(TRHDT)(J/M^2)    1 D1
      name_j(k) = 'RHDT'
      lname_j(k) = 'RHDT=A1BYA2*AJ(6)*DTSRCE+AJ(9)'
      units_j(k) = 'J/m^2'
      ia_j(k) = ia_src
c  AJ13
      k=k+1
      J_SHDT  = k ! SHEATDT (J/M**2)                              1 SF
      name_j(k) = 'SHDT'
      lname_j(k) = 'SENS HEAT FLUX'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_src
c  AJ14
      k=k+1
      J_EVHDT = k ! EVHDT (J/M**2)                                1 SF
      name_j(k) = 'EVHDT'
      lname_j(k) = 'LATENT HEAT FLUX'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_src
c  AJ15
      k=k+1
      J_F2DT  = k ! F2DT (J/M**2)                                 1 GD
      name_j(k) = 'F2DT'
      lname_j(k) = 'F2DT'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_src
c  AJ16
      k=k+1
      J_HZ1   = k ! HEATZ1=AJ(EDIFS)+AJ(F1DT)                     1 D1
      name_j(k) = 'HZ1'
      lname_j(k) = 'HEATZ1=AJ(41)+AJ(42)'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_src
c  AJ17
      k=k+1
      J_TG2   = k ! TG2 (K-TF)                                    1 GD
      name_j(k) = 'TG2'
      lname_j(k) = 'TG2'
      units_j(k) = 'degC'
      ia_j(k) = ia_src
c  AJ18
      k=k+1
      J_TG1   = k ! TG1 (K-TF)                                    1 GD
      name_j(k) = 'TG1'
      lname_j(k) = 'TG1'
      units_j(k) = 'degC'
      ia_j(k) = ia_src
c  AJ19
      k=k+1
      J_EVAP  = k ! EVAP (KG/M**2)                                1 GD
      name_j(k) = 'EVAP'
      lname_j(k) = 'EVAP'
      units_j(k) = 'mm/day'
      ia_j(k) = ia_src
c  AJ20
      k=k+1
      J_PRCP  = k ! PRCP=AJ(PRCPSS)+AJ(PRCPMC) (100 PA)           1 D1
      name_j(k) = 'PRCP'
      lname_j(k) = 'PREC'
      units_j(k) = 'mm/day'
      ia_j(k) = ia_src
c  AJ21
      k=k+1
      J_TX    = k ! TX (K-TF)  (INTEGRAL OVER ATMOSPHERE OF)      4 DA
      name_j(k) = 'TX'
      lname_j(k) = 'TX'
      units_j(k) = 'degC'
      ia_j(k) = ia_dga
c  AJ22
      k=k+1
      J_TX1   = k ! TX1 (K-TF)                                    4 DA
      name_j(k) = 'TX1'
      lname_j(k) = 'TX1'
      units_j(k) = 'degC'
      ia_j(k) = ia_dga
c  AJ23
      k=k+1
      J_TSRF  = k ! TS (K-TF)                                     3 SF
      name_j(k) = 'TSRF'
      lname_j(k) = 'TS'
      units_j(k) = 'degC'
      ia_j(k) = ia_srf
c  AJ24
      k=k+1
      J_DTSGST= k ! DTH/DPHI  (STRATOSPHERE)                      4 DA
      name_j(k) = 'DTSGST'
      lname_j(k) = 'DTH/DPHI  (STRATOSPHERE)'
      units_j(k) = 'unknown'
      ia_j(k) = ia_dga
c  AJ25
      k=k+1
      J_DTDGTR= k ! DTH/DPHI  (TROPOSPHERE)                       4 DA
      name_j(k) = 'DTDGTR'
      lname_j(k) = 'DTH/DPHI  (TROPOSPHERE)'
      units_j(k) = 'unknown'
      ia_j(k) = ia_dga
c  AJ26
      k=k+1
      J_RICST = k ! .0625*DTH*DLNP/(DU*DU+DV*DV)  (STRATOSPHERE)  4 DA
      name_j(k) = 'RICST'
      lname_j(k) = 'RICH NUM (STRATOSPHERE)'
      units_j(k) = '1'
      ia_j(k) = ia_dga
c  AJ27
      k=k+1
      J_RICTR = k ! .0625*DTH*DLNP/(DU*DU+DV*DV)  (TROPOSPHERE)   4 DA
      name_j(k) = 'RICTR'
      lname_j(k) = 'RICH NUM (TROPOSPHERE)'
      units_j(k) = '1'
      ia_j(k) = ia_dga
c  AJ28
      k=k+1
      J_ROSST = k ! 4*UMAX/(DX*SINJ)  (STRATOSPHERE)              4 DA
      name_j(k) = 'ROSST'
      lname_j(k) = 'ROSS NUM (STRATOSPHERE)'
      units_j(k) = '1'
      ia_j(k) = ia_dga
c  AJ29
      k=k+1
      J_ROSTR = k ! 4*UMAX/(DX*SINJ)  (TROPOSPHERE)               4 DA
      name_j(k) = 'ROSTR'
      lname_j(k) = 'ROSS NUM (TROPOSPHERE)'
      units_j(k) = '1'
      ia_j(k) = ia_dga
c  AJ30
      k=k+1
      J_RSI   = k ! RSI (1)                                       1 GD
      name_j(k) = 'RSI'
      lname_j(k) = 'POICE'
      units_j(k) = '1'
      ia_j(k) = ia_src
c  AJ31
      k=k+1
      J_RSNOW = k ! PSNOW (1)                                     4 DA
      name_j(k) = 'RSNOW'
      lname_j(k) = 'PSNOW'
      units_j(k) = '1'
      ia_j(k) = ia_dga
c  AJ32
      k=k+1
      J_SWCOR = k ! SW CORRECTION   (obsolete)                    2 RD
      name_j(k) = 'SWCOR'
      lname_j(k) = 'SW CORRECTION'
      units_j(k) = 'unknown'
      ia_j(k) = ia_rad
c  AJ33
      k=k+1
      J_OHT   = k ! OCEAN TRANSPORT                               1 GD
      name_j(k) = 'OHT'
      lname_j(k) = 'OCEAN TRANSPORT'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_src
c  AJ34
      k=k+1
      J_OMLT  = k ! OCEAN TEMPERATURE AT MAX. MIXED LAYER DEPTH   1 GD
      name_j(k) = 'OMLT'
      lname_j(k) = 'OCEAN TEMP AT MAX. MIXED LAYER DEPTH'
      units_j(k) = 'unknown'
      ia_j(k) = ia_src
c  AJ35
      k=k+1
      J_DTDJS = k ! T(J+1)-T(J-1)  (SUM OVER STRATOSPHERE OF)     4 DA
      name_j(k) = 'DTDJS'
      lname_j(k) = 'T(J+1)-T(J-1)  (SUM OVER STRATOSPHERE OF)'
      units_j(k) = 'K'
      ia_j(k) = ia_dga
c  AJ36
      k=k+1
      J_DTDJT = k ! T(J+1)-T(J-1)  (SUM OVER TROPOSPHERE OF)      4 DA
      name_j(k) = 'DTDJT'
      lname_j(k) = 'T(J+1)-T(J-1)  (SUM OVER TROPOSPHERE OF)'
      units_j(k) = 'K'
      ia_j(k) = ia_dga
c  AJ37
      k=k+1
      J_LSTR  = k ! SQRT(DTH/DLNP)/SINJ  (STRATOSPHERE)           4 DA
      name_j(k) = 'LSTR'
      lname_j(k) = 'SQRT(DTH/DLNP)/SINJ  (STRATOSPHERE)'
      units_j(k) = 'unknown'
      ia_j(k) = ia_dga
c  AJ38
      k=k+1
      J_LTRO  = k ! SQRT(DTH/DLNP)/SINJ  (TROPOSPHERE)            4 DA
      name_j(k) = 'LTRO'
      lname_j(k) = 'SQRT(DTH/DLNP)/SINJ'
      units_j(k) = 'unknown'
      ia_j(k) = ia_dga
c  AJ39
      k=k+1
      J_EPRCP = k ! ENERGP (J/M**2)                               1 CN
      name_j(k) = 'EPRCP'
      lname_j(k) = 'ENERGP'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_src
c  AJ40
      k=k+1
      J_ERUN1 = k ! ERUN1 (J/M**2)                                1 GP
      name_j(k) = 'ERUN1'
      lname_j(k) = 'ERUN1'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_src
c  AJ41
      k=k+1
      J_EDIFS = k ! EDIFS (J/M**2)                                1 GP
      name_j(k) = 'EDIFS'
      lname_j(k) = 'EDIFS'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_src
c  AJ42
      k=k+1
      J_F1DT  = k ! F1DT (J/M**2)                                 1 GD
      name_j(k) = 'F1DT'
      lname_j(k) = 'F1DT'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_src
c  AJ43
      k=k+1
      J_ERUN2 = k ! ERUN2 (J/M**2)                                1 GP
      name_j(k) = 'ERUN2'
      lname_j(k) = 'ERUN2'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_src
c  AJ44
      k=k+1
      J_HZ0   = k ! HEATZ0=RHDT+SHDT+EVHDT+EPRCP-EDIFS (J/M**2)   1 D1
      name_j(k) = 'HZ0'
      lname_j(k) = 'NET HT Z0=AJ(12)+AJ(13)+AJ(14)+AJ(39)-AJ(40)'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_src
c  AJ45
      k=k+1
      J_DIFS  = k ! DIFS (KG/M**2)                                1 GP
      name_j(k) = 'DIFS'
      lname_j(k) = 'DIFS'
      units_j(k) = 'mm/day'
      ia_j(k) = ia_src
c  AJ46
      k=k+1
      J_IMELT = k ! AIFO ; BRUN2 ; CRUN2+CIFI                     1 GP
      name_j(k) = 'IMELT'
      lname_j(k) = 'AIFO; BRUN2; CRUN2+CIFI'
      units_j(k) = 'mm/day'
      ia_j(k) = ia_src
c  AJ47
      k=k+1
      J_RUN2  = k ! RUN2 (KG/M**2)                                1 GP
      name_j(k) = 'RUN2'
      lname_j(k) = 'RUN2'
      units_j(k) = 'mm/day'
      ia_j(k) = ia_src
c  AJ48
      k=k+1
      J_DWTR2 = k ! DWTR2=AJ(DIFS)-AJ(RUN2) (KG/M**2)             1 D1
      name_j(k) = 'DWTR2'
      lname_j(k) = 'HEAT RUNOFF THROUGH THE MIXED LAYER DEPTH'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_src
c  AJ49
      k=k+1
      J_WTR1  = k ! WTR1 (KG/M**2)                                1 GD
      name_j(k) = 'WTR1'
      lname_j(k) = 'WTR1'
      units_j(k) = 'kg/m^2'
      ia_j(k) = ia_src
c  AJ50
      k=k+1
      J_ACE1  = k ! ACE1 (KG/M**2)                                1 GD
      name_j(k) = 'ACE1'
      lname_j(k) = 'ACE1'
      units_j(k) = 'kg/m^2'
      ia_j(k) = ia_src
c  AJ51
      k=k+1
      J_WTR2  = k ! WTR2 (KG/M**2)                                1 GD
      name_j(k) = 'WTR2'
      lname_j(k) = 'WTR2'
      units_j(k) = 'kg/m^2'
      ia_j(k) = ia_src
c  AJ52
      k=k+1
      J_ACE2  = k ! ACE2 (KG/M**2)                                1 GD
      name_j(k) = 'ACE2'
      lname_j(k) = 'ACE2'
      units_j(k) = 'kg/m^2'
      ia_j(k) = ia_src
c  AJ53
      k=k+1
      J_SNOW  = k ! SNOW (KG/M**2)                                1 GD
      name_j(k) = 'SNOW'
      lname_j(k) = 'SNOW'
      units_j(k) = 'kg/m^2'
      ia_j(k) = ia_src
c  AJ54
      k=k+1
      J_RUN1  = k ! RUN1 (KG/M**2)                                1 GP
      name_j(k) = 'RUN1'
      lname_j(k) = 'RUN1'
      units_j(k) = 'mm/day'
      ia_j(k) = ia_src
c  AJ55
      k=k+1
      J_BRTEMP= k ! BTEMPW-TF                                     2 RD
      name_j(k) = 'BRTEMP'
      lname_j(k) = 'BTEMPW'
      units_j(k) = 'degC'
      ia_j(k) = ia_rad
c  AJ56
      k=k+1
      J_HZ2   = k ! HEATZ2=AJ(F2DT)+AJ(ERUN2) (J/M**2)            1 D1
      name_j(k) = 'HZ2'
      lname_j(k) = 'HEATZ2=AJ(15)+AJ(43)'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_src
c  AJ57
      k=k+1
      J_PCLDSS= k ! PCLDSS (1)  (COMPOSITE OVER ATMOSPHERE)       2 RD
      name_j(k) = 'PCLDSS'
      lname_j(k) = 'PCLDSS (COMPOSITE OVER ATMOSPHERE)'
      units_j(k) = '1'
      ia_j(k) = ia_rad
c  AJ58
      k=k+1
      J_PCLDMC= k ! PCLDMC (1)  (COMPOSITE OVER ATMOSPHERE)       2 RD
      name_j(k) = 'PCLDMC'
      lname_j(k) = 'PCLDMC (COMPOSITE OVER ATMOSPHERE)'
      units_j(k) = '1'
      ia_j(k) = ia_rad
c  AJ59
      k=k+1
      J_PCLD  = k ! PCLD (1)  (COMPOSITE OVER ATMOSPHERE)         2 RD
      name_j(k) = 'PCLD'
      lname_j(k) = 'PCLD (COMPOSITE OVER ATMOSPHERE)'
      units_j(k) = '1'
      ia_j(k) = ia_rad
c  AJ60
      k=k+1
      J_CTOPP = k ! CLDTOPMC=AJ(CDLDEP)/AJ(PCLDMC) (100 PA)       0 D1
      name_j(k) = 'CTOPP'
      lname_j(k) = 'CLDTOPMC=AJ(80)/AJ(58)'
      units_j(k) = '100 PA'
      ia_j(k) = ia_src
c  AJ61
      k=k+1
      J_PRCPSS= k ! PRCPSS (100 PA)                               1 CN
      name_j(k) = 'PRCPSS'
      lname_j(k) = 'PRCPSS'
      units_j(k) = 'mm/day'
      ia_j(k) = ia_src
c  AJ62
      k=k+1
      J_PRCPMC= k ! PRCPMC (100 PA)                               1 CN
      name_j(k) = 'PRCPMC'
      lname_j(k) = 'PRCPMC'
      units_j(k) = '100 PA'
      ia_j(k) = ia_src
c  AJ63
      k=k+1
      J_QP    = k ! Q*P (100 PA)  (INTEGRAL OVER ATMOSPHERE OF)   4 DA
      name_j(k) = 'QP'
      lname_j(k) = 'Q*P (INTEGRAL OVER ATMOSPHERE OF)'
      units_j(k) = 'kg/m^2'
      ia_j(k) = ia_dga
c  AJ64
      k=k+1
      J_GAM   = k ! GAM  (K/M)  (*SIG(TROPOSPHERE)/GRAV)          4 DA
      name_j(k) = 'GAM'
      lname_j(k) = 'GAM  (*SIG(TROPOSPHERE)/GRAV)'
      units_j(k) = 'K/m'
      ia_j(k) = ia_dga
c  AJ65
      k=k+1
      J_GAMM  = k ! GAMM  (K-S**2/M**2)  (SIG(TROPOSPHERE)/GAMD)  4 DA
      name_j(k) = 'GAMM'
      lname_j(k) = 'GAMM (SIG(TROPOSPHERE)/GAMD)'
      units_j(k) = 'K s^2/m^2'
      ia_j(k) = ia_dga
c  AJ66
      k=k+1
      J_GAMC  = k ! GAMC  (K/M)                                   4 DA
      name_j(k) = 'GAMC'
      lname_j(k) = 'GAMC'
      units_j(k) = 'K/m'
      ia_j(k) = ia_dga
c  AJ67
      k=k+1
      J_TRINCG= k ! TRINCG (W/M**2)                               2 RD
      name_j(k) = 'TRINCG'
      lname_j(k) = 'LW INC ON ZO'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_rad
c  AJ68
      k=k+1
      J_FTHERM= k ! ENERGY DIFFUSION INTO THERMOCLINE (W/M**2) .5*9 MN
      name_j(k) = 'FTHERM'
      lname_j(k) = 'ENERGY DIFFUSION INTO THERMOCLINE'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_12hr
c  AJ69
      k=k+1
      J_TYPE  = k ! PTYPE                                         1 GD
      name_j(k) = 'PTYPE'
      lname_j(k) = 'sum of surface type fraction'
      units_j(k) = 'unknown'
      ia_j(k) = ia_src
c  AJ70
      k=k+1
      J_HSURF = k ! TRNFP0-TRNFG (W/M**2)                         2 RD
      name_j(k) = 'HSURF'
      lname_j(k) = 'TRNFP0-TRNFG'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_rad
c  AJ71
      k=k+1
      J_HATM  = k ! TRNFP1-TRNFG (W/M**2)                         2 RD
      name_j(k) = 'HATM'
      lname_j(k) = 'TRNFP1-TRNFG'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_rad
c  AJ72
      k=k+1
C**** Note: next eight diagnostics must remain in order
      J_PLAVIS= k ! PLAVIS*S0*COSZ (W/M**2)                       2 RD
      name_j(k) = 'PLAVIS'
      lname_j(k) = 'PLAVIS*S0*COSZ'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_rad
c  AJ73
      k=k+1
      J_PLANIR= k ! PLANIR*S0*COSZ (W/M**2)                       2 RD
      name_j(k) = 'PLANIR'
      lname_j(k) = 'PLANIR*S0*COSZ'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_rad
c  AJ74
      k=k+1
      J_ALBVIS= k ! ALBVIS*S0*COSZ (W/M**2)                       2 RD
      name_j(k) = 'ALBVIS'
      lname_j(k) = 'ALBVIS*S0*COSZ'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_rad
c  AJ75
      k=k+1
      J_ALBNIR= k ! ALBNIR*S0*COSZ (W/M**2)                       2 RD
      name_j(k) = 'ALBNIR'
      lname_j(k) = 'ALBNIR*S0*COSZ'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_rad
c  AJ76
      k=k+1
      J_SRRVIS= k ! SRRVIS*S0*COSZ (W/M**2)                       2 RD
      name_j(k) = 'SRRVIS'
      lname_j(k) = 'SRRVIS*S0*COSZ'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_rad
c  AJ77
      k=k+1
      J_SRRNIR= k ! SRRNIR*S0*COSZ (W/M**2)                       2 RD
      name_j(k) = 'SRRNIR'
      lname_j(k) = 'SRRNIR*S0*COSZ'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_rad
c  AJ78
      k=k+1
      J_SRAVIS= k ! SRAVIS*S0*COSZ (W/M**2)                       2 RD
      name_j(k) = 'SRAVIS'
      lname_j(k) = 'SRAVIS*S0*COSZ'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_rad
c  AJ79
      k=k+1
      J_SRANIR= k ! SRANIR*S0*COSZ (W/M**2)                       2 RD
      name_j(k) = 'SRANIR'
      lname_j(k) = 'SRANIR*S0*COSZ'
      units_j(k) = 'W/m^2'
      ia_j(k) = ia_rad
c  AJ80
      k=k+1
      J_CDLDEP= k ! PBOTMC-PTOPMC (100 PA)                        2 RD
      name_j(k) = 'CDLDEP'
      lname_j(k) = 'PBOTMC-PTOPMC'
      units_j(k) = '100 PA'
      ia_j(k) = ia_rad
c
      k=k+1
      name_j(k) = 'AJ81'
      lname_j(k) = 'unknown'
      units_j(k) = 'unknown'
      ia_j(k) = 0
c
      k=k+1
      name_j(k) = 'AJ82'
      lname_j(k) = 'unknown'
      units_j(k) = 'unknown'
      ia_j(k) = 0
c
      k=k+1
      name_j(k) = 'AJ83'
      lname_j(k) = 'unknown'
      units_j(k) = 'unknown'
      ia_j(k) = 0
c
      k=k+1
      name_j(k) = 'AJ84'
      lname_j(k) = 'unknown'
      units_j(k) = 'unknown'
      ia_j(k) = 0
c
      k=k+1
      name_j(k) = 'AJ85'
      lname_j(k) = 'unknown'
      units_j(k) = 'unknown'
      ia_j(k) = 0
c
      k=k+1
      name_j(k) = 'AJ86'
      lname_j(k) = 'unknown'
      units_j(k) = 'unknown'
      ia_j(k) = 0
c
      k=k+1
      name_j(k) = 'AJ87'
      lname_j(k) = 'unknown'
      units_j(k) = 'unknown'
      ia_j(k) = 0
c
      k=k+1
      name_j(k) = 'AJ88'
      lname_j(k) = 'unknown'
      units_j(k) = 'unknown'
      ia_j(k) = 0
c
      k=k+1
      name_j(k) = 'AJ89'
      lname_j(k) = 'unknown'
      units_j(k) = 'unknown'
      ia_j(k) = 0
c
      k=k+1
      name_j(k) = 'AJ90'
      lname_j(k) = 'unknown'
      units_j(k) = 'unknown'
      ia_j(k) = 0
c
      k=k+1
      name_j(k) = 'AJ91'
      lname_j(k) = 'unknown'
      units_j(k) = 'unknown'
      ia_j(k) = 0
c
      k=k+1
      name_j(k) = 'AJ92'
      lname_j(k) = 'unknown'
      units_j(k) = 'unknown'
      ia_j(k) = 0
c
      k=k+1
      name_j(k) = 'AJ93'
      lname_j(k) = 'unknown'
      units_j(k) = 'unknown'
      ia_j(k) = 0
c
      k=k+1
      name_j(k) = 'AJ94'
      lname_j(k) = 'unknown'
      units_j(k) = 'unknown'
      ia_j(k) = 0
c
      return
      end subroutine j_defs

      subroutine ij_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
C**** AIJ diagnostic names:
C**** NAME     NO.    DESCRIPTION   (SCALE)*IDACC  LOCATION
C**** ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
      k=k+1 !  'AIJ001'
      IJ_RSOI = k ! POICE (1)            1 GD
      lname_ij(k) = 'ACCUMULATED OCEAN ICE FRACTION'
      units_ij(k) = '1'
      name_ij(k) = 'RSOI'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ002'
      IJ_RSNW = k ! PSNOW (1)            4 DA
      lname_ij(k) = 'ACCUMULATED SNOW FRACTION'
      units_ij(k) = '1'
      name_ij(k) = 'RSNW'
      ia_ij(k) = ia_dga
c
      k=k+1 !  'AIJ003'
      IJ_SNOW = k ! SNOW (KG/M**2)       4 DA
      lname_ij(k) = 'ACCUMULATED SNOW MASS'
      units_ij(k) = 'kg/m^2'
      name_ij(k) = 'SNOW'
      ia_ij(k) = ia_dga
c
      k=k+1 !  'AIJ004'
      IJ_SHDT = k ! SHDT (J/M**2)        1 SF
      lname_ij(k) = 'INTEGRATED SENSIBLE HEAT FLUX'
      units_ij(k) = 'J/m^2'
      name_ij(k) = 'SHDT'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ005'
      IJ_PREC = k ! PREC (KG/M**2)       1 CN
      lname_ij(k) = 'ACCUMULATED PRECIPITATION'
      units_ij(k) = 'kg/m^2'
      name_ij(k) = 'PREC'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ006'
      IJ_EVAP = k ! EVAP (KG/M**2)       1 SF
      lname_ij(k) = 'ACCUMULATED EVAPORATION'
      units_ij(k) = 'kg/m^2'
      name_ij(k) = 'EVAP'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ007'
      IJ_BETA = k ! BETA (1)             1 GD
      lname_ij(k) = 'ACCUMULATED EVAPORATION EFFICIENCY'
      units_ij(k) = '1'
      name_ij(k) = 'BETA'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ008'
      IJ_PRES = k ! PIJ (100 PA)  (NO PRINTOUT)  4 DA
      lname_ij(k) = 'ACCUMULATED SURFACE PRESSURE'
      units_ij(k) = '100 PA'
      name_ij(k) = 'PRES'
      ia_ij(k) = ia_dga
c
      k=k+1 !  'AIJ009'
      IJ_PHI1K = k ! PHI1000 (M**2/S**2) 4 DA
      lname_ij(k) = ' PHI1000'
      units_ij(k) = 'm^2/s^2'
      name_ij(k) = 'PHI1K'
      ia_ij(k) = ia_dga
c
      k=k+1 !  'AIJ010'
      IJ_PHI850 = k ! PHI850 (M**2/S**2-1500*GRAV) 4 DA
      lname_ij(k) = ' PHI850'
      units_ij(k) = 'm^2/s^2'
      name_ij(k) = 'PHI850'
      ia_ij(k) = ia_dga
c
      k=k+1 !  'AIJ011'
      IJ_PHI700 = k ! PHI700-3000*GRAV  4 DA
      lname_ij(k) = 'PHI700'
      units_ij(k) = 'm^2/s^2'
      name_ij(k) = 'PHI700'
      ia_ij(k) = ia_dga
c
      k=k+1 !  'AIJ012'
      IJ_PHI500 = k ! PHI500-5600*GRAV  4 DA
      lname_ij(k) = 'PHI500'
      units_ij(k) = 'm^2/s^2'
      name_ij(k) = 'PHI500'
      ia_ij(k) = ia_dga
c
      k=k+1 !  'AIJ013'
      IJ_PHI300 = k ! PHI300-9500*GRAV  4 DA
      lname_ij(k) = 'PHI300'
      units_ij(k) = 'm^2/s^2'
      name_ij(k) = 'PHI300'
      ia_ij(k) = ia_dga
c
      k=k+1 !  'AIJ014'
      IJ_PHI100 = k ! PHI100-16400*GRAV 4 DA
      lname_ij(k) = 'PHI100'
      units_ij(k) = 'm^2/s^2'
      name_ij(k) = 'PHI100'
      ia_ij(k) = ia_dga
c
      k=k+1 !  'AIJ015'
      IJ_PHI30 = k ! PHI30-24000*GRAV   4 DA
      lname_ij(k) = 'PHI30'
      units_ij(k) = 'm^2/s^2'
      name_ij(k) = 'PHI30'
      ia_ij(k) = ia_dga
c
      k=k+1 !  'AIJ016'
      IJ_T850 = k ! T850-TF (K-TF)*GRAV) (NO PRT) 4 DA
      lname_ij(k) = 'T850'
      units_ij(k) = 'degC'
      name_ij(k) = 'T850'
      ia_ij(k) = ia_dga
c
      k=k+1 !  'AIJ017'
      IJ_PMCCLD = k ! PCLDMC (1)  (COMPOSITE OVER ATMOSPHERE) 2 RD
      lname_ij(k) = 'ACCUMULATED CONVECTIVE CLOUD COVER'
      units_ij(k) = '1'
      name_ij(k) = 'PMCCLD'
      ia_ij(k) = ia_rad
c
      k=k+1 !  'AIJ018'
      IJ_CLDTPPR = k ! P-CLOUD TOP   (100 PA)                  2 RD
      lname_ij(k) = 'ACCUMULATED CLOUD TOP PRESSURE'
      units_ij(k) = '100 PA'
      name_ij(k) = 'CLDTPPR'
      ia_ij(k) = ia_rad
c
      k=k+1 !  'AIJ019'
      IJ_CLDCV = k ! PCLD (1)  (COMPOSITE OVER ATMOSPHERE)   2 RD
      lname_ij(k) = 'ACCUMULATED CLOUD COVER'
      units_ij(k) = '1'
      name_ij(k) = 'CLDCV'
      ia_ij(k) = ia_rad
c
      k=k+1 !  'AIJ020'
      IJ_PEV  = k
       ! 16*P4*(SHA*T4+Z4)*V1*DSIG*DXV (100 W*M/S**2) (UV GRID) 4 DA
      lname_ij(k) = '16*P4*(SHA*T4+Z4)*V1*DSIG*DXV'
      units_ij(k) = '100 W*m/s^2'
      name_ij(k) = 'PEV'
      ia_ij(k) = ia_dga
c
      k=k+1 !  'AIJ021'
      IJ_TRNFP0 = k ! TRNFP0 (W/M**2)                         2 RS
      lname_ij(k) = 'ACCUMULATED NET LONGWAVE FLUX, TOA'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'TRNFP0'
      ia_ij(k) = ia_rad
c
      k=k+1 !  'AIJ022'
      IJ_SRTR = k ! SRHDT+TRHDT (J/M**2)                    1 SF
      lname_ij(k) = 'ACCUMULATED RADIATION ABS BY SURF'
      units_ij(k) = 'J/m^2'
      name_ij(k) = 'SRTR'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ023'
      IJ_NETH = k    ! SRHDT+TRHDT+SHDT+EVHDT+ENRGP (J/M**2)   1 SC
      lname_ij(k) = 'ACCUMLATED NET SURFACE HEATING'
      units_ij(k) = 'J/m^2'
      name_ij(k) = 'NETH'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ024'
      IJ_SRNFP0 = k ! SRNFP0 (W/M**2)                         2 RD
      lname_ij(k) = 'ACCUMULATED NET SHORTWAVE FLUX, TOA'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SRNFP0'
      ia_ij(k) = ia_rad
c
      k=k+1 !  'AIJ025'
      IJ_SRINCP0 = k ! SRINCP0 (W/M**2)                        2 RD
      lname_ij(k) = 'ACCUMULATED INCIDENT SHORTWAVE FLUX, TOA'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SRINCP0'
      ia_ij(k) = ia_rad
c
      k=k+1 !  'AIJ026'
      IJ_SRNFG = k   ! SRNFG (W/M**2)                          2 RD
      lname_ij(k) = 'ACCUMULATED NET SHORTWAVE FLUX, SURF'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SRNFG'
      ia_ij(k) = ia_rad
c
      k=k+1 !  'AIJ027'
      IJ_SRINCG = k ! SRINCG (W/M**2)                         2 RD
      lname_ij(k) = 'ACCUMULATED INCIDENT SHORTWAVE FLUX, SURF'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SRINCG'
      ia_ij(k) = ia_rad
c
      k=k+1 !  'AIJ028'
      IJ_TG1  = k ! TG1 (K-TF)                                1 GD
      lname_ij(k) = 'ACCUMULATED GROUND TEMPERATURE'
      units_ij(k) = 'degC'
      name_ij(k) = 'TG1'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ029'
      IJ_RSIT = k ! POICE+PLICE+(IF SNOW)PEARTH               4 DA
      lname_ij(k) = 'POICE+PLICE+(IF SNOW)PEARTH'
      units_ij(k) = '1'
      name_ij(k) = 'RSIT'
      ia_ij(k) = ia_dga
c
      k=k+1 !  'AIJ030'
      IJ_TDSL = k ! DIURNAL DELTA TS (K) OVER SOIL (NO PRT) .5*9 MN
      lname_ij(k) = 'DIURNAL DELTA TS (k) OVER SOIL'
      units_ij(k) = 'K'
      name_ij(k) = 'TDSL'
      ia_ij(k) = ia_12hr
c
      k=k+1 !  'AIJ031'
      IJ_DTDP = k ! DTHETA/DPHI (K S**2/M**2) IN TROPOSPHERE  4 DA
      lname_ij(k) = 'DTHETA/DPHI IN TROPOSPHERE'
      units_ij(k) = 'K s^2/m^2'
      name_ij(k) = 'DTDP'
      ia_ij(k) = ia_dga
c
      k=k+1 !  'AIJ032'
      IJ_RUNE = k ! RUN1 OVER EARTH  (KG/M**2)                1 PG
      lname_ij(k) = 'RUN1 OVER EARTH'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'RUNE'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ033'
      IJ_RUNLI = k ! RUN1 OVER LAND ICE  (KG/M**2) (NO PRT)    1 PG
      lname_ij(k) = 'RUN1 OVER LAND ICE'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'RUNLI'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ034'
      IJ_WS   = k ! SURFACE WIND SPEED (M/S)                  3 SF
      lname_ij(k) = 'ACCUMULATED SURFACE WIND SPEED'
      units_ij(k) = 'm/s'
      name_ij(k) = 'WS'
      ia_ij(k) = ia_srf
c
      k=k+1 !  'AIJ035'
      IJ_TS   = k ! TS (K-TF)                                 3 SF
      lname_ij(k) = 'ACCUMULATED SURFACE AIR TEMPERATURE'
      units_ij(k) = 'degC'
      name_ij(k) = 'TS'
      ia_ij(k) = ia_srf
c
      k=k+1 !  'AIJ036'
      IJ_US   = k ! US (M/S)                                  3 SF
      lname_ij(k) = 'ACCUMULATED EASTWARD COMP. OF SURFACE WIND'
      units_ij(k) = 'm/s'
      name_ij(k) = 'US'
      ia_ij(k) = ia_srf
c
      k=k+1 !  'AIJ037'
      IJ_VS   = k ! VS (M/S)                                  3 SF
      lname_ij(k) = 'ACCUMULATED NORTHWARD COMP. OF SURFACE WIND'
      units_ij(k) = 'm/s'
      name_ij(k) = 'VS'
      ia_ij(k) = ia_srf
c
      k=k+1 !  'AIJ038'
      IJ_SLP  = k ! PSL (100 PA-1000)  (USING TS)             4 DA
      lname_ij(k) = 'PSL USING TS'
      units_ij(k) = 'mb-1000'
      name_ij(k) = 'SLP'
      ia_ij(k) = ia_dga
c
      k=k+1 !  'AIJ039'
      IJ_UJET = k ! UJET (M/S)                                4 DA
      lname_ij(k) = 'UJET'
      units_ij(k) = 'm/s'
      name_ij(k) = 'UJET'
      ia_ij(k) = ia_dga
c
      k=k+1 !  'AIJ040'
      IJ_VJET = k ! VJET (M/S)                                4 DA
      lname_ij(k) = 'VJET'
      units_ij(k) = 'm/s'
      name_ij(k) = 'VJET'
      ia_ij(k) = ia_dga
c
      k=k+1 !  'AIJ041'
      IJ_PCLDL = k ! PCLD(LOW) (1)                            2 RD
      lname_ij(k) = 'ACCUMULATED LOW CLOUD COVER'
      units_ij(k) = '1'
      name_ij(k) = 'PCLDL'
      ia_ij(k) = ia_rad
c
      k=k+1 !  'AIJ042'
      IJ_PCLDM = k ! PCLD(MID) (1)                            2 RD
      lname_ij(k) = 'ACCUMULATED MIDDLE CLOUD COVER'
      units_ij(k) = '1'
      name_ij(k) = 'PCLDM'
      ia_ij(k) = ia_rad
c
      k=k+1 !  'AIJ043'
      IJ_PCLDH = k ! PCLD(HIGH) (1)                           2 RD
      lname_ij(k) = 'ACCUMULATED HIGH CLOUD COVER'
      units_ij(k) = '1'
      name_ij(k) = 'PCLDH'
      ia_ij(k) = ia_rad
c
      k=k+1 !  'AIJ044'
      IJ_BTMPW = k ! BTEMPW-TF (K-TF)                         2 RD
      lname_ij(k) = 'ACCUMULATED BRIGHTNESS TEMPERATURE IN WINDOW'
      units_ij(k) = 'degC'
      name_ij(k) = 'BTMPW'
      ia_ij(k) = ia_rad
c
      k=k+1 !  'AIJ045'
      IJ_SRREF = k ! PLAVIS*S0*COSZ (W/M**2)                  2 RD
      lname_ij(k) = 'ACCUMULATED INC. SHORTWAVE x VISIBLE ALBEDO'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SRREF'
      ia_ij(k) = ia_rad
c
      k=k+1 !  'AIJ046'
      IJ_TOC2 = k ! TGO2= TOCEAN(2)  (C)                   .5*9 MN
      lname_ij(k) = 'TGO2=ODATA(4)'
      units_ij(k) = 'degC'
      name_ij(k) = 'TOC2'
      ia_ij(k) = ia_12hr
c
      k=k+1 !  'AIJ047'
      IJ_TAUS = k ! TAUS  (MOM. SURF. DRAG) (kg/m**2) (NO PRT)  3 SF
      lname_ij(k) = 'ACCUMULATED MAG. OF SURFACE STRESS'
      units_ij(k) = 'kg/m*s^2'
      name_ij(k) = 'TAUS'
      ia_ij(k) = ia_srf
c
      k=k+1 !  'AIJ048'
      IJ_TAUUS = k ! TAUUS (MOM. SURF. DRAG) (kg/m**2) (NO PRT)  3 SF
      lname_ij(k) = 'ACCUMULATED EASTWARD COMP. OF SURFACE STRESS'
      units_ij(k) = 'kg/m*s^2'
      name_ij(k) = 'TAUUS'
      ia_ij(k) = ia_srf
c
      k=k+1 !  'AIJ049'
      IJ_TAUVS = k ! TAUVS (MOM. SURF. DRAG) (kg/m**2) (NO PRT)  3 SF
      lname_ij(k) = 'ACCUMULATED NORTHWARD COMP. OF SURFACE STRESS'
      units_ij(k) = 'kg/m*s^2'
      name_ij(k) = 'TAUVS'
      ia_ij(k) = ia_srf
c
      k=k+1 !  'AIJ050'
      IJ_GWTR = k ! WATER1+WATER2+ICE1+ICE2 (EARTH POINTS ONLY) 1 GD
      lname_ij(k) = 'WATER1+WATER2+ICE1+ICE2'
      units_ij(k) = 'kg/m^2'
      name_ij(k) = 'GWTR'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ051'
      IJ_QS   = k ! QS                                (NO PRT)  3 SF
      lname_ij(k) = 'ACCUMULATED SURFACE AIR SPECIFIC HUMIDITY'
      units_ij(k) = '1'
      name_ij(k) = 'QS'
      ia_ij(k) = ia_srf
c
      k=k+1 !  'AIJ052'
      IJ_STRNGTS = k ! MAX(0,33-1.8*DAILY MEAN ON TS IN C)    .5*9 MN
      lname_ij(k) = 'MAX(0,33-1.8*DAILY MEAN ON TS IN C)'
      units_ij(k) = 'Fahrenheit*day'
      name_ij(k) = 'STRNGTS'
      ia_ij(k) = ia_12hr
c
      k=k+1 !  'AIJ053'
      IJ_ARUNU = k ! 40.6+.72*(2TS(C)-(QSATS-QS)*LHA/SHA)       3 SF
      lname_ij(k) = 'UNDERGROUND RUNOFF'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'ARUNU'
      ia_ij(k) = ia_srf
c
      k=k+1 !  'AIJ054'
      IJ_DTGDTS = k ! 18*(DEL(TG)/DEL(TS)-1),DEL= DIURN MX-MN .5*9 MN
      lname_ij(k) = '18*(DEL(TG)/DEL(TS)-1)'
      units_ij(k) = 'unknown'
      name_ij(k) = 'DTGDTS'
      ia_ij(k) = ia_12hr
c
      k=k+1 !  'AIJ055'
      IJ_PUQ  = k ! 8*P*U*Q (VERT. INTEGRATED) (12.5 PA*M/S) 4 DA
      lname_ij(k) = '8*P*U*Q (VERTICALLY INTEGRATED)'
      units_ij(k) = '12.5 PA*m/s'
      name_ij(k) = 'PUQ'
      ia_ij(k) = ia_dga
c
      k=k+1 !  'AIJ056'
      IJ_PVQ  = k ! 8*P*V*Q (VERT. INTEGRATED) (12.5 PA*M/S) 4 DA
      lname_ij(k) = '8*P*V*Q (VERTICALLY INTEGRATED)'
      units_ij(k) = '12.5 PA*m/s'
      name_ij(k) = 'PVQ'
      ia_ij(k) = ia_dga
c
      k=k+1 !  'AIJ057'
      IJ_TGO  = k ! TGO= TOCEAN(1)  (C)                      1 GD
      lname_ij(k) = 'TGO=ODATA(1)'
      units_ij(k) = 'degC'
      name_ij(k) = 'TGO'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ058'
      IJ_MSI2 = k ! ACE2OI= MSI2*POICE  (KG/M**2)        1 GD
      lname_ij(k) = 'ACE2OI=ODATA(3)*POICE'
      units_ij(k) = 'kg/m^2'
      name_ij(k) = 'MSI2'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ059'
      IJ_WLM  = k ! WIND SPEED IN TOP LAYER (M/S) before SDRAG 1 SD
      lname_ij(k) = 'WIND SPEED IN TOP LAYER'
      units_ij(k) = 'm/s'
      name_ij(k) = 'WLM'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ060'
      IJ_TGO2 = k ! TGO12= TOCEAN(3) (C)                  .5*9 MN
      lname_ij(k) = 'TGO12=ODATA(5)'
      units_ij(k) = 'degC'
      name_ij(k) = 'TGO2'
      ia_ij(k) = ia_12hr
c
      k=k+1 !  'AIJ061'
      IJ_EVAPO = k ! EVAP*POCEAN  (KG/M**2)                  1 GD
      lname_ij(k) = 'EVAP*POCEAN'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'EVAPO'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ062'
      IJ_EVAPI = k ! EVAP*POICE  (KG/M**2)                   1 GD
      lname_ij(k) = 'EVAP*POICE'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'EVAPI'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ063'
      IJ_EVAPLI = k ! EVAP OVER LAND ICE  (KG/M**2)          1 GD
      lname_ij(k) = 'EVAP OVER LAND ICE'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'EVAPLI'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ064'
      IJ_EVAPE = k ! EVAP OVER EARTH  (KG/M**2)              1 GD
      lname_ij(k) = 'EVAP OVER EARTH'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'EVAPE'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ065'
      IJ_F0OC = k ! F0DT*POCEAN, NET HEAT AT Z0  (J/M**2)    1 GD
      lname_ij(k) = 'F0DT*POCEAN, NET HEAT AT Z0'
      units_ij(k) = 'J/m^2'
      name_ij(k) = 'F0OC'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ066'
      IJ_F0OI = k ! F0DT*POICE, NET HEAT AT Z0  (J/M**2)     1 GD
      lname_ij(k) = 'F0DT*POICE, NET HEAT AT Z0'
      units_ij(k) = 'J/m^2'
      name_ij(k) = 'F0OI'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ067'
      IJ_F0LI = k ! F0DT, NET HEAT AT Z0 OVER LAND ICE  (J/M**2) 1 GD
      lname_ij(k) = 'F0DT, NET HEAT AT Z0 OVER LAND ICE'
      units_ij(k) = 'J/m^2'
      name_ij(k) = 'F0LI'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ068'
      IJ_F0E  = k ! F0DT, NET HEAT AT Z0 OVER EARTH  (J/M**2) 1 GD
      lname_ij(k) = ' F0DT, NET HEAT AT Z0 OVER EARTH'
      units_ij(k) = 'J/m^2'
      name_ij(k) = 'F0E'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ069'
      IJ_F1LI = k ! F1DT OVER LAND ICE  (J/M**2)             1 PG
      lname_ij(k) = 'F1DT OVER LAND ICE'
      units_ij(k) = 'J/m^2'
      name_ij(k) = 'F1LI'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ070'
      IJ_SNWF = k ! SNOW FALL  (KG/M**2)                     1 PR
      lname_ij(k) = 'SNOW FALL (H2O EQUIV)'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'SNWF'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ071'
      IJ_TSLI = k ! SURF AIR TEMP OVER LAND ICE  (C)  NISURF*1 SF
      lname_ij(k) = 'SURF AIR TEMP OVER LAND ICE'
      units_ij(k) = 'degC'
      name_ij(k) = 'TSLI'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ072'
      IJ_ERUN2 = k ! F2DT OVER LAND ICE  (J/M**2)            1 PG
      lname_ij(k) = 'F2DT OVER LAND ICE'
      units_ij(k) = 'J/m^2'
      name_ij(k) = 'ERUN2'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ073'
      IJ_SHDTLI = k ! SHDT OVER LAND ICE  (J/M**2)           3 SF
      lname_ij(k) = 'SENS HEAT FLUX OVER LAND ICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SHDTLI'
      ia_ij(k) = ia_srf
c
      k=k+1 !  'AIJ074'
      IJ_EVHDT = k ! EVHDT OVER LAND ICE  (J/M**2)           3 SF
      lname_ij(k) = 'LATENT HEAT FLUX OVER LAND ICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'EVHDT'
      ia_ij(k) = ia_srf
c
      k=k+1 !  'AIJ075'
      IJ_TRHDT = k ! TRHDT OVER LAND ICE  (J/M**2)           3 SF
      lname_ij(k) = 'TRHDT OVER LAND ICE'
      units_ij(k) = 'J/m^2'
      name_ij(k) = 'TRHDT'
      ia_ij(k) = ia_srf
c
      k=k+1 !  'AIJ076'
      IJ_TMAX  = k ! MAX(COMPOSITE TS)                      12 SF
      lname_ij(k) = 'MAX(COMPOSITE TS)'
      units_ij(k) = 'degC'
      name_ij(k) = 'TMAX'
      ia_ij(k) = ia_nmo
c
      k=k+1 !  'AIJ077'
      IJ_TMIN  = k ! MIN(COMPOSITE TS)                      12 SF
      lname_ij(k) = 'MIN(COMPOSITE TS)'
      units_ij(k) = 'degC'
      name_ij(k) = 'TMIN'
      ia_ij(k) = ia_nmo
c
      k=k+1 !  'AIJ078'
      IJ_TMNMX  = k ! MIN(DIURNAL MAX OF COMPOSITE TS)      12 MN
      lname_ij(k) = 'MIN(DIURNAL MAX OF COMPOSITE TS)'
      units_ij(k) = 'degC'
      name_ij(k) = 'TMNMX'
      ia_ij(k) = ia_nmo
c
      k=k+1 !  'AIJ079'
      IJ_PEVAP = k ! POTENTIAL EVAPORATION (KG/M**2)         1 EA
      lname_ij(k) = 'POTENTIAL EVAPORATION'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'PEVAP'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ080'
      IJ_TMAXE = k ! MAX TS OVER EARTH FOR CURRENT DAY (K).5*9 MN
      lname_ij(k) = 'MAX TS OVER EARTH FOR CURRENT DAY'
      units_ij(k) = 'K'
      name_ij(k) = 'TMAXE'
      ia_ij(k) = ia_12hr
c
      k=k+1 !  'AIJ081'
      IJ_WMSUM = k ! LIQUID WATER PATH (kg/M**2)             1 CL
      lname_ij(k) = 'LIQUID WATER PATH'
      units_ij(k) = 'kg/m^2'
      name_ij(k) = 'WMSUM'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ082'
      IJ_PSCLD = k ! SHALLOW CONVECTIVE CLOUD COVER  (1)     1 CL
      lname_ij(k) = 'ACCUMULATED SHALLOW CONVECTIVE CLOUD COVER'
      units_ij(k) = '1'
      name_ij(k) = 'PSCLD'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ083'
      IJ_PDCLD = k ! DEEP CONVECTIVE CLOUD COVER     (1)     1 CL
      lname_ij(k) = 'ACCUMULATED DEEP CONVECTIVE CLOUD COVER'
      units_ij(k) = '1'
      name_ij(k) = 'PDCLD'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ084'
      IJ_DCNVFRQ = k ! DEEP CONVECTIVE CLOUD FREQUENCY (1)   1 CL
      lname_ij(k) = 'ACCUMULATED DEEP CONVECTIVE CLOUD OCCURENCE'
      units_ij(k) = '1'
      name_ij(k) = 'DCNVFRQ'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ085'
      IJ_SCNVFRQ = k ! SHALLOW CONVECTIVE CLOUD FREQUENCY (1) 1 CL
      lname_ij(k) = 'ACCUMULATED SHALLOW CONVECTIVE CLOUD OCCURENCE'
      units_ij(k) = '1'
      name_ij(k) = 'SCNVFRQ'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ086'
      IJ_EMTMOM = k ! INCIDENT MTN EAST MOM. FLUX (MB-M/S**2)  1 SD
      lname_ij(k) = 'INCIDENT MTN EAST MOMENTUM FLUX'
      units_ij(k) = 'mb m/s^2'
      name_ij(k) = 'EMTMOM'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ087'
      IJ_SMTMOM = k ! INCIDENT MTN SOUTH MOM. FLUX (MB-M/S**2) 1 SD
      lname_ij(k) = 'INCIDENT MTN SOUTH MOMENTUM FLUX'
      units_ij(k) = 'mb m/s^2'
      name_ij(k) = 'SMTMOM'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ088'
      IJ_FPEU = k ! EAST-WEST POT. ENTHALPY FLUX (W)   /3600.*1 DY
      lname_ij(k) = 'EAST-WEST POTENTIAL ENTHALPY FLUX'
      units_ij(k) = 'W'
      name_ij(k) = 'FPEU'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ089'
      IJ_FPEV = k ! NORTH-SOUTH POT. ENTHALPY FLUX (W) /3600.*1 DY
      lname_ij(k) = 'NORTH-SOUTH POTENTIAL ENTHALPY FLUX'
      units_ij(k) = 'W'
      name_ij(k) = 'FPEV'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ090'
      IJ_FMU  = k ! EAST-WEST MASS FLUX (KG/S) 100./GRAV/3600.*1 DY
      lname_ij(k) = 'EAST-WEST MASS FLUX'
      units_ij(k) = 'kg/s'
      name_ij(k) = 'FMU'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ091'
      IJ_FMV  = k ! NORTH-SOUTH MASS FLUX (KG/S) 100./GRAV/3600.*1 DY
      lname_ij(k) = 'NORTH-SOUTH MASS FLUX'
      units_ij(k) = 'kg/s'
      name_ij(k) = 'FMV'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ092'
      IJ_FQU = k ! EAST-WEST WATER VAPOR FLUX (KG/S)   /3600.*1 DY
      lname_ij(k) = 'EAST-WEST WATER VAPOR FLUX'
      units_ij(k) = 'kg/s'
      name_ij(k) = 'FQU'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ093'
      IJ_FQV = k ! NORTH-SOUTH WATER VAPOR FLUX (KG/S) /3600.*1 DY
      lname_ij(k) = 'NORTH-SOUTH WATER VAPOR FLUX'
      units_ij(k) = 'kg/s'
      name_ij(k) = 'FQV'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ094'
      IJ_FGZU = k ! EAST-WEST GEOPOTENTIAL FLUX (W)    /3600.*1 DY
      lname_ij(k) = 'EAST-WEST GEOPOTENTIAL FLUX'
      units_ij(k) = 'W'
      name_ij(k) = 'FGZU'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ095'
      IJ_FGZV = k ! NORTH-SOUTH GEOPOTENTIAL FLUX (W)  /3600.*1 DY
      lname_ij(k) = 'NORTH-SOUTH GEOPOTENTIAL FLUX'
      units_ij(k) = 'W'
      name_ij(k) = 'FGZV'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ096'
      IJ_ERVR = k ! Energy Outflow by Rivers (10**10 W) E-10/DTS*1 RV
      lname_ij(k) = 'unknown'
      units_ij(k) = 'unknown'
      name_ij(k) = 'ERVR'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ097'
      IJ_MRVR = k ! Mass Outflow by Rivers (10**5 kg/s)  E-5/DTS*1 RV
      lname_ij(k) = 'unknown'
      units_ij(k) = 'unknown'
      name_ij(k) = 'MRVR'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ098'
      IJ_SDRAG = k ! DU/DT BY SDRAG (M S-2)                       1 SD
      lname_ij(k) = 'unknown'
      units_ij(k) = 'unknown'
      name_ij(k) = 'SDRAG'
      ia_ij(k) = ia_src
c
      k=k+1 !  'AIJ099'
      IJ_LKON = k
      lname_ij(k) = 'LAST DAY OF ICE-FREE LAKE'
      units_ij(k) = 'JULIAN DAY'
      name_ij(k) = 'LKONDAY'
      ia_ij(k) = ia_nmo
c
      k=k+1 !  'AIJ100'
      IJ_LKOFF = k 
      lname_ij(k) = 'LAST DAY OF ICED-UP LAKE'
      units_ij(k) = 'JULIAN DAY'
      name_ij(k) = 'LKOFFDAY'
      ia_ij(k) = ia_nmo
c
      return
      end subroutine ij_defs

      subroutine il_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      name_il(k) = 'AIL01'
      lname_il(k) = 'U (SUM FOR J=JEQ+1,JEQ,JEQ-1,JEQ-2) (PU GRID)'
      units_il(k) = 'm/s'
c
      k=k+1
      name_il(k) = 'AIL02'
      lname_il(k) = 'U (SUM FOR J=JEQ+1,JEQ,JEQ-1,JEQ-2) (PU GRID)'
      units_il(k) = 'm/s'
c
      k=k+1
      name_il(k) = 'AIL03'
      lname_il(k) = 'SD (SUM FOR J=JEQ,JEQ-1,JEQ-2)'
      units_il(k) = '100 N/s'
c
      k=k+1
      name_il(k) = 'AIL04'
      lname_il(k) = 'TX (SUM FOR J=JEQ,JEQ-1,JEQ-2)'
      units_il(k) = 'degC'
c
      k=k+1
      name_il(k) = 'AIL05'
      lname_il(k) = 'RH (SUM FOR J=JEQ,JEQ-1,JEQ-2)'
      units_il(k) = '1'
c
      k=k+1
      name_il(k) = 'AIL06'
      lname_il(k) = 'DTX(MC)*P*DA (SUM FOR J=JEQ,JEQ-1,JEQ-2)'
      units_il(k) = '100 K*N'
c
      k=k+1
      name_il(k) = 'AIL07'
      lname_il(k) = '(SRHR+TRHR)*DA (SUM FOR J=JEQ,JEQ-1,JEQ-2)'
      units_il(k) = 'W'
c
      k=k+1
      name_il(k) = 'AIL08'
      lname_il(k) = 'unknown'
      units_il(k) = 'unknown'
c
      k=k+1
      name_il(k) = 'AIL09'
      lname_il(k) = 'SD (AT LAT 50 N) (COMMENTED OUT)'
      units_il(k) = '100 N/s'
c
      k=k+1
      name_il(k) = 'AIL10'
      lname_il(k) = 'TX-273.16  (AT LAT 50 N)'
      units_il(k) = 'unknown'
c
      k=k+1
      name_il(k) = 'AIL11'
      lname_il(k) = 'SR+TR  (AT LAT 50 N)'
      units_il(k) = 'unknown'
c
      k=k+1
      name_il(k) = 'AIL12'
      lname_il(k) = '2*U  (AT LAT 50 N)'
      units_il(k) = 'unknown'
c
      k=k+1
      name_il(k) = 'AIL13'
      lname_il(k) = 'SD  (AT LAT 70 N) (COMMENTED OUT)'
      units_il(k) = 'unknown'
c
      k=k+1
      name_il(k) = 'AIL14'
      lname_il(k) = 'TX-273.16  (AT LAT 70 N)  (COMMENTED OUT)'
      units_il(k) = 'unknown'
c
      k=k+1
      name_il(k) = 'AIL15'
      lname_il(k) = 'SR+TR  (AT LAT 70 N)'
      units_il(k) = 'unknown'
c
      k=k+1
      name_il(k) = 'AIL16'
      lname_il(k) = '2*U  (AT LAT 70 N)        (COMMENTED OUT)'
      units_il(k) = 'unknown'
c
      return
      end subroutine il_defs

      subroutine jl_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      name_jl(k) = 'AJL01'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL02'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL03'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL04'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL05'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL06'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL07'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL08'
      lname_jl(k) = 'FMX(MC)*P'
      units_jl(k) = '100 PA'
c
      k=k+1
      name_jl(k) = 'AJL09'
      lname_jl(k) = 'SRHR'
      units_jl(k) = 'W/m^2'
c
      k=k+1
      name_jl(k) = 'AJL10'
      lname_jl(k) = 'TRHR'
      units_jl(k) = 'W/m^2'
c
      k=k+1
      name_jl(k) = 'AJL11'
      lname_jl(k) = 'DTX(SS)*P'
      units_jl(k) = '100 K*PA'
c
      k=k+1
      name_jl(k) = 'AJL12'
      lname_jl(k) = 'DT(DC)*P'
      units_jl(k) = '100 K*PA'
c
      k=k+1
      name_jl(k) = 'AJL13'
      lname_jl(k) = 'DT(MC)*P  DRY HEATING'
      units_jl(k) = '100 PA*K'
c
      k=k+1
      name_jl(k) = 'AJL14'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL15'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL16'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL17'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL18'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL19'
      lname_jl(k) = 'PCLD*P (TOTAL)'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL20'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL21'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL22'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL23'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL24'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL25'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL26'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL27'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL28'
      lname_jl(k) = 'PCLD*P (SS)'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL29'
      lname_jl(k) = 'PCLD*P (MC)'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL30'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL31'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL32'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL33'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL34'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL35'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL36'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL37'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL38'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL39'
      lname_jl(k) = 'DU(MC)*P (UV GRID)'
      units_jl(k) = '100 N/m/s'
c
      k=k+1
      name_jl(k) = 'AJL40'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL41'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL42'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL43'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL44'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL45'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL46'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL47'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL48'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL49'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL50'
      lname_jl(k) = 'DT(MC)*P  CHANGE OF PHASE'
      units_jl(k) = '100 PA*K'
c
      k=k+1
      name_jl(k) = 'AJL51'
      lname_jl(k) = 'CLHE*DQ(MC BEFORE COND)*P'
      units_jl(k) = '100 PA*K'
c
      k=k+1
      name_jl(k) = 'AJL52'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL53'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL54'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL55'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL56'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      name_jl(k) = 'AJL57'
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      return
      end subroutine jl_defs

      subroutine sjl_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      name_sjl(k) = 'ASJL01'
      lname_sjl(k) = 'TX'
      units_sjl(k) = 'degC'
c
      k=k+1
      name_sjl(k) = 'ASJL02'
      lname_sjl(k) = 'PHI'
      units_sjl(k) = 'm^2/s^2'
c
      k=k+1
      name_sjl(k) = 'ASJL03'
      lname_sjl(k) = 'SRHR'
      units_sjl(k) = 'W/m^2'
c
      k=k+1
      name_sjl(k) = 'ASJL04'
      lname_sjl(k) = 'TRHR'
      units_sjl(k) = 'W/m^2'
c
      return
      end subroutine sjl_defs

      subroutine jk_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      name_jk(k) = 'AJK01'
      lname_jk(k) = 'DP=(PDN-PM(K+1);PDN=MAX(PM(K+1),PS)'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK02'
      lname_jk(k) = 'DP4  (PDN-PM(K+1))   (UV GRID)'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK03'
      lname_jk(k) = '(TX-273.16)*DP'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK04'
      lname_jk(k) = 'PHI*DP'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK05'
      lname_jk(k) = 'Q*DP'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK06'
      lname_jk(k) = 'TH*DP'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK07'
      lname_jk(k) = 'RH*DP'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK08'
      lname_jk(k) = 'U*DP4  (UV GRID)'
      units_jk(k) = '100 PA*m/s'
c
      k=k+1
      name_jk(k) = 'AJK09'
      lname_jk(k) = 'V*DP4  (UV GRID)'
      units_jk(k) = '100 PA*m/s'
c
      k=k+1
      name_jk(k) = 'AJK10'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK11'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK12'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK13'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK14'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK15'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK16'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK17'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK18'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK19'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK20'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK21'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK22'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK23'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK24'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK25'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK26'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK27'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK28'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK29'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK30'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK31'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK32'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK33'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK34'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK35'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK36'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK37'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK38'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK39'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK40'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK41'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK42'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK43'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK44'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK45'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK46'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK47'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK48'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK49'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK50'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      name_jk(k) = 'AJK51'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      return
      end subroutine jk_defs

      subroutine ijg_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      name_ijg(k) = 'AIJG01'
      lname_ijg(k) = 'LAYER 1 REL SAT OF BARE AMD VEG. SOIL'
      units_ijg(k) = '%'
c
      k=k+1
      name_ijg(k) = 'AIJG02'
      lname_ijg(k) = 'LAYER 2 REL SAT OF BARE AMD VEG. SOIL'
      units_ijg(k) = '%'
c
      k=k+1
      name_ijg(k) = 'AIJG03'
      lname_ijg(k) = 'LAYER 3 REL SAT OF BARE AMD VEG. SOIL'
      units_ijg(k) = '%'
c
      k=k+1
      name_ijg(k) = 'AIJG04'
      lname_ijg(k) = 'LAYER 4 REL SAT OF BARE AMD VEG. SOIL'
      units_ijg(k) = '%'
c
      k=k+1
      name_ijg(k) = 'AIJG05'
      lname_ijg(k) = 'BETA (BARE SOIL)'
      units_ijg(k) = '%'
c
      k=k+1
      name_ijg(k) = 'AIJG06'
      lname_ijg(k) = 'BETA (PENMAN)'
      units_ijg(k) = '%'
c
      k=k+1
      name_ijg(k) = 'AIJG07'
      lname_ijg(k) = 'CANOPY  REL SATURATION'
      units_ijg(k) = '%'
c
      k=k+1
      name_ijg(k) = 'AIJG08'
      lname_ijg(k) = 'LAYER 1 REL SATURATION OF VEG. SOIL'
      units_ijg(k) = '%'
c
      k=k+1
      name_ijg(k) = 'AIJG09'
      lname_ijg(k) = 'LAYER 2 REL SATURATION OF VEG. SOIL'
      units_ijg(k) = '%'
c
      k=k+1
      name_ijg(k) = 'AIJG10'
      lname_ijg(k) = 'LAYER 3 REL SATURATION OF VEG. SOIL'
      units_ijg(k) = '%'
c
      k=k+1
      name_ijg(k) = 'AIJG11'
      lname_ijg(k) = 'BETA (BARE SOIL & VEGETATION)'
      units_ijg(k) = '%'
c
      k=k+1
      name_ijg(k) = 'AIJG12'
      lname_ijg(k) = 'CONDUCTANCE OF ATMOSPHERE'
      units_ijg(k) = '.01 m/s'
c
      k=k+1
      name_ijg(k) = 'AIJG13'
      lname_ijg(k) = 'CONDUCTANCE OF CANOPY'
      units_ijg(k) = '.01 m/s'
c
      k=k+1
      name_ijg(k) = 'AIJG14'
      lname_ijg(k) = 'PENMAN POTENTIAL EVAPORATION'
      units_ijg(k) = 'mm/day'
c
      k=k+1
      name_ijg(k) = 'AIJG15'
      lname_ijg(k) = 'TEMP OF LAYER 1 BARE SOIL AND SNOW'
      units_ijg(k) = 'degC'
c
      k=k+1
      name_ijg(k) = 'AIJG16'
      lname_ijg(k) = 'TEMP OF SOIL LAYER 2 - BARE SOIL'
      units_ijg(k) = 'degC'
c
      k=k+1
      name_ijg(k) = 'AIJG17'
      lname_ijg(k) = 'TEMP OF SOIL LAYER 3 - BARE SOIL'
      units_ijg(k) = 'degC'
c
      k=k+1
      name_ijg(k) = 'AIJG18'
      lname_ijg(k) = 'BARE SOIL EVAPORATION'
      units_ijg(k) = 'mm/day'
c
      k=k+1
      name_ijg(k) = 'AIJG19'
      lname_ijg(k) = 'DRY CANOPY EVAPORATION'
      units_ijg(k) = 'mm/day'
c
      k=k+1
      name_ijg(k) = 'AIJG20'
      lname_ijg(k) = 'WET CANOPY EVAPORATION'
      units_ijg(k) = 'mm/day'
c
      k=k+1
      name_ijg(k) = 'AIJG21'
      lname_ijg(k) = 'TEMP OF CANOPY AND SNOW'
      units_ijg(k) = 'degC'
c
      k=k+1
      name_ijg(k) = 'AIJG22'
      lname_ijg(k) = 'TEMP OF SOIL LAYER 1 - VEGETATED SOIL'
      units_ijg(k) = 'degC'
c
      k=k+1
      name_ijg(k) = 'AIJG23'
      lname_ijg(k) = 'TEMP OF SOIL LAYER 2 - VEGETATED SOIL'
      units_ijg(k) = 'degC'
c
      k=k+1
      name_ijg(k) = 'AIJG24'
      lname_ijg(k) = 'TEMP OF SOIL LAYER 3 - VEGETATED SOIL'
      units_ijg(k) = 'degC'
c
      k=k+1
      name_ijg(k) = 'AIJG25'
      lname_ijg(k) = 'AVERAGE WATER TABLE'
      units_ijg(k) = 'm'
c
      k=k+1
      name_ijg(k) = 'AIJG26'
      lname_ijg(k) = 'BETAV,OVER VEGETATION'
      units_ijg(k) = '%'
c
      k=k+1
      name_ijg(k) = 'AIJG27'
      lname_ijg(k) = 'BETAT,TRANSPIRATION'
      units_ijg(k) = '%'
c
      k=k+1
      name_ijg(k) = 'AIJG28'
      lname_ijg(k) = 'SNOW DEPTH OVER BARE SOIL'
      units_ijg(k) = 'MM H2O'
c
      k=k+1
      name_ijg(k) = 'AIJG29'
      lname_ijg(k) = 'SNOW DEPTH OVER VEG SOIL'
      units_ijg(k) = 'MM H2O'
c
      return
      end subroutine ijg_defs

      subroutine consrv_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      name_consrv(k) = 'CONSRV01'
      lname_consrv(k) = 'INSTANTANEOUS AM'
      units_consrv(k) = '10**9 J*s/m^2'
c
      k=k+1
      name_consrv(k) = 'CONSRV02'
      lname_consrv(k) = 'CHANGE OF AM BY ADVECTION'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV03'
      lname_consrv(k) = 'CHANGE OF AM BY CORIOLIS FORCE'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV04'
      lname_consrv(k) = 'CHANGE OF AM BY ADVEC + COR'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV05'
      lname_consrv(k) = 'CHANGE OF AM BY PRESSURE GRAD'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV06'
      lname_consrv(k) = 'CHANGE OF AM BY DYNAMICS'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV07'
      lname_consrv(k) = 'CHANGE OF AM BY SURFACE FRIC'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV08'
      lname_consrv(k) = 'CHANGE OF AM BY STRATOS DRAG'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV09'
      lname_consrv(k) = 'CHANGE OF AM BY FILTER'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV10'
      lname_consrv(k) = 'CHANGE OF AM BY DAILY RESTOR'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV11'
      lname_consrv(k) = 'SUM OF AM CHANGES'
      units_consrv(k) = '10**2 J/m^2'
c
      k=k+1
      name_consrv(k) = 'CONSRV12'
      lname_consrv(k) = 'INSTANTANEOUS KE'
      units_consrv(k) = '10**3 J/m^2'
c
      k=k+1
      name_consrv(k) = 'CONSRV13'
      lname_consrv(k) = 'CHANGE OF KE BY ADVECTION'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV14'
      lname_consrv(k) = 'CHANGE OF KE BY CORIOLIS FORCE'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV15'
      lname_consrv(k) = 'CHANGE OF KE BY ADVEC + COR'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV16'
      lname_consrv(k) = 'CHANGE OF KE BY PRESSURE GRAD'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV17'
      lname_consrv(k) = 'CHANGE OF KE BY DYNAMICS'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV18'
      lname_consrv(k) = 'CHANGE OF KE BY MOIST CONVEC'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV19'
      lname_consrv(k) = 'CHANGE OF KE BY SURF + DRY CONV'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV20'
      lname_consrv(k) = 'CHANGE OF KE BY STRATOS DRAG'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV21'
      lname_consrv(k) = 'CHANGE OF KE BY FILTER'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV22'
      lname_consrv(k) = 'CHANGE OF KE BY DAILY RESTOR'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV23'
      lname_consrv(k) = 'SUM OF KE CHANGES'
      units_consrv(k) = '10**-3 W/m^2'
c
      k=k+1
      name_consrv(k) = 'CONSRV24'
      lname_consrv(k) = 'INSTANTANEOUS MASS'
      units_consrv(k) = 'kg/m^2'
c
      k=k+1
      name_consrv(k) = 'CONSRV25'
      lname_consrv(k) = 'CHANGE OF MASS BY DYNAMICS'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV26'
      lname_consrv(k) = 'CHANGE OF MASS BY FILTER'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV27'
      lname_consrv(k) = 'CHANGE OF MASS BY DAILY RESTOR'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV28'
      lname_consrv(k) = 'SUM MASS CHANGES'
      units_consrv(k) = '10**-8 kg/m^2/s'
c
      k=k+1
      name_consrv(k) = 'CONSRV29'
      lname_consrv(k) = 'INSTANTANE TPE'
      units_consrv(k) = '10**5 J/m^2'
c
      k=k+1
      name_consrv(k) = 'CONSRV30'
      lname_consrv(k) = 'CHANGE OF TPE BY DYNAMICS'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV31'
      lname_consrv(k) = 'CHANGE OF TPE BY CONDENSATION'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV32'
      lname_consrv(k) = 'CHANGE OF TPE BY RADIATION'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV33'
      lname_consrv(k) = 'CHANGE OF TPE BY SURFACE INTER'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV34'
      lname_consrv(k) = 'CHANGE OF TPE BY FILTER'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV35'
      lname_consrv(k) = 'CHANGE OF TPE BY DAILY RESTOR'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV36'
      lname_consrv(k) = 'SUM OF TPE CHANGES'
      units_consrv(k) = '10**-2 W/m^2'
c
      k=k+1
      name_consrv(k) = 'CONSRV37'
      lname_consrv(k) = 'INSTANT WATER'
      units_consrv(k) = '10**-2 kg/m^2'
c
      k=k+1
      name_consrv(k) = 'CONSRV38'
      lname_consrv(k) = 'CHANGE OF WATER BY DYNAMICS'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV39'
      lname_consrv(k) = 'CHANGE OF WATER BY CONDENSATION'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV40'
      lname_consrv(k) = 'CHANGE OF WATER BY SURFACE EVAP'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV41'
      lname_consrv(k) = 'CHANGE OF WATER BY DAILY RESTOR'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV42'
      lname_consrv(k) = 'SUM WATER CHANGES'
      units_consrv(k) = '10**-8 kg/m^2/s'
c
      k=k+1
      name_consrv(k) = 'CONSRV43'
      lname_consrv(k) = 'unknown'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV44'
      lname_consrv(k) = 'unknown'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV45'
      lname_consrv(k) = 'unknown'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV46'
      lname_consrv(k) = 'unknown'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV47'
      lname_consrv(k) = 'unknown'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV48'
      lname_consrv(k) = 'unknown'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV49'
      lname_consrv(k) = 'unknown'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV50'
      lname_consrv(k) = 'unknown'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV51'
      lname_consrv(k) = 'unknown'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV52'
      lname_consrv(k) = 'unknown'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV53'
      lname_consrv(k) = 'unknown'
      units_consrv(k) = 'unknown'
c
      k=k+1
      name_consrv(k) = 'CONSRV54'
      lname_consrv(k) = 'unknown'
      units_consrv(k) = 'unknown'
c
      return
      end subroutine consrv_defs

      subroutine ijk_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      name_ijk(k) = 'UDPB'
      lname_ijk(k) = 'u-wind x delta p, b-grid'
      units_ijk(k) = '100 PA*m/s'
c
      k=k+1
      name_ijk(k) = 'VDPB'
      lname_ijk(k) = 'v-wind x delta p, b-grid'
      units_ijk(k) = '100 PA*m/s'
c
      k=k+1
      name_ijk(k) = 'DSEDPBx4'
      lname_ijk(k) = 'dry stat. energy x delta p x 4, b-grid'
      units_ijk(k) = '100 N/s^2'
c
      k=k+1
      name_ijk(k) = 'DPB'
      lname_ijk(k) = 'delta p, b-grid'
      units_ijk(k) = '100 PA'
c
      k=k+1
      name_ijk(k) = 'TDPBx4'
      lname_ijk(k) = 'temperature x delta p x 4, b-grid'
      units_ijk(k) = '100 K*PA'
c
      k=k+1
      name_ijk(k) = 'QDPBx4'
      lname_ijk(k) = 'spec. humidity x delta p x 4, b-grid'
      units_ijk(k) = '100 PA'
c
      return
      end subroutine ijk_defs

      subroutine ijl_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      name_ijl(k) = 'AIJL1'
      lname_ijl(k) = 'unknown'
      units_ijl(k) = 'unknown'
c
      k=k+1
      name_ijl(k) = 'AIJL2'
      lname_ijl(k) = 'unknown'
      units_ijl(k) = 'unknown'
c
      k=k+1
      name_ijl(k) = 'AIJL3'
      lname_ijl(k) = 'unknown'
      units_ijl(k) = 'unknown'
c
      k=k+1
      name_ijl(k) = 'AIJL4'
      lname_ijl(k) = 'unknown'
      units_ijl(k) = 'unknown'
c
      k=k+1
      name_ijl(k) = 'AIJL5'
      lname_ijl(k) = 'unknown'
      units_ijl(k) = 'unknown'
c
      return
      end subroutine ijl_defs

      subroutine wave_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      name_wave(k) = 'U850_EQ'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'V850_EQ'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'U300_EQ'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'V300_EQ'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'U050_EQ'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'V050_EQ'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'Z922_50N'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'Z700_50N'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'Z500_50N'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'Z300_50N'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'Z100_50N'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'Z010_50N'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      return
      end subroutine wave_defs

      subroutine jlsp_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      name_jlsp(k) = 'AJLSP01'
      lname_jlsp(k) = 'unknown'
      units_jlsp(k) = 'unknown'
c
      k=k+1
      name_jlsp(k) = 'AJLSP02'
      lname_jlsp(k) = 'unknown'
      units_jlsp(k) = 'unknown'
c
      k=k+1
      name_jlsp(k) = 'AJLSP03'
      lname_jlsp(k) = 'unknown'
      units_jlsp(k) = 'unknown'
c
      return
      end subroutine jlsp_defs

      subroutine pj_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      name_pj(k) = 'APJ1'
      lname_pj(k) = 'P'
      units_pj(k) = '100 PA'
c
      k=k+1
      name_pj(k) = 'APJ2'
      lname_pj(k) = '4*P4I (UV GRID)'
      units_pj(k) = '100 PA'
c
      k=k+1
      name_pj(k) = 'APJ3'
      lname_pj(k) = 'unknown'
      units_pj(k) = 'unknown'
c
      return
      end subroutine pj_defs

      subroutine tsf_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      name_tsf(k) = 'TSFREZ1'
      lname_tsf(k) = 'FIRST DAY OF GROWING SEASON'
      units_tsf(k) = 'JULIAN DAY'
      tf_day1 = k
c
      k=k+1
      name_tsf(k) = 'TSFREZ2'
      lname_tsf(k) = 'LAST DAY OF GROWING SEASON'
      units_tsf(k) = 'JULIAN DAY'
      tf_last = k
c
      k=k+1
      name_tsf(k) = 'LKICEON'
      lname_tsf(k) = 'LAST DAY OF ICE-FREE LAKE'
      units_tsf(k) = 'JULIAN DAY'
      tf_lkon = k
c
      k=k+1
      name_tsf(k) = 'LKICEOFF'
      lname_tsf(k) = 'LAST DAY OF ICED-UP LAKE'
      units_tsf(k) = 'JULIAN DAY'
      tf_lkoff = k
c
      return
      end subroutine tsf_defs
