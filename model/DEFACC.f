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
      call aj_defs
      dj_name=aj_name
      do k=1,kaj ! to avoid naming conflicts, put a prefix
         aj_name(k) = 'aj_'//trim(aj_name(k))
         dj_name(k) = 'reg_'//trim(dj_name(k))
      enddo
      call apj_defs
      call ajl_defs
      call asjl_defs
      call aij_defs
      call ail_defs
      call aijg_defs
      call consrv_defs
      call wave_defs
      call ajk_defs
      call aijk_defs
      call aijl_defs
      call ajlsp_defs
      return
      end subroutine def_acc

      subroutine iparm_defs
      use E001M12_COM
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
      use E001M12_COM
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

      subroutine aj_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c  AJ01
      k=k+1
      aj_name(k) = 'SRINCP0'
      aj_lname(k) = 'INC SW'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 2
c  AJ02
      k=k+1
      aj_name(k) = 'SRNFP0'
      aj_lname(k) = 'PLANETARY ALBEDO'
      aj_units(k) = '1'
      aj_ia(k) = 2
c  AJ03
      k=k+1
      aj_name(k) = 'SRNFP1'
      aj_lname(k) = 'SRNFP1'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 2
c  AJ04
      k=k+1
      aj_name(k) = 'SRABS'
      aj_lname(k) = 'SRABSATM'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 2
c  AJ05
      k=k+1
      aj_name(k) = 'SRINCG'
      aj_lname(k) = 'SW INC ON ZO'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 2
c  AJ06
      k=k+1
      aj_name(k) = 'SRNFG'
      aj_lname(k) = 'SW ABS AT ZO'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 2
c  AJ07
      k=k+1
      aj_name(k) = 'TRNFP0'
      aj_lname(k) = 'TRNFP0=AJ(74)+A2BYA1*AJ(9)/DTSRCE'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 2
c  AJ08
      k=k+1
      aj_name(k) = 'TRNFP1'
      aj_lname(k) = 'TRNFP1=AJ(75)+A2BYA1*AJ(9)/DTSRCE'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 2
c  AJ09
      k=k+1
      aj_name(k) = 'TRHDT'
      aj_lname(k) = 'NET LW AT ZO (positive downward)'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 1
c  AJ10
      k=k+1
      aj_name(k) = 'RNFP0'
      aj_lname(k) = 'RNFP0=AJ(2)+AJ(7)'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 2
c  AJ11
      k=k+1
      aj_name(k) = 'RNFP1'
      aj_lname(k) = 'RNFP1=AJ(3)+AJ(8)'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 2
c  AJ12
      k=k+1
      aj_name(k) = 'RHDT'
      aj_lname(k) = 'RHDT=A1BYA2*AJ(6)*DTSRCE+AJ(9)'
      aj_units(k) = 'J/m^2'
      aj_ia(k) = 1
c  AJ13
      k=k+1
      aj_name(k) = 'SHDT'
      aj_lname(k) = 'SENS HEAT FLUX'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 1
c  AJ14
      k=k+1
      aj_name(k) = 'EVHDT'
      aj_lname(k) = 'LATENT HEAT FLUX'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 1
c  AJ15
      k=k+1
      aj_name(k) = 'F2DT'
      aj_lname(k) = 'F2DT'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 1
c  AJ16
      k=k+1
      aj_name(k) = 'HZ1'
      aj_lname(k) = 'HEATZ1=AJ(41)+AJ(42)'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 1
c  AJ17
      k=k+1
      aj_name(k) = 'TG2'
      aj_lname(k) = 'TG2'
      aj_units(k) = 'degC'
      aj_ia(k) = 1
c  AJ18
      k=k+1
      aj_name(k) = 'TG1'
      aj_lname(k) = 'TG1'
      aj_units(k) = 'degC'
      aj_ia(k) = 1
c  AJ19
      k=k+1
      aj_name(k) = 'EVAP'
      aj_lname(k) = 'EVAP'
      aj_units(k) = 'mm/day'
      aj_ia(k) = 1
c  AJ20
      k=k+1
      aj_name(k) = 'PRCP'
      aj_lname(k) = 'PREC'
      aj_units(k) = 'mm/day'
      aj_ia(k) = 1
c  AJ21
      k=k+1
      aj_name(k) = 'TX'
      aj_lname(k) = 'TX'
      aj_units(k) = 'degC'
      aj_ia(k) = 4
c  AJ22
      k=k+1
      aj_name(k) = 'TX1'
      aj_lname(k) = 'TX1'
      aj_units(k) = 'degC'
      aj_ia(k) = 4
c  AJ23
      k=k+1
      aj_name(k) = 'TSRF'
      aj_lname(k) = 'TS'
      aj_units(k) = 'degC'
      aj_ia(k) = 3
c  AJ24
      k=k+1
      aj_name(k) = 'DTSGST'
      aj_lname(k) = 'DTH/DPHI  (STRATOSPHERE)'
      aj_units(k) = 'unknown'
      aj_ia(k) = 4
c  AJ25
      k=k+1
      aj_name(k) = 'DTDGTR'
      aj_lname(k) = 'DTH/DPHI  (TROPOSPHERE)'
      aj_units(k) = 'unknown'
      aj_ia(k) = 4
c  AJ26
      k=k+1
      aj_name(k) = 'RICST'
      aj_lname(k) = 'RICH NUM (STRATOSPHERE)'
      aj_units(k) = '1'
      aj_ia(k) = 4
c  AJ27
      k=k+1
      aj_name(k) = 'RICTR'
      aj_lname(k) = 'RICH NUM (TROPOSPHERE)'
      aj_units(k) = '1'
      aj_ia(k) = 4
c  AJ28
      k=k+1
      aj_name(k) = 'ROSST'
      aj_lname(k) = 'ROSS NUM (STRATOSPHERE)'
      aj_units(k) = '1'
      aj_ia(k) = 4
c  AJ29
      k=k+1
      aj_name(k) = 'ROSTR'
      aj_lname(k) = 'ROSS NUM (TROPOSPHERE)'
      aj_units(k) = '1'
      aj_ia(k) = 4
c  AJ30
      k=k+1
      aj_name(k) = 'RSI'
      aj_lname(k) = 'POICE'
      aj_units(k) = '1'
      aj_ia(k) = 1
c  AJ31
      k=k+1
      aj_name(k) = 'RSNOW'
      aj_lname(k) = 'PSNOW'
      aj_units(k) = '1'
      aj_ia(k) = 4
c  AJ32
      k=k+1
      aj_name(k) = 'SWCOR'
      aj_lname(k) = 'SW CORRECTION'
      aj_units(k) = 'unknown'
      aj_ia(k) = 2
c  AJ33
      k=k+1
      aj_name(k) = 'OHT'
      aj_lname(k) = 'OCEAN TRANSPORT'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 1
c  AJ34
      k=k+1
      aj_name(k) = 'OMLT'
      aj_lname(k) = 'OCEAN TEMP AT MAX. MIXED LAYER DEPTH'
      aj_units(k) = 'unknown'
      aj_ia(k) = 1
c  AJ35
      k=k+1
      aj_name(k) = 'DTDJS'
      aj_lname(k) = 'T(J+1)-T(J-1)  (SUM OVER STRATOSPHERE OF)'
      aj_units(k) = 'K'
      aj_ia(k) = 4
c  AJ36
      k=k+1
      aj_name(k) = 'DTDJT'
      aj_lname(k) = 'T(J+1)-T(J-1)  (SUM OVER TROPOSPHERE OF)'
      aj_units(k) = 'K'
      aj_ia(k) = 4
c  AJ37
      k=k+1
      aj_name(k) = 'LSTR'
      aj_lname(k) = 'SQRT(DTH/DLNP)/SINJ  (STRATOSPHERE)'
      aj_units(k) = 'unknown'
      aj_ia(k) = 4
c  AJ38
      k=k+1
      aj_name(k) = 'LTRO'
      aj_lname(k) = 'SQRT(DTH/DLNP)/SINJ'
      aj_units(k) = 'unknown'
      aj_ia(k) = 4
c  AJ39
      k=k+1
      aj_name(k) = 'EPRCP'
      aj_lname(k) = 'ENERGP'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 1
c  AJ40
      k=k+1
      aj_name(k) = 'ERUN1'
      aj_lname(k) = 'ERUN1'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 1
c  AJ41
      k=k+1
      aj_name(k) = 'EDIFS'
      aj_lname(k) = 'EDIFS'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 1
c  AJ42
      k=k+1
      aj_name(k) = 'F1DT'
      aj_lname(k) = 'F1DT'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 1
c  AJ43
      k=k+1
      aj_name(k) = 'ERUN2'
      aj_lname(k) = 'ERUN2'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 1
c  AJ44
      k=k+1
      aj_name(k) = 'HZ0'
      aj_lname(k) = 'NET HT Z0=AJ(12)+AJ(13)+AJ(14)+AJ(39)-AJ(40)'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 1
c  AJ45
      k=k+1
      aj_name(k) = 'DIFS'
      aj_lname(k) = 'DIFS'
      aj_units(k) = 'mm/day'
      aj_ia(k) = 1
c  AJ46
      k=k+1
      aj_name(k) = 'IMELT'
      aj_lname(k) = 'AIFO; BRUN2; CRUN2+CIFI'
      aj_units(k) = 'mm/day'
      aj_ia(k) = 1
c  AJ47
      k=k+1
      aj_name(k) = 'RUN2'
      aj_lname(k) = 'RUN2'
      aj_units(k) = 'mm/day'
      aj_ia(k) = 1
c  AJ48
      k=k+1
      aj_name(k) = 'DWTR2'
      aj_lname(k) = 'HEAT RUNOFF THROUGH THE MIXED LAYER DEPTH'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 1
c  AJ49
      k=k+1
      aj_name(k) = 'WTR1'
      aj_lname(k) = 'WTR1'
      aj_units(k) = 'kg/m^2'
      aj_ia(k) = 1
c  AJ50
      k=k+1
      aj_name(k) = 'ACE1'
      aj_lname(k) = 'ACE1'
      aj_units(k) = 'kg/m^2'
      aj_ia(k) = 1
c  AJ51
      k=k+1
      aj_name(k) = 'WTR2'
      aj_lname(k) = 'WTR2'
      aj_units(k) = 'kg/m^2'
      aj_ia(k) = 1
c  AJ52
      k=k+1
      aj_name(k) = 'ACE2'
      aj_lname(k) = 'ACE2'
      aj_units(k) = 'kg/m^2'
      aj_ia(k) = 1
c  AJ53
      k=k+1
      aj_name(k) = 'SNOW'
      aj_lname(k) = 'SNOW'
      aj_units(k) = 'kg/m^2'
      aj_ia(k) = 1
c  AJ54
      k=k+1
      aj_name(k) = 'RUN1'
      aj_lname(k) = 'RUN1'
      aj_units(k) = 'mm/day'
      aj_ia(k) = 1
c  AJ55
      k=k+1
      aj_name(k) = 'BRTEMP'
      aj_lname(k) = 'BTEMPW'
      aj_units(k) = 'degC'
      aj_ia(k) = 2
c  AJ56
      k=k+1
      aj_name(k) = 'HZ2'
      aj_lname(k) = 'HEATZ2=AJ(15)+AJ(43)'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 1
c  AJ57
      k=k+1
      aj_name(k) = 'PCLDSS'
      aj_lname(k) = 'PCLDSS (COMPOSITE OVER ATMOSPHERE)'
      aj_units(k) = '1'
      aj_ia(k) = 2
c  AJ58
      k=k+1
      aj_name(k) = 'PCLDMC'
      aj_lname(k) = 'PCLDMC (COMPOSITE OVER ATMOSPHERE)'
      aj_units(k) = '1'
      aj_ia(k) = 2
c  AJ59
      k=k+1
      aj_name(k) = 'PCLD'
      aj_lname(k) = 'PCLD (COMPOSITE OVER ATMOSPHERE)'
      aj_units(k) = '1'
      aj_ia(k) = 2
c  AJ60
      k=k+1
      aj_name(k) = 'CTOPP'
      aj_lname(k) = 'CLDTOPMC=AJ(80)/AJ(58)'
      aj_units(k) = '100 PA'
      aj_ia(k) = 0
c  AJ61
      k=k+1
      aj_name(k) = 'PRCPSS'
      aj_lname(k) = 'PRCPSS'
      aj_units(k) = 'mm/day'
      aj_ia(k) = 1
c  AJ62
      k=k+1
      aj_name(k) = 'PRCPMC'
      aj_lname(k) = 'PRCPMC'
      aj_units(k) = '100 PA'
      aj_ia(k) = 1
c  AJ63
      k=k+1
      aj_name(k) = 'QP'
      aj_lname(k) = 'Q*P (INTEGRAL OVER ATMOSPHERE OF)'
      aj_units(k) = 'kg/m^2'
      aj_ia(k) = 4
c  AJ64
      k=k+1
      aj_name(k) = 'GAM'
      aj_lname(k) = 'GAM  (*SIG(TROPOSPHERE)/GRAV)'
      aj_units(k) = 'K/m'
      aj_ia(k) = 4
c  AJ65
      k=k+1
      aj_name(k) = 'GAMM'
      aj_lname(k) = 'GAMM (SIG(TROPOSPHERE)/GAMD)'
      aj_units(k) = 'K s^2/m^2'
      aj_ia(k) = 4
c  AJ66
      k=k+1
      aj_name(k) = 'GAMC'
      aj_lname(k) = 'GAMC'
      aj_units(k) = 'K/m'
      aj_ia(k) = 4
c  AJ67
      k=k+1
      aj_name(k) = 'TRINCG'
      aj_lname(k) = 'LW INC ON ZO'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 2
c  AJ68
      k=k+1
      aj_name(k) = 'FTHERM'
      aj_lname(k) = 'ENERGY DIFFUSION INTO THERMOCLINE'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 9
c  AJ69
      k=k+1
      aj_name(k) = 'SPTYPE'
      aj_lname(k) = 'sum of surface type fraction'
      aj_units(k) = 'unknown'
      aj_ia(k) = 1
c  AJ70
      k=k+1
      aj_name(k) = 'HSURF'
      aj_lname(k) = 'TRNFP0-TRNFG'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 2
c  AJ71
      k=k+1
      aj_name(k) = 'HATM'
      aj_lname(k) = 'TRNFP1-TRNFG'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 2
c  AJ72
      k=k+1
      aj_name(k) = 'PLAVIS'
      aj_lname(k) = 'PLAVIS*S0*COSZ'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 2
c  AJ73
      k=k+1
      aj_name(k) = 'PLANIR'
      aj_lname(k) = 'PLANIR*S0*COSZ'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 2
c  AJ74
      k=k+1
      aj_name(k) = 'ALBVIS'
      aj_lname(k) = 'ALBVIS*S0*COSZ'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 2
c  AJ75
      k=k+1
      aj_name(k) = 'ALBNIR'
      aj_lname(k) = 'ALBNIR*S0*COSZ'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 2
c  AJ76
      k=k+1
      aj_name(k) = 'SRRVIS'
      aj_lname(k) = 'SRRVIS*S0*COSZ'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 2
c  AJ77
      k=k+1
      aj_name(k) = 'SRRNIR'
      aj_lname(k) = 'SRRNIR*S0*COSZ'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 2
c  AJ78
      k=k+1
      aj_name(k) = 'SRAVIS'
      aj_lname(k) = 'SRAVIS*S0*COSZ'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 2
c  AJ79
      k=k+1
      aj_name(k) = 'SRANIR'
      aj_lname(k) = 'SRANIR*S0*COSZ'
      aj_units(k) = 'W/m^2'
      aj_ia(k) = 2
c  AJ80
      k=k+1
      aj_name(k) = 'CDLDEP'
      aj_lname(k) = 'PBOTMC-PTOPMC'
      aj_units(k) = '100 PA'
      aj_ia(k) = 2
c
      k=k+1
      aj_name(k) = 'AJ81'
      aj_lname(k) = 'unknown'
      aj_units(k) = 'unknown'
      aj_ia(k) = 0
c
      k=k+1
      aj_name(k) = 'AJ82'
      aj_lname(k) = 'unknown'
      aj_units(k) = 'unknown'
      aj_ia(k) = 0
c
      k=k+1
      aj_name(k) = 'AJ83'
      aj_lname(k) = 'unknown'
      aj_units(k) = 'unknown'
      aj_ia(k) = 0
c
      k=k+1
      aj_name(k) = 'AJ84'
      aj_lname(k) = 'unknown'
      aj_units(k) = 'unknown'
      aj_ia(k) = 0
c
      k=k+1
      aj_name(k) = 'AJ85'
      aj_lname(k) = 'unknown'
      aj_units(k) = 'unknown'
      aj_ia(k) = 0
c
      k=k+1
      aj_name(k) = 'AJ86'
      aj_lname(k) = 'unknown'
      aj_units(k) = 'unknown'
      aj_ia(k) = 0
c
      k=k+1
      aj_name(k) = 'AJ87'
      aj_lname(k) = 'unknown'
      aj_units(k) = 'unknown'
      aj_ia(k) = 0
c
      k=k+1
      aj_name(k) = 'AJ88'
      aj_lname(k) = 'unknown'
      aj_units(k) = 'unknown'
      aj_ia(k) = 0
c
      k=k+1
      aj_name(k) = 'AJ89'
      aj_lname(k) = 'unknown'
      aj_units(k) = 'unknown'
      aj_ia(k) = 0
c
      k=k+1
      aj_name(k) = 'AJ90'
      aj_lname(k) = 'unknown'
      aj_units(k) = 'unknown'
      aj_ia(k) = 0
c
      k=k+1
      aj_name(k) = 'AJ91'
      aj_lname(k) = 'unknown'
      aj_units(k) = 'unknown'
      aj_ia(k) = 0
c
      k=k+1
      aj_name(k) = 'AJ92'
      aj_lname(k) = 'unknown'
      aj_units(k) = 'unknown'
      aj_ia(k) = 0
c
      k=k+1
      aj_name(k) = 'AJ93'
      aj_lname(k) = 'unknown'
      aj_units(k) = 'unknown'
      aj_ia(k) = 0
c
      k=k+1
      aj_name(k) = 'AJ94'
      aj_lname(k) = 'unknown'
      aj_units(k) = 'unknown'
      aj_ia(k) = 0
c
      return
      end subroutine aj_defs

      subroutine aij_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1 !  'AIJ001'
      aij_lname(k) = 'ACCUMULATED OCEAN ICE FRACTION'
      aij_units(k) = '1'
      aij_name(k) = 'RSOI'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ002'
      aij_lname(k) = 'ACCUMULATED SNOW FRACTION'
      aij_units(k) = '1'
      aij_name(k) = 'RSNW'
      aij_ia(k) = 4
c
      k=k+1 !  'AIJ003'
      aij_lname(k) = 'ACCUMULATED SNOW MASS'
      aij_units(k) = 'kg/m^2'
      aij_name(k) = 'SNOW'
      aij_ia(k) = 4
c
      k=k+1 !  'AIJ004'
      aij_lname(k) = 'INTEGRATED SENSIBLE HEAT FLUX'
      aij_units(k) = 'J/m^2'
      aij_name(k) = 'SHDT'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ005'
      aij_lname(k) = 'ACCUMULATED PRECIPITATION'
      aij_units(k) = 'kg/m^2'
      aij_name(k) = 'PREC'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ006'
      aij_lname(k) = 'ACCUMULATED EVAPORATION'
      aij_units(k) = 'kg/m^2'
      aij_name(k) = 'EVAP'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ007'
      aij_lname(k) = 'ACCUMULATED EVAPORATION EFFICIENCY'
      aij_units(k) = '1'
      aij_name(k) = 'BETA'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ008'
      aij_lname(k) = 'ACCUMULATED SURFACE PRESSURE'
      aij_units(k) = '100 PA'
      aij_name(k) = 'PRES'
      aij_ia(k) = 4
c
      k=k+1 !  'AIJ009'
      aij_lname(k) = ' PHI1000'
      aij_units(k) = 'm^2/s^2'
      aij_name(k) = 'PHI1K'
      aij_ia(k) = 4
c
      k=k+1 !  'AIJ010'
      aij_lname(k) = ' PHI850'
      aij_units(k) = 'm^2/s^2'
      aij_name(k) = 'PHI850'
      aij_ia(k) = 4
c
      k=k+1 !  'AIJ011'
      aij_lname(k) = 'PHI700'
      aij_units(k) = 'm^2/s^2'
      aij_name(k) = 'PHI700'
      aij_ia(k) = 4
c
      k=k+1 !  'AIJ012'
      aij_lname(k) = 'PHI500'
      aij_units(k) = 'm^2/s^2'
      aij_name(k) = 'PHI500'
      aij_ia(k) = 4
c
      k=k+1 !  'AIJ013'
      aij_lname(k) = 'PHI300'
      aij_units(k) = 'm^2/s^2'
      aij_name(k) = 'PHI300'
      aij_ia(k) = 4
c
      k=k+1 !  'AIJ014'
      aij_lname(k) = 'PHI100'
      aij_units(k) = 'm^2/s^2'
      aij_name(k) = 'PHI100'
      aij_ia(k) = 4
c
      k=k+1 !  'AIJ015'
      aij_lname(k) = 'PHI30'
      aij_units(k) = 'm^2/s^2'
      aij_name(k) = 'PHI30'
      aij_ia(k) = 4
c
      k=k+1 !  'AIJ016'
      aij_lname(k) = 'T850'
      aij_units(k) = 'degC'
      aij_name(k) = 'T850'
      aij_ia(k) = 4
c
      k=k+1 !  'AIJ017'
      aij_lname(k) = 'ACCUMULATED CONVECTIVE CLOUD COVER'
      aij_units(k) = '1'
      aij_name(k) = 'PCLDMC'
      aij_ia(k) = 2
c
      k=k+1 !  'AIJ018'
      aij_lname(k) = 'ACCUMULATED CLOUD TOP PRESSURE'
      aij_units(k) = '100 PA'
      aij_name(k) = 'CLDTPPR'
      aij_ia(k) = 2
c
      k=k+1 !  'AIJ019'
      aij_lname(k) = 'ACCUMULATED CLOUD COVER'
      aij_units(k) = '1'
      aij_name(k) = 'PCLD'
      aij_ia(k) = 2
c
      k=k+1 !  'AIJ020'
      aij_lname(k) = '16*P4*(SHA*T4+Z4)*V1*DSIG*DXV'
      aij_units(k) = '100 W*m/s^2'
      aij_name(k) = 'PEV'
      aij_ia(k) = 4
c
      k=k+1 !  'AIJ021'
      aij_lname(k) = 'ACCUMULATED NET LONGWAVE FLUX, TOA'
      aij_units(k) = 'W/m^2'
      aij_name(k) = 'TRNFP0'
      aij_ia(k) = 2
c
      k=k+1 !  'AIJ022'
      aij_lname(k) = 'ACCUMULATED RADIATION ABS BY SURF'
      aij_units(k) = 'J/m^2'
      aij_name(k) = 'SRTR'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ023'
      aij_lname(k) = 'ACCUMLATED NET SURFACE HEATING'
      aij_units(k) = 'J/m^2'
      aij_name(k) = 'NETH'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ024'
      aij_lname(k) = 'ACCUMULATED NET SHORTWAVE FLUX, TOA'
      aij_units(k) = 'W/m^2'
      aij_name(k) = 'SRNFP0'
      aij_ia(k) = 2
c
      k=k+1 !  'AIJ025'
      aij_lname(k) = 'ACCUMULATED INCIDENT SHORTWAVE FLUX, TOA'
      aij_units(k) = 'W/m^2'
      aij_name(k) = 'SRINCP0'
      aij_ia(k) = 2
c
      k=k+1 !  'AIJ026'
      aij_lname(k) = 'ACCUMULATED NET SHORTWAVE FLUX, SURF'
      aij_units(k) = 'W/m^2'
      aij_name(k) = 'SRNFG'
      aij_ia(k) = 2
c
      k=k+1 !  'AIJ027'
      aij_lname(k) = 'ACCUMULATED INCIDENT SHORTWAVE FLUX, SURF'
      aij_units(k) = 'W/m^2'
      aij_name(k) = 'SRINCG'
      aij_ia(k) = 2
c
      k=k+1 !  'AIJ028'
      aij_lname(k) = 'ACCUMULATED GROUND TEMPERATURE'
      aij_units(k) = 'degC'
      aij_name(k) = 'TG1'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ029'
      aij_lname(k) = 'POICE+PLICE+(IF SNOW)PEARTH'
      aij_units(k) = '1'
      aij_name(k) = 'RSIT'
      aij_ia(k) = 4
c
      k=k+1 !  'AIJ030'
      aij_lname(k) = 'DIURNAL DELTA TS (k) OVER SOIL'
      aij_units(k) = 'K'
      aij_name(k) = 'TDSL'
      aij_ia(k) = 9
c
      k=k+1 !  'AIJ031'
      aij_lname(k) = 'DTHETA/DPHI IN TROPOSPHERE'
      aij_units(k) = 'K s^2/m^2'
      aij_name(k) = 'DTDP'
      aij_ia(k) = 4
c
      k=k+1 !  'AIJ032'
      aij_lname(k) = 'RUN1 OVER EARTH'
      aij_units(k) = 'mm/day'
      aij_name(k) = 'RUNE'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ033'
      aij_lname(k) = 'RUN1 OVER LAND ICE'
      aij_units(k) = 'mm/day'
      aij_name(k) = 'RUNLI'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ034'
      aij_lname(k) = 'ACCUMULATED SURFACE WIND SPEED'
      aij_units(k) = 'm/s'
      aij_name(k) = 'WS'
      aij_ia(k) = 3
c
      k=k+1 !  'AIJ035'
      aij_lname(k) = 'ACCUMULATED SURFACE AIR TEMPERATURE'
      aij_units(k) = 'degC'
      aij_name(k) = 'TS'
      aij_ia(k) = 3
c
      k=k+1 !  'AIJ036'
      aij_lname(k) = 'ACCUMULATED EASTWARD COMP. OF SURFACE WIND'
      aij_units(k) = 'm/s'
      aij_name(k) = 'US'
      aij_ia(k) = 3
c
      k=k+1 !  'AIJ037'
      aij_lname(k) = 'ACCUMULATED NORTHWARD COMP. OF SURFACE WIND'
      aij_units(k) = 'm/s'
      aij_name(k) = 'VS'
      aij_ia(k) = 3
c
      k=k+1 !  'AIJ038'
      aij_lname(k) = 'PSL USING TS'
      aij_units(k) = 'mb-1000'
      aij_name(k) = 'SLP'
      aij_ia(k) = 4
c
      k=k+1 !  'AIJ039'
      aij_lname(k) = 'UJET'
      aij_units(k) = 'm/s'
      aij_name(k) = 'UJET'
      aij_ia(k) = 4
c
      k=k+1 !  'AIJ040'
      aij_lname(k) = 'VJET'
      aij_units(k) = 'm/s'
      aij_name(k) = 'VJET'
      aij_ia(k) = 4
c
      k=k+1 !  'AIJ041'
      aij_lname(k) = 'ACCUMULATED LOW CLOUD COVER'
      aij_units(k) = '1'
      aij_name(k) = 'PCLDL'
      aij_ia(k) = 2
c
      k=k+1 !  'AIJ042'
      aij_lname(k) = 'ACCUMULATED MIDDLE CLOUD COVER'
      aij_units(k) = '1'
      aij_name(k) = 'PCLDM'
      aij_ia(k) = 2
c
      k=k+1 !  'AIJ043'
      aij_lname(k) = 'ACCUMULATED HIGH CLOUD COVER'
      aij_units(k) = '1'
      aij_name(k) = 'PCLDH'
      aij_ia(k) = 2
c
      k=k+1 !  'AIJ044'
      aij_lname(k) = 'ACCUMULATED BRIGHTNESS TEMPERATURE IN WINDOW'
      aij_units(k) = 'degC'
      aij_name(k) = 'BTMPW'
      aij_ia(k) = 2
c
      k=k+1 !  'AIJ045'
      aij_lname(k) = 'ACCUMULATED INC. SHORTWAVE x VISIBLE ALBEDO'
      aij_units(k) = 'W/m^2'
      aij_name(k) = 'SRREF'
      aij_ia(k) = 2
c
      k=k+1 !  'AIJ046'
      aij_lname(k) = 'TGO2=ODATA(4)'
      aij_units(k) = 'degC'
      aij_name(k) = 'TOC2'
      aij_ia(k) = 9
c
      k=k+1 !  'AIJ047'
      aij_lname(k) = 'ACCUMULATED MAG. OF SURFACE STRESS'
      aij_units(k) = 'kg/m*s^2'
      aij_name(k) = 'TAUS'
      aij_ia(k) = 3
c
      k=k+1 !  'AIJ048'
      aij_lname(k) = 'ACCUMULATED EASTWARD COMP. OF SURFACE STRESS'
      aij_units(k) = 'kg/m*s^2'
      aij_name(k) = 'TAUUS'
      aij_ia(k) = 3
c
      k=k+1 !  'AIJ049'
      aij_lname(k) = 'ACCUMULATED NORTHWARD COMP. OF SURFACE STRESS'
      aij_units(k) = 'kg/m*s^2'
      aij_name(k) = 'TAUVS'
      aij_ia(k) = 3
c
      k=k+1 !  'AIJ050'
      aij_lname(k) = 'WATER1+WATER2+ICE1+ICE2'
      aij_units(k) = 'kg/m^2'
      aij_name(k) = 'GWTR'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ051'
      aij_lname(k) = 'ACCUMULATED SURFACE AIR SPECIFIC HUMIDITY'
      aij_units(k) = '1'
      aij_name(k) = 'QS'
      aij_ia(k) = 3
c
      k=k+1 !  'AIJ052'
      aij_lname(k) = 'MAX(0,33-1.8*DAILY MEAN ON TS IN C)'
      aij_units(k) = 'Fahrenheit*day'
      aij_name(k) = 'STRNGTS'
      aij_ia(k) = 9
c
      k=k+1 !  'AIJ053'
      aij_lname(k) = 'UNDERGROUND RUNOFF'
      aij_units(k) = 'mm/day'
      aij_name(k) = 'ARUNU'
      aij_ia(k) = 3
c
      k=k+1 !  'AIJ054'
      aij_lname(k) = '18*(DEL(TG)/DEL(TS)-1)'
      aij_units(k) = 'unknown'
      aij_name(k) = 'DTGDTS'
      aij_ia(k) = 9
c
      k=k+1 !  'AIJ055'
      aij_lname(k) = '8*P*U*Q (VERTICALLY INTEGRATED)'
      aij_units(k) = '12.5 PA*m/s'
      aij_name(k) = 'PUQ'
      aij_ia(k) = 4
c
      k=k+1 !  'AIJ056'
      aij_lname(k) = '8*P*V*Q (VERTICALLY INTEGRATED)'
      aij_units(k) = '12.5 PA*m/s'
      aij_name(k) = 'PVQ'
      aij_ia(k) = 4
c
      k=k+1 !  'AIJ057'
      aij_lname(k) = 'TGO=ODATA(1)'
      aij_units(k) = 'degC'
      aij_name(k) = 'TGO'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ058'
      aij_lname(k) = 'ACE2OI=ODATA(3)*POICE'
      aij_units(k) = 'kg/m^2'
      aij_name(k) = 'MSI2'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ059'
      aij_lname(k) = 'WIND SPEED IN TOP LAYER'
      aij_units(k) = 'm/s'
      aij_name(k) = 'WLM'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ060'
      aij_lname(k) = 'TGO12=ODATA(5)'
      aij_units(k) = 'degC'
      aij_name(k) = 'TGO2'
      aij_ia(k) = 9
c
      k=k+1 !  'AIJ061'
      aij_lname(k) = 'EVAP*POCEAN'
      aij_units(k) = 'mm/day'
      aij_name(k) = 'EVAPO'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ062'
      aij_lname(k) = 'EVAP*POICE'
      aij_units(k) = 'mm/day'
      aij_name(k) = 'EVAPI'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ063'
      aij_lname(k) = 'EVAP OVER LAND ICE'
      aij_units(k) = 'mm/day'
      aij_name(k) = 'EVAPLI'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ064'
      aij_lname(k) = 'EVAP OVER EARTH'
      aij_units(k) = 'mm/day'
      aij_name(k) = 'EVAPE'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ065'
      aij_lname(k) = 'F0DT*POCEAN, NET HEAT AT Z0'
      aij_units(k) = 'J/m^2'
      aij_name(k) = 'F0OC'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ066'
      aij_lname(k) = 'F0DT*POICE, NET HEAT AT Z0'
      aij_units(k) = 'J/m^2'
      aij_name(k) = 'F0OI'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ067'
      aij_lname(k) = 'F0DT, NET HEAT AT Z0 OVER LAND ICE'
      aij_units(k) = 'J/m^2'
      aij_name(k) = 'F0LI'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ068'
      aij_lname(k) = ' F0DT, NET HEAT AT Z0 OVER EARTH'
      aij_units(k) = 'J/m^2'
      aij_name(k) = 'F0E'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ069'
      aij_lname(k) = 'F1DT OVER LAND ICE'
      aij_units(k) = 'J/m^2'
      aij_name(k) = 'F1LI'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ070'
      aij_lname(k) = 'SNOW FALL (H2O EQUIV)'
      aij_units(k) = 'mm/day'
      aij_name(k) = 'SNWF'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ071'
      aij_lname(k) = 'SURF AIR TEMP OVER LAND ICE'
      aij_units(k) = 'degC'
      aij_name(k) = 'TSLI'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ072'
      aij_lname(k) = 'F2DT OVER LAND ICE'
      aij_units(k) = 'J/m^2'
      aij_name(k) = 'ERUN2'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ073'
      aij_lname(k) = 'SENS HEAT FLUX OVER LAND ICE'
      aij_units(k) = 'W/m^2'
      aij_name(k) = 'SHDTLI'
      aij_ia(k) = 3
c
      k=k+1 !  'AIJ074'
      aij_lname(k) = 'LATENT HEAT FLUX OVER LAND ICE'
      aij_units(k) = 'W/m^2'
      aij_name(k) = 'EVHDT'
      aij_ia(k) = 3
c
      k=k+1 !  'AIJ075'
      aij_lname(k) = 'TRHDT OVER LAND ICE'
      aij_units(k) = 'J/m^2'
      aij_name(k) = 'TRHDT'
      aij_ia(k) = 3
c
      k=k+1 !  'AIJ076'
      aij_lname(k) = 'MAX(COMPOSITE TS)'
      aij_units(k) = 'degC'
      aij_name(k) = 'TMAX'
      aij_ia(k) = 12
c
      k=k+1 !  'AIJ077'
      aij_lname(k) = 'MIN(COMPOSITE TS)'
      aij_units(k) = 'degC'
      aij_name(k) = 'TMIN'
      aij_ia(k) = 12
c
      k=k+1 !  'AIJ078'
      aij_lname(k) = 'MIN(DIURNAL MAX OF COMPOSITE TS)'
      aij_units(k) = 'degC'
      aij_name(k) = 'TMNMX'
      aij_ia(k) = 12
c
      k=k+1 !  'AIJ079'
      aij_lname(k) = 'POTENTIAL EVAPORATION'
      aij_units(k) = 'mm/day'
      aij_name(k) = 'PEVAP'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ080'
      aij_lname(k) = 'MAX TS OVER EARTH FOR CURRENT DAY'
      aij_units(k) = 'K'
      aij_name(k) = 'TMAXE'
      aij_ia(k) = 9
c
      k=k+1 !  'AIJ081'
      aij_lname(k) = 'LIQUID WATER PATH'
      aij_units(k) = 'kg/m^2'
      aij_name(k) = 'WMSUM'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ082'
      aij_lname(k) = 'ACCUMULATED SHALLOW CONVECTIVE CLOUD COVER'
      aij_units(k) = '1'
      aij_name(k) = 'PSCLD'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ083'
      aij_lname(k) = 'ACCUMULATED DEEP CONVECTIVE CLOUD COVER'
      aij_units(k) = '1'
      aij_name(k) = 'PDCLD'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ084'
      aij_lname(k) = 'ACCUMULATED DEEP CONVECTIVE CLOUD OCCURENCE'
      aij_units(k) = '1'
      aij_name(k) = 'DCNVFRQ'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ085'
      aij_lname(k) = 'ACCUMULATED SHALLOW CONVECTIVE CLOUD OCCURENCE'
      aij_units(k) = '1'
      aij_name(k) = 'SCNVFRQ'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ086'
      aij_lname(k) = 'INCIDENT MTN EAST MOMENTUM FLUX'
      aij_units(k) = 'mb m/s^2'
      aij_name(k) = 'EMTMOM'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ087'
      aij_lname(k) = 'INCIDENT MTN SOUTH MOMENTUM FLUX'
      aij_units(k) = 'mb m/s^2'
      aij_name(k) = 'SMTMOM'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ088'
      aij_lname(k) = 'EAST-WEST POTENTIAL ENTHALPY FLUX'
      aij_units(k) = 'W'
      aij_name(k) = 'FPEU'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ089'
      aij_lname(k) = 'NORTH-SOUTH POTENTIAL ENTHALPY FLUX'
      aij_units(k) = 'W'
      aij_name(k) = 'FPEV'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ090'
      aij_lname(k) = 'EAST-WEST MASS FLUX'
      aij_units(k) = 'kg/s'
      aij_name(k) = 'FMU'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ091'
      aij_lname(k) = 'NORTH-SOUTH MASS FLUX'
      aij_units(k) = 'kg/s'
      aij_name(k) = 'FMV'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ092'
      aij_lname(k) = 'EAST-WEST WATER VAPOR FLUX'
      aij_units(k) = 'kg/s'
      aij_name(k) = 'FQU'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ093'
      aij_lname(k) = 'NORTH-SOUTH WATER VAPOR FLUX'
      aij_units(k) = 'kg/s'
      aij_name(k) = 'FQV'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ094'
      aij_lname(k) = 'EAST-WEST GEOPOTENTIAL FLUX'
      aij_units(k) = 'W'
      aij_name(k) = 'FGZU'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ095'
      aij_lname(k) = 'NORTH-SOUTH GEOPOTENTIAL FLUX'
      aij_units(k) = 'W'
      aij_name(k) = 'FGZV'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ096'
      aij_lname(k) = 'unknown'
      aij_units(k) = 'unknown'
      aij_name(k) = 'ERVR'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ097'
      aij_lname(k) = 'unknown'
      aij_units(k) = 'unknown'
      aij_name(k) = 'MRVR'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ098'
      aij_lname(k) = 'unknown'
      aij_units(k) = 'unknown'
      aij_name(k) = 'SDRAG'
      aij_ia(k) = 1
c
      k=k+1 !  'AIJ099'
      aij_lname(k) = 'unknown'
      aij_units(k) = 'unknown'
      aij_name(k) = 'AIJ099'
c
      k=k+1 !  'AIJ100'
      aij_lname(k) = 'unknown'
      aij_units(k) = 'unknown'
      aij_name(k) = 'AIJ100'
c
      return
      end subroutine aij_defs

      subroutine ail_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      ail_name(k) = 'AIL01'
      ail_lname(k) = 'U (SUM FOR J=JEQ+1,JEQ,JEQ-1,JEQ-2) (PU GRID)'
      ail_units(k) = 'm/s'
c
      k=k+1
      ail_name(k) = 'AIL02'
      ail_lname(k) = 'U (SUM FOR J=JEQ+1,JEQ,JEQ-1,JEQ-2) (PU GRID)'
      ail_units(k) = 'm/s'
c
      k=k+1
      ail_name(k) = 'AIL03'
      ail_lname(k) = 'SD (SUM FOR J=JEQ,JEQ-1,JEQ-2)'
      ail_units(k) = '100 N/s'
c
      k=k+1
      ail_name(k) = 'AIL04'
      ail_lname(k) = 'TX (SUM FOR J=JEQ,JEQ-1,JEQ-2)'
      ail_units(k) = 'degC'
c
      k=k+1
      ail_name(k) = 'AIL05'
      ail_lname(k) = 'RH (SUM FOR J=JEQ,JEQ-1,JEQ-2)'
      ail_units(k) = '1'
c
      k=k+1
      ail_name(k) = 'AIL06'
      ail_lname(k) = 'DTX(MC)*P*DA (SUM FOR J=JEQ,JEQ-1,JEQ-2)'
      ail_units(k) = '100 K*N'
c
      k=k+1
      ail_name(k) = 'AIL07'
      ail_lname(k) = '(SRHR+TRHR)*DA (SUM FOR J=JEQ,JEQ-1,JEQ-2)'
      ail_units(k) = 'W'
c
      k=k+1
      ail_name(k) = 'AIL08'
      ail_lname(k) = 'unknown'
      ail_units(k) = 'unknown'
c
      k=k+1
      ail_name(k) = 'AIL09'
      ail_lname(k) = 'SD (AT LAT 50 N) (COMMENTED OUT)'
      ail_units(k) = '100 N/s'
c
      k=k+1
      ail_name(k) = 'AIL10'
      ail_lname(k) = 'TX-273.16  (AT LAT 50 N)'
      ail_units(k) = 'unknown'
c
      k=k+1
      ail_name(k) = 'AIL11'
      ail_lname(k) = 'SR+TR  (AT LAT 50 N)'
      ail_units(k) = 'unknown'
c
      k=k+1
      ail_name(k) = 'AIL12'
      ail_lname(k) = '2*U  (AT LAT 50 N)'
      ail_units(k) = 'unknown'
c
      k=k+1
      ail_name(k) = 'AIL13'
      ail_lname(k) = 'SD  (AT LAT 70 N) (COMMENTED OUT)'
      ail_units(k) = 'unknown'
c
      k=k+1
      ail_name(k) = 'AIL14'
      ail_lname(k) = 'TX-273.16  (AT LAT 70 N)  (COMMENTED OUT)'
      ail_units(k) = 'unknown'
c
      k=k+1
      ail_name(k) = 'AIL15'
      ail_lname(k) = 'SR+TR  (AT LAT 70 N)'
      ail_units(k) = 'unknown'
c
      k=k+1
      ail_name(k) = 'AIL16'
      ail_lname(k) = '2*U  (AT LAT 70 N)        (COMMENTED OUT)'
      ail_units(k) = 'unknown'
c
      return
      end subroutine ail_defs

      subroutine ajl_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      ajl_name(k) = 'AJL01'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL02'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL03'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL04'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL05'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL06'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL07'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL08'
      ajl_lname(k) = 'FMX(MC)*P'
      ajl_units(k) = '100 PA'
c
      k=k+1
      ajl_name(k) = 'AJL09'
      ajl_lname(k) = 'SRHR'
      ajl_units(k) = 'W/m^2'
c
      k=k+1
      ajl_name(k) = 'AJL10'
      ajl_lname(k) = 'TRHR'
      ajl_units(k) = 'W/m^2'
c
      k=k+1
      ajl_name(k) = 'AJL11'
      ajl_lname(k) = 'DTX(SS)*P'
      ajl_units(k) = '100 K*PA'
c
      k=k+1
      ajl_name(k) = 'AJL12'
      ajl_lname(k) = 'DT(DC)*P'
      ajl_units(k) = '100 K*PA'
c
      k=k+1
      ajl_name(k) = 'AJL13'
      ajl_lname(k) = 'DT(MC)*P  DRY HEATING'
      ajl_units(k) = '100 PA*K'
c
      k=k+1
      ajl_name(k) = 'AJL14'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL15'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL16'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL17'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL18'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL19'
      ajl_lname(k) = 'PCLD*P (TOTAL)'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL20'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL21'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL22'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL23'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL24'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL25'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL26'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL27'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL28'
      ajl_lname(k) = 'PCLD*P (SS)'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL29'
      ajl_lname(k) = 'PCLD*P (MC)'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL30'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL31'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL32'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL33'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL34'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL35'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL36'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL37'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL38'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL39'
      ajl_lname(k) = 'DU(MC)*P (UV GRID)'
      ajl_units(k) = '100 N/m/s'
c
      k=k+1
      ajl_name(k) = 'AJL40'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL41'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL42'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL43'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL44'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL45'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL46'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL47'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL48'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL49'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL50'
      ajl_lname(k) = 'DT(MC)*P  CHANGE OF PHASE'
      ajl_units(k) = '100 PA*K'
c
      k=k+1
      ajl_name(k) = 'AJL51'
      ajl_lname(k) = 'CLHE*DQ(MC BEFORE COND)*P'
      ajl_units(k) = '100 PA*K'
c
      k=k+1
      ajl_name(k) = 'AJL52'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL53'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL54'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL55'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL56'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      k=k+1
      ajl_name(k) = 'AJL57'
      ajl_lname(k) = 'unknown'
      ajl_units(k) = 'unknown'
c
      return
      end subroutine ajl_defs

      subroutine asjl_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      asjl_name(k) = 'ASJL01'
      asjl_lname(k) = 'TX'
      asjl_units(k) = 'degC'
c
      k=k+1
      asjl_name(k) = 'ASJL02'
      asjl_lname(k) = 'PHI'
      asjl_units(k) = 'm^2/s^2'
c
      k=k+1
      asjl_name(k) = 'ASJL03'
      asjl_lname(k) = 'SRHR'
      asjl_units(k) = 'W/m^2'
c
      k=k+1
      asjl_name(k) = 'ASJL04'
      asjl_lname(k) = 'TRHR'
      asjl_units(k) = 'W/m^2'
c
      return
      end subroutine asjl_defs

      subroutine ajk_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      ajk_name(k) = 'AJK01'
      ajk_lname(k) = 'DP=(PDN-PM(K+1);PDN=MAX(PM(K+1),PS)'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK02'
      ajk_lname(k) = 'DP4  (PDN-PM(K+1))   (UV GRID)'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK03'
      ajk_lname(k) = '(TX-273.16)*DP'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK04'
      ajk_lname(k) = 'PHI*DP'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK05'
      ajk_lname(k) = 'Q*DP'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK06'
      ajk_lname(k) = 'TH*DP'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK07'
      ajk_lname(k) = 'RH*DP'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK08'
      ajk_lname(k) = 'U*DP4  (UV GRID)'
      ajk_units(k) = '100 PA*m/s'
c
      k=k+1
      ajk_name(k) = 'AJK09'
      ajk_lname(k) = 'V*DP4  (UV GRID)'
      ajk_units(k) = '100 PA*m/s'
c
      k=k+1
      ajk_name(k) = 'AJK10'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK11'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK12'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK13'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK14'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK15'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK16'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK17'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK18'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK19'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK20'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK21'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK22'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK23'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK24'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK25'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK26'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK27'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK28'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK29'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK30'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK31'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK32'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK33'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK34'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK35'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK36'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK37'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK38'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK39'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK40'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK41'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK42'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK43'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK44'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK45'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK46'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK47'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK48'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK49'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK50'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      k=k+1
      ajk_name(k) = 'AJK51'
      ajk_lname(k) = 'unknown'
      ajk_units(k) = 'unknown'
c
      return
      end subroutine ajk_defs

      subroutine aijg_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      aijg_name(k) = 'AIJG01'
      aijg_lname(k) = 'LAYER 1 REL SAT OF BARE AMD VEG. SOIL'
      aijg_units(k) = '%'
c
      k=k+1
      aijg_name(k) = 'AIJG02'
      aijg_lname(k) = 'LAYER 2 REL SAT OF BARE AMD VEG. SOIL'
      aijg_units(k) = '%'
c
      k=k+1
      aijg_name(k) = 'AIJG03'
      aijg_lname(k) = 'LAYER 3 REL SAT OF BARE AMD VEG. SOIL'
      aijg_units(k) = '%'
c
      k=k+1
      aijg_name(k) = 'AIJG04'
      aijg_lname(k) = 'LAYER 4 REL SAT OF BARE AMD VEG. SOIL'
      aijg_units(k) = '%'
c
      k=k+1
      aijg_name(k) = 'AIJG05'
      aijg_lname(k) = 'BETA (BARE SOIL)'
      aijg_units(k) = '%'
c
      k=k+1
      aijg_name(k) = 'AIJG06'
      aijg_lname(k) = 'BETA (PENMAN)'
      aijg_units(k) = '%'
c
      k=k+1
      aijg_name(k) = 'AIJG07'
      aijg_lname(k) = 'CANOPY  REL SATURATION'
      aijg_units(k) = '%'
c
      k=k+1
      aijg_name(k) = 'AIJG08'
      aijg_lname(k) = 'LAYER 1 REL SATURATION OF VEG. SOIL'
      aijg_units(k) = '%'
c
      k=k+1
      aijg_name(k) = 'AIJG09'
      aijg_lname(k) = 'LAYER 2 REL SATURATION OF VEG. SOIL'
      aijg_units(k) = '%'
c
      k=k+1
      aijg_name(k) = 'AIJG10'
      aijg_lname(k) = 'LAYER 3 REL SATURATION OF VEG. SOIL'
      aijg_units(k) = '%'
c
      k=k+1
      aijg_name(k) = 'AIJG11'
      aijg_lname(k) = 'BETA (BARE SOIL & VEGETATION)'
      aijg_units(k) = '%'
c
      k=k+1
      aijg_name(k) = 'AIJG12'
      aijg_lname(k) = 'CONDUCTANCE OF ATMOSPHERE'
      aijg_units(k) = '.01 m/s'
c
      k=k+1
      aijg_name(k) = 'AIJG13'
      aijg_lname(k) = 'CONDUCTANCE OF CANOPY'
      aijg_units(k) = '.01 m/s'
c
      k=k+1
      aijg_name(k) = 'AIJG14'
      aijg_lname(k) = 'PENMAN POTENTIAL EVAPORATION'
      aijg_units(k) = 'mm/day'
c
      k=k+1
      aijg_name(k) = 'AIJG15'
      aijg_lname(k) = 'TEMP OF LAYER 1 BARE SOIL AND SNOW'
      aijg_units(k) = 'degC'
c
      k=k+1
      aijg_name(k) = 'AIJG16'
      aijg_lname(k) = 'TEMP OF SOIL LAYER 2 - BARE SOIL'
      aijg_units(k) = 'degC'
c
      k=k+1
      aijg_name(k) = 'AIJG17'
      aijg_lname(k) = 'TEMP OF SOIL LAYER 3 - BARE SOIL'
      aijg_units(k) = 'degC'
c
      k=k+1
      aijg_name(k) = 'AIJG18'
      aijg_lname(k) = 'BARE SOIL EVAPORATION'
      aijg_units(k) = 'mm/day'
c
      k=k+1
      aijg_name(k) = 'AIJG19'
      aijg_lname(k) = 'DRY CANOPY EVAPORATION'
      aijg_units(k) = 'mm/day'
c
      k=k+1
      aijg_name(k) = 'AIJG20'
      aijg_lname(k) = 'WET CANOPY EVAPORATION'
      aijg_units(k) = 'mm/day'
c
      k=k+1
      aijg_name(k) = 'AIJG21'
      aijg_lname(k) = 'TEMP OF CANOPY AND SNOW'
      aijg_units(k) = 'degC'
c
      k=k+1
      aijg_name(k) = 'AIJG22'
      aijg_lname(k) = 'TEMP OF SOIL LAYER 1 - VEGETATED SOIL'
      aijg_units(k) = 'degC'
c
      k=k+1
      aijg_name(k) = 'AIJG23'
      aijg_lname(k) = 'TEMP OF SOIL LAYER 2 - VEGETATED SOIL'
      aijg_units(k) = 'degC'
c
      k=k+1
      aijg_name(k) = 'AIJG24'
      aijg_lname(k) = 'TEMP OF SOIL LAYER 3 - VEGETATED SOIL'
      aijg_units(k) = 'degC'
c
      k=k+1
      aijg_name(k) = 'AIJG25'
      aijg_lname(k) = 'AVERAGE WATER TABLE'
      aijg_units(k) = 'm'
c
      k=k+1
      aijg_name(k) = 'AIJG26'
      aijg_lname(k) = 'BETAV,OVER VEGETATION'
      aijg_units(k) = '%'
c
      k=k+1
      aijg_name(k) = 'AIJG27'
      aijg_lname(k) = 'BETAT,TRANSPIRATION'
      aijg_units(k) = '%'
c
      k=k+1
      aijg_name(k) = 'AIJG28'
      aijg_lname(k) = 'SNOW DEPTH OVER BARE SOIL'
      aijg_units(k) = 'MM H2O'
c
      k=k+1
      aijg_name(k) = 'AIJG29'
      aijg_lname(k) = 'SNOW DEPTH OVER VEG SOIL'
      aijg_units(k) = 'MM H2O'
c
      return
      end subroutine aijg_defs

      subroutine consrv_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      consrv_name(k) = 'CONSRV01'
      consrv_lname(k) = 'INSTANTANEOUS AM'
      consrv_units(k) = '10**9 J*s/m^2'
c
      k=k+1
      consrv_name(k) = 'CONSRV02'
      consrv_lname(k) = 'CHANGE OF AM BY ADVECTION'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV03'
      consrv_lname(k) = 'CHANGE OF AM BY CORIOLIS FORCE'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV04'
      consrv_lname(k) = 'CHANGE OF AM BY ADVEC + COR'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV05'
      consrv_lname(k) = 'CHANGE OF AM BY PRESSURE GRAD'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV06'
      consrv_lname(k) = 'CHANGE OF AM BY DYNAMICS'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV07'
      consrv_lname(k) = 'CHANGE OF AM BY SURFACE FRIC'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV08'
      consrv_lname(k) = 'CHANGE OF AM BY STRATOS DRAG'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV09'
      consrv_lname(k) = 'CHANGE OF AM BY FILTER'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV10'
      consrv_lname(k) = 'CHANGE OF AM BY DAILY RESTOR'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV11'
      consrv_lname(k) = 'SUM OF AM CHANGES'
      consrv_units(k) = '10**2 J/m^2'
c
      k=k+1
      consrv_name(k) = 'CONSRV12'
      consrv_lname(k) = 'INSTANTANEOUS KE'
      consrv_units(k) = '10**3 J/m^2'
c
      k=k+1
      consrv_name(k) = 'CONSRV13'
      consrv_lname(k) = 'CHANGE OF KE BY ADVECTION'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV14'
      consrv_lname(k) = 'CHANGE OF KE BY CORIOLIS FORCE'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV15'
      consrv_lname(k) = 'CHANGE OF KE BY ADVEC + COR'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV16'
      consrv_lname(k) = 'CHANGE OF KE BY PRESSURE GRAD'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV17'
      consrv_lname(k) = 'CHANGE OF KE BY DYNAMICS'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV18'
      consrv_lname(k) = 'CHANGE OF KE BY MOIST CONVEC'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV19'
      consrv_lname(k) = 'CHANGE OF KE BY SURF + DRY CONV'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV20'
      consrv_lname(k) = 'CHANGE OF KE BY STRATOS DRAG'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV21'
      consrv_lname(k) = 'CHANGE OF KE BY FILTER'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV22'
      consrv_lname(k) = 'CHANGE OF KE BY DAILY RESTOR'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV23'
      consrv_lname(k) = 'SUM OF KE CHANGES'
      consrv_units(k) = '10**-3 W/m^2'
c
      k=k+1
      consrv_name(k) = 'CONSRV24'
      consrv_lname(k) = 'INSTANTANEOUS MASS'
      consrv_units(k) = 'kg/m^2'
c
      k=k+1
      consrv_name(k) = 'CONSRV25'
      consrv_lname(k) = 'CHANGE OF MASS BY DYNAMICS'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV26'
      consrv_lname(k) = 'CHANGE OF MASS BY FILTER'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV27'
      consrv_lname(k) = 'CHANGE OF MASS BY DAILY RESTOR'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV28'
      consrv_lname(k) = 'SUM MASS CHANGES'
      consrv_units(k) = '10**-8 kg/m^2/s'
c
      k=k+1
      consrv_name(k) = 'CONSRV29'
      consrv_lname(k) = 'INSTANTANE TPE'
      consrv_units(k) = '10**5 J/m^2'
c
      k=k+1
      consrv_name(k) = 'CONSRV30'
      consrv_lname(k) = 'CHANGE OF TPE BY DYNAMICS'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV31'
      consrv_lname(k) = 'CHANGE OF TPE BY CONDENSATION'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV32'
      consrv_lname(k) = 'CHANGE OF TPE BY RADIATION'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV33'
      consrv_lname(k) = 'CHANGE OF TPE BY SURFACE INTER'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV34'
      consrv_lname(k) = 'CHANGE OF TPE BY FILTER'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV35'
      consrv_lname(k) = 'CHANGE OF TPE BY DAILY RESTOR'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV36'
      consrv_lname(k) = 'SUM OF TPE CHANGES'
      consrv_units(k) = '10**-2 W/m^2'
c
      k=k+1
      consrv_name(k) = 'CONSRV37'
      consrv_lname(k) = 'INSTANT WATER'
      consrv_units(k) = '10**-2 kg/m^2'
c
      k=k+1
      consrv_name(k) = 'CONSRV38'
      consrv_lname(k) = 'CHANGE OF WATER BY DYNAMICS'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV39'
      consrv_lname(k) = 'CHANGE OF WATER BY CONDENSATION'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV40'
      consrv_lname(k) = 'CHANGE OF WATER BY SURFACE EVAP'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV41'
      consrv_lname(k) = 'CHANGE OF WATER BY DAILY RESTOR'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV42'
      consrv_lname(k) = 'SUM WATER CHANGES'
      consrv_units(k) = '10**-8 kg/m^2/s'
c
      k=k+1
      consrv_name(k) = 'CONSRV43'
      consrv_lname(k) = 'unknown'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV44'
      consrv_lname(k) = 'unknown'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV45'
      consrv_lname(k) = 'unknown'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV46'
      consrv_lname(k) = 'unknown'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV47'
      consrv_lname(k) = 'unknown'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV48'
      consrv_lname(k) = 'unknown'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV49'
      consrv_lname(k) = 'unknown'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV50'
      consrv_lname(k) = 'unknown'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV51'
      consrv_lname(k) = 'unknown'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV52'
      consrv_lname(k) = 'unknown'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV53'
      consrv_lname(k) = 'unknown'
      consrv_units(k) = 'unknown'
c
      k=k+1
      consrv_name(k) = 'CONSRV54'
      consrv_lname(k) = 'unknown'
      consrv_units(k) = 'unknown'
c
      return
      end subroutine consrv_defs

      subroutine aijk_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      aijk_name(k) = 'UDPB'
      aijk_lname(k) = 'u-wind x delta p, b-grid'
      aijk_units(k) = '100 PA*m/s'
c
      k=k+1
      aijk_name(k) = 'VDPB'
      aijk_lname(k) = 'v-wind x delta p, b-grid'
      aijk_units(k) = '100 PA*m/s'
c
      k=k+1
      aijk_name(k) = 'DSEDPBx4'
      aijk_lname(k) = 'dry stat. energy x delta p x 4, b-grid'
      aijk_units(k) = '100 N/s^2'
c
      k=k+1
      aijk_name(k) = 'DPB'
      aijk_lname(k) = 'delta p, b-grid'
      aijk_units(k) = '100 PA'
c
      k=k+1
      aijk_name(k) = 'TDPBx4'
      aijk_lname(k) = 'temperature x delta p x 4, b-grid'
      aijk_units(k) = '100 K*PA'
c
      k=k+1
      aijk_name(k) = 'QDPBx4'
      aijk_lname(k) = 'spec. humidity x delta p x 4, b-grid'
      aijk_units(k) = '100 PA'
c
      return
      end subroutine aijk_defs

      subroutine aijl_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      aijl_name(k) = 'AIJL1'
      aijl_lname(k) = 'unknown'
      aijl_units(k) = 'unknown'
c
      k=k+1
      aijl_name(k) = 'AIJL2'
      aijl_lname(k) = 'unknown'
      aijl_units(k) = 'unknown'
c
      k=k+1
      aijl_name(k) = 'AIJL3'
      aijl_lname(k) = 'unknown'
      aijl_units(k) = 'unknown'
c
      k=k+1
      aijl_name(k) = 'AIJL4'
      aijl_lname(k) = 'unknown'
      aijl_units(k) = 'unknown'
c
      k=k+1
      aijl_name(k) = 'AIJL5'
      aijl_lname(k) = 'unknown'
      aijl_units(k) = 'unknown'
c
      return
      end subroutine aijl_defs

      subroutine wave_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      wave_name(k) = 'U850_EQ'
      wave_lname(k) = 'unknown'
      wave_units(k) = 'unknown'
c
      k=k+1
      wave_name(k) = 'V850_EQ'
      wave_lname(k) = 'unknown'
      wave_units(k) = 'unknown'
c
      k=k+1
      wave_name(k) = 'U300_EQ'
      wave_lname(k) = 'unknown'
      wave_units(k) = 'unknown'
c
      k=k+1
      wave_name(k) = 'V300_EQ'
      wave_lname(k) = 'unknown'
      wave_units(k) = 'unknown'
c
      k=k+1
      wave_name(k) = 'U050_EQ'
      wave_lname(k) = 'unknown'
      wave_units(k) = 'unknown'
c
      k=k+1
      wave_name(k) = 'V050_EQ'
      wave_lname(k) = 'unknown'
      wave_units(k) = 'unknown'
c
      k=k+1
      wave_name(k) = 'Z922_50N'
      wave_lname(k) = 'unknown'
      wave_units(k) = 'unknown'
c
      k=k+1
      wave_name(k) = 'Z700_50N'
      wave_lname(k) = 'unknown'
      wave_units(k) = 'unknown'
c
      k=k+1
      wave_name(k) = 'Z500_50N'
      wave_lname(k) = 'unknown'
      wave_units(k) = 'unknown'
c
      k=k+1
      wave_name(k) = 'Z300_50N'
      wave_lname(k) = 'unknown'
      wave_units(k) = 'unknown'
c
      k=k+1
      wave_name(k) = 'Z100_50N'
      wave_lname(k) = 'unknown'
      wave_units(k) = 'unknown'
c
      k=k+1
      wave_name(k) = 'Z010_50N'
      wave_lname(k) = 'unknown'
      wave_units(k) = 'unknown'
c
      return
      end subroutine wave_defs

      subroutine ajlsp_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      ajlsp_name(k) = 'AJLSP01'
      ajlsp_lname(k) = 'unknown'
      ajlsp_units(k) = 'unknown'
c
      k=k+1
      ajlsp_name(k) = 'AJLSP02'
      ajlsp_lname(k) = 'unknown'
      ajlsp_units(k) = 'unknown'
c
      k=k+1
      ajlsp_name(k) = 'AJLSP03'
      ajlsp_lname(k) = 'unknown'
      ajlsp_units(k) = 'unknown'
c
      return
      end subroutine ajlsp_defs

      subroutine apj_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      apj_name(k) = 'APJ1'
      apj_lname(k) = 'P'
      apj_units(k) = '100 PA'
c
      k=k+1
      apj_name(k) = 'APJ2'
      apj_lname(k) = '4*P4I (UV GRID)'
      apj_units(k) = '100 PA'
c
      k=k+1
      apj_name(k) = 'APJ3'
      apj_lname(k) = 'unknown'
      apj_units(k) = 'unknown'
c
      return
      end subroutine apj_defs

      subroutine tsf_defs
      use DAGCOM
      implicit none
      integer :: k
c
      k=0
c
      k=k+1
      tsf_name(k) = 'TSFREZ1'
      tsf_lname(k) = 'FIRST DAY OF GROWING SEASON'
      tsf_units(k) = 'JULIAN DAY'
c
      k=k+1
      tsf_name(k) = 'TSFREZ2'
      tsf_lname(k) = 'LAST DAY OF GROWING SEASON'
      tsf_units(k) = 'JULIAN DAY'
c
      return
      end subroutine tsf_defs
