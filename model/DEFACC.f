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
      call wave_defs
      call jk_defs
      call ijk_defs
      call diurn_defs
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
      call set_iparm('LS1',LS1)
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
      call set_iparm('MDYN',MDYN)
      call set_iparm('MCNDS',MCNDS)
      call set_iparm('MRAD',MRAD)
      call set_iparm('MSURF',MSURF)
      call set_iparm('MDIAG',MDIAG)
      call set_iparm('MELSE',MELSE)
      call set_iparm('MODRD',MODRD)
      call set_iparm('MODD5K',MODD5K)
      call set_iparm('MODD5S',MODD5S)
      call set_iparm('IYEAR1',IYEAR1)
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
      use CLOUDS, only : u00wtr,u00ice
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
!!    call set_dparm('SKIPSE',SKIPSE)
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
      use CONSTANT, only : grav,sday,shw,rgas,omega,bygrav,gamd
      use MODEL_COM, only : jm,lm,ls1,dtsrc,fim,sige, DO_GWDRAG,
     *  iDO_GWDRAG
      use DAGCOM
      use DAGPCOM, only : p1000k
      implicit none
      integer :: k
c
      do k=1,kaj
         write(name_j(k),'(a2,i3.3)') 'AJ',k
         lname_j(k) = 'unused'
         units_j(k) = 'unused'
         stitle_j(k)= 'no output'
         scale_j(k) = 1.
         ia_j(k)    = 1.
      enddo
c
      k=0
c  AJ01
      k=k+1
      J_SRINCP0= k ! SRINCP0 (W/M**2)                              2 RD
      name_j(k) = 'inc_sw'
      lname_j(k) = 'SOLAR RADIATION INCIDENT ON PLANET'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' INC SW(WT/M**2)'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ02
      k=k+1
      J_SRNFP0=  k ! SRNFP0 (W/M**2)                               2 RD
      name_j(k) = 'sw_abs_p0'
      lname_j(k) = ' SOLAR RADIATION ABSORBED BY PLANET'
      units_j(k) = 'W/m^2'
      stitle_j(k)= '0SW ABS BELOW P0'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ03
      k=k+1
      J_SRNFP1=  k ! SRNFP1 (W/M**2)                               2 RD
      name_j(k) = 'sw_abs_p1'
      lname_j(k) = 'SOLAR RADIATION ABSORBED BELOW PTOP'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' SW ABS BELOW P1'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ04
      k=k+1
      J_SRABS =  k ! SRABSATM=AJ(SRNFP0)-AJ(SRNFG) (W/M**2)        2 D1
      name_j(k) = 'sw_abs_atm'
      lname_j(k) = 'SOLAR RADIATION ABSORBED BY ATMOSPHERE'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' SW ABS BY ATMOS'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ05
      k=k+1
      J_SRINCG=  k ! SRINCG (W/M**2)                               2 RD
      name_j(k) = 'sw_inc_z0'
      lname_j(k) = 'SOLAR RADIATION INCIDENT ON GROUND'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' SW INC ON Z0   '
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ06
      k=k+1
      J_SRNFG =  k ! SRNFG (W/M**2)                                2 RD
      name_j(k) = 'sw_abs_z0'
      lname_j(k) = 'SOLAR RADIATION ABSORBED BY GROUND'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' SW ABS AT Z0   '
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ07
      k=k+1
      J_TRNFP0=  k ! TRNFP0=AJ(ALBVIS)+A2BYA1*AJ(TRHDT)/DTS(W/M**2)2 D1
      name_j(k) = 'net_lw_p0'
      lname_j(k) = 'THERMAL RADIATION EMITTED BY PLANET'
      units_j(k) = 'W/m^2'
      stitle_j(k)= '0NET LW AT P0   '
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ08
      k=k+1
      J_TRNFP1=  k ! TRNFP1=AJ(ALBNIR)+A2BYA1*AJ(TRHDT)/DTS(W/M**2)2 D1
      name_j(k) = 'net_lw_p1'
      lname_j(k) = 'NET THERMAL RADIATION AT PTOP'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET LW AT P1   '
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ09
      k=k+1
      J_TRHDT =  k ! TRHDT (J/M**2)                                1 SF
      name_j(k) = 'net_lw_z0'
      lname_j(k) = 'NET THERMAL RADIATION AT GROUND'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET LW AT Z0   '
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c  AJ10
      k=k+1
      J_RNFP0 = k ! RNFP0=AJ(SRNFP0)+AJ(TRNFP0) (W/M**2)          2 D1
      name_j(k) = 'net_rad_p0'
      lname_j(k) = 'NET RADIATION OF PLANET'
      units_j(k) = 'W/m^2'
      stitle_j(k)= '0NET RAD AT P0  '
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ11
      k=k+1
      J_RNFP1 = k ! RNFP1=AJ(SRNFP1)+AJ(TRNFP1) (W/M**2)          2 D1
      name_j(k) = 'net_rad_p1'
      lname_j(k) = 'NET RADIATION BELOW PTOP'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET RAD AT P1  '
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ12
      k=k+1
      J_RHDT  = k ! RHDT=A1BYA2*AJ(SRNFG)*DTS+AJ(TRHDT)(J/M^2)    1 D1
      name_j(k) = 'net_rad_z0'
      lname_j(k) = 'NET RADIATION ABSORBED BY GROUND'
      units_j(k) = 'J/m^2'
      stitle_j(k)= ' NET RAD AT Z0  '
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c  AJ13
      k=k+1
      J_SHDT  = k ! SHEATDT (J/M**2)                              1 SF
      name_j(k) = 'snsht_flx'
      lname_j(k) = 'SENSIBLE HEAT FLUX INTO THE GROUND'
      units_j(k) = 'W/m^2'
      stitle_j(k)= '0SENSBL HEAT FLX'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c  AJ14
      k=k+1
      J_EVHDT = k ! EVHDT (J/M**2)                                1 SF
      name_j(k) = 'evht_flx'
      lname_j(k) = 'LATENT HEAT FLUX INTO THE GROUND'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' EVAPOR HEAT FLX'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c  AJ15
      k=k+1
      J_F2DT  = k ! F2DT (J/M**2)                                 1 GD
      name_j(k) = 'ht_cond_z1z2'
      lname_j(k) = 'CONDUCTION AT BOTTOM OF GROUND LAYER 2'
      units_j(k) = 'W/m^2'
      stitle_j(k)= '0CONDC AT -Z1-Z2'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c  AJ16
      k=k+1
      J_HZ1   = k ! HEATZ1=AJ(EDIFS)+AJ(F1DT)                     1 D1
      name_j(k) = 'net_ht_z1'
      lname_j(k) = 'NET HEAT FLUX BETWEEN GROUND LAYERS'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET HEAT AT -Z1'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c  AJ17
      k=k+1
      J_TG2   = k ! TG2 (K-TF)                                    1 GD
      name_j(k) = 'tg2'
      lname_j(k) = 'TEMPERATURE OF GROUND LAYER 2'
      units_j(k) = '.1 C'
      stitle_j(k)= ' TG2 (.1 C)     '
      scale_j(k) = 10.
      ia_j(k) = ia_src
c  AJ18
      k=k+1
      J_TG1   = k ! TG1 (K-TF)                                    1 GD
      name_j(k) = 'tg1'
      lname_j(k) = 'TEMPERATURE OF GROUND LAYER 1'
      units_j(k) = '.1 C'
      stitle_j(k)= '1TG1 (.1 C)     '
      scale_j(k) = 10.
      ia_j(k) = ia_src
c  AJ19
      k=k+1
      J_EVAP  = k ! EVAP (KG/M**2)                                1 GD
      name_j(k) = 'evap'
      lname_j(k) = 'EVAPORATION'
      units_j(k) = 'mm/day'
      stitle_j(k)= ' EVAPOR (MM/DAY)'
      scale_j(k) = SDAY/DTSRC
      ia_j(k) = ia_src
c  AJ20
      k=k+1
      J_PRCP  = k ! PRCP=AJ(PRCPSS)+AJ(PRCPMC) (100 PA)           1 D1
      name_j(k) = 'prec'
      lname_j(k) = 'PRECIPITATION'
      units_j(k) = 'mm/day'
      stitle_j(k)= ' PRECIP (MM/DAY)'
      scale_j(k) = 100.*SDAY/(DTsrc*GRAV)
      ia_j(k) = ia_src
c  AJ21
      k=k+1
      J_TX    = k ! TX (K-TF)  (INTEGRAL OVER ATMOSPHERE OF)      4 DA
      name_j(k) = 'tair'
      lname_j(k) = 'AIR TEMPERATURE'
      units_j(k) = '.1 C'
      stitle_j(k)= ' T AIR (.1 C)   '
      scale_j(k) = 10.
      ia_j(k) = ia_dga
c  AJ22
      k=k+1
      J_TX1   = k ! TX1 (K-TF)                                    4 DA
      name_j(k) = 't1'
      lname_j(k) = 'TEMPERATURE OF AIR LAYER 1'
      units_j(k) = '.1 C'
      stitle_j(k)= ' T1 (.1 C)      '
      scale_j(k) = 10.
      ia_j(k) = ia_dga
c  AJ23
      k=k+1
      J_TSRF  = k ! TS (K-TF)                                     3 SF
      name_j(k) = 'tsurf'
      lname_j(k) = 'SURFACE AIR TEMPERATURE'
      units_j(k) = '.1 C'
      stitle_j(k)= '0T SURF (.1 C)  '
      scale_j(k) = 10.
      ia_j(k) = ia_srf
c  AJ24
      k=k+1
      J_DTSGST= k ! DTH/DPHI  (STRATOSPHERE)                      4 DA
      name_j(k) = 'sstab_strat'
      lname_j(k) = 'STRATOSPHERIC STATIC STABILITY'
      units_j(k) = 'C/km'
      stitle_j(k)= '0STAT STB(STRAT)'
      scale_j(k) = 1.D3*GRAV*P1000K
      ia_j(k) = ia_dga
c  AJ25
      k=k+1
      J_DTDGTR= k ! DTH/DPHI  (TROPOSPHERE)                       4 DA
      name_j(k) = 'sstab_trop'
      lname_j(k) = 'TROPOSPHERIC STATIC STABILITY'
      units_j(k) = 'C/km'
      stitle_j(k)= ' STAT STB(TROPO)'
      scale_j(k) = 1.D3*GRAV*P1000K
      ia_j(k) = ia_dga
c  AJ26
      k=k+1
      J_RICST = k ! .0625*DTH*DLNP/(DU*DU+DV*DV)  (STRATOSPHERE)  4 DA
      name_j(k) = 'rich_num_strat'
      lname_j(k) = 'STRATOSPHERIC RICHARDSON NUMBER'
      units_j(k) = '1'
      stitle_j(k)= '0RICH NUM(STRAT)'
      scale_j(k) = 16.*RGAS
      ia_j(k) = ia_dga
c  AJ27
      k=k+1
      J_RICTR = k ! .0625*DTH*DLNP/(DU*DU+DV*DV)  (TROPOSPHERE)   4 DA
      name_j(k) = 'rich_num_trop'
      lname_j(k) = 'TROPOSPHERIC RICHARDSON NUMBER'
      units_j(k) = '1'
      stitle_j(k)= ' RICH NUM(TROPO)'
      scale_j(k) = 16.*RGAS
      ia_j(k) = ia_dga
c  AJ28
      k=k+1
      J_ROSST = k ! 4*UMAX/(DX*SINJ)  (STRATOSPHERE)              4 DA
      name_j(k) = 'ross_num_strat'
      lname_j(k) = 'STRATOSPHERIC ROSSBY NUMBER'
      units_j(k) = '1'
      stitle_j(k)= ' ROSS NUM(STRAT)'
      scale_j(k) = .5/(2.*OMEGA*FIM)
      ia_j(k) = ia_dga
c  AJ29
      k=k+1
      J_ROSTR = k ! 4*UMAX/(DX*SINJ)  (TROPOSPHERE)               4 DA
      name_j(k) = 'ross_num_trop'
      lname_j(k) = 'TROPOSPHERIC ROSSBY NUMBER'
      units_j(k) = '1'
      stitle_j(k)= ' ROSS NUM(TROPO)'
      scale_j(k) = .5/(2.*OMEGA*FIM)
      ia_j(k) = ia_dga
c  AJ30
      k=k+1
      J_RSI   = k ! RSI (1)                                       1 GD
      name_j(k) = 'ocn_lak_ice_frac'
      lname_j(k) = 'OCEAN/LAKE ICE COVER'
      units_j(k) = '%'
      stitle_j(k)= ' OC/LK ICE COVER'
      scale_j(k) = 100.
      ia_j(k) = ia_src
c  AJ31
      k=k+1
      J_RSNOW = k ! PSNOW (1)                                     4 DA
      name_j(k) = 'snow_cover'
      lname_j(k) = 'SNOW COVER'
      units_j(k) = '%'
      stitle_j(k)= '0SNOW COVER     '
      scale_j(k) = 100.
      ia_j(k) = ia_src
c  AJ32
      k=k+1
      J_SWCOR = k ! SW CORRECTION   (obsolete)                    2 RD
      name_j(k) = 'SWCOR'
      lname_j(k) = 'SW CORRECTION'
      units_j(k) = 'W/M^2'
      stitle_j(k)= ' SW CORRECTION  '
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ33
      k=k+1
      J_OHT   = k ! OCEAN TRANSPORT                               1 GD
      name_j(k) = 'ocn_ht_trans'
      lname_j(k) = 'CONVERGED OCEAN HEAT TRANSPORT'
      units_j(k) = 'W/m^2'
      stitle_j(k)= '0OCEAN TRANSPORT'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c  AJ34    ! this includes OCEAN TEMPERATURE AT BASE MIXED LAYER
      k=k+1
      J_TG3  = k ! TEMPERATURE GROUND LAYER 3                     1 GD
      name_j(k) = 'tg3'
      lname_j(k) = 'TEMPERATURE GROUND LAYER 3'
      units_j(k) = '.1 C'
      stitle_j(k)= ' TG3 (.1 C)     '
      scale_j(k) = 10.
      ia_j(k) = ia_src
c  AJ35
      k=k+1
      J_DTDJS = k ! T(J+1)-T(J-1)  (SUM OVER STRATOSPHERE OF)     4 DA
      name_j(k) = 'dtdlat_strat'
      lname_j(k) = 'STRATO TEMP CHANGE PER DEGREE LATITUDE'
      units_j(k) = 'deg C/deg lat'
      stitle_j(k)= '0DT/DLAT(STRAT) '
      scale_j(k) = .5D2*(JM-1.)/((SIGE(LS1)-SIGE(LSTR+1)+1d-12)*180.)
      ia_j(k) = ia_dga
c  AJ36
      k=k+1
      J_DTDJT = k ! T(J+1)-T(J-1)  (SUM OVER TROPOSPHERE OF)      4 DA
      name_j(k) = 'dtdlat_trop'
      lname_j(k) = 'TROPO TEMP CHANGE PER DEGREE LATITUDE'
      units_j(k) = 'deg C/deg lat'
      stitle_j(k)= ' DT/DLAT(TROPO) '
      scale_j(k) = .5d2*(JM-1.)/((SIGE(1)-SIGE(LS1))*180.)
      ia_j(k) = ia_dga
c  AJ37
      k=k+1
      J_LSTR  = k ! SQRT(DTH/DLNP)/SINJ  (STRATOSPHERE)           4 DA
      name_j(k) = 'ross_radius_strat'
      lname_j(k) = 'ROSSBY RADIUS IN THE STRATOSPHERE'
      units_j(k) = '10**5 m'
      stitle_j(k)= ' L(STRAT)(10**5)'
      scale_j(k) = 1d-5*SQRT(RGAS)/(2.*OMEGA)
      ia_j(k) = ia_dga
c  AJ38
      k=k+1
      J_LTRO  = k ! SQRT(DTH/DLNP)/SINJ  (TROPOSPHERE)            4 DA
       name_j(k) = 'ross_radius_trop'
      lname_j(k) = 'ROSSBY RADIUS IN THE TROPOSPHERE'
      units_j(k) = '10**5 m'
      stitle_j(k)=  ' L(TROP) (10**5)'
      scale_j(k) = 1d-5*SQRT(RGAS)/(2.*OMEGA)
      ia_j(k) = ia_dga
c  AJ39
      k=k+1
      J_EPRCP = k ! ENERGP (J/M**2)                               1 CN
      name_j(k) = 'prec_ht_flx'
      lname_j(k) = 'PRECIPITATION HEAT FLUX INTO THE GROUND'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' PRECIP HEAT FLX'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c  AJ40
      k=k+1
      J_ERUN1 = k ! ERUN1 (J/M**2)                                1 GP
      name_j(k) = 'ht_runoff_z0'
      lname_j(k) = 'HEAT RUNOFF Z0'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' HEAT RUNOFF Z0 '
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c  AJ41
      k=k+1
      J_EDIFS = k ! EDIFS (J/M**2)                                1 GP
      name_j(k) = 'ht_wtr_difs_z1'
      lname_j(k) = 'HEAT DIFFUSION OF WATER OR ICE AT -Z1'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' HT WTR DIFS -Z1'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c  AJ42
      k=k+1
      J_F1DT  = k ! F1DT (J/M**2)                                 1 GD
      name_j(k) = 'ht_cond_z1'
      lname_j(k) = 'CONDUCTION AT BOTTOM OF GROUND LAYER 1'
      units_j(k) = 'W/m^2'
      stitle_j(k)= '0CONDUCTN AT -Z1'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c  AJ43
      k=k+1
      J_ERUN2 = k ! ERUN2 (J/M**2)                                1 GP
      name_j(k) = 'ht_ice_z1z2'
      lname_j(k) = 'ENERGY OF MELTING (OR TRANS) AT -Z1-Z2'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' ICE ENRG -Z1-Z2'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c  AJ44
      k=k+1
      J_HZ0   = k ! HEATZ0=RHDT+SHDT+EVHDT+EPRCP-EDIFS (J/M**2)   1 D1
      name_j(k) = 'nt_ht_z0'
      lname_j(k) = 'NET HEATING AT GROUND SURFACE'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET HEAT AT Z0 '
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c  AJ45
      k=k+1
      J_DIFS  = k ! DIFS (KG/M**2)                                1 GP
      name_j(k) = 'h2o_difs_z1'
      lname_j(k) = 'WATER OR ICE DIFFUSION AT -Z1'
      units_j(k) = 'mm/day'
      stitle_j(k)= ' H2O DIFS AT -Z1'
      scale_j(k) = SDAY/DTSRC
      ia_j(k) = ia_src
c  AJ46
      k=k+1
      J_IMELT = k ! AIFO ; BRUN2 ; CRUN2+CIFI                     1 GP
      name_j(k) = 'ice_thru_z1z2'
      lname_j(k) = 'ICE MELTING (OR TRANSPORT) AT -Z1-Z2'
      units_j(k) = 'mm/day'
      stitle_j(k)= ' ICE THRU -Z1-Z2'
      scale_j(k) = SDAY/DTSRC
      ia_j(k) = ia_src
c  AJ47
      k=k+1
      J_RUN2  = k ! RUN2 (KG/M**2)                                1 GP
      name_j(k) = 'h2o_runoff_mld'
      lname_j(k) = 'WATER RUNOFF THROUGH MIXED LAYER DEPTH'
      units_j(k) = 'mm/day'
      stitle_j(k)= ' WATR RUNOFF MLD'
      scale_j(k) = SDAY/DTSRC
      ia_j(k) = ia_src
c  AJ48
      k=k+1
      J_DWTR2 = k ! DWTR2=AJ(DIFS)-AJ(RUN2) (KG/M**2)             1 D1
      name_j(k) = 'ht_runoff_mld'
      lname_j(k) = 'HEAT RUNOFF THROUGH MIXED LAYER DEPTH'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' HEAT RUNOFF MLD'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c  AJ49
      k=k+1
      J_WTR1  = k ! WTR1 (KG/M**2)                                1 GD
      name_j(k) = 'wat_g1'
      lname_j(k) = 'WATER IN GROUND LAYER 1'
      units_j(k) = 'kg/m^2'
      stitle_j(k)= '0WATER IN G1    '
      scale_j(k) = 1.
      ia_j(k) = ia_src
c  AJ50
      k=k+1
      J_ACE1  = k ! ACE1 (KG/M**2)                                1 GD
      name_j(k) = 'ice_g1'
      lname_j(k) = 'ICE IN GROUND LAYER 1'
      units_j(k) = 'kg/m^2'
      stitle_j(k)= ' ICE IN G1      '
      scale_j(k) = 1.
      ia_j(k) = ia_src
c  AJ51
      k=k+1
      J_WTR2  = k ! WTR2 (KG/M**2)                                1 GD
      name_j(k) = 'wat_g2'
      lname_j(k) = 'WATER IN GROUND LAYER 2'
      units_j(k) = 'kg/m^2'
      stitle_j(k)= ' WATER IN G2    '
      scale_j(k) = 1.
      ia_j(k) = ia_src
c  AJ52
      k=k+1
      J_ACE2  = k ! ACE2 (KG/M**2)                                1 GD
      name_j(k) = 'ice_g2'
      lname_j(k) = 'ICE IN GROUND LAYER 2'
      units_j(k) = 'kg/m^2'
      stitle_j(k)= ' ICE IN G2      '
      scale_j(k) = 1.
      ia_j(k) = ia_src
c  AJ53
      k=k+1
      J_SNOW  = k ! SNOW (KG/M**2)                                1 GD
      name_j(k) = 'snowdp'
      lname_j(k) = 'SNOW DEPTH'
      units_j(k) = 'kg/m^2'
      stitle_j(k)= ' SNOW DEPTH     '
      scale_j(k) = 1.
      ia_j(k) = ia_src
c  AJ54
      k=k+1
      J_RUN1  = k ! RUN1 (KG/M**2)                                1 GP
      name_j(k) = 'wat_runoff_z0'
      lname_j(k) = 'WATER RUNOFF AT GROUND SURFACE'
      units_j(k) = 'mm/day'
      stitle_j(k)= '0WATER RUNOFF Z0'
      scale_j(k) = SDAY/DTSRC
      ia_j(k) = ia_src
c  AJ55
      k=k+1
      J_BRTEMP= k ! BTEMPW-TF                                     2 RD
      name_j(k) = 'btemp_window'
      lname_j(k) = 'BRIGHTNESS TEMP THROUGH WINDOW REGION'
      units_j(k) = 'degC'
      stitle_j(k)= ' LW WINDOW BTEMP'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ56
      k=k+1
      J_HZ2   = k ! HEATZ2=AJ(F2DT)+AJ(ERUN2) (J/M**2)            1 D1
      name_j(k) = 'net_ht_z1z2'
      lname_j(k) = 'NET HEAT FLUX AT BOTTOM OF GRND LAY 2'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET HEAT -Z1-Z2'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c  AJ57
      k=k+1
      J_PCLDSS= k ! PCLDSS (1)  (COMPOSITE OVER ATMOSPHERE)       2 RD
      name_j(k) = 'sscld'
      lname_j(k) = 'SUPER SATURATION CLOUD COVER'
      units_j(k) = '%'
      stitle_j(k)= '0TOT SUP SAT CLD'
      scale_j(k) = 100.
      ia_j(k) = ia_rad
c  AJ58
      k=k+1
      J_PCLDMC= k ! PCLDMC (1)  (COMPOSITE OVER ATMOSPHERE)       2 RD
      name_j(k) = 'mccld'
      lname_j(k) = 'MOIST CONVECTIVE CLOUD COVER'
      units_j(k) = '%'
      stitle_j(k)= ' TOT MST CNV CLD'
      scale_j(k) = 100.
      ia_j(k) = ia_rad
c  AJ59
      k=k+1
      J_PCLD  = k ! PCLD (1)  (COMPOSITE OVER ATMOSPHERE)         2 RD
      name_j(k) = 'totcld'
      lname_j(k) = 'TOTAL CLOUD COVER'
      units_j(k) = '%'
      stitle_j(k)= ' TOTAL CLD COVER'
      scale_j(k) = 100.
      ia_j(k) = ia_rad
c  AJ60
      k=k+1
      J_CTOPP = k ! CLDTOPMC=AJ(CDLDEP)/AJ(PCLDMC) (100 PA)       0 D1
      name_j(k) = 'mc_clddp'
      lname_j(k) = 'MOIST CONVECTIVE CLOUD DEPTH'
      units_j(k) = 'mb'
      stitle_j(k)= ' MC CLD DPTH(MB)'
      scale_j(k) = 1.
      ia_j(k) = ia_src
c  AJ61
      k=k+1
      J_PRCPSS= k ! PRCPSS (100 PA)                               1 CN
      name_j(k) = 'ssprec'
      lname_j(k) = 'SUPER SATURATION PRECIPITATION'
      units_j(k) = 'mm/day'
      stitle_j(k)= '0SS PRECIP(MM/D)'
      scale_j(k) = 100.*SDAY/(DTsrc*GRAV)
      ia_j(k) = ia_src
c  AJ62
      k=k+1
      J_PRCPMC= k ! PRCPMC (100 PA)                               1 CN
      name_j(k) = 'mcprec'
      lname_j(k) = 'MOIST CONVECTIVE PRECIPITATION'
      units_j(k) = 'mm/day'
      stitle_j(k)= ' MC PRECIP(MM/D)'
      scale_j(k) = 100.*SDAY/(DTsrc*GRAV)
      ia_j(k) = ia_src
c  AJ63
      k=k+1
      J_QP    = k ! Q*P (100 PA)  (INTEGRAL OVER ATMOSPHERE OF)   4 DA
      name_j(k) = 'atmh2o'
      lname_j(k) = 'WATER CONTENT OF ATMOSPHERE'
      units_j(k) = 'mm'
      stitle_j(k)= ' H2O OF ATM (MM)'
      scale_j(k) = 100.*BYGRAV
      ia_j(k) = ia_dga
c  AJ64
      k=k+1
      J_GAM   = k ! GAM  (K/M)  (*SIG(TROPOSPHERE)/GRAV)          4 DA
      name_j(k) = 'lapse_rate'
      lname_j(k) = 'MEAN LAPSE RATE'
      units_j(k) = 'K/km'
      stitle_j(k)= '0GAM(K/KM)      '
      scale_j(k) = 1d3*GRAV
      ia_j(k) = ia_dga
c  AJ65
      k=k+1
      J_GAMM  = k ! GAMM  (K-S**2/M**2)  (SIG(TROPOSPHERE)/GAMD)  4 DA
      name_j(k) = 'lapse_rate_m'
      lname_j(k) = 'MOIST ADIABATIC LAPSE RATE'
      units_j(k) = 'K/km'
      stitle_j(k)= ' GAMM(K/KM)     '
      scale_j(k) = 1.D3*GAMD/(SIGE(1)-SIGE(LS1))
      ia_j(k) = ia_dga
c  AJ66
      k=k+1
      J_GAMC  = k ! GAMC  (K/M)                                   4 DA
      name_j(k) = 'lapse_rate_c'
      lname_j(k) = 'GAMC'
      units_j(k) = 'K/Km'
      stitle_j(k)= ' GAMC(K/KM)     '
      scale_j(k) = 1d3
      ia_j(k) = ia_dga
c  AJ67
      k=k+1
      J_TRINCG= k ! TRINCG (W/M**2)                               2 RD
      name_j(k) = 'lw_inc_z0'
      lname_j(k) = 'THERMAL RADIATION INCIDENT ON GROUND'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' LW INC ON Z0   '
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ68
      k=k+1
      J_FTHERM= k ! ENERGY DIFFUSION INTO THERMOCLINE (W/M**2) .5*9 MN
      name_j(k) = 'ht_thermocline'
      lname_j(k) = 'ENERGY DIFFUSION INTO THE THERMOCLINE'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' HT INTO THRMOCL'
      scale_j(k) = 2d3*SHW/SDAY
      ia_j(k) = ia_12hr
c  AJ69
      k=k+1
      J_TYPE  = k ! PTYPE                                         1 GD
      name_j(k) = 'surf_type_frac'
      lname_j(k) = 'SURF TYPE FRACT'
      units_j(k) = '%'
      stitle_j(k)= ' SURF TYPE FRACT'
      scale_j(k) = 100.
      ia_j(k) = ia_src
c  AJ70
      k=k+1
      J_HSURF = k ! TRNFP0-TRNFG (W/M**2)                         2 RD
      name_j(k) = 'HSURF'
      lname_j(k) = 'SURFACE THERMAL HEATING'
      units_j(k) = 'W/m^2'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ71
      k=k+1
      J_HATM  = k ! TRNFP1-TRNFG (W/M**2)                         2 RD
      name_j(k) = 'HATM'
      lname_j(k) = 'ATMOSPHERIC THERMAL HEATING'
      units_j(k) = 'W/m^2'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ72
      k=k+1
C**** Note: next eight diagnostics must remain in order
      J_PLAVIS= k ! PLAVIS*S0*COSZ (W/M**2)                       2 RD
      name_j(k) = 'plan_refl_vis'
      lname_j(k) = 'PLANETARY REFLECTED RAD. IN VISUAL'
      units_j(k) = 'W/M^2'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ73
      k=k+1
      J_PLANIR= k ! PLANIR*S0*COSZ (W/M**2)                       2 RD
      name_j(k) = 'plan_refl_nir'
      lname_j(k) = 'PLANETARY REFLECTED RAD. IN NEAR IR'
      units_j(k) = 'W/M^2'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ74
      k=k+1
      J_ALBVIS= k ! ALBVIS*S0*COSZ (W/M**2)                       2 RD
      name_j(k) = 'surf_refl_vis'
      lname_j(k) = 'GROUND REFLECTED RAD. IN VISUAL'
      units_j(k) = 'W/M^2'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ75
      k=k+1
      J_ALBNIR= k ! ALBNIR*S0*COSZ (W/M**2)                       2 RD
      name_j(k) = 'surf_refl_nir'
      lname_j(k) = 'GROUND REFLECTED RAD. IN NEAR IR'
      units_j(k) = 'W/M^2'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ76
      k=k+1
      J_SRRVIS= k ! SRRVIS*S0*COSZ (W/M**2)                       2 RD
      name_j(k) = 'atm_refl_vis'
      lname_j(k) = 'ATMOSPHERIC REFLECTED RAD. IN VISUAL'
      units_j(k) = 'W/M^2'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ77
      k=k+1
      J_SRRNIR= k ! SRRNIR*S0*COSZ (W/M**2)                       2 RD
      name_j(k) = 'atm_refl_nir'
      lname_j(k) = 'ATMOSPHERIC REFLECTED RAD. IN NEAR IR'
      units_j(k) = 'W/M^2'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ78
      k=k+1
      J_SRAVIS= k ! SRAVIS*S0*COSZ (W/M**2)                       2 RD
      name_j(k) = 'atm_abs_vis'
      lname_j(k) = 'ATMOSPHERIC ABSORPTION IN VISUAL'
      units_j(k) = 'W/M^2'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ79
      k=k+1
      J_SRANIR= k ! SRANIR*S0*COSZ (W/M**2)                       2 RD
      name_j(k) = 'atm_abs_nir'
      lname_j(k) = 'ATMOSPHERIC ABSORPTION IN NEAR IR'
      units_j(k) = 'W/M^2'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c  AJ80
      k=k+1
      J_CDLDEP= k ! PBOTMC-PTOPMC (100 PA)                        2 RD
      name_j(k) = 'CDLDEP'
      lname_j(k) = 'MOIST CONVECTIVE CLOUD DEPTH'
      units_j(k) = 'mb'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_CLRTOA= k ! clear sky radiative forcing top of atmosphere
      name_j(k) = 'net_clr_toa'
      lname_j(k) = 'NET CLEAR SKY RADIATION AT P0'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET CLR RAD  P0'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_CLRTRP= k ! clear sky radiative forcing tropopause
      name_j(k) = 'net_clr_trp'
      lname_j(k) = 'NET CLEAR SKY RADIATION AT TROPOPAUSE (WMO)'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET CLR RAD TRP'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c AJ83
      k=k+1
      J_TOTTRP= k ! total radiative forcing tropopause
      name_j(k) = 'net_tot_trp'
      lname_j(k) = 'NET RADIATION AT TROPOPAUSE (WMO)'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET RAD (TROPP)'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      return
      end subroutine j_defs

      subroutine ij_defs
      use DAGCOM
      use BDIJ
      use MODEL_COM
      use constant
      use DAGPCOM
      implicit none
      integer :: k
c
      do k=1,kaijx
         write(name_ij(k),'(a3,i3.3)') 'AIJ',k
         lname_ij(k) = 'unused'
         units_ij(k) = 'unused'
         ia_ij(k) = ia_src
         scale_ij(k) = 1.
         jgrid_ij(k) = 1
         iw_ij(k) = iw_all
         ir_ij(k) = ir_pct
      enddo
c
      k=0
C**** AIJ diagnostic names:
C**** NAME     NO.    DESCRIPTION   (SCALE)*IDACC  LOCATION
C**** ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
      k=k+1 !  'AIJ001'
      IJ_RSOI = k ! POICE (1)            1 GD
      lname_ij(k) = 'OCEAN/LAKE ICE COVERAGE'
      units_ij(k) = '%'
      name_ij(k) = 'RSOI'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !  'AIJ002'
      IJ_RSNW = k ! PSNOW (1)            1 GD
      lname_ij(k) = 'SNOW COVERAGE'
      units_ij(k) = '%'
      name_ij(k) = 'RSNW'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !  'AIJ003'
      IJ_SNOW = k ! SNOW (KG/M**2)       1 GD
      lname_ij(k) = 'SNOW DEPTH'    ! 'SNOW MASS'
      units_ij(k) = 'mm H2O'
      name_ij(k) = 'SNOW'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
cc    iw_ij(k) = iw_firm
      ir_ij(k) = ir_0_26_150
c
      k=k+1 !  'AIJ004'
      IJ_SHDT = k ! SHDT (J/M**2)        1 SF
      lname_ij(k) = 'SENSIBLE HEAT FLUX'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SHDT'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      ir_ij(k) = ir_m265_95
c
      k=k+1 !  'AIJ005'
      IJ_PREC = k ! PREC (KG/M**2)       1 CN
      lname_ij(k) = 'PRECIPITATION'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'PREC'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
      ir_ij(k) = ir_0_3_15
c
      k=k+1 !  'AIJ006'
      IJ_EVAP = k ! EVAP (KG/M**2)       1 SF
      lname_ij(k) = 'EVAPORATION'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'EVAP'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
      ir_ij(k) = ir_0_3_15
c
      k=k+1 !  'AIJ007'
      IJ_BETA = k ! BETA (1)             1 GD
      lname_ij(k) = 'GROUND WETNESS (VEG ROOTS)'
      units_ij(k) = '%'
      name_ij(k) = 'BETA'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      iw_ij(k) = iw_veg
c
      k=k+1 !  'AIJ008'
      IJ_PRES = k ! PIJ (100 PA)  (NO PRINTOUT)  4 DA
      lname_ij(k) = 'SURFACE PRESSURE  - PTOP'
      units_ij(k) = 'mb'
      name_ij(k) = 'PRES'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_1775
c
      k=k+1 !  'AIJ009'
      IJ_PHI1K = k ! PHI1000 (M**2/S**2) 4 DA
      lname_ij(k) = '1000mb HEIGHT'
      units_ij(k) = 'm'
      name_ij(k) = 'PHI1K'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m190_530
c
      k=k+1 !  'AIJ010'
      IJ_PHI850 = k ! PHI850 (M**2/S**2-1500*GRAV) 4 DA
      lname_ij(k) = '850 mb HEIGHT'
      units_ij(k) = 'm-1500'
      name_ij(k) = 'PHI850'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m265_95
c
      k=k+1 !  'AIJ011'
      IJ_PHI700 = k ! PHI700-3000*GRAV  4 DA
      lname_ij(k) = '700 mb HEIGHT'
      units_ij(k) = 'm-3000'
      name_ij(k) = 'PHI700'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m530_190
c
      k=k+1 !  'AIJ012'
      IJ_PHI500 = k ! PHI500-5600*GRAV  4 DA
      lname_ij(k) = '500 mb HEIGHT'
      units_ij(k) = 'm-5600'
      name_ij(k) = 'PHI500'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m1325_475
c
      k=k+1 !  'AIJ013'
      IJ_PHI300 = k ! PHI300-9500*GRAV  4 DA
      lname_ij(k) = '300 mb HEIGHT'
      units_ij(k) = 'm-9500'
      name_ij(k) = 'PHI300'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m2650_950
c
      k=k+1 !  'AIJ014'
      IJ_PHI100 = k ! PHI100-16400*GRAV 4 DA
      lname_ij(k) = '100 mb HEIGHT'
      units_ij(k) = 'm-16400'
      name_ij(k) = 'PHI100'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m2650_950
c
      k=k+1 !  'AIJ015'
      IJ_PHI30 = k ! PHI30-24000*GRAV   4 DA
      lname_ij(k) = '30 mb HEIGHT'
      units_ij(k) = 'm-24000'
      name_ij(k) = 'PHI30'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m3975_1425
c
      if (kgz_max.gt.k-IJ_PHI1K+1) then
      k=k+1 
      IJ_PHI10 = k ! PHI10-30000*GRAV   4 DA
      lname_ij(k) = '10 mb HEIGHT'
      units_ij(k) = 'm-30000'
      name_ij(k) = 'PHI10'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m5300_1900
      end if
c
      if (kgz_max.gt.k-IJ_PHI1K+1) then
      k=k+1 
      IJ_PHI3p4 = k ! PHI3.4-40000*GRAV   4 DA
      lname_ij(k) = '3.4 mb HEIGHT'
      units_ij(k) = 'm-40000'
      name_ij(k) = 'PHI3.4'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m5300_1900
      end if
c
      if (kgz_max.gt.k-IJ_PHI1K+1) then
      k=k+1 
      IJ_PHI0p7 = k ! PHI0.7-50000*GRAV   4 DA
      lname_ij(k) = '0.7 mb HEIGHT'
      units_ij(k) = 'm-50000'
      name_ij(k) = 'PHI0.7'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m5300_1900
      end if
c
      if (kgz_max.gt.k-IJ_PHI1K+1) then
      k=k+1 
      IJ_PHI0p16 = k ! PHI0.16-61000*GRAV   4 DA
      lname_ij(k) = '0.16 mb HEIGHT'
      units_ij(k) = 'm-61000'
      name_ij(k) = 'PHI0.16'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m5300_1900
      end if
c
      if (kgz_max.gt.k-IJ_PHI1K+1) then
      k=k+1 
      IJ_PHI0p07 = k ! PHI0.07-67000*GRAV   4 DA
      lname_ij(k) = '0.07 mb HEIGHT'
      units_ij(k) = 'm-67000'
      name_ij(k) = 'PHI0.07'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m5300_1900
      end if
c
      if (kgz_max.gt.k-IJ_PHI1K+1) then
      k=k+1 
      IJ_PHI0p03 = k ! PHI0.03-72000*GRAV   4 DA
      lname_ij(k) = '0.03 mb HEIGHT'
      units_ij(k) = 'm-72000'
      name_ij(k) = 'PHI0.03'
      ia_ij(k) = ia_dga
      scale_ij(k) = BYGRAV
      ir_ij(k) = ir_m5300_1900
      end if
c
      k=k+1 !  'AIJ016'
      IJ_T850 = k ! T850-TF (K-TF)*GRAV) (NO PRT) 4 DA
      lname_ij(k) = 'TEMPERATURE AT 850mb'
      units_ij(k) = 'degC'
      name_ij(k) = 'T850'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
c
      k=k+1 !  'AIJ017'
      IJ_PMCCLD = k ! PCLDMC (1)  (COMPOSITE OVER ATMOSPHERE) 2 RD
      lname_ij(k) = 'CONVECTIVE CLOUD COVER'
      units_ij(k) = '%'
      name_ij(k) = 'PMCCLD'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.
c
      k=k+1 !  'AIJ018'
      IJ_CLDTPPR = k ! P-CLOUD TOP   (100 PA)                  2 RD
      lname_ij(k) = 'CLOUD TOP PRESSURE x TOTAL CLOUD COVER'
      units_ij(k) = 'mb'
      name_ij(k) = 'CLDTPPR'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
cc    iw_ij(k) = iw_cldcv  ! built in
      ir_ij(k) = ir_0_1775
c
      k=k+1 !  'AIJ019'
      IJ_CLDCV = k ! PCLD (1)  (COMPOSITE OVER ATMOSPHERE)   2 RD
      lname_ij(k) = 'TOTAL CLOUD COVER'
      units_ij(k) = '%'
      name_ij(k) = 'CLDCV'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.
c
      k=k+1 !  'AIJ020'
      IJ_DSEV  = k
       ! 16*P4*(SHA*T4+Z4)*V1*DSIG*DXV (100 W*M/S**2) (UV GRID) 4 DA
      lname_ij(k) = 'TOTAL NT DRY STAT ENRGY' !  NT: NORTHWARD TRANSPORT
      units_ij(k) = '10^14 W'
      name_ij(k) = 'DSEV'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.d-14*100.*BYGRAV/16.
      ir_ij(k) = ir_m95_265
c
      k=k+1 !  'AIJ021'
      IJ_TRNFP0 = k ! TRNFP0 (W/M**2)                         2 RS
      lname_ij(k) = 'NET THERMAL RADIATION, TOA'   ! >0 if down !
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'TRNFP0'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_m530_190
c
      k=k+1 !  'AIJ022'
      IJ_SRTR = k ! SRHDT+TRHDT (J/M**2)                    1 RD/SF
      lname_ij(k) = 'NET RADIATION AT GROUND'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SRTR'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      ir_ij(k) = ir_m95_265
c
      k=k+1 !  'AIJ023'
      IJ_NETH = k    ! SRHDT+TRHDT+SHDT+EVHDT+ENRGP (J/M**2)   1 SC
      lname_ij(k) = 'NET HEATING AT GROUND'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'NETH'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      ir_ij(k) = ir_m530_190
c
      k=k+1 !  'AIJ024'
      IJ_SRNFP0 = k ! SRNFP0 (W/M**2)                         2 RD
      lname_ij(k) = 'NET SOLAR RADIATION, TOA'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SRNFP0'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
c
      k=k+1 !  'AIJ025'
      IJ_SRINCP0 = k ! SRINCP0 (W/M**2)                        2 RD
      lname_ij(k) = 'INCIDENT SOLAR RADIATION, TOA'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SRINCP0'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
c
      k=k+1 !  'AIJ026'
      IJ_SRNFG = k   ! SRNFG (W/M**2)                          2 RD
      lname_ij(k) = 'NET SOLAR RADIATION, SURF'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SRNFG'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
c
      k=k+1 !  'AIJ027'
      IJ_SRINCG = k ! SRINCG (W/M**2)                         2 RD
      lname_ij(k) = 'INCIDENT SOLAR RADIATION, SURF'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SRINCG'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
c
      k=k+1 !  'AIJ028'
      IJ_TG1  = k ! TG1 (K-TF)                                1 GD
      lname_ij(k) = 'GROUND TEMPERATURE'
      units_ij(k) = 'degC'
      name_ij(k) = 'TG1'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
c
      k=k+1 !  'AIJ029'
      IJ_RSIT = k ! POICE+PLICE+(IF SNOW)PEARTH               4 DA
      lname_ij(k) = 'SNOW AND ICE COVERAGE'
      units_ij(k) = '%'
      name_ij(k) = 'RSIT'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !  'AIJ030'
      IJ_TDSL = k ! DIURNAL DELTA TS (K) OVER SOIL (NO PRT) .5*9 MN
      lname_ij(k) = 'DIURNAL SURF AIR TEMP RANGE OVER SOIL'
      units_ij(k) = 'K'
      name_ij(k) = 'TDSL'
      ia_ij(k) = ia_12hr
      scale_ij(k) = 2.
      iw_ij(k) = iw_soil
      ir_ij(k) = ir_0_18
c
      k=k+1 !  'AIJ031'
      IJ_DTDP = k ! DTHETA/DPHI (K S**2/M**2) IN TROPOSPHERE  4 DA
      lname_ij(k) = 'TROP STATIC STABILITY'
      units_ij(k) = 'degC/km'
      name_ij(k) = 'DTDP'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1000.*GRAV*P1000K
      ir_ij(k) = ir_0_18
c
      k=k+1 !  'AIJ032'
      IJ_RUNE = k ! RUN1 OVER EARTH  (KG/M**2)                1 PG
      lname_ij(k) = 'GROUND RUNOFF OVER SOIL'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'RUNE'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
      iw_ij(k) = iw_soil
      ir_ij(k) = ir_m1_3
c
      k=k+1 !  'AIJ033'
      IJ_RUNLI = k ! RUN1 OVER LAND ICE  (KG/M**2) (NO PRT)    1 PG
      lname_ij(k) = 'SURFACE RUNOFF OVER LAND ICE'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'RUNLI'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
      iw_ij(k) = iw_lice
      ir_ij(k) = ir_m1_3
c
      k=k+1 !  'AIJ034'
      IJ_WS   = k ! SURFACE WIND SPEED (M/S)                  3 SF
      lname_ij(k) = 'SURFACE WIND SPEED'
      units_ij(k) = 'm/s'
      name_ij(k) = 'WS'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_18
c
      k=k+1 !  'AIJ035'
      IJ_TS   = k ! TS (K-TF)                                 3 SF
      lname_ij(k) = 'SURFACE AIR TEMPERATURE'
      units_ij(k) = 'degC'
      name_ij(k) = 'TS'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
c
      k=k+1 !  'AIJ036'
      IJ_US   = k ! US (M/S)                                  3 SF
      lname_ij(k) = 'U COMPONENT OF SURFACE AIR WIND'
      units_ij(k) = 'm/s'
      name_ij(k) = 'US'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m9_26
c
      k=k+1 !  'AIJ037'
      IJ_VS   = k ! VS (M/S)                                  3 SF
      lname_ij(k) = 'V COMPONENT OF SURFACE AIR WIND'
      units_ij(k) = 'm/s'
      name_ij(k) = 'VS'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m9_26
c
      k=k+1 !  'AIJ038'
      IJ_SLP  = k ! PSL (100 PA-1000)  (USING TS)             4 DA
      lname_ij(k) = 'SEA LEVEL PRESSURE'
      units_ij(k) = 'mb-1000'
      name_ij(k) = 'SLP'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.
      ir_ij(k) = ir_m9_26
c
      k=k+1 !  'AIJ039'
      IJ_UJET = k ! UJET (M/S)                                4 DA
      lname_ij(k) = 'U COMPONENT OF JET WINDS'
      units_ij(k) = 'm/s'
      name_ij(k) = 'UJET'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.
      jgrid_ij(k) = 2
      ir_ij(k) = ir_m38_106
c
      k=k+1 !  'AIJ040'
      IJ_VJET = k ! VJET (M/S)                                4 DA
      lname_ij(k) = 'V COMPONENT OF JET WINDS'
      units_ij(k) = 'm/s'
      name_ij(k) = 'VJET'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.
      jgrid_ij(k) = 2
      ir_ij(k) = ir_m38_106
c
      k=k+1 !  'AIJ041'
      IJ_PCLDL = k ! PCLD(LOW) (1)                            2 RD
      lname_ij(k) = 'LOW LEVEL CLOUDINESS'
      units_ij(k) = '%'
      name_ij(k) = 'PCLDL'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.
c
      k=k+1 !  'AIJ042'
      IJ_PCLDM = k ! PCLD(MID) (1)                            2 RD
      lname_ij(k) = 'MIDDLE LEVEL CLOUDINESS'
      units_ij(k) = '%'
      name_ij(k) = 'PCLDM'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.
c
      k=k+1 !  'AIJ043'
      IJ_PCLDH = k ! PCLD(HIGH) (1)                           2 RD
      lname_ij(k) = 'HIGH LEVEL CLOUDINESS'
      units_ij(k) = '%'
      name_ij(k) = 'PCLDH'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.
c
      k=k+1 !  'AIJ044'
      IJ_BTMPW = k ! BTEMPW-TF (K-TF)                         2 RD
      lname_ij(k) = 'BRIGHTNESS TEMP THRU WNDW' ! window region
      units_ij(k) = 'degC'
      name_ij(k) = 'BTMPW'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
c
      k=k+1 !  'AIJ045'
      IJ_SRREF = k ! PLAVIS*S0*COSZ (W/M**2)                  2 RD
      lname_ij(k) = 'REFLECTED SOLAR RADIATION IN VISUAL'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SRREF'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
c
      k=k+1 !  'AIJ046'
      IJ_TOC2 = k ! TGO2= TOCEAN(2)  (C)                   .5*9 MN
      lname_ij(k) = 'OCEAN TEMPERATURE BELOW MIXED LAYER' ! lyr 2
      units_ij(k) = 'degC'
      name_ij(k) = 'TOC2'
      ia_ij(k) = ia_12hr
      scale_ij(k) = 2.
      iw_ij(k) = iw_ocn
      ir_ij(k) = ir_m9_26
c
      k=k+1 !  'AIJ047'
      IJ_TAUS = k ! TAUS  (MOM. SURF. DRAG) (kg/m**2) (NO PRT)  3 SF
      lname_ij(k) = 'MAG OF MOMENTUM SURFACE DRAG'
      units_ij(k) = 'g/m*s^2'
      name_ij(k) = 'TAUS'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1000.
      ir_ij(k) = ir_0_710
c
      k=k+1 !  'AIJ048'
      IJ_TAUUS = k ! TAUUS (MOM. SURF. DRAG) (kg/m**2) (NO PRT)  3 SF
      lname_ij(k) = 'U COMPON OF MOMENTUM SRF DRAG'
      units_ij(k) = 'g/m*s^2'
      name_ij(k) = 'TAUUS'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1000.
      ir_ij(k) = ir_m2650_950
c
      k=k+1 !  'AIJ049'
      IJ_TAUVS = k ! TAUVS (MOM. SURF. DRAG) (kg/m**2) (NO PRT)  3 SF
      lname_ij(k) = 'V COMPON OF MOMENTUM SRF DRAG'
      units_ij(k) = 'g/m*s^2'
      name_ij(k) = 'TAUVS'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1000.
      ir_ij(k) = ir_m2650_950
c
      k=k+1 !  'AIJ050'
      IJ_GWTR = k ! WATER1+WATER2+ICE1+ICE2 (EARTH POINTS ONLY) 1 GD
      lname_ij(k) = 'TOTAL EARTH WATER' ! includes ice
      units_ij(k) = 'kg/m^2'
      name_ij(k) = 'GWTR'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      iw_ij(k) = iw_soil
      ir_ij(k) = ir_0_710
c
      k=k+1 !  'AIJ051'
      IJ_QS   = k ! QS                                (NO PRT)  3 SF
      lname_ij(k) = 'SURFACE AIR SPECIFIC HUMIDITY'
      units_ij(k) = '10^-4 g/g'
      name_ij(k) = 'QS'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.d4
      ir_ij(k) = ir_0_180
c
      k=k+1 !  'AIJ052'
      IJ_STRNGTS = k ! MAX(0,65F-TS_daily_avg in F)          .5*9 MN
      lname_ij(k) = 'MONTHLY HEATING' ! monthly heating need ?
      units_ij(k) = 'degF days'
      name_ij(k) = 'STRNGTS'
      ia_ij(k) = ia_12hr
      scale_ij(k) = 2.*30.
      ir_ij(k) = ir_0_3550
c
      k=k+1 !  'AIJ053'
      IJ_ARUNU = k ! ARUNU                                      1 EA
      lname_ij(k) = 'UNDERGROUND RUNOFF OVER SOIL'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'ARUNU'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
      iw_ij(k) = iw_soil
      ir_ij(k) = ir_m1_3
c
      k=k+1 !  'AIJ054'
      IJ_DTGDTS = k ! 18*(DEL(TG)/DEL(TS)-1),DEL= DIURN MX-MN .5*9 MN
      lname_ij(k) = 'PLANT WATER STRESS'
      units_ij(k) = '1'
      name_ij(k) = 'DTGDTS'
      ia_ij(k) = ia_12hr
      scale_ij(k) = 2.*30.
      ir_ij(k) = ir_m190_530
c
      k=k+1 !  'AIJ055'
      IJ_PUQ  = k ! 8*P*U*Q (VERT. INTEGRATED) (12.5 PA*M/S) 4 DA
      lname_ij(k) = 'EAST-WEST HUMIDITY FLUX (VERT SUM)'
      units_ij(k) = 'mb*m/s'
      name_ij(k) = 'PUQ'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1./8.
      ir_ij(k) = ir_m45_130
c
      k=k+1 !  'AIJ056'
      IJ_PVQ  = k ! 8*P*V*Q (VERT. INTEGRATED) (12.5 PA*M/S) 4 DA
      lname_ij(k) = 'NORTH-SOUTH HUMIDITY FLUX (VERT SUM)'
      units_ij(k) = 'mb*m/s'
      name_ij(k) = 'PVQ'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1./8.
      ir_ij(k) = ir_m45_130
c
      k=k+1 !  'AIJ057'
      IJ_TGO  = k ! TGO= TOCEAN(1)  (C)                      1 GD
      lname_ij(k) = 'MIXED-LAYER OCEAN TEMPERATURE'    ! layer 1
      units_ij(k) = 'degC'
      name_ij(k) = 'TGO'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      iw_ij(k) = iw_ocn
      ir_ij(k) = ir_m9_26
c
      k=k+1 !  'AIJ058'
      IJ_MSI2 = k ! ACE2OI= MSI2*POICE  (KG/M**2)        1 GD
      lname_ij(k) = 'LAYER 2 OCEAN ICE MASS x POICE'
      units_ij(k) = 'tons/m^2'
      name_ij(k) = 'MSI2'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-3
cc    iw_ij(k) = iw_oic !! built in
      ir_ij(k) = ir_0_4
c
      k=k+1 !  'AIJ059'
      IJ_WLM  = k ! WIND SPEED IN TOP LAYER (M/S) before SDRAG 1 SD
      lname_ij(k) = 'WIND SPEED IN TOP LAYER'
      units_ij(k) = 'm/s'
      name_ij(k) = 'WLM'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      jgrid_ij(k) = 2
      ir_ij(k) = ir_0_26_150
c
      k=k+1 !  'AIJ060'
      IJ_TGO2 = k ! TGO12= TOCEAN(3) (C)                  .5*9 MN
      lname_ij(k) = 'OCEAN TEMPERATURE AT ANN-MAX MIXED-LAYER' ! layer 3
      units_ij(k) = 'degC'
      name_ij(k) = 'TGO2'
      ia_ij(k) = ia_12hr
      scale_ij(k) = 2.
      iw_ij(k) = iw_ocn
      ir_ij(k) = ir_m9_26
c
      k=k+1 !  'AIJ061'
      IJ_EVAPO = k ! EVAP*POCEAN  (KG/M**2)                  1 GD
      lname_ij(k) = 'OCEAN EVAPORATION x POCEAN'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'EVAPO'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
c     iw built-in
      ir_ij(k) = ir_0_3_15
c
      k=k+1 !  'AIJ062'
      IJ_EVAPI = k ! EVAP*POICE  (KG/M**2)                   1 GD
      lname_ij(k) = 'OCEAN ICE EVAPORATION x POICE'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'EVAPI'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
c     iw built-in
      ir_ij(k) = ir_0_3_15
c
      k=k+1 !  'AIJ063'
      IJ_EVAPLI = k ! EVAP OVER LAND ICE  (KG/M**2)          1 GD
      lname_ij(k) = 'LAND ICE EVAPORATION x PLICE'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'EVAPLI'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
c     iw built-in
      ir_ij(k) = ir_0_3_15
c
      k=k+1 !  'AIJ064'
      IJ_EVAPE = k ! EVAP OVER EARTH  (KG/M**2)              1 GD
      lname_ij(k) = 'SOIL EVAPORATION x PSOIL'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'EVAPE'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
c     iw built-in
      ir_ij(k) = ir_0_3_15
c
      k=k+1 !  'AIJ065'
      IJ_F0OC = k ! F0DT*POCEAN, NET HEAT AT Z0  (J/M**2)    1 GD
      lname_ij(k) = 'NET HEAT INTO OCEAN x POCEAN'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'F0OC'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
c     iw built-in
      ir_ij(k) = ir_m530_190
c
      k=k+1 !  'AIJ066'
      IJ_F0OI = k ! F0DT*POICE, NET HEAT AT Z0  (J/M**2)     1 GD
      lname_ij(k) = 'NET HEAT INTO OCEAN ICE x POICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'F0OI'
      ia_ij(k) = ia_src
      scale_ij(k) = (1./DTsrc)
c     iw built-in
      ir_ij(k) = ir_m530_190
c
      k=k+1 !  'AIJ067'
      IJ_F0LI = k ! F0DT, NET HEAT AT Z0 OVER LAND ICE  (J/M**2) 1 GD
      lname_ij(k) = 'NET HEAT INTO LAND ICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'F0LI'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      iw_ij(k) = iw_lice
      ir_ij(k) = ir_m530_190
c
      k=k+1 !  'AIJ068'
      IJ_F0E  = k ! F0DT, NET HEAT AT Z0 OVER EARTH  (J/M**2) 1 GD
      lname_ij(k) = 'NET HEAT INTO SOIL'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'F0E'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      iw_ij(k) = iw_soil
      ir_ij(k) = ir_m530_190
c
      k=k+1 !  'AIJ069'
      IJ_F1LI = k ! F1DT OVER LAND ICE  (J/M**2)             1 PG
      lname_ij(k) = 'CONDUCTION AT LYR 1 BOTTOM OVER LAND ICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'F1LI'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      iw_ij(k) = iw_lice
      ir_ij(k) = ir_m530_190
c
      k=k+1 !  'AIJ070'
      IJ_SNWF = k ! SNOW FALL  (KG/M**2)                     1 PR
      lname_ij(k) = 'SNOW FALL (H2O EQUIV)'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'SNWF'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
      ir_ij(k) = ir_0_3_15
c
      k=k+1 !  'AIJ071'
      IJ_TSLI = k ! SURF AIR TEMP OVER LAND ICE  (C)  NISURF*1 SF
      lname_ij(k) = 'SURF AIR TEMP OVER LAND ICE'
      units_ij(k) = 'degC'
      name_ij(k) = 'TSLI'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d0/NIsurf
      iw_ij(k) = iw_lice
      ir_ij(k) = ir_m80_28
c
      k=k+1 !  'AIJ072'
      IJ_ERUN2 = k ! F2DT OVER LAND ICE  (J/M**2)            1 PG
      lname_ij(k) = 'TRANSPORT OF ENERGY AT LYR 2 BOTTOM OVER LAND ICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'ERUN2'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      iw_ij(k) = iw_lice
      ir_ij(k) = ir_m80_28
c
      k=k+1 !  'AIJ073'
      IJ_SHDTLI = k ! SHDT OVER LAND ICE  (J/M**2)           1 SF
      lname_ij(k) = 'SENS HEAT FLUX OVER LAND ICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'SHDTLI'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      iw_ij(k) = iw_lice
      ir_ij(k) = ir_m38_106
c
      k=k+1 !  'AIJ074'
      IJ_EVHDT = k ! EVHDT OVER LAND ICE  (J/M**2)           1 SF
      lname_ij(k) = 'LATENT HEAT FLUX OVER LAND ICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'EVHDT'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      iw_ij(k) = iw_lice
      ir_ij(k) = ir_m38_106
c
      k=k+1 !  'AIJ075'
      IJ_TRHDT = k ! TRHDT OVER LAND ICE  (J/M**2)           1 SF
      lname_ij(k) = 'NET THERMAL RADIATION INTO LAND ICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'TRHDT'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      iw_ij(k) = iw_lice
      ir_ij(k) = ir_m38_106
c
c     k=k+1 !  'AIJ076'
cfree IJ_TMAX  = k ! MAX(COMPOSITE TS)                      12 SF
cfree lname_ij(k) = 'MAX(COMPOSITE TS)'
cfree units_ij(k) = 'degK'
cfree name_ij(k) = 'TMAX'
cfree ia_ij(k) = ia_inst
cfree scale_ij(k) = 1.
c
c     k=k+1 !  'AIJ077'
cfree IJ_TMIN  = k ! MIN(COMPOSITE TS)                      12 SF
cfree lname_ij(k) = 'MIN(COMPOSITE TS)'
cfree units_ij(k) = 'degK'
cfree name_ij(k) = 'TMIN'
cfree ia_ij(k) = ia_inst
cfree scale_ij(k) = 1.
c
      k=k+1 !  'AIJ078'
      IJ_TMNMX  = k ! MIN(DIURNAL MAX OF COMPOSITE TS)      12 MN
      lname_ij(k) = 'SURFC AIR TEMPERATURE: LOWEST DIURNAL HIGH  + TF'
      units_ij(k) = 'degC' ! after offset (TF=273.16)
      name_ij(k) = 'TMNMX'
      ia_ij(k) = ia_inst
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
c
      k=k+1 !  'AIJ079'
      IJ_PEVAP = k ! POTENTIAL EVAPORATION (KG/M**2)         1 EA
      lname_ij(k) = 'POTENTIAL EVAPORATION'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'PEVAP'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
      ir_ij(k) = ir_0_26_150
c
      k=k+1 !  'AIJ080'
      IJ_TMAXE = k ! MAX TS OVER EARTH FOR CURRENT DAY (C).5*9 MN
      lname_ij(k) = 'SURFACE AIR TEMPERATURE: DIURNAL HIGH/SOIL'
      units_ij(k) = 'degC'
      name_ij(k) = 'TMAXE'
      ia_ij(k) = ia_12hr
      scale_ij(k) = 2.
      iw_ij(k) = iw_soil
      ir_ij(k) = ir_m80_28
c
      k=k+1 !  'AIJ081'
      IJ_WMSUM = k ! LIQUID WATER PATH (kg/M**2)             1 CL
      lname_ij(k) = 'LIQUID WATER PATH'
      units_ij(k) = '.1 kg/m^2'
      name_ij(k) = 'WMSUM'
      ia_ij(k) = ia_src
      scale_ij(k) = 10.
      ir_ij(k) = ir_0_18
c
      k=k+1 !  'AIJ082'
      IJ_PSCLD = k ! SHALLOW CONVECTIVE CLOUD COVER  (1)     1 CL
      lname_ij(k) = 'SHALLOW CONVECTIVE CLOUD COVER'
      units_ij(k) = '%'
      name_ij(k) = 'PSCLD'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !  'AIJ083'
      IJ_PDCLD = k ! DEEP CONVECTIVE CLOUD COVER     (1)     1 CL
      lname_ij(k) = 'DEEP CONVECTIVE CLOUD COVER'
      units_ij(k) = '%'
      name_ij(k) = 'PDCLD'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !  'AIJ084'
      IJ_DCNVFRQ = k ! DEEP CONVECTIVE CLOUD OCCURRENCE (1)   1 CL
      lname_ij(k) = 'DEEP CONVECTIVE CLOUD FREQUENCY'
      units_ij(k) = '%'
      name_ij(k) = 'DCNVFRQ'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !  'AIJ085'
      IJ_SCNVFRQ = k ! SHALLOW CONVECTIVE CLOUD OCCURRENCE (1) 1 CL
      lname_ij(k) = 'SHALLOW CONV CLOUD FREQUENCY'
      units_ij(k) = '%'
      name_ij(k) = 'SCNVFRQ'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
c     k=k+1 !  'AIJ086'
cfree IJ_EMTMOM = k ! INCIDENT MTN EAST MOM. FLUX (MB-M/S**2)  1 SD
cfree lname_ij(k) = 'INCIDENT MTN EAST MOMENTUM FLUX'
cfree units_ij(k) = 'mb m/s^2'
cfree name_ij(k) = 'EMTMOM'
cfree ia_ij(k) = ia_src
cfree scale_ij(k) = 1./DTsrc
c
c     k=k+1 !  'AIJ087'
cfree IJ_SMTMOM = k ! INCIDENT MTN SOUTH MOM. FLUX (MB-M/S**2) 1 SD
cfree lname_ij(k) = 'INCIDENT MTN SOUTH MOMENTUM FLUX'
cfree units_ij(k) = 'mb m/s^2'
cfree name_ij(k) = 'SMTMOM'
cfree ia_ij(k) = ia_src
cfree scale_ij(k) = 1./DTsrc
c
      k=k+1 !  'AIJ088'
      IJ_FPEU = k ! EAST-WEST POT. ENTHALPY FLUX (W)   /3600.*1 DY
      lname_ij(k) = 'EAST-WEST POTENTIAL ENTHALPY FLUX'
      units_ij(k) = '10^10 W'
      name_ij(k) = 'FPEU'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-10/DTsrc
      ir_ij(k) = ir_m38_106
c
      k=k+1 !  'AIJ089'
      IJ_FPEV = k ! NORTH-SOUTH POT. ENTHALPY FLUX (W) /3600.*1 DY
      lname_ij(k) = 'NORTH-SOUTH POTENTIAL ENTHALPY FLUX'
      units_ij(k) = '10^10 W'
      name_ij(k) = 'FPEV'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-10/DTsrc
      ir_ij(k) = ir_m38_106
c
      k=k+1 !  'AIJ090'
      IJ_FMU  = k ! EAST-WEST MASS FLUX (KG/S) 100./GRAV/3600.*1 DY
      lname_ij(k) = 'EAST-WEST MASS FLUX'
      units_ij(k) = '10^10 kg/s'
      name_ij(k) = 'FMU'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-10*100.*BYGRAV/DTsrc
      ir_ij(k) = ir_m38_106
c
      k=k+1 !  'AIJ091'
      IJ_FMV  = k ! NORTH-SOUTH MASS FLUX (KG/S) 100./GRAV/3600.*1 DY
      lname_ij(k) = 'NORTH-SOUTH MASS FLUX'
      units_ij(k) = '10^10 kg/s'
      name_ij(k) = 'FMV'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-10*100.*BYGRAV/DTsrc
c
      k=k+1 !  'AIJ092'
      IJ_FQU = k ! EAST-WEST WATER VAPOR FLUX (KG/S)   /3600.*1 DY
      lname_ij(k) = 'EAST-WEST WATER VAPOR FLUX'
      units_ij(k) = '10^6 kg/s'
      name_ij(k) = 'FQU'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-6/DTsrc
      ir_ij(k) = ir_m80_28
c
      k=k+1 !  'AIJ093'
      IJ_FQV = k ! NORTH-SOUTH WATER VAPOR FLUX (KG/S) /3600.*1 DY
      lname_ij(k) = 'NORTH-SOUTH WATER VAPOR FLUX'
      units_ij(k) = '10^6 kg/s'
      name_ij(k) = 'FQV'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-6/DTsrc
      ir_ij(k) = ir_m80_28
c
      k=k+1 !  'AIJ094'
      IJ_FGZU = k ! EAST-WEST GEOPOTENTIAL FLUX (W)    /3600.*1 DY
      lname_ij(k) = 'EAST-WEST GEOPOTENTIAL FLUX'
      units_ij(k) = '10^10 W'
      name_ij(k) = 'FGZU'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-10/DTsrc
      ir_ij(k) = ir_m1325_475
c
      k=k+1 !  'AIJ095'
      IJ_FGZV = k ! NORTH-SOUTH GEOPOTENTIAL FLUX (W)  /3600.*1 DY
      lname_ij(k) = 'NORTH-SOUTH GEOPOTENTIAL FLUX'
      units_ij(k) = '10^10 W'
      name_ij(k) = 'FGZV'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-10/DTsrc
      ir_ij(k) = ir_m1325_475
c
      k=k+1 !  'AIJ096'
      IJ_ERVR = k ! Energy Outflow by Rivers (10**10 W) E-10/DTS*1 RV
      lname_ij(k) = 'Energy Outflow by Rivers'
      units_ij(k) = '10^10 W'
      name_ij(k) = 'ERVR'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-10/DTsrc
      ir_ij(k) = ir_m1325_475
c
      k=k+1 !  'AIJ097'
      IJ_MRVR = k ! Mass Outflow by Rivers (10**5 kg/s)  E-5/DTS*1 RV
      lname_ij(k) = 'Mass Outflow by Rivers'
      units_ij(k) = '100 tons/s'
      name_ij(k) = 'MRVR'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-5/DTsrc
      ir_ij(k) = ir_m1325_475
c
      k=k+1 !  'AIJ098'
      IJ_SDRAG = k ! DU/DT BY SDRAG (M S-2)                       1 SD
      lname_ij(k) = 'STRATOSPHERIC DRAG'
      units_ij(k) = '10^-6 m/s/s'
      name_ij(k) = 'SDRAG'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d6/DTsrc
      ir_ij(k) = ir_m190_530
c
      k=k+1 !  'AIJ099'
      IJ_LKON = k
      lname_ij(k) = 'LAST DAY OF ICE-FREE LAKE'
      units_ij(k) = 'JULIAN DAY'
      name_ij(k) = 'LKONDAY'
      ia_ij(k) = ia_inst        ! only works for ann.means w/ NMONAV=1 ?
      iw_ij(k) = iw_lake
      ir_ij(k) = ir_0_180
c
      k=k+1 !  'AIJ100'
      IJ_LKOFF = k
      lname_ij(k) = 'LAST DAY OF ICED-UP LAKE'
      units_ij(k) = 'JULIAN DAY'
      name_ij(k) = 'LKOFFDAY'
      ia_ij(k) = ia_inst        ! only works for ann.means w/ NMONAV=1 ?
      iw_ij(k) = iw_lake
      ir_ij(k) = ir_0_180
c
C**** Here I am adding all the previous AIJG to AIJ
C**** actual number for Gxx = 100 + xx
      k=k+1
      IJ_G01 = k
      name_ij(k) = 'AIJG01'
      lname_ij(k) = 'LAYER 1 BARE SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      iw_ij(k) = iw_bare
      ir_ij(k) = ir_0_710
c
      k=k+1
      IJ_G02 = k
      name_ij(k) = 'AIJG02'
      lname_ij(k) = 'LAYER 2 BARE SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      iw_ij(k) = iw_bare
      ir_ij(k) = ir_0_710
c
      k=k+1
      IJ_G03 = k
      name_ij(k) = 'AIJG03'
      lname_ij(k) = 'LAYER 3 BARE SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      iw_ij(k) = iw_bare
      ir_ij(k) = ir_0_710
c
      k=k+1
      IJ_G04 = k
      name_ij(k) = 'AIJG04'
      lname_ij(k) = 'LAYER 4 BARE SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      iw_ij(k) = iw_bare
      ir_ij(k) = ir_0_710
c
      k=k+1
      IJ_G05 = k
      name_ij(k) = 'AIJG05'
      lname_ij(k) = 'BARE SOIL WETNESS, BETA'
      units_ij(k) = '%'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      iw_ij(k) = iw_bare
c
      k=k+1
      IJ_G06 = k
      name_ij(k) = 'AIJG06'
      lname_ij(k) = 'PENMAN SOIL WETNESS, BETA'
      units_ij(k) = '%'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      iw_ij(k) = iw_soil
c
      k=k+1
      IJ_G07 = k
      name_ij(k) = 'AIJG07'
      lname_ij(k) = 'VEGETATION CANOPY SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_0_710
c
      k=k+1
      IJ_G08 = k
      name_ij(k) = 'AIJG08'
      lname_ij(k) = 'LAYER 1 VEGETATED SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_0_710
c
      k=k+1
      IJ_G09 = k
      name_ij(k) = 'AIJG09'
      lname_ij(k) = 'LAYER 2 VEGETATED SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_0_710
c
      k=k+1
      IJ_G10 = k
      name_ij(k) = 'AIJG10'
      lname_ij(k) = 'LAYER 3 VEGETATED SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_0_710
c
      k=k+1
      IJ_G11 = k
      name_ij(k) = 'AIJG11'
      lname_ij(k) = 'BARE & VEGETATED SOIL WETNESS'
      units_ij(k) = '%'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      iw_ij(k) = iw_soil
c
      k=k+1
      IJ_G12 = k
      name_ij(k) = 'AIJG12'
      lname_ij(k) = 'CONDUCTANCE OF ATMOSPHERE'
      units_ij(k) = 'm/s'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      iw_ij(k) = iw_soil
      ir_ij(k) = ir_0_71
c
      k=k+1
      IJ_G13 = k
      name_ij(k) = 'AIJG13'
      lname_ij(k) = 'CONDUCTANCE OF CANOPY'
      units_ij(k) = '.01 m/s'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_0_71
c
      k=k+1
      IJ_G14 = k
      name_ij(k) = 'AIJG14'
      lname_ij(k) = 'PENMAN POTENTIAL EVAPORATION'
      units_ij(k) = 'mm/day'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
      iw_ij(k) = iw_soil
      ir_ij(k) = ir_0_26_150
c
      k=k+1
      IJ_G15 = k
      name_ij(k) = 'AIJG15'
      lname_ij(k) = 'BARE SOIL LAYER 1 TEMPERATURE'
      units_ij(k) = 'degC'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./NIsurf
      iw_ij(k) = iw_bare
      ir_ij(k) = ir_m80_28
c
      k=k+1
      IJ_G16 = k
      name_ij(k) = 'AIJG16'
      lname_ij(k) = 'BARE SOIL LAYER 1 TEMPERATURE'
      units_ij(k) = 'degC'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./NIsurf
      iw_ij(k) = iw_bare
      ir_ij(k) = ir_m80_28
c
      k=k+1
      IJ_G17 = k
      name_ij(k) = 'AIJG17'
      lname_ij(k) = 'BARE SOIL LAYER 3 TEMPERATURE'
      units_ij(k) = 'degC'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./NIsurf
      iw_ij(k) = iw_bare
      ir_ij(k) = ir_m80_28
c
      k=k+1
      IJ_G18 = k
      name_ij(k) = 'AIJG18'
      lname_ij(k) = 'BARE SOIL EVAPORATION'
      units_ij(k) = 'mm/day'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
      iw_ij(k) = iw_bare
      ir_ij(k) = ir_0_3_15
c
      k=k+1
      IJ_G19 = k
      name_ij(k) = 'AIJG19'
      lname_ij(k) = 'DRY CANOPY EVAPORATION'
      units_ij(k) = 'mm/day'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_0_3_15
c
      k=k+1
      IJ_G20 = k
      name_ij(k) = 'AIJG20'
      lname_ij(k) = 'WET CANOPY EVAPORATION'
      units_ij(k) = 'mm/day'
      ia_ij(k) = ia_src
      scale_ij(k) = SDAY/DTsrc
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_0_3_15
c
      k=k+1
      IJ_G21 = k
      name_ij(k) = 'AIJG21'
      lname_ij(k) = 'CANOPY TEMPERATURE'
      units_ij(k) = 'degC'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./NIsurf
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_m80_28
c
      k=k+1
      IJ_G22 = k
      name_ij(k) = 'AIJG22'
      lname_ij(k) = 'VEGETATED SOIL LAYER 1 TEMPERATURE'
      units_ij(k) = 'degC'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./NIsurf
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_m80_28
c
      k=k+1
      IJ_G23 = k
      name_ij(k) = 'AIJG23'
      lname_ij(k) = 'VEGETATED SOIL LAYER 2 TEMPERATURE'
      units_ij(k) = 'degC'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./NIsurf
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_m80_28
c
      k=k+1
      IJ_G24 = k
      name_ij(k) = 'AIJG24'
      lname_ij(k) = 'VEGETATED SOIL LAYER 3 TEMPERATURE'
      units_ij(k) = 'degC'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./NIsurf
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_m80_28
c
      k=k+1
      IJ_G25 = k
      name_ij(k) = 'AIJG25'
      lname_ij(k) = 'AVERAGE WATER TABLE DEPTH'
      units_ij(k) = 'm'
      ia_ij(k) = ia_src
      scale_ij(k) = -1./NIsurf
      iw_ij(k) = iw_soil
      ir_ij(k) = ir_0_3_15
c
      k=k+1
      IJ_G26 = k
      name_ij(k) = 'AIJG26'
      lname_ij(k) = 'VEGETATED SOIL WETNESS'
      units_ij(k) = '%'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      iw_ij(k) = iw_veg
c
      k=k+1
      IJ_G27 = k
      name_ij(k) = 'AIJG27'
      lname_ij(k) = 'TRANSPIRATION EFFICIENCY, BETAT'
      units_ij(k) = '%'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      iw_ij(k) = iw_veg
c
      k=k+1
      IJ_G28 = k
      name_ij(k) = 'AIJG28'
      lname_ij(k) = 'SNOW DEPTH OVER BARE SOIL'
      units_ij(k) = 'mm H2O'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_0_3550
c
      k=k+1
      IJ_G29 = k
      name_ij(k) = 'AIJG29'
      lname_ij(k) = 'SNOW DEPTH OVER VEG SOIL'
      units_ij(k) = 'mm H2O'
      scale_ij(k) = 1000.
      ia_ij(k) = ia_src
      iw_ij(k) = iw_veg
      ir_ij(k) = ir_0_3550
c
c Gravity Wave diagnostics   
      iDO_GWDRAG = 0
      if (DO_GWDRAG) then
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW1 = k 
      name_ij(k) = 'ij_def_drag_mom_flux'
      lname_ij(k) = 'DEFORM. DRAG MOM FLUX'
      units_ij(k) = '.1 N/M**2'
      scale_ij(k) = 10.*100.*BYGRAV
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m1_3
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW2 = k
      name_ij(k) = 'ij_mtn_wave_mom_flux'
      lname_ij(k) = 'MTN WAVE MOM. FLUX'
      units_ij(k) = '.1 N/M**2'  ! dynes/cm^2
      scale_ij(k) = 10.*100.*BYGRAV
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m1_3
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW3 = k
      name_ij(k) = 'ij_shr_wave_mom_flux'
      lname_ij(k) = 'SHEAR WAVE MOM. FLUX'
      units_ij(k) = '.001 N/M**2'
      scale_ij(k) = 1000.*100.*BYGRAV
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m1_3
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW4 = k
      name_ij(k) = 'ij_mc_c_m10r_mom_flux'
      lname_ij(k) = 'MC C=-10R MOM. FLUX'
      units_ij(k) = '.001 N/M**2'
      scale_ij(k) = 1000.*100.*BYGRAV
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m1_3
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW5 = k
      name_ij(k) = 'ij_mc_c_m20r_mom_flux'
      lname_ij(k) = 'MC C=-20R MOM. FLUX'
      units_ij(k) = '.001 N/M**2'
      scale_ij(k) = 1000.*100.*BYGRAV
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m1_3
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW6 = k
      name_ij(k) = 'ij_mc_c_m40r_mom_flux'
      lname_ij(k) = 'MC C=-40R MOM. FLUX'
      units_ij(k) = '.001 N/M**2'
      scale_ij(k) = 1000.*100.*BYGRAV
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m1_3
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW7 = k
      name_ij(k) = 'ij_phase_speed_of_wind_shear'
      lname_ij(k) = 'PHASE SPEED OF SHEAR WAVE'
      units_ij(k) = 'M/S'
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m45_130
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW8 = k
      name_ij(k) = 'ij_source_speed_of_mc'
      lname_ij(k) = 'MC SOURCE WIND SPEED'
      units_ij(k) = 'M/S'
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m45_130
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW9 = k
      name_ij(k) = 'ij_exit_tot_mom_flux'
      lname_ij(k) = 'EXIT TOT. MOM. FLUX'
      units_ij(k) = '.0001 N/M**2'
      scale_ij(k) = 10000.*100.*BYGRAV
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m1_3
      iDO_GWDRAG = k-IJ_GW1+1
      END IF

      if (k .gt. kaij) then
        write (6,*) 'ij_defs: Increase kaij=',kaij,' to at least ',k
        stop 'kaij too small'
      end if

      return
      end subroutine ij_defs

      subroutine il_defs
      USE CONSTANT, only : grav,rgas,by3,sha,bygrav
      USE MODEL_COM, only : dtsrc,jeq
      USE GEOM, only : dxyp
      use DAGCOM
      implicit none
      real*8 :: bydj,bydjuv,daeq
      integer :: k,j
c
      do k=1,kail
         write(name_il(k),'(a3,i3.3)') 'AIL',k
         lname_il(k) = 'unused'
         units_il(k) = 'unused'
         scale_il(k) = 1.
         ia_il(k)    = 0.
      enddo

C**** some scaling numbers for the equatorial diags.
      bydj   = 1./dble(j5n-j5s+1)
      bydjuv = 1./dble(j5nuv-j5suv+1)
      daeq=0.
      do j=j5s,j5n
        daeq=daeq+DXYP(J)
      end do
C****
      k=0
c
      k = k + 1
      IL_UEQ=k
      name_il(k) = 'u_equator'
      lname_il(k) = 'ZONAL WIND (U COMPONENT) AROUND +/- 5 DEG'
      units_il(k) = 'METERS/SECOND'
      scale_il(k) = bydjuv
      ia_il(k)    = ia_dga
      k = k + 1
      IL_VEQ=k
      name_il(k) = 'v_equator'
      lname_il(k) = 'MERIDIONAL WIND (V COMPONENT) AROUND +/- 5 DEG'
      units_il(k) = 'METERS/SECOND'
      scale_il(k) = bydjuv
      ia_il(k)    = ia_dga
      k = k + 1
      IL_WEQ=k
      name_il(k) = 'vvel_equator'
      lname_il(k) = 'VERTICAL VELOCITY AROUND +/- 5 DEG'
      units_il(k) = '10**-4 METERS/SECOND'
      scale_il(k) = -1d4*RGAS*BYGRAV/daeq
      ia_il(k)    = ia_dga
      k = k + 1
      IL_TEQ=k
      name_il(k) = 'temp_equator'
      lname_il(k) = 'TEMPERATURE AROUND +/- 5 DEG'
      units_il(k) = 'DEGREES CENTIGRADE'
      scale_il(k) = bydj
      ia_il(k)    = ia_dga
      k = k + 1
      IL_QEQ=k
      name_il(k) = 'rh_equator'
      lname_il(k) = 'RELATIVE HUMIDITY AROUND +/- 5 DEG'
      units_il(k) = 'PERCENT'
      scale_il(k) = 1d2*bydj
      ia_il(k)    = ia_dga
      k = k + 1
      IL_MCEQ=k
      name_il(k) = 'mcheat_equator'
      lname_il(k) = 'MOIST CONVECTIVE HEATING AROUND +/- 5 DEG'
      units_il(k) = '10**13 WATTS/DSIG'
      scale_il(k) = 100d-13*SHA/(GRAV*DTsrc)
      ia_il(k)    = ia_src
      k = k + 1
      IL_REQ=k
      name_il(k) = 'rad_cool_equator'
      lname_il(k) = 'TOTAL RADIATIVE COOLING AROUND +/- 5 DEG'
      units_il(k) = '10**13 WATTS/DSIG'
      scale_il(k) = -1d-13
      ia_il(k)    = ia_rad
      k = k + 1
c      name_il(k) = ''
c      lname_il(k) = ''
c      units_il(k) = ''
c      scale_il(k) = 1.
c      ia_il(k)    = 0
      k = k + 1
      IL_W50N=k
      name_il(k) = 'vvel_50N'
      lname_il(k) = 'VERTICAL VELOCITY AT 50 N'
      units_il(k) = '10**-4 METERS/SECOND'
      scale_il(k) = 1d4*RGAS/(GRAV*DXYP(J50N))
      ia_il(k)    = ia_dga
      k = k + 1
      IL_T50N=k
      name_il(k) = 'temp_50N'
      lname_il(k) = 'TEMPERATURE AT 50 N'
      units_il(k) = 'DEGREES CENTIGRADE'
      scale_il(k) = 1.
      ia_il(k)    = ia_dga
      k = k + 1
      IL_R50N=k
      name_il(k) = 'rad_cool_50N'
      lname_il(k) = 'TOTAL RADIATIVE COOLING AT 50 N'
      units_il(k) = '10**13 WATTS/UNIT SIGMA'
      scale_il(k) = 1d-13
      ia_il(k)    = ia_rad
      k = k + 1
      IL_U50N=k
      name_il(k) = 'u_50N'
      lname_il(k) = 'ZONAL WIND AT 50 N'
      units_il(k) = 'METERS/SECOND'
      scale_il(k) = 0.5
      ia_il(k)    = ia_dga
      k = k + 1
      IL_W70N=k
      name_il(k) = 'vvel_70N'
      lname_il(k) = 'VERTICAL VELOCITY AT 70 N'
      units_il(k) = '10**-4 METERS/SECOND'
      scale_il(k) = -1d4*RGAS/(GRAV*DXYP(J70N))
      ia_il(k)    = ia_dga
      k = k + 1
      IL_T70N=k
      name_il(k) = 'temp_70N'
      lname_il(k) = 'TEMPERATURE AT 70 N'
      units_il(k) = 'DEGREES CENTIGRADE'
      scale_il(k) = 1.
      ia_il(k)    = ia_dga
      k = k + 1
      IL_R70N=k
      name_il(k) = 'rad_cool_70N'
      lname_il(k) = 'TOTAL RADIATIVE COOLING AT 70 N'
      units_il(k) = '10**13 WATTS/UNIT SIGMA'
      scale_il(k) = -1d-13
      ia_il(k)    = ia_rad
      k = k + 1
      IL_U70N=k
      name_il(k) = 'u_70N'
      lname_il(k) = 'ZONAL WIND AT 70 N'
      units_il(k) = 'METERS/SECOND'
      scale_il(k) = 0.5
      ia_il(k)    = ia_dga
c
      return
      end subroutine il_defs


      subroutine jl_defs
      use CONSTANT, only : sday,grav,twopi,sha,rgas,bygrav,radius
      use MODEL_COM, only : fim,dtsrc,nidyn,byim
      use GEOM, only : dlon
      use DAGCOM
      implicit none
      integer :: k
c
      do k=1,kajlx
         write(sname_jl(k),'(a3,i3.3)') 'AJL',k
         lname_jl(k) = 'unused'
         units_jl(k) = 'unused'
      enddo
c
      k=0
c
      k=k+1
      jl_mcmflx = k
      sname_jl(k) = 'mc_mflx' !'AJL08'                  'FMX(MC)*P'
      lname_jl(k) = 'VERTICAL MASS EXCHANGE FROM MOIST CONVECTION'
      units_jl(k) = '10**9 KG/SECOND' !'100 PA'
      scale_jl(k) = 100.D-9*TWOPI/(GRAV*DTsrc*DLON*FIM)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_srhr = k
      sname_jl(k) = 'srad_heat' !'AJL09'
      lname_jl(k) = 'SOLAR RADIATION HEATING RATE' !'SRHR'
      units_jl(k) = 'HUNDREDTHS OF DEGREES KELVIN/DAY' !'W/m^2'
      scale_jl(k) = 100.D-2*GRAV*SDAY/SHA
      ia_jl(k) = ia_rad
      jgrid_jl(k) = 1
c
      k=k+1
      jl_trcr = k
      sname_jl(k) = 'trad_cool' !'AJL10'
      lname_jl(k) = 'THERMAL RADIATION COOLING RATE' !'TRHR'
      units_jl(k) = 'HUNDREDTHS OF DEGREES K/DAY' !'W/m^2'
      scale_jl(k) = -100.D-2*GRAV*SDAY/SHA
      ia_jl(k) = ia_rad
      jgrid_jl(k) = 1
c
      k=k+1
      jl_sshr = k
      sname_jl(k) = 'lscond_heat' !'AJL11'
      lname_jl(k) = 'HEATING BY LARGE SCALE CONDENSATION' !'DTX(SS)*P'
      units_jl(k) = '10**13 WATTS/UNIT SIGMA' !'100 K*PA'
      scale_jl(k) = 100.D-13*SHA*TWOPI/(GRAV*DTsrc*DLON*FIM)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_trbhr = k
      sname_jl(k) = 'turb_heat' !'AJL12'
      lname_jl(k) = 'HEATING BY TURBULENCE' !'DT(DC)*P'
      units_jl(k) = '10**13 WATTS/UNIT SIGMA' !'100 K*PA'
      scale_jl(k) = 100.D-13*SHA*TWOPI/(GRAV*DTsrc*DLON*FIM)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_mchr = k
      sname_jl(k) = 'AJL13'
      lname_jl(k) = 'DT(MC)*P  DRY HEATING'
      units_jl(k) = '100 PA*K'
c
!     k=k+1
!     jl_free14 = k
!     sname_jl(k) = 'AJL14'
!     lname_jl(k) = 'unknown'
!     units_jl(k) = 'unknown'
c
!     k=k+1
!     jl_free15 = k
!     sname_jl(k) = 'AJL15'
!     lname_jl(k) = 'unknown'
!     units_jl(k) = 'unknown'
c
      k=k+1
      jl_ape = k
      sname_jl(k) = 'avail_pe' !'AJL16'
      lname_jl(k) = 'AVAILABLE POTENTIAL ENERGY'
      units_jl(k) = '10**5 JOULES/M**2/UNIT SIGMA'
      scale_jl(k) = 50.D-5*RGAS*BYIM*BYGRAV
      ia_jl(k) = ia_dga
      jgrid_jl(k) = 1
c
      k=k+1
      jl_dtdyn = k
      sname_jl(k) = 'DT_DYNAMICS' !'AJL17'
      lname_jl(k) = 'DTEMP/DT BY DYNAMICS'
      units_jl(k) = '10**-1 DEG-K/DAY'
      scale_jl(k) = 1.E1*SDAY*NIDYN/(FIM*7200.)
      ia_jl(k) = ia_dga
      jgrid_jl(k) = 1
c
      k=k+1
      jl_totcld = k
      sname_jl(k) = 'totcld' !'AJL19'
      lname_jl(k) = 'TOTAL CLOUD COVER' !'PCLD*P (TOTAL)'
      units_jl(k) = 'PERCENT'
      scale_jl(k) = 100.*BYIM
      ia_jl(k) = ia_rad
      jgrid_jl(k) = 1
c
      k=k+1
      jl_gwFirst = k   ! The next consececutive n are Gravity Wave
      jl_dudfmdrg = k
      sname_jl(k) = 'dudt_dfmdrg' !'AJL18'
      lname_jl(k) = 'DU/DT BY STRAT DEFORM DRAG'
      units_jl(k) = '10**-6 M/S/S'
      scale_jl(k) = 1.D6/(FIM*DTsrc)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_dumtndrg = k
      sname_jl(k) = 'dudt_mtndrg' !'AJL20'
      lname_jl(k) = 'DU/DT BY STRAT MTN DRAG'
      units_jl(k) = '10**-6 M S-2' ! ??? 'M/S'
      scale_jl(k) = 1.D6/(FIM*DTsrc)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_dushrdrg = k
      sname_jl(k) = 'dudt_shrdrg'
      lname_jl(k) = 'DU/DT BY STRAT SHR DRAG'
      units_jl(k) = '10**-6 M/S/S'
      scale_jl(k) = 1.D6/(FIM*DTsrc)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_dumcdrgm10 = k
      sname_jl(k) = 'dudt_mcdrgm10' !'AJL22'
      lname_jl(k) = 'DU/DT BY STRAT MC DRAG C=-10'
      units_jl(k) = '10**-6 M/S/S'
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_dumcdrgp10 = k
      sname_jl(k) = 'dudt_mcdrgp10' !'AJL23'
      lname_jl(k) = 'DU/DT BY STRAT MC DRAG C=+10'
      units_jl(k) = '10**-6 M/S/S'
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_dumcdrgm40 = k
      sname_jl(k) = 'dudt_mcdrgm40' !'AJL24'
      lname_jl(k) = 'DU/DT BY STRAT MC DRAG C=-40'
      units_jl(k) = '10**-6 M/S/S'
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_dumcdrgp40 = k
      sname_jl(k) = 'dudt_mcdrgp40' !'AJL25'
      lname_jl(k) = 'DU/DT BY STRAT MC DRAG C=+40'
      units_jl(k) = '10**-6 M/S/S'
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_dumcdrgm20 = k
      sname_jl(k) = 'dudt_mcdrgm20' !'AJL26'
      lname_jl(k) = 'DU/DT BY STRAT MC DRAG C=-20'
      units_jl(k) = '10**-6 M/S/S'
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c Last of the Gravity Wave JL's
      k=k+1
      jl_dumcdrgp20 = k
      sname_jl(k) = 'dudt_mcdrgp20' !'AJL27'
      lname_jl(k) = 'DU/DT BY STRAT MC DRAG C=+20'
      units_jl(k) = '10**-6 M/S/S'
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_sscld = k
      sname_jl(k) = 'sscld' !'AJL28'
      lname_jl(k) = 'SUPER SATURATION CLOUD COVER' !'PCLD*P (SS)'
      units_jl(k) = 'PERCENT'
      scale_jl(k) = 100.*BYIM
      ia_jl(k) = ia_rad
      jgrid_jl(k) = 1
c
      k=k+1
      jl_mccld = k
      sname_jl(k) = 'mccld' !'AJL29'
      lname_jl(k) = 'MOIST CONVECTIVE CLOUD COVER' !'PCLD*P (MC)'
      units_jl(k) = 'PERCENT'
      scale_jl(k) = 100.*BYIM
      ia_jl(k) = ia_rad
      jgrid_jl(k) = 1
c
!     k=k+1
!     jl_free30 = k
!     sname_jl(k) = 'AJL30'
!     lname_jl(k) = 'unknown'
!     units_jl(k) = 'unknown'
c
      k=k+1
      jl_sdifcoef = k
      sname_jl(k) = 'strat_diff_coeff' !'AJL31'
      lname_jl(k) = 'STRAT. DIFFUSION COEFF'
      units_jl(k) = 'M*M/S'
c
      k=k+1
      jl_dudtsdif = k
      sname_jl(k) = 'dudt_sdiff' !'AJL32'     ! gwdrag
      lname_jl(k) = 'DU/DT  STRATOSPHERIC DIFFUSION'
      units_jl(k) = 'M/S'
c
      k=k+1
      jl_dtdtsdrg = k
      sname_jl(k) = 'DT_SDRAG' !'AJL33'
      lname_jl(k) = 'DTEMP/DT BY STRATOSPHERIC DRAG'
      units_jl(k) = '10**-1 DEG-K/DAY'
      scale_jl(k) = 1.E1*SDAY/(FIM*DTsrc)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
!     k=k+1
!     jl_free34 = k
!     sname_jl(k) = 'AJL34'
!     lname_jl(k) = 'unknown'
!     units_jl(k) = 'unknown'
c
!     k=k+1
!     jl_free35 = k
!     sname_jl(k) = 'AJL35'
!     lname_jl(k) = 'unknown'
!     units_jl(k) = 'unknown'
c
      k=k+1
      jl_epflxv = k
      sname_jl(k) = 'epflx_vert' !'AJL36'
      lname_jl(k) = 'VERTICAL ELIASSEN-PALM FLUX'
      units_jl(k) = '10**17 JOULES'
      scale_jl(k) = .125*100.D-17*TWOPI*RADIUS*BYGRAV/(DLON*FIM)
      ia_jl(k) = ia_dga
      jgrid_jl(k) = 1
c
      k=k+1
      jl_epflxn = k
      sname_jl(k) = 'epflx_north' !'AJL37'
      lname_jl(k) = 'NORTHWARD ELIASSEN-PALM FLUX'
      units_jl(k) = '10**17 JOULES/UNIT SIGMA'
      scale_jl(k) = 100.D-17*TWOPI*RADIUS*BYGRAV/(DLON*FIM)
      ia_jl(k) = ia_dga
      jgrid_jl(k) = 2
c
      k=k+1
      jl_damdc = k
      sname_jl(k) = 'del_am_dc' !'AJL38'
      lname_jl(k) = 'CHANGE OF ANG. MOMENTUM BY TURBULENCE'
      units_jl(k) = '10**18 JOULE/UNIT SIGMA'
      scale_jl(k) = 100.D-18*TWOPI*RADIUS/(GRAV*DTsrc*DLON*FIM)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_dammc = k
      sname_jl(k) = 'del_am_mc' !'AJL39'          'DU(MC)*P (UV GRID)'
      lname_jl(k) = 'CHANGE OF ANG. MOMENTUM BY MOIST CONV' ! UV GRID
      units_jl(k) = '10**18 JOULE/UNIT SIGMA' !'100 N/m/s'
      scale_jl(k) = 100.D-18*TWOPI*RADIUS/(GRAV*DTsrc*DLON*FIM)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_40 = k
      sname_jl(k) = 'AJL40' !DU(ED)*P*(DTSURF*DSIG*ED/DZ**2)  (UV GRID)
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      jl_uepac = k
      sname_jl(k) = 'u_epac' !'AJL41'
      lname_jl(k) = 'U WIND AVERAGED OVER I,EAST PACIFIC'
      units_jl(k) = 'TENTHS OF METERS/SECOND'
      scale_jl(k) = .2E+1
      ia_jl(k) = ia_dga
      jgrid_jl(k) = 2
c
      k=k+1
      jl_vepac = k
      sname_jl(k) = 'v_epac' !'AJL42'
      lname_jl(k) = 'V WIND AVERAGED OVER EAST PACIFIC'
      units_jl(k) = 'TENTHS OF METERS/SECOND'
      scale_jl(k) = .2E+1
      ia_jl(k) = ia_dga
      jgrid_jl(k) = 2
c
      k=k+1
      jl_wepac = k
      sname_jl(k) = 'vvel_epac' !'AJL43'
      lname_jl(k) = 'VERTICAL VELOCITY FOR EAST PACIFIC'
      units_jl(k) = '10**-5 METERS/SECOND'
      scale_jl(k) = -1.D5*RGAS/(5.*GRAV)
      ia_jl(k) = ia_dga
      jgrid_jl(k) = 1
c
      k=k+1
      jl_uwpac = k
      sname_jl(k) = 'u_wpac' !'AJL44'
      lname_jl(k) = 'U WIND AVERAGED OVER I=WEST PACIFIC'
      units_jl(k) = 'TENTHS OF METERS/SECOND'
      scale_jl(k) = .2E+1
      ia_jl(k) = ia_dga
      jgrid_jl(k) = 2
c
      k=k+1
      jl_vwpac = k
      sname_jl(k) = 'v_wpac' !'AJL45'
      lname_jl(k) = 'V WIND AVERAGED OVER I=WEST PACIFIC'
      units_jl(k) = 'TENTHS OF METERS/SECOND'
      scale_jl(k) = .2E+1
      ia_jl(k) = ia_dga
      jgrid_jl(k) = 2
c
      k=k+1
      jl_wwpac = k
      sname_jl(k) = 'vvel_wpac' !'AJL46'
      lname_jl(k) = 'VERTICAL VELOCITY FOR WEST PACIFIC'
      units_jl(k) = '10**-5 METERS/SECOND'
      scale_jl(k) = -1.D5*RGAS/(5.*GRAV)
      ia_jl(k) = ia_dga
      jgrid_jl(k) = 1
c
      k=k+1
      jl_47 = k
      sname_jl(k) = 'AJL47' !V-V*  =D((V-VI)*(T-TI)/DTHDP)/DP
      lname_jl(k) = 'unknown'
      units_jl(k) = 'unknown'
c
      k=k+1
      jl_zmfntmom = k
      sname_jl(k) = 'zmf_nt_mom' !'AJL48'
      lname_jl(k) = 'NORTH TRANS ZON. MOM. BY ZON. MEAN FLOW'
      units_jl(k) = 'unknown'
c
      k=k+1
      jl_totntmom = k
      sname_jl(k) = 'tot_nt_mom' !'AJL49'
      lname_jl(k) = 'TOTAL NORTH TRANS ZON. MOM.'
      units_jl(k) = 'unknown'
c
      k=k+1
      jl_mchphas = k
      sname_jl(k) = 'AJL50'
      lname_jl(k) = 'DT(MC)*P  CHANGE OF PHASE'
      units_jl(k) = '100 PA*K'
c
      k=k+1
      jl_mcdtotw = k
      sname_jl(k) = 'mc_del_tot_wat' !'AJL51' 'CLHE*DQ(MC BEFORE COND)*P'
      lname_jl(k) = 'CHANGE IN TOTAL WATER BY MOIST CONV'
      units_jl(k) = '100 PA*K'
c
      k=k+1
      jl_dudtsdrg = k
      sname_jl(k) = 'dudt_sdrag' !'AJL52'
      lname_jl(k) = 'DU/DT BY SDRAG'
      units_jl(k) = '10**-6 M S-2'
      scale_jl(k) = 1.D6/(FIM*DTsrc)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 2
c
      k=k+1
      jl_mcdlht = k
      sname_jl(k) = 'moist_lat' !'AJL53'
      lname_jl(k) = 'CHANGE OF LATENT HEAT BY MOIST CONV.'
      units_jl(k) = '10**14 W/UNIT SIGMA'
      scale_jl(k) = 100.D-14*SHA*TWOPI/(GRAV*DTsrc*DLON*FIM)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_trbke = k
      sname_jl(k) = 'tke' !'AJL54'
      lname_jl(k) = 'TURBULENT KINETIC ENERGY'
      units_jl(k) = 'W/M^2'
      scale_jl(k) = 1.
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_trbdlht = k
      sname_jl(k) = 'turb_lat' !'AJL55'
      lname_jl(k) = 'CHANGE OF LATENT HEAT BY TURBULENCE'
      units_jl(k) = '10**14 W/UNIT SIGMA'
      scale_jl(k) = 100.D-14*SHA*TWOPI/(GRAV*DTsrc*DLON*FIM)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_mcheat = k
      sname_jl(k) = 'tot_ht_mc' !'AJL56'
      lname_jl(k) = 'TOTAL HEATING BY MOIST CONVECTION (Q1)'
      units_jl(k) = '10**14 WATTS/DSIG'
      scale_jl(k) = 100.D-14*SHA*TWOPI/(GRAV*DTsrc*DLON*FIM)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      k=k+1
      jl_mcdry = k
      sname_jl(k) = 'tot_dry_mc' !'AJL57'
      lname_jl(k) = 'TOTAL DRYING BY MOIST CONVECTION (Q2)'
      units_jl(k) = '10**14 WATTS/DSIG'
      scale_jl(k) = 100.D-14*SHA*TWOPI/(GRAV*DTsrc*DLON*FIM)
      ia_jl(k) = ia_src
      jgrid_jl(k) = 1
c
      return
      end subroutine jl_defs

      subroutine sjl_defs
      use CONSTANT, only : grav,sday,sha
      use MODEL_COM, only : byim
      use DAGCOM
      implicit none
      integer :: k
c
      do k=1,kasjl
         write(name_sjl(k),'(a4,i3.3)') 'ASJL',k
         lname_sjl(k) = 'unused'
         units_sjl(k) = 'unused'
      enddo
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
      name_sjl(k) = 'srad_heat' !'ASJL03'
      lname_sjl(k) = 'SOLAR RADIATION HEATING RATE' !'SRHR'
      units_sjl(k) = 'HUNDREDTHS OF DEGREES KELVIN/DAY' !'W/m^2'
      scale_sjl(k) = 100.D-2*GRAV*SDAY*BYIM/SHA
      ia_sjl(k) = ia_rad
c
      k=k+1
      name_sjl(k) = 'trad_cool' !'ASJL04'
      lname_sjl(k) = 'THERMAL RADIATION COOLING RATE' !'TRHR'
      units_sjl(k) = 'HUNDREDTHS OF DEGREES K/DAY' !'W/m^2'
      scale_sjl(k) = -100.D-2*GRAV*SDAY*BYIM/SHA
      ia_sjl(k) = ia_rad
c
      return
      end subroutine sjl_defs

      subroutine jk_defs
      use CONSTANT, only : twopi
      use MODEL_COM, only : fim,byim
      use GEOM, only : dlon
      use DAGCOM
      implicit none
      integer :: k
c
      do k=1,kajkx
         write(sname_jk(k),'(a3,i3.3)') 'AJK',k
         lname_jk(k) = 'unused'
         units_jk(k) = 'unused'
      enddo
c
      k=0
c
      k=k+1
      jk_dpa = k             !'AJK01'
      sname_jk(k) = 'dp_cp1' !   DP=MIN(PM(K),PS)-MIN(PM(K+1),PS)
      lname_jk(k) =  'PRESSURE DIFFERENCES (CP,PT)' ! DP (PT GRID)
      units_jk(k) = 'mb'
      scale_jk(k) = byim
      ia_jk(k) = ia_dga
      jgrid_jk(k) = 1
c
      k=k+1
      jk_dpb = k
      sname_jk(k) = 'dp_cp2' !'AJK02'
      lname_jk(k) = 'PRESSURE DIFFERENCES (CP,UV)' ! DP (UV GRID)
      units_jk(k) = 'mb'
      scale_jk(k) = byim
      ia_jk(k) = ia_dga
      jgrid_jk(k) = 2
c
      k=k+1
      jk_temp = k
      sname_jk(k) = 'temp' !'AJK03'
      lname_jk(k) = 'TEMPERATURE' !'(TX-273.16)*DP'
      units_jk(k) = 'DEGREES CENTIGRADE'
      jgrid_jk(k) = 1
c
      k=k+1
      jk_hght = k
      sname_jk(k) = 'height' !'AJK04'
      lname_jk(k) = 'HEIGHT' !'PHI*DP'
      units_jk(k) = 'HUNDREDS OF METERS'
      jgrid_jk(k) = 1
c
      k=k+1
      jk_q = k
      sname_jk(k) = 'q' !'AJK05'
      lname_jk(k) = 'SPECIFIC HUMIDITY' !'Q*DP'
      units_jk(k) = '10**-5 KG H2O/KG AIR'
      jgrid_jk(k) = 1
c
      k=k+1
      jk_theta = k
      sname_jk(k) = 'pot_temp' !'AJK06'
      lname_jk(k) = 'POTENTIAL TEMPERATURE' !'TH*DP'
      units_jk(k) = 'DEGREES KELVIN'
      jgrid_jk(k) = 1
c
      k=k+1
      jk_rh = k
      sname_jk(k) = 'rh' !'AJK07'
      lname_jk(k) = 'RELATIVE HUMIDITY' !'RH*DP'
      units_jk(k) = 'PERCENT'
      jgrid_jk(k) = 1
c
      k=k+1
      jk_u = k
      sname_jk(k) = 'u' !'AJK08'
      lname_jk(k) = 'ZONAL WIND (U COMPONENT)' !'U*DP4  (UV GRID)'
      units_jk(k) = 'TENTHS OF METERS/SECOND' !'100 PA*m/s'
      jgrid_jk(k) = 2
c
      k=k+1
      jk_v = k
      sname_jk(k) = 'v' !'AJK09'
      lname_jk(k) = 'MERIDIONAL WIND (V COMPONENT)' !'V*DP4  (UV GRID)'
      units_jk(k) = 'HUNDREDTHS OF METERS/SECOND' !'100 PA*m/s'
      jgrid_jk(k) = 2
c
      k=k+1
      jk_zmfke = k
      sname_jk(k) = 'zon_ke' !'AJK10'
      lname_jk(k) = 'KINETIC ENERGY OF ZONAL MEAN FLOW'
      units_jk(k) = 'm2/s2'
c
      k=k+1
      jk_totke = k
      sname_jk(k) = 'tot_ke' !'AJK11'
      lname_jk(k) = 'TOTAL KINETIC ENERGY'
      units_jk(k) = '10**4 JOULES/M**2/UNIT SIGMA'
      jgrid_jk(k) = 2
c
      k=k+1
      jk_zmfntsh = k
      sname_jk(k) = 'zmf_nt_sh' !'AJK12'
      lname_jk(k) = 'NORTH. TRANS. SENS. HT. BY ZON. MEAN FLOW'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_totntsh = k
      sname_jk(k) = 'tot_nt_sh' !'AJK13'
      lname_jk(k) = 'TOT NORTH. TRANS. SENS. HT'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_zmfntgeo = k
      sname_jk(k) = 'zmf_nt_geo' !'AJK14'
      lname_jk(k) = 'NORTH. TRANS. GEOPOT. BY ZON. MEAN FLOW'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_totntgeo = k
      sname_jk(k) = 'tot_nt_geo' !'AJK15'
      lname_jk(k) = 'TOT NORTH. TRANS. GEOPOT.'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_zmfntlh = k
      sname_jk(k) = 'zmf_nt_lh' !'AJK16'
      lname_jk(k) = 'NORTH TRANS LAT HT BY ZON. MEAN FLOW'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_totntlh = k
      sname_jk(k) = 'tot_nt_lh' !'AJK17'
      lname_jk(k) = 'TOTAL NORTHWARD TRANSPORT OF LATENT HEAT'
      units_jk(k) = '10**14 WATTS/UNIT SIG'
      jgrid_jk(k) = 2
c
      k=k+1
      jk_zmfntke = k
      sname_jk(k) = 'zmf_nt_ke' !'AJK18'
      lname_jk(k) = 'NORTH TRANS KE BY ZON. MEAN FLOW'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_totntke = k
      sname_jk(k) = 'tot_nt_ke' !'AJK19'
      lname_jk(k) = 'TOTAL NORTHWARD TRANSPORT OF KINETIC ENERGY'
      units_jk(k) = '10**12 WATTS/DSIG'
      jgrid_jk(k) = 2
c
      k=k+1
      jk_zmfntmom = k
      sname_jk(k) = 'zmf_nt_mom' !'AJK20'
      lname_jk(k) = 'NORTH TRANS ZON. MOM. BY ZON. MEAN FLOW'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_totntmom = k
      sname_jk(k) = 'tot_nt_mom' !'AJK21'
      lname_jk(k) = 'TOT NORTH. TRANS. ZON. MOM.'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_p2kedpgf = k
      sname_jk(k) = 'p2k_eddy_pgf' !'AJK22'
      lname_jk(k) = 'P-K BY EDDY PRESSURE GRADIENT FORCE'
      units_jk(k) = '10**-1 W/M**2/UNIT SIGMA'
      jgrid_jk(k) = 2
c
      k=k+1
      jk_dpsqr = k
      sname_jk(k) = 'dp_sqr' !'AJK23'
      lname_jk(k) = 'SQUARE OF PRESSURE DIFFERENCES'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_nptsavg = k
      sname_jk(k) = 'npts_avg' !'AJK24'
      lname_jk(k) = 'NUMBER OF GRIDPOINTS INCLUDED IN AVERAGE (CP)'
      units_jk(k) = '1'
      scale_jk(k) = TWOPI/(DLON*FIM)
      ia_jk(k) = ia_dga
      jgrid_jk(k) = 2
c
      k=k+1
      jk_vvel = k
      sname_jk(k) = 'vvel' !'AJK25'
      lname_jk(k) = 'VERTICAL VELOCITY'
      units_jk(k) = '10**-5 MILLIBARS/SECOND'
      jgrid_jk(k) = 1
c
      k=k+1
      jk_zmfvtdse = k
      sname_jk(k) = 'zmf_vt_dse' !'AJK26'
      lname_jk(k) = 'VERT TRANS DRY STAT. ENER. BY ZON. MEAN FLOW'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_totvtdse = k
      sname_jk(k) = 'tot_vt_dse' !'AJK27'
      lname_jk(k)='TOTAL LGE SCALE VERT. TRANS. OF DRY STAT. ENER. (CP)'
      units_jk(k) = '10**14 WATTS'
      jgrid_jk(k) = 1
c
      k=k+1
      jk_zmfvtlh = k
      sname_jk(k) = 'zmf_vt_lh' !'AJK28'
      lname_jk(k) = 'VERT TRANS LATENT HEAT BY ZON. MEAN FLOW'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_totvtlh = k
      sname_jk(k) = 'tot_vt_lh' !'AJK29'
      lname_jk(k) = 'TOTAL LARGE SCALE VERT. TRANS. OF LATENT HEAT (CP)'
      units_jk(k) = '10**13 WATTS'
      jgrid_jk(k) = 1
c
      k=k+1
      jk_vtgeoeddy = k
      sname_jk(k) = 'vt_geopot_eddy' !'AJK30'
      lname_jk(k) = 'VERT. TRANS. OF GEOPOTENTIAL ENERGY BY EDDIES (CP)'
      units_jk(k) = '10**12 WATTS'
      jgrid_jk(k) = 1
c
      k=k+1
      jk_barekegen = k
      sname_jk(k) = 'baroc_eddy_ke_gen' !'AJK31'
      lname_jk(k) = 'BAROCLINIC EDDY KINETIC ENERGY GEN.'
      units_jk(k) = '10**-1 WATTS/M**2/SIGMA'
      jgrid_jk(k) = 1
c
      k=k+1
      jk_potvort = k
      sname_jk(k) = 'pot_vort' !'AJK32'
      lname_jk(k) = 'POTENTIAL VORTICITY (CP)'
      units_jk(k) = '10**-6 K/(MB-S)'
      jgrid_jk(k) = 1
c
      k=k+1
      jk_vtpv = k
      sname_jk(k) = 'vt_pv' !'AJK33'
      lname_jk(k) = 'VERT. TRANSPORT OF POTENTIAL VORTICITY (CP)'
      units_jk(k) = '10**4 KG-DEG K/MB/S/S'
      jgrid_jk(k) = 1
c
      k=k+1
      jk_vtpveddy = k
      sname_jk(k) = 'vt_pv_eddy' !'AJK34'
      lname_jk(k) = 'VERT. TRANS. OF POT. VORT. BY EDDIES (CP)'
      units_jk(k) = '10**4 KG-DEG K/MB/S/S'
      jgrid_jk(k) = 1
c
      k=k+1
      jk_nptsavg1 = k
      sname_jk(k) = 'npts_avg1' !'AJK35'
      lname_jk(k) = 'NUMBER OF GRIDPOINTS IN AVERAGE (PT GRID)'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_totvtke = k
      sname_jk(k) = 'tot_vt_ke' !'AJK36'
      lname_jk(k) ='TOTAL LGE SCALE VERT. TRANS. OF KINETIC ENERGY (CP)'
      units_jk(k) = '10**11 WATTS'
      jgrid_jk(k) = 2
c
      k=k+1
      jk_vtameddy = k
      sname_jk(k) = 'vt_am_eddy' !'AJK37'
      lname_jk(k) = 'VERT. TRANS. OF ANG. MOMENTUM BY EDDIES (CP)'
      units_jk(k) = '10**16 JOULES'
      jgrid_jk(k) = 2
c
      k=k+1
      jk_totvtam = k
      sname_jk(k) = 'tot_vt_am' !'AJK38'
      lname_jk(k) = 'TOTAL LGE SCALE VERT. TRANS. OF ANG. MOMENTUM (CP)'
      units_jk(k) = '10**18 JOULES'
      jgrid_jk(k) = 2
c
      k=k+1
      jk_sheth = k
      sname_jk(k) = 'sheth' !'AJK39'
      lname_jk(k) = 'unknown'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_dudtmadv = k
      sname_jk(k) = 'dudt_mean_advec' !'AJK40'
      lname_jk(k) = 'DU/DT BY MEAN ADVECTION (CP)'
      units_jk(k) = '10**-6 M/S/S'
      jgrid_jk(k) = 1
c
      k=k+1
      jk_dtdtmadv = k
      sname_jk(k) = 'dtempdt_mean_advec' !'AJK41'
      lname_jk(k) = 'DTEMP/DT BY MEAN ADVECTION (CP)'
      units_jk(k) = '10**-1 DEG-K/DAY'
      jgrid_jk(k) = 1
c
      k=k+1
      jk_dudttem = k
      sname_jk(k) = 'dudt_advec_tem' !'AJK42'
      lname_jk(k) = 'DU/DT BY TRANSFORMED ADVECTION (CP)'
      units_jk(k) = '10**-6 M/S/S'
      jgrid_jk(k) = 1
c
      k=k+1
      jk_dtdttem = k
      sname_jk(k) = 'dtempdt_advec_tem' !'AJK43'
      lname_jk(k) = 'DTEMP/DT BY TRANSFORMED ADVECTION (CP)'
      units_jk(k) = '10**-1 DEG-K/DAY'
      jgrid_jk(k) = 1
c
      k=k+1
      jk_epflxncp = k
      sname_jk(k) = 'epflx_north_cp' !'AJK44'
      lname_jk(k) = 'NORTHWARD COMP. OF ELIASSEN-PALM FLUX (CP)'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_epflxvcp = k
      sname_jk(k) = 'epflx_vert_cp' !'AJK45'
      lname_jk(k) = 'NORTHWARD COMP. OF ELIASSEN-PALM FLUX (CP)'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_uinst = k
      sname_jk(k) = 'u_inst' !'AJK46'
      lname_jk(k) = 'INSTANTANEOUS ZONAL AVERAGE OF ZONAL WIND'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_totdudt = k
      sname_jk(k) = 'tot_dudt' !'AJK47'
      lname_jk(k) = 'DU/DT   TOTAL CHANGE (CP)'
      units_jk(k) = '10**-6 M/S/S'
      jgrid_jk(k) = 2
c
      k=k+1
      jk_tinst = k
      sname_jk(k) = 't_inst' !'AJK48'
      lname_jk(k) = 'INSTANTANEOUS ZONAL AVERAGE OF TEMPERATURE'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_totdtdt = k
      sname_jk(k) = 'dtempdt' !'AJK49'
      lname_jk(k) = 'DTEMP/DT   TOTAL CHANGE (CP)'
      units_jk(k) = '10**-1 DEG-K/DAY'
      jgrid_jk(k) = 1
c
      k=k+1
      jk_eddvtpt = k
      sname_jk(k) = 'edd_vt_pt'
      lname_jk(k) = 'EDDY VERTICAL TRANSPORT OF POT. TEMP.'
      units_jk(k) = 'unknown'
c
      k=k+1
      jk_cldh2o = k
      sname_jk(k) = 'cldh2o' !'AJK51'
      lname_jk(k) = 'TOTAL CLOUD WATER CONTENT'
      units_jk(k) = '10**-6 KG/KG'
      jgrid_jk(k) = 1
c
      return
      end subroutine jk_defs

      subroutine ijk_defs
      use CONSTANT, only : bygrav,tf
      use DAGCOM
      implicit none
      integer :: k
c
      do k=1,kaijkx
         write(name_ijk(k),'(a4,i3.3)') 'AIJK',k
         lname_ijk(k) = 'unused'
         units_ijk(k) = 'unused'
         scale_ijk(k) = 1.
         off_ijk(k)   = 0.
      enddo
c
      k=0
c
      k=k+1
      IJK_U=k
      name_ijk(k) = 'u' !'UDPB'
      lname_ijk(k) = 'U-WIND            x delta p, b-grid'
      units_ijk(k) = 'm/s'
      scale_ijk(k) = 1.
      off_ijk(k)   = 0.
c
      k=k+1
      IJK_V=k
      name_ijk(k) = 'v' !'VDPB'
      lname_ijk(k) = 'V-WIND            x delta p, b-grid'
      units_ijk(k) = 'm/s'
      scale_ijk(k) = 1.
      off_ijk(k)   = 0.
c
      k=k+1
      IJK_DSE=k
      name_ijk(k) = 'dse' !'DSEDPB'
      lname_ijk(k) = 'DRY STAT. ENERGY  x delta p x 4, b-grid'
      units_ijk(k) = 'm^2/s^2'
      scale_ijk(k) = .25
      off_ijk(k)   = 0.
c
      k=k+1
      IJK_DP=k
      name_ijk(k) = 'dp' !'DPB'
      lname_ijk(k) = 'DELTA-P           b-grid'
      units_ijk(k) = '100 PA'
      scale_ijk(k) = 1.
      off_ijk(k)   = 0.
c
      k=k+1
      IJK_T=k
      name_ijk(k) = 't' !'TDPB'
      lname_ijk(k) = 'TEMPERATURE       x delta p x 4, b-grid'
      units_ijk(k) = 'C'
      scale_ijk(k) = 0.25
      off_ijk(k)   = -TF
c
      k=k+1
      IJK_Q=k
      name_ijk(k) = 'q' !'QDPB'
      lname_ijk(k) = 'SPECIFIC HUMIDITY x delta p x 4, b-grid'
      units_ijk(k) = '10**-5'
      scale_ijk(k) = 0.25*1d5
      off_ijk(k)   = 0.
c
      return
      end subroutine ijk_defs


      subroutine wave_defs
      use DAGCOM
      implicit none
      integer :: k
c
      do k=1,kwp
         write(name_wave(k),'(a5,i3.3)') 'WAVEP',k
         lname_wave(k) = 'unused'
         units_wave(k) = 'unused'
      enddo
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


      subroutine pj_defs
      use DAGCOM
      implicit none
      integer :: k
c
      do k=1,kapj
         write(name_pj(k),'(a3,i3.3)') 'APJ',k
         lname_pj(k) = 'unused'
         units_pj(k) = 'unused'
      enddo
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
      do k=1,ktsf
         write(name_tsf(k),'(a6,i3.3)') 'TSFREZ',k
         lname_tsf(k) = 'unused'
         units_tsf(k) = 'unused'
      enddo
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

      subroutine diurn_defs
!@sum  diurn_defs definitions for diurnal diagnostic accumulated arrays
!@auth G. Schmidt
!@ver  1.0
      use CONSTANT, only : sha,rgas,twopi,sday,grav
      use MODEL_COM, only : dtsrc,nisurf
      use DAGCOM
      implicit none
      integer :: k
c
      do k=1,NDIUVAR
         write(name_dd(k),'(a5,i3.3)') 'DIURN',k
         lname_dd(k) = 'unused'
         units_dd(k) = 'unused'
         scale_dd(k) = 1.
      enddo
c
      k=0
c
      k=k+1
      IDD_ISW=k
      name_dd(k)='INCSW'
      units_dd(k)='W/M^2'
      scale_dd(k)=1.
      lname_dd(k)=' INC SW RADIATION'
c
      k=k+1
      IDD_PALB=k
      name_dd(k)='PLALB'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' P ALBD '
c
      k=k+1
      IDD_GALB=k
      name_dd(k)='GRDALB'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' G ALBD '
c
      k=k+1
      IDD_ABSA=k
      name_dd(k)='ABSATM'
      units_dd(k)='W/M^2'
      scale_dd(k)=1.
      lname_dd(k)=' ABS ATM'
c
      k=k+1
      IDD_ECND=k
      name_dd(k)='ENRGCND'
      units_dd(k)='W/M^2'
      scale_dd(k)=100.*SHA/(GRAV*DTsrc)
      lname_dd(k)=' E CNDS '
c
      k=k+1
      IDD_SPR=k
      name_dd(k)='SRFPRS'
      units_dd(k)='mb'
      scale_dd(k)=1.
      lname_dd(k)=' SRF PRS'
c
      k=k+1
      IDD_PT5=k
      name_dd(k)='POTT5'
      units_dd(k)='K'
      scale_dd(k)=1.
      lname_dd(k)=' PT 5   '
c
      k=k+1
      IDD_PT4=k
      name_dd(k)='POTT4'
      units_dd(k)='K'
      scale_dd(k)=1.
      lname_dd(k)=' PT 4   '
c
      k=k+1
      IDD_PT3=k
      name_dd(k)='POTT3'
      units_dd(k)='K'
      scale_dd(k)=1.
      lname_dd(k)=' PT 3   '
c
      k=k+1
      IDD_PT2=k
      name_dd(k)='POTT2'
      units_dd(k)='K'
      scale_dd(k)=1.
      lname_dd(k)=' PT 2   '
c
      k=k+1
      IDD_PT1=k
      name_dd(k)='POTT1'
      units_dd(k)='K'
      scale_dd(k)=1.
      lname_dd(k)=' PT 1   '
c
      k=k+1
      IDD_TS=k
      name_dd(k)='TSURF'
      units_dd(k)='K'
      scale_dd(k)=1.
      lname_dd(k)=' TS     '
c
      k=k+1
      IDD_TG1=k
      name_dd(k)='TGRND'
      units_dd(k)='K'
      scale_dd(k)=1.
      lname_dd(k)=' TG1    '
c
      k=k+1
      IDD_Q5=k
      name_dd(k)='Q5'
      units_dd(k)='1d-5 kg/kg'
      scale_dd(k)=1d5
      lname_dd(k)=' Q 5    '
c
      k=k+1
      IDD_Q4=k
      name_dd(k)='Q4'
      units_dd(k)='1d-5 kg/kg'
      scale_dd(k)=1d5
      lname_dd(k)=' Q 4    '
c
      k=k+1
      IDD_Q3=k
      name_dd(k)='Q3'
      units_dd(k)='1d-5 kg/kg'
      scale_dd(k)=1d5
      lname_dd(k)=' Q 3    '
c
      k=k+1
      IDD_Q2=k
      name_dd(k)='Q2'
      units_dd(k)='1d-5 kg/kg'
      scale_dd(k)=1d5
      lname_dd(k)=' Q 2    '
c
      k=k+1
      IDD_Q1=k
      name_dd(k)='Q1'
      units_dd(k)='1d-5 kg/kg'
      scale_dd(k)=1d5
      lname_dd(k)=' Q 1    '
c
      k=k+1
      IDD_QS=k
      name_dd(k)='QSURF'
      units_dd(k)='1d-5 kg/kg'
      scale_dd(k)=1d5
      lname_dd(k)=' QS     '
c
      k=k+1
      IDD_QG=k
      name_dd(k)='QGRND'
      units_dd(k)='1d-5 kg/kg'
      scale_dd(k)=1d5
      lname_dd(k)=' QG     '
c
      k=k+1
c      IDD_=k
c      name_dd(k)=''
c      units_dd(k)=''
c      scale_dd(k)=
c      lname_dd(k)='        '
c
      k=k+1
c      IDD_=k
c      name_dd(k)=''
c      units_dd(k)=''
c      scale_dd(k)=''
c      lname_dd(k)='        '
c
      k=k+1
c      IDD_=k
c      name_dd(k)=''
c      units_dd(k)=''
c      scale_dd(k)=''
c      lname_dd(k)='        '
c
      k=k+1
c      IDD_=k
c      name_dd(k)=''
c      units_dd(k)=''
c      scale_dd(k)=''
c      lname_dd(k)='        '
c
      k=k+1
c      IDD_=k
c      name_dd(k)=''
c      units_dd(k)=''
c      scale_dd(k)=''
c      lname_dd(k)='        '
c
      k=k+1
c      IDD_=k
c      name_dd(k)=''
c      units_dd(k)=''
c      scale_dd(k)=''
c      lname_dd(k)='        '
c
      k=k+1
c      IDD_=k
c      name_dd(k)=''
c      units_dd(k)=''
c      scale_dd(k)=''
c      lname_dd(k)='        '
c
      k=k+1
      IDD_SWG=k
      name_dd(k)='SWGRND'
      units_dd(k)='W/M^2'
      scale_dd(k)=NISURF/DTsrc
      lname_dd(k)=' SW ON G'
c
      k=k+1
      IDD_LWG=k
      name_dd(k)='LWGRND'
      units_dd(k)='W/M^2'
      scale_dd(k)=NISURF/DTsrc
      lname_dd(k)=' LW AT G'
c
      k=k+1
      IDD_SH=k
      name_dd(k)='SENSHT'
      units_dd(k)='W/M^2'
      scale_dd(k)=NISURF/DTsrc
      lname_dd(k)=' SNSB HT'
c
      k=k+1
      IDD_LH=k
      name_dd(k)='LATHT'
      units_dd(k)='W/M^2'
      scale_dd(k)=NISURF/DTsrc
      lname_dd(k)=' LAT HT '
c
      k=k+1
      IDD_HZ0=k
      name_dd(k)='NETHT'
      units_dd(k)='W/M^2'
      scale_dd(k)=NISURF/DTsrc
      lname_dd(k)=' HEAT Z0'
c
      k=k+1
      IDD_UG=k
      name_dd(k)='UGEOS'
      units_dd(k)='0.1 m/s'
      scale_dd(k)=10.
      lname_dd(k)=' UG*10  '
c
      k=k+1
      IDD_VG=k
      name_dd(k)='VGEOS'
      units_dd(k)='0.1 m/s'
      scale_dd(k)=10.
      lname_dd(k)=' VG*10  '
c
      k=k+1
      IDD_WG=k
      name_dd(k)='WGEOS'
      units_dd(k)='0.1 m/s'
      scale_dd(k)=10.
      lname_dd(k)=' WG*10  '
c
      k=k+1
      IDD_US=k
      name_dd(k)='USURF'
      units_dd(k)='0.1 m/s'
      scale_dd(k)=10.
      lname_dd(k)=' US*10  '
c
      k=k+1
      IDD_VS=k
      name_dd(k)='VSURF'
      units_dd(k)='m/s'
      scale_dd(k)=10.
      lname_dd(k)=' VS*10  '
c
      k=k+1
      IDD_WS=k
      name_dd(k)='WSURF'
      units_dd(k)='m/s'
      scale_dd(k)=10.
      lname_dd(k)=' WS*10  '
c
      k=k+1
      IDD_CIA=k
      name_dd(k)='CRISOANG'
      units_dd(k)=''
      scale_dd(k)=360./TWOPI
      lname_dd(k)=' ALPHA0 '
c
      k=k+1
      IDD_RIS=k
      name_dd(k)='RIS1'
      units_dd(k)='0.01'
      scale_dd(k)=100.
      lname_dd(k)=' RIS1*E2'
c
      k=k+1
      IDD_RIG=k
      name_dd(k)='RIGS'
      units_dd(k)='0.01'
      scale_dd(k)=100.
      lname_dd(k)=' RIGS*E2'
c
      k=k+1
      IDD_CM=k
      name_dd(k)='CM'
      units_dd(k)='1d-4'
      scale_dd(k)=1d4
      lname_dd(k)=' CM*E4  '
c
      k=k+1
      IDD_CH=k
      name_dd(k)='CH'
      units_dd(k)='1d-4'
      scale_dd(k)=1d4
      lname_dd(k)=' CH*E4  '
c
      k=k+1
      IDD_CQ=k
      name_dd(k)='CQ'
      units_dd(k)='1d-4'
      scale_dd(k)=1d4
      lname_dd(k)=' CQ*E4  '
c
      k=k+1
      IDD_EDS=k
      name_dd(k)='EDS'
      units_dd(k)='0.1'
      scale_dd(k)=10.
      lname_dd(k)=' EDS1*10'
c
      k=k+1
      IDD_DBL=k
      name_dd(k)='DBL'
      units_dd(k)='m'
      scale_dd(k)=1.
      lname_dd(k)=' DBL    '
c
      k=k+1
      IDD_DCF=k
      name_dd(k)='DCFREQ'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' DC FREQ'
c
      k=k+1
      IDD_LDC=k
      name_dd(k)='LDC'
      units_dd(k)='0.1 L'
      scale_dd(k)=10.
      lname_dd(k)=' LDC*10 '
c
      k=k+1
      IDD_PR=k
      name_dd(k)='PREC'
      units_dd(k)='0.01 mm/day'
      scale_dd(k)=100.*100.*SDAY/(DTsrc*GRAV)
      lname_dd(k)=' PRC*100'  ! check scale
c
      k=k+1
      IDD_EV=k
      name_dd(k)='EVAP'
      units_dd(k)='0.01 mm/day'
      scale_dd(k)=100.*SDAY*NISURF/DTsrc
      lname_dd(k)=' EVP*100'
c
      k=k+1
      IDD_DMC=k
      name_dd(k)='DEEPMC'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' DEEP MC'
c
      k=k+1
      IDD_SMC=k
      name_dd(k)='SHLWMC'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' SHLW MC'
c
      k=k+1
      IDD_CL7=k
      name_dd(k)='CLD7'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' CLD 7  '
c
      k=k+1
      IDD_CL6=k
      name_dd(k)='CLD6'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' CLD 6  '
c
      k=k+1
      IDD_CL5=k
      name_dd(k)='CLD5'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' CLD 5  '
c
      k=k+1
      IDD_CL4=k
      name_dd(k)='CLD4'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' CLD 4  '
c
      k=k+1
      IDD_CL3=k
      name_dd(k)='CLD3'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' CLD 3  '
c
      k=k+1
      IDD_CL2=k
      name_dd(k)='CLD2'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' CLD 2  '
c
      k=k+1
      IDD_CL1=k
      name_dd(k)='CLD1'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' CLD 1  '
c
      k=k+1
      IDD_W=k
      name_dd(k)='VERTVEL'
      units_dd(k)='1d-5 m/s'
      scale_dd(k)=1.
      lname_dd(k)=' W TO-5 '
c
      k=k+1
      IDD_CCV=k
      name_dd(k)='CLDCOV'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' C COVER'
c
      k=k+1
      IDD_SSP=k
      name_dd(k)='SSPREC'
      units_dd(k)='0.01 mm/day'
      scale_dd(k)=100.*100.*SDAY/(DTsrc*GRAV)
      lname_dd(k)=' SSP*100'
c
      k=k+1
      IDD_MCP=k
      name_dd(k)='MCPREC'
      units_dd(k)='0.01 mm/day'
      scale_dd(k)=100.*100.*SDAY/(DTsrc*GRAV)
      lname_dd(k)=' MCP*100'

      return
      end subroutine diurn_defs
