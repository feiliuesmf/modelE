      MODULE BDJ
!@sum  stores information for outputting latitude diagnostics
!@auth M. Kelley

      IMPLICIT NONE

!@param nj_out number of j-format output fields = (kd1m+10)
      integer, parameter :: nj_out=79

!@var units string containing output field units
      CHARACTER(LEN=50), DIMENSION(nj_out) :: UNITS
!@var lname string describing output field
      CHARACTER(LEN=50), DIMENSION(nj_out) :: LNAME
!@var sname string referencing output field in self-desc. output file
      CHARACTER(LEN=30), DIMENSION(nj_out) :: SNAME

!@var nt_j expanded version of diagj: ndex that includes albedo-indices
      integer, dimension(nj_out) :: nt_j
!@var nstype_out set to diagj: ntype+1 since indexing in diagj is 0:ntype
!@var iotype current output surf type (1 through nstype_out)
      integer :: nstype_out,iotype
      END MODULE BDJ

      SUBROUTINE J_TITLES
      USE BDJ
      IMPLICIT NONE
      INTEGER :: K
c
      k = 0
c
      k = k + 1
      sname(k) = 'inc_sw'
      lname(k) = 'SOLAR RADIATION INCIDENT ON PLANET'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'sw_abs_p0'
      lname(k) = ' SOLAR RADIATION ABSORBED BY PLANET'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'sw_abs_p1'
      lname(k) = 'SOLAR RADIATION ABSORBED BELOW PTOP'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'sw_abs_atm'
      lname(k) = 'SOLAR RADIATION ABSORBED BY ATMOSPHERE'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'sw_inc_z0'
      lname(k) = 'SOLAR RADIATION INCIDENT ON GROUND'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'sw_abs_z0'
      lname(k) = 'SOLAR RADIATION ABSORBED BY GROUND'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'net_lw_p0'
      lname(k) = 'THERMAL RADIATION EMITTED BY PLANET'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'net_lw_p1'
      lname(k) = 'NET THERMAL RADIATION AT PTOP'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'net_lw_z0'
      lname(k) = 'NET THERMAL RADIATION AT GROUND'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'net_rad_p0'
      lname(k) = 'NET RADIATION OF PLANET'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'net_rad_p1'
      lname(k) = 'NET RADIATION BELOW PTOP'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'net_rad_z0'
      lname(k) = 'NET RADIATION ABSORBED BY GROUND'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'net_clr_toa'
      lname(k) = 'NET CLEAR SKY RADIATION AT P0'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'net_clr_trp'
      lname(k) = 'NET CLEAR SKY RADIATION AT TROPOPAUSE (WMO)'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'net_tot_trp'
      lname(k) = 'NET RADIATION AT TROPOPAUSE (WMO)'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'snsht_flx'
      lname(k) = 'SENSIBLE HEAT FLUX INTO THE GROUND'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'evht_flx'
      lname(k) = 'EVAPORATION HEAT FLUX INTO THE GROUND'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'ht_cond_z1z2'
      lname(k) = 'CONDUCTION AT BOTTOM OF GROUND LAYER 2'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'net_ht_z1'
      lname(k) = 'NET HEAT FLUX BETWEEN GROUND LAYERS'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'tg2'
      lname(k) = 'TEMPERATURE OF GROUND LAYER 2'
      units(k) = '.1 C'
      k = k + 1
      sname(k) = 'tg1'
      lname(k) = 'TEMPERATURE OF GROUND LAYER 1'
      units(k) = '.1 C'
      k = k + 1
      sname(k) = 'evap'
      lname(k) = 'EVAPORATION'
      units(k) = 'MM/DAY'
      k = k + 1
      sname(k) = 'prec'
      lname(k) = 'PRECIPITATION'
      units(k) = 'MM/DAY'
      k = k + 1
      sname(k) = 'tair'
      lname(k) = 'AIR TEMPERATURE'
      units(k) = '.1 C'
      k = k + 1
      sname(k) = 't1'
      lname(k) = 'TEMPERATURE OF AIR LAYER 1'
      units(k) = '.1 C'
      k = k + 1
      sname(k) = 'tsurf'
      lname(k) = 'SURFACE AIR TEMPERATURE'
      units(k) = '.1 C'
      k = k + 1
      sname(k) = 'sstab_strat'
      lname(k) = 'STRATOSPHERIC STATIC STABILITY'
      units(k) = 'C/km'
      k = k + 1
      sname(k) = 'sstab_trop'
      lname(k) = 'TROPOSPHERIC STATIC STABILITY'
      units(k) = 'C/km'
      k = k + 1
      sname(k) = 'rich_num_strat'
      lname(k) = 'STRATOSPHERIC RICHARDSON NUMBER'
      units(k) = '1'
      k = k + 1
      sname(k) = 'rich_num_trop'
      lname(k) = 'TROPOSPHERIC RICHARDSON NUMBER'
      units(k) = '1'
      k = k + 1
      sname(k) = 'ross_num_strat'
      lname(k) = 'STRATOSPHERIC ROSSBY NUMBER'
      units(k) = '1'
      k = k + 1
      sname(k) = 'ross_num_trop'
      lname(k) = 'TROPOSPHERIC ROSSBY NUMBER'
      units(k) = '1'
      k = k + 1
      sname(k) = 'ocn_lak_ice_frac'
      lname(k) = 'OCEAN/LAKE ICE COVER'
      units(k) = '%'
      k = k + 1
      sname(k) = 'snow_cover'
      lname(k) = 'SNOW COVER'
      units(k) = '%'
c      k = k + 1
c      sname(k) = 'sw_correc'
c      lname(k) = 'SOLAR RADIATION CORRECTION'
c      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'ocn_ht_trans'
      lname(k) = 'CONVERGED OCEAN HEAT TRANSPORT'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'tg3'
      lname(k) = 'OCEAN TEMP AT MAXIMUM MIXED LAYER DEPTH'
      units(k) = '.1 C'
      k = k + 1
      sname(k) = 'dtdlat_strat'
      lname(k) = 'STRATO TEMP CHANGE PER DEGREE LATITUDE'
      units(k) = 'deg C/deg lat'
      k = k + 1
      sname(k) = 'dtdlat_trop'
      lname(k) = 'TROPO TEMP CHANGE PER DEGREE LATITUDE'
      units(k) = 'deg C/deg lat'
      k = k + 1
      sname(k) = 'ross_radius_strat'
      lname(k) = 'ROSSBY RADIUS IN THE STRATOSPHERE'
      units(k) = '10**5 m'
      k = k + 1
      sname(k) = 'ross_radius_trop'
      lname(k) = 'ROSSBY RADIUS IN THE TROPOSPHERE'
      units(k) = '10**5 m'
      k = k + 1
      sname(k) = 'prec_ht_flx'
      lname(k) = 'PRECIPITATION HEAT FLUX INTO THE GROUND'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'ht_runoff_z0'
      lname(k) = 'HEAT RUNOFF Z0'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'ht_wtr_difs_z1'
      lname(k) = 'HEAT DIFFUSION OF WATER OR ICE AT -Z1'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'ht_cond_z1'
      lname(k) = 'CONDUCTION AT BOTTOM OF GROUND LAYER 1'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'ht_ice_z1z2'
      lname(k) = 'ENERGY OF MELTING (OR TRANS) AT -Z1-Z2'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'nt_ht_z0'
      lname(k) = 'NET HEATING AT GROUND SURFACE'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'h2o_difs_z1'
      lname(k) = 'WATER OR ICE DIFFUSION AT -Z1'
      units(k) = 'mm/day'
      k = k + 1
      sname(k) = 'ice_thru_z1z2'
      lname(k) = 'ICE MELTING (OR TRANSPORT) AT -Z1-Z2'
      units(k) = 'mm H2O/day'
      k = k + 1
      sname(k) = 'h2o_runoff_mld'
      lname(k) = 'WATER RUNOFF THROUGH MIXED LAYER DEPTH'
      units(k) = 'mm/day'
      k = k + 1
      sname(k) = 'ht_runoff_mld'
      lname(k) = 'HEAT RUNOFF THROUGH MIXED LAYER DEPTH'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'wat_g1'
      lname(k) = 'WATER IN GROUND LAYER 1'
      units(k) = 'kg/m**2'
      k = k + 1
      sname(k) = 'ice_g1'
      lname(k) = 'ICE IN GROUND LAYER 1'
      units(k) = 'kg/m**2'
      k = k + 1
      sname(k) = 'wat_g2'
      lname(k) = 'WATER IN GROUND LAYER 2'
      units(k) = 'kg/m**2'
      k = k + 1
      sname(k) = 'ice_g2'
      lname(k) = 'ICE IN GROUND LAYER 2'
      units(k) = 'kg/m**2'
      k = k + 1
      sname(k) = 'snowdp'
      lname(k) = 'SNOW DEPTH'
      units(k) = 'kg/m**2'
      k = k + 1
      sname(k) = 'wat_runoff_z0'
      lname(k) = 'WATER RUNOFF AT GROUND SURFACE'
      units(k) = 'mm/day'
      k = k + 1
      sname(k) = 'btemp_window'
      lname(k) = 'BRIGHTNESS TEMP THROUGH WINDOW REGION'
      units(k) = 'deg C'
      k = k + 1
      sname(k) = 'net_ht_z1z2'
      lname(k) = 'NET HEAT FLUX AT BOTTOM OF GRND LAY 2'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'sscld'
      lname(k) = 'SUPER SATURATION CLOUD COVER'
      units(k) = '%'
      k = k + 1
      sname(k) = 'mccld'
      lname(k) = 'MOIST CONVECTIVE CLOUD COVER'
      units(k) = '%'
      k = k + 1
      sname(k) = 'totcld'
      lname(k) = 'TOTAL CLOUD COVER'
      units(k) = '%'
      k = k + 1
      sname(k) = 'mc_clddp'
      lname(k) = 'MC CLD DPTH'
      units(k) = 'MB'
      k = k + 1
      sname(k) = 'ssprec'
      lname(k) = 'SUPER SATURATION PRECIPITATION'
      units(k) = 'MM/DAY'
      k = k + 1
      sname(k) = 'mcprec'
      lname(k) = 'MOIST CONVECTIVE PRECIPITATION'
      units(k) = 'MM/DAY'
      k = k + 1
      sname(k) = 'atmh2o'
      lname(k) = 'WATER CONTENT OF ATMOSPHERE'
      units(k) = 'MM'
      k = k + 1
      sname(k) = 'lapse_rate'
      lname(k) = 'MEAN LAPSE RATE'
      units(k) = 'K/KM'
      k = k + 1
      sname(k) = 'lapse_rate_m'
      lname(k) = 'MOIST ADIABATIC LAPSE RATE'
      units(k) = 'K/KM'
      k = k + 1
      sname(k) = 'lapse_rate_c'
      lname(k) = 'GAMC'
      units(k) = 'K/KM'
      k = k + 1
      sname(k) = 'lw_inc_z0'
      lname(k) = 'THERMAL RADIATION INCIDENT ON GROUND'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'ht_thermocline'
      lname(k) = 'ENERGY DIFFUSION INTO THE THERMOCLINE'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'surf_type_frac'
      lname(k) = 'SURF TYPE FRACT'
      units(k) = '%'
      k = k + 1
      sname(k) = 'plan_alb'
      lname(k) = 'PLANETARY ALBEDO'
      units(k) = '%'
      k = k + 1
      sname(k) = 'plan_alb_vis'
      lname(k) = 'PLANETARY ALBEDO IN VISUAL'
      units(k) = '%'
      k = k + 1
      sname(k) = 'plan_alb_nir'
      lname(k) = 'PLANETARY ALBEDO IN NEAR IR'
      units(k) = '%'
      k = k + 1
      sname(k) = 'surf_alb'
      lname(k) = 'GROUND ALBEDO'
      units(k) = '%'
      k = k + 1
      sname(k) = 'surf_alb_vis'
      lname(k) = 'GROUND ALBEDO IN VISUAL'
      units(k) = '%'
      k = k + 1
      sname(k) = 'surf_alb_nir'
      lname(k) = 'GROUND ALBEDO IN NEAR IR'
      units(k) = '%'
      k = k + 1
      sname(k) = 'atm_alb_vis'
      lname(k) = 'ATMOSPHERIC ALBEDO IN VISUAL'
      units(k) = '%'
      k = k + 1
      sname(k) = 'atm_alb_nir'
      lname(k) = 'ATMOSPHERIC ALBEDO IN NEAR IR'
      units(k) = '%'
      k = k + 1
      sname(k) = 'atm_abs_vis'
      lname(k) = 'ATMOSPHERIC ABSORPTION IN VISUAL'
      units(k) = 'WT/M**2'
      k = k + 1
      sname(k) = 'atm_abs_nir'
      lname(k) = 'ATMOSPHERIC ABSORPTION IN NEAR IR'
      units(k) = 'WT/M**2'

      RETURN
      END SUBROUTINE J_TITLES
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
C***7   7  THERMAL RADIATION EMITTED BY PLANET (W/M**2)     
C***8   8  THERMAL RADIATION AT PTOP (W/M**2)          
C***9   9  THERMAL RADIATION EMITTED BY GROUND (W/M**2)
C**10  67  THERMAL RADIATION INCIDENT ON GROUND (W/M**2)
C****
C**11  55  BRIGHTNESS TEMPERATURE THROUGH WINDOW REGION (K-TF)
C****  10  NET RADIATION ABSORBED BY PLANET (W/M**2)
C****  11  NET RADIATION ABSORBED BELOW PTOP (W/M**2)
C****  12  NET RADIATION ABSORBED BY GROUND (W/M**2)
C****  81  NET CLEAR SKY RADIATION AT P0 W/M**2
C****  82  NET CLEAR SKY RADIATION AT TROPOPAUSE (WMO) W/M**2
C****  83  NET RADIATION AT TROPOPAUSE (WMO) W/M**2
C****  13  SENSIBLE HEAT FLUX INTO THE GROUND (W/M**2)
C****  14  EVAPORATION HEAT FLUX INTO THE GROUND (W/M**2)
C****  39  PRECIPITATION HEAT FLUX INTO THE GROUND (W/M**2)
C****
C**21  40  HEAT RUNOFF FROM FIRST GROUND LAYER (W/M**2)
C****  44  NET HEATING AT Z0 (W/M**2)
C****  42  CONDUCTION AT -Z1 (W/M**2)
C****  41  HEAT OF WATER OR ICE DUFFUSION AT -Z1 (W/M**2)
C****  16  NET HEATING AT -Z1 (W/M**2)
C****  15  CONDUCTION AT -Z1-Z2 (W/M**2)
C****  43  ENERGY OF ICE MELTING (OR TRANSPORTING) AT -Z1-Z2 (W/M**2)
C****  56  NET HEATING AT -Z1-Z2 (W/M**2)
C****  33  OCEAN TRANSPORT (W/M**2)
C****  48  HEAT RUNOFF THROUGH THE MIXED LAYER DEPTH (W/M**2)
C****
C**31  68  ENERGY DIFFUSION INTO THE THERMOCLINE (W/M**2)
C****  18  MEAN TEMPERATURE OF FIRST GROUND LAYER (.1 K-TF)
C****  17  MEAN TEMPERATURE OF SECOND GROUND LAYER (.1 K-TF)
C****  34  OCEAN TEMPERATURE AT THE MAXIMUM MIXED LAYER DEPTH
C****  23  SURFACE AIR TEMPERATURE (.1 K-TF)
C****  22  FIRST LAYER AIR TEMPERATURE (.1 K-TF)
C****  21  COMPOSITE AIR TEMPERATURE (.1 K-TF)
C****  35  STRATO TEMPERATURE CHANGE PER DEGREE LATITUDE (10**-2 K)
C****  36  TROPO TEMPERATURE CHANGE PER DEGREE LATITUDE (10**-2 K)
C****  24  STRATOSPHERIC STATIC STABILITY (10**-3 K/M)
C****
C**41  25  TROPOSPHERIC STATIC STABILITY (10**-3 K/M)
C****  26  STRATOSPHERIC RICHARDSON NUMBER (1)
C****  27  TROPOSPHERIC RICHARDSON NUMBER (1)
C****  28  STRATOSPHERIC ROSSBY NUMBER (1)
C****  29  TROPOSPHERIC ROSSBY NUMBER (1)
C****  37  L IN THE STRATOSPHERE (10**5 M)
C****  38  L IN THE TROPOSPHERE (10**5 M)
C****  64  GAM  (10**-3 K/M)
C****  65  GAMM  (10**-3 K/M)
C****  66  GAMC  (10**-3 K/M)
C****
C**51  57  INTEGRATED SUPER-SATURATION CLOUD COVER (10**-2)
C****  58  INTEGRATED MOIST CONVECTIVE CLOUD COVER (10**-2)
C****  59  INTEGRATED TOTAL CLOUD COVER (10**-2)
C****  60  MOIST CONVECTIVE CLOUD DEPTH (100 N)
C****  61  SUPER SATURATION PRECIPITATION (KG/M**2/86400 S)
C****  62  MOIST CONVECTIVE PRECIPITATION (KG/M**2/86400 S)
C****  20  PRECIPITATION (KG/M**2/86400 S)
C****  19  EVAPORATION (KG/M**2/86400 S)
C****  63  WATER CONTENT OF ATMOSPHERE (KG/M**2)
C****  54  WATER RUNOFF AT Z0 (KG/M**2/86400 S)
C****
C**61  45  WATER OR ICE DIFFUSION AT -Z1 (KG/M**2/86400 S)
C****  46  ICE MELTING (OR TRANSPORTING) AT -Z1-Z2 (KG/M**2/86400 S)
C****  47  WATER RUNOFF THROUGH MIXED LAYER DEPTH (KG/M**2/86400 S)
C****  49  WATER CONTAINED IN FIRST GROUND LAYER (KG/M**2)
C****  50  ICE CONTAINED IN FIRST GROUND LAYER (KG/M**2)
C****  51  WATER CONTAINED IN SECOND GROUND LAYER (KG/M**2)
C****  52  ICE CONTAINED IN SECOND GROUND LAYER (KG/M**2)
C****  53  SNOW DEPTH (KG/M**2)
C****  31  SNOW COVER (10**-2)
C****  30  OCEAN ICE COVER (10**-2)
C**71  69  SURFACE TYPE FRACTION (1)
C****
c      USE PRTCOM, only :
      USE CONSTANT, only :
     &     grav,rgas,sday,twopi,omega,kapa,bygrav
      USE MODEL_COM, only :
     &     im,jm,lm,fim,flice,
     &     dtsrc,fland,idacc,jhour,jhour0,jdate,jdate0,amon,amon0,
     &     jyear,jyear0,ls1,sige,itime,itime0,nday,xlabel,lrunid,ntype
      USE GEOM, only :
     &     dxyp,lat,LAT_DG
      USE DAGCOM, only :
     &     QCHECK,acc_period,aj,areg,jreg,kdiag,namreg,nreg,kaj,ia_j,
     &     j_ctopp,j_cdldep,j_pcldmc,j_srabs,j_srnfp0,j_srnfg,j_trnfp0,
     &     j_hsurf,j_trhdt,j_trnfp1,j_hatm,j_rnfp0,j_rnfp1,j_srnfp1,
     &     j_rhdt,j_hz1,j_edifs,j_f1dt,j_prcp,j_prcpss,j_prcpmc,j_hz0,
     &     j_shdt,j_evhdt,j_eprcp,j_erun1,j_hz2,j_f2dt,j_erun2,j_type
      USE BDJ, only :
     &     nt_j,nstype_out,iotype
      IMPLICIT NONE
      SAVE
      DOUBLE PRECISION, DIMENSION(JM) ::
     &     CONTJ,CONTO,CONTL,CONTOI,S1,FLAT
      DOUBLE PRECISION, DIMENSION(NTYPE,JM) :: SPTYPE
      INTEGER, DIMENSION(JM) :: MLAT
      DOUBLE PRECISION, DIMENSION(NREG) :: SAREA
      INTEGER, DIMENSION(2) :: MHEM
      DOUBLE PRECISION, DIMENSION(2) :: FHEM
      INTEGER III
!@param NTYPE_OUT number of output budgets pages
      INTEGER, PARAMETER :: NTYPE_OUT=NTYPE+2  ! to include composites
C****
!@var TERRAIN name of surface type
c      CHARACTER*16, DIMENSION(0:NTYPE_OUT) :: TERRAIN = (/
c     *     '    (GLOBAL)','(OPEN OCEAN)',' (OCEAN ICE)','      (LAND)',
c     *     '  (LAND ICE)',' (OPEN LAKE)','  (LAKE ICE)'/)
C**** old version
c      CHARACTER*16 :: TERRAIN(0:NTYPE_OUT) = (/'    (GLOBAL)',
c     *     '      (LAND)','     (OCEAN)',' (OCEAN ICE)'/)
C**** possible expanded version (including composites)
      CHARACTER*16, DIMENSION(0:NTYPE_OUT) :: TERRAIN = (/
     *     '    (GLOBAL)','(OPEN OCEAN)',' (OCEAN ICE)','     (OCEAN)',
     *     '      (LAND)','  (LAND ICE)',' (OPEN LAKE)','  (LAKE ICE)',
     *     '     (LAKES)'/)
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
     8  ' ROSS NUM(STRAT)', ' ROSS NUM(TROPO)', ' OC/LK ICE COVER',
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
     1  '0SS PRECIP(MM/D)', ' MC PRECIP(MM/D)', ' H2O OF ATM (MM)',
     4  '0GAM(K/KM)      ', ' GAMM(K/KM)     ', ' GAMC(K/KM)     ',
     7  ' LW INC ON Z0   ', ' HT INTO THRMOCL', ' SURF TYPE FRACT',
     & ('                ',iii=1,11),
     *  ' NET CLR RAD  P0', ' NET CLR RAD TRP', ' NET RAD (TROPP)',
     & ('                ',iii=1,11)/)
      DATA TITLEA/' PLANETARY ALBDO',' PLAN ALB VISUAL',
     *  ' PLAN ALB NEARIR', ' SURFACE G ALBDO', ' SURF ALB VISUAL',
     *  ' SURF ALB NEARIR', '0ATMO ALB VISUAL', ' ATMO ALB NEARIR',
     *  ' ATMO ABS VISUAL', ' ATMO ABS NEARIR'/
C**** Arrays needed for full output
      REAL*8, DIMENSION(JM+3,KAJ) :: BUDG
      CHARACTER*16,DIMENSION(KAJ) :: TITLEO
C**** weighting functions for surface types
      DOUBLE PRECISION, DIMENSION(0:NTYPE_OUT,NTYPE) ::
c     *     WT=RESHAPE(         ! old version with NTYPE=3
c     *        (/1.,1.,0.,0., 1.,0.,1.,0.,  1.,0.,0.,1./),
c     *     (/NTYPE_OUT+1,NTYPE/) )
c     *     WT=RESHAPE(         ! seperate types
c     *        (/1.,1.,0.,0.,0.,0.,0., 1.,0.,1.,0.,0.,0.,0.,
c     *          1.,0.,0.,1.,0.,0.,0., 1.,0.,0.,0.,1.,0.,0.,
c     *          1.,0.,0.,0.,0.,1.,0., 1.,0.,0.,0.,0.,0.,1./),
c     *     (/NTYPE_OUT+1,NTYPE/) )
     *     WT=RESHAPE(          ! seperate types + composites
     *     (/1.,1.,0.,1.,0.,0.,0.,0.,0., 1.,0.,1.,1.,0.,0.,0.,0.,0.,
     *       1.,0.,0.,0.,1.,0.,0.,0.,0., 1.,0.,0.,0.,0.,1.,0.,0.,0.,
     *       1.,0.,0.,0.,0.,0.,1.,0.,1., 1.,0.,0.,0.,0.,0.,0.,1.,1./),
     *     (/NTYPE_OUT+1,NTYPE/) )
      INTEGER, DIMENSION(71) :: NDEX = (/
     *  1,2,3,4,5,6,7,8,9,67,           55,10,11,12,81,82,83,13,14,39,
     *     40,44,
     *  42,41,16,15,43,56,33,48,68,18,  17,34,23,22,21,35,36,24,25,26,
     *  27,28,29,37,38,64,65,66,57,58,  59,60,61,62,20,19,63,54,45,46,
     *  47,49,50,51,52,53,31,30,69/)
      INTEGER, DIMENSION(10) :: INNUM,INDEN
      DATA INNUM/2,72,73,6,74,75,76,77,78,79/, INDEN/3*1,5,6*1/
      DOUBLE PRECISION, DIMENSION(KAJ) :: SCALE
      DATA SCALE/6*1.,  6*1.,  4*1.,2*10.,  2*1.,4*10.,  6*100.,
     *  100.,2*1.,10.,2*100.,  6*1.,  6*1.,  6*1.,  2*1.,3*100.,1.,
     *  6*1., 2*1.,100.,3*1., 22*1./

      DOUBLE PRECISION :: A1BYA2,A2BYA1,AMULT,BYA1,BYIACC,
     &     FGLOB,GSUM,GSUM2,GWT,HSUM,HSUM2,HWT,QDEN,QJ,QNUM,DAYS,WTX

      INTEGER :: I,IACC,INC,J,JH,JHEMI,JMHALF,JR,K,KA,KD1M,M,MD,N,ND,NN
     *     ,IT
      INTEGER :: IFIRST = 1
      DOUBLE PRECISION, PARAMETER :: P1000=1000.
      IF (IFIRST.EQ.1) THEN
      IFIRST=0
C**** INITIALIZE CERTAIN QUANTITIES  (KD1M LE 69)
      nstype_out = ntype_out+1 ! must do this before calling open_j
      call j_titles
      KD1M=71  ! was 68, now includes J_TYPE,J_CLRTOA,J_CLRTPR,J_TOTTRP
               ! removed J_SWCOR
      INC=1+(JM-1)/24
      JMHALF=JM/2
      SAREA=0.
      DO J=1,JM
        S1(J)=IM
        DO I=1,IM
          JR=JREG(I,J)
          SAREA(JR)=SAREA(JR)+DXYP(J)
        END DO
      END DO
      S1(1)=1.
      S1(JM)=1.
      SCALE(9)=1./DTSRC
      SCALE(12)=1./DTSRC
      SCALE(13)=1./DTSRC
      SCALE(14)=1./DTSRC
      SCALE(15)=1./DTSRC
      SCALE(16)=1./DTSRC
      SCALE(19)=SDAY/DTSRC
      SCALE(20)=100.*SDAY/(DTsrc*GRAV)
      SCALE(24)=1.D3*GRAV*(P1000**KAPA)
      SCALE(25)=SCALE(24)
      SCALE(26)=16.*RGAS
      SCALE(27)=16.*RGAS
      SCALE(28)=.5/(2.*OMEGA*FIM)
      SCALE(29)=.5/(2.*OMEGA*FIM)
      SCALE(33)=1./DTSRC
      SCALE(35)=.5D2*(JM-1.)/((SIGE(LS1)-SIGE(LM+1)+1.D-12)*180.)
      SCALE(36)=.5E2*(JM-1.)/((SIGE(1)-SIGE(LS1))*180.)
      SCALE(37)=1.D-5*SQRT(RGAS)/(2.*OMEGA)
      SCALE(38)=SCALE(37)
      SCALE(39)=1./DTSRC
      SCALE(40)=1./DTSRC
      SCALE(41)=1./DTSRC
      SCALE(42)=1./DTSRC
      SCALE(43)=1./DTSRC
      SCALE(44)=1./DTSRC
      SCALE(45)=SDAY/DTSRC
      SCALE(46)=SDAY/DTSRC
      SCALE(47)=SDAY/DTSRC
      SCALE(48)=1./DTSRC
      SCALE(54)=SDAY/DTSRC
      SCALE(56)=1./DTSRC
      SCALE(61)=SCALE(20)
      SCALE(62)=SCALE(20)
      SCALE(63)=100.*BYGRAV
      SCALE(64)=1.D3*GRAV
      SCALE(65)=1.D3*.0098/(SIGE(1)-SIGE(LS1))
      SCALE(66)=1.D3
      SCALE(68)=2.E3*4185./SDAY
      END IF
C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      iotype = 0
      IF (QCHECK) call open_j(trim(acc_period)//'.j'//XLABEL(1:LRUNID))

C**** CALCULATE THE DERIVED QUANTITIES
      BYA1=1./(IDACC(1)+1.D-20)
      A2BYA1=DFLOAT(IDACC(2))/DFLOAT(IDACC(1))
      A1BYA2=IDACC(1)/(IDACC(2)+1.D-20)
      DO JR=1,23 ! only 23 will fit on a green sheet
        AREG(JR,J_TRNFP0)=AREG(JR,J_HSURF)+A2BYA1*AREG(JR,J_TRHDT)/DTSRC
        AREG(JR,J_TRNFP1)=AREG(JR,J_HATM)+A2BYA1*AREG(JR,J_TRHDT)/DTSRC
        AREG(JR,J_SRABS) =AREG(JR,J_SRNFP0)-AREG(JR,J_SRNFG)
        AREG(JR,J_RNFP0) =AREG(JR,J_SRNFP0)+AREG(JR,J_TRNFP0)
        AREG(JR,J_RNFP1) =AREG(JR,J_SRNFP1)+AREG(JR,J_TRNFP1)
        AREG(JR,J_RHDT)  =A1BYA2*AREG(JR,J_SRNFG)*DTSRC+AREG(JR,J_TRHDT)
        AREG(JR,J_PRCP)  =AREG(JR,J_PRCPSS)+AREG(JR,J_PRCPMC)
        AREG(JR,J_HZ1)   =AREG(JR,J_EDIFS)+AREG(JR,J_F1DT)
        AREG(JR,J_HZ0)   =AREG(JR,J_RHDT)+AREG(JR,J_SHDT)+
     *       AREG(JR,J_EVHDT)+AREG(JR,J_EPRCP)-AREG(JR,J_ERUN1)
        AREG(JR,J_HZ2)   =AREG(JR,J_F2DT)+AREG(JR,J_ERUN2)
        AREG(JR,J_CTOPP) =IDACC(2)*SAREA(JR)*AREG(JR,J_CDLDEP)/
     *       (AREG(JR,J_PCLDMC)+1.D-20)
      END DO
      DO J=1,JM
      DO IT=1,NTYPE
        SPTYPE(IT,J) =AJ(J,J_TYPE,IT)*BYA1
        AJ(J,J_TRNFP0,IT)=AJ(J,J_HSURF,IT)+A2BYA1*AJ(J,J_TRHDT,IT)/DTSRC
        AJ(J,J_TRNFP1,IT)=AJ(J,J_HATM,IT) +A2BYA1*AJ(J,J_TRHDT,IT)/DTSRC
c       AJ(J,J_SWCOR,IT) =(1.-SRCOR)*AJ(J,J_SRNFG,IT)   !obsolete
        AJ(J,J_SRABS ,IT)=AJ(J,J_SRNFP0,IT)-AJ(J,J_SRNFG,IT)
        AJ(J,J_RNFP0 ,IT)=AJ(J,J_SRNFP0,IT)+AJ(J,J_TRNFP0,IT)
        AJ(J,J_RNFP1 ,IT)=AJ(J,J_SRNFP1,IT)+AJ(J,J_TRNFP1,IT)
        AJ(J,J_RHDT  ,IT)=A1BYA2*AJ(J,J_SRNFG,IT)*DTSRC+AJ(J,J_TRHDT,IT)
        AJ(J,J_PRCP  ,IT)=AJ(J,J_PRCPSS,IT)+AJ(J,J_PRCPMC,IT)
        AJ(J,J_HZ1   ,IT)=AJ(J,J_EDIFS,IT)+AJ(J,J_F1DT,IT)
        AJ(J,J_HZ0   ,IT)=AJ(J,J_RHDT,IT)+AJ(J,J_SHDT,IT)+
     *       AJ(J,J_EVHDT,IT)+AJ(J,J_EPRCP,IT)-AJ(J,J_ERUN1,IT)
        AJ(J,J_HZ2,   IT)=AJ(J,J_F2DT,IT)+AJ(J,J_ERUN2,IT)
C**** this is incorrect: ratios should not be summed over type
        AJ(J,J_CTOPP,IT) =IDACC(2)*SPTYPE(IT,J)*AJ(J,J_CDLDEP,IT)/
     *       (AJ(J,J_PCLDMC,IT)+1.D-20)
      END DO
      END DO
      DAYS=(Itime-Itime0)/DFLOAT(nday)
C****
C**** LOOP OVER SURFACE TYPES: 1 TO NTYPE
C****
      IF (KDIAG(1).GT.7) GO TO 510
      DO M=0,NTYPE_OUT
      WRITE (6,901) XLABEL
      WRITE (6,902) TERRAIN(M),JYEAR0,AMON0,JDATE0,JHOUR0,
     *  JYEAR,AMON,JDATE,JHOUR,ITIME,DAYS

      WRITE (6,903) (NINT(LAT_DG(J,1)),J=JM,INC,-INC)
      WRITE (6,905)
      DO K=1,KD1M
      N=NDEX(K)
      IACC=IDACC(IA_J(N))
C**** set weighting for denominator (different only for J_TYPE)
      MD=M
      IF (N.EQ.J_TYPE) MD=0
      GSUM=0.
      GWT=0.
      DO JHEMI=1,2
        HSUM=0.
        HWT=0.
        DO JH=1,JMHALF
          J=(JHEMI-1)*JMHALF+JH
C**** Sum over types
          QJ=0
          WTX=0
          DO IT=1,NTYPE
            QJ =QJ +WT(M ,IT)*AJ(J,N,IT)
            WTX=WTX+WT(MD,IT)*SPTYPE(IT,J)
          END DO
          QJ=QJ*SCALE(N)
          WTX=WTX*IACC
          FLAT(J)=QJ/(WTX+1.D-20)
          MLAT(J)=NINT(FLAT(J))
          HSUM=HSUM+QJ*DXYP(J)*(FIM+1.-S1(J))
          HWT=HWT+WTX*DXYP(J)*(FIM+1.-S1(J))
        END DO
        FHEM(JHEMI)=HSUM/(HWT+1.D-20)
        GSUM=GSUM+HSUM
        GWT=GWT+HWT
      END DO
      FGLOB=GSUM/(GWT+1.D-20)
      IF (M.EQ.0) CALL KEYDJ (N,FGLOB,FHEM(2))
C**** Save BUDG for full output
      BUDG(1:JM,K)=FLAT(1:JM)
      BUDG(JM+1,K)=FHEM(1)
      BUDG(JM+2,K)=FHEM(2)
      BUDG(JM+3,K)=FGLOB
      TITLEO(K)=TITLE(N)
      NT_J(K)=N
c     GO TO (350,350,350,350,350,350,  350,350,350,350,350,350,
c    *       350,350,348,348,350,350,  345,345,350,350,350,350,
c    *       340,350,350,345,345,350,  350,340,348,350,350,350,
c    *       345,345,348,345,348,348,  348,348,345,345,345,345,
c    *       350,350,350,350,350,345,  350,348,350,350,350,350,
c    *       345,345,350,345,345,345,  350,345,350,350,350,350),N
      SELECT CASE (N)           ! output format
      CASE (25,32)
        WRITE (6,906) TITLE(N),FGLOB,FHEM(2),FHEM(1),
     *       (FLAT(J),J=JM,INC,-INC)
      CASE (19,20,28,29,37,38,40,45:48,54,61,62,64:66,68)
        WRITE (6,911) TITLE(N),FGLOB,FHEM(2),FHEM(1),
     *       (FLAT(J),J=JM,INC,-INC)
      CASE (15,16,33,39,41:44,56)
        WRITE (6,912) TITLE(N),FGLOB,FHEM(2),FHEM(1),
     *       (MLAT(J),J=JM,INC,-INC)
      CASE DEFAULT
        WRITE (6,907) TITLE(N),FGLOB,FHEM(2),FHEM(1),
     *       (MLAT(J),J=JM,INC,-INC)
      END SELECT
      IF (N.EQ.1) THEN
C**** CALCULATE AND PRINT ALBEDOS
      DO KA=1,10
        NN=INNUM(KA)
        ND=INDEN(KA)
        AMULT=1.
        IF (KA.LE.1.OR.KA.EQ.4) AMULT=-1.
        GSUM=0.
        GSUM2=0.
        DO JHEMI=1,2
          HSUM=0.
          HSUM2=0.
          DO JH=1,JMHALF
            J=(JHEMI-1)*JMHALF+JH
C**** Sum over types
            QNUM=0
            QDEN=0
            DO IT=1,NTYPE
              QNUM=QNUM+WT(M,IT)*AJ(J,NN,IT)
              QDEN=QDEN+WT(M,IT)*AJ(J,ND,IT)
            END DO
            FLAT(J)=AMULT*(100.*     QNUM/(QDEN    +1.D-20)-50.)+50.
            MLAT(J)=FLAT(J)+.5
            HSUM=HSUM+QNUM*DXYP(J)*(FIM+1.-S1(J))
            HSUM2=HSUM2+QDEN*DXYP(J)*(FIM+1.-S1(J))
          END DO
          FHEM(JHEMI)=50.+AMULT*(100.*HSUM/(HSUM2+1.D-20)-50.)
          GSUM=GSUM+HSUM
          GSUM2=GSUM2+HSUM2
        END DO
        FGLOB=50.+AMULT*(100.*GSUM/(GSUM2+1.D-20)-50.)
        IF (M.EQ.0.AND.KA.EQ.1) CALL KEYDJA (FGLOB)
C**** Save BUDG for full output
      BUDG(1:JM,KA+KD1M)=FLAT(1:JM)
      BUDG(JM+1,KA+KD1M)=FHEM(1)
      BUDG(JM+2,KA+KD1M)=FHEM(2)
      BUDG(JM+3,KA+KD1M)=FGLOB
      TITLEO(KA+KD1M)=TITLE(KA)
      NT_J(KA+KD1M)=KA+KD1M
      WRITE (6,912) TITLEA(KA),FGLOB,FHEM(2),FHEM(1),
     *     (MLAT(J),J=JM,INC,-INC)
      END DO
      END IF
      END DO
      WRITE (6,903) (NINT(LAT_DG(J,1)),J=JM,INC,-INC)
      WRITE (6,905)
      iotype = iotype + 1
      IF (QCHECK) CALL POUT_J(TITLEO,BUDG,KD1M+10,TERRAIN(M))
      IF (KDIAG(1).GT.1) RETURN
      END DO
      if(qcheck) call close_j
  510 IF (KDIAG(1).GT.0.AND.KDIAG(1).NE.8) RETURN
C****
C**** PRODUCE REGIONAL STATISTICS
C****
      WRITE (6,901) XLABEL
      WRITE (6,902) '   (REGIONS)    ',JYEAR0,AMON0,JDATE0,JHOUR0,
     *  JYEAR,AMON,JDATE,JHOUR,ITIME,DAYS
      WRITE(6,918)(NAMREG(1,K),K=1,23),(NAMREG(2,K),K=1,23)
      DO K=1,KD1M
      N=NDEX(K)
      BYIACC=1./(IDACC(IA_J(N))+1.D-20)
      DO JR=1,23
        FLAT(JR)=AREG(JR,N)*SCALE(N)*BYIACC/SAREA(JR)
        MLAT(JR)=NINT(FLAT(JR))
      END DO
CF       DO 523 J=1,23
CF523    RBUDG(J,K)=FLAT(J)
c     GO TO (550,550,550,550,550,550,  550,550,550,550,550,550,
c    *       550,550,550,550,550,550,  540,540,550,550,550,550,
c    *       540,550,550,540,540,540,  540,540,550,550,550,550,
c    *       540,540,550,540,550,550,  550,550,540,540,540,540,
c    *       540,540,550,550,540,540,  550,550,550,550,550,550,
c    *       540,540,540,540,540,540,  550,540,550,550,550,550),N
      SELECT CASE (N) ! output format
      CASE (19,20,25,28:32,37,38,40,45:50,53,54,61:66,68)
        WRITE (6,910) TITLE(N),(FLAT(JR),JR=1,23)
      CASE DEFAULT
        WRITE (6,909) TITLE(N),(MLAT(JR),JR=1,23)
      END SELECT
      IF (N.EQ.1) THEN ! do albedo printout next
C**** CALCULATE AND PRINT ALBEDOS FOR REGIONAL STATISTICS
      DO KA=1,10
        NN=INNUM(KA)
        ND=INDEN(KA)
        AMULT=1.
        IF (KA.LE.1.OR.KA.EQ.4) AMULT=-1.
        DO JR=1,23
          FLAT(JR)=AMULT*(100.*AREG(JR,NN)/(AREG(JR,ND)+1.D-20)-50.)+50.
          MLAT(JR)=FLAT(JR)+.5
        END DO
CF       DO 613 J=1,23
CF613    RBUDG(J,KA+KD1M)=FLAT(J)
        WRITE (6,909) TITLEA(KA),(MLAT(JR),JR=1,23)
      END DO
      END IF
      END DO
      WRITE (6,905)
      WRITE(6,918)(NAMREG(1,K),K=1,23),(NAMREG(2,K),K=1,23)
      RETURN
C****
  901 FORMAT ('1',A)
  902 FORMAT ('0** BUDGETS',A16,'**   From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,'  Model-Time:',I9,5X,
     *  'Dif:',F7.2,' Days')
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

      MODULE BDJK
!@sum  stores information for outputting lat-pressure diagnostics
!@auth M. Kelley

      IMPLICIT NONE

!@param njk_out number of jk-format output fields (not all used)
      integer, parameter :: njk_out=44

!@var title string, formed as concatentation of lname//units
      CHARACTER(LEN=64), DIMENSION(njk_out) :: TITLE
!@var units string containing output field units
      CHARACTER(LEN=50), DIMENSION(njk_out) :: UNITS
!@var lname string describing output field
      CHARACTER(LEN=50), DIMENSION(njk_out) :: LNAME
!@var sname string referencing output field in self-desc. output file
      CHARACTER(LEN=30), DIMENSION(njk_out) :: SNAME

!@var in_jkmap flag telling pout_jl to use jk titles
      logical :: in_jkmap
!@var nt_jk index telling pout_jl which field is being output
      integer :: nt_jk

      END MODULE BDJK

      SUBROUTINE JK_TITLES
      USE BDJK
      IMPLICIT NONE
      INTEGER K
c
      k = 0
c
      k = k + 1 ! 01
      sname(k) = 'temp'
      lname(k) = 'TEMPERATURE'
      units(k) = 'DEGREES CENTIGRADE'
      k = k + 1 ! 02
      sname(k) = 'height'
      lname(k) = 'HEIGHT'
      units(k) = 'HUNDREDS OF METERS'
      k = k + 1 ! 03
      sname(k) = 'q'
      lname(k) = 'SPECIFIC HUMIDITY'
      units(k) = '10**-5 KG H2O/KG AIR'
      k = k + 1 ! 04
      sname(k) = 'rh'
      lname(k) = 'RELATIVE HUMIDITY'
      units(k) = 'PERCENT'
      k = k + 1 ! 05
      sname(k) = 'u'
      lname(k) = 'ZONAL WIND (U COMPONENT)'
      units(k) = 'TENTHS OF METERS/SECOND'
      k = k + 1 ! 06
      sname(k) = 'v'
      lname(k) = 'MERIDIONAL WIND (V COMPONENT)'
      units(k) = 'HUNDREDTHS OF METERS/SECOND'
      k = k + 1 ! 07
      sname(k) = ''
      lname(k) = ''
      units(k) = ''
      k = k + 1 ! 08
      sname(k) = ''
      lname(k) = ''
      units(k) = ''
      k = k + 1 ! 09
      sname(k) = 'baroc_eddy_ke_gen'
      lname(k) = 'BAROCLINIC EDDY KINETIC ENERGY GEN.'
      units(k) = '10**-1 WATTS/M**2/SIGMA'
      k = k + 1 ! 10
      sname(k) = 'nt_eddy_qgpv'
      lname(k) = 'NORTH. TRANS. OF EDDY Q-G POT. VORTICITY'
      units(k) = '10**-6 M/S**2'
      k = k + 1 ! 11
      sname(k) = 'p2k_eddy_pgf'
      lname(k) = 'P-K BY EDDY PRESSURE GRADIENT FORCE'
      units(k) = '10**-1 W/M**2/UNIT SIGMA'
      k = k + 1 ! 12
      sname(k) = 'dyn_conv_eddy_geop'
      lname(k) = 'DYNAMIC CONVERGENCE OF EDDY GEOPOTENTIAL'
      units(k) = '.1 WATTS/M**2/DSIGMA'
      k = k + 1 ! 13
      sname(k) = 'nt_sheat_eddy'
      lname(k) = 'NORTH. TRANS. OF SENSIBLE HEAT BY EDDIES'
      units(k) = '10**14 WATTS/DSIGMA'
      k = k + 1 ! 14
      sname(k) = 'dyn_conv_dse'
      lname(k) = 'DYNAMIC CONVERGENCE OF DRY STATIC ENERGY'
      units(k) = '10 WATTS/M**2/DSIGMA'
      k = k + 1 ! 15
      sname(k) = ''
      lname(k) = ''
      units(k) = ''
      k = k + 1 ! 16
      sname(k) = ''
      lname(k) = ''
      units(k) = ''
      k = k + 1 ! 17
      sname(k) = 'stand_eddy_ke'
      lname(k) = 'STANDING EDDY KINETIC ENERGY'
      units(k) = '10**4 JOULES/M**2/UNIT SIGMA'
      k = k + 1 ! 18
      sname(k) = 'eddy_ke'
      lname(k) = 'EDDY KINETIC ENERGY'
      units(k) = '10**4 JOULES/M**2/UNIT SIGMA'
      k = k + 1 ! 19
      sname(k) = 'tot_ke'
      lname(k) = 'TOTAL KINETIC ENERGY'
      units(k) = '10**4 JOULES/M**2/UNIT SIGMA'
      k = k + 1 ! 20
      sname(k) = ''
      lname(k) = ''
      units(k) = ''
      k = k + 1 ! 21
      sname(k) = 'pot_temp'
      lname(k) = 'POTENTIAL TEMPERATURE'
      units(k) = 'DEGREES KELVIN'
      k = k + 1 ! 22
      sname(k) = 'nt_dse_stand_eddy'
      lname(k) = 'NOR. TRANS. OF DRY STAT. ENERGY BY STAND. EDDIES'
      units(k) = '10**14 W/DSIG'
      k = k + 1 ! 23
      sname(k) = 'nt_dse_eddy'
      lname(k) = 'NORTH. TRANS. OF DRY STATIC ENERGY BY EDDIES'
      units(k) = '10**14 WATTS/DSIG'
      k = k + 1 ! 24
      sname(k) = 'tot_nt_dse'
      lname(k) = 'TOTAL NORTH. TRANSPORT OF DRY STATIC ENERGY'
      units(k) = '10**15 WATTS/DSIG'
      k = k + 1 ! 25
      sname(k) = 'nt_lh_eddy'
      lname(k) = 'NORTHWARD TRANSPORT OF LATENT HEAT BY EDDIES'
      units(k) = '10**13 WATTS/DSIG'
      k = k + 1 ! 26
      sname(k) = 'tot_nt_lh'
      lname(k) = 'TOTAL NORTHWARD TRANSPORT OF LATENT HEAT'
      units(k) = '10**14 WATTS/UNIT SIG'
      k = k + 1 ! 27
      sname(k) = 'nt_se_eddy'
      lname(k) = 'NORTH.TRANSPORT OF STATIC ENERGY BY EDDIES'
      units(k) = '10**14 WATTS/DSIGMA'
      k = k + 1 ! 28
      sname(k) = 'tot_nt_se'
      lname(k) = 'TOTAL NORTHWARD TRANSPORT OF STATIC ENERGY'
      units(k) = '10**15 WATTS/DSIGMA'
      k = k + 1 ! 29
      sname(k) = ''
      lname(k) = ''
      units(k) = ''
      k = k + 1 ! 30
      sname(k) = 'tot_nt_ke'
      lname(k) = 'TOTAL NORTHWARD TRANSPORT OF KINETIC ENERGY'
      units(k) = '10**12 WATTS/DSIG'
      k = k + 1 ! 31
      sname(k) = 'nt_am_stand_eddy'
      lname(k) = 'NORTH. TRANS. OF ANG. MOMENTUM BY STAND. EDDIES'
      units(k) = '10**18 J/DSIG'
      k = k + 1 ! 32
      sname(k) = 'nt_am_eddy'
      lname(k) = 'NORTH. TRANS. OF ANG. MOMENTUM BY EDDIES'
      units(k) = '10**18 JOULES/DSIGMA'
      k = k + 1 ! 33
      sname(k) = 'tot_nt_am'
      lname(k) = 'TOTAL NORTHWARD TRANSPORT OF ANG. MOMENTUM'
      units(k) = '10**19 JOULES/DSIG'
      k = k + 1 ! 34
      sname(k) = 'we_flx_nor'
      lname(k) = 'NORTHWARD WAVE ENERGY FLUX'
      units(k) = '10**11 JOULES/METER/UNIT SIGMA'
      k = k + 1 ! 35
      sname(k) = '' !'we_flx_vert'
      lname(k) = '' !'VERTICAL WAVE ENERGY FLUX'
      units(k) = '' !'10**11 JOULES/METER'
      k = k + 1 ! 36
      sname(k) = 'we_flx_div'
      lname(k) = 'DIVERGENCE OF THE WAVE ENERGY FLUX'
      units(k) = '10**-5 M/S**2'
      k = k + 1 ! 37
      sname(k) = 'refr_ind_wave1'
      lname(k) = 'REFRACTION INDEX FOR WAVE NUMBER 1'
      units(k) = '10**-8 PER METER**2'
      k = k + 1 ! 38
      sname(k) = 'refr_ind_wave2'
      lname(k) = 'REFRACTION INDEX FOR WAVE NUMBER 2'
      units(k) = '10**-8 PER METER**2'
      k = k + 1 ! 39
      sname(k) = 'refr_ind_wave3'
      lname(k) = 'REFRACTION INDEX FOR WAVE NUMBER 3'
      units(k) = '10**-8 PER METER**2'
      k = k + 1 ! 40
      sname(k) = 'refr_ind_wave6'
      lname(k) = 'REFRACTION INDEX FOR WAVE NUMBER 6'
      units(k) = '10**-8 PER METER**2'
      k = k + 1 ! 41
      sname(k) = 'refr_ind_wave9'
      lname(k) = 'REFRACTION INDEX FOR WAVE NUMBER 9'
      units(k) = '10**-8 PER METER**2'
      k = k + 1 ! 42
      sname(k) = 'del_qgpv'
      lname(k) = 'Q-G POT. VORTICITY CHANGE OVER LATITUDES'
      units(k) = '10**-12 1/(SEC-M)'
      k = k + 1 ! 43
      sname(k) = 'cldh2o'
      lname(k) = 'TOTAL CLOUD WATER CONTENT'
      units(k) = '10**-6 KG/KG'
      k = k + 1 ! 44
      sname(k) = 'nt_lh_stand_eddy'
      lname(k) = 'N. TRANSPORT OF LATENT HEAT BY STAND. EDDIES'
      units(k) = '10**13 WATTS/DSIG'

c create titles by concatenating long names with units
c no checks whether total length of lname+units exceeds length of title
      do k=1,njk_out
         title(k)=''
         if(lname(k).ne.'')
     &        title(k) = trim(lname(k))//' ('//trim(units(k))//')'
      enddo
      RETURN
      END SUBROUTINE JK_TITLES

      SUBROUTINE DIAGJK
c      USE PRTCOM, only :
      USE CONSTANT, only :
     &     grav,rgas,kapa,sday,lhe,twopi,omega,sha,bygrav,tf
      USE MODEL_COM, only :
     &     im,jm,lm,fim, xlabel,lrunid,
     &     BYIM,DSIG,DT,DTsrc,IDACC,IMH,LS1,NDAA,nidyn,
     &     PTOP,PMTOP,PSFMPT,SIG,SIGE,JHOUR    ! ,skipse
      USE GEOM, only :
     &     AREAG,BYDXYP,COSP,COSV,DLON,DXV,DXYP,DXYV,DYP,FCOR,RADIUS,WTJ
      USE DAGPCOM, only :
     &     PLM
      USE DAGCOM, only :
     &     ajk,ajl,asjl,ajlsp,kdiag,aijl,aijk,nwav_dag,kajlsp,LM_REQ
     &     ,qcheck, acc_period,ijk_u,ijk_v,ijk_t,ijk_q,ijk_dp,ijk_dse
     *     ,kep
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

      INTEGER :: LINECT,JMHALF,INC
      COMMON/DJLCOM/LINECT,JMHALF,INC

      DOUBLE PRECISION, DIMENSION(JM,LM,0:NWAV_DAG,KAJLSP):: SPSTAD
      DOUBLE PRECISION, DIMENSION(JM,LM,0:NWAV_DAG,KAJLSP):: SPTRAN
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

C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      IF(QCHECK) call open_jl(trim(acc_period)//'.jk'//XLABEL(1:LRUNID))

C**** INITIALIZE CERTAIN QUANTITIES
      call jk_titles
      call jl_titles

      IF (KDIAG(2).GE.8) GO TO 120
      XWON=TWOPI/(DLON*FIM)
      LMP1=LM+1
      KM=LM
      KMM1=KM-1
      INC=1+(JM-1)/24
      JMHALF=JM/2
      BY100G=.01*BYGRAV
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
      SCALE=100.D-9*XWON*BYIADA*BYGRAV
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
c     IF (SKIPSE.EQ.1.) GO TO 180
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
C     DEK(I,J)=DEK(I,J)+AIJK(I,J,K,IJK_DSE)    
C     ELK(I,J)=ELK(I,J)+AIJK(I,J,K,IJK_Q)
C 125 AMK(I,J)=AMK(I,J)+AIJK(I,J,K,IJK_U)
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
      DO 156 N=0,NWAV_DAG
  156 SPSTAD(J,K,N,1)=.5*FIM*(FCPVA(N,J,K)*FCJKA(N,J,K)+
     *   FCPVB(N,J,K)*FCJKB(N,J,K))/DSSIG
      CALL FFT(AIJL(1,J,K,4),FCJKA(0,J,K),FCJKB(0,J,K))
      DO 157 N=0,NWAV_DAG
  157 SPSTAD(J,K,N,2)=.5*FIM*(FCPVA(N,J,K)*FCJKA(N,J,K)+
     *   FCPVB(N,J,K)*FCJKB(N,J,K))/DSSIG
      CALL FFT(AIJL(1,J,K,1),FCJKA(0,J,K),FCJKB(0,J,K))
      DO 158 N=0,NWAV_DAG
  158 SPSTAD(J,K,N,3)=.5*FIM*(FCPVA(N,J,K)*FCJKA(N,J,K)+
     *   FCPVB(N,J,K)*FCJKB(N,J,K))/DSSIG
      DO 159 N=0,NWAV_DAG
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
      IF (AIJK(I,J,K,IJK_DP).EQ.0.) GO TO 160
      DPTI=DPTI+AIJK(I,J,K,IJK_DP)
      BYDPK=1./(AIJK(I,J,K,IJK_DP)+1.D-20)
      PUTI=PUTI+AIJK(I,J,K,IJK_U)
      PVTI=PVTI+AIJK(I,J,K,IJK_V)
      DE4TI=DE4TI+AIJK(I,J,K,IJK_DSE)
      EL4QI=EL4QI+AIJK(I,J,K,IJK_Q)
      SKEI=SKEI+(AIJK(I,J,K,IJK_U)*AIJK(I,J,K,IJK_U)
     *            +AIJK(I,J,K,IJK_V)*AIJK(I,J,K,IJK_V))*BYDPK
      SNDEGI=SNDEGI+(AIJK(I,J,K,IJK_DSE)*AIJK(I,J,K,IJK_V)*BYDPK)
      SENDEG(I,J)=SENDEG(I,J)
     *  +(AIJK(I,J,K,IJK_DSE)*AIJK(I,J,K,IJK_V)*BYDPK)
      SNELGI=SNELGI+(AIJK(I,J,K,IJK_Q)*AIJK(I,J,K,IJK_V)*BYDPK)
      SNAMI=SNAMI+AIJK(I,J,K,IJK_U)*AIJK(I,J,K,IJK_V)*BYDPK
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
  180 SCALE=50.D-4*BYIMDA*BYGRAV
c     IF (SKIPSE.EQ.1.) GO TO 190
      CALL JKMAP (17,PLM,AX,SCALE,ONES,ONES,KM,2,2)
  190 DO 200 K=1,KM
      DO 200 J=2,JM
  200 AX(J,K)=AJK(J,K,11)-AJK(J,K,10)
      CALL JKMAP (18,PLM,AX,SCALE,ONES,ONES,KM,2,2)
      CALL JKMAP (19,PLM,AJK(1,1,11),SCALE,ONES,ONES,KM,2,2)
C**** POTENTIAL TEMPERATURE, POTENTIAL VORTICITY
      DO 205 LR=1,LM_REQ
      DO 205 J=1,JM
  205 ARQX(J,LR)=ASJL(J,LR,1)*BYIMDA*ONESPO(J)+TF
      CALL JKMAPS (21,PLM,AJK(1,1,6),P1000K,BYP,ONES,KM,2,1,
     *  ARQX,P1000K,ONES,BYPKS)
      SCALE=1.D6*BYIMDA*P1000K
      CALL JLMAP (63,PLM,AJK(1,1,32),SCALE,BYDXYP,ONES,KM,2,1)
C****
C**** NORTHWARD TRANSPORTS AT CONSTANT PRESSURE
C****
C**** NORTHWARD TRANSPORT OF SENSIBLE HEAT BY EDDIES
      SCALE=25.D-14*XWON*SHA*BYIADA*BYGRAV
      DO 210 K=1,KM
      DO 210 J=2,JM
  210 AX(J,K)=AJK(J,K,13)-AJK(J,K,12)
      CALL JKMAP (13,PLM,AX,SCALE,DXV,ONES,KM,2,2)
C**** NORTHWARD TRANSPORT OF DRY STATIC ENERGY BY STANDING EDDIES,
C****   EDDIES, AND TOTAL
C**** Individual wave transports commented out. (gas - 05/2001)
      SCALE=25.D-14*XWON*BYIADA*BYGRAV
c     IF (SKIPSE.EQ.1.) GO TO 220
      CALL JKMAP (22,PLM,BX,SCALE,DXV,ONES,KM,2,2)
C      DO 219 N=1,5
C  219 CALL JLMAP(114+N,PLM,SPSTAD(1,1,N,1),SCALE,DXV,ONES,LM,2,2)
  220 DO 230 K=1,KM
      DO 230 J=2,JM
      AX(J,K)=SHA*(AJK(J,K,13)-AJK(J,K,12))+(AJK(J,K,15)-AJK(J,K,14))
  230 BX(J,K)=SHA*AJK(J,K,13)+AJK(J,K,15)
      CALL JKMAP (23,PLM,AX,SCALE,DXV,ONES,KM,2,2)
C      DO 231 N=1,9
C  231 CALL JLMAP(119+N,PLM,SPTRAN(1,1,N,1),SCALE,DXV,ONES,LM,2,2)
      SCALE=SCALE*.1
      CALL JKMAP (24,PLM,BX,SCALE,DXV,ONES,KM,2,2)
C**** NORTHWARD TRANSPORT OF LATENT HEAT BY STAND. EDDY, EDDIES AND TOTA
      DO 240 K=1,KM
      DO 240 J=2,JM
      DX(J,K)=AJK(J,K,17)-AJK(J,K,16)
  240 AX(J,K)=AX(J,K)+LHE*DX(J,K)
      SCALE=25.D-13*XWON*LHE*BYIADA*BYGRAV
c     IF(SKIPSE.EQ.1.) GO TO 242
      CALL JKMAP(44,PLM,EX,SCALE,DXV,ONES,KM,2,2)
C      DO 241 N=1,5
C  241 CALL JLMAP(128+N,PLM,SPSTAD(1,1,N,2),SCALE,DXV,ONES,LM,2,2)
  242 CONTINUE
      CALL JKMAP (25,PLM,DX,SCALE,DXV,ONES,KM,2,2)
C      DO 243 N=1,9
C  243 CALL JLMAP(133+N,PLM,SPTRAN(1,1,N,2),SCALE,DXV,ONES,LM,2,2)
      SCALE=SCALE*.1
      CALL JKMAP (26,PLM,AJK(1,1,17),SCALE,DXV,ONES,KM,2,2)
C**** NORTHWARD TRANSPORT OF STATIC ENERGY BY EDDIES AND TOTAL
      DO 245 K=1,KM
      DO 245 J=2,JM
  245 DX(J,K)=BX(J,K)+LHE*AJK(J,K,17)
      SCALE=25.D-14*XWON*BYIADA*BYGRAV
      CALL JKMAP (27,PLM,AX,SCALE,DXV,ONES,KM,2,2)
      SCALE=SCALE*.1
      CALL JKMAP (28,PLM,DX,SCALE,DXV,ONES,KM,2,2)
C**** NORTHWARD TRANSPORT OF KINETIC ENERGY
      SCALE=50.D-12*XWON*BYIADA*BYGRAV
      CALL JKMAP (30,PLM,AJK(1,1,19),SCALE,DXV,ONES,KM,2,2)
C**** NOR. TRANS. OF ANG. MOMENTUM BY STANDING EDDIES, EDDIES AND TOTAL
      SCALE=100.D-18*XWON*RADIUS*BYIADA*BYGRAV
c     IF (SKIPSE.EQ.1.) GO TO 250
      CALL JKMAP (31,PLM,CX,SCALE,DXCOSV,ONES,KM,2,2)
C      DO 249 N=1,5
C  249 CALL JLMAP(142+N,PLM,SPSTAD(1,1,N,3),SCALE,DXCOSV,ONES,LM,2,2)
  250 DO 260 K=1,KM
      DO 260 J=2,JM
      CX(J,K)=AJK(J,K,21)-AJK(J,K,20)
  260 DX(J,K)=AJK(J,K,21)+RADIUS*OMEGA*COSV(J)*AJK(J,K,9)
      CALL JKMAP (32,PLM,CX,SCALE,DXCOSV,ONES,KM,2,2)
C      DO 261 N=1,9
C  261 CALL JLMAP(147+N,PLM,SPTRAN(1,1,N,3),SCALE,DXCOSV,ONES,LM,2,2)
      SCALE=.1*SCALE
      CALL JKMAP (33,PLM,DX,SCALE,DXCOSV,ONES,KM,2,2)
C****
C**** DYNAMIC CONVERGENCE OF ENERGY
C****
      SCALE=25.D-1*BYIMDA*BYGRAV
      DO 370 K=1,KM
C     CX(1,K)=-BX(2,K)*SCALE*DXV(2)
      CX(1,K)=0.
      CX(JM,K)=BX(JM,K)*SCALE*DXV(JM)
C     DX(1,K)=-(AJK(2,K,15)-AJK(2,K,14))*SCALE*DXV(2)
      DX(1,K)=0.
      DX(JM,K)=(AJK(JM,K,15)-AJK(JM,K,14))*SCALE*DXV(JM)
      DO 370 J=2,JM-1
      CX(J,K)=(BX(J,K)*DXV(J)-BX(J+1,K)*DXV(J+1))*SCALE
      DX(J,K)=((AJK(J,K,15)-AJK(J,K,14))*DXV(J) -
     *  (AJK(J+1,K,15)-AJK(J+1,K,14))*DXV(J+1))*SCALE
  370 CONTINUE
      SCALE=-100.D-1*BYIMDA*BYGRAV
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
      SCALE=50.E1*BYIMDA*RGAS*BYGRAV
      CALL JKMAP (9,PLM,AJK(1,1,31),SCALE,BYDXYP,ONES,KM,2,1)
      SCALE=100.E1*BYIMDA/(GRAV*2.*DT)
      CALL JKMAP (11,PLM,AJK(1,1,22),SCALE,ONES,ONES,KM,2,2)
C****
C**** VERTICAL TRANSPORTS
C****
C**** VERTICAL TRANSPORT OF GEOPOTENTIAL ENERGY BY EDDIES
      SCALE=-100.D-12*XWON*BYIADA*BYGRAV
      CALL JLMAP (99,PM,AJK(1,1,30),SCALE,ONES,ONES,KMM1,1,1)
C**** VERTICAL TRANSPORT OF DRY STATIC ENERGY BY EDDIES AND TOTAL
      DO 390 K=1,KMM1
      DO 390 J=1,JM
      AX(J,K)=AJK(J,K,27)-AJK(J,K,26)
  390 BX(J,K)=AJK(J,K,29)-AJK(J,K,28)
      SCALE=-100.D-12*XWON*BYIADA*BYGRAV
      CALL JLMAP (100,PM,AX,SCALE,ONES,ONES,KMM1,1,1)
      SCALE=SCALE*.01
      CALL JLMAP (101,PM,AJK(1,1,27),SCALE,ONES,ONES,KMM1,1,1)
C**** VERTICAL TRANSPORT OF LATENT HEAT BY EDDIES AND TOTAL
      SCALE=-100.D-12*XWON*LHE*BYIADA*BYGRAV
      CALL JLMAP (102,PM,BX,SCALE,ONES,ONES,KMM1,1,1)
      SCALE=SCALE*.1
      CALL JLMAP (103,PM,AJK(1,1,29),SCALE,ONES,ONES,KMM1,1,1)
C**** VERTICAL TRANSPORT OF STATIC ENERGY BY EDDIES AND TOTAL
      DO 420 K=1,KMM1
      DO 420 J=1,JM
      AX(J,K)=AX(J,K)+LHE*BX(J,K)
  420 BX(J,K)=AJK(J,K,27)+LHE*AJK(J,K,29)
      SCALE=-100.D-13*XWON*BYIADA*BYGRAV
      CALL JLMAP (104,PM,AX,SCALE,ONES,ONES,KMM1,1,1)
      SCALE=SCALE*.1
      CALL JLMAP (105,PM,BX,SCALE,ONES,ONES,KMM1,1,1)
C**** VERTICAL TRANSPORT OF KINETIC ENERGY
      SCALE=-12.5E-11*XWON*BYIADA*BYGRAV
      CALL JLMAP (110,PM,AJK(1,1,36),SCALE,ONES,ONES,KMM1,1,2)
C**** VERTICAL TRANSPORT OF ANGULAR MOMENTUM BY LARE SCALE MOTIONS
      SCALE=-25.D-16*XWON*RADIUS*BYIADA*BYGRAV
      CALL JLMAP (111,PM,AJK(1,1,37),SCALE,COSV,ONES,KMM1,1,2)
      SCALE=1.D-2*SCALE
      CALL JLMAP (112,PM,AJK(1,1,38),SCALE,COSV,ONES,KMM1,1,2)
C**** VERTICAL TRANSPORT OF POTENTIAL VORTICITY TOTAL AND BY EDDIES
      SCALE=-25.D-4*XWON*P1000K*BYIADA*BYGRAV
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
      SCALE=25.D-11*XWON*BYIADA*BYGRAV
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
      SCALE=-100.D-11*XWON*BYIADA*BYGRAV
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
      DALPHA=(AJK(J,KUP,3)/(AJK(J,KUP,1)+1.D-20)+TF)/PUP-
     *  (AJK(J,KDN,3)/(AJK(J,KDN,1)+1.D-20)+TF)/PDN
      DTHETA=AJK(J,KUP,6)/(AJK(J,KUP,1)+1.D-20)-
     *  AJK(J,KDN,6)/(AJK(J,KDN,1)+1.D-20)
      THETA=AJK(J,K,6)/(AJK(J,K,1)+1.D-20)
      TX=AJK(J,K,3)/(AJK(J,K,1)+1.D-20)+TF
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
      TX=BYDP2*(AJK(J-1,K,3)+AJK(J,K,3))+TF
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
      SCALE=1.D6/((IDACC(4)-1)*(DTsrc*NDAA+DT+DT))
      CALL JLMAP (40,PLM,AJK(1,1,47),SCALE,ONES,ONES,KM,2,2)
  730 CONTINUE
C**** Depending on whether EP fluxes have been specially calcualted
C**** output full or approximate version
      IF (KEP.gt.0) THEN
        CALL EPFLXP
      ELSE ! these are not very good 
        SCALE=1.D6*BYIADA
        CALL JLMAP (59,PLM,AJK(1,1,40),SCALE,ONES,ONES,KM,2,1)
        SCALE=1.D6
        CALL JLMAP (60,PLM,AX,SCALE,ONES,ONES,KM,2,1)
C**** WIND: TRANSFORMED ADVECTION, LAGRANGIAN CONVERGENCE (DEL.F)
        SCALE=1.D6*BYIADA
        CALL JLMAP (64,PLM,AJK(1,1,42),SCALE,ONES,ONES,KM,2,1)
        CALL JLMAP (65,PLM,BX,SCALE,ONES,ONES,KM,2,1)
      END IF
C**** WIND: DU/DT BY STRAT. DRAG
      SCALE=1.D6/(FIM*IDACC(1)*DTsrc+1.E-20)
      CALL JLMAP (1,PLM,AJL(1,1,20),SCALE,ONES,ONES,LM,2,2)
C**** DU/DT BY SDRAG
      SCALE=1.D6/(FIM*IDACC(1)*DTsrc+1.D-20)
      CALL JLMAP (53,PLM,AJL(1,1,52),SCALE,ONES,ONES,LM,2,2)
C**** TEMPERATURE: RATE OF CHANGE, ADVECTION, EDDY CONVERGENCE
      IF (IDACC(4).LE.1) GO TO 750
      SCALE=1.D1*SDAY/((IDACC(4)-1)*(DTsrc*NDAA+DT+DT))
      CALL JLMAP (66,PLM,AJK(1,1,49),SCALE,ONES,PKM,KM,2,1)
  750 SCALE=1.D1*SDAY*BYIADA
      CALL JLMAP (67,PLM,AJK(1,1,41),SCALE,ONES,PKM,KM,2,1)
      SCALE=1.D1*SDAY
      CALL JLMAP (69,PLM,CX,SCALE,ONES,ONES,KM,2,1)
C**** TEMPERATURE: TRANSFORMED ADVECTION
      SCALE=1.D1*SDAY*BYIADA
      CALL JLMAP (70,PLM,AJK(1,1,43),SCALE,ONES,PKM,KM,2,1)
C**** CHANGE IN TEMPERATURE BY STRATOSPHERIC DRAG
      SCALE=1.E1*SDAY/(FIM*IDACC(1)*DTsrc+1.D-20)
      CALL JLMAP (29,PLM,AJL(1,1,33),SCALE,ONES,PKM,KM,2,1)
C**** CHANGE IN TEMPERATURE BY DYNAMICS
      SCALE=1.E1*SDAY*NIDYN/(FIM*IDACC(4)*7200.+1.E-20)
      CALL JLMAP (21,PLM,AJL(1,1,17),SCALE,ONESPO,PKM,KM,2,1)
  799 CONTINUE
      RETURN
  901 FORMAT (
     *  '010**14 WATTS = .2067 * 10**19 CALORIES/DAY'/
     *  ' 10**18 JOULES = .864 * 10**30 GM*CM**2/SEC/DAY')
      END SUBROUTINE DIAGJK

      SUBROUTINE JKMAP (NT,PM,AX,SCALE,SCALEJ,SCALEK,KMAX,JWT,J1)
      USE DAGCOM, only : QCHECK,acc_period,iu_jl,lm_req
      USE MODEL_COM, only :
     &     jm,lm,JDATE,JDATE0,JMON0,JMON,AMON0,AMON,JYEAR,JYEAR0,XLABEL
      USE GEOM, only :
     &     LAT_DG,WTJ
      USE BDJK, only :
     &     title,in_jkmap,nt_jk
      IMPLICIT NONE

      INTEGER, DIMENSION(JM) :: MLAT
      DOUBLE PRECISION, DIMENSION(JM) :: FLAT,ASUM
      DOUBLE PRECISION, DIMENSION(2) :: AHEM
      DOUBLE PRECISION, DIMENSION(JM,LM) :: CX

      DOUBLE PRECISION, DIMENSION(JM,LM,2) :: DSJK
      DOUBLE PRECISION, DIMENSION(2,LM,2) :: DSHEM
      DOUBLE PRECISION, DIMENSION(LM,2) :: DSGLOB
      COMMON/WORK5/DSJK,DSHEM,DSGLOB

      INTEGER :: LINECT,JMHALF,INC
      COMMON/DJLCOM/LINECT,JMHALF,INC

      INTEGER :: J1,JWT,KMAX,NT
      DOUBLE PRECISION :: SCALE,SCALER
      DOUBLE PRECISION, DIMENSION(JM,LM) :: AX
      DOUBLE PRECISION, DIMENSION(JM,LM_REQ) :: ARQX
      DOUBLE PRECISION, DIMENSION(JM) :: SCALEJ,SCALJR
      DOUBLE PRECISION, DIMENSION(LM) :: SCALEK
      DOUBLE PRECISION, DIMENSION(LM_REQ) :: SCALLR
      DOUBLE PRECISION, DIMENSION(LM+LM_REQ) :: PM

      CHARACTER*4 DASH,WORD(4)
      DATA DASH/'----'/,WORD/'SUM','MEAN',' ','.1*'/

      INTEGER :: IWORD,J,J0,JH,JHEMI,K,L ,ksx,klmax
      DOUBLE PRECISION :: AGLOB,FGLOB,FLATJ,G1,H1,H2,SUMFAC

      REAL*8, DIMENSION(JM+3,LM+4) :: XJL ! for binary output
      CHARACTER XLB*16,CLAT*16,CPRES*16,CBLANK*16,TITLEO*80
      DATA CLAT/'LATITUDE'/,CPRES/'PRESSURE (MB)'/,CBLANK/' '/
C****
C**** PRODUCE A LATITUDE BY LAYER TABLE OF THE ARRAY A
C****
   10 LINECT=LINECT+KMAX+7
      IF (LINECT.LE.60) GO TO 20
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
      LINECT=KMAX+8
   20 WRITE (6,901) TITLE(NT),(DASH,J=J1,JM,INC)
      WRITE (6,904) WORD(JWT),(NINT(LAT_DG(J,J1)),J=JM,J1,-INC)
      WRITE (6,905) (DASH,J=J1,JM,INC)
      J0=J1-1
         DO 40 L=1,LM+LM_REQ+1
         DO 40 J=1,JM+3
   40    XJL(J,L) = -1.E30
         KSX = 0            ! KSX = LAYERS GENERATED AT ENTRY
  100 DO 110 J=J1,JM
      DO 110 K=1,KMAX
  110 CX(J,K)=AX(J,K)*SCALE*SCALEJ(J)*SCALEK(K)
         KLMAX = KMAX+KSX
C**** HORIZONTAL SUMS AND TABLE ENTRIES
      DO 140 K=KMAX,1,-1
      AGLOB=0.
      DO 130 JHEMI=1,2
      AHEM(JHEMI)=0.
      DO 120 JH=1,JMHALF
      J=(JHEMI-1)*(JMHALF-J0)+JH+J0
      FLAT(J)=CX(J,K)/(DSJK(J,K,J1)+1.D-20)
         XJL(J,K) = FLAT(J)
         IF (DSJK(J,K,J1).EQ.0.) XJL(J,K) = -1.E30
      MLAT(J)=NINT(FLAT(J))
  120 AHEM(JHEMI)=AHEM(JHEMI)+CX(J,K)*WTJ(J,JWT,J1)
  130 AGLOB=AGLOB+AHEM(JHEMI)/JWT
      H1=AHEM(1)/(DSHEM(1,K,J1)+1.D-20)
      H2=AHEM(2)/(DSHEM(2,K,J1)+1.D-20)
      G1=AGLOB/(DSGLOB(K,J1)+1.D-20)
         XJL(JM+3,K)=H1   ! SOUTHERN HEM
         XJL(JM+2,K)=H2   ! NORTHERN HEM
         XJL(JM+1,K)=G1   ! GLOBAL
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
         XJL(J,LM+LM_REQ+1)=ASUM(J)
  180 MLAT(J)=NINT(ASUM(J)*SUMFAC)
      AGLOB=0.
      DO 200 JHEMI=1,2
      AHEM(JHEMI)=0.
      DO 190 JH=1,JMHALF
      J=(JHEMI-1)*(JMHALF-J0)+JH+J0
      AHEM(JHEMI)=AHEM(JHEMI)+ASUM(J)*WTJ(J,JWT,J1)*SUMFAC
  190 CONTINUE
  200 AGLOB=AGLOB+AHEM(JHEMI)/JWT
         XJL(JM+3,LM+LM_REQ+1)=AHEM(1)   ! SOUTHERN HEM
         XJL(JM+2,LM+LM_REQ+1)=AHEM(2)   ! NORTHERN HEM
         XJL(JM+1,LM+LM_REQ+1)=AGLOB     ! GLOBAL
         XLB=' '//acc_period(1:3)//' '//acc_period(4:12)//'  '
         TITLEO=TITLE(NT)//XLB
         in_jkmap=.true.
         nt_jk=nt
         IF(QCHECK) CALL POUT_JL(TITLEO,J1,KLMAX,XJL,PM,CLAT
     *        ,CPRES)
         in_jkmap=.false.
      WRITE (6,903) WORD(IWORD),AGLOB,AHEM(2),AHEM(1),
     *  (MLAT(J),J=JM,J1,-INC)
         IF (NT.EQ.1) CALL KEYJKT (AGLOB,ASUM)
         IF (NT.EQ.18.OR.NT.EQ.19) CALL KEYJKE (NT,AHEM,ASUM)
         IF (NT.GE.22.AND.NT.LE.33) CALL KEYJKN (NT,ASUM,SUMFAC)
      RETURN
C****
      ENTRY JKMAPS (NT,PM,AX,SCALE,SCALEJ,SCALEK,KMAX,JWT,J1,
     *  ARQX,SCALER,SCALJR,SCALLR)
         KSX = 3
         DO 205 L=1,LM+LM_REQ+1
         DO 205 J=1,JM+3
  205    XJL(J,L) = -1.E30
      LINECT=LINECT+KMAX+10
      IF (LINECT.LE.60) GO TO 230
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
      LINECT=KMAX+11
  230 J0=J1-1
C**** PRODUCE UPPER STRATOSPHERE NUMBERS FIRST
      WRITE (6,901) TITLE(NT),(DASH,J=J1,JM,INC)
      WRITE (6,904) WORD(JWT),(NINT(LAT_DG(J,J1)),J=JM,J1,-INC)
      WRITE (6,905) (DASH,J=J1,JM,INC)
      DO 260 L=LM_REQ,1,-1
      FGLOB=0.
      DO 250 JHEMI=1,2
      AHEM(JHEMI)=0.
      DO 240 JH=1,JMHALF
      J=(JHEMI-1)*(JMHALF-J0)+JH-J0
      FLATJ=ARQX(J,L)*SCALER*SCALJR(J)*SCALLR(L)
         XJL(J,L+KMAX) = FLATJ
      MLAT(J)=NINT(FLATJ)
  240 AHEM(JHEMI)=AHEM(JHEMI)+FLATJ*WTJ(J,JWT,J1)
  250 FGLOB=FGLOB+AHEM(JHEMI)/JWT
         XJL(JM+3,L+KMAX)=AHEM(1)   ! SOUTHERN HEM
         XJL(JM+2,L+KMAX)=AHEM(2)   ! NORTHERN HEM
         XJL(JM+1,L+KMAX)=FGLOB     ! GLOBAL
  260 WRITE (6,902) PM(L+LM),FGLOB,AHEM(2),AHEM(1),
     *  (MLAT(J),J=JM,J1,-INC)
      GO TO 100
  901 FORMAT ('0',30X,A64,'  CP'/1X,30('-'),24A4)
  902 FORMAT (1X,F7.3,3F8.1,1X,24I4)
  903 FORMAT (A6,2X,3F8.1,1X,24I4)
  904 FORMAT (' P(MB)   ',A4,' G      NH      SH  ',24I4)
  905 FORMAT (1X,32('-'),24A4)
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
      END SUBROUTINE JKMAP

      MODULE BDJL
!@sum  stores information for outputting lat-sigma diagnostics
!@auth M. Kelley

      IMPLICIT NONE

!@param njl_out number of jl-format output fields (not all used)
      integer, parameter :: njl_out=156

!@var title string, formed as concatentation of lname//units
      CHARACTER(LEN=64), DIMENSION(njl_out) :: TITLE
!@var units string containing output field units
      CHARACTER(LEN=50), DIMENSION(njl_out) :: UNITS
!@var lname string describing output field
      CHARACTER(LEN=50), DIMENSION(njl_out) :: LNAME
!@var sname string referencing output field in self-desc. output file
      CHARACTER(LEN=30), DIMENSION(njl_out) :: SNAME

!@var in_jlmap flag telling pout_jl to use jl titles
      logical :: in_jlmap
!@var nt_jl index telling pout_jl which field is being output
      integer :: nt_jl

      END MODULE BDJL

      SUBROUTINE JL_TITLES
      USE BDJL
      IMPLICIT NONE
      INTEGER :: K
c
      k = 0
c
      k = k + 1 ! 001
      sname(k) = 'del_u_sdrag'
      lname(k) = 'ZONAL WIND CHANGE BY STRATOSPHERIC DRAG'
      units(k) = '10**-6 M S-2'
      k = k + 1 ! 002
      sname(k) = '' !'height_jl'
      lname(k) = '' !'HEIGHT'
      units(k) = '' !'HUNDREDS OF METERS'
      k = k + 1 ! 003
      sname(k) = '' !'q_jl'
      lname(k) = '' !'SPECIFIC HUMIDITY'
      units(k) = '' !'10**-5 KG H2O/KG AIR'
      k = k + 1 ! 004
      sname(k) = '' !'rh_jl'
      lname(k) = '' !'RELATIVE HUMIDITY'
      units(k) = '' !'PERCENT'
      k = k + 1 ! 005
      sname(k) = '' !'u_jl'
      lname(k) = '' !'ZONAL WIND (U COMPONENT)'
      units(k) = '' !'TENTHS OF METERS/SECOND'
      k = k + 1 ! 006
      sname(k) = '' !'v_jl'
      lname(k) = '' !'MERIDIONAL WIND (V COMPONENT)'
      units(k) = '' !'HUNDREDTHS OF METERS/SECOND'
      k = k + 1 ! 007
      sname(k) = '' !'psi'
      lname(k) = '' !'STREAM FUNCTION'
      units(k) = '' !'10**9 KILOGRAMS/SECOND'
      k = k + 1 ! 008
      sname(k) = '' !'vvel'
      lname(k) = '' !'VERTICAL VELOCITY'
      units(k) = '' !'10**-5 MILLIBARS/SECOND'
      k = k + 1 ! 009
      sname(k) = '' !'baroc_eddy_ke_gen_jl'
      lname(k) = '' !'BAROCLINIC EDDY KINETIC ENERGY GEN.'
      units(k) = '' !'10**-1 WATTS/M**2/SIGMA'
      k = k + 1 ! 010
      sname(k) = 'mc_mflx'
      lname(k) = 'VERTICAL MASS EXCHANGE FROM MOIST CONVECTION'
      units(k) = '10**9 KG/SECOND'
      k = k + 1 ! 011
      sname(k) = 'srad_heat'
      lname(k) = 'SOLAR RADIATION HEATING RATE'
      units(k) = 'HUNDREDTHS OF DEGREES KELVIN/DAY'
      k = k + 1 ! 012
      sname(k) = 'trad_cool'
      lname(k) = 'THERMAL RADIATION COOLING RATE'
      units(k) = 'HUNDREDTHS OF DEGREES K/DAY'
      k = k + 1 ! 013
      sname(k) = 'rad_cool'
      lname(k) = 'TOTAL RADIATION COOLING RATE'
      units(k) = '10**13 WATTS/UNIT SIGMA'
      k = k + 1 ! 014
      sname(k) = 'lscond_heat'
      lname(k) = 'HEATING BY LARGE SCALE CONDENSATION'
      units(k) = '10**13 WATTS/UNIT SIGMA'
      k = k + 1 ! 015
      sname(k) = 'turb_heat'
      lname(k) = 'HEATING BY TURBULENCE'
      units(k) = '10**13 WATTS/UNIT SIGMA'
      k = k + 1 ! 016
      sname(k) = 'turb_lat'
      lname(k) = 'CHANGE OF LATENT HEAT BY TURBULENCE'
      units(k) = '10**14 W/UNIT SIGMA'
      k = k + 1 ! 017
      sname(k) = 'moist_lat'
      lname(k) = 'CHANGE OF LATENT HEAT BY MOIST CONV.'
      units(k) = '10**14 W/UNIT SIGMA'
      k = k + 1 ! 018
      sname(k) = '' !'eddy_ke'
      lname(k) = '' !'EDDY KINETIC ENERGY'
      units(k) = '' !'10**4 JOULES/M**2/UNIT SIGMA'
      k = k + 1 ! 019
      sname(k) = '' !'tot_ke'
      lname(k) = '' !'TOTAL KINETIC ENERGY'
      units(k) = '' !'10**4 JOULES/M**2/UNIT SIGMA'
      k = k + 1 ! 020
      sname(k) = 'avail_pe'
      lname(k) = 'AVAILABLE POTENTIAL ENERGY'
      units(k) = '10**5 JOULES/M**2/UNIT SIGMA'
      k = k + 1 ! 021
      sname(k) = 'DT_DYNAMICS' 
      lname(k) = 'DTEMP/DT BY DYNAMICS' 
      units(k) = '10**-1 DEG-K/DAY'
      k = k + 1 ! 022
      sname(k) = '' !'nt_dse_stand_eddy'
      lname(k) = '' !'NOR. TRANS. OF DRY STAT. ENERGY BY STAND. EDDIES'
      units(k) = '' !'10**14 W/DSIG'
      k = k + 1 ! 023
      sname(k) = '' !'nt_dse_eddy'
      lname(k) = '' !'NORTH. TRANS. OF DRY STATIC ENERGY BY EDDIES'
      units(k) = '' !'10**14 WATTS/DSIG'
      k = k + 1 ! 024
      sname(k) = '' !'tot_nt_dse'
      lname(k) = '' !'TOTAL NORTH. TRANSPORT OF DRY STATIC ENERGY'
      units(k) = '' !'10**15 WATTS/DSIG'
      k = k + 1 ! 025
      sname(k) = '' !'nt_lh_eddy'
      lname(k) = '' !'NORTHWARD TRANSPORT OF LATENT HEAT BY EDDIES'
      units(k) = '' !'10**13 WATTS/DSIG'
      k = k + 1 ! 026
      sname(k) = '' !'tot_nt_lh'
      lname(k) = '' !'TOTAL NORTHWARD TRANSPORT OF LATENT HEAT'
      units(k) = '' !'10**14 WATTS/UNIT SIG'
      k = k + 1 ! 027
      sname(k) = '' !'nt_se_eddy'
      lname(k) = '' !'NORTH.TRANSPORT OF STATIC ENERGY BY EDDIES'
      units(k) = '' !'10**14 WATTS/DSIGMA'
      k = k + 1 ! 028
      sname(k) = '' !'tot_nt_se'
      lname(k) = '' !'TOTAL NORTHWARD TRANSPORT OF STATIC ENERGY'
      units(k) = '' !'10**15 WATTS/DSIGMA'
      k = k + 1 ! 029
      sname(k) = 'DT_SDRAG' 
      lname(k) = 'DTEMP/DT BY STRATOSPHERIC DRAG' 
      units(k) = '10**-1 DEG-K/DAY' 
      k = k + 1 ! 030
      sname(k) = '' !'tot_nt_ke'
      lname(k) = '' !'TOTAL NORTHWARD TRANSPORT OF KINETIC ENERGY'
      units(k) = '' !'10**12 WATTS/DSIG'
      k = k + 1 ! 031
      sname(k) = '' !'nt_am_stand_eddy'
      lname(k) = '' !'NORTH. TRANS. OF ANG. MOMENTUM BY STAND. EDDIES'
      units(k) = '' !'10**18 J/DSIG'
      k = k + 1 ! 032
      sname(k) = '' !'nt_am_eddy'
      lname(k) = '' !'NORTH. TRANS. OF ANG. MOMENTUM BY EDDIES'
      units(k) = '' !'10**18 JOULES/DSIGMA'
      k = k + 1 ! 033
      sname(k) = '' !'tot_nt_am'
      lname(k) = '' !'TOTAL NORTHWARD TRANSPORT OF ANG. MOMENTUM'
      units(k) = '' !'10**19 JOULES/DSIG'
      k = k + 1 ! 034
      sname(k) = '' !'vt_dse_eddy'
      lname(k) = '' !'VERT. TRANS. OFDRY STATIC ENERGY BY EDDIES'
      units(k) = '' !'10**12 WATTS'
      k = k + 1 ! 035
      sname(k) = '' !'tot_vt_dse'
      lname(k) = '' !'TOT. LARGE SCALE VERT. TRANS. OF DRY STAT. ENER.'
      units(k) = '' !'10**14 WATTS'
      k = k + 1 ! 036
      sname(k) = '' !'vt_lh_eddy'
      lname(k) = '' !'VERTICAL TRANSPORT OF LATENT HEAT BY EDDIES'
      units(k) = '' !'10**12 WATTS'
      k = k + 1 ! 037
      sname(k) = '' !'tot_vt_lh'
      lname(k) = '' !'TOTAL LARGE SCALE VERT. TRANS. OF LATENT HEAT'
      units(k) = '' !'10**13 WATTS'
      k = k + 1 ! 038
      sname(k) = '' !'vt_se_eddy'
      lname(k) = '' !'VERTICAL TRANSPORT OF STATIC ENERGY BY EDDIES'
      units(k) = '' !'10**13 WATTS'
      k = k + 1 ! 039
      sname(k) = '' !'tot_vt_se'
      lname(k) = '' !'TOTAL LARGE SCALE VERT. TRANS. OF STATIC ENERGY'
      units(k) = '' !'10**14 WATTS'
      k = k + 1 ! 040
      sname(k) = 'tot_dudt'
      lname(k) = 'DU/DT   TOTAL CHANGE (CP)'
      units(k) = '10**-6 M/S/S'
      k = k + 1 ! 041
      sname(k) = '' !'tot_vt_ke'
      lname(k) = '' !'TOTAL LARGE SCALE VERT. TRANS. OF KINETIC ENERGY'
      units(k) = '' !'10**11 WATTS'
      k = k + 1 ! 042
      sname(k) = '' !'vt_am_eddy'
      lname(k) = '' !'VERT. TRANS. OFANG. MOMENTUM BY EDDIES'
      units(k) = '' !'10**16JOULES'
      k = k + 1 ! 043
      sname(k) = '' !'tot_vt_am'
      lname(k) = '' !'TOTAL LARGE SCALE VERT. TRANS. OF ANG. MOMENTUM'
      units(k) = '' !'10**18 JOULES'
      k = k + 1 ! 044
      sname(k) = 'del_am_dc'
      lname(k) = 'CHANGE OF ANG. MOMENTUM BY TURBULENCE'
      units(k) = '10**18 JOULE/UNIT SIGMA'
      k = k + 1 ! 045
      sname(k) = 'del_am_mc'
      lname(k) = 'CHANGE OF ANG. MOMENTUM BY MOIST CONV'
      units(k) = '10**18 JOULE/UNIT SIGMA'
      k = k + 1 ! 046
      sname(k) = '' !'del_am_diff'
      lname(k) = '' !'CHANGE OF ANG. MOMENTUM BY DIFFUSION'
      units(k) = '' !'10**18 JOULES/UNIT SIGMA'
      k = k + 1 ! 047
      sname(k) = 'u_epac'
      lname(k) = 'U WIND AVERAGED OVER I,EAST PACIFIC'
      units(k) = 'TENTHS OF METERS/SECOND'
      k = k + 1 ! 048
      sname(k) = 'v_epac'
      lname(k) = 'V WIND AVERAGED OVER EAST PACIFIC'
      units(k) = 'TENTHS OF METERS/SECOND'
      k = k + 1 ! 049
      sname(k) = 'vvel_epac'
      lname(k) = 'VERTICAL VELOCITY FOR EAST PACIFIC'
      units(k) = '10**-5 METERS/SECOND'
      k = k + 1 ! 050
      sname(k) = 'u_wpac'
      lname(k) = 'U WIND AVERAGED OVER I=WEST PACIFIC'
      units(k) = 'TENTHS OF METERS/SECOND'
      k = k + 1 ! 051
      sname(k) = 'v_wpac'
      lname(k) = 'V WIND AVERAGED OVER I=WEST PACIFIC'
      units(k) = 'TENTHS OF METERS/SECOND'
      k = k + 1 ! 052
      sname(k) = 'vvel_wpac'
      lname(k) = 'VERTICAL VELOCITY FOR WEST PACIFIC'
      units(k) = '10**-5 METERS/SECOND'
      k = k + 1 ! 053
      sname(k) = 'dudt_sdrag'
      lname(k) = 'DU/DT BY SDRAG'
      units(k) = '10**-6 M S-2'
      k = k + 1 ! 054
      sname(k) = 'nt_qgpv'
      lname(k) = 'NORTHWARD TRANSPORT OF Q-G POT. VORTICITY'
      units(k) = '10**18 JOULES/DSIG'
      k = k + 1 ! 055
      sname(k) = '' !'p2k_eddy_pgf'
      lname(k) = '' !'P-K BY EDDY PRESSURE GRADIENT FORCE'
      units(k) = '' !'10**-1 W/M**2/UNIT SIGMA'
      k = k + 1 ! 056
      sname(k) = '' !'del_qgpv'
      lname(k) = '' !'Q-G POT. VORTICITY CHANGE OVER LATITUDES'
      units(k) = '' !'10**-12 1/(SEC-M)'
      k = k + 1 ! 057
      sname(k) = '' !'psi_tem'
      lname(k) = '' !'TRANSFORMED STREAM FUNCTION'
      units(k) = '' !'10**9 KILOGRAMS/SECOND'
      k = k + 1 ! 058
      sname(k) = '' !'dyn_conv_eddy_geopot'
      lname(k) = '' !'DYNAMIC CONVERGENCE OF EDDY GEOPOTENTIAL'
      units(k) = '' !'.1 WATTS/M**2/DSIGMA'
      k = k + 1 ! 059
      sname(k) = 'dudt_mean_advec'
      lname(k) = 'DU/DT BY MEAN ADVECTION (CP)'
      units(k) = '10**-6 M/S/S'
      k = k + 1 ! 060
      sname(k) = 'dudt_eddy_conv'
      lname(k) = 'DU/DT BY EDDY CONVERGENCE (CP)'
      units(k) = '10**-6 M/S/S'
      k = k + 1 ! 061
      sname(k) = 'psi_cp'
      lname(k) = 'STREAM FUNCTION (CP)'
      units(k) = '10**9 KILOGRAMS/SECOND'
      k = k + 1 ! 062
      sname(k) = 'vvel'
      lname(k) = 'VERTICAL VELOCITY'
      units(k) = '10**-5 MILLIBARS/SECOND'
      k = k + 1 ! 063
      sname(k) = 'pot_vort'
      lname(k) = 'POTENTIAL VORTICITY (CP)'
      units(k) = '10**-6 K/(MB-S)'
      k = k + 1 ! 064
      sname(k) = 'dudt_advec_tem'
      lname(k) = 'DU/DT BY TRANSFORMED ADVECTION (CP)'
      units(k) = '10**-6 M/S/S'
      k = k + 1 ! 065
      sname(k) = 'dudt_epdiv'
      lname(k) = 'DU/DT BY ELIASSEN-PALM DIVERGENCE (CP)'
      units(k) = '10**-6 M/S/S'
      k = k + 1 ! 066
      sname(k) = 'dtempdt'
      lname(k) = 'DTEMP/DT   TOTAL CHANGE (CP)'
      units(k) = '10**-1 DEG-K/DAY'
      k = k + 1 ! 067
      sname(k) = 'dtempdt_mean_advec'
      lname(k) = 'DTEMP/DT BY MEAN ADVECTION (CP)'
      units(k) = '10**-1 DEG-K/DAY'
      k = k + 1 ! 068
      sname(k) = 'stdev_dp'
      lname(k) = 'STANDARD DEVIATION OF PRESSURE DIFFERENCES'
      units(k) = 'MB'
      k = k + 1 ! 069
      sname(k) = 'dtempdt_eddy_conv'
      lname(k) = 'DTEMP/DT BY EDDY CONVERGENCE (CP)'
      units(k) = '10**-1 DEG-K/DAY'
      k = k + 1 ! 070
      sname(k) = 'dtempdt_advec_tem'
      lname(k) = 'DTEMP/DT BY TRANSFORMED ADVECTION (CP)'
      units(k) = '10**-1 DEG-K/DAY'
      k = k + 1 ! 071
      sname(k) = '' !'refr_ind_wave1'
      lname(k) = '' !'REFRACTION INDEX FOR WAVE NUMBER 1'
      units(k) = '' !'10**-8 PER METER**2'
      k = k + 1 ! 072
      sname(k) = '' !'refr_ind_wave2'
      lname(k) = '' !'REFRACTION INDEX FOR WAVE NUMBER 2'
      units(k) = '' !'10**-8 PER METER**2'
      k = k + 1 ! 073
      sname(k) = '' !'refr_ind_wave3'
      lname(k) = '' !'REFRACTION INDEX FOR WAVE NUMBER 3'
      units(k) = '' !'10**-8 PER METER**2'
      k = k + 1 ! 074
      sname(k) = '' !'refr_ind_wave6'
      lname(k) = '' !'REFRACTION INDEX FOR WAVE NUMBER 6'
      units(k) = '' !'10**-8 PER METER**2'
      k = k + 1 ! 075
      sname(k) = '' !'refr_ind_wave9'
      lname(k) = '' !'REFRACTION INDEX FOR WAVE NUMBER 9'
      units(k) = '' !'10**-8 PER METER**2'
      k = k + 1 ! 076
      sname(k) = ''
      lname(k) = ''
      units(k) = ''
      k = k + 1 ! 077
      sname(k) = 'totcld'
      lname(k) = 'TOTAL CLOUD COVER'
      units(k) = 'PERCENT'
      k = k + 1 ! 078
      sname(k) = 'sscld'
      lname(k) = 'SUPER SATURATION CLOUD COVER'
      units(k) = 'PERCENT'
      k = k + 1 ! 079
      sname(k) = 'mccld'
      lname(k) = 'MOIST CONVECTIVE CLOUD COVER'
      units(k) = 'PERCENT'
      k = k + 1 ! 080
      sname(k) = 'phi_amp_wave1'
      lname(k) = 'AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 1'
      units(k) = 'METERS'
      k = k + 1 ! 081
      sname(k) = 'phi_amp_wave2'
      lname(k) = 'AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 2'
      units(k) = 'METERS'
      k = k + 1 ! 082
      sname(k) = 'phi_amp_wave3'
      lname(k) = 'AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 3'
      units(k) = 'METERS'
      k = k + 1 ! 083
      sname(k) = 'phi_amp_wave4'
      lname(k) = 'AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 4'
      units(k) = 'METERS'
      k = k + 1 ! 084
      sname(k) = 'phi_phase_wave1'
      lname(k) = 'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 1'
      units(k) = 'DEG WEST LONG'
      k = k + 1 ! 085
      sname(k) = 'phi_phase_wave2'
      lname(k) = 'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 2'
      units(k) = 'DEG WEST LONG'
      k = k + 1 ! 086
      sname(k) = 'phi_phase_wave3'
      lname(k) = 'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 3'
      units(k) = 'DEG WEST LONG'
      k = k + 1 ! 087
      sname(k) = 'phi_phase_wave4'
      lname(k) = 'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 4'
      units(k) = 'DEG WEST LONG'
      k = k + 1 ! 088
      sname(k) = '' !'nt_sensht_eddy'
      lname(k) = '' !'NORTH. TRANS. OF SENSIBLE HEAT BY EDDIES'
      units(k) = '' !'10**14 WATTS/DSIGMA'
      k = k + 1 ! 089
      sname(k) = 'npts_avg'
      lname(k) = 'NUMBER OF GRIDPOINTS INCLUDED IN AVERAGE (CP)'
      units(k) = '1'
      k = k + 1 ! 090
      sname(k) = '' !'vt_geopot_eddy'
      lname(k) = '' !'VERT. TRANS. OF GEOPOTENTIAL ENERGY BY EDDIES'
      units(k) = '' !'10**12 WATTS'
      k = k + 1 ! 091
      sname(k) = 'dp_cp'
      lname(k) = 'PRESSURE DIFFERENCES (CP)'
      units(k) = 'MB'
      k = k + 1 ! 092
      sname(k) = 'tke'
      lname(k) = 'TURBULENT KINETIC ENERGY'
      units(k) = 'W/M^2'
      k = k + 1 ! 093
      sname(k) = '' !'dyn_conv_dse'
      lname(k) = '' !'DYNAMIC CONVERGENCE OF DRY STATIC ENERGY'
      units(k) = '' !'10 WATTS/M**2/DSIGMA'
      k = k + 1 ! 094
      sname(k) = 'epflx_div'
      lname(k) = 'DIVERGENCE OF THE ELIASSEN-PALM FLUX'
      units(k) = '10**17 JOULES/UNIT SIGMA'
      k = k + 1 ! 095
      sname(k) = 'epflx_north'
      lname(k) = 'NORTHWARD ELIASSEN-PALM FLUX'
      units(k) = '10**17 JOULES/UNIT SIGMA'
      k = k + 1 ! 096
      sname(k) = 'epflx_vert'
      lname(k) = 'VERTICAL ELIASSEN-PALM FLUX'
      units(k) = '10**17 JOULES'
      k = k + 1 ! 097
      sname(k) = 'tot_dry_mc'
      lname(k) = 'TOTAL DRYING BY MOIST CONVECTION (Q2)'
      units(k) = '10**14 WATTS/DSIG'
      k = k + 1 ! 098
      sname(k) = 'tot_ht_mc'
      lname(k) = 'TOTAL HEATING BY MOIST CONVECTION (Q1)'
      units(k) = '10**14 WATTS/DSIG'
      k = k + 1 ! 099
      sname(k) = 'vt_geopot_eddy'
      lname(k) = 'VERT. TRANS. OF GEOPOTENTIAL ENERGY BY EDDIES (CP)'
      units(k) = '10**12 WATTS'
      k = k + 1 ! 100
      sname(k) = 'vt_dse_eddy'
      lname(k) = 'VERT. TRANS. OF DRY STATIC ENERGY BY EDDIES (CP)'
      units(k) = '10**12 WATTS'
      k = k + 1 ! 101
      sname(k) = 'tot_vt_dse'
      lname(k) = 'TOTAL LGE SCALE VERT. TRANS. OF DRY STAT. ENER. (CP)'
      units(k) = '10**14 WATTS'
      k = k + 1 ! 102
      sname(k) = 'vt_lh_eddy'
      lname(k) = 'VERTICAL TRANSPORT OF LATENT HEAT BY EDDIES (CP)'
      units(k) = '10**12 WATTS'
      k = k + 1 ! 103
      sname(k) = 'tot_vt_lh'
      lname(k) = 'TOTAL LARGE SCALE VERT. TRANS. OF LATENT HEAT (CP)'
      units(k) = '10**13 WATTS'
      k = k + 1 ! 104
      sname(k) = 'vt_se_eddy'
      lname(k) = 'VERTICAL TRANSPORT OF STATIC ENERGY BY EDDIES (CP)'
      units(k) = '10**13 WATTS'
      k = k + 1 ! 105
      sname(k) = 'tot_vt_se'
      lname(k) = 'TOTAL LARGE SCALE VERT. TRANS. OF STATIC ENERGY (CP)'
      units(k) = '10**14 WATTS'
      k = k + 1 ! 106
      sname(k) = 'psi_tem'
      lname(k) = 'TRANSFORMED STREAM FUNCTION (CP)'
      units(k) = '10**9 KG/SEC'
      k = k + 1 ! 107
      sname(k) = 'vt_pv'
      lname(k) = 'VERT. TRANSPORT OF POTENTIAL VORTICITY (CP)'
      units(k) = '10**4 KG-DEG K/MB/S/S'
      k = k + 1 ! 108
      sname(k) = 'vt_pv_eddy'
      lname(k) = 'VERT. TRANS. OF POT. VORT. BY EDDIES (CP)'
      units(k) = '10**4 KG-DEG K/MB/S/S'
      k = k + 1 ! 109
      sname(k) = ''
      lname(k) = ''
      units(k) = ''
      k = k + 1 ! 110
      sname(k) = 'tot_vt_ke'
      lname(k) = 'TOTAL LGE SCALE VERT. TRANS. OF KINETIC ENERGY (CP)'
      units(k) = '10**11 WATTS'
      k = k + 1 ! 111
      sname(k) = 'vt_am_eddy'
      lname(k) = 'VERT. TRANS. OF ANG. MOMENTUM BY EDDIES (CP)'
      units(k) = '10**16 JOULES'
      k = k + 1 ! 112
      sname(k) = 'tot_vt_am'
      lname(k) = 'TOTAL LGE SCALE VERT. TRANS. OF ANG. MOMENTUM (CP)'
      units(k) = '10**18 JOULES'
      k = k + 1 ! 113
      sname(k) = 'epflx_vert_cp'
      lname(k) = 'VERTICAL ELIASSEN-PALM FLUX (CP)'
      units(k) = '10**11 JOULES/METER'
      k = k + 1 ! 114
      sname(k) = ''
      lname(k) = ''
      units(k) = ''
      k = k + 1 ! 115
      sname(k) = 'stand_eddy_nt_dse_wave1'
      lname(k) = 'STAND. N. TRANS.OF D. STATIC ENERGY, WAVE #1'
      units(k) = '10**14 WATTS/DSIG'
      k = k + 1 ! 116
      sname(k) = 'stand_eddy_nt_dse_wave2'
      lname(k) = 'STAND. N. TRANS.OF D. STATIC ENERGY, WAVE #2'
      units(k) = '10**14 WATTS/DSIG'
      k = k + 1 ! 117
      sname(k) = 'stand_eddy_nt_dse_wave3'
      lname(k) = 'STAND. N. TRANS.OF D. STATIC ENERGY, WAVE #3'
      units(k) = '10**14 WATTS/DSIG'
      k = k + 1 ! 118
      sname(k) = 'stand_eddy_nt_dse_wave4'
      lname(k) = 'STAND. N. TRANS.OF D. STATIC ENERGY, WAVE #4'
      units(k) = '10**14 WATTS/DSIG'
      k = k + 1 ! 119
      sname(k) = 'stand_eddy_nt_dse_wave5'
      lname(k) = 'STAND. N. TRANS.OF D. STATIC ENERGY, WAVE #5'
      units(k) = '10**14 WATTS/DSIG'
      k = k + 1 ! 120
      sname(k) = 'trans_eddy_nt_dse_wave1'
      lname(k) = 'TRNSNT N. TRANS.OF D. STATIC ENERGY, WAVE #1'
      units(k) = '10**14 WATTS/DSIG'
      k = k + 1 ! 121
      sname(k) = 'trans_eddy_nt_dse_wave2'
      lname(k) = 'TRNSNT N. TRANS.OF D. STATIC ENERGY, WAVE #2'
      units(k) = '10**14 WATTS/DSIG'
      k = k + 1 ! 122
      sname(k) = 'trans_eddy_nt_dse_wave3'
      lname(k) = 'TRNSNT N. TRANS.OF D. STATIC ENERGY, WAVE #3'
      units(k) = '10**14 WATTS/DSIG'
      k = k + 1 ! 123
      sname(k) = 'trans_eddy_nt_dse_wave4'
      lname(k) = 'TRNSNT N. TRANS.OF D. STATIC ENERGY, WAVE #4'
      units(k) = '10**14 WATTS/DSIG'
      k = k + 1 ! 124
      sname(k) = 'trans_eddy_nt_dse_wave5'
      lname(k) = 'TRNSNT N. TRANS.OF D. STATIC ENERGY, WAVE #5'
      units(k) = '10**14 WATTS/DSIG'
      k = k + 1 ! 125
      sname(k) = 'trans_eddy_nt_dse_wave6'
      lname(k) = 'TRNSNT N. TRANS.OF D. STATIC ENERGY, WAVE #6'
      units(k) = '10**14 WATTS/DSIG'
      k = k + 1 ! 126
      sname(k) = 'trans_eddy_nt_dse_wave7'
      lname(k) = 'TRNSNT N. TRANS.OF D. STATIC ENERGY, WAVE #7'
      units(k) = '10**14 WATTS/DSIG'
      k = k + 1 ! 127
      sname(k) = 'trans_eddy_nt_dse_wave8'
      lname(k) = 'TRNSNT N. TRANS.OF D. STATIC ENERGY, WAVE #8'
      units(k) = '10**14 WATTS/DSIG'
      k = k + 1 ! 128
      sname(k) = 'trans_eddy_nt_dse_wave9'
      lname(k) = 'TRNSNT N. TRANS.OF D. STATIC ENERGY, WAVE #9'
      units(k) = '10**14 WATTS/DSIG'
      k = k + 1 ! 129
      sname(k) = 'stand_eddy_nt_lh_wave1'
      lname(k) = 'STAND. N. TRANS. OF LATENT HEAT, WAVE #1'
      units(k) = '10**13 WATTS/DSIG'
      k = k + 1 ! 130
      sname(k) = 'stand_eddy_nt_lh_wave2'
      lname(k) = 'STAND. N. TRANS. OF LATENT HEAT, WAVE #2'
      units(k) = '10**13 WATTS/DSIG'
      k = k + 1 ! 131
      sname(k) = 'stand_eddy_nt_lh_wave3'
      lname(k) = 'STAND. N. TRANS. OF LATENT HEAT, WAVE #3'
      units(k) = '10**13 WATTS/DSIG'
      k = k + 1 ! 132
      sname(k) = 'stand_eddy_nt_lh_wave4'
      lname(k) = 'STAND. N. TRANS. OF LATENT HEAT, WAVE #4'
      units(k) = '10**13 WATTS/DSIG'
      k = k + 1 ! 133
      sname(k) = 'stand_eddy_nt_lh_wave5'
      lname(k) = 'STAND. N. TRANS. OF LATENT HEAT, WAVE #5'
      units(k) = '10**13 WATTS/DSIG'
      k = k + 1 ! 134
      sname(k) = 'trans_eddy_nt_lh_wave1'
      lname(k) = 'TRANSIENT N. TRANS. OF LATENT HEAT, WAVE #1'
      units(k) = '10**13 WATTS/DSIG'
      k = k + 1 ! 135
      sname(k) = 'trans_eddy_nt_lh_wave2'
      lname(k) = 'TRANSIENT N. TRANS. OF LATENT HEAT, WAVE #2'
      units(k) = '10**13 WATTS/DSIG'
      k = k + 1 ! 136
      sname(k) = 'trans_eddy_nt_lh_wave3'
      lname(k) = 'TRANSIENT N. TRANS. OF LATENT HEAT, WAVE #3'
      units(k) = '10**13 WATTS/DSIG'
      k = k + 1 ! 137
      sname(k) = 'trans_eddy_nt_lh_wave4'
      lname(k) = 'TRANSIENT N. TRANS. OF LATENT HEAT, WAVE #4'
      units(k) = '10**13 WATTS/DSIG'
      k = k + 1 ! 138
      sname(k) = 'trans_eddy_nt_lh_wave5'
      lname(k) = 'TRANSIENT N. TRANS. OF LATENT HEAT, WAVE #5'
      units(k) = '10**13 WATTS/DSIG'
      k = k + 1 ! 139
      sname(k) = 'trans_eddy_nt_lh_wave6'
      lname(k) = 'TRANSIENT N. TRANS. OF LATENT HEAT, WAVE #6'
      units(k) = '10**13 WATTS/DSIG'
      k = k + 1 ! 140
      sname(k) = 'trans_eddy_nt_lh_wave7'
      lname(k) = 'TRANSIENT N. TRANS. OF LATENT HEAT, WAVE #7'
      units(k) = '10**13 WATTS/DSIG'
      k = k + 1 ! 141
      sname(k) = 'trans_eddy_nt_lh_wave8'
      lname(k) = 'TRANSIENT N. TRANS. OF LATENT HEAT, WAVE #8'
      units(k) = '10**13 WATTS/DSIG'
      k = k + 1 ! 142
      sname(k) = 'trans_eddy_nt_lh_wave9'
      lname(k) = 'TRANSIENT N. TRANS. OF LATENT HEAT, WAVE #9'
      units(k) = '10**13 WATTS/DSIG'
      k = k + 1 ! 143
      sname(k) = 'stand_eddy_nt_am_wave1'
      lname(k) = 'STAND. N. TRANS. OF ANG. MOMENT., WAVE #1'
      units(k) = '10**18 JOULES/DSIG'
      k = k + 1 ! 144
      sname(k) = 'stand_eddy_nt_am_wave2'
      lname(k) = 'STAND. N. TRANS. OF ANG. MOMENT., WAVE #2'
      units(k) = '10**18 JOULES/DSIG'
      k = k + 1 ! 145
      sname(k) = 'stand_eddy_nt_am_wave3'
      lname(k) = 'STAND. N. TRANS. OF ANG. MOMENT., WAVE #3'
      units(k) = '10**18 JOULES/DSIG'
      k = k + 1 ! 146
      sname(k) = 'stand_eddy_nt_am_wave4'
      lname(k) = 'STAND. N. TRANS. OF ANG. MOMENT., WAVE #4'
      units(k) = '10**18 JOULES/DSIG'
      k = k + 1 ! 147
      sname(k) = 'stand_eddy_nt_am_wave5'
      lname(k) = 'STAND. N. TRANS. OF ANG. MOMENT., WAVE #5'
      units(k) = '10**18 JOULES/DSIG'
      k = k + 1 ! 148
      sname(k) = 'trans_eddy_nt_am_wave1'
      lname(k) = 'TRNSNT N. TRANS. OF ANG. MOMENT., WAVE #1'
      units(k) = '10**18 JOULES/DS'
      k = k + 1 ! 149
      sname(k) = 'trans_eddy_nt_am_wave2'
      lname(k) = 'TRNSNT N. TRANS. OF ANG. MOMENT., WAVE #2'
      units(k) = '10**18 JOULES/DS'
      k = k + 1 ! 150
      sname(k) = 'trans_eddy_nt_am_wave3'
      lname(k) = 'TRNSNT N. TRANS. OF ANG. MOMENT., WAVE #3'
      units(k) = '10**18 JOULES/DS'
      k = k + 1 ! 151
      sname(k) = 'trans_eddy_nt_am_wave4'
      lname(k) = 'TRNSNT N. TRANS. OF ANG. MOMENT., WAVE #4'
      units(k) = '10**18 JOULES/DS'
      k = k + 1 ! 152
      sname(k) = 'trans_eddy_nt_am_wave5'
      lname(k) = 'TRNSNT N. TRANS. OF ANG. MOMENT., WAVE #5'
      units(k) = '10**18 JOULES/DS'
      k = k + 1 ! 153
      sname(k) = 'trans_eddy_nt_am_wave6'
      lname(k) = 'TRNSNT N. TRANS. OF ANG. MOMENT., WAVE #6'
      units(k) = '10**18 JOULES/DS'
      k = k + 1 ! 154
      sname(k) = 'trans_eddy_nt_am_wave7'
      lname(k) = 'TRNSNT N. TRANS. OF ANG. MOMENT., WAVE #7'
      units(k) = '10**18 JOULES/DS'
      k = k + 1 ! 155
      sname(k) = 'trans_eddy_nt_am_wave8'
      lname(k) = 'TRNSNT N. TRANS. OF ANG. MOMENT., WAVE #8'
      units(k) = '10**18 JOULES/DS'
      k = k + 1 ! 156
      sname(k) = 'trans_eddy_nt_am_wave9'
      lname(k) = 'TRNSNT N. TRANS. OF ANG. MOMENT., WAVE #9'
      units(k) = '10**18 JOULES/DS'

c create titles by concatenating long names with units
c no checks whether total length of lname+units exceeds length of title
      do k=1,njl_out
         title(k)=''
         if(lname(k).ne.'')
     &        title(k) = trim(lname(k))//' ('//trim(units(k))//')'
      enddo

      RETURN
      END SUBROUTINE JL_TITLES

      SUBROUTINE DIAGJL
c      USE PRTCOM, only :
      USE CONSTANT, only :
     &     grav,rgas,kapa,sday,twopi,sha,bygrav
      USE MODEL_COM, only :
     &     im,jm,lm,FIM,
     &     BYDSIG,BYIM,DSIG,DTsrc,IDACC,IMH,
     &     PTOP,PMTOP,SIG,SIGE,JHOUR,PSFMPT
      USE GEOM, only :
     &     AREAG,BYDXYP,COSP,COSV,DLON,DXV,DXYP,DXYV,FCOR,RADIUS,WTJ
      USE DAGPCOM, only :
     &     PLE,PLM
      USE DAGCOM, only :
     &     ajl,apj,asjl,kdiag,aij,LM_REQ, qcheck,ij_phi1k
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(JM) ::
     &     BYP,BYPV,BYDAPO,BYPDA,DXCOSV,DACOSV,DXYPPO,ONES,ONESPO
      DOUBLE PRECISION, DIMENSION(JM,LM) :: AX,BX,CX,DX
      DOUBLE PRECISION, DIMENSION(JM,LM_REQ) :: ARQX
      DOUBLE PRECISION, DIMENSION(LM_REQ) :: BYDPS,BYPKS
      DOUBLE PRECISION, DIMENSION(0:IMH) :: AN,BN
      DOUBLE PRECISION, DIMENSION(JM,8,4) :: AMPLTD,PHASE

      INTEGER :: LINECT,JMHALF,INC
      COMMON/DJLCOM/LINECT,JMHALF,INC

      DOUBLE PRECISION, DIMENSION(LM) :: PKM,BYD2SG
      DOUBLE PRECISION, DIMENSION(LM+LM_REQ) :: PL

      DOUBLE PRECISION, DIMENSION(7), PARAMETER ::
     &     PMB=(/999.9,850.,700.,500.,300.,100.,30./)
      DOUBLE PRECISION, PARAMETER :: P1000=1000.

      INTEGER :: J,K,KM,L,LDN,LS,LUP,N

      DOUBLE PRECISION ::
     &     BDN,BUP,BY100G,BYIACN,BYIADA,BYIARD,BYIMDA,
     &     DAM4,DXCVN,DXCVS,ELOFIM,FIMDA,P1000K,
     &     SCALE,SCALE2,SCALES,SCALEV,XWON

C**** INITIALIZE CERTAIN QUANTITIES
      XWON=TWOPI/(DLON*FIM)
      INC=1+(JM-1)/24
      JMHALF=JM/2
      BY100G=.01*BYGRAV
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
C     SCALE=25.D-9*BYIADA*BYGRAV
C     CALL JLMAP (7,PLE,AX,SCALE,DXV,ONES,LM,2,2)
C     CALL JLMAP (57,PLE,BX,SCALE,DXV,ONES,LM,2,2)
C**** VERTICAL VELOCITY AND MASS FLUX MOIST CONVECTION
C     SCALE=-1.D5*BYIMDA
C     CALL JLMAP (8,PLE,AJL(1,1,6),SCALE,BYDAPO,ONES,LM-1,2,1)
      SCALE=100.D-9*XWON*BYIACN/(GRAV*DTsrc)
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
      SCALES=SCALE*PSFMPT
      CALL JLMAPS (13,PLM,AX,SCALE,DXYPPO,BYDSIG,LM,1,1,
     *  ARQX,SCALES,DXYPPO,BYDPS)
C**** TOTAL, SUPER SATURATION, AND CONVECTIVE CLOUD COVER
      SCALE=100.*BYIARD*BYIM
      CALL JLMAP (77,PLM,AJL(1,1,19),SCALE,ONESPO,ONES,LM,2,1)
      CALL JLMAP (78,PLM,AJL(1,1,28),SCALE,ONESPO,ONES,LM,2,1)
      CALL JLMAP (79,PLM,AJL(1,1,29),SCALE,ONESPO,ONES,LM,2,1)
C**** TURBULENT KINETIC ENERGY
      SCALE=BYIACN
      CALL JLMAP (92,PLM,AJL(1,1,54),SCALE,ONES,ONES,LM,2,1)
C**** HEATING BY LARGE SCALE COND., MOIST CONVECTION AND TURBULENCE
      SCALE=100.D-13*XWON*SHA*BYIACN/(GRAV*DTsrc)
      CALL JLMAP (14,PLM,AJL(1,1,11),SCALE,DXYPPO,BYDSIG,LM,1,1)
      CALL JLMAP (15,PLM,AJL(1,1,12),SCALE,DXYPPO,ONES,LM,1,1)
      SCALE=0.1*SCALE
      CALL JLMAP (16,PLM,AJL(1,1,55),SCALE,DXYPPO,BYDSIG,LM,1,1)
      CALL JLMAP (17,PLM,AJL(1,1,53),SCALE,DXYPPO,BYDSIG,LM,1,1)
      CALL JLMAP (98,PLM,AJL(1,1,56),SCALE,DXYPPO,ONES,LM,1,1)
      CALL JLMAP (97,PLM,AJL(1,1,57),SCALE,DXYPPO,ONES,LM,1,1)
C****
C**** ENERGY
C****
C**** AVAILABLE POTENTIAL ENERGY
      SCALE=50.D-5*RGAS*BYIMDA*BYGRAV
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
      SCALE=25.D-18*XWON*BYIADA*RADIUS*BYGRAV * 0. ! AJL47 not done
      CALL JLMAP (54,PLM,CX,SCALE,ONES,ONES,LM,1,1)
C****
C**** VERTICAL TRANSPORTS
C****
C**** VERTICAL TRANSPORT OF ANGULAR MOMENTUM BY SMALL SCALE MOTIONS
      SCALE=100.D-18*XWON*RADIUS*BYIACN/(GRAV*DTsrc)
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
      SCALE=100.D-17*XWON*BYIADA*RADIUS*BYGRAV
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
      CALL FFT (AIJ(1,J,IJ_PHI1K-1+K),AN,BN)
      DO 620 N=1,4
      AMPLTD(J,K,N)=SQRT(AN(N)*AN(N)+BN(N)*BN(N))
      PHASE(J,K,N)=(ATAN2(BN(N),AN(N))-TWOPI)/N+ELOFIM
      IF (PHASE(J,K,N).LE.-.5*TWOPI) PHASE(J,K,N)=PHASE(J,K,N)+TWOPI
      PHASE(J,K,N)=-PHASE(J,K,N)
  620 CONTINUE
      SCALE=BYIADA*BYGRAV
      DO 630 N=1,4
  630 CALL JLMAP (N+79,PMB,AMPLTD(1,1,N),SCALE,ONES,ONES,KM,2,1)
      SCALE=360./TWOPI
      DO 640 N=1,4
  640 CALL JLMAP (N+83,PMB,PHASE(1,1,N),SCALE,ONES,ONES,KM,2,1)
      if(qcheck) call close_jl
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
      USE DAGCOM, only : QCHECK,acc_period,iu_jl,LM_REQ
      USE MODEL_COM, only :
     &     jm,lm,DSIG,JDATE,JDATE0,AMON,AMON0,JYEAR,JYEAR0,SIGE,XLABEL
      USE GEOM, only :
     &     LAT_DG,WTJ
      USE BDJL, only :
     &     title,in_jlmap,nt_jl
      IMPLICIT NONE

      INTEGER, DIMENSION(JM) :: MLAT
      DOUBLE PRECISION, DIMENSION(JM) :: FLAT,ASUM
      DOUBLE PRECISION, DIMENSION(2) :: FHEM,HSUM

      INTEGER :: LINECT,JMHALF,INC
      COMMON/DJLCOM/LINECT,JMHALF,INC

      INTEGER :: J1,JWT,LMAX,NT
      DOUBLE PRECISION :: SCALE,SCALER
      DOUBLE PRECISION, DIMENSION(JM,LM) :: AX
      DOUBLE PRECISION, DIMENSION(JM,LM_REQ) :: ARQX
      DOUBLE PRECISION, DIMENSION(JM) :: SCALEJ,SCALJR
      DOUBLE PRECISION, DIMENSION(LM) :: SCALEL
      DOUBLE PRECISION, DIMENSION(LM_REQ) :: SCALLR
      DOUBLE PRECISION, DIMENSION(LM+LM_REQ) :: PL

      CHARACTER*4 DASH,WORD(4)
      DATA DASH/'----'/,WORD/'SUM','MEAN',' ','.1*'/

      INTEGER :: IWORD,J,J0,JH,JHEMI,K,L  ,ksx,klmax
      DOUBLE PRECISION :: FGLOB,FLATJ,GSUM,SDSIG,SUMFAC

      REAL*8, DIMENSION(JM+3,LM+LM_REQ+1) :: XJL ! for binary output
      CHARACTER XLB*16,CLAT*16,CPRES*16,CBLANK*16,TITLEO*80
      DATA CLAT/'LATITUDE'/,CPRES/'PRESSURE (MB)'/,CBLANK/' '/
C****
C**** PRODUCE A LATITUDE BY LAYER TABLE OF THE ARRAY A
C****
   10 LINECT=LINECT+LMAX+7
      IF (LINECT.LE.60) GO TO 20
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
      LINECT=LMAX+8
   20 WRITE (6,901) TITLE(NT),(DASH,J=J1,JM,INC)
      WRITE (6,904) WORD(JWT),(NINT(LAT_DG(J,J1)),J=JM,J1,-INC)
      WRITE (6,905) (DASH,J=J1,JM,INC)
      J0=J1-1
         DO 40 L=1,LM+LM_REQ+1
         DO 40 J=1,JM+3
   40    XJL(J,L) = -1.E30
         KSX = 0            ! KSX = LAYERS GENERATED AT ENTRY
  100 SDSIG=1.-SIGE(LMAX+1)
         KLMAX = LMAX+KSX
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
         XJL(J,L) = FLAT(J)
      MLAT(J)=NINT(FLAT(J))
  115 ASUM(J)=ASUM(J)+FLAT(J)*DSIG(L)/SDSIG
  120 FHEM(JHEMI)=FHEM(JHEMI)+FLAT(J)*WTJ(J,JWT,J1)
  130 FGLOB=FGLOB+FHEM(JHEMI)/JWT
         XJL(JM+3,L)=FHEM(1)   ! SOUTHERN HEM
         XJL(JM+2,L)=FHEM(2)   ! NORTHERN HEM
         XJL(JM+1,L)=FGLOB     ! GLOBAL
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
         DO 180 J=J1,JM
  180    XJL(J   ,LM+LM_REQ+1)=ASUM(J)
         XJL(JM+3,LM+LM_REQ+1)=HSUM(1)   ! SOUTHERN HEM
         XJL(JM+2,LM+LM_REQ+1)=HSUM(2)   ! NORTHERN HEM
         XJL(JM+1,LM+LM_REQ+1)=GSUM      ! GLOBAL
         XLB=' '//acc_period(1:3)//' '//acc_period(4:12)//'  '
         TITLEO=TITLE(NT)//XLB
         in_jlmap=.true.
         nt_jl=nt
         IF(QCHECK) CALL POUT_JL(TITLEO,J1,KLMAX,XJL,PL,CLAT
     *        ,CPRES)
         in_jlmap=.false.
      WRITE (6,903) WORD(IWORD),GSUM,HSUM(2),HSUM(1),
     *  (MLAT(J),J=JM,J1,-INC)
C        IF (NT.EQ.1) CALL KEYJLT (GSUM,ASUM)
C        IF (NT.EQ.18.OR.NT.EQ.19) CALL KEYJLE (NT,HSUM,ASUM)
C        IF (NT.GE.22.AND.NT.LE.33) CALL KEYJLN (NT,ASUM,SUMFAC)
      RETURN
C****
      ENTRY JLMAPS (NT,PL,AX,SCALE,SCALEJ,SCALEL,LMAX,JWT,J1,
     *  ARQX,SCALER,SCALJR,SCALLR)
         KSX = 3
         DO 205 L=1,LM+3
         DO 205 J=1,JM
  205    XJL(J,L) = -1.E30
      LINECT=LINECT+LMAX+10
      IF (LINECT.LE.60) GO TO 200
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
      LINECT=LMAX+11
  200 J0=J1-1
C**** PRODUCE UPPER STRATOSPHERE NUMBERS FIRST
      WRITE (6,901) TITLE(NT),(DASH,J=J1,JM,INC)
      WRITE (6,904) WORD(JWT),(NINT(LAT_DG(J,J1)),J=JM,J1,-INC)
      WRITE (6,905) (DASH,J=J1,JM,INC)
      DO 230 L=LM_REQ,1,-1
      FGLOB=0.
      DO 220 JHEMI=1,2
      FHEM(JHEMI)=0.
      DO 210 JH=1,JMHALF
      J=(JHEMI-1)*(JMHALF-J0)+JH-J0
      FLATJ=ARQX(J,L)*SCALER*SCALJR(J)*SCALLR(L)
         XJL(J,L+LMAX) = FLATJ
      MLAT(J)=NINT(FLATJ)
  210 FHEM(JHEMI)=FHEM(JHEMI)+FLATJ*WTJ(J,JWT,J1)
  220 FGLOB=FGLOB+FHEM(JHEMI)/JWT
         XJL(JM+3,L+LMAX)=FHEM(1)   ! SOUTHERN HEM
         XJL(JM+2,L+LMAX)=FHEM(2)   ! NORTHERN HEM
         XJL(JM+1,L+LMAX)=FGLOB     ! GLOBAL
  230 WRITE (6,902) PL(L+LM),FGLOB,FHEM(2),FHEM(1),
     *  (MLAT(J),J=JM,J1,-INC)
      GO TO 100
  901 FORMAT ('0',30X,A64/2X,32('-'),24A4)
  902 FORMAT (1X,F8.3,3F8.1,1X,24I4)
  903 FORMAT (1X,A6,2X,3F8.1,1X,24I4)
  904 FORMAT ('  P(MB)   ',A4,' G      NH      SH  ',24I4)
  905 FORMAT (2X,32('-'),24A4)
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
      END SUBROUTINE JLMAP

      MODULE BDIL
!@sum  stores information for outputting lon-height diagnostics
!@auth M. Kelley

      IMPLICIT NONE

!@param nil_out number of il-format output fields
      integer, parameter :: nil_out=16

!@var title string, formed as concatentation of lname//units
      CHARACTER(LEN=64), DIMENSION(nil_out) :: TITLE
!@var units string containing output field units
      CHARACTER(LEN=50), DIMENSION(nil_out) :: UNITS
!@var lname string describing output field
      CHARACTER(LEN=50), DIMENSION(nil_out) :: LNAME
!@var sname string referencing output field in self-desc. output file
      CHARACTER(LEN=30), DIMENSION(nil_out) :: SNAME

!@var nt_il index telling pout_il which field is being output
      integer :: nt_il

      END MODULE BDIL

      SUBROUTINE IL_TITLES
      USE BDIL
      IMPLICIT NONE
      INTEGER :: K
c
      k = 0
c
      k = k + 1
      sname(k) = 'u_equator'
      lname(k) = 'ZONAL WIND (U COMPONENT) AROUND +/- 5 DEG'
      units(k) = 'METERS/SECOND'
      k = k + 1
      sname(k) = 'v_equator'
      lname(k) = 'MERIDIONAL WIND (V COMPONENT) AROUND +/- 5 DEG'
      units(k) = 'METERS/SECOND'
      k = k + 1
      sname(k) = 'vvel_equator'
      lname(k) = 'VERTICAL VELOCITY AROUND +/- 5 DEG'
      units(k) = '10**-4 METERS/SECOND'
      k = k + 1
      sname(k) = 'temp_equator'
      lname(k) = 'TEMPERATURE AROUND +/- 5 DEG'
      units(k) = 'DEGREES CENTIGRADE'
      k = k + 1
      sname(k) = 'rh_equator'
      lname(k) = 'RELATIVE HUMIDITY AROUND +/- 5 DEG'
      units(k) = 'PERCENT'
      k = k + 1
      sname(k) = 'mcheat_equator'
      lname(k) = 'MOIST CONVECTIVE HEATING AROUND +/- 5 DEG'
      units(k) = '10**13 WATTS/DSIG'
      k = k + 1
      sname(k) = 'rad_cool_equator'
      lname(k) = 'TOTAL RADIATIVE COOLING AROUND +/- 5 DEG'
      units(k) = '10**13 WATTS/DSIG'
      k = k + 1
      sname(k) = ''
      lname(k) = ''
      units(k) = ''
      k = k + 1
      sname(k) = 'vvel_50N'
      lname(k) = 'VERTICAL VELOCITY AT 50 N'
      units(k) = '10**-4 METERS/SECOND'
      k = k + 1
      sname(k) = 'temp_50N'
      lname(k) = 'TEMPERATURE AT 50 N'
      units(k) = 'DEGREES CENTIGRADE'
      k = k + 1
      sname(k) = 'rad_cool_50N'
      lname(k) = 'TOTAL RADIATIVE COOLING AT 50 N'
      units(k) = '10**13 WATTS/UNIT SIGMA'
      k = k + 1
      sname(k) = 'u_50N'
      lname(k) = 'ZONAL WIND AT 50 N'
      units(k) = 'METERS/SECOND'
      k = k + 1
      sname(k) = 'vvel_70N'
      lname(k) = 'VERTICAL VELOCITY AT 70 N'
      units(k) = '10**-4 METERS/SECOND'
      k = k + 1
      sname(k) = 'temp_70N'
      lname(k) = 'TEMPERATURE AT 70 N'
      units(k) = 'DEGREES CENTIGRADE'
      k = k + 1
      sname(k) = 'rad_cool_70N'
      lname(k) = 'TOTAL RADIATIVE COOLING AT 70 N'
      units(k) = '10**13 WATTS/UNIT SIGMA'
      k = k + 1
      sname(k) = 'u_70N'
      lname(k) = 'ZONAL WIND AT 70 N'
      units(k) = 'METERS/SECOND'

c create titles by concatenating long names with units
      do k=1,nil_out
         title(k)=''
         if(lname(k).ne.'')
     &        title(k) = trim(lname(k))//' ('//trim(units(k))//')'
      enddo

      return
      END SUBROUTINE IL_TITLES

      SUBROUTINE DIAGIL
      USE CONSTANT, only : grav,rgas,sha,bygrav
      USE MODEL_COM, only : im,jm,lm,fim,bydsig,dsig,dtsrc,idacc,jeq,psf
     *     ,ptop,sig,sige,xlabel,lrunid
      USE GEOM, only :  areag,dxyp,dxyv
      USE DAGPCOM, only : plm,ple
      USE DAGCOM, only : ail,lm_req,j50n,j70n,acc_period, qcheck
      USE BDIL, only :
     &     title,nt_il
      IMPLICIT NONE

      INTEGER :: LINECT,JMHALF,INC
      COMMON/DILCOM/LINECT,JMHALF,INC
      DOUBLE PRECISION, DIMENSION(LM) :: ONES
      DOUBLE PRECISION, DIMENSION(LM+LM_REQ) :: PL
      INTEGER :: J,L,K
      DOUBLE PRECISION :: BYIACN,BYIADA,BYIARD,SCALE
      INTEGER, PARAMETER :: K_PIL=16
      INTEGER, DIMENSION(K_PIL) ::
     *  KNDEX=(/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/),
     *  CASEIL=(/1,1,2,3,3, 4,5,0,2,3, 0,1,2,3,0, 1/)
c    *  CASEIL=(/1,1,2,3,3, 4,5,0,2,3, 5,1,2,3,5, 1/)
      DOUBLE PRECISION SIL(K_PIL)

C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      IF(QCHECK) call open_il(trim(acc_period)//'.il'//XLABEL(1:LRUNID))

C**** INITIALIZE CERTAIN QUANTITIES
      call il_titles
      INC=1+(JM-1)/24
      JMHALF=JM/2
      DO 20 L=1,LM
      ONES(L)=1.
   20 CONTINUE
      LINECT=65
      BYIACN=1./(IDACC(1)+1.D-20)
      BYIARD=1./(IDACC(2)+1.D-20)
      BYIADA=1./(IDACC(4)+1.D-20)
      SIL( 1)=BYIADA/3.
      SIL( 2)=BYIADA/3.
      SIL( 3)=
     *  -1.D4*BYIADA*RGAS*BYGRAV/(DXYP(JEQ)+DXYP(JEQ-1)+DXYP(JEQ-2))
      SIL( 4)=BYIADA/3.
      SIL( 5)=1.D2*BYIADA/3.
      SIL( 6)=100.D-13*SHA*BYIACN/(GRAV*DTsrc)
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
      nt_il = k
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
      if(qcheck) call close_il
      RETURN
      END SUBROUTINE DIAGIL

      SUBROUTINE ILMAP (TITLE,PL,AX,SCALE,SCALEL,LMAX,JWT,ISHIFT)
      USE DAGCOM, only : qcheck,acc_period,iu_il
      USE CONSTANT, only : twopi
      USE MODEL_COM, only : im,jm,lm,dsig,jdate,jdate0,amon,amon0,jyear
     *     ,jyear0,sige,xlabel
      USE GEOM, only : dlon,lon_dg
      IMPLICIT NONE

      REAL*8 :: XIL(IM,LM),ZONAL(LM) ! used for post-proc
      CHARACTER XLB*80,CWORD*8
      CHARACTER*16, PARAMETER :: CBLANK=' ',
     &  CLAT='LONGITUDE',CPRES='PRESSURE (MB)'

      DOUBLE PRECISION :: FGLOB,FLON,GSUM,SDSIG, SCALE
      INTEGER :: I,K,L, LMAX,JWT,ISHIFT

      INTEGER, DIMENSION(IM) :: MLON
      DOUBLE PRECISION, DIMENSION(IM) :: ASUM

      INTEGER :: LINECT,JMHALF,INC
      COMMON/DILCOM/LINECT,JMHALF,INC

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
        WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
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
         XLB(65:80)=' '//acc_period(1:3)//' '//acc_period(4:12)//'  '
         IF(QCHECK) CALL POUT_IL(XLB,ISHIFT,LMAX,XIL,PL,CLAT,CPRES
     *        ,ASUM,GSUM,ZONAL)

      RETURN
  901 FORMAT ('0',30X,A64/1X,14('-'),36A3)
  902 FORMAT (F6.1,F8.1,1X,36I3)
  903 FORMAT (F14.1,1X,36I3)
  904 FORMAT (' P(MB)',4X,A4,1X,36I3)  ! U-grid (i.e., centers)
  905 FORMAT (1X,14('-'),36A3)
  906 FORMAT (' P(MB)',4X,A4,36I3)     ! V-grid (i.e., edges)
  907 FORMAT ('1',27A4,I4,1X,A3,I3,' TO ',I3,1X,A3,I5)
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

      SUBROUTINE DIAG7P
C****
C**** THIS ENTRY PRINTS THE TABLES
C****
c      USE PRTCOM, only :
      USE DAGCOM, only :
     &     nwav_dag,wave,Max12HR_sequ,Min12HR_sequ
      USE MODEL_COM, only :
     &     im,IDACC,JDATE,JDATE0,AMON,AMON0,JYEAR,JYEAR0,XLABEL
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
     &     IDACC9,K,KPAGE,KQ,KTABLE,
     &     M,MMAXP1,N,NMAX,NS,NUA,NX

      NMAX=NWAV_DAG
      IDACC9=IDACC(9)
      IF (IDACC9.LE.MMAX) RETURN
C**** PATCH NEEDED IF SEVERAL RESTART FILES WERE ACCUMULATED
      IF (IDACC(12).LE.1) GO TO 320
      IDACC9=Min12HR_sequ           ! in case a February was included
      BYIA12=1./IDACC(12)
      WAVE(:,:,:,:)=WAVE(:,:,:,:)*BYIA12
  320 CONTINUE
      IF (IDACC9.GT.Max12HR_sequ) IDACC9=Max12HR_sequ
C****
C**** OUTPUT WAVE POWER AT THE EQUATOR
C****
      MMAXP1=MMAX+1
      DO 400 KPAGE=1,2
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
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
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
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
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
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

      MODULE BDIJ
!@sum  stores information for outputting lon-lat diagnostics
!@auth M. Kelley

      IMPLICIT NONE

!@param nij_out number of ij-format output fields
      integer, parameter :: nij_out=67

!@var title string, formed as concatentation of lname//units
      CHARACTER(LEN=48), DIMENSION(nij_out) :: TITLE
!@var units string containing output field units
      CHARACTER(LEN=50), DIMENSION(nij_out) :: UNITS
!@var lname string describing output field
      CHARACTER(LEN=50), DIMENSION(nij_out) :: LNAME
!@var sname string referencing output field in self-desc. output file
      CHARACTER(LEN=30), DIMENSION(nij_out) :: SNAME

!@var nt_ij index telling pout_ij which field is being output
      integer :: nt_ij

!@param LEGEND "contour levels" for ij-maps
      CHARACTER(LEN=40), DIMENSION(24), PARAMETER :: LEGEND=(/
     *  '0=0,1=5...9=45,A=50...K=100             ',
     *  '0=0...9=90,A=100...I=180...R=270        ',
     *  '1=.5...9=4.5,A=5...Z=17.5,+=MORE        ',
     *  '1=.1...9=.9,A=1...Z=3.5,+=MORE          ',
     *  '1=2...9=18,A=20...Z=70,+=MORE           ',
     *  '1=50...9=450,A=500...Z=1750,+=MORE      ',
     *  '1=100...9=900,A=1000...Z=3500,+=MORE    ',
     *  '1=20...9=180,A=200...Z=700,+=MORE       ',
     *  'A=1...Z=26,3=30...9=90,+=100-150,*=MORE ',
     *  '0=0,A=.1...Z=2.6,3=3...9=9,+=10-15      ',
     *  '-=LESS,Z=-78...0=0...9=27,+=MORE        ',
     *  '-=LESS,Z=-260...0=0...9=90,+=MORE       ',
     *  '-=LESS,Z=-520...0=0...9=180,+=MORE      ',
     *  '-=LESS,Z=-1300...0=0...9=450,+=MORE     ',
     *  '-=LESS,Z=-2600...0=0...9=900,+=MORE     ',
     *  '-=LESS,Z=-3900...0=0...9=1350,+=MORE    ',
     *  '-=LESS,Z=-5200...0=0...9=1800,+=MORE    ',
     *  '-=LESS,9=-.9...0=0,A=.1...Z=2.6,+=MORE  ',
     *  '-=LESS,9=-45...0=0,A=5...K=45...+=MORE  ',
     *  '-=LESS,9=-90...0=0,A=10...Z=260,+=MORE  ',
     *  '-=LESS,9=-180...A=20...Z=520,+=MORE     ',
     *  '-=LESS,9=-9...0=0,A=1...Z=26,+=MORE     ',
     *  '-=LESS,9=-36...0=0,A=4...Z=104,+=MOR    ',
     *  '1=5...9=45,A=50...Z=175,+=MORE          ' /)

!@param [ABCDE]CHAR "color bars" for ij-maps
      CHARACTER(LEN=38), PARAMETER :: ACHAR=
     &     ' 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ+'
      CHARACTER(LEN=23), PARAMETER :: BCHAR=
     &     ' 0123456789ABCDEFGHIJKX'
      CHARACTER(LEN=38), PARAMETER :: CCHAR=
     &     '-9876543210ABCDEFGHIJKLMNOPQRSTUVWXYZ+'
      CHARACTER(LEN=37), PARAMETER :: DCHAR=
     &     ' 0ABCDEFGHIJKLMNOPQRSTUVWXYZ3456789+*'
      CHARACTER(LEN=38), PARAMETER :: ECHAR=
     &     '-ZYXWVUTSRQPONMLKJIHGFEDCBA0123456789+'

      END MODULE BDIJ

      SUBROUTINE IJ_TITLES
      USE BDIJ
      IMPLICIT NONE
      INTEGER :: k,row,col
c
      k = 0
c
      k = k + 1
      sname(k) = 'topog'
      lname(k) = 'TOPOGRAPHY'
      units(k) = 'METERS'
      k = k + 1
      sname(k) = 'frac_land'
      lname(k) = 'LAND COVERAGE'
      units(k) = '%'
      k = k + 1
      sname(k) = 'frac_ocnice'
      lname(k) = 'OCEAN ICE COVERAGE'
      units(k) = '%'
      k = k + 1
      sname(k) = 'frac_snow'
      lname(k) = 'SNOW COVERAGE'
      units(k) = '%'
      k = k + 1
      sname(k) = 'snowdp'
      lname(k) = 'SNOW DEPTH'
      units(k) = 'MM H2O'
      k = k + 1
      sname(k) = 'frac_snow_and_ice'
      lname(k) = 'SNOW AND ICE COVERAGE'
      units(k) = '%'
      k = k + 1
      sname(k) = 'prec'
      lname(k) = 'PRECIPITATION'
      units(k) = 'MM/DAY'
      k = k + 1
      sname(k) = 'evap'
      lname(k) = 'EVAPORATION'
      units(k) = 'MM/DAY'
      k = k + 1
      sname(k) = 'snsht_flx'
      lname(k) = 'SENSIBLE HEAT FLUX'
      units(k) = 'WATTS/M**2'
      k = k + 1
      sname(k) = 'grnd_wetness'
      lname(k) = 'GROUND WETNESS'
      units(k) = '%'
      k = k + 1
      sname(k) = 'grnd_runoff'
      lname(k) = 'GROUND RUNOFF'
      units(k) = 'MM/DAY'
      k = k + 1
      sname(k) = 'grnd_temp'
      lname(k) = 'GROUND TEMPERATURE'
      units(k) = 'DEGREES C'
      k = k + 1
      sname(k) = 'srf_wspeed'
      lname(k) = 'SURFACE WIND SPEED'
      units(k) = 'METERS/SEC'
      k = k + 1
      sname(k) = 'jet_speed'
      lname(k) = 'JET SPEED'
      units(k) = 'METERS/SEC'
      k = k + 1
      sname(k) = 'srf_wspeed_bgrid'
      lname(k) = 'SURF WIND SPEED FROM U,V'
      units(k) = 'M/S'
      k = k + 1
      sname(k) = 'mtn_wave_mom_flx'
      lname(k) = 'MTN WAVE MOM. FLUX'
      units(k) = 'D/CM**2'
      k = k + 1
      sname(k) = 'jet_dir'
      lname(k) = 'JET DIRECTION'
      units(k) = 'CW NOR'
      k = k + 1
      sname(k) = 'srf_wind_dir'
      lname(k) = 'SURFACE WIND DIRECTION'
      units(k) = 'CW NOR'
      k = k + 1
      sname(k) = 'totcld'
      lname(k) = 'TOTAL CLOUD COVER'
      units(k) = '%'
      k = k + 1
      sname(k) = 'mccld'
      lname(k) = 'CONVECTIVE CLOUD COVER'
      units(k) = '%'
      k = k + 1
      sname(k) = 'pcldtop'
      lname(k) = 'CLOUD TOP PRESSURE'
      units(k) = 'MB'
      k = k + 1
      sname(k) = 'pcld_low'
      lname(k) = 'LOW LEVEL CLOUDINESS'
      units(k) = '%'
      k = k + 1
      sname(k) = 'pcld_mid'
      lname(k) = 'MIDDLE LEVEL CLOUDINESS'
      units(k) = '%'
      k = k + 1
      sname(k) = 'pcld_hi'
      lname(k) = 'HIGH LEVEL CLOUDINESS'
      units(k) = '%'
      k = k + 1
      sname(k) = 'net_rad_planet'
      lname(k) = 'NET RAD. OF PLANET'
      units(k) = 'WATTS/M**2'
      k = k + 1
      sname(k) = 'net_rad_z0'
      lname(k) = 'NET RADIATION AT Z0'
      units(k) = 'WATTS/M**2'
      k = k + 1
      sname(k) = 'btemp_wndw'
      lname(k) = 'BRIGHTNESS TEMP THRU WNDW'
      units(k) = 'DEG C'
      k = k + 1
      sname(k) = 'plan_alb'
      lname(k) = 'PLANETARY ALBEDO'
      units(k) = '%'
      k = k + 1
      sname(k) = 'grnd_alb'
      lname(k) = 'GROUND ALBEDO'
      units(k) = '%'
      k = k + 1
      sname(k) = 'vis_alb'
      lname(k) = 'VISUAL ALBEDO'
      units(k) = '%'
      k = k + 1
      sname(k) = 'net_thrml_rad'
      lname(k) = 'NET THRML RADIATION'
      units(k) = 'WATTS/M**2'
      k = k + 1
      sname(k) = 'net_ht_z0'
      lname(k) = 'NET HEAT AT Z0'
      units(k) = 'WATTS/M**2'
      k = k + 1
      sname(k) = 'trop_st_stab'
      lname(k) = 'TROP STATIC STABILITY'
      units(k) = 'DEG K/KM'
      k = k + 1
      sname(k) = 'tot_nt_dse'
      lname(k) = 'TOTAL NT DRY STAT ENR'
      units(k) = '10**14 WT'
      k = k + 1
      sname(k) = 'stand_eddy_nt_dse'
      lname(k) = 'NT DRY STAT ENR BY ST ED'
      units(k) = 'E14 WT'
      k = k + 1
      sname(k) = 'trans_eddy_nt_dse'
      lname(k) = 'NT DRY STAT ENR BY TR ED'
      units(k) = 'E14 WT'
      k = k + 1
      sname(k) = 'z850'
      lname(k) = '850 MB HEIGHT'
      units(k) = 'METERS-1500'
      k = k + 1
      sname(k) = 'z700'
      lname(k) = '700 MB HEIGHT'
      units(k) = 'METERS-3000'
      k = k + 1
      sname(k) = 'z500'
      lname(k) = '500 MB HEIGHT'
      units(k) = 'METERS-5600'
      k = k + 1
      sname(k) = 'z300'
      lname(k) = '300 MB HEIGHT'
      units(k) = 'METERS-9500'
      k = k + 1
      sname(k) = 'z100'
      lname(k) = '100 MB HEIGHT'
      units(k) = 'METERS-16400'
      k = k + 1
      sname(k) = 'z030'
      lname(k) = ' 30 MB HEIGHT'
      units(k) = 'METERS-24000'
      k = k + 1
      sname(k) = 'dztemp_1000-850'
      lname(k) = 'THICKNESS TEMP 1000-850'
      units(k) = 'C'
      k = k + 1
      sname(k) = 'dztemp_850-700'
      lname(k) = 'THICKNESS TEMP 850-700'
      units(k) = 'C'
      k = k + 1
      sname(k) = 'dztemp_700-500'
      lname(k) = 'THICKNESS TEMP 700-500'
      units(k) = 'C'
      k = k + 1
      sname(k) = 'dztemp_500-300'
      lname(k) = 'THICKNESS TEMP 500-300'
      units(k) = 'C'
      k = k + 1
      sname(k) = 'dztemp_300-100'
      lname(k) = 'THICKNESS TEMP 300-100'
      units(k) = 'C'
      k = k + 1
      sname(k) = 'dztemp_100-030'
      lname(k) = 'THICKNESS TEMP 100-30'
      units(k) = 'C'
      k = k + 1
      sname(k) = 'grndh2o'
      lname(k) = 'TOTAL EARTH WATER'
      units(k) = 'KG/M**2'
      k = k + 1
      sname(k) = 'cld_wat_path'
      lname(k) = 'LIQUID WATER PATH '
      units(k) = '.1 KG/M**2'
      k = k + 1
      sname(k) = 'dcfreq'
      lname(k) = 'DEEP CONV CLOUD FREQUENCY'
      units(k) = '%'
      k = k + 1
      sname(k) = 'scfreq'
      lname(k) = 'SHALLOW CONV CLOUD FREQUENCY'
      units(k) = '%'
      k = k + 1
      sname(k) = 'dccld'
      lname(k) = 'DEEP CONVECTIVE CLOUD COVER'
      units(k) = '%'
      k = k + 1
      sname(k) = 'sccld'
      lname(k) = 'SHALLOW CONVECTIVE CLOUD COVER'
      units(k) = '%'
      k = k + 1
      sname(k) = 'dudt_sdrag'
      lname(k) = 'DU/DT BY SDRAG'
      units(k) = '10**-6 M S-2'
      k = k + 1
      sname(k) = ''
      lname(k) = ''
      units(k) = ''
      k = k + 1
      sname(k) = ''
      lname(k) = ''
      units(k) = ''
c imported from ijmap
      k = k + 1
      sname(k) = 'slp'
      lname(k) = 'SEA LEVEL PRESSURE'
      units(k) = 'MB-1000'
      k = k + 1
      sname(k) = 'tsurf'
      lname(k) = 'SURFACE TEMPERATURE'
      units(k) = 'DEGREES C'
      k = k + 1
      sname(k) = 'z850_instant'
      lname(k) = 'INSTANTANEOUS 850 MB HEIGHTS'
      units(k) = 'DEKAMETERS-100'
      k = k + 1
      sname(k) = 'slp_using_t1'
      lname(k) = 'SEA LEVEL PRESSURE (USING T1)'
      units(k) = 'MB-1000'
      k = k + 1
      sname(k) = 'tsurf_using_t1'
      lname(k) = 'SURFACE TEMPERATURE (LAPSE RATE FROM T1)'
      units(k) = 'DEG C'
c imported from sijmap
      k = k + 1
      sname(k) = 'usurf'
      lname(k) = 'U COMP. OF SURFACE AIR WIND'
      units(k) = 'm/s'
      k = k + 1
      sname(k) = 'vsurf'
      lname(k) = 'V COMP. OF SURFACE AIR WIND'
      units(k) = 'm/s'
      k = k + 1
      sname(k) = 'taumag'
      lname(k) = 'MAG OF COMP. MOM. SURF. DRAG'
      units(k) = 'kg/m*s**2'
      k = k + 1
      sname(k) = 'tauus'
      lname(k) = 'U COMP. OF MOM. SRF DRAG'
      units(k) = 'kg/m*S**2'
      k = k + 1
      sname(k) = 'tauvs'
      lname(k) = 'V COMP. OF MOM. SRF DRAG'
      units(k) = 'kg/m*S**2'

c create titles by concatenating long names with units
c no checks whether total length of lname+units exceeds length of title
      do k=1,nij_out
         title(k)=''
         if(lname(k).ne.'')
     &        title(k) = trim(lname(k))//' ('//trim(units(k))//')'
      enddo

      END SUBROUTINE IJ_TITLES

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
C**12  28  FIRST LAYER GROUND TEMPERATURE (K-TF)                  1
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
C**27  44  BRIGHTNESS TEMPERATURE THROUGH WINDOW REGION (K-TF)
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
C**49  50  TOTAL EARTH WATER (KG H2O/M**2)
C****  81  LIQUID WATER PATH (.1 kg/M**2)
C****  84  DEEP CONVECTIVE CLOUD FREQUENCY (%)
C****  85  SHALLOW CONVECTIVE CLOUD FREQUENCY (%)
C****  83  DEEP CONVECTIVE CLOUD COVER (%)
C****  82  SHALLOW CONVECTIVE CLOUD COVER (%)
C****
C**55  98  DU/DT BY SDRAG (10**-6 M S-2)
C****
      USE DAGCOM, only : QCHECK,acc_period,iu_ij
      USE CONSTANT, only :
     &     grav,rgas,sday,twopi,sha,kapa,bygrav,tf
      USE MODEL_COM, only :
     &     im,jm,lm,fim,jeq,
     &     BYIM,DTsrc,FLAND,IDACC,JHOUR,JHOUR0,JDATE,JDATE0,
     &     AMON,AMON0,JYEAR,JYEAR0,NDAY,           !   SKIPSE,
     &     Itime,Itime0,XLABEL,LRUNID,ZATMO
      USE GEOM, only :
     &     AREAG,DXV,DXYP,DXYV
      USE DAGCOM, only :
     &     ajk,kdiag,aij,kaij,aijk,tsfrez,ia_ij,
     *     ijk_v,ijk_dse,ijk_u,ijk_dp,
     &     IJ_PEV,IJ_TRNFP0,IJ_SRNFP0,IJ_SLP,IJ_TS !not a generic subr.
      USE BDIJ, only :
     &     title,nt_ij,legend,achar,bchar,cchar,dchar,echar
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(IM,JM,2) :: ENDE16
      common/ntdse/ende16 ! aka sendeg
      INTEGER, DIMENSION(3) :: MLAT,MGLOBE
      DOUBLE PRECISION, DIMENSION(3) :: FLAT,FNH,FGLOBE,GNUM,GDEN
      REAL*8, DIMENSION(IM,JM,3) :: SMAP
      REAL*8, DIMENSION(JM,3) :: SMAPJ
      CHARACTER XLB*32

      CHARACTER*1 LINE(IM,3),LONGTD(36)
      DATA LONGTD/'+',35*' '/

      INTEGER, DIMENSION(60), PARAMETER :: IND=(/
     &   1, 1, 1, 2, 3,29,    5, 6, 4, 7,32,28,   34,39,36,86,39,36,
     &  19,17,18,41,42,43,   21,22,44,24,26,45,   21,23,31,20, 1, 2,
     &  10,11,12,13,14,15,    9,10,11,12,13,14,   50,81,84,85,83,82,
     &  98,    1,1,1,1,1/)
C**** IA now set from DEFACC

      INTEGER, DIMENSION(3,24) :: ILEG
      DATA ILEG/7,3*1,9,1,   10,10,12, 1,18,11,    3, 5, 3, 4, 2, 2,
     *   1, 1, 6, 1, 1, 1,   13,20,11, 1, 1, 1,   13,13, 3,20,20,18,
     *  12,13,14,15,15,16,   11,11,11,11,11,11,    8, 3, 1, 1, 1, 1,
     *  21, 5*1/
      DOUBLE PRECISION, DIMENSION(60) :: SCALE
      DATA SCALE/1.,3*100.,1.,100.,  3*1.,100.,2*1.,  6*1.,
     *  2*100.,1.,3*100.,  3*1.,3*100.,  2*1.,2.,21*1.,
     *  6*1./
      DOUBLE PRECISION, DIMENSION(60) :: FAC
      DATA FAC/.01,3*.2,1.,.2,  2*10.,.1,.2,10.,.3333333,
     *  2.,.5,2.,10.,2*.1,  2*.2,.02,3*.2,  .05,.1,.3333333,3*.2,
     *  2*.05,2.,2*.1,10.,  .1,.05,.02,.01,.01,.006666667,  6*.3333333,
     *  .05,2.0,.2,.2,.2,.2,  .05,5*0./
      INTEGER, DIMENSION(KAIJ) :: JGRID
      DATA JGRID/19*1,2, 18*1,2*2, 40*1, 5*1,2*2,13*1/
      INTEGER, DIMENSION(3) :: JGRID_THIS_ROW ! for green sheet printout
      DOUBLE PRECISION, DIMENSION(7), PARAMETER :: PMB=(/
     &     1000.,850.,700.,500.,300.,100.,30. /)
      DOUBLE PRECISION, DIMENSION(7), PARAMETER :: GHT=(/
     &     0.,1500.,3000.,5600.,9500.,16400.,24000. /)
      DOUBLE PRECISION, PARAMETER :: P1000=1000.

      INTEGER ::
     &     I,IFRSTP,ILINE,INC,IQ1,IQ2,IQ3,J,K,
     &     KC,KCOLMN,KPAGE,KR,KROW,KT,L,LASTP,M,N,NDEX,NDEX2,KM

      DOUBLE PRECISION ::
     &     A,BYIACC,BYIACN,BYIADA,DAREA,
     &     FDEN,FLATK,FNUM,PLAND,DAYS,ZNDE16,ZS,DPTI,PVTI,
     &     DE4TI,BYDPK,SZNDEG

C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      IF(QCHECK) call open_ij(trim(acc_period)//'.ij'//XLABEL(1:LRUNID))

C**** INITIALIZE CERTAIN QUANTITIES
      call ij_titles
      XLB=acc_period(1:3)//' '//acc_period(4:12)//' '//XLABEL(1:LRUNID)
      INC=1+(JM-1)/24
      ILINE=36*INC
      IQ1=1+IM/(4*INC)
      LONGTD(IQ1)=LONGTD(1)
      IQ2=1+IM/(2*INC)
      LONGTD(IQ2)=LONGTD(1)
      IQ3=1+3*IM/(4*INC)
      LONGTD(IQ3)=LONGTD(1)
      SCALE(7)=SDAY/DTsrc
      SCALE(8)=SDAY/DTSRC
      SCALE(9)=1./DTSRC
      SCALE(11)=SDAY/DTSRC
      SCALE(13)=1.
      SCALE(16)=1E3*BYGRAV
      SCALE(26)=1./DTSRC
      SCALE(32)=1./DTSRC
      SCALE(33)=1.D3*GRAV*(P1000**KAPA)
      SCALE(34)=6.25E-14*BYGRAV
      SCALE(35)=SCALE(34)
      SCALE(36)=SCALE(34)
      DO 70 M=37,42
   70 SCALE(M)=1.*BYGRAV
      DO 80 M=43,48
   80 SCALE(M)=1./(RGAS*LOG(PMB(M-42)/PMB(M-41)))
      SCALE(50)=10.
      SCALE(51)=100.
      SCALE(52)=100.
      SCALE(53)=100.
      SCALE(54)=100.
      SCALE(55)=1.D6/(DTsrc+1.D-20)
      IF (IDACC(12).LT.1) IDACC(12)=1
C****
      DAYS=(Itime-Itime0)/DFLOAT(nday)
CF*** NO PALMER INDEX FOR FINE GRID RUNS
      BYIADA=1./(IDACC(4)+1.D-20)
c     IF (SKIPSE.EQ.1.) GO TO 160
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
      IF (AIJK(I,J,K,IJK_DP).EQ.0.) GO TO 105
      DPTI=DPTI+AIJK(I,J,K,IJK_DP)
      BYDPK=1./(AIJK(I,J,K,IJK_DP)+1.D-20)
      PVTI=PVTI+AIJK(I,J,K,IJK_V)
      DE4TI=DE4TI+AIJK(I,J,K,IJK_DSE)
      ENDE16(I,J,1)=ENDE16(I,J,1)
     *  +(AIJK(I,J,K,3)*AIJK(I,J,K,IJK_V)*BYDPK)
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
  160 CONTINUE
      DO 180 N=1,KAIJ
      IF (JGRID(N).EQ.2) GO TO 180
      DO 170 I=1,IM
      AIJ(I,1,N)=AIJ(1,1,N)
  170 AIJ(I,JM,N)=AIJ(1,JM,N)
  180 CONTINUE
      IFRSTP=1
      LASTP=10
      IF (KDIAG(3).GT.0) LASTP=9-KDIAG(3)
      IF (KDIAG(3).LT.0) IFRSTP=-KDIAG(3)
      IF (KDIAG(3).LT.0) LASTP=IFRSTP
      DO 610 KPAGE=IFRSTP,LASTP
      WRITE (6,901) XLABEL
      WRITE (6,902) JYEAR0,AMON0,JDATE0,JHOUR0,
     *  JYEAR,AMON,JDATE,JHOUR,ITIME,DAYS
      DO 610 KROW=1,2
      IF (KPAGE.EQ.LASTP .AND. KROW.EQ.2) GO TO 610
      KR=2*(KPAGE-1)+KROW
c      WRITE (6,903) (TITLE(K,KR),K=1,3)
      WRITE (6,903) (TITLE(3*(KR-1)+K)(1:32),K=1,3)
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
      BYIACC=1./(IDACC(ia_ij(IND(K)))+1.D-20)
      GO TO (320,340,400,400,440,400, 440,440,460,400,420,460,
     *       380,300,300,475,240,240, 400,400,260,400,400,400,
     *       220,420,460,260,260,260, 460,460,380,420,280,280,
     *       460,460,460,460,460,460, 360,360,360,360,360,360,
     *       380,380,400,400,400,400, 420,215,215,215,215,215),K
C**** Blanks
  215 CONTINUE
      DO I=1,IM
        LINE(I,KCOLMN)=' '
      END DO
      FLATK=0.
      GO TO 500
C**** SUM OF TWO ARRAYS
  220 DO 230 I=1,IM
      A=(AIJ(I,J,IJ_TRNFP0)+AIJ(I,J,IJ_SRNFP0))*SCALE(K)*BYIACC
          SMAP(I,J,KCOLMN)=A
      FLATK=FLATK+A
      N=28.5+A*FAC(K)
      IF (N.LT.1 ) N=1
      IF (N.GT.38) N=38
  230 LINE(I,KCOLMN)=ECHAR(N:N)
      GO TO 500
C**** WIND DIRECTION
  240 IF (J.EQ.1) GO TO 500
      DO 250 I=1,IM
      A=360.*ATAN2(AIJ(I,J,NDEX)+1.D-20,AIJ(I,J,NDEX+1)+1.D-20)/TWOPI
          SMAP(I,J,KCOLMN)=A
      FLATK=FLATK+A
      N=2.5+A*FAC(K)
      IF (N.LT.2) N=N+36
  250 LINE(I,KCOLMN)=ACHAR(N:N)
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
          SMAP(I,J,KCOLMN)=A
  270 LINE(I,KCOLMN)=ACHAR(N:N)
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
  280 CONTINUE
c     IF (SKIPSE.EQ.1.) GO TO 510
      DO 290 I=1,IM
      A=ENDE16(I,J,NDEX)*SCALE(K)/(IDACC(4)+1d-20)
          SMAP(I,J,KCOLMN)=A
      FLATK=FLATK+A
      N=11.5+A*FAC(K)
      IF (N.LT.1) N=1
      IF (N.GT.38) N=38
  290 LINE(I,KCOLMN)=CCHAR(N:N)
      FLAT(KCOLMN)=FLATK
      DAREA=DXYV(J)
      GO TO 505
C**** MAGNITUDE OF TWO PERPENDICULAR COMPONENTS
  300 IF (K.NE.15.AND.J.EQ.1) GO TO 500
      DO 310 I=1,IM
      A=SQRT(AIJ(I,J,NDEX)**2+AIJ(I,J,NDEX+1)**2)*SCALE(K)*BYIACC
          SMAP(I,J,KCOLMN)=A
      FLATK=FLATK+A
      N=2.5+A*FAC(K)
      IF (N.GT.38) N=38
  310 LINE(I,KCOLMN)=ACHAR(N:N)
      GO TO 500
C**** SURFACE TOPOGRAPHY
  320 DO 330 I=1,IM
      ZS=ZATMO(I,J)*BYGRAV
          SMAP(I,J,KCOLMN)=ZS
      FLATK=FLATK+ZS
      N=2.5+.01*ZS
      IF (ZS.LE.0.) N=1
      IF (N.GT.38) N=38
  330 LINE(I,KCOLMN)=ACHAR(N:N)
      GO TO 500
C**** LAND COVERAGE
  340 DO 350 I=1,IM
      PLAND=FLAND(I,J)*100.
          SMAP(I,J,KCOLMN)=PLAND
      FLATK=FLATK+PLAND
      N=2.5+PLAND*.2
      IF (PLAND.LE.0.) N=1
      IF (PLAND.GE.100.) N=23
  350 LINE(I,KCOLMN)=BCHAR(N:N)
      GO TO 500
C**** THICKNESS TEMPERATURES
  360 DO 370 I=1,IM
      A=((AIJ(I,J,NDEX+1)-AIJ(I,J,NDEX))*BYIACC
     *  +(GHT(NDEX-7)-GHT(NDEX-8))*GRAV)*SCALE(K)-TF
          SMAP(I,J,KCOLMN)=A
      FLATK=FLATK+A
      N=28.5+A*FAC(K)
      IF (N.LT.1) N=1
      IF (N.GT.38) N=38
  370 LINE(I,KCOLMN)=ECHAR(N:N)
      GO TO 500
C**** POSITIVE QUANTITIES UNIFORMLY SCALED
  380 DO 390 I=1,IM
      A=AIJ(I,J,NDEX)*SCALE(K)*BYIACC
          SMAP(I,J,KCOLMN)=A
      FLATK=FLATK+A
      N=2.5+A*FAC(K)
      IF (A.EQ.0.) N=1
      IF (N.GT.38) N=38
  390 LINE(I,KCOLMN)=ACHAR(N:N)
      GO TO 500
C**** PERCENTAGES
  400 DO 410 I=1,IM
      A=AIJ(I,J,NDEX)*SCALE(K)*BYIACC
          SMAP(I,J,KCOLMN)=A
      FLATK=FLATK+A
      N=2.5+A*FAC(K)
      IF (A.LE.0.) N=1
      IF (A*FAC(K).GE.20.) N=23
  410 LINE(I,KCOLMN)=BCHAR(N:N)
      GO TO 500
C**** SIGNED QUANTITIES UNIFORMLY SCALED (LETTERS +, NUMBERS -)
  420 DO 430 I=1,IM
      A=AIJ(I,J,NDEX)*SCALE(K)*BYIACC
          SMAP(I,J,KCOLMN)=A
      FLATK=FLATK+A
      N=11.5+A*FAC(K)
      IF (N.LT.1) N=1
      IF (N.GT.38) N=38
  430 LINE(I,KCOLMN)=CCHAR(N:N)
      IF (K.EQ.34) FLATK=FLATK*FIM
      GO TO 500
C**** PRECIPITATION AND EVAPORATION
  440 DO 450 I=1,IM
      A=AIJ(I,J,NDEX)*SCALE(K)*BYIACC
          SMAP(I,J,KCOLMN)=A
      FLATK=FLATK+A
      N=1
      IF (A.LE.0.) GO TO 450
      N=2.5+A*FAC(K)
      IF (N.GT.28) N=(N+263)/10
      IF (N.GT.35) N=(N+180)/6
      IF (N.GT.37) N=37
  450 LINE(I,KCOLMN)=DCHAR(N:N)
      GO TO 500
C**** SIGNED QUANTITIES UNIFORMLY SCALED (NUMBERS +, LETTERS -)
  460 DO 470 I=1,IM
      A=AIJ(I,J,NDEX)*SCALE(K)*BYIACC
          SMAP(I,J,KCOLMN)=A
      FLATK=FLATK+A
      N=28.5+A*FAC(K)
      IF (N.LT.1 ) N=1
      IF (N.GT.38) N=38
  470 LINE(I,KCOLMN)=ECHAR(N:N)
      GO TO 500
C**** ABSOLUTE VALUE OF QUANTITIES UNIFORMLY SCALED
  475 DO 477 I=1,IM
      A=ABS(AIJ(I,J,NDEX)*SCALE(K)*BYIACC)
          SMAP(I,J,KCOLMN)=A
      FLATK=FLATK+A
      N=2.5+A*FAC(K)
      IF(A.EQ.0.) N=1
      IF (N.GT.38) N=38
  477 LINE(I,KCOLMN)=ACHAR(N:N)
      GO TO 500
C**** POSITIVE QUANTITIES NON-UNIFORMLY SCALED
  480 DO 490 I=1,IM
      A=AIJ(I,J,NDEX)*SCALE(K)*BYIACC
          SMAP(I,J,KCOLMN)=A
      FLATK=FLATK+A
      N=2.5+A*FAC(K)
      IF (N.GE.13) N=(N+123)/10
      IF (N.GT.38) N=38
  490 LINE(I,KCOLMN)=ACHAR(N:N)
      GO TO 500
C**** LENGTH OF GROWING SEASON
  491 DO 492 I=1,IM
      A=TSFREZ(I,J,2)-TSFREZ(I,J,1)
      IF (A.LT.0.) A=A+365.
          SMAP(I,J,KCOLMN)=A
      FLATK=FLATK+A
      N=2.5+A*FAC(K)
      IF (A.LE.0.) N=1
      IF (N.GT.38) N=38
  492 LINE(I,KCOLMN)=ACHAR(N:N)
      GO TO 500
C**** PALMER DROUGHT INDEX
  493 DO 494 I=1,IM
      A=0.
      FLATK=FLATK+A
      N=11.5+A*FAC(K)
      IF (N.LT.1 ) N=1
      IF (N.GT.38) N=38
  494 LINE(I,KCOLMN)=CCHAR(N:N)
  500 FLAT(KCOLMN)=FLATK*BYIM
      DAREA=DXYP(J)
      IF (JGRID(NDEX).EQ.2) DAREA=DXYV(J)
      JGRID_THIS_ROW(KCOLMN) = JGRID(NDEX)
      IF (J.GE.JEQ) FNH(KCOLMN)=FNH(KCOLMN)+FLAT(KCOLMN)*DAREA
  505 FGLOBE(KCOLMN)=FGLOBE(KCOLMN)+FLAT(KCOLMN)*DAREA
CB       DIJMAP(IM+1,J,K)=FLAT(KCOLMN)
      MLAT(KCOLMN)=NINT(FLAT(KCOLMN))
  510 CONTINUE
          DO KCOLMN=1,3
             SMAPJ(J,KCOLMN) = FLAT(KCOLMN)
          END DO
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
  600 WRITE (6,909) (LEGEND(ILEG(KCOLMN,KR))(1:40),KCOLMN=1,2),
     *  LEGEND(ILEG(3,KR))(1:36)
      IF(QCHECK) THEN
        DO KC=1,3
          nt_ij = 3*(KR-1)+KC
          IF(TITLE(3*(KR-1)+KC).ne.' ')
     &          CALL POUT_IJ(TITLE(3*(KR-1)+KC)//XLB,
     *          SMAP(1,1,KC),SMAPJ(1,KC),FGLOBE(KC),JGRID_THIS_ROW(KC))
c     *        WRITE(iu_ij) TITLE(KC,KR),XLB,((SMAP(I,J,KC),I=1,IM),
c     *               J=1,JM),(SMAPJ(J,KC),J=1,JM),SNGL(FGLOBE(KC))
            END DO
         END IF
  610 CONTINUE
  690 CONTINUE
C****
C**** PRODUCE FULL PAGE I,J MAPS
C****
      nt_ij = 58
      CALL IJMAP (TITLE(nt_ij),AIJ(1,1,IJ_SLP),BYIADA)
      BYIACN=1./(IDACC(3)+1.D-20)
      nt_ij = 59
      CALL IJMAP (TITLE(nt_ij),AIJ(1,1,IJ_TS),BYIACN)
C     CALL IJMAP (TITLE(61),AIJ(1,1,IJ_SLP1),BYIADA)
C     CALL IJMAP (TITLE(62),AIJ(1,1,IJ_TS1),BYIADA)
      IF(QCHECK) THEN
         CALL SIJMAP
         call close_ij
         CALL IJKMAP
      ENDIF
      RETURN
C****
  901 FORMAT ('1',A)
  902 FORMAT ('0',15X,'From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,'  Model-Time:',I9,5X,
     *  'Dif:',F7.2,' Days')
  903 FORMAT ('0',6X,A32,13X,A32,13X,A32)
  906 FORMAT ('+',6X,36A1,9X,36A1,9X,36A1)
  909 FORMAT (7X,A40,5X,A40,5X,A36)
  910 FORMAT (1X,F5.1,1X,36A1,F8.1,1X,36A1,F8.1,1X,36A1)
  911 FORMAT (1X,F5.1,1X,36A1,F8.1,1X,36A1,I8,1X,36A1)
  914 FORMAT (1X,I5,1X,36A1,F8.1,1X,36A1,F8.1,1X,36A1)
  916 FORMAT (1X,I5,1X,36A1,I8,1X,36A1,F8.1,1X,36A1)
  917 FORMAT (1X,I5,1X,36A1,I8,1X,36A1,I8,1X,36A1)
      END SUBROUTINE DIAGIJ

      SUBROUTINE IJMAP (TITLE,ARRAY,BYIACC)
      USE DAGCOM, only : QCHECK,acc_period,iu_ij
      USE MODEL_COM, only :
     &     im,jm,BYIM,FLAND,NDAY,JHOUR,JHOUR0,JDATE,JDATE0,AMON,AMON0,
     &     JYEAR,JYEAR0,Itime,Itime0,XLABEL,lrunid
      USE GEOM, only :
     &     LAT_DG,LON_DG
      IMPLICIT NONE

      DOUBLE PRECISION :: BYIACC
      DOUBLE PRECISION, DIMENSION(IM,JM) :: ARRAY
      CHARACTER*48 TITLE

      CHARACTER*1 IDX(12),BLANK,xlb*32
      CHARACTER(LEN=3), DIMENSION(IM) :: LINE
      CHARACTER(LEN=9) :: AVG
      DATA IDX/'0','1','2','3','4','5','6','7','8','9','-','*'/
      DATA BLANK/' '/

      REAL*8, DIMENSION(IM,JM) :: SMAP
      REAL*8, DIMENSION(JM) :: SMAPJ
      DOUBLE PRECISION :: A,FLAT,DAYS
      INTEGER :: I,IA,INC,J,JA,JX,K,LD

C****
C**** INITIALIZE CERTAIN QUANTITIES
C****
c      BYIM=1./IM
      XLB=acc_period(1:3)//' '//acc_period(4:12)//' '//XLABEL(1:LRUNID)
      INC=(IM+35)/36
      LD=360/IM
      DO 40 I=1,IM
      WRITE(LINE(I),'(I3)') I
   40 CONTINUE
      AVG='     MEAN'
      DAYS=(Itime-Itime0)/DFLOAT(nday)
      WRITE(6,901)XLABEL
      WRITE (6,902) JYEAR0,AMON0,JDATE0,JHOUR0,
     *  JYEAR,AMON,JDATE,JHOUR,ITIME,DAYS
      WRITE(6,900) TITLE
      WRITE (6,910) (LINE(I),I=1,IM,INC),AVG
      WRITE(6,940)
      WRITE(6,940)
C**** OUTSIDE J LOOP
      DO 300 J=JM,1,-1
      FLAT=0.
      DO 250 I=1,IM
      A=ARRAY(I,J)*BYIACC
           SMAP(I,J)=A
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
         SMAPJ(J)=FLAT
      WRITE(AVG,'(F9.2)') FLAT
      WRITE (6,920) NINT(LAT_DG(J,1)),J,(LINE(I),I=1,IM,INC),AVG
      DO 260 I=1,IM
      IF (FLAND(I,J).GE..5) GO TO 260
      DO 255 K=1,3
  255 LINE(I)(K:K)=BLANK
  260 CONTINUE
C     WRITE (6,925) (LINE(I),I=1,IM,INC)
      WRITE (6,925) (LINE(I),I=1,IM,INC)
      WRITE (6,925) (LINE(I),I=1,IM,INC)
  300 IF (JM.LE.24) WRITE (6,940)
      WRITE (6,930) (LON_DG(I,1),I=1,IM,INC*2)
         IF(QCHECK) CALL POUT_IJ(TITLE//XLB,SMAP,SMAPJ,-1.D30,1)
c                   WRITE(iu_ij) TITLE,XLB,SMAP,SMAPJ,-1.E30
      RETURN
C****
  900 FORMAT('0',45X,A48)
  901 FORMAT ('1',A)
  902 FORMAT ('0',15X,'From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,'  Model-Time:',I9,5X,
     *  'Dif:',F7.2,' Days')
  910 FORMAT('0LAT  J/I  ',36A3,A9)
  920 FORMAT(2I4,3X,36A3,A9)
  925 FORMAT('+',10X,36A3)
  930 FORMAT('0  LONG',2X,19F6.1)
  940 FORMAT(' ')
      END SUBROUTINE IJMAP

      SUBROUTINE DIAG9P
C****
C**** THIS ENTRY PRODUCES TABLES OF CONSERVATION QUANTITIES
C****
c      USE PRTCOM, only :
      USE CONSTANT, only :
     &     twopi
      USE MODEL_COM, only :
     &     jm,fim,idacc,jhour,jhour0,jdate,jdate0,amon,amon0,
     &     jyear,jyear0,nday,jeq,itime,itime0,xlabel
      USE GEOM, only :
     &     areag,dlon,dxyp,dxyv,LAT_DG
      USE DAGCOM, only :
     &     consrv,kcon,scale_con,title_con,nsum_con,ia_con,kcmx
      IMPLICIT NONE

      INTEGER, DIMENSION(JM) :: MAREA
      DOUBLE PRECISION, DIMENSION(KCON) :: FGLOB
      DOUBLE PRECISION, DIMENSION(2,KCON) :: FHEM
      INTEGER, DIMENSION(JM,KCON) :: MLAT
      DOUBLE PRECISION, DIMENSION(JM+3,KCON) :: CNSLAT
      CHARACTER*4, PARAMETER :: HEMIS(2) = (/' SH ',' NH '/),
     *     DASH = ('----')
      DOUBLE PRECISION, PARAMETER :: XWON = TWOPI/(DLON*FIM)

      INTEGER :: inc,j,jhemi,jnh,jp1,jpm,jsh,jv1,jvm,jx,n
      DOUBLE PRECISION :: aglob,ahem,feq,fnh,fsh,days
C**** CALCULATE SCALING FACTORS
      IF (IDACC(12).LT.1) IDACC(12)=1
C**** CALCULATE SUMMED QUANTITIES
C**** LOOP BACKWARDS SO THAT INITIALISATION IS DONE BEFORE SUMMATION!
      DO J=1,JM
        DO N=KCMX,1,-1
          IF (NSUM_CON(N).eq.0) THEN
            CONSRV(J,N)=0.
          ELSEIF (NSUM_CON(N).gt.0) THEN
            CONSRV(J,NSUM_CON(N))=CONSRV(J,NSUM_CON(N))+CONSRV(J,N)
     *           *SCALE_CON(N)*IDACC(12)/(IDACC(IA_CON(N))+1d-20)
          END IF
        END DO
      END DO
C**** CALCULATE FINAL ANGULAR MOMENTUM + KINETIC ENERGY ON VELOCITY GRID
      DO N=1,23
        FEQ=CONSRV(JEQ,N)*SCALE_CON(N)/(IDACC(IA_CON(N))+1d-20)
        FGLOB(N)=FEQ
        FHEM(1,N)=.5*FEQ
        FHEM(2,N)=.5*FEQ
        CNSLAT(JEQ,N)=FEQ/(FIM*DXYV(JEQ))
        DO JSH=2,JEQ-1
          JNH=2+JM-JSH
          FSH=CONSRV(JSH,N)*SCALE_CON(N)/(IDACC(IA_CON(N))+1d-20)
          FNH=CONSRV(JNH,N)*SCALE_CON(N)/(IDACC(IA_CON(N))+1d-20)
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
C**** CALCULATE ALL OTHER CONSERVED QUANTITIES ON TRACER GRID
      DO N=24,KCMX
        FGLOB(N)=0.
        FHEM(1,N)=0.
        FHEM(2,N)=0.
        DO JSH=1,JEQ-1
          JNH=1+JM-JSH
          FSH=CONSRV(JSH,N)*SCALE_CON(N)/(IDACC(IA_CON(N))+1d-20)
          FNH=CONSRV(JNH,N)*SCALE_CON(N)/(IDACC(IA_CON(N))+1d-20)
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
      DAYS=(Itime-Itime0)/DFLOAT(nday)
      DO N=1,KCMX
        DO J=1,JM
          MLAT(J,N)=NINT(CNSLAT(J,N))
        END DO
        CNSLAT(JM+1,N)=FHEM(1,N)
        CNSLAT(JM+2,N)=FHEM(2,N)
        CNSLAT(JM+3,N)=FGLOB(N)
      END DO
      DO JHEMI=2,1,-1
        WRITE (6,901) XLABEL
        WRITE (6,902) JYEAR0,AMON0,JDATE0,JHOUR0,
     *       JYEAR,AMON,JDATE,JHOUR,ITIME,DAYS
        JP1=1+(JHEMI-1)*(JEQ-1)
        JPM=JHEMI*(JEQ-1)
        JV1=2+(JHEMI-1)*(JEQ-2)
        JVM=JEQ+(JHEMI-1)*(JEQ-2)
C**** PRODUCE TABLES FOR ANGULAR MOMENTUM AND KINETIC ENERGY
        WRITE (6,903) (DASH,J=JV1,JVM,INC)
        WRITE (6,904) HEMIS(JHEMI),(NINT(LAT_DG(JX,2)),JX=JVM,JV1,-INC)
        WRITE (6,903) (DASH,J=JV1,JVM,INC)
        DO N=1,23
          WRITE (6,905) TITLE_CON(N),FGLOB(N),FHEM(JHEMI,N),
     *         (MLAT(JX,N),JX=JVM,JV1,-INC)
        END DO
        DO J=JV1,JVM
          MAREA(J)=1.D-10*XWON*FIM*DXYV(J)+.5
        END DO
        WRITE (6,906) AGLOB,AHEM,(MAREA(JX),JX=JVM,JV1,-INC)
C**** PRODUCE TABLES FOR OTHER CONSERVED QUANTITIES
        WRITE (6,907)
        WRITE (6,903) (DASH,J=JP1,JPM,INC)
        WRITE (6,904) HEMIS(JHEMI),(NINT(LAT_DG(JX,1)),JX=JPM,JP1,-INC)
        WRITE (6,903) (DASH,J=JP1,JPM,INC)
        DO N=24,KCMX
          WRITE (6,905) TITLE_CON(N),FGLOB(N),FHEM(JHEMI,N),
     *         (MLAT(JX,N),JX=JPM,JP1,-INC)
        END DO
        DO J=JP1,JPM
          MAREA(J)=1.D-10*XWON*FIM*DXYP(J)+.5
        END DO
        WRITE (6,906) AGLOB,AHEM,(MAREA(JX),JX=JPM,JP1,-INC)
      END DO
      RETURN
C****
  901 FORMAT ('1',A)
  902 FORMAT ('0Conservation Quantities       From:',
     *  I6,A6,I2,',  Hr',I3,  6X,  'To:',I6,A6,I2,', Hr',I3,
     *  '  Model-Time:',I9,5X,'Dif:',F7.2,' Days')
  903 FORMAT (1X,25('--'),13(A4,'--'))
  904 FORMAT (35X,'GLOBAL',A7,2X,13I6)
  905 FORMAT (A32,2F9.2,1X,13I6)
  906 FORMAT ('0AREA (10**10 M**2)',F22.1,F9.1,1X,13I6)
  907 FORMAT ('0')
      END SUBROUTINE DIAG9P

      SUBROUTINE DIAG5P
C**** THIS ENTRY PRINTS THE SPECTRAL ANALYSIS TABLES
C****
c      USE PRTCOM, only :
      USE CONSTANT, only :
     &     grav,rgas,twopi
      USE MODEL_COM, only :
     &     im,jm,lm,fim,
     &     DT,IDACC,JHOUR,JHOUR0,JDATE,JDATE0,
     &     AMON,AMON0,JYEAR,JYEAR0,LS1,JEQ,XLABEL,istrat
!!!  &     ,SKIPSE
      USE GEOM, only :
     &     DLON,DXYV
      USE DAGCOM, only :
     &     speca,atpe,ajk,aijk,kspeca,ktpe,nhemi,nspher,ijk_u,klayer
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
c      CHARACTER*16 :: SPHERE(2) = (/
c     *     'STRATOSPHERE    ','TROPOSPHERE     '/)
      CHARACTER*16 :: SPHERE(4)=
     *     (/'TROPOSPHERE     ','LOW STRATOSPHERE',
     *       'MID STRATOSPHERE','UPP STRATOSPHERE'/)

      INTEGER ::
     &     I,IUNITJ,IUNITW,J,J45N,
     &     K,KPAGE,KROW,KSPHER,L,
     &     M,MAPE,MTPE,N,NM,NM1

      DOUBLE PRECISION ::
     &     FACTOR,FNM

      NM=1+IM/2
      IF (IDACC(12).LT.1) IDACC(12)=1
c     IF (SKIPSE.GE.1.) GO TO 600
C      JEQ=1+JM/2
      J45N=2.+.75*(JM-1.)
C****
C**** STANDING KINETIC ENERGY
C****
      DO 710 K=1,NSPHER
      DO 710 N=1,NM
  710 SPECA(N,1,K)=0.
      DO 770 L=1,LM
        KSPHER=KLAYER(L)
c      KSPHER=2
c      IF (L.GE.LS1) KSPHER=1
      DO 770 J=2,JM
      IF (AJK(J,L,2).LE.1.D-20) GO TO 770
      FACTOR=FIM*DXYV(J)/AJK(J,L,2)
      DO 769 K=IZERO,LM,LM
      DO 720 I=1,IM
  720 X(I)=AIJK(I,J,L+K,IJK_U)
      CALL FFTE (X,X)
      IF (J.EQ.JEQ) GO TO 750
      DO 730 N=1,NM
  730 SPECA(N,1,KSPHER)=SPECA(N,1,KSPHER)+X(N)*FACTOR
      IF (J.NE.J45N) GO TO 769
      DO 740 N=1,NM
  740 SPECA(N,1,KSPHER+2)=SPECA(N,1,KSPHER+2)+X(N)*FACTOR
      GO TO 769
  750 DO 760 N=1,NM
      SPECA(N,1,KSPHER+2)=SPECA(N,1,KSPHER+2)+X(N)*FACTOR
      SPECA(N,1,KSPHER)=SPECA(N,1,KSPHER)+.5*X(N)*FACTOR
  760 SPECA(N,1,KSPHER+1)=SPECA(N,1,KSPHER+1)+.5*X(N)*FACTOR
      IF (K.EQ.LM) KSPHER=KSPHER+1
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
      SCALE(16)=100.D-12/(GRAV*DT*(.5*IDACC(9)+1.D-20))
      SCALE(17)=SCALE(16)*RGAS
      SCALE(18)=100.D-17/(GRAV*IDACC(4)+1.D-20)
      DO 605 K=1,KSPECA
  605 SCALE(K)=(TWOPI/(DLON*FIM))*SCALE(K)
      IUNITJ=17
      IUNITW=12
      DO 690 KPAGE=1,4   ! one for each loc. SH/NH/EQ/45N
C**** WRITE HEADINGS
      WRITE (6,901) XLABEL
      WRITE (6,902) JYEAR0,AMON0,JDATE0,JHOUR0,JYEAR,AMON,JDATE,JHOUR,
     *  IUNITJ,IUNITW
      DO 670 KROW=1,2+ISTRAT
      IF (JM.GE.25.AND.KROW.EQ.2) WRITE (6,901)
      WRITE (6,903) LATITD(KPAGE),SPHERE(KROW)
c      KSPHER=2*(KPAGE-1)+KROW
      KSPHER=4*(KROW-1)+KPAGE
C**** WRITE KINETIC AND AVAILABLE POTENTIAL ENERGY BY WAVE NUMBER
      DO 610 M=1,KSPECA
      F0(M)=SPECA(1,M,KSPHER)*SCALE(M)
      MN(M)=NINT(F0(M))
  610 FNSUM(M)=0.
      WRITE (6,904) MN
      DO 630 N=2,NM
c      KSPHER=2*(KPAGE-1)+KROW
      KSPHER=4*(KROW-1)+KPAGE
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
  902 FORMAT ('0**  Spectral Analysis **      From:',
     *  I6,A6,I2,',  Hr',I3,  6X,  'To:',I6,A6,I2,', Hr',I3,
     *  '       UNITS 10**',I2,' JOULES AND 10**',I2,' WATTS')
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
     *  '        ','        ','        ','        ','        ',
     *  '        ','        ','0SW ON G',' LW AT G',' SNSB HT',
     *  ' LAT HT ',' HEAT Z0','0UG*10  ',' VG*10  ',' WG*10  ',
     *  ' US*10  ',' VS*10  ',' WS*10  ',' ALPHA0 ','0RIS1*E2',
     *  ' RIGS*E2',' CM*E4  ',' CH*E4  ',' CQ*E4  ',' EDS1*10',
     *  '0DBL    ',' DC FREQ',' LDC*10 ','0PRC*10 ',' EVP*10 ',
     *  ' DEEP MC',' SHLW MC',' CLD 7  ',' CLD 6  ',' CLD 5  ',
     *  ' CLD 4  ',' CLD 3  ',' CLD 2  ',' CLD 1  ',' W TO-5 ',
     *  ' C COVER',' SS P*10',' MC P*10'/
      END BLOCK DATA BDDLY

      SUBROUTINE DIAG6
C****
C**** THIS SUBROUTINE PRINTS THE DIURNAL CYCLE OF SOME QUANTITIES
C****
c      USE PRTCOM, only :
      USE CONSTANT, only :
     &     grav,sha,sday,twopi
      USE MODEL_COM, only :
     &     dtsrc,idacc,
     &     IJD6,JDATE,JDATE0,AMON,AMON0,JYEAR,JYEAR0,NAMD6,
     &     NISURF,XLABEL
      USE DAGCOM, only :
     &     adaily,kdiag,ndlyvar,hr_in_day
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(NDLYVAR) :: SCALE
      DATA SCALE/1.,2*100.,2*1.,  5*1.,  3*1.,2*1.D5,  5*1.D5,
     *  5*100.,  2*100.,3*1.,  2*1.,3*10.,  3*10.,1.,100.,
     *  100.,3*1.D4,1*10.,  1.,100.,10.,2*1.,
     *  2*100., 7*100., 1.,100., 2*10./
      DOUBLE PRECISION, DIMENSION(HR_IN_DAY+1) :: XHOUR
      INTEGER, DIMENSION(HR_IN_DAY+1) :: MHOUR

      CHARACTER*8 TITLE
      COMMON/D6COM/TITLE(NDLYVAR)

      DOUBLE PRECISION :: AVE,AVED,AVEN,BYIDAC,DTSURF
      INTEGER :: I,IH,IREGF,IREGL,IS,K,KQ,KR,NDAYS

C****
      NDAYS=IDACC(9)/2
      IF (NDAYS.LE.0) RETURN
      DTSURF=DTsrc/NISURF
      BYIDAC=1./NDAYS
      SCALE(5)=100.*SHA/(GRAV*DTsrc)
      SCALE(28)=1./DTSURF
      SCALE(29)=1./DTSURF
      SCALE(30)=1./DTSURF
      SCALE(31)=1./DTSURF
      SCALE(32)=1./DTSURF
      SCALE(39)=360./TWOPI
      SCALE(49)=100.*100.*SDAY/(DTsrc*GRAV)
      SCALE(62)=100.*100.*SDAY/(DTsrc*GRAV)
      SCALE(63)=SCALE(62)
      SCALE(50)=100.*SDAY/DTSURF
C****
      IREGF=1
      IREGL=4
      IF (KDIAG(6).GT.0) IREGL=4-KDIAG(6)
      IF (KDIAG(6).LT.0.AND.KDIAG(6).GT.-5) IREGF=-KDIAG(6)
      IF (KDIAG(6).LT.0) IREGL=IREGF
      DO 500 KR=IREGF,IREGL
      WRITE (6,901) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
      WRITE (6,903) NAMD6(KR),IJD6(1,KR),IJD6(2,KR),(I,I=1,HR_IN_DAY)
      DO 500 KQ=1,NDLYVAR
      IF (KQ.EQ.48) GO TO 200
      IF (len_trim(TITLE(KQ)).eq.0) GO TO 500
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
  901 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
  903 FORMAT ('0',A4,I2,',',I2,' ',I2,23I5,'  AVE')
  904 FORMAT (A8,25I5)
      END SUBROUTINE DIAG6

      SUBROUTINE DIAG4
C**** THIS ENTRY PRODUCES A TIME HISTORY TABLE OF ENERGIES
C****
c      USE PRTCOM, only :
      USE CONSTANT, only :
     &     grav,rgas,twopi,bygrav
      USE MODEL_COM, only :
     &     im,jm,lm,fim,
     &     IDACC,JHOUR,JHOUR0,JDATE,JDATE0,AMON,AMON0,
     &     JYEAR,JYEAR0,NDA4,NDAY,Itime0,XLABEL,istrat
      USE GEOM, only :
     &     DLON
      USE DAGCOM, only :
     &     energy,nehist,hist_days
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(NEHIST) :: SUM,SCALE
      INTEGER, DIMENSION(NEHIST) :: IK
      DOUBLE PRECISION, DIMENSION(NEHIST,HIST_DAYS+1) :: EHIST

      INTEGER ::
     &     I,IDACC5,ItimeX,IDAYX,IDAYXM,K,K0
      DOUBLE PRECISION ::
     &     TOFDYX, FAC

      IDACC5=IDACC(5)
      IF (IDACC5.LE.0) RETURN
      IF (IDACC(12).LT.1) IDACC(12)=1
      SCALE(1)=100.D-18*BYGRAV
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
      DO K0=1,MIN(1+ISTRAT,2)
        WRITE (6,901) XLABEL
        IF (K0.eq.1) THEN
          FAC = 1.
          WRITE (6,902) JYEAR0,AMON0,JDATE0,JHOUR0,JYEAR,AMON,JDATE
     *         ,JHOUR
          WRITE (6,903)
        ELSE
          FAC = 10.
          WRITE (6,906) JYEAR0,AMON0,JDATE0,JHOUR0,JYEAR,AMON,JDATE
     *         ,JHOUR
          WRITE (6,907)
        END IF
        SUM(:)=0.
        DO I=1,IDACC5
          ItimeX=Itime0+I*NDA4-1
          IDAYX=1+ItimeX/NDAY
          IDAYXM=MOD(IDAYX,100000)
          TOFDYX=MOD(ItimeX,NDAY)*24./NDAY
          DO K=1,MIN((K0-1)*20,NEHIST)
            EHIST(K,I)=ENERGY(K,I)*SCALE(K)*FAC
            IK(K)=EHIST(K,I)+.5
            SUM(K)=SUM(K)+ENERGY(K,I)
          END DO 
          WRITE (6,904) IDAYXM,TOFDYX,IK
        END DO
        DO K=1,NEHIST
          EHIST(K,101)=SUM(K)*SCALE(K)*FAC/IDACC5
          IK(K)=EHIST(K,101)+.5
        END DO
        WRITE (6,905) IK
        IF (K0.eq.1) CALL KEYD4 (IK)
      END DO
      RETURN
C****
  901 FORMAT ('1',A)
  902 FORMAT ('0** ENERGY HISTORY **      From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,
     *  '    UNITS OF 10**18 JOULES')
  903 FORMAT ('0',15X,21('-'),' TROPOSPHERE ',22('-'),5X,21('-'),
     *  '  LOW STRAT. * 10 ',19('-')/8X,2(11X,'ZKE',8X,'EKE',7X,
     *     'SEKE',9X,
     * 'ZPE',10X,'EPE')/3X,'DAY  HOUR     SH   NH    SH   NH     1    2
     *    SH    NH     SH    NH      SH   NH    SH   NH    SH   NH     S
     *H    NH     SH    NH'/1X,132('='))
  904 FORMAT (I6,F6.1,1X,3(I6,I5),2(I7,I6),2X,3(I6,I5),2(I7,I6))
  905 FORMAT (1X,132('=')/8X,'MEAN ',3(I6,I5),2(I7,I6),2X,3(I6,I5),
     *  2(I7,I6))
  906 FORMAT ('0** ENERGY HISTORY **   DAY',I6,', HR',I3,' (',I2,A5,I5,
     *  ')    TO    DAY',I5,', HR',I3,' (',I2,A5,I5,
     *  ')    UNITS OF 10**17 JOULES')
  907 FORMAT ('0',15X,19('-'),' MID STRATOSPHERE ',19('-'),5X,18('-'),
     *  ' HIGH STRAT. * 10  ',19('-')/8X,2(11X,'ZKE',8X,'EKE',7X,
     *   'NHKE',9X,
     *  'ZPE',10X,'EPE')/3X,'DAY  HOUR     SH   NH    SH   NH    1    2
     *    SH    NH     SH    NH      SH   NH    SH   NH    1    2      S
     *H    NH     SH    NH'/1X,132('='))
  920 FORMAT (1X)
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
c      USE PRTCOM, only :
      USE CONSTANT, only :
     &     twopi
      USE MODEL_COM, only :
     &     jm,lm,jeq, JHOUR,JHOUR0,
     &     JDATE,JDATE0,JMON,JMON0,AMON,AMON0,JYEAR,JYEAR0,
     &     Itime,ItimeI,Itime0,PTOP,SIG,XLABEL,
     &     PSFMPT,AMONTH,nday
      USE GEOM, only :
     &     DLAT,DXYP,LAT_DG
      USE DAGCOM, only :
     &     keyct,keynr,nehist,nkeynr
      USE PARAM
      IMPLICIT NONE
      SAVE
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
     &     I,I35,I70,IEND,ISIGN,
     &     J,J60,JMAX,JEQP1,JNDEX,JSTART,
     &     K,KEYMAX,KNDEX,
     &     L,LL,LMAX,LNLM,LNM,LSLM,LSM,
     &     N,NT,NTDIF

      DOUBLE PRECISION ::
     &     A,BIG,CEPT,CHECK,
     &     FGLOB,FNH,GSUM,HN,HS,PISG,PISN,
     &     SAVE,SUMFAC,DAYS,TEQ,TNOR,TSOU,
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
      DO 510 I=JEQ+1,JM
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
  394 KEYNR(38,KEYCT)=NINT(LAT_DG(JNDEX,2))
  396 KEYNR(NT+9,KEYCT)=NINT(BIG)
      RETURN
  398 KEYNR(42,KEYCT)=NINT(LAT_DG(JNDEX,2))
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
      DAYS=(Itime-Itime0)/DFLOAT(nday)
      KEYNR(1,KEYCT)=JMON0
      IF (Itime.eq.ItimeI+1) KEYNR(1,KEYCT)=0
      IF (KEYCT.GE.2.AND.KEYNR(1,KEYCT-1).EQ.JMON0) KEYCT=KEYCT-1
      WRITE(6,901) XLABEL
      WRITE(6,910) JYEAR0,AMON0,JDATE0,JHOUR0,
     *  JYEAR,AMON,JDATE,JHOUR,ITIME,DAYS
      WRITE(6,902)
      DO 810 K=1,KEYCT
      IF (KEYNR(1,K).EQ.1) WRITE (6,905)
  810 WRITE(6,905) AMONTH(KEYNR(NDEX(1),K)),(KEYNR(NDEX(I),K),I=2,42)
      WRITE (6,915)
CB       DO 815 I=1,42
CB815    FKEYDS(I)=KEYNR(I,KEYCT)
      KEYCT=KEYCT+1
      KEYMAX=49
      IF (KEYNR(1,1).NE.0) KEYMAX=48
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
  910 FORMAT ('0',15X,'From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,
     *  '  Model-Time:',I9,5X,'Dif:',F7.2,' Days')
  915 FORMAT('0')
      END SUBROUTINE DIAGKS

      SUBROUTINE SIJMAP
C****
C**** CALCULATE EXTRA QUANTITIES FOR SAVING
C****
      USE MODEL_COM, only :
     &     im,jm,IDACC,BYIM,XLABEL
      USE DAGCOM, only :
     &     aij,acc_period,ia_ij,kaij,
     &     IJ_US,IJ_VS,IJ_TAUS,IJ_TAUUS,IJ_TAUVS !not a generic subr.
      USE BDIJ, only :
     &     title,nt_ij
      IMPLICIT NONE

      CHARACTER XLB*32
      REAL*8 SMAP(IM,JM),SMAPJ(JM),UNDEF,flat
      INTEGER i,j,k,kxlb,nt(kaij)
C****
C**** INITIALIZE CERTAIN QUANTITIES
C****
      UNDEF = -1.E30
      KXLB = INDEX(XLABEL(1:11),'(')-1
      IF(KXLB.le.0) KXLB = 10
      XLB = ' '
      XLB(1:13)=acc_period(1:3)//' '//acc_period(4:12)
      XLB(15:14+KXLB) = XLABEL(1:KXLB)
C**** Set up index for output index (temporary expedient)
      nt(ij_us)=63
      nt(ij_vs)=64
      nt(ij_taus)=65
      nt(ij_tauus)=66
      nt(ij_tauvs)=67
C****
C**** Print simple selected AIJ diagnostics
C****
      DO K=1,KAIJ
        IF (K.eq.IJ_US .or. K.eq.IJ_VS .or. K.eq.IJ_TAUS .or. K.eq
     *       .IJ_TAUUS .or. K.eq.IJ_TAUVS) THEN
          DO J=1,JM
            FLAT = 0.
            DO I=1,IM
              SMAP(I,J) =AIJ(I,J,K)/(IDACC(ia_ij(K))+1.E-20)
              FLAT =FLAT+AIJ(I,J,K)/(IDACC(ia_ij(K))+1.E-20)
            END DO
            SMAPJ(J) = FLAT*byIM
          END DO
          nt_ij=nt(K)
          CALL POUT_IJ(TITLE(nt_ij)//XLB,SMAP,SMAPJ,UNDEF,1)
        END IF
      END DO
C****
      RETURN
      END SUBROUTINE SIJMAP

      MODULE BDIJK
!@sum  stores information for outputting lon-lat-press diagnostics
!@auth M. Kelley

      IMPLICIT NONE

!@param nijk_out number of ijk-format output fields
      integer, parameter :: nijk_out=5

!@var title string, formed as concatentation of lname//units
      CHARACTER(LEN=48), DIMENSION(nijk_out) :: TITLE
!@var units string containing output field units
      CHARACTER(LEN=50), DIMENSION(nijk_out) :: UNITS
!@var lname string describing output field
      CHARACTER(LEN=50), DIMENSION(nijk_out) :: LNAME
!@var sname string referencing output field in self-desc. output file
      CHARACTER(LEN=30), DIMENSION(nijk_out) :: SNAME

!@var nt_ijk index telling pout_ijk which field is being output
      integer :: nt_ijk

      END MODULE BDIJK

      SUBROUTINE IJK_TITLES
      USE BDIJK
      IMPLICIT NONE
      INTEGER :: K
c
      k = 0
c
      k = k + 1
      sname(k) = 'u'
      lname(k) = 'ZONAL WIND (U COMPONENT)'
      units(k) = 'METERS/SECOND'
      k = k + 1
      sname(k) = 'v'
      lname(k) = 'MERIDIONAL WIND (V COMPONENT)'
      units(k) = 'METERS/SECOND'
      k = k + 1
      sname(k) = 'temp'
      lname(k) = 'TEMPERATURE'
      units(k) = 'deg C'
      k = k + 1
      sname(k) = 'q'
      lname(k) = 'SPECIFIC HUMIDITY'
      units(k) = '10**-5'
      k = k + 1
      sname(k) = 'height'
      lname(k) = 'HEIGHT'
      units(k) = 'm'

c create titles by concatenating long names with units
c no checks whether total length of lname+units exceeds length of title
      do k=1,nijk_out
         title(k)=''
         if(lname(k).ne.'')
     &        title(k) = trim(lname(k))//' ('//trim(units(k))//')'
      enddo

      return
      END SUBROUTINE IJK_TITLES

      SUBROUTINE IJKMAP
C****
C**** CALCULATE 3-D output fields on B-grid
C****
      USE CONSTANT, only :
     &     grav,sha
      USE MODEL_COM, only :
     &     im,jm,lm,
     &     PTOP,SIG,PSFMPT,XLABEL,LRUNID
      USE DAGCOM, only :
     &     aijk,acc_period,ijk_u,ijk_v,ijk_t,ijk_q,ijk_dp,ijk_dse
      USE BDIJK, only :
     &     nt_ijk
      IMPLICIT NONE

      CHARACTER XLB*24,TITLEX*56
      CHARACTER*80 TITLEL(LM)
      REAL*8 SMAP(IM,JM-1,LM),SMAPJK(JM-1,LM),SMAPK(LM)
      REAL*8 UNDEF,flat,press,dp
      CHARACTER*4 CPRESS(LM)
      INTEGER i,j,l,kxlb,ni

C**** OPEN PLOTTABLE OUTPUT FILE
      call open_ijk(trim(acc_period)//'.ijk'//XLABEL(1:LRUNID))

C****
C**** INITIALIZE CERTAIN QUANTITIES
C****
      call ijk_titles
      UNDEF = -1.E30
      SMAPK(:) = UNDEF
      KXLB = INDEX(XLABEL(1:11),'(')-1
      IF(KXLB.le.0) KXLB = 10
      XLB = ' '
      XLB(1:13)=acc_period(1:3)//' '//acc_period(4:12)
      XLB(15:14+KXLB) = XLABEL(1:KXLB)
C****
C**** Complete 3D-field titles  and select fields
C****
      DO L=1,LM
      PRESS = PTOP+PSFMPT*SIG(L)
      IF (PRESS.GE.1.) THEN
        WRITE(CPRESS(L),'(I4)') NINT(PRESS)
      ELSE
        WRITE(CPRESS(L),'(F4.3)') PRESS
      END IF
      END DO
C****
      TITLEX = '     U-WIND           at        mb (m/s, UV G)'
      DO 715 L=1,LM

      DO J=1,JM-1
      DO I=1,IM
      DP=AIJK(I,J+1,L,IJK_DP)
      SMAP(I,J,L) = UNDEF
      IF(DP.GT.0.) SMAP(I,J,L) = AIJK(I,J+1,L,IJK_U)/DP
      ENDDO
      ENDDO

      DO J=1,JM-1
        SMAPJK(J,L) = UNDEF
        NI = 0
        FLAT = 0.
        DO I=1,IM
        IF (SMAP(I,J,L).NE.UNDEF) THEN
          FLAT = FLAT+SMAP(I,J,L)
          NI = NI+1
        END IF
        END DO
        IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
      END DO
      WRITE(TITLEX(27:30),'(A)') CPRESS(L)
      TITLEL(L) = TITLEX//XLB
  715 CONTINUE
      nt_ijk = 1
      CALL POUT_IJK(TITLEL,SMAP,SMAPJK,SMAPK)
C****
      TITLEX = '     V-WIND           at        mb (m/s, UV G)'
      DO 725 L=1,LM

      DO J=1,JM-1
      DO I=1,IM
      DP=AIJK(I,J+1,L,IJK_DP)
      SMAP(I,J,L) = UNDEF
      IF(DP.GT.0.) SMAP(I,J,L) = AIJK(I,J+1,L,IJK_V)/DP
      ENDDO
      ENDDO

      DO J=1,JM-1
        SMAPJK(J,L) = UNDEF
        NI = 0
        FLAT = 0.
        DO I=1,IM
        IF (SMAP(I,J,L).NE.UNDEF) THEN
          FLAT = FLAT+SMAP(I,J,L)
          NI = NI+1
        END IF
        END DO
        IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
      END DO
      WRITE(TITLEX(27:30),'(A)') CPRESS(L)
      TITLEL(L) = TITLEX//XLB
  725 CONTINUE
      nt_ijk = 2
      CALL POUT_IJK(TITLEL,SMAP,SMAPJK,SMAPK)
C****
      TITLEX = '     TEMPERATURE      at        mb (C, UV Grid)'
      DO 735 L=1,LM

      DO J=1,JM-1
      DO I=1,IM
      DP=AIJK(I,J+1,L,IJK_DP)
      SMAP(I,J,L) = UNDEF
      IF(DP.GT.0.) SMAP(I,J,L) = (.25*AIJK(I,J+1,L,IJK_T)/DP - 273.16)
      ENDDO
      ENDDO

      DO J=1,JM-1
        SMAPJK(J,L) = UNDEF
        NI = 0
        FLAT = 0.
        DO I=1,IM
        IF (SMAP(I,J,L).NE.UNDEF) THEN
          FLAT = FLAT+SMAP(I,J,L)
          NI = NI+1
        END IF
        END DO
        IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
      END DO
      WRITE(TITLEX(27:30),'(A)') CPRESS(L)
      TITLEL(L) = TITLEX//XLB
  735 CONTINUE
      nt_ijk = 3
      CALL POUT_IJK(TITLEL,SMAP,SMAPJK,SMAPK)
C****
      TITLEX = '    SPECIFIC HUMIDITY at        mb (10**-5, UV G)'
      DO 745 L=1,LM

      DO J=1,JM-1
      DO I=1,IM
      DP=AIJK(I,J+1,L,IJK_DP)
      SMAP(I,J,L) = UNDEF
      IF(DP.GT.0.) SMAP(I,J,L) = .25*1.E5*AIJK(I,J+1,L,IJK_Q)/DP
      ENDDO
      ENDDO

      DO J=1,JM-1
        SMAPJK(J,L) = UNDEF
        NI = 0
        FLAT = 0.
        DO I=1,IM
        IF (SMAP(I,J,L).NE.UNDEF) THEN
          FLAT = FLAT+SMAP(I,J,L)
          NI = NI+1
        END IF
        END DO
        IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
      END DO
      WRITE(TITLEX(27:30),'(A)') CPRESS(L)
      TITLEL(L) = TITLEX//XLB
  745 CONTINUE
      nt_ijk = 4
      CALL POUT_IJK(TITLEL,SMAP,SMAPJK,SMAPK)
C****
      TITLEX = '       HEIGHT         at        mb (m, UV Grid)'
      DO 755 L=1,LM

      DO J=1,JM-1
      DO I=1,IM
      DP=AIJK(I,J+1,L,IJK_DP)
      SMAP(I,J,L) = UNDEF
      IF(DP.GT.0.) SMAP(I,J,L) = .25*(AIJK(I,J+1,L,IJK_DSE)-
     *   SHA*AIJK(I,J+1,L,IJK_T))/(GRAV*DP)
      ENDDO
      ENDDO

      DO J=1,JM-1
        SMAPJK(J,L) = UNDEF
        NI = 0
        FLAT = 0.
        DO I=1,IM
        IF (SMAP(I,J,L).NE.UNDEF) THEN
          FLAT = FLAT+SMAP(I,J,L)
          NI = NI+1
        END IF
        END DO
        IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
      END DO
      WRITE(TITLEX(27:30),'(A)') CPRESS(L)
      TITLEL(L) = TITLEX//XLB
  755 CONTINUE
      nt_ijk = 5
      CALL POUT_IJK(TITLEL,SMAP,SMAPJK,SMAPK)

      call close_ijk

      RETURN
      END SUBROUTINE IJKMAP
