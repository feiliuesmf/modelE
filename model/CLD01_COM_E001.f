      MODULE CLD01_COM_E001
!@sum  CLD01_COM_E001 model variables for moist convction and
!@sum          large-scale condensation
!@auth M.S.Yao/T. Del Genio (modularisation by Gavin Schmidt)
!@ver  1.0 (taken from CB265)
      USE E001M12_COM, only : IM,JM,LM
C**** Note that we USE the NAMELIST set-able parameters from CLD01 only
C**** to be able to pass them to INPUT. and possibly write out to rsf
C**** and acc files
      USE CLD01, only : U00wtr,U00ice,LMCM

      IMPLICIT NONE
      SAVE
!@var PREC precipitation (kg/m^2)
      REAL*8, DIMENSION(IM,JM) :: PREC
!@var TPREC temperature of preciptiation (C)
      REAL*8, DIMENSION(IM,JM) :: TPREC
!@var PRECSS precipitation from super-saturation (kg/m^2)
      REAL*8, DIMENSION(IM,JM) :: PRECSS
!@var TTOLD,QTOLD previous potential temperature, humidity
      REAL*8, DIMENSION(IM,JM,LM) :: TTOLD,QTOLD
!@var SVLHX,SVLAT previous latent heat of evaporation
      REAL*8, DIMENSION(IM,JM,LM) :: SVLHX,SVLAT
!@var RHSAV previous relative humidity
      REAL*8, DIMENSION(IM,JM,LM) :: RHSAV
!@var CLDSAV previous cloud cover area (percent)
      REAL*8, DIMENSION(IM,JM,LM) :: CLDSAV
!@var PBLTOP height of PBL (m) (NOT USED)
      REAL*8, DIMENSION(IM,JM) :: PBLTOP

C**** variables saved for radiation calculations
!@var TAUSS optical depth from super-saturated clouds
      REAL*8, DIMENSION(LM,IM,JM) :: TAUSS
!@var TAUMC optical depth from moist-convective clouds
      REAL*8, DIMENSION(LM,IM,JM) :: TAUMC
!@var CLDSS super-saturated cloud cover area (percent)
      REAL*8, DIMENSION(LM,IM,JM) :: CLDSS
!@var CLDMC moist convective cloud cover area (percent)
      REAL*8, DIMENSION(LM,IM,JM) :: CLDMC
!@var CSIZMC,CSIZSS mc,ss effective cloud droplet radius (microns)
      REAL*8, DIMENSION(LM,IM,JM) :: CSIZMC,CSIZSS

      END MODULE CLD01_COM_E001


