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
!@var TAUSS optical depth from super-saturated clouds
      REAL*8, DIMENSION(IM,JM,LM) :: TAUSS
!@var TAUMC optical depth from moist-convective clouds
      REAL*8, DIMENSION(IM,JM,LM) :: TAUMC
!@var TTOLD,QTOLD previous potential temperature, humidity
      REAL*8, DIMENSION(IM,JM,LM) :: TTOLD,QTOLD
!@var SVLHX,SVLAT previous latent heat of evaporation
      REAL*8, DIMENSION(IM,JM,LM) :: SVLHX,SVLAT
!@var RHSAV previous relative humidity
      REAL*8, DIMENSION(IM,JM,LM) :: RHSAV
!@var SVWMX addition of liquid water from MC
      REAL*8, DIMENSION(IM,JM,LM) :: SVWMX
!@var CLDSAV previous cloud cover area (percent)
      REAL*8, DIMENSION(IM,JM,LM) :: CLDSAV
!@var PBLTOP height of PBL (m) (NOT USED)
      REAL*8, DIMENSION(IM,JM) :: PBLTOP

!@var CLDSS super-saturated cloud cover area (percent)
      REAL*8, DIMENSION(IM,JM,LM) :: CLDSS
!@var CLDMC moist convective cloud cover area (percent)
      REAL*8, DIMENSION(IM,JM,LM) :: CLDMC
!@var CSIZE effective cloud droplet radius (microns) 1:mc 2:ss
      REAL*8, DIMENSION(IM,JM,LM,2) :: CSIZE

!@var AQ1,AQ2 unused diagnostics?
      REAL*8 :: AQ1,AQ2
      DIMENSION AQ1(IM,JM,LM),AQ2(IM,JM,LM)

      END MODULE CLD01_COM_E001


