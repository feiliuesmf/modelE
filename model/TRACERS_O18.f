#include "rundeck_opts.h"

!@sum  TRACERS_O18 water isotope specific routines and functions
!@auth Gavin Schmidt

      FUNCTION FRACVL(TEMP,trname)
!@sum FRACVL Calculate vapour --> liquid equilibrium fractionation factor
!@+          by a quadratic fit to Majoube'71
!@auth Gavin Schmidt
      IMPLICIT NONE
!@var TEMP  temperature (deg C)
      REAL*8, INTENT(IN) :: TEMP
      CHARACTER, INTENT(IN) :: trname*8
      INTEGER, PARAMETER :: NTSPM=4
      REAL*8, PARAMETER ::
     *     A(NTSPM) = (/  0d0, -3.75d-7, -6.375d-6,  0d0 /),
     *     B(NTSPM) = (/  0d0,  1.025d-4, 1.2475d-3, 0d0 /),
     *     C(NTSPM) = (/  1d0,  0.9884d0, 0.9001d0 , 1d0 /)
!@var ITR   species number of tracer
C****       1: Fresh Water 2: o18 3: Deu  4: Tritium
      INTEGER ITR
      REAL*8 FRACVL
C****
      select case (trname)
      case ('Water')
        ITR=1
      case ('H2O18')
        ITR=2
      case ('HDO')
        ITR=3
      case ('HTO')   ! tritium data not yet confirmed
        ITR=4
      case default
        write(6,*) "Tracer name ",trname," not defined in FRACVL"
        stop
      end select
      FRACVL=C(ITR) + TEMP*(B(ITR) + TEMP*A(ITR))
C****
      RETURN
      END FUNCTION FRACVL

      FUNCTION FRACVS(TEMP,trname)
!@sum FRACVS Calculate vapour --> solid (ice) equilibrium fractionation
!@+          factor
!@auth Gavin Schmidt
      IMPLICIT NONE
!@var TEMP  temperature (deg C)
      REAL*8, INTENT(IN) :: TEMP
      CHARACTER, INTENT(IN) :: trname*8
      INTEGER, PARAMETER :: NTSPM=4
      REAL*8, PARAMETER ::
     *     A(NTSPM) = (/  0d0,  1.36d-4,   1.46d-3, 0d0/),
     *     B(NTSPM) = (/  1d0,  0.9850d0, 0.8834d0, 1d0/)
!@var ITR   species number of tracer
C****       1: Fresh Water 2: o18 3: Deu  4: Tritium
      INTEGER ITR
      REAL*8 FRACVS
C****
      select case (trname)
      case ('Water')
        ITR=1
      case ('H2O18')
        ITR=2
      case ('HDO')
        ITR=3
      case ('HTO')   ! tritium data not yet confirmed
        ITR=4
      case default
        write(6,*) "Tracer name ",trname," not defined in FRACVS"
        stop
      end select
      FRACVS=B(ITR) + A(ITR)*TEMP
C****
      RETURN
      END FUNCTION FRACVS

       FUNCTION FRACLS(trname)
!@sum FRACLS Calculate liquid --> solid equilibrium fractionation factor
!@auth Gavin Schmidt
      IMPLICIT NONE
      CHARACTER, INTENT(IN) :: trname*8
      INTEGER, PARAMETER :: NTSPM=4
      REAL*8, PARAMETER ::
     *     A(NTSPM) = (/ 1d0, 1.0035d0, 1.0208d0,  1.04d0/)
!@var ITR   species number of tracer
C****       1: Fresh Water 2: o18 3: Deu  4: Tritium
      INTEGER ITR
      REAL*8 FRACLS
C****
      select case (trname)
      case ('Water')
        ITR=1
      case ('H2O18')
        ITR=2
      case ('HDO')
        ITR=3
      case ('HTO')   ! tritium data not yet confirmed
        ITR=4
      case default
        FRACLS=1.    ! no fractionation
        RETURN
      end select
      FRACLS=A(ITR)
C****
      RETURN
      END FUNCTION FRACLS

      FUNCTION FRACLK(WS,trname)
!@sum FRACLK calculates the liquid/vapor kinetic fractionation factor
!@+          from either (Merlivat and Jouzel,1972) or as a constant
!@auth Gavin Schmidt
      IMPLICIT NONE
!@var WS surface wind speed (m/s)
      REAL*8, INTENT(IN) :: WS
      CHARACTER, INTENT(IN) :: trname*8
      INTEGER, PARAMETER :: NTSPM=4
      REAL*8, PARAMETER ::
     *     A(NTSPM) = (/ 1d0, 0.994d0,  0.99472d0,  0.98944d0  /),
     *     B(NTSPM) = (/ 0d0, 0.285d-3, 0.2508d-3,  0.5016d-3  /),
     *     C(NTSPM) = (/ 0d0, 0.82d-3,  0.7216d-3,  0.14432d-2 /)
!@var ITR   species number of tracer
C****       1: Fresh Water 2: o18 3: Deu  4: Tritium
      INTEGER ITR
      REAL*8 FRACLK
C****
C**** Calculate kinetic fractionation factor:
C****
      select case (trname)
      case ('Water')
        ITR=1
      case ('H2O18')
        ITR=2
      case ('HDO')
        ITR=3
      case ('HTO')   ! tritium data not yet confirmed
        ITR=4
      case default
        write(6,*) "Tracer name ",trname," not defined in FRACLK"
        stop
      end select
      FRACLK=A(ITR)
      IF(WS.GE.7.) FRACLK=1d0-(B(ITR)*WS+C(ITR))
C****
      RETURN
      END FUNCTION FRACLK

#ifdef TRACERS_SPECIAL_O18
      SUBROUTINE ISOEQUIL(N,TEMP,QMV,QML,TRMV,TRML)
!@sum ISOEQUIL equilibrates isotopes in vapour and liquid/solid reservoirs
!@auth Gavin Schmidt/Georg Hoffmann
      USE CONSTANT, only : tf
      USE TRACER_COM, only: trname,tr_wd_TYPE,nWATER
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      REAL*8, INTENT(IN) :: TEMP  !@var TEMP temperature (K)
!@var QMV,QML vapour and liquid water masses (kg or kg/m^2)
      REAL*8, INTENT(IN) :: QMV,QML  
!@var TRMV,TRML vapour and liquid tracer mass (kg or kg/m^2)
      REAL*8, INTENT(INOUT) :: TRMV,TRML
      REAL*8 TDEGC,ZALPH,ZDELEQU,ZXFAC,FRACVL,FRACVS

      SELECT CASE(tr_wd_TYPE(N))
        CASE(nWATER)
          TDEGC = TEMP - TF
          IF (QMV.GT.0.) THEN
            IF (TDEGC.GE.0.) THEN
              ZALPH = 1./FRACVL(TDEGC,trname(N)) 
            ELSE
              ZALPH = 1./FRACVS(TDEGC,trname(N))
            END IF
            ZXFAC = ZALPH*QML/QMV
            ZDELEQU = (TRML - ZXFAC*TRMV)/(1.+ZXFAC)
            TRML = TRML - ZDELEQU
            TRMV = TRMV + ZDELEQU
          END IF
      END SELECT

      RETURN
C****
      END SUBROUTINE ISOEQUIL
#endif
