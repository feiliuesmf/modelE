      MODULE LAKES_COM
!@sum  LAKES_COM model variables for Lake/Rivers module
!@auth Gavin Schmidt
!@ver  1.0
      USE E001M12_COM, only : IM,JM

      IMPLICIT NONE

!@var MWL mass of lake water (kg)
      REAL*8, DIMENSION(IM,JM) :: MWL 
!@var GML total enthalpy of lake (J)
      REAL*8, DIMENSION(IM,JM) :: GML 
!@var TLAKE temperature of lake (C)
      REAL*8, DIMENSION(IM,JM) :: TLAKE 
!@var TFL freezing temperature for lakes (=0 C)
      REAL*8, PARAMETER :: TFL = 0
!@var FLEADLK lead fraction for lakes
      REAL*8, PARAMETER :: FLEADLK=.05d0   ! = 0?
!@var T50 50 day mean temperature (used for estimating lake ice cover)
      REAL*8, SAVE,DIMENSION(IM,JM) :: T50

      END MODULE LAKES_COM
