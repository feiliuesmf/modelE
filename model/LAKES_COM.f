      MODULE LAKES_COM
!@sum  LAKES_COM model variables for Lake/Rivers module
!@auth Gavin Schmidt
!@ver  1.0
      USE E001M12_COM, only : IM,JM
c      USE LAKES, only : LMLI

      IMPLICIT NONE

!@var MWL mass of lake water (kg)
      REAL*8, DIMENSION(IM,JM) :: MWL 
!@var GML total enthalpy of lake (J)
      REAL*8, DIMENSION(IM,JM) :: GML 
!@var TSL surface temp. of lake (C)
      REAL*8, DIMENSION(IM,JM) :: TSL 
!@var TFL freezing temperature for lakes (=0 C)
      REAL*8, PARAMETER :: TFL = 0
!@var FLEADLK lead fraction for lakes
      REAL*8, PARAMETER :: FLEADLK=.05d0   ! = 0?
!@var RLI lake ice fraction
c      REAL*8, DIMENSION(IM,JM) :: RLI 
!@var MLI mass of lake ice (kg/m^2)
c      REAL*8, DIMENSION(IM,JM,2) :: MLI 
!@var HLI enthalpy of lake ice (J/m^2)
c      REAL*8, DIMENSION(IM,JM,LMLI) :: HLI 

      END MODULE LAKES_COM
