      MODULE RADNCB
      
      USE E001M12_COM, only : im,jm,lm

      IMPLICIT NONE
      SAVE

      INTEGER, PARAMETER :: LM_REQ=3 ! # radiation-only layers
      DOUBLE PRECISION, DIMENSION(IM,JM,LM_REQ) :: RQT
      DOUBLE PRECISION, DIMENSION(IM,JM,LM+1) :: SRHR,TRHR
      DOUBLE PRECISION, DIMENSION(IM,JM,4) :: FSF

      DOUBLE PRECISION, DIMENSION(IM,JM) :: COSZ1

!@var S0X solar constant multiplication factor (set-able in NAMELIST)
      REAL*8 :: S0X = 1.
!@var CO2 carbon dioxide multiplication factor (set-able in NAMELIST)
      REAL*8 :: CO2 = 1.

      END MODULE RADNCB
