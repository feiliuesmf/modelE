!@sum  RES_M12 Resolution info for 12 layer, 4x5 non-strat model
!@auth Original Development Team
!@ver  1.0

      MODULE RESOLUTION
!@sum  RESOLUTION contains horiz/vert resolution variables
!@auth Original Development Team
!@ver  1.0
      IMPLICIT NONE
      SAVE
!@var IM,JM longitudinal and latitudinal number of grid boxes
!@var LM number of vertical levels 
      INTEGER, PARAMETER :: IM=72,JM=46,LM=12

      END MODULE RESOLUTION

C**** The vertical resolution also determines whether
C**** stratospheric wave drag will be applied or not.
C**** Hence also included here are some dummy routines for non-strat
C**** models. 

      SUBROUTINE DUMMY_STRAT
!@sum DUMMY dummy routines for non-stratospheric models
      ENTRY GWDRAG
      ENTRY VDIFF
      ENTRY EPFLUX
      ENTRY EPFLXI
      ENTRY diaga0
      ENTRY io_strat
      RETURN
      END SUBROUTINE DUMMY_STRAT
