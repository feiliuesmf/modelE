!@sum  RES_M23 Resolution info for 23 layer, 4x5 strat model
!@auth Jean Lerner
!@ver  1.0

      MODULE RESOLUTION
!@sum  RESOLUTION contains horiz/vert resolution variables
!@auth Original Development Team
!@ver  1.0
      IMPLICIT NONE
      SAVE
!@var IM,JM longitudinal and latitudinal number of grid boxes
!@var LM number of vertical levels 
      INTEGER, PARAMETER :: IM=72,JM=46,LM=23

      END MODULE RESOLUTION

