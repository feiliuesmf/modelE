#include "rundeck_opts.h"

      MODULE STRAITS
!@sum  STRAITS ocean strait related variables
!@+    RESOLUTION DEPENDENT: This version is for 144x90 - F
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
      USE TRACER_COM, only : ntm
#endif
      USE SEAICE, only : lmi
      USE OCEAN, only : lmo
      IMPLICIT NONE
      SAVE
C**** These values are highly resolution dependent
C****
C****     Strait           From         To       LM    Width
C****     ------           ----         --       --    -----
C****  1  Fury & Hecla   38,80 ES    40,80 WN     2    20000
C****  2  Nares          45,85 EN    47,86 WS     5     5000
C****  3  Belle Isle     49,71 EN    51,71 WS     4    40000
C****  4  Gibraltar      70,63 EN    72,64 WS     5    25000
C****  5  Irish Sea      70,72 NW    70,74 ES     3    30000
C****  6  Kattegat       77,74 EN    78,73 WS     2    60000
C****  7  Bosporous      83,65 EN    84,66 WS     2     6000
C****  8  White Sea      88,78 EN    89,79 WS     2    40000
C****  9  Bab al Mandab  89,53 ES    90,52 WN     6    25000
C**** 10  Hormuz         94,58 ES    96,58 WN     2   100000  ! 50000 
C**** 11  Malacca       112,48 EN   114,46 WS     2    50000
C**** 12  Soya-kaikyo   128,68 EN   129,69 WS     2    40000
C****
      INTEGER, PARAMETER :: NMST=12  !@param NMST no. of ocean straits
!@var MMST mass of water in strait (kg)
!@var MUST mass flux of water in strait (kg/s)
!@var G0MST,GXMST,GYMST pot. enthalpy of water in strait (+ moments) (J)
!@var S0MST,SXMST,SYMST salinity of water in strait (+ moments) (kg)
      REAL*8, DIMENSION(LMO,NMST) :: MMST,MUST,G0MST,GXMST,GZMST,S0MST
     *     ,SXMST,SZMST

!@var WIST width of strait (m)
!@var DIST distance along strait (m)
!@var DISTPG distance between centre points of adjoining ocean boxes (m)
      REAL*8, DIMENSION(NMST) :: DIST,DISTPG,
     *     WIST = (/  2d4,  5d3,  4d4, 2.5d4,   3d4,  6d4,
     *                6d3,  4d4,2.5d4,   1d5,   5d4,  4d4 /)

!@var XST,YST local coordinates [-1,1] for strait entry/exit points
      REAL*8, DIMENSION(NMST,2) :: XST = RESHAPE( (/
     *     .6d0,   0d0,   1d0,  .8d0,  .6d0,  .8d0,  .8d0,  .8d0,
     *     .8d0,   1d0,  .8d0,  .8d0, -.6d0,  -1d0,  -1d0, -.8d0,
     *      0d0, -.8d0, -.8d0, -.8d0, -.8d0,  -1d0, -.8d0, -.8d0
     *     /), (/NMST,2/) ),
     *     YST = RESHAPE( (/
     *     .8d0,   1d0,   0d0,  .8d0,   1d0, -.8d0,  .8d0,  .8d0,
     *    -.8d0,   0d0, -.8d0,  .8d0,  .8d0,   0d0,  .6d0, -.8d0,
     *     -1d0,  .8d0, -.8d0, -.8d0,  .8d0,   0d0,  .8d0, -.8d0
     *     /), (/NMST,2/) )

!@var IST,JST i,j coordinates of ends of straits
      INTEGER, DIMENSION(NMST,2) ::
     *     IST = RESHAPE( (/
     *     38, 45, 49, 70, 70, 77, 83, 88, 89, 94,112,128,
     *     40, 47, 51, 72, 70, 78, 84, 89, 90, 96,114,129/),
     *     (/NMST,2/) ),
     *     JST = RESHAPE( (/
     *     80, 85, 71, 63, 72, 74, 65, 78, 53, 58, 48, 68,
     *     80, 86, 71, 64, 74, 73, 66, 79, 52, 58, 46, 69/),
     *     (/NMST,2/) )

!@var LMST no. of levels in strait
      INTEGER, DIMENSION(NMST) :: LMST = (/
     *      2,  5,  4,  5,  3,  2,  2,  2,  6,  2,  2,  2/)

!@var name_st Names of straits
      CHARACTER*20, DIMENSION(NMST) :: NAME_ST = (/
     *     'Fury & Hecla  ', 'Nares         ', 'Belle Isle    ',
     *     'Gibraltar     ', 'Irish Sea     ', 'Kattegat      ',
     *     'Bosporous     ', 'White Sea     ', 'Bab al Mandab ',
     *     'Hormuz        ', 'Malacca       ', 'Soya-kaikyo   '/)

!@var RSIST Sea ice fraction in strait
!@var RSIXST Center of sea ice in strait (m)
!@var MSIST Mass of ice within strait (kg)
!@var HSIST Enthalpy of ice within strait (J)
!@var SSIST Salinity of ice within strait (kg)
      REAL*8, DIMENSION(NMST) :: RSIST,RSIXST
      REAL*8, DIMENSION(2,NMST) :: MSIST
      REAL*8, DIMENSION(LMI,NMST) :: HSIST,SSIST

!@param USIFAC ratio of strait sea ice velocity to current
      REAL*8 :: USIFAC = 0.1d0   ! used to be 1. (too much)

#ifdef TRACERS_OCEAN
!@var TRMST,TXMST,TZMST tracer amount in strait (+ moments) (kg)
      REAL*8, DIMENSION(LMO,NMST,NTM) :: TRMST, TXMST, TZMST
#endif
#ifdef TRACERS_WATER
!@var TRSIST tracer amount in with strait (kg)
      REAL*8, DIMENSION(NTM,LMI,NMST) :: TRSIST
#endif

      END MODULE STRAITS

C**** Resolution dependent strait diagnostic routines
C**** Called from routines in ocean diagnostics

      SUBROUTINE STRMJL_STRAITS(L,SF,OLNST,FACST)
!@sum Add strait flow from to the Stream Function (from STRMJL)
      USE OCEAN, only : jm,lmo
      USE STRAITS, only : nmst
      IMPLICIT NONE
      REAL*8 SF(JM-1,0:LMO,0:4),FACST,OLNST(LMO,NMST)
      INTEGER :: L

C**** Nares: (45,85) to (47,86)
      SF(85,L,1) = SF(85,L,1) + OLNST(L+1,2)*FACST
C**** Gibraltar: (70,63) to (72,64)
      SF(63,L,1) = SF(63,L,1) + OLNST(L+1,4)*FACST
C**** Irish Sea: (70,72) to (70,74)
      SF(72,L,1) = SF(72,L,1) + OLNST(L+1,5)*FACST
      SF(73,L,1) = SF(73,L,1) + OLNST(L+1,5)*FACST
C**** Kattegut (77,74) to (78,73)
      SF(73,L,1) = SF(73,L,1) - OLNST(L+1,6)*FACST
C**** Bosporous: (83,65) to (84,66)
      SF(65,L,1) = SF(65,L,1) + OLNST(L+1,7)*FACST
C**** White Sea: (88,78) to (89,79)
      SF(78,L,1) = SF(78,L,1) + OLNST(L+1,8)*FACST
C**** Bab-al-Mandab: (89,53) to (90,52)
      SF(52,L,3) = SF(52,L,3) - OLNST(L+1,9)*FACST
C**** Soya: (128,68) to (129,69)
      SF(68,L,2) = SF(68,L,2) + OLNST(L+1,12)*FACST
C**** Malacca: (112,48) to (114,46), from Indian Ocean to Pacific Ocean
      SF(46,L,2) = SF(46,L,2) - OLNST(L+1,11)*FACST   !*.5
c      DO 510 J=1,46
c      SF(J,L,2) = SF(J,L,2) - OLNST(L+1,11)*FACST*.5
c  510 SF(J,L,3) = SF(J,L,3) + OLNST(L+1,11)*FACST*.5
C****
      RETURN
      END

      SUBROUTINE STRMIJ_STRAITS(J,SF,OLNST,FACST)
!@sum Add strait flow from to the Stream Function (from STRMIJ)
      USE OCEAN, only : im,jm,lmo
      USE STRAITS, only : nmst,lmst
      IMPLICIT NONE
      REAL*8 SF(IM,JM),FACST,OLNST(LMO,NMST)
      INTEGER :: J

C**** Fury & Hecla: (38,80) to (40,80)
      IF (J.eq.80) SF(38,80) = SF(38,80) + SUM(OLNST(1:LMST(1),1))*FACST
      IF (J.eq.80) SF(39,80) = SF(39,80) + SUM(OLNST(1:LMST(1),1))*FACST
C**** Nares: (45,85) to (47,86)
      IF (J.eq.86) SF(45,86) = SF(45,86) + SUM(OLNST(1:LMST(2),2))*FACST
      IF (J.eq.86) SF(46,86) = SF(46,86) + SUM(OLNST(1:LMST(2),2))*FACST
C**** Gibraltar: (70,63) to (72,64)
      IF (J.eq.64) SF(70,64) = SF(70,64) + SUM(OLNST(1:LMST(4),4))*FACST
      IF (J.eq.64) SF(71,64) = SF(71,64) + SUM(OLNST(1:LMST(4),4))*FACST
C**** Kattegat: (77,74) to (78,73)
      IF (J.eq.74) SF(77,74) = SF(77,74) + SUM(OLNST(1:LMST(6),6))*FACST
C**** Bosporous: (83,65) to (84,66)
      IF (J.eq.66) SF(83,66) = SF(83,66) + SUM(OLNST(1:LMST(7),7))*FACST
C**** White Sea: (88,78) to (89,79)
      IF (J.eq.79) SF(88,79) = SF(88,79) + SUM(OLNST(1:LMST(8),8))*FACST
C**** Bab-al-Mandab: (89,53) to (90,52)
      IF (J.eq.53) SF(89,53) = SF(89,53) + SUM(OLNST(1:LMST(9),9))*FACST
C**** Hormuz: (94,58) to (96,58)
      IF (J.eq.58) SF(94,58) = SF(94,58)+SUM(OLNST(1:LMST(10),10))*FACST
      IF (J.eq.58) SF(95,58) = SF(95,58)+SUM(OLNST(1:LMST(10),10))*FACST
C**** Soya: (128,68) to (129,69)
      IF (J.eq.69) SF(128,69)=SF(128,69)+SUM(OLNST(1:LMST(12),12))*FACST
C**** Malacca: (112,48) to (114,46)
      IF (J.eq.48) SF(112,48)=SF(112,48)+SUM(OLNST(1:LMST(11),11))*FACST
      IF (J.eq.48) SF(113,48)=SF(113,48)+SUM(OLNST(1:LMST(11),11))*FACST
C****
      RETURN
      END

      SUBROUTINE OTJ_STRAITS(X,SOLNST,SCALE,KQ)
!@sum Calculate transport through straits from latitude to another
!@+   within the same basin (from OTJOUT)
      USE OCEAN, only : jm
      USE STRAITS, only : nmst
      IMPLICIT NONE
      REAL*8, INTENT(INOUT) :: X(0:JM,4,3)
      REAL*8, INTENT(IN) :: SCALE,SOLNST(NMST)
      INTEGER, INTENT(IN) ::  KQ

C**** Nares: (45,85) to (47,86)
      X(85,1,KQ) = X(85,1,KQ) + SOLNST(2)*SCALE
      X(85,4,KQ) = X(85,4,KQ) + SOLNST(2)*SCALE
C**** Gibraltar: (70,63) to (72,64)
      X(63,1,KQ) = X(63,1,KQ) + SOLNST(4)*SCALE
      X(63,4,KQ) = X(63,4,KQ) + SOLNST(4)*SCALE
C**** Irish Sea: (70,72) to (70,74)
      X(72,1,KQ) = X(72,1,KQ) + SOLNST(5)*SCALE
      X(72,4,KQ) = X(72,4,KQ) + SOLNST(5)*SCALE
      X(73,1,KQ) = X(73,1,KQ) + SOLNST(5)*SCALE
      X(73,4,KQ) = X(73,4,KQ) + SOLNST(5)*SCALE
C**** Kattegut (77,74) to (78,73)
      X(73,1,KQ) = X(73,1,KQ) - SOLNST(6)*SCALE
      X(73,4,KQ) = X(73,4,KQ) - SOLNST(6)*SCALE
C**** Bosporous: (83,65) to (84,66)
      X(65,1,KQ) = X(65,1,KQ) + SOLNST(7)*SCALE
      X(65,4,KQ) = X(65,4,KQ) + SOLNST(7)*SCALE
C**** White Sea: (88,78) to (89,79)
      X(78,1,KQ) = X(78,1,KQ) + SOLNST(8)*SCALE
      X(78,4,KQ) = X(78,4,KQ) + SOLNST(8)*SCALE
C**** Bab-al-Mandab: (89,53) to (90,52)
      X(52,3,KQ) = X(52,3,KQ) - SOLNST(9)*SCALE
      X(52,4,KQ) = X(52,4,KQ) - SOLNST(9)*SCALE
C**** Soya: (128,68) to (129,69)
      X(68,2,KQ) = X(68,2,KQ) + SOLNST(12)*SCALE
      X(68,4,KQ) = X(68,4,KQ) + SOLNST(12)*SCALE
C****
C**** Calculate transport through straits from one basin to another
C****
C**** Malacca: (112,48) to (114,46)
      X(46,2,KQ) = X(46,2,KQ) - SOLNST(11)*SCALE   !*.5
      X(46,4,KQ) = X(46,4,KQ) - SOLNST(11)*SCALE   !*.5
      X(47,2,KQ) = X(47,2,KQ) - SOLNST(11)*SCALE   !*.5
      X(47,4,KQ) = X(47,4,KQ) - SOLNST(11)*SCALE   !*.5
c      DO 510 J=1,46
c      X(J,2,KQ) = X(J,2,KQ) - SOLNST(11)*SCALE*.5
c      X(J,4,KQ) = X(J,4,KQ) - SOLNST(11)*SCALE*.5
c      X(J,3,KQ) = X(J,3,KQ) + SOLNST(11)*SCALE*.5
c  510 X(J,4,KQ) = X(J,4,KQ) + SOLNST(11)*SCALE*.5
C****
      RETURN
      END
