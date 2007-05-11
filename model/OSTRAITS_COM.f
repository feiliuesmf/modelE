#include "rundeck_opts.h"

      MODULE STRAITS
!@sum  STRAITS ocean strait related variables 
!@+    RESOLUTION DEPENDENT: This version is for 72x46 - M
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
C****  1  Fury & Hecla   19,42 ES    20,40 WN     2    20000
C****  2  Nares          22,43 EN    24,44 WS     5     5000
C****  3  Gibraltar      35,32 EN    37,33 WS     5    25000
C****  4  English        36,36 EN    37,37 WS     2    35000
C****  5  Kattegat       38,38 EN    40,38 WS     2    60000
C****  6  Bosporous      42,33 EN    43,34 WS     2     6000
C****  7  Red Sea        44,29 ES    45,28 WN     6   250000
C****  8  Bab al Mandab  45,28 ES    46,27 WN     6    25000
C****  9  Hormuz         47,30 ES    49,29 WN     2   100000  ! 50000
C**** 10  Malacca        56,25 EN    58,24 WS     3    50000
C**** 11  Korea          62,32 EN    63,33 WS     4   170000
C**** 12  Soya-kaikyo    64,34 EN    65,35 WS     2    40000
C**** 12alt Tsugaru      64,34 EN    66,34 WS     4    20000
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
     *     WIST = (/  2d4,   5d3, 2.5d4, 3.5d4,   6d4, 6d3,
     *              2.5d5, 2.5d4,   1d5,   5d4, 1.7d5, 4d4/)

!@var XST,YST local coordinates [-1,1] for strait entry/exit points
      REAL*8, DIMENSION(NMST,2) :: XST = RESHAPE( (/
     *     0d0  , 6d-1, 6d-1,  1d0,  1d0,  0d0, 6d-1, 6d-1,  1d0, 6d-1,
     *     0d0  , 6d-1,
     *     0d0  ,-8d-1,-8d-1,  0d0,-6d-1,-8d-1,-8d-1, -1d0,-7d-1, -1d0,
     *     -6d-1, -1d0 /), (/NMST,2/) ),
     *     YST = RESHAPE( (/
     *     -1d0 , 8d-1, 8d-1,  0d0,  0d0,  1d0,-8d-1,-8d-1,  0d0,-8d-1,
     *      1d0 , 8d-1,
     *      1d0 ,-6d-1,-6d-1, -1d0,-8d-1,-6d-1, 6d-1,  0d0, 7d-1,  0d0,
     *     -8d-1,  0d0 /), (/NMST,2/) )

!@var IST,JST i,j coordinates of ends of straits
      INTEGER, DIMENSION(NMST,2) ::
     *     IST = RESHAPE( (/
     *     19, 22, 35, 36, 38, 42, 44, 45, 47, 56, 62, 64,
     *     20, 24, 37, 37, 40, 43, 45, 46, 49, 58, 63, 65/),
     *     (/NMST,2/) ),
     *     JST = RESHAPE( (/
     *     42, 43, 32, 36, 38, 33, 29, 28, 30, 25, 32, 34,
     *     40, 44, 33, 37, 38, 34, 28, 27, 29, 24, 33, 35/),
     *     (/NMST,2/) )

!@var LMST no. of levels in strait
      INTEGER, DIMENSION(NMST) :: LMST = (/
     *     2,  5,  5,  2,  2,  2,  6,  6,  2,  3,  4,  2/)

!@var name_st Names of straits
      CHARACTER*20, DIMENSION(NMST) :: NAME_ST = (/
     *     'Fury & Hecla  ', 'Nares         ', 'Gibraltar     ',
     *     'English       ', 'Kattegat      ', 'Bosporous     ',
     *     'Red Sea       ', 'Bab al Mandab ', 'Hormuz        ',
     *     'Malacca       ', 'Korea         ', 'Soya-kaikyo   '/)

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

C**** Fury & Hecla: (19,42) to (20,40)
      SF(41,L,1) = SF(41,L,1) - OLNST(L+1,1)*FACST
      SF(40,L,1) = SF(40,L,1) - OLNST(L+1,1)*FACST
C**** Nares: (22,43) to (24,44)
      SF(43,L,1) = SF(43,L,1) + OLNST(L+1,2)*FACST
C**** Gibralter: (35,32) to (37,33)
      SF(32,L,1) = SF(32,L,1) + OLNST(L+1,3)*FACST
C**** English: (36,36) to (37,37)
      SF(36,L,1) = SF(36,L,1) + OLNST(L+1,4)*FACST
C**** Bosporous: (42,33) to (43,34)
      SF(33,L,1) = SF(33,L,1) + OLNST(L+1,6)*FACST
C**** Red Sea: (44,29) to (45,28)
      SF(28,L,3) = SF(28,L,3) - OLNST(L+1,7)*FACST
C**** Bab-al-Mandab: (45,28) to (46,27)
      SF(27,L,3) = SF(27,L,3) - OLNST(L+1,8)*FACST
C**** Hormuz: (47,30) to (49,29)
      SF(29,L,3) = SF(29,L,3) - OLNST(L+1,9)*FACST
C**** Korea: (62,32) to (63,33)
      SF(32,L,2) = SF(32,L,2) + OLNST(L+1,11)*FACST
C**** Soya: (64,34) to (65,35)
      SF(34,L,2) = SF(34,L,2) + OLNST(L+1,12)*FACST
C**** Malacca: (56,25) to (58,24), from Indian Ocean to Pacific Ocean
      SF(24,L,2) = SF(24,L,2) - OLNST(L+1,10)*FACST  !*.5
c      DO 510 J=1,23
c      SF(J,L,2) = SF(J,L,2) - OLNST(L+1,10)*FACST*.5
c  510 SF(J,L,3) = SF(J,L,3) + OLNST(L+1,10)*FACST*.5
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

C**** Fury & Hecla: (19,42) to (20,40)
      IF (J.eq.41) SF(19,41) = SF(19,41) + SUM(OLNST(1:LMST(1),1))*FACST
C**** Nares: (22,43) to (24,44)
      IF (J.eq.44) SF(22,44) = SF(22,44) + SUM(OLNST(1:LMST(2),2))*FACST
      IF (J.eq.44) SF(23,44) = SF(23,44) + SUM(OLNST(1:LMST(2),2))*FACST
C**** Gibrater: (35,32) to (37,33)
      IF (J.eq.33) SF(35,33) = SF(35,33) + SUM(OLNST(1:LMST(3),3))*FACST
      IF (J.eq.33) SF(36,33) = SF(36,33) + SUM(OLNST(1:LMST(3),3))*FACST
C**** Engish: (36,36) to (37,37)
      IF (J.eq.37) SF(36,37) = SF(36,37) + SUM(OLNST(1:LMST(4),4))*FACST
C**** Bosporous: (42,33) to (43,34)
      IF (J.eq.34) SF(42,34) = SF(42,34) + SUM(OLNST(1:LMST(6),6))*FACST
C**** Red Sea: (44,29) to (45,28)
      IF (J.eq.29) SF(44,29) = SF(44,29) + SUM(OLNST(1:LMST(7),7))*FACST
C**** Bab-al-Mandab: (45,28) to (46,27)
      IF (J.eq.28) SF(45,28) = SF(45,28) + SUM(OLNST(1:LMST(8),8))*FACST
C**** Hormuz: (47,30) to (49,29)
      IF (J.eq.30) SF(47,30) = SF(47,30) + SUM(OLNST(1:LMST(9),9))*FACST
      IF (J.eq.30) SF(48,30) = SF(48,30) + SUM(OLNST(1:LMST(9),9))*FACST
C**** Korea: (62,32) to (63,33)
      IF (J.eq.33) SF(62,33) = SF(62,33)+SUM(OLNST(1:LMST(11),11))*FACST
C**** Soya: (64,34) to (65,35)
      IF (J.eq.35) SF(64,35) = SF(64,35)+SUM(OLNST(1:LMST(12),12))*FACST
C**** Malacca: (56,25) to (58,24),
      IF (J.eq.25) SF(56,25) = SF(56,25)+SUM(OLNST(1:LMST(10),10))*FACST
      IF (J.eq.25) SF(57,25) = SF(57,25)+SUM(OLNST(1:LMST(10),10))*FACST
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

C**** Fury & Hecla: (19,42) to (20,40)
      X(40,1,KQ) = X(40,1,KQ) - SOLNST(1)*SCALE
      X(40,4,KQ) = X(40,4,KQ) - SOLNST(1)*SCALE
      X(41,1,KQ) = X(41,1,KQ) - SOLNST(1)*SCALE
      X(41,4,KQ) = X(41,4,KQ) - SOLNST(1)*SCALE
C**** Nares: (22,43) to (24,44)
      X(43,1,KQ) = X(43,1,KQ) + SOLNST(2)*SCALE
      X(43,4,KQ) = X(43,4,KQ) + SOLNST(2)*SCALE
C**** Gibralter: (35,32) to (37,33)
      X(32,1,KQ) = X(32,1,KQ) + SOLNST(3)*SCALE
      X(32,4,KQ) = X(32,4,KQ) + SOLNST(3)*SCALE
C**** English: (36,36) to (37,37)
      X(36,1,KQ) = X(36,1,KQ) + SOLNST(4)*SCALE
      X(36,4,KQ) = X(36,4,KQ) + SOLNST(4)*SCALE
C**** Bosporous: (42,33) to (43,34)
      X(33,1,KQ) = X(33,1,KQ) + SOLNST(6)*SCALE
      X(33,4,KQ) = X(33,4,KQ) + SOLNST(6)*SCALE
C**** Red Sea: (44,29) to (45,28)
      X(28,3,KQ) = X(28,3,KQ) - SOLNST(7)*SCALE
      X(28,4,KQ) = X(28,4,KQ) - SOLNST(7)*SCALE
C**** Bab-al-Mandab: (45,28) to (46,27)
      X(27,3,KQ) = X(27,3,KQ) - SOLNST(8)*SCALE
      X(27,4,KQ) = X(27,4,KQ) - SOLNST(8)*SCALE
C**** Hormuz: (47,30) to (49,29)
      X(29,3,KQ) = X(29,3,KQ) - SOLNST(9)*SCALE
      X(29,4,KQ) = X(29,4,KQ) - SOLNST(9)*SCALE
C**** Korea: (62,32) to (63,33)
      X(32,2,KQ) = X(32,2,KQ) + SOLNST(11)*SCALE
      X(32,4,KQ) = X(32,4,KQ) + SOLNST(11)*SCALE
C**** Soya: (64,34) to (65,35)
      X(34,2,KQ) = X(34,2,KQ) + SOLNST(12)*SCALE
      X(34,4,KQ) = X(34,4,KQ) + SOLNST(12)*SCALE
C****
C**** Calculate transport through straits from one basin to another
C****
C**** Malacca: (56,25) to (58,24)
c      DO 510 J=1,23
c      X( J,2,KQ) = X( J,2,KQ) + SOLNST(10)*SCALE
c  510 X( J,3,KQ) = X( J,3,KQ) - SOLNST(10)*SCALE
      X(24,2,KQ) = X(24,2,KQ) - SOLNST(10)*SCALE
      X(24,4,KQ) = X(24,4,KQ) - SOLNST(10)*SCALE

      RETURN
      END
