#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_WATER
#endif
      MODULE STRAITS
!@sum  STRAITS ocean strait related variables 
!@+    RESOLUTION DEPENDENT: This version is for 72x46 - M
!@auth Gary Russell/Gavin Schmidt/Allegra LeGrande
!@ver  1.0
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
      USE OCN_TRACER_COM, only : ntm
#endif
      USE SEAICE, only : lmi
      USE OCEAN, only : lmo
      IMPLICIT NONE
      SAVE
C**** These values are highly resolution dependent
C****
C**** LGM 21 kyr 
C****     Strait           From         To       LM    Width
C****     ------           ----         --       --    -----
C****  X  Fury & Hecla  Above ground/under ice
C****  X  Nares         Above ground/under ice
C****  1  Gibraltar      35,32 EN    37,33 WS     5    25000
C****  X  English       Above ground/under ice
C****  X  Kattegat      Above ground/under ice
C****  X  Bosporous     Black Sea -> Lake
C****  2  Red Sea        44,29 ES    45,28 WN     5   250000
C****  3  Bab al Mandab  45,28 ES    46,27 WN     2    25000
C****  X  Hormuz         Above ground
C****  X  Malacca        Above ground
C**** 4  Korea          63,31 EN    63,33 WS     4    17000!85000 makes deep salty (Tsushima gone)
C****  X  Soya-kaikyo    Above ground

      INTEGER, PARAMETER :: NMST=4  !@param NMST no. of ocean straits
      integer, parameter :: n_STgbrltr=1 ,n_STrs=2,n_STbam=3,n_STkor=4
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
     *     WIST = (/  2.5d4, 2.5d5, 2.5d4, 1.7d5/)

!@var XST,YST local coordinates [-1,1] for strait entry/exit points
      REAL*8, DIMENSION(NMST,2) :: XST = RESHAPE( (/
     *     6d-1,  6d-1, 6d-1, -1d0, 
     *     -8d-1,-8d-1, -1d0, -6d1 /), (/NMST,2/) ),
     *     YST = RESHAPE( (/
     *     8d-1,  -8d-1,-8d-1, 1d0,
     *     -6d-1, 6d-1,  0d0, -8d-1 /), (/NMST,2/) )

!@var IST,JST i,j coordinates of ends of straits
      INTEGER, DIMENSION(NMST,2) ::
     *     IST = RESHAPE( (/
     *     35, 44, 45, 63,
     *     37, 45, 46, 63 /), (/NMST,2/) )
     *     , JST = RESHAPE( (/
     *     32, 29, 28, 31,
     *     33, 28, 27, 33 /), (/NMST,2/) )

!@var LMST no. of levels in strait
      INTEGER, DIMENSION(NMST) :: LMST = (/
     *     5,   5,  2,  4 /) ! 180m,144m,17m,100m ... given 120m slr to 0k

!@var name_st Names of straits
      CHARACTER*20, DIMENSION(NMST) :: NAME_ST = (/
     *     'Gibraltar     ', 'Red Sea       '
     *     ,'Bab al Mandab ', 'Korea         '/)

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
      USE STRAITS, only : nmst,n_STgbrltr,n_STrs,n_STbam,n_STkor
      IMPLICIT NONE
      REAL*8 SF(JM-1,0:LMO,0:4),FACST,OLNST(LMO,NMST)
      INTEGER :: L

C**** Gibralter: (35,32) to (37,33)
      SF(32,L,1) = SF(32,L,1) + OLNST(L+1,n_STgbrltr)*FACST
C**** Red Sea: (44,29) to (45,28)
      SF(28,L,3) = SF(28,L,3) - OLNST(L+1,n_STrs)*FACST
C**** Bab-al-Mandab: (45,28) to (46,27)
      SF(27,L,3) = SF(27,L,3) - OLNST(L+1,n_STbam)*FACST
C**** Korea: (63,31) to (63,33)
      SF(31,L,2) = SF(31,L,2) + OLNST(L+1,n_STkor)*FACST
      SF(32,L,2) = SF(32,L,2) + OLNST(L+1,n_STkor)*FACST

      RETURN
      END

      SUBROUTINE STRMIJ_STRAITS(J,SF,OLNST,FACST)
!@sum Add strait flow from to the Stream Function (from STRMIJ)
      USE OCEAN, only : im,jm,lmo
      USE STRAITS, only : nmst,lmst,n_STgbrltr,n_STrs,n_STbam,n_STkor
      IMPLICIT NONE
      REAL*8 SF(IM,JM),FACST,OLNST(LMO,NMST)
      INTEGER :: J

C**** Gibrater: (35,32) to (37,33)
      IF (J.eq.33) SF(35,33) = SF(35,33) + SUM(OLNST(1:LMST(n_STgbrltr)
     *     ,n_STgbrltr))*FACST
      IF (J.eq.33) SF(36,33) = SF(36,33) + SUM(OLNST(1:LMST(n_STgbrltr)
     *     ,n_STgbrltr))*FACST
C**** Red Sea: (44,29) to (45,28)
      IF (J.eq.29) SF(44,29) = SF(44,29) + SUM(OLNST(1:LMST(n_STrs)
     *     ,n_STrs))*FACST
C**** Bab-al-Mandab: (45,28) to (46,27)
      IF (J.eq.28) SF(45,28) = SF(45,28) + SUM(OLNST(1:LMST(n_STbam)
     *     ,n_STbam))*FACST

      RETURN
      END

      SUBROUTINE OTJ_STRAITS(X,SOLNST,SCALE,KQ)
!@sum Calculate transport through straits from latitude to another
!@+   within the same basin (from OTJOUT)
      USE OCEAN, only : jm
      USE STRAITS, only : nmst,n_STgbrltr,n_STrs,n_STbam,n_STkor
      IMPLICIT NONE
      REAL*8, INTENT(INOUT) :: X(0:JM,4,3)
      REAL*8, INTENT(IN) :: SCALE,SOLNST(NMST)
      INTEGER, INTENT(IN) ::  KQ

C**** Gibralter: (35,32) to (37,33)
      X(32,1,KQ) = X(32,1,KQ) + SOLNST(n_STgbrltr)*SCALE
      X(32,4,KQ) = X(32,4,KQ) + SOLNST(n_STgbrltr)*SCALE
C**** Red Sea: (44,29) to (45,28)
      X(28,3,KQ) = X(28,3,KQ) - SOLNST(n_STrs)*SCALE
      X(28,4,KQ) = X(28,4,KQ) - SOLNST(n_STrs)*SCALE
C**** Bab-al-Mandab: (45,28) to (46,27)
      X(27,3,KQ) = X(27,3,KQ) - SOLNST(n_STbam)*SCALE
      X(27,4,KQ) = X(27,4,KQ) - SOLNST(n_STbam)*SCALE
C**** Korea: (63,31) to (63,33)
      X(31,2,KQ) = X(31,2,KQ) + SOLNST(n_STkor)*SCALE
      X(31,4,KQ) = X(31,4,KQ) + SOLNST(n_STkor)*SCALE
      X(32,2,KQ) = X(32,2,KQ) + SOLNST(n_STkor)*SCALE
      X(32,4,KQ) = X(32,4,KQ) + SOLNST(n_STkor)*SCALE

C**** Calculate transport through straits from one basin to another
C****

      RETURN
      END
