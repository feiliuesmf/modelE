#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_WATER
#endif

      Module STRAITS
C****
!@sum   STRAITS ocean strait related variables
!@+     RESOLUTION DEPENDENT: This version is for 288x180 - S
!@auth  Gary Russell/Gavin Schmidt
!@ver   1.0
C****
      Use OCEAN,  Only: LMO
      Use SEAICE, Only: LMI

#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
      Use OCN_TRACER_COM, Only: NTM
#endif

      Implicit None
      Save
C****
C****     Strait              From         To         LM    Width
C****     ------              ----         --         --    -----
C****  1  Dolphin & Union   51,160 ES    53,159  N     3    32000
C****  2  Dease             57,159 EN    61,159 WN     3    25000
C****  3  Fury & Hecla      76,160 EN    79,160 W      2    15000
C****  4  Nares             91,171 EN    93,172 WS     6    30000
C****  5  Gibraltar        139,126 EN   141,126 WS     7    15000
C****  6  English          145,141 EN   146,142 WS     2    35000
C****  7  Dardanelles      165,130 EN   167,131 WS     2     1000
C****  8  Bosporous        167,131 EN   168,132 WS     2     2000
C****  9  Bab al Mandab    179,104 ES   180,103 WN     5    25000
C**** 10  Malacca          225, 93 ES   227, 92 WN     2    40000
C**** 11  Salat Sunda      229, 84 EN   230, 85 WS     2    25000
C****
C**** Fixed parameters of straits
C****      NMST    = number of ocean straits
C****  IST(NMST,2) = I index of west and east ocean cells
C****  JST(NMST,2) = J index of west and east ocean cells
C**** LMST(NMST)   = number of layers in strait
C****  XST(NMST,2) = X location of strait end in ocean cell [-1:1]
C****  YST(NMST,2) = Y location of strait end in ocean cell [-1:1]
C**** WIST(NMST)   = width of strait (m)
C****
      Integer*4,Parameter :: NMST =11,
     *  IST(NMST,2) = Reshape ( (/
     *     51, 57, 76, 91,139,145,165,167,179,225,229,
     *     53, 61, 79, 93,141,146,167,168,180,227,230/),(/NMST,2/) ),
     *  JST(NMST,2) = Reshape ( (/
     *    160,159,160,171,126,141,130,131,104, 93, 84,
     *    159,159,160,172,126,142,131,132,103, 92, 85/),(/NMST,2/) ),
     *  LMST(NMST) = (/ 3, 2, 2, 6, 7, 2, 2, 2, 4, 2, 2/)
      Real*8,Parameter ::
     *  XST(NMST,2) = Reshape ( (/
     *     .7d0, .7d0, .7d0, .7d0, .7d0, .7d0,
     *           .7d0, .7d0, .7d0, .7d0, .7d0,
     *     .0d0,-.7d0,-1.d0,-.7d0,-.7d0,-.7d0,
     *          -.7d0,-.7d0,-.7d0,-.7d0,-.7d0 /), (/ NMST,2 /) ),
     *  YST(NMST,2) = Reshape ( (/
     *    -.7d0, .7d0, .7d0, .7d0, .7d0, .7d0,
     *           .7d0, .7d0,-.7d0,-.7d0, .7d0,
     *     1.d0, .7d0, .0d0,-.7d0,-.7d0,-.7d0,
     *          -.7d0,-.7d0, .7d0, .7d0,-.7d0 /), (/ NMST,2 /) ),
     *  WIST(NMST) = (/ 32d3, 25d3, 15d3, 30d3, 15d3, 35d3,
     *                   1d3,  2d3, 25d3, 40d3, 25d3 /)
      Character*20,Parameter :: NAME_ST(NMST) = (/
     *  'Dolphin & Union', 'Dease  ', 'Fury & Hecla', 'Nares',
     *  'Gibraltar'      , 'English', 'Dardanelles' , 'Bosporous',
     *  'Bab al Mandab'  , 'Malacca', 'Salat Sunda' /)

      Real*8 ::
     *   MMST(LMO,NMST), !  mass of water in strait (kg)
     *   MUST(LMO,NMST), !  mass flux of water in strait (kg/s)
     *  G0MST(LMO,NMST), !  pot. enthalpy of water in strait (J)
     *  GXMST(LMO,NMST), !  west-east gradient of pot. enthalpy (J)
     *  GZMST(LMO,NMST), !  down. vert. gradient of pot. enthalpy (J)
     *  S0MST(LMO,NMST), !  mass of salt in strait (kg)
     *  SXMST(LMO,NMST), !  west-east gradient of salt (kg)
     *  SZMST(LMO,NMST), !  down. vert. gradient of salt (kg)
     *       DIST(NMST), !  fixed length of strait (m)
     *     DISTPG(NMST), !  fixed distance between ocean cell centers
     *      RSIST(NMST), !  horizontal sea ice cover in strait
     *     RSIXST(NMST), !  center of sea ice cover in strait
     *    MSIST(2,NMST), !  2 mass layers of sea ice in strait (kg)
     *  HSIST(LMI,NMST), !  LMI layers of heat content in strait (J)
     *  SSIST(LMI,NMST), !  LMI layers of salt in strait (kg)
     *    USIFAC = .1d0  !  ratio of strait sea ice velocity to current

#ifdef TRACERS_OCEAN
      Real*8 TRMST(LMO,NMST,NTM),TXMST(LMO,NMST,NTM),TZMST(LMO,NMST,NTM)
#endif
#ifdef TRACERS_WATER
      Real*8 TRSIST(NTM,LMI,NMST)
#endif

      EndModule STRAITS

C**** Resolution dependent strait diagnostic routines
C**** Called from routines in ocean diagnostics

      Subroutine STRMJL_STRAITS (L,SF,OLNST,FACST)
!@sum Add strait flow from to the Stream Function (from STRMJL)
      Use OCEAN,   Only: JM,LMO
      Use STRAITS, Only: NMST
      Implicit  None
      Real*8    SF(JM-1,0:LMO,0:4), FACST, OLNST(LMO,NMST)
      Integer*4 L
C**** Dolphin & Union: (51,160) to (53,159)
      SF(159,L,1) = SF(159,L,1) - OLNST(L+1,1)*FACST
C**** Nares: (91,171) to (93,172)
      SF(171,L,1) = SF(171,L,1) + OLNST(L+1,4)*FACST
C**** English: (145,141) to (146,142)
      SF(141,L,1) = SF(141,L,1) + OLNST(L+1,6)*FACST
C**** Dardanelles: (165,130) to (167,131)
      SF(130,L,1) = SF(130,L,1) + OLNST(L+1,7)*FACST
C**** Bosporous: (167,131) to (168,132)
      SF(131,L,1) = SF(131,L,1) + OLNST(L+1,8)*FACST
C**** Bab al Mandab: (179,104) to (180,103)
      SF(103,L,3) = SF(103,L,3) - OLNST(L+1,9)*FACST
C**** Malacca: (225,93) to (227,92), from Indian to Pacific Ocean
      SF( 92,L,2) = SF( 92,L,2) - OLNST(L+1,10)*FACST*.5
      SF( 92,L,3) = SF( 92,L,3) - OLNST(L+1,10)*FACST*.5
C**** Salat Sunda: (229,84) to (230,85), from Indian to Pacific Ocean
      SF( 84,L,2) = SF( 84,L,2) + OLNST(L+1,11)*FACST*.5
      SF( 84,L,3) = SF( 84,L,3) + OLNST(L+1,11)*FACST*.5
      Return
      End

      Subroutine STRMIJ_STRAITS (J,SF,OLNST,FACST)
!@sum Add strait flow from to the Stream Function (from STRMIJ)
      Use OCEAN,   Only: IM,JM,LMO
      Use STRAITS, Only: NMST,LMST
      Implicit  None
      Real*8    SF(IM,JM),FACST,OLNST(LMO,NMST)
      Integer*4 J
C**** Dolphin & Union: (51,160) to (53,159)
      If (J ==159)  Then
        SF(51,159) = SF(51,159) + Sum(OLNST(1:LMST(1),1))*FACST
        SF(52,159) = SF(52,159) + Sum(OLNST(1:LMST(1),1))*FACST
C**** Dease: (57,159) to (61,159)
        SF(57,159) = SF(57,159) + Sum(OLNST(1:LMST(2),2))*FACST
        SF(58,159) = SF(58,159) + Sum(OLNST(1:LMST(2),2))*FACST
        SF(59,159) = SF(59,159) + Sum(OLNST(1:LMST(2),2))*FACST
        SF(60,159) = SF(60,159) + Sum(OLNST(1:LMST(2),2))*FACST ; EndIf
C**** Fury & Hecla: (76,160) to (79,160)
      If (J ==160)  Then
        SF(76,160) = SF(76,160) + Sum(OLNST(1:LMST(3),3))*FACST
        SF(77,160) = SF(77,160) + Sum(OLNST(1:LMST(3),3))*FACST
        SF(78,160) = SF(78,160) + Sum(OLNST(1:LMST(3),3))*FACST ; EndIf
C**** Nares: (91,171) to (93,172)
      If (J ==172)  Then
        SF(91,172) = SF(91,172) + Sum(OLNST(1:LMST(4),4))*FACST
        SF(92,172) = SF(92,172) + Sum(OLNST(1:LMST(4),4))*FACST ; EndIf
C**** Gibraltar: (139,126) to (141,126)
      If (J == 126)  Then
        SF(139,126) = SF(139,122) + Sum(OLNST(1:LMST(5),5))*FACST
        SF(140,126) = SF(140,122) + Sum(OLNST(1:LMST(5),5))*FACST ;EndIf
C**** English: (145,141) to (146,142)
      If (J == 142)
     *  SF(145,142) = SF(145,142) + Sum(OLNST(1:LMST(6),6))*FACST
C**** Dardanelles: (165,130) to (167,131)
      If (J == 131)  Then
        SF(165,131) = SF(165,131) + Sum(OLNST(1:LMST(7),7))*FACST
        SF(166,131) = SF(166,131) + Sum(OLNST(1:LMST(7),7))*FACST ;EndIf
C**** Bosporous: (167,131) to (168,132)
      If (J == 132)
     *  SF(167,132) = SF(167,132) + Sum(OLNST(1:LMST(8),8))*FACST
C**** Bab-al-Mandab: (179,104) to (180,103)
      If (J == 103)
     *  SF(179,103) = SF(179,103) + Sum(OLNST(1:LMST(9),9))*FACST
C**** Malacca: (225,93) to (227,92)
      If (J == 92)  Then
        SF(225,92) = SF(225,92) + Sum(OLNST(1:LMST(10),10))*FACST
        SF(226,92) = SF(226,92) + Sum(OLNST(1:LMST(10),10))*FACST ;EndIf
C**** Salat Sunda: (229,84) to (230,85)
      If (J == 85)
     *  SF(229,85) = SF(229,85) + Sum(OLNST(1:LMST(11),11))*FACST
      Return
      End

      Subroutine OTJ_STRAITS (X,SOLNST,SCALE,KQ)
!@sum Calculate transport through straits from latitude to another
!@+   within the same basin (from OTJOUT)
      Use OCEAN,   Only: JM
      Use STRAITS, Only: NMST
      Implicit NONE
      Real*8,Intent(InOut) :: X(0:JM,4,3)
      Real*8,   Intent(In) :: SCALE, SOLNST(NMST)
      Integer*4,Intent(In) :: KQ
C**** Dolphin & Union (51,160) to (53,159)
      X(159,1,KQ) = X(159,1,KQ) - SOLNST(1)*SCALE
      X(159,4,KQ) = X(159,4,KQ) - SOLNST(1)*SCALE
C**** Nares: (91,171) to (93,172)
      X(171,1,KQ) = X(171,1,KQ) + SOLNST(4)*SCALE
      X(171,4,KQ) = X(171,4,KQ) + SOLNST(4)*SCALE
C**** English: (145,141) to (146,142)
      X(141,1,KQ) = X(141,1,KQ) + SOLNST(6)*SCALE
      X(141,4,KQ) = X(141,4,KQ) + SOLNST(6)*SCALE
C**** Dardanelles: (165,130) to (167,131)
      X(130,1,KQ) = X(130,1,KQ) + SOLNST(7)*SCALE
      X(130,4,KQ) = X(130,4,KQ) + SOLNST(7)*SCALE
C**** Bosporous: (167,131) to (168,132)
      X(131,1,KQ) = X(131,1,KQ) + SOLNST(8)*SCALE
      X(131,4,KQ) = X(131,4,KQ) + SOLNST(8)*SCALE
C**** Bab-al-Mandab: (179,104) to (180,103)
      X(103,3,KQ) = X(103,3,KQ) - SOLNST(9)*SCALE
      X(103,4,KQ) = X(103,4,KQ) - SOLNST(9)*SCALE
C**** Malacca: (225,93) to (227,92)
      X(92,2,KQ) = X(92,2,KQ) - SOLNST(10)*SCALE*.5
      X(92,3,KQ) = X(92,3,KQ) - SOLNST(10)*SCALE*.5
      X(92,4,KQ) = X(92,4,KQ) - SOLNST(10)*SCALE
C**** Salat Sunda: (229,84) to (230,85)
      X(84,2,KQ) = X(84,2,KQ) + SOLNST(11)*SCALE*.5
      X(84,3,KQ) = X(84,3,KQ) + SOLNST(11)*SCALE*.5
      X(84,4,KQ) = X(84,4,KQ) + SOLNST(11)*SCALE
      Return
C****
!@param NMST no. of ocean straits
!@var IST,JST i,j coordinates of ends of straits
!@var LMST no. of levels in strait
!@var XST,YST local coordinates [-1,1] for strait entry/exit points
!@var WIST width of strait (m)
!@var NAME_ST Names of straits
!@var MMST mass of water in strait (kg)
!@var MUST mass flux of water in strait (kg/s)
!@var G0MST,GXMST,GZMST pot. enthalpy of water in strait (+ moments) (J)
!@var S0MST,SXMST,SZMST salinity of water in strait (+ moments) (kg)
!@var DIST distance along strait (m)
!@var DISTPG distance between centre points of adjoining ocean boxes (m)
!@var RSIST Sea ice fraction in strait
!@var RSIXST Center of sea ice in strait (m)
!@var MSIST Mass of ice within strait (kg)
!@var HSIST Enthalpy of ice within strait (J)
!@var SSIST Salinity of ice within strait (kg)
!@param USIFAC ratio of strait sea ice velocity to current
!@var TRMST,TXMST,TZMST tracer amount in strait (+ moments) (kg)
!@var TRSIST tracer amount in with strait (kg)
      End
