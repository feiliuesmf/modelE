!@sum  RES_M24S Resolution info for 24 layer, 4x5 strat model
!@+    (M20b with 4 more layers between 984mb and 150mb)
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
!@var LS1 Layers LS1->LM: constant pressure levels, L<LS1: sigma levels
      INTEGER, PARAMETER :: IM=72,JM=46,LM=24, LS1=17

!@var PSF,PMTOP global mean surface, model top pressure  (mb)
!@var PTOP pressure at interface level sigma/const press coord syst (mb)
      REAL*8, PARAMETER :: PSF=984.d0, PTOP = 150.d0, PMTOP = .3d0
!@var PSFMPT,PSTRAT pressure due to troposhere,stratosphere
      REAL*8, PARAMETER :: PSFMPT=PSF-PTOP, PSTRAT=PTOP-PMTOP

!@var SIGE sigma levels at layer interfaces (1)
      REAL*8, PARAMETER, DIMENSION(LM+1) :: SIGE = ((/
     t     PSF, 974.d0, 959.d0, 934.d0, 900.d0, 857.d0, ! Pbot L=1,..
     t  804.d0, 744.d0, 677.d0, 607.d0, 535.d0, 463.d0, !      L=...
     t  393.d0, 325.d0, 260.d0, 200.d0,                 !      L=...
     1    PTOP,                                         !      L=LS1
     s  110.d0,  80.d0,   55.d0,   35.d0,    20.d0,     !      L=...
     s  10.d0,    3.d0,   PMTOP /)                      !      L=..,LM+1
     *  - PTOP)/(PSF-PTOP)
!!!!  Note:   sige(1)=1,  sige(ls1)=0,  sige(lm+1)=-pstrat/psfmpt

!@var SIG,DSIG,byDSIG mid point, depth, 1/depth of sigma levels (1)
      REAL*8, dimension(lm), parameter ::
     &     SIG    = (sige(1:lm)+sige(2:lm+1))*0.5d0,
     &     DSIG   =  sige(1:lm)-sige(2:lm+1),
     &     byDSIG =  1./DSIG

C**** KEP depends on whether stratos. EP flux diagnostics are calculated
C**** If dummy EPFLUX is used set KEP=0, otherwise KEP=21
!@param KEP number of lat/height E-P flux diagnostics
      INTEGER, PARAMETER :: KEP=21

C**** Based on model top, determine how much of stratosphere is resolved
C****         PMTOP >= 10 mb,    ISTRAT = 0
C**** 1 mb <= PMTOP <  10 mb,    ISTRAT = 1
C****         PMTOP <   1 mb,    ISTRAT = 2
      INTEGER, PARAMETER :: ISTRAT = 2

      END MODULE RESOLUTION

C**** The vertical resolution also determines whether
C**** stratospheric wave drag will be applied or not.
C**** This resolution is for a stratospheric model and so must be used
C**** in conjunction with the strat. modules
