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
!@var LS1 Layers LS1->LM: constant pressure levels, L<LS1: sigma levels
      INTEGER, PARAMETER :: IM=72,JM=46,LM=12, LS1=9

!@var PSF,PMTOP global mean surface, model top pressure  (mb)
!@var PTOP pressure at interface level sigma/const press coord syst (mb)
      REAL*8, PARAMETER :: PSF=984.d0, PMTOP = 10.d0, PTOP = 150.d0
!@var PSFMPT,PSTRAT pressure due to troposhere,stratosphere
      REAL*8, PARAMETER :: PSFMPT=PSF-PTOP, PSTRAT=PTOP-PMTOP

!@var SIGE sigma levels at layer interfaces (1)
      REAL*8, PARAMETER, DIMENSION(LM+1) :: SIGE = ((/
     t     PSF,   934.d0,  854.d0,  720.d0,  550.d0,    ! Pbot L=1,5
     t  390.d0,   285.d0,  210.d0,                      !      L=...
     1    PTOP,                                         !      L=LS1
     s  100.d0,    60.d0,   30.d0,   PMTOP /)           !      L=..,LM+1
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
      INTEGER, PARAMETER :: KEP=0

C**** Based on model top, determine how much of stratosphere is resolved
C****         PMTOP >= 10 mb,    ISTRAT = 0
C**** 1 mb <= PMTOP <  10 mb,    ISTRAT = 1
C****         PMTOP <   1 mb,    ISTRAT = 2
      INTEGER, PARAMETER :: ISTRAT = 0

      END MODULE RESOLUTION

C**** The vertical resolution also determines whether
C**** stratospheric wave drag will be applied or not.
C**** Hence also included here are some dummy routines for non-strat
C**** models.

      SUBROUTINE DUMMY_STRAT
!@sum DUMMY dummy routines for non-stratospheric models
C**** Dummy routines in place of STRATDYN
      ENTRY init_GWDRAG
      ENTRY GWDRAG
      ENTRY VDIFF
      ENTRY io_strat
C**** Dummy routines in place of STRAT_DIAG (EP flux calculations)
C**** Note that KEP=0 is set to zero above for the dummy versions.
      ENTRY EPFLUX
      ENTRY EPFLXI
      ENTRY EPFLXP
C****
      RETURN
      END SUBROUTINE DUMMY_STRAT
