!@sum  RES_M23 Resolution info for 23 layer, 4x5 strat model
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
      INTEGER, PARAMETER :: IM=72,JM=46,LM=23, LS1=12

!@var PSF,PMTOP global mean surface, model top pressure  (mb)
!@var PTOP pressure at interface level sigma/const press coord syst (mb)
      REAL*8, PARAMETER :: PSF=984.d0, PTOP = 150.d0, PMTOP=2.0576514d-3
!@var PSFMPT,PSTRAT pressure due to troposhere,stratosphere
      REAL*8, PARAMETER :: PSFMPT=PSF-PTOP, PSTRAT=PTOP-PMTOP

!@var SIGE sigma levels at layer interfaces (1)
      REAL*8, PARAMETER, DIMENSION(LM+1) :: SIGE = ((/
     t     PSF,   960.d0,  929.d0,  884.d0,  819.d0,    ! Pbot L=1,5
     t  710.d0,   570.d0,  425.d0,  314.d0,  245.d0,    !      L=6,10
     t  192.d0,                                         !      L=11
     1    PTOP,                                         !      L=LS1
     s  117.d0,   86.2d0,  56.2d0,  31.6d0,  17.8d0,    !      L=13-17
     s  10.0d0,   4.63d0,  1.46d0,  0.460995258d0,      !      L=18-21
     s  0.1449994968d0,    3.12002802d-2,    PMTOP /)   !      L=..,LM+1
     *  - PTOP)/(PSF-PTOP)
!!!!  Note:   sige(1)=1,  sige(ls1)=0,  sige(lm)=-pstrat/psfmpt

!@var SIG,DSIG,byDSIG mid point, depth, 1/depth of sigma levels (1)
      REAL*8, dimension(lm), parameter ::
     &     SIG    = (sige(1:lm)+sige(2:lm+1))*0.5d0,
     &     DSIG   =  sige(1:lm)-sige(2:lm+1),
     &     byDSIG =  1./DSIG

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
