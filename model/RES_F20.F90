!@sum RES_F20.F90   Resolution file, 2x2.5 Lat-Lon Grid, 20 layers, top at .1 mb, no GWDRAG
!@ver 2013/03/21
!@auth Original Development Team

      Module RESOLUTION
      Implicit None
!@var IM,JM = longitudinal and latitudinal number of grid cells
!@var LM    = number of dynamical layers
!@var LS1   = lowest layer of strtosphere
      Integer*4,Parameter :: IM=144,JM=90,LM=20, LS1=11

!@var MDRYA = dry atmospheric mass (kg/m^2) = 100*PSF/GRAV
!@var MTOP  = mass above dynamical top (kg/m^2) = 100*PMTOP/GRAV
!@var MFIXs = summation of MFIX (kg/m^2) = 100*(PTOP-PMTOP)/GRAV
!@var MVAR  = spatially and temporally varying column mass (kg/m^2)
!@var MSURF = MTOP + MFIXs + MVAR (kg/m^2)
!@var AM(L) = MFIX(L) + MVAR*MFRAC(L) (kg/m^2)
!@var MFIX(L)  = fixed mass in each layer (kg/m^2) = 100*[PLBOT(L)-PLBOT(L+1)]/GRAV
!@var MFRAC(L) = fraction of variable mass in each layer = DSIG(L)
      Real*8,Parameter :: MDRYA = 98400/9.80665d0, MTOP = 10/9.80665d0, MFIXs = 14990/9.80665d0, &
         MFIX(LM) = (/ 0d0,0d0,0d0,0d0,0d0, 0d0,0d0,0d0,0d0,0d0,  &
                       4000/9.80665d0, 3000/9.80665d0, 2500/9.80665d0, 2000/9.80665d0, 1500/9.80665d0, &
                      1000/9.80665d0,  700/9.80665d0,  200/9.80665d0,   70/9.80665d0,   20/9.80665d0 /), &
         MFRAC(LM) = (/ 20/834d0, 30/834d0, 50/834d0, 74/834d0,100/834d0, &
                       160/834d0,160/834d0,105/834d0, 75/834d0, 60/834d0, & 
                       0d0,0d0,0d0,0d0,0d0, 0d0,0d0,0d0,0d0,0d0 /)      

!**** Vertival resolution
!****                         ---MSURF=10034.0---    ---MSURF=5781.8----
!**** Layer   MFIX   MFRAC    MVAR     AM     SUM    MVAR     AM     SUM
!**** =====   ====   =====    ====     ==     ===    ====     ==     ===
!****   Top                                   1.0                    1.0
!****    20    2.0    0        0.0    2.0             0.0    2.0
!****    19    7.1    0        0.0    7.1             0.0    7.1
!****    18   20.4    0        0.0   20.4             0.0   20.4
!****
!****    12  305.9    0        0.0  203.9             0.0  305.9
!****    11  407.9    0        0.0  224.3  1529.6     0.0  407.9  1529.6
!****    10    0.0  60/834   611.8  611.8           305.9  305.9
!****     9    0.0  70/834
!****
!****     3    0.0  50/834
!****     2    0.0  30/834
!****     1    0.0  20/834   203.9  203.9 10034.0   102.0  102.0  5781.8
!****         ----  ------   -----  -----           -----  -----
!****       1528.6 834/834  8504.4 10033.0         4252.2 5780.8  

!@var PSF,PMTOP global mean surface, model top pressure  (mb)
!@var PTOP pressure at interface level sigma/const press coord syst (mb)
!@var PSFMPT,PSTRAT pressure due to troposhere,stratosphere
!@var PLbot pressure levels at bottom of layers (mb)
      Real*8,Parameter :: PSF=984.d0, PTOP = 150.d0, PMTOP = .1d0, PSFMPT = PSF-PTOP, &
                          PSTRAT = PTOP-PMTOP, &
         PLBOT(1:LM+1) = (/ PSF,   964d0, 934d0, 884d0, 810d0, &  ! Pbot L=1,..
                            710d0, 550d0, 390d0, 285d0, 210d0, &  !      L=...
                            PTOP,                              &  !      L=LS1
                            110d0,  80d0,  55d0,  35d0,  20d0, &  !      L=...
                             10d0,   3d0,   1d0,  .3d0, PMTOP /)  !      L=..,LM+1

!**** KEP depends on whether stratos. EP flux diagnostics are calculated
!**** If dummy EPFLUX is used set KEP=0, otherwise KEP=21
!@param KEP number of lat/height E-P flux diagnostics
      Integer*4,Parameter :: KEP = 0

!**** Based on model top, determine how much of stratosphere is resolved
!**** ISTRAT = 2:          PMTOP <   1 mb
!**** ISTRAT = 1:  1 mb <= PMTOP <  10 mb
!**** ISTRAT = 0: 10 mb <= PMTOP
      Integer*4,Parameter :: ISTRAT = 2

      EndModule RESOLUTION


      Subroutine DUMMY_STRAT
!**** Dummy routines in place of STRATDYN
!**** The vertical resolution also determines whether stratospheric wave drag will be applied or not.
!**** Hence also included here are some dummy routines for non-stratospheric models.
!@sum DUMMY dummy routines for non-stratospheric models
      Entry INIT_GWDRAG
      Entry GWDRAG
      Entry VDIFF
      Entry io_strat
      Entry ALLOC_STRAT_COM
!**** Dummy routines in place of STRAT_DIAG (EP flux calculations)
!**** Note that KEP=0 is set to zero above for the dummy versions.
      Entry EPFLUX
      Entry EPFLXI
      Entry EPFLXP
      Return
      EndSubroutine DUMMY_STRAT

