      MODULE DYNAMICS
!@sum  DYNAMICS contains all the routines for calulating dynamics
!@sum           related variables
!@auth Original development team
!@ver  1.0
!@cont ADVECM,ADVECV,AFLUX,PGF,SHAP1D,AVRX,FILTER,FLTRUV,CALC_AMPK

      USE E001M12_COM, only : im,jm,lm,imh,sig,sige,dsig,psf,ptop,ls1,u
     *     ,v,t,q,p,wm,mfiltr,zatmo,fim,mrch,modd5k,psfmpt,bydsig,byim
      USE CONSTANT, only : grav,rgas,kapa,sday,lhm,lhe,lhs,twopi,omega,
     *     bbyg,gbyrb,bykapa,bykapap1,bykapap2
      USE SOMTQ_COM, only : tmom,qmom
      USE PBLCOM, only : tsavg

      USE GEOM
      IMPLICIT NONE
      SAVE
C**** Some helpful arrays (arrays should be L first)
!@var  PLIJ  Surface pressure: P(I,J) or PSF-PTOP (mb)
      REAL*8, DIMENSION(LM,IM,JM) :: PLIJ
!@var  PDSIG  Surface pressure * DSIG(L) (mb)
      REAL*8, DIMENSION(LM,IM,JM) :: PDSIG
!@var  AM  Air mass of each box (kg/m^2)
!      REAL*8, DIMENSION(LM,IM,JM) :: AM     ! PLIJ*DSIG(L)*100/grav
!@var  BYAM  1/Air mass (m^2/kg)
!      REAL*8, DIMENSION(LM,IM,JM) :: BYAM
!@var  PMID  Pressure at mid point of box (mb)
      REAL*8, DIMENSION(LM,IM,JM) :: PMID    ! SIG(L)*PLIJ+PTOP
!@var  PK   PMID**KAPA
      REAL*8, DIMENSION(LM,IM,JM) :: PK
!@var  PEUP  Pressure at lower edge of box (incl. surface) (mb)
      REAL*8, DIMENSION(LM+1,IM,JM) :: PEDN  ! SIGE(L)*PLIJ+PTOP
!@var  PEK  PEUP**KAPA
      REAL*8, DIMENSION(LM+1,IM,JM) :: PEK
!@var  SQRTP  square root of P (used in diagnostics)
      REAL*8, DIMENSION(IM,JM) :: SQRTP

C**** module should own dynam variables used by other routines
!@var PTOLD pressure at beginning of dynamic time step (for clouds)
      REAL*8, DIMENSION(IM,JM)    :: PTOLD
!@var SD_CLOUDS vert. integrated horizontal convergence (for clouds)
      REAL*8, DIMENSION(IM,JM,LM) :: SD_CLOUDS
!@var GZ geopotential height (for Clouds and Diagnostics)
      REAL*8, DIMENSION(IM,JM,LM) :: GZ
!@var DPDX,DPDY surface pressure gradients (for PBL)
c      REAL*8, SAVE,DIMENSION(IM,JM)    :: DPDX,DPDY


      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: PU,PV,CONV
      DOUBLE PRECISION, DIMENSION(IM,JM,LM-1) :: SD
      DOUBLE PRECISION, DIMENSION(IM,JM) :: PIT
      EQUIVALENCE (SD(1,1,1),CONV(1,1,2))
      EQUIVALENCE (PIT(1,1),CONV(1,1,1))

      END MODULE DYNAMICS

