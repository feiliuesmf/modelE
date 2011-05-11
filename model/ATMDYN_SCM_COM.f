      MODULE DYNAMICS
!@sum  DYNAMICS contains all the pressure and momentum related variables
!@auth Original development team
      USE DOMAIN_DECOMP_ATM, ONLY : grid
      USE RESOLUTION , ONLY : LS1, LM
      IMPLICIT NONE
      SAVE

!**** Vertical resolution dependent variables (set in INPUT)
!@var SIGE sigma levels at layer interfaces (1)
!!!!  Note:   sige(1)=1,  sige(ls1)=0,  sige(lm+1)=-pstrat/psfmpt
      REAL*8, DIMENSION(LM+1) :: SIGE
!@var SIG,DSIG,byDSIG mid point, depth, 1/depth of sigma levels (1)
      REAL*8, DIMENSION(LM) ::
     &     SIG,    ! = (sige(1:lm)+sige(2:lm+1))*0.5d0,
     &     DSIG,   ! =  sige(1:lm)-sige(2:lm+1),
     &     byDSIG  ! =  1./DSIG

!@dbparam DT (atmospheric) dynamics time step (s)
      REAL*8 :: DT    =  450.         ! DT = DTdyn_atm
!@var     NIdyn:  DT atm_dyn  =  DTsrc/NIdyn     (NIdyn=DTsrc/DT)
      INTEGER :: NIdyn

!@dbparam NFILTR: DT_filter   =  NFILTR*DTsrc
      INTEGER :: NFILTR = 1

!@var NSTEP number of dynamics steps since start of run
!@var MRCH  flags position in dynamics cycle (>0 fw, <0 bw step)
      INTEGER :: NSTEP,MRCH

!     Controls for FLTRUV (momentum/velocity filter)
!@dbparam DT_XUfilter dU is multiplied by dt/DT_XUfilter in E-W
!@dbparam DT_XVfilter dV is multiplied by dt/DT_XVfilter in E-W
!@dbparam DT_YUfilter dU is multiplied by dt/DT_YUfilter in N-S
!@dbparam DT_YVfilter dV is multiplied by dt/DT_YVfilter in N-S
      REAL*8 :: ! a value of 0. switches off a given filter
     &     DT_XUfilter=0.,DT_XVfilter=0.,DT_YUfilter=0.,DT_YVfilter=0.
!@var QUVfilter: True if any of DT_[XY][UV]filter are not=0
      LOGICAL :: QUVfilter
!@dbparam do_polefix if =1 : u,v tendencies are corrected near the pole
      INTEGER :: do_polefix=1     ! default is to enable corrections

!@dbparam MFILTR: if 1 => SLP, if 2 => T, if 3 => SLP&T is filtered
      integer :: MFILTR = 1
!@dbparam ang_uv =1 to conserve ang mom in UVfilter
      INTEGER :: ang_uv = 1 ! UV-filter conserves ang mom

!**** Stratospheric drag related parameters
!@dbparam X_SDRAG.  SDRAG ~X_SDRAG(1)+X_SDRAG(2)*wind_magnitude
      REAL*8, DIMENSION(2) :: X_SDRAG = (/2.5D-4,2.5D-5/)
!@dbparam C_SDRAG.  SDRAG=C_SDRAG (const.) above PTOP
      REAL*8 :: C_SDRAG = 2.5D-5
      REAL*8, DIMENSION(LS1:LM) :: CSDRAGL
!@dbparam P_CSDRAG pressure level above which const.drag is increased
      REAL*8 :: P_CSDRAG=0.
!@dbparam P(P)_SDRAG pressure level above which SDRAG is applied (mb)
      REAL*8 :: P_SDRAG=0., PP_SDRAG = 1.d0 ! (PP_... near poles)
!@var L(P)SDRAG level above which SDRAG is applied (near pole)
      INTEGER :: LSDRAG=LM, LPSDRAG=LM  ! non-polar, polar limit
!@var ANG_SDRAG if =1: ang.momentum lost by SDRAG is added in below PTOP
      INTEGER :: ANG_SDRAG=1  ! default: SDRAG does conserve ang.mom
!@dbparam Wc_JDRAG critical velocity for J.Hansen/Judith Perlwitz drag
      REAL*8 :: Wc_JDRAG=30.d0  !  if 0.: no JDRAG-feature in Sdrag
!@dbparam wmax imposed limit for stratospheric winds (m/s) in SDRAG
      real*8 :: wmax=200.d0
!@dbparam VSDRAGL layer dependent tuning factor for stratospheric drag
!@+   (not =1 e.g. if used with explicit grav.wave drag scheme)
      real*8 :: VSDRAGL(LS1:LM) = 1d0
!@dbparam  USE_UNR_DRAG   if 1 =>SDRAG is turned off and GWD is applied.
!@+    if 0 => SDRAG is kept intact and alternative GWD is not employed.
      INTEGER :: USE_UNR_DRAG=0  ! default: SDRAG is kept intact.

C**** Variables specific for stratosphere and/or strat diagnostics
!@var DO_GWDRAG when true, prints Gravity Wave diagnostics
      LOGICAL :: DO_GWDRAG = .false.
!@var iDO_GWDRAG number if AIJ Gravity wave diagnostics
      INTEGER :: iDO_GWDRAG = 0


!@var  PK   PMID**KAPA
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PK
!@var GZ geopotential height (for Clouds and Diagnostics)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: GZ

C**** Some helpful arrays (arrays should be L first)
!@var  PLIJ  Surface pressure: P(I,J) or PSF-PTOP (mb)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PLIJ
!@var  PDSIG  Surface pressure * DSIG(L) (mb)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PDSIG
!@var  AM  Air mass of each box (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: AM     ! PLIJ*DSIG(L)*100/grav
!@var  BYAM  1/Air mass (m^2/kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: BYAM
!@var  PMID  Pressure at mid point of box (mb)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PMID    ! SIG(L)*PLIJ+PTOP
!@var  PEUP  Pressure at lower edge of box (incl. surface) (mb)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PEDN  ! SIGE(L)*PLIJ+PTOP
!@var  PEK  PEUP**KAPA
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PEK
!@var  SQRTP  square root of P (used in diagnostics)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SQRTP
!@var  PTROPO  Pressure at mid point of tropopause level (mb)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: PTROPO
!@var  LTROPO  Tropopause layer
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: LTROPO

C**** module should own dynam variables used by other routines
!@var PTOLD pressure at beginning of dynamic time step (for clouds)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)    :: PTOLD
!@var SD_CLOUDS vert. integrated horizontal convergence (for clouds)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SD_CLOUDS
!@var DPDX_BY_RHO,DPDY_BY_RHO (pressure gradients)/density at L=1
      REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: DPDX_BY_RHO,DPDY_BY_RHO
!@var DPDX_BY_RHO_0,DPDY_BY_RHO_0 surface (pressure gradients)/density
      REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: DPDX_BY_RHO_0
      REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: DPDY_BY_RHO_0

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PU,PV
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), TARGET :: CONV
cgsfc      REAL*8, DIMENSION(IM,JM,LM-1) :: SD
      REAL*8, ALLOCATABLE :: SD(:,:,:)
!@var PIT  pressure tendency (mb m^2/s)
      REAL*8, ALLOCATABLE :: PIT(:,:)
cgsfc      EQUIVALENCE (SD(1,1,1),CONV(1,1,2))
cgsfc      EQUIVALENCE (PIT(1,1),CONV(1,1,1))

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PHI,SPA
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DUT,DVT
!!@var xAVRX scheme-depend. coefficient for AVRX: 1,byrt2 (2nd,4th order)
!      REAL*8 xAVRX

!@var PUA,PVA,SDA,PS save PU,PV,SD,P for hourly tracer advection
!@var MB Air mass array for tracers (before advection)
!@var MA Air mass array for tracers (updated during advection)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PUA,PVA,SDA,MB,MA
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: PS

!@var WCP vertical mass flux in a constant-pressure vertical
!@+   coordinate whose pressure levels are the global means of
!@+   each layer.
!@var WCPsig: WCP interpolated to terrain-following coordinate
!@+   surfaces (sigma levels)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: WCP,WCPsig

!@var DKE change in KE due to dissipation (SURF/DC/MC) (m^2/s^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DKE
!@var KEA KE on the A grid (m^2/s^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: KEA ! ke on A grid
!@var UALIJ,VALIJ U,V on the A grid (m/s)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: UALIJ,VALIJ
!@var WSAVE vertical velocity (m/s)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: WSAVE
!@var SMASS local but "SAVE"d array ADVECV in MOMEN2ND made global
!@    here since its use does not go beyond ATMDYN that calls ADVECV
      REAL*8, ALLOCATABLE:: SMASS(:)

c     INTEGER, DIMENSION(2) :: t_dyn_a, t_dyn_b, t_dyn_c, t_dyn_d
c     INTEGER, DIMENSION(2) :: t_aflux, t_advecm, t_advecv
c     INTEGER, DIMENSION(2) :: t_PGF, t_filter, t_fltruv, t_calc_pijl
c     INTEGER, DIMENSION(2) :: t_sdrag, t_calc_trop
      INTEGER, DIMENSION(2) :: t_advecv

      REAL*8, PARAMETER :: COS_LIMIT = 0.15d0

      END MODULE DYNAMICS


      SUBROUTINE ALLOC_DYNAMICS(grid)
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID
      USE RESOLUTION , ONLY : LM
      USE DYNAMICS, ONLY : PLIJ,PDSIG,AM,BYAM,PMID,PK,PEDN,PEK,
     $                     SD_CLOUDS,GZ,PU,PV,CONV,PHI,SPA,DUT,
     $                     DVT,PUA,PVA,SDA,MB,MA,DKE,WSAVE,SD,PIT,
     $                     SQRTP,PTROPO,LTROPO,PTOLD,DPDX_BY_RHO,
     $                     DPDY_BY_RHO,DPDX_BY_RHO_0,DPDY_BY_RHO_0,
     $                     PS,SMASS,KEA,UALIJ,VALIJ
c     USE DYNAMICS, only : t_dyn_a, t_dyn_b, t_dyn_c, t_dyn_d
c     USE DYNAMICS, only : t_aflux, t_advecm, t_advecv

      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER

      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO
      
      ! K-I-J arrays
      ALLOCATE ( PLIJ(LM,I_0H:I_1H,J_0H:J_1H), 
     $           PDSIG(LM,I_0H:I_1H,J_0H:J_1H),
     $             AM(LM,I_0H:I_1H,J_0H:J_1H),  
     $           BYAM(LM,I_0H:I_1H,J_0H:J_1H),
     $           PMID(LM,I_0H:I_1H,J_0H:J_1H),    
     $           PK(LM,I_0H:I_1H,J_0H:J_1H),
     $           PEDN(LM+1,I_0H:I_1H,J_0H:J_1H), 
     $            PEK(LM+1,I_0H:I_1H,J_0H:J_1H),
     $   STAT = IER)
      ! I-J-K arrays
      ALLOCATE( SD_CLOUDS(I_0H:I_1H,J_0H:J_1H,LM),  
     $                 GZ(I_0H:I_1H,J_0H:J_1H,LM), 
     $                 PU(I_0H:I_1H,J_0H:J_1H,LM),  
     $                 PV(I_0H:I_1H,J_0H:J_1H,LM), 
     $               CONV(I_0H:I_1H,J_0H:J_1H,LM), 
     $               PHI(I_0H:I_1H,J_0H:J_1H,LM), 
     $                SPA(I_0H:I_1H,J_0H:J_1H,LM), 
     $                DUT(I_0H:I_1H,J_0H:J_1H,LM), 
     $                DVT(I_0H:I_1H,J_0H:J_1H,LM), 
     $                PUA(I_0H:I_1H,J_0H:J_1H,LM), 
     $                PVA(I_0H:I_1H,J_0H:J_1H,LM), 
     $                SDA(I_0H:I_1H,J_0H:J_1H,LM),  
     $                MB(I_0H:I_1H,J_0H:J_1H,LM), 
     $                MA(I_0H:I_1H,J_0H:J_1H,LM), 
     $                DKE(I_0H:I_1H,J_0H:J_1H,LM), 
     $                KEA(I_0H:I_1H,J_0H:J_1H,LM),
     $                UALIJ(LM,I_0H:I_1H,J_0H:J_1H),
     $                VALIJ(LM,I_0H:I_1H,J_0H:J_1H),
     $              WSAVE(I_0H:I_1H,J_0H:J_1H,LM-1), 
     $              PIT(I_0H:I_1H,J_0H:J_1H),
     $              SD(I_0H:I_1H,J_0H:J_1H,LM-1),
     $   STAT = IER)

      !hack to remove NaNs from "SD" - possible mistake somewhere
      CONV(:,:,:) = 0.d0

      ! I-J arrays
      ALLOCATE(  SQRTP(I_0H:I_1H,J_0H:J_1H), 
     $          PTROPO(I_0H:I_1H,J_0H:J_1H),
     $          LTROPO(I_0H:I_1H,J_0H:J_1H),  
     $          PTOLD(I_0H:I_1H,J_0H:J_1H),
     $          DPDX_BY_RHO(I_0H:I_1H,J_0H:J_1H), 
     $          DPDY_BY_RHO(I_0H:I_1H,J_0H:J_1H),
     $          DPDX_BY_RHO_0(I_0H:I_1H,J_0H:J_1H), 
     $          DPDY_BY_RHO_0(I_0H:I_1H,J_0H:J_1H),
     $          PS(I_0H:I_1H,J_0H:J_1H),
     $   STAT = IER)

      ! J arrays
      ALLOCATE(  SMASS(J_0H:J_1H), 
     $   STAT = IER)

! correct or wrong, but being static all arrays were initialized
! to zero by default. They have to be initialized to something now 
! to avoid floating point exceptions...
      DPDX_BY_RHO(I_0H:I_1H,J_0H:J_1H) = 0.d0
      DPDY_BY_RHO(I_0H:I_1H,J_0H:J_1H) = 0.d0
      DPDX_BY_RHO_0(I_0H:I_1H,J_0H:J_1H) = 0.d0
      DPDY_BY_RHO_0(I_0H:I_1H,J_0H:J_1H) = 0.d0

      PU(I_0H:I_1H,J_0H:J_1H,1:LM) = 0.d0
      PV(I_0H:I_1H,J_0H:J_1H,1:LM) = 0.d0

      SD_CLOUDS(I_0H:I_1H,J_0H:J_1H,1:LM) = 0.d0
      PIT(:,:) = 0.d0

      END SUBROUTINE ALLOC_DYNAMICS

      subroutine perturb_temps
C**** Perturb tropospheric temperatures by at most 1 degree C
      use resolution, only : ls1,lm
      use atm_com, only : t,pk
      use random
      use domain_decomp_atm, only : grid,get
      implicit none
      integer :: I,J,L,I_0,I_1,J_0,J_1
      real*8 :: tijl,x
      integer :: nij_before_j0,nij_after_j1,nij_after_i1

      CALL GET(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1)

      CALL CALC_AMPK(LM)
      DO L=1,LS1-1
        call burn_random(nij_before_j0(J_0))
        DO J=J_0,J_1
        call burn_random((I_0-1))
        DO I=I_0,I_1
           TIJL=T(I,J,L)*PK(L,I,J)-1.+2*RANDU(X)
           T(I,J,L)=TIJL/PK(L,I,J)
        END DO
        call burn_random(nij_after_i1(I_1))
        END DO
        call burn_random(nij_after_j1(J_1))
      END DO

      return
      end subroutine perturb_temps

      SUBROUTINE CALC_VERT_AMP(P0,LMAX,PL,AM,PDSIG,PEDN,PMID)
!@sum  CALC_VERT_AMPK calculates air mass and pressure vertical arrays
!@auth Jean Lerner/Gavin Schmidt
      USE CONSTANT, only : bygrav
      USE RESOLUTION, only : ls1,ptop,psfmpt,pmtop
      USE RESOLUTION, only : lm
      USE ATM_COM, only : lm_req,req_fac,req_fac_m,req_fac_d
      USE DYNAMICS, only : dsig,sig,sige
      IMPLICIT NONE

      REAL*8, INTENT(IN) :: P0 !@var P0 surface pressure (-PTOP) (mb)
      INTEGER, INTENT(IN) :: LMAX !@var LMAX max level for calculation
!@var AM mass at each level (kg/m2)
!@var PDSIG pressure interval at each level (mb)
!@var PMID mid-point pressure (mb)
      REAL*8, INTENT(OUT), DIMENSION(LMAX) :: AM,PDSIG,PMID,PL
!@var PEDN edge pressure (top of box) (mb)
      REAL*8, INTENT(OUT), DIMENSION(LMAX+1) :: PEDN
      INTEGER :: L  !@var L  loop variables

C**** Calculate air mass, layer pressures
C**** Note that only layers LS1 and below vary as a function of surface
C**** pressure.
C**** Note Air mass is calculated in (kg/m^2)

      DO L=1,LS1-1
        PL(L)   = P0
        PDSIG(L)= P0*DSIG(L)
        PMID(L) = SIG(L)*P0+PTOP
        PEDN(L) = SIGE(L)*P0+PTOP
        AM  (L) = PDSIG(L)*1d2*BYGRAV
      END DO
      DO L=LS1,MIN(LMAX,LM)
        PL(L)   = PSFMPT
        PDSIG(L)= PSFMPT*DSIG(L)
        PMID(L) = SIG(L)*PSFMPT+PTOP
        PEDN(L) = SIGE(L)*PSFMPT+PTOP
        AM  (L) = PDSIG(L)*1d2*BYGRAV
      END DO
      IF (LMAX.ge.LM) PEDN(LM+1) = SIGE(LM+1)*PSFMPT+PTOP
C**** Rad. equ. layers if necessary (only PEDN,AM,PMID)
      IF (LMAX.eq.LM+LM_REQ) THEN
        PMID(LM+1:LM+LM_REQ) = REQ_FAC_M(1:LM_REQ)*PMTOP
          AM(LM+1:LM+LM_REQ) = REQ_FAC_D(1:LM_REQ)*PMTOP*1d2*BYGRAV
        PEDN(LM+2:LM+LM_REQ) = REQ_FAC(1:LM_REQ-1)*PEDN(LM+1)
        PEDN(LM+LM_REQ+1)=0.    ! 1d-5  ! why not zero?
      END IF

      RETURN
      END SUBROUTINE CALC_VERT_AMP

      SUBROUTINE DAILY_atmdyn(end_of_day)
!@sum  DAILY performs daily tasks at end-of-day and maybe at (re)starts
!@auth Original Development Team
!@calls constant:orbit, calc_ampk
      USE RESOLUTION, only : ls1,ptop,psf
      USE RESOLUTION, only : im,jm,lm
      USE ATM_COM, only : p
      USE MODEL_COM, only : itime,itimei
      USE GEOM, only : areag,axyp
      USE DOMAIN_DECOMP_ATM, only : grid, GET, GLOBALSUM, AM_I_ROOT
c      USE ATMDYN, only : CALC_AMPK
      IMPLICIT NONE
      REAL*8 :: DELTAP,PBAR,SMASS
      REAL*8 :: CMASS(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                grid%J_STRT_HALO:grid%J_STOP_HALO)
      INTEGER i,j,l
      LOGICAL, INTENT(IN) :: end_of_day
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, I_0,I_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

#ifndef SCM

      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      IF (.not.(end_of_day.or.itime.eq.itimei)) RETURN

C**** Tasks to be done at end of day and at initial starts only
C****
C**** THE GLOBAL MEAN PRESSURE IS KEPT CONSTANT AT PSF MILLIBARS
C****
C**** CALCULATE THE CURRENT GLOBAL MEAN PRESSURE
      DO J=J_0,J_1
      DO I=I_0,I_1
        CMASS(I,J)=P(I,J)*AXYP(I,J)
      END DO
      END DO
      if(have_south_pole) cmass(2:im,1)=cmass(1,1)
      if(have_north_pole) cmass(2:im,jm)=cmass(1,jm)
      CALL GLOBALSUM(grid, CMASS, SMASS, ALL=.TRUE.)
      PBAR=SMASS/AREAG+PTOP
C**** CORRECT PRESSURE FIELD FOR ANY LOSS OF MASS BY TRUNCATION ERROR
C****   except if it was just done (restart from itime=itimei)
      DELTAP=PSF-PBAR
      if(itime.eq.itimei .and. abs(deltap).lt.1.d-10) return
      P=P+DELTAP

      CALL CALC_AMPK(LS1-1)

      if (AM_I_ROOT()) then
         IF (ABS(DELTAP).gt.1d-6)
     *      WRITE (6,'(A25,F10.6/)') '0PRESSURE ADDED IN GMP IS',DELTAP
      end if
#endif

      RETURN
      END SUBROUTINE DAILY_atmdyn

      subroutine init_sdrag
      USE RESOLUTION, only : ls1,pstrat
      USE RESOLUTION, only : LM
      USE ATM_COM, only : pednl00,pmidl00
      use DYNAMICS, only : 
     &     USE_UNR_DRAG
     &     ,X_SDRAG,C_SDRAG,LSDRAG,P_SDRAG,LPSDRAG,PP_SDRAG,ang_sdrag
     &     ,P_CSDRAG,CSDRAGL,Wc_Jdrag,wmax,VSDRAGL
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      USE Dictionary_mod
      implicit none
      INTEGER L,LCSDRAG

      call sync_param( "X_SDRAG", X_SDRAG, 2 )
      call sync_param( "C_SDRAG", C_SDRAG )
      call sync_param( "P_CSDRAG", P_CSDRAG )
      call sync_param( "P_SDRAG", P_SDRAG )
      call sync_param( "PP_SDRAG", PP_SDRAG )
      call sync_param( "ANG_SDRAG", ANG_SDRAG )
      call sync_param( "Wc_Jdrag", Wc_Jdrag )
      call sync_param( "VSDRAGL", VSDRAGL, lm-ls1+1 )
      call sync_param( "wmax", wmax )

C**** Calculate levels for application of SDRAG: LSDRAG,LPSDRAG->LM i.e.
C**** all levels above and including P_SDRAG mb (PP_SDRAG near poles)
C**** If P is the edge between 2 levels, take the higher level.
C**** Also find CSDRAGL, the coefficients of C_Sdrag as a function of L

      LSDRAG=LM ; LPSDRAG=LM ; LCSDRAG=LM ; CSDRAGL=C_SDRAG
      DO L=1,LM
        IF (PEDNL00(L+1)-1d-5.lt.P_SDRAG .and.
     *      PEDNL00(L)  +1d-5.gt.P_SDRAG)         LSDRAG=L
        IF (PEDNL00(L+1)-1d-5.lt.PP_SDRAG .and.
     *      PEDNL00(L)  +1d-5.gt.PP_SDRAG)        LPSDRAG=L
        IF (PEDNL00(L+1)-1d-5.lt.P_CSDRAG .and.
     *      PEDNL00(L)  +1d-5.gt.P_CSDRAG)        LCSDRAG=L
      END DO
      DO L=LCSDRAG,LSDRAG-1
         CSDRAGL(L) = C_SDRAG + max( 0.d0 , (X_SDRAG(1)-C_SDRAG) *
     *     LOG(P_CSDRAG/(PMIDL00(L))) / LOG(P_CSDRAG/P_SDRAG) )
      END DO
      if (AM_I_ROOT()) then
         WRITE(6,*) "Levels for  LSDRAG =",LSDRAG ,"->",LM
         WRITE(6,*) "Levels for LPSDRAG =",LPSDRAG,"->",LM," near poles"
         WRITE(6,*) "C_SDRAG coefficients:",CSDRAGL(LS1:LSDRAG-1)
      end if

      RETURN
C****
      end subroutine init_sdrag

      subroutine read_nmc
C**** read atmospheric initial conditions file
      use resolution, only : ptop
      use resolution, only : im,jm,lm
      use atm_com, only : u,v,t,p,q
! move this back to atm_com.f after moving these arrays to atm_com
      use atm_com, only : pk,pmid,pedn,ualij,valij
      use domain_decomp_atm, only : get,grid,readt_parallel
      use filemanager
      implicit none
      integer :: I,J,L,iu_AIC,I_0,I_1,J_0,J_1
      logical :: have_south_pole,have_north_pole

      CALL GET(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      call openunit("AIC",iu_AIC,.true.,.true.)

      CALL READT_PARALLEL(grid,iu_AIC,NAMEUNIT(iu_AIC),P,1) ! Psurf
      DO J=J_0,J_1
        DO I=I_0,I_1
          P(I,J)=P(I,J)-PTOP    ! Psurf -> P
        END DO
      END DO
      CALL CALC_AMPK(LM)
      DO L=1,LM
        CALL READT_PARALLEL(grid,iu_AIC,NAMEUNIT(iu_AIC),U(:,:,L),1) ! U
      END DO
      DO L=1,LM
        CALL READT_PARALLEL(grid,iu_AIC,NAMEUNIT(iu_AIC),V(:,:,L),1) ! V
      END DO
      DO L=1,LM
        CALL READT_PARALLEL(grid,iu_AIC,NAMEUNIT(iu_AIC),T(:,:,L),1) ! Temperature
        DO J=J_0,J_1
        DO I=I_0,I_1
C**** REPLACE TEMPERATURE BY POTENTIAL TEMPERATURE
          T(I,J,L)=T(I,J,L)/PK(L,I,J)
        END DO
        END DO
      END DO
      DO L=1,LM        ! alternatively, only read in L=1,LS1 ; skip rest
        CALL READT_PARALLEL(grid,iu_AIC,NAMEUNIT(iu_AIC),Q(:,:,L),1) ! Q
      END DO

C**** Close "AIC"
      call closeunit(iu_AIC)

C**** INITIALIZE VERTICAL SLOPES OF T,Q
      call tq_zmom_init(t,q,PMID,PEDN)

#if defined(SCM) || defined(CUBED_SPHERE)
c in these cases, assume input U/V are on the A grid
      DO J=J_0,J_1
      DO I=I_0,I_1
        ualij(:,i,j) = u(i,j,:)
        valij(:,i,j) = v(i,j,:)
      END DO
      END DO
#else
c assume input U/V are on the B grid.  Need to calculate A-grid winds.
      call recalc_agrid_uv
c the latlon version of recalc_agrid_uv does not fill the poles.
c replicate polar data to avoid compiler traps in INPUT only.
      if(have_south_pole) then
        ualij(1,2:im,1) = ualij(1,1,1)
        valij(1,2:im,1) = valij(1,1,1)
      endif
      if(have_north_pole) then
        ualij(1,2:im,jm) = ualij(1,1,jm)
        valij(1,2:im,jm) = valij(1,1,jm)
      endif
#endif

      return
      end subroutine read_nmc

