#include "rundeck_opts.h"

      MODULE DYNAMICS
!@sum  DYNAMICS contains all the pressure and momentum related variables
!@auth Original development team
      USE DOMAIN_DECOMP_ATM, ONLY : grid
      USE RESOLUTION , ONLY : LS1,LM
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

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PU,PV
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:), TARGET :: CONV
      REAL*8, ALLOCATABLE :: SD(:,:,:)
!@var PIT  pressure tendency (mb m^2/s)
      REAL*8, ALLOCATABLE :: PIT(:,:)

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DUT,DVT,SPA

!@var WCP vertical mass flux in a constant-pressure vertical
!@+   coordinate whose pressure levels are the global means of
!@+   each layer.
!@var WCPsig: WCP interpolated to terrain-following coordinate
!@+   surfaces (sigma levels)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: WCP,WCPsig

!@var SMASS local but "SAVE"d array ADVECV in MOMEN2ND made global
!@    here since its use does not go beyond ATMDYN that calls ADVECV
      REAL*8, ALLOCATABLE:: SMASS(:)

      REAL*8, PARAMETER :: COS_LIMIT = 0.15d0

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

      END MODULE DYNAMICS


!#ifndef SCM
      SUBROUTINE ALLOC_DYNAMICS(grid)
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID,AM_I_ROOT
      USE RESOLUTION , ONLY : LM
      USE DYNAMICS, ONLY : SIGE,SIG,DSIG,BYDSIG
      USE DYNAMICS, ONLY :
     $                     PU,PV,CONV,SPA,DUT,
     $                     DVT,SD,PIT,
     $                     SMASS,WCP,WCPsig
      USE ATM_COM, only : lm_req
      USE ATM_COM, only : PL00,PMIDL00,PDSIGL00,AML00,BYAML00,PEDNL00
      USE RESOLUTION, only : ls1,psfmpt,plbot,ptop
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: LMR,IER

      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO
      
C****
C**** Set dependent vertical resolution variables
C****
      SIGE(:) = (PLbot(:)-PTOP)/PSFMPT
      SIG(:)  = (sige(1:lm)+sige(2:lm+1))*0.5d0
      DSIG(:) =  sige(1:lm)-sige(2:lm+1)
      byDSIG  =  1./DSIG
C**** Check the vertical layering defined in RES_ (is sige(ls1)=0 ?)
      IF (SIGE(LS1).ne.0.) then
        if (AM_I_ROOT())
     *       write(6,*) 'bad vertical layering: ls1,sige(ls1)',
     &       ls1,sige(ls1)
        call stop_model('INPUT: ls1 incorrectly set in RES_',255)
      END IF
C**** Calculate default vertical arrays (including rad. eq. layers)
      LMR=LM+LM_REQ
      CALL CALC_VERT_AMP(PSFMPT,LMR,PL00,AML00,PDSIGL00,PEDNL00,PMIDL00)
      BYAML00(:)=1./AML00(:)

      ! I-J-K arrays
      ALLOCATE(
     $                 PU(I_0H:I_1H,J_0H:J_1H,LM),  
     $                 PV(I_0H:I_1H,J_0H:J_1H,LM), 
     $               CONV(I_0H:I_1H,J_0H:J_1H,LM), 
     $                SPA(I_0H:I_1H,J_0H:J_1H,LM), 
     $                DUT(I_0H:I_1H,J_0H:J_1H,LM), 
     $                DVT(I_0H:I_1H,J_0H:J_1H,LM), 
     $              PIT(I_0H:I_1H,J_0H:J_1H),
     $              SD(I_0H:I_1H,J_0H:J_1H,LM-1),
     $   STAT = IER)

      !hack to remove NaNs from "SD" - possible mistake somewhere
      CONV(:,:,:) = 0.d0

      ! I-J arrays
      ALLOCATE(
     $          WCP(I_0H:I_1H,J_0H:J_1H,LM),
     $          WCPsig(I_0H:I_1H,J_0H:J_1H,LM),
     $   STAT = IER)

      ! J arrays
      ALLOCATE(  SMASS(J_0H:J_1H), 
     $   STAT = IER)

! correct or wrong, but being static all arrays were initialized
! to zero by default. They have to be initialized to something now 
! to avoid floating point exceptions...

      PU(I_0H:I_1H,J_0H:J_1H,1:LM) = 0.d0
      PV(I_0H:I_1H,J_0H:J_1H,1:LM) = 0.d0

      PIT(:,:) = 0.d0

      END SUBROUTINE ALLOC_DYNAMICS
!#endif

      SUBROUTINE CALC_VERT_AMP (P0,LMAX,PL,MA,PDSIG,PEDN,PMID)
!@sum  CALC_VERT_AMPK calculates air mass and pressure vertical arrays
!@vers 2013/03/26
!@auth Jean Lerner/Gavin Schmidt
      USE CONSTANT, only : bygrav
      USE RESOLUTION, only : ls1,ptop,psfmpt,pmtop
      USE RESOLUTION, only : lm
      USE ATM_COM, only : lm_req,req_fac,req_fac_m,req_fac_d
      USE DYNAMICS, only : dsig,sig,sige
      IMPLICIT NONE

      REAL*8, INTENT(IN) :: P0 !@var P0 surface pressure (-PTOP) (mb)
      INTEGER, INTENT(IN) :: LMAX !@var LMAX max level for calculation
!@var MA mass per unit area for each layer (kg/m^2)
!@var PDSIG pressure interval at each level (mb)
!@var PMID mid-point pressure (mb)
      REAL*8, INTENT(OUT), DIMENSION(LMAX) :: MA,PDSIG,PMID,PL
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
        MA  (L) = PDSIG(L)*1d2*BYGRAV
      END DO
      DO L=LS1,MIN(LMAX,LM)
        PL(L)   = PSFMPT
        PDSIG(L)= PSFMPT*DSIG(L)
        PMID(L) = SIG(L)*PSFMPT+PTOP
        PEDN(L) = SIGE(L)*PSFMPT+PTOP
        MA  (L) = PDSIG(L)*1d2*BYGRAV
      END DO
      IF (LMAX.ge.LM) PEDN(LM+1) = SIGE(LM+1)*PSFMPT+PTOP
C**** Rad. equ. layers if necessary (only PEDN,MA,PMID)
      IF (LMAX.eq.LM+LM_REQ) THEN
        PMID(LM+1:LM+LM_REQ) = REQ_FAC_M(1:LM_REQ)*PMTOP
          MA(LM+1:LM+LM_REQ) = REQ_FAC_D(1:LM_REQ)*PMTOP*1d2*BYGRAV
        PEDN(LM+2:LM+LM_REQ) = REQ_FAC(1:LM_REQ-1)*PEDN(LM+1)
        PEDN(LM+LM_REQ+1)=0.    ! 1d-5  ! why not zero?
      END IF

      RETURN
      END SUBROUTINE CALC_VERT_AMP

      subroutine read_nmc
C**** read atmospheric initial conditions file
      use resolution, only : ptop
      use resolution, only : im,jm,lm
      use atm_com, only : u,v,t,p,q
! move this back to atm_com.f after moving these arrays to atm_com
      use atm_com, only : pk,pmid,pedn,ualij,valij
      use domain_decomp_atm, only : getDomainBounds,grid,readt_parallel
      use filemanager
      implicit none
      integer :: I,J,L,iu_AIC,I_0,I_1,J_0,J_1
      logical :: have_south_pole,have_north_pole

      call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1, 
     &               J_STRT=J_0, J_STOP=J_1,
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

      subroutine perturb_temps
C**** Perturb tropospheric temperatures by at most 1 degree C
      use resolution, only : ls1,lm
      use atm_com, only : t,pk
      use random
      use domain_decomp_atm, only : grid,getDomainBounds
      implicit none
      integer :: I,J,L,I_0,I_1,J_0,J_1
      real*8 :: tijl,x
      integer :: nij_before_j0,nij_after_j1,nij_after_i1

      call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1, 
     &     J_STRT=J_0, J_STOP=J_1)

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

      SUBROUTINE DAILY_atmdyn(end_of_day)
!@sum  DAILY performs daily tasks at end-of-day and maybe at (re)starts
!@auth Original Development Team
!@calls constant:orbit, calc_ampk
      USE RESOLUTION, only : ls1,ptop,psf
      USE RESOLUTION, only : im,jm,lm
      USE ATM_COM, only : p
      USE MODEL_COM, only : itime,itimei
      USE GEOM, only : areag,axyp
      USE DOMAIN_DECOMP_ATM, only : grid, getDomainBounds
      USE DOMAIN_DECOMP_ATM, only : GLOBALSUM, AM_I_ROOT
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

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
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
