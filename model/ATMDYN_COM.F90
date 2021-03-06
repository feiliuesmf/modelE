#include "rundeck_opts.h"

      Module DYNAMICS
!@sum  DYNAMICS contains all the pressure and momentum related variables
!@vers 2013/03/29
!@auth Original development team
      Use DOMAIN_DECOMP_ATM, Only: GRID
      Use RESOLUTION,        Only: LM,LS1
      Implicit  None

!**** Vertical resolution dependent variables (set in INPUT)
!@var SIGE sigma levels at layer interfaces (1)
!@var SIG,DSIG,byDSIG mid point, depth, 1/depth of sigma levels (1)
      Real*8 :: SIGE(LM+1), &  !  sige(1)=1,  sige(ls1)=0,  sige(lm+1)=-pstrat/psfmpt
                 SIG(LM),   &  ! = (sige(1:lm)+sige(2:lm+1))*0.5d0,
                DSIG(LM),   &  ! =  sige(1:lm)-sige(2:lm+1),
              byDSIG(LM)       ! =  1./DSIG

!@var MU,MV,MW,CONV (kg/s) = mass fluxes
      Real*8,Allocatable :: MU(:,:,:), MV(:,:,:), MW(:,:,:),CONV(:,:,:), &
                            PU(:,:,:), PV(:,:,:), SD(:,:,:), &
                           DUT(:,:,:),DVT(:,:,:),SPA(:,:,:)

!@var WCP = vertical mass flux in a constant-pressure vertical
!@+   coordinate whose pressure levels are the global means of each layer.
!@var WCPsig: WCP interpolated to terrain-following coordinate surfaces (sigma levels)
      Real*8,Allocatable :: WCP(:,:,:),WCPsig(:,:,:)

!@var SMASS = local but "SAVE"d array ADVECV in MOMEN2ND made global
!@    here since its use does not go beyond ATMDYN that calls ADVECV
      Real*8,Allocatable :: SMASS(:)

      Real*8,Parameter :: COS_LIMIT = 0.15d0

!@dbparam DT = (atmospheric) dynamics time step (s)
!@var NIDYN = DTsrc / DT
!@var NSTEP = number of DT steps during dynamics
!@var MRCH  = kind of step: 0 = initial forward, -1 = backward, 2 = even leap-frog, -2 = odd leap-frog
      Real*8  :: DT=450
      Integer :: NIDYN,NSTEP,MRCH

!@dbparam MFILTR: if 1 => PSL, if 2 => T, if 3 => PSL&T is filtered
!@dbparam NFILTR = DT_filter / DTsrc
!@dbparam DT_XUfilter dU is multiplied by dt/DT_XUfilter in E-W, value of 0 switches off filter
!@dbparam DT_XVfilter dV is multiplied by dt/DT_XVfilter in E-W
!@dbparam DT_YUfilter dU is multiplied by dt/DT_YUfilter in N-S
!@dbparam DT_YVfilter dV is multiplied by dt/DT_YVfilter in N-S
!@dbparam DO_POLEFIX = 1 to correct u,v tendencies near poles
!@dbparam ANG_UV     = 1 to conserve ang mom in UVfilter
!**** Controls for FLTRUV (momentum/velocity filter)
      Integer :: MFILTR=1, NFILTR=1, DO_POLEFIX=1, ANG_UV=1
      Real*8  :: DT_XUfilter=0, DT_XVfilter=0, DT_YUfilter=0, DT_YVfilter=0
!@var QUVfilter = True if any of DT_[XY][UV]filter are not 0
      Logical :: QUVfilter

!**** Stratospheric drag related parameters
!@dbparam X_SDRAG.  SDRAG ~X_SDRAG(1)+X_SDRAG(2)*wind_magnitude
!@dbparam C_SDRAG.  SDRAG=C_SDRAG (const.) above PTOP
!@dbparam P_CSDRAG pressure level above which const.drag is increased
!@dbparam P_SDRAG = pressure level above which SDRAG is applied (mb), PP_SDRAG = near poles
!@dbparam Wc_JDRAG = critical velocity for J.Hansen/Judith Perlwitz drag; if 0 no JDRAG feature in SDRAG
!@dbparam WMAX = imposed limit for stratospheric winds (m/s) in SDRAG
!@dbparam VSDRAGL = tuning factor for stratospheric drag (not =1 e.g. if used with explicit grav.wave drag scheme)
!@dbparam USE_UNR_DRAG: if 1 => SDRAG is turned off and GWD is applied
!@+                     if 0 => SDRAG is kept intact and alternative GWD is not employed
      Real*8  :: X_SDRAG(2) = (/2.5d-4,2.5d-5/), C_SDRAG = 2.5d-5, &
                 CSDRAGL(LS1:LM), &
                 P_CSDRAG=0, P_SDRAG=0, PP_SDRAG=1, &
                 Wc_JDRAG=30, WMAX=200, VSDRAGL(LS1:LM)=1
      Integer :: USE_UNR_DRAG=0
!@var LSDRAG = level above which SDRAG is applied, LPSDRAG = near pole
!@var ANG_SDRAG if =1: angular momentum lost by SDRAG is added in below PTOP
      Integer :: LSDRAG=LM, LPSDRAG=LM, ANG_SDRAG=1


!@var linear_sdrag flag whether to use a condition-independent
!@+   rayleigh friction timescale (rtau) rather than the default
!@+   inverse timescale of cdn*|u|/deltaz, where
!@+   cdn is a drag coefficient that may depend on wind speed
!@+       and other factors
!@+   |u| is the wind speed
!@+   deltaz is local layer thickness
      logical :: linear_sdrag
!@var l1_rtau lowest layer at which to apply linear_sdrag
      integer :: l1_rtau
!@dbparam rtau rayleigh friction timescale (seconds) as a function of layer
!@+            (only used with linear sdrag scheme)
      real*8, allocatable :: rtau(:) ! allocated over l1_rtau:lm

!@var mincolmass, maxcolmass minimum/maximum allowed column mass (kg/m2)
      real*8 :: mincolmass, maxcolmass

!**** Variables specific for stratosphere and/or strat diagnostics
!@var DO_GWDRAG when true, prints Gravity Wave diagnostics
!@var iDO_GWDRAG number if AIJ Gravity wave diagnostics
      Logical :: DO_GWDRAG = .False.
      Integer :: iDO_GWDRAG = 0

      EndModule DYNAMICS


!!!#ifndef SCM
      Subroutine ALLOC_DYNAMICS (GRID)
      Use DOMAIN_DECOMP_ATM, Only: DIST_GRID, AM_I_ROOT
      Use RESOLUTION, Only: LM,LS1, PSFMPT,PLBOT,PTOP
      Use DYNAMICS, Only: SIGE,SIG,DSIG,BYDSIG, MU,MV,MW,CONV, PU,PV,SD, DUT,DVT,SPA, SMASS,WCP,WCPsig
      Use ATM_COM,  Only: LM_REQ, PL00,PMIDL00,PDSIGL00,AML00,BYAML00,PEDNL00
      Implicit  None
      TYPE (DIST_GRID), Intent(In) :: GRID
      Integer :: I1H,INH,J1H,JNH, LMR,IER

      I1H = GRID%I_STRT_HALO  ;  INH = GRID%I_STOP_HALO  ;  J1H = GRID%J_STRT_HALO  ;  JNH = GRID%J_STOP_HALO

!****
!**** Set dependent vertical resolution variables
!****
      SIGE(:) = (PLBOT(:)-PTOP)/PSFMPT
      SIG(:)  = (sige(1:lm)+sige(2:lm+1))*0.5d0
      DSIG(:) =  sige(1:lm)-sige(2:lm+1)
    byDSIG(:) =  1 / DSIG(:)
!**** Check the vertical layering defined in RES_ (is sige(ls1)=0 ?)
      IF (SIGE(LS1).ne.0.) then
         If (AM_I_ROOT())  Write (6,*) 'bad vertical layering: ls1,sige(ls1)',ls1,sige(ls1)
         call stop_model('INPUT: ls1 incorrectly set in RES_',255)  ;  END IF
!**** Calculate default vertical arrays (including rad. eq. layers)
      LMR = LM + LM_REQ
      Call CALC_VERT_AMP (PSFMPT,LMR,PL00,AML00,PDSIGL00,PEDNL00,PMIDL00)
      BYAML00(:) = 1 / AML00(:)

      Allocate (MU(I1H:INH,J1H:JNH,LM),  MV(I1H:INH,J1H:JNH,LM),  MW(I1H:INH,J1H:JNH,LM-1), &
              CONV(I1H:INH,J1H:JNH,LM), &
                PU(I1H:INH,J1H:JNH,LM),  PV(I1H:INH,J1H:JNH,LM),  SD(I1H:INH,J1H:JNH,LM-1), &
               DUT(I1H:INH,J1H:JNH,LM), DVT(I1H:INH,J1H:JNH,LM), SPA(I1H:INH,J1H:JNH,LM), &
            WCPsig(I1H:INH,J1H:JNH,LM), WCP(I1H:INH,J1H:JNH,LM), &
                     SMASS(J1H:JNH),  Stat=IER)

! correct or wrong, but being static all arrays were initialized
! to zero by default. They have to be initialized to something now
! to avoid floating point exceptions...
      MU(:,:,:) = 0  ;  MV(:,:,:) = 0  ;  CONV(:,:,:) = 0
      PU(:,:,:) = 0  ;  PV(:,:,:) = 0
      EndSubroutine ALLOC_DYNAMICS
!!!#endif


      Subroutine MAtoP (MA)                                        
!@sum MAtoP calculates haloed pressure arrays PEDN, PMID, PDSIG and PK from haloed air mass MA
      Use CONSTANT,   Only: kg2mb,KAPA
      Use RESOLUTION, Only: LM, MTOP
      Use ATM_COM,    Only: PEDN,PMID,PDSIG,PK
      Use DOMAIN_DECOMP_ATM, Only: GRID
      Use DOMAIN_DECOMP_1D,  Only: GetDomainBounds, HALO_UPDATE_COLUMN, SOUTH
      Implicit  None
      Real*8,Intent(In) :: MA(LM, GRID%I_STRT_HALO:GRID%I_STOP_HALO, GRID%J_STRT_HALO:GRID%J_STOP_HALO)
      Real*8  :: M
      Integer :: I,J,L, I1,IN,J1,JN

#ifndef CUBED_SPHERE                                   !  Lat-Lon Grid
      I1 = GRID%I_STRT_HALO  ;  IN = GRID%I_STOP_HALO  !  1:IM
      J1 = GRID%J_STRT_HALO  ;  JN = GRID%J_STOP       !  primary row limits used by velocity
#endif

#ifdef CUBED_SPHERE                                    !  Cube-Sphere grid
      I1 = GRID%I_STRT_HALO  ;  IN = GRID%I_STOP_HALO  !  haloed primary column limits
      J1 = GRID%J_STRT_HALO  ;  JN = GRID%J_STOP_HALO  !  haloed primary row limits
#endif

      Do J=J1,JN  ;  Do I=I1,IN
         M = MTOP
         Do L=LM,1,-1
            PEDN(L,I,J) = kg2mb * (M + MA(L,I,J))
            PMID(L,I,J) = kg2mb * (M + MA(L,I,J)*.5)
           PDSIG(L,I,J) = kg2mb * MA(L,I,J)
              PK(L,I,J) = PMID(L,I,J)**KAPA
            M = M + MA(L,I,J)  ;  EndDo  ;  EndDo  ;  EndDo

      Return
      EndSubroutine MAtoP


      Subroutine CALC_VERT_AMP (P0,LMAX,PL,MA,PDSIG,PEDN,PMID)
!@sum  CALC_VERT_AMPK calculates air mass and pressure vertical arrays
!@auth Jean Lerner/Gavin Schmidt
      Use CONSTANT,   Only: byGRAV
      Use RESOLUTION, Only: LM,LS1, PTOP,PSFMPT,PMTOP
      Use ATM_COM,    Only: LM_REQ, REQ_FAC,REQ_FAC_M,REQ_FAC_D
      Use DYNAMICS,   Only: dsig,sig,sige
      Implicit  None

!@var LMAX = max level for calculation
!@var P0   = surface pressure - PTOP (mb)
!@var MA mass per unit area for each layer (kg/m^2)
!@var PDSIG pressure interval at each level (mb)
!@var PMID mid-point pressure (mb)
!@var PEDN edge pressure (top of box) (mb)
      Integer,Intent(In)  :: LMAX
      Real*8, Intent(In)  :: P0
      Real*8, Intent(Out) :: MA(LMAX),PDSIG(LMAX),PMID(LMAX),PL(LMAX),PEDN(LMAX+1)
      Integer :: L

!**** Calculate air mass, layer pressures
!**** Note that only layers LS1 and below vary as a function of surface
!**** pressure.
      Do L=1,LS1-1
           PL(L) = P0
        PDSIG(L) = P0*DSIG(L)
         PMID(L) = SIG(L)*P0 + PTOP
         PEDN(L) = SIGE(L)*P0 + PTOP
           MA(L) = PDSIG(L)*1d2*BYGRAV  ;  EndDo
      Do L=LS1,Min(LMAX,LM)
           PL(L) = PSFMPT
        PDSIG(L) = PSFMPT*DSIG(L)
         PMID(L) = SIG(L)*PSFMPT + PTOP
         PEDN(L) = SIGE(L)*PSFMPT + PTOP
           MA(L) = PDSIG(L)*1d2*BYGRAV  ;  EndDo
      If (LMAX >= LM)  PEDN(LM+1) = SIGE(LM+1)*PSFMPT + PTOP
!*** Radiation equilibrium layers if necessary
      If (LMAX == LM+LM_REQ)  Then
         PMID(LM+1:LM+LM_REQ) = REQ_FAC_M(1:LM_REQ)*PMTOP
           MA(LM+1:LM+LM_REQ) = REQ_FAC_D(1:LM_REQ)*PMTOP*1d2*BYGRAV
         PEDN(LM+2:LM+LM_REQ) = REQ_FAC(1:LM_REQ-1)*PEDN(LM+1)
         PEDN(LM+LM_REQ+1) = 0  ;  EndIf

      Return
      EndSubroutine CALC_VERT_AMP


      Subroutine READ_NMC
!**** read atmospheric initial conditions file
      Use CONSTANT,   Only: mb2kg,areag,rgas
      Use RESOLUTION, Only: IM,JM,LM, MTOP,MFIX,MFIXs,MFRAC, PSF,PTOP
      Use ATM_COM,    Only: MA,U,V,T,P,Q, PK,PMID,PEDN,UALIJ,VALIJ, ZATMO
      Use DOMAIN_DECOMP_ATM, Only: GRID, GetDomainBounds, globalsum
      use pario, only : par_open,par_close,read_dist_data
      use GEOM, only : axyp
      use Dictionary_mod
      Implicit none
      Integer :: I,J,L,fid, I1,IN,J1,JN
      Logical :: QSP,QNP
      Real*8  :: MVAR
      integer :: initial_psurf_from_topo
      real*8, dimension(:,:), allocatable :: expz,aexpz
      real*8 :: aexpz_sum

      Call GetDomainBounds (GRID, I_STRT=I1, I_STOP=IN, J_STRT=J1, J_STOP=JN, &
                                  HAVE_SOUTH_POLE=QSP, HAVE_NORTH_POLE=QNP)


      fid = par_open(grid,'AIC','read')
      call read_dist_data(grid,fid,'p',p)
      call read_dist_data(grid,fid,'u',u)
      call read_dist_data(grid,fid,'v',v)
      call read_dist_data(grid,fid,'t',t)
      call read_dist_data(grid,fid,'q',q)
      call par_close(grid,fid)

      if(is_set_param('initial_psurf_from_topo')) &
           call get_param('initial_psurf_from_topo',initial_psurf_from_topo)
      if(initial_psurf_from_topo==1) then
!**** Reset initial surface pressure to be approximately hydrostatically consistent
!**** with the orography.  Regional lapse rates are not taken into account (yet),
!**** as this degree of precision is likely not necessary for the cold-start
!**** scenarios for which this option was created.
        allocate(expz(grid%i_strt_halo:grid%i_stop_halo, &
                      grid%j_strt_halo:grid%j_stop_halo), &
                aexpz(grid%i_strt_halo:grid%i_stop_halo, &
                      grid%j_strt_halo:grid%j_stop_halo) )
        do J=J1,JN
        do I=I1,IN
          ! note zatmo is actually gravity times surface elevation
          expz(i,j) = exp(-zatmo(i,j)/(rgas*t(i,j,1)))
          aexpz(i,j) = axyp(i,j)*expz(i,j)
        enddo
        enddo
        call globalsum(grid,aexpz,aexpz_sum,all=.true.)
        do J=J1,JN
        do I=I1,IN
          p(i,j) = psf*expz(i,j)/(aexpz_sum/areag)
        enddo
        enddo
        deallocate(expz,aexpz)
      endif

      Do J=J1,JN
      Do I=I1,IN
         MVAR = P(I,J)*mb2kg - MFIXs - MTOP
         MA(:,I,J) = MFIX(:) + MVAR*MFRAC(:)
         P(I,J) = P(I,J) - PTOP  !  Psurf -> P
      EndDo
      EndDo
      Call CALC_AMPK (LM)

!**** Convert Temperature to Potential Temperature
      Do L=1,LM
        T(I1:IN,J1:JN,L) = T(I1:IN,J1:JN,L) / PK(L,I1:IN,J1:JN)
      EndDo

!**** INITIALIZE VERTICAL SLOPES OF T,Q
      Call TQ_ZMOM_INIT (T,Q,PMID,PEDN)

#if defined(SCM) || defined(CUBED_SPHERE)
! in these cases, assume input U/V are on the A grid
      Do J=J1,JN  ;  Do I=I1,IN
         UALIJ(:,I,J) = U(I,J,:)
         VALIJ(:,I,J) = V(I,J,:)  ;  EndDo  ;  EndDo
#else
! assume input U/V are on the B grid.  Need to calculate A-grid winds.
      Call RECALC_AGRID_UV
! the latlon version of recalc_agrid_uv does not fill the poles.
! replicate polar data to avoid compiler traps in INPUT only.
      If (QSP)  Then
         UALIJ(1,2:IM,1) = UALIJ(1,1,1)
         VALIJ(1,2:IM,1) = VALIJ(1,1,1)  ;  EndIf
      If (QNP)  Then
         UALIJ(1,2:IM,JM) = UALIJ(1,1,JM)
         VALIJ(1,2:IM,JM) = VALIJ(1,1,JM)  ;  EndIf
#endif

      Return
      EndSubroutine READ_NMC


      Subroutine PERTURB_TEMPS
!**** Perturb tropospheric temperatures by at most 1 degree C
      Use RESOLUTION, Only: LM,LS1
      Use ATM_COM,    Only: T,PK
      Use RANDOM
      Use domain_decomp_atm, only : grid,getDomainBounds
      Implicit None
      Integer :: I,J,L, I1,IN,J1,JN
      Real*8  :: TIJL,X
      Integer :: nij_before_j0,nij_after_j1,nij_after_i1

      Call GetDomainBounds (GRID, I_STRT=I1, I_STOP=IN, J_STRT=J1, J_STOP=JN)

      Call CALC_AMPK (LM)
      Do L=1,LS1-1
         Call BURN_RANDOM (nij_before_j0(J1))
         Do J=J1,JN
            Call BURN_RANDOM ((I1-1))
            Do I=I1,IN
               TIJL = T(I,J,L)*PK(L,I,J) - 1 + 2*RANDU(X)
               T(I,J,L) = TIJL/PK(L,I,J)  ;  EndDo
            Call BURN_RANDOM (nij_after_i1(IN))  ;  EndDo
         Call BURN_RANDOM (nij_after_j1(JN))  ;  EndDo

      Return
      EndSubroutine PERTURB_TEMPS


      Subroutine INIT_SDRAG
      Use RESOLUTION, Only: LM,LS1, PSTRAT
      Use ATM_COM,    Only: PEDNL00,PMIDL00
      Use DYNAMICS,   Only: LSDRAG,LPSDRAG,ANG_SDRAG,USE_UNR_DRAG, &
                            X_SDRAG,C_SDRAG,P_SDRAG,PP_SDRAG,P_CSDRAG,CSDRAGL,Wc_JDRAG,WMAX,VSDRAGL
      use dynamics, only : l1_rtau,rtau,linear_sdrag
      Use DOMAIN_DECOMP_ATM, Only: AM_I_ROOT
      Use Dictionary_mod
      Implicit None
      Integer :: L,LCSDRAG,nrtau
      character(len=1) :: partype

      linear_sdrag = is_set_param('rtau')

      if(linear_sdrag) then

        call query_param('rtau',nrtau,partype)
        l1_rtau = 1 + lm - nrtau
        allocate(rtau(l1_rtau:lm))
        call get_param('rtau',rtau,nrtau)

      else

        Call sync_param ("X_SDRAG",  X_SDRAG, 2 )
        Call sync_param ("C_SDRAG",  C_SDRAG )
        Call sync_param ("P_CSDRAG", P_CSDRAG )
        Call sync_param ("P_SDRAG",  P_SDRAG )
        Call sync_param ("PP_SDRAG", PP_SDRAG )
        Call sync_param ("ANG_SDRAG",ANG_SDRAG )
        Call sync_param ("Wc_Jdrag", Wc_Jdrag )
        Call sync_param ("VSDRAGL",  VSDRAGL, LM-LS1+1 )
        Call sync_param ("wmax",     WMAX )

!**** Calculate levels for application of SDRAG: LSDRAG,LPSDRAG->LM i.e.
!**** all levels above and including P_SDRAG mb (PP_SDRAG near poles)
!**** If P is the edge between 2 levels, take the higher level.
!**** Also find CSDRAGL, the coefficients of C_Sdrag as a function of L

        LSDRAG=LM ; LPSDRAG=LM ; LCSDRAG=LM ; CSDRAGL=C_SDRAG
        DO L=1,LM
         If (PEDNL00(L+1)-1d-5 <  P_SDRAG .and. PEDNL00(L)+1d-5 >  P_SDRAG)  LSDRAG  = L
         If (PEDNL00(L+1)-1d-5 < PP_SDRAG .and. PEDNL00(L)+1d-5 > PP_SDRAG)  LPSDRAG = L
         If (PEDNL00(L+1)-1d-5 < P_CSDRAG .and. PEDNL00(L)+1d-5 > P_CSDRAG)  LCSDRAG = L
        EndDo
        DO L=LCSDRAG,LSDRAG-1
         CSDRAGL(L) = C_SDRAG + Max(0d0, (X_SDRAG(1)-C_SDRAG)*Log(P_CSDRAG/(PMIDL00(L))) / Log(P_CSDRAG/P_SDRAG))
        EndDo
        If (AM_I_ROOT()) then
         Write (6,*) "Levels for  LSDRAG =",LSDRAG ,"->",LM
         Write (6,*) "Levels for LPSDRAG =",LPSDRAG,"->",LM," near poles"
         Write (6,*) "C_SDRAG coefficients:",CSDRAGL(LS1:LSDRAG-1)
        EndIf

      endif  ! linear_drag or not

      Return
      EndSubroutine INIT_SDRAG


      Subroutine DAILY_ATMDYN (end_of_day)
!@sum  DAILY performs daily tasks at end-of-day and maybe at (re)starts
!@auth Original Development Team
!@calls constant:orbit, calc_ampk
      Use RESOLUTION, Only: IM,JM,LM,LS1, PTOP,PSF
      Use ATM_COM,    Only: P
      Use MODEL_COM,  Only: ITIME,ITIMEI
      Use GEOM,       Only: AREAG,AXYP
      Use DOMAIN_DECOMP_ATM, Only: GRID, GetDomainBounds, GLOBALSUM, AM_I_ROOT
!      USE ATMDYN, only : CALC_AMPK
      Implicit None
      Logical,Intent(In) :: END_of_DAY
      Real*8  :: DELTAP,PBAR,SMASS, CMASS(GRID%I_STRT_HALO:GRID%I_STOP_HALO,GRID%J_STRT_HALO:GRID%J_STOP_HALO)
      Integer :: I,J,L, I1,IN,J1,JN
      Logical :: QSP,QNP

#ifndef SCM
      If (.not.(END_of_DAY .or. ITIME==ITIMEI))  Return
      Call GetDomainBounds (GRID, I_STRT=I1, I_STOP=IN, J_STRT=J1, J_STOP=JN, &
                                  HAVE_SOUTH_POLE=QSP, HAVE_NORTH_POLE=QNP)

!**** Tasks to be done at end of day and at initial starts only
!****
!**** THE GLOBAL MEAN PRESSURE IS KEPT CONSTANT AT PSF MILLIBARS
!****
!**** CALCULATE THE CURRENT GLOBAL MEAN PRESSURE
      CMASS(I1:IN,J1:JN) = P(I1:IN,J1:JN) * AXYP(I1:IN,J1:JN)
      If (QSP)  CMASS(2:IM,1)  = CMASS(1,1)
      If (QNP)  CMASS(2:IM,JM) = CMASS(1,JM)
      Call GLOBALSUM (GRID, CMASS, SMASS, ALL=.TRUE.)
      PBAR = SMASS/AREAG + PTOP
!**** CORRECT PRESSURE FIELD FOR ANY LOSS OF MASS BY TRUNCATION ERROR
!****   except if it was just done (restart from itime=itimei)
      DELTAP = PSF-PBAR
      If (ITIME==ITIMEI .and. Abs(DELTAP) < 1d-10)  Return
      P(:,:) = P(:,:) + DELTAP
      Call CALC_AMPK (LS1-1)

      If (AM_I_ROOT() .and. Abs(DELTAP) > 1d-6) &
         Write (6,'(A25,F10.6/)') '0PRESSURE ADDED IN GMP IS',DELTAP
#endif
      Return
      EndSubroutine DAILY_ATMDYN
