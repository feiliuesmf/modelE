#define VERIFY_(rc) If (rc /= ESMF_SUCCESS) Call abort_core(__LINE__,rc)

module FV_CS_Mod

  use ESMF_MOD
  use FV_UTILS

  implicit none
  private

  public :: createInternalRestart
  public :: gridCompInit
  public :: setupForESMF
  public :: addQfieldToFVImport
  public :: hydrostatic
  public :: accumulate_mass_fluxes

  logical, parameter :: hydrostatic = .true.

  ! private data of convenience
  Integer :: rc ! return code from ESMF
  integer, parameter :: prec = kind(0.0D0) ! double precision

  contains

!-----------------------------------------------------------------------------
  subroutine setupForESMF(fv, vm, grid, cf, config_file)
!-------------------------------------------------------------------------------
    use DOMAIN_DECOMP_ATM, only : AM_I_ROOT
#ifdef CUBED_SPHERE
    use fvdycorecubed_gridcomp, only: SetServices
#else
    use fvdycore_gridcompmod, only: SetServices
#endif

    Type(FV_CORE), intent(inout) :: fv
    type (ESMF_vm),   intent(in) :: vm
    type (ESMF_grid), intent(inout) :: grid      ! fv grid
    type (ESMF_config), intent(inout) :: cf      ! config object
    character(len=*),  intent(in) :: config_file ! filename for resource file

! CC : Changing CUBED_SPHERE gridCompName causes the model to crash. Why?
! e.g. character(len=*), parameter :: gridCompName = 'FVCS-DYCORE'
#ifdef CUBED_SPHERE
    character(len=*), parameter :: gridCompName = 'FVCORE'
#else
    character(len=*), parameter :: gridCompName = 'FV dynamics'
#endif

    fv%vm  = vm
    fv%grid = grid

    ! Load configuration information from resource file
    !  ------------------------------------------------
    cf = load_configuration(config_file)

    !  Create the dynamics component, using same layout as application
    !  ------------------------------------------------------------------
    fv%gc = ESMF_GridCompCreate ( name=gridCompName, &
         & grid=grid, gridcomptype=ESMF_ATM,         &
         & config=cf, rc=rc)
    VERIFY_(rc)

    ! Create couplings
    fv%import = ESMF_StateCreate('fv dycore imports', ESMF_STATE_IMPORT, rc=rc)
    VERIFY_(rc)
    fv%export = ESMF_StateCreate('fv dycore exports', ESMF_STATE_EXPORT, rc=rc)
    VERIFY_(rc)

    !  Register services for components
    !  --------------------------------
    call ESMF_GridCompSetServices ( fv%gc, SetServices, rc )
    VERIFY_(rc)

    ! Create layout resource file - independent of actual restart file.
    If (AM_I_ROOT()) Then
       Call Write_Layout(FVCORE_LAYOUT, fv, hydrostatic)
    End If

  contains
      
    !---------------------------------------------------------------------------
    Subroutine write_layout(fname, fv, hydrostatic)
    !---------------------------------------------------------------------------
      use dist_grid_mod, only: grid, axisIndex, getAxisIndex
#ifdef CUBED_SPHERE
      Use MAPL_IOMod, only: GETFILE, Free_file
#else
      Use GEOS_IOMod, only: GETFILE, Free_file
#endif
      USE RESOLUTION, only: IM, JM, LM, LS1
      Use MODEL_COM,  only: DT=>DTsrc
      character(len=*), intent(in) :: fname
      type (fv_core)  , intent(in) :: fv
      logical         , intent(in) :: hydrostatic

      Type (axisIndex), pointer :: AI(:,:)
      Integer :: unit
      Integer :: npes, my_pet
      Integer :: mppnx, mppny

      call ESMF_VMget(fv%vm, petcount=npes, rc=rc)
      VERIFY_(rc)
      Allocate(AI(npes,3))
      call getAxisIndex(grid, AI)
      VERIFY_(rc)

      unit = GetFile(fname, form="formatted", rc=rc)
      VERIFY_(rc)
    
      mppnx = 0
      mppny = 0
      if(mod(NPES,6) == 0) then
         mppnx = int(floor(sqrt(real(NPES/6))))
         mppny = (NPES / mppnx) / 6
      endif

      ! Lat-Lon parameters
      write(unit,*)' # empty line #'
      Write(unit,*)'xy_yz_decomp:',1,npes,npes,1
      Write(unit,*)'    im: ',IM
      Write(unit,*)'    jm: ',JM
      Write(unit,*)'    km: ',LM
      Write(unit,*)'    dt: ',DT
      Write(unit,*)' ntotq: ',1
      Write(unit,*)'  imxy: ',IM
      Write(unit,'(a,100(1x,I3))')'  jmxy: ',AI(:,2)%MAX-AI(:,2)%MIN+1
      Write(unit,'(a,100(1x,I3))')'  jmyz: ',AI(:,2)%MAX-AI(:,2)%MIN+1
      Write(unit,*)'   kmyz: ',LM

      ! Common parameters
      Write(unit,*)'nsplit: ',0
      Write(unit,*)'    nq: ',1

      ! Cubed-Sphere parameters
      Write(unit,*)'      npx: ',IM
      Write(unit,*)'      npy: ',JM
      Write(unit,*)'      npz: ',LM
      Write(unit,*)'       dt: ',DT
      Write(unit,*)'   npes_x: ',mppnx
      Write(unit,*)'   npes_y: ',mppny
      Write(unit,*)' fv_debug: ',.false.
      Write(unit,*)' hydrostatic: ',.true.
      Write(unit,*)' Make_NH: ',.false.

      Call Free_File(unit)

      Deallocate(AI)

    End Subroutine write_layout

  end subroutine setupForESMF

!-------------------------------------------------------------------------------
  subroutine gridCompInit(fv, clock)
!-------------------------------------------------------------------------------
    use FV_StateMod, only:FV_RESET_CONSTANTS
    use CONSTANT, only: pi, omega, sha, radius, rvap, grav, lhe, rgas, kapa

    Type(FV_CORE), intent(inout) :: fv
    type (ESMF_clock), intent(inout) :: clock
    logical :: z_tracer

    call ESMF_GridCompInitialize ( fv%gc, importState=fv%import, &
         & exportState=fv%export, clock=clock,                   &
         & phase=ESMF_SINGLEPHASE, rc=rc )
    VERIFY_(rc)

! z_tracer=.false. is necessary to get correct mass flux exports.
! Eventually, fix the FV code that is used when z_tracer=.true. 
    z_tracer = .false.

! The FV dycore default compilation if real(4) and is controlled by the
! SINGLE_FV flag in FVdycoreCubed_GridComp/fvdycore/GNUmakefile.
    call  FV_RESET_CONSTANTS( FV_PI=REAL(pi), &
                            & FV_OMEGA=REAL(omega) ,&
                            & FV_CP=REAL(rgas/kapa) ,&
                            & FV_RADIUS=REAL(radius) ,&
                            & FV_RGAS=REAL(rgas) ,&
                            & FV_RVAP=REAL(rvap) ,&
                            & FV_KAPPA=REAL(kapa) ,&
                            & FV_GRAV=REAL(grav) ,&
                            & FV_HLV=REAL(lhe) ,&
                            & FV_ZVIR=REAL(rvap/rgas-1)  )

    call allocateFvExport3D ( fv%export,'U_DGRID' )
    call allocateFvExport3D ( fv%export,'V_DGRID' )
    call allocateFvExport3D ( fv%export,'TH' )
    call allocateFvExport3D ( fv%export,'PLE' )
    call allocateFvExport3D ( fv%export,'Q' )
    call allocateFvExport3D ( fv%export,'MFX' )
    call allocateFvExport3D ( fv%export,'MFY' )
    call allocateFvExport3D ( fv%export,'MFZ' )

  end subroutine gridCompInit

!-------------------------------------------------------------------------------
  subroutine addQfieldToFVImport(fv)
!-------------------------------------------------------------------------------
    use ESMFL_MOD, Only: ESMFL_StateGetPointerToData
    Type(FV_CORE), intent(inout) :: fv

    type (ESMF_FieldBundle)        :: bundle
    type (ESMF_Field)              :: Qfield

    ! Specific Humidity - need to reserve space in FV import state
    call ESMFL_StateGetPointerToData (fv%export, fv%q, 'Q', rc=rc)
    VERIFY_(rc)
    fv%q = 0.0

    Qfield = ESMF_FieldCreate( fv%grid, fv%q, name='Q', rc=rc)
    VERIFY_(rc)

    ! First obtain a reference field from the state (avoids tedious details)
    call ESMF_StateGet(fv%import, 'TRADV'  , bundle, rc=rc)
    VERIFY_(rc)

    Call ESMF_FieldBundleAdd(bundle, Qfield, rc=rc)
    VERIFY_(rc)

  end subroutine addQfieldToFVImport

!-------------------------------------------------------------------------------
  Subroutine createInternalRestart(fv, istart, cf, clock)
!-------------------------------------------------------------------------------
    USE DOMAIN_DECOMP_ATM, ONLY: GRID, GET, AM_I_ROOT
    Use MAPL_IOMod, only: GETFILE, Free_file, GEOS_VarWrite=>MAPL_VarWrite, Write_parallel
    USE RESOLUTION, only: IM, JM, LM, LS1, PMTOP, Ptop, PSFMPT
    Use DYNAMICS, only: sige, sig
    Use MODEL_COM, only: DT=>DTsrc
    Use ATM_COM, only: U, V, T, P, Q
    Use Constant, only: omega, radius, grav, rgas, kapa, deltx

    Type (FV_Core), Intent(InOut) :: fv
    integer, intent(in) :: istart
    Type (ESMF_Config), Intent(InOut) :: cf
    Type (ESMF_Clock),  Intent(In) :: clock

    Character(Len=ESMF_MAXSTR) :: rst_file
    Integer :: unit

    Integer, Parameter :: N_TRACERS = 0
    Integer :: I_0, I_1, j_0, j_1, j_0h, j_1h, L
    Logical :: exist

    Real(kind=prec) :: ak(size(sige)), bk(size(sige))
    real(kind=prec), allocatable, dimension(:,:,:) :: U_d
    real(kind=prec), allocatable, dimension(:,:,:) :: V_d
    real(kind=prec), allocatable, dimension(:,:,:) :: PE, PKZ, PT
    real(kind=prec), allocatable, dimension(:,:,:) :: W
    real(kind=prec), allocatable, dimension(:,:,:) :: DZ

    Call ESMF_ConfigGetAttribute(cf, value=rst_file, &
         & label='FVCORE_INTERNAL_RESTART_FILE:',    &
         & default=FVCORE_INTERNAL_RESTART,rc=rc)

    if(istart .ge. extend_run) then
    ! Check to see if restart file exists
       inquire(file=FVCORE_INTERNAL_RESTART, EXIST=exist)
       if (exist) then
          if (AM_I_ROOT()) then
             print*,'Using checkpoint file: ',FVCORE_INTERNAL_RESTART
          end if
          return
       end if
       ! Uh oh
       call stop_model('fv part of restart file not found',255)
       return
    end if

    ! If we got to here, then this means then we'll have to create a restart file
    ! from scratch.
    unit = GetFile(rst_file, form="unformatted", rc=rc)
    VERIFY_(rc)

    ! 1) Start date
    Call write_start_date(clock, unit)

    ! 2) Grid size
    Call WRITE_PARALLEL( (/ IM, JM, LM, LM+1-LS1, N_TRACERS /), unit )

    ! 3) Pressure coordinates
    ! Keep in mind that L is reversed between these two models

    Call Compute_ak_bk(ak, bk, sige, Ptop, PSFMPT, unit)
    Call WRITE_PARALLEL(ak, unit)
    Call WRITE_PARALLEL(bk, unit)

    Call GET(grid, i_strt=I_0, i_stop=I_1, &
                   j_strt=j_0, j_stop=j_1, &
                   j_strt_halo=j_0h, j_stop_halo=j_1h)

    ! 4) 3D fields velocities
    Allocate(U_d(I_0:I_1, J_0:J_1, LM))
    Allocate(V_d(I_0:I_1, J_0:J_1, LM))
    U_d = 0.d0
    V_d = 0.d0
    Call GEOS_VarWrite(unit, grid%ESMF_GRID, U_d(I_0:I_1,J_0:J_1,:))
    Call GEOS_VarWrite(unit, grid%ESMF_GRID, V_d(I_0:I_1,J_0:J_1,:))
    Deallocate(V_d)
    Deallocate(U_d)

    ! Compute potential temperature from modelE (1 mb -> 1 pa ref)
    Allocate(PT(I_0:I_1, J_0:J_1, LM))
    Call ConvertPotTemp_GISS2FV(getVirtTemp(T(I_0:I_1,J_0:J_1,:), Q(I_0:I_1,J_0:J_1,:)), PT)
    Call GEOS_VarWrite(unit, grid%ESMF_GRID, PT)
    Deallocate(PT)

    ! Compute PE, PKZ from modelE
    Allocate(PKZ(I_0:I_1, J_0:J_1, LM))
    Allocate(PE(I_0:I_1, J_0:J_1, LM+1))
    Call ComputePressureLevels(unit, grid, getvirtTemp(T, Q), P, SIG, SIGE, Ptop, KAPA, PE, PKZ )

    Call GEOS_VarWrite(unit, grid%ESMF_GRID, PE)
    Call GEOS_VarWrite(unit, grid%ESMF_GRID, PKZ)

    Deallocate(PE)
    Deallocate(PKZ)

!    if (.not. hydrostatic) then
       Allocate(W(I_0:I_1, J_0:J_1, LM))
       Allocate(DZ(I_0:I_1, J_0:J_1, LM))
       W = 0.d0
       DZ = 0.d0
       Call GEOS_VarWrite(unit, grid%ESMF_GRID, DZ(I_0:I_1,J_0:J_1,:))
       Call GEOS_VarWrite(unit, grid%ESMF_GRID, W(I_0:I_1,J_0:J_1,:))
       Deallocate(W)
       Deallocate(DZ)
!    endif

    Call Free_File(unit)

  CONTAINS

    ! Computes virtual pot. temp. from pot. temp. and specific humidity
    !------------------------------------------------------------------
    Function getVirtTemp(T, Q) Result(T_virt)
    !------------------------------------------------------------------
      Use Constant, only: DELTX
      Real(kind=prec), Intent(In) :: T(:,:,:)
      Real(kind=prec), Intent(In) :: Q(:,:,:)
      Real(kind=prec)             :: T_virt(size(T,1), size(T,2), size(T,3))

      T_virt = T * (1 + deltx * Q)

    End Function getVirtTemp

    !------------------------------------------------------------------
    Subroutine ComputePressureLevels(unit,grid,T_virt,P,sig,sige,ptop,kapa,PE,PKZ)
    !------------------------------------------------------------------
      USE DOMAIN_DECOMP_ATM, only: dist_grid, get
      USE RESOLUTION, only: IM, LM
      Integer, intent(in) :: unit
      type (dist_grid) :: grid
      real(kind=prec), dimension(grid%i_strt_halo:,grid%j_strt_halo:,:) :: T_virt
      real(kind=prec), dimension(grid%i_strt_halo:,grid%j_strt_halo:) :: P
      real(kind=prec), dimension(grid%i_strt:grid%i_stop,grid%j_strt:grid%j_stop,LM) :: PKZ
      real(kind=prec), dimension(grid%i_strt:grid%i_stop,grid%j_strt:grid%j_stop,LM+1) :: PE
      real(kind=prec) :: sig(:), sige(:)
      real(kind=prec) :: ptop, kapa

      Integer :: I_0, I_1, i_0h, i_1h, j_0, j_1, j_0h, j_1h
      Integer :: I,J,L,L_fv
      Real(kind=prec), Allocatable :: PK(:,:,:), PELN(:,:,:), PE_trans(:,:,:)

      !    Request local bounds from modelE grid.
      Call GET(grid, i_strt=I_0, i_stop=I_1, j_strt=j_0, j_stop=j_1, &
            & i_strt_halo=i_0h, i_stop_halo=i_1h, j_strt_halo=j_0h, j_stop_halo=j_1h)

      PE = -99999
      PKZ = -99999

      Call ConvertPressure_GISS2FV(EdgePressure_GISS(), PE(i_0:i_1,j_0:j_1,:))

      Allocate(pk(i_0h:i_1h,j_0h:j_1h,LM+1))
      PK(I_0:I_1,J_0:J_1,:) = PE**KAPA
      do l=1,LM
        do j=j_0,j_1
          do i=I_0,I_1
            if (PE(i,j,l+1)-PE(i,j,l) /= 0.0) then
               PKZ(i,j,l) = ( PK(i,j,l+1)-PK(i,j,l) ) / &
                            ( KAPA*log( PE(i,j,l+1)/PE(i,j,l) ) )
            endif
          enddo
        enddo
      enddo
      deallocate(pk)

    End Subroutine ComputePressureLevels

    !------------------------------------------------------------------
    Subroutine Compute_ak_bk(ak, bk, sige, Ptop, PSFMPT, unit)
    !------------------------------------------------------------------
      USE RESOLUTION, only: LM, LS1
      Real(kind=prec) :: sige(:)
      Real(kind=prec) :: Ptop, PSFMPT
      Integer :: unit

      Real(kind=prec), intent(out) :: ak(size(sige)), bk(size(sige))
      Integer :: L, L_fv,k

      Do L = 1, LM+1
         L_fv = LM+2-L
         Select Case(L)
         Case (:LS1-1)
            ak(L_fv) = Ptop*(1-sige(L)) * PRESSURE_UNIT_RATIO ! convert from hPa
            bk(L_fv) = sige(L)
         Case (LS1:)
            ak(L_fv)   = (sige(L)*PSFMPT + Ptop) * PRESSURE_UNIT_RATIO! convert from hPa
            bk(L_fv)   = 0
         End Select
      End Do

    End Subroutine Compute_ak_bk

    !------------------------------------------------------------------
    Subroutine write_start_date(clock, unit)
    !------------------------------------------------------------------
      Type (ESMF_Clock), Intent(In) :: clock
      Integer          , Intent(In) :: unit

      Type (ESMF_Time)  :: currentTime
      Integer :: int_pack(6), YEAR, MONTH, DAY, HOUR, MINUTE, SECOND

      call ESMF_ClockGet(clock, currTime=currentTime, rc=rc)

      call ESMF_TimeGet(currentTime, YY=year, MM= month, &
           DD=day, H=hour, M=minute, &
           S=second, rc=rc)

      INT_PACK(1) = YEAR
      INT_PACK(2) = MONTH
      INT_PACK(3) = DAY
      INT_PACK(4) = HOUR
      INT_PACK(5) = MINUTE
      INT_PACK(6) = SECOND
      Call WRITE_PARALLEL(INT_PACK(1:6), unit=UNIT)
    End Subroutine write_start_date

  End Subroutine createInternalRestart

!-------------------------------------------------------------------------------
  subroutine accumulate_mass_fluxes(fv)
!-------------------------------------------------------------------------------
    use ESMFL_MOD, Only: ESMFL_StateGetPointerToData
    Use Resolution, only: LM
    USE GEOM, ONLY: AXYP
    USE ATM_COM, ONLY: PUA,PVA,SDA
    USE DYNAMICS, ONLY: PU,PV,CONV,SD,PIT
    Use MODEL_COM, only: DT=>DTsrc
    USE DOMAIN_DECOMP_ATM, only: get, grid
    Use Constant, only: grav
    implicit none
    type (FV_core) :: fv
    real*4, Dimension(:,:,:), Pointer :: mfx_X, mfx_Y, mfx_Z
    integer :: I_0, I_1, J_0, J_1
    integer :: i,j,l,k
    integer :: rc
    real(kind=prec) :: DTfac

    DTfac = DT

    Call Get(grid, i_strt=I_0, i_stop=I_1, j_strt=j_0, j_stop=j_1)

    ! Horizontal and Vertical mass fluxes
    !---------------
    !     \item {\tt MFX}:       Mass-Weighted U-Wind on C-Grid (Pa m^2/s)
    !     \item {\tt MFY}:       Mass-Weighted V-wind on C-Grid (Pa m^2/s)
    !     \item {\tt MFZ}:       Vertical mass flux (kg/(m^2*s))

    call ESMFL_StateGetPointerToData ( fv%export,mfx_X,'MFX',rc=rc)
    VERIFY_(rc)
    call ESMFL_StateGetPointerToData ( fv%export,mfx_Y,'MFY',rc=rc)
    VERIFY_(rc)
    call ESMFL_StateGetPointerToData ( fv%export,mfx_Z,'MFZ',rc=rc)
    VERIFY_(rc)

#ifdef NO_MASS_FLUX
    mfx_X = 0
    mfx_Y = 0
    mfx_Z = 0
#endif

    PU(I_0:I_1,J_0:J_1,1:LM) = ReverseLevels(mfx_X)
    PV(I_0:I_1,J_0:J_1,1:LM) = ReverseLevels(mfx_Y)
    mfx_Z = ReverseLevels(mfx_Z)
    SD(I_0:I_1,J_0:J_1,1:LM-1) = mfx_Z(:,:,1:LM-1)

! Add missing factor of grav to FV exports, change units

    pu = pu*grav/PRESSURE_UNIT_RATIO
    pv = pv*grav/PRESSURE_UNIT_RATIO

! Change Units of vertical mass fluxes to mb m^2/s
    do l=1,lm-1
       do j=j_0,j_1
          do i=I_0,I_1
             sd(i,j,l) = grav*grav*AXYP(i,j)*sd(i,j,l)/PRESSURE_UNIT_RATIO
          enddo
       enddo
    enddo

    ! Surface Pressure tendency - vert integral of horizontal convergence
    PIT(I_0:I_1,J_0:J_1) = mfx_Z(:,:,0) + sum(SD(I_0:I_1,J_0:J_1,1:LM-1),3)

    ! Recopy into CONV to support prior usage
    CONV(I_0:I_1,J_0:J_1,1) = PIT(I_0:I_1,J_0:J_1)
    CONV(I_0:I_1,J_0:J_1,2:LM) = SD(I_0:I_1,J_0:J_1,1:LM-1)

    PUA(I_0:I_1,J_0:J_1,:) = PUA(I_0:I_1,J_0:J_1,:) + PU(I_0:I_1,J_0:J_1,:)*DTfac
    PVA(I_0:I_1,J_0:J_1,:) = PVA(I_0:I_1,J_0:J_1,:) + PV(I_0:I_1,J_0:J_1,:)*DTfac
    SDA(I_0:I_1,J_0:J_1,1:LM-1) = SDA(I_0:I_1,J_0:J_1,1:LM-1) + SD(I_0:I_1,J_0:J_1,1:LM-1)*DTfac

  end subroutine accumulate_mass_fluxes

end module FV_CS_Mod
