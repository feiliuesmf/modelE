#define VERIFY_(rc) If (rc /= ESMF_SUCCESS) Call abort_core(__LINE__,rc)
#define RETURN_(status) If (Present(rc)) rc=status; return


!#define NO_FORCING


!---------------------------------------------------------------------------------------------
! The following module provides an interface for modelE to use the ESMF-wrapped version of
! the FV dynamical core from GEOS-5.  Many elements closely resemble an ESMF coupler component,
! which is to be expected given the intended role for this layer.
!---------------------------------------------------------------------------------------------

module FV_INTERFACE_MOD
  use esmf_mod
!!$  use GEOS_Mod, Only: ESMFL_StateGetPointerToData
  use ESMFL_MOD, Only: ESMFL_StateGetPointerToData
  implicit none
  private

  ! except for

  public :: FV_CORE                ! Derived type to encapsulate FV + modelE data
  public :: Init_app_clock         ! modelE does not use ESMF clocks, but ESMF requires one
  public :: Initialize             ! Uses modelE data to initialize the FV gridded component
  public :: Finalize               ! Uses modelE data to finalize the FV gridded component
  public :: Compute_Tendencies     ! Temporarily a separate method, but will move within Run() soon.
  public :: Run                    ! Execute the FV method to integrate the core forward in time
  public :: Checkpoint             ! Unimplemented
  public :: edgepressure_giss
  public :: DryTemp_GISS
  ! This data structure is o convenient entity for storing persistent data between
  ! calls to this module.  In addition to
  Type FV_CORE
     PRIVATE

     type (esmf_gridcomp) :: gc   ! This is the handle for the fv dynamical core

     type (esmf_grid)     :: grid ! Although modelE is not an ESMF component, it does have an ESMF_Grid
     type (esmf_vm)       :: vm   ! Should be eliminated ... Only used to obtain NPES for config file.

     ! Import and Export states for FV dycore
     type(esmf_state) :: import   ! Allocated within FV component
     type(esmf_state) :: export   ! Allocated within FV component

     ! The following pointers can be re-extracted from the import state at each iteration,
     ! but it is convenient to have a simpler means of access.
     real*4, pointer, dimension(:,:,:) :: dudt, dvdt, dtdt, dpedt  ! Tendencies
     real*4, pointer, dimension(:,:,:) :: Q  ! Humidity
     real*4, pointer, dimension(:,:,:) :: Qtr  ! other tracers
     real*4, pointer, dimension(:,:)   :: phis


     ! modelE does not work directly with tendencies.  Instead, tendencies are derived
     ! by differencing before and after physics.  Therefore, the final dynamical state
     ! must be preserved for the following
     ! of modelE fields
     real*8, pointer, dimension(:,:,:) :: U_old, V_old, dPT_old, PE_old, dT_old

  END Type FV_CORE

  ! private data of convenience
  Integer :: rc ! return code from ESMF

  Interface Initialize
     module procedure initialize_fv
  End Interface

  Interface Run
     module procedure run_fv
  End Interface

  Interface reverse
     Module Procedure reverse_3d_r8
     Module Procedure reverse_3d_r4
  End Interface

  Interface ConvertPotTemp_GISS2FV
     Module Procedure CnvPotTemp_GISS2FV_r8
     Module Procedure CnvPotTemp_GISS2FV_r4
  End Interface

  Interface ConvertPotTemp_FV2GISS
     Module Procedure CnvPotTemp_FV2GISS_r4
  End Interface

  Interface ConvertPressure_GISS2FV
     Module Procedure ConvertPressure_GISS2FV_r4
     Module Procedure ConvertPressure_GISS2FV_r8
  End Interface

  ! The following parameters address the fact that FV and modelE use
  ! different units and different reference pressures for potential temperature.
  ! Superficially there is redundancy between the two sets, but in some sense
  ! this is merely coincidental.

  Real*8, parameter :: REF_PRESSURE_GISS = 100 ! 1 mb = 100 pa
  Real*8, parameter :: REF_PRESSURE_FV   =   1 ! pa
  Real*8, parameter :: REF_RATIO = REF_PRESSURE_FV / REF_PRESSURE_GISS

  Real*8, parameter :: PRESSURE_UNIT_GISS  =  100 ! 1 mb
  Real*8, parameter :: PRESSURE_UNIT_FV    =    1 ! 1 pa
  Real*8, parameter :: PRESSURE_UNIT_RATIO = PRESSURE_UNIT_GISS/PRESSURE_UNIT_FV

  character(len=*), parameter :: FVCORE_INTERNAL_RESTART = 'fv_internal_restart.dat'
  character(len=*), parameter :: FVCORE_IMPORT_RESTART   = 'fv_import_restart.dat'
  character(len=*), parameter :: FVCORE_LAYOUT           = 'fv_layout.rc'

  character(len=*), parameter :: TENDENCIES_FILE = 'tendencies_checkpoint'

  integer, parameter :: FROM_OBSERVED_DATA = 2
  integer, parameter :: FROM_LATEST_SAVE_FILE = 10

contains

  subroutine Initialize_fv(fv, istart, vm, grid, clock, config_file)
    use resolution, only: IM, LM
    use fvdycore_gridcompmod, only: SetServices
    use GEOS_BaseMod, only: GEOS_FieldCreate
    use DOMAIN_DECOMP, only : modelE_grid => grid, get

    type (fv_core),    intent(inout) :: fv
    integer,           intent(in) :: istart
    type (esmf_vm),    intent(in) :: vm
    type (esmf_grid),  intent(in) :: grid
    type (esmf_clock), intent(in) :: clock
    character(len=*),  intent(in) :: config_file ! filename for resource file
    type (esmf_config) :: cf
    Integer :: L
    Integer :: j_0, j_1
    type (ESMF_Bundle)               :: bundle
    type (ESMF_Field)                :: Qfield
    type (ESMF_Array)                :: Qarray
    type (ESMF_FieldDataMap)         :: datamap

    fv % vm  =vm
    fv % grid=grid

    ! Load configuration information from resource file
    !  ------------------------------------------------
    cf = load_configuration(config_file)

    !  Create the dynamics component, using same layout as application
    !  ------------------------------------------------------------------
    fv % gc = ESMF_GridCompCreate ( vm=vm, name='FV dynamics', &
         & grid=grid, gridcomptype=ESMF_ATM,              &
         & config=cf, rc=rc)
    VERIFY_(rc)

    ! Create couplings
    fv % import = ESMF_StateCreate ( 'fv dycore imports', ESMF_STATE_IMPORT, rc=rc )
    VERIFY_(rc)
    fv % export = ESMF_StateCreate ( 'fv dycore exports', ESMF_STATE_EXPORT, rc=rc )
    VERIFY_(rc)

    !  Register services for components
    !  --------------------------------
!   print*,'calling set services'
    call ESMF_GridCompSetServices ( fv % gc, SetServices, rc )
    VERIFY_(rc)


    ! The FV components requires its own restart file for managing
    ! its internal state.  We check to see if the file already exists, and if not
    ! create one based upon modelE's internal state.
    call create_restart_file(fv, istart, cf, clock)

!   print*,'calling fv initialize'
    call ESMF_GridCompInitialize ( fv % gc, importState=fv % import, exportState=fv % export, clock=clock, &
         & phase=ESMF_SINGLEPHASE, rc=rc )
    VERIFY_(rc)

    ! Initialize component (and import/export states)
    !  ----------------------------------------------
    Call allocate_tendency_storage(fv, istart)

    call allocateFvExport3D ( fv % export,'U' )
    call allocateFvExport3D ( fv % export,'V' )
    call allocateFvExport3D ( fv % export,'TH' )
    call allocateFvExport3D ( fv % export,'PLE' )
    call allocateFvExport3D ( fv % export,'Q' )
    call allocateFvExport3D ( fv % export,'MFX_A' )
    call allocateFvExport3D ( fv % export,'MFY_A' )
    call allocateFvExport3D ( fv % export,'MFZ' )


!   print*,'called fv initialize'

    ! Specific Humidity - need to reserve space in FV import state
    call ESMFL_StateGetPointerToData ( fv % export,fv % q,'Q',rc=rc)
    VERIFY_(rc)
    Call get(modelE_grid, J_STRT=J_0, J_STOP=J_1)
!!$    Allocate(fv % Qtr(IM,J_0:J_1,LM))


    Qarray = ESMF_ArrayCreate(fv % q, ESMF_DATA_REF, RC=rc)
    VERIFY_(rc)
        call ESMF_FieldDataMapSetDefault(datamap, dataRank=3, rc=rc)
    Qfield = ESMF_FieldCreate( grid, Qarray, horzRelloc=ESMF_CELL_CENTER, &
         datamap=datamap, name="Q", rc=rc)
    VERIFY_(rc)
    ! First obtain a reference field from the state (avoids tedious details)
    call ESMF_StateGetBundle (fv % import, 'QTR'  , bundle, rc=rc)
    VERIFY_(rc)
    Call ESMF_BundleAddField(bundle, Qfield, rc=rc)
    VERIFY_(rc)

    contains

    subroutine allocateFvExport3D ( state, name )
    type(ESMF_State),  intent(INOUT) :: state
    character(len=*),  intent(IN   ) :: name

    real, pointer                    :: ptr(:,:,:)
    logical :: alloc
    integer :: status

    alloc = .true.
    call ESMFL_StateGetPointerToData ( state, ptr , name , alloc , status )
    VERIFY_(status)

    end subroutine allocateFvExport3D

  end subroutine Initialize_fv

  function EdgePressure_GISS() Result(PE)
    USE RESOLUTION, only: IM, LM, LS1
    Use MODEL_COM, only : SIG, SIGE, Ptop, PSFMPT, P
    use DOMAIN_DECOMP, only: grid, get
    USE CONSTANT, only: KAPA

    REAL*8 :: PE(IM,grid % J_STRT:grid % J_STOP,LM+1)

    INTEGER :: L, j_0, j_1

    call get(grid, J_STRT=J_0, J_STOP=J_1)

    Do L = 1, LM+1

       If (L < LS1) THEN
          PE(:,:,L) = SIGE(L)*P(:,J_0:J_1) + Ptop
       Else
          PE(:,:,L)  = SIGE(L)*PSFMPT + Ptop
       End IF

    End Do

  end function EdgePressure_GISS

  ! Compute Delta-pressure for GISS model
  function DeltPressure_GISS() Result(dP)
    USE RESOLUTION, only: IM, LM, LS1
    Use MODEL_COM, only: T
    USE DOMAIN_DECOMP, only: grid, GET

    REAL*8 :: dP(IM,grid % J_STRT:grid % J_STOP,LM)
    REAL*8 :: PE(IM,grid % J_STRT:grid % J_STOP,LM+1)

    INTEGER :: J_0,J_1,k
    Call Get(grid, J_STRT=J_0, J_STOP=J_1)

    PE = EdgePressure_GISS()
    do k = 1, LM
       dP(:,:,k) = (PE(:,:,k)-PE(:,:,k+1))
    end do

  end function DeltPressure_GISS


  ! Convert Potential Temperature into (dry) Temperature
  function DeltPressure_DryTemp_GISS() Result(dPT)
    USE RESOLUTION, only: IM, LM, LS1
    Use MODEL_COM, only: T
    USE DOMAIN_DECOMP, only: grid, GET

    REAL*8 :: dPT(IM,grid % J_STRT:grid % J_STOP,LM)
    REAL*8 :: PE(IM,grid % J_STRT:grid % J_STOP,LM+1)
    REAL*8 :: T_dry(IM,grid % J_STRT:grid % J_STOP,LM)

    INTEGER :: J_0,J_1,k
    Call Get(grid, J_STRT=J_0, J_STOP=J_1)

    T_dry = DryTemp_GISS()
    PE = EdgePressure_GISS()
    do k = 1, LM
       dPT(:,:,k) = (PE(:,:,k)-PE(:,:,k+1)) * T_dry(:,:,k)
    end do

  end function DeltPressure_DryTemp_GISS


  Subroutine allocate_tendency_storage(fv, istart)
    Use DOMAIN_DECOMP, only: GRID, GET, AM_I_ROOT
    USE RESOLUTION, only: IM, LM, LS1
    USE MODEL_COM, only: U, V, T
    use FILEMANAGER
    Implicit None
    Type (FV_Core) :: fv
    integer, intent(in) :: istart

    Integer :: J_0H, J_1H, J_0, J_1
    integer :: iunit

    Call Get(grid, J_strt_halo=J_0H, J_stop_halo=J_1H, &
         & J_STRT=J_0, J_STOP=J_1)

    ! 1) Link/copy modelE data to import state
    call ESMFL_StateGetPointerToData ( fv % import,fv % dudt,'DUDT',rc=rc)
    VERIFY_(rc)
    call ESMFL_StateGetPointerToData ( fv % import,fv % dvdt,'DVDT',rc=rc)
    VERIFY_(rc)
    call ESMFL_StateGetPointerToData ( fv % import,fv % dtdt,'DTDT',rc=rc)
    VERIFY_(rc)
    call ESMFL_StateGetPointerToData ( fv % import,fv % dpedt,'DPEDT',rc=rc)
    VERIFY_(rc)
    call ESMFL_StateGetPointerToData ( fv % import,fv % phis, 'PHIS', rc=rc)
    VERIFY_(rc)

    ! 2) Allocate space for storing old values from which to compute tendencies
    Allocate(fv % U_old(IM,J_0:J_1,LM), &
         &   fv % V_old(IM,J_0:J_1,LM), &
         &   fv % dPT_old(IM,J_0:J_1,LM), &
         &   fv % dT_old(IM,J_0:J_1,LM), &
         &   fv % PE_old(IM,J_0:J_1,LM+1))

    select case (istart)
    case (FROM_OBSERVED_DATA)
       ! Do a cold start.  Set Old = Current.
       fv % dudt=0
       fv % dvdt=0
       fv % dtdt=0
       fv % dpedt=0

       fv % U_old = U(:,J_0:J_1,:)
       fv % V_old = V(:,J_0:J_1,:)
       fv % dPT_old = DeltPressure_DryTemp_GISS()
       fv % dT_old = DryTemp_GISS()
       fv % PE_old   = EdgePressure_GISS()
    case (FROM_LATEST_SAVE_FILE)
       ! input couplings are in a file somewhere and already read in
       if (AM_I_ROOT()) call openunit(TENDENCIES_FILE, iunit, qbin=.true.,qold=.true.)
       call readArr(iunit, fv % U_old)
       call readArr(iunit, fv % V_old)
       call readArr(iunit, fv % dPT_old)
       call readArr(iunit, fv % PE_old)
       call readArr(iunit, fv % dT_old)
       if (AM_I_ROOT()) call closeunit(iunit)
       call compute_tendencies(fv)
       if (AM_I_ROOT()) call closeUnit(iunit)
    case default
       call stop_model('ISTART option not supported',istart)
    end select


  contains

    subroutine readArr(iunit, arr)
      use resolution, only: IM, JM
      use DOMAIN_DECOMP, only: grid, unpack_data, get
      integer, intent(in) :: iunit
      real*8, intent(out) :: arr(:,:,:)
      real*8, allocatable :: padArr(:,:,:)
      real*8, allocatable :: globalArr(:,:,:)

      allocate(globalArr(IM,JM,size(arr,3)))
      allocate(padArr(1:IM,j_0h:j_1h, size(arr,3)))

      if (AM_I_ROOT()) read(iunit) globalArr
      call unpack_data(grid, globalArr, padArr)
      arr(:,:,:) = padArr(:,j_0:j_1,:)

      deallocate(padArr)
      deallocate(globalArr)

    end subroutine readArr

  End Subroutine allocate_tendency_storage

  subroutine finalize(fv, clock)
    Type (FV_Core) :: fv
    type (esmf_clock), intent(in) :: clock
    integer :: rc

    call saveTendencies(fv)
    call ESMF_GridCompFinalize ( fv % gc, fv % import, fv % export, clock, rc )
    Deallocate(fv % U_old, fv % V_old, fv % dPT_old, fv % dT_old, fv % PE_old)

  end subroutine finalize


  ! Compute tendencies
  ! ------------------
  subroutine compute_tendencies(fv)
    USE RESOLUTION, only: IM, LM
    USE MODEL_COM, Only: U, V, T
    USE DOMAIN_DECOMP, only: grid, get
    USE DYNAMICS, only: DUT, DVT
    USE CONSTANT, only: KAPA
    Implicit None
    Type (FV_CORE) :: fv

    Integer :: J_0, J_1
    Integer :: J_0h, J_1h
    Call Get(grid, j_strt=j_0, j_stop=j_1)

    ! U, V
    DUT(:,J_0:J_1,:) = Tendency(U(:,J_0:J_1,:), fv % U_old(:,J_0:J_1,:))
    DVT(:,J_0:J_1,:) = Tendency(V(:,J_0:J_1,:), fv % V_old(:,J_0:J_1,:))
    call ConvertUV_GISS2FV(DUT(:,J_0:J_1,:), DVT(:,J_0:J_1,:), fv % dudt, fv % dvdt)

    ! delta pressure weighted Temperature
    fv  %  dtdt = reverse(DeltPressure_GISS() * Tendency(DryTemp_GISS(), fv % dT_old)) * &
         & (PRESSURE_UNIT_RATIO)

    ! Edge Pressure
       Call ConvertPressure_GISS2FV( Tendency(EdgePressure_GISS(), fv % PE_old), fv % dpedt)

#ifdef NO_FORCING
    fv % dudt = 0
    fv % dvdt = 0
    fv % dTdt = 0
    fv % dPEdt = 0
#endif

!!$$    call write_profile(real(fv % dudt,8),'Tendency U')
!!$$    call write_profile(real(fv % dVdt,8),'Tendency V')
!!$$    call write_profile(real(fv % dTdt,8),'Tendency T')
!!$$    call write_profile(real(fv % dPEdt,8),'Tendency PE')

  end subroutine compute_tendencies

  !-----------
  ! d/dt A
  ! Note that FV needs tendency with the order of vertical levels reversed.
  !-----------
  Function Tendency(A, A_old) Result(tend)
    USE MODEL_COM, Only: DTsrc
    REAL*8, INTENT(IN) :: A(:,:,:)
    REAL*8, INTENT(IN) :: A_old(:,:,:)
    REAL*8 :: tend(size(a,1),size(a,2),size(a,3))

    tend = (A - A_old)/DTsrc

  End Function Tendency

  subroutine run_fv(fv, clock)
    USE DOMAIN_DECOMP, only: grid, halo_update, get, grid, NORTH
    USE MODEL_COM, Only : U, V, T, P, IM, JM, LM, ZATMO
    USE MODEL_COM, only : NIdyn, DT, DTSRC
    USE SOMTQ_COM, only: TMOM, MZ
    USE ATMDYN, only: CALC_AMP, CALC_PIJL, AFLUX, COMPUTE_MASS_FLUX_DIAGS
    USE DYNAMICS, only: MA, PHI, GZ
    USE DYNAMICS, ONLY: PU, PV, CONV
    USE DYNAMICS, ONLY: SD, AFLUX

    Type (FV_CORE)    :: fv
    Type (ESMF_Clock) :: clock
    type (ESMF_TimeInterval) :: timeInterval

    REAL*8, DIMENSION(IM,grid % J_STRT_HALO:grid % J_STOP_HALO,LM) :: PIJL
    integer :: istep, NS, NIdyn_fv
    integer :: rc

!@sum  CALC_AMP Calc. AMP: kg air*grav/100, incl. const. pressure strat
    call calc_amp(P, MA)
    CALL CALC_PIJL(LM,P,PIJL)

    Call Copy_modelE_to_FV_import(fv)

    call clear_accumulated_mass_fluxes()
    ! Run dycore
    NIdyn_fv = DTsrc / (DT)
    do istep = 1, NIdyn_fv

       call ESMF_GridCompRun ( fv % gc, fv % import, fv % export, clock, 91, rc=rc )
       call ESMF_GridCompRun ( fv % gc, fv % import, fv % export, clock, rc )

       call ESMF_TimeIntervalSet(timeInterval, s = nint(DT), rc=rc)
       call ESMF_ClockAdvance(clock, timeInterval, rc=rc)

       call accumulate_mass_fluxes(fv)
       TMOM = 0 ! for now
       call Copy_FV_export_to_modelE(fv) ! inside loop to accumulate PUA,PVA,SDA
       phi = compute_phi(P, T, TMOM(MZ,:,:,:), ZATMO)
       call compute_mass_flux_diags(phi, pu, pv, dt)
    end do

    gz  = phi
  end subroutine run_fv

  subroutine checkpoint(fv, clock)
    use GEOS_mod, only: GEOS_RecordPhase
    use DOMAIN_DECOMP, only: AM_I_ROOT
    Type (FV_Core),    intent(inout) :: fv
    type (esmf_clock), intent(in) :: clock
    integer :: rc

    character(len=*), parameter :: SUFFIX_TEMPLATE = '.YYYYMMDD_HHMMz.bin'
    character(len=len(SUFFIX_TEMPLATE)) :: suffix
    Type (ESMF_Time)  :: currentTime
    integer :: year, month, day, hour, minute, second

    call ESMF_GridCompFinalize( fv % gc, fv % import, fv % export, clock, GEOS_RecordPhase, rc=rc)
    ! Now move the file into a more useful name
    if (AM_I_ROOT()) then
    ! 1) FV names restarts as: fvcore_internal_checkoint.YYYYMMDD_HHMMz.bin
       call ESMF_ClockGet(clock, currTime=currentTime, rc=rc)
       call ESMF_TimeGet(currentTime, YY=year, MM= month, &
            DD=day, H=hour, M=minute, &
            S=second, rc=rc)
       write(suffix,'(".",i4.4,i2.2,i2.2,"_",i2.2,i2.2,"z.bin")'), &
            & year, month, day, hour, minute
       call system('mv ' // FVCORE_INTERNAL_RESTART // suffix // ' ' // FVCORE_INTERNAL_RESTART)
    end if

    call saveTendencies(fv)
       

  end subroutine checkpoint

  subroutine saveTendencies(fv)
    use DOMAIN_DECOMP, only: AM_I_ROOT
    use FILEMANAGER
    type (fv_core),    intent(inout) :: fv
    integer :: iunit

    if (AM_I_ROOT()) call openunit(TENDENCIES_FILE, iunit, qbin=.true.)

    call saveArr(iunit, fv % U_old)
    call saveArr(iunit, fv % V_old)
    call saveArr(iunit, fv % dPT_old)
    call saveArr(iunit, fv % PE_old)
    call saveArr(iunit, fv % dT_old)

    if (am_i_root()) call closeunit(iunit)

  contains

    subroutine saveArr(iunit, arr)
      use resolution, only: IM, JM
      use DOMAIN_DECOMP, only: grid, pack_data, get
      integer, intent(in) :: iunit
      real*8, intent(in) :: arr(:,:,:)
      real*8, allocatable :: padArr(:,:,:)
      real*8, allocatable :: globalArr(:,:,:)
      integer :: j_0, j_1
      integer :: j_0h, j_1h

      Call Get(grid, j_strt=j_0, j_stop=j_1, &
           & j_strt_halo = j_0h, j_stop_halo = j_1h)
      allocate(globalArr(IM,JM,size(arr,3)))
      allocate(padArr(1:IM,j_0h:j_1h, size(arr,3)))
      
      padArr(:,j_0:j_1,:) = arr(:,:,:)
      call pack_data(grid, padArr, globalArr)
      if (AM_I_ROOT()) write(iunit) globalArr

      deallocate(padArr)
      deallocate(globalArr)

    end subroutine saveArr

  end subroutine saveTendencies

  !----------------------------
  !  Internal routines
  !----------------------------

  !----------------------------
  function init_app_clock(start_time, end_time, interval) Result(clock)
    integer :: start_time(6)
    integer :: end_time(6)
    integer :: interval
    type (esmf_clock)              :: clock

    type (esmf_time) :: startTime
    type (esmf_time) :: stopTime
    type (esmf_timeinterval) :: timeStep
    type (esmf_calendar) :: gregorianCalendar

    ! initialize calendar to be Gregorian type
    gregorianCalendar = esmf_calendarcreate("GregorianCalendar", ESMF_CAL_GREGORIAN, rc)
    VERIFY_(rc)

    ! initialize start time
    write(*,*)'Time Set Start: ',START_TIME
    call esmf_timeset(startTime, YY=START_TIME(1), MM= START_TIME(2), &
         & DD=START_TIME(3), H=START_TIME(4),                         &
         & M=START_TIME(5),  S=START_TIME(6),                         &
         & calendar=gregorianCalendar, rc=rc)
    VERIFY_(rc)

    ! initialize stop time
    write(*,*)'Time Set End: ',END_TIME
    call esmf_timeset(stopTime, YY=END_TIME(1), MM= END_TIME(2), &
         & DD=END_TIME(3), H=END_TIME(4),                        &
         & M=END_TIME(5),  S=END_TIME(6),                        &
         & calendar=gregorianCalendar, rc=rc)
    VERIFY_(rc)

    ! initialize time interval
    call esmf_timeintervalset(timeStep, &
         & S=INT(interval), rc=rc)
    VERIFY_(rc)

    ! initialize the clock with the above values
    clock = esmf_clockcreate("ApplClock", timeStep, startTime, stopTime, rc=rc)
    VERIFY_(rc)

  end function init_app_clock
  !----------------------------

  !----------------------------
  function load_configuration(config_file) result( config )
    use FILEMANAGER
    Use MODEL_COM,  only: DT
    character(len=*), parameter :: Iam="FV_INTERFACE::loadconfiguration"
    character(len=*), intent(in) :: config_file
    type (esmf_config)           :: config

    integer :: iunit

    config = esmf_configcreate(rc=rc)
    VERIFY_(rc)

    call openunit(config_file, iunit, qbin=.false., qold=.false.)
    write(iunit,*)'FVCORE_INTERNAL_CHECKPOINT_FILE:  ', FVCORE_INTERNAL_RESTART
    write(iunit,*)'FVCORE_INTERNAL_RESTART_FILE:     ', FVCORE_INTERNAL_RESTART
!!$    write(iunit,*)'FVCORE_IMPORT_CHECKPOINT_FILE:    ', TENDENCIES_FILE
!!$    write(iunit,*)'FVCORE_IMPORT_RESTART_FILE:       ', TENDENCIES_FILE
    write(iunit,*)'FVCORE_LAYOUT_FILE:               ', FVCORE_LAYOUT
    write(iunit,*)'RUN_DT:                           ', DT
    call closeUnit(iunit)

    call esmf_configloadfile(config, config_file, rc=rc)
    VERIFY_(rc)

  end function load_configuration
  !----------------------------

  Subroutine Create_Restart_File(fv, istart, cf, clock)
    USE DOMAIN_DECOMP, ONLY: GRID, GET, AM_I_ROOT
    Use GEOS_IOMod, only: GETFILE, Free_file, GEOS_VarWrite, Write_parallel
    USE RESOLUTION, only: IM, JM, LM, LS1
    Use MODEL_COM, only: sige, sig, Ptop, DT, PMTOP
    Use MODEL_COM, only: U, V, T, P, PSFMPT, Q
    USE GEOM, only: AREAG, DXYP
    Use Constant, only: omega, radius, grav, rgas, kapa, deltx

    Type (FV_Core), Intent(InOut) :: fv
    integer, intent(in) :: istart
    Type (ESMF_Config), Intent(InOut) :: cf
    Type (ESMF_Clock),  Intent(In) :: clock

    Character(Len=ESMF_MAXSTR) :: rst_file
    Integer :: unit

    Integer, Parameter :: N_TRACERS = 0
    Integer :: j_0, j_1, j_0h, j_1h, L
    Logical :: exist

    Real*8 :: ak(size(sige)), bk(size(sige))

    real*8, allocatable, dimension(:,:,:) :: U_d
    real*8, allocatable, dimension(:,:,:) :: V_d
    real*8, allocatable, dimension(:,:,:) :: PE, PKZ, PT


    ! 1) Create layout resource file - independent of actual restart file.
    If (AM_I_ROOT()) Then
       Call Write_Layout(FVCORE_LAYOUT, fv)
    End If

    Call ESMF_ConfigGetAttribute(cf, value=rst_file, label='FVCORE_INTERNAL_RESTART_FILE:', &
         & default=FVCORE_INTERNAL_RESTART,rc=rc)

    ! Check to see if restart file already exists
    print*,'Looking for an FV restart ...',istart
    select case (istart)
    case (FROM_LATEST_SAVE_FILE)
       inquire(file=FVCORE_INTERNAL_RESTART, EXIST=exist)
       if (exist) then
          if (AM_I_ROOT()) then
             print*,'Using checkpoint file: ',FVCORE_INTERNAL_RESTART
          end if
          return
       end if
    case (FROM_OBSERVED_DATA)
       ! Forcing creation of FV restart from GISS data for now
!!$       ! possibly exists from somewhere else?
!!$       inquire(file=rst_file,EXIST=exist)
!!$       if (exist) then
!!$          if (AM_I_ROOT()) Then
!!$             print*,'Apparently a restart file for FV already exists: ',trim(rst_file)
!!$             print*,'Returning to main program'
!!$             print*,' '
!!$          end if
!!$          return
!!$       end if
    case default
       ! Uh oh
       call stop_model('ISTART value not supported for FV restart.',ISTART)
       return
    end select

    ! If we got to here, then this means then we'll have to create a restart file
    ! from scratch.
    unit = GetFile(rst_file, form="unformatted", rc=rc)
    VERIFY_(rc)
    
    ! 1) Start date
    Call write_start_date(clock, unit)
    
    ! 2) Grid size
    Call WRITE_PARALLEL( (/ IM, JM, LM, LS1-1, N_TRACERS /), unit )
    
    ! 3) Pressure coordinates
    ! Keep in mind that L is reversed between these two models
    
    Call Compute_ak_bk(ak, bk, sige, Ptop, PSFMPT, unit)
    Call WRITE_PARALLEL( ak, unit)
    Call WRITE_PARALLEL( bk, unit)
    
    Call GET(grid, j_strt=j_0, j_stop=j_1, j_strt_halo=j_0h, j_stop_halo=j_1h)
    
    ! 4) 3D fields velocities
    Allocate(U_d(IM, J_0:J_1, LM))
    Allocate(V_d(IM, J_0:J_1, LM))
    
    write(*,*)'Calling ComputeRestartVelocities()'
    Call ComputeRestartVelocities(unit, grid, U, V, U_d, V_d)
    
!!$      call set_zonal_flow(U_d, V_d, j_0, j_1)
    
    Call GEOS_VarWrite(unit, grid % ESMF_GRID, U_d(:,J_0:J_1,:))
    Call GEOS_VarWrite(unit, grid % ESMF_GRID, V_d(:,J_0:J_1,:))
    
    Deallocate(V_d)
    Deallocate(U_d)
    
    ! Compute potential temperature from modelE (1 mb -> 1 pa ref)
    Allocate(PT(IM, J_0:J_1, LM))
    Call ConvertPotTemp_GISS2FV(VirtualTemp(T(:,J_0:J_1,:), Q(:,J_0:J_1,:)), PT)
    Call GEOS_VarWrite(unit, grid % ESMF_GRID, PT)
    Deallocate(PT)
    
    ! Compute PE, PKZ from modelE
    Allocate(PKZ(IM, J_0:J_1, LM))
    Allocate(PE(IM, J_0:J_1, LM+1))
    Call ComputePressureLevels(unit, grid, VirtualTemp(T, Q), P, SIG, SIGE, Ptop, KAPA, PE, PKZ )
    
    Call GEOS_VarWrite(unit, grid % ESMF_GRID, PE)
    Call GEOS_VarWrite(unit, grid % ESMF_GRID, PKZ)
    
    Deallocate(PE)
    Deallocate(PKZ)
    
    Call Free_File(unit)
          
  CONTAINS

    ! Computes virtual pot. temp. from pot. temp. and specific humidity
    !------------------------------------------------------------------
    Function VirtualTemp(T, Q) Result(T_virt)
      Use Constant, only: DELTX
      Real*8, Intent(In) :: T(:,:,:)
      Real*8, Intent(In) :: Q(:,:,:)
      Real*8             :: T_virt(size(T,1), size(T,2), size(T,3))

      T_virt = T * (1 + deltx * Q)

    End Function VirtualTemp

    Subroutine ComputePressureLevels(unit, grid, T_virt, P, sig, sige, ptop, kapa, PE, PKZ)
      USE DOMAIN_DECOMP, only: dist_grid, get
      USE RESOLUTION, only: IM, LM, LS1
      Integer, intent(in) :: unit
      type (dist_grid) :: grid
      real*8, dimension(:,grid % j_strt_halo:,:) :: T_virt
      real*8, dimension(:,grid % j_strt_halo:) :: P
      real*8, dimension(IM,grid % j_strt:grid % j_stop,LM) :: PKZ
      real*8, dimension(IM,grid % j_strt:grid % j_stop,LM+1) :: PE
      real*8 :: sig(:), sige(:)
      real*8 :: ptop, kapa

      Integer :: L, L_fv
      Integer :: j_0, j_1, j_0h, j_1h
      Integer :: J
      Real*8, Allocatable :: PK(:,:,:), PELN(:,:,:), PE_trans(:,:,:)

      !    Request local bounds from modelE grid.
      Call GET(grid, j_strt=j_0, j_stop=j_1, j_strt_halo=j_0h, j_stop_halo=j_1h)

      PE = -99999
      Allocate(pe_trans(IM,LM+1,j_0h:j_1h),pk(IM,j_0h:j_1h,LM+1), peln(IM,LM+1,j_0h:j_1h))

      Call ConvertPressure_GISS2FV(EdgePressure_GISS(), PE(:,j_0:j_1,:))

      ! Compute PKZ - Consistent mean for p^kappa
      ! use pkez() routine within the FV dycore component.

      Do j = j_0, j_1
         PE_trans(:,:,j) = PE(:,j,:)
      END DO
      pk=0
      peln=0
      pkz=0

      CALL pkez(1, IM, LM, J_0, J_1, 1, LM, 1, IM, PE_trans(:,:,j_0:j_1), &
           &  PK(:,j_0:j_1,:), KAPA, LS1-1, PELN(:,:,j_0:j_1), pkz(:,j_0:j_1,:), .true.)

      deallocate(pe_trans, pk, peln)

    End Subroutine ComputePressureLevels

    Subroutine ComputeRestartVelocities(unit, grid, U_b, V_b, U_d, V_d)
      use domain_decomp, only: DIST_GRID, NORTH, HALO_UPDATE

      Integer, intent(in) :: unit
      Type (Dist_Grid) :: grid
      real*8, dimension(:,grid % j_strt_halo:,:) :: U_b, V_b

      real*8, intent(out) :: U_d(:,:,:), V_d(:,:,:)

      Call HALO_UPDATE(grid, U_b, FROM=NORTH)
      Call Regrid_B_to_D(Reverse(U_b), Reverse(V_b), U_d, V_d)

    End Subroutine ComputeRestartVelocities

    Subroutine Compute_ak_bk(ak, bk, sige, Ptop, PSFMPT, unit)
      USE RESOLUTION, only: LM, LS1
      Real*8 :: sige(:)
      Real*8 :: Ptop, PSFMPT
      Integer :: unit

      Real*8, intent(out) :: ak(size(sige)), bk(size(sige))
      Integer :: L, L_fv

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

    Subroutine write_layout(fname, fv)
      USE RESOLUTION, only: IM, JM, LM, LS1
      Use MODEL_COM,  only: DT
      character(len=*), intent(in) :: fname
      type (fv_core), intent(in)   :: fv

      Type (esmf_axisindex), pointer :: AI(:,:)
      Integer :: unit
      Integer :: npes

      call esmf_vmget(fv % vm, petcount = npes, rc=rc)
      Allocate(AI(npes,3))
      call esmf_gridgetallaxisindex(fv % grid, globalai=ai, horzrelloc=ESMF_CELL_CENTER,  &
           &       vertRelLoc=ESMF_CELL_CELL, rc=rc)

      unit = GetFile(fname, form="formatted", rc=rc)
      write(unit,*)' # empty line #'
      Write(unit,*)'xy_yz_decomp:',1,npes,npes,1
      Write(unit,*)'    im: ',IM
      Write(unit,*)'    jm: ',JM
      Write(unit,*)'    km: ',LM
      Write(unit,*)'    dt: ',DT
      Write(unit,*)'nsplit: ',0
      Write(unit,*)' ntotq: ',1
      Write(unit,*)'    nq: ',1
      Write(unit,*)'  imxy: ',IM
      Write(unit,'(a,100(1x,I4))')'  jmxy: ',AI(:,2) % MAX-AI(:,2) % MIN+1
      Write(unit,'(a,100(1x,I4))')'  jmyz: ',AI(:,2) % MAX-AI(:,2) % MIN+1
      Write(unit,*)'   kmyz: ',LM
      Call Free_File(unit)

      Deallocate(AI)

    End Subroutine write_layout

    Subroutine write_start_date(clock, unit)
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

  End Subroutine Create_Restart_File


  function PKZ_GISS() Result(PKZ)
    USE RESOLUTION, only: IM, LM, LS1
    Use MODEL_COM, only : SIG, SIGE, Ptop, PSFMPT, P
    use DOMAIN_DECOMP, only: grid, get
    USE CONSTANT, only: KAPA

    REAL*8 :: PE
    REAL*8 :: PKZ(IM,grid % J_STRT:grid % J_STOP,LM)

    INTEGER :: L, j_0, j_1

    call get(grid, J_STRT=J_0, J_STOP=J_1)

    Do L = 1, LM

       If (L < LS1) THEN
          PKZ(:,:,L) = (SIG(L)*P(:,J_0:J_1) + Ptop) ** KAPA
       Else
          PKZ(:,:,L) = (SIG(L)*PSFMPT + Ptop) ** KAPA
       End IF

    End Do

  end function PKZ_GISS

  ! Convert Potential Temperature into (dry) Temperature
  function DryTemp_GISS() Result(T_dry)
    USE RESOLUTION, only: IM, LM, LS1
    Use MODEL_COM, only: T
    USE DOMAIN_DECOMP, only: grid, GET

    REAL*8 :: T_dry(IM,grid % J_STRT:grid % J_STOP,LM)
    REAL*8 :: PKZ(IM,grid % J_STRT:grid % J_STOP,LM)

    INTEGER :: J_0,J_1
    Call Get(grid, J_STRT=J_0, J_STOP=J_1)

    PKZ = PKZ_GISS()
    T_dry = PKZ * T(:,J_0:J_1,:)

  end function DryTemp_GISS

  Subroutine Copy_modelE_to_FV_import(fv)
    USE MODEL_COM, only:  Q     ! Secific Humidity
    USE MODEL_COM, only:  ZATMO ! Geopotential Height?
    USE MODEL_COM, Only : U, V, T
    Use DOMAIN_DECOMP, Only: grid, Get
    Type (FV_CORE) :: fv

    Integer :: nq

    Integer :: j_0, j_1

    Call Get(grid, j_strt=j_0, j_stop=j_1)

    ! 1) Link/copy modelE data to import state
    fv % phis=ZATMO(:,j_0:j_1)
#ifdef NO_FORCING
    fv % phis = 0
#endif

    ! Moisture
    fv % Q = Reverse(Q(:,j_0:j_1,:))
#ifdef NO_FORCING
    fv % Q = 0
#endif

  End Subroutine Copy_modelE_to_FV_import


  Subroutine Copy_FV_export_to_modelE(fv)
    Use Resolution, only: IM,JM,LM,LS1
    USE DYNAMICS, ONLY: PUA,PVA,SDA, PU, PV, SD
    Use MODEL_COM, only: U, V, T, P, PSFMPT, Q
    Use MODEL_COM, only : SIGE, Ptop, P
    USE DOMAIN_DECOMP, only: grid, GET, SOUTH, NORTH, HALO_UPDATE
    Use ATMDYN, only: CALC_AMPK
    USE GEOM
    Type (FV_CORE) :: fv
    real*4, Dimension(:,:,:), Pointer :: U_a, V_a, T_fv, PLE

    Integer :: unit

    Integer :: i,j,k
    Integer :: i_0, i_1, j_0, j_1

    Call Get(grid, i_strt=i_0, i_stop=i_1, j_strt=j_0, j_stop=j_1)

    ! First compute updated values for modelE.  Then capture the
    ! new state in fv % *_old for computing tendencies.
    ! ----------------------------------------------------------

    ! Velocity field
    !---------------
    call ESMFL_StateGetPointerToData ( fv % export,U_a,'U',rc=rc)
      VERIFY_(rc)
    call ESMFL_StateGetPointerToData ( fv % export,V_a,'V',rc=rc)
      VERIFY_(rc)
    call Regrid_A_to_B(Reverse(U_a), Reverse(V_a), U(:,J_0:J_1,:), V(:,J_0:J_1,:))
    fv % U_old = U(:,J_0:J_1,:)
    fv % V_old = V(:,J_0:J_1,:)

    ! Potential temperature (save dry Temperature for computing tendencies)
    !----------------------------------------------------------------------
    call ESMFL_StateGetPointerToData ( fv % export,T_fv,'TH',rc=rc)
      VERIFY_(rc)
    call ConvertPotTemp_FV2GISS(T_fv, T)

    ! Use edge pressure export to update surface pressure
    !----------------------------------------------------------------------
    call ESMFL_StateGetPointerToData ( fv % export,PLE,'PLE',rc=rc)
      VERIFY_(rc)

    call ConvertPressure_FV2GISS(PLE, fv % PE_old)
    ! Just need surface pressure - Ptop
    P(:,J_0:J_1) = (fv % PE_old(:,:,1) - Ptop)/SIGE(1)
    CALL CALC_AMPK(LS1-1)

    ! Preserve state information for later computation of tendencies.
    fv % dT_old = DryTemp_GISS()
    fv % dPT_old = DeltPressure_DryTemp_GISS()

  End Subroutine Copy_FV_export_to_modelE

  subroutine ConvertPressure_GISS2FV_r4(P_giss, P_fv)
    Real*8, intent(in) :: P_giss(:,:,:)
    Real*4, intent(out) :: P_fv(:,:,:)   ! no halo in this case

    P_fv = Reverse(P_giss * PRESSURE_UNIT_RATIO)

  end Subroutine ConvertPressure_GISS2FV_r4

  subroutine ConvertPressure_GISS2FV_r8(P_giss, P_fv)
    Real*8, intent(in) :: P_giss(:,:,:)
    Real*8, intent(out) :: P_fv(:,:,:)   ! no halo in this case

    P_fv = Reverse(P_giss * PRESSURE_UNIT_RATIO)

  end Subroutine ConvertPressure_GISS2FV_r8

  ! Convert pressure from GISS representation to FV representation.
  ! Both the order of levels and the units must be adjusted.
  subroutine ConvertPressure_FV2GISS(P_fv, P_giss)
    Real*4, intent(in) :: P_fv(:,:,:)
    Real*8, intent(out) :: P_giss(:,:,:)

    P_giss = Reverse(P_fv / PRESSURE_UNIT_RATIO)

  end Subroutine ConvertPressure_FV2GISS

  ! Convert potential temperature between the two representations.
  subroutine CnvPotTemp_GISS2FV_r8(PT_giss, PT_fv)
    Use CONSTANT, only : KAPA
    Real*8, intent(in) :: PT_giss(:,:,:)
    Real*8, intent(out) :: PT_fv(:,:,:)

    PT_fv = Reverse(PT_giss * REF_RATIO ** KAPA)

  end Subroutine CnvPotTemp_GISS2FV_r8

  ! Convert potential temperature between the two representations.
  subroutine CnvPotTemp_GISS2FV_r4(PT_giss, PT_fv)
    Use CONSTANT, only : KAPA
    Real*8, intent(in) :: PT_giss(:,:,:)
    Real*4, intent(out) :: PT_fv(:,:,:)

    PT_fv = Reverse(PT_giss * REF_RATIO ** KAPA)

  end Subroutine CnvPotTemp_GISS2FV_r4

  !------------------------------------------------------------------------
  ! Convert potential temperature as exported by FV to the corresponding
  ! potential temperauture within modelE.  Note that FV export uses
  ! a reference pressure of 10^5 pa.
  !------------------------------------------------------------------------
  subroutine CnvPotTemp_FV2GISS_r4(PT_fv, PT_giss)
    Use CONSTANT, only : KAPA
    Real*4, intent(in) :: PT_fv(:,:,:)
    Real*8, intent(out) :: PT_giss(:,:,:)

    ! As an export, FV provides Pot Temp with a reference
    ! pressure of 10^5 Pa, whereas GISS uses 100 Pa (1 mb)
    REAL*8, PARAMETER :: REF_RATIO = 100000 / 100
    integer :: n

    ! PT_GISS has a halo, while PT_fv does not.
    n = size(PT_GISS,2) ! exclude halo
    PT_giss(:,2:n-1,:) = Reverse(PT_fv / REF_RATIO ** KAPA)

  end Subroutine CnvPotTemp_FV2GISS_r4

  !------------------------------------------------------------------------
  ! Given velocities from modelE (B grid, k=1 bottom), produce
  ! velocities for import to FV (A grid, k=1 top)
  !------------------------------------------------------------------------
  subroutine ConvertUV_GISS2FV(U_giss, V_giss, U_fv, V_fv)
    Real*8, intent(in) :: U_giss(:,:,:)
    Real*8, intent(in) :: V_giss(:,:,:)
    Real*4, intent(out) :: U_fv(:,:,:)
    Real*4, intent(out) :: V_fv(:,:,:)

    Call Regrid_B_to_A(Reverse(U_giss), Reverse(V_giss), U_fv, V_fv)

  end subroutine ConvertUV_GISS2FV

  !------------------------------------------------------------------------
  ! This following routine interpolates U, V from the Arakawa A grid
  ! to the Arakawa B grid.  B-grid velocities correspond to GISS modelE,
  ! while A-grid velocities are used for import/export of the FV component.
  !------------------------------------------------------------------------
  subroutine Regrid_A_to_B(U_a, V_a, U_b, V_b)
    Use Resolution, only : IM, LM
    Use DOMAIN_DECOMP, only : grid, get, SOUTH, HALO_UPDATE
    Real*4, intent(in), Dimension(:,grid % J_STRT:,:) :: U_a, V_a
    Real*8, intent(out),Dimension(:,grid % J_STRT:,:) :: U_b, V_b

    Real*8, allocatable, dimension(:,:,:) :: Ua_halo, Va_halo
    Integer :: i,j,k,im1
    integer :: j_0stgr, j_1stgr,J_0h,J_1h,J_0,J_1

    Call Get(grid, J_STRT_STGR=j_0stgr, J_STOP_STGR=j_1stgr, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H, &
         & J_STRT=J_0, J_STOP=J_1)
    Allocate(Ua_halo(IM,J_0h:J_1h,LM), Va_halo(IM,J_0h:J_1h,LM))

    Ua_halo(:,J_0:J_1,:) = U_a
    Va_halo(:,J_0:J_1,:) = V_a
    Call Halo_Update(grid, Ua_halo,FROM=SOUTH)
    Call Halo_Update(grid, Va_halo,FROM=SOUTH)

    Do k = 1, LM
       Do j = j_0STGR, j_1STGR
          im1 = IM
          Do i = 1, IM

             u_b(i,j,k) = (Ua_halo(im1,j-1,k) + Ua_halo(i,j-1,k) + Ua_halo(im1,j,k) + Ua_halo(i,j,k))/4
             v_b(i,j,k) = (Va_halo(im1,j-1,k) + Va_halo(i,j-1,k) + Va_halo(im1,j,k) + Va_halo(i,j,k))/4
             im1 = i
 
          End do
       end do
    end do

    Deallocate(Ua_halo, Va_halo)

  end subroutine Regrid_A_to_B

  !------------------------------------------------------------------------
  ! This following routine interpolates U, V from the Arakawa B grid
  ! to the Arakawa A grid.  B-grid velocities correspond to GISS modelE,
  ! while A-grid velocities are used for import/export of the FV component.
  !------------------------------------------------------------------------
  subroutine Regrid_B_to_A(U_b, V_b, U_a, V_a)
    Use Resolution, only : IM, JM, LM
    Use DOMAIN_DECOMP, only : grid, get, NORTH, HALO_UPDATE
    Real*8, intent(in), Dimension(:,grid % J_STRT:,:) :: U_b, V_b
    Real*4, intent(out), Dimension(:,grid % J_STRT:,:) :: U_a, V_a

    Real*8, allocatable, dimension(:,:,:) :: Ub_halo, Vb_halo

    Integer :: i,j,k,im1
    integer :: j_0, j_1, j_0s, j_1s, j_0h, j_1h
    logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

    Call Get(grid, J_STRT=J_0, J_STOP=J_1, J_STRT_SKP=J_0s, J_STOP_SKP=J_1S, &
         & HAVE_SOUTH_POLE=HAVE_SOUTH_POLE, HAVE_NORTH_POLE=HAVE_NORTH_POLE, &
         & J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

    Allocate(Ub_halo(IM,J_0h:J_1h,LM), Vb_halo(IM,J_0h:J_1h,LM))

    Ub_halo(:,J_0:J_1,:) = U_b
    Vb_halo(:,J_0:J_1,:) = V_b

    Call Halo_Update(grid, Ub_halo,FROM=NORTH)
    Call Halo_Update(grid, Vb_halo,FROM=NORTH)

    Do k = 1, LM
       Do j = j_0s,j_1s
          im1 = IM
          Do i = 1, IM

             u_a(im1,j,k) = (ub_halo(im1,j,k) + ub_halo(i,j,k) + ub_halo(im1,j+1,k) + ub_halo(i,j+1,k))/4
             v_a(im1,j,k) = (vb_halo(im1,j,k) + vb_halo(i,j,k) + vb_halo(im1,j+1,k) + vb_halo(i,j+1,k))/4
             im1 = i

          End do
       end do
    end do

    deallocate(ub_halo, Vb_halo)

    ! Polar conditions are a bit more complicated
    ! First determine an "absolute" U, V at the pole, then use sin/cos to
    ! map to the indidual longitudes.

    If (HAVE_SOUTH_POLE) then
       Call FixPole(U_b(:,2,:), V_b(:,2,:), U_a(:,1,:), V_a(:,1,:))
    End If

    If (HAVE_NORTH_POLE) then
       Call FixPole(U_b(:,JM,:), V_b(:,JM,:), U_a(:,JM,:), V_a(:,JM,:))
    End If

  Contains
    ! This routine works for both poles.
    ! The sign is correct for the NP, but (-1) appears as a factor twice
    ! in the derivation for the south pole.
    !----------------------------------------
    Subroutine FixPole(Ub, Vb, Ua, Va)
      USE GEOM, only : SINIV, COSIV, SINIP, COSIP ! trig of lon and stgr lon
      Real*8, intent(in) :: Ub(IM, LM), Vb(IM,LM)
      Real*4, intent(out) :: Ua(IM, LM), Va(IM,LM)

      Real*8 :: us, vs
      Integer :: k

      do k = 1, LM

         us = Sum( -ub(:,k) * SINIP - vb(:,k) * COSIP) / IM
         vs = Sum(  ub(:,k) * COSIP - vb(:,k) * SINIP) / IM

         ua(:,k) = -us * SINIV + vs * COSIV
         va(:,k) = -us * COSIV - vs * SINIV

      end do

    End Subroutine FixPole
  end subroutine Regrid_B_to_A

  !--------------------------------------------------------------------
  ! This following routine interpolates U, V from the Arakawa B grid
  ! to the Arakawa D grid.  B-grid velocities correspond to GISS modelE,
  ! while D-grid velocities are used for the initial conditions of FV
  !--------------------------------------------------------------------
  Subroutine regrid_B_to_D(u_b, v_b, u_d, v_d)
    USE RESOLUTION, only: IM, JM, LM
    Use DOMAIN_DECOMP, only: grid, get
    Implicit None
    Real*8, intent(in) :: U_b(:,grid % j_strt_halo:,:)
    Real*8, intent(in) :: V_b(:,grid % j_strt_halo:,:)
    Real*8, intent(out) :: U_d(:,grid % j_strt:,:)
    Real*8, intent(out) :: V_d(:,grid % j_strt:,:)

    Integer :: i, j, k, ip1
    Integer :: j_0, j_1, j_0s, j_1s
    Logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

    Call Get(grid, j_strt=j_0, j_stop=j_1, j_strt_skp=j_0s, j_stop_skp=j_1s, &
         & HAVE_SOUTH_POLE=HAVE_SOUTH_POLE, HAVE_NORTH_POLE=HAVE_NORTH_POLE)

    If (HAVE_SOUTH_POLE) Then
       ! 1st interpolate to 'A' grid
       Call FixPole(U_b(:,2,:), V_b(:,2,:), U_d(:,1,:))
     ! V_d(:,1,:) = (V_b(:,2,:) + CSHIFT(V_b(:,2,:),1,1))/2
       U_d(:,1,:)=0 ! not used (but needs legal value)
    End If

    If (HAVE_NORTH_POLE) Then
       ! 1st interpolate to 'A' grid
       Call FixPole(U_b(:,JM,:), V_b(:,JM,:), U_d(:,JM,:))
     ! V_d(:,JM,:) = (V_b(:,JM,:) + CSHIFT(V_b(:,JM,:),1,1))/2
    End If

! V-grid
    Do k = 1, LM
       Do j = j_0s, j_1s
          Do i = 1, IM
             V_d(i,j,k) = (V_b(i,j,k) + V_b(i,j+1,k)) /2
          End Do
       !  V_d(:,j,:) = (V_b(:,j,:) + CSHIFT(V_b(:,j,:),1,1))/2
       End Do
    End Do
! U-grid
    Do k = 1, LM
       Do j = j_0s, j_1
          i = IM
          Do ip1 = 1, IM
             U_d(i,j,k) = (U_b(i,j,k) + U_b(ip1,j,k)) /2
             i = ip1
          End Do
       !  U_d(:,j,:) = (U_b(:,j,:) + CSHIFT(U_b(:,j,:),1,1))/2
       End Do
    End Do

  Contains

    Subroutine FixPole(Ub, Vb, Ud)
      USE GEOM, only : SINIP, COSIP ! trig of lon and stgr lon
      Real*8, intent(in) :: Ub(IM, LM), Vb(IM,LM)
      Real*8, intent(out) :: Ud(IM, LM)

      Real*8 :: US, VS
      Integer :: k

      do k = 1, LM

         us = Sum( -ub(:,k) * SINIP - vb(:,k) * COSIP) / IM
         vs = Sum(  ub(:,k) * COSIP - vb(:,k) * SINIP) / IM

         ud(:,k) = (ub(:,k) + (-us * SINIP + vs * COSIP))/2

      end do

    End Subroutine FixPole

  End Subroutine regrid_B_to_D

  ! Reverse the order of vertical levels.
  Function reverse_3d_r8(A) Result(B)
    Real*8, Intent(In) :: A(:,:,:)
    Real*8             :: B(Size(A,1),Size(A,2),size(A,3))

    Integer, parameter :: K_IDX = 3
    Integer :: k, n


    n = Size(A, K_IDX)

    Do k = 1, n
       B(:,:,k) = A(:,:,1+n-k)
    End Do

  End Function reverse_3d_r8

  ! Single precision variant of Reverse()
  Function reverse_3d_r4(A) Result(B)
    Real*4, Intent(In) :: A(:,:,:)
    Real*4             :: B(Size(A,1),Size(A,2),size(A,3))

    Integer, parameter :: K_IDX = 3
    Integer :: k, n

    n = Size(A, K_IDX)

    Do k = 1, n
       B(:,:,k) = A(:,:,1+n-k)
    End Do

  End Function reverse_3d_r4

  Subroutine Write_Profile(arr, name)
    Use RESOLUTION,    Only: IM, JM, LM
    Use DOMAIN_DECOMP, Only: grid, PACK_DATA, AM_I_ROOT
    Real*8, intent(in) :: arr(:,:,:)
    character(len=*), intent(in) :: name

    Integer :: k, km
    Real*8 :: rng(3,LM)
    Real*8 :: arr_global(IM,JM,size(arr,3))
    Real*8 :: arr_tmp(IM,grid % j_STRT_HALO:grid % J_stop_HALO,size(arr,3))

    arr_tmp(:,grid % J_strt:grid % J_STOP,:)=arr

    Call PACK_DATA(grid, arr_tmp, arr_global)

    IF (AM_I_ROOT()) Then
       rng(1,:) = MINVAL(MINVAL(arr_global,DIM=1),DIM=1)
       rng(2,:) = MAXVAL(MAXVAL(arr_global,DIM=1),DIM=1)
       rng(3,:) = SUM(SUM(arr_global,DIM=1),DIM=1)/(IM*JM)

       print*,'***********'
       print*,'stats for ',trim(name)
       km = size(arr,3)

       Do k = 1, km
          Write(*,'(a,i4.0,3(f21.9,1x))')'k:',k,rng(:,k)
       End Do
       print*,'***********'
       print*,' '
    End IF

  End Subroutine Write_Profile

  function compute_phi(P, T, SZ, zatmo) result(phi)
    USE CONSTANT, only: KAPA, BYKAPA, BYKAPAP1, BYKAPAP2, RGAS
    USE MODEL_COM, only: LS1, LM, DSIG, SIG, SIGE, PTOP, PSFMPT
    USE MODEL_COM, only: IM, JM, LM
    USE DOMAIN_DECOMP, Only: grid, Get
    USE GEOM, only: IMAXJ
    implicit none

    real*8, intent(in) :: P(:,grid % J_STRT_HALO:)
    real*8, intent(in) :: T(:,grid % J_STRT_HALO:,:)
    real*8, intent(in) :: SZ(:,grid % J_STRT_HALO:,:)
    real*8, intent(in) :: ZATMO(:,grid % J_STRT_HALO:)
    real*8 :: phi(IM,grid % J_STRT_HALO:grid % J_STOP_HALO,LM)

    REAL*8 :: PKE(LS1:LM+1)
    REAL*8 :: PIJ, PHIDN
    REAL*8 :: DP, BYDP, P0, TZBYDP, X
    REAL*8 :: PUP, PKUP, PKPUP, PKPPUP
    REAL*8 :: PDN, PKDN, PKPDN, PKPPDN
    INTEGER :: I, J, L

    INTEGER :: J_0, J_1
    LOGICAL :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE

    call get(grid, J_STRT=J_0, J_STOP=J_1, HAVE_NORTH_POLE=HAVE_NORTH_POLE, &
         & HAVE_SOUTH_POLE = HAVE_SOUTH_POLE)

    DO L=LS1,LM+1
       PKE(L)=(PSFMPT*SIGE(L)+PTOP)**KAPA
    END DO

!$OMP  PARALLEL DO PRIVATE(I,J,L,DP,P0,PIJ,PHIDN,TZBYDP,X,
!$OMP*             BYDP,PDN,PKDN,PKPDN,PKPPDN,PUP,PKUP,PKPUP,PKPPUP)
    DO J=J_0,J_1
       DO I=1,IMAXJ(J)

          PIJ=P(I,J)

          PDN=PIJ+PTOP
          PKDN=PDN**KAPA
          PHIDN=ZATMO(I,J)

          !**** LOOP OVER THE LAYERS
          DO L=1,LM
             PKPDN=PKDN*PDN
             PKPPDN=PKPDN*PDN
             IF(L.GE.LS1) THEN
                DP=DSIG(L)*PSFMPT
                BYDP=1./DP
                P0=SIG(L)*PSFMPT+PTOP
                TZBYDP=2.*SZ(I,J,L)*BYDP
                X=T(I,J,L)+TZBYDP*P0
                PUP=SIGE(L+1)*PSFMPT+PTOP
                PKUP=PKE(L+1)
                PKPUP=PKUP*PUP
                PKPPUP=PKPUP*PUP
             ELSE
                DP=DSIG(L)*PIJ
                BYDP=1./DP
                P0=SIG(L)*PIJ+PTOP
                TZBYDP=2.*SZ(I,J,L)*BYDP
                X=T(I,J,L)+TZBYDP*P0
                PUP=SIGE(L+1)*PIJ+PTOP
                PKUP=PUP**KAPA
                PKPUP=PKUP*PUP
                PKPPUP=PKPUP*PUP
             END IF
             !**** CALCULATE PHI, MASS WEIGHTED THROUGHOUT THE LAYER
             PHI(I,J,L)=PHIDN+RGAS*(X*PKDN*BYKAPA-TZBYDP*PKPDN*BYKAPAP1 &
                  &      -(X*(PKPDN-PKPUP)*BYKAPA-TZBYDP*(PKPPDN-PKPPUP)*BYKAPAP2) &
                  &      *BYDP*BYKAPAP1)
             !**** CALULATE PHI AT LAYER TOP (EQUAL TO BOTTOM OF NEXT LAYER)
             PHIDN=PHIDN+RGAS*(X*(PKDN-PKUP)*BYKAPA-TZBYDP*(PKPDN-PKPUP) &
                  &     *BYKAPAP1)
             PDN=PUP
             PKDN=PKUP
          END DO
       END DO
    END DO
    !$OMP END PARALLEL DO

    !**** SET POLAR VALUES FROM THOSE AT I=1
    IF (HAVE_SOUTH_POLE) THEN
       DO L=1,LM
          PHI(2:IM,1,L)=PHI(1,1,L)
       END DO
    END IF

    IF (HAVE_NORTH_POLE) THEN
       DO L=1,LM
          PHI(2:IM,JM,L)=PHI(1,JM,L)
       END DO
    END IF

  end function compute_phi

  subroutine clear_accumulated_mass_fluxes()
    USE DYNAMICS, ONLY: PUA,PVA,SDA

    PUA(:,:,:) = 0.0
    PVA(:,:,:) = 0.0
    SDA(:,:,:) = 0.0

  end subroutine clear_accumulated_mass_fluxes

  subroutine accumulate_mass_fluxes(fv)
    Use Resolution, only: IM,JM,LM,LS1
    USE DYNAMICS, ONLY: PUA,PVA,SDA
    USE DYNAMICS, ONLY: PU,PV,CONV,SD,PIT
    USE MODEL_COM, only: DTsrc,DT,DSIG
    USE DOMAIN_DECOMP, only: get, grid, NORTH, SOUTH, HALO_UPDATE
    Use Constant, only: radius,pi,grav
    implicit none
    type (FV_core) :: fv
    real*4, Dimension(:,:,:), Pointer :: PLE, mfx_X, mfx_Y, mfx_Z
    integer :: J_0, J_1
    integer :: J_0S, J_1S
    integer :: J_0H, J_1H
    integer :: i,im1,j,l
    integer :: rc
    logical :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE
    real*8 :: DTLF, DTfac
    real*8 :: area,dlon,dlat,LAT

    DTfac = DT

    Call Get(grid, j_strt=j_0, j_stop=j_1, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S, &
         & J_STRT_HALO=J_0H, J_STOP_HALO = J_1H, &
         & HAVE_NORTH_POLE = HAVE_NORTH_POLE, HAVE_SOUTH_POLE = HAVE_SOUTH_POLE)

    ! Horizontal and Vertical mass fluxes
    !---------------

    call ESMFL_StateGetPointerToData ( fv % export,mfx_X,'MFX_A',rc=rc)
    VERIFY_(rc)
    call ESMFL_StateGetPointerToData ( fv % export,mfx_Y,'MFY_A',rc=rc)
    VERIFY_(rc)
    call ESMFL_StateGetPointerToData ( fv % export,mfx_Z,'MFZ',rc=rc)
    VERIFY_(rc)
    call ESMFL_StateGetPointerToData ( fv % export,PLE,'PLE',rc=rc)
    VERIFY_(rc)
   
!     \item {\tt MFX}:       Mass-Weighted U-Wind on C-Grid (Pa m^2/s)
!     \item {\tt MFY}:       Mass-Weighted V-wind on C-Grid (Pa m^2/s)
!     \item {\tt MFZ}:       Vertical mass flux (kg/(m^2*s))
#ifdef NO_MASS_FLUX
    mfx_X = 0
    mfx_Y = 0
    mfx_Z = 0
#endif
    call Regrid_A_to_B(Reverse(mfx_X), Reverse(mfx_Y), PU(:,J_0:J_1,:), PV(:,J_0:J_1,:))
    PU = PU/PRESSURE_UNIT_RATIO
    PV = PV/PRESSURE_UNIT_RATIO
 
    mfx_Z = Reverse(mfx_Z)/PRESSURE_UNIT_RATIO
    dlon = 2.0*pi/im
    dlat = pi/(jm-1)
    do l=1,lm
       do j=j_0,j_1
          LAT = ((-pi/2.0) + (j-0.5)*dlat )
          do i=1,im
             area = radius**2 * dlon * COS(LAT) * 2 * sin(dlat/2)
             mfx_Z(i,j-j_0+1,l) = grav*area*mfx_Z(i,j-j_0+1,l) ! convert to (Pa m^2/s)
          enddo
       enddo
    enddo
    SD(:,J_0:J_1,1:LM-1) = (mfx_Z(:,:,1:LM-1)) ! SD only goes up to LM-1

    ! Surface Pressure tendency - vert integral of horizontal convergence
    PIT(:,J_0:J_1) = mfx_Z(:,:,1) + sum(SD(:,J_0:J_1,1:LM-1),3)

    ! Recopy into CONV to support prior usage
    CONV(:,J_0:J_1,1) = PIT(:,J_0:J_1)
    CONV(:,J_0:J_1,2:LM) = SD(:,J_0:J_1,1:LM-1)

    PUA(:,J_0:J_1,:) = PUA(:,J_0:J_1,:) + PU(:,J_0:J_1,:)*DTfac
    PVA(:,J_0:J_1,:) = PVA(:,J_0:J_1,:) + PV(:,J_0:J_1,:)*DTfac

    SDA(:,J_0:J_1,1:LM-1) = SDA(:,J_0:J_1,1:LM-1) + SD(:,J_0:J_1,1:LM-1)*DTfac

  end subroutine accumulate_mass_fluxes

  subroutine set_zonal_flow(U_d, V_d, j0, j1)
    use GEOM, only: cosp
    use DOMAIN_DECOMP, only: grid
    Real*8, intent(out) :: U_d(:,grid % j_strt:,:)
    Real*8, intent(out) :: V_d(:,grid % j_strt:,:)
    integer :: j0, j1
    integer :: j

    do j = j0, j1
       U_d(:,j,:) = cosp(j)
       V_d(:,j,:) = 0
    end do
  end subroutine set_zonal_flow

  subroutine printSnapshot(fv)
    use MODEL_COM, only: U, V, T, Q, P
    type (fv_core) :: fv
    integer, save :: counter = 0

    counter = counter + 1
    write(40,*)'Iteration: ',counter
    write(40,*)'U : ', U(4, 3, 3)
    write(40,*)'V : ', V(4, 3, 3)
    write(40,*)'T : ', T(4, 3, 3)
    write(40,*)'Q : ', Q(4, 3, 3)
    write(40,*)'P : ', P(4, 3)
    write(40,*)'dudt  : ', fv % dudt(4, 3, 3), sum(fv% dudt)
    write(40,*)'dvdt  : ', fv % dvdt(4, 3, 3), sum(fv% dvdt)
    write(40,*)'dtdt  : ', fv % dtdt(4, 3, 3), sum(fv% dtdt)
    write(40,*)'dpedt : ', fv % dpedt(4, 3, 14), sum(fv% dpedt)
    write(40,*)'Q     : ', fv % q(4, 3, 3)
    write(40,*)'PHIS  : ', fv % phis(4, 3)
    write(40,*)'*******'
    write(40,*)' '
  end subroutine printSnapshot

end module FV_INTERFACE_MOD

!----------------------------------------------------------------
! The following routine interfaces the VERIFY_ macro with the
! GISS termination routine.   The line number and return code are
! written into a buffer which is passed to stop_model().
!----------------------------------------------------------------
Subroutine abort_core(line,rc)
  Implicit None
  Integer, Intent(In) :: line
  Integer, Intent(In) :: rc

  Character(len=100) :: buf

  Write(buf,*)'FV core failure at line',line,'\n error code:', rc
  Call stop_model(Trim(buf),99)

End Subroutine abort_core

module dynamics_save

  implicit none
  private

  public :: initialize
  public :: save_dynamics_state
  public :: restore_dynamics_state
  public :: write_dynamics_state

  real*8, allocatable, dimension(:,:,:), save :: U_save, V_save, T_save
  real*8, allocatable, dimension(:,:), save :: P_save
  real*8, allocatable, dimension(:,:,:), save :: Q_save, WM_save
  real*8, allocatable, dimension(:,:,:), save :: gz_save, phi_save
  real*8, allocatable, dimension(:,:,:,:), save :: tmom_save
  real*8, allocatable, dimension(:,:,:), save :: pdsig_save

contains

  subroutine initialize()
    use resolution, only : IM, JM, LM
    use somtq_com, only : tmom, NMOM

    allocate(U_save(IM,1:JM,LM))
    allocate(V_save(IM,1:JM,LM))
    allocate(T_save(IM,1:JM,LM))
    allocate(P_save(IM,1:JM))
    allocate(Q_save(IM,1:JM,LM))
    allocate(WM_save(IM,1:JM,LM))

    allocate(gz_save(IM,1:JM,LM))
    allocate(phi_save(IM,1:JM,LM))
    allocate(tmom_save(NMOM,IM,1:JM,LM))
    allocate(pdsig_save(LM,IM,1:JM))

  end subroutine initialize

  subroutine save_dynamics_state()
    use resolution, only : IM, JM, LM
    use model_com, only: U, V, T, P, Q, WM
    use dynamics, only: gz, phi, pdsig
    use somtq_com, only : tmom

    U_save = U(:,1:JM,:)
    V_save = V(:,1:JM,:)
    T_save = T(:,1:JM,:)
    P_save = P(:,1:JM)
    Q_save = Q(:,1:JM,:)
    WM_save = WM(:,1:JM,:)

    gz_save = gz(:,1:JM,:)
    phi_save = phi(:,1:JM,:)
    tmom_save = tmom(:,:,1:JM,:)
    pdsig_save = pdsig(:,:,1:JM)

  end subroutine save_dynamics_state

  subroutine restore_dynamics_state()
    use resolution, only : IM, JM, LM
    use model_com, only: U, V, T, P, Q, WM
    use dynamics, only: gz, phi, pdsig
    use somtq_com, only : tmom

    U(:,1:JM,:) = U_save
    V(:,1:JM,:) = V_save
    T(:,1:JM,:) = T_save
    P(:,1:JM)   = P_save
    Q(:,1:JM,:) = Q_save
    WM(:,1:JM,:) = WM_save

    gz(:,1:JM,:)   = gz_save
    phi(:,1:JM,:)  = phi_save
    tmom(:,:,1:JM,:) = tmom_save
    pdsig(:,:,1:JM) = pdsig_save

  end subroutine restore_dynamics_state

  subroutine write_dynamics_state(dyn)
    use resolution, only : IM, JM, LM
    use model_com, only: U, V, T, P
    character(len=*), intent(in) :: dyn
    integer :: unit

    select case(dyn)
    case ('FV')
       unit = 21
    case ('MODELE A')
       unit = 22
    case ('MODELE B')
       unit = 23
    case ('MODELE C')
       unit = 24
    case default
    end select

    write(unit) 'U',U(:,1:JM,:)
    write(unit) 'V',V(:,1:JM,:)
    write(unit) 'T',T(:,1:JM,:)
    write(unit) 'P',P(:,1:JM)


    select case(dyn)
    case ('FV')
       unit = 31
    case ('MODELE A')
       unit = 32
    case ('MODELE B')
       unit = 33
    case ('MODELE C')
       unit = 34
    case default
    end select

    write(*,*) 'DYNAMICS for ' // trim(dyn)
    call write_surface_zone('U',U(:,:,1))
    call write_surface_zone('V',V(:,:,1))
    call write_surface_zone('T',T(:,:,1))
    call write_surface_zone('P',P(:,:))

    call write_avg('U',U)
    call write_avg('V',V)
    call write_avg('T',T)

  contains

    subroutine write_surface_zone(name, array)
      character(len=*), intent(in) :: name
      real*8, intent(in) :: array(:,:)
      integer :: j

      write(unit,*)' '
      write(unit,*)'*********************************'
      write(unit,*)'Surface Zonal Average for ' // trim(name)
      write(unit,*)'*********************************'
      do j = 1, JM
         write(unit,*) j, sum(array(:,j))/IM
      end do
      write(unit,*)' '
    end subroutine write_surface_zone

    subroutine write_avg(name, array)
      use geom, only: AREAG, dxyp
      character(len=*), intent(in) :: name
      real*8, intent(in) :: array(:,:,:)
      integer :: k

      write(unit,*)' '
      write(unit,*)'*********************************'
      write(unit,*)' Level Average for ' // trim(name)
      write(unit,*)'*********************************'
      do k = 1, LM
         write(unit,*) k, sum(sum(array(:,1:JM,k),1)*dxyp(1:jm))/AREAG
      end do
      write(unit,*)' '
    end subroutine write_avg

  end subroutine write_dynamics_state


end module dynamics_save
