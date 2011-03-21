#define VERIFY_(rc) If (rc /= ESMF_SUCCESS) Call abort_core(__LINE__,rc)

module FV_UTILS

  USE ESMF_MOD
  USE CONSTANT, only: KAPA

  implicit none
  private

  public :: FV_CORE
  public :: abort_core
  public :: ReverseLevels
  public :: ConvertPotTemp_GISS2FV
  public :: ConvertPotTemp_FV2GISS
  public :: ConvertPressure_GISS2FV
  public :: ConvertPressure_FV2GISS
  public :: ConvertUV_GISS2FV
  public :: EdgePressure_GISS
  public :: DryTemp_GISS
  public :: DeltPressure_DryTemp_GISS
  public :: DeltPressure_GISS
  public :: compute_phi
  public :: MoveFile
  public :: allocateFvExport3D
  public :: DumpState
  public :: load_configuration
  public :: Copy_modelE_to_FV_import
  public :: Copy_FV_export_to_modelE
  public :: Tendency
  public :: ClearTendencies
  public :: SaveTendencies
  public :: clear_accumulated_mass_fluxes
  public :: writeFVCSstate
  public :: StateWriteToFile

  Interface ReverseLevels
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

  ! Public parameters
  character(len=*), parameter, public :: FVCORE_INTERNAL_RESTART = 'fvcore_internal_rst'
  character(len=*), parameter, public :: FVCORE_IMPORT_RESTART   = 'fvcore_import_rst'
  character(len=*), parameter, public :: FVCORE_LAYOUT           = 'fvcore_layout.rc'

  character(len=*), parameter, public :: TENDENCIES_FILE = 'tendencies_checkpoint'

  Real*8, parameter, public :: PRESSURE_UNIT_GISS  =  100 ! 1 mb
  Real*8, parameter, public :: PRESSURE_UNIT_FV    =    1 ! 1 pa
  Real*8, parameter, public :: PRESSURE_UNIT_RATIO = PRESSURE_UNIT_GISS/PRESSURE_UNIT_FV

  integer, parameter, public :: INITIAL_START = 2 ! ISTART=2: cold start, no FV restart files
  integer, parameter, public :: EXTEND_RUN = 3    ! ISTART>2: FV restart files are read


  ! This data structure is o convenient entity for storing persistent data between
  ! calls to this module.  In addition to
  Type FV_CORE


     type (ESMF_gridcomp) :: gc   ! This is the handle for the fv dynamical core

     type (ESMF_grid)     :: grid ! Although modelE is not an ESMF component, it does have an ESMF_Grid
     type (ESMF_vm)       :: vm   ! Should be eliminated ... Only used to obtain NPES for config file.

     ! Import and Export states for FV dycore
     type(ESMF_state) :: import   ! Allocated within FV component
     type(ESMF_state) :: export   ! Allocated within FV component

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

  ! The following parameters address the fact that FV and modelE use
  ! different units and different reference pressures for potential temperature.
  ! Superficially there is redundancy between the two sets, but in some sense
  ! this is merely coincidental.
  Real*8, parameter :: REF_PRESSURE_GISS = 100 ! 1 mb = 100 pa
  Real*8, parameter :: REF_PRESSURE_FV   =   1 ! pa
  Real*8, parameter :: REF_RATIO = REF_PRESSURE_FV / REF_PRESSURE_GISS

contains

  !----------------------------------------------------------------
  ! The following routine interfaces the VERIFY_ macro with the
  ! GISS termination routine.   The line number and return code are
  ! written into a buffer which is passed to stop_model().
  !----------------------------------------------------------------
!-------------------------------------------------------------------------------
  Subroutine abort_core(line,rc)
!-------------------------------------------------------------------------------
    Implicit None
    Integer, Intent(In) :: line
    Integer, Intent(In) :: rc

    Character(len=100) :: buf

    Write(buf,*)'FV core failure at line',line,'\n error code:', rc
    Call stop_model(Trim(buf),99)

  End Subroutine abort_core

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

  function EdgePressure_GISS() Result(PE)
    USE RESOLUTION, only: IM, LM, LS1
    Use MODEL_COM, only : SIG, SIGE, Ptop, PSFMPT, P
    use domain_decomp_atm, only: grid, get

    REAL*8 :: PE(grid%I_STRT:grid%I_STOP,grid%J_STRT:grid%J_STOP,LM+1)

    INTEGER :: L, i_0, i_1, j_0, j_1

    call get(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1)

    Do L = 1, LM+1

       If (L < LS1) THEN
          PE(:,:,L) = SIGE(L)*P(I_0:I_1,J_0:J_1) + Ptop
       Else
          PE(:,:,L)  = SIGE(L)*PSFMPT + Ptop
       End IF

    End Do

  end function EdgePressure_GISS

  ! Compute Delta-pressure for GISS model
  function DeltPressure_GISS() Result(dP)
    USE RESOLUTION, only: IM, LM, LS1
    Use MODEL_COM, only: T
    USE DOMAIN_DECOMP_ATM, only: grid, GET

    REAL*8 :: dP(grid%I_STRT:grid%I_STOP,grid%J_STRT:grid%J_STOP,LM)
    REAL*8 :: PE(grid%I_STRT:grid%I_STOP,grid%J_STRT:grid%J_STOP,LM+1)

    INTEGER :: k

    PE = EdgePressure_GISS()
    do k = 1, LM
       dP(:,:,k) = (PE(:,:,k)-PE(:,:,k+1))
    end do

  end function DeltPressure_GISS

  ! Convert Potential Temperature into (dry) Temperature
  function DeltPressure_DryTemp_GISS() Result(dPT)
    USE RESOLUTION, only: IM, LM, LS1
    Use MODEL_COM, only: T
    USE DOMAIN_DECOMP_ATM, only: grid, GET

    REAL*8 :: dPT(grid%I_STRT:grid%I_STOP,grid%J_STRT:grid%J_STOP,LM)
    REAL*8 :: PE(grid%I_STRT:grid%I_STOP,grid%J_STRT:grid%J_STOP,LM+1)
    REAL*8 :: T_dry(grid%I_STRT:grid%I_STOP,grid%J_STRT:grid%J_STOP,LM)

    INTEGER :: J_0,J_1,k
    Call Get(grid, J_STRT=J_0, J_STOP=J_1)

    T_dry = DryTemp_GISS()
    PE = EdgePressure_GISS()
    do k = 1, LM
       dPT(:,:,k) = (PE(:,:,k)-PE(:,:,k+1)) * T_dry(:,:,k)
    end do

  end function DeltPressure_DryTemp_GISS

  ! Convert Potential Temperature into (dry) Temperature
  function DryTemp_GISS() Result(T_dry)
    USE RESOLUTION, only: IM, LM, LS1
    Use MODEL_COM, only: T
    USE DOMAIN_DECOMP_ATM, only: grid, GET

    REAL*8 :: T_dry(grid%I_STRT:grid%I_STOP,grid%J_STRT:grid%J_STOP,LM)
    REAL*8 :: PKZ(grid%I_STRT:grid%I_STOP,grid%J_STRT:grid%J_STOP,LM)

    INTEGER :: I_0, I_1, J_0,J_1
    Call Get(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1)

    PKZ = PKZ_GISS()
    T_dry = PKZ * T(I_0:I_1,J_0:J_1,:)

  end function DryTemp_GISS

  function PKZ_GISS() Result(PKZ)
    USE RESOLUTION, only: IM, LM, LS1
    Use MODEL_COM, only : SIG, Ptop, PSFMPT, P
    use domain_decomp_atm, only: grid, get

    REAL*8 :: PKZ(grid%I_STRT:grid%I_STOP,grid%J_STRT:grid%J_STOP,LM)

    INTEGER :: L
    INTEGER :: I_0, I_1, J_0,J_1
    Call Get(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1)

    Do L = 1, LM

       If (L < LS1) THEN
          PKZ(:,:,L) = (SIG(L)*P(I_0:I_1,J_0:J_1) + Ptop) ** KAPA
       Else
          PKZ(:,:,L) = (SIG(L)*PSFMPT + Ptop) ** KAPA
       End IF

    End Do
  end function PKZ_GISS

!-------------------------------------------------------------------------------
  subroutine ConvertUV_GISS2FV(U_orig, V_orig, U_d, V_d)
!-------------------------------------------------------------------------------
! Cubed-sphere modelE works with native-grid winds, so this routine is trivial.
    Use Resolution, only : LM
    Use Domain_decomp_atm, only : grid, get
    Real*8, intent(in), Dimension(grid%I_STRT:,grid%J_STRT:,:) :: U_orig, V_orig
    Real*4, intent(out), Dimension(grid%I_STRT:,grid%J_STRT:,:) :: U_d, V_d

    Integer :: i,j,k
    integer :: I_0, I_1, j_0, j_1

    Call Get(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1)

    Do k = 1, LM
       Do j = j_0,j_1
          Do i = i_0,i_1
             U_d(i,j,k) = U_orig(i,j,k)
             V_d(i,j,k) = V_orig(i,j,k)
          End do
       end do
    end do
  end subroutine ConvertUV_GISS2FV

!-------------------------------------------------------------------------------
  subroutine ConvertPressure_GISS2FV_r4(P_giss, P_fv)
!-------------------------------------------------------------------------------
    Real*8, intent(in) :: P_giss(:,:,:)
    Real*4, intent(out) :: P_fv(:,:,:)   ! no halo in this case

    P_fv = ReverseLevels(P_giss * PRESSURE_UNIT_RATIO)

  end Subroutine ConvertPressure_GISS2FV_r4

!-------------------------------------------------------------------------------
  subroutine ConvertPressure_GISS2FV_r8(P_giss, P_fv)
!-------------------------------------------------------------------------------
    Real*8, intent(in) :: P_giss(:,:,:)
    Real*8, intent(out) :: P_fv(:,:,:)   ! no halo in this case

    P_fv = ReverseLevels(P_giss * PRESSURE_UNIT_RATIO)

  end Subroutine ConvertPressure_GISS2FV_r8

  ! Convert pressure from GISS representation to FV representation.
  ! Both the order of levels and the units must be adjusted.
!-------------------------------------------------------------------------------
  subroutine ConvertPressure_FV2GISS(P_fv, P_giss)
!-------------------------------------------------------------------------------
    Real*4, intent(in) :: P_fv(:,:,:)
    Real*8, intent(out) :: P_giss(:,:,:)

    P_giss = ReverseLevels(P_fv / PRESSURE_UNIT_RATIO)

  end Subroutine ConvertPressure_FV2GISS

  ! Convert potential temperature between the two representations.
!-------------------------------------------------------------------------------
  subroutine CnvPotTemp_GISS2FV_r8(PT_giss, PT_fv)
!-------------------------------------------------------------------------------
    Real*8, intent(in) :: PT_giss(:,:,:)
    Real*8, intent(out) :: PT_fv(:,:,:)

    PT_fv = ReverseLevels(PT_giss * REF_RATIO ** KAPA)

  end Subroutine CnvPotTemp_GISS2FV_r8

  ! Convert potential temperature between the two representations.
!-------------------------------------------------------------------------------
  subroutine CnvPotTemp_GISS2FV_r4(PT_giss, PT_fv)
!-------------------------------------------------------------------------------
    Real*8, intent(in) :: PT_giss(:,:,:)
    Real*4, intent(out) :: PT_fv(:,:,:)

    PT_fv = ReverseLevels(PT_giss * REF_RATIO ** KAPA)

  end Subroutine CnvPotTemp_GISS2FV_r4

  !------------------------------------------------------------------------
  ! Convert potential temperature as exported by FV to the corresponding
  ! potential temperauture within modelE.  Note that FV export uses
  ! a reference pressure of 10^5 pa.
  !------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  subroutine CnvPotTemp_FV2GISS_r4(PT_fv, PT_giss)
!-------------------------------------------------------------------------------
    Real*4, intent(in) :: PT_fv(:,:,:)
    Real*8, intent(out) :: PT_giss(:,:,:)

    ! As an export, FV provides Pot Temp with a reference
    ! pressure of 10^5 Pa, whereas GISS uses 100 Pa (1 mb)
    REAL*8, PARAMETER :: REF_RATIO = 100000 / 100
    integer :: m,n

    ! PT_GISS has a halo, while PT_fv does not.
    m = size(PT_GISS,1) ! exclude halo
    n = size(PT_GISS,2) ! exclude halo
    PT_giss(2:m-1,2:n-1,:) = ReverseLevels(PT_fv / REF_RATIO ** KAPA)

  end Subroutine CnvPotTemp_FV2GISS_r4

!-------------------------------------------------------------------------------
  function compute_phi(P, T, SZ, zatmo) result(phi)
!-------------------------------------------------------------------------------
    USE CONSTANT, only : rgas,bykapa,bykapap1,bykapap2
    USE MODEL_COM, only: LS1, LM, DSIG, SIG, SIGE, PTOP, PSFMPT
    USE MODEL_COM, only: IM, JM, LM
    USE DOMAIN_DECOMP_ATM, Only: grid, Get
    implicit none

    real*8, intent(in), dimension(grid%i_strt_halo:grid%i_stop_halo, &
                                  grid%j_strt_halo:grid%j_stop_halo) :: p,zatmo
    real*8, intent(in), dimension(grid%i_strt_halo:grid%i_stop_halo, &
                                  grid%j_strt_halo:grid%j_stop_halo,lm) :: t,sz
    real*8, dimension(grid%i_strt_halo:grid%i_stop_halo, &
                      grid%j_strt_halo:grid%j_stop_halo,lm) :: phi

    REAL*8 :: PKE(LS1:LM+1)
    REAL*8 :: PIJ, PHIDN
    REAL*8 :: DP, BYDP, P0, TZBYDP, X
    REAL*8 :: PUP, PKUP, PKPUP, PKPPUP
    REAL*8 :: PDN, PKDN, PKPDN, PKPPDN
    INTEGER :: I, J, L

    INTEGER :: I_0, I_1, J_0, J_1
    LOGICAL :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE

    call get(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1, &
         & HAVE_NORTH_POLE = HAVE_NORTH_POLE, &
         & HAVE_SOUTH_POLE = HAVE_SOUTH_POLE)

    DO L=LS1,LM+1
       PKE(L)=(PSFMPT*SIGE(L)+PTOP)**KAPA
    END DO

!$OMP*             BYDP,PDN,PKDN,PKPDN,PKPPDN,PUP,PKUP,PKPUP,PKPPUP)
    DO J=J_0,J_1

       DO I=I_0,I_1

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

!-------------------------------------------------------------------------------
  subroutine MoveFile(fv_fname, fv_dfname, suffix)
!-------------------------------------------------------------------------------
    character(len=*), intent(in) :: fv_fname, fv_dfname, suffix

    call system('mv ' // FVCORE_INTERNAL_RESTART // suffix // ' ' // trim(fv_fname) )
    call system('mv ' // FVCORE_IMPORT_RESTART // ' ' // trim(fv_dfname) )
  end subroutine MoveFile

!-------------------------------------------------------------------------------
  subroutine allocateFvExport3D ( state, name )
!-------------------------------------------------------------------------------
    use ESMFL_MOD, Only: ESMFL_StateGetPointerToData
    type(ESMF_State),  intent(INOUT) :: state
    character(len=*),  intent(IN   ) :: name

    real, pointer :: ptr(:,:,:)
    logical       :: alloc

    alloc = .true.
    call ESMFL_StateGetPointerToData ( state, ptr , name , alloc, rc=rc )
    VERIFY_(rc)

  end subroutine allocateFvExport3D

#ifdef CUBED_SPHERE
#define USE_MAPL
#endif

#ifdef USE_MAPL
!-------------------------------------------------------------------------------
  subroutine DumpState(fv, clock, fv_fname, fv_dfname, suffix, isFinalize)
!-------------------------------------------------------------------------------
    use MAPL_mod, only: MAPL_MetaComp
    use MAPL_mod, only: MAPL_GetObjectFromGC
    use MAPL_mod, only: MAPL_Get
    use MAPL_mod, only: MAPL_GetResource
    use domain_decomp_atm, only: AM_I_ROOT

    Type (FV_CORE),    intent(inout) :: fv
    type (ESMF_clock), intent(inout) :: clock
    character(len=*), intent(in) :: fv_fname, fv_dfname, suffix
    logical, intent(in) :: isFinalize

    type (MAPL_MetaComp), pointer :: internalState
    type (ESMF_State) :: ESMFInternalState
    integer :: hdr
    integer :: rc

    call SaveTendencies(fv, FVCORE_IMPORT_RESTART)

    if(isFinalize) then
       call ESMF_GridCompFinalize( fv%gc, importState=fv%import, exportState=fv%export, &
                                   clock=clock, rc=rc)
    else
       ! workaround for RecordPhase since modelE is not using MAPL interface 
       ! for alarms
       call MAPL_GetObjectFromGC( fv%gc, internalSTate)
       call MAPL_Get(internalState, internal_ESMF_state=ESMFInternalState, rc=rc)
       VERIFY_(rc)
       call MAPL_GetResource( internalState   , hdr,         &
                               default=0, &
                               LABEL="INTERNAL_HEADER:", &
                               RC=rc)
       VERIFY_(rc)
       call StateWriteToFile(ESMFInternalState, clock, &
            & FVCORE_INTERNAL_RESTART // suffix,       &
            & 'binary',internalState, hdr==1, rc=rc)
       VERIFY_(rc)
    endif

    ! Now move the file into a more useful name
    if (AM_I_ROOT()) then
      call MoveFile(fv_fname, fv_dfname, suffix)
    end if

  end subroutine DumpState
#else

!-------------------------------------------------------------------------------
  subroutine DumpState(fv, clock, fv_fname, fv_dfname, suffix, isFinalize)
!-------------------------------------------------------------------------------
    use GEOS_mod, only: RecordPhase=>GEOS_RecordPhae
    use domain_decomp_atm, only: AM_I_ROOT

    Type (FV_CORE),    intent(inout) :: fv
    type (ESMF_clock), intent(in) :: clock
    character(len=*), intent(in) :: fv_fname, fv_dfname, suffix
    logical, intent(in) :: isFinalize

    call SaveTendencies(fv, FVCORE_IMPORT_RESTART)

    if(isFinalize) then
       call ESMF_GridCompFinalize( fv%gc, fv%import, fv%export, clock, rc=rc)
       VERIFY_(rc)
    else
       call ESMF_GridCompFinalize( fv%gc, fv%import, fv%export, clock, &
            &  phase=RecordPhase, rc=rc)
       VERIFY_(rc)
    endif

    ! Now move the file into a more useful name
    if (AM_I_ROOT()) then
      call MoveFile(fv_fname, fv_dfname, suffix)
    end if

  end subroutine DumpState
#endif

!-------------------------------------------------------------------------------
  function load_configuration(config_file) result( config )
!-------------------------------------------------------------------------------
    use ESMF_mod
    use FILEMANAGER
    Use MODEL_COM,  only: DT
    character(len=*), parameter :: Iam="FV_INTERFACE::loadconfiguration"
    character(len=*), intent(in) :: config_file
    type (ESMF_config)           :: config

    integer :: iunit
    type (ESMF_VM) :: vm

    config = ESMF_configcreate(rc=rc)
    VERIFY_(rc)

    call openunit(config_file, iunit, qbin=.false., qold=.false.)
    write(iunit,*)'FVCORE_INTERNAL_CHECKPOINT_FILE:  ', FVCORE_INTERNAL_RESTART
    write(iunit,*)'FVCORE_INTERNAL_RESTART_FILE:     ', FVCORE_INTERNAL_RESTART
    write(iunit,*)'FVCORE_LAYOUT:                    ', FVCORE_LAYOUT
    write(iunit,*)'RUN_DT:                           ', DT
    call closeUnit(iunit)

    Call ESMF_VMGetGlobal(vm, rc)
    call ESMF_VMbarrier(vm, rc)
    call ESMF_configloadfile(config, config_file, rc=rc)
    VERIFY_(rc)

  end function load_configuration


!-------------------------------------------------------------------------------
  Subroutine Copy_FV_export_to_modelE(fv)
!-------------------------------------------------------------------------------
    use ESMFL_MOD, Only: ESMFL_StateGetPointerToData
    Use Resolution, only: IM,JM,LM,LS1
    USE DYNAMICS, ONLY: PUA,PVA,SDA, PU, PV, SD
    Use MODEL_COM, only: U, V, T, P, PSFMPT, Q
    Use MODEL_COM, only : Ptop, P
    USE DOMAIN_DECOMP_ATM, only: grid, GET
    USE GEOM
    
    Type (FV_CORE) :: fv
    real*4, Dimension(:,:,:), Pointer :: T_fv, PLE, U_d, V_d

    Integer :: unit

    Integer :: i,j,k
    Integer :: i_0, i_1, j_0, j_1

    Call Get(grid, i_strt=i_0, i_stop=i_1, j_strt=j_0, j_stop=j_1)

    ! First compute updated values for modelE.  Then capture the
    ! new state in fv%*_old for computing tendencies.
    ! ----------------------------------------------------------
    call ESMFL_StateGetPointerToData ( fv%export,U_d,'U_DGRID',rc=rc)
    VERIFY_(rc)
    call ESMFL_StateGetPointerToData ( fv%export,V_d,'V_DGRID',rc=rc)
    VERIFY_(rc)
    U(I_0:I_1,J_0:J_1,:) = ReverseLevels(U_d)
    V(I_0:I_1,J_0:J_1,:) = ReverseLevels(V_d)
    fv%U_old = U(I_0:I_1,J_0:J_1,:)
    fv%V_old = V(I_0:I_1,J_0:J_1,:)

    ! Potential temperature (save dry Temperature for computing tendencies)
    !----------------------------------------------------------------------
    call ESMFL_StateGetPointerToData ( fv%export,T_fv,'TH',rc=rc)
    VERIFY_(rc)
    call ConvertPotTemp_FV2GISS(T_fv, T)

    ! Use edge pressure export to update surface pressure
    !----------------------------------------------------------------------
    call ESMFL_StateGetPointerToData ( fv%export,PLE,'PLE',rc=rc)
    VERIFY_(rc)

    call ConvertPressure_FV2GISS(PLE, fv%PE_old)
    ! Just need surface pressure - Ptop
    P(I_0:I_1,J_0:J_1) = fv%PE_old(:,:,1) - Ptop
    CALL CALC_AMPK(LS1-1)

    ! Preserve state information for later computation of tendencies.
    fv%dT_old = DryTemp_GISS()
    fv%dPT_old = DeltPressure_DryTemp_GISS()

#if defined(USE_FV_Q)
    Q(I_0:I_1,j_0:j_1,:) = ReverseLevels(fv%Q)
#endif

  End Subroutine Copy_FV_export_to_modelE

!-------------------------------------------------------------------------------
  Subroutine Copy_modelE_to_FV_import(fv)
!-------------------------------------------------------------------------------
    USE MODEL_COM, only:  Q     ! Secific Humidity
    USE MODEL_COM, only:  ZATMO ! Geopotential Height?
    USE MODEL_COM, Only : U, V, T
    Use Domain_decomp_atm, Only: grid, Get
    Type (FV_CORE) :: fv

    Integer :: nq

    Integer :: I_0, I_1, j_0, j_1

    Call Get(grid, i_strt=I_0, i_stop=I_1, j_strt=j_0, j_stop=j_1)

    ! 1) Link/copy modelE data to import state
    fv%phis=ZATMO(I_0:I_1,j_0:j_1)
#ifdef NO_FORCING
    fv%phis = 0
#endif

    ! Moisture
#ifndef ADIABATIC
    fv%Q = ReverseLevels(Q(I_0:I_1,j_0:j_1,:))
#ifdef NO_FORCING
    fv%Q = 0
#endif
#endif

  End Subroutine Copy_modelE_to_FV_import

  !-----------
  ! d/dt A
  ! Note that FV needs tendency with the order of vertical levels reversed.
  !-----------
!-------------------------------------------------------------------------------
  Function Tendency(A, A_old) Result(tend)
!-------------------------------------------------------------------------------
    USE MODEL_COM, Only: DTsrc
    REAL*8, INTENT(IN) :: A(:,:,:)
    REAL*8, INTENT(IN) :: A_old(:,:,:)
    REAL*8 :: tend(size(a,1),size(a,2),size(a,3))

    tend = (A - A_old)/DTsrc

  End Function Tendency

!-------------------------------------------------------------------------------
  subroutine ClearTendencies(fv)
!-------------------------------------------------------------------------------
    Type (FV_CORE) :: fv

    fv%dudt = 0
    fv%dvdt = 0
    fv%dtdt = 0
    fv%dpedt = 0

  end subroutine ClearTendencies

!-------------------------------------------------------------------------------
  subroutine SaveTendencies(fv, fv_dfname)
!-------------------------------------------------------------------------------
    use domain_decomp_atm, only: AM_I_ROOT
    use FILEMANAGER

    type (fv_core),    intent(inout) :: fv
    integer :: iunit
    character(len=*), intent(in) :: fv_dfname

    call OpenUnit(trim(fv_dfname) , iunit, qbin=.true.)

    call SaveArr(iunit, fv%U_old)
    call SaveArr(iunit, fv%V_old)
    call SaveArr(iunit, fv%dPT_old)
    call SaveArr(iunit, fv%PE_old)
    call SaveArr(iunit, fv%dT_old)

    call closeunit(iunit)

  contains

    !---------------------------------------------------------------------------
    subroutine SaveArr(iunit, arr)
    !---------------------------------------------------------------------------
      use domain_decomp_atm, only: grid, dwrite8_parallel, get
      integer, intent(in) :: iunit
      real*8, intent(in) :: arr(:,:,:)
      real*8, allocatable :: padArr(:,:,:)
      integer :: I_0,  I_1,  J_0,  J_1
      integer :: I_0H, I_1H, J_0H, J_1H

      Call Get(grid, i_strt=I_0, i_stop=I_1, j_strt=J_0, j_stop=J_1, &
           & i_strt_halo = I_0H, i_stop_halo = I_1H, &
           & j_strt_halo = J_0H, j_stop_halo = J_1H)
      allocate(padArr(I_0H:I_1H,J_0H:J_1H, size(arr,3)))

      padArr(I_0:I_1,J_0:J_1,:) = arr(:,:,:)
      call dwrite8_parallel(grid, iunit, nameunit(iunit), padArr)

      deallocate(padArr)

    end subroutine SaveArr

  end subroutine SaveTendencies

!-------------------------------------------------------------------------------
  subroutine clear_accumulated_mass_fluxes()
!-------------------------------------------------------------------------------
    USE DYNAMICS, ONLY: PUA,PVA,SDA
    PUA(:,:,:) = 0.0
    PVA(:,:,:) = 0.0
    SDA(:,:,:) = 0.0
  end subroutine clear_accumulated_mass_fluxes

!-------------------------------------------------------------------------------
  subroutine writeFVCSstate
!-------------------------------------------------------------------------------
! This is a temporary routine to write the Cubed Sphere state to a file
! that can be read with grads (thus real*4 fields are written out).
! Unfortunately the fields are represented in the CS-grid and so a conversion
! to lat-lon may be desirable.
    Use MAPL_IOMod, only: GETFILE, Free_file, MAPL_VarWrite
    USE DOMAIN_DECOMP_ATM, only: grid
    USE MODEL_COM, Only : U, V, T, P, IM, JM, LM
    real*4 u4(im,jm,lm), v4(im,jm,lm)
    real*4 pt4(im,jm,lm), p4(im,jm)
    integer :: unit, rc

    u4=u; v4=v; pt4=t; p4=p
    unit = GetFile('fvout.dat', form="unformatted", rc=rc)
    VERIFY_(rc)
    Call MAPL_VarWrite(unit, grid%ESMF_GRID, p4)
    Call MAPL_VarWrite(unit, grid%ESMF_GRID, u4)
    Call MAPL_VarWrite(unit, grid%ESMF_GRID, v4)
    Call MAPL_VarWrite(unit, grid%ESMF_GRID, pt4)
    Call Free_File(unit)

  end subroutine writeFVCSstate

!-------------------------------------------------------------------------------
  subroutine StateWriteToFile(STATE,CLOCK,FILENAME,FILETYPE,MPL,HDR,RC)
!-------------------------------------------------------------------------------
    use MAPL_mod, only: MAPL_MetaComp
    Use MAPL_IOMod, only: GETFILE, Free_file, MAPL_VarWrite, Write_parallel
    USE RESOLUTION, only: IM, JM, LM
    type(ESMF_State),                 intent(INOUT) :: STATE
    type(ESMF_Clock),                 intent(IN   ) :: CLOCK
    character(len=*),                 intent(IN   ) :: FILENAME
    character(LEN=*),                 intent(IN   ) :: FILETYPE
    type(MAPL_MetaComp),              intent(INOUT) :: MPL
    logical,                          intent(IN   ) :: HDR
    integer, optional,                intent(  OUT) :: RC

    character(len=ESMF_MAXSTR), parameter :: IAm="StateWriteToFile"

    type (ESMF_StateItemType), pointer    :: ITEMTYPES(:)
    character(len=ESMF_MAXSTR ), pointer  :: ITEMNAMES(:)
    integer                               :: ITEMCOUNT
    integer                               :: UNIT
    integer                               :: I, J, K, L, N
    integer                               :: YYYY, MM, DD, H, M, S
    type(ESMF_Time)                       :: currentTime
    integer                               :: HEADER(6)

! Get information from state
!---------------------------

    call ESMF_StateGet(STATE,ITEMCOUNT=ITEMCOUNT,RC=rc)
    VERIFY_(rc)

    allocate(ITEMNAMES(ITEMCOUNT),STAT=rc)
    VERIFY_(rc)
    allocate(ITEMTYPES(ITEMCOUNT),STAT=rc)
    VERIFY_(rc)

    call ESMF_StateGet(STATE,ITEMNAMELIST=ITEMNAMES, &
         & STATEITEMTYPELIST=ITEMTYPES,RC=rc)
    VERIFY_(rc)

! Open file
!----------

    if (filetype == 'binary' .or. filetype == 'BINARY') then
       UNIT = GETFILE(FILENAME, form="unformatted", rc=rc)
       VERIFY_(rc)
    elseif(filetype=="formatted".or.filetype=="FORMATTED") then
       UNIT = GETFILE(FILENAME, form="formatted", rc=rc)
       VERIFY_(rc)
    else
       UNIT=0
    end if

! Write data
!-----------

    call ESMF_ClockGet (clock, currTime=currentTime, rc=rc)
    VERIFY_(rc)
    call ESMF_TimeGet(CurrentTime, &
      &   YY=YYYY, MM=MM, DD=DD,   &
      &   H=H, M=M, S=S, rc=rc)
    VERIFY_(rc)

    HEADER(1) = YYYY
    HEADER(2) = MM
    HEADER(3) = DD
    HEADER(4) = H
    HEADER(5) = M
    HEADER(6) = S
    
    call Write_Parallel(HEADER, UNIT, RC=rc)
    VERIFY_(rc)
    
    HEADER(1) = IM
    HEADER(2) = JM
    HEADER(3) = LM
    HEADER(4) = 0
    HEADER(5) = 0
    
    call Write_Parallel(HEADER(1:5), UNIT, RC=rc)
    VERIFY_(rc)
    
    if(UNIT/=0) then
       do J = 1, ITEMCOUNT
          call MAPL_VarWrite(UNIT=UNIT, STATE=STATE, NAME=ITEMNAMES(J), rc=rc)
          VERIFY_(rc)
       end do
       call FREE_FILE(UNIT)
    else
       rc = -1  ! not yet
       VERIFY_(rc)
    endif
    
    deallocate(ITEMNAMES) 
    deallocate(ITEMTYPES)

  end subroutine StateWriteToFile

end module FV_UTILS
