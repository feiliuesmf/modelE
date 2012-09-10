#include "rundeck_opts.h"
module newTracer_COM
  implicit none
  private

  public :: getTracerNames
  public :: MAXLEN_TRACER_NAME

  integer, parameter :: MAXLEN_TRACER_NAME = 8

contains

!**** Each tracer has a variable name and a unique index
!**** The chemistry species need to be declared first, until the
!**** do igas=1,ntm_chem instances get corrected.

  subroutine printTracerNames(tracerNames)
    use domain_decomp_atm, only: am_i_root

    character(len=MAXLEN_TRACER_NAME), allocatable :: tracerNames(:)
    integer :: ntm       ! # of tracers
    integer :: i

    NTM = size(tracerNames)
    
    if (am_i_root()) then
      do i=1,NTM
        write(6,*) 'TRACER',i,trim(tracerNames(i))
      end do
    end if
  end subroutine printTracerNames
    
 
!TODO in F2003 this should be a function
  subroutine getTracerNames(tracerNames)  

!TODO use "ONLY" clause to specify which controls are used.
    use RunTimeControls_mod
    character(len=MAXLEN_TRACER_NAME), allocatable :: tracerNames(:)

    allocate(tracerNames(0))

    if (tracers_special_shindell) call appendShindellTracers()

    if ((.not. tracers_amp) .and. tracers_water) call appendNames(['Water   '])

    if (tracers_special_o18) then
      call appendNames(['H2O18   ','HDO     '])
!!$      call appendNames(['HTO     ','H2O17   '])
    end if

    if (tracers_gasexch_ocean_cfc) call appendNames(['CFCn    '])

    if (tracers_gasexch_ocean_co2 .or. tracers_gasexch_land_co2) then
      call appendNames(['CO2n    '])
    end if

    if (tracers_special_lerner) then
      call appendNames([ &
           & 'SF6     ','Rn222   ','CO2     ','N2O     ', &
           & 'CFC11   ','14CO2   ','CH4     ','O3      ','SF6_c   '])
    end if

    if (tracers_aerosols_koch) call appendKochTracers()

    if (tracers_aerosols_ocean) call appendNames(['OCocean '])

    if (tracers_dust) then
      call appendNames(['Clay    ','Silt1   ','Silt2   ','Silt3   '])
      if (tracers_dust_silt4) call appendNames(['Silt4   '])
    end if

    if (tracers_nitrate) call appendNames(['NH3     ','NH4     ','NO3p    '])
    if (tracers_hetchem) then
      call appendNames(['SO4_d1  ','SO4_d2  ','SO4_d3  '])
      if (tracers_nitrate) call appendNames(['N_d1    ','N_d2    ','N_d3    '])
    end if

    if (tracers_cosmo) then
      if (tracers_radon) call appendNames(['Pb210   '])
      call appendNames(['Be7     ','Be10    '])
      if (tracers_radon) call appendNames(['Rn222   '])
    end if

    if (tracers_minerals) then
      call appendNames([ &
           &     'ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar', &
           &     'Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps', &
           &     'Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps', &
           &     'Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps'])
    end if

    if (tracers_quarzhem) call appendNames(['Sil1QuHe','Sil2QuHe','Sil3QuHe'])

    if (tracers_air .or. htap_like_diags) call appendNames(['Air     '])

    if (tracers_amp) call appendAMP_tracers()

    if (tracers_tomas) call appendTomasTracers()

    ! Print it all out so we can see what's what!
    call printTracerNames(tracerNames)

  contains

    subroutine appendShindellTracers()
      call appendNames([ &
           &    'Ox      ','NOx     ','ClOx    ','BrOx    ',            &
           &                          'N2O5    ','HNO3    ','H2O2    ', &
           &    'CH3OOH  ','HCHO    ','HO2NO2  ','CO      ','CH4     ', &
           &    'PAN     ','Isoprene','AlkylNit','Alkenes ','Paraffin' ] )

      if (tracers_terp) call appendNames([ 'Terpenes' ])

      if (tracers_aerosols_soa) then
        call appendNames(['isopp1g ','isopp1a ','isopp2g ','isopp2a '])
        if (tracers_terp) then
          call appendNames(['apinp1g ','apinp1a ','apinp2g ','apinp2a '])
        end if
      end if

      call appendNames([ &
           &                       'HCl     ','HOCl    ','ClONO2  ', &
           & 'HBr     ','HOBr    ','BrONO2  ','N2O     ','CFC     ' ])

      if (shindell_strat_extra) then
        if (accmip_like_diags) then
          call appendNames(['codirect','stratOx ','GLT     '])
        else
          call appendNames(['GLT     '])  ! used to also do Be7, Be10
        end if
      end if
    end subroutine appendShindellTracers

    subroutine appendKochTracers()
      call appendNames(['DMS     ','MSA     ','SO2     ','SO4     ','H2O2_s  '])
      if (.not. sulf_only_aerosols) then
        call appendNames(['seasalt1','seasalt2','BCII    ','BCIA    ','BCB     '])
        if (tracers_aerosols_vbs) then
          call appendNames([ &
               & 'vbsGm2  ','vbsGm1  ','vbsGz   ','vbsGp1  ','vbsGp2  ', &
               & 'vbsGp3  ','vbsGp4  ','vbsGp5  ','vbsGp6  ', &
               & 'vbsAm2  ','vbsAm1  ','vbsAz   ','vbsAp1  ','vbsAp2  ', &
               &    'vbsAp3  ','vbsAp4  ','vbsAp5  ','vbsAp6  '])
        else
          call appendNames(['OCII    ','OCIA    ','OCB     '])
        end if
      end if
    end subroutine appendKochTracers

    ! The order of the ntmAMP aerosols matters!!!
    subroutine appendAMP_tracers()
      call appendNames(['M_NO3   ','M_NH4   ','M_H2O   '])

      if ( tracers_amp_m1 .or. tracers_amp_m2 .or. tracers_amp_m3 &
           .or. tracers_amp_m5 .or. tracers_amp_m6 .or. tracers_amp_m7) then
        call appendNames(['M_AKK_SU','N_AKK_1 '])                   !AKK
      end if
      call appendNames([ &
           & 'M_ACC_SU','N_ACC_1 ', &                                !ACC
           & 'M_DD1_SU','M_DD1_DU','N_DD1_1 ', &                     !DD1  
           & 'M_DS1_SU','M_DS1_DU','N_DS1_1 '])                      !DS1
      if (tracers_amp_m1 .or. tracers_amp_m2 .or. tracers_amp_m3 &
           & .or. tracers_amp_m4) then
        call appendNames([ &
             & 'M_DD2_SU','M_DD2_DU','N_DD2_1 ', &                   !DD2
             & 'M_DS2_SU','M_DS2_DU','N_DS2_1 '])                    !DS2
      end if

      if (tracers_amp_m1 .or.  tracers_amp_m2 .or. tracers_amp_m3 &
           &  .or. tracers_amp_m5 .or.  tracers_amp_m6 .or. tracers_amp_m7) then
        call appendNames([ &
             &    'M_SSA_SU','M_SSA_SS', &                           !SSA
             &    'M_SSC_SS'])                                       !SSC  
      end if
      if ( tracers_amp_m4 .or. tracers_amp_m8) then
        call appendNames(['M_SSS_SU','M_SSS_SS'])                    !SSS
      end if

      call appendNames([ &
           & 'M_OCC_SU','M_OCC_OC','N_OCC_1 ', &                     !OCC
           & 'M_BC1_SU','M_BC1_BC','N_BC1_1 ', &                     !BC1
           & 'M_BC2_SU','M_BC2_BC','N_BC2_1 '])                      !BC2

      if ( tracers_amp_m1 .or. tracers_amp_m5) then
        call appendNames(['M_BC3_SU','M_BC3_BC','N_BC3_1 '])          !BC3
      end if
      if ( tracers_amp_m2 .or. tracers_amp_m6) then
        call appendNames(['M_OCS_SU','M_OCS_OC','N_OCS_1 '])          !OCS
      end if

      if ( tracers_amp_m1 .or. tracers_amp_m2 .or. tracers_amp_m6) then
        call appendNames(['M_DBC_SU','M_DBC_BC','M_DBC_DU','N_DBC_1 ']) !DBC
      end if
      if (tracers_amp_m1 .or. tracers_amp_m2 .or. tracers_amp_m3 &
           & .or. tracers_amp_m6 .or. tracers_amp_m7) then
        call appendNames(['M_BOC_SU','M_BOC_BC','M_BOC_OC','N_BOC_1 ']) !BOC
      end if

      if ( tracers_amp_m1 .or. tracers_amp_m2 .or. tracers_amp_m5 &
           & .or. TRACERS_AMP_M6) then
        call appendNames(['M_BCS_SU','M_BCS_BC','N_BCS_1 '])          !BCS
      end if

      call appendNames([ &
           &    'M_MXX_SU', &                                      !MXX
           &    'M_MXX_BC','M_MXX_OC','M_MXX_DU','M_MXX_SS','N_MXX_1 ', &
           &    'H2SO4   ','DMS     ','SO2     ','H2O2_s  ','NH3     '])
    end subroutine appendAMP_tracers

    subroutine appendTomasTracers()
      call appendNames([ &
           &    'H2SO4   ','DMS     ','SO2     ','SOAgas  ','H2O2_s  ', &
           &    'NH3     ','NH4     '])
           ! H202_s should be used if gas chemistry is off.  &
      if (tomas_12_10nm) then
      call appendNames([ &
           &    'ASO4__01','ASO4__02','ASO4__03','ASO4__04','ASO4__05', &
           &    'ASO4__06','ASO4__07','ASO4__08','ASO4__09','ASO4__10', &
           &    'ASO4__11','ASO4__12', &
           &    'ANACL_01','ANACL_02','ANACL_03','ANACL_04','ANACL_05', &
           &    'ANACL_06','ANACL_07','ANACL_08','ANACL_09','ANACL_10', &
           &    'ANACL_11','ANACL_12', &
           &    'AECOB_01','AECOB_02','AECOB_03','AECOB_04','AECOB_05', &
           &    'AECOB_06','AECOB_07','AECOB_08','AECOB_09','AECOB_10', &
           &    'AECOB_11','AECOB_12', &
           &    'AECIL_01','AECIL_02','AECIL_03','AECIL_04','AECIL_05', &
           &    'AECIL_06','AECIL_07','AECIL_08','AECIL_09','AECIL_10', &
           &    'AECIL_11','AECIL_12', &
           &    'AOCOB_01','AOCOB_02','AOCOB_03','AOCOB_04','AOCOB_05', &
           &    'AOCOB_06','AOCOB_07','AOCOB_08','AOCOB_09','AOCOB_10', &
           &    'AOCOB_11','AOCOB_12', &
           &    'AOCIL_01','AOCIL_02','AOCIL_03','AOCIL_04','AOCIL_05', &
           &    'AOCIL_06','AOCIL_07','AOCIL_08','AOCIL_09','AOCIL_10', &
           &    'AOCIL_11','AOCIL_12', &
           &    'ADUST_01','ADUST_02','ADUST_03','ADUST_04','ADUST_05', &
           &    'ADUST_06','ADUST_07','ADUST_08','ADUST_09','ADUST_10', &
           &    'ADUST_11','ADUST_12', &
           &    'ANUM__01','ANUM__02','ANUM__03','ANUM__04','ANUM__05', &
           &    'ANUM__06','ANUM__07','ANUM__08','ANUM__09','ANUM__10', &
           &    'ANUM__11','ANUM__12', &
           &    'AH2O__01','AH2O__02','AH2O__03','AH2O__04','AH2O__05', &
           &    'AH2O__06','AH2O__07','AH2O__08','AH2O__09','AH2O__10', &
           &    'AH2O__11','AH2O__12' &
           & ])
      end if
      if (tomas_12_3nm) then
      call appendNames([ &
           &    'ASO4__01','ASO4__02','ASO4__03','ASO4__04','ASO4__05', &
           &    'ASO4__06','ASO4__07','ASO4__08','ASO4__09','ASO4__10', &
           &    'ASO4__11','ASO4__12','ASO4__13','ASO4__14','ASO4__15', &
           &    'ANACL_01','ANACL_02','ANACL_03','ANACL_04','ANACL_05', &
           &    'ANACL_06','ANACL_07','ANACL_08','ANACL_09','ANACL_10', &
           &    'ANACL_11','ANACL_12','ANACL_13','ANACL_14','ANACL_15', &
           &    'AECOB_01','AECOB_02','AECOB_03','AECOB_04','AECOB_05', &
           &    'AECOB_06','AECOB_07','AECOB_08','AECOB_09','AECOB_10', &
           &    'AECOB_11','AECOB_12','AECOB_13','AECOB_14','AECOB_15', &
           &    'AECIL_01','AECIL_02','AECIL_03','AECIL_04','AECIL_05', &
           &    'AECIL_06','AECIL_07','AECIL_08','AECIL_09','AECIL_10', &
           &    'AECIL_11','AECIL_12','AECIL_13','AECIL_14','AECIL_15', &
           &    'AOCOB_01','AOCOB_02','AOCOB_03','AOCOB_04','AOCOB_05', &
           &    'AOCOB_06','AOCOB_07','AOCOB_08','AOCOB_09','AOCOB_10', &
           &    'AOCOB_11','AOCOB_12','AOCOB_13','AOCOB_14','AOCOB_15', &
           &    'AOCIL_01','AOCIL_02','AOCIL_03','AOCIL_04','AOCIL_05', &
           &    'AOCIL_06','AOCIL_07','AOCIL_08','AOCIL_09','AOCIL_10', &
           &    'AOCIL_11','AOCIL_12','AOCIL_13','AOCIL_14','AOCIL_15', &
           &    'ADUST_01','ADUST_02','ADUST_03','ADUST_04','ADUST_05', &
           &    'ADUST_06','ADUST_07','ADUST_08','ADUST_09','ADUST_10', &
           &    'ADUST_11','ADUST_12','ADUST_13','ADUST_14','ADUST_15', &
           &    'ANUM__01','ANUM__02','ANUM__03','ANUM__04','ANUM__05', &
           &    'ANUM__06','ANUM__07','ANUM__08','ANUM__09','ANUM__10', &
           &    'ANUM__11','ANUM__12','ANUM__13','ANUM__14','ANUM__15', &
           &    'AH2O__01','AH2O__02','AH2O__03','AH2O__04','AH2O__05', &
           &    'AH2O__06','AH2O__07','AH2O__08','AH2O__09','AH2O__10', &
           &    'AH2O__11','AH2O__12','AH2O__13','AH2O__14','AH2O__15'&
           & ])
      endif

    end subroutine appendTomasTracers

    subroutine appendNames(names)
      character(len=MAXLEN_TRACER_NAME), intent(in) :: names(:)
      
      character(len=MAXLEN_TRACER_NAME), allocatable :: tmpNames(:)
      integer :: n
      
      ! more elegant implementation available for true F2003
      n = size(tracerNames)
      
      allocate(tmpNames(n))
      tmpNames = tracerNames
      
      deallocate(tracerNames)
      allocate(tracerNames(n + size(names)))
      tracerNames(1:n) = tmpNames
      tracerNames(n+1:) = names
      
      deallocate(tmpNames)
      
    end subroutine appendNames

  end subroutine getTracerNames

end module newTracer_COM
