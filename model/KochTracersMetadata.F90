!------------------------------------------------------------------------------
module KochTracersMetadata_mod
!------------------------------------------------------------------------------
!@sum  KochTracersMetadata_mod encapsulates the KOCH tracers metadata
!@auth NCCS ASTG
  use sharedTracersMetadata_mod, only: DMS_setspec, &
    SO2_setspec, H2O2_s_setspec
  use sharedTracersMetadata_mod, only: convert_HSTAR
  use OldTracer_mod, only: oldAddTracer
  use OldTracer_mod, only: nPart, nGAS
  use OldTracer_mod, only: set_tr_mm
  use OldTracer_mod, only: set_ntm_power
  use OldTracer_mod, only: set_trpdens
  use OldTracer_mod, only: set_trradius
  use OldTracer_mod, only: set_fq_aer
  use OldTracer_mod, only: set_tr_wd_type
  use OldTracer_mod, only: set_HSTAR
  use OldTracer_mod, only: set_tr_RKD
  use OldTracer_mod, only: set_tr_DHD
  use OldTracer_mod, only: tr_RKD 
  use OldTracer_mod, only: set_ntisurfsrc
  use OldTracer_mod, only: set_needtrs
  use TRACER_COM, only:  n_MSA, n_SO2,  n_SO4, n_DMS, n_seasalt1,  n_seasalt2, &
    n_BCII,  n_BCIA,  n_BCB, n_OCII,  n_OCIA,  n_OCB, n_H2O2_s
  use TRACER_COM, only: aer_int_yr
  use TRACER_COM, only: offline_dms_ss, offline_ss
  use TRACER_COM, only: set_ntsurfsrc
  use Dictionary_mod, only: sync_param
  use RunTimeControls_mod, only: tracers_drydep
  use RunTimeControls_mod, only: sulf_only_aerosols
  use RunTimeControls_mod, only: tracers_special_shindell
  use Tracer_mod, only: Tracer
#ifdef TRACERS_AEROSOLS_VBS
  use TRACERS_VBS, only: ivbs_m2,ivbs_m1,ivbs_m0,ivbs_p1,ivbs_p2,ivbs_p3,&
                         ivbs_p4,ivbs_p5,ivbs_p6
#endif /* TRACERS_AEROSOLS_VBS */

  implicit none
  private

  public Koch_initMetadata

  integer :: n ! class scoped temporary tracer index

!------------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  subroutine KOCH_InitMetadata(pTracer)
!------------------------------------------------------------------------------
    class (Tracer), pointer :: pTracer

    call  DMS_setSpec('DMS')
    call  MSA_setSpec('MSA')
    call  SO2_setSpec('SO2')
    call  SO4_setSpec('SO4')
    if (.not. tracers_special_shindell) then
      call  H2O2_s_setSpec('H2O2_s')
    end if
    if (.not. sulf_only_aerosols) then
      call  seasalt1_setSpec('seasalt1')
      call  seasalt2_setSpec('seasalt2')
      call  BCII_setSpec('BCII')
      call  BCIA_setSpec('BCIA')
      call  BCB_setSpec('BCB')
    end if
#ifdef TRACERS_AEROSOLS_VBS
      call  VBS_setSpec('vbsGm2', ivbs_m2,'igas')
      call  VBS_setSpec('vbsGm1', ivbs_m1,'igas')
      call  VBS_setSpec('vbsGz', ivbs_m0,'igas')
      call  VBS_setSpec('vbsGp1', ivbs_p1,'igas')
      call  VBS_setSpec('vbsGp2', ivbs_p2,'igas')
      call  VBS_setSpec('vbsGp3', ivbs_p3,'igas')
      call  VBS_setSpec('vbsGp4', ivbs_p4,'igas')
      call  VBS_setSpec('vbsGp5', ivbs_p5,'igas')
      call  VBS_setSpec('vbsGp6', ivbs_p6,'igas')

      call  VBS_setSpec('vbsAm2', ivbs_m2,'iaer')
      call  VBS_setSpec('vbsAm1', ivbs_m1,'iaer')
      call  VBS_setSpec('vbsAz', ivbs_m0,'iaer')
      call  VBS_setSpec('vbsAp1', ivbs_p1,'iaer')
      call  VBS_setSpec('vbsAp2', ivbs_p2,'iaer')
      call  VBS_setSpec('vbsAp3', ivbs_p3,'iaer')
      call  VBS_setSpec('vbsAp4', ivbs_p4,'iaer')
      call  VBS_setSpec('vbsAp5', ivbs_p5,'iaer')
      call  VBS_setSpec('vbsAp6', ivbs_p6,'iaer')
#else
      call  OCII_setSpec('OCII')   !Insoluble industrial organic mass
      call  OCIA_setSpec('OCIA')   !Aged industrial organic mass
      call  OCB_setSpec('OCB')     !Biomass organic mass
#endif /* TRACERS_AEROSOLS_VBS */
     
!------------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------------

    subroutine MSA_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_MSA = n
      call set_ntm_power(n, -13)
      call set_tr_mm(n, 96.d+0) !(H2O2 34;SO2 64)
      call set_trpdens(n, 1.7d3) !kg/m3 this is sulfate value
      call set_trradius(n, 5.d-7 ) !m (SO4 3;BC 1;OC 3)
      call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
    end subroutine MSA_setSpec

    subroutine SO4_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_SO4 = n 
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 96.d+0)
      call set_trpdens(n, 1.7d3) !kg/m3 this is sulfate value
      call set_trradius(n, 3.d-7 ) !m
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
    end subroutine SO4_setSpec

    subroutine seasalt1_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_seasalt1 = n
      call set_ntsurfsrc(n,  0) ! ocean bubbles
      call set_ntisurfsrc(n, 1)
      call set_ntm_power(n, -10)
      call set_tr_mm(n, 75.d0)  !Na x 3.256
      call set_trpdens(n, 2.2d3) !kg/m3 This is for non-hydrated
      call set_trradius(n, 4.4d-7 ) ! This is non-hydrated
      call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
    end subroutine seasalt1_setSpec

    subroutine seasalt2_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_seasalt2 = n
      call set_ntsurfsrc(n,  0) ! ocean bubbles
      call set_ntisurfsrc(n, 1)
      call set_ntm_power(n, -9)
      call set_tr_mm(n, 75.d0)  !Na x 3.256
      call set_trpdens(n, 2.2d3) !kg/m3 This is for non-hydrated
      call set_trradius(n, 5.0d-6) ! This is non-hydrated
      if (OFFLINE_DMS_SS.ne.1 .and. OFFLINE_SS.ne.1) then
        call set_trradius(n, 1.7d-6 ) ! This is non-hydrated
      end if
      call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
    end subroutine seasalt2_setSpec

    subroutine BCII_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_BCII = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, 1.3d3) !kg/m3
      call set_trradius(n, 1.d-7 ) !m
      call set_fq_aer(n, 0.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
    end subroutine BCII_setSpec

    subroutine BCIA_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_BCIA = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, 1.3d3) !kg/m3
      call set_trradius(n, 1.d-7 ) !m
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
    end subroutine BCIA_setSpec

    subroutine BCB_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_BCB = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, 1.3d3) !kg/m3
      call set_trradius(n, 1.d-7 ) !m
      call set_fq_aer(n, 0.6d0 ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
#ifdef DYNAMIC_BIOMASS_BURNING
      ! 12 below are the 12 VDATA veg types or Ent remapped to them,
      ! from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
      emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 8.2731990d-09,&
        1.7817767d-07, 8.5456378d-08, 2.5662467d-07, 3.3909114d-07,&
        1.1377826d-07, 2.9145593d-07, 0.0000000d+00, 0.0000000d+00,&
        0.0000000d+00, 0.0000000d+00/)
#endif
    end subroutine BCB_setSpec

#ifdef TRACERS_AEROSOLS_VBS
    subroutine VBS_setSpec(name, n, index, type) result(label)
      character(len=*), intent(in) :: name
      integer, intent(in) :: index
      character(len=4), intent(in) :: type

      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      label = n

      select case (type)
      case ('igas')
        vbs_tr%igas(index) = n
        call num_srf_sources(n,.false.)
      case ('iaer')
        vbs_tr%iaer(index) = n
        call num_srf_sources(n,.true.)
      end select

      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      select case(name)
      case ('vbsGm2', 'vbsGm1', 'vbsGz',  'vbsGp1', 'vbsGp2', &! VBS gas-phase
        'vbsGp3', 'vbsGp4', 'vbsGp5', 'vbsGp6')
        call set_tr_wd_type(n, ngas)
        call set_tr_RKD(n, 1.d4 / convert_HSTAR ) !Henry; from mole/(L atm) to mole/J
      if (tracers_drydep) call set_HSTAR(n, tr_RKD(n)*convert_HSTAR)
      case ('vbsAm2', 'vbsAm1', 'vbsAz',  'vbsAp1', 'vbsAp2', ! VBS aerosol-phase
        &      'vbsAp3', 'vbsAp4', 'vbsAp5', 'vbsAp6')
#ifdef DYNAMIC_BIOMASS_BURNING
        ! 12 below are the 12 VDATA veg types or Ent remapped to them,
        ! from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
        emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 6.4230818d-07,&
          1.6844633d-06, 8.7537586d-07, 2.5902559d-06, 4.3200689d-06,&
          2.6824284d-06, 3.9549395d-06, 0.0000000d+00, 0.0000000d+00,&
          0.0000000d+00, 0.0000000d+00/)
#endif
        call set_tr_wd_type(n, npart)
        call set_trpdens(n, 1.5d3) !kg/m3
        call set_trradius(n, 3.d-7 ) !m
        call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      end select
    end subroutine VBS_setSpec
#endif /* TRACERS_AEROSOLS_VBS */

    subroutine OCII_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_OCII = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7 ) !m
      call set_fq_aer(n, 0.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
    end subroutine OCII_setSpec

    subroutine OCIA_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_OCIA = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7 ) !m
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
    end subroutine OCIA_setSpec

    subroutine OCB_setSpec(name)
      use OldTracer_mod, only: om2oc, set_om2oc
      character(len=*), intent(in) :: name
      real*8 :: tmp
      n = oldAddTracer(name)
      n_OCB = n
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      tmp = om2oc(n_OCB)
      call sync_param("OCB_om2oc",tmp)
      call set_om2oc(n_OCB, tmp)
#endif
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 15.6d0)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7 ) !m
      call set_fq_aer(n, 0.8d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
#ifdef DYNAMIC_BIOMASS_BURNING
      ! 12 below are the 12 VDATA veg types or Ent remapped to them,
      ! from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
      emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 6.4230818d-07,&
        1.6844633d-06, 8.7537586d-07, 2.5902559d-06, 4.3200689d-06,&
        2.6824284d-06, 3.9549395d-06, 0.0000000d+00, 0.0000000d+00,&
        0.0000000d+00, 0.0000000d+00/)
#endif
    end subroutine OCB_setSpec

  end subroutine KOCH_InitMetadata

end module KochTracersMetadata_mod



