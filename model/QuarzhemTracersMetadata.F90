!------------------------------------------------------------------------------
module QuarzhemTracersMetadata_mod
!------------------------------------------------------------------------------
!@sum  QuarzhemTracersMetadata_mod encapsulates the Quarzhem tracers metadata
!@auth NCCS ASTG
  use sharedTracersMetadata_mod
  use TRACER_COM, only: n_sil1quhe,n_sil2quhe,n_sil3quhe
  use TRACER_COM, only: set_ntsurfsrc
  use OldTracer_mod, only: oldAddTracer
  use OldTracer_mod, only: set_tr_mm, set_ntm_power
  use OldTracer_mod, only: set_trpdens
  use OldTracer_mod, only: set_trradius
  use OldTracer_mod, only: set_tr_wd_TYPE
  use OldTracer_mod, only: set_isDust
  use OldTracer_mod, only: set_fq_aer
  use OldTracer_mod, only: nPart
  use tracers_dust,only : &
    DensityHematite, DensityQuartz, FreeFe, frHemaInQuarAggr, &
    pureByTotalHematite
  use tracers_dust,only : DensityHematite, DensityQuartz
  use RunTimeControls_mod, only: tracers_drydep
  use Tracer_mod, only: Tracer_type
  implicit none
  private

  public Quarzhem_initMetadata

  integer :: n ! class scoped temporary tracer index

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  subroutine Quarzhem_InitMetadata(pTracer)
!------------------------------------------------------------------------------
    class (Tracer), pointer :: pTracer

    call  Sil1quhe_setSpec('Sil1quhe')
    call  Sil2quhe_setSpec('Sil2quhe')
    call  Sil3quhe_setSpec('Sil3quhe')
   
!------------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------------

      subroutine Sil1quhe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil1quhe=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frHemaInQuarAggr) * DensityQuartz &
          + frHemaInQuarAggr * DensityHematite)
      if (tracers_drydep) call set_trradius(n, 1.47D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Sil1quhe_setSpec

      subroutine Sil2quhe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil2quhe=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frHemaInQuarAggr) * DensityQuartz &
          + frHemaInQuarAggr * DensityHematite)
      if (tracers_drydep) call set_trradius(n, 2.94D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Sil2quhe_setSpec

      subroutine Sil3quhe_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil3quhe=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, (1.d0 - frHemaInQuarAggr) * DensityQuartz &
          + frHemaInQuarAggr * DensityHematite)
      if (tracers_drydep) call set_trradius(n, 5.88D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
      end subroutine Sil3quhe_setSpec

  end subroutine Quarzhem_InitMetadata

end module QuarzhemTracersMetadata_mod



