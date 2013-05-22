!------------------------------------------------------------------------------
module MineralsTracersMetadata_mod
  !------------------------------------------------------------------------------
!@sum  MineralsTracersMetadata_mod encapsulates the Minerals tracers metadata
!@auth NCCS ASTG
  use sharedTracersMetadata_mod
  use TRACER_COM, only: &
    n_clayilli,n_claykaol,n_claysmec,n_claycalc, &
    n_clayquar,n_sil1quar,n_sil1feld,n_sil1calc, &
    n_sil1hema,n_sil1gyps,n_sil2quar,n_sil2feld, &
    n_sil2calc,n_sil2hema,n_sil2gyps,n_sil3quar, &
    n_sil3feld,n_sil3calc,n_sil3hema,n_sil3gyps, &
    n_sil1quhe,n_sil2quhe,n_sil3quhe
  use TRACER_COM, only: set_ntsurfsrc
  use OldTracer_mod, only: oldAddTracer
  use OldTracer_mod, only: set_tr_mm, set_ntm_power
  use OldTracer_mod, only: set_trpdens
  use OldTracer_mod, only: set_trradius
  use OldTracer_mod, only: set_tr_wd_TYPE
  use OldTracer_mod, only: set_isDust
  use OldTracer_mod, only: set_fq_aer
  use OldTracer_mod, only: nPart
  use RunTimeControls_mod, only: tracers_drydep
  use tracers_dust,only : DensityHematite, DensityQuartz
  use Tracer_mod, only: Tracer_type

  implicit none
  private

  public Minerals_initMetadata

  integer :: n ! class scoped temporary tracer index

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  subroutine Minerals_InitMetadata(pTracer)
!------------------------------------------------------------------------------
    class (Tracer), pointer :: pTracer

    call  Clayilli_setSpec('Clayilli') ! http://webmineral.com/data/Illite.shtml
    call  Claykaol_setSpec('Claykaol') ! http://www.webmineral.com/data/Kaolinite.shtml
    call  Claysmec_setSpec('Claysmec') ! http://www.webmineral.com/data/Rectorite.shtml
    call  Claycalc_setSpec('Claycalc') ! http://www.webmineral.com/data/Calcite.shtml
    call  Clayquar_setSpec('Clayquar') ! http://www.webmineral.com/data/Quartz.shtml
    call  Sil1quar_setSpec('Sil1quar') ! http://www.webmineral.com/data/Quartz.shtml
    call  Sil1feld_setSpec('Sil1feld') ! http://www.mindat.org/min-1624.html
    call  Sil1calc_setSpec('Sil1calc') ! http://www.webmineral.com/data/Calcite.shtml
    call  Sil1hema_setSpec('Sil1hema') ! http://www.webmineral.com/data/Hematite.shtml
    call  Sil1gyps_setSpec('Sil1gyps') ! http://www.webmineral.com/data/Gypsum.shtml
    call  Sil2quar_setSpec('Sil2quar') ! http://www.webmineral.com/data/Quartz.shtml
    call  Sil2feld_setSpec('Sil2feld') ! http://www.mindat.org/min-1624.html
    call  Sil2calc_setSpec('Sil2calc') ! http://www.webmineral.com/data/Calcite.shtml
    call  Sil2hema_setSpec('Sil2hema') ! http://www.webmineral.com/data/Hematite.shtml
    call  Sil2gyps_setSpec('Sil2gyps') ! http://www.webmineral.com/data/Gypsum.shtml
    call  Sil3quar_setSpec('Sil3quar') ! http://www.webmineral.com/data/Quartz.shtml
    call  Sil3feld_setSpec('Sil3feld') ! http://www.mindat.org/min-1624.html
    call  Sil3calc_setSpec('Sil3calc') ! http://www.webmineral.com/data/Calcite.shtml
    call  Sil3hema_setSpec('Sil3hema') ! http://www.webmineral.com/data/Hematite.shtml
    call  Sil3gyps_setSpec('Sil3gyps') ! http://www.webmineral.com/data/Gypsum.shtml

!------------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------------

    subroutine Clayilli_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Clayilli=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.795d3) ! measured; http://www.mindat.org/min-2011.html
      if (tracers_drydep) call set_trradius(n, 0.46D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)

    end subroutine Clayilli_setSpec

    subroutine Claykaol_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Claykaol=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.63d3) ! calculated; http://www.mindat.org/min-2011.html
      if (tracers_drydep) call set_trradius(n, 0.46D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)

    end subroutine Claykaol_setSpec

    subroutine Claysmec_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Claysmec=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.35D3) ! for Montmorillonite
      if (tracers_drydep) call set_trradius(n, 0.46D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)

    end subroutine Claysmec_setSpec

    subroutine Claycalc_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Claycalc=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.71d3) ! measured; http://www.mindat.org/min-859.html
      if (tracers_drydep) call set_trradius(n, 0.46D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Claycalc_setSpec

    subroutine Clayquar_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Clayquar=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, DensityQuartz)
      if (tracers_drydep) call set_trradius(n, 0.46D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Clayquar_setSpec

    subroutine Sil1quar_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil1quar=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, DensityQuartz)
      if (tracers_drydep) call set_trradius(n, 1.47D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil1quar_setSpec

    subroutine Sil1feld_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil1feld=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.68d3) ! average Plagioclase
      if (tracers_drydep) call set_trradius(n, 1.47D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil1feld_setSpec

    subroutine Sil1calc_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil1calc=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.71d3) ! measured; http://www.mindat.org/min-859.html
      if (tracers_drydep) call set_trradius(n, 1.47D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil1calc_setSpec

    subroutine Sil1hema_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil1hema=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, DensityHematite)
      if (tracers_drydep) call set_trradius(n, 1.47D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil1hema_setSpec

    subroutine Sil1gyps_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil1gyps=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.312d3) ! measured; http://www.mindat.org/min-1784.html
      if (tracers_drydep) call set_trradius(n, 1.47D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil1gyps_setSpec

    subroutine Sil2quar_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil2quar=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, DensityQuartz)
      if (tracers_drydep) call set_trradius(n, 2.94D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil2quar_setSpec

    subroutine Sil2feld_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil2feld=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.68d3) ! average Plagioclase
      if (tracers_drydep) call set_trradius(n, 2.94D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil2feld_setSpec

    subroutine Sil2calc_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil2calc=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.71d3) ! measured; http://www.mindat.org/min-859.html
      if (tracers_drydep) call set_trradius(n, 2.94D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil2calc_setSpec

    subroutine Sil2hema_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil2hema=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, DensityHematite)
      if (tracers_drydep) call set_trradius(n, 2.94D-06)

      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil2hema_setSpec

    subroutine Sil2gyps_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil2gyps=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.312d3) ! measured; http://www.mindat.org/min-1784.html
      if (tracers_drydep) call set_trradius(n, 2.94D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil2gyps_setSpec

    subroutine Sil3quar_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil3quar=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, DensityQuartz)
      if (tracers_drydep) call set_trradius(n, 5.88D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil3quar_setSpec

    subroutine Sil3feld_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil3feld=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.68d3) ! average Plagioclase
      if (tracers_drydep) call set_trradius(n, 5.88D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil3feld_setSpec

    subroutine Sil3calc_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil3calc=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.71d3) ! measured; http://www.mindat.org/min-859.html
      if (tracers_drydep) call set_trradius(n, 5.88D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil3calc_setSpec

    subroutine Sil3hema_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil3hema=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, DensityHematite)
      if (tracers_drydep) call set_trradius(n, 5.88D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil3hema_setSpec

    subroutine Sil3gyps_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Sil3gyps=n
      call set_ntm_power(n, -9)
      call set_trpdens(n, 2.312d3) ! measured; http://www.mindat.org/min-1784.html
      if (tracers_drydep) call set_trradius(n, 5.88D-06)
      call set_fq_aer(n, 5.D-1)
      call set_rc_washt(n, 5.D-1)
      call set_tr_wd_type(n, nPART)
      call set_tr_mm(n, 1.d+0)
      call set_isdust(n, 1)
    end subroutine Sil3gyps_setSpec

  end subroutine Minerals_InitMetadata

end module MineralsTracersMetadata_mod



