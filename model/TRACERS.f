#include "rundeck_opts.h"

#ifdef TRACERS_ON
!@sum  TRACERS: generic tracer routines used for all tracers
!@+    Routines included:
!@+      Generic diags: set_generic_tracer_diags
!@+      Apply previously set sources: apply_tracer_sources
!@+      Radioactive Decay: tdecay
!@+      Gravitaional Settling: trgrav
!@+      Check routine: checktr
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0

      SUBROUTINE set_generic_tracer_diags
!@sum set_generic_tracer_diags init trace gas attributes and diagnostics
!@auth J. Lerner
!@calls sync_param
      USE CONSTANT, only: mair,sday
      USE MODEL_COM, only: dtsrc,nisurf
      USE DIAG_COM, only: ia_src,ia_12hr,ir_log2,ir_0_71
      USE TRACER_COM
      USE TRDIAG_COM
      USE PARAM
      USE DOMAIN_DECOMP, only: AM_I_ROOT
      implicit none
      integer :: l,k,n
      character*10, DIMENSION(NTM) :: CMR,CMRWT
      logical :: T=.TRUE. , F=.FALSE.
      character*50 :: unit_string

C**** Get factor to convert from mass to volume mr if necessary
      do n=1,ntm
        if (to_volume_MixRat(n) .eq.1) then
          MMR_to_VMR(n) = mass2vol(n)
          cmr(n) = 'V/V air'
        else
          MMR_to_VMR(n) = 1.d0
          cmr(n) = 'kg/kg air'
        endif
#ifdef TRACERS_WATER
        if (to_per_mil(n) .eq.1) then
          cmrwt(n) = 'per mil'
          cmr(n) = 'per mil'     ! this overrides to_volume_ratio
          MMR_to_VMR(n) = 1.
        else
          cmrwt(n) = 'kg/m^2'
        end if
#endif
      end do
C****
C**** TAJLN(J,L,KQ,N)  (SUM OVER LONGITUDE AND TIME OF)
C****
C**** jlq_power Exponent associated with a physical process
C****      (for printing tracers).   (ntm_power+jlq_power)
C**** jls_power Exponent associated with a source/sink (for printing)
C**** Defaults for JLN
      scale_jlq(:) = 1./DTsrc
      jgrid_jlq(:) = 1
      ia_jlq(:) = ia_src

C**** Tracer concentration
      do n=1,ntm
        k = 1        ! <<<<< Be sure to do this
        jlnt_conc = k
        sname_jln(k,n) = trim(trname(n))//'_CONCENTRATION'
        lname_jln(k,n) = trim(trname(n))//' CONCENTRATION'
        jlq_power(k) = 0.
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),cmr(n))
        scale_jlq(k) = 1.d0
        scale_jln(n) = MMR_to_VMR(n)
#ifdef TRACERS_WATER
        if (to_per_mil(n).gt.0) units_jln(k,n) = unit_string(0,cmr(n))
#endif
C**** Tracer mass
        k = k + 1
        jlnt_mass = k
        sname_jln(k,n) = trim(trname(n))//'_MASS'
        lname_jln(k,n) = trim(trname(n))//' MASS'
        jlq_power(k) = 4
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_OM_SP) ||\
    (defined TRACERS_AMP)
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k)+13
     *       ,'kg')
#else
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k)
     *       ,'kg/m^2')
#endif
        scale_jlq(k) = 1.d0
#ifdef TRACERS_WATER
C****   TRACER CONCENTRATION IN CLOUD WATER
        k = k + 1
        jlnt_cldh2o = k
        sname_jln(k,n) = trim(trname(n))//'_WM_CONC'
        lname_jln(k,n) = trim(trname(n))//' CLOUD WATER CONCENTRATION'
        jlq_power(k) = 4
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k)
     *       ,'kg/kg water')
        scale_jlq(k) = 1.d0
        if (to_per_mil(n).gt.0) units_jln(k,n) = unit_string(0,cmr(n))
#endif
C**** Physical processes affecting tracers
C****   F (TOTAL NORTHWARD TRANSPORT OF TRACER MASS)  (kg)
        k = k + 1
        jlnt_nt_tot = k
        sname_jln(k,n) = 'tr_nt_tot_'//trname(n)
        lname_jln(k,n) = 'TOTAL NORTHWARD TRANSPORT OF '//
     &     trim(trname(n))//' MASS'
        jlq_power(k) = 11
        jgrid_jlq(k) = 2
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),'kg/s')
C****   STM/SM (MEAN MERIDIONAL N.T. OF TRACER MASS)  (kg)
        k = k + 1
        jlnt_nt_mm = k
        sname_jln(k,n) = 'tr_nt_mm_'//trname(n)
        lname_jln(k,n) = 'NORTHWARD TRANS. OF '//
     &     trim(trname(n))//' MASS BY MERIDIONAL CIRC.'
        jlq_power(k) = 10
        jgrid_jlq(k) = 2
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),'kg/s')
C****   F (TOTAL VERTICAL TRANSPORT OF TRACER MASS)  (kg)
        k = k + 1
        jlnt_vt_tot = k
        sname_jln(k,n) = 'tr_vt_tot_'//trname(n)
        lname_jln(k,n) = 'TOTAL VERTICAL TRANSPORT OF '//
     &     trim(trname(n))//' MASS'
        jlq_power(k) = 11
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),'kg/s')
C****   STM/SM (MEAN MERIDIONAL V.T. OF TRACER MASS)  (kg)
        k = k + 1
        jlnt_vt_mm = k
        sname_jln(k,n) = 'tr_vt_mm_'//trname(n)
        lname_jln(k,n) = 'VERTICAL TRANS. OF '//
     &     trim(trname(n))//' MASS BY MERIDIONAL CIRC.'
        jlq_power(k) = 10.
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),'kg/s')
C****   TMBAR-TM (CHANGE OF TRACER MASS BY MOIST CONVEC)(kg)
        k = k + 1
        jlnt_mc = k
        sname_jln(k,n) = 'tr_mc_'//trname(n)
        lname_jln(k,n) = 'CHANGE OF '//
     &     trim(trname(n))//' MASS BY MOIST CONVECTION'
        jlq_power(k) = 10
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),'kg/s')
C****   TMBAR-TM (CHANGE OF TRACER MASS BY Large-scale CONDENSE)  (kg)
        k = k + 1
        jlnt_lscond = k
        sname_jln(k,n) = 'tr_lscond'//trname(n)
        lname_jln(k,n) ='CHANGE OF '//
     &     trim(trname(n))//' MASS BY LARGE-SCALE CONDENSE'
        jlq_power(k) = 10.
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),'kg/s')
C****   TMBAR-TM (CHANGE OF TRACER MASS BY DRY CONVEC)  (kg)
        k = k + 1
        jlnt_turb = k
        sname_jln(k,n) = 'tr_turb_'//trname(n)
        lname_jln(k,n) = 'CHANGE OF '//
     &     trim(trname(n))//' MASS BY TURBULENCE/DRY CONVECTION'
        jlq_power(k) = 10
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),'kg/s')

      end do

      if (k.gt. ktajl) then
        if (AM_I_ROOT()) write (6,*)
     &   'tjl_defs: Increase ktajl=',ktajl,' to at least ',k
        call stop_model('ktajl too small',255)
      end if

C**** CONTENTS OF TAIJLN(I,J,LM,N)  (SUM OVER TIME OF)
C****        TML (M*M * KG TRACER/KG AIR)
C**** Set defaults that are true for all tracers and layers
      ia_ijt    = ia_src
      ir_ijt(:) = ir_log2   !n
      do n=1,ntm
        ijtm_power(n) = ntm_power(n)+4 !n for integrated mass
        ijtc_power(n) = ntm_power(n)+1 !n for concentration
#ifdef TRACERS_WATER
        if (to_per_mil(n) .eq.1) ijtc_power(n) = 0
#endif
      end do
C**** Tracer concentrations (TAIJLN)
      do n=1,ntm
      do l=1,lm
        write(sname_ijt(l,n),'(a,i2.2)') trim(TRNAME(n))//'_L_',l
        write(lname_ijt(l,n),'(a,i2)')   trim(TRNAME(n))//' L ',l
        units_ijt(l,n) = unit_string(ijtc_power(n),cmr(n))
        scale_ijt(l,n) = MMR_to_VMR(n)*10.**(-ijtc_power(n))
      end do
      end do

C**** AIJN
C****     1  TM (SUM OVER ALL LAYERS) (M*M * KG TRACER/KG AIR)
C****     2  TRS (SURFACE TRACER CONC.) (M*M * KG TRACER/KG AIR)
C****     3  TM (SUM OVER ALL LAYERS) (M*M * KG TRACER)
      do n=1,ntm
C**** Summation of mass over all layers
      k = 1        ! <<<<< Be sure to do this
      tij_mass = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_Total_Mass'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//' Total Mass'
        units_tij(k,n) = unit_string(ijtm_power(n),'kg/m^2')
        scale_tij(k,n) = 10.**(-ijtm_power(n))
C**** Average concentration over layers
      k = k+1
      tij_conc = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_Average'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//' Average'
        units_tij(k,n) = unit_string(ijtc_power(n),cmr(n))
        scale_tij(k,n) = MMR_to_VMR(n)*10.**(-ijtc_power(n))
C**** Surface concentration
      k = k+1
      tij_surf = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_At_Surface'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//' At Surface'
        units_tij(k,n) = unit_string(ijtc_power(n),cmr(n))
        scale_tij(k,n)=MMR_to_VMR(n)*10.**(-ijtc_power(n))/
     *                 REAL(NIsurf,KIND=8)
C**** Surface concentration by volume (units kg/m^3)
      k = k+1
      tij_surfbv = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *     '_byVol_At_Surface'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' By Vol At Surface'
        units_tij(k,n) = unit_string(ijtc_power(n),'kg/m^3')
        scale_tij(k,n)=MMR_to_VMR(n)*10.**(-ijtc_power(n))/
     *                 REAL(NIsurf,KIND=8)
#ifdef TRACERS_GASEXCH_Natassa
C**** Gas Exchange Solubility coefficient
      k = k+1
      tij_alpha = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_Solubility'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//' Solubility'
        ijtc_power(n) = 0
        units_tij(k,n) = unit_string(ijtc_power(n),'mol/l/atm')
        scale_tij(k,n) = 10.**(-ijtc_power(n))
C**** Gas Exchange Coefficient (piston velocity)
      k = k+1
      tij_kw = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_Piston_Veloc'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//' Piston Veloc'
        ijtc_power(n) = -4
        units_tij(k,n) = unit_string(ijtc_power(n),'m/s')
        scale_tij(k,n) = 10.**(-ijtc_power(n))
C**** Gas Exchange flux at ocean surface
      k = k+1
      tij_gasx = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_Gas_Exch_Flux'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//' Gas Exch Flux'
        ijtc_power(n) = -6
        units_tij(k,n) = unit_string(ijtc_power(n),'kg/m^2')
        scale_tij(k,n) = 10.**(-ijtc_power(n))

      print*,'TRACERS.f Natassa: ktaij=',ktaij,', k=',k
#endif
#ifdef TRACERS_WATER
C**** the following diagnostics are set assuming that the particular
C**** tracer exists in water.
C**** Tracers in precipitation (=Wet deposition)
      k = k+1
      tij_prec = k
      if (dowetdep(n)) then
        if (to_per_mil(n) .eq.1) then
          write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_in_prec'
          write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *         ' in Precip'
          units_tij(k,n)=cmrwt(n)
          scale_tij(k,n)=10.**(-ijtc_power(n))/dtsrc
        else
          write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_wet_dep'
          write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *         ' Wet Deposition'
          units_tij(k,n)=unit_string(ijtc_power(n)-5,trim(cmrwt(n))
     *         //'/s')
          scale_tij(k,n)=10.**(-ijtc_power(n)+5)/dtsrc
        end if
      end if
C**** Tracers in evaporation
      k = k+1
      tij_evap = k
      if (tr_wd_type(n).eq.nWater) then
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_in_evap'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' in Evaporation'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(ijtc_power(n),cmrwt(n))
        else
          units_tij(k,n)=unit_string(ijtc_power(n),trim(cmrwt(n))//'/s')
        end if
        scale_tij(k,n)=10.**(-ijtc_power(n))/dtsrc
      endif
C**** Tracers in river runoff
      k = k+1
      tij_rvr = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_in_rvr'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' in River Outflow'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          scale_tij(k,n)=1.
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg')
          scale_tij(k,n)=10.**(-ijtc_power(n)-3)
        end if
C**** Tracers in sea ice
      k = k+1
      tij_seaice = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_in_ice'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' in Sea Ice  x POICE'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          scale_tij(k,n)=1.
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,cmrwt(n))
          scale_tij(k,n)=10.**(-ijtc_power(n)-3)
        end if
C**** Tracers conc. in ground component (ie. water or ice surfaces)
      k = k+1
      tij_grnd = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_at_Grnd'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' at Ground'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg wat')
        end if
        scale_tij(k,n)=10.**(-ijtc_power(n)-3)/REAL(NIsurf,KIND=8)
C**** Tracers conc. in lakes (layer 1)
      k = k+1
      tij_lk1 = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_Lake1'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' Lakes layer 1'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg wat')
        end if
        scale_tij(k,n)=10.**(-ijtc_power(n)-3)/REAL(NIsurf,KIND=8)
C**** Tracers conc. in lakes (layer 2)
      k = k+1
      tij_lk2 = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_Lake2'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' Lakes layer 2'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg wat')
        end if
        scale_tij(k,n)=10.**(-ijtc_power(n)-3)/REAL(NIsurf,KIND=8)
C**** Tracers conc. in soil water
      k = k+1
      tij_soil = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_in_Soil'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' Soil Water'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg wat')
        end if
        scale_tij(k,n)=10.**(-ijtc_power(n)-3)/REAL(NIsurf,KIND=8)
C**** Tracers conc. in land snow water
      k = k+1
      tij_snow = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_in_Snow'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' Land Snow Water'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg wat')
        end if
        scale_tij(k,n)=10.**(-ijtc_power(n)-3)/REAL(NIsurf,KIND=8)
C**** Tracer ice-ocean flux
      k = k+1
      tij_icocflx = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_ic_oc_flx'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' Ice-Ocean Flux'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
        else
          units_tij(k,n)=unit_string(ijtc_power(n)-5,'kg/m^2/s')
        end if
        scale_tij(k,n)=10.**(-ijtc_power(n)+5)/DTsrc
C**** Tracers integrated E-W atmospheric flux
      k = k+1
      tij_uflx = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_uflx'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' E-W Atmos Flux'
        units_tij(k,n)=unit_string(ijtc_power(n)+10,'kg/s')
        scale_tij(k,n)=10.**(-ijtc_power(n)-10)/DTsrc
C**** Tracers integrated E-W atmospheric flux
      k = k+1
      tij_vflx = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_vflx'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' N-S Atmos Flux'
        units_tij(k,n)=unit_string(ijtc_power(n)+10,'kg/s')
        scale_tij(k,n)=10.**(-ijtc_power(n)-10)/DTsrc
#endif
#ifdef TRACERS_DRYDEP
C**** Tracers dry deposition flux.
      k = k+1
      tij_drydep = k
      if (dodrydep(n)) then
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_dry_dep'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' Dry Deposition'
        units_tij(k,n)=unit_string(ijtc_power(n)-5,'kg/m^2/s')
        scale_tij(k,n)=10.**(-ijtc_power(n)+5)/DTsrc
      end if
      k = k+1
      tij_gsdep = k
      if (trradius(n).gt.0) then
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_gs_dep'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' Gravitional Settling'
        units_tij(k,n)=unit_string(ijtc_power(n)-5,'kg/m^2/s')
        scale_tij(k,n)=10.**(-ijtc_power(n)+5)/DTsrc
      end if
#endif

      if (k .gt. ktaij) then
        if (AM_I_ROOT()) write (6,*)
     &   'tij_defs: Increase ktaij=',ktaij,' to at least ',k
        call stop_model('ktaij too small',255)
      end if

      end do

      END SUBROUTINE set_generic_tracer_diags


      SUBROUTINE apply_tracer_2Dsource(dtstep)
!@sum apply_tracer_2Dsource adds surface sources to tracers
!@auth Jean Lerner/Gavin Schmidt
      USE MODEL_COM, only : jm
      USE GEOM, only : imaxj
      USE QUSDEF, only : mz,mzz
      USE TRACER_COM, only : ntm,trm,trmom,ntsurfsrc,ntisurfsrc,trname
      USE FLUXES, only : trsource,trflux1,trsrfflx
      USE TRDIAG_COM, only : taijs=>taijs_loc
      USE TRDIAG_COM, only : tajls=>tajls_loc,ijts_source,jls_source
     *     ,itcon_surf
      USE DOMAIN_DECOMP, ONLY : GRID, GET
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: dtstep
      INTEGER n,ns,naij,najl,j,i
      REAL*8, DIMENSION(grid%J_STRT_HALO:grid%J_STOP_HALO) :: dtracer
      REAL*8 ftr1

      INTEGER :: J_0, J_1

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)
C**** This is tracer independent coding designed to work for all
C**** surface sources.
C**** Note that tracer flux is added to first layer either implicitly
C**** in ATURB or explicitly in 'apply_fluxes_to_atm' call in SURFACE.

      do n=1,ntm
        trflux1(:,:,n) = 0.
C**** Non-interactive sources
        do ns=1,ntsurfsrc(n)
C**** diagnostics
          naij = ijts_source(ns,n)
          IF (naij > 0) THEN
          taijs(:,:,naij) = taijs(:,:,naij) + trsource(:,:,ns,n)*dtstep
          ENDIF
          najl = jls_source(ns,n)
          IF (najl > 0) THEN
            DO J=J_0,J_1
              tajls(j,1,najl) = tajls(j,1,najl)+
     *             sum(trsource(1:imaxj(j),j,ns,n))*dtstep
            END  DO
          END IF
          dtracer(:)=0.
          DO J=J_0,J_1
            dtracer(j)=0.
            do i=1,imaxj(j)
              dtracer(j)=dtracer(j)+trsource(i,j,ns,n)*dtstep
            end do
          end do
          if (itcon_surf(ns,n).gt.0)
     *         call DIAGTCB(dtracer(:),itcon_surf(ns,n),n)
C**** trflux1 is total flux into first layer
          trflux1(:,:,n) = trflux1(:,:,n)+trsource(:,:,ns,n)
        end do

C**** modify vertical moments (only from non-interactive sources)
       do j=j_0,j_1
        trmom( mz,:,j,1,n) = trmom( mz,:,j,1,n)-1.5*trflux1(:,j,n)
     *       *dtstep
        trmom(mzz,:,j,1,n) = trmom(mzz,:,j,1,n)+0.5*trflux1(:,j,n)
     *       *dtstep

C**** Accumulate interactive sources as well
        trflux1(:,j,n) = trflux1(:,j,n)+trsrfflx(:,j,n)
       end do

C**** Technically speaking the vertical moments should be modified here
C**** as well. But for consistency with water vapour we only modify
C**** moments for dew.
        do j=J_0,J_1
          do i=1,imaxj(j)
            if (trsrfflx(i,j,n).lt.0 .and. trm(i,j,1,n).gt.0) then
              ftr1=-trsrfflx(i,j,n)*dtstep/trm(i,j,1,n)
              trmom(:,i,j,1,n)=trmom(:,i,j,1,n)*(1.-ftr1)
            end if
          end do
        end do
      end do
C****
      RETURN
      END SUBROUTINE apply_tracer_2Dsource

      MODULE apply3d
!@sum apply3d is used simply so that I can get optional arguments
!@+   to work. If anyone can some up with something neater, let me know.
      CONTAINS
      SUBROUTINE apply_tracer_3Dsource( ns , n , momlog )
!@sum apply_tracer_3Dsource adds 3D sources to tracers
!@auth Jean Lerner/Gavin Schmidt
      USE CONSTANT, only : teeny
      USE MODEL_COM, only : jm,im,lm,dtsrc
      USE GEOM, only : imaxj,bydxyp,lat_dg,lon_dg
      USE QUSDEF, only: nmom
      USE TRACER_COM, only : ntm,trm,trmom,trname,alter_sources,
     * num_regions,reg_S,reg_N,reg_E,reg_W,ef_FACT3d,num_tr_sectors
     * ,tr_sect_index3D,num_tr_sectors3d
      USE FLUXES, only : tr3Dsource
      USE TRDIAG_COM, only : tajls=>tajls_loc,jls_3Dsource,itcon_3Dsrc
     *     ,ijts_3Dsource,taijs=>taijs_loc
      USE DOMAIN_DECOMP, only : GRID, GET, write_parallel
      IMPLICIT NONE
!@var MOM true (default) if moments are to be modified
      logical, optional, intent(in) :: momlog
      integer, intent(in) :: n,ns
      real*8 fr3d,taijsum(im,grid%j_strt_halo:grid%j_stop_halo,lm)
      logical :: domom
      integer najl,i,j,l,naij,kreg,nsect,nn
      INTEGER :: J_0, J_1

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

C**** Ensure that this is a valid tracer and source
      if (n.eq.0 .or. ns.eq.0) return
C**** parse options
      domom=.true.
      if( present(momlog) ) then
        domom=momlog
      end if
C**** This is tracer independent coding designed to work for all
C**** 3D sources.
C**** Modify tracer amount, moments, and diagnostics
      najl = jls_3Dsource(ns,n)
      naij = ijts_3Dsource(ns,n)

C**** apply tracer source alterations if requested in rundeck:
      if(alter_sources)then
      do nsect=1,num_tr_sectors3d(n,ns)
      do l=1,lm ! making non-openMP suspecting bug...
      do j=j_0,j_1
      do i=1,imaxj(j)
        do kreg=1,num_regions
          if(lat_dg(j,1) >= REG_S(kreg) .and. lat_dg(j,1)
     &    <= REG_N(kreg) .and. lon_dg(i,1) >= REG_W(kreg)
     &    .and. lon_dg(i,1) < REG_E(kreg) )then
            nn=tr_sect_index3D(n,ns,nsect)
            if(ef_FACT3d(nn,kreg) > -1.e20)
     &      tr3Dsource(i,j,l,ns,n) = tr3Dsource(i,j,l,ns,n) *
     &      ef_FACT3d(nn,kreg) 
          endif 
        enddo
      enddo
      enddo 
      enddo
      enddo
      endif

!$OMP PARALLEL DO PRIVATE (L,I,J,fr3d)
      do l=1,lm
      do j=j_0,j_1
      do i=1,imaxj(j)
C**** calculate fractional loss
        fr3d=0.
        if (tr3Dsource(i,j,l,ns,n).lt.0.) then
          fr3d = -tr3Dsource(i,j,l,ns,n)*dtsrc/(trm(i,j,l,n)+
     *         tiny(trm(i,j,l,n)))
          fr3d=min(1d0,fr3d)
          if (domom) trmom(1:nmom,i,j,l,n) =
     *         trmom(1:nmom,i,j,l,n)*(1.-fr3d)
        end if
C**** update tracer mass and diagnostics
        trm(i,j,l,n) = trm(i,j,l,n)+tr3Dsource(i,j,l,ns,n)*dtsrc
        if (1.-fr3d.le.1d-16) trm(i,j,l,n) = 0.
        IF (jls_3Dsource(ns,n) > 0) tajls(j,l,najl)
     &       =tajls(j,l,najl)+tr3Dsource(i,j,l,ns,n)*dtsrc
        taijsum(i,j,l)=tr3Dsource(i,j,l,ns,n)*dtsrc
      end do; end do; end do
!$OMP END PARALLEL DO
      if (naij.gt.0) then
      do j=j_0,j_1
        do i=1,imaxj(j)
          taijs(i,j,naij) = taijs(i,j,naij) + sum(taijsum(i,j,:))
        end do
      end do
      end if
      if (itcon_3Dsrc(ns,n).gt.0)
     *  call DIAGTCA(itcon_3Dsrc(ns,n),n)

C****
      RETURN
      END SUBROUTINE apply_tracer_3Dsource

      end module apply3d
      SUBROUTINE TDECAY
!@sum TDECAY decays radioactive tracers every source time step
!@auth Gavin Schmidt/Jean Lerner
      USE MODEL_COM, only : im,jm,lm,itime,dtsrc
      USE FLUXES, only : tr3Dsource
      USE GEOM, only : imaxj
      USE TRACER_COM, only : ntm,trm,trmom,trdecay,itime_tr0,n_Pb210,
     &     trname, n_Rn222
#ifdef TRACERS_WATER
     *     ,trwm
      USE SEAICE_COM, only : trsi
      USE LAKES_COM, only : trlake
      USE LANDICE_COM, only : trlndi,trsnowli
      USE GHY_COM, only : tr_w_ij,tr_wsn_ij
#endif
      USE TRDIAG_COM, only : tajls=>tajls_loc,jls_decay,itcon_decay
      USE DOMAIN_DECOMP, only : GRID, GET
      IMPLICIT NONE
      real*8, save, dimension(ntm) :: expdec = 1.
      real*8, dimension(im,grid%J_STRT_HALO:grid%J_STOP_HALO,lm) :: told
      logical, save :: ifirst=.true.
      integer n,najl,j,l
      integer :: J_0, J_1

      CALL GET(grid, J_STRT = J_0, J_STOP = J_1)

      if (ifirst) then
        do n=1,ntm
          if (trdecay(n).gt.0.0) expdec(n)=exp(-trdecay(n)*dtsrc)
        end do
        ifirst = .false.
      end if

      do n=1,ntm
        if (trdecay(n).gt.0. .and. itime.ge.itime_tr0(n)) then
C**** Atmospheric decay
          told(:,:,:)=trm(:,:,:,n)

#ifdef TRACERS_WATER
     *               +trwm(:,:,:,n)
          trwm(:,:,:,n)   = expdec(n)*trwm(:,:,:,n)
#endif
          if (trname(n) .eq. "Rn222" .and. n_Pb210.gt.0) then
            tr3Dsource(:,:,:,1,n_Pb210)= trm(:,:,:,n)*(1-expdec(n))*210.
     *           /222./dtsrc
          end if

          trm(:,:,:,n)    = expdec(n)*trm(:,:,:,n)
          trmom(:,:,:,:,n)= expdec(n)*trmom(:,:,:,:,n)

#ifdef TRACERS_WATER
C**** Note that ocean tracers are dealt with by separate ocean code.
C**** Decay sea ice tracers
          trsi(n,:,:,:)   = expdec(n)*trsi(n,:,:,:)
C**** ...lake tracers
          trlake(n,:,:,:) = expdec(n)*trlake(n,:,:,:)
C**** ...land surface tracers
          tr_w_ij(n,:,:,:,:) = expdec(n)*tr_w_ij(n,:,:,:,:)
          tr_wsn_ij(n,:,:,:,:)= expdec(n)*tr_wsn_ij(n,:,:,:,:)
          trsnowli(n,:,:) = expdec(n)*trsnowli(n,:,:)
          trlndi(n,:,:)   = expdec(n)*trlndi(n,:,:)
#endif
C**** atmospheric diagnostics
          najl = jls_decay(n)
          do l=1,lm
          do j=J_0,J_1
            tajls(j,l,najl)=tajls(j,l,najl)+sum(trm(1:imaxj(j),j,l,n)
#ifdef TRACERS_WATER
     *           +trwm(1:imaxj(j),j,l,n)
#endif
     *           -told(1:imaxj(j),j,l))
          enddo
          enddo


          call DIAGTCA(itcon_decay(n),n)
        end if
      end do
C****
      return
      end subroutine tdecay


      SUBROUTINE TRGRAV
!@sum TRGRAV gravitationally settles particular tracers
!@auth Gavin Schmidt/Reha Cakmur
      USE CONSTANT, only : grav,deltx,lhe,rgas
      USE MODEL_COM, only : im,jm,lm,itime,dtsrc,zatmo,t,q
      USE GEOM, only : imaxj,bydxyp
      USE SOMTQ_COM, only : mz,mzz,mzx,myz,zmoms
      USE DYNAMICS, only : gz,pmid,pk
      USE TRACER_COM, only : ntm,trm,trmom,itime_tr0,trradius
     *     ,trname,trpdens
#ifdef TRACERS_AMP
     *     ,AMP_MODES_MAP,ntmAMP
      USE AMP_AEROSOL, only : DIAM, AMP_dens
#endif
      USE TRDIAG_COM, only : tajls=>tajls_loc,jls_grav
      USE DOMAIN_DECOMP, only : GRID, GET
      IMPLICIT NONE
      real*8 :: stokevdt,fgrfluxd,fluxd,fluxu,press,airden,temp,rh,qsat
     *     ,vgs,tr_radius,tr_dens
      real*8, dimension(im,grid%J_STRT_HALO:grid%J_STOP_HALO,lm) :: told
      integer n,najl,i,j,l
      integer :: J_0, J_1
      logical :: hydrate

      CALL GET(grid, J_STRT = J_0, J_STOP = J_1)

      do n=1,ntm
        if (trradius(n).gt.0. .and. itime.ge.itime_tr0(n)) then
C**** need to hydrate the sea salt before determining settling
          hydrate = (trname(n).eq.'seasalt1'.or.trname(n).eq.'seasalt2'
     *         .or.trname(n).eq.'M_SSA_SS'.or. trname(n).eq.'M_SSC_SS'
     *         .or.trname(n).eq.'M_SSS_SS')

C**** Gravitational settling
!$OMP  PARALLEL DO PRIVATE (l,i,j,press,airden,temp,rh,stokevdt,
!$OMP* fgrfluxd,fluxd,fluxu)

          do j=J_0,J_1
          do i=1,imaxj(j)
            fluxd=0.
            do l=lm,1,-1        ! loop down

               tr_dens = trpdens(n)
               tr_radius = trradius(n)

#ifdef TRACERS_AMP
      if (n.le.ntmAMP) then
        if(AMP_MODES_MAP(n).gt.0) then
        if(DIAM(i,j,l,AMP_MODES_MAP(n)).gt.0.) 
     &      tr_radius =DIAM(i,j,l,AMP_MODES_MAP(n)) *0.5

          call AMPtrdens(i,j,l,n)

          tr_dens =AMP_dens(i,j,l,AMP_MODES_MAP(n))

        endif   
      endif 
#endif  
              told(i,j,l)=trm(i,j,l,n)
C**** air density + relative humidity (wrt water)
              press=pmid(l,i,j)
              temp=pk(l,i,j)*t(i,j,l)
              airden=100.d0*press/(rgas*temp*(1.+q(i,j,l)*deltx))
              rh=q(i,j,l)/qsat(temp,lhe,press)

C**** calculate stoke's velocity (including possible hydration effects
C**** and slip correction factor)
              stokevdt=dtsrc*vgs(airden,rh,tr_radius,tr_dens,hydrate)

              fluxu=fluxd       ! from previous level

C**** Calculate height differences using geopotential
              if (l.eq.1) then  ! layer 1
                fgrfluxd=0.     ! calc now done in PBL
              else              ! above layer 1
                fgrfluxd=stokevdt*grav/(gz(i,j,l)-gz(i,j,l-1))
              end if
              fluxd = trm(i,j,l,n)*fgrfluxd ! total flux down
              trm(i,j,l,n) = trm(i,j,l,n)*(1.-fgrfluxd)+fluxu
              if (1.-fgrfluxd.le.1d-16) trm(i,j,l,n) = fluxu
              trmom(zmoms,i,j,l,n) = trmom(zmoms,i,j,l,n)*(1.-fgrfluxd)
              
            end do
          end do
          end do
!$OMP END PARALLEL DO
          najl = jls_grav(n)
          IF (najl > 0) THEN
            do l=1,lm
              do j=J_0,J_1
                tajls(j,l,najl)=tajls(j,l,najl)
     &               +sum(trm(1:imaxj(j),j,l,n)-told(1:imaxj(j),j,l))
              enddo
            enddo
          END IF
        end if
      end do
C****
      return
      end subroutine trgrav

      REAL*8 FUNCTION vgs(airden,rh1,tr_radius,tr_dens,hydrate)
!@sum vgs returns settling velocity for tracers (m/s)
!@auth Gavin Schmidt/Reha Cakmur
      USE CONSTANT, only : visc_air,by3,pi,gasc,avog,rt2,deltx
     *     ,mair,grav
      IMPLICIT NONE
      real*8, intent(in) ::  airden,rh1,tr_radius,tr_dens
      logical, intent(in) :: hydrate
      real*8  wmf,frpath
      real*8, parameter :: dair=3.65d-10 !m diameter of air molecule
C**** coefficients for slip factor calculation
      real*8, parameter :: s1=1.247d0, s2=0.4d0, s3=1.1d0
C**** coefficients for hydration radius if required 
      real*8, parameter :: c1=0.7674d0, c2=3.079d0, c3=2.573d-11,
     *     c4=-1.424d0
      real*8 dens, rad, rh

C**** Determine if hydration of aerosols is required
      if (hydrate) then
        rh=max(0.01d0,min(rh1,0.99d0))
C**** hydrated radius (reference?)
        rad=(c1*tr_radius**(c2)/(c3*tr_radius**(c4)-log10(rh))
     *       + tr_radius**3)**by3
C**** hydrated density
        dens=((rad**3 - tr_radius**3)*1000.d0
     *       + tr_radius**3*tr_dens)/rad**3
      else                      ! take dry values
        dens = tr_dens
        rad = tr_radius 
      end if

C**** calculate stoke's velocity
      vgs=2.*grav*dens*rad**2/(9.*visc_air)

C**** slip correction factor
c wmf is the additional velocity if the particle size is small compared
c   to the mean free path of the air; important in the stratosphere
      frpath=1d-3*mair/(pi*rt2*avog*airden*(dair)**2.)
      wmf=frpath/tr_radius*(s1+s2*exp(-s3*tr_radius/frpath))
      vgs=(1.d0+wmf)*vgs
C****
      return
      end function vgs

      SUBROUTINE read_monthly_sources(iu,jdlast,tlca,tlcb,data,
     *  frac,imon)
!@sum Read in monthly sources and interpolate to current day
!@+   Calling routine must have the lines:
!@+      real*8 tlca(im,jm,nm),tlcb(im,jm,nm)
!@+      integer imon(nm)   ! nm=number of files that will be read
!@+      data jdlast /0/
!@+      save jdlast,tlca,tlcb,imon
!@+   Input: iu, the fileUnit#; jdlast
!@+   Output: interpolated data array + two monthly data arrays
!@auth Jean Lerner and others
      USE FILEMANAGER, only : NAMEUNIT
      USE DOMAIN_DECOMP, only : GRID, GET, AM_I_ROOT
      USE DOMAIN_DECOMP, only : READT_PARALLEL, REWIND_PARALLEL
      USE MODEL_COM, only: jday,im,jm,idofm=>JDmidOfM
      implicit none
      real*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     tlca,tlcb,data
      real*8 :: frac
      integer imon,iu,jdlast

      integer :: J_0, J_1

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

      if (jdlast.EQ.0) then ! NEED TO READ IN FIRST MONTH OF DATA
        imon=1          ! imon=January
        if (jday.le.16)  then ! JDAY in Jan 1-15, first month is Dec
          CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),0,tlca,12)
          CALL REWIND_PARALLEL( iu )
        else            ! JDAY is in Jan 16 to Dec 16, get first month
  120     imon=imon+1
          if (jday.gt.idofm(imon) .AND. imon.le.12) go to 120
          CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),0,tlca,imon-1)
          if (imon.eq.13)  CALL REWIND_PARALLEL( iu )
        end if
      else              ! Do we need to read in second month?
        if (jday.ne.jdlast+1) then ! Check that data is read in daily
          if (jday.ne.1 .OR. jdlast.ne.365) then
            if (AM_I_ROOT()) write(6,*)
     *      'Incorrect values in Tracer Source:JDAY,JDLAST=',JDAY,JDLAST
            call stop_model('stopped in TRACERS.f',255)
          end if
          imon=imon-12  ! New year
          go to 130
        end if
        if (jday.le.idofm(imon)) go to 130
        imon=imon+1     ! read in new month of data
        tlca(:,:) = tlcb(:,:)
        if (imon.eq.13) CALL REWIND_PARALLEL( iu  )
      end if
      CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),0,tlcb,1)
  130 continue
c**** Interpolate two months of data to current day
      frac = float(idofm(imon)-jday)/(idofm(imon)-idofm(imon-1))
      data(:,J_0:J_1) = tlca(:,J_0:J_1)*frac + tlcb(:,J_0:J_1)*(1.-frac)
      return
      end subroutine read_monthly_sources
#endif

      subroutine checktr(subr)
!@sum  CHECKTR Checks whether atmos tracer variables are reasonable
!@auth Gavin Schmidt
!@ver  1.0
#ifdef TRACERS_ON
      USE CONSTANT, only : teeny
      USE MODEL_COM, only : ls1,im,jm,lm,q,wm
      USE GEOM, only : dxyp,imaxj
      USE SOMTQ_COM, only : qmom
      USE DYNAMICS, only : am
      USE FLUXES, only : gtracer
      USE TRACER_COM
      USE DOMAIN_DECOMP, ONLY: GRID, GET, HERE, AM_I_ROOT
      IMPLICIT NONE
      LOGICAL QCHECKT
      INTEGER I,J,L,N,m, imax,jmax,lmax
      REAL*8 relerr, errmax,errsc
!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR
      INTEGER :: J_0, J_1, nj
      CALL GET(GRID, J_STRT=J_0, J_STOP=J_1)
      nj = J_1 - J_0 + 1

      CALL CHECK4(gtracer(1,1,1,J_0),NTM,4,IM,nJ,SUBR,'GTRACE')
      do n=1,ntm
        CALL CHECK4(trmom(:,:,J_0:J_1,:,n),NMOM,IM,nJ,LM,SUBR,
     *       'X'//trname(n))
        CALL CHECK3(trm(:,J_0:J_1,:,n),IM,nJ,LM,SUBR,trname(n))
#ifdef TRACERS_WATER
        CALL CHECK3(trwm(:,J_0:J_1,:,n),IM,nJ,LM,SUBR,'WM'//trname(n))
#endif

C**** check for negative tracer amounts (if t_qlimit is set)
        if (t_qlimit(n)) then
          QCHECKT=.false.
          do l=1,lm
          do j=j_0,j_1
          do i=1,imaxj(j)
            if (trm(i,j,l,n).lt.0) then
              if (AM_I_ROOT())
     *         write(6,*) "Negative mass for ",trname(n),i,j,l,trm(i,j,l
     *             ,n)," after ",SUBR,"."
              QCHECKT=.true.
            end if
          end do
          end do
          end do
          if (QCHECKT)
     &         call stop_model("CHECKTR: Negative tracer amount",255)
        end if

C**** check whether air mass is conserved

        if (trname(n).eq.'Air') then
          errmax = 0. ; lmax=1 ; imax=1 ; jmax=1
          do l=1,lm
          do j=j_0,j_1
          do i=1,imaxj(j)
            relerr=abs(trm(i,j,l,n)-am(l,i,j)*dxyp(j))/
     *           (am(l,i,j)*dxyp(j))
            if (relerr.gt.errmax) then
              lmax=l ; imax=i ; jmax=j ; errmax=relerr
            end if
          end do
          end do
          end do
          print*,"Relative error in air mass after ",trim(subr),":",imax
     *         ,jmax,lmax,errmax,trm(imax,jmax,lmax,n),am(lmax,imax
     *         ,jmax)*dxyp(jmax)
        end if

#ifdef TRACERS_WATER
        if (trname(n).eq.'Water') then
          errmax = 0. ; lmax=1 ; imax=1 ; jmax=1
          do l=1,lm
          do j=j_0,j_1
          do i=1,imaxj(j)
            errsc=(q(i,j,l)+sum(abs(qmom(:,i,j,l))))*am(l,i,j)*dxyp(j)
            if (errsc.eq.0.) errsc=1.
            relerr=abs(trm(i,j,l,n)-q(i,j,l)*am(l,i,j)*dxyp(j))/errsc
            if (wm(i,j,l).gt.0 .and. trwm(i,j,l,n).gt.1.) relerr
     *           =max(relerr,(trwm(i,j,l,n)-wm(i,j,l)*am(l,i,j)*dxyp(j))
     *           /(wm(i,j,l)*am(l,i,j)*dxyp(j)))
            if ((wm(i,j,l).eq.0 .and.trwm(i,j,l,n).gt.1) .or. (wm(i,j,l)
     *           .gt.teeny .and.trwm(i,j,l,n).eq.0))
     *           print*,"Liquid water mismatch: ",subr,i,j,l,trwm(i,j,l
     *           ,n),wm(i,j,l)*am(l,i,j)*dxyp(j)
            do m=1,nmom
              relerr=max(relerr,(trmom(m,i,j,l,n)-qmom(m,i,j,l)*am(l,i,j
     *             )*dxyp(j))/errsc)
            end do
            if (relerr.gt.errmax) then
              lmax=l ; imax=i ; jmax=j ; errmax=relerr
            end if
          end do
          end do
          end do
          print*,"Relative error in water mass after ",trim(subr),":"
     *         ,imax,jmax,lmax,errmax,trm(imax,jmax,lmax,n),q(imax,jmax
     *         ,lmax)*am(lmax,imax,jmax)*dxyp(jmax),trwm(imax,jmax,lmax
     *         ,n),wm(imax,jmax,lmax)*am(lmax,imax,jmax)*dxyp(jmax)
     *         ,(trmom(m,imax,jmax,lmax,n),qmom(m,imax,jmax,lmax)
     *         *am(lmax,imax,jmax)*dxyp(jmax),m=1,nmom)
        end if
#endif
      end do
#endif
      return
      end subroutine checktr


      SUBROUTINE io_tracer(kunit,iaction,ioerr)
!@sum  io_tracer reads and writes tracer variables to file
!@auth Jean Lerner
!@ver  1.0
#ifdef TRACERS_ON
      USE MODEL_COM, only: ioread,iowrite,irsfic,irsficno,irerun,lhead
      USE DOMAIN_DECOMP, only : grid, AM_I_ROOT, PACK_DATA, UNPACK_DATA
     &,PACK_DATAj, UNPACK_DATAj, PACK_BLOCK, UNPACK_BLOCK, PACK_COLUMN
     &,UNPACK_COLUMN, esmf_bcast
      USE TRACER_COM
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM, only: yNO3,pHOx,pNOx,pOx,yCH3O2,yC2O3,
     & yROR,yXO2,yAldehyde,yXO2N,yRXPAR,ss,corrOx,JPPJ
#ifdef SHINDELL_STRAT_CHEM
     & ,SF3,SF2,pClOx,pClx,pOClOx,pBrOx,yCl2,yCl2O2
#endif
#ifdef INTERACTIVE_WETLANDS_CH4 
      use TRACER_SOURCES, only: day_ncep,DRA_ch4,sum_ncep,PRS_ch4,
     & HRA_ch4,iday_ncep,i0_ncep,iHch4,iDch4,i0ch4,first_ncep,first_mod
     & ,max_days,nra_ncep,nra_ch4,maxHR_ch4
#endif
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      USE fluxes,ONLY : nstype,pprec,pevap
      USE tracers_dust,ONLY : hbaij,ricntd
#endif
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records

      REAL*8, DIMENSION(NMOM,IM,JM,LM,NTM) :: TRMOM_GLOB
      REAL*8, DIMENSION(     IM,JM,LM,NTM) :: TRM_GLOB
            
#ifdef TRACERS_SPECIAL_Shindell
      REAL*8, DIMENSION(IM,JM,LM) :: yNO3_glob,pHOx_glob,pNOx_glob,
     &pOx_glob,yCH3O2_glob,yC2O3_glob,yROR_glob,yXO2_glob,
     &yAldehyde_glob,yXO2N_glob,yRXPAR_glob
#ifdef SHINDELL_STRAT_CHEM
     &,SF3_glob,SF2_glob,pClOx_glob,pClx_glob,pOClOx_glob,
     &pBrOx_glob,yCl2_glob,yCl2O2_glob
#endif     
      REAL*8, DIMENSION(JM,LM,12) :: corrOx_glob
      REAL*8, DIMENSION(JPPJ,LM,IM,JM) :: ss_glob     
#ifdef INTERACTIVE_WETLANDS_CH4 
      REAL*8, dimension(IM,JM,max_days,nra_ncep) :: day_ncep_glob
      REAL*8, dimension(IM,JM,max_days,nra_ch4)  :: DRA_ch4_glob
      REAL*8, dimension(IM,JM,nra_ncep)          :: sum_ncep_glob
      REAL*8, dimension(IM,JM,nra_ch4)           :: PRS_ch4_glob
      REAL*8, dimension(IM,JM,maxHR_ch4,nra_ch4) :: HRA_ch4_glob
      INTEGER, dimension(IM,JM,nra_ch4) :: iHch4_glob,iDch4_glob,
     &                                     i0ch4_glob,first_mod_glob
#endif
#endif     
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      REAL*8,DIMENSION(Im,Jm) :: pprec_glob,ricntd_glob,hbaij_glob
      REAL*8,DIMENSION(Im,Jm,Nstype) :: pevap_glob
#endif
      INTEGER :: ITM,ITM1,ITM2
#ifdef TRACERS_WATER
      CHARACTER*80 :: HEADER, MODULE_HEADER = "TRACERW01"
      REAL*8, DIMENSION(     IM,JM,LM,NTM) :: TRWM_GLOB

      write (MODULE_HEADER(lhead+1:80),'(a,i2,a,a,i1,a,i2,a,i2,a)')
     *     'R8 TRM(im,jm,lm,',NTM,')',
     *     ',TRmom(',NMOM,',im,jm,lm,',NTM,'),trwm(im,jm,lm,',NTM,')'
#else
      CHARACTER*80 :: HEADER, MODULE_HEADER = "TRACER01"

      write (MODULE_HEADER(lhead+1:80),'(a,i2,a,a,i1,a,i2,a)')
     *           'R8 TRM(im,jm,lm,',NTM,')',
     *  ',TRmom(',NMOM,',im,jm,lm,',NTM,')'
#endif

      SELECT CASE (IACTION)

      CASE (:IOWRITE) ! output to end-of-month restart file
      
#ifdef TRACERS_SPECIAL_Shindell
       CALL PACK_DATA(grid, yNO3     , yNO3_glob)
       CALL PACK_DATA(grid, pHOx     , pHOx_glob)
       CALL PACK_DATA(grid, pNOx     , pNOx_glob)
       CALL PACK_DATA(grid, pOx      , pOx_glob)
       CALL PACK_DATA(grid, yCH3O2   , yCH3O2_glob)
       CALL PACK_DATA(grid, yC2O3    , yC2O3_glob)
       CALL PACK_DATA(grid, yROR     , yROR_glob)
       CALL PACK_DATA(grid, yXO2     , yXO2_glob)
       CALL PACK_DATA(grid, yAldehyde, yAldehyde_glob)
       CALL PACK_DATA(grid, yXO2N    , yXO2N_glob)
       CALL PACK_DATA(grid, yRXPAR   , yRXPAR_glob)  
       DO ITM=1,12
         CALL PACK_DATAj(grid, corrOx(:,:,ITM), corrOx_glob(:,:,ITM))
       END DO
       CALL PACK_BLOCK(grid, ss(:,:,:,:), ss_glob(:,:,:,:))
#ifdef SHINDELL_STRAT_CHEM
       CALL PACK_DATA(grid, SF3      , SF3_glob)
       CALL PACK_DATA(grid, SF2      , SF2_glob)
       CALL PACK_DATA(grid, pClOx    , pClOx_glob)
       CALL PACK_DATA(grid, pClx     , pClx_glob)
       CALL PACK_DATA(grid, pOClOx   , pOClOx_glob)
       CALL PACK_DATA(grid, pBrOx    , pBrOx_glob)
       CALL PACK_DATA(grid, yCl2     , yCl2_glob)
       CALL PACK_DATA(grid, yCl2O2   , yCl2O2_glob)
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
       DO ITM=1,max_days
        DO ITM2=1,nra_ncep
         CALL PACK_DATA(grid,day_ncep(:,:,ITM,ITM2),
     &   day_ncep_glob(:,:,ITM,ITM2))
        END DO
        DO ITM2=1,nra_ch4
         CALL PACK_DATA(grid,DRA_ch4(:,:,ITM,ITM2),
     &   DRA_ch4_glob(:,:,ITM,ITM2))
        END DO
       END DO 
       DO ITM1=1,nra_ncep
        CALL PACK_DATA(grid,sum_ncep(:,:,ITM1),sum_ncep_glob(:,:,ITM1))
       END DO 
       DO ITM1=1,nra_ch4
        CALL PACK_DATA(grid,PRS_ch4(:,:,ITM1),PRS_ch4_glob(:,:,ITM1))
       END DO
       DO ITM=1,maxHR_ch4
        DO ITM2=1,nra_ch4
         CALL PACK_DATA(grid,HRA_ch4(:,:,ITM,ITM2),
     &   HRA_ch4_glob(:,:,ITM,ITM2))
        END DO
       END DO       
       CALL PACK_DATA(grid,iHch4(:,:,:),iHch4_glob(:,:,:)) !integer
       CALL PACK_DATA(grid,iDch4(:,:,:),iDch4_glob(:,:,:)) !integer
       CALL PACK_DATA(grid,i0ch4(:,:,:),i0ch4_glob(:,:,:)) !integer
       CALL PACK_DATA(grid,first_mod(:,:,:),first_mod_glob(:,:,:)) !integer
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
       CALL pack_data(grid,pprec,pprec_glob)
       CALL pack_data(grid,pevap,pevap_glob)
       CALL pack_data(grid,hbaij,hbaij_glob)
       CALL pack_data(grid,ricntd,ricntd_glob)
#endif
#endif
        DO ITM=1,NTM
          CALL PACK_DATA(grid, TRM(:,:,:,ITM), TRM_GLOB(:,:,:,ITM))
          CALL PACK_COLUMN(grid, TRmom(:,:,:,:,ITM),
     &         TRmom_GLOB(:,:,:,:,ITM))
#ifdef TRACERS_WATER
          CALL PACK_DATA(grid, TRWM(:,:,:,ITM), TRWM_GLOB(:,:,:,ITM))
#endif
        END DO
        IF (AM_I_ROOT()) THEN
          WRITE (kunit,err=10) MODULE_HEADER,TRM_glob,TRmom_glob
#ifdef TRACERS_WATER
     *       ,TRWM_glob
#endif
#ifdef TRACERS_SPECIAL_Shindell
     *       ,yNO3_glob,pHOx_glob,pNOx_glob,pOx_glob,yCH3O2_glob
     *       ,yC2O3_glob,yROR_glob,yXO2_glob,yAldehyde_glob
     *       ,yXO2N_glob,yRXPAR_glob,corrOx_glob,ss_glob
#ifdef SHINDELL_STRAT_CHEM
     *       ,SF3_glob,SF2_glob,pClOx_glob,pClx_glob,pOClOx_glob
     *       ,pBrOx_glob,yCl2_glob,yCl2O2_glob
#endif
#ifdef INTERACTIVE_WETLANDS_CH4 
     *       ,day_ncep_glob,DRA_ch4_glob,sum_ncep_glob,PRS_ch4_glob
     *       ,HRA_ch4_glob,iday_ncep    ,i0_ncep     ,iHch4_glob
     *       ,iDch4_glob,i0ch4_glob,first_ncep     ,first_mod_glob
#endif
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,hbaij_glob,ricntd_glob,pprec_glob,pevap_glob
#endif
       END IF     !only root processor writes

      CASE (IOREAD:)          ! input from restart file
        SELECT CASE (IACTION)
        CASE (ioread,irerun,irsfic,irsficno) ! restarts
          if ( AM_I_ROOT() ) then
            READ (kunit,err=10) HEADER,TRM_glob,TRmom_glob
#ifdef TRACERS_WATER
     *       ,TRWM_glob
#endif
#ifdef TRACERS_SPECIAL_Shindell
     *       ,yNO3_glob,pHOx_glob,pNOx_glob,pOx_glob,yCH3O2_glob
     *       ,yC2O3_glob,yROR_glob,yXO2_glob,yAldehyde_glob
     *       ,yXO2N_glob,yRXPAR_glob,corrOx_glob,ss_glob
#ifdef SHINDELL_STRAT_CHEM
     *       ,SF3_glob,SF2_glob,pClOx_glob,pClx_glob,pOClOx_glob
     *       ,pBrOx_glob,yCl2_glob,yCl2O2_glob
#endif
#ifdef INTERACTIVE_WETLANDS_CH4 
     *       ,day_ncep_glob,DRA_ch4_glob,sum_ncep_glob,PRS_ch4_glob
     *       ,HRA_ch4_glob,iday_ncep     ,i0_ncep     ,iHch4_glob
     *       ,iDch4_glob,i0ch4_glob,first_ncep     ,first_mod_glob
#endif
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     ,hbaij_glob,ricntd_glob,pprec_glob,pevap_glob
#endif
            IF (HEADER(1:lhead).ne.MODULE_HEADER(1:lhead)) THEN
              PRINT*,"Discrepancy in module version ",HEADER,
     &                MODULE_HEADER
              GO TO 10
            END IF
          end if      ! AM_I_ROOT
C**** ESMF: Copy global read data into the corresponding local (distributed) arrays.
       DO ITM=1,NTM
         CALL UNPACK_DATA  (grid,   TRM_GLOB(:,:,:,  itm),
     &                              TRM(:,:,:,  itm))
         CALL UNPACK_COLUMN(grid, TRMOM_GLOB(:,:,:,:,itm),
     &                              TRmom(:,:,:,:,itm))
#ifdef TRACERS_WATER
         CALL UNPACK_DATA  (grid,  TRWM_GLOB(:,:,:,  itm),
     &                              TRWM(:,:,:,  itm))
#endif
        END DO
#ifdef TRACERS_SPECIAL_Shindell
       CALL UNPACK_DATA(grid, yNO3_glob     , yNO3)
       CALL UNPACK_DATA(grid, pHOx_glob     , pHOx)
       CALL UNPACK_DATA(grid, pNOx_glob     , pNOx)
       CALL UNPACK_DATA(grid, pOx_glob      , pOx)
       CALL UNPACK_DATA(grid, yCH3O2_glob   , yCH3O2)
       CALL UNPACK_DATA(grid, yC2O3_glob    , yC2O3)
       CALL UNPACK_DATA(grid, yROR_glob     , yROR)
       CALL UNPACK_DATA(grid, yXO2_glob     , yXO2)
       CALL UNPACK_DATA(grid, yAldehyde_glob, yAldehyde)
       CALL UNPACK_DATA(grid, yXO2N_glob    , yXO2N)
       CALL UNPACK_DATA(grid, yRXPAR_glob   , yRXPAR)  
       DO ITM=1,12
         CALL UNPACK_DATAj(grid, corrOx_glob(:,:,ITM), corrOx(:,:,ITM))
       END DO
       CALL UNPACK_BLOCK(grid, ss_glob(:,:,:,:), ss(:,:,:,:))
#ifdef SHINDELL_STRAT_CHEM
       CALL UNPACK_DATA(grid, SF3_glob      , SF3)
       CALL UNPACK_DATA(grid, SF2_glob      , SF2)
       CALL UNPACK_DATA(grid, pClOx_glob    , pClOx)
       CALL UNPACK_DATA(grid, pClx_glob     , pClx)
       CALL UNPACK_DATA(grid, pOClOx_glob   , pOClOx)
       CALL UNPACK_DATA(grid, pBrOx_glob    , pBrOx)
       CALL UNPACK_DATA(grid, yCl2_glob     , yCl2)
       CALL UNPACK_DATA(grid, yCl2O2_glob   , yCl2O2)
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
       DO ITM=1,max_days
        DO ITM2=1,nra_ncep
         CALL UNPACK_DATA(grid,day_ncep_glob(:,:,ITM,ITM2),
     &   day_ncep(:,:,ITM,ITM2))
        END DO
        DO ITM2=1,nra_ch4
         CALL UNPACK_DATA(grid,DRA_ch4_glob(:,:,ITM,ITM2),
     &   DRA_ch4(:,:,ITM,ITM2))
        END DO
       END DO 
       DO ITM1=1,nra_ncep
        CALL UNPACK_DATA(grid,sum_ncep_glob(:,:,ITM1),
     &  sum_ncep(:,:,ITM1))
       END DO 
       DO ITM1=1,nra_ch4
        CALL UNPACK_DATA(grid,PRS_ch4_glob(:,:,ITM1),PRS_ch4(:,:,ITM1))
       END DO
       DO ITM=1,maxHR_ch4
        DO ITM2=1,nra_ch4
         CALL UNPACK_DATA(grid,HRA_ch4_glob(:,:,ITM,ITM2),
     &   HRA_ch4(:,:,ITM,ITM2))
        END DO
       END DO       
       CALL UNPACK_DATA(grid,iHch4_glob(:,:,:),iHch4(:,:,:)) !integer
       CALL UNPACK_DATA(grid,iDch4_glob(:,:,:),iDch4(:,:,:)) !integer
       CALL UNPACK_DATA(grid,i0ch4_glob(:,:,:),i0ch4(:,:,:)) !integer
       CALL UNPACK_DATA(grid,first_mod_glob(:,:,:),first_mod(:,:,:)) !integer
#endif
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
        CALL unpack_data(grid,hbaij_glob,hbaij)
        CALL unpack_data(grid,ricntd_glob,ricntd)
        CALL unpack_data(grid,pprec_glob,pprec)
        CALL unpack_data(grid,pevap_glob,pevap)
#endif
C**** ESMF: Broadcast all non-distributed read arrays.
#ifdef INTERACTIVE_WETLANDS_CH4
        call ESMF_BCAST( grid, iday_ncep )
        call ESMF_BCAST( grid, i0_ncep   )
        call ESMF_BCAST( grid, first_ncep)
#endif
        END SELECT
      END SELECT

      RETURN
 10   IOERR=1
#endif
      RETURN
      END SUBROUTINE io_tracer


      subroutine num_srf_sources(nt,nsrc)
!@sum reads headers from emission files to return
!@+ source names and determine the number of sources
!@+ from the number of files in the rundeck of the form:
!@+ trname_##. Then assigns each source to sector(s),
!@+ based on definitions in the rundeck.
!@auth Greg Faluvegi

      use TRACER_COM, only : ntsurfsrcmax,trname,
     & num_tr_sectors,n_max_sect,tr_sect_index,num_sectors,
     & tr_sect_name, sect_name
      USE DOMAIN_DECOMP, only: GRID, GET, write_parallel
      USE FILEMANAGER, only: openunit,closeunit
      USE PARAM, only : sync_param
    
      implicit none

!@var nsrc number of source to define ntsurfsrc(n)
      integer, intent(out) :: nsrc
      integer, intent(in) :: nt
      integer :: n,iu,i,j,nsect,nn
      character*2 :: fnum
      character*80 :: fname
      character*32 :: pname
      character*124 :: tr_sectors_are
      character(len=300) :: out_line
      logical :: qexist

! loop through potential number of surface sources, checking if
! those files exist. If they do, obtain the source name by reading
! the header. If not, the number of sources for this tracer has 
! been reached.

      nsrc=0
      loop_n: do n=1,ntsurfsrcmax
        if(n < 10) then ; write(fnum,'(a1,I1)')'0',n
        else ; write(fnum,'(I2)')n ; endif
        fname=trim(trname(nt))//'_'//fnum
        inquire(file=fname,exist=qexist)
        if(qexist) then
          nsrc=nsrc+1
          call openunit(fname,iu,.true.)
          call read_emis_header(nt,n,iu)
          call closeunit(iu)
! -- begin sector  stuff --
          tr_sectors_are = ' '
          pname=trim(trim(fname)//'_sect')
          call sync_param(pname,tr_sectors_are)
          num_tr_sectors(nt,n)=0
          i=1
          do while(i < len(tr_sectors_are))
            j=index(tr_sectors_are(i:len(tr_sectors_are))," ")
            if (j > 1) then
              num_tr_sectors(nt,n)=num_tr_sectors(nt,n)+1
              i=i+j
            else
              i=i+1
            end if
          enddo
          if(num_tr_sectors(nt,n) > n_max_sect) call stop_model
     &    ("num_tr_sectors problem",255)
          if(num_tr_sectors(nt,n) > 0)then
            read(tr_sectors_are,*)
     &      tr_sect_name(nt,n,1:num_tr_sectors(nt,n))
            do nsect=1,num_tr_sectors(nt,n)
              tr_sect_index(nt,n,nsect)=0
              loop_nn: do nn=1,num_sectors
                if(trim(tr_sect_name(nt,n,nsect)) ==
     &          trim(sect_name(nn))) then
                  tr_sect_index(nt,n,nsect)=nn
                  exit loop_nn
                endif
              enddo loop_nn
            enddo
          endif
! -- end sector stuff --
        else
          exit loop_n
        endif
      enddo loop_n

! and make sure there isn't a skip:

      n=n+1
      if(n < 10) then ; write(fnum,'(a1,I1)')'0',n
      else ; write(fnum,'(I2)')n ; endif
      fname=trim(trname(nt))//'_'//fnum
      inquire(file=fname,exist=qexist)
      if(qexist) then
        write(out_line,*)'problem with num_srf_sources'
        call write_parallel(trim(out_line))
        call stop_model(trim(out_line),255)
      endif

      end subroutine num_srf_sources


      subroutine read_sfc_sources(n,nsrc)
!@sum reads surface (2D generally non-interactive) sources
!@auth Jean Lerner/Greg Faluvegi
      USE MODEL_COM, only: itime,jday,jyear,im,jm
      USE DOMAIN_DECOMP, only: GRID, GET, readt_parallel, write_parallel
      USE FILEMANAGER, only: openunit,closeunit, nameunit
      USE TRACER_COM,only:itime_tr0,trname,sfc_src,ntm,ntsurfsrcmax,
     & freq,nameT,ssname,ty_start,ty_end,kstep
   
      implicit none
      
      integer :: iu,ns,k,ipos
      integer, intent(in) :: nsrc,n
      character*80 :: fname
      character*2 :: fnum
      character(len=300) :: out_line
      logical,dimension(ntm,ntsurfsrcmax) :: ifirst2=.true.
      real*8 :: alpha
      real*8, dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     & sfc_a,sfc_b
      
      save ifirst2

      INTEGER :: J_1, J_0, J_0H, J_1H

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

      if (itime < itime_tr0(n)) return
      if (nsrc <= 0) return
   
      call GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

      do ns=1,nsrc

! open file and read its header:

        if(ns < 10) then
          write(fnum,'(a1,I1)')'0',ns
        else
          write(fnum,'(I2)')ns
        endif
        fname=trim(trname(n))//'_'//fnum
        call openunit(fname,iu,.true.)
        call read_emis_header(n,ns,iu) ! maybe just ROOT should do?

! now read the data: (should be in kg/m2/s please)
 
! -------------- non-transient emissions ----------------------------!
        if(ty_start(n,ns)==ty_end(n,ns))then 

          select case(freq(n,ns))
          case('a')        ! annual file, only read first time
            if(ifirst2(n,ns)) then
              call readt_parallel(grid,iu,fname,0,sfc_src(:,:,n,ns),1)
              write(out_line,*)trim(nameT(n,ns)),
     &        ' ',trim(ssname(n,ns)),' ann source read.'
              call write_parallel(trim(out_line))
              ifirst2(n,ns) = .false.
            endif
          case('m')        ! monthly file, interpolate to now
            call read_mon_src_2(n,ns,iu,sfc_src(:,:,n,ns))
            ifirst2(n,ns) = .false. ! needed?
          end select

! --------------- transient emissions -------------------------------!
        else                
          select case(freq(n,ns))
          case('a')        ! annual file, only read first time + new steps
            ipos=1
            alpha=0.d0 ! before start year, use start year value
            if(jyear>ty_end(n,ns).or.
     &      (jyear==ty_end(n,ns).and.jday>=183))then
              alpha=1.d0 ! after end year, use end year value     
              ipos=(ty_end(n,ns)-ty_start(n,ns))/kstep
            endif
            do k=ty_start(n,ns),ty_end(n,ns)-kstep,kstep
!should do!   if(jyear==k .and. jday==183)ifirst2(n,ns)=.true.
              if(jyear>k .or. (jyear==k.and.jday>=183)) then
                if(jyear<k+kstep.or.(jyear==k+kstep.and.jday<183))then
                  ipos=1+(k-ty_start(n,ns))/kstep ! (integer artithmatic)
                  alpha=(365.d0*(0.5+real(jyear-1-k))+jday) / 
     &                  (365.d0*real(kstep))
                  exit
                endif
              endif
            enddo
!should do! if(ifirst2(n,ns)) then
              call readt_parallel(grid,iu,fname,0,sfc_a(:,:),ipos)
              call readt_parallel(grid,iu,fname,0,sfc_b(:,:),1)
!should do! endif
            sfc_src(:,:,n,ns)=sfc_a(:,:)*(1.d0-alpha)+sfc_b(:,:)*alpha
            if(alpha>0.d0)then 
              write(out_line,*)
     &        trim(nameT(n,ns)),' ',trim(ssname(n,ns)),' at ',
     &        100.d0*alpha,' % of period ',k,' to ',k+kstep
              call write_parallel(trim(out_line))
            endif
            ifirst2(n,ns) = .false.

          case('m')        ! monthly file, interpolate to now
            call read_mon_src_2(n,ns,iu,sfc_src(:,:,n,ns))
            ifirst2(n,ns) = .false. ! needed?
          end select

        endif
        call closeunit(iu)

      enddo

      return
      end subroutine read_sfc_sources


      subroutine read_emis_header(n,ns,iu)
!@sum read_emis_header reads the emissions file's header and
!@+   reports back the meta-data. 
!@auth Greg Faluvegi
     
      use TRACER_COM, only : trname,nameT
      USE DOMAIN_DECOMP, only: GRID, GET, write_parallel
      USE FILEMANAGER, only: openunit,closeunit

      implicit none

      integer, intent(in) :: iu,n,ns
      integer :: error
      character*80 :: message,header
      character(len=300) :: out_line

      error=0
      read(iu)header
      if(len_trim(header) < 1)error=1

      call parse_header(n,ns,trim(header),error)
  
      select case(error)
      case default ! nothing
      case(1) ; message='read_emis_header: missing header'
      case(2) ; message='read_emis_header: problem with freq'
      case(3) ; message='read_emis_header: a and m are choices for freq'
      case(4) ; message='read_emis_header: problem with tracer name'
      case(5) ; message='read_emis_header: tracer name mismatch'
      case(6) ; message='read_emis_header: problem with source'
      case(7) ; message='read_emis_header: problem with res'
      case(8) ; message='read_emis_header: M and F are choices for res'
      case(9) ; message='read_emis_header: transient years seem wrong'
      case(10); message='read_emis_header: trans yrs not 10 years apart'
      end select
      if(error > 0) then
        write(out_line,*) trim(header)
        call write_parallel(trim(out_line))
        write(out_line,*) trim(message)
        call write_parallel(trim(out_line))
        call stop_model('problem reading emissions',255)
      endif

      end subroutine read_emis_header


      subroutine parse_header(n,ns,str,error)
!@sum parse_header gets the informantion from emissions file's
!@+  header and reports back the meta-data. 
!@auth Greg Faluvegi
     
      use TRACER_COM, only : trname,freq,nameT,res,ssname,Tyears,
     & ty_start,ty_end

      implicit none

      integer, intent(in) :: n,ns
      integer :: n1,n2,error
      character*(*) str

      n1 = scan( str,'=')         ! freq
      if(str(1:n1-1) /= 'freq')error=2
      read(str(n1+1:n1+1),*)freq(n,ns)
      if(freq(n,ns) /= 'a' .and. freq(n,ns) /= 'm')error=3

      str = str(n1+3:)            ! tracer name
      n1 = scan( str, '=')
      if(str(1:n1-1) /= 'name')error=4
      n2 = scan( str,' ')
      read(str(n1+1:n2-1),*)nameT(n,ns)
      if(trim(trname(n)) /= trim(nameT(n,ns)))error=5

      str = str(n2+1:)            ! source name
      n1 = scan( str,'=')
      if(str(1:n1-1) /= 'source')error=6
      n2 = scan( str,' ')
      read(str(n1+1:n2-1),*)ssname(n,ns)

      str = str(n2+1:)            ! resolution
      n1 = scan( str,'=')
      if(str(1:n1-1) /= 'res')error=7
      read(str(n1+1:n1+1),*)res(n,ns)
      if(res(n,ns) /= 'M' .and. res(n,ns) /= 'F')error=8

      str = str(n1+3:)            ! optional transient years
      n1 = scan( str, '=')
      if(str(1:n1-1) /= 'y')then
        ty_start(n,ns)=0; ty_end(n,ns)=0
      else 
        read(str(n1+1:n1+9),*)Tyears(n,ns)
        read(Tyears(n,ns)(1:4),'(I4)')ty_start(n,ns)
        read(Tyears(n,ns)(6:9),'(I4)')ty_end(n,ns)
        if(ty_start(n,ns) /= ty_end(n,ns))then
          if(ty_start(n,ns) < 0 .or. ty_start(n,ns) > 3000)error=9
          if(ty_end(n,ns)   < 0 .or. ty_end(n,ns)   > 3000)error=9
          if(ty_end(n,ns)-ty_start(n,ns) < 10)error=10 ! decades expected
        endif
      endif

      end subroutine parse_header


      SUBROUTINE read_mon_src_2(n,ns,iu,data)
!@sum Read in monthly sources and interpolate to current day
!@+   Calling routine must have the lines:
!@+      integer imon(nm)   ! nm=number of files that will be read
!@+      data jdlast /0/
!@+      save jdlast,imon
!@+   Input: iu, the fileUnit#; jdlast
!@+   Output: interpolated data array + two monthly data arrays
!@+ Note: this started out very similar to read_monthly_sources, 
!@+ but I had problems with that routine. This one is less effecient,
!@+ in that it doesn't keep file open, and has to read monthly files
!@+ each day. My intent is to eventually fix this problem and merge
!@+ with read_monthly_sources. --gsf
!@auth Greg Faluvegi, Jean Lerner and others

      USE FILEMANAGER, only : NAMEUNIT
      USE DOMAIN_DECOMP, only : GRID, GET, AM_I_ROOT, write_parallel,
     & READT_PARALLEL, REWIND_PARALLEL, BACKSPACE_PARALLEL
      USE MODEL_COM, only: jday,im,jm,idofm=>JDmidOfM,jyear
      USE TRACER_COM, only: ssname,nameT,ty_start,ty_end,kstep

      implicit none

      real*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     tlca,tlcb,data
      real*8 :: frac,alpha
      integer ::  imon,iu,n,ns,ipos,k,nn
      character*80 :: junk
      character(len=300) :: out_line
      real*8, dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     & sfc_a,sfc_b

      integer :: J_0, J_1

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

! -------------- non-transient emissions ----------------------------!
      if(ty_start(n,ns)==ty_end(n,ns))then
        imon=1
        if (jday <= 16)  then ! JDAY in Jan 1-15, first month is Dec
          call readt_parallel(grid,iu,nameunit(iu),0,tlca,12)
          call rewind_parallel( iu );if (AM_I_ROOT()) read( iu ) junk
        else            ! JDAY is in Jan 16 to Dec 16, get first month
          do while(jday > idofm(imon) .AND. imon <= 12)
            imon=imon+1
          enddo
          call readt_parallel(grid,iu,nameunit(iu),0,tlca,imon-1)
          if (imon == 13)then
            call rewind_parallel( iu );if (AM_I_ROOT()) read( iu ) junk
          endif
        end if
        call readt_parallel(grid,iu,nameunit(iu),0,tlcb,1)

c****   Interpolate two months of data to current day
        frac = float(idofm(imon)-jday)/(idofm(imon)-idofm(imon-1))
        data(:,J_0:J_1)=tlca(:,J_0:J_1)*frac+tlcb(:,J_0:J_1)*(1.-frac)
        write(out_line,*)
     &  trim(nameT(n,ns)),' ',trim(ssname(n,ns)),' interp to now ',frac
        call write_parallel(trim(out_line))
! --------------- transient emissions -------------------------------!
      else 
        ipos=1
        alpha=0.d0 ! before start year, use start year value
        if(jyear>ty_end(n,ns).or.
     &  (jyear==ty_end(n,ns).and.jday>=183))then
          alpha=1.d0 ! after end year, use end year value
          ipos=(ty_end(n,ns)-ty_start(n,ns))/kstep
        endif
        do k=ty_start(n,ns),ty_end(n,ns)-kstep,kstep
          if(jyear>k .or. (jyear==k.and.jday>=183)) then
            if(jyear<k+kstep .or. (jyear==k+kstep.and.jday<183))then
              ipos=1+(k-ty_start(n,ns))/kstep ! (integer artithmatic)
              alpha=(365.d0*(0.5+real(jyear-1-k))+jday) /
     &              (365.d0*real(kstep))
              exit
            endif
          endif
        enddo
!
! read the two necessary months from the first decade:
!
        imon=1
        if (jday <= 16)  then ! JDAY in Jan 1-15, first month is Dec
         call readt_parallel(grid,iu,nameunit(iu),0,tlca,(ipos-1)*12+12)
         do nn=1,12; call backspace_parallel(iu); enddo
        else            ! JDAY is in Jan 16 to Dec 16, get first month
         do while(jday > idofm(imon) .AND. imon <= 12)
           imon=imon+1
         enddo
         call  readt_parallel
     &   (grid,iu,nameunit(iu),0,tlca,(ipos-1)*12+imon-1)
         if (imon == 13)then
           do nn=1,12; call backspace_parallel(iu); enddo               
         endif
        end if
        call readt_parallel(grid,iu,nameunit(iu),0,tlcb,1)
c****   Interpolate two months of data to current day
        frac = float(idofm(imon)-jday)/(idofm(imon)-idofm(imon-1))
        sfc_a(:,J_0:J_1)=tlca(:,J_0:J_1)*frac+tlcb(:,J_0:J_1)*(1.-frac)
        call rewind_parallel( iu );if (AM_I_ROOT()) read( iu ) junk

        ipos=ipos+1
        imon=1
        if (jday <= 16)  then ! JDAY in Jan 1-15, first month is Dec
         call readt_parallel(grid,iu,nameunit(iu),0,tlca,(ipos-1)*12+12)
         do nn=1,12; call backspace_parallel(iu); enddo
        else            ! JDAY is in Jan 16 to Dec 16, get first month
         do while(jday > idofm(imon) .AND. imon <= 12)
           imon=imon+1
         enddo
         call readt_parallel
     &   (grid,iu,nameunit(iu),0,tlca,(ipos-1)*12+imon-1)
         if (imon == 13)then
           do nn=1,12; call backspace_parallel(iu); enddo
         endif
        end if
        call readt_parallel(grid,iu,nameunit(iu),0,tlcb,1)
c****   Interpolate two months of data to current day
        frac = float(idofm(imon)-jday)/(idofm(imon)-idofm(imon-1))
        sfc_b(:,J_0:J_1)=tlca(:,J_0:J_1)*frac+tlcb(:,J_0:J_1)*(1.-frac)

! now interpolate between the two time periods:

        data(:,:)=sfc_a(:,:)*(1.d0-alpha)+sfc_b(:,:)*alpha

        if(alpha>0.d0)then 
          write(out_line,*)
     &    trim(nameT(n,ns)),' ',trim(ssname(n,ns)),' at ',
     &    100.d0*alpha,' % of period ',k,' to ',k+kstep,
     &    ' and monthly fraction= ',frac
          call write_parallel(trim(out_line))
        endif
 
      endif

      return
      end subroutine read_mon_src_2


      subroutine setup_emis_sectors_regions
!@sum setup_emis_sectors_regions reads from the rundeck the 
!@+ geographic regions and sectors associated with tracer
!@+ emissions and saves names.
!@+ Also reads the factors associated with each sector and
!@+ region. Output IJ map of regions.
!@auth Greg Faluvegi

      use TRACER_COM, only : n_max_sect,reg_n,reg_s,reg_e,reg_w,
     & n_max_reg,alter_sources,ef_REG_IJ,ef_REG_IJ_glob,
     & ef_REG_IJ_glob4,ef_fact,num_regions,num_sectors,sect_name
      USE DOMAIN_DECOMP, only:GRID,GET,AM_I_ROOT,PACK_DATA,UNPACK_DATA,
     & UNPACK_DATAj,write_parallel
      USE GEOM, only: lat_dg, lon_dg, imaxj
      USE FILEMANAGER, only: openunit,closeunit
      USE PARAM, only : sync_param

      implicit none

      integer :: i,j,n,iu
      character*80 :: title
      character*2 :: fnum
      character*124 :: sectors_are,regions_are

      INTEGER J_0, J_1
      CALL GET(grid, J_STRT=J_0,J_STOP=J_1)

      sectors_are=' '
      regions_are=' '
      call sync_param("sectors_are",sectors_are)
      call sync_param("regions_are",regions_are)

! see how many sectors there are, save names in array:

      num_sectors=0
      i=1
      do while(i < len(sectors_are))
        j=index(sectors_are(i:len(sectors_are))," ")
        if (j > 1) then
          num_sectors=num_sectors+1
          i=i+j
        else
          i=i+1
        end if
      enddo
      if (num_sectors > n_max_sect) call stop_model
     &("n_max_sect must be increased",255)
      if(num_sectors > 0 ) read(sectors_are,*)
     & sect_name(1:num_sectors)

! see how many regions there are, save names in array:

      num_regions=0
      i=1
      do while(i < len(regions_are))
        j=index(regions_are(i:len(regions_are))," ")
        if (j > 1) then
          num_regions=num_regions+1
          i=i+j
        else
          i=i+1
        end if
      enddo
      if (num_regions > n_max_reg) call stop_model
     &("n_max_reg must be increased",255)
      if(num_regions>0) then
        REG_N(:)=-1.e30 ! initialize to undefined
        REG_S(:)=-1.e30; REG_E(:)=-1.e30; REG_W(:)=-1.e30
        call sync_param("REG_N",REG_N,num_regions)
        call sync_param("REG_S",REG_S,num_regions)
        call sync_param("REG_E",REG_E,num_regions)
        call sync_param("REG_W",REG_W,num_regions)
      end if 

! read the actual emission altering factors:
      do n=1,num_sectors
        if(n < 10) then ; write(fnum,'(a1,I1)')'0',n
        else ; write(fnum,'(I2)')n ; endif
        ef_fact(n,:)=1.d0
        call sync_param("SECT_"//fnum,ef_fact(n,:),num_regions)
      enddo

! check if there is any altering requested:
      alter_sources = .false.
      do i=1,num_sectors; do j=1,num_regions
        if(ef_fact(i,j)/=1.0 .and. ef_fact(i,j)/=-1.d30)
     &  alter_sources = .true.
      enddo ; enddo

! write out regions to GISS-format IJ file, each restart:
      if(alter_sources)then
        ef_REG_IJ(:,J_0:J_1)=0.d0
        do j=J_0,J_1; do i=1,imaxj(j); do n=1,num_regions
          if(lat_dg(j,1) >= REG_S(n) .and. lat_dg(j,1)
     &    <= REG_N(n) .and. lon_dg(i,1) >= REG_W(n)
     &    .and. lon_dg(i,1) < REG_E(n) ) ef_REG_IJ(i,j)=
     &    max(ef_REG_IJ(i,j),dble(n))
        enddo; enddo; enddo
        call pack_data(grid,ef_REG_IJ,ef_REG_IJ_glob)
        if(am_i_root( ))then
          title='Regions defined in rundeck for altering tracer sources'
          call openunit('EF_REG',iu,.true.)
          ef_REG_IJ_glob4(:,:)=sngl(ef_REG_IJ_glob(:,:))
          write(iu)title,ef_REG_IJ_glob4
          call closeunit(iu)
        endif
      endif

      end subroutine setup_emis_sectors_regions
