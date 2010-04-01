#include "rundeck_opts.h"

      MODULE TRACER_GASEXCH_COM

      USE TRACER_COM, only : ntm    !tracers in air-sea gas exch

      implicit none

#include "dimension2.h"

      private

      public alloc_gasexch_com
      public gather_gasexch_com_arrays
      public scatter_gasexch_com_arrays

      public
     .  tracflx !  tracer flux at air-sea intfc

      public 
     .  tracflx_glob

      public
     .  tracflx1d

      public atracflx
      public atrac

      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: atracflx
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: atrac

      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: tracflx !  tracer flux at air-sea intfc
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: tracflx_glob
  
      real*8 tracflx1d(ntm)
      common /gasexch3/tracflx1d
!$OMP THREADPRIVATE(/gasexch3/)


      contains

!------------------------------------------------------------------------------
      subroutine alloc_gasexch_com

      USE TRACER_COM, only : ntm    !tracers in air-sea gas exch

#ifdef OBIO_ON_GARYocean
      USE OCEANRES, only   : idm=>imo,jdm=>jmo
      USE MODEL_COM, only : iia=>im,jja=>jm
      USE OCEANR_DIM, only : ogrid
#else
      USE hycom_dim_glob
      USE HYCOM_DIM, only : ogrid
#endif
      integer i_0,i_1,j_0,j_1
      integer j_0h,j_1h

      i_0=ogrid%I_STRT
      i_1=ogrid%I_STOP
      j_0=ogrid%J_STRT
      j_1=ogrid%J_STOP
      j_0h=ogrid%J_STRT_HALO
      j_1h=ogrid%J_STOP_HALO

      print*,'alloc_gasexch_com:',i_0,i_1,j_0,j_1
      print*, ntm

      ALLOCATE(tracflx(i_0:i_1,j_0h:j_1h,ntm))
      ALLOCATE(tracflx_glob(idm,jdm,ntm))

      ALLOCATE(atracflx(iia,jja,ntm))
      ALLOCATE(atrac(iia,jja,ntm))

      end subroutine alloc_gasexch_com


!------------------------------------------------------------------------------
      subroutine gather_gasexch_com_arrays

#ifdef OBIO_ON_GARYocean
      USE OCEANR_DIM, only : ogrid
#else
      USE HYCOM_DIM, only : ogrid
#endif
      USE DOMAIN_DECOMP_1D, ONLY: PACK_DATA

      call pack_data( ogrid, tracflx,tracflx_glob )

      end subroutine gather_gasexch_com_arrays

!------------------------------------------------------------------------------
      subroutine scatter_gasexch_com_arrays

#ifdef OBIO_ON_GARYocean
      USE OCEANR_DIM, only : ogrid
#else
      USE HYCOM_DIM, only : ogrid
#endif
      USE DOMAIN_DECOMP_1D, ONLY: UNPACK_DATA

      call unpack_data( ogrid, tracflx_glob,tracflx )

      end subroutine scatter_gasexch_com_arrays

      END MODULE TRACER_GASEXCH_COM

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
         !ONLY FOR HYCOM
#ifndef OBIO_ON_GARYocean
      subroutine init_gasexch_co2

!this routine is called from inside OCEAN_hycom.f and only from ROOT 
!therefore arrays here have to be global


      USE HYCOM_DIM_GLOB, only : kk,iia,jja,kdm,idm,jdm
      USE HYCOM_DIM, only : ogrid
      USE hycom_atm, only : gtracer => gtracer_glob,gtracer_loc,focean
      USE HYCOM_SCALARS, only : nstep

      USE PARAM, only: get_param


      USE TRACER_COM, only : n_CO2n,ntm  !tracers involved in air-sea gas exch

      USE TRACER_GASEXCH_COM, only : atrac

      USE obio_com, only : pCO2,pCO2_glob
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT, agrid => grid,
     &     pack_data, unpack_data, pack_block, unpack_block
      implicit none
      integer nt,i,j

      call pack_data(ogrid, pCO2, pCO2_glob)
      call pack_block( agrid,GTRACER_loc,GTRACER)

      if (AM_I_ROOT()) then

      do j=1,jja
      do i=1,iia
      if (focean(i,j).gt.0.) then
c         do nt=1,ntm
            GTRACER(n_CO2n,1,i,j)=atrac(i,j,n_CO2n)
c         enddo
      endif
      enddo
      enddo

      endif ! i am root

      call unpack_data(ogrid, pCO2_glob, pCO2)
      call unpack_block( agrid,GTRACER,GTRACER_loc)

      end subroutine init_gasexch_co2
#endif

c ---------------------------------------------------------------------

      SUBROUTINE TRACERS_GASEXCH_ocean_CO2_PBL(tg1,ws,
     . alati,psurf,itr,trconstflx,byrho,Kw_gas,alpha_gas,
     . beta_gas,trsf,trcnst,ilong,jlat)
  

      USE CONSTANT, only:    rhows,mair
      USE TRACER_COM, only : trname,tr_mm,vol2mass,n_CO2n
      USE obio_incom, only : awan
#ifdef OBIO_ON_GARYocean
      USE MODEL_COM,  only : nstep=>itime
#else
      USE HYCOM_SCALARS, only : nstep
#endif
      
      implicit none

      integer :: ilong,jlat,itr
      real*8  :: tg1,ws,psurf,trconstflx,byrho,trsf,trcnst
      real*8  :: alati,Kw_gas,alpha_gas,beta_gas
      real*8  :: Sc_gas
      real*8, external :: sc_co2,sol_co2

!@var  psurf surface pressure
!@var  alati SSS at i,j
!@var  Kw_gas  gas exchange transfer velocity at i,j only over ocean
!@var  alpha_gas  solubility of gas
!@var  beta_gas  conversion term  that includes solubility
!@var  trconstflx  constant component of surface tracer flux
!      byrho   1/surface air density

!routine for PBL calculations of CO2 gas exchange

      !---------------------------------------------------------------
      !TRANSFER VELOCITY
      !---------------------------------------------------------------
      !Schmidt number for gas
       Sc_gas=sc_co2(tg1)

      !wind speed ws: magn. of surf. wind modified by buoyancy flux (m/s)
      !compute transfer velocity Kw only over ocean
      if (Sc_gas .le. 0.) then
        write(*,'(a,2i4,a,2f9.3)')
     .          'warning: Sc_gas negtv, at ',ilong,jlat,
     .          ', Sc_gas,temp_c=',Sc_gas,tg1
         Kw_gas=1.e-10
      else
         Kw_gas=(Sc_gas/660.d0)**(-0.5d0) * ws * ws * awan !units of m/s
      endif

      !---------------------------------------------------------------
      !gas SOLUBILITY
      !---------------------------------------------------------------
      !alpha --solubility of CO2 in seawater
      alpha_gas=sol_co2(tg1,alati)    !mol/m^3/picoatm
      alpha_gas = alpha_gas * 1.d-6   !mol,CO2/m3/uatm
      alpha_gas = alpha_gas * 1024.5  !Watson Gregg
      !---------------------------------------------------------------
      !psurf is in mb. 
      beta_gas = alpha_gas * psurf/1013.25      !stdslp and psurf in mb, no need to change units

      !trsf is really sfac = Kw_gas * beta_gas
      !the term 1.0d6 / vol2mass(n_CO2n) is needed to convert uatm -> kg,co2/kg,air 
      !in the denominator of alpha
      trsf = Kw_gas * beta_gas * 1.0d6 / vol2mass(n_CO2n)

      !trconstflx comes in from SURFACE.f and has units kg,co2/kg,air/m2
      !therefore trcnst needs to be myltiplied by byrho before it is sent to  PBL.f
      trcnst = Kw_gas * alpha_gas * trconstflx * byrho     
     .                * 1.0d6 / vol2mass(n_CO2n)    

        if (ilong.eq.1. .and. jlat.eq.45) then
        write(*,'(a,3i7,11e12.4)')'PBL, TRACER_GASEXCH_CO2:',
!       write(*,'(a,3i7,11e12.4)')'44444444444444444444444',
     .   nstep,ilong,jlat,tg1,(Sc_gas/660.d0)**(-0.5d0),ws*ws,
     .   Kw_gas,alpha_gas,beta_gas,trsf,trcnst,trconstflx,byrho,rhows
        endif


      RETURN
      END SUBROUTINE TRACERS_GASEXCH_ocean_CO2_PBL

c ---------------------------------------------------------------------
c ---------------------------------------------------------------------
c used with TRACERS_GASEXCH_ocean_CO2 to compute transfer velocity for CO2
c
c taken from Watsons' code
c ---------------------------------------------------------------------
      REAL*8 FUNCTION sol_co2(pt,ps)
c-------------------------------------------------------------------
c
c     CO2 Solubility in seawater
c
c     pt:       temperature (degrees Celcius)
c     ps:       salinity    (o/oo)
c     sol_co2:  in mol/m3/pptv
c               1 pptv = 1 part per trillion = 10^-12 atm = 1 picoatm
c-------------------------------------------------------------------

      REAL*8    pt,ps,ta,tk100,tk1002


      ta = (pt + 273.16d0)
      tk100 = ta*0.01d0
      tk1002 = tk100*tk100
      sol_co2 = exp(-162.8301d0 + 218.2968d0/tk100  +
     .          90.9241d0*log(tk100) - 1.47696d0*tk1002 +
     .          ps * (.025695d0 - .025225d0*tk100 +
     .          0.0049867d0*tk1002))

      END

      REAL*8 FUNCTION alpha_gas2(pt,ps)
c-------------------------------------------------------------------
c
c     CO2 Solubility in seawater
c
c     pt:       temperature (degrees Celcius)
c     ps:       salinity    (o/oo)
c     sol_co2:  in mol/m3/pptv
c               1 pptv = 1 part per trillion = 10^-12 atm = 1 picoatm
c     Scaled by 1024.5 (density of sea water) and 10^6 to get
c     alpha_gas in  mol/m3/uatm
c-------------------------------------------------------------------

      REAL*8    pt,ps,sol_CO2

      alpha_gas2=sol_CO2(pt,ps) * 1.d-6 * 1024.5  !mol,CO2/m3/uatm
      
      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      REAL*8 FUNCTION sc_co2(t)
c---------------------------------------------------
c     CO2 Schmidt number 
c     as a function of temperature. 
c
c     t: temperature (degree Celcius)
c---------------------------------------------------
      IMPLICIT NONE 
      REAL*8  a1 ( 11: 12), a2 ( 11: 12), a3 ( 11: 12), a4 ( 11: 12)
      REAL*8 t
c
      sc_co2 = 2073.1 - 125.62*t + 3.6276*t**2 - 0.043219*t**3
c
      RETURN 
      END 

