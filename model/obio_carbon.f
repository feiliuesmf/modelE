#include "rundeck_opts.h"
      subroutine obio_carbon(gro,vrbos,kmax,i,j)
c
c  Computes carbon cycling.  Uses Aumont et al (2002; JGR) for
c  semi-labile DOC (because of basic similarities in model
c  formulation, specifically, presense of grazers and particulate
c  detritus).
c  ncar(1) = DOC semi-labile (uM(C))
c  ncar(2) = DIC  (uM(C))
c

      USE MODEL_COM, only: dtsrc

      USE obio_dim
      USE obio_incom, only : cnratio,rlamdoc,rkdoc1,rkdoc2
     .                      ,rlampoc,uMtomgm3,Pzo,awan,stdslp
     .                      ,excz,resz,remin,excp,resp
      USE obio_forc, only: wind,atmCO2
      USE obio_com, only : bn,C_tend,obio_P,P_tend,car
     .                    ,tfac,det,D_tend,tzoo,pnoice,pCO2_ij
     .                    ,temp1d,saln1d,dp1d,rhs,alk1d

#ifdef OBIO_ON_GARYocean
#ifndef TRACER_GASEXCH_ocean_CO2
#ifdef TRACERS_OceanBiology
     .                    ,ao_co2flux
#endif
#endif

      USE MODEL_COM, only : nstep=>itime
      USE OCEANRES, only : kdm=>lmo
#else
      USE hycom_dim_glob, only : kdm
      USE hycom_scalars, only : nstep
#endif

#ifdef TRACERS_GASEXCH_ocean_CO2
      USE TRACER_COM, only : ntm,tr_mm    !tracers involved in air-sea gas exch
      USE TRACER_GASEXCH_COM, only : tracflx1d
#endif
      

      implicit none


      integer :: i,j,k

      integer :: nt,kmax
      real  :: cchlratio,mgchltouMC,rmmzoo,docexcz,rndep,docdep
      real  :: docbac,docdet,dicresz,sumdoc,sumutk,sumres,totgro
      real  :: docexcp,dicresp,scco2,scco2arg,wssq,rkwco2
      real  :: Ts,tk,tk100,tk1002,ff,xco2,deltco2,flxmolm3,flxmolm3h
      real  :: gro(kdm,nchl)
      real  :: pHsfc
      real term

      logical vrbos


      do nt=1,ncar
       do k=1,kdm
        C_tend(k,nt) = 0.0
       enddo
      enddo

! Detrital, bacterial, and grazing components
      do k = 1,kmax

        cchlratio = bn(k)*cnratio
        mgchltouMC = cchlratio/uMtomgm3

        !---------------------------------------------------------------
        !DOC
        rmmzoo = obio_P(k,ntyp)/(Pzo+obio_P(k,ntyp))

        docexcz = excz*rmmzoo*obio_P(k,ntyp)  !zoopl production DOC
        term = - docexcz*pnoice
        rhs(k,ntyp,14) = term
        P_tend(k,ntyp) = P_tend(k,ntyp) + term

        rndep = rlamdoc*(obio_P(k,1)/(rkdoc1 + obio_P(k,1)))
        docdep = car(k,1)/(rkdoc2+car(k,1))
        docbac = tfac(k)*rndep*docdep*car(k,1)   !bacterial loss DOC
        docdet = tfac(k)*rlampoc*det(k,1)        !detrital production DOC

        term = (docexcz*mgchltouMC
     .                   +  docdet/uMtomgm3-docbac)*pnoice
        rhs(k,13,14) = term
        C_tend(k,1) = term

        !adjust detritus
        term = - docdet !carbon/nitrogen detritus
        rhs(k,10,14) = term
        D_tend(k,1) = D_tend(k,1) + term  

        !equivalent amount of DOC from detritus goes into 
        !nitrate, bypassing DON
        term = docdet/cnratio
        rhs(k,1,14) = term
        P_tend(k,1) = P_tend(k,1) + term 

        !---------------------------------------------------------------
        !DIC
        dicresz = tzoo*resz*obio_P(k,ntyp) !zoopl production DIC (resp)
        term = - dicresz*pnoice
        rhs(k,ntyp,15) =  term
        P_tend(k,ntyp) = P_tend(k,ntyp) + term

        term = dicresz*mgchltouMC * pnoice
        rhs(k,14,15) = term
        C_tend(k,2) = term

!       if(k.eq.1)
!    .     write(*,'(a,3i5,e12.4)')'dicterm1',
!    .     nstep,i,j,term

        term = docbac * pnoice
        rhs(k,14,14) = term
        C_tend(k,2) = C_tend(k,2) + term
     
!       if(k.eq.1)
!    .     write(*,'(a,3i5,e12.4)')'dicterm2',
!    .     nstep,i,j,term

        term = tfac(k)*remin(1)*det(k,1)/uMtomgm3 * pnoice
        rhs(k,14,10) = term
        C_tend(k,2) = C_tend(k,2) + term

!       if(k.eq.1)
!    .     write(*,'(a,3i5,e12.4)')'dicterm3',
!    .     nstep,i,j,term

!     if(k.eq.1)write(*,'(a,3i5,11e12.4)')'dic_carbon:',
!    . nstep,i,j,tzoo,resz,obio_P(k,ntyp),dicresz,pnoice,
!    . mgchltouMC,docbac,tfac(k),remin(1),det(k,1),uMtomgm3

      enddo  !k=1,kmax

cdiag   if (vrbos) write(*,'(a,i7,e12.4)')
cdiag.        'obio_carbon1: ', nstep,C_tend(1,2)

! Phytoplankton components related to growth
      do k = 1,kmax

        cchlratio = bn(k)*cnratio
        mgchltouMC = cchlratio/uMtomgm3
        sumdoc = 0.0
        sumutk = 0.0
        sumres = 0.0

        do nt = 1,nchl
         !!totgro = gro(k,nt)*obio_P(k,nt+nnut)
         totgro = gro(k,nt)
          docexcp = excp*totgro   !phyto production DOC
           dicresp = resp*totgro   !phyto production DIC

            term = - docexcp * pnoice
            rhs(k,nt+nnut,14) = term
            P_tend(k,nt+nnut) = P_tend(k,nt+nnut) + term

            term = - dicresp * pnoice
            rhs(k,nt+nnut,15) = term
            P_tend(k,nt+nnut) = P_tend(k,nt+nnut) + term

           sumdoc = sumdoc + docexcp
          sumutk = sumutk + totgro
         sumres = sumres + dicresp

        enddo !nt

        term = sumdoc * mgchltouMC * pnoice     !phyto prod DOC
        rhs(k,13,5) = term        
        C_tend(k,1) = C_tend(k,1) + term

        term = ((sumres-sumutk)*mgchltouMC) * pnoice     !phyto prod DIC
        rhs(k,14,5) = term
        C_tend(k,2) = C_tend(k,2) + term

!       if(k.eq.1)
!    .     write(*,'(a,3i5,e12.4)')'dicterm4',
!    .     nstep,i,j,term


      enddo !k=1,kmax

cdiag if (vrbos) write(*,'(a,i7,e12.4)')
cdiag.    'obio_carbon2: ', nstep,C_tend(1,2)

c pCO2
#ifdef pCO2_ONLINE
      !this ppco2 routine comes from OCMIP. I am not using psurf
      !and thus not compute dtco2 because these are computed in PBL
      call ppco2(temp1d(1),saln1d(1),car(1,2),alk1d(1),
     .           obio_P(1,1),obio_P(1,3),atmCO2,
     .           pCO2_ij,pHsfc)

!note: pco2 is computed as if it is 100% open ocean cell. This is why
!in the flux computation below we need to take into account pnoice
!also in the diagnostics

!     !limits on pco2 ---more work needed
! ppco2 does not handle well the extreme salinity cases, such as
! when ice melts/forms, in river outflows.
      if (saln1d(1).ge.40. .and. pCO2_ij.lt.100.)pCO2_ij=100.
      if (saln1d(1).le.31. .and. pCO2_ij.gt.800.)pCO2_ij=800.
      if (pCO2_ij .lt. 100.) pCO2_ij=100.
      if (pCO2_ij .gt.1000.) pCO2_ij=1000.

      if(vrbos)then
        write(*,'(a,3i5,9e12.4)')
     .    'carbon: ONLINE',nstep,i,j,temp1d(1),saln1d(1),
     .               car(1,2),alk1d(1),
     .               obio_P(1,1),obio_P(1,3),pCO2_ij,
     .               pHsfc,pnoice
      endif

#else

      call ppco2tab(temp1d(1),saln1d(1),car(1,2),alk1d(1),pCO2_ij)
      if (vrbos) then
         write(*,'(a,3i5,6e12.4)')
     .    'carbon: OFFLINE',nstep,i,j,temp1d(1),saln1d(1),
     .                      car(1,2),alk1d(1),pCO2_ij,pHsfc
      endif
#endif

c Update DIC for sea-air flux of CO2

!this is for gas exchange + ocean biology
#if defined(TRACERS_GASEXCH_ocean) && defined(TRACERS_GASEXCH_ocean_CO2)
      k = 1
      do nt=1,ntm
      term = tracflx1d(nt)           ! mol/m2/s
     .     * 3600.D0                 ! mol/m2/hr
     .     /dp1d(k)                  ! mol/m3/hr
!    .     / tr_mm(nt)/1.D-3         ! mol/m3/hr
     .     * 1000.D0*pnoice          !units of uM/hr (=mili-mol/m3/hr)
      rhs(k,14,16) = term
      C_tend(k,2) = C_tend(k,2) + term

         if (vrbos) then
         write(*,'(a,3i7,i3,4e12.4)')
     .     'obio_carbon (coupled):',
     .     nstep,i,j,nt,tr_mm(nt),dp1d(1),tracflx1d(nt),term
         endif

      enddo

#else

!this is for only ocean biology but no gas exchange: 
#ifdef TRACERS_OceanBiology
      !when ocean biology but no CO2 gas exch
      !atmco2 is set to constant
      k = 1
      Ts = temp1d(k)
      scco2 = 2073.1 - 125.62*Ts + 3.6276*Ts**2 - 0.043219*Ts**3
      wssq = wind*wind
      if (scco2.lt.0.) then
        scco2arg=1.d-10
        rkwco2=1.d-10
      else
        scco2arg = (scco2/660.D0)**(-0.5)      !Schmidt number
        rkwco2 = awan*wssq*scco2arg           !transfer coeff (units of m/s)
      endif
      tk = 273.15+Ts
      tk100 = tk*0.01
      tk1002 = tk100*tk100
      ff = exp(-162.8301 + 218.2968/tk100  +    !solubility
     .         90.9241*log(tk100) - 1.47696*tk1002 +
     .         saln1d(k) * (.025695 - .025225*tk100 +
     .         0.0049867*tk1002))

      xco2 = atmCO2*1013.D0/stdslp
      deltco2 = (xco2-pCO2_ij)*ff*1.0245D-3
      flxmolm3 = (rkwco2*deltco2/dp1d(k))   !units of mol/m3/s
      flxmolm3h = flxmolm3*3600.D0          !units of mol/m3/hr
      term = flxmolm3h*1000.D0*pnoice       !units of uM/hr (=mili-mol/m^3/hr)
      rhs(k,14,16) = term
      C_tend(k,2) = C_tend(k,2) + term

#ifndef TRACER_GASEXCH_ocean_CO2
#ifdef TRACERS_OceanBiology
      ao_co2flux= rkwco2*(xco2-pCO2_ij)*ff*1.0245D-3*pnoice  ! air-sea co2 flux
     .            *3600.D0                                   ! mol/m2/hr
     .            *44.d0*24.d0*365.d0                        ! grC/m2/yr
            write(*,'(a,3i5,7e12.4)')'obio_carbon, fluxdiag:',
     .      nstep,i,j,scco2arg,wssq,rkwco2,
     .      xco2,pCO2_ij,ff*1.0245D-3,ao_co2flux
#endif
#endif

!       if(k.eq.1)
!    .     write(*,'(a,3i5,e12.4)')'dicterm5',
!    .     nstep,i,j,term

      if (vrbos) then
       write(6,'(a,3i7,9e12.4)')'obio_carbon(watson):',
     .   nstep,i,j,Ts,scco2arg,wssq,rkwco2,ff,xco2,pCO2_ij,
     .   rkwco2*(xco2-pCO2_ij)*ff*1.0245D-3,term
      endif

#endif

#endif

      return
      end

c ---------------------------------------------------------------------------
#ifndef pCO2_ONLINE
      subroutine ppco2tab(T,S,car1D,TA,pco21D)
 
c  Computes pCO2 in the surface layer and delta pCO2 with the
c  atmosphere using OCMIP protocols.  Uses pre-computed look-up
c  table to increase computational speed.
c  Var  min     max     increment
c  T0   -2      37.5    0.5
c  sal  30      39.5    0.5
c  DIC  1800    2450    2
c  TA   2000    2500    2
c

      USE obio_dim, only: ALK_CLIM
      USE obio_incom, only : it0inc,idicinc,itainc,pco2tab
     .                      ,nt0,nsal,ndic,nta

      implicit none

      integer  :: it0,isal,idic,ita

      real, parameter :: tabar=2310.0 !mean total alkalinity uE/kg; OCMIP
      real, parameter :: Sbar=34.836  !global mean annual salinity (area-wt)
      real :: T,S,car1D,pco21D,Tsfc,sal,DIC,TA

c Get pco2
      pco21D = 0.0
       Tsfc = T
        sal = S
        DIC = car1D     !uM
        if (ALK_CLIM.eq.0) TA = tabar*S/Sbar  !adjust alk for salinity
        it0 = nint((Tsfc+2.0)*2.0)/it0inc + 1
       isal = nint((sal-30.0)*2.0) + 1
       idic = (nint(DIC-1800.0))/idicinc + 1
        ita = (nint(TA-2000.0))/itainc + 1

cdiag  write(*,'(a,4f9.3,4i6)')
cdiag.         'ppco2: ',Tsfc,sal,DIC,TA,it0,isal,idic,ita

       !test within bounds
       if ( it0.gt.nt0)  then
           it0= min( it0,nt0 )
c          write(*,'(a,2i5,f9.3)')
c    .             'correction in ppco2tab, it0=',it0,nt0,Tsfc
       endif
       if (isal.gt.nsal) then
           isal= min(isal,nsal)
c          write(*,'(a,2i5,f9.3)')
c    .             'correction in ppco2tab, isal=',isal,nsal,sal
       endif
       if (idic.gt.ndic) then
           idic= min(idic,ndic)
c          write(*,'(a,2i5,f9.3)')
c    .             'correction in ppco2tab, idic=',idic,ndic,DIC
       endif
       if ( ita.gt.nta) then
            ita = min(ita,nta)
c          write(*,'(a,2i5,f9.3)')
c    .             'correction in ppco2tab, ita =',ita,nta,TA
       endif

       if (it0.lt.  1) ita = max(it0,1)
       if (isal.lt. 1) isal= max(isal,1)
       if (idic.lt. 1) idic= max(idic,1)
       if (ita.lt.  1) ita = max(ita,1)

       pco21D = pco2tab(it0,isal,idic,ita)

      return
      end
#endif

c ---------------------------------------------------------------------------
      subroutine ppco2(T,S,car1D,TA,nitr,silic,atmCO2,
     .                 pco2,pHsfc)
c
c  Computes pCO2 in the surface layer and delta pCO2 with the 
c  atmosphere using OCMIP protocols.
c

      USE obio_dim, only: ALK_CLIM
      USE obio_incom, only: pHmin,pHmax

      implicit none

      real*8, parameter :: tabar=2310.0D0    !mean total alkalinity uE/kg; OCMIP
      real*8, parameter :: stdslp=1013.25D0  !standard sea level pressure in mb
      real*8, parameter :: Sbar=34.836D0     !global mean annual salinity (area-wt)

      real*8 dtco2,atmCO2
      real*8 T,S,car1D,TA,nitr,silic,pco2
      real*8 phlo,phhi,ph,atmpres,pHsfc
      real*8 Tsfc,sal,DIC,PO4,Si,dic_in,ta_in,pt_in,sit_in,xco2_in
      real*8 co2star,pco2surf,dpco2

c
c  Constants
      dtco2 = 0.0
      pco2 = 0.0
      phlo = pHmin
      phhi = pHmax
      phlo = 1.D0   !range for pH
      phhi =16.D0   !range for pH
      ph = 8.0
c
c  Use OCMIP subroutines
       Tsfc = T
       sal = S
       DIC = car1D             !uM
       PO4 = nitr*0.1          !uM phosphate, converted from NO3/PO4
c                               ratio from Conkright et al 1994, 
c                               global using 1st 3 standard depths 
c                               (0, 10, and 20m)
        Si= silic              !uM Si
c
c  Convert to units for co2calc
       dic_in = dic*1.0E-3  !uM to mol/m3
       if (ALK_CLIM.eq.0) TA = tabar*S/Sbar  !adjust alk for salinity
       ta_in = ta*1024.5*1.0E-6  !uE/kg to E/m3
       pt_in = PO4*1.0E-3   !uM to mol/m3
       sit_in = Si*1.0E-3   !uM to mol/m3
       xco2_in = atmco2
       !!!atmpres = slp/stdslp  --not done here; do it in PBL
c
!      print*,Tsfc,sal,dic_in,ta_in,pt_in,sit_in,
!    *      phlo,phhi,ph,xco2_in

       call co2calc_SWS(Tsfc,sal,dic_in,ta_in,pt_in,sit_in,
     *      phlo,phhi,ph,xco2_in,co2star,pco2surf)
        pco2 = pco2surf
       !!!dtco2 = dco2star
       pHsfc = pH
c
!      print*, 'pco2=',pco2

      return
      end



c_ ---------------------------------------------------------------------
c_ RCS lines preceded by "c_ "
c_ --------------------------------------------------------------------
c_
c_ $Source: /home/ialeinov/GIT_transition/cvsroot_fixed/modelE/model/obio_carbon.f,v $ 
c_ $Revision: 2.26 $
c_ $Date: 2009/06/04 22:33:55 $   ;  $State: Exp $
c_ $Author: aromanou $ ;  $Locker:  $
c_
c_ ---------------------------------------------------------------------
c_ $Log: obio_carbon.f,v $
c_ Revision 2.26  2009/06/04 22:33:55  aromanou
c_
c_ correcting bug
c_ adding diagnostics for gas exchange on ocean grid (for when TRASERS_GASEXCH_ocean is not defined)
c_ atmos-only-over-ocean diagnostics
c_
c_ Revision 2.25  2009/06/01 14:14:44  aromanou
c_
c_ new diagnostic for gas exchange flux on the ocean grid (Watson's parameterization)
c_
c_ Revision 2.24  2009/05/22 09:19:45  aromanou
c_
c_ alternative definition of the compensation depth (still a little too shallow
c_ and too uniform).  might improve later
c_
c_ new alkalinity diagnostics
c_
c_ Revision 2.23  2009/05/15 07:51:34  aromanou
c_
c_ removing bug from obio_carbon tracer moments calculation
c_ pco2 diagnostics only over open ocean (*(1-oRSI))
c_ comments in albedo
c_
c_ Revision 2.22  2009/04/10 17:01:39  aromanou
c_ testing chl=0 in albedo.f and obio_ocalbedo. Rick, you do not need to restart your runs.
c_
c_ correcting bug in pco2 diagnostic
c_
c_ rearranging terms in alkalinity formulation
c_
c_ Revision 2.21  2009/03/27 09:40:16  aromanou
c_
c_ CO2 gas exchange in Russell ocean. case of constant atmosph CO2.
c_
c_ Revision 2.20  2009/03/22 19:01:07  aromanou
c_ CO2 gas exchange in Russell ocean. Needs checking, particularly wrt units.
c_ Renamed some of the fossil rundeck options, set atrac=pco2 and interpolate onto
c_ the atmospheric grid.
c_
c_ Updated rundecks. Cases with prognostic or seawifs chlorophyl,
c_ obio-radiation coupling, gas exchange on the ocean grid (ie non-interactive
c_ atmsopheric CO2 concentrations) are working
c_ properly (as best as I can tell, anyway).
c_ Interactive gas exchange needs a bit more testing.
c_
c_ Please remember to add #include "rundeck_opts.h" when you create new routines.
c_
c_ Revision 2.19  2009/01/08 14:33:26  aromanou
c_ corrected pressure and in situ temp calculation in obio_model.
c_
c_ Revision 2.18  2008/12/29 20:12:48  aromanou
c_
c_ "fishes in another ocean"
c_
c_ ocean carbon bio-geo-chemistry module (NOBM; Watson Gregg) coupled to Gary
c_ Russell ocean. Gas exchange still on ocean grid. To be changed soon...
c_
c_ Relevant rundeck is E1arobio_g.R based on Larissa's E1F40o13.R
c_
c_ Tested results so as not to cause problems to other rundecks. (Larissa, I did
c_ my best...)
c_
c_ More testing needed to
c_ -- verify there are no restart problems/multiple processor problems
c_ -- few points remain unclear till ocean geometry is settled in Gary's ocean model
c_
c_ Revision 2.17  2008/10/06 13:39:57  aromanou
c_
c_ small changes to correct deep fluxes.
c_ pH is now passed from main model.
c_
c_ Revision 2.16  2008/08/27 02:13:51  ialeinov
c_ Some fixes needed to run the Dynamic Vegetation together with OBIO
c_ (and HYCOM). Mostly related to partial cells.
c_ Aslo used some quick hacks to prevent zero divisions and uninitialized
c_ variables in obio code - Natassa, please check.
c_
c_ Revision 2.15  2008/08/26 00:00:44  aromanou
c_
c_ changed names of precompiler options (CPP's) to prepare for full carbon
c_ (ocean-land-atmosphere) coupling.
c_
c_ Ocean CPPs are:
c_ TRACERS_GASEXCH_ocean         CPP for all ocean-atmos gas exch
c_ TRACERS_GASEXCH_ocean_CO2     CPP for CO2 exchange btw ocean-atmos
c_ TRACERS_GASEXCH_ocean_CFC     CPP for CFC exchange btw ocean-atmos
c_ Please note that presently TRACERS_GASEXCH_ocean_CO2 and
c_ TRACERS_GASEXCH_ocean_CFC cannot be used together (although it should be
c_ easy to code).
c_
c_ Atmospheric concentration of CO2 (atmCO2) is set in the rundeck. This
c_ is the case where atmCO2 is held constant.
c_
c_ Updated rundeck E1arobio1.R  to reflect all these changes
c_ Also removed obsolete rundecks
c_
c_ Revision 2.14  2008/08/07 20:59:17  aromanou
c_
c_ pco2 is now computed online instead of using a look-up table
c_ to do this need to use pco2_ONLINE precompiler option.
c_
c_ E1arobio1.R is updated to show all flags that can be used with the ocean biology code.
c_
c_ Revision 1.5  2004/07/02 15:29:01  orr
c_ Added missing parenthesis at end of "kf" equation (bug identified by Anne Mouchet)
c_
c_ Revision 1.4  2004/06/16 13:18:29  orr
c_ Initial Revision
c_
c_ ---------------------------------------------------------------------
c_ 
      subroutine co2calc_SWS(t,s,dic_in,ta_in,pt_in,sit_in
     &                  ,phlo,phhi,ph,xco2_in
     &                  ,co2star,pCO2surf)
C
C-------------------------------------------------------------------------
C
C Modified from co2calc.f (RCS version 1.8, OCMIP-2) 
C - by A. Mouchet, 2004:
C
C NOW; "All" constants are given on seawater H scale (hSWS) 
C - HOWEVER, dissociation constants for S and F are on the 'free' H scale 
C            (this is necessary to work with hSWS)
C - this routine corrects for inconsistencies in the previous version.
C
C - Other changes:
C   * use ta_iter_SWS instead of ta_iter_1;
C   * hSWS replaces htotal;
C   * moved upward the calculation of bt, st and ft 
C     (needed in evaluation of kb);
C   * added various comments and references.
C
C
C SUBROUTINE CO2CALC_SWS
C
C PURPOSE
C	Calculate delta co2* from total alkalinity and total CO2 at
C temperature (t), salinity (s) and "atmpres" atmosphere total pressure. 
C
C USAGE
C       call co2calc_SWS(t,s,dic_in,ta_in,pt_in,sit_in
C    &                  ,phlo,phhi,ph,xco2_in,atmpres
C    &                  ,co2star,dco2star,pCO2surf,dpco2)
C
C INPUT
C	dic_in = total inorganic carbon (mol/m^3) 
C                where 1 T = 1 metric ton = 1000 kg
C	ta_in  = total alkalinity (eq/m^3) 
C	pt_in  = inorganic phosphate (mol/m^3) 
C	sit_in = inorganic silicate (mol/m^3) 
C	t      = temperature (degrees C)
C	s      = salinity (PSU)
C	phlo   = lower limit of pH range
C	phhi   = upper limit of pH range
C	xco2_in=atmospheric mole fraction CO2 in dry air (ppmv) 
C	atmpres= atmospheric pressure in atmospheres (1 atm==1013.25mbar)
C
C       Note: arguments dic_in, ta_in, pt_in, sit_in, and xco2_in are 
C             used to initialize variables dic, ta, pt, sit, and xco2.
C             * Variables dic, ta, pt, and sit are in the common block 
C               "species".
C             * Variable xco2 is a local variable.
C             * Variables with "_in" suffix have different units 
C               than those without.

C OUTPUT
C	co2star  = CO2*water (mol/m^3)
C	dco2star = delta CO2 (mol/m^3)
c       pco2surf = oceanic pCO2 (ppmv)
c       dpco2    = Delta pCO2, i.e, pCO2ocn - pCO2atm (ppmv)
C
C IMPORTANT: Some words about units - (JCO, 4/4/1999)
c     - Models carry tracers in mol/m^3 (on a per volume basis)
c     - Conversely, this routine, which was written by observationalists 
c       (C. Sabine and R. Key), passes input arguments in umol/kg  
c       (i.e., on a per mass basis)
c     - I have changed things slightly so that input arguments are in mol/m^3,
c     - Thus, all input concentrations (dic_in, ta_in, pt_in, and st_in) 
c       should be given in mol/m^3; output arguments "co2star" and "dco2star"  
c       are likewise in mol/m^3.

C FILES and PROGRAMS NEEDED
C	drtsafe
C	ta_iter_SWS
C
C--------------------------------------------------------------------------
C
        implicit none

        real*8 t,s,dic_in,ta_in,pt_in,sit_in,phlo,phhi,
     .         ph,xco2_in,co2star,pCO2surf
     
        real*8 st,ft,ff,x1,x2,xacc,hSWS,hSWS2,bt,scl,s15,s2,sqrtis,
     .         tk,dic,permeg,xco2,tk100,tk1002,sqrts,permil,pt,sit,ta,
     .         dlogtk

        real*8 drtsafe

        real*8 invtk,is,is2
        real*8 k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi
        common /const/k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi,ff,hSWS
        common /species/bt,st,ft,sit,pt,dic,ta
        external ta_iter_SWS
C

c       ---------------------------------------------------------------------
C       Change units from the input of mol/m^3 -> mol/kg:
c       (1 mol/m^3)  x (1 m^3/1024.5 kg)
c       where the ocean's mean surface density is 1024.5 kg/m^3
c       Note: mol/kg are actually what the body of this routine uses 
c       for calculations.  
c       ---------------------------------------------------------------------
        permil = 1.d0 / 1024.5d0
c       To convert input in mol/m^3 -> mol/kg 

!       print*,permil,pt_in,sit_in,ta_in,dic_in

        pt=pt_in*permil
        sit=sit_in*permil
        ta=ta_in*permil
        dic=dic_in*permil
      
!       print*,pt,sit,ta,dic

c       ---------------------------------------------------------------------
C       Change units from uatm to atm. That is, atm is what the body of 
c       this routine uses for calculations.
c       ---------------------------------------------------------------------
        permeg=1.e-6
c       To convert input in uatm -> atm
        xco2=xco2_in*permeg
c       ---------------------------------------------------------------------
C
C*************************************************************************
C Calculate all constants needed to convert between various measured
C carbon species. References for each equation are noted in the code. 
C Once calculated, the constants are
C stored and passed in the common block "const". The original version of this
C code was based on the code by Dickson in Version 2 of "Handbook of Methods
C for the Analysis of the Various Parameters of the Carbon Dioxide System
C in Seawater", DOE, 1994 (SOP No. 3, p25-26). 
C
C Derive simple terms used more than once
C
      tk = 273.15 + t
      tk100 = tk/100.0
      tk1002=tk100*tk100
      invtk=1.0/tk
      dlogtk=log(tk)
      is=19.924*s/(1000.-1.005*s)
      is2=is*is
      sqrtis=sqrt(is)
      s2=s*s
      sqrts=sqrt(s)
      s15=s**1.5
      scl=s/1.80655
C
C------------------------------------------------------------------------
C Calculate concentrations for borate, sulfate, and fluoride
C
C Uppstrom (1974)
      bt = 0.000232 * scl/10.811
C Morris & Riley (1966)
      st = 0.14 * scl/96.062
C Riley (1965)
      ft = 0.000067 * scl/18.9984
C
C------------------------------------------------------------------------
C f = k0(1-pH2O)*correction term for non-ideality
C
C Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
C
      ff = exp(-162.8301 + 218.2968/tk100  +
     & 90.9241*log(tk100) - 1.47696*tk1002 +
     & s * (.025695 - .025225*tk100 + 
     & 0.0049867*tk1002))
C
C K0 from Weiss 1974
C
      k0 = exp(93.4517/tk100 - 60.2409 + 23.3585 * log(tk100) +
     & s * (.023517 - 0.023656 * tk100 + 0.0047036 * tk1002))

C
C------------------------------------------------------------------------
C k1 = [H][HCO3]/[H2CO3]
C k2 = [H][CO3]/[HCO3]     on hSWS
C
C Millero p.664 (1995) using Mehrbach et al. data on SEAWATER scale 
C (Original reference: Dickson and Millero, DSR, 1987)
C
      k1=10**(-1*(3670.7*invtk - 62.008 + 9.7944*dlogtk -
     & 0.0118 * s + 0.000116*s2))
C
      k2=10**(-1*(1394.7*invtk + 4.777 - 
     & 0.0184*s + 0.000118*s2))
C
C------------------------------------------------------------------------
C k1p = [H][H2PO4]/[H3PO4] on hSWS
C
C Millero p.670 (1995)
C
      k1p = exp(-4576.752*invtk + 115.540 - 18.453 * dlogtk +
     & (-106.736*invtk + 0.69171) * sqrts +
     & (-0.65643*invtk - 0.01844) * s)
C
C------------------------------------------------------------------------
C k2p = [H][HPO4]/[H2PO4] on hSWS
C
C Millero p.670 (1995)
C
      k2p = exp(-8814.715*invtk + 172.1033 - 27.927 * dlogtk +
     & (-160.340*invtk + 1.3566) * sqrts +
     & (0.37335*invtk - 0.05778) * s)
C
C------------------------------------------------------------------------
C k3p = [H][PO4]/[HPO4] on hSWS
C
C Millero p.670 (1995)
C
      k3p = exp(-3070.75*invtk - 18.126 + 
     & (17.27039*invtk + 2.81197) *
     & sqrts + (-44.99486*invtk - 0.09984) * s)
C
C------------------------------------------------------------------------
C ksi = [H][SiO(OH)3]/[Si(OH)4] on hSWS
C
C Millero p.671 (1995) using data from Yao and Millero (1995)
C change to (mol/ kg soln)
C
      ksi = exp(-8904.2*invtk + 117.400 - 19.334 * dlogtk +
     & (-458.79*invtk + 3.5913) * sqrtis +
     & (188.74*invtk - 1.5998) * is +
     & (-12.1652*invtk + 0.07871) * is2 +
     & log(1.0-0.001005*s))
C
C------------------------------------------------------------------------
C kw = [H][OH] on hSWS
C
C Millero p.670 (1995) using composite data
C
      kw = exp(-13847.26*invtk + 148.9802 - 23.6521 * dlogtk +
     & (118.67*invtk - 5.977 + 1.0495 * dlogtk) *
     & sqrts - 0.01615 * s)
C
C------------------------------------------------------------------------
C ks = [H][SO4]/[HSO4] on free H scale
C
C Dickson (1990, J. chem. Thermodynamics 22, 113)
C change to (mol/ kg soln)
C
      ks=exp(-4276.1*invtk + 141.328 - 23.093*dlogtk +
     & (-13856*invtk + 324.57 - 47.986*dlogtk) * sqrtis +
     & (35474*invtk - 771.54 + 114.723*dlogtk) * is -
     & 2698*invtk*is**1.5 + 1776*invtk*is2 +
     & log(1.0 - 0.001005*s))
C
C------------------------------------------------------------------------
C kf = [H][F]/[HF] on free H scale
C
C Dickson and Riley (1979)
C change to (mol/ kg soln)
C
      kf=exp(1590.2*invtk - 12.641 + 1.525*sqrtis +
     & log(1.0 - 0.001005*s)) 
C
C------------------------------------------------------------------------
C kb = [H][BO2]/[HBO2] on hSWS
C
C Dickson p.673 (1990)
C change from htotal to hSWS
C
      kb=exp( (-8966.90 - 2890.53*sqrts - 77.942*s +
     & 1.728*s15 - 0.0996*s2)*invtk +
     & (148.0248 + 137.1942*sqrts + 1.62142*s) +
     & (-24.4344 - 25.085*sqrts - 0.2474*s) *
     & dlogtk + 0.053105*sqrts*tk +
     & log((1+(st/ks)+(ft/kf))/(1+(st/ks))) )
C
C*************************************************************************
C
C Calculate [H+] SWS when DIC and TA are known at T, S and 1 atm.
C The solution converges to err of xacc. The solution must be within
C the range x1 to x2.
C
C If DIC and TA are known then either a root finding or iterative method
C must be used to calculate hSWS. In this case we use the Newton-Raphson
C "safe" method taken from "Numerical Recipes" (function "rtsafe.f" with
C error trapping removed).
C
C As currently set, this procedure iterates about 12 times. The x1 and x2
C values set below will accomodate ANY oceanographic values. If an initial
C guess of the pH is known, then the number of iterations can be reduced to
C about 5 by narrowing the gap between x1 and x2. It is recommended that
C the first few time steps be run with x1 and x2 set as below. After that,
C set x1 and x2 to the previous value of the pH +/- ~0.5. The current
C setting of xacc will result in co2star accurate to 3 significant figures
C (xx.y). Making xacc bigger will result in faster convergence also, but this
C is not recommended (xacc of 10**-9 drops precision to 2 significant figures).
C
C Parentheses added around negative exponents (Keith Lindsay)
C
      x1 = 10.0**(-phhi)
      x2 = 10.0**(-phlo)
c	xacc = 1.e-10
      xacc = 1.e-14
      hSWS = drtsafe(ta_iter_SWS,x1,x2,xacc)
C
C Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2, 
C ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
C
      hSWS2=hSWS*hSWS
      co2star=dic*hSWS2/(hSWS2 + k1*hSWS + k1*k2)
      !!!co2starair=xco2*ff*atmpres  !not needed here
      !!!dco2star=co2starair-co2star
      ph=-log10(hSWS)

c
c       ---------------------------------------------------------------
cc      Add two output arguments for storing pCO2surf
cc      Should we be using K0 or ff for the solubility here?
c       ---------------------------------------------------------------
        pCO2surf = co2star / ff
        !!!dpCO2    = pCO2surf - xco2*atmpres
C
C  Convert units of output arguments
c      Note: co2star and dco2star are calculated in mol/kg within this routine 
c      Thus Convert now from mol/kg -> mol/m^3
       co2star  = co2star / permil
       !!!dco2star = dco2star / permil

c      Note: pCO2surf and dpCO2 are calculated in atm above. 
c      Thus convert now to uatm
       pCO2surf = pCO2surf / permeg
       !!!dpCO2    = dpCO2 / permeg
C
      return
      end

c_ ---------------------------------------------------------------------
c_ RCS lines preceded by "c_ "
c_ ---------------------------------------------------------------------
c_
c_ $Source: /home/ialeinov/GIT_transition/cvsroot_fixed/modelE/model/obio_carbon.f,v $ 
c_ $Revision: 2.26 $
c_ $Date: 2009/06/04 22:33:55 $   ;  $State: Exp $
c_ $Author: aromanou $ ;  $Locker:  $
c_
c_ ---------------------------------------------------------------------
c_ $Log: obio_carbon.f,v $
c_ Revision 2.26  2009/06/04 22:33:55  aromanou
c_
c_ correcting bug
c_ adding diagnostics for gas exchange on ocean grid (for when TRASERS_GASEXCH_ocean is not defined)
c_ atmos-only-over-ocean diagnostics
c_
c_ Revision 2.25  2009/06/01 14:14:44  aromanou
c_
c_ new diagnostic for gas exchange flux on the ocean grid (Watson's parameterization)
c_
c_ Revision 2.24  2009/05/22 09:19:45  aromanou
c_
c_ alternative definition of the compensation depth (still a little too shallow
c_ and too uniform).  might improve later
c_
c_ new alkalinity diagnostics
c_
c_ Revision 2.23  2009/05/15 07:51:34  aromanou
c_
c_ removing bug from obio_carbon tracer moments calculation
c_ pco2 diagnostics only over open ocean (*(1-oRSI))
c_ comments in albedo
c_
c_ Revision 2.22  2009/04/10 17:01:39  aromanou
c_ testing chl=0 in albedo.f and obio_ocalbedo. Rick, you do not need to restart your runs.
c_
c_ correcting bug in pco2 diagnostic
c_
c_ rearranging terms in alkalinity formulation
c_
c_ Revision 2.21  2009/03/27 09:40:16  aromanou
c_
c_ CO2 gas exchange in Russell ocean. case of constant atmosph CO2.
c_
c_ Revision 2.20  2009/03/22 19:01:07  aromanou
c_ CO2 gas exchange in Russell ocean. Needs checking, particularly wrt units.
c_ Renamed some of the fossil rundeck options, set atrac=pco2 and interpolate onto
c_ the atmospheric grid.
c_
c_ Updated rundecks. Cases with prognostic or seawifs chlorophyl,
c_ obio-radiation coupling, gas exchange on the ocean grid (ie non-interactive
c_ atmsopheric CO2 concentrations) are working
c_ properly (as best as I can tell, anyway).
c_ Interactive gas exchange needs a bit more testing.
c_
c_ Please remember to add #include "rundeck_opts.h" when you create new routines.
c_
c_ Revision 2.19  2009/01/08 14:33:26  aromanou
c_ corrected pressure and in situ temp calculation in obio_model.
c_
c_ Revision 2.18  2008/12/29 20:12:48  aromanou
c_
c_ "fishes in another ocean"
c_
c_ ocean carbon bio-geo-chemistry module (NOBM; Watson Gregg) coupled to Gary
c_ Russell ocean. Gas exchange still on ocean grid. To be changed soon...
c_
c_ Relevant rundeck is E1arobio_g.R based on Larissa's E1F40o13.R
c_
c_ Tested results so as not to cause problems to other rundecks. (Larissa, I did
c_ my best...)
c_
c_ More testing needed to
c_ -- verify there are no restart problems/multiple processor problems
c_ -- few points remain unclear till ocean geometry is settled in Gary's ocean model
c_
c_ Revision 2.17  2008/10/06 13:39:57  aromanou
c_
c_ small changes to correct deep fluxes.
c_ pH is now passed from main model.
c_
c_ Revision 2.16  2008/08/27 02:13:51  ialeinov
c_ Some fixes needed to run the Dynamic Vegetation together with OBIO
c_ (and HYCOM). Mostly related to partial cells.
c_ Aslo used some quick hacks to prevent zero divisions and uninitialized
c_ variables in obio code - Natassa, please check.
c_
c_ Revision 2.15  2008/08/26 00:00:44  aromanou
c_
c_ changed names of precompiler options (CPP's) to prepare for full carbon
c_ (ocean-land-atmosphere) coupling.
c_
c_ Ocean CPPs are:
c_ TRACERS_GASEXCH_ocean         CPP for all ocean-atmos gas exch
c_ TRACERS_GASEXCH_ocean_CO2     CPP for CO2 exchange btw ocean-atmos
c_ TRACERS_GASEXCH_ocean_CFC     CPP for CFC exchange btw ocean-atmos
c_ Please note that presently TRACERS_GASEXCH_ocean_CO2 and
c_ TRACERS_GASEXCH_ocean_CFC cannot be used together (although it should be
c_ easy to code).
c_
c_ Atmospheric concentration of CO2 (atmCO2) is set in the rundeck. This
c_ is the case where atmCO2 is held constant.
c_
c_ Updated rundeck E1arobio1.R  to reflect all these changes
c_ Also removed obsolete rundecks
c_
c_ Revision 2.14  2008/08/07 20:59:17  aromanou
c_
c_ pco2 is now computed online instead of using a look-up table
c_ to do this need to use pco2_ONLINE precompiler option.
c_
c_ E1arobio1.R is updated to show all flags that can be used with the ocean biology code.
c_
c_ Revision 1.4  2004/06/16 13:18:29  orr
c_ Initial Revision
c_
c_ ---------------------------------------------------------------------
c_ 
        subroutine ta_iter_SWS(x,fn,df)

        implicit none

        real*8 x,fn,df,b2,db,dic,bt,pt,sit,ta,x3,x2,c,a,a2,da,b,
     .         st,ft,ff,hSWS

        real*8 k12,k12p,k123p
        real*8 k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi
        common /const/k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi,ff,hSWS
        common /species/bt,st,ft,sit,pt,dic,ta
C
C Modified from ta_iter_1.f (RCS version 1.2, OCMIP-2)
C - by A. Mouchet, 2004:
C Fixed Problems w/ version of ta_iter_1.f used in OCMIP-2 (vers. 1.2)
C  1) fixed errors in signs, parenthesis and coefficient c in derivative
C  2) changed from Total to Seawater Scale 
C     * c defined for seawater H scale; 
C     * fn and df adapted to KF on free H scale
C     * comments have been adapted
C

C
C This routine expresses TA as a function of DIC, hSWS and constants.
C It also calculates the derivative of this function with respect to 
C hSWS. It is used in the iterative solution for hSWS. In the call
C "x" is the input value for hSWS, "fn" is the calculated value for TA
C and "df" is the value for dTA/dhSWS
C
      x2=x*x
      x3=x2*x
      k12 = k1*k2
      k12p = k1p*k2p
      k123p = k12p*k3p
      c = 1.0 + st/ks + ft/kf
      a = x3 + k1p*x2 + k12p*x + k123p
      a2=a*a
      da = 3.0*x2 + 2.0*k1p*x + k12p
      b = x2 + k1*x + k12
      b2=b*b
      db = 2.0*x + k1
C
C	fn = hco3+co3+borate+oh+hpo4+2*po4+silicate-hfree-hso4-hf-h3po4-ta
C===========================================================================
C
      fn = k1*x*dic/b +
     &       2.0*dic*k12/b +
     &       bt/(1.0 + x/kb) +
     &       kw/x +
     &       pt*k12p*x/a +
     &       2.0*pt*k123p/a +
     &       sit/(1.0 + x/ksi) -
     &       x/c -
     &       st/(1.0 + ks/(x/c)) -
     &       ft/(1.0 + kf/(x/c)) -
     &       pt*x3/a -
     &       ta
C
C	df = dfn/dx
C
       df = ((k1*dic*b) - k1*x*dic*db)/b2 -
     &       2.0*dic*k12*db/b2 -
     &       bt/kb/(1.0+x/kb)**2. -
     &       kw/x2 +
     &       (pt*k12p*(a - x*da))/a2 -
     &       2.0*pt*k123p*da/a2 -
     &       sit/ksi/(1.0+x/ksi)**2. -
     &       1.0/c -
     &       st *(1.0 + ks/(x/c))**(-2.0) *(ks*c/x2) -
     &       ft*(1.0 + kf/(x/c))**(-2.0) *(kf*c/x2) -
     &       pt*x2*(3.0*a-x*da)/a2
C
      return
      end

c_ ---------------------------------------------------------------------
c_ RCS lines preceded by "c_ "
c_ ---------------------------------------------------------------------
c_
c_ $Source: /home/ialeinov/GIT_transition/cvsroot_fixed/modelE/model/obio_carbon.f,v $ 
c_ $Revision: 2.26 $
c_ $Date: 2009/06/04 22:33:55 $   ;  $State: Exp $
c_ $Author: aromanou $ ;  $Locker:  $
c_
c_ ---------------------------------------------------------------------
c_ $Log: obio_carbon.f,v $
c_ Revision 2.26  2009/06/04 22:33:55  aromanou
c_
c_ correcting bug
c_ adding diagnostics for gas exchange on ocean grid (for when TRASERS_GASEXCH_ocean is not defined)
c_ atmos-only-over-ocean diagnostics
c_
c_ Revision 2.25  2009/06/01 14:14:44  aromanou
c_
c_ new diagnostic for gas exchange flux on the ocean grid (Watson's parameterization)
c_
c_ Revision 2.24  2009/05/22 09:19:45  aromanou
c_
c_ alternative definition of the compensation depth (still a little too shallow
c_ and too uniform).  might improve later
c_
c_ new alkalinity diagnostics
c_
c_ Revision 2.23  2009/05/15 07:51:34  aromanou
c_
c_ removing bug from obio_carbon tracer moments calculation
c_ pco2 diagnostics only over open ocean (*(1-oRSI))
c_ comments in albedo
c_
c_ Revision 2.22  2009/04/10 17:01:39  aromanou
c_ testing chl=0 in albedo.f and obio_ocalbedo. Rick, you do not need to restart your runs.
c_
c_ correcting bug in pco2 diagnostic
c_
c_ rearranging terms in alkalinity formulation
c_
c_ Revision 2.21  2009/03/27 09:40:16  aromanou
c_
c_ CO2 gas exchange in Russell ocean. case of constant atmosph CO2.
c_
c_ Revision 2.20  2009/03/22 19:01:07  aromanou
c_ CO2 gas exchange in Russell ocean. Needs checking, particularly wrt units.
c_ Renamed some of the fossil rundeck options, set atrac=pco2 and interpolate onto
c_ the atmospheric grid.
c_
c_ Updated rundecks. Cases with prognostic or seawifs chlorophyl,
c_ obio-radiation coupling, gas exchange on the ocean grid (ie non-interactive
c_ atmsopheric CO2 concentrations) are working
c_ properly (as best as I can tell, anyway).
c_ Interactive gas exchange needs a bit more testing.
c_
c_ Please remember to add #include "rundeck_opts.h" when you create new routines.
c_
c_ Revision 2.19  2009/01/08 14:33:26  aromanou
c_ corrected pressure and in situ temp calculation in obio_model.
c_
c_ Revision 2.18  2008/12/29 20:12:48  aromanou
c_
c_ "fishes in another ocean"
c_
c_ ocean carbon bio-geo-chemistry module (NOBM; Watson Gregg) coupled to Gary
c_ Russell ocean. Gas exchange still on ocean grid. To be changed soon...
c_
c_ Relevant rundeck is E1arobio_g.R based on Larissa's E1F40o13.R
c_
c_ Tested results so as not to cause problems to other rundecks. (Larissa, I did
c_ my best...)
c_
c_ More testing needed to
c_ -- verify there are no restart problems/multiple processor problems
c_ -- few points remain unclear till ocean geometry is settled in Gary's ocean model
c_
c_ Revision 2.17  2008/10/06 13:39:57  aromanou
c_
c_ small changes to correct deep fluxes.
c_ pH is now passed from main model.
c_
c_ Revision 2.16  2008/08/27 02:13:51  ialeinov
c_ Some fixes needed to run the Dynamic Vegetation together with OBIO
c_ (and HYCOM). Mostly related to partial cells.
c_ Aslo used some quick hacks to prevent zero divisions and uninitialized
c_ variables in obio code - Natassa, please check.
c_
c_ Revision 2.15  2008/08/26 00:00:44  aromanou
c_
c_ changed names of precompiler options (CPP's) to prepare for full carbon
c_ (ocean-land-atmosphere) coupling.
c_
c_ Ocean CPPs are:
c_ TRACERS_GASEXCH_ocean         CPP for all ocean-atmos gas exch
c_ TRACERS_GASEXCH_ocean_CO2     CPP for CO2 exchange btw ocean-atmos
c_ TRACERS_GASEXCH_ocean_CFC     CPP for CFC exchange btw ocean-atmos
c_ Please note that presently TRACERS_GASEXCH_ocean_CO2 and
c_ TRACERS_GASEXCH_ocean_CFC cannot be used together (although it should be
c_ easy to code).
c_
c_ Atmospheric concentration of CO2 (atmCO2) is set in the rundeck. This
c_ is the case where atmCO2 is held constant.
c_
c_ Updated rundeck E1arobio1.R  to reflect all these changes
c_ Also removed obsolete rundecks
c_
c_ Revision 2.14  2008/08/07 20:59:17  aromanou
c_
c_ pco2 is now computed online instead of using a look-up table
c_ to do this need to use pco2_ONLINE precompiler option.
c_
c_ E1arobio1.R is updated to show all flags that can be used with the ocean biology code.
c_
c_ Revision 1.1  1999/04/03 22:00:42  orr
c_ Initial revision
c_
c_ ---------------------------------------------------------------------
c_ 
      REAL*8 FUNCTION DRTSAFE(FUNCD,X1,X2,XACC)
C
C	File taken from Numerical Recipes. Modified  R.M.Key 4/94
C
      implicit none

      integer j,maxit
      real*8 xacc,dx,temp,f,df,dxold,swap,xh,xl,x1,fl,x2,fh

      MAXIT=100
      CALL FUNCD(X1,FL,DF)
      CALL FUNCD(X2,FH,DF)
      IF(FL .LT. 0.0) THEN
        XL=X1
        XH=X2
      ELSE
        XH=X1
        XL=X2
        SWAP=FL
        FL=FH
        FH=SWAP
      END IF
      DRTSAFE=.5*(X1+X2)
      DXOLD=ABS(X2-X1)
      DX=DXOLD
      CALL FUNCD(DRTSAFE,F,DF)
      DO 100, J=1,MAXIT
        IF(((DRTSAFE-XH)*DF-F)*((DRTSAFE-XL)*DF-F) .GE. 0.0 .OR.
     &     ABS(2.0*F) .GT. ABS(DXOLD*DF)) THEN
          DXOLD=DX
          DX=0.5*(XH-XL)
          DRTSAFE=XL+DX
          IF(XL .EQ. DRTSAFE)RETURN
        ELSE
          DXOLD=DX
          DX=F/DF
          TEMP=DRTSAFE
          DRTSAFE=DRTSAFE-DX
          IF(TEMP .EQ. DRTSAFE)RETURN
        END IF
        IF(ABS(DX) .LT. XACC)RETURN
        CALL FUNCD(DRTSAFE,F,DF)
        IF(F .LT. 0.0) THEN
          XL=DRTSAFE
          FL=F
        ELSE
          XH=DRTSAFE
          FH=F
        END IF
  100  CONTINUE
      RETURN
      END

