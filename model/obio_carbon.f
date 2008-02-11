#include "rundeck_opts.h"
      subroutine obio_carbon(vrbos,kmax,i,j)
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

      USE obio_forc, only: wind,atmCO2

      USE obio_incom, only : cnratio,rlamdoc,rkdoc1,rkdoc2
     .                      ,rlampoc,uMtomgm3,Pzo,awan,stdslp
     .                      ,excz,resz,remin,excp,resp

      USE obio_com, only : bn,C_tend,obio_P,P_tend,car
     .                    ,tfac,det,D_tend,tzoo,gro,pnoice,pCO2_ij
     .                    ,temp1d,saln1d,dp1d,pnoice2,rhs,alk1d

#ifdef TRACERS_GASEXCH_CO2_Natassa

      USE TRACER_COM, only : ntm,tr_mm    !tracers involved in air-sea gas exch

      USE TRACER_GASEXCH_COM, only : tracflx1d

#endif
      
      implicit none

#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"


      integer :: nt,kmax
      real  :: cchlratio,mgchltouMC,rmmzoo,docexcz,rndep,docdep
      real  :: docbac,docdet,dicresz,sumdoc,sumutk,sumres,totgro
      real  :: docexcp,dicresp,scco2,scco2arg,wssq,rkwco2
      real  :: Ts,tk,tk100,tk1002,ff,xco2,deltco2,flxmolm3,flxmolm3h

      logical vrbos

      real term



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

        term = docbac * pnoice
        rhs(k,14,14) = term
        C_tend(k,2) = C_tend(k,2) + term
     
        term = tfac(k)*remin(1)*det(k,1)/uMtomgm3 * pnoice
        rhs(k,14,10) = term
        C_tend(k,2) = C_tend(k,2) + term

cdiag   if (vrbos .and. k.eq.1)then
cdiag      if (nstep.eq.120)write(909,'(a)')
cdiag.       'nstep,k,dicresz,mgchltouMC,docbac,tfac(k),det(k,1)'
cdiag      write(909,'(2i7,5e12.4)')
cdiag.        nstep,k,dicresz,mgchltouMC,docbac,tfac(k),det(k,1)

cdiag      if (nstep.eq.120)write(910,'(a)')
cdiag.'nstep,k,tzoo,resz,P(k,ntyp),tfac(k),sst,P(k,1),car(k,1),det(k,1)'
cdiag      write(910,'(2i7,8e12.4)')
cdiag.   nstep,k,tzoo,resz,obio_P(k,ntyp),tfac(k),temp1d(k)
cdiag.  ,obio_P(k,1),car(k,1),det(k,1)
cdiag   endif

      enddo  !k=1,kmax

cdiag   if (vrbos) write(908,'(a,i7,e12.4)')'1: ', nstep,C_tend(1,2)

! Phytoplankton components related to growth
      do k = 1,kmax

        cchlratio = bn(k)*cnratio
        mgchltouMC = cchlratio/uMtomgm3
        sumdoc = 0.0
        sumutk = 0.0
        sumres = 0.0

        do nt = 1,nchl
         totgro = gro(k,nt)*obio_P(k,nt+nnut)
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

cdiag    if (vrbos .and. k.eq.1)
!        if(nstep.ge.48.and.nstep.lt.50)write(*,*)
!    .   'obio_carbon1:',
!    .   nstep,nt+nnut,i,j,k,gro(k,nt),obio_P(k,nt+nnut),totgro,
!    .   excp,docexcp,resp,dicresp,sumdoc,sumutk,sumres

        enddo !nt

        term = sumdoc * mgchltouMC * pnoice     !phyto prod DOC
        rhs(k,13,5) = term        
        C_tend(k,1) = C_tend(k,1) + term

        term = ((sumres-sumutk)*mgchltouMC) * pnoice     !phyto prod DIC
        rhs(k,14,5) = term
        C_tend(k,2) = C_tend(k,2) + term

      enddo !k=1,kmax

cdiag   if (vrbos) write(908,'(a,i7,e12.4)')'2: ', nstep,C_tend(1,2)

c pCO2
      call ppco2tab(temp1d(1),saln1d(1),car(1,2),alk1d(1),pCO2_ij)

cdiag   if (vrbos) then
cdiag   do nt = 1,nchl
cdiag     write(802,'(i7,2e12.4)')
cdiag.          nstep,gro(1,nt),obio_P(1,nt+nnut)
cdiag   enddo
cdiag   endif



!!!!!!this is all needs to be done in the gas exchange routine

!this part needs pco2_ij to be passed to atmosphere to compute the 
!flux. prolly we need to abandon this step here, compute surface tracer,
!pass to the atmosphere, compute flux and then use kpp to change the surface
!flux of the tracer? no? why not?


c Update DIC for sea-air flux of CO2

!this is for gas exchange + ocean biology: MAKE SURE!!!!
#if defined(TRACERS_GASEXCH_Natassa) && defined(TRACERS_GASEXCH_CO2_Natassa)
      !when ocean biology + CO2 gas exch
      k = 1
      do nt=1,ntm
      term = tracflx1d(nt)           ! kg/m2/s
     .     * 3600.0                  ! kg/m2/hr
     .     /dp1d(k)
!    .     / tr_mm(nt)/1.e-3/dp1d(k) ! mol/m3/hr
!    .     * 1.e6                    ! micro-mol/m3/hr
     .     * 1000.0*pnoice           !units of uM/hr (=mili-mol/lt/hr)
      rhs(k,14,16) = term
      C_tend(k,2) = C_tend(k,2) + term

         if (vrbos)
!    .   write(6,'(a,3i7,i3,4e12.4)')
     .   write(6,*)
     .     'obio_carbon (coupled):',
     .     nstep,i,j,nt,tr_mm(nt),dp1d(1),tracflx1d(nt),term

      enddo

#else

!this is for only ocean biology but no gas exchange: MAKE SURE!!!!
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
        scco2arg = (scco2/660.0)**(-0.5)      !Schmidt number
        rkwco2 = awan*wssq*scco2arg           !transfer coeff (units of m/s)
      endif
      tk = 273.15+Ts
      tk100 = tk*0.01
      tk1002 = tk100*tk100
      ff = exp(-162.8301 + 218.2968/tk100  +    !solubility
     .         90.9241*log(tk100) - 1.47696*tk1002 +
     .         saln1d(k) * (.025695 - .025225*tk100 +
     .         0.0049867*tk1002))

      xco2 = atmCO2*1013.0/stdslp
      deltco2 = (xco2-pCO2_ij)*ff*1.0245E-3
      flxmolm3 = (rkwco2*deltco2/dp1d(k))   !units of mol/m3/s
      flxmolm3h = flxmolm3*3600.0           !units of mol/m3/hr
      term = flxmolm3h*1000.0*pnoice        !units of uM/hr (=mili-mol/lt/hr)
      rhs(k,14,16) = term
      C_tend(k,2) = C_tend(k,2) + term

      if (vrbos)
     . write(6,'(a,3i7,8e12.4)')'obio_carbon(watson):',
     .   nstep,i,j,Ts,scco2arg,wssq,rkwco2,ff,xco2,pCO2_ij,term

#endif

#endif

cdiag   if (vrbos) write(908,'(a,i7,e12.4)')'3: ', nstep,C_tend(1,2)

      return
      end

c ---------------------------------------------------------------------------
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

      USE obio_incom, only : it0inc,idicinc,itainc,pco2tab
     .                      ,nt0,nsal,ndic,nta
      USE obio_dim, only: ALK_CLIM

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
