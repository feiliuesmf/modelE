#include "rundeck_opts.h"
      subroutine obio_carbon(vrbos,kmax)
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

      USE obio_com, only : bn,C_tend,cnratio,obio_P,P_tend,car
     .                    ,tfac,det,D_tend,tzoo,gro,pnoice,pCO2_ij
     .                    ,temp1d,saln1d,dp1d,pnoice2

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
        P_tend(k,ntyp) = P_tend(k,ntyp) - docexcz*pnoice
        rndep = rlamdoc*(obio_P(k,1)/(rkdoc1 + obio_P(k,1)))
        docdep = car(k,1)/(rkdoc2+car(k,1))
        docbac = tfac(k)*rndep*docdep*car(k,1)   !bacterial loss DOC
        docdet = tfac(k)*rlampoc*det(k,1)        !detrital production DOC

        C_tend(k,1) = (docexcz*mgchltouMC
     .                   +  docdet/uMtomgm3-docbac)*pnoice

        !adjust detritus
        D_tend(k,1) = D_tend(k,1) - docdet  !carbon/nitrogen detritus

        !equivalent amount of DOC from detritus goes into 
        !nitrate, bypassing DON
        P_tend(k,1) = P_tend(k,1) + docdet/cnratio

        !---------------------------------------------------------------
        !DIC
        dicresz = tzoo*resz*obio_P(k,ntyp) !zoopl production DIC (resp)
        P_tend(k,ntyp) = P_tend(k,ntyp) - dicresz*pnoice
        C_tend(k,2) = (dicresz*mgchltouMC + docbac
     .                  + tfac(k)*remin(1)*det(k,1)/uMtomgm3)
     .                  * pnoice


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
            P_tend(k,nt+nnut) = P_tend(k,nt+nnut) -
     .                         (docexcp+dicresp)*pnoice
           sumdoc = sumdoc + docexcp
          sumutk = sumutk + totgro
         sumres = sumres + dicresp
cdiag        if (vrbos .and. k.eq.1)then
cdiag           if (nstep.eq.120)write(911,'(a)')
cdiag. 'nstep,k,nt,gro(k,nt),P(k,nt),resp,sumres,sumutk,mgchltouMC'
cdiag         write(911,'(3i7,6e12.4)')
cdiag.        nstep,k,nt,gro(k,nt),obio_P(k,nt+nnut),resp,
cdiag.        sumres,sumutk,mgchltouMC
cdiag        endif
        enddo !nt
        C_tend(k,1) = C_tend(k,1) 
     .              + sumdoc * mgchltouMC * pnoice
        C_tend(k,2) = C_tend(k,2) 
     .              + ((sumres-sumutk)*mgchltouMC) * pnoice

      enddo !k=1,kmax

cdiag   if (vrbos) write(908,'(a,i7,e12.4)')'2: ', nstep,C_tend(1,2)

c pCO2
      call ppco2tab(temp1d(1),saln1d(1),car(1,2),pCO2_ij)

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

#if defined(TRACERS_GASEXCH_Natassa) && defined(TRACERS_GASEXCH_CO2_Natassa)
      !when ocean biology + CO2 gas exch
      k = 1
      do nt=1,ntm
      C_tend(k,2) = C_tend(k,2)
     .            + tracflx1d(nt) * trcfrq  ! kg/m2/s
     .            * 3600.0                  ! kg/m2/hr
     .            / tr_mm(nt)/1.e-3/dp1d(k) ! mol/m3/hr
     .            * 1000.0*pnoice        !units of uM/hr

       if (vrbos)
     . write(6,'(a,i7,i3,4e12.4)')
     .  'obio_carbon:',nstep,nt,tr_mm(nt),dp1d(1),tracflx1d(nt),
     .   tracflx1d(nt)*trcfrq*3600./tr_mm(nt)/1.e-3/dp1d(1)

      enddo

#else

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
        scco2arg = (scco2/660.0)**(-0.5)
        rkwco2 = awan*wssq*scco2arg           !units of m/s
      endif
      tk = 273.15+Ts
      tk100 = tk*0.01
      tk1002 = tk100*tk100
      ff = exp(-162.8301 + 218.2968/tk100  +
     .         90.9241*log(tk100) - 1.47696*tk1002 +
     .         saln1d(k) * (.025695 - .025225*tk100 +
     .         0.0049867*tk1002))

      xco2 = atmCO2*1013.0/stdslp
      deltco2 = (xco2-pCO2_ij)*ff*1.0245E-3
      flxmolm3 = (rkwco2*deltco2/dp1d(k))   !units of mol/m3/s
      flxmolm3h = flxmolm3*3600.0                !units of mol/m3/hr

      if (vrbos)
     . write(6,'(a,i7,6e12.4)')
     .  'obio_carbon:',nstep,Ts,scco2arg,wssq,rkwco2,xco2,flxmolm3h

      C_tend(k,2) = C_tend(k,2) + flxmolm3h*1000.0*pnoice !units of uM/hr
#endif

#endif

cdiag   if (vrbos) write(908,'(a,i7,e12.4)')'3: ', nstep,C_tend(1,2)

      return
      end

c ---------------------------------------------------------------------------
      subroutine ppco2tab(T,S,car1D,pco21D)
 
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
         TA = tabar*S/Sbar  !adjust alk for salinity
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
