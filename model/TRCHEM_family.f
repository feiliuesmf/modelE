#include "rundeck_opts.h"
c     Family chemistry calculations:
c     equil values are production/loss from rxns within family only
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
      SUBROUTINE Oxinit(lmax,I,J)
!@sum Oxinit Find O,O1D and Ox initial conc assuming equilibrium with O3
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on famchem0C8_M23p)
c
C**** GLOBAL parameters and variables:
c
      USE TRACER_COM, only : n_CH4, n_Ox
      USE TRCHEM_Shindell_COM, only:ss,rr,y,nO2,nM,nH2O,nO,nO1D,nO3,pOx
C
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
C
!@var az,bz,P1 dummy variables
!@var L dummy loop variable
!@var I,J passed horizontal position indicies
!@var lmax maximum altitude for chemistry
!@var iO3form O3 formation reaction from O + O2
      real*8 az, bz, P1
      integer, intent(IN) :: lmax,I,J
      integer L,iO3form
c
#ifdef SHINDELL_STRAT_CHEM
      iO3form=94
#else
      iO3form=47
#endif
C
      do L=1,lmax
c       for concentration of O:
        az=(ss(2,L,I,J)+ss(3,L,I,J))/(rr(iO3form,L)*y(nO2,L))
c       for concentration of O(1D):
        bz=ss(2,L,I,J)/
     &  (rr(8,L)*y(nO2,L)+rr(9,L)*y(nM,L)+
     &  rr(10,L)*y(nH2O,L)+rr(11,L)*y(n_CH4,L))
        P1=1./(1.+az+bz)
        y(nO,L)=P1*az*y(n_Ox,L)
        y(nO1D,L)=P1*bz*y(n_Ox,L)
        y(nO3,L)=y(n_Ox,L)-y(nO,L)-y(nO1D,L)
        if(y(nO,L).lt.0.)y(nO,L)=0.d0
        if(y(nO1D,L).lt.0.)y(nO1D,L)=0.d0
        if(y(nO3,L).lt.1.)y(nO3,L)=1.d0
        if(y(n_Ox,L).lt.1.)y(n_Ox,L)=1.d0
        pOx(I,J,L)=y(nO3,L)/y(n_Ox,L)
      enddo
c
      return
      END SUBROUTINE OXinit
C
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
      SUBROUTINE NOxfam(lmax,I,J)
!@sum NOxfam Find NOx family (NO,NO2,NO3,HONO) partitioning assuming
!@+   equilibrium at a given concentration of NOx. Only called during
!@+   daylight, then assume NO3=HONO=1.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on famchem0C8_M23p)
c
C**** GLOBAL parameters and variables:
c
      USE DYNAMICS, only: LTROPO
      USE TRACER_COM, only : n_NOx
      USE TRCHEM_Shindell_COM, only:rr,y,yNO3,nO3,nHO2,yCH3O2,nO,nC2O3,
     &                        ta,nXO2,ss,nNO,nNO2,pNOx,nNO3,nHONO
#ifdef SHINDELL_STRAT_CHEM     
     &                        ,nClO,nOClO,nBrO
#endif
C
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var b,c,p1,p2 dummy variables
!@var L dummy loop variable
!@var I,J passed horizontal position indicies
!@var lmax maximum altitude for chemistry
!@var iNO2form NO2 formation reaction from O + NO
      real*8 b,c,p1,p2
      integer L,iNO2form
      integer, intent(IN) :: lmax,I,J
c
#ifdef SHINDELL_STRAT_CHEM
      iNO2form=95
#else
      iNO2form=48
#endif
c
      do L=1,lmax
c       If dawn then set NO3 back to zero:
        IF(yNO3(I,J,L).GT.0.)yNO3(I,J,L)=0.d0
c       B is for NO->NO2 reactions :
        B=rr(5,L)*y(nO3,L)+rr(6,L)*y(nHO2,L)
     &   +rr(iNO2form,L)*y(nO,L)

c       Troposphere
        if(l.lt.LTROPO(I,J)+1)then
         B=B+rr(20,L)*yCH3O2(I,J,L)
     &   +rr(39,L)*y(nC2O3,L)+4.2d-12*exp(180./ta(L))*y(nXO2,L)
        else
c       Stratosphere
#ifdef SHINDELL_STRAT_CHEM
c
c      calculate NO3 abundance (works fine, but negligible contribution)
c       Aqq=rr(5,L)*y(nO3,L)
c       Bqq=rr(17,L)*y(nNO,L)+rr(66,L)*
c     *    y(nCl,alt)+ss(6,L,i,j)+ss(5,L,i,j)
c       if(y(n_NOx,L).gt.0)then
c        Gqq=(rr(16,L)*y(nOH,L)*y(n_HNO3,L)+
c     *   rr(65,L)*y(n_ClONO2,L)*y(nO,L)+
c     *   rr(92,L)*y(n_N2O5,L)*y(nM,L)+ss(11,L,i,j)*
c     *   y(n_HO2NO2,L)+ss(7,L,i,j)*y(n_N2O5,L)+ss(22,L,i,j)*
c     *   y(n_ClONO2,L))/y(n_NOx,L)
c       else
c        Gqq=0.d0
c       endif
c       y(nNO3,L)=(Aqq*y(nNO2,L)+Gqq*y(n_NOx,L))/
c     *   (Bqq+rr(99,L)*y(nNO2,L))
c
         B=B+rr(64,L)*y(nClO,L)+rr(67,L)*y(nOClO,L)
     &   +rr(71,L)*y(nBrO,L)
c     &   +rr(17,L)*y(nNO3,L)
#endif
        endif
c
C       C is for NO2->NO reactions :
        C=ss(1,L,I,J)+rr(26,L)*y(nO,L)
        if(l.lt.LTROPO(I,J)+1)C=C
     &   +rr(7,L)*y(nO3,L)*0.25 !forms NO3, assume some goes to NO
c        most likely rxns: NO2+NO3->NO+NO2, J5:NO3->NO+O2, J6:NO3->NO2+O
c     &   +rr(24,L)*y(nNO3,L) !no NO3 during day
c
        p2=B/(B+C)
        p1=1-p2
        y(nNO,L)=p1*y(n_NOx,L)
        y(nNO2,L)=p2*y(n_NOx,L)
C       Set limits on NO, NO2, NOx:
        if(y(nNO,L).lt.1.) y(nNO,L) = 1.d0
        if(y(nNO2,L).lt.1.)y(nNO2,L)= 1.d0
        if(y(n_NOx,L).lt.1.)y(n_NOx,L)= 1.d0
        pNOx(I,J,L)=y(nNO2,L)/y(n_NOx,L)
        y(nNO3,L) =1.d0
        y(nHONO,L)=1.d0
      enddo

      return
      END SUBROUTINE NOxfam
C
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
      SUBROUTINE HOxfam(lmax,I,J)
!@sum HOxfam Find HOx family (OH,HO2) partitioning assuming equilibrium
!@+   concentration of HOx.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on famchem0C8_M23p)
c
C**** GLOBAL parameters and variables:
c
      USE MODEL_COM, only : LM
#ifdef SHINDELL_STRAT_CHEM
     &      ,ptop,psf,sig
#endif
      USE DYNAMICS, only: LTROPO
      USE TRACER_COM, only : n_CH4,n_HNO3,n_CH3OOH,n_H2O2,n_HCHO,n_CO,
     &                       n_Paraffin,n_Alkenes,n_Isoprene,n_AlkylNit,
     &                      rsulf1,rsulf2,rsulf4,n_SO2,n_DMS
#ifdef SHINDELL_STRAT_CHEM
     &                       ,n_HBr,n_HOCl,n_HCl
#endif

      USE TRCHEM_Shindell_COM, only:pHOx,rr,y,nNO2,nNO,yCH3O2,nH2O,nO3,
     &                           nO2,nM,nHO2,nOH,nH2,nAldehyde,nXO2,
     &                           nXO2N,ta,ss,nC2O3,nROR,
     &                           yso2,ydms
#ifdef SHINDELL_STRAT_CHEM
     &                           ,nBrO,nClO,nOClO,nBr,nCl,SF3,nO
#endif
C
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var aqqz,bqqz,cqqz,cz,dz,sqroot,temp_yHOx,rcqqz,ratio dummy vars
!@var L dummy loop variable
!@var I,J passed horizontal position indicies
!@var lmax maximum altitude for chemistry
!@var iH2O2form H2O2 formation reaction from OH + OH
!@var iHNO3form HNO3 formation reaction from OH + NO2
!@var iHONOform HONO formation reaction from NO + OH
!@var rHprod,rHspecloss,rkzero,rktot temporary var during OH->H rxns
      real*8 aqqz, bqqz, cqqz, cz, dz, sqroot, temp_yHOx,
     *rcqqz,ratio,rHprod,rHspecloss,rkzero,rktot
      integer L
      integer, intent(IN) :: lmax,I,J
c
#ifdef SHINDELL_STRAT_CHEM
      integer, parameter :: iH2O2form=96,iHNO3form=97,iHONOform=100
#else
      integer, parameter :: iH2O2form=49,iHNO3form=50,iHONOform=53
#endif
#ifdef SHINDELL_STRAT_CHEM
!@var PRES local nominal pressure for regional Ox tracers
      REAL*8, DIMENSION(LM) :: PRES
#endif
c
#ifdef SHINDELL_STRAT_CHEM
      do L=1,LM
       PRES(L)=SIG(L)*(PSF-PTOP)+PTOP
      end do
#endif
C
cc    Troposphere
      do L=1,LTROPO(I,J)   ! >> beginning of altitude loop <<
c
c      First calculate equilibrium amount of HOx
c      A: loss rxns with HOx**2, B: loss rxns linear in HOx, C: prod
c       equations are in terms of HO2 (so *pHOx when OH is reactant)
C
       aqqz=2.*(pHOx(I,J,L)*rr(1,L) + pHOx(I,J,L)*pHOx(I,J,L)*
     & (rr(3,L)+rr(iH2O2form,L)) + rr(15,L))
C
       bqqz=pHOx(I,J,L)*(rr(12,L)*y(n_CH4,L)+rr(16,L)*
     & y(n_HNO3,L)+rr(23,L)*y(n_CH3OOH,L)+rr(iHNO3form,L)
     & *y(nNO2,L)+rr(iHONOform,L)*y(nNO,L))+rr(22,L)*yCH3O2(I,J,L)
     & +pHOx(I,J,L)*(rr(38,L)*y(nAldehyde,L)+rr(37,L)
     & *y(n_Paraffin,L)*0.89d0+rr(34,L)*y(n_Alkenes,L)
     & +rr(30,L)*y(n_Isoprene,L)*0.15d0+rr(33,L)*y(n_AlkylNit,L))
     & +rr(43,L)*y(nXO2,L)+y(nXO2N,L)*
     & (rr(44,L)*rr(43,L)/(4.2d-12*exp(180./ta(L))))
C oxidation of DMS,SO2 
     & +pHOx(I,J,L)*(rsulf1(i,j,l)*ydms(i,j,l) + 
     & rsulf2(i,j,l)*ydms(i,j,l))

       cqqz=(2.d0*ss(4,L,I,J)*y(n_H2O2,L)+ss(9,L,I,J)*y(n_HNO3,L)+
     & ss(13,L,I,J)*y(n_HCHO,L)+ss(14,L,I,J)*y(n_CH3OOH,L)+
     & (rr(20,L)*y(nNO,L)+0.66d0*rr(27,L)*yCH3O2(I,J,L))
     & *yCH3O2(I,J,L))
C
       cqqz=cqqz+
     & ((2.d0*rr(10,L)*y(nH2O,L)+rr(11,L)*y(n_CH4,L))*
     & ss(2,L,I,J)*y(nO3,L))/
     & (rr(8,L)*y(nO2,L)+rr(9,L)*y(nM,L)+
     & rr(10,L)*y(nH2O,L)+rr(11,L)*y(n_CH4,L))
     & +ss(16,L,I,J)*y(nAldehyde,L)*2.d0+(rr(39,L)*y(nNO,L)
     & +rr(40,L)*y(nC2O3,L)*2.d0)*y(nC2O3,L)
     & +(rr(42,L)*0.94d0+1.6d3)*y(nROR,L)+rr(35,L)*y(n_Alkenes,L)
     & *y(nO3,L)*0.65d0+rr(31,L)*y(n_Isoprene,L)*y(nO3,L)*0.58d0
c
       sqroot=sqrt(bqqz*bqqz+4.d0*aqqz*cqqz)
       y(nHO2,L)=(sqroot-bqqz)/(2.d0*aqqz)
       y(nOH,L)=pHOx(I,J,L)*y(nHO2,L)
       temp_yHOx=y(nOH,L)+y(nHO2,L)
c
c      Now partition HOx into OH and HO2:
c      CZ: OH->HO2 reactions :
       cz=rr(2,L)*y(nO3,L)+rr(13,L)*y(n_CO,L)
     & +rr(14,L)*y(n_H2O2,L)+rr(19,L)*y(nH2,L)
     & +rr(21,L)*y(n_HCHO,L)+rr(37,L)*y(n_Paraffin,L)*
     & *0.11d0+rr(30,L)*y(n_Isoprene,L)*0.85d0
     & +rr(34,L)*y(n_Alkenes,L)
C SO2 oxidation 
     & + rsulf4(i,j,l)*yso2(i,j,l)
C
C      DZ: HO2->OH reactions :
       dz=rr(4,L)*y(nO3,L)+rr(6,L)*y(nNO,L)
     & +rr(41,L)*0.79d0*y(nC2O3,L)
       if((2.d0*ss(4,L,I,J)+y(nOH,L)*rr(14,L)).ne.0.)then
         dz=dz+rr(15,L)*y(nHO2,L)*(2.d0*ss(4,L,I,J)
     &   /(2.d0*ss(4,L,I,J)+y(nOH,L)*rr(14,L)))
       end if
c      Previous few lines additional OH production via rxn 41,
c      which also produces HO2, and R15 then S4/(S4+S14) fraction
C
       if(cz+dz.gt.0.)then
        y(nOH,L)=(dz/(cz+dz))*temp_yHOx
        if(y(nOH,L).gt.temp_yHOx)y(nOH,L)=temp_yHOx-1.d0
c---->  warning: OH caps follow
        if(j.le.3.and.y(nOH,L).ge.3.d5) y(nOH,L)=3.d5
        if(j.ge.44.and.y(nOH,L).ge.3.d5)y(nOH,L)=3.d5
       else
        y(nOH,L)=1.d0
       endif
       y(nHO2,L)=(temp_yHOx-y(nOH,L))
C      Some limits on OH, HO2:
       if(y(nOH,L).lt.1.d0)y(nOH,L)=1.d0
       if(y(nHO2,L).lt.1.d0)y(nHO2,L)=1.d0
       if(y(nHO2,L).gt.1.d9)y(nHO2,L)=1.d9
       pHOx(I,J,L)=y(nOH,L)/y(nHO2,L)
C
      enddo  ! >> end of altitude loop <<
c
#ifdef SHINDELL_STRAT_CHEM
cc    Stratosphere
      do L=LTROPO(I,J)+1,lmax
c
c      First calculate equilibrium amount of HOx
c      A: loss rxns with HOx**2, B: loss rxns linear in HOx, C: prod
c       equations are in terms of HO2 (so *pHOx when OH is reactant)
C
       aqqz=2.d0*(pHOx(I,J,L)*rr(1,L) + pHOx(I,J,L)*pHOx(I,J,L)*
     & (rr(3,L)+rr(iH2O2form,L)) + rr(15,L))
C
       bqqz=pHOx(I,J,L)*(rr(12,L)*y(n_CH4,L)+rr(16,L)*
     & y(n_HNO3,L)+rr(23,L)*y(n_CH3OOH,L)+rr(iHNO3form,L)
     & *y(nNO2,L)+rr(iHONOform,L)*y(nNO,L))+rr(22,L)*yCH3O2(I,J,L)
     &  +rr(52,L)*y(n_HCl,L)*pHOx(I,J,L)+rr(53,L)*y(n_HOCl,L)
     &  *pHOx(I,J,L)+rr(56,L)*y(nOClO,L)*pHOx(I,J,L)+
     &  rr(59,L)*y(nCl,L)
     &  +rr(62,L)*y(nClO,L)*pHOx(I,J,L)+rr(63,L)*y(nClO,L)
     &  +rr(68,L)*y(n_HBr,L)*pHOx(I,J,L)+rr(72,L)*y(nBr,L)
     &  +rr(73,L)*y(nBrO,L)+rr(81,L)*y(nBrO,L)*pHOx(I,J,L)
c     &  +pHOx(I,J,L)*rr(18,L)*y(n_HO2NO2,L)+rr(98,L)*y(nNO2,L)
c
c     Use OH production without O1D explicitly
       cqqz=2.d0*ss(4,L,i,j)*y(n_H2O2,L)+ss(9,L,i,j)*y(n_HNO3,L)
c     &	+(2*rr(10,L)*y(nH2O,L)+rr(11,L)*y(n_CH4,L))*y(nO1D,L)
c     &	((2*rr(10,L)*y(nH2O,L)+rr(11,L)*y(n_CH4,L))*
c     &  ss(2,L,i,j)*y(nO3,L))/
c     &  (rr(8,L)*y(n_O2,L)+rr(9,L)*y(nM,L)+
c     &  rr(10,L)*y(nH2O,L)+rr(11,L)*y(n_CH4,L))
     &  +rr(54,L)*y(n_HCl,L)*y(nO,L)+rr(55,L)*y(n_HOCl,L)*y(nO,L)
     &  +rr(57,L)*y(n_HOCl,L)*y(nCl,L)+rr(58,L)*y(nCl,L)*
     &  y(n_H2O2,L)+rr(79,L)*y(nBr,L)*y(n_H2O2,L)
     &  +rr(84,L)*y(n_HBr,L)*y(nO,L)
     &  +SF3(I,J,L)*pHOx(I,J,L) !water vapor photolysis in SRBs

c      production from O1D limited to O1D amount
       rcqqz=rr(8,L)*y(nO2,L)+rr(9,L)*y(nM,L)+
     & rr(10,L)*y(nH2O,L)+rr(11,L)*y(n_CH4,L)
       if(rcqqz.gt.1)then
        ratio=1.d0/rcqqz
       else
        ratio=1.d0
       endif
       cqqz=cqqz+ratio*        
     & ((2.d0*rr(10,L)*y(nH2O,L)+rr(11,L)*y(n_CH4,L))*
     & ss(2,L,I,J)*y(nO3,L))/
     & (rr(8,L)*y(nO2,L)+rr(9,L)*y(nM,L)+
     & rr(10,L)*y(nH2O,L)+rr(11,L)*y(n_CH4,L))
c     &  +ss(8,L,i,j)*y(nHONO,L)
c     &  +rr(28,L)*y(nNO3,L)*y(n_HCHO,L) !no NO3 or HONO during day
c
c      if(J.eq.jprn.and.I.eq.iprn.and.L.eq.lprn)
c     & write(*,*) 'HOxfam: a,b,c,p = ',aqqz,bqqz,cqqz,pHOx(I,J,L)
c
       sqroot=sqrt(bqqz*bqqz+4.d0*aqqz*cqqz)
       y(nHO2,L)=(sqroot-bqqz)/(2.d0*aqqz)
       y(nOH,L)=pHOx(I,J,L)*y(nHO2,L)
       temp_yHOx=y(nOH,L)+y(nHO2,L)
c
c     if(J.eq.jprn.and.I.eq.iprn.and.L.eq.lprn)
c    &write(*,*)'eq end pHOx, aqqz, bqqz, cqqz, HO2, HOx, OH = ',
c    &pHOx(I,J,L),aqqz,bqqz,cqqz,y(nHO2,L),temp_yHOx,y(nOH,L)
c
c      Now partition HOx into OH and HO2:
c      CZ: OH->HO2 reactions :
       cz=rr(2,L)*y(nO3,L)+rr(13,L)*y(n_CO,L)
     & +rr(14,L)*y(n_H2O2,L)+rr(19,L)*y(nH2,L)
     &  +rr(21,L)*y(n_HCHO,L)
     &  +rr(61,L)*y(nClO,L)+rr(80,L)*y(nBrO,L)
c
       dz=rr(4,L)*y(nO3,L)+rr(6,L)*y(nNO,L)
c     &+rr(15,L)*y(nHO2,L)
c     &*(2.d0*ss(4,L,i,j)/(2.d0*ss(4,L,i,j)+y(nOH,L)*rr(14,L)))
c      Previous two lines additional OH production via 
c      R15 then S4/(S4+S14) fraction
     &  +rr(60,L)*y(nCl,L)+rr(90,L)*y(nO,L)
c
       if(cz+dz.gt.0)then
        y(nOH,L)=(dz/(cz+dz))*temp_yHOx
        if(y(nOH,L).gt.temp_yHOx)y(nOH,L)=temp_yHOx-1.d0
       else
        y(nOH,L)=1.d0
       endif
       y(nHO2,L)=(temp_yHOx-y(nOH,L))

c      At low pressures, include loss of OH into H
       if(PRES(L).lt.2.d0)then
        rHprod=rr(89,L)*y(nOH,L)*y(nO,L)
        rHspecloss=y(nO3,L)*1.4d-10*exp(-470./ta(L))
        rkzero=y(nM,L)*5.7d-32*((ta(L)/300.d0)**(-1.6))
        rktot=(rkzero/(1+(rkzero/7.5d-11)))
        rHspecloss=rHspecloss+y(nO2,L)*rktot
        y(nOH,L)=y(nOH,L)-rHprod/rHspecloss
       endif
c
        if(y(nOH,L).lt.1)y(nOH,L)=1.d0
        if(y(nHO2,L).lt.1)y(nHO2,L)=1.d0
        if(y(nHO2,L).gt.1.d9)y(nHO2,L)=1.d9
       pHOx(I,J,L)=y(nOH,L)/y(nHO2,L)
c      if(J.eq.jprn.and.I.eq.iprn.and.L.eq.lprn)
c    & write(*,*)'part end HO2, HOx, OH, cz, dz = ',
c    & y(nHO2,L),temp_yHOx,y(nOH,L),cz,dz
c
      enddo  ! end of altitude loop
#endif
c
      return
      END SUBROUTINE HOxfam
c
C
      SUBROUTINE ClOxfam(lmax,I,J,ClOx_old)
c
#ifdef SHINDELL_STRAT_CHEM
c
!@sum ClOxfam Find ClOx family (Cl,ClO,OClO,Cl2,Cl2O2) partitioning 
!@+   assuming equilibrium within ClOx.
!@auth Drew Shindell
!@ver  1.0 (based on ds4p_famchem_M23)
c
C**** GLOBAL parameters and variables:
c
      USE DYNAMICS, only: LTROPO
      USE MODEL_COM, only :  LM
      USE TRACER_COM, only : n_ClOx,n_HOCl,n_ClONO2,n_HCl,n_H2O2,n_CH4
      USE TRCHEM_Shindell_COM, only:pClOx,rr,y,nClO,nOClO,nCl,nCl2O2,
     &    ta,ss,nO3,nHO2,nNO3,nO,nNO,nBr,nOH,nBrO,nCH3O2,nM,nCl2,nH2,
     &    SZA,dt2,pClx,pOClOx,nNO2
C
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var A,B,C,D,F,G,Q,V,X,YY,dClOx,ww,p1,p2,p3,dOClO,ratioc dummy vars
!@var ClOx_old total ClOx at start of chemical timestep
!@var rnormnum is temporary var for conservation of Cl
!@var destCl, prodCl temporary vars for simple steady state Cl calc
      REAL*8, DIMENSION(LM)      :: ClOx_old 
      real*8 A,B,C,D,F,G,Q,V,X,YY,dClOx,ww,p1,p2,p3,dOClO,ratioc,
     &rnormnum,destCl,prodCl
      integer L
      integer, intent(IN) :: lmax,I,J

      do L=LTROPO(I,J),lmax  !  LS1,lmax

c      full ClOxfam code from offline photochemistry
c       goto 10
c       if(y(n_ClOx,L).gt.0.6d-9*y(nM,L))y(n_ClOx,L)=0.6d-9*y(nM,L)
       y(nClO,L)=y(n_ClOx,L)*pClOx(I,J,L)
       y(nOClO,L)=y(n_ClOx,L)*pOClOx(I,J,L)
       y(nCl,L)=y(n_ClOx,L)*pClx(I,J,L)
       if(y(n_ClOx,L).eq.0)goto50

c      low temperature stabilizes ClO dimer, use [Cl2O2] only for
c      calculating Cl amount, otherwise ignore
       if(ta(L).le.210)then
        y(nCl2O2,L)=(rr(102,L)*y(nClO,L)*y(nClO,L))/
     * (rr(93,L)+ss(20,L,i,j))
        if(y(nCl2O2,L).gt.y(nClO,L))y(nCl2O2,L)=y(nClO,L)
        y(nClO,L)=y(nClO,L)-y(nCl2O2,L)
       endif

       A=y(nO3,L)*rr(45,L)+y(nOClO,L)*rr(47,L)+y(nHO2,L)*
     *   rr(60,L)+y(nNO3,L)*rr(66,L)
       B=y(nCl,L)*rr(47,L)+y(nO,L)*rr(50,L)+y(nNO,L)*
     *   rr(67,L)+y(nBr,L)*rr(74,L)+ss(19,L,i,j)
       C=y(nO,L)*rr(46,L)+y(nO3,L)*(rr(48,L)+rr(49,L))+
     *   y(nOH,L)*rr(61,L)+y(nNO,L)*rr(64,L)+y(nBrO,L)*
     *   (rr(75,L)+rr(76,L))+ss(17,L,i,j)+rr(85,L)*
     *   y(nCH3O2,L)
       D=y(nO3,L)*rr(49,L)+y(nBrO,L)*rr(75,L)
       F=(rr(53,L)*y(nOH,L)*y(n_HOCl,L)+rr(55,L)*y(nO,L)*
     *   y(n_HOCl,L)+rr(65,L)*y(n_ClONO2,L)*y(nO,L)+rr(93,L)*
     *   y(nCl2O2,L)*y(nM,L)*2.d0)/y(n_ClOx,L)
       G=rr(62,L)*y(nOH,L)+rr(63,L)*y(nHO2,L)+rr(77,L)*
     *   y(nBrO,L)+2.d0*rr(102,L)*y(nClO,L)+rr(103,L)*y(nNO2,L)
       Q=rr(56,L)*y(nOH,L)
       V=C-D
       X=rr(51,L)*y(nOH,L)*y(nCl2,L)+rr(54,L)*y(nO,L)*
     *   y(n_HCl,L)+ 2.d0*
     *   ss(18,L,i,j)*y(nCl2,L)+2.d0*ss(20,L,i,j)*y(nCl2O2,L)+
     *   ss(21,L,i,j)*y(n_HOCl,L)+ss(22,L,i,j)*y(n_ClONO2,L)
       X=X/y(n_ClOx,L)
       YY=rr(57,L)*y(n_HOCl,L)+rr(58,L)*y(n_H2O2,L)+rr(59,L)*
     *   y(nHO2,L)+rr(82,L)*y(n_CH4,L)+rr(83,L)*y(nH2,L)
       if((dt2*y(n_ClOx,L)).ne.0)then
        dClOx=(y(n_ClOx,L)-ClOx_old(L))/(dt2*y(n_ClOx,L))
       else
        dClOx=0.d0
       endif
       ww=A+yy+dClOx

         if(v.ne.0)then 
          p1=(-(ww/v)*(B+C+G+dClOx))+A-B
          if(p1.ne.0)then
           p1=abs(-((x/v)*(B+C+G+dClOx)+F+B)/p1)
           if((B+D-Q-dClOx).ne.0)then
            p3=abs((-p1*D+D)/(B+D-Q-dClOx))
           else 
            p3=0.d0
           endif   
           p2=1-p1-p3
          else 
           p2=1.d0
           p1=0.d0
           p3=0.d0  
          endif  
         else 
          p2=1.d0 
          p1=0.d0
          p3=0.d0  
         endif          
       
       if(SZA.gt.90.)then
        dOClO=(y(nClO,L)*D-y(nOClO,L)*(B+Q))*dt2
        y(nOClO,L)=y(nOClO,L)+dOClO
        if(y(nOClO,L).lt.0.)y(nOClO,L)=0.d0
       endif 
       
c       if(SZA.lt.93)then       
          y(nCl,L)=p1*y(n_ClOx,L)
c         faction of Cl recycled within ClOx vs other paths
c          ratioC=(y(nO3,L)*rr(45,L)+y(nOClO,L)*rr(47,L)+
c     *    y(nHO2,L)*rr(60,L))
c          ratioC=ratioC/(ratioC+y(n_HOCl,L)*rr(57,L)+
c     *    y(n_H2O2,L)*rr(58,L)+y(nHO2,L)*rr(59,L)
c     *    +y(n_CH4,L)*rr(82,L)+y(nH2,L)*rr(83,L))
c          y(nCl,L)=y(nCl,L)*ratioC
          y(nClO,L)=p2*y(n_ClOx,L)
c         else 
cc          y(nCl,L)=0.d0
c        dClO=(
cc     *   y(nCl,L)*A+
c     *   y(nOClO,L)*B+y(nClOx,L)*F-y(nClO,L)*(C+G))*dt2
c          y(nClO,L)=y(nClO,L)+dClO   
c          y(nCl,L)=y(n_ClOx,L)-(y(nClO,L)+y(nOClO,L))
c         endif   
       if(SZA.lt.90.)y(nOClO,L)=p3*y(n_ClOx,L)  
         if(y(nCl,L).lt.0)y(nCl,L)=0.d0
         if(y(nClO,L).lt.1)y(nClO,L)=0.d0
         if(y(nOClO,L).lt.1)y(nOClO,L)=0.d0
         if(y(n_ClOx,L).lt.1)y(n_ClOx,L)=1.0
           
c      write(*,*) ' L, Cl, ClO, OClO, ClOx =',L,y(nCl,L),
c     * y(nClO,L),y(nOClO,L),y(nClOx,L)

c      Normalize so that amount of ClOx doesn't change
       if((y(nCl,L)+y(nClO,L)+y(nOClO,L)).gt.0.)then
          rnormnum=y(n_ClOx,L)/(y(nCl,L)+y(nClO,L)+y(nOClO,L))
          y(nCl,L)=y(nCl,L)*rnormnum
          y(nClO,L)=y(nClO,L)*rnormnum
          y(nOClO,L)=y(nOClO,L)*rnormnum
       endif

c      if(J.eq.jprn.and.I.eq.iprn.and.L.eq.lprn)write(*,*)
c     & 'ClOx fam: Cl, ClO, NO, NOx, rxn coeff ClO+NO',
c     & y(nCl,L),y(nClO,L),y(nNO,L),y(nNOx,L),
c     & rr(64,L),ww,v,B,C,G,dClOx,A,x,F,p1
c      if(J.eq.jprn.and.I.eq.iprn.and.L.eq.lprn)write(*,*)
c     & 'A: ',y(nO3,L),rr(45,L),y(nOClO,L),rr(47,L),y(nHO2,L),
c     *   rr(60,L),y(nNO3,L),rr(66,L)

       pClOx(I,J,L)=y(nClO,L)/y(n_ClOx,L)
       pClx(I,J,L)=y(nCl,L)/y(n_ClOx,L)
       pOClOx(I,J,L)=y(nOClO,L)/y(n_ClOx,L)

c      Temporarily, no het chem, Cl2 and Cl2O2 both zero.
c       y(nCl2O2,L)=0.d0
       y(nCl2,L)=0.d0 !het rxns set to produce Cl instead of Cl2

       goto 50 !skip alternative method

 10    continue
cc     simple Cl/ClO calculation
       if(y(n_ClOx,L).eq.0)goto50
       destCl=y(nO3,L)*rr(45,L)+y(nHO2,L)*rr(60,L)
     *    +y(n_HOCl,L)*rr(57,L)+
     *    y(n_H2O2,L)*rr(58,L)+y(nHO2,L)*rr(59,L)
     *    +y(n_CH4,L)*rr(82,L)+y(nH2,L)*rr(83,L)
       prodCl=(y(nO,L)*rr(46,L)+y(nO3,L)*rr(48,L)+
     *   y(nOH,L)*rr(61,L)+y(nNO,L)*rr(64,L)+ss(17,L,i,j))*
     *   y(n_ClOx,L)

c      if(J.eq.jprn.and.I.eq.iprn.and.L.eq.lprn)write(*,*)
c     * 'p diags: ',y(nO,L)*rr(46,L),y(nO3,L)*rr(48,L),
c     *   y(nOH,L)*rr(61,L),y(nNO,L)*rr(64,L),ss(17,L,i,j)
c      if(J.eq.jprn.and.I.eq.iprn.and.L.eq.lprn)write(*,*)
c     * 'l diags: ',y(nO3,L)*rr(45,L),y(nHO2,L)*rr(60,L)
c     *    ,y(n_HOCl,L)*rr(57,L),
c     *    y(n_H2O2,L)*rr(58,L),y(nHO2,L)*rr(59,L)
c     *    ,y(n_CH4,L)*rr(82,L),y(nH2,L)*rr(83,L)

       y(nCl,L)=prodCl/destCl
       y(nClO,L)=y(n_ClOx,L)-y(nCl,L)

c       p2=1.d0/((C/A)+1.d0)
c       p1=1.d0-p2
c       y(nCl,L)=p1*y(n_ClOx,L)
c       y(nClO,L)=p2*y(n_ClOx,L)

c      Normalize so that amount of ClOx doesn't change
       if((y(nCl,L)+y(nClO,L)).gt.0.)then
          rnormnum=y(n_ClOx,L)/(y(nCl,L)+y(nClO,L))
          y(nCl,L)=y(nCl,L)*rnormnum
          y(nClO,L)=y(nClO,L)*rnormnum
       endif

       pClOx(I,J,L)=y(nClO,L)/y(n_ClOx,L)
       pClx(I,J,L)=y(nCl,L)/y(n_ClOx,L)
       pOClOx(I,J,L)=0.d0

c      Temporarily, no het chem, Cl2 and Cl2O2 both zero.
       y(nCl2O2,L)=0.d0
       y(nCl2,L)=0.d0
       y(nOClO,L)=0.d0

 50   enddo  ! end of altitude loop
#endif
      return
      END SUBROUTINE ClOxfam
C
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
      SUBROUTINE BrOxfam(lmax,I,J)
c
#ifdef SHINDELL_STRAT_CHEM
c
!@sum BrOxfam Find BrOx family (Br,BrO) partitioning 
!@+   assuming equilibrium within BrOx.
!@auth Drew Shindell
!@ver  1.0 (based on ds4p_famchem_M23)
c
C**** GLOBAL parameters and variables:
c
      USE DYNAMICS, only: LTROPO
      USE TRACER_COM, only : n_BrOx,n_H2O2,n_HBr,n_HOBr,n_BrONO2
      USE TRCHEM_Shindell_COM, only:rr,y,nO3,nClO,nOClO,nNO,nO,nBr,nOH,
     &    nBrO,ss,nHO2,nNO2,pBrOx
C
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var a,b,c,d,eq,f,p2,p1 dummy vars
      real*8 a,b,c,d,eq,f,p2,p1
      integer L
      integer, intent(IN) :: lmax,I,J

      do L=LTROPO(I,J),lmax  ! LS1,lmax
         a=y(nO3,L)*rr(70,L)+y(nOClO,L)*rr(74,L)
         b=y(nO,L)*rr(69,L)+y(nNO,L)*rr(71,L)+y(nClO,L)*
     *   (rr(75,L)+rr(76,L))+2*y(nBrO,L)*rr(78,L)+y(nOH,L)*
     *   rr(80,L)+ss(25,L,i,j)
       c=rr(73,L)*y(nHO2,L)+rr(77,L)*y(nClO,L)+rr(81,L)*
     *   y(nOH,L)+rr(104,L)*y(nNO2,L)    
       d=rr(72,L)*y(nHO2,L)+rr(79,L)*y(n_H2O2,L)
       eq=rr(68,L)*y(n_HBr,L)*y(nOH,L)+rr(84,L)*y(n_HBr,L)*
     *  y(nO,L)+ss(24,L,i,j)*y(n_HOBr,L)
       f=ss(23,L,i,j)*y(n_BrONO2,L)
       if(a+b.ne.0)then
        p2=a/(a+b)
       else
        p2=0.d0
       endif
       if(p2.lt.0)p2=0.d0
       if(p2.gt.1)p2=1.d0
         p1=1-p2    
         y(nBr,L)=p1*y(n_BrOx,L)
         y(nBrO,L)=p2*y(n_BrOx,L)
         if(y(nBr,L).lt.1)y(nBr,L)=0.d0
         if(y(nBrO,L).lt.1)y(nBrO,L)=0.d0
         if(y(nBr,L).gt.1d9)y(nBr,L)=0.d0
         if(y(nBrO,L).gt.1d9)y(nBrO,L)=0.d0
         if(y(n_BrOx,L).lt.1)y(n_BrOx,L)=1.d0
c       write(*,*) ' L, Br, BrO, BrOx =',L,y(nBr,L),y(nBrO,L),
c     * y(n_BrOx,L)

       pBrOx(I,J,L)=y(nBrO,L)/y(n_BrOx,L)
      enddo ! end of altitude loop
#endif
c
      return
      END SUBROUTINE BrOxfam


