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
#ifdef Shindell_Strat_chem
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
        if(y(nO,L).lt.0.)y(nO,L)=0.0
        if(y(nO1D,L).lt.0.)y(nO1D,L)=0.0
        if(y(nO3,L).lt.1.)y(nO3,L)=1.0
        if(y(n_Ox,L).lt.1.)y(n_Ox,L)=1.0
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
      USE TRACER_COM, only : n_NOx
      USE TRCHEM_Shindell_COM, only:rr,y,yNO3,nO3,nHO2,yCH3O2,nO,nC2O3,
     &                        ta,nXO2,ss,nNO,nNO2,pNOx,nNO3,nHONO,LS1J
#ifdef Shindell_Strat_chem     
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
#ifdef Shindell_Strat_chem
      iNO2form=95
#else
      iNO2form=48
#endif
c
      do L=1,lmax
c       If dawn then set NO3 back to zero:
        IF(yNO3(I,J,L).GT.0.)yNO3(I,J,L)=0.
c       B is for NO->NO2 reactions :
        B=rr(5,L)*y(nO3,L)+rr(6,L)*y(nHO2,L)
     &   +rr(iNO2form,L)*y(nO,L)

c       Troposphere
        if(l.lt.LS1J(J))then
         B=B+rr(20,L)*yCH3O2(I,J,L)
     &   +rr(39,L)*y(nC2O3,L)+4.2E-12*exp(180./ta(L))*y(nXO2,L)
        else
c       Stratosphere
#ifdef Shindell_Strat_chem
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
c        Gqq=0.
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
        if(l.lt.LS1J(J))C=C
     &   +rr(7,L)*y(nO3,L)*0.25 !forms NO3, assume some goes to NO
c        most likely rxns: NO2+NO3->NO+NO2, J5:NO3->NO+O2, J6:NO3->NO2+O
c     &   +rr(24,L)*y(nNO3,L) !no NO3 during day
c
        p2=B/(B+C)
        p1=1-p2
        y(nNO,L)=p1*y(n_NOx,L)
        y(nNO2,L)=p2*y(n_NOx,L)
C       Set limits on NO, NO2, NOx:
        if(y(nNO,L).lt.1.) y(nNO,L) = 1.0
        if(y(nNO2,L).lt.1.)y(nNO2,L)= 1.0
        if(y(n_NOx,L).lt.1.)y(n_NOx,L)= 1.0
        pNOx(I,J,L)=y(nNO2,L)/y(n_NOx,L)
        y(nNO3,L) =1.0
        y(nHONO,L)=1.0
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
      USE TRACER_COM, only : n_CH4,n_HNO3,n_CH3OOH,n_H2O2,n_HCHO,n_CO,
     &                       n_Paraffin,n_Alkenes,n_Isoprene,n_AlkylNit
#ifdef Shindell_Strat_chem
     &                       ,n_HBr,n_HOCl,n_HCl
#endif

      USE TRCHEM_Shindell_COM, only:pHOx,rr,y,nNO2,nNO,yCH3O2,nH2O,nO3,
     &                           nO2,nM,nHO2,nOH,nH2,nAldehyde,nXO2,
     &                           nXO2N,ta,ss,nC2O3,nROR,LS1J
#ifdef Shindell_Strat_chem
     &                           ,nBrO,nClO,nOClO,nBr,nCl,SF3
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
#ifdef Shindell_Strat_chem
      integer, parameter :: iH2O2form=96,iHNO3form=97,iHONOform=100
#else
      integer, parameter :: iH2O2form=49,iHNO3form=50,iHONOform=53
#endif
C
cc    Troposphere
      do L=1,LS1J(J)-1   ! >> beginning of altitude loop <<
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
     & *y(n_Paraffin,L)*0.89+rr(34,L)*y(n_Alkenes,L)
     & +rr(30,L)*y(n_Isoprene,L)*0.15+rr(33,L)*y(n_AlkylNit,L))
     & +rr(43,L)*y(nXO2,L)+y(nXO2N,L)*
     & (rr(44,L)*rr(43,L)/(4.2E-12*exp(180./ta(L))))
C
       cqqz=(2.*ss(4,L,I,J)*y(n_H2O2,L)+ss(9,L,I,J)*y(n_HNO3,L)+
     & ss(13,L,I,J)*y(n_HCHO,L)+ss(14,L,I,J)*y(n_CH3OOH,L)+
     & (rr(20,L)*y(nNO,L)+0.66*rr(27,L)*yCH3O2(I,J,L))
     & *yCH3O2(I,J,L))
C
       cqqz=cqqz+
     & ((2.*rr(10,L)*y(nH2O,L)+rr(11,L)*y(n_CH4,L))*
     & ss(2,L,I,J)*y(nO3,L))/
     & (rr(8,L)*y(nO2,L)+rr(9,L)*y(nM,L)+
     & rr(10,L)*y(nH2O,L)+rr(11,L)*y(n_CH4,L))
     & +ss(16,L,I,J)*y(nAldehyde,L)*2.+(rr(39,L)*y(nNO,L)
     & +rr(40,L)*y(nC2O3,L)*2.)*y(nC2O3,L)
     & +(rr(42,L)*0.94+1.6E3)*y(nROR,L)+rr(35,L)*y(n_Alkenes,L)
     & *y(nO3,L)*0.65+rr(31,L)*y(n_Isoprene,L)*y(nO3,L)*0.58
c
       sqroot=sqrt(bqqz*bqqz+4.*aqqz*cqqz)
       y(nHO2,L)=(sqroot-bqqz)/(2.*aqqz)
       y(nOH,L)=pHOx(I,J,L)*y(nHO2,L)
       temp_yHOx=y(nOH,L)+y(nHO2,L)
c
c      Now partition HOx into OH and HO2:
c      CZ: OH->HO2 reactions :
       cz=rr(2,L)*y(nO3,L)+rr(13,L)*y(n_CO,L)
     & +rr(14,L)*y(n_H2O2,L)+rr(19,L)*y(nH2,L)
     & +rr(21,L)*y(n_HCHO,L)+rr(37,L)*y(n_Paraffin,L)*
     & *0.11+rr(30,L)*y(n_Isoprene,L)*0.85
     & +rr(34,L)*y(n_Alkenes,L)
C
C      DZ: HO2->OH reactions :
       dz=rr(4,L)*y(nO3,L)+rr(6,L)*y(nNO,L)
     & +rr(41,L)*0.79*y(nC2O3,L)+rr(15,L)*y(nHO2,L)
     & *(2.*ss(4,L,I,J)/(2.*ss(4,L,I,J)+y(nOH,L)*rr(14,L)))
c      Previous two lines additional OH production via rxn 41,
c      which also produces HO2, and R15 then S4/(S4+S14) fraction
C
       if(cz+dz.gt.0.)then
        y(nOH,L)=(dz/(cz+dz))*temp_yHOx
        if(y(nOH,L).gt.temp_yHOx)y(nOH,L)=temp_yHOx-1.0
c---->  warning: OH caps follow
        if(j.le.3.and.y(nOH,L).ge.3.E5) y(nOH,L)=3.E5
        if(j.ge.44.and.y(nOH,L).ge.3.E5)y(nOH,L)=3.E5
       else
        y(nOH,L)=1.0
       endif
       y(nHO2,L)=(temp_yHOx-y(nOH,L))
C      Some limits on OH, HO2:
       if(y(nOH,L).lt.1.)y(nOH,L)=1.0
       if(y(nHO2,L).lt.1.)y(nHO2,L)=1.0
       if(y(nHO2,L).gt.1.E9)y(nHO2,L)=1.E9
       pHOx(I,J,L)=y(nOH,L)/y(nHO2,L)
C
      enddo  ! >> end of altitude loop <<
c
#ifdef Shindell_Strat_chem
cc    Stratosphere
      do L=LS1J(J),lmax
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
     &  +rr(52,L)*y(n_HCl,L)*pHOx(I,J,L)+rr(53,L)*y(n_HOCl,L)
     &  *pHOx(I,J,L)+rr(56,L)*y(nOClO,L)*pHOx(I,J,L)+
     &  rr(59,L)*y(nCl,L)+
     &  +rr(62,L)*y(nClO,L)*pHOx(I,J,L)+rr(63,L)*y(nClO,L)
     &  +rr(68,L)*y(n_HBr,L)*pHOx(I,J,L)+rr(72,L)*y(nBr,L)
     &  +rr(73,L)*y(nBrO,L)+rr(81,L)*y(nBrO,L)*pHOx(I,J,L)
c     &  +pHOx(I,J,L)*rr(18,L)*y(n_HO2NO2,L)+rr(98,L)*y(nNO2,L)
c
c     Use OH production without O1D explicitly
       cqqz=2*ss(4,L,i,j)*y(n_H2O2,L)+ss(9,L,i,j)*y(n_HNO3,L)
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
       rcqqz=rr(8,L)*y(n_O2,L)+rr(9,L)*y(nM,L)+
     & rr(10,L)*y(nH2O,L)+rr(11,L)*y(n_CH4,L)
       if(rcqqz.gt.1)then
        ratio=1./rcqqz
       else
        ratio=1.
       endif
       cqqz=cqqz+ratio*        
     & ((2.*rr(10,L)*y(nH2O,L)+rr(11,L)*y(n_CH4,L))*
     & ss(2,L,I,J)*y(nO3,L))/
     & (rr(8,L)*y(nO2,L)+rr(9,L)*y(nM,L)+
     & rr(10,L)*y(nH2O,L)+rr(11,L)*y(n_CH4,L))
c     &  +ss(8,L,i,j)*y(nHONO,L)
c     &  +rr(28,L)*y(nNO3,L)*y(n_HCHO,L) !no NO3 or HONO during day
c
c      if(J.eq.jprn.and.I.eq.iprn.and.L.eq.lprn)
c     & write(*,*) 'HOxfam: a,b,c,p = ',aqqz,bqqz,cqqz,pHOx(I,J,L)
c
       sqroot=sqrt(bqqz*bqqz+4*aqqz*cqqz)
       y(nHO2,L)=(sqroot-bqqz)/(2*aqqz)
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
c     &*(2*ss(4,L,i,j)/(2*ss(4,L,i,j)+y(nOH,L)*rr(14,L)))
c      Previous two lines additional OH production via 
c      R15 then S4/(S4+S14) fraction
     &  +rr(60,L)*y(nCl,L)+rr(90,L)*y(nO,L)
c
       if(cz+dz.gt.0)then
        y(nOH,L)=(dz/(cz+dz))*temp_yHOx
        if(y(nOH,L).gt.temp_yHOx)y(nOH,L)=temp_yHOx-1.0
       else
        y(nOH,L)=1.0
       endif
       y(nHO2,L)=(temp_yHOx-y(nOH,L))

c      At low pressures, include loss of OH into H
       if(L.gt.LM-4)then
        rHprod=rr(89,L)*y(nOH,L)*y(nO,L)
        rHspecloss=y(nO3,L)*1.4E-10*exp(-470./ta(L))
        rkzero=y(nM,L)*5.7d-32*((ta(L)/300)**-1.6)
        rktot=(rkzero/(1+(rkzero/7.5d-11)))
        rHspecloss=rHspecloss+y(nO2,L)*rktot
        y(nOH,L)=y(nOH,L)-rHprod/rHspecloss
       endif
c
        if(y(nOH,L).lt.1)y(nOH,L)=1.0
        if(y(nHO2,L).lt.1)y(nHO2,L)=1.0
        if(y(nHO2,L).gt.1E9)y(nHO2,L)=1.E9
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
