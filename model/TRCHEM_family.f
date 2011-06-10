#include "rundeck_opts.h"
c Family chemistry calculations: Equilibrium values are production/loss
c from reactions *within* family only:


      SUBROUTINE Oxinit(lmax,I,J)
!@sum Oxinit Find O,O1D & Ox initial conc assuming equilibrium with O3
!@auth Drew Shindell (modelEifications by Greg Faluvegi)

C**** GLOBAL parameters and variables:
      USE RESOLUTION, only : ls1,ptop,psf
      USE RESOLUTION, only : LM
      USE DYNAMICS, only : sig
      USE TRACER_COM, only : n_CH4, n_Ox
      USE TRCHEM_Shindell_COM, only:ss,rr,y,nO2,nM,nH2O,nO,nO1D,nO3,pOx

      IMPLICIT NONE

C**** Local parameters and variables and arguments:

!@var az,bz,P1 dummy variables
!@var L dummy loop variable
!@var I,J passed horizontal position indicies
!@var lmax maximum altitude for chemistry
!@var iO3form O3 formation reaction from O + O2
!@var PRES local nominal pressure
      integer, intent(IN)   :: lmax,I,J
      integer               :: L,iO3form
      REAL*8, DIMENSION(LM) :: PRES
      real*8                :: az, bz, P1

      PRES(1:LM)=SIG(1:LM)*(PSF-PTOP)+PTOP
#ifdef TRACERS_TERP
      iO3form=98
#else
      iO3form=95
#endif  /* TRACERS_TERP */

      do L=1,lmax
c       for concentration of O:
        az=(ss(2,L,I,J)+ss(3,L,I,J))/(rr(iO3form,L)*y(nO2,L))
c       for concentration of O(1D):
        bz=ss(2,L,I,J)/(rr(8,L)*y(nO2,L)+rr(9,L)*y(nM,L)+
     &  rr(10,L)*y(nH2O,L)+rr(11,L)*y(n_CH4,L))
        if(PRES(L) < 50.) then
          bz=bz*2.5d0
!test   else if(PRES(L) > 100.) then
!test     bz=bz*0.9d0
        endif
        P1=1.d0/(1.d0+az+bz)
        y(nO,L)=P1*az*y(n_Ox,L)
        y(nO1D,L)=P1*bz*y(n_Ox,L)
        y(nO3,L)=y(n_Ox,L)-y(nO,L)-y(nO1D,L)
        if(y(nO,L) < 0.)  y(nO,L)  =0.d0
        if(y(nO1D,L) < 0.)y(nO1D,L)=0.d0
        if(y(nO3,L) < 1.) y(nO3,L) =1.d0
        if(y(n_Ox,L) < 1.)y(n_Ox,L)=1.d0
        pOx(I,J,L)=y(nO3,L)/y(n_Ox,L)
      enddo
c
      return
      END SUBROUTINE OXinit



      SUBROUTINE NOxfam(lmax,I,J)
!@sum NOxfam Find NOx family (NO,NO2,NO3,HONO) partitioning assuming
!@+   equilibrium at a given concentration of NOx. Only called during
!@+   daylight, then assume NO3=HONO=1.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)

C**** GLOBAL parameters and variables:

      USE RESOLUTION, only         : LS1
      USE ATM_COM, only            : LTROPO
      USE TRACER_COM, only         : n_NOx
      USE TRCHEM_Shindell_COM, only:rr,y,yNO3,nO3,nHO2,yCH3O2,nO,nC2O3,
     &                  ta,nXO2,ss,nNO,nNO2,pNOx,nNO3,nHONO,which_trop
     &                  ,nClO,nOClO,nBrO

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var b,c,p1,p2 dummy variables
!@var L dummy loop variable
!@var I,J passed horizontal position indicies
!@var lmax maximum altitude for chemistry
!@var iNO2form NO2 formation reaction from O + NO
!@var maxl LTROPO(I,J) or LS1-1, depending upon what_trop variable
      integer             :: L,iNO2form,maxl
      integer, intent(IN) :: lmax,I,J
      real*8              :: b,c,p1,p2
      
#ifdef TRACERS_TERP
      iNO2form=99
#else
      iNO2form=96
#endif  /* TRACERS_TERP */

      select case(which_trop)
      case(0); maxl=ltropo(I,J)
      case(1); maxl=ls1-1
      case default; call stop_model('which_trop problem 8',255)
      end select

      do L=1,lmax
c       If dawn then set NO3 back to zero:
        IF(yNO3(I,J,L) > 0.) yNO3(I,J,L)=0.d0
c       B is for NO->NO2 reactions :
        B=rr(5,L)*y(nO3,L)+rr(6,L)*y(nHO2,L)+rr(iNO2form,L)*y(nO,L)

        if(l <= maxl)then  ! Troposphere:
          B=B+rr(20,L)*yCH3O2(I,J,L)
     &    +rr(39,L)*y(nC2O3,L)+4.2d-12*exp(180./ta(L))*y(nXO2,L)
        else               ! Stratosphere:
          B=B+rr(64,L)*y(nClO,L)+rr(67,L)*y(nOClO,L)
     &    +rr(71,L)*y(nBrO,L)
        endif

C       C is for NO2->NO reactions :
        C=ss(1,L,I,J)+rr(26,L)*y(nO,L)
        !forms NO3, assume some goes to NO:
        if(l <= maxl) C = C + rr(7,L)*y(nO3,L)*0.25d0 
        p2=B/(B+C)
        p1=1-p2
        y(nNO,L)= p1*y(n_NOx,L)
        y(nNO2,L)=p2*y(n_NOx,L)
C       Set limits on NO, NO2, NOx:
        if(y(nNO,L)   < 1.)   y(nNO,L) = 1.d0
        if(y(nNO2,L)  < 1.)  y(nNO2,L) = 1.d0
        if(y(n_NOx,L) < 1.) y(n_NOx,L) = 1.d0
        pNOx(I,J,L)=y(nNO2,L)/y(n_NOx,L)
        y(nNO3,L) =1.d0
        y(nHONO,L)=1.d0
      enddo

      return
      END SUBROUTINE NOxfam



      SUBROUTINE HOxfam(lmax,I,J)
!@sum HOxfam Find HOx family (OH,HO2) partitioning assuming equilibrium
!@+   concentration of HOx.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)

C**** GLOBAL parameters and variables:
      USE RESOLUTION, only : LS1,ptop,psf
      USE RESOLUTION, only : LM
      USE DYNAMICS, only : sig
      USE GEOM, only : LAT2D_DG
      USE ATM_COM, only: LTROPO
      USE TRACER_COM, only : n_CH4,n_HNO3,n_CH3OOH,n_H2O2,n_HCHO,n_CO,
     &                       n_Paraffin,n_Alkenes,n_Isoprene,n_AlkylNit,
#ifdef TRACERS_TERP
     &                       n_Terpenes,
#endif  /* TRACERS_TERP */
     &                       rsulf1,rsulf2,rsulf4,n_SO2,n_DMS
     &                       ,n_HBr,n_HOCl,n_HCl
      USE TRCHEM_Shindell_COM, only:pHOx,rr,y,nNO2,nNO,yCH3O2,nH2O,nO3,
     &                        nO2,nM,nHO2,nOH,nH2,nAldehyde,nXO2,nXO2N,
     &                        ta,ss,nC2O3,nROR,yso2,ydms,which_trop
     &                        ,nBrO,nClO,nOClO,nBr,nCl,SF3,nO,nCH3O2

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var aqqz,bqqz,cqqz,cz,dz,sqroot,temp_yHOx,rcqqz,ratio dummy vars
!@var L dummy loop variable
!@var I,J passed horizontal position indicies
!@var lmax maximum altitude for chemistry
!@var iH2O2form H2O2 formation reaction from OH + OH
!@var iHNO3form HNO3 formation reaction from OH + NO2
!@var iHONOform HONO formation reaction from NO + OH
!@var rHprod,rHspecloss,rkzero,rktot temporary var during OH->H rxns
!@var maxl either LTROPO(I,J) or LS1-1 depending on which_trop dbparam
!@var PRES local nominal pressure for regional Ox tracers

#ifdef TRACERS_TERP
      integer, parameter :: iH2O2form=100,iHNO3form=101,iHONOform=104
     &                     ,iTerpenesOH=92,iTerpenesO3=93
#else
      integer, parameter :: iH2O2form=97,iHNO3form=98,iHONOform=101
#endif  /* TRACERS_TERP */
      integer             :: L, maxl 
      integer, intent(IN) :: lmax,I,J
      real*8              :: aqqz, bqqz, cqqz, cz, dz, sqroot, 
     &   temp_yHOx,rcqqz,ratio,rHprod,rHspecloss,rkzero,rktot
      REAL*8, DIMENSION(LM) :: PRES

      PRES(1:LM)=SIG(1:LM)*(PSF-PTOP)+PTOP

      select case(which_trop)
      case(0); maxl=ltropo(I,J)
      case(1); maxl=ls1-1
      case default; call stop_model('which_trop problem 9',255)
      end select

C ------------ Troposphere ------------
      do L=1,maxl   ! >> beginning of altitude loop <<
c First calculate equilibrium amount of HOx:
c A: loss rxns with HOx**2
c B: loss rxns linear in HOx
c C: prod equations in terms of HO2 (so *pHOx when OH is reactant)

       aqqz=2.d0*(pHOx(I,J,L)*rr(1,L) + pHOx(I,J,L)*pHOx(I,J,L)*
     & (rr(3,L)+rr(iH2O2form,L)) + rr(15,L))

       bqqz=pHOx(I,J,L)*(rr(12,L)*y(n_CH4,L)+rr(16,L)*
     & y(n_HNO3,L)+rr(23,L)*y(n_CH3OOH,L)+rr(iHNO3form,L)
     & *y(nNO2,L)+rr(iHONOform,L)*y(nNO,L))+rr(22,L)*yCH3O2(I,J,L)
     & +pHOx(I,J,L)*(rr(38,L)*y(nAldehyde,L)+rr(37,L)
     & *y(n_Paraffin,L)*0.89d0+rr(34,L)*y(n_Alkenes,L)
     & +rr(30,L)*y(n_Isoprene,L)*0.15d0+rr(33,L)*y(n_AlkylNit,L)
#ifdef TRACERS_TERP
     & +rr(iTerpenesOH,L)*y(n_Terpenes,L)*0.15d0
#endif  /* TRACERS_TERP */
     & )
     & +rr(43,L)*y(nXO2,L)+y(nXO2N,L)*
     & (rr(44,L)*rr(43,L)/(4.2d-12*exp(180./ta(L))))
       ! oxidation of DMS,SO2: 
     & +pHOx(I,J,L)*(rsulf1(i,j,l)*ydms(i,j,l) + 
     & rsulf2(i,j,l)*ydms(i,j,l))

       cqqz=(2.d0*ss(4,L,I,J)*y(n_H2O2,L)+ss(9,L,I,J)*y(n_HNO3,L)+
     & ss(13,L,I,J)*y(n_HCHO,L)+ss(14,L,I,J)*y(n_CH3OOH,L)+
     & (rr(20,L)*y(nNO,L)+0.66d0*rr(27,L)*yCH3O2(I,J,L))
     & *yCH3O2(I,J,L))

       cqqz=cqqz+
     & ((2.d0*rr(10,L)*y(nH2O,L)+rr(11,L)*y(n_CH4,L))*
     & ss(2,L,I,J)*y(nO3,L))/
     & (rr(8,L)*y(nO2,L)+rr(9,L)*y(nM,L)+
     & rr(10,L)*y(nH2O,L)+rr(11,L)*y(n_CH4,L))
     & +ss(16,L,I,J)*y(nAldehyde,L)*2.d0+(rr(39,L)*y(nNO,L)
     & +rr(40,L)*y(nC2O3,L)*2.d0)*y(nC2O3,L)
     & +(rr(42,L)*0.94d0+1.6d3)*y(nROR,L)+rr(35,L)*y(n_Alkenes,L)
     & *y(nO3,L)*0.65d0+rr(31,L)*y(n_Isoprene,L)*y(nO3,L)*0.58d0
#ifdef TRACERS_TERP
     & +rr(iTerpenesO3,L)*y(n_Terpenes,L)*y(nO3,L)*0.58d0
#endif  /* TRACERS_TERP */

       sqroot=sqrt(bqqz*bqqz+4.d0*aqqz*cqqz)
       y(nHO2,L)=(sqroot-bqqz)/(2.d0*aqqz)
       y(nOH,L)=pHOx(I,J,L)*y(nHO2,L)
       temp_yHOx=y(nOH,L)+y(nHO2,L)

c Now partition HOx into OH and HO2:
       ! CZ: OH->HO2 reactions :
       cz=rr(2,L)*y(nO3,L)+rr(13,L)*y(n_CO,L)
     & +rr(14,L)*y(n_H2O2,L)+rr(19,L)*y(nH2,L)
     & +rr(21,L)*y(n_HCHO,L)+rr(37,L)*y(n_Paraffin,L)
     & *0.11d0+rr(30,L)*y(n_Isoprene,L)*0.85d0
#ifdef TRACERS_TERP
     & +rr(iTerpenesOH,L)*y(n_Terpenes,L)*0.85d0
#endif  /* TRACERS_TERP */
     & +rr(34,L)*y(n_Alkenes,L)
       ! SO2 oxidation: 
     & + rsulf4(i,j,l)*yso2(i,j,l)

       ! DZ: HO2->OH reactions :
       dz=rr(4,L)*y(nO3,L)+rr(6,L)*y(nNO,L)
     & +rr(41,L)*0.79d0*y(nC2O3,L)

c Previous few lines represent additional OH production via reaction 41
c which also produces HO2 and R15 then S4/(S4+S14) fraction.

       if(cz+dz > 0.)then
         y(nOH,L)=(dz/(cz+dz))*temp_yHOx
         if(y(nOH,L) > temp_yHOx) y(nOH,L)=temp_yHOx-1.d0
c---->   warning: OH caps follow   <----
!4x5hard if(j <= 3 .and. y(nOH,L) >= 3.d5) y(nOH,L)=3.d5
!4x5hard if(j >= 44 .and. y(nOH,L) >= 3.d5)y(nOH,L)=3.d5
         if(lat2d_dg(i,j) <= -80. .or. lat2d_dg(i,j) >= 80.)then
           y(nOH,L)=min(y(nOH,L),3.d5)
         endif
       else
         y(nOH,L)=1.d0
       endif
       y(nHO2,L)=(temp_yHOx-y(nOH,L))
       ! Some limits on OH, HO2:
       if(y(nOH,L) < 1.d0)y(nOH,L)=1.d0
       if(y(nHO2,L) < 1.d0)y(nHO2,L)=1.d0
       if(y(nHO2,L) > 1.d9)y(nHO2,L)=1.d9
       pHOx(I,J,L)=y(nOH,L)/y(nHO2,L)

      enddo  ! >> end of troposphere loop <<

C ------------ Stratosphere ------------
      do L=maxl+1,lmax
c First calculate equilibrium amount of HOx:
c A: loss rxns with HOx**2,
c B: loss rxns linear in HOx,
c C: prod equations in terms of HO2 (so *pHOx when OH is reactant):

       aqqz=2.d0*(pHOx(I,J,L)*rr(1,L) + pHOx(I,J,L)*pHOx(I,J,L)*
     & (rr(3,L)+rr(iH2O2form,L)) + rr(15,L))

       bqqz=pHOx(I,J,L)*(rr(12,L)*y(n_CH4,L)+rr(16,L)*
     & y(n_HNO3,L)+rr(iHNO3form,L)*y(nNO2,L)+rr(iHONOform,L)*y(nNO,L))
     & +rr(52,L)*y(n_HCl,L)*pHOx(I,J,L)+rr(53,L)*y(n_HOCl,L)
     & *pHOx(I,J,L)+rr(56,L)*y(nOClO,L)*pHOx(I,J,L)+
     & rr(59,L)*y(nCl,L)
     & +rr(62,L)*y(nClO,L)*pHOx(I,J,L)+rr(63,L)*y(nClO,L)
     & +rr(68,L)*y(n_HBr,L)*pHOx(I,J,L)+rr(72,L)*y(nBr,L)
     & +rr(73,L)*y(nBrO,L)+rr(81,L)*y(nBrO,L)*pHOx(I,J,L)

       ! Use OH production without O1D explicitly:
       cqqz=2.d0*ss(4,L,i,j)*y(n_H2O2,L)+ss(9,L,i,j)*y(n_HNO3,L)
     & +ss(21,L,i,j)*y(n_HOCl,L) 
     & +rr(54,L)*y(n_HCl,L)*y(nO,L)+rr(55,L)*y(n_HOCl,L)*y(nO,L)
     & +rr(57,L)*y(n_HOCl,L)*y(nCl,L)+rr(58,L)*y(nCl,L)*
     & y(n_H2O2,L)+rr(79,L)*y(nBr,L)*y(n_H2O2,L)
     & +rr(84,L)*y(n_HBr,L)*y(nO,L)
     
       ! water vapor photolysis in SRBs:
       if(PRES(L) < 10.) cqqz = cqqz + 0.5d0*SF3(I,J,L)*y(nH2O,L) 

       ! production from O1D limited to O1D amount:
       rcqqz=rr(8,L)*y(nO2,L)+rr(9,L)*y(nM,L)+
     & rr(10,L)*y(nH2O,L)+rr(11,L)*y(n_CH4,L)
       if(rcqqz > 1)then
         ratio=1.d0/rcqqz
       else
         ratio=1.d0
       endif
       cqqz=cqqz+ratio*        
     & ((2.d0*rr(10,L)*y(nH2O,L)+rr(11,L)*y(n_CH4,L))*
     & ss(2,L,I,J)*y(nO3,L))/
     & (rr(8,L)*y(nO2,L)+rr(9,L)*y(nM,L)+
     & rr(10,L)*y(nH2O,L)+rr(11,L)*y(n_CH4,L))

       sqroot=sqrt(bqqz*bqqz+4.d0*aqqz*cqqz)
       y(nHO2,L)=(sqroot-bqqz)/(2.d0*aqqz)
       y(nOH,L)=pHOx(I,J,L)*y(nHO2,L)
       temp_yHOx=y(nOH,L)+y(nHO2,L)

c Now partition HOx into OH and HO2:
c CZ: OH->HO2 reactions :
       cz=rr(2,L)*y(nO3,L)+rr(13,L)*y(n_CO,L)
     & +rr(14,L)*y(n_H2O2,L)+rr(19,L)*y(nH2,L)
     & +rr(21,L)*y(n_HCHO,L)
     & +rr(61,L)*y(nClO,L)+rr(80,L)*y(nBrO,L)

       dz=rr(4,L)*y(nO3,L)+rr(6,L)*y(nNO,L)
     & +rr(60,L)*y(nCl,L)+rr(90,L)*y(nO,L)

       if(cz+dz > 0)then
         y(nOH,L)=(dz/(cz+dz))*temp_yHOx
         if(y(nOH,L) > temp_yHOx)y(nOH,L)=temp_yHOx-1.d0
       else
         y(nOH,L)=1.d0
       endif
       y(nHO2,L)=(temp_yHOx-y(nOH,L))

c At low pressures, include loss of OH into atomic H using production
c via OH + O -> O2 + H, loss via H + O3 -> OH + O2 and
c H + O2 + M -> HO2 + M :
       if(PRES(L) < 2.d0)then
         rHprod=rr(89,L)*y(nOH,L)*y(nO,L)
         rHspecloss=y(nO3,L)*1.4d-10*exp(-470./ta(L))
         rkzero=y(nM,L)*5.7d-32*((ta(L)/300.d0)**(-1.6))
         rktot=(rkzero/(1+(rkzero/7.5d-11)))
         rHspecloss=rHspecloss+y(nO2,L)*rktot
         y(nOH,L)=y(nOH,L)-rHprod/rHspecloss
       endif

       if(y(nOH,L)  < 1)   y(nOH,L)  = 1.d0
       if(y(nHO2,L) < 1)   y(nHO2,L) = 1.d0
       if(y(nHO2,L) > 1.d9)y(nHO2,L) = 1.d9
       pHOx(I,J,L)=y(nOH,L)/y(nHO2,L)

      enddo  ! end of stratossphere loop

      return
      END SUBROUTINE HOxfam



      SUBROUTINE ClOxfam(lmax,I,J,ClOx_old)
!@sum ClOxfam Find ClOx family (Cl,ClO,OClO,Cl2,Cl2O2) partitioning 
!@+   assuming equilibrium within ClOx.
!@auth Drew Shindell

C**** GLOBAL parameters and variables:

      USE ATM_COM, only   : LTROPO
      USE RESOLUTION, only : LS1
      USE RESOLUTION, only : LM
      USE TRACER_COM, only : n_ClOx,n_HOCl,n_ClONO2,n_HCl,n_H2O2,n_CH4
      USE TRCHEM_Shindell_COM, only:pClOx,rr,y,nClO,nOClO,nCl,nCl2O2,
     &    ta,ss,nO3,nHO2,nNO3,nO,nNO,nBr,nOH,nBrO,nCH3O2,nM,nCl2,nH2,
     &    SZA,dt2,pClx,pOClOx,nNO2,which_trop,yCl2,yCl2O2

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var A,B,C,D,F,G,Q,V,X,YY,dClOx,ww,p1,p2,p3,dOClO,ratioc dummy vars
!@var ClOx_old total ClOx at start of chemical timestep
!@var rnormnum is temporary var for conservation of Cl
!@var destCl, prodCl temporary vars for simple steady state Cl calc
!@var maxl either LTROPO(I,J) or LS1-1 depending on which_trop dbparam
      integer               :: L,maxl
      integer, intent(IN)   :: lmax,I,J
      REAL*8, DIMENSION(LM) :: ClOx_old 
      real*8                :: A,B,C,D,F,G,Q,V,X,YY,dClOx,ww,p1,p2,p3,
     &                         dOClO,ratioc,rnormnum,destCl,prodCl
#ifdef TRACERS_TERP
      integer, parameter :: iClOplusClO=106,iCl2O2decomp=97,
     &                      iClOplusNO2=107
#else
      integer, parameter :: iClOplusClO=103,iCl2O2decomp=94,
     &                      iClOplusNO2=104
#endif  /* TRACERS_TERP */

      select case(which_trop)
      case(0); maxl=ltropo(I,J)
      case(1); maxl=ls1-1
      case default; call stop_model('which_trop problem 10',255)
      end select

      do L=maxl+1,lmax ! stratosphere
c Set Cl2 and default Cl2O2:
       y(nCl2,L)=yCl2(I,J,L)       ! set in chemstep
       yCl2O2(I,J,L)=0.d0          ! non-zero at low temp, see below
       y(nCl2O2,L)=yCl2O2(I,J,L)   ! non-zero at low temp, see below

c Full ClOxfam code from offline photochemistry:
       y(nClO,L)=y(n_ClOx,L)*pClOx(I,J,L)
       y(nOClO,L)=y(n_ClOx,L)*pOClOx(I,J,L)
       y(nCl,L)=y(n_ClOx,L)*pClx(I,J,L)
       if(y(n_ClOx,L) == 0) CYCLE

c Low temperature stabilizes ClO dimer, use [Cl2O2] only for
c calculating Cl amount, otherwise ignore:
       if(ta(L) <= 220.)then
         y(nCl2O2,L)=(rr(iClOplusClO,L)*y(nClO,L)*y(nClO,L))/
     &   (rr(iCl2O2decomp,L)+ss(20,L,i,j))
         if(y(nCl2O2,L) > y(nClO,L)) y(nCl2O2,L)=y(nClO,L)
         y(nClO,L)=y(nClO,L)-y(nCl2O2,L)
       endif

       A=y(nO3,L)*rr(45,L)+y(nOClO,L)*rr(47,L)+y(nHO2,L)*rr(60,L)
       B=y(nCl,L)*rr(47,L)+y(nO,L)*rr(50,L)+y(nNO,L)*
     &   rr(67,L)+y(nBr,L)*rr(74,L)+ss(19,L,i,j)
       C=y(nO,L)*rr(46,L)+y(nO3,L)*(rr(48,L)+rr(49,L))+
     &   y(nOH,L)*rr(61,L)+y(nNO,L)*rr(64,L)+y(nBrO,L)*
     &   (rr(75,L)+rr(76,L))+ss(17,L,i,j)+rr(85,L)*
     &   y(nCH3O2,L)
     &   +y(nClO,L)*(1.d-12*exp(-1590./TA(L))+3.d-11*exp(-2450./TA(L))
     &   + 3.5d-13*exp(-1370./TA(L))) 
       D=y(nO3,L)*rr(49,L)+y(nBrO,L)*rr(75,L)
       F=(rr(53,L)*y(nOH,L)*y(n_HOCl,L)+rr(55,L)*y(nO,L)*
     &   y(n_HOCl,L)+rr(65,L)*y(n_ClONO2,L)*y(nO,L))/y(n_ClOx,L)
       G=rr(62,L)*y(nOH,L)+rr(63,L)*y(nHO2,L)+rr(77,L)*
     &   y(nBrO,L)+2.d0*rr(iClOplusClO,L)*y(nClO,L)+
     &   rr(iClOplusNO2,L)*y(nNO2,L)
       Q=rr(56,L)*y(nOH,L)
       V=C-D
       X=rr(51,L)*y(nOH,L)*y(nCl2,L)+rr(54,L)*y(nO,L)*
     &   y(n_HCl,L)+ 2.d0*
     &   ss(18,L,i,j)*y(nCl2,L)+2.d0*ss(20,L,i,j)*y(nCl2O2,L)+
     &   ss(21,L,i,j)*y(n_HOCl,L)+ss(22,L,i,j)*y(n_ClONO2,L)
       X=X/y(n_ClOx,L)
       YY=rr(57,L)*y(n_HOCl,L)+rr(58,L)*y(n_H2O2,L)+rr(59,L)*
     &   y(nHO2,L)+rr(82,L)*y(n_CH4,L)+rr(83,L)*y(nH2,L)
       if((dt2*y(n_ClOx,L)) /= 0)then
         dClOx=(y(n_ClOx,L)-ClOx_old(L))/(dt2*y(n_ClOx,L))
       else
         dClOx=0.d0
       endif
       ww=A+yy+dClOx

       if(v /= 0)then 
         p1=(-(ww/v)*(B+C+G+dClOx))+A-B
         if(p1 /= 0)then
           p1=abs(-((x/v)*(B+C+G+dClOx)+F+B)/p1)
           if((B+D-Q-dClOx) /= 0)then
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
       
       if(SZA > 90.)then
         dOClO=(y(nClO,L)*D-y(nOClO,L)*(B+Q))*dt2
         y(nOClO,L)=y(nOClO,L)+dOClO
         if(y(nOClO,L) < 0.) y(nOClO,L)=0.d0
       endif 
       
       y(nCl,L)=p1*y(n_ClOx,L)
       y(nClO,L)=p2*y(n_ClOx,L)
  
       if(SZA < 90.) y(nOClO,L)=p3*y(n_ClOx,L)  
       if(y(nCl,L) < 0)   y(nCl,L)   =0.d0
       if(y(nClO,L) < 1)  y(nClO,L)  =0.d0
       if(y(nOClO,L) < 1) y(nOClO,L) =0.d0
       if(y(n_ClOx,L) < 1)y(n_ClOx,L)=1.0
           

c Normalize so that amount of ClOx doesn't change:
       if((y(nCl,L)+y(nClO,L)+y(nOClO,L)) > 0.)then
          rnormnum=y(n_ClOx,L)/(y(nCl,L)+y(nClO,L)+y(nOClO,L))
          y(nCl,L)=y(nCl,L)*rnormnum
          y(nClO,L)=y(nClO,L)*rnormnum
          y(nOClO,L)=y(nOClO,L)*rnormnum
       endif

       pClOx(I,J,L)=y(nClO,L)/y(n_ClOx,L)
       pClx(I,J,L)=y(nCl,L)/y(n_ClOx,L)
       pOClOx(I,J,L)=y(nOClO,L)/y(n_ClOx,L)

      enddo  ! end of altitude loop
      yCl2O2(I,J,:)=y(nCl2O2,:)
        
      return
      END SUBROUTINE ClOxfam



      SUBROUTINE BrOxfam(lmax,I,J)
!@sum BrOxfam Find BrOx family (Br,BrO) partitioning 
!@+   assuming equilibrium within BrOx.
!@auth Drew Shindell

C**** GLOBAL parameters and variables:
      USE RESOLUTION, only : LS1
      USE ATM_COM, only    : LTROPO
      USE TRACER_COM, only : n_BrOx,n_H2O2,n_HBr,n_HOBr,n_BrONO2
      USE TRCHEM_Shindell_COM, only:rr,y,nO3,nClO,nOClO,nNO,nO,nBr,nOH,
     &    nBrO,ss,nHO2,nNO2,pBrOx,which_trop

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var a,b,c,d,eq,f,p2,p1 dummy vars
!@var maxL either LTROPO(I,J) or LS1-1 depending on which_trop dbparam
      integer             :: L,maxl
      integer, intent(IN) :: lmax,I,J
      real*8              :: a,b,c,d,eq,f,p2,p1
#ifdef TRACERS_TERP
      integer, parameter :: iBrOplusNO2=108
#else
      integer, parameter :: iBrOplusNO2=105
#endif  /* TRACERS_TERP */
      
      select case(which_trop)
      case(0); maxl=ltropo(I,J)
      case(1); maxl=ls1-1
      case default; call stop_model('which_trop problem 11',255)
      end select

      do L=maxl+1,lmax  ! stratosphere
        a=y(nO3,L)*rr(70,L)+y(nOClO,L)*rr(74,L)
        b=y(nO,L)*rr(69,L)+y(nNO,L)*rr(71,L)+y(nClO,L)*
     &  (rr(75,L)+rr(76,L))+2*y(nBrO,L)*rr(78,L)+y(nOH,L)*
     &  rr(80,L)+ss(25,L,i,j)
        c=rr(73,L)*y(nHO2,L)+rr(77,L)*y(nClO,L)+rr(81,L)*
     &  y(nOH,L)+rr(iBrOplusNO2,L)*y(nNO2,L)    
        d=rr(72,L)*y(nHO2,L)+rr(79,L)*y(n_H2O2,L)
        eq=rr(68,L)*y(n_HBr,L)*y(nOH,L)+rr(84,L)*y(n_HBr,L)*
     &  y(nO,L)+ss(24,L,i,j)*y(n_HOBr,L)
        f=ss(23,L,i,j)*y(n_BrONO2,L)
        if(a+b /= 0)then
          p2=a/(a+b)
        else
          p2=0.d0
        endif
        if(p2 < 0)p2=0.d0
        if(p2 > 1)p2=1.d0
        p1=1-p2    
        y(nBr,L)=p1*y(n_BrOx,L)
        y(nBrO,L)=p2*y(n_BrOx,L)
        if(y(nBr,L) < 1)   y(nBr,L)   =0.d0
        if(y(nBrO,L) < 1)  y(nBrO,L)  =0.d0
        if(y(nBr,L) > 1d9) y(nBr,L)   =0.d0
        if(y(nBrO,L) > 1d9)y(nBrO,L)  =0.d0
        if(y(n_BrOx,L) < 1)y(n_BrOx,L)=1.d0
        pBrOx(I,J,L)=y(nBrO,L)/y(n_BrOx,L)
      enddo ! end of altitude loop

      return
      END SUBROUTINE BrOxfam


