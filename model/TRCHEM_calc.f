#include "rundeck_opts.h"
      SUBROUTINE chemstep(I,J,changeL)
!@sum chemstep Calculate new concentrations after photolysis & chemistry
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on chemcalc0C5.4_M23p.f from model II)
!@calls rates,chem1,chem1prn
c
C**** GLOBAL parameters and variables:
C
      USE DOMAIN_DECOMP, only : GRID, GET
      USE MODEL_COM, only       : im,jm,lm
#if (defined regional_Ox_tracers) || (defined SHINDELL_STRAT_CHEM)
     &                            ,ptop,psf,sig
#endif
      USE DYNAMICS, only        : am, byam,LTROPO
      USE GEOM, only            : BYDXYP,dxyp
#ifdef regional_Ox_tracers
     &                            ,LAT_DG,LON_DG
#endif
      USE TRACER_DIAG_COM, only : jls_OHcon,jls_H2Omr,jls_day,tajls
C    &                            ! ,ijs_OxL1
#ifdef regional_Ox_tracers
     &  ,jls_Oxloss,jls_Oxprod,ijs_Oxprod,ijs_Oxloss
#endif
      USE TRACER_COM, only: n_CH4,n_CH3OOH,n_Paraffin,n_PAN,n_Isoprene,
     &                   n_AlkylNit,n_Alkenes,n_N2O5,n_NOx,n_HO2NO2,
     &                   n_Ox,n_HNO3,n_H2O2,n_CO,n_HCHO,trm,ntm,n_N2O,
     &                   n_ClOx,n_BrOx,n_HCl,n_HOCl,n_ClONO2,n_HBr,
     &                   n_HOBr,n_BrONO2,n_CFC,ntm_chem
#ifdef regional_Ox_tracers
     &         ,NregOx,regOx_t,regOx_b,regOx_n,regOx_s,regOx_e,regOx_w
#endif  
      USE TRCHEM_Shindell_COM, only: chemrate,photrate,mass2vol,
     &                   yCH3O2,yC2O3,yXO2,yXO2N,yRXPAR,yAldehyde,
     &                   yROR,nCH3O2,nC2O3,nXO2,nXO2N,nRXPAR,
     &                   nAldehyde,nROR,nr,nn,dt2,nss,ks,dest,prod,
     &                   ny,rr,nO1D,nOH,nNO,nHO2,ta,nM,ss,
     &                   nO3,nNO2,nNO3,prnrts,jprn,iprn,lprn,ay,
     &                   prnchg,y,bymass2vol,kss,nps,kps,nds,kds,
     &                   npnr,nnr,ndnr,kpnr,kdnr,nH2O
#ifdef SHINDELL_STRAT_CHEM
     &                   ,SF3,ratioNs,ratioN2,rNO2frac,nO,nClO,nBrO
     &                   ,rNOfrac,rNOdenom,nOClO,nCl,nBr,changeCl
#endif
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var changeL change due to chemistry in mass/time
!@var I,J passed horizontal spatial indicies
!@var L,iter,Lz dummy loop variable
!@var maxl highest level with chemistry, maxT top of troposphere
!@var qqqCH3O2,CH3O2loss,XO2_NO,XO2N_HO2,RXPAR_PAR,ROR_CH2,C2O3prod,
!@+   C2O3dest,XO2prod,XO2dest,XO2_XO2,XO2Nprod,XO2Ndest,RXPARprod,
!@+   RXPARdest,Aldehydeprod,Aldehydedest,RORprod,RORdest,total,
!@+   rnewval,dNOx,ratio,sumD,newD,ratioD,newP,ratioP,changeA
!@+   sumP dummy temp variables
!@+   sumN,sumC,sumH,sumB,sumO,sumA variables for O3 catalytic diags
!@var tempiter,tempiter2 temp vars for equilibrium calcs iterations
!@var changeX temporary variable for equil calcs
!@var iHO2NO2form HO2NO2 formation reaction
!@var iN2O5form N2O5 formation reaction
!@var iPANform PAN formation reaction
!@var iHO2NO2_OH HO2NO2 oxidation by OH reaction
!@var iHO2NO2decomp HO2NO2 decomposition reaction
!@var iN2O5decomp N2O5 decomposition reaction
!@var iPANdecomp PAN decomposition reaction
!@var rMAbyM is airmass over air concentration
!@var dxbym2v is dxyp over mass2volume
!@var sv_changeN2O N2O change without portion making N2 (for N cons)
!@var vClONO2, vBrONO2 temporary vars within N conservation
!@var changeH2O chemical change in H2O
!@var Oxcorr account for Ox change from within NOx partitioning
!@var rNO3prod,rNO2prod,rNOprod to acct for dOx from NOx partitioning
      REAL*8, DIMENSION(LM,ntm)     :: changeL
      REAL*8, DIMENSION(LM) :: rMAbyM,sv_changeN2O,changeH2O,Oxcorr
      INTEGER, INTENT(IN) :: I,J
      INTEGER L,iter,maxl,igas,maxT,Lz
#ifdef SHINDELL_STRAT_CHEM
      INTEGER, PARAMETER :: iHO2NO2form=98,iN2O5form=99,
     &iPANform=101,iHO2NO2_OH=18,iHO2NO2decomp=91,iN2O5decomp=92
     &,iPANdecomp=29
!@param JN J around 30 N
!@param JS J around 30 S
!@param JNN,JSS Js for "high-lat" definition
      INTEGER, PARAMETER :: JS = JM/3 + 1, JN = 2*JM/3
      INTEGER, PARAMETER :: JNN = 5*JM/6, JSS= JM/6 + 1
#else
      INTEGER, PARAMETER :: iHO2NO2form=51,iN2O5form=52,
     &iPANform=54,iHO2NO2_OH=18,iHO2NO2decomp=45,iN2O5decomp=46
     &,iPANdecomp=29
#endif
      REAL*8 qqqCH3O2,CH3O2loss,XO2_NO,XO2N_HO2,RXPAR_PAR,ROR_CH2,
     & C2O3prod,C2O3dest,XO2prod,XO2dest,XO2_XO2,XO2Nprod,XO2Ndest,
     & RXPARprod,RXPARdest,Aldehydeprod,Aldehydedest,RORprod,RORdest,
     & total,rnewval,dNOx,ratio,sumD,newD,ratioD,newP,ratioP,
     & changeA,sumP,tempiter,tempiter2,sumC,sumN,sumH,sumB,sumO,sumA,
     & dxbym2v,changeX,vClONO2,vBrONO2,conc2mass,rNO3prod,rNO2prod,
     & rNOprod,changeAldehyde
#if (defined regional_Ox_tracers) || (defined SHINDELL_STRAT_CHEM)
!@var PRES local nominal pressure for regional Ox tracers
      REAL*8, DIMENSION(LM) :: PRES
#endif
#ifdef regional_Ox_tracers
!@var Oxloss Ox chemical loss for use with regional Ox tracers
!@var Oxprod Ox chemical production for use with regional Ox tracers
!@var nREG index of regional Ox tracer number
      REAL*8, DIMENSION(LM) :: Oxloss, Oxprod
      INTEGER nREG
#endif
!@param nlast either ntm or ntm-NregOx for chemistry loops
      INTEGER, PARAMETER :: nlast =
#ifdef regional_Ox_tracers
     &                              ntm_chem-NregOx
#else
     &                              ntm_chem
#endif

      INTEGER J_0H, J_1H

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )

C
C     TROPOSPHERIC CHEMISTRY ONLY or TROP+STRAT:
#ifdef SHINDELL_STRAT_CHEM
      maxl=LM
#else
      maxl=LTROPO(I,J)
#endif
      maxT=LTROPO(I,J)

#if (defined regional_Ox_tracers) || (defined SHINDELL_STRAT_CHEM)
      do L=1,maxL
       PRES(L)=SIG(L)*(PSF-PTOP)+PTOP
      end do
#endif
      do L=1,maxT
       y(nCH3O2,L)   =yCH3O2(I,J,L)
       y(nC2O3,L)    =yC2O3(I,J,L)
       y(nXO2,L)     =yXO2(I,J,L)
       y(nXO2N,L)    =yXO2N(I,J,L)
       y(nRXPAR,L)   =yRXPAR(I,J,L)
       y(nAldehyde,L)=yAldehyde(I,J,L)
       y(nROR,L)     =yROR(I,J,L)
      enddo
#ifdef SHINDELL_STRAT_CHEM
      do L=maxT+1,LM
       y(nCH3O2,L)=yCH3O2(I,J,L)
      enddo
#endif
C
cc    calculate reaction rates with present concentrations
      call rates(maxl,I,J)
c
c     chem1 call sample:
c     (klist,l,numel,nlist,ndlist,rate,change,multip)
c              numel=number of elements in reaction list nlist (1 or 2)
c              change=dest or prod array
c              multip=1(prod) or -1(dest)
c     chemical destruction:
      call chem1(kdnr,maxl,2,nn,ndnr,chemrate,dest,-1)
c     chemical production:
      call chem1(kpnr,maxl,2,nnr,npnr,chemrate,prod,1)
c     photolytic destruction:
      call chem1(kds,maxl,1,ks,nds,photrate,dest,-1)
c     photolytic production:
      call chem1(kps,maxl,2,kss,nps,photrate,prod,1)
C
c     Oxidation of Isoprene and Alkenes produces less than one
c      HCHO, Alkenes, and CO per rxn, correct here following Houweling
      do L=1,maxl
       prod(n_CO,L)=prod(n_CO,L)-0.63d0*chemrate(35,L)
       prod(n_HCHO,L)=prod(n_HCHO,L)-0.36d0*chemrate(35,L)
       prod(n_HCHO,L)=prod(n_HCHO,L)-0.39d0*chemrate(30,L)
       prod(n_Alkenes,L)=prod(n_Alkenes,L)-0.42d0*chemrate(30,L)
       prod(n_HCHO,L)=prod(n_HCHO,L)-0.10d0*chemrate(31,L)
       prod(n_Alkenes,L)=prod(n_Alkenes,L)-0.45d0*chemrate(31,L)
      enddo
c
#ifdef SHINDELL_STRAT_CHEM
c     modify Cl source from CFC photolysis to account for other
c      sources (e.g. other CFCs, HCFCs, methyl chloride)
      do L=maxT+1,LM
       prod(n_ClOx,L)=prod(n_ClOx,L)+ss(26,L,I,J)*y(n_CFC,L)*dt2*0.75d0
c     add Br source using CFC photolysis as proxy
       prod(n_HBr,L)=prod(n_HBr,L)+ss(26,L,I,J)*y(n_CFC,L)*dt2*1.d-3
      enddo
c
c     Remove some lower strat polar HNO3 in PSCs
      do L=1,maxL
       IF(PRES(L).lt.245.d0 .and. PRES(L).ge.31.6d0)THEN
         if(ta(L).lt.198.d0)then
C this was:if(J.le.17.or.J.ge.39)dest(n_HNO3,L)=
           if(J.le.JS+1.or.J.ge.JNN+1)dest(n_HNO3,L)=
     *     dest(n_HNO3,L)-0.05d0*y(n_HNO3,L)
         endif
       END IF
      enddo
#endif
c
c     set CH3O2 values (concentration = production/specific loss)
      do L=1,maxl
        iter=1
        qqqCH3O2=(rr(11,L)*y(nO1D,L)+rr(12,L)*y(nOH,L))
     &  *y(n_CH4,L)+rr(23,L)*y(n_CH3OOH,L)*y(nOH,L)

       tempiter=rr(20,L)*y(nNO,L)+rr(22,L)*y(nHO2,L)
 10    CH3O2loss=tempiter+rr(27,L)*yCH3O2(I,J,L)
       if(CH3O2loss.gt.1d-7)then
          y(nCH3O2,L)=qqqCH3O2/CH3O2loss
        else
          y(nCH3O2,L)=1.d0
        endif
        yCH3O2(I,J,L)=y(nCH3O2,L)
        iter=iter+1
        if(iter.le.7)goto10 ! replace with while loop?
      enddo
c
      do L=1,maxT
c
c       Set C2O3, XO2, XO2N, RXPAR, Aldehyde & ROR values
c
c        First set various specific loss rates
         XO2_NO=y(nNO,L)*4.2d-12*exp(180.d0/ta(L))
         XO2N_HO2=y(nHO2,L)*y(nNO,L)*rr(44,L)*
     &   rr(43,L)/XO2_NO
         RXPAR_PAR=y(n_Paraffin,L)*8.d-11
         ROR_CH2=1.6d3
c
c       Set value for C2O3
        iter=1
        C2O3prod=rr(38,L)*yAldehyde(I,J,L)*y(nOH,L)+
     &  (rr(29,L)*y(nM,L)+ss(15,L,I,J))*y(n_PAN,L)+0.15d0*
     &  rr(31,L)*y(nO3,L)*y(n_Isoprene,L)
        tempiter=rr(39,L)*y(nNO,L)+rr(iPANform,L)*y(nNO2,L)+
     &   rr(41,L)*y(nHO2,L)
 20     C2O3dest=tempiter+rr(40,L)*yC2O3(I,J,L)
        if(C2O3dest.gt.1d-7)then
          y(nC2O3,L)=(C2O3prod/C2O3dest)
        else
          y(nC2O3,L)=1.d0
        endif
        yC2O3(I,J,L)=y(nC2O3,L)
        iter=iter+1
        if(iter.le.7)goto20     ! replace with while loop?
c
c       Set value for XO2
        iter=1
        XO2prod=ss(16,L,I,J)*yAldehyde(I,J,L)+
     &  y(nC2O3,L)*(rr(39,L)*y(nNO2,L)+rr(40,L)*
     &  y(nC2O3,L)*2.d0+rr(41,L)*y(nHO2,L))
     &  +rr(42,L)*yROR(I,J,L)*0.96d0
     &  +y(nOH,L)*(rr(37,L)*y(n_Paraffin,L)*0.87d0+rr(34,L)*
     &  y(n_Alkenes,L)+rr(30,L)*y(n_Isoprene,L)*0.85d0+
     &  rr(33,L)*y(n_AlkylNit,L))+
     &  y(nO3,L)*(rr(35,L)*y(n_Alkenes,L)*0.29d0+
     &  rr(31,L)*y(n_Isoprene,L)*0.18d0)
        tempiter=XO2_NO+rr(43,L)*y(nHO2,L)
        tempiter2=1.7d-14*exp(1300.d0/ta(L))
 30     XO2_XO2=yXO2(I,J,L)*tempiter2
        XO2dest=tempiter+XO2_XO2
        if(XO2dest.gt.1.d-7.and.ss(16,L,I,J).gt.5.d-5)then
          y(nXO2,L)=(XO2prod/XO2dest)
        else
          y(nXO2,L)=1.d0
        endif
        yXO2(I,J,L)=y(nXO2,L)
        iter=iter+1
        if(iter.le.7)goto30   ! replace with while loop?
c
c       Set value for XO2N
        XO2Nprod=rr(37,L)*y(n_Paraffin,L)*y(nOH,L)*0.13d0+
     &  rr(42,L)*yROR(I,J,L)*0.04d0+rr(30,L)*y(n_Isoprene,L)*
     &  y(nOH,L)*0.15d0
        XO2Ndest=XO2N_HO2+rr(44,L)*y(nNO,L)
        if(XO2Ndest.gt.1.d-7)then
          y(nXO2N,L)=(XO2Nprod/XO2Ndest)
        else
          y(nXO2N,L)=1.d0
        endif
        yXO2N(I,J,L)=y(nXO2N,L)
c
c       Set value for RXPAR
        RXPARprod=rr(37,L)*y(n_Paraffin,L)*y(nOH,L)*0.11d0+
     &  rr(34,L)*yROR(I,J,L)*2.1d0+rr(35,L)*y(n_Alkenes,L)*
     &  y(nO3,L)*0.9d0
        RXPARdest=RXPAR_PAR
        if(RXPARdest.gt.0.d0)then
          y(nRXPAR,L)=(RXPARprod/RXPARdest)
        else
          y(nRXPAR,L)=1.d0
        endif
        yRXPAR(I,J,L)=y(nRXPAR,L)
c
c       Set value for Aldehyde
        Aldehydeprod=rr(37,L)*y(n_Paraffin,L)*y(nOH,L)*0.11d0+
     &  rr(34,L)*y(n_Alkenes,L)*y(nOH,L)+
     &  rr(42,L)*yROR(I,J,L)*1.1d0+rr(35,L)*y(n_Alkenes,L)*
     &  y(nO3,L)*0.44d0
        Aldehydedest=rr(38,L)*y(nOH,L)+ss(16,L,I,J)
c
c       Check for equilibrium
        if(Aldehydedest*y(nAldehyde,L)*dt2.lt.y(nAldehyde,L))then
         changeAldehyde=
     &   (Aldehydeprod-y(nAldehyde,L)*Aldehydedest)*dt2
         if(changeAldehyde.gt.y(nAldehyde,L))
     &   changeAldehyde=y(nAldehyde,L)
         y(nAldehyde,L)=y(nAldehyde,L)+changeAldehyde
         if(y(nAldehyde,L).lt.0.d0)y(nAldehyde,L)=1.d0
        else
          y(nAldehyde,L)=(Aldehydeprod/(Aldehydedest+0.5d-5))
        endif
        yAldehyde(I,J,L)=y(nAldehyde,L)
c
c       Set value for ROR
        RORprod=rr(37,L)*y(n_Paraffin,L)*y(nOH,L)*0.76d0+
     &  rr(42,L)*yROR(I,J,L)*0.02d0
        RORdest=rr(42,L)+ROR_CH2
        if(RORdest.gt.0.d0)then
          y(nROR,L)=(RORprod/RORdest)
        else
          y(nROR,L)=1.d0
        endif
        yROR(I,J,L)=y(nROR,L)
c
c       Add in parrafin loss term via rxpar reaction
        dest(n_Paraffin,L)=dest(n_Paraffin,L)-y(nRXPAR,L)*RXPAR_PAR*dt2
c
c       Add in CH3OOH production via XO2N + HO2
        prod(n_CH3OOH,L)=prod(n_CH3OOH,L)+XO2N_HO2*y(nXO2N,L)*dt2
c
      end do  ! L
c
c     If NOx in equil with N2O5, HO2NO2 or PAN, remove from changes
      do L=1,maxl
       if(-dest(n_N2O5,L).ge.y(n_N2O5,L).or.
     &chemrate(iN2O5form,L).gt.y(n_NOx,L))then
        dest(n_NOx,L)=dest(n_NOx,L)+2.d0*chemrate(iN2O5form,L)
        prod(n_NOx,L)=prod(n_NOx,L)-2.d0*(chemrate(iN2O5decomp,L)
     &   +photrate(7,L))
       endif
       if(-dest(n_HO2NO2,L).ge.y(n_HO2NO2,L).or.
     &chemrate(iHO2NO2form,L).gt.y(n_NOx,L))then
        dest(n_NOx,L)=dest(n_NOx,L)+chemrate(iHO2NO2form,L)
        prod(n_NOx,L)=prod(n_NOx,L)-(chemrate(iHO2NO2_OH,L)+
     &   chemrate(iHO2NO2decomp,L)+photrate(10,L)+photrate(11,L))
       endif
       if(-dest(n_PAN,L).ge.y(n_PAN,L).or.
     &chemrate(iPANform,L).gt.y(n_NOx,L))then
        dest(n_NOx,L)=dest(n_NOx,L)+chemrate(iPANform,L)
        prod(n_NOx,L)=prod(n_NOx,L)-(chemrate(iPANdecomp,L)+
     &   photrate(15,L))
       endif
#ifdef SHINDELL_STRAT_CHEM
c     If BrOx in equil with HOBr or BrONO2, remove from changes
       if(-dest(n_HOBr,L).ge.y(n_HOBr,L).or.
     &chemrate(73,L).gt.0.5d0*y(n_BrOx,L))then
        dest(n_BrOx,L)=dest(n_BrOx,L)+chemrate(73,L)
        prod(n_BrOx,L)=prod(n_BrOx,L)-photrate(24,L)
       endif
       if(-dest(n_BrONO2,L).ge.y(n_BrONO2,L).or.
     &chemrate(104,L).gt.0.5d0*y(n_BrOx,L))then
        dest(n_BrOx,L)=dest(n_BrOx,L)+chemrate(104,L)
        prod(n_BrOx,L)=prod(n_BrOx,L)-photrate(23,L)
        dest(n_NOx,L)=dest(n_NOx,L)+chemrate(104,L)
        prod(n_NOx,L)=prod(n_NOx,L)-photrate(23,L)
       endif
c     If ClOx in equil with HOCl or ClONO2, remove from changes
       if(-dest(n_HOCl,L).ge.y(n_HOCl,L).or.
     &chemrate(63,L).gt.y(n_ClOx,L))then
        dest(n_ClOx,L)=dest(n_ClOx,L)+chemrate(63,L)
        prod(n_ClOx,L)=prod(n_ClOx,L)-(photrate(21,L)+chemrate(55,L))
       endif
       if(-dest(n_ClONO2,L).ge.y(n_ClONO2,L).or.
     &chemrate(103,L).gt.0.8d0*y(n_ClOx,L))then
        dest(n_ClOx,L)=dest(n_ClOx,L)+chemrate(103,L)
        prod(n_ClOx,L)=prod(n_ClOx,L)-(photrate(22,L)+chemrate(65,L))
        dest(n_NOx,L)=dest(n_NOx,L)+chemrate(103,L)
        prod(n_NOx,L)=prod(n_NOx,L)-(photrate(22,L)+chemrate(65,L))
       endif
#endif
      enddo
c
#ifdef SHINDELL_STRAT_CHEM
c      Calculate water vapor change
       do L=maxT+1,maxl
        if(y(nO1D,L).gt.0.d0)changeH2O(L)=(2*y(n_CH4,L)*
     *   (rr(11,L)*y(nO1D,L)+rr(12,L)*y(nOH,L))-SF3(I,J,L))*dt2
       enddo
c
c     Calculate ozone change due to within NOx partitioning
      do L=1,LM
c      account for NO2 and NO ozone destruction
       rNO2prod=rr(18,L)*y(nOH,L)*y(n_HO2NO2,L)+
     *  rr(91,L)*y(n_HO2NO2,L)+ss(9,L,I,J)*y(n_HNO3,L)+
     *  ss(10,L,I,J)*y(n_HO2NO2,L)+ss(23,L,I,J)*y(n_BrONO2,L)
       rNOprod=rr(87,L)*y(n_N2O,L)*y(nO1D,L)
       rNO3prod=rr(65,L)*y(nO,L)*y(n_ClONO2,L)+
     *  ss(7,L,I,J)*y(n_N2O5,L)+ss(11,L,I,J)*y(n_HO2NO2,L)+
     *  ss(22,L,I,J)*y(n_ClONO2,L)
c      add in production of NO and NO2 from NO3
       rNO3prod=rNO3prod*ss(6,L,I,J)/(ss(5,L,I,J)+ss(6,L,I,J)+1.d0)
       rNO2prod=rNO2prod+rNO3prod
       rNOprod=rNOprod+rNO3prod
c
       ratioNs=rNO2prod/rNOprod
       ratioN2=y(nNO2,L)/y(nNO,L)
       if(ratioNs.gt.ratioN2)then !excess NO2 production
c
c      account for NO2 that then goes via NO2+O->NO+O2, NO2->NO+O
       rNO2frac=(rr(26,L)*y(nO,L)-ss(1,L,I,J))/
     *  (rr(97,L)*y(nOH,L)+
     *  rr(98,L)*y(nHO2,L)+rr(99,L)*y(nNO3,L)+
     *  rr(103,L)*y(nClO,L)+rr(104,L)*y(nBrO,L)+
     *  rr(26,L)*y(nO,L)+ss(1,L,I,J))
       Oxcorr(L)=(rNO2prod-rNOprod)*rNO2frac*dt2*y(nNO,L)/y(n_NOx,L)
       if(Oxcorr(L).gt.-1.d18.and.Oxcorr(L).lt.1.d18)then
        dest(n_Ox,L)=dest(n_Ox,L)-Oxcorr(L)
       else
        write(6,*) 'Oxcorr fault NO2:',ratioNs,ratioN2,rNO2frac,
     &   rNO2prod,rNOprod
       endif
c
       else !excess NO prodcution
c
c      account for NO that then goes via NO+O3->NO2+O2 or NO+O+M->NO2+M
       rNOfrac=(rr(5,L)*y(nO3,L)+rr(95,L)*y(nO,L))
       rNOdenom=(rr(5,L)*y(nO3,L)+rr(95,L)*y(nO,L)+
     *  rr(6,L)*y(nHO2,L)+rr(44,L)*y(nXO2N,L)+1.d0)
        if(l.le.maxT)then
c        Troposphere
         rNOdenom=rNOdenom+rr(20,L)*yCH3O2(I,J,L)
     &   +rr(39,L)*y(nC2O3,L)+4.2d-12*exp(180/ta(L))*y(nXO2,L)
        else
c        Stratosphere
         rNOdenom=rNOdenom+rr(64,L)*y(nClO,L)+
     &   rr(67,L)*y(nOClO,L)+rr(71,L)*y(nBrO,L)
        endif
       rNOfrac=rNOfrac/rNOdenom
       Oxcorr(L)=(rNOprod-rNO2prod)*rNOfrac*dt2*y(nNO2,L)/y(n_NOx,L)
       if(Oxcorr(L).gt.-1.d18.and.Oxcorr(L).lt.1.d18)then
        dest(n_Ox,L)=dest(n_Ox,L)-Oxcorr(L)
       else
        write(6,*) 'Oxcorr fault NO:',I,J,L,ratioNs,ratioN2,rNOfrac,
     &   rNO2prod,rNOprod,y(nNO2,L),y(nNO,L),rNOdenom,y(nO,L),y(nO3,L)
       endif
c
       endif
      enddo
#endif
c
c     print rxtn rates if desired (chem1prn : argument before multip is
c     index = number of call):
c
       if(prnrts.and.J.eq.jprn.and.I.eq.iprn)then
        do igas=1,nlast
         total=0.d0
         write(6,108) ' Species: ',ay(igas)
c
         call chem1prn
     *   (kdnr,2,nn,ndnr,chemrate,1,-1,igas,total,maxl,I,J)
c
         if(igas.eq.n_NOx)then
         if(-dest(n_HO2NO2,lprn).ge.y(n_HO2NO2,lprn).or.
     &chemrate(iHO2NO2form,lprn).gt.y(n_NOx,lprn))
     *    write(6,110)'loss by reaction    (HO2NO2 formation) removed',
     *    chemrate(iHO2NO2form,lprn)
         if(-dest(n_N2O5,lprn).ge.y(n_N2O5,lprn).or.
     &chemrate(iN2O5form,lprn).gt.y(n_NOx,lprn))
     *    write(6,110)'losses by reaction    (N2O5 formation) removed',
     *    2*chemrate(iN2O5form,lprn)
         if(-dest(n_PAN,lprn).ge.y(n_PAN,lprn).or.
     &chemrate(iPANform,lprn).gt.y(n_NOx,lprn))
     *    write(6,110)'losses by reaction     (PAN formation) removed',
     *    chemrate(iPANform,lprn)
         endif
         call chem1prn
     &   (kpnr,2,nnr,npnr,chemrate,2,1,igas,total,maxl,I,J)
         if(igas.eq.n_NOx)then
         if(-dest(n_HO2NO2,lprn).ge.y(n_HO2NO2,lprn).or.
     &chemrate(iHO2NO2form,lprn).gt.y(n_NOx,lprn))
     *    write(6,110)'gain by reactions destroying HO2NO2) rmoved  ',
     *    (rr(iHO2NO2_OH,lprn)*y(nOH,L)+
     *    rr(iHO2NO2decomp,lprn)*y(nM,lprn)+
     *    ss(10,lprn,I,J)+ss(11,lprn,I,J))*y(n_HO2NO2,lprn)*dt2
         if(-dest(n_N2O5,lprn).ge.y(n_N2O5,lprn).or.
     &chemrate(iN2O5form,lprn).gt.y(n_NOx,lprn))
     *        write(6,110)
     *        'gains by reaction    (N2O5 decomposition) removed',
     *        2.d0*chemrate(iN2O5decomp,lprn)
         if(-dest(n_PAN,lprn).ge.y(n_PAN,lprn).or.
     &chemrate(iPANform,lprn).gt.y(n_NOx,lprn))
     *    write(6,110)'gain by reaction    (from PAN) removed',
     *    chemrate(iPANdecomp,lprn)
        endif
c
         call chem1prn(kds,1,ks,nds,photrate,3,-1,igas,total,maxl,I,J)
c
         call chem1prn(kps,2,kss,nps,photrate,4,1,igas,total,maxl,I,J)
#ifdef SHINDELL_STRAT_CHEM
       if(igas.eq.n_Ox)write(6,110)'Ox change due to within NOx rxns  ',
     *    -Oxcorr(lprn)
#endif
c
         if(igas.eq.n_NOx)then
         if(-dest(n_N2O5,lprn).ge.y(n_N2O5,lprn).or.
     &chemrate(iN2O5form,lprn).gt.y(n_NOx,lprn))
     *    write(6,110)'gains by reaction 7 (N2O5 photolysis) removed',
     *    ss(7,lprn,I,J)*y(n_N2O5,lprn)*dt2
         if(-dest(n_N2O5,lprn).ge.y(n_N2O5,lprn).or.
     &chemrate(iN2O5form,lprn).gt.y(n_NOx,lprn))
     *        write(6,110)
     *        'net change due to N2O5 is ',
     *        2.d0*(y(n_N2O5,lprn)-(rr(iN2O5form,lprn)*y(nNO3,lprn)*
     *        y(nNO2,lprn))/
     *        (rr(iN2O5decomp,lprn)*y(nM,lprn)+ss(7,lprn,I,J)))
         if(-dest(n_HO2NO2,lprn).ge.y(n_HO2NO2,lprn).or.
     &chemrate(iHO2NO2form,lprn).gt.y(n_NOx,lprn))
     *    write(6,110)'gain by rxns 10 & 11 (HO2NO2 photolysis) rmoved'
     *    ,(ss(10,lprn,I,J)+ss(11,lprn,I,J))*y(n_HO2NO2,lprn)*dt2
         if(-dest(n_HO2NO2,lprn).ge.y(n_HO2NO2,lprn).or.
     &chemrate(iHO2NO2form,lprn).gt.y(n_NOx,lprn))
     *        write(6,110)'net change due to HO2NO2 is ',
     *        y(n_HO2NO2,lprn)-((rr(iHO2NO2form,lprn)*y(nHO2,lprn)*
     *       y(nNO2,lprn))/(rr(iHO2NO2_OH,lprn)*
     *     y(nOH,lprn)+rr(iHO2NO2decomp,lprn)*y(nM,lprn)
     *     +ss(10,lprn,I,J)+ss(11,lprn,I,J)))
         if(-dest(n_PAN,lprn).ge.y(n_PAN,lprn).or.
     &chemrate(iPANform,lprn).gt.y(n_NOx,lprn))
     *        write(6,110)'net change due to PAN is ',
     *        y(n_PAN,lprn)-((rr(iPANform,lprn)*y(nC2O3,lprn)*
     *        y(nNO2,lprn))/
     *        (rr(iPANdecomp,lprn)*y(nM,lprn)+ss(15,lprn,I,J)))
         endif
         if(igas.eq.n_Ox.or.igas.eq.n_NOx)total=
     *    100.d0*(dest(igas,lprn)+prod(igas,lprn))/y(igas,lprn)
#ifdef SHINDELL_STRAT_CHEM
        if(igas.eq.n_BrOx)then
         if(-dest(n_HOBr,lprn).ge.y(n_HOBr,lprn).or.
     &chemrate(73,lprn).gt.0.5d0*y(n_BrOx,lprn))then
         write(6,110)'gain by rxns 24 (HOBr photolysis) removed'
     *    ,ss(24,lprn,i,j)*y(n_HOBr,lprn)*dt2
         write(6,110)'prod by rxn 73 removed'
     *    ,chemrate(73,lprn)
         endif
         if(-dest(n_BrONO2,lprn).ge.y(n_BrONO2,lprn).or.
     &chemrate(104,lprn).gt.0.5d0*y(n_BrOx,lprn))then
         write(6,110)'gain by rxns 23 (BrONO2 photolysis) removed'
     *    ,ss(23,lprn,i,j)*y(n_BrONO2,lprn)*dt2
         write(6,110)'prod by rxn 104 removed'
     *    ,chemrate(104,lprn)
         endif
        endif
        if(igas.eq.n_NOx)then
         if(-dest(n_BrONO2,lprn).ge.y(n_BrONO2,lprn).or.
     &chemrate(104,lprn).gt.0.5d0*y(n_BrOx,lprn))then
         write(6,110)'gain by rxns 23 (BrONO2 photolysis) removed'
     *    ,ss(23,lprn,i,j)*y(n_BrONO2,lprn)*dt2
         write(6,110)'prod by rxn 104 removed'
     *    ,chemrate(104,lprn)
         endif
        endif
        if(igas.eq.n_ClOx)then
         if(-dest(n_HOCl,lprn).ge.y(n_HOCl,lprn).or.
     &chemrate(63,lprn).gt.y(n_ClOx,lprn))then
         write(6,110)'gain by rxn 21 (HOCl photolysis) removed'
     *    ,ss(21,lprn,i,j)*y(n_HOCl,lprn)*dt2
         write(6,110)'gain by rxn 55 removed'
     *    ,chemrate(55,lprn)
         write(6,110)'prod by rxn 63 removed'
     *    ,chemrate(63,lprn)
         endif
         if(-dest(n_ClONO2,lprn).ge.y(n_ClONO2,lprn).or.
     &chemrate(103,lprn).gt.0.8d0*y(n_ClOx,lprn))then
         write(6,110)'gain by rxn 22 (ClONO2 photolysis) removed'
     *    ,ss(22,lprn,i,j)*y(n_ClONO2,lprn)*dt2
         write(6,110)'gain by rxn 65 removed'
     *    ,chemrate(65,lprn)
         write(6,110)'prod by rxn 103 removed'
     *    ,chemrate(103,lprn)
         endif
        endif
        if(igas.eq.n_NOx)then
         if(-dest(n_ClONO2,lprn).ge.y(n_ClONO2,lprn).or.
     &chemrate(103,lprn).gt.0.8d0*y(n_ClOx,lprn))then
         write(6,110)'gain by rxn 22 (ClONO2 photolysis) removed'
     *    ,ss(22,lprn,i,j)*y(n_ClONO2,lprn)*dt2
         write(6,110)'gain by rxn 65 removed'
     *    ,chemrate(65,lprn)
         write(6,110)'prod by rxn 103 removed'
     *    ,chemrate(103,lprn)
         endif
        endif
#endif
         if(igas.eq.n_CH3OOH) write(6,'(a48,a6,e10.3)')
     *    'production from XO2N + HO2 ','dy = ',
     *    y(nHO2,lprn)*y(nNO,lprn)*rr(44,lprn)*rr(43,lprn)/
     *    (y(nNO,lprn)*4.2d-12*exp(180.d0/ta(lprn)))
     *    *y(nXO2N,lprn)*dt2
         if(igas.eq.n_Paraffin) write(6,'(a48,a6,e10.3)')
     *    'destruction from RXPAR ',
     *    'dy = ',-y(nRXPAR,lprn)*y(n_Paraffin,lprn)*8.d-11*dt2
         write(6,118) ' Total change in ',ay(igas),
     *  ' is ',total,' percent; dy= ',dest(igas,lprn)+prod(igas,lprn)
         write(6,*)
        enddo ! igas
       endif  ! chem diags
 108  format(a10,2x,a8)
 110  format(a68,e10.3)
 118  format(a17,a8,a4,f10.0,a14,e12.3)
c
      if(prnchg.and.J.eq.jprn.and.I.eq.iprn)then
#ifdef SHINDELL_STRAT_CHEM
       write(*,*) 'Percentage ozone loss per cycle at I,J:',I,J
       write(*,'(a41,a56)') '  L   ClOx    NOx     HOx     BrOx    Ox ',
     *    '   NO2+O     NO+O3     ClO+O     Cl+O3    NO2+hv    net'
       do Lz=LM,12,-1
        sumC=rr(45,Lz)*y(nCl,Lz)*y(nO3,Lz)+
     *   rr(46,Lz)*y(nClO,Lz)*y(nO,Lz)+
     *   rr(49,Lz)*y(nClO,Lz)*y(nO3,Lz)+
     *   rr(50,Lz)*y(nOClO,Lz)*y(nO,Lz)
c     *   -ss(17,Lz,i,j)*y(nClO,Lz)
        sumN=rr(5,Lz)*y(nNO,Lz)*y(nO3,Lz)+
     *   rr(26,Lz)*y(nNO2,Lz)*y(nO,Lz)+
     *   rr(7,Lz)*y(nNO2,Lz)*y(nO3,Lz)+
     *   rr(95,Lz)*y(nNO,Lz)*y(nO,Lz)
     *   -ss(1,Lz,i,j)*y(nNO2,Lz)
        sumH=rr(2,Lz)*y(nOH,Lz)*y(nO3,Lz)+
     *   rr(4,Lz)*y(nHO2,Lz)*y(nO3,Lz)+
     *   rr(89,Lz)*y(nOH,Lz)*y(nO,Lz)+
     *   rr(90,Lz)*y(nHO2,Lz)*y(nO,Lz)
        sumB=rr(69,Lz)*y(nBrO,Lz)*y(nO,Lz)+
     *   rr(70,Lz)*y(nBr,Lz)*y(nO3,Lz)
c     *   -ss(25,Lz,i,j)*y(nBrO,Lz)
        sumO=2*rr(88,Lz)*y(nO,Lz)*y(nO3,Lz)
        sumA=sumC+sumN+sumH+sumB+sumO
        write(*,'(i3,1x,5(f7.2,1x),6(e9.2,1x))') Lz,100.d0*sumC/sumA,
     *   100.d0*sumN/sumA,100.d0*sumH/sumA,100.d0*sumB/sumA,
     *   100.d0*sumO/sumA,rr(26,Lz)*y(nNO2,Lz)*y(nO,Lz),
     *   rr(5,Lz)*y(nNO,Lz)*y(nO3,Lz),
     *   rr(46,Lz)*y(nClO,Lz)*y(nO,Lz),
     *   rr(45,Lz)*y(nCl,Lz)*y(nO3,Lz),ss(1,Lz,i,j)*y(nNO2,Lz),sumA
       enddo
       write(*,*) ' '
#endif
      write(6,'(a35,3(2x,i2))')
     * ' Total change by species at I, J, L',i,j,lprn
      endif
c
cc    Loop to calculate tracer changes
      do L=1,maxl
        rMAbyM(L)=AM(L,I,J)/y(nM,L)
      enddo


      do igas=1,nlast
       dxbym2v=dxyp(J)/mass2vol(igas)
       do L=1,maxl
         conc2mass=rMAbyM(L)*dxbym2v
c
         changeL(L,igas)=
     &   (dest(igas,L)+prod(igas,L))*conc2mass
c
#ifdef SHINDELL_STRAT_CHEM
c         For regional tracers
c         if(igas.eq.n_Ox)then
c           Oxloss(L)=dest(igas,L)*conc2mass
c           Oxprod(L)=prod(igas,L)*conc2mass
c         end if
#endif
c
#ifdef regional_Ox_tracers     
         if(igas.eq.n_Ox)then
           Oxloss(L)=dest(igas,L)*dxyp(J)*AM(L,I,J)*bymass2vol(igas)
     *     /y(nM,L)
           Oxprod(L)=prod(igas,L)*dxyp(J)*AM(L,I,J)*bymass2vol(igas)
     *     /y(nM,L)
           TAJLS(J,L,jls_Oxprod)=TAJLS(J,L,jls_Oxprod)+prod(igas,L)
           TAJLS(J,L,jls_Oxloss)=TAJLS(J,L,jls_Oxloss)+dest(igas,L)
           TAIJS(I,J,ijs_Oxprod)=TAIJS(I,J,ijs_Oxprod)+prod(igas,L)
           TAIJS(I,J,ijs_Oxloss)=TAIJS(I,J,ijs_Oxloss)+dest(igas,L)
         end if 
#endif    
c        Set N2O5 to equilibrium when necessary (near ground,
c        N2O5 is thermally unstable, has a very short lifetime)
         if(igas.eq.n_N2O5.and.(-dest(igas,L).ge.y(n_N2O5,L)*0.75d0.or.
     &    chemrate(iN2O5form,L).gt.y(n_NOx,L)))then
           rnewval=(rr(iN2O5form,L)*y(nNO3,L)*y(nNO2,L))/
     &     (rr(iN2O5decomp,L)*y(nM,L)+ss(7,L,I,J)+1.d-12)
           if(rnewval.lt.1.d0)rnewval=1.d0
           changeL(L,igas)=(rnewval-y(n_N2O5,L))
           if(changeL(L,igas).gt.0.33d0*y(nNO2,L))changeL(L,igas)=
     &      0.33d0*y(nNO2,L)
           changeL(L,igas)=changeL(L,igas)*conc2mass
         endif
c
c        Conserve NOx with respect to N2O5
         if(igas.eq.n_NOx.and.(-dest(n_N2O5,L).ge.y(n_N2O5,L).or.
     &    chemrate(iN2O5form,L).gt.y(n_NOx,L)))then
          rnewval=(rr(iN2O5form,L)*y(nNO3,L)*y(nNO2,L))/
     &    (rr(iN2O5decomp,L)*y(nM,L)+ss(7,L,I,J)+1.d-12)
          if(rnewval.lt.1.d0)rnewval=1.d0
           changeX=(rnewval-y(n_N2O5,L))
           if(changeX.gt.0.33d0*y(nNO2,L))changeX=0.33d0*y(nNO2,L)
          changeL(L,igas)=
     &    changeL(L,igas)-changeX*conc2mass
         endif
c
c        Set HO2NO2 to equil when necessary
         if(igas.eq.n_HO2NO2.and.(-dest(igas,L).ge.y(n_HO2NO2,L).or.
     &    chemrate(iHO2NO2form,L).gt.y(n_NOx,L)))then
          rnewval=(rr(iHO2NO2form,L)*y(nHO2,L)*y(nNO2,L))/
     *    (rr(iHO2NO2_OH,L)*y(nOH,L)+rr(iHO2NO2decomp,L)*y(nM,L)
     *    +ss(10,L,I,J)+ss(11,L,I,J)+1.d-12)
          if(rnewval.lt.1.d0)rnewval=1.d0
           changeL(L,igas)=(rnewval-y(n_HO2NO2,L))
           if(changeL(L,igas).gt.0.33d0*y(nNO2,L))changeL(L,igas)=
     *      0.33d0*y(nNO2,L)
           changeL(L,igas)=changeL(L,igas)*conc2mass
         endif
c
c        Conserve NOx with respect to HO2NO2
         if(igas.eq.n_NOx.and.(-dest(n_HO2NO2,L).ge.y(n_HO2NO2,L).or.
     &    chemrate(iHO2NO2form,L).gt.y(n_NOx,L)))then
          rnewval=(rr(iHO2NO2form,L)*y(nHO2,L)*y(nNO2,L))/
     *    (rr(iHO2NO2_OH,L)*y(nOH,L)+rr(iHO2NO2decomp,L)*y(nM,L)
     *    +ss(10,L,I,J)+ss(11,L,I,J)+1.d-12)
          if(rnewval.lt.1.d0)rnewval=1.d0
          changeX=(rnewval-y(n_HO2NO2,L))
          if(changeX.gt.0.33d0*y(nNO2,L))changeX=0.33d0*y(nNO2,L)
          changeL(L,igas)=changeL(L,igas)-
     *     changeX*conc2mass
         endif
c
c        Set PAN to equilibrium when necessary (near ground,
c         PAN is thermally unstable, has a very short lifetime)
         if(igas.eq.n_PAN.and.(-dest(igas,L).ge.y(n_PAN,L).or.
     &    chemrate(iPANform,L).gt.y(n_NOx,L)))then
          rnewval=(rr(iPANform,L)*y(nC2O3,L)*y(nNO2,L))/
     *    (rr(iPANdecomp,L)*y(nM,L)+ss(15,L,I,J)+1.d-12)
          if(rnewval.lt.1.d0)rnewval=1.d0
          changeL(L,igas)=(rnewval-y(n_PAN,L))
          if(changeL(L,igas).gt.0.33d0*y(nNO2,L))changeL(L,igas)=
     &     0.33d0*y(nNO2,L)
          changeL(L,igas)=changeL(L,igas)*conc2mass
         endif
c
c        Conserve NOx with respect to PAN
         if(igas.eq.n_NOx.and.(-dest(n_PAN,L).ge.y(n_PAN,L).or.
     &    chemrate(iPANform,L).gt.y(n_NOx,L)))then
          rnewval=(rr(iPANform,L)*y(nC2O3,L)*y(nNO2,L))/
     *    (rr(iPANdecomp,L)*y(nM,L)+ss(15,L,I,J)+1.d-12)
          if(rnewval.lt.1.d0)rnewval=1.d0
          changeX=(rnewval-y(n_PAN,L))
          if(changeX.gt.0.33d0*y(nNO2,L))changeX=0.33d0*y(nNO2,L)
          changeL(L,igas)=changeL(L,igas)-changeX*
     *     conc2mass
         endif
c
#ifdef SHINDELL_STRAT_CHEM
c        Set HOBr to equilibrium when necessary
         if(igas.eq.n_HOBr.and.(-dest(igas,L).ge.y(n_HOBr,L).or.
     &    chemrate(73,L).gt.0.5d0*y(n_BrOx,L)))then
          rnewval=(rr(73,L)*y(nBrO,L)*y(nHO2,L))/
     *     (ss(24,L,i,j)+1.d-12)
          if(rnewval.lt.1.d0)rnewval=1.d0
          changeL(L,igas)=(rnewval-y(n_HOBr,L))
          if(changeL(L,igas).gt.0.5d0*y(nBrO,L))changeL(L,igas)=
     *     0.5d0*y(nBrO,L)
          changeL(L,igas)=changeL(L,igas)*conc2mass
         endif
c
c        Conserve BrOx with respect to HOBr
         if(igas.eq.n_BrOx.and.(-dest(n_HOBr,L).ge.y(n_HOBr,L).or.
     &    chemrate(73,L).gt.0.5d0*y(n_BrOx,L)))then
          rnewval=(rr(73,L)*y(nBrO,L)*y(nHO2,L))/
     *     (ss(24,L,i,j)+1.d-12)
          if(rnewval.lt.1.d0)rnewval=1.d0
          changeX=(rnewval-y(n_HOBr,L))
          if(changeX.gt.0.5d0*y(nBrO,L))changeX=0.5d0*y(nBrO,L)
          changeL(L,igas)=changeL(L,igas)-
     *     changeX*conc2mass
         endif
c
c        Set BrONO2 to equilibrium when necessary
         if(igas.eq.n_BrONO2.and.(-dest(igas,L).ge.y(n_BrONO2,L).or.
     &    chemrate(104,L).gt.0.5d0*y(n_BrOx,L)))then
          rnewval=(rr(104,L)*y(nBrO,L)*y(nNO2,L))/
     *     (ss(23,L,i,j)+1.d-12)
          if(rnewval.lt.1.d0)rnewval=1.d0
          changeL(L,igas)=(rnewval-y(n_BrONO2,L))
          if(changeL(L,igas).gt.0.5d0*y(nBrO,L))changeL(L,igas)=
     *     0.5d0*y(nBrO,L)
          changeL(L,igas)=changeL(L,igas)*conc2mass
         endif
c
c        Conserve BrOx with respect to BrONO2
         if(igas.eq.n_BrOx.and.(-dest(n_BrONO2,L).ge.y(n_BrONO2,L).or.
     &    chemrate(104,L).gt.0.5d0*y(n_BrOx,L)))then
          rnewval=(rr(104,L)*y(nBrO,L)*y(nNO2,L))/
     *     (ss(23,L,i,j)+1.d-12)
          if(rnewval.lt.1.d0)rnewval=1.d0
          changeX=(rnewval-y(n_BrONO2,L))
          if(changeX.gt.0.5d0*y(nBrO,L))changeX=0.5d0*y(nBrO,L)
          changeL(L,igas)=changeL(L,igas)-changeX*
     *     conc2mass
         endif
c
c        Conserve NOx with respect to BrONO2
         if(igas.eq.n_NOx.and.(-dest(n_BrONO2,L).ge.y(n_BrONO2,L).or.
     &    chemrate(104,L).gt.0.5d0*y(n_BrOx,L)))then
          rnewval=(rr(104,L)*y(nBrO,L)*y(nNO2,L))/
     *     (ss(23,L,i,j)+1.d-12)
          if(rnewval.lt.1.d0)rnewval=1.d0
          changeX=(rnewval-y(n_BrONO2,L))
          if(changeX.gt.0.5d0*y(nBrO,L))changeX=0.5d0*y(nBrO,L)
          changeL(L,igas)=changeL(L,igas)-changeX*
     *     conc2mass
         endif
c
c        Set ClONO2 to equilibrium when necessary
         if(igas.eq.n_ClONO2.and.(-dest(igas,L).ge.y(n_ClONO2,L).or.
     &    chemrate(103,L).gt.0.8d0*y(n_ClOx,L)))then
          rnewval=(rr(103,L)*y(nClO,L)*y(nNO2,L))/
     *     (ss(22,L,i,j)+rr(65,L)*y(nO,L)+1.d-12)
          if(rnewval.lt.1.d0)rnewval=1.d0
          changeL(L,igas)=(rnewval-y(n_ClONO2,L))
          if(changeL(L,igas).gt.0.3d0*y(nClO,L))changeL(L,igas)=
     *     0.3d0*y(nClO,L)
          if(-changeL(L,igas).gt.0.8d0*y(n_ClONO2,L))
     *     changeL(L,igas)=-0.8d0*y(n_ClONO2,L)
          changeL(L,igas)=changeL(L,igas)*conc2mass
         endif
c
c        Conserve ClOx with respect to ClONO2
         if(igas.eq.n_ClOx.and.(-dest(n_ClONO2,L).ge.y(n_ClONO2,L).or.
     &    chemrate(103,L).gt.0.8d0*y(n_ClOx,L)))then
          rnewval=(rr(103,L)*y(nClO,L)*y(nNO2,L))/
     *     (ss(22,L,i,j)+rr(65,L)*y(nO,L)+1.d-12)
          if(rnewval.lt.1.d0)rnewval=1.d0
          changeX=(rnewval-y(n_ClONO2,L))
          if(changeX.gt.0.3d0*y(nClO,L))changeX=0.3d0*y(nClO,L)
          if(-changeX.gt.0.8d0*y(n_ClONO2,L))changeX=
     &    -0.8d0*y(n_ClONO2,L)
          changeL(L,igas)=changeL(L,igas)-changeCl*
     *     conc2mass
         endif
c
c        Conserve NOx with respect to ClONO2
         if(igas.eq.n_NOx.and.(-dest(n_ClONO2,L).ge.y(n_ClONO2,L).or.
     &    chemrate(103,L).gt.0.8d0*y(n_ClOx,L)))then
          rnewval=(rr(103,L)*y(nClO,L)*y(nNO2,L))/
     *     (ss(22,L,i,j)+rr(65,L)*y(nO,L)+1.d-12)
          if(rnewval.lt.1.d0)rnewval=1.d0
          changeX=(rnewval-y(n_ClONO2,L))
          if(changeX.gt.0.3d0*y(nClO,L))changeX=0.3d0*y(nClO,L)
          if(-changeX.gt.0.8d0*y(n_ClONO2,L))changeX=
     &     -0.8d0*y(n_ClONO2,L)
          changeL(L,igas)=changeL(L,igas)-changeCl*
     *     conc2mass
         endif
c
c        Set HOCl to equilibrium when necessary
         if(igas.eq.n_HOCl.and.(-dest(igas,L).ge.y(n_HOCl,L).or.
     &    chemrate(63,L).gt.y(n_ClOx,L)))then
          rnewval=(rr(63,L)*y(nClO,L)*y(nHO2,L))/
     *     (ss(21,L,i,j)+rr(55,L)*y(nO,L)+1.d-12)
          if(rnewval.lt.1.d0)rnewval=1.d0
          changeL(L,igas)=(rnewval-y(n_HOCl,L))
          if(changeL(L,igas).gt.0.3d0*y(nClO,L))changeL(L,igas)=
     *     0.3d0*y(nClO,L)
          changeL(L,igas)=changeL(L,igas)*conc2mass
         endif
c
c        Conserve ClOx with respect to HOCl
         if(igas.eq.n_ClOx.and.(-dest(n_HOCl,L).ge.y(n_HOCl,L).or.
     &    chemrate(63,L).gt.y(n_ClOx,L)))then
          rnewval=(rr(63,L)*y(nClO,L)*y(nHO2,L))/
     *     (ss(21,L,i,j)+rr(55,L)*y(nO,L)+1.d-12)
          if(rnewval.lt.1.d0)rnewval=1.d0
          changeX=(rnewval-y(n_HOCl,L))
          if(changeX.gt.0.3d0*y(nClO,L))changeX=0.3d0*y(nClO,L)
          changeL(L,igas)=changeL(L,igas)-changeX*
     *     conc2mass
         endif
#endif
c
       end do ! L
      end do  ! igas
c
#ifdef SHINDELL_STRAT_CHEM
c     separate N2O change for N cons, leave out N2O->N2+O fm cons
      do L=1,maxl
      sv_changeN2O(L)=-chemrate(87,L)*dxyp(J)*rMAbyM(L)/mass2vol(n_N2O)
      enddo
#endif
c
cc    Ensure nitrogen conservation
c     (since equilibration of short lived gases may alter this)
      if(prnchg.and.J.eq.jprn.and.I.eq.iprn)then
       write(*,*) 'changes (mass) before nitrogen conservation routine'
       write(*,*) 'NOx, N2O5, HO2NO2, HNO3, PAN, AlkylNit, N2O'
#ifdef SHINDELL_STRAT_CHEM
       write(*,*) 'ClONO2, BrONO2'
#endif
       write(*,*) changeL(lprn,n_NOx),changeL(lprn,n_N2O5),
     & changeL(lprn,n_HO2NO2),changeL(lprn,n_HNO3),
     & changeL(lprn,n_PAN),changeL(lprn,n_AlkylNit)
#ifdef SHINDELL_STRAT_CHEM
     & ,changeL(lprn,n_N2O)
     & ,changeL(lprn,n_ClONO2),changeL(lprn,n_BrONO2)
       write(*,*) 'N2O change w/o rxns forming N2',sv_changeN2O(lprn)
#endif
      endif
c
      do L=1,maxl
c
c       First check for nitrogen loss > 100%
        if(-changeL(L,n_NOx).gt.trm(I,J,L,n_NOx))
     &  changeL(L,n_NOx)=1.d0-trm(I,J,L,n_NOx)
        if(-changeL(L,n_N2O5).gt.trm(I,J,L,n_N2O5))
     &  changeL(L,n_N2O5)=1.d0-trm(I,J,L,n_N2O5)
        if(-changeL(L,n_HO2NO2).gt.trm(I,J,L,n_HO2NO2))
     &  changeL(L,n_HO2NO2)=1.d0-trm(I,J,L,n_HO2NO2)
        if(-changeL(L,n_HNO3).gt.trm(I,J,L,n_HNO3))
     &  changeL(L,n_HNO3)=1.d0-trm(I,J,L,n_HNO3)
        if(-changeL(L,n_PAN).gt.trm(I,J,L,n_PAN))
     &  changeL(L,n_PAN)=1.d0-trm(I,J,L,n_PAN)
        if(-changeL(L,n_AlkylNit).gt.trm(I,J,L,n_AlkylNit))
     &  changeL(L,n_AlkylNit)=1.d0-trm(I,J,L,n_AlkylNit)
#ifdef SHINDELL_STRAT_CHEM
        if(-changeL(L,n_ClONO2).gt.trm(I,J,L,n_ClONO2))
     &  changeL(L,n_ClONO2)=1.d0-trm(I,J,L,n_ClONO2)
        if(-changeL(L,n_BrONO2).gt.trm(I,J,L,n_BrONO2))
     &  changeL(L,n_BrONO2)=1.d0-trm(I,J,L,n_BrONO2)
#endif
c
c       Next insure balance between dNOx and sum of dOthers
        sumN=(2.d0*changeL(L,n_N2O5))*mass2vol(n_N2O5)+
     *  (changeL(L,n_HNO3))*mass2vol(n_HNO3)+
     *  (changeL(L,n_HO2NO2))*mass2vol(n_HO2NO2)+
     *  (changeL(L,n_PAN))*mass2vol(n_PAN)+
     *  (changeL(L,n_AlkylNit))*mass2vol(n_AlkylNit)
#ifdef SHINDELL_STRAT_CHEM
          if(L.ge.maxT+1)sumN=sumN+
     *     changeL(L,n_ClONO2)*mass2vol(n_ClONO2)+
     *     changeL(L,n_BrONO2)*mass2vol(n_BrONO2)
         dNOx=changeL(L,n_NOx)*mass2vol(n_NOx)+
     *    2*sv_changeN2O(L)*mass2vol(n_N2O)
#else
        dNOx=changeL(L,n_NOx)*mass2vol(n_NOx)
#endif
        if(prnchg.and.J.eq.jprn.and.I.eq.iprn.and.L.eq.lprn)
     &       write(*,*)
     *       'other N changes, dNOx (less prod fm N2O) = (molec) ',
     *       sumN,dNOx
        ratio=-sumN/dNOx
        IF(ratio.le.0.999d0 .OR. ratio.ge.1.001d0) THEN
         if(dNOx.gt.0.d0)then
c         NOx being produced (net positive change)
          if (ratio.gt.1.d0)then
           sumD=0.d0
c          reduce N destruction to match NOx prodcution
           if(changeL(L,n_N2O5).lt.0.d0)   sumD=sumD+
     *     2.d0*changeL(L,n_N2O5)*mass2vol(n_N2O5)
           if(changeL(L,n_HO2NO2).lt.0.d0) sumD=sumD+
     *     changeL(L,n_HO2NO2)*mass2vol(n_HO2NO2)
           if(changeL(L,n_HNO3).lt.0.d0)   sumD=sumD+
     *     changeL(L,n_HNO3)*mass2vol(n_HNO3)
           if(changeL(L,n_PAN).lt.0.d0)    sumD=sumD+
     *     changeL(L,n_PAN)*mass2vol(n_PAN)
           if(changeL(L,n_AlkylNit).lt.0.d0)sumD=sumD+
     *     changeL(L,n_AlkylNit)*mass2vol(n_AlkylNit)
#ifdef SHINDELL_STRAT_CHEM
           if(changeL(L,n_ClONO2).lt.0.d0)sumD=sumD+
     *      changeL(L,n_ClONO2)*mass2vol(n_ClONO2)
           if(changeL(L,n_BrONO2).lt.0.d0)sumD=sumD+
     *      changeL(L,n_BrONO2)*mass2vol(n_BrONO2)
#endif
           newD=(sumN/ratio)+sumD-sumN
           ratioD=newD/sumD
           if(changeL(L,n_N2O5).lt.0.d0)    changeL(L,n_N2O5)=
     *     changeL(L,n_N2O5)    *ratioD
           if(changeL(L,n_HO2NO2).lt.0.d0)  changeL(L,n_HO2NO2)=
     *     changeL(L,n_HO2NO2)  *ratioD
           if(changeL(L,n_HNO3).lt.0.d0)    changeL(L,n_HNO3)=
     *     changeL(L,n_HNO3)    *ratioD
           if(changeL(L,n_PAN).lt.0.d0)     changeL(L,n_PAN)=
     *     changeL(L,n_PAN)     *ratioD
           if(changeL(L,n_AlkylNit).lt.0.d0)changeL(L,n_AlkylNit)=
     *     changeL(L,n_AlkylNit)*ratioD
#ifdef SHINDELL_STRAT_CHEM
           vClONO2=changeL(L,n_ClONO2)*(1.d0-ratioD)
           if(changeL(L,n_ClONO2).lt.0.d0)changeL(L,n_ClONO2)=
     *      changeL(L,n_ClONO2)*ratioD
           changeL(L,n_ClOx)=changeL(L,n_ClOx)+vClONO2*
     *      (mass2vol(n_ClONO2)/mass2vol(n_ClOx)) !ensure Cl cons
           vBrONO2=changeL(L,n_BrONO2)*(1.d0-ratioD)
           if(changeL(L,n_BrONO2).lt.0.d0)changeL(L,n_BrONO2)=
     *      changeL(L,n_BrONO2)*ratioD
           changeL(L,n_BrOx)=changeL(L,n_BrOx)+vBrONO2*
     *      (mass2vol(n_BrONO2)/mass2vol(n_BrOx)) !ensure Br cons
#endif
          endif
c
          if (ratio.le.1.d0.and.ratio.gt.0.d0)then
c          reduce NOx production to match N loss
           changeL(L,n_NOx)=changeL(L,n_NOx)*ratio
#ifdef SHINDELL_STRAT_CHEM
           changeL(L,n_NOx)=changeL(L,n_NOx)-
     *      2.d0*sv_changeN2O(L)*mass2vol(n_N2O)/mass2vol(n_NOx)
#endif
          endif
         else
c         NOx being destroyed (net change is negative)
          if (ratio.gt.1.d0)then
           sumP=0
c          reduce N production to match NOx loss
           if(changeL(L,n_N2O5).gt.0.d0)    sumP=sumP+
     *     2.d0*changeL(L,n_N2O5)*mass2vol(n_N2O5)
           if(changeL(L,n_HO2NO2).gt.0.d0)  sumP=sumP+
     *     changeL(L,n_HO2NO2)*mass2vol(n_HO2NO2)
           if(changeL(L,n_HNO3).gt.0.d0)    sumP=sumP+
     *     changeL(L,n_HNO3)*mass2vol(n_HNO3)
           if(changeL(L,n_PAN).gt.0.d0)     sumP=sumP+
     *     changeL(L,n_PAN)*mass2vol(n_PAN)
           if(changeL(L,n_AlkylNit).gt.0.d0)sumP=sumP+
     *     changeL(L,n_AlkylNit)*mass2vol(n_AlkylNit)
#ifdef SHINDELL_STRAT_CHEM
           if(changeL(L,n_ClONO2).gt.0.d0)sumP=sumP+
     *      changeL(L,n_ClONO2)*mass2vol(n_ClONO2)
           if(changeL(L,n_BrONO2).gt.0.d0)sumP=sumP+
     *      changeL(L,n_BrONO2)*mass2vol(n_BrONO2)
#endif
           newP=(sumN/ratio)+sumP-sumN
           ratioP=newP/sumP
           if(changeL(L,n_N2O5).gt.0.d0)     changeL(L,n_N2O5)=
     *        changeL(L,n_N2O5)*ratioP
           if(changeL(L,n_HO2NO2).gt.0.d0)  changeL(L,n_HO2NO2)=
     *        changeL(L,n_HO2NO2)*ratioP
           if(changeL(L,n_HNO3).gt.0.d0)    changeL(L,n_HNO3)=
     *        changeL(L,n_HNO3)*ratioP
           if(changeL(L,n_PAN).gt.0.d0)     changeL(L,n_PAN)=
     *        changeL(L,n_PAN)*ratioP
           if(changeL(L,n_AlkylNit).gt.0.d0)changeL(L,n_AlkylNit)=
     *        changeL(L,n_AlkylNit)*ratioP
#ifdef SHINDELL_STRAT_CHEM
           vClONO2=changeL(L,n_ClONO2)*(1.d0-ratioP)
           if(changeL(L,n_ClONO2).gt.0.d0)changeL(L,n_ClONO2)=
     *      changeL(L,n_ClONO2)*ratioP
           changeL(L,n_ClOx)=changeL(L,n_ClOx)+vClONO2*
     *      (mass2vol(n_ClONO2)/mass2vol(n_ClOx)) !ensure Cl cons
           vBrONO2=changeL(L,n_BrONO2)*(1.d0-ratioP)
           if(changeL(L,n_BrONO2).gt.0.d0)changeL(L,n_BrONO2)=
     *      changeL(L,n_BrONO2)*ratioP
           changeL(L,n_BrOx)=changeL(L,n_BrOx)+vBrONO2*
     *      (mass2vol(n_BrONO2)/mass2vol(n_BrOx)) !ensure Br cons
#endif
          endif
          if (ratio.le.1.d0.and.ratio.gt.0.d0)then
c          increase N production to match NOx loss (12/4/2001)
           sumP=0.d0
           if(changeL(L,n_N2O5).gt.0.d0)    sumP=sumP+
     *     2.d0*changeL(L,n_N2O5)*mass2vol(n_N2O5)
           if(changeL(L,n_HO2NO2).gt.0.d0)  sumP=sumP+
     *     changeL(L,n_HO2NO2)*mass2vol(n_HO2NO2)
           if(changeL(L,n_HNO3).gt.0.d0)    sumP=sumP+
     *     changeL(L,n_HNO3)*mass2vol(n_HNO3)
           if(changeL(L,n_PAN).gt.0.d0)     sumP=sumP+
     *     changeL(L,n_PAN)*mass2vol(n_PAN)
           if(changeL(L,n_AlkylNit).gt.0.d0)sumP=sumP+
     *     changeL(L,n_AlkylNit)*mass2vol(n_AlkylNit)
           newP=(sumN/ratio)+sumP-sumN
           ratioP=newP/sumP
           if(changeL(L,n_N2O5).gt.0.d0)    changeL(L,n_N2O5)=
     *     changeL(L,n_N2O5)*ratioP
           if(changeL(L,n_HO2NO2).gt.0.d0)  changeL(L,n_HO2NO2)=
     *     changeL(L,n_HO2NO2)*ratioP
           if(changeL(L,n_HNO3).gt.0.d0)    changeL(L,n_HNO3)=
     *     changeL(L,n_HNO3)*ratioP
           if(changeL(L,n_PAN).gt.0.d0)     changeL(L,n_PAN)=
     *     changeL(L,n_PAN)*ratioP
           if(changeL(L,n_AlkylNit).gt.0.d0)changeL(L,n_AlkylNit)=
     *     changeL(L,n_AlkylNit)*ratioP
          endif
         endif

        END IF ! skipped section above if ratio very close to one
C
        if(prnchg.and.J.eq.jprn.and.I.eq.iprn.and.L.eq.lprn)
     &  write(*,*) 'ratio for conservation =',ratio
C
      end do ! L
c
#ifdef SHINDELL_STRAT_CHEM
c     Remove some of the HNO3 formed heterogeneously, as it doesn't come
c     back to the gas phase
      do L=maxT+1,LM
       changeL(L,n_HNO3)=changeL(L,n_HNO3)-
     *  0.005d0*rr(105,L)*y(n_HNO3,L)*dt2
     *  *(dxyp(J)*rMAbyM(L))/(mass2vol(n_HNO3))
      enddo
#endif
C
cc    Print chemical changes in a particular grid box if desired
      if(prnchg.and.J.eq.jprn.and.I.eq.iprn)then

        do igas=1,nlast
         changeA=changeL(Lprn,igas)*y(nM,lprn)*mass2vol(igas)*
     *   bydxyp(J)*byam(lprn,I,J)
         if(y(igas,lprn).eq.0.d0)then
            write(6,156) ay(igas),': ',changeA,' molecules;  y=0'
         else
            write(6,155) ay(igas),': ',changeA
     *           ,' molecules produced; ',
     *   (100.d0*changeA)/y(igas,lprn),' percent of'
     *   ,y(igas,lprn),'(',1.d9*y(igas,lprn)/y(nM,lprn),' ppbv)'
         endif
c
         if(igas.eq.nlast)then
#ifdef SHINDELL_STRAT_CHEM
         if(LPRN.GE.maxT+1)then
            write(6,155) ay(nH2O),': ',
     *      changeH2O(lprn),' molecules produced; ',
     *      (100*changeH2O(lprn))/y(nH2O,lprn),' percent of',
     *      y(nH2O,lprn),'(',1.d6*y(nH2O,lprn)/y(nM,lprn),' ppmv)'
         else
          write(6,'(a10,58x,e13.3,6x,f10.3,a5)')
     *    ' H2O     :',y(nH2O,LPRN),(y(nH2O,LPRN)/
     *    y(nM,LPRN))*1.d6,' ppmv'
         endif
#else
          write(6,'(a10,58x,e13.3,6x,f10.3,a5)')
     *    ' H2O     :',y(nH2O,LPRN),(y(nH2O,LPRN)/
     *    y(nM,LPRN))*1.d6,' ppmv'
#endif
          write(6,'(a10,58x,e13.3,6x,f10.3,a5)')
     *    ' CH3O2   :',yCH3O2(I,J,LPRN),(yCH3O2(I,J,LPRN)/
     *    y(nM,LPRN))*1.d9,' ppbv'
          write(6,'(a10,58x,e13.3,6x,f10.3,a5)')
     *    ' C2O3    :',y(nC2O3,LPRN),(y(nC2O3,LPRN)/y(nM,LPRN))*1.d9,
     *    ' ppbv'
          write(6,'(a10,58x,e13.3,6x,f10.3,a5)')
     *    ' XO2     :',y(nXO2,LPRN),(y(nXO2,LPRN)/y(nM,LPRN))*1.d9,
     *    ' ppbv'
          write(6,'(a10,58x,e13.3,6x,f10.3,a5)')
     *    ' XO2N    :',y(nXO2N,LPRN),(y(nXO2N,LPRN)/
     *    y(nM,LPRN))*1.d9,' ppbv'
          write(6,'(a10,58x,e13.3,6x,f10.3,a5)')
     *    ' RXPAR   :',y(nRXPAR,LPRN),(y(nRXPAR,LPRN)/
     *    y(nM,LPRN))*1.d9,' ppbv'
          write(6,'(a10,58x,e13.3,6x,f10.3,a5)')
     *    ' Aldehyde:',y(nAldehyde,LPRN),(y(nAldehyde,LPRN)/
     *    y(nM,LPRN))*1.d9,' ppbv'
          write(6,'(a10,58x,e13.3,6x,f10.3,a5)')
     *    ' ROR     :',y(nROR,LPRN),(y(nROR,LPRN)/
     *    y(nM,LPRN))*1.d9,' ppbv'
         endif
        enddo
      endif  !end of prnchg loop
c
c*** tracer masses & slopes are now updated in apply_tracer_3Dsource ***
      do igas=1,ntm_chem
       do L=1,maxl
#ifdef regional_Ox_tracers
c*** calculate chemical changes for regional Ox tracers ***
        if(igas.gt.nlast) then
          nREG=igas-nlast
          changeL(L,igas)=trm(I,J,L,igas)*Oxloss(L)/trm(I,J,L,n_Ox)
C ----- in region criterion -----
          if(lat_dg(J,1).ge.regOx_s(nREG).and.(lat_dg(J,1).lt.
     &    regOx_n(nREG).or.(lat_dg(J,1).eq.regOx_n(nREG).and.J.eq.JM))
     &                        .and.
     &    PRES(L).le.regOx_b(nREG).and.(PRES(L).gt.regOx_t(nREG)
     &    .or.(PRES(L).eq.regOx_t(nREG).and.L.eq.maxl))
     &                        .and.
     &    lon_dg(I,1).ge.regOx_w(nREG).and.(lon_dg(I,1).lt.
     &    regOx_e(nREG).or.(lon_dg(I,1).eq.regOx_e(nREG).and.I.eq.IM)))
C -------------------------------
     &    changeL(L,igas)=changeL(L,igas)+Oxprod(L)    
        end if
#endif
c*** limit the change due to chemistry ***
        if(changeL(L,igas).gt.1.d20) then
           WRITE(99,*)'change set to 0 in chemstep: I,J,L,igas,change'
     &    ,I,J,L,igas,changeL(L,igas)
           changeL(L,igas) = 0.d0
        endif
        if(-changeL(L,igas).gt.trm(I,J,L,igas)) THEN
          if(prnchg)
     &    WRITE(99,*)'change .gt. mass, so use 95%: I,J,L,igas,change'
     &    ,I,J,L,igas,changeL(L,igas)
          changeL(L,igas) = -0.95d0*trm(I,J,L,igas)
        endif
C surface Ox change diagnostic:
C       if(L.eq.1.and.igas.eq.n_Ox.and.changeL(L,igas).le.1.d20)then
C          changeA=(changeL(L,igas)*y(nM,L)*mass2vol(igas))*
C    *      bydxyp(J)*byam(L,I,J)
C          TAIJS(I,J,ijs_OxL1)=TAIJS(I,J,ijs_OxL1)+1.d9*changeA/y(nM,L)
C       end if

#ifdef  SHINDELL_STRAT_CHEM
c     No bromine chemistry or HOCl chemistry in troposphere
         if(L.le.maxT.and.igas.eq.n_BrOx.or.igas.eq.n_BrONO2)
     *    changeL(L,igas)=0.d0
         if(L.le.maxT.and.igas.eq.n_HBr.or.igas.eq.n_HOBr)
     *    changeL(L,igas)=0.d0
         if(PRES(L).lt.300.d0.and.igas.eq.n_HOCl)changeL(L,igas)=0.d0
#endif
       end do    ! L
      end do     ! igas

#ifdef  SHINDELL_STRAT_CHEM
c     No ozone chemistry (for now) at mesospheric levels
      do L=1,maxl
        if(PRES(L).lt.2.d-1) changeL(L,n_Ox)=0.d0
      end do
#endif

C**** special diags not associated with a particular tracer
      DO L=1,maxl
         if (y(nOH,L).gt.0.d0 .and. y(nOH,L).lt.1.d20) 
     &   TAJLS(J,L,jls_OHcon)=TAJLS(J,L,jls_OHcon)+y(nOH,L)
         TAJLS(J,L,jls_H2Omr)=TAJLS(J,L,jls_H2Omr)+(y(nH2O,L)/y(nM,L))
      END DO
      TAJLS(J,1,jls_day)=TAJLS(J,1,jls_day)+1.d0 
CCCC  need to save 3D OH here too...

 155  format(1x,a8,a2,e13.3,a21,f10.0,a11,2x,e13.3,3x,a1,f12.5,a6)
 156  format(1x,a8,a2,e13.3,a16)
c
      return
      end SUBROUTINE chemstep
cc    __________________________________________________________________
      SUBROUTINE rates(maxl,I,J)
!@sum rates calculate reaction rates with present concentrations
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on chemcalc0C5.4_M23p)
c
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: nr,chemrate,photrate,rr,y,nn,dt2,
     &                          ss,ks,ny,dest,prod,JPPJ,nhet
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var kalt local dummy L-loop variable
!@var maxl passed highest chemistry level
!@var ireac,igas dummy loop variables
!@var I,J passed horizontal spatial indicies
      INTEGER kalt, ireac, igas, maxl
      INTEGER, INTENT(IN) :: I,J
c
      do kalt=1,maxl
        do ireac=1,nr-nhet       ! non-heterogeneous
          chemrate(ireac,kalt)=rr(ireac,kalt)*y(nn(1,ireac),kalt)*
     &    y(nn(2,ireac),kalt)*dt2
        enddo
        do ireac=nr-nhet+1,nr    ! heterogeneous
          chemrate(ireac,kalt)=rr(ireac,kalt)*y(nn(1,ireac),kalt)*dt2
        enddo
        do ireac=1,JPPJ          ! photolysis
          photrate(ireac,kalt)=ss(ireac,kalt,I,J)*y(ks(ireac),kalt)*dt2
        enddo
c       Initialize change arrays
        do igas=1,ny
          dest(igas,kalt)=0.d0
          prod(igas,kalt)=0.d0
        enddo
      enddo
      return
      end SUBROUTINE rates

cc    __________________________________________________________________
      SUBROUTINE chem1(kdnr,maxl,numel,nn,ndnr,chemrate,dest,multip)
!@sum chem1 calculate chemical destruction/production
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on chemcalc0C5.4_M23p)
c
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: p_2, p_3, p_4, ny, numfam,nfam
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var maxl passed highest chemistry level
!@var numel first index of nn array (shouldn't code like this anymore?)
!@var kdnr kdnr,kpnr,kds, or kps    passed from chemstep
!@var nn nn,nnr,ks, or kss          passed from chemstep
!@var ndnr ndnr,npnr,nds, or nps    passed from chemstep
!@var chemrate chemrate or photrate passed from chemstep
!@var dest dest or prod             passed from chemstep
!@var multip -1 for destruction, +1 for production
!@var i,ireac,igas,ial,nbeg,nend dummy loop variables
!@var dk dummy variable
      INTEGER ireac,igas,ial,i,dk,nbeg,nend
      INTEGER, INTENT(IN)            :: maxl, numel,multip
      INTEGER, DIMENSION(p_4)        :: kdnr
      INTEGER, DIMENSION(numel,p_2)  :: nn       ! can't do anymore?
      INTEGER, DIMENSION(p_3)        :: ndnr
      REAL*8,  DIMENSION(p_2,maxl)   :: chemrate ! can't do anymore?
      REAL*8,  DIMENSION(ny,maxl)    :: dest
C
      ireac=0
c     Reactive families
      do igas=1,numfam
        dk=kdnr(igas+1)-kdnr(igas)
        if(dk.ge.1) then
          do i=1,dk
           ireac=ireac+1
           do ial=1,maxl
            if(nn(1,ndnr(ireac)).ge.nfam(igas).and.nn(1,ndnr(ireac)).lt.
     *       nfam(igas+1))then
             dest(igas,ial)=dest(igas,ial)+multip
     &       *chemrate(ndnr(ireac),ial)
c            Save change array for individual family elements
             dest(nn(1,ndnr(ireac)),ial)=dest(nn(1,ndnr(ireac)),ial)+
     *          multip*chemrate(ndnr(ireac),ial)
            endif
            if(numel.eq.2)then
             if(nn(2,ndnr(ireac)).ge.nfam(igas).and.nn(2,ndnr(ireac))
     *        .lt.nfam(igas+1))then
              dest(igas,ial)=dest(igas,ial)+
     *        multip*chemrate(ndnr(ireac),ial)
              dest(nn(2,ndnr(ireac)),ial)=dest(nn(2,ndnr(ireac)),ial)+
     *        multip*chemrate(ndnr(ireac),ial)
             endif
            endif
           end do ! ial
          end do  ! i
        end if
      end do      ! igas
C
c     INDIVIDUAL SPECIES
      nbeg=numfam+1
      nend=nfam(1)-1
      do igas=nbeg,nend
        dk=kdnr(igas+1)-kdnr(igas)
        if(dk.ge.1) then
          do i=1,dk
           ireac=ireac+1
           do ial=1,maxl
            dest(igas,ial)=dest(igas,ial)+multip*
     &      chemrate(ndnr(ireac),ial)
           end do
          end do
        end if
      end do
      return
      end SUBROUTINE chem1

cc    __________________________________________________________________
      SUBROUTINE chem1prn(kdnr,numel,nn,ndnr,chemrate,
     &                    index,multip,igas,total,maxl,I,J)
!@sum chem1prn for printing out the chemical reactions
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on chemcalc0C5.4_M23p)
c
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: ay, lprn, nfam, p_4, numfam, y,
     &                              p_2, p_3
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var kdnr kdnr,kpnr,kds, or kps from chemstep
!@var numel first index of nn array (shouldn't code like this anymore?)
!@var nn nn,nnr,ks, or kss from chemstep
!@var ndnr ndnr,npnr,nds, or nps from chemstep
!@var chemrate chemrate or photrate from chemstep
!@var index passed index to know which call this is... {1,2,3,4}
!@var multip 1 for production, -1 for destruction
!@var igas passed index for gas number
!@var total dummy summation
!@var maxl highest chemistry level
!@var I,J passed horizontal spatial indicies
!@var label character string for printing
!@var irec dummy loop variables
!@var per dummy temp variable
      INTEGER, INTENT(IN) :: igas,I,J,maxl,multip,index,numel
      INTEGER ireac
      REAL*8 total,per
      REAL*8,  DIMENSION(p_2,maxl)   :: chemrate !can't do this anymore?
      INTEGER, DIMENSION(p_3)        :: ndnr
      INTEGER, DIMENSION(numel,p_2)  :: nn       !can't do this anymore?
      INTEGER, DIMENSION(p_4)        :: kdnr
      character*17 label

c     skip this section during (within) family chemisty calls
      if(igas.le.numfam) then
c
c     FAMILIES
      if(kdnr(igas+1)-kdnr(igas).lt.1)goto40
      do ireac=kdnr(igas),kdnr(igas+1)-1
         if(index.le.2)label=' chem reaction # '
         if(index.gt.2)label=' phot reaction # '
         if(nn(1,ndnr(ireac)).ge.nfam(igas).and.nn(1,ndnr(ireac)).lt.
     *     nfam(igas+1))then
           per=0.d0
           if(y(igas,lprn).ne.0.d0)per=multip*100.d0*
     *     chemrate(ndnr(ireac),lprn)/y(igas,lprn)
           write(6,177) label,ndnr(ireac),' percent change from ',
     *     ay(nn(1,ndnr(ireac))),' = ',per,
     *     ' dy=',multip*chemrate(ndnr(ireac),lprn)
           total=total+per
         endif
         if(numel.eq.2)then
          if(nn(2,ndnr(ireac)).ge.nfam(igas).and.nn(2,ndnr(ireac)).lt.
     *     nfam(igas+1))then
           per=0.d0
           if(y(igas,lprn).ne.0.d0)per=multip*100.d0*
     *     chemrate(ndnr(ireac),lprn)/y(igas,lprn)
           write(6,177) label,ndnr(ireac),' percent change from ',
     *     ay(nn(2,ndnr(ireac))),' = ',per,
     *     ' dy=',multip*chemrate(ndnr(ireac),lprn)
           total=total+per
          endif
         endif  ! end of numel=2 loop
      end do  !end of ireac loop
      goto40

      end if ! non family chem only loop

c     INDIVIDUAL SPECIES
      if(kdnr(igas+1)-kdnr(igas).lt.1)goto40
      do ireac=kdnr(igas),kdnr(igas+1)-1
        if(index.le.2)label=' chem reaction # '
        if(index.gt.2)label=' phot reaction # '
c       skip same reaction if written twice
        if (ireac.gt.1) then
          if (ndnr(ireac).eq.ndnr(ireac-1)) goto 30 
        end if
         if(nn(1,ndnr(ireac)).eq.igas)then
           per=0.d0
           if(y(igas,lprn).ne.0.d0)per=100.d0*multip*
     *     chemrate(ndnr(ireac),lprn)/y(igas,lprn)
           write(6,106) label,ndnr(ireac),' percent change = ',per,
     *     ' dy=',multip*chemrate(ndnr(ireac),lprn)
           total=total+per
         endif
         if(numel.eq.2)then
          if(nn(2,ndnr(ireac)).eq.igas)then
           per=0.d0
           if(y(igas,lprn).ne.0.d0)per=100.d0*multip*
     *     chemrate(ndnr(ireac),lprn)/y(igas,lprn)
           write(6,106) label,ndnr(ireac),' percent change = ',per,
     *     ' dy=',multip*chemrate(ndnr(ireac),lprn)
           total=total+per
          endif
         endif  !end numel=2 loop
 106     format(a17,i3,a18,f10.0,a4,e12.3)
 177     format(a17,i3,a21,a8,a3,f10.0,a4,e12.3)
  30  end do ! ireac
  40  continue
      return
      end SUBROUTINE chem1prn
