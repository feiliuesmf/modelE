#include "rundeck_opts.h"
      SUBROUTINE chemstep(I,J,changeL,ierr_loc)
!@sum chemstep Calculate new concentrations after photolysis & chemistry
!@vers 2013/03/27
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@calls rates,chem1,chem1prn
c
C**** GLOBAL parameters and variables:
C
      USE SOMTQ_COM, only       : qmom
      USE RAD_COM, only         : clim_interact_chem
      USE RESOLUTION, only      : ls1
      USE RESOLUTION, only      : im,jm,lm
      USE ATM_COM, only         : Q
      USE DOMAIN_DECOMP_ATM,only : grid,getDomainBounds,write_parallel
      USE ATM_COM, only         : MA, byMA,ltropo
      USE GEOM, only            : byaxyp,axyp
      USE TRDIAG_COM, only : taijls=>taijls_loc,jls_OHcon,jls_day
     &     ,jls_OxpT,jls_OxdT,jls_Oxp,jls_Oxd,jls_COp,jls_COd,ijlt_OH
     &     ,ijlt_HO2,ijlt_COp,ijlt_COd,ijlt_Oxd,ijlt_Oxp,ijlt_CH4d
     &     ,ijlt_OxpRO2
     &     ,jls_ClOcon,jls_H2Ocon,jls_H2Ochem
      use OldTracer_mod, only: vol2mass, mass2vol
      USE TRACER_COM, only: n_CH4,n_CH3OOH,n_Paraffin,n_PAN,n_Isoprene,
     &                   n_stratOx
#ifdef TRACERS_TERP
      USE TRACER_COM, only: n_Terpenes
#endif  /* TRACERS_TERP */
      USE TRACER_COM, only: n_AlkylNit,n_Alkenes,n_N2O5,n_NOx,n_HO2NO2
#ifdef TRACERS_AEROSOLS_SOA
      USE TRACER_COM, only: n_isopp1g,n_isopp1a,n_isopp2g,n_isopp2a
#ifdef TRACERS_TERP
      USE TRACER_COM, only: n_apinp1g,n_apinp1a,n_apinp2g,n_apinp2a
#endif  /* TRACERS_TERP */
#endif  /* TRACERS_AEROSOLS_SOA */
      USE TRACER_COM, only: n_Ox,n_HNO3,n_H2O2,n_CO,n_HCHO,trm,NTM,
     &                  n_N2O,n_ClOx,n_BrOx,n_HCl,n_HOCl,n_ClONO2,n_HBr,
     &                  n_HOBr,n_BrONO2,n_CFC,ntm_chem
#ifdef TRACERS_WATER
      use OldTracer_mod, only: tr_wd_type, nWater, tr_H2ObyCH4
      USE TRACER_COM, only: trm, trmom 
#endif
#ifdef TRACERS_HETCHEM
      USE TRACER_COM, only: krate,n_N_d1,n_N_d2,n_N_d3
#endif
      USE TRCHEM_Shindell_COM, only: chemrate,photrate,cpd,
     &                   yCH3O2,yC2O3,yXO2,yXO2N,yRXPAR,yAldehyde,
     &                   yROR,nCH3O2,nC2O3,nXO2,nXO2N,nRXPAR,
     &                   nAldehyde,nROR,nr,nn,dt2,nss,ks,dest,prod,
     &                   ny,rr,nO1D,nOH,nNO,nHO2,ta,nM,ss,
     &                   nO3,nNO2,nNO3,prnrts,jprn,iprn,lprn,ay,
     &                   prnchg,y,kss,nps,kps,nds,kds,
     &                   npnr,nnr,ndnr,kpnr,kdnr,nH2O,which_trop,
     &                   Jacet,acetone
     &                   ,SF3,ratioNs,ratioN2,rNO2frac,nO,nClO,nBrO
     &                   ,rNOfrac,rNOdenom,nOClO,nCl,nBr
     &                   ,nCl2,yCl2,SF2,nO2,MWabyMWw,yCl2O2,pscX
#ifdef TRACERS_AEROSOLS_SOA
       USE TRACERS_SOA, only: apartmolar,whichsoa,soa_apart,LM_soa
#endif  /* TRACERS_AEROSOLS_SOA */
      USE DIAG_COM, only : aj,j_h2och4,ftype,ntype
      USE ATM_COM, only : pmidl00
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var changeL change due to chemistry (kg)
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
!@var dxbym2v is axyp over mass2volume
!@var sv_changeN2O N2O change without portion making N2 (for N cons)
!@var vClONO2, vBrONO2 temporary vars within N conservation
!@var changeH2O chemical change in H2O
!@var Oxcorr account for Ox change from within NOx partitioning
!@var rNO3prod,rNO2prod,rNOprod to acct for dOx from NOx partitioning
!@var PRES local nominal pressure for regional Ox tracers
      INTEGER, INTENT(IN) :: I,J
      INTEGER, INTENT(INOUT) :: ierr_loc
      INTEGER :: L,iter,maxl,igas,maxT,Lz,it,n
      INTEGER :: J_0, J_1
#ifdef TRACERS_TERP
      INTEGER, PARAMETER :: iHO2NO2form=102,iN2O5form=103,
     &iPANform=105,iHO2NO2_OH=18,iHO2NO2decomp=95,iN2O5decomp=96
     &,iPANdecomp=29,iClOplusNO2=107,iBrOplusNO2=108,iClOplusClO=106
     &,iOHplusNO2=101,iNOplusO=99
     &,iTerpenesOH=92,iTerpenesO3=93
#else
      INTEGER, PARAMETER :: iHO2NO2form=99,iN2O5form=100,
     &iPANform=102,iHO2NO2_OH=18,iHO2NO2decomp=92,iN2O5decomp=93
     &,iPANdecomp=29,iClOplusNO2=104,iBrOplusNO2=105,iClOplusClO=103
     &,iOHplusNO2=98,iNOplusO=96
#endif  /* TRACERS_TERP */
      character(len=300) :: out_line
      logical            :: jay
      REAL*8, DIMENSION(LM,ntm) :: changeL
      REAL*8, DIMENSION(LM) :: rMAbyM,sv_changeN2O,changeH2O,Oxcorr,
     & PRES,dQ,dQM,fraQ2,c2ml,conOH,conClO,conH2O,
     &     NprodOx_pos,NprodOx_neg
      REAL*8 qqqCH3O2,CH3O2loss,XO2_NO,XO2N_HO2,RXPAR_PAR,ROR_CH2,
     & C2O3prod,C2O3dest,XO2prod,XO2dest,XO2_XO2,XO2Nprod,XO2Ndest,
     & RXPARprod,RXPARdest,Aldehydeprod,Aldehydedest,RORprod,RORdest,
     & total,rnewval,dNOx,ratio,sumD,newD,ratioD,newP,ratioP,
     & changeA,sumP,tempiter,tempiter2,sumC,sumN,sumH,sumB,sumO,sumA,
     & dxbym2v,changeX,vClONO2,vBrONO2,conc2mass,rNO3prod,rNO2prod,
     & rNOprod,changeAldehyde,rxnN2,rxnN3,rxnN4,NprodOx,NlossNOx,byta,
     & diffCH3O2,tempAcet,prodCH3O2,dQMsum

      call getDomainBounds(grid, J_STRT    =J_0,  J_STOP    =J_1)
      
      jay = (J >= J_0 .and. J <= J_1) 
     
C Pick top level for chemistry, and tropopause layer:  
      select case(which_trop)
      case(0); maxT=ltropo(I,J)
      case(1); maxT=ls1-1
      case default; call stop_model('which_trop problem 1',255)
      end select 
      maxl=LM

      PRES(1:maxL)=PMIDL00(1:maxL)   !SIG(1:maxL)*(PSF-PTOP)+PTOP
      
      do L=1,maxT
        y(nCH3O2,L)   =    yCH3O2(I,J,L)
        y(nC2O3,L)    =     yC2O3(I,J,L)
        y(nXO2,L)     =      yXO2(I,J,L)
        y(nXO2N,L)    =     yXO2N(I,J,L)
        y(nRXPAR,L)   =    yRXPAR(I,J,L)
        y(nAldehyde,L)= yAldehyde(I,J,L)
        y(nROR,L)     =      yROR(I,J,L)
      enddo
      do L=maxT+1,LM
       y(nCH3O2,L)=yCH3O2(I,J,L)
      enddo
C
C Calculate reaction rates with present concentrations:
      call rates(maxl,I,J)
c
c chem1 call sample:
c (klist,l,numel,nlist,ndlist,rate,change,multip)
c numel=number of elements in reaction list nlist (1 or 2)
c change=dest or prod array
c multip=1(prod) or -1(dest)

c chemical destruction:
      call chem1(kdnr,maxl,2,nn,ndnr,chemrate,dest,-1)
c chemical production:
      call chem1(kpnr,maxl,2,nnr,npnr,chemrate,prod,1)
c photolytic destruction:
      call chem1(kds,maxl,1,ks,nds,photrate,dest,-1)
c photolytic production:
      call chem1(kps,maxl,2,kss,nps,photrate,prod,1)

c Add additional Cl from CFC photolysis + background :
      do L=1,maxl
        prod(n_ClOx,L)=prod(n_ClOx,L)+0.33d0*photrate(26,L)+
     &  7.5d-3*photrate(28,L)
        prod(n_BrOx,L)=prod(n_BrOx,L)+5.55d-4*photrate(26,L)+
     &  5.2d-6*photrate(28,L)
      enddo

c Oxidation of Isoprene and Alkenes produces less than one
c HCHO, Alkenes, and CO per rxn, correct here following Houweling:
      do L=1,maxl
        prod(n_CO,L)=prod(n_CO,L)-0.63d0*chemrate(35,L)
        prod(n_HCHO,L)=prod(n_HCHO,L)-0.36d0*chemrate(35,L)
        prod(n_HCHO,L)=prod(n_HCHO,L)-0.39d0*chemrate(30,L)
#ifdef TRACERS_TERP
     &                               -0.39d0*chemrate(iTerpenesOH,L)
#endif  /* TRACERS_TERP */
        prod(n_Alkenes,L)=prod(n_Alkenes,L)-0.42d0*chemrate(30,L)
#ifdef TRACERS_TERP
     &                               -0.42d0*chemrate(iTerpenesOH,L)
#endif  /* TRACERS_TERP */
        prod(n_HCHO,L)=prod(n_HCHO,L)-0.10d0*chemrate(31,L)
#ifdef TRACERS_TERP
     &                               -0.10d0*chemrate(iTerpenesO3,L)
#endif  /* TRACERS_TERP */
        prod(n_Alkenes,L)=prod(n_Alkenes,L)-0.45d0*chemrate(31,L)
#ifdef TRACERS_TERP
     &                               -0.45d0*chemrate(iTerpenesO3,L)
#endif  /* TRACERS_TERP */
#ifdef TRACERS_HETCHEM
        dest(n_HNO3,l)=dest(n_HNO3,l)-krate(i,j,l,1,1)*y(n_HNO3,l)*dt2
#endif
      enddo
#ifdef TRACERS_AEROSOLS_SOA
      call soa_apart ! calculate current apartmolar factors
#ifdef SOA_DIAGS
     &              (I,J)
#endif  /* SOA_DIAGS */
      do L=1,LM_soa
        prod(n_isopp1g,L)=prod(n_isopp1g,L)+
     &                    apartmolar(L,whichsoa(n_isopp1a))*
     &                    (chemrate(30,L)+chemrate(31,L))
        prod(n_isopp2g,L)=prod(n_isopp2g,L)+
     &                    apartmolar(L,whichsoa(n_isopp2a))*
     &                    (chemrate(30,L)+chemrate(31,L))
#ifdef TRACERS_TERP
        prod(n_apinp1g,L)=prod(n_apinp1g,L)+
     &                    apartmolar(L,whichsoa(n_apinp1a))*
     &                    chemrate(iTerpenesO3,L)
        prod(n_apinp2g,L)=prod(n_apinp2g,L)+
     &                    apartmolar(L,whichsoa(n_apinp2a))*
     &                    chemrate(iTerpenesO3,L)
#endif  /* TRACERS_TERP */
      enddo
#endif  /* TRACERS_AEROSOLS_SOA */

c Set CH3O2 values (concentration = production/specific loss):
      do L=1,maxT
        iter=1
        qqqCH3O2=(rr(11,L)*y(nO1D,L)+rr(12,L)*y(nOH,L))
     &  *y(n_CH4,L)+rr(23,L)*y(n_CH3OOH,L)*y(nOH,L)
        tempAcet=2.d0*Jacet(L)*acetone(I,J,L)
        prodCH3O2=qqqCH3O2+tempAcet
        tempiter=rr(20,L)*y(nNO,L)+rr(22,L)*y(nHO2,L)
        do while(iter <= 7)
          CH3O2loss=tempiter+rr(27,L)*yCH3O2(I,J,L)
          if(CH3O2loss > 1.d-7)then
            y(nCH3O2,L)=prodCH3O2/CH3O2loss
          else
            y(nCH3O2,L)=1.d0
          endif
          iter=iter+1
        end do
        
c Conserve carbon wrt CH3O2 changes:
        diffCH3O2=y(nCH3O2,L)-yCH3O2(I,J,L)
        if(diffCH3O2 > tempAcet)then
c         reduce non-acetone source gases (CH4 and CH3OOH):
          dest(n_CH4,L)=dest(n_CH4,L)-(diffCH3O2-tempAcet)
     &    *(qqqCH3O2-rr(23,L)*y(n_CH3OOH,L)*y(nOH,L))/qqqCH3O2
          dest(n_CH3OOH,L)=dest(n_CH3OOH,L)-(diffCH3O2-tempAcet)
     &    *(rr(23,L)*y(n_CH3OOH,L)*y(nOH,L))/qqqCH3O2
        else if(diffCH3O2 < tempAcet)then
c         increase non-acetone product gases:
          prod(n_HCHO,L)=prod(n_HCHO,L)-(diffCH3O2-tempAcet)
     &    *(CH3O2loss-rr(22,L)*y(nHO2,L))/CH3O2loss
          prod(n_CH3OOH,L)=prod(n_CH3OOH,L)-(diffCH3O2-tempAcet)
     &    *(rr(22,L)*y(nHO2,L))/CH3O2loss
        end if
        yCH3O2(I,J,L)=y(nCH3O2,L)
      enddo
      
      do L=maxT+1,maxl
        iter=1
        qqqCH3O2=(rr(11,L)*y(nO1D,L)+rr(12,L)*y(nOH,L))
     &  *y(n_CH4,L)+rr(23,L)*y(n_CH3OOH,L)*y(nOH,L)
     &  +rr(82,l)*y(nCl,L)
        tempiter=rr(20,L)*y(nNO,L)+rr(22,L)*y(nHO2,L)
     &  +rr(85,l)*y(nClO,l)
        do while (iter <= 7)
          CH3O2loss=tempiter+rr(27,L)*yCH3O2(I,J,L)
          if(CH3O2loss > 1.d-7)then
            y(nCH3O2,L)=qqqCH3O2/CH3O2loss
          else
            y(nCH3O2,L)=1.d-5
          endif
          iter=iter+1
        end do            
            
c Conserve carbon wrt CH3O2 changes:
        diffCH3O2=y(nCH3O2,L)-yCH3O2(I,J,L)
        if(diffCH3O2 > 0.d0)then
c         reduce source gases (CH4 and CH3OOH):
          dest(n_CH4,L)=dest(n_CH4,l)-diffCH3O2
     &    *(qqqCH3O2-rr(23,L)*y(n_CH3OOH,L)*y(nOH,L))/qqqCH3O2
          dest(n_CH3OOH,L)=dest(n_CH3OOH,l)-diffCH3O2
     &    *(rr(23,L)*y(n_CH3OOH,L)*y(nOH,L))/qqqCH3O2
        else if(diffCH3O2 < 0.d0)then
c         increase product gases:
          prod(n_HCHO,l)=prod(n_HCHO,l)-diffCH3O2
     &    *(CH3O2loss-rr(22,L)*y(nHO2,L))/CH3O2loss
          prod(n_CH3OOH,l)=prod(n_CH3OOH,l)-diffCH3O2
     &    *(rr(22,L)*y(nHO2,L))/CH3O2loss
        endif
        yCH3O2(I,J,L)=y(nCH3O2,L)
      enddo

      do L=1,maxT ! ---------- troposphere loop ---------

c Set C2O3, XO2, XO2N, RXPAR, Aldehyde & ROR values:

c        First set various specific loss rates:
         XO2_NO=y(nNO,L)*4.2d-12*exp(180.d0/ta(L))
         XO2N_HO2=y(nHO2,L)*y(nNO,L)*rr(44,L)*
     &   rr(43,L)/XO2_NO
         RXPAR_PAR=y(n_Paraffin,L)*8.d-11
         ROR_CH2=1.6d3

c       Set value for C2O3:
        iter=1
        C2O3prod=rr(38,L)*yAldehyde(I,J,L)*y(nOH,L)+
     &  (rr(29,L)*y(nM,L)+ss(15,L,I,J))*y(n_PAN,L)+
     &  0.15d0*rr(31,L)*y(nO3,L)*y(n_Isoprene,L)
#ifdef TRACERS_TERP
     & +0.15d0*rr(iTerpenesO3,L)*y(nO3,L)*y(n_Terpenes,L)
#endif  /* TRACERS_TERP */
        tempiter=rr(39,L)*y(nNO,L)+rr(iPANform,L)*y(nNO2,L)+
     &  rr(41,L)*y(nHO2,L)
        do while (iter <= 7)
          C2O3dest=tempiter+rr(40,L)*yC2O3(I,J,L)
          if(C2O3dest > 1.d-7)then
            y(nC2O3,L)=(C2O3prod/C2O3dest)
          else
            y(nC2O3,L)=1.d0
          endif
          yC2O3(I,J,L)=y(nC2O3,L)
          iter=iter+1
        end do

c       Set value for XO2:
! remember to update voc2nox if you update any of the following XO2 loss reactions
        iter=1
        XO2prod=ss(16,L,I,J)*yAldehyde(I,J,L)+
     &  y(nC2O3,L)*(rr(39,L)*y(nNO2,L)+rr(40,L)*
     &  y(nC2O3,L)*2.d0+rr(41,L)*y(nHO2,L))
     &  +rr(42,L)*yROR(I,J,L)*0.96d0
     &  +y(nOH,L)*(rr(37,L)*y(n_Paraffin,L)*0.87d0+rr(34,L)*
     &  y(n_Alkenes,L)+rr(30,L)*y(n_Isoprene,L)*0.85d0+
#ifdef TRACERS_TERP
     &  rr(iTerpenesOH,L)*y(n_Terpenes,L)*0.85d0+
#endif  /* TRACERS_TERP */
     &  rr(33,L)*y(n_AlkylNit,L))+
     &  y(nO3,L)*(rr(35,L)*y(n_Alkenes,L)*0.29d0+
     &  rr(31,L)*y(n_Isoprene,L)*0.18d0
#ifdef TRACERS_TERP
     &  +rr(iTerpenesO3,L)*y(n_Terpenes,L)*0.18d0
#endif  /* TRACERS_TERP */
     &           )
        tempiter=XO2_NO+rr(43,L)*y(nHO2,L)
        tempiter2=1.7d-14*exp(1300.d0/ta(L))
        do while (iter <= 7)
          XO2_XO2=yXO2(I,J,L)*tempiter2
          XO2dest=tempiter+XO2_XO2
          if(XO2dest > 1.d-7.and.ss(16,L,I,J) > 1.d-6)then
            y(nXO2,L)=(XO2prod/XO2dest)
          else
            y(nXO2,L)=1.d0
          endif
          yXO2(I,J,L)=y(nXO2,L)
          iter=iter+1
        end do

c       Set value for XO2N:
        XO2Nprod=rr(37,L)*y(n_Paraffin,L)*y(nOH,L)*0.13d0+
     &  rr(42,L)*yROR(I,J,L)*0.04d0+rr(30,L)*y(n_Isoprene,L)*
     &  y(nOH,L)*0.15d0
#ifdef TRACERS_TERP
     & +rr(iTerpenesOH,L)*y(n_Terpenes,L)*y(nOH,L)*0.15d0
#endif  /* TRACERS_TERP */
        XO2Ndest=XO2N_HO2+rr(44,L)*y(nNO,L)
        if(XO2Ndest > 1.d-7)then
          y(nXO2N,L)=(XO2Nprod/XO2Ndest)
        else
          y(nXO2N,L)=1.d0
        endif
        yXO2N(I,J,L)=y(nXO2N,L)

#ifdef ACCMIP_LIKE_DIAGS
        TAIJLS(I,J,L,ijlt_OxpRO2)=TAIJLS(I,J,L,ijlt_OxpRO2)+
     &  (y(nXO2,L)*XO2_NO+y(nXO2N,L)*y(nNO,L)*rr(44,L))*cpd
#endif

c       Set value for RXPAR:
        RXPARprod=rr(37,L)*y(n_Paraffin,L)*y(nOH,L)*0.11d0+
     &  rr(34,L)*yROR(I,J,L)*2.1d0+rr(35,L)*y(n_Alkenes,L)*
     &  y(nO3,L)*0.9d0
        RXPARdest=RXPAR_PAR
        if(RXPARdest > 0.d0)then
          y(nRXPAR,L)=(RXPARprod/RXPARdest)
        else
          y(nRXPAR,L)=1.d0
        endif
        yRXPAR(I,J,L)=y(nRXPAR,L)

c       Set value for Aldehyde:
        Aldehydeprod=rr(37,L)*y(n_Paraffin,L)*y(nOH,L)*0.11d0+
     &  rr(34,L)*y(n_Alkenes,L)*y(nOH,L)+
     &  rr(42,L)*yROR(I,J,L)*1.1d0+rr(35,L)*y(n_Alkenes,L)*
     &  y(nO3,L)*0.44d0
        Aldehydedest=rr(38,L)*y(nOH,L)+ss(16,L,I,J)
c       Check for equilibrium:
        if(Aldehydedest*y(nAldehyde,L)*dt2 < y(nAldehyde,L))then
          changeAldehyde=
     &    (Aldehydeprod-y(nAldehyde,L)*Aldehydedest)*dt2
          if(changeAldehyde > y(nAldehyde,L))
     &    changeAldehyde=y(nAldehyde,L)
          y(nAldehyde,L)=y(nAldehyde,L)+changeAldehyde
          if(y(nAldehyde,L) < 0.d0) y(nAldehyde,L)=1.d0
        else
          y(nAldehyde,L)=(Aldehydeprod/(Aldehydedest+0.5d-5))
        endif
        yAldehyde(I,J,L)=y(nAldehyde,L)

c       Set value for ROR:
        RORprod=rr(37,L)*y(n_Paraffin,L)*y(nOH,L)*0.76d0+
     &  rr(42,L)*yROR(I,J,L)*0.02d0
        RORdest=rr(42,L)+ROR_CH2
        if(RORdest > 0.d0)then
          y(nROR,L)=(RORprod/RORdest)
        else
          y(nROR,L)=1.d0
        endif
        yROR(I,J,L)=y(nROR,L)

c       Add parrafin loss term via rxpar reaction and
c       prod term via isoprene rxns:
        dest(n_Paraffin,L)=dest(n_Paraffin,L)-y(nRXPAR,L)*RXPAR_PAR*dt2
        prod(n_Paraffin,L)=prod(n_Paraffin,L)+0.63d0*y(n_Isoprene,L)
     &  *(rr(30,L)*y(nOH,L)+rr(31,L)*y(nO3,L))*dt2
#ifdef TRACERS_TERP
     &  +5.0d0*0.63d0*y(n_Terpenes,L)
     &  *(rr(iTerpenesOH,L)*y(nOH,L)+rr(iTerpenesO3,L)*y(nO3,L))*dt2
#endif  /* TRACERS_TERP */

c       Add CH3OOH production via XO2N + HO2:
        prod(n_CH3OOH,L)=prod(n_CH3OOH,L)+XO2N_HO2*y(nXO2N,L)*dt2
c
      end do  ! --------------------------------------

c If NOx in equil with N2O5, HO2NO2, or PAN, remove from changes:
      do L=1,maxl
        if(-dest(n_N2O5,L) >= y(n_N2O5,L) .or.
     &  chemrate(iN2O5form,L) > y(n_NOx,L))then
          dest(n_NOx,L)=dest(n_NOx,L)+2.d0*chemrate(iN2O5form,L)
          prod(n_NOx,L)=prod(n_NOx,L)-2.d0*(chemrate(iN2O5decomp,L)
     &    +photrate(7,L))
        endif
        if(-dest(n_HO2NO2,L) >= y(n_HO2NO2,L) .or.
     &  chemrate(iHO2NO2form,L) > y(n_NOx,L))then
          dest(n_NOx,L)=dest(n_NOx,L)+chemrate(iHO2NO2form,L)
          prod(n_NOx,L)=prod(n_NOx,L)-(chemrate(iHO2NO2_OH,L)+
     &    chemrate(iHO2NO2decomp,L)+photrate(10,L)+photrate(11,L))
        endif
        if(-dest(n_PAN,L) >= y(n_PAN,L) .or.
     &  chemrate(iPANform,L) > y(n_NOx,L))then
          dest(n_NOx,L)=dest(n_NOx,L)+chemrate(iPANform,L)
          prod(n_NOx,L)=prod(n_NOx,L)-(chemrate(iPANdecomp,L)+
     &    photrate(15,L))
        endif
        
c If BrOx in equil with HOBr or BrONO2, remove from changes:
        if(-dest(n_HOBr,L) >= y(n_HOBr,L).or.
     &  chemrate(73,L) > 0.5d0*y(n_BrOx,L))then
          dest(n_BrOx,L)=dest(n_BrOx,L)+chemrate(73,L)
          prod(n_BrOx,L)=prod(n_BrOx,L)-photrate(24,L)
        endif
        if(-dest(n_BrONO2,L) >= y(n_BrONO2,L).or.
     &  chemrate(iBrOplusNO2,L) > 0.5d0*y(n_BrOx,L))then
          dest(n_BrOx,L)=dest(n_BrOx,L)+chemrate(iBrOplusNO2,L)
          prod(n_BrOx,L)=prod(n_BrOx,L)-photrate(23,L)
          dest(n_NOx,L)=dest(n_NOx,L)+chemrate(iBrOplusNO2,L)
          prod(n_NOx,L)=prod(n_NOx,L)-photrate(23,L)
        endif
        
c If ClOx in equil with HOCl or ClONO2, remove from changes:
        if(-dest(n_HOCl,L) >= y(n_HOCl,L) .or.
     &  chemrate(63,L) > y(n_ClOx,L))then
          dest(n_ClOx,L)=dest(n_ClOx,L)+chemrate(63,L)
          prod(n_ClOx,L)=prod(n_ClOx,L)-(photrate(21,L)+chemrate(55,L))
        endif
        if(-dest(n_ClONO2,L) >= y(n_ClONO2,L) .or.
     &  chemrate(iClOplusNO2,L) > 0.8d0*y(n_ClOx,L))then
          dest(n_ClOx,L)=dest(n_ClOx,L)+chemrate(iClOplusNO2,L)
          prod(n_ClOx,L)=prod(n_ClOx,L)-(photrate(22,L)+chemrate(65,L))
          dest(n_NOx,L)=dest(n_NOx,L)+chemrate(iClOplusNO2,L)
          prod(n_NOx,L)=prod(n_NOx,L)-(photrate(22,L)+chemrate(65,L))
        endif
      enddo

c Calculate water vapor change AND APPLY TO MODEL Q VARIABLE:
      do L=1,maxl !! for a long time, this was: maxT+1,maxl
        changeH2O(L)=(2.d0*y(n_CH4,L)*
     *  (rr(11,L)*y(nO1D,L)+rr(12,L)*y(nOH,L)+rr(82,L)*y(nCl,L))
     *  -2.0d0*SF3(I,J,L)*y(nH2O,L))*dt2  
C       And apply that change here and accumulate a diagnostic:
C       --- y --- :
        y(nH2O,L)=y(nH2O,L)+changeH2O(L)
C       --- Q --- :
        dQ(L) = changeH2O(L)/(y(nM,L)*MWabyMWw)
        dQM(L) = dQ(L)*MA(L,I,J)*axyp(I,J)
        if(clim_interact_chem > 0)then
          fraQ2(l)=(Q(I,J,L)+changeH2O(L)/(y(nM,L)*MWabyMWw))/Q(I,J,L)
          Q(I,J,L) = Q(I,J,L) + dQ(L)
C       -- Qmom --:
          if(changeH2O(L) < 0.)then
            qmom(:,i,j,l)=qmom(:,i,j,l)*fraQ2(l)
            if(fraQ2(l) <= 0.98)then
              write(out_line,*)'> 2% Q change in calc IJL,change='
     &        ,I,J,L,fraQ2(l)
              call write_parallel(trim(out_line),crit=.true.)
              call stop_model('big Q change in calc',255)
            endif
          endif
        endif
      enddo

C     -- diags --:
      call inc_tajls2_column(i,j,1,maxl,lm,jls_H2Ochem,dQM)
      if(clim_interact_chem > 0)then
        dQMsum = sum(dQM(1:maxl))/axyp(i,j)
        do it=1,ntype
          call inc_aj(i,j,it,j_h2och4,dQMsum*ftype(it,i,j))
        enddo
#ifdef TRACERS_WATER
C     -- water tracers --:
        do n=1,ntm
          select case (tr_wd_type(n))
          case (nWater)           ! water: add CH4-sourced water to tracers
            do l=1,maxl
              trm(i,j,l,n) = trm(i,j,l,n) + tr_H2ObyCH4(n)*dQM(l)
              if(changeH2O(L) < 0.) trmom(:,i,j,l,n) = trmom(:,i,j,l,n)
     *             *fraQ2(l)
            enddo
          end select
        end do
#endif
      endif

c Calculate ozone change due to within-NOx partitioning:
      do L=1,LM
        if(y(nO1D,L) == 0.) CYCLE
c       account for NO2 and NO ozone destruction:
        rNO2prod=rr(18,L)*y(nOH,L)*y(n_HO2NO2,L)+
     &  rr(iHO2NO2decomp,L)*y(n_HO2NO2,L)+ss(9,L,I,J)*y(n_HNO3,L)+
     &  ss(10,L,I,J)*y(n_HO2NO2,L)+ss(23,L,I,J)*y(n_BrONO2,L)
        rNOprod=rr(87,L)*y(n_N2O,L)*y(nO1D,L)
        rNO3prod=rr(65,L)*y(nO,L)*y(n_ClONO2,L)+
     &  ss(7,L,I,J)*y(n_N2O5,L)+ss(11,L,I,J)*y(n_HO2NO2,L)+
     &  ss(22,L,I,J)*y(n_ClONO2,L)
c       add production of NO and NO2 from NO3:
        rNO3prod=rNO3prod*ss(6,L,I,J)/(ss(5,L,I,J)+ss(6,L,I,J)+1.d0)
        rNO2prod=rNO2prod+rNO3prod
        rNOprod=rNOprod+rNO3prod
        ratioNs=rNO2prod/rNOprod
        ratioN2=y(nNO2,L)/y(nNO,L)
        
        if(ratioNs > ratioN2)then !excess NO2 production
        
c         account for NO2 that then goes via NO2+O->NO+O2, NO2->NO+O:
          rNO2frac=(rr(26,L)*y(nO,L)-ss(1,L,I,J))/
     &    (rr(iOHplusNO2,L)*y(nOH,L)+
     &    rr(iHO2NO2form,L)*y(nHO2,L)+rr(iN2O5form,L)*y(nNO3,L)+
     &    rr(iClOplusNO2,L)*y(nClO,L)+rr(iBrOplusNO2,L)*y(nBrO,L)+
     &    rr(26,L)*y(nO,L)+ss(1,L,I,J))
          Oxcorr(L)=(rNO2prod-rNOprod)*rNO2frac*dt2*y(nNO,L)/y(n_NOx,L)
          if(Oxcorr(L) > -1.d18 .and. Oxcorr(L) < 1.d18)then
            dest(n_Ox,L)=dest(n_Ox,L)-Oxcorr(L)
          else
            ierr_loc=ierr_loc+1 ! will stop model in masterchem
            write(out_line,'(a17,5(1X,E11.4))')
     &      'Oxcorr fault NO2:',
     &      ratioNs,ratioN2,rNO2frac,rNO2prod,rNOprod
            call write_parallel(trim(out_line),crit=.true.)      
            return
          endif

        else                      !excess NO prodcution

c         account for NO that then goes via NO+O3->NO2+O2
c         or NO+O+M->NO2+M:
          rNOfrac=(rr(5,L)*y(nO3,L)+rr(iNOplusO,L)*y(nO,L))
          rNOdenom=(rr(5,L)*y(nO3,L)+rr(iNOplusO,L)*y(nO,L)+
     &    rr(6,L)*y(nHO2,L)+rr(44,L)*y(nXO2N,L)+1.d0)+
     &    rr(20,L)*yCH3O2(I,J,L)+
     &    rr(39,L)*y(nC2O3,L)+4.2d-12*exp(180/ta(L))*y(nXO2,L)+
     &    rr(64,L)*y(nClO,L)+
     &    rr(67,L)*y(nOClO,L)+rr(71,L)*y(nBrO,L)

          rNOfrac=rNOfrac/rNOdenom
          Oxcorr(L)=(rNOprod-rNO2prod)*rNOfrac*dt2*y(nNO2,L)/y(n_NOx,L)
          if(Oxcorr(L) > -1.d18 .and. Oxcorr(L) < 1.d18)then
            dest(n_Ox,L)=dest(n_Ox,L)-Oxcorr(L)
          else
            ierr_loc=ierr_loc+1 ! will stop model in masterchem
            write(out_line,'(a16,3I4,10(1X,E11.4))')'Oxcorr fault NO:',
     &      I,J,L,ratioNs,ratioN2,rNOfrac,rNO2prod,rNOprod,y(nNO2,L),
     &      y(nNO,L),rNOdenom,y(nO,L),y(nO3,L)
            call write_parallel(trim(out_line),crit=.true.)    
            return
          endif
c
        endif
      enddo ! 1->LM loop

c Calculate ozone change due to Cl2O2 cycling:
      do L=1,LM
        if(yCl2O2(I,J,L) > 1d1)dest(n_Ox,L)=dest(n_Ox,L) - 0.75d0*
     &  rr(45,L)*y(nCl,L)*y(nO3,L)*dt2*yCl2O2(I,J,L)*1.5d9/y(nM,L)
      enddo

c Include oxidation of CO by O(1D)
      do L=1,LM
        dest(n_CO,L)=dest(n_CO,L)-1.0d-9*y(n_CO,L)*y(nO1D,L)*dt2
      end do

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c           Print chemistry diagnostics if desired :
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C             REACTION RATES, CHEMICAL CHANGES
c (chem1prn: argument before multip is index = number of call):
      
      if(prnrts .and. J==jprn .and. I==iprn)then
        do igas=1,ntm_chem
          total=0.d0
          write(out_line,108)' Species: ',ay(igas)
          call write_parallel(trim(out_line),crit=jay)

          call chem1prn
     &    (kdnr,2,nn,ndnr,chemrate,1,-1,igas,total,maxl,I,J,jay)

          if(igas == n_NOx)then
            if(-dest(n_HO2NO2,lprn) >= y(n_HO2NO2,lprn) .or.
     &      chemrate(iHO2NO2form,lprn) > y(n_NOx,lprn)) then
              write(out_line,110)
     &        'loss by reaction    (HO2NO2 formation) removed',
     &        chemrate(iHO2NO2form,lprn)
              call write_parallel(trim(out_line),crit=jay)
            endif
            if(-dest(n_N2O5,lprn) >= y(n_N2O5,lprn) .or.
     &      chemrate(iN2O5form,lprn) > y(n_NOx,lprn)) then     
              write(out_line,110)
     &        'losses by reaction    (N2O5 formation) removed',
     &        2.d0*chemrate(iN2O5form,lprn)
              call write_parallel(trim(out_line),crit=jay)
            endif
            if(-dest(n_PAN,lprn) >= y(n_PAN,lprn) .or.
     &      chemrate(iPANform,lprn) > y(n_NOx,lprn)) then
              write(out_line,110)
     &        'losses by reaction     (PAN formation) removed',
     &        chemrate(iPANform,lprn)
              call write_parallel(trim(out_line),crit=jay)
            endif
          endif
          
          call chem1prn
     &    (kpnr,2,nnr,npnr,chemrate,2,1,igas,total,maxl,I,J,jay)
     
          if(igas == n_NOx)then
            if(-dest(n_HO2NO2,lprn) >= y(n_HO2NO2,lprn) .or.
     &      chemrate(iHO2NO2form,lprn) > y(n_NOx,lprn)) then
              write(out_line,110)
     &        'gain by reactions destroying HO2NO2) rmoved  ',
     &        (rr(iHO2NO2_OH,lprn)*y(nOH,L)+
     &        rr(iHO2NO2decomp,lprn)*y(nM,lprn)+
     &        ss(10,lprn,I,J)+ss(11,lprn,I,J))*y(n_HO2NO2,lprn)*dt2
              call write_parallel(trim(out_line),crit=jay)     
            endif
            if(-dest(n_N2O5,lprn) >= y(n_N2O5,lprn).or.
     &      chemrate(iN2O5form,lprn) > y(n_NOx,lprn)) then
              write(out_line,110)
     &        'gains by reaction    (N2O5 decomposition) removed',
     &        2.d0*chemrate(iN2O5decomp,lprn)
              call write_parallel(trim(out_line),crit=jay)
            endif
            if(-dest(n_PAN,lprn) >= y(n_PAN,lprn).or.
     &      chemrate(iPANform,lprn) > y(n_NOx,lprn)) then
              write(out_line,110)
     &        'gain by reaction    (from PAN) removed',
     &        chemrate(iPANdecomp,lprn)
              call write_parallel(trim(out_line),crit=jay)
            endif
          endif

          call chem1prn
     &    (kds,1,ks,nds,photrate,3,-1,igas,total,maxl,I,J,jay)
          call chem1prn
     &    (kps,2,kss,nps,photrate,4,1,igas,total,maxl,I,J,jay)

          if(igas == n_Ox) then
            write(out_line,110)'Ox change due to within NOx rxns  ',
     &      -Oxcorr(lprn)
            call write_parallel(trim(out_line),crit=jay)
          endif
          if(igas == n_NOx)then
            if(-dest(n_N2O5,lprn) >= y(n_N2O5,lprn) .or.
     &      chemrate(iN2O5form,lprn) > y(n_NOx,lprn)) then
              write(out_line,110)'gains by reaction 7'//
     &        ' (N2O5 photolysis) removed',ss(7,lprn,I,J)*
     &        y(n_N2O5,lprn)*dt2
              call write_parallel(trim(out_line),crit=jay)
            endif
            if(-dest(n_N2O5,lprn) >= y(n_N2O5,lprn) .or.
     &      chemrate(iN2O5form,lprn) > y(n_NOx,lprn)) then 
              write(out_line,110)'net change due to N2O5 is ',
     &        2.d0*(y(n_N2O5,lprn)-(rr(iN2O5form,lprn)*y(nNO3,lprn)*
     &        y(nNO2,lprn))/
     &        (rr(iN2O5decomp,lprn)*y(nM,lprn)+ss(7,lprn,I,J)))
              call write_parallel(trim(out_line),crit=jay)
            endif
            if(-dest(n_HO2NO2,lprn) >= y(n_HO2NO2,lprn) .or.
     &      chemrate(iHO2NO2form,lprn) > y(n_NOx,lprn)) then
              write(out_line,110)'gain by rxns 10 & 11 (HO2NO2'
     &        //' photolysis) rmoved',(ss(10,lprn,I,J)+
     &        ss(11,lprn,I,J))*y(n_HO2NO2,lprn)*dt2
              call write_parallel(trim(out_line),crit=jay)
            endif
            if(-dest(n_HO2NO2,lprn) >= y(n_HO2NO2,lprn) .or.
     &      chemrate(iHO2NO2form,lprn) > y(n_NOx,lprn)) then
              write(out_line,110)'net change due to HO2NO2 is ',
     &        y(n_HO2NO2,lprn)-((rr(iHO2NO2form,lprn)*y(nHO2,lprn)*
     &        y(nNO2,lprn))/(rr(iHO2NO2_OH,lprn)*
     &        y(nOH,lprn)+rr(iHO2NO2decomp,lprn)*y(nM,lprn)
     &        +ss(10,lprn,I,J)+ss(11,lprn,I,J)))
              call write_parallel(trim(out_line),crit=jay)
            endif
            if(-dest(n_PAN,lprn) >= y(n_PAN,lprn) .or.
     &      chemrate(iPANform,lprn) > y(n_NOx,lprn)) then
              write(out_line,110)'net change due to PAN is ',
     &        y(n_PAN,lprn)-((rr(iPANform,lprn)*y(nC2O3,lprn)*
     &        y(nNO2,lprn))/
     &        (rr(iPANdecomp,lprn)*y(nM,lprn)+ss(15,lprn,I,J)))              
              call write_parallel(trim(out_line),crit=jay)    
            endif
          endif
                
          if(igas == n_Ox .or. igas == n_NOx) total=
     &    100.d0*(dest(igas,lprn)+prod(igas,lprn))/y(igas,lprn)
     
          if(igas == n_BrOx)then
            if(-dest(n_HOBr,lprn) >= y(n_HOBr,lprn).or.
     &      chemrate(73,lprn) > 0.5d0*y(n_BrOx,lprn))then
              write(out_line,110)
     &        'gain by rxns 24 (HOBr photolysis) removed'
     &        ,ss(24,lprn,i,j)*y(n_HOBr,lprn)*dt2
              call write_parallel(trim(out_line),crit=jay)
              write(out_line,110)'loss by rxn 73 removed'
     &        ,chemrate(73,lprn)
              call write_parallel(trim(out_line),crit=jay)
            endif
            if(-dest(n_BrONO2,lprn) >= y(n_BrONO2,lprn) .or.
     &      chemrate(iClOplusNO2,lprn) > 0.5d0*y(n_BrOx,lprn))then
              write(out_line,110)
     &        'gain by rxns 23 (BrONO2 photolysis) removed'
     &        ,ss(23,lprn,i,j)*y(n_BrONO2,lprn)*dt2
              call write_parallel(trim(out_line),crit=jay)
               write(out_line,110)'loss by rxn iClOplusNO2 removed'
     &        ,chemrate(iClOplusNO2,lprn)
              call write_parallel(trim(out_line),crit=jay)
            endif
          endif
          
          if(igas == n_NOx)then
            if(-dest(n_BrONO2,lprn) >= y(n_BrONO2,lprn) .or.
     &      chemrate(iClOplusNO2,lprn) > 0.5d0*y(n_BrOx,lprn))then
              write(out_line,110)
     &        'gain by rxns 23 (BrONO2 photolysis) removed'
     &        ,ss(23,lprn,i,j)*y(n_BrONO2,lprn)*dt2
              call write_parallel(trim(out_line),crit=jay)
              write(out_line,110)'loss by rxn iClOplusNO2 removed'
     &        ,chemrate(iClOplusNO2,lprn)
              call write_parallel(trim(out_line),crit=jay)     
            endif
          endif
          
          if(igas == n_ClOx)then
            if(-dest(n_HOCl,lprn) >= y(n_HOCl,lprn) .or.
     &      chemrate(63,lprn) > y(n_ClOx,lprn))then
              write(out_line,110)
     &        'gain by rxn 21 (HOCl photolysis) removed'
     &        ,ss(21,lprn,i,j)*y(n_HOCl,lprn)*dt2
              call write_parallel(trim(out_line),crit=jay)
              write(out_line,110)'gain by rxn 55 removed'
     &        ,chemrate(55,lprn)
              call write_parallel(trim(out_line),crit=jay)
              write(out_line,110)'loss by rxn 63 removed'
     &        ,chemrate(63,lprn)
              call write_parallel(trim(out_line),crit=jay)
            endif 
            if(-dest(n_ClONO2,lprn) >= y(n_ClONO2,lprn) .or.
     &      chemrate(iClOplusClO,lprn) > 0.8d0*y(n_ClOx,lprn))then
              write(out_line,110)
     &        'gain by rxn 22 (ClONO2 photolysis) removed'
     &        ,ss(22,lprn,i,j)*y(n_ClONO2,lprn)*dt2
              call write_parallel(trim(out_line),crit=jay)
              write(out_line,110)'gain by rxn 65 removed'
     &        ,chemrate(65,lprn)
              call write_parallel(trim(out_line),crit=jay)
              write(out_line,110)'loss by rxn iClOplusClO removed'
     &        ,chemrate(iClOplusClO,lprn)
              call write_parallel(trim(out_line),crit=jay)
            endif
          endif
        
          if(igas == n_NOx)then
            if(-dest(n_ClONO2,lprn) >= y(n_ClONO2,lprn) .or.
     &      chemrate(iClOplusClO,lprn) > 0.8d0*y(n_ClOx,lprn))then
              write(out_line,110)
     &        'gain by rxn 22 (ClONO2 photolysis) removed'
     &        ,ss(22,lprn,i,j)*y(n_ClONO2,lprn)*dt2
              call write_parallel(trim(out_line),crit=jay)
              write(out_line,110)'gain by rxn 65 removed'
     &        ,chemrate(65,lprn)
              call write_parallel(trim(out_line),crit=jay)
              write(out_line,110)'loss by rxn iClOplusClO removed'
     &        ,chemrate(iClOplusClO,lprn)
              call write_parallel(trim(out_line),crit=jay)
            endif
          endif

          if(igas == n_CH3OOH) then
            write(out_line,'(a48,a6,e10.3)')
     &      'production from XO2N + HO2 ','dy = ',
     &      y(nHO2,lprn)*y(nNO,lprn)*rr(44,lprn)*rr(43,lprn)/
     &      (y(nNO,lprn)*4.2d-12*exp(180.d0/ta(lprn)))
     &      *y(nXO2N,lprn)*dt2            
            call write_parallel(trim(out_line),crit=jay)
          endif

#ifdef TRACERS_HETCHEM
          if(igas == n_HNO3) then
            write(out_line,'(a48,a6,e10.3)')
     &      'destruction from HNO3 +dust ','dy = ',
     &      -y(n_HNO3,lprn)*krate(iprn,jprn,lprn,1,1)*dt2
            call write_parallel(trim(out_line),crit=jay)
          endif
#endif
          if(igas == n_Paraffin) then
            write(out_line,'(a48,a6,e10.3)')'destruction from RXPAR ',
     &      'dy = ',-y(nRXPAR,lprn)*y(n_Paraffin,lprn)*8.d-11*dt2
            call write_parallel(trim(out_line),crit=jay)
          endif
          
          write(out_line,118) ' Total change in ',ay(igas),
     &    ' is ',total,' percent; dy= ',dest(igas,lprn)+prod(igas,lprn)
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,*) ' '
          call write_parallel(trim(out_line),crit=jay)
        enddo ! igas
      endif  ! chem diags
 108  format(a10,2x,a8)
 110  format(a68,e10.3)
 118  format(a17,a8,a4,f10.0,a14,e12.3)

      if(prnchg .and. J == jprn .and. I == iprn) then
      
        write(out_line,*)
     &  'Percentage ozone loss per cycle at I,J:',I,J
        call write_parallel(trim(out_line),crit=jay)
        write(out_line,'(a41,a56,a15)')
     &  '  L   ClOx    NOx     HOx     BrOx    Ox ',
     &  '   NO2+O     NO+O3     ClO+O     Cl+O3    NO2+hv    net',
     &  '    JO2+hv   JNO'     
        call write_parallel(trim(out_line),crit=jay)

        do Lz=LM,LS1,-1
          sumC=rr(45,Lz)*y(nCl,Lz)*y(nO3,Lz)+
     &    rr(46,Lz)*y(nClO,Lz)*y(nO,Lz)+
     &    rr(49,Lz)*y(nClO,Lz)*y(nO3,Lz)+
     &    rr(50,Lz)*y(nOClO,Lz)*y(nO,Lz) ! -ss(17,Lz,i,j)*y(nClO,Lz)
          sumN=rr(5,Lz)*y(nNO,Lz)*y(nO3,Lz)+
     &    rr(26,Lz)*y(nNO2,Lz)*y(nO,Lz)+
     &    rr(7,Lz)*y(nNO2,Lz)*y(nO3,Lz)+
     &    rr(iNOplusO,Lz)*y(nNO,Lz)*y(nO,Lz)
     &    -ss(1,Lz,i,j)*y(nNO2,Lz)
          sumH=rr(2,Lz)*y(nOH,Lz)*y(nO3,Lz)+
     &    rr(4,Lz)*y(nHO2,Lz)*y(nO3,Lz)+
     &    rr(89,Lz)*y(nOH,Lz)*y(nO,Lz)+
     &    rr(90,Lz)*y(nHO2,Lz)*y(nO,Lz)
          sumB=rr(69,Lz)*y(nBrO,Lz)*y(nO,Lz)+
     &    rr(70,Lz)*y(nBr,Lz)*y(nO3,Lz) ! -ss(25,Lz,i,j)*y(nBrO,Lz)
          sumO=2*rr(88,Lz)*y(nO,Lz)*y(nO3,Lz)
          sumA=sumC+sumN+sumH+sumB+sumO
          write(out_line,'(i3,1x,5(f7.2,1x),8(e9.2,1x))')
     &    Lz,100.d0*sumC/sumA,
     &    100.d0*sumN/sumA,100.d0*sumH/sumA,100.d0*sumB/sumA,
     &    100.d0*sumO/sumA,rr(26,Lz)*y(nNO2,Lz)*y(nO,Lz),
     &    rr(5,Lz)*y(nNO,Lz)*y(nO3,Lz),
     &    rr(46,Lz)*y(nClO,Lz)*y(nO,Lz),
     &    rr(45,Lz)*y(nCl,Lz)*y(nO3,Lz),ss(1,Lz,i,j)*y(nNO2,Lz),sumA
     &    ,ss(27,Lz,i,j),SF2(i,j,Lz)
          call write_parallel(trim(out_line),crit=jay)
        enddo
        write(out_line,*) ' '
        call write_parallel(trim(out_line),crit=jay)
        write(out_line,'(a35,3(2x,i2))')
     &  ' Total change by species at I, J, L',i,j,lprn
        call write_parallel(trim(out_line),crit=jay)
      endif ! end of chemistry diagnostics ----------------------------

c Loops to calculate tracer changes:

      rMAbyM(1:maxl) = MA(1:maxl,I,J) / y(nM,1:maxl)

      do igas=1,ntm_chem ! TRACER LOOP -----------------
      
       dxbym2v=axyp(I,J)*vol2mass(igas)
       do L=1,maxl
         conc2mass=rMAbyM(L)*dxbym2v
         c2ml(l) = conc2mass
         changeL(L,igas)=
     &   (dest(igas,L)+prod(igas,L))*conc2mass
         if(igas == n_CO)then
#ifdef HTAP_LIKE_DIAGS
           TAIJLS(I,J,L,ijlt_COp)=TAIJLS(I,J,L,ijlt_COp)+prod(igas,L)
     *          *cpd
           TAIJLS(I,J,L,ijlt_COd)=TAIJLS(I,J,L,ijlt_COd)+dest(igas,L)
     *          *cpd
#endif
         endif
         if(igas == n_Ox)then
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
!NEED      if(trm(i,j,l,n_Ox)==0.)call stop_model('zero ozone',255)
!NEED      changeL(L,n_stratOx)=dest(igas,L)*conc2mass*
!NEED&     trm(i,j,l,n_stratOx)/trm(i,j,l,n_Ox)
!NEED      if(L>maxT)changeL(L,n_stratOx)=changeL(L,n_stratOx)+ 
!NEED&     prod(igas,L)*conc2mass*trm(i,j,l,n_stratOx)/trm(i,j,l,n_Ox)
!NEED      if((trm(i,j,l,n_stratOx)+changeL(l,n_stratOx)) < 1.d0)
!NEED&     changeL(l,n_stratOx) = 1.d0 - trm(i,j,l,n_stratOx)
#endif
#ifdef HTAP_LIKE_DIAGS
           TAIJLS(I,J,L,ijlt_Oxp)=TAIJLS(I,J,L,ijlt_Oxp)+prod(igas,L)
     *          *cpd
           TAIJLS(I,J,L,ijlt_Oxd)=TAIJLS(I,J,L,ijlt_Oxd)+dest(igas,L)
     *          *cpd
         else if(igas==n_CH4)then
           ! destruction only
           TAIJLS(I,J,L,ijlt_CH4d)=
     &          TAIJLS(I,J,L,ijlt_CH4d)+dest(igas,L)*cpd
#endif
         endif
         
c Set N2O5 to equilibrium when necessary (near ground,
c N2O5 is thermally unstable, has a very short lifetime):
         if(igas == n_N2O5.and.(-dest(igas,L) >= y(n_N2O5,L)*0.75d0.or.
     &   chemrate(iN2O5form,L) > y(n_NOx,L)))then
           rnewval=(rr(iN2O5form,L)*y(nNO3,L)*y(nNO2,L))/
     &     (rr(iN2O5decomp,L)*y(nM,L)+ss(7,L,I,J)+1.d-12)
           if(rnewval < 1.d0)rnewval=1.d0
           changeL(L,igas)=(rnewval-y(n_N2O5,L))
           if(changeL(L,igas) > 0.33d0*y(nNO2,L))changeL(L,igas)=
     &     0.33d0*y(nNO2,L)
           changeL(L,igas)=changeL(L,igas)*conc2mass
         endif

c Conserve NOx with respect to N2O5:
         if(igas == n_NOx.and.(-dest(n_N2O5,L) >= y(n_N2O5,L).or.
     &   chemrate(iN2O5form,L) > y(n_NOx,L)))then
           rnewval=(rr(iN2O5form,L)*y(nNO3,L)*y(nNO2,L))/
     &     (rr(iN2O5decomp,L)*y(nM,L)+ss(7,L,I,J)+1.d-12)
           if(rnewval < 1.d0)rnewval=1.d0
           changeX=(rnewval-y(n_N2O5,L))
           if(changeX > 0.33d0*y(nNO2,L))changeX=0.33d0*y(nNO2,L)
           changeL(L,igas)=
     &     changeL(L,igas)-changeX*conc2mass
         endif

c Set HO2NO2 to equil when necessary:
         if(igas == n_HO2NO2.and.(-dest(igas,L) >= y(n_HO2NO2,L).or.
     &   chemrate(iHO2NO2form,L) > y(n_NOx,L)))then
           rnewval=(rr(iHO2NO2form,L)*y(nHO2,L)*y(nNO2,L))/
     &     (rr(iHO2NO2_OH,L)*y(nOH,L)+rr(iHO2NO2decomp,L)*y(nM,L)
     &     +ss(10,L,I,J)+ss(11,L,I,J)+1.d-12)
           if(rnewval < 1.d0)rnewval=1.d0
           changeL(L,igas)=(rnewval-y(n_HO2NO2,L))
           if(changeL(L,igas) > 0.33d0*y(nNO2,L))changeL(L,igas)=
     &     0.33d0*y(nNO2,L)
           changeL(L,igas)=changeL(L,igas)*conc2mass
         endif

c Conserve NOx with respect to HO2NO2:
         if(igas == n_NOx.and.(-dest(n_HO2NO2,L) >= y(n_HO2NO2,L).or.
     &   chemrate(iHO2NO2form,L) > y(n_NOx,L)))then
           rnewval=(rr(iHO2NO2form,L)*y(nHO2,L)*y(nNO2,L))/
     &     (rr(iHO2NO2_OH,L)*y(nOH,L)+rr(iHO2NO2decomp,L)*y(nM,L)
     &     +ss(10,L,I,J)+ss(11,L,I,J)+1.d-12)
           if(rnewval < 1.d0)rnewval=1.d0
           changeX=(rnewval-y(n_HO2NO2,L))
           if(changeX > 0.33d0*y(nNO2,L))changeX=0.33d0*y(nNO2,L)
           changeL(L,igas)=changeL(L,igas)-
     &     changeX*conc2mass
         endif

c Set PAN to equilibrium when necessary (near ground,
c PAN is thermally unstable, has a very short lifetime):
         if(igas == n_PAN.and.(-dest(igas,L) >= y(n_PAN,L).or.
     &   chemrate(iPANform,L) > y(n_NOx,L)))then
           rnewval=(rr(iPANform,L)*y(nC2O3,L)*y(nNO2,L))/
     &     (rr(iPANdecomp,L)*y(nM,L)+ss(15,L,I,J)+1.d-12)
           if(rnewval < 1.d0)rnewval=1.d0
           changeL(L,igas)=(rnewval-y(n_PAN,L))
           if(changeL(L,igas) > 0.33d0*y(nNO2,L))changeL(L,igas)=
     &     0.33d0*y(nNO2,L)
           changeL(L,igas)=changeL(L,igas)*conc2mass
         endif

c Conserve NOx with respect to PAN:
         if(igas == n_NOx.and.(-dest(n_PAN,L) >= y(n_PAN,L).or.
     &   chemrate(iPANform,L) > y(n_NOx,L)))then
           rnewval=(rr(iPANform,L)*y(nC2O3,L)*y(nNO2,L))/
     &     (rr(iPANdecomp,L)*y(nM,L)+ss(15,L,I,J)+1.d-12)
           if(rnewval < 1.d0)rnewval=1.d0
           changeX=(rnewval-y(n_PAN,L))
           if(changeX > 0.33d0*y(nNO2,L))changeX=0.33d0*y(nNO2,L)
           changeL(L,igas)=changeL(L,igas)-changeX*
     &     conc2mass
         endif

c Cacluate Cl2 amount to P/L:
         if((ss(18,L,I,J)+rr(51,L)*y(nOH,L)) > 0.)then
           y(nCl2,L)=rr(57,L)*y(n_HOCl,L)*y(nCl,L) / 
     &     (ss(18,L,I,J)+rr(51,L)*y(nOH,L) + 1.d-12)
         else
           y(nCl2,L)=0.d0
         endif
         yCl2(I,J,L)=y(nCl2,L)

c Set HOBr to equilibrium when necessary:
         if(igas == n_HOBr.and.(-dest(igas,L) >= y(n_HOBr,L).or.
     &     chemrate(73,L) > 0.5d0*y(n_BrOx,L)))then
           rnewval=(rr(73,L)*y(nBrO,L)*y(nHO2,L))/
     &     (ss(24,L,i,j)+1.d-12)
           if(rnewval < 1.d0)rnewval=1.d0
           changeL(L,igas)=(rnewval-y(n_HOBr,L))
           if(changeL(L,igas) > 0.5d0*y(nBrO,L))changeL(L,igas)=
     &     0.5d0*y(nBrO,L)
           changeL(L,igas)=changeL(L,igas)*conc2mass
         endif

c Conserve BrOx with respect to HOBr:
         if(igas == n_BrOx.and.(-dest(n_HOBr,L) >= y(n_HOBr,L).or.
     &   chemrate(73,L) > 0.5d0*y(n_BrOx,L)))then
           rnewval=(rr(73,L)*y(nBrO,L)*y(nHO2,L))/
     &     (ss(24,L,i,j)+1.d-12)
           if(rnewval < 1.d0)rnewval=1.d0
           changeX=(rnewval-y(n_HOBr,L))
           if(changeX > 0.5d0*y(nBrO,L))changeX=0.5d0*y(nBrO,L)
           changeL(L,igas)=changeL(L,igas)-
     &     changeX*conc2mass
         endif

c Set BrONO2 to equilibrium when necessary:
         if(igas == n_BrONO2.and.(-dest(igas,L) >= y(n_BrONO2,L).or.
     &   chemrate(iBrOplusNO2,L) > 0.5d0*y(n_BrOx,L)))then
           rnewval=(rr(iBrOplusNO2,L)*y(nBrO,L)*y(nNO2,L))/
     &     (ss(23,L,i,j)+1.d-12)
           if(rnewval < 1.d0)rnewval=1.d0
           changeL(L,igas)=(rnewval-y(n_BrONO2,L))
           if(changeL(L,igas) > 0.5d0*y(nBrO,L))changeL(L,igas)=
     &     0.5d0*y(nBrO,L)
           changeL(L,igas)=changeL(L,igas)*conc2mass
         endif

c Conserve BrOx with respect to BrONO2:
         if(igas == n_BrOx.and.(-dest(n_BrONO2,L) >= y(n_BrONO2,L).or.
     &   chemrate(iBrOplusNO2,L) > 0.5d0*y(n_BrOx,L)))then
           rnewval=(rr(iBrOplusNO2,L)*y(nBrO,L)*y(nNO2,L))/
     &     (ss(23,L,i,j)+1.d-12)
           if(rnewval < 1.d0)rnewval=1.d0
           changeX=(rnewval-y(n_BrONO2,L))
           if(changeX > 0.5d0*y(nBrO,L))changeX=0.5d0*y(nBrO,L)
           changeL(L,igas)=changeL(L,igas)-changeX*
     &     conc2mass
         endif

c Conserve NOx with respect to BrONO2:
         if(igas == n_NOx.and.(-dest(n_BrONO2,L) >= y(n_BrONO2,L).or.
     &   chemrate(iBrOplusNO2,L) > 0.5d0*y(n_BrOx,L)))then
           rnewval=(rr(iBrOplusNO2,L)*y(nBrO,L)*y(nNO2,L))/
     &     (ss(23,L,i,j)+1.d-12)
           if(rnewval < 1.d0)rnewval=1.d0
           changeX=(rnewval-y(n_BrONO2,L))
           if(changeX > 0.5d0*y(nBrO,L))changeX=0.5d0*y(nBrO,L)
           changeL(L,igas)=changeL(L,igas)-changeX*
     &     conc2mass
         endif

c Set ClONO2 to equilibrium when necessary:
         if(igas == n_ClONO2.and.(-dest(igas,L) >= y(n_ClONO2,L).or.
     &   chemrate(iClOplusNO2,L) > 0.8d0*y(n_ClOx,L)))then
           rnewval=(rr(iClOplusNO2,L)*y(nClO,L)*y(nNO2,L))/
     &     (ss(22,L,i,j)+rr(65,L)*y(nO,L)+1.d-12)
           if(rnewval < 1.d0)rnewval=1.d0
           changeL(L,igas)=(rnewval-y(n_ClONO2,L))
           if(changeL(L,igas) > 0.3d0*y(nClO,L))changeL(L,igas)=
     &     0.3d0*y(nClO,L)
           if(-changeL(L,igas) > 0.8d0*y(n_ClONO2,L))
     &     changeL(L,igas)=-0.8d0*y(n_ClONO2,L)
           changeL(L,igas)=changeL(L,igas)*conc2mass
         endif

c Conserve ClOx with respect to ClONO2:
         if(igas == n_ClOx.and.(-dest(n_ClONO2,L) >= y(n_ClONO2,L).or.
     &   chemrate(iClOplusNO2,L) > 0.8d0*y(n_ClOx,L)))then
           rnewval=(rr(iClOplusNO2,L)*y(nClO,L)*y(nNO2,L))/
     &     (ss(22,L,i,j)+rr(65,L)*y(nO,L)+1.d-12)
           if(rnewval < 1.d0)rnewval=1.d0
           changeX=(rnewval-y(n_ClONO2,L))
           if(changeX > 0.3d0*y(nClO,L))changeX=0.3d0*y(nClO,L)
           if(-changeX > 0.8d0*y(n_ClONO2,L))changeX=
     &     -0.8d0*y(n_ClONO2,L)
           changeL(L,igas)=changeL(L,igas)-changeX*
     &     conc2mass
         endif

c Conserve NOx with respect to ClONO2:
         if(igas == n_NOx.and.(-dest(n_ClONO2,L) >= y(n_ClONO2,L).or.
     &   chemrate(iClOplusNO2,L) > 0.8d0*y(n_ClOx,L)))then
           rnewval=(rr(iClOplusNO2,L)*y(nClO,L)*y(nNO2,L))/
     &     (ss(22,L,i,j)+rr(65,L)*y(nO,L)+1.d-12)
           if(rnewval < 1.d0)rnewval=1.d0
           changeX=(rnewval-y(n_ClONO2,L))
           if(changeX > 0.3d0*y(nClO,L))changeX=0.3d0*y(nClO,L)
           if(-changeX > 0.8d0*y(n_ClONO2,L))changeX=
     &     -0.8d0*y(n_ClONO2,L)
           changeL(L,igas)=changeL(L,igas)-changeX*
     &     conc2mass
         endif

c Set HOCl to equilibrium when necessary:
         if(igas == n_HOCl.and.(-dest(igas,L) >= y(n_HOCl,L).or.
     &   chemrate(63,L) > y(n_ClOx,L)))then
           rnewval=(rr(63,L)*y(nClO,L)*y(nHO2,L) + 
     &     rr(51,L)*y(nCl2,L)*y(nOH,L)) /
     &     (ss(21,L,i,j)+rr(55,L)*y(nO,L)+rr(51,L)*y(nCl2,L)+1.d-12)
           if(rnewval < 1.d0)rnewval=1.d0
           changeL(L,igas)=(rnewval-y(n_HOCl,L))
           if(changeL(L,igas) > 0.3d0*y(nClO,L))changeL(L,igas)=
     &     0.3d0*y(nClO,L)
           changeL(L,igas)=changeL(L,igas)*conc2mass
         endif

c Conserve ClOx with respect to HOCl:
         if(igas == n_ClOx.and.(-dest(n_HOCl,L) >= y(n_HOCl,L).or.
     &   chemrate(63,L) > y(n_ClOx,L)))then
           rnewval=(rr(63,L)*y(nClO,L)*y(nHO2,L) + 
     &     rr(51,L)*y(nCl2,L)*y(nOH,L)) /
     &     (ss(21,L,i,j)+rr(55,L)*y(nO,L)+rr(51,L)*y(nCl2,L)+1.d-12)
           if(rnewval < 1.d0)rnewval=1.d0
           changeX=(rnewval-y(n_HOCl,L))
           if(changeX > 0.3d0*y(nClO,L))changeX=0.3d0*y(nClO,L)
           changeL(L,igas)=changeL(L,igas)-changeX*
     &     conc2mass
         endif
       end do ! L

       if(igas == n_CO)then
         call inc_tajls_column(i,j,1,maxl,maxl,jls_COp,
     &        prod(igas,1:maxl)*c2ml(1:maxl))
         call inc_tajls_column(i,j,1,maxl,maxl,jls_COd,
     &        dest(igas,1:maxl)*c2ml(1:maxl))
       endif
       if(igas == n_Ox)then
         call inc_tajls_column(i,j,1,maxl,maxl,jls_Oxp ,
     &        prod(igas,1:maxl)*c2ml(1:maxl))
         call inc_tajls_column(i,j,1,maxT,maxT,jls_OxpT,
     &        prod(igas,1:maxT)*c2ml(1:maxT))
         call inc_tajls_column(i,j,1,maxl,maxl,jls_Oxd ,
     &        dest(igas,1:maxl)*c2ml(1:maxl))
         call inc_tajls_column(i,j,1,maxT,maxT,jls_OxdT,
     &        dest(igas,1:maxT)*c2ml(1:maxT))
       endif
      end do  ! igas ! end of TRACER LOOP -----------------

c Separate N2O change for N cons, leave out N2O->N2+O fromm cons:
      sv_changeN2O(1:maxl)=
     &-chemrate(87,1:maxl)*axyp(i,j)*rMAbyM(1:maxl)*vol2mass(n_N2O)

c Ensure nitrogen conservation,
c (since equilibration of short lived gases may alter this):

      if(prnchg .and. J == jprn .and. I == iprn)then
        write(out_line,*)
     &  'changes (mass) before nitrogen conservation routine'
        call write_parallel(trim(out_line),crit=jay)
        write(out_line,*) 'NOx, N2O5, HO2NO2, HNO3, PAN, AlkylNit, N2O'
        call write_parallel(trim(out_line),crit=jay)
        write(out_line,*) 'ClONO2, BrONO2'
        call write_parallel(trim(out_line),crit=jay)
        write(out_line,*) changeL(lprn,n_NOx),changeL(lprn,n_N2O5),
     &  changeL(lprn,n_HO2NO2),changeL(lprn,n_HNO3),
     &  changeL(lprn,n_PAN),changeL(lprn,n_AlkylNit)
     &  ,changeL(lprn,n_N2O)
     &  ,changeL(lprn,n_ClONO2),changeL(lprn,n_BrONO2)
        call write_parallel(trim(out_line),crit=jay)
        write(out_line,*)
     &  'N2O change w/o rxns forming N2',sv_changeN2O(lprn)
        call write_parallel(trim(out_line),crit=jay)
#ifdef TRACERS_HETCHEM
        write(out_line,*) 'HNO3 loss on dust replaced for cons ',
     &  (krate(i,j,lprn,1,1)*y(n_HNO3,lprn)*dt2)*rMAbyM(lprn)*axyp(I,J)
        call write_parallel(trim(out_line),crit=jay)
#endif
      endif

      do L=1,maxl ! start big L-LOOP ---------------

c First check for nitrogen loss > 100% :
        if(-changeL(L,n_NOx) > trm(I,J,L,n_NOx))
     &  changeL(L,n_NOx)=1.d0-trm(I,J,L,n_NOx)
        if(-changeL(L,n_N2O5) > trm(I,J,L,n_N2O5))
     &  changeL(L,n_N2O5)=1.d0-trm(I,J,L,n_N2O5)
        if(-changeL(L,n_HO2NO2) > trm(I,J,L,n_HO2NO2))
     &  changeL(L,n_HO2NO2)=1.d0-trm(I,J,L,n_HO2NO2)
        if(-changeL(L,n_HNO3) > trm(I,J,L,n_HNO3))
     &  changeL(L,n_HNO3)=1.d0-trm(I,J,L,n_HNO3)
        if(-changeL(L,n_PAN) > trm(I,J,L,n_PAN))
     &  changeL(L,n_PAN)=1.d0-trm(I,J,L,n_PAN)
        if(-changeL(L,n_AlkylNit) > trm(I,J,L,n_AlkylNit))
     &  changeL(L,n_AlkylNit)=1.d0-trm(I,J,L,n_AlkylNit)
        if(-changeL(L,n_ClONO2) > trm(I,J,L,n_ClONO2))
     &  changeL(L,n_ClONO2)=1.d0-trm(I,J,L,n_ClONO2)
        if(-changeL(L,n_BrONO2) > trm(I,J,L,n_BrONO2))
     &  changeL(L,n_BrONO2)=1.d0-trm(I,J,L,n_BrONO2)
#ifdef TRACERS_HETCHEM
        changeL(L,n_HNO3)=changeL(L,n_HNO3)+(krate(i,j,l,1,1)
     &  *y(n_HNO3,l)*dt2)*rMAbyM(L)*axyp(i,j)*vol2mass(n_HNO3)
        if(prnchg .and. i == iprn .and. j == jprn) then
          write(out_line,*)
     &    changeL(L,n_HNO3),krate(i,j,l,1,1),y(n_HNO3,l)
          call write_parallel(trim(out_line),crit=jay)
        endif   
#endif

c Next insure balance between dNOx and sum of dOthers:
        sumN=(2.d0*changeL(L,n_N2O5))*mass2vol(n_N2O5)+
     &  (changeL(L,n_HNO3))*mass2vol(n_HNO3)+
     &  (changeL(L,n_HO2NO2))*mass2vol(n_HO2NO2)+
     &  (changeL(L,n_PAN))*mass2vol(n_PAN)+
     &  (changeL(L,n_AlkylNit))*mass2vol(n_AlkylNit)

        if(prnchg.and.J==jprn.and.I==iprn.and.L==lprn) then
          write(out_line,*) 'ratio for conservation =',ratio
          call write_parallel(trim(out_line),crit=jay)
        endif
        
        sumN=sumN+
     &  changeL(L,n_ClONO2)*mass2vol(n_ClONO2)+
     &  changeL(L,n_BrONO2)*mass2vol(n_BrONO2)
        dNOx=changeL(L,n_NOx)*mass2vol(n_NOx)+
     &  2.d0*sv_changeN2O(L)*mass2vol(n_N2O)
        if(prnchg.and.J==jprn.and.I==iprn.and.L==lprn) then
          write(out_line,*)
     &    'other N changes, dNOx (less prod fm N2O) = (molec) ',
     &    sumN,dNOx
          call write_parallel(trim(out_line),crit=jay)
        endif
        ratio=-sumN/dNOx
        if(ratio <= 0.999d0 .or. ratio >= 1.001d0) then
         if(dNOx > 0.d0)then ! NOx produced (net positive change)
          if (ratio > 1.d0)then
           sumD=0.d0
c          reduce N destruction to match NOx prodcution:
           if(changeL(L,n_N2O5) < 0.d0)   sumD=sumD+
     &     2.d0*changeL(L,n_N2O5)*mass2vol(n_N2O5)
           if(changeL(L,n_HO2NO2) < 0.d0) sumD=sumD+
     &     changeL(L,n_HO2NO2)*mass2vol(n_HO2NO2)
           if(changeL(L,n_HNO3) < 0.d0)   sumD=sumD+
     &     changeL(L,n_HNO3)*mass2vol(n_HNO3)
           if(changeL(L,n_PAN) < 0.d0)    sumD=sumD+
     &     changeL(L,n_PAN)*mass2vol(n_PAN)
           if(changeL(L,n_AlkylNit) < 0.d0)sumD=sumD+
     &     changeL(L,n_AlkylNit)*mass2vol(n_AlkylNit)
           if(changeL(L,n_ClONO2) < 0.d0)sumD=sumD+
     &     changeL(L,n_ClONO2)*mass2vol(n_ClONO2)
           if(changeL(L,n_BrONO2) < 0.d0)sumD=sumD+
     &     changeL(L,n_BrONO2)*mass2vol(n_BrONO2)
           newD=(sumN/ratio)+sumD-sumN
           ratioD=newD/sumD
           if(changeL(L,n_N2O5) < 0.d0)    changeL(L,n_N2O5)=
     &     changeL(L,n_N2O5)    *ratioD
           if(changeL(L,n_HO2NO2) < 0.d0)  changeL(L,n_HO2NO2)=
     &     changeL(L,n_HO2NO2)  *ratioD
           if(changeL(L,n_HNO3) < 0.d0)    changeL(L,n_HNO3)=
     &     changeL(L,n_HNO3)    *ratioD
           if(changeL(L,n_PAN) < 0.d0)     changeL(L,n_PAN)=
     &     changeL(L,n_PAN)     *ratioD
           if(changeL(L,n_AlkylNit) < 0.d0)changeL(L,n_AlkylNit)=
     &     changeL(L,n_AlkylNit)*ratioD
           vClONO2=changeL(L,n_ClONO2)*(1.d0-ratioD)
           if(changeL(L,n_ClONO2) < 0.d0)changeL(L,n_ClONO2)=
     &     changeL(L,n_ClONO2)*ratioD
           changeL(L,n_ClOx)=changeL(L,n_ClOx)+vClONO2*
     &     (mass2vol(n_ClONO2)*vol2mass(n_ClOx)) !ensure Cl cons
           vBrONO2=changeL(L,n_BrONO2)*(1.d0-ratioD)
           if(changeL(L,n_BrONO2) < 0.d0)changeL(L,n_BrONO2)=
     &     changeL(L,n_BrONO2)*ratioD
           changeL(L,n_BrOx)=changeL(L,n_BrOx)+vBrONO2*
     &     (mass2vol(n_BrONO2)*vol2mass(n_BrOx)) !ensure Br cons
          endif

          if (ratio <= 1.d0 .and. ratio > 0.d0)then
c          reduce NOx production to match N loss:
           changeL(L,n_NOx)=changeL(L,n_NOx)*ratio
           changeL(L,n_NOx)=changeL(L,n_NOx)-
     &     2.d0*sv_changeN2O(L)*mass2vol(n_N2O)*vol2mass(n_NOx)
          endif
         else       ! NOx destroyed (net change is negative):
          if (ratio > 1.d0)then
           sumP=0.d0
c          reduce N production to match NOx loss:
           if(changeL(L,n_N2O5) > 0.d0)    sumP=sumP+
     &     2.d0*changeL(L,n_N2O5)*mass2vol(n_N2O5)
           if(changeL(L,n_HO2NO2) > 0.d0)  sumP=sumP+
     &     changeL(L,n_HO2NO2)*mass2vol(n_HO2NO2)
           if(changeL(L,n_HNO3) > 0.d0)    sumP=sumP+
     &     changeL(L,n_HNO3)*mass2vol(n_HNO3)
           if(changeL(L,n_PAN) > 0.d0)     sumP=sumP+
     &     changeL(L,n_PAN)*mass2vol(n_PAN)
           if(changeL(L,n_AlkylNit) > 0.d0)sumP=sumP+
     &     changeL(L,n_AlkylNit)*mass2vol(n_AlkylNit)
           if(changeL(L,n_ClONO2) > 0.d0)sumP=sumP+
     &      changeL(L,n_ClONO2)*mass2vol(n_ClONO2)
           if(changeL(L,n_BrONO2) > 0.d0)sumP=sumP+
     &      changeL(L,n_BrONO2)*mass2vol(n_BrONO2)
           newP=(sumN/ratio)+sumP-sumN
           if (sumP == 0.) then
             write(out_line,*)'SUMP = 0***', sumP, L,i,j,
     &       changeL(L,n_HNO3)*mass2vol(n_HNO3)
             call write_parallel(trim(out_line),crit=.true.)
           endif
           ratioP=newP/sumP
           if(changeL(L,n_N2O5) > 0.d0)     changeL(L,n_N2O5)=
     &     changeL(L,n_N2O5)*ratioP
           if(changeL(L,n_HO2NO2) > 0.d0)  changeL(L,n_HO2NO2)=
     &     changeL(L,n_HO2NO2)*ratioP
           if(changeL(L,n_HNO3) > 0.d0)    changeL(L,n_HNO3)=
     &     changeL(L,n_HNO3)*ratioP
           if(changeL(L,n_PAN) > 0.d0)     changeL(L,n_PAN)=
     &     changeL(L,n_PAN)*ratioP
           if(changeL(L,n_AlkylNit) > 0.d0)changeL(L,n_AlkylNit)=
     &     changeL(L,n_AlkylNit)*ratioP
           vClONO2=changeL(L,n_ClONO2)*(1.d0-ratioP)
           if(changeL(L,n_ClONO2) > 0.d0)changeL(L,n_ClONO2)=
     &     changeL(L,n_ClONO2)*ratioP
           changeL(L,n_ClOx)=changeL(L,n_ClOx)+vClONO2*
     &     (mass2vol(n_ClONO2)*vol2mass(n_ClOx)) !ensure Cl cons
           vBrONO2=changeL(L,n_BrONO2)*(1.d0-ratioP)
           if(changeL(L,n_BrONO2) > 0.d0)changeL(L,n_BrONO2)=
     &     changeL(L,n_BrONO2)*ratioP
           changeL(L,n_BrOx)=changeL(L,n_BrOx)+vBrONO2*
     &     (mass2vol(n_BrONO2)*vol2mass(n_BrOx)) !ensure Br cons
          endif
          if (ratio <= 1.d0 .and. ratio > 0.d0)then
c          reduce NOx destruction to match N production:
           changeL(L,n_NOx)=changeL(L,n_NOx)*ratio
           changeL(L,n_NOx)=changeL(L,n_NOx)-
     &     2.d0*sv_changeN2O(L)*mass2vol(n_N2O)*vol2mass(n_NOx)
          endif
         endif
#ifdef TRACERS_HETCHEM
         changeL(L,n_HNO3)=changeL(L,n_HNO3)-(krate(i,j,l,1,1)
     &   *y(n_HNO3,l)*dt2)*rMAbyM(L)*axyp(I,J)*vol2mass(n_HNO3)
         changeL(L,n_N_d1)=changeL(L,n_N_d1)+(krate(i,j,l,2,1)
     &   *y(n_HNO3,l)*dt2)*rMAbyM(L)*axyp(I,J)*vol2mass(n_HNO3)
         changeL(L,n_N_d2)=changeL(L,n_N_d2)+(krate(i,j,l,3,1)
     &   *y(n_HNO3,l)*dt2)*rMAbyM(L)*axyp(I,J)*vol2mass(n_HNO3)
         changeL(L,n_N_d3)=changeL(L,n_N_d3)+(krate(i,j,l,4,1)
     &   *y(n_HNO3,l)*dt2)*rMAbyM(L)*axyp(I,J)*vol2mass(n_HNO3)
#endif

        end if ! skipped section above if ratio very close to one

        if(prnchg.and.J==jprn.and.I==iprn.and.L==lprn) then
          write(out_line,*) 'ratio for conservation =',ratio
          call write_parallel(trim(out_line),crit=jay)
        endif

c       Calculate NOx and Ox changes due to atomic nitrogen
c       produced by SRB photlysis (SF2 is NO + hv rate) :
        byta=1.d0/ta(L)
c       rxnN1=3.8d-11*exp(85d0*byta)*y(nOH,L)
        ! that's N+OH->NO+H, not in JPL (rates from IUPAC 1989)
        rxnN2=1.5d-11*exp(-3600.d0*byta)*y(nO2,L) ! N+O2->NO+O
          rxnN3=5.8d-12*exp(220.d0*byta)*y(nNO2,L)  ! N+O2->N2O+O
          rxnN4=2.1d-11*exp(100.d0*byta)*y(nNO,L)   ! N+O2->N2+O
          NprodOx=2.0d0*SF2(I,J,L)*y(nNO,L)*dt2               
          NlossNOx=3.0d1*NprodOx*(rxnN3+rxnN4)/(rxnN2+rxnN3+rxnN4)
        changeL(L,n_NOx)=changeL(L,n_NOx)-NlossNOx
     &  *(axyp(I,J)*rMAbyM(L))*vol2mass(n_NOx)
        conc2mass=axyp(I,J)*rMAbyM(L)*vol2mass(n_Ox)
        changeL(L,n_Ox)=changeL(L,n_Ox)+NprodOx*conc2mass
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
!NEED   if(L>maxT .or. NprodOx<0.)then
!NEED     changeL(L,n_stratOx)=changeL(L,n_stratOx)+
!NEED&    NprodOx*conc2mass*trm(i,j,l,n_stratOx)/trm(i,j,l,n_Ox)
!NEED     if((trm(i,j,l,n_stratOx)+changeL(l,n_stratOx)) < 1.d0)
!NEED&    changeL(l,n_stratOx) = 1.d0 - trm(i,j,l,n_stratOx)
!NEED   endif
#endif
        if(NprodOx <  0.) then ! necessary?
          NprodOx_pos(l) = 0.
          NprodOx_neg(l) = NprodOx*conc2mass
#ifdef HTAP_LIKE_DIAGS
          TAIJLS(I,J,L,ijlt_Oxd)=TAIJLS(I,J,L,ijlt_Oxd)+NprodOx*cpd
#endif
        else 
          NprodOx_neg(l) = 0.
          NprodOx_pos(l) = NprodOx*conc2mass
#ifdef HTAP_LIKE_DIAGS
          TAIJLS(I,J,L,ijlt_Oxp)=TAIJLS(I,J,L,ijlt_Oxp)+NprodOx*cpd
#endif
        endif 
        if(prnchg.and.J==jprn.and.I==iprn.and.l==lprn) then
          write(out_line,*) 'NOx loss & Ox gain due to rxns  w/ N '
     &    ,NlossNOx,NprodOx
          call write_parallel(trim(out_line),crit=jay)
        endif
      end do ! end big L loop -----------------

      call inc_tajls_column(i,j,1,maxl,lm,jls_Oxd ,NprodOx_neg)
      call inc_tajls_column(i,j,1,maxT,lm,jls_OxdT,NprodOx_neg)
      call inc_tajls_column(i,j,1,maxl,lm,jls_Oxp ,NprodOx_pos)
      call inc_tajls_column(i,j,1,maxT,lm,jls_OxpT,NprodOx_pos)

! Remove some of the HNO3 formed heterogeneously, as it doesn't come
! back to the gas phase:
      do L=1,maxL
        if(pscX(L)) changeL(L,n_HNO3)=changeL(L,n_HNO3)-
     &  2.0d-3*y(n_HNO3,L)*(axyp(i,j)*rMAbyM(L))*vol2mass(n_HNO3)
      enddo

c Print chemical changes in a particular grid box if desired:
      if(prnchg .and. J==jprn .and. I==iprn)then
       do igas=1,ntm_chem
         changeA=changeL(Lprn,igas)*y(nM,lprn)*mass2vol(igas)*
     &   byaxyp(i,J)*byMA(lprn,I,J)
         if(y(igas,lprn) == 0.d0)then
           write(out_line,156) ay(igas),': ',changeA,' molecules;  y=0'
           call write_parallel(trim(out_line),crit=jay)
         else
           write(out_line,155)ay(igas),': ',changeA
     &     ,' molecules produced; ',
     &     (100.d0*changeA)/y(igas,lprn),' percent of'
     &     ,y(igas,lprn),'(',1.d9*y(igas,lprn)/y(nM,lprn),' ppbv)'
           call write_parallel(trim(out_line),crit=jay)
         endif

         if(igas == ntm_chem)then
         if(LPRN >= maxT+1)then
            write(out_line,155) ay(nH2O),': ',
     &      changeH2O(lprn),' molecules produced; ',
     &      (100*changeH2O(lprn))/y(nH2O,lprn),' percent of',
     &      y(nH2O,lprn),'(',1.d6*y(nH2O,lprn)/y(nM,lprn),' ppmv)'
            call write_parallel(trim(out_line),crit=jay)
         else
          write(out_line,'(a10,58x,e13.3,6x,f10.3,a5)')
     &    ' H2O     :',y(nH2O,LPRN),(y(nH2O,LPRN)/
     &    y(nM,LPRN))*1.d6,' ppmv'
          call write_parallel(trim(out_line),crit=jay)
         endif
          write(out_line,'(a10,58x,e13.3,6x,f10.3,a5)')
     &    ' CH3O2   :',yCH3O2(I,J,LPRN),(yCH3O2(I,J,LPRN)/
     &    y(nM,LPRN))*1.d9,' ppbv'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,'(a10,58x,e13.3,6x,f10.3,a5)')
     &    ' C2O3    :',y(nC2O3,LPRN),(y(nC2O3,LPRN)/y(nM,LPRN))*1.d9,
     &    ' ppbv'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,'(a10,58x,e13.3,6x,f10.3,a5)')
     &    ' XO2     :',y(nXO2,LPRN),(y(nXO2,LPRN)/y(nM,LPRN))*1.d9,
     &    ' ppbv'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,'(a10,58x,e13.3,6x,f10.3,a5)')
     &    ' XO2N    :',y(nXO2N,LPRN),(y(nXO2N,LPRN)/
     &    y(nM,LPRN))*1.d9,' ppbv'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,'(a10,58x,e13.3,6x,f10.3,a5)')
     &    ' RXPAR   :',y(nRXPAR,LPRN),(y(nRXPAR,LPRN)/
     &    y(nM,LPRN))*1.d9,' ppbv'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,'(a10,58x,e13.3,6x,f10.3,a5)')
     &    ' Aldehyde:',y(nAldehyde,LPRN),(y(nAldehyde,LPRN)/
     &    y(nM,LPRN))*1.d9,' ppbv'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,'(a10,58x,e13.3,6x,f10.3,a5)')
     &    ' ROR     :',y(nROR,LPRN),(y(nROR,LPRN)/
     &    y(nM,LPRN))*1.d9,' ppbv'
          call write_parallel(trim(out_line),crit=jay)
         endif
       enddo
      endif  !end this section of chem diags 

C Tracer masses & slopes are now updated in apply_tracer_3Dsource,
C so here, just saved in changeL:
      do igas=1,ntm_chem
       do L=1,maxl
c Limit the change due to chemistry:
        if(changeL(L,igas) > 1.d20) then
          WRITE(out_line,*)
     &    'change set to 0 in chemstep: I,J,L,igas,change'
     &    ,I,J,L,igas,changeL(L,igas)
          call write_parallel(trim(out_line),unit=99,crit=.true.)
          changeL(L,igas) = 0.d0
        endif
        if(-changeL(L,igas) > trm(I,J,L,igas)) THEN
          if(prnchg)then
            WRITE(out_line,*)
     &      'change > mass, so use 95%: I,J,L,igas,change'
     &      ,I,J,L,igas,changeL(L,igas)
            call write_parallel(trim(out_line),unit=99,crit=.true.)
          endif
          changeL(L,igas) = -0.95d0*trm(I,J,L,igas)
        endif
       end do    ! L
      end do     ! igas

C**** special diags not associated with a particular tracer
      
      DO L=1,maxl
        conOH(l) = 0.
        if (y(nOH,L) > 0.d0 .and. y(nOH,L) < 1.d20)then
          conOH(l) = y(nOH,L)
          TAIJLS(I,J,L,ijlt_OH)=TAIJLS(I,J,L,ijlt_OH)+y(nOH,L)
#ifdef HTAP_LIKE_DIAGS
     &                                             /y(nM,L)
#endif
        end if
        if (y(nHO2,L) > 0.d0 .and. y(nHO2,L) < 1.d20)
     &       TAIJLS(I,J,L,ijlt_HO2)=TAIJLS(I,J,L,ijlt_HO2)+y(nHO2,L)
        conClO(l) = 0.
        if (y(nClO,L) > 0.d0 .and. y(nClO,L) < 1.d20)
     &       conClO(l) = y(nClO,L)/y(nM,L)
        conH2O(l) = 0.
        if (y(nH2O,L) > 0.d0 .and. y(nH2O,L) < 1.d20)
     &       conH2O(l) = y(nH2O,L)/y(nM,L)
      END DO
      call inc_tajls2_column(i,j,1,maxl,lm,jls_OHcon,conOH)
      call inc_tajls2_column(i,j,1,maxl,lm,jls_ClOcon,conClO)
      call inc_tajls2_column(i,j,1,maxl,lm,jls_H2Ocon,conH2O)
      CALL INC_TAJLS2(I,J,1,jls_day,1.d0)

 155  format(1x,a8,a2,e13.3,a21,f10.0,a11,2x,e13.3,3x,a1,f12.5,a6)
 156  format(1x,a8,a2,e13.3,a16)

      return
      end SUBROUTINE chemstep



      SUBROUTINE rates(maxl,I,J)
!@sum rates calculate reaction rates with present concentrations
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
c
C**** GLOBAL parameters and variables:

      USE TRCHEM_Shindell_COM, only: nr,chemrate,photrate,rr,y,nn,dt2,
     &                          ss,ks,ny,dest,prod,JPPJ,nhet

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var kalt local dummy L-loop variable
!@var maxl passed highest chemistry level
!@var ireac,igas dummy loop variables
!@var I,J passed horizontal spatial indicies
      INTEGER :: kalt, ireac, igas, maxl
      INTEGER, INTENT(IN) :: I,J

C Set up rates:
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
        
c Initialize change arrays:
        do igas=1,ny
          dest(igas,kalt)=0.d0
          prod(igas,kalt)=0.d0
        enddo
      enddo
      return
      end SUBROUTINE rates



      SUBROUTINE chem1(kdnr,maxl,numel,nn,ndnr,chemrate,dest,multip)
!@sum chem1 calculate chemical destruction/production
!@auth Drew Shindell (modelEifications by Greg Faluvegi)

C**** GLOBAL parameters and variables:

      USE TRCHEM_Shindell_COM, only: p_2, p_3, p_4, ny, numfam,nfam

      IMPLICIT NONE

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

      ireac=0
      
c Reactive families:

      do igas=1,numfam
        dk=kdnr(igas+1)-kdnr(igas)
        if(dk >= 1) then
          do i=1,dk
            ireac=ireac+1
            do ial=1,maxl
              if(nn(1,ndnr(ireac)) >= nfam(igas) .and. 
     &        nn(1,ndnr(ireac)) < nfam(igas+1))then
                dest(igas,ial)=dest(igas,ial)+multip
     &          *chemrate(ndnr(ireac),ial)
c               Save change array for individual family elements:
                dest(nn(1,ndnr(ireac)),ial)=dest(nn(1,ndnr(ireac)),ial)
     &          + multip*chemrate(ndnr(ireac),ial)
              endif
              if(numel == 2)then
                if(nn(2,ndnr(ireac)) >= nfam(igas) .and. 
     &          nn(2,ndnr(ireac)) < nfam(igas+1))then
                  dest(igas,ial)=dest(igas,ial)+
     &            multip*chemrate(ndnr(ireac),ial)
                  dest(nn(2,ndnr(ireac)),ial)=
     &            dest(nn(2,ndnr(ireac)),ial) + 
     &            multip*chemrate(ndnr(ireac),ial)
                endif
              endif
            end do ! ial
          end do  ! i
        end if
      end do      ! igas

c Individual Species:

      nbeg=numfam+1
      nend=nfam(1)-1
      do igas=nbeg,nend
        dk=kdnr(igas+1)-kdnr(igas)
        if(dk >= 1) then
          do i=1,dk
            ireac=ireac+1
            do ial=1,maxl
              dest(igas,ial)=dest(igas,ial)+multip*
     &        chemrate(ndnr(ireac),ial)
            end do
          end do
        end if
      end do
      return
      end SUBROUTINE chem1



      SUBROUTINE chem1prn(kdnr,numel,nn,ndnr,chemrate,
     &                    index,multip,igas,total,maxl,I,J,jay)
!@sum chem1prn for printing out the chemical reactions
!@auth Drew Shindell (modelEifications by Greg Faluvegi)

C**** GLOBAL parameters and variables:

      USE DOMAIN_DECOMP_ATM, only : write_parallel
      USE TRCHEM_Shindell_COM, only: ay, lprn, nfam, p_4, numfam, y,
     &                              p_2, p_3

      IMPLICIT NONE

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
      INTEGER, DIMENSION(p_3)        :: ndnr
      INTEGER, DIMENSION(numel,p_2)  :: nn
      INTEGER, DIMENSION(p_4)        :: kdnr      
      INTEGER                        :: ireac
      character*17                   :: label
      character(len=300)             :: out_line
      logical                        :: jay
      REAL*8                         :: total,per
      REAL*8, DIMENSION(p_2,maxl)    :: chemrate

c FAMILIES ONLY:

      if(igas <= numfam) then
        if(kdnr(igas+1)-kdnr(igas) < 1) return
        do ireac=kdnr(igas),kdnr(igas+1)-1
          if(index <= 2)then 
            label=' chem reaction # '
          else
            label=' phot reaction # '
          endif
          if(nn(1,ndnr(ireac)) >= nfam(igas) .and. 
     &    nn(1,ndnr(ireac)) < nfam(igas+1))then
            per=0.d0
            if(y(igas,lprn) /= 0.d0) per=multip*100.d0*
     &      chemrate(ndnr(ireac),lprn)/y(igas,lprn)
            write(out_line,177) label,ndnr(ireac),' percent change'
     &      //' from ',ay(nn(1,ndnr(ireac))),' = ',per,
     &      ' dy=',multip*chemrate(ndnr(ireac),lprn)
            call write_parallel(trim(out_line),crit=jay)
            total=total+per
          endif
          if(numel == 2)then
            if(nn(2,ndnr(ireac)) >= nfam(igas) .and. 
     &      nn(2,ndnr(ireac)) < nfam(igas+1))then
              per=0.d0
              if(y(igas,lprn) /= 0.d0) per=multip*100.d0*
     &        chemrate(ndnr(ireac),lprn)/y(igas,lprn)
              write(out_line,177) label,ndnr(ireac),' percent change'
     &        //' from ',ay(nn(2,ndnr(ireac))),' = ',per,
     &        ' dy=',multip*chemrate(ndnr(ireac),lprn)
              call write_parallel(trim(out_line),crit=jay)
              total=total+per
            endif
          endif  
        end do 
        return
      end if

c INDIVIDUAL SPECIES:

      if(kdnr(igas+1)-kdnr(igas) < 1)return
      do ireac=kdnr(igas),kdnr(igas+1)-1
        if(index <= 2) then
          label=' chem reaction # '
        else
          label=' phot reaction # '
        endif
c       skip same reaction if written twice:
        if ((ireac > 1) .and. (ndnr(ireac) == ndnr(ireac-1))) CYCLE
        if(nn(1,ndnr(ireac)) == igas)then
          per=0.d0
          if(y(igas,lprn) /= 0.d0) per=100.d0*multip*
     &    chemrate(ndnr(ireac),lprn)/y(igas,lprn)
          write(out_line,106) label,ndnr(ireac),' percent change = '
     &    ,per,' dy=',multip*chemrate(ndnr(ireac),lprn)
          call write_parallel(trim(out_line),crit=jay)
          total=total+per
        endif
        if(numel == 2)then
          if(nn(2,ndnr(ireac)) == igas)then
            per=0.d0
            if(y(igas,lprn) /= 0.d0) per=100.d0*multip*
     &      chemrate(ndnr(ireac),lprn)/y(igas,lprn)
            write(out_line,106) label,ndnr(ireac),' percent change = '
     &      ,per,' dy=',multip*chemrate(ndnr(ireac),lprn)
            call write_parallel(trim(out_line),crit=jay)
            total=total+per
          endif
        endif 
      end do
 106  format(a17,i3,a18,f10.0,a4,e12.3)
 177  format(a17,i3,a21,a8,a3,f10.0,a4,e12.3)
  
      return
      end SUBROUTINE chem1prn
