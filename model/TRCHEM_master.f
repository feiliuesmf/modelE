#include "rundeck_opts.h"
      SUBROUTINE masterchem
!@sum masterchem main chemistry routine
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on masterchem000_M23p)   
!@calls photoj,checktracer,Crates,Oxinit,HOxfam,NOxfam,chemstep
C
C PLEASE SEE THE WARNING ABOUT CHANGEL VARIABLE IN THE
C STRATOSPHERIC OVERWRITE SECTION.
c
C**** GLOBAL parameters and variables:
c
      USE MODEL_COM, only: q,JDAY,IM,JM,sig,ptop,psf,byim,ls1,
     & COUPLED_CHEM
      USE DYNAMICS, only: pedn,LTROPO
      USE RADNCB, only : COSZ1,salbfj=>salb,rcloudfj=>rcld,O3_rad_save
      USE GEOM, only : BYDXYP, DXYP, LAT_DG, IMAXJ
      USE FLUXES, only : tr3Dsource
      USE TRACER_COM, only: n_Ox,n_NOx,n_N2O5,n_HNO3,n_H2O2,n_CH3OOH,
     &                  n_HCHO,n_HO2NO2,n_CO,n_CH4,n_PAN,n_Isoprene,
     &                  n_AlkylNit,n_Alkenes,n_Paraffin,ntm_chem,
     &                  n_DMS, n_MSA, n_SO2,n_SO4,n_H2O2_s,
     &                  oh_live,no3_live,nChemistry,nStratwrite
#ifdef SHINDELL_STRAT_CHEM
     &                  ,n_HBr,n_HOCl,n_HCl,n_ClONO2,n_ClOx,n_BrOx,
     &                  n_BrONO2,n_CFC,n_HOBR
#endif
#ifdef regional_Ox_tracers
     &                  ,NregOx,n_OxREG1
#endif
      USE CONSTANT, only: radian,gasc,mair,mb2kg,pi
      USE TRACER_DIAG_COM, only : jls_N2O5sulf,tajls
      USE TRCHEM_Shindell_COM
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@param by35 1/35 used for spherical geometry constant
!@param JN J around 30 N
!@param JS J around 30 S
!@param nlast either ntm_chem or ntm_chem-NregOx for chemistry loops
      INTEGER, PARAMETER :: JS = JM/3 + 1, JN = 2*JM/3
      INTEGER, PARAMETER :: nlast =
#ifdef regional_Ox_tracers
     &                              ntm_chem-NregOx
#else
     &                              ntm_chem
#endif
      REAL*8, PARAMETER  :: by35=1.d0/35.d0
!@var FASTJ_PFACT temp factor for vertical pressure-weighting
!@var FACT1 temp variable for start overwrite
!@var bydtsrc reciprocal of the timestep dtsrc
!@var local logical for error checking 
!@var byam75 the reciprocal air mass near 75 hPa level
!@var average tropospheric ch4 value near 569 hPa level
!@var PRES2 local nominal pressure for verticle interpolations
!@var thick thickness of each layer in km
!@var tempO2, photO2 and sHerzberg for O2->Ox diagnostic
!@var ClOx_old total ClOx at start of chemical timestep
!@var ClTOT total chlorine in all forms (reactive and reservoir)
!@var colmO2 is overhead oxygen column
!@var CH4FACT, FACTJ, r179 for setting CH4 ICs and strat distribution
!@var changeClONO2,changeClOx,changeHOCl,changeHCl nighttime changes
!@var changehetClONO2 nighttime het change in ClONO2 (on sulfate)
!@var pscX flag for psc existence
!@var chg106,chg107,chg108,chg109 reaction rates for het rxns on pscs
!@var rmrClOx,rmrBrOx dummy vars with mixing ratios of halogens
!@var rmv dummy variable for halogne removal in trop vs height
!@var airvol Air volume in grid cell
!@var changeL 2D array holds the local change due to chem or strat 
!@+   overwrite until adding to tr3Dsource (replaces 4D "change")
      REAL*8, DIMENSION(LM,NTM) :: changeL
      REAL*8, DIMENSION(LM) :: airvol
      REAL*8 :: CH4FACT,FACTj,r179
      REAL*8, DIMENSION(LM) :: PRES2
      REAL*8 FASTJ_PFACT, FACT1, bydtsrc, byam75
#ifdef SHINDELL_STRAT_CHEM
      REAL*8, DIMENSION(LM)      :: ClOx_old 
      REAL*8, DIMENSION(LM), PARAMETER      :: thick = (/
     & 0.3d0,0.36d0,0.44d0,0.65d0,1.2d0,1.6d0,2.0d0,2.0d0,1.6d0,
     & 1.6d0,1.5d0,1.6d0,1.9d0,2.5d0,3.7d0,3.6d0,4.0d0,5.2d0,
     & 7.4d0,8.6d0,8.6d0,9.2d0,12.4d0/)
      REAL*8, DIMENSION(JM)      :: tempO2 
      REAL*8, DIMENSION(JM,LM)   :: photO2 
      REAL*8 sHerzberg,CLTOT,colmO2,changeClONO2,changeClOx,
     *changeHOCl,changeHCl,changehetClONO2,pscX,chg106,chg107,
     *chg108,chg109,rmrClOx,rmrBrOx,rmv,rmrOx
#endif
      REAL*8 :: CH4_569
      LOGICAL error
!@var I,J,L,N,igas,inss,LL,Lqq,JJ,J3,L2,n2 dummy loop variables
      INTEGER igas,LL,I,J,L,N,inss,Lqq,JJ,J3,L2,n2,Jqq,Iqq
!@var imonth,m only needed to choose Ox strat correction factors
      INTEGER imonth,m
#ifdef regional_Ox_tracers
!@var sumOx for summing regional Ox tracers
!@var bysumOx reciprocal of sum of regional Ox tracers
      REAL*8 sumOx, bysumOx
#endif
c
C++++ First, some INITIALIZATIONS :
      bydtsrc = 1.d0/dtsrc
      BYFJM=1.d0/float(JM)
      PRES2(:)=SIG(:)*(PSF-PTOP)+PTOP
C
c Set "chemical time step". Really this is a method of applying only
c a fraction of the chemistry change to the tracer mass for the first
c 30 hours.  That fraction is: dt2/dtscr.  E.g. in the first hour it
c is (dtsrc/24)/dtsrc = 1/24th of the chemistry change is applied.
c This is to work around initial instabilities.
c
      if(Itime-ItimeI.le.3)then
        dt2=dtsrc/24.d0          ! 150.
      elseif(Itime-ItimeI.gt.3.and.Itime-ItimeI.le.6)then
        dt2=dtsrc/12.d0          ! 300.
      elseif(Itime-ItimeI.gt.6.and.Itime-ItimeI.le.11)then
        dt2=dtsrc/6.d0           ! 600.
      elseif(Itime-ItimeI.gt.11.and.Itime-ItimeI.le.30)then
        dt2=dtsrc/2.4d0          ! 1500.
      elseif(Itime-ItimeI.gt.30)then
        dt2=dtsrc                ! 3600
      endif
c
c     Calculate new photolysis rates every n_phot main timesteps
      MODPHOT=MOD(Itime-ItimeI,n_phot)
C****
C**** CALCULATE TX, THE REAL TEMPERATURE
C**** (note this section is already done in DIAG.f)
      DO L=1,LM
        TX(1,1,L)=T(1,1,L)*PK(L,1,1)
        TX(1,JM,L)=T(1,JM,L)*PK(L,1,JM)
        DO I=2,IM
          TX(I,1,L)=TX(1,1,L)
          TX(I,JM,L)=TX(1,JM,L)
        END DO
        DO J=2,JM-1
        DO I=1,IM
          TX(I,J,L)=T(I,J,L)*PK(L,I,J)
        END DO
        END DO
      END DO
C
!$OMP  PARALLEL DO PRIVATE (airvol, changeL, FASTJ_PFACT, 
#ifdef regional_Ox_tracers
!$OMP* bysumOx, sumOx,
#endif
#ifdef SHINDELL_STRAT_CHEM
!$OMP* CLTOT, ClOx_old, COLMO2, changeClONO2, changeClOx,
!$OMP* changehetClONO2, changeHOCl, changeHCl,
!$OMP* chg106, chg107, chg108, chg109,
!$OMP* pscx, rmrclox, rmrbrox, rmrox, rmv,
#endif
#ifdef regional_Ox_tracers
!$OMP* n2,
#endif
!$OMP* LL, I, igas, inss, J, L, Lqq, N, error )
c
      DO J=1,JM                          ! >>>> MAIN J LOOP BEGINS <<<<
#ifdef SHINDELL_STRAT_CHEM
      DU_O3(J)=0.d0
#endif
      DO I=1,IMAXJ(J)                    ! >>>> MAIN I LOOP BEGINS <<<<
C
      DO L=1,LM
#ifdef regional_Ox_tracers
C As per Volker:
C make sure the regional Ox tracers didn't diverge from total Ox:
       sumOx=0.d0
       bysumOx=0.d0
       do n2=n_OxREG1,ntm_chem
         sumOx=sumOx+trm(i,j,l,n2)
       end do
       sumOx=MAX(sumOx,1.d-9) ! minimum to prevent NaN's as per Drew
       bysumOx=1.d0/sumOx 
       do n2=n_OxREG1,ntm_chem
         trm(i,j,l,n2)=trm(i,j,l,n2)*trm(i,j,l,n_Ox)*bysumOx
       end do
#endif
c      initialize the 2D change variable:
       changeL(L,:)=0.d0
c      save presure and temperature in local arrays:
       pres(L)=PMID(L,I,J)
       ta(L)=TX(I,J,L)
c
c Air volume for sulfate SA calculation (m3)
       airvol(L)=(((DXYP(J)*AM(L,I,J)*1.d3)/mair)*gasc
     &      *ta(L))/pres(L)*100.d0 
c
c      calculate M and set fixed ratios for O2 & H2:
       y(nM,L)=pres(L)/(ta(L)*1.38d-19)
       y(nO2,L)=y(nM,L)*pfix_O2
#ifdef SHINDELL_STRAT_CHEM
       if(pres2(l).gt.20.d0)then
         y(nH2,L)=y(nM,L)*pfix_H2
       else
C was:   y(nH2,L)=y(nM,L)*pfix_H2*7.d1/(7.d1+L-LTROPO(I,J)+1)
C        Now: a drop of 0.3 molec/cm3 per 12 hPa decrease:
         y(nH2,L)=y(nM,L)*pfix_H2 + 2.5d-2*(pres2(L)-20.d0)
       endif
       CLTOT=0.d0
#else
       y(nH2,L)=y(nM,L)*pfix_H2
#endif
c
c      Tracers (converted from mass to number density)
c
       do igas=1,nlast
         y(igas,L)=trm(I,J,L,igas)*y(nM,L)*mass2vol(igas)*
     *   BYDXYP(J)*BYAM(L,I,J)
       enddo
c
c      If desired, fix the methane concentration used in chemistry
       if(fix_CH4_chemistry.eq.1) THEN
         if(J.lt.JEQ)then ! SH
           y(n_CH4,L)=y(nM,L)*pfix_CH4_S
         else             ! NH
           y(n_CH4,L)=y(nM,L)*pfix_CH4_N
         endif
       end if
#ifdef SHINDELL_STRAT_CHEM
c
c Set CLTOT based on CFCs (3.3 ppbv yield from complete oxidation of
c 1.7 ppbv CFC plus 0.5 ppbv background
      if(L.ge.LTROPO(I,J)+1)then
        CLTOT=((y(n_CFC,1)/y(nM,1)-y(n_CFC,L)/y(nM,L))*(3.3d0/1.8d0)*
     *  y(n_CFC,1)/(1.8d-9*y(nM,1)))
        CLTOT=CLTOT+0.5d-9
        CLTOT=CLTOT*y(nM,L)/
     *  (y(n_ClOx,L)+y(n_HCl,L)+y(n_HOCl,L)+y(n_ClONO2,L))
        y(n_ClOx,L)=y(n_ClOx,L)*CLTOT
        changeL(L,n_ClOx)=trm(I,J,L,n_ClOx)*(CLTOT-1.D0)
        y(n_HCl,L)=y(n_HCl,L)*CLTOT
        changeL(L,n_HCl)=trm(I,J,L,n_HCl)*(CLTOT-1.D0)
        y(n_HOCl,L)=y(n_HOCl,L)*CLTOT
        changeL(L,n_HOCl)=trm(I,J,L,n_HOCl)*(CLTOT-1.D0)
        y(n_ClONO2,L)=y(n_ClONO2,L)*CLTOT
        changeL(L,n_ClONO2)=trm(I,J,L,n_ClONO2)*(CLTOT-1.D0)
      endif
c     Save initial ClOx amount for use in ClOxfam
      ClOx_old(L)=trm(I,J,L,n_ClOx)*y(nM,L)*mass2vol(n_ClOx)*
     *  BYDXYP(J)*BYAM(L,I,J)
#endif
c
c Limit N2O5 number density:
       if(y(n_N2O5,L).lt.1.d0)y(n_N2O5,L)=1.d0
c Set H2O, based on Q:
       y(nH2O,L)=Q(I,J,L)*MWabyMWw*y(nM,L)
c
       if(Itime.eq.ItimeI)then !initial startup
        y(nAldehyde,L)=y(nM,L)*pfix_Aldehyde
       else
        y(nAldehyde,L)=yAldehyde(I,J,L)
       endif
c      set [NO]=0 for first HOx calc, NO2 = NOx
c      set reactive species for use in family chemistry & nighttime NO2
       y(nNO2,L)     =y(n_NOx,L)*pNOx(I,J,L)
       y(nNO,L)      =y(n_NOx,L)*(1.-pNOx(I,J,L))
       y(nO3,L)      =pOx(I,J,L)*y(n_Ox,L)
       y(nCH3O2,L)   =yCH3O2(I,J,L)
       y(nC2O3,L)    =yC2O3(I,J,L)
       y(nXO2,L)     =yXO2(I,J,L)
       y(nXO2N,L)    =yXO2N(I,J,L)
       y(nRXPAR,L)   =yRXPAR(I,J,L)
       y(nROR,L)     =yROR(I,J,L)
#ifdef SHINDELL_STRAT_CHEM
       y(nOClO,L)=y(n_ClOx,L)*pOClOx(I,J,L)
       y(nClO,L)=y(n_ClOx,L)*pClOx(I,J,L)
       y(nCl,L)=y(n_ClOx,L)*pClx(I,J,L)
       y(nBr,L)=y(n_BrOx,L)*(1.d0-pBrOx(I,J,L))
       y(nBrO,L)=y(n_BrOx,L)*pBrOx(I,J,L)
#endif
      END DO ! L
c
C For solar zenith angle, we now use the arccosine of the COSZ1
C from the radiation code, which is the cosine of the solar zenith 
C angle averaged over the physics time step.
C If the solar zenith angle (sza) from the radiation code is > 90 deg,
C (and hence COSZ1 is set to 0), recalculate it with get_sza routine:
      IF(COSZ1(I,J).eq.0.d0) THEN
        call get_sza(I,J,sza)
      ELSE
        sza = acos(COSZ1(I,J))*byradian
      END IF
c
c     update radiation and temperatures, call PHOTOLYSIS every
c     [desired number] of hours:
      if(MODPHOT.eq.0)then             !  >>>> PHOTOLYSIS IF BEGIN <<<<
c
c      additional SUNLIGHT criterion (see also fam chem criterion):
       if((SALBFJ(I,J).ne.0.d0).AND.(sza.lt.szamax))then!>>SUNLIGHT<<<
c
c       define column temperatures to be sent to FASTJ:
        TFASTJ = ta
c
#ifdef TRACERS_SPECIAL_Shindell
#ifdef SHINDELL_STRAT_CHEM
c       define pressures to be sent to FASTJ (centers):
        DO LL=1,LM
          PFASTJ2(LL) = PMID(LL,I,J)
        END DO

      PFASTJ2(LM+1)=0.00206d0  ! P at SIGE(LM+1)
      PFASTJ2(LM+2)=0.00058d0  ! unknown, so fudged for now.
      PFASTJ2(LM+3)=0.00028d0  ! unknown, so fudged for now.
#else
c       define pressures to be sent to FASTJ (centers):
        DO LL=2,2*LM,2
          PFASTJ(LL) = PMID(LL/2,I,J)
        END DO
c       define pressures to be sent to FASTJ (edges):
        PFASTJ(1)=PEDN(1,I,J)
        DO LL=3,(2*LM)+1,2
          PFASTJ(LL) = PEDN((LL+1)/2,I,J)
        END DO
        PFASTJ((2*LM)+2) = 0.1d0*PFASTJ((2*LM)+1)
C       This is a fudge, so that we don't have to get mesosphere data:
        PFASTJ((2*LM)+2) = 0.00058d0 ! for 23 layer model...
#endif
#endif
c
#ifdef SHINDELL_STRAT_CHEM
c Pass O3 array (in ppmv) to fastj.
c Above these levels fastj uses climatological (Nagatani) O3
c read in by chem_init. lg.feb99., gsf.apr01. dts.aug.02
        DO LL=1,LM
          O3_FASTJ(LL)=y(nO3,LL)/y(nM,LL)
        END DO
#else
c Interpolate O3 (in ppm) from bottom LTROPO(I,J) model sigma levels
c (ie those levels normally in troposphere) onto bottom 2*LTROPO(I,J)
c FASTJ levels. Above these levels fastj uses climatological (Nagatani)
c O3 read in by chem_init.
c
c       define O3 to be sent to FASTJ (centers):
        DO LL=2,2*LTROPO(I,J),2
          O3_FASTJ(LL)=y(nO3,LL/2)
        ENDDO
c       define O3 to be sent to FASTJ (edges):
        O3_FASTJ(1)=y(nO3,1)*O3_1_fact ! see parameter declaration...
        DO LL=3,2*(LTROPO(I,J)-1)-1,2
c         interpolation factor, based on pressure:
          FASTJ_PFACT=(PFASTJ(LL)-PFASTJ(LL-1))/
     &    (PFASTJ(LL+1)-PFASTJ(LL-1))
          O3_FASTJ(LL)=y(nO3,(LL-1)/2)+
     &    (y(nO3,((LL-1)/2)+1)-y(nO3,(LL-1)/2))*FASTJ_PFACT
c         lower limit on O3:
          IF(O3_FASTJ(LL).LT.0.d0) O3_FASTJ(LL)=-O3_FASTJ(LL)
        ENDDO
c       Convert to ppm (units used in fastj)
        DO LL=1,2*LTROPO(I,J)
          O3_FASTJ(LL)=O3_FASTJ(LL)/(PFASTJ(LL)*2.55d10)
        ENDDO
#endif
C
C CALL THE PHOTOLYSIS SCHEME:
C
        call photoj(I,J)
C
C And fill in the photolysis coefficients: ZJ --> ss:
#ifdef SHINDELL_STRAT_CHEM
        colmO2=5.6d20 
#endif
        DO L=JPNL,1,-1
          do inss=1,JPPJ
           ss(inss,L,I,J)=zj(L,inss)*
     &     by35*SQRT(1.224d3*(cos(ABS(LAT_DG(J,1))*radian))**2.+1.d0)
          enddo
#ifdef SHINDELL_STRAT_CHEM
          colmO2=colmO2+y(nO2,L)*thick(L)*1.d5
          if(L.ge.LTROPO(I,J)+1)then
            SF3(I,J,L)=1.3d-6*EXP(-1.d-7*colmO2**.35)
          else
            SF3(I,J,L)=0.d0
          endif
          SF3(I,J,L)=SF3(I,J,L)*
     &    by35*SQRT(1.224d3*(cos(ABS(LAT_DG(J,1))*radian))**2.+1.d0)
#endif
        END DO

       endif                               ! >>>> SUNLIGHT IF END <<<<
      endif                             !  >>>> PHOTOLYSIS IF END <<<<
c
c     calculate the chemical reaction rates:
      call Crates (I,J)      
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
cc    Main chemistry calculations    cc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c     Family partitioning:
c
      if((SALBFJ(I,J).ne.0.d0).AND.(sza.lt.szamax))then ! --SUNLIGHT--
#ifdef SHINDELL_STRAT_CHEM
       call Oxinit(LM,I,J)
       call HOxfam(LM,I,J)
       call NOxfam(LM,I,J)
       call BrOxfam(LM,I,J)
#else
       call Oxinit(LTROPO(I,J),I,J)
       call HOxfam(LTROPO(I,J),I,J)
       call NOxfam(LTROPO(I,J),I,J)
#endif
c
cc    Non-family chemistry:
c
      call chemstep(I,J,changeL)
c
C Accumulate 3D radical arrays to pass to aerosol code
      if(coupled_chem.eq.1) then
        do l=1,LTROPO(I,J)
          oh_live(i,j,l)=y(nOH,L)
          no3_live(i,j,l)=yNO3(i,j,l)
        end do
      endif

#ifdef SHINDELL_STRAT_CHEM
      call ClOxfam(LM,I,J,ClOx_old)
#endif
c
C Some Chemistry Diagnostics:
      if(prnchg.and.J.eq.jprn.and.I.eq.iprn)then
       l=lprn
       write(6,*) ' '
       write(6,*) 'Family ratios at I,J,L: ',i,j,l
       write(6,*) 'OH/HO2 = ',y(nOH,l)/y(nHO2,l)
       write(6,*) 'O/O3 = ',y(nO,l)/y(nO3,l)
       write(6,*) 'O1D/O3 = ',y(nO1D,l)/y(nO3,l),
     *  '  J(O1D) = ',ss(2,l,I,J)
c     *  ,' J(O3) = ',ss(3,l,I,J)
       write(6,*) 'NO/NO2 = ',y(nNO,l)/y(nNO2,l),
     *  '   J(NO2) = ',ss(1,l,I,J)
c       write(6,*) 'pHOx, pNOx, pOx = ',pHOx(I,J,l),
c     * pNOx(I,J,l),pOx(I,J,l)
       write(6,*) 'conc OH = ',y(nOH,l)
#ifdef SHINDELL_STRAT_CHEM
       write(*,*) 'Cl,ClO,Cl2O2,OClO,Cl2 = ',y(nCl,l),
     *  y(nClO,l),y(nCl2O2,l),y(nOClO,l),y(nCl2,l)
       write(*,*) 'Br,BrO = ',y(nBr,l),y(nBrO,l)
       write(*,*) 'pCl,pClO,pOClO,pBrO = ',pClx(I,J,l),
     *  pClOx(I,J,l),pOClOx(I,J,l),pBrOx(I,J,l)
#endif
       write(6,*)'sun, SALBFJ,sza,I,J,Itime= ',SALBFJ(I,J),sza,I,J,Itime
       do inss=1,JPPJ
        write(6,195) ' J',inss,ay(ks(inss)),' = ',
     *   (ss(inss,Lqq,I,J),Lqq=1,11)
       enddo
       write(6,196) ' RCloud',(RCLOUDFJ(Lqq,I,J),Lqq=1,11)
       write(6,196) ' Ozone ',(y(nO3,Lqq),Lqq=1,11)
       write(6,*) ' '
      endif
 195  format(a2,i2,1x,a8,a3,11(1x,e9.2))
 196  format(a7,9x,11(1x,e9.2))
c
      else                           !             D A R K N E S S
C
C*****************************************************************
C      >> N2O5 sink on sulfate aerosol :: <<
C      REACTION PROBABLITY FORMULATION
C
C      To evaluate 'k' for N2O5 + H2O -> 2HNO3,
C       assume k = GAMMA.SA.v / 4 (Dentener and Crutzen, 1993)
C
C      GAMMA = 0.1, SA = given, v = molecular velocity (cm/s)
C      where v  = SQRT (8.Kb.T / PI*M); Kb = 1.38062E-23;
C      T = Temp (K); M = mass N2O5 (Kg)
C
C      Off-line sulfate fields to run in uncoupled mode 
C      give SO4 in cm2/cm3 'surface area density'.
C
C      On-line sulfate in coupled mode must be converted to cm2/cm3
C      Assume a monodispersed aerosol with radius 0.078 microns
C
C*****************************************************************
c
       do L=1,LTROPO(I,J) ! (troposphere)

         if (coupled_chem.eq.1) then
C Convert units of kg sulfate to cm2 sulfate /cm3 air
C Number of particles * mean surface area
C Mean surface area of sulfate particle = 7.6E-10 cm2
C Total surface area per grid box (cm2)
           sulfate(I,J,L)=(trm(I,J,L,n_SO4)*1.d3)/tr_mm(n_SO4)
     &     *6.022d23*7.6d-10
C Divide by grid box volume (cm3)
           sulfate(I,J,L)=sulfate(I,J,L)/(airvol(L)*1.d6)
         endif

         pfactor=dxyp(J)*AM(L,I,J)/y(nM,L)
         bypfactor=1.D0/pfactor
         RVELN2O5=SQRT(TX(I,J,L)*RKBYPIM)*100.d0
C        Calculate sulfate sink, and cap it at 90% of N2O5:
         wprod_sulf=
     &   DT2*sulfate(I,J,L)*y(n_N2O5,L)*RGAMMASULF*RVELN2O5*0.25d0
         if(wprod_sulf.gt.0.9d0*y(n_N2O5,L))wprod_sulf=0.9d0*y(n_N2O5,L)
         prod_sulf=wprod_sulf*pfactor
         TAJLS(J,L,jls_N2O5sulf)=TAJLS(J,L,jls_N2O5sulf)
     &     -(prod_sulf*bymass2vol(n_N2O5))
C
C*****************************************************************
c        g signifies gas phase
c        while prod_sulf and wprod_sulf are sulfate rxn
c        wprods are in molecules/cm3/s
c        prods/mass2vol are in mass units to add to tracers
c
c        NO3 amounts are a function of reaction 7 (NO2 + O3 -> NO3),
c        24, 25 (leave out 28, 0.9*32, 46, 52 outside NOx family)
c        NO2, similarly leave out 29, 45, and 46.
c        Keep NOx unchanged as this is only intrafamily
C*****************************************************************
C
c        define O3 from Ox:
         y(nO3,L)=y(n_Ox,L)*pOx(I,J,L)
c
c        calculate change in NO3:
         dNO3=rr(7,L)*y(nNO2,L)*y(n_Ox,L)-(rr(24,L)*y(nNO2,L)+
     &        2*rr(25,L)*yNO3(I,J,L))*yNO3(I,J,L) -(rr(36,L)
     &        *y(n_Alkenes,L)+rr(32,L)*y(n_Isoprene,L))*yNO3(I,J,L)
#ifdef SHINDELL_STRAT_CHEM
         dNO3=dNO3-(rr(28,L)*y(n_HCHO,L)+rr(99,L)*y(nNO2,L))
     &        *yNO3(I,J,L)+rr(92,L)*y(n_N2O5,L)
#else
         dNO3=dNO3-(rr(28,L)*y(n_HCHO,L)+rr(52,L)*y(nNO2,L))
     &        *yNO3(I,J,L)+rr(46,L)*y(n_N2O5,L)
#endif
         dNO3=dNO3*dt2!portions of dNO3 reflect fixes by dts 12/19/01
c
c        limit the change in NO3:
         if(-dNO3.gt.0.66d0*yNO3(I,J,L))dNO3=-0.66d0*yNO3(I,J,L)
c
c        apply the NO3 change limit the value to positive and 1/2 NOx:
         yNO3(I,J,L)=yNO3(I,J,L)+dNO3
         if(yNO3(I,J,L).lt.0.d0)yNO3(I,J,L)=0.d0
         if(yNO3(I,J,L).gt.y(n_NOx,L)*0.5d0)
     &   yNO3(I,J,L)=y(n_NOx,L)*0.5d0
c
c        calculate and limit NO2
         y(nNO2,L)=y(n_NOx,L)-yNO3(I,J,L)
         pNOx(I,J,L)=y(nNO2,L)/(y(n_NOx,L)+1.d-10)
         if(pNOx(I,J,L).gt.1.d0)pNOx(I,J,L)=1.d0
         if(pNOx(I,J,L).lt.0.5d0)pNOx(I,J,L)=0.5d0
C
C        LOWER LIMIT ON N2O5:
         if(y(n_N2O5,L).le.1.d0)y(n_N2O5,L)=1.d0
C
C Calculate and limit gaseous changes to HNO3, HCHO, N2O5, Aldehyde,
C Alkenes, Isoprene, and AlkylNit:
C
         gwprodHNO3=(y(n_HCHO,L)*rr(28,L)+2.5d-15*
     &              yAldehyde(I,J,L))*yNO3(I,J,L)*dt2
         if(gwprodHNO3.gt.0.5d0*y(n_NOx,L))gwprodHNO3=0.5d0*y(n_NOx,L)
         if(gwprodHNO3.gt.y(n_HCHO,L))gwprodHNO3=y(n_HCHO,L)
         gprodHNO3=gwprodHNO3*pfactor
C
#ifdef SHINDELL_STRAT_CHEM
         gwprodN2O5=(yNO3(I,J,L)*y(nNO2,L)*rr(99,L)-y(n_N2O5,L)
     &              *rr(92,L))*dt2
#else
         gwprodN2O5=(yNO3(I,J,L)*y(nNO2,L)*rr(52,L)-y(n_N2O5,L)
     &              *rr(46,L))*dt2
#endif
         if(gwprodN2O5.gt.0.25d0*y(n_NOx,L))gwprodN2O5=0.25d0*y(n_NOx,L)
         if(-gwprodN2O5.gt.0.5d0*y(n_N2O5,L))
     *        gwprodN2O5=-0.49d0*y(n_N2O5,L)
C
         changeAldehyde=(rr(36,L)*y(n_Alkenes,L)+rr(32,L)*
     &               y(n_Isoprene,L)*0.12d0-2.5d-15*yAldehyde(I,J,L))*
     &               yNO3(I,J,L)*dt2
         if(-changeAldehyde.gt.0.75d0*yAldehyde(I,J,L))changeAldehyde=
     &   -0.75d0*yAldehyde(I,J,L)
C
         changeAlkenes=(rr(32,L)*y(n_Isoprene,L)*0.45d0-rr(36,L)*
     &                y(n_Alkenes,L))*yNO3(I,J,L)*dt2+(rr(31,L)*
     &                y(n_Isoprene,L)-rr(35,L)*y(n_Alkenes,L))*
     &                y(nO3,L)*dt2
         if(-changeAlkenes.gt.0.75d0*y(n_Alkenes,L))changeAlkenes=
     &   -0.75d0*y(n_Alkenes,L)
C
         changeIsoprene=-(rr(32,L)*yNO3(I,J,L)
     &                  +rr(31,L)*y(nO3,L))*y(n_Isoprene,L)*dt2
         if(-changeIsoprene.gt.0.75d0*y(n_Isoprene,L))changeIsoprene=
     &   -0.75d0*y(n_Isoprene,L)
C
         changeHCHO=(rr(36,L)*y(n_Alkenes,L)+
     &              rr(32,L)*y(n_Isoprene,L)*0.03d0)*yNO3(I,J,L)*dt2
     &              -gwprodHNO3+(rr(31,L)*y(n_Isoprene,L)*0.9d0
     &              +rr(35,L)*y(n_Alkenes,L))*y(nO3,L)*0.64d0*dt2
C
         changeAlkylNit=rr(32,L)*y(n_Isoprene,L)*
     &   yNO3(I,J,L)*dt2*0.9d0
         if(-changeAlkylNit.gt.0.75d0*y(n_AlkylNit,L))changeAlkylNit=
     &   -0.75d0*y(n_AlkylNit,L)
C
c        Convert some changes to molecules/cm3/s:
         changeHNO3=gwprodHNO3+2*wprod_sulf  !always positive
         changeNOx=-gwprodHNO3-2*gwprodN2O5-(0.9d0*rr(32,L)*
     &    y(n_Isoprene,L)+2.5d-15*yAldehyde(I,J,L))*yNO3(I,J,L)*dt2
         if(-changeNOx.gt.y(n_NOx,L))changeNOx=-.95d0*y(n_NOx,L)!dts9/01
         changeN2O5=gwprodN2O5-wprod_sulf
c
c Ensure nitrogen conservation (presumably dNOx<0, others >0):
c
         rlossN=0.d0
         rprodN=0.d0
         if(changeNOx.lt.0.d0)then
          rlossN=rlossN+changeNOx
         else
          rprodN=rprodN+changeNOx
         endif
         if(changeHNO3.lt.0.d0)then
          rlossN=rlossN+changeHNO3
         else
          rprodN=rprodN+changeHNO3
         endif
         if(changeN2O5.lt.0.d0)then
          rlossN=rlossN+2*changeN2O5
         else
          rprodN=rprodN+2*changeN2O5
         endif
         if(changeAlkylNit.lt.0.d0)then
          rlossN=rlossN+changeAlkylNit
         else
          rprodN=rprodN+changeAlkylNit
         endif
         if(rprodN.gt.rlossN)then
           ratioN=-rlossN/rprodN
           if(changeNOx.gt.0.d0)     changeNOx    =changeNOx   *ratioN
           if(changeHNO3.gt.0.d0)    changeHNO3   =changeHNO3  *ratioN
           if(changeN2O5.gt.0.d0)    changeN2O5   =changeN2O5  *ratioN
           if(changeAlkylNit.gt.0.d0)changeAlkylNit=
     &                                           changeAlkylNit*ratioN
         else
           ratioN=rprodN/(-rlossN)
           if(changeNOx.lt.0.d0)     changeNOx    =changeNOx   *ratioN
           if(changeHNO3.lt.0.d0)    changeHNO3   =changeHNO3  *ratioN
           if(changeN2O5.lt.0.d0)    changeN2O5   =changeN2O5  *ratioN
           if(changeAlkylNit.lt.0.d0)changeAlkylNit=
     &                                           changeAlkylNit*ratioN
         endif
C
C Apply Alkenes, AlkyNit, and Aldehyde changes here:
C
         y(n_Alkenes,L)=y(n_Alkenes,L)+changeAlkenes
         y(n_AlkylNit,L)=y(n_AlkylNit,L)+changeAlkylNit
         yAldehyde(I,J,L)=yAldehyde(I,J,L)+changeAldehyde
C
C Note: there is a lower limit of 1 placed on the resulting tracer mass
C from the following changes. This is to prevent negative tracer mass:
C
C -- HCHO --
c        Gas phase NO3 + HCHO -> HNO3 + CO yield of HCHO & CO
         changeL(L,n_HCHO)=changeHCHO*pfactor*bymass2vol(n_HCHO)
         if(-changeL(L,n_HCHO).gt.trm(I,J,L,n_HCHO))then
           changeL(L,n_HCHO)=-.95d0*trm(I,J,L,n_HCHO)
           changeHCHO=changeL(L,n_HCHO)*mass2vol(n_HCHO)*bypfactor
         endif
         IF((trm(i,j,l,n_HCHO)+changeL(l,n_HCHO)).lt.1.d0) THEN
           changeL(l,n_HCHO) = 1.d0 - trm(i,j,l,n_HCHO)
           changeHCHO=changeL(L,n_HCHO)*mass2vol(n_HCHO)*bypfactor
         ENDIF
         wprodHCHO=changeHCHO
C -- CO --
         changeL(L,n_CO)=gprodHNO3*bymass2vol(n_CO)
         IF((trm(i,j,l,n_CO)+changeL(l,n_CO)).lt.1.d0)
     &   changeL(l,n_CO) = 1.d0 - trm(i,j,l,n_CO)
         wprodCO=gwprodHNO3   ! <<< note
C -- HNO3 --  (HNO3 from gas and het phase rxns )
         changeL(L,n_HNO3)=changeHNO3*pfactor*bymass2vol(n_HNO3)
         IF((trm(i,j,l,n_HNO3)+changeL(l,n_HNO3)).lt.1.d0) THEN
           changeL(l,n_HNO3) = 1.d0 - trm(i,j,l,n_HNO3)
           changeHNO3=changeL(L,n_HNO3)*mass2vol(n_HNO3)*bypfactor
         END IF
C -- N2O5 --  (N2O5 from gas and het phase rxns)
         changeL(L,n_N2O5)=changeN2O5*pfactor*bymass2vol(n_N2O5)
         IF((trm(i,j,l,n_N2O5)+changeL(l,n_N2O5)).lt.1.d0) THEN
           changeL(l,n_N2O5) = 1.d0 - trm(i,j,l,n_N2O5)
           changeN2O5=changeL(L,n_N2O5)*mass2vol(n_N2O5)*bypfactor
         END IF
c -- NOx --   (NOx from gas phase rxns)
         changeL(L,n_NOx)=changeNOx*pfactor*bymass2vol(n_NOx)
         IF((trm(i,j,l,n_NOx)+changeL(l,n_NOx)).lt.1.d0) THEN
           changeL(l,n_NOx) = 1.d0 - trm(i,j,l,n_NOx)
           changeNOx=changeL(L,n_NOx)*mass2vol(n_NOx)*bypfactor
         END IF
C -- Alkenes --  (Alkenes from gas phase rxns)
         changeL(L,n_Alkenes)=
     &   changeAlkenes*pfactor*bymass2vol(n_Alkenes)
         IF((trm(i,j,l,n_Alkenes)+changeL(l,n_Alkenes)).lt.1.d0)THEN
           changeL(l,n_Alkenes) = 1.d0 - trm(i,j,l,n_Alkenes)
           changeAlkenes=changeL(L,n_Alkenes)*mass2vol(n_Alkenes)
     &     *bypfactor
         END IF
c -- Isoprene -- (Isoprene from gas phase rxns)
         changeL(L,n_Isoprene)=
     &   changeIsoprene*pfactor*bymass2vol(n_Isoprene)
         IF((trm(i,j,l,n_Isoprene)+changeL(l,n_Isoprene)).lt.1.d0)
     &   THEN
           changeL(l,n_Isoprene) = 1.d0 - trm(i,j,l,n_Isoprene)
           changeIsoprene=changeL(L,n_Isoprene)*mass2vol(n_Isoprene)
     &     *bypfactor
         END IF
c -- AlkylNit -- (AlkylNit from gas phase rxns)
         changeL(L,n_AlkylNit)=
     &   changeAlkylNit*pfactor*bymass2vol(n_AlkylNit)
         IF((trm(i,j,l,n_AlkylNit)+changeL(l,n_AlkylNit)).lt.1.d0)
     &   THEN
           changeL(l,n_AlkylNit) = 1.d0 - trm(i,j,l,n_AlkylNit)
           changeAlkylNit=changeL(L,n_AlkylNit)*mass2vol(n_AlkylNit)
     &     *bypfactor
         END IF

C Accumulate 3D radical arrays to pass to aerosol code
C Make sure we get the nightime values
C Set OH to zero for now

         if(coupled_chem.eq.1) then
           oh_live(i,j,l)=0.d0
           no3_live(i,j,l)=yNO3(i,j,l)
         endif

C
c Some More Chemistry Diagnostics:
         if(prnchg.and.J.eq.jprn.and.I.eq.iprn.and.L.eq.lprn)then
          write(6,*) 'dark, SALBFJ,sza,I,J,L,Itime= ',SALBFJ(I,J),sza,I
     &          ,J,L,Itime
          write(6,198) ay(n_NOx),': ',
     *         changeNOx,' molecules produced; ',
     *      100.d0*(changeNOx)/y(n_NOx,L),' percent of'
     *     ,y(n_NOx,L),'(',1.d9*y(n_NOx,L)/y(nM,L),' ppbv)'
          write(6,198) ay(n_HNO3),': ',
     *         changeHNO3,' molecules produced; ',
     *     100.d0*(changeHNO3)/y(n_HNO3,L),' percent of'
     *     ,y(n_HNO3,L),'(',1.d9*y(n_HNO3,L)/y(nM,L),' ppbv)'
          write(6,198) ay(n_N2O5),': ',
     *         changeN2O5,' net molec produced; ',
     *      100.d0*(changeN2O5)/y(n_N2O5,L),' percent of'
     *     ,y(n_N2O5,L),'(',1.d9*y(n_N2O5,L)/y(nM,L),' ppbv)'
          write(6,198) ay(n_N2O5),': ',
     *            gwprodN2O5,' molec prod fm gas;  ',
     *      100.d0*(gwprodN2O5)/y(n_N2O5,L),' percent of'
     *     ,y(n_N2O5,L),'(',1.d9*y(n_N2O5,L)/y(nM,L),' ppbv)'
          write(6,198) ay(n_HCHO),': ',
     *         wprodHCHO,' molecules produced; ',
     *      100.d0*(wprodHCHO)/y(n_HCHO,L),' percent of'
     *     ,y(n_HCHO,L),'(',1.d9*y(n_HCHO,L)/y(nM,L),' ppbv)'
          write(6,198) 'Aldehyde',': ',
     *         changeAldehyde,' molecules produced; ',
     *      100.d0*(changeAldehyde)/yAldehyde(I,J,L),' percent of'
     *     ,yAldehyde(I,J,L),'(',1.d9*yAldehyde(I,J,L)/y(nM,L),' ppbv)'
          write(6,198) 'Alkenes ',': ',
     *         changeAlkenes,' molecules produced; ',
     *      100.d0*(changeAlkenes)/y(n_Alkenes,L),' percent of'
     *     ,y(n_Alkenes,L),'(',1.d9*y(n_Alkenes,L)/y(nM,L),' ppbv)'
          write(6,198) 'Isoprene',': ',
     *         changeIsoprene,' molecules produced; ',
     *      100.d0*(changeIsoprene)/y(n_Isoprene,L),' percent of'
     *     ,y(n_Isoprene,L),'(',1.d9*y(n_Isoprene,L)/y(nM,L),' ppbv)'
          write(6,198) 'AlkylNit',': ',
     *         changeAlkylNit,' molecules produced; ',
     *      100.d0*(changeAlkylNit)/y(n_AlkylNit,L),' percent of'
     *     ,y(n_AlkylNit,L),'(',1.d9*y(n_AlkylNit,L)/y(nM,L),' ppbv)'
          write(6,199) 'NO2, NO3  = ',y(nNO2,L),yNO3(I,J,L)
         endif
 198     format(1x,a8,a2,e13.3,a21,f10.0,a11,2x,e13.3,3x,a1,f12.5,a6)
 199     format(1x,a20,2(2x,e13.3))
C
C Make sure nighttime chemistry changes are not too big:
C
         error=.false.
         if(changeNOx.lt.-1.d15.OR.changeNOx.gt.1.d15) then
          write(*,*) 'Big chg@ Itime,I,J,L,NOx ',Itime,I,J,L,changeNOx
          error=.true.
         end if
         if(changeHNO3.lt.-1.d15.OR.changeHNO3.gt.1.d15) then
          write(*,*) 'Big chg@ Itime,I,J,L,HNO3',Itime,I,J,L,changeHNO3
          error=.true.
         end if
         if(changeN2O5.lt.-1.d15.OR.changeN2O5.gt.1.d15) then
          write(*,*) 'Big chg@ Itime,I,J,L,N2O5',Itime,I,J,L,changeN2O5
          error=.true.
         end if
         if(wprodHCHO.lt.-1.d15.OR.wprodHCHO.gt.1.d15) then
          write(*,*)'Big chg@ Itime,I,J,L,HCHO',Itime,I,J,L,wprodHCHO
          error=.true.
         endif
         if(error)call stop_model('nighttime chem: big changes',255)
C
       enddo  ! troposphere loop
C
#ifdef SHINDELL_STRAT_CHEM
C
cc     Nighttime stratospheric chemistry
       do L=LTROPO(I,J)+1,LM
c
         pfactor=dxyp(J)*AM(L,I,J)/y(nM,L)
         wprod_sulf=DT2*y(n_N2O5,L)*rr(105,L)
c        define O3 from Ox:
         y(nO3,L)=y(n_Ox,L)*pOx(I,J,L)
c
c        calculate change in NO3:
         dNO3=rr(7,L)*y(nNO2,L)*y(n_Ox,L)-(rr(24,L)*y(nNO2,L)+
     &        2.d0*rr(25,L)*yNO3(I,J,L))*yNO3(I,J,L)
         dNO3=dNO3-(rr(28,L)*y(n_HCHO,L)+rr(99,L)*y(nNO2,L))
     &        *yNO3(I,J,L)+rr(92,L)*y(n_N2O5,L)
         dNO3=dNO3*dt2
c
c        limit the change in NO3:
         if(-dNO3.gt.0.66d0*yNO3(I,J,L))dNO3=-0.66d0*yNO3(I,J,L)
c
c        apply the NO3 change limit the value to positive and 1/2 NOx:
         yNO3(I,J,L)=yNO3(I,J,L)+dNO3
         if(yNO3(I,J,L).lt.0.d0)yNO3(I,J,L)=0.d0
         if(yNO3(I,J,L).gt.y(n_NOx,L)*0.5d0)yNO3(I,J,L)=
     &    y(n_NOx,L)*0.5d0
c
c        calculate and limit NO2
         y(nNO2,L)=y(n_NOx,L)-yNO3(I,J,L)
         pNOx(I,J,L)=y(nNO2,L)/(y(n_NOx,L)+1.d-10)
         if(pNOx(I,J,L).gt.1.d0)pNOx(I,J,L)=1.d0
         if(pNOx(I,J,L).lt.0.5d0)pNOx(I,J,L)=0.5d0
C
C        LOWER LIMIT ON N2O5:
         if(y(n_N2O5,L).le.1.d0)y(n_N2O5,L)=1.d0
C
C Calculate and limit gaseous changes to HNO3, HCHO, N2O5, Aldehyde,
C Alkenes, Isoprene, and AlkylNit:
C
         gwprodHNO3=y(n_HCHO,L)*rr(28,L)*yNO3(I,J,L)*dt2
         if(gwprodHNO3.gt.0.5d0*y(n_NOx,L))gwprodHNO3=0.5d0*y(n_NOx,L)
         if(gwprodHNO3.gt.y(n_HCHO,L))gwprodHNO3=y(n_HCHO,L)
         gprodHNO3=gwprodHNO3*pfactor
         gwprodN2O5=(yNO3(I,J,L)*y(nNO2,L)*rr(99,L)-y(n_N2O5,L)
     &              *rr(92,L))*dt2
         if(gwprodN2O5.gt.0.25d0*y(n_NOx,L))gwprodN2O5=0.25d0*y(n_NOx,L)
         if(-gwprodN2O5.gt.0.5d0*y(n_N2O5,L))
     *        gwprodN2O5=-0.49d0*y(n_N2O5,L)
         changeClONO2=y(nClO,L)*rr(103,L)*y(nNO2,L)*dt2
         if(changeClONO2.ge.y(nClO,L))changeClONO2=0.8d0*y(nClO,L)
         if(changeClONO2.ge.y(nNO2,L))changeClONO2=0.8d0*y(nNO2,L)
         changeClOx=-changeClONO2
c
c        Convert some changes to molecules/cm3/s:(w/o ClONO2->HNO3 het or PSC)
         changeHNO3=gwprodHNO3+2.d0*wprod_sulf  !always positive
         changeNOx=-gwprodHNO3-2.d0*gwprodN2O5
         if(-changeNOx.gt.y(n_NOx,L))changeNOx=-.95d0*y(n_NOx,L)

         changeN2O5=gwprodN2O5-wprod_sulf
c
c Ensure nitrogen conservation (presumably dNOx<0, others >0):
c
         rlossN=0.d0
         rprodN=0.d0
         if(changeNOx.lt.0.d0)then
          rlossN=rlossN+changeNOx
         else
          rprodN=rprodN+changeNOx
         endif
         if(changeHNO3.lt.0.d0)then
          rlossN=rlossN+changeHNO3
         else
          rprodN=rprodN+changeHNO3
         endif
         if(changeN2O5.lt.0.d0)then
          rlossN=rlossN+2.d0*changeN2O5
         else
          rprodN=rprodN+2.d0*changeN2O5
         endif
         if(rprodN.gt.rlossN)then
          ratioN=-rlossN/rprodN
          if(changeNOx.gt.0.d0)     changeNOx     =changeNOx     *ratioN
          if(changeHNO3.gt.0.d0)    changeHNO3    =changeHNO3    *ratioN
          if(changeN2O5.gt.0.d0)    changeN2O5    =changeN2O5    *ratioN
         else
          ratioN=rprodN/(-rlossN)
          if(changeNOx.lt.0.d0)     changeNOx     =changeNOx     *ratioN
          if(changeHNO3.lt.0.d0)    changeHNO3    =changeHNO3    *ratioN
          if(changeN2O5.lt.0.d0)    changeN2O5    =changeN2O5    *ratioN
         endif
c
         changeNOx=changeNOx-changeClONO2
         if(-changeNOx.gt.y(n_NOx,L))changeNOx=-.95d0*y(n_NOx,L)

cc     Het rxn ClONO2+H2O on sulfate only if no PSCs
         if(ta(l).gt.198.d0.and.rr(106,L).gt.2.d-35)then
          changehetClONO2=-(rr(106,L)*y(n_ClONO2,L))*dt2
          if(changehetClONO2.ge.y(n_ClONO2,L))changehetClONO2=
     &     -0.8d0*y(n_ClONO2,L)
          changeClONO2=changeClONO2+changehetClONO2
          changeHOCl=-changehetClONO2
          changeHNO3=changeHNO3-changehetClONO2
         else
          changeHOCl=0.d0
         endif

         changeHCl=0.d0

cc     PSC chemistry
c 105 N2O5    +H2O     -->HNO3    +HNO3  (calculated above)
c 106 ClONO2  +H2O     -->HOCl    +HNO3
c 107 ClONO2  +HCl     -->Cl      +HNO3        !really makes Cl2
c 108 HOCl    +HCl     -->Cl      +H2O         !raeally makes Cl2
c 109 N2O5    +HCl     -->Cl      +HNO3        !really makes ClNO2  2

         pscX=0.d0
         if(pres2(l).le.150.d0.and.pres2(l).ge.35.d0)then
         if((J.LE.JS).OR.(J.GT.JN))then !not in tropics
         if(ta(l).le.198.d0)then
           pscX=1.d0
           chg106=rr(106,L)*y(n_ClONO2,L)*y(nH2O,L)*dt2 !inlc H2O?
           if(chg106.ge.0.4d0*y(n_ClONO2,L))chg106=0.4d0*y(n_ClONO2,L)
           chg107=rr(107,L)*y(n_ClONO2,L)*y(n_HCl,L)*dt2
           if(chg107.ge.0.4d0*y(n_ClONO2,L))chg107=0.4d0*y(n_ClONO2,L)
           if(chg107.ge.0.3d0*y(n_HCl,L))chg107=0.3d0*y(n_HCl,L)
           chg108=rr(108,L)*y(n_HOCl,L)*y(n_HCl,L)*dt2
           if(chg108.ge.0.8d0*y(n_HOCl,L))chg108=0.8d0*y(n_HOCl,L)
           if(chg108.ge.0.3d0*y(n_HCl,L))chg108=0.3d0*y(n_HCl,L)
           chg109=rr(109,L)*y(n_N2O5,L)*y(n_HCl,L)*dt2
           if(chg109.ge.0.5d0*y(n_N2O5,L))chg109=0.5d0*y(n_N2O5,L)
           if(chg109.ge.0.3d0*y(n_HCl,L))chg109=0.3d0*y(n_HCl,L)
           changeClONO2=changeClONO2-chg106-chg107
           changeHOCl=chg106-chg108
           changeN2O5=changeN2O5-chg109
           changeHCl=changeHCl-chg107-chg108-chg109
           changeHNO3=changeHNO3+chg106+chg107+chg109
           changeClOx=changeClOx+chg107+chg108+chg109
c Note that really the last 3 produce Cl2, not ClOx, and Cl2 at
c night is stable and doesn't go back into ClONO2, so should
c eventually keep track of Cl2/ClOx partitioning
         endif
         endif
         endif
C
C Make sure nighttime chemistry changes are not too big:
C
         error=.false.
         if(changeNOx.lt.-1.d15.OR.changeNOx.gt.1.d15) then
          write(*,*) 'Big chg@ Itime,I,J,L,NOx ',Itime,I,J,L,changeNOx
          error=.true.
         end if
         if(changeHNO3.lt.-1.d15.OR.changeHNO3.gt.1.d15) then
          write(*,*) 'Big chg@ Itime,I,J,L,HNO3',Itime,I,J,L,changeHNO3
          error=.true.
         end if
         if(changeN2O5.lt.-1.d15.OR.changeN2O5.gt.1.d15) then
          write(*,*) 'Big chg@ Itime,I,J,L,N2O5',Itime,I,J,L,changeN2O5
          error=.true.
         end if
         if(wprodHCHO.lt.-1.d15.OR.wprodHCHO.gt.1.d15) then
          write(*,*)'Big chg@ Itime,I,J,L,HCHO',Itime,I,J,L,wprodHCHO
          error=.true.
         endif
         if(error)call stop_model('nighttime chem: big changes',255)
c
C -- HNO3 --  (HNO3 from gas and het phase rxns )
         changeL(L,n_HNO3)=changeHNO3*pfactor*bymass2vol(n_HNO3)
         IF((trm(i,j,l,n_HNO3)+changeL(l,n_HNO3)).lt.1.d0) THEN
           changeL(l,n_HNO3) = 1.d0 - trm(i,j,l,n_HNO3)
           changeHNO3=changeL(L,n_HNO3)*mass2vol(n_HNO3)*bypfactor
         END IF
C -- N2O5 --  (N2O5 from gas and het phase rxns)
         changeL(L,n_N2O5)=changeN2O5*pfactor*bymass2vol(n_N2O5)
         IF((trm(i,j,l,n_N2O5)+changeL(l,n_N2O5)).lt.1.d0) THEN
           changeL(l,n_N2O5) = 1.d0 - trm(i,j,l,n_N2O5)
           changeN2O5=changeL(L,n_N2O5)*mass2vol(n_N2O5)*bypfactor
         END IF
c -- NOx --   (NOx from gas phase rxns)
         changeL(L,n_NOx)=changeNOx*pfactor*bymass2vol(n_NOx)
         IF((trm(i,j,l,n_NOx)+changeL(l,n_NOx)).lt.1.d0) THEN
           changeL(l,n_NOx) = 1.d0 - trm(i,j,l,n_NOx)
           changeNOx=changeL(L,n_NOx)*mass2vol(n_NOx)*bypfactor
         END IF
c -- ClONO2 --   (ClONO2 from gas and het phase rxns)
         changeL(L,n_ClONO2)=changeClONO2*pfactor*
     *    bymass2vol(n_ClONO2)
         IF((trm(i,j,l,n_ClONO2)+changeL(l,n_ClONO2)).lt.1.d0) THEN
           changeL(l,n_ClONO2) = 1.d0 - trm(i,j,l,n_ClONO2)
           changeClONO2=changeL(L,n_ClONO2)*mass2vol(n_ClONO2)*
     *      bypfactor
         END IF
c -- ClOx --   (ClOx from gas and het phase rxns)
         changeL(L,n_ClOx)=changeClOx*pfactor*bymass2vol(n_ClOx)
         IF((trm(i,j,l,n_ClOx)+changeL(l,n_ClOx)).lt.1.d0) THEN
           changeL(l,n_ClOx) = 1.d0 - trm(i,j,l,n_ClOx)
           changeClOx=changeL(L,n_ClOx)*mass2vol(n_ClOx)*bypfactor
         END IF
         if(pscX.gt.0.9d0)then
c -- HOCl --   (HOCl from het phase rxns)
         changeL(L,n_HOCl)=changeHOCl*pfactor*bymass2vol(n_HOCl)
         IF((trm(i,j,l,n_HOCl)+changeL(l,n_HOCl)).lt.1.d0) THEN
           changeL(l,n_HOCl) = 1.d0 - trm(i,j,l,n_HOCl)
           changeHOCl=changeL(L,n_HOCl)*mass2vol(n_HOCl)*bypfactor
         END IF
c -- HCl --   (HCl from het phase rxns)
         changeL(L,n_HCl)=changeHCl*pfactor*bymass2vol(n_HCl)
         IF((trm(i,j,l,n_HCl)+changeL(l,n_HCl)).lt.1.d0) THEN
           changeL(l,n_HCl) = 1.d0 - trm(i,j,l,n_HCl)
           changeHCl=changeL(L,n_HCl)*mass2vol(n_HCl)*bypfactor
         END IF
         endif  ! pscX > 0.9
C
c Some More Chemistry Diagnostics:
         if(prnchg.and.J.eq.jprn.and.I.eq.iprn.and.L.eq.lprn)then
          write(6,*) 'dark, SALBFJ,sza,I,J,L,Itime= ',SALBFJ(I,J),sza,I
     &,J,L,Itime
          if(pscX.gt.0.9d0)then
           write(6,*) 'There are PSCs, T =',ta(L)
          else
           write(6,*) 'There are no PSCs, T =',ta(L)
          endif
          write(6,198) ay(n_NOx),': ',
     *         changeNOx,' molecules produced; ',
     *      100.d0*(changeNOx)/y(n_NOx,L),' percent of'
     *     ,y(n_NOx,L),'(',1.d9*y(n_NOx,L)/y(nM,L),' ppbv)'
          write(6,198) ay(n_HNO3),': ',
     *         changeHNO3,' molecules produced; ',
     *     100.d0*(changeHNO3)/y(n_HNO3,L),' percent of'
     *     ,y(n_HNO3,L),'(',1.d9*y(n_HNO3,L)/y(nM,L),' ppbv)'
          write(6,198) ay(n_N2O5),': ',
     *         changeN2O5,' net molec produced; ',
     *      100.d0*(changeN2O5)/y(n_N2O5,L),' percent of'
     *     ,y(n_N2O5,L),'(',1.d9*y(n_N2O5,L)/y(nM,L),' ppbv)'
          write(6,198) ay(n_N2O5),': ',
     *            gwprodN2O5,' molec prod fm gas;  ',
     *      100.d0*(gwprodN2O5)/y(n_N2O5,L),' percent of'
     *     ,y(n_N2O5,L),'(',1.d9*y(n_N2O5,L)/y(nM,L),' ppbv)'
          write(6,198) ay(n_N2O5),': ',
     *         -wprod_sulf,' molec prod fm sulf; ',
     *     -100.d0*(wprod_sulf)/y(n_N2O5,L),' percent of'
     *    ,y(n_N2O5,L),'(',1.d9*y(n_N2O5,L)/y(nM,L),' ppbv)'
          write(6,198) ay(n_ClONO2),': ',
     *         changeClONO2,' molecules produced; ',
     *    100.d0*(changeClONO2)/y(n_ClONO2,L),' percent of'
     *    ,y(n_ClONO2,L),'(',1.d9*y(n_ClONO2,L)/y(nM,L),' ppbv)'
          write(6,198) ay(n_ClOx),': ',
     *         changeClOx,' molecules produced; ',
     *    100.d0*(changeClOx)/y(n_ClOx,L),' percent of'
     *    ,y(n_ClOx,L),'(',1.d9*y(n_ClOx,L)/y(nM,L),' ppbv)'
          write(6,198) ay(n_HOCl),': ',
     *         changeHOCl,' molecules produced; ',
     *    100.d0*(changeHOCl)/y(n_HOCl,L),' percent of'
     *    ,y(n_HOCl,L),'(',1.d9*y(n_HOCl,L)/y(nM,L),' ppbv)'
           write(6,198) ay(n_HCl),': ',
     *         changeHCl,' molecules produced; ',
     *    100.d0*(changeHCl)/y(n_HCl,L),' percent of'
     *    ,y(n_HCl,L),'(',1.d9*y(n_HCl,L)/y(nM,L),' ppbv)'
          write(6,199) 'NO2, NO3  = ',y(nNO2,L),yNO3(I,J,L)
         endif
       enddo ! stratosphere L loop

#endif
c
      endif                         ! >>> END sunlight/darkness IF <<<
C

#ifdef SHINDELL_STRAT_CHEM
      LL=LM
#else
      LL=LTROPO(I,J)
#endif
      DO L=1,LL  ! loop over troposphere again (or trop+strat)
C       >> Lower limit on HO2NO2 of 1.0 <<
        if(trm(i,j,l,n_HO2NO2)+changeL(l,n_HO2NO2).lt.1.d0)
     &  changeL(l,n_HO2NO2) = 1.d0 - trm(i,j,l,n_HO2NO2)
c
#ifdef SHINDELL_STRAT_CHEM
C       Limit ClOx to 0.1 or 0.2 ppbv at highest levels:
        if(pres2(L).le.0.05d0)then
          rmrClOx=(trm(I,J,L,n_ClOx)+changeL(L,n_ClOx))*
     *    mass2vol(n_ClOx)*BYDXYP(J)*BYAM(L,I,J)
          if(rmrClOx.gt.1.d-10)changeL(L,n_ClOx)=
     *    (1.d-10-trm(I,J,L,n_ClOx))
     *    /(mass2vol(n_ClOx)*BYDXYP(J)*BYAM(L,I,J))
        endif
        if(pres2(L).le.0.2d0 .and. pres2(L).gt.0.05d0)then
          rmrClOx=(trm(I,J,L,n_ClOx)+changeL(L,n_ClOx))
     *    *mass2vol(n_ClOx)*BYDXYP(J)*BYAM(L,I,J)
          if(rmrClOx.gt.2.d-10)changeL(L,n_ClOx)=
     *    (2.d-10-trm(I,J,L,n_ClOx))
     *     /(mass2vol(n_ClOx)*BYDXYP(J)*BYAM(L,I,J))
        endif
C        
C       Limit BrOx to 1. pptv above 5 hPa:
        if(pres2(L).le.5.d0.and.J.ne.1.and.J.ne.JM)then
           rmrBrOx=(trm(I,J,L,n_BrOx)+changeL(L,n_BrOx))*
     *     mass2vol(n_BrOx)*BYDXYP(J)*BYAM(L,I,J)
           if(rmrBrOx.gt.1.d-12)changeL(L,n_BrOx)=
     *     (1.d-12-trm(I,J,L,n_BrOx))/
     *     (mass2vol(n_BrOx)*BYDXYP(J)*BYAM(L,I,J))
        endif
C
C       Limit Ox to 4. ppmv above 0.2 hPa:
        if(pres2(L).le.0.2d0)then
          rmrOx=(trm(I,J,L,n_Ox)+changeL(L,n_Ox))*
     *     mass2vol(n_Ox)*BYDXYP(J)*BYAM(L,I,J)
          if(rmrOx.gt.4.d-6)changeL(L,n_Ox)=(4.d-6-trm(I,J,L,n_Ox))
     *    /(mass2vol(n_Ox)*BYDXYP(J)*BYAM(L,I,J))
        endif
c
cc      Troposphere halogen sink (Br) & (Cl)
        IF (L.lt.LTROPO(I,J)) THEN 
Corig     rmv=(1.d0-0.95d0**(LTROPO(I,J)-L)) ! had level dependence
          rmv=max(0.d0,    ! remove between 0 and 
     &        min(0.99d0,  ! 99% only...
     &        2.389d-1 * LOG(pres2(L))  - 1.05d0
     &        ))           ! about 5% at 100mb and 60% at 1000mb
          changeL(L,n_ClOx)=changeL(L,n_ClOx)-
     *               (trm(I,J,L,n_ClOx)*rmv)
          changeL(L,n_HCl)=changeL(L,n_HCl)-
     *               (trm(I,J,L,n_HCl)*rmv)
          changeL(L,n_HOCl)=changeL(L,n_HOCl)-
     *               (trm(I,J,L,n_HOCl)*0.85d0)
          changeL(L,n_ClONO2)=changeL(L,n_ClONO2)-
     *               (trm(I,J,L,n_ClONO2)*rmv)
          changeL(L,n_BrOx)=changeL(L,n_BrOx)-
     *               (trm(I,J,L,n_BrOx)*0.9d0)
          changeL(L,n_HBr)=changeL(L,n_HBr)-
     *               (trm(I,J,L,n_HBr)*0.9d0)
          changeL(L,n_HOBr)=changeL(L,n_HOBr)-
     *               (trm(I,J,L,n_HOBr)*0.9d0)
          changeL(L,n_BrONO2)=changeL(L,n_BrONO2)-
     *               (trm(I,J,L,n_BrONO2)*0.9d0)
        END IF
#endif
        if(checktracer_on) call checktracer(I,J)
C
C Save chemistry changes for updating tracers in apply_tracer_3Dsource.
        DO N=1,NTM_CHEM
          tr3Dsource(i,j,l,nChemistry,n) = changeL(l,n) * bydtsrc
        END DO
c Reset radiation O3 values here, if interactive ozone:
        !
        !
        !
#ifdef SHINDELL_STRAT_CHEM
        DU_O3(J)=0.d0! this used to be updated with radiation O3 value
#endif
c
      END DO       ! end current altitude loop
C
CC    if(checktracer_on) call checktracer(I,J)
c
      END DO ! >>>> MAIN I LOOP ENDS <<<<
c
      END DO ! >>>> MAIN J LOOP ENDS <<<<
!$OMP END PARALLEL DO  

#ifdef SHINDELL_STRAT_CHEM
      do j=1,jm
       DU_O3(J)=1000.d0*DU_O3(J)*BYIM
      enddo
c      
      if(MOD(Itime,24).eq.0)then
       write(*,*) 'Ozone column fm -90 to +90'
       write(*,'(46(f4.0,1x))') (DU_O3(J),J=1,JM)
      endif
c
      if(prnchg)then
       write(*,*) 'Map of O3 production from O2 (Herz & SRB + NO SRB)'
       write(*,*)
     & 'NOTE: lower limit of strat is actually LTROPO(I,J), however!'
       write(*,'(a4,7(i10))') 'Jqq:',(Jqq,Jqq=3,44,6)
       do Lqq=LM,LS1,-1 ! inconvenient to print down to LTROPO(I,J)
        pres(Lqq)=(PSF-PTOP)*SIG(Lqq)+PTOP
        do jqq=1,JM
         tempO2(Jqq)=pres(Lqq)/(TX(1,Jqq,Lqq)*1.38d-19)*0.209476d0
         photO2(Jqq,Lqq)=0.d0
         do Iqq=1,IMAXJ(jqq)
       photO2(Jqq,Lqq)=photO2(Jqq,Lqq)+2*ss(27,Lqq,Iqq,Jqq)*tempO2(Jqq)
         enddo
        enddo
        write(6,'(i2,7(1x,E10.3))')
     *   Lqq,(photO2(Jqq,Lqq)/IMAXJ(Jqq),Jqq=3,44,6) !makes 2 Os
       enddo
      endif
#endif
C
C If not doing stratospheric chemistry, overwrite stratosphere:
#ifndef SHINDELL_STRAT_CHEM
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
cc    Stratospheric Overwrite of tracers                cc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C W A R N I N G :If there is ever stratospheric chemistry (i.e. the
C               'changeL' variable for L>LTROPO(I,J) is non-zero at this
C               point in the code), then the stratospheric changes below
C               should be altered.  Currently, they are functions of
C               tracer mass UNCHANGED by chemistry !
C 
C Calculate an average tropical CH4 value near 569 hPa:
C
      CH4_569=0.d0   
      DO J=JS+1,JN
       FACTj= 1.d6*bydxyp(j)
       CH4_569=CH4_569+FACTJ*SUM(F569M*trm(:,J,L569M,n_CH4)*
     & byam(L569M,:,J)+F569P*trm(:,J,L569P,n_CH4)*byam(L569P,:,J))*BYIM
      END DO
      CH4_569=CH4_569/(JN-JS)
c
C     0.55866= 1/1.79 is 1/(obs. tropsph. CH4):
      r179=1.d0/1.79d0 
C
      if(correct_strat_Ox) then
C       Use Ox correction factor for the proper month:
        imonth= 1
        DO m=2,12
         IF((JDAY.LE.MDOFM(m)).AND.(JDAY.GT.MDOFM(m-1))) THEN
          imonth=m
          GOTO 217
         END IF
        END DO
 217    CONTINUE
        if(MOD(Itime,24).eq.0)
     &  write(6,*)'Using Ox strat correction factors for month ',imonth
      end if
C
      do j=1,jm
       J3=MAX(1,NINT(float(j)*float(JCOlat)*BYFJM))! index for CO
       do i=1,IM
         changeL(:,:)=0.d0 ! initilize the change
         do L=LTROPO(I,J)+1,LM    ! >> BEGIN LOOP OVER STRATOSPHERE <<
c         Update stratospheric ozone to amount set in radiation:
          if(correct_strat_Ox) then
            changeL(L,n_Ox)=O3_rad_save(L,I,J)*DXYP(J)*O3MULT
     &      *corrOx(J,L,imonth) - trm(I,J,L,n_Ox)
          else
            changeL(L,n_Ox)=O3_rad_save(L,I,J)*DXYP(J)*O3MULT
     &      - trm(I,J,L,n_Ox)
          end if
          byam75=F75P*byam(L75P,I,J)+F75M*byam(L75M,I,J)
          FACT1=2.0d-9*DXYP(J)*am(L,I,J)*byam75
C         We think we have too little stratospheric NOx, so, to
C         increase the flux into the troposphere, increasing
C         previous stratospheric value by 70% here: GSF/DTS 9.15.03: 
          changeL(L,n_NOx)=
     &    trm(I,J,L,n_Ox)*2.3d-4*1.7d0 - trm(I,J,L,n_NOx)
          changeL(L,n_N2O5)=  FACT1            - trm(I,J,L,n_N2O5)
          changeL(L,n_HNO3)=trm(I,J,L,n_Ox)*4.2d-3-trm(I,J,L,n_HNO3)
          changeL(L,n_H2O2)=  FACT1            - trm(I,J,L,n_H2O2)
          changeL(L,n_CH3OOH)=FACT1            - trm(I,J,L,n_CH3OOH)
          changeL(L,n_HCHO)=  FACT1            - trm(I,J,L,n_HCHO)
          changeL(L,n_HO2NO2)=FACT1*70.d0      - trm(I,J,L,n_HO2NO2)
C         note above: 70. = 1.4E-7/2.0E-9
          changeL(L,n_CO)=(COlat(J3)*COalt(L)*1.D-9*bymass2vol(n_CO)
     &    *AM(L,I,J)*DXYP(J))                   - trm(I,J,L,n_CO)
          changeL(L,n_PAN)     = FACT1*1.d-4 - trm(I,J,L,n_PAN)
          changeL(L,n_Isoprene)= FACT1*1.d-4 - trm(I,J,L,n_Isoprene)
          changeL(L,n_AlkylNit)= FACT1*1.d-4 - trm(I,J,L,n_AlkylNit)
          changeL(L,n_Alkenes) = FACT1*1.d-4 - trm(I,J,L,n_Alkenes)
          changeL(L,n_Paraffin)= FACT1*1.d-4 - trm(I,J,L,n_Paraffin)
c
c         Overwrite stratospheric ch4 based on HALOE obs for tropics
c         and extratropics and scale by the ratio of near-569hPa
c         mixing ratios to 1.79:
c
          CH4FACT=CH4_569*r179
          IF((J.LE.JS).OR.(J.GT.JN)) THEN                ! extratropics
            DO L2=L,LTROPO(I,J)+1,-1
              IF(CH4altX(L2).ne.0.d0) THEN
                CH4FACT=CH4FACT*CH4altX(L2)
                EXIT
              END IF
            END DO
          ELSE IF((J.GT.JS).AND.(J.LE.JN)) THEN           ! tropics
            DO L2=L,LTROPO(I,J)+1,-1
              IF(CH4altT(L2).ne.0.d0) THEN
                CH4FACT=CH4FACT*CH4altT(L2)
                EXIT
              END IF
            END DO
          END IF
C**** ensure that strat overwrite is only a sink
          changeL(L,n_CH4)=-MAX(0d0,trm(I,J,L,n_CH4)-
     &         (AM(L,I,J)*DXYP(J)*CH4FACT*1.d-6))
C
C Save stratosphic change for updating tracer in apply_tracer_3Dsource
          DO N=1,NTM_CHEM
            tr3Dsource(i,j,l,nStratwrite,n) = changeL(l,n)*bydtsrc
          END DO
C
        end do              ! >> END LOOP OVER STRATOSPHERE <<
       end do !i
      end do  !j
#endif
c
      RETURN
      END SUBROUTINE masterchem
c
      SUBROUTINE Crates(I,J)
!@sum Crates calculate chemical reaction rates for each altitude,
!@+   using JPL 00.  Includes special calculations for pressure
!@+   dependent reactions. Specifically:
!@+   #13 CO+OH->HO2+CO2, #15 HO2+HO2->H2O2+O2, #16 OH+HNO3->H2O+NO3,
!@+   and reactions #29, and #42.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on masterchem000_M23p)
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only: LM,JM
#ifdef SHINDELL_STRAT_CHEM
     &     ,ptop,psf,sig
#endif
      USE CONSTANT, only: PI
      USE DYNAMICS, only: LTROPO
      USE TRCHEM_Shindell_COM, only: nr2,nr3,nmm,nhet,ta,ea,rr,pe,
     &                          r1,sb,nst,y,nM,nH2O,ro,sn
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
C
!@param JN J around 30 N
!@param JS J around 30 S
!@param JNN,JSS Js for "high-lat" definition
!@param nlast either ntm_chem or ntm_chem-NregOx for chemistry loops
      INTEGER, PARAMETER :: JS = JM/3 + 1, JN = 2*JM/3
      INTEGER, PARAMETER :: JNN = 5*JM/6, JSS= JM/6 + 1
!@var I,J passed horizontal position indicies
!@var dd,pp,fw,rkp,rk2,rk3M,nb,rrrr,temp dummy "working" variables
!@var L,jj dummy loop variables
!@var byta reciprocal of the local temperature
!@var rkext aerosol extinction from SAGE obs
!@var pscEx NAT PSC surface conc per unit volume (cm^2/cm^3)
!@var Ltop is number of levels with chemistry
      REAL*8 byta, dd, pp, fw, rkp, rk2, rk3M, rrrr,temp
      INTEGER L,jj,nb,Ltop
      INTEGER, INTENT(IN) :: I,J
#ifdef SHINDELL_STRAT_CHEM
!@var PRES local nominal pressure
      REAL*8, DIMENSION(LM) :: PSCEX,PRES,rkext
#endif
C
#ifdef SHINDELL_STRAT_CHEM
      Ltop=LM
      PRES(:)=SIG(:)*(PSF-PTOP)+PTOP
      rkext(:)=0.d0 ! initialize over L
#else
      Ltop=LTROPO(I,J)
#endif
      do L=1,Ltop            !  >>> BEGIN ALTITUDE LOOP <<<
        byta=1.d0/ta(L)
        do jj=1,nr2             ! bimolecular rates start
          IF(ea(jj).ne.0.d0) THEN
            rr(jj,L)=pe(jj)*exp(-ea(jj)*byta)
          ELSE
            rr(jj,L)=pe(jj)
          END IF
c         for #13, k=pe*(1+0.6*(Patm/1013)) Patm=[M]*(T*1.38E-19)
          if(jj.eq.13) rr(jj,L) =
     &    pe(jj)*(1.d0+0.6d0*((y(nM,L)*ta(L)*1.38d-19)/1013.d0))
c         for reaction #15, k=(kc+kp)fw, kc=rr
          if(jj.eq.15)then
            rkp=1.7d-33*y(nM,L)*exp(1000.d0*byta)
            fw=(1.d0+1.4d-21*y(nH2O,L)*exp(2200.d0*byta))
            rr(jj,L)=(rr(jj,L)+rkp)*fw
          endif
c         for #16, k=[pe*exp(-e(jj)/ta(l))]+k3[M]/(1+k3[M]/k2)
          if(jj.eq.16)then
            rk3M=y(nM,l)*1.90d-33*exp(725.d0*byta)
            rk2=4.10d-16*exp(1440.d0*byta)
            rr(jj,L)=rr(jj,L)+rk3M/(1.d0+(rk3M/rk2))
          endif
          if(jj.eq.29)rr(jj,L)=rr(jj,L)/y(nM,L)!PAN+M really PAN
          if(jj.eq.42)rr(jj,L)=rr(jj,L)/y(nM,L)!ROR+M really ROR
        end do                ! bimolecular rates end
c
        do jj=1,nr3           ! trimolecular rates start
         rr(nr2+jj,L)=y(nM,L)*ro(jj)*(300.d0*byta)**sn(jj)
c        if(r1(jj).ge.1.d-30)then
         if(sb(jj).ge.0.01d0)then !dts 3/29/02 alteration of above line
           dd=rr(nr2+jj,L)/(r1(jj)*(300.d0*byta)**sb(jj))
           pp=0.6d0**(1.d0/(1.d0+(log10(dd))**2.))
           rr(nr2+jj,L)=(rr(nr2+jj,L)/(1.d0+dd))*pp
         endif
        end do                ! trimolecular rates end
c
        nb=nr2-nmm
        if(nmm.ge.1) then
          do jj=1,nmm         ! monomolecular rates start
           rrrr=exp(0.5d0*ea(jj+nb)*byta)!0.5 for precision,cor next lin
           rr(jj+nb,L)=rr(nst(jj),L)/
     &     (rrrr*pe(jj+nb)*rrrr*y(nM,l))
          end do              ! monomolecular rates end
        end if
c
#ifdef SHINDELL_STRAT_CHEM
c     Calculate rates for heterogeneous reactions (Divided by solid in Chem1)
c       1=N2O5 + H2O --> 2HNO3          gamma=0.1 (aero), 0.0003 (PSC)
c       2=ClONO2 + H2O --> HOCl + HNO3  gamma=0.001 (aero), 0.001 (PSC)
c       3=ClONO2 + HCl --> Cl2 + HNO3   gamma=0.1
c       4=HOCl + HCl --> Cl2 + H2O      gamma=0.1
c       5=N2O5 + HCl --> ClNO2 + HNO3   gamma=0.003
C       if(l.lt.12.or.l.gt.18)then   !aerosols (14-33 km) & PSCs 14-22 km
        if(pres(l).ge.150.d0 .or. pres(l).le.5.d0)then 
          do jj=nr2+nr3+1,nr2+nr3+nhet
            rr(jj,L)=1.0d-35
          enddo 
          goto5
        else  
cc        Aerosol profiles and latitudinal distribution of extnct coeffs
cc        (in km**-1) from SAGE II data on GISS web site
          if(pres(l).le.150.d0.and.pres(l).gt.90.d0)then
            if(J.GT.JS.AND.J.LE.JN)then !tropics
              rkext(l)=300.d-5
            else
              rkext(l)=10.d-5
            endif
          endif
          if(pres(l).le.90.d0.and.pres(l).gt.56.2d0)then
            if(J.GT.JS.AND.J.LE.JN)then !tropics
              rkext(l)=60.d-5
            elseif(j.le.JSS.or.j.ge.JNN)then !high-lats
              rkext(l)=20.d-5
            else
              rkext(l)=40.d-5
            endif
          endif
          if(pres(l).le.56.2d0.and.pres(l).ge.31.6d0)then
            if(J.GT.JS.AND.J.LE.JN)then !tropics
              rkext(l)=40.d-5
            elseif(j.le.JSS.or.j.ge.JNN)then !high-lats
              rkext(l)=3.d-5
            else
              rkext(l)=13.d-5
            endif
          endif
          if(pres(l).le.31.6d0.and.pres(l).ge.17.8d0)then
            if(J.GT.JS.AND.J.LE.JN)then !tropics
              rkext(l)=20.d-5
            elseif(j.le.JSS.or.j.ge.JNN)then !high-lats
              rkext(l)=1.3d-5
            else
              rkext(l)=6.d-5
            endif
          endif
          if(pres(l).le.17.8d0.and.pres(l).ge.10.0d0)then
            if(J.GT.JS.AND.J.LE.JN)then !tropics
              rkext(l)=1.8d-5
            elseif(j.le.JSS.or.j.ge.JNN)then !high-lats
              rkext(l)=1.2d-6
            else
              rkext(l)=4.8d-6
            endif
          endif
          if(pres(l).le.10.0d0.and.pres(l).ge.4.6d0)then
            if(J.GT.JS.AND.J.LE.JN)then !tropics
              rkext(l)=1.2d-5
            elseif(j.le.JSS.or.j.ge.JNN)then !high-lats
              rkext(l)=0.8d-6
            else
              rkext(l)=3.2d-6
            endif
          endif

cc        PSCs (pressure < 150mb and pressure > 5mb here)
          pscEx(l)=0.d0!NAT PSC surface conc per unit volume (cm^2/cm^3)
          if(pres(l).ge.31.6d0.and.ta(l).le.198.d0)then
            IF((J.LE.JS).OR.(J.GT.JN))pscEx(l)=1.0d-8 !not in tropics
          endif

c         Reaction 1 on sulfate and PSCs      
          temp=sqrt(8.d0*1.38d-16*ta(l)*6.02d23/(PI*108.d0))
          rr(nr2+nr3+1,L)=0.5d0*rkext(l)*1.d-5*temp*0.5d0
          if(pres(l).le.86.2d0 .and. pres(l).gt.4.6d0)
     &    rr(nr2+nr3+1,L)=rr(nr2+nr3+1,L)*.45d0
          if(pres(l).gt.31.6d0) rr(nr2+nr3+1,L)=
     &    rr(nr2+nr3+1,L)+0.25d0*pscEx(l)*temp*0.0009d0

c         Reaction 2 on sulfate and PSCs      
          temp=sqrt(8.d0*1.38d-16*ta(l)*6.02d23/(PI*97.d0))
          rr(nr2+nr3+2,L)=0.5d0*rkext(l)*1.d-5*rr(nr2+nr3+2,L)*1.d-3
          if(pres(l).gt.31.6d0) rr(nr2+nr3+2,L)=
     &    rr(nr2+nr3+2,L)+0.25d0*pscEx(l)*temp*0.001d0

          if(pres(l).gt.31.6d0) then
c           rr(nr2+nr3+3,L)=sqrt(8.d0*1.38d-16*ta(l)*6.02d23/(PI*97.d0))
            rr(nr2+nr3+3,L)=0.25d0*pscEx(l)*temp*0.1d0
            rr(nr2+nr3+4,L)=
     &      sqrt(8.d0*1.38d-16*ta(l)*6.02d23/(PI*52.d0))
            rr(nr2+nr3+4,L)=0.25d0*pscEx(l)*rr(nr2+nr3+4,L)*0.1d0
            rr(nr2+nr3+5,L)=
     &      sqrt(8.d0*1.38d-16*ta(l)*6.02d23/(PI*108.d0))
            rr(nr2+nr3+5,L)=0.25d0*pscEx(l)*rr(nr2+nr3+5,L)*0.003d0
          endif
        endif  
   5    continue       
#endif
      end do                  !  >>> END ALTITUDE LOOP <<<
 
      RETURN
      END SUBROUTINE Crates
c
c
#ifdef TRACERS_SPECIAL_Shindell
      SUBROUTINE checktracer(I,J)
!@sum checktracer for various debugging of tracer chemistry
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on masterchem000_M23p)
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only  : Itime, LM
      USE DYNAMICS, only   : LTROPO
      USE TRACER_COM, only : ntm, trname, n_Ox, ntm_chem
#ifdef regional_Ox_tracers
     & ,NregOx
#endif
      USE TRCHEM_Shindell_COM, only: y, nM
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var I,J passed horizontal position indicies
!@var L,igas dummy loop variables
!@var tlimit if tracer goes above this limit, model stops
!@var checkOx logical: should I check for large tropospheric Ox?
!@var checkmax logical: should I check for large tracers throughout?
!@var checkNeg logical: should I check for negative tracers?
!@var checkNaN logical: should I check for unreal tracers?
!@var nlast either ntm or ntm-NregOx
C
      INTEGER, PARAMETER :: nlast=
#ifdef regional_Ox_tracers
     & ntm_chem-NregOx
#else
     & ntm_chem
#endif
      INTEGER L, igas
      INTEGER, INTENT(IN) :: I,J
      REAL*8, DIMENSION(nlast) :: tlimit
      DATA tlimit/
     &9.d-5,1.d-5,1.d-7,3.d-6,1.d-1,1.d-6,3.d-6,1.d-1,1.d-1,1.d-1,
     &1.d-1, 1.d-1, 1.d-1, 1.d-1, 1.d-1
#ifdef SHINDELL_STRAT_CHEM
     &,1.d-1, 1.d-1, 1.d-1, 1.d-1, 1.d-1
     &,1.d-1, 1.d-1, 1.d-1, 1.d-1, 1.d-1/
#else
     & /
#endif
    
      LOGICAL checkOx, checkmax, checkNeg, checkNan
      DATA checkNeg /.true./
      DATA checkNan /.true./
      DATA checkOx  /.true./
      DATA checkmax /.false./
c
      IF(i.eq.1.and.j.eq.1)
     & WRITE(6,*) 'WARNING: checktracer call is active.'
      IF(checkmax)
     & call stop_model('checktracer: set tlimit for tracers 11->25',255)
C     please (re)set tlimit values for tracers 11 through 15 in the
C     data statement above. Then, please delete the above stop.
C
C check if ozone gets really big in the troposphere:
       IF(checkOx) THEN
       do L=1,LTROPO(I,J)
         if(y(n_Ox,L)/y(nM,L).gt.1.d-5) then
           write(6,*)'Ox @ I,J,L,Ox,tau:',I,J,L,y(n_Ox,L),Itime
           call stop_model('checktracer: Ox too big in tropo.',255)
         end if
       end do
       END IF
c general check on maximum of tracers:
      IF(checkmax) THEN
      do L=1,LM
       do igas=1,nlast
        if(y(igas,L)/y(nM,L).gt.tlimit(igas)) then
          write(6,*) trname(igas),'@ I,J,L,Ox :',I,J,L,y(igas,L)
          call stop_model('checktracer: tracer upper limit',255)
        end if
       end do
      end do
      END IF
c check for negative tracers:
      IF(checkNeg) THEN
      do L=1,LM
      do igas=1,nlast
       if(y(igas,L).lt.0.d0) THEN
         write(6,*)trname(igas),
     &   'negative @ tau,I,J,L,y:',Itime,I,J,L,y(igas,L)
         call stop_model('checktracer: tracer is negative',255)
       end if
      enddo
      end do
      END IF
c check for unreal (not-a-number) tracers:
      IF(checkNaN) THEN
      do L=1,LM
      do igas=1,nlast
        if(.NOT.(y(igas,L).gt.0.d0.OR.y(igas,L).le.0.d0)) THEN
         write(6,*)trname(igas),
     &   'is not a number @ tau,I,J,L,y:',Itime,I,J,L,y(igas,L)
         call stop_model('checktracer: tracer is NaN',255)
        end if
      enddo
      end do
      END IF
      RETURN

      END SUBROUTINE checktracer
#endif
      
      

