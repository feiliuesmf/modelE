#include "rundeck_opts.h"
      SUBROUTINE masterchem
!@sum masterchem main chemistry routine
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on masterchem000_M23p)   
!@calls photoj,checktracer,Crates,Oxinit,HOxfam,NOxfam,chemstep
C
C IF ALTERING THIS ROUTINE, PLEASE SEE THE WARNING ABOUT THE CHANGEL
C VARIABLE IN THE STRATOSPHERIC OVERWRITE SECTION.
c
C**** GLOBAL parameters and variables:
c
      USE SOMTQ_COM, only   : qmom
      USE DOMAIN_DECOMP,only: GRID,GET,AM_I_ROOT,PACK_COLUMN,
     &                        GLOBALSUM, PACK_DATA, PACK_DATAj,
     &                        write_parallel
      USE MODEL_COM, only   : Q,JDAY,IM,JM,sig,ptop,psf,ls1,JYEAR,
     &                        COUPLED_CHEM
      USE CONSTANT, only    : radian,gasc,mair,mb2kg,pi,avog
      USE DYNAMICS, only    : pedn,LTROPO
      USE FILEMANAGER, only : openunit,closeunit
      USE RAD_COM, only     : COSZ1,alb,rcloudfj=>rcld,
     &                        rad_to_chem,O3_tracer_save,H2ObyCH4,
     &                        SRDN,rad_to_file
      USE GEOM, only        : BYDXYP, DXYP, LAT_DG, IMAXJ
      USE FLUXES, only      : tr3Dsource
      USE TRACER_COM, only  : n_Ox,n_NOx,n_N2O5,n_HNO3,n_H2O2,n_CH3OOH,
     &                        n_HCHO,n_HO2NO2,n_CO,n_CH4,n_PAN,
     &                        n_Isoprene,n_AlkylNit,n_Alkenes,
     &                        n_Paraffin,ntm_chem,n_DMS,n_MSA,n_SO2,
     &                        n_SO4,n_H2O2_s,oh_live,no3_live,
     &                        nChemistry,nStratwrite,rsulf1,rsulf2,
     &                        rsulf3,rsulf4,TR_MM,trname
#ifdef SHINDELL_STRAT_CHEM
     &                        ,n_HBr,n_HOCl,n_HCl,n_ClONO2,n_ClOx,
     &                        n_BrOx,n_BrONO2,n_CFC,n_N2O,n_HOBR
#endif
#ifdef regional_Ox_tracers
     &                        ,NregOx,n_OxREG1
#endif
#ifdef TRACERS_HETCHEM
     &                        ,krate,n_N_d1,n_N_d2,n_N_d3
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      USE TRACER_SOURCES, only: avg_model,n__sw
#endif
      USE TRDIAG_COM, only    : jls_N2O5sulf,tajls=>tajls_loc,
     & taijs=>taijs_loc,ijs_JH2O2,ijs_NO3,jls_COp,jls_COd,jls_Oxp,
     & jls_Oxd
#ifdef SHINDELL_STRAT_CHEM
     &                           ,jls_ClOcon,jls_H2Ocon
#endif
      USE TRCHEM_Shindell_COM

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@param by35 1/35 used for spherical geometry constant
!@param JN J around 30 N
!@param JS J around 30 S
!@param JNN,JSS Js for "high-lat" definition
!@param nlast either ntm_chem or ntm_chem-NregOx for chemistry loops
      INTEGER, PARAMETER :: JS = JM/3 + 1, JN = 2*JM/3
      INTEGER, PARAMETER :: JNN = 5*JM/6, JSS= JM/6 + 1
      INTEGER, PARAMETER :: nlast =
#ifdef regional_Ox_tracers
     &                              ntm_chem-NregOx
#else
     &                              ntm_chem
#endif
      REAL*8, PARAMETER  :: by35=1.d0/35.d0
      REAL*8, PARAMETER  :: bymair = 1.d0/mair
      REAL*8, PARAMETER, DIMENSION(LM) :: thick = (/
     & 0.3d0,0.36d0,0.44d0,0.65d0,1.2d0,1.6d0,2.0d0,2.0d0,1.6d0,
     & 1.6d0,1.5d0,1.6d0,1.9d0,2.5d0,3.7d0,3.6d0,4.0d0,5.2d0,
     & 7.4d0,8.6d0,8.6d0,9.2d0,12.4d0/)
!@var FASTJ_PFACT temp factor for vertical pressure-weighting
!@var FACT1,2,3 temp variable for start overwrite
!@var bydtsrc reciprocal of the timestep dtsrc
!@var local logical for error checking 
!@var byam75 the reciprocal air mass near 75 hPa level
!@var average tropospheric ch4 value near 569 hPa level
!@var PRES2 local nominal pressure for verticle interpolations
!@var thick thickness of each layer in km
!@var photO2_glob for O2->Ox diagnostic
!@var ClOx_old total ClOx at start of chemical timestep
!@var ClTOT total chlorine in all forms (reactive and reservoir)
!@var colmO2, colmO3 are overhead oxygen and ozone columns
!@var CH4FACT, r179 for setting CH4 ICs and strat distribution
!@var changeClONO2,changeClOx,changeHOCl,changeHCl nighttime changes
!@var changehetClONO2 nighttime het change in ClONO2 (on sulfate)
!@var pscX flag for psc existence
!@var chg106,chg107,chg108,chg109 reaction rates for het rxns on pscs
!@var rmrClOx,rmrBrOx dummy vars with mixing ratios of halogens
!@var rmv dummy variable for halogne removal in trop vs height
!@var changeL 2D array holds the local change due to chem or strat 
!@+   overwrite until adding to tr3Dsource (replaces 4D "change")
!@var PIfact strat-overwrite adjustment for preindustrial runs
!@var pfactor to convert units on species chemical changes
!@var bypfactor to convert units on species chemical changes
!@var dNO3,gwprodHNO3,gprodHNO3,gwprodN2O5,changeAldehyde,
!@+   changeAlkenes,changeIsoprene,changeHCHO,changeAlkylNit,
!@+   changeHNO3,changeNOx,changeN2O5,wprodHCHO working variables to 
!@+   calculate nighttime chemistry changes
!@var rlossN,rprodN,ratioN variables for nitrogen conservation
!@var I,J,L,N,igas,inss,LL,Lqq,JJ,J3,L2,n2 dummy loop variables
!@var avgTT_CH4 Itime avg CH4 # density at LTROPO between 20N and 20S
!@var avgTT_H2O Itime avg H2O # density at LTROPO between 20N and 20S
!@var countTT # of points between 20N and 20S on LTROPO plane
!@var aero yes(1) or no(0) tag of non-zero rkext from Crates
!@var imonth,m only needed to choose Ox strat correction factors
!@var maxl chosen tropopause 0=LTROPO(I,J), 1=LS1-1
!@var sumOx for summing regional Ox tracers
!@var bysumOx reciprocal of sum of regional Ox tracers
      REAL*8, DIMENSION(LM,NTM) :: changeL
      REAL*8, DIMENSION(NTM)    :: PIfact
      REAL*8, DIMENSION(LM)     :: PRES2
      REAL*8 :: FACT1,FACT2,FACT3,FACT4,FACT5,FACT6,FACT7,fact_so4,
     &  FASTJ_PFACT,bydtsrc,byam75,byavog,CH4FACT,r179,rlossN,
     &  rprodN,ratioN,pfactor,bypfactor,gwprodHNO3,gprodHNO3,
     &  gwprodN2O5,wprod_sulf,wprodCO,dNO3,wprodHCHO,prod_sulf,
     &  RVELN2O5,changeAldehyde,changeAlkenes,changeAlkylNit,
     &  changeIsoprene,changeHCHO,changeHNO3,changeNOx,changeN2O5,
     &  changeOx,fraQ,CH4_569,count_569
#ifdef TRACERS_HETCHEM
      REAL*8 :: changeN_d1,changeN_d2,changeN_d3
#endif
#ifdef regional_Ox_tracers
      REAL*8 :: sumOx, bysumOx
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      REAL*8 :: temp_SW
#endif
#ifdef SHINDELL_STRAT_CHEM
      REAL*8, DIMENSION(LM)     :: ClOx_old  
      REAL*8 :: CLTOT,colmO2,colmO3,changeClONO2,changeClOx,
     & changeHOCl,changeHCl,changehetClONO2,pscX,chg106,chg107,
     & chg108,chg109,rmrClOx,rmrBrOx,rmv,rmrOx,avgTT_H2O,avgTT_CH4,
     & countTT
      INTEGER, DIMENSION(LM)    :: aero
#endif
      INTEGER                   :: igas,LL,I,J,L,N,inss,Lqq,J3,L2,n2,
     &                             Jqq,Iqq,imonth,m,maxl,iu
      LOGICAL                   :: error, jay
      CHARACTER*4               :: ghg_name
      CHARACTER*80              :: ghg_file
      character(len=300)        :: out_line

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Some variables defined especially for MPI compliance :         C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8, dimension(LM,IM,JM,5) :: rad_to_file_glob
      real*8, dimension(LM,IM,JM)   :: ss27_glob
      real*8 :: ss27(LM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO)
      real*8, dimension(JM)         :: DU_O3_glob
#ifdef SHINDELL_STRAT_CHEM
      real*8, dimension(JM,LM)      :: photO2_glob
#endif
      real*8 :: avgTT_CH4_part(GRID%J_STRT_HALO:GRID%J_STOP_HALO),
     &          avgTT_H2O_part(GRID%J_STRT_HALO:GRID%J_STOP_HALO),
     &            CH4_569_part(GRID%J_STRT_HALO:GRID%J_STOP_HALO),
     &            countTT_part(GRID%J_STRT_HALO:GRID%J_STOP_HALO),
     &          count_569_part(GRID%J_STRT_HALO:GRID%J_STOP_HALO)
      
      INTEGER :: J_0, J_1, J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE     
      
      CALL GET(grid, J_STRT    =J_0,  J_STOP    =J_1,
     &               J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      
      byavog = 1.d0/avog
#ifdef INITIAL_GHG_SETUP
C-------- special section for ghg runs ---------
      write(out_line,*)'Warning: INITIAL_GHG_SETUP is on!'
      call write_parallel(trim(out_line))
      if(use_rad_ch4>0 .or. use_rad_n2o>0 .or. use_rad_cfc>0)then
        rad_to_file(:,:,J_0:J_1,1)=rad_to_chem(:,:,J_0:J_1,1)
        rad_to_file(:,:,J_0:J_1,2)=rad_to_chem(:,:,J_0:J_1,2)
        do j=J_0,J_1
          rad_to_file(:,:,j,3)=rad_to_chem(:,:,j,3)*2.69e20*byavog*
     &    dxyp(j)*tr_mm(n_N2O) ! i.e. in trm units now!
          rad_to_file(:,:,j,4)=rad_to_chem(:,:,j,4)*2.69e20*byavog*
     &    dxyp(j)*tr_mm(n_CH4) ! i.e. in trm units now!
          rad_to_file(:,:,j,5)=rad_to_chem(:,:,j,5)*2.69e20*byavog*
     &    dxyp(j)*tr_mm(n_CFC)*fact_CFC ! i.e. in trm units now!
        enddo 
        do m=1,5
          call PACK_COLUMN
     &    (grid, rad_to_file(:,:,:,m), rad_to_file_glob(:,:,:,m))
        end do
        if(AM_I_ROOT( ))then
          write(ghg_name,'(I4)') JYEAR  
          ghg_file='GHG_IC_'//ghg_name
          call openunit(ghg_file,iu,.true.,.false.)
          do j=1,5
            write(iu)ghg_file,rad_to_file_glob(:,:,:,j)
          enddo
          call closeunit(iu)          
          write(6,*)'Stopping in masterchem to output inital'
          write(6,*)'conditions for ghgs. You must now recompile'
          write(6,*)'(rerun setup) after removing the'
          write(6,*)'INITIAL_GHG_SETUP directive from your rundeck'
          write(6,*)' & linking the GHGic file to: ',trim(ghg_file)
          write(6,*)'Address questions to Greg Faluvegi. Thanks.'
          call stop_model('Normal INITIAL_GHG_SETUP stop.',255)
        endif
      endif 
#endif

C Some INITIALIZATIONS :
      byavog  = 1.d0/avog
      bydtsrc = 1.d0/dtsrc
      BYFJM   = 1.d0/real(JM)
      PRES2(:)= SIG(:)*(PSF-PTOP)+PTOP
#ifdef SHINDELL_STRAT_CHEM
      if(H2ObyCH4 /= 0.)               
     &call stop_model('H2ObyCH4 not = 0 and strat chem is on.',255)
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
C Not really related to chemistry, but convenient place to update
C running-averages for interactive wetlands CH4:
      do J=J_0,J_1; do I=1,IMAXJ(J)
        temp_SW=ALB(I,J,1)*(SRDN(I,J)+1.d-20)*COSZ1(I,J)
        call running_average(temp_SW,I,J,1.d0,n__sw)
      end do      ; end do
#endif

C Calculation of gas phase reaction rates for sulfur chemistry:
      CALL GET_SULF_GAS_RATES
      
#ifdef TRACERS_HETCHEM
c Calculation of removal rates on dust surfaces:
      CALL HETCDUST
#endif

c Set "chemical time step". Really this is a method of applying only
c a fraction of the chemistry change to the tracer mass for the first
c 30 hours.  That fraction is: dt2/dtscr.  E.g. in the first hour it
c is (dtsrc/24)/dtsrc = 1/24th of the chemistry change is applied.
c This is to work around initial instabilities.

      if(Itime-ItimeI <= 3)then
        dt2=dtsrc/24.d0          ! e.g. 150.
      elseif(Itime-ItimeI > 3 .and. Itime-ItimeI <= 6)then
        dt2=dtsrc/12.d0          ! e.g. 300.
      elseif(Itime-ItimeI > 6 .and. Itime-ItimeI <= 11)then
        dt2=dtsrc/6.d0           ! e.g. 600.
      elseif(Itime-ItimeI > 11 .and. Itime-ItimeI <= 30)then
        dt2=dtsrc/2.4d0          ! e.g. 1500.
      elseif(Itime-ItimeI > 30)then
        dt2=dtsrc                ! e.g. 3600
      endif

c Calculate new photolysis rates every n_phot main timesteps:
      MODPHOT= 0 !!!!!!!!!!!!MOD(Itime-ItimeI,n_phot)

C CALCULATE TX, THE REAL TEMPERATURE:
C (note this section is already done in DIAG.f)
      IF(HAVE_SOUTH_POLE) THEN
        DO L=1,LM
          TX(1,1,L)=T(1,1,L)*PK(L,1,1)
          TX(2:IM,1,L)=TX(1,1,L)
        END DO
      ENDIF  
      IF(HAVE_NORTH_POLE) THEN
        DO L=1,LM
          TX(1,JM,L)=T(1,JM,L)*PK(L,1,JM)
          TX(2:IM,JM,L)=TX(1,JM,L)
        END DO
      ENDIF
      DO L=1,LM
        DO J=J_0S,J_1S
          TX(1:IM,J,L)=T(1:IM,J,L)*PK(L,1:IM,J)
        END DO
      END DO

#ifdef SHINDELL_STRAT_CHEM
C info to set strat H2O based on tropical tropopause H2O and CH4:
      if(Itime == ItimeI)then
        avgTT_H2O_part(:)=0.d0
        avgTT_CH4_part(:)=0.d0
        countTT_part(:)=0.d0
        do J=J_0S,J_1S
          if(LAT_DG(J,1) >= -20. .and. LAT_DG(J,1) <= 20.)then
            do I=1,IMAXJ(J)
              avgTT_H2O_part(J) = avgTT_H2O_part(J) + 
     &                            Q(I,J,LTROPO(I,J))*MWabyMWw
              if(use_rad_ch4 > 0) then
                avgTT_CH4_part(J) = avgTT_CH4_part(J) + 
     &          rad_to_chem(LTROPO(I,J),I,J,4)
     &          *2.69d20*byavog*mair*BYAM(LTROPO(I,J),I,J)
              else
                avgTT_CH4_part(J) = avgTT_CH4_part(J) +
     &          trm(I,J,LTROPO(I,J),n_CH4)
     &          *mass2vol(n_CH4)*BYDXYP(J)*BYAM(LTROPO(I,J),I,J)
              endif
              countTT_part(J) = countTT_part(J) + 1.d0
            end do
          end if
        end do
        CALL GLOBALSUM(grid, avgTT_CH4_part, avgTT_CH4, all=.true.)
        CALL GLOBALSUM(grid, avgTT_H2O_part, avgTT_H2O, all=.true.)
        CALL GLOBALSUM(grid, countTT_part,   countTT,   all=.true.)
        if(countTT <= 0.)call stop_model('countTT.le.0',255)
      end if
#endif


CCCC!$OMP  PARALLEL DO PRIVATE (changeL, FASTJ_PFACT,
CCCC!$OMP* rlossN,rprodN,ratioN,pfactor,bypfactor,gwprodHNO3,gprodHNO3,
CCCC!$OMP* gwprodN2O5,wprod_sulf,wprodCO,dNO3,wprodHCHO,prod_sulf,rveln2o5,
CCCC!$OMP* changeAldehyde,changeAlkenes,changeIsoprene,changeHCHO,
CCCC!$OMP* changeAlkylNit,changeHNO3,changeNOx,changeN2O5,changeOx,FACT_SO4,
CCCC#ifdef TRACERS_HETCHEM
CCCC!$OMP* changeN_d1,changeN_d2,changeN_d3,
CCCC#endif
CCCC#ifdef SHINDELL_STRAT_CHEM
CCCC!$OMP* aero, CLTOT, ClOx_old, COLMO2, COLMO3, changeClONO2, changeClOx,
CCCC!$OMP* changehetClONO2, changeHOCl, changeHCl,
CCCC!$OMP* chg106, chg107, chg108, chg109, fraQ,
CCCC!$OMP* pscx, rmrclox, rmrbrox, rmrox, rmv,
CCCC#endif
CCCC#ifdef regional_Ox_tracers
CCCC!$OMP* bysumOx, sumOx,     n2,
CCCC#endif
CCCC!$OMP* igas, inss, J, jay, L, LL, Lqq, maxl, N, error )
CCCC!$OMP* SHARED (N_NOX,N_HNO3,N_N2O5,N_HCHO,N_ALKENES,N_ISOPRENE,
CCCC!$OMP* N_ALKYLNIT)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO J=J_0,J_1                  ! >>>> MAIN J LOOP BEGINS <<<<

      DO I=1,IMAXJ(J)               ! >>>> MAIN I LOOP BEGINS <<<<
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if(checktracer_on) call checktracer(I,J)
      select case(which_trop)
      case(0); maxl=ltropo(I,J)
      case(1); maxl=ls1-1
      case default; call stop_model('which_trop problem 4',255)
      end select

      DO L=1,LM
#ifdef regional_Ox_tracers
C      Per Volker: make sure the regional Ox tracers didn't
C      diverge from total Ox:
       sumOx  =0.d0
       bysumOx=0.d0
       do n2=n_OxREG1,ntm_chem
         sumOx=sumOx+trm(i,j,l,n2)
       end do
       sumOx=MAX(sumOx,1.d-9) ! Per Drew: minimum to prevent NaN's
       bysumOx=1.d0/sumOx 
       do n2=n_OxREG1,ntm_chem
         trm(i,j,l,n2)=trm(i,j,l,n2)*trm(i,j,l,n_Ox)*bysumOx
       end do
#endif
c Initialize the 2D change variable:
       changeL(L,:)=0.d0
c Save presure and temperature in local arrays:
       pres(L)=PMID(L,I,J)
       ta(L)  =TX(I,J,L)
c Calculate M and set fixed ratios for O2 & H2:
       y(nM,L)=pres(L)/(ta(L)*1.38d-19)
       y(nO2,L)=y(nM,L)*pfix_O2
#ifdef SHINDELL_STRAT_CHEM
       if(pres2(l) > 20.d0)then
         y(nH2,L)=y(nM,L)*pfix_H2
       else
         ! Was: y(nH2,L)=y(nM,L)*pfix_H2*7.d1/(7.d1+L-maxl+1)
         ! Now: a drop of 0.3 molec/cm3 per 12 hPa decrease:
         y(nH2,L)=y(nM,L)*pfix_H2 + 2.5d-2*(pres2(L)-20.d0)
       endif
       CLTOT=0.d0
#else
       y(nH2,L)=y(nM,L)*pfix_H2
#endif

c Tracers (converted from mass to number density):
       do igas=1,nlast
         y(igas,L)=trm(I,J,L,igas)*y(nM,L)*mass2vol(igas)*
     &   BYDXYP(J)*BYAM(L,I,J)
       enddo

C Concentrations of DMS and SO2 for sulfur chemistry:
       if (coupled_chem == 1) then
         ydms(i,j,l)=trm(i,j,l,n_dms)*y(nM,L)*(28.0D0/62.0D0)*
     &   BYDXYP(J)*BYAM(L,I,J)
         yso2(i,j,l)=trm(i,j,l,n_so2)*y(nM,L)*(28.0D0/64.0D0)*
     &   BYDXYP(J)*BYAM(L,I,J)
       else
         ! Convert from pptv to molecule cm-3:
         ydms(i,j,l)=dms_offline(i,j,l)*1.0D-12*y(nM,L)
         yso2(i,j,l)=so2_offline(i,j,l)*1.0D-12*y(nM,L)
       endif

c If desired, fix the methane concentration used in chemistry:
C WHY IS THIS NECESSARY ANY MORE, NOW THAT GET_CH4_IC IS CALLED?
#ifndef SHINDELL_STRAT_CHEM
       if(fix_CH4_chemistry == 1) THEN
         if(J < JEQ)then ! SH
           y(n_CH4,L)=y(nM,L)*pfix_CH4_S
         else            ! NH
           y(n_CH4,L)=y(nM,L)*pfix_CH4_N
         endif
       end if
#endif
#ifdef SHINDELL_STRAT_CHEM
c Save initial ClOx amount for use in ClOxfam:
       ClOx_old(L)=trm(I,J,L,n_ClOx)*y(nM,L)*mass2vol(n_ClOx)*
     & BYDXYP(J)*BYAM(L,I,J)
#endif

c Limit N2O5 number density:
       if(y(n_N2O5,L) < 1.) y(n_N2O5,L)=1.d0
c Set H2O, based on Q:
       y(nH2O,L)=Q(I,J,L)*MWabyMWw*y(nM,L)
#ifdef SHINDELL_STRAT_CHEM
c Initialize stratospheric y(H2O) & GCM Q variable (!),
c based on tropical tropopause H2O and CH4:
       if(Itime == ItimeI .and. L > LTROPO(I,J)) then
         y(nH2O,L) =  y(nM,L)*(avgTT_H2O/countTT +
     &   2.d0*(avgTT_CH4/countTT-y(n_CH4,L)/y(nM,L)))
         fraQ=(y(nH2O,L)/(y(nM,L)*MWabyMWw))/Q(I,J,L)
         Q(I,J,L)=y(nH2O,L)/(y(nM,L)*MWabyMWw)
         if(fraQ < 1.)qmom(:,i,j,l)=qmom(:,i,j,l)*fraQ
       end if 
#endif

c Initialize various other species:
c - set [NO]=0 (where?) for first HOx calc, NO2 = NOx:
c - set reactive species for use in family chemistry & nighttime NO2:

       if(Itime == ItimeI)then 
         y(nAldehyde,L)=y(nM,L)*pfix_Aldehyde
       else
         y(nAldehyde,L)=yAldehyde(I,J,L)
       endif
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
       y(nCl2,L)     =yCl2(I,J,L)
       y(nCl2O2,L)   =yCl2O2(I,J,L)
       y(nOClO,L)    =y(n_ClOx,L)*pOClOx(I,J,L)
       y(nClO,L)     =y(n_ClOx,L)*pClOx(I,J,L)
       y(nCl,L)      =y(n_ClOx,L)*pClx(I,J,L)
       y(nBr,L)      =y(n_BrOx,L)*(1.d0-pBrOx(I,J,L))
       y(nBrO,L)     =y(n_BrOx,L)*pBrOx(I,J,L)
#endif
      END DO ! L

C For solar zenith angle, we use the arccosine of the COSZ1
C from the radiation code, which is the cosine of the solar zenith 
C angle averaged over the physics time step.
C If the solar zenith angle (sza) from the radiation code is > 90 deg,
C (and hence COSZ1 is set to 0), recalculate it with get_sza routine:
      IF(COSZ1(I,J) == 0.d0) THEN
        call get_sza(I,J,sza)
      ELSE
        sza = acos(COSZ1(I,J))*byradian
      END IF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                 BEGIN PHOTOLYSIS                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if(MODPHOT == 0)then ! i.e. not every time step

       ! additional SUNLIGHT criterion (see also fam chem criterion):
       if((ALB(I,J,1) /= 0.d0).AND.(sza < szamax))then
       
c       define column temperatures to be sent to FASTJ:
        TFASTJ = ta

#ifdef SHINDELL_STRAT_CHEM
c       define pressures to be sent to FASTJ (centers):
        PFASTJ2(1:LM) = PMID(1:LM,I,J)
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
        if(LM /= 23) call stop_model('check PFASTJ2(LM+X).',255)
        
#ifdef SHINDELL_STRAT_CHEM
c Pass O3 array (in ppmv) to fastj. Above these levels fastj2 uses
C Nagatani climatological O3, read in by chem_init: 
        if(LM /= 23)call stop_model('check O3_FASTJ(near LM)',255)
        DO LL=1,LM
          if(LL >= LM-1)y(nO3,LL)=y(n_Ox,LL)
          O3_FASTJ(LL)=y(nO3,LL)/y(nM,LL)
        END DO
#else
c Interpolate O3 (in ppmv) from bottom maxl model sigma levels
c (i.e. those levels normally in troposphere) onto bottom 2*maxl
c FASTJ levels. Above these levels fastj uses Nagatani climatological
c O3, read in by chem_init:
c
c       define O3 to be sent to FASTJ (centers):
        DO LL=2,2*maxl,2
          O3_FASTJ(LL)=y(nO3,LL/2)
        ENDDO
c       define O3 to be sent to FASTJ (edges):
        O3_FASTJ(1)=y(nO3,1)*O3_1_fact ! see parameter declaration...
        DO LL=3,2*maxl-1,2
c         interpolation factor, based on pressure:
          FASTJ_PFACT=(PFASTJ(LL)-PFASTJ(LL-1))/
     &    (PFASTJ(LL+1)-PFASTJ(LL-1))
          O3_FASTJ(LL)=y(nO3,(LL-1)/2)+
     &    (y(nO3,((LL-1)/2)+1)-y(nO3,(LL-1)/2))*FASTJ_PFACT
c         lower limit on O3:
          IF(O3_FASTJ(LL) < 0.d0) O3_FASTJ(LL)=-O3_FASTJ(LL) ! weird
        ENDDO
c       Convert to ppm (units used in fastj)
        DO LL=1,2*maxl
          O3_FASTJ(LL)=O3_FASTJ(LL)/(PFASTJ(LL)*2.55d10)
        ENDDO
#endif

        call photoj(I,J) ! CALL THE PHOTOLYSIS SCHEME

C Define and alter resulting photolysis coefficients (zj --> ss):
#ifdef SHINDELL_STRAT_CHEM
        colmO2=5.6d20 
        colmO3=5.0d16 
#endif
        DO L=JPNL,1,-1
          do inss=1,JPPJ
            ss(inss,L,I,J)=zj(L,inss)
            if(inss /= 22)ss(inss,L,I,J)=ss(inss,L,I,J)*
     &      (by35*SQRT(1.224d3*(cos(ABS(LAT_DG(J,1))*radian))**2.
     &      +1.d0)) ! spherical correction but excluding ClONO2.
#ifdef SHINDELL_STRAT_CHEM
c           reduce rates for gases that photolyze in window region
c           (~200nm):
            if(inss == 28)ss(inss,L,I,J)=ss(inss,L,I,J)*1.25d-1 !N2O
            if(inss == 26)ss(inss,L,I,J)=ss(inss,L,I,J)*1.25d-1 !CFC
#endif
          enddo
          taijs(i,j,ijs_JH2O2(l))=taijs(i,j,ijs_JH2O2(l))+ss(4,l,i,j)
#ifdef SHINDELL_STRAT_CHEM
          colmO2=colmO2+y(nO2,L)*thick(L)*1.d5
          colmO3=colmO3+y(nO3,L)*thick(L)*1.d5
c SF3 is photolysis of water in Schumann-Runge bands based on:
c Nicolet, Pl. Space Sci., p 871, 1983.
C SF3_fact is, if x[ ] = bin4_flux[ ]:
C {(x[present] - x[1988]) / (x[1991] - x[1988])} * 0.1E-6
C This gets ADDED to the 1.3E-6 factor in the SF3 calculation. Here,
C bin4_flux is a proxy for the flux from all 175-200nm bins. 
C (Drew says the ratio would be the same.)
          if(SF2_fact == 0.)call stop_model('SF2_fact=0 in master',255)
          if(pres2(L) <= 10.)then
            if((SF3_FACT+1.3d-6) < 0.)call stop_model
     &      ('(SF3_FACT+1.3d-6) < 0 in master',255)
            SF3(I,J,L)=6.d0*(SF3_FACT+1.3d-6)*EXP(-1.d-7*colmO2**.35)
     &      *by35*SQRT(1.224d3*(cos(ABS(LAT_DG(J,1))*radian))**2.+1.d0)
            SF3(I,J,L)=SF3(I,J,L)*5.d-2
          else
            SF3(I,J,L)=0.d0
          endif
C SF2 is photlysis of NO in bands (0-0) and (1-0) based on Nicolet,
C Pl. Space Sci., p 111, 1980. SF2_fact is a ratio representative of
C bands (0-0) and (1-0); =  bin5_flux[present] / bin5_flux[1988] :
          if(L > maxl)then
            if(colmO2 > 2.d19)then
              SF2(I,J,L)=4.5d-6*EXP(-(1.d-8*colmO2**.38+5.d-19*colmO3))
            else
              SF2(I,J,L)=4.75d-6*EXP(-1.5d-20*colmO2)
            endif
            SF2(I,J,L)=SF2(I,J,L)*SF2_fact*
     &      by35*SQRT(1.224d3*(cos(ABS(LAT_DG(J,1))*radian))**2.+1.d0)
          else
            SF2(I,J,L)=0.d0
          endif
#endif
        END DO

       endif ! (sunlight)
      endif  ! (modphot)                           
      
CCCCCCCCCCCCCCCCC END PHOTOLYSIS SECTION CCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                 BEGIN CHEMISTRY                                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


c Calculate the chemical reaction rates:
      call Crates (I,J
#ifdef SHINDELL_STRAT_CHEM
     &                ,aero
#endif
     &                     )      


      if((ALB(I,J,1) /= 0.d0).AND.(sza < szamax))then
CCCCCCCCCCCCCCCCCCCC   SUNLIGHT   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCC FAMILY PARTITIONING CCCCCCCCCCCCCCCCCCCCCCCCCC

#ifdef SHINDELL_STRAT_CHEM
       call Oxinit(LM,I,J)
       call HOxfam(LM,I,J)
       call NOxfam(LM,I,J)
       call BrOxfam(LM,I,J)
#else
       call Oxinit(maxl,I,J)
       call HOxfam(maxl,I,J)
       call NOxfam(maxl,I,J)
#endif
CCCCCCCCCCCCCCCCC NON-FAMILY CHEMISTRY CCCCCCCCCCCCCCCCCCCCCCCC

      call chemstep(I,J,changeL)

C Save 3D radical arrays to pass to aerosol code:
      if(coupled_chem == 1) then
        do l=1,maxl
          oh_live(i,j,l)=y(nOH,L)
          no3_live(i,j,l)=yNO3(i,j,l)
        end do
      endif

#ifdef SHINDELL_STRAT_CHEM
      call ClOxfam(LM,I,J,ClOx_old) ! needed something from chemstep.
#endif

CCCCCCCCCCCCC PRINT SOME CHEMISTRY DIAGNOSTICS CCCCCCCCCCCCCCCC
      if(prnchg .and. J == jprn .and. I == iprn) then
       jay = (J >= J_0 .and. J <= J_1) 
       l=lprn
       write(out_line,*) ' '
       call write_parallel(trim(out_line),crit=jay)
       write(out_line,*) 'Family ratios at I,J,L: ',i,j,l
       call write_parallel(trim(out_line),crit=jay)
       write(out_line,*) 'OH/HO2 = ',y(nOH,l)/y(nHO2,l)
       call write_parallel(trim(out_line),crit=jay)
       write(out_line,*) 'O/O3 = ',y(nO,l)/y(nO3,l)
       call write_parallel(trim(out_line),crit=jay)
       write(out_line,*) 'O1D/O3 = ',y(nO1D,l)/y(nO3,l),
     &  '  J(O1D) = ',ss(2,l,I,J)
       call write_parallel(trim(out_line),crit=jay)
       write(out_line,*) 'NO/NO2 = ',y(nNO,l)/y(nNO2,l),
     &  '   J(NO2) = ',ss(1,l,I,J)
       call write_parallel(trim(out_line),crit=jay)
       write(out_line,*) 'conc OH = ',y(nOH,l)
       call write_parallel(trim(out_line),crit=jay)
#ifdef SHINDELL_STRAT_CHEM
       write(out_line,*) 'Cl,ClO,Cl2O2,OClO,Cl2 = ',y(nCl,l),
     &  y(nClO,l),y(nCl2O2,l),y(nOClO,l),y(nCl2,l)
       call write_parallel(trim(out_line),crit=jay)
       write(out_line,*) 'Br,BrO = ',y(nBr,l),y(nBrO,l)
       call write_parallel(trim(out_line),crit=jay)
       write(out_line,*) 'pCl,pClO,pOClO,pBrO = ',pClx(I,J,l),
     &  pClOx(I,J,l),pOClOx(I,J,l),pBrOx(I,J,l)
       call write_parallel(trim(out_line),crit=jay)
#endif
       write(out_line,*)
     & 'sun, SALBFJ,sza,I,J,Itime= ',ALB(I,J,1),sza,I,J,Itime
       call write_parallel(trim(out_line),crit=jay)
       do inss=1,JPPJ
        write(out_line,195) ' J',inss,ay(ks(inss)),' = ',
     &  (ss(inss,Lqq,I,J),Lqq=1,LS1-1)
        call write_parallel(trim(out_line),crit=jay)
       enddo
       write(out_line,196) ' RCloud',(RCLOUDFJ(Lqq,I,J),Lqq=1,LS1-1)
       call write_parallel(trim(out_line),crit=jay)
       write(out_line,196) ' Ozone ',(y(nO3,Lqq),Lqq=1,LS1-1)
       call write_parallel(trim(out_line),crit=jay)
       write(out_line,*) ' '
       call write_parallel(trim(out_line),crit=jay)
      endif
 195  format(a2,i2,1x,a8,a3,11(1x,e9.2))
 196  format(a7,9x,11(1x,e9.2))
CCCCCCCCCCCCCCCCCCCC END CHEM DIAG SECT CCCCCCCCCCCCCCCCCCCCCCC
      
      else

CCCCCCCCCCCCCCCCCCCC END SUNLIGHT CCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCC   DARKNESS  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C*****************************************************************
C               ABOUT: N2O5 sink on sulfate aerosol 
C                REACTION PROBABLITY FORMULATION:
C
C To evaluate 'k' for N2O5 + H2O -> 2HNO3, assume k = GAMMA*SA*v/4
C (Dentener and Crutzen, 1993). GAMMA = 0.1, SA = given,
C v = molecular velocity (cm/s) where v = SQRT (8.Kb.T / PI*M);
C Kb = 1.38062E-23; T = Temp (K); M = mass N2O5 (Kg)
C
C Off-line sulfate fields to run in uncoupled mode give SO4 in
C cm2/cm3 'surface area density'.
C
C On-line sulfate in coupled mode must be converted to aerosol
C surface (cm2 aerosol/cm3 air) via aerosol volume fraction
C (cm3 aerosol/cm3 air). Assume a monodispersed aerosol with diameter
C of 0.078 microns. Specific aerosol density = 1.1g/cm3 (1.7g/cm3 ?)
C See Dentener and Crutzen, 1993 for details. Mr[sulfate] = 96.0g;
C Mr[(NH4)HSO4] = 115.0gC
C*****************************************************************

CCCCCCCCCCCCCCCC NIGHTTIME  TROPOSPHERE  CCCCCCCCCCCCCCCCCCCCCC

      do L=1,maxl ! (troposphere)

        if (coupled_chem == 1) then
          ! Convert SO4 from mass (kg) to aerosol surface per grid box:
          fact_so4 = BYDXYP(J)*BYAM(L,I,J)*28.0D0*y(nM,L) ! *96/96
          ! sulfate(I,J,L)=(trm(I,J,L,n_SO4)*fact_so4*3.0D0)
     &    ! /(3.9D-6*avog*1.1D0)
          sulfate(I,J,L)=trm(I,J,L,n_SO4)*fact_so4*byavog*1.063636d7
        endif

        pfactor=dxyp(J)*AM(L,I,J)/y(nM,L)
        bypfactor=1.D0/pfactor
        RVELN2O5=SQRT(TX(I,J,L)*RKBYPIM)*100.d0
C       Calculate sulfate sink, and cap it at 90% of N2O5:
        wprod_sulf=
     &  dt2*sulfate(I,J,L)*y(n_N2O5,L)*RGAMMASULF*RVELN2O5*0.25d0
        if(wprod_sulf > 0.9d0*y(n_N2O5,L))wprod_sulf=0.9d0*y(n_N2O5,L)
        prod_sulf=wprod_sulf*pfactor
        TAJLS(J,L,jls_N2O5sulf)=TAJLS(J,L,jls_N2O5sulf)
     &  -(prod_sulf*bymass2vol(n_N2O5))

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

c       define O3 from Ox:
        y(nO3,L)=y(n_Ox,L)*pOx(I,J,L)

c       calculate change in NO3:
        dNO3=rr(7,L)*y(nNO2,L)*y(n_Ox,L)-(rr(24,L)*y(nNO2,L)+
     &       2.d0*rr(25,L)*yNO3(I,J,L))*yNO3(I,J,L) -(rr(36,L)
     &       *y(n_Alkenes,L)+rr(32,L)*y(n_Isoprene,L))*yNO3(I,J,L)
C       including DMS+NO3 :
     &       - ydms(i,j,l)*rsulf3(i,j,l)*yNO3(I,J,L)
#ifdef SHINDELL_STRAT_CHEM
        dNO3=dNO3-(rr(28,L)*y(n_HCHO,L)+rr(99,L)*y(nNO2,L))
     &       *yNO3(I,J,L)+rr(92,L)*y(n_N2O5,L)
#else
        dNO3=dNO3-(rr(28,L)*y(n_HCHO,L)+rr(52,L)*y(nNO2,L))
     &       *yNO3(I,J,L)+rr(46,L)*y(n_N2O5,L)
#endif
        dNO3=dNO3*dt2

c       limit the change in NO3:
        if(-dNO3 > 0.66d0*yNO3(I,J,L)) dNO3=-0.66d0*yNO3(I,J,L)

c       apply the NO3 change; limit the value to positive & 1/2 NOx:
        yNO3(I,J,L)=yNO3(I,J,L)+dNO3
        if(yNO3(I,J,L) < 0.d0) yNO3(I,J,L)=0.d0
        if(yNO3(I,J,L) > y(n_NOx,L)*0.5d0)
     &  yNO3(I,J,L)=y(n_NOx,L)*0.5d0

c       calculate and limit NO2:
        y(nNO2,L)=y(n_NOx,L)-yNO3(I,J,L)
        pNOx(I,J,L)=y(nNO2,L)/(y(n_NOx,L)+1.d-10)
        if(pNOx(I,J,L) > 1.d0) pNOx(I,J,L)=1.d0
        if(pNOx(I,J,L) < 0.5d0)pNOx(I,J,L)=0.5d0

C       LOWER LIMIT ON N2O5:
        if(y(n_N2O5,L) <= 1.d0) y(n_N2O5,L)=1.d0

C Calculate and limit gaseous changes to HNO3, HCHO, N2O5, Aldehyde,
C Alkenes, Isoprene, and AlkylNit:

        gwprodHNO3=(y(n_HCHO,L)*rr(28,L)+2.5d-15*
     &             yAldehyde(I,J,L))*yNO3(I,J,L)*dt2
        if(gwprodHNO3 > 0.5d0*y(n_NOx,L))gwprodHNO3=0.5d0*y(n_NOx,L)
        if(gwprodHNO3 > y(n_HCHO,L))gwprodHNO3=y(n_HCHO,L)
        gprodHNO3=gwprodHNO3*pfactor

#ifdef SHINDELL_STRAT_CHEM
        gwprodN2O5=(yNO3(I,J,L)*y(nNO2,L)*rr(99,L)-y(n_N2O5,L)
     &             *rr(92,L))*dt2
#else
        gwprodN2O5=(yNO3(I,J,L)*y(nNO2,L)*rr(52,L)-y(n_N2O5,L)
     &             *rr(46,L))*dt2
#endif
        if(gwprodN2O5 > 0.25d0*y(n_NOx,L))gwprodN2O5=0.25d0*y(n_NOx,L)
        if(-gwprodN2O5 > 0.5d0*y(n_N2O5,L))
     &       gwprodN2O5=-0.49d0*y(n_N2O5,L)

        changeAldehyde=(rr(36,L)*y(n_Alkenes,L)+rr(32,L)*
     &              y(n_Isoprene,L)*0.12d0-2.5d-15*yAldehyde(I,J,L))*
     &              yNO3(I,J,L)*dt2
        if(-changeAldehyde > 0.75d0*yAldehyde(I,J,L))changeAldehyde=
     &  -0.75d0*yAldehyde(I,J,L)

        changeAlkenes=(rr(32,L)*y(n_Isoprene,L)*0.45d0-rr(36,L)*
     &               y(n_Alkenes,L))*yNO3(I,J,L)*dt2+(rr(31,L)*
     &               y(n_Isoprene,L)-rr(35,L)*y(n_Alkenes,L))*
     &               y(nO3,L)*dt2
        if(-changeAlkenes > 0.75d0*y(n_Alkenes,L))changeAlkenes=
     &  -0.75d0*y(n_Alkenes,L)

        changeIsoprene=-(rr(32,L)*yNO3(I,J,L)
     &                 +rr(31,L)*y(nO3,L))*y(n_Isoprene,L)*dt2
        if(-changeIsoprene > 0.75d0*y(n_Isoprene,L))changeIsoprene=
     &  -0.75d0*y(n_Isoprene,L)

        changeHCHO=(rr(36,L)*y(n_Alkenes,L)+
     &             rr(32,L)*y(n_Isoprene,L)*0.03d0)*yNO3(I,J,L)*dt2
     &             -gwprodHNO3+(rr(31,L)*y(n_Isoprene,L)*0.9d0
     &             +rr(35,L)*y(n_Alkenes,L))*y(nO3,L)*0.64d0*dt2

        changeAlkylNit=rr(32,L)*y(n_Isoprene,L)*
     &                yNO3(I,J,L)*dt2*0.9d0
        if(-changeAlkylNit > 0.75d0*y(n_AlkylNit,L))changeAlkylNit=
     &  -0.75d0*y(n_AlkylNit,L)

c Convert some changes to molecules/cm3/s:
        changeHNO3=gwprodHNO3+2*wprod_sulf  !always positive
        changeNOx=-gwprodHNO3-2*gwprodN2O5-(0.9d0*rr(32,L)*
     &  y(n_Isoprene,L)+2.5d-15*yAldehyde(I,J,L))*yNO3(I,J,L)*dt2
        if(-changeNOx > y(n_NOx,L))changeNOx=-.95d0*y(n_NOx,L)
        changeN2O5=gwprodN2O5-wprod_sulf

c Ensure nitrogen conservation (presumably dNOx<0, others >0):
        rlossN=0.d0
        rprodN=0.d0
        if(changeNOx < 0.d0)then
          rlossN=rlossN+changeNOx
        else
          rprodN=rprodN+changeNOx
        endif
        if(changeHNO3 < 0.d0)then
          rlossN=rlossN+changeHNO3
        else
          rprodN=rprodN+changeHNO3
        endif
        if(changeN2O5 < 0.d0)then
          rlossN=rlossN+2*changeN2O5
        else
          rprodN=rprodN+2*changeN2O5
        endif
        if(changeAlkylNit < 0.d0)then
          rlossN=rlossN+changeAlkylNit
        else
          rprodN=rprodN+changeAlkylNit
        endif
        if(rprodN > rlossN)then
          ratioN=-rlossN/rprodN
          if(changeNOx > 0.d0)     changeNOx    =changeNOx   *ratioN
          if(changeHNO3 > 0.d0)    changeHNO3   =changeHNO3  *ratioN
          if(changeN2O5 > 0.d0)    changeN2O5   =changeN2O5  *ratioN
          if(changeAlkylNit > 0.d0)changeAlkylNit=
     &                                           changeAlkylNit*ratioN
        else
          ratioN=rprodN/(-rlossN)
          if(changeNOx < 0.d0)     changeNOx    =changeNOx   *ratioN
          if(changeHNO3 < 0.d0)    changeHNO3   =changeHNO3  *ratioN
          if(changeN2O5 < 0.d0)    changeN2O5   =changeN2O5  *ratioN
          if(changeAlkylNit < 0.d0)changeAlkylNit=
     &                                           changeAlkylNit*ratioN
        endif
#ifdef TRACERS_HETCHEM
C       Include reactions on dust for HNO3:
        changeHNO3 = changeHNO3 - krate(i,j,l,1,1)*y(n_HNO3,l)*dt2
        changeN_d1 = krate(i,j,l,2,1) * y(n_HNO3,l) *dt2
        changeN_d2 = krate(i,j,l,3,1) * y(n_HNO3,l) *dt2
        changeN_d3 = krate(i,j,l,4,1) * y(n_HNO3,l) *dt2
#endif

C Apply Alkenes, AlkyNit, and Aldehyde changes here:
        y(n_Alkenes,L)  =y(n_Alkenes,L)  +changeAlkenes
        y(n_AlkylNit,L) =y(n_AlkylNit,L) +changeAlkylNit
        yAldehyde(I,J,L)=yAldehyde(I,J,L)+changeAldehyde

C Note: the lower limit of 1 placed on the resulting tracer mass
C from the following changes is to prevent negative tracer mass:

C -- HCHO --
c       Gas phase NO3 + HCHO -> HNO3 + CO yield of HCHO & CO:
        changeL(L,n_HCHO)=changeHCHO*pfactor*bymass2vol(n_HCHO)
        if(-changeL(L,n_HCHO) > trm(I,J,L,n_HCHO))then
          changeL(L,n_HCHO)=-.95d0*trm(I,J,L,n_HCHO)
          changeHCHO=changeL(L,n_HCHO)*mass2vol(n_HCHO)*bypfactor
        endif
        IF((trm(i,j,l,n_HCHO)+changeL(l,n_HCHO)) < 1.d0) THEN
          changeL(l,n_HCHO) = 1.d0 - trm(i,j,l,n_HCHO)
          changeHCHO=changeL(L,n_HCHO)*mass2vol(n_HCHO)*bypfactor
        ENDIF
        wprodHCHO=changeHCHO
C -- CO --
        changeL(L,n_CO)=gprodHNO3*bymass2vol(n_CO)
        IF((trm(i,j,l,n_CO)+changeL(l,n_CO)) < 1.d0)
     &  changeL(l,n_CO) = 1.d0 - trm(i,j,l,n_CO)
        wprodCO=gwprodHNO3   ! <<< note
        if(changeL(L,n_CO) >= 0.) then  
          TAJLS(J,L,jls_COp)=TAJLS(J,L,jls_COp)+changeL(L,n_CO)
        else
          TAJLS(J,L,jls_COd)=TAJLS(J,L,jls_COd)+changeL(L,n_CO)
        endif       
C -- HNO3 --  (HNO3 from gas and het phase rxns )
        changeL(L,n_HNO3)=changeHNO3*pfactor*bymass2vol(n_HNO3)
        IF((trm(i,j,l,n_HNO3)+changeL(l,n_HNO3)) < 1.d0) THEN
          changeL(l,n_HNO3) = 1.d0 - trm(i,j,l,n_HNO3)
          changeHNO3=changeL(L,n_HNO3)*mass2vol(n_HNO3)*bypfactor
        END IF
#ifdef TRACERS_HETCHEM
        changeL(L,n_N_d1)=changeN_d1*pfactor*bymass2vol(n_N_d1)
        if(i==36.and.j==28.and.l==1) then
          write(out_line,*)'Mchange L 2 ', changeL(L,n_N_d1),changeN_d1
          call write_parallel(trim(out_line),crit=.true.)
        endif
        IF((trm(i,j,l,n_N_d1)+changeL(l,n_N_d1)) < 1.d0) THEN
          changeL(l,n_N_d1) = 1.d0 - trm(i,j,l,n_N_d1)
          changeN_d1=changeL(L,n_N_d1)*mass2vol(n_N_d1)*bypfactor
        END IF
        changeL(L,n_N_d2)=changeN_d2*pfactor*bymass2vol(n_N_d2)
        IF((trm(i,j,l,n_N_d2)+changeL(l,n_N_d2)) < 1.d0) THEN
          changeL(l,n_N_d2) = 1.d0 - trm(i,j,l,n_N_d2)
          changeN_d2=changeL(L,n_N_d2)*mass2vol(n_N_d2)*bypfactor
        END IF
        changeL(L,n_N_d3)=changeN_d3*pfactor*bymass2vol(n_N_d3)
        IF((trm(i,j,l,n_N_d3)+changeL(l,n_N_d3)) < 1.d0) THEN
          changeL(l,n_N_d3) = 1.d0 - trm(i,j,l,n_N_d3)
          changeN_d3=changeL(L,n_N_d3)*mass2vol(n_N_d3)*bypfactor
        END IF
#endif
C -- N2O5 --  (N2O5 from gas and het phase rxns)
        changeL(L,n_N2O5)=changeN2O5*pfactor*bymass2vol(n_N2O5)
        IF((trm(i,j,l,n_N2O5)+changeL(l,n_N2O5)) < 1.d0) THEN
          changeL(l,n_N2O5) = 1.d0 - trm(i,j,l,n_N2O5)
          changeN2O5=changeL(L,n_N2O5)*mass2vol(n_N2O5)*bypfactor
        END IF
c -- NOx --   (NOx from gas phase rxns)
        changeL(L,n_NOx)=changeNOx*pfactor*bymass2vol(n_NOx)
        IF((trm(i,j,l,n_NOx)+changeL(l,n_NOx)) < 1.d0) THEN
          changeL(l,n_NOx) = 1.d0 - trm(i,j,l,n_NOx)
          changeNOx=changeL(L,n_NOx)*mass2vol(n_NOx)*bypfactor
        END IF
C -- Alkenes --  (Alkenes from gas phase rxns)
        changeL(L,n_Alkenes)=
     &  changeAlkenes*pfactor*bymass2vol(n_Alkenes)
        IF((trm(i,j,l,n_Alkenes)+changeL(l,n_Alkenes)) < 1.d0)THEN
          changeL(l,n_Alkenes) = 1.d0 - trm(i,j,l,n_Alkenes)
          changeAlkenes=changeL(L,n_Alkenes)*mass2vol(n_Alkenes)
     &    *bypfactor
        END IF
c -- Isoprene -- (Isoprene from gas phase rxns)
        changeL(L,n_Isoprene)=
     &  changeIsoprene*pfactor*bymass2vol(n_Isoprene)
        IF((trm(i,j,l,n_Isoprene)+changeL(l,n_Isoprene)) < 1.d0)
     &  THEN
          changeL(l,n_Isoprene) = 1.d0 - trm(i,j,l,n_Isoprene)
          changeIsoprene=changeL(L,n_Isoprene)*mass2vol(n_Isoprene)
     &    *bypfactor
        END IF
c -- AlkylNit -- (AlkylNit from gas phase rxns)
        changeL(L,n_AlkylNit)=
     &  changeAlkylNit*pfactor*bymass2vol(n_AlkylNit)
        IF((trm(i,j,l,n_AlkylNit)+changeL(l,n_AlkylNit)) < 1.d0)
     &  THEN
          changeL(l,n_AlkylNit) = 1.d0 - trm(i,j,l,n_AlkylNit)
          changeAlkylNit=changeL(L,n_AlkylNit)*mass2vol(n_AlkylNit)
     &    *bypfactor
        END IF

C Save 3D radical arrays to pass to aerosol code:
C Make sure we get the nightime values; Set OH to zero for now:
        if(coupled_chem == 1) then
          oh_live(i,j,l)=0.d0
          no3_live(i,j,l)=yNO3(i,j,l)
        endif

CCCCCCCCCCCCC PRINT SOME CHEMISTRY DIAGNOSTICS CCCCCCCCCCCCCCCC
        if(prnchg.and.J == jprn.and.I == iprn.and.L == lprn)then
          jay = (J >= J_0 .and. J <= J_1)
          write(out_line,*)
     &    'dark, SALBFJ,sza,I,J,L,Itime= ',ALB(I,J,1),sza,I,J,L,Itime
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(n_NOx),': ',
     &    changeNOx,' molecules produced; ',
     &    100.d0*(changeNOx)/y(n_NOx,L),' percent of'
     &    ,y(n_NOx,L),'(',1.d9*y(n_NOx,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(n_HNO3),': ',
     &    changeHNO3,' molecules produced; ',
     &    100.d0*(changeHNO3)/y(n_HNO3,L),' percent of'
     &    ,y(n_HNO3,L),'(',1.d9*y(n_HNO3,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
#ifdef TRACERS_HETCHEM
          write(out_line,198) ay(n_HNO3),': ',
     &    (-krate(i,j,l,1,1)*y(n_HNO3,l)*dt2),' molecules dest dust ',
     &    (100.d0*(-krate(i,j,l,1,1)*y(n_HNO3,l)*dt2))/y(n_HNO3,L),
     &    ' percent of'
     &    ,y(n_HNO3,L),'(',1.d9*y(n_HNO3,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
#endif
          write(out_line,198) ay(n_N2O5),': ',
     &    changeN2O5,' net molec produced; ',
     &    100.d0*(changeN2O5)/y(n_N2O5,L),' percent of'
     &    ,y(n_N2O5,L),'(',1.d9*y(n_N2O5,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(n_N2O5),': ',
     &    gwprodN2O5,' molec prod fm gas;  ',
     &    100.d0*(gwprodN2O5)/y(n_N2O5,L),' percent of'
     &    ,y(n_N2O5,L),'(',1.d9*y(n_N2O5,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(n_HCHO),': ',
     &    wprodHCHO,' molecules produced; ',
     &    100.d0*(wprodHCHO)/y(n_HCHO,L),' percent of'
     &    ,y(n_HCHO,L),'(',1.d9*y(n_HCHO,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) 'Aldehyde',': ',
     &    changeAldehyde,' molecules produced; ',
     &    100.d0*(changeAldehyde)/yAldehyde(I,J,L),' percent of'
     &    ,yAldehyde(I,J,L),'(',1.d9*yAldehyde(I,J,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) 'Alkenes ',': ',
     &    changeAlkenes,' molecules produced; ',
     &    100.d0*(changeAlkenes)/y(n_Alkenes,L),' percent of'
     &    ,y(n_Alkenes,L),'(',1.d9*y(n_Alkenes,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) 'Isoprene',': ',
     &    changeIsoprene,' molecules produced; ',
     &    100.d0*(changeIsoprene)/y(n_Isoprene,L),' percent of'
     &    ,y(n_Isoprene,L),'(',1.d9*y(n_Isoprene,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) 'AlkylNit',': ',
     &    changeAlkylNit,' molecules produced; ',
     &    100.d0*(changeAlkylNit)/y(n_AlkylNit,L),' percent of'
     &    ,y(n_AlkylNit,L),'(',1.d9*y(n_AlkylNit,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,199) 'NO2, NO3  = ',y(nNO2,L),yNO3(I,J,L)
          call write_parallel(trim(out_line),crit=jay)
        endif
 198    format(1x,a8,a2,e13.3,a21,f10.0,a11,2x,e13.3,3x,a1,f12.5,a6)
 199    format(1x,a20,2(2x,e13.3))
CCCCCCCCCCCCCCCCCCCC END CHEM DIAG SECT CCCCCCCCCCCCCCCCCCCCCCC

C Make sure nighttime chemistry changes are not too big:
        error=.false.
        if(changeNOx < -1.d15.OR.changeNOx > 1.d15) then
          write(6,*) 'Big chg@ Itime,I,J,L,NOx ',Itime,I,J,L,changeNOx
          error=.true.
        end if
        if(changeHNO3 < -1.d15.OR.changeHNO3 > 1.d15) then
          write(6,*) 'Big chg@ Itime,I,J,L,HNO3',Itime,I,J,L,changeHNO3
          error=.true.
        end if
        if(changeN2O5 < -1.d15.OR.changeN2O5 > 1.d15) then
          write(6,*) 'Big chg@ Itime,I,J,L,N2O5',Itime,I,J,L,changeN2O5
          error=.true.
        end if
        if(wprodHCHO < -1.d15.OR.wprodHCHO > 1.d15) then
          write(6,*)'Big chg@ Itime,I,J,L,HCHO',Itime,I,J,L,wprodHCHO
          error=.true.
        endif
        if(error)call stop_model('nighttime chem: big changes',255)

C       ACCUMULATE 3D NO3 diagnostic: 
        if (yNO3(I,J,L) > 0.d0 .and. yNO3(I,J,L) < 1.d20)
     &  taijs(i,j,ijs_NO3(l))=taijs(i,j,ijs_NO3(l))+yNO3(i,j,l)
     
       enddo  ! troposphere loop

CCCCCCCCCCCCCCCC END NIGHTTIME TROPOSPHERE CCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCC NIGHTTIME STRATOSPHERE CCCCCCCCCCCCCCCCCCCCCCC

#ifdef SHINDELL_STRAT_CHEM
       do L=maxl+1,LM ! (stratosphere)

         pfactor=dxyp(J)*AM(L,I,J)/y(nM,L)
         wprod_sulf=DT2*y(n_N2O5,L)*rr(105,L) !rxn on sulfate & PSCs
         if(wprod_sulf >= 0.5d0*y(n_N2O5,L))wprod_sulf=0.5d0*y(n_N2O5,L)

c        calculate and limit change in NO3:
         dNO3=rr(7,L)*y(nNO2,L)*y(n_Ox,L)-(rr(24,L)*y(nNO2,L)+
     &        2.d0*rr(25,L)*yNO3(I,J,L))*yNO3(I,J,L)
         dNO3=dNO3-(rr(28,L)*y(n_HCHO,L)+rr(99,L)*y(nNO2,L))
     &        *yNO3(I,J,L)+rr(92,L)*y(n_N2O5,L)
         dNO3=dNO3*dt2
         if(-dNO3 > 0.66d0*yNO3(I,J,L))dNO3=-0.66d0*yNO3(I,J,L)

c        apply the NO3 change; limit value to positive and 1/2 NOx:
         yNO3(I,J,L)=yNO3(I,J,L)+dNO3
         if(yNO3(I,J,L) < 0.d0)yNO3(I,J,L)=0.d0
         if(yNO3(I,J,L) > y(n_NOx,L)*0.5d0)yNO3(I,J,L)=
     &   y(n_NOx,L)*0.5d0

c        calculate and limit NO2:
         y(nNO2,L)=y(n_NOx,L)-yNO3(I,J,L)
         pNOx(I,J,L)=y(nNO2,L)/(y(n_NOx,L)+1.d-10)
         if(pNOx(I,J,L) > 1.d0)pNOx(I,J,L)=1.d0
         if(pNOx(I,J,L) < 0.5d0)pNOx(I,J,L)=0.5d0

C        LOWER LIMIT ON N2O5:
         if(y(n_N2O5,L) <= 1.d0)y(n_N2O5,L)=1.d0

C Calculate and limit gaseous changes to HNO3, HCHO, N2O5, Aldehyde,
C Alkenes, Isoprene, and AlkylNit:
         gwprodHNO3=y(n_HCHO,L)*rr(28,L)*yNO3(I,J,L)*dt2
         if(gwprodHNO3 > 0.5d0*y(n_NOx,L))gwprodHNO3=0.5d0*y(n_NOx,L)
         if(gwprodHNO3 > y(n_HCHO,L))gwprodHNO3=y(n_HCHO,L)
         gprodHNO3=gwprodHNO3*pfactor
         gwprodN2O5=(yNO3(I,J,L)*y(nNO2,L)*rr(99,L)-y(n_N2O5,L)
     &              *rr(92,L))*dt2
         if(gwprodN2O5 > 0.25d0*y(n_NOx,L))gwprodN2O5=0.25d0*y(n_NOx,L)
         if(-gwprodN2O5 > 0.5d0*y(n_N2O5,L))
     &   gwprodN2O5=-0.49d0*y(n_N2O5,L)
         changeClONO2=y(nClO,L)*rr(103,L)*y(nNO2,L)*dt2
         if(changeClONO2 >= y(nClO,L))changeClONO2=0.8d0*y(nClO,L)
         if(changeClONO2 >= y(nNO2,L))changeClONO2=0.8d0*y(nNO2,L)
         changeClOx=-changeClONO2

c Convert some changes to molecules/cm3/s:(w/o ClONO2->HNO3 het or PSC)
         changeHNO3=gwprodHNO3+2.d0*wprod_sulf  !always positive
         changeNOx=-gwprodHNO3-2.d0*gwprodN2O5
         if(-changeNOx > y(n_NOx,L))changeNOx=-.95d0*y(n_NOx,L)

         changeN2O5=gwprodN2O5-wprod_sulf

c Ensure nitrogen conservation (presumably dNOx<0, others >0):
         rlossN=0.d0
         rprodN=0.d0
         if(changeNOx < 0.d0)then
           rlossN=rlossN+changeNOx
         else
           rprodN=rprodN+changeNOx
         endif
         if(changeHNO3 < 0.d0)then
           rlossN=rlossN+changeHNO3
         else
           rprodN=rprodN+changeHNO3
         endif
         if(changeN2O5 < 0.d0)then
           rlossN=rlossN+2.d0*changeN2O5
         else
           rprodN=rprodN+2.d0*changeN2O5
         endif
         if(rprodN > rlossN)then
           ratioN=-rlossN/rprodN
           if(changeNOx > 0.d0)   changeNOx     =changeNOx     *ratioN
           if(changeHNO3 > 0.d0)  changeHNO3    =changeHNO3    *ratioN
           if(changeN2O5 > 0.d0)  changeN2O5    =changeN2O5    *ratioN
         else
           ratioN=rprodN/(-rlossN)
           if(changeNOx < 0.d0)   changeNOx     =changeNOx     *ratioN
           if(changeHNO3 < 0.d0)  changeHNO3    =changeHNO3    *ratioN
           if(changeN2O5 < 0.d0)  changeN2O5    =changeN2O5    *ratioN
         endif

         changeNOx=changeNOx-changeClONO2
         if(-changeNOx > y(n_NOx,L))changeNOx=-.95d0*y(n_NOx,L)

c Heterogeneous reaction ClONO2+H2O on sulfate (and PSCs if present):
         if(rr(106,L) > 2.d-35)then
           changehetClONO2=-(rr(106,L)*y(n_ClONO2,L))*dt2
           if(changehetClONO2 >= 0.5*y(n_ClONO2,L))changehetClONO2=
     &     -0.5d0*y(n_ClONO2,L)
           changeClONO2=changeClONO2+changehetClONO2
           changeHOCl=-changehetClONO2
           changeHNO3=changeHNO3-changehetClONO2
         else
           changeHOCl=0.d0
         endif
         
         changeHCl=0.d0

c Polar Stratospheric Clouds (PSC) Chemistry:
c 105 N2O5    +H2O     -->HNO3    +HNO3  (calculated above)
c 106 ClONO2  +H2O     -->HOCl    +HNO3  (calculated above)
c 107 ClONO2  +HCl     -->Cl      +HNO3  !really makes Cl2
c 108 HOCl    +HCl     -->Cl      +H2O   !raeally makes Cl2
c 109 N2O5    +HCl     -->Cl      +HNO3  !really makes ClNO2  2

         pscX=0.d0
         if((pres2(l) <= 150.d0 .and. pres2(l) >= 35.d0) ! certain pres
     &   .and.((lat_dg(J,1)<=-66.).OR.(lat_dg(J,1)>=66.))! near poles
     &   .and. (ta(l) <= T_thresh))then                  ! cold enough
           pscX=1.d0
           chg107=rr(107,L)*y(n_ClONO2,L)*dt2
           if(chg107 >= 0.4d0*y(n_ClONO2,L))chg107=0.4d0*y(n_ClONO2,L)
           if(chg107 >= 0.3d0*y(n_HCl,L))chg107=0.3d0*y(n_HCl,L)
           chg108=rr(108,L)*y(n_HOCl,L)*dt2
           if(chg108 >= 0.8d0*y(n_HOCl,L))chg108=0.8d0*y(n_HOCl,L)
           if(chg108 >= 0.3d0*y(n_HCl,L))chg108=0.3d0*y(n_HCl,L)
           chg109=rr(109,L)*y(n_N2O5,L)*dt2
           if(chg109 >= 0.5d0*y(n_N2O5,L))chg109=0.5d0*y(n_N2O5,L)
           if(chg109 >= 0.3d0*y(n_HCl,L))chg109=0.3d0*y(n_HCl,L)
           changeClONO2=changeClONO2-chg107
           changeHOCl=changeHOCl-chg108
           changeN2O5=changeN2O5-chg109
           changeHCl=changeHCl-chg107-chg108-chg109
           changeHNO3=changeHNO3+chg107+chg109
           changeClOx=changeClOx+chg107+chg108+chg109
c          Note that really the last 3 produce Cl2, not ClOx, and Cl2
C          at night is stable and doesn't go back into ClONO2, so
C          should eventually keep track of Cl2/ClOx partitioning!
         endif

c Remove some of the HNO3 formed heterogeneously, as it doesn't come
c back to the gas phase:
         if(pres2(L) < 245.d0 .and. pres2(L) >= 31.6d0 .and.
     &   (J <= JSS+1.or.J >= JNN+1) .and. ta(L) <= T_thresh)then
           changeHNO3 = changeHNO3 - 3.0d-3*y(n_HNO3,L)
         endif
         if(aero(L)  /=  0)
     &   changeHNO3 = changeHNO3 -1.2d-4*y(n_HNO3,L)
         if(aero(L)  /=  0.and.pres2(L) <= 50.d0.and.pres2(L) > 38.d0)
     &   changeHNO3 = changeHNO3 +9.6d-5*y(n_HNO3,L)
         if(aero(L)  /=  0.and.pres2(L) < 245.d0.and.pres2(L) > 50.d0)
     &   changeHNO3 = changeHNO3 - 1.5d-5*(pres2(L)-30.d0)*y(n_HNO3,L)

C Make sure nighttime chemistry changes are not too big:
         error=.false.
         if(changeNOx < -1.d15.OR.changeNOx > 1.d15) then
           write(*,*) 'Big chg@ Itime,I,J,L,NOx ',Itime,I,J,L,changeNOx
           error=.true.
         end if
         if(changeHNO3 < -1.d15.OR.changeHNO3 > 1.d15) then
           write(*,*) 'Big chg@ Itime,I,J,L,HNO3',Itime,I,J,L,changeHNO3
           error=.true.
         end if
         if(changeN2O5 < -1.d15.OR.changeN2O5 > 1.d15) then
           write(*,*) 'Big chg@ Itime,I,J,L,N2O5',Itime,I,J,L,changeN2O5
           error=.true.
         end if
         if(wprodHCHO < -1.d15.OR.wprodHCHO > 1.d15) then
           write(*,*)'Big chg@ Itime,I,J,L,HCHO',Itime,I,J,L,wprodHCHO
           error=.true.
         endif
         if(error)call stop_model('nighttime chem: big changes',255)

C -- HNO3 --  (HNO3 from gas and het phase rxns )
         changeL(L,n_HNO3)=changeHNO3*pfactor*bymass2vol(n_HNO3)
         IF((trm(i,j,l,n_HNO3)+changeL(l,n_HNO3)) < 1.d0) THEN
           changeL(l,n_HNO3) = 1.d0 - trm(i,j,l,n_HNO3)
           changeHNO3=changeL(L,n_HNO3)*mass2vol(n_HNO3)*bypfactor
         END IF
C -- N2O5 --  (N2O5 from gas and het phase rxns)
         changeL(L,n_N2O5)=changeN2O5*pfactor*bymass2vol(n_N2O5)
         IF((trm(i,j,l,n_N2O5)+changeL(l,n_N2O5)) < 1.d0) THEN
           changeL(l,n_N2O5) = 1.d0 - trm(i,j,l,n_N2O5)
           changeN2O5=changeL(L,n_N2O5)*mass2vol(n_N2O5)*bypfactor
         END IF
c -- NOx --   (NOx from gas phase rxns)
         changeL(L,n_NOx)=changeNOx*pfactor*bymass2vol(n_NOx)
         IF((trm(i,j,l,n_NOx)+changeL(l,n_NOx)) < 1.d0) THEN
           changeL(l,n_NOx) = 1.d0 - trm(i,j,l,n_NOx)
           changeNOx=changeL(L,n_NOx)*mass2vol(n_NOx)*bypfactor
         END IF
c --  Ox --   ( Ox from gas phase rxns)
         if(pres2(l) > 1.0.and.pres2(l) < 25.0)then
           changeOx=-(rr(7,L)*y(nNO2,L) + y(n_H2O2,L)*0.7d0*1.4d-11*
     &     exp(-2000./TA(L)))*y(n_Ox,L)*dt2*2.5d0
           changeL(L,n_Ox)=changeOx*pfactor*bymass2vol(n_Ox)
           IF((trm(i,j,l,n_Ox)+changeL(l,n_Ox)) < 1.d0) THEN
             changeL(l,n_Ox) = 1.d0 - trm(i,j,l,n_Ox)
             changeOx=changeL(L,n_Ox)*mass2vol(n_Ox)*bypfactor
           END IF
           if(changeL(L,n_Ox) >= 0.) then  
             TAJLS(J,L,jls_Oxp)=TAJLS(J,L,jls_Oxp)+changeL(L,n_Ox)
           else
             TAJLS(J,L,jls_Oxd)=TAJLS(J,L,jls_Oxd)+changeL(L,n_Ox)
           endif 
         endif
c -- ClONO2 --   (ClONO2 from gas and het phase rxns)
         changeL(L,n_ClONO2)=changeClONO2*pfactor*
     &   bymass2vol(n_ClONO2)
         IF((trm(i,j,l,n_ClONO2)+changeL(l,n_ClONO2)) < 1.d0) THEN
           changeL(l,n_ClONO2) = 1.d0 - trm(i,j,l,n_ClONO2)
           changeClONO2=changeL(L,n_ClONO2)*mass2vol(n_ClONO2)*
     &     bypfactor
         END IF
c -- ClOx --   (ClOx from gas and het phase rxns)
         changeL(L,n_ClOx)=changeClOx*pfactor*bymass2vol(n_ClOx)
         IF((trm(i,j,l,n_ClOx)+changeL(l,n_ClOx)) < 1.d0) THEN
           changeL(l,n_ClOx) = 1.d0 - trm(i,j,l,n_ClOx)
           changeClOx=changeL(L,n_ClOx)*mass2vol(n_ClOx)*bypfactor
         END IF
         if(pscX > 0.9d0)then
c -- HOCl --   (HOCl from het phase rxns)
           changeL(L,n_HOCl)=changeHOCl*pfactor*bymass2vol(n_HOCl)
           IF((trm(i,j,l,n_HOCl)+changeL(l,n_HOCl)) < 1.d0) THEN
             changeL(l,n_HOCl) = 1.d0 - trm(i,j,l,n_HOCl)
             changeHOCl=changeL(L,n_HOCl)*mass2vol(n_HOCl)*bypfactor
           END IF
c -- HCl --   (HCl from het phase rxns)
           changeL(L,n_HCl)=changeHCl*pfactor*bymass2vol(n_HCl)
           IF((trm(i,j,l,n_HCl)+changeL(l,n_HCl)) < 1.d0) THEN
             changeL(l,n_HCl) = 1.d0 - trm(i,j,l,n_HCl)
             changeHCl=changeL(L,n_HCl)*mass2vol(n_HCl)*bypfactor
           END IF
         endif  ! pscX > 0.9

CCCCCCCCCCCCC PRINT SOME CHEMISTRY DIAGNOSTICS CCCCCCCCCCCCCCCC
         if(prnchg .and. J == jprn .and. I == iprn .and. L == lprn)then
           jay = (J >= J_0 .and. J <= J_1)
           write(out_line,*) 'dark, SALBFJ,sza,I,J,L,Itime= ',
     &     ALB(I,J,1),sza,I,J,L,Itime
           call write_parallel(trim(out_line),crit=jay)
           if(pscX > 0.9d0)then
             write(out_line,*) 'There are PSCs, T =',ta(L)
             call write_parallel(trim(out_line),crit=jay)
           else
             write(out_line,*) 'There are no PSCs, T =',ta(L)
             call write_parallel(trim(out_line),crit=jay)
           endif
           write(out_line,198) ay(n_NOx),': ',
     &     changeNOx,' molecules produced; ',
     &     100.d0*(changeNOx)/y(n_NOx,L),' percent of'
     &     ,y(n_NOx,L),'(',1.d9*y(n_NOx,L)/y(nM,L),' ppbv)'
           call write_parallel(trim(out_line),crit=jay)
           write(out_line,198) ay(n_HNO3),': ',
     &     changeHNO3,' molecules produced; ',
     &     100.d0*(changeHNO3)/y(n_HNO3,L),' percent of'
     &     ,y(n_HNO3,L),'(',1.d9*y(n_HNO3,L)/y(nM,L),' ppbv)'
           call write_parallel(trim(out_line),crit=jay)
           write(out_line,198) ay(n_N2O5),': ',
     &     changeN2O5,' net molec produced; ',
     &     100.d0*(changeN2O5)/y(n_N2O5,L),' percent of'
     &     ,y(n_N2O5,L),'(',1.d9*y(n_N2O5,L)/y(nM,L),' ppbv)'
           call write_parallel(trim(out_line),crit=jay)
           write(out_line,198) ay(n_N2O5),': ',
     &     gwprodN2O5,' molec prod fm gas;  ',
     &     100.d0*(gwprodN2O5)/y(n_N2O5,L),' percent of'
     &     ,y(n_N2O5,L),'(',1.d9*y(n_N2O5,L)/y(nM,L),' ppbv)'
           call write_parallel(trim(out_line),crit=jay)
           write(out_line,198) ay(n_N2O5),': ',
     &     -wprod_sulf,' molec prod fm sulf; ',
     &     -100.d0*(wprod_sulf)/y(n_N2O5,L),' percent of'
     &     ,y(n_N2O5,L),'(',1.d9*y(n_N2O5,L)/y(nM,L),' ppbv)'
           call write_parallel(trim(out_line),crit=jay)
           write(out_line,198) ay(n_ClONO2),': ',
     &     changeClONO2,' molecules produced; ',
     &     100.d0*(changeClONO2)/y(n_ClONO2,L),' percent of'
     &     ,y(n_ClONO2,L),'(',1.d9*y(n_ClONO2,L)/y(nM,L),' ppbv)'
           call write_parallel(trim(out_line),crit=jay)
           write(out_line,198) ay(n_ClOx),': ',
     &     changeClOx,' molecules produced; ',
     &     100.d0*(changeClOx)/y(n_ClOx,L),' percent of'
     &     ,y(n_ClOx,L),'(',1.d9*y(n_ClOx,L)/y(nM,L),' ppbv)'
           call write_parallel(trim(out_line),crit=jay)
           write(out_line,198) ay(n_HOCl),': ',
     &     changeHOCl,' molecules produced; ',
     &     100.d0*(changeHOCl)/y(n_HOCl,L),' percent of'
     &     ,y(n_HOCl,L),'(',1.d9*y(n_HOCl,L)/y(nM,L),' ppbv)'
           call write_parallel(trim(out_line),crit=jay)
           write(out_line,198) ay(n_HCl),': ',
     &     changeHCl,' molecules produced; ',
     &     100.d0*(changeHCl)/y(n_HCl,L),' percent of'
     &     ,y(n_HCl,L),'(',1.d9*y(n_HCl,L)/y(nM,L),' ppbv)'
           call write_parallel(trim(out_line),crit=jay)
           write(out_line,199) 'NO2, NO3  = ',y(nNO2,L),yNO3(I,J,L)
           call write_parallel(trim(out_line),crit=jay)
           write(out_line,198) ay(n_Ox),': ',
     &     changeOx,' molecules produced; ',
     &     100.d0*(changeOx)/y(n_Ox,L),' percent of'
     &     ,y(n_Ox,L),'(',1.d9*y(n_Ox,L)/y(nM,L),' ppbv)'
           call write_parallel(trim(out_line),crit=jay)
         endif
CCCCCCCCCCCCCCCCCCCC END CHEM DIAG SECT CCCCCCCCCCCCCCCCCCCCCCC

         if (y(nClO,L) > 0.d0 .and. y(nClO,L) < 1.d20)
     &   TAJLS(J,L,jls_ClOcon)=TAJLS(J,L,jls_ClOcon)+y(nClO,L)/y(nM,L)
         if (y(nH2O,L) > 0.d0 .and. y(nH2O,L) < 1.d20)
     &   TAJLS(J,L,jls_H2Ocon)=TAJLS(J,L,jls_H2Ocon)+y(nH2O,L)/y(nM,L)

       enddo ! end stratosphere loop
CCCCCCCCCCCCCCCC END NIGHTTIME STRATOSPHERE CCCCCCCCCCCCCCCCCCC
#endif
      endif
CCCCCCCCCCCCCCCCCCCC END DARKNESS CCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#ifdef SHINDELL_STRAT_CHEM
      LL=LM
#else
      LL=maxl
#endif
      DO L=1,LL  

C Lower limit on HO2NO2 : 
        if(trm(i,j,l,n_HO2NO2)+changeL(l,n_HO2NO2) < 1.d0)
     &  changeL(l,n_HO2NO2) = 1.d0 - trm(i,j,l,n_HO2NO2)

#ifdef SHINDELL_STRAT_CHEM
c Tropospheric halogen sink Br & Cl :
        IF (L <= maxl) THEN 
          rmv=max(0.d0,    ! remove between 0 and 
     &        min(0.99d0,  ! 99% only...
     &        2.389d-1 * LOG(pres2(L))  - 1.05d0
     &        ))           ! about 5% at 100mb and 60% at 1000mb
          changeL(L,n_ClOx)  =-trm(I,J,L,n_ClOx)  *rmv
          changeL(L,n_HCl)   =-trm(I,J,L,n_HCl)   *rmv
          changeL(L,n_HOCl)  =-trm(I,J,L,n_HOCl)  *.85d0
          changeL(L,n_ClONO2)=-trm(I,J,L,n_ClONO2)*rmv
          changeL(L,n_BrOx)  =-trm(I,J,L,n_BrOx)  *0.9d0
          changeL(L,n_HBr)   =-trm(I,J,L,n_HBr)   *0.9d0
          changeL(L,n_HOBr)  =-trm(I,J,L,n_HOBr)  *0.9d0
          changeL(L,n_BrONO2)=-trm(I,J,L,n_BrONO2)*0.9d0
        END IF

c Set CLTOT based on CFCs (3.3 ppbv yield from complete oxidation of
c 1.8 ppbv CFC plus 0.5 ppbv background) :
        if(L > maxl)then
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !WARNING: RESETTING SOME Y's HERE; SO DON'T USE THEM BELOW!     
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          y(n_ClOx,L)=(trm(I,J,L,n_ClOx)+changeL(L,n_ClOx))*y(nM,L)*
     &    mass2vol(n_ClOx)*BYDXYP(J)*BYAM(L,I,J)
          y(n_HCl,L)= (trm(I,J,L,n_HCl)+changeL(L,n_HCl))*y(nM,L)*
     &    mass2vol(n_HCl)*BYDXYP(J)*BYAM(L,I,J)
          y(n_HOCl,L)=(trm(I,J,L,n_HOCl)+changeL(L,n_HOCl))*y(nM,L)*
     &    mass2vol(n_HOCl)*BYDXYP(J)*BYAM(L,I,J)
          y(n_ClONO2,L)=(trm(I,J,L,n_ClONO2)+changeL(L,n_ClONO2))*
     &    y(nM,L)*mass2vol(n_ClONO2)*BYDXYP(J)*BYAM(L,I,J)

          CLTOT=((y(n_CFC,1)/y(nM,1)-y(n_CFC,L)/y(nM,L))*(3.3d0/1.8d0)*
     &    y(n_CFC,1)/(1.8d-9*y(nM,1)))
          CLTOT=CLTOT+0.5d-9
          CLTOT=CLTOT*y(nM,L)/
     &    (y(n_ClOx,L)+y(n_HCl,L)+y(n_HOCl,L)+y(n_ClONO2,L))
          if(prnchg.and.J == jprn.and.I == iprn.and.L == lprn)then  
            write(out_line,66) CLTOT
            call write_parallel(trim(out_line),crit=jay)
 66         format ('CLTOT = ',F20.5)
          endif
          IF(CLTOT <= 0.999d0 .OR. CLTOT >= 1.001d0) THEN
            changeL(L,n_ClOx)=changeL(L,n_ClOx)*CLTOT+
     &      trm(I,J,L,n_ClOx)*(CLTOT-1.D0)
            changeL(L,n_HCl)=changeL(L,n_HCl)*CLTOT+
     &      trm(I,J,L,n_HCl)*(CLTOT-1.D0)
            changeL(L,n_HOCl)=changeL(L,n_HOCl)*CLTOT+
     &      trm(I,J,L,n_HOCl)*(CLTOT-1.D0)
c           Conserve N wrt ClONO2 once inital Cl changes past:
            if(Itime-ItimeI >= 6)then
              changeL(L,n_NOx)=changeL(L,n_NOx)-
     &        (trm(I,J,L,n_ClONO2)+changeL(L,n_ClONO2))*
     &        (CLTOT-1.D0)*tr_mm(n_NOx)/tr_mm(n_ClONO2)
              if(-changeL(L,n_NOx) > trm(I,J,L,n_NOx))changeL(L,n_NOx)=
     &        -0.8d0*trm(I,J,L,n_NOx)
            endif
            changeL(L,n_ClONO2)=changeL(L,n_ClONO2)*CLTOT+
     &      trm(I,J,L,n_ClONO2)*(CLTOT-1.D0)
          ENDIF

c Set Total Bromine(using CLTOT name!) based on CFCs (1.0 pptv yield
C from complete oxidation of 1.8 ppbv CFC plus 0.5 pptv background) :
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !WARNING: RESETTING SOME Y's HERE; SO DON'T USE THEM BELOW!     
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          y(n_BrOx,L)=(trm(I,J,L,n_BrOx)+changeL(L,n_BrOx))*y(nM,L)*
     &    mass2vol(n_BrOx)*BYDXYP(J)*BYAM(L,I,J)
          y(n_HBr,L)= (trm(I,J,L,n_HBr)+changeL(L,n_HBr))*y(nM,L)*
     &    mass2vol(n_HBr)*BYDXYP(J)*BYAM(L,I,J)
          y(n_HOBr,L)=(trm(I,J,L,n_HOBr)+changeL(L,n_HOBr))*y(nM,L)*
     &    mass2vol(n_HOBr)*BYDXYP(J)*BYAM(L,I,J)
          y(n_BrONO2,L)=(trm(I,J,L,n_BrONO2)+changeL(L,n_BrONO2))*
     &    y(nM,L)*mass2vol(n_BrONO2)*BYDXYP(J)*BYAM(L,I,J)
     
          CLTOT=((y(n_CFC,1)/y(nM,1)-y(n_CFC,L)/y(nM,L))*(1.0d-3/1.8d0)
     &    *y(n_CFC,1)/(1.8d-9*y(nM,1)))
          CLTOT=CLTOT+0.5d-12
          CLTOT=CLTOT*y(nM,L)/
     &    (y(n_BrOx,L)+y(n_HBr,L)+y(n_HOBr,L)+y(n_BrONO2,L))
          if(prnchg.and.J == jprn.and.I == iprn.and.L == lprn)then  
            write(out_line,67) CLTOT
            call write_parallel(trim(out_line),crit=jay)
 67         format ('BrTOT = ',F20.5)
          endif
          IF(CLTOT <= 0.999d0 .OR. CLTOT >= 1.001d0) THEN
            changeL(L,n_BrOx)=changeL(L,n_BrOx)*CLTOT+
     &      trm(I,J,L,n_BrOx)*(CLTOT-1.D0)
            changeL(L,n_HBr)=changeL(L,n_HBr)*CLTOT+
     &      trm(I,J,L,n_HBr)*(CLTOT-1.D0)
            changeL(L,n_HOBr)=changeL(L,n_HOBr)*CLTOT+
     &      trm(I,J,L,n_HOBr)*(CLTOT-1.D0)
c           Conserve N wrt BrONO2 once inital Br changes past:
            if(Itime-ItimeI >= 6)then
              changeL(L,n_NOx)=changeL(L,n_NOx)-
     &        (trm(I,J,L,n_BrONO2)+changeL(L,n_BrONO2))*
     &        (CLTOT-1.D0)*tr_mm(n_NOx)/tr_mm(n_BrONO2)
              if(-changeL(L,n_NOx) > trm(I,J,L,n_NOx))changeL(L,n_NOx)=
     &        -0.8d0*trm(I,J,L,n_NOx)
            endif
            changeL(L,n_BrONO2)=changeL(L,n_BrONO2)*CLTOT+
     &      trm(I,J,L,n_BrONO2)*(CLTOT-1.D0)
          ENDIF
        endif ! L > maxl
#endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Save chemistry changes for applying in apply_tracer_3Dsource.  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        DO N=1,NTM_CHEM
          tr3Dsource(i,j,l,nChemistry,n) = changeL(l,n) * bydtsrc
        END DO
#ifdef TRACERS_HETCHEM
        tr3Dsource(i,j,l,nChemistry,n_N_d1) = changeL(l,n_N_d1) *bydtsrc
        tr3Dsource(i,j,l,nChemistry,n_N_d2) = changeL(l,n_N_d2) *bydtsrc
        tr3Dsource(i,j,l,nChemistry,n_N_d3) = changeL(l,n_N_d3) *bydtsrc
#endif

      END DO ! end current altitude loop

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      END DO ! >>>> MAIN I LOOP ENDS <<<<

      END DO ! >>>> MAIN J LOOP ENDS <<<<
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCC!$OMP END PARALLEL DO  

#ifdef SHINDELL_STRAT_CHEM
      if(prnchg)then
       if(LM /= 23) call stop_model('alter O3 map for new LM',255)
       DO J=J_0,J_1; DO I=1,IMAXJ(J)
         ss27(1:LM,I,J)=ss(27,1:LM,I,J)
       ENDDO; ENDDO
       CALL PACK_COLUMN(grid,ss27,ss27_glob)
       IF(AM_I_ROOT()) THEN
         write(6,*) 'Map of O3 production from O2 (Herz & SRB)'
         if(which_trop == 0)write(6,*)
     &   'NOTE: lower limit of strat is actually LTROPO(I,J), however!'
         write(6,'(a4,7(i10))') 'Jqq:',(Jqq,Jqq=3,44,6)
         do Lqq=LM,LS1,-1 ! inconvenient to print down to LTROPO(I,J)
           do jqq=1,JM
             photO2_glob(Jqq,Lqq)=0.d0
             do Iqq=1,IMAXJ(jqq)
               photO2_glob(Jqq,Lqq)=photO2_glob(Jqq,Lqq)
     &         +2.d0*ss27_glob(Lqq,Iqq,Jqq) ! makes 2 Os
             enddo
           enddo
           write(6,'(i2,7(1x,E10.3))') 
     &     Lqq,(photO2_glob(Jqq,Lqq)/REAL(IMAXJ(Jqq)),Jqq=3,44,6)
         enddo
       END IF ! end AM_I_ROOT
      endif
#endif

CCCCCCCCCCCCCCCCCC END CHEMISTRY SECTION CCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                 BEGIN OVERWRITING                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C If fix_CH4_chemistry is turned on, reset the CH4 tracer everywhere
C to initial conditions (in mixing ratio units) :
      if(fix_CH4_chemistry == 1) call get_CH4_IC(1)
   
C If desired, use correction factors on stratospheric Ox:   
      if(correct_strat_Ox) then
        imonth= 1
        DO m=2,12
          IF((JDAY <= MDOFM(m)).AND.(JDAY > MDOFM(m-1))) THEN
            imonth=m
            EXIT
          END IF
        END DO
        if(Itime == ItimeI) then
          write(out_line,*)'Warning: Please remember that Ox'
     &    //' stratospheric correction factors are on and change with'
     &    //' time.'
          call write_parallel(trim(out_line))
        end if
      end if

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Special cases of overwriting, when doing stratospheric chemistry C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
#ifdef SHINDELL_STRAT_CHEM      
C N2O, CFC, and optional CH4 L=1 overwriting: with all these "fact"s
C this looks complicated, but basically, you are either converting
C from mixing ratio to KG (normal case) or from cm*atm to KG  
C (interactive radiation case - for more on that conversion, see
C the notes on O3MULT in the TRCHEM_Shindell_COM program):
      PIfact(:)=1.d0     
      if(PI_run == 1)then
        PIfact(n_NOx)=PIratio_N
        if(use_rad_n2o > 0 .or. use_rad_cfc > 0) then
          write(out_line,*)
     &    'warning: use_rad_{cfc,n2o} overrides PIfact({cfc,n2o})' 
          call write_parallel(trim(out_line))
        endif     
        if(use_rad_n2o == 0) PIfact(n_N2O)=PIratio_N2O
        if(use_rad_cfc == 0) PIfact(n_CFC)=PIratio_CFC
      endif
      fact2=n2o_pppv  ! default N2O mixing ratio overwrite
      fact3=cfc_pppv  ! default CFC mixing ratio overwrite
      fact7=fact_cfc
      if(use_rad_cfc == 0)fact7=1.d0
      do j=J_0,J_1
       fact6=2.69d20*dxyp(j)*byavog
       fact5=fact6 
       fact4=fact6
       do i=1,IMAXJ(j)
        fact1=bymair*am(1,i,j)*dxyp(j)
        if(use_rad_n2o == 0)fact4=fact1 
        if(use_rad_cfc == 0)fact5=fact1
        if(use_rad_n2o > 0)fact2=rad_to_chem(1,i,j,3)
        if(use_rad_cfc > 0)fact3=rad_to_chem(1,i,j,5)
        tr3Dsource(i,j,1,nChemistry,n_N2O)=0.d0
        tr3Dsource(i,j,1,nChemistry,n_CFC)=0.d0
        tr3Dsource(i,j,1,nStratwrite,n_N2O)=(fact2*fact4*
     &  tr_mm(n_N2O)*PIfact(n_N2O) - trm(i,j,1,n_N2O))*bydtsrc
        tr3Dsource(i,j,1,nStratwrite,n_CFC)=(fact3*fact5*fact7*
     &  tr_mm(n_CFC)*PIfact(n_CFC) - trm(i,j,1,n_CFC))*bydtsrc
        if(use_rad_ch4 > 0)then
          tr3Dsource(i,j,1,nChemistry,n_CH4)=0.d0
          tr3Dsource(i,j,1,nStratwrite,n_CH4)=(rad_to_chem(1,i,j,4)*
     &    fact6*tr_mm(n_CH4)-trm(i,j,1,n_CH4))*bydtsrc
        endif
       end do
      end do

C For Ox, NOx, BrOx, and ClOx, we have overwriting where P < 0.1mb:
      !(Interpolate BrOx & ClOx altitude-dependence to model resolution)
      CALL LOGPINT(LCOalt,PCOalt,BrOxaltIN,LM,PRES2,BrOxalt,.true.)
      CALL LOGPINT(LCOalt,PCOalt,ClOxaltIN,LM,PRES2,ClOxalt,.true.)
      do L=LS1,LM
       if(pres2(L) < 0.1)then
        do j=J_0,J_1         
          do i=1,IMAXJ(j)
            ! -- Ox --
            tr3Dsource(i,j,L,nChemistry,n_Ox)=0.d0
            if(correct_strat_Ox) then
             tr3Dsource(i,j,L,nStratwrite,n_Ox)=(rad_to_chem(L,i,j,1)*
     &       dxyp(j)*O3MULT*corrOx(J,L,imonth)-trm(i,j,L,n_Ox))*bydtsrc
            else
             tr3Dsource(i,j,L,nStratwrite,n_Ox)=(rad_to_chem(L,i,j,1)*
     &       dxyp(j)*O3MULT                   -trm(i,j,L,n_Ox))*bydtsrc
            end if
            ! -- ClOx --
            tr3Dsource(i,j,L,nChemistry,n_ClOx)=0.d0
            tr3Dsource(i,j,L,nStratwrite,n_ClOx)=(1.d-11*ClOxalt(l)
     &      *tr_mm(n_ClOx)*bymair*am(L,i,j)*dxyp(j)
     &      -trm(i,j,L,n_ClOx))*bydtsrc    
            ! -- BrOx --
            tr3Dsource(i,j,L,nChemistry,n_BrOx)=0.d0
            tr3Dsource(i,j,L,nStratwrite,n_BrOx)=(1.d-11*BrOxalt(l)
     &      *tr_mm(n_BrOx)*bymair*am(L,i,j)*dxyp(j)
     &      -trm(i,j,L,n_BrOx))*bydtsrc
            ! -- NOx --
            if(PIfact(n_NOx) /= 1.) ! what's this for? I forget.
     &      tr3Dsource(i,j,L,nChemistry,n_NOx)=0.d0 !75=1*300*2.5*.1
            tr3Dsource(i,j,L,nStratwrite,n_NOx)=(75.d-11
     &      *am(L,i,j)*dxyp(j)*PIfact(n_NOx)-trm(i,j,L,n_NOx))*bydtsrc 
          end do ! I 
        end do   ! J
       end if    ! pressure
      end do     ! L

#else

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C No stratospheric chemistry; overwrite all tracers in stratosphere C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C W A R N I N G : This section should never be used if there is
C chemistry done in the stratosphere, because the stratospheric
C changes below assume that the changeL variable is zero for L>maxl
C at this point in the code. Not the case if chemistry was done.
C To put it another way, the overwritings below are explicitly 
C functions of tracer mass UNCHANGED by chemistry !

C determine pre-industrial factors, if any:
      PIfact(:)=1.d0
      if(PI_run == 1) then
        do N=1,NTM
          select case(trname(n))
          case('NOx','HNO3','N2O5','HO2NO2')
            PIfact(n)=PIratio_N
          case('CO')
            PIfact(n)=PIratio_CO_S
          case('PAN','Isoprene','AlkylNit','Alkenes','Paraffin')
            PIfact(n)=PIratio_other
          end select
        end do
      endif

C Calculate an average tropical CH4 value near 569 hPa::
      CH4_569_part(:)=0.d0   
      count_569_part(:)=0.d0
      DO J=J_0,J_1
        if(LAT_DG(J,1) >= -30. .and. LAT_DG(J,1) <= 30.)then
          do I=1,IMAXJ(J)
            count_569_part(J)=count_569_part(J)+1.d0
            CH4_569_part(J)=CH4_569_part(J)+1.d6*bydxyp(j)*
     &      (F569M*trm(I,J,L569M,n_CH4)*byam(L569M,I,J)+
     &      F569P*trm(I,J,L569P,n_CH4)*byam(L569P,I,J))
          end do
        end if
      END DO
      CALL GLOBALSUM(grid,  CH4_569_part,  CH4_569, all=.true.)
      CALL GLOBALSUM(grid,count_569_part,count_569, all=.true.)
      if(count_569 <= 0.)call stop_model('count_569.le.0',255)
      CH4_569 = CH4_569 / count_569

      r179=1.d0/1.79d0 ! 1.79 is observed trop. CH4

      do j=J_0,J_1
       J3=MAX(1,NINT(float(j)*float(JCOlat)*BYFJM))! index for CO
       do i=1,IMAXJ(J)
         select case(which_trop)
         case(0); maxl=ltropo(I,J)
         case(1); maxl=ls1-1
         case default; call stop_model('which_trop problem 5',255)
         end select
         changeL(:,:)=0.d0 ! initilize the change
         
         do L=maxl+1,LM    ! >> BEGIN LOOP OVER STRATOSPHERE <<
c         Update stratospheric ozone to amount set in radiation:
          if(correct_strat_Ox) then
            changeL(L,n_Ox)=rad_to_chem(L,I,J,1)*DXYP(J)*O3MULT
     &      *corrOx(J,L,imonth) - trm(I,J,L,n_Ox)
          else
            changeL(L,n_Ox)=rad_to_chem(L,I,J,1)*DXYP(J)*O3MULT
     &                          - trm(I,J,L,n_Ox)
          end if
          byam75=F75P*byam(L75P,I,J)+F75M*byam(L75M,I,J)
          FACT1=2.0d-9*DXYP(J)*am(L,I,J)*byam75
C         We think we have too little stratospheric NOx, so, to
C         increase the flux into the troposphere, increasing
C         previous stratospheric value by 70% here: GSF/DTS 9.15.03: 
          changeL(L,n_NOx)=trm(I,J,L,n_Ox)*2.3d-4*1.7d0*PIfact(n_NOx)
     &                                         - trm(I,J,L,n_NOx)
          changeL(L,n_N2O5)=  FACT1*PIfact(n_N2O5)- trm(I,J,L,n_N2O5)
          changeL(L,n_HNO3)=trm(I,J,L,n_Ox)*4.2d-3*PIfact(n_HNO3)
     &                                         - trm(I,J,L,n_HNO3)
          changeL(L,n_H2O2)=  FACT1            - trm(I,J,L,n_H2O2)
          changeL(L,n_CH3OOH)=FACT1            - trm(I,J,L,n_CH3OOH)
          changeL(L,n_HCHO)=  FACT1            - trm(I,J,L,n_HCHO)
          ! here 70. = 1.4E-7/2.0E-9
          changeL(L,n_HO2NO2)=FACT1*70.d0*PIfact(n_HO2NO2)
     &                                         - trm(I,J,L,n_HO2NO2)
          changeL(L,n_CO)=(COlat(J3)*COalt(L)*1.d-9*bymass2vol(n_CO)
     &    *AM(L,I,J)*DXYP(J))*0.4d0*PIfact(n_CO)-trm(I,J,L,n_CO)!<<40%
          changeL(L,n_PAN)     = FACT1*1.d-4*PIfact(n_PAN)
     &                           - trm(I,J,L,n_PAN)
          changeL(L,n_Isoprene)= FACT1*1.d-4*PIfact(n_Isoprene)
     &                           - trm(I,J,L,n_Isoprene)
          changeL(L,n_AlkylNit)= FACT1*1.d-4*PIfact(n_AlkylNit)
     &                           - trm(I,J,L,n_AlkylNit)
          changeL(L,n_Alkenes) = FACT1*1.d-4*PIfact(n_Alkenes)
     &                           - trm(I,J,L,n_Alkenes)
          changeL(L,n_Paraffin)= FACT1*1.d-4*PIfact(n_Paraffin)
     &                           - trm(I,J,L,n_Paraffin)

c Overwrite stratospheric ch4 based on HALOE obs for tropics and
C extratropics and scale by the ratio of near-569hPa mixing ratios
C to 1.79:
          if (fix_CH4_chemistry == 0) then ! -------------------------
          CH4FACT=CH4_569*r179
          IF((J <= JS).OR.(J > JN)) THEN               ! extratropics
            DO L2=L,maxl+1,-1
              IF(CH4altX(L2) /= 0.d0) THEN
                CH4FACT=CH4FACT*CH4altX(L2)
                EXIT
              END IF
            END DO
          ELSE IF((J > JS).AND.(J <= JN)) THEN         ! tropics
            DO L2=L,maxl+1,-1
              IF(CH4altT(L2) /= 0.d0) THEN
                CH4FACT=CH4FACT*CH4altT(L2)
                EXIT
              END IF
            END DO
          END IF
          select case(PI_run)
          case(1) ! preindustrial
            if(j <= jm/2)then; PIfact(n_CH4)=pfix_CH4_S
            else             ; PIfact(n_CH4)=pfix_CH4_N
            end if
            changeL(L,n_CH4)=am(l,i,j)*dxyp(j)*TR_MM(n_CH4)*bymair
     &      *PIfact(n_CH4)  - trm(I,J,L,n_CH4)
          case default
            ! also ensure that strat overwrite is only a sink:
            changeL(L,n_CH4)=-MAX(0d0,trm(I,J,L,n_CH4)-
     &      (AM(L,I,J)*DXYP(J)*CH4FACT*1.d-6))
          end select
          else   ! ------ fixed CH4 set in get_CH4_IC(1) -------------
            changeL(L,n_CH4)=0.d0
          end if ! ---------------------------------------------------

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Save overwrite changes for applying in apply_tracer_3Dsource.  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

          DO N=1,NTM_CHEM
            tr3Dsource(i,j,l,nStratwrite,n) = changeL(l,n)*bydtsrc
          END DO

        end do ! >> END LOOP OVER STRATOSPHERE <<
       end do  ! i
      end do   ! j
#endif

CCCCCCCCCCCCCCCCCC END OVERWRITE SECTION CCCCCCCCCCCCCCCCCCCCCC

c Save new tracer Ox field for use in radiation or elsewhere:
      do j=J_0,J_1
        DU_O3(J)=0.d0 ! Drew's diagnostic...
        do i=1,imaxj(j) 
#ifdef SHINDELL_STRAT_CHEM
         maxl = LM
#else
         select case(which_trop)
         case(0); maxl=ltropo(I,J)
         case(1); maxl=ls1-1
         case default; call stop_model('maxl error',255)
         end select
#endif
         do l=1,maxl
           o3_tracer_save(l,i,j)=(trm(i,j,l,n_Ox) +
     &     (tr3Dsource(i,j,l,nChemistry,n_Ox) + 
     &     tr3Dsource(i,j,l,nStratwrite,n_Ox))*dtsrc)
     &     *bydxyp(j)*byO3MULT
           DU_O3(J)=DU_O3(J)+o3_tracer_save(l,i,j)
         end do
         if(maxl < LM) then
           do l=maxl+1,LM
             o3_tracer_save(l,i,j)=rad_to_chem(l,i,j,1)
             DU_O3(J)=DU_O3(J)+o3_tracer_save(l,i,j)
           end do
         end if
        end do ! i
        DU_O3(J)=1.d3*DU_O3(J)/IMAXJ(J)
      end do   ! j
      
      if(MOD(Itime,24) == 0)then
       call PACK_DATA( grid, DU_O3, DU_O3_glob )
       IF(AM_I_ROOT()) THEN
         write(6,*) 'Ozone column fm -90 to +90'
         write(6,'(46(f4.0,1x))') (DU_O3_glob(J),J=1,JM)
       END IF
      endif

      RETURN
      END SUBROUTINE masterchem



      SUBROUTINE Crates(I,J
#ifdef SHINDELL_STRAT_CHEM
     &   ,aero
#endif
     & )
!@sum Crates calculate chemical reaction rates for each altitude,
!@+   using JPL 00.  Includes special calculations for pressure
!@+   dependent reactions. Specifically:
!@+   #13 CO+OH->HO2+CO2, #15 HO2+HO2->H2O2+O2, #16 OH+HNO3->H2O+NO3,
!@+   and reactions #29, and #42.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on masterchem000_M23p)

C**** GLOBAL parameters and variables:
      USE MODEL_COM, only: LM,JM,LS1,JEQ,ptop,psf,sig,Itime,ItimeI
      USE RAD_COM, only  : rad_to_chem
      USE GEOM, only     : LAT_DG
      USE CONSTANT, only : PI
      USE DYNAMICS, only : LTROPO
      USE TRCHEM_Shindell_COM, only: nr2,nr3,nmm,nhet,ta,ea,rr,pe,
     &               r1,sb,nst,y,nM,nH2O,ro,sn,which_trop,T_thresh

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
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
      REAL*8                :: byta,dd,pp,fw,rkp,rk2,rk3M,rrrr,temp
      INTEGER               :: L,jj,nb,Ltop
      INTEGER, INTENT(IN)   :: I,J
#ifdef SHINDELL_STRAT_CHEM
!@var PRES local nominal pressure
!@var LAXt,LAXb lowest and highest levels to have nonzero 
!@+   RAD-code aerosol extinction 
!@var aero array =1 for nonzero rkext, otherwise 0.
      REAL*8, DIMENSION(LM) :: PSCEX,PRES,rkext
      INTEGER               :: LAXt,LAXb
      INTEGER, INTENT(OUT), dimension(LM) :: aero
#endif

#ifdef SHINDELL_STRAT_CHEM
      aero(:)=0
      Ltop=LM
      PRES(:)=SIG(:)*(PSF-PTOP)+PTOP
      rkext(:)=0.d0 ! initialize over L
      if(rad_to_chem(1,I,J,2) /= 0.)call stop_model('kext prob 0',255)
      do L=2,Ltop
        if(rad_to_chem(L,I,J,2) /= 0..and.rad_to_chem(L-1,I,J,2) == 0.)
     &  LAXb=L
        if(rad_to_chem(L,I,J,2) == 0..and.rad_to_chem(L-1,I,J,2) /= 0.)
     &  LAXt=L-1
      end do
#else
      select case(which_trop)
      case(0); Ltop=ltropo(I,J)
      case(1); Ltop=ls1-1
      case default; call stop_model('which_trop problem 6',255)
      end select
#endif
      do L=1,Ltop            !  >>> BEGIN ALTITUDE LOOP <<<
        byta=1.d0/ta(L)
        do jj=1,nr2             ! bimolecular rates start
          IF(ea(jj) /= 0.d0) THEN
            rr(jj,L)=pe(jj)*exp(-ea(jj)*byta)
          ELSE
            rr(jj,L)=pe(jj)
          END IF
c         for #13, k=pe*(1+0.6*(Patm/1013)) Patm=[M]*(T*1.38E-19)
          if(jj == 13) rr(jj,L) =
     &    pe(jj)*(1.d0+0.6d0*((y(nM,L)*ta(L)*1.38d-19)/1013.d0))
c         for reaction #15, k=(kc+kp)fw, kc=rr
          if(jj == 15)then
            rkp=1.7d-33*y(nM,L)*exp(1000.d0*byta)
            fw=(1.d0+1.4d-21*y(nH2O,L)*exp(2200.d0*byta))
            rr(jj,L)=(rr(jj,L)+rkp)*fw
          endif
c         for #16, k=[pe*exp(-e(jj)/ta(l))]+k3[M]/(1+k3[M]/k2)
          if(jj == 16)then
            rk3M=y(nM,l)*1.90d-33*exp(725.d0*byta)
            rk2=4.10d-16*exp(1440.d0*byta)
            rr(jj,L)=rr(jj,L)+rk3M/(1.d0+(rk3M/rk2))
          endif
          if(jj == 29)rr(jj,L)=rr(jj,L)/y(nM,L)!PAN+M really PAN
          if(jj == 42)rr(jj,L)=rr(jj,L)/y(nM,L)!ROR+M really ROR
        end do                ! bimolecular rates end

#ifdef SHINDELL_STRAT_CHEM
        rr(87,L)=rr(87,L)*1.76d0  ! N2O+O(1D)-->NO+NO 
#endif
        do jj=1,nr3           ! trimolecular rates start
          rr(nr2+jj,L)=y(nM,L)*ro(jj)*(300.d0*byta)**sn(jj)
          if(sb(jj) >= 0.01d0)then 
            dd=rr(nr2+jj,L)/(r1(jj)*(300.d0*byta)**sb(jj))
            pp=0.6d0**(1.d0/(1.d0+(log10(dd))**2.))
            rr(nr2+jj,L)=(rr(nr2+jj,L)/(1.d0+dd))*pp
          endif
        end do                ! trimolecular rates end

        nb=nr2-nmm
        if(nmm >= 1) then
          do jj=1,nmm         ! monomolecular rates start
           ! 0.5 for precision,correct following line:
           rrrr=exp(0.5d0*ea(jj+nb)*byta)
           rr(jj+nb,L)=rr(nst(jj),L)/(rrrr*pe(jj+nb)*rrrr*y(nM,l))     
          end do              ! monomolecular rates end
        end if

#ifdef SHINDELL_STRAT_CHEM
c Calculate rates for heterogeneous reactions (Divided by solid
C in Chem1). sticking coefficients from JPL '02:
c       1=N2O5 + H2O --> 2HNO3          gamma=0.2 (aero),   0.0004 (PSC)
c       2=ClONO2 + H2O --> HOCl + HNO3  gamma=1.8d-4 (aero), 0.004 (PSC)
c       3=ClONO2 + HCl --> Cl2 + HNO3   gamma=0.2
c       4=HOCl + HCl --> Cl2 + H2O      gamma=0.1
c       5=N2O5 + HCl --> ClNO2 + HNO3   gamma=0.003
C
C Aerosols (14-33 km) & PSCs 14-22 km.
C
c Aerosol profiles and latitudinal distribution of extinction 
c coefficients(in km**-1) are from SAGE II data on GISS web site:

        if(pres(l) >= 245.d0 .or. pres(l) <= 5.d0)then 
          do jj=nr2+nr3+1,nr2+nr3+nhet
            rr(jj,L)=1.0d-35
          enddo 
          CYCLE
        else  
          if((pres(l) < 245.d0.and.pres(l) > 150.d0) .or. 
     &    LAXb < 1.or.LAXb > ltop.or.LAXt < 1.or.LAXt > ltop)then 
            rkext(l)=0.d0
          else
            if(pres(l) <= 150..and.pres(l) > 31.60)then
              if(l < LAXb) then
                rkext(l)=5.d-2*rad_to_chem(LAXb,i,j,2)
              else if(l > LAXt) then
                rkext(l)=0.33d0*rkext(l-1)
              else
                rkext(l)=5.d-2*rad_to_chem(l,i,j,2)
              endif
            endif
            if(pres(l) <= 31.6d0.and.pres(l) >= 17.8d0)then
              if(l < LAXb) then
                call stop_model('kext problem 1',255)
              else if(l > LAXt) then
                rkext(l)=2.0d0*rkext(l-1)
              else
                rkext(l)=5.d-2*rad_to_chem(l,i,j,2)
              endif
            endif
            if(pres(l) <= 17.8d0.and.pres(l) >= 10.0d0)then
              if(l < LAXb) then
                call stop_model('kext problem 2',255)
              else if(l > LAXt) then
                rkext(l)=16.d0*8.33333d-2*rkext(l-1)
              else
                rkext(l)=5.d-2*rad_to_chem(l,i,j,2)
              endif
            endif
            if(pres(l) <= 10.0d0.and.pres(l) >= 4.6d0)then
              if(l < LAXb) then
                call stop_model('kext problem 3',255)
              else if(l > LAXt) then
                rkext(l)=0.4d0*6.6667d-1*rkext(l-1)
              else
                rkext(l)=0.5d-2*rad_to_chem(l,i,j,2)
              endif
            endif
          endif

          IF(PRES(L) < 90. .and. J > JS .and. J <= JN)
     &    rkext(l)=rkext(l)*0.1d0   !<<<<<<<<<<< NOTE <<<<<<<<<<<
   
          rkext(l)=rkext(l)*2.0d0   !<<<<<<<<<<< NOTE <<<<<<<<<<<
           
          if(rkext(l) /= 0.)aero(l) = 1

c         PSCs: (pressure < 245mb and pressure > 5mb here):
          !NAT PSC surface conc per unit volume (cm^2/cm^3)
          pscEx(l)=0.d0
          if(pres(l) >= 31.6d0.and.ta(l) <= T_thresh .and.
     &    ((lat_dg(J,1) <= -66.).OR.(lat_dg(J,1) >= 66.)))then
            pscEx(l)=2.d-6 
            if(ta(l) <= T_thresh-7.d0)pscEx(l)=pscEx(l)*1.d1
          endif

c         Reaction 1 on sulfate and PSCs:      
          temp=sqrt(8.d0*1.38d-16*ta(l)*6.02d23/(PI*108.d0))
          rr(nr2+nr3+1,L)=0.5d0*rkext(l)*1.d-5*temp*0.2d0
          if(pres(l) > 31.6d0) rr(nr2+nr3+1,L)=
     &    rr(nr2+nr3+1,L)+0.25d0*pscEx(l)*temp*0.0008d0

c         Reaction 2 on sulfate and PSCs:      
          temp=sqrt(8.d0*1.38d-16*ta(l)*6.02d23/(PI*97.d0))
          rr(nr2+nr3+2,L)=0.5d0*rkext(l)*1.d-5*temp*1.8d-4
          if(pres(l) > 31.6d0) rr(nr2+nr3+2,L)=
     &    rr(nr2+nr3+2,L)+0.25d0*pscEx(l)*temp*8.d-3 

          if(pres(l) > 31.6d0) then
            rr(nr2+nr3+3,L)=0.25d0*pscEx(l)*temp*0.2d0
            rr(nr2+nr3+4,L)=
     &      sqrt(8.d0*1.38d-16*ta(l)*6.02d23/(PI*52.d0))
            rr(nr2+nr3+4,L)=0.25d0*pscEx(l)*rr(nr2+nr3+4,L)*0.1d0
            rr(nr2+nr3+5,L)=
     &      sqrt(8.d0*1.38d-16*ta(l)*6.02d23/(PI*108.d0))
            rr(nr2+nr3+5,L)=0.25d0*pscEx(l)*rr(nr2+nr3+5,L)*0.003d0
          endif
        endif  
#endif
      end do                  !  >>> END ALTITUDE LOOP <<<
 
      RETURN
      END SUBROUTINE Crates



#ifdef TRACERS_SPECIAL_Shindell
      SUBROUTINE checktracer(I,J)
!@sum checktracer for various debugging of tracer chemistry
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on masterchem000_M23p)

C**** GLOBAL parameters and variables:
      USE MODEL_COM, only  : Itime, LM, LS1
      USE DYNAMICS, only   : LTROPO
      USE TRACER_COM, only : ntm, trname, n_Ox, ntm_chem
#ifdef regional_Ox_tracers
     & ,NregOx
#endif
      USE TRCHEM_Shindell_COM, only: y, nM, which_trop

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var I,J passed horizontal position indicies
!@var L,igas dummy loop variables
!@var tlimit if tracer goes above this limit, model stops
!@var checkOx logical: should I check for large tropospheric Ox?
!@var checkmax logical: should I check for large tracers throughout?
!@var checkNeg logical: should I check for negative tracers?
!@var checkNaN logical: should I check for unreal tracers?
!@var nlast either ntm or ntm-NregOx
!@var maxL LTROPO(I,J) or LS1-1, depending upon which_trop variable

      INTEGER, PARAMETER :: nlast=
#ifdef regional_Ox_tracers
     & ntm_chem-NregOx
#else
     & ntm_chem
#endif
      INTEGER                  :: L, igas, maxL
      INTEGER, INTENT(IN)      :: I,J
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
    
      LOGICAL :: checkOx, checkmax, checkNeg, checkNan
      DATA checkNeg /.true./
      DATA checkNan /.true./
      DATA checkOx  /.true./
      DATA checkmax /.false./

      IF(i == 1.and.j == 1)
     & WRITE(6,*) 'WARNING: checktracer call is active.'
      select case(which_trop)
      case(0); maxl=ltropo(I,J)
      case(1); maxl=ls1-1
      case default; call stop_model('which_trop problem 7',255)
      end select

      IF(checkmax)
     &call stop_model('checktracer: set tlimit for tracers 11->25',255)
C     please (re)set tlimit values for tracers 11 through 15 in the
C     data statement above. Then, please delete the above stop.

C check if ozone gets really big:
       IF(checkOx) THEN
       do L=1,maxL
         if(y(n_Ox,L)/y(nM,L) > 1.d-5) then
           write(6,*)'Ox @ I,J,L,Ox,Itime:',I,J,L,y(n_Ox,L),Itime
           call stop_model('checktracer: Ox too big in tropo.',255)
         end if
       end do
#ifdef SHINDELL_STRAT_CHEM
       do L=maxL+1,LM
         if(y(n_Ox,L)/y(nM,L) > 1.5d-5) then
           write(6,*)'Ox @ I,J,L,Ox,Itime:',I,J,L,y(n_Ox,L),Itime
           call stop_model('checktracer: Ox too big in strato.',255)
         end if
       end do
#endif
       END IF
       
c general check on maximum of tracers:
      IF(checkmax) THEN
      do L=1,LM
       do igas=1,nlast
        if(y(igas,L)/y(nM,L) > tlimit(igas)) then
          write(6,*) trname(igas),'@ I,J,L,Y :',I,J,L,y(igas,L)
          call stop_model('checktracer: tracer upper limit',255)
        end if
       end do
      end do
      END IF

c check for negative tracers:
      IF(checkNeg) THEN
      do L=1,LM
       do igas=1,nlast
        if(y(igas,L) < 0.d0) THEN
          write(6,*)trname(igas),
     &    'negative @ tau,I,J,L,y:',Itime,I,J,L,y(igas,L)
          call stop_model('checktracer: tracer is negative',255)
        end if
       enddo
      end do
      END IF

c check for unreal (not-a-number) tracers (maybe SGI only?):
      IF(checkNaN) THEN
      do L=1,LM
       do igas=1,nlast
        if(.NOT.(y(igas,L) > 0.d0.OR.y(igas,L) <= 0.d0)) THEN
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
      
      
