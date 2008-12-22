#include "rundeck_opts.h"
      MODULE AMP_AEROSOL
!@sum Driver for Aerosol Microphysics
!@auth Susanne Bauer
      USE TRACER_COM
      USE AERO_CONFIG, ONLY: NMODES
      USE AERO_PARAM,  ONLY: NEMIS_SPCS
      USE MODEL_COM,   ONLY: LM
      IMPLICIT NONE
      SAVE

C**************  Latitude-Dependant (allocatable) *******************
      ! Mie lookup tables
      REAL*8, DIMENSION(17,15,6,23)      :: AMP_EXT, AMP_ASY, AMP_SCA   !(17,15,6,12) (IM,RE,lambda,size)
      REAL*8, DIMENSION(17,15,23)        :: AMP_Q55
      ! 1 Dim arrays for Radiation
      REAL*8, DIMENSION(LM,nmodes)       :: Reff_LEV, NUMB_LEV
      COMPLEX*8, DIMENSION(LM,nmodes,6)  :: RindexAMP
      REAL*8, DIMENSION(LM,nmodes,7)      :: dry_Vf_LEV

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)       :: AQsulfRATE !(i,j,l)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: DIAM       ![m](i,j,l,nmodes)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: MASSH2O    ![ug/m^3](i,j,l,nmodes)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: AMP_dens   !density(i,j,l,nmodes)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: AMP_TR_MM  !molec. mass(i,j,l,nmodes)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: NACTV      != 1.0D-30  ![#/m^3](i,j,l,nmodes)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: VDDEP_AERO != 1.0D-30  ![m/s](i,j,nmodes,2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)         :: DTR_AMPe   !(jm,ntmAMP) ! Emission diagnostic - hardcoded to 10 in TRDIAG
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)       :: DTR_AMP    !(7,jm,ntmAMP)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)       :: DTR_AMPm   !(2,jm,ntmAMP)	  
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: NUMB_SS  ! Sea salt number concentration [#/gb]
#ifdef BLK_2MOM
      REAL*8, ALLOCATABLE, DIMENSION(:,:)         :: NACTC      ! = 1.0D-30  ![#/m3](l,nmodes)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)         :: NAERC      ! = 1.0D-30  ![#/m3](l,nmodes)
#endif


!-------------------------------------------------------------------------------------------------------------------------
!     The array VDDEP_AERO(X,Y,Z,I,1) contains current values for the dry deposition velocities 
!     for aerosol number concentrations for mode I. 
!     The array VDDEP_AERO(X,Y,Z,I,2) contains current values for the dry deposition velocities 
!     for aerosol mass   concentrations for mode I. 
!     Values in VDDEP_AERO are saved in subr. MATRIX at each time step. 
!-------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------
!     The array NACTV(X,Y,Z,I) contains current values of the number of aerosol particles 
!     activated in clouds for each mode I for use outside of the MATRIX microphysical module.
!     Values in NACTV are saved in subr. MATRIX at each time step. 
!-------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------      
!     1 - BC  2-BCmix 3-OC 4-OCmix 5-SS1  6-SS2 7-D1 8-D2
!-------------------------------------------------------------------------------------------------------------------------      
!-------------------------------------------------------------------------------------------------------------------------
!     The array DIAM(x,y,z) contains current values of some measure of average ambient mode diameter for each mode  
!       for use outside of the MATRIX microphysical module where it is calculated. 
!       Values in DIAM are saved at the top of the subr. MATRIX before microphysical evolution 
!       for the current time step is done. 
!
!     The current measure of particle diameter is the diameter of average mass:
! 
!        DIAM(x,y,z) = [ (6/pi) * (Mi/Ni) * (1/D) ]^(1/3)
!
!     with Mi the total mass concentration (including water) in mode i, Ni the number concentration in mode i, and
!     D a constant ambient particle density, currently set to D = DENSP = 1.4 g/cm^3. 
!-------------------------------------------------------------------------------------------------------------------------      
      END MODULE AMP_AEROSOL

      SUBROUTINE MATRIX_DRV
      USE TRACER_COM
      USE TRDIAG_COM, only : taijs=>taijs_loc,ijts_AMPp,ijts_AMPe
     $             ,ijts_AMPm,ijts_AMPext,ijts_AMPpdf,
     $             tajln=>tajln_loc,itcon_AMP,itcon_AMPe,itcon_AMPm
     $             ,itcon_surf
      USE AMP_AEROSOL
      USE AEROSOL_SOURCES, only: off_HNO3

      USE MODEL_COM, only : im,jm,lm     ! dimensions
     $                     ,t            ! potential temperature (C)
     $                     ,q            ! saturatered pressure
     $                     ,dtsrc
      USE GEOM, only: dxyp,imaxj,BYDXYP
      USE CONSTANT,   only:  lhe,mair,gasc   
      USE FLUXES, only: tr3Dsource,trsource,trsrfflx,trflux1
      USE DYNAMICS,   only: pmid,pk,byam,gz, am   ! midpoint pressure in hPa (mb)
!                                           and pk is t mess up factor
!                                           BYAM  1/Air mass (m^2/kg)
      USE AERO_CONFIG
      USE AERO_INIT
      USE AERO_PARAM, only: IXXX, IYYY, ILAY, NEMIS_SPCS
      USE AERO_SETUP
      USE PBLCOM,     only: EGCM !(LM,IM,JM) 3-D turbulent kinetic energy [m^2/s^2]
      USE DOMAIN_DECOMP,only: GRID, GET

      IMPLICIT NONE

      REAL(8):: TK,RH,PRES,TSTEP,AQSO4RATE
      REAL(8):: AERO(NAEROBOX)     ! aerosol conc. [ug/m^3] or [#/m^3]
      REAL(8):: GAS(NGASES)        ! gas-phase conc. [ug/m^3]
      REAL(8):: EMIS_MASS(NEMIS_SPCS) ! mass emission rates [ug/m^3]
      REAL(8):: SPCMASS(NMASS_SPCS+2)
      REAL(8):: DT_AERO(NDIAG_AERO,NAEROBOX) !NDIAG_AERO=15
      REAL(8):: yS, yM, ZHEIGHT1,WUP,AVOL 
      REAL(8) :: PDF1(NBINS)               ! number or mass conc. at each grid point [#/m^3] or [ug/m^3]       
      REAL(8) :: PDF2(NBINS)               ! number or mass conc. at each grid point [#/m^3] or [ug/m^3]       
      INTEGER:: j,l,i,n,J_0, J_1
C**** functions
      REAL(8):: QSAT
c***  daily output
      character*30 diam_1

      CALL GET(grid, J_STRT =J_0, J_STOP =J_1)
#ifndef  TRACERS_SPECIAL_Shindell
      CALL READ_OFFHNO3(OFF_HNO3)
#endif
      DTR_AMP(:,J_0:J_1,:)      = 0.d0
      DTR_AMPm(:,J_0:J_1,:)     = 0.d0
      NACTV(:,J_0:J_1,:,:)      = 0.d0 
      VDDEP_AERO(:,J_0:J_1,:,:) = 0.d0 
      DIAM(:,J_0:J_1,:,:)       = 0.d0
      AMP_dens(:,J_0:J_1,:,:)   = 0.d0
      AMP_TR_MM(:,J_0:J_1,:,:)  = 0.d0
      NUMB_SS(:,J_0:J_1,:,:)    = 0.d0
#ifdef BLK_2MOM
      NACTC(:,:)=0.d0
      NAERC(:,:)=0.d0
#endif




      DO L=1,LM                            
      DO J=J_0,J_1                          
      DO I=1,IM                     

      IXXX = I
      IYYY = J
      ILAY = L
      DT_AERO(:,:) = 0.d0
      EMIS_MASS(:) = 0.d0
      AERO(:)      = 0.d0
! meteo
      TK = pk(l,i,j)*t(i,j,l)           !should be in [K]
      RH = MIN(1.,q(i,j,l)/QSAT(TK,lhe,pmid(l,i,j))) ! rH [0-1]
      PRES= pmid(l,i,j)*100.                  ! pmid in [hPa]
      TSTEP=dtsrc
      ZHEIGHT1 = GZ(i,j,l) /1000./9.81
      WUP = SQRT(.6666667*EGCM(l,i,j))  ! updraft velocity

c avol [m3/gb] mass of air pro m3      
      AVOL = am(l,i,j)*dxyp(j)/mair*1000.d0*gasc*tk/pres 
! in-cloud SO4 production rate [ug/m^3/s] ::: AQsulfRATE [kg] 
      AQSO4RATE = AQsulfRATE (i,j,l)* 1.d9  / AVOL /dtsrc
c conversion trm [kg/gb] -> [ug /m^3]
      GAS(1) = trm(i,j,l,n_H2SO4)* 1.d9 / AVOL! [ug H2SO4/m^3]
c conversion trm [kg/kg] -> [ug /m^3]
      GAS(2) = off_HNO3(i,j,l)*1.d9 * 1.292 !   [ug HNO3/m^3]
c conversion trm [kg/gb] -> [ug /m^3]
      GAS(3) = trm(i,j,l,n_NH3)* 1.d9 / AVOL!   [ug NH3 /m^3]
!  [kg/s] -> [ug/m3/s]

       DO n=1,ntmAMP 
c conversion trm [kg/gb] -> AERO [ug/m3]
         if(AMP_NUMB_MAP(n).eq. 0) then
       AERO(AMP_AERO_MAP(n)) =trm(i,j,l,n)*1.d9 / AVOL ! ug/m3
          else

       AERO(AMP_AERO_MAP(n)) =trm(i,j,l,n)/ AVOL       !  #/m3
          endif
       ENDDO

       if (L.eq.1) then     
!      Emis Mass [ug/m3/s] <-- trflux1[kg/s]
#ifdef TRACERS_AMP_M4
      EMIS_MASS(2) =  MAX(trflux1(i,j,n_M_ACC_SU)*1.d9 / AVOL,0.d0)
      EMIS_MASS(3) =  MAX(trflux1(i,j,n_M_BC1_BC)*1.d9 / AVOL,0.d0)
      EMIS_MASS(4) =  MAX(trflux1(i,j,n_M_OCC_OC)*1.d9 / AVOL,0.d0)
      EMIS_MASS(5) =  MAX(trflux1(i,j,n_M_DD1_DU)*1.d9 / AVOL,0.d0)
      EMIS_MASS(6) =  MAX(trflux1(i,j,n_M_SSS_SS)*1.d9 / AVOL,0.d0)
      EMIS_MASS(10)=  MAX(trflux1(i,j,n_M_DD2_DU)*1.d9 / AVOL,0.d0)
#else
      EMIS_MASS(1) =  MAX(trflux1(i,j,n_M_AKK_SU)*1.d9 / AVOL,0.d0)
      EMIS_MASS(2) =  MAX(trflux1(i,j,n_M_ACC_SU)*1.d9 / AVOL,0.d0)
      EMIS_MASS(3) =  MAX(trflux1(i,j,n_M_BC1_BC)*1.d9 / AVOL,0.d0)
      EMIS_MASS(4) =  MAX(trflux1(i,j,n_M_OCC_OC)*1.d9 / AVOL,0.d0)
      EMIS_MASS(5) =  MAX(trflux1(i,j,n_M_DD1_DU)*1.d9 / AVOL,0.d0)
      EMIS_MASS(6) =  MAX(trflux1(i,j,n_M_SSA_SS)*1.d9 / AVOL,0.d0)
      EMIS_MASS(7) =  MAX(trflux1(i,j,n_M_SSC_SS)*1.d9 / AVOL,0.d0)
      EMIS_MASS(10)=  MAX(trflux1(i,j,n_M_DD2_DU)*1.d9 / AVOL,0.d0)
#endif
        endif
!      Emis Mass [ug/m3/s] <-- trflux1[kg/s]
#ifdef TRACERS_AMP_M4
      EMIS_MASS(2) =  EMIS_MASS(2) + (tr3Dsource(i,j,l,2,n_M_ACC_SU)*1.d9 / AVOL)
      EMIS_MASS(3) =  EMIS_MASS(3) + (tr3Dsource(i,j,l,2,n_M_BC1_BC)*1.d9 / AVOL)
      EMIS_MASS(9) =  EMIS_MASS(9) + (tr3Dsource(i,j,l,2,n_M_OCC_OC)*1.d9 / AVOL)
#else
      EMIS_MASS(1) =  EMIS_MASS(1) + (tr3Dsource(i,j,l,2,n_M_AKK_SU)*1.d9 / AVOL)
      EMIS_MASS(2) =  EMIS_MASS(2) + (tr3Dsource(i,j,l,2,n_M_ACC_SU)*1.d9 / AVOL)
      EMIS_MASS(3) =  EMIS_MASS(3) + (tr3Dsource(i,j,l,2,n_M_BC1_BC)*1.d9 / AVOL)
      EMIS_MASS(8) =  EMIS_MASS(8) + (tr3Dsource(i,j,l,2,n_M_BOC_BC)*1.d9 / AVOL)
      EMIS_MASS(9) =  EMIS_MASS(9) + (tr3Dsource(i,j,l,2,n_M_BOC_OC)*1.d9 / AVOL)
#endif
       CALL SPCMASSES(AERO,GAS,SPCMASS)

       CALL MATRIX(AERO,GAS,EMIS_MASS,TSTEP,TK,RH,PRES,AQSO4RATE,WUP,DT_AERO) 
       CALL SIZE_PDFS(AERO,PDF1,PDF2)
 
       DO n=1, ntmAMP
          if(AMP_NUMB_MAP(n).eq. 0) then
      tr3Dsource(i,j,l,1,n) =((AERO(AMP_AERO_MAP(n)) *AVOL *1.d-9)
     *        -trm(i,j,l,n)) /dtsrc 
          else
      tr3Dsource(i,j,l,1,n) =((AERO(AMP_AERO_MAP(n)) *AVOL)
     *        -trm(i,j,l,n)) /dtsrc
          endif   
       ENDDO

       NUMB_SS(i,j,l,1) = AERO(22) *AVOL
       NUMB_SS(i,j,l,2) = AERO(25) *AVOL ! but has only tiny number in it
    

      tr3Dsource(i,j,l,1,n_H2SO4) =((GAS(1)*AVOL *1.d-9)
     *        -trm(i,j,l,n_H2SO4)) /dtsrc 
      tr3Dsource(i,j,l,1,n_NH3)   =((GAS(3)*AVOL *1.d-9)
     *        -trm(i,j,l,n_NH3)) /dtsrc

#ifdef  TRACERS_SPECIAL_Shindell
      tr3Dsource(i,j,l,3,n_HNO3)  =((GAS(2)/1.292 * 1.d-9)
     *        -trm(i,j,l,n_HNO3))/dtsrc
#endif
c       DT_AERO(:,:) = DT_AERO(:,:) * dtsrc !DT_AERO [# or ug/m3/s] , taijs [kg m2/kg(air)], byam [kg/m2]

c Update physical properties per mode
       do n=1,ntmAMP
c Diagnostic of Processes - Sources and Sincs - timestep included
          if(AMP_NUMB_MAP(n).eq. 0) then  !taijs [kg/s] -> in acc [kg/m2*s]
        taijs(i,j,ijts_AMPp(1,n)) =taijs(i,j,ijts_AMPp(1,n))+(DT_AERO(9,AMP_AERO_MAP(n))* AVOL / 1.d9)
      DTR_AMP(1,j,n)=DTR_AMP(1,j,n)+(DT_AERO(9,AMP_AERO_MAP(n))* AVOL /1.d9)
        taijs(i,j,ijts_AMPp(2,n)) =taijs(i,j,ijts_AMPp(2,n))+(DT_AERO(10,AMP_AERO_MAP(n))* AVOL/ 1.d9)
      DTR_AMP(2,j,n)=DTR_AMP(2,j,n)+(DT_AERO(10,AMP_AERO_MAP(n))* AVOL /1.d9)
        taijs(i,j,ijts_AMPp(3,n)) =taijs(i,j,ijts_AMPp(3,n))+(DT_AERO(11,AMP_AERO_MAP(n))* AVOL/ 1.d9)
      DTR_AMP(3,j,n)=DTR_AMP(3,j,n)+(DT_AERO(11,AMP_AERO_MAP(n))* AVOL /1.d9)
        taijs(i,j,ijts_AMPp(4,n)) =taijs(i,j,ijts_AMPp(4,n))+(DT_AERO(12,AMP_AERO_MAP(n))* AVOL/ 1.d9)
      DTR_AMP(4,j,n)=DTR_AMP(4,j,n)+(DT_AERO(12,AMP_AERO_MAP(n))* AVOL /1.d9)
        taijs(i,j,ijts_AMPp(5,n)) =taijs(i,j,ijts_AMPp(5,n))+(DT_AERO(13,AMP_AERO_MAP(n))* AVOL/ 1.d9)
      DTR_AMP(5,j,n)=DTR_AMP(5,j,n)+(DT_AERO(13,AMP_AERO_MAP(n))* AVOL /1.d9)
        taijs(i,j,ijts_AMPp(6,n)) =taijs(i,j,ijts_AMPp(6,n))+(DT_AERO(14,AMP_AERO_MAP(n))* AVOL/ 1.d9)
      DTR_AMP(6,j,n)=DTR_AMP(6,j,n)+(DT_AERO(14,AMP_AERO_MAP(n))* AVOL /1.d9)
        taijs(i,j,ijts_AMPp(7,n)) =taijs(i,j,ijts_AMPp(7,n))+(DT_AERO(15,AMP_AERO_MAP(n))* AVOL/ 1.d9)
      DTR_AMP(7,j,n)=DTR_AMP(7,j,n)+(DT_AERO(15,AMP_AERO_MAP(n))* AVOL /1.d9)
          else
        taijs(i,j,ijts_AMPp(1,n)) =taijs(i,j,ijts_AMPp(1,n))+(DT_AERO(2,AMP_AERO_MAP(n))* AVOL)
      DTR_AMP(1,j,n)=DTR_AMP(1,j,n)+(DT_AERO(2,AMP_AERO_MAP(n)) * AVOL)
        taijs(i,j,ijts_AMPp(2,n)) =taijs(i,j,ijts_AMPp(2,n))+(DT_AERO(3,AMP_AERO_MAP(n))* AVOL)
      DTR_AMP(2,j,n)=DTR_AMP(2,j,n)+(DT_AERO(3,AMP_AERO_MAP(n)) * AVOL)
       taijs(i,j,ijts_AMPp(3,n)) =taijs(i,j,ijts_AMPp(3,n))+(DT_AERO(1,AMP_AERO_MAP(n))* AVOL)
      DTR_AMP(3,j,n)=DTR_AMP(3,j,n)+(DT_AERO(1,AMP_AERO_MAP(n)) * AVOL)
        taijs(i,j,ijts_AMPp(4,n)) =taijs(i,j,ijts_AMPp(4,n))+(DT_AERO(4,AMP_AERO_MAP(n))* AVOL)
      DTR_AMP(4,j,n)=DTR_AMP(4,j,n)+(DT_AERO(4,AMP_AERO_MAP(n)) * AVOL)
        taijs(i,j,ijts_AMPp(5,n)) =taijs(i,j,ijts_AMPp(5,n))+(DT_AERO(5,AMP_AERO_MAP(n))* AVOL)
      DTR_AMP(5,j,n)=DTR_AMP(5,j,n)+(DT_AERO(5,AMP_AERO_MAP(n)) * AVOL)
        taijs(i,j,ijts_AMPp(6,n)) =taijs(i,j,ijts_AMPp(6,n))+(DT_AERO(6,AMP_AERO_MAP(n))* AVOL)
      DTR_AMP(6,j,n)=DTR_AMP(6,j,n)+(DT_AERO(6,AMP_AERO_MAP(n)) * AVOL)
        taijs(i,j,ijts_AMPp(7,n)) =taijs(i,j,ijts_AMPp(7,n))+(DT_AERO(7,AMP_AERO_MAP(n))* AVOL)
      DTR_AMP(7,j,n)=DTR_AMP(7,j,n)+(DT_AERO(7,AMP_AERO_MAP(n)) * AVOL)
          endif
       select case (trname(n)) !taijs [kg * m2/kg air] -> in acc [kg/kg air]
      CASE('N_AKK_1 ','N_ACC_1 ','N_DD1_1 ','N_DS1_1 ','N_DD2_1 ','N_DS2_1 ','N_OCC_1 ','N_BC1_1 ',
     *     'N_BC2_1 ','N_BC3_1 ','N_DBC_1 ','N_BOC_1 ','N_BCS_1 ','N_MXX_1 ','N_OCS_1 ')
c - 3d acc output
        taijs(i,j,ijts_AMPm(l,1,n))=taijs(i,j,ijts_AMPm(l,1,n)) + DIAM(i,j,l,AMP_MODES_MAP(n))
        taijs(i,j,ijts_AMPm(l,2,n))=taijs(i,j,ijts_AMPm(l,2,n)) + (NACTV(i,j,l,AMP_MODES_MAP(n))*AVOL*byam(l,i,j))
C** Surabi's addition to save as /cm-3 or just mult. std o/p by 1.2*1.e-16 to get units in cm-3
c       taijs(i,j,ijts_AMPm(l,2,n))=taijs(i,j,ijts_AMPm(l,2,n)) + (NACTV(i,j,l,AMP_MODES_MAP(n))*1.e-6)
c       if(NACTV(i,j,l,AMP_MODES_MAP(n)) .gt. 10.)
c    *  write(6,*)"Aerosolconv",AVOL*byam(l,i,j),NACTV(i,j,l,AMP_MODES_MAP(n)),
c    *l,i,j,n,dxyp(j),mair,gasc,tk,pres,am(l,i,j)

c - 2d PRT Diagnostic
         DTR_AMPm(1,j,n)=DTR_AMPm(1,j,n)+(DIAM(i,j,l,AMP_MODES_MAP(n))*1d6*1d18)
         DTR_AMPm(2,j,n)=DTR_AMPm(2,j,n)+NACTV(i,j,l,AMP_MODES_MAP(n))*AVOL
       end select

      enddo !n
c - special diag: Size distribution pdfs
       if (l.eq.1) taijs(i,j,ijts_AMPpdf(l,:))=taijs(i,j,ijts_AMPpdf(l,:)) + (PDF1(:)*AVOL*byam(l,i,j))
c - N_SSA, N_SSC, M_SSA_SU
        taijs(i,j,ijts_AMPext(l,1))=taijs(i,j,ijts_AMPext(l,1)) + (NACTV(i,j,l,SEAS_MODE_MAP(1))*AVOL*byam(l,i,j))
        taijs(i,j,ijts_AMPext(l,2))=taijs(i,j,ijts_AMPext(l,2)) + (NACTV(i,j,l,SEAS_MODE_MAP(2))*AVOL*byam(l,i,j))
        taijs(i,j,ijts_AMPext(l,3))=taijs(i,j,ijts_AMPext(l,3)) +  DIAM(i,j,l,SEAS_MODE_MAP(1)) 
        taijs(i,j,ijts_AMPext(l,4))=taijs(i,j,ijts_AMPext(l,4)) +  DIAM(i,j,l,SEAS_MODE_MAP(2)) 
        taijs(i,j,ijts_AMPext(l,5))=taijs(i,j,ijts_AMPext(l,5)) + (AERO(22) *AVOL *byam(l,i,j)) 
        taijs(i,j,ijts_AMPext(l,6))=taijs(i,j,ijts_AMPext(l,6)) + (AERO(25) *AVOL *byam(l,i,j)) 
#ifdef BLK_2MOM
      NACTC(l,:) = NACTV(i,j,l,:)
c     NAERC(l,:) = AERO(AMP_AERO_MAP(n)) !trm(i,j,l,:)/AVOL
      NAERC(l,:) = trm(i,j,l,:)/AVOL
c     do n=1, ntmAMP
c     if (NAERC(l,n).gt.0.0001) write(6,*)"ACTIV",NAERC(l,n),
c    *AERO(AMP_AERO_MAP(n)),
c    *trm(i,j,l,n),AVOL,i,j,l,n,NACTC(l,n)
c     enddo
#endif
      ENDDO !i
      ENDDO !j
      ENDDO !l

      do n=1,ntmAMP 

      CALL DIAGTCB(DTR_AMP(1,:,n),itcon_amp(1,n),n)
      CALL DIAGTCB(DTR_AMP(2,:,n),itcon_amp(2,n),n)
      CALL DIAGTCB(DTR_AMP(3,:,n),itcon_amp(3,n),n)
      CALL DIAGTCB(DTR_AMP(4,:,n),itcon_amp(4,n),n)
      CALL DIAGTCB(DTR_AMP(5,:,n),itcon_amp(5,n),n)
      CALL DIAGTCB(DTR_AMP(6,:,n),itcon_amp(6,n),n)
      CALL DIAGTCB(DTR_AMP(7,:,n),itcon_amp(7,n),n)

       select case (trname(n))       
      CASE('N_AKK_1 ','N_ACC_1 ','N_DD1_1 ','N_DS1_1 ','N_DD2_1 ',
     *     'N_DS2_1 ','N_OCC_1 ','N_BC1_1 ',
     *     'N_BC2_1 ','N_BC3_1 ','N_DBC_1 ','N_BOC_1 ','N_BCS_1 ','N_MXX_1 ','N_OCS_1 ')
      CALL DIAGTCB(DTR_AMPm(1,:,n),itcon_AMPm(1,n),n)
      CALL DIAGTCB(DTR_AMPm(2,:,n),itcon_AMPm(2,n),n)
      CASE('M_DD1_DU ','M_DD2_DU ','M_DDD_DU','M_SSA_SS ','M_SSC_SS ','M_SSS_SS')
      CALL DIAGTCB(DTR_AMPe(:,n),itcon_surf(1,n),n)
       end select
      enddo

c      diam_1='diam_1'
c      call INST_ncoutd(diam_1,DIAM(:,:,:,1))

      RETURN
      END SUBROUTINE MATRIX_DRV
c -----------------------------------------------------------------

c -----------------------------------------------------------------
      SUBROUTINE AMPtrdens(i,j,l,n)
!----------------------------------------------------------------------------------------------------------------------
!     Routine to calculate the actual density per mode
!----------------------------------------------------------------------------------------------------------------------
      USE TRACER_COM
      USE AMP_AEROSOL, only : AMP_dens
      USE AERO_CONFIG, ONLY: NMODES

      IMPLICIT NONE
      Integer :: i,j,l,n

 
         if(AMP_MODES_MAP(n).gt.0)
     &   AMP_dens(i,j,l,AMP_MODES_MAP(n)) = 
     &   sum(trpdens(AMP_trm_nm1(n):AMP_trm_nm2(n)) * trm(i,j,l,AMP_trm_nm1(n):AMP_trm_nm2(n))) 
     & / sum(trm(i,j,l,AMP_trm_nm1(n):AMP_trm_nm2(n)))
         if (AMP_dens(i,j,l,AMP_MODES_MAP(n)).le.0) AMP_dens(i,j,l,AMP_MODES_MAP(n)) = trpdens(AMP_MODES_MAP(n))

      RETURN
      END SUBROUTINE AMPtrdens
c -----------------------------------------------------------------
c -----------------------------------------------------------------
      SUBROUTINE AMPtrmass(i,j,l,n)
!----------------------------------------------------------------------------------------------------------------------
!     Routine to calculate the actual molecular mass per mode
!----------------------------------------------------------------------------------------------------------------------
      USE TRACER_COM
      USE AMP_AEROSOL, only : AMP_TR_MM
      USE AERO_CONFIG, ONLY: NMODES

      IMPLICIT NONE
      Integer :: i,j,l,n

         if(AMP_MODES_MAP(n).gt.0)
     &   AMP_TR_MM(i,j,l,AMP_MODES_MAP(n)) = 
     &   sum(tr_mm(AMP_trm_nm1(n):AMP_trm_nm2(n)) * trm(i,j,l,AMP_trm_nm1(n):AMP_trm_nm2(n))) 
     & / sum(trm(i,j,l,AMP_trm_nm1(n):AMP_trm_nm2(n)))
         if (AMP_TR_MM(i,j,l,AMP_MODES_MAP(n)).le.0) AMP_TR_MM(i,j,l,AMP_MODES_MAP(n)) = tr_mm(AMP_MODES_MAP(n))


      RETURN
      END SUBROUTINE AMPtrmass
c -----------------------------------------------------------------
      SUBROUTINE SPCMASSES(AERO,GAS,SPCMASS)
!----------------------------------------------------------------------------------------------------------------------
!     Routine to calculate the total mass concentration of each model species:
!     SULF, BCAR, OCAR, DUST, SEAS, NO3, NH4. Aerosol water is not treated. 
!----------------------------------------------------------------------------------------------------------------------
      USE AERO_SETUP, ONLY: SULF_MAP, BCAR_MAP, OCAR_MAP, DUST_MAP, SEAS_MAP
      USE AERO_PARAM
      USE AERO_CONFIG
      IMPLICIT NONE
      REAL(8) :: AERO(NAEROBOX)
      REAL(8) :: GAS(NGASES)    
      REAL(8) :: SPCMASS(NMASS_SPCS+2)
      SPCMASS(1) = SUM( AERO( SULF_MAP(:) ) ) + GAS( GAS_H2SO4 )
      SPCMASS(2) = SUM( AERO( BCAR_MAP(:) ) )
      SPCMASS(3) = SUM( AERO( OCAR_MAP(:) ) )
      SPCMASS(4) = SUM( AERO( DUST_MAP(:) ) )
      SPCMASS(5) = SUM( AERO( SEAS_MAP(:) ) )
      SPCMASS(6) = AERO( MASS_NO3 )           + GAS( GAS_HNO3 )
      SPCMASS(7) = AERO( MASS_NH4 )           + GAS( GAS_NH3  )

 
      RETURN
      END SUBROUTINE SPCMASSES
     
      SUBROUTINE SIZE_PDFS(AERO,PDF1,PDF2)
      USE AERO_PARAM, ONLY: PI6, DENSP, IXXX, IYYY, ILAY
      USE AERO_CONFIG, ONLY: NMODES, NAEROBOX,NBINS
      USE AERO_SETUP, ONLY: SIG0, CONV_DPAM_TO_DGN, NUMB_MAP, MODE_NAME
      USE AMP_AEROSOL, only: DIAM
      IMPLICIT NONE

      ! Arguments.
       REAL(8), INTENT(IN) :: AERO(NAEROBOX)! aerosol conc. [ug/m^3] or [#/m^3]

      ! Local variables. 

      INTEGER :: I, N
!      INTEGER, PARAMETER :: NBINS = 30! 200    ! number of bins [1]     defined in config 
      REAL(8) :: DGRID(NBINS)              ! fixed diameter        grid [um]
      REAL(8) :: MGRID(NBINS)              ! fixed mass/particle   grid [ug/particle]
      REAL(8) :: DLOWER(NBINS)             ! lower boundary fixed diameter grid [um]
      REAL(8) :: DUPPER(NBINS)             ! upper boundary fixed diameter grid [um]
      REAL(8) :: NTOT(NMODES)              ! number concentration for each mode [#/m^3]
      REAL(8) :: PDF(NBINS,2,NMODES)       ! number or mass conc. at each grid point [#/m^3] or [ug/m^3]       
      REAL(8) :: PDF1(NBINS)               ! number or mass conc. at each grid point [#/m^3] or [ug/m^3]       
      REAL(8) :: PDF2(NBINS)               ! number or mass conc. at each grid point [#/m^3] or [ug/m^3]       
      REAL(8) :: DNDLOGD(NMODES)           ! dN/dlog10(Dp) [ #/m^3]
      REAL(8) :: DMDLOGD(NMODES)           ! dM/dlog10(Dp) [ug/m^3]
      REAL(8) :: RDMIN                     ! reciprocal of DMIN to optimize coagulation [1/um]
      REAL(8) :: RDLOGDSC                  ! reciprocal of log10 of the grid spacing [1]
      REAL(8) :: SCALE, F, SUM1, SUM2      ! scratch variables 
      REAL(8) :: DMINL, DMAXL, DG          ! diameters [um]  
      REAL(8) :: FLN                       ! function for lognormal distribution [1]  
      REAL(8), PARAMETER :: DMIN =  0.001D+00   ! smallest particle diameter of the discrete grid [um]
      REAL(8), PARAMETER :: DMAX = 20.000D+00   ! largest  particle diameter of the discrete grid [um]


        DMAXL = DMAX
        DMINL = DMIN

      SCALE    = ( DMAXL / DMINL )**(1.0D+00/REAL(NBINS-1))
      RDLOGDSC = 1.0D+00 / LOG10( SCALE )
      RDMIN    = 1.0D+00 / DMINL
      DO I=1, NBINS
        DGRID(I)  = DMINL * SCALE**(I-1)                  ! [um]
        DLOWER(I) = DGRID(I) / SCALE**0.5D+00             ! [um]
        DUPPER(I) = DGRID(I) * SCALE**0.5D+00             ! [um]
        MGRID(I)  = 1.0D-06 * DENSP * PI6 * DGRID(I)**3   ! [ug/particle]
        DO N=1, NMODES
          DG = 1.0D+06 * DIAM(IXXX,IYYY,ILAY,N) * CONV_DPAM_TO_DGN(N)   ! convert [m] to [um] and Dbar to Dg
          NTOT(N) = AERO( NUMB_MAP(N) )
          F = NTOT(N) * FLN( DGRID(I), DG, SIG0(N) )
          PDF(I,1,N) = F * ( DUPPER(I) - DLOWER(I) )
          PDF(I,2,N) = PDF(I,1,N) * MGRID(I)
          DNDLOGD(N) = PDF(I,1,N) * RDLOGDSC * 1.0D-06           ! convert from [#/m^3] to [#/cm^3]
          DNDLOGD(N) = MAX( DNDLOGD(N), 1.0D-30 )
          DMDLOGD(N) = PDF(I,2,N) * RDLOGDSC                     ! [ug/m^3]
          DMDLOGD(N) = MAX( DMDLOGD(N), 1.0D-30 )
        ENDDO
c        WRITE(IUNIT,91) I, DGRID(I), DNDLOGD(:)
c        WRITE(JUNIT,91) I, DGRID(I), DMDLOGD(:)
      ENDDO

        PDF1(:) = 0.0D+00
        PDF2(:) = 0.0D+00
      DO N=1, NMODES
        DO I=1, NBINS
          PDF1(I) = PDF1(I)  + PDF(I,1,N)
          PDF2(I) = PDF2(I)  + PDF(I,2,N)
          SUM1 = SUM1 + PDF(I,1,N)
          SUM2 = SUM2 + PDF(I,2,N)
        ENDDO
      ENDDO

      RETURN 
      END SUBROUTINE SIZE_PDFS


      REAL(8) FUNCTION FLN(X,XG,SIGMAG)
      REAL(8) :: X      ! particle radius or diameter [any units]
      REAL(8) :: XG     ! geometric mean radius or diameter [any units]
      REAL(8) :: SIGMAG ! geometric standard deviation [monodisperse = 1.0]
      REAL(8), PARAMETER :: SQRTTWOPI = 2.506628275D+00
      FLN = EXP(-0.5D+00*(LOG(X/XG)/LOG(SIGMAG))**2) / (X*LOG(SIGMAG)*SQRTTWOPI)
      RETURN
      END FUNCTION FLN
      subroutine alloc_tracer_amp_com(grid)
!@SUM  To alllocate arrays whose sizes now need to be determined
!@+    at run-time
!@auth Susanne Bauer
!@ver  1.0
      use domain_decomp, only : dist_grid, get
      use model_com, only     : im,lm
      use tracer_com, only    : ntmAMP
      use amp_aerosol
      use aero_config, only   : nmodes

      IMPLICIT NONE

      type (dist_grid), intent(in) :: grid
      integer :: ier, J_1H, J_0H
      logical :: init = .false.

      if(init)return
      init=.true.
    
      call get( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
 

! I,J,L
      allocate(  AQsulfRATE(IM,J_0H:J_1H,LM)   )
! other dimensions
      allocate(  DTR_AMPe(J_0H:J_1H,ntmAMP)    )
      allocate(  DTR_AMPm(2,J_0H:J_1H,ntmAMP)  )
      allocate(  DTR_AMP(7,J_0H:J_1H,ntmAMP)   )
      allocate(  DIAM(IM,J_0H:J_1H,LM,nmodes)  )
      allocate(  MASSH2O(IM,J_0H:J_1H,LM,nmodes)  )
      allocate(  AMP_TR_MM(IM,J_0H:J_1H,LM,nmodes)  )
      allocate(  AMP_dens(IM,J_0H:J_1H,LM,nmodes)  )
      allocate(  NACTV(IM,J_0H:J_1H,LM,nmodes) )
      allocate(  VDDEP_AERO(IM,J_0H:J_1H,nmodes,2))
      allocate(  NUMB_SS(IM,J_0H:J_1H,LM,2))

      NACTV   = 1.0D-30
      DIAM    = 1.0D-30
      MASSH2O = 1.0D-30
             
#ifdef BLK_2MOM
      NACTV      = 1.0D-30
      allocate(  NACTC(LM,nmodes) )
      NACTC      = 1.0D-30
      allocate(  NAERC(LM,nmodes) )
      NAERC      = 1.0D-30
#endif    
      return
      end subroutine alloc_tracer_amp_com
      
