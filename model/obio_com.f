#include "rundeck_opts.h"
      MODULE obio_com
!@sum  obio_com contains the parameters, arrays and definitions
!@+    necessary for the OceanBiology routines
!@auth NR
!@ver  1.0e-11

      USE obio_dim

      USE hycom_dim_glob
      implicit none

!this part is taken out of common_blocks.h in hycom
c --- dobio       activate Watson Gregg's ocean biology code
      logical dobio
      data dobio/.true./
c

      real, ALLOCATABLE, DIMENSION(:,:)    :: tzoo2d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: tfac3d,bn3d,wshc3d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: Fescav3d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: obio_wsd2d,obio_wsh2d
      real, ALLOCATABLE, DIMENSION(:,:,:,:):: rmuplsr3d,rikd3d
      real, ALLOCATABLE, DIMENSION(:,:,:,:):: acdom3d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: gcmax         !cocco max growth rate
      real, ALLOCATABLE, DIMENSION(:,:)    :: pCO2          !partial pressure of CO2
      real, ALLOCATABLE, DIMENSION(:,:)    :: pCO2_glob
      real, ALLOCATABLE, DIMENSION(:,:)    :: pp2tot_day    !net pp total per day
      real, ALLOCATABLE, DIMENSION(:,:)    :: tot_chlo      !tot chlorophyl at surf. layer

      real, parameter :: obio_deltath = 1.0  !time step in hours
      !!real, parameter :: obio_deltat = obio_deltath*3600.0 !time step in seconds
      real, parameter :: obio_deltat = obio_deltath    !time step in hrs 
                                                       !because all rates are in hrs
 

      integer, parameter :: EUZ_DEFINED=1


      real, parameter :: rlamz=1.0,greff=0.25 !other zoopl. parameters
      real, parameter :: drate=0.05/24.0      !phytoplankton death rate/hr
      real, parameter :: dratez1=0.1/24.0     !zooplankton death rate/hr
      real, parameter :: dratez2=0.5/24.0     !zooplankton death rate/hr
      real, parameter :: regen=0.25           !regeneration fraction

      real co2mon(26,12)        !26 years 1979-2004, 12 months

      integer npst,npnd   !starting and ending array index for PAR
      data npst,npnd /3,17/

      real WtoQ           !Watts/m2 to quanta/m2/s conversion
      common /bwq/ WtoQ(nlt)

! reduced rank arrays for obio_model calculations
      integer ihra_ij
      common /reducarr2/ihra_ij
!$OMP THREADPRIVATE(/reducarr2/)

      real temp1d,dp1d,obio_P,det,car,avgq1d,gcmax1d
     .    ,saln1d,p1d,alk1d
      common /reducarr1/temp1d(kdm),dp1d(kdm),obio_P(kdm,ntyp+n_inert)
     .                 ,det(kdm,ndet),car(kdm,ncar),avgq1d(kdm)
     .                 ,gcmax1d(kdm),saln1d(kdm),p1d(kdm+1)
     .                 ,alk1d(kdm)
!$OMP THREADPRIVATE(/reducarr1/)

      real atmFe_ij,covice_ij
      common /reducarr3/atmFe_ij,covice_ij
!$OMP THREADPRIVATE(/reducarr3/)


      integer inwst,inwnd,jnwst,jnwnd     !starting and ending indices 
                                          !for daylight 
                                          !in i and j directions
      common /bnwi/  inwst,inwnd,jnwst,jnwnd

      real    acdom
      common /bcdom/ acdom(kdm,nlt)       !absorptio coefficient of CDOM
!$OMP THREADPRIVATE(/bcdom/)


      real P_tend                             !bio tendency (dP/dt)
      common /bkpt/ P_tend(kdm,ntyp+n_inert)
!$OMP THREADPRIVATE(/bkpt/)

      real rmuplsr                            !growth+resp rate
      common /bgro/ rmuplsr(kdm,nchl)   
!$OMP THREADPRIVATE(/bgro/)

      real D_tend                             !detrtial tendency
      common /bdtend/ D_tend(kdm,ndet)     
!$OMP THREADPRIVATE(/bdtend/)

      real obio_ws                            !phyto sinking rate
      common /bkws/ obio_ws(kdm+1,nchl)   
!$OMP THREADPRIVATE(/bkws/)

      real bn                                 !N:chl ratio
      common /bbn/ bn(kdm)              
!$OMP THREADPRIVATE(/bbn/)

      real tfac                               !phyto T-dependence
      common /btfac/ tfac(kdm)          
!$OMP THREADPRIVATE(/btfac/)

      real pnoice                    !pct ice-free
      common /bpnoice/ pnoice  
!$OMP THREADPRIVATE(/bpnoice/)

      real wsdet                              !detrital sinking rate
      common /bkwsdet2/ wsdet(kdm+1,ndet) 
!$OMP THREADPRIVATE(/bkwsdet2/)

      real rikd                               !photoadaption state
      common /bikd/ rikd(kdm,nchl)      
!$OMP THREADPRIVATE(/bikd/)

      real tzoo                               !herbivore T-dependence
      common /bzoo/ tzoo                
!$OMP THREADPRIVATE(/bzoo/)

      real Fescav
      common /bscav/ Fescav(kdm)           !iron scavenging rate
!$OMP THREADPRIVATE(/bscav/)


C if NCHL_DEFINED > 3
      real wshc                            !cocco sinking rate
      common /bwshc/ wshc(kdm)             
!$OMP THREADPRIVATE(/bwshc/)
C endif


      real :: C_tend                        !carbon tendency
      common /bctend/ C_tend(kdm,ncar)      
!$OMP THREADPRIVATE(/bctend/)

      real :: pCO2_ij                 !partial pressure of CO2
      common /bco2ij/ pCO2_ij               
!$OMP THREADPRIVATE(/bco2ij/)

      real :: gro                           !realized growth rate
      common /bgro2D/ gro(kdm,nchl)         
!$OMP THREADPRIVATE(/bgro2D/)

      integer :: day_of_month, hour_of_day

      real :: rhs
      common /brhs/ rhs(kdm,14,16)
!$OMP THREADPRIVATE(/brhs/)

      real :: pp2_1d          
      common /bpp2/ pp2_1d(kdm,nchl)           !net primary production
!$OMP THREADPRIVATE(/bpp2/)



      contains

      subroutine alloc_obio_com

      USE obio_dim
      USE hycom_dim_glob
      USE hycom_dim, only : j_0h,j_1h

      ALLOCATE(tzoo2d(idm,j_0h:j_1h))
      ALLOCATE(tfac3d(idm,j_0h:j_1h,kdm),bn3d(idm,j_0h:j_1h,kdm))
      ALLOCATE(wshc3d(idm,j_0h:j_1h,kdm),Fescav3d(idm,j_0h:j_1h,kdm))
      ALLOCATE(obio_wsd2d(idm,j_0h:j_1h,nchl),
     &     obio_wsh2d(idm,j_0h:j_1h,nchl))
      ALLOCATE(rmuplsr3d(idm,j_0h:j_1h,kdm,nchl),
     &     rikd3d(idm,j_0h:j_1h,kdm,nchl))
      ALLOCATE(acdom3d(idm,j_0h:j_1h,kdm,nlt))
      ALLOCATE(gcmax(idm,j_0h:j_1h,kdm))
      ALLOCATE(pCO2(idm,j_0h:j_1h))           
      ALLOCATE(pCO2_glob(idm,jdm))           
      ALLOCATE(pp2tot_day(idm,j_0h:j_1h))
      ALLOCATE(tot_chlo(idm,j_0h:j_1h))

      end subroutine alloc_obio_com

      subroutine gather_pCO2
      USE HYCOM_DIM, only : ogrid
      USE DOMAIN_DECOMP, ONLY: PACK_DATA
 
      call pack_data( ogrid, pCO2, pCO2_glob )

      end subroutine gather_pCO2


      END MODULE obio_com
