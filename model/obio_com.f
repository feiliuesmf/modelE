      MODULE obio_com
!@sum  obio_com contains the parameters, arrays and definitions
!@+    necessary for the OceanBiology routines
!@auth NR
!@ver  1.0e-11

      USE obio_dim

      implicit none

      SAVE


#include "dimensions.h"
#include "dimension2.h"


!this part is taken out of common_blocks.h in hycom
c --- dobio       activate Watson Gregg's ocean biology code
      logical dobio
      data dobio/.true./

      !these are defined in the rundeck for the coupled runs
!     character*60 flnmoas,flnmsolz,flnmbio1,flnmbio2,flnmbio3

!     data flnmoas/
!    .   '/g6/aromanou/NewGrid/20w/oasim181x180_20w'/
!     data flnmsolz/
!    .   '/g6/aromanou/NewGrid/20w/solz_forbio_20w.dat'/
!     data flnmbio1/
!    .   '/g6/aromanou/NewGrid/20w/nitoa181x180_20w.asc'/
!     data flnmbio2/
!    .   '/g6/aromanou/NewGrid/20w/siloa181x180_20w.asc'/
!     data flnmbio3/
!    .   '/g6/aromanou/NewGrid/20w/glodap181x180_20w.asc'/
c




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

      real rmu4     !growth on ammonium, NH4
      real rmu3     !growth on nitrate, NO4
      real rmu5     !growth on silica
      real rmuf     !growth on iron
      real zoo      !herbivores (zooplankton)
      real dphy     !death rate of phytoplankton
      real co2mon(26,12)        !26 years 1979-2004, 12 months
      real viscfac
      common /bptend1/ rmu4(nchl),rmu3(nchl),rmu5(nchl),rmuf(nchl)
      common /bptend2/ zoo(ntyp),dphy(ntyp),viscfac(kdm)
!$OMP THREADPRIVATE(/bptend1/)
!$OMP THREADPRIVATE(/bptend2/)

      integer npst,npnd   !starting and ending array index for PAR
      real WtoQ           !Watts/m2 to quanta/m2/s conversion
      data npst,npnd /3,17/
      common /bwq/ WtoQ(nlt)

! reduced rank arrays for obio_model calculations
      integer ihra_ij
      real temp1d,dp1d,obio_P,det,car,avgq1d,gcmax1d,atmFe_ij,covice_ij
     .    ,saln1d,p1d,alk1d
      common /reducarr1/temp1d(kdm),dp1d(kdm),obio_P(kdm,ntyp+n_inert)
     .                 ,det(kdm,ndet),car(kdm,ncar),avgq1d(kdm)
     .                 ,gcmax1d(kdm),saln1d(kdm),p1d(kdm+1)
     .                 ,alk1d(kdm)
      common /reducarr2/ihra_ij
      common /reducarr3/atmFe_ij,covice_ij
!$OMP THREADPRIVATE(/reducarr1/)
!$OMP THREADPRIVATE(/reducarr2/)
!$OMP THREADPRIVATE(/reducarr3/)


      integer inwst,inwnd,jnwst,jnwnd     !starting and ending indices 
                                          !for daylight 
                                          !in i and j directions
      real    acdom

      common /bnwi/  inwst,inwnd,jnwst,jnwnd
      common /bcdom/ acdom(kdm,nlt)       !absorptio coefficient of CDOM


      real P_tend,rmuplsr,D_tend
      real obio_ws,bn,tfac,pnoice,pnoice2
      real wsdet,rikd
      real Fescav
      real tzoo
      common /bkpt/ P_tend(kdm,ntyp+n_inert)!bio tendency (dP/dt)
      common /bgro/ rmuplsr(kdm,nchl)   !growth+resp rate
      common /bkws/ obio_ws(kdm,nchl)   !phyto sinking rate
      common /bzoo/ tzoo                !herbivore T-dependence
      common /bkwsdet2/ wsdet(kdm+1,ndet) !detrital sinking rate
      common /btfac/ tfac(kdm)          !phyto T-dependence
      common /bpnoice/ pnoice,pnoice2   !pct ice-free
      common /bbn/ bn(kdm)              !N:chl ratio
      common /bikd/ rikd(kdm,nchl)      !photoadaption state

!$OMP THREADPRIVATE(/bkpt/)
!$OMP THREADPRIVATE(/bkws/)
!$OMP THREADPRIVATE(/bkwsdet2/)
!$OMP THREADPRIVATE(/bbn/)
!$OMP THREADPRIVATE(/bpnoice/)
!$OMP THREADPRIVATE(/btfac/)
!$OMP THREADPRIVATE(/bikd/)
!$OMP THREADPRIVATE(/bgro/)
!$OMP THREADPRIVATE(/bzoo/)

C if NCHL_DEFINED > 3
      real wshc,gcmax

      common /bwshc/ wshc(kdm)             !cocco sinking rate
      common /bgcmax/ gcmax(idm,jdm,kdm)   !cocco max growth rate
      common /bscav/ Fescav(kdm)           !iron scavenging rate
!$OMP THREADPRIVATE(/bscav/)
!$OMP THREADPRIVATE(/bwshc/)

C endif

      common /bdtend/ D_tend(kdm,ndet)     !detrtial tendency
!$OMP THREADPRIVATE(/bdtend/)


      real :: C_tend
      real :: pCO2, pCO2_ij
      common /bctend/ C_tend(kdm,ncar)      !carbon tendency
      common /bco2/ pCO2(idm,jdm)           !partial pressure of CO2
      common /bco2ij/ pCO2_ij               !partial pressure of CO2
!$OMP THREADPRIVATE(/bctend/)
!$OMP THREADPRIVATE(/bco2ij/)

      real :: gro
      common /bgro2D/ gro(kdm,nchl)         !realized growth rate
!$OMP THREADPRIVATE(/bgro2D/)

      integer :: day_of_month, hour_of_day

      real :: rhs
      common /brhs/ rhs(kdm,14,16)
!$OMP THREADPRIVATE(/brhs/)

      END MODULE obio_com
