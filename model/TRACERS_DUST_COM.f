#include "rundeck_opts.h"
      MODULE tracers_dust_com
!@sum  dust tracer parameters and variables
!@auth Reha Cakmur, Jan Perlwitz, Ina Tegen

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      USE constant,ONLY : By6
      USE resolution,ONLY : Im,Jm,Lm
      USE tracer_com,ONLY : Ntm_dust

!@param nDustEmij index of dust emission in ijts_source
!@param nDustTurbij index of dust dry turbulent deposition in ijts_source
!@param nDustGravij index of dust gravitational settling in ijts_source
!@param nDustWetij index of dust wet deposition in ijts_source
      INTEGER,PARAMETER :: nDustEmij=1,
     &                     nDustTurbij=2,
     &                     nDustGravij=3,
     &                     nDustWetij=4
!@param nDustEmjl index of dust emission in jls_source
!@param nDustTurbjl index of dust dry turbulent deposition in jls_source
!@param nDustGrav3Djl index of dust gravitational settling in jls_3Dsource
!@param nDustWetjl index of dust wet deposition in jls_3Dsource
      INTEGER,PARAMETER :: nDustEmjl=1,
     &                     nDustTurbjl=2,
     &                     nDustGrav3Djl=1,
     &                     nDustWet3Djl=2
!@param Z scavenging ratio for wet deposition of dust      
      REAL*8,PARAMETER :: Z=700.
!@param Zld assumed thickness of layers [m]
      REAL*8,PARAMETER :: Zld(12)=(/400.,800.,1400.,2100.,2600.,2200.,
     .     2000.,2100.,2300.,2700.,3300.,4300./)
!@param dradius soil dust particle radius [m]
#ifdef TRACERS_DUST
      REAL*8,PARAMETER :: Dradius(Ntm_dust)=(/0.75D-06,2.2D-06,4.4D-06,
     .     6.7D-06/) !n=1: clay, n=2,3,4: silt
#else
#ifdef TRACERS_MINERALS
      REAL*8,PARAMETER :: Dradius(Ntm_dust)=(/0.75D-06,0.75D-06,
     &     0.75D-06,0.75D-06,0.75D-06,2.2D-06,2.2D-06,2.2D-06,2.2D-06,
     &     2.2D-06,4.4D-06,4.4D-06,4.4D-06,4.4D-06,4.4D-06,6.7D-06,
     &     6.7D-06,6.7D-06,6.7D-06,6.7D-06/)
#endif
#endif
!@param uplfac uplift factor for each size class of soil dust [kg*s**2/m**5]
#ifdef TRACERS_DUST
      REAL*8,PARAMETER :: Uplfac(Ntm_dust)=(/52.D-9,52.D-9,52.D-9,
     &     52.D-9/)
#else
#ifdef TRACERS_MINERALS
      REAL*8,PARAMETER :: Uplfac(Ntm_dust)=(/52.D-9,52.D-9,52.D-9,
     &     52.D-9,52.D-9,52.D-9,52.D-9,52.D-9,52.D-9,52.D-9,52.D-9,
     &     52.D-9,52.D-9,52.D-9,52.D-9,52.D-9,52.D-9,52.D-9,52.D-9,
     &     52.D-9/)
#endif
#endif
!@param By8 0.25d0/2d0
      REAL*8,PARAMETER :: By8=0.25D0/2D0
!@param fracn fraction of uplifted soil for each size class of dust [1]
#ifdef TRACERS_DUST
      REAL*8 :: Fracn(Ntm_dust)=(/By6,By8,By8,By8/)
#else
#ifdef TRACERS_MINERALS
      REAL*8 :: Fracn(Ntm_dust)=(/By6,By6,By6,By6,By6,By8,By8,By8,By8,
     &     By8,By8,By8,By8,By8,By8,By8,By8,By8,By8,By8/)
#endif
#endif
!@var hbaij  accumulated precipitation - evaporation balance
      REAL*8 :: hbaij(im,jm),ricntd(im,jm),wsubtke(Im,Jm),wsubwd(Im,Jm),
     &     wsubwm(Im,Jm)
!@var dryhr  number of hours with evaporation-precipitation greater Zero
!@var dryhr  to allow dust emission
!@var frclay fraction of clay
!@var frsilt fraction of silt
!@var vtrsh  threshold wind speed above which dust emission is allowed
      REAL*4 :: dryhr(im,jm),frclay(im,jm),frsilt(im,jm),vtrsh(im,jm)
#endif

      END MODULE tracers_dust_com
