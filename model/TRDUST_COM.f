#include "rundeck_opts.h"
      MODULE tracers_dust
!@sum  dust tracer parameters and variables
!@auth Reha Cakmur, Jan Perlwitz, Ina Tegen

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      USE constant,ONLY : By6
      USE resolution,ONLY : Im,Jm,Lm
      USE model_com,ONLY : JMperY,JDperY
      USE tracer_com,ONLY : Ntm_dust

      IMPLICIT NONE

!@param nDustTurbij index of dust dry turbulent deposition in ijts_source
!@param nDustWetij index of dust wet deposition in ijts_source
      INTEGER,PARAMETER :: nDustTurbij=2,
     &                     nDustWetij=4
!@param nDustTurbjl index of dust dry turbulent deposition in jls_source
!@param nDustWetjl index of dust wet deposition in jls_3Dsource
      INTEGER,PARAMETER :: nDustTurbjl=2,
     &                     nDustWet3Djl=2
!@param uplfac uplift factor for each size class of soil dust [kg*s**2/m**5]
#ifdef TRACERS_DUST_CUB_SAH
      REAL*8,PARAMETER :: Upclsi=52.D-9
#else !default case
      REAL*8,PARAMETER :: Upclsi=2.2D-9
#endif
#ifdef TRACERS_DUST_CUB_SAH
!@param By8 0.25d0/2d0
      REAL*8,PARAMETER :: By8=0.25D0/2D0
#else !default case
      REAL*8,PARAMETER :: By4=1D0/4D0
#endif
!@param fracn fraction of uplifted soil for each size class of dust [1]
#ifdef TRACERS_DUST_CUB_SAH
      REAL*8 :: Fracl=By6,Frasi=By8
c      REAL*8 :: Fracn(Ntm_dust)=(/By6,By8,By8,By8/)
#else !default case
      REAL*8 :: Fracl=0.18D0,Frasi=By4
c      REAL*8 :: Fracn(Ntm_dust)=(/0.18D0,By4,By4,By4/)
#endif
!@var hbaij  accumulated precipitation - evaporation balance
      REAL*8 :: hbaij(im,jm),ricntd(im,jm)
!@var dryhr  number of hours with evaporation-precipitation greater Zero
!@var dryhr  to allow dust emission
!@var frclay fraction of clay
!@var frsilt fraction of silt
!@var vtrsh  threshold wind speed above which dust emission is allowed
      REAL*4 :: dryhr(im,jm),frclay(im,jm),frsilt(im,jm),vtrsh(im,jm)
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
!@param Mtrac number of different fields with tracer fractions in grid box
!@param Mtrac 5 clay; 5 silt
      INTEGER,PARAMETER :: Mtrac=10
!@var minfr distribution of tracer fractions in grid box
      REAL*4 :: minfr(Im,Jm,Mtrac)
#endif
!@var qdust flag whether conditions for dust emission are fulfilled
      LOGICAL :: qdust(Im,Jm)
      REAL*4 :: ers_data(im,jm,JMperY),gin_data(im,jm)
      INTEGER,PARAMETER :: lim=294,ljm=244,lkm=9
      REAL*8 :: curint(Im,Jm),table(lim,ljm,lkm),x1(lim),x2(ljm),x3(lkm)
      REAL*8 :: wsubtke_com(Im,Jm),wsubwd_com(Im,Jm),wsubwm_com(Im,Jm)
      REAL*8 :: prelay(im,jm,Lm)
!@var d_dust Prescribed daily dust emissions
      REAL*8 :: d_dust(Im,Jm,Ntm_dust,JDperY)
#endif

      END MODULE tracers_dust
