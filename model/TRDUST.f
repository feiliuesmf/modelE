#include "rundeck_opts.h"
      MODULE tracers_dust
!@sum  dust tracer parameters and variables
!@auth Reha Cakmur, Jan Perlwitz, Ina Tegen

#ifdef TRACERS_DUST
      USE constant,ONLY : By6
      USE resolution,ONLY : Im,Jm
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
      REAL*8,PARAMETER :: Uplfac(Ntm_dust)=(/52.D-9,52.D-9,52.D-9,
     &     52.D-9/)
!@param By8 0.25d0/2d0
      REAL*8,PARAMETER :: By8=0.25D0/2D0
!@param fracn fraction of uplifted soil for each size class of dust [1]
      REAL*8 :: Fracn(Ntm_dust)=(/By6,By8,By8,By8/)
!@var hbaij  accumulated precipitation - evaporation balance
      REAL*8 :: hbaij(im,jm),ricntd(im,jm)
!@var dryhr  number of hours with evaporation-precipitation greater Zero
!@var dryhr  to allow dust emission
!@var frclay fraction of clay
!@var frsilt fraction of silt
!@var vtrsh  threshold wind speed above which dust emission is allowed
      REAL*4 :: dryhr(im,jm),frclay(im,jm),frsilt(im,jm),vtrsh(im,jm)
!@var qdust  flag whether conditions for dust emission are fulfilled
      LOGICAL :: qdust(Im,Jm)
#endif

      CONTAINS

      SUBROUTINE dust_emission_constraints(i,j)

#ifdef TRACERS_DUST
      USE resolution,ONLY : Jm
      USE model_com,ONLY : dtsrc,fearth,nisurf
      USE fluxes,ONLY : prec,evapor
      USE geom,ONLY : imaxj
      USE ghycom,ONLY : snowe !earth snow amount
      USE pblcom,ONLY : wsavg

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: i,j

      REAL*8 :: hbaijold,hbaijd
      LOGICAL :: pmei

c     Checking whether accumulated precipitation - evaporation
c     less/equal than Zero for a succeeding number of hours greater/equal
c     than threshold dryhr to permit dust emission
      hbaij(i,j)=hbaijold+prec(i,j)*fearth(i,j)-evapor(i,j,4)
      hbaijd=hbaij(i,j)-hbaijold
      IF (hbaijd <= 0) THEN
        ricntd(i,j)=ricntd(i,j)+Dtsrc/3600./nisurf
        IF (ricntd(i,j) >= dryhr(i,j) .AND. dryhr(i,j) /= 0) THEN
          pmei=.TRUE.
        ELSE
          pmei=.FALSE.
        END IF
      ELSE
        ricntd(i,j)=0.
        pmei=.FALSE.
      END IF

      IF (pmei) THEN
        IF (fearth(i,j) > 0. .AND. snowe(i,j) <= 1 .AND.
     &       vtrsh(i,j) > 0. .AND. wsavg(i,j) > vtrsh(i,j)) THEN
          qdust(i,j)=.TRUE.
        ELSE
          qdust(i,j)=.FALSE.
        END IF
      END IF
#endif

      RETURN
      END SUBROUTINE dust_emission_constraints

      SUBROUTINE local_dust_emission(i,j,n)
!@sum  dust source flux
!@auth Jan Perlwitz, Reha Cakmur, Ina Tegen

#ifdef TRACERS_DUST
      USE model_com,ONLY : fearth
      USE tracer_com,ONLY : trname
      USE fluxes,ONLY : dustflux
      USE geom,ONLY : dxyp
      USE pblcom,ONLY : wsavg

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: i,j,n

      REAL*8 :: frtrac

#ifdef TRACERS_DUST_MINERAL8
#ifdef TRACERS_DUST_TURB
      CALL loc_dustflux_turb_min8
#else
      CALL loc_dustflux_cub_min8
#endif
#else
#ifdef TRACERS_DUST_TURB
      CALL loc_dustflux_turb_sah
#else ! default case
      SELECT CASE(trname(n))
      CASE ('Clay')
        frtrac=frclay(i,j)
      CASE ('Silt1','Silt2','Silt3')
        frtrac=frsilt(i,j)
      END SELECT
      CALL loc_dustflux_cub_sah(n,qdust(i,j),dxyp(j),wsavg(i,j),
     &     vtrsh(i,j),frtrac,fearth(i,j),dustflux(i,j,n))
#endif
#endif
#endif
            
      RETURN
      END SUBROUTINE local_dust_emission

      SUBROUTINE loc_dustflux_cub_sah(n,qdustij,dxypj,wsavgij,vtrshij,
     &     frtracij,fearthij,dustfluxij)
!@sum  local dust source flux physics according to Ina's old cubic scheme
!@auth Ina Tegen, Jan Perlwitz, Reha Cakmur

#ifdef TRACERS_DUST

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: n
      REAL*4,INTENT(IN) :: vtrshij
c     vtrshij  wind speed threshold in grid cell i,j [m/s]
      REAL*8,INTENT(IN) :: dxypj,fearthij,frtracij,wsavgij
c     dxypj    area of grid cell at j [m**2]
c     fearthij  fraction of land area in grid cell i,j [1]
c     frtrac fraction of dust tracer n in grid cell i,j [1]
c     wsavgij  wind speed at surface in grid cell i,j [m/s]
      LOGICAL,INTENT(IN) :: qdustij
c     qdustij local flag whether conditions for dust emission are fulfilled
      REAL*8,INTENT(OUT) :: dustfluxij
c     dustfluxij  local source flux of dust tracer n [kg/s]

      IF (.NOT. qdustij) THEN
        dustfluxij=0
      ELSE
        dustfluxij=Uplfac(n)*frtracij*Fracn(n)*dxypj*fearthij*
     &       (wsavgij-vtrshij)*wsavgij**2
      END IF
#endif

      RETURN
      END SUBROUTINE loc_dustflux_cub_sah

      SUBROUTINE loc_dustflux_cub_min8
!@sum  local dust source flux physics with Ina's cubic scheme and 8 minerals
!@auth Jan Perlwitz, ...

      RETURN
      END SUBROUTINE loc_dustflux_cub_min8

      SUBROUTINE loc_dustflux_turb_sah
!@sum  local dust source flux physics with turbulent fluxes and Sahara dust
!@auth Reha Cakmur, ...

      RETURN
      END SUBROUTINE loc_dustflux_turb_sah

      SUBROUTINE loc_dustflux_turb_min8
!@sum  local dust source flux physics with turbulent fluxes and 8 minerals
!@auth Jan Perlwitz, Reha Cakmur

      RETURN
      END SUBROUTINE loc_dustflux_turb_min8

      END MODULE tracers_dust
