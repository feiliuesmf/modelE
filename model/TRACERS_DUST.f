#include "rundeck_opts.h"
      SUBROUTINE dust_wet(i,j)
!@sum  Simple scheme for wet deposition of dust/mineral tracers
!@auth Ina Tegen, Reha Cakmur, Jan Perlwitz

#ifndef TRACERS_WATER
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      USE constant,ONLY : Grav
      USE resolution,ONLY : Jm,Lm
      USE model_com,ONLY : zatmo
      USE fluxes,ONLY : prec
      USE clouds,ONLY : tm_dust,tmom_dust,trprc_dust
      USE tracer_com,ONLY : Ntm_dust,trname
      USE tracers_dust,ONLY : prelay
      USE dynamics,ONLY : gz

      IMPLICIT NONE

      REAL*8,PARAMETER :: Z=700.

      INTEGER,INTENT(IN) :: i,j

      INTEGER :: l,layer,n
      INTEGER,DIMENSION(Jm) :: lwdep
      REAL*8,DIMENSION(Jm) :: h
      REAL*8 :: y
      REAL*8 :: height
      REAL*8,DIMENSION(Lm,Ntm_dust) :: tmold

      COMMON /dustprv/ l,layer,y,height,lwdep,h,tmold
!$OMP  THREADPRIVATE (/dustprv/)

#ifdef WET_DEPO_Ina
      SELECT CASE(j)
      CASE (1:6)
        lwdep(j)=3
        h(j)=2800
        lwdep(jm+1-j)=3
        h(jm+1-j)=2800
      CASE (7:12)
        lwdep(j)=4
        h(j)=4900
        lwdep(jm+1-j)=4
        h(jm+1-j)=4900
      CASE (13:16)
        LWDEP(J)=5
        h(j)=7400
        lwdep(jm+1-j)=5
        h(jm+1-j)=7400
      CASE (17:23)
        lwdep(j)=6
        h(j)=10300
        lwdep(jm+1-j)=6
        h(jm+1-j)=10300
      END SELECT
      
      y = Z*(prec(i,j)/h(j))
      IF (y > 1.) y=1.

#else

      layer=0
      DO l=Lm,1,-1
        IF (prelay(i,j,l) /= 0.) THEN 
          layer=l
          EXIT
        END IF
      END DO

      IF (layer == 1) THEN
        height=(gz(i,j,layer)-zatmo(i,j))/Grav
      ELSE IF (layer /= 0) THEN
        height=(gz(i,j,layer)-gz(i,j,1))/Grav
      END IF
      IF (layer /= 0) THEN
        y=Z*(prec(i,j)/height)
        IF (y > 1.) y=1.
      ELSE
        y=0.D0
      END IF

#endif

c**** Wet Deposition

      tmold=tm_dust
      DO n=1,Ntm_dust
        DO l=1,layer
          tm_dust(l,n)=tm_dust(l,n)*(1-y)
          tmom_dust(:,l,n)=tmom_dust(:,l,n)*(1-y)
        END DO
      END DO
      trprc_dust=tmold-tm_dust

#endif
#endif

      RETURN
      END SUBROUTINE dust_wet
