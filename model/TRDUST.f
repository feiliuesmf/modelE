#include "rundeck_opts.h"
      SUBROUTINE dust_emission_constraints(i,j,itype,ptype,wsm)
!@sum  local constrainsts for dust tracer emission valid for all dust bins
!@auth Jan Perlwitz, Reha Cakmur, Ina Tegen

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      USE model_com,ONLY : dtsrc,nisurf,jmon,wfcs
      USE tracer_com,ONLY : imDUST
      USE fluxes,ONLY : prec,pprec,pevap
      USE ghycom,ONLY : snowe,wearth,aiearth
      USE tracers_dust,ONLY : curint,dryhr,ers_data,hbaij,lim,ljm,qdust,
     &     ricntd,table,vtrsh,x1,x2,x3,wsubtke_com,wsubwd_com,wsubwm_com
      USE pbl_drv,ONLY : wsubtke,wsubwd,wsubwm

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: i,j,itype
      REAL*8,INTENT(IN) :: ptype,wsm

      REAL*8 :: hbaijold,hbaijd
      REAL*8 :: soilwet
      LOGICAL :: pmei
      REAL*8 :: sigma,ans,dy
      REAL*8 :: soilvtrsh,workij1,workij2

      IF (imDUST == 0) THEN

#ifndef DUST_EMISSION_EXTERN

#ifdef TRACERS_DUST_CUB_SAH
c     Checking whether accumulated precipitation - evaporation
c     less/equal than Zero for a succeeding number of hours greater/equal
c     than threshold dryhr to permit dust emission
      hbaijold=hbaij(i,j)
      hbaij(i,j)=hbaijold+pprec(i,j)*ptype/nisurf-pevap(i,j,itype)
      hbaijd=hbaij(i,j)-hbaijold
      IF (itype == 4 .AND. hbaijd <= 0) THEN
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

      IF (pmei .AND. snowe(i,j) <= 1 .AND. vtrsh(i,j) > 0. .AND.
     &     wsm > vtrsh(i,j)) THEN
        qdust(i,j)=.TRUE.
      ELSE
        qdust(i,j)=.FALSE.
      END IF
#else !default case
      IF (itype == 4 .AND. snowe(i,j) <= 1
     &     .AND. ers_data(i,j,jmon) < -13. .AND.
     &     ers_data(i,j,jmon) /= -99.) THEN
        qdust(i,j)=.TRUE.
      ELSE
        qdust(i,j)=.FALSE.
      END IF

      IF (qdust(i,j)) THEN

        soilwet=(WEARTH(I,J)+AIEARTH(I,J))/(WFCS(I,J)+1.D-20)
        if (soilwet.gt.1.) soilwet=1.d0
        soilvtrsh=8.d0*(exp(0.25d0*soilwet))

        curint(i,j)=0.d0
        workij1=0.d0
        workij2=0.d0

c     There is no moist convection, sigma is composed of TKE and DRY
c     convective velocity scale

        if (wsubwm == 0.) then
          sigma=wsubtke+wsubwd
          if (sigma > 0.1 .OR. wsm > 1.) then
            call ratint2(x1,x2,x3,table,lim,ljm,wsm,sigma,soilvtrsh,ans,
     &           dy)
            curint(i,j)=exp(ans)
          endif
        endif

c     When there is moist convection, the sigma is the combination of
c     all three subgrid scale parameters (i.e. independent or dependent)
c     Takes into account that the moist convective velocity scale acts
c     only over 5% of the area.

        if (wsubwm /= 0.) then
          sigma=wsubtke+wsubwd+wsubwm
          if (sigma > 0.1 .OR. wsm > 1.) then
            call ratint2(x1,x2,x3,table,lim,ljm,wsm,sigma,soilvtrsh,ans,
     &           dy)
            workij1=exp(ans)*0.05 !!!0.05 for the MC area
          endif

          sigma=wsubtke+wsubwd
          if (sigma > 0.1 .OR. wsm > 1.) then
            call ratint2(x1,x2,x3,table,lim,ljm,wsm,sigma,soilvtrsh,ans,
     &           dy)
            workij2=exp(ans)*0.95 !!!0.95 for the rest
          endif
          curint(i,j)=workij1+workij2
        endif

        if (sigma == 0.) then
          if (wsm > soilvtrsh) then
            curint(i,j)=(wsm-soilvtrsh)*wsm**2
          else
            curint(i,j)=0D0
          endif
        endif

      END IF
#endif
#else

      IF (itype == 4) THEN
        wsubtke_com(i,j)=wsubtke
        wsubwd_com(i,j)=wsubwd
        wsubwm_com(i,j)=wsubwm
      END IF
      
#endif

      ELSE IF (imDUST == 1) THEN
        IF (itype == 4) THEN
          qdust(i,j)=.TRUE.
        ELSE
          qdust(i,j)=.FALSE.
        END IF
      END IF

#endif

      RETURN
      END SUBROUTINE dust_emission_constraints

      SUBROUTINE local_dust_emission(i,j,n,wsm,ptype,dsrcflx)
!@sum  selects routine for calculating local dust source flux
!@auth Jan Perlwitz, Reha Cakmur, Ina Tegen

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      USE constant,ONLY : sday
      USE model_com,ONLY : jday
      USE geom,ONLY : dxyp
      USE tracer_com,ONLY : trname,imDUST,n_clayilli
      USE tracers_dust,ONLY : curint,fracn,frclay,frsilt,gin_data,qdust,
     &     uplfac,vtrsh,d_dust
#ifdef TRACERS_MINERALS
     &     ,minfr
#endif

#ifndef DUST_EMISSION_EXTERN
      IMPLICIT NONE

      INTEGER,INTENT(IN) :: i,j,n
      REAL*8,INTENT(IN) :: ptype,wsm
      REAL*8,INTENT(OUT) :: dsrcflx

      INTEGER :: n1
      REAL*8 :: frtrac

      IF (imDUST == 0) THEN
c     Interactive dust emission

      IF (.NOT. qdust(i,j)) THEN
        dsrcflx=0D0
      ELSE
#ifdef TRACERS_DUST
#ifdef TRACERS_DUST_CUB_SAH
        SELECT CASE(trname(n))
        CASE ('Clay')
          frtrac=frclay(i,j)
        CASE ('Silt1','Silt2','Silt3')
          frtrac=frsilt(i,j)
        END SELECT

        dsrcflx=Uplfac(n)*frtrac*Fracn(n)*(wsm-vtrsh(i,j))*wsm**2

#else ! default case

        dsrcflx=Uplfac(n)*Fracn(n)*gin_data(i,j)*curint(i,j)
#endif
#else
#ifdef TRACERS_MINERALS

        SELECT CASE(trname(n))
        CASE ('ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar',
     &        'Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps')
          n1=n-n_clayilli+1
        CASE ('Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps')
          n1=n-n_clayilli-4
        CASE ('Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps')
          n1=n-n_clayilli-9
        END SELECT

#ifdef TRACERS_DUST_CUB_SAH
        SELECT CASE(trname(n))
        CASE ('ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar')
          frtrac=frclay(i,j)
        CASE ('Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps',
     &        'Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps',
     &        'Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps')
          frtrac=frsilt(i,j)
        END SELECT

        dsrcflx=minfr(i,j,n1)*Uplfac(n)*frtrac*Fracn(n)*
     &       (wsm-vtrsh(i,j))*wsm**2

#else ! default case

        dsrcflx=minfr(i,j,n1)*Uplfac(n)*Fracn(n)*gin_data(i,j)*
     &       curint(i,j)

#endif
#endif
#endif

      END IF

      ELSE IF (imDUST == 1) THEN
c     prescribed AEROCOM dust emission

      IF (.NOT. qdust(i,j)) THEN
        dsrcflx=0D0
      ELSE

        dsrcflx=d_dust(i,j,n,jday)/Sday/dxyp(j)/ptype

#ifdef TRACERS_MINERALS
      
        SELECT CASE(trname(n))
        CASE ('ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar',
     &        'Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps')
          n1=n-n_clayilli+1
        CASE ('Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps')
          n1=n-n_clayilli-4
        CASE ('Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps')
          n1=n-n_clayilli-9
        END SELECT

        dsrcflx=dsrcflx*minfr(i,j,n1)

#endif

      END IF

      END IF

#endif
#endif
            
      RETURN
      END SUBROUTINE local_dust_emission

      SUBROUTINE ratint2(x1a,x2a,x3a,ya,m,n,x1,x2,x3,y,dy)

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      implicit none
      INTEGER, INTENT(IN) :: m,n
      INTEGER NMAX,MMAX
      REAL*8, INTENT(IN) :: x1,x2,x3,x1a(m),x2a(n),x3a(9),ya(m,n,9)
      REAL*8, INTENT(OUT):: y,dy
      PARAMETER (NMAX=4,MMAX=4)
      INTEGER i,j,k,jjj,iii,xx,yy,zz,kkk
      REAL*8 ymtmp(MMAX),yntmp(nmax),x11(4),x22(4),x33(4)
      real*8 yotmp(MMAX)

      call locate(x1a,m,x1,xx)
      call locate(x2a,n,x2,yy)
      call locate(x3a,9,x3,zz)

      do i=1,4
         if (zz.eq.1) then
            kkk=i
         else
            kkk=i-1
         endif
         x33(i)=x3a(zz+kkk-1)
         do k=1,4
            iii=k-1
            x22(k)=x2a(yy+iii-1)
            do j=1,4
               jjj=j-1
               x11(j)=x1a(xx+jjj-1)
               yntmp(j)=ya(xx+jjj-1,yy+iii-1,zz+kkk-1)
c               print *,yntmp(j)
            enddo
            if (yntmp(1).eq.-1000.or.yntmp(2).eq.-1000.) then
               ymtmp(k)=-1000.
            else
               call ratint(x11,yntmp,4,x1,ymtmp(k),dy)
            endif
         enddo
         if (ymtmp(1).eq.-1000.or.ymtmp(2).eq.-1000.) then
            yotmp(i)=-1000.
         else
            call ratint(x22,ymtmp,4,x2,yotmp(i),dy)
         endif
c       print *,yotmp(i)
      enddo
      if (yotmp(1).eq.-1000.or.yotmp(2).eq.-1000.)  then
         y=-1000.
      else
         if (x3.ge.2.*x1) then
            call polint(x33,yotmp,4,x3,y,dy)
         else
            call ratint(x33,yotmp,4,x3,y,dy)
         endif
      endif
#endif

      return
      END SUBROUTINE RATINT2
C  (C) Copr. 1986-92 Numerical Recipes Software 'W3.

      SUBROUTINE ratint(xa,ya,n,x,y,dy)

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      IMPLICIT NONE

      INTEGER NMAX
      INTEGER, INTENT(IN) :: n
      REAL*8, INTENT(IN) :: x,xa(n),ya(n)
      real*8 TINY
      REAL*8, INTENT(OUT) :: y,dy
      PARAMETER (NMAX=4,TINY=1.d-25)
      INTEGER i,m,ns
      REAL*8 dd,h,hh,t,w,c(NMAX),d(NMAX)
      ns=1
      hh=abs(x-xa(1))
      do 11 i=1,n
        h=abs(x-xa(i))
        if (h.eq.0.d0)then
          y=ya(i)
          dy=0.0d0
          return
        else if (h.lt.hh) then
          ns=i
          hh=h
        endif
        c(i)=ya(i)
        d(i)=ya(i)+TINY
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          w=c(i+1)-d(i)
          h=xa(i+m)-x
          t=(xa(i)-x)*d(i)/h
          dd=t-c(i+1)
          if(dd.eq.0.d0) CALL stop_model('failure in ratint',255)
          dd=w/dd
          d(i)=c(i+1)*dd
          c(i)=t*dd
 12     continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
        ns=ns-1
        endif
        y=y+dy
13    continue
#endif

      return
      END SUBROUTINE RATINT

      SUBROUTINE polint(xa,ya,n,x,y,dy)

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      INTEGER NMAX
      REAL*8, INTENT(OUT) :: y,dy
      REAL*8, INTENT(IN) :: x,xa(n),ya(n)
      PARAMETER (NMAX=300)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.d0) CALL stop_model('failure in polint',255)
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
#endif

      return
      END SUBROUTINE POLINT



      SUBROUTINE locate(xx,n,x,j)

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      implicit none
      INTEGER, INTENT(IN):: n
      INTEGER, INTENT(OUT):: j
      REAL*8, INTENT(IN):: x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
#endif

      return
      END SUBROUTINE LOCATE
