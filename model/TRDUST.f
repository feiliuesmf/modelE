#include "rundeck_opts.h"
      SUBROUTINE dust_emission_constraints(i,j,itype,ptype,wsgcm,
     &     dsteve1,dsteve2,soilvtrsh,
     &     wsubtke,wsubwd,wsubwm)
!@sum  local constrainsts for dust tracer emission valid for all dust bins
!@auth Jan Perlwitz, Reha Cakmur, Ina Tegen

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      USE model_com,ONLY : dtsrc,nisurf,wfcs
      USE tracer_com,ONLY : imDUST
      USE fluxes,ONLY : prec,pprec,pevap
      USE ghy_com,ONLY : snowe,wearth,aiearth
      USE tracers_dust,ONLY : pdfint,dryhr,hbaij,lim,ljm,lkm,qdust,
     &     ricntd,table,vtrsh,x1,x2,x3,wsubtke_com,wsubwd_com,wsubwm_com

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: i,j,itype
      REAL*8,INTENT(IN) :: ptype,wsgcm
      REAL*8,INTENT(OUT) :: dsteve1,dsteve2,soilvtrsh
      REAL*8,INTENT(IN) :: wsubtke,wsubwd,wsubwm

      REAL*8 :: mcfrac=0.05
      REAL*8 :: hbaijold,hbaijd
      REAL*8 :: soilwet
      LOGICAL :: pmei
      REAL*8 :: sigma,ans,dy
      REAL*8 :: workij1,workij2,wsgcm1

      IF (imDUST == 0) THEN

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
     &     wsgcm > vtrsh(i,j)) THEN
        qdust(i,j)=.TRUE.
      ELSE
        qdust(i,j)=.FALSE.
      END IF
#else !default case
      IF (itype == 4 .AND. snowe(i,j) <= 1) THEN
        qdust(i,j)=.TRUE.
      ELSE
        qdust(i,j)=.FALSE.
        dsteve1=0.D0
        dsteve2=0.D0
        soilvtrsh=0.D0
      END IF

      IF (qdust(i,j)) THEN

        soilwet=(WEARTH(I,J)+AIEARTH(I,J))/(WFCS(I,J)+1.D-20)
        if (soilwet.gt.1.) soilwet=1.d0
        soilvtrsh=8.d0*(exp(0.7d0*soilwet))

        pdfint(i,j)=0.d0
        workij1=0.d0
        workij2=0.d0
        wsgcm1=wsgcm

c     There is no moist convection, sigma is composed of TKE and DRY
c     convective velocity scale
        IF (wsubwm == 0.) THEN
          sigma=wsubtke+wsubwd
c     No need to calculate the emission below these values since
c     the emission is zero
          IF (sigma > 0.1 .OR. wsgcm1 > 1.) THEN
c     This is the case when sigma is very small and we approximate
c     the function by standard dust emission formula
            IF (sigma < 0.0005 .AND. wsgcm1 > 1.) THEN
              IF (wsgcm1 > soilvtrsh) THEN
                pdfint(i,j)=(wsgcm1-soilvtrsh)*wsgcm1**2.D0
              ELSE
                pdfint(i,j)=0.d0
              END IF
c     This is the case when wsgcm1 is very small and we set it
c     equal to one of the smallest values in the table index
            ELSE IF (sigma > 0.1 .AND. wsgcm1 < 0.0005) THEN
              wsgcm1=0.0005d0
c     Linear Polynomial fit (Default)
              CALL polint3dlin(x1,x2,x3,table,lim,ljm,lkm,wsgcm1,sigma,
     &             soilvtrsh,ans,dy)
c     Cubic Polynomial fit (Not Used, Optional)
c              CALL polint3dlicub(x1,x2,x3,table,lim,ljm,lkm,wsgcm1, 
c     &             sigma,soilvtrsh,ans,dy)
              pdfint(i,j)=exp(ans)
            ELSE
c     Linear Polynomial fit (Default)
              CALL polint3dlin(x1,x2,x3,table,lim,ljm,lkm,wsgcm1,sigma,
     &             soilvtrsh,ans,dy)
c     Cubic Polynomial fit (Not Used, Optional) 
c              CALL polint3dlicub(x1,x2,x3,table,lim,ljm,lkm,wsgcm1,
c     &             sigma,soilvtrsh,ans,dy)
              pdfint(i,j)=exp(ans)
            END IF
          END IF

        ELSE

c     When there is moist convection, the sigma is the combination of
c     all three subgrid scale parameters (i.e. independent or dependent)
c     Takes into account that the moist convective velocity scale acts
c     only over 5% of the area.

          sigma=wsubtke+wsubwd+wsubwm
c     No need to calculate the emission below these values since
c     the emission is Zero
          IF (sigma > 0.1 .OR. wsgcm1 > 1.) THEN
c     This is the case when sigma is very small and we approximate
c     the function by standard dust emission formula
            IF (sigma < 0.0005 .AND. wsgcm1 > 1.) THEN
              IF (wsgcm1 > soilvtrsh) THEN
                workij1=mcfrac*(wsgcm1-soilvtrsh)*wsgcm1**2.D0
              ELSE
                workij1=0.d0
              END IF
c     This is the case when wsgcm1 is very small and we set it
c     equal to one of the smallest values in the table index
            ELSE IF (sigma > 0.1 .AND. wsgcm1 < 0.0005) THEN
              wsgcm1=0.0005d0
c     Linear Polynomial fit (Default)
              CALL polint3dlin(x1,x2,x3,table,lim,ljm,lkm,wsgcm1,sigma,
     &             soilvtrsh,ans,dy)
c     Cubic Polynomial fit (Not Used, Optional)
c               call polint3dlicub(x1,x2,x3,table,lim,ljm,lkm,wsgcm1,
c     &              sigma,soilvtrsh,ans,dy)
              workij1=mcfrac*exp(ans)
            ELSE
c     Linear Polynomial fit (Default)
              CALL polint3dlin(x1,x2,x3,table,lim,ljm,lkm,wsgcm1,sigma,
     &             soilvtrsh,ans,dy)
c     Cubic Polynomial fit (Not Used, Optional)
c               CALL polint3dlicub(x1,x2,x3,table,lim,ljm,lkm,wsgcm1,
c     &              sigma,soilvtrsh,ans,dy)
              workij1=mcfrac*exp(ans)
            END IF
          END IF

          sigma=wsubtke+wsubwd
c     No need to calculate the emission below these values since
c     the emission is Zero
          IF (sigma > 0.1 .OR. wsgcm1 > 1.) THEN
c     This is the case when sigma is very small and we approximate
c     the function by standard dust emission formula
            IF (sigma < 0.0005 .AND. wsgcm1 > 1.) THEN
              IF (wsgcm1 > soilvtrsh) THEN
                workij2=(1.d0-mcfrac)*(wsgcm1-soilvtrsh)*wsgcm1**2.D0
              ELSE
                workij2=0.d0
              END IF
c     This is the case when wsgcm1 is very small and we set it
c     equal to one of the smallest values in the table index
            ELSE IF (sigma > 0.1 .AND. wsgcm1 < 0.0005) THEN
              wsgcm1=0.0005d0
c     Linear Polynomial fit (Default)
              CALL polint3dlin(x1,x2,x3,table,lim,ljm,lkm,wsgcm1,sigma,
     &             soilvtrsh,ans,dy)
c     Cubic Polynomial fit (Not Used, Optional)
c               CALL polint3dlicub(x1,x2,x3,table,lim,ljm,lkm,wsgcm1,
c     &              sigma,soilvtrsh,ans,dy)
              workij2=(1.d0-mcfrac)*exp(ans)
            ELSE
c     Linear Polynomial fit (Default)
              CALL polint3dlin(x1,x2,x3,table,lim,ljm,lkm,wsgcm1,sigma,
     &             soilvtrsh,ans,dy)
c     Cubic Polynomial fit (Not Used, Optional)
c               CALL polint3dlicub(x1,x2,x3,table,lim,ljm,lkm,wsgcm1,
c     &              sigma,soilvtrsh,ans,dy)
              workij2=(1.d0-mcfrac)*exp(ans)
            END IF
          END IF
          pdfint(i,j)=workij1+workij2
        END IF

        IF (sigma == 0.) THEN
          IF (wsgcm1 > soilvtrsh) THEN
            pdfint(i,j)=(wsgcm1-soilvtrsh)*wsgcm1**2.D0
          ELSE
            pdfint(i,j)=0.D0
          END IF
        END IF

        IF (pdfint(i,j) > 0.D0) THEN
          dsteve1=1.D0
        ELSE
          dsteve1=0.D0
        END IF

        IF (vtrsh(i,j) > 0. .AND. wsgcm1 > vtrsh(i,j)) THEN
          dsteve2=1.D0
        ELSE
          dsteve2=0.D0
        END IF

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

      SUBROUTINE local_dust_emission(i,j,n,wsgcm,ptype,dsrcflx,dsrcflx2)
!@sum  selects routine for calculating local dust source flux
!@auth Jan Perlwitz, Reha Cakmur, Ina Tegen
#ifndef TRACERS_AMP
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      USE constant,ONLY : sday
      USE model_com,ONLY : jday,jmon
      USE geom,ONLY : dxyp
      USE tracer_com,ONLY : trname,imDUST,n_clay,n_clayilli
#ifdef TRACERS_QUARZHEM
     &     ,FreeFe
#endif
      USE tracers_dust,ONLY : pdfint,Fracl,Frasi,frclay,frsilt,gin_data,
     &     qdust,upclsi,vtrsh,d_dust,ers_data
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
     &     ,minfr
#endif

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: i,j,n
      REAL*8,INTENT(IN) :: ptype,wsgcm
      REAL*8,INTENT(OUT) :: dsrcflx,dsrcflx2

      INTEGER :: n1
      REAL*8 :: frtrac,zfac

      IF (imDUST == 0) THEN
c     Interactive dust emission

      IF (.NOT. qdust(i,j)) THEN
        dsrcflx=0.D0
        dsrcflx2=0.D0
      ELSE
#ifdef TRACERS_DUST

#ifdef TRACERS_DUST_CUB_SAH

        SELECT CASE(trname(n))
        CASE ('Clay')
          frtrac=frclay(i,j)*Fracl
        CASE ('Silt1','Silt2','Silt3','Silt4')
          frtrac=frsilt(i,j)*Frasi
        END SELECT

        dsrcflx=Upclsi*frtrac*(wsgcm-vtrsh(i,j))*wsgcm**2

#else ! default case

        SELECT CASE(trname(n))
        CASE ('Clay')
          frtrac=Fracl
        CASE ('Silt1','Silt2','Silt3','Silt4')
          frtrac=Frasi
        END SELECT

c ..........
c dust emission above threshold from sub grid scale wind fluctuations
c ..........
        dsrcflx=Upclsi*frtrac*ers_data(i,j,jmon)*gin_data(i,j)
     &       *pdfint(i,j)
c ..........
c emission according to cubic scheme
c ..........
        IF (vtrsh(i,j) > 0. .AND. wsgcm > vtrsh(i,j)) THEN
          dsrcflx2=Upclsi*frtrac*gin_data(i,j)*ers_data(i,j,jmon)
     &         *(wsgcm-vtrsh(i,j))*wsgcm**2
        ELSE
          dsrcflx2=0.D0
        END IF

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

#endif

#ifdef TRACERS_DUST_CUB_SAH

        SELECT CASE(trname(n))
#ifdef TRACERS_MINERALS
        CASE ('ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar')
          frtrac=frclay(i,j)*Fracl*minfr(i,j,n1)
        CASE ('Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps',
     &        'Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps',
     &        'Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps')
          frtrac=frsilt(i,j)*Frasi*minfr(i,j,n1)
        CASE ('Sil1Quar','Sil2Quar','Sil3Quar')
          zfac=frsilt(i,j)*Frasi
          frtrac=zfac*minfr(i,j,n1)
#ifdef TRACERS_QUARZHEM
     &         -zfac*MIN((1.D0-FreeFe)*minfr(i,j,9),DBLE(minfr(i,j,n1)))
#endif
#endif
#ifdef TRACERS_QUARZHEM
        CASE ('Sil1QuHe','Sil2QuHe','Sil3QuHe')
          frtrac=frsilt(i,j)*Frasi
     &         *MIN((1.D0-FreeFe)*minfr(i,j,9),DBLE(minfr(i,j,6)))
#endif
        END SELECT

        dsrcflx=Upclsi*frtrac*(wsgcm-vtrsh(i,j))*wsgcm**2

#else ! default case

        SELECT CASE(trname(n))
#ifdef TRACERS_MINERALS
        CASE ('ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar')
          frtrac=Fracl*minfr(i,j,n1)
        CASE ('Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps',
     &        'Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps',
     &        'Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps')
          frtrac=Frasi*minfr(i,j,n1)
        CASE ('Sil1Quar','Sil2Quar','Sil3Quar')
          frtrac=Frasi*minfr(i,j,n1)
#ifdef TRACERS_QUARZHEM
     &         -Frasi
     &         *MIN((1.D0-FreeFe)*minfr(i,j,9),DBLE(minfr(i,j,n1)))
#endif
#endif
#ifdef TRACERS_QUARZHEM
        CASE ('Sil1QuHe','Sil2QuHe','Sil3QuHe')
          frtrac=Frasi
     &         *MIN((1.D0-FreeFe)*minfr(i,j,9),DBLE(minfr(i,j,6)))
#endif
        END SELECT

        dsrcflx=Upclsi*frtrac*ers_data(i,j,jmon)*gin_data(i,j)
     &       *pdfint(i,j)

#endif
#endif

      END IF

      ELSE IF (imDUST == 1) THEN
c     prescribed AEROCOM dust emission

      IF (.NOT. qdust(i,j)) THEN
        dsrcflx=0.D0
      ELSE

#ifdef TRACERS_DUST

        SELECT CASE(trname(n))
        CASE ('Clay','Silt1','Silt2','Silt3')
          n1=n-n_clay+1
          dsrcflx=d_dust(i,j,n1,jday)/Sday/dxyp(j)/ptype
        CASE DEFAULT
          dsrcflx=0.D0
        END SELECT

#else
#ifdef TRACERS_MINERALS
      
        SELECT CASE(trname(n))
        CASE ('ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar')
          n1=n-n_clayilli+1
          dsrcflx=d_dust(i,j,1,jday)
        CASE ('Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps')
          n1=n-n_clayilli+1
          dsrcflx=d_dust(i,j,2,jday)
        CASE ('Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps')
          n1=n-n_clayilli-4
          dsrcflx=d_dust(i,j,3,jday)
        CASE ('Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps')
          n1=n-n_clayilli-9
          dsrcflx=d_dust(i,j,4,jday)
        END SELECT

#endif

#ifdef TRACERS_QUARZHEM

        SELECT CASE(trname(n))
        CASE ('Sil1QuHe')
          dsrcflx=d_dust(i,j,2,jday)
        CASE ('Sil2QuHe')
          dsrcflx=d_dust(i,j,3,jday)
        CASE ('Sil3QuHe')
          dsrcflx=d_dust(i,j,4,jday)
        END SELECT

#endif

        SELECT CASE(trname(n))
#ifdef TRACERS_MINERALS
        CASE ('Sil1Quar','Sil2Quar','Sil3Quar')
          frtrac=minfr(i,j,n1)
#ifdef TRACERS_QUARZHEM
     &         -MIN((1.D0-FreeFe)*minfr(i,j,9),DBLE(minfr(i,j,n1)))
#endif
#endif
#ifdef TRACERS_QUARZHEM
        CASE ('Sil1QuHe','Sil2QuHe','Sil3QuHe')
          frtrac=MIN((1.D0-FreeFe)*minfr(i,j,9),DBLE(minfr(i,j,6)))
#endif
#ifdef TRACERS_MINERALS
        CASE DEFAULT
          frtrac=minfr(i,j,n1)
#endif
        END SELECT

        dsrcflx=dsrcflx*frtrac/Sday/dxyp(j)/ptype

#endif

      END IF

      END IF

#endif
# else ! TRACERS_AMP
     
      USE constant,ONLY : sday
      USE model_com,ONLY : jday,jmon
      USE geom,ONLY : dxyp
      USE tracer_com,ONLY : trname,imDUST,n_clay,n_clayilli
      USE tracers_dust,ONLY : pdfint,Fracl,Frasi,frclay,frsilt,gin_data,
     &     qdust,upclsi,vtrsh,d_dust,ers_data

      IMPLICIT NONE
      INTEGER,INTENT(IN) :: i,j,n
      REAL*8,INTENT(IN) :: ptype,wsgcm
      REAL*8,INTENT(OUT) :: dsrcflx,dsrcflx2

      INTEGER :: n1
      REAL*8 :: frtrac,zfac

      IF (imDUST == 0) THEN !Interactive dust emission
      IF (.NOT. qdust(i,j)) THEN
        dsrcflx=0.D0
        dsrcflx2=0.D0
      ELSE
        if (N.gt.4) PRINT*,'WRONG DUST EMISSION TRACER' 
        if (N.gt.4) stop 
#ifdef TRACERS_DUST_CUB_SAH
        if (n.eq.1) then
          frtrac=frclay(i,j)*Fracl
        else  
          frtrac=frsilt(i,j)*Frasi
        endif
        dsrcflx=Upclsi*frtrac*(wsgcm-vtrsh(i,j))*wsgcm**2
#else 

        if (n.eq.1) then   
          frtrac=Fracl
        else  
          frtrac=Frasi
        endif
c ..........
c dust emission above threshold from sub grid scale wind fluctuations
c ..........
        dsrcflx=Upclsi*frtrac*ers_data(i,j,jmon)*gin_data(i,j)
     &       *pdfint(i,j)
c ..........
c emission according to cubic scheme
c ..........
        IF (vtrsh(i,j) > 0. .AND. wsgcm > vtrsh(i,j)) THEN
          dsrcflx2=Upclsi*frtrac*gin_data(i,j)*ers_data(i,j,jmon)
     &         *(wsgcm-vtrsh(i,j))*wsgcm**2
        ELSE
          dsrcflx2=0.D0
        END IF
#endif

#ifdef TRACERS_DUST_CUB_SAH
        dsrcflx=Upclsi*frtrac*(wsgcm-vtrsh(i,j))*wsgcm**2
#else
        dsrcflx=Upclsi*frtrac*ers_data(i,j,jmon)*gin_data(i,j)
     &       *pdfint(i,j)
#endif

      END IF
      ELSE IF (imDUST == 1) THEN !prescribed AEROCOM dust emission
      IF (.NOT. qdust(i,j)) THEN
        dsrcflx=0.D0
      ELSE
!        n1=n-n_clay+1 n stuff needs to be sorted out SUSA
!        dsrcflx=d_dust(i,j,n1,jday)/Sday/dxyp(j)/ptype
!        dsrcflx=dsrcflx*frtrac/Sday/dxyp(j)/ptype

      END IF
      END IF

#endif ! TRACERS_AMP  

            
      RETURN
      END SUBROUTINE local_dust_emission

      SUBROUTINE polint3dlin(x1a,x2a,x3a,ya,m,n,lkm,x1,x2,x3,y,dy) 
 
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)  || (defined TRACERS_AMP)

      implicit none 
      INTEGER, INTENT(IN) :: m,n,lkm 
      REAL*8, INTENT(IN) :: x1,x2,x3,x1a(m),x2a(n),x3a(lkm),ya(m,n,lkm) 
      REAL*8, INTENT(OUT):: y,dy 
      INTEGER, PARAMETER :: nmax=2
      INTEGER i,j,k,jjj,iii,xx,yy,zz,kkk 
      REAL*8 ymtmp(nmax),yntmp(nmax),x11(nmax),x22(nmax),x33(nmax) 
      real*8 yotmp(nmax) 
 
      call locate(x1a,m,x1,xx) 
      call locate(x2a,n,x2,yy) 
      call locate(x3a,lkm,x3,zz) 

      do i=1,nmax
         kkk=i
         x33(i)=x3a(zz+kkk-1)
         do k=1,nmax
            iii=k
            x22(k)=x2a(yy+iii-1)  
            do j=1,nmax
               jjj=j
               x11(j)=x1a(xx+jjj-1)
               yntmp(j)=ya(xx+jjj-1,yy+iii-1,zz+kkk-1)
            enddo
            if (yntmp(1).eq.-1000) then
               ymtmp(k)=-1000.
            else
               call polint(x11,yntmp,nmax,x1,ymtmp(k),dy)
            endif
         enddo
         if (ymtmp(1).eq.-1000) then
            yotmp(i)=-1000.
         else
            call polint(x22,ymtmp,nmax,x2,yotmp(i),dy)
         endif
      enddo
      if (yotmp(2).eq.-1000)  then
         y=-1000.
      else
         call polint(x33,yotmp,nmax,x3,y,dy)
      endif
#endif
      return
      END SUBROUTINE POLINT3DLIN
C  (C) Copr. 1986-92 Numerical Recipes Software 'W3.

      SUBROUTINE polint3dcub(x1a,x2a,x3a,ya,m,n,lkm,x1,x2,x3,y,dy)

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      implicit none
      INTEGER, INTENT(IN) :: m,n,lkm
      REAL*8, INTENT(IN) :: x1,x2,x3,x1a(m),x2a(n),x3a(lkm),ya(m,n,lkm)
      REAL*8, INTENT(OUT):: y,dy
      INTEGER, PARAMETER :: nmax=4
      INTEGER i,j,k,jjj,iii,xx,yy,zz,kkk
      REAL*8 ymtmp(nmax),yntmp(nmax),x11(nmax),x22(nmax),x33(nmax)
      real*8 yotmp(nmax)

      call locate(x1a,m,x1,xx)
      call locate(x2a,n,x2,yy)
      call locate(x3a,lkm,x3,zz)

      do i=1,nmax
         kkk=i-1
         x33(i)=x3a(zz+kkk-1)
         do k=1,nmax
            iii=k-1
            x22(k)=x2a(yy+iii-1)
            do j=1,nmax
               jjj=j-1
               x11(j)=x1a(xx+jjj-1)
               yntmp(j)=ya(xx+jjj-1,yy+iii-1,zz+kkk-1)
            enddo
            if (yntmp(1).eq.-1000.or.yntmp(2).eq.-1000.) then
               ymtmp(k)=-1000.
            else
               call polint(x11,yntmp,nmax,x1,ymtmp(k),dy)
            endif
         enddo
         if (ymtmp(1).eq.-1000.or.ymtmp(2).eq.-1000.) then
            yotmp(i)=-1000.
         else
            call polint(x22,ymtmp,nmax,x2,yotmp(i),dy)
         endif
      enddo
      if (yotmp(3).eq.-1000.or.yotmp(4).eq.-1000.)  then
         y=-1000.
      else
         call polint(x33,yotmp,nmax,x3,y,dy)
      endif
#endif
      return
      END SUBROUTINE POLINT3DCUB
C  (C) Copr. 1986-92 Numerical Recipes Software 'W3.

      SUBROUTINE polint(xa,ya,n,x,y,dy)

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n
      REAL*8, INTENT(OUT) :: y,dy
      REAL*8, INTENT(IN) :: x,xa(n),ya(n)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(n),d(n)
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

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
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
