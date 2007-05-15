#include "rundeck_opts.h"
      SUBROUTINE dust_emission_constraints(itype,ptype,wsgcm,pbl_args)
!@sum  local constrainsts for dust tracer emission valid for all dust bins
!@auth Jan Perlwitz, Reha Cakmur, Ina Tegen

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      USE model_com,ONLY : dtsrc,nisurf
      USE socpbl,ONLY : t_pbl_args
      USE tracer_com,ONLY : imDust
      USE tracers_dust,ONLY : lim,ljm,lkm,table,x1,x2,x3

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: itype
      REAL*8,INTENT(IN) :: ptype,wsgcm

      type(t_pbl_args),INTENT(INOUT) :: pbl_args

      REAL*8 :: snowe,vtrsh
      REAL*8 :: dsteve1,dsteve2
      REAL*8 :: soilvtrsh
      LOGICAL :: qdust
#ifdef TRACERS_DUST_CUB_SAH
      REAL*8 :: dryhr,hbaij,pprec,pevap,ricntd
      REAL*8 :: hbaijd,hbaijold
      LOGICAL :: pmei
#else
      REAL*8 :: wearth,aiearth,wfcs,pdfint,wsubtke,wsubwd,wsubwm
      REAL*8 :: mcfrac=0.05
      REAL*8 :: soilwet,sigma,ans,dy,workij1,workij2,wsgcm1
#endif

c**** input
      snowe=pbl_args%snowe
      vtrsh=pbl_args%vtrsh
#ifdef TRACERS_DUST_CUB_SAH
      dryhr=pbl_args%dryhr
      pprec=pbl_args%pprec
      pevap=pbl_args%pevap
      hbaij=pbl_args%hbaij
      ricntd=pbl_args%ricntd
#else
      wearth=pbl_args%wearth
      aiearth=pbl_args%aiearth
      wfcs=pbl_args%wfcs
      wsubtke=pbl_args%wsubtke
      wsubwd=pbl_args%wsubwd
      wsubwm=pbl_args%wsubwm
#endif

      IF (imDUST == 0) THEN

#ifdef TRACERS_DUST_CUB_SAH
c     Checking if accumulated precipitation - evaporation
c     less/equal than Zero for a succeeding number of hours greater/equal
c     than threshold dryhr to permit dust emission

      hbaijold=hbaij
      hbaij=hbaijold+pprec*ptype/nisurf-pevap
      hbaijd=hbaij-hbaijold
      IF (itype == 4 .AND. hbaijd <= 0.D0) THEN
        ricntd=ricntd+Dtsrc/3600./nisurf
        IF (ricntd >= dryhr .AND. dryhr /= 0.D0) THEN
          pmei=.TRUE.
        ELSE
          pmei=.FALSE.
        END IF
      ELSE
        ricntd=0.D0
        pmei=.FALSE.
      END IF

      IF (vtrsh > 0.D0 .AND. wsgcm > vtrsh) THEN
        dsteve2=1.D0
      ELSE
        dsteve2=0.D0
      END IF
      IF (pmei .AND. snowe <= 1 .AND. vtrsh > 0.D0 .AND. wsgcm > vtrsh)
     &     THEN
        dsteve1=1.D0
        qdust=.TRUE.
      ELSE
        dsteve1=0.D0
        qdust=.FALSE.
      END IF

#else
c**** default case

      IF (itype == 4 .AND. snowe <= 1.D0) THEN
        qdust=.TRUE.
      ELSE
        qdust=.FALSE.
        dsteve1=0.D0
        dsteve2=0.D0
        soilvtrsh=0.D0
      END IF

      IF (qdust) THEN

        soilwet=(wearth+aiearth)/(wfcs+1.D-20)
        if (soilwet.gt.1.D0) soilwet=1.d0
        soilvtrsh=8.d0*(exp(0.7d0*soilwet))

        pdfint=0.d0
        workij1=0.d0
        workij2=0.d0
        wsgcm1=wsgcm

c     There is no moist convection, sigma is composed of TKE and DRY
c     convective velocity scale
        IF (wsubwm == 0.D0) THEN
          sigma=wsubtke+wsubwd
c     No need to calculate the emission below these values since
c     the emission is zero
          IF (sigma > 0.1D0 .OR. wsgcm1 > 1.D0) THEN
c     This is the case when sigma is very small and we approximate
c     the function by standard dust emission formula
            IF (sigma < 0.0005D0 .AND. wsgcm1 > 1.D0) THEN
              IF (wsgcm1 > soilvtrsh) THEN
                pdfint=(wsgcm1-soilvtrsh)*wsgcm1**2.D0
              ELSE
                pdfint=0.d0
              END IF
c     This is the case when wsgcm1 is very small and we set it
c     equal to one of the smallest values in the table index
            ELSE IF (sigma > 0.1D0 .AND. wsgcm1 < 0.0005D0) THEN
              wsgcm1=0.0005d0
c     Linear Polynomial fit (Default)
              CALL polint3dlin(x1,x2,x3,table,lim,ljm,lkm,wsgcm1,sigma,
     &             soilvtrsh,ans,dy)
c     Cubic Polynomial fit (Not Used, Optional)
c              CALL polint3dlicub(x1,x2,x3,table,lim,ljm,lkm,wsgcm1, 
c     &             sigma,soilvtrsh,ans,dy)
              pdfint=exp(ans)
            ELSE
c     Linear Polynomial fit (Default)
              CALL polint3dlin(x1,x2,x3,table,lim,ljm,lkm,wsgcm1,sigma,
     &             soilvtrsh,ans,dy)
c     Cubic Polynomial fit (Not Used, Optional) 
c              CALL polint3dlicub(x1,x2,x3,table,lim,ljm,lkm,wsgcm1,
c     &             sigma,soilvtrsh,ans,dy)
              pdfint=exp(ans)
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
          IF (sigma > 0.1D0 .OR. wsgcm1 > 1.D0) THEN
c     This is the case when sigma is very small and we approximate
c     the function by standard dust emission formula
            IF (sigma < 0.0005D0 .AND. wsgcm1 > 1.D0) THEN
              IF (wsgcm1 > soilvtrsh) THEN
                workij1=mcfrac*(wsgcm1-soilvtrsh)*wsgcm1**2.D0
              ELSE
                workij1=0.d0
              END IF
c     This is the case when wsgcm1 is very small and we set it
c     equal to one of the smallest values in the table index
            ELSE IF (sigma > 0.1D0 .AND. wsgcm1 < 0.0005D0) THEN
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
          IF (sigma > 0.1D0 .OR. wsgcm1 > 1.D0) THEN
c     This is the case when sigma is very small and we approximate
c     the function by standard dust emission formula
            IF (sigma < 0.0005D0 .AND. wsgcm1 > 1.D0) THEN
              IF (wsgcm1 > soilvtrsh) THEN
                workij2=(1.d0-mcfrac)*(wsgcm1-soilvtrsh)*wsgcm1**2.D0
              ELSE
                workij2=0.d0
              END IF
c     This is the case when wsgcm1 is very small and we set it
c     equal to one of the smallest values in the table index
            ELSE IF (sigma > 0.1D0 .AND. wsgcm1 < 0.0005D0) THEN
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
          pdfint=workij1+workij2
        END IF

        IF (sigma == 0.D0) THEN
          IF (wsgcm1 > soilvtrsh) THEN
            pdfint=(wsgcm1-soilvtrsh)*wsgcm1**2.D0
          ELSE
            pdfint=0.D0
          END IF
        END IF

        IF (pdfint > 0.D0) THEN
          dsteve1=1.D0
        ELSE
          dsteve1=0.D0
        END IF

        IF (vtrsh > 0.D0 .AND. wsgcm1 > vtrsh) THEN
          dsteve2=1.D0
        ELSE
          dsteve2=0.D0
        END IF

      END IF

#endif

      ELSE IF (imDUST == 1) THEN
        IF (itype == 4) THEN
          qdust=.TRUE.
        ELSE
          qdust=.FALSE.
        END IF
        soilvtrsh=0.D0
      END IF

c**** output
      pbl_args%dust_event1=dsteve1
      pbl_args%dust_event2=dsteve2
      pbl_args%qdust=qdust
#ifdef TRACERS_DUST_CUB_SAH
      pbl_args%hbaij=hbaij
      pbl_args%ricntd=ricntd
#else
      pbl_args%wtrsh=soilvtrsh
      pbl_args%pdfint=pdfint
#endif

#endif

      RETURN
      END SUBROUTINE dust_emission_constraints

      SUBROUTINE local_dust_emission(n,ptype,wsgcm,pbl_args,dsrcflx,
     &     dsrcflx2)
!@sum  selects routine for calculating local dust source flux
!@auth Jan Perlwitz, Reha Cakmur, Ina Tegen

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      USE socpbl,ONLY : t_pbl_args
      USE tracer_com,ONLY : Ntm_dust,trname,imDUST,n_clay,n_clayilli
#ifdef TRACERS_QUARZHEM
     &     ,FreeFe
#endif
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
      USE tracers_dust,ONLY : Mtrac
#endif
      USE tracers_dust,ONLY : Fracl,Frasi,upclsi

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: n
      REAL*8,INTENT(IN) :: ptype,wsgcm
      TYPE(t_pbl_args),INTENT(IN) :: pbl_args

      REAL*8,INTENT(OUT) :: dsrcflx,dsrcflx2

      INTEGER :: n1
      REAL*8 :: vtrsh
      REAL*8 :: d_dust(Ntm_dust)
      REAL*8 :: frtrac
      LOGICAL :: qdust
#ifdef TRACERS_DUST_CUB_SAH
      REAL*8 :: frclay,frsilt
#else
      REAL*8 :: ers_data,src_fnct,pdfint
#endif
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
      REAL*8 :: minfr(Mtrac)
#endif

c**** input
      qdust=pbl_args%qdust
      vtrsh=pbl_args%vtrsh
      IF (imDust == 1) d_dust(:)=pbl_args%d_dust(:)
#ifdef TRACERS_DUST_CUB_SAH
      frclay=pbl_args%frclay
      frsilt=pbl_args%frsilt
#else
      ers_data=pbl_args%ers_data
      src_fnct=pbl_args%src_fnct
      pdfint=pbl_args%pdfint
#endif
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
      minfr(:)=pbl_args%minfr(:)
#endif

c**** initialize
      dsrcflx=0.D0
      dsrcflx2=0.D0
      IF (qdust) THEN

      IF (imDUST == 0) THEN
c**** Interactive dust emission

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)

        SELECT CASE(trname(n))
        CASE ('Clay','ClayIlli','ClayKaol','ClaySmec','ClayCalc',
     &        'ClayQuar')
          frtrac=Fracl
#ifdef TRACERS_DUST_CUB_SAH
     &         *frclay
#endif
        CASE ('Silt1','Silt2','Silt3','Silt4','Sil1Quar','Sil1Feld',
     &        'Sil1Calc','Sil1Hema','Sil1Gyps','Sil2Quar','Sil2Feld',
     &        'Sil2Calc','Sil2Hema','Sil2Gyps','Sil3Quar','Sil3Feld',
     &        'Sil3Calc','Sil3Hema','Sil3Gyps','Sil1QuHe','Sil2QuHe',
     &        'Sil3QuHe')
          frtrac=Frasi
#ifdef TRACERS_DUST_CUB_SAH
     &         *frsilt
#endif
        END SELECT

#ifdef TRACERS_MINERALS
        SELECT CASE(trname(n))
        CASE ('ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar')
          n1=n-n_clayilli+1
        CASE ('Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps')
          n1=n-n_clayilli+1
        CASE ('Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps')
          n1=n-n_clayilli-4
        CASE ('Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps')
          n1=n-n_clayilli-9
        END SELECT
#endif

        SELECT CASE(trname(n))
#ifdef TRACERS_MINERALS
        CASE ('Sil1Quar','Sil2Quar','Sil3Quar')
          frtrac=frtrac*minfr(n1)
#ifdef TRACERS_QUARZHEM
     &         -frtrac*MIN((1.D0-FreeFe)*minfr(9),minfr(n1))
#endif
#endif
#ifdef TRACERS_QUARZHEM
        CASE ('Sil1QuHe','Sil2QuHe','Sil3QuHe')
          frtrac=frtrac*MIN((1.D0-FreeFe)*minfr(9),minfr(6))
#endif
#ifdef TRACERS_MINERALS
        CASE DEFAULT
          frtrac=frtrac*minfr(n1)
#endif
        END SELECT

#else
#ifdef TRACERS_AMP
        SELECT CASE (n)
        CASE (1)
          frtrac=Fracl
#ifdef TRACERS_DUST_CUB_SAH
     &         *frclay
#endif
        CASE (2,3,4)
          frtrac=Frasi
#ifdef TRACERS_DUST_CUB_SAH
     &         *frsilt
#endif
        END SELECT
#endif

#endif

#ifdef TRACERS_DUST_CUB_SAH
c**** legacy dust emission scheme
        dsrcflx=Upclsi*frtrac*(wsgcm-vtrsh)*wsgcm**2
#else
c**** default case
c ..........
c dust emission above threshold from sub grid scale wind fluctuations
c ..........
        dsrcflx=Upclsi*frtrac*ers_data*src_fnct*pdfint
c ..........
c emission according to cubic scheme (diagnostics only)
c ..........
        IF (vtrsh > 0. .AND. wsgcm > vtrsh) THEN
          dsrcflx2=Upclsi*frtrac*src_fnct*ers_data*(wsgcm-vtrsh)
     &         *wsgcm**2
        END IF
#endif


      ELSE IF (imDUST == 1) THEN
c**** prescribed AEROCOM dust emission

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)

#ifdef TRACERS_DUST

        SELECT CASE(trname(n))
        CASE ('Clay','Silt1','Silt2','Silt3')
          n1=n-n_clay+1
          dsrcflx=d_dust(n1)
        END SELECT

#else
#ifdef TRACERS_MINERALS
      
        SELECT CASE(trname(n))
        CASE ('ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar')
          n1=n-n_clayilli+1
          dsrcflx=d_dust(1)
        CASE ('Sil1Quar','Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps')
          n1=n-n_clayilli+1
          dsrcflx=d_dust(2)
        CASE ('Sil2Quar','Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps')
          n1=n-n_clayilli-4
          dsrcflx=d_dust(3)
        CASE ('Sil3Quar','Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps')
          n1=n-n_clayilli-9
          dsrcflx=d_dust(4)
        END SELECT

#endif

#ifdef TRACERS_QUARZHEM

        SELECT CASE(trname(n))
        CASE ('Sil1QuHe')
          dsrcflx=d_dust(2)
        CASE ('Sil2QuHe')
          dsrcflx=d_dust(3)
        CASE ('Sil3QuHe')
          dsrcflx=d_dust(4)
        END SELECT

#endif

        SELECT CASE(trname(n))
#ifdef TRACERS_MINERALS
        CASE ('Sil1Quar','Sil2Quar','Sil3Quar')
          frtrac=minfr(n1)
#ifdef TRACERS_QUARZHEM
     &         -MIN((1.D0-FreeFe)*minfr(9),minfr(n1))
#endif
#endif
#ifdef TRACERS_QUARZHEM
        CASE ('Sil1QuHe','Sil2QuHe','Sil3QuHe')
          frtrac=MIN((1.D0-FreeFe)*minfr(9),minfr(6))
#endif
#ifdef TRACERS_MINERALS
        CASE DEFAULT
          frtrac=minfr(n1)
#endif
        END SELECT

        dsrcflx=dsrcflx*frtrac

#endif

#else

#ifdef TRACERS_AMP
        dsrcflx=d_dust(n)
#endif

#endif

      END IF

      END IF

#endif
            
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
