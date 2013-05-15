#include "rundeck_opts.h"
      SUBROUTINE dust_emission_constraints(itype,wsgcm,pbl_args)
!@sum  local constrainsts for dust tracer emission valid for all dust bins
!@auth Jan Perlwitz, Reha Cakmur, Ina Tegen

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      USE model_com,ONLY : dtsrc
      use fluxes, only : nisurf
      USE socpbl,ONLY : t_pbl_args
      USE tracers_dust,ONLY : imDust
      use wspdf_mod, only : lim,ljm,lkm,table,x1,x2,x3
      use polint_mod
      IMPLICIT NONE

      INTEGER,INTENT(IN) :: itype
      REAL*8,INTENT(IN) :: wsgcm

      type(t_pbl_args),INTENT(INOUT) :: pbl_args

      REAL*8 :: snowe,vtrsh
      REAL*8 :: dsteve1,dsteve2
      REAL*8 :: soilvtrsh
      LOGICAL :: qdust
      REAL*8 :: dryhr,hbaij,pprec,pevap,ricntd
      REAL*8 :: hbaijd,hbaijold
      LOGICAL :: pmei
      REAL*8 :: wearth,aiearth,wfcs,pdfint,wsubtke,wsubwd,wsubwm
      REAL(KIND=8) :: soilwet,sigma,ans,dy,workij1,workij2,wsgcm1,mcfrac
      CHARACTER(17) :: fname='WARNING_in_TRDUST'
      CHARACTER(25) :: subr='dust_emission_constraints'
      CHARACTER(5) :: vname1='wsgcm',vname2='sigma'
      CHARACTER(9) :: vname3='soilvtrsh'

c**** input
      snowe=pbl_args%snow
      vtrsh=pbl_args%vtrsh
      dryhr=pbl_args%dryhr
      pprec=pbl_args%pprec
      pevap=pbl_args%pevap
      hbaij=pbl_args%hbaij
      ricntd=pbl_args%ricntd
      wearth=pbl_args%wearth
      aiearth=pbl_args%aiearth
      wfcs=pbl_args%wfcs
      wsubtke=pbl_args%wsubtke
      wsubwd=pbl_args%wsubwd
      wsubwm=pbl_args%wsubwm
      mcfrac=pbl_args%mcfrac

      qdust = .false.
      dsteve1 = 0.D0
      dsteve2 = 0.D0
      soilvtrsh = 0.D0
      pdfint = 0.D0

      IF (imDUST == 2) THEN

c**** legacy dust emission
c     Checking if accumulated precipitation - evaporation
c     less/equal than Zero for a succeeding number of hours greater/equal
c     than threshold dryhr to permit dust emission

      hbaijold=hbaij
      hbaij=hbaijold+pprec/nisurf-pevap
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

      ELSE IF (imDUST == 0) THEN

c**** dust emission using probability density function of wind speed
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
              CALL check_upper_limit(sigma,x2(Ljm),fname,subr,vname2)
              CALL check_upper_limit(soilvtrsh,x3(Lkm),fname,subr
     &             ,vname3)
c     Linear Polynomial fit (Default)
              CALL polint3dlin(x1,x2,x3,table,lim,ljm,lkm,wsgcm1,sigma,
     &             soilvtrsh,ans,dy)
c     Cubic Polynomial fit (Not Used, Optional)
c              CALL polint3dlicub(x1,x2,x3,table,lim,ljm,lkm,wsgcm1, 
c     &             sigma,soilvtrsh,ans,dy)
              pdfint=exp(ans)
            ELSE
              CALL check_upper_limit(wsgcm1,x1(Lim),fname,subr,vname1)
              CALL check_upper_limit(sigma,x2(Ljm),fname,subr,vname2)
              CALL check_upper_limit(soilvtrsh,x3(Lkm),fname,subr
     &             ,vname3)
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
c     only over the area with downdrafts (mcfrac).

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
              CALL check_upper_limit(sigma,x2(Ljm),fname,subr,vname2)
              CALL check_upper_limit(soilvtrsh,x3(Lkm),fname,subr
     &             ,vname3)
c     Linear Polynomial fit (Default)
              CALL polint3dlin(x1,x2,x3,table,lim,ljm,lkm,wsgcm1,sigma,
     &             soilvtrsh,ans,dy)
c     Cubic Polynomial fit (Not Used, Optional)
c               call polint3dlicub(x1,x2,x3,table,lim,ljm,lkm,wsgcm1,
c     &              sigma,soilvtrsh,ans,dy)
              workij1=mcfrac*exp(ans)
            ELSE
              CALL check_upper_limit(wsgcm1,x1(Lim),fname,subr,vname1)
              CALL check_upper_limit(sigma,x2(Ljm),fname,subr,vname2)
              CALL check_upper_limit(soilvtrsh,x3(Lkm),fname,subr
     &             ,vname3)
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
              CALL check_upper_limit(sigma,x2(Ljm),fname,subr,vname2)
              CALL check_upper_limit(soilvtrsh,x3(Lkm),fname,subr
     &             ,vname3)
c     Linear Polynomial fit (Default)
              CALL polint3dlin(x1,x2,x3,table,lim,ljm,lkm,wsgcm1,sigma,
     &             soilvtrsh,ans,dy)
c     Cubic Polynomial fit (Not Used, Optional)
c               CALL polint3dlicub(x1,x2,x3,table,lim,ljm,lkm,wsgcm1,
c     &              sigma,soilvtrsh,ans,dy)
              workij2=(1.d0-mcfrac)*exp(ans)
            ELSE
              CALL check_upper_limit(wsgcm1,x1(Lim),fname,subr,vname1)
              CALL check_upper_limit(sigma,x2(Ljm),fname,subr,vname2)
              CALL check_upper_limit(soilvtrsh,x3(Lkm),fname,subr
     &             ,vname3)
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
      pbl_args%hbaij=hbaij
      pbl_args%ricntd=ricntd
      pbl_args%wtrsh=soilvtrsh
      pbl_args%pdfint=pdfint

#endif

      RETURN
      END SUBROUTINE dust_emission_constraints

      SUBROUTINE local_dust_emission(n,wsgcm,pbl_args,dsrcflx,
     &     dsrcflx2)
!@sum  selects routine for calculating local dust source flux
!@auth Jan Perlwitz, Reha Cakmur, Ina Tegen

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)  ||\
    (defined TRACERS_TOMAS)
      USE socpbl,ONLY : t_pbl_args
      use tracer_com, only: Ntm_dust
      use OldTracer_mod, only: trname
      use tracers_dust,only : nAerocomDust,CWiCub,FClWiCub,FSiWiCub,
     &     CWiPdf,FracClayPDFscheme,FracSiltPDFscheme,imDust
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
     &     , n_soildust
#endif

      IMPLICIT NONE

      INTEGER,INTENT(IN) :: n
      REAL*8,INTENT(IN) :: wsgcm
      TYPE(t_pbl_args),INTENT(IN) :: pbl_args

      REAL*8,INTENT(OUT) :: dsrcflx,dsrcflx2

      integer :: n1
      REAL*8 :: vtrsh
      real(kind=8) :: d_dust(nAerocomDust)
      REAL*8 :: frtrac
      LOGICAL :: qdust
      REAL*8 :: frclay,frsilt
      real(kind=8) :: ers_data,dustSourceFunction,soilvtrsh,pdfint
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
      real(kind=8) :: mineralFractions( Ntm_dust )
#endif

c**** input
      qdust=pbl_args%qdust
      vtrsh=pbl_args%vtrsh
      IF (imDust == 1) d_dust(:)=pbl_args%d_dust(:)
      frclay=pbl_args%frclay
      frsilt=pbl_args%frsilt
      ers_data=pbl_args%ers_data
      dustSourceFunction = pbl_args%dustSourceFunction
      soilvtrsh=pbl_args%wtrsh
      pdfint=pbl_args%pdfint
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
      mineralFractions( : ) = pbl_args%mineralFractions( : )
#endif

c**** initialize
      dsrcflx=0.D0
      dsrcflx2=0.D0
      IF (qdust) THEN

      IF (imDUST /= 1) THEN
c**** Interactive dust emission

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)

        SELECT CASE(trname(n))
        CASE ('Clay','ClayIlli','ClayKaol','ClaySmec','ClayCalc',
     &        'ClayQuar')
          IF (imDust == 0) THEN
            frtrac = FracClayPDFscheme
          ELSE IF (imDust == 2) THEN
            frtrac=FClWiCub*frclay
          END IF
        CASE ('Silt1','Silt2','Silt3','Silt4','Sil1Quar','Sil1Feld',
     &        'Sil1Calc','Sil1Hema','Sil1Gyps','Sil2Quar','Sil2Feld',
     &        'Sil2Calc','Sil2Hema','Sil2Gyps','Sil3Quar','Sil3Feld',
     &        'Sil3Calc','Sil3Hema','Sil3Gyps','Sil1QuHe','Sil2QuHe',
     &        'Sil3QuHe')
          IF (imDust == 0) THEN
            frtrac = FracSiltPDFscheme
          ELSE IF (imDust == 2) THEN
            frtrac=FSiWiCub*frsilt
          END IF
        case default
          return
        END SELECT

#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
        frtrac = frtrac * mineralFractions( n - n_soildust + 1 )
#endif

#else /* TRACERS_DUST || TRACERS_MINERALS || TRACERS_QUARZHEM */
#if (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
        SELECT CASE (n)
        CASE (1)
          IF (imDust == 0) THEN
            frtrac = FracClayPDFscheme
          ELSE IF (imDust == 2) THEN
            frtrac=FClWiCub*frclay
          END IF
        CASE (2,3,4)
          IF (imDust == 0) THEN
            frtrac = FracSiltPDFscheme
          ELSE IF (imDust == 2) THEN
            frtrac=FSiWiCub*frsilt
          END IF
        END SELECT
#endif /* TRACERS_AMP || TRACERS_TOMAS */

#endif

c**** legacy dust emission scheme
        IF (imDust == 2) THEN
          dsrcflx=CWiCub*frtrac*(wsgcm-vtrsh)*wsgcm**2

c**** default case
c ..........
c dust emission above threshold from sub grid scale wind fluctuations
c ..........
        ELSE IF (imDust == 0) THEN
          dsrcflx = CWiPdf*frtrac*ers_data*dustSourceFunction*pdfint
c ..........
c emission according to cubic scheme, but with pdf sheme parameters
c (only used as diagnostic variable)
c ..........
          IF (soilvtrsh > 0. .AND. wsgcm > soilvtrsh) THEN
            dsrcflx2 = CWiPdf*frtrac*dustSourceFunction*ers_data
     &           *(wsgcm-soilvtrsh)*wsgcm**2
          END IF
        END IF

      ELSE IF (imDUST == 1) THEN
c**** prescribed AEROCOM dust emission

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)

        SELECT CASE(trname(n))
        CASE ('Clay','Silt1','Silt2','Silt3')
          dsrcflx = d_dust( n - n_soildust + 1 )
        CASE ('ClayIlli','ClayKaol','ClaySmec','ClayCalc','ClayQuar')
          dsrcflx=d_dust(1)
        case( 'Sil1Quar', 'Sil1Feld', 'Sil1Calc', 'Sil1Hema', 'Sil1Gyps'
     &         , 'Sil1QuHe' )
          dsrcflx=d_dust(2)
        case( 'Sil2Quar', 'Sil2Feld', 'Sil2Calc', 'Sil2Hema', 'Sil2Gyps'
     &         , 'Sil2QuHe' )
          dsrcflx=d_dust(3)
       CASE( 'Sil3Quar', 'Sil3Feld', 'Sil3Calc', 'Sil3Hema', 'Sil3Gyps'
     &         , 'Sil3QuHe')
          dsrcflx=d_dust(4)
        END SELECT

#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
        dsrcflx = dsrcflx * mineralFractions( n - n_soildust + 1 )
#endif

#else /* TRACERS_DUST || TRACERS_MINERALS || TRACERS_QUARZHEM */

#if (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
        dsrcflx=d_dust(n)
#endif

#endif

      END IF

      END IF

#endif
            
      return

      END SUBROUTINE local_dust_emission

c****
c**** This is legacy code for the old wet deposition scheme.
c**** To use the code the comiler directive TRACERS_WATER must
c**** not be defined and WET_DEPO_Ina must be defined in rundeck.
c****
      SUBROUTINE dust_wet(i,j)
!@sum  Simple scheme for wet deposition of dust/mineral tracers
!@auth Ina Tegen, Reha Cakmur, Jan Perlwitz

#ifndef TRACERS_WATER
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
      USE constant,ONLY : Grav
      USE resolution,ONLY : Jm,Lm
      USE atm_com,ONLY : zatmo,gz
      USE fluxes,ONLY : prec
      USE clouds,ONLY : tm_dust,tmom_dust,trprc_dust
      USE tracer_com,ONLY : Ntm_dust,trname
      USE tracers_dust,ONLY : prelay

      IMPLICIT NONE

      REAL*8,PARAMETER :: Z=700.

      INTEGER,INTENT(IN) :: i,j

      INTEGER :: l,layer,n
      INTEGER,DIMENSION(Jm) :: lwdep
      REAL*8,DIMENSION(Jm) :: h
      REAL*8 :: y
      REAL*8 :: height
      REAL*8,DIMENSION(Lm,Ntm_dust) :: tmold

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
