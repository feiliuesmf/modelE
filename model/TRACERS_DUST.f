#include "rundeck_opts.h"
      SUBROUTINE tracers_dust_old
!@sum soil dust sources and sinks
!auth Reha Cakmur, Jan Perlwitz, Ina Tegen

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)

c      do j=1,jm
c         do i=1,im
c            if (bldata(i,j,1) > vtrsh(i,j)tadd(i,j,1)=tadd(i,j,1)+1.
c            if (frclay(i,j)*pmei(i,j) > 0)tadd(i,j,3)=tadd(i,j,3)+1.
c         enddo
c      enddo

#ifdef DUST_EMISSION_EXTERN
      CALL dust_emission
#endif
#ifndef TRACERS_WATER
      CALL dust_wet
#endif
#ifndef TRACERS_DRYDEP
      CALL dust_turb
      CALL dust_grav
#endif
#endif

      RETURN
      END SUBROUTINE tracers_dust_old

      SUBROUTINE dust_emission
!@sum  dust source flux
!@auth Jan Perlwitz, Reha Cakmur, Ina Tegen

#ifdef TRACERS_DUST
      USE model_com,ONLY : Dtsrc,fearth,im,jm,jmon,wfcs,itime
      USE fluxes,ONLY : prec,evapor,trsrfflx
      USE tracer_com,ONLY : n_clay,Ntm_dust,trm
      USE trdiag_com,ONLY : ijts_source,jls_source,taijs,tajls
      USE geom,ONLY : dxyp,imaxj
      USE ghy_com,ONLY : snowe !earth snow amount
      USE pblcom,ONLY : wsavg
      USE tracers_dust,ONLY : dryhr,frclay,frsilt,hbaij,ricntd,vtrsh,
     &     ers_data,gin_data,wsubtke=>wsubtke_com,wsubwd=>wsubwd_com,
     &     wsubwm=>wsubwm_com
      USE trdiag_com,ONLY : nDustEmij,nDustEmjl
      USE ghy_com,ONLY : wearth,aiearth

      IMPLICIT NONE

      INTEGER :: i,j,n,n1,naij,najl
      REAL*8 :: hbaijold,hbaijd,workijn,soilwet,soilvtrsh
      REAL*8 :: dsrcflx(Ntm_dust)
c     dsrcflx  dust source flux for Ntm_dust tracers [kg/s]
      LOGICAL :: pmei(Im,Jm)

#ifdef TRACERS_DUST_CUB_SAH
c     Checking whether accumulated precipitation - evaporation
c     less/equal than Zero for a succeeding number of hours greater/equal
c     than threshold dryhr to permit dust emission
      DO j=1,Jm
        DO i=1,Imaxj(j)
          hbaijold=hbaij(i,j)
          hbaij(i,j)=hbaijold+prec(i,j)*fearth(i,j)-evapor(i,j,4)
          hbaijd=hbaij(i,j)-hbaijold
          IF (hbaijd <= 0) THEN
            ricntd(i,j)=ricntd(i,j)+1.
            IF (ricntd(i,j) >= dryhr(i,j) .AND. dryhr(i,j) /= 0) THEN
              pmei(i,j)=.TRUE.
            ELSE
              pmei(i,j)=.FALSE.
            END IF
          ELSE
            ricntd(i,j)=0.
            pmei(i,j)=.FALSE.
          END IF
        END DO
      END DO
#endif

c     Loop for calculating dust source flux for each tracer
      DO n=1,Ntm_dust
        n1=n_clay+n-1
        trsrfflx(:,:,n1)=0D0
      END DO

      DO j=1,Jm
        DO i=1,Imaxj(j)

          dsrcflx=0D0

#ifdef TRACERS_DUST_CUB_SAH

          IF (fearth(i,j) > 0. .AND. snowe(i,j) <= 1 .AND.
     &         vtrsh(i,j) > 0. .AND. wsavg(i,j) > vtrsh(i,j) .AND.
     &         pmei(i,j)) THEN

#ifdef TRACERS_DUST_MINERAL8

            CALL loc_dsrcflx_cub_min8

#else

            CALL loc_dsrcflx_cub_sah(dxyp(j),wsavg(i,j),vtrsh(i,j),
     &           frclay(i,j),frsilt(i,j),fearth(i,j),dsrcflx)

#endif

          ENDIF

#else

          IF (fearth(i,j) > 0. .AND. snowe(i,j) <= 1
     &         .AND. ers_data(i,j,jmon) < -13. .AND.
     &         ers_data(i,j,jmon) /= -99.) THEN

            soilwet=(WEARTH(I,J)+AIEARTH(I,J))/(WFCS(I,J)+1.D-20)
            if (soilwet.gt.1.) soilwet=1.d0
            soilvtrsh=8.d0*(exp(0.25d0*soilwet))

c            WRITE(*,*) 'In dust_emission: itime,i,j,wearth,aiearth,',
c     &           'wfcs,soilwet:',itime,i,j,wearth(i,j),aiearth(i,j),
c     &           wfcs(i,j),soilwet

#ifdef TRACERS_DUST_MINERAL8

            CALL loc_dsrcflx_turb_min8

#else
c     default case

            CALL loc_dsrcflx_turb_sah(dxyp(j),wsavg(i,j),soilvtrsh,
     &           fearth(i,j),gin_data(i,j),wsubtke(i,j),wsubwd(i,j),
     &           wsubwm(i,j),dsrcflx)

c            WRITE(*,*) 'In dust_emission: itime,i,j,dsrcflx:',itime,i,j,
c     &           dsrcflx

          END IF

#endif
#endif

          DO n=1,Ntm_dust
            n1=n_clay+n-1
            trsrfflx(i,j,n1)=dsrcflx(n)
          END DO

        END DO
      END DO

      DO n=1,Ntm_dust
        n1=n_clay+n-1
        naij=ijts_source(nDustEmij,n1)
        najl=jls_source(nDustEmjl,n1)
c        WRITE(*,*) 'naij,najl:',naij,najl
        DO j=1,Jm
          DO i=1,Imaxj(j)
            IF (trsrfflx(i,j,n1) > 0.) THEN
              workijn=trsrfflx(i,j,n1)*Dtsrc
              trm(i,j,1,n1)=trm(i,j,1,n1)+workijn

              taijs(i,j,naij)=taijs(i,j,naij)+workijn
            END IF
          END DO
          tajls(j,1,najl)=tajls(j,1,najl)+SUM(trsrfflx(:,j,n1))*Dtsrc
        END DO
        trsrfflx(:,:,n1)=0D0
      END DO
#endif
            
      RETURN
      END SUBROUTINE dust_emission

      SUBROUTINE loc_dsrcflx_cub_sah(dxypj,wsavgij,vtrshij,frclayij,
     &     frsiltij,fearthij,dsrcflx)
!@sum  local dust source flux physics according to Ina's old cubic scheme
!@auth Ina Tegen, Jan Perlwitz, Reha Cakmur

#ifdef TRACERS_DUST
      USE tracer_com,ONLY : Ntm_dust
      USE tracers_dust_com,ONLY : Fracn,Uplfac

      IMPLICIT NONE

      REAL*4,INTENT(IN) :: frclayij,frsiltij,vtrshij
c     frclayij fraction of clay in grid cell i,j [1]
c     frsiltij fraction of silt in grid cell i,j [1]
c     vtrshij  wind speed threshold in grid cell i,j [m/s]
c     dxypj    area of grid cell at j [m**2]
      REAL*8,INTENT(IN) :: dxypj,fearthij,wsavgij
c     fearthij  fraction of land area in grid cell i,j [1]
c     wsavgij  wind speed at surface in grid cell i,j [m/s]
      REAL*8,INTENT(OUT) :: dsrcflx(Ntm_dust)
c     dsrcflx  dust source flux for Ntm_dust tracers [kg/s]
      INTEGER :: n
      REAL*8 :: workij

      DO n=1,Ntm_dust

        workij=dxypj*fearthij*(wsavgij-vtrshij)*wsavgij**2
        SELECT CASE(n)
        CASE(1)
          dsrcflx(n)=Uplfac(n)*frclayij*Fracn(n)*workij
        CASE(2:Ntm_dust)
          dsrcflx(n)=Uplfac(n)*frsiltij*Fracn(n)*workij
        END SELECT

      END DO
#endif

      RETURN
      END SUBROUTINE loc_dsrcflx_cub_sah

      SUBROUTINE loc_dsrcflx_cub_min8
!@sum  local dust source flux physics with Ina's cubic scheme and 8 minerals
!@auth Jan Perlwitz, ...

      RETURN
      END SUBROUTINE loc_dsrcflx_cub_min8

      SUBROUTINE loc_dsrcflx_turb_sah(dxypj,wsavgij,vtrshij,
     &     fearthij,gin_dataij,wsubtkeij,wsubwdij,wsubwmij,
     &     dsrcflx)
!@sum  local dust source flux physics with turbulent fluxes and Sahara dust
!@auth Reha Cakmur, Jan Perlwitz

#ifdef TRACERS_DUST
      USE tracer_com,ONLY : Ntm_dust
      USE tracers_dust,ONLY : lim,ljm,Fracn,Uplfac,table,x1,x2,x3
      USE model_com, ONLY : itime
      INTEGER :: i,j,n
      REAL*8 sigma,ans,dy
      REAL*8,INTENT(OUT) :: dsrcflx(Ntm_dust)
      REAL*4,INTENT(IN) :: gin_dataij
      REAL*8,INTENT(IN) :: dxypj,fearthij,wsavgij,vtrshij
      REAL*8,INTENT(IN) :: wsubtkeij,wsubwdij,wsubwmij
      REAL*8 :: workij,workij1,workij2

      do kk=1,9
        x3(kk)=5.d0+1.d0*kk
      enddo


      do j=1,ljm
        if (j.le.5) then
          x2(j)=.001d0*j
        endif
        if (j.gt.5.and.j.le.104) then
          x2(j)=0.005d0*j-0.02d0
        endif
        if (j.gt.104.and.j.le.194) then
          x2(j)=0.05d0*j-4.7d0
        endif
        if (j.gt.194) then
          x2(j)=0.5d0*j-92.d0
        endif
      enddo
      do i=1,lim
        if (i.le.51) then
          x1(i)=.0001d0*(i-1)
        endif
        if (i.gt.51.and.i.le.59) then
          x1(i)=0.005d0*i-.25d0
        endif
        if (i.gt.59.and.i.le.258) then
          x1(i)=0.05d0*i-2.95d0
        endif
        if (i.gt.258.and.i.le.278) then
          x1(i)=0.5d0*i-119.5d0
        endif
        if (i.gt.278) then
          x1(i)=1.d0*i-259.d0
        endif
      enddo

      workij=0.d0
      workij1=0.d0
      workij2=0.d0
c
c       There is no moist convection, sigma is composed of TKE and DRY convective
c       velocity scale
      if (wsubwmij == 0.) then
        sigma=wsubtkeij+wsubwdij
        if (sigma > 0.1 .OR. wsavgij > 1.) then
          call ratint2(x1,x2,x3,table,lim,ljm,wsavgij,sigma,vtrshij,ans,
     &         dy)
          workij=exp(ans)
        endif
      endif

c       When there is moist convection, the sigma is the combination of all
c       three subgrid scale parameters (i.e. independent or dependent)
c       Takes into account that the moist convective velocity scale acts
c       only over 5% of the area.

      if (wsubwmij /= 0.) then
        sigma=wsubtkeij+wsubwdij+wsubwmij
        if (sigma > 0.1 .OR. wsavgij > 1.) then
          call ratint2(x1,x2,x3,table,lim,ljm,wsavgij,sigma,vtrshij,ans,
     &         dy)
          workij1=exp(ans)*0.05 !!!0.05 for the MC area
        endif

        sigma=wsubtkeij+wsubwdij
        if (sigma > 0.1 .OR. wsavgij > 1.) then
          call ratint2(x1,x2,x3,table,lim,ljm,wsavgij,sigma,vtrshij,ans,
     &         dy)
          workij2=exp(ans)*0.95 !!!0.95 for the rest
        endif
        workij=workij1+workij2
      endif

      if (sigma == 0.) then
        if (wsavgij > vtrshij) then
          workij=(wsavgij-vtrshij)*wsavgij**2
        else
          workij=0.
        endif
      endif


      DO n=1,Ntm_dust
        SELECT CASE(n)
        CASE(1)
          dsrcflx(n)=dxypj*gin_dataij*fearthij*Uplfac(n)*Fracn(n)*workij
        CASE(2:Ntm_dust)
          dsrcflx(n)=dxypj*gin_dataij*fearthij*Uplfac(n)*Fracn(n)*workij
        END SELECT
      END DO
c      WRITE(*,*) 'loc_dsrcflx_turb_sah: sigma,vtrshij,wsavgij,curint:',
c     &       sigma,vtrshij,wsavgij,workij
#endif

      RETURN
      END SUBROUTINE loc_dsrcflx_turb_sah

      SUBROUTINE loc_dsrcflx_turb_min8
!@sum  local dust source flux physics with turbulent fluxes and 8 minerals
!@auth Jan Perlwitz, Reha Cakmur

      RETURN
      END SUBROUTINE loc_dsrcflx_turb_min8

      SUBROUTINE dust_wet
!@sum  Computes dust wet deposition
!@auth Reha Cakmur, Jan Perlwitz, Ina Tegen

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      USE constant,ONLY : Grav
      USE resolution,ONLY : Im,Jm,Lm
      USE model_com,ONLY : Dtsrc,zatmo
      USE fluxes,ONLY : prec,trsrfflx,tr3Dsource
      USE tracer_com,ONLY : n_clay,Ntm_dust,trm,trmom
      USE trdiag_com,ONLY : ijts_source,jls_3Dsource,taijs,tajls
      USE tracers_dust_com,ONLY : Z
      USE tracers_dust,ONLY : nDustWetij,nDustWet3Djl,prelay
      USE CLOUDS_COM, ONLY : cldmc,cldss
      USE DYNAMICS, ONLY : gz

      IMPLICIT NONE

      INTEGER :: i,j,l,n,n1,naij,najl,layer
      INTEGER,DIMENSION(Jm) :: lwdep
      REAL*8,DIMENSION(jm) :: h
      REAL*8 :: y
      REAL*8 :: cloudlay,sum1,cloudfrac(im,jm),height,cldmc1(im,jm,lm)

#ifdef WET_DEPO_Ina
      DO J=1,6
         LWDEP(J) = 3
         H(J) = 2800
         LWDEP(JM+1-J) = 3
         H(JM+1-J) = 2800
      ENDDO
      DO J=7,12
         LWDEP(J) = 4
         H(J) = 4900
         LWDEP(JM+1-J) = 4
         H(JM+1-J) = 4900
      ENDDO
      DO J=13,16
         LWDEP(J) = 5
         H(J) = 7400
         LWDEP(JM+1-J) = 5
         H(JM+1-J) = 7400
      ENDDO
      DO J=17,23
         LWDEP(J) = 6
         H(J) = 10300
         LWDEP(JM+1-J) = 6
         H(JM+1-J) = 10300
      ENDDO
      
c**** Wet Deposition
      DO n=1,Ntm_dust
        n1=n_clay+n-1
        tr3Dsource(:,:,:,nDustWet3Djl,n1)=0D0
        trsrfflx(:,:,n1)=0D0
        DO j=1,Jm
          DO l=1,lwdep(j)
            DO i=1,Im
              y = Z*(prec(i,j)/h(j))
              IF (y > 1.) y=1.
              tr3Dsource(i,j,l,nDustWet3Djl,n1)=y*trm(i,j,l,n1)/Dtsrc
              trsrfflx(i,j,n1)=trsrfflx(i,j,n1)+
     &             tr3Dsource(i,j,l,nDustWet3Djl,n1)
              trm(i,j,l,n1)=trm(i,j,l,n1)*(1.-y)
              trmom(:,i,j,l,n1)=trmom(:,i,j,l,n1)*(1.-y)
            END DO
          END DO
        END DO
      END DO

      DO n=1,Ntm_dust
        n1=n_clay+n-1
        naij=ijts_source(nDustWetij,n1)
        najl=jls_3Dsource(nDustWet3Djl,n1)
        taijs(:,:,naij)=taijs(:,:,naij)+trsrfflx(:,:,n1)*Dtsrc
        DO j=1,Jm
          DO l=1,lwdep(j)
            tajls(j,l,najl)=tajls(j,l,najl)+
     &           SUM(tr3Dsource(:,j,l,nDustWet3Djl,n1))*Dtsrc
          END DO
        END DO
      END DO
#else

      DO n=1,Ntm_dust
        n1=n_clay+n-1
        tr3Dsource(:,:,:,nDustWet3Djl,n1)=0D0
        trsrfflx(:,:,n1)=0D0
        DO j=1,Jm
          DO i=1,im

            layer=0
            do l=lm,1,-1
              if (prelay(i,j,l).ne.0.) then 
                layer=l
                exit
              endif
            enddo

            if (layer.eq.1) then
              height=(gz(i,j,layer)-zatmo(i,j))/Grav
              y = Z*(prec(i,j)/height)
              IF (y > 1.) y=1.
              DO l=1,layer
                tr3Dsource(i,j,l,nDustWet3Djl,n1)=y*trm(i,j,l,n1)/Dtsrc
                trsrfflx(i,j,n1)=trsrfflx(i,j,n1)+
     &               tr3Dsource(i,j,l,nDustWet3Djl,n1)
                trm(i,j,l,n1)=trm(i,j,l,n1)*(1.-y)
                trmom(:,i,j,l,n1)=trmom(:,i,j,l,n1)*(1.-y)
              END DO

              naij=ijts_source(nDustWetij,n1)   
              najl=jls_3Dsource(nDustWet3Djl,n1)
              taijs(i,j,naij)=taijs(i,j,naij)+trsrfflx(i,j,n1)*Dtsrc
              DO l=1,Lm
                tajls(j,l,najl)=tajls(j,l,najl)+
     &               tr3Dsource(i,j,l,nDustWet3Djl,n1)*Dtsrc
              END DO

            else if (layer.ne.0) then
              height=(gz(i,j,layer)-gz(i,j,1))/Grav
              y = Z*(prec(i,j)/height)
              IF (y > 1.) y=1.
              DO l=1,layer
                tr3Dsource(i,j,l,nDustWet3Djl,n1)=y*trm(i,j,l,n1)/Dtsrc
                trsrfflx(i,j,n1)=trsrfflx(i,j,n1)+
     &               tr3Dsource(i,j,l,nDustWet3Djl,n1)
                trm(i,j,l,n1)=trm(i,j,l,n1)*(1.-y)
                trmom(:,i,j,l,n1)=trmom(:,i,j,l,n1)*(1.-y)
              END DO
              naij=ijts_source(nDustWetij,n1)   
              najl=jls_3Dsource(nDustWet3Djl,n1)
              taijs(i,j,naij)=taijs(i,j,naij)+trsrfflx(i,j,n1)*Dtsrc
              DO l=1,Lm
                tajls(j,l,najl)=tajls(j,l,najl)+
     &               tr3Dsource(i,j,l,nDustWet3Djl,n1)*Dtsrc
              END DO
            endif

          END DO
        END DO
      END DO

#endif

#endif

      RETURN
      END SUBROUTINE dust_wet

      SUBROUTINE dust_grav
!@sum  Computes dust gravitational deposition

#ifdef TRACERS_DUST
      USE constant,ONLY : visc_air,grav
      USE resolution,ONLY : Im,Jm,Lm
      USE model_com,ONLY : Dtsrc
      USE geom, ONLY : dxyp,bydxyp
      USE qusdef,ONLY : zmoms
      USE fluxes,ONLY : tr3Dsource
      USE tracer_com,ONLY : n_clay,Ntm_dust,trpdens,trm,trmom
      USE trdiag_com,ONLY : jls_grav,jls_3Dsource,taijn,
     &     tajls
#ifdef TRACERS_DRYDEP
     &     ,tij_gsdep
#endif
      USE tracers_dust_com,ONLY : zld,dradius,nDustGrav3Djl

      IMPLICIT NONE

      INTEGER :: j,l,n,n1,naij,najl
      REAL*8,DIMENSION(Ntm_dust) :: stokevdt
      REAL*8 :: stokefac1,stokefac2
      REAL*8 :: work(Im,Jm,Ntm_dust)

      IF (Lm .NE. 12)
     &     CALL stop_model
     & ('Stopped in dust_grav: Wrong vertical resolution (.ne. 12)',255)

#ifdef TRACERS_DUST_MINERAL8
      CALL stoke_mindust8(stokevdt)
#else
c     default case
      CALL stoke_sahdust(stokevdt)
#endif

      DO n=1,Ntm_dust
        n1=n_clay+n-1
        DO j=1,Jm
          stokefac1=stokevdt(n)/zld(1)*Dtsrc
          work(:,j,n)=stokefac1*trm(:,j,1,n1)*bydxyp(j)
c          trdrydep(n1,1,:,j)=trdrydep(n1,1,:,j)+
c     *         (1.-RSI(I,J))*(FOCEAN(I,J)+FLAKE(I,J))*work(:,j,n)
c          trdrydep(n1,2,:,j)=trdrydep(n1,2,:,j)+
c     *         RSI(I,J))*(FOCEAN(I,J)+FLAKE(I,J))*work(:,j,n)
c          trdrydep(n1,3,:,j)=trdrydep(n1,3,:,j)+
c     *         FLICE(I,J)*work(:,j,n)
c          trdrydep(n1,4,:,j)=trdrydep(n1,3,:,j)+
c     *         fearth(i,j)*work(:,j,n)
        END DO
        DO l=1,Lm-1
          stokefac1=stokevdt(n)/zld(l)
          stokefac2=stokevdt(n)/zld(l+1)
          tr3Dsource(:,:,l,nDustGrav3Djl,n1)=stokefac1*trm(:,:,l,n1)-
     &         stokefac2*trm(:,:,l+1,n1)
          stokefac1=stokefac1*Dtsrc
          stokefac2=stokefac2*Dtsrc
          trm(:,:,l,n1)=(1-stokefac1)*trm(:,:,l,n1)+
     &         stokefac2*trm(:,:,l+1,n1)
          trmom(zmoms,:,:,l,n1)=(1-stokefac1)*trmom(zmoms,:,:,l,n1)
        END DO
        stokefac1=stokevdt(n)/zld(Lm)
        tr3Dsource(:,:,Lm,nDustGrav3Djl,n1)=stokefac1*trm(:,:,Lm,n1)
        stokefac1=stokefac1*Dtsrc
        trm(:,:,Lm,n1)=(1-stokefac1)*trm(:,:,Lm,n1)
        trmom(zmoms,:,:,Lm,n1)=(1-stokefac1)*trmom(zmoms,:,:,Lm,n1)
      END DO

      DO n=1,Ntm_dust
        n1=n_clay+n-1
        najl=jls_grav(n1)
#ifdef TRACERS_DRYDEP
        taijn(:,:,tij_gsdep,n1)=taijn(:,:,tij_gsdep,n1) +
     &         work(:,:,n)
#endif
        DO l=1,Lm
          DO j=1,Jm
            tajls(j,l,najl)=tajls(j,l,najl)-
     &           SUM(tr3Dsource(:,j,l,nDustGrav3Djl,n1))*Dtsrc
          END DO
        END DO
      END DO
#endif

      RETURN
      END SUBROUTINE dust_grav

      SUBROUTINE stoke_sahdust(stokevdt)
!@sum  gravitational settling velocity for Sahara dust case
!@auth Jan Perlwitz, Reha Cakmur, Ina Tegen

#ifdef TRACERS_DUST
      USE constant,ONLY : Grav,Visc_air
      USE model_com,ONLY : Dtsrc
      USE tracer_com,ONLY : Ntm_dust,trpdens
      USE tracers_dust_com,ONLY : dradius

      IMPLICIT NONE

      REAL*8,INTENT(OUT) :: stokevdt(Ntm_dust)
      INTEGER :: n

      DO n=1,Ntm_dust
        stokevdt(n)=3.398D8*dradius(n)**2D0
c        stokevdt(n)=2D0*Grav*trpdens(n)*dradius(n)**2D0/(9D0*Visc_air)
      END DO
c      WRITE(*,*) 'stokevdt:',stokevdt
#endif

      RETURN
      END SUBROUTINE stoke_sahdust


      SUBROUTINE stoke_mindust8(stokevdt)
!@sum  gravitational settling velocity for 8 minerals case

#ifdef TRACERS_DUST
      USE tracer_com,ONLY : Ntm_dust

      IMPLICIT NONE

      REAL*8 :: stokevdt(Ntm_dust)
#endif

      RETURN
      END SUBROUTINE stoke_mindust8

      SUBROUTINE dust_turb
!@sum  Computes turbulent dust deposition
c*****dry deposition (turbulent mixing) 1cm/s (Giorgi (86), JGR)

#ifdef TRACERS_DUST
      USE resolution,ONLY : Jm
      USE MODEL_COM, only: dtsrc
      USE fluxes,ONLY : trsrfflx
      USE tracer_com, ONLY : n_clay,Ntm_dust,trm,trmom
      USE TRDIAG_COM, only : ijts_source,jls_source,taijs,tajls
      USE tracers_dust_com, only: nDustTurbij,nDustTurbjl

      IMPLICIT NONE

!param Thickness of first layer [m]
      REAL,PARAMETER :: Zld1=400.
      INTEGER :: j,n,n1,naij,najl
      REAL*8 :: dratio

      DO n=1,Ntm_dust
        n1=n_clay+n-1
        dratio=0.001D0*dtsrc/Zld1
        IF (dratio > 1.) dratio=1.
        trsrfflx(:,:,n1)=dratio*trm(:,:,1,n1)/Dtsrc
        trm(:,:,1,n1)=(1-dratio)*trm(:,:,1,n1)
        trmom(:,:,:,1,n1)=(1-dratio)*trmom(:,:,:,1,n1)

        naij=ijts_source(nDustTurbij,n1)
        taijs(:,:,naij)=taijs(:,:,naij)+dratio*trm(:,:,1,n1)
        najl=jls_source(nDustTurbjl,n1)
        DO j=1,Jm
          tajls(j,1,najl)=tajls(j,1,najl)+SUM(trsrfflx(:,j,n1))*Dtsrc
        END DO
        trsrfflx(:,:,n1)=0D0
      END DO
#endif

      RETURN
      END SUBROUTINE dust_turb
