#include "rundeck_opts.h"
      SUBROUTINE tracers_dust
!@sum soil dust sources and sinks
!auth Reha Cakmur, Jan Perlwitz, Ina Tegen
      USE resolution,ONLY : Im,Jm
      USE tracers_dust_com,ONLY : dryhr,frclay,frsilt,vtrsh,hbaij
      USE filemanager,ONLY : openunit,closeunit

      IMPLICIT NONE

      INTEGER :: i,j
      INTEGER :: io_data
      LOGICAL :: ifirst=.true.

      if (ifirst) then
        hbaij=0D0
c**** Read input: threshold speed
        call openunit('VTRSH',io_data,.true.,.true.)
        read (io_data) vtrsh
        call closeunit(io_data)
c**** Read input: fraction clay
        call openunit('FRCLAY',io_data,.true.,.true.)
        read (io_data) frclay
        call closeunit(io_data)
c**** Read input: fraction silt
        call openunit('FRSILT',io_data,.true.,.true.)
        read (io_data) frsilt
        call closeunit(io_data)
c**** Read input: prec-evap data
        call openunit('DRYHR',io_data,.true.,.true.)
        read (io_data) dryhr
        call closeunit(io_data)
#ifdef TRACERS_DUST_MINERAL8
        CALL openunit('MINFR',io_data,.true.,.true.)
        CALL closeunit(io_data)
#endif
        ifirst=.false.
      endif

c      do j=1,jm
c         do i=1,im
c            if (bldata(i,j,1) > vtrsh(i,j)tadd(i,j,1)=tadd(i,j,1)+1.
c            if (frclay(i,j)*pmei(i,j) > 0)tadd(i,j,3)=tadd(i,j,3)+1.
c         enddo
c      enddo

      CALL dust_emission
      CALL dust_grav
      CALL dust_wet
      CALL dust_turb

      RETURN
      END SUBROUTINE tracers_dust

      SUBROUTINE dust_emission
!@sum  dust source flux
!@auth Jan Perlwitz, Reha Cakmur, Ina Tegen
      USE resolution, ONLY : Im,Jm
      USE model_com,ONLY : Dtsrc,fearth
      USE qusdef,ONLY : mz
      USE fluxes,ONLY : prec,evapor,trsrfflx
      USE tracer_com,ONLY : n_clay,Ntm_dust,trm,trmom
      USE tracer_diag_com,ONLY : ijts_source,jls_source,taijs,tajls
      USE geom,ONLY : dxyp
      USE ghycom,ONLY : snowe !earth snow amount
      USE pblcom,ONLY : wsavg
      USE tracers_dust_com,ONLY : dryhr,frclay,frsilt,hbaij,ricntd,
     &     vtrsh,nDustEmij,nDustEmjl

      IMPLICIT NONE

      INTEGER :: i,j,n,n1,naij,najl
      REAL*8 :: hbaijold,hbaijd,workijn
      REAL*8 :: dsrcflx(Ntm_dust)
c     dsrcflx  dust source flux for Ntm_dust tracers [kg/s]
      LOGICAL :: pmei(Im,Jm)

c     Checking whether accumulated precipitation - evaporation
c     less/equal than Zero for a succeeding number of hours greater/equal
c     than threshold dryhr to permit dust emission
      DO j=1,Jm
        DO i=1,Im
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

c     Loop for calculating dust source flux for each tracer
      DO n=1,Ntm_dust
        n1=n_clay+n-1
        trsrfflx(:,:,n1)=0D0
      END DO

      DO j=1,Jm
        DO i=1,Im

          IF (fearth(i,j) > 0. .AND. snowe(i,j) <= 1 .AND.
     &         vtrsh(i,j) > 0. .AND. wsavg(i,j) > vtrsh(i,j) .AND.
     &         pmei(i,j)) THEN

#ifdef TRACERS_DUST_MINERAL8
#ifdef TRACERS_DUST_TURB
            CALL loc_dsrcflx_turb_min8
#else
            CALL loc_dsrcflx_cub_min8
#endif
#else
#ifdef TRACERS_DUST_TURB
            CALL loc_dsrcflx_turb_sah
#else
c     default case
            CALL loc_dsrcflx_cub_sah(dxyp(j),wsavg(i,j),vtrsh(i,j),
     &           frclay(i,j),frsilt(i,j),fearth(i,j),dsrcflx)
#endif
#endif

            DO n=1,Ntm_dust
              n1=n_clay+n-1
              trsrfflx(i,j,n1)=dsrcflx(n)
            END DO

          END IF

        END DO
      END DO

      DO n=1,Ntm_dust
        n1=n_clay+n-1
        naij=ijts_source(nDustEmij,n1)
        najl=jls_source(nDustEmjl,n1)
c        WRITE(*,*) 'naij,najl:',naij,najl
        DO j=1,Jm
          DO i=1,Im
            IF (trsrfflx(i,j,n1) > 0.) THEN
              workijn=trsrfflx(i,j,n1)*Dtsrc
              trm(i,j,1,n1)=trm(i,j,1,n1)+workijn
              trmom(mz,i,j,1,n1)=trmom(mz,i,j,1,n1)-workijn

              taijs(i,j,naij)=taijs(i,j,naij)+workijn
            END IF
          END DO
          tajls(j,1,najl)=tajls(j,1,najl)+SUM(trsrfflx(:,j,n1))*Dtsrc
        END DO
      END DO
            
      RETURN
      END SUBROUTINE dust_emission

      SUBROUTINE loc_dsrcflx_cub_sah(dxypj,wsavgij,vtrshij,frclayij,
     &     frsiltij,fearthij,dsrcflx)
!@sum  local dust source flux physics according to Ina's old cubic scheme
!@auth Ina Tegen, Jan Perlwitz, Reha Cakmur

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

      RETURN
      END SUBROUTINE loc_dsrcflx_cub_sah

      SUBROUTINE loc_dsrcflx_cub_min8
!@sum  local dust source flux physics with Ina's cubic scheme and 8 minerals
!@auth Jan Perlwitz, ...

      RETURN
      END SUBROUTINE loc_dsrcflx_cub_min8

      SUBROUTINE loc_dsrcflx_turb_sah
!@sum  local dust source flux physics with turbulent fluxes and Sahara dust
!@auth Reha Cakmur, ...

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
      USE resolution,ONLY : Im,Jm,Lm
      USE model_com,ONLY : Dtsrc
      USE qusdef,ONLY : mx,my,mz
      USE fluxes,ONLY : prec,trsrfflx,tr3Dsource
      USE tracer_com,ONLY : n_clay,Ntm_dust,trm,trmom
      USE tracer_diag_com,ONLY : ijts_source,jls_3Dsource,taijs,tajls
      USE tracers_dust_com,ONLY : Z,nDustWetij,nDustWet3Djl

      IMPLICIT NONE

      INTEGER :: i,j,l,n,n1,naij,najl
      INTEGER,DIMENSION(Jm) :: lwdep
      REAL*8,DIMENSION(jm) :: h
      REAL*8 :: y

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
              trmom(mx,i,j,l,n1)=trmom(mx,i,j,l,n1)*(1.-y)
              trmom(my,i,j,l,n1)=trmom(my,i,j,l,n1)*(1.-y)
              trmom(mz,i,j,l,n1)=trmom(mz,i,j,l,n1)*(1.-y)
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

      END SUBROUTINE dust_wet

      SUBROUTINE dust_grav
!@sum  Computes dust gravitational deposition
      USE CONSTANT,ONLY : visc_air,grav
      USE resolution,ONLY : Jm,Lm
      USE MODEL_COM,ONLY: Dtsrc
      USE qusdef,ONLY : mz
      USE fluxes,ONLY : tr3Dsource,trgrdep
      USE TRACER_COM,ONLY : n_clay,Ntm_dust,trpdens,trm,trmom
      USE TRACER_DIAG_COM, only : ijts_source,jls_3Dsource,taijs,tajls
      USE tracers_dust_com, only: zld,dradius,nDustGravij,nDustGrav3Djl

      IMPLICIT NONE

      INTEGER :: j,l,n,n1,naij,najl
      REAL*8,DIMENSION(Ntm_dust) :: stokevdt
      REAL*8 :: stokefac1,stokefac2

#ifdef TRACERS_DUST_MINERAL8
      CALL stoke_mindust8(stokevdt)
#else
c     default case
      CALL stoke_sahdust(stokevdt)
#endif

      DO n=1,Ntm_dust
        n1=n_clay+n-1
        stokefac1=stokevdt(n)/zld(1)*Dtsrc
        trgrdep(n1,:,:)=stokefac1*trm(:,:,1,n1)
c        WRITE(*,*) 'n1,stokefac1,trgrdep(n1,:,:):',n1,stokefac1,
c     &       trgrdep(n1,:,:)
        DO l=1,Lm-1
          stokefac1=stokevdt(n)/zld(l)
          stokefac2=stokevdt(n)/zld(l+1)
          tr3Dsource(:,:,l,nDustGrav3Djl,n1)=stokefac1*trm(:,:,l,n1)-
     &         stokefac2*trm(:,:,l+1,n1)
          stokefac1=stokefac1*Dtsrc
          stokefac2=stokefac2*Dtsrc
          trm(:,:,l,n1)=(1-stokefac1)*trm(:,:,l,n1)+
     &         stokefac2*trm(:,:,l+1,n1)
          trmom(mz,:,:,l,n1)=(1-stokefac1)*trmom(mz,:,:,l,n1)+stokefac2*
     &         trmom(mz,:,:,l+1,n1)
        END DO
        stokefac1=stokevdt(n)/zld(Lm)
        tr3Dsource(:,:,Lm,nDustGrav3Djl,n1)=stokefac1*trm(:,:,Lm,n1)
        stokefac1=stokefac1*Dtsrc
        trm(:,:,Lm,n1)=(1-stokefac1)*trm(:,:,Lm,n1)
        trmom(mz,:,:,Lm,n1)=(1-stokefac1)*trmom(mz,:,:,Lm,n1)
      END DO

      DO n=1,Ntm_dust
        n1=n_clay+n-1
        naij=ijts_source(nDustGravij,n1)
        najl=jls_3Dsource(nDustGrav3Djl,n1)
c        WRITE(*,*) 'naij,najl:',naij,najl
        taijs(:,:,naij)=taijs(:,:,naij)+trgrdep(n1,:,:)
c        WRITE(*,*) 'n1,taijs(:,:,naij):',n1,taijs(:,:,naij)
        DO l=1,Lm
          DO j=1,Jm
            tajls(j,l,najl)=tajls(j,l,najl)+
     &           SUM(tr3Dsource(:,j,l,nDustGrav3Djl,n1))*Dtsrc
          END DO
        END DO
c        WRITE(*,*) 'n1,tajls(:,:,najl):',n1,tajls(:,:,najl)
      END DO

      RETURN
      END SUBROUTINE dust_grav

      SUBROUTINE stoke_sahdust(stokevdt)
!@sum  gravitational settling velocity for Sahara dust case
!@auth Jan Perlwitz, Reha Cakmur, Ina Tegen
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

      RETURN
      END SUBROUTINE stoke_sahdust


      SUBROUTINE stoke_mindust8(stokevdt)
!@sum  gravitational settling velocity for 8 minerals case
      USE tracer_com,ONLY : Ntm_dust

      IMPLICIT NONE

      REAL*8 :: stokevdt(Ntm_dust)

      RETURN
      END SUBROUTINE stoke_mindust8

      SUBROUTINE dust_turb
!@sum  Computes turbulent dust deposition
c*****dry deposition (turbulent mixing) 1cm/s (Giorgi (86), JGR)
      USE resolution,ONLY : Jm
      USE MODEL_COM, only: dtsrc
      USE qusdef,ONLY : mz
      USE fluxes,ONLY : trsrfflx
      USE tracer_com, ONLY : n_clay,Ntm_dust,trm,trmom
      USE TRACER_DIAG_COM, only : ijts_source,jls_source,taijs,tajls
      USE tracers_dust_com, only: zld,nDustTurbij,nDustTurbjl

      IMPLICIT NONE

      INTEGER :: j,n,n1,naij,najl
      REAL*8 :: dratio

      DO n=1,Ntm_dust
        n1=n_clay+n-1
        dratio=0.001D0*dtsrc/zld(1)
        IF (dratio > 1.) dratio=1.
        trsrfflx(:,:,n1)=dratio*trm(:,:,1,n1)/Dtsrc
        trm(:,:,1,n1)=(1-dratio)*trm(:,:,1,n1)
        trmom(mz,:,:,1,n1)=(1-dratio)*trmom(mz,:,:,1,n1)

        naij=ijts_source(nDustTurbij,n1)
        taijs(:,:,naij)=taijs(:,:,naij)+dratio*trm(:,:,1,n1)
        najl=jls_source(nDustTurbjl,n1)
        DO j=1,Jm
          tajls(j,1,najl)=tajls(j,1,najl)+SUM(trsrfflx(:,j,n1))*Dtsrc
        END DO
        trsrfflx(:,:,n1)=0D0
      END DO

      RETURN
      END SUBROUTINE dust_turb
