#include "rundeck_opts.h"
      MODULE tracers_dust
!@sum  dust tracer parameters and variables
!@auth Reha Cakmur, Jan Perlwitz, Ina Tegen
!@ver 1.0

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)
      USE constant,ONLY : By6
      USE resolution,ONLY : Im,Jm,Lm
      USE model_com,ONLY : JMperY,JDperY
#ifndef TRACERS_AMP
      USE tracer_com,ONLY : Ntm_dust
#endif
      IMPLICIT NONE

!@param uplfac uplift factor for each size class of soil dust [kg*s**2/m**5]
#ifdef TRACERS_DUST_CUB_SAH
      REAL*8,PARAMETER :: Upclsi=52.D-9
#else
c**** default case
      REAL*8,PARAMETER :: Upclsi=12.068996D-9
#endif
#ifdef TRACERS_DUST_CUB_SAH
!@param By8 0.25d0/2d0
      REAL*8,PARAMETER :: By8=0.25D0/2D0
#else
c**** default case
      REAL*8,PARAMETER :: By4=1D0/4D0
#endif
!@param fracn fraction of uplifted soil for each size class of dust [1]
#ifdef TRACERS_DUST_CUB_SAH
      REAL*8 :: Fracl=By6,Frasi=By8
#else
c**** default case
      REAL*8 :: Fracl=0.092335D0,Frasi=0.226916D0
#endif
!@var hbaij accumulated precipitation - evaporation balance  [kg/m^2]
!@var ricntd no. of hours with negative precipitation - evaporation balance [1]
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: hbaij,ricntd
!@var dryhr  number of hours with evaporation-precipitation greater Zero
!@+          to allow dust emission
!@var frclay fraction of clay
!@var frsilt fraction of silt
!@var vtrsh  threshold wind speed above which dust emis. is allowed [m/s]
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: dryhr,frclay,frsilt,vtrsh
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
!@param Mtrac number of different fields with tracer fractions in grid box
!@+           5 clay; 5 silt
      INTEGER,PARAMETER :: Mtrac=10
!@var minf distribution of tracer fractions in grid box
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: minfr
#endif
!@var ers_data field of ERS data
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: ers_data
!@var src_fnct distribution of preferred sources
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: src_fnct
      INTEGER,PARAMETER :: Lim=220,Ljm=220,Lkm=22
!@var x1,x2,x3 indices of lock up table for emission
      REAL*8 :: x1(Lim),x2(Ljm),x3(Lkm)
!@var table array of lock up table for emission local to each grid box
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: table
!@var wsubtke_com distributed array of subscale turbulent term
!@var wsubwd_com distributed array of subscale dry convective term
!@var wsubwm_com distributed array of subscale moist convective term
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: wsubtke_com,wsubwd_com,
     &     wsubwm_com
!@var prelay distributed array with some prec info needed for simple wet
!@+          deposition scheme
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: prelay
!@var d_dust prescribed daily dust emissions [kg] (e.g. AEROCOM)
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: d_dust
#endif

      END MODULE tracers_dust

      SUBROUTINE alloc_dust(grid)
!@sum  allocates dust/mineral tracer arrays
!@auth Jan Perlwitz
!@ver  1.0

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)

      USE domain_decomp,ONLY : dist_grid
      USE resolution,ONLY : Lm
      USE model_com,ONLY : JMperY,JDperY
      USE tracer_com,ONLY : Ntm_dust
      USE tracers_dust,ONLY : hbaij,ricntd,dryhr,frclay,frsilt,vtrsh,
     &     src_fnct,ers_data,wsubtke_com,wsubwd_com,wsubwm_com,prelay,
     &     d_dust,lim,ljm,lkm,table
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
     &     ,Mtrac,minfr
#endif

      IMPLICIT NONE

      TYPE(DIST_GRID),INTENT(IN) :: grid

      INTEGER :: i_0h,i_1h,j_1h,j_0h
      INTEGER :: ier
      LOGICAL,SAVE :: qfirst=.TRUE.

      IF (.NOT. qfirst) RETURN
      qfirst=.FALSE.

      i_0h=grid%i_strt_halo
      i_1h=grid%i_stop_halo
      j_0h=grid%j_strt_halo
      j_1h=grid%j_stop_halo

      ALLOCATE(hbaij(i_0h:i_1h,j_0h:j_1h),ricntd(i_0h:i_1h,j_0h:j_1h),
     &     dryhr(i_0h:i_1h,j_0h:j_1h),frclay(i_0h:i_1h,j_0h:j_1h),
     &     frsilt(i_0h:i_1h,j_0h:j_1h),vtrsh(i_0h:i_1h,j_0h:j_1h),
     &     src_fnct(i_0h:i_1h,j_0h:j_1h),
     &     ers_data(i_0h:i_1h,j_0h:j_1h,JMperY),
     &     wsubtke_com(i_0h:i_1h,j_0h:j_1h),
     &     wsubwd_com(i_0h:i_1h,j_0h:j_1h),
     &     wsubwm_com(i_0h:i_1h,j_0h:j_1h),
     &     prelay(i_0h:i_1h,j_0h:j_1h,LM),
     &     d_dust(i_0h:i_1h,j_0h:j_1h,Ntm_dust,JDperY),
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
     &     minfr(i_0h:i_1h,j_0h:j_1h,Mtrac),
#endif
     &     STAT=ier)

      ALLOCATE(table(Lim,Ljm,Lkm),STAT=ier)

      d_dust(i_0h:i_1h,j_0h:j_1h,:,:)=0.D0

#endif

      RETURN
      END SUBROUTINE alloc_dust

      SUBROUTINE init_dust
!@sum  reads in source and parameter files for dust/mineral tracer at startup
!@auth Jan Perlwitz
!@ver  1.0
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)

      USE filemanager,ONLY: nameunit,openunit,closeunit
      USE domain_decomp,ONLY: am_i_root,dread_parallel,esmf_bcast,grid,
     &     unpack_data,write_parallel
      USE resolution, ONLY : Im,Jm
      USE model_com,ONLY : JDperY,JMperY
      USE tracer_com,ONLY : imDUST,kim,kjm,table1,x11,x21,Ntm_dust
      USE tracers_dust,ONLY : dryhr,frclay,frsilt,vtrsh,ers_data,
     &   src_fnct,table,x1,x2,x3,lim,ljm,lkm,d_dust
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
     &   ,Mtrac,minfr
#endif

      IMPLICIT NONE

      INCLUDE 'netcdf.inc'

      INTEGER :: i,ierr,j,io_data,k
      INTEGER startd(3),countd(3),statusd
      INTEGER idd1,idd2,idd3,idd4,ncidd1,ncidd2,ncidd3,ncidd4
      REAL*8 :: zsum,tabsum
c**** temporary arrays to read in global data; deallocated after read in
      REAL*4,ALLOCATABLE,DIMENSION(:,:) :: work
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: d_dust_glob

      LOGICAL,SAVE :: qfirst=.TRUE.
      CHARACTER :: cierr*3,name*256

      IF (.NOT. qfirst) RETURN
      qfirst=.FALSE.

c**** read in lookup table for calculation of mean surface wind speed from PDF
      IF (am_i_root()) THEN
        CALL openunit('LKTAB1',io_data,.TRUE.,.TRUE.)
        READ (io_data,IOSTAT=ierr) table1
        name=TRIM(nameunit(io_data))
        CALL closeunit(io_data)
      END IF
      IF (ierr == 0) THEN
        CALL write_parallel(' Read from file '//TRIM(name),unit=6)
      ELSE
        WRITE(cierr,'(I2)') ierr
        CALL write_parallel(' READ ERROR ON FILE '//TRIM(name)//
     &       ': IOSTAT='//cierr,unit=6)
        CALL stop_model('init_dust: READ ERROR',255)
      END IF
      CALL esmf_bcast(grid,table1)
c**** index of table for sub grid scale velocity (sigma) from 0.0001 to 50 m/s
      zsum=0.D0
      DO j=1,Kjm
        IF (j <= 30) THEN
          zsum=zsum+0.0001d0+FLOAT(j-1)*0.00008d0
          x21(j)=zsum
        ELSE IF (j > 30) THEN
          zsum=zsum-0.055254d0+0.005471d0*FLOAT(j)-
     &         1.938365d-4*FLOAT(j)**2.d0+
     &         3.109634d-6*FLOAT(j)**3.d0-
     &         2.126684d-8*FLOAT(j)**4.d0+
     &         5.128648d-11*FLOAT(j)**5.d0
          x21(j)=zsum
        END IF
      END DO
c**** index of table for GCM surface wind speed from 0.0001 to 50 m/s
      x11(:)=x21(:)

c**** prescribed AEROCOM dust emissions
      IF (imDust == 1) THEN

        ALLOCATE(work(Im,Jm),d_dust_glob(Im,Jm,Ntm_dust,JDperY))

        IF (am_i_root()) THEN

          statusd=NF_OPEN('dust_bin1',NCNOWRIT,ncidd1)
          statusd=NF_OPEN('dust_bin2',NCNOWRIT,ncidd2)
          statusd=NF_OPEN('dust_bin3',NCNOWRIT,ncidd3)
          statusd=NF_OPEN('dust_bin4',NCNOWRIT,ncidd4)

          statusd=NF_INQ_VARID(ncidd1,'dust',idd1)
          statusd=NF_INQ_VARID(ncidd2,'dust',idd2)
          statusd=NF_INQ_VARID(ncidd3,'dust',idd3)
          statusd=NF_INQ_VARID(ncidd4,'dust',idd4)

          startd(1)=1
          startd(2)=1

          countd(1)=Im
          countd(2)=Jm
          countd(3)=1

          DO k=1,JDperY

            IF (k > 59) THEN
              startd(3)=k+1
            ELSE
              startd(3)=k
            END IF

            statusd=NF_GET_VARA_REAL(ncidd1,idd1,startd,countd,work)
            d_dust_glob(:,:,1,k)=DBLE(work(:,:))
            statusd=NF_GET_VARA_REAL(ncidd2,idd2,startd,countd,work)
            d_dust_glob(:,:,2,k)=DBLE(work(:,:))
            statusd=NF_GET_VARA_REAL(ncidd3,idd3,startd,countd,work)
            d_dust_glob(:,:,3,k)=DBLE(work(:,:))
            statusd=NF_GET_VARA_REAL(ncidd4,idd4,startd,countd,work)
            d_dust_glob(:,:,4,k)=DBLE(work(:,:))

          END DO

          statusd=NF_CLOSE('dust_bin1',NCNOWRIT,ncidd1)
          statusd=NF_CLOSE('dust_bin2',NCNOWRIT,ncidd2)
          statusd=NF_CLOSE('dust_bin3',NCNOWRIT,ncidd3)
          statusd=NF_CLOSE('dust_bin4',NCNOWRIT,ncidd4)

        END IF                  ! end am_i_root

        CALL write_parallel(' Read from file dust_bin[1-4]',unit=6)
        CALL unpack_data(grid,d_dust_glob,d_dust)

        DEALLOCATE(work,d_dust_glob)

      ELSE IF (imDUST == 0) THEN
c**** interactive dust emissions

c**** Read input: threshold speed
        CALL openunit('VTRSH',io_data,.TRUE.,.TRUE.)
        CALL dread_parallel(grid,io_data,nameunit(io_data),vtrsh)
        CALL closeunit(io_data)

c**** Read input: fraction clay
        CALL openunit('FRCLAY',io_data,.TRUE.,.TRUE.)
        CALL dread_parallel(grid,io_data,nameunit(io_data),frclay)
        CALL closeunit(io_data)
        
c**** Read input: fraction silt
        CALL openunit('FRSILT',io_data,.TRUE.,.TRUE.)
        CALL dread_parallel(grid,io_data,nameunit(io_data),frsilt)
        CALL closeunit(io_data)

c**** Read input: prec-evap data
        CALL openunit('DRYHR',io_data,.TRUE.,.TRUE.)
        CALL dread_parallel(grid,io_data,nameunit(io_data),dryhr)
        CALL closeunit(io_data)

c**** Read input: ERS data
        CALL openunit('ERS',io_data,.TRUE.,.TRUE.)
        DO k=1,JMperY
          CALL dread_parallel(grid,io_data,nameunit(io_data),
     &         ers_data(:,:,k))
        END DO
        CALL closeunit(io_data)

c**** Read input: source function data
        CALL openunit('GIN',io_data,.TRUE.,.TRUE.)
        CALL dread_parallel(grid,io_data,nameunit(io_data),src_fnct)
        CALL closeunit(io_data)

c**** Read input: EMISSION LOOKUP TABLE data
        IF (am_i_root()) THEN
          CALL openunit('LKTAB',io_data,.TRUE.,.TRUE.)
          DO k=1,Lkm
            READ(io_data,IOSTAT=ierr) ((table(i,j,k),i=1,Lim),j=1,Ljm)
          END DO
          name=nameunit(io_data)
          CALL closeunit(io_data)
        END IF
        IF (ierr == 0) THEN
          CALL write_parallel(' Read from file '//TRIM(name),unit=6)
        ELSE
          WRITE(cierr,'(I2)') ierr
          CALL write_parallel(' READ ERROR ON FILE '//TRIM(name)//
     &         ': IOSTAT='//cierr,unit=6)
          CALL stop_model('init_dust: READ ERROR',255)
        END IF
        CALL esmf_bcast(grid,table)
          
c**** index of table for threshold velocity from 6.5 to 17 m/s
        DO k=1,Lkm
          x3(k)=6.d0+0.5d0*k
        END DO

c**** index of table for sub grid scale velocity (sigma) from .0001 to 30 m/s
        zsum=0.d0
        DO j=1,Ljm
          IF (j <= 30) THEN
            zsum=zsum+0.0001d0+FLOAT(j-1)*0.00008d0
            x2(j)=zsum
          ELSE IF (j > 30) THEN
            zsum=zsum-0.055254d0+0.005471d0*FLOAT(j)-
     &           1.938365d-4*FLOAT(j)**2.d0+
     &           3.109634d-6*FLOAT(j)**3.d0-
     &           2.126684d-8*FLOAT(j)**4.d0+
     &           5.128648d-11*FLOAT(j)**5.d0
            x2(j)=zsum
          END IF
        END DO
c**** index of table for GCM surface wind speed from 0.0001 to 30 m/s
        x1(:)=x2(:)
      ELSE
        CALL stop_model
     &     ('Stopped in tracer_IC: parameter imDUST must be 0 or 1',255)
      END IF

#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
        CALL openunit('MINFR',io_data,.TRUE.,.TRUE.)
        CALL dread_parallel(grid,io_data,nameunit(io_data),minfr)
        CALL closeunit(io_data)
#endif

#endif
      RETURN
      END SUBROUTINE init_dust
