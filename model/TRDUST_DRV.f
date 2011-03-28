#include "rundeck_opts.h"
      module trdust_drv
!@sum  trdust_drv routines with resolution dependent variables for
!@+               soil dust aerosols
!@auth Jan Perlwitz
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)

      use filemanager,only: nameunit,openunit,closeunit
      use constant, only: rgas
      use resolution, only: im,jm,lm
      use domain_decomp_atm, only: am_i_root,grid,dread_parallel
     &     ,esmf_bcast,write_parallel,get
      use model_com, only: ioread,iowrite,irsfic,irsficno
     &     ,irerun,JDperY,JMperY,itime
      use fluxes, only: dust_flux_glob
#if (defined TRACERS_DUST) || (defined TRACERS_AMP)
     &     ,dust_flux2_glob
#endif
#ifdef TRACERS_DRYDEP
     &     ,depo_turb_glob,depo_grav_glob
#endif
#ifdef TRACERS_WATER
     &     ,trprec
#else
     &     ,trprec_dust
#endif
      use tracer_com,only: n_clay,n_clayilli,n_sil1quhe,ntm_dust,trm
     &     ,trname,coupled_chem
#ifdef TRACERS_DRYDEP
     &     ,dodrydep
#endif
#ifdef TRACERS_WATER
     &     ,dowetdep
#endif
      use trdiag_com, only: trcsurf,trcSurfByVol
      use tracers_dust
#ifdef NEW_IO
      use pario, only: defvar,read_dist_data,write_dist_data
#endif

      implicit none

      integer :: i_0,i_1,j_0,j_1

#endif /*TRACERS_DUST || TRACERS_MINERALS || TRACERS_QUARZHEM || TRACERS_AMP*/

      contains

c init_dust
      subroutine init_dust
!@sum  init_dust reads in source and parameter files for dust/mineral tracer at startup
!@auth Jan Perlwitz
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP)

      IMPLICIT NONE

      INCLUDE 'netcdf.inc'

      integer :: i,ierr,j,io_data,k,ires,n,n1
      INTEGER startd(3),countd(3),statusd
      INTEGER idd1,idd2,idd3,idd4,ncidd1,ncidd2,ncidd3,ncidd4
      REAL*8 :: zsum,tabsum
c**** temporary array to read in data
      REAL*4,DIMENSION(grid%i_strt:grid%i_stop,
     &                 grid%j_strt:grid%j_stop) :: work ! no halo

      LOGICAL,SAVE :: qfirst=.TRUE.
      CHARACTER :: cierr*3,name*256
      CHARACTER(50) :: OptModelVers='No Entry'

      IF (.NOT. qfirst) RETURN
      qfirst=.FALSE.

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1, I_STRT=I_0, I_STOP=I_1)

#ifdef TRACERS_DUST
      n_soilDust = n_clay
#else
#ifdef TRACERS_MINERALS
      n_soilDust = n_clayilli
#else
#ifdef TRACERS_QUARZHEM
      n_soilDust = n_sil1quhe
#endif
#endif
#endif

#ifndef TRACERS_AMP
c**** initialize dust names
      do n=1,Ntm_dust
        n1=n_soilDust+n-1
        dust_names(n) = trname(n1)
      end do
#endif

c**** read in lookup table for calculation of mean surface wind speed from PDF
      IF (am_i_root()) THEN
        CALL openunit('LKTAB1',io_data,.TRUE.,.TRUE.)
        READ (io_data,IOSTAT=ierr) table1
        name=TRIM(nameunit(io_data))
        CALL closeunit(io_data)
        IF (ierr == 0) then
          write(6,*) ' Read from file '//TRIM(name)
        ELSE
          WRITE(cierr,'(I2)') ierr
          write(6,*) ' READ ERROR ON FILE '//TRIM(name)//' rc='//cierr
        END IF
      END IF
      CALL esmf_bcast(grid,ierr)
      IF (ierr.ne.0) CALL stop_model('init_dust: READ ERROR',255)
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

c these netcdf reads are still latlon-specific.
c will call read_dist_data for cubed sphere compatibility

        statusd=NF_OPEN('dust_bin1',NCNOWRIT,ncidd1)
        statusd=NF_OPEN('dust_bin2',NCNOWRIT,ncidd2)
        statusd=NF_OPEN('dust_bin3',NCNOWRIT,ncidd3)
        statusd=NF_OPEN('dust_bin4',NCNOWRIT,ncidd4)

        statusd=NF_INQ_VARID(ncidd1,'dust',idd1)
        statusd=NF_INQ_VARID(ncidd2,'dust',idd2)
        statusd=NF_INQ_VARID(ncidd3,'dust',idd3)
        statusd=NF_INQ_VARID(ncidd4,'dust',idd4)

        startd(1)=i_0
        startd(2)=j_0

        countd(1)=1+(i_1-i_0)
        countd(2)=1+(j_1-j_0)
        countd(3)=1

        DO k=1,JDperY

          IF (k > 59) THEN
            startd(3)=k+1
          ELSE
            startd(3)=k
          END IF

          statusd=NF_GET_VARA_REAL(ncidd1,idd1,startd,countd,work)
          d_dust(i_0:i_1,j_0:j_1,1,k)=DBLE(work(i_0:i_1,j_0:j_1))
          statusd=NF_GET_VARA_REAL(ncidd2,idd2,startd,countd,work)
          d_dust(i_0:i_1,j_0:j_1,2,k)=DBLE(work(i_0:i_1,j_0:j_1))
          statusd=NF_GET_VARA_REAL(ncidd3,idd3,startd,countd,work)
          d_dust(i_0:i_1,j_0:j_1,3,k)=DBLE(work(i_0:i_1,j_0:j_1))
          statusd=NF_GET_VARA_REAL(ncidd4,idd4,startd,countd,work)
          d_dust(i_0:i_1,j_0:j_1,4,k)=DBLE(work(i_0:i_1,j_0:j_1))

        END DO

        statusd=NF_CLOSE(ncidd1)
        statusd=NF_CLOSE(ncidd2)
        statusd=NF_CLOSE(ncidd3)
        statusd=NF_CLOSE(ncidd4)

        CALL write_parallel(' Read from file dust_bin[1-4]',unit=6)

      ELSE IF (imDust == 2) THEN
c**** legacy emission scheme
        IF (Im /= 72 .OR. Jm /= 46) CALL stop_model
     & ('Stopped in init_dust: imDust=2 works only for Im=72 and Jm=46'
     &       ,255)

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

      ELSE IF (imDust == 0) THEN
c**** Probability density function scheme for dust emission

c**** Read input: ERS data
        CALL openunit('ERS',io_data,.TRUE.,.TRUE.)
        DO k=1,JMperY
          CALL dread_parallel(grid,io_data,nameunit(io_data),
     &         ers_data(:,:,k))
        END DO
        CALL closeunit(io_data)

c**** Read input: source function data
        call openunit('DSRC',io_data,.true.,.true.)
        call dread_parallel(grid,io_data,nameunit(io_data)
     &       ,dustSourceFunction)
        CALL closeunit(io_data)

c**** set parameters depending on the preferred sources chosen
        select case(prefDustSources)
        case(0)                 ! Ginoux 2001 sources w/ vegetation mask
          select case(im)
          case(72)              ! uses old values for Ginoux 2001 source file
            fracClayPDFscheme = 0.092335D0 ! not optimized yet
            fracSiltPDFscheme = 0.226916D0 ! not optimized yet
            ires=1
          case(144)
            if (coupled_chem == 1) then
#ifdef NUDGE_ON
              fracClayPDFscheme = 0.036636533D0
              fracSiltPDFscheme = 0.023985135D0
#else
              fracClayPDFscheme = 0.091393298D0
              fracSiltPDFscheme = 0.10973922D0
#endif
              OptModelVers='02/20/2010, 23:30 EST'
            else
#ifdef NUDGE_ON
              fracClayPDFscheme = 0.039218312D0
              fracSiltPDFscheme = 0.021643261D0
#else
              fracClayPDFscheme = 0.10309873D0
              fracSiltPDFscheme = 0.080987688D0
#endif
              OptModelVers='02/20/2010, 23:30 EST'
            end if
            ires=2
          case(288)             ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.091393298D0  ! not optimized yet
              fracSiltPDFscheme = 0.10973922D0   ! not optimized yet
            else
              fracClayPDFscheme = 0.10309873D0   ! not optimized yet
              fracSiltPDFscheme = 0.080987688D0  ! not optimized yet
            end if
            ires=3
          case(360)             ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.091393298D0  ! not optimized yet
              fracSiltPDFscheme = 0.10973922D0   ! not optimized yet
            else
              fracClayPDFscheme = 0.10309873D0   ! not optimized yet
              fracSiltPDFscheme = 0.080987688D0  ! not optimized yet
            end if
            ires=4
          case(90)              ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.091393298D0  ! not optimized yet
              fracSiltPDFscheme = 0.10973922D0   ! not optimized yet
            else
              fracClayPDFscheme = 0.10309873D0   ! not optimized yet
              fracSiltPDFscheme = 0.080987688D0  ! not optimized yet
            end if
            ires=5
          end select
        case(1)                 ! Ginoux 2009 sources w/ vegetation mask
          select case(im)
          case(72)              ! uses old values for Ginoux 2001 source file
            fracClayPDFscheme = 0.092335D0 ! not optimized yet
            fracSiltPDFscheme = 0.226916D0 ! not optimized yet
            ires=1
          case(144)
            if (coupled_chem == 1) then
#ifdef NUDGE_ON
              fracClayPDFscheme = 0.034805282D0
              fracSiltPDFscheme = 0.0307729D0
#else
              fracClayPDFscheme = 0.086874723D0
              fracSiltPDFscheme = 0.097082074D0
#endif
              OptModelVers='02/20/2010, 23:30 EST'
            else
#ifdef NUDGE_ON
              fracClayPDFscheme = 0.038017106D0
              fracSiltPDFscheme = 0.029505731D0
#else
              fracClayPDFscheme = 0.091387274D0
              fracSiltPDFscheme = 0.15582714D0
#endif
              OptModelVers='02/20/2010, 23:30 EST'
            end if
            ires=2
          case(288)             ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.086874723D0 ! not optimized yet
              fracSiltPDFscheme = 0.097082074D0 ! not optimized yet
            else
              fracClayPDFscheme = 0.091387274D0 ! not optimized yet
              fracSiltPDFscheme = 0.15582714D0  ! not optimized yet
            end if
            ires=3
          case(360)             ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.086874723D0 ! not optimized yet
              fracSiltPDFscheme = 0.097082074D0 ! not optimized yet
            else
              fracClayPDFscheme = 0.091387274D0 ! not optimized yet
              fracSiltPDFscheme = 0.15582714D0  ! not optimized yet
            end if
            ires=4
          case(90)              ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.086874723D0 ! not optimized yet
              fracSiltPDFscheme = 0.097082074D0 ! not optimized yet
            else
              fracClayPDFscheme = 0.091387274D0 ! not optimized yet
              fracSiltPDFscheme = 0.15582714D0  ! not optimized yet
            end if
            ires=5
          end select
        case(2)                 ! Ginoux 2009 sources w/o vegetation mask
          select case(im)
          case(72)              ! uses old values for Ginoux 2001 source file
            fracClayPDFscheme = 0.092335D0 ! not optimized yet
            fracSiltPDFscheme = 0.226916D0 ! not optimized yet
            ires=1
          case(144)
            if (coupled_chem == 1) then
#ifdef NUDGE_ON
              fracClayPDFscheme = 0.028216341D0
              fracSiltPDFscheme = 0.010625339D0
#else
              fracClayPDFscheme = 0.051027254D0
              fracSiltPDFscheme = 0.047415049D0
#endif
              OptModelVers='02/20/2010, 23:30 EST'
            else
#ifdef NUDGE_ON
              fracClayPDFscheme = 0.03028031D0
              fracSiltPDFscheme = 0.097915021D-1
#else
              fracClayPDFscheme = 0.059176377D0
              fracSiltPDFscheme = 0.056770932D0
#endif
              OptModelVers='02/20/2010, 23:30 EST'
            end if
            ires=2
          case(288)             ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.051027254D0 ! not optimized yet
              fracSiltPDFscheme = 0.047415049D0 ! not optimized yet
            else
              fracClayPDFscheme = 0.059176377D0 ! not optimized yet
              fracSiltPDFscheme = 0.056770932D0 ! not optimized yet
            end if
            ires=3
          case(360)             ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.051027254D0 ! not optimized yet
              fracSiltPDFscheme = 0.047415049D0 ! not optimized yet
            else
              fracClayPDFscheme = 0.059176377D0 ! not optimized yet
              fracSiltPDFscheme = 0.056770932D0 ! not optimized yet
            end if
            ires=4
          case(90)              ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.051027254D0 ! not optimized yet
              fracSiltPDFscheme = 0.047415049D0 ! not optimized yet
            else
              fracClayPDFscheme = 0.059176377D0 ! not optimized yet
              fracSiltPDFscheme = 0.056770932D0 ! not optimized yet
            end if
            ires=5
          end select
        case(3)                 ! Grini/Zender sources
          select case(im)
          case(72)              ! uses old values for Ginoux 2001 source file
            fracClayPDFscheme = 0.092335D0 ! not optimized yet
            fracSiltPDFscheme = 0.226916D0 ! not optimized yet
            ires=1
          case(144)
            if (coupled_chem == 1) then
#ifdef NUDGE_ON
              fracClayPDFscheme = 0.015809495D0
              fracSiltPDFscheme = 0.059532525D-1
#else
              fracClayPDFscheme = 0.03586315D0
              fracSiltPDFscheme = 0.028140008D0
#endif
              OptModelVers='02/20/2010, 23:30 EST'
            else
#ifdef NUDGE_ON
              fracClayPDFscheme = 0.016681537D0
              fracSiltPDFscheme = 0.054484947D-1
#else
              fracClayPDFscheme = 0.036893354D0
              fracSiltPDFscheme = 0.039649002D0
#endif
              OptModelVers='02/20/2010, 23:30 EST'
            end if
            ires=2
          case(288)             ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.03586315D0  ! not optimized yet
              fracSiltPDFscheme = 0.028140008D0 ! not optimized yet
            else
              fracClayPDFscheme = 0.036893354D0 ! not optimized yet
              fracSiltPDFscheme = 0.039649002D0 ! not optimized yet
            end if
            ires=3
          case(360)             ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.03586315D0  ! not optimized yet
              fracSiltPDFscheme = 0.028140008D0 ! not optimized yet
            else
              fracClayPDFscheme = 0.036893354D0 ! not optimized yet
              fracSiltPDFscheme = 0.039649002D0 ! not optimized yet
            end if
            ires=4
          case(90)              ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.03586315D0  ! not optimized yet
              fracSiltPDFscheme = 0.028140008D0 ! not optimized yet
            else
              fracClayPDFscheme = 0.036893354D0 ! not optimized yet
              fracSiltPDFscheme = 0.039649002D0 ! not optimized yet
            end if
            ires=5
          end select
        case(4)                 ! Tegen sources
          select case(im)
          case(72)              ! uses old values for Ginoux 2001 source file
            fracClayPDFscheme = 0.092335D0 ! not optimized yet
            fracSiltPDFscheme = 0.226916D0 ! not optimized yet
            ires=1
          case(144)
            if (coupled_chem == 1) then
#ifdef NUDGE_ON
              fracClayPDFscheme = 0.046381759D0
              fracSiltPDFscheme = 0.040315694D-1
#else
              fracClayPDFscheme = 0.083940461D0
              fracSiltPDFscheme = 0.06811547D0
#endif
              OptModelVers='02/20/2010, 23:30 EST'
            else
#ifdef NUDGE_ON
              fracClayPDFscheme = 0.048860416D0
              fracSiltPDFscheme = 0.029009232D-1
#else
              fracClayPDFscheme = 0.09620785D0
              fracSiltPDFscheme = 0.058139336D0
#endif
              OptModelVers='02/20/2010, 23:30 EST'
            end if
            ires=2
          case(288)             ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.083940461D0 ! not optimized yet
              fracSiltPDFscheme = 0.06811547D0  ! not optimized yet
            else
              fracClayPDFscheme = 0.09620785D0  ! not optimized yet
              fracSiltPDFscheme = 0.058139336D0 ! not optimized yet
            end if
            ires=3
          case(360)             ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.083940461D0 ! not optimized yet
              fracSiltPDFscheme = 0.06811547D0  ! not optimized yet
            else
              fracClayPDFscheme = 0.09620785D0  ! not optimized yet
              fracSiltPDFscheme = 0.058139336D0 ! not optimized yet
            end if
            ires=4
          case(90)              ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.083940461D0 ! not optimized yet
              fracSiltPDFscheme = 0.06811547D0  ! not optimized yet
            else
              fracClayPDFscheme = 0.09620785D0  ! not optimized yet
              fracSiltPDFscheme = 0.058139336D0 ! not optimized yet
            end if
            ires=5
          end select
        end select

        if (am_i_root()) then
          write(6,*) ' Actually used parameters for soil dust emission:'
          write(6,'(1x,a28,f12.9)') '  Clay: fracClayPDFscheme = '
     &         ,fracClayPDFscheme
          write(6,'(1x,a28,f12.9)') '  Silt: fracSiltPDFscheme = '
     &         ,fracSiltPDFscheme
          if (prefDustSources <= numDustSourceOpt-1) then
            write(6,*) '  For prefDustSources =',prefDustSources,','
            write(6,*) '  these parameters are the optimized values for'
            write(6,*) '  following file with preferred dust sources:'
            write(6,*) '  >> '
     &           ,trim(dustSourceFile(prefDustSources,ires)),' <<'
            write(6,*) '  for model version from ',trim(OptModelVers)
            if (coupled_chem == 1) then
              write(6,*) '  for all aerosols and chemistry'
            else
              write(6,*) '  for dust aerosols w/o other aerosols and'
     &             ,' chemistry'
            end if
            write(6,*) '  For a free choice of emission parameters set'
            write(6,*) '  prefDustSources >',numDustSourceOpt-1
            write(6,*) '  and set fracClayPDFscheme and'
     &           ,' fracSiltPDFscheme in rundeck'
          end if
        end if

c**** Read input: EMISSION LOOKUP TABLE data
        IF (am_i_root()) THEN
          CALL openunit('LKTAB',io_data,.TRUE.,.TRUE.)
          DO k=1,Lkm
            READ(io_data,IOSTAT=ierr) ((table(i,j,k),i=1,Lim),j=1,Ljm)
          END DO
          name=nameunit(io_data)
          CALL closeunit(io_data)
          IF (ierr == 0) THEN
            write(6,*) ' Read from file '//TRIM(name)
          ELSE
            WRITE(cierr,'(I2)') ierr
            write(6,*) ' READ ERROR ON FILE '//TRIM(name)//' rc='//cierr
          END IF
        END IF
        CALL esmf_bcast(grid,ierr)
        if(ierr.ne.0) CALL stop_model('init_dust: READ ERROR',255)
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
     &     ('Stopped in init_dust: parameter imDUST must be <= 2',255)
      END IF

#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
        CALL openunit('MINFR',io_data,.TRUE.,.TRUE.)
        CALL dread_parallel(grid,io_data,nameunit(io_data),minfr)
        CALL closeunit(io_data)
#endif

#endif /*TRACERS_DUST || TRACERS_MINERALS || TRACERS_QUARZHEM*/
      RETURN
      END SUBROUTINE init_dust

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
c io_trDust
      subroutine io_trDust(kunit,iaction)
!@sum  io_trDust I/O of specific dust aerosol diagnostics (not netcdf)
!@auth Jan Perlwitz
      use domain_decomp_1d, only: pack_data,unpack_data
      implicit none

!@var kunit unit number of read/write
!@var iaction flag for reading or writing to file
      integer,intent(in) :: iaction,kunit

!@var iostat I/O status
      integer :: iostat
!@var header text with information on variables
      character(len=80) :: header

!@var dustDiagSubdd_glob global structured variable of dust aerosol diagnostics
      type(dustDiagSubdd) :: dustDiagSubdd_glob

      allocate(dustDiagSubdd_glob%dustEmission(im,jm,Ntm_dust))
      allocate(dustDiagSubdd_glob%dustEmission2(im,jm,Ntm_dust))
      allocate(dustDiagSubdd_glob%dustDepoTurb(im,jm,Ntm_dust))
      allocate(dustDiagSubdd_glob%dustDepoGrav(im,jm,Ntm_dust))
      allocate(dustDiagSubdd_glob%dustMassInPrec(im,jm,Ntm_dust))
      allocate(dustDiagSubdd_glob%dustSurfMixR(im,jm,Ntm_dust))
      allocate(dustDiagSubdd_glob%dustSurfConc(im,jm,Ntm_dust))
      allocate(dustDiagSubdd_glob%dustMass(im,jm,lm,Ntm_dust))
      allocate(dustDiagSubdd_glob%dustConc(im,jm,lm,Ntm_dust))

      select case(iaction)

      case(:iowrite)            ! write to restart file
        call pack_data(grid,dustDiagSubdd_acc%dustEmission
     &       ,dustDiagSubdd_glob%dustEmission)
        call pack_data(grid,dustDiagSubdd_acc%dustEmission2
     &       ,dustDiagSubdd_glob%dustEmission2)
        call pack_data(grid,dustDiagSubdd_acc%dustDepoTurb
     &       ,dustDiagSubdd_glob%dustDepoTurb)
        call pack_data(grid,dustDiagSubdd_acc%dustDepoGrav
     &       ,dustDiagSubdd_glob%dustDepoGrav)
        call pack_data(grid,dustDiagSubdd_acc%dustMassInPrec
     &       ,dustDiagSubdd_glob%dustMassInPrec)
        call pack_data(grid,dustDiagSubdd_acc%dustSurfMixR
     &       ,dustDiagSubdd_glob%dustSurfMixR)
        call pack_data(grid,dustDiagSubdd_acc%dustSurfConc
     &       ,dustDiagSubdd_glob%dustSurfConc)
        call pack_data(grid,dustDiagSubdd_acc%dustMass
     &       ,dustDiagSubdd_glob%dustMass)
        call pack_data(grid,dustDiagSubdd_acc%dustConc
     &       ,dustDiagSubdd_glob%dustConc)
        header='For subdaily dust tracers diagnostics: dustDiagSubdd'
        if (am_i_root()) then
          write(kunit,iostat=iostat) header
     &         ,dustDiagSubdd_glob%dustEmission
     &         ,dustDiagSubdd_glob%dustEmission2
     &         ,dustDiagSubdd_glob%dustDepoTurb
     &         ,dustDiagSubdd_glob%dustDepoGrav
     &         ,dustDiagSubdd_glob%dustMassInPrec
     &         ,dustDiagSubdd_glob%dustSurfMixR
     &         ,dustDiagSubdd_glob%dustSurfConc
     &         ,dustDiagSubdd_glob%dustMass
     &         ,dustDiagSubdd_glob%dustConc
          if (iostat > 0) call stop_model
     &         ('In io_trdust_drv: Restart file write error',255)
        end if

      case(ioread:)
        select case(iaction)    ! read from restart file
        case(ioread,irerun,irsfic,irsficno) ! restarts
          if (am_i_root()) then
            read(kunit,iostat=iostat) header
     &         ,dustDiagSubdd_glob%dustEmission
     &         ,dustDiagSubdd_glob%dustEmission2
     &         ,dustDiagSubdd_glob%dustDepoTurb
     &         ,dustDiagSubdd_glob%dustDepoGrav
     &         ,dustDiagSubdd_glob%dustMassInPrec
     &         ,dustDiagSubdd_glob%dustSurfMixR
     &         ,dustDiagSubdd_glob%dustSurfConc
     &         ,dustDiagSubdd_glob%dustMass
     &         ,dustDiagSubdd_glob%dustConc
            if (iostat > 0) call stop_model
     &           ('In io trdust_drv: Restart file read error',255)
          end if
          call unpack_data(grid,dustDiagSubdd_glob%dustEmission
     &         ,dustDiagSubdd_acc%dustEmission)
          call unpack_data(grid,dustDiagSubdd_glob%dustEmission2
     &         ,dustDiagSubdd_acc%dustEmission2)
          call unpack_data(grid,dustDiagSubdd_glob%dustDepoTurb
     &         ,dustDiagSubdd_acc%dustDepoTurb)
          call unpack_data(grid,dustDiagSubdd_glob%dustDepoGrav
     &         ,dustDiagSubdd_acc%dustDepoGrav)
          call unpack_data(grid,dustDiagSubdd_glob%dustMassInPrec
     &         ,dustDiagSubdd_acc%dustMassInPrec)
          call unpack_data(grid,dustDiagSubdd_glob%dustSurfMixR
     &         ,dustDiagSubdd_acc%dustSurfMixR)
          call unpack_data(grid,dustDiagSubdd_glob%dustSurfConc
     &         ,dustDiagSubdd_acc%dustSurfConc)
          call unpack_data(grid,dustDiagSubdd_glob%dustMass
     &         ,dustDiagSubdd_acc%dustMass)
          call unpack_data(grid,dustDiagSubdd_glob%dustConc
     &         ,dustDiagSubdd_acc%dustConc)
        end select

      end select

      deallocate(dustDiagSubdd_glob%dustEmission)
      deallocate(dustDiagSubdd_glob%dustEmission2)
      deallocate(dustDiagSubdd_glob%dustDepoTurb)
      deallocate(dustDiagSubdd_glob%dustDepoGrav)
      deallocate(dustDiagSubdd_glob%dustMassInPrec)
      deallocate(dustDiagSubdd_glob%dustSurfMixR)
      deallocate(dustDiagSubdd_glob%dustSurfConc)
      deallocate(dustDiagSubdd_glob%dustMass)
      deallocate(dustDiagSubdd_glob%dustConc)

      return
      end subroutine io_trDust
#endif /*TRACERS_DUST || TRACERS_MINERALS || TRACERS_QUARZHEM*/

#ifdef NEW_IO
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
c def_rsf_trdust
      subroutine def_rsf_trdust(fid)
!@sum def_rsf_trdust defines control info in restart files specifically for
!@+                  soil dust aerosols
!@auth Jan Perlwitz

      implicit none

!@var fid file id
      integer,intent(in) :: fid

      call defvar(grid,fid,dustDiagSubdd_acc%dustEmission(:,:,:)
     &     ,'dustEmission'//'(dist_im,dist_jm,Ntm_dust)')
      call defvar(grid,fid,dustDiagSubdd_acc%dustEmission2(:,:,:)
     &     ,'dustEmission2'//'(dist_im,dist_jm,Ntm_dust)')
      call defvar(grid,fid,dustDiagSubdd_acc%dustDepoTurb(:,:,:)
     &     ,'dustDepoTurb'//'(dist_im,dist_jm,Ntm_dust)')
      call defvar(grid,fid,dustDiagSubdd_acc%dustDepoGrav(:,:,:)
     &     ,'dustDepoGrav'//'(dist_im,dist_jm,Ntm_dust)')
      call defvar(grid,fid,dustDiagSubdd_acc%dustMassInPrec(:,:,:)
     &     ,' dustMassInPrec'//'(dist_im,dist_jm,Ntm_dust)')
      call defvar(grid,fid,dustDiagSubdd_acc%dustSurfMixR(:,:,:)
     &     ,'dustSurfMixR'//'(dist_im,dist_jm,Ntm_dust)')
      call defvar(grid,fid,dustDiagSubdd_acc%dustSurfConc(:,:,:),'
     &     dustSurfConc'//'(dist_im,dist_jm,Ntm_dust)')
      call defvar(grid,fid,dustDiagSubdd_acc%dustMass(:,:,:,:)
     &     ,'dustMass'//'(dist_im,dist_jm,lm,Ntm_dust)')
      call defvar(grid,fid,dustDiagSubdd_acc%dustConc(:,:,:,:)
     &     ,'dustConc'//'(dist_im,dist_jm,lm,Ntm_dust)')

      return
      end subroutine def_rsf_trdust

c new_io_trdust
      subroutine new_io_trdust(fid,iaction)
!@sum  new_io_trdust netcdf I/O specifically for soil dust aerosols
!@auth Jan Perlwitz

      implicit none

!@var fid file id
!@var iaction flag for reading or writing to file
      integer,intent(in) :: fid,iaction

      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid,fid,'dustEmission'
     &       ,dustDiagSubdd_acc%dustEmission(:,:,:))
        call write_dist_data(grid,fid,'dustEmission2'
     &       ,dustDiagSubdd_acc%dustEmission2(:,:,:))
        call write_dist_data(grid,fid,'dustDepoTurb'
     &       ,dustDiagSubdd_acc%dustDepoTurb(:,:,:))
        call write_dist_data(grid,fid,'dustDepoGrav'
     &       ,dustDiagSubdd_acc%dustDepoGrav(:,:,:))
        call write_dist_data(grid,fid,'dustMassInPrec'
     &       ,dustDiagSubdd_acc%dustMassInPrec(:,:,:))
        call write_dist_data(grid,fid,'dustSurfMixR'
     &       ,dustDiagSubdd_acc%dustSurfMixR(:,:,:))
        call write_dist_data(grid,fid,'dustSurfConc'
     &       ,dustDiagSubdd_acc%dustSurfConc(:,:,:))
        call write_dist_data(grid,fid,'dustMass'
     &       ,dustDiagSubdd_acc%dustMass(:,:,:,:))
        call write_dist_data(grid,fid,'dustConc'
     &       ,dustDiagSubdd_acc%dustConc(:,:,:,:))
      case (ioread)            ! input from restart file
        call read_dist_data(grid,fid,'dustEmission'
     &       ,dustDiagSubdd_acc%dustEmission(:,:,:))
        call read_dist_data(grid,fid,'dustEmission2'
     &       ,dustDiagSubdd_acc%dustEmission2(:,:,:))
        call read_dist_data(grid,fid,'dustDepoTurb'
     &       ,dustDiagSubdd_acc%dustDepoTurb(:,:,:))
        call read_dist_data(grid,fid,'dustDepoGrav'
     &       ,dustDiagSubdd_acc%dustDepoGrav(:,:,:))
        call read_dist_data(grid,fid,'dustMassInPrec'
     &       ,dustDiagSubdd_acc%dustMassInPrec(:,:,:))
        call read_dist_data(grid,fid,'dustSurfMixR'
     &       ,dustDiagSubdd_acc%dustSurfMixR(:,:,:))
        call read_dist_data(grid,fid,'dustSurfConc'
     &       ,dustDiagSubdd_acc%dustSurfConc(:,:,:))
        call read_dist_data(grid,fid,'dustMass'
     &       ,dustDiagSubdd_acc%dustMass(:,:,:,:))
        call read_dist_data(grid,fid,'dustConc'
     &       ,dustDiagSubdd_acc%dustConc(:,:,:,:))
      end select

      return
      end subroutine new_io_trdust
#endif /*TRACERS_DUST || TRACERS_MINERALS || TRACERS_QUARZHEM*/
#endif /*NEW_IO*/

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
c accSubddDust
      subroutine accSubddDust(dustDiagSubdd_acc)
!@sum  accSubddDust accumulates specific soil dust aerosol variables for
!@+                 subdaily diagnostics
!@auth Jan Perlwitz
      use atm_com, only : t,byam,pk,pmid
      implicit none

      type(dustDiagSubdd),intent(inout) :: dustDiagSubdd_acc

      integer :: l,n,n1

      do n=1,Ntm_dust
        n1=n_soilDust+n-1
        dustDiagSubdd_acc%dustEmission(:,:,n)
     &       =dustDiagSubdd_acc%dustEmission(:,:,n)+dust_flux_glob(:,:,n
     &       )
#if (defined TRACERS_DUST) || (defined TRACERS_AMP)
        dustDiagSubdd_acc%dustEmission2(:,:,n)
     &       =dustDiagSubdd_acc%dustEmission(:,:,n)+dust_flux2_glob(:,:
     &       ,n)
#endif
#ifdef TRACERS_DRYDEP
        if (dodrydep(n1)) then
          dustDiagSubdd_acc%dustDepoTurb(:,:,n)
     &         =dustDiagSubdd_acc%dustDepoTurb(:,:,n)
     &         +sum(depo_turb_glob(:,:,:,n1),dim=3)
          dustDiagSubdd_acc%dustDepoGrav(:,:,n)
     &         =dustDiagSubdd_acc%dustDepoGrav(:,:,n)
     &         +sum(depo_grav_glob(:,:,:,n1),dim=3)
        end if
#endif
#ifdef TRACERS_WATER
        if (dowetdep(n1)) then
          dustDiagSubdd_acc%dustMassInPrec(:,:,n)
     &         =dustDiagSubdd_acc%dustMassInPrec(:,:,n)+trprec(n1,:,:)
        end if
#else
        dustDiagSubdd_acc%dustMassInPrec(:,:,n)
     &       =dustDiagSubdd_acc%dustMassInPrec(:,:,n)+trprec_dust(n,:,:)
#endif
        dustDiagSubdd_acc%dustSurfMixR(:,:,n)
     &       =dustDiagSubdd_acc%dustSurfMixR(:,:,n)+trcSurf(:,:,n1)
        dustDiagSubdd_acc%dustSurfConc(:,:,n)
     &       =dustDiagSubdd_acc%dustSurfConc(:,:,n)+trcSurfByVol(:,:,n1)
        dustDiagSubdd_acc%dustMass(:,:,:,n)=dustDiagSubdd_acc%dustMass(:
     &       ,:,:,n)+trm(:,:,:,n1)
        do l=1,lm
          dustDiagSubdd_acc%dustConc(:,:,l,n) =
     &         dustDiagSubdd_acc%dustConc(:,:,l,n) + trm(:,:,l,n1)
     &         *byam(l,:,:)*1d2*pmid(l,:,:)/(rgas*t(:,:,l)*pk(l,:,:))
        end do
      end do

      return
      end subroutine accSubddDust
#endif /*TRACERS_DUST || TRACERS_MINERALS || TRACERS_QUARZHEM*/

      end module trdust_drv
