#include "rundeck_opts.h"
      MODULE obio_forc

      USE obio_dim

      USE hycom_dim_glob
      implicit none


      integer, ALLOCATABLE, DIMENSION(:,:) :: ihra            !counter for daylight hours

      !real, ALLOCATABLE, DIMENSION(:,:)    :: asolz
      !real, ALLOCATABLE, DIMENSION(:,:)    :: awind           !wind speed from modelE (see hycom2.f)
      real, ALLOCATABLE, DIMENSION(:,:)    :: osolz
      real, ALLOCATABLE, DIMENSION(:,:)    :: owind           !wind speed in hycom grid (see hycom2.f)
      real, ALLOCATABLE, DIMENSION(:,:)    :: osolz_glob
      real, ALLOCATABLE, DIMENSION(:,:)    :: owind_glob           !wind speed in hycom grid (see hycom2.f)

      real, ALLOCATABLE, DIMENSION(:,:,:)  :: tirrq3d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: avgq            !mean daily irradiance in quanta
      real, ALLOCATABLE, DIMENSION(:,:)    :: atmFe
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: atmFe_all       !surface iron deposition
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: alk             !alkalinity from climatology in 'umol/kg'

#ifdef OBIO_RAD_coupling
!      real*8, ALLOCATABLE, DIMENSION(:,:)    :: avisdir,avisdif
!     .                                         ,anirdir,anirdif
      real*8, ALLOCATABLE, DIMENSION(:,:)    :: ovisdir,ovisdif
     .                                         ,onirdir,onirdif
      real*8, ALLOCATABLE, DIMENSION(:,:)   :: ovisdir_glob,ovisdif_glob
     .                                        ,onirdir_glob,onirdif_glob
#endif
#ifndef OBIO_RAD_coupling
      real, ALLOCATABLE, DIMENSION(:,:,:,:,:):: Eda,Esa       !direct,diffuse downwelling irradiance
#endif


      real solz               !mean cosine solar zenith angle
      real sunz               !solar zenith angle
      common /brod1/ solz,sunz
!$OMP THREADPRIVATE(/brod1/)

#ifdef OBIO_RAD_coupling 
      real eda_frac,esa_frac
      common /frac_oasim/eda_frac(nlt),esa_frac(nlt)

      real ovisdir_ij,ovisdif_ij,onirdir_ij,onirdif_ij
      common /rada2o_ij/ ovisdir_ij,ovisdif_ij,onirdir_ij,onirdif_ij
!$OMP THREADPRIVATE(/rada2o_ij/)
#else
      real Eda2,Esa2
      common /beda2/ Eda2(nlt,nhn),Esa2(nlt,nhn)
!$OMP THREADPRIVATE(/beda2/)
#endif

      real Ed,Es 
      common /beds/  Ed(nlt),Es(nlt)
!$OMP THREADPRIVATE(/beds/)

      real wind               !surface wind from atmos
      common /bwind/ wind
!$OMP THREADPRIVATE(/bwind/)

      real tirrq                   !total mean irradiance in quanta
      common /blte/ tirrq(kdm)
!$OMP THREADPRIVATE(/blte/)

      real rmud                    !downwelling irradiance average cosine
      common /bmud /rmud 
!$OMP THREADPRIVATE(/bmud/)

      real atmCO2
!     real, parameter :: atmCO2=368.6          !uatm for year 2000
!                  !     atmCO2=280.           !uatm for preindustial runs
!                  !     atmCO2=371.3          !uatm or ppmv (equivalent);
!                                              !global mean
!                                              !2000-2003 from OCMIP

      contains

      subroutine alloc_obio_forc

      USE obio_dim
      USE hycom_dim_glob
      USE hycom_dim, only : j_0h,j_1h

      !ALLOCATE(asolz(iia,jja))
      !ALLOCATE(awind(iia,jja))
      ALLOCATE(osolz(idm,j_0h:j_1h))
      ALLOCATE(owind(idm,j_0h:j_1h))
      ALLOCATE(osolz_glob(idm,jdm))
      ALLOCATE(owind_glob(idm,jdm))

      ALLOCATE(tirrq3d(idm,j_0h:j_1h,kdm))
      ALLOCATE(ihra(idm,j_0h:j_1h))
      ALLOCATE(avgq(idm,j_0h:j_1h,kdm))
!      ALLOCATE(atmFe(idm,jdm),atmFe_all(idm,jdm,12))
!      ALLOCATE(alk(idm,jdm,kdm))
      ALLOCATE(atmFe(idm,jdm),atmFe_all(idm,j_0h:j_1h,12))
      ALLOCATE(alk(idm,j_0h:j_1h,kdm))


#ifdef OBIO_RAD_coupling
!      ALLOCATE(avisdir(iia,jja),avisdif(iia,jja)
!     .        ,anirdir(iia,jja),anirdif(iia,jja))
      ALLOCATE(ovisdir(idm,j_0h:j_1h),ovisdif(idm,j_0h:j_1h)
     .        ,onirdir(idm,j_0h:j_1h),onirdif(idm,j_0h:j_1h))
      ALLOCATE(ovisdir_glob(idm,jdm),ovisdif_glob(idm,jdm)
     .        ,onirdir_glob(idm,jdm),onirdif_glob(idm,jdm))
#endif
cddd#ifndef OBIO_RAD_coupling
cddd      ALLOCATE(Eda(idm,jdm,nlt,nhn,12))
cddd      ALLOCATE(Esa(idm,jdm,nlt,nhn,12))
cddd#endif

      end subroutine alloc_obio_forc

      subroutine scatter_obio_forc_arrays

      USE HYCOM_DIM, only : ogrid
      USE DOMAIN_DECOMP, ONLY: UNPACK_DATA

      call unpack_data( ogrid, owind_glob, owind )
      call unpack_data( ogrid, osolz_glob, osolz )
#ifdef OBIO_RAD_coupling
      call unpack_data( ogrid, ovisdir_glob, ovisdir )
      call unpack_data( ogrid, ovisdif_glob, ovisdif )
      call unpack_data( ogrid, onirdir_glob, onirdir )
      call unpack_data( ogrid, onirdif_glob, onirdif )
#endif

      end subroutine scatter_obio_forc_arrays

      END MODULE obio_forc
