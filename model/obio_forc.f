#include "rundeck_opts.h"
      MODULE obio_forc

      USE obio_dim

      USE hycom_dim_glob
      implicit none

!!#include "dimensions.h"
#include "dimension2.h"

      integer, ALLOCATABLE, DIMENSION(:,:) :: ihra            !counter for daylight hours

      real, ALLOCATABLE, DIMENSION(:,:)    :: osolz
      real, ALLOCATABLE, DIMENSION(:,:)    :: owind           !wind speed in hycom grid (see hycom2.f)
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: tirrq3d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: avgq            !mean daily irradiance in quanta
      real, ALLOCATABLE, DIMENSION(:,:)    :: atmFe
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: atmFe_all       !surface iron deposition
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: alk             !alkalinity from climatology in 'umol/kg'

      real, ALLOCATABLE, DIMENSION(:,:)    :: asolz
      real, ALLOCATABLE, DIMENSION(:,:)    :: awind           !wind speed from modelE (see hycom2.f)
#ifdef OBIO_RAD_coupling
      real*8, ALLOCATABLE, DIMENSION(:,:)    :: ovisdir,ovisdif
     .                                         ,onirdir,onirdif
      real*8, ALLOCATABLE, DIMENSION(:,:)    :: avisdir,avisdif
     .                                         ,anirdir,anirdif
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

!     real    rod,ros       !surface reflectance for
!                           !direct (rod) and diffuse (ros)
!                           !components separately
!     common /brod2/ rod(nlt),ros(nlt)
!!$OMP THREADPRIVATE(/brod2/)

      real wind               !surface wind from atmos
      common /bwind/ wind
!$OMP THREADPRIVATE(/bwind/)

      real tirrq                   !total mean irradiance in quanta
      common /blte/ tirrq(kdm)
!$OMP THREADPRIVATE(/blte/)

      real rmud                    !downwelling irradiance average cosine
      common /bmud /rmud 
!$OMP THREADPRIVATE(/bmud/)

      real, parameter :: atmCO2=368.6          !uatm for year 2000
                   !     atmCO2=280.           !uatm for preindustial runs
                   !     atmCO2=371.3          !uatm or ppmv (equivalent);
                                               !global mean
                                               !2000-2003 from OCMIP

      contains

      subroutine alloc_obio_forc

      USE obio_dim
      USE hycom_dim_glob

      ALLOCATE(osolz(idm,jdm))
      ALLOCATE(owind(idm,jdm))
      ALLOCATE(tirrq3d(idm,jdm,kdm))
      ALLOCATE(ihra(idm,jdm))
      ALLOCATE(avgq(idm,jdm,kdm))
      ALLOCATE(atmFe(idm,jdm),atmFe_all(idm,jdm,12))
      ALLOCATE(alk(idm,jdm,kdm))

      ALLOCATE(asolz(iia,jja))
      ALLOCATE(awind(iia,jja))

#ifdef OBIO_RAD_coupling
      ALLOCATE(ovisdir(idm,jdm),ovisdif(idm,jdm)
     .        ,onirdir(idm,jdm),onirdif(idm,jdm))
      ALLOCATE(avisdir(iia,jja),avisdif(iia,jja)
     .        ,anirdir(iia,jja),anirdif(iia,jja))
#endif
#ifndef OBIO_RAD_coupling
      ALLOCATE(Eda(idm,jdm,nlt,nhn,12))
      ALLOCATE(Esa(idm,jdm,nlt,nhn,12))
#endif

      end subroutine alloc_obio_forc


      END MODULE obio_forc
