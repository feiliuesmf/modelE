#include "rundeck_opts.h"
      MODULE obio_forc

      USE obio_dim

      USE hycom_dim_glob
      implicit none

!!#include "dimensions.h"
#include "dimension2.h"


!     real sunz               !solar zenith angle
!     real solz_all,solz      !mean cosine solar zenith angle
!     real solz2,sunz2
!     common /bsolz/solz_all(idm,jdm,12,12),solz2(12),sunz2(12)

      real asolz(iia,jja)
      real osolz(idm,jdm)
      real solz               !mean cosine solar zenith angle
      real sunz               !solar zenith angle
      common /brod1/ solz,sunz
!$OMP THREADPRIVATE(/brod1/)

#ifdef OBIO_RAD_coupling 
      real eda_frac,esa_frac
      common /frac_oasim/eda_frac(nlt),esa_frac(nlt)

      real avisdir,avisdif,anirdir,anirdif
      real ovisdir,ovisdif,onirdir,onirdif
      common /rada2o/  avisdir(iia,jja),avisdif(iia,jja)
     .                ,anirdir(iia,jja),anirdif(iia,jja)
     .                ,ovisdir(idm,jdm),ovisdif(idm,jdm)
     .                ,onirdir(idm,jdm),onirdif(idm,jdm)
      real ovisdir_ij,ovisdif_ij,onirdir_ij,onirdif_ij
      common /rada2o_ij/ ovisdir_ij,ovisdif_ij,onirdir_ij,onirdif_ij
!$OMP THREADPRIVATE(/rada2o_ij/)

#ifdef CHL_from_SeaWIFs
      real, ALLOCATABLE :: chl_3d(:,:,:)
      real achl(iia,jja)
      real chl
      common /alb_rad_chl/chl
!$OMP THREADPRIVATE(/alb_rad_chl/)
#endif

!      real obio_bocvn,obio_xocvn   !albedo bands that are used in modelE
!      common /alb_rad/ obio_bocvn(6),obio_xocvn(6)
!!$OMP THREADPRIVATE(/alb_rad/)

#else
      !real Eda(idm,jdm,nlt,nhn,12)    !direct downwelling irradiance
      !real Esa(idm,jdm,nlt,nhn,12)    !diffuse downwelling irradiance
      !common /beda/  Eda(idm,jdm,nlt,nhn,12),Esa(idm,jdm,nlt,nhn,12)
      real, ALLOCATABLE :: Eda(:,:,:,:,:),Esa(:,:,:,:,:)
      real Eda2,Esa2,Ed,Es 
      common /beda2/ Eda2(nlt,nhn),Esa2(nlt,nhn)
!$OMP THREADPRIVATE(/beda2/)
#endif

      real Ed,Es 
      common /beds/  Ed(nlt),Es(nlt)
!$OMP THREADPRIVATE(/beds/)

      real awind      !wind speed from modelE (see hycom2.f)
      real owind      !wind speed in hycom grid (see hycom2.f)
      common /owind/  awind(iia,jja),owind(idm,jdm)
     
      real    rod,ros       !surface reflectance for
                            !direct (rod) and diffuse (ros)
                            !components separately
      common /brod2/ rod(nlt),ros(nlt)
!$OMP THREADPRIVATE(/brod2/)

      real wind               !surface wind from atmos
      common /bwind/ wind
!$OMP THREADPRIVATE(/bwind/)

      real tirrq3d
      common /btirrq/ tirrq3d(idm,jdm,kdm)
      real tirrq                   !total mean irradiance in quanta
      common /blte/ tirrq(kdm)
!$OMP THREADPRIVATE(/blte/)

      integer ihra                 !counter for daylight hours
      common /bhra/ ihra(idm,jdm)

      real rmud                    !downwelling irradiance average cosine
      common /bmud /rmud 
!$OMP THREADPRIVATE(/bmud/)

      real avgq    !mean daily irradiance in quanta
      common /bavq/ avgq(idm,jdm,kdm)

      real atmFe,atmFe_all
      common /biron/atmFe(idm,jdm),atmFe_all(idm,jdm,12) !surface iron deposition

      real alk
      common /balk/alk(idm,jdm,kdm)    !alkalinity from climatology in 'umol/kg'

      real, parameter :: atmCO2=368.6          !uatm for year 2000
                   !     atmCO2=280.           !uatm for preindustial runs
                   !     atmCO2=371.3          !uatm or ppmv (equivalent);
                                               !global mean
                                               !2000-2003 from OCMIP

      END MODULE obio_forc
