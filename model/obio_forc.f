      MODULE obio_forc

      USE obio_dim

      implicit none

      SAVE 

#include "dimensions.h"
#include "dimension2.h"


!     real sunz               !solar zenith angle
!     real solz_all,solz      !mean cosine solar zenith angle
!     real solz2,sunz2
!     common /bsolz/solz_all(idm,jdm,12,12),solz2(12),sunz2(12)

      real asolz(iia,jja)
      real osolz(idm,jdm)
      real sunz               !solar zenith angle
      real solz               !mean cosine solar zenith angle

      !real Eda(idm,jdm,nlt,nhn,12)    !direct downwelling irradiance
      !real Esa(idm,jdm,nlt,nhn,12)    !diffuse downwelling irradiance
      !common /beda/  Eda(idm,jdm,nlt,nhn,12),Esa(idm,jdm,nlt,nhn,12)
      real, ALLOCATABLE :: Eda(:,:,:,:,:),Esa(:,:,:,:,:)

      real Eda2,Esa2,Ed,Es 
      common /beda2/ Eda2(nlt,nhn),Esa2(nlt,nhn)
      common /beds/  Ed(nlt),Es(nlt)
!$OMP THREADPRIVATE(/beds/)

      real awind      !wind speed from modelE (see hycom2.f)
      real owind      !wind speed in hycom grid (see hycom2.f)
      common /owind/  awind(iia,jja),owind(idm,jdm)
     

      real    rod,ros       !surface reflectance for
                            !direct (rod) and diffuse (ros)
                            !components separately
      common /brod/ rod(nlt),ros(nlt),solz,sunz
!$OMP THREADPRIVATE(/brod/)

      real wind               !surface wind from atmos
      common /bwind/wind
!$OMP THREADPRIVATE(/bwind/)

      real tirrq  !total mean irradiance in quanta
      common /blte/ tirrq(kdm)
!$OMP THREADPRIVATE(/blte/)

      integer ihra !counter for daylight hours
      common /bhra/ ihra(idm,jdm)

      real rmud          !downwelling irradiance average cosine
      common /bmud /rmud 

      real avgq    !mean daily irradiance in quanta
      common /bavq/ avgq(idm,jdm,kdm)

      real atmFe,atmFe_all
      common /biron/atmFe(idm,jdm),atmFe_all(idm,jdm,12) !surface iron deposition

      real alk
      common /balk/alk(idm,jdm,kdm)    !alkalinity from climatology in 'umol/kg'


      real, parameter :: atmCO2=368.6          !uatm for year 2000
                   !     atmCO2=371.3          !uatm or ppmv (equivalent); 
                                               !global mean
                                               !2000-2003 from OCMIP

      END MODULE obio_forc
