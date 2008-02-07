      module hycom_cpler
      USE HYCOM_DIM, only : iia,jja,iio,jjo,isp,ifp,ilp,ii,jj,ip
      implicit none
!!#include "a2o.h"

      integer nwgta2o,nwgta2o2,nwgto2a
      parameter (nwgta2o=18,nwgta2o2=18,nwgto2a=39)
c
      real*8 wlista2o(iio,jjo,nwgta2o),wtaua2o(iio,jjo,nwgta2o2)
     .    ,wlisto2a(iia,jja,nwgto2a)
     .    ,wlisto2a_e(iia,jja,nwgto2a)
     .    ,wlisto2a_n(iia,jja,nwgto2a)
     .    ,ofrac(iia,jja),agcmgdsz(jja)
     .    ,ocellsz(iio,jjo),coso(iio,jjo),sino(iio,jjo)
      integer ilista2o(iio,jjo,nwgta2o),jlista2o(iio,jjo,nwgta2o)
     .                                 ,nlista2o(iio,jjo)
     .       ,itaua2o(iio,jjo,nwgta2o2), jtaua2o(iio,jjo,nwgta2o2)
     .                                 , ntaua2o(iio,jjo)
      integer ilisto2a(iia,jja,nwgto2a),jlisto2a(iia,jja,nwgto2a)
     .                                 ,nlisto2a(iia,jja)
     .     ,ilisto2a_e(iia,jja,nwgto2a),jlisto2a_e(iia,jja,nwgto2a)
     .                                 ,nlisto2a_e(iia,jja)
     .     ,ilisto2a_n(iia,jja,nwgto2a),jlisto2a_n(iia,jja,nwgto2a)
     .                                 ,nlisto2a_n(iia,jja)

      contains

      subroutine ssto2a(fldo,flda)
c --- mapping sst from 'o' grid to 'a' grd
c
c --- fldo:  input field from ogcm grid
c     flda: output field onto agcm grid
c
      !USE HYCOM_DIM
      implicit none
!!#include "dimensions.h"
#include "dimension2.h"
!!#include "a2o.h"
      real*8 flda(iia,jja),fldo(iio,jjo),tto,tta
c
c$OMP PARALLEL DO
      do 16 ja=1,jja
      do 16 ia=1,iia
      flda(ia,ja)=0.
c
      do 17 n=1,nlisto2a(ia,ja)
 17   flda(ia,ja)=flda(ia,ja)+fldo(ilisto2a(ia,ja,n),jlisto2a(ia,ja,n))
     .                                              *wlisto2a(ia,ja,n)
 16   continue
c$OMP END PARALLEL DO             
c
      return
      end subroutine ssto2a
c
      subroutine veca2o(tauxa,tauya,tauxo,tauyo)
c --- tauxa/tauya: input taux (E-ward)/tauy (N-ward) on agcm grid (N/m*m)
c --- tauxo/tauyo:output taux (S-ward)/tauy (E-ward) on ogcm grid (N/m*m)
c    
      !USE HYCOM_DIM
      implicit none
!!#include "dimensions.h"
#include "dimension2.h"
!!#include "a2o.h"
      real*8 tauxa(iia,jja),tauya(iia,jja),tauxo(iio,jjo),tauyo(iio,jjo)
     .    ,sward(iio,jjo),eward(iio,jjo),tta,tto
c
cdiag tto=0.
c --- mapping tauxa/tauya to ogcm grid
c$OMP PARALLEL DO
      do 6 j=1,jj               
      do 6 l=1,isp(j)           
      do 6 i=ifp(j,l),ilp(j,l)
      eward(i,j)=0.
      sward(i,j)=0.
c
      do 7 n=1,ntaua2o(i,j)
      eward(i,j)=eward(i,j)+tauxa(itaua2o(i,j,n),jtaua2o(i,j,n))
     .                                          *wtaua2o(i,j,n)
 7    sward(i,j)=sward(i,j)-tauya(itaua2o(i,j,n),jtaua2o(i,j,n))
     .                                          *wtaua2o(i,j,n)
 6    continue
c$OMP END PARALLEL DO
c
c --- rotate sward/eward to fit onto Panam grid
c$OMP PARALLEL DO
      do 9 j=1,jj
      do 9 l=1,isp(j)           
      do 9 i=ifp(j,l),ilp(j,l)
      tauxo(i,j)= sward(i,j)*coso(i,j)+eward(i,j)*sino(i,j)
      tauyo(i,j)= eward(i,j)*coso(i,j)-sward(i,j)*sino(i,j)
 9    continue
c$OMP END PARALLEL DO           
      return
      end subroutine veca2o
c
      subroutine flxa2o(flda,fldo)
c --- mapping flux-like field from agcm to ogcm
c     input: flda (W/m*m), output: fldo (W/m*m)
c
      !USE HYCOM_DIM
      implicit none
!!#include "dimensions.h"
#include "dimension2.h"
!!#include "a2o.h"
      real*8 flda(iia,jja),fldo(iio,jjo),tta,tto
c
c$OMP PARALLEL DO
      do 8 j=1,jj
      do 8 l=1,isp(j)
      do 8 i=ifp(j,l),ilp(j,l)
      fldo(i,j)=0.
c
      do 9 n=1,nlista2o(i,j)
      fldo(i,j)=fldo(i,j)+flda(ilista2o(i,j,n),jlista2o(i,j,n))
     .                        *wlista2o(i,j,n)
 9    continue
 8    continue
c$OMP END PARALLEL DO
c
      return
      end subroutine flxa2o
c
      subroutine veco2a(tauxo,tauyo,tauxa,tauya)
c --- tauxo/tauyo: input taux (S-ward)/tauy (E-ward) on ogcm grid (N/m*m)
c --- tauxa/tauya:output taux (E-ward)/tauy (N-ward) on agcm grid (N/m*m)
c    
      !USE HYCOM_DIM
      implicit none
!!#include "dimensions.h"
#include "dimension2.h"
!!#include "a2o.h"
      real*8 tauxa(iia,jja),tauya(iia,jja),tauxo(iio,jjo),tauyo(iio,jjo)
     .      ,nward(iio,jjo),eward(iio,jjo),tta,tto,sine
c
c$OMP PARALLEL DO
      do 10 j=1,jj
      do 10 i=1,ii
      nward(i,j)=0.
 10   eward(i,j)=0.
c$OMP END PARALLEL DO
c
c --- rotate taux/tauy to n/e orientation at local p point on panam grid
c --- check velocity bounds
c$OMP PARALLEL DO PRIVATE(jb,sine)
      do 12 j=1,jj
      jb=mod(j,jj)+1
      do 12 l=1,isp(j)           
      do 12 i=ifp(j,l),ilp(j,l)
      if (ip(i,j).eq.1) then
      sine=sino(i,j)*sino(i,j)+coso(i,j)*coso(i,j)
      nward(i,j)=((tauyo(i,j)+tauyo(i ,jb))*sino(i,j)
     .           -(tauxo(i,j)+tauxo(i+1,j))*coso(i,j))/(2.*sine)
      eward(i,j)=((tauyo(i,j)+tauyo(i ,jb))*coso(i,j)
     .           +(tauxo(i,j)+tauxo(i+1,j))*sino(i,j))/(2.*sine)
      endif
 12   continue
c$OMP END PARALLEL DO           
c
c --- mapping nward/eward from ogcm grid to agcm grid
c
c$OMP PARALLEL DO
      do 16 ja=1,jja
      do 16 ia=1,iia
      tauxa(ia,ja)=0.
      tauya(ia,ja)=0.
c
      do 17 n=1,nlisto2a_e(ia,ja)
 17   tauxa(ia,ja)=tauxa(ia,ja)+eward(ilisto2a_e(ia,ja,n)
     .          ,jlisto2a_e(ia,ja,n))*wlisto2a_e(ia,ja,n)
c
      do 18 n=1,nlisto2a_n(ia,ja)
 18   tauya(ia,ja)=tauya(ia,ja)+nward(ilisto2a_n(ia,ja,n)
     .          ,jlisto2a_n(ia,ja,n))*wlisto2a_n(ia,ja,n)
 16   continue
c$OMP END PARALLEL DO
c
      return
      end subroutine veco2a
c
      subroutine tempro2a(fldo,flda)
c --- mapping sqrt(sqrt(temp**4)) from 'o' grid to 'a' grid, unit: K
c
c --- fldo:  input field from ogcm grid
c     flda: output field onto agcm grid
c
      !USE HYCOM_DIM
      implicit none
!!#include "dimensions.h"
#include "dimension2.h"
!!#include "a2o.h"
      real*8 flda(iia,jja),fldo(iio,jjo),tto,tta
c
c$OMP PARALLEL DO
      do 16 ja=1,jja
      do 16 ia=1,iia
      flda(ia,ja)=0.
c
      do 17 n=1,nlisto2a(ia,ja)
 17   flda(ia,ja)=flda(ia,ja)+(fldo(ilisto2a(ia,ja,n),jlisto2a(ia,ja,n))
     .    +273.16d0)**4*wlisto2a(ia,ja,n)
      flda(ia,ja)=sqrt(sqrt(flda(ia,ja)))       ! Kelvin for radiation
 16   continue
c$OMP END PARALLEL DO             
c
      return
      end subroutine tempro2a

      end module hycom_cpler
