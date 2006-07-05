      subroutine ssto2a(fldo,flda)
c --- mapping sst from 'o' grid to 'a' grd
c
c --- fldo:  input field from ogcm grid
c     flda: output field onto agcm grid
c
      implicit none
#include "dimensions.h"
#include "dimension2.h"
#include "a2o.h"
      real*8 flda(iia,jja),fldo(iio,jjo),tto,tta
      integer io,jo
c
cdiag tta=0.
cdiag tto=0.
c$OMP PARALLEL DO PRIVATE(io,jo)
      do 16 ja=1,jja
      do 16 ia=1,iia
      flda(ia,ja)=0.
c
      do 17 n=1,nlisto2a(ia,ja)
      io=ilisto2a(ia,ja,n) 
      jo=jlisto2a(ia,ja,n) 
 17   flda(ia,ja)=flda(ia,ja)+fldo(io,jo)*wlisto2a(ia,ja,n)
cdiag tta=tta+flda(ia,ja)*ofrac(ia,ja)*agcmgdsz(ja)
 16   continue
c$OMP END PARALLEL DO             
c
cdiag      do 18 io=1,iio
cdiag      do 18 jo=1,jjo
cdiag 18   tto=tto+fldo(io,jo)*ocellsz(io,jo)
cdiag      if (tto.ne.0.) then
cdiag        write(*,'(a,2e16.4,f12.8)')'chk bud ssto2a=',tta,tto,tta/tto
cdiag      else
cdiag        write(*,'(a,2e16.4)') 'chk bud ssto2a=',tta,tto
cdiag      endif
cdiag      if (abs(tta-tto)/tta.gt. 0.1) stop 'budget'
c
      return
      end
c
      subroutine veca2o(tauxa,tauya,tauxo,tauyo)
c --- tauxa/tauya: input taux (E-ward)/tauy (N-ward) on agcm grid (N/m*m)
c --- tauxo/tauyo:output taux (S-ward)/tauy (E-ward) on ogcm grid (N/m*m)
c    
      implicit none
#include "dimensions.h"
#include "dimension2.h"
#include "a2o.h"
      real*8 tauxa(iia,jja),tauya(iia,jja),tauxo(iio,jjo),tauyo(iio,jjo)
     .    ,sward(iio,jjo),eward(iio,jjo),tta,tto
c
cdiag tto=0.
c --- mapping tauxa/tauya to ogcm grid
c$OMP PARALLEL DO PRIVATE(ia,ja) !!! REDUCTION(+:tto)
      do 6 j=1,jj               
      do 6 l=1,isp(j)           
      do 6 i=ifp(j,l),ilp(j,l)
      eward(i,j)=0.
      sward(i,j)=0.
c
      do 7 n=1,ntaua2o(i,j)
      ia=itaua2o(i,j,n) 
      ja=jtaua2o(i,j,n) 
      eward(i,j)=eward(i,j)+tauxa(ia,ja)*wtaua2o(i,j,n)
 7    sward(i,j)=sward(i,j)-tauya(ia,ja)*wtaua2o(i,j,n)
cdiag tto=tto+eward(i,j)*ocellsz(i,j)
 6    continue
c$OMP END PARALLEL DO
c
cdiag tta=0.
cc$OMP PARALLEL DO REDUCTION(+:tta)
cdiag do 8 ia=1,iia
cdiag do 8 ja=1,jja
cdiag 8    tta=tta+tauxa(ia,ja)*agcmocn(ia,ja)
cc$OMP END PARALLEL DO
cdiag  if (tto.ne.0.)
cdiag  .write(*,'(a,2e16.4,f12.8)')'chk bud veca2o=',tta,tto,tta/tto
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
      end
c
      subroutine flxa2o(flda,fldo)
c --- mapping flux-like field from agcm to ogcm
c     input: flda (W/m*m), output: fldo (W/m*m)
c
      implicit none
#include "dimensions.h"
#include "dimension2.h"
#include "a2o.h"
      real*8 flda(iia,jja),fldo(iio,jjo),tta,tto
c
cdiag      tto=0.
cdiag      tta=0.
cdiag      do 18 ia=1,iia
cdiag      do 18 ja=1,jja
cdiag      tta=tta+flda(ia,ja)*agcmocn(ia,ja)
cdiag 18   continue
c
c$OMP PARALLEL DO PRIVATE(ia,ja)
      do 8 j=1,jj
      do 8 l=1,isp(j)
      do 8 i=ifp(j,l),ilp(j,l)
      fldo(i,j)=0.
c
      do 9 n=1,nlista2o(i,j)
      ia=ilista2o(i,j,n)
      ja=jlista2o(i,j,n)
      fldo(i,j)=fldo(i,j)+flda(ia,ja)*wlista2o(i,j,n)
 9    continue
cdiag      tto=tto+fldo(i,j)*ocellsz(i,j)
 8    continue
c$OMP END PARALLEL DO
c
cdiag      if (tto.ne.0.) then
cdiag        write(*,'(a,2e16.4,a,f12.8)')'chk bud flxa2m=',tta,tto
cdiag     .                                          ,' %=',tta/tto
cdiag      else
cdiag        write(*,'(a,2e16.4)')'chk bud flxa2m=',tta,tto
cdiag      endif
      return
      end
c
      subroutine veco2a(tauxo,tauyo,tauxa,tauya)
c --- tauxo/tauyo: input taux (S-ward)/tauy (E-ward) on ogcm grid (N/m*m)
c --- tauxa/tauya:output taux (E-ward)/tauy (N-ward) on agcm grid (N/m*m)
c    
      implicit none
#include "dimensions.h"
#include "dimension2.h"
#include "a2o.h"
      real*8 tauxa(iia,jja),tauya(iia,jja),tauxo(iio,jjo),tauyo(iio,jjo)
     .      ,nward(iio,jjo),eward(iio,jjo),tta,tto,sin
      integer io,jo
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
c$OMP PARALLEL DO PRIVATE(jb,sin)
      do 12 j=1,jj
      jb=mod(j,jj)+1
      do 12 l=1,isp(j)           
      do 12 i=ifp(j,l),ilp(j,l)
      if (iv(i,j)+iv(i,jb)+iu(i,j)+iu(i+1,j).eq.4) then
      sin=sino(i,j)*sino(i,j)+coso(i,j)*coso(i,j)
      nward(i,j)=((tauyo(i,j)+tauyo(i ,jb))*sino(i,j)
     .           -(tauxo(i,j)+tauxo(i+1,j))*coso(i,j))/(2.*sin)
      eward(i,j)=((tauyo(i,j)+tauyo(i ,jb))*coso(i,j)
     .           +(tauxo(i,j)+tauxo(i+1,j))*sino(i,j))/(2.*sin)
      endif
 12   continue
c$OMP END PARALLEL DO           
c
c --- mapping nward/eward from ogcm grid to agcm grid
c
c$OMP PARALLEL DO PRIVATE(io,jo)
      do 16 ja=1,jja
      do 16 ia=1,iia
      tauxa(ia,ja)=0.
      tauya(ia,ja)=0.
c
      do 17 n=1,nlisto2a(ia,ja)
      io=ilisto2a(ia,ja,n)
      jo=jlisto2a(ia,ja,n)
      tauxa(ia,ja)=tauxa(ia,ja)+eward(io,jo)*wlisto2a(ia,ja,n)
 17   tauya(ia,ja)=tauya(ia,ja)+nward(io,jo)*wlisto2a(ia,ja,n)
 16   continue
c$OMP END PARALLEL DO
c
      return
      end
