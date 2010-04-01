#include "rundeck_opts.h"
      module hycom_cpler
      USE HYCOM_DIM_GLOB, only : iia,jja,iio,jjo,isp,ifp,ilp,ii,jj,ip
      USE HYCOM_SCALARS, only : flnma2o,flnma2o_s,flnmo2a,flnmo2a_f
     &   ,flnma2o_tau,flnmcoso,lp
      USE HYCOM_DIM, only : agrid,ogrid
     &    ,aJ_0, aJ_1, aJ_0H, aJ_1H,
     &      J_0,  J_1,  J_0H,  J_1H
      USE DOMAIN_DECOMP_1D, only : am_i_root,pack_data,unpack_data
      USE HYCOM_DIM, only : agrid,ogrid
c
      implicit none
      private

      public ssta2o,ssto2a,veca2o,flxa2o,flxo2a,veco2a,tempro2a,cpl_wgt
      public ssto2a_global,flxa2o_global

      public nwgta2o,nwgto2a

      public wlista2o, wtaua2o, wlista2o_s, wlisto2a, wlisto2a_f
     .    ,ilista2o_s, jlista2o_s,nlista2o_s, coso, sino
     .     ,ilista2o,  jlista2o,  nlista2o
     .     ,itaua2o,   jtaua2o,   ntaua2o
     .     ,ilisto2a,  jlisto2a,  nlisto2a
     .     ,ilisto2a_f,jlisto2a_f,nlisto2a_f

      integer nwgta2o,nwgto2a
#ifdef ATM4x5_HYCOM2deg
      parameter (nwgta2o=18,nwgto2a=39)
#endif
#ifdef ATM2x2h_HYCOM2deg
      parameter (nwgta2o=36,nwgto2a=19)
#endif
#ifdef ATM2x2h_HYCOM1deg
      parameter (nwgta2o=37,nwgto2a=48)
#endif
c
      real*8 wlista2o(iio,jjo,nwgta2o),wtaua2o(iio,jjo,nwgta2o)
     .    ,wlista2o_s(iio,jjo,nwgta2o)
     .    ,wlisto2a(iia,jja,nwgto2a), wlisto2a_f(iia,jja,nwgto2a)
     .    ,coso(iio,jjo),sino(iio,jjo)
      integer ilista2o_s(iio,jjo,nwgta2o),jlista2o_s(iio,jjo,nwgta2o)
     .                                 ,nlista2o_s(iio,jjo)
     .       ,ilista2o(iio,jjo,nwgta2o),jlista2o(iio,jjo,nwgta2o)
     .                                 ,nlista2o(iio,jjo)
     .       ,itaua2o(iio,jjo,nwgta2o), jtaua2o(iio,jjo,nwgta2o)
     .                                 , ntaua2o(iio,jjo)
      integer ilisto2a(iia,jja,nwgto2a),jlisto2a  (iia,jja,nwgto2a)
     .                                 ,nlisto2a    (iia,jja)
     .       ,ilisto2a_f(iia,jja,nwgto2a),jlisto2a_f(iia,jja,nwgto2a)
     .                                 ,nlisto2a_f  (iia,jja)

      contains
      subroutine ssta2o(flda,fldo)
c --- mapping scalar-like field from agcm A grid to ogcm A grid
c     input: flda, output: fldo 
c
      implicit none
      integer i,j,l,n
      real*8, intent(in)  :: flda(iia,jja)
      real*8, intent(out) :: fldo(iio,jjo)
c
c$OMP PARALLEL DO
      do 8 j=1,jj
      do 8 l=1,isp(j)
      do 8 i=ifp(j,l),ilp(j,l)
      fldo(i,j)=0.
c
      do 9 n=1,nlista2o_s(i,j)
      fldo(i,j)=fldo(i,j)+flda(ilista2o_s(i,j,n),jlista2o_s(i,j,n))
     .                        *wlista2o_s(i,j,n)
 9    continue
 8    continue
c$OMP END PARALLEL DO
c
      return
      end subroutine ssta2o

      subroutine ssto2a(fldo_loc,flda_loc)
c --- mapping sst from ogcm A grid to agcm A grid
c     input: fldo_loc, output: flda_loc
c
      implicit none
      integer n,ia,ja
      real*8 fldo_loc(iio,J_0H:J_1H),flda_loc(iia,aJ_0H:aJ_1H)
      real*8, allocatable :: flda(:,:),fldo(:,:)
      if(am_i_root()) allocate(flda(iia,jja),fldo(iio,jjo))
      call pack_data(ogrid,fldo_loc,fldo)
      if(am_i_root()) then
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
      endif ! am_i_root
      call unpack_data(agrid,flda,flda_loc)
      if(am_i_root()) deallocate(flda,fldo)
      return
      end subroutine ssto2a

      subroutine ssto2a_global(fldo,flda)
c --- mapping sst from ogcm A grid to agcm A grid
c     input: fldo, output: flda
c
      implicit none
      integer n,ia,ja
      real*8 flda(iia,jja),fldo(iio,jjo)
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
      end subroutine ssto2a_global
c
c
      subroutine veca2o(tauxa_loc,tauya_loc,tauxo_loc,tauyo_loc)
c --- mapping vector from agcm A grid to ogcm A grid
c --- input  tauxa/tauya (N/m2): E_/N_ward on agcm A grid
c --- output tauxo/tauyo (N/m2): +i_/+j_ward on ogcm A grid (S_/E_ward in Mercador domain)
c
      implicit none
      integer i,j,l,n
      real*8, dimension(iia,aJ_0H:aJ_1H) :: tauxa_loc,tauya_loc
      real*8, dimension(iio, J_0H: J_1H) :: tauxo_loc,tauyo_loc
      real*8, dimension(:,:), allocatable :: tauxa,tauya,tauxo,tauyo,
     &     sward,eward
      if(am_i_root()) then
        allocate(
     &       tauxa(iia,jja),tauya(iia,jja),
     &       tauxo(iio,jjo),tauyo(iio,jjo),
     &       sward(iio,jjo),eward(iio,jjo)
     &       )
      endif
      call pack_data(agrid,tauxa_loc,tauxa)
      call pack_data(agrid,tauya_loc,tauya)
      if(am_i_root()) then
c
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
      endif ! am_i_root
      call unpack_data(ogrid,tauxo,tauxo_loc)
      call unpack_data(ogrid,tauyo,tauyo_loc)
      if(am_i_root()) then
        deallocate(tauxa,tauya,tauxo,tauyo,sward,eward)
      endif
      return
      end subroutine veca2o
c
      subroutine flxa2o(flda_loc,fldo_loc)
c --- mapping flux-like field from agcm A grid to ogcm A grid
c     input: flda (W/m*m), output: fldo (W/m*m)
c
      implicit none
      integer i,j,l,n
      real*8 flda_loc(iia,aJ_0H:aJ_1H),fldo_loc(iio,J_0H:J_1H)
      real*8, allocatable :: flda(:,:),fldo(:,:)
      if(am_i_root()) allocate(flda(iia,jja),fldo(iio,jjo))
      call pack_data(agrid,flda_loc,flda)
c
      if(am_i_root()) then
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
      endif ! am_i_root
c
      call unpack_data(ogrid,fldo,fldo_loc)
      if(am_i_root()) deallocate(flda,fldo)
      return
      end subroutine flxa2o

      subroutine flxa2o_global(flda,fldo)
c --- mapping flux-like field from agcm A grid to ogcm A grid
c     input: flda (W/m*m), output: fldo (W/m*m)
c
      implicit none
      integer i,j,l,n
      real*8, intent(in)  :: flda(iia,jja)
      real*8, intent(out) :: fldo(iio,jjo)
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
      end subroutine flxa2o_global
c
      subroutine flxo2a(fldo,flda)
c --- mapping flux-like field from ogcm A grid to agcm A grid
c     input: fldo (W/m*m), output: flda (W/m*m)
c
      implicit none
      real*8, intent(in)  :: fldo(iio,jjo)
      real*8, intent(out) :: flda(iia,jja)
      integer n,ia,ja
c
c$OMP PARALLEL DO
      do 8 ja=1,jja
      do 8 ia=1,iia
      flda(ia,ja)=0.
c
      do 9 n=1,nlisto2a_f(ia,ja)
      flda(ia,ja)=flda(ia,ja)+ 
     . fldo(ilisto2a_f(ia,ja,n),jlisto2a_f(ia,ja,n))*wlisto2a_f(ia,ja,n)
 9    continue
 8    continue
c$OMP END PARALLEL DO
c
      return
      end subroutine flxo2a
c
      subroutine veco2a(tauxo_loc,tauyo_loc,tauxa_loc,tauya_loc)
c --- mapping vector from ogcm C grid to agcm A grid
c --- input  tauxo/tauyo (N/m2): +i_/+j_ward on ogcm C grid (S_/E_ward in Mercador domain)
c --- output tauxa/tauya (N/m2): E_/N_ward on agcm A grid
c
      implicit none
      integer i,j,l,n,ia,ja,jb
      real*8, dimension(iia,aJ_0H:aJ_1H) :: tauxa_loc,tauya_loc
      real*8, dimension(iio, J_0H: J_1H) :: tauxo_loc,tauyo_loc
      real*8, dimension(:,:), allocatable :: tauxa,tauya,tauxo,tauyo,
     &     nward,eward
      real*8 sine
      if(am_i_root()) then
        allocate(
     &       tauxa(iia,jja),tauya(iia,jja),
     &       tauxo(iio,jjo),tauyo(iio,jjo),
     &       nward(iio,jjo),eward(iio,jjo)
     &       )
      endif
      call pack_data(ogrid,tauxo_loc,tauxo)
      call pack_data(ogrid,tauyo_loc,tauyo)
      if(am_i_root()) then
c
c$OMP PARALLEL DO
      do 10 j=1,jj
      do 10 i=1,ii
      nward(i,j)=0.
 10   eward(i,j)=0.
c$OMP END PARALLEL DO
c
c --- rotate taux/tauy to n/e orientation at ocean A grid
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
c --- mapping nward/eward from ogcm to agcm, both on A grid
c
c$OMP PARALLEL DO
      do 16 ja=1,jja
      do 16 ia=1,iia
      tauxa(ia,ja)=0.
      tauya(ia,ja)=0.
c
      do 17 n=1,nlisto2a(ia,ja)
      tauxa(ia,ja)=tauxa(ia,ja)+eward(ilisto2a(ia,ja,n)
     .            ,jlisto2a(ia,ja,n))*wlisto2a(ia,ja,n)
 17   tauya(ia,ja)=tauya(ia,ja)+nward(ilisto2a(ia,ja,n)
     .            ,jlisto2a(ia,ja,n))*wlisto2a(ia,ja,n)
 16   continue
c$OMP END PARALLEL DO
c
      endif ! am_i_root
      call unpack_data(agrid,tauxa,tauxa_loc)
      call unpack_data(agrid,tauya,tauya_loc)
      if(am_i_root()) then
        deallocate(tauxa,tauya,tauxo,tauyo,nward,eward)
      endif
      return
      end subroutine veco2a
c
      subroutine tempro2a(fldo_loc,flda_loc)
c --- mapping sqrt(sqrt(temp**4)) from ogcm A grid to agcm A grid
c --- input: fldo in deg C; outout: flda in deg K
c
      implicit none
      integer n,ia,ja
      real*8 fldo_loc(iio,J_0H:J_1H),flda_loc(iia,aJ_0H:aJ_1H)
      real*8, allocatable :: flda(:,:),fldo(:,:)
      if(am_i_root()) allocate(flda(iia,jja),fldo(iio,jjo))
      call pack_data(ogrid,fldo_loc,fldo)
      if(am_i_root()) then
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
      endif ! am_i_root
      call unpack_data(agrid,flda,flda_loc)
      if(am_i_root()) deallocate(flda,fldo)
      return
      end subroutine tempro2a
c
      subroutine cpl_wgt
      implicit none
      integer :: iz,jz
c
#ifdef  ATM4x5_HYCOM2deg
       integer, parameter :: nsize1=10249200, nsize2=2079936
#endif     
#ifdef ATM2x2h_HYCOM2deg
       integer, parameter :: nsize1=20358000, nsize2=3991680
#endif
#ifdef ATM2x2h_HYCOM1deg
       integer, parameter :: nsize1=83034720, nsize2=10005120
#endif
c --- read in all weights
      if (iio*jjo*((nwgta2o*2+1)*4+nwgta2o*8).ne.nsize1 .or.
     .    iia*jja*((nwgto2a*2+1)*4+nwgto2a*8).ne.nsize2) then
        write(lp,'(a,2i12,a,2i12)') 'wrong size in cpler '
     .  ,iio*jjo*((nwgta2o*2+1)*4+nwgta2o*8)
     .  ,iia*jja*((nwgto2a*2+1)*4+nwgto2a*8)
     .  ,' should be ',nsize1,nsize2
        stop ' wrong size in cpler'
      endif
c
      open(21,file=flnma2o,form='unformatted',status='old',
     .  access='direct',recl=nsize1)
      read(21,rec=1) ilista2o,jlista2o,wlista2o,nlista2o
      close(21)
c
      open(22,file=flnma2o_tau,form='unformatted',status='old',
     .  access='direct',recl=nsize1)
      read(22,rec=1) itaua2o,jtaua2o,wtaua2o,ntaua2o
      close(22)
c
      open(23,file=flnmo2a,form='unformatted',status='old',
     .  access='direct',recl=nsize2)
      read(23,rec=1) ilisto2a,jlisto2a,wlisto2a,nlisto2a
      close(23)
c
      open(24,file=flnmcoso,form='unformatted',status='old')
      read(24) iz,jz,coso,sino
      close(24)
      if (iz.ne.iio .or. jz.ne.jjo) then
        write(lp,*) ' iz,jz=',iz,jz
        stop '(wrong iz/jz in cososino.8bin)'
      endif
c
#ifdef ATM2x2h_HYCOM1deg
      open(25,file=flnmo2a_f,form='unformatted',status='old',     ! TNL
     .  access='direct',recl=nsize2)
      read(25,rec=1) ilisto2a_f,jlisto2a_f,wlisto2a_f,nlisto2a_f
      close(25)
c
      open(26,file=flnma2o_s,form='unformatted',status='old',   ! TNL
     .  access='direct',recl=nsize1)
      read(26,rec=1) ilista2o_s,jlista2o_s,wlista2o_s,nlista2o_s
      close(26)
#endif
c
      return
      end subroutine cpl_wgt

      end module hycom_cpler
