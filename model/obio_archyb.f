!#include "hycom_mpi_hacks.h"
      subroutine obio_archyb
c
c --- write archive file for time level n to flnm ( b i n a r y  hycom fmt)
c
      USE MODEL_COM, only :
     *  itime,iyear1,nday,jdendofm,jyear,jmon,jday,jdate,jhour,aMON
     * ,xlabel,lrunid

      USE HYCOM_SCALARS, only : nstep,time,lp,theta,onem
     &     ,thref
      USE HYCOM_DIM_GLOB, only : jj,JDM,kk,ntrcr,ii,idm,kdm,iia,jja
      USE HYCOM_ARRAYS_GLOB, only: tracer,temp,saln,dp
      USE HYCOM_CPLER, only: ssto2a
      USE hycom_atm, only: focean
c
      implicit none
      integer i,j,k,l
c
      integer no,nop,nt
      character flnm*40,intvl*3,title*80
      real*4 array4(iia,jja)
      real*8 array8(iia,jja)
c
      call getdte(Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)
c --- check if ogcm date matches agcm date
      if (nstep.eq.1) then
        write(flnm,'(a3,i4.4,2a)') amon,0,'.obio',xlabel(1:lrunid)
      elseif (abs((itime+1.)/nday-time).gt.1.e-5) then
        write(*,*) 'mismatching archive date in agcm/ogcm=',
     .     (itime+1.)/nday,time
        stop 'mismatching archive date'
      else
        write(flnm,'(a3,i4.4,2a)') amon,Jyear,'.obio',xlabel(1:lrunid)
      endif
c
      if (jdate.lt.100) then
        write (intvl,'(i3.3)') jdate
      else
        stop 'jdate >100'
      endif
c
      nop=12
      write (lp,'(a/9x,a)') 'storing history data in',flnm
c
      open (unit=nop,file=flnm,status='unknown',
     .      form='unformatted')

      do nt=1,ntrcr
      do k=1,kk

        !convert to atmosgrid
        call ssto2a(tracer(:,:,k,nt),array8)
        array4=array8

        do j=1,jja
        do i=1,iia
         if (focean(i,j).le.0.) array4(i,j)=1.e33
        enddo
        enddo

        write(title,'(a,i4,a,i4)')'tracer, nt=',nt,', k=',k
        write (nop)title,array4 

      enddo
      enddo
!      dp, densities,sst, sss, mld, trac1-16,pco2,gasexchflx
 
      return
      end subroutine obio_archyb
