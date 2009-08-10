!#include "hycom_mpi_hacks.h"
      subroutine obio_archyb(nn,dpav,temav,salav,th3av,
     .                       dpmxav,oiceav)
c
c --- write archive file for time level n to flnm ( b i n a r y  hycom fmt)
c
      USE MODEL_COM, only :
     *  itime,iyear1,nday,jdendofm,jyear,jmon,jday,jdate,jhour,aMON
     * ,xlabel,lrunid

      USE HYCOM_SCALARS, only : nstep,time,lp,theta,onem
     &     ,thref,baclin
      USE HYCOM_DIM_GLOB, only : jj,jdm,kk,ntrcr,ii,idm,kdm,iia,jja
      USE HYCOM_ARRAYS_GLOB, only: tracer,temp,saln,p,dpmixl
      USE obio_com, only : pCO2_glob,ao_co2flux_glob,pCO2av, 
     .                     ao_co2fluxav,diag_counter,tracav,plevav
c
      implicit none
      real*8, dimension(idm,jdm) :: dpmxav,oiceav
      real*8, dimension(idm,jdm,kdm) :: dpav,temav,salav,th3av
      real factor
      integer i,j,k,l,kn,nn
     
c
      integer no,nop,nt
      character flnm*40,intvl*3,title*80
c
      call getdte(Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)
c --- check if ogcm date matches agcm date
      if (nstep.eq.1) then
        write(flnm,'(a3,i4.4,2a)') amon,0,'.obio',xlabel(1:lrunid)
!     elseif (abs((itime+1.)/nday-time).gt.1.e-5) then
!       write(0,*) 'obio: mismatching archive date in agcm/ogcm=',
!    .     (itime+1.)/nday,time
!       stop 'obio_archyb: mismatching archive date'
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

      factor=baclin/(jdate*86400.)
c
      open (unit=nop,file=flnm,status='unknown',
     .      form='unformatted')

      ! ocean layer depth
      do k=1,kk
        write(title,'(a,i4)')'depth, k=',k
        call write2giss(nop,p(:,:,k)/onem,title)
      enddo

      ! ocean temperature
      do k=1,kk
      kn=k+nn
        write(title,'(a,i4)')'temp, k=',k
        call write2giss(nop,temp(:,:,kn),title)
      enddo

      ! ocean salinity
      do k=1,kk
      kn=k+nn
        write(title,'(a,i4)')'saln, k=',k
        call write2giss(nop,saln(:,:,kn),title)
      enddo

      !mld
        write(title,'(a)')'MLD (m)'
        call write2giss(nop,dpmixl(:,:,2)/onem,title)

      !tracer
      do nt=1,ntrcr
      do k=1,kk
        write(title,'(a,i4,a,i4)')'tracer, nt=',nt,', k=',k
        call write2giss(nop,tracer(:,:,k,nt),title)
      enddo
      enddo

      !pco2
        write(title,'(a)')'pCO2, water'
        call write2giss(nop,pCO2_glob,title)

      !ao co2 flux
        write(title,'(a)')'AO CO2 flux'
        call write2giss(nop,ao_co2flux_glob,title)


! time averaged arrays

      !dpav
      do k=1,kk
        write(title,'(a,i4)')'dpav, k=',k
        call write2giss(nop,dpav(:,:,k)/onem,title)
      enddo
   
      !tempav
      do k=1,kk
        write(title,'(a,i4)')'tempav, k=',k
        call write2giss(nop,temav(:,:,k),title)
      enddo
   
      !salav
      do k=1,kk
        write(title,'(a,i4)')'salav, k=',k
        call write2giss(nop,salav(:,:,k),title)
      enddo
   
      !th3av
      do k=1,kk
        write(title,'(a,i4)')'th3av, k=',k
        call write2giss(nop,th3av(:,:,k),title)
      enddo
   
      !dpmxav
        write(title,'(a)')'dpmxav'
        call write2giss(nop,dpmxav(:,:)/onem,title)
   
      !oiceav
        write(title,'(a)')'oiceav'
        call write2giss(nop,oiceav(:,:),title)
   
      !tracav
      !no need to divide by diag_counter because pressumably
      !this is already done when divide by plevav
      do nt=1,ntrcr
      tracav(:,:,:,nt)=tracav(:,:,:,nt)/plevav(:,:,:)
      do k=1,kk
        write(title,'(a,i4,a,i4)')'tracav, nt=',nt,', k=',k
        call write2giss(nop,tracav(:,:,k,nt),title)
      enddo
      enddo
   
      !pco2av
        pco2av=pco2av/diag_counter
        write(title,'(a)')'pCO2av'
        call write2giss(nop,pCO2av,title)

      !pco2av
        ao_co2fluxav=ao_co2fluxav/diag_counter
        write(title,'(a)')'ao_co2fluxav'
        call write2giss(nop,ao_co2fluxav,title)

!zero out for next diagnostic period
      diag_counter= 0.
      tracav = 0.
      plevav = 0.
      pco2av = 0.
      ao_co2fluxav = 0.

      return
      end subroutine obio_archyb

      subroutine write2giss(nop,array_o,title)

      USE HYCOM_DIM_GLOB, only : iia,jja,iio,jjo
      USE HYCOM_CPLER, only: ssto2a
      USE HYCOM_SCALARS, only : onem
      USE hycom_atm, only: focean

      integer nop
      real*4 array4(iia,jja)
      real*8 array8(iia,jja),array_o(iio,jjo),focean_glob(iia,jja)
      character title*80

        !convert to atmosgrid
        call ssto2a(array_o,array8)
        array4=array8
        do j=1,jja
        do i=1,iia
         if (focean(i,j).le.0.) array4(i,j)=1.e33
        enddo
        enddo
        write (nop)title,array4 

      end subroutine write2giss
