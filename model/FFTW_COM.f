      module fftw_com
!@sum This module provides variables and routines 
!@+   to compute FFTs using the externam FFTW library
!@auth Denis Gueyffier

      implicit none
      include "fftw3.f"
      save
      integer*8 :: plan
      integer, parameter :: N=72
      real *8 :: infftw(N)
      double complex ::  outfftw(N/2 + 1)
      contains


      subroutine fft0
      implicit none
      include "fftw3.f"

c***  Create plan 
      call dfftw_plan_dft_r2c_1d(plan,N,infftw,outfftw,FFTW_ESTIMATE)

      end subroutine fft0
c*

      subroutine fftend
      implicit none
      call dfftw_destroy_plan(plan)
      end subroutine fftend
c*


      subroutine fft(F,A,B)
!@sum Computing FFT using FFTW
!@+   Input/Output is reorder to conform to G. Russell's convention
!@+   Gary Russell samples from 1 (->2 pi/ N)  to N (->2 pi )
!@+   FFTW uses a more widespread convention (same as Matlab for example): 
!@+   sampling from  0 (->0 )  to N-1 (-> 2  pi (1 -1/ N) ) 
      use FFT72
      implicit none
      real*8, intent(in) :: F(1:N)
      real*8, intent(out) :: A(0:N/2),B(0:N/2)          
      real*8 :: Atemp(1:N/2+1),Btemp(1:N/2+1),
     &     Ar(0:N/2),Br(0:N/2) 
      integer :: i, k, K_(N)

c***  Converting input from Gary Russell's convention to FFTW convention
      infftw(2:N)=F(1:N-1)/N
      infftw(1)=F(N)/N

c***  Compute FFT
      call dfftw_execute_(plan)

      write(*,*) 'Afttw',real(outfftw)
      write(*,*) 'Bfttw',imag(outfftw)

c***  Converting back output from FFTW convention to Gary Russell's convention 
      Atemp(:)=2.d0*real(outfftw)
      Atemp(1)=Atemp(1)/2.d0
      Atemp(N/2+1)=Atemp(N/2+1)/2.d0

      Btemp(:)=-2.d0*imag(outfftw)
      Btemp(1)=Btemp(1)/2.d0
      Btemp(N/2+1)=Btemp(N/2+1)/2.d0

      A(0:N/2)=Atemp(1:N/2+1)
      B(0:N/2)=Btemp(1:N/2+1)

      end subroutine fft

      end module fftw_com
