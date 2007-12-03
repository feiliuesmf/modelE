      subroutine obio_update(vrbos,kmax,errcon)
 
c  Performs updating of biological particles, uses mid-point
c  leap frog method.

      USE obio_dim
      USE obio_incom,only: wsdeth,obio_wsh
      USE obio_com, only: P_tend,obio_deltat,D_tend,C_tend
     .                   ,obio_P,det,car
     .                   ,dp1d,wsdet,p1d,obio_ws
     .                   ,rhs


      implicit none
 
#include "dimensions.h"
#include "dimension2.h" 
#include "common_blocks.h"

      integer :: nt,kmax
      real    :: Pnew,Dnew,Cnew
      logical :: vrbos,errcon
 
c  Loop to update
c   indexes are mixed because new H has already been computed
c   in update.F, but P has not been updated yet

      do 1000 k = 1,kmax

        do nt = 1,ntyp+n_inert
         Pnew = (obio_P(k ,nt) +  P_tend(k,nt)*obio_deltat)
         obio_P(k,nt) = Pnew
        enddo
 
        do nt = 1,ndet
         Dnew = (det(k,nt)     +  D_tend(k,nt)*obio_deltat)
         det(k,nt) = Dnew
        enddo

        do nt = 1,ncar
         Cnew = (car(k,nt) +  C_tend(k,nt)*obio_deltat)
         car(k,nt) = Cnew
        enddo

 1000 continue

!!!   return
!---------------------------------------------------------------
! --- phyto sinking and detrital settling
!---------------------------------------------------------------
      if (kmax.le.1) return

   
       !phyto sinking
       do nt=1,nchl

          do k=kmax+1,2,-1
            !need to multiply by timestep in hrs
            !(which we do not do here becz timestep=1hr)
            obio_ws(k,nt)=.5*(obio_ws(k,nt)+obio_ws(k-1,nt))
          end do
          obio_ws(     1,nt)=0.                !  no flux through sea surface


cdiag      if (vrbos) then
cdiag        do k=1,kmax
cdiag        write(*,'(a,3i5,3e12.4)')'befr sinking',
cdiag.       kmax,nt,k,obio_P(k,nnut+nt),p1d(k),obio_ws(k,nt)
cdiag        enddo
cdiag      endif

           do k=1,kmax
            rhs(k,nnut+nt,16)=obio_P(k,nnut+nt)
           enddo
           call advc1d(kmax,obio_P(1,nnut+nt),p1d,obio_ws(1,nt),
     .                 vrbos,errcon)
           if (errcon) write(*,*)'error in phyto sinking: nt=',nt

           do k=1,kmax
            rhs(k,nnut+nt,16)=
     ,          (obio_P(k,nnut+nt)-rhs(k,nnut+nt,16))/obio_deltat
           enddo


cdiag      if (vrbos) then
cdiag        do k=1,kmax
cdiag        write(*,'(a,3i5,3e12.4)')'aftr sinking',
cdiag.       kmax,nt,k,obio_P(k,nnut+nt),p1d(k),obio_ws(k,nt)
cdiag        enddo
cdiag      endif


       enddo


       !detrital settling
       do nt=1,ndet
          do k=kmax+1,2,-1
            !detritus
            !need to multiply by timestep in hrs
            !(which we do not do here becz timestep=1hr)
            wsdet(k,nt)=.5*(wsdet(k,nt)+wsdet(k-1,nt))
            !for the moment let wsdet be constant and not depending on T
            !when I tried to change this I would get crossover errors
            !for layers of near-zero thickness
            wsdet(k,nt)=wsdeth(nt)

          end do
          wsdet(     1,nt)=0.                !  no flux through sea surface

cdiag     if (vrbos.and.nt.eq.1) then
cdiag       do k=1,kdm
cdiag          write(lp,'(a,3i5,3e12.4)')
cdiag.         'bfre bcdond: ',nt,kmax,k,wsdet(k,nt),p1d(k),det(k,nt)
cdiag       enddo
cdiag       write(400,'(a,3i5,2e12.4)')
cdiag.     'bfre bcdond: ',nstep,kmax,k,wsdet(k,nt),p1d(k)
cdiag     endif

          !no flux through the sea floor
!         do k=1,kmax
!            wsdet(k,nt)=min(wsdet(k,nt),p1d(kmax+1)-p1d(k))
!         enddo
!for the moment let flux go through the sea floor.
!need to change that and create an array that will actually 
!hold the excess stuff (sediment array) to be used in 
!sequastration studies. 
!(sediment=stuff(bottom layer)-stuff(layer above)

cdiag     if (vrbos) then
cdiag       do k=1,kdm
cdiag          write(400,'(a,3i5,2e12.4)')
cdiag.         'bfre advc1d: ',nt,kmax,k,wsdet(k,nt),det(k,nt)
cdiag       enddo
cdiag     endif

          do k=1,kmax
           rhs(k,nnut+nchl+nzoo+nt,16)=det(k,nt)
          enddo
          call advc1d(kmax,det(1,nt),p1d,wsdet(1,nt),vrbos,errcon)
          if (errcon) write(*,*)'error in detritus component: nt=',nt
          do k=1,kmax
           rhs(k,nnut+nchl+nzoo+nt,16)=
     .            (det(k,nt)-rhs(k,nnut+nchl+nzoo+nt,16))/obio_deltat
          enddo

cdiag     if (vrbos) then
cdiag       do k=1,kdm
cdiag          write(400,'(a,3i5,2e12.4)')
cdiag.         'aftr advc1d: ',nt,kmax,k,wsdet(k,nt),det(k,nt)
cdiag       enddo
cdiag     endif

       end do


      return
      end
