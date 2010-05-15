#include "rundeck_opts.h"
      subroutine obio_sinksettl(vrbos,kmax,errcon,i,j)

      USE obio_dim
      USE obio_incom,only: wsdeth,mgchltouMC
      USE obio_com, only: P_tend,obio_deltat,D_tend,C_tend
     .                   ,obio_P,det,car
     .                   ,dp1d,wsdet,p1d,obio_ws
     .                   ,rhs,zc
#ifdef OBIO_ON_GARYocean
      USE MODEL_COM,  only : nstep=>itime
      USE ODIAG, only: ij_cexp,oij=>oij_loc
      USE OCEAN, only: lmm,dxypo
#else
      USE hycom_scalars, only: nstep,baclin
#endif

      implicit none

      integer :: i,j,k,nt,kmax,kzc
      real    :: trnd, cexp
      logical :: vrbos,errcon


!---------------------------------------------------------------
! --- phyto sinking and detrital settling
!---------------------------------------------------------------

#ifdef OBIO_ON_GARYocean

!the sinking term is given in units (m/hr)*(mgr,chl/m3)
!in order to be converted into mgr,chl/m3/hr as the tendency
!terms are in the phytoplankton equations, 
!we need to multiply by dz of each layer:
!  dz(k  ) * P_tend(k  ) = dz(k  ) * P_tend(k  ) - trnd
!  dz(k+1) * P_tend(k+1) = dz(k+1) * P_tend(k+1) + trnd
!this way we ensure conservation of tracer after vertical adjustment
!the /hr factor is bcz the obio timestep is in hrs.

      !phyto sinking
      do nt = nnut+1,ntyp-nzoo
        do k = 1,kmax
        rhs(k,nt,16) = 0.
        enddo
        do k = 1,kmax-1
         trnd = obio_P(k,nt)*obio_ws(k,nt-nnut)
         P_tend(k  ,nt) = P_tend(k,  nt) - trnd/dp1d(k  )
         P_tend(k+1,nt) = P_tend(k+1,nt) + trnd/dp1d(k+1)

         rhs(k  ,nt,16) = rhs(k  ,nt,16) - trnd/dp1d(k  )
         rhs(k+1,nt,16) = rhs(k+1,nt,16) + trnd/dp1d(k+1)
        enddo  ! k
!let phytoplankton that reaches the bottom, disappear in the sediment
!        k = kmax
!        trnd = obio_P(k,nt)*obio_ws(k,nt-nnut)
!        P_tend(k,nt)   = P_tend(k,nt)   - trnd/dp1d(k)
!        rhs(k,nt,16)= - trnd/dp1d(k)
      enddo ! n

     
      !diagnostic for total carbon export at compensation depth
      !total carbon = sinking phyto + settling C detritus
      !term1: sinking phytoplankton
      !find layer index for zc
      kzc = 1
      do k=kmax+1,1,-1
           if (p1d(k).gt.zc) kzc = k
      enddo
      if (kzc.lt.1) kzc=1
      if (kzc.gt.lmm(i,j)) kzc=lmm(i,j)

      cexp = 0.
      do nt=nnut+1,nnut+nchl
        k = kzc    !dont integrate over depth
        cexp = cexp
     .        + obio_P(k,nt)*obio_ws(k,nt-nnut)/dp1d(k)    ! mgm3/hr 
     .        * mgchltouMC * 1.d-3                         ! -> uM,C/hr=mili-mol,C/m3/hr
     .        * 24.d0 * 365.d0 *dp1d(k)                    ! -> mol,C/m3/hr
     .        * 12.d0 * 1.d-15                             ! -> Pg,C/m2/yr
     .        * dxypo(j)                                   ! -> Pg,C/yr

cdiag write(*,'(a,5i5,4e12.4)')'sinksettl, cexp1:',
cdiag.   nstep,i,j,k,nt,obio_P(k,nt),obio_ws(k,nt-nnut),dp1d(k),cexp
      enddo

      OIJ(I,J,IJ_cexp) = OIJ(I,J,IJ_cexp)  + cexp 

cdiag write(*,'(a,3i5,e12.4)')'sinksettl, cexp1a:',
cdiag.   nstep,i,j,cexp


      !detritus settling
      do nt = 1,ndet
        do k=1,kmax
        rhs(k,nnut+nchl+nzoo+nt,16) = 0.
        enddo
        do k = 1,kmax-1
         trnd = det(k,nt)*wsdet(k,nt)
         D_tend(k  ,nt) = D_tend(k  ,nt) - trnd/dp1d(k  )
         D_tend(k+1,nt) = D_tend(k+1,nt) + trnd/dp1d(k+1)

         rhs(k  ,nnut+nchl+nzoo+nt,16)= rhs(k  ,nnut+nchl+nzoo+nt,16) 
     .                                - trnd/dp1d(k  )
         rhs(k+1,nnut+nchl+nzoo+nt,16)= rhs(k+1,nnut+nchl+nzoo+nt,16) 
     .                                + trnd/dp1d(k+1)

        enddo  ! k
!let detritus that reaches the bottom, disappear in the sediment
!        k = kmax
!        trnd = det(k,nt)*wsdet(k,nt)
!        D_tend(k,nt)   = D_tend(k,nt)   - trnd/dp1d(k)
!        rhs(k,nnut+nchl+nzoo+nt,16)= - trnd/dp1d(k)
      enddo ! nt

      !diagnostic for carbon export at compensation depth
      !term2: settling C detritus contribution
      cexp =  0.
        k = kzc
        nt= 1            !only the for carbon detritus
        cexp = cexp 
     .        + det(k,nt)*wsdet(k,nt)/dp1d(k)      ! micro-grC/lt/hr == mili,grC/m3/hr
     .        * 1.d-3 /12.d0                       ! -> mol,C/m3/hr
     .        * 24.d0 * 365.d0 *dp1d(k)            ! -> mol,C/m2/yr
     .        * 12.d0 * 1.d-15                     ! -> Pgr,C/m2/yr
     .        * dxypo(j)                           ! -> Pg,C/yr

cdiag write(*,'(a,4i5,4e12.4)')'sinksettl, cexp2:',
cdiag.   nstep,i,j,k,det(k,nt),wsdet(k,nt),dp1d(k),cexp

      OIJ(I,J,IJ_cexp) = OIJ(I,J,IJ_cexp) + cexp

cdiag write(*,'(a,3i5,e12.4)')'sinksettl, cexp2a:',
cdiag.   nstep,i,j,cexp



#else     /* HYCOM */
#ifndef noBIO            /****** for all runs except noBIO tests ********/
      if (kmax.le.1) return


       !phyto sinking
       do nt=1,nchl

          do k=kmax+1,2,-1
            !obio_ws is in m/hr
            obio_ws(k,nt)=.5*(obio_ws(k,nt)+obio_ws(k-1,nt))
          end do
          obio_ws(     1,nt)=0.                !  no flux through sea surface


cdiag      if (vrbos) then
cdiag        do k=1,kmax
cdiag        write(*,'(a,6i5,3e12.4)')'befr sinking',
cdiag.       nstep,kmax,nt,i,j,k,obio_P(k,nnut+nt),
cdiag.       p1d(k),obio_ws(k,nt)*baclin/3600.
cdiag        enddo
cdiag      endif

           do k=1,kmax
            rhs(k,nnut+nt,16)=obio_P(k,nnut+nt)
           enddo
           call advc1d(kmax,obio_P(1,nnut+nt),p1d,
     .                 obio_ws(1,nt)*baclin/3600.,
     .                 vrbos,errcon)
           if (errcon) write(*,*)'error in phyto sinking: nt=',nt
     .                        ,i,j

           do k=1,kmax
            rhs(k,nnut+nt,16)=
     ,          (obio_P(k,nnut+nt)-rhs(k,nnut+nt,16))/obio_deltat
           enddo


cdiag      if (vrbos) then
cdiag        do k=1,kmax
cdiag        write(*,*)'aftr sinking',
cdiag.       nstep,kmax,nt,i,j,k,obio_P(k,nnut+nt),
cdiag.       p1d(k),obio_ws(k,nt)
cdiag        enddo
cdiag      endif


       enddo   !nchl
       !detrital settling
       do nt=1,ndet
          do k=kmax+1,2,-1
            !detritus settling rates in m/hr
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
cdiag      if (vrbos) then
cdiag        do k=1,kmax
cdiag        write(*,*)'befr settling: ',
cdiag.       nstep,kmax,nt,i,j,k,wsdet(k,nt),p1d(k),det(k,nt)
cdiag        enddo
cdiag      endif

          !no flux through the sea floor
          do k=1,kmax
             wsdet(k,nt)=min(wsdet(k,nt),p1d(kmax+1)-p1d(k))
          enddo
!need to change that (?) and create an array that will actually
!hold the excess stuff (sediment array) to be used in
!sequastration studies.
!(sediment=stuff(bottom layer)-stuff(layer above)

          do k=1,kmax
           !this is not really rhs, just a temp storage space
           rhs(k,nnut+nchl+nzoo+nt,16)=det(k,nt)
          enddo
          call advc1d(kmax,det(1,nt),p1d,
     .      wsdet(1,nt)*baclin/3600.,vrbos,errcon)
          if (errcon) write(*,*)
     .     'error in detritus component: nt=',nt,i,j,kmax
          do k=1,kmax

           !this is really rhs
           rhs(k,nnut+nchl+nzoo+nt,16)=
     .            (det(k,nt)-rhs(k,nnut+nchl+nzoo+nt,16))/obio_deltat
          enddo

cdiag      if (vrbos) then
cdiag        do k=1,kmax
cdiag        write(*,*)'aftr settling: ',
cdiag.       nstep,kmax,nt,i,j,k,wsdet(k,nt),p1d(k),det(k,nt)
cdiag        enddo
cdiag      endif


       end do  !ndet

#else
       errcon = .false.
#endif
#endif

      end subroutine obio_sinksettl
