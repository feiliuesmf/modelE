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
      USE OCEAN, only: dxypo,lmm
#else
      USE hycom_scalars, only: nstep,baclin
#endif

      implicit none

      integer :: i,j,k,nt,kmax,kzc
      real    :: trnd
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
      do k=kmax+1,1,-1
           if (p1d(k).gt.zc) kzc = k
      enddo
      if (kzc.lt.1) kzc=1
      if (kzc.gt.lmm(i,j)) kzc=lmm(i,j)

      do nt=nnut+1,nnut+nchl
         OIJ(I,J,IJ_cexp) = OIJ(I,J,IJ_cexp) 
                          !mgm3/hr -> uM/hr=mili,molC/m3/hr-> mol,C/m3/hr
     .                    + rhs(kzc,nt,16) * mgchltouMC * 1.d-3    
                          ! -> Pgr,C/yr
     .                    * 24.d0 * 365.d0 *dxypo(J)*75.d0
     .                    * 12.d0 * 1.d-15
      enddo

cdiag write(*,'(a,3i5,8e12.4)')'carbon export1:',
cdiag.   nstep,i,j,obio_P(kzc,5),obio_ws(kzc,1),
cdiag.             obio_P(kzc,6),obio_ws(kzc,2),
cdiag.             obio_P(kzc,7),obio_ws(kzc,3),
cdiag.             obio_P(kzc,8),obio_ws(kzc,4) 
cdiag write(*,'(a,3i5,4e12.4)')'carbon export1a:',
cdiag.   nstep,i,j,       rhs(kzc,5,16),
cdiag.                    rhs(kzc,6,16),
cdiag.                    rhs(kzc,7,16),
cdiag.                    rhs(kzc,8,16)

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
      OIJ(I,J,IJ_cexp) = OIJ(I,J,IJ_cexp) 
                       !micro-grC/lt/hr -> mili,grC/m3/hr -> mol,C/m3/hr
     .                 + rhs(kzc,nnut+nchl+nzoo+1,16)*1.d-3 /12.d0   
                          ! -> Pgr,C/yr
     .                    * 24.d0 * 365.d0 *dxypo(J)*75.d0
     .                    * 12.d0 * 1.d-15

cdiag write(*,'(a,4i5,5e12.4)')'carbon export2:',
cdiag.   nstep,i,j,kzc,zc,p1d(kzc),det(kzc,1),wsdet(kzc,1),
cdiag.                    rhs(kzc,nnut+nchl+nzoo+1,16)

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
cdiag        write(*,*)'befr sinking',
cdiag.       nstep,kmax,nt,i,j,k,obio_P(k,nnut+nt),
cdiag.       p1d(k),obio_ws(k,nt)
cdiag        enddo
cdiag      endif

           do k=1,kmax
            rhs(k,nnut+nt,16)=obio_P(k,nnut+nt)
           enddo
           call advc1d(kmax,obio_P(1,nnut+nt),p1d,
     .                 obio_ws(1,nt)*baclin/3600.,
     .                 vrbos,errcon)
           if (errcon) write(*,*)'error in phyto sinking: nt=',nt

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
