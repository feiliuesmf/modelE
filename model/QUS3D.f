#include "rundeck_opts.h"
      MODULE TRACER_ADV
!@sum MODULE TRACER_ADV arrays needed for tracer advection

      USE MODEL_COM, ONLY : IM,JM,LM,byim
      SAVE
      INTEGER, PARAMETER :: ncmax=10
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NSTEPX
      INTEGER :: NCYC
      real*8, dimension(:,:), allocatable :: pv_south
      integer, dimension(:), allocatable :: ncycxy

C**** zonal mean diags
      REAL*8,  ALLOCATABLE, DIMENSION(:,:)   :: sfbm,sbm,sbf,
     *                                          sfcm,scm,scf
C**** vertically integrated horizontal fluxes
      REAL*8,  ALLOCATABLE, DIMENSION(:,:)   :: safv,sbfv

      real*8, parameter :: mrat_lim=0.375

      contains

      SUBROUTINE AADVQ(RM,RMOM,qlimit,tname)
      USE DOMAIN_DECOMP, only :
     &     grid, get, HALO_UPDATE,HALO_UPDATE_COLUMN, NORTH,SOUTH
      USE QUSDEF
      USE QUSCOM, ONLY : IM,JM,LM
      USE DYNAMICS, ONLY: pu=>pua, pv=>pva, sd=>sda, mb, ma
      IMPLICIT NONE
      character*8 tname          !tracer name
      logical :: qlimit
      REAL*8, dimension(im,grid%J_STRT_HALO:grid%J_STOP_HALO,lm) :: 
     &                  rm
      REAL*8, dimension(NMOM,IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) 
     &               :: rmom
      REAL*8, dimension(im,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     mflx,mwdn,fdn,fdn0
      REAL*8, dimension(nmom,im,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     fmomdn

      INTEGER :: I,J,L,N,nc,ncxy
      integer :: jmin_x,jmax_x
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0H, J_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)


      DO L=1,LM
         MA(:,:,L) = MB(:,:,L) ! fill in halo lats
      ENDDO

      do nc=1,ncyc

        if(nc.gt.1) CALL HALO_UPDATE(grid, ma, FROM=NORTH+SOUTH)
        CALL HALO_UPDATE(grid, rm,    FROM=NORTH+SOUTH)
        CALL HALO_UPDATE_COLUMN(grid, rmom, FROM=NORTH+SOUTH)

        do j=j_0,j_1
          do i=1,im
            mwdn(i,j) = 0.
            fdn(i,j) = 0.
            fdn0(i,j) = 0.
            fmomdn(:,i,j) = 0.
          enddo
        enddo

        DO L=1,LM+1             ! L=LM+1 is for the last z pass

          if(l.le.lm) then

            do ncxy=1,ncycxy(l) ! loop over horizontal cycles

              if(ncxy.gt.1) then ! update boundaries if necessary
                CALL HALO_UPDATE(grid, ma(:,:,l), FROM=NORTH+SOUTH)
                CALL HALO_UPDATE(grid, rm(:,:,l), FROM=NORTH+SOUTH)
                CALL HALO_UPDATE_COLUMN(
     &                           grid, rmom(:,:,:,l), FROM=NORTH+SOUTH)
              endif

C**** Fill in values at the poles
C**** SOUTH POLE:
              if (HAVE_SOUTH_POLE) then
                DO I=2,IM
                  RM(I,1 ,L) =   RM(1,1 ,L)
                  DO N=1,NMOM
                    RMOM(N,I,1 ,L) =  RMOM(N,1,1 ,L)
                  enddo
                enddo
              endif             !SOUTH POLE

c**** NORTH POLE:
              if (HAVE_NORTH_POLE) then
                DO I=2,IM
                  RM(I,JM,L) =   RM(1,JM,L)
                  DO N=1,NMOM
                    RMOM(N,I,JM,L) =  RMOM(N,1,JM,L)
                  enddo
                enddo
              endif             ! NORTH POLE

              jmin_x = j_0s; jmax_x = j_1s
              if(qlimit) then
c when flow out both sides would cause negative tracer mass, modify moments
                if(.not.HAVE_SOUTH_POLE) then
                  j=j_0h
                  do i=1,im
                    if(pv(i,j,l).gt.0. .and. pv_south(i,l).lt.0.) then
                      call checkflux(pv_south(i,l),pv(i,j,l),ma(i,j,l),
     &                     rm(i,j,l),rmom(my,i,j,l),rmom(myy,i,j,l))
                    endif
                  enddo
                endif
                do j=max(2,j_0),min(jm-1,j_1h)
                  do i=1,im
                    if(pv(i,j,l).gt.0. .and. pv(i,j-1,l).lt.0.) then
                      call checkflux(pv(i,j-1,l),pv(i,j,l),ma(i,j,l),
     &                     rm(i,j,l),rmom(my,i,j,l),rmom(myy,i,j,l))
                    endif
                  enddo
                enddo
                CALL AADVQY(RM(1,j_0h,l),RMOM(1,1,j_0h,l),
     &               MA(1,j_0h,l),pv(1,j_0h,l)
     &               ,sbf(j_0h,l),sbm(j_0h,l),sfbm(j_0h,l),sbfv)
                CALL AADVQX(RM(1,j_0h,l),RMOM(1,1,j_0h,l),
     &               MA(1,j_0h,l),pu(1,j_0h,l),jmin_x,jmax_x,
     &               nstepx(j_0h,l),safv)
              else
                CALL AADVQY2(RM(1,j_0h,l),RMOM(1,1,j_0h,l),
     &               MA(1,j_0h,l),pv(1,j_0h,l)
     &               ,sbf(j_0h,l),sbm(j_0h,l),sfbm(j_0h,l),sbfv)
                CALL AADVQX2(RM(1,j_0h,l),RMOM(1,1,j_0h,l),
     &               MA(1,j_0h,l),pu(1,j_0h,l),jmin_x,jmax_x,
     &               nstepx(j_0h,l),safv)
              endif

            enddo               ! ncxy
          endif                 ! l.le.lm

c when flow out both sides would cause negative tracer mass, modify moments
          if(qlimit .and. l.gt.1 .and. l.lt.lm) then
            do j=j_0,j_1
              do i=1,im
                if(sd(i,j,l).gt.0. .and. sd(i,j,l-1).lt.0.) then
                  call checkflux(sd(i,j,l-1),sd(i,j,l),ma(i,j,l),
     &                 rm(i,j,l),rmom(mz,i,j,l),rmom(mzz,i,j,l))
                endif
              enddo
            enddo
          endif
          if(l.gt.1) then
            if(qlimit) then
              CALL AADVQZ(RM(1,j_0h,l-1),RMOM(1,1,j_0h,l-1),
     &             MA(1,j_0h,l-1),SD(1,j_0h,L-1),mwdn,fdn,fmomdn,fdn0
     &             ,scf(j_0h,l-1),scm(j_0h,l-1),sfcm(j_0h,l-1))
            else
              CALL AADVQZ2(RM(1,j_0h,l-1),RMOM(1,1,j_0h,l-1),
     &             MA(1,j_0h,l-1),SD(1,j_0h,L-1),mwdn,fdn,fmomdn!,fdn0
     &             ,scf(j_0h,l-1),scm(j_0h,l-1),sfcm(j_0h,l-1))
            endif
          endif

        ENDDO                   ! l
      enddo                     ! ncyc

      RETURN
      END SUBROUTINE AADVQ

      SUBROUTINE AADVQ0(dt_dummy)
!@sum AADVQ0 initialises advection of tracer.
!@+   Decide how many cycles to take such that mass does not become
!@+   too small during any of the operator splitting steps of each cycle
!@auth Maxwell Kelley
      USE DYNAMICS, ONLY: mu=>pua, mv=>pva, mw=>sda, mb, ma
      USE DOMAIN_DECOMP, ONLY : GRID, GET, GLOBALSUM, HALO_UPDATE
     &     ,halo_updatej, globalmax
      USE DOMAIN_DECOMP, ONLY : NORTH, SOUTH, AM_I_ROOT
      USE QUSCOM, ONLY : IM,JM,LM
      IMPLICIT NONE
      real*8, intent(in) :: dt_dummy
      INTEGER :: i,j,l,n,nc,nbad,nbad_loc,ierr_loc,ierr,nc3d,ncxy
     &     ,nstepj_dum,ncycxy_loc
      REAL*8 :: byn,ssp,snp,mvbyn,mwbyn,mpol,byn3d
      real*8, dimension(im) :: mubyn,am,mi
      integer, dimension(grid%j_strt_halo:grid%j_stop_halo) :: nstep_dum
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) ::
     &     ma2d,mb2d

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0H, J_1H
      INTEGER :: J_0S, J_1S, J_0STG, J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_HALO=J_0H,   J_STOP_HALO=J_1H,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      DO L=1,LM
        IF (HAVE_SOUTH_POLE) MU(:,1,L) = 0.
        PV_SOUTH(:,L) = MV(:,J_0H,L)
        DO J=J_0H,J_1S
          MV(:,J,L) = MV(:,J+1,L)
        END DO
        IF (HAVE_NORTH_POLE) MU(:,JM,L) = 0.
        IF (HAVE_NORTH_POLE) MV(:,JM,L) = 0.
        MW(:,j_0:j_1,L) = -MW(:,j_0:j_1,L) ! accum phase should do the switch
      ENDDO
      MW(:,j_0:j_1,LM) = 0.
      CALL HALO_UPDATE(grid, MV, FROM=NORTH)

c
c first determine the horizontal-vertical ncyc
c
      nbad = 1
      ncyc = 0
      do while(nbad.gt.0)
        ncyc = ncyc + 1
        if(ncyc.gt.ncmax) then
          if (AM_I_ROOT()) then
            write(6,*) 'stop: ncyc>ncmax in AADVQ0'
          end if
          call stop_model('AADVQ0: ncyc>ncmax',255)
        end if
        byn = 1./ncyc
        nbad_loc = 0
        do l=1,lm
          ma(:,j_0h:j_1h,l) = mb(:,j_0h:j_1h,l)
        enddo
        do nc=1,ncyc

C**** check whether horizontal fluxes reduce mass too much
          do l=1,lm
            do j=J_0S,J_1S
              i=1
              ma(i,j,l) = ma(i,j,l) + 
     &             byn*(mu(im ,j,l)-mu(i,j,l)+mv(i,j-1,l)-mv(i,j,l))
              if (ma(i,j,l).lt.mrat_lim*mb(i,j,l)) then
                nbad_loc = nbad_loc + 1
              endif
              do i=2,im
                ma(i,j,l) = ma(i,j,l) +
     &               byn*(mu(i-1,j,l)-mu(i,j,l)+mv(i,j-1,l)-mv(i,j,l))
                if (ma(i,j,l).lt.mrat_lim*mb(i,j,l)) then
                  nbad_loc = nbad_loc + 1
                endif
              enddo
            enddo
            if (HAVE_SOUTH_POLE) then
              ssp = sum(ma(:, 1,l)-mv(:,   1,l)*byn)*byim
              ma(:,1 ,l) = ssp
              if (ma(1,1,l).lt.mrat_lim*mb(1,1,l)) then
                nbad_loc = nbad_loc + 1
              endif
            endif
            if (HAVE_NORTH_POLE) then
              snp = sum(ma(:,jm,l)+mv(:,jm-1,l)*byn)*byim
              ma(:,jm,l) = snp
              if (ma(1,jm,l).lt.mrat_lim*mb(1,jm,l)) then
                nbad_loc = nbad_loc + 1
              endif
            endif
          end do                ! l
          CALL GLOBALSUM(grid, nbad_loc, nbad, all=.true.)
          IF(NBAD.GT.0) exit    ! nc loop

c check courant numbers in the z direction
          do l=1,lm
            if(l.lt.lm) then
              do j=j_0,j_1
                do i=1,im
                  mwbyn = mw(i,j,l)*byn
                  if((ma(i,j,l)-mwbyn)*(ma(i,j,l+1)+mwbyn).lt.0.) then
                    nbad_loc = nbad_loc + 1
                  endif
                enddo
              enddo
            endif
c update mass from z fluxes
            if(nc.lt.ncyc) then
              if(l.eq.1) then   ! lowest layer
                do j=J_0,J_1
                  do i=1,im
                    ma(i,j,l) = ma(i,j,l)-mw(i,j,l)*byn
                  enddo
                enddo
              elseif(l.eq.lm) then ! topmost layer
                do j=J_0,J_1
                  do i=1,im
                    ma(i,j,l) = ma(i,j,l)+mw(i,j,l-1)*byn
                  enddo
                enddo
              else              ! interior layers
                do j=J_0,J_1
                  do i=1,im
                    ma(i,j,l) = ma(i,j,l)+(mw(i,j,l-1)-mw(i,j,l))*byn
                  end do
                end do
              endif
            endif               ! nc.lt.ncyc
          enddo                 ! l
          CALL GLOBALSUM(grid, nbad_loc, nbad, all=.true.)
          IF(NBAD.GT.0) exit    ! nc loop

        enddo                   ! nc loop
      enddo                     ! nbad .gt. 0

C**** Divide the mass fluxes by the number of 3D cycles
      if(ncyc.gt.1) then
        byn = 1./ncyc
        pv_south(:,:) = pv_south(:,:)*byn
        DO L=1,LM
          mu(:,:,l)=mu(:,:,l)*byn
          mv(:,:,l)=mv(:,:,l)*byn
          mw(:,:,l)=mw(:,:,l)*byn
        ENDDO
      endif

c
c now determine the ncycxy for each level, and nstepx for each lat/level
c
      ierr_loc = 0  ! for xstep errors
      do l=1,lm
        if(ncyc.gt.1) mb2d(:,j_0h:j_1h) = mb(:,j_0h:j_1h,l)
        nstepx(j_0:j_1,l) = 1
        ncycxy(l) = 1
        nc3dloop: do nc3d=1,ncyc
          nbad = 1
          do while(nbad.gt.0)
            if(ncycxy(l).gt.ncmax) exit nc3dloop
            if(nc3d.eq.1) then
              ma2d(:,j_0h:j_1h) = mb(:,j_0h:j_1h,l)
            else
              ma2d(:,j_0:j_1) = mb2d(:,j_0:j_1)
              CALL HALO_UPDATE(grid, MA2D, FROM=SOUTH+NORTH) ! rarely needed
            endif
            byn = 1./ncycxy(l)
            nbad = 0
            nstep_dum(j_0:j_1) = 1
            do nc=1,ncycxy(l)

c check y direction courant numbers
              do j=max(2,j_0),min(j_1,jm-2)
                do i=1,im
                  mvbyn = mv(i,j,l)*byn
                  if((ma2d(i,j)-mvbyn)*(ma2d(i,j+1)+mvbyn).lt.0.) then
                    nbad = nbad + 1
                  endif
                enddo
              enddo
              if(have_south_pole) then
                j=1
                mpol = ma2d(1,j)*im
                do i=1,im
                  mvbyn = mv(i,j,l)*byn
                  if((mpol-mvbyn)*(ma2d(i,j+1)+mvbyn).lt.0.) then
                    nbad = nbad + 1
                  endif
                enddo
              endif
              if(have_north_pole) then
                j=jm-1
                mpol = ma2d(1,j+1)*im
                do i=1,im
                  mvbyn = mv(i,j,l)*byn
                  if((ma2d(i,j)-mvbyn)*(mpol+mvbyn).lt.0.) then
                    nbad = nbad + 1
                  endif
                enddo
              endif
c check mass ratios after y direction
              do j=J_0S,J_1S
                do i=1,im
                  ma2d(i,j) = ma2d(i,j) + (mv(i,j-1,l)-mv(i,j,l))*byn
                  if (ma2d(i,j).lt.mrat_lim*mb(i,j,l)) then
                    nbad = nbad + 1
                  endif
                end do
              end do
              if (HAVE_SOUTH_POLE) then
                ssp = sum(ma2d(:, 1)-mv(:,   1,l)*byn)*byim
                ma2d(:,1 ) = ssp
                if (ma2d(1,1).lt.mrat_lim*mb(1,1,l)) then
                  nbad = nbad + 1
                endif
              endif
              if (HAVE_NORTH_POLE) then
                snp = sum(ma2d(:,jm)+mv(:,jm-1,l)*byn)*byim
                ma2d(:,jm) = snp
                if (ma2d(1,jm).lt.mrat_lim*mb(1,jm,l)) then
                  nbad = nbad + 1
                endif
              endif
              if(NBAD.GT.0) then
                ncycxy(l) = ncycxy(l) + 1
                exit            ! nc loop
              endif

C**** determine nstepx
              do j=J_0S,J_1S
                mubyn(:) = mu(:,j,l)*byn
                call XSTEP(j,l,ierr_loc,
     &               ma2d(1,j),mubyn,nstepj_dum,am,mi)
                nstep_dum(j) = max(nstep_dum(j),nstepj_dum)
                if(nc.lt.ncycxy(l)) ma2d(:,j) = mi(:)
              enddo

c update boundary airmasses to avoid layerwise mpi communication during iteration
              if(nc.lt.ncycxy(l)) then
                if(.not.have_south_pole) then
                  j=j_0h
                  i=1
                  ma2d(i,j) = ma2d(i,j) + byn*
     &                 (mu(im,j,l)-mu(i,j,l)+pv_south(i,l)-mv(i,j,l))
                  do i=2,im
                    ma2d(i,j) = ma2d(i,j) + byn*
     &                   (mu(i-1,j,l)-mu(i,j,l)+pv_south(i,l)-mv(i,j,l))
                  enddo
                endif
                if(.not.have_north_pole) then
                  j=j_1h
                  i=1
                  ma2d(i,j) = ma2d(i,j) + byn*
     &                 (mu(im,j,l)-mu(i,j,l)+mv(i,j-1,l)-mv(i,j,l))
                  do i=2,im
                    ma2d(i,j) = ma2d(i,j) + byn*
     &                   (mu(i-1,j,l)-mu(i,j,l)+mv(i,j-1,l)-mv(i,j,l))
                  enddo
                endif
              endif
            
            enddo               ! nc loop
          enddo                 ! nbad.gt.0

c now add the z mass tendency at this level
          if(nc3d.lt.ncyc) then
            if(l.eq.1) then     ! lowest layer
              do j=J_0,J_1
                do i=1,im
                  mb2d(i,j) = ma2d(i,j)-mw(i,j,l)
                enddo
              enddo
            else if(l.eq.lm) then ! topmost layer
              do j=J_0,J_1
                do i=1,im
                  mb2d(i,j) = ma2d(i,j)+mw(i,j,l-1)
                end do
              end do
            else                ! interior layers
              do j=J_0,J_1
                do i=1,im
                  mb2d(i,j) = ma2d(i,j)+(mw(i,j,l-1)-mw(i,j,l))
                enddo
              enddo
            endif
          endif

          nstepx(j_0:j_1,l) = max(nstepx(j_0:j_1,l),nstep_dum(j_0:j_1))
          
        enddo nc3dloop

c globalmax of ncycxy. DOMAIN_DECOMP needs a vector version of globalmax
        ncycxy_loc = ncycxy(l)
        CALL GLOBALMAX(grid, ncycxy_loc, ncycxy(l))
        if(ncycxy(l).gt.ncmax) then
          if (AM_I_ROOT()) then
            write(6,*) 'stop: ncycxy>ncmax in AADVQ0'
          endif
          call stop_model('AADVQ0: ncycxy>ncmax',255)
        endif

C**** Further divide the xy mass fluxes by the number of xy cycles at each level
        if(ncycxy(l).gt.1) then
          byn = 1./ncycxy(l)
          pv_south(:,l) = pv_south(:,l)*byn
          mu(:,:,l)=mu(:,:,l)*byn
          mv(:,:,l)=mv(:,:,l)*byn
        endif

      enddo ! l

      CALL GLOBALSUM(grid, ierr_loc, ierr, all=.true.)
      IF(ierr.GT.0)
     &     call stop_model('too many steps in xstep',255)

      if(AM_I_ROOT()) then
        if(ncyc.gt.1) write(6,*) 'AADVQ0: ncyc>1',ncyc
        do l=1,lm
          if(ncycxy(l).gt.1)
     &         write(6,*) 'AADVQ0: ncycxy>1 at l= ',l,ncycxy(l)
        enddo
      endif
      
      RETURN
      END subroutine AADVQ0

      end MODULE TRACER_ADV

      SUBROUTINE XSTEP(jprt,lprt,ierr,M,MU,NSTEP,am,mi)
!@sum XSTEP determines the number of X timesteps for tracer dynamics
!@+    using Courant limits
!@auth J. Lerner and M. Kelley
!@ver  1.0
      USE QUSCOM, ONLY : IM
      IMPLICIT NONE
      integer :: jprt,lprt,ierr
      REAL*8, dimension(im) :: m,mu
      REAL*8, dimension(im) :: am ! workspace
      REAL*8, dimension(im) :: mi ! output mass
      integer :: nstep
      integer :: i,ns
      REAL*8 :: courmax
      integer, parameter :: nstepmax=20 ! 60 needed for 0.5 degree resolution
c      ierr=0
      nstep=0
      courmax = 2.
      do while(courmax.gt.1.)
        nstep = nstep+1   !(1+int(courmax))
        if(nstep.eq.nstepmax)  then
          write(6,*) 'xstep: j,l,nstep,courmax=',
     &         jprt,lprt,nstep,courmax
          ierr=1
          return
        endif
        do i=1,im
          am(i) = mu(i)/nstep
          mi(i) = m(i)
        enddo
        courmax = 0.
        do ns=1,nstep
          i=1
          if(am(i).gt.0.) then
            courmax = max(courmax,+am(i)/mi(i))
          else
            courmax = max(courmax,-am(i)/mi(i+1))
          endif
          i=im
          if(am(i).gt.0.) then
            courmax = max(courmax,+am(i)/mi(i))
          else
            courmax = max(courmax,-am(i)/mi(1))
          endif
          do i=2,im-1
            if(am(i).gt.0.) then
              courmax = max(courmax,+am(i)/mi(i))
            else
              courmax = max(courmax,-am(i)/mi(i+1))
            endif
            mi(i) = mi(i)+am(i-1)-am(i)
          enddo
          i=1
          mi(i) = mi(i)+am(im)-am(i)
          i=im
          mi(i) = mi(i)+am(i-1)-am(i)
        enddo ! ns=1,nstep
      enddo      ! while(courmax.gt.1.)
      RETURN
      END SUBROUTINE XSTEP

      SUBROUTINE ALLOC_TRACER_ADV(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
!@ver  1.0
      USE TRACER_ADV
      USE DOMAIN_DECOMP, ONLY : DIST_GRID, GET
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H
      INTEGER :: IER

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

      ALLOCATE(NSTEPX(J_0H:J_1H,LM))

      ALLOCATE( PV_SOUTH(IM,LM) )
      ALLOCATE( NCYCXY(LM) )

      ALLOCATE( sfbm(J_0H:J_1H,LM),
     *           sbm(J_0H:J_1H,LM),
     *           sbf(J_0H:J_1H,LM),
     *          sfcm(J_0H:J_1H,LM),
     *           scm(J_0H:J_1H,LM),
     *           scf(J_0H:J_1H,LM) )

      ALLOCATE( safv(IM,J_0H:J_1H),
     *          sbfv(IM,J_0H:J_1H) )
  
      END SUBROUTINE ALLOC_TRACER_ADV

      subroutine aadvqx(rm,rmom,mass,mu,jmin,jmax,nstep,safv)
!@sum  AADVQX advection driver for x-direction
!@auth Maxwell Kelley
      use DOMAIN_DECOMP, only : grid, GET
      use QUSDEF
      use QUSCOM, only : im,jm
      implicit none
      integer :: jmin,jmax
      REAL*8, dimension(im,grid%J_STRT_HALO:grid%J_STOP_HALO) :: 
     &                  rm,mass,mu
      REAL*8, dimension(NMOM,IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &                  rmom
      integer, dimension(grid%J_STRT_HALO:grid%J_STOP_HALO) :: nstep
      REAL*8, dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: safv
      integer :: i,ii,j,ns
      real*8 :: frac1,fracm,fw,fe,feim,amw,dm2,mold,mnew,bymnew
      real*8 :: fe0,fe_pass,feim0,fw0,rm0,fex_pass,fexx_pass
      real*8, dimension(nmom) :: fmomw,fmome,fmomeim
      real*8, dimension(im) :: am

c**** Get useful local parameters for domain decomposition
      integer :: J_0, J_1, J_0S, J_1S
      CALL GET(grid, J_STRT = J_0 , J_STOP=J_1,
     &             J_STRT_SKP=J_0S,J_STOP_SKP=J_1S )

      do j=jmin,jmax
      am(:) = mu(:,j)/nstep(j)
c****
c**** loop over timesteps
c****
      do ns=1,nstep(j)

c when flow out both sides would cause negative tracer mass, modify moments
      i=1
        if(am(i).gt.0. .and. am(im).lt.0.) then
          call checkflux(am(im),am(i),mass(i,j),
     &         rm(i,j),rmom(mx,i,j),rmom(mxx,i,j))
        endif
      do i=2,im
        if(am(i).gt.0. .and. am(i-1).lt.0.) then
          call checkflux(am(i-1),am(i),mass(i,j),
     &         rm(i,j),rmom(mx,i,j),rmom(mxx,i,j))
        endif
      enddo

      i=im
      if(am(i).lt.0.) then  ! air mass flux is negative
        ii=1
        frac1=+1.
      else                      ! air mass flux is positive
        ii=i
        frac1=-1.
      endif
      fracm=am(i)/mass(ii,j)
      frac1=fracm+frac1
      fe=fracm*(rm(ii,j)-frac1*(rmom(mx,ii,j)-
     &     (frac1+fracm)*rmom(mxx,ii,j)))
      fmome(mx)=am(i)*(fracm*fracm*(rmom(mx,ii,j)
     &     -3.*frac1*rmom(mxx,ii,j))-3.*fe)
      fmome(mxx)=am(i)*(am(i)*fracm**3 *rmom(mxx,ii,j)
     &     -5.*(am(i)*fe+fmome(mx)))

      ! cross moments
      fmome(my)  = fracm*(rmom(my,ii,j)-frac1*rmom(mxy,ii,j))
      fmome(mxy) = am(i)*(fracm*fracm*rmom(mxy,ii,j)-3.*fmome(my))
      fmome(mz)  = fracm*(rmom(mz,ii,j)-frac1*rmom(mzx,ii,j))
      fmome(mzx) = am(i)*(fracm*fracm*rmom(mzx,ii,j)-3.*fmome(mz))
      fmome(myy) = fracm*rmom(myy,ii,j)
      fmome(mzz) = fracm*rmom(mzz,ii,j)
      fmome(myz) = fracm*rmom(myz,ii,j)

! flux limitations
      fe0 = fe
      fe_pass = fe
      fex_pass = fmome(mx)
      fexx_pass = fmome(mxx)
      if(am(i).gt.0.) then
        if(fe.lt.0.) then
          fe=0.
          fe_pass=0.
          fex_pass=0
          fexx_pass=0.
        elseif(fe.gt.rm(i,j)) then
          fe=rm(i,j)
          fe_pass = fe
          fex_pass=am(i)*(-3.*fe)
          fexx_pass=am(i)*(-5.*(am(i)*fe+fex_pass))
        endif
      else
        if(fe.gt.0.) then
          fe=0.
          fe0=0.
          fmome((/mx,mxx/))=0.
        elseif(fe.lt.-rm(1,j)) then
          fe=-rm(1,j)
          fe0 = fe
          fmome(mx)=am(i)*(-3.*fe)
          fmome(mxx)=am(i)*(-5.*(am(i)*fe+fmome(mx)))
        endif
      endif
      feim = fe
      feim0 = fe0
      fmomeim(:) = fmome(:)
      amw = am(im)
      fw = fe
      fw0 = fe_pass
      fmomw(:) = fmome(:)
      fmomw(mx) = fex_pass
      fmomw(mxx) = fexx_pass
      safv(i,j)  = safv(i,j) + fe

      do i=1,im-1
         if(am(i).lt.0.) then ! air mass flux is negative
            ii=i+1
            frac1=+1.
         else                 ! air mass flux is positive
            ii=i
            frac1=-1.
         endif
         fracm=am(i)/mass(ii,j)
         frac1=fracm+frac1
         fe=fracm*(rm(ii,j)-frac1*(rmom(mx,ii,j)-
     &        (frac1+fracm)*rmom(mxx,ii,j)))
         fmome(mx)=am(i)*(fracm*fracm*(rmom(mx,ii,j)
     &        -3.*frac1*rmom(mxx,ii,j))-3.*fe)
         fmome(mxx)=am(i)*(am(i)*fracm**3 *rmom(mxx,ii,j)
     &        -5.*(am(i)*fe+fmome(mx)))
      ! cross moments
         fmome(my)  = fracm*(rmom(my,ii,j)-frac1*rmom(mxy,ii,j))
         fmome(mxy) = am(i)*(fracm*fracm*rmom(mxy,ii,j)-3.*fmome(my))
         fmome(mz)  = fracm*(rmom(mz,ii,j)-frac1*rmom(mzx,ii,j))
         fmome(mzx) = am(i)*(fracm*fracm*rmom(mzx,ii,j)-3.*fmome(mz))
         fmome(myy) = fracm*rmom(myy,ii,j)
         fmome(mzz) = fracm*rmom(mzz,ii,j)
         fmome(myz) = fracm*rmom(myz,ii,j)

! flux limitations
         fe0 = fe
         fe_pass = fe
         fex_pass = fmome(mx)
         fexx_pass = fmome(mxx)
         if(am(i).gt.0.) then
           if(fe.lt.0.) then
             fe=0.
             fe_pass=0.
             fex_pass=0
             fexx_pass=0.
           elseif(fe.gt.rm(i,j)) then
             fe=rm(i,j)
             fe_pass = fe
             fex_pass=am(i)*(-3.*fe)
             fexx_pass=am(i)*(-5.*(am(i)*fe+fex_pass))
           endif
         else
           if(fe.gt.0.) then
             fe=0.
             fe0=0.
             fmome((/mx,mxx/))=0.
           elseif(fe.lt.-rm(i+1,j)) then
             fe=-rm(i+1,j)
             fe0 = fe
             fmome(mx)=am(i)*(-3.*fe)
             fmome(mxx)=am(i)*(-5.*(am(i)*fe+fmome(mx)))
           endif
         endif

         mold=mass(i,j)
         mnew=mold+amw-am(i)
         bymnew = 1./mnew
         dm2=amw+am(i)
         rm0=rm(i,j)+fw0-fe0
         rm(i,j)=rm(i,j)+fw-fe
      !
         rmom(mx,i,j)=(rmom(mx,i,j)*mold-3.*(-dm2*rm0
     &     +mold*(fw0+fe0))+(fmomw(mx)-fmome(mx)))*bymnew
         rmom(mxx,i,j) = (rmom(mxx,i,j)*mold*mold
     &     +2.5*rm0*(mold*mold-mnew*mnew-3.*dm2*dm2)
     &     +5.*(mold*(mold*(fw0-fe0)-fmomw(mx)
     &     -fmome(mx))+dm2*rmom(mx,i,j)*mnew)
     &     +(fmomw(mxx)-fmome(mxx))) * (bymnew*bymnew)
      ! cross moments
         rmom(my,i,j)=rmom(my,i,j)+fmomw(my)-fmome(my)
         rmom(mxy,i,j)=(rmom(mxy,i,j)*mold-3.*(-dm2*rmom(my,i,j) +
     &        mold*(fmomw(my)+fmome(my))) +
     &        (fmomw(mxy)-fmome(mxy)))*bymnew
         rmom(mz,i,j)=rmom(mz,i,j)+fmomw(mz)-fmome(mz)
         rmom(mzx,i,j)=(rmom(mzx,i,j)*mold-3.*(-dm2*rmom(mz,i,j) +
     &        mold*(fmomw(mz)+fmome(mz))) +
     &        (fmomw(mzx)-fmome(mzx)))*bymnew
      !
         rmom(myy,i,j)=rmom(myy,i,j)+fmomw(myy)-fmome(myy)
         rmom(mzz,i,j)=rmom(mzz,i,j)+fmomw(mzz)-fmome(mzz)
         rmom(myz,i,j)=rmom(myz,i,j)+fmomw(myz)-fmome(myz)

         mass(i,j) = mnew

! clean up roundoff errors
         if(rm(i,j).le.0d0) then
           rm(i,j)=0d0; rmom(:,i,j)=0d0
         endif

         amw = am(i)
         fw = fe
         fmomw(:) = fmome(:)
         fw0 = fe_pass
         fmomw(mx) = fex_pass
         fmomw(mxx) = fexx_pass
         safv(i,j)  = safv(i,j) + fe

      enddo ! i

      i = im
      fe = feim
      fe0 = feim0
      fmome(:) = fmomeim(:)
      mold=mass(i,j)
      mnew=mold+amw-am(i)
      bymnew = 1./mnew
      dm2=amw+am(i)
      rm0=rm(i,j)+fw0-fe0
      rm(i,j)=rm(i,j)+fw-fe
      !
      rmom(mx,i,j)=(rmom(mx,i,j)*mold-3.*(-dm2*rm0
     &     +mold*(fw0+fe0))+(fmomw(mx)-fmome(mx)))*bymnew
      rmom(mxx,i,j) = (rmom(mxx,i,j)*mold*mold
     &     +2.5*rm0*(mold*mold-mnew*mnew-3.*dm2*dm2)
     &     +5.*(mold*(mold*(fw0-fe0)-fmomw(mx)
     &     -fmome(mx))+dm2*rmom(mx,i,j)*mnew)
     &     +(fmomw(mxx)-fmome(mxx))) * (bymnew*bymnew)
      ! cross moments
      rmom(my,i,j)=rmom(my,i,j)+fmomw(my)-fmome(my)
      rmom(mxy,i,j)=(rmom(mxy,i,j)*mold-3.*(-dm2*rmom(my,i,j) +
     &     mold*(fmomw(my)+fmome(my))) +
     &     (fmomw(mxy)-fmome(mxy)))*bymnew
      rmom(mz,i,j)=rmom(mz,i,j)+fmomw(mz)-fmome(mz)
      rmom(mzx,i,j)=(rmom(mzx,i,j)*mold-3.*(-dm2*rmom(mz,i,j) +
     &     mold*(fmomw(mz)+fmome(mz))) +
     &     (fmomw(mzx)-fmome(mzx)))*bymnew
      !
      rmom(myy,i,j)=rmom(myy,i,j)+fmomw(myy)-fmome(myy)
      rmom(mzz,i,j)=rmom(mzz,i,j)+fmomw(mzz)-fmome(mzz)
      rmom(myz,i,j)=rmom(myz,i,j)+fmomw(myz)-fmome(myz)

      mass(i,j) = mnew

! clean up roundoff errors
      if(rm(i,j).le.0d0) then
        rm(i,j)=0d0; rmom(:,i,j)=0d0
      endif

      enddo ! ns
      enddo ! j

      return
c****
      end subroutine aadvqx

      subroutine aadvqy(rm,rmom,mass,mv  ,sbf,sbm,sfbm,sbfv)
!@sum  AADVQY advection driver for y-direction
!@auth Maxwell Kelley
      use DOMAIN_DECOMP, only : grid, get
      use DOMAIN_DECOMP, only : halo_update,halo_update_column
      use DOMAIN_DECOMP, only : NORTH, SOUTH, AM_I_ROOT
      use QUSDEF
      use QUSCOM, only : im,jm,byim
      implicit none
      REAL*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) :: 
     &                  rm,mass,mv
      REAL*8, dimension(NMOM,IM,grid%J_STRT_HALO:
     &                          grid%J_STOP_HALO) :: rmom
      REAL*8, dimension(GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: 
     &                                         sfbm,sbm,sbf
      REAL*8, dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: sbfv
      integer :: i,j,jj
      REAL*8 :: m_sp,m_np,rm_sp,rm_np,rzm_sp,rzm_np,rzzm_sp,rzzm_np
      real*8, dimension(im) :: mvs,fs
      real*8, dimension(nmom,im) :: fmoms
      real*8, dimension(nmom) :: fmomn
      real*8 :: frac1,fracm,fn,mold,mnew,bymnew,dm2
      real*8 :: fn0,fn_pass,rm0,fny_pass,fnyy_pass
      real*8, dimension(im) :: fs0

c****Get relevant local distributed parameters
      INTEGER J_0,J_1,J_0H,J_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0,
     &               J_STOP = J_1,
     &               J_STRT_HALO = J_0H,
     &               J_STOP_HALO = J_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

c      if(rehalo_mom) then
c        CALL HALO_UPDATE_COLUMN(grid, rmom, FROM=NORTH+SOUTH)
c      endif

c**** scale polar boxes to their full extent
! set horizontal moments to zero at pole
      if (HAVE_SOUTH_POLE) then
        mass(:,1)=mass(:,1)*im
        m_sp = mass(1,1 )
        rm(:,1)=rm(:,1)*im
        rm_sp = rm(1,1 )
        do i=1,im
           rmom(:,i,1 )=rmom(:,i,1 )*im
        enddo
        rzm_sp  = rmom(mz ,1,1 )
        rzzm_sp = rmom(mzz,1,1 )
        rmom(ihmoms,:,1) = 0.
c        rmom((/mx,mxx,myy,mxy,mzx/),:,1) = 0.
      end if                       !SOUTH POLE

      if (HAVE_NORTH_POLE) then
        mass(:,jm)=mass(:,jm)*im
        m_np = mass(1,jm)
        rm(:,jm)=rm(:,jm)*im
        rm_np = rm(1,jm)
        do i=1,im
           rmom(:,i,jm)=rmom(:,i,jm)*im
        enddo
        rzm_np  = rmom(mz ,1,jm)
        rzzm_np = rmom(mzz,1,jm)
        rmom(ihmoms,:,jm) = 0.
c        rmom((/mx,mxx,myy,mxy,mzx/),:,jm) = 0.
      end if                       !NORTH POLE

      if(have_north_pole) then
        mv(:,jm) = 0.
      endif
      if(have_south_pole) then
        mvs(:) = 0.
        fs(:) = 0.
        fmoms(:,:) = 0.
      else
        j=j_0-1
        do i=1,im
          if(mv(i,j).lt.0.) then ! air mass flux is negative
            jj=j+1
            frac1=+1.
          else                  ! air mass flux is positive
            jj=j
            frac1=-1.
          endif
          fracm=mv(i,j)/mass(i,jj)
          frac1=fracm+frac1
          fn=fracm*(rm(i,jj)-frac1*(rmom(my,i,jj)-
     &         (frac1+fracm)*rmom(myy,i,jj)))
          fmomn(my)=mv(i,j)*(fracm*fracm*(rmom(my,i,jj)
     &         -3.*frac1*rmom(myy,i,jj))-3.*fn)
          fmomn(myy)=mv(i,j)*(mv(i,j)*fracm**3 *rmom(myy,i,jj)
     &         -5.*(mv(i,j)*fn+fmomn(my)))
          fmomn(mz)  = fracm*(rmom(mz,i,jj)-frac1*rmom(myz,i,jj))
          fmomn(myz) = mv(i,j)*
     &         (fracm*fracm*rmom(myz,i,jj)-3.*fmomn(mz))
          fmomn(mx)  = fracm*(rmom(mx,i,jj)-frac1*rmom(mxy,i,jj))
          fmomn(mxy) = mv(i,j)*
     &         (fracm*fracm*rmom(mxy,i,jj)-3.*fmomn(mx))
          fmomn(mzz) = fracm*rmom(mzz,i,jj)
          fmomn(mxx) = fracm*rmom(mxx,i,jj)
          fmomn(mzx) = fracm*rmom(mzx,i,jj)

! flux limitations
          fn0 = fn
          fn_pass = fn
          fny_pass = fmomn(my)
          fnyy_pass = fmomn(myy)
          if(mv(i,j).gt.0.) then
            if(fn.lt.0.) then
              fn=0.
              fn_pass=0.
              fny_pass=0
              fnyy_pass=0.
            elseif(fn.gt.rm(i,j)) then
              fn=rm(i,j)
              fn_pass = fn
              fny_pass=mv(i,j)*(-3.*fn)
              fnyy_pass=mv(i,j)*(-5.*(mv(i,j)*fn+fny_pass))
            endif
          else
            if(fn.gt.0.) then
              fn=0.
              fn0=0.
              fmomn((/my,myy/))=0.
            elseif(fn.lt.-rm(i,j+1)) then
              fn=-rm(i,j+1)
              fn0 = fn
              fmomn(my)=mv(i,j)*(-3.*fn)
              fmomn(myy)=mv(i,j)*(-5.*(mv(i,j)*fn+fmomn(my)))
            endif
          endif

          mvs(i) = mv(i,j)
          fs(i) = fn
          fmoms(:,i) = fmomn(:)

          fs0(i) = fn_pass
          fmoms(my,i) = fny_pass
          fmoms(myy,i) = fnyy_pass

        enddo ! i
      endif
      do j=j_0,j_1
      do i=1,im
         if(mv(i,j).lt.0.) then ! air mass flux is negative
            jj=j+1
            frac1=+1.
          else                   ! air mass flux is positive
            jj=j
            frac1=-1.
         endif
         fracm=mv(i,j)/mass(i,jj)
         frac1=fracm+frac1
         fn=fracm*(rm(i,jj)-frac1*(rmom(my,i,jj)-
     &        (frac1+fracm)*rmom(myy,i,jj)))
         fmomn(my)=mv(i,j)*(fracm*fracm*(rmom(my,i,jj)
     &        -3.*frac1*rmom(myy,i,jj))-3.*fn)
         fmomn(myy)=mv(i,j)*(mv(i,j)*fracm**3 *rmom(myy,i,jj)
     &        -5.*(mv(i,j)*fn+fmomn(my)))
         fmomn(mz)  = fracm*(rmom(mz,i,jj)-frac1*rmom(myz,i,jj))
         fmomn(myz) = mv(i,j)*
     &        (fracm*fracm*rmom(myz,i,jj)-3.*fmomn(mz))
         fmomn(mx)  = fracm*(rmom(mx,i,jj)-frac1*rmom(mxy,i,jj))
         fmomn(mxy) = mv(i,j)*
     &        (fracm*fracm*rmom(mxy,i,jj)-3.*fmomn(mx))
         fmomn(mzz) = fracm*rmom(mzz,i,jj)
         fmomn(mxx) = fracm*rmom(mxx,i,jj)
         fmomn(mzx) = fracm*rmom(mzx,i,jj)

! flux limitations
         fn0 = fn
         fn_pass = fn
         fny_pass = fmomn(my)
         fnyy_pass = fmomn(myy)
         if(mv(i,j).gt.0.) then
           if(fn.lt.0.) then
             fn=0.
             fn_pass=0.
             fny_pass=0
             fnyy_pass=0.
           elseif(fn.gt.rm(i,j)) then
             fn=rm(i,j)
             fn_pass = fn
             fny_pass=mv(i,j)*(-3.*fn)
             fnyy_pass=mv(i,j)*(-5.*(mv(i,j)*fn+fny_pass))
           endif
         elseif(mv(i,j).lt.0.) then
           if(fn.gt.0.) then
             fn=0.
             fn0=0.
             fmomn((/my,myy/))=0.
           elseif(fn.lt.-rm(i,j+1)) then
             fn=-rm(i,j+1)
             fn0 = fn
             fmomn(my)=mv(i,j)*(-3.*fn)
             fmomn(myy)=mv(i,j)*(-5.*(mv(i,j)*fn+fmomn(my)))
           endif
         endif

         mold=mass(i,j)
         mnew=mold+mvs(i)-mv(i,j)
         bymnew = 1./mnew
         dm2=mvs(i)+mv(i,j)
         rm0=rm(i,j)+fs0(i)-fn0
         rm(i,j)=rm(i,j)+fs(i)-fn
      !
         rmom(my,i,j)=(rmom(my,i,j)*mold-3.*(-dm2*rm0
     &     +mold*(fs0(i)+fn0))+(fmoms(my,i)-fmomn(my)))*bymnew
         rmom(myy,i,j) = (rmom(myy,i,j)*mold*mold
     &     +2.5*rm0*(mold*mold-mnew*mnew-3.*dm2*dm2)
     &     +5.*(mold*(mold*(fs0(i)-fn0)-fmoms(my,i)
     &     -fmomn(my))+dm2*rmom(my,i,j)*mnew)
     &     +(fmoms(myy,i)-fmomn(myy))) * (bymnew*bymnew)
      ! cross moments
         rmom(mz,i,j)=rmom(mz,i,j)+fmoms(mz,i)-fmomn(mz)
         rmom(myz,i,j)=(rmom(myz,i,j)*mold-3.*(-dm2*rmom(mz,i,j) +
     &        mold*(fmoms(mz,i)+fmomn(mz))) +
     &        (fmoms(myz,i)-fmomn(myz)))*bymnew
         rmom(mx,i,j)=rmom(mx,i,j)+fmoms(mx,i)-fmomn(mx)
         rmom(mxy,i,j)=(rmom(mxy,i,j)*mold-3.*(-dm2*rmom(mx,i,j) +
     &        mold*(fmoms(mx,i)+fmomn(mx))) +
     &        (fmoms(mxy,i)-fmomn(mxy)))*bymnew
      !
         rmom(mzz,i,j)=rmom(mzz,i,j)+fmoms(mzz,i)-fmomn(mzz)
         rmom(mxx,i,j)=rmom(mxx,i,j)+fmoms(mxx,i)-fmomn(mxx)
         rmom(mzx,i,j)=rmom(mzx,i,j)+fmoms(mzx,i)-fmomn(mzx)

         mass(i,j) = mnew

! clean up roundoff errors
         if(rm(i,j).le.0d0) then
           rm(i,j)=0d0; rmom(:,i,j)=0d0
         endif

         mvs(i) = mv(i,j)
         fs(i) = fn
         fmoms(:,i) = fmomn(:)

         fs0(i) = fn_pass
         fmoms(my,i) = fny_pass
         fmoms(myy,i) = fnyy_pass

         sbfv(i,j) = sbfv(i,j) + fn
         sbf(j) = sbf(j) + fn
         sbm(j) = sbm(j) + mv(i,j)
         if(mv(i,j).ne.0.) sfbm(j) = sfbm(j) + fn/mv(i,j)

      enddo ! i
      enddo ! j

c**** average and unscale polar boxes
! horizontal moments are zero at pole
      if (HAVE_SOUTH_POLE) then
        mass(:,1 ) = (m_sp + sum(mass(:,1 )-m_sp))*byim
        rm(:,1 ) = (rm_sp + sum(rm(:,1 )-rm_sp))*byim
        rmom(mz ,:,1 ) = (rzm_sp + 
     *       sum(rmom(mz ,:,1 )-rzm_sp ))*byim
        rmom(mzz,:,1 ) = (rzzm_sp+ 
     *       sum(rmom(mzz,:,1 )-rzzm_sp))*byim
        rmom(ihmoms,:,1) = 0
      end if   !SOUTH POLE

      if (HAVE_NORTH_POLE) then
        mass(:,jm) = (m_np + sum(mass(:,jm)-m_np))*byim
        rm(:,jm) = (rm_np + sum(rm(:,jm)-rm_np))*byim
        rmom(mz ,:,jm) = (rzm_np + 
     &       sum(rmom(mz ,:,jm)-rzm_np ))*byim
        rmom(mzz,:,jm) = (rzzm_np+ 
     &       sum(rmom(mzz,:,jm)-rzzm_np))*byim
        rmom(ihmoms,:,jm) = 0.
      end if  !NORTH POLE

      return
      end subroutine aadvqy

      subroutine aadvqz(rm,rmom,mass,mw,mwdn,fdn,fmomdn,fdn0
     &     ,scf,scm,sfcm)
!@sum  AADVQZ advection driver for z-direction
!@auth Maxwell Kelley
      use DOMAIN_DECOMP, only : grid, GET
      use QUSDEF
      use QUSCOM, only : im
      implicit none
      REAL*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,2) 
     &        :: rm,mass
      REAL*8, dimension(NMOM,IM,grid%j_strt_halo:grid%j_stop_halo,2) 
     &        :: rmom
      REAL*8, dimension(NMOM,IM,grid%j_strt_halo:grid%j_stop_halo) 
     &        :: fmomdn
      REAL*8, dimension(IM,grid%j_strt_halo:grid%j_stop_halo) ::
     &     mw,mwdn,fdn,fdn0
      REAL*8, dimension(GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     sfcm,scm,scf
      real*8, dimension(nmom) :: fmomup
      integer :: i,j,l,ll
      real*8 :: frac1,fracm,fup,mold,mnew,bymnew,dm2
      real*8 :: fup0,fup_pass,rm0,fupz_pass,fupzz_pass

c**** Get useful local parameters for domain decomposition
      integer :: J_0, J_1
      CALL GET( grid, J_STRT=J_0 , J_STOP=J_1 )

      do j=j_0,j_1
      do i=1,im
         if(mw(i,j).lt.0.) then ! air mass flux is negative
            ll=2
            frac1=+1.
          else                   ! air mass flux is positive
            ll=1
            frac1=-1.
         endif
         fracm=mw(i,j)/mass(i,j,ll)
         frac1=fracm+frac1
         fup=fracm*(rm(i,j,ll)-frac1*(rmom(mz,i,j,ll)-
     &        (frac1+fracm)*rmom(mzz,i,j,ll)))
         fmomup(mz)=mw(i,j)*(fracm*fracm*(rmom(mz,i,j,ll)
     &        -3.*frac1*rmom(mzz,i,j,ll))-3.*fup)
         fmomup(mzz)=mw(i,j)*(mw(i,j)*fracm**3 *rmom(mzz,i,j,ll)
     &        -5.*(mw(i,j)*fup+fmomup(mz)))
         fmomup(my)  = fracm*(rmom(my,i,j,ll)-frac1*rmom(myz,i,j,ll))
         fmomup(myz) = mw(i,j)*
     &        (fracm*fracm*rmom(myz,i,j,ll)-3.*fmomup(my))
         fmomup(mx)  = fracm*(rmom(mx,i,j,ll)-frac1*rmom(mzx,i,j,ll))
         fmomup(mzx) = mw(i,j)*
     &        (fracm*fracm*rmom(mzx,i,j,ll)-3.*fmomup(mx))
         fmomup(myy) = fracm*rmom(myy,i,j,ll)
         fmomup(mxx) = fracm*rmom(mxx,i,j,ll)
         fmomup(mxy) = fracm*rmom(mxy,i,j,ll)

! flux limitations
         fup0 = fup
         fup_pass = fup
         fupz_pass = fmomup(mz)
         fupzz_pass = fmomup(mzz)
         if(mw(i,j).gt.0.) then
           if(fup.lt.0.) then
             fup=0.
             fup_pass=0.
             fupz_pass=0
             fupzz_pass=0.
           elseif(fup.gt.rm(i,j,1)) then
             fup=rm(i,j,1)
             fup_pass = fup
             fupz_pass=mw(i,j)*(-3.*fup)
             fupzz_pass=mw(i,j)*(-5.*(mw(i,j)*fup+fupz_pass))
           endif
         elseif(mw(i,j).lt.0.) then
           if(fup.gt.0.) then
             fup=0.
             fup0=0.
             fmomup((/mz,mzz/))=0.
           elseif(fup.lt.-rm(i,j,2)) then
             fup=-rm(i,j,2)
             fup0 = fup
             fmomup(mz)=mw(i,j)*(-3.*fup)
             fmomup(mzz)=mw(i,j)*(-5.*(mw(i,j)*fup+fmomup(mz)))
           endif
         endif

         mold=mass(i,j,1)
         mnew=mold+mwdn(i,j)-mw(i,j)
         bymnew = 1./mnew
         dm2=mwdn(i,j)+mw(i,j)
         rm0=rm(i,j,1)+fdn0(i,j)-fup0
         rm(i,j,1)=rm(i,j,1)+fdn(i,j)-fup
      !
         rmom(mz,i,j,1)=(rmom(mz,i,j,1)*mold-3.*(-dm2*rm0
     &     +mold*(fdn0(i,j)+fup0))+(fmomdn(mz,i,j)-fmomup(mz)))*bymnew
         rmom(mzz,i,j,1) = (rmom(mzz,i,j,1)*mold*mold
     &     +2.5*rm0*(mold*mold-mnew*mnew-3.*dm2*dm2)
     &     +5.*(mold*(mold*(fdn0(i,j)-fup0)-fmomdn(mz,i,j)
     &     -fmomup(mz))+dm2*rmom(mz,i,j,1)*mnew)
     &     +(fmomdn(mzz,i,j)-fmomup(mzz))) * (bymnew*bymnew)
      ! cross moments
         rmom(my,i,j,1)=rmom(my,i,j,1)+fmomdn(my,i,j)-fmomup(my)
         rmom(myz,i,j,1)=(rmom(myz,i,j,1)*mold-3.*(-dm2*rmom(my,i,j,1) +
     &        mold*(fmomdn(my,i,j)+fmomup(my))) +
     &        (fmomdn(myz,i,j)-fmomup(myz)))*bymnew
         rmom(mx,i,j,1)=rmom(mx,i,j,1)+fmomdn(mx,i,j)-fmomup(mx)
         rmom(mzx,i,j,1)=(rmom(mzx,i,j,1)*mold-3.*(-dm2*rmom(mx,i,j,1) +
     &        mold*(fmomdn(mx,i,j)+fmomup(mx))) +
     &        (fmomdn(mzx,i,j)-fmomup(mzx)))*bymnew
      !
         rmom(myy,i,j,1)=rmom(myy,i,j,1)+fmomdn(myy,i,j)-fmomup(myy)
         rmom(mxx,i,j,1)=rmom(mxx,i,j,1)+fmomdn(mxx,i,j)-fmomup(mxx)
         rmom(mxy,i,j,1)=rmom(mxy,i,j,1)+fmomdn(mxy,i,j)-fmomup(mxy)

         mass(i,j,1) = mnew

! clean up roundoff errors
         if(rm(i,j,1).le.0d0) then
           rm(i,j,1)=0d0; rmom(:,i,j,1)=0d0
         endif

         mwdn(i,j) = mw(i,j)
         fdn(i,j) = fup
         fmomdn(:,i,j) = fmomup(:)

         fdn0(i,j) = fup_pass
         fmomdn(mz,i,j) = fupz_pass
         fmomdn(mzz,i,j) = fupzz_pass

         scf(j) = scf(j) + fup
         scm(j) = scm(j) + mw(i,j)
         if(mw(i,j).ne.0.) sfcm(j) = sfcm(j) + fup/mw(i,j)

      enddo ! i
      enddo ! j

      return
c****
      end subroutine aadvqz

      subroutine checkflux(aml,amr,m,rm,rxm,rxxm)
      implicit none
      real*8 :: aml,amr,m,rm,rxm,rxxm
      real*8 :: a,fl,fr
c flux out the right side
      A = AMR / M
      FR = A*(RM + (1.-A)*(RXM + (1.-2.*A)*RXXM))
c flux out the left side
      A = AML / M
      FL = A*(RM - (1.+A)*(RXM - (1.+2.*A)*RXXM))
c
      if(rm+fl-fr.le.0.) then
        rxm = 0.
        rxxm = 0.
      endif
      return
      end subroutine checkflux

      subroutine aadvqx2(rm,rmom,mass,mu,jmin,jmax,nstep,safv)
!@sum  AADVQX2 version of AADVQX without flux limiter
!@auth Maxwell Kelley
      use DOMAIN_DECOMP, only : grid, GET
      use QUSDEF
      use QUSCOM, only : im,jm
      implicit none
      integer :: jmin,jmax
      REAL*8, dimension(im,grid%J_STRT_HALO:grid%J_STOP_HALO) :: 
     &                  rm,mass,mu
      REAL*8, dimension(NMOM,IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &                  rmom
      integer, dimension(grid%J_STRT_HALO:grid%J_STOP_HALO) :: nstep
      REAL*8, dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: safv

      integer :: i,ii,j,ns
      real*8 :: frac1,fracm,fw,fe,feim,amw,dm2,mold,mnew,bymnew
      real*8, dimension(nmom) :: fmomw,fmome,fmomeim
      real*8, dimension(im) :: am

c**** Get useful local parameters for domain decomposition
      integer :: J_0, J_1, J_0S, J_1S
      CALL GET(grid, J_STRT = J_0 , J_STOP=J_1,
     &             J_STRT_SKP=J_0S,J_STOP_SKP=J_1S )

      do j=jmin,jmax
      am(:) = mu(:,j)/nstep(j)
c****
c**** loop over timesteps
c****
      do ns=1,nstep(j)

      i=im
      if(am(i).lt.0.) then  ! air mass flux is negative
        ii=1
        frac1=+1.
      else                      ! air mass flux is positive
        ii=i
        frac1=-1.
      endif
      fracm=am(i)/mass(ii,j)
      frac1=fracm+frac1
      fw=fracm*(rm(ii,j)-frac1*(rmom(mx,ii,j)-
     &     (frac1+fracm)*rmom(mxx,ii,j)))
      fmomw(mx)=am(i)*(fracm*fracm*(rmom(mx,ii,j)
     &     -3.*frac1*rmom(mxx,ii,j))-3.*fw)
      fmomw(mxx)=am(i)*(am(i)*fracm**3 *rmom(mxx,ii,j)
     &     -5.*(am(i)*fw+fmomw(mx)))
      ! cross moments
      fmomw(my)  = fracm*(rmom(my,ii,j)-frac1*rmom(mxy,ii,j))
      fmomw(mxy) = am(i)*(fracm*fracm*rmom(mxy,ii,j)-3.*fmomw(my))
      fmomw(mz)  = fracm*(rmom(mz,ii,j)-frac1*rmom(mzx,ii,j))
      fmomw(mzx) = am(i)*(fracm*fracm*rmom(mzx,ii,j)-3.*fmomw(mz))
      fmomw(myy) = fracm*rmom(myy,ii,j)
      fmomw(mzz) = fracm*rmom(mzz,ii,j)
      fmomw(myz) = fracm*rmom(myz,ii,j)
      feim = fw
      fmomeim(:) = fmomw(:)
      amw = am(im)
      safv(i,j)  = safv(i,j) + fw

      do i=1,im-1
         if(am(i).lt.0.) then ! air mass flux is negative
            ii=i+1
            frac1=+1.
         else                 ! air mass flux is positive
            ii=i
            frac1=-1.
         endif
         fracm=am(i)/mass(ii,j)
         frac1=fracm+frac1
         fe=fracm*(rm(ii,j)-frac1*(rmom(mx,ii,j)-
     &        (frac1+fracm)*rmom(mxx,ii,j)))
         fmome(mx)=am(i)*(fracm*fracm*(rmom(mx,ii,j)
     &        -3.*frac1*rmom(mxx,ii,j))-3.*fe)
         fmome(mxx)=am(i)*(am(i)*fracm**3 *rmom(mxx,ii,j)
     &        -5.*(am(i)*fe+fmome(mx)))
      ! cross moments
         fmome(my)  = fracm*(rmom(my,ii,j)-frac1*rmom(mxy,ii,j))
         fmome(mxy) = am(i)*(fracm*fracm*rmom(mxy,ii,j)-3.*fmome(my))
         fmome(mz)  = fracm*(rmom(mz,ii,j)-frac1*rmom(mzx,ii,j))
         fmome(mzx) = am(i)*(fracm*fracm*rmom(mzx,ii,j)-3.*fmome(mz))
         fmome(myy) = fracm*rmom(myy,ii,j)
         fmome(mzz) = fracm*rmom(mzz,ii,j)
         fmome(myz) = fracm*rmom(myz,ii,j)

         mold=mass(i,j)
         mnew=mold+amw-am(i)
         bymnew = 1./mnew
         dm2=amw+am(i)
         rm(i,j)=rm(i,j)+fw-fe
      !
         rmom(mx,i,j)=(rmom(mx,i,j)*mold-3.*(-dm2*rm(i,j)
     &     +mold*(fw+fe))+(fmomw(mx)-fmome(mx)))*bymnew
         rmom(mxx,i,j) = (rmom(mxx,i,j)*mold*mold
     &     +2.5*rm(i,j)*(mold*mold-mnew*mnew-3.*dm2*dm2)
     &     +5.*(mold*(mold*(fw-fe)-fmomw(mx)
     &     -fmome(mx))+dm2*rmom(mx,i,j)*mnew)
     &     +(fmomw(mxx)-fmome(mxx))) * (bymnew*bymnew)
      ! cross moments
         rmom(my,i,j)=rmom(my,i,j)+fmomw(my)-fmome(my)
         rmom(mxy,i,j)=(rmom(mxy,i,j)*mold-3.*(-dm2*rmom(my,i,j) +
     &        mold*(fmomw(my)+fmome(my))) +
     &        (fmomw(mxy)-fmome(mxy)))*bymnew
         rmom(mz,i,j)=rmom(mz,i,j)+fmomw(mz)-fmome(mz)
         rmom(mzx,i,j)=(rmom(mzx,i,j)*mold-3.*(-dm2*rmom(mz,i,j) +
     &        mold*(fmomw(mz)+fmome(mz))) +
     &        (fmomw(mzx)-fmome(mzx)))*bymnew
      !
         rmom(myy,i,j)=rmom(myy,i,j)+fmomw(myy)-fmome(myy)
         rmom(mzz,i,j)=rmom(mzz,i,j)+fmomw(mzz)-fmome(mzz)
         rmom(myz,i,j)=rmom(myz,i,j)+fmomw(myz)-fmome(myz)

         mass(i,j) = mnew

         amw = am(i)
         fw = fe
         fmomw(:) = fmome(:)

         safv(i,j)  = safv(i,j) + fe

      enddo ! i

      i = im
      fe = feim
      fmome(:) = fmomeim(:)
      mold=mass(i,j)
      mnew=mold+amw-am(i)
      bymnew = 1./mnew
      dm2=amw+am(i)
      rm(i,j)=rm(i,j)+fw-fe
      !
      rmom(mx,i,j)=(rmom(mx,i,j)*mold-3.*(-dm2*rm(i,j)
     &     +mold*(fw+fe))+(fmomw(mx)-fmome(mx)))*bymnew
      rmom(mxx,i,j) = (rmom(mxx,i,j)*mold*mold
     &     +2.5*rm(i,j)*(mold*mold-mnew*mnew-3.*dm2*dm2)
     &     +5.*(mold*(mold*(fw-fe)-fmomw(mx)
     &     -fmome(mx))+dm2*rmom(mx,i,j)*mnew)
     &     +(fmomw(mxx)-fmome(mxx))) * (bymnew*bymnew)
      ! cross moments
      rmom(my,i,j)=rmom(my,i,j)+fmomw(my)-fmome(my)
      rmom(mxy,i,j)=(rmom(mxy,i,j)*mold-3.*(-dm2*rmom(my,i,j) +
     &     mold*(fmomw(my)+fmome(my))) +
     &     (fmomw(mxy)-fmome(mxy)))*bymnew
      rmom(mz,i,j)=rmom(mz,i,j)+fmomw(mz)-fmome(mz)
      rmom(mzx,i,j)=(rmom(mzx,i,j)*mold-3.*(-dm2*rmom(mz,i,j) +
     &     mold*(fmomw(mz)+fmome(mz))) +
     &     (fmomw(mzx)-fmome(mzx)))*bymnew
      !
      rmom(myy,i,j)=rmom(myy,i,j)+fmomw(myy)-fmome(myy)
      rmom(mzz,i,j)=rmom(mzz,i,j)+fmomw(mzz)-fmome(mzz)
      rmom(myz,i,j)=rmom(myz,i,j)+fmomw(myz)-fmome(myz)

      mass(i,j) = mnew

      enddo ! ns
      enddo ! j

      return
c****
      end subroutine aadvqx2

      subroutine aadvqy2(rm,rmom,mass,mv  ,sbf,sbm,sfbm,sbfv)
!@sum  AADVQY2 version of AADVQY without flux limiter
!@auth Maxwell Kelley
      use DOMAIN_DECOMP, only : grid, get
      use DOMAIN_DECOMP, only : halo_update,halo_update_column
      use DOMAIN_DECOMP, only : NORTH, SOUTH, AM_I_ROOT
      use QUSDEF
      use QUSCOM, only : im,jm,byim
      implicit none
      REAL*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) :: 
     &                  rm,mass,mv
      REAL*8, dimension(NMOM,IM,grid%J_STRT_HALO:
     &                          grid%J_STOP_HALO) :: rmom
      REAL*8, dimension(GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: 
     &                                         sfbm,sbm,sbf
      REAL*8, dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: sbfv
      integer :: i,j,jj
      REAL*8 :: m_sp,m_np,rm_sp,rm_np,rzm_sp,rzm_np,rzzm_sp,rzzm_np
      real*8, dimension(im) :: mvs,fs
      real*8, dimension(nmom,im) :: fmoms
      real*8, dimension(nmom) :: fmomn
      real*8 :: frac1,fracm,fn,mold,mnew,bymnew,dm2

c****Get relevant local distributed parameters
      INTEGER J_0,J_1,J_0H,J_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0,
     &               J_STOP = J_1,
     &               J_STRT_HALO = J_0H,
     &               J_STOP_HALO = J_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

c**** scale polar boxes to their full extent
! set horizontal moments to zero at pole
      if (HAVE_SOUTH_POLE) then
        mass(:,1)=mass(:,1)*im
        m_sp = mass(1,1 )
        rm(:,1)=rm(:,1)*im
        rm_sp = rm(1,1 )
        do i=1,im
           rmom(:,i,1 )=rmom(:,i,1 )*im
        enddo
        rzm_sp  = rmom(mz ,1,1 )
        rzzm_sp = rmom(mzz,1,1 )
c        rmom(ihmoms,:,1) = 0.
        rmom((/mx,mxx,myy,mxy,mzx/),:,1) = 0.
      end if                       !SOUTH POLE

      if (HAVE_NORTH_POLE) then
        mass(:,jm)=mass(:,jm)*im
        m_np = mass(1,jm)
        rm(:,jm)=rm(:,jm)*im
        rm_np = rm(1,jm)
        do i=1,im
           rmom(:,i,jm)=rmom(:,i,jm)*im
        enddo
        rzm_np  = rmom(mz ,1,jm)
        rzzm_np = rmom(mzz,1,jm)
c        rmom(ihmoms,:,jm) = 0.
        rmom((/mx,mxx,myy,mxy,mzx/),:,jm) = 0.
      end if                       !NORTH POLE

      if(have_north_pole) then
        mv(:,jm) = 0.
      endif
      if(have_south_pole) then
        mvs(:) = 0.
        fs(:) = 0.
        fmoms(:,:) = 0.
      else
        j=j_0-1
        do i=1,im
          if(mv(i,j).lt.0.) then ! air mass flux is negative
            jj=j+1
            frac1=+1.
          else                  ! air mass flux is positive
            jj=j
            frac1=-1.
          endif
          fracm=mv(i,j)/mass(i,jj)
          frac1=fracm+frac1
          fn=fracm*(rm(i,jj)-frac1*(rmom(my,i,jj)-
     &         (frac1+fracm)*rmom(myy,i,jj)))
          fmomn(my)=mv(i,j)*(fracm*fracm*(rmom(my,i,jj)
     &         -3.*frac1*rmom(myy,i,jj))-3.*fn)
          fmomn(myy)=mv(i,j)*(mv(i,j)*fracm**3 *rmom(myy,i,jj)
     &         -5.*(mv(i,j)*fn+fmomn(my)))
          fmomn(mz)  = fracm*(rmom(mz,i,jj)-frac1*rmom(myz,i,jj))
          fmomn(myz) = mv(i,j)*
     &         (fracm*fracm*rmom(myz,i,jj)-3.*fmomn(mz))
          fmomn(mx)  = fracm*(rmom(mx,i,jj)-frac1*rmom(mxy,i,jj))
          fmomn(mxy) = mv(i,j)*
     &         (fracm*fracm*rmom(mxy,i,jj)-3.*fmomn(mx))
          fmomn(mzz) = fracm*rmom(mzz,i,jj)
          fmomn(mxx) = fracm*rmom(mxx,i,jj)
          fmomn(mzx) = fracm*rmom(mzx,i,jj)
          mvs(i) = mv(i,j)
          fs(i) = fn
          fmoms(:,i) = fmomn(:)
        enddo ! i
      endif
      do j=j_0,j_1
      do i=1,im
         if(mv(i,j).lt.0.) then ! air mass flux is negative
            jj=j+1
            frac1=+1.
          else                   ! air mass flux is positive
            jj=j
            frac1=-1.
         endif
         fracm=mv(i,j)/mass(i,jj)
         frac1=fracm+frac1
         fn=fracm*(rm(i,jj)-frac1*(rmom(my,i,jj)-
     &        (frac1+fracm)*rmom(myy,i,jj)))
         fmomn(my)=mv(i,j)*(fracm*fracm*(rmom(my,i,jj)
     &        -3.*frac1*rmom(myy,i,jj))-3.*fn)
         fmomn(myy)=mv(i,j)*(mv(i,j)*fracm**3 *rmom(myy,i,jj)
     &        -5.*(mv(i,j)*fn+fmomn(my)))
         fmomn(mz)  = fracm*(rmom(mz,i,jj)-frac1*rmom(myz,i,jj))
         fmomn(myz) = mv(i,j)*
     &        (fracm*fracm*rmom(myz,i,jj)-3.*fmomn(mz))
         fmomn(mx)  = fracm*(rmom(mx,i,jj)-frac1*rmom(mxy,i,jj))
         fmomn(mxy) = mv(i,j)*
     &        (fracm*fracm*rmom(mxy,i,jj)-3.*fmomn(mx))
         fmomn(mzz) = fracm*rmom(mzz,i,jj)
         fmomn(mxx) = fracm*rmom(mxx,i,jj)
         fmomn(mzx) = fracm*rmom(mzx,i,jj)

         mold=mass(i,j)
         mnew=mold+mvs(i)-mv(i,j)
         bymnew = 1./mnew
         dm2=mvs(i)+mv(i,j)
         rm(i,j)=rm(i,j)+fs(i)-fn
      !
         rmom(my,i,j)=(rmom(my,i,j)*mold-3.*(-dm2*rm(i,j)
     &     +mold*(fs(i)+fn))+(fmoms(my,i)-fmomn(my)))*bymnew
         rmom(myy,i,j) = (rmom(myy,i,j)*mold*mold
     &     +2.5*rm(i,j)*(mold*mold-mnew*mnew-3.*dm2*dm2)
     &     +5.*(mold*(mold*(fs(i)-fn)-fmoms(my,i)
     &     -fmomn(my))+dm2*rmom(my,i,j)*mnew)
     &     +(fmoms(myy,i)-fmomn(myy))) * (bymnew*bymnew)
      ! cross moments
         rmom(mz,i,j)=rmom(mz,i,j)+fmoms(mz,i)-fmomn(mz)
         rmom(myz,i,j)=(rmom(myz,i,j)*mold-3.*(-dm2*rmom(mz,i,j) +
     &        mold*(fmoms(mz,i)+fmomn(mz))) +
     &        (fmoms(myz,i)-fmomn(myz)))*bymnew
         rmom(mx,i,j)=rmom(mx,i,j)+fmoms(mx,i)-fmomn(mx)
         rmom(mxy,i,j)=(rmom(mxy,i,j)*mold-3.*(-dm2*rmom(mx,i,j) +
     &        mold*(fmoms(mx,i)+fmomn(mx))) +
     &        (fmoms(mxy,i)-fmomn(mxy)))*bymnew
      !
         rmom(mzz,i,j)=rmom(mzz,i,j)+fmoms(mzz,i)-fmomn(mzz)
         rmom(mxx,i,j)=rmom(mxx,i,j)+fmoms(mxx,i)-fmomn(mxx)
         rmom(mzx,i,j)=rmom(mzx,i,j)+fmoms(mzx,i)-fmomn(mzx)

         mass(i,j) = mnew

         mvs(i) = mv(i,j)
         fs(i) = fn
         fmoms(:,i) = fmomn(:)

         sbfv(i,j) = sbfv(i,j) + fn
         sbf(j) = sbf(j) + fn
         sbm(j) = sbm(j) + mv(i,j)
         if(mv(i,j).ne.0.) sfbm(j) = sfbm(j) + fn/mv(i,j)

      enddo ! i
      enddo ! j

c**** average and unscale polar boxes
! horizontal moments are zero at pole
      if (HAVE_SOUTH_POLE) then
        mass(:,1 ) = (m_sp + sum(mass(:,1 )-m_sp))*byim
        rm(:,1 ) = (rm_sp + sum(rm(:,1 )-rm_sp))*byim
        rmom(mz ,:,1 ) = (rzm_sp + 
     *       sum(rmom(mz ,:,1 )-rzm_sp ))*byim
        rmom(mzz,:,1 ) = (rzzm_sp+ 
     *       sum(rmom(mzz,:,1 )-rzzm_sp))*byim
        rmom(ihmoms,:,1) = 0
      end if   !SOUTH POLE

      if (HAVE_NORTH_POLE) then
        mass(:,jm) = (m_np + sum(mass(:,jm)-m_np))*byim
        rm(:,jm) = (rm_np + sum(rm(:,jm)-rm_np))*byim
        rmom(mz ,:,jm) = (rzm_np + 
     &       sum(rmom(mz ,:,jm)-rzm_np ))*byim
        rmom(mzz,:,jm) = (rzzm_np+ 
     &       sum(rmom(mzz,:,jm)-rzzm_np))*byim
        rmom(ihmoms,:,jm) = 0.
      end if  !NORTH POLE

      return
      end subroutine aadvqy2

      subroutine aadvqz2(rm,rmom,mass,mw,mwdn,fdn,fmomdn
     &     ,scf,scm,sfcm)
!@sum  AADVQZ2 version of AADVQZ without flux limiter
!@auth Maxwell Kelley
      use DOMAIN_DECOMP, only : grid, GET
      use QUSDEF
      use QUSCOM, only : im
      implicit none
      REAL*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,2) 
     &        :: rm,mass
      REAL*8, dimension(NMOM,IM,grid%j_strt_halo:grid%j_stop_halo,2) 
     &        :: rmom
      REAL*8, dimension(NMOM,IM,grid%j_strt_halo:grid%j_stop_halo) 
     &        :: fmomdn
      REAL*8, dimension(IM,grid%j_strt_halo:grid%j_stop_halo) ::
     &     mw,mwdn,fdn
      REAL*8, dimension(GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     sfcm,scm,scf
      real*8, dimension(nmom) :: fmomup
      integer :: i,j,l,ll
      real*8 :: frac1,fracm,fup,mold,mnew,bymnew,dm2

c**** Get useful local parameters for domain decomposition
      integer :: J_0, J_1
      CALL GET( grid, J_STRT=J_0 , J_STOP=J_1 )

      do j=j_0,j_1
      do i=1,im
         if(mw(i,j).lt.0.) then ! air mass flux is negative
            ll=2
            frac1=+1.
          else                   ! air mass flux is positive
            ll=1
            frac1=-1.
         endif
         fracm=mw(i,j)/mass(i,j,ll)
         frac1=fracm+frac1
         fup=fracm*(rm(i,j,ll)-frac1*(rmom(mz,i,j,ll)-
     &        (frac1+fracm)*rmom(mzz,i,j,ll)))
         fmomup(mz)=mw(i,j)*(fracm*fracm*(rmom(mz,i,j,ll)
     &        -3.*frac1*rmom(mzz,i,j,ll))-3.*fup)
         fmomup(mzz)=mw(i,j)*(mw(i,j)*fracm**3 *rmom(mzz,i,j,ll)
     &        -5.*(mw(i,j)*fup+fmomup(mz)))
         fmomup(my)  = fracm*(rmom(my,i,j,ll)-frac1*rmom(myz,i,j,ll))
         fmomup(myz) = mw(i,j)*
     &        (fracm*fracm*rmom(myz,i,j,ll)-3.*fmomup(my))
         fmomup(mx)  = fracm*(rmom(mx,i,j,ll)-frac1*rmom(mzx,i,j,ll))
         fmomup(mzx) = mw(i,j)*
     &        (fracm*fracm*rmom(mzx,i,j,ll)-3.*fmomup(mx))
         fmomup(myy) = fracm*rmom(myy,i,j,ll)
         fmomup(mxx) = fracm*rmom(mxx,i,j,ll)
         fmomup(mxy) = fracm*rmom(mxy,i,j,ll)

         mold=mass(i,j,1)
         mnew=mold+mwdn(i,j)-mw(i,j)
         bymnew = 1./mnew
         dm2=mwdn(i,j)+mw(i,j)
         rm(i,j,1)=rm(i,j,1)+fdn(i,j)-fup
      !
         rmom(mz,i,j,1)=(rmom(mz,i,j,1)*mold-3.*(-dm2*rm(i,j,1)
     &     +mold*(fdn(i,j)+fup))+(fmomdn(mz,i,j)-fmomup(mz)))*bymnew
         rmom(mzz,i,j,1) = (rmom(mzz,i,j,1)*mold*mold
     &     +2.5*rm(i,j,1)*(mold*mold-mnew*mnew-3.*dm2*dm2)
     &     +5.*(mold*(mold*(fdn(i,j)-fup)-fmomdn(mz,i,j)
     &     -fmomup(mz))+dm2*rmom(mz,i,j,1)*mnew)
     &     +(fmomdn(mzz,i,j)-fmomup(mzz))) * (bymnew*bymnew)
      ! cross moments
         rmom(my,i,j,1)=rmom(my,i,j,1)+fmomdn(my,i,j)-fmomup(my)
         rmom(myz,i,j,1)=(rmom(myz,i,j,1)*mold-3.*(-dm2*rmom(my,i,j,1) +
     &        mold*(fmomdn(my,i,j)+fmomup(my))) +
     &        (fmomdn(myz,i,j)-fmomup(myz)))*bymnew
         rmom(mx,i,j,1)=rmom(mx,i,j,1)+fmomdn(mx,i,j)-fmomup(mx)
         rmom(mzx,i,j,1)=(rmom(mzx,i,j,1)*mold-3.*(-dm2*rmom(mx,i,j,1) +
     &        mold*(fmomdn(mx,i,j)+fmomup(mx))) +
     &        (fmomdn(mzx,i,j)-fmomup(mzx)))*bymnew
      !
         rmom(myy,i,j,1)=rmom(myy,i,j,1)+fmomdn(myy,i,j)-fmomup(myy)
         rmom(mxx,i,j,1)=rmom(mxx,i,j,1)+fmomdn(mxx,i,j)-fmomup(mxx)
         rmom(mxy,i,j,1)=rmom(mxy,i,j,1)+fmomdn(mxy,i,j)-fmomup(mxy)

         mass(i,j,1) = mnew

         mwdn(i,j) = mw(i,j)
         fdn(i,j) = fup
         fmomdn(:,i,j) = fmomup(:)

         scf(j) = scf(j) + fup
         scm(j) = scm(j) + mw(i,j)
         if(mw(i,j).ne.0.) sfcm(j) = sfcm(j) + fup/mw(i,j)

      enddo ! i
      enddo ! j

      return
c****
      end subroutine aadvqz2
