C**** QUSEM12 E001M12 SOMTQ QUSB261AM12
C**** QUSBM9=QUSB140M9 with correction to second order moment calc.
C**** Changes for constant pressure above LS1 + double precision
C**** QUS   is Russell quadratic upstream scheme for temperature
C**** and water vapor advection, with limits applied to water vapor.
C**** Changes for constant pressure above LS1
C**** FQU,FQV for additional diagnostics
C**** Routines included: AADVT, AADVTX, AADVTY, AADVTZ

      MODULE QUSCOM
!@sum  QUSCOM contains gcm-specific advection parameters/workspace
!@auth Maxwell Kelley
      USE QUSDEF
      IMPLICIT NONE
      SAVE
      INTEGER :: IM,JM,LM
      INTEGER :: XSTRIDE,YSTRIDE,ZSTRIDE
      DOUBLE PRECISION :: BYIM
C**** AIR MASS FLUXES
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: MFLX
C**** WORKSPACE FOR AADVTX
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: AM,F_I
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FMOM_I
C**** WORKSPACE FOR AADVTY
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: BM,F_J
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FMOM_J
C**** WORKSPACE FOR AADVTZ
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: CM,F_L
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FMOM_L

      END MODULE QUSCOM

      SUBROUTINE init_QUS(IM_GCM,JM_GCM,LM_GCM)
!@sum  init_QUS sets gcm-specific advection parameters/workspace
!@auth Maxwell Kelley
      use QUSCOM
      INTEGER, INTENT(IN) :: IM_GCM,JM_GCM,LM_GCM
C**** SET RESOLUTION
      IM = IM_GCM
      JM = JM_GCM
      LM = LM_GCM
      BYIM = 1.D0/DBLE(IM)
      XSTRIDE = 1
      YSTRIDE = IM
      ZSTRIDE = IM*JM
C**** ALLOCATE SPACE FOR AIR MASS FLUXES
      ALLOCATE(MFLX(IM,JM,LM))
C**** ALLOCATE WORKSPACE FOR AADVTX
      ALLOCATE(AM(IM),F_I(IM),FMOM_I(NMOM,IM))
C**** ALLOCATE WORKSPACE FOR AADVTY
      ALLOCATE(BM(JM),F_J(JM),FMOM_J(NMOM,JM))
C**** ALLOCATE WORKSPACE FOR AADVTZ
      ALLOCATE(CM(LM),F_L(LM),FMOM_L(NMOM,LM))
      RETURN
      END SUBROUTINE init_QUS

      SUBROUTINE AADVT (MA,RM,RMOM,SD,PU,PV,DT,QLIMIT,FQU,FQV)
!@sum  AADVT advection driver
!@auth G. Russell, modified by Maxwell Kelley
c****
c**** AADVT advects tracers using the Quadradic Upstream Scheme.
c****
c**** input:
c****  pu,pv,sd (kg/s) = east-west,north-south,vertical mass fluxes
c****      qlimit = whether moment limitations should be used
C****         DT (s) = time step
c****
c**** input/output:
c****     rm = tracer concentration
c****   rmom = moments of tracer concentration
c****     ma (kg) = fluid mass
c****
      USE QUSCOM, ONLY : IM,JM,LM, MFLX
      USE QUSDEF
      IMPLICIT NONE

      double precision, dimension(im,jm,lm) :: rm,ma
      double precision, dimension(NMOM,IM,JM,LM) :: rmom

      REAL*8, INTENT(IN) :: DT
      double precision, dimension(im,jm,lm), intent(in) :: pu,pv
      double precision, dimension(im,jm,lm-1), intent(in) :: sd
      LOGICAL, INTENT(IN) :: QLIMIT

      double precision, dimension(im,jm), intent(inout) :: fqu,fqv

      INTEGER :: I,J,L
      DOUBLE PRECISION :: BYMA

C**** Initialise diagnostics
      FQU=0.  ; FQV=0.

C**** Fill in values at the poles
      DO L=1,LM
         RM(2:IM,1 ,L) =   RM(1,1 ,L)
         RM(2:IM,JM,L) =   RM(1,JM,L)
         DO I=2,IM
            RMOM(:,I,1 ,L) =  RMOM(:,1,1 ,L)
            RMOM(:,I,JM,L) =  RMOM(:,1,JM,L)
         enddo
      enddo
C****
C**** convert from concentration to mass units
C****
      DO L=1,LM
      DO J=1,JM
      DO I=1,IM
         RM(I,J,L)=RM(I,J,L)*MA(I,J,L)
         RMOM(:,I,J,L)=RMOM(:,I,J,L)*MA(I,J,L)
      enddo
      enddo
      enddo
C****
C**** Advect the tracer using the quadratic upstream scheme
C****
      mflx(:,:,:)=pu(:,:,:)*(.5*dt)
      CALL AADVTX (RM,RMOM,MA,MFLX,QLIMIT,FQU)
      mflx(:,1:jm-1,:)=pv(:,2:jm,:)*dt
      mflx(:,jm,:)=0.
      CALL AADVTY (RM,RMOM,MA,MFLX,QLIMIT,FQV)
      mflx(:,:,1:lm-1)=sd(:,:,1:lm-1)*(-dt)
      mflx(:,:,lm)=0.
      CALL AADVTZ (RM,RMOM,MA,MFLX,QLIMIT)
      mflx(:,:,:)=pu(:,:,:)*(.5*dt)
      CALL AADVTX (RM,RMOM,MA,MFLX,QLIMIT,FQU)
C****
C**** convert from mass to concentration units
C****
      DO L=1,LM
      DO J=1,JM
      DO I=1,IM
         BYMA = 1.D0/MA(I,J,L)
         RM(I,J,L)=RM(I,J,L)*BYMA
         RMOM(:,I,J,L)=RMOM(:,I,J,L)*BYMA
      enddo
      enddo
      enddo
      RETURN
      END

      subroutine aadvtx(rm,rmom,mass,mu,qlimit,fqu)
!@sum  AADVTX advection driver for x-direction
!@auth Maxwell Kelley
c****
c**** aadvtx advects tracers in the west to east direction using the
c**** quadratic upstream scheme.  if qlimit is true, the moments are
c**** limited to prevent the mean tracer from becoming negative.
c****
c**** input:
c****     mu (kg) = west-east mass flux, positive eastward
c****      qlimit = whether moment limitations should be used
c****
c**** input/output:
c****     rm (kg) = tracer mass
c****   rmom (kg) = moments of tracer mass
c****   mass (kg) = fluid mass
c****
      use QUSCOM, only : im,jm,lm, xstride,am,f_i,fmom_i
      use QUSDEF
      implicit none
      double precision, dimension(im,jm,lm) :: rm,mass,mu
      double precision, dimension(NMOM,IM,JM,LM) :: rmom
      logical ::  qlimit
      DOUBLE PRECISION, INTENT(OUT), DIMENSION(IM,JM) :: FQU
      integer :: i,j,l,ierr,nerr
c**** loop over layers and latitudes
      do l=1,lm
      do j=2,jm-1
      am(:) = mu(:,j,l)
c****
c**** call 1-d advection routine
c****
      call adv1d(rm(1,j,l),rmom(1,1,j,l), f_i,fmom_i, mass(1,j,l),
     &        am, im, qlimit,xstride,xdir,ierr,nerr)
      if (ierr.gt.0) then
        write(6,*) "Error in aadvtx: i,j,l=",nerr,j,l
        if (ierr.eq.2) then
          write(0,*) "Error in qlimit: abs(a) > 1"
          call exit_rc(11)
        end if
      end if
c****
c**** store tracer flux in fqu array
c****
      fqu(:,j)  = fqu(:,j) + f_i(:)
      enddo ! j
      enddo ! l
      return
c****
      end subroutine aadvtx

      subroutine aadvty(rm,rmom,mass,mv,qlimit,fqv)
!@sum  AADVTY advection driver for y-direction
!@auth Maxwell Kelley
c****
c**** aadvty advects tracers in the south to north direction using the
c**** quadratic upstream scheme.  if qlimit is true, the moments are
c**** limited to prevent the mean tracer from becoming negative.
c****
c**** input:
c****     mv (kg) = north-south mass flux, positive northward
c****      qlimit = whether moment limitations should be used
c****
c**** input/output:
c****     rm (kg) = tracer mass
c****   rmom (kg) = moments of tracer mass
c****   mass (kg) = fluid mass
c****
      use QUSCOM, only : im,jm,lm, ystride,bm,f_j,fmom_j, byim
      use QUSDEF
      implicit none
      double precision, dimension(im,jm,lm) :: rm,mass,mv
      double precision, dimension(NMOM,IM,JM,LM) :: rmom
      logical ::  qlimit
      double precision, intent(out), dimension(im,jm) :: fqv
      integer :: i,j,l,ierr,nerr
      double precision ::
     &     m_sp,m_np,rm_sp,rm_np,rzm_sp,rzm_np,rzzm_sp,rzzm_np
c**** loop over layers
      do l=1,lm
c**** scale polar boxes to their full extent
      mass(:,1:jm:jm-1,l)=mass(:,1:jm:jm-1,l)*im
      m_sp = mass(1,1 ,l)
      m_np = mass(1,jm,l)
      rm(:,1:jm:jm-1,l)=rm(:,1:jm:jm-1,l)*im
      rm_sp = rm(1,1 ,l)
      rm_np = rm(1,jm,l)
      do i=1,im
         rmom(:,i,1 ,l)=rmom(:,i,1 ,l)*im
         rmom(:,i,jm,l)=rmom(:,i,jm,l)*im
      enddo
      rzm_sp  = rmom(mz ,1,1 ,l)
      rzzm_sp = rmom(mzz,1,1 ,l)
      rzm_np  = rmom(mz ,1,jm,l)
      rzzm_np = rmom(mzz,1,jm,l)
c**** loop over longitudes
      do i=1,im
c****
c**** load 1-dimensional arrays
c****
      bm   (:) = mv(i,:,l) !/nstep
      bm(jm)= 0.
      rmom(ihmoms,i,1 ,l) = 0.! horizontal moments are zero at pole
      rmom(ihmoms,i,jm,l) = 0.
c****
c**** call 1-d advection routine
c****
      call adv1d(rm(i,1,l),rmom(1,i,1,l), f_j,fmom_j, mass(i,1,l),
     &     bm, jm,qlimit,ystride,ydir,ierr,nerr)
      if (ierr.gt.0) then
        write(6,*) "Error in aadvty: i,j,l=",i,nerr,l
        if (ierr.eq.2) then
          write(0,*) "Error in qlimit: abs(b) > 1"
          call exit_rc(11)
        endif
      end if
c**** store tracer flux in fqv array
      fqv(i,:) = fqv(i,:) + f_j(:)
      fqv(i,jm) = 0.   ! play it safe
      rmom(ihmoms,i,1 ,l) = 0.! horizontal moments are zero at pole
      rmom(ihmoms,i,jm,l) = 0.
      enddo ! end loop over longitudes
c**** average and unscale polar boxes
      mass(:,1 ,l) = (m_sp + sum(mass(:,1 ,l)-m_sp))*byim
      mass(:,jm,l) = (m_np + sum(mass(:,jm,l)-m_np))*byim
      rm(:,1 ,l) = (rm_sp + sum(rm(:,1 ,l)-rm_sp))*byim
      rm(:,jm,l) = (rm_np + sum(rm(:,jm,l)-rm_np))*byim
      rmom(mz ,:,1 ,l) = (rzm_sp  + sum(rmom(mz ,:,1 ,l)-rzm_sp ))*byim
      rmom(mzz,:,1 ,l) = (rzzm_sp + sum(rmom(mzz,:,1 ,l)-rzzm_sp))*byim
      rmom(mz ,:,jm,l) = (rzm_np  + sum(rmom(mz ,:,jm,l)-rzm_np ))*byim
      rmom(mzz,:,jm,l) = (rzzm_np + sum(rmom(mzz,:,jm,l)-rzzm_np))*byim
      enddo ! end loop over levels
      return
c****
      end subroutine aadvty


      subroutine aadvtz(rm,rmom,mass,mw,qlimit)
!@sum  AADVTZ advection driver for z-direction
!@auth Maxwell Kelley
c****
c**** aadvtz advects tracers in the upward vertical direction using the
c**** quadratic upstream scheme.  if qlimit is true, the moments are
c**** limited to prevent the mean tracer from becoming negative.
c****
c**** input:
c****     mw (kg) = vertical mass flux, positive upward
c****      qlimit = whether moment limitations should be used
c****
c**** input/output:
c****     rm (kg) = tracer mass
c****   rmom (kg) = moments of tracer mass
c****   mass (kg) = fluid mass
c****
      use QUSCOM, only : im,jm,lm, zstride,cm,f_l,fmom_l
      use QUSDEF
      implicit none
      double precision, dimension(im,jm,lm) :: rm,mass,mw
      double precision, dimension(NMOM,IM,JM,LM) :: rmom
      logical ::  qlimit
      integer :: i,j,l,ierr,nerr
c**** loop over latitudes and longitudes
      do j=1,jm
      do i=1,im
      cm(:) = mw(i,j,:)
      cm(lm)= 0.
c****
c**** call 1-d advection routine
c****
      call adv1d(rm(i,j,1),rmom(1,i,j,1),f_l,fmom_l,mass(i,j,1),
     &        cm,lm,qlimit,zstride,zdir,ierr,nerr)
      if (ierr.gt.0) then
        write(6,*) "Error in aadvtz: i,j,l=",i,j,nerr
        if (ierr.eq.2) then
          write(0,*) "Error in qlimit: abs(c) > 1"
          call exit_rc(11)
        endif
      end if
      enddo ! i
      enddo ! j
      return
c****
      end subroutine aadvtz
