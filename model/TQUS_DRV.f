      MODULE TRACER_ADV
!@sum MODULE TRACER_ADV arrays needed for tracer advection

      USE MODEL_COM, ONLY : IM,JM,LM,byim
      SAVE
      INTEGER, PARAMETER :: ncmax=10
      INTEGER NSTEPX1(JM,LM,NCMAX),NSTEPX2(JM,LM,NCMAX),
     *        NSTEPY(LM,NCMAX),NSTEPZ(IM*JM,NCMAX), NCYC
      double precision, dimension(jm,lm) :: sfbm,sbm,sbf,sfcm,scm,scf

      contains

      SUBROUTINE AADVQ (RM,RMOM,QLIMIT,tname)
!@sum  AADVQ advection driver
!@auth G. Russell, modified by Maxwell Kelley
!@+        Jean Lerner modified this for tracers in mass units
c****
c**** AADVQ advects tracers using the Quadradic Upstream Scheme.
c****
c**** input:
c****  pu,pv,sd (kg/s) = east-west,north-south,vertical mass fluxes
c****      qlimit = whether moment limitations should be used
C****         DT (s) = time step
c****
c**** input/output:
c****     rm = tracer mass
c****   rmom = moments of tracer mass
c****     ma (kg) = fluid mass
c****
      USE QUSCOM, ONLY : MFLX,nmom
      USE DYNAMICS, ONLY: pu=>pua, pv=>pva, sd=>sda, mb,ma
      IMPLICIT NONE

      double precision, dimension(im,jm,lm) :: rm
      double precision, dimension(nmom,im,jm,lm) :: rmom
      logical, intent(in) :: qlimit
      character*8 tname          !tracer name
      integer :: I,J,L,n,nx

C**** Fill in values at the poles
C$OMP  PARALLEL DO PRIVATE (I,L)
      do l=1,lm
         rm(2:im,1 ,l) =   rm(1,1 ,l)
         rm(2:im,jm,l) =   rm(1,jm,l)
         do i=2,im
            rmom(:,i,1 ,l) =  rmom(:,1,1 ,l)
            rmom(:,i,jm,l) =  rmom(:,1,jm,l)
         enddo
      enddo
C$OMP  END PARALLEL DO
C****
C**** Load mass after advection from mass before advection
C****
ccc   ma(:,:,:) = mb(:,:,:)
C$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         MA(:,:,L) = MB(:,:,L)
      ENDDO
C$OMP  END PARALLEL DO
C****
C**** Advect the tracer using the quadratic upstream scheme
C****
C**** loop over cycles
      do n=1,ncyc
ccc   mflx(:,:,:)=pu(:,:,:)
C$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         MFLX(:,:,L) = PU(:,:,L)
      ENDDO
C$OMP  END PARALLEL DO

      call aadvqx (rm,rmom,ma,mflx,qlimit,tname,nstepx1(1,1,n))

ccc   mflx(:,:,:)=pv(:,:,:)
C$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         MFLX(:,:,L) = PV(:,:,L)
      ENDDO
C$OMP  END PARALLEL DO

      call aadvqy (rm,rmom,ma,mflx,qlimit,tname,nstepy(1,n),
     &    sbf,sbm,sfbm)

ccc   mflx(:,:,:)=sd(:,:,:)
C$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         MFLX(:,:,L) = SD(:,:,L)
      ENDDO
C$OMP  END PARALLEL DO

      call aadvqz (rm,rmom,ma,mflx,qlimit,tname,nstepz(1,n),
     &    scf,scm,sfcm)

ccc   mflx(:,:,:)=pu(:,:,:)
C$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         MFLX(:,:,L) = PU(:,:,L)
      ENDDO
C$OMP  END PARALLEL DO

      call aadvqx (rm,rmom,ma,mflx,qlimit,tname,nstepx2(1,1,n))
      end do

      return
      end SUBROUTINE AADVQ


      SUBROUTINE AADVQ0(DT)
!@sum AADVQ0 initialises advection of tracer.
!@+   Decide how many cycles to take such that mass does not become
!@+   negative during any of the operator splitting steps of each cycle
c****
C**** The MA array space is temporarily put to use in this section
      USE DYNAMICS, ONLY: mu=>pua, mv=>pva, mw=>sda, mb, ma
      USE QUSCOM, ONLY : IM,JM,LM
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: DT
      LOGICAL :: mneg
      INTEGER :: i,j,l,n,nc,im1
      double precision :: byn,ssp,snp

ccc   mu(:,:,:) = mu(:,:,:)*(.5*dt)
ccc   mv(:,1:jm-1,:) = mv(:,2:jm,:)*dt
ccc   mv(:,jm,:) = 0.
ccc   mw(:,:,1:lm-1) = mw(:,:,1:lm-1)*(-dt)
C
C$OMP  PARALLEL DO PRIVATE (J,L)
      DO L=1,LM
         MU(:,1,L) = 0.
         DO J=2,JM-1
            MU(:,J,L) = MU(:,J,L)*(.5*DT)
         ENDDO
         MU(:,JM,L) = 0.
      ENDDO
C$OMP  END PARALLEL DO
C
C$OMP  PARALLEL DO PRIVATE (J,L)
      DO L=1,LM
         DO J=1,JM-1
            MV(:,J,L) = MV(:,J+1,L)*DT
         ENDDO
         MV(:,JM,L) = 0.
      ENDDO
C$OMP  END PARALLEL DO
C
C$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM-1
         MW(:,:,L) = MW(:,:,L)*(-DT)
      ENDDO
C$OMP  END PARALLEL DO
C
c     for some reason mu is not zero at the poles...
ccc   mu(:,1,:) = 0.
ccc   mu(:,jm,:) = 0.
C**** Set things up
      mneg=.true.
      ncyc = 0
      do while(mneg)
      ncyc = ncyc + 1
      byn = 1./ncyc
      mneg=.false.
      ma(:,:,:) = mb(:,:,:)
      do nc=1,ncyc
C****     1/2 x-direction
        do l=1,lm
        do j=1,jm
          im1 = im
          do i=1,im
            ma(i,j,l) = ma(i,j,l) + (mu(im1,j,l)-mu(i,j,l))*byn
            if (ma(i,j,l)/mb(i,j,l).lt.0.5) then
    !         write(6,900)'ma(i,j,l) x1: nc=',i,j,l,ma(i,j,l)/mb(i,j,l)
    !*             ,nc
              mneg=.true.
              go to 595
            endif
            im1 = i
          end do
        end do
        end do
C****         y-direction
        do l=1,lm              !Interior
        do j=2,jm-1
        do i=1,im
          ma(i,j,l) = ma(i,j,l) + (mv(i,j-1,l)-mv(i,j,l))*byn
          if (ma(i,j,l)/mb(i,j,l).lt.0.5) then
    !       write(6,900)'ma(i,j,l) y: nc=',i,j,l,ma(i,j,l)/mb(i,j,l),nc
            mneg=.true.
            go to 595
          endif
        end do
        end do
        end do
        do l=1,lm              !Poles
          ssp = sum(ma(:, 1,l)-mv(:,   1,l)*byn)*byim
          snp = sum(ma(:,jm,l)+mv(:,jm-1,l)*byn)*byim
          ma(:,1 ,l) = ssp
          ma(:,jm,l) = snp
          if (ma(1,1,l)/mb(1,1,l).lt.0.5) then
    !       write(6,910)'ma(1,1,l) ysp: nc=',l,ma(1,1,l)/mb(1,1,l),nc
            mneg=.true.
            go to 595
          endif
          if (ma(1,jm,l)/mb(1,jm,l).lt.0.5) then
    !       write(6,910)'ma(1,jm,l) ynp: nc=',l,ma(1,jm,l)/mb(1,jm,l),nc
            mneg=.true.
            go to 595
          endif
        end do
C****         z-direction
        do l=2,lm-1
        do j=1,jm
        do i=1,im
          ma(i,j,l) = ma(i,j,l)+(mw(i,j,l-1)-mw(i,j,l))*byn
          if (ma(i,j,l)/mb(i,j,l).lt.0.5) then
    !       write(6,900)'ma(i,j,l) z : nc=',i,j,l,ma(i,j,l)/mb(i,j,l),nc
            mneg=.true.
            go to 595
          endif
        end do
        end do
        end do
        l = 1
        do j=1,jm
        do i=1,im
          ma(i,j,l) = ma(i,j,l)-mw(i,j,l)*byn
          if (ma(i,j,l)/mb(i,j,l).lt.0.5) then
    !       write(6,900)'ma(i,j,l) z1 : nc=',i,j,l,ma(i,j,l)/mb(i,j,l)
    !*           ,nc
            mneg=.true.
            go to 595
          endif
        end do
        end do
        l = lm
        do j=1,jm
        do i=1,im
          ma(i,j,l) = ma(i,j,l)+mw(i,j,l-1)*byn
          if (ma(i,j,l)/mb(i,j,l).lt.0.5) then
    !       write(6,900)'ma(i,j,l) zlm: nc=',i,j,l,ma(i,j,l)/mb(i,j,l)
    !*           ,nc
            mneg=.true.
            go to 595
          endif
        end do
        end do
C****     1/2 x-direction
        do l=1,lm
        do j=1,jm
          im1 = im
          do i=1,im
            ma(i,j,l) = ma(i,j,l) + (mu(im1,j,l)-mu(i,j,l))*byn
            if (ma(i,j,l)/mb(i,j,l).lt.0.5) then
    !         write(6,900)'ma(i,j,l) x2: nc=',i,j,l,ma(i,j,l)/mb(i,j,l)
    !*             ,nc
              mneg=.true.
              go to 595
            endif
            im1 = i
          end do
        end do
        end do
      end do
 595  continue
      if(ncyc.ge.10) then
        write(6,*) 'stop: ncyc=10 in AADVQ0'
        call exit_rc(11)
      end if
 600  enddo
      if(ncyc.gt.1) write(6,*) 'AADVQ0: ncyc>1',ncyc
C**** Divide the mass fluxes by the number of cycles
      byn = 1./ncyc
ccc   mu(:,:,:)=mu(:,:,:)*byn
ccc   mv(:,1:jm-1,:)=mv(:,1:jm-1,:)*byn
ccc   mw(:,:,:)=mw(:,:,:)*byn
C$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         mu(:,:,l)=mu(:,:,l)*byn
      ENDDO
C$OMP  END PARALLEL DO
C$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         mv(:,1:jm-1,l)=mv(:,1:jm-1,l)*byn
      ENDDO
C$OMP  END PARALLEL DO
C$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         mw(:,:,l)=mw(:,:,l)*byn
      ENDDO
C$OMP  END PARALLEL DO
C****
C**** Decide how many timesteps to take by computing Courant limits
C****
ccc   MA(:,:,:) = MB(:,:,:)
C$OMP  PARALLEL DO PRIVATE (L)
      DO L=1,LM
         MA(:,:,L) = MB(:,:,L)
      ENDDO
C$OMP  END PARALLEL DO
      do n=1,ncyc
        call xstep (MA,nstepx1(1,1,n))
        call ystep (MA,nstepy(1,n))
        call zstep (MA,nstepz(1,n))
        call xstep (MA,nstepx2(1,1,n))
      end do
      RETURN
  900 format (1x,a,3i4,f10.4,i5)
  910 format (1x,a,i4,f10.4,i5)
      END subroutine AADVQ0

      end MODULE TRACER_ADV


      subroutine aadvQx(rm,rmom,mass,mu,qlimit,tname,nstep)
!@sum  AADVQX advection driver for x-direction
!@auth Maxwell Kelley; modified by J. Lerner
c****
c**** aadvtQ advects tracers in the west to east direction using the
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
ccc   use QUSCOM, only : im,jm,lm, xstride,am,f_i,fmom_i
      use QUSCOM, only : im,jm,lm, xstride
      use QUSDEF
      implicit none
      double precision, dimension(im,jm,lm) :: rm,mass,mu
      double precision, dimension(nmom,im,jm,lm) :: rmom
      logical ::  qlimit
      DOUBLE PRECISION  AM(IM), F_I(IM), FMOM_I(NMOM,IM)
      character*8 tname
      integer :: i,j,l,ierr,nerr,ns,nstep(jm,lm),ICKERR

c**** loop over layers and latitudes
      ICKERR=0
C$OMP  PARALLEL DO PRIVATE (J,L,NS,AM,F_I,FMOM_I,IERR,NERR)
      do l=1,lm
      do j=2,jm-1
      am(:) = mu(:,j,l)/nstep(j,l)
c****
c**** call 1-d advection routine
c****
      do ns=1,nstep(j,l)
      call adv1d(rm(1,j,l),rmom(1,1,j,l), f_i,fmom_i, mass(1,j,l),
     &        am, im, qlimit,xstride,xdir,ierr,nerr)
      if (ierr.gt.0) then
        write(6,*) "Error in aadvQx: i,j,l=",nerr,j,l,' ',tname
        if (ierr.eq.2) write(6,*) "Error in qlimit: abs(a) > 1"
        if (ierr.eq.2) ICKERR=1
      end if
      enddo ! ns
      enddo ! j
      enddo ! l
C$OMP  END PARALLEL DO
C
      IF(ICKERR.GT.0)  CALL EXIT_RC(11)
C
      return
c****
      end subroutine aadvQx


      subroutine aadvQy(rm,rmom,mass,mv,qlimit,tname,nstep,
     &   sbf,sbm,sfbm)
!@sum  AADVQY advection driver for y-direction
!@auth Maxwell Kelley; modified by J. Lerner
c****
c**** aadvQy advects tracers in the south to north direction using the
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
      use CONSTANT, only : teeny
ccc   use QUSCOM, only : im,jm,lm, ystride,bm,f_j,fmom_j, byim
      use QUSCOM, only : im,jm,lm, ystride,               byim
      use QUSDEF
      implicit none
      double precision, dimension(im,jm,lm) :: rm,mass,mv
      double precision, dimension(nmom,im,jm,lm) :: rmom
      logical ::  qlimit
      double precision, dimension(im,jm) :: fqv
      double precision, intent(out), dimension(jm,lm) :: sfbm,sbm,sbf
      DOUBLE PRECISION  BM(JM),F_J(JM),FMOM_J(NMOM,JM)
      character*8 tname
      integer :: i,j,l,ierr,nerr,ns,nstep(lm),ICKERR
      double precision ::
     &     m_sp,m_np,rm_sp,rm_np,rzm_sp,rzm_np,rzzm_sp,rzzm_np

c**** loop over layers
      ICKERR=0
C$OMP  PARALLEL DO PRIVATE (I,J,L,M_SP,M_NP,RM_SP,RM_NP,RZM_SP,RZM_NP,
C$OMP*                RZZM_SP,RZZM_NP,BM,F_J,FMOM_J,FQV,NS,IERR,NERR)
      do l=1,lm
      fqv(:,:) = 0.

c**** loop over timesteps
      do ns=1,nstep(l)

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
      bm (:) = mv(i,:,l)/nstep(l)
      bm(jm) = 0.
      rmom(ihmoms,i,1 ,l) = 0.! horizontal moments are zero at pole
      rmom(ihmoms,i,jm,l) = 0.
c****
c**** call 1-d advection routine
c****
      call adv1d(rm(i,1,l),rmom(1,i,1,l), f_j,fmom_j, mass(i,1,l),
     &     bm, jm,qlimit,ystride,ydir,ierr,nerr)
      if (ierr.gt.0) then
        write(6,*) "Error in aadvQy: i,j,l=",i,nerr,l,' ',tname
        if (ierr.eq.2) write(6,*) "Error in qlimit: abs(b) > 1"
        if (ierr.eq.2) ICKERR=1
      end if
      fqv(i,:) = fqv(i,:) + f_j(:)  !store tracer flux in fqv array
      fqv(i,jm) = 0.   ! play it safe
      rmom(ihmoms,i,1 ,l) = 0.! horizontal moments are zero at pole
      rmom(ihmoms,i,jm,l) = 0.
c     sbfijl(i,:,l) = sbfijl(i,:,l)+f_j(:)
      enddo  ! end loop over longitudes
c**** average and unscale polar boxes
      mass(:,1 ,l) = (m_sp + sum(mass(:,1 ,l)-m_sp))*byim
      mass(:,jm,l) = (m_np + sum(mass(:,jm,l)-m_np))*byim
      rm(:,1 ,l) = (rm_sp + sum(rm(:,1 ,l)-rm_sp))*byim
      rm(:,jm,l) = (rm_np + sum(rm(:,jm,l)-rm_np))*byim
      rmom(mz ,:,1 ,l) = (rzm_sp  + sum(rmom(mz ,:,1 ,l)-rzm_sp ))*byim
      rmom(mzz,:,1 ,l) = (rzzm_sp + sum(rmom(mzz,:,1 ,l)-rzzm_sp))*byim
      rmom(mz ,:,jm,l) = (rzm_np  + sum(rmom(mz ,:,jm,l)-rzm_np ))*byim
      rmom(mzz,:,jm,l) = (rzzm_np + sum(rmom(mzz,:,jm,l)-rzzm_np))*byim

      enddo  ! end loop over timesteps

      do j=1,jm-1   !diagnostics
        sfbm(j,l) = sfbm(j,l) + sum(fqv(:,j)/(mv(:,j,l)+teeny))
        sbm (j,l) = sbm (j,l) + sum(mv(:,j,l))
        sbf (j,l) = sbf (j,l) + sum(fqv(:,j))
      enddo

      enddo  ! end loop over levels
C$OMP  END PARALLEL DO
C
      IF(ICKERR.NE.0)  call exit_rc(11)
C
      return
c****
      end subroutine aadvQy


      subroutine aadvQz(rm,rmom,mass,mw,qlimit,tname,nstep,
     &  scf,scm,sfcm)
!@sum  AADVQZ advection driver for z-direction
!@auth Maxwell Kelley; modified by J. Lerner
c****
c**** aadvQz advects tracers in the upward vertical direction using the
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
      use CONSTANT, only : teeny
ccc   use QUSCOM, only : im,jm,lm, zstride,cm,f_l,fmom_l
      use QUSCOM, only : im,jm,lm, zstride
      use QUSDEF
      use GEOM, only : imaxj
      use MODEL_COM, only : fim
      implicit none
      double precision, dimension(im,jm,lm) :: rm,mass,mw
      double precision, dimension(nmom,im,jm,lm) :: rmom
      integer, dimension(im,jm) :: nstep
      logical ::  qlimit
      double precision, dimension(lm) :: fqw
      double precision, intent(out), dimension(jm,lm) :: sfcm,scm,scf
      character*8 tname
      DOUBLE PRECISION  CM(LM),F_L(LM),FMOM_L(NMOM,LM)
      integer :: i,j,l,ierr,nerr,ns,ICKERR

c**** loop over latitudes and longitudes
      ICKERR=0.
C$OMP  PARALLEL DO PRIVATE (I,J,L,NS,FQW,CM,F_L,FMOM_L,FQW,IERR,NERR)
      do j=1,jm
      do i=1,imaxj(j)
      fqw(:) = 0.
      cm(:) = mw(i,j,:)/nstep(i,j)
      cm(lm)= 0.
c****
c**** call 1-d advection routine
c****
      do ns=1,nstep(i,j)
      call adv1d(rm(i,j,1),rmom(1,i,j,1),f_l,fmom_l,mass(i,j,1),
     &        cm,lm,qlimit,zstride,zdir,ierr,nerr)
      if (ierr.gt.0) then
        write(6,*) "Error in aadvQz: i,j,l=",i,j,nerr,' ',tname
        if (ierr.eq.2) write(6,*) "Error in qlimit: abs(c) > 1"
        if (ierr.eq.2) ICKERR=1
      end if
      fqw(:)  = fqw(:) + f_l(:) !store tracer flux in fqw array
      enddo ! ns
      do l=1,lm-1   !diagnostics
        sfcm(j,l) = sfcm(j,l) + fqw(l)/(mw(i,j,l)+teeny)
        scm (j,l) = scm (j,l) + mw(i,j,l)
        scf (j,l) = scf (j,l) + fqw(l)
      enddo
      enddo ! i
      if (j.eq.1.or.j.eq.jm) then
        do l=1,lm
          do i=2,im
            rm(i,j,l)=rm(1,j,l)
            rmom(:,i,j,l)=rmom(:,1,j,l)
            mass(i,j,l)=mass(1,j,l)
          end do
        end do
        sfcm(j,:) = fim*sfcm(j,:)
        scm(j,:)  = fim*scm(j,:)
        scf(j,:)  = fim*scf(j,:)
      end if
      enddo ! j
C$OMP  END PARALLEL DO
C
      IF(ICKERR.GT.0)  call exit_rc(11)
C
      return
c****
      end subroutine aadvQz



      SUBROUTINE XSTEP (M,NSTEPX)
!@sum XSTEP determines the number of X timesteps for tracer dynamics
!@+    using Courant limits
!@auth J. Lerner and M. Kelley
!@ver  1.0
      USE QUSCOM, ONLY : IM,JM,LM,byim
      USE DYNAMICS, ONLY: mu=>pua
      IMPLICIT NONE
      double precision, dimension(im,jm,lm) :: m
      double precision, dimension(im) :: a,am,mi
      integer, dimension(jm,lm) :: nstepx
      integer :: l,j,i,ip1,im1,nstep,ns,ICKERR
      double precision :: courmax

C**** Decide how many timesteps to take by computing Courant limits
C
      ICKERR = 0
C$OMP  PARALLEL DO PRIVATE (I,IP1,IM1,J,L,NSTEP,NS,COURMAX,A,AM,MI)
      DO 420 L=1,LM
      DO 420 J=2,JM-1
      nstep=0
      courmax = 2.
      do while(courmax.gt.1.)
        nstep = nstep+1   !(1+int(courmax))
        am(:) = mu(:,j,l)/nstep
        mi(:) = m (:,j,l)
        courmax = 0.
        do ns=1,nstep
          i = im
          do ip1=1,im
            if(am(i).gt.0.) then
               a(i) = am(i)/mi(i)
               courmax = max(courmax,+a(i))
            else
               a(i) = am(i)/mi(ip1)
               courmax = max(courmax,-a(i))
            endif
          i = ip1
          enddo  ! ip1=1,im
C**** Update air mass
          im1 = im
          do i=1,im
            mi(i) = mi(i)+am(im1)-am(i)
          im1 = i
          enddo
        enddo    ! ns=1,nstep
        if(nstep.ge.20) write(6,*) 'aadvqx: nstep.ge.20'
        if(nstep.ge.20)  then
           write(6,*) 'aadvqx: j,l,nstep,courmax=',j,l,nstep,courmax
           courmax=-1.
           ICKERR=1
        end if
      enddo      ! while(courmax.gt.1.)
C**** Correct air mass
      M(:,J,L) = MI(:)
      NSTEPX(J,L) = NSTEP
c     if(nstep.gt.2 .and. nx.eq.1)
c    *  write(6,'(a,3i3,f7.4)')
c    *  'aadvqx: j,l,nstep,courmax=',j,l,nstep,courmax
  420 CONTINUE
C$OMP  END PARALLEL DO
C
      IF(ICKERR.GT.0)  call exit_rc(11)
C
      RETURN
      END


      SUBROUTINE YSTEP (M,NSTEPY)
!@sum YSTEP determines the number of Y timesteps for tracer dynamics
!@+    using Courant limits
!@auth J. Lerner and M. Kelley
!@ver  1.0
      USE QUSCOM, ONLY : IM,JM,LM,byim
      USE DYNAMICS, ONLY: mv=>pva
      IMPLICIT NONE
      double precision, dimension(im,jm,lm) :: m
      double precision, dimension(im,jm) :: mij
      double precision, dimension(jm) :: b,bm
      integer, dimension(lm) :: nstepy
      integer :: jprob,iprob,nstep,ns,i,j,l,ICKERR
      double precision :: courmax,byn,sbms,sbmn

C**** decide how many timesteps to take (all longitudes at this level)
      ICKERR=0
C$OMP  PARALLEL DO PRIVATE (I,J,L,NS,NSTEP,COURMAX,BYN,B,BM,MIJ,
C$OMP*          IPROB,JPROB,SBMS,SBMN)
      DO 440 L=1,LM
C**** Scale poles
      m(:, 1,l) =   m(:, 1,l)*im !!!!! temporary
      m(:,jm,l) =   m(:,jm,l)*im !!!!! temporary
C**** begin computation
      nstep=0
      courmax = 2.
      do while(courmax.gt.1.)
        nstep = nstep+1   !(1+int(courmax))
        byn = 1./nstep
        courmax = 0.
        mij(:,:) = m(:,:,l)
        do ns=1,nstep
          do j=1,jm-1
          do i=1,im
            bm(j) = mv(i,j,l)*byn
            if(bm(j).gt.0.) then
               b(j) = bm(j)/mij(i,j)
               courmax = max(courmax,+b(j))
               if (courmax.eq.b(j)) jprob = j
               if (courmax.eq.b(j)) iprob = i
            else
               b(j) = bm(j)/mij(i,j+1)
               courmax = max(courmax,-b(j))
               if (courmax.eq.-b(j)) jprob = j
               if (courmax.eq.-b(j)) iprob = i
            endif
          enddo
          enddo
C**** Update air mass at poles
          sbms = sum(mv(:,   1,l))*byn
          sbmn = sum(mv(:,jm-1,l))*byn
          mij(:, 1) = mij(:, 1)-sbms
          mij(:,jm) = mij(:,jm)+sbmn
C**** Update air mass in the interior
          do j=2,jm-1
            mij(:,j) = mij(:,j)+(mv(:,j-1,l)-mv(:,j,l))*byn
          enddo
        enddo    ! ns=1,nstep
        if(nstep.ge.20) then
           write(6,*) 'courmax=',courmax,l,iprob,jprob
           write(6,*) 'aadvqy: nstep.ge.20'
           ICKERR=1
           courmax = -1.
        endif
      enddo      ! while(courmax.gt.1.)
C**** Correct air mass
      m(:,:,l) = mij(:,:)
      NSTEPY(L) = nstep
c       if(nstep.gt.1. and. nTRACER.eq.1) write(6,'(a,2i3,f7.4)')
c    *    'aadvqy: l,nstep,courmax=',l,nstep,courmax
C**** Unscale poles
      m(:, 1,l) =   m(:, 1,l)*byim !!! undo temporary
      m(:,jm,l) =   m(:,jm,l)*byim !!! undo temporary
  440 CONTINUE
C$OMP  END PARALLEL DO
C
      IF(ICKERR.GT.0)  call exit_rc(11)
C
      RETURN
      END


      SUBROUTINE ZSTEP (M,NSTEPZ)
!@sum ZSTEP determines the number of Z timesteps for tracer dynamics
!@+    using Courant limits
!@auth J. Lerner and M. Kelley
!@ver  1.0
      USE QUSCOM, ONLY : IM,JM,LM,byim
      USE DYNAMICS, ONLY: mw=>sda
      IMPLICIT NONE
      double precision, dimension(im,jm,lm) :: m
      double precision, dimension(lm) :: ml
      double precision, dimension(0:lm) :: c,cm
      integer, dimension(im*jm) :: nstepz
      integer :: nstep,ns,l,i,ICKERR
      double precision :: courmax,byn

C**** decide how many timesteps to take
      ICKERR=0
C$OMP  PARALLEL DO PRIVATE (I,L,NS,NSTEP,COURMAX,BYN,C,CM,ML)
      DO I=1,IM*JM
      nstep=0
      courmax = 2.
      do while(courmax.gt.1.)
        nstep = nstep+1   !(1+int(courmax))
        byn = 1./nstep
        cm(1:lm) = mw(i,1,1:lm)*byn
        ml(:)  = m(i,1,:)
        CM(LM)= 0. ! VERY IMPORTANT TO SET THIS TO ZERO
        CM( 0)= 0. ! VERY IMPORTANT TO SET THIS TO ZERO
        courmax = 0.
        do ns=1,nstep
          do l=1,lm-1
            if(cm(l).gt.0.) then
               c(l) = cm(l)/ml(l)
               courmax = max(courmax,+c(l))
            else
               c(l) = cm(l)/ml(l+1)
               courmax = max(courmax,-c(l))
            endif
          enddo
          do l=1,lm
            ml(l) = ml(l)+(cm(l-1)-cm(l))
          enddo
        enddo    ! ns=1,nstep
        if(nstep.ge.20) write(6,*) 'aadvqz: nstep.ge.20'
        if(nstep.ge.20)  then
           write(6,*)  'aadvqz: nstep.ge.20'
           ICKERR=1
           courmax = -1.
        end if
      enddo      ! while(courmax.gt.1.)
C**** Correct air mass
      m(i,1,:) = ml(:)
      NSTEPZ(I) = NSTEP
c     if(nstep.gt.1 .and. nTRACER.eq.1) write(6,'(a,2i7,f7.4)')
c    *   'aadvqz: i,nstep,courmax=',i,nstep,courmax
      END DO
C$OMP  END PARALLEL DO
C
      IF(ICKERR.GT.0)  call exit_rc(11)
C
      RETURN
      END
