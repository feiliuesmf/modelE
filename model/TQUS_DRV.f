      MODULE TRACER_ADV
!@sum MODULE TRACER_ADV arrays needed for tracer advection

      USE MODEL_COM, ONLY : IM,JM,LM,byim
      SAVE
      INTEGER, PARAMETER :: ncmax=10
      INTEGER NSTEPX1(JM,LM,NCMAX),NSTEPX2(JM,LM,NCMAX),
     *        NSTEPY(LM,NCMAX),NSTEPZ(IM*JM,NCMAX), NCYC
      DOUBLE PRECISION, DIMENSION(IM,JM,LM) :: MA
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
      USE DYNAMICS, ONLY: pu=>pua, pv=>pva, sd=>sda, mb
      IMPLICIT NONE

      double precision, dimension(im,jm,lm) :: rm
      double precision, dimension(nmom,im,jm,lm) :: rmom
      logical, intent(in) :: qlimit
      character*8 tname          !tracer name
      integer :: I,J,L,n,nx


C**** Fill in values at the poles
      do l=1,lm
         rm(2:im,1 ,l) =   rm(1,1 ,l)
         rm(2:im,jm,l) =   rm(1,jm,l)
         do i=2,im
            rmom(:,i,1 ,l) =  rmom(:,1,1 ,l)
            rmom(:,i,jm,l) =  rmom(:,1,jm,l)
         enddo
      enddo
C****
C**** Load mass after advection from mass before advection
C****
      ma(:,:,:) = mb(:,:,:)
C****
C**** Advect the tracer using the quadratic upstream scheme
C****
C**** loop over cycles
      do n=1,ncyc
c     write(*,*) ' Processing cycle X',n
      mflx(:,:,:)=pu(:,:,:)
      call aadvqx (rm,rmom,ma,mflx,qlimit,tname,nstepx1(1,1,n))

c     write(*,*) ' Processing cycle Y',n
      mflx(:,:,:)=pv(:,:,:)
      call aadvqy (rm,rmom,ma,mflx,qlimit,tname,nstepy(1,n),
     &    sbf,sbm,sfbm)

c     write(*,*) ' Processing cycle Z',n
      mflx(:,:,:)=sd(:,:,:)
      call aadvqz (rm,rmom,ma,mflx,qlimit,tname,nstepz(1,n),
     &    scf,scm,sfcm)

c     write(*,*) ' Processing cycle 2X',n
      mflx(:,:,:)=pu(:,:,:)
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
      USE DYNAMICS, ONLY: mu=>pua, mv=>pva, mw=>sda, mb
      USE QUSCOM, ONLY : IM,JM,LM
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: DT
      LOGICAL :: mneg
      INTEGER :: i,j,l,n,nc,im1
      double precision :: byn,ssp,snp

      mu(:,:,:) = mu(:,:,:)*(.5*dt)
      mv(:,1:jm-1,:) = mv(:,2:jm,:)*dt
      mv(:,jm,:) = 0.
      mw(:,:,1:lm-1) = mw(:,:,1:lm-1)*(-dt)
c     for some reason mu is not zero at the poles...
      mu(:,1,:) = 0.
      mu(:,jm,:) = 0.
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
              write(6,900)'ma(i,j,l) x1: nc=',i,j,l,ma(i,j,l)/mb(i,j,l)
     *             ,nc
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
            write(6,900)'ma(i,j,l) y: nc=',i,j,l,ma(i,j,l)/mb(i,j,l),nc
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
            write(6,910)'ma(1,1,l) ysp: nc=',l,ma(1,1,l)/mb(1,1,l),nc
            mneg=.true.
            go to 595
          endif
          if (ma(1,jm,l)/mb(1,jm,l).lt.0.5) then
            write(6,910)'ma(1,jm,l) ynp: nc=',l,ma(1,jm,l)/mb(1,jm,l),nc
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
            write(6,900)'ma(i,j,l) z : nc=',i,j,l,ma(i,j,l)/mb(i,j,l),nc
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
            write(6,900)'ma(i,j,l) z1 : nc=',i,j,l,ma(i,j,l)/mb(i,j,l)
     *           ,nc
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
            write(6,900)'ma(i,j,l) zlm: nc=',i,j,l,ma(i,j,l)/mb(i,j,l)
     *           ,nc
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
              write(6,900)'ma(i,j,l) x2: nc=',i,j,l,ma(i,j,l)/mb(i,j,l)
     *             ,nc
              mneg=.true.
              go to 595
            endif
            im1 = i
          end do
        end do
        end do
      end do
 595  continue
      if(ncyc.ge.10) stop 'ncyc=10 in AADVQ0'
 600  enddo
      if(ncyc.gt.1) write(6,*) 'AADVQ0: ncyc>1',ncyc
C**** Divide the mass fluxes by the number of cycles
      byn = 1./ncyc
      mu(:,:,:)=mu(:,:,:)*byn
      mv(:,1:jm-1,:)=mv(:,1:jm-1,:)*byn
      mw(:,:,:)=mw(:,:,:)*byn
C****
C**** Decide how many timesteps to take by computing Courant limits
C****
      MA(:,:,:) = MB(:,:,:)
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
      use QUSCOM, only : im,jm,lm, xstride,am,f_i,fmom_i
      use QUSDEF
      implicit none
      double precision, dimension(im,jm,lm) :: rm,mass,mu
      double precision, dimension(nmom,im,jm,lm) :: rmom
      logical ::  qlimit
      character*8 tname
      integer :: i,j,l,ierr,nerr,ns,nstep(jm,lm)

c**** loop over layers and latitudes
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
        if (ierr.eq.2) stop "Error in qlimit: abs(a) > 1"
      end if
      enddo ! ns
      enddo ! j
      enddo ! l
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
      use QUSCOM, only : im,jm,lm, ystride,bm,f_j,fmom_j, byim
      use QUSDEF
      implicit none
      double precision, dimension(im,jm,lm) :: rm,mass,mv
      double precision, dimension(nmom,im,jm,lm) :: rmom
      logical ::  qlimit
      double precision, dimension(im,jm) :: fqv
      double precision, intent(out), dimension(jm,lm) :: sfbm,sbm,sbf
      character*8 tname
      integer :: i,j,l,ierr,nerr,ns,nstep(lm)
      double precision ::
     &     m_sp,m_np,rm_sp,rm_np,rzm_sp,rzm_np,rzzm_sp,rzzm_np

c**** loop over layers
      do l=1,lm
      fqv(:,:) = 0.
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
c**** loop over timesteps
      do ns=1,nstep(l)
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
        if (ierr.eq.2) stop "Error in qlimit: abs(b) > 1"
      end if
      fqv(i,:) = fqv(i,:) + f_j(:)  !store tracer flux in fqv array
      fqv(i,jm) = 0.   ! play it safe
      rmom(ihmoms,i,1 ,l) = 0.! horizontal moments are zero at pole
      rmom(ihmoms,i,jm,l) = 0.
c     sbfijl(i,:,l) = sbfijl(i,:,l)+f_j(:)
      enddo  ! end loop over longitudes
      enddo  ! end loop over timesteps
      do j=1,jm-1   !diagnostics
        sfbm(j,l) = sfbm(j,l) + sum(fqv(:,j)/mv(:,j,l))
        sbm (j,l) = sbm (j,l) + sum(mv(:,j,l))
        sbf (j,l) = sbf (j,l) + sum(fqv(:,j))
      enddo

c**** average and unscale polar boxes
      mass(:,1 ,l) = (m_sp + sum(mass(:,1 ,l)-m_sp))*byim
      mass(:,jm,l) = (m_np + sum(mass(:,jm,l)-m_np))*byim
      rm(:,1 ,l) = (rm_sp + sum(rm(:,1 ,l)-rm_sp))*byim
      rm(:,jm,l) = (rm_np + sum(rm(:,jm,l)-rm_np))*byim
      rmom(mz ,:,1 ,l) = (rzm_sp  + sum(rmom(mz ,:,1 ,l)-rzm_sp ))*byim
      rmom(mzz,:,1 ,l) = (rzzm_sp + sum(rmom(mzz,:,1 ,l)-rzzm_sp))*byim
      rmom(mz ,:,jm,l) = (rzm_np  + sum(rmom(mz ,:,jm,l)-rzm_np ))*byim
      rmom(mzz,:,jm,l) = (rzzm_np + sum(rmom(mzz,:,jm,l)-rzzm_np))*byim
      enddo  ! end loop over levels
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
      use QUSCOM, only : im,jm,lm, zstride,cm,f_l,fmom_l
      use QUSDEF
      implicit none
      double precision, dimension(im,jm,lm) :: rm,mass,mw
      double precision, dimension(nmom,im,jm,lm) :: rmom
      integer, dimension(im,jm) :: nstep
      logical ::  qlimit
      double precision, dimension(lm) :: fqw
      double precision, intent(out), dimension(jm,lm) :: sfcm,scm,scf
      character*8 tname
      integer :: i,j,l,ierr,nerr,ns

c**** loop over latitudes and longitudes
      do j=1,jm
      do i=1,im
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
        if (ierr.eq.2) stop "Error in qlimit: abs(c) > 1"
      end if
      fqw(:)  = fqw(:) + f_l(:) !store tracer flux in fqw array
      enddo ! ns
      do l=1,lm-1   !diagnostics
        sfcm(j,l) = sfcm(j,l) + (fqw(l)/mw(i,j,l))
        scm (j,l) = scm (j,l) + mw(i,j,l)
        scf (j,l) = scf (j,l) + fqw(l)
      enddo
      enddo ! i
      enddo ! j
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
      integer :: l,j,i,ip1,im1,nstep,ns
      double precision :: courmax

C**** Decide how many timesteps to take by computing Courant limits
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
        if(nstep.ge.20) stop 'aadvqx: nstep.ge.20' !for debuging only
      enddo      ! while(courmax.gt.1.)
C**** Correct air mass
      M(:,J,L) = MI(:)
      NSTEPX(J,L) = NSTEP
c     if(nstep.gt.2 .and. nx.eq.1)
c    *  write(6,'(a,3i3,f7.4)')
c    *  'aadvqx: j,l,nstep,courmax=',j,l,nstep,courmax
  420 CONTINUE
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
      integer :: jprob,iprob,nstep,ns,i,j,l
      double precision :: courmax,byn,sbms,sbmn

C**** decide how many timesteps to take (all longitudes at this level)
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
           stop 'aadvqy: nstep.ge.20'
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
      integer :: nstep,ns,l,i
      double precision :: courmax,byn

C**** decide how many timesteps to take
      DO 430 I=IM,IM*(JM-1)+1
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
        if(nstep.ge.20) stop 'aadvqz: nstep.ge.20'
      enddo      ! while(courmax.gt.1.)
C**** Correct air mass
      m(i,1,:) = ml(:)
      NSTEPZ(I) = NSTEP
c     if(nstep.gt.1 .and. nTRACER.eq.1) write(6,'(a,2i7,f7.4)')
c    *   'aadvqz: i,nstep,courmax=',i,nstep,courmax
  430 CONTINUE
      RETURN
      END
