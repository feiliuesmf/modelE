c
c July 2008 version of the Quadratic Upstream Scheme for a
c cubed-sphere grid.  Development/testing took place within
c the GSFC FVcubed dynamical core environment.  ModelE-ification,
c removal of debugging code, other changes for production
c mode coming soon.
c
c M. Kelley 12/2008
c

c
c TODO: make a calc_ncyc which reflects our particular operator-splitting sequence.
c
      module tqus_com
c      use constant, only : grav ! mfx,mfy have units of delp, so grav not needed
      use fv_mp_mod, only : mp_reduce_max,tile
      use fv_control_mod, only : npx,npy,npz
      use fv_grid_utils_mod, only: sina_u, sina_v
      use fv_grid_tools_mod,  only: area
      use tqus_mp_mod, only : gid,domain,domain_decomp,
     &     is,ie,js,je, isd, ied, jsd, jed
      use mpp_domains_mod
      implicit none
      save
      private
      integer, parameter :: nmom=9,
     &     mx=1,my=2,mz=3, mxx=4,myy=5,mzz=6, mxy=7,myz=8,mzx=9
      real*8 :: xtoy_w,xtoy_e,xtoy_s,xtoy_n
      real*8, dimension(:,:,:,:,:), allocatable :: trmom
      real*8, dimension(:,:,:), allocatable :: ma,mb,mu,mv,mw
      real*8, dimension(:,:,:), allocatable :: rm ! temporary for testing
      integer :: nxg,nyg,ntm,ncyc
      integer :: itr
      public :: tqus_init,trmom,aadvq0,aadvq
      contains

      subroutine tqus_init(nq)
      implicit none
      integer, intent(in) :: nq
      call domain_decomp
c npx,npy are the number of A-grid cells + 1.
      nxg = npx-1
      nyg = npy-1
      ntm = nq
      allocate(trmom(nmom,isd:ied,jsd:jed,npz,ntm))
      trmom(:,:,:,:,:) = -1d30
      allocate(
     &     ma(isd:ied,jsd:jed,npz)
     &    ,mb(isd:ied,jsd:jed,npz)
     &    ,mu(isd:ied+1,jsd:jed,npz)
     &    ,mv(isd:ied,jsd:jed+1,npz)
     &    ,mw(isd:ied,jsd:jed,npz)
     &    ,rm(isd:ied,jsd:jed,npz) ! temporary for testing
     &     )
      ma(:,:,:) = -1d30
      mb(:,:,:) = -1d30
      mu(:,:,:) = -1d30
      mv(:,:,:) = -1d30
      mw(:,:,:) = -1d30
      rm(:,:,:) = -1d30

      xtoy_n = 0
      if(je==nyg.and.mod(tile,2).eq.1) xtoy_n = -1
      xtoy_s = 0
      if(js==1  .and.mod(tile,2).eq.0) xtoy_s = -1
      xtoy_w = 0
      if(is==1  .and.mod(tile,2).eq.1) xtoy_w = +1
      xtoy_e = 0
      if(ie==nxg.and.mod(tile,2).eq.0) xtoy_e = +1

      return
      end subroutine tqus_init

ctst      subroutine calc_ncyc(q_split,cx,cy)
ctst! copied from fv_tracer2d.  Assume for now that cx,cy correspond to
ctst! "mass" courant numbers, and that vertical courant numbers will
ctst! never be limiting.
ctst      implicit none
ctst! inputs:
ctst      integer :: q_split
ctst      real*8, intent(INOUT) ::  cx(is:ie+1,jsd:jed  ,npz) ! Courant Number X-Dir
ctst      real*8, intent(INOUT) ::  cy(isd:ied,js :je +1,npz) ! Courant Number Y-Dir
ctst! local vars
ctst      real*8 :: cmax(npz)
ctst      real*8 :: c_global
ctst      integer :: i,j,k
ctst!
ctst      if ( q_split == 0 ) then
ctst        do k=1,npz
ctst          cmax(k) = 0.
ctst          do j=js,je
ctst            do i=is,ie
ctst              cmax(k) = max(abs(cx(i,j,k))+1.-sina_u(i,j),
ctst     &             abs(cy(i,j,k))+1.-sina_v(i,j), cmax(k))
ctst            enddo
ctst          enddo
ctst        enddo
ctst        call mp_reduce_max(cmax,npz)
ctst! find global max courant number and define nsplt to scale cx,cy,mfx,mfy
ctst        c_global = maxval(cmax)
ctst        ncyc = int(1. + c_global)
ctst        if ( gid == 0 .and. ncyc > 5 )  write(6,*) 'Tracer_2d_split=',
ctst     &       ncyc, c_global
ctst      else
ctst        ncyc = q_split
ctst      endif
ctst
ctst      return
ctst      end subroutine calc_ncyc

      subroutine aadvq0(mfx,mfy,delp1,delp2)
c this version is for testing within fvcubed.  mfx and mfy have
c already been rescaled by tracer2d, so no division
c of mass fluxes by the number of cycles
      real*8, dimension(is:ie+1,js:je,npz) :: mfx
      real*8, dimension(is:ie,js:je+1,npz) :: mfy
      real*8, dimension(is:ie,js:je,npz) :: delp1,delp2
      integer :: i,j,l
      real*8 :: byn

      itr = 0

      do l=1,npz
      do j=js,je
      do i=is,ie
        mb(i,j,l) = area(i,j)*delp1(i,j,npz-l+1) !*grav not needed
        ma(i,j,l) = area(i,j)*delp2(i,j,npz-l+1) !*grav not needed
      enddo
      enddo
      do j=js,je
        mu(is:ie+1,j,l) = mfx(is:ie+1,j,npz-l+1)
      enddo
      do j=js,je+1
        mv(is:ie,j,l) = mfy(is:ie,j,npz-l+1)
      enddo
      enddo ! l

      call mpp_update_domains(mb,domain)
      call mpp_update_domains(mu, mv, domain, gridtype=CGRID_SW)

      do j=js,je
      do i=is,ie
        mw(i,j,npz) = 0d0
      enddo
      enddo
      do l=npz,2,-1
      do j=js,je
      do i=is,ie
        mw(i,j,l-1) = mw(i,j,l) + ma(i,j,l) - mb(i,j,l)
     &       - (mu(i,j,l)-mu(i+1,j,l)+mv(i,j,l)-mv(i,j+1,l))
      enddo
      enddo
      enddo

      ncyc = 2

      byn = 1./ncyc
      DO L=1,NPZ
         mu(:,:,l)=mu(:,:,l)*byn
         mv(:,:,l)=mv(:,:,l)*byn
         mw(:,:,l)=mw(:,:,l)*byn
      ENDDO

      if(gid.eq.3) then
        i=ie
        j=je
        l=12
        write(6,*) 'cour max ',
     &       maxval(mu(is+1:ie+1,js:je,:)/mb(is:ie,js:je,:))
c     &      ,maxval(mv(is:ie,js+1:je+1,:)/mb(is:ie,js:je,:))
c     &      ,maxval(mw(is:ie,js:je,:)/mb(is:ie,js:je,:))
        write(6,*) 'cour check ',
     &        mu(i+1,j,l)/mb(i,j,l)
     &       ,mv(i,j+1,l)/mb(i,j,l)
      endif

      return
      end subroutine aadvq0

      subroutine do_edges_and_corners(ma,rm,rmom,mu,mv)
      real*8, dimension(isd:ied,jsd:jed) :: ma,rm
      real*8, dimension(isd:ied+1,jsd:jed) :: mu ! check these bounds
      real*8, dimension(isd:ied,jsd:jed+1) :: mv ! check these bounds
      real*8, dimension(nmom,isd:ied,jsd:jed) :: rmom
      integer :: i,j,iup,iam,jup,jbm
      real*8 :: xsigni,xsigno,ysigni,ysigno,am,bm,mnewx,mnewy

      call rotate_edges(rmom)

c
c east-west edges (excluding corners)
c
      do i=1,nxg,nxg-1
        if(i.lt.is .or. i.gt.ie) cycle
        if(i.eq.1) then
          iup=is-1
          iam=is
          xsigno = -1.
        else
          iup=ie+1
          iam=ie+1
          xsigno = +1.
        endif
        xsigni = -xsigno
        do j=max(2,jsd),min(jed,nyg-1)
          am=mu(iam,j)*xsigno
          if(am.ge.0.) then
            call qusout(am,ma(i,j),rm(i,j),
     &           rmom(mx,i,j),rmom(mxx,i,j),xsigno)
            call lusout(am,ma(i,j),rmom(my,i,j),rmom(mxy,i,j),xsigno)
            call lusout(am,ma(i,j),rmom(mz,i,j),rmom(mzx,i,j),xsigno)
            rmom((/myy,myz,mzz/),i,j) =
     &           rmom((/myy,myz,mzz/),i,j)*(1.-am/ma(i,j))
          else
            am = -am
            call qusin(am,ma(i,j),rm(i,j),rmom(mx,i,j),rmom(mxx,i,j),
     &        ma(iup,j),rm(iup,j),rmom(mx,iup,j),rmom(mxx,iup,j),xsigni)
            call lusin(am,ma(i,j),rmom(my,i,j),rmom(mxy,i,j),
     &           ma(iup,j),rmom(my,iup,j),rmom(mxy,iup,j),xsigni)
            call lusin(am,ma(i,j),rmom(mz,i,j),rmom(mzx,i,j),
     &           ma(iup,j),rmom(mz,iup,j),rmom(mzx,iup,j),xsigni)
            call susin(am,rmom(myy,iup,j),ma(iup,j),rmom(myy,iup,j))
            call susin(am,rmom(myz,iup,j),ma(iup,j),rmom(myz,iup,j))
            call susin(am,rmom(mzz,iup,j),ma(iup,j),rmom(mzz,iup,j))
          endif
          ma(i,j) = ma(i,j) + mu(iam,j)*xsigni
        enddo
      enddo

c
c north-south edges (excluding corners)
c
      do j=1,nyg,nyg-1
        if(j.lt.js .or. j.gt.je) cycle
        if(j.eq.1) then
          jup=js-1
          jbm=js
          ysigno = -1.
        else
          jup=je+1
          jbm=je+1
          ysigno = +1.
        endif
        ysigni = -ysigno
        do i=max(2,isd),min(ied,nxg-1)
          bm=mv(i,jbm)*ysigno
          if(bm.ge.0.) then
            call qusout(bm,ma(i,j),rm(i,j),
     &           rmom(my,i,j),rmom(myy,i,j),ysigno)
            call lusout(bm,ma(i,j),rmom(mx,i,j),rmom(mxy,i,j),ysigno)
            call lusout(bm,ma(i,j),rmom(mz,i,j),rmom(myz,i,j),ysigno)
            rmom((/mxx,mzx,mzz/),i,j) =
     &           rmom((/mxx,mzx,mzz/),i,j)*(1.-bm/ma(i,j))
          else
            bm = -bm
            call qusin(bm,ma(i,j),rm(i,j),rmom(my,i,j),rmom(myy,i,j),
     &        ma(i,jup),rm(i,jup),rmom(my,i,jup),rmom(myy,i,jup),ysigni)
            call lusin(bm,ma(i,j),rmom(mx,i,j),rmom(mxy,i,j),
     &           ma(i,jup),rmom(mx,i,jup),rmom(mxy,i,jup),ysigni)
            call lusin(bm,ma(i,j),rmom(mz,i,j),rmom(myz,i,j),
     &           ma(i,jup),rmom(mz,i,jup),rmom(myz,i,jup),ysigni)
            call susin(bm,rmom(mxx,i,jup),ma(i,jup),rmom(mxx,i,jup))
            call susin(bm,rmom(mzx,i,jup),ma(i,jup),rmom(mzx,i,jup))
            call susin(bm,rmom(mzz,i,jup),ma(i,jup),rmom(mzz,i,jup))
          endif
          ma(i,j) = ma(i,j) + mv(i,jbm)*ysigni
        enddo
      enddo

c
c corners:
c
      do j=1,nyg,nyg-1
        if(j.lt.js .or. j.gt.je) cycle
        if(j.eq.1) then
          jup=js-1
          jbm=js
          ysigno = -1.
        else
          jup=je+1
          jbm=je+1
          ysigno = +1.
        endif
        ysigni = -ysigno
        do i=1,nxg,nxg-1
          if(i.lt.is .or. i.gt.ie) cycle
          if(i.eq.1) then
            iup=is-1
            iam=is
            xsigno = -1.
          else
            iup=ie+1
            iam=ie+1
            xsigno = +1.
          endif
          xsigni = -xsigno
          am=mu(iam,j)*xsigno
          bm=mv(i,jbm)*ysigno
          rmom((/mxx,myy,mxy/),i,j) = 0d0
          if(am.ge.0. .and. bm.ge.0.) then ! flow exiting both x and y
            call lusout_2sides(am,bm,ma(i,j),rm(i,j),
     &           rmom(mx,i,j),rmom(my,i,j),xsigno,ysigno)
            call lusout_2sides(am,bm,ma(i,j),rmom(mz,i,j),
     &           rmom(mzx,i,j),rmom(myz,i,j),xsigno,ysigno)
            rmom(mzz,i,j) = rmom(mzz,i,j)*(1.-(am+bm)/ma(i,j))
          elseif(am.lt.0. .and. bm.lt.0.) then ! flow entering both x and y
            am = -am
            bm = -bm
            call lusin_2sides(am,bm
     &           ,ma(i,j),rm(i,j),rmom(mx,i,j),rmom(my,i,j)
     &           ,ma(iup,j),rm(iup,j),rmom(mx,iup,j),rmom(my,iup,j)
     &           ,ma(i,jup),rm(i,jup),rmom(mx,i,jup),rmom(my,i,jup)
     &           ,xsigni,ysigni)
            call lusin_2sides(am,bm
     &           ,ma(i,j),rmom(mz,i,j),rmom(mzx,i,j),rmom(myz,i,j)
     &         ,ma(iup,j),rmom(mz,iup,j),rmom(mzx,iup,j),rmom(myz,iup,j)
     &         ,ma(i,jup),rmom(mz,i,jup),rmom(mzx,i,jup),rmom(myz,i,jup)
     &           ,xsigni,ysigni)
            call susin_2sides(am,bm
     &           ,ma(i,j),rmom(mzz,i,j)
     &           ,ma(iup,j),rmom(mzz,iup,j)
     &           ,ma(i,jup),rmom(mzz,i,jup)
     &           )
          elseif(am.le.0.) then ! flow exits in y, then enters in x
            call lusout(bm,ma(i,j),rm(i,j),rmom(my,i,j),ysigno)
            call lusout(bm,ma(i,j),rmom(mz,i,j),rmom(myz,i,j),ysigno)
            rmom((/mx,mzx,mzz/),i,j) =
     &           rmom((/mx,mzx,mzz/),i,j)*(1.-bm/ma(i,j))
            mnewy = ma(i,j)-bm
            am = -am
            call lusin(am,mnewy,rm(i,j),rmom(mx,i,j),
     &           ma(iup,j),rm(iup,j),rmom(mx,iup,j),xsigni)
            call lusin(am,mnewy,rmom(mz,i,j),rmom(mzx,i,j),
     &           ma(iup,j),rmom(mz,iup,j),rmom(mzx,iup,j),xsigni)
            call susin(am,rmom(my,i,j),ma(iup,j),rmom(my,iup,j))
            call susin(am,rmom(myz,i,j),ma(iup,j),rmom(myz,iup,j))
            call susin(am,rmom(mzz,i,j),ma(iup,j),rmom(mzz,iup,j))
          elseif(bm.le.0.) then ! flow exits in x, then enters in y
            call lusout(am,ma(i,j),rm(i,j),rmom(mx,i,j),xsigno)
            call lusout(am,ma(i,j),rmom(mz,i,j),rmom(mzx,i,j),xsigno)
            rmom((/my,myz,mzz/),i,j) =
     &           rmom((/my,myz,mzz/),i,j)*(1.-am/ma(i,j))
            mnewx = ma(i,j)-am
            bm = -bm
            call lusin(bm,mnewx,rm(i,j),rmom(my,i,j),
     &           ma(i,jup),rm(i,jup),rmom(my,i,jup),ysigni)
            call lusin(bm,mnewx,rmom(mz,i,j),rmom(myz,i,j),
     &           ma(i,jup),rmom(mz,i,jup),rmom(myz,i,jup),ysigni)
            call susin(bm,rmom(mx,i,j),ma(i,jup),rmom(mx,i,jup))
            call susin(bm,rmom(mzx,i,j),ma(i,jup),rmom(mzx,i,jup))
            call susin(bm,rmom(mzz,i,j),ma(i,jup),rmom(mzz,i,jup))
          endif
          ma(i,j) = ma(i,j) + mu(iam,j)*xsigni + mv(i,jbm)*ysigni
        enddo
      enddo

      return
      end subroutine do_edges_and_corners

      subroutine check_edges_and_corners(ma,rm,rmom,mu,mv)
c HALO ROWS NOT DONE HERE
      real*8, dimension(isd:ied,jsd:jed) :: ma,rm
      real*8, dimension(isd:ied+1,jsd:jed) :: mu ! check these bounds
      real*8, dimension(isd:ied,jsd:jed+1) :: mv ! check these bounds
      real*8, dimension(nmom,isd:ied,jsd:jed) :: rmom
      integer :: i,j,iam,jbm
      real*8 :: xsigno,ysigno,am,bm
c
c east-west edges (excluding corners)
c
      do i=1,nxg,nxg-1
        if(i.lt.is .or. i.gt.ie) cycle
        if(i.eq.1) then
          iam=is
          xsigno = -1.
        else
          iam=ie+1
          xsigno = +1.
        endif
        do j=max(2,js),min(je,nyg-1)
          am=mu(iam,j)*xsigno
          if(am.gt.0.) then
            call check_qusout(am,ma(i,j),rm(i,j),
     &           rmom(mx,i,j),rmom(mxx,i,j),xsigno)
          endif
        enddo
      enddo

c
c north-south edges (excluding corners)
c
      do j=1,nyg,nyg-1
        if(j.lt.js .or. j.gt.je) cycle
        if(j.eq.1) then
          jbm=js
          ysigno = -1.
        else
          jbm=je+1
          ysigno = +1.
        endif
        do i=max(2,is),min(ie,nxg-1)
          bm=mv(i,jbm)*ysigno
          if(bm.gt.0.) then
            call check_qusout(bm,ma(i,j),rm(i,j),
     &           rmom(my,i,j),rmom(myy,i,j),ysigno)
          endif
        enddo
      enddo

c
c corners:
c
      do j=1,nyg,nyg-1
        if(j.lt.js .or. j.gt.je) cycle
        if(j.eq.1) then
          jbm=js
          ysigno = -1.
        else
          jbm=je+1
          ysigno = +1.
        endif
        do i=1,nxg,nxg-1
          if(i.lt.is .or. i.gt.ie) cycle
          if(i.eq.1) then
            iam=is
            xsigno = -1.
          else
            iam=ie+1
            xsigno = +1.
          endif
          am=mu(iam,j)*xsigno
          bm=mv(i,jbm)*ysigno
          if(am.gt.0. .and. bm.gt.0.) then ! flow exiting both x and y
            call check_lusout_2sides(am,bm,ma(i,j),rm(i,j),
     &           rmom(mx,i,j),rmom(my,i,j),xsigno,ysigno)
          elseif(am.lt.0. .and. bm.lt.0.) then ! flow entering both x and y
          elseif(am.le.0.) then ! flow exits in y, then enters in x
            call check_lusout(bm,ma(i,j),rm(i,j),rmom(my,i,j),ysigno)
          elseif(bm.le.0.) then ! flow exits in x, then enters in y
            call check_lusout(am,ma(i,j),rm(i,j),rmom(mx,i,j),xsigno)
          endif
        enddo
      enddo
      return
      end subroutine check_edges_and_corners

      subroutine rotate_edges(rmom)
c xtoy =  0 means no change of orientation
c xtoy = +1 or -1 means x' = xtoy*y, y' = -xtoy*x
c xtoy for a direction defines the transformation for data coming from
c       that direction, so this routine is called after a halo_update.
      real*8, dimension(nmom,isd:ied,jsd:jed) :: rmom
      integer :: i,j
      real*8 :: xtoy,tmp
c
c east-west edges (including corners)
c
      do i=0,nxg+1,nxg+1
        if(i.lt.isd .or. i.gt.ied) cycle
        if(i.eq.0) then
          xtoy = xtoy_w
        else
          xtoy = xtoy_e
        endif
        if(xtoy.eq.0) cycle
        do j=max(1,jsd),min(jed,nyg)
          call rotxy(rmom(mx,i,j),rmom(my,i,j),xtoy)
          call rotxy(rmom(mzx,i,j),rmom(myz,i,j),xtoy)
          rmom(mxy,i,j) = -rmom(mxy,i,j)
          call swapxy(rmom(mxx,i,j),rmom(myy,i,j))
        enddo
      enddo

c
c north-south edges (including corners)
c
      do j=0,nyg+1,nyg+1
        if(j.lt.jsd .or. j.gt.jed) cycle
        if(j.eq.0) then
          xtoy = xtoy_s
        else
          xtoy = xtoy_n
        endif
        if(xtoy.eq.0) cycle
        do i=max(1,isd),min(ied,nxg)
          call rotxy(rmom(mx,i,j),rmom(my,i,j),xtoy)
          call rotxy(rmom(mzx,i,j),rmom(myz,i,j),xtoy)
          rmom(mxy,i,j) = -rmom(mxy,i,j)
          call swapxy(rmom(mxx,i,j),rmom(myy,i,j))
        enddo
      enddo
      return
      end subroutine rotate_edges

      subroutine rotxy(rx,ry,xtoy)
      real*8 :: rx,ry,xtoy
      real*8 :: tmp
      tmp = ry
      ry = -xtoy*rx
      rx = +xtoy*tmp
      return
      end subroutine rotxy

      subroutine swapxy(rx,ry)
      real*8 :: rx,ry
      real*8 :: tmp
      tmp = ry
      ry = rx
      rx = tmp
      return
      end subroutine swapxy

      subroutine checkfobs_x(ma,mu,rm,rmom)
      real*8, dimension(isd:ied,jsd:jed) :: ma,rm
      real*8, dimension(nmom,isd:ied,jsd:jed) :: rmom
      real*8, dimension(isd:ied+1,jsd:jed) :: mu
      integer :: i,j
      logical :: changed_mom
      do j=max(1,jsd),min(jed,nyg)
      do i=max(1,isd),min(ied,nxg)
        if(mu(i+1,j).gt.0. .and. mu(i,j).lt.0.) then
          call checkfobs(mu(i,j),mu(i+1,j),ma(i,j),
     &         rm(i,j),rmom(mx,i,j),rmom(mxx,i,j),
     &         changed_mom)
        endif
      enddo
      enddo
      return
      end subroutine checkfobs_x

      subroutine checkfobs_y(ma,mv,rm,rmom)
      real*8, dimension(isd:ied,jsd:jed) :: ma,rm
      real*8, dimension(nmom,isd:ied,jsd:jed) :: rmom
      real*8, dimension(isd:ied,jsd:jed+1) :: mv
      integer :: i,j
      logical :: changed_mom
      do j=max(1,jsd),min(jed,nyg)
      do i=is,ie
        if(mv(i,j+1).gt.0. .and. mv(i,j).lt.0.) then
          call checkfobs(mv(i,j),mv(i,j+1),ma(i,j),
     &         rm(i,j),rmom(my,i,j),rmom(myy,i,j),
     &         changed_mom)
        endif
      enddo
      enddo
      return
      end subroutine checkfobs_y

      subroutine checkfobs_z(ma,mw,rm,rmom)
      real*8, dimension(isd:ied,jsd:jed) :: ma,rm
      real*8, dimension(nmom,isd:ied,jsd:jed) :: rmom
      real*8, dimension(isd:ied,jsd:jed,2) :: mw
      integer :: i,j
      logical :: changed_mom
      do j=js,je
      do i=is,ie
        if(mw(i,j,2).gt.0. .and. mw(i,j,1).lt.0.) then
          call checkfobs(mw(i,j,1),mw(i,j,2),ma(i,j),
     &         rm(i,j),rmom(mz,i,j),rmom(mzz,i,j),
     &         changed_mom)
        endif
      enddo
      enddo
      return
      end subroutine checkfobs_z

      subroutine checkfobs(aml,amr,m,rm,rxm,rxxm,changed)
      implicit none
      real*8 :: aml,amr,m,rm,rxm,rxxm
      logical :: changed
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
        changed=.true.
      else
        changed=.false.
      endif
      return
      end subroutine checkfobs

      SUBROUTINE AADVQ(RM_from_fv,RMOM)!,tname
      IMPLICIT NONE
c      character*8 tname          !tracer name
      REAL*8, dimension(is:ie,js:je,npz) :: rm_from_fv
      REAL*8, dimension(nmom,isd:ied,jsd:jed,npz) :: rmom
c locals
      REAL*8, dimension(isd:ied,jsd:jed) :: mflx,mwdn,fdn,fdn0
      REAL*8, dimension(nmom,isd:ied,jsd:jed) :: fmomdn
      real*8, dimension(jsd:jed) :: mu_west,mu_east
      real*8, dimension(isd:ied) :: mv_south,mv_north

      INTEGER :: I,J,L,nc

      itr = itr + 1

C****
C**** Load mass after advection from mass before advection
C****
      DO L=1,NPZ
         MA(:,:,L) = MB(:,:,L) ! fill in halo lats
      ENDDO

C****
C**** convert from concentration to mass units and flip vertical ordering
C****
      do l=1,npz
      do j=js,je
      do i=is,ie
        RM(I,J,L)=RM_from_fv(I,J,NPZ-L+1)*MA(I,J,L)
      enddo
      enddo
      enddo

      if(rmom(1,isd,jsd,1).eq.-1d30) call init_rmom(rm,rmom,ma)

      do nc=1,ncyc

c horizontal flux limits at edges
      do l=1,npz
        call check_edges_and_corners(ma(isd,jsd,l),rm(isd,jsd,l),
     &       rmom(1,isd,jsd,l),mu(isd,jsd,l),mv(isd,jsd,l))
      enddo

c 3D halo updates
      call mpp_update_domains(ma,domain)
      call mpp_update_domains(rm,domain)
      call mpp_update_domains_column(rmom,domain)

      do j=js,je
      do i=is,ie
        mwdn(i,j) = 0.
        fdn(i,j) = 0.
        fdn0(i,j) = 0.
        fmomdn(:,i,j) = 0.
      enddo
      enddo

      DO L=1,NPZ+1 ! npz+1 is for the uppermost aadvqz call

        if(l.le.npz) then

c horz transport across the edges of faces
          call do_edges_and_corners(ma(isd,jsd,l),rm(isd,jsd,l),
     &         rmom(1,isd,jsd,l),mu(isd,jsd,l),mv(isd,jsd,l))

c testing
c          rm(is:ie,js:je,l) =
c     *         RM_from_fv(is:ie,js:je,NPZ-L+1)*MB(is:ie,js:je,l)
c          call zero_corners_and_edges(rmom(1,isd,jsd,l))
c          MA(is:ie,js:je,l) = MB(is:ie,js:je,l)
c          cycle

c horz transport in the interior of faces.
c mu,mv temporarily set to zero at edges.
          if(is.eq.1) then
            do j=jsd,jed
              mu_west(j) = mu(is  ,j,l)
              mu(is  ,j,l) = 0d0
            enddo
          endif
          if(ie.eq.nxg) then
            do j=jsd,jed
              mu_east(j) = mu(ie+1,j,l)
              mu(ie+1,j,l) = 0d0
            enddo
          endif
          CALL AADVQX(RM(isd,jsd,l),RMOM(1,isd,jsd,l),
     &         MA(isd,jsd,l),mu(isd,jsd,l))
          if(is.eq.1) then
            do j=jsd,jed
              mu(is  ,j,l) = mu_west(j)
            enddo
          endif
          if(ie.eq.nxg) then
            do j=jsd,jed
              mu(ie+1,j,l) = mu_east(j)
            enddo
          endif

c testing
c          rm(is:ie,js:je,l) =
c     *         RM_from_fv(is:ie,js:je,NPZ-L+1)*MB(is:ie,js:je,l)
c          call zero_corners_and_edges(rmom(1,isd,jsd,l))
c          ma(:,:,l) = mb(:,:,l)
c          MA(is:ie,js:je,l) = MB(is:ie,js:je,l)
c          cycle

          if(js.eq.1) then
            do i=is,ie
              mv_south(i) = mv(i,js  ,l)
              mv(i,js  ,l) = 0d0
            enddo
          endif
          if(je.eq.nyg) then
            do i=is,ie
              mv_north(i) = mv(i,je+1,l)
              mv(i,je+1,l) = 0d0
            enddo
          endif
          CALL AADVQY(RM(isd,jsd,l),RMOM(1,isd,jsd,l),
     &         MA(isd,jsd,l),mv(isd,jsd,l))
          if(js.eq.1) then
            do i=is,ie
              mv(i,js  ,l) = mv_south(i)
            enddo
          endif
          if(je.eq.nyg) then
            do i=is,ie
              mv(i,je+1,l) = mv_north(i)
            enddo
          endif

c testing
c          rm(is:ie,js:je,l) =
c     *         RM_from_fv(is:ie,js:je,NPZ-L+1)*MB(is:ie,js:je,l)
c          call zero_corners_and_edges(rmom(1,isd,jsd,l))
c          ma(:,:,l) = mb(:,:,l)
c          MA(is:ie,js:je,l) = MB(is:ie,js:je,l)
c          cycle

        endif

c vertical transport
        if(l.gt.1 .and. l.lt.npz) then
          call checkfobs_z(mw(isd,jsd,l-1),
     &         ma(isd,jsd,l),rm(isd,jsd,l),rmom(1,isd,jsd,l))
        endif
        if(l.gt.1) then
          CALL AADVQZ(RM(isd,jsd,l-1),RMOM(1,isd,jsd,l-1),
     &         MA(isd,jsd,l-1),mw(isd,jsd,l-1),
     &         mwdn,fdn,fmomdn,fdn0)
        endif

c testing
        if(l.gt.1) then
c          rm(is:ie,js:je,l-1) =
c     *         RM_from_fv(is:ie,js:je,NPZ-(L-1)+1)*MB(is:ie,js:je,l-1)
c          call zero_corners_and_edges(rmom(1,isd,jsd,l-1))
c          MA(is:ie,js:je,l-1) = MB(is:ie,js:je,l-1)
        endif

      ENDDO ! end loop over l

      enddo ! nc=1,ncyc

C****
C**** convert from mass to concentration units
C****
c      if(.false.) then ! testing
      do l=1,npz
      do j=js,je
      do i=is,ie
        RM_from_fv(I,J,NPZ-L+1)=RM(I,J,L)/MA(I,J,L)
      enddo
      enddo
      enddo
c      endif
c      rmom = 0d0 ! testing

c      if(itr.eq.5 .and. gid.eq.3) then
cc        rm(is:ie,js:je,:) = mb(is:ie,js:je,:)+ncyc*(
cc     &       mu(is:ie,js:je,:)-mu(is+1:ie+1,js:je,:)
cc     &       +mv(is:ie,js:je,:)-mv(is:ie,js+1:je+1,:)
cc     &       )
cc        rm(is:ie,js:je,2:npz) = rm(is:ie,js:je,2:npz)
cc     &       +ncyc*mw(is:ie,js:je,1:npz-1)
cc        rm(is:ie,js:je,1:npz) = rm(is:ie,js:je,1:npz)
cc     &       -ncyc*mw(is:ie,js:je,1:npz)
c        l=12
c        i=ie; j=je
c        write(6,*) 'check ne',1d0-ma(i,j,l)/rm(i,j,l)
c        i=ie; j=js
c        write(6,*) 'check se',1d0-ma(i,j,l)/rm(i,j,l)
c        i=is; j=je
c        write(6,*) 'check nw',1d0-ma(i,j,l)/rm(i,j,l)
c        i=is; j=js
c        write(6,*) 'check sw',1d0-ma(i,j,l)/rm(i,j,l)
c
c        i=(is+ie)/2; j=js
c        write(6,*) 'check so',1d0-ma(i,j,l)/rm(i,j,l)
c        i=(is+ie)/2; j=je
c        write(6,*) 'check no',1d0-ma(i,j,l)/rm(i,j,l)
c
c        i=is; j=(js+je)/2
c        write(6,*) 'check we',1d0-ma(i,j,l)/rm(i,j,l)
c        i=ie; j=(js+je)/2
c        write(6,*) 'check ea',1d0-ma(i,j,l)/rm(i,j,l)
c
c        i=(is+ie)/2; j=(js+je)/2
c        write(6,*) 'check ce',1d0-ma(i,j,l)/rm(i,j,l)
c      endif

c      if(itr.eq.5 .and. gid.eq.47) then
c        write(6,*) 'nbad ',
c     &       count(abs(1d0-rm(is:ie,js:je,1)/ma(is:ie,js:ie,1)).gt.1d-6)
c        write(6,*) 'check glb',gid,
c     &       minval(1d0-rm(is:ie,js:je,:)/ma(is:ie,js:ie,:))
c     &      ,maxval(1d0-rm(is:ie,js:je,:)/ma(is:ie,js:ie,:))
c        write(6,*) 'minloc glb',gid,
c     &       minloc(1d0-rm(is:ie,js:je,:)/ma(is:ie,js:ie,:))
c        write(6,*) 'maxloc glb',gid,
c     &       maxloc(1d0-rm(is:ie,js:je,:)/ma(is:ie,js:ie,:))
c        do j=js,je
c        do i=is,ie
c          if(abs(1d0-rm(i,j,1)/ma(i,j,1)).gt.1d-6) then
c            write(6,*) 'badpt' ,i,j,1d0-rm(i,j,1)/ma(i,j,1)
c            write(6,*) 'mu ',mu(i:i+1,j,1)
c            write(6,*) 'mv ',mv(i,j:j+1,1)
c          endif
c        enddo
c        enddo
c      endif

c      if(itr.eq.4 .and. tile.eq.2 .and. ie.eq.nxg .and. je.eq.nyg) then
c        l=1
c        i=ie
c        j=je
c        write(6,*) 'ne rm ',rm(i,j,l)/ma(i,j,l)
c        write(6,*) 'ne rx ',rmom((/mx,mzx/),i,j,l)/ma(i,j,l)
c        write(6,*) 'ne ry ',rmom((/my,myz/),i,j,l)/ma(i,j,l)
c      endif
c      if(itr.eq.4 .and. tile.eq.3 .and. ie.eq.nxg .and. js.eq.1) then
c        l=1
c        i=ie
c        j=1
c        write(6,*) 'se rm ',rm(i,j,l)/ma(i,j,l)
c        write(6,*) 'se rx ',rmom((/mx,mzx/),i,j,l)/ma(i,j,l)
c        write(6,*) 'se ry ',rmom((/my,myz/),i,j,l)/ma(i,j,l)
c      endif
c      if(itr.eq.4 .and. tile.eq.4 .and. is.eq.1 .and. js.eq.1) then
c        l=1
c        i=1
c        j=1
c        write(6,*) 'sw rm ',rm(i,j,l)/ma(i,j,l)
c        write(6,*) 'sw rx ',rmom((/mx,mzx/),i,j,l)/ma(i,j,l)
c        write(6,*) 'sw ry ',rmom((/my,myz/),i,j,l)/ma(i,j,l)
c      endif

c      l=1
c      if(itr.eq.4) then
c        if(is.eq.1 .and. js.eq.1)
c     &       write(6,*) 'sw rm ',tile,rm(is+1,js,l)/ma(is+1,js,l)
c        if(ie.eq.nxg .and. js.eq.1)
c     &       write(6,*) 'se rm ',tile,rm(ie,js+1,l)/ma(ie,js+1,l)
c        if(is.eq.1 .and. je.eq.nyg)
c     &       write(6,*) 'nw rm ',tile,rm(is+1,je,l)/ma(is+1,je,l)
c        if(ie.eq.nxg .and. je.eq.nyg)
c     &       write(6,*) 'ne rm ',tile,rm(ie,je-1,l)/ma(ie,je-1,l)
c      endif

      RETURN
      END SUBROUTINE AADVQ

      subroutine zero_corners_and_edges(rmom)
      REAL*8, dimension(nmom,isd:ied,jsd:jed) :: rmom

c      rmom(:,isd:ied,jsd:jed) = 0d0

c      if(is.eq.  1) rmom(:,is,js:je) = 0d0
c      if(ie.eq.nxg) rmom(:,ie,js:je) = 0d0
c      if(js.eq.  1) rmom(:,is:ie,js) = 0d0
c      if(je.eq.nyg) rmom(:,is:ie,je) = 0d0

c      if(is.eq.  1 .and. js.eq.  1) rmom(:,is,js) = 0d0
c      if(is.eq.  1 .and. je.eq.nyg) rmom(:,is,je) = 0d0
c      if(ie.eq.nxg .and. js.eq.  1) rmom(:,ie,js) = 0d0
c      if(ie.eq.nxg .and. je.eq.nyg) rmom(:,ie,je) = 0d0

      return
      end subroutine zero_corners_and_edges

      subroutine aadvqx(rm,rmom,mass,mu)
!@sum  AADVQX advection driver for x-direction
!@auth Maxwell Kelley
      REAL*8, dimension(isd:ied,jsd:jed) :: rm,mass
      REAL*8, dimension(isd:ied+1,jsd:jed) :: mu
      REAL*8, dimension(nmom,isd:ied,jsd:jed) :: rmom
      integer :: i,ii,j,ns
      real*8 :: am,frac1,fracm,fw,fe,dm2,mold,mnew,bymnew
      real*8 :: fe0,fe_pass,fw0,rm0,fex_pass,fexx_pass
      real*8, dimension(nmom) :: fmomw,fmome

      call checkfobs_x(mass,mu,rm,rmom)

      do j=max(1,jsd),min(jed,nyg) ! calculate halo rows as well

c-----------------------------------------------------------
      ! calculate tracer mass flux f
c-----------------------------------------------------------
c--------------------------------------------------------------------
      ! calculate tracer fluxes of slopes and curvatures
c--------------------------------------------------------------------
c-------------------------------------------------------------------
c update tracer mass, moments of tracer mass, air mass distribution
c-------------------------------------------------------------------
      i=is-1
      am = mu(i+1,j)
      if(am.le.0.) then      ! air mass flux is negative
        ii=i+1
        frac1=+1.
      else                      ! air mass flux is positive
        ii=i
        frac1=-1.
      endif
      fracm=am/mass(ii,j)
      frac1=fracm+frac1
      fe=fracm*(rm(ii,j)-frac1*(rmom(mx,ii,j)-
     &     (frac1+fracm)*rmom(mxx,ii,j)))
      fmome(mx)=am*(fracm*fracm*(rmom(mx,ii,j)
     &     -3.*frac1*rmom(mxx,ii,j))-3.*fe)
      fmome(mxx)=am*(am*fracm**3 *rmom(mxx,ii,j)
     &     -5.*(am*fe+fmome(mx)))

      ! cross moments
      fmome(my)  = fracm*(rmom(my,ii,j)-frac1*rmom(mxy,ii,j))
      fmome(mxy) = am*(fracm*fracm*rmom(mxy,ii,j)-3.*fmome(my))
      fmome(mz)  = fracm*(rmom(mz,ii,j)-frac1*rmom(mzx,ii,j))
      fmome(mzx) = am*(fracm*fracm*rmom(mzx,ii,j)-3.*fmome(mz))
      fmome(myy) = fracm*rmom(myy,ii,j)
      fmome(mzz) = fracm*rmom(mzz,ii,j)
      fmome(myz) = fracm*rmom(myz,ii,j)

! flux limitations
      fe0 = fe
      fe_pass = fe
      fex_pass = fmome(mx)
      fexx_pass = fmome(mxx)
      if(am.gt.0.) then
        if(fe.lt.0.) then
          fe=0.
          fe_pass=0.
          fex_pass=0
          fexx_pass=0.
        elseif(fe.gt.rm(i,j)) then
          fe=rm(i,j)
          fe_pass = fe
          fex_pass=am*(-3.*fe)
          fexx_pass=am*(-5.*(am*fe+fex_pass))
        endif
      else if(am.lt.0.) then
        if(fe.gt.0.) then
          fe=0.
          fe0=0.
          fmome((/mx,mxx/))=0.
        elseif(fe.lt.-rm(i+1,j)) then
          fe=-rm(i+1,j)
          fe0 = fe
          fmome(mx)=am*(-3.*fe)
          fmome(mxx)=am*(-5.*(am*fe+fmome(mx)))
        endif
      endif
      fw = fe
      fw0 = fe_pass
      fmomw(:) = fmome(:)
      fmomw(mx) = fex_pass
      fmomw(mxx) = fexx_pass

      do i=is,ie
         am = mu(i+1,j)
         if(am.lt.0.) then ! air mass flux is negative
            ii=i+1
            frac1=+1.
         else                 ! air mass flux is positive
            ii=i
            frac1=-1.
         endif
         fracm=am/mass(ii,j)
         frac1=fracm+frac1
         fe=fracm*(rm(ii,j)-frac1*(rmom(mx,ii,j)-
     &        (frac1+fracm)*rmom(mxx,ii,j)))
         fmome(mx)=am*(fracm*fracm*(rmom(mx,ii,j)
     &        -3.*frac1*rmom(mxx,ii,j))-3.*fe)
         fmome(mxx)=am*(am*fracm**3 *rmom(mxx,ii,j)
     &        -5.*(am*fe+fmome(mx)))
      ! cross moments
         fmome(my)  = fracm*(rmom(my,ii,j)-frac1*rmom(mxy,ii,j))
         fmome(mxy) =am*(fracm*fracm*rmom(mxy,ii,j)-3.*fmome(my))
         fmome(mz)  = fracm*(rmom(mz,ii,j)-frac1*rmom(mzx,ii,j))
         fmome(mzx) =am*(fracm*fracm*rmom(mzx,ii,j)-3.*fmome(mz))
         fmome(myy) = fracm*rmom(myy,ii,j)
         fmome(mzz) = fracm*rmom(mzz,ii,j)
         fmome(myz) = fracm*rmom(myz,ii,j)

! flux limitations
         fe0 = fe
         fe_pass = fe
         fex_pass = fmome(mx)
         fexx_pass = fmome(mxx)
         if(am.gt.0.) then
           if(fe.lt.0.) then
             fe=0.
             fe_pass=0.
             fex_pass=0
             fexx_pass=0.
           elseif(fe.gt.rm(i,j)) then
             fe=rm(i,j)
             fe_pass = fe
             fex_pass=am*(-3.*fe)
             fexx_pass=am*(-5.*(am*fe+fex_pass))
           endif
         elseif(am.lt.0.) then
           if(fe.gt.0.) then
             fe=0.
             fe0=0.
             fmome((/mx,mxx/))=0.
           elseif(fe.lt.-rm(i+1,j)) then
             fe=-rm(i+1,j)
             fe0 = fe
             fmome(mx)=am*(-3.*fe)
             fmome(mxx)=am*(-5.*(am*fe+fmome(mx)))
           endif
         endif

         mold=mass(i,j)
         mnew=mold+mu(i,j)-am
         bymnew = 1./mnew
         dm2=mu(i,j)+am
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

         fw = fe
         fmomw(:) = fmome(:)
         fw0 = fe_pass
         fmomw(mx) = fex_pass
         fmomw(mxx) = fexx_pass

      enddo ! i
      enddo ! j

      return
c****
      end subroutine aadvqx

      subroutine aadvqy(rm,rmom,mass,mv)
!@sum  AADVQY advection driver for y-direction
!@auth Maxwell Kelley
      implicit none
      REAL*8, dimension(isd:ied,jsd:jed) :: rm,mass
      REAL*8, dimension(isd:ied,jsd:jed+1) :: mv
      REAL*8, dimension(nmom,isd:ied,jsd:jed) :: rmom
      integer :: i,j,jj
      real*8, dimension(isd:ied) :: mvs,fs
      real*8, dimension(nmom,isd:ied) :: fmoms
      real*8, dimension(nmom) :: fmomn
      real*8 :: frac1,fracm,fn,mold,mnew,bymnew,dm2,am
      real*8 :: fn0,fn_pass,rm0,fny_pass,fnyy_pass
      real*8, dimension(isd:ied) :: fs0

      call checkfobs_y(mass,mv,rm,rmom)

c-----------------------------------------------------------
      ! calculate tracer mass flux f
      ! and fluxes of slopes and curvatures fmom
      ! update tracer mass, moments of tracer mass, air mass distribution
c--------------------------------------------------------------------
      j=js-1
      do i=is,ie
        am = mv(i,j+1)
        if(am.le.0.) then ! air mass flux is negative
          jj=j+1
          frac1=+1.
        else                    ! air mass flux is positive
          jj=j
          frac1=-1.
        endif
        fracm=am/mass(i,jj)
        frac1=fracm+frac1
        fn=fracm*(rm(i,jj)-frac1*(rmom(my,i,jj)-
     &       (frac1+fracm)*rmom(myy,i,jj)))
        fmomn(my)=am*(fracm*fracm*(rmom(my,i,jj)
     &       -3.*frac1*rmom(myy,i,jj))-3.*fn)
        fmomn(myy)=am*(am*fracm**3 *rmom(myy,i,jj)
     &       -5.*(am*fn+fmomn(my)))
        fmomn(mz)  = fracm*(rmom(mz,i,jj)-frac1*rmom(myz,i,jj))
        fmomn(myz) = am*
     &       (fracm*fracm*rmom(myz,i,jj)-3.*fmomn(mz))
        fmomn(mx)  = fracm*(rmom(mx,i,jj)-frac1*rmom(mxy,i,jj))
        fmomn(mxy) = am*
     &       (fracm*fracm*rmom(mxy,i,jj)-3.*fmomn(mx))
        fmomn(mzz) = fracm*rmom(mzz,i,jj)
        fmomn(mxx) = fracm*rmom(mxx,i,jj)
        fmomn(mzx) = fracm*rmom(mzx,i,jj)

! flux limitations
        fn0 = fn
        fn_pass = fn
        fny_pass = fmomn(my)
        fnyy_pass = fmomn(myy)
        if(am.gt.0.) then
          if(fn.lt.0.) then
            fn=0.
            fn_pass=0.
            fny_pass=0
            fnyy_pass=0.
          elseif(fn.gt.rm(i,j)) then
            fn=rm(i,j)
            fn_pass = fn
            fny_pass=am*(-3.*fn)
            fnyy_pass=am*(-5.*(am*fn+fny_pass))
          endif
        elseif(am.lt.0.) then
          if(fn.gt.0.) then
            fn=0.
            fn0=0.
            fmomn((/my,myy/))=0.
          elseif(fn.lt.-rm(i,j+1)) then
            fn=-rm(i,j+1)
            fn0 = fn
            fmomn(my)=am*(-3.*fn)
            fmomn(myy)=am*(-5.*(am*fn+fmomn(my)))
          endif
        endif

        mvs(i) = am
        fs(i) = fn
        fmoms(:,i) = fmomn(:)

        fs0(i) = fn_pass
        fmoms(my,i) = fny_pass
        fmoms(myy,i) = fnyy_pass
        
      enddo                     ! i

      do j=js,je
      do i=is,ie
         am = mv(i,j+1)
         if(am.lt.0.) then ! air mass flux is negative
            jj=j+1
            frac1=+1.
          else                   ! air mass flux is positive
            jj=j
            frac1=-1.
         endif
         fracm=am/mass(i,jj)
         frac1=fracm+frac1
         fn=fracm*(rm(i,jj)-frac1*(rmom(my,i,jj)-
     &        (frac1+fracm)*rmom(myy,i,jj)))
         fmomn(my)=am*(fracm*fracm*(rmom(my,i,jj)
     &        -3.*frac1*rmom(myy,i,jj))-3.*fn)
         fmomn(myy)=am*(am*fracm**3 *rmom(myy,i,jj)
     &        -5.*(am*fn+fmomn(my)))
         fmomn(mz)  = fracm*(rmom(mz,i,jj)-frac1*rmom(myz,i,jj))
         fmomn(myz) = am*
     &        (fracm*fracm*rmom(myz,i,jj)-3.*fmomn(mz))
         fmomn(mx)  = fracm*(rmom(mx,i,jj)-frac1*rmom(mxy,i,jj))
         fmomn(mxy) = am*
     &        (fracm*fracm*rmom(mxy,i,jj)-3.*fmomn(mx))
         fmomn(mzz) = fracm*rmom(mzz,i,jj)
         fmomn(mxx) = fracm*rmom(mxx,i,jj)
         fmomn(mzx) = fracm*rmom(mzx,i,jj)

! flux limitations
         fn0 = fn
         fn_pass = fn
         fny_pass = fmomn(my)
         fnyy_pass = fmomn(myy)
         if(am.gt.0.) then
           if(fn.lt.0.) then
             fn=0.
             fn_pass=0.
             fny_pass=0
             fnyy_pass=0.
           elseif(fn.gt.rm(i,j)) then
             fn=rm(i,j)
             fn_pass = fn
             fny_pass=am*(-3.*fn)
             fnyy_pass=am*(-5.*(am*fn+fny_pass))
           endif
         elseif(am.lt.0.) then
           if(fn.gt.0.) then
             fn=0.
             fn0=0.
             fmomn((/my,myy/))=0.
           elseif(fn.lt.-rm(i,j+1)) then
             fn=-rm(i,j+1)
             fn0 = fn
             fmomn(my)=am*(-3.*fn)
             fmomn(myy)=am*(-5.*(am*fn+fmomn(my)))
           endif
         endif

         mold=mass(i,j)
         mnew=mold+mvs(i)-am
         bymnew = 1./mnew
         dm2=mvs(i)+am
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

         mvs(i) = am
         fs(i) = fn
         fmoms(:,i) = fmomn(:)

         fs0(i) = fn_pass
         fmoms(my,i) = fny_pass
         fmoms(myy,i) = fnyy_pass

      enddo ! i
      enddo ! j

      return
      end subroutine aadvqy

      subroutine aadvqz(rm,rmom,mass,mw,mwdn,fdn,fmomdn,fdn0)
      REAL*8, dimension(isd:ied,jsd:jed,2) :: rm,mass
      REAL*8, dimension(nmom,isd:ied,jsd:jed,2) :: rmom
      REAL*8, dimension(nmom,isd:ied,jsd:jed) :: fmomdn
      REAL*8, dimension(isd:ied,jsd:jed) :: mw,mwdn,fdn,fdn0

      real*8, dimension(nmom) :: fmomup
      integer :: i,j,l,ll
      real*8 :: frac1,fracm,fup,mold,mnew,bymnew,dm2,am
      real*8 :: fup0,fup_pass,rm0,fupz_pass,fupzz_pass

c-----------------------------------------------------------
      ! calculate tracer mass flux f
      ! and fluxes of slopes and curvatures fmom
      ! update tracer mass, moments of tracer mass, air mass distribution
c--------------------------------------------------------------------
      do j=js,je
      do i=is,ie
         am = mw(i,j)
         if(am.lt.0.) then ! air mass flux is negative
            ll=2
            frac1=+1.
          else                   ! air mass flux is positive
            ll=1
            frac1=-1.
         endif
         fracm=am/mass(i,j,ll)
         frac1=fracm+frac1
         fup=fracm*(rm(i,j,ll)-frac1*(rmom(mz,i,j,ll)-
     &        (frac1+fracm)*rmom(mzz,i,j,ll)))
         fmomup(mz)=am*(fracm*fracm*(rmom(mz,i,j,ll)
     &        -3.*frac1*rmom(mzz,i,j,ll))-3.*fup)
         fmomup(mzz)=am*(am*fracm**3 *rmom(mzz,i,j,ll)
     &        -5.*(am*fup+fmomup(mz)))
         fmomup(my)  = fracm*(rmom(my,i,j,ll)-frac1*rmom(myz,i,j,ll))
         fmomup(myz) = am*
     &        (fracm*fracm*rmom(myz,i,j,ll)-3.*fmomup(my))
         fmomup(mx)  = fracm*(rmom(mx,i,j,ll)-frac1*rmom(mzx,i,j,ll))
         fmomup(mzx) = am*
     &        (fracm*fracm*rmom(mzx,i,j,ll)-3.*fmomup(mx))
         fmomup(myy) = fracm*rmom(myy,i,j,ll)
         fmomup(mxx) = fracm*rmom(mxx,i,j,ll)
         fmomup(mxy) = fracm*rmom(mxy,i,j,ll)

! flux limitations
         fup0 = fup
         fup_pass = fup
         fupz_pass = fmomup(mz)
         fupzz_pass = fmomup(mzz)
         if(am.gt.0.) then
           if(fup.lt.0.) then
             fup=0.
             fup_pass=0.
             fupz_pass=0
             fupzz_pass=0.
           elseif(fup.gt.rm(i,j,1)) then
             fup=rm(i,j,1)
             fup_pass = fup
             fupz_pass=am*(-3.*fup)
             fupzz_pass=am*(-5.*(am*fup+fupz_pass))
           endif
         elseif(am.lt.0.) then
           if(fup.gt.0.) then
             fup=0.
             fup0=0.
             fmomup((/mz,mzz/))=0.
           elseif(fup.lt.-rm(i,j,2)) then
             fup=-rm(i,j,2)
             fup0 = fup
             fmomup(mz)=am*(-3.*fup)
             fmomup(mzz)=am*(-5.*(am*fup+fmomup(mz)))
           endif
         endif

         mold=mass(i,j,1)
         mnew=mold+mwdn(i,j)-am
         bymnew = 1./mnew
         dm2=mwdn(i,j)+am
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

         mwdn(i,j) = am
         fdn(i,j) = fup
         fmomdn(:,i,j) = fmomup(:)

         fdn0(i,j) = fup_pass
         fmomdn(mz,i,j) = fupz_pass
         fmomdn(mzz,i,j) = fupzz_pass

      enddo ! i
      enddo ! j

      return
c****
      end subroutine aadvqz

      subroutine init_rmom(rm,rmom,ma)
      REAL*8, dimension(isd:ied,jsd:jed,npz) :: rm,ma
      REAL*8, dimension(nmom,isd:ied,jsd:jed,npz) :: rmom
      integer :: i,j,l
      real*8 :: r1,r2
      call mpp_update_domains(ma,domain)
      call mpp_update_domains(rm,domain)
      rmom(:,:,:,:) = 0d0
      do l=1,npz
      do j=js,je
      do i=is,ie
        r1 = rm(i-1,j,l)/ma(i-1,j,l)
        r2 = rm(i+1,j,l)/ma(i+1,j,l)
        rmom(mx,i,j,l) = .25*(r2-r1)
        r1 = rm(i,j-1,l)/ma(i,j-1,l)
        r2 = rm(i,j+1,l)/ma(i,j+1,l)
        rmom(my,i,j,l) = .25*(r2-r1)
      enddo ! i
      if(is.eq.1) then
        i=is
        r1 = rm(i,j,l)/ma(i,j,l)
        r2 = rm(i+1,j,l)/ma(i+1,j,l)
        rmom(mx,i,j,l) = .5*(r2-r1)
      endif
      if(ie.eq.nxg) then
        i=ie
        r1 = rm(i-1,j,l)/ma(i-1,j,l)
        r2 = rm(i,j,l)/ma(i,j,l)
        rmom(mx,i,j,l) = .5*(r2-r1)
      endif
      enddo ! j
      if(js.eq.1) then
        j=js
        do i=is,ie
          r1 = rm(i,j,l)/ma(i,j,l)
          r2 = rm(i,j+1,l)/ma(i,j+1,l)
          rmom(my,i,j,l) = .5*(r2-r1)
        enddo ! i
      endif
      if(je.eq.nyg) then
        j=je
        do i=is,ie
          r1 = rm(i,j-1,l)/ma(i,j-1,l)
          r2 = rm(i,j,l)/ma(i,j,l)
          rmom(my,i,j,l) = .5*(r2-r1)
        enddo ! i
      endif
      do j=js,je
      do i=is,ie
        rmom(:,i,j,l) = rmom(:,i,j,l)*ma(i,j,l)
      enddo
      enddo
      enddo ! l
      return
      end subroutine init_rmom

      subroutine qusin(am,m,rm,rxm,rxxm,m_up,rm_up,rxm_up,rxxm_up,xsign)
c update rxm,rxxm as if upwind box had zero moments and r = fm/am
      implicit none
      real*8 :: m,m_up,am,rm,rxm,rxxm,rm_up,rxm_up,rxxm_up,xsign
      real*8 :: mnew,a,fm,dr,ddr,mrat
      A    = AM / M_UP
      FM   = A*(RM_UP + (1.-A)*(XSIGN*RXM_UP + (1.-2.*A)*RXXM_UP))
      MNEW   =  M + AM
      rxm = rxm*xsign
      dr  = rm/m - fm/am
      ddr = rxm/m-dr*(1d0-am/m)
      mrat = m/mnew
      RM   =  RM   + FM
      RXM  = (RXM  + 3.*AM*DR )*MRAT
      RXXM = (RXXM + 5.*AM*DDR)*MRAT*MRAT
      rxm = rxm*xsign
      return
      end subroutine qusin

      subroutine lusin(am,m,rm,rxm,m_up,rm_up,rxm_up,xsign)
c update rxm as if upwind box had zero moments and r = fm/am
      implicit none
      real*8 :: m,m_up,am,rm,rxm,rm_up,rxm_up,xsign
      real*8 :: mnew,a,fm,amdr
      A    = AM / M_UP
      FM   = A*(RM_UP + (1.-A)*(XSIGN*RXM_UP))
      MNEW   =  M + AM
      amdr = am*rm/m - fm
      RM  = RM + FM
      RXM = (RXM + 3.*XSIGN*AMDR) * M / MNEW
      return
      end subroutine lusin

      subroutine susin(am,rm,m_up,rm_up)
      implicit none
      real*8 :: m_up,am,rm,rm_up
      RM = RM + AM*RM_UP/M_UP
      return
      end subroutine susin

      subroutine qusout(am,m,rm,rxm,rxxm,xsign)
      implicit none
      real*8 :: am,m,rm,rxm,rxxm,xsign
      real*8 :: bym,mnew,a,rl,rr,x
      real*8, parameter :: by3=1d0/3d0
      RXM = RXM*XSIGN
      bym = 1./m
      A  = AM*bym
      RL = (RM - RXM + RXXM)*bym
      x = 1.-2.*a
      RR = (RM + RXM*x + RXXM*1.5*(X*X-by3))*bym
      MNEW = M - AM
      RM = RM - A*(RM + (1.-A)*(RXM + (1.-2.*A)*RXXM))
      RXM = MNEW*.5*(RR-RL)
      RXXM = MNEW*.5*(RR+RL)-RM
      RXM = RXM*XSIGN
      return
      end subroutine qusout

      subroutine check_qusout(am,m,rm,rxm,rxxm,xsign)
      implicit none
      real*8 :: am,m,rm,rxm,rxxm,xsign
      real*8 :: a,fxtra,rm0
      A      = AM / M
      fxtra = A*(1.-A)*(RXM*XSIGN + (1.-2.*A)*RXXM)
      RM0 = RM*(1.-A)
      if(RM0 .lt. fxtra) then
        rxm = rxm*rm0/fxtra
        rxxm = rxxm*rm0/fxtra
      endif
      return
      end subroutine check_qusout

      subroutine lusout(am,m,rm,rxm,xsign)
      implicit none
      real*8 :: am,m,rm,rxm,xsign
      real*8 :: mnew,a
      A      = AM / M
      MNEW   =  M-AM
      RM  = RM - A*(RM + (1.-A)*RXM*XSIGN)
      rxm = rxm*(mnew/m)**2
      return
      end subroutine lusout

      subroutine check_lusout(am,m,rm,rxm,xsign)
      implicit none
      real*8 :: am,m,rm,rxm,rxxm,xsign
      real*8 :: a,fxtra,rm0
      A      = AM / M
      fxtra = A*(1.-A)*RXM*XSIGN
      RM0 = RM*(1.-A)
      if(RM0 .lt. fxtra) then
        rxm = rxm*rm0/fxtra
      endif
      return
      end subroutine check_lusout

      subroutine lusout_2sides(am,bm,m,rm,rxm,rym,xsign,ysign)
      implicit none
      real*8 :: am,bm,m,rm,rxm,rym,xsign,ysign
      real*8 :: mnew,mnewx,mnewy,a,b
      A      = AM / M
      B      = BM / M
      MNEWX   =  M-AM
      MNEWY   =  M-BM
      MNEW    =  M-AM-BM
      RM  = RM -A*(RM + (1.-A)*RXM*XSIGN) -B*(RM + (1.-B)*RYM*YSIGN)
      rxm = rxm*(mnew*mnewx/(m*m))
      rym = rym*(mnew*mnewy/(m*m))
      return
      end subroutine lusout_2sides

      subroutine check_lusout_2sides(am,bm,m,rm,rxm,rym,xsign,ysign)
      implicit none
      real*8 :: am,bm,m,rm,rxm,rym,xsign,ysign
      real*8 :: a,b,fxtra,rm0
      A      = AM / M
      B      = BM / M
      fxtra = +A*(1.-A)*RXM*XSIGN + B*(1.-B)*RYM*YSIGN
      RM0 = RM*(1.-A-B) 
      if(RM0 .lt. fxtra) then
        rxm = rxm*rm0/fxtra
        rym = rym*rm0/fxtra
      endif
      return
      end subroutine check_lusout_2sides

      subroutine lusin_2sides(am,bm
     &     ,m,rm,rxm,rym
     &     ,m_w,rm_w,rxm_w,rym_w
     &     ,m_s,rm_s,rxm_s,rym_s
     &     ,xsign,ysign)
      implicit none
      real*8 :: am,bm
     &     ,m,rm,rxm,rym
     &     ,m_w,rm_w,rxm_w,rym_w
     &     ,m_s,rm_s,rxm_s,rym_s
     &     ,xsign,ysign
      real*8 :: mnew,mnewx,mnewy,a,b,fm_w,fm_s,r_oldx,r_oldy
      RXM = RXM*XSIGN
      RYM = RYM*YSIGN
      A      = AM / M_W
      B      = BM / M_S
      MNEWX   =  M+AM
      MNEWY   =  M+BM
      MNEW    =  M+AM+BM
      FM_W   = A*(RM_W + (1.-A)*RXM_W*XSIGN)
      FM_S   = B*(RM_S + (1.-B)*RYM_S*YSIGN)
      R_OLDX = RM/M
      R_OLDY = RM/M
      R_OLDX = .5*(R_OLDX + (RM+FM_S)/MNEWY)
      R_OLDY = .5*(R_OLDY + (RM+FM_W)/MNEWX)
      RM  = RM + FM_W + FM_S
      RXM = RXM + B*RXM_S*XSIGN
      RYM = RYM + A*RYM_W*YSIGN
      RXM = (RXM + 3.*(AM*R_OLDX - FM_W)) * MNEWY/MNEW
      RYM = (RYM + 3.*(BM*R_OLDY - FM_S)) * MNEWX/MNEW
      RXM = RXM*XSIGN
      RYM = RYM*YSIGN
      return
      end subroutine lusin_2sides

      subroutine susin_2sides(am,bm,m,rm,m_w,rm_w,m_s,rm_s)
      implicit none
      real*8 :: am,bm,m,rm,m_w,rm_w,m_s,rm_s
      real*8 :: a,b,fm_w,fm_s
      A      = AM / M_W
      B      = BM / M_S
      RM  = RM + A*RM_W + B*RM_S
      return
      end subroutine susin_2sides

      end module tqus_com
