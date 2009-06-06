      module cs2ll_utils
!@auth M. Kelley
!@ver 1.0
!@sum  cs2ll_utils provides procedures for interpolating
!@+    from a cubed-sphere to a lat-lon grid, and MPI-based
!@+    communication procedures for averaging and gathering
!@+    data around latitude circles.
!@+    Currently, only the gnomonic cubed-sphere geometry is
!@+    implemented.
!@+    The communication procedures in this module are applicable
!@+    to grids other than the cubed sphere.
!@+    This module depends upon the dd2d_utils module.
!@+    
!@usage Access is through the following types/routines.
!@+     Below, [ij][se]d refer to the minimum/maximum ilon/jlat
!@+        indices lying in the cubed-sphere domain of the PE.
!@+        They are elements of cs2ll_type.
!@+
!@+     cs2ll_type: a derived type containing communication and
!@+     interpolation info, initialized by calling subroutine
!@+        init_cs2ll_type(
!@+     &     grid,        ! an instance of dist_grid for the cubed sphere
!@+     &     imlon,jmlat, ! number of longitudes and latitudes
!@+     &     cs2ll,       ! the instance of cs2ll_type to be initialized
!@+     &     lons,lats    ! 1D arrays containing longitudes and latitudes
!@+     &     )
!@+
!@+     pack_zonal(cs2ll,jlat,nl,arr,arr_il)
!@+        gathers local array arr(isd:ied,1:nl) to a "global" array
!@+        arr_il(1:imlon,1:nl) on the root processor for latitude index
!@+        jlat.  Data is taken only between the first and last
!@+        valid longitudes at jlat.
!@+
!@+     sum_zonal(cs2ll,jlat,nl,arr,arrsum)
!@+        sums local array arr(isd:ied,1:nl) over longitude and
!@+        combines the partial sum with those from other processors
!@+        having latitude index jlat, placing the result arrsum(1:nl)
!@+        on every processor at jlat. Partial sums are only between the
!@+        first and last valid longitudes at jlat.
!@+
!@+     interp_xy_4D(grid,cs2ll,arrcs,arrll,nl,nk)
!@+     interp_xy_3D(grid,cs2ll,arrcs,arrll,nk)
!@+        Interpolate a cubed-sphere array arrcs(nl,isd:ied,jsd:jed,nk)
!@+             to a       lat-lon    array arrll(nl,isd:ied,jsd:jed,nk).
!@+        The 3D version omits the nl dimension.
!@+        For arrcs, [ij][se]d are the grid%  instances.
!@+        For arrll, [ij][se]d are the cs2ll% instances.
!@+        Interpolation is bilinear in nondimensional cube coords x,y.
!@+        These routines have been coded but not yet tested.
!@+
!@+     Parallelization efficiency is greatly improved by employing
!@+      the following information in cs2ll_type:
!@+     1. Rather than looping from jlat=jsd to jlat=jed,
!@+        follow the order given by jlat_sched:
!@+        do jj=1,jmlat
!@+           jlat = cs2ll_type%jlat_sched(jj)
!@+           if(cs2ll_type%ni(jlat).eq.0) cycle ! no valid lons at this lat
!@+        enddo
!@+        This allows latitudes with more participating PEs to
!@+        be scheduled later, maximizing the initial degree of
!@+        PE independence.
!@+     2. Root status is assigned in a round-robin fashion for each
!@+        band of latitudes having the same set of participating PEs.
!@+        Serial calculations (e.g. FFTs) performed on data obtained
!@+        via pack_zonal can be deferred until all possible PEs
!@+        for a given band of latitudes have received data upon which to
!@+        perform the serial operation.  The test jlat_calc(jlat)>0
!@+        indicates whether to proceed with the serial calculation
!@+        for the latitude at which pack_zonal sent data to my PE.
!@+        The value of jlat_calc(jlat) contains the latitude index
!@+        at which my PE received data.

      implicit none
      private
      include 'mpif.h'

      public :: cs2ll_type
      type cs2ll_type
        integer :: imlon ! number of longitudes
        integer :: jmlat ! number of latitudes
        integer :: isd ! minimum ilon in my domain
        integer :: ied ! maximum ilon in my domain
        integer :: jsd ! minimum jlat in my domain
        integer :: jed ! maximum jlat in my domain
! is(jlat) = first valid ilon at latitude jlat
        integer, dimension(:), allocatable :: is
! ie(jlat) = last valid ilon at latitude jlat
        integer, dimension(:), allocatable :: ie
! ni(jlat) = number of valid ilons at latitude jlat
! NOTE WRAPAROUND ISSUES!
        integer, dimension(:), allocatable :: ni
! comm(jlat) = MPI communicator for latitude jlat
        integer, dimension(:), allocatable :: comm
! nproc_j(jlat) = number of processors at latitude jlat
        integer, dimension(:), allocatable :: nproc_j
! am_i_rootj(jlat) = am i the root PE for latitude jlat
        logical, dimension(:), allocatable :: am_i_rootj
! jlat_sched(j) the jth jlat to be done
        integer, dimension(:), allocatable :: jlat_sched
! jlat_calc(jlat) = upon reaching jlat, do calculations for jlat_calc
        integer, dimension(:), allocatable :: jlat_calc
! isr(iproc,jlat) = first valid ilon at latitude jlat for PE iproc
! ier(iproc,jlat) = last valid ilon at latitude jlat for PE iproc
! these are only set for latitudes at which I am root.
        integer, dimension(:,:), allocatable :: isr
        integer, dimension(:,:), allocatable :: ier
!
        integer :: nsendmax ! = 1+ied-isd
        integer :: nrecvmax ! maximum number of lons root will receive

! ix,jy(ilon,jlat) = i,j position of lon(ilon),lat(jlat) on my grid
! these are real*8 numbers used for interpolations
! ix = i.x means a fraction x of the distance between i and i+1
! jy = j.y means a fraction y of the distance between j and j+1
        real*8, dimension(:,:), allocatable :: ix
        real*8, dimension(:,:), allocatable :: jy

      end type cs2ll_type

      public :: init_cs2ll_type,pack_zonal,sum_zonal,sumxpe_zonal

      contains

      subroutine init_cs2ll_type(grid,imlon,jmlat,cs2ll,lons,lats)
      use dd2d_utils, only : dist_grid
      use geom, only : ll2csxy
      type(dist_grid), intent(in) :: grid
      type(cs2ll_type), intent(out) :: cs2ll
      integer, intent(in) :: imlon,jmlat
      real*8, dimension(imlon), intent(in) :: lons
      real*8, dimension(jmlat), intent(in) :: lats
      integer :: ilon,jlat,jj,tile,n,nmax,jlat1,jlat2,j_sv,nrecv
      integer :: ierr,group_world,group_jlat,comm_jlat
      real*8 :: cs_xmin,cs_xmax,cs_ymin,cs_ymax,x,y
      integer, dimension(:), allocatable :: pelist
      integer, dimension(:,:), allocatable :: pelist_sv
      real*8, dimension(:), allocatable :: hash_pelist

      allocate(pelist(grid%nproc),pelist_sv(grid%nproc,jmlat))
      allocate(hash_pelist(jmlat))

c compute the x,y bounds of this processor
      cs_xmin = -1d0 + 2d0*(grid%is-1)/grid%npx
      cs_xmax = -1d0 + 2d0*(grid%ie  )/grid%npx
      cs_ymin = -1d0 + 2d0*(grid%js-1)/grid%npx
      cs_ymax = -1d0 + 2d0*(grid%je  )/grid%npx

      cs2ll%imlon = imlon
      cs2ll%jmlat = jmlat

c allocate cs2ll arrays
      allocate(cs2ll%comm(jmlat)
     &        ,cs2ll%ni(jmlat)
     &        ,cs2ll%is(jmlat)
     &        ,cs2ll%ie(jmlat)
     &        ,cs2ll%nproc_j(jmlat)
     &        ,cs2ll%am_i_rootj(jmlat)
     &        ,cs2ll%jlat_sched(jmlat)
     &        ,cs2ll%jlat_calc(jmlat)
     &     )

c
c find which latlon points are in our domain
c
      cs2ll%jsd = jmlat+1
      cs2ll%jed= -1
      do jlat=1,jmlat
        cs2ll%ni(jlat) = 0
        cs2ll%is(jlat) = imlon+1 
        cs2ll%ie(jlat) = -1
        do ilon=1,imlon
          call ll2csxy(lons(ilon),lats(jlat),x,y,tile)
          if(tile.ne.grid%tile) cycle
          if(x.lt.cs_xmin .or. x.gt.cs_xmax) cycle
          if(y.lt.cs_ymin .or. y.gt.cs_ymax) cycle
          cs2ll%ni(jlat) = cs2ll%ni(jlat) + 1
          cs2ll%is(jlat)=min(cs2ll%is(jlat),ilon)
          cs2ll%ie(jlat)=max(cs2ll%ie(jlat),ilon)
        enddo
        if(cs2ll%ni(jlat).gt.0) then
          cs2ll%jsd=min(cs2ll%jsd,jlat)
          cs2ll%jed=max(cs2ll%jed,jlat)
        endif
      enddo

c
c allocate and fill ix,jy now that we know the ilon,jlat extent of our domain
c
      cs2ll%isd = minval(cs2ll%is)
      cs2ll%ied = maxval(cs2ll%ie)
      allocate(cs2ll%ix(cs2ll%isd:cs2ll%ied,cs2ll%jsd:cs2ll%jed))
      allocate(cs2ll%jy(cs2ll%isd:cs2ll%ied,cs2ll%jsd:cs2ll%jed))
      do jlat=cs2ll%jsd,cs2ll%jed
      do ilon=cs2ll%isd,cs2ll%ied
        cs2ll%ix(ilon,jlat) = 1d30
        cs2ll%jy(ilon,jlat) = 1d30
        call ll2csxy(lons(ilon),lats(jlat),x,y,tile)
        if(tile.ne.grid%tile) cycle
        if(x.lt.cs_xmin .or. x.gt.cs_xmax) cycle
        if(y.lt.cs_ymin .or. y.gt.cs_ymax) cycle
        cs2ll%ix(ilon,jlat) = .5d0*(1d0+x)*grid%npx+.5d0
        cs2ll%jy(ilon,jlat) = .5d0*(1d0+y)*grid%npx+.5d0
      enddo
      enddo

c
c count how many processors are at each jlat
c
      pelist_sv = -1
      do jlat=1,jmlat
        if(cs2ll%ni(jlat).gt.0) then
          n = grid%gid
        else
          n = -1
        endif
        pelist = -1
        call mpi_allgather(n,1,MPI_INTEGER,pelist,1,
     &       MPI_INTEGER,MPI_COMM_WORLD,ierr)
        n = count(pelist.ge.0)
        cs2ll%nproc_j(jlat) = n
        pelist_sv(1:n,jlat) = pack(pelist,pelist.ge.0)
        hash_pelist(jlat) = r8hash(pelist_sv(1:n,jlat),n)
      enddo

c
c allocate isr,ier now that we know nproc_j
c only allocate over the latitudes in our domain
c
      nmax = maxval(cs2ll%nproc_j(cs2ll%jsd:cs2ll%jed))
      allocate(cs2ll%isr(nmax,cs2ll%jsd:cs2ll%jed))
      allocate(cs2ll%ier(nmax,cs2ll%jsd:cs2ll%jed))

c
c schedule the jlat sequence in order of increasing nproc_j
c
      jj = 0
      do n=1,maxval(cs2ll%nproc_j)
        do jlat=1,jmlat
          if(cs2ll%nproc_j(jlat).eq.n) then
            jj = jj + 1
            cs2ll%jlat_sched(jj) = jlat
          endif
        enddo
      enddo

c
c rotate the PE ordering in each band of latitudes having the same PEs.
c decide the latitudes at which the PEs in each group will do calculations
c
      cs2ll%jlat_calc(:) = -1
      jlat1 = 1
      do while(jlat1.le.jmlat)
        jlat2 = jlat1
        do while(jlat2.lt.jmlat)
          if(hash_pelist(jlat2+1).ne.hash_pelist(jlat1)) exit
          jlat2 = jlat2 + 1
        enddo
        n = cs2ll%nproc_j(jlat1)
        j_sv = -1
        do jlat=jlat1,jlat2
          if(jlat.gt.jlat1 .and. n.gt.1) then
            pelist_sv(2:n,jlat) = pelist_sv(1:n-1,jlat-1)
            pelist_sv(1,jlat) = pelist_sv(n,jlat-1)
          endif
          if(grid%gid.eq.pelist_sv(1,jlat)) then
            cs2ll%am_i_rootj(jlat) = .true.
            j_sv = jlat
          else
            cs2ll%am_i_rootj(jlat) = .false.
          endif
          if(mod(1+jlat-jlat1,n).eq.0 .or. jlat.eq.jlat2) then
            cs2ll%jlat_calc(jlat) = j_sv
            j_sv = -1
          endif
        enddo
        jlat1 = jlat2 + 1
      enddo

c
c create the communicator for each jlat
c
      call mpi_comm_group(MPI_COMM_WORLD,group_world,ierr)
      do jlat=1,jmlat
        n = cs2ll%nproc_j(jlat)
        pelist(1:n) = pelist_sv(1:n,jlat)

        call mpi_group_incl(group_world,n,pelist,group_jlat,ierr)
        call mpi_comm_create(MPI_COMM_WORLD,group_jlat,
     &         comm_jlat,ierr)

        if(cs2ll%ni(jlat).gt.0) then
          cs2ll%comm(jlat) = comm_jlat
          call mpi_gather(cs2ll%is(jlat),1,MPI_INTEGER,
     &         cs2ll%isr(1,jlat),  1,MPI_INTEGER,0,comm_jlat,ierr)
          call mpi_gather(cs2ll%ie(jlat),1,MPI_INTEGER,
     &         cs2ll%ier(1,jlat),  1,MPI_INTEGER,0,comm_jlat,ierr)
        else
          cs2ll%comm(jlat) = MPI_COMM_NULL
        endif

      enddo

c
c some useful send and receive info
c
      cs2ll%nsendmax = 1+cs2ll%ied-cs2ll%isd
      cs2ll%nrecvmax = 0
      do jlat=cs2ll%jsd,cs2ll%jed
        if(cs2ll%am_i_rootj(jlat)) then
          nrecv = 0
          do n=1,cs2ll%nproc_j(jlat)
            nrecv = nrecv + 1+cs2ll%ier(n,jlat)-cs2ll%isr(n,jlat)
          enddo
          cs2ll%nrecvmax=max(cs2ll%nrecvmax,nrecv)
        endif
      enddo

      deallocate(pelist,pelist_sv,hash_pelist)

      return
      end subroutine init_cs2ll_type

      function r8hash(ivec,n)
c hash a sequence of n nonnegative integers
c into a single real*8 value
c use the fortran random_number machinery
c the particular value of the hash will be system-dependent,
c but equality of hashes for different inputs is system-independent
      implicit none
      integer :: n
      integer :: ivec(n)
      real*8 :: r8hash
      integer :: i,j,seed_size
      integer, dimension(5) :: seeds=(/ 1, 3, 5, 7, 11 /)
      real*8 :: randn
      call random_seed(size=seed_size)
      call random_seed(put=seeds(1:seed_size))
      r8hash = 0d0
      do i=1,n
        do j=1,1+ivec(i)
          call random_number(randn)
        enddo
        r8hash = r8hash + randn
      enddo
      return
      end function r8hash

      subroutine pack_zonal(cs2ll,jlat,nl,arr,arr_il)
      implicit none
      type(cs2ll_type) :: cs2ll
      integer :: jlat,nl
      real*8 arr(cs2ll%isd:cs2ll%ied,nl)
      real*8 arr_il(cs2ll%imlon,nl)
      real*8 bufsend(cs2ll%nsendmax*nl)
      real*8 bufrecv(cs2ll%nrecvmax*nl)
      integer :: i,l,m,iproc,nsend,ierr,nproc
      integer, dimension(:), allocatable :: cnts,displs
      nproc = cs2ll%nproc_j(jlat)
      if(nproc.eq.1) then
        arr_il(:,:) = arr(:,:)
        return
      endif
c  pack the send message into contiguous memory
      m = 0
      do l=1,nl
      do i=cs2ll%is(jlat),cs2ll%ie(jlat)
        m = m + 1
        bufsend(m) = arr(i,l)
      enddo
      enddo
      nsend = m
      if(cs2ll%am_i_rootj(jlat)) then
        allocate(cnts(nproc),displs(nproc))
        do iproc=1,nproc
          cnts(iproc)=1+cs2ll%ier(iproc,jlat)-cs2ll%isr(iproc,jlat)
          cnts(iproc) = cnts(iproc)*nl
        enddo
        displs(1)=0
        do iproc=2,nproc
          displs(iproc) = cnts(iproc-1) + displs(iproc-1)
        enddo
      endif
      call mpi_gatherv(bufsend,nsend,MPI_DOUBLE_PRECISION,
     &     bufrecv,cnts,displs,MPI_DOUBLE_PRECISION,0,
     &     cs2ll%comm(jlat),ierr)
      if(cs2ll%am_i_rootj(jlat)) then
        deallocate(cnts,displs)
c  unpack the messages from the processors at this latitude.
c  because of wraparound, summation must be used.
        arr_il(:,:) = 0d0
        m = 0
        do iproc=1,nproc
          do l=1,nl
          do i=cs2ll%isr(iproc,jlat),cs2ll%ier(iproc,jlat)
            m = m + 1
            arr_il(i,l) = arr_il(i,l) + bufrecv(m)
          enddo
          enddo
        enddo
      endif
      return
      end subroutine pack_zonal

      subroutine sum_zonal(cs2ll,jlat,nl,arr,arrsum)
      implicit none
      type(cs2ll_type) :: cs2ll
      integer :: jlat,nl
      real*8 arr(cs2ll%isd:cs2ll%ied,nl)
      real*8 arrsum(nl)
      real*8 sum_partial(nl)
      integer :: ierr
      if(cs2ll%nproc_j(jlat).eq.1) then
c no communications required at this latitude
        arrsum(:)=sum(arr(cs2ll%is(jlat):cs2ll%ie(jlat),:),1)
      else
c compute the partial sum on my processor
        sum_partial(:)=sum(arr(cs2ll%is(jlat):cs2ll%ie(jlat),:),1)
c sum over the processors at this latitude
        call mpi_allreduce(sum_partial,arrsum,nl,MPI_DOUBLE_PRECISION,
     &       MPI_SUM,cs2ll%comm(jlat),ierr)
      endif
      return
      end subroutine sum_zonal

      subroutine sumxpe_zonal(cs2ll,jlat,nl,arr,arrsum)
      implicit none
      type(cs2ll_type) :: cs2ll
      integer :: jlat,nl
      real*8, intent(inout) ::  arr(nl)
      real*8, intent(out), optional :: arrsum(nl)
      real*8 tmparr(nl)
      integer :: ierr
      if(cs2ll%nproc_j(jlat).eq.1) then
c no communications required at this latitude
        if(present(arrsum)) arrsum(:)=arr(:)
      else
        if(present(arrsum)) then
          call mpi_allreduce(arr,arrsum,nl,MPI_DOUBLE_PRECISION,
     &         MPI_SUM,cs2ll%comm(jlat),ierr)
        else
          tmparr = arr
          call mpi_allreduce(tmparr,arr,nl,MPI_DOUBLE_PRECISION,
     &         MPI_SUM,cs2ll%comm(jlat),ierr)
        endif
      endif
      return
      end subroutine sumxpe_zonal

      end module cs2ll_utils

      subroutine interp_to_jlat_4D(grid,cs2ll,arrcs,arrll,nl,nk,jlat)
c note: halo cells of arrcs are assumed to have been filled!
      use dd2d_utils, only : dist_grid
      use cs2ll_utils, only : cs2ll_type
      implicit none
      type(dist_grid), intent(in) :: grid
      type(cs2ll_type), intent(in) :: cs2ll
      integer :: nl,nk,jlat
      real*8, dimension(nl,grid%isd:grid%ied,grid%jsd:grid%jed,nk)
     &     :: arrcs
      real*8, dimension(nl,cs2ll%isd:cs2ll%ied,nk) :: arrll
      integer :: l,ilon,k,i,j
      real*8 :: x,y,wti,wtj

      if(jlat.lt.cs2ll%jsd) return
      if(jlat.gt.cs2ll%jed) return

      do k=1,nk
      do ilon=cs2ll%isd,cs2ll%ied
        x = cs2ll%ix(ilon,jlat)
        y = cs2ll%jy(ilon,jlat)
        if(x.eq.1d30 .or. y.eq.1d30) then
          arrll(:,ilon,k) = 0d0
          cycle
        endif
        i = x
        j = y
        wti = x-i
        wtj = y-j
        do l=1,nl
          arrll(l,ilon,k) =
     &         wtj*(wti*arrcs(l,i+1,j+1,k)+(1.-wti)*arrcs(l,i,j+1,k))
     &   +(1.-wtj)*(wti*arrcs(l,i+1,j  ,k)+(1.-wti)*arrcs(l,i,j  ,k))
        enddo
      enddo
      enddo
      return
      end subroutine interp_to_jlat_4D

      subroutine interp_to_jlat_3D(grid,cs2ll,arrcs,arrll,nk,jlat)
c note: halo cells of arrcs are assumed to have been filled!
      use dd2d_utils, only : dist_grid
      use cs2ll_utils, only : cs2ll_type
      implicit none
      type(dist_grid), intent(in) :: grid
      type(cs2ll_type), intent(in) :: cs2ll
      integer :: nk,jlat
      real*8, dimension(grid%isd:grid%ied,grid%jsd:grid%jed,nk)
     &     :: arrcs
      real*8, dimension(cs2ll%isd:cs2ll%ied,nk) :: arrll
      integer :: ilon,k,i,j
      real*8 :: x,y,wti,wtj

      if(jlat.lt.cs2ll%jsd) return
      if(jlat.gt.cs2ll%jed) return

      do k=1,nk
      do ilon=cs2ll%isd,cs2ll%ied
        x = cs2ll%ix(ilon,jlat)
        y = cs2ll%jy(ilon,jlat)
        if(x.eq.1d30 .or. y.eq.1d30) then
          arrll(ilon,k) = 0d0
          cycle
        endif
        i = x
        j = y
        wti = x-i
        wtj = y-j
        arrll(ilon,k) =
     &         wtj*(wti*arrcs(i+1,j+1,k)+(1.-wti)*arrcs(i,j+1,k))
     &   +(1.-wtj)*(wti*arrcs(i+1,j  ,k)+(1.-wti)*arrcs(i,j  ,k))
      enddo
      enddo
      return
      end subroutine interp_to_jlat_3D

      subroutine corner_fill_4D(grid,arr,nl,nk)
c NOTE: valid halo cells of arr are assumed to have been filled!
c Fill imaginary halo cells at the corners of the cube with values qf
c such that bilinear interpolation in cube coordinates x,y produces
c results at cube corners that are equal to the average of the 3
c surrounding cells.  An example for the upper right corner of a
c cube face is drawn here.
c
c     ------                Value at cube corner C =
c    |      |               (q1+q2+q3)/3 = (q1+q2+q3+qf)/4
c    |  q1  |  qf           -> qf = (q1+q2+q3)/3
c    |      |
c     ------C------ 
c    |      |      |
c    |  q2  |  q3  |
c    |      |      |
c     ------ ------
c
      use dd2d_utils, only : dist_grid
      implicit none
      type(dist_grid), intent(in) :: grid
      integer :: nl,nk
      real*8, dimension(nl,grid%isd:grid%ied,grid%jsd:grid%jed,nk)
     &     :: arr
      integer :: i,j,k,l,m,n
      integer :: ncor
      integer, dimension(4) :: ifill,jfill
      integer, dimension(3,4) :: icor,jcor

      call corner_fill_info(grid,ncor,icor,jcor,ifill,jfill)
      if(ncor.eq.0) return

      do k=1,nk
        do m=1,ncor
          i = ifill(m); j = jfill(m)
          do l=1,nl
            arr(l,i,j,k) = 0.
            do n=1,3
              arr(l,i,j,k) = arr(l,i,j,k) + arr(l,icor(n,m),jcor(n,m),k)
            enddo
            arr(l,i,j,k) = arr(l,i,j,k)/3d0
          enddo
        enddo
      enddo

      return
      end subroutine corner_fill_4D

      subroutine corner_fill_3D(grid,arr,nk)
c like corner_fill_4D without the "l" dimension.
      use dd2d_utils, only : dist_grid
      implicit none
      type(dist_grid), intent(in) :: grid
      integer :: nk
      real*8, dimension(grid%isd:grid%ied,grid%jsd:grid%jed,nk) :: arr
      integer :: i,j,k,m,n
      integer :: ncor
      integer, dimension(4) :: ifill,jfill
      integer, dimension(3,4) :: icor,jcor

      call corner_fill_info(grid,ncor,icor,jcor,ifill,jfill)
      if(ncor.eq.0) return

      do k=1,nk
        do m=1,ncor
          i = ifill(m); j = jfill(m)
          arr(i,j,k) = 0.
          do n=1,3
            arr(i,j,k) = arr(i,j,k) + arr(icor(n,m),jcor(n,m),k)
          enddo
          arr(i,j,k) = arr(i,j,k)/3d0
        enddo
      enddo

      return
      end subroutine corner_fill_3D

      subroutine corner_fill_info(grid,ncor,icor,jcor,ifill,jfill)
      use dd2d_utils, only : dist_grid
      implicit none
      type(dist_grid), intent(in) :: grid
      integer, intent(out) :: ncor
      integer, dimension(3,4), intent(out) :: icor,jcor
      integer, dimension(4), intent(out) :: ifill,jfill

      ncor = 0
      if(grid%js.eq.1) then
        if(grid%is.eq.1) then
          ncor = ncor + 1
          ifill(ncor) = 0; jfill(ncor) = 0;
          icor(:,ncor) = (/ 1, 0, 1 /)
          jcor(:,ncor) = (/ 0, 1, 1 /)
        endif
        if(grid%ie.eq.grid%npx) then
          ncor = ncor + 1
          ifill(ncor) = grid%ie+1; jfill(ncor) = 0;
          icor(:,ncor) = (/ grid%ie, grid%ie, grid%ie+1 /)
          jcor(:,ncor) = (/ 0, 1, 1 /)
        endif
      endif
      if(grid%je.eq.grid%npy) then
        if(grid%is.eq.1) then
          ncor = ncor + 1
          ifill(ncor) = 0; jfill(ncor) = grid%npy+1;
          icor(:,ncor) = (/ 0, 1, 1 /)
          jcor(:,ncor) = (/ grid%npy, grid%npy, grid%npy+1 /)
        endif
        if(grid%ie.eq.grid%npx) then
          ncor = ncor + 1
          ifill(ncor) = grid%npx+1; jfill(ncor) = grid%npy+1;
          icor(:,ncor) = (/ grid%npx, grid%npx+1, grid%npx /)
          jcor(:,ncor) = (/ grid%npy, grid%npy, grid%npy+1 /)
        endif
      endif

      return
      end subroutine corner_fill_info

