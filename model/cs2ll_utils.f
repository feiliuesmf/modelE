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
!@+    This module is also the TEMPORARY home of the
!@+    geometry-independent xgridremap procedure.
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
!@+     interp_to_jlat_4D(grid,cs2ll,arrcs,arrll,nl,nk,jlat)
!@+     interp_to_jlat_3D(grid,cs2ll,arrcs,arrll,nk,jlat)
!@+        Interpolate a cubed-sphere array arrcs(nl,isd:ied,jsd:jed,nk)
!@+             to a      longitude   array arrll(nl,isd:ied,nk) for
!@+             lat index jlat.
!@+        The 3D version omits the nl dimension.
!@+        For arrcs, [ij][se]d are the grid%  instances.
!@+        For arrll, [i][se]d  are the cs2ll% instances.
!@+        Interpolation is bilinear in nondimensional cube coords x,y.
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

! ix,jy(ilon,jlat) = i,j position of lon(ilon),lat(jlat) on CS grid
! these are real*8 numbers used for interpolations
! ix = i.x means a fraction x of the distance between i and i+1
! jy = j.y means a fraction y of the distance between j and j+1
        real*8, dimension(:,:), allocatable :: ix
        real*8, dimension(:,:), allocatable :: jy

      end type cs2ll_type

      public :: init_cs2ll_type,pack_zonal,sum_zonal,sumxpe_zonal

      type redist_type
        integer :: npsrc,npdst,nproc,nproc_tile
        integer :: isd,ied,jsd,jed,npx,npy,ntiles
        integer, dimension(:), allocatable :: js_pack,je_pack
        integer, dimension(:,:), allocatable :: nb_pack
        integer, dimension(:,:,:), allocatable :: is_pack,ie_pack
        integer, dimension(:), allocatable :: js_unpack,je_unpack
        integer, dimension(:,:), allocatable :: nb_unpack
        integer, dimension(:,:,:), allocatable :: is_unpack,ie_unpack
        integer, dimension(:), allocatable :: send_cnts,send_displs
        integer, dimension(:), allocatable :: recv_cnts,recv_displs
      end type redist_type

c documentation on these to be added when interfaces are finalized
      public :: xgridremap_type
      public :: init_xgridremap_type,xgridremap_ij,xgridremap_lij
      type xgridremap_type
        type(redist_type) :: redist
        integer :: is,ie,js,je
        integer :: isd,ied,jsd,jed
        integer :: mpoly
        integer, dimension(:,:), allocatable :: kpoly
        integer, dimension(:), allocatable :: recv_index
        real*8, dimension(:), allocatable :: wts
      end type xgridremap_type

c documentation on these to be added when interfaces are finalized
      public :: cs2llint_type
      public :: init_cs2llint_type,cs2llint_ij,cs2llint_lij
      type cs2llint_type
        integer :: isdcs,iedcs,jsdcs,jedcs
        integer ::       imlon,jsdll,jedll
        integer :: nproc
        integer :: npts_interp,npts_unpack
        real*8, dimension(:), allocatable :: ix,jy
        integer, dimension(:), allocatable :: i_unpack,j_unpack
        integer, dimension(:), allocatable :: send_cnts,send_displs
        integer, dimension(:), allocatable :: recv_cnts,recv_displs
      end type cs2llint_type

      public :: ll2csint_type
      public :: init_ll2csint_type,ll2csint_ij,ll2csint_lij
      type ll2csint_type
        integer :: iscs,iecs,jscs,jecs
        integer :: imlon,jsdll,jedll,jlat_min,jlat_max
        integer :: nproc
        integer :: npts_interp,npts_unpack
        real*8, dimension(:), allocatable :: ix,jy
        integer, dimension(:), allocatable :: i_unpack,j_unpack
        integer, dimension(:), allocatable :: send_cnts,send_displs
        integer, dimension(:), allocatable :: recv_cnts,recv_displs
      end type ll2csint_type

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

      subroutine init_redist_type(grid_src,grid_dst,
     &     ndepend, isrc,jsrc,tsrc, idst,jdst,tdst,
     &     redist)
      use dd2d_utils, only : dist_grid
      type(dist_grid), intent(in) :: grid_src,grid_dst
      integer, intent(in) :: ndepend
      integer, intent(in), dimension(ndepend) ::
     &     isrc,jsrc,tsrc, idst,jdst,tdst
      type(redist_type), intent(out) :: redist
c local vars
      integer :: is,ie,js,je,npsrc,npdst,nproc,npx,npy,nbmax
      integer, dimension(:), allocatable ::
     &     isdst,iedst,jsdst,jedst,tiledst
      integer, dimension(:,:,:), allocatable :: itmp
      integer :: b,i,j,m,n,p,ierr
      logical, dimension(:,:,:), allocatable :: qmask

c for now, assume src pe set same as dst pe set
      redist%nproc  = grid_src%nproc
      redist%npsrc  = grid_src%nproc
      redist%npdst  = grid_dst%nproc
      redist%nproc_tile  = grid_src%nproc_tile
      redist%ntiles = grid_src%ntiles
      redist%npx    = grid_src%npx
      redist%npy    = grid_src%npy
      redist%isd    = grid_src%isd
      redist%ied    = grid_src%ied
      redist%jsd    = grid_src%jsd
      redist%jed    = grid_src%jed

      is = grid_src%is
      ie = grid_src%ie
      js = grid_src%js
      je = grid_src%je
      npsrc = redist%npsrc
      npdst = redist%npdst
      nproc = redist%nproc
      npx = redist%npx
      npy = redist%npy

c collect info about the domain decomp of grid_dst
c for now, assume src pe set same as dst pe set
      allocate(isdst(npdst),iedst(npdst),jsdst(npdst),jedst(npdst),
     &     tiledst(npdst))

      call mpi_allgather(grid_dst%is,1,MPI_INTEGER,isdst,1,
     &     MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call mpi_allgather(grid_dst%ie,1,MPI_INTEGER,iedst,1,
     &     MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call mpi_allgather(grid_dst%js,1,MPI_INTEGER,jsdst,1,
     &     MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call mpi_allgather(grid_dst%je,1,MPI_INTEGER,jedst,1,
     &     MPI_INTEGER,MPI_COMM_WORLD,ierr)
      call mpi_allgather(grid_dst%tile,1,MPI_INTEGER,tiledst,1,
     &     MPI_INTEGER,MPI_COMM_WORLD,ierr)

c
c tabulate which local gridcells to send to each destination pe,
c and the 2D-array counts/displacements for the send buffer
c
      allocate(redist%js_pack(npdst), redist%je_pack(npdst))
      allocate(redist%nb_pack(js:je,npdst))
      if(nproc == 6) then
        nbmax = 4 ! some LL PEs need 4 if npes=6 for ll2cs dir
      else
        nbmax = 2  
      endif
      allocate(redist%is_pack(nbmax,js:je,npdst),
     &         redist%ie_pack(nbmax,js:je,npdst))
      allocate(redist%send_cnts(nproc), redist%send_displs(nproc))
      redist%js_pack(:)   = redist%npy+1
      redist%je_pack(:)   = -1
      allocate(qmask(is:ie,js:je,npdst))
      qmask(:,:,:) = .false.
      do n=1,ndepend
        if(tsrc(n) .ne. grid_src%tile) cycle
        i = isrc(n)
        j = jsrc(n)
        if(i .lt. is .or. i .gt. ie) cycle
        if(j .lt. js .or. j .gt. je) cycle
        do p=1,npdst
          if(tdst(n) .ne. tiledst(p)) cycle
          if(idst(n) .lt. isdst(p) .or. idst(n) .gt. iedst(p)) cycle
          if(jdst(n) .lt. jsdst(p) .or. jdst(n) .gt. jedst(p)) cycle
          redist%js_pack(p) = min(redist%js_pack(p),j)
          redist%je_pack(p) = max(redist%je_pack(p),j)
          qmask(i,j,p) = .true.
          exit
        enddo
      enddo
      do p=1,npdst
        do j=redist%js_pack(p),redist%je_pack(p)
          call get_i1i2(qmask(is,j,p),1+ie-is,
     &         n,
     &         redist%is_pack(1,j,p),
     &         redist%ie_pack(1,j,p),nbmax)
          redist%nb_pack(j,p) = n
          redist%is_pack(1:n,j,p) = redist%is_pack(1:n,j,p) + is-1
          redist%ie_pack(1:n,j,p) = redist%ie_pack(1:n,j,p) + is-1
        enddo
      enddo
      deallocate(qmask,isdst,iedst,jsdst,jedst,tiledst)
      m = 0
      do p=1,npdst
        redist%send_displs(p) = m ! assuming src pe set same as dst pe set
        do j=redist%js_pack(p),redist%je_pack(p)
          do b=1,redist%nb_pack(j,p)
            do i=redist%is_pack(b,j,p),redist%ie_pack(b,j,p)
              m = m + 1
            enddo
          enddo
        enddo
      enddo
      do p=1,nproc-1
        redist%send_cnts(p) =
     &       redist%send_displs(p+1)-redist%send_displs(p)
      enddo
      redist%send_cnts(nproc) = m-redist%send_displs(nproc)

c
c set up the counts/displacements and unpack info for the recv buffer
c for now, assume src pe set same as dst pe set
c
      allocate(redist%js_unpack(npsrc), redist%je_unpack(npsrc))
      allocate(redist%nb_unpack(npy,npsrc))
      allocate(redist%is_unpack(nbmax,npy,npsrc),
     &         redist%ie_unpack(nbmax,npy,npsrc))
      allocate(redist%recv_cnts(nproc), redist%recv_displs(nproc))

      call mpi_alltoall(redist%send_cnts,1,MPI_INTEGER,
     &                  redist%recv_cnts,1,MPI_INTEGER,
     &                  MPI_COMM_WORLD,ierr)
      redist%recv_displs(1) = 0
      do p=2,npsrc
        redist%recv_displs(p) =
     &       redist%recv_displs(p-1)+redist%recv_cnts(p-1)
      enddo
      call mpi_alltoall(redist%js_pack,  1,MPI_INTEGER,
     &                  redist%js_unpack,1,MPI_INTEGER,
     &                  MPI_COMM_WORLD,ierr)
      call mpi_alltoall(redist%je_pack,  1,MPI_INTEGER,
     &                  redist%je_unpack,1,MPI_INTEGER,
     &                  MPI_COMM_WORLD,ierr)

      n = nbmax*npy
      allocate(itmp(nbmax,npy,npdst))
      do p=1,npdst
        itmp(:,:,p) = npx+1
        itmp(:,js:je,p) = redist%is_pack(:,js:je,p)
      enddo
      call mpi_alltoall(itmp,            n,MPI_INTEGER,
     &                  redist%is_unpack,n,MPI_INTEGER,
     &                  MPI_COMM_WORLD,ierr)
      do p=1,npdst
        itmp(:,:,p) = -1
        itmp(:,js:je,p) = redist%ie_pack(:,js:je,p)
      enddo
      call mpi_alltoall(itmp,            n,MPI_INTEGER,
     &                  redist%ie_unpack,n,MPI_INTEGER,
     &                  MPI_COMM_WORLD,ierr)
      deallocate(itmp)

      allocate(itmp(1,npy,npdst))
      do p=1,npdst
        itmp(:,:,p) = 0
        itmp(1,js:je,p) = redist%nb_pack(js:je,p)
      enddo
      call mpi_alltoall(itmp,            npy,MPI_INTEGER,
     &                  redist%nb_unpack,npy,MPI_INTEGER,
     &                  MPI_COMM_WORLD,ierr)
      deallocate(itmp)

      return
      end subroutine init_redist_type

      subroutine get_i1i2(q,im,n,i1,i2,nmax)
c Finds the number of intervals over which q==true, and the
c beginning/ending indices of each interval.
c Wraparound is disabled.
c Originally written for ocean basins.  Move to generic utilities file.
      implicit none
      integer, intent(in) :: im  ! number of points
      logical, dimension(im), intent(in) :: q ! logical mask
      integer, intent(in) :: nmax ! max number of contiguous intervals
      integer, intent(out) :: n ! number of contiguous intervals
      integer, dimension(nmax) :: i1,i2 ! start,end of each interval
      integer :: i,nn
c first find the number of intervals
      n = 0
      do i=1,im-1
        if(q(i) .and. .not.q(i+1)) n = n + 1
      enddo
      if(q(im)) n = n + 1
      if(n.gt.nmax) then
        write(6,*) 'for get_i1i2, increase nmax to ',n
        call stop_model('get_i1i2: n>nmax',255)
      endif
c find the start/end of each interval
      i = 1
      do nn=1,n
        do while(.not.q(i))
          i = i + 1
        enddo
        i1(nn) = i
        do while(q(i))
          i = i + 1
          if(i.gt.im) exit
        enddo
        i2(nn) = i-1
      enddo
      return
      end subroutine get_i1i2

      subroutine redist_data_ij(redist,arr,arr_glob)
      type(redist_type), intent(in) :: redist
      real*8, dimension(redist%isd:redist%ied,redist%jsd:redist%jed),
     &     intent(in) :: arr
      real*8, dimension(redist%npx,redist%npy,redist%ntiles),
     &     intent(out) :: arr_glob
c local vars
      integer :: b,i,j,m,n,p,tile,nproc,ierr
      real*8, dimension(:), allocatable :: bufsend,bufrecv
      integer, dimension(:), allocatable :: scnts,sdspl,rcnts,rdspl

c
c allocate send/recv buffers
c
      nproc = redist%nproc
      allocate(scnts(nproc),sdspl(nproc),rcnts(nproc),rdspl(nproc))
      scnts(:) = redist%send_cnts
      sdspl(:) = redist%send_displs
      rcnts(:) = redist%recv_cnts
      rdspl(:) = redist%recv_displs

      allocate(bufsend(sum(scnts)),bufrecv(sum(rcnts)))

c
c pack data into 1D send buffer
c
      m = 0
      do p=1,redist%npdst
        do j=redist%js_pack(p),redist%je_pack(p)
          do b=1,redist%nb_pack(j,p)
            do i=redist%is_pack(b,j,p),redist%ie_pack(b,j,p)
              m = m + 1
              bufsend(m) = arr(i,j)
            enddo
          enddo
        enddo
      enddo

c
c exchange data
c
      call mpi_alltoallv(bufsend, scnts, sdspl, MPI_DOUBLE_PRECISION,
     &                   bufrecv, rcnts, rdspl, MPI_DOUBLE_PRECISION,
     &                   MPI_COMM_WORLD, ierr)

c
c unpack 1D recv buffer into global array
c
      m = 0
      p = 0
      do tile=1,redist%ntiles
        do n=1,redist%nproc_tile
          p = p + 1
          do j=redist%js_unpack(p),redist%je_unpack(p)
            do b=1,redist%nb_unpack(j,p)
              do i=redist%is_unpack(b,j,p),redist%ie_unpack(b,j,p)
                m = m + 1
                arr_glob(i,j,tile) = bufrecv(m)
              enddo
            enddo
          enddo
        enddo
      enddo

c
c deallocate workspace
c
      deallocate(bufsend,bufrecv)
      deallocate(scnts,sdspl,rcnts,rdspl)

      return
      end subroutine redist_data_ij

      subroutine init_xgridremap_type(grid_src,grid_dst,
     &     npoly, isrc,jsrc,tsrc, idst,jdst,tdst, area,
     &     remap)
c
c Note: this routine needs i,j,t,area lists to be GLOBAL
c
      use dd2d_utils, only : dist_grid
      type(dist_grid), intent(in) :: grid_src,grid_dst
      integer, intent(in) :: npoly
      integer, intent(in), dimension(npoly) ::
     &     isrc,jsrc,tsrc, idst,jdst,tdst
      real*8, intent(in), dimension(npoly) :: area
      type(xgridremap_type), intent(out) :: remap
c local vars
      integer :: isd,ied,jsd,jed ! dst grid
      integer :: is,ie,js,je ! dst grid
      integer :: i,j,k,m,mm,n,p,npoly_loc, bb,ii,jj,kk,tt
      integer, dimension(:), allocatable ::
     &     isrc_loc,jsrc_loc,tsrc_loc, idst_loc,jdst_loc,
     &     isrc_srt,jsrc_srt,tsrc_srt
      real*8, dimension(:), allocatable :: area_loc,area_srt
      integer, parameter :: maxovlap=10000
      integer, dimension(maxovlap) :: isrc_ij,jsrc_ij,tsrc_ij
      real*8,  dimension(maxovlap) :: area_ij
      real*8 :: area_sum

c
c initialize the redist object for mpi send/recv info
c
      call init_redist_type(grid_src,grid_dst,
     &     npoly, isrc,jsrc,tsrc, idst,jdst,tdst,
     &     remap%redist)

      is = grid_dst%is
      ie = grid_dst%ie
      js = grid_dst%js
      je = grid_dst%je
      isd = grid_dst%isd
      ied = grid_dst%ied
      jsd = grid_dst%jsd
      jed = grid_dst%jed

      remap%is = grid_dst%is
      remap%ie = grid_dst%ie
      remap%js = grid_dst%js
      remap%je = grid_dst%je
      remap%isd = grid_dst%isd
      remap%ied = grid_dst%ied
      remap%jsd = grid_dst%jsd
      remap%jed = grid_dst%jed

c
c count the number of xgrid polygons in the domain of this destination
c PE and allocate the remap arrays having this size
c
      allocate(isrc_loc(npoly),jsrc_loc(npoly),tsrc_loc(npoly))
      allocate(idst_loc(npoly),jdst_loc(npoly),area_loc(npoly))
      m = 0
      do n=1,npoly
        if(tdst(n) .ne. grid_dst%tile) cycle
        i = idst(n)
        j = jdst(n)
        if(i .lt. is .or. i .gt. ie) cycle
        if(j .lt. js .or. j .gt. je) cycle
        m = m + 1
        isrc_loc(m) = isrc(n)
        jsrc_loc(m) = jsrc(n)
        tsrc_loc(m) = tsrc(n)
        idst_loc(m) = idst(n)
        jdst_loc(m) = jdst(n)
        area_loc(m) = area(n)
      enddo
      npoly_loc = m
      remap%mpoly = npoly_loc
      allocate(isrc_srt(npoly_loc),jsrc_srt(npoly_loc))
      allocate(tsrc_srt(npoly_loc),area_srt(npoly_loc))

c
c For reproducibility on different numbers of PEs:
c Sort the xgrid polygons into dst i,j order.
c At each dst i,j sort the polygons into src i,j,t order.
c Lazy sort OK during initialization.
c
      allocate(remap%kpoly(is:ie,js:je))
      m = 0
      do j=js,je
        do i=is,ie
          k = 0
          do n=1,npoly_loc
            if(j .ne. jdst_loc(n)) cycle
            if(i .ne. idst_loc(n)) cycle
            k = k + 1
            isrc_ij(k) = isrc_loc(n)
            jsrc_ij(k) = jsrc_loc(n)
            tsrc_ij(k) = tsrc_loc(n)
            area_ij(k) = area_loc(n)
          enddo
          remap%kpoly(i,j) = k
          do tt=minval(tsrc_ij(1:k)),maxval(tsrc_ij(1:k))
            do jj=1,grid_src%npy
              do ii=1,grid_src%npx
                do kk=1,k
                  if(tt .ne. tsrc_ij(kk)) cycle
                  if(jj .ne. jsrc_ij(kk)) cycle
                  if(ii .ne. isrc_ij(kk)) cycle
                  m = m + 1
                  isrc_srt(m) = isrc_ij(kk)
                  jsrc_srt(m) = jsrc_ij(kk)
                  tsrc_srt(m) = tsrc_ij(kk)
                  area_srt(m) = area_ij(kk)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

c
c Compute area weights.
c Find where each src i,j,t will be in the 1D recv buffer after
c mpi_alltoallv(). Lazy search mimics an unpack operation.
c
      allocate(remap%recv_index(npoly_loc))
      remap%recv_index(:) = -1
      allocate(remap%wts(npoly_loc))
      m = 0
      do j=js,je
      do i=is,ie
      area_sum = sum(area_srt(m+1:m+remap%kpoly(i,j)))
      do k=1,remap%kpoly(i,j)
        m = m + 1
        remap%wts(m) = area_srt(m)/area_sum
        mm = 0
        p = 0
        search_loop: do tt=1,remap%redist%ntiles
        do n=1,remap%redist%nproc_tile
        p = p + 1
        do jj=remap%redist%js_unpack(p),remap%redist%je_unpack(p)
        do bb=1,remap%redist%nb_unpack(jj,p)
        do ii=remap%redist%is_unpack(bb,jj,p),
     &        remap%redist%ie_unpack(bb,jj,p)
          mm = mm + 1
          if(tt .eq. tsrc_srt(m) .and.
     &       jj .eq. jsrc_srt(m) .and.
     &       ii .eq. isrc_srt(m)) then
            remap%recv_index(m) = mm
            exit search_loop
          endif
        enddo ! ii
        enddo ! bb
        enddo ! jj
        enddo ! n
        enddo search_loop ! tt
      enddo
      enddo
      enddo

c
c deallocate workspace
c
      deallocate(isrc_loc,jsrc_loc,tsrc_loc)
      deallocate(idst_loc,jdst_loc,area_loc)
      deallocate(isrc_srt,jsrc_srt,tsrc_srt,area_srt)

      return
      end subroutine init_xgridremap_type

      subroutine xgridremap_ij(remap,arr,arr_out)
      type(xgridremap_type), intent(in) :: remap
      real*8, dimension(remap%redist%isd:remap%redist%ied,
     &                  remap%redist%jsd:remap%redist%jed),
     &     intent(in) :: arr
      real*8, dimension(remap%isd:remap%ied,remap%jsd:remap%jed),
     &     intent(out) :: arr_out
c local vars
      integer :: b,i,j,k,m,n,p,nproc,ierr
      real*8, dimension(:), allocatable :: bufsend,bufrecv
      integer, dimension(:), allocatable :: scnts,sdspl,rcnts,rdspl

c
c allocate send/recv buffers
c
      nproc = remap%redist%nproc
      allocate(scnts(nproc),sdspl(nproc),rcnts(nproc),rdspl(nproc))
      scnts(:) = remap%redist%send_cnts
      sdspl(:) = remap%redist%send_displs
      rcnts(:) = remap%redist%recv_cnts
      rdspl(:) = remap%redist%recv_displs

      allocate(bufsend(sum(scnts)),bufrecv(sum(rcnts)))

c
c pack data into 1D send buffer
c
      m = 0
      do p=1,remap%redist%npdst
        do j=remap%redist%js_pack(p),remap%redist%je_pack(p)
          do b=1,remap%redist%nb_pack(j,p)
            do i=remap%redist%is_pack(b,j,p),remap%redist%ie_pack(b,j,p)
              m = m + 1
              bufsend(m) = arr(i,j)
            enddo
          enddo
        enddo
      enddo

c
c exchange data
c
      call mpi_alltoallv(bufsend, scnts, sdspl, MPI_DOUBLE_PRECISION,
     &                   bufrecv, rcnts, rdspl, MPI_DOUBLE_PRECISION,
     &                   MPI_COMM_WORLD, ierr)

c
c remap using the contents of the 1D receive buffer
c
      m = 0
      do j=remap%js,remap%je
        do i=remap%is,remap%ie
          arr_out(i,j) = 0.
          do k=1,remap%kpoly(i,j)
            m = m + 1
            n = remap%recv_index(m)
            arr_out(i,j) = arr_out(i,j) + bufrecv(n)*remap%wts(m)
          enddo
        enddo
      enddo
c
c deallocate workspace
c
      deallocate(bufsend,bufrecv)
      deallocate(scnts,sdspl,rcnts,rdspl)

      return
      end subroutine xgridremap_ij

      subroutine xgridremap_lij(remap,arr,arr_out)
      type(xgridremap_type), intent(in) :: remap
      real*8, dimension(:,remap%redist%isd:,remap%redist%jsd:),
     &     intent(in) :: arr
      real*8, dimension(:,remap%isd:,remap%jsd:),
     &     intent(out) :: arr_out
c local vars
      integer :: b,i,j,k,l,lm,m,n,p,nproc,ierr
      real*8, dimension(:), allocatable :: bufsend,bufrecv
      integer, dimension(:), allocatable :: scnts,sdspl,rcnts,rdspl

      lm = size(arr,1)
c
c allocate send/recv buffers
c
      nproc = remap%redist%nproc
      allocate(scnts(nproc),sdspl(nproc),rcnts(nproc),rdspl(nproc))
      scnts(:) = lm*(remap%redist%send_cnts)
      sdspl(:) = lm*(remap%redist%send_displs)
      rcnts(:) = lm*(remap%redist%recv_cnts)
      rdspl(:) = lm*(remap%redist%recv_displs)

      allocate(bufsend(sum(scnts)),bufrecv(sum(rcnts)))

c
c pack data into 1D send buffer
c
      m = 0
      do p=1,remap%redist%npdst
        do j=remap%redist%js_pack(p),remap%redist%je_pack(p)
          do b=1,remap%redist%nb_pack(j,p)
            do i=remap%redist%is_pack(b,j,p),remap%redist%ie_pack(b,j,p)
              do l=1,lm
                m = m + 1
                bufsend(m) = arr(l,i,j)
              enddo
            enddo
          enddo
        enddo
      enddo

c
c exchange data
c
      call mpi_alltoallv(bufsend, scnts, sdspl, MPI_DOUBLE_PRECISION,
     &                   bufrecv, rcnts, rdspl, MPI_DOUBLE_PRECISION,
     &                   MPI_COMM_WORLD, ierr)

c
c remap using the contents of the 1D receive buffer
c
      m = 0
      do j=remap%js,remap%je
        do i=remap%is,remap%ie
          arr_out(:,i,j) = 0.
          do k=1,remap%kpoly(i,j)
            m = m + 1
            n = lm*(remap%recv_index(m)-1)
            do l=1,lm
              n = n + 1
              arr_out(l,i,j) = arr_out(l,i,j) + bufrecv(n)*remap%wts(m)
            enddo
          enddo
        enddo
      enddo

c
c deallocate workspace
c
      deallocate(bufsend,bufrecv)
      deallocate(scnts,sdspl,rcnts,rdspl)

      return
      end subroutine xgridremap_lij

      subroutine xgridremap_ijl(remap,arr,arr_out)
      type(xgridremap_type), intent(in) :: remap
      real*8, dimension(remap%redist%isd:,remap%redist%jsd:,:),
     &     intent(in) :: arr
      real*8, dimension(remap%isd:,remap%jsd:,:),
     &     intent(out) :: arr_out
c local vars
      integer :: b,i,j,k,l,lm,m,n,p,nproc,ierr
      real*8, dimension(:), allocatable :: bufsend,bufrecv,arrl
      integer, dimension(:), allocatable :: scnts,sdspl,rcnts,rdspl

      lm = size(arr,3)
      allocate(arrl(lm))
c
c allocate send/recv buffers
c
      nproc = remap%redist%nproc
      allocate(scnts(nproc),sdspl(nproc),rcnts(nproc),rdspl(nproc))
      scnts(:) = lm*(remap%redist%send_cnts)
      sdspl(:) = lm*(remap%redist%send_displs)
      rcnts(:) = lm*(remap%redist%recv_cnts)
      rdspl(:) = lm*(remap%redist%recv_displs)

      allocate(bufsend(sum(scnts)),bufrecv(sum(rcnts)))

c
c pack data into 1D send buffer
c
      m = 0
      do p=1,remap%redist%npdst
        do j=remap%redist%js_pack(p),remap%redist%je_pack(p)
          do b=1,remap%redist%nb_pack(j,p)
            do i=remap%redist%is_pack(b,j,p),remap%redist%ie_pack(b,j,p)
              do l=1,lm
                m = m + 1
! against the grain for simplicity later
                bufsend(m) = arr(i,j,l)
              enddo
            enddo
          enddo
        enddo
      enddo

c
c exchange data
c
      call mpi_alltoallv(bufsend, scnts, sdspl, MPI_DOUBLE_PRECISION,
     &                   bufrecv, rcnts, rdspl, MPI_DOUBLE_PRECISION,
     &                   MPI_COMM_WORLD, ierr)

c
c remap using the contents of the 1D receive buffer
c
      m = 0
      do j=remap%js,remap%je
        do i=remap%is,remap%ie
          arrl(:) = 0.
          do k=1,remap%kpoly(i,j)
            m = m + 1
            n = lm*(remap%recv_index(m)-1)
            do l=1,lm
              n = n + 1
              arrl(l) = arrl(l) + bufrecv(n)*remap%wts(m)
            enddo
          enddo
          arr_out(i,j,:) = arrl(:)
        enddo
      enddo

c
c deallocate workspace
c
      deallocate(bufsend,bufrecv,arrl)
      deallocate(scnts,sdspl,rcnts,rdspl)

      return
      end subroutine xgridremap_ijl

      subroutine init_cs2llint_type(grid_cs,grid_ll, lons,lats,
     &     cs2llint)
      use dd2d_utils, only : dist_grid
      use geom, only : ll2csxy,csxy2ll,shiftwest
      use constant, only : pi,twopi
      type(dist_grid), intent(in) :: grid_cs,grid_ll
      real*8, dimension(grid_ll%npx) :: lons
      real*8, dimension(grid_ll%npy) :: lats
      type(cs2llint_type), intent(out) :: cs2llint
c local vars
      integer :: imlon,jmlat,ilon,jlat,i,j
      integer :: npts,tile
      integer :: nproc,llproc,csproc
      integer :: ierr
      integer, dimension(:), allocatable :: jell
      integer, dimension(:), allocatable :: ilist,jlist
      real*8 :: dxh,xn,x,y,cs_xmin,cs_xmax,cs_ymin,cs_ymax,lonpass,
     &     londum,angdum,corang
      real*8, dimension(:), allocatable :: ix,jy,edgang
      imlon = grid_ll%npx
      jmlat = grid_ll%npy

      cs2llint%isdcs = grid_cs%isd
      cs2llint%iedcs = grid_cs%ied
      cs2llint%jsdcs = grid_cs%jsd
      cs2llint%jedcs = grid_cs%jed

      cs2llint%imlon = grid_ll%npx
      cs2llint%jsdll = grid_ll%jsd
      cs2llint%jedll = grid_ll%jed

      nproc = grid_cs%nproc
      cs2llint%nproc = nproc

      allocate(cs2llint%send_cnts(nproc),cs2llint%send_displs(nproc))
      allocate(cs2llint%recv_cnts(nproc),cs2llint%recv_displs(nproc))

c collect info about the domain decomp of grid_ll
c for now, assume cs pe set same as ll pe set and
c that grid_ll has a j-only decomp
      allocate(jell(nproc))

      call mpi_allgather(grid_ll%je,1,MPI_INTEGER,jell,1,
     &     MPI_INTEGER,MPI_COMM_WORLD,ierr)

c
c Set up the list of ll points that will be handled by this CS PE.
c
      allocate(ix(imlon*jmlat/2), jy(imlon*jmlat/2))
      allocate(ilist(imlon*jmlat/2), jlist(imlon*jmlat/2))
      cs_xmin = -1d0 + 2d0*(grid_cs%is-1)/grid_cs%npx
      cs_xmax = -1d0 + 2d0*(grid_cs%ie  )/grid_cs%npx
      cs_ymin = -1d0 + 2d0*(grid_cs%js-1)/grid_cs%npx
      cs_ymax = -1d0 + 2d0*(grid_cs%je  )/grid_cs%npx
      dxh = 1d0/grid_cs%npx
      xn = 1d0 - dxh
      allocate(edgang(grid_cs%npx))
      do j=1,grid_cs%npx ! compute pseudolatitudes
        y = -1d0 + 2d0*(dble(j)-.5d0)/grid_cs%npx
        call csxy2ll(xn,y,1,londum,edgang(j)) ! use tile 1 as reference
      enddo
      corang = edgang(grid_cs%npx)
      npts = 0
      cs2llint%send_cnts(:) = 0
      llproc = 1
      do jlat=1,jmlat
        if(jlat .gt. jell(llproc)) llproc=llproc+1
        do ilon=1,imlon
          lonpass = lons(ilon) + shiftwest
          if(lonpass.gt.pi) lonpass=lonpass-twopi
          call ll2csxy(lonpass,lats(jlat),x,y,tile)
          if(tile.ne.grid_cs%tile) cycle
          if(x.le.cs_xmin .or. x.gt.cs_xmax) cycle
          if(y.le.cs_ymin .or. y.gt.cs_ymax) cycle
          npts = npts + 1
          cs2llint%send_cnts(llproc) = cs2llint%send_cnts(llproc) + 1
          ix(npts) = .5d0*(1d0+x)*grid_cs%npx+.5d0
          jy(npts) = .5d0*(1d0+y)*grid_cs%npx+.5d0
          ilist(npts) = ilon
          jlist(npts) = jlat
c adjust for grid discontinuity near cube edges using angle coordinate
          if(abs(x) .gt. xn) then ! x-edge.  adjust jy
            call csxy2ll(abs(x),y,1,londum,angdum)
            if(abs(angdum).le.corang) then ! avoid corner triangles
              j = jy(npts)
              if(y.gt.0.) then
                if(angdum.lt.edgang(j  )) j = j-1
              else
                if(angdum.gt.edgang(j+1)) j = j+1
              endif
              jy(npts) = j +
     &           (angdum-edgang(j))/(edgang(j+1)-edgang(j))
            endif
          endif
          if(abs(y) .gt. xn) then ! y-edge. adjust ix
            call csxy2ll(abs(y),x,1,londum,angdum)
            if(abs(angdum).le.corang) then ! avoid corner triangles
              i = ix(npts)
              if(x.gt.0.) then
                if(angdum.lt.edgang(i  )) i = i-1
              else
                if(angdum.gt.edgang(i+1)) i = i+1
              endif
              ix(npts) = i +
     &             (angdum-edgang(i))/(edgang(i+1)-edgang(i))
            endif
          endif
        enddo
      enddo
      deallocate(edgang)
      cs2llint%npts_interp = npts
      allocate(cs2llint%ix(npts),cs2llint%jy(npts))
      cs2llint%ix(1:npts) = ix(1:npts)
      cs2llint%jy(1:npts) = jy(1:npts)
      cs2llint%send_displs(1) = 0
      do llproc=2,nproc
        cs2llint%send_displs(llproc) =
     &       cs2llint%send_displs(llproc-1)+cs2llint%send_cnts(llproc-1)
      enddo

c
c CS PEs inform LL PEs about the interpolation pts they will be sending
c
      npts = imlon*(1+grid_ll%je-grid_ll%js)
      cs2llint%npts_unpack = npts
      call mpi_alltoall(cs2llint%send_cnts,1,MPI_INTEGER,
     &                  cs2llint%recv_cnts,1,MPI_INTEGER,
     &                  MPI_COMM_WORLD,ierr)
      if(sum(cs2llint%recv_cnts) .ne. npts) then
        write(6,*)
     &       'bad counts in init_cs2llint_type for LL PE ',grid_ll%gid
        call stop_model('bad counts in init_cs2llint_type',255)
      endif
      cs2llint%recv_displs(1) = 0
      do csproc=2,nproc
        cs2llint%recv_displs(csproc) =
     &       cs2llint%recv_displs(csproc-1)+cs2llint%recv_cnts(csproc-1)
      enddo
      allocate(cs2llint%i_unpack(npts),cs2llint%j_unpack(npts))
      call mpi_alltoallv(
     &     ilist,             cs2llint%send_cnts,cs2llint%send_displs,
     &          MPI_INTEGER,
     &     cs2llint%i_unpack, cs2llint%recv_cnts,cs2llint%recv_displs,
     &          MPI_INTEGER,
     &     MPI_COMM_WORLD, ierr)
      call mpi_alltoallv(
     &     jlist,             cs2llint%send_cnts,cs2llint%send_displs,
     &          MPI_INTEGER,
     &     cs2llint%j_unpack, cs2llint%recv_cnts,cs2llint%recv_displs,
     &          MPI_INTEGER,
     &     MPI_COMM_WORLD, ierr)

c
c deallocate workspace
c
      deallocate(ix,jy,ilist,jlist,jell)

      return
      end subroutine init_cs2llint_type

      subroutine cs2llint_ij(grid_cs,cs2llint,arrcs,arrll)
      use dd2d_utils, only : dist_grid,halo_update
      type(dist_grid), intent(in) :: grid_cs
      type(cs2llint_type), intent(in) :: cs2llint
      real*8, dimension(cs2llint%isdcs:cs2llint%iedcs,
     &                  cs2llint%jsdcs:cs2llint%jedcs) :: arrcs
c      intent(in) :: arrcs
      real*8, dimension(cs2llint%imlon,cs2llint%jsdll:cs2llint%jedll),
     &     intent(out) :: arrll
c local vars
      integer :: i,j,n,nproc,ierr
      real*8, dimension(:), allocatable :: bufsend,bufrecv
      integer, dimension(:), allocatable :: scnts,sdspl,rcnts,rdspl
      real*8 :: wti,wtj

c
c allocate send/recv buffers
c
      nproc = cs2llint%nproc
      allocate(scnts(nproc),sdspl(nproc),rcnts(nproc),rdspl(nproc))
      scnts(:) = cs2llint%send_cnts
      sdspl(:) = cs2llint%send_displs
      rcnts(:) = cs2llint%recv_cnts
      rdspl(:) = cs2llint%recv_displs

      allocate(bufsend(cs2llint%npts_interp),
     &         bufrecv(cs2llint%npts_unpack))

c
c Interpolate and store in 1D send buffer
c
      call halo_update(grid_cs,arrcs)
      call corner_fill_3D(grid_cs,arrcs,1)
      do n=1,cs2llint%npts_interp
        wti = cs2llint%ix(n)
        wtj = cs2llint%jy(n)
        i = wti
        j = wtj
        wti = wti-i
        wtj = wtj-j
        bufsend(n) = 
     &         wtj*(wti*arrcs(i+1,j+1)+(1.-wti)*arrcs(i,j+1))
     &   +(1.-wtj)*(wti*arrcs(i+1,j  )+(1.-wti)*arrcs(i,j  ))
      enddo

c
c send interpolants to LL PEs
c
      call mpi_alltoallv(bufsend, scnts, sdspl, MPI_DOUBLE_PRECISION,
     &                   bufrecv, rcnts, rdspl, MPI_DOUBLE_PRECISION,
     &                   MPI_COMM_WORLD, ierr)

c
c unpack into latlon array
c
      do n=1,cs2llint%npts_unpack
        arrll(cs2llint%i_unpack(n),cs2llint%j_unpack(n)) = bufrecv(n)
      enddo

c
c deallocate workspace
c
      deallocate(bufsend,bufrecv)
      deallocate(scnts,sdspl,rcnts,rdspl)

      return
      end subroutine cs2llint_ij

      subroutine cs2llint_lij(grid_cs,cs2llint,arrcs,arrll)
      use dd2d_utils, only : dist_grid,halo_update
      type(dist_grid), intent(in) :: grid_cs
      type(cs2llint_type), intent(in) :: cs2llint
      real*8, dimension(:,cs2llint%isdcs:,cs2llint%jsdcs:) :: arrcs
c      intent(in) :: arrcs
      real*8, dimension(:,:,cs2llint%jsdll:), intent(out) :: arrll
c local vars
      integer :: i,j,l,lm,m,n,nproc,ierr
      real*8, dimension(:), allocatable :: bufsend,bufrecv
      integer, dimension(:), allocatable :: scnts,sdspl,rcnts,rdspl
      real*8 :: wti,wtj

      lm = size(arrcs,1)
c
c allocate send/recv buffers
c
      nproc = cs2llint%nproc
      allocate(scnts(nproc),sdspl(nproc),rcnts(nproc),rdspl(nproc))
      scnts(:) = lm*(cs2llint%send_cnts)
      sdspl(:) = lm*(cs2llint%send_displs)
      rcnts(:) = lm*(cs2llint%recv_cnts)
      rdspl(:) = lm*(cs2llint%recv_displs)

      allocate(bufsend(lm*cs2llint%npts_interp),
     &         bufrecv(lm*cs2llint%npts_unpack))

c
c Interpolate and store in 1D send buffer
c
      call halo_update(grid_cs,arrcs,jdim=3)
      call corner_fill_4D(grid_cs,arrcs,lm,1)
      m = 0
      do n=1,cs2llint%npts_interp
        wti = cs2llint%ix(n)
        wtj = cs2llint%jy(n)
        i = wti
        j = wtj
        wti = wti-i
        wtj = wtj-j
        do l=1,lm
          m = m + 1
          bufsend(m) = 
     &           wtj*(wti*arrcs(l,i+1,j+1)+(1.-wti)*arrcs(l,i,j+1))
     &     +(1.-wtj)*(wti*arrcs(l,i+1,j  )+(1.-wti)*arrcs(l,i,j  ))
        enddo
      enddo

c
c send interpolants to LL PEs
c
      call mpi_alltoallv(bufsend, scnts, sdspl, MPI_DOUBLE_PRECISION,
     &                   bufrecv, rcnts, rdspl, MPI_DOUBLE_PRECISION,
     &                   MPI_COMM_WORLD, ierr)

c
c unpack into latlon array
c
      m = 0
      do n=1,cs2llint%npts_unpack
        i = cs2llint%i_unpack(n)
        j = cs2llint%j_unpack(n)
        do l=1,lm
          m = m + 1
          arrll(l,i,j) = bufrecv(m)
        enddo
      enddo

c
c deallocate workspace
c
      deallocate(bufsend,bufrecv)
      deallocate(scnts,sdspl,rcnts,rdspl)

      return
      end subroutine cs2llint_lij

      subroutine init_ll2csint_type(grid_ll,grid_cs,
     &     lons,lats, jlat_min,jlat_max,
     &     iscs,iecs,jscs,jecs,loncs,latcs,
     &     ll2csint)
      use constant, only : pi,twopi
      use dd2d_utils, only : dist_grid
      type(dist_grid), intent(in) :: grid_ll,grid_cs
!@var lons,lats global coordinates of LL grid
      real*8, dimension(grid_ll%npx) :: lons
      real*8, dimension(grid_ll%npy) :: lats
!@var jlat_min,jlat_max min/max valid lats of global grid
      integer, intent(in) :: jlat_min,jlat_max
!@var loncs,latcs coordinate pairs requested by this CS PE
      integer, intent(in) :: iscs,iecs,jscs,jecs
      real*8, dimension(iscs:iecs,jscs:jecs), intent(in) :: loncs,latcs
      type(ll2csint_type), intent(out) :: ll2csint
c local vars
      integer :: i,j,k,imlon,jmlat,n,npts_glob,npts,nllreq_cs,
     &     jmin_search,jmax_search
      integer :: nproc,llproc,csproc
      integer :: ierr
      integer, dimension(:), allocatable :: req_cnts,displs,ilist,jlist,
     &     ilist_loc,jlist_loc,req_i,req_j
      real*8, dimension(:), allocatable :: ix,jy,req_lons,req_lats
      real*8 :: latmin,latmax

      ll2csint%iscs = iscs
      ll2csint%iecs = iecs
      ll2csint%jscs = jscs
      ll2csint%jecs = jecs

      imlon = grid_ll%npx
      jmlat = grid_ll%npy
      ll2csint%imlon = imlon
      ll2csint%jsdll = grid_ll%jsd
      ll2csint%jedll = grid_ll%jed
      ll2csint%jlat_min = jlat_min
      ll2csint%jlat_max = jlat_max

      nproc = grid_cs%nproc
      ll2csint%nproc = nproc

      allocate(ll2csint%send_cnts(nproc),ll2csint%send_displs(nproc))
      allocate(ll2csint%recv_cnts(nproc),ll2csint%recv_displs(nproc))

c
c CS PEs inform LL PEs about the interpolation pts they want
c
      allocate(req_cnts(nproc),displs(nproc))
      nllreq_cs = (1+iecs-iscs)*(1+jecs-jscs)
      ll2csint%npts_unpack = nllreq_cs
      allocate(ilist_loc(nllreq_cs),jlist_loc(nllreq_cs))
      n = 0
      do j=jscs,jecs
        do i=iscs,iecs
          n = n + 1
          ilist_loc(n) = i
          jlist_loc(n) = j
        enddo
      enddo

      call mpi_allgather(nllreq_cs,1,MPI_INTEGER,req_cnts,1,
     &     MPI_INTEGER,MPI_COMM_WORLD,ierr)
      npts_glob = sum(req_cnts)
      displs(1) = 0
      do csproc=2,nproc
        displs(csproc) = displs(csproc-1)+req_cnts(csproc-1)
      enddo
      allocate(req_lons(npts_glob),req_lats(npts_glob))
      allocate(req_i(npts_glob),req_j(npts_glob))
      allocate(ilist(npts_glob),jlist(npts_glob))
      call mpi_allgatherv(loncs,nllreq_cs,MPI_DOUBLE_PRECISION,
     &          req_lons,req_cnts,displs,MPI_DOUBLE_PRECISION,
     &          MPI_COMM_WORLD,ierr)
      call mpi_allgatherv(latcs,nllreq_cs,MPI_DOUBLE_PRECISION,
     &          req_lats,req_cnts,displs,MPI_DOUBLE_PRECISION,
     &          MPI_COMM_WORLD,ierr)
      call mpi_allgatherv(ilist_loc,nllreq_cs,MPI_INTEGER,
     &          req_i,req_cnts,displs,MPI_INTEGER,
     &          MPI_COMM_WORLD,ierr)
      call mpi_allgatherv(jlist_loc,nllreq_cs,MPI_INTEGER,
     &          req_j,req_cnts,displs,MPI_INTEGER,
     &          MPI_COMM_WORLD,ierr)

c
c LL PE checks which requested points lie in its interp domain,
c which ranges from lats(js) (included) to lats(je+1) (not included)
c
      allocate(ix(npts_glob), jy(npts_glob))
      if(grid_ll%have_south_pole) then
        latmin=-pi/2.
      else
        latmin = lats(grid_ll%js)
      endif
      if(grid_ll%have_north_pole) then
        latmax=+pi/2.
      else
        latmax = lats(grid_ll%je+1)
      endif
      jmax_search = min(grid_ll%je,jlat_max)
      jmin_search = max(grid_ll%js,jlat_min)
      n = 0
      npts = 0
      do csproc=1,nproc
        ll2csint%send_cnts(csproc) = 0
        do k=1,req_cnts(csproc)
          n = n + 1
          if(req_lats(n) .lt. latmin) cycle
          if(req_lats(n) .ge. latmax) cycle
          ll2csint%send_cnts(csproc) = ll2csint%send_cnts(csproc) + 1
          npts = npts + 1
          do j=jmax_search,jmin_search,-1
            if(req_lats(n) .ge. lats(j)) exit
          enddo
          if(j.eq.jlat_max) then
            jy(npts) = j + (req_lats(n)-lats(j))/(pi/2.-lats(j))
          elseif(j.lt.jlat_min) then
            jy(npts) = j + (req_lats(n)+pi/2.)/(lats(j+1)+pi/2.)
          else
            jy(npts) = j + (req_lats(n)-lats(j))/(lats(j+1)-lats(j))
          endif
          do i=imlon,1,-1
            if(req_lons(n) .ge. lons(i)) exit
          enddo
          if(i.eq.imlon) then
            ix(npts) = i + (req_lons(n)-lons(i))/(twopi+lons(1)-lons(i))
          elseif(i.eq.0) then
            ix(npts) = i + (twopi+req_lons(n)-lons(imlon))/
     &                  (twopi+lons(1)-lons(imlon))
          else
            ix(npts) = i + (req_lons(n)-lons(i))/(lons(i+1)-lons(i))
          endif
          ilist(npts) = req_i(n)
          jlist(npts) = req_j(n)
        enddo
      enddo
      ll2csint%npts_interp = npts
      allocate(ll2csint%ix(npts),ll2csint%jy(npts))
      ll2csint%ix(1:npts) = ix(1:npts)
      ll2csint%jy(1:npts) = jy(1:npts)
      ll2csint%send_displs(1) = 0
      do csproc=2,nproc
        ll2csint%send_displs(csproc) =
     &       ll2csint%send_displs(csproc-1)+ll2csint%send_cnts(csproc-1)
      enddo

c
c LL PEs inform CS PEs about the interpolation pts they will be sending
c

      call mpi_alltoall(ll2csint%send_cnts,1,MPI_INTEGER,
     &                  ll2csint%recv_cnts,1,MPI_INTEGER,
     &                  MPI_COMM_WORLD,ierr)
      if(sum(ll2csint%recv_cnts) .ne. nllreq_cs) then
        write(6,*)
     &       'bad counts in init_ll2csint_type for CS PE ',grid_cs%gid
        call stop_model('bad counts in init_ll2csint_type',255)
      endif

      ll2csint%recv_displs(1) = 0
      do llproc=2,nproc
        ll2csint%recv_displs(llproc) =
     &       ll2csint%recv_displs(llproc-1)+ll2csint%recv_cnts(llproc-1)
      enddo

      allocate(ll2csint%i_unpack(nllreq_cs),
     &         ll2csint%j_unpack(nllreq_cs))
      call mpi_alltoallv(
     &     ilist,             ll2csint%send_cnts,ll2csint%send_displs,
     &          MPI_INTEGER,
     &     ll2csint%i_unpack, ll2csint%recv_cnts,ll2csint%recv_displs,
     &          MPI_INTEGER,
     &     MPI_COMM_WORLD, ierr)
      call mpi_alltoallv(
     &     jlist,             ll2csint%send_cnts,ll2csint%send_displs,
     &          MPI_INTEGER,
     &     ll2csint%j_unpack, ll2csint%recv_cnts,ll2csint%recv_displs,
     &          MPI_INTEGER,
     &     MPI_COMM_WORLD, ierr)


c
c deallocate workspace
c
      deallocate(req_cnts,displs)
      deallocate(ilist,jlist,ilist_loc,jlist_loc,req_i,req_j)
      deallocate(ix,jy,req_lons,req_lats)

      return
      end subroutine init_ll2csint_type

      subroutine ll2csint_ij(grid_ll,ll2csint,arrll,arrcs)
      use dd2d_utils, only : dist_grid
      use domain_decomp_1d, only : halo_update,north
      type(dist_grid), intent(in) :: grid_ll
      type(ll2csint_type), intent(in) :: ll2csint
      real*8, dimension(grid_ll%npx,grid_ll%jsd:grid_ll%jed)
     &     :: arrll
c      intent(in) :: arrll
      real*8, dimension(ll2csint%iscs:ll2csint%iecs,
     &                  ll2csint%jscs:ll2csint%jecs), 
     &     intent(out) :: arrcs
c local vars
      integer :: i,j,im,n,nproc,ierr
      real*8, dimension(:), allocatable :: bufsend,bufrecv
      integer, dimension(:), allocatable :: scnts,sdspl,rcnts,rdspl
      real*8, dimension(:,:), allocatable :: arr_pad
      real*8 :: wti,wtj

c
c allocate send/recv buffers
c
      nproc = ll2csint%nproc
      allocate(scnts(nproc),sdspl(nproc),rcnts(nproc),rdspl(nproc))
      scnts(:) = ll2csint%send_cnts
      sdspl(:) = ll2csint%send_displs
      rcnts(:) = ll2csint%recv_cnts
      rdspl(:) = ll2csint%recv_displs

      allocate(bufsend(ll2csint%npts_interp),
     &         bufrecv(ll2csint%npts_unpack))
      im = grid_ll%npx
      allocate(arr_pad(0:im+1,grid_ll%jsd:grid_ll%jed))

c
c Interpolate and store in 1D send buffer
c
      call halo_update(grid_ll,arrll,from=north)
      do j=grid_ll%js,min(grid_ll%jed,grid_ll%npy)
        do i=1,im
          arr_pad(i,j) = arrll(i,j)
        enddo
        arr_pad(0,j) = arr_pad(im,j)
        arr_pad(im+1,j) = arr_pad(1,j)
      enddo
      if(grid_ll%have_south_pole) then ! define a single value at the SP
        j = ll2csint%jlat_min
        arr_pad(:,j-1) = sum(arrll(:,j))/im
      endif
      if(grid_ll%have_north_pole) then ! define a single value at the NP
        j = ll2csint%jlat_max
        arr_pad(:,j+1) = sum(arrll(:,j))/im
      endif
      do n=1,ll2csint%npts_interp
        wti = ll2csint%ix(n)
        wtj = ll2csint%jy(n)
        i = wti
        j = wtj
        wti = wti-i
        wtj = wtj-j
        bufsend(n) = 
     &         wtj*(wti*arr_pad(i+1,j+1)+(1.-wti)*arr_pad(i,j+1))
     &   +(1.-wtj)*(wti*arr_pad(i+1,j  )+(1.-wti)*arr_pad(i,j  ))
      enddo

c
c send interpolants to CS PEs
c
      call mpi_alltoallv(bufsend, scnts, sdspl, MPI_DOUBLE_PRECISION,
     &                   bufrecv, rcnts, rdspl, MPI_DOUBLE_PRECISION,
     &                   MPI_COMM_WORLD, ierr)

c
c unpack into CS array
c
      do n=1,ll2csint%npts_unpack
        arrcs(ll2csint%i_unpack(n),ll2csint%j_unpack(n)) = bufrecv(n)
      enddo

c
c deallocate workspace
c
      deallocate(bufsend,bufrecv,arr_pad)
      deallocate(scnts,sdspl,rcnts,rdspl)

      return
      end subroutine ll2csint_ij

      subroutine ll2csint_lij(grid_ll,ll2csint,arrll,arrcs,
     &     is_ll_vector)
      use dd2d_utils, only : dist_grid
      use domain_decomp_1d, only : halo_update_column,north
      type(dist_grid), intent(in) :: grid_ll
      type(ll2csint_type), intent(in) :: ll2csint
      real*8, dimension(:,:,grid_ll%jsd:) :: arrll
c      intent(in) :: arrll
      real*8, dimension(:,ll2csint%iscs:,ll2csint%jscs:), 
     &     intent(out) :: arrcs
      logical, intent(in), optional :: is_ll_vector
c local vars
      integer :: i,j,l,lm,im,m,n,nproc,ierr
      real*8, dimension(:), allocatable :: bufsend,bufrecv
      integer, dimension(:), allocatable :: scnts,sdspl,rcnts,rdspl
      real*8, dimension(:,:,:), allocatable :: arr_pad
      real*8 :: wti,wtj
      logical :: is_ll_vector_

      lm = size(arrll,1)

      is_ll_vector_ = .false.
      if(present(is_ll_vector)) is_ll_vector_=is_ll_vector
      if(is_ll_vector_ .and. lm.ne.2) call stop_model(
     &     'll2csint_lij: is_ll_vector needs size==2',255)

c
c allocate send/recv buffers
c
      nproc = ll2csint%nproc
      allocate(scnts(nproc),sdspl(nproc),rcnts(nproc),rdspl(nproc))

      scnts(:) = lm*(ll2csint%send_cnts)
      sdspl(:) = lm*(ll2csint%send_displs)
      rcnts(:) = lm*(ll2csint%recv_cnts)
      rdspl(:) = lm*(ll2csint%recv_displs)

      allocate(bufsend(lm*ll2csint%npts_interp),
     &         bufrecv(lm*ll2csint%npts_unpack))

      im = grid_ll%npx
      allocate(arr_pad(lm,0:im+1,grid_ll%jsd:grid_ll%jed))

c
c Interpolate and store in 1D send buffer
c
      call halo_update_column(grid_ll,arrll,from=north)
      do j=grid_ll%js,min(grid_ll%jed,grid_ll%npy)
        do i=1,im
          arr_pad(:,i,j) = arrll(:,i,j)
        enddo
        arr_pad(:,0,j) = arr_pad(:,im,j)
        arr_pad(:,im+1,j) = arr_pad(:,1,j)
      enddo
      if(grid_ll%have_south_pole) then ! define a single value at the SP
        j = ll2csint%jlat_min
c        if(is_ll_vector_) then ! todo: need sin(lon),cos(lon)
c        else
        do l=1,lm
          arr_pad(l,:,j-1) = sum(arrll(l,:,j))/im
        enddo
c        endif
      endif
      if(grid_ll%have_north_pole) then ! define a single value at the NP
        j = ll2csint%jlat_max
c        if(is_ll_vector_) then ! todo: need sin(lon),cos(lon)
c        else
        do l=1,lm
          arr_pad(l,:,j+1) = sum(arrll(l,:,j))/im
        enddo
c        endif
      endif
      m = 0
      do n=1,ll2csint%npts_interp
        wti = ll2csint%ix(n)
        wtj = ll2csint%jy(n)
        i = wti
        j = wtj
        wti = wti-i
        wtj = wtj-j
        do l=1,lm
          m = m + 1
          bufsend(m) = 
     &         wtj*(wti*arr_pad(l,i+1,j+1)+(1.-wti)*arr_pad(l,i,j+1))
     &   +(1.-wtj)*(wti*arr_pad(l,i+1,j  )+(1.-wti)*arr_pad(l,i,j  ))
        enddo
      enddo

c
c send interpolants to CS PEs
c
      call mpi_alltoallv(bufsend, scnts, sdspl, MPI_DOUBLE_PRECISION,
     &                   bufrecv, rcnts, rdspl, MPI_DOUBLE_PRECISION,
     &                   MPI_COMM_WORLD, ierr)

c
c unpack into CS array
c
      m = 0
      do n=1,ll2csint%npts_unpack
        i = ll2csint%i_unpack(n)
        j = ll2csint%j_unpack(n)
        do l=1,lm
          m = m + 1
          arrcs(l,i,j) = bufrecv(m)
        enddo
      enddo

c
c deallocate workspace
c
      deallocate(bufsend,bufrecv,arr_pad)
      deallocate(scnts,sdspl,rcnts,rdspl)

      return
      end subroutine ll2csint_lij

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

