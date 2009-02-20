!@sum  DIAG_ZONAL defines the resolution and array bounds for zonal
!@+    diagnostics, provides gather/scatter procedures for
!@+    accumulation arrays, and will soon include other code
!@+    and better documentation.
!@auth Cubed Sphere Development Team
!@ver  1.0
      module diag_zonal
      USE CONSTANT, only : twopi
      USE MODEL_COM, only : im,jm
      USE DOMAIN_DECOMP_ATM, only : dist_grid,sumxpe,am_i_root
      use precision_mod, only : reduce_precision
      implicit none
      private

!@param JM_BUDG grid size for budget page diags
!@param IMLON,JMLAT latlon grid sizes giving approx. equivalent res.
      INTEGER, PARAMETER, public :: JM_BUDG=46,JMLAT=2*JM
      INTEGER, PARAMETER :: IMLON=4*IM,IMH=2*IM
!@var XWON scale factor for diag. printout needed for Wonderland model
      REAL*8, public :: XWON = 1d0

c      public :: imlon,imh

      public :: pack_lc,unpack_lc,get_alloc_bounds

      interface pack_lc
        module procedure pack_lc_2d
        module procedure pack_lc_3d
c        module procedure pack_lc_4d
      end interface pack_lc

      interface unpack_lc
        module procedure unpack_lc_2d
        module procedure unpack_lc_3d
c        module procedure unpack_lc_4d
      end interface unpack_lc

      contains

c      subroutine get_bounds(grid,
c     &     j_strt_budg,j_stop_budg,
c     &     j_strt_jk,j_stop_jk
c     &     )
c      type(dist_grid), intent(in) :: grid
c      integer, optional, intent(out) :: j_strt_budg,j_stop_budg
c      integer, optional, intent(out) :: j_strt_jk,j_stop_jk
c      if(present(j_strt_budg)) j_strt_budg = grid%j_strt
c      if(present(j_stop_budg)) j_stop_budg = grid%j_stop
c      if(present(j_strt_jk)) j_strt_jk = grid%j_strt
c      if(present(j_stop_jk)) j_stop_jk = grid%j_stop
c      return
c      end subroutine get_bounds

      subroutine get_alloc_bounds(grid,
     &     j_strt_budg,j_stop_budg,
     &     j_strt_jk,j_stop_jk
     &     )
      type(dist_grid), intent(in) :: grid
      integer, optional, intent(out) :: j_strt_budg,j_stop_budg
      integer, optional, intent(out) :: j_strt_jk,j_stop_jk
      if(present(j_strt_budg)) j_strt_budg = 1
      if(present(j_stop_budg)) j_stop_budg = jm_budg
      if(present(j_strt_jk)) j_strt_jk = 1
      if(present(j_stop_jk)) j_stop_jk = jmlat
      return
      end subroutine get_alloc_bounds

      subroutine pack_lc_2d(grid,arr_loc,arr_glob)
      type(dist_grid), intent(in) :: grid
      real*8, intent(inout) :: arr_loc(:,:)
      real*8, intent(inout) :: arr_glob(:,:)
      call sumxpe(arr_loc, arr_glob, increment=.true.)
      if(am_i_root()) call reduce_precision(arr_glob,1d-8)
      arr_loc=0
      return
      end subroutine pack_lc_2d
      subroutine pack_lc_3d(grid,arr_loc,arr_glob)
      type(dist_grid), intent(in) :: grid
      real*8, intent(inout) :: arr_loc(:,:,:)
      real*8, intent(inout) :: arr_glob(:,:,:)
      call sumxpe(arr_loc, arr_glob, increment=.true.)
      if(am_i_root()) call reduce_precision(arr_glob,1d-8)
      arr_loc=0
      return
      end subroutine pack_lc_3d
c      subroutine pack_lc_4d(grid,arr_loc,arr_glob)
c      type(dist_grid), intent(in) :: grid
c      real*8, intent(inout) :: arr_loc(:,:,:,:)
c      real*8, intent(inout) :: arr_glob(:,:,:,:)
c      call sumxpe(arr_loc, arr_glob, increment=.true.)
c      if(am_i_root()) call reduce_precision(arr_glob,1d-8)
c      arr_loc=0
c      return
c      end subroutine pack_lc_4d

      subroutine unpack_lc_2d(grid,arr_glob,arr_loc)
      type(dist_grid), intent(in) :: grid
      real*8, intent(in)   :: arr_glob(:,:)
      real*8, intent(out)  :: arr_loc(:,:)
      arr_loc=0
      return
      end subroutine unpack_lc_2d
      subroutine unpack_lc_3d(grid,arr_glob,arr_loc)
      type(dist_grid), intent(in) :: grid
      real*8, intent(in)   :: arr_glob(:,:,:)
      real*8, intent(out)  :: arr_loc(:,:,:)
      arr_loc=0
      return
      end subroutine unpack_lc_3d
      subroutine unpack_lc_4d(grid,arr_glob,arr_loc)
      type(dist_grid), intent(in) :: grid
      real*8, intent(in)   :: arr_glob(:,:,:,:)
      real*8, intent(out)  :: arr_loc(:,:,:,:)
      arr_loc=0
      return
      end subroutine unpack_lc_4d

      end module diag_zonal
