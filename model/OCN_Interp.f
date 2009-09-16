#include "rundeck_opts.h"

!#define AG2OG_PRECIP_BUNDLE
!#define OG2AG_TOC2SST_BUNDLE
!#define AG2OG_OCEANS_BUNDLE
!#define OG2AG_OCEANS_BUNDLE      

!#define DEBUG_OG2AG

!#define BUNDLE_INTERP
!#define BUNDLE_INTERP_NEW

#ifdef BUNDLE_INTERP
#ifdef BUNDLE_INTERP_NEW
      module bundle_arrays
      implicit none
      private

      public lookup_str
      public ba_init, ba_add, ba_bundle, ba_unbundle,
     &     ba_add_new, ba_copy

      integer, parameter :: N_LOOKUP_MAX=256

      type lookup_record
        integer :: lm,km
        real*8, pointer :: src(:), dest(:)
        real*8, pointer :: src_w(:,:), dest_w(:,:)
      end type lookup_record

      type lookup_str
        integer :: si_0, si_1, sj_0, sj_1  ! source bounds
        integer :: di_0, di_1, dj_0, dj_1  ! destination bounds
        integer  :: n_lookup=0
        type (lookup_record) :: lr(N_LOOKUP_MAX)
      end type lookup_str

      interface ba_add
         module procedure ba_add_2
         module procedure ba_add_3
      end interface


      contains

      subroutine ba_init( lstr, si_0, si_1, sj_0, sj_1,
     &     di_0, di_1, dj_0, dj_1 )
      implicit none
      type (lookup_str) :: lstr
      integer :: si_0, si_1, sj_0, sj_1,
     &     di_0, di_1, dj_0, dj_1 
      
      lstr%si_0 = si_0
      lstr%si_1 = si_1
      lstr%sj_0 = sj_0
      lstr%sj_1 = sj_1
      lstr%di_0 = di_0
      lstr%di_1 = di_1
      lstr%dj_0 = dj_0
      lstr%dj_1 = dj_1

      lstr%n_lookup = 0
      end subroutine ba_init


      subroutine ba_add_2( lstr, src, dest, src_w, dest_w )
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid  !remove : for debugging only
      implicit none
      type (lookup_str) :: lstr
      real*8, dimension(:,:), target :: src, dest
      real*8, dimension(:,:), target, optional :: src_w, dest_w

      call ba_add_new( lstr, src, dest, shape(src), 'ij', src_w, dest_w)

      end subroutine ba_add_2


      subroutine ba_add_3( lstr, src, dest, src_w, dest_w )
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid  !remove : for debugging only
      implicit none
      type (lookup_str) :: lstr
      real*8, dimension(:,:,:), target :: src, dest
      real*8, dimension(:,:), target, optional :: src_w, dest_w

      call ba_add_new( lstr, src, dest, shape(src), 'lij', src_w,dest_w)

      end subroutine ba_add_3

      subroutine ba_add_new( lstr, src, dest, shp, flag, src_w, dest_w )
      implicit none
      type (lookup_str) :: lstr
      real*8, dimension(*), target :: src, dest
      integer :: shp(:)
      character*(*) :: flag
      real*8, dimension(:,:), target, optional :: src_w, dest_w

      integer :: si_0, si_1, sj_0, sj_1,
     &     di_0, di_1, dj_0, dj_1 
      integer :: sim, sjm, dim, djm, lm, km

      lstr%n_lookup = lstr%n_lookup + 1
      if ( lstr%n_lookup > N_LOOKUP_MAX )
     &     call stop_model("ba_add: increase N_LOOKUP_MAX", 255)

      print *,"shp=" ,shp

      si_0 = lstr%si_0
      si_1 = lstr%si_1
      sj_0 = lstr%sj_0
      sj_1 = lstr%sj_1
      di_0 = lstr%di_0
      di_1 = lstr%di_1
      dj_0 = lstr%dj_0
      dj_1 = lstr%dj_1

      sim = si_1 - si_0 + 1
      sjm = sj_1 - sj_0 + 1
      dim = di_1 - di_0 + 1
      djm = dj_1 - dj_0 + 1

      select case( flag )
      case('ij')
        lm = 1
        km = 1
      case('lij')
        lm = shp(1)
        km = 1
      case('ijk')
        lm = 1
        km = shp(3)
      case('lijk')
        lm = shp(1)
        km = shp(4)
      case default
        call stop_model("ba_add: unexpected flag", 255)
      end select

      lstr%lr(lstr%n_lookup)%lm = lm
      lstr%lr(lstr%n_lookup)%km = km

      lstr%lr(lstr%n_lookup)%src => src(1:lm*sim*sjm*km)
      lstr%lr(lstr%n_lookup)%dest => dest(1:lm*dim*djm*km)

      if ( present(src_w) .and. present(dest_w) ) then
        lstr%lr(lstr%n_lookup)%src_w => src_w(:,:)
        lstr%lr(lstr%n_lookup)%dest_w => dest_w(:,:)
      else
        if ( present(src_w) ) call stop_model(
     &       "ba_add: use both weights or none", 255)
        nullify( lstr%lr(lstr%n_lookup)%src_w )
        nullify( lstr%lr(lstr%n_lookup)%dest_w )
      endif

      end subroutine ba_add_new

      subroutine ba_bundle( lstr, buf_s, buf_d )
      !USE DOMAIN_DECOMP_ATM, only : agrid=>grid  !remove : for debugging only
      implicit none
      type (lookup_str) :: lstr
      real*8, dimension(:,:,:), pointer :: buf_s, buf_d

      integer :: si_0, si_1, sj_0, sj_1,
     &     di_0, di_1, dj_0, dj_1 
      integer im,jm,km,lm
      integer i,j,k,l,m,n,ind
      

      si_0 = lstr%si_0
      si_1 = lstr%si_1
      sj_0 = lstr%sj_0
      sj_1 = lstr%sj_1
      di_0 = lstr%di_0
      di_1 = lstr%di_1
      dj_0 = lstr%dj_0
      dj_1 = lstr%dj_1

      im = si_1 - si_0 + 1
      jm = sj_1 - sj_0 + 1

      n = 0
      do k=1,lstr%n_lookup
        n = n+lstr%lr(k)%km*lstr%lr(k)%lm
      enddo

      allocate( buf_s(n, si_0:si_1, sj_0:sj_1) )
      allocate( buf_d(n, di_0:di_1, dj_0:dj_1) )

      n = 0
      do m = 1,lstr%n_lookup
        km = lstr%lr(m)%km
        lm = lstr%lr(m)%lm
        do k=0,km-1
          do l=0,lm-1
            n = n+1

            do j=0,jm-1
              do i=0,im-1
                ind = l + i*lm + j*im*lm + k*jm*im*lm + 1
                buf_s(n,i+si_0,j+sj_0) = lstr%lr(m)%src(ind)
              enddo
            enddo

            if ( associated(lstr%lr(m)%src_w) ) then
              do j=0,jm-1
                do i=0,im-1
                  ind = l + i*lm + j*im*lm + k*jm*im*lm + 1
                  buf_s(n,i+si_0,j+sj_0) = buf_s(n,i+si_0,j+sj_0)
     &                 * lstr%lr(m)%src_w(i+1,j+1)
                enddo
              enddo
            endif

          enddo
        enddo
      enddo

      end subroutine ba_bundle


      subroutine ba_unbundle( lstr, buf_s, buf_d )
      implicit none
      type (lookup_str) :: lstr
      real*8, dimension(:,:,:), pointer :: buf_s, buf_d

      integer :: di_0, di_1, dj_0, dj_1
      integer im,jm,km,lm
      integer i,j,k,l,m,n,ind

      di_0 = lstr%di_0
      di_1 = lstr%di_1
      dj_0 = lstr%dj_0
      dj_1 = lstr%dj_1

      im = di_1 - di_0 + 1
      jm = dj_1 - dj_0 + 1

      n = 0
      do m = 1,lstr%n_lookup
        km = lstr%lr(m)%km
        lm = lstr%lr(m)%lm
        do k=0,km-1
          do l=0,lm-1
            n = n+1

            do j=0,jm-1
              do i=0,im-1
                ind = l + i*lm + j*im*lm + k*jm*im*lm + 1
                lstr%lr(m)%dest(ind) = buf_d(n,i+di_0,j+dj_0)
              enddo
            enddo

          enddo
        enddo
      enddo

      n = 0
      do m = 1,lstr%n_lookup
        km = lstr%lr(m)%km
        lm = lstr%lr(m)%lm
        do k=0,km-1
          do l=0,lm-1
            n = n+1
            
            if ( associated(lstr%lr(m)%dest_w ) ) then
              do j=0,jm-1
                do i=0,im-1
                  ind = l + i*lm + j*im*lm + k*jm*im*lm + 1
                  if ( lstr%lr(m)%dest_w(i+1,j+1) .ne. 0.d0 ) then
                    lstr%lr(m)%dest(ind) = lstr%lr(m)%dest(ind)
     &                   / lstr%lr(m)%dest_w(i+1,j+1)
                  else
                    lstr%lr(m)%dest(ind) = 0.d0
                  endif
                enddo
              enddo
            endif

          enddo
        enddo
      enddo

      deallocate( buf_s )
      deallocate( buf_d )     

      end subroutine ba_unbundle

      subroutine ba_copy( lstr )
      implicit none
      type (lookup_str) :: lstr

      integer :: si_0, si_1, sj_0, sj_1,
     &     di_0, di_1, dj_0, dj_1
      integer im,jm,km,lm,nm
      integer m

      si_0 = lstr%si_0
      si_1 = lstr%si_1
      sj_0 = lstr%sj_0
      sj_1 = lstr%sj_1
      di_0 = lstr%di_0
      di_1 = lstr%di_1
      dj_0 = lstr%dj_0
      dj_1 = lstr%dj_1

      im = si_1 - si_0 + 1
      jm = sj_1 - sj_0 + 1

      ! make sure we copy arrays of the same size
      if ( im .ne. di_1 - di_0 + 1 .or. jm .ne. dj_1 - dj_0 + 1 )
     &     call stop_model("ba_copy: dims of arrays differ", 255)

      do m = 1,lstr%n_lookup
        km = lstr%lr(m)%km
        lm = lstr%lr(m)%lm
        nm = lm*im*jm*km
        lstr%lr(m)%dest(1:nm) = lstr%lr(m)%src(1:nm)
      enddo

      end subroutine ba_copy

      end module bundle_arrays

#else /* BUNDLE_INTERP_NEW */

      module bundle_arrays
      implicit none
      private

      public lookup_str
      public ba_init, ba_add, ba_bundle, ba_unbundle

      interface ba_add
         module procedure ba_add_2
         module procedure ba_add_3
      end interface

      integer, parameter :: N_LOOKUP_MAX=256

      type lookup_record
        integer :: rank,km
        real*8, pointer :: src2(:,:), dest2(:,:)
        real*8, pointer :: src_w(:,:), dest_w(:,:)
        real*8, pointer :: src3(:,:,:), dest3(:,:,:)
      end type lookup_record

      type lookup_str
        integer :: si_0, si_1, sj_0, sj_1  ! source bounds
        integer :: di_0, di_1, dj_0, dj_1  ! destination bounds
        integer  :: n_lookup=0
        type (lookup_record) :: lr(N_LOOKUP_MAX)
      end type lookup_str


      contains

      subroutine ba_init( lstr, si_0, si_1, sj_0, sj_1,
     &     di_0, di_1, dj_0, dj_1 )
      implicit none
      type (lookup_str) :: lstr
      integer :: si_0, si_1, sj_0, sj_1,
     &     di_0, di_1, dj_0, dj_1 
      
      lstr%si_0 = si_0
      lstr%si_1 = si_1
      lstr%sj_0 = sj_0
      lstr%sj_1 = sj_1
      lstr%di_0 = di_0
      lstr%di_1 = di_1
      lstr%dj_0 = dj_0
      lstr%dj_1 = dj_1

      lstr%n_lookup = 0
      end subroutine ba_init


      subroutine ba_add_2( lstr, src, dest, src_w, dest_w )
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid  !remove : for debugging only
      implicit none
      type (lookup_str) :: lstr
      real*8, dimension(:,:), target :: src, dest
      real*8, dimension(:,:), target, optional :: src_w, dest_w

      lstr%n_lookup = lstr%n_lookup + 1
      if ( lstr%n_lookup > N_LOOKUP_MAX )
     &     call stop_model("ba_add: increase N_LOOKUP_MAX", 255)

      lstr%lr(lstr%n_lookup)%src2 => src(:,:)
      lstr%lr(lstr%n_lookup)%dest2 => dest(:,:)
      lstr%lr(lstr%n_lookup)%rank = 2
      lstr%lr(lstr%n_lookup)%km = 1

      if ( present(src_w) .and. present(dest_w) ) then
        lstr%lr(lstr%n_lookup)%src_w => src_w(:,:)
        lstr%lr(lstr%n_lookup)%dest_w => dest_w(:,:)
      else
        if ( present(src_w) ) call stop_model(
     &       "ba_add: use both weights or none", 255)
        nullify( lstr%lr(lstr%n_lookup)%src_w )
        nullify( lstr%lr(lstr%n_lookup)%dest_w )
      endif

      if (agrid%gid .eq. 0) write(*,*) "nlookup, KM=",
     & lstr%n_lookup,lstr%lr(lstr%n_lookup)%km

      end subroutine ba_add_2


      subroutine ba_add_3( lstr, src, dest, src_w, dest_w )
      implicit none
      type (lookup_str) :: lstr
      real*8, dimension(:,:,:), target :: src, dest
      real*8, dimension(:,:), target, optional :: src_w, dest_w

      lstr%n_lookup = lstr%n_lookup + 1
      if ( lstr%n_lookup > N_LOOKUP_MAX )
     &     call stop_model("ba_add: increase N_LOOKUP_MAX", 255)

      lstr%lr(lstr%n_lookup)%src3 => src(:,:,:)
      lstr%lr(lstr%n_lookup)%dest3 => dest(:,:,:)
      lstr%lr(lstr%n_lookup)%rank = 3
      lstr%lr(lstr%n_lookup)%km = size(src,1)

      if ( present(src_w) .and. present(dest_w) ) then
        lstr%lr(lstr%n_lookup)%src_w => src_w(:,:)
        lstr%lr(lstr%n_lookup)%dest_w => dest_w(:,:)
      else
        if ( present(src_w) ) call stop_model(
     &       "ba_add: use both weights or none", 255)
        nullify( lstr%lr(lstr%n_lookup)%src_w )
        nullify( lstr%lr(lstr%n_lookup)%dest_w )
      endif

      end subroutine ba_add_3


      subroutine ba_bundle( lstr, buf_s, buf_d )
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid  !remove : for debugging only
      implicit none
      type (lookup_str) :: lstr
      real*8, dimension(:,:,:), pointer :: buf_s, buf_d

      integer :: si_0, si_1, sj_0, sj_1,
     &     di_0, di_1, dj_0, dj_1 
      integer k, l, n

      si_0 = lstr%si_0
      si_1 = lstr%si_1
      sj_0 = lstr%sj_0
      sj_1 = lstr%sj_1
      di_0 = lstr%di_0
      di_1 = lstr%di_1
      dj_0 = lstr%dj_0
      dj_1 = lstr%dj_1

      n = 1
      do k=1,lstr%n_lookup
        n = n+lstr%lr(k)%km
      enddo

      allocate( buf_s(n, si_0:si_1, sj_0:sj_1) )
      allocate( buf_d(n, di_0:di_1, dj_0:dj_1) )

      n = 1
      do k=1,lstr%n_lookup
        select case(lstr%lr(k)%rank)
        case(2)
          buf_s(n,:,:) = lstr%lr(k)%src2(:,:)
          if ( associated(lstr%lr(k)%src_w ) )
     &         buf_s(n,:,:) = buf_s(n,:,:) * lstr%lr(k)%src_w(:,:)
        case(3)
          buf_s(n:n+lstr%lr(k)%km-1,:,:) = lstr%lr(k)%src3(:,:,:)
          if ( associated(lstr%lr(k)%src_w ) ) then
             do l=n,n+lstr%lr(k)%km-1
                buf_s(l,:,:) =
     &               buf_s(l,:,:) * lstr%lr(k)%src_w(:,:)
             enddo
          endif
        end select
        n = n+lstr%lr(k)%km
        if (agrid%gid .eq. 0) write(*,*),"index n=",n
      enddo

      end subroutine ba_bundle


      subroutine ba_unbundle( lstr, buf_s, buf_d )
      implicit none
      type (lookup_str) :: lstr
      real*8, dimension(:,:,:), pointer :: buf_s, buf_d

      integer i,j,k,l,n

      n = 1
      do k=1,lstr%n_lookup
        select case(lstr%lr(k)%rank)
        case(2)
           lstr%lr(k)%dest2(:,:) = buf_d(n,:,:)
        case(3)
           lstr%lr(k)%dest3(:,:,:) = buf_d(n:n+lstr%lr(k)%km-1,:,:)
        end select
      n = n+lstr%lr(k)%km
      enddo
 	 
      do k=1,lstr%n_lookup
         select case(lstr%lr(k)%rank)
         case(2)
            if ( associated(lstr%lr(k)%dest_w ) ) then
               where ( lstr%lr(k)%dest_w(:,:) .ne. 0.d0 )
                  lstr%lr(k)%dest2(:,:) =
     &                 lstr%lr(k)%dest2(:,:) / lstr%lr(k)%dest_w(:,:)
               elsewhere
                  lstr%lr(k)%dest2(:,:) = 0.d0
               endwhere
            endif
         case(3)
         if ( associated(lstr%lr(k)%dest_w ) ) then
           do l=1,size(lstr%lr(k)%dest3,1)
            where ( lstr%lr(k)%dest_w(:,:) .ne. 0.d0 )
              lstr%lr(k)%dest3(l,:,:) =
     &             lstr%lr(k)%dest3(l,:,:) / lstr%lr(k)%dest_w(:,:)
            elsewhere
              lstr%lr(k)%dest3(l,:,:) = 0.d0
            endwhere
           enddo
         endif
      end select
      enddo

      deallocate( buf_s )
      deallocate( buf_d )     

      end subroutine ba_unbundle

      end module bundle_arrays
#endif /* BUNDLE_INTERP_NEW */

      subroutine bundle_interpolation(remap,lstr)
      USE bundle_arrays
      USE cs2ll_utils, only : xgridremap_type,xgridremap_lij
      implicit none
      type (xgridremap_type) :: remap
      type (lookup_str) :: lstr
      real*8, pointer :: buf_s(:,:,:), buf_d(:,:,:)

      call ba_bundle( lstr, buf_s, buf_d )

      call xgridremap_lij( remap, buf_s, buf_d)

      call ba_unbundle( lstr, buf_s, buf_d )

      end subroutine bundle_interpolation

#endif /* BUNDLE_INTERP */

#ifdef AG2OG_PRECIP_BUNDLE
      SUBROUTINE AG2OG_precip

!@sum  INT_AG2OG_precip is for interpolation of precipitation
!!       arrays from atmospheric grid to the ocean grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>im, aJM=>jm
      USE OCEAN,      only : oIM=>im, oJM=>jm

      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only : am_i_root,ocn_pack=>pack_data
      USE OCEANR_DIM, only : ogrid

      USE OCEAN, only : oDXYPO=>DXYPO,oDLATM=>DLATM, OXYP
      Use GEOM,  only : AXYP

#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : NTM
#ifdef TRACERS_WATER
      USE FLUXES,  only : aTRPREC=>TRPREC, aTRUNPSI=>TRUNPSI
      USE OFLUXES, only : oTRPREC, oTRUNPSI
#endif
#endif
      USE SEAICE_COM, only : aRSI=>RSI
      USE FLUXES, only : aPREC=>PREC, aEPREC=>EPREC
     *     , aRUNPSI=>RUNPSI, aSRUNPSI=>SRUNPSI, aERUNPSI=>ERUNPSI

      USE OFLUXES, only : oRSI, oPREC, oEPREC
     *     , oRUNPSI, oSRUNPSI, oERUNPSI
      USE regrid_com, only : remap_A2O
      USE bundle_arrays

      IMPLICIT NONE
      INTEGER N
      REAL*8, allocatable :: aWEIGHT(:,:),oWEIGHT(:,:)
      REAL*8, allocatable :: aWEIGHT1(:,:),oWEIGHT1(:,:)
      REAL*8, allocatable :: o8_glob(:,:)
      REAL*4, allocatable :: o4_glob(:,:)

      character*80 :: name,title
      type (lookup_str) :: lstr
      integer :: i,j,l

      allocate(aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

      !-- need to allocate separate array for different weights
      allocate(aweight1(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight1(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

      !-- initializing "lstr"
      call ba_init( lstr, aGRID%I_STRT_HALO, aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO, aGRID%J_STOP_HALO,
     &     oGRID%I_STRT_HALO, oGRID%I_STOP_HALO,
     &     oGRID%J_STRT_HALO, oGRID%J_STOP_HALO )

      aWEIGHT(:,:) = 1.-aRSI(:,:) !!  open ocean fraction
      call ba_add(lstr, aWEIGHT,oWEIGHT )  ! must be called before calling ba_add(aQuantity,oQuantity,aWeight,oWeight)
      call ba_add(lstr,aPREC,oPREC,aWEIGHT,oWEIGHT )
      call ba_add(lstr,aEPREC,oEPREC,aWEIGHT,oWEIGHT )
#if (defined TRACERS_OCEAN) && (defined TRACERS_WATER)
      DO N=1,NTM
        aTRPREC(N,:,:) = aTRPREC(N,:,:)/AXYP(:,:)
      END DO
      call ba_add(lstr, aTRPREC, oTRPREC, aWEIGHT, oWEIGHT )
#endif

      aWEIGHT1(:,:) = aRSI(:,:)   
!     oWEIGHT1(:,:) = 1.-oWEIGHT(:,:)            ! using the property REGRID(1-RSI)=1-REGRID(RSI)  
      call ba_add(lstr,aWEIGHT1,oWEIGHT1)        ! this line should be removed when previous line is uncommented 
      call ba_add(lstr,aRUNPSI,oRUNPSI,aWEIGHT1,oWEIGHT1)
      call ba_add(lstr,aSRUNPSI,oSRUNPSI,aWEIGHT1,oWEIGHT1)
      call ba_add(lstr,aERUNPSI,oERUNPSI,aWEIGHT1,oWEIGHT1)

#if (defined TRACERS_OCEAN) && (defined TRACERS_WATER)
      call ba_add( lstr,aTRUNPSI,oTRUNPSI,aWEIGHT1,oWEIGHT1)
#endif

      call ba_add( lstr, aRSI, oRSI)


c*   actual interpolation here
      call bundle_interpolation(remap_A2O,lstr)

      oEPREC(:,:) = oEPREC(:,:)*OXYP(:,:)
      oSRUNPSI(:,:) = oSRUNPSI(:,:)*OXYP(:,:)
      oERUNPSI(:,:) = oERUNPSI(:,:)*OXYP(:,:)

#if (defined TRACERS_OCEAN) && (defined TRACERS_WATER)
      DO N=1,NTM
        oTRPREC(N,:,:) = oTRPREC(N,:,:)*OXYP(:,:)
      END DO
      DO N=1,NTM
        oTRUNPSI(N,:,:) = oTRUNPSI(N,:,:)*OXYP(:,:)
      END DO
#endif

      deallocate(aweight,oweight)
      deallocate(aweight1,oweight1)

      RETURN
      END SUBROUTINE AG2OG_precip

#else
      SUBROUTINE AG2OG_precip

!@sum  INT_AG2OG_precip is for interpolation of precipitation
!!       arrays from atmospheric grid to the ocean grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>im, aJM=>jm
      USE OCEAN,      only : oIM=>im, oJM=>jm

      USE DOMAIN_DECOMP_ATM, only : agrid=>grid

      USE OCEAN, only : oDXYPO=>DXYPO,oDLATM=>DLATM, OXYP
      Use GEOM,  only : AXYP

#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : NTM
#ifdef TRACERS_WATER
      USE FLUXES,  only : aTRPREC=>TRPREC, aTRUNPSI=>TRUNPSI
      USE OFLUXES, only : oTRPREC, oTRUNPSI
#endif
#endif
      USE SEAICE_COM, only : aRSI=>RSI
      USE FLUXES, only : aPREC=>PREC, aEPREC=>EPREC
     *     , aRUNPSI=>RUNPSI, aSRUNPSI=>SRUNPSI, aERUNPSI=>ERUNPSI

      USE OFLUXES, only : oRSI, oPREC, oEPREC
     *     , oRUNPSI, oSRUNPSI, oERUNPSI

      USE INT_AG2OG_MOD, only : INT_AG2OG

      IMPLICIT NONE

      INTEGER N

      REAL*8, allocatable :: aWEIGHT(:,:)

      character*80 :: name

      allocate(aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))

      aWEIGHT(:,:) = 1.- aRSI(:,:)  !!  open ocean fraction
      CALL INT_AG2OG(aPREC,oPREC, aWEIGHT)

      CALL INT_AG2OG(aEPREC,oEPREC, aWEIGHT)
      oEPREC(:,:) = oEPREC(:,:)*OXYP(:,:)

      aWEIGHT(:,:) = aRSI(:,:)
      CALL INT_AG2OG(aRUNPSI,oRUNPSI, aWEIGHT)

      CALL INT_AG2OG(aSRUNPSI,oSRUNPSI, aWEIGHT)
      oSRUNPSI(:,:) = oSRUNPSI(:,:)*OXYP(:,:)

      CALL INT_AG2OG(aERUNPSI,oERUNPSI, aWEIGHT)
      oERUNPSI(:,:) = oERUNPSI(:,:)*OXYP(:,:)

      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(aRSI,oRSI, aWEIGHT)

#if (defined TRACERS_OCEAN) && (defined TRACERS_WATER)

      aWEIGHT(:,:) = 1.- aRSI(:,:)
      DO N=1,NTM
        aTRPREC(N,:,:) = aTRPREC(N,:,:)/AXYP(:,:)
      END DO
      CALL INT_AG2OG(aTRPREC,oTRPREC, aWEIGHT, NTM)
      DO N=1,NTM
        oTRPREC(N,:,:) = oTRPREC(N,:,:)*OXYP(:,:)
      END DO

      aWEIGHT(:,:) = aRSI(:,:)
      CALL INT_AG2OG(aTRUNPSI,oTRUNPSI, aWEIGHT, NTM)
      DO N=1,NTM
        oTRUNPSI(N,:,:) = oTRUNPSI(N,:,:)*OXYP(:,:)
      END DO
#endif

      deallocate(aweight)

      RETURN
      END SUBROUTINE AG2OG_precip
#endif /* ag2og_precip_bundle */

#ifdef OG2AG_TOC2SST_BUNDLE
      SUBROUTINE OG2AG_TOC2SST
!@sum  OG2AG_oceans: ocean arrays used in the subr. TOC2SST are gathered
!       on the ocean grid, interpolated to the atm. grid, and scattered
!!      on the atm. grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>IM, aJM=>JM

      USE OCEAN, only : oIM=>IM, oJM=>JM, oLM=>LMO, oFOCEAN=>FOCEAN
     *                , oDXYPO=>DXYPO, OXYP, oIMAXJ=>IMAXJ
     *                , oCOSI=>COSIC,oSINI=>SINIC
     *                , IVSPO=>IVSP,IVNPO=>IVNP
     *                , sinpo, sinvo
#ifndef CUBE_GRID
      Use GEOM,  only : aCOSI=>COSIP,aSINI=>SINIP
#endif
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid, am_i_root
      use dd2d_utils, only : atm_pack=>pack_data
      USE DOMAIN_DECOMP_1D, only : OCN_UNPACK=>UNPACK_DATA
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE,SOUTH
      USE OCEANR_DIM, only : ogrid

      USE OCEAN, only : MO, UO,VO, G0M
     *     , S0M, OGEOZ,OGEOZ_SV
#ifdef TRACERS_OCEAN
     *     , TRMO
#endif
      USE AFLUXES, only : aMO, aUO1,aVO1, aG0
     *     , aS0, aOGEOZ,aOGEOZ_SV
#ifdef TRACERS_OCEAN
     *     , aTRAC
#endif
      USE OFLUXES, only : oRSI

#ifdef TRACERS_OCEAN
      USE TRACER_COM, only: NTM
      USE OCN_TRACER_COM, only: conc_from_fw
#endif
#ifdef TRACERS_OceanBiology
!only for TRACERS_OceanBiology and not for seawifs
!/bc we interpolate an internal field
      USE obio_com, only: tot_chlo
      USE FLUXES, only : CHL
#endif
#ifdef TRACERS_GASEXCH_ocean_CO2
      USE MODEL_COM, only: nstep=>itime
      USE obio_com, only: pCO2
      USE TRACER_COM, only : vol2mass
      USE FLUXES, only : gtracer
#endif
      USE regrid_com, only : remap_O2A
      USE bundle_arrays
      USE INT_OG2AG_MOD, only : INT_OG2AG

      IMPLICIT NONE
      INTEGER N
      INTEGER IER, I,J,K,L, NT
      INTEGER oJ_0,oJ_1, oI_0,oI_1, oJ_0S,oJ_1S
      REAL*8, allocatable :: oFOCEAN_loc(:,:)
      REAL*8 :: UNP,VNP,AWT1,AWT2
      REAL*8, ALLOCATABLE :: oG0(:,:,:), oS0(:,:,:)
     *                     , oUO1(:,:), oVO1(:,:), oTRAC(:,:,:)
     *                     , oTOT_CHLO_loc(:,:),opCO2_loc(:,:)
     *                     , oMOtmp(:,:,:),OGEOZtmp(:,:)
     *                     , OGEOZ_SVtmp(:,:)
      REAL*8, allocatable :: aWEIGHT(:,:),oWEIGHT(:,:)
      REAL*8, allocatable :: aWEIGHT1(:,:),oWEIGHT1(:,:)
      REAL*8, allocatable :: aones(:,:,:),oones(:,:,:),
     &    aones_glob(:,:,:,:)
      REAL*8, allocatable :: atwos(:,:,:),otwos(:,:,:),
     &    atwos_glob(:,:,:)
      REAL*8, allocatable :: athrees(:,:),othrees(:,:),
     &    athrees_glob(:,:,:) 
      REAL*4, allocatable :: a4_glob(:,:,:)
      character*80 :: name,title
      type (lookup_str) :: lstr

      allocate(aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

      !-- need to allocate separate array for different weights
      allocate(aweight1(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight1(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

      write(*,*) "OG2AG_TOC2SST_BUNDLE"

      !-- initializing "lstr"
      call ba_init( lstr, 1, oIM,
     &     oGRID%J_STRT_HALO, oGRID%J_STOP_HALO,
     &     aGRID%I_STRT_HALO, aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO, aGRID%J_STOP_HALO)

      oI_0 = oGRID%I_STRT
      oI_1 = oGRID%I_STOP
      oJ_0 = oGRID%j_STRT
      oJ_1 = oGRID%j_STOP
      oJ_0S = oGRID%j_STRT_SKP
      oJ_1S = oGRID%j_STOP_SKP

      ALLOCATE
     *  (oG0(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,2), STAT = IER)
      ALLOCATE
     *  (oS0(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,2), STAT = IER)
      ALLOCATE
     *  (oUO1(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO), STAT = IER)
      ALLOCATE
     *  (oVO1(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO), STAT = IER)
#ifdef TRACERS_OCEAN
      ALLOCATE
     *  (oTRAC(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,NTM), STAT = IER)
#endif
      ALLOCATE
     *  (oTOT_CHLO_loc(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
     * , STAT = IER)
      ALLOCATE
     *  (opCO2_loc(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) ,STAT = IER)

      allocate(oFOCEAN_loc(oIM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

      allocate(oMOtmp(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,2),
     &     OGEOZtmp(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
     &     OGEOZ_SVtmp(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))


      CALL OCN_UNPACK (oGRID,oFOCEAN,oFOCEAN_loc)


#ifdef DEBUG_OG2AG
      allocate(aones(4,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oones(4,oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(atwos(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO,5))
      allocate(otwos(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,5))
      allocate(athrees(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(othrees(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))


      do l=1,4
      oones(l,:,:)=5.*l
      enddo
      do k=1,5
      otwos(:,:,k)=4.*k+agrid%gid
      enddo
      othrees=3.+agrid%gid

      oweight= oFOCEAN_loc(:,:)
      oweight1= oFOCEAN_loc(:,:)
c      CALL INT_OG2AG(otwos,atwos,oWEIGHT,5,5,.FALSE.)
      call ba_add(lstr,oWEIGHT,aWEIGHT)
      call ba_add_new(lstr,oones,aones,shape(oones),'lij',
     &   oWEIGHT,aWEIGHT )
      call ba_add_new(lstr,otwos,atwos,shape(otwos),'ijk',
     &   oWEIGHT,aWEIGHT )
      call ba_add(lstr,oWEIGHT1,aWEIGHT1)
      call ba_add(lstr,othrees,athrees,oWEIGHT1,aWEIGHT1 )
#else /* not DEBUG_OG2AG */

      oWEIGHT(:,:) = oFOCEAN_loc(:,:)
      call ba_add( lstr, oWEIGHT, aWEIGHT)

      do L=1,2
      oMOtmp(:,:,L) = MO(:,:,L)
      if (ogrid%HAVE_NORTH_POLE) 
     &     oMOtmp(2:oIM,oJM,L) = oMOtmp(1,oJM,L)
      enddo
      call ba_add_new( lstr, oMOtmp, aMO, shape(oMOtmp), 
     &     'ijk', oWEIGHT, aWEIGHT) 
c      CALL INT_OG2AG(MO,aMO,oWEIGHT,oLM,2,.FALSE.)

      OGEOZtmp=OGEOZ
      if (ogrid%HAVE_NORTH_POLE) 
     &     OGEOZtmp(2:oIM,oJM) = OGEOZtmp(1,oJM)
      call ba_add( lstr, OGEOZtmp, aOGEOZ, oWEIGHT, aWEIGHT)   
c      CALL INT_OG2AG(OGEOZ,aOGEOZ, oWEIGHT, .FALSE.)

      OGEOZ_SVtmp=OGEOZ_SV
      if (ogrid%HAVE_NORTH_POLE) 
     &     OGEOZ_SVtmp(2:oIM,oJM) = OGEOZ_SVtmp(1,oJM)
      call ba_add( lstr, OGEOZ_SVtmp,aOGEOZ_SV, oWEIGHT, aWEIGHT)
      CALL INT_OG2AG(OGEOZ_SV,aOGEOZ_SV, oWEIGHT, .FALSE.)

      oWEIGHT1(:,:) = MO(:,:,1)*oFOCEAN_loc(:,:)
      if (ogrid%HAVE_NORTH_POLE) 
     &     oWEIGHT1(2:oIM,oJM) = oWEIGHT1(1,oJM)
      call ba_add( lstr, oWEIGHT1, aWEIGHT1)

      oG0(:,:,:) = 0.d0
      DO L = 1,2
        DO J=oJ_0,oJ_1
          DO I=oI_0,oIMAXJ(J)
            IF (oFOCEAN_loc(I,J).gt.0.) THEN
              oG0(I,J,L) = G0M(I,J,L)/(MO(I,J,L)*OXYP(I,J))
            END IF
          END DO
        END DO
      END DO
      do L=1,2
      if (ogrid%HAVE_NORTH_POLE) oG0(2:oIM,oJM,L) = oG0(1,oJM,L)
      enddo
      call ba_add_new( lstr, oG0, aG0, shape(oG0),  'ijk', 
     &     oWEIGHT1, aWEIGHT1) 
c      CALL INT_OG2AG(oG0,aG0, oWEIGHT1, 2,2,.TRUE.)

      oS0(:,:,:) = 0.d0
      DO L = 1,2
        DO J=oJ_0,oJ_1
          DO I=oI_0,oIMAXJ(J)
            IF (oFOCEAN_loc(I,J).gt.0.) THEN
              oS0(I,J,L) = S0M(I,J,L)/(MO(I,J,L)*OXYP(I,J))
            END IF
          END DO
        END DO
      END DO
      do L=1,2
      if (ogrid%HAVE_NORTH_POLE) oS0(2:oIM,oJM,L) = oS0(1,oJM,L)
      enddo
      call ba_add_new( lstr, oS0, aS0, shape(oS0), 'ijk',
     &     oWEIGHT1, aWEIGHT1) 
c      CALL INT_OG2AG(oS0,aS0, oWEIGHT1, 2,2,.TRUE.)

c Discontinued method for ocean C-grid -> atm A-grid:
c use a variant of INT_OG2AG aware of C-grid staggering.
c May be reinstated later when cubed-sphere INT_OG2AG
c has this capability.
c      oUO1(:,:) = UO(:,:,1)
c      oVO1(:,:) = VO(:,:,1)
c      oWEIGHT(:,:) = 1.d0
c      CALL INT_OG2AG(oUO1,oVO1,aUO1,aVO1, oWEIGHT
c     *             , IVSPO,IVNPO)
c
c ocean C-grid -> atm A-grid method requiring fewer INT_OG2AG variants:
c ocean C -> ocean A followed by ocean A -> atm A via INT_OG2AG
c
      call halo_update(ogrid,vo(:,:,1),from=south)
      do j=oJ_0S,oJ_1S
c area weights that would have been used by HNTRP for ocean C -> ocean A
        awt1 = (sinpo(j)-sinvo(j-1))/(sinvo(j)-sinvo(j-1))
        awt2 = 1.-awt1
        i=1
          oUO1(i,j) = .5*(UO(i,j,1)+UO(oIM,j,1))
          oVO1(i,j) = VO(i,j-1,1)*awt1+VO(i,j,1)*awt2
        do i=2,oIM
          oUO1(i,j) = .5*(UO(i,j,1)+UO(i-1,j,1))
          oVO1(i,j) = VO(i,j-1,1)*awt1+VO(i,j,1)*awt2
        enddo
      enddo
      if(oGRID%have_south_pole) then
        oUO1(:,1) = 0.; oVO1(:,1) = 0.
      endif
      if(oGRID%have_north_pole) then ! NP U,V from prognostic polar U,V
        oUO1(:,oJM) = UO(oIM,oJM,1)*oCOSI(:) + UO(IVNPO,oJM,1)*oSINI(:)
! oVO1 currently has no effect when atm is lat-lon
        oVO1(:,oJM) = UO(IVNPO,oJM,1)*oCOSI(:) - UO(oIM,oJM,1)*oSINI(:)
      endif

      if (ogrid%HAVE_NORTH_POLE) oUO1(2:oIM,oJM) = oUO1(1,oJM)
      if (ogrid%HAVE_NORTH_POLE) oVO1(2:oIM,oJM) = oVO1(1,oJM)

      call ba_add(lstr, oUO1, aUO1)
      call ba_add(lstr, oVO1, aVO1)

#ifndef CUBE_GRID   /* this portion will disappear if bundling is used for CS only */
      if(aGRID%have_north_pole) then ! latlon atm needs single polar vector
        UNP = SUM(aUO1(:,aJM)*aCOSI(:))*2/aIM
        VNP = SUM(aUO1(:,aJM)*aSINI(:))*2/aIM
        aUO1(1,aJM) = UNP
        aVO1(1,aJM) = VNP
      endif
#endif

#ifdef TRACERS_OCEAN
C**** surface tracer concentration
      oWEIGHT(:,:) = MO(:,:,1)*oFOCEAN_loc(:,:)
      DO NT = 1,NTM
        if (conc_from_fw(nt)) then  ! define conc from fresh water
        DO J=oJ_0,oJ_1
          DO I=oI_0,oIMAXJ(J)
            IF (oFOCEAN_loc(I,J).gt.0.) THEN
              oTRAC(I,J,NT)=TRMO(I,J,1,NT)/(MO(I,J,1)*OXYP(I,J)
     *             -S0M(I,J,1))
            ELSE
              oTRAC(I,J,NT)=0.
            END IF
          END DO
        END DO
        else  ! define conc from total sea water mass
        DO J=oJ_0,oJ_1
          DO I=oI_0,oIMAXJ(J)
            IF (oFOCEAN_loc(I,J).gt.0.) THEN
              oTRAC(I,J,NT)=TRMO(I,J,1,NT)/(MO(I,J,1)*OXYP(I,J))
            ELSE
              oTRAC(I,J,NT)=0.
            END IF
          END DO
        END DO
        end if
      END DO
      CALL INT_OG2AG(oTRAC,aTRAC, oWEIGHT, NTM,NTM,.TRUE.)

#ifdef TRACERS_OceanBiology
!total ocean chlorophyll. Units are kg,chlorophyll/m3 of seawater
!tot_chlo is defined over all ocean points. Here only use open water
!chorophyll, because that is what is seen by radiation
      oWEIGHT(:,:) = oFOCEAN_loc(:,:) * (1.d0 - oRSI(:,:))
      DO J=oJ_0,oJ_1
        DO I=oI_0,oIMAXJ(J)
          IF (oFOCEAN_loc(I,J).gt.0.) THEN
            oTOT_CHLO_loc(I,J) = tot_chlo(I,J)
          ELSE
            oTOT_CHLO_loc(I,J)=0.
          END IF
        END DO
      END DO
      CALL INT_OG2AG(oTOT_CHLO_loc,CHL, oWEIGHT, .FALSE.)
#endif

#ifdef TRACERS_GASEXCH_ocean_CO2
!partial CO2 pressure in seawater. Units are uatm.
!defined only over open ocean cells, because this is what is
!involved in gas exchage with the atmosphere.
      oWEIGHT(:,:) = oFOCEAN_loc(:,:)*(1.d0-oRSI(:,:))
      DO J=oJ_0,oJ_1
        DO I=oI_0,oIMAXJ(J)
          IF (oFOCEAN_loc(I,J).gt.0.) THEN
            opCO2_loc(I,J) = pCO2(I,J)
          ELSE
            opCO2_loc(I,J)=0.
          END IF
        END DO
      END DO
      DO NT = 1,NTM
         CALL INT_OG2AG(opCO2_loc,aTRAC(:,:,NT), oWEIGHT, .FALSE.)

         !opco2 is in uatm, convert to kg,CO2/kg,air
         aTRAC(:,:,NT) = aTRAC(:,:,NT) * vol2mass(nt)* 1.d-6 ! ppmv (uatm) -> kg,CO2/kg,air

         !gtracer is first set in TRACER_DRV, then atrac is interpolated
         !here from pco2 in the ocean and later OCNDYN sets gtracer=atrac
         !Therefore in timesetep 0 (cold start) pco2 has not yet been defined
         !and atrac has to be hard coded in here, so that we do not have
         !urealistic tracer flux at the air-sea interface during step0.

         if (nstep.eq.0) aTRAC(:,:,NT) = gtracer(NT,1,:,:)

      END DO
#endif
#endif

#endif /* DEBUG_OG2AG */

c*    actual interpolation here
      call bundle_interpolation(remap_O2A,lstr)


      DEALLOCATE(oG0, oS0, oUO1,oVO1, oTOT_CHLO_loc, opCO2_loc)
#ifdef TRACERS_OCEAN
      DEALLOCATE(oTRAC)
#endif

#ifdef DEBUG_OG2AG 
c      allocate(aones_glob(aIM,aJM,5,6))
      allocate(atwos_glob(aIM,aJM,6))
c      allocate(athrees_glob(aIM,aJM,6))
      allocate(a4_glob(aIM,aJM,6))
c      call atm_pack(agrid,aones,aones_glob)
      call atm_pack(agrid,atwos(:,:,5),atwos_glob)
c      call atm_pack(agrid,athrees,athrees_glob)

      if (am_i_root()) then
        write(2500,200) atwos_glob
        open(2300,FILE="O2A",FORM='unformatted',
     &       STATUS='unknown')
        title="test"
        a4_glob=atwos_glob
        write(2300) title,a4_glob
        close(2300)
        open(2600,FILE="O2Aa8",FORM='unformatted',
     &       STATUS='unknown')
        title="testa8"
        write(2600) title,atwos_glob
        close(2600)

      endif
 200        format(f25.22)
      deallocate(aones,oones,atwos,otwos,athrees,othrees)
      deallocate(atwos_glob,a4_glob)

      call stop_model("debug bundle",255)

#endif /* DEBUG_OG2AG */


      deallocate(oweight,aweight,
     &     oweight1,aweight1,
     &     oMOtmp,OGEOZtmp,OGEOZ_SVtmp,
     &     ofocean_loc)



      RETURN
      END SUBROUTINE OG2AG_TOC2SST
#else

      SUBROUTINE OG2AG_TOC2SST
!@sum  OG2AG_oceans: ocean arrays used in the subr. TOC2SST are gathered
!       on the ocean grid, interpolated to the atm. grid, and scattered
!!      on the atm. grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>IM, aJM=>JM

      USE OCEAN, only : oIM=>IM, oJM=>JM, oLM=>LMO, oFOCEAN=>FOCEAN
     *                , oDXYPO=>DXYPO, OXYP, oIMAXJ=>IMAXJ
     *                , oCOSI=>COSIC,oSINI=>SINIC
     *                , IVSPO=>IVSP,IVNPO=>IVNP
     *                , sinpo, sinvo
#ifndef CUBE_GRID
      Use GEOM,  only : aCOSI=>COSIP,aSINI=>SINIP
#endif
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only : OCN_UNPACK=>UNPACK_DATA
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE,SOUTH
      USE OCEANR_DIM, only : ogrid

      USE OCEAN, only : MO, UO,VO, G0M
     *     , S0M, OGEOZ,OGEOZ_SV
#ifdef TRACERS_OCEAN
     *     , TRMO
#endif
      USE AFLUXES, only : aMO, aUO1,aVO1, aG0
     *     , aS0, aOGEOZ,aOGEOZ_SV
#ifdef TRACERS_OCEAN
     *     , aTRAC
#endif
      USE OFLUXES, only : oRSI

#ifdef TRACERS_OCEAN
      USE TRACER_COM, only: NTM
      USE OCN_TRACER_COM, only: conc_from_fw
#endif
#ifdef TRACERS_OceanBiology
!only for TRACERS_OceanBiology and not for seawifs
!/bc we interpolate an internal field
      USE obio_com, only: tot_chlo
      USE FLUXES, only : CHL
#endif
#ifdef TRACERS_GASEXCH_ocean_CO2
      USE MODEL_COM, only: nstep=>itime
      USE obio_com, only: pCO2
      USE TRACER_COM, only : vol2mass
      USE FLUXES, only : gtracer
#endif

      USE INT_OG2AG_MOD, only : INT_OG2AG

      IMPLICIT NONE

      INTEGER IER, I,J, L, NT
      INTEGER oJ_0,oJ_1, oI_0,oI_1, oJ_0S,oJ_1S
      REAL*8, allocatable ::
     * oWEIGHT(:,:), oFOCEAN_loc(:,:)
      REAL*8 :: UNP,VNP,AWT1,AWT2
      REAL*8, ALLOCATABLE :: oG0(:,:,:), oS0(:,:,:)
     *                     , oUO1(:,:), oVO1(:,:), oTRAC(:,:,:)
     *                     , oTOT_CHLO_loc(:,:),opCO2_loc(:,:)
     *                     ,atest(:,:)
      character*80 :: name
      oI_0 = oGRID%I_STRT
      oI_1 = oGRID%I_STOP
      oJ_0 = oGRID%j_STRT
      oJ_1 = oGRID%j_STOP
      oJ_0S = oGRID%j_STRT_SKP
      oJ_1S = oGRID%j_STOP_SKP


c      write(*,*) "OG2AG_TOC2SST old version"

      ALLOCATE
     *  (oG0(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,2), STAT = IER)
      ALLOCATE
     *  (oS0(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,2), STAT = IER)
      ALLOCATE
     *  (oUO1(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO), STAT = IER)
      ALLOCATE
     *  (oVO1(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO), STAT = IER)
#ifdef TRACERS_OCEAN
      ALLOCATE
     *  (oTRAC(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,NTM), STAT = IER)
#endif
      ALLOCATE
     *  (oTOT_CHLO_loc(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
     * , STAT = IER)
      ALLOCATE
     *  (opCO2_loc(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) ,STAT = IER)

      allocate(oFOCEAN_loc(oIM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
     &         oWEIGHT(oIM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) )

      CALL OCN_UNPACK (oGRID,oFOCEAN,oFOCEAN_loc)

      oWEIGHT(:,:) = oFOCEAN_loc(:,:)
      CALL INT_OG2AG(MO,aMO,oWEIGHT,oLM,2,.FALSE.)

      oG0(:,:,:) = 0.d0
      oWEIGHT(:,:) = MO(:,:,1)*oFOCEAN_loc(:,:)
      DO L = 1,2
        DO J=oJ_0,oJ_1
          DO I=oI_0,oIMAXJ(J)
            IF (oFOCEAN_loc(I,J).gt.0.) THEN
              oG0(I,J,L) = G0M(I,J,L)/(MO(I,J,L)*OXYP(I,J))
            END IF
          END DO
        END DO
      END DO
      CALL INT_OG2AG(oG0,aG0, oWEIGHT, 2,2,.TRUE.)

      oS0(:,:,:) = 0.d0
      DO L = 1,2
        DO J=oJ_0,oJ_1
          DO I=oI_0,oIMAXJ(J)
            IF (oFOCEAN_loc(I,J).gt.0.) THEN
              oS0(I,J,L) = S0M(I,J,L)/(MO(I,J,L)*OXYP(I,J))
            END IF
          END DO
        END DO
      END DO
      CALL INT_OG2AG(oS0,aS0, oWEIGHT, 2,2,.TRUE.)

      oWEIGHT(:,:) = oFOCEAN_loc(:,:)
      CALL INT_OG2AG(OGEOZ,aOGEOZ, oWEIGHT, .FALSE.)

      CALL INT_OG2AG(OGEOZ_SV,aOGEOZ_SV, oWEIGHT, .FALSE.)

c Discontinued method for ocean C-grid -> atm A-grid:
c use a variant of INT_OG2AG aware of C-grid staggering.
c May be reinstated later when cubed-sphere INT_OG2AG
c has this capability.
c      oUO1(:,:) = UO(:,:,1)
c      oVO1(:,:) = VO(:,:,1)
c      oWEIGHT(:,:) = 1.d0
c      CALL INT_OG2AG(oUO1,oVO1,aUO1,aVO1, oWEIGHT
c     *             , IVSPO,IVNPO)
c
c ocean C-grid -> atm A-grid method requiring fewer INT_OG2AG variants:
c ocean C -> ocean A followed by ocean A -> atm A via INT_OG2AG
c
      call halo_update(ogrid,vo(:,:,1),from=south)
      do j=oJ_0S,oJ_1S
c area weights that would have been used by HNTRP for ocean C -> ocean A
        awt1 = (sinpo(j)-sinvo(j-1))/(sinvo(j)-sinvo(j-1))
        awt2 = 1.-awt1
        i=1
          oUO1(i,j) = .5*(UO(i,j,1)+UO(oIM,j,1))
          oVO1(i,j) = VO(i,j-1,1)*awt1+VO(i,j,1)*awt2
        do i=2,oIM
          oUO1(i,j) = .5*(UO(i,j,1)+UO(i-1,j,1))
          oVO1(i,j) = VO(i,j-1,1)*awt1+VO(i,j,1)*awt2
        enddo
      enddo
      if(oGRID%have_south_pole) then
        oUO1(:,1) = 0.; oVO1(:,1) = 0.
      endif
      if(oGRID%have_north_pole) then ! NP U,V from prognostic polar U,V
        oUO1(:,oJM) = UO(oIM,oJM,1)*oCOSI(:) + UO(IVNPO,oJM,1)*oSINI(:)
! oVO1 currently has no effect when atm is lat-lon
        oVO1(:,oJM) = UO(IVNPO,oJM,1)*oCOSI(:) - UO(oIM,oJM,1)*oSINI(:)
      endif
      oWEIGHT(:,:) = 1.d0
      CALL INT_OG2AG(oUO1, aUO1, oWEIGHT, .FALSE., AvgPole=.FALSE.)
      CALL INT_OG2AG(oVO1, aVO1, oWEIGHT, .FALSE., AvgPole=.FALSE.)
#ifndef CUBE_GRID
      if(aGRID%have_north_pole) then ! latlon atm needs single polar vector
        UNP = SUM(aUO1(:,aJM)*aCOSI(:))*2/aIM
        VNP = SUM(aUO1(:,aJM)*aSINI(:))*2/aIM
        aUO1(1,aJM) = UNP
        aVO1(1,aJM) = VNP
      endif
#endif

#ifdef TRACERS_OCEAN
C**** surface tracer concentration
      oWEIGHT(:,:) = MO(:,:,1)*oFOCEAN_loc(:,:)
      DO NT = 1,NTM
        if (conc_from_fw(nt)) then  ! define conc from fresh water
        DO J=oJ_0,oJ_1
          DO I=oI_0,oIMAXJ(J)
            IF (oFOCEAN_loc(I,J).gt.0.) THEN
              oTRAC(I,J,NT)=TRMO(I,J,1,NT)/(MO(I,J,1)*OXYP(I,J)
     *             -S0M(I,J,1))
            ELSE
              oTRAC(I,J,NT)=0.
            END IF
          END DO
        END DO
        else  ! define conc from total sea water mass
        DO J=oJ_0,oJ_1
          DO I=oI_0,oIMAXJ(J)
            IF (oFOCEAN_loc(I,J).gt.0.) THEN
              oTRAC(I,J,NT)=TRMO(I,J,1,NT)/(MO(I,J,1)*OXYP(I,J))
            ELSE
              oTRAC(I,J,NT)=0.
            END IF
          END DO
        END DO
        end if
      END DO
      CALL INT_OG2AG(oTRAC,aTRAC, oWEIGHT, NTM,NTM,.TRUE.)

#ifdef TRACERS_OceanBiology
!total ocean chlorophyll. Units are kg,chlorophyll/m3 of seawater
!tot_chlo is defined over all ocean points. Here only use open water
!chorophyll, because that is what is seen by radiation
      oWEIGHT(:,:) = oFOCEAN_loc(:,:) * (1.d0 - oRSI(:,:))
      DO J=oJ_0,oJ_1
        DO I=oI_0,oIMAXJ(J)
          IF (oFOCEAN_loc(I,J).gt.0.) THEN
            oTOT_CHLO_loc(I,J) = tot_chlo(I,J)
          ELSE
            oTOT_CHLO_loc(I,J)=0.
          END IF
        END DO
      END DO
      CALL INT_OG2AG(oTOT_CHLO_loc,CHL, oWEIGHT, .FALSE.)
#endif

#ifdef TRACERS_GASEXCH_ocean_CO2
!partial CO2 pressure in seawater. Units are uatm.
!defined only over open ocean cells, because this is what is
!involved in gas exchage with the atmosphere.
      oWEIGHT(:,:) = oFOCEAN_loc(:,:)*(1.d0-oRSI(:,:))
      DO J=oJ_0,oJ_1
        DO I=oI_0,oIMAXJ(J)
          IF (oFOCEAN_loc(I,J).gt.0.) THEN
            opCO2_loc(I,J) = pCO2(I,J)
          ELSE
            opCO2_loc(I,J)=0.
          END IF
        END DO
      END DO
      DO NT = 1,NTM
         CALL INT_OG2AG(opCO2_loc,aTRAC(:,:,NT), oWEIGHT, .FALSE.)

         !opco2 is in uatm, convert to kg,CO2/kg,air
         aTRAC(:,:,NT) = aTRAC(:,:,NT) * vol2mass(nt)* 1.d-6 ! ppmv (uatm) -> kg,CO2/kg,air

         !gtracer is first set in TRACER_DRV, then atrac is interpolated
         !here from pco2 in the ocean and later OCNDYN sets gtracer=atrac
         !Therefore in timesetep 0 (cold start) pco2 has not yet been defined
         !and atrac has to be hard coded in here, so that we do not have
         !urealistic tracer flux at the air-sea interface during step0.

         if (nstep.eq.0) aTRAC(:,:,NT) = gtracer(NT,1,:,:)

      END DO
#endif
#endif

      DEALLOCATE(oG0, oS0, oUO1,oVO1, oTOT_CHLO_loc, opCO2_loc)
#ifdef TRACERS_OCEAN
      DEALLOCATE(oTRAC)
#endif

      deallocate(oweight,ofocean_loc)

      RETURN
      END SUBROUTINE OG2AG_TOC2SST

#endif /* OG2AG_TOC2SST_BUNDLE */

#ifdef AG2OG_OCEANS_BUNDLE

      SUBROUTINE AG2OG_oceans
!@sum  AG2OG_oceans: all atmospheric arrays used in the subr. OCEANS are gathered
!       on the atmospheric grid, interpolated to the ocean grid, and scattered
!!      on the ocean grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>IM, aJM=>JM

      USE OCEAN, only : oIM=>IM, oJM=>JM, oFOCEAN=>FOCEAN
     &                , oDXYPO=>DXYPO, OXYP
     &                , IVSPO=>IVSP,IVNPO=>IVNP
     &                , oSINI=>SINIC, oCOSI=>COSIC 
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
#if defined (TRACERS_OceanBiology) && !defined (TRACERS_GASEXCH_ocean)
      USE OCN_TRACER_COM, only: NTM
#else
      USE TRACER_COM, only: NTM
#endif
#endif

      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      use domain_decomp_1d, only : pack_data
      USE OCEANR_DIM, only : ogrid

      USE SEAICE_COM, only : aRSI=>RSI
      USE FLUXES, only : aSOLAR=>SOLAR, aE0=>E0, aEVAPOR=>EVAPOR
     *     , aRUNOSI=>RUNOSI,aERUNOSI=>ERUNOSI,aSRUNOSI=>SRUNOSI
     *     , aFLOWO=>FLOWO,aEFLOWO=>EFLOWO, aAPRESS=>APRESS
     *     , aMELTI=>MELTI,aEMELTI=>EMELTI,aSMELTI=>SMELTI
     *     , aDMUA=>DMUA, aDMVA=>DMVA
     *     , aGMELT=>GMELT, aEGMELT=>EGMELT
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
     *     , aTRFLOWO=>TRFLOWO, aTREVAPOR=>TREVAPOR
     *     , aTRUNOSI=>TRUNOSI, aTRMELTI=>TRMELTI
     *     , aTRGMELT=>TRGMELT
#ifdef TRACERS_DRYDEP
     *     , aTRDRYDEP=>TRDRYDEP
#endif
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      USE FLUXES, only : aTRGASEX=>TRGASEX
      USE TRACER_GASEXCH_COM, only: tracflx_glob
#endif
#ifdef OBIO_RAD_coupling
      USE RAD_COM, only : avisdir    => FSRDIR
     *                   ,asrvissurf => SRVISSURF
     *                   ,avisdif    => FSRDIF
     *                   ,anirdir    => DIRNIR
     *                   ,anirdif    => DIFNIR
#endif
#ifdef TRACERS_OceanBiology
      USE RAD_COM, only : aCOSZ1=>COSZ1
      USE PBLCOM, only : awind=>wsavg
#endif

      USE OFLUXES, only : oRSI, oSOLAR, oE0, oEVAPOR
     *     , oRUNOSI, oERUNOSI, oSRUNOSI
     *     , oFLOWO, oEFLOWO, oAPRESS
     *     , oMELTI, oEMELTI, oSMELTI
     *     , oDMUA,oDMVA
     *     , oGMELT, oEGMELT
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
     *     , oTRFLOWO, oTREVAPOR
     *     , oTRUNOSI, oTRMELTI
     *     , oTRGMELT
#ifdef TRACERS_DRYDEP
     *     , oTRDRYDEP
#endif
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      USE MODEL_COM, only: nstep=>itime
      USE OFLUXES, only : oTRGASEX
#endif
#ifdef OBIO_RAD_coupling
      USE obio_forc, only : ovisdir,ovisdif,onirdir,onirdif
#endif
#ifdef TRACERS_OceanBiology
      USE obio_forc, only : osolz
      USE obio_forc, only : owind
#endif

      Use GEOM,  only : AXYP,aIMAXJ=>IMAXJ

      USE MODEL_COM, only : aFOCEAN_loc=>FOCEAN

      USE regrid_com, only : remap_A2O
      USE bundle_arrays

      IMPLICIT NONE
      REAL*8, allocatable :: aWEIGHT(:,:),oWEIGHT(:,:)
      REAL*8, allocatable :: aWEIGHT1(:,:),oWEIGHT1(:,:)
      REAL*8, allocatable :: aWEIGHT2(:,:),oWEIGHT2(:,:)
      REAL*8, allocatable :: aE0tmp(:,:),oE0tmp(:,:),
     &     aEVAPORtmp(:,:),oEVAPORtmp(:,:),
     &     aSOLAR1tmp(:,:),oSOLAR1tmp(:,:),
     &     aSOLAR3tmp(:,:),oSOLAR3tmp(:,:),
     &     aDMUA1tmp(:,:),oDMUA1tmp(:,:),
     &     aDMVA1tmp(:,:),oDMVA1tmp(:,:)
      REAL*8, allocatable :: o8_glob(:,:)
      REAL*4, allocatable :: o4_glob(:,:)
      INTEGER, PARAMETER :: NSTYPE=4
      INTEGER I,J,N
      INTEGER aJ_0,aJ_1, aI_0,aI_1
      INTEGER oJ_0,oJ_1
      REAL*8,
     *     DIMENSION(aGRID%I_STRT:aGRID%I_STOP,
     *     aGRID%J_STRT:aGRID%J_STOP)::
     *     aFact
      REAL*8 :: oDMUA1sp,oDMVA1sp,oDMUA1np,oDMVA1np
      character*80 :: name,title
      type (lookup_str) :: lstr


      allocate(aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(aweight1(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight1(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(aweight2(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight2(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

      !-- initializing "lstr"
      call ba_init( lstr, aGRID%I_STRT_HALO, aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO, aGRID%J_STOP_HALO,
     &     oGRID%I_STRT_HALO, oGRID%I_STOP_HALO,
     &     oGRID%J_STRT_HALO, oGRID%J_STOP_HALO )

      aJ_0 = aGRID%j_STRT
      aJ_1 = aGRID%j_STOP
      aI_0 = aGRID%I_STRT
      aI_1 = aGRID%I_STOP
      oJ_0 = oGRID%j_STRT
      oJ_1 = oGRID%j_STOP


      call ba_add( lstr, aRSI, oRSI)
      
      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aFact(I,J) = 1.d0/(AXYP(I,J)*aFOCEAN_loc(I,J))
            aFLOWO(I,J) = aFLOWO(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      call ba_add( lstr, aFLOWO,oFLOWO)

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aEFLOWO(I,J) = aEFLOWO(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      call ba_add( lstr, aEFLOWO,oEFLOWO)

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aMELTI(I,J) = aMELTI(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      call ba_add( lstr, aMELTI,oMELTI)

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aEMELTI(I,J) = aEMELTI(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      call ba_add( lstr, aEMELTI,oEMELTI)

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aSMELTI(I,J) = aSMELTI(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      call ba_add( lstr, aSMELTI,oSMELTI)

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aGMELT(I,J) = aGMELT(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      call ba_add(lstr, aGMELT,oGMELT)

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aFact(I,J) = 1.d0/AXYP(I,J)
            aEGMELT(I,J) = aEGMELT(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      call ba_add( lstr, aEGMELT,oEGMELT)

      call ba_add( lstr, aAPRESS,oAPRESS)

      aWEIGHT(:,:) = aRSI(:,:)
      call ba_add( lstr, aWEIGHT, oWEIGHT)  
      call ba_add( lstr, aRUNOSI, oRUNOSI, aWEIGHT, oWEIGHT)
      call ba_add( lstr, aERUNOSI, oERUNOSI, aWEIGHT, oWEIGHT)
      call ba_add( lstr, aSRUNOSI, oSRUNOSI, aWEIGHT, oWEIGHT)

      aWEIGHT1(:,:) = 1.d0 - aRSI(:,:)
      call ba_add( lstr, aWEIGHT1, oWEIGHT1)

!     allocate temporary arrays 
      allocate(aE0tmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oE0tmp(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(aEVAPORtmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oEVAPORtmp(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(aSOLAR1tmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oSOLAR1tmp(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(aSOLAR3tmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oSOLAR3tmp(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(aDMUA1tmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oDMUA1tmp(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(aDMVA1tmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oDMVA1tmp(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

!     copy E0(:,:,1) in temporary 2d array on atm grid, 
!     and return temporary 2d array on ocean grid
      aE0tmp(:,:)=aE0(:,:,1)
      call ba_add(lstr, aE0tmp, oE0tmp, aWEIGHT1, oWEIGHT1)

!     copy EVAPOR(:,:,1) in temporary 2d array on atm grid, 
!     and return temporary 2d array on ocean grid
      aEVAPORtmp(:,:)=aEVAPOR(:,:,1)
      call ba_add(lstr, aEVAPORtmp, oEVAPORtmp, aWEIGHT1, oWEIGHT1)

!     copy SOLAR(1,:,:) in temporary 2d array on atm grid, 
!     and return temporary 2d array on ocean grid
      aSOLAR1tmp(:,:)=aSOLAR(1,:,:)
      call ba_add( lstr, aSOLAR1tmp, oSOLAR1tmp, aWEIGHT1, oWEIGHT1)

      aWEIGHT2(:,:) = aRSI(:,:)
      call ba_add( lstr, aWEIGHT2, oWEIGHT2)

!     copy SOLAR(1,:,:) in temporary 2d array on atm grid, 
!     and return temporary 2d array on ocean grid
      aSOLAR3tmp(:,:)=aSOLAR(3,:,:)
      call ba_add( lstr, aSOLAR3tmp, oSOLAR3tmp, aWEIGHT2, oWEIGHT2)

#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
      aWEIGHT(:,:) = 1.d0
      DO N=1,NTM
        DO J=aJ_0,aJ_1
          DO I=aI_0,aIMAXJ(J)
            IF (aFOCEAN_loc(I,J).gt.0.) THEN
              aFact(I,J) = 1.d0/(AXYP(I,J)*aFOCEAN_loc(I,J))
              aTRFLOWO(N,I,J) = aTRFLOWO(N,I,J)*aFact(I,J)
            END IF
          END DO
        END DO
      END DO
      CALL INT_AG2OG(aTRFLOWO,oTRFLOWO, aWEIGHT, NTM)

      DO N=1,NTM
        DO J=aJ_0,aJ_1
          DO I=aI_0,aIMAXJ(J)
            IF (aFOCEAN_loc(I,J).gt.0.) THEN
              aTRMELTI(N,I,J) = aTRMELTI(N,I,J)*aFact(I,J)
            END IF
          END DO
        END DO
      END DO
      CALL INT_AG2OG(aTRMELTI,oTRMELTI, aWEIGHT, NTM)

      aWEIGHT(:,:) = aRSI(:,:)
      CALL INT_AG2OG(aTRUNOSI,oTRUNOSI, aWEIGHT, NTM)

      aWEIGHT(:,:) = 1.d0
      DO N=1,NTM
        DO J=aJ_0,aJ_1
          DO I=aI_0,aIMAXJ(J)
            IF (aFOCEAN_loc(I,J).gt.0.) THEN
              aFact(I,J) = 1.d0/AXYP(I,J)
              aTRGMELT(N,I,J) = aTRGMELT(N,I,J)*aFact(I,J)
            END IF
          END DO
        END DO
      END DO
      CALL INT_AG2OG(aTRGMELT,oTRGMELT, aWEIGHT, NTM)
      DO N=1,NTM
        oTRGMELT(N,:,:) = oTRGMELT(N,:,:)*OXYP(:,:)
      END DO

      aWEIGHT(:,:) = 1.d0 - aRSI(:,:)
      CALL INT_AG2OG(aTREVAPOR,oTREVAPOR, aWEIGHT, NTM, NSTYPE,1)

#ifdef TRACERS_DRYDEP
      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(aTRDRYDEP,oTRDRYDEP, aWEIGHT, NTM, NSTYPE,1)
#endif
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(aTRGASEX,oTRGASEX, tracflx_glob, aWEIGHT
     *             , NTM, NSTYPE,1)
#endif

#ifdef OBIO_RAD_coupling
      aWEIGHT(:,:) = aFOCEAN_loc(:,:)
      CALL INT_AG2OG(aVISDIR,aSRVISSURF,oVISDIR, aWEIGHT)

      CALL INT_AG2OG(aVISDIF,oVISDIF, aWEIGHT)

      CALL INT_AG2OG(aNIRDIR,oNIRDIR, aWEIGHT)

      CALL INT_AG2OG(aNIRDIF,oNIRDIF, aWEIGHT)
#endif

#ifdef TRACERS_OceanBiology
      aWEIGHT(:,:) = aFOCEAN_loc(:,:)
      CALL INT_AG2OG(aCOSZ1,oSOLZ, aWEIGHT)

      CALL INT_AG2OG(aWIND,oWIND, aWEIGHT)
#endif

!      aWEIGHT(:,:) = 1.d0
      
      do j=aGRID%J_STRT_HALO,aGRID%J_STOP_HALO
         do i=aGRID%I_STRT_HALO,aGRID%I_STOP_HALO
            aDMUA1tmp(i,j) = 0.
            if (aFOCEAN_loc(i,j).gt.0.) then
               aDMUA1tmp(i,j) = aDMUA(i,j,1) 
            endif
         enddo
      enddo

      call ba_add( lstr, aDMUA1tmp, oDMUA1tmp)

      do j=aGRID%J_STRT_HALO,aGRID%J_STOP_HALO
         do i=aGRID%I_STRT_HALO,aGRID%I_STOP_HALO
            aDMVA1tmp(i,j) = 0.
            if (aFOCEAN_loc(i,j).gt.0.) then
               aDMVA1tmp(i,j) = aDMVA(i,j,1) 
            endif
         enddo
      enddo

      call ba_add( lstr, aDMVA1tmp, oDMVA1tmp)

c*   actual interpolation here
      call bundle_interpolation(remap_A2O,lstr)
c*

      oEGMELT(:,:) = oEGMELT(:,:)*OXYP(:,:)

      oE0(:,:,1)=oE0tmp(:,:)  

      oEVAPOR(:,:,1)=oEVAPORtmp(:,:)  

      oSOLAR(1,:,:)=oSOLAR1tmp(:,:)  

      oSOLAR(3,:,:)=oSOLAR3tmp(:,:) 

      if (ogrid%HAVE_SOUTH_POLE) then
         oDMUA1sp =  SUM(oDMUA1tmp(:,  1)*oCOSI(:))*2/oIM
         oDMVA1sp = -SUM(oDMUA1tmp(:,  1)*oSINI(:))*2/oIM
         oDMUA1tmp(1,  1) = oDMUA1sp
         oDMVA1tmp(1,  1) = oDMVA1sp
      endif
      if (ogrid%HAVE_NORTH_POLE) then
         oDMUA1np =  SUM(oDMUA1tmp(:,oJM)*oCOSI(:))*2/oIM
         oDMVA1np =  SUM(oDMUA1tmp(:,oJM)*oSINI(:))*2/oIM
         oDMUA1tmp(1,oJM) = oDMUA1np
         oDMVA1tmp(1,oJM) = oDMVA1np
      endif

      oDMUA(:,:,1)=oDMUA1tmp(:,:)
      oDMVA(:,:,1)=oDMVA1tmp(:,:)
 
      deallocate(aweight,oweight,
     &     aweight1,oweight1,
     &     aweight2,oweight2)
      deallocate(aE0tmp,oE0tmp,
     &     aEVAPORtmp,oEVAPORtmp,
     &     aSOLAR1tmp,oSOLAR1tmp,
     &     aSOLAR3tmp,oSOLAR3tmp,
     &     aDMUA1tmp,oDMUA1tmp,
     &     aDMVA1tmp,oDMVA1tmp)


      END SUBROUTINE AG2OG_oceans
#else

      SUBROUTINE AG2OG_oceans
!@sum  AG2OG_oceans: all atmospheric arrays used in the subr. OCEANS are gathered
!       on the atmospheric grid, interpolated to the ocean grid, and scattered
!!      on the ocean grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>IM, aJM=>JM

      USE OCEAN, only : oIM=>IM, oJM=>JM, oFOCEAN=>FOCEAN
     *                , oDXYPO=>DXYPO, OXYP
     *                , IVSPO=>IVSP,IVNPO=>IVNP

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
#if defined (TRACERS_OceanBiology) && !defined (TRACERS_GASEXCH_ocean)
      USE OCN_TRACER_COM, only: NTM
#else
      USE TRACER_COM, only: NTM
#endif
#endif

      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      use domain_decomp_1d, only : pack_data
      USE OCEANR_DIM, only : ogrid

      USE SEAICE_COM, only : aRSI=>RSI
      USE FLUXES, only : aSOLAR=>SOLAR, aE0=>E0, aEVAPOR=>EVAPOR
     *     , aRUNOSI=>RUNOSI,aERUNOSI=>ERUNOSI,aSRUNOSI=>SRUNOSI
     *     , aFLOWO=>FLOWO,aEFLOWO=>EFLOWO, aAPRESS=>APRESS
     *     , aMELTI=>MELTI,aEMELTI=>EMELTI,aSMELTI=>SMELTI
     *     , aDMUA=>DMUA, aDMVA=>DMVA
     *     , aGMELT=>GMELT, aEGMELT=>EGMELT
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
     *     , aTRFLOWO=>TRFLOWO, aTREVAPOR=>TREVAPOR
     *     , aTRUNOSI=>TRUNOSI, aTRMELTI=>TRMELTI
     *     , aTRGMELT=>TRGMELT
#ifdef TRACERS_DRYDEP
     *     , aTRDRYDEP=>TRDRYDEP
#endif
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      USE FLUXES, only : aTRGASEX=>TRGASEX
      USE TRACER_GASEXCH_COM, only: tracflx_glob
#endif
#ifdef OBIO_RAD_coupling
      USE RAD_COM, only : avisdir    => FSRDIR
     *                   ,asrvissurf => SRVISSURF
     *                   ,avisdif    => FSRDIF
     *                   ,anirdir    => DIRNIR
     *                   ,anirdif    => DIFNIR
#endif
#ifdef TRACERS_OceanBiology
      USE RAD_COM, only : aCOSZ1=>COSZ1
      USE PBLCOM, only : awind=>wsavg
#endif

      USE OFLUXES, only : oRSI, oSOLAR, oE0, oEVAPOR
     *     , oRUNOSI, oERUNOSI, oSRUNOSI
     *     , oFLOWO, oEFLOWO, oAPRESS
     *     , oMELTI, oEMELTI, oSMELTI
     *     , oDMUA,oDMVA
     *     , oGMELT, oEGMELT
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
     *     , oTRFLOWO, oTREVAPOR
     *     , oTRUNOSI, oTRMELTI
     *     , oTRGMELT
#ifdef TRACERS_DRYDEP
     *     , oTRDRYDEP
#endif
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      USE MODEL_COM, only: nstep=>itime
      USE OFLUXES, only : oTRGASEX
#endif
#ifdef OBIO_RAD_coupling
      USE obio_forc, only : ovisdir,ovisdif,onirdir,onirdif
#endif
#ifdef TRACERS_OceanBiology
      USE obio_forc, only : osolz
      USE obio_forc, only : owind
#endif

      Use GEOM,  only : AXYP,aIMAXJ=>IMAXJ

      USE MODEL_COM, only : aFOCEAN_loc=>FOCEAN

      USE INT_AG2OG_MOD, only : INT_AG2OG

      IMPLICIT NONE

      INTEGER, PARAMETER :: NSTYPE=4
      INTEGER I,J, N
      INTEGER aJ_0,aJ_1, aI_0,aI_1
      INTEGER oJ_0,oJ_1

      REAL*8,
     * DIMENSION(aGRID%I_STRT:aGRID%I_STOP,aGRID%J_STRT:aGRID%J_STOP)::
     * aFact

      REAL*8, allocatable :: aweight(:,:)
      real*8, dimension(oIM,oJM) :: otest_glob
      real*4, dimension(oIM,oJM) :: tr4

      character*80 :: title,name

      allocate (aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO) )

      aJ_0 = aGRID%j_STRT
      aJ_1 = aGRID%j_STOP
      aI_0 = aGRID%I_STRT
      aI_1 = aGRID%I_STOP

      oJ_0 = oGRID%j_STRT
      oJ_1 = oGRID%j_STOP

      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(aRSI,oRSI, aWEIGHT)
      
      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aFact(I,J) = 1.d0/(AXYP(I,J)*aFOCEAN_loc(I,J))
            aFLOWO(I,J) = aFLOWO(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(aFLOWO,oFLOWO, aWEIGHT)

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aEFLOWO(I,J) = aEFLOWO(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(aEFLOWO,oEFLOWO, aWEIGHT)

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aMELTI(I,J) = aMELTI(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(aMELTI,oMELTI, aWEIGHT)

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aEMELTI(I,J) = aEMELTI(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(aEMELTI,oEMELTI, aWEIGHT)

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aSMELTI(I,J) = aSMELTI(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(aSMELTI,oSMELTI, aWEIGHT)

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aGMELT(I,J) = aGMELT(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(aGMELT,oGMELT, aWEIGHT)

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (aFOCEAN_loc(I,J).gt.0.) THEN
            aFact(I,J) = 1.d0/AXYP(I,J)
            aEGMELT(I,J) = aEGMELT(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(aEGMELT,oEGMELT, aWEIGHT)
      oEGMELT(:,:) = oEGMELT(:,:)*OXYP(:,:)

      aWEIGHT(:,:) = aRSI(:,:)
      CALL INT_AG2OG(aRUNOSI,oRUNOSI, aWEIGHT)

      CALL INT_AG2OG(aERUNOSI,oERUNOSI, aWEIGHT)

      CALL INT_AG2OG(aSRUNOSI,oSRUNOSI, aWEIGHT)

      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(aAPRESS,oAPRESS, aWEIGHT)

      aWEIGHT(:,:) = 1.d0 - aRSI(:,:)
      CALL INT_AG2OG(aE0,oE0, aWEIGHT, NSTYPE,1)

      CALL INT_AG2OG(aEVAPOR,oEVAPOR, aWEIGHT, NSTYPE,1)

      CALL INT_AG2OG(aSOLAR,oSOLAR, aWEIGHT, 3,3,1)

      aWEIGHT(:,:) = aRSI(:,:)
      CALL INT_AG2OG(aSOLAR,oSOLAR, aWEIGHT, 3,3,3)

#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
      aWEIGHT(:,:) = 1.d0
      DO N=1,NTM
        DO J=aJ_0,aJ_1
          DO I=aI_0,aIMAXJ(J)
            IF (aFOCEAN_loc(I,J).gt.0.) THEN
              aFact(I,J) = 1.d0/(AXYP(I,J)*aFOCEAN_loc(I,J))
              aTRFLOWO(N,I,J) = aTRFLOWO(N,I,J)*aFact(I,J)
            END IF
          END DO
        END DO
      END DO
      CALL INT_AG2OG(aTRFLOWO,oTRFLOWO, aWEIGHT, NTM)

      DO N=1,NTM
        DO J=aJ_0,aJ_1
          DO I=aI_0,aIMAXJ(J)
            IF (aFOCEAN_loc(I,J).gt.0.) THEN
              aTRMELTI(N,I,J) = aTRMELTI(N,I,J)*aFact(I,J)
            END IF
          END DO
        END DO
      END DO
      CALL INT_AG2OG(aTRMELTI,oTRMELTI, aWEIGHT, NTM)

      aWEIGHT(:,:) = aRSI(:,:)
      CALL INT_AG2OG(aTRUNOSI,oTRUNOSI, aWEIGHT, NTM)

      aWEIGHT(:,:) = 1.d0
      DO N=1,NTM
        DO J=aJ_0,aJ_1
          DO I=aI_0,aIMAXJ(J)
            IF (aFOCEAN_loc(I,J).gt.0.) THEN
              aFact(I,J) = 1.d0/AXYP(I,J)
              aTRGMELT(N,I,J) = aTRGMELT(N,I,J)*aFact(I,J)
            END IF
          END DO
        END DO
      END DO
      CALL INT_AG2OG(aTRGMELT,oTRGMELT, aWEIGHT, NTM)
      DO N=1,NTM
        oTRGMELT(N,:,:) = oTRGMELT(N,:,:)*OXYP(:,:)
      END DO

      aWEIGHT(:,:) = 1.d0 - aRSI(:,:)
      CALL INT_AG2OG(aTREVAPOR,oTREVAPOR, aWEIGHT, NTM, NSTYPE,1)

#ifdef TRACERS_DRYDEP
      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(aTRDRYDEP,oTRDRYDEP, aWEIGHT, NTM, NSTYPE,1)
#endif
#endif
#endif
#ifdef TRACERS_GASEXCH_ocean
      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(aTRGASEX,oTRGASEX, tracflx_glob, aWEIGHT
     *             , NTM, NSTYPE,1)
#endif

#ifdef OBIO_RAD_coupling
      aWEIGHT(:,:) = aFOCEAN_loc(:,:)
      CALL INT_AG2OG(aVISDIR,aSRVISSURF,oVISDIR, aWEIGHT)

      CALL INT_AG2OG(aVISDIF,oVISDIF, aWEIGHT)

      CALL INT_AG2OG(aNIRDIR,oNIRDIR, aWEIGHT)

      CALL INT_AG2OG(aNIRDIF,oNIRDIF, aWEIGHT)
#endif

#ifdef TRACERS_OceanBiology
      aWEIGHT(:,:) = aFOCEAN_loc(:,:)
      CALL INT_AG2OG(aCOSZ1,oSOLZ, aWEIGHT)

      CALL INT_AG2OG(aWIND,oWIND, aWEIGHT)
#endif
      aWEIGHT(:,:) = 1.d0
      
      CALL INT_AG2OG(aDMUA,aDMVA,oDMUA,oDMVA, aWEIGHT,aFOCEAN_loc
     *              ,NSTYPE,1)

      deallocate(aweight)

      END SUBROUTINE AG2OG_oceans
#endif /*ag2og_oceans_bundle */


c      subroutine debug_acoupling(name,ar8)
c!@auth D. Gueyffier
c      USE RESOLUTION, only : aIM=>IM, aJM=>JM
c      use domain_decomp_atm, only : agrid=>grid,pack_data,am_i_root
c      real*8, intent(in) :: ar8(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
c     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO)
c      real*8, allocatable :: ar8glob(:,:,:)
c      real*4, allocatable :: ar4(:,:,:)
c      character*80 :: name,title,fname
c
c      allocate(ar8glob(aIM,aJM,6),ar4(aIM,aJM,6))
c
c      call pack_data(agrid,ar8,ar8glob)
c
c      if (am_i_root()) then
c         ar4=ar8glob
c         title="testa-"//trim(name)
c         fname=title
c         open(2000,FILE=fname,FORM='unformatted', STATUS='unknown')
c         write(2000) title,ar4
c         close(2000)
c         name=""
c      endif
c
c      deallocate(ar8glob,ar4)
c
c      end subroutine debug_acoupling


c      subroutine debug_ocoupling(name,or8)
c!@auth D. Gueyffier
c      USE OCEAN, only : oIM=>IM, oJM=>JM, oFOCEAN=>FOCEAN
c      use domain_decomp_1d, only : pack_data,am_i_root
c      USE OCEANR_DIM, only : ogrid
c      real*8, intent(in) :: or8(oIM,
c     &     oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
c      real*8, allocatable :: or8glob(:,:)
c      real*4, allocatable :: or4(:,:)
c      character*80 :: title,name,fname
c
c      allocate(or8glob(oIM,oJM),or4(oIM,oJM))
c
cc      write(*,*) "max o dist=",maxval(
cc     &     or8(:,oGRID%J_STRT:oGRID%J_STOP)
cc     &                                ),
cc     &           "min odist=",minval(
cc     &     or8(:,oGRID%J_STRT:oGRID%J_STOP)
cc     &                               )
c
c      call pack_data(ogrid,or8,or8glob)
c
c      if (am_i_root()) then
c         or4=or8glob
c         title="testo-"//trim(name)
c         fname=title
c         open(30,FILE=fname,FORM='unformatted', STATUS='unknown')
c         write(30) title,or4
c         write(*,*) "max "//trim(name),"=",maxval(or8glob),
c     &        "min "//trim(name),"=",minval(or8glob),
c     &        "maxloc "//trim(name),"=",maxloc(or8glob),
c     &        "minloc "//trim(name),"=",minloc(or8glob)
c 
c         close(30)
c         name=""
c      endif
c
c      deallocate(or8glob,or4)
c
c      end subroutine debug_ocoupling
c

#ifdef OG2AG_OCEANS_BUNDLE      
      SUBROUTINE OG2AG_oceans
!@sum  OG2AG_oceans: ocean arrays for sea ice formation calculated in the
!       subr. OCEANS are gathered on the ocean grid, interpolated to the
!!      atmospheric grid, and scattered on the atmospheric ocean grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>IM, aJM=>JM

      USE OCEAN, only : oIM=>IM,oJM=>JM, oFOCEAN=>FOCEAN

#ifdef TRACERS_ON
      USE TRACER_COM, only: ntm_atm=>ntm
#endif
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only: ntm
#endif

      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only : UNPACK_DATA
      USE OCEANR_DIM, only : ogrid

      USE MODEL_COM, ONLY : aFOCEAN_loc=>FOCEAN

      USE FLUXES, only : aDMSI=>DMSI,aDHSI=>DHSI,aDSSI=>DSSI
#ifdef TRACERS_ON
     *     , aDTRSI=>DTRSI
#endif

      USE OFLUXES, only : oDMSI, oDHSI, oDSSI
#ifdef TRACERS_OCEAN
     *     , oDTRSI
#endif

      USE regrid_com, only : remap_O2A
      USE bundle_arrays
      IMPLICIT NONE

      REAL*8,
     * DIMENSION(oIM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)::
     *     oFOCEAN_loc

      REAL*8, allocatable :: aWEIGHT(:,:),oWEIGHT(:,:)
      REAL*8, allocatable :: oDMSItmp(:,:,:),aDMSItmp(:,:,:)
      REAL*8, allocatable :: oDHSItmp(:,:,:),aDHSItmp(:,:,:)
      REAL*8, allocatable :: oDSSItmp(:,:,:),aDSSItmp(:,:,:)
      integer :: i,j,l
      type (lookup_str) :: lstr


      allocate(aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

      allocate(oDMSItmp(2,oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(oDHSItmp(2,oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(oDSSItmp(2,oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

      allocate(aDMSItmp(2,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(aDHSItmp(2,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(aDSSItmp(2,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))


      write(*,*) "OG2AG_OCEANS_BUNDLE"

      !-- initializing "lstr"
      call ba_init( lstr, 1, oIM,
     &     oGRID%J_STRT_HALO, oGRID%J_STOP_HALO,
     &     aGRID%I_STRT_HALO, aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO, aGRID%J_STOP_HALO)


      CALL UNPACK_DATA (oGRID,oFOCEAN,oFOCEAN_loc)

      oWEIGHT(:,:) = oFOCEAN_loc(:,:)
      call ba_add(lstr, oWEIGHT, aWEIGHT)

      oDMSItmp=oDMSI
      do L=1,2
         if (ogrid%HAVE_NORTH_POLE) 
     &        oDMSItmp(L,2:oIM,oJM) = oDMSItmp(L,1,oJM)
      enddo
      call ba_add_new(lstr, oDMSItmp, aDMSItmp, shape(oDMSItmp),
     &     'lij', oWEIGHT, aWEIGHT) 
c      CALL INT_OG2AG(oDMSI,aDMSI, oWEIGHT, 2)

      oDHSItmp=oDHSI
      do L=1,2
         if (ogrid%HAVE_NORTH_POLE) 
     &        oDHSItmp(L,2:oIM,oJM) = oDHSItmp(L,1,oJM)
      enddo
      call ba_add_new(lstr, oDHSItmp, aDHSItmp, shape(oDHSItmp),
     &     'lij', oWEIGHT, aWEIGHT) 
c      CALL INT_OG2AG(oDHSI,aDHSI, oWEIGHT, 2)

      oDSSItmp=oDSSI
      do L=1,2
         if (ogrid%HAVE_NORTH_POLE) 
     &        oDSSItmp(L,2:oIM,oJM) = oDSSItmp(L,1,oJM)
      enddo
      call ba_add_new(lstr, oDSSItmp, aDSSItmp, shape(oDSSItmp),
     &     'lij', oWEIGHT, aWEIGHT) 
c      CALL INT_OG2AG(oDSSI,aDSSI, oWEIGHT, 2)

#if (defined TRACERS_OCEAN) && (defined TRACERS_ON)

      IF (NTM == NTM_ATM) THEN

        oWEIGHT(:,:) = oFOCEAN_loc(:,:)
        CALL INT_OG2AG(oDTRSI,aDTRSI, oWEIGHT, 2, NTM)

      ELSE
C**** if number of ocean and atm. tracers are not the same
C**** do something in here

      END IF
#endif

c*   actual interpolation here
      call bundle_interpolation(remap_O2A,lstr)
c*


      do j=agrid%j_strt,agrid%j_stop
      do i=agrid%i_strt,agrid%i_stop
        if(aFOCEAN_loc(i,j) > 0.) then
           aDMSI(:,i,j) = aDMSItmp(:,i,j)
           aDHSI(:,i,j) = aDHSItmp(:,i,j)
           aDSSI(:,i,j) = aDSSItmp(:,i,j)
        endif
      enddo
      enddo

      deallocate(oWEIGHT, aWEIGHT)
      deallocate(oDMSItmp,aDMSItmp)
      deallocate(oDHSItmp,aDHSItmp)
      deallocate(oDSSItmp,aDSSItmp)

      RETURN
      END SUBROUTINE OG2AG_oceans

#else
      SUBROUTINE OG2AG_oceans
!@sum  OG2AG_oceans: ocean arrays for sea ice formation calculated in the
!       subr. OCEANS are gathered on the ocean grid, interpolated to the
!!      atmospheric grid, and scattered on the atmospheric ocean grid
!@auth Larissa Nazarenko
!@ver  1.0

      USE RESOLUTION, only : aIM=>IM, aJM=>JM

      USE OCEAN, only : oIM=>IM,oJM=>JM, oFOCEAN=>FOCEAN

#ifdef TRACERS_ON
      USE TRACER_COM, only: ntm_atm=>ntm
#endif
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only: ntm
#endif

      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only : UNPACK_DATA
      USE OCEANR_DIM, only : ogrid

      USE MODEL_COM, ONLY : aFOCEAN_loc=>FOCEAN

      USE FLUXES, only : aDMSI=>DMSI,aDHSI=>DHSI,aDSSI=>DSSI
#ifdef TRACERS_ON
     *     , aDTRSI=>DTRSI
#endif

      USE OFLUXES, only : oDMSI, oDHSI, oDSSI
#ifdef TRACERS_OCEAN
     *     , oDTRSI
#endif

      USE INT_OG2AG_MOD, only : INT_OG2AG

      IMPLICIT NONE

      REAL*8,
     * DIMENSION(oIM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)::
     * oWEIGHT, oFOCEAN_loc

      CALL UNPACK_DATA (oGRID,oFOCEAN,oFOCEAN_loc)

      oWEIGHT(:,:) = oFOCEAN_loc(:,:)
      CALL INT_OG2AG(oDMSI,aDMSI, oWEIGHT, 2)

      CALL INT_OG2AG(oDHSI,aDHSI, oWEIGHT, 2)

      CALL INT_OG2AG(oDSSI,aDSSI, oWEIGHT, 2)

#if (defined TRACERS_OCEAN) && (defined TRACERS_ON)

      IF (NTM == NTM_ATM) THEN

        oWEIGHT(:,:) = oFOCEAN_loc(:,:)
        CALL INT_OG2AG(oDTRSI,aDTRSI, oWEIGHT, 2, NTM)

      ELSE
C**** if number of ocean and atm. tracers are not the same
C**** do something in here

      END IF
#endif
      RETURN
      END SUBROUTINE OG2AG_oceans

#endif /* OG2AG_OCEANS_BUNDLE */

      MODULE IGOG_regrid_info
!@sum IGOG_info saves instances of hntrp_type for use in
!@+   ice <-> ocean regrids
      use hntrp_mod, only : hntrp_type
      implicit none
      save

      type(hntrp_type) :: hntrp_i2o_u ! ice u C -> ocn u C
      type(hntrp_type) :: hntrp_i2o_v ! ice v C -> ocn v C
      logical :: hntrp_i2o_uv_need_init = .true.

      type(hntrp_type) :: hntrp_o2i_u ! ocn u C -> ice u A
      type(hntrp_type) :: hntrp_o2i_v ! ocn v C -> ice v A
      logical :: hntrp_o2i_uv_need_init = .true.

      END MODULE IGOG_regrid_info


      SUBROUTINE IG2OG_oceans
!@sum IG2OG_oceans interpolates relevant DYNSI outputs to the
!@+   ocean grid
!@auth M. Kelley
      use hntrp_mod
      use IGOG_regrid_info
      use icedyn_com, only : iDMUI=>DMUI, iDMVI=>DMVI
      use ofluxes,    only : oDMUI,oDMVI
      USE ICEDYN, only : iIM=>imicdyn,iJM=>jmicdyn
      USE OCEAN,  only : oIM=>im,oJM=>jm
      USE ICEDYN, only : iDLATM=>DLATM
      USE OCEAN,  only : oDLATM=>DLATM
      USE ICEDYN,     only : iGRID=>grid_icdyn
      USE OCEANR_DIM, only : oGRID
      USE DOMAIN_DECOMP_1D, only : BAND_PACK
      implicit none
      real*8, dimension(:,:), allocatable ::
     &     ones_band,idmui_band,idmvi_band
      integer :: jmin,jmax

      if(iIM.eq.oIM .and. iJM.eq.oJM) then
        oDMUI(:,:) = iDMUI(:,:)
        oDMVI(:,:) = iDMVI(:,:)
        return
      endif

      if(hntrp_i2o_uv_need_init) then
        call Init_Hntrp_Type(hntrp_i2o_u,
     &     iGRID, .5d0,iDLATM,
     &     oGRID, .5d0,oDLATM,
     &     0.d0)
        call Init_Hntrp_Type(hntrp_i2o_v,
     &     iGRID, 0.d0,iDLATM,
     &     oGRID, 0.d0,oDLATM,
     &     0.d0,
     &     JMA_4interp=iJM-1,JMB_4interp=oJM-1) ! for secondary lats
        hntrp_i2o_uv_need_init = .false.
      endif

c regrid DMUI from ice C to ocn C
      jmin = hntrp_i2o_u%bpack%jband_strt
      jmax = hntrp_i2o_u%bpack%jband_stop
      ALLOCATE(iDMUI_band(iIM,jmin:jmax),ones_band(iIM,jmin:jmax))
      ones_band(:,:) = 1d0
      call BAND_PACK (hntrp_i2o_u%bpack, iDMUI, iDMUI_band)
      if(iGRID%have_north_pole) then
        iDMUI_band(:,iJM) = 0.
      endif
      call HNTR8_band (ones_band, iDMUI_band, hntrp_i2o_u, oDMUI)
      deallocate(iDMUI_band,ones_band)

c regrid DMVI from ice C to ocn C
      jmin = hntrp_i2o_v%bpack%jband_strt
      jmax = hntrp_i2o_v%bpack%jband_stop
      ALLOCATE(iDMVI_band(iIM,jmin:jmax),ones_band(iIM,jmin:jmax))
      ones_band(:,:) = 1d0
      call BAND_PACK (hntrp_i2o_v%bpack, iDMVI, iDMVI_band)
      if(iGRID%have_north_pole) then
        iDMVI_band(:,iJM) = 0.
      endif
      call HNTR8_band (ones_band, iDMVI_band, hntrp_i2o_v, oDMVI)
      deallocate(iDMVI_band,ones_band)

      if(oGRID%have_north_pole) then ! INT_AG2OG_Vector2 set them to zero
        oDMUI(:,oJM) = 0.0
        oDMVI(:,oJM) = 0.0
      endif

      return
      END SUBROUTINE IG2OG_oceans

      SUBROUTINE OG2IG_uvsurf
!@sum OG2IG_uvsurf interpolates ocean surface velocity to the
!@+   DYNSI A-grid
!@auth M. Kelley
      use hntrp_mod
      use IGOG_regrid_info
      USE RESOLUTION, only : aIM=>im,aJM=>jm
      USE ICEDYN, only     : iIM=>imicdyn,iJM=>jmicdyn
      USE OCEAN,  only     : oIM=>im,oJM=>jm
      USE FLUXES, only : atm_uosurf=>uosurf,atm_vosurf=>vosurf
      USE ICEDYN_COM, only : uosurf,vosurf
      USE ICEDYN, only : iDLATM=>DLATM
      USE OCEAN,  only : oDLATM=>DLATM
      USE ICEDYN,     only : iGRID=>grid_icdyn
      USE OCEANR_DIM, only : oGRID
      USE DOMAIN_DECOMP_1D, only : BAND_PACK
      USE OCEAN, only : UO,VO
      implicit none
      real*8, dimension(:,:), allocatable ::
     &     ones_band,ocnu_band,ocnv_band
      integer :: jmin,jmax

c If DYNSI grid == ATM grid, simply replicate ATM copy
c (Note this requires that TOC2SST has been called first)
      if(aIM.ne.aJM .and. ! extra check if the atm is lat-lon
     &     aIM.eq.iIM .and. aJM.eq.iJM) then
        uosurf(:,:) = atm_uosurf(:,:)
        vosurf(:,:) = atm_vosurf(:,:)
        return
      endif

c call HNTRP
      if(hntrp_o2i_uv_need_init) then
        call Init_Hntrp_Type(hntrp_o2i_u,
     &       oGRID, .5d0,oDLATM,
     &       iGRID, 0.d0,iDLATM,
     &       0.d0)
        call Init_Hntrp_Type(hntrp_o2i_v,
     &       oGRID, 0.d0,oDLATM,
     &       iGRID, 0.d0,iDLATM,
     &       0.d0,
     &       JMA_4interp=oJM-1) ! secondary ocean lats
        hntrp_o2i_uv_need_init = .false.
      endif
c regrid uosurf
      jmin = hntrp_o2i_u%bpack%jband_strt
      jmax = hntrp_o2i_u%bpack%jband_stop
      ALLOCATE(ocnu_band(oIM,jmin:jmax),ones_band(oIM,jmin:jmax))
      ones_band(:,:) = 1d0
      call BAND_PACK (hntrp_o2i_u%bpack, uo(:,:,1), ocnu_band)
      call HNTR8_band (ones_band, ocnu_band, hntrp_o2i_u, uosurf)
      deallocate(ocnu_band,ones_band)
c regrid vosurf
      jmin = hntrp_o2i_v%bpack%jband_strt
      jmax = hntrp_o2i_v%bpack%jband_stop
      ALLOCATE(ocnv_band(oIM,jmin:jmax),ones_band(oIM,jmin:jmax))
      ones_band(:,:) = 1d0
      call BAND_PACK (hntrp_o2i_v%bpack, vo(:,:,1), ocnv_band)
      call HNTR8_band (ones_band, ocnv_band, hntrp_o2i_v, vosurf)
      deallocate(ocnv_band,ones_band)
      if(igrid%have_north_pole) then
c Eventually, uosurf/vosurf will be interpolated directly to
c the ice B grid, which has no polar point.  For now, setting
c A-grid values at the polar point to zero
        uosurf(:,iJM) = 0.
        vosurf(:,iJM) = 0.
      endif
      return
      END SUBROUTINE OG2IG_uvsurf
