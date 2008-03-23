#include "mpi_defs.h"
module SparseCommunicator_mod
   
   implicit none
#ifdef USE_ESMF
   include 'mpif.h'
#endif
   private

   public :: SparseCommunicator_type
   public :: SparseCommunicator
   public :: clean
   public :: getOffset
   public :: getOffsets
   public :: countInDomain
   public :: gatherIJ
   public :: scatterIJ

   public :: NOT_FOUND

   type SparseCommunicator_type
      private
      integer :: mpiCommmunicator
      integer :: mpiRank
      integer, pointer :: mpiCounts(:)
      integer, pointer :: mpiDispls(:)
      integer, pointer :: localOffsets(:)
      integer, pointer :: globalOffsets(:)
   end type SparseCommunicator_type

   integer, parameter :: NOT_FOUND = -1

   interface gatherIJ
      module procedure gather_IJ
      module procedure gather_IJx
      module procedure gather_IJxx
   end interface

   interface scatterIJ
      module procedure scatter_IJ
      module procedure scatter_IJx
      module procedure scatter_IJxx
   end interface

   interface clean
      module procedure cleanSparseCommunicator
   end interface

!   include 'mpif.h'

   integer, parameter :: ROOT = 0

contains

   function SparseCommunicator(points, lboundLocal, uboundLocal, lboundGlobal, uboundGlobal, comm) result(sparseComm)
      integer, intent(in) :: points(:,:)
      integer, intent(in), optional :: comm
      integer, intent(in) :: lboundLocal(2)
      integer, intent(in) :: uboundLocal(2)
      integer, intent(in) :: lboundGlobal(2)
      integer, intent(in) :: uboundGlobal(2)
      type (SparseCommunicator_type) :: sparseComm

      integer :: numPointsLocal
      integer :: numPointsGlobal
      integer, allocatable :: localOffsets(:)
      integer, allocatable :: globalOffsets(:)
      integer :: p, npes, ier, comm_

#ifdef USE_ESMF
      numPointsLocal  = countInDomain(points, lBoundLocal, uBoundLocal, .true. )
      numPointsGlobal = countInDomain(points, lboundglobal, uboundglobal, .false. )

      comm_ = MPI_COMM_WORLD
      if(present(comm)) comm_ = comm
      sparseComm % mpiCommmunicator = comm_
      call mpi_comm_size(comm_, npes, ier)
      call mpi_comm_rank(comm_, sparseComm % mpiRank, ier)

      allocate(sparseComm % mpiCounts(0: npes-1))
      allocate(sparseComm % mpiDispls(0: npes-1))

      ! determine contribution from each processor
      call mpi_gather(numPointsLocal, 1, MPI_INTEGER, &
           & sparseComm % mpiCounts, 1, MPI_INTEGER,  &
           & ROOT, comm_, ier)

      ! determine offsets in combined list
      sparseComm % mpiDispls(0) = 0
      do p = 1, npes - 1
         sparseComm % mpiDispls(p) = sparseComm % mpiDispls(p-1) + sparseComm % mpiCounts(p-1)
      end do

      allocate(sparseComm % localOffsets(numPointsLocal))

      allocate(localOffsets(size(points,2)))
      localOffsets = getOffsets(points, lboundlocal, uboundlocal, .true.)

      ! only keep the ones that are local - use F90 pack to eliminate others
      sparseComm % localOffsets  = pack(localOffsets, localOffsets /= NOT_FOUND)

      allocate(sparseComm % globalOffsets(numPointsGlobal))
      allocate(globalOffsets(size(sparseComm % localOffsets)))

      ! Rather than force a global ordering of the points, we just gather the
      ! global offsets for local points rather than just computing global offsets
      ! on the root processor.
      globalOffsets = pack(getOffsets(points, lboundGlobal, uboundGlobal, .false.), localOffsets /= NOT_FOUND)

      call mpi_gatherV(globalOffsets, size(globalOffsets), MPI_INTEGER, &
           & sparseComm % globalOffsets, sparseComm % mpiCounts, sparseComm % mpiDispls, MPI_INTEGER, &
           & ROOT, comm_, ier)
#else
      npes=1
      sparseComm % mpiCommmunicator = 0
      sparseComm % mpiRank          = 0
      allocate(sparseComm % mpiCounts(0: npes-1))
      allocate(sparseComm % mpiDispls(0: npes-1))

      numPointsLocal  = countInDomain(points, lBoundLocal, uBoundLocal, .true. )
      numPointsGlobal = countInDomain(points, lboundglobal, uboundglobal, .false. )
      if(numPointsLocal.ne.numPointsGlobal) write(*,*)'error mkb'

      sparseComm % mpiCounts(0) = numPointsLocal
      sparseComm % mpiDispls(0) = 0

      allocate(sparseComm % localOffsets(numPointsLocal))

      allocate(localOffsets(size(points,2)))
      localOffsets = getOffsets(points, lboundlocal, uboundlocal, .true.)

      ! only keep the ones that are local - use F90 pack to eliminate others
      sparseComm % localOffsets  = pack(localOffsets, localOffsets /= NOT_FOUND)

      allocate(sparseComm % globalOffsets(numPointsGlobal))
      allocate(globalOffsets(size(sparseComm % localOffsets)))

      ! Rather than force a global ordering of the points, we just gather the
      ! global offsets for local points rather than just computing global offsets
      ! on the root processor.
      sparseComm % globalOffsets = &
        & pack(getOffsets(points, lboundGlobal, uboundGlobal, .false.), localOffsets /= NOT_FOUND)
#endif

      deallocate(globalOffsets)
      deallocate(localOffsets)

   end function SparseCommunicator

   integer function getOffset(point, lowerBounds, upperBounds, addHalo) result(offset)
      integer, intent(in) :: point(2)
      integer, intent(in) :: lowerBounds(2)
      integer, intent(in) :: upperBounds(2)
      logical, intent(in) :: addHalo

      logical :: inBoundsI, inBoundsJ
      integer, parameter :: HALO_WIDTH = 1

      inBoundsI = (point(1) >= lowerBounds(1) .and. point(1) <= upperBounds(1))
      if ( addHalo ) then
        inBoundsJ = (point(2) >= lowerBounds(2)+HALO_WIDTH .and. point(2) <= upperBounds(2)-HALO_WIDTH)
      else
        inBoundsJ = (point(2) >= lowerBounds(2)            .and. point(2) <= upperBounds(2)           )
      endif

      if (inBoundsI .and. inBoundsJ) then
         offset = (point(1) - lowerBounds(1))  + (point(2) - lowerBounds(2)) * (upperBounds(1) - lowerBounds(1) + 1)
      else
         offset = NOT_FOUND
      end if

   end function getOffset

   function getOffsets(points, lowerBounds, upperBounds, addHalo) result(offsets)
      integer, intent(in) :: points(:,:)  ! i,j in first index
      integer, intent(in) :: lowerBounds(2)
      integer, intent(in) :: upperBounds(2)
      logical, intent(in) :: addHalo

      integer :: offsets(size(points,2))
      integer :: iPoint

      do iPoint = 1, size(points,2)
         offsets(iPoint) = getOffset(points(:,iPoint), lowerBounds, upperBounds, addHalo)
      end do

   end function getOffsets

   integer function countInDomain(points, lowerBounds, upperBounds, addHalo)
      integer, intent(in) :: points(:,:)  ! i,j in first index
      integer, intent(in) :: lowerBounds(2)
      integer, intent(in) :: upperBounds(2)
      logical, intent(in) :: addHalo

      integer :: offsets(size(points,2))
      integer :: iPoint

      countInDomain = count (getOffsets(points, lowerBounds, upperBounds, addHalo) /= NOT_FOUND)

   end function countInDomain

   subroutine cleanSparseCommunicator(this)
      type (SparseCommunicator_type) :: this

      deallocate(this % globalOffsets)
      deallocate(this % localOffsets)
      deallocate(this % mpiDispls)
      deallocate(this % mpiCounts)

   end subroutine cleanSparseCommunicator

   logical function amRoot(this)
      type (SparseCommunicator_type) :: this
      amRoot = (this % mpiRank == ROOT)
   end function amRoot

   subroutine gather(this, localArray, globalArray, numBlks, lstride, gstride)
      type (SparseCommunicator_type) :: this
      REAL*8, intent(in)    :: localArray(*)
      REAL*8, intent(inOut) :: globalArray(*)
      integer, intent(in)    :: numBlks, lstride, gstride

      integer :: ier, i, j, bloc, loc

      REAL*8, allocatable :: localBuffer(:)
      REAL*8, allocatable :: globalBuffer(:)

      allocate(localBuffer(size(this % localOffsets)*numBlks))
      allocate(globalBuffer(size(this % globalOffsets)*numBlks))
      
      !pack
      if(size(this % localOffsets).gt.0)then
         bloc=1
         do i = 1, size(this % localOffsets)
            loc = 1 + this % localOffsets(i)
            do j = 1, numBlks
               localBuffer(bloc)=localArray(loc)
               loc = loc+lstride
               bloc = bloc+1
            end do
         end do
      endif

#ifdef USE_ESMF
      call mpi_GatherV(localBuffer, size(localBuffer), MPI_DOUBLE_PRECISION, &
           & globalBuffer, this % mpiCounts*numBlks, this % mpiDispls*numBlks, MPI_DOUBLE_PRECISION, &
           & ROOT, this % mpiCommmunicator, ier)
#else
      globalBuffer(:) = localBuffer(:)
#endif

      !unpack at root
      if (amRoot(this)) then
         bloc=1
         do i = 1, size(this % globalOffsets)
            loc = 1 + this % globalOffsets(i)
            do j = 1, numBlks
               globalArray(loc) = globalBuffer(bloc)
               loc = loc + gstride
               bloc = bloc + 1
            end do
         end do
      endif

      deallocate(globalBuffer)
      deallocate(localBuffer)

   end subroutine gather

   subroutine scatter(this, localArray, globalArray, numBlks, lstride, gstride)
      type (SparseCommunicator_type) :: this
      REAL*8, intent(inOut)    :: localArray(*)
      REAL*8, intent(in) :: globalArray(*)
      integer, intent(in)    :: numBlks, lstride, gstride

      integer :: ier, i, j, bloc, loc

      REAL*8, allocatable :: localBuffer(:)
      REAL*8, allocatable :: globalBuffer(:)

      allocate(localBuffer(size(this % localOffsets)*numBlks))
      allocate(globalBuffer(size(this % globalOffsets)*numBlks))

      !pack at root
      if (amRoot(this)) then
         bloc=1
         do i = 1, size(this % globalOffsets)
            loc = 1 + this % globalOffsets(i)
            do j = 1, numBlks
               globalBuffer(bloc) = globalArray(loc)
               loc = loc + gstride
               bloc = bloc + 1
            end do
         end do
      endif

#ifdef USE_ESMF
      call mpi_ScatterV(globalBuffer, this % mpiCounts*numBlks, this % mpiDispls*numBlks,  &
           &  MPI_DOUBLE_PRECISION, localBuffer, size(localBuffer),  MPI_DOUBLE_PRECISION, &
           & ROOT, this % mpiCommmunicator, ier)
#else
      localBuffer(:) = globalBuffer(:)
#endif
      
      !unpack
      if(size(this % localOffsets).gt.0)then
         bloc=1
         do i = 1, size(this % localOffsets)
            loc = 1 + this % localOffsets(i)
            do j = 1, numBlks
               localArray(loc) = localBuffer(bloc)
               loc = loc+lstride
               bloc = bloc+1
            end do
         end do
      endif

      deallocate(globalBuffer)
      deallocate(localBuffer)

   end subroutine scatter

   SUBROUTINE  gather_IJ(this,ARR,ARR_GLOB)
   IMPLICIT NONE
   type (SparseCommunicator_type), intent(in) :: this
   REAL*8, INTENT(IN)  :: ARR(:,:)
   REAL*8, INTENT(INOUT) :: ARR_GLOB(:,:)
   !local variables
   Integer :: arrRank, numBlks, lstride, gstride

   arrRank=size(shape(ARR_GLOB))
   numBlks=1
   lstride=0
   gstride=0

   call gather(this, ARR, ARR_GLOB, numBlks, lstride, gstride)

   RETURN
   END SUBROUTINE gather_IJ

   SUBROUTINE  gather_IJx(this,ARR,ARR_GLOB)
   IMPLICIT NONE
   type (SparseCommunicator_type), intent(in) :: this
   REAL*8, INTENT(IN)  :: ARR(:,:,:)
   REAL*8, INTENT(INOUT) :: ARR_GLOB(:,:,:)
   !local variables
   Integer :: arrRank, numBlks, lstride, gstride

   arrRank=size(shape(ARR_GLOB))
   numBlks=size(arr_glob,3)
   lstride=ubound(ARR_GLOB,1)*(ubound(ARR,2)-lbound(ARR,2)+1)
   gstride=ubound(ARR_GLOB,1)*ubound(ARR_GLOB,2)

   call gather(this, ARR, ARR_GLOB, numBlks, lstride, gstride)

   RETURN
   END SUBROUTINE gather_IJx

   SUBROUTINE  gather_IJxx(this,ARR,ARR_GLOB)
   IMPLICIT NONE
   type (SparseCommunicator_type), intent(in) :: this
   REAL*8, INTENT(IN)  :: ARR(:,:,:,:)
   REAL*8, INTENT(INOUT) :: ARR_GLOB(:,:,:,:)
   !local variables
   Integer :: arrRank, numBlks, lstride, gstride

   arrRank=size(shape(ARR_GLOB))
   numBlks=size(arr_glob,3)*size(arr_glob,4)
   lstride=ubound(ARR_GLOB,1)*(ubound(ARR,2)-lbound(ARR,2)+1)
   gstride=ubound(ARR_GLOB,1)*ubound(ARR_GLOB,2)

   call gather(this, ARR, ARR_GLOB, numBlks, lstride, gstride)

   RETURN
   END SUBROUTINE gather_IJxx

   SUBROUTINE  scatter_IJ(this,ARR_GLOB,ARR)
   IMPLICIT NONE
   type (SparseCommunicator_type), intent(in) :: this
   REAL*8, INTENT(IN) :: ARR_GLOB(:,:)
   REAL*8, INTENT(INOUT)  :: ARR(:,:)
   !local variables
   Integer :: arrRank, numBlks, lstride, gstride

   arrRank=size(shape(ARR_GLOB))
   numBlks=1
   lstride=0
   gstride=0

   call scatter(this, ARR, ARR_GLOB, numBlks, lstride, gstride)

   RETURN
   END SUBROUTINE scatter_IJ

   SUBROUTINE  scatter_IJx(this,ARR_GLOB,ARR)
   IMPLICIT NONE
   type (SparseCommunicator_type), intent(in) :: this
   REAL*8, INTENT(IN) :: ARR_GLOB(:,:,:)
   REAL*8, INTENT(INOUT)  :: ARR(:,:,:)
   !local variables
   Integer :: arrRank, numBlks, lstride, gstride

   arrRank=size(shape(ARR_GLOB))
   numBlks=size(arr_glob,3)
   lstride=ubound(ARR_GLOB,1)*(ubound(ARR,2)-lbound(ARR,2)+1)
   gstride=ubound(ARR_GLOB,1)*ubound(ARR_GLOB,2)

   call scatter(this, ARR, ARR_GLOB, numBlks, lstride, gstride)

   RETURN
   END SUBROUTINE scatter_IJx

   SUBROUTINE  scatter_IJxx(this,ARR_GLOB,ARR)
   IMPLICIT NONE
   type (SparseCommunicator_type), intent(in) :: this
   REAL*8, INTENT(IN) :: ARR_GLOB(:,:,:,:)
   REAL*8, INTENT(INOUT)  :: ARR(:,:,:,:)
   !local variables
   Integer :: arrRank, numBlks, lstride, gstride

   arrRank=size(shape(ARR_GLOB))
   numBlks=size(arr_glob,3)*size(arr_glob,4)
   lstride=ubound(ARR_GLOB,1)*(ubound(ARR,2)-lbound(ARR,2)+1)
   gstride=ubound(ARR_GLOB,1)*ubound(ARR_GLOB,2)

   call scatter(this, ARR, ARR_GLOB, numBlks, lstride, gstride)

   RETURN
   END SUBROUTINE scatter_IJxx

end module SparseCommunicator_mod
