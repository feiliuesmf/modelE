!***********************************************************
!* WARNING                                                 *
!*                                                         *
!* This file has been automatically generated via a rather *
!* perverse set of preprocessors.  Do _NOT_ directly edit  *
!* this file.  Instead, on should obtain the $?.nw (noweb) *
!* and make changes there.                                 *
!*                                                         *
!* For more information, contact:                          *
!*        Tom Clune <Thomas.L.Clune@nasa.gov>              *
!*                                                         *
!***********************************************************
#include "assert.h"

Module FieldComm
  Implicit None
  Private

  Public :: Field_Halo
  !Public :: Field_Redist
  Public :: Field_Reduce
  Public :: Field_Gather
  !Public :: Field_Scatter

  Integer, Parameter, Public :: FIELD_SUM = 1
  Integer, Parameter, Public :: FIELD_MIN = 2
  Integer, Parameter, Public :: FIELD_MAX = 3
  Integer, Parameter, Public :: FIELD_REDUCE_OPS(3) = (/ FIELD_SUM, FIELD_MIN, FIELD_MAX /)

#ifdef USE_ESMF
  include 'mpif.h'
#endif

  Integer, Parameter :: J_GRID = 2
  Integer, Parameter :: halowidth = 1


  Interface Field_Reduce
     Module Procedure field_reduce_Real4_D
     Module Procedure field_reduce_Real8_D
     Module Procedure field_reduce_Integer_D
  End Interface


Contains

  Subroutine Field_Halo(field, direction, rc)
    Use ESMF_MOD, Only: ESMF_Field
    Use ESMF_MOD, Only: ESMF_DataType
    Use ESMF_MOD, Only: ESMF_DataKind
    Use ESMF_MOD, Only: ESMF_SUCCESS
#ifdef USE_ESMF
    USE ESMF_MOD ! only way to trick absoft into allowing "==" to be brought from ESMF_Base
    Use ESMF_MOD, Only: ESMF_FieldGet
    Use ESMF_MOD, Only: ESMF_FieldGetArray
    Use ESMF_MOD, Only: ESMF_Array
    Use ESMF_MOD, Only: ESMF_ArrayCreate
    Use ESMF_MOD, Only: ESMF_ArrayGetData
    Use ESMF_MOD, Only: ESMF_ArrayF90Allocate
    USE ESMF_Mod, Only: ESMF_KIND_R4
    USE ESMF_Mod, Only: ESMF_KIND_R8
    USE ESMF_Mod, Only: ESMF_DATA_REAL
    Use ESMF_MOD, Only: ESMF_FieldDataMap
    Use ESMF_MOD, Only: ESMF_FieldDataMapGet
    Use ESMF_MOD, Only: ESMF_FieldGetDataPointer
    Use ESMF_MOD, Only: ESMF_DATA_REF
    Use ESMF_MOD_private, ONLY : NORTH, SOUTH
#endif
  ! Arguments
    Type (ESMF_Field), Intent(InOut) :: field
    Integer, Intent(In) :: direction
    Integer, Optional, Intent(Out) :: rc
    Type (ESMF_DataType) :: type
    Type (ESMF_DataKind) :: kind
#ifdef USE_ESMF
    ! Local variables
    Type (ESMF_Array) :: f_array
    Integer :: mpi_type
    Integer :: comm, root_pet = 0, npes, mpi_rank, root, pet
    Integer :: ier
    Integer :: rank, dist_idx
    Integer, Allocatable :: lb(:), ub(:) ! interior bounds
    Integer, Allocatable :: lbh(:), ubh(:) ! halo bounds
    Integer, Allocatable :: loc_cnt(:)

    Type (ESMF_FieldDataMap) :: map

    Logical :: from_south, from_north
    Integer :: pe_south, pe_north
    Integer :: stag, rtag
    Integer :: status(MPI_STATUS_SIZE)
    
    Real (Kind=ESMF_KIND_R4), Pointer :: f_ptr_Real4_3D (:,:,:) 

    
    Real (Kind=ESMF_KIND_R8), Pointer :: f_ptr_Real8_1D (:) 

    
    Real (Kind=ESMF_KIND_R8), Pointer :: f_ptr_Real8_2D (:,:) 

    
    Real (Kind=ESMF_KIND_R8), Pointer :: f_ptr_Real8_3D (:,:,:) 

    
    Real (Kind=ESMF_KIND_R8), Pointer :: f_ptr_Real8_4D (:,:,:,:) 

    
    Real (Kind=ESMF_KIND_R8), Pointer :: f_ptr_Real8_5D (:,:,:,:,:) 

    
    Integer , Pointer :: f_ptr_Integer_1D (:) 

    
    Integer , Pointer :: f_ptr_Integer_2D (:,:) 

    
    Integer , Pointer :: f_ptr_Integer_3D (:,:,:) 

    comm = Field_Get_MPI_comm(field, comm_size=npes, comm_rank=mpi_rank, comm_root=root, rc=rc)
    ASSERT(rc == ESMF_SUCCESS,'find_index')

    ! Determine rank of field array.
    !-------------------------------
    Call ESMF_FieldGet(field, array=f_array, rc=rc)
    ASSERT(rc == ESMF_SUCCESS,'fieldget')
    Call ESMF_ArrayGet(f_array, rank=rank, type=type, kind=kind, rc=rc)
    ASSERT(rc == ESMF_SUCCESS,'arrayget')

    ! create a new mpi type for use in communication
    !-------------------------------
    Allocate(loc_cnt(rank), STAT=ier)
    Call Field_AnalyzeCommPattern(field, mpi_type, dist_idx=dist_idx, counts=loc_cnt, rc=rc)
    ASSERT(rc == ESMF_SUCCESS,'analyze')
    Deallocate(loc_cnt, STAT=ier)

    ! Determine neigboring processes
    !-------------------------------
    Call GetNeighbors(mpi_rank, npes, pe_south, pe_north)

    Allocate(lb(rank), ub(rank), lbh(rank), ubh(rank), STAT = ier)

    ! Define a template for haloing various data patterns (types, kinds,
    ! and ranks). A generic pointer would make this only need a single
    ! instantiation, but this mechanism should be more portable.
    ! -------------------------------
    
     
    If (rank == 3 .and. type == ESMF_DATA_REAL .and. kind == ESMF_R4) Then 
       Call ESMF_FieldGetDataPointer(field, ptr=f_ptr_Real4_3D, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'getdata pointer') 
 
       lbh = LBound( f_ptr_Real4_3D ) 
       ubh = lbh; ubh(dist_idx) = Ubound( f_ptr_Real4_3D, DIM=dist_idx ) 
 
       lb  = lbh; lb(dist_idx) = lb(dist_idx) + haloWidth 
       ub  = ubh; ub(dist_idx) = ub(dist_idx) - haloWidth 
        
       from_north = (IAND(direction, NORTH) > 0) 
       from_south = (IAND(direction, SOUTH) > 0) 
 
       If (from_north) Then 
          stag = 100 + mpi_rank 
          rtag = 100 + pe_north 
 
          Call MPI_Sendrecv( f_ptr_Real4_3D(lb(1),lb(2),lb(3)),  1, mpi_type, pe_south, stag, & 
               &             f_ptr_Real4_3D(ubh(1),ubh(2),ubh(3)), 1, mpi_type, pe_north, rtag, & 
               &             comm, status, ier) 
       End If 
 
       If (from_south) Then 
          stag = 100 + mpi_rank 
          rtag = 100 + pe_south 
 
          Call MPI_Sendrecv( f_ptr_Real4_3D(ub(1),ub(2),ub(3)),  1, mpi_type, pe_north, stag, & 
               &             f_ptr_Real4_3D(lbh(1),lbh(2),lbh(3)), 1, mpi_type, pe_south, rtag, & 
               &             comm, status, ier) 
       End If 
 
       Nullify(f_ptr_Real4_3D) 
    End If 
     

    
     
    If (rank == 1 .and. type == ESMF_DATA_REAL .and. kind == ESMF_R8) Then 
       Call ESMF_FieldGetDataPointer(field, ptr=f_ptr_Real8_1D, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'getdata pointer') 
 
       lbh = LBound( f_ptr_Real8_1D ) 
       ubh = lbh; ubh(dist_idx) = Ubound( f_ptr_Real8_1D, DIM=dist_idx ) 
 
       lb  = lbh; lb(dist_idx) = lb(dist_idx) + haloWidth 
       ub  = ubh; ub(dist_idx) = ub(dist_idx) - haloWidth 
        
       from_north = (IAND(direction, NORTH) > 0) 
       from_south = (IAND(direction, SOUTH) > 0) 
 
       If (from_north) Then 
          stag = 100 + mpi_rank 
          rtag = 100 + pe_north 
 
          Call MPI_Sendrecv( f_ptr_Real8_1D(lb(1)),  1, mpi_type, pe_south, stag, & 
               &             f_ptr_Real8_1D(ubh(1)), 1, mpi_type, pe_north, rtag, & 
               &             comm, status, ier) 
       End If 
 
       If (from_south) Then 
          stag = 100 + mpi_rank 
          rtag = 100 + pe_south 
 
          Call MPI_Sendrecv( f_ptr_Real8_1D(ub(1)),  1, mpi_type, pe_north, stag, & 
               &             f_ptr_Real8_1D(lbh(1)), 1, mpi_type, pe_south, rtag, & 
               &             comm, status, ier) 
       End If 
 
       Nullify(f_ptr_Real8_1D) 
    End If 
     

    
     
    If (rank == 2 .and. type == ESMF_DATA_REAL .and. kind == ESMF_R8) Then 
       Call ESMF_FieldGetDataPointer(field, ptr=f_ptr_Real8_2D, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'getdata pointer') 
 
       lbh = LBound( f_ptr_Real8_2D ) 
       ubh = lbh; ubh(dist_idx) = Ubound( f_ptr_Real8_2D, DIM=dist_idx ) 
 
       lb  = lbh; lb(dist_idx) = lb(dist_idx) + haloWidth 
       ub  = ubh; ub(dist_idx) = ub(dist_idx) - haloWidth 
        
       from_north = (IAND(direction, NORTH) > 0) 
       from_south = (IAND(direction, SOUTH) > 0) 
 
       If (from_north) Then 
          stag = 100 + mpi_rank 
          rtag = 100 + pe_north 
 
          Call MPI_Sendrecv( f_ptr_Real8_2D(lb(1),lb(2)),  1, mpi_type, pe_south, stag, & 
               &             f_ptr_Real8_2D(ubh(1),ubh(2)), 1, mpi_type, pe_north, rtag, & 
               &             comm, status, ier) 
       End If 
 
       If (from_south) Then 
          stag = 100 + mpi_rank 
          rtag = 100 + pe_south 
 
          Call MPI_Sendrecv( f_ptr_Real8_2D(ub(1),ub(2)),  1, mpi_type, pe_north, stag, & 
               &             f_ptr_Real8_2D(lbh(1),lbh(2)), 1, mpi_type, pe_south, rtag, & 
               &             comm, status, ier) 
       End If 
 
       Nullify(f_ptr_Real8_2D) 
    End If 
     

    
     
    If (rank == 3 .and. type == ESMF_DATA_REAL .and. kind == ESMF_R8) Then 
       Call ESMF_FieldGetDataPointer(field, ptr=f_ptr_Real8_3D, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'getdata pointer') 
 
       lbh = LBound( f_ptr_Real8_3D ) 
       ubh = lbh; ubh(dist_idx) = Ubound( f_ptr_Real8_3D, DIM=dist_idx ) 
 
       lb  = lbh; lb(dist_idx) = lb(dist_idx) + haloWidth 
       ub  = ubh; ub(dist_idx) = ub(dist_idx) - haloWidth 
        
       from_north = (IAND(direction, NORTH) > 0) 
       from_south = (IAND(direction, SOUTH) > 0) 
 
       If (from_north) Then 
          stag = 100 + mpi_rank 
          rtag = 100 + pe_north 
 
          Call MPI_Sendrecv( f_ptr_Real8_3D(lb(1),lb(2),lb(3)),  1, mpi_type, pe_south, stag, & 
               &             f_ptr_Real8_3D(ubh(1),ubh(2),ubh(3)), 1, mpi_type, pe_north, rtag, & 
               &             comm, status, ier) 
       End If 
 
       If (from_south) Then 
          stag = 100 + mpi_rank 
          rtag = 100 + pe_south 
 
          Call MPI_Sendrecv( f_ptr_Real8_3D(ub(1),ub(2),ub(3)),  1, mpi_type, pe_north, stag, & 
               &             f_ptr_Real8_3D(lbh(1),lbh(2),lbh(3)), 1, mpi_type, pe_south, rtag, & 
               &             comm, status, ier) 
       End If 
 
       Nullify(f_ptr_Real8_3D) 
    End If 
     

    
     
    If (rank == 4 .and. type == ESMF_DATA_REAL .and. kind == ESMF_R8) Then 
       Call ESMF_FieldGetDataPointer(field, ptr=f_ptr_Real8_4D, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'getdata pointer') 
 
       lbh = LBound( f_ptr_Real8_4D ) 
       ubh = lbh; ubh(dist_idx) = Ubound( f_ptr_Real8_4D, DIM=dist_idx ) 
 
       lb  = lbh; lb(dist_idx) = lb(dist_idx) + haloWidth 
       ub  = ubh; ub(dist_idx) = ub(dist_idx) - haloWidth 
        
       from_north = (IAND(direction, NORTH) > 0) 
       from_south = (IAND(direction, SOUTH) > 0) 
 
       If (from_north) Then 
          stag = 100 + mpi_rank 
          rtag = 100 + pe_north 
 
          Call MPI_Sendrecv( f_ptr_Real8_4D(lb(1),lb(2),lb(3),lb(4)),  1, mpi_type, pe_south, stag, & 
               &             f_ptr_Real8_4D(ubh(1),ubh(2),ubh(3),ubh(4)), 1, mpi_type, pe_north, rtag, & 
               &             comm, status, ier) 
       End If 
 
       If (from_south) Then 
          stag = 100 + mpi_rank 
          rtag = 100 + pe_south 
 
          Call MPI_Sendrecv( f_ptr_Real8_4D(ub(1),ub(2),ub(3),ub(4)),  1, mpi_type, pe_north, stag, & 
               &             f_ptr_Real8_4D(lbh(1),lbh(2),lbh(3),lbh(4)), 1, mpi_type, pe_south, rtag, & 
               &             comm, status, ier) 
       End If 
 
       Nullify(f_ptr_Real8_4D) 
    End If 
     

    
     
    If (rank == 5 .and. type == ESMF_DATA_REAL .and. kind == ESMF_R8) Then 
       Call ESMF_FieldGetDataPointer(field, ptr=f_ptr_Real8_5D, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'getdata pointer') 
 
       lbh = LBound( f_ptr_Real8_5D ) 
       ubh = lbh; ubh(dist_idx) = Ubound( f_ptr_Real8_5D, DIM=dist_idx ) 
 
       lb  = lbh; lb(dist_idx) = lb(dist_idx) + haloWidth 
       ub  = ubh; ub(dist_idx) = ub(dist_idx) - haloWidth 
        
       from_north = (IAND(direction, NORTH) > 0) 
       from_south = (IAND(direction, SOUTH) > 0) 
 
       If (from_north) Then 
          stag = 100 + mpi_rank 
          rtag = 100 + pe_north 
 
          Call MPI_Sendrecv( f_ptr_Real8_5D(lb(1),lb(2),lb(3),lb(4),lb(5)),  1, mpi_type, pe_south, stag, & 
               &             f_ptr_Real8_5D(ubh(1),ubh(2),ubh(3),ubh(4),ubh(5)), 1, mpi_type, pe_north, rtag, & 
               &             comm, status, ier) 
       End If 
 
       If (from_south) Then 
          stag = 100 + mpi_rank 
          rtag = 100 + pe_south 
 
          Call MPI_Sendrecv( f_ptr_Real8_5D(ub(1),ub(2),ub(3),ub(4),ub(5)),  1, mpi_type, pe_north, stag, & 
               &             f_ptr_Real8_5D(lbh(1),lbh(2),lbh(3),lbh(4),lbh(5)), 1, mpi_type, pe_south, rtag, & 
               &             comm, status, ier) 
       End If 
 
       Nullify(f_ptr_Real8_5D) 
    End If 
     

    
     
    If (rank == 1 .and. type == ESMF_DATA_INTEGER .and. kind == ESMF_I4) Then 
       Call ESMF_FieldGetDataPointer(field, ptr=f_ptr_Integer_1D, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'getdata pointer') 
 
       lbh = LBound( f_ptr_Integer_1D ) 
       ubh = lbh; ubh(dist_idx) = Ubound( f_ptr_Integer_1D, DIM=dist_idx ) 
 
       lb  = lbh; lb(dist_idx) = lb(dist_idx) + haloWidth 
       ub  = ubh; ub(dist_idx) = ub(dist_idx) - haloWidth 
        
       from_north = (IAND(direction, NORTH) > 0) 
       from_south = (IAND(direction, SOUTH) > 0) 
 
       If (from_north) Then 
          stag = 100 + mpi_rank 
          rtag = 100 + pe_north 
 
          Call MPI_Sendrecv( f_ptr_Integer_1D(lb(1)),  1, mpi_type, pe_south, stag, & 
               &             f_ptr_Integer_1D(ubh(1)), 1, mpi_type, pe_north, rtag, & 
               &             comm, status, ier) 
       End If 
 
       If (from_south) Then 
          stag = 100 + mpi_rank 
          rtag = 100 + pe_south 
 
          Call MPI_Sendrecv( f_ptr_Integer_1D(ub(1)),  1, mpi_type, pe_north, stag, & 
               &             f_ptr_Integer_1D(lbh(1)), 1, mpi_type, pe_south, rtag, & 
               &             comm, status, ier) 
       End If 
 
       Nullify(f_ptr_Integer_1D) 
    End If 
     

    
     
    If (rank == 2 .and. type == ESMF_DATA_INTEGER .and. kind == ESMF_I4) Then 
       Call ESMF_FieldGetDataPointer(field, ptr=f_ptr_Integer_2D, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'getdata pointer') 
 
       lbh = LBound( f_ptr_Integer_2D ) 
       ubh = lbh; ubh(dist_idx) = Ubound( f_ptr_Integer_2D, DIM=dist_idx ) 
 
       lb  = lbh; lb(dist_idx) = lb(dist_idx) + haloWidth 
       ub  = ubh; ub(dist_idx) = ub(dist_idx) - haloWidth 
        
       from_north = (IAND(direction, NORTH) > 0) 
       from_south = (IAND(direction, SOUTH) > 0) 
 
       If (from_north) Then 
          stag = 100 + mpi_rank 
          rtag = 100 + pe_north 
 
          Call MPI_Sendrecv( f_ptr_Integer_2D(lb(1),lb(2)),  1, mpi_type, pe_south, stag, & 
               &             f_ptr_Integer_2D(ubh(1),ubh(2)), 1, mpi_type, pe_north, rtag, & 
               &             comm, status, ier) 
       End If 
 
       If (from_south) Then 
          stag = 100 + mpi_rank 
          rtag = 100 + pe_south 
 
          Call MPI_Sendrecv( f_ptr_Integer_2D(ub(1),ub(2)),  1, mpi_type, pe_north, stag, & 
               &             f_ptr_Integer_2D(lbh(1),lbh(2)), 1, mpi_type, pe_south, rtag, & 
               &             comm, status, ier) 
       End If 
 
       Nullify(f_ptr_Integer_2D) 
    End If 
     

    
     
    If (rank == 3 .and. type == ESMF_DATA_INTEGER .and. kind == ESMF_I4) Then 
       Call ESMF_FieldGetDataPointer(field, ptr=f_ptr_Integer_3D, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'getdata pointer') 
 
       lbh = LBound( f_ptr_Integer_3D ) 
       ubh = lbh; ubh(dist_idx) = Ubound( f_ptr_Integer_3D, DIM=dist_idx ) 
 
       lb  = lbh; lb(dist_idx) = lb(dist_idx) + haloWidth 
       ub  = ubh; ub(dist_idx) = ub(dist_idx) - haloWidth 
        
       from_north = (IAND(direction, NORTH) > 0) 
       from_south = (IAND(direction, SOUTH) > 0) 
 
       If (from_north) Then 
          stag = 100 + mpi_rank 
          rtag = 100 + pe_north 
 
          Call MPI_Sendrecv( f_ptr_Integer_3D(lb(1),lb(2),lb(3)),  1, mpi_type, pe_south, stag, & 
               &             f_ptr_Integer_3D(ubh(1),ubh(2),ubh(3)), 1, mpi_type, pe_north, rtag, & 
               &             comm, status, ier) 
       End If 
 
       If (from_south) Then 
          stag = 100 + mpi_rank 
          rtag = 100 + pe_south 
 
          Call MPI_Sendrecv( f_ptr_Integer_3D(ub(1),ub(2),ub(3)),  1, mpi_type, pe_north, stag, & 
               &             f_ptr_Integer_3D(lbh(1),lbh(2),lbh(3)), 1, mpi_type, pe_south, rtag, & 
               &             comm, status, ier) 
       End If 
 
       Nullify(f_ptr_Integer_3D) 
    End If 
     

    Deallocate(lb, ub, lbh, ubh)
    Call MPI_Type_Free(mpi_type, ier)

#else
  ! do nothing if no ESMF
#endif
    If (Present(rc)) rc = ESMF_SUCCESS
  End Subroutine Field_Halo

  !< [[Field_Redist]] >>
  
  Subroutine Field_Reduce_Real4_D (field, result, rtype, hemi_result, AI_mask, broadcast, rc) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Use ESMF_MOD, Only: ESMF_Field 
    Use ESMF_MOD, Only: ESMF_AxisIndex 
    Use ESMF_MOD, Only: ESMF_FieldGetDataPointer 
    Use ESMF_MOD, Only: ESMF_SUCCESS 
     
  !*********************** 
  ! Arguments 
    Type (ESMF_Field),               Intent(In)  :: field 
    Real (Kind=ESMF_KIND_R4),                   Intent(Out) :: result 
    Integer,                         Intent(In)  :: rtype              ! Sum, Max, Min, ... 
    Real (Kind=ESMF_KIND_R4),         Optional, Intent(Out) :: hemi_result(2)     ! Reduce on each hemisphere 
    Type (ESMF_AxisIndex), Optional, Intent(In)  :: AI_Mask(:)         ! restrict reduction to sub domain 
    Logical,               Optional, Intent(In)  :: broadcast          ! Default puts result only on root PET 
    Integer,               Optional, Intent(Out) :: rc 
 
  ! Local variables 
    Real (Kind=ESMF_KIND_R4), Pointer :: fp_5d(:,:,:,:,:)  ! pointer to field data for each rank 
    Real (Kind=ESMF_KIND_R4), Pointer :: fp_4d(:,:,:,:) 
    Real (Kind=ESMF_KIND_R4), Pointer :: fp_3d(:,:,:) 
    Real (Kind=ESMF_KIND_R4), Pointer :: fp_2d(:,:) 
    Real (Kind=ESMF_KIND_R4), Pointer :: fp_1d(:) 
 
    Real (Kind=ESMF_KIND_R4), Allocatable :: x(:), psum(:) ! temporary buffer 
 
    Integer :: dist_idx 
    Integer :: rank_result 
    Integer :: ni, nj, nx1, nx2 
    Integer, Allocatable :: shp_mask(:), lb(:), ub(:) 
    Integer, Allocatable :: cnts(:), dspls(:) 
    Integer :: n_loc, n_glob, lb_glob, ub_glob 
    Integer :: field_rank, rank 
    Integer :: ldx, m, n 
    Integer :: comm, root_pet, mype, npes, pe, ier 
    Integer :: eq, eq_sth, eq_nth, j_equator 
    Logical :: field_is_staggered 
    Logical :: bcast 
    Integer, Parameter :: NH = 2, SH = 1 
 
    bcast =.false. 
    IF (Present(broadcast)) bcast = broadcast 
 
    dist_idx = Field_FindIndex(field, J_GRID, rc=rc) 
    ASSERT(rc == ESMF_SUCCESS,'find_index') 
 
    field_rank = Field_GetRank(field, rc=rc) 
    ASSERT(rc == ESMF_SUCCESS,'get rank') 
    Allocate(lb(field_rank), ub(field_rank), shp_mask(field_rank), STAT=ier) 
    ASSERT_ALWAYS(ier==0,'allocate()') 
 
    Call Field_GetBounds(field, lb, ub, rc=rc) 
 
    n = Product (ub(1:field_rank) - lb(1:field_rank) + 1) 
    Allocate(x(n), STAT=ier) 
    ASSERT_ALWAYS(ier== 0,'allocate()') 
 
    If (Present(AI_Mask)) Then 
       ASSERT(Size(AI_Mask) == field_rank, 'inconsistent ranks') 
       lb = Max(lb, AI_Mask%min) 
       ub = Min(ub, AI_Mask%max) 
    End If 
 
    shp_mask = ub - lb +1 
 
    Select Case (field_rank) 
    Case (1) 
       Call ESMF_FieldGetDataPointer(field, fp_1d, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'fp_1d') 
       Call CopyWindow_Real4_D(fp_1d(lb(1):ub(1)), x, n) 
       Nullify(fp_1d) 
    Case (2) 
       Call ESMF_FieldGetDataPointer(field, fp_2d, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'fp_2d') 
       Call CopyWindow_Real4_D(fp_2d(lb(1):ub(1),lb(2):ub(2)), x, n) 
       Nullify(fp_2d) 
    Case (3) 
       Call ESMF_FieldGetDataPointer(field, fp_3d, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'fp_3d') 
       Call CopyWindow_Real4_D(fp_3d(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)), x, n) 
       Nullify(fp_3d) 
    Case (4) 
       Call ESMF_FieldGetDataPointer(field, ptr=fp_4d, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'fp_4d') 
       Call CopyWindow_Real4_D(fp_4d(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)), x, n) 
       Nullify(fp_4d) 
    Case (5) 
       Call ESMF_FieldGetDataPointer(field, fp_5d, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'fp_5d') 
       Call CopyWindow_Real4_D(fp_5d(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5)), x, n) 
       Nullify(fp_5d) 
    End Select 
 
    Do rank = 1, field_rank 
 
       If (rank == dist_idx) Cycle 
 
       ldx = Product(shp_mask(:rank-1)) 
       m   = shp_mask(rank) 
       If (m == 1) cycle 
       n   = Product(shp_mask(rank+1:field_rank)) 
 
       Call LocalReduce_Real4_D(x(1), ldx, m, n, rtype) 
       shp_mask(rank) = 1 ! new shape 
 
    End Do 
 
 
 
    n_loc  = Max(0, shp_mask(dist_idx)) 
 
#ifdef USE_ESMF 
    comm = Field_Get_MPI_comm(field, npes, mype, root_pet) 
 
    Allocate(cnts(0:npes-1), dspls(0:npes-1), STAT=ier) 
    cnts=0 
 
    Call MPI_AllGather(n_loc, 1, MPI_INTEGER, cnts, 1, MPI_INTEGER, comm, ier) 
    Call MPI_AllReduce(lb(dist_idx), lb_glob, 1, MPI_Integer, MPI_MIN, comm, ier) 
    Call MPI_AllReduce(ub(dist_idx), ub_glob, 1, MPI_Integer, MPI_MAX, comm, ier) 
 
    dspls(0) = 0 
    Do pe = 1, npes-1 
       dspls(pe) = dspls(pe-1) + cnts(pe-1) 
    End Do 
 
    n_glob=Sum(cnts) 
    ASSERT(n_glob == ub_glob - lb_glob + 1,'inconsistent') 
    Allocate(psum(lb_glob:ub_glob), STAT=ier) 
 
    If (bcast) Then 
       Call MPI_AllGatherV(x, n_loc, MPI_REAL, psum, cnts, dspls, MPI_REAL, comm, ier) 
    Else 
       Call MPI_GatherV(x, n_loc, MPI_REAL, psum, cnts, dspls, MPI_REAL, root_pet, comm, ier) 
    End If 
#else 
    Allocate(psum(n_loc)) 
    psum = x(1:n_loc) 
#endif 
 
 
    ASSERT_ALWAYS(ANY(rtype == FIELD_REDUCE_OPS),'undefined reduction operation') 
 
    Select Case (rtype) 
    Case (FIELD_SUM) 
       result = Sum(psum) 
    Case (FIELD_MAX) 
       result = MaxVal(psum) 
    Case (FIELD_MIN) 
       result = MinVal(psum) 
    End Select 
 
 
    If (Present(hemi_result)) Then 
       !Call Grid_GetExtents(grid, equator=equator) 
 
       hemi_result=0 
       !!$$If (field_is_staggered) Then 
       !!$$   eq_sth = eq 
       !!$$   eq_nth = eq + 2 
       !!$$Else 
       !!$$   eq_sth = eq 
       !!$$   eq_nth = eq 
       !!$$End If 
       !!$$ 
       !!$$Select Case (rtype) 
       !!$$Case (FIELD_SUM) 
       !!$$   hemi_result(SH) = Sum(psum(:eq_sth)) 
       !!$$   hemi_result(NH) = Sum(psum(eq_nth:)) 
       !!$$   If (field_is_staggered) hemi_result = hemi_result + 0.5 * psum(eq+1) 
       !!$$Case (FIELD_MAX) 
       !!$$   hemi_result(SH) = MaxVal(psum(:eq_sth)) 
       !!$$   hemi_result(NH) = MaxVal(psum(eq_nth:)) 
       !!$$   If (field_is_staggered) hemi_result = Max(hemi_result, psum(eq+1)) 
       !!$$Case (FIELD_MIN) 
       !!$$   hemi_result(SH) = MinVal(psum(:eq_sth)) 
       !!$$   hemi_result(NH) = MinVal(psum(eq_nth:)) 
       !!$$   If (field_is_staggered) hemi_result = Min(hemi_result, psum(eq+1)) 
       !!$$End Select 
 
 
    End If 
 
    Deallocate(psum) 
#ifdef USE_ESMF 
    Deallocate(cnts, dspls) 
#endif 
    Deallocate(lb, ub, shp_mask) 
    Deallocate(x) 
 
  End Subroutine Field_Reduce_Real4_D 

  
  Subroutine Field_Reduce_Real8_D (field, result, rtype, hemi_result, AI_mask, broadcast, rc) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Use ESMF_MOD, Only: ESMF_Field 
    Use ESMF_MOD, Only: ESMF_AxisIndex 
    Use ESMF_MOD, Only: ESMF_FieldGetDataPointer 
    Use ESMF_MOD, Only: ESMF_SUCCESS 
     
  !*********************** 
  ! Arguments 
    Type (ESMF_Field),               Intent(In)  :: field 
    Real (Kind=ESMF_KIND_R8),                   Intent(Out) :: result 
    Integer,                         Intent(In)  :: rtype              ! Sum, Max, Min, ... 
    Real (Kind=ESMF_KIND_R8),         Optional, Intent(Out) :: hemi_result(2)     ! Reduce on each hemisphere 
    Type (ESMF_AxisIndex), Optional, Intent(In)  :: AI_Mask(:)         ! restrict reduction to sub domain 
    Logical,               Optional, Intent(In)  :: broadcast          ! Default puts result only on root PET 
    Integer,               Optional, Intent(Out) :: rc 
 
  ! Local variables 
    Real (Kind=ESMF_KIND_R8), Pointer :: fp_5d(:,:,:,:,:)  ! pointer to field data for each rank 
    Real (Kind=ESMF_KIND_R8), Pointer :: fp_4d(:,:,:,:) 
    Real (Kind=ESMF_KIND_R8), Pointer :: fp_3d(:,:,:) 
    Real (Kind=ESMF_KIND_R8), Pointer :: fp_2d(:,:) 
    Real (Kind=ESMF_KIND_R8), Pointer :: fp_1d(:) 
 
    Real (Kind=ESMF_KIND_R8), Allocatable :: x(:), psum(:) ! temporary buffer 
 
    Integer :: dist_idx 
    Integer :: rank_result 
    Integer :: ni, nj, nx1, nx2 
    Integer, Allocatable :: shp_mask(:), lb(:), ub(:) 
    Integer, Allocatable :: cnts(:), dspls(:) 
    Integer :: n_loc, n_glob, lb_glob, ub_glob 
    Integer :: field_rank, rank 
    Integer :: ldx, m, n 
    Integer :: comm, root_pet, mype, npes, pe, ier 
    Integer :: eq, eq_sth, eq_nth, j_equator 
    Logical :: field_is_staggered 
    Logical :: bcast 
    Integer, Parameter :: NH = 2, SH = 1 
 
    bcast =.false. 
    IF (Present(broadcast)) bcast = broadcast 
 
    dist_idx = Field_FindIndex(field, J_GRID, rc=rc) 
    ASSERT(rc == ESMF_SUCCESS,'find_index') 
 
    field_rank = Field_GetRank(field, rc=rc) 
    ASSERT(rc == ESMF_SUCCESS,'get rank') 
    Allocate(lb(field_rank), ub(field_rank), shp_mask(field_rank), STAT=ier) 
    ASSERT_ALWAYS(ier==0,'allocate()') 
 
    Call Field_GetBounds(field, lb, ub, rc=rc) 
 
    n = Product (ub(1:field_rank) - lb(1:field_rank) + 1) 
    Allocate(x(n), STAT=ier) 
    ASSERT_ALWAYS(ier== 0,'allocate()') 
 
    If (Present(AI_Mask)) Then 
       ASSERT(Size(AI_Mask) == field_rank, 'inconsistent ranks') 
       lb = Max(lb, AI_Mask%min) 
       ub = Min(ub, AI_Mask%max) 
    End If 
 
    shp_mask = ub - lb +1 
 
    Select Case (field_rank) 
    Case (1) 
       Call ESMF_FieldGetDataPointer(field, fp_1d, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'fp_1d') 
       Call CopyWindow_Real8_D(fp_1d(lb(1):ub(1)), x, n) 
       Nullify(fp_1d) 
    Case (2) 
       Call ESMF_FieldGetDataPointer(field, fp_2d, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'fp_2d') 
       Call CopyWindow_Real8_D(fp_2d(lb(1):ub(1),lb(2):ub(2)), x, n) 
       Nullify(fp_2d) 
    Case (3) 
       Call ESMF_FieldGetDataPointer(field, fp_3d, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'fp_3d') 
       Call CopyWindow_Real8_D(fp_3d(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)), x, n) 
       Nullify(fp_3d) 
    Case (4) 
       Call ESMF_FieldGetDataPointer(field, ptr=fp_4d, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'fp_4d') 
       Call CopyWindow_Real8_D(fp_4d(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)), x, n) 
       Nullify(fp_4d) 
    Case (5) 
       Call ESMF_FieldGetDataPointer(field, fp_5d, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'fp_5d') 
       Call CopyWindow_Real8_D(fp_5d(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5)), x, n) 
       Nullify(fp_5d) 
    End Select 
 
    Do rank = 1, field_rank 
 
       If (rank == dist_idx) Cycle 
 
       ldx = Product(shp_mask(:rank-1)) 
       m   = shp_mask(rank) 
       If (m == 1) cycle 
       n   = Product(shp_mask(rank+1:field_rank)) 
 
       Call LocalReduce_Real8_D(x(1), ldx, m, n, rtype) 
       shp_mask(rank) = 1 ! new shape 
 
    End Do 
 
 
 
    n_loc  = Max(0, shp_mask(dist_idx)) 
 
#ifdef USE_ESMF 
    comm = Field_Get_MPI_comm(field, npes, mype, root_pet) 
 
    Allocate(cnts(0:npes-1), dspls(0:npes-1), STAT=ier) 
    cnts=0 
 
    Call MPI_AllGather(n_loc, 1, MPI_INTEGER, cnts, 1, MPI_INTEGER, comm, ier) 
    Call MPI_AllReduce(lb(dist_idx), lb_glob, 1, MPI_Integer, MPI_MIN, comm, ier) 
    Call MPI_AllReduce(ub(dist_idx), ub_glob, 1, MPI_Integer, MPI_MAX, comm, ier) 
 
    dspls(0) = 0 
    Do pe = 1, npes-1 
       dspls(pe) = dspls(pe-1) + cnts(pe-1) 
    End Do 
 
    n_glob=Sum(cnts) 
    ASSERT(n_glob == ub_glob - lb_glob + 1,'inconsistent') 
    Allocate(psum(lb_glob:ub_glob), STAT=ier) 
 
    If (bcast) Then 
       Call MPI_AllGatherV(x, n_loc, MPI_DOUBLE_PRECISION, psum, cnts, dspls, MPI_DOUBLE_PRECISION, comm, ier) 
    Else 
       Call MPI_GatherV(x, n_loc, MPI_DOUBLE_PRECISION, psum, cnts, dspls, MPI_DOUBLE_PRECISION, root_pet, comm, ier) 
    End If 
#else 
    Allocate(psum(n_loc)) 
    psum = x(1:n_loc) 
#endif 
 
 
    ASSERT_ALWAYS(ANY(rtype == FIELD_REDUCE_OPS),'undefined reduction operation') 
 
    Select Case (rtype) 
    Case (FIELD_SUM) 
       result = Sum(psum) 
    Case (FIELD_MAX) 
       result = MaxVal(psum) 
    Case (FIELD_MIN) 
       result = MinVal(psum) 
    End Select 
 
 
    If (Present(hemi_result)) Then 
       !Call Grid_GetExtents(grid, equator=equator) 
 
       hemi_result=0 
       !!$$If (field_is_staggered) Then 
       !!$$   eq_sth = eq 
       !!$$   eq_nth = eq + 2 
       !!$$Else 
       !!$$   eq_sth = eq 
       !!$$   eq_nth = eq 
       !!$$End If 
       !!$$ 
       !!$$Select Case (rtype) 
       !!$$Case (FIELD_SUM) 
       !!$$   hemi_result(SH) = Sum(psum(:eq_sth)) 
       !!$$   hemi_result(NH) = Sum(psum(eq_nth:)) 
       !!$$   If (field_is_staggered) hemi_result = hemi_result + 0.5 * psum(eq+1) 
       !!$$Case (FIELD_MAX) 
       !!$$   hemi_result(SH) = MaxVal(psum(:eq_sth)) 
       !!$$   hemi_result(NH) = MaxVal(psum(eq_nth:)) 
       !!$$   If (field_is_staggered) hemi_result = Max(hemi_result, psum(eq+1)) 
       !!$$Case (FIELD_MIN) 
       !!$$   hemi_result(SH) = MinVal(psum(:eq_sth)) 
       !!$$   hemi_result(NH) = MinVal(psum(eq_nth:)) 
       !!$$   If (field_is_staggered) hemi_result = Min(hemi_result, psum(eq+1)) 
       !!$$End Select 
 
 
    End If 
 
    Deallocate(psum) 
#ifdef USE_ESMF 
    Deallocate(cnts, dspls) 
#endif 
    Deallocate(lb, ub, shp_mask) 
    Deallocate(x) 
 
  End Subroutine Field_Reduce_Real8_D 

  
  Subroutine Field_Reduce_Integer_D (field, result, rtype, hemi_result, AI_mask, broadcast, rc) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Use ESMF_MOD, Only: ESMF_Field 
    Use ESMF_MOD, Only: ESMF_AxisIndex 
    Use ESMF_MOD, Only: ESMF_FieldGetDataPointer 
    Use ESMF_MOD, Only: ESMF_SUCCESS 
     
  !*********************** 
  ! Arguments 
    Type (ESMF_Field),               Intent(In)  :: field 
    Integer ,                   Intent(Out) :: result 
    Integer,                         Intent(In)  :: rtype              ! Sum, Max, Min, ... 
    Integer ,         Optional, Intent(Out) :: hemi_result(2)     ! Reduce on each hemisphere 
    Type (ESMF_AxisIndex), Optional, Intent(In)  :: AI_Mask(:)         ! restrict reduction to sub domain 
    Logical,               Optional, Intent(In)  :: broadcast          ! Default puts result only on root PET 
    Integer,               Optional, Intent(Out) :: rc 
 
  ! Local variables 
    Integer , Pointer :: fp_5d(:,:,:,:,:)  ! pointer to field data for each rank 
    Integer , Pointer :: fp_4d(:,:,:,:) 
    Integer , Pointer :: fp_3d(:,:,:) 
    Integer , Pointer :: fp_2d(:,:) 
    Integer , Pointer :: fp_1d(:) 
 
    Integer , Allocatable :: x(:), psum(:) ! temporary buffer 
 
    Integer :: dist_idx 
    Integer :: rank_result 
    Integer :: ni, nj, nx1, nx2 
    Integer, Allocatable :: shp_mask(:), lb(:), ub(:) 
    Integer, Allocatable :: cnts(:), dspls(:) 
    Integer :: n_loc, n_glob, lb_glob, ub_glob 
    Integer :: field_rank, rank 
    Integer :: ldx, m, n 
    Integer :: comm, root_pet, mype, npes, pe, ier 
    Integer :: eq, eq_sth, eq_nth, j_equator 
    Logical :: field_is_staggered 
    Logical :: bcast 
    Integer, Parameter :: NH = 2, SH = 1 
 
    bcast =.false. 
    IF (Present(broadcast)) bcast = broadcast 
 
    dist_idx = Field_FindIndex(field, J_GRID, rc=rc) 
    ASSERT(rc == ESMF_SUCCESS,'find_index') 
 
    field_rank = Field_GetRank(field, rc=rc) 
    ASSERT(rc == ESMF_SUCCESS,'get rank') 
    Allocate(lb(field_rank), ub(field_rank), shp_mask(field_rank), STAT=ier) 
    ASSERT_ALWAYS(ier==0,'allocate()') 
 
    Call Field_GetBounds(field, lb, ub, rc=rc) 
 
    n = Product (ub(1:field_rank) - lb(1:field_rank) + 1) 
    Allocate(x(n), STAT=ier) 
    ASSERT_ALWAYS(ier== 0,'allocate()') 
 
    If (Present(AI_Mask)) Then 
       ASSERT(Size(AI_Mask) == field_rank, 'inconsistent ranks') 
       lb = Max(lb, AI_Mask%min) 
       ub = Min(ub, AI_Mask%max) 
    End If 
 
    shp_mask = ub - lb +1 
 
    Select Case (field_rank) 
    Case (1) 
       Call ESMF_FieldGetDataPointer(field, fp_1d, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'fp_1d') 
       Call CopyWindow_Integer_D(fp_1d(lb(1):ub(1)), x, n) 
       Nullify(fp_1d) 
    Case (2) 
       Call ESMF_FieldGetDataPointer(field, fp_2d, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'fp_2d') 
       Call CopyWindow_Integer_D(fp_2d(lb(1):ub(1),lb(2):ub(2)), x, n) 
       Nullify(fp_2d) 
    Case (3) 
       Call ESMF_FieldGetDataPointer(field, fp_3d, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'fp_3d') 
       Call CopyWindow_Integer_D(fp_3d(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)), x, n) 
       Nullify(fp_3d) 
    Case (4) 
       Call ESMF_FieldGetDataPointer(field, ptr=fp_4d, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'fp_4d') 
       Call CopyWindow_Integer_D(fp_4d(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4)), x, n) 
       Nullify(fp_4d) 
    Case (5) 
       Call ESMF_FieldGetDataPointer(field, fp_5d, rc=rc) 
       ASSERT(rc == ESMF_SUCCESS,'fp_5d') 
       Call CopyWindow_Integer_D(fp_5d(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5)), x, n) 
       Nullify(fp_5d) 
    End Select 
 
    Do rank = 1, field_rank 
 
       If (rank == dist_idx) Cycle 
 
       ldx = Product(shp_mask(:rank-1)) 
       m   = shp_mask(rank) 
       If (m == 1) cycle 
       n   = Product(shp_mask(rank+1:field_rank)) 
 
       Call LocalReduce_Integer_D(x(1), ldx, m, n, rtype) 
       shp_mask(rank) = 1 ! new shape 
 
    End Do 
 
 
 
    n_loc  = Max(0, shp_mask(dist_idx)) 
 
#ifdef USE_ESMF 
    comm = Field_Get_MPI_comm(field, npes, mype, root_pet) 
 
    Allocate(cnts(0:npes-1), dspls(0:npes-1), STAT=ier) 
    cnts=0 
 
    Call MPI_AllGather(n_loc, 1, MPI_INTEGER, cnts, 1, MPI_INTEGER, comm, ier) 
    Call MPI_AllReduce(lb(dist_idx), lb_glob, 1, MPI_Integer, MPI_MIN, comm, ier) 
    Call MPI_AllReduce(ub(dist_idx), ub_glob, 1, MPI_Integer, MPI_MAX, comm, ier) 
 
    dspls(0) = 0 
    Do pe = 1, npes-1 
       dspls(pe) = dspls(pe-1) + cnts(pe-1) 
    End Do 
 
    n_glob=Sum(cnts) 
    ASSERT(n_glob == ub_glob - lb_glob + 1,'inconsistent') 
    Allocate(psum(lb_glob:ub_glob), STAT=ier) 
 
    If (bcast) Then 
       Call MPI_AllGatherV(x, n_loc, MPI_INTEGER, psum, cnts, dspls, MPI_INTEGER, comm, ier) 
    Else 
       Call MPI_GatherV(x, n_loc, MPI_INTEGER, psum, cnts, dspls, MPI_INTEGER, root_pet, comm, ier) 
    End If 
#else 
    Allocate(psum(n_loc)) 
    psum = x(1:n_loc) 
#endif 
 
 
    ASSERT_ALWAYS(ANY(rtype == FIELD_REDUCE_OPS),'undefined reduction operation') 
 
    Select Case (rtype) 
    Case (FIELD_SUM) 
       result = Sum(psum) 
    Case (FIELD_MAX) 
       result = MaxVal(psum) 
    Case (FIELD_MIN) 
       result = MinVal(psum) 
    End Select 
 
 
    If (Present(hemi_result)) Then 
       !Call Grid_GetExtents(grid, equator=equator) 
 
       hemi_result=0 
       !!$$If (field_is_staggered) Then 
       !!$$   eq_sth = eq 
       !!$$   eq_nth = eq + 2 
       !!$$Else 
       !!$$   eq_sth = eq 
       !!$$   eq_nth = eq 
       !!$$End If 
       !!$$ 
       !!$$Select Case (rtype) 
       !!$$Case (FIELD_SUM) 
       !!$$   hemi_result(SH) = Sum(psum(:eq_sth)) 
       !!$$   hemi_result(NH) = Sum(psum(eq_nth:)) 
       !!$$   If (field_is_staggered) hemi_result = hemi_result + 0.5 * psum(eq+1) 
       !!$$Case (FIELD_MAX) 
       !!$$   hemi_result(SH) = MaxVal(psum(:eq_sth)) 
       !!$$   hemi_result(NH) = MaxVal(psum(eq_nth:)) 
       !!$$   If (field_is_staggered) hemi_result = Max(hemi_result, psum(eq+1)) 
       !!$$Case (FIELD_MIN) 
       !!$$   hemi_result(SH) = MinVal(psum(:eq_sth)) 
       !!$$   hemi_result(NH) = MinVal(psum(eq_nth:)) 
       !!$$   If (field_is_staggered) hemi_result = Min(hemi_result, psum(eq+1)) 
       !!$$End Select 
 
 
    End If 
 
    Deallocate(psum) 
#ifdef USE_ESMF 
    Deallocate(cnts, dspls) 
#endif 
    Deallocate(lb, ub, shp_mask) 
    Deallocate(x) 
 
  End Subroutine Field_Reduce_Integer_D 

  Subroutine Field_Gather(field, array, rc)
    USE ESMF_MOD ! only way to trick absoft into allowing "==" to be brought from ESMF_Base
    Use ESMF_MOD, Only: ESMF_Field
    Use ESMF_MOD, Only: ESMF_Array
    USE ESMF_Mod, Only: ESMF_KIND_R4
    USE ESMF_Mod, Only: ESMF_KIND_R8
    USE ESMF_Mod, Only: ESMF_DATA_REAL
    Use ESMF_MOD, Only: ESMF_DataType
    Use ESMF_MOD, Only: ESMF_DataKind
#ifdef USE_ESMF
    Use ESMF_MOD, Only: ESMF_FieldGet
    Use ESMF_MOD, Only: ESMF_FieldGetArray
    Use ESMF_MOD, Only: ESMF_ArrayCreate
    Use ESMF_MOD, Only: ESMF_ArrayGetData
    Use ESMF_MOD, Only: ESMF_ArrayF90Allocate
    Use ESMF_MOD, Only: ESMF_FieldDataMap
    Use ESMF_MOD, Only: ESMF_FieldDataMapGet
    Use ESMF_MOD, Only: ESMF_FieldGetDataPointer
#endif

    ! Arguments
    Type (ESMF_Field), Intent(InOut) :: field ! should be intent in (ESMF bug)
    Type (ESMF_Array), Intent(Out) :: array
    Integer, Optional, Intent(Out) :: rc

    Type (ESMF_DataType) :: type
    Type (ESMF_DataKind) :: kind

#ifdef USE_ESMF
    ! Local variables
    Type (ESMF_Array) :: f_array
    Integer :: mpi_type
    Integer :: base(1), buf_send(1), buf_recv(1)
    Integer :: loc_cnt
    Integer, Allocatable :: cnts(:), displs(:)
    Integer :: comm, root_pet = 0, npes, mpi_rank, root, pet
    Integer :: ier
    Integer :: rank, dist_idx
    Integer, Allocatable :: local_shape(:), global_shape(:) ! allocate to size "rank"
    Integer, Allocatable :: lb(:)
    Type (ESMF_FieldDataMap) :: map
    
    Real (Kind=ESMF_KIND_R4), Pointer :: f_ptr_Real4_3D (:,:,:) 
    Real (Kind=ESMF_KIND_R4), Pointer :: a_ptr_Real4_3D (:,:,:) 

    
    Real (Kind=ESMF_KIND_R8), Pointer :: f_ptr_Real8_1D (:) 
    Real (Kind=ESMF_KIND_R8), Pointer :: a_ptr_Real8_1D (:) 

    
    Real (Kind=ESMF_KIND_R8), Pointer :: f_ptr_Real8_2D (:,:) 
    Real (Kind=ESMF_KIND_R8), Pointer :: a_ptr_Real8_2D (:,:) 

    
    Real (Kind=ESMF_KIND_R8), Pointer :: f_ptr_Real8_3D (:,:,:) 
    Real (Kind=ESMF_KIND_R8), Pointer :: a_ptr_Real8_3D (:,:,:) 

    
    Real (Kind=ESMF_KIND_R8), Pointer :: f_ptr_Real8_4D (:,:,:,:) 
    Real (Kind=ESMF_KIND_R8), Pointer :: a_ptr_Real8_4D (:,:,:,:) 

    
    Real (Kind=ESMF_KIND_R8), Pointer :: f_ptr_Real8_5D (:,:,:,:,:) 
    Real (Kind=ESMF_KIND_R8), Pointer :: a_ptr_Real8_5D (:,:,:,:,:) 

    
    Integer , Pointer :: f_ptr_Integer_1D (:) 
    Integer , Pointer :: a_ptr_Integer_1D (:) 

    
    Integer , Pointer :: f_ptr_Integer_2D (:,:) 
    Integer , Pointer :: a_ptr_Integer_2D (:,:) 

    
    Integer , Pointer :: f_ptr_Integer_3D (:,:,:) 
    Integer , Pointer :: a_ptr_Integer_3D (:,:,:) 

    comm = Field_Get_MPI_comm(field, comm_size=npes, comm_rank=mpi_rank, comm_root=root, rc=rc)

    ! Determine rank of field array.
    !-------------------------------
    Call ESMF_FieldGet(field, array=f_array, rc=rc)
    Call ESMF_ArrayGet(f_array, rank=rank, type=type, kind=kind, rc=rc)
    Allocate(local_shape(rank), global_shape(rank), lb(rank), STAT = ier);

    ! create a new mpi type for use in communication
    !-------------------------------
    Call Field_AnalyzeCommPattern(field, mpi_type, dist_idx=dist_idx, counts=local_shape, rc=rc)

    ! account for halo width
    !-------------------------------
    loc_cnt = local_shape(dist_idx) - 2 * haloWidth

    ! Gather the local extents in distributed direction.
    ! Probably could get this by some non-communication means for optimization.
    !-------------------------------
    Allocate(cnts(0:npes-1), displs(0:npes-1), STAT=ier)
    Call MPI_Gather(loc_cnt, 1, MPI_INTEGER, cnts, 1, MPI_INTEGER, root_pet, comm, ier)

    ! Determine global shape for result array, as well as displs based on cnts.
    !-------------------------------
    If (mpi_rank == root_pet) Then
       displs(0) = 0
       Do pet = 1, npes - 1
          displs(pet) = displs(pet-1) + cnts(pet-1)
       End Do
       global_shape = local_shape ! Except at dist_idx, it is the total domain.
       global_shape(dist_idx) = Sum(cnts)
    Else
       global_shape = 1 ! non empty for MPI
    End If

    ! Create an array from global shape.
    ! To do - worry about arrays that do not start at "1".
    !-------------------------------
    array = ESMF_ArrayCreate( rank=rank, type=type, kind=kind, counts=global_shape, rc=rc)

    ! Define a template for gathering various data types, kinds, and ranks.
    ! A generic pointer would make this only need a single instantiation,
    ! but this mechanism should be more portable.
    !-------------------------------
    
     
    If (rank == 3 .and. type == ESMF_DATA_REAL .and. kind == ESMF_R4) Then 
       Call ESMF_FieldGetDataPointer(field, ptr=f_ptr_Real4_3D, rc=rc) 
       Call ESMF_ArrayGetData( array, fptr=a_ptr_Real4_3D, rc=rc) 
        
       lb = LBound( f_ptr_Real4_3D) 
       lb(dist_idx) = lb(dist_idx) + haloWidth 
        
       Call MPI_GatherV(f_ptr_Real4_3D (lb(1),lb(2),lb(3)), loc_cnt, mpi_type, & 
            &           a_ptr_Real4_3D,       cnts,  displs, mpi_type, & 
            & root_pet, comm, ier) 
       Nullify(f_ptr_Real4_3D) 
    End If 
     

    
     
    If (rank == 1 .and. type == ESMF_DATA_REAL .and. kind == ESMF_R8) Then 
       Call ESMF_FieldGetDataPointer(field, ptr=f_ptr_Real8_1D, rc=rc) 
       Call ESMF_ArrayGetData( array, fptr=a_ptr_Real8_1D, rc=rc) 
        
       lb = LBound( f_ptr_Real8_1D) 
       lb(dist_idx) = lb(dist_idx) + haloWidth 
        
       Call MPI_GatherV(f_ptr_Real8_1D (lb(1)), loc_cnt, mpi_type, & 
            &           a_ptr_Real8_1D,       cnts,  displs, mpi_type, & 
            & root_pet, comm, ier) 
       Nullify(f_ptr_Real8_1D) 
    End If 
     

    
     
    If (rank == 2 .and. type == ESMF_DATA_REAL .and. kind == ESMF_R8) Then 
       Call ESMF_FieldGetDataPointer(field, ptr=f_ptr_Real8_2D, rc=rc) 
       Call ESMF_ArrayGetData( array, fptr=a_ptr_Real8_2D, rc=rc) 
        
       lb = LBound( f_ptr_Real8_2D) 
       lb(dist_idx) = lb(dist_idx) + haloWidth 
        
       Call MPI_GatherV(f_ptr_Real8_2D (lb(1),lb(2)), loc_cnt, mpi_type, & 
            &           a_ptr_Real8_2D,       cnts,  displs, mpi_type, & 
            & root_pet, comm, ier) 
       Nullify(f_ptr_Real8_2D) 
    End If 
     

    
     
    If (rank == 3 .and. type == ESMF_DATA_REAL .and. kind == ESMF_R8) Then 
       Call ESMF_FieldGetDataPointer(field, ptr=f_ptr_Real8_3D, rc=rc) 
       Call ESMF_ArrayGetData( array, fptr=a_ptr_Real8_3D, rc=rc) 
        
       lb = LBound( f_ptr_Real8_3D) 
       lb(dist_idx) = lb(dist_idx) + haloWidth 
        
       Call MPI_GatherV(f_ptr_Real8_3D (lb(1),lb(2),lb(3)), loc_cnt, mpi_type, & 
            &           a_ptr_Real8_3D,       cnts,  displs, mpi_type, & 
            & root_pet, comm, ier) 
       Nullify(f_ptr_Real8_3D) 
    End If 
     

    
     
    If (rank == 4 .and. type == ESMF_DATA_REAL .and. kind == ESMF_R8) Then 
       Call ESMF_FieldGetDataPointer(field, ptr=f_ptr_Real8_4D, rc=rc) 
       Call ESMF_ArrayGetData( array, fptr=a_ptr_Real8_4D, rc=rc) 
        
       lb = LBound( f_ptr_Real8_4D) 
       lb(dist_idx) = lb(dist_idx) + haloWidth 
        
       Call MPI_GatherV(f_ptr_Real8_4D (lb(1),lb(2),lb(3),lb(4)), loc_cnt, mpi_type, & 
            &           a_ptr_Real8_4D,       cnts,  displs, mpi_type, & 
            & root_pet, comm, ier) 
       Nullify(f_ptr_Real8_4D) 
    End If 
     

    
     
    If (rank == 5 .and. type == ESMF_DATA_REAL .and. kind == ESMF_R8) Then 
       Call ESMF_FieldGetDataPointer(field, ptr=f_ptr_Real8_5D, rc=rc) 
       Call ESMF_ArrayGetData( array, fptr=a_ptr_Real8_5D, rc=rc) 
        
       lb = LBound( f_ptr_Real8_5D) 
       lb(dist_idx) = lb(dist_idx) + haloWidth 
        
       Call MPI_GatherV(f_ptr_Real8_5D (lb(1),lb(2),lb(3),lb(4),lb(5)), loc_cnt, mpi_type, & 
            &           a_ptr_Real8_5D,       cnts,  displs, mpi_type, & 
            & root_pet, comm, ier) 
       Nullify(f_ptr_Real8_5D) 
    End If 
     

    
     
    If (rank == 1 .and. type == ESMF_DATA_INTEGER .and. kind == ESMF_I4) Then 
       Call ESMF_FieldGetDataPointer(field, ptr=f_ptr_Integer_1D, rc=rc) 
       Call ESMF_ArrayGetData( array, fptr=a_ptr_Integer_1D, rc=rc) 
        
       lb = LBound( f_ptr_Integer_1D) 
       lb(dist_idx) = lb(dist_idx) + haloWidth 
        
       Call MPI_GatherV(f_ptr_Integer_1D (lb(1)), loc_cnt, mpi_type, & 
            &           a_ptr_Integer_1D,       cnts,  displs, mpi_type, & 
            & root_pet, comm, ier) 
       Nullify(f_ptr_Integer_1D) 
    End If 
     

    
     
    If (rank == 2 .and. type == ESMF_DATA_INTEGER .and. kind == ESMF_I4) Then 
       Call ESMF_FieldGetDataPointer(field, ptr=f_ptr_Integer_2D, rc=rc) 
       Call ESMF_ArrayGetData( array, fptr=a_ptr_Integer_2D, rc=rc) 
        
       lb = LBound( f_ptr_Integer_2D) 
       lb(dist_idx) = lb(dist_idx) + haloWidth 
        
       Call MPI_GatherV(f_ptr_Integer_2D (lb(1),lb(2)), loc_cnt, mpi_type, & 
            &           a_ptr_Integer_2D,       cnts,  displs, mpi_type, & 
            & root_pet, comm, ier) 
       Nullify(f_ptr_Integer_2D) 
    End If 
     

    
     
    If (rank == 3 .and. type == ESMF_DATA_INTEGER .and. kind == ESMF_I4) Then 
       Call ESMF_FieldGetDataPointer(field, ptr=f_ptr_Integer_3D, rc=rc) 
       Call ESMF_ArrayGetData( array, fptr=a_ptr_Integer_3D, rc=rc) 
        
       lb = LBound( f_ptr_Integer_3D) 
       lb(dist_idx) = lb(dist_idx) + haloWidth 
        
       Call MPI_GatherV(f_ptr_Integer_3D (lb(1),lb(2),lb(3)), loc_cnt, mpi_type, & 
            &           a_ptr_Integer_3D,       cnts,  displs, mpi_type, & 
            & root_pet, comm, ier) 
       Nullify(f_ptr_Integer_3D) 
    End If 
     

    ! Clean up
    !-------------------------------
    Deallocate(cnts, displs)
    Deallocate(local_shape, global_shape, lb)
    Call MPI_Type_Free(mpi_type, ier)

#else

    ! Gather just is a copy of the existing array.
    Call ESMF_FieldGetArray(field, anArray=array, rc=rc)
#endif

  End Subroutine Field_Gather

  !< [[Field_Scatter]] >>

  ! Internal methods
#ifdef USE_ESMF
  Subroutine Field_AnalyzeCommPattern(field, mpi_field_type, dist_idx, counts, rc)
    Use ESMF_MOD, Only: ESMF_SUCCESS
    Use ESMF_MOD, Only: ESMF_Field
    Use ESMF_MOD, Only: ESMF_FieldGet
    Use ESMF_MOD, Only: ESMF_FieldGetArray
    Use ESMF_MOD, Only: ESMF_FieldDataMap
    Use ESMF_MOD, Only: ESMF_Array
    Use ESMF_MOD, Only: ESMF_ArrayGet
    Use ESMF_MOD, Only: ESMF_DataType
    Use ESMF_MOD, Only: ESMF_DataKind
  ! Arguments
    Type (ESMF_Field), Intent(In) :: field
    Integer, Intent(Out) :: mpi_field_type
    Integer, Intent(Out) :: dist_idx
    Integer, Intent(Out) :: counts(:)
    Integer, Optional, Intent(Out) :: rc

  ! Local variables

    Integer :: mpi_base_type
    Type (ESMF_DataType) :: esmf_type
    Type (ESMF_DataKind) :: esmf_kind
    Type (ESMF_Array) :: array_temp
    Type (ESMF_FieldDataMap) :: map
    Integer :: idx
    Integer :: rank
    Integer :: ier

    Call ESMF_FieldGetArray(field, array=array_temp, rc=rc)
    Call ESMF_ArrayGet(array_temp, rank=rank, type=esmf_type, kind=esmf_kind, counts=counts, rc=rc)

    dist_idx = Field_FindIndex(field, J_GRID, rc = rc)
    mpi_base_type = Get_MPI_Type(esmf_type, esmf_kind)
    mpi_field_type = CreateDist_MPI_Type(mpi_base_type, counts, dist_idx)

    If (present(rc)) rc = ESMF_SUCCESS
    Return

  End Subroutine Field_AnalyzeCommPattern

  Function CreateDist_MPI_Type(base_type, counts, dist_idx) Result(new_type)
    Integer, Intent(In) :: base_type
    Integer, Intent(In) :: counts(:)
    Integer, Intent(In) :: dist_idx
    Integer :: new_type

    Integer :: stride, n_blocks, blocklen
    Integer :: ext_lb
    Integer :: base_byte_len, new_len
    Integer :: vector_type
    Integer :: ier

    n_blocks = Product(counts(dist_idx+1:))
    blocklen = Product(counts(:dist_idx-1))
    stride = counts(dist_idx) * blocklen

    Call MPI_Type_vector(n_blocks, blocklen, stride, base_type, vector_type, ier)
    Call MPI_Type_extent(base_type, base_byte_len, ier)
    new_len = base_byte_len * blocklen
    Call MPI_Type_struct(2, (/ 1, 1 /), (/ 0, new_len /), &
         & (/ vector_type, MPI_UB /), vector_type, ier)
    Call MPI_Type_Commit(vector_type, ier)

    new_type = vector_type

  End Function CreateDist_MPI_Type


  Function Field_Get_MPI_comm(field, comm_size, comm_rank, comm_root, rc) Result(comm)
    Use ESMF_MOD, Only: ESMF_Field
    Use ESMF_MOD, Only: ESMF_FieldGet
    Use ESMF_MOD, Only: ESMF_VM
    Use ESMF_MOD, Only: ESMF_VMGetGlobal
    Use ESMF_MOD, Only: ESMF_VMGet

    ! Arguments
    Type (ESMF_Field), Intent(In) :: field ! ignored for now
    Integer, Intent(Out) :: comm_size
    Integer, Intent(Out) :: comm_rank
    Integer, Intent(Out) :: comm_root
    Integer, Intent(Out), Optional :: rc
    Integer :: comm ! the mpi communicator
  !--------------------------------------------

    Type (ESMF_VM) :: vm

    Call ESMF_VMGetGlobal(vm, rc)
    Call ESMF_VMGet(vm, petCount=comm_size, localPet=comm_rank, mpiCommunicator=comm, rc=rc)
    comm_root = 0

  End Function Field_Get_MPI_comm

  Function Get_MPI_Type(esmf_type, esmf_kind) Result(mpi_type)
    Use ESMF_MOD ! Cannot get ".eq" for type and kind with "Only" on Absoft
    Use ESMF_MOD, Only: ESMF_R4
    Use ESMF_MOD, Only: ESMF_DataType
    Use ESMF_MOD, Only: ESMF_DataKind
    Use ESMF_MOD, Only: ESMF_DATA_REAL
  ! Arguments
    Type (ESMF_DataType) :: esmf_type
    Type (ESMF_DataKind) :: esmf_kind
    Integer :: mpi_type
  !----------------------------------------

    If (esmf_type == ESMF_DATA_REAL) Then
       If (esmf_kind == ESMF_R4) Then
          mpi_type = MPI_REAL
       Else
          mpi_type = MPI_DOUBLE_PRECISION
       End If
    Else
       mpi_type = MPI_INTEGER
    End If

  End Function Get_MPI_Type

  Subroutine GetNeighbors(rank, npes, pe_south, pe_north)
    Integer, Intent(In) :: rank
    Integer, Intent(In) :: npes
    Integer, Intent(Out) :: pe_south
    Integer, Intent(Out) :: pe_north

    If (rank > 0) Then
       pe_south = rank - 1
    Else
       pe_south = MPI_PROC_NULL
    End If

    If (rank < npes-1) Then
       pe_north = rank + 1
    Else
       pe_north = MPI_PROC_NULL
    End If

  End Subroutine GetNeighbors

#endif
  Function Field_GetRank(field, rc) Result(rank)
    Use ESMF_MOD, Only: ESMF_Field
    Use ESMF_MOD, Only: ESMF_SUCCESS
#ifdef USE_ESMF
    Use ESMF_MOD, Only: ESMF_FieldGet
    Use ESMF_MOD, Only: ESMF_FieldDataMap
    Use ESMF_MOD, Only: ESMF_FieldDataMapGet
#else
    USE ESMF_MOD_private, Only: GetRank_priv => Field_GetRank
#endif
    Type (ESMF_Field), Intent(In) :: field
    Integer, Optional, Intent(Out) :: rc
    Integer :: rank ! result


#ifdef USE_ESMF
    Type (ESMF_FieldDataMap) :: map

    Call ESMF_FieldGet(field, datamap=map, rc=rc)
    Call ESMF_FieldDataMapGet(map, datarank=rank, rc=rc)

#else
    rank = GetRank_priv(field)
#endif

    If (Present(rc)) rc = ESMF_SUCCESS
  End Function Field_GetRank

  Subroutine Field_GetBounds(field, lb, ub, rc)
    Use ESMF_MOD, Only: ESMF_Field
    USE ESMF_MOD, Only: ESMF_SUCCESS
#ifdef USE_ESMF
    Use ESMF_MOD, Only: ESMF_FieldGet
    Use ESMF_MOD, Only: ESMF_FieldDataMap
    Use ESMF_MOD, Only: ESMF_FieldDataMapGet
    Use ESMF_MOD, Only: ESMF_Array
    Use ESMF_MOD, Only: ESMF_ArrayGet
#else
    Use ESMF_MOD_private, Only: Field_GetShape
#endif

    Type (ESMF_Field), Intent(In) :: field
    Integer, Intent(Out) :: lb(:)
    Integer, Intent(Out) :: ub(:)
    Integer, Intent(Out) :: rc
    Integer :: rank ! result
#ifdef USE_ESMF
    Type (ESMF_Array) :: array
#endif
    Integer :: dist_idx

    rank = Field_GetRank(field,rc=rc)

#ifdef USE_ESMF
    Call ESMF_FieldGet(field, array=array, rc=rc)
    Call ESMF_ArrayGet(array, lbounds=lb, ubounds=ub, rc=rc)

    dist_idx = Field_FindIndex(field, J_GRID, rc=rc)

    lb(dist_idx) = lb(dist_idx) + halowidth
    ub(dist_idx) = ub(dist_idx) - halowidth
#else
    lb = 1
    Call Field_GetShape(field, ub)
    rc = ESMF_SUCCESS
#endif

  End Subroutine Field_GetBounds



  Integer Function Field_FindIndex(field, idx, rc)
    Use ESMF_MOD, Only: ESMF_Field
    Use ESMF_MOD, Only: ESMF_SUCCESS
#ifdef USE_ESMF
    Use ESMF_MOD, Only: ESMF_FieldGet
    Use ESMF_MOD, Only: ESMF_FieldDataMap
    Use ESMF_MOD, Only: ESMF_FieldDataMapGet
    Use ESMF_MOD, Only: ESMF_FAILURE
#else
    Use ESMF_MOD_private, Only: Field_GetDistRank
#endif
    Type (ESMF_Field), Intent(In) :: field
    Integer, Intent(In) :: idx
    Integer, Optional, Intent(Out) :: rc

#ifdef USE_ESMF
    Integer :: rank
    Integer :: ier
    Type (ESMF_FieldDataMap) :: map
    Integer, Allocatable :: dataIndices(:)

    Call ESMF_FieldGet(field, datamap=map, rc=rc)

    Call ESMF_FieldDataMapGet(map, datarank=rank, rc=rc)
    Allocate(dataindices(rank), STAT=ier)
    ASSERT_ALWAYS(ier == 0, 'allocate()')
    Call ESMF_FieldDataMapGet(map, dataIndexList=dataIndices, rc=rc)

    Field_FindIndex = FindIndex(dataIndices, idx, rc=rc)
    Deallocate(dataindices)
#else
    Field_FindIndex = Field_GetDistRank(field)
#endif
    If (Present(rc)) rc = ESMF_SUCCESS
  End Function Field_FindIndex

  Integer Function FindIndex(list, idx, rc)
    Use ESMF_MOD, Only: ESMF_SUCCESS
    Use ESMF_MOD, Only: ESMF_FAILURE
    Integer, Intent(In) :: list(:)
    Integer, Intent(In) :: idx
    Integer, Optional, Intent(Out) :: rc

    Integer :: i
    Logical :: found

    found =.false.
    Do i = 1, Size(list)
       If (list(i) == idx) Then
          found = .true.
          Exit
       End If
    End Do

    If (found) Then
       FindIndex = i
       If (Present(rc)) rc = ESMF_SUCCESS
    Else
       FindIndex = -1
       If (Present(rc)) rc = ESMF_FAILURE
    End If

  End Function FindIndex
  
  Subroutine CopyWindow_Real4_D(x, y, n) 
    USE ESMF_MOD, Only: ESMF_KIND_R4 
    USE ESMF_MOD, Only: ESMF_KIND_R8 
    Integer :: n 
    Real (Kind=ESMF_KIND_R4), Intent(In)  :: x(n) 
    Real (Kind=ESMF_KIND_R4), Intent(Out) :: y(n) 
 
    y = x 
 
  End Subroutine CopyWindow_Real4_D 

  
  Subroutine CopyWindow_Real8_D(x, y, n) 
    USE ESMF_MOD, Only: ESMF_KIND_R4 
    USE ESMF_MOD, Only: ESMF_KIND_R8 
    Integer :: n 
    Real (Kind=ESMF_KIND_R8), Intent(In)  :: x(n) 
    Real (Kind=ESMF_KIND_R8), Intent(Out) :: y(n) 
 
    y = x 
 
  End Subroutine CopyWindow_Real8_D 

  
  Subroutine CopyWindow_Integer_D(x, y, n) 
    USE ESMF_MOD, Only: ESMF_KIND_R4 
    USE ESMF_MOD, Only: ESMF_KIND_R8 
    Integer :: n 
    Integer , Intent(In)  :: x(n) 
    Integer , Intent(Out) :: y(n) 
 
    y = x 
 
  End Subroutine CopyWindow_Integer_D 

  
  Subroutine LocalReduce_Real4_D(x, ldx, n, m, rtype) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Integer, Intent(In) :: ldx, n, m 
    Real (Kind=ESMF_KIND_R4) :: x(ldx * n * m) 
    Integer, Intent(In) :: rtype 
 
    Integer :: i, j, k 
    Integer :: n1, n2, n3 
    Real (Kind=ESMF_KIND_R4) :: s 
 
    Select Case (rtype) 
    Case (FIELD_SUM) 
       n2 = 1 
       n3 = 1 
       Do k = 1, m 
          Do i = 1, ldx 
             n1 = n3 
             s = 0 
             Do j = 1, n 
                s = s + x(n1) 
                n1 = n1 + ldx 
             End Do 
             x(n2) = s 
             n2 = n2 + 1 
             n3 = n3 + 1 
          End Do 
          n3 = ldx*n*k + 1 
       End Do 
    Case (FIELD_MAX) 
       n2 = 1 
       n3 = 1 
       Do k = 1, m 
          Do i = 1, ldx 
             n1 = n3 
             s = x(n1) 
             Do j = 1, n 
                s = Max(s, x(n1)) 
                n1 = n1 + ldx 
             End Do 
             x(n2) = s 
             n2 = n2 + 1 
             n3 = n3 + 1 
          End Do 
          n3 = ldx*n*k + 1 
       End Do 
    Case (FIELD_MIN) 
       n2 = 1 
       n3 = 1 
       Do k = 1, m 
          Do i = 1, ldx 
             n1 = n3 
             s = x(n1) 
             Do j = 1, n 
                s = Min(s, x(n1)) 
                n1 = n1 + ldx 
             End Do 
             x(n2) = s 
             n2 = n2 + 1 
             n3 = n3 + 1 
          End Do 
          n3 = ldx*n*k + 1 
       End Do 
    End Select 
        
  End Subroutine LocalReduce_Real4_D 

  
  Subroutine LocalReduce_Real8_D(x, ldx, n, m, rtype) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Integer, Intent(In) :: ldx, n, m 
    Real (Kind=ESMF_KIND_R8) :: x(ldx * n * m) 
    Integer, Intent(In) :: rtype 
 
    Integer :: i, j, k 
    Integer :: n1, n2, n3 
    Real (Kind=ESMF_KIND_R8) :: s 
 
    Select Case (rtype) 
    Case (FIELD_SUM) 
       n2 = 1 
       n3 = 1 
       Do k = 1, m 
          Do i = 1, ldx 
             n1 = n3 
             s = 0 
             Do j = 1, n 
                s = s + x(n1) 
                n1 = n1 + ldx 
             End Do 
             x(n2) = s 
             n2 = n2 + 1 
             n3 = n3 + 1 
          End Do 
          n3 = ldx*n*k + 1 
       End Do 
    Case (FIELD_MAX) 
       n2 = 1 
       n3 = 1 
       Do k = 1, m 
          Do i = 1, ldx 
             n1 = n3 
             s = x(n1) 
             Do j = 1, n 
                s = Max(s, x(n1)) 
                n1 = n1 + ldx 
             End Do 
             x(n2) = s 
             n2 = n2 + 1 
             n3 = n3 + 1 
          End Do 
          n3 = ldx*n*k + 1 
       End Do 
    Case (FIELD_MIN) 
       n2 = 1 
       n3 = 1 
       Do k = 1, m 
          Do i = 1, ldx 
             n1 = n3 
             s = x(n1) 
             Do j = 1, n 
                s = Min(s, x(n1)) 
                n1 = n1 + ldx 
             End Do 
             x(n2) = s 
             n2 = n2 + 1 
             n3 = n3 + 1 
          End Do 
          n3 = ldx*n*k + 1 
       End Do 
    End Select 
        
  End Subroutine LocalReduce_Real8_D 

  
  Subroutine LocalReduce_Integer_D(x, ldx, n, m, rtype) 
    Use ESMF_MOD, Only: ESMF_KIND_R4 
    Use ESMF_MOD, Only: ESMF_KIND_R8 
    Integer, Intent(In) :: ldx, n, m 
    Integer  :: x(ldx * n * m) 
    Integer, Intent(In) :: rtype 
 
    Integer :: i, j, k 
    Integer :: n1, n2, n3 
    Integer  :: s 
 
    Select Case (rtype) 
    Case (FIELD_SUM) 
       n2 = 1 
       n3 = 1 
       Do k = 1, m 
          Do i = 1, ldx 
             n1 = n3 
             s = 0 
             Do j = 1, n 
                s = s + x(n1) 
                n1 = n1 + ldx 
             End Do 
             x(n2) = s 
             n2 = n2 + 1 
             n3 = n3 + 1 
          End Do 
          n3 = ldx*n*k + 1 
       End Do 
    Case (FIELD_MAX) 
       n2 = 1 
       n3 = 1 
       Do k = 1, m 
          Do i = 1, ldx 
             n1 = n3 
             s = x(n1) 
             Do j = 1, n 
                s = Max(s, x(n1)) 
                n1 = n1 + ldx 
             End Do 
             x(n2) = s 
             n2 = n2 + 1 
             n3 = n3 + 1 
          End Do 
          n3 = ldx*n*k + 1 
       End Do 
    Case (FIELD_MIN) 
       n2 = 1 
       n3 = 1 
       Do k = 1, m 
          Do i = 1, ldx 
             n1 = n3 
             s = x(n1) 
             Do j = 1, n 
                s = Min(s, x(n1)) 
                n1 = n1 + ldx 
             End Do 
             x(n2) = s 
             n2 = n2 + 1 
             n3 = n3 + 1 
          End Do 
          n3 = ldx*n*k + 1 
       End Do 
    End Select 
        
  End Subroutine LocalReduce_Integer_D 

End Module FieldComm
