      MODULE TRIDIAG_MOD
!@sum TRIDIAG_MOD contains subroutine TRIDIAG
      PRIVATE
      PUBLIC TRIDIAG
      PUBLIC TRIDIAG_NEW

      Interface Tridiag_new
      Module Procedure tridiag
      Module procedure tridiag_2d_glob
      Module procedure tridiag_2d_dist
      End Interface

      contains

      SUBROUTINE TRIDIAG(A,B,C,R,U,N)
!@sum  TRIDIAG  solves a tridiagonal matrix equation (A,B,C)U=R
!@auth Numerical Recipes
!@ver  1.0
      IMPLICIT NONE
      INTEGER, PARAMETER :: NMAX = 5000  !@var NMAX workspace
      INTEGER, INTENT(IN):: N         !@var N    dimension of arrays
      REAL*8, INTENT(IN) :: A(N)   !@var A    coefficients of u_i-1
      REAL*8, INTENT(IN) :: B(N)   !@var B    coefficients of u_i
      REAL*8, INTENT(IN) :: C(N)   !@var C    coefficients of u_i+1
      REAL*8, INTENT(IN) :: R(N)   !@var R    RHS vector
      REAL*8, INTENT(OUT):: U(N)   !@var U    solution vector
      REAL*8 :: BET                   !@var BET  work variable
      REAL*8 :: GAM(NMAX)             !@var GAM  work array
      INTEGER :: J                    !@var J    loop variable

      IF ( N > NMAX )
     &     call stop_model("TRIDIAG: N > NMAX, increase NMAX",255)
      BET=B(1)
      IF (BET.eq.0) call stop_model("TRIDIAG: DENOMINATOR = ZERO",255)
      U(1)=R(1)/BET
      DO J=2,N
        GAM(J)=C(J-1)/BET
        BET=B(J)-A(J)*GAM(J)
        IF (BET.eq.0) call stop_model("TRIDIAG: DENOMINATOR = ZERO",255)
        U(J)=(R(J)-A(J)*U(J-1))/BET
      END DO
      DO J=N-1,1,-1
        U(J)=U(J)-GAM(J+1)*U(J+1)
      END DO
      RETURN
      END SUBROUTINE TRIDIAG

      SUBROUTINE TRIDIAG_2D_GLOB(A, B, C, R, U)
!@sum  TRIDIAG  solves an array of tridiagonal matrix equations (A,B,C)U=R
!@auth Numerical Recipes
!@ver  1.0
      IMPLICIT NONE

      REAL*8, INTENT(INOUT) :: A(:,:),B(:,:),C(:,:),R(:,:)
      REAL*8, INTENT(OUT)   :: U(:,:)

      REAL*8 :: BET                   !@var BET  work variable
      REAL*8 :: GAM(Size(A,2))        !@var GAM  work array

      Integer :: i, j
      Integer :: N, N_arr

      N = SIZE(A,2)
      n_arr = Size(A,1)
      do i=1,n_arr
        BET=B(i,1)
        IF (BET.eq.0) then
         print*, "TRIDIAG_2D_GLOB: DENOMINATOR = ZERO"
         stop
        end if
        U(i,1)=R(i,1)/BET
        DO J=2,N
          GAM(J)=C(i,J-1)/BET
          BET=B(i,J)-A(i,J)*GAM(J)
          IF (BET.eq.0) then
           print*, "TRIDIAG: DENOMINATOR = ZERO"
           stop
          end if
          U(i,J)=( R(i,J) - A(i,J)*U(i,J-1) )/BET
        END DO
        DO J=N-1,1,-1
          U(i,J)=U(i,J)-GAM(J+1)*U(i,J+1)
        END DO
      end do

      RETURN
      END SUBROUTINE TRIDIAG_2D_GLOB


      SUBROUTINE TRIDIAG_2D_DIST(A_dist, B_dist, C_dist, R_dist,
     &                           U_dist,grid, j_lower, j_upper )
!@sum  TRIDIAG  solves an array of tridiagonal matrix equations (A,B,C)U=R
!@auth Numerical Recipes
!@ver  1.0
      USE DOMAIN_DECOMP, ONLY : DIST_GRID
      USE DOMAIN_DECOMP, ONLY : TRANSP
      IMPLICIT NONE

      Type (DIST_GRID), Intent(IN) :: grid
      REAL*8, INTENT(INOUT) :: A_dist(:,grid%j_strt_halo:)
      REAL*8, INTENT(INOUT) :: B_dist(:,grid%j_strt_halo:)
      REAL*8, INTENT(INOUT) :: C_dist(:,grid%j_strt_halo:)
      REAL*8, INTENT(INOUT) :: R_dist(:,grid%j_strt_halo:)
      REAL*8, INTENT(OUT)   :: U_dist(:,grid%j_strt_halo:)
      INTEGER, INTENT(IN)   :: J_LOWER, J_UPPER

      REAL*8, ALLOCATABLE :: A_tr(:,:), B_tr(:,:),
     &                       C_tr(:,:), R_tr(:,:),
     &                       U_tr(:,:)

      REAL*8 :: BET                   !@var BET  work variable
      REAL*8 :: GAM(grid%jm_world)        !@var GAM  work array

      Integer :: i, j
      Integer :: N, N_i


! Determine the size of the global arrays
      N = grid%jm_world
      n_i = grid%ni_loc

! Matrix size consistent with array size?
      if ( J_upper > N ) then
        print*, 'TRIDIAG: upper bound of matrix arrays is too large'
        print*, 'j_upper = ', j_upper, 'jm =', n, ' ( need j_upper<=jm)'
        call stop_model('TRIDIAG: j_upper argument too large', 255)
      end if

! Allocate the transposed arrays
      allocate( a_tr(n_i,n), b_tr(n_i,n), c_tr(n_i,n), r_tr(n_i,n) )
      allocate( u_tr(n_i,n) )

! First get the transpose of A,B,C,R
      call transp(grid, a_dist, a_tr)
      call transp(grid, b_dist, b_tr)
      call transp(grid, c_dist, c_tr)
      call transp(grid, r_dist, r_tr)

! Solve
        do i=1,n_i
          BET=B_tr(i,j_lower)
          IF (BET.eq.0) then
           print*, "TRIDIAG_2D_DIST: DENOM. = ZERO  i,j= ", i,' 1'
           stop
          end if
          U_tr(i,j_lower)=R_tr(i,j_lower)/BET
          DO J=j_lower+1, j_upper
            GAM(J)=C_tr(i,J-1)/BET
            BET=B_tr(i,J)-A_tr(i,J)*GAM(J)
            IF (BET.eq.0) then
             print*, "TRIDIAG_2D_DIST: DENOM. = ZERO i,j= ", i, j
             stop
            end if
          U_tr(i,J)=( R_tr(i,J) - A_tr(i,J)*U_tr(i,J-1) )/BET
          END DO
          DO J=j_upper-1,j_lower,-1
            U_tr(i,J)=U_tr(i,J)-GAM(J+1)*U_tr(i,J+1)
          END DO
        end do


! Transfer the solution to the j-distributed array
        call transp( grid, u_dist, u_tr, reverse=.true.)

      RETURN
      END SUBROUTINE TRIDIAG_2D_DIST



      END MODULE TRIDIAG_MOD
