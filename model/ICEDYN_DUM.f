      MODULE ICEDYN
      integer :: imic=0
      contains
      subroutine alloc_icedyn(grid)
      use DOMAIN_DECOMP, only : DYN_GRID
      TYPE (DYN_GRID), INTENT(IN) :: grid
      return
      end subroutine alloc_icedyn
      END MODULE ICEDYN

      MODULE ICEDYN_COM
      integer :: kticij=0
      real*8,dimension(2) :: rsix=0,rsiy=0,usi=0,vsi=0
      real*8,dimension(2) :: USIDT=0,VSIDT=0,RSISAVE=0
      real*8,dimension(2) :: icij=0,ticij=0
      contains
      subroutine alloc_icedyn_com(grid)
      use DOMAIN_DECOMP, only : DYN_GRID
      TYPE (DYN_GRID), INTENT(IN) :: grid
      return
      end subroutine alloc_icedyn_com
      END MODULE ICEDYN_COM

      SUBROUTINE ICEDYN_DUM
!@sum ICEDYN_DUM dummy routines to replace ice dynamics
      ENTRY alloc_icedyn
      ENTRY alloc_icedyn_com
      ENTRY io_icedyn
      ENTRY io_icdiag
      ENTRY reset_icdiag
      ENTRY ADVSI
      ENTRY init_icedyn
      ENTRY diag_ICEDYN
      RETURN
      END SUBROUTINE ICEDYN_DUM

      SUBROUTINE DYNSI
!@sum DYNSI simple coding to estimate ice-ocean friction velocity
!@auth Gavin Schmidt
      USE CONSTANT, only : rhows
      USE MODEL_COM, only : im,jm,kocean,focean,dtsrc
      USE GEOM, only : imaxj
      USE SEAICE, only : oi_ustar0
      USE SEAICE_COM, only : rsi
      USE FLUXES, only : UI2rho,dmua,dmva

      IMPLICIT NONE
      INTEGER I,J
      REAL*8 ustar1

      IF (KOCEAN.ge.1) THEN
        DO J=1,JM
        DO I=1,IMAXJ(J)
c          UI2rho(I,J) = rhows*(oi_ustar0)**2  ! default
C**** with wind stress dependence
          if (rsi(i,j)*focean(i,j).gt.0) then
            ustar1= SQRT(SQRT(DMUA(I,J,2)**2+DMVA(I,J,2)**2)/(RSI(i,j)
     *           *focean(i,j)*DTSRC*RHOWS))
            UI2rho(I,J)=rhows*(oi_ustar0*max(1d0,1d3*ustar1))**2
          end if
        END DO
        END DO
      END IF

      RETURN
      END SUBROUTINE DYNSI

