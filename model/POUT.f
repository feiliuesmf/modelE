!@sum  POUT default output routines for standard GISS formats
!@auth Gavin Schmidt
!@ver  1.0

C****
C**** If other formats are desired, please replace these routines
C**** with ones appropriate for your chosen output format, but with the
C**** same interface
C****
C**** Note: it would be nice to amalgamate IL and JL, but that will
C**** have to wait.

      subroutine open_ij(filename)
!@sum  OPEN_IJ opens the lat-lon binary output file
!@auth M. Kelley
!@ver  1.0
      USE DAGCOM, only : iu_ij
      USE FILEMANAGER
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename

      call openunit(filename,iu_ij,.true.,.false.)
      return
      end subroutine open_ij

      subroutine close_ij
!@sum  CLOSE_IJ closes the lat-lon binary output file
!@auth M. Kelley
!@ver  1.0
      USE DAGCOM, only : iu_ij
      USE FILEMANAGER
      IMPLICIT NONE
      call closeunit(iu_ij)
      return
      end subroutine close_ij

      subroutine POUT_IJ(TITLE,XIJ,XJ,XSUM,IJGRID)
!@sum  POUT_IJ output lat-lon binary records
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : IM,JM
      USE DAGCOM, only : iu_ij
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, INTENT(IN) :: TITLE*80
!@var XIJ lat/lon output field
      REAL*8, DIMENSION(IM,JM), INTENT(IN) :: XIJ
!@var XJ lat sum/mean of output field
      REAL*8, DIMENSION(JM), INTENT(IN) :: XJ
!@var XSUM global sum/mean of output field
      REAL*8, INTENT(IN) :: XSUM
!@var IJGRID = 1 for primary lat-lon grid, 2 for secondary lat-lon grid
      INTEGER, INTENT(IN) :: IJGRID

      WRITE(iu_ij) TITLE,SNGL(XIJ),SNGL(XJ),SNGL(XSUM)
      return
      end

      subroutine open_jl(filename)
!@sum  OPEN_JL opens the lat-height binary output file
!@auth M. Kelley
!@ver  1.0
      USE DAGCOM, only : iu_jl
      USE FILEMANAGER
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename

      call openunit(filename,iu_jl,.true.,.false.)
      return
      end subroutine open_jl

      subroutine close_jl
!@sum  CLOSE_JL closes the lat-height binary output file
!@auth M. Kelley
!@ver  1.0
      USE DAGCOM, only : iu_jl
      IMPLICIT NONE
      close(iu_jl)
      return
      end subroutine close_jl

      subroutine POUT_JL(TITLE,J1,KLMAX,XJL,PM,CX,CY)
!@sum  POUT_JL output lat-height binary records
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : JM,LM
      USE GEOM, only : lat_dg
      USE DAGCOM, only : lm_req,iu_jl
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, INTENT(IN) :: TITLE*80
!@var KLMAX max level to output
!@var J1 minimum j value to output (needed for secondary grid fields)
      INTEGER, INTENT(IN) :: KLMAX,J1
!@var XJL output field
!@+       (J1:JM,1:KLMAX) is field
!@+       (JM+1:JM+3,1:KLMAX) are global/NH/SH average over L
!@+       (J1:JM+3,LM+LM_REQ+1) are averages over J
      REAL*8, DIMENSION(JM+3,LM+LM_REQ+1), INTENT(IN) :: XJL
!@var PM pressure levels (MB)
      REAL*8, DIMENSION(LM+LM_REQ), INTENT(IN) :: PM

      CHARACTER*16, INTENT(IN) :: CX,CY
      CHARACTER*16, PARAMETER :: CBLANK = ' '
      REAL*8 XCOOR(JM)
      INTEGER J,L,JXMAX

      XCOOR(1:JM-J1) = LAT_DG(J1:JM,J1)
      JXMAX = JM-J1+1

      WRITE (iu_jl) TITLE,JXMAX,KLMAX,1,1,
     *     ((SNGL(XJL(J1+J-1,L)),J=1,JXMAX),L=1,KLMAX)
     *     ,(SNGL(XCOOR(J)),J=1,JXMAX)
     *     ,(SNGL(PM(L)),L=1,KLMAX),1.,1.,CX,CY,CBLANK,CBLANK,'NASAGISS'
     *     ,(SNGL(XJL(J,LM+LM_REQ+1)),J=J1,JM+3)
     *     ,((SNGL(XJL(J,L)),J=JM+1,JM+3),L=1,KLMAX)

      return
      end

      subroutine open_il(filename)
!@sum  OPEN_IL opens the lon-height binary output file
!@auth M. Kelley
!@ver  1.0
      USE DAGCOM, only : iu_il
      USE FILEMANAGER
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename

      call openunit(filename,iu_il,.true.,.false.)
      return
      end subroutine open_il

      subroutine close_il
!@sum  CLOSE_IL closes the lon-height binary output file
!@auth M. Kelley
!@ver  1.0
      USE DAGCOM, only : iu_il
      USE FILEMANAGER
      IMPLICIT NONE
      call closeunit(iu_il)
      return
      end subroutine close_il

      subroutine POUT_IL(TITLE,ISHIFT,KLMAX,XIL,PM,CX,CY,
     *     ASUM,GSUM,ZONAL)
!@sum  POUT_IL output lon-height binary records
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : IM,LM
      USE GEOM, only : lon_dg
      USE DAGCOM, only : lm_req,iu_il
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, INTENT(IN) :: TITLE*80
!@var KLMAX max level to output
!@var ISHIFT flag for secondary grid
      INTEGER, INTENT(IN) :: KLMAX,ISHIFT
!@var XIL output field
      REAL*8, DIMENSION(IM,LM+LM_REQ+1), INTENT(IN) :: XIL
!@var PM pressure levels (MB)
      REAL*8, DIMENSION(LM+LM_REQ), INTENT(IN) :: PM
!@var ASUM vertical mean/sum
      REAL*8, DIMENSION(IM), INTENT(IN) :: ASUM
!@var GSUM total sum/mean
      REAL*8, INTENT(IN) :: GSUM
!@var ZONAL zonal sum/mean
      REAL*8, DIMENSION(LM+LM_REQ), INTENT(IN) :: ZONAL

      CHARACTER*16, INTENT(IN) :: CX,CY
      CHARACTER*16, PARAMETER :: CBLANK = ' '
      REAL*8 XCOOR(IM)
      INTEGER I,L

      XCOOR(1:IM) = LON_DG(1:IM,ISHIFT)
      WRITE (iu_il) TITLE,IM,KLMAX,1,1,
     *     ((SNGL(XIL(I,L)),I=1,IM),L=1,KLMAX),(SNGL(XCOOR(I)),I=1,IM)
     *     ,(SNGL(PM(L)),L=1,KLMAX),0.,0.,CX,CY,CBLANK,CBLANK,'NASAGISS'
     *     ,(SNGL(ASUM(I)),I=1,IM),SNGL(GSUM),(SNGL(ZONAL(L)),L=1,KLMAX)

      return
      end

      subroutine close_j
!@sum  CLOSE_J closes the latitudinal budget-page ascii output file
!@auth M. Kelley
!@ver  1.0
      USE DAGCOM, only : iu_j
      USE FILEMANAGER
      IMPLICIT NONE
      call closeunit(iu_j)
      return
      end subroutine close_j

      subroutine open_j(filename)
!@sum  OPEN_J opens the latitudinal budget-page ascii output file
!@auth M. Kelley
!@ver  1.0
      USE DAGCOM, only : iu_j
      USE FILEMANAGER
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename

      call openunit(filename,iu_j,.false.,.false.)
      return
      end subroutine open_j

      subroutine POUT_J(TITLE,BUDG,KMAX,TERRAIN)
!@sum  POUT_J output zonal budget ascii file (aplot format)
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : JM
      USE DAGCOM, only : KAJ,iu_j
      USE GEOM, only : lat_dg
      IMPLICIT NONE
      CHARACTER*16, DIMENSION(KAJ),INTENT(INOUT) :: TITLE
      CHARACTER*16, INTENT(IN) :: TERRAIN
      REAL*8, DIMENSION(JM+3,KAJ), INTENT(IN) :: BUDG
      INTEGER, INTENT(IN) :: KMAX
      INTEGER K,N,J

C**** Convert spaces in TITLE to underscore
      DO K=1,KMAX
        do n=2,len_trim(title(K))    ! skip control character
          if (title(K)(n:n).eq.' ') title(K)(n:n)='_'
        end do
        title(K)(1:1)=' '
      END DO

      WRITE(iu_j,*) "Zonal Budgets for surface type ",TERRAIN
      WRITE(iu_j,*) "Latitude"
      WRITE(iu_j,*) "Zonal Average"
      WRITE(iu_j,*) "Lat",(TRIM(TITLE(K)),K=1,KMAX)

      DO J=1,JM
        WRITE(iu_j,'(I4,100(1X,F8.3))') NINT(LAT_DG(J,1)),
     *       (BUDG(J,K),K=1,KMAX)
      END DO
      WRITE(iu_j,*)
C**** output hemispheric and global means
      WRITE(iu_j,'(A4,100F8.3)') "NH",(BUDG(JM+1,K),K=1,KMAX)
      WRITE(iu_j,'(A4,100F8.3)') "SH",(BUDG(JM+2,K),K=1,KMAX)
      WRITE(iu_j,'(A4,100F8.3)') "GLOB",(BUDG(JM+3,K),K=1,KMAX)
      WRITE(iu_j,*)

      return
      end

      subroutine open_ijk(filename)
!@sum  OPEN_IJK opens the lat-lon-height binary output file
!@auth M. Kelley
!@ver  1.0
      USE DAGCOM, only : iu_ijk
      USE FILEMANAGER
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename

      call openunit(filename,iu_ijk,.true.,.false.)
      return
      end subroutine open_ijk

      subroutine close_ijk
!@sum  CLOSE_IJK closes the lat-lon-height binary output file
!@auth M. Kelley
!@ver  1.0
      USE DAGCOM, only : iu_ijk
      USE FILEMANAGER
      IMPLICIT NONE
      call closeunit(iu_ijk)
      return
      end subroutine close_ijk

      subroutine POUT_IJK(TITLE,XIJK,XJK,XK)
!@sum  POUT_IJK outputs lat-lon-height binary records
!@auth M. Kelley
!@ver  1.0
      USE MODEL_COM, only : IM,JM,LM
      USE DAGCOM, only : iu_ijk
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, DIMENSION(LM), INTENT(IN) :: TITLE*80
!@var XIJK lat/lon/height output field
      REAL*8, DIMENSION(IM,JM-1,LM), INTENT(IN) :: XIJK
!@var XJK lat sum/mean of output field
      REAL*8, DIMENSION(JM-1,LM), INTENT(IN) :: XJK
!@var XK global sum/mean of output field
      REAL*8, DIMENSION(LM), INTENT(IN) :: XK

      INTEGER :: I,K

      DO K=1,LM
         WRITE(iu_ijk) TITLE(K),
c fill in missing first row for GISS format(should put in missing value)
     &        ((0.,i=1,im)),SNGL(XIJK(:,:,K)),
     &        0.           , SNGL(XJK(:,K)),
     &        SNGL(XK(K))
      ENDDO
      return
      end
