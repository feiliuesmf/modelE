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

      subroutine POUT_IJ(TITLE,XIJ,XJ,XSUM)
!@sum  POUT_IJ output lat-lon binary records
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : IM,JM
      USE DAGCOM, only : iu_ij
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, INTENT(IN) :: TITLE*80
!@var XIJ lat/lon output field 
      REAL*4, DIMENSION(IM,JM), INTENT(IN) :: XIJ
!@var XJ lat sum/mean of output field 
      REAL*4, DIMENSION(JM), INTENT(IN) :: XJ
!@var XSUM global sum/mean of output field 
      REAL*4, INTENT(IN) :: XSUM

      WRITE(iu_ij) TITLE,XIJ,XJ,XSUM
      return
      end

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
      REAL*4, DIMENSION(JM+3,LM+LM_REQ+1), INTENT(IN) :: XJL
!@var PM pressure levels (MB)
      REAL*8, DIMENSION(LM+LM_REQ), INTENT(IN) :: PM 
      
      CHARACTER*16, INTENT(IN) :: CX,CY
      CHARACTER*16, PARAMETER :: CBLANK = ' '
      REAL*4 XCOOR(JM)
      INTEGER J,L,JXMAX
      
      XCOOR(1:JM-J1) = LAT_DG(J1:JM,J1)
      JXMAX = JM-J1+1

      WRITE (iu_jl) TITLE,JXMAX,KLMAX,1,1,
     *     ((XJL(J1+J-1,L),J=1,JXMAX),L=1,KLMAX),(XCOOR(J),J=1,JXMAX)
     *     ,(SNGL(PM(L)),L=1,KLMAX),1.,1.,CX,CY,CBLANK,CBLANK,'NASAGISS'
     *     ,(XJL(J,LM+LM_REQ+1),J=J1,JM+3),((XJL(J,L),J=JM+1,JM+3),L=1
     *     ,KLMAX)

      return
      end

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
      REAL*4, DIMENSION(IM,LM+LM_REQ+1), INTENT(IN) :: XIL
!@var PM pressure levels (MB)
      REAL*8, DIMENSION(LM+LM_REQ), INTENT(IN) :: PM 
!@var ASUM vertical mean/sum
      REAL*4, DIMENSION(IM), INTENT(IN) :: ASUM 
!@var GSUM total sum/mean 
      REAL*4, INTENT(IN) :: GSUM
!@var ZONAL zonal sum/mean 
      REAL*4, DIMENSION(LM+LM_REQ), INTENT(IN) :: ZONAL 
      
      CHARACTER*16, INTENT(IN) :: CX,CY
      CHARACTER*16, PARAMETER :: CBLANK = ' '
      REAL*4 XCOOR(IM)
      INTEGER I,L
      
      XCOOR(1:IM) = LON_DG(1:IM,ISHIFT)
      WRITE (iu_il) TITLE,IM,KLMAX,1,1,
     *     ((XIL(I,L),I=1,IM),L=1,KLMAX),(XCOOR(I),I=1,IM)
     *     ,(SNGL(PM(L)),L=1,KLMAX),0.,0.,CX,CY,CBLANK,CBLANK,'NASAGISS'
     *     ,(ASUM(I),I=1,IM),GSUM,(ZONAL(L),L=1,KLMAX)

      return
      end

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
      REAL*4, DIMENSION(JM+3,KAJ), INTENT(IN) :: BUDG
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

