      MODULE SOMTQ_COM
!@sum  SOMTQ_COM contains the arrays containing second order moments
!@auth Gary Russell
!@ver  1.0
      USE QUSDEF
      USE MODEL_COM, only : im,jm,lm
      USE DOMAIN_DECOMP, only : grid
      IMPLICIT NONE
      SAVE
!     REAL*8, DIMENSION(NMOM,IM, GRID%J_STRT_HALO:GRID%J_STOP_HALO ,LM)  
!    &        :: TMOM,QMOM
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: TMOM,QMOM

      END MODULE SOMTQ_COM

      SUBROUTINE ALLOC_SMOMTQ(grid)
!@sum  init_smomtq allocates the arrays in this module which
!@+    must now be dynamic for the distributed memory implementation.
!@auth Rosalinda de Fainchtein
!@ver  1.0
      USE DOMAIN_DECOMP, ONLY : DYN_GRID
      USE QUSDEF, ONLY : NMOM
      USE MODEL_COM, ONLY : LM
      USE SOMTQ_COM, ONLY : TMOM,QMOM
      IMPLICIT NONE
      TYPE (DYN_GRID), INTENT(IN) :: grid

      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER

      INTEGER :: I,J,L

      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO

      ALLOCATE ( TMOM(NMOM , I_0H:I_1H , J_0H:J_1H , LM) )
      ALLOCATE ( QMOM(NMOM , I_0H:I_1H , J_0H:J_1H , LM) )

      END SUBROUTINE ALLOC_SMOMTQ

      SUBROUTINE io_somtq(kunit,iaction,ioerr)
!@sum  io_somtq reads and writes second order moments to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,lhead
      USE DOMAIN_DECOMP, only : grid
      USE SOMTQ_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "QUS01"

      write (MODULE_HEADER(lhead+1:80),'(a7,i2,a)')
     * 'R8 dim(',nmom,',im,jm,lm):Tmom,Qmom'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)           ! output to standard restart file
        WRITE (KUNIT,ERR=10) MODULE_HEADER,TMOM,QMOM
      CASE (IOREAD:)            ! input from restart file
        READ (KUNIT,ERR=10) HEADER,TMOM,QMOM
        IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
          PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
          GO TO 10
        END IF
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_somtq

      subroutine tq_zmom_init(t,q)
      USE MODEL_COM, only : im,jm,lm,sige,sig
      USE DOMAIN_DECOMP, ONLY: grid
      USE SOMTQ_COM
      implicit none
      REAL*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lm) :: t,q
      integer :: i,j,l
      REAL*8 :: rdsig

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0S, J_1S, J_0STG, J_1STG
C****
C**** Extract useful local domain parameters from "grid"
C****
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      J_0 = grid%J_STRT
      J_1 = grid%J_STOP
      J_0S = grid%J_STRT_SKP
      J_1S = grid%J_STOP_SKP
      J_0STG = grid%J_STRT_STGR
      J_1STG = grid%J_STOP_STGR

C**** INITIALIZES VERTICAL SLOPES OF T,Q
      DO J=J_0,J_1
      DO I=1,IM
         RDSIG=(SIG(1)-SIGE(2))/(SIG(1)-SIG(2))
         TMOM(MZ,I,J,1)=(T(I,J,2)-T(I,J,1))*RDSIG
         QMOM(MZ,I,J,1)=(Q(I,J,2)-Q(I,J,1))*RDSIG
         IF(Q(I,J,1)+QMOM(MZ,I,J,1).LT.0.) QMOM(MZ,I,J,1)=-Q(I,J,1)
         DO L=2,LM-1
            RDSIG=(SIG(L)-SIGE(L+1))/(SIG(L-1)-SIG(L+1))
            TMOM(MZ,I,J,L)=(T(I,J,L+1)-T(I,J,L-1))*RDSIG
            QMOM(MZ,I,J,L)=(Q(I,J,L+1)-Q(I,J,L-1))*RDSIG
            IF(Q(I,J,L)+QMOM(MZ,I,J,L).LT.0.) QMOM(MZ,I,J,L)=-Q(I,J,L)
         END DO
         RDSIG=(SIG(LM)-SIGE(LM+1))/(SIG(LM-1)-SIG(LM))
         TMOM(MZ,I,J,LM)=(T(I,J,LM)-T(I,J,LM-1))*RDSIG
         QMOM(MZ,I,J,LM)=(Q(I,J,LM)-Q(I,J,LM-1))*RDSIG
         IF(Q(I,J,LM)+QMOM(MZ,I,J,LM).LT.0.) QMOM(MZ,I,J,LM)=-Q(I,J,LM)
      END DO
      END DO
      return
      end subroutine tq_zmom_init
