      module GHYCOM

      USE E001M12_COM
     &     , only : im,jm
      USE SLE001
     &     , only : ngm,imt
      IMPLICIT NONE

C bare/veg not in merged array because WBARE does not contain
C 0 index for legacy reasons
      DOUBLE PRECISION, DIMENSION(IM,JM,  NGM) :: WBARE
      DOUBLE PRECISION, DIMENSION(IM,JM,0:NGM) :: WVEGE
      DOUBLE PRECISION, DIMENSION(IM,JM,0:NGM) :: HTBARE
      DOUBLE PRECISION, DIMENSION(IM,JM,0:NGM) :: HTVEGE
      DOUBLE PRECISION, DIMENSION(IM,JM,2) :: SNOWBV
      common/soils3/ wbare,wvege,htbare,htvege,snowbv
C handle for reading/writing ghdata
      DOUBLE PRECISION, DIMENSION(IM,JM,4*NGM+5) :: GHDATA
      equivalence (ghdata(1,1,1),wbare(1,1,1))
      DOUBLE PRECISION, DIMENSION(IM,JM,NGM) :: DZ_IJ
      DOUBLE PRECISION, DIMENSION(IM,JM,IMT,NGM) :: Q_IJ,QK_IJ
      DOUBLE PRECISION, DIMENSION(IM,JM) :: SL_IJ
C this common only for purpose of reading its contents from a
C file opened in fortran unformatted sequential access mode
C containing its contents in a contiguous real*4 block
      COMMON/SDATA/ DZ_IJ,Q_IJ,QK_IJ,SL_IJ

      DOUBLE PRECISION, DIMENSION(IM,JM,NGM) :: AFR
      DOUBLE PRECISION, DIMENSION(3,IM,JM) :: ALA,ACS
      DOUBLE PRECISION, DIMENSION(IM,JM) :: AFB,AVH

      end module GHYCOM
