      program vegtypet2b
!@sum converting txt file vegtype.global to binary. 
!@+ In addtion we adopt a new format which is compatible with 
!@+ what's used elsewhere in the model:
!@+ do K-1:VEGTYPE
!@+ read(iu_data) TITLE_VEG_TYPE(K), IUSE_glob(1:IM,1:JM,K)
!@! enddo 
!@  note that the upper bound is now NVEGTYPE. 
!@+ The main change is that now the file contains many entries where IUSE_glob=0.
!@+
!@+ The arrays IREG_glob and ILAND_GLOB are deduced from
!@+ IUSE_glob by checking which values of IREG_glob are nonzero
!@+ 
!@+ The new binary file vegtype.global is used in a separate routine (regridinput) to regrid 
!@+ the fractions IUSE_glob to the cubed-sphere in such a way that sum(IUSE(I0,J0,:),3) = 1000 
!@+ for any (I0,J0) on the cubed-sphere surface.
!@+ 
!@auth Denis Gueyffier
!@+ 
      integer, parameter :: IM=72, JM=46
      INTEGER, PARAMETER :: NPOLY   = 20,
     &                      NTYPE   = 16,
     &                      NVEGTYPE= 74
      INTEGER, PARAMETER :: NWAT=6
      integer, dimension(IM,JM) :: IREG_glob
      integer, dimension(IM,JM,NTYPE) :: ILAND_glob,IUSE_glob
      real*4, dimension(IM,JM,NVEGTYPE) :: FUSE_glob   !fraction btw 0. and 1.
      CHARACTER*1, dimension(70) :: COM
      character(len=300) :: out_line
      INTEGER, DIMENSION(NVEGTYPE) :: IZO
      INTEGER :: I,J,K,L,NWAT2,IDUMMY,iu_data,iols,iveg
      INTEGER, DIMENSION(NWAT) :: IWATER
      character*80 :: fname,oname,TITLE_VEG_TYPE,TITLE
      character(len=2) :: vegtype 

      write(6,*) 'READING land types and fractions'
      write(6,*) 'and converting to binary'
      fname='vegtype.global'
      open(20,FILE=fname,FORM="FORMATTED",status="UNKNOWN")

 100  READ(20,'(20I4)',end=110) I,J,IREG_glob(I,J),
     &     (ILAND_glob(I,J,K),K=1,IREG_glob(I,J)),
     &     (IUSE_glob(I,J,K),K=1,IREG_glob(I,J))
      GO TO 100
 110  CONTINUE
      close(20)

      FUSE_glob(:,:,:)=0

c***  Warning : ILAND_GLOB takes values between 0 and 73, and 0 isn't a friendly bound 
c***  for fortran arrays
      do I=1,IM
         do J=1,JM
            write(*,*) "ireg=",IREG_glob(I,J)
            do K=1,IREG_glob(I,J)
               FUSE_glob(I,J,ILAND_GLOB(I,J,K)+1)=
     &              IUSE_glob(I,J,K)/1000.
            enddo
         enddo
      enddo
      
      
      oname=trim(fname)//".bin"
      write(*,*) oname

      open(unit=3200, FILE=oname,FORM='unformatted',STATUS='unknown')

c      TITLE="vegetation fractions, 4x5 grid, please add more info here"
c      write(3200) TITLE

      do iveg=0,NVEGTYPE-1
         write(vegtype,'(i2)') iveg
         TITLE_VEG_TYPE="veg type="//vegtype
         write(3200) TITLE_VEG_TYPE,FUSE_glob(:,:,iveg+1)
      enddo

c***  Check that fractions at a particular cell add up to 1
      if (any(abs(sum(FUSE_glob,3)-1.0) .gt. 0.001)) then 
         write(6,*) "Warning veg. fractions at cell:",i,j,
     &        "do not add up to 1. Difference is more than 1 per mil"
      endif

c      write(6,*) sum(FUSE_glob,3)
 
      close(unit=3200)

      end program vegtypet2b
