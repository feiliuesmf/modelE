      program glmeltt2b
      implicit none
      integer, parameter :: IM=72,JM=46,ncols=72
      character*72 :: TITLE
      character*80 :: TITLE80,fname,oname
      CHARACTER*1, DIMENSION(IM,JM) :: CGLM   ! global array
      integer :: iu_GL,Irow,i,j,len
      real*4, allocatable :: RGLM(:,:)

      iu_GL=20
      fname="GLMELT_4X5.OCN"

      open(unit=iu_GL,FILE=fname,FORM ="FORMATTED",
     &     status ="UNKNOWN")

      READ  (iu_GL,100) TITLE
  100 FORMAT(A<ncols>)

      WRITE (6,*) 'Read on unit ',iu_GL,': ',TITLE
      READ  (iu_GL,*)

C**** assumes a fix width slab - will need adjusting for CS
      DO Irow=1,1+(IM-1)/ncols
        DO J=JM,1,-1
          READ (iu_GL,'(72A1)')
     *         (CGLM(I,J),I=ncols*(Irow-1)+1,MIN(IM,Irow*ncols))
        END DO
      END DO
  200 FORMAT(<ncols>A1)
C****

      allocate(RGLM(IM,JM))

      do i=1,IM
         do j=1,JM
            if (CGLM(i,j) .eq. "1") then
               RGLM(i,j)=1.0
            else
               RGLM(i,j)=0
            endif
         enddo
      enddo

      write(*,*) "RGLM=",RGLM
      close(iu_GL)

      oname=trim(fname)//".bin"
      write(*,*) oname
      write(*,*) "size iglm",size(RGLM)
      open(unit=3200, FILE=oname,FORM='unformatted',STATUS='unknown')


      TITLE80=TITLE    
      TITLE80=trim(TITLE80)
      write(3200) TITLE80,RGLM
      close(unit=3200)
      deallocate(RGLM)

      end program glmeltt2b
