      program vegtypet2b
!@sum converting txt file vegtype.global to binary. 
!@+ In addtion we adopt a new format which is compatible with 
!@+ what's used elsewhere in the model:
!@+ do K-1:VEGTYPE
!@+ read(iu_data) TITLE_VEG_TYPE(K), IUSE_glob(1:IM,1:JM,K)
!@! enddo 
!@  note that the upper bound is now NVEGTYPE. 
!@+ The main change is that now the file contains many entries where IUSE_glob=0.
!@+ Add explanation about ordering of veg. types
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
      real*4, dimension(IM,JM,NVEGTYPE) :: FUSE_glob,order   !fraction btw 0. and 1., index for veg. type ordering
      real*8, dimension(im,jm,ntype) :: XLAI_glob
      real*4, dimension(im,jm,nvegtype) :: XLAI_out

      INTEGER :: I,J,K,L,iu_data,iveg
      character*80 :: fname,oname,TITLE_VEG_TYPE,TITLE
      character(len=2) :: vegtype,c2month 
      character(len=1) :: c1month 

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
            do K=1,IREG_glob(I,J)
               FUSE_glob(I,J,ILAND_GLOB(I,J,K)+1)=
     &              IUSE_glob(I,J,K)/1000.
c               order(I,J,ILAND_GLOB(I,J,K)+1)=K
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
c      do iveg=0,NVEGTYPE-1
c         write(vegtype,'(i2)') iveg
c         TITLE_VEG_TYPE="order type="//vegtype
c         write(3200) TITLE_VEG_TYPE,order(:,:,iveg+1)
c      enddo


c***  Check that fractions at a particular cell add up to 1
      if (any(abs(sum(FUSE_glob,3)-1.0) .gt. 0.001)) then 
         write(6,*) "Warning veg. fractions at cell:",i,j,
     &        "do not add up to 1. Difference is more than 1 per mil"
      endif

c      write(6,*) sum(FUSE_glob,3)
 
      close(unit=3200)

C     Read lai:
      do imonth=1,12
         if (imonth .lt. 10) then
            write(c1month,'(i1)') imonth
            c2month="0"//c1month
         else
            write(c2month,'(i2)') imonth
         endif
         fname='lai'//c2month//'.global'
         write(*,*) fname
         open(20,FILE=fname,FORM="FORMATTED",status="UNKNOWN")
 10      READ(20,"(3I3,20F5.1)",END=20) I,J,INDEX,
     &        (XLAI_glob(I,J,K),K=1,INDEX)
         GOTO 10
 20     close(20)
        
        XLAI_out(:,:,:)=0.
        
        do I=1,IM
           do J=1,JM
              do K=1,IREG_glob(I,J)
                 XLAI_out(I,J,ILAND_GLOB(I,J,K)+1)=
     &                XLAI_glob(I,J,K)
              enddo
           enddo
        enddo
        
        oname=trim(fname)//".bin"
        write(*,*) oname
        open(unit=3200, FILE=oname,FORM='unformatted',
     &       STATUS='unknown')
        
        do iveg=0,NVEGTYPE-1
           write(vegtype,'(i2)') iveg
           TITLE_VEG_TYPE="LAI corresponding to veg type="//vegtype
           write(3200) TITLE_VEG_TYPE,XLAI_out(:,:,iveg+1)
        enddo
        
        close(3200)

      enddo

      end program vegtypet2b
