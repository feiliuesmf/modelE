      module regrid
      private
      public :: x_2grids, init_regrid, do_regrid

      type x_grid      ! stores exchange grid information 
      integer  :: ncells
      integer, allocatable, dimension(:,:) :: ijcub
      integer, allocatable, dimension(:,:) :: ijlatlon
      real*8, allocatable, dimension(:) :: xgrid_area
      integer, allocatable, dimension(:) :: tile
      integer :: maxkey      
      end type x_grid

      type x_2grids   ! stores x-grid info plus source and target grid info 
      type (x_grid) :: xgrid
      integer :: ntilessource  ! #tiles of source grid (1 for latlon, 6 for cubed sphere)
      integer :: ntilestarget  ! #tiles of target grid (1 for latlon, 6 for cubed sphere)
      integer :: imsource      ! im for source grid
      integer :: jmsource      ! jm for source grid
      integer :: imtarget      ! im for target grid
      integer :: jmtarget      ! jm for target grid
      end type x_2grids

      contains 


      subroutine init_regrid(fname,x2grids,imsource,jmsource,
     &     ntilessource,imtarget,jmtarget,ntilestarget)

!@sum  Reads regriding file then instanciates the x_2grids derived type (x_2grids type=x_grid plus info about
!@+    source and target grids). 
!@auth Denis Gueyffier
      implicit none
      include 'netcdf.inc'

      type (x_2grids) :: x2grids
      type (x_grid) :: xgrid
      real*8,  allocatable, dimension(:) :: xgrid_area  !local variable
      integer, allocatable, dimension(:) :: tile        !local variable
      integer, allocatable, dimension(:,:) :: ijcub     !local variable
      integer, allocatable, dimension(:,:) :: ijlatlon  !local variable
      integer :: ncells                                 !local variable
      integer, intent(in) :: imsource,jmsource,imtarget,jmtarget
      integer, intent(in) :: ntilessource,ntilestarget
      integer :: status,fid,n,vid,ikey,jlat
      integer :: itile,j,idomain,iic,jjc,index,indexc,nc2
      integer :: ierr,i
      character(len=30) :: fname

      x2grids%imsource=imsource
      x2grids%jmsource=jmsource
      x2grids%ntilessource=ntilessource
      x2grids%imtarget=imtarget
      x2grids%jmtarget=jmtarget
      x2grids%ntilestarget=ntilestarget

c     
c     Read weights
c     
      status = nf_open(trim(fname),nf_nowrite,fid)
      
      if (status .ne. NF_NOERR) write(*,*) 
     *     "UNABLE TO OPEN REMAP FILE",trim(fname)
      
      status = nf_inq_dimid(fid,'ncells',vid)
      status = nf_inq_dimlen(fid,vid,ncells)

     

      allocate(xgrid_area(ncells))
      allocate(ijcub(2,ncells))
      allocate(ijlatlon(2,ncells))
      allocate(tile(ncells))

      status = nf_inq_varid(fid,'xgrid_area',vid)
      status = nf_get_var_double(fid,vid,xgrid_area)
      status = nf_inq_varid(fid,'tile1',vid)
      status = nf_get_var_int(fid,vid,tile)
      status = nf_inq_varid(fid,'tile1_cell',vid)
      status = nf_get_var_int(fid,vid,ijcub)
      status = nf_inq_varid(fid,'tile2_cell',vid)
      status = nf_get_var_int(fid,vid,ijlatlon)
      status = nf_close(fid)
    
      allocate(x2grids%xgrid%xgrid_area(ncells))
      allocate(x2grids%xgrid%ijcub(2,ncells))
      allocate(x2grids%xgrid%ijlatlon(2,ncells))
      allocate(x2grids%xgrid%tile(ncells))
      
      x2grids%xgrid%xgrid_area(:)=xgrid_area(:)
      x2grids%xgrid%ijcub(:,:)=ijcub(:,:)
      x2grids%xgrid%ijlatlon(:,:)=ijlatlon(:,:)
      x2grids%xgrid%tile(:)=tile(:)
      x2grids%xgrid%ncells=ncells

      deallocate(xgrid_area)
      deallocate(ijcub)
      deallocate(ijlatlon)
      deallocate(tile)

      end subroutine init_regrid

      subroutine do_regrid(x2grids,tsource,ttarget)
!@sum  regrid data from source grid to target grid
!@auth Denis Gueyffier
      implicit none
      type (x_2grids), intent(in) :: x2grids
      real*8 :: tsource(x2grids%imsource,x2grids%jmsource,
     &     x2grids%ntilessource)
      real*8 :: ttarget(x2grids%imtarget,x2grids%jmtarget,
     &     x2grids%ntilestarget)
     &     ,atarget(x2grids%imtarget,x2grids%jmtarget,
     &     x2grids%ntilestarget)
      integer :: n,icub,jcub,i,j,itile,icc,jcc,il,jl


      if ((x2grids%ntilessource .eq. 1) .and.
     &     (x2grids%ntilestarget .eq. 6)) then

c     ll2cs

         atarget(:,:,:) = 0.d0
         ttarget(:,:,:) = 0.d0


         do n=1,x2grids%xgrid%ncells
            itile=x2grids%xgrid%tile(n)
            icc=x2grids%xgrid%ijcub(1,n)
            jcc=x2grids%xgrid%ijcub(2,n)
            il=x2grids%xgrid%ijlatlon(1,n)
            jl=x2grids%xgrid%ijlatlon(2,n)

            atarget(icc,jcc,itile) = atarget(icc,jcc,itile)
     &           + x2grids%xgrid%xgrid_area(n)
            ttarget(icc,jcc,itile) = ttarget(icc,jcc,itile)
     &           + x2grids%xgrid%xgrid_area(n)
     &           *tsource(il,jl,1)
         enddo

         do itile=1,x2grids%ntilestarget
            do j=1,x2grids%jmtarget
               do i=1,x2grids%imtarget
                  ttarget(i,j,itile) = ttarget(i,j,itile)
     &                 /atarget(i,j,itile)
               enddo
            enddo
         enddo

      endif

      end subroutine do_regrid

      end module regrid
 
      program ll2cs
!@sum Program to regrid latlon input files to the cubed sphere grid. 
!@+   Program is standalone and runs in serial mode only 
!@+
!@+   The latlon input file has to conform to the GISS file format:
!@+   char*80 title, real*4 data
!@+
!@+
!@auth Denis Gueyffier
      use regrid, only : init_regrid, x_2grids, do_regrid
      implicit none
      character*150 :: filesource,filetarget,regridfile,cims,cjms,cnts,
     &     cimt,cjmt,cntt,fsource,ftarget
      character*80 :: title
      real*4, allocatable :: tsource4(:,:,:), ttarget4(:,:,:)
      real*8, allocatable :: tsource(:,:,:),ttarget(:,:,:)
      integer :: ims,jms,nts,imt,jmt,ntt,nargs,maxrec,iuin,iuout
      type (x_2grids) :: x2grids

      NARGS = IARGC()
      IF(NARGS.lt.9) write(*,*) "function needs 9 arguments";

      call getarg(1,filesource)
      call getarg(2,filetarget)
      call getarg(3,regridfile)
      call getarg(4,cims)
      call getarg(5,cjms)
      call getarg(6,cnts)
      call getarg(7,cimt)
      call getarg(8,cjmt)
      call getarg(9,cntt)
      read(cims,'(I4)') ims
      read(cjms,'(I4)') jms
      read(cnts,'(I4)') nts
      read(cimt,'(I4)') imt
      read(cjmt,'(I4)') jmt
      read(cntt,'(I4)') ntt


c***  Initialize exchange grid
      call init_regrid(regridfile,x2grids,ims,jms,nts,imt,jmt,ntt)


c***  Read data (must be in GISS format)
      fsource=trim(filesource)
      ftarget=trim(filetarget)

      iuin=200
      iuout=300
      maxrec=0

      open(iuin,FILE=fsource,FORM='unformatted', STATUS='old')
      open(iuout,FILE=ftarget,FORM='unformatted', STATUS='unknown')

      allocate( tsource(ims,jms,nts), tsource4(ims,jms,nts))
      allocate( ttarget(imt,jmt,ntt), ttarget4(imt,jmt,ntt))

      do 
         read(unit=iuin,end=30) title,tsource4
         tsource=tsource4
         call do_regrid(x2grids,tsource,ttarget)
         ttarget4=ttarget
         write(unit=iuout) title,ttarget4
         maxrec=maxrec+1
      enddo
 30   continue

      write(6,*) "We have remapped ",maxrec," records"

      deallocate(tsource4, tsource) 
      deallocate(ttarget4, ttarget)

      close(iuin)
      close(iuout)

      end program ll2cs



