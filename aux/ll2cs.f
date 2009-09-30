      module regrid
!@auth D,Gueyffier
!@ver 1.0
!@sum    ll2cs is a generic tool to convert latlon input files to the cubed sphere
!@+  it is used together with companion script remap.pl and regridding parameter file
!@+
!@usage  to compile: cd aux; gmake regrid 
!@+  to run : ./remap.pl -par parameter_file -in input_file -out output_file
!@+  where parameter_file is e.g. regrid.par
!@+  and input_file is e.g. top_index_360x180.ij.ext
!@+  output_file is the resulting cubed sphere file
!@+
!@+  The parameter file contains the following entries
!@+  filedir      =  path to input file (e.g. /discover/nobackup/projects/giss/prod_input_files/)
!@+  regridfile   =  remapping file containing interpolation weights (e.g. remap360-180C90-90.nc)
!@+  imsource     =  im resolution of source (latlon) grid
!@+  jmsource     =  jm resolution of source (latlon) grid
!@+  ntilessource =  number of faces of source (latlon) grid (=1)
!@+  imtarget     =  im resolution of target (cubed sphere) grid
!@+  jmtarget     =  jm resolution of target (cubed sphere) grid
!@+  ntilestarget =  number of faces of target (cubed sphere) grid (=6)
!@+  format       =  'ij' : each record contains title, data(imsource,jmsource)
!@+                  'ijl': each record contains title, data(imsource,jmsource,levels)
!@+                  'lij': each record contains title, data(levels,imsouce,jmsource)
!@+  nfields      =  sets the number of fields per record
!@+                         i.e.  title, data_1(im,jm), data_2(im,jm)..data_nfields(im,jm)
!@+  levels (optional) = for 'ijl' or 'lij' arrays only, sets the number of "l" levels 
!@+  title        =  defines if each record contains a title or not
!@+  maintile     =  defines if the first record is a title with no data
!@+
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
     &     cimt,cjmt,cntt,cformat,cfields,fsource,ftarget,ctitle,
     &     clevels,cmaintitle
      character*80 :: title
      real*4, allocatable :: sij4(:,:,:,:) , tij4(:,:,:,:)   ! ij arrays
      real*8, allocatable :: sij(:,:,:,:)  , tij(:,:,:,:)
      real*4, allocatable :: sijl4(:,:,:,:) , tijl4(:,:,:,:)  !ijl arrays
      real*8, allocatable :: sijl(:,:,:,:)  , tijl(:,:,:,:)
      real*4, allocatable :: slij4(:,:,:,:) , tlij4(:,:,:,:)  !lij arrays
      real*8, allocatable :: slij(:,:,:,:)  , tlij(:,:,:,:)
      integer :: ims,jms,nts,imt,jmt,ntt,nfields,nlevels,n,
     &     nargs,maxrec,iuin,iuout
      type (x_2grids) :: x2grids

      nargs = IARGC()
      IF(nargs.lt.14) write(*,*) "ll2cs needs 14 arguments";

      call getarg(1,filesource)
      call getarg(2,filetarget)
      call getarg(3,regridfile)
      call getarg(4,cims)
      call getarg(5,cjms)
      call getarg(6,cnts)
      call getarg(7,cimt)
      call getarg(8,cjmt)
      call getarg(9,cntt)
      call getarg(10,cformat)
      call getarg(11,cfields)
      call getarg(12,ctitle)
      call getarg(13,clevels)
      call getarg(14,cmaintitle)
      read(cims,'(I4)') ims
      read(cjms,'(I4)') jms
      read(cnts,'(I4)') nts
      read(cimt,'(I4)') imt
      read(cjmt,'(I4)') jmt
      read(cntt,'(I4)') ntt
      read(cfields,'(I4)') nfields
      read(clevels,'(I4)') nlevels


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

      cformat = trim(cformat)
      ctitle = trim(ctitle)
      cmaintitle=trim(cmaintitle)

      if (cmaintitle .eq. 'yes' .or. cmaintitle .eq. 'YES'
     &        .or. cmaintitle .eq. 'y' .or. cmaintitle .eq. 'Y') then
         read(unit=iuin) title
         write(unit=iuout) title
      endif

      if (cformat .eq. 'ij' .or. cformat .eq. 'IJ' 
     &     .or. cformat .eq. 'giss' .or. cformat .eq. 'GISS') then

         allocate(sij(ims,jms,nts,nfields),sij4(ims,jms,nts,nfields))
         allocate(tij(imt,jmt,ntt,nfields),tij4(imt,jmt,ntt,nfields))         
         
         if (ctitle .eq. 'yes' .or. ctitle .eq. 'YES'
     &        .or. ctitle .eq. 'y' .or. ctitle .eq. 'Y') then
            do 
               read(unit=iuin,end=30) title,sij4
               sij=sij4
               do n=1,nfields
                  call do_regrid(x2grids,sij(:,:,:,n),tij(:,:,:,n))
               enddo
               tij4=tij
               write(unit=iuout) title,tij4
               maxrec=maxrec+1
            enddo
 30      continue

         else if (ctitle .eq. 'no' .or. ctitle .eq. 'NO'
     &        .or. ctitle .eq. 'n' .or. ctitle .eq. 'N') then
            do 
               read(unit=iuin,end=40) sij4
               sij=sij4
               do n=1,nfields
                  call do_regrid(x2grids,sij(:,:,:,n),tij(:,:,:,n))
               enddo
               tij4=tij
               write(unit=iuout) tij4
               maxrec=maxrec+1
            enddo
 40      continue
         endif

         deallocate(sij4, sij) 
         deallocate(tij4, tij)

         
      else if (cformat .eq. 'ijl' .or. cformat .eq. 'IJL') then
         write(*,*) "IJL"            
         allocate(sijl(ims,jms,nlevels,nts),sijl4(ims,jms,nlevels,nts))
         allocate(tijl(imt,jmt,nlevels,ntt),tijl4(imt,jmt,nlevels,ntt))         
         
         if (ctitle .eq. 'yes' .or. ctitle .eq. 'YES'
     &        .or. ctitle .eq. 'y' .or. ctitle .eq. 'Y') then
            do 
               read(unit=iuin,end=50) title,sijl4
               sijl=sijl4
               do n=1,nlevels
                  call do_regrid(x2grids,sijl(:,:,n,:),tijl(:,:,n,:))
               enddo
               tijl4=tijl
               write(unit=iuout) title,tijl4
               maxrec=maxrec+1
            enddo
 50      continue

         else if (ctitle .eq. 'no' .or. ctitle .eq. 'NO'
     &        .or. ctitle .eq. 'n' .or. ctitle .eq. 'N') then
            do 
               read(unit=iuin,end=60) sijl4
               sijl=sijl4
               do n=1,nlevels
                  call do_regrid(x2grids,sijl(:,:,n,:),tijl(:,:,n,:))
               enddo
               tijl4=tijl
               write(unit=iuout) tijl4
               maxrec=maxrec+1
            enddo
 60      continue
         endif

         deallocate(sijl4, sijl) 
         deallocate(tijl4, tijl)
         
      else if (cformat .eq. 'lij' .or. cformat .eq. 'LIJ') then
         
         allocate(slij(nlevels,ims,jms,nts),slij4(nlevels,ims,jms,nts))
         allocate(tlij(nlevels,imt,jmt,ntt),tlij4(nlevels,imt,jmt,ntt))         
         
         if (ctitle .eq. 'yes' .or. ctitle .eq. 'YES'
     &        .or. ctitle .eq. 'y' .or. ctitle .eq. 'Y') then
            do 
               read(unit=iuin,end=70) title,slij4
               slij=slij4
               do n=1,nlevels
                  call do_regrid(x2grids,slij(n,:,:,:),tlij(n,:,:,:))
               enddo
               tlij4=tlij
               write(unit=iuout) title,tlij4
               maxrec=maxrec+1
            enddo
 70         continue

         else if (ctitle .eq. 'no' .or. ctitle .eq. 'NO'
     &           .or. ctitle .eq. 'n' .or. ctitle .eq. 'N') then
            do 
               read(unit=iuin,end=80) slij4
               slij=slij4
               do n=1,nlevels
                  call do_regrid(x2grids,slij(n,:,:,:),tlij(n,:,:,:))
               enddo
               tlij4=tlij
               write(unit=iuout) tlij4
               maxrec=maxrec+1
            enddo
 80         continue
         endif
         
         deallocate(slij4, slij) 
         deallocate(tlij4, tlij)
         
      endif

      close(iuin)
      close(iuout)

 
      write(6,*) "remapped file contains:"
      write(6,100) maxrec," records"
      write(6,200) "record format ",trim(cformat)
      write(6,200) trim(cfields)," fields per record"
      write(6,200) trim(clevels)," levels per array"
      write(6,200) "records contains title? ", ctitle
      write(6,200) "file contains main title? ", cmaintitle

 100  format(3X, I2, A)
 200  format(4X, A, A)


      end program ll2cs



