      module regrid
!@auth D,Gueyffier
!@ver 1.0
!@sum    ll2csnc is a generic tool to convert netcdf latlon input files to the cubed sphere
!@+  it is used together with companion script remap.pl and regridding parameter file
!@+
!@usage  to compile: cd aux; gmake regrid 
!@+  to run : ./ncremap.pl -par parameter_file -in input_file -out output_file.nc
!@+  where parameter_file is e.g. nvregrid.par
!@+  and input_file is e.g. 
!@+  output_file.nc is the resulting cubed sphere netcdf file
!@+
!@+  The parameter file contains the following entries
!@+  filedir      =  path to input file (e.g. /discover/nobackup/projects/giss/prod_input_files/)
!@+  regridfile   =  remapping file containing interpolation weights (e.g. remap720-360C90-90.nc)
!@+  ntilessource =  number of faces of source (latlon) grid (=1)
!@+  ntilestarget =  number of faces of target (cubed sphere) grid (=6)
!@+  imtarget     =  im resolution of target (cubed sphere) grid
!@+  jmtarget     =  jm resolution of target (cubed sphere) grid
!@+  format       =  format of the variable that the user wants to regrid (e.g. time,lat,lon)
!@+                  note that with netcdf the dimensions are declared in reverse order 
!@+                  ( == C ordering == row major) with respect to fortran ( == column major)
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

      write(*,*) "initialized regridding weights"

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
 
      program ncll2cs
!@sum Program to regrid latlon input files to the cubed sphere grid. 
!@+   Program is standalone and runs in serial mode only 
!@+
!@auth Denis Gueyffier
      use regrid, only : init_regrid, x_2grids, do_regrid
      use ncio, only : defvar, write_data
      implicit none
      include 'netcdf.inc'
      character*150 :: filesource,filetarget,regridfile,cnts,
     &     cimt,cjmt,cntt,cformat,ccopy,cidl
      character*80 :: title
      real*4, allocatable :: sij4(:,:,:) , sij4shift(:,:,:) 
     &     , tij4(:,:,:)   ! ij arrays
      real*8, allocatable :: sij(:,:,:)  , tij(:,:,:)
      real*4, allocatable :: sijl4(:,:,:,:), sijl4shift(:,:,:,:) 
     &     , tijl4(:,:,:,:)     !ijl arrays
      real*8, allocatable :: sijl(:,:,:,:)  , tijl(:,:,:,:)
      real*4, allocatable :: slij4(:,:,:,:) , tlij4(:,:,:,:)  !lij arrays
      real*8, allocatable :: slij(:,:,:,:)  , tlij(:,:,:,:)
      integer :: ims,jms,lms,kms,nts,imt,jmt,ntt,nfields,nlevels,n,
     &     nargs,maxrec,iuin,iuout,i,imax,status,fid,fidt,vid
      type (x_2grids) :: x2grids
      integer :: offset(10),dim(10)
      integer :: ndims, nvars, ngatts, itype, ndim, natt, UNLIMDIMID,
     &     idlms,idjms,idims,is
      character*20 :: cdim(10),cform
      character*80 :: cval,vname
      integer ::  ishape(3)   ! variable shape
      integer :: dtype,nd,shp(7)
      logical :: r4_on_disk !real*8 stored as real*4
      

      nargs = IARGC()
      IF(nargs.lt.9) write(*,*) "ncll2cs needs 9 arguments";

      call getarg(1,filesource)
      call getarg(2,filetarget)
      call getarg(3,regridfile)
      call getarg(4,cnts)
      call getarg(5,cimt)
      call getarg(6,cjmt)
      call getarg(7,cntt)
      call getarg(8,cformat)
      call getarg(9,cidl)
      read(cnts,'(I4)') nts
      read(cimt,'(I4)') imt
      read(cjmt,'(I4)') jmt
      read(cntt,'(I4)') ntt

      write(*,*) "netcdf format=",cformat

c*** Parse format
      i=1
      offset(1)=1
      cformat=adjustl(trim(cformat))  ! remove leading and trailing blank spaces
      ccopy=cformat

c    search for commas and separate different dimensions
      write(*,*) "Parsing dimensions"      
      do while (index(ccopy,',') .ne. 0) 
         i=i+1
         offset(i)=index(ccopy,',')
         cdim(i-1)=ccopy(:offset(i)-1)
         ccopy=ccopy(offset(i)+1:)
         write(*,*) "dim(",i-1,")=",cdim(i-1)
      enddo
      cdim(i)=ccopy
      write(*,*) "dim(",i,")=",cdim(i)
      imax=i

      status = nf_open(trim(filesource),nf_nowrite,fid)
      if (status .ne. nf_noerr) write(*,*) 
     *     "UNABLE TO OPEN SOURCE FILE ",trim(filesource)

      status = nf_create(trim(filetarget),nf_clobber,fidt)
      if(status .ne. nf_noerr) write(*,*)
     &     "ERROR CREATING ",trim(filetarget)

c     read dimensions from input netcdf file     
      do i=1,imax
         status = nf_inq_dimid(fid,cdim(i),vid)
         status = nf_inq_dimlen(fid,vid,dim(i))
      enddo
     
c     find if format is 'ij', 'ijl', 'kij'
      if (imax .eq. 3) then
         if (cdim(imax) .eq. 'lon' .and. cdim(imax-1) .eq. 'lat') then 
            cform='ijl'
            ims=dim(3)
            jms=dim(2)
            lms=dim(1)
         else
            cform='kij'
            kms=dim(3)
            ims=dim(2)
            jms=dim(1)
         endif
      elseif (imax .eq. 2) then
         cform='ij'
         ims=dim(2)
         jms=dim(1)
      endif

      write(*,*) "array format=",cform
      if (cform .eq. 'ijl') write(*,*) "i j l=",ims,jms,lms
      if (cform .eq. 'ij') write(*,*) "i j=",ims,jms

c***  Initialize exchange grid
      call init_regrid(regridfile,x2grids,ims,jms,nts,imt,jmt,ntt)

c***  first loop to find ids of the different dimensions
      idims=0;idjms=0;idlms=0
      status = nf_inq(fid, ndims, nvars, ngatts, UNLIMDIMID)
      do i = 1, nvars
         status = nf_inq_var(fid, i, cval, itype, ndim, ishape, natt)
         if (cform .eq. 'ijl') then
            if (cval .eq. 'lon')   idims=ishape(1)
            if (cval .eq. 'lat')   idjms=ishape(1)
            if (cval .eq. cdim(1)) idlms=ishape(1)
         endif
         if (cform .eq. 'ij') then
            if (cval .eq. 'lon')   idims=ishape(1)
            if (cval .eq. 'lat')   idjms=ishape(1)
         endif
      enddo
      write(*,*) "idims idjms idlms",idims,idjms,idlms

c***  2nd loop to define variables which have format == cformat 

      do i = 1, nvars
         status = nf_inq_var(fid, i, cval, itype, ndim, ishape, natt)
         
         if (cform .eq. 'ijl') then
c     extract values of variables having correct format 
            if (ishape(1) .eq. idims .and. 
     &          ishape(2) .eq. idjms .and. 
     &          ishape(3) .eq. idlms .and. ndim .eq. 3) then
               
               allocate(tijl(imt,jmt,lms,ntt))
               status = nf_inq_varid(fid,cval,vid)
               
c     define each remaped variable to netcdf output file
               vname=trim(cval)//'(im,jm,lm,tile)'
               call defvar(fidt,imt,jmt,ntt,tijl,trim(vname))
               deallocate(tijl)
            endif
         endif

         if (cform .eq. 'ij') then
c     extract values of variables having correct format 
            if (ishape(1) .eq. idims .and. 
     &          ishape(2) .eq. idjms .and. ndim .eq. 2) then

               allocate(tij(imt,jmt,ntt))
               status = nf_inq_varid(fid,cval,vid)
               
c     define each remaped variable to netcdf output file
               vname=trim(cval)//'(im,jm,tile)'
               call defvar(fidt,imt,jmt,ntt,tij,trim(vname))
               deallocate(tij)
            endif
         endif       
      enddo      !end of first loop

      status = nf_close(fid)
      status = nf_enddef(fidt)

      if (status .ne. NF_NOERR) 
     &     write(*,*) "Problem with enddef"

c***  3rd loop to write variables which have format == cformat 
      status = nf_open(trim(filesource),nf_nowrite,fid)
      status = nf_inq(fid, ndims, nvars, ngatts, UNLIMDIMID)

      do i = 1, nvars
         status = nf_inq_var(fid, i, cval, itype, ndim, ishape, natt)
         if (cform .eq. 'ijl') then
c     find ids of the different dimensions
         if (cval .eq. 'lon')   idims=ishape(1)
         if (cval .eq. 'lat')   idjms=ishape(1)
         if (cval .eq. cdim(1)) idlms=ishape(1)

c     extract values of variables having correct format 
            if (ishape(1) .eq. idims .and. 
     &          ishape(2) .eq. idjms .and. 
     &          ishape(3) .eq. idlms .and. ndim .eq. 3) then

               allocate(sijl(ims,jms,lms,nts),sijl4(ims,jms,lms,nts),
     &              sijl4shift(ims,jms,lms,nts),
     &              tijl(imt,jmt,lms,ntt),tijl4(imt,jmt,lms,ntt))
               
               status = nf_inq_varid(fid,cval,vid)
               status = nf_get_var_real(fid,vid,sijl4)
               
c     shift by 180 degrees if data aligned on Greenwich meridian
               if (trim(cidl) .ne. 'y' .and. trim(cidl) .ne. 'Y') then
                  is=ims/2
                  do n=1,is
                     sijl4shift(n,:,:,:)=sijl4(n+is,:,:,:)
                     sijl4shift(n+is,:,:,:)=sijl4(n,:,:,:)
                  enddo
                  sijl4=sijl4shift
               endif
c     remap variables with correct format
               sijl=sijl4

               do n=1,lms
                  call do_regrid(x2grids,sijl(:,:,n,:),tijl(:,:,n,:))
               enddo
               tijl4=tijl

c     write each remaped variable to netcdf output file  
               write(*,*) "write ",adjustl(trim(cval))," in output file"
               call write_data(fidt,trim(cval),tijl4)

               deallocate(sijl4, sijl4shift, sijl) 
               deallocate(tijl4, tijl)
            endif
         endif

         if (cform .eq. 'ij') then
c     find ids of the different dimensions
            if (cval .eq. 'lon')   idims=ishape(1)
            if (cval .eq. 'lat')   idjms=ishape(1)

c     extract values of variables having correct format 
            if (ishape(1) .eq. idjms .and. 
     &           ishape(2) .eq. idims .and. ndim .eq. 2) then

               allocate(sij(ims,jms,nts),sij4(ims,jms,nts),
     &              sij4shift(ims,jms,nts),
     &          tij(imt,jmt,ntt),tij4(imt,jmt,ntt))    
               
               status = nf_inq_varid(fid,cval,vid)
               status = nf_get_var_real(fid,vid,sij4)

c     shift by 180 degrees if data aligned on Greenwich meridian
               if (trim(cidl) .ne. 'y' .and. trim(cidl) .ne. 'Y') then
                  is=ims/2
                  do n=1,is
                     sij4shift(n,:,:)=sij4(n+is,:,:)
                     sij4shift(n+is,:,:)=sij4(n,:,:)
                  enddo
                  sij4=sij4shift
               endif

c     remap variables with correct format
               sij=sij4
               call do_regrid(x2grids,sij(:,:,:),tij(:,:,:))
               tij4=tij

c     write each remaped variable to netcdf output file  
               write(*,*) "write ",adjustl(trim(cval))," in output file"
               call write_data(fidt,trim(cval),tij4)

               deallocate(sij4, sij,sij4shift) 
               deallocate(tij4, tij)
            endif
         endif       

         if (cformat .eq. 'kij') then
c     find ids of the different dimensions
            if (cval .eq. 'lon')   idims=ishape(1)
            if (cval .eq. 'lat')   idjms=ishape(1)
            if (cval .eq. cdim(1)) idlms=ishape(1)

c     extract values of variables having correct format 
            if (ishape(1) .eq. idlms .and. 
     &           ishape(2) .eq. idims .and. 
     &           ishape(3) .eq. idjms .and. ndim .eq. 3) then
               
               allocate(slij(lms,ims,jms,nts),slij4(lms,ims,jms,nts))
               allocate(tlij(lms,imt,jmt,ntt),tlij4(lms,imt,jmt,ntt))    
               
               status = nf_inq_varid(fid,cval,vid)
               status = nf_get_var_real(fid,vid,slij4)
               
c     remap variables with correct format
               slij=slij4
               
               do n=1,lms
                  call do_regrid(x2grids,slij(n,:,:,:),tlij(n,:,:,:))
               enddo
               tlij4=tlij
               
c     write each remaped variable to netcdf output file  
               write(*,*) "write ",adjustl(trim(cval))," in output file"
               call write_data(fidt,trim(cval),tlij4)
               
               deallocate(slij4, slij) 
               deallocate(tlij4, tlij)
               
            endif
         endif
      enddo    !end of 2nd loop

      status = nf_close(fidt)
      status = nf_close(fid)
      
      write(6,*) "wrote:",trim(filetarget),">>>>>>>>>>>>"
      write(6,*) ""
      
      end program ncll2cs
