      subroutine init_xgrid_zonal(x2grids,imsource,jmsource,
     &     ntilessource,imtarget,jmtarget,ntilestarget)

!@sum Initialize exchange grid for constant latitude band zonal means (no zigzag)
!@auth Denis Gueyffier
      use regrid_com
      implicit none
      include 'netcdf.inc'
      include 'mpif.h'

      type (x_2grids), intent(inout) :: x2grids
      integer, intent(in) :: imsource,jmsource,imtarget,jmtarget
      integer, intent(in) :: ntilessource,ntilestarget
      real*8,  allocatable, dimension(:) :: xgrid_area  !local variable
      integer, allocatable, dimension(:) :: tile        !local variable
      integer, allocatable, dimension(:,:) :: ijcub     !local variable
      integer, allocatable, dimension(:,:) :: ijlatlon  !local variable
      integer :: ncells                                 !local variable
      integer :: maxkey                                 !local variable

      integer :: status,fid,n,vid,ikey,jlat,nc,nc2,
     *     ierr,icc,jcc,jlatm,ic,jc
      integer :: itile,j,idomain,iic,jjc,index,indexc
      character*200 :: exchfile
      real*8 :: checkarea
      character(len=10) :: imch,jmch,icch,jcch 


      x2grids%imsource=imsource
      x2grids%jmsource=jmsource
      x2grids%ntilessource=ntilessource
      x2grids%imtarget=imtarget
      x2grids%jmtarget=jmtarget
      x2grids%ntilestarget=ntilestarget

      write(imch,'(i10)') imsource
      write(jmch,'(i10)') jmsource
      write(icch,'(i10)') imtarget
      write(jcch,'(i10)') jmtarget
      imch=trim(adjustl(imch))
      jmch=trim(adjustl(jmch))
      icch=trim(adjustl(icch))
      jcch=trim(adjustl(jcch))
      exchfile="remap"//trim(imch)//"-"//trim(jmch)
     *     //"C"//trim(icch)//"-"//trim(jcch)//".nc"
      write(*,*) "filename=",exchfile
      imch=trim(adjustl(imch))
      jmch=trim(adjustl(jmch))
      icch=trim(adjustl(icch))
      jcch=trim(adjustl(jcch))
      
      exchfile="remap"//trim(imch)//"-"//trim(jmch)
     *     //"C"//trim(icch)//"-"//trim(jcch)//".nc"
      write(*,*) "-->filename=",exchfile

c      exchfile="exexch.nc"

c      
c Read weights
c
      status = nf_open(trim(exchfile),nf_nowrite,fid)
      if (status .ne. NF_NOERR) write(*,*) "UNABLE TO OPEN REMAP FILE"
 
      status = nf_inq_dimid(fid,'ncells',vid)
      status = nf_inq_dimlen(fid,vid,ncells)

c
c Allocate arrays with size depending on ncells
c
      allocate(xgrid_area(ncells))
      allocate(ijcub(2,ncells))
      allocate(ijlatlon(2,ncells))
      allocate(tile(ncells))
      
      if (AM_I_ROOT()) then   
         status = nf_inq_varid(fid,'xgrid_area',vid)
         status = nf_get_var_double(fid,vid,xgrid_area)

         
         status = nf_inq_varid(fid,'tile1',vid)
         status = nf_get_var_int(fid,vid,tile)
         
         status = nf_inq_varid(fid,'tile1_cell',vid)
         status = nf_get_var_int(fid,vid,ijcub)
         
         status = nf_inq_varid(fid,'tile2_cell',vid)
         status = nf_get_var_int(fid,vid,ijlatlon)
         
         status = nf_close(fid)
      endif

c
c     Broadcast x_grid area and indices to all procs
c
      nc2=2*ncells
     
      write(6,*) "BEFMPI ROLLED"

      call MPI_BCAST( xgrid_area, ncells, MPI_DOUBLE_PRECISION,
     *     0, MPI_COMM_WORLD, ierr ) 

      call MPI_BCAST( ijcub, nc2, MPI_INTEGER, 0, 
     *     MPI_COMM_WORLD, ierr ) 

      call MPI_BCAST( ijlatlon, nc2, MPI_INTEGER, 0, 
     *     MPI_COMM_WORLD, ierr ) 

      call MPI_BCAST( tile, ncells, MPI_INTEGER, 0, 
     *     MPI_COMM_WORLD, ierr ) 

      write(6,*) "AFT BCAST ROLLED"
c****
c*     ikey=ikey+1
c*     index=(iic-1)*jcmax+jjc
c*     array index_key : ikey->index
c*     array jlat_key : ikey->jlat
c*     array xarea_key : ikey->xgrid_area
c****

      
c     calculate maxkey = size of azonal 
 
      ikey=1
      jlatm=x2grids%jmtarget
      do jlat=1,jlatm
         do n=1,ncells
            itile=tile(n)
            if (itile .eq. mytile) then
               j=ijlatlon(2,n)
               if (j .eq. jlat) then
                  ikey=ikey+1
               endif
            endif
         enddo
      enddo
      
      maxkey=ikey-1
      x2grids%xgrid%maxkey=maxkey
      x2grids%xgrid%ncells=ncells
    
      allocate(x2grids%xgrid%index_key(maxkey))
      allocate(x2grids%xgrid%jlat_key(maxkey))
      allocate(x2grids%xgrid%xarea_key(maxkey))

      write(6,*) "AFT ALOC ROLLED, maxkey",maxkey

      x2grids%xgrid%index_key(:)=0
      x2grids%xgrid%jlat_key(:)=0
      x2grids%xgrid%xarea_key(:)=0.d0

      checkarea=0.d0
      ikey=1
      ic=x2grids%imsource
      jc=x2grids%jmsource
      do iic=1,ic
         do jjc=1,jc
            index=(iic-1)*jc+jjc
            do n=1,ncells
               itile=tile(n)
               if (itile .eq. mytile) then
                  icc=ijcub(1,n)
                  jcc=ijcub(2,n)
                  indexc=(icc-1)*jc+jcc
c               write(6,*) "index indexc",index,indexc
                  if ( index .eq. indexc ) then
                     jlat=ijlatlon(2,n)
                     x2grids%xgrid%index_key(ikey)=index
                     x2grids%xgrid%jlat_key(ikey)=jlat
                     x2grids%xgrid%xarea_key(ikey)=xgrid_area(n)    
c                     write(6,*) "ikey=",ikey
                     ikey=ikey+1
                     checkarea=checkarea+xgrid_area(n)
                  endif
               endif
            enddo
         enddo
      enddo
      write(*,*) "ikey=",ikey-1
      write(*,*) "->checkarea=",checkarea

      deallocate(xgrid_area)
      deallocate(ijcub)
      deallocate(ijlatlon)
      deallocate(tile)
      
      end subroutine init_xgrid_zonal
c****


      subroutine zonalmean_band(x2grids,value,iinput,jinput,
     *     zonal_mean,area_latband)

!@sum  Calculate zonal mean using latitute band algorithm (not zigzags on budget grid)
!@auth Denis Gueyffier

      use regrid_com
      implicit none
      include 'netcdf.inc'
      type (x_2grids), intent(in) :: x2grids
      integer, intent(in) :: iinput,jinput
      real*8 :: value,zonal_mean(x2grids%jmtarget),
     &     area_latband(x2grids%jmtarget),
     &     zstore,astore
      integer :: ikey,jlat,current_index,index

      if ( ( (iinput .le. ie) .and. (iinput .ge. is)) .and.
     &     ( (jinput .le. je) .and. (jinput .ge. js)) ) then

      current_index=x2grids%jmsource*(iinput-1)+jinput

      do ikey=1,x2grids%xgrid%maxkey
         index=x2grids%xgrid%index_key(ikey)            
         if (index .eq. current_index) then
            jlat=x2grids%xgrid%jlat_key(ikey) 
            zstore=zonal_mean(jlat)
            zonal_mean(jlat)=zstore
     &           +x2grids%xgrid%xarea_key(ikey)*value
            astore=area_latband(jlat)
            area_latband(jlat)=astore
     &           +x2grids%xgrid%xarea_key(ikey)
         endif
      enddo

      endif

      end subroutine zonalmean_band
c*


      subroutine parallel_regrid(x2grids,tsource,   
     &     ttarget,atarget)

!@sum  Parallel regriding from source to target grid
!@+    x_2grids is used to determine type of source and target 
!@+    grids (cubed-sphere or lat-lon)
!@auth Denis Gueyffier
      use regrid_com
      implicit none
      include 'netcdf.inc'
      type (x_2grids), intent(in) :: x2grids
      real*8 :: tsource(isd:ied,jsd:jed)
      real*8 :: ttarget(x2grids%imtarget,x2grids%jmtarget)
     &     ,atarget(x2grids%imtarget,x2grids%jmtarget)
      integer :: n,icub,jcub,i,j,itile,ilon,jlat,ikey
      character*120:: ofi
      integer ::  status, fid, vid
       
      if ( (x2grids%ntilessource .eq. 6) .and.   !cs2ll
     &     (x2grids%ntilestarget .eq. 1) ) then   
         
         atarget(:,:) = 0.d0
         ttarget(:,:) = 0.d0

         do ikey=1,x2grids%xgrid%maxkey
            icub=x2grids%xgrid%icub_key(ikey)
            jcub=x2grids%xgrid%jcub_key(ikey)
            ilon=x2grids%xgrid%ilon_key(ikey)
            jlat=x2grids%xgrid%jlat_key(ikey)
            itile=x2grids%xgrid%itile_key(ikey)
            
            atarget(ilon,jlat) = atarget(ilon,jlat) + 
     &           x2grids%xgrid%xarea_key(ikey)
            ttarget(ilon,jlat) = ttarget(ilon,jlat) + 
     &           x2grids%xgrid%xarea_key(ikey)
     &           * tsource(icub,jcub) 
         enddo
         
c
c     sum all contributions
c     
         call SUMXPE(ttarget)
         call SUMXPE(atarget)
         
c     
c     root proc section
c     
         if (AM_I_ROOT()) then   
            ttarget(:,:) = ttarget(:,:)/atarget(:,:)

            ofi='tstout.nc'
            status = nf_open(trim(ofi),nf_write,fid)
            if (status .ne. NF_NOERR) write(*,*) NF_STRERROR(status)
            status = nf_inq_varid(fid,'lwup_sfc',vid)
            write(*,*) NF_STRERROR(status)
            status = nf_put_var_double(fid,vid,ttarget)
            write(*,*) "STATUS",NF_STRERROR(status),"<<"
            status = nf_close(fid)
         endif
         
      endif  !cs2ll

      end subroutine parallel_regrid
c*


      subroutine regrid_exact(x2grids,tsource,ttarget,atarget)
 
!@sum  Parallel regriding using double-double arithmetics for reproducible results
!@+    independent of #procs
!@auth Denis Gueyffier
      use regrid_com
      implicit none
      include 'netcdf.inc'
      type (x_2grids), intent(in) :: x2grids
      real*8 :: tsource(isd:ied,jsd:jed)
      real*8 :: ttarget(x2grids%imtarget,x2grids%jmtarget),
     &     atarget(x2grids%imtarget,x2grids%jmtarget)
      complex*16 :: ttarget_dd(x2grids%imtarget,x2grids%jmtarget),
     &     atarget_dd(x2grids%imtarget,x2grids%jmtarget)
      complex*16 :: ttarget_tmp(1),tmp_dd,atarget_tmp(1),xarea_dd(1),
     &     tmp_dd_vec(1)
      integer :: n,icub,jcub,i,j,itile,ilon,jlat,ikey,imlon,jmlat
      character*120:: ofi
      integer ::  status, fid, vid
     

      if ( (x2grids%ntilessource .eq. 6) .and. !cs2ll
     &     (x2grids%ntilestarget .eq. 1) ) then   

         imlon=x2grids%imtarget 
         jmlat=x2grids%jmtarget
         
         do j=1,jmlat
            do i=1,imlon
               atarget_dd(i,j) = dcmplx(0.d0,0.d0)
               ttarget_dd(i,j) = dcmplx(0.d0,0.d0)
            enddo
         enddo

         write(*,*) "maxkey==",x2grids%xgrid%maxkey

         do ikey=1,x2grids%xgrid%maxkey

            icub=x2grids%xgrid%icub_key(ikey)
            jcub=x2grids%xgrid%jcub_key(ikey)
            ilon=x2grids%xgrid%ilon_key(ikey)
            jlat=x2grids%xgrid%jlat_key(ikey)
            itile=x2grids%xgrid%itile_key(ikey)
          
            atarget_tmp(1)=atarget_dd(ilon,jlat)
            xarea_dd(1)=dcmplx(x2grids%xgrid%xarea_key(ikey)
     &           , 0.d0)
            call add_dd(xarea_dd,atarget_tmp,1) 
            atarget_dd(ilon,jlat)=atarget_tmp(1) !OK

            ttarget_tmp(1)=ttarget_dd(ilon,jlat)
            call dxd(x2grids%xgrid%xarea_key(ikey),
     &           tsource(icub,jcub),tmp_dd)
            tmp_dd_vec(1)=tmp_dd
            call add_dd(tmp_dd_vec,ttarget_tmp,1)   
            ttarget_dd(ilon,jlat)=ttarget_tmp(1)

         enddo
      
         call sumxpe2d_exact(ttarget_dd,imlon,jmlat)
         call sumxpe2d_exact(atarget_dd,imlon,jmlat)
         
         write(*,*) "MYTILE",mytile
c     
c     root proc section
c     

         if (AM_I_ROOT()) then   
            do j=1,jmlat
               do i=1,imlon
                  ttarget(i,j)=real(ttarget_dd(i,j))
                  atarget(i,j)=real(atarget_dd(i,j))
               enddo
            enddo
            
            ttarget(:,:) = ttarget(:,:)/atarget(:,:)
            
            ofi='tstout.nc'
            status = nf_open(trim(ofi),nf_write,fid)
            if (status .ne. NF_NOERR) write(*,*) NF_STRERROR(status)
            status = nf_inq_varid(fid,'lwup_sfc',vid)
            write(*,*) NF_STRERROR(status)
            status = nf_put_var_double(fid,vid,ttarget)
            write(*,*) "STATUS",NF_STRERROR(status),"<<"
            status = nf_close(fid)
         endif
      endif
      
      end subroutine regrid_exact
c*


      subroutine root_regrid(x2gridsroot,tsource,ttarget)

!@sum  Root processor regrids data from source grid to target grid
!@auth Denis Gueyffier
      use regrid_com
      implicit none
      type (x_2gridsroot), intent(in) :: x2gridsroot
      real*8 :: tsource(x2gridsroot%imsource,x2gridsroot%jmsource,
     &     x2gridsroot%ntilessource)
      real*8 :: ttarget(x2gridsroot%imtarget,x2gridsroot%jmtarget,
     &     x2gridsroot%ntilestarget)
     &     ,atarget(x2gridsroot%imtarget,x2gridsroot%jmtarget,
     &     x2gridsroot%ntilestarget)
      integer :: n,icub,jcub,i,j,itile,icc,jcc,il,jl


      if ((x2gridsroot%ntilessource .eq. 1) .and. 
     &     (x2gridsroot%ntilestarget .eq. 6)) then
         
c     ll2cs
         
         if (AM_I_ROOT()) then   
            write(6,*) "ROOT REGRID LL2CS"
            atarget(:,:,:) = 0.d0
            ttarget(:,:,:) = 0.d0
            
            do n=1,x2gridsroot%xgridroot%ncells
               itile=x2gridsroot%xgridroot%tile(n)
               icc=x2gridsroot%xgridroot%ijcub(1,n)
               jcc=x2gridsroot%xgridroot%ijcub(2,n)
               il=x2gridsroot%xgridroot%ijlatlon(1,n)
               jl=x2gridsroot%xgridroot%ijlatlon(2,n)
               
               atarget(icc,jcc,itile) = atarget(icc,jcc,itile) 
     &              + x2gridsroot%xgridroot%xgrid_area(n)
               ttarget(icc,jcc,itile) = ttarget(icc,jcc,itile) 
     &              + x2gridsroot%xgridroot%xgrid_area(n)
     &              *tsource(il,jl,1)
            enddo
            
            do itile=1,x2gridsroot%ntilestarget
               do j=1,x2gridsroot%jmtarget
                  do i=1,x2gridsroot%imtarget
                     ttarget(i,j,itile) = ttarget(i,j,itile)
     &                    /atarget(i,j,itile)
                  enddo
               enddo
            enddo
            write(6,*) "END ROOT REGRID LL2CS"
         endif
         
      endif
      
      end subroutine root_regrid
c*

         
      subroutine init_regrid(x2grids,imsource,jmsource,
     &     ntilessource,imtarget,jmtarget,ntilestarget)

!@sum  Reads regriding file on root proc, broadcasts the 
!@+    x_grid data to all processes then instanciates locally
!@+    the x_2grids derived type (x_2grids type=x_grid plus info about
!@+    source and target grids). It also initializes domain decomposition
!@+    variables through dist_grid derived type. Will soon use dd2d 
!@+    derived type in place of dist_grid.
!@auth Denis Gueyffier
      use regrid_com
      implicit none
      include 'netcdf.inc'
      include 'mpif.h'
c      type (dist_grid), intent(in) :: grid
      type (x_2grids), intent(inout) :: x2grids
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
      integer :: ierr,i,icub,jcub,maxkey
      character*200 :: exchfile
      character(len=10) :: imch,jmch,icch,jcch


c modify lines below when we begin using derived type dd2d 
#ifndef CUBE_GRID  
c     set variables ("Constructor")
      is=grid%I_STRT
      ie=grid%I_STOP
      isd=grid%I_STRT_HALO
      ied=grid%I_STOP_HALO

      js=grid%J_STRT
      je=grid%J_STOP
      jsd=grid%J_STRT_HALO
      jed=grid%J_STOP_HALO

      call MPI_COMM_RANK(MPI_COMM_WORLD,gid,ierr) ! will soon be replaced by gid = grid%dd2d%gid
c      call MPI_COMM_SIZE(MPI_COMM_WORLD,totPEs,ierr)

      mytile=(gid/dom_per_tile)+1   ! will soon be mytile=grid%dd2d%mytile
#endif

      x2grids%imsource=imsource
      x2grids%jmsource=jmsource
      x2grids%ntilessource=ntilessource
      x2grids%imtarget=imtarget
      x2grids%jmtarget=jmtarget
      x2grids%ntilestarget=ntilestarget

c
      if (AM_I_ROOT()) then   

         write(imch,'(i10)') imsource
         write(jmch,'(i10)') jmsource
         write(icch,'(i10)') imtarget
         write(jcch,'(i10)') jmtarget
         imch=trim(adjustl(imch))
         jmch=trim(adjustl(jmch))
         icch=trim(adjustl(icch))
         jcch=trim(adjustl(jcch))
         exchfile="remap"//trim(imch)//"-"//trim(jmch)
     *        //"C"//trim(icch)//"-"//trim(jcch)//".nc"
         write(*,*) "filename=",exchfile
c     
c     Read weights
c     
         status = nf_open(trim(exchfile),nf_nowrite,fid)
         if (status .ne. NF_NOERR) write(*,*) 
     *        "UNABLE TO OPEN REMAP FILE"
         
         status = nf_inq_dimid(fid,'ncells',vid)
         status = nf_inq_dimlen(fid,vid,ncells)
      endif
            
c     Broadcast value of ncells & Allocate arrays with size 
c     depending on ncells on each processor

      call MPI_BCAST( ncells, 1, MPI_INTEGER, 0, 
     *     MPI_COMM_WORLD, ierr ) 

      allocate(xgrid_area(ncells))
      allocate(ijcub(2,ncells))
      allocate(ijlatlon(2,ncells))
      allocate(tile(ncells))


      if (AM_I_ROOT()) then   
         status = nf_inq_varid(fid,'xgrid_area',vid)
         status = nf_get_var_double(fid,vid,xgrid_area)
         status = nf_inq_varid(fid,'tile1',vid)
         status = nf_get_var_int(fid,vid,tile)
         status = nf_inq_varid(fid,'tile1_cell',vid)
         status = nf_get_var_int(fid,vid,ijcub)
         status = nf_inq_varid(fid,'tile2_cell',vid)
         status = nf_get_var_int(fid,vid,ijlatlon)
         status = nf_close(fid)
      endif
      
c
c     Broadcast x_grid area and indices to all procs
c
      nc2=2*ncells
     
      call MPI_BCAST( xgrid_area, ncells, MPI_DOUBLE_PRECISION,
     *     0, MPI_COMM_WORLD, ierr ) 
      call MPI_BCAST( ijcub, nc2, MPI_INTEGER, 0, 
     *     MPI_COMM_WORLD, ierr ) 
      call MPI_BCAST( ijlatlon, nc2, MPI_INTEGER, 0, 
     *     MPI_COMM_WORLD, ierr ) 
      call MPI_BCAST( tile, ncells, MPI_INTEGER, 0, 
     *     MPI_COMM_WORLD, ierr ) 

      
      
      if ((ntilessource .eq. 6) .and. (ntilestarget .eq. 1)) then   !cs2ll
         
c     fill x2grids with values local to the domain
c     first calculate maxkey
         
         ikey=1
         do n=1,ncells
            icub=ijcub(1,n)
            jcub=ijcub(2,n)
            itile=tile(n)
            if (itile .eq. mytile) then
               if ( ( (icub .le. ie) .and. (icub .ge. is)) .and.
     &              ( (jcub .le. je) .and. (jcub .ge. js)) ) then
                  ikey=ikey+1
               endif
            endif
         enddo
         
         maxkey=ikey-1
         
         write(*,*) "maxkey=",maxkey

         allocate(x2grids%xgrid%icub_key(maxkey),
     &        x2grids%xgrid%jcub_key(maxkey),
     &        x2grids%xgrid%ilon_key(maxkey),
     &        x2grids%xgrid%jlat_key(maxkey),
     &        x2grids%xgrid%xarea_key(maxkey),
     &        x2grids%xgrid%itile_key(maxkey)
     &        ) 
         
         ikey=1
         
         do n=1,ncells
            icub=ijcub(1,n)
            jcub=ijcub(2,n)
            itile=tile(n)
            
            if (itile .eq. mytile) then
               if ( ( (icub .le. ie) .and. (icub .ge. is)) .and.
     &              ( (jcub .le. je) .and. (jcub .ge. js)) ) then
                  x2grids%xgrid%icub_key(ikey)=icub
                  x2grids%xgrid%jcub_key(ikey)=jcub
                  x2grids%xgrid%ilon_key(ikey)=ijlatlon(1,n)
                  x2grids%xgrid%jlat_key(ikey)=ijlatlon(2,n)
                  x2grids%xgrid%xarea_key(ikey)=xgrid_area(n)
                  x2grids%xgrid%itile_key(ikey)=tile(n)
                  ikey=ikey+1
               endif
            endif
         enddo

         x2grids%xgrid%maxkey=maxkey
         
      else if ((ntilessource .eq. 1) .and. (ntilestarget .eq. 6)) then !ll2cs
         ikey=1

         
      endif    
         deallocate(xgrid_area)
         deallocate(ijcub)
         deallocate(ijlatlon)
         deallocate(tile)
         
     
      
      end subroutine init_regrid
c*


      subroutine init_regrid_root(x2gridsroot,imsource,jmsource,
     &     ntilessource,imtarget,jmtarget,ntilestarget)

!@sum  Reads regriding file on root proc, broadcasts the 
!@+    x_grid data to all processes then instanciates locally
!@+    the x_2grids derived type (x_2grids type=x_grid plus info about
!@+    source and target grids). It also initializes domain decomposition
!@+    variables through dist_grid derived type. Will soon use dd2d 
!@+    derived type in place of dist_grid.
!@auth Denis Gueyffier
      use regrid_com
      implicit none
      include 'netcdf.inc'
      include 'mpif.h'
      type (x_2gridsroot), intent(inout) :: x2gridsroot
      type (x_gridroot) :: xgridroot
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
      character*200 :: exchfile
      character(len=10) :: imch,jmch,icch,jcch


c modify lines below when we begin using derived type dd2d 
#ifndef CUBE_GRID  
c     set variables ("Constructor")
      is=grid%I_STRT
      ie=grid%I_STOP
      isd=grid%I_STRT_HALO
      ied=grid%I_STOP_HALO

      js=grid%J_STRT
      je=grid%J_STOP
      jsd=grid%J_STRT_HALO
      jed=grid%J_STOP_HALO

      call MPI_COMM_RANK(MPI_COMM_WORLD,gid,ierr) ! will soon be replaced by gid = grid%dd2d%gid
c      call MPI_COMM_SIZE(MPI_COMM_WORLD,totPEs,ierr)

      mytile=(gid/dom_per_tile)+1   ! will soon be mytile=grid%dd2d%mytile
#endif

      x2gridsroot%imsource=imsource
      x2gridsroot%jmsource=jmsource
      x2gridsroot%ntilessource=ntilessource
      x2gridsroot%imtarget=imtarget
      x2gridsroot%jmtarget=jmtarget
      x2gridsroot%ntilestarget=ntilestarget


c
      if (AM_I_ROOT()) then   

         write(imch,'(i10)') imsource
         write(jmch,'(i10)') jmsource
         write(icch,'(i10)') imtarget
         write(jcch,'(i10)') jmtarget
         imch=trim(adjustl(imch))
         jmch=trim(adjustl(jmch))
         icch=trim(adjustl(icch))
         jcch=trim(adjustl(jcch))
         exchfile="remap"//trim(imch)//"-"//trim(jmch)
     *        //"C"//trim(icch)//"-"//trim(jcch)//".nc"
         write(*,*) "filename=",exchfile
c     
c     Read weights
c     
         status = nf_open(trim(exchfile),nf_nowrite,fid)
         if (status .ne. NF_NOERR) write(*,*) 
     *        "UNABLE TO OPEN REMAP FILE"
         
         status = nf_inq_dimid(fid,'ncells',vid)
         status = nf_inq_dimlen(fid,vid,ncells)
      endif
            
c     Broadcast value of ncells & Allocate arrays with size 
c     depending on ncells on each processor

      call MPI_BCAST( ncells, 1, MPI_INTEGER, 0, 
     *     MPI_COMM_WORLD, ierr ) 

      allocate(xgrid_area(ncells))
      allocate(ijcub(2,ncells))
      allocate(ijlatlon(2,ncells))
      allocate(tile(ncells))


      if (AM_I_ROOT()) then   
         status = nf_inq_varid(fid,'xgrid_area',vid)
         status = nf_get_var_double(fid,vid,xgrid_area)
         status = nf_inq_varid(fid,'tile1',vid)
         status = nf_get_var_int(fid,vid,tile)
         status = nf_inq_varid(fid,'tile1_cell',vid)
         status = nf_get_var_int(fid,vid,ijcub)
         status = nf_inq_varid(fid,'tile2_cell',vid)
         status = nf_get_var_int(fid,vid,ijlatlon)
         status = nf_close(fid)
      endif
      
c
c     Broadcast x_grid area and indices to all procs
c
      nc2=2*ncells
     
      call MPI_BCAST( xgrid_area, ncells, MPI_DOUBLE_PRECISION,
     *     0, MPI_COMM_WORLD, ierr ) 
      call MPI_BCAST( ijcub, nc2, MPI_INTEGER, 0, 
     *     MPI_COMM_WORLD, ierr ) 
      call MPI_BCAST( ijlatlon, nc2, MPI_INTEGER, 0, 
     *     MPI_COMM_WORLD, ierr ) 
      call MPI_BCAST( tile, ncells, MPI_INTEGER, 0, 
     *     MPI_COMM_WORLD, ierr ) 

      allocate(x2gridsroot%xgridroot%xgrid_area(ncells))
      allocate(x2gridsroot%xgridroot%ijcub(2,ncells))
      allocate(x2gridsroot%xgridroot%ijlatlon(2,ncells))
      allocate(x2gridsroot%xgridroot%tile(ncells))

      x2gridsroot%xgridroot%xgrid_area(:)=xgrid_area(:)
      x2gridsroot%xgridroot%ijcub(:,:)=ijcub(:,:)
      x2gridsroot%xgridroot%ijlatlon(:,:)=ijlatlon(:,:)
      x2gridsroot%xgridroot%tile(:)=tile(:)
      x2gridsroot%xgridroot%ncells=ncells

      deallocate(xgrid_area)
      deallocate(ijcub)
      deallocate(ijlatlon)
      deallocate(tile)

      end subroutine init_regrid_root
c*



      subroutine readt_regrid_parallel(x2gridsroot,iunit,name,nskip,
     &     ttargloc,ipos)

!@sum  Read input data on lat-lon grid, regrid to global cubbed sphere grid
!@+    then scatter to all subdomains
!@auth Denis Gueyffier

      use regrid_com, only :gid,AM_I_ROOT,x_2gridsroot
      use gatscat_mod
      implicit none
      type (x_2gridsroot), intent(in) :: x2gridsroot
      integer, intent(in) :: iunit
      character*16, intent(in) :: name
      integer, intent(in) :: nskip

      real*8, intent(inout) :: ttargloc(is:ie,js:je) ! remapped & scattered data 
      real*4 :: tsource4(x2gridsroot%imsource,x2gridsroot%jmsource,   ! real*4 data read from input file
     &     x2gridsroot%ntilessource) 
      real*8 :: tsource(x2gridsroot%imsource,x2gridsroot%jmsource,  
     &     x2gridsroot%ntilessource) 
      real*8 :: ttarget(x2gridsroot%imtarget,x2gridsroot%jmtarget,  
     &     x2gridsroot%ntilestarget)
      real*4 :: X              
      integer, intent(in) :: ipos
      integer :: n,ierr
      character*80 :: TITLE     


      if (AM_I_ROOT()) then   
         do n=1,ipos-1
            read(UNIT=iunit,IOSTAT=ierr)
         enddo
         read(UNIt=iunit,IOSTAT=ierr) TITLE, (X,n=1,nskip), tsource4
c     convert from real*4 to real*8
         tsource=tsource4
c     Regrid, form global array
         call root_regrid(x2gridsroot,tsource,ttarget)
       endif

c     Scatter data to every processor
      call unpack_data(ttarget,ttargloc)
      write(6,*) "TTARGLOC>>>",ttargloc,"<<<<"
      end subroutine readt_regrid_parallel
c*


      subroutine dread_regrid_parallel(x2gridsroot,iunit,name,ttargloc)
     
!@sum  Read input data on lat-lon grid, regrid to global cubbed sphere grid
!@+    then scatter to all subdomains
!@auth Denis Gueyffier

      use regrid_com, only :gid,AM_I_ROOT,x_2gridsroot
      use gatscat_mod
      implicit none
      type (x_2gridsroot), intent(in) :: x2gridsroot
      integer, intent(in) :: iunit
      character*16, intent(in) :: name

      real*8, intent(inout) :: ttargloc(is:ie,js:je) ! remapped & scattered data 
      real*4 :: tsource4(x2gridsroot%imsource,x2gridsroot%jmsource,   ! real*4 data read from input file
     &     x2gridsroot%ntilessource) 
      real*8 :: tsource(x2gridsroot%imsource,x2gridsroot%jmsource,  
     &     x2gridsroot%ntilessource) 
      real*8 :: ttarget(x2gridsroot%imtarget,x2gridsroot%jmtarget,  
     &     x2gridsroot%ntilestarget)
      real*4 :: X              
      integer :: n,ierr



   
      if (AM_I_ROOT()) then   
         read(UNIt=iunit,IOSTAT=ierr) tsource4
c     convert from real*4 to real*8
         tsource=tsource4
c     Regrid, form global array
         call root_regrid(x2gridsroot,tsource,ttarget)
       endif

c     Scatter data to every processor
      call unpack_data(ttarget,ttargloc)
      write(6,*) "TTARGLOC>>>",ttargloc,"<<<<"

      end subroutine dread_regrid_parallel
c*





      subroutine dxd (da, db, ddc)
c  This subroutine computes ddc = da x db.
      use regrid_com
      implicit none
      real*8 da, db
      real*8 a1, a2, b1, b2, con, split, t1, t2
      complex*16 ddc
      parameter (split = 134217729.d0)

c   Split da and db into two parts with at most 26 bits each, using the 
c   Dekker-Veltkamp method.
      con = da * split
      a1 = con - (con - da)
      a2 = da - a1
      con = db * split
      b1 = con - (con - db)
      b2 = db - b1

c   Compute da * db using Dekker's method.
      t1 = da * db
      t2 = (((a1*b1 - t1) + a1*b2) + a2*b1) + a2*b2

      ddc = dcmplx(t1, t2)

      return
      end subroutine dxd
c*

      subroutine add_DD (dda, ddb, len)
!@sum  Compute dda + ddb using double-double arithmetics
!@auth Denis Gueyffier, adapted from  Yun He and Chris Ding, 
!@auth Journal of Supercomputing, 2001
      use regrid_com
      implicit none
      real*8 e, t1, t2
      integer i, len
      complex*16 :: dda(len), ddb(len)

      do i = 1, len
c   Compute dda + ddb using Knuth's trick.
      t1 = real(dda(i)) + real(ddb(i))
      e = t1 - real(dda(i))
      t2 = ((real(ddb(i)) - e) + (real(dda(i)) - (t1 - e)))
     &     +imag(dda(i)) + imag(ddb(i))

c   The result is t1 + t2, after normalization.
      ddb(i) = cmplx ( t1 + t2, t2 - ((t1 + t2) - t1) )
      enddo

      return
      end subroutine add_DD
c*

      subroutine test_globalsum_exact()
!@sum Testing reproducible globalsum 
!@auth Denis Gueyffier
      use regrid_com
      implicit none
      include 'mpif.h'
      integer i, imt, jmt
      integer :: myPE, totPEs, stat(MPI_STATUS_SIZE), ierr
      integer :: start, end, len, MPI_SUMDD, itype
      external add_DD
      real*8, dimension(:), allocatable:: array, local_array
      complex*16 :: global_sum

      imt = 120
      jmt = 64

      call MPI_COMM_RANK( MPI_COMM_WORLD, myPE, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, totPEs, ierr )

      len = imt*jmt/totPEs + 1
      write(*,*) 'myPE=', myPE, ' len=', len
      allocate (local_array(len))

C  read data for root proc 
      if (myPE .eq. 0) then
         allocate(array(len*totPEs))
         open(10, file ='etaana.dat', status='unknown')
         do i= 1, imt*jmt
            read(10,*) array(i)
         enddo
         close(10)
         do i=imt*jmt+1,len*totPEs
            array(i)=0.0d0
         enddo
      endif 

C  scatter array
      call MPI_SCATTER(array, len, MPI_REAL8,
     &     local_array, len, MPI_REAL8, 0,
     &     MPI_COMM_WORLD, ierr)

      write(*,*) "HERE, myPE=",myPE

      if (myPE .eq. 0) deallocate (array) 

      call globalsum_exact(global_sum,local_array,len)

      deallocate (local_array)

      end subroutine test_globalsum_exact
c*     



      subroutine globalsum_exact(global_sum,local_array,len)
!@sum This code calculates a reproducible globalsum 
!@+   independent on number of processors using double-double precision method
!@+   complex->complex*16, cmplx->dcmplx, 
!@+   MPI_REAL->MPI_REAL8, MPI_COMPLEX->MPI_DOUBLE_COMPLEX
!@auth Denis Gueyffier, based on Yun He and Chris Ding, 
!@+   Journal of Supercomputing, 2001

      use regrid_com
      implicit none
      include 'mpif.h'
      integer i, imt, jmt
      integer :: myPE, stat(MPI_STATUS_SIZE), ierr
      integer ::  len, MPI_SUMDD, itype
      external add_DD
      real*8, dimension(len):: local_array
      complex*16 :: local_sum, global_sum

      write(*,*) "BEGIN gsum"

      call MPI_COMM_RANK( MPI_COMM_WORLD, myPE, ierr )

C  operator MPI_SUMDD is created based on an external function add_dd
      call MPI_OP_CREATE(add_DD, .TRUE., MPI_SUMDD, ierr)

      
C  Each processor calculates the local_sum of its own section first. 
C  Complex number is defined to represent local_sum and local_sum err.
      local_sum = 0.0
      do i = 1, len
         call add_DD(dcmplx(local_array(i), 0.0d0), local_sum,1)
      enddo

      write(*,*) "MyPE bef reduce=",myPE

C  add all local_sums on each PE to PE0 with MPI_SUMDD.
C  global_sum is a complex number, represents final (sum, error).
      call MPI_REDUCE (local_sum, global_sum, 1, MPI_DOUBLE_COMPLEX, 
     &     MPI_SUMDD, 0, MPI_COMM_WORLD, ierr)

      write(*,*) "MyPE aft reduce=",myPE

      if (myPE.eq.0) then
         write(*,*)'quad precision sum, error= ', global_sum
      endif

      end subroutine globalsum_exact
c*



      subroutine sumxpe2d_exact(arr,iarr,jarr)
!@sum This code calculates a reproducible non rank reducing sum 
!@+   independent on number of processors using double-double arithmetics
!@auth Denis Gueyffier
      implicit none
      include 'mpif.h'
      integer :: ierr,arr_size,myPE
      integer :: MPI_SUMDD
      external add_DD
      integer, intent(in) :: iarr,jarr
      complex*16, dimension(iarr,jarr) :: arr
      complex*16, dimension(:),allocatable ::arr_rsh,arr_tmp

      call MPI_COMM_RANK( MPI_COMM_WORLD, myPE, ierr )
c     operator MPI_SUMDD is created based on an external function add_dd
      call MPI_OP_CREATE(add_DD, .TRUE., MPI_SUMDD, ierr)
      
c     reduction 
      arr_size = size(arr)
      write(*,*) "arrsize",arr_size

c      write(*,*) "arr=",arr(:,:)

      allocate(arr_tmp(arr_size))
      allocate(arr_rsh(arr_size))
      arr_rsh=reshape(arr,(/arr_size/))

      call MPI_REDUCE (arr_rsh, arr_tmp, arr_size, MPI_DOUBLE_COMPLEX, 
     &     MPI_SUMDD, 0, MPI_COMM_WORLD, ierr)
      arr=reshape(arr_tmp,shape(arr))
      deallocate(arr_tmp)
      deallocate(arr_rsh) 
    
      end subroutine sumxpe2d_exact
c*

      subroutine init_csgrid_debug()

!@sum  temporary instanciation of CS grid and domain decomp. using direct
!@+    access to MPP. Used only for debugging purpose.
!@auth Denis Gueyffier

      use mpp_mod
      use mpp_domains_mod
      use fv_mp_mod, only : mp_start,mp_stop
     &     ,fv_domain_decomp=>domain_decomp
      use fv_mp_mod, only : is, ie, js, je, isd, ied, jsd, jed
      use fv_mp_mod, only : gid, domain, tile, npes_x, npes_y
      use fv_grid_utils_mod, only : cosa_s,sina_s
     &     ,grid_utils_init
c     &     ,sw_corner,se_corner,ne_corner,nw_corner
      use fv_arrays_mod, only: fv_atmos_type
      use fv_grid_tools_mod,  only: init_grid, cosa, sina, area, area_c,
     &     dx, dy, dxa, dya, dxc, dyc, grid_type, dx_const, dy_const
      use fv_control_mod, only : uniform_ppm,c2l_ord
      use gs_domain_decomp, ng=>halo_width
      use gatscat_mod

      implicit none
      integer :: i,j,npz
      integer :: npes,ndims
      integer :: commID
      integer, dimension(:), allocatable :: pelist

      real*8, dimension(:,:), allocatable :: P

      type(fv_atmos_type) :: atm
      character*80 :: grid_name = 'Gnomonic'
      character*120:: grid_file = 'Inline'
      logical :: non_ortho

c***  Temporarily use instanciation of grid through fv_grid_tools_mod's init_grid()

      call mpp_init(MPP_VERBOSE)
c code copied from fv_init:
      npes = mpp_npes()
      allocate(pelist(npes))
      call mpp_get_current_pelist( pelist, commID=commID )
      call mp_start(commID)
c      write(6,*) 'commID ',commID
      npx = 48; npy = npx ! 1x1 resolution
      ng = 3 ! number of ghost zones required

      call fv_domain_decomp(npx+1,npy+1,6,ng,grid_type)
c      call mp_stop()

      ndims = 2
      npz = 5
      call init_grid(atm,grid_name,grid_file,
     &     npx+1, npy+1, npz, ndims, 6, ng)

      non_ortho=.true.
      call grid_utils_init(Atm, npx+1, npy+1, npz, Atm%grid, Atm%agrid,
     &     area, area_c, cosa, sina, dx, dy, dxa, dya, non_ortho,
     &     uniform_ppm, grid_type, c2l_ord)

      end subroutine init_csgrid_debug
c****



