      subroutine init_xgrid_zonal(x2grids,imsource,jmsource,
     &     ntilessource,imtarget,jmtarget,ntilestarget)

!@sum Initialize exchange grid for constant latitude band zonal means (not zigzags)
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
     *     ierr,icc,jcc
      integer :: itile,j,idomain,iic,jjc,index,indexc
      character*200 :: exchfile
      real*8 :: checkarea
      character(len=10) :: imch,jmch,icch,jcch 

c         write(imch,'(i10)') ilonm
c         write(jmch,'(i10)') jlatm
c         write(icch,'(i10)') ic
c         write(jcch,'(i10)') jc

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
c*     array index_from_key : ikey->index
c*     array jlat_from_key : ikey->jlat
c*     array xarea_from_key : ikey->xgrid_area
c****

      
c     calculate maxkey = size of azonal 
 
      ikey=1
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
    
      allocate(x2grids%xgrid%index_from_key(maxkey))
      allocate(x2grids%xgrid%jlat_from_key(maxkey))
      allocate(x2grids%xgrid%xarea_from_key(maxkey))

      write(6,*) "AFT ALOC ROLLED, maxkey",maxkey

      x2grids%xgrid%index_from_key(:)=0
      x2grids%xgrid%jlat_from_key(:)=0
      x2grids%xgrid%xarea_from_key(:)=0.d0

      checkarea=0.d0
      ikey=1
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
                     x2grids%xgrid%index_from_key(ikey)=index
                     x2grids%xgrid%jlat_from_key(ikey)=jlat
                     x2grids%xgrid%xarea_from_key(ikey)=xgrid_area(n)    
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
      real*8 :: value,zonal_mean(jlatm),area_latband(jlatm),
     & zstore,astore
      integer :: ikey,jlat,current_index,index

      if ( ( (iinput .le. ie) .and. (iinput .ge. is)) .and.
     &     ( (jinput .le. je) .and. (jinput .ge. js)) ) then

      current_index=jc*(iinput-1)+jinput

      do ikey=1,x2grids%xgrid%maxkey
         index=x2grids%xgrid%index_from_key(ikey)            
         if (index .eq. current_index) then
            jlat=x2grids%xgrid%jlat_from_key(ikey) 
            zstore=zonal_mean(jlat)
            zonal_mean(jlat)=zstore
     &           +x2grids%xgrid%xarea_from_key(ikey)*value
            astore=area_latband(jlat)
            area_latband(jlat)=astore
     &           +x2grids%xgrid%xarea_from_key(ikey)
         endif
      enddo

      endif

      end subroutine zonalmean_band
c*


      subroutine parallel_regrid_cs2ll(x2grids,tcub_loc,   !make it generic, remove cs2ll from proc. name
     &     tlatlon,alatlon)

!@sum  Parallel regriding from cubed-sphere to lat-lon grid
!@auth Denis Gueyffier
      use regrid_com
      implicit none
      type (x_2grids), intent(in) :: x2grids
      real*8 :: tlatlon(ilonm,jlatm),alatlon(ilonm,jlatm)
      real*8 :: tcub_loc(isd:ied,jsd:jsd)
      integer :: n,icub,jcub,i,j,itile,nij_latlon
       
      alatlon(:,:) = 0.d0
      tlatlon(:,:) = 0.d0

      do n=1,x2grids%xgrid%ncells
        icub=x2grids%xgrid%ijcub(1,n)
        jcub=x2grids%xgrid%ijcub(2,n)
        itile=x2grids%xgrid%tile(n)
c        write(6,*) 'ie, icub, is=',ie,icub,is
c        write(6,*) 'je, jcub, js=',je,jcub,js

        if (itile .eq. mytile) then
           if ( ( (icub .le. ie) .and. (icub .ge. is)) .and.
     &          ( (jcub .le. je) .and. (jcub .ge. js)) ) then
c     write(6,*) 'xg_area',x2grids%xgrid%xgrid_area(n)
              i=x2grids%xgrid%ijlatlon(1,n)
              j=x2grids%xgrid%ijlatlon(2,n)
              alatlon(i,j) = alatlon(i,j) + 
     &             x2grids%xgrid%xgrid_area(n)
              tlatlon(i,j) = tlatlon(i,j) + 
     &             x2grids%xgrid%xgrid_area(n)*tcub_loc(icub,jcub)
           endif
        endif
      enddo

c
c     sum all contributions
c

      call SUMXPE(tlatlon)
      call SUMXPE(alatlon)

c
c     root proc section
c

      if (AM_I_ROOT()) then   
         tlatlon(:,:) = tlatlon(:,:)/alatlon(:,:)

      do i=1,ilonm
         do j=1,jlatm
            write(*,*) "par_regrid",tlatlon(i,j),"<----"
         enddo
      enddo

      endif

      end subroutine parallel_regrid_cs2ll
c*


      subroutine regrid_cs2ll_exact(x2grids,tcub_loc,tlatlon,alatlon)
 
!@sum  Parallel regriding from cubed-sphere to lat-lon grid
!@+    version using double-double arithmetics for reproducible results
!@+    independent of #procs
!@auth Denis Gueyffier
      use regrid_com
      implicit none
      include 'netcdf.inc'
      type (x_2grids), intent(in) :: x2grids
      real*8 :: tlatlon(ilonm,jlatm),alatlon(ilonm,jlatm)
      real*8 :: tcub_loc(isd:ied,jsd:jsd)
      real*8 :: tlatlon_dbg(ilonm,jlatm),
     &     alatlon_dbg(ilonm,jlatm)      
      complex*16 :: tlatlon_dd(ilonm,jlatm),alatlon_dd(ilonm,jlatm)
      complex*16 :: tlatlon_tmp(1),tmp_dd,alatlon_tmp(1),xarea_dd(1),
     &     tmp_dd_vec(1)
      integer :: n,icub,jcub,i,j,itile
      character*120:: ofi
      integer ::  status, fid, vid
     
      write(*,*) "CS2LL REPRO"
      write(*,*) "isd,ied,jsd,jed,ilonm,jlatm=",isd,ied,jsd,jed,ilonm
     & ,jlatm
       
 
      do j=1,jlatm
         do i=1,ilonm
            alatlon_dd(i,j) = dcmplx(0.d0,0.d0)
            tlatlon_dd(i,j) = dcmplx(0.d0,0.d0)
         enddo
      enddo

      alatlon_dbg(:,:) = 0d0
      tlatlon_dbg(:,:) = 0d0

      do n=1,x2grids%xgrid%ncells
        icub=x2grids%xgrid%ijcub(1,n)
        jcub=x2grids%xgrid%ijcub(2,n)
        itile=x2grids%xgrid%tile(n)

        if (itile .eq. mytile) then
           if ( ( (icub .le. ie) .and. (icub .ge. is)) .and.
     &          ( (jcub .le. je) .and. (jcub .ge. js)) ) then
c              write(*,*) "inside ie ,is , je, js"
              i=x2grids%xgrid%ijlatlon(1,n)
              j=x2grids%xgrid%ijlatlon(2,n)
              
              alatlon_tmp(1)=alatlon_dd(i,j)
              xarea_dd(1)=dcmplx(x2grids%xgrid%xgrid_area(n), 0.d0)
              call add_dd(xarea_dd,alatlon_tmp,1) 
              alatlon_dd(i,j)=alatlon_tmp(1)     !OK
              
              tlatlon_tmp(1)=tlatlon_dd(i,j)
              call dxd(x2grids%xgrid%xgrid_area(n),
     &             tcub_loc(icub,jcub),tmp_dd)
              tmp_dd_vec(1)=tmp_dd
              call add_dd(tmp_dd_vec,tlatlon_tmp,1)   
              tlatlon_dd(i,j)=tlatlon_tmp(1)
c
c              alatlon_dbg(i,j) = alatlon_dbg(i,j) + xgrid_area(n)
c              tlatlon_dbg(i,j) = tlatlon_dbg(i,j) + xgrid_area(n)
c     &             *tcub_loc(icub,jcub)
c
c              tlatlon_dd(i,j)=dcmplx(1.0d0,0.0d0) !remove this, debug only
c              alatlon_dd(i,j)=dcmplx(1.0d0,0.0d0) !remove this, debug only
           endif
        endif
      enddo

      call sumxpe2d_exact(tlatlon_dd)
      call sumxpe2d_exact(alatlon_dd)

c      write(*,*) "counting"

c
c      call sumxpe(tlatlon_dbg)
c      call sumxpe(alatlon_dbg)
c

c
c     root proc section
c

      if (AM_I_ROOT()) then   
         write(*,*) "HERE"
         do j=1,jlatm
            do i=1,ilonm
               tlatlon(i,j)=real(tlatlon_dd(i,j))
               alatlon(i,j)=real(alatlon_dd(i,j))
            enddo
         enddo
         tlatlon(:,:) = tlatlon(:,:)/alatlon(:,:)
         
c     tlatlon_dbg(:,:) = tlatlon_dbg(:,:)
c     &          /alatlon_dbg(:,:)
         
         
         do j=1,jlatm
            do i=1,ilonm
               write(*,*) "tlat",tlatlon(i,j)
c     &              ,tlatlon_dbg(i,j),"<<<"
            enddo
         enddo
         
c         ofi='tstout.nc'
c         status = nf_open(trim(ofi),nf_write,fid)
c         if (status .ne. NF_NOERR) write(*,*) NF_STRERROR(status)
c         status = nf_inq_varid(fid,'lwup_sfc',vid)
c         write(*,*) NF_STRERROR(status)
c         status = nf_put_var_double(fid,vid,tlatlon)
c         write(*,*) "STATUS",NF_STRERROR(status),"<<"
c         
c         status = nf_close(fid)
      endif

      end subroutine regrid_cs2ll_exact
c*


      subroutine root_regrid_ll2cs(x2grids,tlatlon,tcubglob)

!@sum  Root processor regrids data from lat-lon -> cubbed sphere 
!@auth Denis Gueyffier
      use regrid_com
      implicit none
      type (x_2grids), intent(in) :: x2grids
      real*8 :: tlatlon(ilonm,jlatm),tcubglob(ic,jc,6),
     &     acubglob(ic,jc,6)
      integer :: n,icub,jcub,i,j,itile,icc,jcc,il,jl


      if (AM_I_ROOT()) then   
         write(6,*) "IN ROOT REGRID LL2CS"
         acubglob(:,:,:) = 0d0
         tcubglob(:,:,:) = 0d0
         
         do n=1,x2grids%xgrid%ncells
            itile=x2grids%xgrid%tile(n)
            icc=x2grids%xgrid%ijcub(1,n)
            jcc=x2grids%xgrid%ijcub(2,n)
            il=x2grids%xgrid%ijlatlon(1,n)
            jl=x2grids%xgrid%ijlatlon(2,n)
            
            acubglob(icc,jcc,itile) = acubglob(icc,jcc,itile) 
     &           + x2grids%xgrid%xgrid_area(n)
            tcubglob(icc,jcc,itile) = tcubglob(icc,jcc,itile) 
     &           + x2grids%xgrid%xgrid_area(n)*tlatlon(il,jl)
         enddo
         
         do itile=1,6
            do j=1,jc
               do i=1,ic
                  tcubglob(i,j,itile) = tcubglob(i,j,itile)
     &                 /acubglob(i,j,itile)
               enddo
            enddo
         enddo
         write(6,*) "END ROOT REGRID LL2CS"
      endif
      
      end subroutine root_regrid_ll2cs
c*

         
      subroutine init_regrid(x2grids,imsource,jmsource,
     &     ntilessource,imtarget,jmtarget,ntilestarget)

!@sum  Reads regriding file on root proc, broadcasts the 
!@+    x_grid data to all processes then instanciates locally
!@+    the regrid derived type (regrid type=x_grid plus extra info about
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
      integer :: ierr
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
c         write(imch,'(i10)') ilonm
c         write(jmch,'(i10)') jlatm
c         write(icch,'(i10)') ic
c         write(jcch,'(i10)') jc

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
c*


      subroutine readt_regrid_parallel(x2grids,iunit,name,nskip,
     &     tcubloc,ipos)

!@sum  Read input data on lat-lon grid, regrid to global cubbed sphere grid
!@+    then scatter to all subdomains
!@auth Denis Gueyffier

      use regrid_com, only :ilonm,jlatm,ic,jc,gid,AM_I_ROOT,x_2grids
      use gatscat_mod
      implicit none
      type (x_2grids), intent(in) :: x2grids
      integer, intent(in) :: iunit
      character*16, intent(in) :: name
      integer, intent(in) :: nskip

      real*8 :: tcubloc(is:ie,js:je) !output regridded and scattered data 
      real*4 :: tllr4(ilonm,jlatm)  ! latlon real*4 data read from input file
      real*8 :: tdatall(ilonm,jlatm)  ! latlon real*8 data 
      real*8 :: tcubglob(ic,jc,6)  ! global array 
      real*4 :: X              !dummy arrays
      integer, intent(in) :: ipos
      integer :: n,ierr
      character*80 :: TITLE     


      if (AM_I_ROOT()) then   
         do n=1,ipos-1
            read(UNIT=iunit,IOSTAT=ierr)
         enddo
         read(UNIt=iunit,IOSTAT=ierr) TITLE, (X,n=1,nskip), tllr4
c     convert from real*4 to real*8
         tdatall=tllr4
c     Regrid from lat-lon to cubbed sphere, form global array
         call root_regrid_ll2cs(x2grids,tdatall,tcubglob)
       endif

c     Scatter data to every processor
c      write(6,*) tcubglob
      call unpack_data(tcubglob,tcubloc)
      write(6,*) tcubloc
      end subroutine readt_regrid_parallel
c*


      subroutine dread_regrid_parallel(x2grids,iunit,name,tcubloc)
     
!@sum  Read input data on lat-lon grid, regrid to global cubbed sphere grid
!@+    then scatter to all subdomains
!@auth Denis Gueyffier

      use regrid_com, only :ilonm,jlatm,ic,jc,AM_I_ROOT,x_2grids
      use gatscat_mod
      implicit none
      type (x_2grids), intent(in) :: x2grids
      integer, intent(in) :: iunit
      character*16, intent(in) :: name
      real*8 :: tcubloc(is:ie,js:je) !output regridded and scattered data 
      real*4 :: tllr4(ilonm,jlatm)  ! latlon real*4 data read from input file
      real*8 :: tdatall(ilonm,jlatm)  ! latlon real*8 data 
      real*8 :: tcubglob(ic,jc,6)  ! global array 
      real*4 :: X              !dummy arrays
      integer :: ierr
   
      if (AM_I_ROOT()) then   
         read(UNIt=iunit,IOSTAT=ierr) tllr4
c     convert from real*4 to real*8
         tdatall=tllr4
c     Regrid from lat-lon to cubbed sphere, form global array
         call root_regrid_ll2cs(x2grids,tdatall,tcubglob)
       endif

c     Scatter data to every processor
c      write(6,*) tcubglob
      call unpack_data(tcubglob,tcubloc)
      write(6,*) tcubloc
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


      subroutine sumxpe2d_exact(arr)
!@sum This code calculates a reproducible non rank reducing sum 
!@+   independent on number of processors using double-double arithmetics
!@auth Denis Gueyffier
      use regrid_com,only: ilonm,jlatm
      implicit none
      include 'mpif.h'
      integer :: ierr,arr_size,myPE
      integer :: MPI_SUMDD
      external add_DD
      complex*16, dimension(ilonm,jlatm) :: arr
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

      USE regrid_com, only:ilonm,jlatm,ic,jc
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
      npx = ic; npy = npx ! 1x1 resolution
      ng = 3 ! number of ghost zones required

      call fv_domain_decomp(npx+1,npy+1,ntiles,ng,grid_type)
c      call mp_stop()

      ndims = 2
      npz = 5
      call init_grid(atm,grid_name,grid_file,
     &     npx+1, npy+1, npz, ndims, ntiles, ng)

      non_ortho=.true.
      call grid_utils_init(Atm, npx+1, npy+1, npz, Atm%grid, Atm%agrid,
     &     area, area_c, cosa, sina, dx, dy, dxa, dya, non_ortho,
     &     uniform_ppm, grid_type, c2l_ord)

      end subroutine init_csgrid_debug
c****



