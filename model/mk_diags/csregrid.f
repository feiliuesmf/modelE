      program csregrid
      implicit none
      character(len=80) :: varname,csfile,llfile,remapfile
      integer :: nargs,iargc
      include 'netcdf.inc'
      integer, dimension(7) :: srt,icnt,ocnt,dids,dsizes,kmod
      integer :: status,ifid,ofid,ivid,ovid,imdid,jmdid,tldid,ndims
      integer :: imcub,imcub2,jmcub2,imlon,imlon2,jmlat,jmlat2,ncells,
     &     ntiles,n,idim,jdim,kacc,k
      integer, dimension(:), allocatable :: tcub
      integer, dimension(:,:), allocatable :: ijcub,ijll
      real*8, dimension(:), allocatable :: xgrid_area,xp,yp
      real*8, dimension(:,:,:), allocatable :: arrcs
      real*8, dimension(:,:), allocatable :: arrll,areall

      nargs = iargc()
      if(nargs.ne.4) then
        write(6,*) 'usage: csregrid remapfile csfile llfile varname'
        stop
      endif

c
c get the command line arguments
c
      call getarg(1,remapfile)
      call getarg(2,csfile)
      call getarg(3,llfile)
      call getarg(4,varname)

c
c read the remap file
c
      call handle_err(nf_open(remapfile,nf_nowrite,ifid),
     &     'opening '//trim(remapfile))
      call get_dimsize(ifid,'imcub',imcub)
      call get_dimsize(ifid,'lon',imlon)
      call get_dimsize(ifid,'lat',jmlat)
      call get_dimsize(ifid,'ncells',ncells)
      allocate(tcub(ncells),ijcub(2,ncells),ijll(2,ncells))
      allocate(xgrid_area(ncells),xp(ncells),yp(ncells))
      call get_var_int(ifid,'tile1',tcub)
      call get_var_int(ifid,'tile1_cell',ijcub)
      call get_var_int(ifid,'tile2_cell',ijll)
      call get_var_r8(ifid,'xgrid_area',xgrid_area)
      call get_var_r8(ifid,'xprime',xp)
      call get_var_r8(ifid,'yprime',yp)
      status = nf_close(ifid)

c
c calculate areas on the latlon grid
c
      allocate(areall(imlon,jmlat))
      do n=1,ncells
        areall(ijll(1,n),ijll(2,n)) = areall(ijll(1,n),ijll(2,n))
     &       + xgrid_area(n)
      enddo

c
c open the input and output files, check consistency of dimensions
c
      call handle_err(nf_open(csfile,nf_nowrite,ifid),
     &     'opening '//trim(csfile))
      call get_dimsize(ifid,'im',imcub2)
      call get_dimsize(ifid,'jm',jmcub2)
      call get_dimsize(ifid,'tile',ntiles)
      if(ntiles.ne.6) stop 'ntiles != 6'
      if(imcub2.ne.jmcub2) then
        write(6,*) trim(csfile),
     &       ' does not appear to be a cubed sphere file: im,jm = ',
     &       imcub2,jmcub2
        stop
      endif
      if(imcub2.ne.imcub) then
        write(6,*) 'mismatched sizes: '
        write(6,*) 'remap cube size = ',imcub
        write(6,*) 'input cube suze = ',imcub2
        stop
      endif

      call handle_err(nf_open(llfile,nf_write,ofid),
     &     'opening '//trim(llfile))
      call get_dimsize(ofid,'lon',imlon2)
      call get_dimsize(ofid,'lat',jmlat2)
      if(imlon2.ne.imlon .or. jmlat2.ne.jmlat) then
        write(6,*) 'mismatched grid sizes: '
        write(6,*) 'remap imlon,jmlat = ',imlon,jmlat
        write(6,*) 'output imlon,jmlat = ',imlon2,jmlat2
        stop
      endif

      call handle_err(nf_inq_varid(ifid,varname,ivid),
     &     trim(varname)//' not present in '//trim(csfile))

      call handle_err(nf_inq_varid(ofid,varname,ovid),
     &     trim(varname)//' not present in '//trim(llfile))

c
c setup for the loop over the fields to be regridded
c
      status = nf_inq_dimid(ifid,'im',imdid)
      status = nf_inq_dimid(ifid,'jm',jmdid)
      status = nf_inq_dimid(ifid,'tile',tldid)
      status = nf_inq_varndims(ifid,ivid,ndims)
      status = nf_inq_vardimid(ifid,ivid,dids)
      if(dids(ndims).ne.tldid) stop 'tile must be the last dimension'
      idim = 0
      jdim = 0
      do n=1,ndims
        if(dids(n).eq.imdid) idim=n
        if(dids(n).eq.jmdid) jdim=n
        status = nf_inq_dimlen(ifid,dids(n),dsizes(n))
      enddo
      if(idim.eq.0 .or. jdim.eq.0 .or. jdim.ne.idim+1)
     &     stop 'array in input file has no im,jm dims'
      srt = 1
      icnt(:) = 1; icnt(ndims) = 6
      icnt(idim:jdim) = (/ imcub, imcub /)
      ocnt(:) = 1
      ocnt(idim:jdim) = (/ imlon, jmlat /)
      kacc = product(dsizes(1:ndims))/(imcub*imcub*6)
      k = 1
      do n=1,ndims-1
        if(n.eq.idim .or. n.eq.jdim) cycle
        kmod(n) = k
        k = k*dsizes(n)
      enddo

c
c loop over fields
c
      allocate(arrcs(imcub,imcub,6),arrll(imlon,jmlat))
      do k=1,kacc
c read the cubed sphere data
        status = nf_get_vara_double(ifid,ivid,srt,icnt,arrcs)
c regrid the data using either first- or second-order method
c todo: take the choice of method as a command-line argument,
c       but need to be able to indicate which quantities require
c       non-negativity constraints for 2nd order method
c first-order method
        call regrid1(arrcs,imcub,arrll,areall,imlon,jmlat,
     &     ijcub,ijll,tcub,xgrid_area,ncells)
c second-order method
c        call regrid2(arrcs,imcub,arrll,areall,imlon,jmlat,
c     &     ijcub,ijll,tcub,xgrid_area,xp,yp,ncells)
c write the latlon data
        status = nf_put_vara_double(ofid,ovid,srt,ocnt,arrll)
        do n=1,ndims-1 ! increment the start vector
          if(n.eq.idim .or. n.eq.jdim) cycle
          if(mod(k,kmod(n)).eq.0) then
            srt(n) = srt(n) + 1
            if(srt(n).gt.dsizes(n)) srt(n)=1
          endif
        enddo
      enddo

c
c close input and output files
c
      status = nf_close(ifid)
      status = nf_close(ofid)


      end program csregrid

      subroutine regrid1(qcs,imcub,qll,areall,imlon,jmlat,
     &     ijcs,ijll,tile,xgrid_area,ncells)
      implicit none
      integer :: imcub,imlon,jmlat,ncells
      real*8, dimension(imcub,imcub,6) :: qcs
      real*8, dimension(imlon,jmlat) :: qll,areall
      integer, dimension(2,ncells) :: ijcs,ijll
      integer, dimension(ncells) :: tile
      real*8, dimension(ncells) :: xgrid_area
      integer :: ics,jcs,ill,jll,n
      qll(:,:) = 0.
      do n=1,ncells
        ics = ijcs(1,n)
        jcs = ijcs(2,n)
        ill = ijll(1,n)
        jll = ijll(2,n)
        qll(ill,jll) = qll(ill,jll) + xgrid_area(n)*
     &       qcs(ics,jcs,tile(n))
      enddo
      qll(:,:) = qll(:,:)/areall(:,:)
      return
      end subroutine regrid1

      subroutine regrid2(qcs,imcub,qll,areall,imlon,jmlat,
     &     ijcs,ijll,tile,xgrid_area,xp,yp,ncells)
      implicit none
      integer :: imcub,imlon,jmlat,ncells
      real*8, dimension(imcub,imcub,6) :: qcs
      real*8, dimension(imlon,jmlat) :: qll,areall
      integer, dimension(2,ncells) :: ijcs,ijll
      integer, dimension(ncells) :: tile
      real*8, dimension(ncells) :: xgrid_area,xp,yp
      real*8, dimension(:,:,:), allocatable :: dqi,dqj
      integer :: ics,jcs,ill,jll,n,itile
      real*8 :: dx,gradfac
      allocate(dqi(imcub,imcub,6),dqj(imcub,imcub,6))
      dx = 2./real(imcub,kind=8)
      gradfac = 1./(2.*dx)
      do itile=1,6
        do jcs=1,imcub
          dqi(1,jcs,itile) = 0.
          do ics=2,imcub-1
            dqi(ics,jcs,itile) = 
     &           (qcs(ics+1,jcs,itile)-qcs(ics-1,jcs,itile))*gradfac
          enddo
          dqi(imcub,jcs,itile) = 0.
        enddo
        dqj(:,1,itile) = 0.
        do jcs=2,imcub-1
          do ics=1,imcub
            dqj(ics,jcs,itile) =
     &           (qcs(ics,jcs+1,itile)-qcs(ics,jcs-1,itile))*gradfac
          enddo
        enddo
        dqj(:,imcub,itile) = 0.
      enddo
      qll(:,:) = 0.
      do n=1,ncells
        ics = ijcs(1,n)
        jcs = ijcs(2,n)
        ill = ijll(1,n)
        jll = ijll(2,n)
c todo: if q is non-negative, do a recon pass over cs cells
c and set dqi=dqj=0 in cells in which any qcs+dqi*xp+dqj*yp < 0
        qll(ill,jll) = qll(ill,jll) + (
     &       qcs(ics,jcs,tile(n))
     &      +dqi(ics,jcs,tile(n))*xp(n)
     &      +dqj(ics,jcs,tile(n))*yp(n)
     &       )*xgrid_area(n)
      enddo
      qll(:,:) = qll(:,:)/areall(:,:)
      deallocate(dqi,dqj)
      return
      end subroutine regrid2

      subroutine get_var_r8(fid,var_name,var)
      implicit none
      include 'netcdf.inc'
      integer :: fid
      character(len=*) :: var_name
      real*8 :: var(1)
      integer :: status,varid
      status = nf_inq_varid(fid,var_name,varid)
      if(status.ne.nf_noerr) then
        write(6,*) 'nonexistent variable ',trim(var_name)
        stop
      endif
      status = nf_get_var_double(fid,varid,var)
      return
      end subroutine get_var_r8
