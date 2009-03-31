      program regridAIJ
!@sum regridding aij diagnostics from CS to latlon grid
!@auth Denis Gueyffier
      use offregrid_com
      implicit none
      include 'netcdf.inc'
      type (x_2gridsoff) :: xcs2ll
      real*4, allocatable :: AIJ_in(:,:,:,:),AIJ_out(:,:,:)
      real*8, allocatable :: tsource(:,:,:)
      real*8, allocatable :: ttargglob(:,:,:)
      integer, allocatable :: lon(:),lat(:)
      character*80 filein,fileout
      integer :: iu_AIJ,iuout,ims,jms,nts,imt,jmt,ntt
      integer :: i,j,k,fid,status,ntiles,im,jm
      integer :: ntilessource,ntilestarget,
     &     kaij,fidin,fidout,nargs
      integer, external :: iargc
      nargs = iargc()
      call getarg(1,fileout)

c***  read aij defined on CS grid
      filein="PARTIAL.accEtstcs3.nc"
      status = nf_open(filein,nf_nowrite,fidin)
      status = nf_open(trim(fileout),nf_write,fidout)

      call get_dimsize(fidin,'kaij',kaij)
      call get_dimsize(fidin,'im',ims)
      call get_dimsize(fidin,'jm',jms)
      write(*,*) "ims, jms, kaij=",ims,jms,kaij

      call get_dimsize(fidout,'lon',imt)
      call get_dimsize(fidout,'lat',jmt)
      write(*,*) "imt, jmt=",imt,jmt

      allocate (AIJ_in(ims,jms,kaij,6),AIJ_out(imt,jmt,kaij))
      allocate (tsource(ims,jms,6),ttargglob(imt,jmt,1) )
      allocate (lon(imt),lat(jmt))

      call get_var_real(fidin,'aij',AIJ_in)
      status= nf_close(fidin)


      call init_offregrid(xcs2ll,ims,jms,6,
     &     imt,jmt,1)

      do k=1,kaij
         tsource(:,:,:)=AIJ_in(:,:,k,:)
         call offline_regrid(xcs2ll,tsource,ttargglob)
         AIJ_out(:,:,k)=ttargglob(:,:,1)
c         write(*,*) "k=",AIJ_out(:,:,k)
      enddo
     
      do i=1,imt
      lon(i)=i
      enddo
      do j=1,jmt
      lat(j)=j
      enddo

c***  Write regridded aij in netcdf file
      call put_var_real(fidout,'aij',AIJ_out)
      call put_var_int(fidout,'lon',lon)
      call put_var_int(fidout,'lat',lat)

      status = nf_close(fidout)

      deallocate (AIJ_in,AIJ_out,tsource,ttargglob)
      
      end program regridAIJ
c*

