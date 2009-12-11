!@sum ncarithm performs the following arithmetic operations on netcdf files:
!@    fileout.nc=a1*file1.nc+a2*file2.nc
!@+   on all variables inside file1.nc and file2.nc which have same name 
!@+   and which have format == cformat  
!@auth D. Gueyffier 
      program ncarithm
      use ncio, only : defvar, write_data
      implicit none
      include 'netcdf.inc'
      character*150 :: file1,file2,fileout,cformat,ca1,ca2,ccopy
      character*80 :: title
      real*4, allocatable :: sij4(:,:), sij4_2(:,:), tij4(:,:)   
      real*8, allocatable :: tij(:,:)
      real*4, allocatable :: sijl4(:,:,:), sijl4_2(:,:,:), tijl4(:,:,:) 
      real*8, allocatable :: tijl(:,:,:)
      real*4 :: a1,a2
      integer :: im,jm,lm,n,maxrec,i,imax,status,fid1,fid2,fidout,vid,
     &     vid2
      integer :: offset(10),dim(10)
      integer :: ndims, nvars, ngatts, itype, ndim, natt, UNLIMDIMID,
     &     idlm,idjm,idim,is
      real*8, allocatable :: lon(:),lat(:),ill(:)
      real*4, allocatable :: lon4(:),lat4(:),ill4(:)
      character*20 :: cdim(10),cform
      character*80 :: cval,vname
      integer ::  ishape(3)   ! variable shape
      integer :: dtype,nd,shp(7)


      if(iargc().ne.6) then
        write(6,*) 
     &        'usage: ncarithm cformat a1 a2 file1 file2 fileout'
        stop
      endif

      call getarg(1,cformat)
      call getarg(2,ca1)
      call getarg(3,ca2)
      call getarg(4,file1)
      call getarg(5,file2)
      call getarg(6,fileout)

      read(ca1,'(F20.10)') a1
      read(ca2,'(F20.10)') a2

      write(*,*) "a1 a2=",a1,a2
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

      status = nf_open(trim(file1),nf_nowrite,fid1)
      if (status .ne. nf_noerr) write(*,*) 
     *     "UNABLE TO OPEN FILE1 ",trim(file1)

      status = nf_open(trim(file2),nf_nowrite,fid2)
      if (status .ne. nf_noerr) write(*,*) 
     *     "UNABLE TO OPEN FILE2 ",trim(file2)

      status = nf_create(trim(fileout),nf_clobber,fidout)
      if(status .ne. nf_noerr) write(*,*)
     &     "ERROR CREATING ",trim(fileout)

c     read dimensions from file1
      do i=1,imax
         status = nf_inq_dimid(fid1,cdim(i),vid)
         status = nf_inq_dimlen(fid1,vid,dim(i))
      enddo
     
c     find if format is 'ij', 'ijl'
      if (imax .eq. 3) then
         if (cdim(imax) .eq. 'lon' .and. cdim(imax-1) .eq. 'lat') then 
            cform='ijl'
            im=dim(3)
            jm=dim(2)
            lm=dim(1)
         endif
      elseif (imax .eq. 2) then
         cform='ij'
         im=dim(2)
         jm=dim(1)
      endif

      write(*,*) "array format=",cform
      if (cform .eq. 'ijl') write(*,*) "i j l=",im,jm,lm
      if (cform .eq. 'ij') write(*,*) "i j=",im,jm

c***  first loop to find ids of the different dimensions
      idim=0;idjm=0;idlm=0
      status = nf_inq(fid1, ndims, nvars, ngatts, UNLIMDIMID)
      do i = 1, nvars
         status = nf_inq_var(fid1, i, cval, itype, ndim, ishape, natt)

         if (cform .eq. 'ijl') then
            if (cval .eq. 'lon')   idim=ishape(1)
            if (cval .eq. 'lat')   idjm=ishape(1)
            if (cval .eq. cdim(1)) idlm=ishape(1)
         endif
         if (cform .eq. 'ij') then
            if (cval .eq. 'lon')   idim=ishape(1)
            if (cval .eq. 'lat')   idjm=ishape(1)
         endif
      enddo
      write(*,*) "idim idjm idlm",idim,idjm,idlm

c***  2nd loop to define variables which have same name and format == cformat 
      do i = 1, nvars
         status = nf_inq_var(fid1, i, cval, itype, ndim, ishape, natt)
         
         if (cform .eq. 'ijl') then
            if (cval .eq. 'lon') then
               allocate(lon(im),lon4(im))
               call defvar(fidout,im,jm,1,lon,'lon(lon)')
            endif
            if (cval .eq. 'lat') then
               allocate(lat(jm),lat4(jm))
               call defvar(fidout,im,jm,1,lat,'lat(lat)')
            endif
            if (cval .eq. cdim(1)) then 
               allocate(ill(lm),ill4(lm))
               vname=trim(cval)//'('//trim(cval)//')'
               call defvar(fidout,im,jm,1,ill,vname)
            endif
         endif
         if (cform .eq. 'ij') then
            if (cval .eq. 'lon') then
               allocate(lon(im),lon4(im))
               call defvar(fidout,im,jm,1,lon,'lon(lon)')
            endif
            if (cval .eq. 'lat') then
               allocate(lat(jm),lat4(jm))
               call defvar(fidout,im,jm,1,lat,'lat(lat)')
            endif
         endif

         if(itype.eq.nf_char) cycle         
         
         status = nf_inq_varid(fid2,cval,vid2)
         
         if(status.ne.nf_noerr) then
            write(6,*) 'file1 array '//trim(cval)//
     &           ' is absent in file2: skipping'
            cycle
         endif

         if (cform .eq. 'ijl') then
            if (ishape(1) .eq. idim .and. 
     &          ishape(2) .eq. idjm .and. 
     &          ishape(3) .eq. idlm .and. ndim .eq. 3) then
               
               allocate(tijl(im,jm,lm))
               
c     define variable to netcdf output file
               vname=trim(cval)//'(lon,lat,'//cdim(1)//')'
               call defvar(fidout,im,jm,1,tijl,trim(vname))
               deallocate(tijl)
            endif
         endif

         if (cform .eq. 'ij') then
            if (ishape(1) .eq. idim .and. 
     &          ishape(2) .eq. idjm .and. ndim .eq. 2) then

               allocate(tij(im,jm))
               
c     define variable to netcdf output file
               vname=trim(cval)//'(lon,lat)'
               call defvar(fidout,im,jm,1,tij,trim(vname))
               deallocate(tij)
            endif
         endif       
      enddo      !end of 2nd loop

      status = nf_close(fid1)
      status = nf_close(fid2)
      status = nf_enddef(fidout)

      if (status .ne. NF_NOERR) 
     &     write(*,*) "Problem with enddef"

c***  3rd loop to perform arithmetic operations and write variables 
      status = nf_open(trim(file1),nf_nowrite,fid1)
      status = nf_open(trim(file2),nf_nowrite,fid2)
      status = nf_inq(fid1, ndims, nvars, ngatts, UNLIMDIMID)

      do i = 1, nvars
         status = nf_inq_var(fid1, i, cval, itype, ndim, ishape, natt)

         if(itype.eq.nf_char) cycle         
         
         status = nf_inq_varid(fid2,cval,vid2)
         
         if(status.ne.nf_noerr) then
            write(6,*) 'file1 array '//trim(cval)//
     &           ' is absent in file2: skipping'
            cycle
         endif

         if (cform .eq. 'ijl') then
c     find ids of the different dimensions
            if (cval .eq. 'lon') then
               status = nf_inq_varid(fid1,cval,vid)
               status = nf_get_var_real(fid1,vid,lon4)
               idim=ishape(1)  
               call write_data(fidout,trim(cval),lon4)
            endif
            if (cval .eq. 'lat') then
               idjm=ishape(1)
               status = nf_inq_varid(fid1,cval,vid)
               status = nf_get_var_real(fid1,vid,lat4)
               call write_data(fidout,trim(cval),lat4)
            endif
            if (cval .eq. cdim(1)) then
               idlm=ishape(1)
               status = nf_inq_varid(fid1,cval,vid)
               status = nf_get_var_real(fid1,vid,ill4)               
               call write_data(fidout,trim(cval),ill4)
            endif
            
c     extract values of variables having correct format 
            if (ishape(1) .eq. idim .and. 
     &           ishape(2) .eq. idjm .and. 
     &           ishape(3) .eq. idlm .and. ndim .eq. 3) then
               
               allocate(sijl4(im,jm,lm),tijl4(im,jm,lm))
               allocate(sijl4_2(im,jm,lm))
               
               status = nf_inq_varid(fid1,cval,vid)
               status = nf_inq_varid(fid2,cval,vid2)
               status = nf_get_var_real(fid1,vid,sijl4)
               status = nf_get_var_real(fid2,vid2,sijl4_2)
               
               tijl4=a1*sijl4+a2*sijl4_2

c     write variable to netcdf output file  
               write(*,*) "write ",adjustl(trim(cval))," in output file"
               call write_data(fidout,trim(cval),tijl4)

               deallocate(sijl4, sijl4_2, tijl4)
            endif
         endif

         if (cform .eq. 'ij') then
c     find ids of the different dimensions
            if (cval .eq. 'lon') then
               idim=ishape(1)
               call write_data(fidout,trim(cval),lon4)
            endif
            if (cval .eq. 'lat') then
               idjm=ishape(1)
               call write_data(fidout,trim(cval),lat4)
            endif
            
c     extract values of variables having correct format 
            if (ishape(1) .eq. idjm .and. 
     &           ishape(2) .eq. idim .and. ndim .eq. 2) then

               allocate(sij4(im,jm),tij4(im,jm))
               allocate(sij4_2(im,jm))
               
               status = nf_inq_varid(fid1,cval,vid)
               status = nf_inq_varid(fid2,cval,vid2)
               status = nf_get_var_real(fid1,vid,sij4)
               status = nf_get_var_real(fid2,vid2,sij4_2)

               tij4=a1*sij4+a2*sij4_2

c     write variable to netcdf output file  

c     write each remaped variable to netcdf output file  
               write(*,*) "write ",adjustl(trim(cval))," in output file"
               call write_data(fidout,trim(cval),tij4)

               deallocate(sij4, sij4_2, tij4)
               deallocate(lon,lat,lon4,lat4)
            endif
         endif       

      enddo    !end of 3rd loop

      if (cform .eq. 'ijl') deallocate(lon,lat,ill,lon4,lat4,ill4)
      if (cform .eq. 'ij') deallocate(lon,lat,lon4,lat4)

      status = nf_close(fidout)
      status = nf_close(fid1)
      status = nf_close(fid2)
      
      write(6,*) "wrote:",trim(fileout),">>>>>>>>>>>>"
      write(6,*) ""

      end program

