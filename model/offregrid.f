c ifort offregrid.f -o offregrid -convert big_endian -m64 -I/usr/local/netcdf-64bits-ifort-gcc/include -L/usr/local/netcdf-64bits-ifort-gcc/lib -lnetcdf

      program offregrid 
      implicit none
      character*200 :: wtfile,ofile,infi,outfile
      character*200 :: outunformat1,outunformat2,outunformat3 
      character*200 :: outunformat4,outunformat5,outunformat6

      character*1 :: ostr
      character*3 :: irstr

      include 'netcdf.inc'
      integer, parameter :: ncells=31344,nmax=500
      integer, parameter :: im=72,jm=46,imc=48,jmc=48
      real*8 :: xgrid_area(ncells),tot_area
      real*8 :: tcub(imc,jmc,6),acub(imc,jmc,6),area_tile(6)
      real*8 :: tll(im,jm,nmax)
      integer :: tile(ncells)
      integer, dimension(2,ncells) :: ijcub,ijlatlon
      integer :: status,vid,fid,itile,srt(2),cnt(2),n,i,j,ic,jc,il,jl
      character*80 TITLE(nmax)
      real*4 data(im,jm,nmax),tout(imc,jmc,6)
      integer iu1,iu2,iu3,iu4,iu5,iu6,irec,recmax,ir,iunit

c read weights
      wtfile='/Users/dgueyffier/fregrid/remapC48N4x5.nc'
      status = nf_open(trim(wtfile),nf_nowrite,fid)
      if (status .ne. NF_NOERR) write(*,*) "UNABLE TO OPEN REMAP FILE"
      status = nf_inq_varid(fid,'xgrid_area',vid)
      status = nf_get_var_double(fid,vid,xgrid_area)
      status = nf_inq_varid(fid,'tile1',vid)
      status = nf_get_var_int(fid,vid,tile)
      status = nf_inq_varid(fid,'tile1_cell',vid)
      status = nf_get_var_int(fid,vid,ijcub)
      status = nf_inq_varid(fid,'tile2_cell',vid)
      status = nf_get_var_int(fid,vid,ijlatlon)
      status = nf_close(fid)


c Read data on latlon grid
c      infi='/Users/dgueyffier/fregrid/Z72X46N.cor4_nocasp'
      infi='/Users/dgueyffier/fregrid/V72X46.1.cor2_no_crops.ext'
      iunit=15

      open( iunit, FILE=infi, 
     &     FORM='unformatted', STATUS='old')

      irec=1

      do
      read(unit=iunit,END=30) TITLE(irec), data(:,:,irec)
 
      write(*,*) TITLE(irec)

      do i=1,im
         do j=1,jm
            tll(i,j,irec)= data(i,j,irec)
         enddo
      enddo

      irec=irec+1

      enddo
 30   continue

      close(iunit)

      write(*,*) "HERE"
      recmax=irec
c
c     prepare unformatted output
c     
      ofile='/Users/dgueyffier/fregrid/output_fregrid/outll4x5-'
      outunformat1=trim(ofile)//'1.dat'            
      outunformat2=trim(ofile)//'2.dat'            
      outunformat3=trim(ofile)//'3.dat'            
      outunformat4=trim(ofile)//'4.dat'            
      outunformat5=trim(ofile)//'5.dat'            
      outunformat6=trim(ofile)//'6.dat'            
      iu1=21
      iu2=22
      iu3=23
      iu4=24
      iu5=25
      iu6=26
      
      open( iu1, FILE=outunformat1, 
     &     FORM='unformatted', STATUS='new')
      open( iu2, FILE=outunformat2, 
     &     FORM='unformatted', STATUS='new')
      open( iu3, FILE=outunformat3, 
     &     FORM='unformatted', STATUS='new')
      open( iu4, FILE=outunformat4, 
     &     FORM='unformatted', STATUS='new')
      open( iu5, FILE=outunformat5, 
     &     FORM='unformatted', STATUS='new')
      open( iu6, FILE=outunformat6, 
     &     FORM='unformatted', STATUS='new')
      
      do ir=1,recmax
         
c     Regrid
         acub(:,:,:) = 0d0
         tcub(:,:,:) = 0d0
         area_tile(:) = 0d0
         do n=1,ncells
            itile=tile(n)
            ic=ijcub(1,n)
            jc=ijcub(2,n)
            il=ijlatlon(1,n)
            jl=ijlatlon(2,n)
c     write(*,*) tll(il,jl)
            area_tile(itile) = area_tile(itile) + xgrid_area(n)
            acub(ic,jc,itile) = acub(ic,jc,itile) + xgrid_area(n)
            tcub(ic,jc,itile) = tcub(ic,jc,itile) 
     &           + xgrid_area(n)*tll(il,jl,ir)
         enddo
         do itile=1,6
            do j=1,jmc
               do i=1,imc
                  tcub(i,j,itile) = tcub(i,j,itile)/acub(i,j,itile)
               enddo
            enddo
         enddo

c     Output area of each tile
         tot_area=0d0
         do itile=1,6
            write(*,*) "area tile",itile,"=",area_tile(itile)
            tot_area=tot_area+area_tile(itile)
         enddo
         
         write(*,*) "tot area=",tot_area
         
c         
c     Write to netcdf output file
c         
         
         do itile=1,6
            write(ostr,'(i1)') itile
            if (ir .lt. 10) then
               write(irstr,'(i1)') ir
            elseif (ir .lt. 100) then
               write(irstr,'(i2)') ir
            else
               write(irstr,'(i3)') ir
            endif
            outfile=trim(ofile)//ostr//'-'//trim(irstr)
     *           //'.nc'
            
            write(*,*) "ni=",im," nj=",jm
            write(*,*) "outfile =",outfile
            status = nf_open(trim(outfile),nf_write,fid)
            if (status .ne. NF_NOERR) write(*,*) NF_STRERROR(status)
            status = nf_inq_varid(fid,'zsurf',vid)
            write(*,*) NF_STRERROR(status)
            status = nf_put_var_double(fid,vid,tcub(:,:,itile))
            write(*,*) NF_STRERROR(status)
            status = nf_close(fid)
            
         enddo

c
c     write to unformated output
c

         tout(:,:,:)=tcub(:,:,:)

         write(unit=iu1) TITLE(ir), tout(:,:,1)
         write(unit=iu2) TITLE(ir), tout(:,:,2)
         write(unit=iu3) TITLE(ir), tout(:,:,3)
         write(unit=iu4) TITLE(ir), tout(:,:,4)
         write(unit=iu5) TITLE(ir), tout(:,:,5)
         write(unit=iu6) TITLE(ir), tout(:,:,6)
         write(*,*) "ir=",ir
                 
      enddo

      close(iu1)
      close(iu2)
      close(iu3)
      close(iu4)
      close(iu5)
      close(iu6)
      

      end program offregrid

