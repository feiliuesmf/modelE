      program ijkdag
C****
C****
      use ncinp
      use ncout
      implicit none
      integer :: nargs,iargc
      character(len=80) :: stop_str,cmd_str
      real, dimension(:), allocatable :: lon_gcm,lat_gcm,plm_gcm
      real, dimension(:,:,:), allocatable :: arrn,arrb,dpb
      character(len=20) :: var_name,acc_name,dim_name

      call getarg(0,cmd_str)
      nargs = iargc()
      if(nargs.ne.2) then
         stop_str = 'usage: '//trim(cmd_str)//' acc_file out_file'
         stop trim(stop_str)
      endif
      call getarg(1,accfile)
      call getarg(2,outfile)

! open the acc file
      call open_acc

! allocate space using these dimensions
      allocate(lon_gcm(im))
      allocate(lat_gcm(jm))
      allocate(plm_gcm(lm))
      allocate(arrn(im,jm,lm))
      allocate(dpb(im,jm,lm))
      allocate(arrb(im,jm-1,lm))
! get lon,lat,p coordinates
      acc_name='lonb'; call getacc(acc_name,lon_gcm)
      acc_name='latb'; call getacc(acc_name,lat_gcm)
      lat_gcm(1:jm-1) = lat_gcm(2:jm) ! shift off nonexistent 1st latitude
      acc_name='p'; call getacc(acc_name,plm_gcm)

! define output file
      call open_out
      dim_name='longitude'; call def_dim_out(dim_name,im)
      dim_name='latitude'; call def_dim_out(dim_name,jm-1)
      dim_name='p'; call def_dim_out(dim_name,lm)

      ndims_out = 1
      dim_name='longitude'; call set_dim_out(dim_name,1)
      units='degrees_east'
      var_name='longitude';call wrtarr(var_name,lon_gcm)
      dim_name='latitude'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='latitude';call wrtarr(var_name,lat_gcm)
      dim_name='p'; call set_dim_out(dim_name,1)
      units='mb'; long_name='Pressure'
      var_name='p';call wrtarr(var_name,plm_gcm)

      ndims_out = 3
      dim_name='longitude'; call set_dim_out(dim_name,1)
      dim_name='latitude'; call set_dim_out(dim_name,2)
      dim_name='p'; call set_dim_out(dim_name,3)

! get DPB array from acc file
      acc_name='DPB'; call getacc(acc_name,dpb)

! u-velocity
      acc_name='UDPB'; call getacc(acc_name,arrn)
      call scale(arrb,arrn,dpb,im,jm,lm,missing,1.)
      var_name='u'; call wrtarr(var_name,arrb)

! v-velocity
      acc_name='VDPB'; call getacc(acc_name,arrn)
      call scale(arrb,arrn,dpb,im,jm,lm,missing,1.)
      var_name='v'; call wrtarr(var_name,arrb)

! temperature
      acc_name='TDPBx4'; call getacc(acc_name,arrn)
      call scale(arrb,arrn,dpb,im,jm,lm,missing,.25)
      var_name='t'; call wrtarr(var_name,arrb)
! and height
      acc_name='DSEDPBx4'; call getacc(acc_name,arrn)
      where(dpb(:,2:jm,:).gt.0.)
         arrb(:,:,:) = (.25*arrn(:,2:jm,:)/dpb(:,2:jm,:) -
     &        sha*arrb(:,:,:))/grav
      elsewhere
         arrb(:,:,:) = missing
      end where
      var_name='z'; call wrtarr(var_name,arrb)

! humidity
      acc_name='QDPBx4'; call getacc(acc_name,arrn)
      call scale(arrb,arrn,dpb,im,jm,lm,missing,.25)
      var_name='q'; call wrtarr(var_name,arrb)

      call close_acc
      call close_out

      end program ijkdag

      subroutine scale(arrb,arrn,dpb,im,jm,lm,missing,fac)
      implicit none
      real, dimension(im,jm-1,lm) :: arrb
      real, dimension(im,jm,lm) :: arrn,dpb
      integer :: im,jm,lm
      real :: missing,fac
      where(dpb(:,2:jm,:).gt.0.)
         arrb(:,:,:) = fac*arrn(:,2:jm,:)/dpb(:,2:jm,:)
      elsewhere
         arrb(:,:,:) = missing
      end where
      return
      end subroutine scale
