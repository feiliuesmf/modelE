c f90 -c -64 NCACC.f
      subroutine write_nc_acc(fileout)
      use E001M12_COM
      use GEOM, lat_radians=>lat,lonab_deg=>lon, dlat_radians=>dlat
      use DAGCOM
      implicit none
      include '/usr/local/netcdf-3.4/include/netcdf.inc'

      integer :: varid,status,ncid
      character(len=80) :: fileout
      integer, dimension(4) :: dimids

      integer ::
     &     lon_did,lat_did,lonb_did,latb_did,sig_did,sige_did,
     &     str_levid,plm_did,ple_did,ntype_did

      character(len=20) :: var_name
      double precision, dimension(jm,ntype) :: aj_tmp

      integer :: i,j,k

      double precision, dimension(jm) :: lat_deg,latb_deg
      double precision, dimension(lm) :: plm_gcm
      double precision, dimension(lm+1) :: ple_gcm

      lat_deg = lat_radians*360./twopi
      latb_deg = lat_deg - dlat_radians*180./twopi
      plm_gcm = (psf-ptop)*sig(1:lm)+ptop
      ple_gcm = (psf-ptop)*sige(1:lm+1)+ptop

      call def_acc

      status = nf_create (trim(fileout), nf_clobber, ncid)
c-----------------------------------------------------------------------
c write global header
c-----------------------------------------------------------------------
      status = nf_put_att_text(ncid,nf_global,'XLABEL',132,xlabel)
      status = nf_put_att_text(ncid,nf_global,'AMON',  4,amon)
      status = nf_put_att_text(ncid,nf_global,'AMON0',  4,amon0)
      call ncwrt_iparm(iparm_name,ncid,iparm,100)
      call ncwrt_dparm(dparm_name,ncid,dparm,100)
c-----------------------------------------------------------------------
c define some popular dimensions
c-----------------------------------------------------------------------
      status = nf_def_dim(ncid, 'longitude', im,          lon_did)
      status = nf_def_dim(ncid, 'latitude' ,jm,          lat_did)
      status = nf_def_dim(ncid, 'lonb', im,          lonb_did)
      status = nf_def_dim(ncid, 'latb' ,jm,        latb_did)
      status = nf_def_dim(ncid, 'sig',   lm,      sig_did)
      status = nf_def_dim(ncid, 'sige', lm+1,      sige_did)
      status = nf_def_dim(ncid, 'p',   lm,      plm_did)
      status = nf_def_dim(ncid, 'ple',   lm+1,      ple_did)
      status = nf_def_dim(ncid, 'NTYPE',ntype, ntype_did)
c-----------------------------------------------------------------------
c define some supplemental variables not currently present in acc-file
c-----------------------------------------------------------------------
c longitudes at a-grid grid box centers
      status = nf_def_var (ncid, 'longitude', nf_real, 1, lon_did,varid)
c latitudes at a-grid grid box centers
      status = nf_def_var (ncid, 'latitude', nf_real, 1, lat_did, varid)
c a-grid grid box areas
      status = nf_def_var (ncid, 'area', nf_real, 1, lat_did, varid)
c longitudes at b-grid grid box centers
      status = nf_def_var (ncid, 'lonb', nf_real, 1, lonb_did, varid)
c latitudes at b-grid grid box centers
      status = nf_def_var (ncid, 'latb', nf_real, 1, latb_did, varid)
c b-grid grid box areas
      status = nf_def_var (ncid, 'areab', nf_real, 1, latb_did, varid)
c b-grid grid box east-west length
      status = nf_def_var (ncid, 'dxv', nf_real, 1, latb_did, varid)
c sigmas
      status = nf_def_var(ncid, 'sig', nf_real, 1, sig_did, varid)
c edge sigmas
      status = nf_def_var(ncid, 'sige', nf_real, 1, sige_did, varid)
c global mean layer midpoint pressures
      status = nf_def_var(ncid, 'p', nf_real, 1, plm_did, varid)
c global mean layer edge pressures
      status = nf_def_var(ncid, 'ple', nf_real, 1, ple_did, varid)

c-----------------------------------------------------------------------
c idacc
c-----------------------------------------------------------------------
      status = nf_def_dim(ncid, 'NDACC', 12,  dimids(1))
      status = nf_def_var(ncid, 'IDACC', nf_int, 1, dimids, varid)
c-----------------------------------------------------------------------
c surface type names
c-----------------------------------------------------------------------
      status=nf_def_dim(ncid,'STYPE_CLEN',len(stype_names(1)),dimids(1))
      dimids(2) = ntype_did
      status = nf_def_var(ncid, 'STYPE_NAMES', nf_char, 2,dimids,varid)
c-----------------------------------------------------------------------
c keynr
c-----------------------------------------------------------------------
      status = nf_def_dim(ncid, 'NKEYNR',42, dimids(1))
      status = nf_def_dim(ncid, 'KEYNR2',50, dimids(2))
      status = nf_def_var(ncid, 'KEYNR', nf_int, 2,dimids,varid)
c-----------------------------------------------------------------------
c tsfrez
c-----------------------------------------------------------------------
      dimids(1) = lon_did
      dimids(2) = lat_did
      call ncdefarr(tsf_name,tsf_lname,tsf_units,
     &     ncid,ktsf,nf_real,2,dimids)
c-----------------------------------------------------------------------
c aj
c-----------------------------------------------------------------------
      dimids(1)=lat_did
      dimids(1)=ntype_did
      call ncdefarr(aj_name,aj_lname,aj_units,
     &     ncid,kaj,nf_real,2,dimids)
      call ncput_idacc(aj_name,aj_ia,ncid,kaj)
      status = nf_def_dim(ncid, 'reg',nreg, dimids(1))
      call ncdefarr(dj_name,aj_lname,aj_units,
     &     ncid,kaj,nf_real,1,dimids)
      call ncput_idacc(dj_name,aj_ia,ncid,kaj)

c-----------------------------------------------------------------------
c apj
c-----------------------------------------------------------------------
      dimids(1)=lat_did
      call ncdefarr(apj_name,apj_lname,apj_units,
     &     ncid,kapj,nf_real,1,dimids)
c-----------------------------------------------------------------------
c ajl
c-----------------------------------------------------------------------
      dimids(1)=lat_did
      dimids(2)=sig_did
      call ncdefarr(ajl_name,ajl_lname,ajl_units,
     &     ncid,kajl,nf_real,2,dimids)
c-----------------------------------------------------------------------
c asjl
c-----------------------------------------------------------------------
      status = nf_def_dim(ncid, 'str_level', lm_req,  str_levid)
      dimids(1)=lat_did
      dimids(2)=str_levid
      call ncdefarr(asjl_name,asjl_lname,asjl_units,
     &     ncid,kasjl,nf_real,2,dimids)

c-----------------------------------------------------------------------
c aij
c-----------------------------------------------------------------------
      dimids(1)=lon_did
      dimids(2)=lat_did
      call ncdefarr(aij_name,aij_lname,aij_units,
     &     ncid,kaij,nf_real,2,dimids)
      call ncput_idacc(aij_name,aij_ia,ncid,kaij)
c-----------------------------------------------------------------------
c ail
c-----------------------------------------------------------------------
      dimids(1)=lon_did
      dimids(2)=sig_did
      call ncdefarr(ail_name,ail_lname,ail_units,
     &     ncid,kail,nf_real,2,dimids)
c-----------------------------------------------------------------------
c aijg
c-----------------------------------------------------------------------
      dimids(1)=lon_did
      dimids(2)=lat_did
      call ncdefarr(aijg_name,aijg_lname,aijg_units,
     &     ncid,kaijg,nf_real,2,dimids)
c-----------------------------------------------------------------------
c energy
c-----------------------------------------------------------------------
      status = nf_def_dim(ncid, 'ENERGY1', 20, dimids(1))
      status = nf_def_dim(ncid, 'ENERGY2',100, dimids(2))
      status = nf_def_var(ncid, 'ENERGY' ,nf_real,2,dimids,varid)
c-----------------------------------------------------------------------
c consrv
c-----------------------------------------------------------------------
      dimids(1)=lat_did
      call ncdefarr(consrv_name,consrv_lname,consrv_units,
     &     ncid,kcon,nf_real,1,dimids)
c-----------------------------------------------------------------------
c speca
c-----------------------------------------------------------------------
      status = nf_def_dim(ncid, 'IM_BY2_PLUS1',im/2+1, dimids(1))
      status = nf_def_dim(ncid, 'SP_QUANT'    ,    20, dimids(2))
      status = nf_def_dim(ncid, 'SP_NSPHER'   ,     8, dimids(3))
      status = nf_def_var(ncid, 'SPECA'       ,nf_real,3,dimids,varid)
c-----------------------------------------------------------------------
c atpe
c-----------------------------------------------------------------------
      status = nf_def_dim(ncid, 'ATPE1', 8, dimids(1)) ! m=1->8 in d5
      status = nf_def_dim(ncid, 'NHEMI', 2, dimids(2))
      status = nf_def_var(ncid, 'ATPE', nf_real, 2, dimids, varid)
c-----------------------------------------------------------------------
c adaily
c-----------------------------------------------------------------------
      status = nf_def_dim(ncid, 'HR_IN_DAY', 24, dimids(1))
      status = nf_def_dim(ncid, 'NDLYVAR',   63, dimids(2))
      status = nf_def_dim(ncid, 'NDLYPT', 4, dimids(3))
      status = nf_def_var(ncid, 'ADAILY',nf_real,3,dimids,varid)
      status = nf_put_att_text(ncid,varid,'NAMD6',16,namd6)
      status = nf_put_att_int(ncid,varid,'IJD6',nf_int,8,jc(81))
c-----------------------------------------------------------------------
c wave power
c-----------------------------------------------------------------------
      status = nf_def_dim(ncid, 'RE_AND_IM', RE_AND_IM, dimids(1))
      status = nf_def_dim(ncid, 'Max12hr_sequ', Max12hr_sequ, dimids(2))
      status = nf_def_dim(ncid, 'NWAV_DAG', NWAV_DAG, dimids(3))
      call ncdefarr(wave_name,wave_lname,wave_units,
     &     ncid,kwp,nf_real,3,dimids)
c-----------------------------------------------------------------------
c ajk
c-----------------------------------------------------------------------
      dimids(1)=latb_did
      dimids(2)=sig_did
      call ncdefarr(ajk_name,ajk_lname,ajk_units,
     &     ncid,kajk,nf_real,2,dimids)
c-----------------------------------------------------------------------
c aijk (need to define the b-grid lon_did,lat_did)
c-----------------------------------------------------------------------
      dimids(1)=lonb_did
      dimids(2)=latb_did
      dimids(3)=sig_did
      call ncdefarr(aijk_name,aijk_lname,aijk_units,
     &     ncid,kaijk,nf_real,3,dimids)
c-----------------------------------------------------------------------
c aijl
c-----------------------------------------------------------------------
      dimids(1)=lon_did
      dimids(2)=lat_did
      dimids(3)=sig_did
      call ncdefarr(aijl_name,aijl_lname,aijl_units,
     &     ncid,kaijl,nf_real,3,dimids)
c-----------------------------------------------------------------------
c ajlsp_[dse,lht,mom]
c-----------------------------------------------------------------------
      dimids(1)=lat_did
      dimids(2)=sig_did
      status = nf_def_dim(ncid,'ZERO_TO_NWAV_DAG',1+NWAV_DAG,dimids(3))
      call ncdefarr(ajlsp_name,ajlsp_lname,ajlsp_units,
     &     ncid,kajlsp,nf_real,3,dimids)

c-----------------------------------------------------------------------
c end netcdf definitions
c-----------------------------------------------------------------------

      status = nf_enddef (ncid)

c-----------------------------------------------------------------------
c write data to netcdf file
c-----------------------------------------------------------------------
      var_name='longitude'; call ncwrtdbl1(var_name,ncid,lonab_deg(1,1))
      var_name='latitude';    call ncwrtdbl1(var_name,ncid,lat_deg)
      var_name='area';   call ncwrtdbl1(var_name,ncid,dxyp)
      var_name='lonb';   call ncwrtdbl1(var_name,ncid,lonab_deg(1,2))
      var_name='latb';   call ncwrtdbl1(var_name,ncid,latb_deg)
      var_name='areab';  call ncwrtdbl1(var_name,ncid,dxyv)
      var_name='dxv';  call ncwrtdbl1(var_name,ncid,dxv)
      var_name='sig';    call ncwrtdbl1(var_name,ncid,sig)
      var_name='sige';   call ncwrtdbl1(var_name,ncid,sige)
      var_name='p';   call ncwrtdbl1(var_name,ncid,plm_gcm)
      var_name='ple';   call ncwrtdbl1(var_name,ncid,ple_gcm)
      var_name='IDACC';  call ncwrtint1(var_name,ncid,idacc)
      var_name='STYPE_NAMES';  call ncwrtchar(var_name,ncid,stype_names)
      var_name='KEYNR';  call ncwrtint1(var_name,ncid,keynr)
      call ncwrtdbl(tsf_name,ncid,im*jm,ktsf,tsfrez)
      do k=1,kaj
         aj_tmp=aj(:,k,:)
         call ncwrtdbl1(aj_name(k),ncid,aj_tmp)
      enddo
      call ncwrtdbl(dj_name,ncid,nreg,kaj,areg)
      call ncwrtdbl(apj_name,ncid,jm,kapj,apj)
      call ncwrtdbl(ajl_name,ncid,jm*lm,kajl,ajl)
      call ncwrtdbl(asjl_name,ncid,jm*lm_req,kasjl,asjl)
      call ncwrtdbl(aij_name,ncid,im*jm,kaij,aij)
      call ncwrtdbl(ail_name,ncid,im*lm,kail,ail)
      call ncwrtdbl(aijg_name,ncid,im*jm,kaijg,aijg)
      var_name='ENERGY'; call ncwrtdbl1(var_name,ncid,energy)
      call ncwrtdbl(consrv_name,ncid,jm,kcon,consrv)
      var_name='SPECA';  call ncwrtdbl1(var_name,ncid,speca)
      var_name='ATPE';   call ncwrtdbl1(var_name,ncid,atpe)
      var_name='ADAILY'; call ncwrtdbl1(var_name,ncid,adaily)
      call ncwrtdbl(wave_name,ncid,re_and_im*Max12hr_sequ*NWAV_DAG,
     &     kwp,wave)
      call ncwrtdbl(ajk_name,ncid,jm*lm,kajk,ajk)
      call ncwrtdbl(aijk_name,ncid,im*jm*lm,kaijk,aijk)
      call ncwrtdbl(aijl_name,ncid,im*jm*lm,kaijl,aijl)
      call ncwrtdbl(ajlsp_name,ncid,jm*lm*(nwav_dag+1),kajlsp,ajlsp)
c-----------------------------------------------------------------------
c close netcdf file
c-----------------------------------------------------------------------
      status = nf_close (ncid)
      return
      end subroutine write_nc_acc

      subroutine ncdefarr(arr_name,arr_lname,arr_units,
     &     ncid,arr_size,nf_type,ndims,dimids)
      implicit none
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
      integer :: ncid,arr_size,nf_type,ndims
      integer, dimension(ndims) :: dimids
      integer :: status
      character(len=20), dimension(arr_size) :: arr_name,arr_units
      character(len=80), dimension(arr_size) :: arr_lname
      integer :: n,varid
      do n=1,arr_size
         status = nf_def_var(
     &        ncid,trim(arr_name(n)),nf_type,ndims,dimids,varid)
         status = nf_put_att_text(ncid,varid,'long_name',
     &        len_trim(arr_lname(n)),trim(arr_lname(n)))
         status = nf_put_att_text(ncid,varid,'units',
     &        len_trim(arr_units(n)),trim(arr_units(n)))
      enddo
      return
      end subroutine ncdefarr

      subroutine ncput_idacc(arr_name,arr_idacc,ncid,arr_size)
      implicit none
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
      integer :: ncid,arr_size
      character(len=20), dimension(arr_size) :: arr_name
      integer, dimension(arr_size) :: arr_idacc
      integer :: status,n,varid
      do n=1,arr_size
         status = nf_inq_varid(ncid,trim(arr_name(n)),varid)
         status = nf_put_att_int(ncid,varid,'idacc',nf_int,1,
     &        arr_idacc(n))
      enddo
      return
      end subroutine ncput_idacc

      subroutine ncwrtdbl(arr_name,ncid,inner_size,arr_size,arr)
      implicit none
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
      integer :: ncid,inner_size,arr_size
      integer :: status
      character(len=20), dimension(arr_size) :: arr_name
      double precision, dimension(inner_size,arr_size) :: arr
      integer :: n,varid
      do n=1,arr_size
         status = nf_inq_varid(ncid,trim(arr_name(n)),varid)
         status = nf_put_var_double(ncid,varid,arr(1,n))
      enddo
      return
      end subroutine ncwrtdbl

      subroutine ncgetdbl(arr_name,ncid,inner_size,arr_size,arr)
      implicit none
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
      integer :: ncid,inner_size,arr_size
      integer :: status
      character(len=20), dimension(arr_size) :: arr_name
      double precision, dimension(inner_size,arr_size) :: arr
      integer :: n,varid
      do n=1,arr_size
         status = nf_inq_varid(ncid,trim(arr_name(n)),varid)
         status = nf_get_var_double(ncid,varid,arr(1,n))
      enddo
      return
      end subroutine ncgetdbl

      subroutine ncwrtdbl1(var_name,ncid,var)
      implicit none
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
      integer :: ncid
      integer :: status,varid
      character(len=20) :: var_name
      double precision :: var
      status = nf_inq_varid(ncid,trim(var_name),varid)
      status = nf_put_var_double(ncid,varid,var)
      return
      end subroutine ncwrtdbl1

      subroutine ncgetdbl1(var_name,ncid,var)
      implicit none
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
      integer :: ncid
      integer :: status,varid
      character(len=20) :: var_name
      double precision :: var
      status = nf_inq_varid(ncid,trim(var_name),varid)
      status = nf_get_var_double(ncid,varid,var)
      return
      end subroutine ncgetdbl1

      subroutine ncwrtint1(var_name,ncid,var)
      implicit none
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
      integer :: ncid
      integer :: status,varid
      character(len=20) :: var_name
      integer :: var
      status = nf_inq_varid(ncid,trim(var_name),varid)
      status = nf_put_var_int(ncid,varid,var)
      return
      end subroutine ncwrtint1

      subroutine ncgetint1(var_name,ncid,var)
      implicit none
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
      integer :: ncid
      integer :: status,varid
      character(len=20) :: var_name
      integer :: var
      status = nf_inq_varid(ncid,trim(var_name),varid)
      status = nf_get_var_int(ncid,varid,var)
      return
      end subroutine ncgetint1

      subroutine ncwrtchar(var_name,ncid,var)
      implicit none
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
      integer :: ncid
      integer :: status,varid
      character(len=20) :: var_name
      character :: var
      status = nf_inq_varid(ncid,trim(var_name),varid)
      status = nf_put_var_text(ncid,varid,var)
      return
      end subroutine ncwrtchar

      subroutine ncgetchar(var_name,ncid,var)
      implicit none
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
      integer :: ncid
      integer :: status,varid
      character(len=20) :: var_name
      character :: var
      status = nf_inq_varid(ncid,trim(var_name),varid)
      status = nf_get_var_text(ncid,varid,var)
      return
      end subroutine ncgetchar

      subroutine ncwrt_iparm(iparm_name,ncid,iparm,kjc)
      implicit none
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
      character(len=20), dimension(kjc) :: iparm_name
      integer, dimension(kjc) :: iparm
      integer :: ncid,kjc
      integer :: k,status
      do k=1,kjc
         if(iparm_name(k).eq.'nowrite') cycle
         status = nf_put_att_int(ncid,nf_global,trim(iparm_name(k)),
     &                           nf_int,1,iparm(k))
      enddo
      return
      end subroutine ncwrt_iparm

      subroutine ncget_iparm(iparm_name,ncid,iparm,kjc)
      implicit none
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
      character(len=20), dimension(kjc) :: iparm_name
      integer, dimension(kjc) :: iparm
      integer :: ncid,kjc
      integer :: k,status
      do k=1,kjc
         if(iparm_name(k).eq.'nowrite') cycle
         status = nf_get_att_int(
     &        ncid,nf_global,trim(iparm_name(k)),iparm(k))
      enddo
      return
      end subroutine ncget_iparm

      subroutine ncwrt_dparm(dparm_name,ncid,dparm,krc)
      implicit none
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
      character(len=20), dimension(krc) :: dparm_name
      double precision, dimension(krc) :: dparm
      integer :: ncid,krc
      integer :: k,status
      do k=1,krc
         if(dparm_name(k).eq.'nowrite') cycle
         status = nf_put_att_double(ncid,nf_global,trim(dparm_name(k)),
     &                           nf_real,1,dparm(k))
      enddo
      return
      end subroutine ncwrt_dparm

      subroutine ncget_dparm(dparm_name,ncid,dparm,krc)
      implicit none
      include '/usr/local/netcdf-3.4/include/netcdf.inc'
      character(len=20), dimension(krc) :: dparm_name
      double precision, dimension(krc) :: dparm
      integer :: ncid,krc
      integer :: k,status
      do k=1,krc
         if(dparm_name(k).eq.'nowrite') cycle
         status=nf_get_att_double(
     &        ncid,nf_global,trim(dparm_name(k)),dparm(k))
      enddo
      return
      end subroutine ncget_dparm
