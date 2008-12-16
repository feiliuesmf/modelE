      subroutine regrid_input(dd2d)
      use regrid_com
      use dd2d_utils
      implicit none
      type (dd2d_grid), intent(in) :: dd2d
      type (x_2gridsroot) :: xll2cs

      call init_regrid_root(xll2cs,72,46,1,48,48,6)
c      call init_regrid(xll2cs,360,180,1,90,90,6)

ccc   regrid boundary condition files 
      write(*,*) "IN REGRID INPUT"
      if (AM_I_ROOT()) then
c         call regridTOPO(xll2cs)
c         call regridOSST(xll2cs)
c         call regridSICE(xll2cs)
c         call regridCDN(xll2cs)
c         call regridVEG(xll2cs)
c         call regridRVR(xll2cs) ! empty for the moment
c         call regridCROPS(xll2cs)
c         call regridTOPINDEX(xll2cs)
c         call regridSOIL(xll2cs)
ccc   Then regrid Initial Condition 
      endif
         call regridGIC(xll2cs,dd2d)
c         call regridAIC(xll2cs)

      end subroutine regrid_input
c     *


      subroutine regridTOPO(x2grids)
c
c     for 1x1 resolution : Jeff posted Z1X1N and Z1X1N_MODELE on Discover
c
      USE FILEMANAGER, only : openunit,closeunit
      use regrid_com
      implicit none
      include 'netcdf.inc'
      character*80 :: TITLE,name,ncfile
      type (x_2gridsroot), intent(inout) :: x2grids
      real*8, allocatable :: ttargglob(:,:,:,:),ones(:,:,:)
      integer :: iu_TOPO,i,j,k,irec,iunit,imt,jmt,ntt,maxrec,status,
     &     vid,fid

      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "imt jmt ntt",imt,jmt,ntt

      allocate( ttargglob(imt,jmt,ntt,nrecmax),ones(imt,jmt,ntt) )

      name="Z72X46N.cor4_nocasp"

      write(*,*) name

      open( iu_TOPO, FILE=name,FORM='unformatted', STATUS='old')

      write(*,*) "iu_TOPO",iu_TOPO

      call read_regrid_4D_1R(x2grids,iu_TOPO,ttargglob,maxrec)

c
c     CONSISTENCY CHECKS: 1) FOCEAN+FLAKE+FGRND+FGICE=1 
c                         2) IF FOCEAN(i,j) > 0 set FGRND=FGRND+FLAKE, FLAKE=0

      ones(:,:,:)=ttargglob(:,:,:,1)
     &     +ttargglob(:,:,:,2)
     &     +ttargglob(:,:,:,3)
     &     +ttargglob(:,:,:,4)
      
      if (any( abs(ones(1:imt,1:jmt,1:ntt) - 1.0) .gt. 1.d-6 ) ) then 
         write(*,*) "WARNING FOCEAN+FLAKE+FGRND+FGICE=1 BROKEN"
      endif

      do k=1,ntt
         do j=1,jmt
            do i=1,imt
               if ( (ttargglob(i,j,k,1) .gt. 1.e-6) .and.
     &              (ttargglob(i,j,k,2) .gt. 1.e-6 ) ) then
                  write(*,*) "FGRND=FGRND+FLAKE, FLAKE=0"
                  ttargglob(i,j,k,3)=ttargglob(i,j,k,2)
     &                 +ttargglob(i,j,k,3)
                  ttargglob(i,j,k,2)=0.0
               endif
            enddo
         enddo
      enddo

      close(iu_TOPO)
 
      ncfile="topo6tiles.nc"

      status = nf_open(trim(ncfile),nf_write,fid)
      if (status .ne. NF_NOERR) write(*,*) "UNABLE TO OPEN FILE"
      status = nf_inq_varid(fid,'zatmo',vid)
      write(*,*) NF_STRERROR(status)
      status = nf_put_var_double(fid,vid,ttargglob(:,:,:,1))
      write(*,*) "STATUS",NF_STRERROR(status),"<<"
      status = nf_close(fid)

      ncfile="topo.nc"

      status = nf_open(trim(ncfile),nf_write,fid)
      if (status .ne. NF_NOERR) write(*,*) "UNABLE TO OPEN FILE"
      status = nf_inq_varid(fid,'zatmo',vid)
      write(*,*) NF_STRERROR(status)
      status = nf_put_var_double(fid,vid,ttargglob(:,:,2,1))
      write(*,*) "STATUS",NF_STRERROR(status),"<<"
      status = nf_close(fid)

      deallocate(ttargglob,ones)
      
      end subroutine regridTOPO 
c*



      subroutine regridOSST(x2grids)
c     for 1x1 resolution :  Gary on Athena clima1/OBS/AMIP/1x1
c     Also SST1x1_HadISST from Hadley (Jeff Jonas) - check if this  
c     has been interpolated from lower resolution 

      USE FILEMANAGER, only : openunit,closeunit
      use regrid_com
      implicit none
      type (x_2gridsroot), intent(in) :: x2grids
      character*80 :: TITLE,name
      real*4 OSTmean(x2grids%imsource,x2grids%jmsource),
     &     OSTend(x2grids%imsource,x2grids%jmsource)
      integer iu_OSST,i,j,k

c      call openunit ("OSST",iu_OSST,.true.,.true.)     

      name="OST4X5.B.1993-2002avg.Hadl1.1"
      
      open(iu_OSST, FILE=name,FORM='unformatted', STATUS='old')

      call read_regrid_write_4D_2R(x2grids,name,iu_OSST)

      close(iu_OSST)

      end subroutine regridOSST
c


      subroutine regridSICE(x2grids)
c     for 1x1 resolution: Gary on Athena clima1/OBS/AMIP/1x1
c     Also ICE_1x1_HadISST from Hadley (Jeff Jonas) - check if this  
c     has been interpolated from lower resolution 
c     /u/cmrun/SICE4X5.B.1993-2002avg.Hadl1.1

      USE FILEMANAGER, only : openunit,closeunit
      use regrid_com
      implicit none
      type (x_2gridsroot), intent(in) :: x2grids
      character*80 TITLE,name
      real*4 SICEmean(x2grids%imsource,x2grids%jmsource),
     &     SICEend(x2grids%imsource,x2grids%jmsource)
      integer iu_SICE,i,j,k

c     call openunit ("SICE",iu_SICE,.true.,.true.) 

      name="SICE4X5.B.1993-2002avg.Hadl1.1"
      open (iu_SICE, FILE=name,FORM='unformatted', STATUS='old')
      read(unit=iu_SICE) TITLE 
      write(*,*) TITLE

      call read_regrid_write_4D_2R(x2grids,name,iu_SICE)
      
      close(iu_SICE)
      
      end subroutine regridSICE



      subroutine regridCDN(x2grids)
c     for 1x1 resolution: Jeff uses CDN=AL30RL360X180N.rep
      USE FILEMANAGER, only : openunit,closeunit
      use regrid_com
      implicit none
      type (x_2gridsroot), intent(in) :: x2grids
      character*80 TITLE,name
      integer iu_CDN

      name="CD4X500S.ext"
      open (iu_CDN, FILE=name,FORM='unformatted', STATUS='old')

      call read_regrid_write_4D_1R(x2grids,name,iu_CDN)
      close(iu_CDN)
      
      end subroutine regridCDN
c*



      subroutine regridVEG(x2grids)
c     for 1x1 resolution: Jeff uses VEG=V360X180_no_crops.rep
c	It is identical to V144X90_no_crops.ext (144X90 data 
c	was just transfered to 360X180 grid without any change)
      use regrid_com
      implicit none
      type (x_2gridsroot), intent(in) :: x2grids
      character*80 TITLE,name
      integer iu_VEG
	
      name="V72X46.1.cor2_no_crops.ext"

      open(iu_VEG,FILE=name,FORM='unformatted', STATUS='old')

      call read_regrid_write_4D_1R_rmax(x2grids,name,iu_VEG,10)
           
      close(iu_VEG)
      
      end subroutine regridVEG
c*


      subroutine regridRVR(x2grids)
c	empty for the moment
      USE FILEMANAGER, only : openunit,closeunit
      use regrid_com
      implicit none
      type (x_2gridsroot), intent(inout) :: x2grids
      character*80 TITLE,name
      integer iu_RVR

ccc   EMPTY FOR THE MOMENT
   
      end subroutine regridRVR
c*


      subroutine regridCROPS(x2grids)
      use regrid_com
      implicit none
      type (x_2gridsroot), intent(inout) :: x2grids
      character*80 name
      integer iu_CROPS

      name="CROPS_72X46N.cor4.ext"

      open(iu_CROPS,FILE=name,FORM='unformatted', STATUS='old')

      call read_regrid_write_4D_1R(x2grids,name,iu_CROPS)

      close(iu_CROPS)

      end subroutine regridCROPS
c*



      subroutine regridTOPINDEX(x2grids)
c     top_index_360x180.ij is not "extended". It has "-1" in ocean 
c     cells (meaning missing data), so you should be careful not to 
c     include those cells in the regridding computations. 
c     Alternatively you can use top_index_360x180.ij.rep which is 
c     direct transfer from top_index_144x90.ij.ext 

      USE FILEMANAGER, only : openunit,closeunit
      use regrid_com
      implicit none
      type (x_2gridsroot), intent(in) :: x2grids
      character*80 TITLE,name
      integer iu_TOP_INDEX
     
      name="top_index_72x46.ij.ext"

      open(iu_TOP_INDEX,FILE=name,FORM='unformatted', STATUS='old')
      
      call read_regrid_write_4D_1R(x2grids,name,iu_TOP_INDEX)
      
      close(iu_TOP_INDEX)
      
      end subroutine regridTOPINDEX
c*


      
      subroutine regridSOIL(x2grids)
c
c     for 1x1 resolution: Jeff uses SOIL=S360X180_0098M.rep
c
      use regrid_com
      implicit none
      type (x_2gridsroot), intent(in) :: x2grids
      real*4, allocatable :: dz(:,:,:),ftext(:,:,:,:),
     &     ftextk(:,:,:,:),sl(:,:)
      real*4, allocatable :: dzout(:,:,:,:),ftextout(:,:,:,:,:),
     &     ftextkout(:,:,:,:,:),slout(:,:,:)
      real*8, allocatable :: tsource(:,:,:)
      real*8, allocatable :: ttargglob(:,:,:)
      character*80 TITLE,name,outunformat
      integer iu_SOIL,iuout,ims,jms,nts,imt,jmt,ntt,i,j,k,l

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "ims,jms,nts,imt,jmt,ntt",ims,jms,nts,imt,jmt,ntt

      name="S4X50093.ext"

      open(iu_SOIL,FILE=name,FORM='unformatted', STATUS='old')

      allocate (dz(ims,jms,6),ftext(ims,jms,6,5),
     &     ftextk(ims,jms,6,5),sl(ims,jms),
     &     dzout(imt,jmt,6,ntt),ftextout(imt,jmt,6,5,ntt),
     &     ftextkout(imt,jmt,6,5,ntt),
     &     slout(imt,jmt,ntt))

      allocate (tsource(ims,jms,nts),
     &     ttargglob(imt,jmt,ntt) )

      read(iu_SOIL) dz,ftext,ftextk,sl

      close(iu_SOIL)
      
      outunformat=trim(name)//".CS"
      
      write(*,*) outunformat

      write(*,*) dz
      iuout=20
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")

      do k=1,6
         do j=1,jms
            do i=1,ims
               tsource(i,j,1)=dz(i,j,k)
            enddo
         enddo
         call root_regrid(x2grids,tsource,ttargglob)
         write(*,*) "k=",k
         dzout(:,:,k,:)=ttargglob(:,:,:)
      enddo

      do k=1,6
         do l=1,5
            tsource(:,:,1)=ftext(:,:,k,l)
            call root_regrid(x2grids,tsource,ttargglob)
            ftextout(:,:,k,l,:)=ttargglob(:,:,:)
         enddo
      enddo


      do k=1,6
         do l=1,5
            tsource(:,:,1)=ftextk(:,:,k,l)
            call root_regrid(x2grids,tsource,ttargglob)
            ftextkout(:,:,k,l,:)=ttargglob(:,:,:)
         enddo
      enddo

      tsource(:,:,1)=sl(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      slout(:,:,:)=ttargglob(:,:,:)
      
      write(unit=iuout) dzout,ftextout,ftextkout,slout

      close(iuout)
      
      end subroutine regridSOIL
c*
      
      


      subroutine regridGIC(x2grids,dd2d)
c
c     for 1x1 resolution: Jeff uses GIC=GIC.360X180.DEC01.1.rep
c
      USE FILEMANAGER, only : openunit,closeunit
      use DOMAIN_DECOMP,only : am_i_root
      use regrid_com
      use dd2d_utils
      use pario, only : defvar,write_data

      implicit none
      include 'netcdf.inc'

      type (x_2gridsroot), intent(in) :: x2grids
      type (dd2d_grid), intent(in) :: dd2d
      real*8, allocatable :: Tocn(:,:,:),MixLD(:,:)        ! OCN01
      real*8, allocatable :: F(:,:),H(:,:,:),snw(:,:),msi(:,:),
     &     ssi(:,:,:),pond_melt(:,:)
      logical, allocatable :: flag_dsws(:,:)               ! SICE02
      real*8, allocatable :: snowe(:,:),Te(:,:),WTRe(:,:),ICEe(:,:),
     &     SNOage(:,:,:),evmax(:,:),fsat(:,:),gq(:,:)      ! EARTH01
      real*8, allocatable :: Wb(:,:,:),Wv(:,:,:),HTb(:,:,:),
     &     HTv(:,:,:),SNWbv(:,:,:)                         ! SOILS02
      real*8, allocatable :: SNOW(:,:),T(:,:,:)            ! GLAIC01

      real*8, allocatable :: Tocn_out(:,:,:,:),MixLD_out(:,:,:)   ! OCN01
      real*8, allocatable :: F_out(:,:,:),H_out(:,:,:,:),
     &     snw_out(:,:,:),msi_out(:,:,:),ssi_out(:,:,:,:),
     &     pond_melt_out(:,:,:)
      integer, allocatable :: flag_dsws_out(:,:,:)                ! SICE02
      real*8, allocatable :: snowe_out(:,:,:),Te_out(:,:,:),
     &     WTRe_out(:,:,:), ICEe_out(:,:,:),SNOage_out(:,:,:,:),
     &     evmax_out(:,:,:),fsat_out(:,:,:),gq_out(:,:,:)         ! EARTH01
      real*8, allocatable :: Wb_out(:,:,:,:),Wv_out(:,:,:,:),
     &     HTb_out(:,:,:,:),HTv_out(:,:,:,:),SNWbv_out(:,:,:,:)   ! SOILS02
      real*8, allocatable :: SNOW_out(:,:,:),T_out(:,:,:,:)       ! GLAIC01

      real*8, allocatable :: tsource(:,:,:)
      real*8, allocatable :: ttargglob(:,:,:)
      character*80 TITLEOCN01,TITLESICE02,TITLEEARTH01,TITLESOILS02,
     &     TITLEGLAIC01,name,outunformat,outnc
      integer :: iu_GIC,iuout,ims,jms,nts,imt,jmt,ntt
      integer :: i,j,k,l,fid,status,ntiles,im,jm,d2,d3
     
      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "ims,jms,nts,imt,jmt,ntt",ims,jms,nts,imt,jmt,ntt

      name="GIC.E046D3M20A.1DEC1955.ext"

      open(iu_GIC,FILE=name,FORM='unformatted', STATUS='old')

      allocate (Tocn(3,ims,jms),MixLD(ims,jms),
     &     F(ims,jms),H(4,ims,jms),snw(ims,jms),msi(ims,jms),
     &     ssi(4,ims,jms),pond_melt(ims,jms),flag_dsws(ims,jms),
     &     snowe(ims,jms),Te(ims,jms),WTRe(ims,jms),ICEe(ims,jms), 
     &     SNOage(3,ims,jms),evmax(ims,jms),fsat(ims,jms),gq(ims,jms),
     &     Wb(6,ims,jms),Wv(7,ims,jms),HTb(7,ims,jms),HTv(7,ims,jms), 
     &     SNWbv(2,ims,jms),
     &     SNOW(ims,jms),T(2,ims,jms) )

      allocate (Tocn_out(3,imt,jmt,ntt),MixLD_out(imt,jmt,ntt),
     &     F_out(imt,jmt,ntt),H_out(4,imt,jmt,ntt),
     &     snw_out(imt,jmt,ntt),msi_out(imt,jmt,ntt),
     &     ssi_out(4,imt,jmt,ntt),pond_melt_out(imt,jmt,ntt),
     &     flag_dsws_out(imt,jmt,ntt),
     &     snowe_out(imt,jmt,ntt),Te_out(imt,jmt,ntt),
     &     WTRe_out(imt,jmt,ntt),ICEe_out(imt,jmt,ntt), 
     &     SNOage_out(3,imt,jmt,ntt),evmax_out(imt,jmt,ntt),
     &     fsat_out(imt,jmt,ntt),gq_out(imt,jmt,ntt),
     &     Wb_out(6,imt,jmt,ntt),Wv_out(7,imt,jmt,ntt),
     &     HTb_out(7,imt,jmt,ntt),HTv_out(7,imt,jmt,ntt), 
     &     SNWbv_out(2,imt,jmt,ntt),
     &     SNOW_out(imt,jmt,ntt),T_out(2,imt,jmt,ntt) )


      allocate (tsource(ims,jms,nts),
     &     ttargglob(imt,jmt,ntt) )

      read(iu_GIC) TITLEOCN01, Tocn,MixLD
      write(*,*) TITLEOCN01
      read(iu_GIC) TITLESICE02, F,H,snw,msi,ssi,pond_melt,flag_dsws
      write(*,*) TITLESICE02
      read(iu_GIC) TITLEEARTH01, snowe,Te,WTRe, ICEe, SNOage,evmax,
     &     fsat,gq
      write(*,*) TITLEEARTH01
      read(iu_GIC) TITLESOILS02, Wb,Wv,HTb,HTv,SNWbv
      write(*,*) TITLESOILS02
      read(iu_GIC) TITLEGLAIC01, SNOW,T
      write(*,*) TITLEGLAIC01

      close(iu_GIC)
      
            
      outunformat=trim(name)//".CS"
      
      write(*,*) outunformat
      
      iuout=20
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")

      do k=1,3
         do j=1,jms
            do i=1,ims
               tsource(i,j,1)=Tocn(k,i,j)
            enddo
         enddo
         call root_regrid(x2grids,tsource,ttargglob)
         write(*,*) "k=",k
         Tocn_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

     
      tsource(:,:,1)=MixLD(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      MixLD_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=F(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      F_out(:,:,:)=ttargglob(:,:,:)

      do k=1,4
         do j=1,jms
            do i=1,ims
               tsource(i,j,1)=H(k,i,j)
            enddo
         enddo
         call root_regrid(x2grids,tsource,ttargglob)
         write(*,*) "k=",k
         H_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

      tsource(:,:,1)=snw(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      snw_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=msi(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      msi_out(:,:,:)=ttargglob(:,:,:)

      do k=1,4
         do j=1,jms
            do i=1,ims
               tsource(i,j,1)=ssi(k,i,j)
            enddo
         enddo
         call root_regrid(x2grids,tsource,ttargglob)
         write(*,*) "k=",k
         ssi_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

      tsource(:,:,1)=pond_melt(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      pond_melt_out(:,:,:)=ttargglob(:,:,:)

      do j=1,jms
         do i=1,ims
            if (flag_dsws(i,j) .eq. .true.) then
               tsource(i,j,1)=1.d0
            else
               tsource(i,j,1)=0.d0
            endif
         enddo
      enddo

      call root_regrid(x2grids,tsource,ttargglob)

c***  Compatibility
      do k=1,6
         do j=1,jms
            do i=1,ims
               if (ttargglob(i,j,k) .ge. 0.5d0) then
                  flag_dsws_out(i,j,k)=1
               else
                  flag_dsws_out(i,j,k)=0
               endif
            enddo
         enddo
      enddo


      tsource(:,:,1)=snowe(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      snowe_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=Te(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      Te_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=WTRe(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      WTRe_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=ICEe(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      ICEe_out(:,:,:)=ttargglob(:,:,:)

      do k=1,3
         do j=1,jms
            do i=1,ims
               tsource(i,j,1)=SNOage(k,i,j)
            enddo
         enddo
         call root_regrid(x2grids,tsource,ttargglob)
         write(*,*) "k=",k
         SNOage_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

      tsource(:,:,1)=evmax(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      evmax_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=fsat(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      fsat_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=gq(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      gq_out(:,:,:)=ttargglob(:,:,:)

      do k=1,6
         do j=1,jms
            do i=1,ims
               tsource(i,j,1)=Wb(k,i,j)
            enddo
         enddo
         call root_regrid(x2grids,tsource,ttargglob)
         write(*,*) "k=",k
         Wb_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

      do k=1,7
         do j=1,jms
            do i=1,ims
               tsource(i,j,1)=Wv(k,i,j)
            enddo
         enddo
         call root_regrid(x2grids,tsource,ttargglob)
         write(*,*) "k=",k
         Wv_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

      do k=1,7
         do j=1,jms
            do i=1,ims
               tsource(i,j,1)=HTb(k,i,j)
            enddo
         enddo
         call root_regrid(x2grids,tsource,ttargglob)
         write(*,*) "k=",k
         HTb_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

      do k=1,7
         do j=1,jms
            do i=1,ims
               tsource(i,j,1)=HTv(k,i,j)
            enddo
         enddo
         call root_regrid(x2grids,tsource,ttargglob)
         write(*,*) "k=",k
         HTv_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

      do k=1,2
         do j=1,jms
            do i=1,ims
               tsource(i,j,1)=SNWbv(k,i,j)
            enddo
         enddo
         call root_regrid(x2grids,tsource,ttargglob)
         write(*,*) "k=",k
         SNWbv_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

      tsource(:,:,1)=SNOW(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      SNOW_out(:,:,:)=ttargglob(:,:,:)

      do k=1,2
         do j=1,jms
            do i=1,ims
               tsource(i,j,1)=T(k,i,j)
            enddo
         enddo
         call root_regrid(x2grids,tsource,ttargglob)
         write(*,*) "k=",k
         T_out(k,:,:,:)=ttargglob(:,:,:)
      enddo


c***  Write Netcdf file
#ifdef TRACERS_WATER
      write(*,*) "STOP TRACERS WATER NOT IMPLEMENTED IN regridinput"
      stop
#endif      
      

      deallocate (Tocn,MixLD,F,H,snw,msi,
     &     ssi,pond_melt,flag_dsws,
     &     snowe,Te,WTRe,ICEe, 
     &     SNOage,evmax,fsat,gq,
     &     Wb,Wv,HTb,HTv, 
     &     SNWbv,SNOW,T )

      deallocate (tsource,ttargglob)     
      
      outnc=trim(name)//"-CS.nc"
      write(*,*) outnc

      if (am_i_root()) then
         status = nf_create(outnc,nf_clobber,fid)
         if (status .ne. NF_NOERR) write(*,*) "UNABLE TO CREATE FILE"
      endif

    
c***  Define OCN variables
      call defvar(dd2d,fid,Tocn_out,'tocean_glob(d3,im,jm,tile)')
      call defvar(dd2d,fid,MixLD_out,'z1o_glob(im,jm,tile)')
c***  Define SICE variables
      call defvar(dd2d,fid,F_out,'rsi_glob(im,jm,tile)')
      call defvar(dd2d,fid,H_out,'hsi_glob(lmi,im,jm,tile)')
      call defvar(dd2d,fid,snw_out,'snowi_glob(im,jm,tile)')
      call defvar(dd2d,fid,msi_out,'msi_glob(im,jm,tile)')
      call defvar(dd2d,fid,ssi_out,'ssi_glob(lmi,im,jm,tile)')
      call defvar(dd2d,fid,pond_melt_out,
     &     'pond_melt_glob(im,jm,tile)')
      call defvar(dd2d,fid,flag_dsws_out,
     &     'flag_dsws_glob(im,jm,tile)')
c***  Define EARTH variables
      call defvar(dd2d,fid,snowe_out,'snowe_glob(im,jm,tile)')
      call defvar(dd2d,fid,Te_out,'tearth_glob(im,jm,tile)')
      call defvar(dd2d,fid,WTRe_out,'wearth_glob(im,jm,tile)')
      call defvar(dd2d,fid,ICEe_out,'aiearth_glob(im,jm,tile)')
      call defvar(dd2d,fid,SNOage_out,'snoage_glob(d3,im,jm,tile)')
      call defvar(dd2d,fid,evmax_out,
     &     'evap_max_ij_glob(im,jm,tile)')
      call defvar(dd2d,fid,fsat_out,'fr_sat_ij_glob(im,jm,tile)')
      call defvar(dd2d,fid,gq_out,'qg_ij_glob(im,jm,tile)')
c***  Define SOIL variables     ! this is the old SOIL02 version, implement the new version
      call defvar(dd2d,fid,Wb_out,'wb_glob(ngm,im,jm,tile)')
      call defvar(dd2d,fid,Wv_out,
     &     'wv_glob(zero_to_ngm,im,jm,tile)')
      call defvar(dd2d,fid,HTb_out,
     &     'htb_glob(zero_to_ngm,im,jm,tile)')
      call defvar(dd2d,fid,HTv_out,
     &     'htv_glob(zero_to_ngm,im,jm,tile)')
      call defvar(dd2d,fid,SNWbv_out,
     &     'snowbv_glob(ls_nfrac,im,jm,tile)')
c***  Define GLAIC variables
      call defvar(dd2d,fid,SNOW_out,'snowli_glob(im,jm,tile)')
      call defvar(dd2d,fid,T_out,'tlandi_glob(d2,im,jm,tile)')
      
      if (am_i_root()) then
         status = nf_enddef(fid)
         if (status .ne. NF_NOERR) write(*,*) "Problem with enddef"
      endif

c***  Write OCN variables
      call write_data(dd2d,fid,'tocean_glob',Tocn_out)
      call write_data(dd2d,fid,'z1o_glob',MixLD_out)
c***  Write SICE variables
      call write_data(dd2d,fid,'rsi_glob',F_out)
      call write_data(dd2d,fid,'hsi_glob',H_out)
      call write_data(dd2d,fid,'snowi_glob',snw_out)
      call write_data(dd2d,fid,'msi_glob',msi_out)
      call write_data(dd2d,fid,'ssi_glob',ssi_out)
      call write_data(dd2d,fid,'pond_melt_glob',pond_melt_out)
      call write_data(dd2d,fid,'flag_dsws_glob',flag_dsws_out)
c***  Write EARTH variables
      call write_data(dd2d,fid,'snowe_glob',snowe_out)
      call write_data(dd2d,fid,'tearth_glob',Te_out)
      call write_data(dd2d,fid,'wearth_glob',WTRe_out)
      call write_data(dd2d,fid,'aiearth_glob',ICEe_out)
      call write_data(dd2d,fid,'snoage_glob',SNOage_out)
      call write_data(dd2d,fid,'evap_max_ij_glob',evmax_out)
      call write_data(dd2d,fid,'fr_sat_ij_glob',fsat_out)
      call write_data(dd2d,fid,'qg_ij_glob',gq_out)
c***  Write SOIL variables     ! this is the old SOIL02 version, implement the new version
      call write_data(dd2d,fid,'wb_glob',Wb_out)
      call write_data(dd2d,fid,'wv_glob',Wv_out)
      call write_data(dd2d,fid,'htb_glob',HTb_out)
      call write_data(dd2d,fid,'htv_glob',HTv_out)
      call write_data(dd2d,fid,'snowbv_glob',SNWbv_out)
c***  Write GLAIC variables
      call write_data(dd2d,fid,'snowli_glob',SNOW_out)
      call write_data(dd2d,fid,'tlandi_glob',T_out)


      deallocate (Tocn_out,MixLD_out,F_out,H_out,
     &     snw_out,msi_out,ssi_out,pond_melt_out,
     &     flag_dsws_out,snowe_out,Te_out,
     &     WTRe_out,ICEe_out,SNOage_out,evmax_out,
     &     fsat_out,gq_out,Wb_out,Wv_out,
     &     HTb_out,HTv_out,SNWbv_out,
     &     SNOW_out,T_out )

      
       if(am_i_root()) status = nf_close(fid)

      end subroutine regridGIC
c*



      subroutine regridAIC(x2grids)
c
c     for 1x1 resolution : Jeff uses AIC=AIC.RES_X40.D771201N.rep
c
      USE FILEMANAGER, only : openunit,closeunit
      use regrid_com
      implicit none
      type (x_2gridsroot), intent(in) :: x2grids
      character*80 TITLE,name
      integer iu_AIC
     
      name="AIC.RES_M20A.D771201"

      open(iu_AIC,FILE=name,FORM='unformatted', STATUS='old')
      
      call read_regrid_write_4D_1R_r8(x2grids,name,iu_AIC)
      
      close(iu_AIC)
      
      end subroutine regridAIC
c*


      subroutine read_recs_1R(tsource,iuin,TITLE,maxrec,im1,jm1,ntl)
      use regrid_com
      implicit none
      integer i,j,k,irec
      integer, intent(in) :: im1,jm1,ntl
      real*4 :: data(im1,jm1,ntl,nrecmax)
      real*8, intent(inout) :: tsource(im1,jm1,ntl,nrecmax)
      integer, intent(in) :: iuin
      character*80, intent(inout) :: TITLE(nrecmax)
      integer, intent(out) :: maxrec
      
      write(*,*) "iuin",iuin
      irec=1

      do
         read(unit=iuin,END=30) TITLE(irec), data(:,:,:,irec)
         write(*,*) "TITLE, irec",TITLE(irec),irec
         tsource(:,:,:,irec)= data(:,:,:,irec)
         irec=irec+1
      enddo

 30   continue

      maxrec=irec-1

      close(iuin)

      end subroutine read_recs_1R
c*


      subroutine read_recs_1R_r8(tsource,iuin,TITLE,maxrec,im1,jm1,
     *     ntl)
      use regrid_com
      implicit none
      integer i,j,k,irec
      integer, intent(in) :: im1,jm1,ntl
      real*8, intent(inout) :: tsource(im1,jm1,ntl,nrecmax)
      integer, intent(in) :: iuin
      character*80, intent(inout) :: TITLE(nrecmax)
      integer, intent(out) :: maxrec
      
      write(*,*) "iuin",iuin
      irec=1

      do
         read(unit=iuin,END=30) TITLE(irec), tsource(:,:,:,irec)
         write(*,*) "TITLE, irec",TITLE(irec),irec
         irec=irec+1
      enddo

 30   continue

      maxrec=irec-1

      close(iuin)

      end subroutine read_recs_1R_r8
c*


      subroutine read_recs_2R(tsource1,tsource2,iuin,
     &     TITLE,maxrec,im1,jm1,ntl)
      use regrid_com
      implicit none
      integer i,j,k,irec
      integer, intent(in) :: im1,jm1,ntl
      real*4 :: data1(im1,jm1,ntl,nrecmax),
     &     data2(im1,jm1,ntl,nrecmax)
      real*8, intent(inout) :: tsource1(im1,jm1,ntl,nrecmax),
     &     tsource2(im1,jm1,ntl,nrecmax)
      integer, intent(in) :: iuin
      character*80, intent(inout) :: TITLE(nrecmax)
      integer, intent(out) :: maxrec
      
      write(*,*) "iuin",iuin
      irec=1

      do
         read(unit=iuin,END=30) TITLE(irec), data1(:,:,:,irec),
     &        data2(:,:,:,irec)

         tsource1(:,:,:,irec)= data1(:,:,:,irec)
         tsource2(:,:,:,irec)= data2(:,:,:,irec)
         irec=irec+1
      enddo

 30   continue

      maxrec=irec-1

      close(iuin)

      end subroutine read_recs_2R
c*


      subroutine read_regrid_write_4D_1R(x2grids,name,iuin)
      use regrid_com
      implicit none
      type(x_2gridsroot), intent(in) :: x2grids
      character*80, intent(in) :: name
      integer, intent(in) :: iuin
      integer :: iuout
      real*8, allocatable :: tsource(:,:,:,:)
      real*8, allocatable :: ttargglob(:,:,:)
      real*4, allocatable :: tout(:,:,:)
      character*80 :: TITLE(nrecmax),outunformat
      integer :: maxrec,irec,ir,ims,jms,nts,imt,jmt,ntt

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "iuin ims,jms,nts,imt,jmt,ntt 4D",iuin,
     &     ims,jms,nts,imt,jmt,ntt
      allocate (tsource(ims,jms,nts,nrecmax),
     &     ttargglob(imt,jmt,ntt),
     &     tout(imt,jmt,ntt))
      tsource(:,:,:,:)=0.0
      
      call read_recs_1R(tsource,iuin,TITLE,
     &        maxrec,ims,jms,nts)
      
      write(*,*) "maxrec",maxrec
      
      outunformat=trim(name)//".CS"
      
      write(*,*) outunformat
      iuout=20
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")

      do ir=1,maxrec
         call root_regrid(x2grids,tsource(:,:,:,ir),ttargglob)
         tout(:,:,:)=ttargglob(:,:,:)
         write(unit=iuout) TITLE(ir),tout(:,:,:)
         write(*,*) "TITLE",TITLE(ir)
      enddo

      close(iuout) 

      deallocate(tsource,ttargglob,tout)

      end subroutine read_regrid_write_4D_1R
c*


      subroutine read_regrid_write_4D_1R_r8(x2grids,name,iuin)
      use regrid_com
      implicit none
      type(x_2gridsroot), intent(in) :: x2grids
      character*80, intent(in) :: name
      integer, intent(in) :: iuin
      integer :: iuout
      real*8, allocatable :: tsource(:,:,:,:)
      real*8, allocatable :: ttargglob(:,:,:)
      character*80 :: TITLE(nrecmax),outunformat
      integer :: maxrec,irec,ir,ims,jms,nts,imt,jmt,ntt

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "iuin ims,jms,nts,imt,jmt,ntt 4D",iuin,
     &     ims,jms,nts,imt,jmt,ntt
      allocate (tsource(ims,jms,nts,nrecmax),
     &     ttargglob(imt,jmt,ntt) )
      tsource(:,:,:,:)=0.0
      
      call read_recs_1R_r8(tsource,iuin,TITLE,
     &        maxrec,ims,jms,nts)
      
      write(*,*) "maxrec",maxrec
      
      outunformat=trim(name)//".CS"
      
      write(*,*) outunformat
      iuout=20
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")

      do ir=1,maxrec
         call root_regrid(x2grids,tsource(:,:,:,ir),ttargglob)
         write(unit=iuout) TITLE(ir),ttargglob(:,:,:)
         write(*,*) "TITLE",TITLE(ir)
      enddo

      close(iuout) 

      deallocate(tsource,ttargglob)

      end subroutine read_regrid_write_4D_1R_r8
c*



      subroutine read_regrid_write_4D_1R_rmax(x2grids,name,iuin,rmax)
      use regrid_com
      implicit none
      type(x_2gridsroot), intent(in) :: x2grids
      character*80, intent(in) :: name
      integer, intent(in) :: iuin
      integer, intent(in):: rmax
      integer :: iuout
      real*8, allocatable :: tsource(:,:,:,:)
      real*8, allocatable :: ttargglob(:,:,:)
      real*4, allocatable :: tout(:,:,:),data(:,:,:,:)
      character*80 :: TITLE(nrecmax),outunformat
      integer :: maxrec,irec,ir,ims,jms,nts,imt,jmt,ntt

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "iuin ims,jms,nts,imt,jmt,ntt 4D",iuin,
     &     ims,jms,nts,imt,jmt,ntt
      allocate (tsource(ims,jms,nts,nrecmax),
     &     ttargglob(imt,jmt,ntt),
     &     tout(imt,jmt,ntt),
     &     data(ims,jms,nts,nrecmax) )
      tsource(:,:,:,:)=0.0
      data(:,:,:,:)=0.0
      

      do irec=1,rmax
         read(unit=iuin) TITLE(irec), data(:,:,:,irec)
         write(*,*) "TITLE, irec",TITLE(irec),irec
         tsource(:,:,:,irec)= data(:,:,:,irec)
      enddo
      maxrec=rmax
            
      write(*,*) "maxrec",maxrec
      
      outunformat=trim(name)//".CS"
      
      write(*,*) outunformat
      iuout=20
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")

      do ir=1,maxrec
         call root_regrid(x2grids,tsource(:,:,:,ir),ttargglob)
         tout(:,:,:)=ttargglob(:,:,:)
         write(unit=iuout) TITLE(ir),tout(:,:,:)
         write(*,*) "TITLE",TITLE(ir)
c         write(*,*) "TOUT>>",tout(:,:,:),"<<<"
      enddo

      close(iuout) 

      deallocate(tsource,ttargglob,tout,data)

      end subroutine read_regrid_write_4D_1R_rmax
c*


      subroutine read_regrid_4D_1R(x2grids,iuin,ttargglob,maxrec)
      use regrid_com
      implicit none
      type(x_2gridsroot), intent(in) :: x2grids
      integer, intent(in) :: iuin
      integer, intent(inout) :: maxrec
      real*8, allocatable :: tsource(:,:,:,:)
      real*8 :: ttargglob(x2grids%imtarget,x2grids%jmtarget,
     &     x2grids%ntilestarget,nrecmax)
      character*80 :: TITLE(nrecmax)
      integer :: ir,ims,jms,nts

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource

      allocate (tsource(ims,jms,nts,nrecmax))
      tsource(:,:,:,:)=0.
      ttargglob(:,:,:,:)=0.

            
      call read_recs_1R(tsource,iuin,TITLE,
     &        maxrec,ims,jms,nts)
      
      write(*,*) "maxrec",maxrec
      
      do ir=1,maxrec
         call root_regrid(x2grids,tsource(:,:,:,ir),
     &        ttargglob(:,:,:,ir))
         write(*,*) "TITLE",TITLE(ir)
      enddo

      deallocate(tsource)

      end subroutine read_regrid_4D_1R
c*



      subroutine read_regrid_write_4D_2R(x2grids,name,iuin)
      use regrid_com
      type(x_2gridsroot), intent(in) :: x2grids
      integer, intent(in) :: iuin
      character*80 name
      integer :: iuout
      real*8, allocatable :: tsource1(:,:,:,:),tsource2(:,:,:,:)
      real*8, allocatable :: ttargglob1(:,:,:),ttargglob2(:,:,:)
      real*4, allocatable :: tout1(:,:,:),tout2(:,:,:)
      character*80 :: TITLE(nrecmax),outunformat
      integer :: maxrec

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "iuin ims,jms,nts,imt,jmt,ntt 4D",iuin,
     &     ims,jms,nts,imt,jmt,ntt
      allocate (tsource1(ims,jms,nts,nrecmax),
     &     tsource2(ims,jms,nts,nrecmax),
     &     ttargglob1(imt,jmt,ntt),
     &     ttargglob2(imt,jmt,ntt),
     &     tout1(imt,jmt,ntt),
     &     tout2(imt,jmt,ntt) )

      tsource1(:,:,:,:)=0.0
      tsource2(:,:,:,:)=0.0
      
      call read_recs_2R(tsource1,tsource2,iuin,TITLE,
     &     maxrec,ims,jms,nts)
      
c      write(*,*) "TITLE RECS",TITLE(:)
      write(*,*) "maxrec",maxrec
      
      outunformat=trim(name)//".CS"
      
      write(*,*) outunformat
      iuout=20
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")

      do ir=1,maxrec
         call root_regrid(x2grids,tsource1(:,:,:,ir),ttargglob1)
         tout1(:,:,:)=ttargglob1(:,:,:)
         call root_regrid(x2grids,tsource2(:,:,:,ir),ttargglob2)
         tout2(:,:,:)=ttargglob2(:,:,:)
         write(unit=iuout) TITLE(ir),tout1(:,:,:),tout2(:,:,:)
         write(*,*) "TITLE",TITLE(ir)
c         write(*,*) "TOUT>>",tout(:,:,:),"<<<"
      enddo
      
      close(iuout) 

      deallocate(tsource1,tsource2,ttargglob1,ttargglob2,tout1,tout2)

      end subroutine read_regrid_write_4D_2R


