      subroutine regrid_input(dd2d)
      use regrid_com
      use dd2d_utils
      implicit none
      type (dd2d_grid), intent(in) :: dd2d
      type (x_2gridsroot) :: xll2cs

      call init_regrid_root(xll2cs,72,46,1,32,32,6)
c      call init_regrid_root(xll144X90C32,72,46,1,32,32,6)
c      call init_regrid_root(xll2cs,360,180,1,32,32,6)
c      call init_regrid_root(xll72X46C32,144,90,1,48,48,6)

ccc   regrid boundary condition files 
      write(*,*) "IN REGRID INPUT"
      if (AM_I_ROOT()) then
c         call regridTOPO(xll72X46C32)
c         call regridOSST(xll72X46C32)
c         call regridSICE(xll72X46C32)
c         call regridCDN(xll72X46C32)
c         call regridVEG(xll72X46C32)
c         call regridRVR(xll72X46C32) ! empty for the moment
c         call regridCROPS(xll72X46C32)
c         call regridTOPINDEX(xll72X46C32)
c         call regridSOIL(xll72X46C32)
      endif
c      call regridGIC(xll2cs,dd2d)
         call regridAIC(xll2cs,dd2d)

      end subroutine regrid_input
c*


      subroutine regridTOPO(x2grids)
c
c     for 1x1 resolution : Jeff posted Z1X1N and Z1X1N_MODELE on Discover
c
      use regrid_com
      implicit none
      include 'netcdf.inc'
      character*80 :: TITLE(nrecmax),name,ncfile
      type (x_2gridsroot), intent(inout) :: x2grids
      real*8, allocatable :: ttargglob(:,:,:,:),ones(:,:,:)
      real*4, allocatable :: ttargr4(:,:,:,:)
      integer :: iu_TOPO,i,j,k,irec,iunit,imt,jmt,ntt,maxrec,status,
     &     vid,fid,ir

      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "imt jmt ntt",imt,jmt,ntt

      allocate( ttargglob(imt,jmt,ntt,nrecmax),
     &     ttargr4(imt,jmt,ntt,nrecmax),
     &     ones(imt,jmt,ntt) )

      name="Z72X46N.cor4_nocasp"

      write(*,*) name

      open( iu_TOPO, FILE=name,FORM='unformatted', STATUS='old')

      write(*,*) "iu_TOPO",iu_TOPO

      call read_regrid_4D_1R(x2grids,iu_TOPO,TITLE,ttargglob,maxrec)

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

      name="Z72X46N.cor4_nocasp.CS"

      write(*,*) name

      open( iu_TOPO, FILE=name,FORM='unformatted', STATUS='unknown')

      ttargr4=ttargglob

      do ir=1,maxrec
         write(unit=iu_TOPO) TITLE(ir), ttargr4(:,:,:,ir)
         write(*,*) TITLE(ir), ttargr4(:,:,:,ir)
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

      use regrid_com
      implicit none
      type (x_2gridsroot), intent(in) :: x2grids
      character*80 :: TITLE,name
      real*4 OSTmean(x2grids%imsource,x2grids%jmsource),
     &     OSTend(x2grids%imsource,x2grids%jmsource)
      integer iu_OSST,i,j,k

c      call openunit ("OSST",iu_OSST,.true.,.true.)     

c      name="OST4X5.B.1993-2002avg.Hadl1.1"
      name="OST4X5.B.1876-85avg.Hadl1.1"
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

      use regrid_com
      implicit none
      type (x_2gridsroot), intent(in) :: x2grids
      character*80 TITLE,name,outunformat, TITLE2(nrecmax)
      real*8, allocatable :: tsource(:,:,:)
      real*8, allocatable :: ttargglob(:,:,:)
      real*4, allocatable :: tin(:,:,:),tout(:,:,:)
      real*8, allocatable :: tsource1(:,:,:,:),tsource2(:,:,:,:)
      real*8, allocatable :: ttargglob1(:,:,:),ttargglob2(:,:,:)
      real*4, allocatable :: tout1(:,:,:),tout2(:,:,:),tbig(:,:,:,:)
      integer iu_SICE,i,j,k,ims,jms,nts,imt,jmt,ntt,iuout,
     &    maxrec,ir

      name="SICE4X5.B.1876-85avg.Hadl1.1"
      write(*,*) name
      open (iu_SICE, FILE=name,FORM='unformatted', STATUS='old')


      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      allocate (tsource(ims,jms,nts),
     &     ttargglob(imt,jmt,ntt),
     &     tin(ims,jms,nts),
     &     tout(imt,jmt,ntt))
      tsource(:,:,:)=0.0
 
      read(unit=iu_SICE) TITLE, tin

      tsource=tin

      call root_regrid(x2grids,tsource(:,:,:),ttargglob)
      tout(:,:,:)=ttargglob(:,:,:)
      
      outunformat=trim(name)//".CS"
      
      write(*,*) outunformat
      iuout=20
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")
      write(unit=iuout) TITLE,tout(:,:,:)
      write(*,*) TITLE

      allocate (tsource1(ims,jms,nts,nrecmax),
     &     tsource2(ims,jms,nts,nrecmax),
     &     ttargglob1(imt,jmt,ntt),
     &     ttargglob2(imt,jmt,ntt),
     &     tout1(imt,jmt,ntt),
     &     tout2(imt,jmt,ntt),tbig(imt,jmt,2,ntt) )

      tsource1(:,:,:,:)=0.0
      tsource2(:,:,:,:)=0.0

      call read_recs_2R(tsource1,tsource2,iu_SICE,TITLE2,
     &     maxrec,ims,jms,nts)

      write(*,*) "maxrec",maxrec

      do ir=1,maxrec
         call root_regrid(x2grids,tsource1(:,:,:,ir),ttargglob1)
         tout1(:,:,:)=ttargglob1(:,:,:)
         call root_regrid(x2grids,tsource2(:,:,:,ir),ttargglob2)
         tout2(:,:,:)=ttargglob2(:,:,:)
         tbig(:,:,1,:)=tout1(:,:,:)
         tbig(:,:,2,:)=tout2(:,:,:)
         write(unit=iuout) TITLE2(ir),tbig
c         write(unit=iuout) TITLE2(ir),tout1(:,:,:),tout2(:,:,:)
         write(*,*) "TITLE",TITLE2(ir)
      enddo

      close(iuout)

      deallocate(tsource1,tsource2,ttargglob1,ttargglob2,tout1,tout2,
     &     tin,tout,tbig)

      close(iu_SICE)
      
      end subroutine regridSICE



      subroutine regridCDN(x2grids)
c     for 1x1 resolution: Jeff uses CDN=AL30RL360X180N.rep

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
c     &     convert='big_endian')

      call read_regrid_write_veg(x2grids,name,iu_VEG)
           
      close(iu_VEG)
      
      end subroutine regridVEG
c*


      subroutine regridRVR(x2grids)
c	empty for the moment

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
     &     ftextkout(:,:,:,:,:),slout(:,:,:),bigarrout(:,:,:,:)
      real*8, allocatable :: tsource(:,:,:)
      real*8, allocatable :: ttargglob(:,:,:)
      character*80 TITLE,name,outunformat
      integer iu_SOIL,iuout,ims,jms,nts,imt,jmt,ntt,i,j,k,l,inds

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "ims,jms,nts,imt,jmt,ntt",ims,jms,nts,imt,jmt,ntt

      name="S4X50093.ext"

      open(iu_SOIL,FILE=name,FORM='unformatted', STATUS='old')
c ,
c     &     convert='big_endian')

      allocate (dz(ims,jms,6),ftext(ims,jms,6,5),
     &     ftextk(ims,jms,6,5),sl(ims,jms),
     &     dzout(imt,jmt,6,ntt),ftextout(imt,jmt,6,5,ntt),
     &     ftextkout(imt,jmt,6,5,ntt),
     &     slout(imt,jmt,ntt),bigarrout(imt,jmt,67,ntt) )

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
c ,
c     &     convert='big_endian')

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
      
      do k=1,6
         bigarrout(:,:,k,:)=dzout(:,:,k,:)
      enddo

      inds=6
      do k=1,6
         do l=1,5
            bigarrout(:,:,inds+k+6*(l-1),:)=ftextout(:,:,k,l,:)
         enddo
      enddo

      inds=inds+30
      do k=1,6
         do l=1,5
            bigarrout(:,:,inds+k+6*(l-1),:)=ftextkout(:,:,k,l,:)
         enddo
      enddo

      inds=inds+30
      bigarrout(:,:,inds+1,:)=slout(:,:,:)

      write(unit=iuout) bigarrout
      
      close(iuout)
      
      end subroutine regridSOIL
c*
      
      


      subroutine regridGIC(x2grids,dd2d)
c
c     for 1x1 resolution: Jeff uses GIC=GIC.360X180.DEC01.1.rep
c
      use DOMAIN_DECOMP_1D,only : am_i_root
      use regrid_com
      use dd2d_utils
      use pario, only : defvar,write_data

      implicit none
      include 'netcdf.inc'

      type (x_2gridsroot), intent(in) :: x2grids
      type (dd2d_grid), intent(in) :: dd2d

c*    read
      real*8, allocatable :: Tocn(:,:,:),MixLD(:,:)        ! OCN01
      real*8, allocatable :: F(:,:),H(:,:,:),snw(:,:),msi(:,:),
     &     ssi(:,:,:),pond_melt(:,:)
      logical, allocatable :: flag_dsws(:,:)               ! SICE02
      real*8, allocatable :: snowe(:,:),Te(:,:),WTRe(:,:),ICEe(:,:),
     &     SNOage(:,:,:),evmax(:,:),fsat(:,:),gq(:,:)      ! EARTH01
      real*8, allocatable :: W(:,:,:,:),HT(:,:,:,:),SNWbv(:,:,:) ! SOILS03
      real*8, allocatable :: SNOW(:,:),T(:,:,:),
     &     MDWN(:,:),EDWN(:,:)
      real*8 ::  ACCPDA, ACCPDG,  EACCPDA, EACCPDG           !GLAI

c*    write
      real*8, allocatable :: Tocn_out(:,:,:,:),MixLD_out(:,:,:)   ! OCN01
      real*8, allocatable :: F_out(:,:,:),H_out(:,:,:,:),
     &     snw_out(:,:,:),msi_out(:,:,:),ssi_out(:,:,:,:),
     &     pond_melt_out(:,:,:), flag_dsws_out(:,:,:)                ! SICE02
      real*8, allocatable :: snowe_out(:,:,:),Te_out(:,:,:),
     &     WTRe_out(:,:,:), ICEe_out(:,:,:),SNOage_out(:,:,:,:),
     &     evmax_out(:,:,:),fsat_out(:,:,:),gq_out(:,:,:)         ! EARTH01
      real*8, allocatable :: W_out(:,:,:,:,:),
     &     HT_out(:,:,:,:,:),SNWbv_out(:,:,:,:)                   ! SOILS02
      real*8, allocatable :: SNOW_out(:,:,:),T_out(:,:,:,:),
     &     MDWN_out(:,:,:),EDWN_out(:,:,:)        !GLAI
      real*8, allocatable :: tsource(:,:,:)
      real*8, allocatable :: ttargglob(:,:,:)
      character*80 TITLEOCN01,TITLESICE02,TITLEEARTH01,TITLESOILS03,
     &     TITLEGLAIC01,name,outnc
      integer :: iu_GIC,iuout,ims,jms,nts,imt,jmt,ntt
      integer :: i,j,k,l,m,fid,status,ntiles,im,jm,d2,d3
     

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "ims,jms,nts,imt,jmt,ntt",ims,jms,nts,imt,jmt,ntt

c      name="GIC.E046D3M20A.1DEC1955.ext"
      name="GIC.360X180.DEC01.1.rep"
c      name="GIC.144X90.DEC01.1.ext"

      open(iu_GIC,FILE=name,FORM='unformatted', STATUS='old')

      allocate (Tocn(3,ims,jms),MixLD(ims,jms),
     &     F(ims,jms),H(4,ims,jms),snw(ims,jms),msi(ims,jms),
     &     ssi(4,ims,jms),pond_melt(ims,jms),flag_dsws(ims,jms),
     &     snowe(ims,jms),Te(ims,jms),WTRe(ims,jms),ICEe(ims,jms), 
     &     SNOage(3,ims,jms),evmax(ims,jms),fsat(ims,jms),gq(ims,jms),
     &     W(7,3,ims,jms),HT(7,3,ims,jms),
     &     SNWbv(2,ims,jms),
     &     SNOW(ims,jms),T(2,ims,jms),
     &     MDWN(ims,jms),EDWN(ims,jms))

      allocate (Tocn_out(3,imt,jmt,ntt),MixLD_out(imt,jmt,ntt),
     &     F_out(imt,jmt,ntt),H_out(4,imt,jmt,ntt),
     &     snw_out(imt,jmt,ntt),msi_out(imt,jmt,ntt),
     &     ssi_out(4,imt,jmt,ntt),pond_melt_out(imt,jmt,ntt),
     &     flag_dsws_out(imt,jmt,ntt),
     &     snowe_out(imt,jmt,ntt),Te_out(imt,jmt,ntt),
     &     WTRe_out(imt,jmt,ntt),ICEe_out(imt,jmt,ntt), 
     &     SNOage_out(3,imt,jmt,ntt),evmax_out(imt,jmt,ntt),
     &     fsat_out(imt,jmt,ntt),gq_out(imt,jmt,ntt),
     &     W_out(7,3,imt,jmt,ntt),
     &     HT_out(7,3,imt,jmt,ntt),
     &     SNWbv_out(3,imt,jmt,ntt),
     &     SNOW_out(imt,jmt,ntt),T_out(2,imt,jmt,ntt),
     &     MDWN_out(imt,jmt,ntt),EDWN_out(imt,jmt,ntt))

      allocate (tsource(ims,jms,nts),
     &     ttargglob(imt,jmt,ntt) )

      read(iu_GIC) TITLEOCN01, Tocn,MixLD
      write(*,*) TITLEOCN01
      read(iu_GIC) TITLESICE02, F,H,snw,msi,ssi,pond_melt,flag_dsws
      write(*,*) TITLESICE02
      read(iu_GIC) TITLEEARTH01, snowe,Te,WTRe, ICEe, SNOage,evmax,
     &     fsat,gq
      write(*,*) TITLEEARTH01
      read(iu_GIC) TITLESOILS03, W,HT,SNWbv
      write(*,*) TITLESOILS03
      read(iu_GIC) TITLEGLAIC01, SNOW,T,MDWN,EDWN,
     &     ACCPDA,ACCPDG,EACCPDA,EACCPDG
      write(*,*) TITLEGLAIC01

      close(iu_GIC)
      

            
      do k=1,3
         tsource(:,:,1)=Tocn(k,:,:)
         call root_regrid(x2grids,tsource,ttargglob)
         Tocn_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

     
      tsource(:,:,1)=MixLD(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      MixLD_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=F(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      F_out(:,:,:)=ttargglob(:,:,:)

      do k=1,4
         tsource(:,:,1)=H(k,:,:)
         call root_regrid(x2grids,tsource,ttargglob)
         H_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

      tsource(:,:,1)=snw(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      snw_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=msi(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      msi_out(:,:,:)=ttargglob(:,:,:)

      do k=1,4
         tsource(:,:,1)=ssi(k,:,:)
         call root_regrid(x2grids,tsource,ttargglob)
         ssi_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

      tsource(:,:,1)=pond_melt(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      pond_melt_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1) = 0.d0

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
      do k=1,ntt
         do j=1,jmt
            do i=1,imt
               if (ttargglob(i,j,k) .ge. 0.5d0) then
                  flag_dsws_out(i,j,k)=1.
               else
                  flag_dsws_out(i,j,k)=0.
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
         tsource(:,:,1)=SNOage(k,:,:)
         call root_regrid(x2grids,tsource,ttargglob)
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
      
      
      do k=1,7
         do m=1,3
            tsource(:,:,1)=W(k,m,:,:)
            call root_regrid(x2grids,tsource,ttargglob)
            W_out(k,m,:,:,:)=ttargglob(:,:,:)
         enddo
      enddo

      do k=1,7
         do m=1,3
            tsource(:,:,1)=HT(k,m,:,:)
            call root_regrid(x2grids,tsource,ttargglob)
            HT_out(k,m,:,:,:)=ttargglob(:,:,:)
         enddo
      enddo

      do k=1,2
         tsource(:,:,1)=SNWbv(k,:,:)
         call root_regrid(x2grids,tsource,ttargglob)
         SNWbv_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

      SNWbv_out(3,:,:,:) = 0.

      tsource(:,:,1)=SNOW(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      SNOW_out(:,:,:)=ttargglob(:,:,:)

      do k=1,2
         tsource(i,j,1)=T(k,i,j)
         call root_regrid(x2grids,tsource,ttargglob)
         T_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

      tsource(:,:,1)=MDWN(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      MDWN_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=EDWN(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      EDWN_out(:,:,:)=ttargglob(:,:,:)
      

c***  Write Netcdf file
#ifdef TRACERS_WATER
      write(*,*) "STOP TRACERS WATER NOT IMPLEMENTED IN regridinput"
      stop
#endif      
      

      deallocate (Tocn,MixLD,F,H,snw,msi,
     &     ssi,pond_melt,flag_dsws,
     &     snowe,Te,WTRe,ICEe, 
     &     SNOage,evmax,fsat,gq,
     &     W,HT,SNWbv,SNOW,T,MDWN,EDWN )

      deallocate (tsource,ttargglob)     
      
      outnc=trim(name)//"-CS.nc"
      write(*,*) outnc
      
      if (am_i_root()) then
         status = nf_create(outnc,nf_clobber,fid)
         if (status .ne. NF_NOERR) write(*,*) "UNABLE TO CREATE FILE"
      endif

    
c***  Define OCN variables
      call defvar(dd2d,fid,Tocn_out,'tocean(d3,im,jm,tile)')
      call defvar(dd2d,fid,MixLD_out,'z1o(im,jm,tile)')
c***  Define SICE variables
      call defvar(dd2d,fid,F_out,'rsi(im,jm,tile)')
      call defvar(dd2d,fid,H_out,'hsi(lmi,im,jm,tile)')
      call defvar(dd2d,fid,snw_out,'snowi(im,jm,tile)')
      call defvar(dd2d,fid,msi_out,'msi(im,jm,tile)')
      call defvar(dd2d,fid,ssi_out,'ssi(lmi,im,jm,tile)')
      call defvar(dd2d,fid,pond_melt_out,
     &     'pond_melt(im,jm,tile)')
      call defvar(dd2d,fid,flag_dsws_out,
     &     'flag_dsws(im,jm,tile)')
c***  Define EARTH variables
      call defvar(dd2d,fid,snowe_out,'snowe(im,jm,tile)')
      call defvar(dd2d,fid,Te_out,'tearth(im,jm,tile)')
      call defvar(dd2d,fid,WTRe_out,'wearth(im,jm,tile)')
      call defvar(dd2d,fid,ICEe_out,'aiearth(im,jm,tile)')
      call defvar(dd2d,fid,SNOage_out,'snoage(d3,im,jm,tile)')
      call defvar(dd2d,fid,evmax_out,
     &     'evap_max_ij(im,jm,tile)')
      call defvar(dd2d,fid,fsat_out,'fr_sat_ij(im,jm,tile)')
      call defvar(dd2d,fid,gq_out,'qg_ij(im,jm,tile)')
c***  Define SOIL variables     ! this is the old SOIL02 version, implement the new version
      call defvar(dd2d,fid,W_out,
     &     'w_ij(zero_to_ngm,ls_nfrac,im,jm,tile)')
      call defvar(dd2d,fid,HT_out,
     &     'ht_ij(zero_to_ngm,ls_nfrac,im,jm,tile)')
      call defvar(dd2d,fid,SNWbv_out,
     &     'snowbv(ls_nfrac,im,jm,tile)')
c***  Define GLAIC variables
      call defvar(dd2d,fid,SNOW_out,'snowli(im,jm,tile)')
      call defvar(dd2d,fid,T_out,'tlandi(d2,im,jm,tile)')
      call defvar(dd2d,fid,MDWN_out,'mdwnimp(im,jm,tile)')
      call defvar(dd2d,fid,EDWN_out,'edwnimp(im,jm,tile)')
      call defvar(dd2d,fid,ACCPDA,'accpda')
      call defvar(dd2d,fid,ACCPDG,'accpdg')
      call defvar(dd2d,fid,EACCPDA,'eaccpda')
      call defvar(dd2d,fid,EACCPDG,'eaccpdg')
      if (am_i_root()) then
         status = nf_enddef(fid)
         if (status .ne. NF_NOERR) write(*,*) "Problem with enddef"
      endif

c***  Write OCN variables
      call write_data(dd2d,fid,'tocean',Tocn_out)
      call write_data(dd2d,fid,'z1o',MixLD_out)
c***  Write SICE variables
      call write_data(dd2d,fid,'rsi',F_out)
      call write_data(dd2d,fid,'hsi',H_out)
      call write_data(dd2d,fid,'snowi',snw_out)
      call write_data(dd2d,fid,'msi',msi_out)
      call write_data(dd2d,fid,'ssi',ssi_out)
      call write_data(dd2d,fid,'pond_melt',pond_melt_out)
      call write_data(dd2d,fid,'flag_dsws',flag_dsws_out)
c***  Write EARTH variables
      call write_data(dd2d,fid,'snowe',snowe_out)
      call write_data(dd2d,fid,'tearth',Te_out)
      call write_data(dd2d,fid,'wearth',WTRe_out)
      call write_data(dd2d,fid,'aiearth',ICEe_out)
      call write_data(dd2d,fid,'snoage',SNOage_out)
      call write_data(dd2d,fid,'evap_max_ij',evmax_out)
      call write_data(dd2d,fid,'fr_sat_ij',fsat_out)
      call write_data(dd2d,fid,'qg_ij',gq_out)
c***  Write SOIL variables     ! this is the old SOIL02 version, implement the new version
      call write_data(dd2d,fid,'w_ij',W_out)
      call write_data(dd2d,fid,'ht_ij',HT_out)
      call write_data(dd2d,fid,'snowbv',SNWbv_out)
c***  Write GLAIC variables
      call write_data(dd2d,fid,'snowli',SNOW_out)
      call write_data(dd2d,fid,'tlandi',T_out)
      call write_data(dd2d,fid,'mdwnimp',MDWN_out)
      call write_data(dd2d,fid,'edwnimp',EDWN_out)
      call write_data(dd2d,fid,'accpda',ACCPDA)
      call write_data(dd2d,fid,'accpdg',ACCPDG)
      call write_data(dd2d,fid,'eaccpda',EACCPDA)
      call write_data(dd2d,fid,'eaccpdg',EACCPDG)


      deallocate (Tocn_out,MixLD_out,F_out,H_out,
     &     snw_out,msi_out,ssi_out,pond_melt_out,
     &     flag_dsws_out,snowe_out,Te_out,
     &     WTRe_out,ICEe_out,SNOage_out,evmax_out,
     &     fsat_out,gq_out,W_out,
     &     HT_out,SNWbv_out,
     &     SNOW_out,T_out,MDWN_out,EDWN_out )

      
      if(am_i_root()) status = nf_close(fid)
      
      write(*,*) "end regrid GIC"

      end subroutine regridGIC
c*



      subroutine regridAIC(x2grids,dd2d)
c
c     for 1x1 resolution : Jeff uses AIC=AIC.RES_X40.D771201N.rep
c
      use regrid_com
      use pario, only : defvar,write_data
      use dd2d_utils
      implicit none
      include 'netcdf.inc'
      type(x_2gridsroot), intent(in) :: x2grids
      type (dd2d_grid), intent(in) :: dd2d
      integer :: iuout
      real*8, allocatable :: tsource(:,:,:)
      real*4, allocatable :: ts4(:,:,:)
      real*8, allocatable :: ttargglob(:,:,:),tcopy(:,:,:)
      character*80 :: TITLE,outunformat,outnc,name
      integer :: maxrec,irec,ir,ims,jms,nts,imt,jmt,ntt,
     &     status,fid,vid
      integer iu_AIC
     

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "ims,jms,nts,imt,jmt,ntt r8",
     &     ims,jms,nts,imt,jmt,ntt
      allocate (tsource(ims,jms,nts),ts4(ims,jms,nts),
     &     ttargglob(imt,jmt,ntt),tcopy(imt,jmt,ntt) )

      tsource(:,:,:)=0.0

      if (am_i_root()) then
      name="AIC.RES_M20A.D771201"
      open(iu_AIC,FILE=name,FORM='unformatted', STATUS='old')

      outunformat=trim(name)//".CS"
      write(*,*) outunformat

      iuout=20
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")

      irec=1

      do
         read(unit=iu_AIC,END=30) TITLE, ts4
         tsource=ts4
         write(*,*) "TITLE, irec",TITLE,irec
         call root_regrid(x2grids,tsource,ttargglob)

         if (irec .eq. 1) tcopy=ttargglob

         write(unit=iuout) TITLE,real(ttargglob,KIND=4)
c          write(unit=iuout) TITLE,ttargglob
         irec=irec+1
      enddo
   
 30   continue

      maxrec=irec-1

      write(*,*) "maxrec",maxrec

      close(iuout) 
      close(iu_AIC)

      endif

      write(*,*) "here w r4"

c      outnc=trim(name)//"-CS.nc"
c      write(*,*) outnc
c
c      if (am_i_root()) then
c         write(*,*) "TCOPY=",tcopy
c         status = nf_create(outnc,nf_clobber,fid)
c         if (status .ne. NF_NOERR) write(*,*) "UNABLE TO CREATE FILE"
c      endif
c
c      call defvar(dd2d,fid,tcopy,'press(im,jm,tile)')
c
c      if (am_i_root()) then
c         status = nf_enddef(fid)
c         if (status .ne. NF_NOERR) write(*,*) "Problem with enddef"
c      endif
c
c      call write_data(dd2d,fid,'press',tcopy)
c
c      if(am_i_root()) status = nf_close(fid)

      deallocate(tsource,ts4,ttargglob,tcopy)
 
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

c      write(*,*) "DATA=",data
      maxrec=irec-1

      close(iuin)

      end subroutine read_recs_1R
c*

      subroutine read_recs_1R_r4_r8(tsource,iuin,TITLE,maxrec,
     *     im1,jm1,ntl)
      use regrid_com
      implicit none
      integer i,j,k,irec
      integer, intent(in) :: im1,jm1,ntl
      real*8, intent(inout) :: tsource(im1,jm1,ntl,nrecmax)
      real*4 :: ts4(im1,jm1,ntl,nrecmax)
      integer, intent(in) :: iuin
      character*80, intent(inout) :: TITLE(nrecmax)
      integer, intent(out) :: maxrec

      write(*,*) "iuin",iuin
      irec=1

      do
         read(unit=iuin,END=30) TITLE(irec), ts4(:,:,:,irec)
         tsource(:,:,:,irec)=ts4(:,:,:,irec)
         write(*,*) "TITLE, irec",TITLE(irec),irec
         irec=irec+1
      enddo

 30   continue

      maxrec=irec-1

      close(iuin)

      end subroutine read_recs_1R_r4_r8
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

c      write(*,*) "TOUT=",tout

      deallocate(tsource,ttargglob,tout)

      end subroutine read_regrid_write_4D_1R
c*


      subroutine read_regrid_write_4D_1R_r8(dd2d,x2grids,name,iuin)
      use regrid_com
      use pario, only : defvar,write_data
      use dd2d_utils
      implicit none
      include 'netcdf.inc'
      type(x_2gridsroot), intent(in) :: x2grids
      type (dd2d_grid), intent(in) :: dd2d
      character*80, intent(in) :: name
      integer, intent(in) :: iuin
      integer :: iuout
      real*8, allocatable :: tsource(:,:,:,:)
      real*8, allocatable :: ttargglob(:,:,:),tcopy(:,:,:)
      character*80 :: TITLE(nrecmax),outunformat,outnc
      integer :: maxrec,irec,ir,ims,jms,nts,imt,jmt,ntt,
     &     status,fid,vid

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "iuin ims,jms,nts,imt,jmt,ntt r8",iuin,
     &     ims,jms,nts,imt,jmt,ntt
      allocate (tsource(ims,jms,nts,nrecmax),
     &     ttargglob(imt,jmt,ntt),tcopy(imt,jmt,ntt) )
      tsource(:,:,:,:)=0.0
      
      call read_recs_1R_r4_r8(tsource,iuin,TITLE,
     &        maxrec,ims,jms,nts)
      
      write(*,*) "maxrec",maxrec
      
      outunformat=trim(name)//".CS"
      
      write(*,*) outunformat
      iuout=20
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")

      do ir=1,maxrec
         call root_regrid(x2grids,tsource(:,:,:,ir),ttargglob)
         if (ir .eq. 1) tcopy=ttargglob
         write(unit=iuout) TITLE(ir),ttargglob(:,:,:)
         write(*,*) "TITLE",TITLE(ir)
      enddo

      close(iuout) 

      write(*,*) "here w r8"

      outnc=trim(name)//"-CS.nc"
      write(*,*) outnc

      if (am_i_root()) then
         write(*,*) "TCOPY=",tcopy
         status = nf_create(outnc,nf_clobber,fid)
         if (status .ne. NF_NOERR) write(*,*) "UNABLE TO CREATE FILE"
      endif

      call defvar(dd2d,fid,tcopy,'press(im,jm,tile)')

      if (am_i_root()) then
         status = nf_enddef(fid)
         if (status .ne. NF_NOERR) write(*,*) "Problem with enddef"
      endif

      call write_data(dd2d,fid,'press',tcopy)

      deallocate(tsource,ttargglob,tcopy)

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
      real*8, allocatable :: ttargglob(:,:,:),arrsum(:,:,:)
      real*4, allocatable :: tout(:,:,:),data(:,:,:,:)
      character*80, allocatable :: TITLE(:)
      character*80 :: outunformat
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
     &     ttargglob(imt,jmt,ntt),arrsum(ims,jmt,ntt),
     &     tout(imt,jmt,ntt),
     &     data(ims,jms,nts,rmax),
     &     TITLE(rmax))
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
c,
c     &     convert='big_endian')

      arrsum(:,:,:)=0.

      do ir=1,maxrec
         call root_regrid(x2grids,tsource(:,:,:,ir),ttargglob)
         arrsum(:,:,:)=arrsum(:,:,:)+ttargglob(:,:,:)
         tout(:,:,:)=ttargglob(:,:,:)
         write(unit=iuout) TITLE(ir),tout(:,:,:)
         write(*,*) "TITLE",TITLE(ir)
c         write(*,*) "TOUT>>",tout,"<<<"
      enddo

      write(*,*) "SUM ARRAY=",arrsum

      close(iuout) 

      deallocate(tsource,ttargglob,arrsum,tout,data,TITLE)

      end subroutine read_regrid_write_4D_1R_rmax
c*

      subroutine read_regrid_write_veg(x2grids,name,iuin)
      use regrid_com
      implicit none
      type(x_2gridsroot), intent(in) :: x2grids
      character*80, intent(in) :: name
      integer, intent(in) :: iuin
      integer :: iuout
      real*8, allocatable :: tsource(:,:,:)
      real*8, allocatable :: ttargglob(:,:,:)
      real*4, allocatable :: tout(:,:,:),data(:,:,:)
c      real*4 :: vadata(11,4,3)
      character*80 :: TITLE
      character*80 :: outunformat
      integer :: ir,ims,jms,nts,imt,jmt,ntt

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "iuin ims,jms,nts,imt,jmt,ntt 4D",iuin,
     &     ims,jms,nts,imt,jmt,ntt

      allocate (tsource(ims,jms,nts),data(ims,jms,nts),
     &     ttargglob(imt,jmt,ntt),tout(imt,jmt,ntt))
c     arrsum(imt,jmt,ntt)
     
      
      outunformat=trim(name)//".CS"
      write(*,*) outunformat

      iuout=20

      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")

      do ir=1,10
         read(unit=iuin) TITLE,data
         write(*,*) "TITLE, ir",TITLE,ir
         tsource= data
         call root_regrid(x2grids,tsource,ttargglob)
c         arrsum(:,:,:)=arrsum(:,:,:)+ttargglob(:,:,:)
         tout=ttargglob
         write(unit=iuout) TITLE,tout
      enddo

      close(iuout) 

      deallocate(tsource,ttargglob,tout,data)

      end subroutine read_regrid_write_veg
c*


      subroutine read_regrid_4D_1R(x2grids,iuin,TITLE,
     &     ttargglob,maxrec)
      use regrid_com
      implicit none
      type(x_2gridsroot), intent(in) :: x2grids
      integer, intent(in) :: iuin
      integer, intent(inout) :: maxrec
      real*8, allocatable :: tsource(:,:,:,:)
      real*8 :: ttargglob(x2grids%imtarget,x2grids%jmtarget,
     &     x2grids%ntilestarget,nrecmax)
      character*80, intent(inout) :: TITLE(nrecmax)
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
      real*4, allocatable :: tout1(:,:,:),tout2(:,:,:),tbig(:,:,:,:)
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
     &     tout2(imt,jmt,ntt),
     &     tbig(imt,jmt,2,ntt) )

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
         tbig(:,:,1,:)=tout1
         tbig(:,:,2,:)=tout2
         write(unit=iuout) TITLE(ir),tbig
         write(*,*) "TITLE",TITLE(ir)
      enddo
      
      close(iuout) 

      deallocate(tsource1,tsource2,ttargglob1,ttargglob2,tout1,tout2,
     &     tbig)

      end subroutine read_regrid_write_4D_2R


