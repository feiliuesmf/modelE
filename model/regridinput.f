      subroutine regrid_input
      use regrid_com
      implicit none
      type (x_2gridsroot) :: xll2cs

      call init_regrid_root(xll2cs,72,46,1,48,48,6)
c      call init_regrid(xll2cs,360,180,1,90,90,6)

ccc   regrid boundary condition files 
      write(*,*) "IN REGRID INPUT"
      if (AM_I_ROOT()) then
c      call regridTOPO(xll2cs)
c      call regridOSST(xll2cs)
      call regridSICE(xll2cs)
c      call regridCDN(xll2cs)
c      call regridVEG(xll2cs)
c      call regridRVR(xll2cs)            ! empty for the moment
c      call regridCROPS(xll2cs)
c      call regridTOPINDEX(xll2cs)
c      call regridSOIL(xll2cs)

ccc   Then regrid Initial Condition 
ccc   If restart != 2 then IC should already be on CS grid
c      if (restart .neq. 2) then
c         call regridGIC(xll2cs)
c         call regridAIC(xll2cs)
c      endif
      endif
      end subroutine regrid_input
c*


      subroutine regridTOPO(x2grids)
c
c     Jeff posted Z1X1N on Discover
c
      USE FILEMANAGER, only : openunit,closeunit
      use regrid_com
      implicit none
      character*80 :: TITLE,name
      type (x_2gridsroot), intent(inout) :: x2grids
      real*4 FOCEAN(x2grids%imsource,x2grids%jmsource),
     &     FLAKE(x2grids%imsource,x2grids%jmsource),
     &     FGRND(x2grids%imsource,x2grids%jmsource),
     &     FGICE(x2grids%imsource,x2grids%jmsource),
     &     ZATMO(x2grids%imsource,x2grids%jmsource)
      integer :: iu_TOPO,i,j,irec,iunit

      
c      call openunit("TOPO",iu_TOPO,.true.,.true.)

      name="Z72X46N.cor4_nocasp"

      open( iu_TOPO, FILE=name,FORM='unformatted', STATUS='old')

      write(*,*) "iu_TOPO",iu_TOPO

      call read_regrid_write_4D_1R(x2grids,name,iu_TOPO)

c     TODO: CONSISTENCY CHECK FOCEAN+FLAKE+FGRND+FGICE=1 .AND. FOCEAN=0 OR 1
c     ALSO IF FOCEAN(i,j) > 0 set FGRND=FGRND+FLAKE, FLAKE=0

      close(iu_TOPO)
      
      end subroutine regridTOPO 
c*


      subroutine regridOSST(x2grids)
c     OSST: 1x1 by Gary on Athena clima1/OBS/AMIP/1x1
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

c      do k=1,12  !12 months
c      read(iu_OSST) TITLE,OSTmean,OSTend 
c      write(*,*) TITLE
c      write(*,*) OSTmean(:,:)
c      enddo

      close(iu_OSST)

      end subroutine regridOSST
c


      subroutine regridSICE(x2grids)
c 	SICE : 1x1 by Gary on Athena clima1/OBS/AMIP/1x1
c 	Also ICE_1x1_HadISST from Hadley (Jeff Jonas) - check if this  
c 	has been interpolated from lower resolution 
c       /u/cmrun/SICE4X5.B.1993-2002avg.Hadl1.1

      USE FILEMANAGER, only : openunit,closeunit
      use regrid_com
      implicit none
      type (x_2gridsroot), intent(in) :: x2grids
      character*80 TITLE,name
      real*4 SICEmean(x2grids%imsource,x2grids%jmsource),
     &     SICEend(x2grids%imsource,x2grids%jmsource)
      integer iu_SICE,i,j,k

c      call openunit ("SICE",iu_SICE,.true.,.true.) 

      name="SICE4X5.B.1993-2002avg.Hadl1.1"
      open (iu_SICE, FILE=name,FORM='unformatted', STATUS='old')
      read(unit=iu_SICE) TITLE 
      write(*,*) TITLE

      call read_regrid_write_4D_2R(x2grids,name,iu_SICE)

c      do k=1,12
c         read(iu_SICE) TITLE,SICEmean,SICEend        
c         write(*,*) TITLE
c	   write(*,*) SICEmean(:,:)
c      enddo
      
      close(iu_SICE)
      
      end subroutine regridSICE


      subroutine regridCDN(x2grids)
c	Jeff uses CDN=AL30RL360X180N.rep
      USE FILEMANAGER, only : openunit,closeunit
      use regrid_com
      implicit none
      type (x_2grids), intent(in) :: x2grids
      character*80 TITLE,name
      real*4 CDN(x2grids%imsource,x2grids%jmsource)
      integer iu_CDN,i,j,k	

c      call openunit ("CDN",iu_CDN,.true.,.true.)

      name="CD4X500S.ext"
      open (iu_CDN, FILE=name,FORM='unformatted', STATUS='old')

      call read_regrid_write_4D_1R(x2grids,name,iu_CDN)
      
      call closeunit(iu_CDN)
      
      end subroutine regridCDN
c


      subroutine regridVEG(x2grids)
c     Jeff uses VEG=V360X180_no_crops.rep
c	It is identical to V144X90_no_crops.ext (144X90 data 
c	was just transfered to 360X180 grid without any change)
      USE FILEMANAGER, only : openunit,closeunit
      use regrid_com
      implicit none
      type (x_2grids), intent(in) :: x2grids
      character*80 TITLE
      real*4 VEG(x2grids%imsource,x2grids%jmsource)
      integer iu_VEG,i,j,k	


      call openunit ("VEG",iu_VEG,.true.,.true.)
      read(iu_veg) TITLE,VEG 
      
      write(*,*) TITLE
      write(*,*) VEG(:,:)
      
      call closeunit(iu_VEG)
      
      end subroutine regridVEG
c


      subroutine regridRVR(x2grids)
c	empty for the moment
      USE FILEMANAGER, only : openunit,closeunit
      use regrid_com
      implicit none
      type (x_2grids), intent(in) :: x2grids
      character*80 TITLE
      real*4 RVR(x2grids%imsource,x2grids%jmsource)
      integer iu_RVR,i,j,k	

      call openunit ("RVR",iu_RVR,.true.,.true.)
      read(iu_RVR) TITLE,RVR 
      
      write(*,*) TITLE
      write(*,*) RVR(:,:)

      call closeunit(iu_RVR)
   
      end subroutine regridRVR
c*


      subroutine regridCROPS(x2grids)
      USE FILEMANAGER, only : openunit,closeunit
      use regrid_com
      implicit none
      type (x_2grids), intent(in) :: x2grids
      character*80 TITLE
      real*4 CROPS(x2grids%imsource,x2grids%jmsource)
      integer iu_CROPS,i,j,k	


      call openunit("CROPS",iu_CROPS,.true.,.true.)
      read(iu_CROPS) TITLE,CROPS 
      
      write(*,*) TITLE
      write(*,*) CROPS(:,:)
      call closeunit(iu_CROPS)

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
      type (x_2grids), intent(in) :: x2grids
      character*80 TITLE
      real*4 TOP_INDEX(x2grids%imsource,x2grids%jmsource)
      integer i,j,k,iu_TOP_INDEX
     
      call openunit("TOP_INDEX",iu_TOP_INDEX,.true.,.true.)
      
      read(iu_TOP_INDEX) TITLE,TOP_INDEX 
      
      write(*,*) TITLE
      write(*,*) TOP_INDEX(:,:)
      
      call closeunit(iu_TOP_INDEX)
      
      end subroutine regridTOPINDEX
c*


      
      subroutine regridSOIL(x2grids)
c     Jeff uses SOIL=S360X180_0098M.rep

      USE FILEMANAGER, only : openunit,closeunit
      use regrid_com
      implicit none
      type (x_2grids), intent(in) :: x2grids
      character*80 TITLE
      real*4 SOIL(x2grids%imsource,x2grids%jmsource)
      integer iu_SOIL,i,j,k	

      call openunit ("SOIL",iu_SOIL,.true.,.true.)
      read(iu_SOIL) TITLE,SOIL 
      
      write(*,*) TITLE
      write(*,*) SOIL(:,:)
      call closeunit(iu_SOIL)
      
      end subroutine regridSOIL
c*
      
      

      subroutine regridGIC(x2grids)
c     Jeff uses GIC=GIC.360X180.DEC01.1.rep

      USE FILEMANAGER, only : openunit,closeunit
      use regrid_com
      implicit none
      type (x_2grids), intent(in) :: x2grids
      character*80 TITLE
      real*4 GIC(x2grids%imsource,x2grids%jmsource)
      integer i,j,k,iu_GIC	

      call openunit ("GIC",iu_GIC,.true.,.true.)
      read(iu_GIC) TITLE,GIC 
      
      write(*,*) TITLE
      write(*,*) GIC(:,:)
      call closeunit(iu_GIC)

      end subroutine regridGIC
c*



      subroutine regridAIC(x2grids)
c     Jeff uses AIC=AIC.RES_X40.D771201N.rep

      USE FILEMANAGER, only : openunit,closeunit
      use regrid_com
      implicit none
      type (x_2grids), intent(in) :: x2grids
      character*80 TITLE
      real*4 AIC(x2grids%imsource,x2grids%jmsource)
      integer i,j,k,iu_AIC

      call openunit("AIC",iu_AIC,.true.,.true.)

      read(iu_AIC) TITLE,AIC 
      
      write(*,*) TITLE
      write(*,*) AIC(:,:)

      call closeunit(iu_AIC)

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
         tsource(:,:,:,irec)= data(:,:,:,irec)
         irec=irec+1
      enddo

 30   continue

      maxrec=irec-1

      close(iuin)

      end subroutine read_recs_1R
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
      type(x_2gridsroot), intent(in) :: x2grids
      integer, intent(in) :: iuin
      character*80 name
      integer :: iuout
      real*8, allocatable :: tsource(:,:,:,:)
      real*8, allocatable :: ttargglob(:,:,:)
      real*4, allocatable :: tout(:,:,:)
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
      allocate (tsource(ims,jms,nts,nrecmax),
     &     ttargglob(imt,jmt,ntt),
     &     tout(imt,jmt,ntt))
      tsource(:,:,:,:)=0.0
      
      call read_recs_1R(tsource,iuin,TITLE,
     &     maxrec,ims,jms,nts)
      
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

      deallocate(tsource,ttargglob,tout)

      end subroutine read_regrid_write_4D_1R
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


