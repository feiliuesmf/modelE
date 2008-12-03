      subroutine regrid_input
      use regrid_com
      implicit none
      type (x_2gridsroot) :: xll2cs

      call init_regrid_root(xll2cs,72,46,1,48,48,6)
c      call init_regrid(xll2cs,360,180,1,90,90,6)

ccc   regrid boundary condition files 
      write(*,*) "IN REGRID INPUT"
      if (AM_I_ROOT()) then
      call regridTOPO(xll2cs)
c      call regridOSST(xll2cs)
c      call regridSICE(xll2cs)
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
      real*8, allocatable :: ttargglob(:,:,:,:),ones(:,:,:)
      integer :: iu_TOPO,i,j,k,irec,iunit,imt,jmt,ntt,maxrec

      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      allocate( ttargglob(imt,jmt,ntt,nrecmax),ones(imt,jmt,ntt) )

      name="Z72X46N.cor4_nocasp"

      open( iu_TOPO, FILE=name,FORM='unformatted', STATUS='old')

      write(*,*) "iu_TOPO",iu_TOPO

      call read_regrid_4D_1R(x2grids,iu_TOPO,ttargglob,maxrec)


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

      deallocate(ttargglob,ones)
      
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
c     SICE : 1x1 by Gary on Athena clima1/OBS/AMIP/1x1
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

c      do k=1,12
c         read(iu_SICE) TITLE,SICEmean,SICEend        
c         write(*,*) TITLE
c	   write(*,*) SICEmean(:,:)
c      enddo
      
      close(iu_SICE)
      
      end subroutine regridSICE



      subroutine regridCDN(x2grids)
c     Jeff uses CDN=AL30RL360X180N.rep
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
c     Jeff uses VEG=V360X180_no_crops.rep
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
c     Jeff uses SOIL=S360X180_0098M.rep
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
      
      

      subroutine regridGIC(x2grids)
c
c     Jeff uses GIC=GIC.360X180.DEC01.1.rep
c
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
c
c     Jeff uses AIC=AIC.RES_X40.D771201N.rep
c
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
         write(*,*) "TITLE, irec",TITLE(irec),irec
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
c         write(*,*) "TOUT>>",tout(:,:,:),"<<<"
      enddo

      close(iuout) 

      deallocate(tsource,ttargglob,tout)

      end subroutine read_regrid_write_4D_1R
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

      deallocate(tsource,ttargglob,tout)

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


      subroutine read_regrid_write_4D_2R(x2grids,name,iuin,rmax)
      use regrid_com
      type(x_2gridsroot), intent(in) :: x2grids
      integer, intent(in) :: iuin
      integer, intent(in), optional :: rmax
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


