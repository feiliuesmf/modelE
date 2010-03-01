      program RDijk2ll_CS
!@sum On CS grid, convert i,j,k (with k=cube face index from 1..6)
!@+   to absolute lat-lon coordinates
!@auth Denis Gueyffier
      use geom, only : geom_cs,lon2d_dg,lat2d_dg
      use regrid_com
      USE DOMAIN_DECOMP_1D, ONLY : init_app
      use domain_decomp_atm, only : grid,init_grid,pack_data 
      USE fms_mod,         only : fms_init, fms_end
      implicit none
      integer, parameter :: imt=90,jmt=90
      integer, parameter :: nrvr = 41 ! # of river mouths
      character*80 :: name,nameout,title,
     &     title1,title2,title3,title4,title5,title6
      real*8,parameter :: undef=-1.d30  !missing value
      real*4, dimension(imt,jmt,6) :: down_lat,down_lon
      real*4, dimension(imt,jmt,6) :: down_lat_911,down_lon_911
      real*8, dimension(imt,jmt,6) :: lon_glob,lat_glob
      integer, dimension(imt,jmt,6) :: idown,jdown,kdown
      integer :: iu_RD,iu_TOPO,iu_MNAME,iu_MIJ,i,j,k
      LOGICAL, dimension(imt,jmt) :: NODIR
      real*4 :: FOCEAN(imt,jmt,6)
      character*8,dimension(nrvr) :: namervr
      character*2,dimension(nrvr) :: mouthI,mouthJ
      character*1,dimension(nrvr) :: mouthK
      integer,dimension(nrvr) :: imouthI,imouthJ,imouthK
      real*4,dimension(nrvr) :: lat_rvr,lon_rvr
      iu_TOPO=30

      call fms_init()

      call init_app()
      call init_grid(grid, imt, jmt, 20, CREATE_CAP=.true.)
      call geom_cs

      if (am_i_root()) then
c*    Read ocean fraction
         if (imt .eq. 32) then
            name="Z_CS32_4X5"
         else if (imt .eq. 90) then
            name="Z_CS90"
         endif
      open(iu_TOPO,FILE=name,FORM='unformatted', STATUS='old')
      read(iu_TOPO) title,FOCEAN
      close(iu_TOPO)

c*    Read names of river mouths
      iu_MNAME=20
      if (imt .eq. 32) then
         name="mouthnames_DtoO_CS32"
      elseif (imt .eq. 90) then
         name="mouthnames_DtoO_CS90"
      endif
      open(iu_MNAME,FILE=name,FORM='formatted', STATUS='old')
      READ (iu_MNAME,'(A8)') (namervr(I),I=1,nrvr) !Read mouths names
      write(*,*) namervr
      close(iu_MNAME)

c*    Read i,j,k coordinates of river mouths (k = index of cube face)
      if (imt .eq. 32) then
         name="mouthij_DtoO_CS32"   
      elseif (imt .eq. 90) then
         name="mouthij_DtoO_CS90"
      endif
      
      open(iu_MIJ,FILE=name,FORM='formatted', STATUS='old')
      READ (iu_MIJ,'(A2,1X,A2,1X,A1)') (      
     &       mouthI(I),mouthJ(I),mouthK(I), I=1,nrvr)
      close(iu_MIJ)

c*     conversion char to int
      do i=1,nrvr
          read(mouthI(i),'(I2)') imouthI(i)
          read(mouthJ(i),'(I2)') imouthJ(i)
          read(mouthK(i),'(I1)') imouthK(i)
          write(*,*) imouthI(i),imouthJ(i),imouthK(i)
      enddo
    
    
      iu_RD=20
      if (imt .eq. 32) then
         name="RDtoO.CS32"
         nameout="RDdistocean_CS32.bin"
      elseif (imt .eq. 90) then
         name="RDtoO.CS90"
         nameout="RDdistocean_CS90.bin"
      endif


      write(*,*) name,imt,jmt

c*    read i,j,k coordinates of downstream cells

         open( iu_RD, FILE=name,FORM='unformatted',
     &        STATUS='unknown')

         read(iu_RD) title,idown,jdown,kdown
         close(iu_RD)
      write(*,*) "read RDijk2ll_CS"
      endif

c*    form global arrays for mapping (i,j,k) -> (lon,lat)
      call pack_data(grid,lon2d_dg,lon_glob)
      call pack_data(grid,lat2d_dg,lat_glob)

      if (am_i_root()) then

c*    convert i,j,k coordinates of river mouths to absolute lat-lon coordinates    
      do k=1,nrvr
            lat_rvr(k)=lat_glob(imouthI(k),imouthJ(k)
     &            ,imouthK(k))
            lon_rvr(k)=lon_glob(imouthI(k),imouthJ(k)
     &            ,imouthK(k))
      enddo

c*    convert i,j,k coordinates of downstream cells to absolute lat-lon coordinates    
      do k=1,6
        do j=1,jmt
          do i=1,imt
            if (FOCEAN(i,j,k) .lt. 1.e-6) then
              down_lat(i,j,k)=lat_glob(idown(i,j,k),jdown(i,j,k)
     &            ,kdown(i,j,k))
              down_lon(i,j,k)=lon_glob(idown(i,j,k),jdown(i,j,k)
     &            ,kdown(i,j,k))

c*    dummy emergency directions
              down_lat_911(i,j,k)=down_lat(i,j,k)
              down_lon_911(i,j,k)=down_lon(i,j,k)
            write(130+k,200) lon_glob(i,j,k),lat_glob(i,j,k),
     &           down_lon(i,j,k)-lon_glob(i,j,k),
     &           down_lat(i,j,k)-lat_glob(i,j,k)
            else
              down_lat(i,j,k)=undef
              down_lon(i,j,k)=undef
              down_lat_911(i,j,k)=undef
              down_lon_911(i,j,k)=undef
            write(130+k,200) lon_glob(i,j,k),lat_glob(i,j,k),
     &           0.,0.
            end if

          enddo
        enddo
      enddo

c*    output everything 
      open(iu_RD,FILE=nameout,FORM='unformatted',
     &        STATUS='unknown')

      if (imt .eq. 32) then
         title1="river directions from dist. to ocean, CS32, April 09"
      elseif (imt .eq. 90) then
         title1="river directions from dist. to ocean, CS90, April 09"
      endif

      title2="Named River Mouths:"
      title3="Latitude of downstream river direction box"
      title4="Longitude of downstream river direction box"
      title5="Latitude of emergency downstream river direction box"
      title6="Longitude of emergency downstream river direction box"

      write(iu_RD) title1   
      write(iu_RD) title2,nrvr,namervr(1:nrvr),lat_rvr(1:nrvr)
     *     ,lon_rvr(1:nrvr)
      write(iu_RD) title3,down_lat
      write(iu_RD) title4,down_lon
      write(iu_RD) title5,down_lat_911
      write(iu_RD) title6,down_lon_911
      close(iu_RD)

      write(*,*) "wrote RD file"

      endif

 200  format(4(1X,f8.3))

      end program RDijk2ll_CS
c*
