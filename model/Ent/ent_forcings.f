      module ent_forcings
      implicit none

      private

      public open_forcings_file, close_forcings_file,  read_forcings

      integer I0_file, I1_file, J0_file, J1_file

      integer iu_SAT  ! - surface air temp
      integer iu_QLAT  ! - latent heat flux
      integer iu_QSEN  ! - sensible heat flux
      integer iu_PREC  ! - precip
      integer iu_MP1  ! - soil matric potential
      integer iu_MP2  ! - soil matric potential
      integer iu_MP3  ! - soil matric potential
      integer iu_MP4  ! - soil matric potential
      integer iu_MP5  ! - soil matric potential
      integer iu_MP6  ! - soil matric potential
      integer iu_SIC1  ! - soil ice
      integer iu_SIC2  ! - soil ice
      integer iu_SIC3  ! - soil ice
      integer iu_SIC4  ! - soil ice
      integer iu_SIC5  ! - soil ice
      integer iu_SIC6  ! - soil ice
      integer iu_VS  ! -  V wind speed
      integer iu_US  ! - U wind speed
      integer iu_CH  ! - heat transf coeff
      integer iu_PS  ! - surface pressure
      integer iu_TP  ! - canopy temperature
      integer iu_DVIS  ! - direct vis
      integer iu_VIS  ! - total vis
      integer iu_QFOL  ! - Qf - foliage surf humidity
      integer iu_ZEN  ! - cos zenith angle
! soil temperature and moisture now from all 6 GISS layers (0-3.5m) -PK 11/29/06
!***note that soil "moisture" is in m, so divide by layer depth to get volum frac*** -PK 
      integer iu_W1   ! - soil water (m), layer 1 (0-0.1 m)
      integer iu_W2   ! - soil water, layer 2 (0.1-0.27 m)
      integer iu_W3   ! - soil water, layer 3 (0.27-0.57 m)
      integer iu_W4   ! - soil water, layer 4 (0.57-1.08 m)
      integer iu_W5   ! - soil water, layer 5 (1.08-1.97 m)
      integer iu_W6   ! - soil water, layer 6 (1.97-3.50 m)
      integer iu_SOILT1   ! - soil temp (C), layer 1
      integer iu_SOILT2   ! - soil temp, layer 2
      integer iu_SOILT3   ! - soil temp, layer 3
      integer iu_SOILT4   ! - soil temp, layer 4
      integer iu_SOILT5   ! - soil temp, layer 5
      integer iu_SOILT6   ! - soil temp, layer 6
!      integer iu_SOILT30cm  !soil temp avg top 30 cm
!      integer iu_SMOIST30cm  !soil volum moist avg top 30 cm

      integer iu_LAI1
      integer iu_LAI2
      integer iu_LAI3
      integer iu_LAI4
      integer iu_LAI5
      integer iu_LAI6
      integer iu_LAI7
      integer iu_LAI8

      integer iu_height1
      integer iu_height2
      integer iu_height3
      integer iu_height4
      integer iu_height5
      integer iu_height6
      integer iu_height7
      integer iu_height8

      contains

      subroutine open_forcings_file( i0f, i1f, j0f, j1f, force_VEG,skip)
      use filemanager
      integer i0f, i1f, j0f, j1f
      logical force_VEG
      integer skip !#HACK, number of records to skip at beginning of forcing file.
      
      ! These indices could be provided in input file. They describe
      ! the dimensions of arrays in input files and should be compatible
      ! with M, JM, I0, I1, J0, J1 in the main program

!      I0_file = 1
!      I1_file = 72
!      J0_file = 1
!      J1_file = 46
!      I0_file = 17 !Ponca
!      I1_file = 17
!      J0_file = 33
!      J1_file = 33

      I0_file = i0f
      I1_file = i1f
      J0_file = j0f
      J1_file = j1f


      call openunit("CH",iu_CH,.true.,.true.)
      call openunit("DVIS",iu_DVIS,.true.,.true.)
      call openunit("MP1",iu_MP1,.true.,.true.)
      call openunit("MP2",iu_MP2,.true.,.true.)
      call openunit("MP3",iu_MP3,.true.,.true.)
      call openunit("MP4",iu_MP4,.true.,.true.)
      call openunit("MP5",iu_MP5,.true.,.true.)
      call openunit("MP6",iu_MP6,.true.,.true.)
      call openunit("PREC",iu_PREC,.true.,.true.)
      call openunit("PS",iu_PS,.true.,.true.)
      call openunit("QFOL",iu_QFOL,.true.,.true.)
      call openunit("QLAT",iu_QLAT,.true.,.true.)
      call openunit("QSEN",iu_QSEN,.true.,.true.)
      call openunit("SAT",iu_SAT,.true.,.true.)
      call openunit("SIC1",iu_SIC1,.true.,.true.)
      call openunit("SIC2",iu_SIC2,.true.,.true.)
      call openunit("SIC3",iu_SIC3,.true.,.true.)
      call openunit("SIC4",iu_SIC4,.true.,.true.)
      call openunit("SIC5",iu_SIC5,.true.,.true.)
      call openunit("SIC6",iu_SIC6,.true.,.true.)
      call openunit("TP",iu_TP,.true.,.true.)
      call openunit("US",iu_US,.true.,.true.)
      call openunit("VIS",iu_VIS,.true.,.true.)
      call openunit("VS",iu_VS,.true.,.true.)
      call openunit("ZEN",iu_ZEN,.true.,.true.)
      call openunit("W1",iu_W1,.true.,.true.)
      call openunit("W2",iu_W2,.true.,.true.)
      call openunit("W3",iu_W3,.true.,.true.)
      call openunit("W4",iu_W4,.true.,.true.)
      call openunit("W5",iu_W5,.true.,.true.)
      call openunit("W6",iu_W6,.true.,.true.)
      call openunit("SOILT1",iu_SOILT1,.true.,.true.)
      call openunit("SOILT2",iu_SOILT2,.true.,.true.)
      call openunit("SOILT3",iu_SOILT3,.true.,.true.)
      call openunit("SOILT4",iu_SOILT4,.true.,.true.)
      call openunit("SOILT5",iu_SOILT5,.true.,.true.)
      call openunit("SOILT6",iu_SOILT6,.true.,.true.)
!      call openunit("SOILT30cm",iu_SOILT30cm,.true.,.true.)
!      call openunit("SMOIST30cm",iu_SMOIST30cm,.true.,.true.)

      if ( force_VEG ) then
        call openunit("LAI1",iu_LAI1,.true.,.true.)
        call openunit("LAI2",iu_LAI2,.true.,.true.)
        call openunit("LAI3",iu_LAI3,.true.,.true.)
        call openunit("LAI4",iu_LAI4,.true.,.true.)
        call openunit("LAI5",iu_LAI5,.true.,.true.)
        call openunit("LAI6",iu_LAI6,.true.,.true.)
        call openunit("LAI7",iu_LAI7,.true.,.true.)
        call openunit("LAI8",iu_LAI8,.true.,.true.)
      endif

      if ( force_VEG ) then
        call openunit("HEIGHT1",iu_height1,.true.,.true.)
        call openunit("HEIGHT2",iu_height2,.true.,.true.)
        call openunit("HEIGHT3",iu_height3,.true.,.true.)
        call openunit("HEIGHT4",iu_height4,.true.,.true.)
        call openunit("HEIGHT5",iu_height5,.true.,.true.)
        call openunit("HEIGHT6",iu_height6,.true.,.true.)
        call openunit("HEIGHT7",iu_height7,.true.,.true.)
        call openunit("HEIGHT8",iu_height8,.true.,.true.)
      endif

      ! skip some records if needed

      call skip_forcings( skip, force_VEG )

      end subroutine open_forcings_file


      subroutine close_forcings_file( force_VEG )
      use filemanager
      logical force_VEG

      call closeunit(iu_CH)
      call closeunit(iu_DVIS)
      call closeunit(iu_MP1)
      call closeunit(iu_MP2)
      call closeunit(iu_MP3)
      call closeunit(iu_MP4)
      call closeunit(iu_MP5)
      call closeunit(iu_MP6)
      call closeunit(iu_PREC)
      call closeunit(iu_PS)
      call closeunit(iu_QFOL)
      call closeunit(iu_QLAT)
      call closeunit(iu_QSEN)
      call closeunit(iu_SAT)
      call closeunit(iu_SIC1)
      call closeunit(iu_SIC2)
      call closeunit(iu_SIC3)
      call closeunit(iu_SIC4)
      call closeunit(iu_SIC5)
      call closeunit(iu_SIC6)
      call closeunit(iu_TP)
      call closeunit(iu_US)
      call closeunit(iu_VIS)
      call closeunit(iu_VS)
      call closeunit(iu_ZEN)
      call closeunit(iu_W1)
      call closeunit(iu_W2)
      call closeunit(iu_W3)
      call closeunit(iu_W4)
      call closeunit(iu_W5)
      call closeunit(iu_W6)
      call closeunit(iu_SOILT1)
      call closeunit(iu_SOILT2)
      call closeunit(iu_SOILT3)
      call closeunit(iu_SOILT4)
      call closeunit(iu_SOILT5)
      call closeunit(iu_SOILT6)
!      call closeunit(iu_SOILT30cm)
!      call closeunit(iu_SMOIST30cm)

       if ( force_VEG ) then
        call closeunit(iu_LAI1)
        call closeunit(iu_LAI2)
        call closeunit(iu_LAI3)
        call closeunit(iu_LAI4)
        call closeunit(iu_LAI5)
        call closeunit(iu_LAI6)
        call closeunit(iu_LAI7)
        call closeunit(iu_LAI8)
       endif

       if ( force_VEG ) then
        call closeunit(iu_height1)
        call closeunit(iu_height2)
        call closeunit(iu_height3)
        call closeunit(iu_height4)
        call closeunit(iu_height5)
        call closeunit(iu_height6)
        call closeunit(iu_height7)
        call closeunit(iu_height8)
       endif
     
      end subroutine close_forcings_file


      subroutine skip_forcings( nskip, force_VEG )
      integer, intent(in) :: nskip
      logical, intent(in) :: force_VEG
      integer i
      
      if (.not.force_VEG) then
      do i=1,nskip
        read(iu_CH)
        read(iu_DVIS)
        read(iu_MP1)
        read(iu_MP2)
        read(iu_MP3)
        read(iu_MP4)
        read(iu_MP5)
        read(iu_MP6)
        read(iu_PREC)
        read(iu_PS)
        read(iu_QFOL)
        read(iu_QLAT)
        read(iu_QSEN)
        read(iu_SAT)
        read(iu_SIC1)
        read(iu_SIC2)
        read(iu_SIC3)
        read(iu_SIC4)
        read(iu_SIC5)
        read(iu_SIC6)
        read(iu_TP)
        read(iu_US)
        read(iu_VIS)
        read(iu_VS)
        read(iu_ZEN)
        read(iu_W1)
        read(iu_W2)
        read(iu_W3)
        read(iu_W4)
        read(iu_W5)
        read(iu_W6)
        read(iu_SOILT1)
        read(iu_SOILT2)
        read(iu_SOILT3)
        read(iu_SOILT4)
        read(iu_SOILT5)
        read(iu_SOILT6)
!        read(iu_SOILT30cm)
!        read(iu_SMOIST30cm)
        if ( force_VEG ) then
          read(iu_LAI1)
          read(iu_LAI2)
          read(iu_LAI3)
          read(iu_LAI4)
          read(iu_LAI5)
          read(iu_LAI6)
          read(iu_LAI7)
          read(iu_LAI8)
        endif
        if ( force_VEG ) then
          read(iu_height1)
          read(iu_height2)
          read(iu_height3)
          read(iu_height4)
          read(iu_height5)
          read(iu_height6)
          read(iu_height7)
          read(iu_height8)
        endif
      enddo
      endif

      end subroutine skip_forcings


      subroutine read_forcings( I0, I1, J0, J1, N_DEPTH, N_CASA_LAYERS,  !added N_CASA_LAYERS -PK 7/23/07
     &     TairC, TcanopyC, Qf, P_mbar, Ca, Ch, U,
     &     IPARdif, IPARdir, CosZen, 
     &     Soiltemp, Soilmoist, Soilmp, fice, LAI, height, 
     &     force_VEG,do_rewind)
      integer, intent(in) :: I0, I1, J0, J1, N_DEPTH, N_CASA_LAYERS
      real*8, dimension(I0:I1,J0:J1) :: TairC, TcanopyC,Qf,P_mbar,Ca,
     &     Ch, U, IPARdif, IPARdir, CosZen  
      real*8, dimension(N_DEPTH,I0:I1,J0:J1) :: Soilmp, fice
      real*8, dimension(N_CASA_LAYERS,I0:I1,J0:J1) :: Soiltemp,Soilmoist  !consistent with other modules -PK
      real*8, dimension(:,:,:),pointer :: LAI, height
      logical :: force_VEG,do_rewind
      !----Local------
      real*4, dimension(I0_file:I1_file,J0_file:J1_file) :: buf, buf1
      integer niter, niter1, i, j


      ! no files for the following data. setting to defaults
      !Qf(I0:I1,J0:J1) = 0.d0
      Ca(I0:I1,J0:J1) = 0.0127609 !(350 ppm @STP in mol m-3)

!      read(iu_CH) niter, buf
!      print *,'Got here in read_forcings.CH'
!      niter1 = niter
!      Ch(I0:I1,J0:J1) = buf(I0:I1,J0:J1)
      read(iu_CH) niter, buf
      niter1 = niter
      Ch(I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_DVIS) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      IPARdir(I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_MP1) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soilmp(1,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_MP2) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soilmp(2,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_MP3) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soilmp(3,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_MP4) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soilmp(4,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_MP5) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soilmp(5,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_MP6) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soilmp(6,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_PREC) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      != buf(I0:I1,J0:J1)

      read(iu_PS) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      P_mbar(I0:I1,J0:J1) = buf(I0:I1,J0:J1)
!      write(669,*) P_mbar(17,33)

      read(iu_QFOL) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Qf(I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_QLAT) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      != buf(I0:I1,J0:J1)

      read(iu_QSEN) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      != buf(I0:I1,J0:J1)

      read(iu_SAT) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      TairC= buf(I0:I1,J0:J1)

      read(iu_SIC1) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      fice(1,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_SIC2) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      fice(2,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_SIC3) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      fice(3,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_SIC4) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      fice(4,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_SIC5) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      fice(5,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_SIC6) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      fice(6,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_TP) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      TcanopyC(I0:I1,J0:J1) = buf(I0:I1,J0:J1) 

      read(iu_US) niter, buf1(I0:I1,J0:J1)
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      != buf(I0:I1,J0:J1)

      read(iu_VS) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      != buf(I0:I1,J0:J1)
      !! really should save absolute value
      do j=J0,J1
        do i=I0,I1
          U(i,j) = sqrt( buf1(i,j)**2 + buf(i,j)**2 )
        enddo
      enddo

      read(iu_VIS) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      IPARdif(I0:I1,J0:J1) = buf(I0:I1,J0:J1) - IPARdir(I0:I1,J0:J1)

      read(iu_ZEN) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      CosZen(I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      !***temporary hack for soilmoist and soiltemp*** -PK 7/23/07 
      read(iu_W1) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soilmoist(1,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

#ifdef NCASA2
      read(iu_W2) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soilmoist(2,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
#endif

!      read(iu_W3) niter, buf
!      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
!      Soilmoist(3,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

!      read(iu_W4) niter, buf
!      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
!      Soilmoist(4,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

!      read(iu_W5) niter, buf
!      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
!      Soilmoist(5,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

!      read(iu_W6) niter, buf
!      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
!      Soilmoist(6,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      read(iu_SOILT1) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soiltemp(1,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

#ifdef NCASA2
      read(iu_SOILT2) niter, buf
      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
      Soiltemp(2,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
#endif

!      read(iu_SOILT3) niter, buf
!      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
!      Soiltemp(3,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

!      read(iu_SOILT4) niter, buf
!      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
!      Soiltemp(4,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

!      read(iu_SOILT5) niter, buf
!      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
!      Soiltemp(5,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

!      read(iu_SOILT6) niter, buf
!      if(niter.ne.niter1) call stop_model("read_forcings: bad rec#",255)
!      Soiltemp(6,I0:I1,J0:J1) = buf(I0:I1,J0:J1)

      if (force_VEG) then

      if ( associated( LAI ) ) then
        read(iu_LAI1) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(1,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_LAI2) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(2,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_LAI3) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(3,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_LAI4) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(4,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_LAI5) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(5,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_LAI6) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(6,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_LAI7) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(7,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_LAI8) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        LAI(8,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
      endif

      if ( associated( height ) ) then
        read(iu_height1) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(1,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_height2) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(2,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_height3) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(3,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_height4) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(4,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_height5) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(5,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_height6) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(6,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_height7) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(7,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
        read(iu_height8) niter, buf
        if(niter.ne.niter1) call stop_model("read_forcings:bad rec",255)
        height(8,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
      endif

      endif !force_VEG

      end subroutine read_forcings
      
      end module ent_forcings
