      module ent_prog_mod
!@sum Utilities for off-line run of Ent model.
!#define OFFLINE TRUE
!#define PRINT_DRIVERS

      use ent_mod
      use FILEMANAGER
      implicit none

      contains

!************************************************************************
      subroutine ent_read_state( cells )
!@sum read ent state from the file
      type(entcelltype_public), intent(out) :: cells(:,:)
      !---
      integer, parameter :: MAX_BUFFER=100000 ! need realistic estimate
      real*8 buffer(MAX_BUFFER)
      integer iu_entstate
      integer ic, jc, i, j

      ic = size(cells,1)
      jc = size(cells,2)

      call openunit('ent_state',iu_entstate,.true.,.true.)
      do j=1,jc
        do i=1,ic
          read(iu_entstate) buffer
          ! check length of buffer : if( buffer(1) > MAX_BUFFER ) ??
          call ent_cell_unpack(buffer, cells(i,j))
        enddo
      enddo
      call closeunit(iu_entstate)

      end subroutine ent_read_state

!************************************************************************
      subroutine ent_write_state( cells )
!@sum write ent state to the file
      type(entcelltype_public), intent(in) :: cells(:,:)
      !---
      real*8, pointer :: buffer(:)
      integer iu_entstate
      integer ic, jc, i, j

      ic = size(cells,1)
      jc = size(cells,2)

      call openunit('ent_state_new',iu_entstate,.true.,.false.)
      do j=1,jc
        do i=1,ic
          call ent_cell_pack(buffer, cells(i,j))
          write(iu_entstate) buffer
          deallocate(buffer)
        enddo
      enddo
      call closeunit(iu_entstate)

      end subroutine ent_write_state

!************************************************************************
      subroutine ent_prescr_init( cells,
     &     IM, JM, I0, I1, J0, J1, jday, year,do_soilinit )
      use ent_prescribed_drv, only : init_canopy_physical,prescr_vegdata
     &                              ,prescr_soilpools  !for prescribing soil C, N pools from external file -PK
      use ent_prescr_veg, only : prescr_calcconst
      use ent_const

      implicit none
      type(entcelltype_public), intent(out) :: cells(:,:)
      integer, intent(in) :: IM, JM, I0, I1, J0, J1, jday, year
      logical, intent(in) :: do_soilinit
      !---Local variables-----
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: vegdata !cohort
      real*8, dimension(N_BANDS,N_COVERTYPES,I0:I1,J0:J1) :: albedodata !patch, NOTE:snow
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: laidata  !cohort
      real*8, dimension(N_COVERTYPES) :: hdata    !cohort
      real*8, dimension(N_COVERTYPES) :: nmdata    !cohort
      real*8, dimension(N_COVERTYPES,N_DEPTH) :: rootprofdata !Root fraction of veg type.
      real*8, dimension(N_COVERTYPES) :: popdata !Dummy population density:  0-bare soil, 1-vegetated
      real*8, dimension(N_COVERTYPES) :: dbhdata !Diameter at breast height for woody veg.(cm)
      real*8, dimension(N_COVERTYPES) :: craddata !Crown radius (m)
      real*8, dimension(N_COVERTYPES,N_BPOOLS,I0:I1,J0:J1) :: cpooldata !Carbon pools in individuals
      integer, dimension(N_COVERTYPES) :: soildata ! soil types 1-bright 2-dark
      real*8, dimension(N_SOIL_TEXTURES,I0:I1,J0:J1) :: soil_texture
      real*8, dimension(I0:I1,J0:J1) :: Ci_ini,CNC_ini,Tcan_ini,Qf_ini
      real*8, dimension(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS,
     &                  I0:I1,J0:J1):: Tpooldata  !in g/m2 -PK
        
      !hemi(J0:JM/2)   = -1    ! S.
      !hemi(JM/2+1:J1) =  1    ! N.

      !Read land surface parameters or use defaults
      !prescr data sets:
      call prescr_calcconst()
 
      call prescr_vegdata(jday, year, 
     &     IM,JM,I0,I1,J0,J1,vegdata,albedodata,laidata,hdata,nmdata,
     &     popdata,dbhdata,craddata,cpooldata,rootprofdata,
     &     soildata,soil_texture,Tpooldata,do_soilinit)
      !print *,"popdata in ent_GISS_init: ",popdata
      print *, "Tpooldata in ent_GISS_init: ",Tpooldata

      call init_canopy_physical(I0, I1, J0, J1,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini)


!#ifdef PRESCR_SOILCARB   !PK  -MOVED TO INSIDE prescr_vegdata
!      if (do_soilinit) then
!     &     call prescr_soilpools(IM,JM,I0,I1,J0,J1,Tpool_ini)
!      endif
!#endif                            

      print *,"soil_texture:",soil_texture
      print *,"vegdata:",vegdata
      !Translate gridded data to Entdata structure
      !GISS data:  a patch per vegetation cover fraction, one cohort per patch
      call ent_cell_set(cells, vegdata, popdata, laidata,
     &     hdata, dbhdata, craddata, cpooldata, nmdata, rootprofdata, 
     &     soildata, albedodata, soil_texture,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini, Tpooldata)  

      end subroutine ent_prescr_init

!************************************************************************
!************************************************************************

      end module ent_prog_mod
!************************************************************************




!************************************************************************
      program ent_prog
      ! this driver just passes the bounds of the region, so that all
      ! arrays can be allocated on a stack
      implicit none
      integer IM, JM, I0, I1, J0, J1, i0f, i1f, j0f, j1f
      integer jday, year, jday2, year2
      logical :: force_VEG
      logical :: do_soilinit,do_soilresp, do_phenology, do_patchdynamics
      logical :: do_spinup
      integer :: skip !#HACK to skip records at beginning of forcing file

      !* Default configuration
      force_VEG = .false.
      do_soilinit = .true.
      do_soilresp = .true.
      do_phenology = .false.
      do_patchdynamics = .false.
      do_spinup = .false.
      skip = 0

! Default values of IM, JM, I0, I1, J0, J1, jday, year, Ci_ini, CNC_ini, 
! Tcan_ini, Qf_ini should be set either here or inside of 
! "read_input_parameters"


      !* Set world bounds (should correspond to format of input files)
      IM = 72; JM = 46

      !* Set bounds of simulated region (i0=i1; j0=j1 for single cell)
      !* (all arrays will be allocated to this size)
      !* Full grid
      !i0=1; i1=IM
      !j0=1; j1=JM

      !* E.g. if one grid cell
      !i0 = 39; i1 = 39
      !j0 = 36; j1 = 36

      !* Single grid cell corresponding to Ponca lat/long
!      i0 = 17; i1 = 17
!      j0 = 33; j1 = 33

      !* Default, entire grid.
      i0 = 1; i1 = 72
      j0 = 1; j1 = 46

      !* dims to default forcings file
      i0f = 1; i1f = 72
      j0f = 1; j1f = 46

      !* Default date start for GISS GCM 10-day forcings
      jday = 152                !June 1
      year = 1980
      jday2 = jday + 10
      year2 = -1 !Initialize

      call read_input_parameters(IM, JM, I0, I1, J0, J1
     &     ,jday, year, jday2, year2
     &     ,i0f, i1f, j0f, j1f
     &     ,force_VEG
     &     ,do_soilinit, do_soilresp, do_phenology, do_patchdynamics
     &     ,skip, do_spinup)
      if (year2.eq.-1) year2=year  !If year2 was not specified, default same as year.

      print *,"ent_input: "
     &     , jday, year, jday2, year2
     &     , i0f, i1f, j0f, j1f
     &     , force_VEG
     &     , do_soilinit, do_soilresp, do_phenology, do_patchdynamics
     &     , do_spinup


      print *,"starting program"
      call run_offline(IM, JM, I0, I1, J0, J1, jday, year,jday2,year2
     &     ,i0f, i1f, j0f, j1f
     &     ,force_VEG
     &     ,do_soilinit, do_soilresp, do_phenology, do_patchdynamics
     &     ,skip,do_spinup)


      print *,"Ent run completed."
      end program ent_prog


!************************************************************************
      subroutine read_input_parameters(IM, JM, I0, I1, J0, J1,
     &     jday, year, jday2, year2, i0f, i1f, j0f, j1f, 
     &     force_VEG,do_soilinit,
     &     do_soilresp, do_phenology, do_patchdynamics,skip,do_spinup)

      use filemanager
      implicit none
      integer IM, JM, I0, I1, J0, J1, i0f, i1f, j0f, j1f
      integer jday, year,jday2, year2
      logical force_VEG
      logical do_soilinit,do_soilresp, do_phenology, do_patchdynamics
      logical do_spinup
      integer skip !#HACK to skip records at beginning of forcing file
      !----
      integer iu_ent_input
      character*80, parameter :: file_ent_input="ent_input"
      namelist /input_parameters/ IM, JM, I0, I1, J0, J1
     &     ,jday, year, jday2, year2
     &     ,i0f, i1f, j0f, j1f, force_VEG, do_soilinit
     &     ,do_soilresp, do_phenology, do_patchdynamics,skip,do_spinup

      print *,"reading input parameters from ", file_ent_input
      call openunit(trim(file_ent_input),iu_ent_input,.false.,.true.)
      read(iu_ent_input, NML=input_parameters, ERR=10)
      call closeunit(iu_ent_input)

      return
 10   continue
      print *,"error reading namelist file:", file_ent_input
      stop 255
      end subroutine read_input_parameters


!************************************************************************
      subroutine run_offline(IM, JM, I0, I1, J0, J1
     &     ,jday, year, jday2, year2
     &     ,i0f, i1f, j0f, j1f
     &     ,force_VEG, do_soilinit
     &     ,do_soilresp, do_phenology, do_patchdynamics
     &     ,skip,do_spinup)
      !****************************************************************
      !* Example program to run Ent coupled to a GCM.
      !* This version assumes an explicit scheme for calculation of
      !* canopy conductance, photosynthesis, and temperature.
      !* - For parallelization, the GCM provides the grid bounds to Ent.
      !* - Ent initializes vegetation structure parameters within these bounds.
      !* - GCM provides initial surface meteorological state variables.
      !* - Ent initializes GCANOPY, Ci, Qf, zeroes GPP
      !* - Simulation loop:
      !*   a. GCM provides meteorological drivers and canopy temperature
      !*   b. Ent updates.
      !*   c. Program gets from Ent: GCANOPY, Ci, Qf, GPP, TRANS_SW, albedo
      !*   d. GCM updates Qf, Tcanopy, given GCANOPY, TRANS_SW, albedo
      !*      and runs tracers on C and N.
      !****************************************************************

      !GCM_coupler:  need to write these routines specific to GCM.
!      use ent_GCM_coupler, only: 
!     &     GCM__get_grid, GCM_get_time, GCM_getdrv_cell, GCM_EWB

      use ent_mod
      !use ent_prescrveg
      use ent_prog_mod
      use ent_forcings
      use filemanager
      use ent_const  !PK 7/07
      !use ent_const, only : JEQUATOR
      
      implicit none
      integer, intent(in) :: IM, JM, I0, I1, J0, J1, i0f, i1f, j0f, j1f
      integer, intent(in) :: jday, year,jday2,year2
      logical, intent(in) :: force_VEG, do_soilinit
      logical, intent(in) :: do_soilresp,do_phenology,do_patchdynamics
      integer, intent(in) :: skip !#HACK to skip records at beginning of forcing file
      logical, intent(in) :: do_spinup
      !---Local----
      integer :: jdaycount  !Only needed for prognostic phenology
      logical :: do_rewind  

      ! Run parameters - NOT NEEDED
!      integer, parameter :: BIOPHYSICS_ONLY = 1
!      integer, parameter :: GISSVEG_CASA = 2
!      integer, parameter :: RUN_TYPE = BIOPHYSICS_ONLY
!      integer, parameter :: RUN_TYPE = GISSVEG_CASA

      ! time settings
      real*8, parameter :: dt = 1800.d0 !
!      real*8, parameter :: max_time = dt * 2!run for x time steps
!      real*8, parameter :: max_time = dt * 48.*10. !run for n days
!      real*8, parameter :: max_time = dt * 48.*365. !run for 1 year
      real*8 :: max_time !In seconds
      real*8, parameter :: save_interval=dt*48.d0*10d0 !save every 10 days

      ! model-dependent aarray sizes (maybe should be USEs from ent_...)
      !integer,parameter :: NSOILLAYERS = 6
      !integer,parameter :: NALBBANDS = 6

      !Coupling variables
      !Cell-level summary values - CALCULATED BY GCM/EWB OR OFF-LINE FILE
      real*8, dimension(I0:I1,J0:J1) ::
     &     TairC                !Air temperature (Celsius) !KIM - for phenology
     &     ,TcanopyC            !Canopy temperature (Celsius)
     &     ,Qf                  !Foliage surface specif humidity (kg vapor/ kg air)
     &     ,P_mbar              !Atmospheric pressure (mb)
     &     ,Ca                  !@Atmos CO2 conc at surface height (mol/m3).
     &     ,Ch                  !Ground to surface heat transfer coefficient 
     &     ,U                   !Surface layer wind speed (m s-1)
      !Radiation may later be broken down into hyperspectral increments.
      ! in an array
     &     ,IPARdif             !Incident diffuse PAR (vis) (W m-2)
     &     ,IPARdir             !Incident direct PAR (vis) (W m-2)
     &     ,CosZen              !cos of solar zenith angle
!***soil temp, moist now explicitly depth-structured*** -PK 7/07
!     &     ,Soiltemp            !soil temp avg top 30 cm (C) -PK 6/27/06
!     &     ,Soilmoist           !soil vol moist avg top 30 cm
      real*8, dimension(N_CASA_LAYERS,I0:I1,J0:J1) ::
     &      Soiltemp           !soil temp 
     &     ,Soilmoist          !soil volum moist 
      real*8, dimension(N_DEPTH,I0:I1,J0:J1) ::
!     &     SoilmoistGCM         !May be an array by depth (units TBA) !omit--see ent_forcings -PK 7/23/07
!     &     ,SoiltempGCM
     &     Soilmp              !Soil matric potential
     &     ,fice                !Fraction of soil water that is ice.
      real*8, dimension(:,:,:), pointer ::
     &     LAI                  ! prescribed LAI (for all pft's in the cell)
     &     ,height              ! prescribed height (for all pft's in the cell)
!      real*8, dimension(:,:,:), pointer ::
!     &     heightm              ! prescribed height (for all pft's in the cell)
      

      !Coupling and diagnostic variables specific to Ent (output)
      real*8, dimension(I0:I1,J0:J1) ::
     &     GCANOPY              !Canopy conductance of water vapor (m s-1). 
     &     ,Ci                  !Internal foliage CO2 concentration (mol/m3)
     &     ,GPP                 !Gross primary productivity (kg[C]/m2/s).
     &     ,TRANS_SW            !Transmittance of shortwave through canopy to soil 
     &     ,z0, CO2flux         !Variables from Ent to GCM
     &     ,C_labile, R_auto    !Variables from Ent to GCM
      !un-comment next one if want to export (see get_ent_exports below) -PK
!     &     ,Soil_resp           !soil respiration, patch level (added for ent_get_exports) -PK
      real*8, dimension(N_DEPTH,I0:I1,J0:J1) ::
     &     betadl
      real*8, dimension(N_BANDS,I0:I1,J0:J1) ::
     &     albedo               !Variables from Ent to GCM
      !un-comment next one if want to export (see get_ent_exports below) -PK
!      real*8, dimension(PTRACE,NPOOLS,N_CASA_LAYERS,I0:I1,J0:J1) ::
!     &     Tpool  !now explicitly depth-structured -PK 7/07
      !---------------------------------------------------------------------

      type(entcelltype_public) cells(I0:I1,J0:J1)
      !integer iu_forcings
      integer iu_results
      real*8 time, time_since_last_save
      integer hemi(I0:I1,J0:J1) !hemisphere flags = 1 for N., =-1 for S.
      logical :: dailyupdate    !For prescribed phenology, litter

      print *,"started run_offline with:", IM, JM, I0, I1, J0, J1

      ! now allocate optional arrays
      if ( force_VEG ) then
        allocate( LAI(N_PFT,I0:I1,J0:J1) )
        allocate( height(N_PFT,I0:I1,J0:J1) )
        !allocate( heightm(N_PFT,I0:I1,J0:J1) )
      else
        nullify( LAI )
        nullify( height )
        !nullify(heightm)
      endif

      ! now Ent should be initialized before any other calls to ent_*
      call ent_initialize(
     &     do_soilresp=do_soilresp
     &     ,do_phenology=do_phenology
     &     ,do_patchdynamics=do_patchdynamics)
 
      !* Set hemisphere flags.
      if ( J0<=JM/2 )   hemi(:,J0:min(JM/2,J1))   = -1    ! S.
      if ( J1>=JM/2+1 ) hemi(:,max(JM/2+1,J0):J1) =  1    ! N.
      print *,"set hemi ok"
      
      !* Initialize ent cells.
      call ent_cell_construct( cells )
      print *,"ent_cell_construct passed ok" 
      !if binary intialization file exists
      ! call ent_read_state( cells )
      !else initialize from scratch
      !jday = 152  !June 1
      !year = 1980


      call ent_prescr_init( cells, IM, JM, I0, I1, J0, J1, jday, year, 
     &     do_soilinit)
      print *,"ent_prescr_init passed ok"
      !print *,"Printing entcells from run_offline"
      !call ent_cell_print( cells )
      !print *,"End print from run_offline"
    
      ! open file with forcings (one record per time step)
      !!call openunit('ent_forcings',iu_forcings,.true.,.true.)
      call open_forcings_file( i0f, i1f, j0f, j1f, force_VEG,skip )


      ! open file for ent results (to save on each time step)
      call openunit('ent_results',iu_results,.true.,.false.)

#ifdef OFFLINE
#ifdef PRINT_DRIVERS
      write(980,*) "I0 I1 J0 J1 N_DEPTH TairC  TcanopyC Qf P_mbar Ca Ch U
     & IPARdif IPARdir CosZen St1 Sm1 
     & Smp1 Smp2 Smp3 Smp4 Smp5 Smp6
     & fice1 fice2 fice3 fice4 fice5 fice6 
     & LAI1 LAI2 LAI3 LAI4 LAI5 LAI6 LAI7 LAI8" 
!     & h1 h2 h3 h4 h5 h6 h7 h8" 

!!!! NEVER ever use explicit number for i/o unit (except 99) !!!
!!! or at least set it to a very big number to avoid conflicts with filemanager
      write(981,*) "GCANOPY Ci Qf GPP TRANS_SW z0 CO2flux 
     & betadl1 betadl2 betadl3 betadl4 betadl5 betadl6
     & albedo1 albedo2 albedo3 albedo4 albedo5 albedo6"
#endif
      !*********** DIAGNOSTICS FOR PLOTTING in ent.f *******************!
      if (N_CASA_LAYERS == 1) then  !labels depend on number of layers -PK 7/07
       write(995,*)"patchnum IPARdir IPARdif coszen pft n lai heightm
     & leaf froot wood surfmet surfstr soilmet                      
     & soilstr cwd surfmic soilmic slow passive
     & C_fol C_w C_froot C_root C_lab C_repro TRANS_SW
     & Ci GPP Rauto Soilresp NPP CO2flux GCANOPY senescefrac"
      else if (N_CASA_LAYERS == 2) then
       write(995,*)"patchnum IPARdir IPARdif coszen pft n lai heightm
     & leaf1 froot1 wood1 surfmet1 surfstr1 soilmet1                    !suffixes 1,2 for 1st, 2nd soil bgc layers -PK  
     & soilstr1 cwd1 surfmic1 soilmic1 slow1 passive1
     & leaf2 froot2 wood2 surfmet2 surfstr2 soilmet2  
     & soilstr2 cwd2 surfmic2 soilmic2 slow2 passive2 
     & C_fol C_w C_froot C_root C_lab C_repro TRANS_SW
     & Ci GPP Rauto Soilresp NPP CO2flux GCANOPY senescefrac"
      end if

!      write(996,*) "area IPARdir IPARdif LAI fv
!     & leaf froot wood surfmet surfstr soilmet                      
!     & soilstr cwd surfmic soilmic slow passive
!     & C_fol  C_w  C_froot  C_root  C_lab 
!     & TRANS_SW Ci  GPP R_auto Soil_resp NPP CO2flux GCANOPY"

!      write(998, *) "TcanopyC cop%GPP cop%R_root cop%R_auto cop%C_fol"

!#endif
!      write(982,*)"CiPa gasc tk vegpar:sigma sqrtexpr kdf rhor kbl alai 
!     &     nm vh Ntot vegalbedo ppar:pcp Kc Ko n1 n2 m1 m2 msat N0 Oi k"
#endif

      
      !* Start loop over time.
      max_time = dt*48*(jday2-jday+1)*(year2-year+1)
      print *,"max_time: ",max_time
      time = 0.d0
      time_since_last_save = 0.d0
      jdaycount = jday
      do_rewind=.false.
      do while( time < max_time )

        print *,"---------------------------------------------------"
        print *,"started time step, time=", time, "dt=", dt
!        write(99,*)"---------------------------------------------------"
!        write(99,*)"started time step, time=", time, "dt=", dt

#ifdef DEBUG        
      call assign_dummyvals(I0,I1,J0,J1,N_DEPTH,
     &  TairC,TcanopyC,Qf,P_mbar,Ca,Ch,U,
     &  IPARdif,IPARdir,CosZen, Soilmoist, Soilmp,fice)
#endif        

      print *,"Got here before read_forcings."
      call read_forcings(I0, I1, J0, J1, N_DEPTH, N_CASA_LAYERS, !added N_CASA_LAYERS -PK 7/23/07
     &       TairC, 
     &       TcanopyC, Qf, P_mbar, Ca, Ch, U,
     &       IPARdif, IPARdir, CosZen,
     &       Soiltemp, Soilmoist, Soilmp, fice, LAI, height, 
     &     force_VEG, do_rewind.and.do_spinup)
        ! hack to use upper layer data only:
        !**omit -- should be done completely separately from Ent** -PK 7/23/07
!        Soiltemp(:,:) = SoiltempGCM(1,:,:)
!        Soilmoist(:,:) = SoilmoistGCM(1,:,:)

        print *, 'call read_forcings ok'  !test -PK 7/11/06

#ifdef OFFLINE
#ifdef PRINT_DRIVERS
        write(980,*) I0, I1, J0, J1, N_DEPTH,
     &       TairC,
     &       TcanopyC, Qf, P_mbar, Ca, Ch, U,
     &       IPARdif, IPARdir, CosZen,
     &       Soiltemp, Soilmoist, Soilmp(:,I0:I1,J0:J1), 
     &       fice(:,I0:I1,J0:J1),
     &       LAI(:,I0:I1,J0:J1)
!     &       heightm(:,I0:I1,J0:J1)
#endif
#endif

        !* Set forcings.

        !#### IF STARTING FROM RESTART, NEED TO INIALIZE ##
        !####     Ci, CNC, Tcan, Qf                      ##

        call ent_set_forcings( cells,
     &       air_temperature=TairC,
     &       canopy_temperature=TcanopyC,
     &       canopy_air_humidity=Qf,
     &       surf_pressure=P_mbar,
     &       surf_CO2=Ca,
     &       heat_transfer_coef=Ch,
     &       wind_speed=U,
     &       total_visible_rad=IPARdir+IPARdif,
     &       direct_visible_rad=IPARdir,
     &       cos_solar_zenith_angle=CosZen,
     &       soil_temp=Soiltemp,  
     &       soil_moist=Soilmoist,
     &       soil_matric_pot=Soilmp,
     &       soil_ice_fraction=fice
     &       )
        print *, 'call ent_set_forcings ok'  !test -PK 7/11/06

      if (mod(time,86400.d0) .EQ. 0.d0) then !temporarily, later use timestruct
         dailyupdate=.true.
      else 
         dailyupdate=.false.
      end if

      !* NEW STREAMLINED CONTROL *!
      if (dailyupdate) then
        if (.not.do_phenology) then
          print *, 'Got here, updating veg.'
          call ent_prescribe_vegupdate(cells,hemi,jdaycount,year,
     &         do_giss_phenology=(.not.do_phenology),
     &         do_giss_lai=(.not.force_VEG),
     &         update_crops=.false.,
     &         laidata=LAI,     ! pft LAI
     &         hdata=height )   ! pft height
        end if                  !.not.do_phenology
      endif                    !dailyupdate

      call ent_run( cells, dt, time) !Change name to ent_processes, calls ent_integrate
      print *, 'call ent_run ok' !test -PK 7/11/06
      !***************************!

      !************************************************************************
!      if (.false.) then
!#ifdef PFT_MODEL_ENT
!        call ent_run (cells, dt, time)
!#endif
!        if (RUN_TYPE.eq.BIOPHYSICS_ONLY) then
!          if (dailyupdate) then
!            if ( force_VEG ) then !if reading in external LAI  -PK 8/16/07 
!              call ent_prescribe_vegupdate(cells,hemi,jday,year,
!     &             do_giss_phenology=.true.,
!     &             do_giss_lai=(.not.force_VEG),
!     &             update_crops=.false.,
!     &             laidata=LAI, ! pft LAI
!     &             hdata=height ) ! pft height
!                                !jday = mod(jday,365) + 1  !Use this if jday is not read in with forcings.
!            else
!              call ent_prescribe_vegupdate(cells,hemi,jday,year,
!     &             do_giss_phenology=.true.,
!     &             do_giss_lai=(.not.force_VEG),
!     &             update_crops=.false.)
!                   !jday = mod(jday,365) + 1  !Use this if jday is not read in with forcings.
!            end if              !force_VEG
!          end if                !dailyupdate
!          call ent_fast_processes( cells, dt )
!        else if (RUN_TYPE.eq.GISSVEG_CASA) THEN
!          if (dailyupdate) then
!            if ( force_VEG ) then !if reading in external LAI  -PK 8/16/07 
!             call ent_prescribe_vegupdate(cells,hemi,jday,year,
!     &           do_giss_phenology=.true.,
!     &           do_giss_lai=(.not.force_VEG),
!     &           update_crops=.false.,
!     &           laidata=LAI,         ! pft LAI
!     &           hdata=height )       ! pft height
!                 !jday = mod(jday,365) + 1  !Use this if jday is not read in with forcings.
!            else
!                 call ent_prescribe_vegupdate(cells,hemi,jday,year,
!     &           do_giss_phenology=.true.,
!     &           do_giss_lai=(.not.force_VEG),
!     &           update_crops=.false.)
!                 !jday = mod(jday,365) + 1  !Use this if jday is not read in with forcings.
!            end if  !force_VEG
!           end if  !dailyupdate
!          !call ent_run( cells, dt)
!          call ent_fast_processes( cells, dt )
!          
!          print *, 'call ent_run ok' !test -PK 7/11/06
!          print *, 'time=',time
!        end if
!
!      end if !.false 
      !************************************************************************

#ifdef DEBUG
        print *,"Got here before ent_get_exports."
#endif
        !* Extract results from Ent structure for GCM or diagnostics.
        call ent_get_exports( cells,
     &       canopy_conductance=GCANOPY,
     &       beta_soil_layers=betadl,
     &       shortwave_transmit=TRANS_SW,
     &       leafinternal_CO2=Ci,
     &       foliage_humidity=Qf,
     &       canopy_gpp=GPP,
     &       roughness_length=z0,
     &       flux_CO2=CO2flux,
     &       C_labile=C_labile,
     &       R_auto=R_auto,
     &       albedo=albedo
        !un-comment next two to write from here rather than ent.f -PK
!     &       ,soilresp=Soil_resp
!     &       ,soilcpools=Tpool
     &     )

        !* Save results to a file.
        write(iu_results) GCANOPY,Ci,Qf,GPP,TRANS_SW,z0,CO2flux,
     &       betadl,albedo


#ifdef OFFLINE
#ifdef PRINT_DRIVERS
        write(981,'(100e16.6)') GCANOPY,Ci,Qf,GPP,TRANS_SW,z0,CO2flux,
     &       betadl,albedo

!        print *,"Got here after iu_results."
!        print *,"canopy_conductance=",GCANOPY(I0,J0)
!        print *,"beta_soil_layers=",betadl(:,I0,J0)
!        print *,"shortwave_transmit=",TRANS_SW(I0,J0)
!        print *,"leafinternal_CO2=",Ci(I0,J0)
!        print *,"foliage_humidity=",Qf(I0,J0)
!        print *,"canopy_gpp=",GPP(I0,J0)
!        print *,"roughness_length=",z0(I0,J0)
!        print *,"flux_CO2=",CO2flux(I0,J0)
!        print *,"albedo=",albedo(:,I0,J0)
#endif
#endif


        !* Write ent state to a restart file.
        time_since_last_save = time_since_last_save + dt
        if ( time_since_last_save > save_interval ) then
          call ent_write_state( cells )
          time_since_last_save = 0.d0
        endif

        if (dailyupdate) then
          jdaycount = jdaycount + 1
          if (jdaycount>365) then !This assumes forcing file is only 1 year long.
            do_rewind = .true.
            jdaycount = 1     !reset
          else
            do_rewind = .false.
          endif
          !jdaycount = mod(jdaycount,365) !mod(365) goes from 0 to 364 for days 365 and 1-364.  Use this if jday is not read in with forcings.
        endif

        time = time + dt
        
      enddo

      call closeunit(iu_results)
      !!call closeunit(iu_forcings)
      call close_forcings_file( force_VEG )

      end subroutine run_offline

!***************************************************************************
      

!***************************************************************************

      subroutine assign_dummyvals(I0,I1,J0,J1,N_DEPTH,
     &     TairC,
     &     TcanopyC,Qf,P_mbar,Ca,Ch,U,
     &     IPARdif,IPARdir,CosZen, Soilmoist, Soilmp,fice)
!! Dummy forcings (in the absence of real ones)
      integer, intent(in) :: I0,I1,J0,J1,N_DEPTH
      real*8,intent(out) :: TairC,TcanopyC,Qf,P_mbar,Ca,Ch,U
      real*8,intent(out) :: IPARdif,IPARdir,CosZen
      real*8, dimension(:,:,:) ::Soilmoist, Soilmp,fice

      TairC = 20.
      TcanopyC = 20.            ! TPJUN1950_datafile - canopy temperature
      Qf = 0.                   !Needs to be passed in from land surface.
      P_mbar = 900.             ! PSJUN1950_datafile - surface pressure
      Ca = 350.0e-6*(P_mbar*100)/(8.3145*298.15) !mol m-3
      Ch = 1.d-2                ! CHJUN1950_datafile - heat transf coeff
      U = 2.                    ! VSJUN1950_datafile -  V wind speed
                                ! USJUN1950_datafile - U wind speed
      IPARdif = 200.            ! VISJUN1950_datafile - total vis
      IPARdir = 50.             ! DVISJUN1950_datafile - direct vis
      CosZen = .5
      Soilmoist(:,:,:) = .01
      Soilmp(:,:,:) = -.1 ! MP1JUN1950_datafile - soil matric potential ???
                                ! should be 6 layers !
      fice(N_DEPTH,I0:I1,J0:J1) = 0. ! SIC1JUN1950_datafile - soil ice ????
                                ! should be 6 layers !
      end subroutine assign_dummyvals

!************************************************************************

  

