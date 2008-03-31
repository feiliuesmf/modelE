#include "rundeck_opts.h" 
      MODULE AEROSOL_SOURCES
!@sum repository for Koch aerosol sources, features, etc.
!@auth Dorothy Koch
      USE TRACER_COM
      IMPLICIT NONE
      SAVE
      INTEGER, PARAMETER :: ndmssrc  = 1
!@param lmAER maximum height for AEROCOM emissions (RESOLUTION DEPENDENT?)
      INTEGER, PARAMETER :: lmAER = 7   
!@param Laircr the number of layers of aircraft data read from file

      INTEGER, PARAMETER ::
     &                      Laircrs      =19
!@param aircrafts_Tyr1, aircrafts_Tyr2 the starting and ending years
!@+     for transient SO2 aircraft emissions (= means non transient)
      integer :: aircrafts_Tyr1=0,aircrafts_Tyr2=0      
!@var DMSinput           DMS ocean source (kg/s/m2)
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: DMSinput ! DMSinput(im,jm,12)
c!@var DMS_AER           DMS prescribed by AERONET (kg S/day/box)
      real*4, ALLOCATABLE, DIMENSION(:,:,:) :: DMS_AER  !(im,jm,366)
c!@var SS1_AER        SALT bin 1 prescribed by AERONET (kg S/day/box)
      real*4, ALLOCATABLE, DIMENSION(:,:,:) :: SS1_AER  !(im,jm,366)
c!@var SS2_AER        SALT bin 2 prescribed by AERONET (kg S/day/box)
      real*4, ALLOCATABLE, DIMENSION(:,:,:) :: SS2_AER  !(im,jm,366)
#ifdef EDGAR_1995
      INTEGER, PARAMETER :: nso2src  = 6
#else
      INTEGER, PARAMETER :: nso2src  = 1
#endif
!@var SO2_src    SO2 industry: surface source (kg/s/box)
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: SO2_src !(im,jm,nso2src)
!@var BCI_src    BC Industrial source (kg/s/box)
      real*8, ALLOCATABLE, DIMENSION(:,:) :: BCI_src !(im,jm)
!@var BCB_src    BC Biomass source (kg/s/box)
      real*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: BCB_src !(im,jm,lmAER,12)
!@var BCBt_src    BC Biomass source (kg/s/box)
      real*8, ALLOCATABLE, DIMENSION(:,:) :: BCBt_src !(im,jm)
!@var OCI_src    OC Industrial source (kg/s/box)
#ifdef TRACERS_OM_SP
      INTEGER, PARAMETER :: nomsrc  = 8
#else
      INTEGER, PARAMETER :: nomsrc = 1
#endif
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: OCI_src !(im,jm,nomsrc)
!@var OCT_src    OC Terpene source (kg/s/box)
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: OCT_src !(im,jm,12)
!@var OCB_src    OC Biomass source (kg/s/box)
      real*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: OCB_src !(im,jm,lmAER,12)
!@var OCBt_src OC trend Biomass source (kg/s/box)
      real*8, ALLOCATABLE, DIMENSION(:,:) :: OCBt_src  !(im,jm)
!@var hBC BC trend source (kg/s/box)
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: hbc  !(im,jm,2)
!@var hOC OC trend source (kg/s/box)
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: hoc  !(im,jm,2)
!@var hso2 so2 trend source (kg/s/box)
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: hso2  !(im,jm,2)
!@var BCII_src_3D  BCI aircraft source (kg/s)
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: BCI_src_3D !(im,jm,lm)
!@var ss_src  Seasalt sources in 2 bins (kg/s/m2)
      INTEGER, PARAMETER :: nsssrc = 2
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: ss_src !(im,jm,nsssrc)
      INTEGER, PARAMETER :: nso2src_3d  = 4
!@var SO2_src_3D SO2 volcanic, aircraft sources (and biomass) (kg/s)
      real*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: SO2_src_3D !(im,jm,lm,nso2src_3d)
!@var SO2_biosrc_3D SO2  biomass(kg/s)
      real*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: SO2_biosrc_3D !(im,jm,lmAER,12)
!@var SO2t_src SO2 biomass trend source
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: SO2t_src !(im,jm,lmAER)
!@var PBLH boundary layer height
!@var MDF is the mass of the downdraft flux
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: 
     *   oh,dho2,perj,tno3,o3_offline  !im,jm,lm
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: ohr,dho2r,perjr,
     *   tno3r,ohsr  !im,jm,lm,12   DMK jmon
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: craft  !(im,jm,lm)
      real*8, ALLOCATABLE, DIMENSION(:,:) :: snosiz
#ifdef TRACERS_RADON
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: rn_src
#endif
!var NH3_src_con Ammonium sources
      REAL*8, ALLOCATABLE, DIMENSION(:,:)       ::
     & NH3_src_con, NH3_src_cyc
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: hnh3_con,hnh3_cyc
!var off_HNO3 off-line HNO3 field, used for nitrate and AMP when gas phase chemistry turned off
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)     ::  off_HNO3, off_SS

      END MODULE AEROSOL_SOURCES

      MODULE LAKI_SOURCE
      IMPLICIT NONE
      SAVE
      INTEGER, DIMENSION(10), PARAMETER :: LAKI_MON = (/6,6,6,
     * 6,7,7,8,9,9,10/)
      INTEGER, DIMENSION(10), PARAMETER :: LAKI_DAY = (/8,11,14,
     * 27,9,29,31,7,26,25/)
      REAL*8, DIMENSION(10), PARAMETER :: LAKI_AMT_T = (/1.98,
     * 3.17,4.41,2.55,2.09,3.10,1.82,1.39,1.03,0.78/)
      REAL*8, DIMENSION(10), PARAMETER :: LAKI_AMT_S = (/8.42,
     * 13.53,18.79,10.85,8.91,13.20,7.78,5.91,4.37,3.32/)
      END MODULE LAKI_SOURCE
      
      subroutine alloc_aerosol_sources(grid)
!@auth D. Koch
      use domain_decomp, only: dist_grid, get
      use AEROSOL_SOURCES, only: DMSinput,DMS_AER,SS1_AER,SS2_AER,
     * nso2src,SO2_src,BCI_src,BCB_src,BCBt_src,lmAER,nomsrc,
     * OCI_src,OCT_src,OCB_src,OCBt_src,BCI_src_3D,imAER,
     * hbc,hoc,hso2,
     * nsssrc,ss_src,nso2src_3d,SO2_src_3D,SO2_biosrc_3D,
     * ohr,dho2r,perjr,
     * tno3r,oh,dho2,perj,tno3,ohsr
     * ,o3_offline
     * ,snosiz
     * ,craft,so2t_src, off_HNO3,off_SS
     * ,NH3_src_con, NH3_src_cyc
     * ,hnh3_con,hnh3_cyc
#ifdef TRACERS_RADON
     * ,rn_src
#endif

      use MODEL_COM, only: im,lm
      
      IMPLICIT NONE
      type (dist_grid), intent(in) :: grid
      integer ::  J_1H, J_0H
      integer :: IER
      logical :: init = .false.

      if(init)return
      init=.true.

      call get( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )

      allocate( craft(IM,J_0H:J_1H,LM)
     * ,STAT=IER )
      allocate( DMSinput(IM,J_0H:J_1H,12) ,STAT=IER)
c     if (imAER.eq.1) then  !don't know if this is OK
      allocate( DMS_AER(IM,J_0H:J_1H,366) ,STAT=IER)
      allocate( SS1_AER(IM,J_0H:J_1H,366) ,STAT=IER)
      allocate( SS2_AER(IM,J_0H:J_1H,366) ,STAT=IER) 
c     endif
      allocate( SO2_src(IM,J_0H:J_1H,nso2src) ,STAT=IER) 
      allocate( BCI_src(IM,J_0H:J_1H) ,STAT=IER)
      allocate( BCB_src(IM,J_0H:J_1H,lmAER,12),STAT=IER )
      allocate( hbc(IM,J_0H:J_1H,2),hoc(IM,J_0H:J_1H,2)
     *  ,hso2(IM,J_0H:J_1H,2) ,STAT=IER)
      allocate( BCBt_src(IM,J_0H:J_1H) ,STAT=IER)
      allocate( OCI_src(IM,J_0H:J_1H,nomsrc) ,STAT=IER)
      allocate( OCB_src(IM,J_0H:J_1H,lmAER,12) ,STAT=IER)
      allocate( OCT_src(IM,J_0H:J_1H,12) ,STAT=IER)
      allocate( OCBt_src(IM,J_0H:J_1H) ,STAT=IER)
      allocate( BCI_src_3D(IM,J_0H:J_1H,lm) ,STAT=IER)
      allocate( ss_src(IM,J_0H:J_1H,nsssrc) ,STAT=IER)
      allocate( SO2_src_3D(IM,J_0H:J_1H,lm,nso2src_3d),STAT=IER )
      allocate( SO2_biosrc_3D(IM,J_0H:J_1H,lmAER,12) ,STAT=IER)
      allocate( SO2t_src(IM,J_0H:J_1H,lmAER) ,STAT=IER)
      allocate( oh(IM,J_0H:J_1H,lm),dho2(IM,J_0H:J_1H,lm),
     * perj(IM,J_0H:J_1H,lm),tno3(IM,J_0H:J_1H,lm)
     * ,o3_offline(IM,J_0H:J_1H,lm),STAT=IER )
      allocate( ohr(IM,J_0H:J_1H,lm),dho2r(IM,J_0H:J_1H,lm),
     * perjr(IM,J_0H:J_1H,lm),tno3r(IM,J_0H:J_1H,lm),
     * ohsr(IM,J_0H:J_1H,lm),STAT=IER )
      allocate( snosiz(IM,J_0H:J_1H) ,STAT=IER)
#ifdef TRACERS_RADON
      allocate( rn_src(IM,J_0H:J_1H,12) ,STAT=IER)
#endif
c Nitrate aerosols
! I,J
      allocate(  NH3_src_con(IM,J_0H:J_1H) )
      allocate(  NH3_src_cyc(IM,J_0H:J_1H) )
      allocate( hnh3_cyc(IM,J_0H:J_1H,2), STAT=IER )
      allocate( hnh3_con(IM,J_0H:J_1H,2), STAT=IER )
c off line 
      allocate(  off_HNO3(IM,J_0H:J_1H,LM)     )
      allocate(  off_SS(IM,J_0H:J_1H,LM)     )

      return
      end subroutine alloc_aerosol_sources      
      subroutine get_O3_offline
!@sum read in ozone fields for aqueous oxidation
c
C**** GLOBAL parameters and variables:
C
      use model_com, only: jday,im,jm,lm
      use filemanager, only: openunit,closeunit
      use aerosol_sources, only: o3_offline
      use domain_decomp, only: grid, get, write_parallel
C     
      implicit none

C**** Local parameters and variables and arguments:
c
!@var nmons: number of monthly input files
      integer, parameter :: nmons=1,levo3=23
      integer, dimension(nmons) :: mon_units,imon
      integer :: i,j,k,l
      integer :: jdlast=0
      logical :: ifirst=.true.
      character*80 title
      character(len=300) :: out_line
      character*10 :: mon_files(nmons) = (/'O3_FIELD'/)
      logical :: mon_bins(nmons)=(/.true./) ! binary file?
      real*8 :: frac
      real*8, dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,levo3,1)
     &        :: src
      real*8, allocatable, dimension(:,:,:,:) :: tlca, tlcb
      save jdlast,mon_units,imon,ifirst,tlca,tlcb
      INTEGER :: J_1, J_0, J_0H, J_1H

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

      if (ifirst) then
        call GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
        allocate(tlca(im,j_0H:j_1H,levo3,nmons),
     &           tlcb(im,j_0H:j_1H,levo3,nmons))
        ifirst=.false.
      endif
      k=1
      call openunit(mon_files(k),mon_units(k),mon_bins(k),.true.)
      call read_mon3Dsources(levo3,mon_units(k),jdlast,
     & tlca(:,:,:,k),tlcb(:,:,:,k),src(:,:,:,k),frac,imon(k))     
      call closeunit(mon_units(k))
      
      jdlast = jday

      if(levo3 /= LM) then ! might be ok, but you should check
        write(out_line,*)
     &  'make sure levels are correct in get_O3_offline'
        call write_parallel(trim(out_line))
        call stop_model('check on get_O3_offline',255)
      endif
      
      o3_offline(:,J_0:J_1,:)=src(:,J_0:J_1,:,k)

      write(out_line,*)
     &'offline ozone interpolated to current day',frac
      call write_parallel(trim(out_line))
      
      return
      END SUBROUTINE get_O3_offline
      
      
      SUBROUTINE read_mon3Dsources(Ldim,iu,jdlast,tlca,tlcb,data1,
     & frac,imon)
! we need to combine this with the one that is used in TRACERS_
! SPECIAL_Shindell, (called read_monthly_3Dsources) and then move it
! into more general tracer code...
!@sum Read in monthly sources and interpolate to current day
!@+   Calling routine must have the lines:
!@+      real*8 tlca(im,jm,Ldim,nm),tlcb(im,jm,Ldim,nm)
!@+      integer imon(nm)   ! nm=number of files that will be read
!@+      data jdlast /0/
!@+      save jdlast,tlca,tlcb,imon
!@+   Input: iu, the fileUnit#; jdlast
!@+   Output: interpolated data array + two monthly data arrays
!@auth Jean Lerner and others / Greg Faluvegi
      USE MODEL_COM, only: jday,im,jm,idofm=>JDmidOfM
      USE FILEMANAGER, only : NAMEUNIT
      USE DOMAIN_DECOMP, only : GRID,GET,READT_PARALLEL,REWIND_PARALLEL
     & ,write_parallel
      implicit none
!@var Ldim how many vertical levels in the read-in file?
!@var L dummy vertical loop variable
      integer Ldim,L,imon,iu,jdlast
      character(len=300) :: out_line
      real*8 :: frac
      real*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     A2D,B2D,dummy
      real*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,Ldim) ::
     &     tlca,tlcb,data1

      integer :: J_0, J_1

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)
C
      if (jdlast == 0) then   ! NEED TO READ IN FIRST MONTH OF DATA
        imon=1                ! imon=January
        if (jday <= 16)  then ! JDAY in Jan 1-15, first month is Dec
          do L=1,Ldim*11
            CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),0,dummy,1)
          end do
          DO L=1,Ldim
            CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),0,A2D,1)
            tlca(:,J_0:J_1,L)=A2D(:,J_0:J_1)
          END DO
          CALL REWIND_PARALLEL( iu )
        else              ! JDAY is in Jan 16 to Dec 16, get first month
  120     imon=imon+1
          if (jday > idofm(imon) .AND. imon <= 12) go to 120
          do L=1,Ldim*(imon-2)
            CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),0,dummy,1)
          end do
          DO L=1,Ldim
            CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),0,A2D,1)
            tlca(:,J_0:J_1,L)=A2D(:,J_0:J_1)
          END DO
          if (imon == 13)  CALL REWIND_PARALLEL( iu )
        end if
      else                         ! Do we need to read in second month?
        if (jday /= jdlast+1) then ! Check that data is read in daily
          if (jday /= 1 .OR. jdlast /= 365) then
            write(out_line,*)'Bad day values in read_monthly_3Dsources'
     &      //': JDAY,JDLAST=',JDAY,JDLAST
            call write_parallel(trim(out_line),crit=.true.)
            call stop_model('Bad values in read_monthly_3Dsources',255)
          end if
          imon=imon-12             ! New year
          go to 130
        end if
        if (jday <= idofm(imon)) go to 130
        imon=imon+1                ! read in new month of data
        if (imon == 13) then
          CALL REWIND_PARALLEL( iu  )
        else
          do L=1,Ldim*(imon-1)
            CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),0,dummy,1)
          end do
        endif
      end if
      DO L=1,Ldim
        CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),0,B2D,1)
        tlcb(:,J_0:J_1,L)=B2D(:,J_0:J_1)
      END DO
 130  continue
c**** Interpolate two months of data to current day
      frac = float(idofm(imon)-jday)/(idofm(imon)-idofm(imon-1))
      data1(:,J_0:J_1,:) =
     & tlca(:,J_0:J_1,:)*frac + tlcb(:,J_0:J_1,:)*(1.-frac)
      return
      end subroutine read_mon3Dsources

      SUBROUTINE READ_OFFHNO3(OUT)
      USE MODEL_COM, only : im,jm,lm,jdate,JDendOFM,jmon
      USE DOMAIN_DECOMP, only : grid,unpack_data,am_i_root
      IMPLICIT NONE
      include 'netcdf.inc'
!@param  nlevnc vertical levels of off-line data  
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM),
     &    intent(out) :: OUT
      INTEGER, PARAMETER :: nlevnc =23
      REAL*4, DIMENSION(IM,JM,nlevnc) :: IN1_glob4, IN2_glob4
      REAL*8, DIMENSION(IM,JM,nlevnc) :: IN1_glob, IN2_glob
      REAL*8, DIMENSION(:,:,:), pointer, save :: IN1, IN2
!@var netcdf integer
      INTEGER :: ncid,id
      INTEGER, save :: step_rea=0, first_call=1
!@var time interpoltation
      REAL*8 :: tau
      integer start(4),count(4),status,l
c -----------------------------------------------------------------
c   Initialisation of the files to be read
c ----------------------------------------------------------------     

      if (first_call==1) then
        first_call=0
        allocate( IN1(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,nlevnc) )
        allocate( IN2(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,nlevnc) )
      endif
      if (step_rea.ne.jmon) then 
        step_rea = JMON
        if ( am_i_root() ) then
          print*,'READING HNO3 OFFLINE ',jmon, step_rea
c -----------------------------------------------------------------
c   Opening of the files to be read
c -----------------------------------------------------------------
          status=NF_OPEN('OFFLINE_HNO3.nc',NCNOWRIT,ncid)
          status=NF_INQ_VARID(ncid,'FELD',id)
C------------------------------------------------------------------
c -----------------------------------------------------------------
c   read
c -----------------------------------------------------------------
          start(1)=1
          start(2)=1
          start(3)=1
          start(4)=step_rea

          count(1)=im
          count(2)=jm
          count(3)=nlevnc
          count(4)=1

          status=NF_GET_VARA_REAL(ncid,id,start,count,IN1_glob4)
          start(4)=step_rea+1
          if (start(4).gt.12) start(4)=1
          status=NF_GET_VARA_REAL(ncid,id,start,count,IN2_glob4)

          status=NF_CLOSE('OFFLINE_HNO3.nc',NCNOWRIT,ncid)

          IN1_glob = IN1_glob4
          IN2_glob = IN2_glob4
        endif ! am_i_root
        call UNPACK_DATA(grid, IN1_glob, IN1)
        call UNPACK_DATA(grid, IN2_glob, IN2)
      endif

C-----------------------------------------------------------------------
      tau=(jdate-.5)/(JDendOFM(jmon)-JDendOFM(jmon-1))
         do l=1,lm
         OUT(:,:,l) = (1.-tau)*IN1(:,:,l)+tau*IN2(:,:,l)  
         enddo
c -----------------------------------------------------------------
      RETURN
      END SUBROUTINE READ_OFFHNO3
c -----------------------------------------------------------------

      SUBROUTINE READ_OFFSS(OUT)
      USE MODEL_COM, only : im,jm,lm,jdate,jmon,JDendOFM
      USE DOMAIN_DECOMP, only : grid,unpack_data,am_i_root
      IMPLICIT NONE
      include 'netcdf.inc'
!@param  nlevnc vertical levels of off-line data  
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM),
     &    intent(out) :: OUT
      INTEGER, PARAMETER :: nlevnc =23
      REAL*4, DIMENSION(IM,JM,nlevnc) :: IN1_glob4, IN2_glob4
      REAL*8, DIMENSION(IM,JM,nlevnc) :: IN1_glob, IN2_glob
      REAL*8, DIMENSION(:,:,:), pointer, save :: IN1_ss, IN2_ss
!@var netcdf integer
      INTEGER :: ncid,id
      INTEGER, save :: step_rea_ss=0, first_call_ss=1
!@var time interpoltation
      REAL*8 :: tau
      integer start(4),count(4),status,l,i,j
c -----------------------------------------------------------------
c   Initialisation of the files to be read
c ----------------------------------------------------------------     

      if (first_call_ss==1) then
        first_call_ss=0
        allocate( IN1_ss(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,nlevnc))
        allocate( IN2_ss(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,nlevnc))
      endif
      if (step_rea_ss.ne.jmon) then 
        step_rea_ss = JMON
        if ( am_i_root() ) then
          print*,'READING SEAS OFFLINE ',jmon, step_rea_ss
c -----------------------------------------------------------------
c   Opening of the files to be read
c -----------------------------------------------------------------
          status=NF_OPEN('OFFLINE_SEAS.nc',NCNOWRIT,ncid)
          status=NF_INQ_VARID(ncid,'FELD',id)
C------------------------------------------------------------------
c -----------------------------------------------------------------
c   read
c -----------------------------------------------------------------
          start(1)=1
          start(2)=1
          start(3)=1
          start(4)=step_rea_ss

          count(1)=im
          count(2)=jm
          count(3)=nlevnc
          count(4)=1

          status=NF_GET_VARA_REAL(ncid,id,start,count,IN1_glob4)
          start(4)=step_rea_ss+1
          if (start(4).gt.12) start(4)=1
          status=NF_GET_VARA_REAL(ncid,id,start,count,IN2_glob4)

          status=NF_CLOSE('OFFLINE_SEAS.nc',NCNOWRIT,ncid)

          IN1_glob = IN1_glob4
          IN2_glob = IN2_glob4
        endif ! am_i_root
        call UNPACK_DATA(grid, IN1_glob, IN1_ss)
        call UNPACK_DATA(grid, IN2_glob, IN2_ss)
      endif

C-----------------------------------------------------------------------
      tau=(jdate-.5)/(JDendOFM(jmon)-JDendOFM(jmon-1))
         do l=1,lm
         OUT(:,:,l) = (1.-tau)*IN1_ss(:,:,l)+tau*IN2_ss(:,:,l)  
         enddo
c -----------------------------------------------------------------
      RETURN
      END SUBROUTINE READ_OFFSS
c -----------------------------------------------------------------


      SUBROUTINE read_E95_SO2_source(nt,iact)
!@sum reads in SO2 surface sources: Edgar 1995 (RIVM)
C**** There are 2 monthly sources and 4 annual sources
C**** Annual sources are read in at start and re-start of run only
C**** Monthly sources are interpolated each day
      USE MODEL_COM, only: itime,jday,JDperY,im,jm
c    * ,DTsrc,ls1,jmon,t,lm
      USE DOMAIN_DECOMP, only : GRID, GET,readt_parallel 
     * ,AM_I_ROOT
      USE GEOM, only: BYDXYP
      USE CONSTANT, only: sday,hrday
      USE FILEMANAGER, only: openunit,closeunit,
     * nameunit
      USE TRACER_COM, only: itime_tr0,trname
      USE AEROSOL_SOURCES, only: src=>SO2_src,nsrc=>nso2src
c    *  ,SO2_biosrc_3D
      implicit none
      character*80 title
      logical :: ifirst=.true.
!@var nanns,nmons: number of annual and monthly input files
      integer, parameter :: nanns=4,nmons=2
      integer ann_units(nanns),mon_units(nmons)
      character*12 :: ann_files(nanns) =
     * (/'SO2_FFNB_E95','SO2_INNB_E95','SO2_WHNB_E95',
     *  'SO2_BFNB_E95'/)
      logical :: ann_bins=.true.
      character*12 :: mon_files(nmons) =
     * (/'SO2_AGNB_E95','SO2_BBNB_E95'/)
      logical :: mon_bins=.true.
      real*8, allocatable, dimension(:,:,:) :: tlca,tlcb  ! for monthly sources
      real*8 frac,bySperHr
      integer :: imon(nmons)
      integer i,j,jj,nt,iact,iu,k,j_0,j_1
      integer :: jdlast=0
      save ifirst,jdlast,tlca,tlcb,mon_units,imon

C Edgar-1995 has its own biomass (2D) source and diagnostic.
C To avoid complications with indexing of 3D sources, leave in the
C 3D biomass array and diagnostic but fill it with Edgar 1995
C emissions. Make all levels other than surface equal to zero (nbell)
       CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

c      so2_biosrc_3D(:,j_0:j_1,:,:)=0.
C Now the Edgar-1995 sources
C
C    K=1    Fossil fuel combustion           (edgar 95 annual)
C    K=2    industrial prod & cons processes (edgar 95 annual)
C    K=3    wastehandling (edgar 95 annual)
C    These 3 are from edgar 95, but with our monthly variation from
C    biomass burning applied to the annual values:
C    K=4    agricultural waste burning
C    K=5    biofuel production, transformation, combustion
C    K=6    biomass burning
C
      if (itime.lt.itime_tr0(nt)) return
      bySperHr = 1.d0/3600.d0

C****
C**** Annual Edgar-1995 Sources
C**** The EDGAR-1995 sources are in KGSO2/4x5grid/HR and need to be
C**** converted to KGSO2/m2/sec:
C****
      if (ifirst) then
       Allocate(tlca(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,nmons),
     &           tlcb(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,nmons))
c      allocate(tlca(im,j_0H:j_1H,nmons),tlcb(im,j_0H:j_1H,nmons))
        call openunit(ann_files,ann_units,ann_bins,.true.)
        k = 0
        do iu = 1,nanns
          k = k+1
c         call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          call readt_parallel (grid,
     &         ann_units(iu),nameunit(ann_units(iu)),0,src(:,:,k),1)
          do j=j_0,j_1
            src(:,j,k) = src(:,j,k)*bydxyp(j)*bySperHr
          end do
        end do
        call closeunit(ann_units)
       
        call openunit(mon_files,mon_units,mon_bins,.true.)
      endif
C****
C**** Monthly sources are interpolated to the current day
C**** The EDGAR-1995 sources are in KG/4x5grid/HR and need to be
C**** converted to KG/m2/sec:

C****
      ifirst = .false.
      jj = 0
      do k=nanns+1,nsrc
        jj = jj+1
        call read_monthly_sources(mon_units(jj),jdlast,
     *    tlca(:,:,jj),tlcb(:,:,jj),src(:,:,k),frac,imon(jj))
        do j=j_0,j_1
C Also fill 3D biomass array
C Units are kgSO2/s
c         if (k.eq.6) then
c         so2_biosrc_3D(:,j,1,jmon)= src(:,j,k)*bySperHr
c         endif
          src(:,j,k) = src(:,j,k)*bydxyp(j)*bySperHr
        end do
      end do
      jdlast = jday
      if (AM_I_ROOT())
     *  write(6,*) trname(nt),'Sources interpolated to current day',frac
      call sys_flush(6)
C****
      END SUBROUTINE read_E95_SO2_source

      SUBROUTINE get_hist_BM(iact)
c historic biomass: linear increase in tropics from 1/2 present day in 1875
      USE AEROSOL_SOURCES, only: BCB_src,OCB_src,BCBt_src,OCBt_src,lmAER
     * ,SO2_biosrc_3D,SO2t_src
      USE MODEL_COM, only: jyear,im,jm,jmon
      USE DOMAIN_DECOMP, only : GRID, GET
      USE GEOM, only: dxyp
      USE TRACER_COM, only: aer_int_yr,imPI,imAER
      USE FILEMANAGER, only: openunit,closeunit
      USE CONSTANT, only: sday
      IMPLICIT NONE
      integer ihyr,l,ii,jj,ll,iuc,mm,iact,j,mmm
      integer J_0,J_1,J_0H,J_1H
      real*8 carbstuff,tfac
      real*8, dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: 
     *  ratb,rato
      
      CALL GET(grid, J_STRT=J_0,       J_STOP=J_1,
     *               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

      if (iact.eq.0) then
      so2_biosrc_3D(:,:,:,:)= 0.d0
      call openunit('SO2_BIOMASS',iuc,.false.,.true.)
      do
      read(iuc,*) ii,jj,mm,ll,carbstuff
      if (ii.eq.0) exit
      if (jj<j_0 .or. jj>j_1) cycle
      if (imPI.eq.1) carbstuff=carbstuff*0.5d0
      SO2_biosrc_3D(ii,jj,ll,mm)=carbstuff/(sday*30.4d0)
      end do
      call closeunit(iuc)
      if (imAER.eq.1) go to 33
      BCB_src(:,:,:,:)=0.d0
      OCB_src(:,:,:,:)=0.d0
      call openunit('BC_BIOMASS',iuc,.false.,.true.)
      do 
      read(iuc,*) mm,ii,jj,carbstuff
      if (ii.eq.0) exit
      if (jj<j_0 .or. jj>j_1) cycle
      if (imPI.eq.1) carbstuff=carbstuff*0.5d0
      carbstuff=carbstuff*dxyp(jj)
      BCB_src(ii,jj,1,mm)=carbstuff
      end do
      call closeunit(iuc)
      if (imAER.eq.2) then
      if (aer_int_yr.lt.1990.or.aer_int_yr.gt.2010) then
      ratb(:,:)=0.d0
      call openunit('BC_BM_RAT',iuc,.false.,.true.)
      do 
      read(iuc,*) ii,jj,carbstuff
      if (ii.eq.0.) exit
      if (jj<j_0 .or. jj>j_1) cycle
      ratb(ii,jj)=carbstuff
      end do
      call closeunit(iuc)
      do mm=1,12
      BCB_src(:,J_0:j_1,1,mm)=BCB_src(:,j_0:j_1,1,mm)
     *  *ratb(:,j_0:j_1)
      end do
      endif
      endif  
      call openunit('OC_BIOMASS',iuc,.false.,.true.)
      do 
      read(iuc,*) mm,ii,jj,carbstuff
      if (ii.eq.0.) exit 
      if (jj<j_0 .or. jj>j_1) cycle
      if (imPI.eq.1) carbstuff=carbstuff*0.5d0
      carbstuff=carbstuff*dxyp(jj)
      OCB_src(ii,jj,1,mm)=carbstuff !this is OM
      end do
      call closeunit(iuc)
      if (imAER.eq.2) then
      if (aer_int_yr.lt.1990.or.aer_int_yr.gt.2010) then
      rato(:,:)=0.d0
      call openunit('OC_BM_RAT',iuc,.false.,.true.)
      do 
      read(iuc,*) ii,jj,carbstuff
      if (ii.eq.0.) exit 
      if (jj<j_0 .or. jj>j_1) cycle
      rato(ii,jj)=carbstuff
      end do
      call closeunit(iuc)
      do mm=1,12
      OCB_src(:,j_0:j_1,1,mm)=OCB_src(:,j_0:j_1,1,mm)
     * *rato(:,j_0:j_1)
      end do
      endif
      endif
      endif

       if (aer_int_yr.eq.0) then
       ihyr=jyear
       else
       ihyr=aer_int_yr
       endif

       BCBt_src(:,:)=BCB_src(:,:,1,jmon)
       OCBt_src(:,:)=OCB_src(:,:,1,jmon)
       SO2t_src(:,:,:)=SO2_biosrc_3D(:,:,:,jmon)
       if (ihyr.le.2000.and.imAER.eq.3) then
       tfac=0.5d0*(1.d0+real(ihyr-1875)/125.d0)
       do j=MAX(j_0,15),MIN(j_1,29)
       BCBt_src(:,j)=BCB_src(:,j,1,jmon)*tfac
       OCBt_src(:,j)=OCB_src(:,j,1,jmon)*tfac
       SO2t_src(:,j,:)=SO2_biosrc_3D(:,j,:,jmon)*tfac
       end do
       endif
 33   continue
      end subroutine get_hist_BM

      SUBROUTINE read_hist_BCOC(iact)
c historic BC emissions
      USE MODEL_COM, only: jyear, jday,im,jm
      USE FILEMANAGER, only: openunit,closeunit
      USE DOMAIN_DECOMP, only :  GRID, GET 
      USE TRACER_COM, only: aer_int_yr
      USE AEROSOL_SOURCES, only: BCI_src,OCI_src,nomsrc
     * ,hbc,hoc
      USE GEOM, only: dxyp
      implicit none
      integer iuc,irr,ihyr,i,j,id,jb1,jb2,iact,ii,jj,nn,
     * iy,ip, j_0,j_1,j_0h,j_1h
      real*8, dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,16) :: 
     * hOC_all,hBC_all
      real*8 d1,d2,d3,xbcff,xbcbm,xombm
      save jb1,jb2
      
      CALL GET(grid, J_STRT=J_0,       J_STOP=J_1,
     *             J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      
      BCI_src(:,:)=0.d0
      OCI_src(:,:,:)=0.d0

      if (aer_int_yr.eq.0) then
      ihyr=jyear
      else
      ihyr=aer_int_yr
      endif
c if run just starting or if it is a new year
c   then open new files
      if (iact.eq.0.or.jday.eq.1) then
      hbc(:,:,:)=0.0
      hoc(:,:,:)=0.0
      hoc_all(:,:,:)=0.0
      hbc_all(:,:,:)=0.0
      if (ihyr.ge.1850.and.ihyr.lt.1860) then
      jb1=1850
      jb2=1860
      irr=1
      else if (ihyr.ge.1860.and.ihyr.lt.1870) then
      jb1=1860
      jb2=1870
      irr=2
      else if (ihyr.ge.1870.and.ihyr.lt.1880) then
      jb1=1870
      jb2=1880
      irr=3
      else if (ihyr.ge.1880.and.ihyr.lt.1890) then
      jb1=1880
      jb2=1890
      irr=4
      else if (ihyr.ge.1890.and.ihyr.lt.1900) then
      jb1=1890
      jb2=1900
      irr=5
      else if (ihyr.ge.1900.and.ihyr.lt.1910) then
      jb1=1900
      jb2=1910
      irr=6
      else if (ihyr.ge.1910.and.ihyr.le.1920) then
      jb1=1910
      jb2=1920
      irr=7
      else if (ihyr.ge.1920.and.ihyr.le.1930) then
      jb1=1920
      jb2=1930
      irr=8
      else if (ihyr.ge.1930.and.ihyr.le.1940) then
      jb1=1930
      jb2=1940
      irr=9
      else if (ihyr.ge.1940.and.ihyr.le.1950) then
      jb1=1940
      jb2=1950
      irr=10
      else if (ihyr.ge.1950.and.ihyr.le.1960) then
      jb1=1950
      jb2=1960
      irr=11
      else if (ihyr.ge.1960.and.ihyr.le.1970) then
      jb1=1960
      jb2=1970
      irr=12
      else if (ihyr.ge.1970.and.ihyr.le.1980) then
      jb1=1970
      jb2=1980
      irr=13
      else if (ihyr.ge.1980.and.ihyr.le.1990) then
      jb1=1980
      jb2=1990
      irr=14
      else if (ihyr.ge.1990.and.ihyr.le.2000) then
      jb1=1990
      jb2=2000
      irr=15
      endif
      call openunit('BC_INDh',iuc, .false.,.true.)
      do ip=1,22754 !14040
      read(iuc,*) ii,jj,iy,xbcff
      if (jj<j_0 .or. jj>j_1) cycle
      hbc_all(ii,jj,iy)=xbcff
      end do
      call closeunit(iuc)
      call openunit('OC_INDh',iuc, .false.,.true.)
      do ip=1,19676 !14040
      read(iuc,*) ii,jj,iy,xbcff
      if (jj<j_0 .or. jj>j_1) cycle
      hoc_all(ii,jj,iy)=xbcff
      end do
      call closeunit(iuc)
      hbc(:,j_0:j_1,1:2)=hbc_all(:,j_0:j_1,irr:irr+1)
      hoc(:,j_0:j_1,1:2)=hoc_all(:,j_0:j_1,irr:irr+1)
c kg/year to kg/s
      hbc(:,j_0:j_1,1:2)=hbc(:,j_0:j_1,1:2)/(365.d0*24.d0*3600.d0)
c OM=1.4 x OC
      hoc(:,j_0:j_1,1:2)=hoc(:,j_0:j_1,1:2)/(365.d0*24.d0*3600.d0)*1.4d0
      endif
c interpolate to model year
c    
      d1=real(ihyr-jb1)
      d2=real(jb2-ihyr)
      d3=real(jb2-jb1)
      BCI_src(:,j_0:j_1)=(d1*hbc(:,j_0:j_1,2)
     * +d2*hbc(:,j_0:j_1,1))/d3
      OCI_src(:,j_0:j_1,1)=(d1*hoc(:,j_0:j_1,2)
     * +d2*hoc(:,j_0:j_1,1))/d3
  
      end subroutine read_hist_BCOC
      
      SUBROUTINE read_hist_SO2(iact)
c historic BC emissions
      USE MODEL_COM, only: jyear, jday,im,jm
      USE FILEMANAGER, only: openunit,closeunit
      USE DOMAIN_DECOMP, only :  GRID, GET 
      USE TRACER_COM, only: aer_int_yr
      USE AEROSOL_SOURCES, only: SO2_src,nomsrc
     * ,hso2
      USE GEOM, only: dxyp
      implicit none
      integer iuc,irr,ihyr,i,j,id,jb1,jb2,iact,ii,jj,nn,
     * iy,ip, j_0,j_1,j_0h,j_1h
      real*8, dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,12) :: 
     * hso2_all
      real*8 d1,d2,d3,xso2ff
      save jb1,jb2
      
      CALL GET(grid, J_STRT=J_0,       J_STOP=J_1,
     *             J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      
      SO2_src(:,:,:)=0.d0

      if (aer_int_yr.eq.0) then
      ihyr=jyear
      else
      ihyr=aer_int_yr
      endif
c if run just starting or if it is a new year
c   then open new files
      if (iact.eq.0.or.jday.eq.1) then
       hso2_all(:,:,:)=0.0
       hso2(:,:,:)=0.0
      if (ihyr.ge.1890.and.ihyr.lt.1900) then
      jb1=1890
      jb2=1900
      irr=1
      else if (ihyr.ge.1900.and.ihyr.lt.1910) then
      jb1=1900
      jb2=1910
      irr=2
      else if (ihyr.ge.1910.and.ihyr.le.1920) then
      jb1=1910
      jb2=1920
      irr=3
      else if (ihyr.ge.1920.and.ihyr.le.1930) then
      jb1=1920
      jb2=1930
      irr=4
      else if (ihyr.ge.1930.and.ihyr.le.1940) then
      jb1=1930
      jb2=1940
      irr=5
      else if (ihyr.ge.1940.and.ihyr.le.1950) then
      jb1=1940
      jb2=1950
      irr=6
      else if (ihyr.ge.1950.and.ihyr.le.1960) then
      jb1=1950
      jb2=1960
      irr=7
      else if (ihyr.ge.1960.and.ihyr.le.1970) then
      jb1=1960
      jb2=1970
      irr=8
      else if (ihyr.ge.1970.and.ihyr.le.1980) then
      jb1=1970
      jb2=1980
      irr=9
      else if (ihyr.ge.1980.and.ihyr.le.1990) then
      jb1=1980
      jb2=1990
      irr=10
      else if (ihyr.ge.1990.and.ihyr.le.2000) then
      jb1=1990
      jb2=2000
      irr=11
      endif
      call openunit('SO2_INDh',iuc, .false.,.true.)
      do ip=1,10661
      read(iuc,*) ii,jj,iy,xso2ff
      if (jj<j_0 .or. jj>j_1) cycle
      hso2_all(ii,jj,iy)=xso2ff
      end do
      call closeunit(iuc)
      hso2(:,j_0:j_1,1:2)=hso2_all(:,j_0:j_1,irr:irr+1)
c kg/year to kg/s and x2 to convert from S to SO2
      hso2(:,j_0:j_1,1:2)=hso2(:,j_0:j_1,1:2)  
     *  *2.d0/(365.d0*24.d0*3600.d0)
      endif
c interpolate to model year
c    
      d1=real(ihyr-jb1)
      d2=real(jb2-ihyr)
      d3=real(jb2-jb1)
      SO2_src(:,j_0:j_1,1)=(d1*hso2(:,j_0:j_1,2)
     * +d2*hso2(:,j_0:j_1,1))/d3
  
      end subroutine read_hist_SO2    
      
      SUBROUTINE read_hist_NH3(iact)
c historic BC emissions
      USE MODEL_COM, only: jyear, jday,im,jm
      USE CONSTANT, only: sday
      USE FILEMANAGER, only: openunit,closeunit
      USE DOMAIN_DECOMP, only :  GRID, GET 
      USE TRACER_COM, only: aer_int_yr
      USE AEROSOL_SOURCES, only: NH3_src_con,
     * NH3_src_cyc,hnh3_cyc,hnh3_con
      USE GEOM, only: dxyp
      implicit none
      integer iuc,irr,ihyr,i,j,id,jb1,jb2,iact,ii,jj,nn,
     * iy,ip, j_0,j_1,j_0h,j_1h
      real*8, dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,12) :: 
     * hnh3_con_all,hnh3_cyc_all
      real*8 d1,d2,d3,xso2ff
      save jb1,jb2
      
      CALL GET(grid, J_STRT=J_0,       J_STOP=J_1,
     *             J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      
      NH3_src_con(:,:)=0.d0
      NH3_src_cyc(:,:)=0.d0

      if (aer_int_yr.eq.0) then
      ihyr=jyear
      else
      ihyr=aer_int_yr
      endif
c if run just starting or if it is a new year
c   then open new files
      if (iact.eq.0.or.jday.eq.1) then
       hnh3_con_all(:,:,:)=0.0
       hnh3_con(:,:,:)=0.0
       hnh3_cyc_all(:,:,:)=0.0
       hnh3_cyc(:,:,:)=0.0
      if (ihyr.ge.1890.and.ihyr.lt.1900) then
      jb1=1890
      jb2=1900
      irr=1
      else if (ihyr.ge.1900.and.ihyr.lt.1910) then
      jb1=1900
      jb2=1910
      irr=2
      else if (ihyr.ge.1910.and.ihyr.le.1920) then
      jb1=1910
      jb2=1920
      irr=3
      else if (ihyr.ge.1920.and.ihyr.le.1930) then
      jb1=1920
      jb2=1930
      irr=4
      else if (ihyr.ge.1930.and.ihyr.le.1940) then
      jb1=1930
      jb2=1940
      irr=5
      else if (ihyr.ge.1940.and.ihyr.le.1950) then
      jb1=1940
      jb2=1950
      irr=6
      else if (ihyr.ge.1950.and.ihyr.le.1960) then
      jb1=1950
      jb2=1960
      irr=7
      else if (ihyr.ge.1960.and.ihyr.le.1970) then
      jb1=1960
      jb2=1970
      irr=8
      else if (ihyr.ge.1970.and.ihyr.le.1980) then
      jb1=1970
      jb2=1980
      irr=9
      else if (ihyr.ge.1980.and.ihyr.le.2000) then
      jb1=1980
      jb2=1990
      irr=10
      else if (ihyr.ge.1990.and.ihyr.le.2000) then
      jb1=1990
      jb2=2000
      irr=11
      endif
      call openunit('NH3SOURCE_CYC',iuc, .false.,.true.)
      do ip=1,8909
      read(iuc,*) ii,jj,iy,xso2ff
      if (jj<j_0 .or. jj>j_1) cycle
      hnh3_cyc_all(ii,jj,iy)=xso2ff
      end do
      call closeunit(iuc)
      call openunit('NH3SOURCE_CON',iuc, .false.,.true.)
      do ip=1,29532
      read(iuc,*) ii,jj,iy,xso2ff
      if (jj<j_0 .or. jj>j_1) cycle
      hnh3_con_all(ii,jj,iy)=xso2ff
      end do
      call closeunit(iuc)
      hnh3_con(:,j_0:j_1,1:2)=hnh3_con_all(:,j_0:j_1,irr:irr+1)
      hnh3_cyc(:,j_0:j_1,1:2)=hnh3_cyc_all(:,j_0:j_1,irr:irr+1)
c conversion [kg N/gb/a] -> [kg NH3 /gb/s]
      hnh3_con(:,j_0:j_1,1:2)=hnh3_con(:,j_0:j_1,1:2)  
     *  *1.2142/(sday*30.4*12.)
      hnh3_cyc(:,j_0:j_1,1:2)=hnh3_cyc(:,j_0:j_1,1:2)  
     *  *1.2142/(sday*30.4*12.)
      endif
c interpolate to model year
c    
      d1=real(ihyr-jb1)
      d2=real(jb2-ihyr)
      d3=real(jb2-jb1)
      NH3_src_con(:,j_0:j_1)=(d1*hnh3_con(:,j_0:j_1,2)
     * +d2*hnh3_con(:,j_0:j_1,1))/d3
      NH3_src_cyc(:,j_0:j_1)=(d1*hnh3_cyc(:,j_0:j_1,2)
     * +d2*hnh3_cyc(:,j_0:j_1,1))/d3
  
      end subroutine read_hist_NH3      
      

      subroutine read_SO2_source(nt)
!@sum reads in  biomass SO2 source
!@auth Koch
c want kg SO2/m2/s
      USE MODEL_COM, only: im,jm,jmon,jday
      USE DOMAIN_DECOMP, only : GRID, GET
      USE TRACER_COM
      USE FILEMANAGER, only: openunit,closeunit
      USE AEROSOL_SOURCES, only: SO2_biosrc_3D,lmAER
      implicit none
c biomass burning parameters:
Cewg Above 56N, allow biomass burning during the months of May(1.5%),
Cewg June(36.9%), July(44.8%), August(16.2%), and September(0.6%).
Cewg Divide monthly totals by number of time steps in the month.
      real*8, parameter :: B60N(12)=(/0.0, 0.0, 0.0, 0.0, 0.015,
     *  0.369, 0.448, 0.162, 0.006, 0.0, 0.0, 0.0/)
Cewg Between 40N and 56N, allow biomass burning during May(25.8%),
Cewg June(36.6%), July(23.5%), August(12.5%), and September(1.7%).
      real*8, parameter :: B40N(12)=(/0.0, 0.0, 0.0, 0.0, 0.258,
     *  0.366, 0.235, 0.125, 0.017, 0.0, 0.0, 0.0/)
Cewg Between 24N and 40N, allow biomass burning from September to April.
      real*8, parameter :: B24N(12)=(/0.125,0.125,0.125,0.125,
     * 0.0, 0.0,  0.0,0.0, 0.125, 0.125, 0.125, 0.125/)
      REAL*8 TB,TBX,TBXX,TBY,TBYY,TBXY,SDDAY, 
     *  endday,fdry,addtc
      integer nt,
     *  iuc2,ib,jb,iburn,
     *  j_0,j_1
      
      CALL GET(grid, J_STRT=J_0,       J_STOP=J_1)
      
C  Read in emissions from biomass burning
       so2_biosrc_3D(:,:,1:lmAER,1:12)= 0.d0
      call openunit('SO2_BIOMASS',iuc2,.false.,.true.)
      DO 
Cewg  SDDAY -- first day (as day of the year) of the 90 driest days
      READ (iuc2,9051) IB,JB,TB,TBX,TBXX,TBY,TBYY,TBXY,SDDAY
 9051  FORMAT(4X,2I5,6E10.3,F5.0)
      IF (IB.EQ.0) exit
      if (jb<j_0 .or. jb>j_1) cycle
C Flux is tons SO2 / 4x5 box / year. Convert to kg SO2 / 4x5 box /sec
c Must be tons SO2/4x5 box/month. Convert to kg SO2/4x5 box/sec
c Must be tons S/4x5 box/month. Convert to kg SO2/4x5 box/sec
       TB = TB * 907.2d0/30.4d0/24.d0/3600.d0/32.d0*64.d0
        IF ((SDDAY+89.) .LE. 365.) THEN
          ENDDAY = SDDAY + 89.
        ELSE
          ENDDAY = (SDDAY + 89.) - 365.
        ENDIF
Cewg Allow burning above 40N only in May, June, July, Aug, and Sept.
      IBURN = 0
      IF (JMON .GE. 5  .AND.  JMON .LE. 9) IBURN = 1
        IF (JB .GE. 34  .AND.  IBURN .EQ. 0) cycle
        IF (JB .GE. 39  .AND.  IBURN .EQ. 1) THEN
          ADDTC = TB*B60N(JMON)
           GO TO 165
        else IF ((JB.GE.34 .AND. JB.LE.38)  .AND.  IBURN .EQ. 1) THEN
          ADDTC = TB*B40N(JMON)
           GO TO 165
Cewg Allow burning between 32N & 40N from September through April
        else IF (JB .GE. 30  .AND.  JB .LE. 33) THEN
          IF (JMON .LE. 4  .OR.  JMON .GE. 9) THEN
             ADDTC = TB*B24N(JMON)
             GO TO 165
          ELSE
            cycle 
          ENDIF
        ENDIF
Cewg Allow burning south of 32N on the 90 driest days of the year
        FDRY = 0.3333333d0
        IF (ENDDAY .LT. SDDAY) THEN
          IF (JDAY.GT.ENDDAY.AND.JDAY.LT.SDDAY) cycle
        ELSE
          IF (JDAY.LT.SDDAY.OR.JDAY.GT.ENDDAY) cycle
        ENDIF
        ADDTC = TB*FDRY
 165    if (imPI.eq.1) ADDTC=0.5d0*ADDTC
       SO2_biosrc_3D(IB,JB,1,jmon) =  ADDTC
       end do
      call closeunit(iuc2)

      end subroutine read_SO2_source

      subroutine read_DMS_sources(swind,itype,i,j,DMS_flux)
!@sum generates DMS ocean source
!@auth Koch
c Monthly DMS ocean concentration sources are read in and combined
c  with wind and ocean temperature functions to get DMS air surface
c  concentrations
c want kg DMS/m2/s
      USE CONSTANT, only: sday
      USE GEOM, only: dxyp
      USE TRACER_COM, only: tr_mm,n_DMS,imAER
      USE MODEL_COM, only: jmon,jday
      USE AEROSOL_SOURCES, only: DMSinput,DMS_AER
      implicit none
      integer jread
      REAL*8 akw,erate
      real*8, INTENT(OUT) :: DMS_flux
      real*8, INTENT(IN) :: swind
      integer, INTENT(IN) :: itype,i,j

      DMS_flux=0.d0
        erate=0.d0
        if (imAER.ne.1) then
        if (itype.eq.1) then
c Nightingale et al
        akw = 0.23d0*swind*swind + 0.1d0 * swind
        akw = akw * 0.24d0
        erate=akw*DMSinput(i,j,jmon)*1.d-9*62.d0 !*tr_mm(nt)
     *       /sday
        endif
        else !AEROCOM run, prescribed flux

c if after Feb 28 skip the leapyear day
         jread=jday
         if (jday.gt.59) jread=jday+1
c         if (j.eq.1.or.j.eq.46) DMS_AER(i,j,jread)
c     *      =DMS_AER(i,j,jread)*72.d0
         erate=DMS_AER(i,j,jread)/sday/dxyp(j)*tr_mm(n_DMS)/32.d0
        endif
        DMS_flux=erate          ! units are kg/m2/s
c
      return
      end subroutine read_DMS_sources

      subroutine read_seasalt_sources(swind,itype,ibin,i,j,ss)
!@sum determines wind-speed dependent oceanic seasalt source
!@auth Koch
c want kg seasalt/m2/s, for now in 2 size bins
      USE TRACER_COM, only: imAER
      USE CONSTANT, only: sday
      USE GEOM, only: dxyp
      USE MODEL_COM, only: jday
      USE AEROSOL_SOURCES, only: SS1_AER,SS2_AER
      implicit none
      REAL*8 erate
      integer jread
      integer, INTENT(IN)::itype,ibin,i,j
      REAL*8, INTENT(IN)::swind
      REAL*8, INTENT(OUT)::ss
c
      ss=0.
        erate=0.d0
       if (imAER.ne.1) then
       if (itype.eq.1) then
c Monahan 1971, bubble source, important for small (<10um) particles
        erate= 1.373d0 * swind**(3.41d0)
        if (ibin.eq.1) then
          ss=erate*2.11d-14     ! submicron (0.1 < r_d < 1.)
        else
          ss=erate*7.78d-14     ! supermicron (1. < r_d < 4.)
        endif
c     units are kg salt/m2/s
       endif
       else
c if after Feb 28 skip the leapyear day
         jread=jday
         if (jday.gt.59) jread=jday+1
        if (ibin.eq.1) then
         ss=SS1_AER(i,j,jread)/sday/dxyp(j)
        else 
         ss=SS2_AER(i,j,jread)/sday/dxyp(j)
        endif
       endif
      return
      end subroutine read_seasalt_sources


      subroutine aerosol_gas_chem
!@sum aerosol gas phase chemistry
!@auth Dorothy Koch
      USE TRACER_COM
      USE TRDIAG_COM, only : tajls=>tajls_loc
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) 
     *     ,jls_OHconk,jls_HO2con,jls_NO3,jls_phot
#endif
#ifdef TRACERS_SPECIAL_Shindell
     &     ,jls_OHcon
#endif
      USE DOMAIN_DECOMP, only: AM_I_ROOT
      USE DOMAIN_DECOMP, only : GRID, GET, UNPACK_DATA, write_parallel
      USE MODEL_COM, only: im,jm,jmon,ls1,lm,dtsrc,t,q,jday,
     * coupled_chem
      USE DYNAMICS, only: pmid,am,pk,LTROPO,byam
      USE GEOM, only: dxyp,imaxj,BYDXYP
      USE FLUXES, only: tr3Dsource
      USE FILEMANAGER, only: openunit,closeunit
      USE AEROSOL_SOURCES, only: ohr,dho2r,perjr,tno3r,oh,
     & dho2,perj,tno3,ohsr,o3_offline
       USE CONSTANT, only : mair
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM, only: which_trop
#endif
c Aerosol chemistry
      implicit none
      logical :: ifirst=.true.
      real*8 ppres,te,tt,mm,dmm,ohmc,r1,d1,r2,d2,ttno3,r3,d3,
     * ddno3,dddms,ddno3a,fmom,dtt
      real*8 rk4,ek4,r4,d4
      real*8 r6,d6,ek9,ek9t,ch2o,eh2o,dho2mc,dho2kg,eeee,xk9,
     * r5,d5,dmssink
#ifdef TRACERS_HETCHEM
     *       ,d41,d42,d43,o3mc,rsulfo3
#endif
      real*8 bciage,ociage
c     real*8, dimension(im,jm,lm,12) :: ohr_glob,dho2r_glob,
c    * perjr_glob,tno3r_glob,ohsr_glob
      real*8, dimension(im,jm,lm) :: 
     * ohr_globm,dho2r_globm,
     * perjr_globm,tno3r_globm,ohsr_glob
      real*4, dimension(im,jm) :: 
     *  ohsr_real
      integer i,j,l,n,iuc,iun,itau,ixx1,ixx2,ichemi,itt,
     * ittime,isp,iix,jjx,llx,ii,jj,ll,iuc2,it,nm,najl,j_0,j_1,
     * j_0s,j_1s,mmm,J_0H,J_1H
#ifdef TRACERS_SPECIAL_Shindell
!@var maxl chosen tropopause 0=LTROPO(I,J), 1=LS1-1
#endif
      integer maxl
      save ifirst
      
      CALL GET(grid, J_STRT=J_0,J_STOP=J_1,
     *    J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,J_STRT_SKP=J_0S,
     * J_STOP_SKP=J_1S)

C**** initialise source arrays
ccOMP PARALLEL DO PRIVATE (L)
        tr3Dsource(:,j_0:j_1,:,1,n_DMS)=0. ! DMS chem sink
#ifndef TRACERS_AMP
        tr3Dsource(:,j_0:j_1,:,1,n_MSA)=0. ! MSA chem sink
        tr3Dsource(:,j_0:j_1,:,1,n_SO4)=0. ! SO4 chem source
#endif
        tr3Dsource(:,j_0:j_1,:,4,n_SO2)=0. ! SO2 chem source
        tr3Dsource(:,j_0:j_1,:,5,n_SO2)=0. ! SO2 chem sink
        tr3Dsource(:,j_0:j_1,:,1,n_H2O2_s)=0. ! H2O2 chem source
        tr3Dsource(:,j_0:j_1,:,2,n_H2O2_s)=0. ! H2O2 chem sink
#ifdef TRACERS_AMP
        tr3Dsource(:,j_0:j_1,:,2,n_H2SO4)=0. ! H2O2 chem sink
#endif
#ifdef TRACERS_HETCHEM
        tr3Dsource(:,j_0:j_1,:,1,n_SO4_d1) =0. ! SO4 on dust
        tr3Dsource(:,j_0:j_1,:,1,n_SO4_d2) =0. ! SO4 on dust
        tr3Dsource(:,j_0:j_1,:,1,n_SO4_d3) =0. ! SO4 on dust
#endif
        if (n_BCII.gt.0) then
          tr3Dsource(:,j_0:j_1,:,1,n_BCII)=0. ! BCII sink
          tr3Dsource(:,j_0:j_1,:,1,n_BCIA)=0. ! BCIA source
        end if
        if (n_OCII.gt.0) then
          tr3Dsource(:,j_0:j_1,:,1,n_OCII)=0. ! OCII sink
          tr3Dsource(:,j_0:j_1,:,1,n_OCIA)=0. ! OCIA source
        end if
        if (n_OCI1.gt.0) then
          tr3Dsource(:,j_0:j_1,:,2,n_OCI1)=0. ! OCI1 sink
          tr3Dsource(:,j_0:j_1,:,1,n_OCA1)=0. ! OCA1 source
        end if
        if (n_OCI2.gt.0) then
          tr3Dsource(:,j_0:j_1,:,2,n_OCI2)=0. ! OCI2 sink
          tr3Dsource(:,j_0:j_1,:,1,n_OCA2)=0. ! OCA2 source
        end if
        if (n_OCI3.gt.0) then
          tr3Dsource(:,j_0:j_1,:,1,n_OCI3)=0. ! OCI3 sink
          tr3Dsource(:,j_0:j_1,:,1,n_OCA3)=0. ! OCA3 source
        end if
ccOMP END PARALLEL DO

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
C Coupled mode: use on-line radical concentrations
      if (coupled_chem.eq.1) then
ccOMP PARALLEL DO PRIVATE (L)
          oh(:,j_0:j_1,:)=oh_live(:,j_0:j_1,:)
          tno3(:,j_0:j_1,:)=no3_live(:,j_0:j_1,:)
c Set h2o2_s =0 and use on-line h2o2 from chemistry
          trm(:,j_0:j_1,:,n_h2o2_s)=0.0
ccOMP END PARALLEL DO
      endif

      if (coupled_chem.eq.0) then
c Use this for chem inputs from B4360C0M23, from Drew
c      if (ifirst) then
         if(AM_I_ROOT( ))then
        call openunit('AER_CHEM',iuc,.true.,.true.)
        do ii=1,jmon
          read(iuc) ichemi
          read(iuc) ohr_globm
          read(iuc) dho2r_globm
          read(iuc) perjr_globm
          read(iuc) tno3r_globm
c         ohr_glob(:,:,:,ii)=ohr_globm(:,:,:)
c         dho2r_glob(:,:,:,ii)=dho2r_globm(:,:,:)
c         perjr_glob(:,:,:,ii)=perjr_globm(:,:,:)
c         tno3r_glob(:,:,:,ii)=tno3r_globm(:,:,:)
        end do
cDMK I could move these loops outside AM_I_ROOT
        call closeunit(iuc)
        endif
        
        call UNPACK_DATA( grid, ohr_globm, ohr )
        call UNPACK_DATA( grid, dho2r_globm, dho2r )
        call UNPACK_DATA( grid, perjr_globm, perjr )
        call UNPACK_DATA( grid, tno3r_globm, tno3r )

         if(AM_I_ROOT( ))then
        call openunit('AER_OH_STRAT',iuc2,.true.,.true.)
        do ii=1,jmon  !12
        do ll=1,lm
         read(iuc2) ohsr_real
          ohsr_glob(:,:,ll)=ohsr_real(:,:)*1.D5
        end do
        end do
        call closeunit(iuc2)
        endif
        call UNPACK_DATA( grid, ohsr_glob, ohsr )
c       ifirst=.false.
c       endif
c I have to read in every timestep unless I can find a better way
c
c skip poles because there was a bug in the input file over the pole
        do j=j_0s,j_1s   
        do i=1,im
        maxl=ltropo(i,j)
        do l=maxl,lm
        ohr(i,j,l)=ohsr(i,j,l)
        end do
        end do
        end do
c impose diurnal variability
        CALL SCALERAD
      endif

C Calculation of gas phase reaction rates
C In coupled chemistry case, routine called from masterchem
c      if (coupled_chem.eq.0) then
#ifndef  TRACERS_SPECIAL_Shindell  
      CALL GET_SULF_GAS_RATES
#endif
c      endif

#ifdef TRACERS_HETCHEM
c calculation of heterogeneous reaction rates: SO2 on dust 
      CALL SULFDUST
c calculation of heterogeneous reaction rates: SO2 on seasalt
c      CALL SULFSEAS 
c     if (COUPLED_CHEM.ne.1) then
c     CALL GET_O3_OFFLINE
c     endif
#endif
#endif
      dtt=dtsrc
C**** THIS LOOP SHOULD BE PARALLELISED
      do 20 l=1,lm
      do 21 j=j_0,j_1
      do 22 i=1,imaxj(j)
C Initialise       
        maxl = ltropo(i,j)
#ifdef TRACERS_SPECIAL_Shindell
        if(which_trop.eq.1)maxl=ls1-1
#endif
      if(l.le.maxl) then
      ppres=pmid(l,i,j)*9.869d-4 !in atm
      te=pk(l,i,j)*t(i,j,l)
      mm=am(l,i,j)*dxyp(j)
      tt = 1.d0/te

c DMM is number density of air in molecules/cm3
      dmm=ppres/(.082d0*te)*6.02d20
      ohmc = oh(i,j,l)          !oh is alread in units of molecules/cm3

      do 23 n=1,ntm

        select case (trname(n))
c    Aging of industrial carbonaceous aerosols 
        case ('BCII')
        bciage=4.3D-6*trm(i,j,l,n) !used this first        
c       bciage=1.0D-6*trm(i,j,l,n)  !2nd
c       bciage=1.0D-7*trm(i,j,l,n)
        tr3Dsource(i,j,l,1,n)=-bciage        
        tr3Dsource(i,j,l,1,n_BCIA)=bciage        

        case ('OCII')
        ociage=7.3D-6*trm(i,j,l,n)       !used this first 
c       ociage=3.6D-6*trm(i,j,l,n)     !2nd
c       ociage=3.D-7*trm(i,j,l,n)
        tr3Dsource(i,j,l,1,n)=-ociage        
        tr3Dsource(i,j,l,1,n_OCIA)=ociage        

        case ('OCI1')
        ociage=4.3D-6*trm(i,j,l,n)
        tr3Dsource(i,j,l,2,n)=-ociage
        tr3Dsource(i,j,l,1,n_OCA1)=ociage
        case ('OCI2')
        ociage=4.3D-6*trm(i,j,l,n)
        tr3Dsource(i,j,l,2,n)=-ociage
        tr3Dsource(i,j,l,1,n_OCA2)=ociage
        case ('OCI3')
        ociage=4.3D-6*trm(i,j,l,n)
        tr3Dsource(i,j,l,2,n)=-ociage
        tr3Dsource(i,j,l,1,n_OCA3)=ociage
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
        case ('DMS')
C***1.DMS + OH -> 0.75SO2 + 0.25MSA
C***2.DMS + OH -> SO2
C***3.DMS + NO3 -> HNO3 + SO2

          r1=rsulf1(i,j,l)*ohmc 
          d1 = exp(-r1*dtsrc)
          r2=rsulf2(i,j,l)*ohmc
          d2 = exp(-r2*dtsrc)

c     NO3 is in mixing ratio: convert to molecules/cm3
c - not necessary for Shindell source
          if (l.gt.8) then
            ttno3=0.d0
            go to 87
          endif
          ttno3 = tno3(i,j,l)   !*6.02d20*ppres/(.082056d0*te)
 87       r3=rsulf3(i,j,l)*ttno3
          d3= exp(-r3*dtsrc)
          ddno3=r3*trm(i,j,l,n)/tr_mm(n)*1000.d0*dtsrc
          dddms=trm(i,j,l,n)/tr_mm(n)*1000.d0
          if (ddno3.gt.dddms) ddno3=dddms

          ddno3=ddno3*0.9
C DMS losses: eqns 1, 2 ,3

          tr3Dsource(i,j,l,1,n) = trm(i,j,l,n)*(d1*d2-1.)/dtsrc

          dmssink=ddno3*tr_mm(n)/1000.d0

          if (dmssink.gt.trm(i,j,l,n)+tr3Dsource(i,j,l,1,n)*dtsrc)
     *         dmssink=trm(i,j,l,n)+tr3Dsource(i,j,l,1,n)*dtsrc
          tr3Dsource(i,j,l,1,n) = tr3Dsource(i,j,l,1,n) - dmssink/dtsrc
          
        case ('MSA')
C MSA gain: eqn 1

          tr3Dsource(i,j,l,1,n) = 0.25d0*Tr_mm(n)/Tr_mm(n_dms)*trm(i,j
     *         ,l,n_dms)*(1.d0 -D1)*SQRT(D2)/dtsrc
          
        case ('SO2')
c SO2 production from DMS
          tr3Dsource(i,j,l,4,n) = (0.75*tr_mm(n)/tr_mm(n_dms)*trm(i,j,l
     *         ,n_dms)*(1.d0 - d1)*sqrt(d2)+ tr_mm(n)/tr_mm(n_dms)*trm(i
     *         ,j,l,n_dms)*(1.d0 - d2)*sqrt(d1)+dmssink*tr_mm(n)
     *         /tr_mm(n_dms))/dtsrc
          
          
          najl = jls_NO3
          tajls(j,l,najl) = tajls(j,l,najl)+ttno3
#endif
        end select
        
 23   CONTINUE
        endif
 22   CONTINUE
 21   CONTINUE
 20   CONTINUE

      do 30 l=1,lm
      do 31 j=j_0,j_1
      do 32 i=1,imaxj(j)

      maxl = ltropo(i,j)
#ifdef TRACERS_SPECIAL_Shindell
      if(which_trop.eq.1)maxl=ls1-1
#endif
      if(l.le.maxl) then

      ppres=pmid(l,i,j)*9.869d-4 !in atm
      te=pk(l,i,j)*t(i,j,l)
      mm=am(l,i,j)*dxyp(j)
      tt = 1.d0/te
      dmm=ppres/(.082d0*te)*6.02d20
      ohmc = oh(i,j,l)          !oh is alread in units of molecules/cm3
#ifdef TRACERS_HETCHEM
      if (COUPLED_CHEM.ne.1) then
      o3mc = o3_offline(i,j,l)*dmm*(28.0D0/48.0D0)*BYDXYP(J)*BYAM(L,I,J)
      else
      o3mc = trm(i,j,l,n_Ox)*dmm*(28.0D0/48.0D0)*BYDXYP(J)*BYAM(L,I,J)
      endif
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)
      do 33 n=1,ntm
        select case (trname(n))
        case ('SO2')
c oxidation of SO2 to make SO4: SO2 + OH -> H2SO4

          r4=rsulf4(i,j,l)*ohmc
          d4 = exp(-r4*dtsrc)

          IF (d4.GE.1.) d4=0.99999d0
#ifdef TRACERS_HETCHEM
      rsulfo3 = 4.39d11*exp(-4131/te)+( 2.56d3*exp(-966/te)) * 10.d5 !assuming pH=5
      rsulfo3 = exp(-rsulfo3*o3mc *dtsrc) !O3 oxidation Maahs '83
       d41 = exp(-rxts1(i,j,l)*dtsrc)     
       d42 = exp(-rxts2(i,j,l)*dtsrc)     
       d43 = exp(-rxts3(i,j,l)*dtsrc)     
       tr3Dsource(i,j,l,5,n) = (-trm(i,j,l,n)*(1.d0-d41)/dtsrc)
     .                       + ( -trm(i,j,l,n)*(1.d0-d4)/dtsrc)
     .                       + ( -trm(i,j,l,n)*(1.d0-d42)/dtsrc)
     .                       + ( -trm(i,j,l,n)*(1.d0-d43)/dtsrc)
#else
       tr3Dsource(i,j,l,5,n) = -trm(i,j,l,n)*(1.d0-d4)/dtsrc 
#ifdef TRACERS_AMP
       tr3Dsource(i,j,l,2,n_H2SO4)=trm(i,j,l,n)*(1.d0-d4)/dtsrc 
#endif        
#endif        
c diagnostics to save oxidant fields
c No need to accumulate Shindell version here because it
c   is done elsewhere
c#ifdef TRACERS_SPECIAL_Shindell
c         najl = jls_OHcon
c#else
          najl = jls_OHconk
c#endif
          tajls(j,l,najl) = tajls(j,l,najl)+oh(i,j,l)
          najl = jls_HO2con
          tajls(j,l,najl) = tajls(j,l,najl)+dho2(i,j,l)

#ifdef TRACERS_HETCHEM
       case ('SO4_d1')
c sulfate production from SO2 on mineral dust aerosol due to O3 oxidation

       tr3Dsource(i,j,l,1,n)=tr3Dsource(i,j,l,1,n)+tr_mm(n)/tr_mm(n_so2)
     *         *(1.d0-d41)*trm(i,j,l,n_so2)            !  SO2
     *         * (1.d0-rsulfo3)                              !+ O3
     *           /dtsrc
       case ('SO4_d2')
c sulfate production from SO2 on mineral dust aerosol

       tr3Dsource(i,j,l,1,n) = tr3Dsource(i,j,l,1,n) +  tr_mm(n)/
     *             tr_mm(n_so2)*(1.d0-d42)*trm(i,j,l,n_so2)
     *         * (1.d0-rsulfo3)                              !+ O3
     *           /dtsrc
       case ('SO4_d3')
c sulfate production from SO2 on mineral dust aerosol

       tr3Dsource(i,j,l,1,n) = tr3Dsource(i,j,l,1,n) +  tr_mm(n)/
     *             tr_mm(n_so2)*(1.d0-d43)*trm(i,j,l,n_so2)
     *         * (1.d0-rsulfo3)                              !+ O3
     *           /dtsrc

#endif
        case('SO4')
C SO4 production
          tr3Dsource(i,j,l,1,n) = tr3Dsource(i,j,l,1,n)+tr_mm(n)
     *         /tr_mm(n_so2)*trm(i,j,l,n_so2)*(1.d0 -d4)/dtsrc
        case('H2O2_s')

          if (coupled_chem.eq.1) go to 140

c hydrogen peroxide formation and destruction:
C***5.H2O2 +hv -> 2OH
C***6.H2O2 + OH -> H2O + HO2
C***7.H2O2 + SO2 -> H2O + SO3 (in-cloud, in CB)
C***9.HO2 + HO2 ->H2O2 + O2
C     HO2 + HO2 + M ->
C     HO2 + HO2 + H2O ->
C     HO2 + HO2 + H2O + M ->

          r6 = 2.9d-12 * exp(-160.d0*tt)*ohmc
          d6 = exp(-r6*dtsrc)
          ek9 = 2.2d-13*exp(600.d0*tt)
          ek9t = 1.9d-20*dmm*0.78d0*exp(980.d0*tt)*1.d-13
          ch2o = q(i,j,l)*6.02d20*28.97d0/18.d0*ppres/(.082d0*te)
          eh2o = 1.+1.4d-21*exp(2200.d0*tt)*ch2o
          dho2mc = dho2(i,j,l)  !/mm*1.292/.033*6.02e17
          dho2kg = dho2(i,j,l)*mm*te*.082056d0/(ppres*28.97d0*6.02d20)
          eeee = eh2o*(ek9+ek9t)*dtt*dho2mc
          xk9 = dho2kg*eeee
c H2O2 production: eqn 9
         
          tr3Dsource(i,j,l,1,n) = tr_mm(n)*xk9/dtsrc
c H2O2 losses:5 and 6
          r5 = perj(i,j,l)
          d5 = exp(-r5*dtsrc)

          tr3Dsource(i,j,l,2,n)=(trm(i,j,l,n))*(d5*d6-1.d0)
     *         /dtsrc
          
          najl = jls_phot
          tajls(j,l,najl) = tajls(j,l,najl)+perj(i,j,l)
        end select

 140    CONTINUE

 33   CONTINUE
#endif
      endif
 32   CONTINUE
 31   CONTINUE
 30   CONTINUE

      RETURN
      END subroutine aerosol_gas_chem

      SUBROUTINE SCALERAD
      use MODEL_COM, only: im,jm,lm,jday,jhour,jmon,itime,nday,dtsrc
      use AEROSOL_SOURCES, only: ohr,dho2r,perjr,tno3r,oh,dho2,perj,tno3
      USE DOMAIN_DECOMP, only:GRID,GET      
      use CONSTANT, only: radian,teeny
c     use RAD_COM, only: cosz1
      implicit none
      real*8 ang1,xnair,vlon,vlat,ctime,timec,p1,p2,p3,fact,rad,
     *  rad1,rad2,rad3,stfac
      real*8, DIMENSION(IM,JM) :: suncos,tczen
      integer i,j,hrstrt,jdstrt,ihr,l,j_0,j_1,j_0h,j_1h,
     * j_0s,j_1s 
      integer, DIMENSION(IM,JM) :: nradn
      LOGICAL HAVE_SOUTH_POLE, HAVE_NORTH_POLE      

      ang1 = 90.d0/91.3125d0
      xnair = 6.022d23 * 1.d3 / 28.9644d0
      jdstrt=0
      hrstrt=0
      CALL GET(grid, J_STRT=J_0,       J_STOP=J_1,
     *               J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     *               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     *               HAVE_NORTH_POLE=HAVE_NORTH_POLE)
C*
      DO 100 j = j_0,j_1
      DO 100 i = 1,im
C*** Calculate cos theta (RAD) at the beginning of the time step (RAD)
      vlon = 180.d0 - (i-1)*5.d0
      vlat=-90.d0+(j-1)*4.d0
      if (HAVE_SOUTH_POLE.and.j.eq.j_0s-1) vlat=-88.d0
      if (HAVE_NORTH_POLE.and.j.eq.j_1s+1) vlat=88.d0      
c      if (j.eq.1) vlat=-88.d0
c      if (j.eq.46) vlat=88.d0
      TIMEC = (mod(itime,365*nday) + 0.5)*DTsrc
      p1 = 15.d0*(timec/3600.d0 + hrstrt - vlon/15.d0 - 12.d0)
      fact = (jdstrt + timec/86400.d0 - 81.1875d0)*ang1
      p2 = 23.5d0*sin(fact*radian)
      p3 = vlat
C*
      rad = (SIN(p3*radian)*SIN(p2*radian)) +
     1            (COS(p1*radian)*COS(p2*radian)*COS(p3*radian))
      if (rad.lt.0.d0) rad = 0.d0
      suncos(I,J) = rad
c     if (i.eq.30.and.j.eq.40) write(6,*) 'ohrad',timec,jday,jhour,
c    * rad,p1,p2,p3
 100  CONTINUE
c Scale OH and PERJ depending on time of day
c   [OH] is approximately proportional to cos(zenith angle),
c   therefore OH is calculated as
c **
c **    OH4HR(i,j,l) = OH5DY(i,j,l) * 6 * SUNCOS(i,j) / TCZEN(i,j)
c **
c ** where TCZEN(i,j) is the total cos(zenith angle) calculated at the
c ** beginning of the day.
c etc for NO3, PERJ. The 6 is because TCZEN is average of 6 times
c    NO3 is different since we want to have it nonzero only during dark
C----------------------------------------------------------------------
        do j = j_0, j_1
           do i = 1, im
             tczen(i,j) = 0.
             nradn(i,j)=0
             vlon = 180.d0 - (i-1)*5.d0
             vlat=-90.d0+(j-1)*4.d0
      if (HAVE_SOUTH_POLE.and.j.eq.j_0s-1) vlat=-88.d0
      if (HAVE_NORTH_POLE.and.j.eq.j_1s+1) vlat=88.d0      
c             if (j.eq.1) vlat=-88.d0
c             if (j.eq.46) vlat=88.d0
              do ihr = 0,20,4
                 ctime = (jday*24.d0) + ihr
                 p1 = 15.d0*(ctime - vlon/15.d0 - 12.d0)
                 fact = (ctime/24.d0 - 81.1875d0)*ang1
                 p2 = 23.5d0 * sin(fact*radian)
                 p3 = vlat
                 rad1 = sin(p3*radian) * sin(p2*radian) +
     m             cos(p3*radian) * cos(p2*radian) * cos(p1*radian)
                 if (rad1.lt.0.d0) rad1 = 0.d0
                 if (rad1.eq.0.0) nradn(i,j)=nradn(i,j)+1
c
                 ctime = (jday*24.d0) + ihr+2
                 p1 = 15.d0*(ctime - vlon/15.d0 - 12.d0)
                 fact = (ctime/24.d0 - 81.1875d0)*ang1
                 p2 = 23.5d0 * sin(fact*radian)
                 p3 = vlat
                 rad2 = sin(p3*radian) * sin(p2*radian) +
     m             cos(p3*radian) * cos(p2*radian) * cos(p1*radian)
                 if (rad2.lt.0.d0) rad2 = 0.d0
                 if (rad2.eq.0.0) nradn(i,j)=nradn(i,j)+2
c
                 ctime = (jday*24) + ihr+4
                 p1 = 15.d0*(ctime - vlon/15.d0 - 12.d0)
                 fact = (ctime/24.d0 - 81.1875d0)*ang1
                 p2 = 23.5d0 * sin(fact*radian)
                 p3 = vlat
                 rad3 = sin(p3*radian) * sin(p2*radian) +
     m              cos(p3*radian) * cos(p2*radian) * cos(p1*radian)
                 if (rad3.lt.0.d0) rad3 = 0.d0
                 if (rad3.eq.0.0) nradn(i,j)=nradn(i,j)+1
c
                 rad = (rad1 + 2.d0*rad2 + rad3)/4.d0
                 tczen(i,j) = tczen(i,j) + rad
                 if(tczen(i,j).eq.0.) tczen(i,j) = teeny
              end do
           end do
        end do
c**************************
          do j = j_0,j_1
           do l = 1,lm
            do i = 1,im
c Get NO3 only if dark, weighted by number of dark hours
c        if (I.EQ.1.AND.L.EQ.1) write(6,*)'NO3R',TAU,J,NRADN(I,J)
            if (suncos(i,j).eq.0.and.nradn(i,j).gt.0) then
            tno3(i,j,l)=tno3r(i,j,l)*24.d0/real(nradn(i,j))  !DMK jmon
            else
            tno3(i,j,l)=0.d0
            endif
  88          if (tczen(i,j).eq.teeny) then
                 oh(i,j,l) = 0.d0
                 perj(i,j,l)=0.d0
                 dho2(i,j,l)=0.d0
              else
                stfac=suncos(i,j)/tczen(i,j)
                if (stfac.gt.1) then
                 stfac=1.
                endif
                 oh(i,j,l)=ohr(i,j,l)*6.d0*stfac    !jmon
                 perj(i,j,l)=perjr(i,j,l)*6.d0*stfac !dmk jmon
                 dho2(i,j,l)=dho2r(i,j,l)*6.d0*stfac !dmk jmon
              end if
           end do
         end do
      end do
      RETURN
      END subroutine SCALERAD


      SUBROUTINE GET_SULFATE(L,temp,fcloud,
     *  wa_vol,wmxtr,sulfin,sulfinc,sulfout,tr_left,
     *  tm,tmcl,airm,LHX,dt_sulf,fcld0)
!@sum  GET_SULFATE calculates formation of sulfate from SO2 and H2O2
!@+    within or below convective or large-scale clouds. Gas
!@+    condensation uses Henry's Law if not freezing.
!@auth Dorothy Koch
!@ver  1.0 (based on CLOUDCHCC and CLOUDCHEM subroutines)
c
C**** GLOBAL parameters and variables:
      USE CONSTANT, only: BYGASC, MAIR,teeny,mb2kg,gasc,LHE
      USE TRACER_COM, only: tr_RKD,tr_DHD,n_H2O2_s,n_SO2
     *     ,trname,ntm,tr_mm,lm,n_SO4,n_H2O2,mass2vol
      USE CLOUDS, only: PL,NTIX,NTX,DXYPJ
      USE MODEL_COM, only: dtsrc,coupled_chem
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
c
!@param BY298K unknown meaning for now (assumed= 1./298K)
!@var Ppas pressure at current altitude (in Pascal=kg/s2/m)
!@var TFAC exponential coeffiecient of tracer condensation temperature
!@+   dependence (mole/joule)
!@var FCLOUD fraction of cloud available for tracer condensation
!@var SSFAC dummy variable (assumed units= kg water?)
!@var L index for altitude loop
!@var N index for tracer number loop
!@var RKD dummy variable (= tr_RKD*EXP[ ])
!@var WA_VOL is the cloud water volume in L
!@var WMX_INC is the change in cloud water ratio
!@var CLWC is the cloud liquid water content: L water/L air
!@var sulfin is the amount of SO2 and H2O2 used to make sulfate
!@var sulfout is the amount of sulfate generated
!@var tr_left is the amount of SO2 and H2O2 left after sulfate is made
!@+    and is now available to condense
!@var PPH is the partial pressure of the gas in M/kg
!@var amass is airmass in kg airm(l)*mb2kg
!@var trd is dissolved portion in moles/L
!@var sulfinc is the change in dissolved SO2 and H2O2 as we use those
!@+  to form sulfate (in addition to what is dissolved)
      REAL*8, PARAMETER :: BY298K=3.3557D-3
      REAL*8 Ppas, tfac, ssfac, RKD, Henry_const(ntm)
      real*8 clwc,rk1f,rkdm(ntm),amass,trd(ntm),trdr(ntm),
     * dso4g,dso4d,pph(ntm),trmol(ntm),trdmol(ntm),dso4gt,dso4dt
      integer n,ih,is,is4,ihx,isx,is4x
      real*8, parameter :: rk1=1.3d-2 !M
      real*8, parameter :: dh1=-1.6736d4 !J/mol
      real*8, parameter :: rk=6.357d14    !1/(M*M*s)
      real*8, parameter :: ea=3.95d4 !J/mol
      REAL*8,  INTENT(IN) :: fcloud,temp,wa_vol,wmxtr,LHX
      real*8, dimension(lm,ntm) :: tm
      real*8, dimension(ntm) :: tmcl,sulfin,sulfout,tr_left
     *  ,sulfinc
      real*8, dimension(lm) :: airm
c     REAL*8,  INTENT(OUT)::
!@var dt_sulf accumulated diagnostic of sulfate chemistry changes
      real*8, dimension(ntm), intent(inout) :: dt_sulf
      real*8 finc,fcld0
      INTEGER, INTENT(IN) :: L
      do n=1,ntx
        sulfin(N)=0.
        sulfinc(N)=0.
        sulfout(N)=0.
        tr_left(N)=1.
      end do
c is this needed?
c#if (defined TRACERS_COSMO) || (defined TRACERS_
c      do n=1,ntx
c       select case (trname(ntix(n)))
c       case ('Pb210   ','Be7     ','Be10    ','Rn222')
c       go to 333
c       end select
c      end do
C#endif
c
c
C**** CALCULATE the fraction of tracer mass that becomes condensate:
c
      if (LHX.NE.LHE.or.fcloud.lt.teeny) go to 333
      finc=(fcloud-fcld0)/fcloud
      if (finc.lt.0d0) finc=0.d0
c First allow for formation of sulfate from SO2 and H2O2. Then remaining
c  gases may be allowed to dissolve (amount given by tr_left)
C H2O2 + SO2 -> H2O + SO3 -> H2SO4
      amass=airm(l)*mb2kg*DXYPJ
      Ppas = PL(L)*1.D2
      tfac = (1./temp - by298k)*bygasc  !mol/J
c  cloud liquid water content
      clwc=wmxtr*mair*ppas/temp*bygasc/1.D6/fcloud
      rk1f=rk1*exp(-dh1*tfac)

      do n=1,ntx
       select case (trname(ntix(n)))
       case('SO2')
       is=ntix(n)
       isx=n
c modified Henry's Law coefficient assuming pH of 4.5
      rkdm(is)=tr_rkd(is)*(1.+ rk1f/3.2d-5)
c mole of tracer, used to limit so4 production
      trmol(is)=1000.*tm(l,isx)/tr_mm(is)*fcloud
c partial pressure of gas x henry's law coefficient
      pph(is)=mass2vol(is)*1.d-3*ppas/amass*
     *   tr_rkd(is)*exp(-tr_dhd(is)*tfac)
c the following is from Phil:
c      reduction in partial pressure as species dissolves
      henry_const(is)=rkdm(is)*exp(-tr_dhd(is)*tfac)
      pph(is)=pph(is)/(1+(henry_const(is)*clwc*gasc*temp))
c again all except tmcl(n)
      trdr(is)=mass2vol(is)*ppas/amass*bygasc
     *    /temp*1.D-3  !M/kg
c dissolved moles
      trdmol(is)=trdr(is)*1000./tr_mm(is)

       case('H2O2','H2O2_s')

         if (trname(ntix(n)).eq."H2O2" .and. coupled_chem.eq.0) goto 400
         if (trname(ntix(n)).eq."H2O2_s" .and. coupled_chem.eq.1) goto
     *        400

       ih=ntix(n)
       ihx=n
c modified Henry's Law coefficient assuming pH of 4.5
      rkdm(ih)=tr_rkd(ih)
c mole of tracer, used to limit so4 production
      trmol(ih)=1000.*tm(l,ihx)/tr_mm(ih)*fcloud
c partial pressure of gas x henry's law coefficient
      pph(ih)=mass2vol(ih)*1.D-3*ppas/amass*
     *   tr_rkd(ih)*exp(-tr_dhd(ih)*tfac)
c the following is from Phil:
c      reduction in partial pressure as species dissolves
      henry_const(ih)=rkdm(ih)*exp(-tr_dhd(ih)*tfac)
      pph(ih)=pph(ih)/(1+(henry_const(ih)*clwc*gasc*temp))
c all except tmcl(n)
      trdr(ih)=mass2vol(ih)*ppas/amass*bygasc/temp*1.D-3  !M/kg
c dissolved moles
      trdmol(ih)=trdr(ih)*1000./tr_mm(ih)

 400   CONTINUE

      end select
      end do
      if (tm(l,ihx).lt.teeny.or.tm(l,isx).lt.teeny) then
      dso4g=0.
      go to 21
      endif
c this part from gas phase:moles/kg/kg
      dso4g=rk*exp(-ea/(gasc*temp))*rk1f
     *    *pph(ih)*pph(is)*dtsrc*wa_vol
c dmk  incremental water volume
      dso4g=dso4g*finc
c should probably be (finc+tr_lef) but then tr_lef has to be saved   
c check to make sure no overreaction: moles of production:
      dso4gt=dso4g*tm(l,ihx)*tm(l,isx)
c can't be more than moles going in:
      if (dso4gt.gt.trmol(is)) then
        dso4g=trmol(is)/(tm(l,ihx)*tm(l,isx))
      endif
      dso4gt=dso4g*tm(l,ihx)*tm(l,isx)
      if (dso4gt.gt.trmol(ih)) then
        dso4g=trmol(ih)/(tm(l,ihx)*tm(l,isx))
      endif
c this part from dissolved gases
 21    dso4d=rk*exp(-ea/(gasc*temp))*rk1f
     *    *trdr(ih)*trdr(is)*dtsrc*wa_vol

      if (trdr(ih).lt.teeny.or.trdr(is).lt.teeny) then
      dso4d=0.
      go to 22
      endif

c check to make sure no overreaction: moles of production:
      dso4dt=dso4d*trdr(ih)*trdr(is)
c can't be more than moles going in:
      if (dso4dt.gt.trdmol(is)) then
        dso4d=trdmol(is)/(trdr(ih)*trdr(is))
      endif
      dso4dt=dso4d*trdr(ih)*trdr(is)
      if (dso4dt.gt.trdmol(ih)) then
        dso4d=trdmol(ih)/(trdr(ih)*trdr(is))
      endif

 22   continue
      do n=1,ntx
       select case (trname(ntix(n)))
       case('SO4','M_ACC_SU')
       is4=ntix(n)
!       is4x=n
       sulfout(is4)=tr_mm(is4)/1000.*(dso4g*tm(l,isx)*tm(l,ihx)
     *  +dso4d*tmcl(isx)*tmcl(ihx)) !kg

       dt_sulf(is4) = dt_sulf(is4) + sulfout(is4)

       case('SO2')
       is=ntix(n)
       isx=n
! is ih/ihx set here, then why isn't is/isx?
       sulfin(is)=-dso4g*tm(l,ihx)*tr_mm(is)/1000. !dimnless
       sulfinc(is)=-dso4d*tmcl(ihx)*tr_mm(is)/1000.
       sulfinc(is)=max(-1d0,sulfinc(is))
       sulfin(is)=max(-1d0,sulfin(is))
       tr_left(isx)=0.
       if (fcloud.gt.abs(sulfin(is))) then
         tr_left(isx)=(fcloud+sulfin(is))
       endif

       dt_sulf(is)=dt_sulf(is)+sulfin(is)*tm(l,isx)+sulfinc(is)*trdr(is)

       case('H2O2','H2O2_s')

         if (trname(ntix(n)).eq."H2O2" .and. coupled_chem.eq.0) goto 401
         if (trname(ntix(n)).eq."H2O2_s" .and. coupled_chem.eq.1) goto
     *        401

       ih=ntix(n)
       ihx=n
       sulfin(ih)=-dso4g*tm(l,isx)*tr_mm(ih)/1000.
       sulfinc(ih)=-dso4d*tmcl(isx)*tr_mm(ih)/1000.
       sulfinc(ih)=max(-1d0,sulfinc(ih))
       sulfin(ih)=max(-1d0,sulfin(ih))
       tr_left(ihx)=0.
       if (fcloud.gt.abs(sulfin(ih))) then
         tr_left(ihx)=fcloud+sulfin(ih)
       endif


 401   CONTINUE

       dt_sulf(ih)=dt_sulf(ih)+sulfin(ih)*tm(l,ihx)+sulfinc(ih)*trdr(ih)

      end select
      END DO

 333  RETURN
      END SUBROUTINE GET_SULFATE

      SUBROUTINE GET_BC_DALBEDO(i,j,bc_dalb)
!@sum Calculates change to albedo of snow on ice and snow on land due
!@+     to BC within the snow.
!@+     Parameterization based on Warren and Wiscombe (1980) (21 inputs)
!@+     or actually on Flanner et al. fig 2 r_e=500 (14 inputs)
!@+     or Warren and Wiscombe (1985) (18 input, old vs new, then I
!@+     continue linearly from 19-29 off the plot)
!@+auth Dorothy Koch
c
      USE CONSTANT, only: rhow
      USE GHY_COM, only: tr_wsn_ij, wsn_ij
      USE SEAICE_COM, only: trsi, snowi
      USE TRACER_COM, only: trname,n_BCB,n_BCII,n_BCIA
      USE VEG_COM, only: afb
      USE RADPAR, only: agesn
      USE FLUXES, only: gtracer,gtemp
      IMPLICIT NONE
c Warren and Wiscombe 1985 includes age dependence
      real*8, parameter :: bc(29)=(/1.d0,2.d0,3.d0,4.d0,5.d0,
     * 6.d0,7.d0,8.d0,9.d0,10.d0,20.d0,30.d0,40.d0,50.d0,60.d0,
     * 70.d0,80.d0,90.d0,100.d0,110.d0,120.d0,130.d0,140.d0,
     * 150.d0,160.d0,170.d0,180.d0,190.d0,200.d0/)
      real*8, parameter :: daln(29)=(/0.d0,0.1d0,0.1d0,0.2d0,
     * 0.2d0,0.2d0,0.2d0,0.3d0,0.3d0,0.4d0,0.7d0,0.9d0,1.1d0,
     * 1.3d0,1.5d0,1.6d0,1.8d0,2.d0,2.2d0,2.4d0,2.6d0,2.8d0,
     * 3.d0,3.2d0,3.4d0,3.6d0,3.8d0,4.d0,4.2d0/)
      real*8, parameter :: dalo(29)=(/0.1d0,0.2d0,0.4d0,0.5d0,
     * 0.6d0,0.7d0,0.8d0,0.9d0,1.d0,1.d0,2.d0,2.6d0,3.2d0,
     * 3.8d0,4.3d0,4.8d0,5.2d0,5.5d0,5.9d0,6.3d0,6.7d0,7.1d0,
     * 7.5d0,7.9d0,8.3d0,8.7d0,9.1d0,9.5d0,9.9d0/)
c Flanner et al
c     real*8, parameter :: bc(14)=(/25.d0,50.d0,100.d0,150.d0,
c    * 200.d0,
c    *250.d0,300.d0,400.d0,500.d0,600.d0,700.d0,800.d0,900.d0,
c    * 1000.d0/)
c     real*8, parameter :: dal(14)=(/1.d0,2.d0,3.d0,4.d0,5.d0,
c    * 6.d0,7.d0,8.d0,9.d0,10.d0,11.d0,11.5d0,12.d0,12.5d0/)
c Warren and Wiscomb 1980
c     real*8, parameter :: bc(21)=(/0.05d0,0.075d0,0.1d0,0.2d0,
c    * 0.3d0,0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0,1.d0,2.d0,
c    * 3.d0,4.d0,5.d0,6.d0,7.d0,8.d0,9.d0,10.d0/)
c     real*8, parameter :: dal(21)=(/2.d0,3.d0,4.d0,6.d0,8.d0,
c    * 9.d0,10.d0,11.d0,12.d0,13.d0,14.d0,16.d0,18.d0,20.d0,
c    * 22.d0,24.d0,26.d0,28.d0,30.d0,32.d0,34.d0/)
      REAL*8  scon,icon,bcsnowb,bcsnowv,bcice,fv,fb,
     * sconb,sconv,bcc,rads
      INTEGER n,ic,ib
      INTEGER, INTENT(IN) :: i,j
      REAL*8, INTENT(OUT) :: bc_dalb
c
c tr_wsn_ij(n,nsl,2,i,j) tracer in snow layer l multiplied by fraction snow, kg/m2
c wsn_ij(nsl,2,i,j)
c trsi(n,nsi,i,j) tracer in sea ice in layer l, kg/m2
c snowi(i,j) snow amount on sea ice, kg/m2
c afb(i,j)=fb, fraction that is bare soil
c fv=1-fb fraction that is vegetated
c rads is the snow grain size determined in GRAINS
c gtracer(n,2,i,j) is tracer concentration in snow on sea ice?
c Maybe I need tracer in snow on sea ice, or mass of sea ice...?
c Does trsi accumulate for ALL tracers?
c fractions??
c
      bc_dalb=0.
      scon=0.
      icon=0.
      bcsnowb=0.
      bcsnowv=0.
      bcice=0.
      fb=afb(i,j)
      fv=1.-fb
       if (wsn_ij(1,1,i,j).gt.0.)  then
       bcsnowb=tr_wsn_ij(n_BCII,1,1,i,j)+
     *  tr_wsn_ij(n_BCIA,1,1,i,j)+tr_wsn_ij(n_BCB,1,1,i,j)
       sconb=bcsnowb/wsn_ij(1,1,i,j)/rhow
       endif
       if (wsn_ij(1,2,i,j).gt.0.)  then
       bcsnowv=tr_wsn_ij(n_BCII,1,2,i,j)+
     *  tr_wsn_ij(n_BCIA,1,2,i,j)+tr_wsn_ij(n_BCB,1,2,i,j)
       sconv=bcsnowv/wsn_ij(1,2,i,j)/rhow
       endif
       scon=(fb*sconb+fv*sconv)*1.D9   !kg/kg to ppmw
       if (snowi(i,j).gt.0.) then
       icon=(gtracer(n_BCII,2,i,j)+gtracer(n_BCIA,2,i,j)+
     *  gtracer(n_BCB,2,i,j))*1.d9
       endif
       bcc=DMAX1(icon,scon)
       call GRAINS(i,j,rads)

       do ib=1,28
       if (bcc.gt.bc(ib).and.bcc.lt.bc(ib+1)) then
       bc_dalb=-(daln(ib)
     *       +(rads-100.d0)/900.d0*(dalo(ib)-daln(ib)))/100.
       go to 33
       endif
       end do
 33    continue
       if (bcc.ge.bc(29)) bc_dalb=-(daln(29)
     *    + (rads-100.d0)/900.d0*(dalo(29)-daln(29)))/100.
c     if (bc_dalb.ne.0.) write(6,*) 'alb_write',i,j,bc_dalb,bcc,rads
      RETURN
      END SUBROUTINE GET_BC_DALBEDO

      SUBROUTINE GRAINS(i,j,rads)
!@sum Estimates snow grain size (microns) based on air temperature
!@+     and snow age. From Susan Marshall's PhD thesis
!@+auth Dorothy Koch
c
      USE CONSTANT, only: pi,gasc
      USE PBLCOM, only: tsavg
      USE GHY_COM, only: snoage
      USE AEROSOL_SOURCES, only: snosiz
      IMPLICIT none
      REAL*8 E,A,age,r0,radmm,ert,
     * tfac,area,delrad
      INTEGER n,ic,ib
      INTEGER, INTENT(IN) :: i,j
      REAL*8, INTENT(OUT) :: rads
      DATA E,A /26020.d0, 29100.d0/
c tsavg(i,j) surface air temperature
c snoage(k,i,j) k=1 ocean ice, k=2 land ice, k=3 land
c   (do these differ within a gridbox?) in days
c RN1 radius from previous timestep
c RADS snow grain radius
c
c Find the age of snow, I assume the age does not
c  vary within the gridbox so just take the max?
       age=DMAX1(snoage(1,i,j),snoage(2,i,j),snoage(3,i,j))
c Use Temperature to check if melting or non-melting snow
       IF (tsavg(i,j).le.273.15) then
c Non-melting snow; distinguish between initial or
c  secondary growth rate
        IF (age.lt.13.5) then
c initial growth
         r0=50.
         radmm=r0+ (0.008d0+age)
         rads=radmm*1000.d0
         snosiz(i,j)=rads+r0
        ELSE
c secondary growth
         r0=150.d0
         ert = dexp(-E/(GASC*TSAVG(I,J)))
         tfac = a*ert
         area = (TFAC/365.d0) * (AGE-12.5d0)
         radmm=dsqrt(area/pi)
         delrad = radmm*1000.d0
         rads=delrad + r0
         snosiz(i,j)=rads
        ENDIF
       ELSE
c melting snow
        radmm=dsqrt(0.05d0)
        rads=snosiz(i,j)+(radmm*100.d0)
        snosiz(i,j)=rads
       ENDIF
       rads=DMIN1(rads,1000.d0)
       rads=DMAX1(rads,100.d0)
      RETURN
      END SUBROUTINE GRAINS
      
            SUBROUTINE get_aircraft_SO2
!@sum  get_aircraft_SO2 to define the 3D source of SO2 from aircraft
!@auth Drew Shindell? / Greg Faluvegi / Jean Learner
!@ver  1.0 (based on DB396Tds3M23)
!! taken from get_aircraft_NOx by Shindell/Faluvegi
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only: itime,jday,JDperY,im,jm,lm,ptop,psf,sig
      USE DOMAIN_DECOMP, only: GRID, GET, write_parallel
      USE CONSTANT, only: sday,hrday
      USE FILEMANAGER, only: openunit,closeunit
      USE FLUXES, only: tr3Dsource
      USE GEOM, only: dxyp
      USE TRACER_COM, only: itime_tr0,trname,n_SO2,n_BCIA
      use AEROSOL_SOURCES, only: Laircrs,aircrafts_Tyr1,aircrafts_Tyr2
      use param, only: sync_param
C
      IMPLICIT NONE
c
!@var nanns,nmons: number of annual and monthly input files
!@var l,ll dummy loop variable
!@var pres local pressure variable
      integer, parameter :: nanns=0,nmons=1
      integer, dimension(nmons) :: mon_units, imon
      integer l,i,j,iu,k,ll,nt,ns,nsect,nn
      character*80 :: title
      character*12, dimension(nmons) :: mon_files=(/'SO2_AIRCRAFT'/)
      character(len=300) :: out_line
      character*124 :: tr_sectors_are
      character*32 :: pname
      logical :: LINJECTS
      logical, dimension(nmons) :: mon_bins=(/.true./) ! binary file?
      real*8 bySperHr
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,Laircrs,1)
     &                                     :: src
      REAL*8, DIMENSION(LM)                :: pres
      REAL*4, PARAMETER, DIMENSION(Laircrs) :: PAIRL =
     & (/1013.25,898.74,794.95,701.08,616.40,540.19,471.81,
     &    410.60,355.99,307.42,264.36,226.32,193.30,165.10,
     &    141.01,120.44,102.87,87.866,75.048/)
      INTEGER :: J_1, J_0, J_0H, J_1H
c
C**** Local parameters and variables and arguments:
      logical :: trans_emis=.false.
      integer :: yr1=0, yr2=0
C**** Aircraft NOx source input is monthly, on 19 levels. Therefore:
C     19x12=228 records. Read it in here and interpolated each day.

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)
      call GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

C****
C**** Monthly sources are interpolated to the current day
C**** The titles of the input files say this is in KG/m2/s, so no
C**** conversion is necessary:
C****
      call sync_param("aircrafts_Tyr1",aircrafts_Tyr1)
      call sync_param("aircrafts_Tyr2",aircrafts_Tyr2)
      k = 1
      if(aircrafts_Tyr1==aircrafts_Tyr2)then
        trans_emis=.false.; yr1=0; yr2=0
      else
        trans_emis=.true.; yr1=aircrafts_Tyr1; yr2=aircrafts_Tyr2
      endif
      call openunit(mon_files(k),mon_units(k),mon_bins(k))
      call read_mon_3D(Laircrs,mon_units(k),
     & src(:,:,:,k),trans_emis,yr1,yr2)
      call closeunit(mon_units(k))
C====
C====   Place aircraft sources onto model levels:
C====
      tr3Dsource(:,J_0:J_1,:,2,n_SO2) = 0.d0
      tr3Dsource(:,J_0:J_1,:,2,n_BCIA) = 0.d0
      PRES(:)=SIG(:)*(PSF-PTOP)+PTOP
      DO J=J_0,J_1
       DO I=1,IM
c multiply N by 3.3 to get NOx
c divide NOx emission by 35 to get S, *2 to get SO2
c divide NOx emission by 350 to get BC
        tr3Dsource(i,j,1,2,n_SO2) = SRC(I,J,1,1)*dxyp(j)*3.3d0/35.d0
     *  *2.d0
        tr3Dsource(i,j,1,2,n_BCIA) = SRC(I,J,1,1)*dxyp(j)*3.3d0/350.d0
        DO LL=2,Laircrs
          LINJECTS=.TRUE.
          DO L=1,LM
           IF(PAIRL(LL) > PRES(L).AND.LINJECTS) THEN
             tr3Dsource(i,j,l,2,n_SO2) =
     &  tr3Dsource(i,j,l,2,n_SO2) + SRC(I,J,LL,1)*dxyp(j)/35.d0*2.d0
     *       *3.3d0
             tr3Dsource(i,j,l,2,n_BCIA) =
     &  tr3Dsource(i,j,l,2,n_BCIA) + SRC(I,J,LL,1)*dxyp(j)/350.d0
     *       *3.3d0
             LINJECTS=.FALSE.
           ENDIF
          ENDDO ! L
        ENDDO   ! LL
       END DO   ! I
      END DO    ! J
      return
      END SUBROUTINE get_aircraft_SO2
      
      SUBROUTINE read_mon_3D
     & (Ldim,iu,data1,trans_emis,yr1,yr2)
!@sum Read in monthly sources and interpolate to current day
!@auth Jean Lerner and others / Greg Faluvegi
! taken from TRACERS_SPECIAL_Shindell, in case we
!  we run aerosols independent of gases
      USE MODEL_COM, only: jday,jyear,im,jm,idofm=>JDmidOfM
      USE TRACER_COM, only: kstep
      USE FILEMANAGER, only : NAMEUNIT
      USE DOMAIN_DECOMP, only : GRID,GET,READT_PARALLEL,REWIND_PARALLEL
     & ,write_parallel,backspace_parallel
      implicit none
!@var Ldim how many vertical levels in the read-in file?
!@var L dummy vertical loop variable
      integer Ldim,L,imon,iu,ipos,k,nn
      character(len=300) :: out_line
      real*8 :: frac, alpha
      real*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     A2D,B2D,dummy
      real*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,Ldim) ::
     &     tlca,tlcb,data1,sfc_a,sfc_b
      logical, intent(in):: trans_emis
      integer, intent(in):: yr1,yr2

      integer :: J_0, J_1

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

C No doubt this code can be combined/compressed, but I am going to
C do the transient and non-transient cases separately for the moment:

! -------------- non-transient emissions ----------------------------!
      if(.not.trans_emis) then
C
      imon=1                ! imon=January
      if (jday <= 16)  then ! JDAY in Jan 1-15, first month is Dec
        write(6,*) 'Not using this first record:'
        call readt_parallel(grid,iu,nameunit(iu),0,dummy,Ldim*11)
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),0,A2D,1)
          tlca(:,J_0:J_1,L)=A2D(:,J_0:J_1)
        enddo
        call rewind_parallel(iu)
      else              ! JDAY is in Jan 16 to Dec 16, get first month
        do while(jday > idofm(imon) .AND. imon <= 12)
          imon=imon+1
        enddo
        write(6,*) 'Not using this first record:'
        call readt_parallel(grid,iu,nameunit(iu),0,dummy,Ldim*(imon-2))
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),0,A2D,1)
          tlca(:,J_0:J_1,L)=A2D(:,J_0:J_1)
        enddo
        if(imon==13) call rewind_parallel(iu)
      end if
      do L=1,Ldim
        call readt_parallel(grid,iu,nameunit(iu),0,B2D,1)
        tlcb(:,J_0:J_1,L)=B2D(:,J_0:J_1)
      enddo
c**** Interpolate two months of data to current day
      frac = float(idofm(imon)-jday)/(idofm(imon)-idofm(imon-1))
      data1(:,J_0:J_1,:) =
     & tlca(:,J_0:J_1,:)*frac + tlcb(:,J_0:J_1,:)*(1.-frac)
      write(out_line,*) '3D source monthly factor=',frac
      call write_parallel(trim(out_line))

! --------------- transient emissions -------------------------------!
      else
        ipos=1
        alpha=0.d0 ! before start year, use start year value
        if(jyear>yr2.or.(jyear==yr2.and.jday>=183))then
          alpha=1.d0 ! after end year, use end year value
          ipos=(yr2-yr1)/kstep
        endif
        do k=yr1,yr2-kstep,kstep
          if(jyear>k .or. (jyear==k.and.jday>=183)) then
            if(jyear<k+kstep .or. (jyear==k+kstep.and.jday<183))then
              ipos=1+(k-yr1)/kstep ! (integer artithmatic)
              alpha=(365.d0*(0.5+real(jyear-1-k))+jday) /
     &              (365.d0*real(kstep))
              exit
            endif
          endif
        enddo
!
! read the two necessary months from the first decade:
!
      imon=1                ! imon=January
      if (jday <= 16)  then ! JDAY in Jan 1-15, first month is Dec
        write(6,*) 'Not using this first record:'
        call readt_parallel
     &  (grid,iu,nameunit(iu),0,dummy,(ipos-1)*12*Ldim+Ldim*11)
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),0,A2D,1)
          tlca(:,J_0:J_1,L)=A2D(:,J_0:J_1)
        enddo
        do nn=1,12*Ldim; call backspace_parallel(iu); enddo
      else              ! JDAY is in Jan 16 to Dec 16, get first month
        do while(jday > idofm(imon) .AND. imon <= 12)
          imon=imon+1
        enddo
        write(6,*) 'Not using this first record:'
        call readt_parallel
     &  (grid,iu,nameunit(iu),0,dummy,(ipos-1)*12*Ldim+Ldim*(imon-2))
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),0,A2D,1)
          tlca(:,J_0:J_1,L)=A2D(:,J_0:J_1)
        enddo
        if(imon==13)then
          do nn=1,12*Ldim; call backspace_parallel(iu); enddo
        endif
      end if
CCCCC write(6,*) 'Not using this first record:'
CCCCC call readt_parallel(grid,iu,nameunit(iu),0,dummy,Ldim*(imon-1))
      do L=1,Ldim
        call readt_parallel(grid,iu,nameunit(iu),0,B2D,1)
        tlcb(:,J_0:J_1,L)=B2D(:,J_0:J_1)
      enddo
      frac = float(idofm(imon)-jday)/(idofm(imon)-idofm(imon-1))
      sfc_a(:,J_0:J_1,:) =
     & tlca(:,J_0:J_1,:)*frac + tlcb(:,J_0:J_1,:)*(1.-frac)
      call rewind_parallel( iu )

      ipos=ipos+1
      imon=1                ! imon=January
      if (jday <= 16)  then ! JDAY in Jan 1-15, first month is Dec
        write(6,*) 'Not using this first record:'
        call readt_parallel
     &  (grid,iu,nameunit(iu),0,dummy,(ipos-1)*12*Ldim+Ldim*11)
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),0,A2D,1)
          tlca(:,J_0:J_1,L)=A2D(:,J_0:J_1)
        enddo
        do nn=1,12*Ldim; call backspace_parallel(iu); enddo
      else              ! JDAY is in Jan 16 to Dec 16, get first month
        do while(jday > idofm(imon) .AND. imon <= 12)
          imon=imon+1
        enddo
        write(6,*) 'Not using this first record:'
        call readt_parallel
     &  (grid,iu,nameunit(iu),0,dummy,(ipos-1)*12*Ldim+Ldim*(imon-2))
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),0,A2D,1)
          tlca(:,J_0:J_1,L)=A2D(:,J_0:J_1)
        enddo
        if(imon==13)then
          do nn=1,12*Ldim; call backspace_parallel(iu); enddo
        endif
      end if
CCCCCCwrite(6,*) 'Not using this first record:'
CCCCCCcall readt_parallel(grid,iu,nameunit(iu),0,dummy,Ldim*(imon-1))
      do L=1,Ldim
        call readt_parallel(grid,iu,nameunit(iu),0,B2D,1)
        tlcb(:,J_0:J_1,L)=B2D(:,J_0:J_1)
      enddo
      frac = float(idofm(imon)-jday)/(idofm(imon)-idofm(imon-1))
      sfc_b(:,J_0:J_1,:) =
     & tlca(:,J_0:J_1,:)*frac + tlcb(:,J_0:J_1,:)*(1.-frac)

! now interpolate between the two time periods:

      data1(:,J_0:J_1,:) =
     & sfc_a(:,J_0:J_1,:)*(1.d0-alpha) + sfc_b(:,J_0:J_1,:)*alpha

      write(out_line,*) '3D source at',
     &100.d0*alpha,' % of period ',k,' to ',k+kstep,
     &' and monthly fraction= ',frac
      call write_parallel(trim(out_line))

      endif ! transient or not

      return
      end subroutine read_mon_3D
      

