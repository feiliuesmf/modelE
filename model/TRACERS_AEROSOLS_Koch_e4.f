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
     * ,craft,so2t_src
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

      return
      end subroutine alloc_aerosol_sources      
      
      subroutine get_O3_offline
!@sum read in ozone fields for aqueous oxidation
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only: jday,im,jm,lm,ptop,psf,sig
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits      
      USE AEROSOL_SOURCES, only: o3_offline
C
      IMPLICIT NONE

C**** Local parameters and variables and arguments:
c
!@var nanns,nmons: number of annual and monthly input files
      integer, parameter :: nanns=0,nmons=1,levo3=23
      integer ann_units(nanns),mon_units(nmons),imon(nmons)
      integer i,j,iu,k,l
      integer      :: jdlast=0  
      logical :: ifirst=.true.
      character*80 title   
      character*10 :: mon_files(nmons) = (/'O3_FIELD'/)
      logical      :: mon_bins(nmons)=(/.true./) ! binary file?
      REAL*8, DIMENSION(IM,JM,levo3,1) :: src
      REAL*8,DIMENSION(Im,Jm,levo3),SAVE :: tlca,tlcb
      save jdlast,mon_units,imon,ifirst
C initialise
c      o3_offline(:,:,:)=0.0d0
c
C     Read it in here and interpolated each day.
C
      if (ifirst) call openunits(mon_files,mon_units,mon_bins,nmons)
      IF (ifirst) jdlast=jday-1
      j = 0
      do k=nanns+1,1
        j = j+1
        call read_monthly_O3_3D(levo3,mon_units(j),jdlast,
     *       tlca,tlcb,src(1,1,1,k),imon(j),ifirst)
      end do
      ifirst=.FALSE.
      jdlast = jday

      call sys_flush(6)

      do k=nanns+1,1; DO l=1,lm; DO J=1,JM; DO I=1,IM
      o3_offline(i,j,l)=src(I,J,L,k)

      end do
      end do
      end do
      end do

      END SUBROUTINE get_O3_offline

      SUBROUTINE read_monthly_O3_3D(Ldim,iu,jdlast,tlca,tlcb,
     *     data1,imon,ifirst)
!@sum Read in monthly sources and interpolate to current day
!@ Author Greg Faluvegi 
!@+   Calling routine must have the lines:
!@+      real*8 tlca(im,jm,Ldim,nm),tlcb(im,jm,Ldim,nm)
!@+      integer imon(nm)   ! nm=number of files that will be read
!@+      save jdlast,tlca,tlcb,imon
!@+   Input: iu, the fileUnit#; jdlast
!@+   Output: interpolated data array + two monthly data arrays
      USE MODEL_COM, only: jday,im,jm,idofm=>JDmidOfM
      implicit none

!@var Ldim how many vertical levels in the read-in file?
      INTEGER,INTENT(IN) :: Ldim,iu,jdlast
      LOGICAL,INTENT(IN) :: ifirst

      INTEGER,INTENT(INOUT) :: imon
      REAL*8,DIMENSION(Im,Jm,Ldim),INTENT(INOUT) :: tlca,tlcb

      REAL*8,INTENT(OUT) :: data1(im,jm,Ldim)

!@var L dummy vertical loop variable
      INTEGER :: L
!@var imom saves month for interpolation
      INTEGER,SAVE :: imom
!@var frac weighting factor for interpolation between months
!@var A2D,B2D 2 dimensional fields to read in data
      REAL*8 :: frac,A2D(im,jm),B2D(im,jm)

      idofm(2)=46

      IF (ifirst) THEN          ! NEED TO READ IN FIRST MONTH OF DATA AT START
        REWIND iu
        imon=1                  ! imon=January
        if (jday < 16)  then    ! JDAY in Jan 1-15, first month is Dec
          do L=1,LDim*11
            read(iu)
          end do
          DO L=1,Ldim
            call readt(iu,0,A2D,im*jm,A2D,1) ! read in December
            tlca(:,:,L)=A2d(:,:)
          END DO
          rewind iu
        else                    ! JDAY is in Jan 16 to Dec 16, get first month
          imon=1
          DO WHILE (jday >= idofm(imon) .AND. imon <= 12)
            imon=imon+1
          END DO
          do L=1,Ldim*(imon-2)
            read(iu)
          end do
          DO L=1,Ldim
            call readt(iu,0,A2D,im*jm,A2D,1) ! read in current month
            tlca(:,:,L)=A2d(:,:)
          END DO
          IF (imon == 13) REWIND iu
        end if
      END IF
c ..........
c Read in next month at start or when middle of the current month is reached
c ..........
      IF (ifirst .OR. (jday == idofm(imon) .AND. jday /= jdlast)) THEN
        IF (.NOT. ifirst) THEN
          tlca(:,:,:) = tlcb(:,:,:)
          IF (imon == 12) REWIND iu
        END IF
        DO L=1,Ldim
          CALL readt(iu,0,B2D,im*jm,B2D,1) ! read in next month
          tlcb(:,:,L)=B2D(:,:)
        END DO
        IF (jday >= idofm(imon)) imon=imon+1
        imom=imon
        IF (imon == 13) imon=1
      END IF
      IF (jday == 1) imom=1

c**** Interpolate two months of data to current day
      frac = float(idofm(imom)-jday)/(idofm(imom)-idofm(imom-1))
      data1(:,:,:) = tlca(:,:,:)*frac + tlcb(:,:,:)*(1.-frac)

      write(6,*) 'Read in ozone offline fields for in cloud oxidation 
     *     interpolated to current day',frac

      END SUBROUTINE read_monthly_O3_3D

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
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits,
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
      logical :: ann_bins(nanns)=(/.true.,.true.,
     *   .true.,.true./) ! binary file?
      character*12 :: mon_files(nmons) =
     * (/'SO2_AGNB_E95','SO2_BBNB_E95'/)
      logical :: mon_bins(nmons)=(/.true.,.true./)
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
        call openunits(ann_files,ann_units,ann_bins,nanns)
        k = 0
        do iu = ann_units(1),ann_units(nanns)
          k = k+1
c         call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          call readt_parallel (grid,iu,nameunit(iu),0,src(:,:,k),1)
          do j=j_0,j_1
            src(:,j,k) = src(:,j,k)*bydxyp(j)*bySperHr
          end do
        end do
        call closeunits(ann_units,nanns)
       
        call openunits(mon_files,mon_units,mon_bins,nmons)
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
      call openunit('SO2_BIOMASS',iuc,.false.)
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
      call openunit('BC_BIOMASS',iuc,.false.)
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
      call openunit('BC_BM_RAT',iuc,.false.)
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
      call openunit('OC_BIOMASS',iuc,.false.)
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
      call openunit('OC_BM_RAT',iuc,.false.)
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
       if (ihyr.le.1990.and.imAER.eq.3) then
       tfac=0.5d0*(1.d0+real(ihyr-1875)/115.d0)
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
      real*8, dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,8) :: 
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
      if (ihyr.ge.1875.and.ihyr.lt.1900) then
      jb1=1875
      jb2=1900
      irr=1
      else if (ihyr.ge.1900.and.ihyr.lt.1925) then
      jb1=1900
      jb2=1925
      irr=2
      else if (ihyr.ge.1925.and.ihyr.lt.1950) then
      jb1=1925
      jb2=1950
      irr=3
      else if (ihyr.ge.1950.and.ihyr.lt.1960) then
      jb1=1950
      jb2=1960
      irr=4
      else if (ihyr.ge.1960.and.ihyr.lt.1970) then
      jb1=1960
      jb2=1970
      irr=5
      else if (ihyr.ge.1970.and.ihyr.lt.1980) then
      jb1=1970
      jb2=1980
      irr=6
      else if (ihyr.ge.1980.and.ihyr.le.1990) then
      jb1=1980
      jb2=1990
      irr=7
      endif
      call openunit('BC_INDh',iuc, .false.)
      do ip=1,7326
      read(iuc,*) ii,jj,iy,xbcff,xbcbm,xombm
      if (jj<j_0 .or. jj>j_1) cycle
      hbc_all(ii,jj,iy)=xbcff+xbcbm
      hoc_all(ii,jj,iy)=xbcff*2.d0+xombm
      end do
      call closeunit(iuc)
      hbc(:,j_0:j_1,1:2)=hbc_all(:,j_0:j_1,irr:irr+1)
      hoc(:,j_0:j_1,1:2)=hoc_all(:,j_0:j_1,irr:irr+1)
c kg/year to kg/s
      hbc(:,j_0:j_1,1:2)=hbc(:,j_0:j_1,1:2)/(365.d0*24.d0*3600.d0)
      hoc(:,j_0:j_1,1:2)=hoc(:,j_0:j_1,1:2)/(365.d0*24.d0*3600.d0)
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
c historic SO2 emissions
      USE MODEL_COM, only: jyear, jday,im,jm
      USE DOMAIN_DECOMP, only : GRID, GET,AM_I_ROOT
     * ,UNPACK_DATA
      USE FILEMANAGER, only: openunit,closeunit
      USE TRACER_COM, only: aer_int_yr
      USE AEROSOL_SOURCES, only: SO2_src,hso2
      USE GEOM, only: dxyp
      implicit none
      integer iuc,irr,ihyr,i,j,id,jb1,jb2,iact,j_0,j_1,nn,j_0h,j_1h
      real*8, DIMENSION(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,8) :: 
     * hso2_all
c     real*4, DIMENSION(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,8) :: 
c    *  hso2_all_read
      real*8, DIMENSION(im,jm,8) :: hso2_all_glob
      real*4, DIMENSION(im,jm,8) :: hso2_all_read
      real*8 d1,d2,d3
      save jb1,jb2

      CALL GET(grid, J_STRT=J_0,       J_STOP=J_1,
     *      J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      
      SO2_src(:,:,:)=0.d0

      if (aer_int_yr.eq.0) then
      ihyr=jyear
      else
      ihyr=aer_int_yr
      endif
c if run just starting or if it is a new year
c   then open new files
      if (iact.eq.0.or.jday.eq.1) then
      if (ihyr.ge.1875.and.ihyr.lt.1900) then
      jb1=1875
      jb2=1900
      irr=1
      else if (ihyr.ge.1900.and.ihyr.lt.1925) then
      jb1=1900
      jb2=1925
      irr=2
      else if (ihyr.ge.1925.and.ihyr.lt.1950) then
      jb1=1925
      jb2=1950
      irr=3
      else if (ihyr.ge.1950.and.ihyr.lt.1960) then
      jb1=1950
      jb2=1960
      irr=4
      else if (ihyr.ge.1960.and.ihyr.lt.1970) then
      jb1=1960
      jb2=1970
      irr=5
      else if (ihyr.ge.1970.and.ihyr.lt.1980) then
      jb1=1970
      jb2=1980
      irr=6
      else if (ihyr.ge.1980.and.ihyr.le.1990) then
      jb1=1980
      jb2=1990
      irr=7
      endif
      if(AM_I_ROOT( ))then      
      call openunit('SO2_INDh',iuc, .true.)
      read(iuc) hso2_all_read
      call closeunit(iuc)
      endif
      hso2_all_glob(:,:,:)=hso2_all_read(:,:,:)
      call UNPACK_DATA( grid, hso2_all_glob, hso2_all )
      hso2(:,j_0:j_1,1:2)=hso2_all(:,j_0:j_1,irr:irr+1)
c kg S/m2/year to kg SO2/s
      do j=j_0,j_1
      hso2(:,j,1:2)=hso2(:,j,1:2)
     *    *dxyp(j)*2.d0
     *    /(365.d0*24.d0*3600.d0)
      end do
      endif
c interpolate to model year
      d1=real(ihyr-jb1)
      d2=real(jb2-ihyr)
      d3=real(jb2-jb1)
      SO2_src(:,:,1)=(d1*hso2(:,:,2)+d2*hso2(:,:,1))/d3

      end subroutine read_hist_SO2

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
      call openunit('SO2_BIOMASS',iuc2,.false.)
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
        tr3Dsource(:,j_0:j_1,:,1,n_H2SO4)=0. ! H2O2 chem sink
        tr3Dsource(:,j_0:j_1,:,1,n_M_ACC_SU)=0. ! ACC SO4 source
        tr3Dsource(:,j_0:j_1,:,1,n_M_AKK_SU)=0. ! AKK SO4 source
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
        call openunit('AER_CHEM',iuc,.true.)
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
        call openunit('AER_OH_STRAT',iuc2,.true.)
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
      CALL GET_O3_OFFLINE
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
      o3mc = o3_offline(i,j,l)*dmm*(28.0D0/48.0D0)*BYDXYP(J)*BYAM(L,I,J)
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
#ifdef TRACERS_SPECIAL_Shindell
          najl = jls_OHcon
#else
          najl = jls_OHconk
#endif
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
     *     ,trname,ntm,tr_mm,lm,n_SO4,n_H2O2
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
#ifdef TRACERS_COSMO
      do n=1,ntx
       select case (trname(ntix(n)))
       case ('Pb210   ','Be7     ','Be10    ','Rn222')
       go to 333
       end select
      end do
#endif
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
      pph(is)=mair*1.d-3*ppas/tr_mm(is)/amass*
     *   tr_rkd(is)*exp(-tr_dhd(is)*tfac)
c the following is from Phil:
c      reduction in partial pressure as species dissolves
      henry_const(is)=rkdm(is)*exp(-tr_dhd(is)*tfac)
      pph(is)=pph(is)/(1+(henry_const(is)*clwc*gasc*temp))
c again all except tmcl(n)
      trdr(is)=mair*ppas/tr_mm(is)/amass*bygasc
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
      pph(ih)=mair*1.D-3*ppas/tr_mm(ih)/amass*
     *   tr_rkd(ih)*exp(-tr_dhd(ih)*tfac)
c the following is from Phil:
c      reduction in partial pressure as species dissolves
      henry_const(ih)=rkdm(ih)*exp(-tr_dhd(ih)*tfac)
      pph(ih)=pph(ih)/(1+(henry_const(ih)*clwc*gasc*temp))
c all except tmcl(n)
      trdr(ih)=mair*ppas/tr_mm(ih)
     *    /amass*bygasc/temp*1.D-3  !M/kg
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


