#include "rundeck_opts.h"

      MODULE TRACER_SOURCES

#ifdef TRACERS_ON

      USE TRACER_COM
      USE MODEL_COM, only : nisurf

      IMPLICIT NONE
      SAVE

!@param Laircr the number of layers of aircraft data read from file
!@param Lsulf the number of layers of sulfate SA data read from file
!@param correct_CO_ind correction factor Lee put in for industrial CO
      INTEGER, PARAMETER :: 
     &                      Laircr      =19,
     &                      Lsulf       =23 ! not LM
      REAL*8, PARAMETER :: correct_CO_ind=520./410.
c
!@var CH4_src           CH4 surface sources and sinks (kg/m2/s)
!@var CO_src             CO surface sources and sinks (kg/m2/s)
!@var Alkenes_src   Alkenes surface sources and sinks (kg/m2/s)
!@var Paraffin_src Paraffin surface sources and sinks (kg/m2/s)
!@var NOx_src           NOx surface sources and sinks (kg/m2/s)
!@var Isoprene_src ISoprene surface sources and sinks (kg/m2/s)
c
c SHOULD PROBABLY USE ntsurfsrc( ) instead of these ...
#ifdef EDGAR_HYDE_SOURCES
      integer, parameter :: nch4src     = 15,
     &                      nCOsrc      =  6,
     &                      nAlkenessrc =  8,
     &                      nParaffinsrc=  8,
     &                      nIsoprenesrc=  1,
     &                      nNOxsrc     =  7
#else
      integer, parameter :: nch4src     = 14,
     &                      nCOsrc      =  2,
     &                      nAlkenessrc =  3,
     &                      nParaffinsrc=  3,
     &                      nIsoprenesrc=  1,
     &                      nNOxsrc     =  3
#endif
      real*8 CO_src(im,jm,nCOsrc),CH4_src(im,jm,nch4src),
     & Alkenes_src(im,jm,nAlkenessrc),Paraffin_src(im,jm,nParaffinsrc),
     & Isoprene_src(im,jm,nIsoprenesrc),NOx_src(im,jm,nNOxsrc)
#ifdef INTERACTIVE_WETLANDS_CH4
!@param nra_ch4 number of running averages needed for int-wetlands
!@param maxHR_ch4 maximum number of sub-daily accumulations
!@param fact_ncep for temp: nothing, for prec: 24 converts from model
!@+     units kg/m2/hr to mm/day
!@param nday_ncep number of days in running average (prec,temp)
!@var by_nday_ncep  1/real(nday_ncep)
!@var day_ncep daily NCEP precipitation (mm/day) (prec,temp)
!@var DRA_ch4 daily running average of model (prec,temp)
!@var avg_ncep running average (prec,temp) over nday_ncep days
!@var avg_modPT equivalent of avg_ncep, but based on model variables
!@var sum_ncep,day_ncep_tmp temp arrays for computing running average
!@var PRS_ch4 period running sum of model (prec,temp)
!@var HRA_ch4 hourly running average of model (prec,temp)
!@var iday_ncep current day (counter) of averaging period (prec,temp)
!@var iHch4 "hourly" index for averages of model (prec,temp)
!@var iDch4 "daily"  index for averages of model (prec,temp)
!@var i0ch4 ponter to current index in running sum of mode (prec,temp)
!@var first_ncep whether in the first ncep averaging period (prec,temp)
!@var first_mod whether in the first model averaging per.   (prec,temp)
      integer, parameter :: nra_ch4 = 2, maxHR_ch4=24*NIsurf
      real*8, parameter, dimension(nra_ch4)  :: fact=(/24.d0,1.d0/)
      integer, parameter, dimension(nra_ch4) :: nday_ncep = (/28,14/)
      real*8, dimension(nra_ch4)             :: by_nday_ncep
      real*8, dimension(im,jm,nday_ncep(1),nra_ch4):: day_ncep,DRA_ch4
      real*8, dimension(im,jm,nra_ch4)       :: avg_ncep,avg_modPT,
     &                                   sum_ncep,day_ncep_tmp,PRS_ch4
      real*8, dimension(im,jm,maxHR_ch4,nra_ch4)   :: HRA_ch4
      integer, dimension(nra_ch4)            :: iday_ncep, i0_ncep
      integer, dimension(im,jm,nra_ch4)      :: iHch4, iDch4, i0ch4  
      logical, dimension(nra_ch4) :: first_ncep = (/.true.,.true./)
     &                               first_mod  = (/.true.,.true./)
#endif
#endif
      END MODULE TRACER_SOURCES
      
      
      MODULE LIGHTNING
!@sum  LIGHTNING_COM model variables lightning parameterization
!@auth Colin Price (modelEification by Greg Faluvegi)
!@ver  1.0 (taken from CB436Tds3M23)

      USE MODEL_COM, only : IM,JM,LM,LS1

      IMPLICIT NONE
      SAVE
      
!@var JN J at 30 N
!@var JS J at 30 S
!@var I dummy

      INTEGER, PARAMETER :: JS = JM/3 + 1, JN = 2*JM/3
      INTEGER I
#ifdef TRACERS_ON
      REAL*8, DIMENSION(IM,JM) :: RNOx_lgt
      REAL*8, DIMENSION(LM) :: SRCLIGHT
#endif
      
!@var HGT Pickering vertical lightning distributions (1998)
      REAL*8 HGT(2,2,16)
      DATA (HGT(1,1,I),I=1,16)/8.2,1.9,2.1,1.6,1.1,1.6,3.0,5.8,
     &  7.6,9.6,10.5,12.3,11.8,12.5,8.1,2.3/
      DATA (HGT(1,2,I),I=1,16)/5.8,2.9,2.6,2.4,2.2,2.1,2.3,6.1,
     & 16.5,14.1,13.7,12.8,12.5,2.8,0.9,0.3/
      DATA (HGT(2,1,I),I=1,16)/20.1,2.3,0.8,1.5,3.4,5.3,3.6,3.8,
     &    5.4,6.6,8.3,9.6,12.8,10.0,6.2,0.3/
      DATA (HGT(2,2,I),I=1,16)/5.8,2.9,2.6,2.4,2.2,2.1,2.3,6.1,
     & 16.5,14.1,13.7,12.8,12.5,2.8,0.9,0.3/ 
 
      END MODULE LIGHTNING 
      
          
#ifdef TRACERS_ON
#ifdef EDGAR_HYDE_SOURCES
!
      subroutine read_CH4_sources(nt,iact)
!@sum reads in CH4 surface sources and sinks: EDGA-HYDE VERSION
!@auth Greg Faluvegi & Jean Lerner
C****
C**** There are 6 monthly sources and 9 annual sources
C**** Annual sources are read in at start and re-start of run only
C**** Monthly sources are interpolated each day
      USE MODEL_COM, only: itime,JDperY,im,jm,jday,focean,fearth,DTsrc
      USE GEOM, only: BYDXYP
      USE LAKES_COM, only: flake
      USE CONSTANT, only: sday
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      USE TRACER_COM, only: itime_tr0,trname
      use TRACER_SOURCES, only: src=>ch4_src,nsrc=>nch4src
      USE TRCHEM_Shindell_COM, only: PI_run,PIratio_indus,PIratio_bburn
 
      implicit none
      character*80 title
!@var adj Factors that tune the total amount of individual sources
      logical :: ifirst=.true.
      logical :: adjust(nsrc) = (/.false.,.false.,.false.,.false.,
     & .true.,.true.,.true.,.true.,.true.,.true.,.false.,.false.,
     & .false.,.false.,.false./)
      real*8 adj(nsrc)
Corig data adj/4*-1.d30, 0.999d0, 1.194d0, 3.5997d-5, 17.330d-4,
      data adj/4*-1.d30, 0.999d0, 1.588d0, 3.5997d-5, 17.330d-4,
     & 5.3558d-5, 0.9818d0, 5*-1.d30/
      INTEGER :: EHorJEAN(nsrc) = (/1,1,1,1,2,2,2,2,2,2,1,1,1,1,1/)
      real*8 ann_conv(2),mon_conv(2)!1=edgar-hyde conversion, 2=Jean's 
      real*8 byDTsrc
C     
C    K=1    Animals (edgar hyde annual)
C    K=2    Fossil Fuel Combustion (edgar hyde annual)
C    K=3    Landfills (edgar hyde annual)
C    K=4    fossil fuel prod, transmis, transform (edgar hyde annual)
C    K=5    termites (from Jean annual)
C    K=6    Soil Sink (from Jean annual)
C    K=7    ocean (from Jean annual)
C    K=8    fresh water lakes (from Jean annual)
C    K=9    Misc. Ground Source (from Jean annual)
C    K=10   wetlands and tundra (from Jean monthly)
C    These 5 are from edgar hyde, but with our monthly variation from
C    biomass burning apllied to the annual values. In the case of
C    agricultural land use, our rice cultivation variation was applied:
C    K=11   agricultural waste burning
C    K=12   deforestation
C    K=13   savannah burning
C    K=14   biofuel production, transformation, combustion
C    K=15   agricultural land activities     
!@var nanns,nmons: number of annual and monthly input files
      integer, parameter :: nanns=9,nmons=6
      integer ann_units(nanns-3),mon_units(nmons)
      character*15 :: ann_files(nanns-3) =
     & (/'CH4_ANIMALS_EH ','CH4_F_FUEL_C_EH','CH4_LANDFILL_EH',
     &   'CH4_F_FUEL_P_EH','CH4_TERMITES   ','CH4_SOIL_ABS   '/)
      logical :: ann_bins(nanns-3)=(/.true.,.true.,.true.,.true.,
     &        .true.,.true./)
      character*9 :: mon_files(nmons) =
     &   (/'CH4_WETL ','CH4_AW_EH','CH4_DF_EH','CH4_SB_EH',
     &     'CH4_BF_EH','CH4_AL_EH'/)
      real*8 adj_wet(jm)
      data adj_wet/15*0.6585,16*1.761,15*0.6585/ !zonal adj FOR WETLANDS
      integer :: kwet=10           !!! position of wetlands array in src
      logical :: mon_bins(nmons)=(/.true.,.true.,.true.,
     &        .true.,.true.,.true./)
      real*8 tlca(im,jm,nmons),tlcb(im,jm,nmons)  ! for monthly sources
      real*8 frac
      integer i,j,jj,nt,iact,iu,k,imon(nmons)
      integer :: jdlast=0
      save ifirst,jdlast,tlca,tlcb,mon_units,imon

      byDTsrc = 1.d0/dtsrc
      if (itime.lt.itime_tr0(nt)) return
      if(PI_run.eq.1)
     &call stop_model('edgar-hyde PI run not set up yet',255)
C****
C**** Annual Sources and sinks
C**** Apply adjustment factors to bring sources into balance.
C**** Also, Jean's sources are in KG CH4/M2/YEAR and
C**** the EDGAR-HYDE sources are in KG/4x5grid/HR. Both need to be
C**** converted to (kg CH4)/m2/sec:
C****
      ann_conv(2)= 1.d0/(sday*JDperY) ! Jean's
      if (ifirst) then
        k = 0
        call openunits(ann_files,ann_units,ann_bins,nanns-3)
        do iu = ann_units(1),ann_units(nanns-3)
          k = k+1
          call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          do j=1,jm
            ann_conv(1)=bydxyp(j)*byDTsrc ! EH
            src(:,j,k) = src(:,j,k)*ann_conv(EHorJEAN(k))
            if(adjust(k)) src(:,j,k) = src(:,j,k)*adj(k)
          end do
        end do
        call closeunits(ann_units,nanns-3)
        ! 3 miscellaneous sources
          k = k+1
          do j=1,jm
            ann_conv(1)=bydxyp(j)*byDTsrc ! EH
            src(:,j,k) = focean(:,j)*ann_conv(EHorJEAN(k))
            if(adjust(k)) src(:,j,k) = src(:,j,k)*adj(k)            
          end do
          k = k+1
          do j=1,jm
            ann_conv(1)=bydxyp(j)*byDTsrc ! EH
            src(:,j,k) =  flake(:,j)*ann_conv(EHorJEAN(k))
            if(adjust(k)) src(:,j,k) = src(:,j,k)*adj(k)
          end do
          k = k+1
          do j=1,jm
            ann_conv(1)=bydxyp(j)*byDTsrc ! EH         
            src(:,j,k) = fearth(:,j)*ann_conv(EHorJEAN(k))
            if(adjust(k)) src(:,j,k) = src(:,j,k)*adj(k)
          end do
        call openunits(mon_files,mon_units,mon_bins,nmons)
      endif
C****
C**** Monthly sources are interpolated to the current day
C**** Apply adjustment factors to bring sources into balance.
C**** Jean's monthly sources are already in correct units.
C**** The EDGAR-HYDE sources are in KG/4x5grid/HR and need to be
C**** converted to (kg CH4)/m2/sec:
C****
      mon_conv(2)= 1.d0 ! Jean's
      ifirst = .false.
      jj = 0
      do k = nanns+1,nsrc
        jj = jj+1
        call read_monthly_sources(mon_units(jj),jdlast,
     *    tlca(1,1,jj),tlcb(1,1,jj),src(1,1,k),frac,imon(jj))
        do j=1,jm
          mon_conv(1)=bydxyp(j)*byDTsrc
          src(:,j,k) = src(:,j,k)*mon_conv(EHorJEAN(k))
          if(adjust(k)) src(:,j,k) = src(:,j,k)*adj(k)
        end do
      end do
      jdlast = jday
      write(6,*) trname(nt),'Sources interpolated to current day',frac
      call sys_flush(6)
C****
C**** Zonal adjustment for combined wetlands and tundra
C****
      do j=1,jm
        src(:,j,kwet) = src(:,j,kwet)*adj_wet(j)
      end do
C****
C**** Also, increase the wetlands + tundra CH4 emissions:
C****
      src(:,:,kwet)=2.2d0*src(:,:,kwet)
C  
      return
      end subroutine read_CH4_sources


      subroutine read_CO_sources(nt,iact)
!@sum reads in CO surface sources and sinks: EDGA-HYDE VERSION
!@auth Jean Lerner/Greg Faluvegi
C****
C**** There is 4 monthly sources and 2 annual source
C**** Annual sources are read in at start and re-start of run only
C**** Monthly sources are interpolated each day
      USE MODEL_COM, only: itime,jday,JDperY,im,jm,DTsrc
      USE GEOM, only: BYDXYP
      USE CONSTANT, only: sday,hrday
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      USE TRACER_COM, only: itime_tr0,trname
      use TRACER_SOURCES, only: src=>CO_src,nsrc=>nCOsrc,correct_CO_ind
      USE TRCHEM_Shindell_COM, only: PI_run,PIratio_indus,PIratio_bburn
      
      implicit none
      character*80 title
      logical :: ifirst=.true.
!@var nanns,nmons: number of annual and monthly input files
      integer, parameter :: nanns=2,nmons=4
      integer ann_units(nanns),mon_units(nmons)
      character*8 :: ann_files(nanns) = (/'CO_FC_EH','CO_IN_EH'/)
      logical :: ann_bins(nanns)=(/.true.,.true./) ! binary file?
      character*10 :: mon_files(nmons) = (/'CO_AW_EH','CO_SB_EH',
     & 'CO_BF_EH','CO_DF_EH'/)
      logical :: mon_bins(nmons)=(/.true.,.true.,.true.,.true./)
      real*8 tlca(im,jm,nmons),tlcb(im,jm,nmons)  ! for monthly sources
      real*8 frac,byDTsrc
      integer i,j,jj,nt,iact,iu,k,imon(nmons)
      integer :: jdlast=0
      save ifirst,jdlast,tlca,tlcb,mon_units,imon
C
C    K=1    Fossil fuel combustion           (edgar hyde annual)
C    K=2    industrial prod & cons processes (edgar hyde annual)
C    These 4 are from edgar hyde, but with our monthly variation from
C    biomass burning apllied to the annual values:
C    K=3    agricultural waste burning
C    K=4    savannah burning
C    K=5    biofuel production, transformation, combustion
C    K=6    deforestation
C
      if (itime.lt.itime_tr0(nt)) return
      byDTsrc = 1.d0/dtsrc
      if(PI_run.eq.1)
     &call stop_model('edgar-hyde PI run not set up yet',255)
C****
C**** Annual Sources and sink
C**** The EDGAR-HYDE sources are in KG/4x5grid/HR and need to be
C**** converted to KG/m2/sec:
C****
      if (ifirst) then
        call openunits(ann_files,ann_units,ann_bins,nanns)
        k = 0
        do iu = ann_units(1),ann_units(nanns)
          k = k+1
          call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          do j=1,jm
            src(:,j,k) = src(:,j,k)*bydxyp(j)*byDTsrc
          end do
        end do
        call closeunits(ann_units,nanns)
        call openunits(mon_files,mon_units,mon_bins,nmons)
      endif
C****
C**** Monthly sources are interpolated to the current day
C**** The EDGAR-HYDE sources are in KG/4x5grid/HR and need to be
C**** converted to KG/m2/sec:
C****
      ifirst = .false.
      jj = 0
      do k=nanns+1,nsrc
        jj = jj+1
        call read_monthly_sources(mon_units(jj),jdlast,
     *    tlca(1,1,jj),tlcb(1,1,jj),src(1,1,k),frac,imon(jj))
        do j=1,jm
          src(:,j,k) = src(:,j,k)*bydxyp(j)*byDTsrc
        end do
      end do
      jdlast = jday
      write(6,*) trname(nt),'Sources interpolated to current day',frac
      call sys_flush(6)
C****
      return
      end subroutine read_CO_sources


      subroutine read_Alkenes_sources(nt,iact)
!@sum reads in Alkenes surface sources and sinks: EDGA-HYDE VERSION
!@auth Jean Lerner/Greg Faluvegi
C****
C**** There are 5 monthly sources and 3 annual source
C**** Annual sources are read in at start and re-start of run only
C**** Monthly sources are interpolated each day
      USE MODEL_COM, only: itime,jday,JDperY,im,jm,DTsrc
      USE CONSTANT, only: sday,hrday
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      USE GEOM, only: BYDXYP
      USE TRACER_COM, only: itime_tr0,trname
      use TRACER_SOURCES, only: src=>Alkenes_src,nsrc=>nAlkenessrc
      USE TRCHEM_Shindell_COM, only: PI_run,PIratio_indus,PIratio_bburn

      implicit none
      character*80 title
      logical :: ifirst=.true.
!@var nanns,nmons: number of annual and monthly input files
      integer, parameter :: nanns=3,nmons=5
      integer ann_units(nanns),mon_units(nmons)
      character*18 :: ann_files(nanns) = (/'Alkenes_FC_EH     ',
     &'Alkenes_IN_EH     ','Alkenes_FP_EH     '/)
      logical :: ann_bins(nanns)=(/.true.,.true.,.true./) ! bin files?
      character*18 :: mon_files(nmons) =      (/'Alkenes_VEGETATION',
     &'Alkenes_AW_EH     ','Alkenes_SB_EH     ','Alkenes_BF_EH     ',
     &'Alkenes_DF_EH     '/)
      logical :: mon_bins(nmons)=(/.true.,.true.,.true.,.true.,.true./)
      real*8 tlca(im,jm,nmons),tlcb(im,jm,nmons)  ! for monthly sources
      real*8 frac,byDTsrc
      integer i,j,jj,nt,iact,iu,k,imon(nmons)
      integer :: jdlast=0
      save ifirst,jdlast,tlca,tlcb,mon_units,imon
C
C    K=1    Fossil fuel combustion           (edgar hyde annual)
C    K=2    industrial prod & cons processes (edgar hyde annual)
C    K=3    fossil fuel prod, transm, transfor      (e-h annual)
C    K=4    Vegetation emissions (from B436000M23 model, monthly)
C    These 4 are from edgar hyde, but with our monthly variation from
C    biomass burning apllied to the annual values:
C    K=5    agricultural waste burning
C    K=6    savannah burning
C    K=7    biofuel production, transformation, combustion
C    K=8    deforestation
      if (itime.lt.itime_tr0(nt)) return
      byDTsrc = 1.d0/dtsrc
      if(PI_run.eq.1)
     &call stop_model('edgar-hyde PI run not set up yet',255)
C****
C**** Annual Sources
C**** The EDGAR-HYDE sources are in KG C/4x5grid/HR and need to be
C**** converted to KG C/m2/sec:
C****
      if (ifirst) then
        call openunits(ann_files,ann_units,ann_bins,nanns)
        k = 0
        do iu = ann_units(1),ann_units(nanns)
          k = k+1
          call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          do j = 1,jm
            src(:,j,k) = src(:,j,k)*BYDXYP(j)*byDTsrc
          end do
        end do
        call closeunits(ann_units,nanns)
        call openunits(mon_files,mon_units,mon_bins,nmons)
      endif
C****
C**** Monthly sources are interpolated to the current day.
C**** The EDGAR-HYDE sources are in KG C/4x5grid/HR and need to be
C**** converted to KG C/m2/sec. Same is true for the B436000M23
C**** vegetation source:
C****
      ifirst = .false.
      jj = 0
      do k=nanns+1,nsrc
        jj = jj+1
        call read_monthly_sources(mon_units(jj),jdlast,
     *    tlca(1,1,jj),tlcb(1,1,jj),src(1,1,k),frac,imon(jj))
          do j = 1,jm
            src(:,j,k) = src(:,j,k)*BYDXYP(j)*byDTsrc
          end do
      end do
      jdlast = jday
      write(6,*) trname(nt),'Sources interpolated to current day',frac
      call sys_flush(6)
C****
      return
      end subroutine read_Alkenes_sources


      subroutine read_Paraffin_sources(nt,iact)
!@sum reads in Paraffin surface sources and sinks: EDGA-HYDE VERSION
!@auth Jean Lerner/Greg Faluvegi
C****
C**** There are 5 monthly sources and 3 annual source
C**** Annual sources are read in at start and re-start of run only
C**** Monthly sources are interpolated each day
      USE MODEL_COM, only: itime,jday,JDperY,im,jm,DTsrc
      USE CONSTANT, only: sday,hrday
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      USE GEOM, only: BYDXYP
      USE TRACER_COM, only: itime_tr0,trname
      use TRACER_SOURCES, only: src=>Paraffin_src,nsrc=>nParaffinsrc
      USE TRCHEM_Shindell_COM, only: PI_run,PIratio_indus,PIratio_bburn
      
      implicit none
      character*80 title
      logical :: ifirst=.true.
!@var nanns,nmons: number of annual and monthly input files
      integer, parameter :: nanns=3,nmons=5
      integer ann_units(nanns),mon_units(nmons)
      character*19 :: ann_files(nanns) = (/'Paraffin_FC_EH     ',
     &'Paraffin_IN_EH     ','Paraffin_FP_EH     '/)
      logical :: ann_bins(nanns)=(/.true.,.true.,.true./) ! bin files?
      character*19 :: mon_files(nmons) =        (/'Paraffin_VEGETATION',
     &'Paraffin_AW_EH     ','Paraffin_SB_EH     ','Paraffin_BF_EH     ',
     &'Paraffin_DF_EH     '/)
      logical :: mon_bins(nmons)=(/.true.,.true.,.true.,.true.,.true./)
      real*8 tlca(im,jm,nmons),tlcb(im,jm,nmons)  ! for monthly sources
      real*8 frac,byDTsrc
      integer i,j,jj,nt,iact,iu,k,imon(nmons)
      integer :: jdlast=0
      save ifirst,jdlast,tlca,tlcb,mon_units,imon
C
C    K=1    Fossil fuel combustion           (edgar hyde annual)
C    K=2    industrial prod & cons processes (edgar hyde annual)
C    K=3    fossil fuel prod, transm, transfor      (e-h annual)
C    K=4    Vegetation emissions (from B436000M23 model, monthly)
C    These 4 are from edgar hyde, but with our monthly variation from
C    biomass burning apllied to the annual values:
C    K=5    agricultural waste burning
C    K=6    savannah burning
C    K=7    biofuel production, transformation, combustion
C    K=8    deforestation
C
      if (itime.lt.itime_tr0(nt)) return
      byDTsrc = 1.d0/dtsrc
      if(PI_run.eq.1)
     &call stop_model('edgar-hyde PI run not set up yet',255)
C****
C**** Annual Sources
C**** The EDGAR-HYDE sources are in KG C/4x5grid/HR and need to be
C**** converted to KG C/m2/sec:
C****
      if (ifirst) then
        call openunits(ann_files,ann_units,ann_bins,nanns)
        k = 0
        do iu = ann_units(1),ann_units(nanns)
          k = k+1
          call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          do j = 1,jm
            src(:,j,k) = src(:,j,k)*BYDXYP(j)*byDTsrc
          end do
        end do
        call closeunits(ann_units,nanns)
        call openunits(mon_files,mon_units,mon_bins,nmons)
      endif
C****
C**** Monthly sources are interpolated to the current day.
C**** The EDGAR-HYDE sources are in KG C/4x5grid/HR and need to be
C**** converted to KG C/m2/sec. Same is true for the B436000M23
C**** vegetation source:
C****
      ifirst = .false.
      jj = 0
      do k=nanns+1,nsrc
        jj = jj+1
        call read_monthly_sources(mon_units(jj),jdlast,
     *  tlca(1,1,jj),tlcb(1,1,jj),src(1,1,k),frac,imon(jj))
        do j = 1,jm
          src(:,j,k) = src(:,j,k)*BYDXYP(j)*byDTsrc
        end do
      end do
      jdlast = jday
      write(6,*) trname(nt),'Sources interpolated to current day',frac
      call sys_flush(6)
C****
      return
      end subroutine read_Paraffin_sources


      subroutine read_NOx_sources(nt,iact)
!@sum reads in NOx surface sources and sinks: EDGA-HYDE VERSION
!@auth Jean Lerner/Greg Faluvegi
C****
C**** There are 5 monthly sources and 2 annual source
C**** Annual sources are read in at start and re-start of run only
C**** Monthly sources are interpolated each day
      USE MODEL_COM, only: itime,jday,JDperY,im,jm,DTsrc
      USE GEOM, only: BYDXYP
      USE CONSTANT, only: sday,hrday
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      USE TRACER_COM, only: itime_tr0,trname
      use TRACER_SOURCES, only: src=>NOx_src,nsrc=>nNOxsrc
      USE TRCHEM_Shindell_COM, only: PI_run,PIratio_indus,PIratio_bburn

      implicit none
      character*80 title
      logical :: ifirst=.true.
!@var nanns,nmons: number of annual and monthly input files
      integer, parameter :: nanns=2,nmons=5
      integer ann_units(nanns),mon_units(nmons)
      character*9 :: ann_files(nanns) = (/'NOx_FC_EH','NOx_IN_EH'/)
      logical :: ann_bins(nanns)=(/.true.,.true./) ! binary files?
      character*9 :: mon_files(nmons) = (/'NOx_AL_EH','NOx_AW_EH',
     &    'NOx_BF_EH','NOx_DF_EH','NOx_SB_EH'/)
      logical :: mon_bins(nmons)=(/.true.,.true.,.true.,.true.,.true./)
      real*8 tlca(im,jm,nmons),tlcb(im,jm,nmons)  ! for monthly sources
      real*8 frac,byDTsrc
      integer i,j,jj,nt,iact,iu,k,imon(nmons)
      integer :: jdlast=0
      save ifirst,jdlast,tlca,tlcb,mon_units,imon

C------ annual -------
C    K=1    Fossil Fuel Combustion (edgar hyde annual)
C    K=2    industrial prod & cons processes (edgar hyde annual)
C------ monthly ------
C    K=3    agricultural land activities (E-H variation as soil NOx)
C    K=4    agricultural waste burning (edgar hyde variation bio.burn)
C    K=5    biofuel production, transformation, combustion (")
C    K=6    deforestation (edgar hyde variation as bio. burn)
C    K=7    savannah burning (edgar hyde variation as bio. burn)

      if (itime.lt.itime_tr0(nt)) return
      byDTsrc = 1.d0/dtsrc
      if(PI_run.eq.1)
     &call stop_model('edgar-hyde PI run not set up yet',255)
C****
C**** Annual Sources
C**** The EDGAR-HYDE sources are in KG N/4x5grid/HR and need to be
C**** converted to KG N/m2/sec:
C****
      if (ifirst) then
        call openunits(ann_files,ann_units,ann_bins,nanns)
        k = 0
        do iu = ann_units(1),ann_units(nanns)
          k = k+1
          call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          do j=1,jm
            src(:,j,k) = src(:,j,k)*BYDXYP(j)*byDTsrc
          end do
        end do
        call closeunits(ann_units,nanns)
        call openunits(mon_files,mon_units,mon_bins,nmons)
      endif
C****
C**** Monthly sources are interpolated to the current day
C**** The EDGAR-HYDE sources are in KG N/4x5grid/HR and need to be
C**** converted to KG N/m2/sec:
C****
      ifirst = .false.
      jj = 0
      do k=nanns+1,nsrc
        jj = jj+1
        call read_monthly_sources(mon_units(jj),jdlast,
     *  tlca(1,1,jj),tlcb(1,1,jj),src(1,1,k),frac,imon(jj))
        do j=1,jm
          src(:,j,k) = src(:,j,k)*BYDXYP(j)*byDTsrc
        end do 
      end do
      jdlast = jday
      write(6,*) trname(nt),'Sources interpolated to current day',frac
      call sys_flush(6)
C****
      return
      end subroutine read_NOx_sources   
!
#else
!
!end of edgar-hyde sources. normal model sources follow:
!  
      subroutine read_CH4_sources(nt,iact)
!@sum reads in CH4 surface sources and sinks
!@auth Jean Lerner
C****
C**** There are 3 monthly sources and 11 annual sources
C**** Annual sources are read in at start and re-start of run only
C**** Monthly sources are interpolated each day
      USE MODEL_COM, only: itime,JDperY,im,jm,jday,focean,fearth
      USE LAKES_COM, only: flake
      USE CONSTANT, only: sday
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      USE TRACER_COM, only: itime_tr0,trname
      use TRACER_SOURCES, only: src=>ch4_src,nsrc=>nch4src
#ifdef INTERACTIVE_WETLANDS_CH4      
     &       ,nday_ncep,by_nday_ncep,first_ncep,iday_ncep,day_ncep
     &       ,i0_ncep,avg_ncep,sum_ncep,fact_ncep,avg_modPT
#endif
      USE TRCHEM_Shindell_COM, only: PI_run,PIratio_indus,PIratio_bburn
      
      implicit none
      character*80 title
!@var adj Factors that tune the total amount of individual sources
!@var PIfact factor for altering source to preindustrial (or not)
      real*8 PIfact
      logical :: ifirst=.true.
      real*8 adj(nsrc)
C  O  data adj/1.3847,1.0285,3.904,1.659,1.233,1.194,0.999,
C! R *  7.2154, 3.7247d-5,3.1399d-4,5.4838d-5,   !Model II prime
C  I *  7.2154, 3.5997d-5,17.330d-4,5.3558d-5,
C  G *  0.4369,0.7533,0.9818/
      data adj/1.0787,0.6891,2.6157,1.659,0.617,1.588,0.999,
     *  7.2154, 3.5997d-5,17.330d-4,5.3558d-5,
     *  0.3277,0.3767,0.9818/
!@var nanns,nmons: number of annual and monthly input files
      integer, parameter :: nanns=11,nmons=3
      integer ann_units(nanns-3),mon_units(nmons)
      character*12 :: ann_files(nanns-3) =
     *  (/'CH4_ANIMALS ','CH4_COALMINE','CH4_GASLEAK ','CH4_GASVENT ',
     *    'CH4_CITYDUMP','CH4_SOIL_ABS','CH4_TERMITES','CH4_COALBURN'/)
      logical :: ann_bins(nanns-3)=(/.true.,.true.,.true.,.true.,
     *        .true.,.true.,.true.,.true./)
      character*8 :: mon_files(nmons) =
     *   (/'CH4_BURN','CH4_RICE','CH4_WETL'/)
      real*8 adj_wet(jm)
      data adj_wet/15*0.6585,16*1.761,15*0.6585/ !zonal adj FOR WETLANDS
      integer :: kwet=14           !!! position of wetlands array in src
      logical :: mon_bins(nmons)=(/.true.,.true.,.true./)
      real*8 tlca(im,jm,nmons),tlcb(im,jm,nmons)  ! for monthly sources
      real*8 frac
      integer i,j,nt,iact,iu,k,imon(nmons)
      integer :: jdlast=0
#ifdef INTERACTIVE_WETLANDS_CH4
      integer n,m
#endif
      save ifirst,jdlast,tlca,tlcb,mon_units,imon

      if (itime.lt.itime_tr0(nt)) return
#ifdef INTERACTIVE_WETLANDS_CH4
      by_nday_ncep(:)=1.d0/real(nday_ncep(:))
#endif
C****
C**** Annual Sources and sinks
C**** Apply adjustment factors to bring sources into balance
C**** Annual sources are in KG C/M2/Y
C**** Sources need to be kg/m^2 s; convert /year to /s
C****
      if (ifirst) then
        k = 0
        call openunits(ann_files,ann_units,ann_bins,nanns-3)
        do iu = ann_units(1),ann_units(nanns-3)
          k = k+1
          PIfact=1.d0
          select case(PI_run)
          case(1) ! pre-industrial
            select case(k)
            case(2,3,4,5,8); PIfact=PIratio_indus
            end select
          end select
          call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          src(:,:,k) = src(:,:,k)*PIfact*adj(k)/(sday*JDperY)
        end do
        call closeunits(ann_units,nanns-3)
        ! 3 miscellaneous sources
          k = k+1
          src(:,:,k) = focean(:,:)*adj(k)/(sday*JDperY)
          k = k+1
          src(:,:,k) =  flake(:,:)*adj(k)/(sday*JDperY)
          k = k+1
          src(:,:,k) = fearth(:,:)*adj(k)/(sday*JDperY)
        call openunits(mon_files,mon_units,mon_bins,nmons)
      endif
C****
C**** Monthly sources are interpolated to the current day
C****
C**** Also, Apply adjustment factors to bring sources into balance
C**** Monthly sources are in KG C/M2/S => src in kg/m^2 s
C****
      ifirst = .false.
      j = 0
      do k = nanns+1,nsrc
        PIfact=1.d0
        select case(PI_run)
        case(1) ! pre-industrial
          select case(k)
          case(nanns+1); PIfact=PIratio_bburn
          end select
        end select
        j = j+1
        call read_monthly_sources(mon_units(j),jdlast,
     *    tlca(1,1,j),tlcb(1,1,j),src(1,1,k),frac,imon(j))
        src(:,:,k) = src(:,:,k)*adj(k)*PIfact
      end do
      jdlast = jday
      write(6,*) trname(nt),'Sources interpolated to current day',frac
      call sys_flush(6)
C****
C**** Zonal adjustment for combined wetlands and tundra
C****
      do j=1,jm
        src(:,j,kwet) = src(:,j,kwet)*adj_wet(j)
      end do
C****
C**** Also, increase the wetlands + tundra CH4 emissions:
C****
      src(:,:,kwet)=2.2d0*src(:,:,kwet)
C  
#ifdef INTERACTIVE_WETLANDS_CH4
C****
C**** Adjust the wetlands+tundra CH4 source based on 1st layer ground
C**** temperature from one week ago, precipitation from 2 weeks ago,
C**** and B. Walter's regression coefficients.  I.e.:
C**** 
C**** CH4emis = CH4emis + (alpha*TempAnom + beta*PrecAnom)
C****
C
C NCEP Precipitation(m=1) and Temperature(m=2) running averages:
C
      do m=1,2
        if(m.gt.2)call stop_model('check on ncep index m',255)
        if(first_ncep(m))then !accumulate first nday_ncep(m) days
          iday_ncep(m) = iday_ncep(m) + 1
          day_ncep(:,:,iday_ncep,m)=PTBA(:,:,m)
          if(iday_ncep(m).eq.nday_ncep(m))then !end of averaging period
            sum_ncep(:,:,m)=0.d0
            do n=1,nday_ncep(m)
              sum_ncep(:,:,m) = sum_ncep(:,:,m) + day_ncep(:,:,n,m)
            end do
            first_ncep(m)=.false.
            iday_ncep(m) = 0
            i0_ncep(m)   = 0
            avg_ncep(:,:,m) = sum_ncep(:,:,m) * by_nday_ncep(m)
          end if
        else                     ! no longer first averaging period
          i0_ncep(m) = i0_ncep(m) + 1
          if(i0_ncep(m).gt.nday_ncep(m)) i0_ncep(m) = 1
C         UPDATE RUNNING AVERAGE:
          day_ncep_tmp(:,:,m) = PTBA(:,:,m)
          sum_ncep(:,:,m)=sum_ncep(:,:,m)-day_ncep(:,:,i0_ncep(m),m)
          day_ncep(:,:,i0_ncep(m),m) = day_ncep_tmp(:,:,m)
          sum_ncep(:,:,m) = sum_ncep(:,:,m) + day_ncep(:,:,i0_ncep(m),m)
          avg_ncep(:,:,m) = sum_ncep(:,:,m) * by_nday_ncep(m)
        endif
      end do ! m
C
C Don't alter source until enough statistics are built up:
      if(first_ncep(1) .or. first_ncep(2)) RETURN
C Otherwise, apply the adjustments, limit change to be positive:
      do m=1,2
        src(:,:,kwet) = src(:,:,kwet) + PTBA(:,:,m+2) * 
     &            (fact_ncep(m)*avg_modPT(:,:,m) - avg_ncep(:,:,m))
      end do
      src(:,:,kwet)=max(src(:,j,kwet),0.d0)   
#endif
      return
      end subroutine read_CH4_sources


      subroutine read_CO_sources(nt,iact)
!@sum reads in CO surface sources and sinks
!@auth Jean Lerner/Greg Faluvegi
C****
C**** There is 1 monthly sources and 1 annual source
C**** Annual sources are read in at start and re-start of run only
C**** Monthly sources are interpolated each day
      USE MODEL_COM, only: itime,jday,JDperY,im,jm
      USE CONSTANT, only: sday,hrday
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      USE TRACER_COM, only: itime_tr0,trname
      use TRACER_SOURCES, only: src=>CO_src,nsrc=>nCOsrc,correct_CO_ind
      USE TRCHEM_Shindell_COM, only: PI_run,PIratio_indus,PIratio_bburn
      implicit none
      character*80 title
      logical :: ifirst=.true.
!@var PIfact factor for altering source to preindustrial (or not)
      real*8 PIfact
!@var nanns,nmons: number of annual and monthly input files
      integer, parameter :: nanns=1,nmons=1
      integer ann_units(nanns),mon_units(nmons)
      character*13 :: ann_files(nanns) = (/'CO_INDUSTRIAL'/)
      logical :: ann_bins(nanns)=(/.true./) ! binary file?
      character*10 :: mon_files(nmons) = (/'CO_BIOMASS'/)
      logical :: mon_bins(nmons)=(/.true./) ! binary file?
      real*8 tlca(im,jm,nmons),tlcb(im,jm,nmons)  ! for monthly sources
      real*8 frac
      integer i,j,nt,iact,iu,k,imon(nmons)
      integer :: jdlast=0
      save ifirst,jdlast,tlca,tlcb,mon_units,imon

      if (itime.lt.itime_tr0(nt)) return
C****
C**** Annual Sources and sink
C**** I believe this source is already in (kg CO)/m2/sec, so no
C**** conversion needed.  However, Lee put in a correction factor of
C**** 520/410 = correct_CO_ind, so put that in here:
C****
      if (ifirst) then
        call openunits(ann_files,ann_units,ann_bins,nanns)
        k = 0
        do iu = ann_units(1),ann_units(nanns)
          k = k+1
          select case(PI_run)
          case(1);      PIfact=PIratio_indus
          case default; PIfact=1.d0
          end select
          call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          src(:,:,k) = src(:,:,k)*correct_CO_ind*PIfact
        end do
        call closeunits(ann_units,nanns)
        call openunits(mon_files,mon_units,mon_bins,nmons)
      endif
C****
C**** Monthly sources are interpolated to the current day
C**** I believe this source is already in (kg CO)/m2/sec, so no
C**** conversion needed:
C****
      ifirst = .false.
      j = 0
      do k=nanns+1,nsrc
        PIfact=1.d0
        select case(PI_run)
        case(1) ! pre-industrial
          select case(k)
          case(nanns+1); PIfact=PIratio_bburn
          end select
        end select
        j = j+1
        call read_monthly_sources(mon_units(j),jdlast,
     *    tlca(1,1,j),tlcb(1,1,j),src(1,1,k),frac,imon(j))
        src(:,:,k) = src(:,:,k)*PIfact
      end do
      jdlast = jday
      write(6,*) trname(nt),'Sources interpolated to current day',frac
      call sys_flush(6)
C****
      return
      end subroutine read_CO_sources


      subroutine read_Alkenes_sources(nt,iact)
!@sum reads in Alkenes surface sources and sinks
!@auth Jean Lerner/Greg Faluvegi
C****
C**** There are 2 monthly sources and 1 annual source
C**** Annual sources are read in at start and re-start of run only
C**** Monthly sources are interpolated each day
      USE MODEL_COM, only: itime,jday,JDperY,im,jm,DTsrc
      USE CONSTANT, only: sday,hrday
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      USE GEOM, only: BYDXYP
      USE TRACER_COM, only: itime_tr0,trname
      use TRACER_SOURCES, only: src=>Alkenes_src,nsrc=>nAlkenessrc
      USE TRCHEM_Shindell_COM, only: PI_run,PIratio_indus,PIratio_bburn
   
      implicit none
      character*80 title
      logical :: ifirst=.true.
!@var nanns,nmons: number of annual and monthly input files
      integer, parameter :: nanns=1,nmons=2
      integer ann_units(nanns),mon_units(nmons)
      character*18 :: ann_files(nanns) = (/'Alkenes_INDUSTRIAL'/)
      logical :: ann_bins(nanns)=(/.true./) ! binary files?
      character*18 :: mon_files(nmons) = (/'Alkenes_BIOMASS   ',
     &                                     'Alkenes_VEGETATION'/)
      logical :: mon_bins(nmons)=(/.true.,.true./) ! binary files?
      real*8 tlca(im,jm,nmons),tlcb(im,jm,nmons)  ! for monthly sources
      real*8 frac,byDTsrc
      integer i,j,jj,nt,iact,iu,k,imon(nmons)
      integer :: jdlast=0
!@var PIfact factor for altering source to preindustrial (or not)
      real*8 PIfact
      save ifirst,jdlast,tlca,tlcb,mon_units,imon

      if (itime.lt.itime_tr0(nt)) return
      byDTsrc = 1.d0/dtsrc
C****
C**** Annual Sources and sink
C**** I believe input file is in (kg C)/4x5 grid/hr. So,
C**** convert to (kg C)/m2/s.
C****
      if (ifirst) then
        call openunits(ann_files,ann_units,ann_bins,nanns)
        k = 0
        do iu = ann_units(1),ann_units(nanns)
          k = k+1
          select case(PI_run)
          case(1);      PIfact=PIratio_indus
          case default; PIfact=1.d0
          end select
          call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          do j = 1,jm
            src(:,j,k) = src(:,j,k)*BYDXYP(j)*byDTsrc*PIfact
          end do
        end do
        call closeunits(ann_units,nanns)
        call openunits(mon_files,mon_units,mon_bins,nmons)
      endif
C****
C**** Monthly sources are interpolated to the current day
C**** I believe input files are in (kg C)/4x5 grid/hr. So,
C**** convert to (kg C)/m2/s.
C****
      ifirst = .false.
      j = 0
      do k=nanns+1,nsrc
        PIfact=1.d0
        select case(PI_run)
        case(1) ! pre-industrial
          select case(k)
          case(nanns+1); PIfact=PIratio_bburn
          end select
        end select
        j = j+1
        call read_monthly_sources(mon_units(j),jdlast,
     *    tlca(1,1,j),tlcb(1,1,j),src(1,1,k),frac,imon(j))
          do jj = 1,jm
           src(:,jj,k)=src(:,jj,k)*BYDXYP(jj)*byDTsrc*PIfact
          end do
      end do
      jdlast = jday
      write(6,*) trname(nt),'Sources interpolated to current day',frac
      call sys_flush(6)
C****
      return
      end subroutine read_Alkenes_sources


      subroutine read_Paraffin_sources(nt,iact)
!@sum reads in Paraffin surface sources and sinks
!@auth Jean Lerner/Greg Faluvegi
C****
C**** There are 2 monthly sources and 1 annual source
C**** Annual sources are read in at start and re-start of run only
C**** Monthly sources are interpolated each day
      USE MODEL_COM, only: itime,jday,JDperY,im,jm,DTsrc
      USE CONSTANT, only: sday,hrday
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      USE GEOM, only: BYDXYP
      USE TRACER_COM, only: itime_tr0,trname
      use TRACER_SOURCES, only: src=>Paraffin_src,nsrc=>nParaffinsrc
      USE TRCHEM_Shindell_COM, only: PI_run,PIratio_indus,PIratio_bburn
     
      implicit none
      character*80 title
      logical :: ifirst=.true.
!@var nanns,nmons: number of annual and monthly input files
      integer, parameter :: nanns=1,nmons=2
      integer ann_units(nanns),mon_units(nmons)
      character*19 :: ann_files(nanns) = (/'Paraffin_INDUSTRIAL'/)
      logical :: ann_bins(nanns)=(/.true./) ! binary files?
      character*19 :: mon_files(nmons) = (/'Paraffin_BIOMASS   ',
     &                                     'Paraffin_VEGETATION'/)
      logical :: mon_bins(nmons)=(/.true.,.true./) ! binary files?
      real*8 tlca(im,jm,nmons),tlcb(im,jm,nmons)  ! for monthly sources
      real*8 frac,byDTsrc
      integer i,j,jj,nt,iact,iu,k,imon(nmons)
      integer :: jdlast=0
!@var PIfact factor for altering source to preindustrial (or not)
      real*8 PIfact
      save ifirst,jdlast,tlca,tlcb,mon_units,imon

      if (itime.lt.itime_tr0(nt)) return
      byDTsrc = 1.d0/dtsrc
C****
C**** Annual Sources and sink
C**** I believe input file is in (kg C)/4x5 grid/hr. So,
C**** convert to (kg C)/m2/s.
C****
      if (ifirst) then
        call openunits(ann_files,ann_units,ann_bins,nanns)
        k = 0
        do iu = ann_units(1),ann_units(nanns)
          k = k+1
          select case(PI_run)
          case(1);      PIfact=PIratio_indus
          case default; PIfact=1.d0
          end select
          call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          do j = 1,jm
            src(:,j,k) = src(:,j,k)*BYDXYP(j)*byDTsrc*PIfact
          end do
        end do
        call closeunits(ann_units,nanns)
        call openunits(mon_files,mon_units,mon_bins,nmons)
      endif
C****
C**** Monthly sources are interpolated to the current day
C**** I believe input files are in (kg C)/4x5 grid/hr. So,
C**** convert to (kg C)/m2/s.
C****
      ifirst = .false.
      j = 0
      do k=nanns+1,nsrc
        PIfact=1.d0
        select case(PI_run)
        case(1) ! pre-industrial
          select case(k)
          case(nanns+1); PIfact=PIratio_bburn
          end select
        end select
        j = j+1
        call read_monthly_sources(mon_units(j),jdlast,
     *  tlca(1,1,j),tlcb(1,1,j),src(1,1,k),frac,imon(j))
        do jj = 1,jm
          src(:,jj,k) = src(:,jj,k)*BYDXYP(jj)*byDTsrc*PIfact
        end do
      end do
      jdlast = jday
      write(6,*) trname(nt),'Sources interpolated to current day',frac
      call sys_flush(6)
C****
      return
      end subroutine read_Paraffin_sources


      subroutine read_NOx_sources(nt,iact)
!@sum reads in NOx surface sources and sinks
!@auth Jean Lerner/Greg Faluvegi
C****
C**** There are 2 monthly sources and 1 annual source
C**** Annual sources are read in at start and re-start of run only
C**** Monthly sources are interpolated each day
      USE MODEL_COM, only: itime,jday,JDperY,im,jm
      USE CONSTANT, only: sday,hrday
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      USE TRACER_COM, only: itime_tr0,trname
      use TRACER_SOURCES, only: src=>NOx_src,nsrc=>nNOxsrc
      USE TRCHEM_Shindell_COM, only: PI_run,PIratio_indus,PIratio_bburn

      implicit none
      character*80 title
      logical :: ifirst=.true.
!@var nanns,nmons: number of annual and monthly input files
      integer, parameter :: nanns=1,nmons=2
      integer ann_units(nanns),mon_units(nmons)
      character*15 :: ann_files(nanns) = (/'NOx_FOSSIL_FUEL'/)
      logical :: ann_bins(nanns)=(/.true./) ! binary files?
      character*11 :: mon_files(nmons) = (/'NOx_BIOMASS',
     &                                     'NOx_SOIL   '/)
      logical :: mon_bins(nmons)=(/.true.,.true./) ! binary files?
      real*8 tlca(im,jm,nmons),tlcb(im,jm,nmons)  ! for monthly sources
      real*8 frac
      integer i,j,nt,iact,iu,k,imon(nmons)
      integer :: jdlast=0
!@var PIfact factor for altering source to preindustrial (or not)
      real*8 PIfact
      save ifirst,jdlast,tlca,tlcb,mon_units,imon

      if (itime.lt.itime_tr0(nt)) return
C****
C**** Annual Sources and sink
C**** The titles of the input file say this is in KG/Y/m2, so convert
C**** yr-1 to s-1:
C****
      if (ifirst) then
        call openunits(ann_files,ann_units,ann_bins,nanns)
        k = 0
        do iu = ann_units(1),ann_units(nanns)
          k = k+1
          select case(PI_run)
          case(1);      PIfact=PIratio_indus
          case default; PIfact=1.d0
          end select
          call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          src(:,:,k) = src(:,:,k)*PIfact/(sday*JDperY)
        end do
        call closeunits(ann_units,nanns)
        call openunits(mon_files,mon_units,mon_bins,nmons)
      endif
C****
C**** Monthly sources are interpolated to the current day
C**** The titles of the input files say this is in KG/m2/s, so no
C**** conversion is necessary:
C****
      ifirst = .false.
      j = 0
      do k=nanns+1,nsrc
        PIfact=1.d0
        select case(PI_run)
        case(1) ! pre-industrial
          select case(k)
          case(nanns+1); PIfact=PIratio_bburn
          end select
        end select
        j = j+1
        call read_monthly_sources(mon_units(j),jdlast,
     *    tlca(1,1,j),tlcb(1,1,j),src(1,1,k),frac,imon(j))
        src(:,:,k) = src(:,:,k) * PIfact
      end do
      jdlast = jday
      write(6,*) trname(nt),'Sources interpolated to current day',frac
      call sys_flush(6)
C****
      return
      end subroutine read_NOx_sources   
#endif


      subroutine read_Isoprene_sources(nt,iact)
!@sum reads in Isoprene surface sources and sinks
!@auth Jean Lerner/Greg Faluvegi
C****
C**** There is 1 monthly source
C**** Monthly sources are interpolated each day
      USE MODEL_COM, only: itime,jday,JDperY,im,jm,DTsrc
      USE CONSTANT, only: sday,hrday
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      USE GEOM, only: BYDXYP
      USE TRACER_COM, only: itime_tr0,trname
      use TRACER_SOURCES, only: src=>Isoprene_src,nsrc=>nIsoprenesrc
     
      implicit none
      character*80 title
      logical :: ifirst=.true.
!@var nanns,nmons: number of annual and monthly input files
      integer, parameter :: nanns=0, nmons=1
      integer mon_units(nmons)
      character*19 :: mon_files(nmons) = (/'Isoprene_VEGETATION'/)
      logical :: mon_bins(nmons)=(/.true./) ! binary file?
      real*8 tlca(im,jm,nmons),tlcb(im,jm,nmons)  ! for monthly sources
      real*8 frac,byDTsrc
      integer i,j,jj,nt,iact,iu,k,imon(nmons)
      integer :: jdlast=0
      save ifirst,jdlast,tlca,tlcb,mon_units,imon

      if (itime.lt.itime_tr0(nt)) return
      byDTsrc = 1.d0/dtsrc
C****
C**** Monthly sources are interpolated to the current day
C****
C**** I believe input files are in (kg C)/4x5 grid/hr. So,
C**** convert to (kg C)/m2/s.
C****
      if(ifirst) call openunits(mon_files,mon_units,mon_bins,nmons)
      ifirst = .false.
      j = 0
      do k=nanns+1,nsrc
        j = j+1
        call read_monthly_sources(mon_units(j),jdlast,
     *    tlca(1,1,j),tlcb(1,1,j),src(1,1,k),frac,imon(j))
        do jj = 1,jm
          src(:,jj,k) = src(:,jj,k)*BYDXYP(jj)*byDTsrc
        end do
      end do
      jdlast = jday
      write(6,*) trname(nt),'Sources interpolated to current day',frac
      call sys_flush(6)
C****
      return
      end subroutine read_Isoprene_sources
c
c
      SUBROUTINE get_lightning_NOx
!@sum  get_lightning_NOx to define the 3D source of NOx from lightning
!@auth Colin Price / Greg Faluvegi
!@ver  1.0 (based on CB436Tds3M23 & DB396Tds3M23)
c
C**** GLOBAL parameters and variables:
      USE FLUXES, only        : tr3Dsource
      USE TRACER_COM, only    : n_NOx,nLightning
      USE LIGHTNING, only     : HGT,JS,JN,SRCLIGHT,RNOx_lgt
      USE CONSTANT, only      : bygrav
      USE MODEL_COM, only     : fland,IM,JM,LS1,LM
      USE DYNAMICS, only      : LTROPO, PHI
      USE GEOM, only          : BYDXYP
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
c
!@var LATINDX Pickering latitude index
!@var LANDINDX Pickering surface type index
!@var IH Pickering altitude index
!@var HEIGHT Pickering HGT variable specific to I,J,L, point
!@var LEVTROP local variable to hold the tropopause level
!@var ALTTROP altitude of tropopause
!@var ALTTOP altitude at the top of ?
!@var I,J,L loop indicies in x,y,z directions
!@param pmin2psec to convert from min-1 to sec-1
      REAL*8, PARAMETER :: pmin2psec = 1.D0/60.D0
      INTEGER LATINDX,LANDINDX,IH,LEVTROP,I,J,L
      REAL*8, DIMENSION(16) :: HEIGHT
      REAL*8 :: ALTTROP,ALTTOP
c
      DO J=1,JM
      DO I=1,IM
C****    Lightning source function altitude dependence:
C****    Determine if latitude is tropical:
         IF(J.GT.JS.AND.J.LE.JN) THEN
           LATINDX=1
         ELSE
           LATINDX=2
         END IF
C****    Determine if over land or ocean:
         IF(fland(I,J).GE.0.5) THEN
           LANDINDX=1
         ELSE
           LANDINDX=2
         END IF
C****    Choose appropriate height file:
         DO IH=1,16
          HEIGHT(IH)=HGT(LATINDX,LANDINDX,IH)*0.01
         ENDDO
C****    Store local tropopause
         LEVTROP=LTROPO(I,J)
         ALTTROP=PHI(I,J,LEVTROP)*BYGRAV*1.D-3
         IF(ALTTROP.EQ.0.) GOTO 1510
C****    Zero source accumulator
         SRCLIGHT = 0. ! 1 to LM levels
C        Determine which GCM level the height belongs to
         DO 1505 IH=1,16
           L=1
 1502      continue
           ALTTOP=
     &     (PHI(I,J,L)+(PHI(I,J,L+1)-PHI(I,J,L))*.5)*BYGRAV*1.D-3
c          RNOx_lgt is gN/min added per grid box (lightning NOx):
           IF(IH.LE.ALTTOP)THEN
             SRCLIGHT(L)=SRCLIGHT(L)+RNOx_lgt(I,J)*HEIGHT(IH)
             goto 1505
           elseif(IH-1.LT.ALTTOP)then
             SRCLIGHT(L)=SRCLIGHT(L)+RNOx_lgt(I,J)*HEIGHT(IH)*
     &       (ALTTOP-IH+1.)
             SRCLIGHT(L+1)=SRCLIGHT(L+1)+RNOx_lgt(I,J)*HEIGHT(IH)*
     &       (IH-ALTTOP)
             goto 1505
           else
             L=L+1
           ENDIF
           goto 1502
 1505    continue
 1510    continue    
C        save tracer 3D source. convert from gN/min to kgN/s :
         DO L=1,LEVTROP
           tr3Dsource(I,J,L,nLightning,n_NOx) = 
     &     SRCLIGHT(L)*pmin2psec*1.d-3
         END DO
         DO L=LEVTROP+1,LM
            tr3Dsource(I,J,LEVTROP,nLightning,n_NOx) = tr3Dsource(I,J
     $       ,LEVTROP,nLightning,n_NOx) + SRCLIGHT(L)*pmin2psec*1.d-3
         END DO
      END DO                    ! I
      END DO ! J
         
C
      END SUBROUTINE get_lightning_NOx
c
c
      SUBROUTINE get_aircraft_NOx
!@sum  get_aircraft_NOx to define the 3D source of NOx from aircraft
!@auth Drew Shindell? / Greg Faluvegi / Jean Learner
!@ver  1.0 (based on DB396Tds3M23)
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only: itime,jday,JDperY,im,jm,lm,ptop,psf,sig
      USE CONSTANT, only: sday,hrday
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      USE FLUXES, only: tr3Dsource
      USE GEOM, only       : dxyp
      USE TRACER_COM, only: itime_tr0,trname,n_NOx,nAircraft
      use TRACER_SOURCES, only: Laircr
C
      IMPLICIT NONE
c
!@var nanns,nmons: number of annual and monthly input files
!@var l,ll dummy loop variable
!@var pres local pressure variable
      character*80 title
      logical LINJECT
      integer, parameter :: nanns=0,nmons=1
      integer ann_units(nanns),mon_units(nmons),imon(nmons)
      character*12 :: mon_files(nmons) = (/'NOx_AIRCRAFT'/)
      logical :: mon_bins(nmons)=(/.true./) ! binary file?
      real*8 tlca(im,jm,Laircr,nmons),tlcb(im,jm,Laircr,nmons)
      real*8 frac
      integer l,i,j,iu,k,ll
      integer :: jdlast=0
      logical :: ifirst = .true.
      save jdlast,tlca,tlcb,mon_units,imon,ifirst
      REAL*8, DIMENSION(IM,JM,Laircr,1)    :: src
      REAL*8, DIMENSION(LM)                :: pres
      REAL*4, PARAMETER, DIMENSION(Laircr) :: PAIRL =
     & (/1013.25,898.74,794.95,701.08,616.40,540.19,471.81,
     &    410.60,355.99,307.42,264.36,226.32,193.30,165.10,
     &    141.01,120.44,102.87,87.866,75.048/)
c
C**** Local parameters and variables and arguments:
c
C**** Aircraft NOx source input is monthly, on 19 levels. Therefore:
C     19x12=228 records. Read it in here and interpolated each day.

      if (itime.lt.itime_tr0(n_NOx)) return

C****
C**** Monthly sources are interpolated to the current day
C**** The titles of the input files say this is in KG/m2/s, so no
C**** conversion is necessary:
C****
      if (ifirst) call openunits(mon_files,mon_units,mon_bins,nmons)
      ifirst=.false.
      j = 0
      do k=nanns+1,1
        j = j+1
        call read_monthly_3Dsources(Laircr,mon_units(j),jdlast,
     *    tlca(1,1,1,j),tlcb(1,1,1,j),src(1,1,1,k),frac,imon(j))
      end do
      jdlast = jday
      write(6,*)trname(n_NOx)
     *     ,'from aircraft interpolated to current day',frac
      call sys_flush(6)
C====
C====   Place aircraft sources onto model levels:
C====
      tr3Dsource(:,:,:,nAircraft,n_NOx) = 0.
      PRES(:)=SIG(:)*(PSF-PTOP)+PTOP
      DO J=1,JM
       DO I=1,IM
        tr3Dsource(i,j,1,nAircraft,n_NOx) = SRC(I,J,1,1)*dxyp(j)
        DO LL=2,Laircr
          LINJECT=.TRUE.
          DO L=1,LM
           IF(PAIRL(LL).GT.PRES(L).AND.LINJECT) THEN
             tr3Dsource(i,j,l,nAircraft,n_NOx) =
     &       tr3Dsource(i,j,l,nAircraft,n_NOx) + SRC(I,J,LL,1)*dxyp(j)
             LINJECT=.FALSE.
           ENDIF
          ENDDO ! L
        ENDDO   ! LL
       END DO   ! I
      END DO    ! J
C
      return
      END SUBROUTINE get_aircraft_NOx
C
C
      SUBROUTINE get_sulfate_N2O5
!@sum  get_sulfate_N2O5 parameterize N2O5 sink on sulfate. sulfate( )
!@+    array is then used in the chemistry to calculate a chemical
!@+    change.  I.e. this source is not applied directly to the tracer
!@+    mass like other 3D sources are.
!@auth Drew Shindell? / Greg Faluvegi
!@ver  1.0 (based on DB396Tds3M23)
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only: jday,im,jm,lm,ptop,psf,sig
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      use TRACER_SOURCES, only: Lsulf
      USE TRCHEM_Shindell_COM, only: sulfate
C
      IMPLICIT NONE

C**** Local parameters and variables and arguments:
c
!@var nanns,nmons: number of annual and monthly input files
!@param Psulf pressure levels of the sulfate SA input file
      real*8, parameter, dimension(Lsulf) :: Psulf = (/
     & 0.9720D+03,0.9445D+03,0.9065D+03,
     & 0.8515D+03,0.7645D+03,0.6400D+03,0.4975D+03,0.3695D+03,
     & 0.2795D+03,0.2185D+03,0.1710D+03,0.1335D+03,0.1016D+03,
     & 0.7120D+02,0.4390D+02,0.2470D+02,0.1390D+02,0.7315D+01,
     & 0.3045D+01,0.9605D+00,0.3030D+00,0.8810D-01,0.1663D-01/)
      integer, parameter :: nanns=0,nmons=1
      integer ann_units(nanns),mon_units(nmons),imon(nmons)
      integer i,j,iu,k,l
      integer      :: jdlast=0  
      logical :: ifirst=.true.
      character*80 title   
      character*10 :: mon_files(nmons) = (/'SULFATE_SA'/)
      logical      :: mon_bins(nmons)=(/.true./) ! binary file?
      real*8 tlca(im,jm,Lsulf,nmons),tlcb(im,jm,Lsulf,nmons),frac
      REAL*8, DIMENSION(LM)    :: pres,srcLout
      REAL*8, DIMENSION(Lsulf) :: srcLin
      REAL*8, DIMENSION(IM,JM,Lsulf,1) :: src
      save jdlast,tlca,tlcb,mon_units,imon,ifirst
c
C**** Sulfate surface area input file is monthly, on LM levels.
C     Read it in here and interpolated each day.
C     I belive no conversion of the sulfate input file is expected,
C     other than interpolation in the vertical.
C
      if (ifirst) call openunits(mon_files,mon_units,mon_bins,nmons)
      ifirst = .false.
      j = 0
      do k=nanns+1,1
        j = j+1
        call read_monthly_3Dsources(Lsulf,mon_units(j),jdlast,
     *    tlca(1,1,1,j),tlcb(1,1,1,j),src(1,1,1,k),frac,imon(j))
      end do
      jdlast = jday
      write(6,*) 'Sulfate surface area interpolated to current day',frac
      call sys_flush(6) 
C====
C====   Place sulfate onto model levels ("sulfate" array to be used in
C====   the chemistry):
C====              
      PRES(:)=SIG(:)*(PSF-PTOP)+PTOP
      do k=nanns+1,1; DO J=1,JM; DO I=1,IM
        DO L=1,Lsulf
          srcLin(L)=src(I,J,L,k)
        END DO
        call LOGPINT(Lsulf,Psulf,srcLin,LM,PRES,srcLout,.true.)
        DO L=1,LM
          sulfate(I,J,L)=srcLout(L)
        ENDDO         
      end do        ; END DO   ; END DO     
C
      return
      END SUBROUTINE get_sulfate_N2O5
C
C
      SUBROUTINE get_dms_offline
!@sum get_dms_offline read in monthly mean fields of DMS concentration
!@+ for use in chemistry. This source is not applied directly to the tracer
!@+    mass like other 3D sources are but used for HOx chemistry.
!@auth Drew Shindell / Greg Faluvegi
!@ver  1.0 
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only: jday,im,jm,lm,ptop,psf,sig
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      use TRACER_SOURCES, only: Lsulf
      USE TRCHEM_Shindell_COM, only: dms_offline
C
      IMPLICIT NONE

C**** Local parameters and variables and arguments:
c
!@var nanns,nmons: number of annual and monthly input files
!@param Psulf pressure levels of the DMS input file
      real*8, parameter, dimension(Lsulf) :: Psulf = (/
     & 0.9720D+03,0.9445D+03,0.9065D+03,
     & 0.8515D+03,0.7645D+03,0.6400D+03,0.4975D+03,0.3695D+03,
     & 0.2795D+03,0.2185D+03,0.1710D+03,0.1335D+03,0.1016D+03,
     & 0.7120D+02,0.4390D+02,0.2470D+02,0.1390D+02,0.7315D+01,
     & 0.3045D+01,0.9605D+00,0.3030D+00,0.8810D-01,0.1663D-01/)
      integer, parameter :: nanns=0,nmons=1
      integer ann_units(nanns),mon_units(nmons),imon(nmons)
      integer i,j,iu,k,l
      integer      :: jdlast=0  
      logical :: ifirst=.true.
      character*80 title   
      character*10 :: mon_files(nmons) = (/'DMS_FIELD'/)
      logical      :: mon_bins(nmons)=(/.true./) ! binary file?
      real*8 tlca(im,jm,Lsulf,nmons),tlcb(im,jm,Lsulf,nmons),frac
      REAL*8, DIMENSION(LM)    :: pres,srcLout
      REAL*8, DIMENSION(Lsulf) :: srcLin
      REAL*8, DIMENSION(IM,JM,Lsulf,1) :: src
      save jdlast,tlca,tlcb,mon_units,imon,ifirst
c
C**** DMS input file is monthly, on LM levels.
C     Read it in here and interpolated each day.
C     Interpolation in the vertical.
C
      if (ifirst) call openunits(mon_files,mon_units,mon_bins,nmons)
      ifirst = .false.
      j = 0
      do k=nanns+1,1
        j = j+1
        call read_monthly_3Dsources(Lsulf,mon_units(j),jdlast,
     *    tlca(1,1,1,j),tlcb(1,1,1,j),src(1,1,1,k),frac,imon(j))

      end do
      jdlast = jday
      write(6,*) 'DMS concentration interpolated to current day',frac
      call sys_flush(6) 
C====
C====   Place DMS onto model levels
C====              
      PRES(:)=SIG(:)*(PSF-PTOP)+PTOP
      do k=nanns+1,1; DO J=1,JM; DO I=1,IM
        DO L=1,Lsulf
          srcLin(L)=src(I,J,L,k)
        END DO
        call LOGPINT(Lsulf,Psulf,srcLin,LM,PRES,srcLout,.true.)
        DO L=1,LM
          dms_offline(I,J,L)=srcLout(L)

        ENDDO         
      end do        ; END DO   ; END DO     
C
      return
      END SUBROUTINE get_dms_offline

      SUBROUTINE get_so2_offline
!@sum get_so2_offline read in monthly mean fields of SO2 concentration
!@+ for use in chemistry. This source is not applied directly to the tracer
!@+    mass like other 3D sources are but used for HOx chemistry.
!@auth Drew Shindell / Greg Faluvegi
!@ver  1.0 
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only: jday,im,jm,lm,ptop,psf,sig
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      use TRACER_SOURCES, only: Lsulf
      USE TRCHEM_Shindell_COM, only: so2_offline
C
      IMPLICIT NONE

C**** Local parameters and variables and arguments:
c
!@var nanns,nmons: number of annual and monthly input files
!@param Psulf pressure levels of the SO2 input file
      real*8, parameter, dimension(Lsulf) :: Psulf = (/
     & 0.9720D+03,0.9445D+03,0.9065D+03,
     & 0.8515D+03,0.7645D+03,0.6400D+03,0.4975D+03,0.3695D+03,
     & 0.2795D+03,0.2185D+03,0.1710D+03,0.1335D+03,0.1016D+03,
     & 0.7120D+02,0.4390D+02,0.2470D+02,0.1390D+02,0.7315D+01,
     & 0.3045D+01,0.9605D+00,0.3030D+00,0.8810D-01,0.1663D-01/)
      integer, parameter :: nanns=0,nmons=1
      integer ann_units(nanns),mon_units(nmons),imon(nmons)
      integer i,j,iu,k,l
      integer      :: jdlast=0  
      logical :: ifirst=.true.
      character*80 title   
      character*10 :: mon_files(nmons) = (/'SO2_FIELD'/)
      logical      :: mon_bins(nmons)=(/.true./) ! binary file?
      real*8 tlca(im,jm,Lsulf,nmons),tlcb(im,jm,Lsulf,nmons),frac
      REAL*8, DIMENSION(LM)    :: pres,srcLout
      REAL*8, DIMENSION(Lsulf) :: srcLin
      REAL*8, DIMENSION(IM,JM,Lsulf,1) :: src
      save jdlast,tlca,tlcb,mon_units,imon,ifirst
c
C**** Sulfate surface area input file is monthly, on LM levels.
C     Read it in here and interpolated each day.
C     I belive no conversion of the sulfate input file is expected,
C     other than interpolation in the vertical.
C
      if (ifirst) call openunits(mon_files,mon_units,mon_bins,nmons)
      ifirst = .false.
      j = 0
      do k=nanns+1,1
        j = j+1
        call read_monthly_3Dsources(Lsulf,mon_units(j),jdlast,
     *    tlca(1,1,1,j),tlcb(1,1,1,j),src(1,1,1,k),frac,imon(j))

      end do
      jdlast = jday
      write(6,*) 'SO2 concentration interpolated to current day',frac
      call sys_flush(6) 
C====
C====   Place SO2 onto model levels
C====              
      PRES(:)=SIG(:)*(PSF-PTOP)+PTOP
      do k=nanns+1,1; DO J=1,JM; DO I=1,IM
        DO L=1,Lsulf
          srcLin(L)=src(I,J,L,k)
        END DO
        call LOGPINT(Lsulf,Psulf,srcLin,LM,PRES,srcLout,.true.)
        DO L=1,LM
          so2_offline(I,J,L)=srcLout(L)

        ENDDO         
      end do        ; END DO   ; END DO     
C
      return
      END SUBROUTINE get_so2_offline
C
C
      SUBROUTINE read_monthly_3Dsources(Ldim,iu,jdlast,tlca,tlcb,data1,
     & frac,imon)
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
      implicit none
!@var Ldim how many vertical levels in the read-in file?
!@var L dummy vertical loop variable
      integer Ldim,L
      real*8 frac, A2D(im,jm), B2D(im,jm)
      real*8 tlca(im,jm,Ldim),tlcb(im,jm,Ldim),data1(im,jm,Ldim)
      integer imon,iu,jdlast
C
      if (jdlast.EQ.0) then   ! NEED TO READ IN FIRST MONTH OF DATA
        imon=1                ! imon=January
        if (jday.le.16)  then ! JDAY in Jan 1-15, first month is Dec
          do L=1,LDim*11
            read(iu)
          end do
          DO L=1,Ldim
            call readt(iu,0,A2D,im*jm,A2D,1)
            tlca(:,:,L)=A2d(:,:)
          END DO
          rewind iu
        else              ! JDAY is in Jan 16 to Dec 16, get first month
  120     imon=imon+1
          if (jday.gt.idofm(imon) .AND. imon.le.12) go to 120
          do L=1,Ldim*(imon-2)
            read(iu)
          end do
          DO L=1,Ldim
            call readt(iu,0,A2D,im*jm,A2D,1)
            tlca(:,:,L)=A2d(:,:)
          END DO
          if (imon.eq.13)  rewind iu
        end if
      else                         ! Do we need to read in second month?
        if (jday.ne.jdlast+1) then ! Check that data is read in daily
          if (jday.ne.1 .OR. jdlast.ne.365) then
            write(6,*)
     *      'Bad values in Tracer 3D Source:JDAY,JDLAST=',JDAY,JDLAST
            call stop_model('Bad values in Tracer 3D Source.',255)
          end if
          imon=imon-12             ! New year
          go to 130
        end if
        if (jday.le.idofm(imon)) go to 130
        imon=imon+1                ! read in new month of data
        tlca(:,:,:) = tlcb(:,:,:)
        if (imon.eq.13) rewind iu
      end if
      DO L=1,Ldim
        call readt(iu,0,B2D,im*jm,B2D,1)
        tlcb(:,:,L)=B2D(:,:)
      END DO
  130 continue
c**** Interpolate two months of data to current day
      frac = float(idofm(imon)-jday)/(idofm(imon)-idofm(imon-1))
      data1(:,:,:) = tlca(:,:,:)*frac + tlcb(:,:,:)*(1.-frac)
      return
      end subroutine read_monthly_3Dsources


      subroutine get_CH4_IC
!@sum get_CH4_IC to generate initial conditions for methane.
!@auth Greg Faluvegi/Drew Shindell
!@ver  1.0 (based on DB396Tds3M23)
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only  : im,jm,lm,ls1,JEQ
      USE GEOM, only       : dxyp
      USE DYNAMICS, only   : am
      USE CONSTANT, only: mair
      USE TRACER_COM, only : trm, TR_MM, n_CH4
      USE TRCHEM_Shindell_COM, only: CH4altT, CH4altX, ch4_init_sh,
     *     ch4_init_nh
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
c
!@var CH4INIT temp variable for ch4 initial conditions
!@var I,J,L dummy loop variables
!@param bymair 1/molecular wt. of air = 1/mair
!@var JN J around 30 N
!@var JS J arount 30 S
      INTEGER, PARAMETER:: JS = JM/3 + 1, JN = 2*JM/3
      REAL*8, PARAMETER :: bymair = 1.d0/mair
      REAL*8 CH4INIT
      INTEGER I, J, L

C     First, the troposphere:
      DO L=1,LS1-1
      DO J=1,JM
C       Initial latitudinal gradient for CH4:
        IF(J.LT.JEQ)THEN ! Southern Hemisphere
          CH4INIT=ch4_init_sh*TR_MM(n_CH4)*bymair*1.d-6*DXYP(J)
        ELSE             ! Northern Hemisphere
          CH4INIT=ch4_init_nh*TR_MM(n_CH4)*bymair*1.d-6*DXYP(J)
        ENDIF
        DO I=1,IM
          trm(i,j,l,n_CH4)=am(L,I,J)*CH4INIT
        END DO
      END DO
      END DO
c
C     Now, the stratosphere:
      do L=LS1,LM
      do j=1,jm
c
c     Define stratospheric ch4 based on HALOE obs for tropics
c     and extratropics and scale by the ratio of initial troposphere
c     mixing ratios to 1.79 (observed):
        IF(J.LT.JEQ)THEN ! Southern Hemisphere
          CH4INIT=ch4_init_sh/1.79d0*TR_MM(n_CH4)*bymair*1.E-6*DXYP(J)
        ELSE             ! Northern Hemisphere
          CH4INIT=ch4_init_nh/1.79d0*TR_MM(n_CH4)*bymair*1.E-6*DXYP(J)
        ENDIF
        IF((J.LE.JS).OR.(J.GT.JN)) THEN                 ! extratropics
          CH4INIT=CH4INIT*CH4altX(L)
        ELSE IF((J.GT.JS).AND.(J.LE.JN)) THEN           ! tropics
          CH4INIT=CH4INIT*CH4altT(L)
        END IF
        do i=1,im  
          trm(i,j,l,n_CH4)= am(L,I,J)*CH4INIT
        end do ! i
      end do   ! j
      end do   ! l
C
      RETURN
      end subroutine get_CH4_IC
C  
C     
      SUBROUTINE calc_lightning(I,J,LMAX,LFRZ)
!@sum calc_lightning to calculate lightning flash amount and cloud-
!@+   to-ground amount, based on cloud top height. WARNING: this 
!@+   routine is apparently resolution dependant.  See the comments.
!@auth Colin Price (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436Tds3M23)
C
C**** GLOBAL parameters and variables:
C
      USE LIGHTNING, only : JN,JS,RNOx_lgt
      USE MODEL_COM, only : fland
      USE GEOM,      only : bydxyp
      USE CONSTANT,  only : bygrav
      USE DYNAMICS,  only : gz
      USE TRACER_DIAG_COM, only : ijs_CtoG,ijs_flash,taijs
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
C
!@var LMAX highest layer of current convective cloud event
!@var lmax_temp local copy of LMAX (alterable)
!@var LFRZ freezing level
!@var HTCON height of convection?
!@var HTFRZ height of freezing?
!@var lightning flashes per minute
!@var th,th2,th3,th4 thickness of cold sector, squared, cubed, etc.
!@var CG fraction of lightning that is cloud-to-ground
!@var zlt ?
!@var I model longitude index
!@var J model latitude index
!@param tune_land multiplier of flash rate over land to match observ.
!@param tune_ocean multiplier of flash rate over land to match observ.
!@param tune_NOx   multiplier of NOx production rate from lightning
      INTEGER, INTENT(IN) :: LMAX,LFRZ,I,J
      INTEGER lmax_temp
      REAL*8 HTCON,HTFRZ,flash,th,th2,th3,th4,zlt,CG
      REAL*8,PARAMETER::tune_land=2.2d0,tune_ocean=3.9d0,tune_NOx=0.67d0
C
c The folowing simple algorithm calculates the lightning
c frequency in each gridbox using the moist convective cloud
c scheme.  The algorithm is based on the Price and Rind (1992)
c JGR paper, and differentiates between oceanic and continental
c convective clouds.  Only the non-entraining plume is used for
c the lightning.  However, I [Greg] have removed the IC variable
C from this routine because Gavin says: "IC distinguishes between
C entraining and non-entraining, but since non-entraining always
C go higher, those are the ones represented by LMCMAX."
C
c The lightning calculation uses the maximum height of the cloud,
c or the maximum depth of the mass flux penetration, LMAX, at
c every timestep.  In each timestep there may be a number of
c different values of LMAX, associated with clouds originating
c from different levels.  However, there are always less values
c of LMAX then the number of vertical layers. gz is the
c geopotential height, used from DYNAMICS module.

      lmax_temp=lmax
      if(j.lt.JS.or.j.gt.JN)lmax_temp=lmax+1 !30 S to 30 N
      HTCON=gz(i,j,lmax_temp)*bygrav*1.d-3
      HTFRZ=gz(i,j,LFRZ)*bygrav*1.d-3

c IF the gridbox is over land or coastline use the continental
c parameterization.  The factor 2.17 is derived from Fig. 1
c of the Price and Rind (1994) Mon. Wea. Rev. paper, and is
c related to the horizontal resolution of the GCM (in this case
c 4x5 degrees).  We use a threshold of 1% for continental
c gridboxes. The units are flashes per minute.

      If (fland(i,j).le.0.01) then
        flash=tune_ocean*2.17d0*6.2d-4*(htcon**1.73d0)  ! ocean
      else
        flash= tune_land*2.17d0*3.44d-5*(htcon**4.92d0) ! continent
      end if

c The formulation by Price and Rind (1993) in Geophysical Res.
c Letters is used to calculate the fraction of total lightning
c that is cloud-to-ground lightning.  This uses the thickness of
c cold sector of the cloud (Hmax-Hzero) as the determining
c parameter.  The algorithm is only valid for thicknesses greater
c than 5.5km and less than 14km.

      th=(HTCON-HTFRZ)
      th=min(max(th,5.5d0),14.d0)
      th2=th*th
      th3=th2*th
      th4=th3*th
      zlt = 0.021d0*th4 - 0.648d0*th3 + 7.493d0*th2 
     &    - 36.544d0*th + 63.088d0
      CG=flash/(1.+zlt)

c Given the number of cloud-to-ground flashes, we can now calculate
c the NOx production rate based on the Price et al. (1997) JGR paper.
c The units are grams of Nitrogen per minute:

      RNOx_lgt(i,j)=15611.d0*(CG + 0.1d0*(flash-CG))*tune_NOx

C If flash is indeed in flashes/min, accumulate it in flashes:
       TAIJS(I,J,ijs_flash)=TAIJS(I,J,ijs_flash) + flash*60.d0
       TAIJS(I,J,ijs_CtoG) =TAIJS(I,J,ijs_CtoG)  +    CG*60.d0

      END SUBROUTINE calc_lightning
c
c
      SUBROUTINE get_sza(I,J,tempsza)
!@sum get_sza calculates the solar angle.  The intention is
!@+   that this routine will only be used when the COSZ1 from the 
!@+   radiation code is < 0, i.e. SZA > 90 deg.
!@auth Greg Faluvegi, based on the Harvard CTM routine SCALERAD
!@ver 1.0

c
C**** GLOBAL parameters and variables:  
C
      USE MODEL_COM, only: IM,JM,nday,Itime,DTsrc 
      USE CONSTANT, only: PI, radian
      USE TRCHEM_Shindell_COM, only: byradian
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments
C
!@param ANG1 ?
!@var DX degree width of a model grid cell (360/IM)
!@var DY degree width of a model grid cell ~(180/JM)
!@var tempsza the solar zenith angle, returned in degrees
!@var P1,P2,P3 ? angles needed to compute COS(SZA) in degrees    
!@var VLAT,VLOM latitude and longitude in degrees
!@var current julian time in seconds
!@var FACT,temp are temp variables
      REAL*8, PARAMETER ::  ANG1 = 90.d0/91.3125d0
      REAL*8 P1,P2,P3,VLAT,VLON,TIMEC,FACT,temp
      REAL*8, INTENT(OUT) :: tempsza
      INTEGER, INTENT(IN) :: I,J
      INTEGER DX,DY

      DX   = NINT(360./REAL(IM))
      DY   = NINT(180./REAL(JM))       
      vlat = -90. + REAL((J-1)*DY)
      if (J.eq.1)  vlat= -90. + 0.5*REAL(DY)
      if (j.eq.JM) vlat=  90. - 0.5*REAL(DY)
      VLON = 180. - REAL((I-1)*DX)
C     This added 0.5 is to make the instantaneous zenith angle
C     more representative throughout the 1 hour time step:
!old  TIMEC = ((JDAY*24.0) + JHOUR  + 0.5)*3600.
      TIMEC = (mod(itime,365*nday) + nday + 0.5)*DTsrc  
!     if(DTsrc.ne.3600.)call stop_model
!    &('DTsrc has changed. Please check on get_sza subroutine',255)
      P1 = 15.*(TIMEC/3600. - VLON/15. - 12.)
      FACT = (TIMEC/86400. - 81.1875)*ANG1
      P2 = 23.5*SIN(FACT*radian)
      P3 = VLAT
      temp = (SIN(P3*radian)*SIN(P2*radian)) +
     &(COS(P1*radian)*COS(P2*radian)*COS(P3*radian))
      tempsza = acos(temp)*byradian

      RETURN
      END SUBROUTINE get_sza
#endif
C
C
      SUBROUTINE LOGPINT(LIN,PIN,AIN,LOUT,POUT,AOUT,min_zero)
!@sum LOGPINT does vertical interpolation of column variable,
!@+   linearly in ln(P).
!@auth Greg Faluvegi
!@ver 1.0
C
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments
C
!@var LIN number of levels for initial input variable
!@var LOUT number of levels for output variable
!@var PIN pressures at LIN levels
!@var POUT pressures at LOUT levels
!@var AIN initial input column variable
!@var AOUT output (interpolated) column variable
!@var LNPIN natural log of PIN
!@var LNPOUT natural log of POUT
!@var min_zero if true, don't allow negatives
!@var slope slope of line used for extrapolations
C
      LOGICAL, INTENT(IN)                 :: min_zero
      INTEGER, INTENT(IN)                 :: LIN, LOUT
      REAL*8, INTENT(IN), DIMENSION(LIN)  :: PIN, AIN
      REAL*8, INTENT(IN),DIMENSION(LOUT)  :: POUT
      REAL*8, INTENT(OUT),DIMENSION(LOUT) :: AOUT 
      REAL*8, DIMENSION(LIN)              :: LNPIN
      REAL*8, DIMENSION(LOUT)             :: LNPOUT       
      INTEGER L1,L2
      REAL*8 slope
C      
      LNPIN(:) = LOG(PIN(:))                   ! take natural log
      LNPOUT(:)= LOG(POUT(:))                  ! of pressures
C      
      DO L1=1,LOUT
       IF (LNPOUT(L1).gt.LNPIN(1)) THEN        ! extrapolate
         slope=(AIN(2)-AIN(1))/(LNPIN(2)-LNPIN(1))
         AOUT(L1)=AIN(1)-slope*(LNPIN(1)-LNPOUT(L1))
       ELSE IF (LNPOUT(L1).lt.LNPIN(LIN)) THEN ! extrapolate
         slope=(AIN(LIN)-AIN(LIN-1))/(LNPIN(LIN)-LNPIN(LIN-1))
         AOUT(L1)=AIN(LIN)+slope*(LNPOUT(L1)-LNPIN(LIN))
       ELSE                                    ! interpolate
        DO L2=1,LIN-1
         IF(LNPOUT(L1).eq.LNPIN(L2)) THEN
           AOUT(L1)=AIN(L2)
         ELSE IF(LNPOUT(L1).eq.LNPIN(L2+1)) THEN
           AOUT(L1)=AIN(L2+1)
         ELSE IF(LNPOUT(L1).lt.LNPIN(L2) .and.
     &   LNPOUT(L1).gt.LNPIN(L2+1)) THEN
           AOUT(L1)=(AIN(L2)*(LNPIN(L2+1)-LNPOUT(L1)) 
     &     +AIN(L2+1)*(LNPOUT(L1)-LNPIN(L2)))/(LNPIN(L2+1)-LNPIN(L2))
         END IF
        END DO
       END IF
      END DO          
C      
C If necessary: limit interpolated array to positive numbers:
C
      IF(min_zero)AOUT(:)=MAX(0.d0,AOUT(:))
C           
      RETURN
      END SUBROUTINE LOGPINT
C
C
#ifdef TRACERS_ON
      SUBROUTINE special_layers_init
!@sum special_layers_init determined special altitude levels that
!@+ are important to Drew Shindell's tracer code.  This is done to
!@+ so that hardcoded levels are not needed when switching
!@+ vertical resolutions.
!@auth Greg Faluvegi
!@ver 1.0
C
C**** Global variables:
c
      USE MODEL_COM, only : LM,SIG,PSF,PTOP
      USE TRCHEM_Shindell_COM, only: L75P,L75M,F75P,F75M, ! FACT1
     &                           L569P,L569M,F569P,F569M  ! CH4
C
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments
C
!@var natural log of nominal pressure for verticle interpolations
!@var log75 natural log of 75 hPa
!@var log569 natural log of 569 hPa
      REAL*8 log75,log569
      REAL*8, DIMENSION(LM) :: LOGP
      INTEGER L
c      
      LOGP(:)=LOG(SIG(:)*(PSF-PTOP)+PTOP)
      log75=LOG(75.d0)
      log569=LOG(569.d0)
c 
      DO L=1,LM-1
        IF(LOGP(L).gt.log75 .and. LOGP(L+1).lt.log75) THEN
          L75P=L+1 ! these are for FACT1 variable in strat overwrite
          L75M=L   ! hence effects several tracers
          F75P=(log75-LOGP(L75M))/(LOGP(L75P)-LOGP(L75M))
          F75M=(LOGP(L75P)-log75)/(LOGP(L75P)-LOGP(L75M))
        END IF
        IF(LOGP(L).gt.log569 .and. LOGP(L+1).lt.log569) THEN
          L569P=L+1 ! these are for the CH4 strat overwrite
          L569M=L
          F569P=(log569-LOGP(L569M))/(LOGP(L569P)-LOGP(L569M))
          F569M=(LOGP(L569P)-log569)/(LOGP(L569P)-LOGP(L569M))
        END IF
      END DO
c      
      RETURN
      END SUBROUTINE special_layers_init
#endif


#ifdef INTERACTIVE_WETLANDS_CH4
      subroutine running_average(var_in,I,J,avg_out,nicall,m)
!@sum running_average keeps a running average of the model 1st layer
!@+ ground temperature and precipitaion for use with the interactive
!@+ wetlands CH4.  I suppose I could generalized this in the future.
!@auth Greg Faluvegi
!@ver 1.0
C
C**** Global variables:
c
      USE TRACER_SOURCES, only: iH=>iHch4,iD=>iDch4,i0=>i0ch4,
     & first_mod,HRA=>HRA_ch4,DRA=>DRA_ch4,PRS=>PRS_ch4,
     & nday_ncep,by_nday_ncep,nra_ch4
C
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments
C     
!@var var_in model variable which will be used in the average
!@var avg_out the updated running average to be returned
!@var nicall number of times routine is called per main timesetep
!@var m for now, 1=precipitation, 2=1st layer ground temperature
!@var temp just for holding current day average for use in avg_out
      real*8, intent(IN) :: var_in
      real*8, intent(OUT):: avg_out
      real*8 temp, by24N
      integer, intent(IN):: nicall, m, I, J
      integer n
            
      if(m.gt.nra_ch4)call stop_model('nra_ch4 problem',255)
      iH(I,J,m) = iH(I,J,m) + 1
      HRA(I,J,iH(I,J,m),m) = var_in
      ! do no more, unless it is the end of the day:
      
      if(iH(I,J,m).eq.24*nicall)then ! end of "day":
        by24N=1.d0/real(24*nicall)
        iH(I,J,m) = 0
        if(first_mod(m))then ! first averaging period only
          iD(I,J,m) = ID(I,J,m) + 1 
          do n=1,24*nicall
            DRA(I,J,iD(I,J,m),m) = DRA(I,J,iD(I,J,m),m) + HRA(I,J,n,m)
          end do
          DRA(I,J,iD(I,J,m),m) = DRA(I,J,iD(I,J,m),m)*by24N
          if(iD(I,J,m).eq.nday_ncep(m))then !end first period
            PRS(I,J,m) = 0.d0
            do n=1,nday_ncep(m)
              PRS(I,J,m) = PRS(I,J,m) + DRA(I,J,iD(I,J,m),m)
            end do
            avg_out = PRS(I,J,m) * by_nday_ncep(m)
            first_mod(m)=.false.
            iD(I,J,m)=0
            i0(I,J,m)=0
          end if
        else ! not first averaging period: update the running average
          i0(I,J,m) = i0(I,J,m) + 1 ! move pointer
          if(i0(I,J,m) .eq. nday_ncep(m)+1) i0(I,J,m)=1
          temp=0.d0
          do n=1,24*nicall
            temp = temp + HRA(I,J,n,m) 
          end do
          temp = temp * by24N ! i.e. today's average
          PRS(I,J,m) = PRS(I,J,m) - DRA(I,J,i0(I,J,m),m)
          DRA(I,J,i0(I,J,m),m) = temp
          PRS(I,J,m) = PRS(I,J,m) + DRA(I,J,i0(I,J,m),m)
          avg_out = PRS(I,J,m) * by_nday_ncep(m)
        end if
      end if
  
      END SUBROUTINE running_average
#endif
