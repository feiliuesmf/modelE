      MODULE TRACER_SOURCES
      USE TRACER_COM
!@param nlightning index number of 3D source for NOx from lightning
!@param naircraft index number of 3D source for NOx from aircraft
!@param nStratwrite idx number of 3D tracer source from strat overwrite
!@param nChemistry index number of 3D tracer source from chemistry
!@param Laircr the number of layers of aircraft data read from file
!@param Lsulf the number of layers of sulfate SA data read from file
!@param correct_CO_ind correction factor Lee put in for industrial CO
      INTEGER, PARAMETER :: nChemistry  = 1,
     &                      nStratwrite = 2,
     &                      nLightning  = 3,
     &                      nAircraft   = 4,
     &                      Laircr      =19,
     &                      Lsulf       =LM
      REAL*8, PARAMETER :: correct_CO_ind=520./410.
c
!@var CH4_src           CH4 surface sources and sinks (kg/s)
!@var CO_src             CO surface sources and sinks (kg/s)
!@var Alkenes_src   Alkenes surface sources and sinks (kg/s)
!@var Paraffin_src Paraffin surface sources and sinks (kg/s)
!@var NOx_src           NOx surface sources and sinks (kg/s)
!@var Isoprene_src ISoprene surface sources and sinks (kg/s)
c
c SHOULD PROBABLY USE ntsurfsrc( ) instead of these ...
      integer, parameter :: nch4src     = 14,
     &                      nCOsrc      =  2,
     &                      nAlkenessrc =  3,
     &                      nParaffinsrc=  3,
     &                      nIsoprenesrc=  1,
     &                      nNOxsrc     =  3
      real*8 CO_src(im,jm,nCOsrc),CH4_src(im,jm,nch4src),
     & Alkenes_src(im,jm,nAlkenessrc),Paraffin_src(im,jm,nParaffinsrc),
     & Isoprene_src(im,jm,nIsoprenesrc),NOx_src(im,jm,nNOxsrc)
      END MODULE TRACER_SOURCES


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
      implicit none
      character*80 title
!@var adj Factors that tune the total amount of individual sources
      logical :: ifirst=.true.
      real*8 adj(nsrc)
      data adj/1.3847,1.0285,3.904,1.659,1.233,1.194,0.999,
!    *  7.2154, 3.7247d-5,3.1399d-4,5.4838d-5,   !Model II prime
     *  7.2154, 3.5997d-5,17.330d-4,5.3558d-5,
     *  0.4369,0.7533,0.9818/
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
      data kwet/14/  !!! position of wetlands array in src
      logical :: mon_bins(nmons)=(/.true.,.true.,.true./)
      real*8 tlca(im,jm,nmons),tlcb(im,jm,nmons)  ! for monthly sources
      real*8 frac
      integer i,j,nt,iact,iu,k,kwet,imon(nmons)
      integer :: jdlast=0
      save ifirst,jdlast,tlca,tlcb,mon_units,imon

      if (itime.lt.itime_tr0(nt)) return
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
          call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          src(:,:,k) = src(:,:,k)*adj(k)/(sday*JDperY)
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
        j = j+1
        call read_monthly_sources(mon_units(j),jdlast,
     *    tlca(1,1,j),tlcb(1,1,j),src(1,1,k),frac,imon(j))
        src(:,:,k) = src(:,:,k)*adj(k)
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
      implicit none
      character*80 title
      logical :: ifirst=.true.
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
          call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          src(:,:,k) = src(:,:,k)*correct_CO_ind
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
        j = j+1
        call read_monthly_sources(mon_units(j),jdlast,
     *    tlca(1,1,j),tlcb(1,1,j),src(1,1,k),frac,imon(j))
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
      USE MODEL_COM, only: itime,jday,JDperY,im,jm
      USE CONSTANT, only: sday,hrday
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      USE GEOM, only: BYDXYP
      USE TRACER_COM, only: itime_tr0,trname
      use TRACER_SOURCES, only: src=>Alkenes_src,nsrc=>nAlkenessrc
      implicit none
      character*80 title
      logical :: ifirst=.true.
!@var nanns,nmons: number of annual and monthly input files
!@param r3600 = 1.d0/3600.
      real*8, parameter :: r3600=1./3600.
      integer, parameter :: nanns=1,nmons=2
      integer ann_units(nanns),mon_units(nmons)
      character*18 :: ann_files(nanns) = (/'Alkenes_INDUSTRIAL'/)
      logical :: ann_bins(nanns)=(/.true./) ! binary files?
      character*18 :: mon_files(nmons) = (/'Alkenes_BIOMASS   ',
     &                                     'Alkenes_VEGETATION'/)
      logical :: mon_bins(nmons)=(/.true.,.true./) ! binary files?
      real*8 tlca(im,jm,nmons),tlcb(im,jm,nmons)  ! for monthly sources
      real*8 frac
      integer i,j,jj,nt,iact,iu,k,imon(nmons)
      integer :: jdlast=0
      save ifirst,jdlast,tlca,tlcb,mon_units,imon

      if (itime.lt.itime_tr0(nt)) return
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
          call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          do j = 1,jm
            src(:,j,k) = src(:,j,k)*BYDXYP(j)*r3600
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
        j = j+1
        call read_monthly_sources(mon_units(j),jdlast,
     *    tlca(1,1,j),tlcb(1,1,j),src(1,1,k),frac,imon(j))
          do jj = 1,jm
            src(:,jj,k) = src(:,jj,k)*BYDXYP(jj)*r3600
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
      USE MODEL_COM, only: itime,jday,JDperY,im,jm
      USE CONSTANT, only: sday,hrday
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      USE GEOM, only: BYDXYP
      USE TRACER_COM, only: itime_tr0,trname
      use TRACER_SOURCES, only: src=>Paraffin_src,nsrc=>nParaffinsrc
      implicit none
      character*80 title
      logical :: ifirst=.true.
!@var nanns,nmons: number of annual and monthly input files
!@var r3600 = 1.d0/3600.
      integer, parameter :: nanns=1,nmons=2
      real*8, parameter :: r3600=1./3600.
      integer ann_units(nanns),mon_units(nmons)
      character*19 :: ann_files(nanns) = (/'Paraffin_INDUSTRIAL'/)
      logical :: ann_bins(nanns)=(/.true./) ! binary files?
      character*19 :: mon_files(nmons) = (/'Paraffin_BIOMASS   ',
     &                                     'Paraffin_VEGETATION'/)
      logical :: mon_bins(nmons)=(/.true.,.true./) ! binary files?
      real*8 tlca(im,jm,nmons),tlcb(im,jm,nmons)  ! for monthly sources
      real*8 frac
      integer i,j,jj,nt,iact,iu,k,imon(nmons)
      integer :: jdlast=0
      save ifirst,jdlast,tlca,tlcb,mon_units,imon

      if (itime.lt.itime_tr0(nt)) return
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
          call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          do j = 1,jm
            src(:,j,k) = src(:,j,k)*BYDXYP(j)*r3600
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
        j = j+1
        call read_monthly_sources(mon_units(j),jdlast,
     *  tlca(1,1,j),tlcb(1,1,j),src(1,1,k),frac,imon(j))
        do jj = 1,jm
          src(:,jj,k) = src(:,jj,k)*BYDXYP(jj)*r3600
        end do
      end do
      jdlast = jday
      write(6,*) trname(nt),'Sources interpolated to current day',frac
      call sys_flush(6)
C****
      return
      end subroutine read_Paraffin_sources


      subroutine read_Isoprene_sources(nt,iact)
!@sum reads in Isoprene surface sources and sinks
!@auth Jean Lerner/Greg Faluvegi
C****
C**** There is 1 monthly source
C**** Monthly sources are interpolated each day
      USE MODEL_COM, only: itime,jday,JDperY,im,jm
      USE CONSTANT, only: sday,hrday
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      USE GEOM, only: BYDXYP
      USE TRACER_COM, only: itime_tr0,trname
      use TRACER_SOURCES, only: src=>Isoprene_src,nsrc=>nIsoprenesrc
      implicit none
      character*80 title
      logical :: ifirst=.true.
!@var nanns,nmons: number of annual and monthly input files
!@var r3600 = 1.d0/3600.
      real*8, parameter :: r3600=1./3600.
      integer, parameter :: nanns=0, nmons=1
      integer mon_units(nmons)
      character*19 :: mon_files(nmons) = (/'Isoprene_VEGETATION'/)
      logical :: mon_bins(nmons)=(/.true./) ! binary file?
      real*8 tlca(im,jm,nmons),tlcb(im,jm,nmons)  ! for monthly sources
      real*8 frac
      integer i,j,jj,nt,iact,iu,k,imon(nmons)
      integer :: jdlast=0
      save ifirst,jdlast,tlca,tlcb,mon_units,imon

      if (itime.lt.itime_tr0(nt)) return
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
          src(:,jj,k) = src(:,jj,k)*BYDXYP(jj)*r3600
        end do
      end do
      jdlast = jday
      write(6,*) trname(nt),'Sources interpolated to current day',frac
      call sys_flush(6)
C****
      return
      end subroutine read_Isoprene_sources


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
          call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          src(:,:,k) = src(:,:,k)/(sday*JDperY)
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
        j = j+1
        call read_monthly_sources(mon_units(j),jdlast,
     *    tlca(1,1,j),tlcb(1,1,j),src(1,1,k),frac,imon(j))
      end do
      jdlast = jday
      write(6,*) trname(nt),'Sources interpolated to current day',frac
      call sys_flush(6)
C****
      return
      end subroutine read_NOx_sources
c
c
      SUBROUTINE get_lightning_NOx
!@sum  get_lightning_NOx to define the 3D source of NOx from lightning
!@auth Colin Price / Greg Faluvegi
!@ver  1.0 (based on CB436Tds3M23 & DB396Tds3M23)
c
C**** GLOBAL parameters and variables:
      USE FLUXES, only: tr3Dsource
      USE TRACER_COM, only: n_NOx
      USE TRACER_SOURCES, only: nLightning
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
c
      tr3Dsource(:,:,:,nLightning,n_NOx) = 0.
      WRITE(6,*) '>>>get_lightning_NOx still a dummy routine<<<'
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
      USE MODEL_COM, only: itime,jday,JDperY,im,jm,lm,ptop
      USE CONSTANT, only: sday,hrday
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      USE FLUXES, only: tr3Dsource
      USE DYNAMICS, only   : PMID
      USE GEOM, only       : dxyp
      USE TRACER_COM, only: itime_tr0,trname,n_NOx
      use TRACER_SOURCES, only: nAircraft,Laircr
C
      IMPLICIT NONE
c
!@var nanns,nmons: number of annual and monthly input files
!@var l,ll dummy loop variable
!@var pres local pressure variable
      character*80 title
      logical LINJECT
      integer, parameter :: nanns=0,nmons=1
      integer ann_units(nanns),mon_units(nmons)
      character*12 :: mon_files(nmons) = (/'NOx_AIRCRAFT'/)
      logical :: mon_bins(nmons)=(/.true./) ! binary file?
      real*8 tlca(im,jm,Laircr,nmons),tlcb(im,jm,Laircr,nmons)
      real*8 frac
      integer l,i,j,iu,k,ll,imon(nmons)
      integer :: jdlast=0
      save jdlast,tlca,tlcb,mon_units,imon
      REAL*8, DIMENSION(IM,JM,Laircr,1) :: src
      REAL*8, DIMENSION(LM)     :: pres
      REAL*4, DIMENSION(Laircr) :: PAIRL
      DATA PAIRL/1013.25,898.74,794.95,701.08,616.40,540.19,471.81,
     &            410.60,355.99,307.42,264.36,226.32,193.30,165.10,
     &            141.01,120.44,102.87,87.866,75.048/
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
      call openunits(mon_files,mon_units,mon_bins,nmons)
      j = 0
      do k=nanns+1,1
        j = j+1
        call read_monthly_3Dsources(Laircr,mon_units(j),jdlast,
     *    tlca(1,1,1,j),tlcb(1,1,1,j),src(1,1,1,k),frac,imon(j))
      end do
      jdlast = jday
      call closeunits(mon_units,nmons)
      write(6,*)trname(n_NOx),'frm aircraft interpd to currnt day',frac
      call sys_flush(6)
C====
C====   Place aircraft sources onto model levels:
C====
      DO J=1,JM
       DO I=1,IM
        DO L=1,LM
          tr3Dsource(i,j,l,nAircraft,n_NOx) = 0.
          if(L.LT.LM) then
            PRES(L)=PMID(L+1,I,J) !Pressure @ layer top?
          else
            PRES(L)=PTOP
          endif
        ENDDO   ! L
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
      SUBROUTINE get_sulfate
!@sum  get_sulfate to parameterize N2O5 sink on sulfate.  The sulfate( )
!@+    array is then used in the chemistry to calculate a chemical
!@+    change.  I.e. this source is not applied directly to the tracer
!@+    mass like other 3D sources are.
!@auth Drew Shindell? / Greg Faluvegi
!@ver  1.0 (based on DB396Tds3M23)
c
C**** GLOBAL parameters and variables:
C
      USE MODEL_COM, only: jday,im,jm,lm
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      use TRACER_SOURCES, only: Lsulf
      USE TRCHEM_Shindell_COM, only: sulfate
C
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
c
!@var nanns,nmons: number of annual and monthly input files
      character*80 title
      integer, parameter :: nanns=0,nmons=1
      integer ann_units(nanns),mon_units(nmons)
      character*10 :: mon_files(nmons) = (/'SULFATE_SA'/)
      logical :: mon_bins(nmons)=(/.true./) ! binary file?
      real*8 tlca(im,jm,Lsulf,nmons),tlcb(im,jm,Lsulf,nmons)
      real*8 frac
      integer i,j,iu,k,imon(nmons)
      integer :: jdlast=0
      REAL*8, DIMENSION(IM,JM,Lsulf,1) :: src
      save jdlast,tlca,tlcb,mon_units,imon
c
C**** Sulfate surface area input file is monthly, on LM levels.
C     Read it in here and interpolated each day.
C****
C**** Monthly sources are interpolated to the current day
C**** I belive no conversion of the sulfate input file is expected:
C****
      call openunits(mon_files,mon_units,mon_bins,nmons)
      j = 0
      do k=nanns+1,1
        j = j+1
        call read_monthly_3Dsources(Lsulf,mon_units(j),jdlast,
     *    tlca(1,1,1,j),tlcb(1,1,1,j),src(1,1,1,k),frac,imon(j))
      end do
      jdlast = jday
      call closeunits(mon_units,nmons)
      write(6,*) 'Sulfate srface area interpolated to current day',frac
      call sys_flush(6)
C
C**** Save in array to be used in chemistry:
C
      sulfate(:,:,:)=src(:,:,:,1)
C
      return
      END SUBROUTINE get_sulfate
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
          DO L=1,Ldim
            call readt(iu,0,A2D,im*jm,A2D,11*Ldim+L)
            tlca(:,:,L)=A2d(:,:)
            rewind iu
          END DO
        else              ! JDAY is in Jan 16 to Dec 16, get first month
  120     imon=imon+1
          if (jday.gt.idofm(imon) .AND. imon.le.12) go to 120
          DO L=1,Ldim-1
            call readt(iu,0,A2D,im*jm,A2D,(imon-1)*Ldim+L)
            tlca(:,:,L)=A2d(:,:)
            rewind iu
          END DO
          if (imon.eq.13)  rewind iu
        end if
      else                         ! Do we need to read in second month?
        if (jday.ne.jdlast+1) then ! Check that data is read in daily
          if (jday.ne.1 .OR. jdlast.ne.365) then
            write(6,*)
     *      'Bad values in Tracer 3D Source:JDAY,JDLAST=',JDAY,JDLAST
            stop
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
!@+ >>WARNING<<: This routine is still mostly 'hard coded' for the
!@+ 4 x 5, 23 layer model. It needs to be generalized!
!@auth Greg Faluvegi/Drew Shindell
!@ver  1.0 (based on DB396Tds3M23)
c
C**** GLOBAL parameters and variables:
C
      USE RESOLUTION, only : im,jm,lm,ls1
      USE MODEL_COM, only  : JEQ
      USE GEOM, only       : dxyp
      USE DYNAMICS, only   : am
      USE CONSTANT, only: mair
      USE TRACER_COM, only : trm, TR_MM, n_CH4
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
c
!@var CH4INIT temp variable for ch4 initial conditions
!@var j2 dummy
!@var I,J,L dummy loop variables
!@param bymair 1/molecular wt. of air = 1/mair
      REAL*8, PARAMETER :: bymair = 1.d0/mair
      REAL*8 CH4INIT
      INTEGER J2, I, J, L

      WRITE(6,*) '>>WARNING : you are using CH4 initial condition in'
      WRITE(6,*) 'the stratosphere which are specific to the 4x5 deg,'
      WRITE(6,*) '23-layer model. Otherwise, subroutine get_CH4_IC'
      WRITE(6,*) 'needs to be generalized. <<'

C     First, the troposphere:
      DO L=1,LS1-1
      DO J=1,JM
C       Initial latitudinal gradient for CH4:
        IF(J.LT.JEQ)THEN ! Southern Hemisphere
          CH4INIT=1.75 *TR_MM(n_CH4)*bymair*1.E-6*DXYP(J)
        ELSE             ! Northern Hemisphere
          CH4INIT=1.855*TR_MM(n_CH4)*bymair*1.E-6*DXYP(J)
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
        j2=NINT((FLOAT(j)-1.)/3.)
         if(j.le.2) j2=1
         if(j.ge.45)j2=14
c
c     Define stratospheric ch4 based on HALOE obs for tropics
c     and extratropics and scale by the ratio of initial troposphere
c     mixing ratios to 1.79 (observed). The 12 "stratospheric" levels
c     are grouped into 4 categories: L={12,13,14} based on 100mb obs.
c     L={15,16,17} on 32mb obs, L={18,19,20} on 3.2mb obs and
c     L={21,22,23} on 0.32mb obs (average of mar,jun,sep,dec).
        IF(J.LT.JEQ)THEN ! Southern Hemisphere
          CH4INIT=1.75/1.79  *TR_MM(n_CH4)*bymair*1.E-6*DXYP(J)
        ELSE             ! Northern Hemisphere
          CH4INIT=1.855/1.79 *TR_MM(n_CH4)*bymair*1.E-6*DXYP(J)
        ENDIF
        IF((J.LE.16).OR.(J.GT.30)) THEN               ! extratropics
          IF((L.GE.LS1).AND.(L.LE.14)) CH4INIT=CH4INIT*1.700!was1.44
          IF((L.GE.15).AND.(L.LE.17))  CH4INIT=CH4INIT*   1.130
          IF((L.GE.18).AND.(L.LE.20))  CH4INIT=CH4INIT*   0.473
          IF((L.GE.21).AND.(L.LE.LM))  CH4INIT=CH4INIT*   0.202
        ELSE IF((J.GT.16).AND.(J.LE.30)) THEN           ! tropics
          IF((L.GE.LS1).AND.(L.LE.14)) CH4INIT=CH4INIT*   1.620
          IF((L.GE.15).AND.(L.LE.17))  CH4INIT=CH4INIT*   1.460
          IF((L.GE.18).AND.(L.LE.20))  CH4INIT=CH4INIT*   0.812
          IF((L.GE.21).AND.(L.LE.LM))  CH4INIT=CH4INIT*   0.230
        END IF
        do i=1,im
          trm(i,j,l,n_CH4)= am(L,I,J)*CH4INIT
        end do ! i
      end do   ! j
      end do   ! l
      j2=0
      RETURN
      end subroutine get_CH4_IC
