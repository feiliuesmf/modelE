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
      real*8 DMSinput(im,jm,12)
c!@var DMS_AER           DMS prescribed by AERONET (kg S/day/box)
      real*4 DMS_AER(im,jm,366)
c!@var SS1_AER        SALT bin 1 prescribed by AERONET (kg S/day/box)
      real*4 SS1_AER(im,jm,366)
c!@var SS2_AER        SALT bin 2 prescribed by AERONET (kg S/day/box)
      real*4 SS2_AER(im,jm,366)
      INTEGER, PARAMETER :: nso2src  = 1
!@var SO2_src    SO2 industry: surface source (kg/s/box)
      real*8 SO2_src(im,jm,nso2src)
!@var BCI_src    BC Industrial source (kg/s/box)
      real*8 BCI_src(im,jm)
!@var BCB_src    BC Biomass source (kg/s/box)
      real*8 BCB_src(im,jm,lmAER,12)
!@var OCI_src    OC Industrial source (kg/s/box)
      real*8 OCI_src(im,jm)
!@var OCT_src    OC Terpene source (kg/s/box)
      real*8 OCT_src(im,jm,12)
!@var OCB_src    OC Biomass source (kg/s/box)
      real*8 OCB_src(im,jm,lmAER,12)
!@var BCII_src_3D  BCI aircraft source (kg/s)
      real*8 BCI_src_3D(im,jm,lm)
!@var ss_src  Seasalt sources in 2 bins (kg/s/m2)
      INTEGER, PARAMETER :: nsssrc = 2
      real*8 ss_src(im,jm,nsssrc)
      INTEGER, PARAMETER :: nso2src_3D  = 3
!@var SO2_src_3D SO2 volcanic, aircraft sources (and biomass) (kg/s)
      real*8 SO2_src_3D(im,jm,lm,nso2src_3D)
!@var SO2_biosrc_3D SO2  biomass(kg/s)
      real*8 SO2_biosrc_3D(im,jm,lmAER,12)
!@var PBLH boundary layer height
!@var MDF is the mass of the downdraft flux
      real*8, DIMENSION(IM,JM):: PBLH = 0,shdtt = 0.   ! ,MDF
      real*8, DIMENSION(IM,JM,LM):: ohr,dho2r,perjr,
     *   tno3r,oh,dho2,perj,tno3,o3_offline
      real*8, DIMENSION(IM,JM,LM,ntm):: aer_tau
      END MODULE AEROSOL_SOURCES
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
      integer, parameter :: nanns=0,nmons=1
      integer ann_units(nanns),mon_units(nmons),imon(nmons)
      integer i,j,iu,k,l
      integer      :: jdlast=0  
      logical :: ifirst=.true.
      character*80 title   
      character*10 :: mon_files(nmons) = (/'O3_FIELD'/)
      logical      :: mon_bins(nmons)=(/.true./) ! binary file?
      real*8 tlca(im,jm,lm,nmons),tlcb(im,jm,lm,nmons),frac
      REAL*8, DIMENSION(IM,JM,LM,1) :: src
      save jdlast,tlca,tlcb,mon_units,imon,ifirst


C initialise
c      o3_offline(:,:,:)=0.0d0
c
C     Read it in here and interpolated each day.
C
      if (ifirst) call openunits(mon_files,mon_units,mon_bins,nmons)
      ifirst = .false.
      j = 0
      do k=nanns+1,1
        j = j+1
        call read_monthly_O3_3D_source(LM,mon_units(j),jdlast,
     *    tlca(1,1,1,j),tlcb(1,1,1,j),src(1,1,1,k),frac,imon(j))
      end do
      jdlast = jday

      write(6,*) 'Read in ozone offline fields for in cloud oxidation 
     *  interpolated to current day',frac
      call sys_flush(6)

      do k=nanns+1,1; DO l=1,lm; DO J=1,JM; DO I=1,IM
      o3_offline(i,j,l)=src(I,J,L,k)

      end do
      end do
      end do
      end do

      END SUBROUTINE get_O3_offline

      SUBROUTINE read_monthly_O3_3D_source(Ldim,iu,jdlast,tlca,
     * tlcb,data1,frac,imon)
!@sum Read in monthly sources and interpolate to current day
!@ Author Greg Faluvegi 
!@+   Calling routine must have the lines:
!@+      real*8 tlca(im,jm,Ldim,nm),tlcb(im,jm,Ldim,nm)
!@+      integer imon(nm)   ! nm=number of files that will be read
!@+      data jdlast /0/
!@+      save jdlast,tlca,tlcb,imon
!@+   Input: iu, the fileUnit#; jdlast
!@+   Output: interpolated data array + two monthly data arrays
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
      end subroutine read_monthly_O3_3D_source

      subroutine read_SO2_source(nt)
!@sum reads in industrial, biomass, volcanic and aircraft SO2 sources
!@auth Koch
c want kg SO2/m2/s
      USE CONSTANT, only: sday
      USE MODEL_COM, only: im,jm,ls1,jmon,jday,dtsrc,zatmo,t,lm
      USE GEOM, only: dxyp
      USE TRACER_COM
      USE FILEMANAGER, only: openunit,closeunit
      USE AEROSOL_SOURCES, only: SO2_src,SO2_src_3D,nso2src,
     *     nso2src_3D,bci_src_3D,SO2_biosrc_3D
      USE DYNAMICS, only: pmid,pk
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

      logical :: ifirst=.true.
      REAL*8 tx,txx,ty,tyy,txy,tmas,
     *  so2_ind_input(im,jm),cfac,
     *  tb,tbx,tbxx,tby,tbyy,tbxy,sdday,endday,fdry,addtc
     *  ,vemis,zg,vht,zh,pl1,pl2,te,hight,vtot
      real*4 craft(im,jm,lm)
      character*56  titleg
      integer i,j,iuc,ij,nt,
     *  iuc2,ib,jb,iburn,
     *  iv,jv,ivht,ir,l,najl
      save ifirst,so2_ind_input
      if (imPI.eq.1) go to 88
      if (ifirst) then
      call openunit('SO2_IND',iuc,.false.)
      so2_ind_input(:,:)=0.   ! initiallise
      DO 10 ij = 1,9999
      read(iuc,902)   I,J,TMAS !,TX,TXX,TY,TYY,TXY
 902  FORMAT(3X,2(I4),E11.3,5F9.2)
      if (i.eq.0) goto 12
      so2_ind_input(i,j)=tmas
 10    continue
 12   continue
      call closeunit(iuc)
      endif
c I think units are: Kg S/box/yr
c We need kg SO2/m2/s
c     cfac=tr_mm(nt)/32.d0/365.d0/sday  !*dtsrc
c Actually they are kg SO2/box/yr:
      cfac=1.d0/365.d0/sday
      do j=1,jm
      do i=1,im
      so2_src(i,j,1)=so2_ind_input(i,j)*cfac
      end do
      end do

c??      DTT=REAL(NDYN)*DT/(86400.*30.) (ADDTC=TB*B60N*DTT)

c Aircraft emissions
      if (ifirst) then
      so2_src_3D(:,:,:,2)=0.
      bci_src_3D(:,:,:)=0.d0
      call openunit('AIRCRAFT',iuc2,.true.)
C Read 13 layer data file for 23 layer model
      if (lm.eq.23) then
      DO L=1,LS1+1
      READ(iuc2) titleg,((craft(i,j,l),i=1,im),j=1,jm)
      END DO
      else
      DO L=1,LM
      READ(iuc2) titleg,((craft(i,j,l),i=1,im),j=1,jm)
      END DO
      endif
      call closeunit(iuc2)
C Set emissions above LS1+1 to zero for 23 layer model
       if (lm.eq.23) then
       DO L=LS1+2,LM
        DO J=1,JM
         DO I=1,IM
       craft(i,j,l)=0.0
       END DO
       END DO
       END DO
       endif
c craft is Kg fuel/day. Convert to Kg SO2/s. 2.3 factor
c adjusts 2015 source to 1990 source.
c 4.d0d-4 converts to kg S
c (for BC multiply this by 0.1)
      so2_src_3D(:,:,:,2)=craft(:,:,:)*4.0d-4*tr_mm(n_SO2)/32.d0
     *  /2.3d0/sday
      bci_src_3D(:,:,:)=so2_src_3D(:,:,:,2)*0.1d0
      ifirst=.false.
      endif
 88   continue
C  Read in emissions from biomass burning
      call openunit('SO2_BIOMASS',iuc2,.false.)
      so2_biosrc_3D(:,:,:,:)=0.   ! initiallise

      DO 154 ij=1,9999
Cewg  SDDAY -- first day (as day of the year) of the 90 driest days
      READ (iuc2,9051) IB,JB,TB,TBX,TBXX,TBY,TBYY,TBXY,SDDAY
 9051  FORMAT(4X,2I5,6E10.3,F5.0)
      IF (IB.EQ.0) GO TO 155
C Flux is tons SO2 / 4x5 box / year. Convert to kg SO2 / 4x5 box /sec
c Must be tons SO2/4x5 box/month. Convert to kg SO2/4x5 box/sec
c Must be tons S/4x5 box/month. Convert to kg SO2/4x5 box/sec
       TB = TB * 907.2d0/30.4d0/24.d0/3600.d0/32.d0*64.d0
        IF ((SDDAY+89.) .LE. 365.) THEN
          ENDDAY = SDDAY + 89.
        ELSE
          ENDDAY = (SDDAY + 89.) - 365.
        ENDIF
 158  CONTINUE
Cewg Allow burning above 40N only in May, June, July, Aug, and Sept.
      IBURN = 0
      IF (JMON .GE. 5  .AND.  JMON .LE. 9) IBURN = 1
        IF (JB .GE. 34  .AND.  IBURN .EQ. 0) GO TO 154
        IF (JB .GE. 39  .AND.  IBURN .EQ. 1) THEN
          ADDTC = TB*B60N(JMON)
           GO TO 165
        ENDIF
        IF ((JB.GE.34 .AND. JB.LE.38)  .AND.  IBURN .EQ. 1) THEN
          ADDTC = TB*B40N(JMON)
           GO TO 165
        ENDIF
Cewg Allow burning between 32N & 40N from September through April
        IF (JB .GE. 30  .AND.  JB .LE. 33) THEN
          IF (JMON .LE. 4  .OR.  JMON .GE. 9) THEN
             ADDTC = TB*B24N(JMON)
             GO TO 165
          ELSE
            GO TO 154
          ENDIF
        ENDIF
Cewg Allow burning south of 32N on the 90 driest days of the year
        FDRY = 1./3.
        IF (ENDDAY .LT. SDDAY) THEN
          IF (REAL(JDAY).GT.ENDDAY.AND.REAL(JDAY).LT.SDDAY) GOTO 154
        ELSE
          IF (REAL(JDAY).LT.SDDAY.OR.REAL(JDAY).GT.ENDDAY) GOTO 154
        ENDIF
        ADDTC = TB*FDRY
 165    if (imPI.eq.1) ADDTC=0.5d0*ADDTC
       SO2_biosrc_3D(IB,JB,1,jmon) =  ADDTC
 400  CONTINUE
 154  CONTINUE
 155  call closeunit(iuc2)

c continuously erupting volcanic emissions
c     from GEIA (Andres and Kasgnoc, 1998)
      so2_src_3D(:,:,:,1)=0.
      call openunit('SO2_VOLCANO',iuc2,.false.)
      do 116 ir=1,49
      read(iuc2,*) iv,jv,ivht,vemis
C Convert emissions from Mg SO2/day to Kg/sec
       vemis = vemis * 1.e3 /sday
C Find layer to put emissions in
C ZG is height of ground above sea level, in meters
       zg = zatmo(iv,jv)
C VHT is height in meters of volcanic plume above sea level
C Assume all emissions occur at top of plume.
       vht = real(ivht)
       zh = zg
       do 21 l=1,lm-1
        pl1=pmid(l,iv,jv)
        pl2=pmid(l+1,iv,jv)
        te=pk(l,iv,jv)*t(iv,jv,l)
        HIGHT = 2.9271d+01*TE*LOG(PL1/PL2)
C ZH is height in meters of top of layer above sea level
        ZH = ZH + HIGHT
        IF (VHT .LT. ZH) GO TO 24
 21    CONTINUE
       GO TO 100
 24    CONTINUE
       so2_src_3D(iv,jv,l,1)=so2_src_3D(iv,jv,l,1)+vemis
 100   CONTINUE
 116  CONTINUE
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
c     include 'netcdf.inc'
c     integer start(3),count(3),status,ncidu,id1
      integer jread
      REAL*8 akw,erate
c     real*4 DMS_AER
      real*8, INTENT(OUT) :: DMS_flux
      real*8, INTENT(IN) :: swind
      integer, INTENT(IN) :: itype,i,j

      DMS_flux=0
c     if (itype.eq.1) then
        erate=0.d0
        if (imAER.eq.0) then
        if (itype.eq.1) then
c Nightingale et al
        akw = 0.23d0*swind*swind + 0.1d0 * swind
        akw = akw * 0.24d0
        erate=akw*DMSinput(i,j,jmon)*1.d-9*62.d0 !*tr_mm(nt)
     *       /sday
        endif
        else !AEROCOM run, prescribed flux

c         status=NF_OPEN('DMS_SEA',NCNOWRIT,ncidu)
c         status=NF_INQ_VARID(ncidu,'dms',id1)

c         start(1)=i
c         start(2)=j
c         start(3)=jday  

c         count(1)=1
c         count(2)=1
c         count(3)=1

c         status=NF_GET_VARA_REAL(ncidu,id1,start,count,DMS_AER)
c         status=NF_CLOSE('DMS.2000_AEROCOM.nc',NCNOWRIT,ncidu)
c        write(6,*)'DMS RRR',i,j,jday,DMS_AER(i,j,jday)
c if after Feb 28 skip the leapyear day
         jread=jday
         if (jday.gt.59) jread=jday+1
         if (j.eq.1.or.j.eq.46) DMS_AER(i,j,jread)
     *      =DMS_AER(i,j,jread)*72.d0
         erate=DMS_AER(i,j,jread)/sday/dxyp(j)*tr_mm(n_DMS)/32.d0
        endif
        DMS_flux=erate          ! units are kg/m2/s
c     endif
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
       if (imAER.eq.0) then
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
c     endif
      return
      end subroutine read_seasalt_sources


      subroutine aerosol_gas_chem
!@sum aerosol gas phase chemistry
!@auth Dorothy Koch
      USE TRACER_COM
      USE TRACER_DIAG_COM, only : tajls   !,jls_3Dsource,itcon_3Dsrc
     *     ,jls_OHconk,jls_HO2con,jls_NO3,jls_phot
     *     , taijs, ijs_dms_dens,ijs_so2_dens,ijs_so4_dens
      USE MODEL_COM, only: im,jm,jmon,ls1,lm,dtsrc,t,q,jday,
     * coupled_chem
      USE DYNAMICS, only: pmid,am,pk,LTROPO,byam
      USE GEOM, only: dxyp,imaxj,BYDXYP
      USE FLUXES, only: tr3Dsource
      USE FILEMANAGER, only: openunit,closeunit
      USE AEROSOL_SOURCES, only: ohr,dho2r,perjr,tno3r,oh,
     & dho2,perj,
     * tno3
c Aerosol chemistry
      implicit none
      logical :: ifirst=.true.
      real*8 ohx(36,24,lm),dho2x(36,24,lm),perjx(36,24,lm)
      real*8 ppres,te,tt,mm,dmm,ohmc,r1,d1,r2,d2,ttno3,r3,d3,
     * ddno3,dddms,ddno3a,fmom,dtt
      real*8 rk4,ek4,r4,d4
      real*8, DIMENSION(IM,JM,LM):: dms_dens,so2_dens,so4_dens
#ifdef TRACERS_HETCHEM
     *       ,d41,d42,d43,d44
#endif
      real*8 r6,d6,ek9,ek9t,ch2o,eh2o,dho2mc,dho2kg,eeee,xk9,
     * r5,d5,dmssink
      real*8 bciage,ociage
      integer i,j,l,n,iuc,iun,itau,ixx1,ixx2,ichemi,itopen,itt,
     * ittime,isp,iix,jjx,llx,ii,jj,ll,iuc2,it,nm,najl
      save ifirst,itopen,iuc

C**** initialise source arrays
!$OMP PARALLEL DO PRIVATE (L)
      do l=1,lm
        tr3Dsource(:,:,l,1,n_DMS)=0. ! DMS chem sink
        tr3Dsource(:,:,l,1,n_MSA)=0. ! MSA chem sink
        tr3Dsource(:,:,l,4,n_SO2)=0. ! SO2 chem source
        tr3Dsource(:,:,l,5,n_SO2)=0. ! SO2 chem sink
        tr3Dsource(:,:,l,1,n_SO4)=0. ! SO4 chem source
        tr3Dsource(:,:,l,1,n_H2O2_s)=0. ! H2O2 chem source
        tr3Dsource(:,:,l,2,n_H2O2_s)=0. ! H2O2 chem sink
#ifdef TRACERS_HETCHEM
        tr3Dsource(:,:,l,1,n_SO4_d1) =0. ! SO4 on dust
        tr3Dsource(:,:,l,1,n_SO4_d2) =0. ! SO4 on dust
        tr3Dsource(:,:,l,1,n_SO4_d3) =0. ! SO4 on dust
        tr3Dsource(:,:,l,1,n_SO4_d4) =0. ! SO4 on dust
#endif
        if (n_BCII.gt.0) then
          tr3Dsource(:,:,l,1,n_BCII)=0. ! BCII sink
          tr3Dsource(:,:,l,1,n_BCIA)=0. ! BCIA source
        end if
        if (n_OCII.gt.0) then
          tr3Dsource(:,:,l,1,n_OCII)=0. ! OCII sink
          tr3Dsource(:,:,l,1,n_OCIA)=0. ! OCIA source
        end if
      end do
!$OMP END PARALLEL DO

C Coupled mode: use on-line radical concentrations
      if (coupled_chem.eq.1) then
!$OMP PARALLEL DO PRIVATE (L)
        DO L=1,LM
          oh(:,:,l)=oh_live(:,:,l)
          tno3(:,:,l)=no3_live(:,:,l)
c Set h2o2_s =0 and use on-line h2o2 from chemistry
          trm(:,:,l,n_h2o2_s)=0.0
        end do
!$OMP END PARALLEL DO
      endif

      if (coupled_chem.eq.0) then
c Use this for chem inputs from B4360C0M23, from Drew

        call openunit('AER_CHEM',iuc,.true.)
        do ii=1,jmon
          read(iuc) ichemi
          read(iuc) ohr
          read(iuc) dho2r
          read(iuc) perjr
          read(iuc) tno3r
        end do
        call closeunit(iuc)

c need to scale TNO3, OH and PERJ using cosine of zenith angle
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
#endif
      dtt=dtsrc
C**** THIS LOOP SHOULD BE PARALLELISED
      do 20 l=1,LM
      do 21 j=1,jm
      do 22 i=1,imaxj(j)
c
C Initialise       
        dms_dens(i,j,l)=0.0D0
        so2_dens(i,j,l)=0.0D0
        so4_dens(i,j,l)=0.0D0

      if(l.le.ltropo(i,j)) then

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
c       bciage=4.3D-6*trm(i,j,l,n) !used this first        
c       bciage=1.0D-6*trm(i,j,l,n)  !2nd
        bciage=1.0D-7*trm(i,j,l,n)
        tr3Dsource(i,j,l,1,n)=-bciage        
        tr3Dsource(i,j,l,1,n_BCIA)=bciage        

        case ('OCII')
c       ociage=7.3D-6*trm(i,j,l,n)       !used this first 
c       ociage=3.6D-6*trm(i,j,l,n)     !2nd
        ociage=3.D-7*trm(i,j,l,n)
        tr3Dsource(i,j,l,1,n)=-ociage        
        tr3Dsource(i,j,l,1,n_OCIA)=ociage        

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
        end select
        
 23   CONTINUE
        endif
 22   CONTINUE
 21   CONTINUE
 20   CONTINUE

cg this is now done automatically in apply_tracer_3Dsource
cg       call DIAGTCA(itcon_3Dsrc(3,n_SO2),n_SO2)


      do 30 l=1,lm
      do 31 j=1,jm
      do 32 i=1,imaxj(j)

      if(l.le.ltropo(i,j)) then

      ppres=pmid(l,i,j)*9.869d-4 !in atm
      te=pk(l,i,j)*t(i,j,l)
      mm=am(l,i,j)*dxyp(j)
      tt = 1.d0/te
      dmm=ppres/(.082d0*te)*6.02d20
      ohmc = oh(i,j,l)          !oh is alread in units of molecules/cm3

      do 33 n=1,ntm

        select case (trname(n))


        case ('SO2')
c oxidation of SO2 to make SO4: SO2 + OH -> H2SO4

          r4=rsulf4(i,j,l)*ohmc
          d4 = exp(-r4*dtsrc)

c     IF (I.EQ.30.AND.J.EQ.30.and.L.EQ.2) WRITE(6,*)'msulf',TE,DMM,
c     *  PPRES,RK4,EK4,R4,D4,ohmc
          IF (d4.GE.1.) d4=0.99999d0
#ifdef TRACERS_HETCHEM
       d41 = exp(-rxts1(i,j,l)*dtsrc)     
       d42 = exp(-rxts2(i,j,l)*dtsrc)     
       d43 = exp(-rxts3(i,j,l)*dtsrc)     
       d44 = exp(-rxts4(i,j,l)*dtsrc)     
       tr3Dsource(i,j,l,5,n) = (-trm(i,j,l,n)*(1.d0-d41)/dtsrc)
     .                       + ( -trm(i,j,l,n)*(1.d0-d4)/dtsrc)
     .                       + ( -trm(i,j,l,n)*(1.d0-d42)/dtsrc)
     .                       + ( -trm(i,j,l,n)*(1.d0-d43)/dtsrc)
     .                       + ( -trm(i,j,l,n)*(1.d0-d44)/dtsrc)
 
#else
       tr3Dsource(i,j,l,5,n) = -trm(i,j,l,n)*(1.d0-d4)/dtsrc 
#endif         
          
c diagnostics to save oxidant fields
          najl = jls_OHconk
          tajls(j,l,najl) = tajls(j,l,najl)+oh(i,j,l)
          najl = jls_HO2con
          tajls(j,l,najl) = tajls(j,l,najl)+dho2(i,j,l)

#ifdef TRACERS_HETCHEM
       case ('SO4_d1')
c sulfate production from SO2 on mineral dust aerosol

       tr3Dsource(i,j,l,1,n) = tr3Dsource(i,j,l,1,n) +  tr_mm(n)/
     *             tr_mm(n_so2)*(1.d0-d41)*trm(i,j,l,n_so2)/dtsrc

       case ('SO4_d2')
c sulfate production from SO2 on mineral dust aerosol

       tr3Dsource(i,j,l,1,n) = tr3Dsource(i,j,l,1,n) +  tr_mm(n)/
     *             tr_mm(n_so2)*(1.d0-d42)*trm(i,j,l,n_so2)/dtsrc

       case ('SO4_d3')
c sulfate production from SO2 on mineral dust aerosol

       tr3Dsource(i,j,l,1,n) = tr3Dsource(i,j,l,1,n) +  tr_mm(n)/
     *             tr_mm(n_so2)*(1.d0-d43)*trm(i,j,l,n_so2)/dtsrc

       case ('SO4_d4')
c sulfate production from SO2 on mineral dust aerosol

       tr3Dsource(i,j,l,1,n) = tr3Dsource(i,j,l,1,n) +  tr_mm(n)/
     *             tr_mm(n_so2)*(1.d0-d44)*trm(i,j,l,n_so2)/dtsrc

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

C Calculate diagnostics for chemistry
C Number density (molecule/cm3) of DMS,SO2,SO4 
          dms_dens(i,j,l)=trm(I,J,L,n_DMS)*dmm*(28.0D0/62.0D0)*
     *   BYDXYP(J)*BYAM(L,I,J)
          so2_dens(i,j,l)=trm(I,J,L,n_SO2)*dmm*(28.0D0/64.0D0)*
     *   BYDXYP(J)*BYAM(L,I,J)
          so4_dens(i,j,l)=trm(I,J,L,n_SO4)*dmm*(28.0D0/96.0D0)*
     *   BYDXYP(J)*BYAM(L,I,J)

        endif           ! end troposphere

C Accumulate diagnostics for chemistry
          taijs(i,j,ijs_dms_dens(l))=taijs(i,j,ijs_dms_dens(l))
     * + dms_dens(i,j,l)

          taijs(i,j,ijs_so2_dens(l))=taijs(i,j,ijs_so2_dens(l))
     * + so2_dens(i,j,l)

          taijs(i,j,ijs_so4_dens(l))=taijs(i,j,ijs_so4_dens(l))
     * + so4_dens(i,j,l)

 32   CONTINUE
 31   CONTINUE
 30   CONTINUE

      RETURN
      END subroutine aerosol_gas_chem


      SUBROUTINE SCALERAD
      use MODEL_COM, only: im,jm,lm,jday,jhour,jmon,itime,nday,dtsrc
      use AEROSOL_SOURCES, only: ohr,dho2r,perjr,tno3r,oh,dho2,perj,tno3
      use CONSTANT, only: radian,teeny
c     use RADNCB, only: cosz1
      implicit none
      real*8 ang1,xnair,vlon,vlat,ctime,timec,p1,p2,p3,fact,rad,
     *  rad1,rad2,rad3,stfac
      real*8 suncos(im,jm),tczen(im,jm)
      integer i,j,hrstrt,jdstrt,ihr,l
      integer nradn(im,jm)

      ang1 = 90.d0/91.3125d0
      xnair = 6.022d23 * 1.d3 / 28.9644d0
      jdstrt=0
      hrstrt=0
C*
      DO 100 j = 1,jm
      DO 100 i = 1,im
C*** Calculate cos theta (RAD) at the beginning of the time step (RAD)
      vlon = 180.d0 - (i-1)*5.d0
      vlat=-90.d0+(j-1)*4.d0
      if (j.eq.1) vlat=-88.d0
      if (j.eq.46) vlat=88.d0
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
        do j = 1, jm
           do i = 1, im
             tczen(i,j) = 0.
             nradn(i,j)=0
             vlon = 180.d0 - (i-1)*5.d0
             vlat=-90.d0+(j-1)*4.d0
             if (j.eq.1) vlat=-88.d0
             if (j.eq.46) vlat=88.d0
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
      do l = 1,lm
          do j = 1,jm
            do i = 1,im
c Get NO3 only if dark, weighted by number of dark hours
c        if (I.EQ.1.AND.L.EQ.1) write(6,*)'NO3R',TAU,J,NRADN(I,J)
            if (suncos(i,j).eq.0.and.nradn(i,j).gt.0) then
            tno3(i,j,l)=tno3r(i,j,l)*24.d0/real(nradn(i,j))
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
c                write (6,*) 'SCALERAD problem',jday,jhour,i,j,
c    *               suncos(i,j),tczen(i,j)
                 stfac=1.
                endif
                 oh(i,j,l)=ohr(i,j,l)*6.d0*stfac
                 perj(i,j,l)=perjr(i,j,l)*6.d0*stfac
c     if (i.eq.30.and.j.eq.40.and.l.eq.2) write(6,*) 'ohrad2',
c    *  suncos(i,j),tczen(i,j)
                 dho2(i,j,l)=dho2r(i,j,l)*6.d0*stfac
              end if
c     if ((i.eq.11.or.i.eq.13.or.i.eq.15).and.j.eq.4.) then
c      if (l.eq.2)
c    *  write(6,*)
c    *  'HO2 SCALE',jday,jhour,
c    * i,l,dho2r(i,j,l),dho2(i,j,l),suncos(i,j),cosz1(i,j),tczen(i,j)
c     endif
           end do
         end do
      end do
      RETURN
      END subroutine SCALERAD


      SUBROUTINE GET_SULFATE(L,temp,fcloud,
     *  wa_vol,wmxtr,sulfin,sulfinc,sulfout,tr_left,
     *  tm,tmcl,airm,LHX,dt_sulf)
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
      real*8 tm(lm,ntm), tmcl(ntm), airm(lm)
c     REAL*8,  INTENT(OUT)::
      real*8 sulfin(ntm),sulfout(ntm),tr_left(ntm)
     *  ,sulfinc(ntm)
!@var dt_sulf accumulated diagnostic of sulfate chemistry changes
      real*8, dimension(ntm), intent(inout) :: dt_sulf
      INTEGER, INTENT(IN) :: L
      do n=1,ntx
        sulfin(N)=0.
        sulfinc(N)=0.
        sulfout(N)=0.
        tr_left(N)=1.
      end do

c
C**** CALCULATE the fraction of tracer mass that becomes condensate:
c
      if (LHX.NE.LHE.or.fcloud.lt.teeny) go to 333
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
c dmk should wa_vol be incremental water volume??
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
       case('SO4')
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

      SUBROUTINE HETER
!@sum heterogeneous production of sulfate
!@auth Dorothy Koch
      USE CONSTANT, only: BYGASC, MAIR,teeny,mb2kg,gasc,pi
      USE TRACER_COM, only:tr_RKD,tr_DHD,n_H2O2_s,n_SO2
     *     ,trname,ntm,tr_mm,n_SO4,trm,trmom,n_H2O2
      USE MODEL_COM, only: im,jm,lm,dtsrc,t,p,coupled_chem
      USE CLOUDS_COM, only:rhsav
      USE DYNAMICS, only: pmid,am,pk,ltropo
      USE GEOM, only: dxyp,imaxj
      USE FLUXES, only : tr3Dsource

      IMPLICIT NONE
      integer i,j,l,n,najl
      REAL*8, PARAMETER :: BY298K=3.3557D-3
      real*8, parameter :: rk1=1.3d-2 !M
      real*8, parameter :: dh1=-1.6736d4 !J/mol
      real*8, parameter :: rk=6.357d14    !1/(M*M*s)
      real*8, parameter :: ea=3.95d4 !J/mol
      real*8 te,amass,ppas,rk1f,tfac,ssfac,RKD,Henry_const(ntm),
     * pph(ntm),r1,A,B,pp,rr,aa,bb,BA,B2,xx,y,pn,tv,avol,
     * dso4g,dso4gt,tso2,th2o2,sulfout,sulfin,wv,ss,clwc,
     * rkdm(ntm),trmol(ntm),tt1,tt2,tt3,ptr


      ptr=0.3d0  !um
      DO 19 L=1,LM
      DO 21 J=1,JM
      DO 20 I=1,IMAXJ(J)

C**** initialise source arrays
        tr3Dsource(i,j,l,2,n_SO4)=0.
        tr3Dsource(i,j,l,6,n_SO2)=0.
        tr3Dsource(i,j,l,3,n_H2O2_s)=0.

      if(l.le.ltropo(i,j)) then

      amass=am(l,i,j)*DXYP(j)   !kg
      ppas=pmid(l,i,j)*100.    !Pa
      te=pk(l,i,j)*t(i,j,l)
      tfac = (1./te - by298k)*bygasc  !mol/J
c  water from deliquescence; assume complete dissociation of (NH)2SO4
c   Assume average radius of 0.1 um
      ss=rhsav(l,i,j)
      if (ss.gt.0.99d0.or.ss.le.0.0d0) go to 88
         R1=log(SS)
         A=3.2d-7/te*1d6
c        B=1.09d-4*7.4 !7.4 is particle mass for 0.1 um particle
         B=1.09d-4*4.d0*pi/3.d0*1760.d0*ptr*ptr*ptr
         pp=-A/R1
         rr=B/R1
         aa=-pp*pp/3.
         bb=(2.*pp*pp*pp+27.*rr)/27.
         BA=(-bb/2.+sqrt(bb*bb/4.+aa*aa*aa/27.))**(1./3.)
         B2=(-(bb/2.+sqrt(bb*bb/4.+aa*aa*aa/27.)))**(1./3.)
         xx=BA+B2
         y=xx-pp/3.
c y is wet diameter. wv is volume of water per particle
       wv=4./3.*pi*(y-ptr)**3.d0*1000.d0/1.E18
c pn is number of particles
       pn=trm(i,j,l,n_so4)/1760.d0*3.d0/(4.d0*pi)*1.E18/(ptr)**3
c tv is total volume
       tv=wv*pn
c clwc is vol water/volume air
c      clwc=wmxtr*mair*ppas/temp*bygasc/1.D6/fcloud
c avol is air volume
      avol=amass/mair*1000.d0*gasc*te/ppas*1000.d0   !L
      clwc=tv/avol
c
      rk1f=rk1*exp(-dh1*tfac)


      if (coupled_chem.eq.1) then
      if (trm(i,j,l,n_so2).lt.teeny.or.
     *  trm(i,j,l,n_h2o2).lt.teeny.or.tv.lt.teeny) then
      go to 88
      endif
      endif

      if (coupled_chem.eq.0) then
      if (trm(i,j,l,n_so2).lt.teeny.or.
     *  trm(i,j,l,n_h2o2_s).lt.teeny.or.tv.lt.teeny) then
      go to 88
      endif
      endif

      DO N=1,NTM
       select case (trname(n))
       case('SO2')
      tso2=trm(i,j,l,n)

c modified Henry's Law coefficient assuming pH of 4.5
      rkdm(n)=tr_rkd(n)*(1.+ rk1f/3.2d-5)
c mole of tracer, used to limit so4 production
      trmol(n)=1000.*trm(i,j,l,n)/tr_mm(n)
c partial pressure of gas x henry's law coefficient
      pph(n)=mair*1.d-3*ppas/tr_mm(n)/amass*
     *   tr_rkd(n)*exp(-tr_dhd(n)*tfac)
c the following is from Phil:
c      reduction in partial pressure as species dissolves
      henry_const(n)=rkdm(n)*exp(-tr_dhd(n)*tfac)
      pph(n)=pph(n)/(1+(henry_const(n)*clwc*gasc*te))

       case('H2O2','H2O2_s')

        if (trname(n).eq."H2O2" .and. coupled_chem.eq.0) goto 402
        if (trname(n).eq."H2O2_s" .and. coupled_chem.eq.1) goto 402

      th2o2=trm(i,j,l,n)
c modified Henry's Law coefficient assuming pH of 4.5
      rkdm(n)=tr_rkd(n)*(1.+ rk1f/3.2d-5)
c mole of tracer, used to limit so4 production
      trmol(n)=1000.*trm(i,j,l,n)/tr_mm(n)
c partial pressure of gas x henry's law coefficient
      pph(n)=mair*1.d-3*ppas/tr_mm(n)/amass*
     *   tr_rkd(n)*exp(-tr_dhd(n)*tfac)
c the following is from Phil:
c      reduction in partial pressure as species dissolves
      henry_const(n)=rkdm(n)*exp(-tr_dhd(n)*tfac)
      pph(n)=pph(n)/(1+(henry_const(n)*clwc*gasc*te))

 402  CONTINUE

      end select
      END DO

c this part from gas phase:moles/kg/kg
        if(coupled_chem.eq.1) then
      dso4g=rk*exp(-ea/(gasc*te))*rk1f
     *    *pph(n_h2o2)*pph(n_so2)*dtsrc*tv
         endif
        if(coupled_chem.eq.0) then
      dso4g=rk*exp(-ea/(gasc*te))*rk1f
     *    *pph(n_h2o2_s)*pph(n_so2)*dtsrc*tv
         endif

c check to make sure no overreaction: moles of production:
      dso4gt=dso4g*tso2*th2o2
c can't be more than moles going in:
      if (dso4gt.gt.trmol(n_so2)) then
        dso4g=trmol(n_so2)/(tso2*th2o2)
      endif
      dso4gt=dso4g*tso2*th2o2

       if (coupled_chem.eq.1) then
      if (dso4gt.gt.trmol(n_h2o2)) then
        dso4g=trmol(n_h2o2)/(tso2*th2o2)
      endif
       endif
       if (coupled_chem.eq.0) then
      if (dso4gt.gt.trmol(n_h2o2_s)) then
        dso4g=trmol(n_h2o2_s)/(tso2*th2o2)
      endif
       endif

      do n=1,ntm
       select case (trname(n))
       case('SO4')
       sulfin=0.
       sulfout=tr_mm(n)/1000.*(dso4g*tso2*th2o2) !kg
c diagnostic
       tr3Dsource(i,j,l,2,n)=sulfout/dtsrc

       sulfin=0.
       case('SO2')
       sulfin=-dso4g*th2o2*tr_mm(n)/1000. !dimnless
       sulfin=max(-1d0,sulfin)

       tr3Dsource(i,j,l,6,n)=trm(i,j,l,n)*sulfin/dtsrc

       case('H2O2','H2O2_s')
        if (trname(n).eq."H2O2" .and. coupled_chem.eq.0) goto 403
        if (trname(n).eq."H2O2_s" .and. coupled_chem.eq.1) goto 403

       sulfin=-dso4g*tso2*tr_mm(n)/1000. !dimnless
       sulfin=max(-1d0,sulfin)

       tr3Dsource(i,j,l,3,n)=trm(i,j,l,n)*sulfin/dtsrc

 403   CONTINUE

      end select
      END DO


  88  continue
       endif

  20  CONTINUE
  21  CONTINUE
  19  CONTINUE

      return
      END subroutine HETER

      SUBROUTINE GET_TAU
!@sum  calculate aerosol optical thickness
!@auth Dorothy Koch
      USE TRACER_COM, only:
     *     trname,ntm,trm
      USE MODEL_COM, only: im,jm,lm,idacc
      USE CLOUDS_COM, only:rhsav
      USE GEOM, only: imaxj,dxyp
      USE AEROSOL_SOURCES, only: aer_tau
      USE TRACER_DIAG_COM, only : taijs,ijts_tau,ia_ijts

      integer i,j,l,n,naijs
      real*8 rhh,tf
c For now we are using the hydrated extinction efficiencies
c   from Chin et al 2002.
c No need to divide by dxyp here because it is done in TRACER_PRT
      DO L=1,LM
      DO J=1,JM
      DO I=1,IMAXJ(J)
      rhh=rhsav(l,i,j)*100.d0
      DO N=1,NTM
       select case (trname(n))
       case('SO4')
       if (rhh.le.10.) then
       tf=5.
       else if (rhh.le.20.) then
       tf=6.
       else if (rhh.le.30.) then
       tf=7.
       else if (rhh.le.40.) then
       tf=8.
       else if (rhh.le.50.) then
       tf=10.
       else if (rhh.le.60.) then
       tf=12.
       else if (rhh.le.70.) then
       tf=14.
       else if (rhh.le.82.) then
       tf=17.
       else if (rhh.le.87.) then
       tf=18.
       else if (rhh.le.93.) then
       tf=20.
       else if (rhh.le.97.) then
       tf=26.
       else
       tf=36.
       endif
       aer_tau(i,j,l,n)=trm(i,j,l,n)*1000.*tf/dxyp(j)

       case ('seasalt1')
       if (rhh.le.10.) then
       tf=1.
       else if (rhh.le.20.) then
       tf=2.
       else if (rhh.le.30.) then
       tf=2.5
       else if (rhh.le.40.) then
       tf=3.
       else if (rhh.le.50.) then
       tf=3.3
       else if (rhh.le.60.) then
       tf=3.6
       else if (rhh.le.70.) then
       tf=4.
       else if (rhh.le.82.) then
       tf=4.5
       else if (rhh.le.87.) then
       tf=5.
       else if (rhh.le.93.) then
       tf=6.
       else if (rhh.le.97.) then
       tf=8.
       else
       tf=20.
       endif
       aer_tau(i,j,l,n)=trm(i,j,l,n)*1000.*tf/dxyp(j)

       case ('seasalt2')
       if (rhh.le.10.) then
       tf=0.1
       else if (rhh.le.30.) then
       tf=.2
       else if (rhh.le.40.) then
       tf=.3
       else if (rhh.le.50.) then
       tf=.35
       else if (rhh.le.60.) then
       tf=.4
       else if (rhh.le.70.) then
       tf=.45
       else if (rhh.le.82.) then
       tf=.5
       else if (rhh.le.87.) then
       tf=.6
       else if (rhh.le.93.) then
       tf=.7
       else if (rhh.le.97.) then
       tf=1.
       else
       tf=2.5
       endif
       aer_tau(i,j,l,n)=trm(i,j,l,n)*1000.*tf/dxyp(j)

       end select
C**** diagnostics
       if (ijts_tau(1,n).gt.0) then
         naijs=ijts_tau(1,n)
         taijs(i,j,naijs)=taijs(i,j,naijs)+aer_tau(i,j,l,n)
       end if
      END DO
      END DO
      END DO
      END DO

c     write(6,*) 'AER_TAU',jhour,
c    * ijts_tau(4),idacc(ia_ijts(ijts_tau(4))),idacc(1)
c     write(6,*) aer_tau(54,14,1,4),aer_tau(20,20,2,6)
c     write(6,*) taijs(20,20,ijts_tau(4))
      return
      END subroutine GET_TAU

