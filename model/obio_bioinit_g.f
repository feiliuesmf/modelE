#include "rundeck_opts.h"

#ifdef OBIO_ON_GARYocean
      subroutine obio_bioinit_g
!note: 
!obio_bioinit is called only for a cold start and reads in 
!INITIAL conditions and interpolates them
!to the ocean grid. Such fields are nitrates,silicate,dic
!obio_init  is called for every start of the run and reads 
!in BOUNDARY conditions and interpolates 
!them to ocean grid. such fields are iron,alkalinity,chlorophyl

c based on /g6/aromanou/Watson_new/BioInit/rstbio.F

c  Makes initialization data files for biological variables.
c  This is for the global model.  To subset, use the routines in
c  /u2/gregg/bio/biodat/subreg.
c  Uses NOAA 2001 atlas for NO3 and SiO2 distributions.
c  Includes initial iron distributions.
c
c  Particle type 1  = nitrate
c  Particle type 2  = ammonium
c  Particle type 3  = silicate
c  Particle type 4  = iron
c  Particle type 5  = diatoms
c  Particle type 6  = chlorophytes
c  Particle type 7  = cyanobacteria
c  Particle type 8  = coccolithophores
c  Particle type 9  = dinoflagellates
c  Particle type 10 = zooplankton
c
c  Detritus type 1  = carbon/nitrogen
c  Detritus type 2  = silica
c  Detritus type 3  = iron
c
c  Carbon type 1    = semi-labile DOC
c  Carbon type 2    = DIC
 
      USE FILEMANAGER, only: openunit,closeunit

      USE obio_dim
      USE obio_incom
      USE obio_forc, only: avgq
      USE obio_com, only: gcmax,tracer_loc,tracer
#ifdef TRACERS_Alkalinity
      USE obio_forc, only: alk
#endif

      USE OCEANRES, only : idm=>imo,jdm=>jmo,kdm=>lmo,dzo
      USE OCEAN, only : LMOM=>LMM,ZOE=>ZE,focean,hocean
      USE OCEANR_DIM, only : ogrid

      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT,pack_data,unpack_data
 
      implicit none

      integer i,j,k,l,nn
      integer i1,i2,j1,j2,k1,k2
      integer lm

      integer nir,nt
      integer iu_bioinit

      real dicmin,dicmax
      real zz

      real fldo2(idm,jdm,kdm)
      real fldoz(idm,jdm,kdm)

      character*2 ntchar
      character*80 filename

      common /bir/ nir(nrg)


      INTEGER :: j_0h,j_1h,i_0h,i_1h

      I_0H = ogrid%I_STRT_HALO
      I_1H = ogrid%I_STOP_HALO
      J_0H = ogrid%J_STRT_HALO
      J_1H = ogrid%J_STOP_HALO


      call alloc_obio_incom

      !call gather_tracer
      call pack_data( ogrid,  tracer_loc, tracer )


      if (AM_I_ROOT()) then
c  Initialize

      tracer(:,:,:,1:ntyp)=0.

      Fer(:,:,:) = 0.d0


      filename='nitrates_inicond'
      call bio_inicond_g(filename,fldo2,fldoz)
      tracer(:,:,:,1)=fldo2         !because 3 is nintrate

      filename='silicate_inicond'
      call bio_inicond_g(filename,fldo2,fldoz)
      tracer(:,:,:,3)=fldo2         !because 3 is silicate

      filename='dic_inicond'
      call bio_inicond_g(filename,fldo2,fldoz)
      dic(:,:,:)=fldo2

#ifdef TRACERS_Alkalinity
      filename='alk_inicond'
      call bio_inicond_g(filename,fldo2,fldoz)
      alk(:,:,:)=fldo2
     
      !remove negative values
      !negs are over land or under ice due to GLODAP missing values in the Arctic Ocean
      !for under ice missing values, use climatological minimums for sets of layers based
      !on GLODAP, rather than setting to the same global min.
      do j=1,jdm
      do i=1,idm
      do k=1,kdm
       if (alk(i,j,k).lt.0. .and. focean(i,j).gt.0) then
          if (zoe(k).le.150.) alk(i,j,k)=2172.      !init neg might be under ice,
          if (zoe(k).gt.150. .and. zoe(k).lt.1200.) alk(i,j,k)=2200. 
          if (zoe(k).ge.1200.) alk(i,j,k)=2300.      
       endif
       if (focean(i,j).le.0) then
         alk(i,j,k)=0.
       endif
      enddo
      enddo
      enddo
      tracer(:,:,:,ntrac)=alk       !because ntrac is alkalinity
#endif

!!these rno3 and so2_init fields are not correct. There are void points due to
!mismatch of the noaa grid and the ocean grid. To fill in have to do
!the interpolation in matlab (furtuna).
!at the same time use dps and interpolate to layer depths from the model

      dicmin = 1.d30
      dicmax =-1.d30

      do k=1,kdm
      do j=1,jdm
       do i=1,idm
         if(tracer(i,j,k,1).le.0.d0)tracer(i,j,k,1)=0.085d0
         if(tracer(i,j,k,3).le.0.d0)tracer(i,j,k,3)=0.297d0
          if (dic(i,j,k).le.0.d0) dic(i,j,k)=1837d0
          dic(i,j,k)=dmax1(dic(i,j,k),1837d0)   !set minimum =1837
          dicmin =dmin1(dicmin,dic(i,j,k))
          dicmax =dmax1(dicmax,dic(i,j,k))
        enddo
       enddo
      enddo
      write(*,'(a,2e12.4)')'BIO: bioinit: dic min-max=',
     .       dicmin,dicmax

c  Obtain region indicators
      write(6,*)'calling fndreg...'
      call fndreg
 
c  Define Fe:NO3 ratios by region, according to Fung et al. (2000)
c  GBC.  Conversion produces nM Fe, since NO3 is as uM

      do j=1,jdm
       do i=1,idm
        do k=1,kdm
          if (ir(i,j) .eq. 1)then        !Antarctic
           Fer(i,j,k) = 3.5E-3
          else if (ir(i,j) .eq. 2)then   !South Indian
           Fer(i,j,k) = 20.0E-3
          else if (ir(i,j) .eq. 3)then   !South Pacific
           Fer(i,j,k) = 4.5E-3
          else if (ir(i,j) .eq. 4)then   !South Atlantic
           Fer(i,j,k) = 20.0E-3
          else if (ir(i,j) .eq. 5)then   !Equatorial Indian
           Fer(i,j,k) = 22.5E-3
          else if (ir(i,j) .eq. 6)then   !Equatorial Pacific
           Fer(i,j,k) = 4.0E-3
          else if (ir(i,j) .eq. 7)then   !Equatorial Atlantic
           Fer(i,j,k) = 20.0E-3
          else if (ir(i,j) .eq. 8)then   !North Indian
           Fer(i,j,k) = 25.0E-3
          else if (ir(i,j) .eq. 9)then   !North Central Pacific
           Fer(i,j,k) = 6.0E-3
          else if (ir(i,j) .eq. 10)then  !North Central Atlantic
           Fer(i,j,k) = 20.0E-3
          else if (ir(i,j) .eq. 11)then  !North Pacific
           Fer(i,j,k) = 6.0E-3
          else if (ir(i,j) .eq. 12)then  !North Atlantic
           Fer(i,j,k) = 10.0E-3
          else if (ir(i,j) .eq. 13)then  !Mediterranean
           Fer(i,j,k) = 20.0E-3
          endif
        enddo  
       enddo  
      enddo  
 
c  Create arrays 
      write(6,*)'Creating bio restart data for ',ntyp,' arrays and'
     . ,kdm,'  layers...'

! the following loop has parallelization problems....

      do 1000 j=1,jdm
      do 1000 i=1,idm
      do 1000 k=1,kdm


          !Nitrate
          !read earlier from file

          !Ammonium
          tracer(i,j,k,2) = 0.5
          if (ZOE(k) .gt. 4000.d0) tracer(i,j,k,2) = Pdeep(2)

          !Silica
          !read earlier from file

          !Iron
          tracer(i,j,k,4) = Fer(i,j,k)*tracer(i,j,k,1)  !Fung et al. 2000
c          if (ir(nw) .eq. 3)then
c           P(i,j,k,4) = 0.04*float(k-1) + 0.2
c           P(i,j,k,4) = 0.06*float(k-1) + 0.2
c           P(i,j,k,4) = 0.08*float(k-1) + 0.2
c           P(i,j,k,4) = min(P(i,j,k,4),0.65)
c           P(i,j,k,4) = min(P(i,j,k,4),0.75)
c          endif
          if (ir(i,j) .eq. 1)then
           tracer(i,j,k,4) = Fer(i,j,k)*0.5*tracer(i,j,k,1)
          endif
          tracer(i,j,k,4) = max(tracer(i,j,k,4),0.01)
          if (ZOE(k) .gt. 4000.d0) tracer(i,j,k,4) = Pdeep(4)

          write(*,'(a,3i5,e12.4,i5,e12.4)')'obio_bioinit, iron:',
     .       i,j,k,Fer(i,j,k),ir(i,j),tracer(i,j,k,4)

          !Herbivores
          do nt = nnut+1,ntyp-nzoo
           tracer(i,j,k,nt) = 0.05
          enddo
          do nt = ntyp-nzoo+1,ntyp
           tracer(i,j,k,nt) = 0.05  !in chl units mg/m3
c          tracer(i,j,k,nt) = 0.05*50.0  !in C units mg/m3
          enddo

          !intert tracer
          do nt = ntyp+1,ntyp+n_inert
           tracer(i,j,k,nt) = tracer(i,j,k,1)
          enddo

          !DIC
          !read earlier from file
          dicmod(i,j,k)=dic(i,j,k)

 1000 continue


c  Detritus (set to 0 for start up)
      write(6,*)'Detritus...'
      cnratio = 106.0/16.0*12.0    !C:N ratio (ugl:uM)
      csratio = 106.0/16.0*12.0    !C:Si ratio (ugl:uM)
      cfratio = 150000.0*12.0*1.0E-3    !C:Fe ratio (ugl:nM)

      do j=1,jdm
      do i=1,idm
      do k=1,kdm
           !only detritus components
           tracer(i,j,k,ntyp+n_inert+1) = tracer(i,j,k,1)*0.25*cnratio !as carbon
           tracer(i,j,k,ntyp+n_inert+2) = tracer(i,j,k,3)*0.1
           tracer(i,j,k,ntyp+n_inert+3) = tracer(i,j,k,4)*0.25
           tracer(i,j,k,ntyp+n_inert+1) = 0.0
           tracer(i,j,k,ntyp+n_inert+2) = 0.0
           tracer(i,j,k,ntyp+n_inert+3) = 0.0
      enddo
      enddo
      enddo

c  Carbon (set to 0 for start up)
c   DIC is derived from GLODAP.  Using mean H from exp601,
c   mean DIC for these values is computed.  Surface DIC is taken
c   as the mean for 020m deeper than the mixed layer, converted from
c   uM/kg to uM
      write(6,*)'Carbon...'
c    conversion from uM to mg/m3

      do j=1,jdm
      do i=1,idm
      do k=1,kdm
          tracer(i,j,k,ntyp+n_inert+ndet+1) = 0.0
          tracer(i,j,k,ntyp+n_inert+ndet+2) = 0.0
      enddo
      enddo
      enddo

      !only carbon components
      do j=1,jdm
      do i=1,idm
      do k=1,kdm
          tracer(i,j,k,ntyp+n_inert+ndet+2) = dicmod(i,j,k)
c         car(i,j,k,1) = 3.0  !from Bissett et al 1999 (uM(C))
c         car(i,j,k,1) = 0.0  !from Walsh et al 1999
      enddo
      enddo
      enddo


      endif   !if am_i_root

      !call scatter_tracer
      call unpack_data( ogrid,  tracer, tracer_loc )

c  Light saturation data
      write(6,*)'Light saturation data...'
      avgq = 0.0d0

      do j=j_0h,j_1h
      do i=i_0h,i_1h
      do k=1,kdm
          avgq(i,j,k) = 25.0
      enddo
      enddo
      enddo

c  Coccolithophore max growth rate
      write(6,*)'Coccolithophore max growth rate...'
      do j=j_0h,j_1h
      do i=i_0h,i_1h
      do k=1,kdm
          gcmax(i,j,k) = 0.0
      enddo
      enddo
      enddo
 
!     !save initialization
!     do nt=1,ntyp+n_inert+ndet+ncar
!       ntchar='00'
!       if(nt.le.9)write(ntchar,'(i1)')nt
!       if(nt.gt.9)write(ntchar,'(i2)')nt
!       print*,'BIO: saving initial tracer fields '
!    .        ,'bioinit_tracer'//ntchar
!       call openunit('bioinit_tracer'//ntchar,iu_bioinit)
!     do k=1,kdm
!     do j=1,jdm
!     do i=1,idm
!          write(iu_bioinit,'(3i4,4e12.4)')
!    .           i,j,k,dzo(k),tracer(i,j,k,nt),focean(i,j),hocean(i,j)
!       enddo
!       enddo
!       enddo
!     call closeunit(iu_bioinit)
!     enddo

      print*,'bioinit: COLD INITIALIZATION'
      call obio_trint

      return
      end subroutine obio_bioinit_g

c------------------------------------------------------------------------------
      subroutine fndreg
 
c  Finds nwater indices corresponding to significant regions,
c  and defines arrays.  Variables representative of the regions will 
c  later be kept in these array indicators.
c  Regions are defined as follows:
c        1 -- Antarctic
c        2 -- South Indian
c        3 -- South Pacific
c        4 -- South Atlantic
c        5 -- Equatorial Indian
c        6 -- Equatorial Pacific
c        7 -- Equatorial Atlantic
c        8 -- North Indian
c        9 -- North Central Pacific
c       10 -- North Central Atlantic
c       11 -- North Pacific
c       12 -- North Atlantic
c       13 -- Mediterranean/Black Seas
 

      USE OCEANR_DIM, only : ogrid   
      USE OCEANRES, only : idm=>imo,jdm=>jmo,kdm=>lmo
      Use GEOM,       only : LON_DG,LAT_DG
      Use OCEAN,      only : LMOM=>LMM, ZOE=>ZE

      USE obio_dim
      USE obio_incom, only: ir

      implicit none


      integer i,j,l,k
      integer iant,isin,ispc,isat,iein,iepc,ieat,incp
     .       ,inca,inat,imed,inin,inpc,nr,ntot

      real rlat,rlon

      integer nir
      common /bir/ nir(nrg)

      real antlat,rnpolat
      data antlat,rnpolat /-40.0, 40.0/
 
c  Set up indicators
      iant = 1   !antarcic region
      isin = 2   !south indian ocean
      ispc = 3   !south pacific
      isat = 4   !south atlantic
      iein = 5   !equatorial indian ocean
      iepc = 6   !equatorial pacific ocean
      ieat = 7   !equatorial atlantic ocean
      inin = 8   !north indian ocean
      incp = 9   !north-central pacific
      inca = 10  !north-central atlantic
      inpc = 11  !north pacific
      inat = 12  !north atlantic
      imed = 13  !mediterranean/black sea
 
c  Initialize region indicator array
      ir = 0
      nir = 0
 
c  Find nwater values corresponding to regions
       do 1000 j=1,jdm
       do 1000 i=1,idm

        rlon=LON_DG(i,1)
        rlat=LAT_DG(j,1)

        if (rlon .gt. 180)rlon = rlon-360.0

        if (ZOE(LMOM(i,j)).le.0) go to 100

c   Antarctic region
        if (rlat .le. antlat)then
         ir(i,j) = iant
         nir(ir(i,j)) = nir(ir(i,j))+1
        endif

c   South Indian region
        if (rlat .le. -30.0 .and. rlat .gt. antlat)then
         if (rlon .gt. 20.0 .and. rlon .lt. 150.0)then
          ir(i,j) = isin
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .le. -10.0 .and. rlat .gt. -30.0)then
         if (rlon .gt. 20.0 .and. rlon .lt. 142.5)then
          ir(i,j) = isin
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif

c   South Pacific region
        if (rlat .le. -10.0 .and. rlat .gt. antlat)then
         if (rlon .le. -70.0 .and. rlon .ge. -180.0)then
          ir(i,j) = ispc
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
         if (rlon .ge. 150.0 .and. rlon .lt. 180.0)then
          ir(i,j) = ispc
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .le. -10.0 .and. rlat .gt. -30.0)then
         if (rlon .ge. 142.5 .and. rlon .lt. 180.0)then
          if (ir(i,j) .eq. 0)then
           ir(i,j) = ispc
           nir(ir(i,j)) = nir(ir(i,j))+1
          endif
         endif
        endif

c   South Atlantic region
        if (rlat .le. -10.0 .and. rlat .gt. antlat)then
         if (rlon .gt. -70.0 .and. rlon .le. 20.0)then
          ir(i,j) = isat
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif

c   Equatorial Indian Ocean region
        if (rlat .le. -8.0 .and. rlat .gt. -10.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 142.5)then
          ir(i,j) = iein
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .le. -6.0 .and. rlat .gt. -8.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 108.0)then
          if (ir(i,j) .eq. 0)then
           ir(i,j) = iein
           nir(ir(i,j)) = nir(ir(i,j))+1
          endif
         endif
        endif
        if (rlat .le. -4.0 .and. rlat .gt. -6.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 105.0)then
          ir(i,j) = iein
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .le. -2.0 .and. rlat .gt. -4.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 103.0)then
          ir(i,j) = iein
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .le. 2.0 .and. rlat .gt. -2.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 104.0)then
          ir(i,j) = iein
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .le. 4.0 .and. rlat .gt. 2.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 103.0)then
          ir(i,j) = iein
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .le. 7.0 .and. rlat .gt. 4.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 102.0)then
          ir(i,j) = iein
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .lt. 10.0 .and. rlat .gt. 7.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 99.0)then
          ir(i,j) = iein
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif

c   Equatorial Pacific region
        if (rlat .lt. 10.0 .and. rlat .gt. -10.0)then
         if (rlon .gt. 90.0)then
          if (ir(i,j) .eq. 0)then
           ir(i,j) = iepc
           nir(ir(i,j)) = nir(ir(i,j))+1
          endif
         endif
         if (rlon .gt. -180.0 .and. rlon .lt. -70.0)then
          if (ir(i,j) .eq. 0)then
           ir(i,j) = iepc
           nir(ir(i,j)) = nir(ir(i,j))+1
          endif
         endif
        endif

c   Equatorial Atlantic region
        if (rlat .lt. 10.0 .and. rlat .gt. -10.0)then
         if (rlon .lt. 20.0 .and. rlon .ge. -74.0)then
          if (ir(i,j) .eq. 0)then
           ir(i,j) = ieat
           nir(ir(i,j)) = nir(ir(i,j))+1
          endif
         endif
        endif

c   North Indian Ocean
        if (rlat .le. 30.0 .and. rlat .ge. 10.0)then
         if (rlon .gt. 20.0 .and. rlon .lt. 99.0)then
          ir(i,j) = inin
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif

c   North Central Pacific region
        if (rlat .ge. 10.0 .and. rlat .le. rnpolat)then
         if (rlon .gt. 99.0)then
          ir(i,j) = incp
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
         if (rlon .le. -100.0)then
          ir(i,j) = incp
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .ge. 10.0 .and. rlat .lt. 18.0)then
         if (rlon .gt. -100.0 .and. rlon .le. -90.0)then
          ir(i,j) = incp
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .ge. 10.0 .and. rlat .lt. 14.0)then
         if (rlon .gt. -90.0 .and. rlon .le. -84.5)then
          ir(i,j) = incp
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif

c   Mediterranean/Black Seas region
        if (rlat .ge. 30.0 .and. rlat .le. 43.0)then
         if (rlon .ge. -5.5 .and. rlon .lt. 60.0)then
          ir(i,j) = imed
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .gt. 43.0 .and. rlat .le. 48.0)then
         if (rlon .ge. 0.0 .and. rlon .lt. 70.0)then
          ir(i,j) = imed
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif

c   North Central Atlantic region
        if (rlat .ge. 7.0 .and. rlat .le. 40.0)then !pickup Carib.Sea
         if (ir(i,j) .eq. 0)then
          ir(i,j) = inca
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif

c  North Pacific and Atlantic
        if (rlat .gt. rnpolat)then
c  North Pacific
         if (rlon .lt. -105.0)then
          ir(i,j) = inpc
          nir(ir(i,j)) = nir(ir(i,j))+1
         else if (rlon .gt. 120.0)then
          ir(i,j) = inpc
          nir(ir(i,j)) = nir(ir(i,j))+1
         else
c  North Atlantic
          if (ir(i,j) .ne. imed)then
           ir(i,j) = inat
           nir(ir(i,j)) = nir(ir(i,j))+1
          endif
         endif
        endif
 100    continue

 1000 continue
 
c  Set nir to minimum 1 value to prevent error in division
      do nr = 1,nrg
       nir(nr) = max(nir(nr),1)
      enddo
 
c  Total up points for check
      ntot = 0
      do nr = 1,nrg
       ntot = ntot + nir(nr)
       write(6,*)'Region, no. points = ',nr,nir(nr)
      enddo
      write(6,*)'Total ocean points = ',ntot
 
      return
      end

c------------------------------------------------------------------------------
      subroutine bio_inicond_g(filename,fldo2,fldoz)

!read in a field and convert to ocean grid (using Gary Russel's routine) 

      USE FILEMANAGER, only: openunit,closeunit
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT

      USE OCEANRES, only : imo,jmo,lmo
      USE OCEAN, only : oDLATM=>DLATM,LMOM=>LMM,ZOE=>ZE,FOCEAN

      implicit none


      integer, parameter :: igrd=360,jgrd=180,kgrd=33
      integer i,j,k,l,n,lm
      integer iu_file,lgth
      real data(igrd,jgrd,kgrd)
      real data_mask(igrd,jgrd)
      real data_min(kgrd),data_max(kgrd)

      real fldo(imo,jmo,kgrd)
      real fldo2(imo,jmo,lmo)
      real fldoz(imo,jmo,lmo)

      real nodc_depths(kgrd)

      integer ii,jj,lcount
      real count_min

      logical vrbos

      character*80 filename

      data nodc_depths/0,  10,  20,  30,  50,  75, 100, 125, 150, 200,
     .          250, 300, 400, 500, 600, 700, 800, 900,1000,1100,1200,
     .    1300,1400,1500,1750,2000,2500,3000,3500,4000,4500,5000,5500/

!--------------------------------------------------------------

      !read no3 files from Watson and convert to ascii
      lgth=len_trim(filename)
      if (AM_I_ROOT())
     .print*, 'bioinit: reading from file...',filename(1:lgth)
      call openunit(filename,iu_file,.false.,.true.)

!NOTE: data starts from Greenwich
!      missing values are -9999
       do k=1,kgrd
       data_min(k)=1.e10
       data_max(k)=-1.e10
         do i=1,igrd
           do j=1,jgrd
             read(iu_file,'(e12.4)')data(i,j,k)
             !preserve the mean for later
              if (data(i,j,k)>0.) then
                data_min(k)=min(data_min(k),data(i,j,k))
                data_max(k)=max(data_max(k),data(i,j,k))
             endif
           enddo
         enddo
      enddo
      call closeunit(iu_file)

!--------------------------------------------------------------

      do k=1,kgrd

      do i=1,igrd
      do j=1,jgrd

        !mask
        data_mask(i,j)=0.d0
        if (data(i,j,k)>=0.d0) data_mask(i,j)=1.d0

cdiag   if (k.eq.1)
cdiag.    write(*,'(a,3i5,2e12.4)')'before hntr80 ',
cdiag.    i,j,1,data(i,j,k),data_mask(i,j)


      enddo    ! i-loop
      enddo    ! j-loop

!compute glb average and replace missing data

      call HNTR80(igrd,jgrd,180.d0,60.d0,
     .             imo,jmo,0.d0,oDLATM,-9999.d0)

      call HNTR8P (data_mask,data(:,:,k),fldo(:,:,k))

      enddo    ! k-loop

cdiag do i=1,imo
cdiag do j=1,jmo
cdiag  write(*,'(a,3i5,e12.4)')'after hntr8p ',
cdiag.    i,j,1,fldo(i,j,1)
cdiag enddo
cdiag enddo

!correction for ocean cells with missing values
      lcount=0
      do j=1,jmo
      do i=1,imo
      !count missing values
      if (lmom(i,j).gt.0.and.fldo(i,j,1).eq.-9999.d0)then
         lcount=lcount+1
      endif
      enddo 
      enddo 
      print*, 'old lcount = ', lcount

      do j=1,jmo
      do i=1,imo

      !replace missing values over ocean
      if (lmom(i,j).gt.0.and.fldo(i,j,1).eq.-9999.d0)then
         do jj=1,10
         do ii=1,10
         if (i+ii.le.imo.and.j+jj.le.jmo) then
         if (lmom(i+ii,j+jj).gt.0.and.
     .       fldo(i+ii,j+jj,1).ne.-9999.d0)then
             count_min=fldo(i+ii,j+jj,1)
             go to 100
         endif
         endif
         enddo
         enddo
 100  continue

      fldo(i,j,1)=count_min
      endif

      enddo 
      enddo 

      !recompute number of missing values over ocean
      lcount=0
      do j=1,jmo
      do i=1,imo
      !count missing values
      if (lmom(i,j).gt.0.and.fldo(i,j,1).eq.-9999.d0)then
         lcount=lcount+1
!        print*,'lcount=',i,j,lcount
      endif
      enddo 
      enddo 
      print*, 'new lcount = ', lcount


!--------------------------------------------------------------
      print*, 'imo,jmo= ',imo, jmo
      do j=1,jmo
      do i=1,imo

      do k=1,lmo
      fldo2(i,j,k)= -9999.
      fldoz(i,j,k)= -9999.
      enddo

cdiag do k=1,kgrd
cdiag write(*,'(a,3i5,e12.4)')'befr VLKtoLZ ',i,j,k,fldo(i,j,k)
cdiag enddo

      IF (FOCEAN(i,j).gt.0) then
         lm=lmom(i,j)
          call VLKtoLZ(kgrd,lm,nodc_depths,ZOE,
     .                 fldo(i,j,:),fldo2(i,j,:),fldoz(i,j,:))   
      ENDIF

cdiag do k=1,lmo
cdiag write(*,'(a,3i5,2e12.4)')'aftr VLKtoLZ ',i,j,k,ZOE(k),fldo2(i,j,k)
cdiag enddo

      enddo
      enddo
  
      return
  
      end subroutine bio_inicond_g



      Subroutine VLKtoLZ (KM,LM, MK,ME, RK, RL,RZ)        !  2008/02/05
C****
C**** VLKtoLZ assumes a continuous piecewise linear tracer distribution,
C**** defined by input tracer concentrations RK at KM specific points.
C**** MK in the downward vertical mass coordinate.
C**** R(M) = {RK(K-1)*[MK(K)-M] + RK(K)*[M-MK(K-1)]} / [MK(K)-MK(K-1)]
C****               when MK(K-1) < M < MK(K).
C**** R(M) = RK(1)  when M < MK(1).
C**** R(M) = is undefined when MK(KM) < M.
C****
C**** VLKtoLZ integrates this tracer distribution over the LM output
C**** layers defined by their layer edges ME, calculating the tracer
C**** mass RM of each layer and the vertical gradient RZ.
C**** RNEW(M) = RL(L) + RZ(L)*[M-MC(L)]/dM(L) when ME(L-1) < M < ME(L)
C**** where MC(L) = .5*[ME(L-1)+ME(L)] and dM(L) = ME(L)-ME(L-1)
C**** Mean concentration of output layers is RL(L) = RM(L)/dM(L).
C****
C**** If ME(L-1) < MK(KM) < ME(L), then RL(L) and RZ(L) are calculated
C**** from the input profile up to MK(KM); RL(L+1:LM) and RZ(L+1:LM)
C**** for deeper layers are undefined, set to DATMIS.
C****
C**** Input:  KM = number of input edges
C****         LM = number of output cells
C****         MK = mass coordinates of input points (kg/m^2)
C****         ME = mass coordinates of output layer edges (kg/m^2)
C****         RK = tracer concentration at input points
C****
C**** Output: RL = mean tracer concentration of each output layer
C****         RZ = vertical gradient of tracer mass of each output layer
C****
C**** Internal: RM = integrated tracer mass of output layers (kg/m^2)
C****           RQ = integrated tracer mass times mass (kg^2/m^4)
C****
      Implicit Real*8 (A-H,M-Z)
      Parameter (DATMIS = -999999)
      Real*8 MK(KM),ME(0:LM), RK(KM), RL(LM),RZ(LM), RM(1024),RQ(1024)



C     If (LM > 1024)  Stop 'LM exceeds internal dimentions in VLKtoLZ'
C****
      RM(1:LM) = 0
      RQ(1:LM) = 0
      K = 1
      L = 1
      MC = .5*(ME(L)+ME(L-1))
C****
C**** Integrate layers with M < MK(1)
C****
      If (ME(0) < MK(1))  GoTo 20
C**** MK(1) <= ME(0), determine K such that MK(K-1) <= ME(0) < MK(K)
   10 If (K == KM)  GoTo 200  ;  K = K+1
      If (MK(K) <= ME(0))  GoTo 10
      GoTo 130  !  MK(K-1) <= ME(0) < MK(K)
C**** ME(0) < MK(1), determine output cell containing MK(1)
   20 If (MK(1) < ME(L))  GoTo 30
C**** ME(L-1) < ME(L) < MK(1), integrate RM from ME(L-1) to ME(L)
      RM(L) = RK(1)*(ME(L)-ME(L-1))
      RQ(L) = 0
      If (L == LM)  GoTo 300  ;  L = L+1  ;  MC = .5*(ME(L)+ME(L-1))
      GoTo 20
C**** ME(L-1) < MK(1) < ME(L), integrate RM from ME(L-1) to MK(1)
   30 RM(L) = RK(1)*(MK(1)-ME(L-1))
      RQ(L) = RK(1)*(MK(1)-ME(L-1))*(.5*(MK(1)+ME(L-1))-MC)
      If (K == KM)  GoTo 220  ;  K = K+1
C****
C**** Integrate layers with MK(1) < M < MK(KM)
C****
  100 If (ME(L) < MK(K))  GoTo 120
C**** ME(L-1) < MK(K-1) < MK(K) < ME(L), integrate from MK(K-1) to MK(K)
      RM(L) = RM(L) + (RK(K)-RK(K-1))*(MK(K)+MK(K-1))/2 +
     +                RK(K-1)*MK(K)-RK(K)*MK(K-1)
      RQ(L) = RQ(L) +
     +  (RK(K)-RK(K-1))*(MK(K)*MK(K)+MK(K)*MK(K-1)+MK(K-1)*MK(K-1))/3 +
     +  (RK(K-1)*(MK(K)+MC)-RK(K)*(MK(K-1)+MC))*(MK(K)+MK(K-1))/2 +
     +  (RK(K)*MK(K-1)-RK(K-1)*MK(K))*MC
      If (K == KM)  GoTo 220  ;  K = K+1
      GoTo 100
C**** ME(L-1) < MK(K-1) < ME(L) < MK(K), integrate from MK(K-1) to ME(L)
  120 RM(L) = RM(L) + ((RK(K)-RK(K-1))*(ME(L)+MK(K-1))/2 +
     +                 (RK(K-1)*MK(K)-RK(K)*MK(K-1))) * (ME(L)-MK(K-1))
     /              / (MK(K)-MK(K-1))
      RQ(L) = RQ(L) +
     +  ((RK(K)-RK(K-1))*(ME(L)*ME(L)+ME(L)*MK(K-1)+MK(K-1)*MK(K-1))/3 +
     +   (RK(K-1)*(MK(K)+MC)-RK(K)*(MK(K-1)+MC))*(ME(L)+MK(K-1))/2 +
     +   (RK(K)*MK(K-1)-RK(K-1)*MK(K))*MC) * (ME(L)-MK(K-1)) /
     /  (MK(K)-MK(K-1))
      If (L == LM)  GoTo 300  ;  L = L+1  ;  MC = .5*(ME(L)+ME(L-1))
  130 If (MK(K) < ME(L))  GoTo 160
C**** MK(K-1) < ME(L-1) < ME(L) < MK(K), integrate from ME(L-1) to ME(L)
  140 RM(L) = ((RK(K)-RK(K-1))*(ME(L)+ME(L-1))/2 +
     +         (RK(K-1)*MK(K)-RK(K)*MK(K-1))) * (ME(L)-ME(L-1)) /
     /        (MK(K)-MK(K-1))
      RQ(L) =
     +  ((RK(K)-RK(K-1))*(ME(L)*ME(L)+ME(L)*ME(L-1)+ME(L-1)*ME(L-1))/3 +
     +   (RK(K-1)*(MK(K)+MC)-RK(K)*(MK(K-1)+MC))*(ME(L)+ME(L-1))/2 +
     +   (RK(K)*MK(K-1)-RK(K-1)*MK(K))*MC) * (ME(L)-ME(L-1)) /
     /  (MK(K)-MK(K-1))
      If (L == LM)  GoTo 300  ;  L = L+1  ;  MC = .5*(ME(L)+ME(L-1))
      If (ME(L) < MK(K))  GoTo 140
C**** MK(K-1) < ME(L-1) < MK(K) < ME(L), integrate from ME(L-1) to MK(K)
  160 RM(L) = RM(L) + ((RK(K)-RK(K-1))*(MK(K)+ME(L-1))/2 +
     +                 (RK(K-1)*MK(K)-RK(K)*MK(K-1))) * (MK(K)-ME(L-1))
     /              / (MK(K)-MK(K-1))
      RQ(L) = RQ(L) +
     +  ((RK(K)-RK(K-1))*(MK(K)*MK(K)+MK(K)*ME(L-1)+ME(L-1)*ME(L-1))/3 +
     +   (RK(K-1)*(MK(K)+MC)-RK(K)*(MK(K-1)+MC))*(MK(K)+ME(L-1))/2 +
     +   (RK(K)*MK(K-1)-RK(K-1)*MK(K))*MC) * (MK(K)-ME(L-1)) /
     /  (MK(K)-MK(K-1))
      If (K == KM)  GoTo 220  ;  K = K+1
      GoTo 100
C****
C**** Calculate RL and RZ from RM and RQ when MK(KM) < ME(LM)
C****
C**** MK(KM) <= ME(0)
  200 RL(:) = DATMIS
      RZ(:) = DATMIS
      Return
C**** ME(L-1) < MK(KM) < ME(L)
  220 Do 230 LL=1,L-1
      RL(LL) =   RM(LL) / (ME(LL)-ME(LL-1))
  230 RZ(LL) = 6*RQ(LL) / (ME(LL)-ME(LL-1))**2
      RL(L)  =   RM(L)  / (MK(KM)-ME(L-1))
      RZ(L)  = 6*(RQ(L) + .5*(ME(L)-MK(KM))*RM(L)) / (MK(KM)-ME(L-1))**2
C**** Vertical gradient is extrapolated half way to .5*[MK(KM)+ME(L)]
      RZ(L)  = RZ(L) * (.5*(MK(KM)+ME(L))-ME(L-1)) / (MK(KM)-ME(L-1))
      RL(L+1:LM) = DATMIS
      RZ(L+1:LM) = DATMIS
      Return
C****
C**** Calculate RL and RZ from RM and RQ when ME(LM) < MK(KM)
C****
  300 Do 310 L=1,LM
      RL(L) =   RM(L) / (ME(L)-ME(L-1))
  310 RZ(L) = 6*RQ(L) / (ME(L)-ME(L-1))**2
      Return
      End 

#endif /*  OBIO_ON_GARYocean */
