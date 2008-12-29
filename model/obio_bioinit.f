#include "rundeck_opts.h"

#ifndef OBIO_ON_GARYocean     !not on Gary's ocean
      subroutine obio_bioinit(nn)
 
!note: 
!obio_bioinit is called only for a cold start and reads in INITIAL conditions and interpolates them
!to the hycom grid. Such fields are nitrates,silicate,dic
!obio_init  is called for every start of the run and reads in BOUNDARY conditions and interpolates 
!them to hycom grid. such fields are iron,alkalinity,chlorophyl

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
      USE obio_com, only: gcmax
 
      USE hycom_dim_glob, only : jj,isp,ifp,ilp,iia,jja,iio,jjo,kdm
     &     ,idm,jdm
      USE hycom_dim, only : isp_l=>isp,ifp_l=>ifp,ilp_l=>ilp
      USE hycom_dim, only : j_0,j_1
      !USE hycom_arrays_glob, only : tracer_glob => tracer
      USE hycom_arrays_glob, only : tracer,dpinit
      USE hycom_scalars, only: onem

      USE DOMAIN_DECOMP, only: AM_I_ROOT
cddd      USE hycom_arrays_glob, only: gather_hycom_arrays,
cddd     &     scatter_hycom_arrays

      implicit none

      integer i,j,k,l,nn


      integer nir,nt
      integer iu_bioinit

      real rlon,rlat,dicmin,dicmax
      real zz

      character*2 ntchar
      character*80 filename

      !real dic_glob(iio,jjo,kdm)

      common /bir/ nir(nrg)

      call alloc_obio_incom

      !call gather_hycom_arrays ! shoul remove later
      call gather_tracer
      call gather_dpinit

      if (AM_I_ROOT()) then
c  Initialize

      tracer(:,:,:,1:ntyp)=0.d0
      Fer(:,:,:) = 0.d0

      filename='nitrates_inicond'
      call bio_inicond(filename,tracer(:,:,:,1))

      filename='silicate_inicond'
      call bio_inicond(filename,tracer(:,:,:,3))

      filename='dic_inicond'
      call bio_inicond(filename,dic(:,:,:))

      !dic(:,j_0h:j_1h,:) = dic_glob(:,:,:)

!     /archive/u/aromanou/Watson_new/BioInit/iron_ron_4x5.asc
!     /archive/u/aromanou/Watson_new/BioInit/CHL_WG_4x5


!!these rno3 and so2_init fields are not correct. There are void points due to
!mismatch of the noaa grid and the hycom grid. To fill in have to do
!the interpolation in matlab (furtuna).
!at the same time use dps and interpolate to layer depths from the model
      do k=1,kdm
       do j=1,jdm
        do i=1,idm
          if(tracer(i,j,k,1).le.0.)tracer(i,j,k,1)=0.085
        enddo
       enddo
      enddo

      do k=1,kdm
       do j=1,jdm
        do i=1,idm
          if(tracer(i,j,k,3).le.0.)tracer(i,j,k,3)=0.297
        enddo
       enddo
      enddo
     
      dicmin = 1.e30
      dicmax =-1.e30
      do k=1,kdm
       do j=1,jdm
        do i=1,idm
          if (dic(i,j,k).le.0.) dic(i,j,k)=1837.
          dic(i,j,k)=max(dic(i,j,k),1837.)   !set minimum =1837
          dicmin =min(dicmin,dic(i,j,k))
          dicmax =max(dicmax,dic(i,j,k))
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

      do j=1,jj
       do k=1,kdm
        do l=1,isp(j)
         do i=ifp(j,l),ilp(j,l)

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
         enddo  !i-loop
        enddo  !l-loop
       enddo  !k-loop
      enddo  !j-loop
 
c  Create arrays 
      write(6,*)'Creating bio restart data for ',ntyp,' arrays and'
     . ,kdm,'  layers...'

! the following loop has parallelization problems....

      do 1000 j=1,jj
      do 1000 l=1,isp(j)
      do 1000 i=ifp(j,l),ilp(j,l)

        zz=0.d0
        do k=1,kdm
          !Nitrate
          !read earlier from file

          zz=zz+dpinit(i,j,k)/onem   !depth in m

          !Ammonium
          tracer(i,j,k,2) = 0.5
          if (zz.gt.4000.)tracer(i,j,k,2) = Pdeep(2)

          !Silica
          !read earlier from file
!         if (zz.gt. 4000.)tracer(i,j,k,2) = Pdeep(3)

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
          if (zz.gt. 4000.)tracer(i,j,k,4) = Pdeep(4)

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

         enddo

 1000 continue


c  Detritus (set to 0 for start up)
      write(6,*)'Detritus...'
      cnratio = 106.0/16.0*12.0    !C:N ratio (ugl:uM)
      csratio = 106.0/16.0*12.0    !C:Si ratio (ugl:uM)
      cfratio = 150000.0*12.0*1.0E-3    !C:Fe ratio (ugl:nM)

      do j=1,jj
       do l=1,isp(j)
        do i=ifp(j,l),ilp(j,l)
         do k=1,kdm
          if (dpinit(i,j,k)/onem .gt. 0.0)then
           !only detritus components
           tracer(i,j,k,ntyp+n_inert+1) = tracer(i,j,k,1)*0.25*cnratio !as carbon
           tracer(i,j,k,ntyp+n_inert+2) = tracer(i,j,k,3)*0.1
           tracer(i,j,k,ntyp+n_inert+3) = tracer(i,j,k,4)*0.25
           tracer(i,j,k,ntyp+n_inert+1) = 0.0
           tracer(i,j,k,ntyp+n_inert+2) = 0.0
           tracer(i,j,k,ntyp+n_inert+3) = 0.0
          endif
         enddo
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
      do j=1,jj
       do l=1,isp(j)
        do i=ifp(j,l),ilp(j,l)
         do k=1,kdm
          tracer(i,j,k,ntyp+n_inert+ndet+1) = 0.0
          tracer(i,j,k,ntyp+n_inert+ndet+2) = 0.0
         enddo
        enddo
       enddo
      enddo

      !only carbon components
      do j=1,jj
       do l=1,isp(j)
        do i=ifp(j,l),ilp(j,l)
         do k = 1,kdm
          tracer(i,j,k,ntyp+n_inert+ndet+2) = dicmod(i,j,k)
         enddo
c         car(i,j,k,1) = 3.0  !from Bissett et al 1999 (uM(C))
c         car(i,j,k,1) = 0.0  !from Walsh et al 1999
        enddo
       enddo
      enddo

      endif
      !call scatter_hycom_arrays
      call scatter_tracer

c  Light saturation data
      if (AM_I_ROOT()) then
      write(6,*)'Light saturation data...'
      endif
      avgq = 0.0

      do j=j_0,j_1
       do k=1,kdm
        do l=1,isp_l(j)
         do i=ifp_l(j,l),ilp_l(j,l)
          avgq(i,j,k) = 25.0
         enddo
        enddo
       enddo
      enddo

c  Coccolithophore max growth rate
      if (AM_I_ROOT()) then
      write(6,*)'Coccolithophore max growth rate...'
      endif
      do j=j_0,j_1
       do k=1,kdm
        do l=1,isp_l(j)
         do i=ifp_l(j,l),ilp_l(j,l)
          gcmax(i,j,k) = 0.0
         enddo
        enddo
       enddo
      enddo
 
      !save initialization
cdiag do nt=1,ntyp+n_inert+ndet+ncar
cdiag   ntchar='00'
cdiag   if(nt.le.9)write(ntchar,'(i1)')nt
cdiag   if(nt.gt.9)write(ntchar,'(i2)')nt
cdiag   print*,'BIO: saving initial tracer fields '
cdiag.        ,'bioinit_tracer'//ntchar
cdiag   call openunit('bioinit_tracer'//ntchar,iu_bioinit)
cdiag   do k=1,kdm
cdiag   do j=j_0,j_1				!  do not parallelize
cdiag   do l=1,isp_l(j)
cdiag   do i=ifp_l(j,l),ilp_l(j,l)
cdiag      write(iu_bioinit,'(3i4,2e12.4)')
cdiag.           i,j,k,dpinit(i,j,k)/onem,tracer(i,j,k,nt)
cdiag   enddo
cdiag   enddo
cdiag   enddo
cdiag   enddo
cdiag call closeunit(iu_bioinit)
cdiag enddo

      if (AM_I_ROOT())
     .print*,'COLD INITIALIZATION'
      call obio_trint(nn)

      return
      end

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
 

      USE obio_dim
      USE obio_incom, only: ir

      USE hycom_dim_glob, only : jj,isp,ifp,ilp
      USE hycom_arrays_glob, only : lonij,latij,dpinit
      USE hycom_scalars, only: onem

      implicit none


      integer i,j,l
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
      do 1000 j=1,jj
      do 1000 l=1,isp(j)
      do 1000 i=ifp(j,l),ilp(j,l)

        rlon=lonij(i,j,3)
        rlat=latij(i,j,3)

        if (rlon .gt. 180)rlon = rlon-360.0

        if (dpinit(i,j,1)/onem .le. 0.0)go to 100

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
100     continue

!       write(121,'(2(i4,1x),4(f8.3,1x),e12.4,i4)')
!    .        i,j,lonij(i,j,3),latij(i,j,3),
!    .        rlon,rlat,dpinit(i,j,1)/onem,ir(i,j)
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
      subroutine bio_inicond(filename,fldo2)

!read in a field and convert to ocean grid (using Shana's routine) 
!convert to the new grid using dps from the model and the remap routine

c --- mapping flux-like field from agcm to ogcm
c     input: flda (W/m*m), output: fldo (W/m*m)
c

      USE FILEMANAGER, only: openunit,closeunit
      USE DOMAIN_DECOMP, only: AM_I_ROOT

      USE hycom_dim_glob, only : jj,isp,ifp,ilp,iia,jja,iio,jjo,kdm
      USE hycom_arrays_glob, only : dpinit
      USE hycom_scalars, only: onem

      USE hycom_cpler, only: wlista2o,ilista2o,jlista2o,nlista2o

      implicit none


      integer, parameter :: igrd=360,jgrd=180,kgrd=33
      integer i,j,k,l,n
      integer iu_file,lgth
      integer i1,j1,iii,jjj,isum,kmax
      integer nodc_kmax

      real data(igrd,jgrd,kgrd)
      real data2(iia,jja,kgrd)
      real data_min(kgrd),data_max(kgrd)
      real sum1
      real dummy1(36,jja,kgrd),dummy2(36,jja,kgrd)
      real fldo(iio,jjo,kgrd)
      real pinit(iio,jjo,kdm+1),fldo2(iio,jjo,kdm)
      real nodc_depths(kgrd),nodc_d(kgrd+1)
      !real dpinit(iio,jjo,kdm)

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
       i1=0
       do i=1,igrd,5
       i1=i1+1
       j1=0
       do j=1,jgrd,4
       j1=j1+1

        sum1=0.
        isum=0
        do iii=1,5
        do jjj=1,4

            if (data(i+iii-1,j+jjj-1,k).ge.0.) then
               sum1=sum1+data(i+iii-1,j+jjj-1,k)
               isum=isum+1
            endif
         enddo
         enddo

          if (isum.gt.0) then
             data2(i1,j1,k)=sum1/float(isum)
          else
             data2(i1,j1,k)=-9999.
          endif
      enddo
      enddo
      enddo

!     !this is needed for dic
!     do k=1,kgrd
!     do i=1,iia
!     do j=1,jja
!     if (data2(i,j,k).gt.0) then
!       if (data2(i,j,k)<data_min(k)) data2(i,j,k)=data_min(k)
!       if (data2(i,j,k)>data_max(k)) data2(i,j,k)=data_max(k)
!     endif
!     enddo
!     enddo
!     enddo

      !set to zero so that interpolation works
      do k=1,kgrd
      do i=1,iia
      do j=1,jja
        if (data2(i,j,k)<0.)data2(i,j,k)=0.
      enddo
      !make it 72x46x33
      data2(i,46,k)=data2(i,45,k)
      enddo
      enddo

      !--------------------------------------------------------
      !***************** important! change to dateline
      !move to dateline
      dummy1=data2(1:36,:,:);
      dummy2=data2(37:72,:,:);
      data2(1:36,:,:)=dummy2;
      data2(37:72,:,:)=dummy1;

      !--------------------------------------------------------
c$OMP PARALLEL DO
      do 8 j=1,jj
      do 8 l=1,isp(j)
      do 8 i=ifp(j,l),ilp(j,l)

      do 9 k=1,kgrd
      fldo(i,j,k)=0.
c
      do 9 n=1,nlista2o(i,j)
      fldo(i,j,k)=fldo(i,j,k)
     .           +data2(ilista2o(i,j,n),jlista2o(i,j,n),k)
     .                       *wlista2o(i,j,n)
 9    continue
 8    continue
c$OMP END PARALLEL DO

!     !this is needed for dic
!     do k=1,kgrd
!     do i=1,iio
!     do j=1,jjo
!       if (fldo(i,j,k)<data_min(k)) fldo(i,j,k)=data_min(k)
!       if (fldo(i,j,k)>data_max(k)) fldo(i,j,k)=data_max(k)
!     enddo
!     enddo
!     enddo

      !--------------------------------------------------------
      !use dpinit(i,j,k)/onem

       pinit(:,:,1)=0.d0
c$OMP PARALLEL DO
       do 10 j=1,jj
       do 10 l=1,isp(j)
       do 10 i=ifp(j,l),ilp(j,l)
       do  k=1,kdm
         pinit(i,j,k+1)=pinit(i,j,k)+dpinit(i,j,k)/onem
       enddo
 10    continue
c$OMP END PARALLEL DO

       fldo2(:,:,:)=-9999.d0
c$OMP PARALLEL DO PRIVATE(kmax,nodc_d,nodc_kmax,vrbos)
       do j=1,jj                       
       do l=1,isp(j)
       do i=ifp(j,l),ilp(j,l)

          !match nodc and model bottom pressure
          !the top match already
          do k=2,kdm+1
           if (pinit(i,j,k) .gt. pinit(i,j,k-1)) kmax=k
          enddo

          !model bottom at kmax+1
          do k=1,kgrd
           if (nodc_depths(k) .le. pinit(i,j,min(21,kmax+1))) then
               nodc_d(k)=nodc_depths(k)
               nodc_kmax=k
           endif
cdiag      write(*,'(a,3i5,2e12.4,i5)')'bioinit: ',
cdiag.               i,j,k,fldo(i,j,k),nodc_d(k),nodc_kmax
          enddo
          nodc_d(nodc_kmax+1)=pinit(i,j,min(21,kmax+1))

!        call remap1d_pcm(fldo(i,j,1:nodc_kmax),nodc_d,nodc_kmax,
!    .             fldo2(i,j,:),pinit(i,j,:),kdm,.false.,i,j)

    
         vrbos=.false.
         call remap1d_plm(fldo(i,j,1:nodc_kmax),nodc_d,nodc_kmax,
     .             fldo2(i,j,1:kdm),pinit(i,j,1:kdm+1),kdm,vrbos,i,j)

       enddo
       enddo
       enddo
c$OMP END PARALLEL DO

      !--------------------------------------------------------

       return
  
       end subroutine bio_inicond
  
      subroutine remap1d_plm(yold,xold,kold,ynew,xnew,knew,vrbos,i,j)
c
c --- consider two stepwise constant functions -yold,ynew- whose
c --- discontinuities are at abscissa values -xold,xnew- respectively.
c --- treat -ynew- as unknown. solve for -ynew- under the condition that
c --- the integral over y*dx is preserved (integration based on PLM).
c
      implicit none
      integer,intent(IN) :: kold,knew
      integer,intent(IN) :: i,j                !  current location in horiz.grid
      real,intent(IN)    :: yold(kold),xold(kold+1),xnew(knew+1)
      real,intent(OUT)   :: ynew(knew)
      logical,intent(IN) :: vrbos        !  if true, print diagnostics
      integer k,ko,n
      real colin,clout,slope,wgta,wgtb,wgtc,yinteg,
     .     ylft(kold),yrgt(kold),xlo,xhi,ra,rb,ya,yb,q,plmslp,
     .     yrka,ylk,yrk,ylkb
      external plmslp
      logical at_top
      real,parameter    :: onemu=1.e-6, acurcy=1.e-6, flag=-999.
      integer,parameter :: iter=2
c
      if (vrbos)
     . write (*,101) i,j,' remap1d -- old profile:     x         y',
     .  (k,xold(k),yold(k),k=1,kold),kold+1,xold(kold+1)
 101  format (2i5,a/(i30,f12.1,f10.2))
c
c --- column integrals (colin/clout) are computed for diagnostic purposes only
      colin=0.
      clout=0.
      do 3 k=1,kold
 3    colin=colin+yold(k)*(xold(k+1)-xold(k))
c
c --- replace each flat segment of stairstep curve by
c --- a slanting segment, using PLM-type limiters.
c
      ylft(   1)=yold(   1)
      yrgt(   1)=yold(   1)
      ylft(kold)=yold(kold)
      yrgt(kold)=yold(kold)
      do 6 n=1,iter                        !  iterate to optimize limiters
      do 2 k=2,kold-1
      if (n.eq.1) then
        yrka=yold(k-1)
        ylk=yold(k)
        yrk=yold(k)
        ylkb=yold(k+1)
      else
        yrka=yrgt(k-1)
        ylk=ylft(k)
        yrk=yrgt(k)
        ylkb=ylft(k+1)
      end if
      wgta=max(onemu,xold(k  )-xold(k-1))
      wgtb=max(onemu,xold(k+1)-xold(k  ))
      wgtc=max(onemu,xold(k+2)-xold(k+1))
      if (k.eq.     1) wgta=onemu
      if (k.eq.kold-1) wgtc=onemu
      slope=plmslp((wgtb*yrka+wgta*ylk)/(wgtb+wgta),
     .     yold(k),(wgtb*ylkb+wgtc*yrk)/(wgtb+wgtc))
      ylft(k)=yold(k)-slope
 2    yrgt(k)=yold(k)+slope
      if (vrbos) print '(8x,a,12x,a,5x,a/(i3,f9.2,5x,2f9.2))',
     .  'y','ylft','yrgt',(ko,yold(ko),ylft(ko),yrgt(ko),ko=1,kold)
 6    continue
      if (vrbos) print '(a/(10f8.2))',
     .  '  target values:',(ynew(ko),ko=1,knew)
c
c --- y in k-th interval now varies from ylft at xold(k) to yrgt at xold(k+1).
c --- find ynew(k) by requiring
c --- that the integral over y*dx from xnew(k) to xnew(k+1) be preserved.
c
      at_top=.true.
      do 4 k=1,knew
      yinteg=0.
      xlo=xnew(k  )
      xhi=xnew(k+1)
ccc      if (vrbos) print '(a,2f9.3)','xlo,xhi =',xlo,xhi
      if (xhi.gt.xlo) then
        at_top=.false.
        do 5 ko=1,kold
        if (xold(ko).ge.xhi) go to 1
c --- integrate over sloping portions of y(x) curve:
        ra=max(xlo,min(xhi,xold(ko  )))
        rb=max(xlo,min(xhi,xold(ko+1)))
cnat    ya=ylft(k)
cnat    yb=yrgt(k)
        ya=ylft(ko)
        yb=yrgt(ko)
        wgta=flag
        wgtb=flag
        if (xold(ko+1).ne.xold(ko)) then
          if (ra.ge.xold(ko).and.ra.le.xold(ko+1)) then
            wgta=(xold(ko+1)-ra)/(xold(ko+1)-xold(ko))
            ya=ylft(ko)*wgta+yrgt(ko)*(1.-wgta)
          end if
          if (rb.ge.xold(ko).and.rb.le.xold(ko+1)) then
            wgtb=(xold(ko+1)-rb)/(xold(ko+1)-xold(ko))
            yb=ylft(ko)*wgtb+yrgt(ko)*(1.-wgtb)
          end if
        end if
        yinteg=yinteg+.5*(ya+yb)*(rb-ra)
ccc        if (vrbos) print '(2i4,4f9.3,3f11.1)',
ccc     .    k,ko,ra,rb,wgta,wgtb,ya,yb,yinteg
 5      continue
        yinteg=yinteg+yb*(xhi-rb)
ccc        if (vrbos) print '(2i4,4f9.3,3f11.1)',
ccc     .    k,0,rb,xhi,wgta,wgtb,yb,yb,yinteg
 1      ynew(k)=yinteg/(xhi-xlo)
      else if (at_top) then
        ynew(k)=yold(   1)
      else                              !  at end
        ynew(k)=yold(kold)
      end if
ccc      if (vrbos) print '(a,f11.1)','ynew =',ynew(k)
      clout=clout+ynew(k)*(xnew(k+1)-xnew(k))
 4    continue
c
      if (abs(clout-colin).gt.acurcy*10.*xold(kold+1))
     .  write (*,100) i,j,' remap1d - column intgl.error',
     .   colin,clout,(clout-colin)/colin
 100  format (2i5,a,2es14.6,es9.1)
c
      if (vrbos)
     . write (*,101) i,j,' remap1d -- new profile:     x         y',
     .  (k,xnew(k),ynew(k),k=1,knew),knew+1,xnew(knew+1)

      return
      end
c
c
      function plmslp(ylft,ymid,yrgt)
c
c --- get slope at point 'ymid' for piecewise linear interpolation
      if (ymid.le.min(ylft,yrgt) .or.
     .    ymid.ge.max(ylft,yrgt)) then
        plmslp=0.
      else if ((yrgt-ylft)*(ylft+yrgt-2.*ymid).gt.0.) then
        plmslp=ymid-ylft
      else
        plmslp=yrgt-ymid
      end if
      return
      end

      subroutine remap1d_pcm(yold,xold,kold,ynew,xnew,knew,vrbos,i,j)
c
c --- consider two stepwise constant functions -yold,ynew- whose
c --- discontinuities are at abscissa values -xold,xnew- respectively.
c --- treat -ynew- as unknown. solve for -ynew- under the condition that
c --- the integral over y*dx is preserved (integration based on PCM).
c
      implicit none
      integer,intent(IN) :: kold,knew
      integer,intent(IN) :: i,j                !  current location in horiz.grid
      real,intent(IN)    :: yold(kold),xold(kold+1),xnew(knew+1)
      real,intent(OUT)   :: ynew(knew)
      logical,intent(IN) :: vrbos        !  if true, print diagnostics
      integer k,ko
      real colin,clout,colmx,yinteg,xlo,xhi,xa,xb
      logical at_top
      real,parameter    :: acurcy=1.e-6
c
      if (vrbos)
     . write (*,101) i,j,' remap1d -- old profile:     x         y',
     .  (k,xold(k),yold(k),k=1,kold),kold+1,xold(kold+1)
 101  format (2i5,a/(i30,f12.1,f10.3))
c
c --- column integrals (colin/clout) are computed for diagnostic purposes only
      colin=0.
      clout=0.
      colmx=0.
      do 3 k=1,kold
      colmx=max(colmx,abs(yold(k)))
 3    colin=colin+yold(k)*(xold(k+1)-xold(k))
c
c --- find ynew(k) by requiring that the integral over y*dx
c --- from xnew(k) to xnew(k+1) be preserved.
c
      at_top=.true.
      do 4 k=1,knew
      yinteg=0.
      xlo=xnew(k  )
      xhi=xnew(k+1)
      if (xhi.gt.xlo) then
        at_top=.false.
        xb=xlo
        do 5 ko=1,kold
        xa=xb
        xb=min(xhi,max(xlo,xold(ko+1)))
        yinteg=yinteg+yold(ko)*(xb-xa)
        if (xa.ge.xhi) go to 1
 5      continue
        xa=xb
        xb=xhi
        yinteg=yinteg+yold(kold)*(xb-xa)
 1      ynew(k)=yinteg/(xhi-xlo)
        clout=clout+yinteg
      else if (at_top) then
        ynew(k)=yold(   1)
      else                              !  at end
        ynew(k)=yold(kold)
      end if
 4    continue
c
      if (abs(clout-colin).gt.acurcy*colmx*xold(kold+1))
     .  write (*,100) i,j,' remap1d - column intgl.error',
     .   colin,clout,(clout-colin)/colin
 100  format (2i5,a,2es14.6,es9.1)
c
      if (vrbos)
     . write (*,101) i,j,' remap1d -- new profile:     x         y',
     .  (k,xnew(k),ynew(k),k=1,knew),knew+1,xnew(knew+1)
      return
      end

#endif    /* not on gary's ocean */
