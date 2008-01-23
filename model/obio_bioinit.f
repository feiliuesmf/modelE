      subroutine obio_bioinit
 
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
 
      implicit none
#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"



      integer igrd,jgrd,kgrd,ir,nir,nt,nerror,ierror,nflg,ntunit
      parameter(igrd=360,jgrd=180,kgrd=nlt)

      integer ichan,ilat,ilon,idic,lgth
      integer kb(kgrd),ipcrs
      integer iu_bioinit,iu_dpinit


      real rno3min,rno3max,so2min,so2max
      real dictot,rlon,rlat,dicavg,startdic
      real rinc
      real Fer(idm,jdm,kdm)
      real rno3(idm,jdm,kgrd),so2(idm,jdm,kgrd)
      real*4 real4(idm,jdm)
      real dicmin,dicmax

      logical save4ini   !switch for when we need to correct ini fields

      real spval
      data spval/-.03125/

      integer jn,jjn,in,iin
      real dic(idm,jdm,kgrd)
      real dicmod(idm,jdm,kdm)
      character*2 ntchar
      character*21 file_dir   !discover
!     character*32 file_dir   !explorer

      common /refin/ ipcrs(idm,jdm)
      common /bir/ ir(idm,jdm),nir(nrg)
      real flag
      !data flag /999.0E9/
      data flag /-99./
      data kb /0,10,20,30,50,75,100,125,150,200,250,300,400,500,600,
     *700,800,900,1000,1100,1200,1300,1400,1500,1750,2000,2500,3000,
     *3500,4000,4500,5000,5500/
 
   
      save4ini=.false.

c  Initialize

c$OMP PARALLEL DO
      do j=1,jj
       do nt= 1,ntyp
        do k = 1,kdm
         do l=1,isp(j)
          do i=ifp(j,l),ilp(j,l)
           tracer(i,j,k,nt) = 0.0
          enddo
         enddo
        enddo
       enddo
       do k = 1,kdm
        do l=1,isp(j)
         do i=ifp(j,l),ilp(j,l)
         Fer(i,j,k) = 0.0
         enddo
        enddo
       enddo
      enddo
c$OMP END PARALLEL DO

      nerror = 0
      ierror = 0
 
      do 14 lgth=60,1,-1
      if (flnmovt(lgth:lgth).ne.' ') go to 15
 14   continue
      stop '(bioinit:flnmovt)'
 15   continue  

c  read NOAA 2001 atlas data
      if (save4ini) then

 
!      open(4,file=flnmbio1,status='unknown')
!      do ichan=1,kgrd
!        do ilon=1,180
!         do ilat=1,181
!          read(4,*)real4(ilat,ilon)
!          util3(ilat,ilon)=1.D0*real4(ilat,ilon)
!         enddo
!        enddo
!        !!call refinp(xpivo,xpivn,util3,ipcrs,180,util1,idm,180)
!
!      if (ichan.eq.1) then
!        do i=1,iold
!        do j=1,jdm 
!        ipo(i,j)=0.
!        if (abs(util3(i,j)-spval).gt.2.) ipo(i,j)=1    ! define old ip
!        enddo
!        enddo
!      endif
!      call refinp(equato,equatn,ipo,util3,util1)
!
!c$OMP PARALLEL DO
!        do ilon=1,jdm
!         do ilat=1,idm
!          rno3(ilat,ilon,ichan)=util1(ilat,ilon)
!         enddo
!        enddo
!c$OMP END PARALLEL DO
!
!      enddo !ichan
!      close(4)
!      write(*,'(a,3i4,e12.4)')
!     . 'BIO: bioinit: rno3 at test point',
!     .    itest,jtest,1,rno3(itest,jtest,1)
!
!      open(4,file=flnmbio2,status='unknown')
!      do ichan=1,kgrd
!        do ilon=1,180
!         do ilat=1,181
!          read(4,*)real4(ilat,ilon)
!          util3(ilat,ilon)=1.D0*real4(ilat,ilon)
!         enddo
!        enddo
!        !!call refinp(xpivo,xpivn,util3,ipcrs,180,util1,idm,180)
!      if (ichan.eq.1) then
!        do i=1,iold
!        do j=1,jdm
!        ipo(i,j)=0.
!        if (abs(util3(i,j)-spval).gt.2.) ipo(i,j)=1    ! define old ip
!        enddo
!        enddo
!      endif
!      call refinp(equato,equatn,ipo,util3,util1)
!
!c$OMP PARALLEL DO
!        do ilon=1,jdm
!         do ilat=1,idm
!          so2(ilat,ilon,ichan)=util1(ilat,ilon)
!         enddo
!        enddo
!c$OMP END PARALLEL DO
!
!      enddo !ichan
!      close(4)
!      write(*,'(a,3i4,e12.4)')
!     .'BIO: bioinit: so2 at test point',itest,jtest,1,so2(itest,jtest,1)
!
!      ! save out rno3 and so2 fields
!      !save out in ascii
!      print*,'saving out in ',flnmovt(1:lgth)//'rno3_so2_init.dat'
!      open(unit=120,file=flnmovt(1:lgth)//'rno3_so2_init.dat',
!     .     status='unknown')
!      rno3min= 1.e30
!      rno3max=-1.e30
!      so2min = 1.e30
!      so2max =-1.e30
!      do ichan=1,kgrd
!       do j=1,jj			!  do not parallelize
!        do l=1,isp(j)
!         do i=ifp(j,l),ilp(j,l)
!          !protect against non-existing values
!          if (rno3(i,j,ichan).lt.0.)rno3(i,j,ichan)=0.
!          if ( so2(i,j,ichan).lt.0.) so2(i,j,ichan)=0.
!
!          rno3min=min(rno3min,rno3(i,j,ichan))
!          rno3max=max(rno3max,rno3(i,j,ichan))
!          so2min =min(so2min,so2(i,j,ichan))
!          so2max =max(so2max,so2(i,j,ichan))
!          write(120,'(3i4,2e12.4)')
!     .     i,j,ichan,rno3(i,j,ichan),so2(i,j,ichan)
!         enddo
!        enddo
!       enddo
!      enddo
!      close(120)
!
!      write(*,'(a,4e12.4)')'BIO: bioinit: rno3,so2 min-max=',
!     .       rno3min,rno3max,so2min,so2max
!
!      !read GLODAP data for dic
!      open(4,file=flnmbio3,status='unknown')
!      do ichan=1,kgrd
!        do ilon=1,180
!         do ilat=1,181
!          read(4,*)real4(ilat,ilon)
!          util3(ilat,ilon)=real4(ilat,ilon)
!         enddo
!        enddo
!        !!call refinp(xpivo,xpivn,util3,ipcrs,180,util1,idm,180)
!        if (ichan.eq.1) then
!          do i=1,iold
!          do j=1,jdm
!          ipo(i,j)=0.
!          if (abs(util3(i,j)-spval).gt.2.) ipo(i,j)=1    ! define old ip
!          enddo
!          enddo
!        endif
!        call refinp(equato,equatn,ipo,util3,util1)
!
!c$OMP PARALLEL DO
!        do ilon=1,jdm
!         do ilat=1,idm
!          dic(ilat,ilon,ichan)=util1(ilat,ilon)
!         enddo
!        enddo
!c$OMP END PARALLEL DO
!
!      enddo !ichan
!      close(4)
!      write(*,'(a,3i4,e12.4)')
!     . 'BIO: bioinit: dic at test point',itest,jtest,1
!     . ,dic(itest,jtest,1)
!
!! save out dic field
!      write (lp,'(a/9x,a)') 'BIO: save dic init field',
!     .   flnmovt(1:lgth)//'dic_init.dat'
!
!      !writeout in ascii
!      print*,'saving out in ',flnmovt(1:lgth)//'dic_init.dat'
!      open(unit=120,file=flnmovt(1:lgth)//'dic_init.dat',
!     .     status='unknown')
!      dicmin = 1.e30
!      dicmax =-1.e30
!      do ichan=1,kgrd
!       do j=1,jj                        !  do not parallelize
!        do l=1,isp(j)
!         do i=ifp(j,l),ilp(j,l)
!          !protect against non-existing values
!          if (dic(i,j,ichan).lt.0.) dic(i,j,ichan)=0.
!          dicmin =min(dicmin,dic(i,j,ichan))
!          dicmax =max(dicmax,dic(i,j,ichan))
!          write(120,'(3i4,e12.4)')i,j,ichan,dic(i,j,ichan)
!         enddo
!        enddo
!       enddo
!      enddo
!      close(120)
!
!      write(*,'(a,2e12.4)')'BIO: bioinit: dic min-max=',
!     .       dicmin,dicmax
!
      !save depths
      print*,'saving out in dp_ini.asc'
      call openunit('dp_ini.asc',iu_dpinit)
       do j=1,jj                        !  do not parallelize
       do l=1,isp(j)
       do i=ifp(j,l),ilp(j,l)
       do k=1,kdm
        write(iu_dpinit,'(4i5,e12.4)')
     .       nstep,i,j,k,dpinit(i,j,k)/onem
       enddo
       enddo
       enddo
       enddo
       call closeunit(iu_dpinit)
 
 
       else   !save4ini is false
!!these rno3 and so2_init fields are not correct. There are void points due to
!mismatch of the noaa grid and the hycom grid. To fill in have to do
!the interpolation in matlab (furtuna).
!at the same time use dps and interpolate to layer depths from the model
!     file_dir='/explore/nobackup/aromanou/RUNS/'
      file_dir='/gpfsm/dnb1/aromanou/'
!     call openunit(file_dir//'rno3_correctini.asc',iu_bioinit)
      call openunit(file_dir//'nitoa195x180_20w.asc',iu_bioinit)
      do ichan=1,kdm
       do j=1,jdm
        do i=1,idm
         read(iu_bioinit,'(f12.8)')tracer(i,j,ichan,1)
          if(tracer(i,j,ichan,1).le.0.)tracer(i,j,ichan,1)=0.085
        enddo
       enddo
      enddo
      call closeunit(iu_bioinit)

!     call openunit(file_dir//'so2_correctini.asc',iu_bioinit)
      call openunit(file_dir//'siloa195x180_20w.asc',iu_bioinit)
      do ichan=1,kdm
       do j=1,jdm
        do i=1,idm
         read(iu_bioinit,'(f12.8)')tracer(i,j,ichan,3)
          if(tracer(i,j,ichan,3).le.0.)tracer(i,j,ichan,3)=0.297
        enddo
       enddo
      enddo
      call closeunit(iu_bioinit)
     
!     call openunit(file_dir//'dic_correctini.asc',iu_bioinit)
      call openunit(file_dir//'dic195x180_20w.asc',iu_bioinit)
      dicmin = 1.e30
      dicmax =-1.e30
      do ichan=1,kdm
       do j=1,jdm
        do i=1,idm
         read(iu_bioinit,'(f12.8)')dic(i,j,ichan)
          if (dic(i,j,ichan).le.0.) dic(i,j,ichan)=1837.
          dic(i,j,ichan)=max(dic(i,j,ichan),1837.)   !set minimum =1837
          dicmin =min(dicmin,dic(i,j,ichan))
          dicmax =max(dicmax,dic(i,j,ichan))
        enddo
       enddo
      enddo
      call closeunit(iu_bioinit)
      write(*,'(a,2e12.4)')'BIO: bioinit: dic min-max=',
     .       dicmin,dicmax

      endif    !save4ini

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

      dictot = 0.0
      idic = 0

! the following loop has parallelization problems....

      do 1000 j=1,jj
      do 1000 l=1,isp(j)
      do 1000 i=ifp(j,l),ilp(j,l)

        do k=1,kdm
          !Nitrate
          !read earlier from file

          !Ammonium
          tracer(i,j,k,2) = 0.5
          if (dpinit(i,j,k)/onem .gt. 4000.) 
     .                      tracer(i,j,k,2) = Pdeep(2)

          !Silica
          !read earlier from file
!         if (dpinit(i,j,k)/onem .gt. 4000.) 
!    .                      tracer(i,j,k,2) = Pdeep(3)

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
          if (dpinit(i,j,k)/onem .gt. 4000.) 
     .                      tracer(i,j,k,4) = Pdeep(4)

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


cdiag do k=1,kdm
cdiag  write(998,'(4i5,4e12.4)')
cdiag.       nstep,i,j,k,dpinit(i,j,k)/onem,
cdiag.       tracer(i,j,k,1),tracer(i,j,k,3),dicmod(i,j,k)
cdiag enddo

 1000 continue


c     dicavg = dictot/float(idic)
c     write(6,*)'avg deep DIC = ',dicavg
 
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
      startdic = 2054.0
      nflg = 0

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

c  Light saturation data
      write(6,*)'Light saturation data...'
      avgq = 0.0

      do j=1,jj
       do k=1,kdm
        do l=1,isp(j)
         do i=ifp(j,l),ilp(j,l)
          avgq(i,j,k) = 25.0
         enddo
        enddo
       enddo
      enddo

c  Coccolithophore max growth rate
      write(6,*)'Coccolithophore max growth rate...'
      do j=1,jj
       do k=1,kdm
        do l=1,isp(j)
         do i=ifp(j,l),ilp(j,l)
          gcmax(i,j,k) = 0.0
         enddo
        enddo
       enddo
      enddo
 
      !save initialization
      do nt=1,ntyp+n_inert+ndet+ncar
!       ntunit=90+nt
        ntchar='00'
        if(nt.le.9)write(ntchar,'(i1)')nt
        if(nt.gt.9)write(ntchar,'(i2)')nt
        print*,'BIO: saving initial tracer fields '
     .        ,'bioinit_tracer'//ntchar
!       open(unit=ntunit,file=flnmovt(1:lgth)//'bioinit_tracer'//ntchar,
!    .       status='unknown')


        call openunit('bioinit_tracer'//ntchar,iu_bioinit)

        do k=1,kdm
        do j=1,jj				!  do not parallelize
        do l=1,isp(j)
        do i=ifp(j,l),ilp(j,l)
           write(iu_bioinit,'(3i4,2e12.4)')
     .           i,j,k,dpinit(i,j,k)/onem,tracer(i,j,k,nt)
        enddo
        enddo
        enddo
        enddo
      call closeunit(iu_bioinit)
      enddo

      if (save4ini) then
      print*, '   '
      stop 'WROTE OUT INI FIELDS. RESTART USING save4ini=.false.'
      endif

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

      implicit none
#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"


      integer iant,isin,ispc,isat,iein,iepc,ieat,incp
     .       ,inca,inat,imed,inin,inpc,nr,ntot

      real rlat,rlon

      integer ir,nir
      common /bir/ ir(idm,jdm),nir(nrg)

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
