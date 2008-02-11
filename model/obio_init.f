c ----------------------------------------------------------------
      subroutine obio_init
c --- biological/light setup
c ----------------------------------------------------------------
c 
      USE FILEMANAGER, only: openunit,closeunit

      USE obio_dim
      USE obio_incom
      USE obio_forc, only : ihra,Eda,Esa,atmFe_all,alk
      USE obio_com, only : npst,npnd,WtoQ,obio_ws,P_tend,D_tend
     .                    ,C_tend,wsdet,gro


      implicit none  

#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"

      integer iu_bio
      integer nt,nl
      integer imon,ihr,nrec,ichan
      integer np,lambda,ic
      integer icd,ntr,ich,ih,iu_fac
      real saw,sbw,sac,sbc
      real*4  facirr4(nh,nch,5,ncd)


      real planck,c,hc,oavo,rlamm,rlam450,Sdom,rlam,hcoavo
     .    ,rnn,rbot,t,tlog,fac,a0,a1,a2,a3,b0,b1,b2,b3

      character*50 title
!     character*50 cfle
      character cacbc*11,cabw*10

      data cacbc,cabw /'acbc25b.dat','abw25b.dat'/

      data a0,a1,a2,a3 /0.9976,0.2194,5.554E-2,6.7E-3/
      data b0,b1,b2,b3 /5.026,-0.01138,9.552E-6,-2.698E-9/

c 
      print*, 'Ocean Biology setup starts'

c  Read in constants, light data
c  Computes constants over entire run of model, reads in required
c  data files, and otherwise obtains one-time-only information
c  necessary for the run.
 
c  Degrees to radians conversion
      pi = dacos(-1.0D0)
      pi2 = pi*2.0
      rad = 180.0D0/pi

      do nt = 1,nchl
       rkn(nt) = 0.0
       rks(nt) = 0.0
       rkf(nt) = 0.0
      enddo
      pHsfc = 8.0
      pHmin = 7.5
      pHmax = 8.6
c
c  Phytoplankton group parameters
      do nt = 1,nchl
       obio_wsd(nt)    = 0.0
       obio_wsh(nt)    = 0.0
      enddo
      do nt = 1,nchl
       rmumax(nt) = 0.0
        rik(1,nt) = 0.0
        rik(2,nt) = 0.0
        rik(3,nt) = 0.0
      enddo
      Pdeep(1) = 32.0    !bottom BC for nitrate
      Pdeep(2) = 0.1     !ammonium
      Pdeep(3) = 60.0    !silica
      Pdeep(4) = 0.6     !iron from Archer and Johnson 2000
      do nt = nnut+1,ntyp
       Pdeep(nt) = 0.0    !chl and herbivores
      enddo
      do nt = 1,ndet
       detdeep(nt) = 0.0 !detritus
      enddo
      cardeep(1) = 0.0   !DOC
      cardeep(2) = 2330.0  !DIC uM(C) from Goyet et al (2000)
c
c  Carbon:chl ratios for different adaptation states
      cchl(1) = 25.0
      cchl(2) = 50.0
      cchl(3) = 80.0
c      cchl(1) = 20.0
c      cchl(2) = 60.0
c      cchl(3) = 100.0
      cnratio = 106.0/16.0*12.0    !C:N ratio (ugl:uM)
      csratio = 106.0/16.0*12.0    !C:Si ratio (ugl:uM)
      cfratio = 150000.0*12.0*1.0E-3    !C:Fe ratio (ugl:nM)
c
!!#if NCHL_DEFINED > 0
      if (nchl > 0) then
c  Diatoms
      nt = 1
      rmumax(nt) = 1.50       !u max in /day at 20C
      obio_wsd(nt)    = 0.75  !sinking rate in m/day
      obio_wsd(nt)    = 0.50  !sinking rate in m/day   !!change Oct27,2008
      rik(1,nt)  = 90.0       !low light-adapted Ik (<50 uE/m2/s)
      rik(2,nt)  = 93.0       !medium light-adapted Ik (50-200 uE/m2/s)
      rik(3,nt)  = 184.0      !high light adapted Ik (>200 uE/m2/s)
      rkn(nt) = 1.0           !M-M half-sat constant for nitrogen
      rks(nt) = 0.2           !M-M half-sat constant for silica
      rkf(nt) = 0.12          !M-M half-sat constant for iron

      endif
!!#endif
!!#if NCHL_DEFINED > 1
      if (nchl > 1) then
c  Chlorophytes
      nt = 2
      rmumax(nt) = rmumax(nt-1)*0.840
      obio_wsd(nt)    = 0.25
      rik(1,nt)  = rik(1,nt-1)*1.077
      rik(2,nt)  = rik(2,nt-1)*0.935
      rik(3,nt)  = rik(3,nt-1)*0.781
      rkn(nt) = rkn(nt-1)*0.75
      rkn(nt) = rkn(nt-1)*0.667   !1/3 distance bet. cocco and dia
      rkf(nt) = rkf(nt-1)*0.835   !midway between cocco's and diatoms
      rkf(nt) = rkf(nt-1)*0.779   !1/3 distance bet. cocco and dia

      endif
!!#endif
!!#if NCHL_DEFINED > 2
      if (nchl > 2) then
c  Cyanobacteria
      nt = 3
      rmumax(nt) = rmumax(nt-2)*0.670
      obio_wsd(nt)    = 0.0085
      rik(1,nt)  = rik(1,nt-2)*0.723
      rik(2,nt)  = rik(2,nt-2)*0.710
      rik(3,nt)  = rik(3,nt-2)*0.256
      rkn(nt) = rkn(nt-2)*0.50
      rkf(nt) = rkf(nt-2)*0.67  !equals cocco

      endif
!!#endif
!!#if NCHL_DEFINED > 3
      if (nchl > 3) then
c  Coccolithophores
      nt = 4
c      rmumax(nt) = rmumax(nt-3)*0.663   !all coccos
       rmumax(nt) = rmumax(nt-3)*0.755   !E. huxleyi only
c      rmumax(nt) = rmumax(nt-3)*0.781   !E. huxleyi only (no Sunda/Hunts)
      obio_wsd(nt)    = 0.82
      obio_wsd(nt)    = 0.648
      rik(1,nt)  = rik(1,nt-3)*0.623
      rik(2,nt)  = rik(2,nt-3)*0.766
      rik(3,nt)  = rik(3,nt-3)*0.899
      rkn(nt) = rkn(nt-3)*0.5
      rkf(nt) = rkf(nt-3)*0.67

      endif
!!#endif
!!#if NCHL_DEFINED > 4
      if (nchl > 4) then
c  Dinoflagellates
      nt = 5
      rmumax(nt) = rmumax(nt-4)*0.335
      obio_wsd(nt)    = 0.0
      rik(1,nt)  = rik(1,nt-4)*1.321
      rik(2,nt)  = rik(2,nt-4)*1.381
      rik(3,nt)  = rik(3,nt-4)*1.463
      rkn(nt) = rkn(nt-4)*1.0
      rkf(nt) = rkf(nt-4)*0.67

      endif
!!#endif
      do nt = 1,nchl
       obio_wsh(nt) = obio_wsd(nt)/24.0  !convert to m/hr
      enddo
c  Detrital sinking rates m/h
      wsdeth(1) = 30.0/24.0     !nitrogen
      wsdeth(2) = 50.0/24.0     !silica
      wsdeth(3) = 20.0/24.0     !iron
c
c  Detrital remineralization rates /hr
      remin(1) = 0.010/24.0            !nitrogen
      remin(2) = 0.0001/24.0           !silica
      remin(3) = 0.020/24.0            !iron
      fescavrate(1) = 2.74E-5/24.0      !low fe scavenging rate/hr
      fescavrate(2) = 50.0*fescavrate(1) !high fe scavenging rate/hr
c
c  (originally done inside lidata subroutine of obio_daysetrad)
c  Reads in radiative transfer data: specifically
c  water data (seawater absorption and total scattering coefficients,
c  and chl-specific absorption and total scattering data for
c  several phytoplankton groups).  PAR (350-700) begins at index 3,
c  and ends at index 17.
c     
c  Water data files
!     cfle = cabw                       
!     open(4,file='/explore/nobackup/aromanou/2.0deg/'//cfle
!    .      ,status='old',form='formatted')

      call openunit('cfle1',iu_bio)
      do ic = 1,6
       read(iu_bio,'(a50)')title
      enddo
      np = 0    
      do nl = 1,nlt
       read(iu_bio,20)lambda,saw,sbw
       lam(nl) = lambda
       aw(nl) = saw
       bw(nl) = sbw
      enddo
      call closeunit(iu_bio)
 20   format(i5,f15.4,f10.4)

c  Phytoplankton group chl-specific absorption and total scattering
c  data.  Chl-specific absorption data is normalized to 440 nm; convert
c  here to actual ac*(440)
!     cfle = cacbc
!     open(4,file='/explore/nobackup/aromanou/2.0deg/'//cfle
!    .      ,status='old',form='formatted')
      call openunit('cfle2',iu_bio)
      do ic = 1,6
       read(iu_bio,'(a50)')title
      enddo
      do nt = 1,nchl
       read(iu_bio,'(a50)')title
       np = 0
       do nl = 1,19
        read(iu_bio,30)lambda,sac,sbc
        ac(nt,nl) = sac
        bc(nt,nl) = sbc
       enddo
       do nl = 20,nlt
        ac(nt,nl) = 0.0
        bc(nt,nl) = 0.0
       enddo
      enddo
      call closeunit(iu_bio)
 30   format(i4,2f10.4)


!ifst part from daysetrad.f
c      h = 6.6256E-34   !Plancks constant J sec
       planck = 6.6256E-34   !Plancks constant J sec
       c = 2.998E8      !speed of light m/sec
c      hc = 1.0/(h*c)
       hc = 1.0/(planck*c)
       oavo = 1.0/6.023E23   ! 1/Avogadros number
       hcoavo = hc*oavo
       do nl = npst,npnd
        rlamm = float(lam(nl))*1.0E-9  !lambda in m
        WtoQ(nl) = rlamm*hcoavo        !Watts to quanta conversion
       enddo
       !CDOM absorption exponent
       rlam450 = 450.0
       Sdom = 0.014
       do nl = 1,nlt
        if (lam(nl) .eq. 450)nl450 = nl
        rlam = float(lam(nl))
        excdom(nl) = exp(-Sdom*(rlam-rlam450))
       enddo
       if (nl450.eq.0) stop 'obio_init: nl450=0'
       !First time thru set ihra to 1 to assure correct value read in
       !from restart file
       do j=1,jj
       do l=1,isp(j)
       do i=ifp(j,l),ilp(j,l)
        ihra(i,j) = 1
       enddo
       enddo
       enddo

!ifst part from edeu.f
       bbw = 0.5            !backscattering to forward scattering ratio
       rmus = 1.0/0.83      !avg cosine diffuse down
       Dmax = 500.0         !depth at which Ed = 0

       rnn = 1.341
       rmuu = 1.0/0.4             !avg cosine diffuse up
       rbot = 0.0                 !bottom reflectance
       rd = 1.5   !these are taken from Ackleson, et al. 1994 (JGR)
       ru = 3.0

c  Read in factors to compute average irradiance
! (this part originally done inside obio_edeu)

      print*,'Reading factors for mean irradiance at depth...'
      print*,'nh,nch,ncd=',nh,nch,ncd

      call openunit('facirr',iu_fac)
      do icd=1,ncd
       do ntr=1,5
        do ich=1,nch
         do  ih=1,nh
          read(iu_fac,*)facirr4(ih,ich,ntr,icd)
           facirr(ih,ich,ntr,icd)=1.D0*facirr4(ih,ich,ntr,icd)
         enddo
        enddo
       enddo
      enddo
      call closeunit(iu_fac)
      print*,'nstep, facirr(10,18,1)=',nstep,facirr4(10,18,1,1)
      print*,'nstep, facirr(10,18,1)=',nstep,facirr(10,18,1,1)


!ifst part from ocalbedo.f
       rn = 1.341        !index of refraction of pure seawater
       roair = 1.2E3     !density of air g/m3
       do nl = 1,nlt
        if (lam(nl) .lt. 900)then
         t = exp(-(aw(nl)+0.5*bw(nl)))
         tlog = alog(1.0E-36+t)
         fac = a0 + a1*tlog + a2*tlog*tlog + a3*tlog*tlog*tlog
         wfac(nl) = min(fac,1.0)
         wfac(nl) = max(fac,0.0)
        else
         fac = b0 + b1*rlam + b2*rlam*rlam + b3*rlam*rlam*rlam
         wfac(nl) = max(fac,0.0)
        endif
       enddo

!ifst part from ptend.f
       do k=1,kdm

         do nt=1,nchl
          obio_ws(k,nt) = 0.0
          gro(k,nt) = 0.0
         enddo
 
         do nt=1,ntyp+n_inert
          P_tend(k,nt) = 0.0
         enddo
 
         do nt=1,ndet
          D_tend(k,nt) = 0.0
          wsdet(k,nt) = 0.0
         enddo
 
         do nt=1,ncar
          C_tend(k,nt) = 0.0
         enddo
       enddo
       do nt=1,nchl
        obio_ws(kdm+1,nt)=0.0
       enddo
       do nt=1,ndet
        wsdet(kdm+1,nt) = 0.0
       enddo
 
    
!ifst part from ppco2tab.f
       ALLOCATE (pco2tab(nt0,nsal,ndic,nta))

!      open(4,file='/explore/nobackup/aromanou/pco2.tbl.asc'
!    .       ,status='old')

       call openunit('pco2table',iu_bio)
       print*, 'obio_init, pco2tbl: ',nta,ndic,nsal,nt0
       do nl=1,nta
        do k=1,ndic
         do j=1,nsal
          do i=1,nt0
           read(iu_bio,'(e12.4)')pco2tab(i,j,k,nl)
          enddo
         enddo
        enddo
       enddo
       call closeunit(iu_bio)
       print*,'BIO: read pCO2 table: ',
     .        pco2tab(1,1,1,1),pco2tab(50,10,100,100)

!read in light (this will be changed later to be passed from atmosphere
      print*, '    '
      print*, 'reading OASIM data.....'

      ALLOCATE (Eda(idm,jdm,nlt,nhn,12),Esa(idm,jdm,nlt,nhn,12))

      open(unit=iu_bio,file='oasimdirect'
     . ,form='unformatted',status='old',access='direct' 
     . ,recl=idm*jdm*8/4)
      nrec=0
      do imon=1,12
      do ihr=1,nhn
      do ichan=1,nlt
        nrec=nrec+1
        read (iu_bio,rec=nrec)
     .       ((Eda(i,j,ichan,ihr,imon),i=1,idm),j=1,jdm)

        nrec=nrec+1
        read (iu_bio,rec=nrec)
     .       ((Esa(i,j,ichan,ihr,imon),i=1,idm),j=1,jdm)
      enddo
      enddo
      enddo
      close(iu_bio)



!read in atmospheric iron deposition (this will also be changed later...)
      print*, '    '
      print*, 'reading iron data.....'

      if (IRON_from.eq.0) then
      open(unit=iu_bio,file='atmFedirect0'
     . ,form='unformatted',status='old',access='direct' 
     . ,recl=idm*jdm*8/4)
      do imon=1,12  !1 year of monthly values
       nrec=imon
       read (iu_bio,rec=nrec)((atmFe_all(i,j,imon),i=1,idm),j=1,jdm)
      enddo
      close(iu_bio)
      endif

      if (IRON_from.eq.1) then
      call openunit('atmFedirect1',iu_bio)
      do imon=1,12  !1 year of monthly values
       do j=1,jdm
        do i=1,idm
        read(iu_bio,'(e12.4)')atmFe_all(i,j,imon)
        enddo
       enddo
      enddo
      call closeunit(iu_bio)
      endif

!read in alkalinity annual mean file
      if (ALK_CLIM.eq.1) then      !read from climatology
      call openunit('alkalindirect',iu_bio)
      do k=1,kdm
       do j=1,jdm
        do i=1,idm
         read(iu_bio,'(f12.8)')alk(i,j,k)
        enddo
       enddo
      enddo
      call closeunit(iu_bio)
      else      !set to zero, obio_carbon sets alk=tabar*sal/sal_mean
      alk = 0.
      endif

! printout some key information
      write(*,*)'**************************************************'
      write(*,*)'**************************************************'
      write(*,*)'**************************************************'
      write(*,*)'           INITIALIZATION                         '

      write(*,*)'ALK_CLIM = ', ALK_CLIM
      if (ALK_CLIM.eq.0) write(*,*) 'ALKALINITY, from SALINITY'
      if (ALK_CLIM.eq.1) write(*,*) 'ALKLNTY, GLODAP annmean'
      if (ALK_CLIM.eq.2) write(*,*) 'ALKALINITY prognostic'

      write(*,*)'IRON_from = ', IRON_from
      if (IRON_from.eq.0) write(*,*) 'IRON FLUXES, GOCART model'
      if (IRON_from.eq.1) write(*,*) 'IRON FLUXES, R.Miller dustfluxes'

      write(*,*)'**************************************************'
      write(*,*)'**************************************************'
      write(*,*)'**************************************************'

      return
      end
