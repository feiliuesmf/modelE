      MODULE obio_incom

! parameters and arrays neccessary for obio_init and obio_bioinit

      USE obio_dim

      implicit none

      SAVE



      real :: rmumax                !max phyto growth rate at 20oC, d/
      real :: obio_wsd,obio_wsh     !phyto sinking rate m/d, m/h
      real :: rik                   !light saturation parameter umol quanta/m2/s
      common /bphy/ rmumax(nchl),obio_wsd(nchl),
     .              obio_wsh(nchl),rik(3,nchl)

      real :: rkn,rks     !half-saturation constants for nitrogen, silica (uM)
      real :: rkf         !half-saturation constant for iron (nM)
      common /brk/ rkn(nchl),rks(nchl),rkf(nchl)

      real rad, pi2
      common /bcnst/ rad,pi2

      real :: pHsfc,pHmin,pHmax   !pH at surface,minimun for iteration,
                                  !max for iteration
      common /cpH/pHsfc,pHmin,pHmax

      real Pdeep,detdeep,cardeep
      common /bpdp/ Pdeep(ntyp)            !deep BC
      common /bddp/ detdeep(ndet)          !detrital deep BC
      common /bcdp/ cardeep(ncar)          !carbon deep BC

      real :: cchl        !C:chl ratio g:g for 3 light states
      real :: cchld       !mean C:chl ratio for the day
      common /bcchl/ cchl(3),cchld

      real :: cnratio     !C:N ratio
      real :: csratio     !C:Si ratio
      real :: cfratio     !C:Fe ratio
      common /beratio/ cnratio,csratio,cfratio

      real :: wsdeth      !sinking rate of detritus (m/d)
      real :: remin       !detrital remineralization rate /d
      common /bkwsdet/ wsdeth(ndet),remin(ndet)

      real :: Fescavrate  !scavenging rate for dissolved iron
      common /bfescav/ Fescavrate(2)

C if CARBON == 1
      real, parameter :: excp=0.05              !excretion of DOC by phyto %growth
      real, parameter :: resp=0.05              !respiration of DIC by phyto %growth
      real, parameter :: excz=0.05/24.0         !excretion of DOC by zoopl/hr
      real, parameter :: resz=0.05/24.0         !respiration of DIC by zoopl/hr
      real, parameter :: phygross=1.0+excp+resp !factor to derive gross PP!

      real, parameter :: rlamdoc=0.017/24.0    !N-dependent DOC remin/hr
      real, parameter :: rkdoc1=0.3*10.0       !N-dep half sat uM(PO4) modified
                                               !for nitrate by mult*10.0,
                                               !based on Conkright et al. 1994
      real, parameter :: rkdoc2=15.0           !DOC-dep half-sat uM(C)
      real, parameter :: rlampoc=0.05/24.0     !detrital breakdown/hr
      real, parameter :: uMtomgm3=12.0         !conversion uM to mg/m3
      real, parameter :: Pzo=1.0*uMtomgm3/50.0 !zoopl half-sat for
                                               !DOC excretion mg/m3(chl,assuming
                                               !C:chl ratio of 50))
      real, parameter :: awan=0.337/(3.6E+5)   !piston vel coeff., from
                                               !Wanninkof 1992, but adjusted
                                               !by OCMIP, and converted from
                                               !cm/hr to m/s
      real, parameter :: stdslp=1013.25        !standard sea level pressure in mb

      real, parameter :: Rm=1.20/24.0          !max zoopl. growth rate/hr
                                               !increase to account for excretion
                                               !and respiration
C if CARBON /=1    parameter(Rm=1.0/24.0)      !max zoopl. growth rate/hr


      integer lam               !wavelength in nm
      real facirr               !array of factors to compute mean irradiance w/in water column
      real aw,bw                !absorption,scattering coefficients of water
      real ac,bc                !absorption and scattering coefficients of chlorophyll
      common /bfac/ facirr(nh,nch,5,ncd)  
      common /bwat/ aw(nlt),bw(nlt),lam(nlt)  
      common /bopt/ ac(nchl,nlt),bc(nchl,nlt)

      real, parameter :: solFe=0.02    !solubility of iron

c     parameter(bn=0.5,bs=0.5)        !N/chl and Si/chl ratios

!save from moved ifst parts
      integer nl450
      real excdom,bbw,Dmax,rd,ru,rmus,rmuu,rn,roair,wfac
       integer, parameter :: it0inc=1,nt0=80/it0inc,isalinc=1,nsal=20
       integer, parameter :: idicinc=2,ndic=(650+idicinc)/idicinc
       integer, parameter :: itainc=2,nta=(500+itainc)/itainc


       !pco2 data file
       real, ALLOCATABLE :: pco2tab(:,:,:,:)
       !real pco2tab
       !common /bpco2tab/pco2tab(nt0,nsal,ndic,nta)

      common/exifst1/nl450
      common/exifst2/excdom(nlt),bbw,Dmax,rd,ru,rmus,rmuu
     .             ,rn,roair,wfac(nlt)

      END MODULE obio_incom
