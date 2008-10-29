#include "rundeck_opts.h"
      module HYCOM_SCALARS
      USE HYCOM_DIM_GLOB

      implicit none

      private

c ---  s w i t c h e s    (if set to .true., then...)
c --- diagno      output model fields and diagnostic messages
c --- thermo      use thermodynamic forcing functions
c --- windf       include wind stress in forcing functions
c --- relax       activate lateral boundary nudging
c --- trcout      advect tracer and save results in history/restart file
c --- dotrcr      perform column physics operations on tracer array(s)
c
      logical, public:: diagno,thermo,windf,relax,trcout,dotrcr
!!      common/swtchs/diagno,thermo,windf,relax,trcout,dotrcr


c     real taux,tauy,wndspd,radflx,airtmp,precip,vapmix,freshw,diafor
c    .    ,pwall,swall,twall
!!      real taux,tauy,oice,oemnp,ustar,ustarb,oflxa2o,osalt
!!      real freshw,diafor
c
!!      common/varbls/time,time0,delt1,dlt,w0,w1,w2,w3,ws0,ws1,ws2,ws3
!!     .     ,area,avgbot,watcum,empcum,slfcum,sala2o,tavini
      real, public :: time,time0,delt1,dlt,w0,w1,w2,w3,ws0,ws1,ws2,ws3,
     .     area,avgbot,watcum,empcum,slfcum,sala2o,tavini

!!      common/varbl2/nstep,nstep0,nstepi,lstep,l0,l1,l2,l3,ls0,ls1
!!     .             ,ls2,ls3,oddev
      integer, public ::  nstep,nstep0,nstepi,lstep,l0,l1,l2,l3,ls0,ls1
     .             ,ls2,ls3,oddev
c
c --- 'baclin' = baroclinic time step
c --- 'batrop' = barotropic time step
c --- 'thkdff' = diffusion velocity (cm/s) for thickness diffusion
c --- 'veldff' = diffusion velocity (cm/s) for momentum dissipation
c --- 'temdff' = diffusion velocity (cm/s) for temp/salin. mixing
c --  'viscos' is nondimensional, used in deformation-dependent viscosity
c --- 'diapyc' = diapycnal diffusivity times buoyancy freq. (cm**2/s**2)
c --- 'trcfrq' = number of time steps between tracer transport calculations
c --- 'h1'     = depth interval used in lateral weighting of hor.pres.grad.
c --- slip = +1  for free-slip boundary cond., slip = -1  for non-slip cond.
c --- 'cbar'   = rms flow speed (cm/s) for linear bottom friction law
c --- 'diagfq' = number of days between model diagnostics (incl.output)
c --- 'ntracr' = number of time steps between tracer transport
c --- 'wuv1/2' = weights for time smoothing of u,v field
c --- 'wts1/2' = weights for time smoothing of t,s field
c --- 'wbaro'  = weight for time smoothing of barotropic u,v,p field
c --- 'thkmin' = minimum mixed-layer thickness
c --- 'thkbot' = thickness of bottom boundary layer
c --- 'botmin' = minimum topo depth
c --- 'ekman'  = thickness of ekman layer
c --- 'sigjmp' = minimum density jump at mixed-layer bottom
c ---' salmin' = minimum salinity allowed in an isopycnic layer
c --- 'acurcy' = permissible roundoff error in column integral calc.
c --- 'nhr   ' = coupling freq. in hours
c
!!      common/parms1/thbase,theta(kdm),baclin,batrop,thkdff,
!!     .              veldff,temdff,viscos,diapyc,vertmx,h1,slip,cbar,
!!     .              diagfq,wuv1,wuv2,wts1,wts2,acurcy,wbaro,thkmin,
!!     .              thkbot,botmin,ekman,sigjmp,salmin(kdm)
      dimension theta(kdm),salmin(kdm)
      real, public ::
     &     theta,thbase,baclin,batrop,thkdff,veldff,temdff,viscos,
     .     diapyc,diapyn,vertmx,h1,slip,cbar,diagfq,wuv1,wuv2,wts1,wts2,
     &     acurcy, wbaro,thkmin,thkbot,botmin,ekman,sigjmp,salmin
c
!!      common/parms2/trcfrq,ntracr,nhr,mixfrq
      integer, public ::       trcfrq,ntracr,nhr,mixfrq
c
c --- 'tenm,onem,...' = pressure thickness values corresponding to 10m,1m,...
c --- 'g'      = gravity acceleration
c --- 'csubp'  = specific heat of air at constant pressure (j/g/deg)
c --- 'spcifh' = specific heat of sea water (j/g/deg)
c --- 'cd'     = drag coefficient
c --- 'ct'     = thermal transfer coefficient
c --- 'airdns' = air density at sea level (g/cm**3)
c --- 'evaplh' = latent heat of evaporation (j/g)
c --- 'thref'  = reference value of specific volume (cm**3/g)
c --- 'epsil'  = small nonzero number used to prevent division by zero
c
!!      common/consts/tenm,onem,tencm,onecm,onemm,g,csubp,spcifh,cd,ct,
!!     .              airdns,evaplh,thref,epsil,huge,radian,pi
c
      real, public ::
     &     tenm,onem,tencm,onecm,onemm,g,csubp,spcifh,cd,ct,airdns,
     .     evaplh,thref,epsil,huge,radian,pi
c
c --- grid point where detailed diagnostics are desired:
      common/testpt/itest,jtest
      integer, public :: itest,jtest
c
c ---      equatn   --  the i index of the equator
c
!!      common/grdparms/equatn
      real, public :: equatn
c
      character*60, public ::
     &     flnmdep,flnmrsi,flnmrso,flnmarc,flnmfor,flnmovt
     .            ,flnmini,flnmriv,flnmbas,flnmdia,flnmlat
     .            ,flnminp,flnmint,flnmins
     .            ,flnmcoso,flnmcosa,flnma2o,flnma2o_tau,flnmo2a
     .            ,flnmo2a_e,flnmo2a_n
!!      common/iovars/flnmdep,flnmrsi,flnmrso,flnmarc,flnmfor,flnmovt
!!     .            ,flnmini,flnmriv,flnmbas,flnmdia,flnmlat
!!     .            ,flnminp,flnmint,flnmins
!!     .            ,flnmcoso,flnmcosa,flnma2o,flnma2o_tau,flnmo2a
!!     .            ,flnmo2a_e,flnmo2a_n


c --- opening the bering strait requires information exchange across a
c --- 'u' face represented in 2 different locations in the tri-pole grid.
c --- the 2 locations are the northern tip of the bering 'inlet' (truncated
c --- bering channel) on the pacific side and the southern (i.e., upper)
c --- tip of the bering inlet in the panam pole patch
c
c --- ipacn,jpac:  grid point north of inlet head on pacific side
c --- ipacs,jpac:  grid point south of inlet head on pacific side
c --- iatln,jatl:  grid point north of inlet head on arctic ocean side
c --- iatls,jatl:  grid point south of inlet head on arctic ocean side
c
c --- thus, the pairs [(ipacn,jpac),(iatln,jatl)],[(ipacs,jpac),(iatls,jatl)]
c --- refer to identical grid cells in physical space.
c
      logical, public, parameter :: beropn=.true.	!  true if bering strait open
#ifdef HYCOM_RESOLUTION_2deg
      integer, public, parameter :: ipacn=67,ipacs=68,jpac= 95
      integer, public, parameter :: iatln= 2,iatls= 1,jatl=156
#endif
#ifdef HYCOM_RESOLUTION_1deg
      integer, public, parameter :: ipacn=137,ipacs=138,jpac=189
      integer, public, parameter :: iatln= 2,iatls= 1,jatl=312
#endif
c-----------------------------------------------------------------------------
c
!!      include 'dimensions.h'
!!!!      include 'common_blocks.h'
c
c --- layer densities (theta units):
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- s i g m a _ 0
ccc      data theta/
ccc     .   24.02, 24.70, 25.28, 25.77, 26.18, 26.52, 26.80, 27.03, 
ccc     .   27.22, 27.38, 27.52, 27.64, 27.74, 27.82, 27.88, 27.92/
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- s i g m a _ 2
c     data theta/30.90,31.87,32.75,33.54,34.24,34.85,35.37,35.80,
c    .           36.15,36.43,36.65,36.82,36.95,37.05,37.13,37.20/    !light
ccc   data theta/31.20,32.51,33.54,34.35,34.99,35.50,35.91,36.24,
ccc  .           36.50,36.70,36.85,36.96,37.04,37.10,37.15,37.20/    !medium
c     data theta/31.85,33.22,34.26,35.04,35.62,36.05,36.37,36.61,
c    .           36.79,36.92,37.01,37.07,37.11,37.14,37.17,37.20/    !heavy
#ifdef HYCOM_RESOLUTION_2deg
      data theta/
     .  28.50,29.71,30.82,31.83,32.74,33.55,34.26,34.87,35.38,35.80,
     .  36.14,36.41,36.62,36.78,36.90,36.99,37.07,37.16,37.28,37.45/  !vhv
#endif
#ifdef HYCOM_RESOLUTION_1deg
      data theta/
     .  28.89,30.07,31.11,32.02,32.81,33.49,34.07,34.56,34.97,35.31
     . ,35.59,35.82,36.01,36.17,36.31,36.44,36.56,36.67,36.77,36.86
     . ,36.94,37.01,37.07,37.12,37.16,37.20/          ! md
#endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c --- 'baclin' = baroclinic time step
c --- 'batrop' = barotropic time step
c --- 'diagfq' = number of days between model diagnostics (incl.output)
c --- 'equatn' = the i index of the equator
#ifdef HYCOM_RESOLUTION_2deg
      data baclin,batrop/3600.,100./,diagfq/365./          ! 2deg full global
      data equatn/122./
#endif
#ifdef HYCOM_RESOLUTION_1deg
      data baclin,batrop/1800., 50./,diagfq/365./          ! 1deg full global
      data equatn/243./
#endif
c
c --- 'thkdff' = diffusion velocity (m/s) for thickness diffusion
c --- 'veldff' = diffusion velocity (m/s) for momentum dissipation
c --- 'temdff' = diffusion velocity (m/s) for temp/salin. mixing
c --- 'viscos' is nondimensional, used in deformation-dependent viscosity
c --- 'vertmx' = scale velocity for vertical momentum mixing (m/s)
c     data thkdff/.005/,veldff/.1/,temdff/.02/,viscos/0.3/,vertmx/0./
      data thkdff/.10/,veldff/.1/,temdff/.02/,viscos/0.3/,vertmx/0./
c
c --- 'diapyc' = diapycnal diffusivity (m^2/s)
c --- 'diapyn' = diapycnal diffusivity times buoyancy freq. (m^2/s^2)
c --- 'h1'     = depth interval used in lateral weighting of hor.pres.grad.
c --- 'thkmin' = minimum mixed-layer thickness (m)
c --- 'acurcy' = permissible roundoff error in column integral calc.
      data diapyn/3.e-7/,diapyc/1.e-2/
      data h1/98060./,thkmin/5./,acurcy/1.e-11/,botmin/30./
c
c --- slip=+1  for free-slip boundary cond., slip=-1  for non-slip cond.
c --- 'cbar'   = rms flow speed (m/s) for linear bottom friction law
c --- 'thkbot' = thickness of bottom boundary layer (m)
c --- 'ekman'  = thickness of ekman layer (m)
c --- 'sigjmp' = minimum density jump at mixed-layer bottom (theta units)
      data slip/-1./,cbar/0.1/,thkbot/10./,ekman/30./,sigjmp/.01/
c
c --- weights for time smoothing
ccc      data wuv1,wuv2/.5,.25/
      data wuv1,wuv2/.75,.125/
ccc      data wts1,wts2/.875,.0625/
ccc      data wts1,wts2/.9375,.03125/
      data wts1,wts2/.96875,.015625/
      data wbaro/.125/
c
c --- layer thicknesses in units of pressure (kg/m/sec^2):
      data tenm,onem,tencm,onecm,onemm/98060.,9806.,980.6,98.06,9.806/
      data radian/57.2957795/,pi/3.1415926536/
c
c --- 'g'      = gravitational acceleration (m/s^2)
c --- 'csubp'  = specific heat of air at constant pressure (j/kg/deg)
c --- 'spcifh' = specific heat of sea water (j/kg/deg)
c --- 'airdns' = air density at sea level (kg/m^3)
c --- 'evaplh' = latent heat of evaporation (j/kg)
c --- 'thref'  = reference value of specific volume (m^3/kg)
c --- 'epsil'  = small nonzero number used to prevent division by zero
      data g/9.806/,csubp/1005.7/,spcifh/4185./
      data airdns/1.2/,evaplh/2.47e6/,thref/1.e-3/,epsil/1.e-11/
c
c --- 'itest,jtest' = grid point where detailed diagnostics are desired
      data itest,jtest/-1,-1/
c
c ---  s w i t c h e s    (if set to .true., then...)
c --- thermo      use thermodynamic forcing functions
c --- windf       use wind stress forcing function
c --- relax       activate lateral boundary nudging
c
#if defined(TRACERS_HYCOM_Ventilation) \
 || defined(TRACERS_GASEXCH_ocean) || defined(TRACERS_OceanBiology)
      data thermo/.true./, windf/.true./,relax/.false./,trcout/.true./
#else     
      data thermo/.true./, windf/.true./,relax/.false./,trcout/.false./
#endif
c
c
c --- use 'huge' to initialize array portions that the code should never access
      data huge/1.e33/
!!! temporary value for debug (for checksums to make sense)
!!!      data huge/0.e0/
      data nhr/1/                        ! couple every nhr hours
      data oddev/-1/
c
c --- i/o file names
c
c     flnmdep = name/location of basin depth array
c     flnmint = name/location of initial -t- field
c     flnmins = name/location of initial -s- field
c     flnminp = name/location of initial -p- field
c     flnmfor = name/location of the forcing functions
c     flnmdia = name/location of the diapycnal flux forcing fields
c     flnmriv = name/location of freshwater (river runoff) fields
c     flnmbas = name/location of basin mask file used in overtn diagno
c     flnmrsi = location (pathname) of restart file (input)
c     flnmrso = location (pathname) of restart file (output)
c     flnmarc = location (pathname) of archive files
c     flnmovt = location (pathname) of ovtn.xxxxxx files
c     flnmlat = location (pathname) of lat/lon at vorticity points
c
      data flnmlat    /'latlonij'/
      data flnmdep    /'hycomtopo'/
      data flnmint    /'temp_ini'/
      data flnmins    /'salt_ini'/
      data flnminp    /'pout_ini'/
      data flnmbas    /'ibasin'/
      data flnma2o    /'flxa2o'/
      data flnma2o_tau/'taua2o'/
      data flnmo2a    /'ssto2a'/
      data flnmo2a_e  /'e_o2a'/
      data flnmo2a_n  /'n_o2a'/
      data flnmcoso   /'cososino'/
      data flnmovt/'./'/

      integer, public :: lp
c
      common /linepr/ lp
c --- 'lp' = logical unit number for printer output
      data lp/6/


      end module HYCOM_SCALARS
