      block data blkdat
c
      include 'dimensions.h'
      include 'common_blocks.h'
c
c --- layer densities (sigma units):
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- s i g m a _ 0
ccc      data sigma/
ccc     .   24.02, 24.70, 25.28, 25.77, 26.18, 26.52, 26.80, 27.03, 
ccc     .   27.22, 27.38, 27.52, 27.64, 27.74, 27.82, 27.88, 27.92/
c     data thbase/25./
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- s i g m a _ 2
c     data sigma/30.90,31.87,32.75,33.54,34.24,34.85,35.37,35.80,
c    .           36.15,36.43,36.65,36.82,36.95,37.05,37.13,37.20/    !light
ccc   data sigma/31.20,32.51,33.54,34.35,34.99,35.50,35.91,36.24,
ccc  .           36.50,36.70,36.85,36.96,37.04,37.10,37.15,37.20/    !medium
      data sigma/31.85,33.22,34.26,35.04,35.62,36.05,36.37,36.61,
     .           36.79,36.92,37.01,37.07,37.11,37.14,37.17,37.20/    !heavy
css   data sigma/
css  .   29.02,29.89,30.71,31.48,32.20,32.87,33.49,34.06,34.58,35.05,
css  .   35.47,35.84,36.16,36.43,36.65,36.82,36.95,37.05,37.13,37.20/
      data thbase/35./
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c --- 'baclin' = baroclinic time step
c --- 'batrop' = barotropic time step
c --- 'diagfq' = number of days between model diagnostics (incl.output)
      data baclin,batrop/3600.,120./,diagfq/365./          ! 2deg full global
c
c --- 'thkdff' = diffusion velocity (cm/s) for thickness diffusion
c --- 'veldff' = diffusion velocity (cm/s) for momentum dissipation
c --- 'temdff' = diffusion velocity (cm/s) for temp/salin. mixing
c --- 'viscos' is nondimensional, used in deformation-dependent viscosity
c --- 'vertmx' = scale velocity for vertical momentum mixing (m/s)
      data thkdff/.002/,veldff/.1/,temdff/.02/,viscos/0.3/,vertmx/0./
css   data thkdff/.005/,veldff/.1/,temdff/.02/,viscos/0.3/,vertmx/0./
c
c --- 'diapyc' = diapycnal diffusivity times buoyancy freq. (m^2/s^2)
c --- 'mixfrq' = number of time steps between diapycnal mixing calculations
c --- 'h1'     = depth interval used in lateral weighting of hor.pres.grad.
c --- 'thkmin' = minimum mixed-layer thickness (m)
c --- 'acurcy' = permissible roundoff error in column integral calc.
      data diapyc/3.e-7/,mixfrq/12/			!  mixfrq=2/day
css   data diapyc/.3e-4/,mixfrq/12/			!  mixfrq=2/day
      data h1/98060./,thkmin/20./,acurcy/1.e-11/
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
c --- 'cd'     = drag coefficient
c --- 'ct'     = thermal transfer coefficient
c --- 'airdns' = air density at sea level (kg/m^3)
c --- 'evaplh' = latent heat of evaporation (j/kg)
c --- 'thref'  = reference value of specific volume (m^3/kg)
c --- 'epsil'  = small nonzero number used to prevent division by zero
      data g/9.806/,csubp/1005.7/,spcifh/4185./,cd/.0013/,ct/.0012/
      data airdns/1.2/,evaplh/2.47e6/,thref/1.e-3/,epsil/1.e-11/
c
c --- 'itest,jtest' = grid point where detailed diagnostics are desired
      data itest,jtest/176,82/
c
c ---      gridn   --  grid size in degrees longitude,
c ---      xpivn   --  the i index of the equator
c
      data xpivn,gridn/115.,2./
c
c ---  s w i t c h e s    (if set to .true., then...)
c --- thermo      use thermodynamic forcing functions
c --- windf       use wind stress forcing function
c --- relax       activate lateral boundary nudging
c
      data thermo/.true./,windf/.true./,relax/.false./
c
c --- 'lp' = logical unit number for printer output
      data lp/6/
c
c --- use 'huge' to initialize array portions that the code should never access
      data huge/1.e33/
      data nhr/4/                        ! couple every nhr hours
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
      data flnmlat/
     .   '/raid10/nick/HYCOM_KPP/cplinput/latlonij.4bin'/
      data flnmdep/
     .   '/raid10/nick/HYCOM_KPP/cplinput/depth181x180c.4bin'/
      data flnmint/
     .   '/raid10/nick/HYCOM_KPP/cplinput/temp181x180jan_hv_SL.asc'/
      data flnmins/
     .   '/raid10/nick/HYCOM_KPP/cplinput/salt181x180jan_hv_SL.asc'/
      data flnminp/
     .   '/raid10/nick/HYCOM_KPP/cplinput/pout181x180jan_hv_SL.asc'/
      data flnmbas/
     .   '/raid10/nick/HYCOM_KPP/cplinput/basinmasks.181x180b'/
      data flnma2o    /'/raid10/nick/HYCOM_KPP/cplinput/flxa2o.8bin'/
      data flnma2o_tau/'/raid10/nick/HYCOM_KPP/cplinput/taua2o.8bin'/
      data flnmo2a    /'/raid10/nick/HYCOM_KPP/cplinput/ssto2a.8bin'/
      data flnmcoso   /'/raid10/nick/HYCOM_KPP/cplinput/cososino.8bin'/
      data flnmcosa   /'/raid10/nick/HYCOM_KPP/cplinput/cosasina.8bin'/
      data flnmijuv   /'/raid10/nick/HYCOM_KPP/cplinput/ijuvo2a.8bin'/
      data flnmovt/'./'/
c
      end
c
c> Revision history:
c>
c> Dec. 2001 - eliminated rotated grid option ('rotat' switch)
