       block data blkkpp
c
#include "dimensions.h"
#include "kpp.h"
c
c --- 'vonk'      = von karman constant
c --- 'zmin,zmax' = zehat limits for velocity scale lookup table
c --- 'umin,umax' = ustar limits for velocity scale lookup table
c --- 'rinfty'    = default = 0.7
c --- 'difm0'     = max viscosity due to shear instability
c --- 'difs0'     = max diffusivity due to shear instability
c --- 'difmiw'    = background/internal wave viscosity (cm^2/s)
c --- 'difsiw'    = background/internal wave diffusivity (cm^2/s)
c --- 'rrho0'     = rp=(alpha*delT)/(beta*delS)
c --- 'dsfmax'    = constant for estimating salt fingering diffusivity
c --- 'ricr'      = critical bulk richardson number
c --- 'epsilon'   = vertical coordinate scale factor
c --- 'cmonob'    = constant for calculating monin-obukov length
c --- 'cekman'    = constant for calculating thickness of ekman layer
c --- 'cs'        = constant for calculating nonlocal flux term
c --- 'cv'        = constant for computing turb shear contributon to bulk ri
c --- 'c11'       = constant for calculating turb velocity scale
c --- 'cstar'     = constant for calculating nonlocal flux term
c --- 'niter'     = no. of iterations in semi-implicit soln. (2 recomended)
c
      data
css  . vonk/.4/,zmin/-.4e-6/,zmax/0./,umin/0./,umax/4.e-2/,rinfty/.7/,
css  . difm0/50.e-4/,difs0/50.e-4/,difmiw/1.e-4/,difsiw/.1e-5/,
     . vonk/.4/,zmin/-.4e-6/,zmax/0./,umin/0./,umax/.16/,rinfty/.7/,
     . difm0/50.e-4/,difs0/50.e-4/,difmiw/1.e-4/,difsiw/1.e-5/,
     . rrho0/1.9/,
     . dsfmax/10.e-4/,ricr/.3/,epsilon/.1/,cmonob/1./,cekman/.7/,
     . cs/98.96/,cv/1.8/,c11/5./,cstar/10./,niter/2/
c
c --- 'jerlv0' = initial jerlov water type (1 to 5)
      data jerlv0/2/
c --- red and blue light extinction coefficients (pressure units)
c --- for jerlov water types 1 to 5 - fraction of penetrating red light
c     betard(1) =  0.35*onem
c     betard(2) =  0.6 *onem
c     betard(3) =  1.0 *onem
c     betard(4) =  1.5 *onem
c     betard(5) =  1.4 *onem
c     betabl(1) = 23.0 *onem
c     betabl(2) = 20.0 *onem
c     betabl(3) = 17.0 *onem
c     betabl(4) = 14.0 *onem
c     betabl(5) =  7.9 *onem
c     redfac(1) = 0.58
c     redfac(2) = 0.62
c     redfac(3) = 0.67
c     redfac(4) = 0.77
c     redfac(5) = 0.78
      data betard/3432.10,5883.60,806.0,14709.0,13728.40/
      data betabl/225538.0,196120.0,166702.0,137284.0,77467.4/
      data redfac/0.58,0.62,0.67,0.77,0.78/
c --- 'tmljmp' = equivalent temperature jump across mixed-layer (degC)
      data tmljmp/0.2/
c
c --- 'mxlkrt' = KT:  HYCOM Kraus-Turner provided by Rainer Bleck
c --- 'mxlgiss' = GISS: activate GISS mixing model
c --- 'mxlkpp' = KPP: activate mixed layer model (mlflag==1)
c --- 'shinst' = KPP: activate shear instability mixing
c --- 'dbdiff' = KPP: activate double diffusion mixing
c --- 'nonloc' = KPP: activate nonlocal b. layer mixing
      data mxlkrt/.false./,mxlkpp/.true./,mxlgiss/.false./
     .    ,shinst/.true./,dbdiff/.true./
     .    ,nonloc/.false./
      end
