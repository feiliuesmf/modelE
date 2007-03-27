c   -----------------------------------------------------------------------------
      include 'bering.h'
c
      c o m m o n
     . u(idm,jdm,2*kdm),v(idm,jdm,2*kdm)      ! velocity components
     .,dp(idm,jdm,2*kdm),dpold(idm,jdm,kdm)   ! layer thickness
     .,dpu(idm,jdm,2*kdm),dpv(idm,jdm,2*kdm)  ! layer thickness at u,v points
     .,p(idm,jdm,kdm+1)                       ! interface pressure
     .,pu(idm,jdm,kdm+1),pv(idm,jdm,kdm+1)    ! interface pres. at u,v points
     .,latij(idm,jdm,4),lonij(idm,jdm,4)      ! latitude/longitude
     .,corio(idm,jdm)                         ! coriolis parameter
     .,potvor(idm,jdm)                        ! potential vorticity
     .,temp(idm,jdm,2*kdm)                    ! temperature
     .,saln(idm,jdm,2*kdm)                    ! salinity
     .,th3d(idm,jdm,2*kdm)                    ! potential density
     .,thstar(idm,jdm,2*kdm)                  ! virtual potential density
     .,thermb(idm,jdm,2*kdm)                  ! difference thstar - th3d
     .,psikk(idm,jdm)                         ! init.montg.pot. in bottom layer
     .,thkk(idm,jdm)                          ! init.thstar in bottom layer
     .,dpmixl(idm,jdm)                        ! Kraus-Turner mixed layer depth
     .,srfhgt(idm,jdm)                        ! sea surface height
c
      real u,v,dp,dpold,dpu,dpv,p,pu,pv,latij,lonij,corio,potvor,
     .     temp,saln,th3d,thstar,thermb,psikk,thkk,dpmixl,srfhgt
c
      c o m m o n
     . montg(idm,jdm,kdm)                     ! montgomery potential
     .,defor1(idm,jdm),defor2(idm,jdm)        ! deformation components
     .,ubavg(idm,jdm,3),vbavg(idm,jdm,3)      ! barotropic velocity
     .,pbavg(idm,jdm,3)                       ! barotropic pressure
     .,ubrhs(idm,jdm),vbrhs(idm,jdm)          ! rhs of barotropic u,v eqns.
     .,utotm(idm,jdm),vtotm(idm,jdm)          ! total (barotrop.+baroclin.)..
     .,utotn(idm,jdm),vtotn(idm,jdm)          ! ..velocities at 2 time levels
     .,uflux(idm,jdm),vflux(idm,jdm)          ! horizontal mass fluxes
     .,uflux1(idm,jdm),vflux1(idm,jdm)        ! more mass fluxes
     .,uflux2(idm,jdm),vflux2(idm,jdm)        ! more mass fluxes
     .,uflux3(idm,jdm),vflux3(idm,jdm)        ! more mass fluxes
     .,uflx(idm,jdm,kdm),vflx(idm,jdm,kdm)    ! more mass fluxes
     .,bolusu(idm,jdm,kdm),bolusv(idm,jdm,kdm)  ! thickness (bolus) fluxes
c
      real montg,defor1,defor2,ubavg,vbavg,pbavg,ubrhs,vbrhs,utotm,
     .     vtotm,utotn,vtotn,uflux,vflux,uflux1,vflux1,uflux2,vflux2,
     .     uflux3,vflux3,uflx,vflx,bolusu,bolusv
c
      c o m m o n /timav/                ! fields needed for time-averaging
     .   uav(idm,jdm,kdm),  vav(idm,jdm,kdm)
     .,dpuav(idm,jdm,kdm),dpvav(idm,jdm,kdm)
     .,temav(idm,jdm,kdm),salav(idm,jdm,kdm)
     .,th3av(idm,jdm,kdm), dpav(idm,jdm,kdm)
     .,ubavav(idm,jdm),vbavav(idm,jdm),pbavav(idm,jdm),sfhtav(idm,jdm)
     .,uflxav(idm,jdm,kdm),vflxav(idm,jdm,kdm)
     .,diaflx(idm,jdm,kdm)                    ! time integral of diapyc.flux
     .,sflxav(idm,jdm),brineav(idm,jdm),eminpav(idm,jdm)
     .,surflav(idm,jdm)
     .,ufxcum(idm,jdm,kdm),vfxcum(idm,jdm,kdm),dpinit(idm,jdm,kdm)
     .,dpmxav(idm,jdm),oiceav(idm,jdm)
c
      real uav,vav,dpuav,dpvav,temav,salav,th3av,dpav,ubavav,vbavav
     .    ,pbavav,sfhtav,uflxav,vflxav,diaflx,sflxav,brineav,eminpav
     .    ,surflav,ufxcum,vfxcum,dpinit
     .    ,dpmxav,oiceav
c
      c o m m o n
     . util1(idm,jdm),util2(idm,jdm)          ! arrays for temporary storage
     .,util3(idm,jdm),util4(idm,jdm)          ! arrays for temporary storage
c
     .,scpx(idm,jdm),scpy(idm,jdm)            ! mesh size at p pts in x,y dir.
     .,scux(idm,jdm),scuy(idm,jdm)            ! mesh size at u pts in x,y dir.
     .,scvx(idm,jdm),scvy(idm,jdm)            ! mesh size at v pts in x,y dir.
     .,scqx(idm,jdm),scqy(idm,jdm)            ! mesh size at q pts in x,y dir.
     .,scu2(idm,jdm),scv2(idm,jdm)            ! grid box size at u,v pts
     .,scp2(idm,jdm),scq2(idm,jdm)            ! grid box size at p,q pts
     .,scuxi(idm,jdm),scvyi(idm,jdm)          ! inverses of scux,scvy
     .,scp2i(idm,jdm),scq2i(idm,jdm)          ! inverses of scp2,scq2
c
     .,pgfx(idm,jdm),pgfy(idm,jdm)            ! horiz. presssure gradient
     .,gradx(idm,jdm),grady(idm,jdm)          ! horiz. presssure gradient
     .,depthu(idm,jdm),depthv(idm,jdm)        ! bottom pres. at u,v points
     .,pvtrop(idm,jdm)                        ! pot.vort. of barotropic flow
     .,depths(idm,jdm)                        ! water depth
     .,drag(idm,jdm)                          ! bottom drag
     .,glue(idm,jdm)                          ! regional viscosity enhancement
     .,dampu(idm,jdm),dampv(idm,jdm)          ! coastal wave damping coeff.
c
      real util1,util2,util3,util4,scpx,scpy,scux,scuy,scvx,scvy,
     .     scqx,scqy,scu2,scv2,scp2,scq2,scuxi,scvyi,scp2i,scq2i,
     .     pgfx,pgfy,gradx,grady,depthu,depthv,pvtrop,depths,drag,
     .     glue,dampu,dampv
c
      c o m m o n
     . uja(idm,jdm),ujb(idm,jdm)              ! velocities at lateral
     .,via(idm,jdm),vib(idm,jdm)              !          neighbor points
     .,pbot(idm,jdm)                          ! bottom pressure at t=0
     .,tracer(idm,jdm,kdm,ntrcr)              ! inert tracer (optional)
     .,tprime(idm,jdm)                        ! temp.change due to surflx
     .,sgain(idm,kdm)                         ! salin.changes from diapyc.mix.
     .,surflx(idm,jdm)                        ! surface thermal energy flux
     .,salflx(idm,jdm)                        ! surface salinity flux
c    .,thkice(idm,jdm)                        ! grid-cell avg. ice thknss (cm)
c    .,covice(idm,jdm)                        ! ice coverage (rel.units)
c    .,temice(idm,jdm)                        ! ice surf.temp.
c    .,odhsi(idm,jdm)                         ! heat borrowed from frozen
     .,odmsi(idm,jdm)                         ! newly formed ice
     .,omlhc(idm,jdm)
     .,dmfz(idm,jdm)                          ! ice mass due to freezing
c
      real uja,ujb,via,vib,pbot,tracer,tprime,sgain,surflx,salflx
c    .   ,thkice,covice,temice,omlhc,dmfz,odhsi
     .   ,odmsi,omlhc,dmfz
c
      integer, dimension (idm,jdm) ::
     .  klist         !k-index of layer below mixl'r
c
      common/int1/klist

c ---  s w i t c h e s    (if set to .true., then...)
c --- diagno      output model fields and diagnostic messages
c --- thermo      use thermodynamic forcing functions
c --- windf       include wind stress in forcing functions
c --- relax       activate lateral boundary nudging
c --- trcout      advect tracer and save results in history/restart file
c --- dotrcr      perform column physics operations on tracer array(s)
c
      logical diagno,thermo,windf,relax,trcout,dotrcr
      common/swtchs/diagno,thermo,windf,relax,trcout,dotrcr
c
      c o m m o n  /frcing/                   !  monthly forcing fields
     . taux(idm,jdm)                          !  wind stress in x direction
     .,tauy(idm,jdm)                          !  wind stress in y direction
c    .,wndspd(idm,jdm,4)                      !  wind speed (tke source)
c    .,airtmp(idm,jdm,4)                      !  pseudo air temperature
c    .,vapmix(idm,jdm,4)                      !  atmosph. vapor mixing ratio
c    .,oprec(idm,jdm)                         !  precipitation
c    .,oevap(idm,jdm)                         !  evaportation
     .,oemnp(idm,jdm)                         !  e - p 
     .,oflxa2o(idm,jdm),oice(idm,jdm)
     .,ustar(idm,jdm)                         ! surface friction velocity
     .,ustarb(idm,jdm)                        ! bottom friction velocity
     .,osalt(idm,jdm)                         ! saltflux from SI(kg/m*m)
c
     .,freshw(idm,jdm)                        !  river & glacier runoff
     .,diafor(idm,jdm)                        !  imposed diapycnal forcing
c
c     real taux,tauy,wndspd,radflx,airtmp,precip,vapmix,freshw,diafor
c    .    ,pwall,swall,twall
      real taux,tauy,oice,oemnp,ustar,ustarb,oflxa2o,osalt
      real freshw,diafor
c
      common/varbls/time,time0,delt1,dlt,w0,w1,w2,w3,ws0,ws1,ws2,ws3
     .     ,area,avgbot,watcum,empcum,slfcum,sala2o,tavini
      real time,time0,delt1,dlt,w0,w1,w2,w3,ws0,ws1,ws2,ws3,
     .     area,avgbot,watcum,empcum,slfcum,sala2o,tavini

      common/varbl2/nstep,nstep0,nstepi,lstep,l0,l1,l2,l3,ls0,ls1
     .             ,ls2,ls3,oddev
      integer       nstep,nstep0,nstepi,lstep,l0,l1,l2,l3,ls0,ls1
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
      common/parms1/thbase,theta(kdm),baclin,batrop,thkdff,
     .              veldff,temdff,viscos,diapyc,vertmx,h1,slip,cbar,
     .              diagfq,wuv1,wuv2,wts1,wts2,acurcy,wbaro,thkmin,
     .              thkbot,botmin,ekman,sigjmp,salmin(kdm)
      real theta,thbase,baclin,batrop,thkdff,veldff,temdff,viscos,
     .     diapyc,vertmx,h1,slip,cbar,diagfq,wuv1,wuv2,wts1,wts2,acurcy,
     .     wbaro,thkmin,thkbot,botmin,ekman,sigjmp,salmin
c
      common/parms2/trcfrq,ntracr,nhr,mixfrq
      integer       trcfrq,ntracr,nhr,mixfrq
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
      common/consts/tenm,onem,tencm,onecm,onemm,g,csubp,spcifh,cd,ct,
     .              airdns,evaplh,thref,epsil,huge,radian,pi
c
      real tenm,onem,tencm,onecm,onemm,g,csubp,spcifh,cd,ct,airdns,
     .     evaplh,thref,epsil,huge,radian,pi
c
c --- grid point where detailed diagnostics are desired:
      common/testpt/itest,jtest
      integer itest,jtest
c
c ---      equatn   --  the i index of the equator
c
      common/grdparms/equatn
      real equatn
c
      character*60 flnmdep,flnmrsi,flnmrso,flnmarc,flnmfor,flnmovt
     .            ,flnmini,flnmriv,flnmbas,flnmdia,flnmlat
     .            ,flnminp,flnmint,flnmins
     .            ,flnmcoso,flnmcosa,flnma2o,flnma2o_tau,flnmo2a
     .            ,flnmo2a_e,flnmo2a_n
      common/iovars/flnmdep,flnmrsi,flnmrso,flnmarc,flnmfor,flnmovt
     .            ,flnmini,flnmriv,flnmbas,flnmdia,flnmlat
     .            ,flnminp,flnmint,flnmins
     .            ,flnmcoso,flnmcosa,flnma2o,flnma2o_tau,flnmo2a
     .            ,flnmo2a_e,flnmo2a_n
c
c> Revision history:
c>
c> July 1997 - eliminated 3-D arrays -uold,vold- (used in time smoothing)
c-----------------------------------------------------------------------------
