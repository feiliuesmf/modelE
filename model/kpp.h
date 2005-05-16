      c o m m o n /kppr1/          !  real kpp variables
cc   . wmt(0:nzehat+1,0:nustar+1)  !  momentum velocity scale table
cc   .,wst(0:nzehat+1,0:nustar+1)  !  scalar velocity scale table
     . zgrid(idm,jdm,kdm+1)        !  grid levels in centimeters
     .,vcty(idm,jdm,kdm+1)         !  vert. viscosity coefficient
     .,dift(idm,jdm,kdm+1)         !  vert. heat diffusivity
     .,difs(idm,jdm,kdm+1)         !  vert. salt diffusivity
     .,ghats(idm,jdm,kdm+1)        !  vert. salt diffusivity
     .,vonk          !  von karman constant
     .,zmin,zmax     !  zehat limits for table
     .,umin,umax     !  ustar limits for table
     .,epsilon       ! vertical coordinate scale factor
     .,cmonob        ! constant for calculating monin-obukov length
     .,cekman        ! constant for calculating thickness of ekman layer
     .,rinfty        ! KPP: value for calculating rshear instability
     .,difm0         ! KPP: max viscosity   due to shear instability
     .,difs0         ! KPP: max diffusivity due to shear instability
     .,difmiw        ! KPP: background/internal wave viscosity   (cm^2/s)
     .,difsiw        ! KPP: background/internal wave diffusivity (cm^2/s)
     .,dsfmax        ! KPP: salt fingering diffusivity factor    (cm^2/s)
     .,rrho0         ! KPP: salt fingering rp=(alpha*delT)/(beta*delS)
     .,ricr          ! KPP: critical bulk richardson number
     .,cs            ! KPP: value for nonlocal flux term
     .,cstar         ! KPP: value for nonlocal flux term
     .,cv            ! KPP: value for turb shear contributon to bulk rich. no.
     .,c11           ! KPP: value for turb velocity scale
     .,deltaz        ! delta zehat in table
     .,deltau        ! delta ustar in table
     .,vtc           ! constant for estimating background shear in rib calc.
     .,cg            ! constant for estimating nonlocal flux term of diff. eq.
     .,dp0enh        ! dist. for tapering diff. enhancement at interface nbl-1
     .,dp00,tmljmp,jerlv0
c                 
      common/kppi/
     . niter         ! KPP: iterations for semi-implicit soln. (2 recomended)
c                 
      real zgrid,vcty,vonk,dift,difs,ghats,
     .     zmin,zmax,umin,umax,rinfty,difm0,difs0,difmiw,difsiw,
     .     rrho0,dsfmax,ricr,epsilon,cmonob,cekman,cs,cstar,cv,c11,
     .     deltaz,deltau,vtc,cg,dp0enh,tmljmp,dp00
      integer niter,jerlv0
c                 
      common/kppr2/
     . betard(5)                              !  red extinction coefficient
     .,betabl(5)                              !  blue extinction coefficient
     .,redfac(5)                              !  fract. of penetr. red light
css  .,clouds(idm,jdm,12)                     !  cloud coverage
     .,clouds(idm)                            !  cloud coverage
     .,ac(idm,jdm)                            !  cloud correction coefs
     .,buoyt(idm,jdm)                         ! turbulent surface buoyancy flux
     .,buoylw(idm,jdm)                        ! longwave surface buoyancy flux
     .,buoysw(idm,jdm)                        ! shortwave surface buoyancy flu
      real betard,betabl,redfac,clouds,ac,buoyt,buoylw,buoysw
c                 
      common/hycom5r/
     . hekman(idm,jdm)                        ! ekman layer thickness
     .,dpbl(idm,jdm)                          ! turbulent boundary layer depth
     .,tmix(idm,jdm)                          ! mixed layer temperature
     .,smix(idm,jdm)                          ! mixed layer salinity
     .,thmix(idm,jdm)                         ! mixed layer theta
     .,umix(idm,jdm)                          ! mixed layer u
     .,vmix(idm,jdm)                          ! mixed layer v
c                 
      common/hycom5i/
     . nmlb(idm,jdm,2)                        ! layer containing mlb.
     .,jerlov(idm,jdm)                        ! jerlov water type 1-5
c                 
      real hekman,dpbl,tmix,smix,thmix,umix,vmix
      integer nmlb,jerlov
c                 
      common/kppswtch/mxlkrt,mxlkpp,mxlgiss,shinst,dbdiff,nonloc
      logical mxlkrt,mxlkpp,mxlgiss,shinst,dbdiff,nonloc
c
