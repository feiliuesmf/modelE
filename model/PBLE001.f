C**** PBLE001 E001M12 SOMTQ PBL_f90 PBLB336EM12
C**** changes for f90
C****
C**** PBLB336EM12 BB285SM12 PBLB295M9     yc 02/11/99
C****
C**** newer SOC pbl
C**** PBLB237JM12                        yc 01/29/99
C**** OPT(3) GOSTMT
C**** All routines necessary for the SOC BL
C**** All lines that are resolution dependent are indicated here
C****  in the update file.
C****
C****
C****
C --------------------------------------------------------------------
C     ABL model using second order closure model -- Greg Hartke
C
C        ocean and ocean ice are treated as rough surfaces
C        roughness lengths from Brutsaert for rough surfaces
C
C --------------------------------------------------------------------
      SUBROUTINE PBL(I,J,ITYPE)
C --------------------------------------------------------------------
C     Variable definitions:
C
C The variables passed thru the common block PBLPAR are parameters
C  necessary to do the PBL solution. These variables have been passed
C  from subroutine SURFCE or subroutine EARTH. The variables are:
C
C     ZGS    = height of the surface layer which is 10 m everywhere.
C     ZS1    = height of the first model layer (m)
C     PIJ    = surface pressure at gridpoint (i,j) (mb)
C     PSK    = surface pressure to the power KAPA
C     TGV    = virtual potential temperature of the ground (K)
C     TKV    = virtual potential temperature of first model layer (K)
C     THV1   = virtual temperature of the first model layer (K)
C     HEMI   = 1 for northern hemisphere, -1 for southern hemisphere
C     SHA    = specific heat at constant pressure (RGAS/KAPA)
C     OMEGA2 = 2.*OMEGA where OMEGA is the angular frequency of
C              the earth (sec)
C     NLAT2  = JM/2, half the total number of latitude grid points
C     JVPO   = 2 at south pole, JM at north pole, otherwise not used
C     IQ1    = IM/4+1
C     IQ2    = IM/2+1
C     IQ3    = 3*IM/4+1
C     IM1    = I except at the poles, where it equals IM
C     POLE   = .TRUE. if at the north or south pole, .FALSE. otherwise
C
C The quantities passed thru common block PBLOUT constitute the output
C  from this PBL subroutine. The variables are:
C
C     US     = x component of surface wind, postive eastward (m/sec)
C     VS     = y component of surface wind, positive northward (m/sec)
C     WS     = magnitude of the surface wind (m/sec)
C     TSV    = virtual potential temperature of the surface (K)
C     QS     = surface value of the specific moisture
C     PSI    = difference in direction between geostrophic and surface
C              winds (radians)
C     DBL    = boundary layer height (m)
C     KM     = momentum transport coefficient charavterizing the
C              boundary layer (m**2/sec)
C     KH     = heat transport coefficient evaluated at ZGS (m**2/sec)
C     USTAR  = friction speed (square root of momentum flux) (m/sec)
C     PPBL   = pressure at DBL (mb)
C     CM     = drag coefficient (dimensionless surface momentum flux)
C     CH     = Stanton number   (dimensionless surface heat flux)
C     CQ     = Dalton number    (dimensionless surface moisture flux)
C     UG     = x component of the geostrophic wind, positive eastward
C              (m/sec)
C     VG     = y component of the geostrophic wind, positive northward
C              (m/sec)
C     WG     = magnitude of the geostrophic wind (m/sec)
C
C --------------------------------------------------------------------
      USE E001M12_COM
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/RDATA/ROUGHL(IM,JM)

      REAL*8 KM,KH,KMSURF,KHSURF
      LOGICAL POLE

      COMMON /PBLPAR/ZGS,ZS1,PIJ,PSK,TGV,TKV,THV1,QG,HEMI,SHA,
     2               OMEGA2,DTSURF,JVPO,IQ1,IQ2,IQ3,IM1,POLE

      COMMON /PBLOUT/US,VS,WS,TSV,QS,PSI,DBL,KM,KH,USTAR,PPBL,
     2               CM,CH,CQ,UG,VG,WG,ZMIX

      PARAMETER (NLAT2=JM/2,EPSLON=1.D-20,radian=3.141592654/180.)

C**** BLDATA 1  COMPOSITE SURFACE WIND MAGNITUDE (M/S)
C****        2  COMPOSITE SURFACE AIR TEMPERATURE (K)
C****        3  COMPOSITE SURFACE AIR SPECIFIC HUMIDITY (1)
C****        4  LAYER TO WHICH DRY CONVECTION MIXES (1)
C****        5  MIXED LAYER DEPTH (Z1O NOT YET PART OF RESTART FILE)
C****        6  COMPOSITE SURFACE U WIND
C****        7  COMPOSITE SURFACE V WIND
C****        8  COMPOSITE SURFACE MOMENTUM TRANSFER (TAU)
C new     9-12  ustar for each ITYPE (sqrt of srfc momentum flux) (m/s)
C****
C**** ROUGHL    LOG10(ZGS/ROUGHNESS LENGTH), prescribed with ZGS=30 m.
C****
      IF (ITYPE.GT.2) THEN
        Z0M=30./(10.**ROUGHL(I,J))
      ENDIF
      ztop=zgs+zs1
      IF (TKV.EQ.TGV) TGV=1.0001*TGV
      IF (TKV.GE.TGV) THEN
C **********************************************************************
C ********** ATMOSPHERE IS STABLE WITH RESPECT TO THE GROUND ***********
C
C  DETERMINE THE VERTICAL LEVEL CORRESPONDING TO THE HEIGHT OF THE PBL:
C   WHEN ATMOSPHERE IS STABLE, CAN COMPUTE DBL BUT DO NOT KNOW THE
C   INDEX OF THE LAYER.
C
        USTAR=BLDATA(I,J,8+ITYPE)
        DBL=0.3*USTAR/OMEGA2
        if (dbl.le.ztop) then
C THE VERTICAL LEVEL FOR WHICH WG IS COMPUTED IS THE FIRST:
          dbl=ztop
          L=1
          ELSE
          if (dbl.gt.3000.) dbl=3000.
C FIND THE VERTICAL LEVEL NEXT HIGHER THAN DBL AND COMPUTE WG THERE:
          zpbl=ztop
          pl1=pij*sig(1)+ptop
          tl1=t(i,j,1)*expbyk(pl1)
          do 100 l=2,ls1
            pl=pij*sig(l)+ptop
            tl=t(i,j,l)*expbyk(pl)
            tbar=thbar(tl1,tl)
            zpbl=zpbl+(rgas/grav)*tbar*(pl1-pl)/pl1
            if (zpbl.ge.dbl) go to 200
            pl1=pl
            tl1=tl
100       continue
200       CONTINUE
        ENDIF
C *********************************************************************
        ELSE
C *********************************************************************
C ********* ATMOSPHERE IS UNSTABLE WITH RESPECT TO THE GROUND *********
C
C  LDC IS THE LEVEL TO WHICH DRY CONVECTION MIXES. IF THE BOUNDARY
C   LAYER HEIGHT IS LESS THAN 3 KM, ASSIGN LDC TO L, OTHERWISE MUST
C   FIND INDEX FOR NEXT MODEL LAYER ABOVE 3 KM:
C
        LDC=BLDATA(I,J,4)
        IF (LDC.EQ.0) LDC=1
        if (ldc.eq.1) then
          dbl=ztop
          l=1
          else
          zpbl=ztop
          pl1=pij*sig(1)+ptop
          tl1=t(i,j,1)*expbyk(pl1)
          zpbl1=ztop
          do 300 l=2,ldc
            pl=pij*sig(l)+ptop
            tl=t(i,j,l)*expbyk(pl)
            tbar=thbar(tl1,tl)
            zpbl=zpbl+(rgas/grav)*tbar*(pl1-pl)/pl1
            if (zpbl.ge.3000.) then
              zpbl=zpbl1
              go to 400
            endif
            pl1=pl
            tl1=tl
            zpbl1=zpbl
300       continue
400       continue
          l=l-1
          dbl=zpbl
        endif
C**********************************************************************
      ENDIF

C *********************************************************************
      ppbl=sige(l)*pij+ptop
      if (l.gt.ls1) ppbl=sige(l)*(psf-ptop)+ptop
      phi=radian*(float(j-1)*180./float(jm-1)-90.)
      coriol=sin(phi)*omega2
      ttop=tkv
      qtop=q(i,j,1)
      tgrnd=tgv
      qgrnd=qg

      if (pole) then
        utop=.25*(u(1,jvpo,1)-u(iq2,jvpo,1)
     2          -(v(iq1,jvpo,1)-v(iq3,jvpo,1))*hemi)
        vtop=.25*(v(1,jvpo,1)-v(iq2,jvpo,1)
     2          +(u(iq1,jvpo,1)-u(iq3,jvpo,1))*hemi)
        ug  =.25*(u(1,jvpo,l)-u(iq2,jvpo,l)
     2          -(v(iq1,jvpo,l)-v(iq3,jvpo,l))*hemi)
        vg  =.25*(v(1,jvpo,l)-v(iq2,jvpo,l)
     2          +(u(iq1,jvpo,l)-u(iq3,jvpo,l))*hemi)
        else
        utop=.25*(u(im1,j,1)+u(i,j,1)+u(im1,j+1,1)+u(i,j+1,1))
        vtop=.25*(v(im1,j,1)+v(i,j,1)+v(im1,j+1,1)+v(i,j+1,1))
        ug  =.25*(u(im1,j,l)+u(i,j,l)+u(im1,j+1,l)+u(i,j+1,l))
        vg  =.25*(v(im1,j,l)+v(i,j,l)+v(im1,j+1,l)+v(i,j+1,l))
      endif

      call advanc(us,vs,tsv,qs,kmsurf,khsurf,ustar,ug,vg,cm,ch,cq,
     2            z0m,z0h,z0q,coriol,utop,vtop,ttop,qtop,tgrnd,
     3            qgrnd,zgs,ztop,zmix,dtsurf,ufluxs,vfluxs,
     4            tfluxs,qfluxs,i,j,itype)

      ws    =sqrt(us*us+vs*vs)
      wg    =sqrt(ug*ug+vg*vg)
      km    =kmsurf
      kh    =khsurf
      psitop=atan2(vg,ug+epslon)
      psisrf=atan2(vs,us+epslon)
      psi   =psisrf-psitop
      bldata(i,j,8+itype)=ustar
C ******************************************************************

      RETURN
      END

      subroutine advanc(us,vs,tsv,qs,kmsurf,khsurf,ustar,ug,vg,cm,ch,cq,
     2                  z0m,z0h,z0q,coriol,utop,vtop,ttop,qtop,tgrnd,
     3                  qgrnd,zgs,ztop,zmix,dtime,ufluxs,vfluxs,
     4                  tfluxs,qfluxs,ilong,jlat,itype)
c ----------------------------------------------------------------------
c this routine time advances the solutions for the wind components u and
c  v, the temperature t, the moisture (or other passive scalar) q, and
c  the turbulence energy e. the quantities u, v, t, and q are computed
c  on the primary grid and e is computed on the secondary grid which is
c  staggered with respect to the primary grid. other quantities computed
c  on the secondary grid are km, kh, ke, lscale (the turbulence length
c  scale) and the dimensionless gradients gm and gh.
c ----------------------------------------------------------------------
c
c internal quantities:
c
c     n           = number of vertical grid points
c     u(n)        = x component of wind
c     v(n)        = y component of wind
c     t(n)        = temperature
c     q(n)        = specific humidity (a passive scalar)
c     e(n-1)      = turbulence energy. computed on the secondary grid.
c     lscale(n-1) = turbulence length scale. computed on secondary grid.
c     z(n)        = altitude of primary vertical grid points
c     zhat(n)     = altitude of secondary vertical grid points
c     dzh(n-1)    = dz evaluated at zhat(i)
c     dz(n)       = dz evaluated at z(i)
c     dxi         = (z(n)-z(1))/(n-1)
c     km(n-1)     = turbulent momentum tranport coefficient. computed
c                   on the secondary grid.
c     kh(n-1)     = turbulent thermometric conductivity. computed
c                   on the secondary grid.
c     ke(n-1)     = tranport coefficient for the turbulence energy.
c                   computed on the secondary grid.
c     ipbl(im,jm,4) array used to keep track of grid points and surface
c                   types for which bl properties were computed at the
c                   last time step. (keeps track of formation and
c                   melting of ice.)
c
c ----------------------------------------------------------------------
c quantities in the call to this routine:
c   output:
c     us          = x component of the surface wind (i.e., due east)
c     vs          = y component of the surface wind (i.e., due north)
c     tsv         = virtual potential surface temperature
c     qs          = surface specific moisture
c     kmsurf      = surface value of km
c     khsurf      = surface value of kh
c     ustar       = friction speed
c     coriol      = 2.*omega*sin(latitude), the coriolis factor
c     zmix        = magic quantity needed in surfce
c     from subroutine dflux:{
c     cm          = dimensionless momentum flux at surface (drag
c                   coefficient)
c     ch          = dimensionless heat flux at surface (stanton number)
c     cq          = dimensionless moisture flux at surface (dalton
c                   number)
c     z0m         = roughness height for momentum (if itype=1 or 2)
c     z0h         = roughness height for heat
c     z0q         = roughness height for moisture }
c   input:
c     z0m         = roughness height for momentum (if itype>2)
c     ug          = x component of the geostrophic wind
c     vg          = y component of the geostrophic wind
c     utop        = x component of wind at the top of the layer
c     vtop        = y component of wind at the top of the layer
c     ttop        = temperature at the top of the layer
c     qtop        = moisture at the top of the layer
c     tgrnd       = temperature of the ground, i.e., at the roughness
c                   height for temperature
c     qgrnd       = moisture at the ground level, i.e., at the roughness
c                   height for temperature
c     zgs         = height of the surface layer (nominally 10 m)
c     ztop        = height of the first model layer, approx 200 m in
c                   the 9 layer model
c     dtime       = time step
c     ilong       = index of the longitude coordinate
c     jlat        = index of the latitude coordinate (these were used
c                   for diagnostics in the development phase.)
c     itype       = 1, ocean
c                 = 2, ocean ice
c                 = 3, land ice
c                 = 4, land
c ----------------------------------------------------------------------
c   value of bgrid (=0.293 based on zs1=200.) is chosen to minimize
c    differences in surface temperature, surface wind direction,
c    surface momentum flux, surface heat flux, and surface moisture flux
c    compared to full simulation for z=[10.,3000.] and n=128.
c ----------------------------------------------------------------------
c ----------------------------------------------------------------------
      implicit real*8 (a-h,k,o-z)
      real*8 lmonin
      parameter (n=8) !n = npbl in the calling routine pbl        .
      parameter (tol=1.e-4)
      PARAMETER (IM=72,JM=46)
      parameter (iprint= 0,jprint=33)  ! set iprint=52 e.g.  to debug
      real*8 lscale(n-1)
      dimension z(n),zhat(n-1),xi(n),xihat(n-1)
      dimension dz(n),dzh(n-1)
      dimension u(n),v(n),t(n),q(n),e(n-1)
      dimension km(n-1),kh(n-1),ke(n-1),gm(n-1),gh(n-1)
      dimension usave(n),vsave(n),tsave(n),qsave(n),esave(n)
      data ifirst/1/
      save bgrid
      common /socabl/uabl(n,im,jm,4),vabl(n,im,jm,4),
     2               tabl(n,im,jm,4),qabl(n,im,jm,4),eabl(n,im,jm,4),
     3               cmgs(im,jm,4),chgs(im,jm,4),cqgs(im,jm,4),
     4               ipbl(im,jm,4)
      common /pblvar/u,v,t,q,e
c
      itmax=1
      if (ifirst.eq.1) then
        call getb(zgs,zdummy,bgrid)
        ifirst=0
      endif

c     if(abs(bgrid-0.29269d0).gt.0.01) then
c     write(99,*) bgrid,ilong,jlat,itype
c     endif

      call griddr(z,zhat,xi,xihat,dz,dzh,zgs,ztop,bgrid,n)
      zmix=dzh(1)+zgs

      do 100 i=1,n-1
        u(i)=uabl(i,ilong,jlat,itype)
        v(i)=vabl(i,ilong,jlat,itype)
        t(i)=tabl(i,ilong,jlat,itype)
        q(i)=qabl(i,ilong,jlat,itype)
        e(i)=eabl(i,ilong,jlat,itype)
        usave(i)=u(i)
        vsave(i)=v(i)
        tsave(i)=t(i)
        qsave(i)=q(i)
        esave(i)=e(i)
100   continue
      u(n)=uabl(n,ilong,jlat,itype)
      v(n)=vabl(n,ilong,jlat,itype)
      t(n)=tabl(n,ilong,jlat,itype)
      q(n)=qabl(n,ilong,jlat,itype)
      usave(n)=u(n)
      vsave(n)=v(n)
      tsave(n)=t(n)
      qsave(n)=q(n)

      cm=cmgs(ilong,jlat,itype)
      ch=chgs(ilong,jlat,itype)
      cq=cqgs(ilong,jlat,itype)
c ----------------------------------------------------------------------
c First, make trial advances of the solutions of the prognostic fields.
c   Must first compute subsidiary quantities such as length scale,
c   transport coefficients, and field scales.
c   NB: ustar0 need be defined only if itmax > 1.
c   Step 1:

      call getl1(e,zhat,dzh,lscale,n)
      call getk(km,kh,ke,gm,gh,u,v,t,e,lscale,dzh,n)
      call stars(ustar,tstar,qstar,lmonin,tgrnd,qgrnd,
     2           u,v,t,q,z,z0m,z0h,z0q,cm,ch,cq,
     3           km,kh,dzh,itype,n)

      call getl2(e,t,zhat,dzh,lscale,ustar,lmonin,n)
      call getk(km,kh,ke,gm,gh,u,v,t,e,lscale,dzh,n)
      call stars(ustar,tstar,qstar,lmonin,tgrnd,qgrnd,
     2           u,v,t,q,z,z0m,z0h,z0q,cm,ch,cq,
     3           km,kh,dzh,itype,n)
c     ustar0=ustar

      call eeqns(esave,e,u,v,t,km,kh,ke,lscale,dz,dzh,
     2               ustar,dtime,n)
      call tqeqns(u,v,tsave,qsave,t,q,z,kh,dz,dzh,
     2            ch,cq,tstar,qstar,z0h,z0q,tgrnd,qgrnd,
     3            ttop,qtop,dtime,n)

      call uveqns(usave,vsave,u,v,z,km,dz,dzh,
     2            ustar,cm,z0m,utop,vtop,dtime,coriol,
     3            ug,vg,ilong,jlat,n)

c     if ((itype.eq.4).or.(itype.eq.3)) then
        if ((ttop.gt.tgrnd).and.(lmonin.lt.0.)) then
          call tfix(t,z,ttop,tgrnd,lmonin,n)
          itmax=2
        endif
c     endif

c     call elevl2(e,u,v,t,km,kh,lscale,dzh,ustar,n)
      call getl2(e,t,zhat,dzh,lscale,ustar,lmonin,n)

c ----------------------------------------------------------------------
c Now iteratively recompute the fields. If itmax > 1, restore commented
c   code that performs test to see if solution has converged. This
c   condition obtains if ustar remains stable to a defined limit
c   between iterations.
c   NB: test need be defined only if itmax > 1. Also, final computation
c       of lscale (via call to getl2 before line 200) should be included
c       if itmax > 1.
c       eeqn is called if level 2.5 is used, elevl2 for level 2 soln.
c   Step 2:

      do 200 iter=1,itmax

        call getk(km,kh,ke,gm,gh,u,v,t,e,lscale,dzh,n)
        call stars(ustar,tstar,qstar,lmonin,tgrnd,qgrnd,
     2             u,v,t,q,z,z0m,z0h,z0q,cm,ch,cq,
     3             km,kh,dzh,itype,n)

c       test=abs(ustar-ustar0)/(ustar+ustar0)
c       if (test.lt.tol) go to 300
c       ustar0=ustar

        call eeqns(esave,e,u,v,t,km,kh,ke,lscale,dz,dzh,
     2                 ustar,dtime,n)
        call tqeqns(u,v,tsave,qsave,t,q,z,kh,dz,dzh,
     2              ch,cq,tstar,qstar,z0h,z0q,tgrnd,qgrnd,
     3              ttop,qtop,dtime,n)

        call uveqns(usave,vsave,u,v,z,km,dz,dzh,
     2              ustar,cm,z0m,utop,vtop,dtime,coriol,
     3              ug,vg,ilong,jlat,n)

        if ((iter.eq.1).and.(itmax.eq.2)) then
c         call elevl2(e,u,v,t,km,kh,lscale,dzh,ustar,n)
          call getl2(e,t,zhat,dzh,lscale,ustar,lmonin,n)
        endif

200   continue
300   continue
c ----------------------------------------------------------------------

      do 400 i=1,n-1
        uabl(i,ilong,jlat,itype)=u(i)
        vabl(i,ilong,jlat,itype)=v(i)
        tabl(i,ilong,jlat,itype)=t(i)
        qabl(i,ilong,jlat,itype)=q(i)
        eabl(i,ilong,jlat,itype)=e(i)
400   continue
      uabl(n,ilong,jlat,itype)=u(n)
      vabl(n,ilong,jlat,itype)=v(n)
      tabl(n,ilong,jlat,itype)=t(n)
      qabl(n,ilong,jlat,itype)=q(n)

      us    = u(1)
      vs    = v(1)
      tsv   = t(1)
      qs    = q(1)
      kmsurf= km(1)
      khsurf= kh(1)
      cmgs(ilong,jlat,itype)=cm
      chgs(ilong,jlat,itype)=ch
      cqgs(ilong,jlat,itype)=cq
      ipbl(ilong,jlat,itype)=1

      ufluxs=km(1)*(u(2)-u(1))/dzh(1)
      vfluxs=km(1)*(v(2)-v(1))/dzh(1)
      tfluxs=kh(1)*(t(2)-t(1))/dzh(1)
      qfluxs=kh(1)*(q(2)-q(1))/dzh(1)

c ----------------------------------------------------------------------
c Diagnostics printed at a selected point:

      if ((ilong.eq.iprint).and.(jlat.eq.jprint)) then
        call output(u,v,t,q,e,lscale,z,zhat,dzh,
     2              km,kh,gm,gh,cm,ch,cq,z0m,z0h,z0q,
     3              ustar,tstar,qstar,lmonin,tgrnd,qgrnd,
     4              utop,vtop,ttop,qtop,
     5              dtime,bgrid,ilong,jlat,iter,itype,n)
      endif
c
c     call check1(ustar,1,ilong,jlat,2)
c     call check1(u,n,ilong,jlat,1)
c     call check1(v,n,ilong,jlat,2)
c     call check1(t,n,ilong,jlat,3)
c     call check1(q,n,ilong,jlat,4)
c ----------------------------------------------------------------------
      return
      end

      subroutine stars(ustar,tstar,qstar,lmonin,tgrnd,qgrnd,
     2                 u,v,t,q,z,z0m,z0h,z0q,cm,ch,cq,
     3                 km,kh,dzh,itype,n)
c ----------------------------------------------------------------------
c Computes the friction speed USTAR, temperature scale TSTAR, and
c  moisture scale QSTAR. The fluxes of momentum, heat, and moisture
c  should be constant in the surface layer, with
c   Momentum flux = USTAR*USTAR
c   Heat flux     = USTAR*TSTAR
c   MOISTURE flux = USTAR*QSTAR
c LMONIN is the Monin-Obukhov length scale.
c The call to dflux computes the drag coefficient, Stanton number, and
c  Dalton number.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,k,o-z)
      parameter (kappa=0.40,grav=9.81)
      parameter (smax=0.25,smin=0.005,cmax=smax*smax,cmin=smin*smin)
      dimension u(n),v(n),t(n),q(n),z(n)
      dimension km(n-1),kh(n-1),dzh(n-1)
      real*8 lmonin

c ----------------------------------------------------------------------
      dz     = dzh(1)
      vel1   = sqrt(u(1)*u(1)+v(1)*v(1))
c     vel2   = sqrt(u(2)*u(2)+v(2)*v(2))
c     dudz   = (vel2-vel1)/dz
      du1=u(2)-u(1)
      dv1=v(2)-v(1)
      dudz=sqrt(du1*du1+dv1*dv1)/dz
      dtdz   = (t(2)-t(1))/dz
      dqdz   = (q(2)-q(1))/dz
      ustar  = sqrt(km(1)*dudz)
      ustar  = max(ustar,1.d-20)
      tstar  = kh(1)*dtdz/ustar
      qstar  = kh(1)*dqdz/ustar
      zgs    = z(1)
      if (ustar.gt.smax*vel1) ustar=smax*vel1
      if (abs(tstar).gt.smax*abs(t(1)-tgrnd)) tstar=smax*(t(1)-tgrnd)
      if (abs(qstar).gt.smax*abs(q(1)-qgrnd)) qstar=smax*(q(1)-qgrnd)
      if (ustar.lt.smin*vel1) ustar=smin*vel1
      if (abs(tstar).lt.smin*abs(t(1)-tgrnd)) tstar=smin*(t(1)-tgrnd)
      if (abs(qstar).lt.smin*abs(q(1)-qgrnd)) qstar=smin*(q(1)-qgrnd)

      lmonin = ustar*ustar*tgrnd/(kappa*grav*tstar)
      call dflux(lmonin,ustar,vel1,z0m,z0h,z0q,zgs,cm,ch,cq,itype)
c ----------------------------------------------------------------------

      return
      end

      subroutine getl1(e,zhat,dzh,lscale,n)
c ----------------------------------------------------------------------
c Computes the master length scale of the turbulence model. Keep in mind
c  that LSCALE is computed on the secondary grid. This routine computes
c  l0 (the asymptotic length scale) according to the usual prescription
c  used in the Mellor and Yamada models.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,k,o-z)
      parameter (alpha=0.20,kappa=0.40)
      real*8 lscale(n-1),l0,l1
      dimension e(n-1),zhat(n-1),dzh(n-1)

      sum1=0.
      sum2=0.
      do 100 j=1,n-1
        sum1=sum1+sqrt(e(j))*zhat(j)*dzh(j)
        sum2=sum2+sqrt(e(j))*dzh(j)
100   continue
      l0=alpha*sum1/sum2

      do 200 i=1,n-1
        l1=kappa*zhat(i)
        lscale(i)=l0*l1/(l0+l1)
c       lscale(i)=l1
200   continue

      return
      end

      subroutine getl2(e,t,zhat,dzh,lscale,ustar,lmonin,n)
c ----------------------------------------------------------------------
c Computes the master length scale of the turbulence model. Note
c  that LSCALE is computed on the secondary grid. l0 in this routine is
c  computed via an analytic approximation to the l0 computed in the
c  full domain simulation.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,k,l,o-z)
      parameter (kappa=0.40,grav=9.81)
      parameter (lcoef=0.060,omega=7.292e-5)
      dimension e(n-1),t(n),zhat(n-1),dzh(n-1),lscale(n-1)

      l0=lcoef*sqrt(ustar*abs(lmonin)/omega)
      if (l0.lt.zhat(1)) l0=zhat(1)

      l1=kappa*zhat(1)
      lscale(1)=l0*l1/(l0+l1)

      do 300 i=2,n-1
        l1=kappa*zhat(i)
        lscale(i)=l0*l1/(l0+l1)
        if (t(i+1).gt.t(i)) then
          bvfrq2=grav*log(t(i+1)/t(i))/dzh(i)
          lmax  =0.75*sqrt(e(i)/(bvfrq2+1.e-40))
          if (lscale(i).gt.lmax) lscale(i)=lmax
        endif
        if (lscale(i).lt.0.5*kappa*zhat(i)) lscale(i)=0.5*kappa*zhat(i)
        if (lscale(i).gt.dzh(i)) lscale(i)=dzh(i)
300   continue

      return
      end

      subroutine dflux(lmonin,ustar,vsurf,z0m,z0h,z0q,zgs,
     2                 cm,ch,cq,itype)
c *********************************************************************
c   Computes the dimensionless surface fluxes of momentum, heat,
c     and moisture (drag coefficient, Stanton number, and Dalton
c     number) for the surface layer.
c
c  inputs are:
c    lmonin = Monin-Obukhov length (m)
c    ustar  = friction speed (sqrt of surface momentum flux) (m/sec)
c    vsurf  = total surface wind speed, used to limit ustar -> cm
c    z0m    = momentum roughness length, prescribed if itype=3 or 4 (m)
c    zgs    = height of the surface layer (m)
c    itype  = integer identifying surface type
c
c  outputs are:
c    cm    = drag coefficient for momentum
c    ch    = Stanton number
c    cq    = Dalton number
c    z0m   = roughness length for momentum, computed if itype=1 or 2
c    z0h   = roughness length for temperature (m)
c    z0q   = roughness length for water vapor (m)
c *********************************************************************
      implicit real*8 (a-h,k,o-z)
      real*8 lmonin,nu,num,nuh,nuq
      parameter (nu=1.5e-5,num=0.135*nu,nuh=0.395*nu,nuq=0.624*nu)
      parameter (epslon=1.e-20)
      parameter (grav=9.81,kappa=0.40)
      parameter (sigma=0.95,sigma1=1.-sigma)
      parameter (gamamu=19.3,gamahu=11.6,gamams=4.8,gamahs=8./sigma)
      parameter (smax=0.25,smin=0.005,cmax=smax*smax,cmin=smin*smin)

      if ((itype.eq.1).or.(itype.eq.2)) then
c *********************************************************************
c  Compute roughness lengths using rough surface formulation:
        z0m=num/ustar+0.018*ustar*ustar/grav
        if (z0m.gt.0.2) z0m=0.2
        if (ustar.le.0.02) then
          z0h=nuh/ustar
          if (z0h.gt.0.5852) z0h=0.5852
          z0q=nuq/ustar
          if (z0q.gt.0.92444) z0q=0.92444
          else
          r0q=(ustar*z0m/nu)**0.25
          z0h=7.4*z0m*exp(-2.4604*r0q)
          z0q=7.4*z0m*exp(-2.2524*r0q)
          if (ustar.lt.0.2) then
            beta=(ustar-0.02)/0.18
            z0h=(1.-beta)*nuh/ustar+beta*z0h
            z0q=(1.-beta)*nuq/ustar+beta*z0q
          endif
        endif
c *********************************************************************
        else
c *********************************************************************
c  For land and land ice, z0m is specified. For z0h and z0q,
c    empirical evidence suggests:
        z0h=z0m*exp(-2.)
        z0q=z0h
c *********************************************************************
      endif

      cmn=kappa*kappa/(log(zgs/z0m)**2)
      chn=kappa*kappa/(log(zgs/z0m)*log(zgs/z0h))
      cqn=kappa*kappa/(log(zgs/z0m)*log(zgs/z0q))

      zgsbyl=zgs/lmonin
      z0mbyl=z0m/lmonin
      z0hbyl=z0h/lmonin
      z0qbyl=z0q/lmonin

c *********************************************************************
c  Now compute DPSI, which is the difference in the psi functions
c    computed at zgs and the relevant roughness height:
c *********************************************************************

      if (lmonin.gt.0.) then
c *********************************************************************
c  Here the atmosphere is stable with respect to the ground:
        dpsim=-gamams*(zgsbyl-z0mbyl)
        dpsih= sigma1*log(zgs/z0h)-sigma*gamahs*(zgsbyl-z0hbyl)
        dpsiq= sigma1*log(zgs/z0q)-sigma*gamahs*(zgsbyl-z0qbyl)
c *********************************************************************
        else
c *********************************************************************
c  Here the atmosphere is unstable with respect to the ground:
        xms  =    (1.-gamamu*zgsbyl)**0.25
        xm0  =    (1.-gamamu*z0mbyl)**0.25
        xhs  =sqrt(1.-gamahu*zgsbyl)
        xh0  =sqrt(1.-gamahu*z0hbyl)
        xqs  =sqrt(1.-gamahu*zgsbyl)
        xq0  =sqrt(1.-gamahu*z0qbyl)
        dpsim=log((1.+xms)*(1.+xms)*(1.+xms*xms)/
     2           ((1.+xm0)*(1.+xm0)*(1.+xm0*xm0)))-
     3        2.*(atan(xms)-atan(xm0))
        dpsih=sigma1*log(zgs/z0h)+2.*sigma*log((1.+xhs)/(1.+xh0))
        dpsiq=sigma1*log(zgs/z0q)+2.*sigma*log((1.+xqs)/(1.+xq0))
c *********************************************************************
      endif

      dm=1./(1.-sqrt(cmn)*dpsim/kappa)**2
      dh=sqrt(dm)/(1.-chn*dpsih/(kappa*sqrt(cmn)))
      dq=sqrt(dm)/(1.-cqn*dpsiq/(kappa*sqrt(cmn)))
      if (dm.lt.1.e-4) dm=1.e-4
      if (dh.lt.1.e-4) dh=1.e-4
      if (dq.lt.1.e-4) dq=1.e-4

      cm=dm*cmn
      ch=dh*chn
      cq=dq*cqn
      if (cm.gt.cmax) cm=cmax
      if (ch.gt.cmax) ch=cmax
      if (cq.gt.cmax) cq=cmax

      if (cm.lt.cmin) cm=cmin
      if (ch.lt.cmin) ch=cmin
      if (cq.lt.cmin) cq=cmin

      return
      end

      subroutine simil(u,t,q,z,ustar,tstar,qstar,
     2                 z0m,z0h,z0q,lmonin,tg,qg)
      implicit real*8 (a-h,k,l,o-z)
      parameter (kappa=0.40)
      parameter (sigma=0.95,sigma1=1.-sigma)
      parameter (gamamu=19.3,gamahu=11.6,gamams=4.8,gamahs=8./sigma)
c *********************************************************************
c  Computes the similarity solutions for wind speed, temperature, and
c   moisture mixing ratio at height z. (Used for diagnostic output.)
c
c  Inputs:
c     z      = height above ground at which soln is being computed (m)
c     ustar  = friction speed (m/sec)
c     tstar  = temperature scale (K)
c     qstar  = moisture scale
c     z0m    = momentum roughness height (m)
c     z0h    = temperature roughness height (m)
c     z0q    = moisture roughness height (m)
c     lmonin = Monin-Obukhov length scale (m)
c     tg     = ground temperature (K)
c     qg     = ground moisture mixing ratio
c
c  Outputs:
c     u      = computed similarity solution for wind speed (m/sec)
c     t      = computed similarity solution for temperature (K)
c     q      = computed similarity solution for moisture mixing ratio
c *********************************************************************

      zbyl  =z  /lmonin
      z0mbyl=z0m/lmonin
      z0hbyl=z0h/lmonin
      z0qbyl=z0q/lmonin

c *********************************************************************
c  Now compute DPSI, which is the difference in the psi functions
c    computed at zgs and the relevant roughness height:
c *********************************************************************

      if (lmonin.gt.0.) then
c *********************************************************************
c  Here the atmosphere is stable with respect to the ground:
        dpsim=-gamams*(zbyl-z0mbyl)
        dpsih= sigma1*log(z/z0h)-sigma*gamahs*(zbyl-z0hbyl)
        dpsiq= sigma1*log(z/z0q)-sigma*gamahs*(zbyl-z0qbyl)
c *********************************************************************
        else
c *********************************************************************
c  Here the atmosphere is unstable with respect to the ground:
        xm   =    (1.-gamamu*  zbyl)**0.25
        xm0  =    (1.-gamamu*z0mbyl)**0.25
        xh   =sqrt(1.-gamahu*  zbyl)
        xh0  =sqrt(1.-gamahu*z0hbyl)
        xq   =sqrt(1.-gamahu*  zbyl)
        xq0  =sqrt(1.-gamahu*z0qbyl)
        dpsim=log((1.+xm )*(1.+xm )*(1.+xm *xm )/
     2           ((1.+xm0)*(1.+xm0)*(1.+xm0*xm0)))-
     3        2.*(atan(xm )-atan(xm0))
        dpsih=sigma1*log(z/z0h)+2.*sigma*log((1.+xh)/(1.+xh0))
        dpsiq=sigma1*log(z/z0q)+2.*sigma*log((1.+xq)/(1.+xq0))
c *********************************************************************
      endif

      u=   (ustar/kappa)*(log(z/z0m)-dpsim)
      t=tg+(tstar/kappa)*(log(z/z0h)-dpsih)
      q=qg+(qstar/kappa)*(log(z/z0q)-dpsiq)

      return
      end

      subroutine getb(zgs,ztop,bgrid)
c ----------------------------------------------------------------------
c This routine computes the value of bgrid to be used in the gridding
c  scheme. This parameter determines the strength of the logarithmic
c  component in the log-linear scheme. This fitting for bgrid was
c  determined by doing a series of off-line runs comparing the reduced
c  domain simulation to the full BL simulation and determining the value
c  of bgrid that gave the best fit for a range of ztop = [50.,200.] m
c  for the reduced domain simulation.
c This form for z1 = zgs + zs1 (in terms of GCM parameters) yields an
c  average value for zs1. The quantity theta was computed on the
c  assumption of zs1=200 m from the original 9-layer model (actually
c  was misconstrued as z1 = 200 m when it should have been zs1 = 200 m)
c  and is then applied to all vertical resolutions.
c
c Input:
c
c    zgs   = The height of the surface layer.
c
c Output:
c
c    ztop  = The height of the top of the BL simulation domain.
c            Corresponds to the height of the middle of the first model
c            layer and is only needed if the BL fields require
c            initialization.
c    bgrid = The parameter that determines the strength of the log
c            term in the log-linear gridding scheme.
c ----------------------------------------------------------------------
      USE E001M12_COM
      IMPLICIT REAL*8 (A-H,O-Z)
      theta=269.0727251
      z1=zgs+0.5*(1.-sige(2))*(psf-ptop)*rgas*theta/(grav*psf)
      x=z1/100.
      ztop=z1
      bgrid=0.177427*x**4 - 1.0504*x**3 + 2.34169*x**2 -
     2      2.4772*x + 1.44509
      return
      end

      subroutine griddr(z,zhat,xi,xihat,dz,dzh,z1,zn,bgrid,n)
c ----------------------------------------------------------------------
c Computes altitudes on the vertical grid. The XI coordinates are
c  uniformly spaced and are mapped in a log-linear fashion onto the Z
c  grid. (The Z's are the physical coords.) Also computes the altitudes
c  on the secondary grid, ZHAT(I), and the derivatives dxi/dz evaluated
c  at both all Z(I) and ZHAT(I). The differentials dz and dzh used
c  throughout the code are
c
c     dz  = deltaxi/(dxi/dz)
c     dzh = deltaxi/(dxi/dzh)
c
c  where deltaxi = dxi = (ztop - zbottom)/(n-1)
c
c The parameter BGRID determines how strongly non-linear the mapping is.
c  BGRID=0 gives linear mapping. Increasing BGRID packs more points into
c  the bottom of the layer.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (tolz=1.e-3)
      dimension z(n),zhat(n-1)
      dimension xi(n),xihat(n-1)
      dimension dz(n),dzh(n-1)
      common /grids_99/z1pass,znpass,b,xipass
      external fgrid

      z1pass=z1
      znpass=zn
      dxi=(zn-z1)/float(n-1)
      b=bgrid
      zmin=z1
      zmax=zn

      do 100 i=1,n-1
        xi(i)=z1+(zn-z1)*float(i-1)/float(n-1)
        xihat(i)=z1+(zn-z1)*(float(i)-0.5)/float(n-1)
100   continue
      xi(n)=zn
      z(1)=z1

      dxidz=1.+bgrid*((zn-z1)/z1-log(zn/z1))
      dz(1)=dxi/dxidz
      xipass=xihat(1)
      zhat(1)=zbrent(fgrid,zmin,zmax,tolz)
      dxidzh=1.+bgrid*((zn-z1)/zhat(1)-log(zn/z1))
      dzh(1)=dxi/dxidzh

      do 200 i=2,n-1
        xipass=xi(i)
        z(i)=zbrent(fgrid,zmin,zmax,tolz)
        xipass=xihat(i)
        zhat(i)=zbrent(fgrid,zmin,zmax,tolz)
        dxidz=1.+bgrid*((zn-z1)/z(i)-log(zn/z1))
        dxidzh=1.+bgrid*((zn-z1)/zhat(i)-log(zn/z1))
        dz(i)=dxi/dxidz
        dzh(i)=dxi/dxidzh
200   continue
      z(n)=zn
      dxidz=1.+bgrid*((zn-z1)/zn-log(zn/z1))
      dz(n)=dxi/dxidz

      return
      end

      function fgrid(z)
c ----------------------------------------------------------------------
c Subsidiary function used for computing the grid points. This defines
c   the functional relationship between z and xi.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /grids_99/z1,zn,bgrid,xi
      fgrid=z+bgrid*((zn-z1)*log(z/z1)-(z-z1)*log(zn/z1))-xi
      return
      end

      function zbrent(func,x1,x2,tol)
c ----------------------------------------------------------------------
c Uses Brent's method to solve FUNC=0. X1 and X2 must bracket the root.
c Taken from Numerical Recipes.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (itmax=100,eps=1.e-3)
      a=x1
      b=x2
      fa=func(a)
      fb=func(b)
      if(fb*fa.gt.0.) then
        write (99,*)  'root must be bracketed for zbrent.'
        stop
      endif
      fc=fb
      do 11 iter=1,itmax
        if(fb*fc.gt.0.) then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*eps*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b)
11    continue
      write (99,*) 'zbrent exceeding maximum iterations.'
      write (99,*) x1,b,x2
      stop 'zbrent'
      return
      end

      subroutine loadbl
c ----------------------------------------------------------------------
c             This routine checks to see if ice has
c              melted or frozen out of a grid box.
c
c For ITYPE=1 (ocean; melted ocean ice since last time step):
c  If there was no computation made for ocean at the last time step,
c  this time step may start from ocean ice result. If there was no
c  ocean nor ocean ice computation at the last time step, nothing
c  need be done.
c
c For ITYPE=2 (ocean ice; frozen from ocean since last time step):
c  If there was no computation made for ocean ice at the last time step,
c  this time step may start from ocean result. If there was no
c  ocean nor ocean ice computation at the last time step, nothing
c  need be done.
c
c For ITYPE=3 (land ice; frozen on land since last time step):
c  If there was no computation made for land ice at the last time step,
c  this time step may start from land result. If there was no
c  land ice nor land computation at the last time step, nothing
c  need be done.
c
c For ITYPE=4 (land; melted land ice since last time step):
c  If there was no computation made for land at the last time step,
c  this time step may start from land ice result. If there was no
c  land nor land ice computation at the last time step, nothing
c  need be done.
c
c In the current version of the GCM, there is no need to check the
c  land or land ice components of the grid box for ice formation and
c  melting because pland and plice are fixed. The source code to do
c  this is retained and deleted in the update deck in the event this
c  capability is added in future versions of the model.
c ----------------------------------------------------------------------
      USE E001M12_COM
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (npbl=8)
      common /socabl/uabl(npbl,im,jm,4),vabl(npbl,im,jm,4),
     2               tabl(npbl,im,jm,4),qabl(npbl,im,jm,4),
     3               eabl(npbl,im,jm,4),
     4               cmgs(im,jm,4),chgs(im,jm,4),cqgs(im,jm,4),
     5               ipbl(im,jm,4)

      do 600 j=1,jm
        if ((j.eq.1).or.(j.eq.jm)) then
          imax=1
          else
          imax=im
        endif
        do 500 i=1,imax

c ******* itype=1: Ocean

          if (ipbl(i,j,1).eq.0) then
            if (ipbl(i,j,2).eq.1) then
              do 100 lpbl=1,npbl-1
                uabl(lpbl,i,j,1)=uabl(lpbl,i,j,2)
                vabl(lpbl,i,j,1)=vabl(lpbl,i,j,2)
                tabl(lpbl,i,j,1)=tabl(lpbl,i,j,2)
                qabl(lpbl,i,j,1)=qabl(lpbl,i,j,2)
                eabl(lpbl,i,j,1)=eabl(lpbl,i,j,2)
  100         continue
              uabl(npbl,i,j,1)=uabl(npbl,i,j,2)
              vabl(npbl,i,j,1)=vabl(npbl,i,j,2)
              tabl(npbl,i,j,1)=tabl(npbl,i,j,2)
              qabl(npbl,i,j,1)=qabl(npbl,i,j,2)
              cmgs(i,j,1)=cmgs(i,j,2)
              chgs(i,j,1)=chgs(i,j,2)
              cqgs(i,j,1)=cqgs(i,j,2)
              bldata(i,j,9)=bldata(i,j,10)
            endif
          endif

c ******* itype=2: Ocean ice

          if (ipbl(i,j,2).eq.0) then
            if (ipbl(i,j,1).eq.1) then
              do 200 lpbl=1,npbl-1
                uabl(lpbl,i,j,2)=uabl(lpbl,i,j,1)
                vabl(lpbl,i,j,2)=vabl(lpbl,i,j,1)
                tabl(lpbl,i,j,2)=tabl(lpbl,i,j,1)
                qabl(lpbl,i,j,2)=qabl(lpbl,i,j,1)
                eabl(lpbl,i,j,2)=eabl(lpbl,i,j,1)
  200         continue
              uabl(npbl,i,j,2)=uabl(npbl,i,j,1)
              vabl(npbl,i,j,2)=vabl(npbl,i,j,1)
              tabl(npbl,i,j,2)=tabl(npbl,i,j,1)
              qabl(npbl,i,j,2)=qabl(npbl,i,j,1)
              cmgs(i,j,2)=cmgs(i,j,1)
              chgs(i,j,2)=chgs(i,j,1)
              cqgs(i,j,2)=cqgs(i,j,1)
              bldata(i,j,10)=bldata(i,j,9)
            endif
          endif

c ******* itype=3: Land ice

          if (ipbl(i,j,3).eq.0) then
            if (ipbl(i,j,4).eq.1) then
              do 300 lpbl=1,npbl-1
                uabl(lpbl,i,j,3)=uabl(lpbl,i,j,4)
                vabl(lpbl,i,j,3)=vabl(lpbl,i,j,4)
                tabl(lpbl,i,j,3)=tabl(lpbl,i,j,4)
                qabl(lpbl,i,j,3)=qabl(lpbl,i,j,4)
                eabl(lpbl,i,j,3)=eabl(lpbl,i,j,4)
  300         continue
              uabl(npbl,i,j,3)=uabl(npbl,i,j,4)
              vabl(npbl,i,j,3)=vabl(npbl,i,j,4)
              tabl(npbl,i,j,3)=tabl(npbl,i,j,4)
              qabl(npbl,i,j,3)=qabl(npbl,i,j,4)
              cmgs(i,j,3)=cmgs(i,j,4)
              chgs(i,j,3)=chgs(i,j,4)
              cqgs(i,j,3)=cqgs(i,j,4)
              bldata(i,j,11)=bldata(i,j,12)
            endif
          endif

c ******* itype=4: Land

          if (ipbl(i,j,4).eq.0) then
            if (ipbl(i,j,3).eq.1) then
              do 400 lpbl=1,npbl-1
                uabl(lpbl,i,j,4)=uabl(lpbl,i,j,3)
                vabl(lpbl,i,j,4)=vabl(lpbl,i,j,3)
                tabl(lpbl,i,j,4)=tabl(lpbl,i,j,3)
                qabl(lpbl,i,j,4)=qabl(lpbl,i,j,3)
                eabl(lpbl,i,j,4)=eabl(lpbl,i,j,3)
  400         continue
              uabl(npbl,i,j,4)=uabl(npbl,i,j,3)
              vabl(npbl,i,j,4)=vabl(npbl,i,j,3)
              tabl(npbl,i,j,4)=tabl(npbl,i,j,3)
              qabl(npbl,i,j,4)=qabl(npbl,i,j,3)
              cmgs(i,j,4)=cmgs(i,j,3)
              chgs(i,j,4)=chgs(i,j,3)
              cqgs(i,j,4)=cqgs(i,j,3)
              bldata(i,j,12)=bldata(i,j,11)
            endif
          endif

  500   continue
  600 continue

      return
      end

      subroutine pgrads1
      USE E001M12_COM
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (iq1=im/4+1,iq2=im/2+1,iq3=3*im/4+1)
      common /pgpass/dpdxr(im,jm),dpdyr(im,jm),phi(im,jm)
     &              ,dpdxr0(im,jm),dpdyr0(im,jm)

c     for gcm main level 1:
      call geopot

      do 200 j=2,jm-1
        do 100 i=1,im

          pij=p(i,j)
          p1=100.*(pij*sig(1)+ptop)
          t1=t(i,j,1)*expbyk(pij*sig(1)+ptop)
          rho1=p1/(rgas*t1)

          index1=i+1
          index2=i-1
          if (i.eq.1) then
            index2=im
          endif
          if (i.eq.im) then
            index1=1
          endif
          dpx=100.*(p(index1,j)-p(index2,j))*sig(1)

          dpdxr(i,j)=0.5*dpx/(dxp(j)*rho1)+0.5*
     2               (phi(index1,j)-phi(index2,j))/dxp(j)

          dpy=100.*(p(i,j+1)-p(i,j-1))*sig(1)

          dpdyr(i,j)=0.5*dpy/(dyp(j)*rho1)+0.5*
     2               (phi(i,j+1)-phi(i,j-1))/dyp(j)

c         at the surface:
          dpx=100.*(p(index1,j)-p(index2,j))

          dpdxr0(i,j)=0.5*dpx/(dxp(j)*rho1)+0.5*
     2               (FDATA(index1,j,1)-FDATA(index2,j,1))/dxp(j)

          dpy=100.*(p(i,j+1)-p(i,j-1))

          dpdyr0(i,j)=0.5*dpy/(dyp(j)*rho1)+0.5*
     2               (FDATA(i,j+1,1)-FDATA(i,j-1,1))/dyp(j)

100     continue
200   continue

      i=1
      j=1
      pij=p(i,j)
      p1=100.*(pij*sig(1)+ptop)
      t1=t(i,j,1)*expbyk(pij*sig(1)+ptop)
      rho1 =p1/(rgas*t1)
      dypsp=2.*dyp(1)

      j=jm
      pij=p(i,j)
      p1=100.*(pij*sig(1)+ptop)
      t1=t(i,j,1)*expbyk(pij*sig(1)+ptop)
      rhojm=p1/(rgas*t1)
      dypnp=2.*dyp(jm)

      dpdxr(1, 1)=0.25*(p(iq1  ,   2)-p(iq3  ,   2)+
     2                  p(iq1+1,   2)-p(iq3+1,   2))*sig(1)*100./
     3                 (dypsp*rho1)+
     4            0.125*(phi(iq1  ,   2)-phi(iq3  ,   2)+
     5                  phi(iq1+1,   2)-phi(iq3+1,   2))/dypsp
      dpdyr(1, 1)=0.25*(p(    1,   2)-p(iq2  ,   2)+
     2                  p(    2,   2)-p(iq2+1,   2))*sig(1)*100./
     3                 (dypsp*rho1)+
     4            0.125*(phi(    1,   2)-phi(iq2  ,   2)+
     5                  phi(    2,   2)-phi(iq2+1,   2))/dypsp

      dpdxr(1,jm)=0.25*(p(iq1  ,jm-1)-p(iq3  ,jm-1)+
     2                  p(iq1+1,jm-1)-p(iq3+1,jm-1))*sig(1)*100./
     3                 (dypnp*rhojm)+
     4            0.125*(phi(iq1  ,jm-1)-phi(iq3  ,jm-1)+
     5                  phi(iq1+1,jm-1)-phi(iq3+1,jm-1))/dypnp
      dpdyr(1,jm)=0.25*(p(iq2  ,jm-1)-p(    1,jm-1)+
     2                  p(iq2+1,jm-1)-p(    2,jm-1))*sig(1)*100./
     3                 (dypnp*rhojm)+
     4            0.125*(phi(iq2  ,jm-1)-phi(    1,jm-1)+
     5                  phi(iq2+1,jm-1)-phi(    2,jm-1))/dypnp

c     at the surface:

      dpdxr0(1, 1)=0.25*(p(iq1  ,   2)-p(iq3  ,   2)+
     2                  p(iq1+1,   2)-p(iq3+1,   2))*100./
     3                 (dypsp*rho1)+
     4            0.125*(FDATA(iq1  ,   2,1)-FDATA(iq3  ,   2,1)+
     5                  FDATA(iq1+1,   2,1)-FDATA(iq3+1,   2,1))/dypsp
      dpdyr0(1, 1)=0.25*(p(    1,   2)-p(iq2  ,   2)+
     2                  p(    2,   2)-p(iq2+1,   2))*100./
     3                 (dypsp*rho1)+
     4            0.125*(FDATA(    1,   2,1)-FDATA(iq2  ,   2,1)+
     5                  FDATA(    2,   2,1)-FDATA(iq2+1,   2,1))/dypsp



      dpdxr0(1,jm)=0.25*(p(iq1  ,jm-1)-p(iq3  ,jm-1)+
     2                  p(iq1+1,jm-1)-p(iq3+1,jm-1))*100./
     3                 (dypnp*rhojm)+
     4            0.125*(FDATA(iq1  ,jm-1,1)-FDATA(iq3  ,jm-1,1)+
     5                  FDATA(iq1+1,jm-1,1)-FDATA(iq3+1,jm-1,1))/dypnp
      dpdyr0(1,jm)=0.25*(p(iq2  ,jm-1)-p(    1,jm-1)+
     2                  p(iq2+1,jm-1)-p(    2,jm-1))*100./
     3                 (dypnp*rhojm)+
     4            0.125*(FDATA(iq2  ,jm-1,1)-FDATA(    1,jm-1,1)+
     5                  FDATA(iq2+1,jm-1,1)-FDATA(    2,jm-1,1))/dypnp

      return
      end

      subroutine geopot
      USE E001M12_COM
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (zgs=10.)
      common /pgpass/dpdxr(im,jm),dpdyr(im,jm),phi(im,jm)
c     note: FDATA(I,J,1) is the geopotential height (9.81*zatm)
c
c     for GCM main level 1:
      do 200 j=2,jm-1
        do 100 i=1,im
          pij=p(i,j)
          p1=pij*sig(1)+ptop
          p1k=expbyk(p1)
          t1=t(i,j,1)*p1k
          z1=zgs+0.5*dsig(1)*rgas*t1*pij/(p1*grav)
          phi(i,j)=grav*z1+FDATA(I,J,1)
100     continue
200   continue

      i=1
      j=1
      pij=p(i,j)
      p1=pij*sig(1)+ptop
      p1k=expbyk(p1)
      t1=t(i,j,1)*p1k
      z1=zgs+0.5*dsig(1)*rgas*t1*pij/(p1*grav)
      phi(i,j)=grav*z1+FDATA(I,J,1)

      j=jm
      pij=p(i,j)
      p1=pij*sig(1)+ptop
      p1k=expbyk(p1)
      t1=t(i,j,1)*p1k
      z1=zgs+0.5*dsig(1)*rgas*t1*pij/(p1*grav)
      phi(i,j)=grav*z1+FDATA(I,J,1)
      return
      end

      subroutine tfix(t,z,ttop,tgrnd,lmonin,n)
      implicit real*8 (a-h,k,o-z)
      real*8 lmonin
      dimension t(n),z(n)
      tsurf=tgrnd+0.2*(ttop-tgrnd)
      t(1)=tsurf
      do 100 i=2,n-1
        t(i)=tsurf+(z(i)-z(1))*(ttop-tsurf)/(z(n)-z(1))
 100  continue
      lmonin=abs(lmonin)
      return
      end

      subroutine ccoeff0
      implicit real*8 (a-h,o-z)
      common/const0/rimax,ghmin,ghmax,gmmax0
      common/const1/d1,d2,d3,d4,d5,s0,s1,s2,s3,s4,s5,s6
      common/const2/c1,c2,c3,c4,c5   ! for level 2 calculation
      common/const3/b1,b123
c----------------------------------------------------------------------
      b1=16.6
      b123=b1**(2./3.)
      g1= 1.87
      g2= 1.63
      g3= 3.14
      g4= 1.33
      g5= 2.31
      g6= 0.4
      g7= 0.0
      g8= 6.73
c
      D1 = -G5*(7./3.*G4+G8)
      D2 = G3**2-1./3.*G2**2-1./4.*G5**2*(G6**2-G7**2)
      D3 = 1./3.*G4*G5**2*(4.*G4+3.*G8)
      D4 = -G5*( 1./3.*G4*(G3**2-G2**2)+G8*(G3**2-1./3.*G2**2)
     &      -G4*G5*(G3*G7-1./3.*G2*G6) )
      D5 = -1./12.*G5**2*(3.*G3**2-G2**2)*(G6**2-G7**2)
      S0 = 1./2.*G1
      S1 = 1./18.*G5*( 3.*G4*G5*(G6+G7)+2.*G4*(G2+3*G3)
     &      -3.*G1*(4.*G4+3.*G8) )
      S2 = -1./8.*G1*G5**2*(G6**2-G7**2)
      S3 = 0.
      S4 = 1./3.*G5
      S5 = -1./3.*G4*G5**2
      S6 = -1./36.*G5*( 6.*G1*(3.*G3-G2)-9.*G1*G5*(G6-G7)
     &      -4.*(3.*G3**2-G2**2) )
c
c      x = -B1**2*Gh
c      y = B1**2*Gm
c    (B1/2*sm_num/den - Sm_num/Den);
c
c      Den  = 1. + D1*Gh + D2*Gm + D3*Gh**2 + D4*Gh*Gm + D5*Gm**2
c      Sm_num = S0 + S1*Gh + S2*Gm
c      Sh_num = S4 + S5*Gh + S6*Gm
c
c     find rimax:
      c1=s5*b1-d3
      c2=d4-(s1+s6)*b1
      c3=s2*b1-d5
      rimax=(-c2-sqrt(c2*c2-4.*c1*c3))/(2*c1)
      rimax=int(rimax*1000.)/1000.
      c4 = d1-s4*b1
      c5 = s0*b1-d2
      write(99,*) "in ccoeff0, rimax=",rimax
c     find ghmax:
      del=(2.*d1+s4)**2-8.*(2.*d3+s5)
      ghmax=( b1*s4-d1-sqrt( (b1*s4-d1)**2-4*(d3-b1*s5) ) )
     &    / (2*(d3-b1*s5))
      ghmax=int(ghmax*10000.)/10000.
      ghmin=-0.281d0
      gmmax0=3.80d0
c
      write(99,*) "g1=",g1
      write(99,*) "g2=",g2
      write(99,*) "g3=",g3
      write(99,*) "g4=",g4
      write(99,*) "g5=",g5
      write(99,*) "g6=",g6
      write(99,*) "g7=",g7
      write(99,*) "g8=",g8
      write(99,*) "rimax=",rimax
      write(99,*) "ghmax=",ghmax
      write(99,*) "ghmin=",ghmin
      write(99,*) "gmmax0=",gmmax0
      return
      end


      subroutine getk(km,kh,ke,gma,gha,u,v,t,e,lscale,dzh,n)
      implicit real*8 (a-h,k,o-z)
      parameter (grav=9.81,sq=0.2)
      dimension km(n-1),kh(n-1),ke(n-1),gma(n-1),gha(n-1)
      dimension u(n),v(n),t(n),e(n-1),dzh(n-1)
      real*8 lscale(n-1)
      common/const0/rimax,ghmin,ghmax,gmmax0
      common/const1/d1,d2,d3,d4,d5,s0,s1,s2,s3,s4,s5,s6
      common/const2/c1,c2,c3,c4,c5   ! for level 2 calculation
      common/const3/b1,b123
      data ifirst/1/
c-----------------------------------------------------------------------
      if(ifirst.eq.1) then  !this loop is called only once
          call ccoeff0
          ifirst=0
      endif
      do 100 i=1,n-1
        an2=2.*grav*(t(i+1)-t(i))/((t(i+1)+t(i))*dzh(i))
        dudz=(u(i+1)-u(i))/dzh(i)
        dvdz=(v(i+1)-v(i))/dzh(i)
        as2=dudz*dudz+dvdz*dvdz
        ell=lscale(i)
        qturb=sqrt(2.*e(i))
        tau=ell/max(qturb,1.d-20)
        gh=-tau*tau*an2
        gm=tau*tau*as2
        if(gh.lt.ghmin) gh=ghmin
        if(gh.gt.ghmax) gh=ghmax
        gmmax=(1.+d1*gh+d3*gh*gh)/(d2+d4*gh)
        gmmax=min(gmmax,gmmax0)
        if(gm.gt.gmmax) gm=gmmax
        Den=1.+D1*Gh+D2*Gm+D3*Gh*Gh+D4*Gh*Gm+D5*Gm*Gm
        sm=(S0+S1*Gh+S2*Gm)/Den
        sh=(S4+S5*Gh+S6*Gm)/Den
        km(i)=min(max(ell*qturb*sm,1.5d-5),100.d0)
        kh(i)=min(max(ell*qturb*sh,2.5d-5),100.d0)
        ke(i)=min(max(ell*qturb*sq,1.5d-5),100.d0)
        gma(i)=gm
        gha(i)=gh
100   continue
c
1004  format(16(1pe16.5))
      return
      end


      subroutine trislv(a,b,c,r,u,n)
c-----to solve the real tridiagonal difference matrix equation
      implicit real*8 (a-h,o-z)
      dimension a(n),b(n),c(n),r(n),u(n),gam(40)
      if(b(1).eq.0.) then
          write(99,*) 'b(1).eq.0.,stop'
          stop
      endif
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
          gam(j)=c(j-1)/bet
          bet=b(j)-a(j)*gam(j)
          if(bet.eq.0.) then
              write(99,*) 'bet.eq.0., stop'
              write(99,*) j,b(j),a(j),gam(j)
              write(99,*) c(j-1),b(1),r(1),u(1)
              stop
          endif
          u(j)=(r(j)-a(j)*u(j-1))/bet
 11   continue
      do 12 j=n-1,1,-1
          u(j)=u(j)-gam(j+1)*u(j+1)
 12   continue
      return
      end

      subroutine eeqns(esave,e,u,v,t,km,kh,ke,lscale,
     &                     dz,dzh,ustar,dtime,n)
      implicit real*8 (a-h,k,o-z)
      parameter (grav=9.81,nm=8)
      dimension esave(n),e(n-1),u(n),v(n),t(n)
      dimension km(n-1),kh(n-1),ke(n-1)
      real*8 lscale(n-1)
      dimension dz(n),dzh(n-1)
      common /trid_b/sub(nm),dia(nm),sup(nm),rhs(nm),rhs1(nm)
      common/const0/rimax,ghmin,ghmax,gmmax0
      common/const3/b1,b123
c ----------------------------------------------------------------------
c     sub(j)*e_jm1_kp1+dia(j)*e_j_kp1+sup(j)*e_jp1_kp1 = rhs(j)
c     from kirk:/u/acyxc/papers/2ndOrder/maple/phik.1,
c       sub(j)=-dtime*aj*P1a_jm1_k/dxi^2
c       sup(j)=-dtime*aj*P1a_jp1_k/dxi^2
c       dia(j)=1-(sub(j)+dia(j))+dtime*P3_j_k
c       rhs(j)=PHI_j_k + dtime*P4_j_k
c       where P1a == P1*a
c       aj=dxi/dz(j) if j refers to the primary grid
c       aj=dxi/dzh(j) if j refers to the secondary grid
c       now for e-eqn j refers to the secondary grid, therefore:
c       aj = dxi/dzh(j)
c       a(j-1/2) = dxi/dzh(j-1/2)=dxi/dz(j)
c       a(j+1/2) = dxi/dzh(j+1/2)=dxi/dz(j+1)
c       p1(j-1/2)=0.5*(ke(j)+ke(j-1))
c       p1(j+1/2)=0.5*(ke(j)+ke(j+1))
c       ke(j)=sq*lscale(j)*qturb, qturb=sqrt(2.*e(j))
c ----------------------------------------------------------------------
c
      do 200 j=2,n-2
          qturb=sqrt(2.*e(j))
          sub(j)=-dtime*0.5*(ke(j)+ke(j-1))/(dzh(j)*dz(j))
          sup(j)=-dtime*0.5*(ke(j)+ke(j+1))/(dzh(j)*dz(j+1))
          dia(j)=1.-(sub(j)+sup(j))+dtime*2*qturb/(b1*lscale(j))
          an2=2*grav*(t(j+1)-t(j))/((t(j+1)+t(j))*dzh(j))
          dudz=(u(j+1)-u(j))/dzh(j)
          dvdz=(v(j+1)-v(j))/dzh(j)
          as2=dudz*dudz+dvdz*dvdz
          rhs(j)=esave(j)+dtime*(km(j)*as2-kh(j)*an2)
200   continue
c
      dia(1)=1.
      sup(1)=0.
      rhs(1)=0.5*b123*ustar*ustar
c
      sub(n-1)=-1.
      dia(n-1)=1.
      rhs(n-1)=0.
c
      call trislv(sub,dia,sup,rhs,e,n-1)
c
      do 400 j=1,n-1
         if(e(j).lt.1.d-20) e(j)=1.d-20
400   continue
c
      return
      end

      subroutine tqeqns(u,v,t0,q0,t,q,z,kh,dz,dzh,
     2                  ch,cq,tstar,qstar,z0h,z0q,tgrnd,qgrnd,
     3                  ttop,qtop,dtime,n)
c ----------------------------------------------------------------------
c this routine computes the matrices for the solutions of the t and q
c  fields as well as appling the boundary conditions.
c  the boundary conditions at the bottom are:
c
c     kh * dt/dz = ch * usurf * (t - tg)
c     kh * dq/dz = cq * usurf * (q - qg), evaluated at the lowest level.
c
c  at the top, the temperature and moisture are prescribed.
c
c the arrays at and aq are the lhs's of the appropriate discretized
c  equation. the vectors bt and bq are similarly the rhs's.
c note that t(i) and q(i) that appear in this routine are the
c  temperature and moisture mixing ratio from the previous time step.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,k,o-z)
      parameter (kappa=0.40,grav=9.81,nm=8)
      dimension u(n),v(n),t0(n),q0(n),t(n),q(n),z(n),kh(n-1)
      dimension dz(n),dzh(n-1)
      common /trid_b/sub(nm),dia(nm),sup(nm),rhs(nm),rhs1(nm)

c ----------------------------------------------------------------------
c compute the lhss of the prognostic equations. the (1,j) and (n,j)
c  components are zeroed to accomodate boundary conditions:

      do 200 i=2,n-1
          sub(i)=-dtime/(dz(i)*dzh(i-1))*kh(i-1)
          sup(i)=-dtime/(dz(i)*dzh(i))*kh(i)
          dia(i)=1.-(sub(i)+sup(i))
200   continue
c ----------------------------------------------------------------------
c now compute the rhss. components (1) and (n) are assigned as bcs:

      do 500 i=2,n-1
        rhs(i)=t0(i)
        rhs1(i)=q0(i)
500   continue
c ----------------------------------------------------------------------
c finally, apply the boundary conditions that the flux be continuous
c  across the boundary at the surface (level 1 of the sub-grid scale
c  model, see comments at beginning of this routine):
c     M.J.Miller et al. 1992:
      wstar3=-1000.*grav*kh(1)*( 2.*(t(2)-t(1))/(t(2)+t(1))
     &                        -(q(2)-q(1)) )/dzh(1)
      if(wstar3.gt.0.) then
        wstar2 = wstar3**(2./3.)
      else
        wstar2 = 0.
      endif
      usurf  = sqrt(u(1)*u(1)+v(1)*v(1)+wstar2)
      facth  = ch*usurf*dzh(1)/kh(1)
      factq  = cq*usurf*dzh(1)/kh(1)

      dia(1) = 1.+facth
      sup(1) = -1.
      rhs(1) = facth*tgrnd

      dia(n)  = 1.
      sub(n)  = 0.
      rhs(n) = ttop

      call trislv(sub,dia,sup,rhs,t,n)
c ----------------------------------------------------------------------
      dia(1) = 1.+factq
      rhs1(1)= factq*qgrnd
      rhs1(n) = qtop
      call trislv(sub,dia,sup,rhs1,q,n)

      return
      end

      subroutine uveqns(u0,v0,u,v,z,km,dz,dzh,
     2                  ustar,cm,z0m,utop,vtop,dtime,coriol,
     3                  ug,vg,ilong,jlat,n)
c ----------------------------------------------------------------------
c this routine computes the matrices for the solutions of the u and v
c  fields as well as appling the boundary conditions.
c  the boundary conditions at the bottom are:
c
c     km * du/dz = cm * usurf * u
c     km * dv/dz = cm * usurf * v, evaluated at the lowest level.
c
c  at the top, the winds are prescribed.
c
c the arrays au and av are the lhs's of the appropriate discretized
c  equation. the vectors bu and bv are similarly the rhs's.
c  nb: u0 and v0 are the velocity components from the previous time
c      step and u and v are (a) the same as u0 and v0 for the
c      computation of the predictors, then are (b) the predictors,
c      i.e., the current iterate of the computed solution.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,k,o-z)
      parameter (kappa=0.40,grav=9.81,epslon=1.e-40,nm=8)
      PARAMETER (IM=72,JM=46)
      dimension u0(n),v0(n),u(n),v(n),z(n)
      dimension km(n-1),dz(n),dzh(n-1)
      common /trid_b/sub(nm),dia(nm),sup(nm),rhs(nm),rhs1(nm)
      common /pgpass/dpdxr(im,jm),dpdyr(im,jm),phi(im,jm)
     &              ,dpdxr0(im,jm),dpdyr0(im,jm)

c ----------------------------------------------------------------------
c compute the lhss of the prognostic equations. the (1,j) and (n,j)
c  components are zeroed to accomodate boundary conditions:

      do 200 i=2,n-1
          sub(i)=-dtime/(dz(i)*dzh(i-1))*km(i-1)
          sup(i)=-dtime/(dz(i)*dzh(i))*km(i)
          dia(i)=1.-(sub(i)+sup(i))
200   continue
c ----------------------------------------------------------------------
c now compute the rhss. components (1) and (n) are assigned as bcs:
      factx=(dpdxr(ilong,jlat)-dpdxr0(ilong,jlat))/(z(n)-z(1))
      facty=(dpdyr(ilong,jlat)-dpdyr0(ilong,jlat))/(z(n)-z(1))
      do 500 i=2,n-1
        dpdx=factx*(z(i)-z(1))+dpdxr0(ilong,jlat)
        dpdy=facty*(z(i)-z(1))+dpdyr0(ilong,jlat)
        rhs(i)=u0(i)+dtime*(coriol*v(i)-dpdx)
        rhs1(i)=v0(i)-dtime*(coriol*u(i)+dpdy)
c       rhs(i)=u0(i)+dtime*coriol*(v(i)-vg)
c       rhs1(i)=v0(i)-dtime*coriol*(u(i)-ug)
500   continue
c ----------------------------------------------------------------------
c finally, apply the boundary conditions that the flux be continuous
c  across the boundary at the surface (level 1 of the sub-grid scale
c  model, see comments at beginning of this routine). at the top, the
c  prognostic fields are set equal to those from the lowest gcm model
c  level:

      usurf  = sqrt(u(1)*u(1)+v(1)*v(1))
      factor = 1.+cm*usurf*dzh(1)/km(1)

      dia(1)= factor
      sup(1)= -1.
      rhs(1)  = 0.
c
      dia(n) = 1.
      sub(n) = 0.
      rhs(n)  = utop
c
      rhs1(1)  = 0.
      rhs1(n)  = vtop
c
      call trislv(sub,dia,sup,rhs,u,n)
      call trislv(sub,dia,sup,rhs1,v,n)
c ----------------------------------------------------------------------

      return
      end


      subroutine tqeqns_sta(u,v,t,q,kh,dz,dzh,
     2            ch,cq,tgrnd,qgrnd,ttop,qtop,n)
c ----------------------------------------------------------------------
c  computes the static solutions of the t and q
c  the boundary conditions at the bottom are:
c     kh * dt/dz = ch * usurf * (t - tg)
c     kh * dq/dz = cq * usurf * (q - qg), evaluated at the lowest level.
c  at the top, the temperature and moisture are prescribed.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,k,o-z)
      parameter (grav=9.81,nm=8)
      dimension u(n),v(n),t(n),q(n),kh(n-1),dz(n),dzh(n-1)
      common /trid_b/sub(nm),dia(nm),sup(nm),rhs(nm),rhs1(nm)
c ----------------------------------------------------------------------
      do 200 i=2,n-1
          sub(i)=-1./(dz(i)*dzh(i-1))*kh(i-1)
          sup(i)=-1./(dz(i)*dzh(i))*kh(i)
          dia(i)=-(sub(i)+sup(i))
200   continue
c ----------------------------------------------------------------------
      do 500 i=2,n-1
        rhs(i)=0.
        rhs1(i)=0.
500   continue
c ----------------------------------------------------------------------
c     M.J.Miller et al. 1992:
      wstar3=-1000.*grav*kh(1)*( 2.*(t(2)-t(1))/(t(2)+t(1))
     &                        -(q(2)-q(1)) )/dzh(1)
      if(wstar3.gt.0.) then
        wstar2 = wstar3**(2./3.)
      else
        wstar2 = 0.
      endif
      usurf  = sqrt(u(1)*u(1)+v(1)*v(1)+wstar2)
      facth  = ch*usurf*dzh(1)/kh(1)
      factq  = cq*usurf*dzh(1)/kh(1)
c
      dia(1) = 1.+facth
      sup(1) = -1.
      rhs(1) = facth*tgrnd

      dia(n)  = 1.
      sub(n)  = 0.
      rhs(n) = ttop
c
      call trislv(sub,dia,sup,rhs,t,n)
c ----------------------------------------------------------------------
      dia(1) = 1.+factq
      rhs1(1)= factq*qgrnd
      rhs1(n) = qtop
      call trislv(sub,dia,sup,rhs1,q,n)
c
      return
      end

      subroutine uveqns_sta(u,v,z,km,dz,dzh,
     2            ustar,cm,utop,vtop,coriol,ilong,jlat,n)
c ----------------------------------------------------------------------
c  computes the static solutions of the u and v
c  the boundary conditions at the bottom are:
c     km * du/dz = cm * usurf * u
c     km * dv/dz = cm * usurf * v, evaluated at the lowest level.
c  at the top, the winds are prescribed.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,k,o-z)
      parameter (nm=8)
      PARAMETER (IM=72,JM=46)
      dimension u(n),v(n),z(n),km(n-1),dz(n),dzh(n-1)
      common /trid_b/sub(nm),dia(nm),sup(nm),rhs(nm),rhs1(nm)
      common /pgpass/dpdxr(im,jm),dpdyr(im,jm),phi(im,jm)
     &              ,dpdxr0(im,jm),dpdyr0(im,jm)
c ----------------------------------------------------------------------
      do 200 i=2,n-1
          sub(i)=-1./(dz(i)*dzh(i-1))*km(i-1)
          sup(i)=-1./(dz(i)*dzh(i))*km(i)
          dia(i)=-(sub(i)+sup(i))
200   continue
c ----------------------------------------------------------------------
      factx=(dpdxr(ilong,jlat)-dpdxr0(ilong,jlat))/(z(n)-z(1))
      facty=(dpdyr(ilong,jlat)-dpdyr0(ilong,jlat))/(z(n)-z(1))
      write(99,*) factx,facty
      do 500 i=2,n-1
        dpdx=factx*(z(i)-z(1))+dpdxr0(ilong,jlat)
        dpdy=facty*(z(i)-z(1))+dpdyr0(ilong,jlat)
        rhs(i)=(coriol*v(i)-dpdx)
        rhs1(i)=-(coriol*u(i)+dpdy)
c       rhs(i)=coriol*(v(i)-vg)
c       rhs1(i)=-coriol*(u(i)-ug)
500   continue
c ----------------------------------------------------------------------
      usurf  = sqrt(u(1)*u(1)+v(1)*v(1))
      factor = 1.+cm*usurf*dzh(1)/km(1)
c
      dia(1)= factor
      sup(1)= -1.
      rhs(1)  = 0.
c
      dia(n) = 1.
      sub(n) = 0.
      rhs(n)  = utop
c
      rhs1(1)  = 0.
      rhs1(n)  = vtop
c
      call trislv(sub,dia,sup,rhs,u,n)
      call trislv(sub,dia,sup,rhs1,v,n)
c
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     from numerical recipes:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getk_old(km,kh,gm,gh,u,v,t,e,lscale,z,zhat,dzh,itype,n)
c ----------------------------------------------------------------------
c This routine computes the turbulent viscosity, KM, and turbulent
c  diffusivity, KH. These coefficients are computed using the second
c  order closure model of Galperin et al. (1988) and are modified under
c  stable conditions to give proper scaling at large Ri. Two grids are
c  used in the computation. The quantities U, V, T, and Z are computed
c  on the primary grid. The quantities E, LSCALE, KM, KH, GM, and GH
c  are computed on the secondary grid, which is staggered with respect
c  to the primary grid.
c
c Requiring GM be positive results in the requirement GH .le. 0.0233.
c  The actual limit applied here is 75% of that too make for a cleaner
c  solution. Actually, this limit is rarely used.
c
c The quantities 1.5E-5 and 2.5E-5 m**2/sec added to km and kh below
c  are the viscosity and thermometric conductivity of air,
c  respectively.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,k,o-z)
      dimension km(n-1),kh(n-1),gm(n-1),gh(n-1)
      dimension u(n),v(n),t(n),z(n),zhat(n-1),e(n-1)
      dimension dzh(n-1)
      real*8 lscale(n-1)
      parameter (grav=9.81)
      parameter (a1=0.92,b1=16.6,c1=0.08,a2=0.74,b2=10.1,sq=0.20)

      parameter (f1=-3.*a1*a2*b1*b2+18.*a2*b2*a1*a1+9.*a1*a2*a2*b1
     2              -54.*a1*a1*a2*a2+9.*a1*a2*b1*b2*c1
     3              +54.*a1*a1*a2*b1*c1)
      parameter (f2= 9.*a1*a2*a2*b1+27.*a1*a2*a2*b2+108.*a1*a1*a2*a2)
      parameter (f3= 21.*a1*a2+3.*a2*b2+a2*b1)
      parameter (f4= 6.*a1*a1-a1*b1+3.*a1*b1*c1)
      parameter (f5= a1*a1*b1*b1-12.*b1*a1**3+36.*a1**4
     2              -6.*a1*a1*b1*b1*c1+9.*a1*a1*b1*b1*c1*c1
     3              +36.*b1*c1*a1**3)
      parameter (f6= 216.*a1*a1*a2*a2-36.*a1*a2*a2*b1+252*a2*a1**3
     2              -90.*a1*a1*a2*b1*c1-30.*a1*a1*a2*b1
     3              -36.*a1*a1*a2*b2+6.*a1*a2*b1*b2
     4              -18.*a1*a2*b1*b2*c1-2.*a1*a2*b1*b1
     5              +6.*a1*a2*b1*b1*c1)
      parameter (f7= 18.*a1*a2*a2*b2+6.*a1*a2*a2*b1+9.*a2*a2*b2*b2
     2              +a2*a2*b1*b1+9.*a1*a1*a2*a2+6.*a2*a2*b1*b2)

      parameter (sm1=1.-3.*c1-6.*a1/b1)
      parameter (sm2=-3.*a2*((b2-3.*a2)*(1.-6.*a1/b1)-3.*c1*(b2+6.*a1)))
      parameter (sm3=9.*a1*a2)
      parameter (sh1=a2*(1.-6.*a1/b1))
      parameter (s1 =3.*a2*(6.*a1+b2))

      parameter (ri0=1.e-7,rimin=-0.44140)
      parameter (ghmin=-0.5*0.75*0.75,ghmax=0.75/(a2*(12.*a1+b1+3.*b2)))
c     parameter (rimaxm=0.045840,rimaxh=0.098800)
      parameter (rimaxm=0.0780715,rimaxh=0.089044)
      parameter (p1=1./3.,p2=4./3.)
      parameter (emin=1.e-2)

      do 100 i=1,n-1

        bvfrq2=grav*log(t(i+1)/t(i))/dzh(i)+1.e-8
        shear2=((u(i+1)-u(i))/dzh(i))**2+
     2         ((v(i+1)-v(i))/dzh(i))**2+1.e-8
        rich=bvfrq2/shear2
        if (abs(rich).lt.ri0) rich=sign(ri0,rich)

        if (rich.lt.rimin) then
          gh(i)=ghmax
          sm=(sm1+gh(i)*sm2)*a1/((1.-s1*gh(i))*(1.-sm3*gh(i)))
          sh=sh1/(1.-s1*gh(i))
        endif

        if ((rich.gt.rimin).and.(rich.lt.rimaxm)) then
          gh(i)=0.5*((f3*rich+f4+
     2             sqrt(f5+f6*rich+f7*rich*rich))/(f1+f2*rich))
          sm=(sm1+gh(i)*sm2)*a1/((1.-s1*gh(i))*(1.-sm3*gh(i)))
          sh=sh1/(1.-s1*gh(i))
        endif

        if (rich.gt.rimaxm) then
          gh(i)=0.5*((f3*rimaxm+f4+
     2             sqrt(f5+f6*rimaxm+f7*rimaxm*rimaxm))/(f1+f2*rimaxm))
          sm=(sm1+gh(i)*sm2)*a1/((1.-s1*gh(i))*(1.-sm3*gh(i)))
c         xm=rimaxm/rich
c         sm=sm*(xm**p1)
          xm=rich/rimaxm
          sm=sm*(2.**p1-1.)/((1.+xm)**p1-1.)
          if (rich.gt.rimaxh) then
            gh(i)=0.5*((f3*rimaxh+f4+
     2             sqrt(f5+f6*rimaxh+f7*rimaxh*rimaxh))/(f1+f2*rimaxh))
            sh=sh1/(1.-s1*gh(i))
c           xh=rimaxh/rich
c           sh=sh*(xh**p2)
            xh=rich/rimaxh
            sh=sh*(2.**p2-1.)/((1.+xh)**p2-1.)
            else
            gh(i)=0.5*((f3*rich+f4+
     2               sqrt(f5+f6*rich+f7*rich*rich))/(f1+f2*rich))
            sh=sh1/(1.-s1*gh(i))
          endif
        endif

        richf=sh*rich/sm
c         if (rich.gt.0.) then
            smmin=2.*emin/(b1*lscale(i)*lscale(i)*shear2*(1.-richf))
            if (sm.lt.smmin) then
              prandtl=sm/sh
              sm=smmin
              sh=sm/prandtl
            endif
c         endif
        gm(i)=-gh(i)/rich

        e(i)=0.5*b1*sm*lscale(i)*lscale(i)*shear2*(1.-richf)

        km(i) = sqrt(b1*sm*(1.-richf)*shear2)*
     2          sm*lscale(i)*lscale(i)+1.5e-5
        kh(i) = sqrt(b1*sm*(1.-richf)*shear2)*
     2          sh*lscale(i)*lscale(i)+2.5e-5

        if (km(i).gt.100.) then
          prndtl=km(i)/kh(i)
          km(i)=100.
          kh(i)=km(i)/prndtl
        endif

100   continue

      do 200 i=2,n-2
        kmtest=0.25*(km(i-1)+km(i+1))
        khtest=0.25*(kh(i-1)+kh(i+1))
        if (km(i).lt.kmtest) km(i)=kmtest
        if (kh(i).lt.khtest) kh(i)=khtest
200   continue

      if (km(1).lt.0.5*km(2)) km(1)=0.5*km(2)
      if (kh(1).lt.0.5*kh(2)) kh(1)=0.5*kh(2)
      if (km(n-1).lt.0.5*km(n-2)) km(n-1)=0.5*km(n-2)
      if (kh(n-1).lt.0.5*kh(n-2)) kh(n-1)=0.5*kh(n-2)

      return
      end

      subroutine elevl2(e,u,v,t,km,kh,lscale,dzh,ustar,n)
c ----------------------------------------------------------------------
c This routine computes the turbulence energy using the
c  Level 2 prescription.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,k,o-z)
      parameter (b1=16.6,grav=9.81,emax=8.)
      real*8 lscale(n-1)
      dimension e(n-1),u(n),v(n),t(n)
      dimension km(n-1),kh(n-1),dzh(n-1)

      do 100 i=1,n-1
        shear2=(((u(i+1)-u(i))/dzh(i))**2+
     2          ((v(i+1)-v(i))/dzh(i))**2)
        tgrad =log(t(i+1)/t(i))/dzh(i)
        extra=b1*lscale(i)*(km(i)*shear2-kh(i)*grav*tgrad)
        if (extra.lt.1.e-3) extra=1.e-3
        e(i)=(0.125*extra*extra)**(1./3.)
        if (e(i).gt.emax) e(i)=emax
100   continue

      return
      end

      subroutine pblini
c -------------------------------------------------------------
c These routines include the array ipbl which indicates if the
c  computation for a particular ITYPE was done last time step.
c Sets up the initialization of wind, temperature, and moisture
c  fields in the boundary layer. The initial values of these
c  fields are obtained by solving the static equations of the
c  Level 2 model. This is used when starting from a restart
c  file that does not have this data stored.
c -------------------------------------------------------------
      USE E001M12_COM
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (zgs=10.,ohmega=7.292e-5,rvx=0.)
      parameter (npbl=8)
      parameter (iq1=im/4+1,iq2=im/2+1,iq3=3*im/4+1)
      dimension tgvdat(im,jm,4)
      common /rdata/roughl(im,jm)
      common /pblvar/uinit(npbl),vinit(npbl),tinit(npbl),qinit(npbl),
     2               einit(npbl)
      common /socabl/uabl(npbl,im,jm,4),vabl(npbl,im,jm,4),
     2               tabl(npbl,im,jm,4),qabl(npbl,im,jm,4),
     3               eabl(npbl,im,jm,4),
     4               cmgs(im,jm,4),chgs(im,jm,4),cqgs(im,jm,4),
     5               ipbl(im,jm,4)
      qsat(tm,pm,qlh)=3.797915*exp(qlh*(7.93252d-6-2.166847d-3/tm))/pm

      call pgrads1   !   added 6/19/00
      do 60 j=1,jm
        do 50 i=1,im
          pland=fdata(i,j,2)
          pwater=1.-pland
          plice=fdata(i,j,3)*pland
          psoil=pland-plice
          poice=odata(i,j,2)*pwater
          pocean=pwater-poice
          tgvdat(i,j,1)=odata(i,j,1) +273.16
          if (pocean.le.0.) tgvdat(i,j,1)=0.
          tgvdat(i,j,2)=gdata(i,j,3) +273.16
          if (poice.le.0.)  tgvdat(i,j,2)=0.
          tgvdat(i,j,3)=gdata(i,j,13)+273.16
          if (plice.le.0.)  tgvdat(i,j,3)=0.
          tgvdat(i,j,4)=gdata(i,j,4) +273.16
          if (psoil.le.0.)  tgvdat(i,j,4)=0.
   50   continue
   60 continue

      call readt (19,0,roughl,im*jm,roughl,1)
      rewind 19
      pi=dacos(-1.d0)
      radian=pi/180.

      call getb(zgs,ztop,bgrid)

      do 400 itype=1,4
        if ((itype.eq.1).or.(itype.eq.4)) then
          elhx=lhe
          else
          elhx=lhs
        endif
        do 300 j=1,jm
          if ((j.eq.1).or.(j.eq.jm)) then
            imax=1
            else
            imax=im
          endif
          jlat=j
          phi=radian*(float(j-1)*180./float(jm-1)-90.)
          coriol=2.*sin(phi)*ohmega

          im1=im
          do 200 i=1,imax
            tgrnd=tgvdat(i,j,itype)
            if (tgrnd.eq.0.) then
              ipbl(i,j,itype)=0
              im1=i
              go to 200
            endif
            ilong=i
            pij=p(i,j)
            ps=pij+ptop
            psk=expbyk(ps)
            qgrnd=qsat(tgrnd,ps,elhx)

            if (j.eq.1) then
c ******************************************************************
c           At the south pole:
              jvpo=2
              hemi=-1.
              utop=.25*(u(1,jvpo,1)-u(iq2,jvpo,1)
     2                -(v(iq1,jvpo,1)-v(iq3,jvpo,1))*hemi)
              vtop=.25*(v(1,jvpo,1)-v(iq2,jvpo,1)
     2                +(u(iq1,jvpo,1)-u(iq3,jvpo,1))*hemi)
c ******************************************************************
            endif

            if (j.eq.jm) then
c ******************************************************************
c     At the north pole:
              jvpo=jm
              hemi=1.
              utop=.25*(u(1,jvpo,1)-u(iq2,jvpo,1)
     2                -(v(iq1,jvpo,1)-v(iq3,jvpo,1))*hemi)
              vtop=.25*(v(1,jvpo,1)-v(iq2,jvpo,1)
     2                +(u(iq1,jvpo,1)-u(iq3,jvpo,1))*hemi)
c ******************************************************************
            endif

            if ((j.gt.1).and.(j.lt.jm)) then
c ******************************************************************
c     Away from the poles:
              utop=.25*(u(im1,j,1)+u(i,j,1)+u(im1,j+1,1)+u(i,j+1,1))
              vtop=.25*(v(im1,j,1)+v(i,j,1)+v(im1,j+1,1)+v(i,j+1,1))
c ******************************************************************
            endif

            qtop=q(i,j,1)
            ttop=t(i,j,1)*(1.+qtop*rvx)*psk
            if (itype.gt.2) then
              zgrnd=30./(10.**roughl(i,j))
              else
              zgrnd=0.1
            endif
            call inits(tgrnd,qgrnd,zgrnd,zgs,ztop,utop,vtop,
     2                 ttop,qtop,coriol,cm,ch,cq,bgrid,ustar,
     3                 ilong,jlat,itype)
            cmgs(i,j,itype)=cm
            chgs(i,j,itype)=ch
            cqgs(i,j,itype)=cq

            do 100 lpbl=1,npbl
              uabl(lpbl,i,j,itype)=uinit(lpbl)
              vabl(lpbl,i,j,itype)=vinit(lpbl)
              tabl(lpbl,i,j,itype)=tinit(lpbl)
              qabl(lpbl,i,j,itype)=qinit(lpbl)
 100        continue

            do 150 lpbl=1,npbl-1
              eabl(lpbl,i,j,itype)=einit(lpbl)
 150        continue

            ipbl(i,j,itype)=1
            bldata(i,j,8+itype)=ustar
            im1=i

 200      continue
 300    continue
c     write (99,1000) itype
 400  continue

      return
 1000 format (1h ,//,1x,'completed initialization, itype = ',i2,//)
      end

      subroutine inits(tgrnd,qgrnd,zgrnd,zgs,ztop,utop,vtop,
     2                 ttop,qtop,coriol,cm,ch,cq,bgrid,ustar,
     3                 ilong,jlat,itype)
c ----------------------------------------------------------------------
c This routine initializes the winds, temperature, and humidity using
c  static solutions of the Level 2 turbulence equations.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,k,o-z)
      parameter (grav=9.81,kappa=0.40,b1=16.6)
      parameter (n=8,itmax=100,w=0.50,tol=1.e-4,epslon=1.e-20)
      parameter (iprint= 0,jprint=25)  ! set iprint=39 e.g..to debug
      dimension u(n),v(n),t(n),q(n),e(n-1)
      dimension km(n-1),kh(n-1),ke(n-1),gm(n-1),gh(n-1)
      dimension usave(n),vsave(n)
      dimension tsave(n),qsave(n),esave(n)
      dimension z(n),zhat(n-1),xi(n),xihat(n-1)
      dimension dz(n),dzh(n-1)
      real*8 lscale(n-1),lmonin
      common /pblvar/u,v,t,q,e
c
c
      pi=dacos(-1.d0)
      radian=pi/180.
      z0m=zgrnd
      z0h=z0m
      z0q=z0m
      write(99,*) "inside inits"
      call griddr(z,zhat,xi,xihat,dz,dzh,zgs,ztop,bgrid,n)

c ----------------------------------------------------------------------
c Initialization for iteration:
      if (coriol.le.0.) then
        hemi=-1.
        else
        hemi= 1.
      endif
      if (tgrnd.gt.ttop) then
        psi1=hemi*5.*radian
        else
        psi1=hemi*15.*radian
      endif
      if (utop.eq.0.) utop=epslon
      if (vtop.eq.0.) vtop=epslon
      if ((utop.gt.0.).and.(vtop.gt.0.)) then
        psi0=atan2(vtop,utop)
      endif
      if ((utop.lt.0.).and.(vtop.gt.0.)) then
        psi0=atan2(abs(utop),vtop)+0.5*pi
      endif
      if ((utop.lt.0.).and.(vtop.lt.0.)) then
        psi0=atan2(abs(vtop),abs(utop))+pi
      endif
      if ((utop.gt.0.).and.(vtop.lt.0.)) then
        psi0=atan2(utop,abs(vtop))+1.5*pi
      endif
      psi=psi0+psi1
      usurf=0.4*sqrt(utop*utop+vtop*vtop)
      u(1)=usurf*cos(psi)
      v(1)=usurf*sin(psi)
      t(1)=tgrnd-0.5*(tgrnd-ttop)
      q(1)=qgrnd-0.5*(qgrnd-qtop)
      e(1)=1.e-2
        usave(1)=u(1)
        vsave(1)=v(1)
        tsave(1)=t(1)
        qsave(1)=q(1)
        esave(1)=e(1)
      u(n)=utop
      v(n)=vtop
      t(n)=ttop
      q(n)=qtop
      e(n-1)=2.e-2

      do 100 i=2,n-1
        u(i)=u(1)+(u(n)  -u(1))*((z(i)-z(1))/(z(n)  -z(1)))
        v(i)=v(1)+(v(n)  -v(1))*((z(i)-z(1))/(z(n)  -z(1)))
        t(i)=t(1)+(t(n)  -t(1))*((z(i)-z(1))/(z(n)  -z(1)))
        q(i)=q(1)+(q(n)  -q(1))*((z(i)-z(1))/(z(n)  -z(1)))
        e(i)=e(1)+(e(n-1)-e(1))*((z(i)-z(1))/(z(n-1)-z(1)))
        usave(i)=u(i)
        vsave(i)=v(i)
        tsave(i)=t(i)
        qsave(i)=q(i)
        esave(i)=e(i)
100   continue

      call getl1(e,zhat,dzh,lscale,n)
      call getk(km,kh,ke,gm,gh,u,v,t,e,lscale,dzh,n)
      call stars(ustar,tstar,qstar,lmonin,tgrnd,qgrnd,
     2           u,v,t,q,z,z0m,z0h,z0q,cm,ch,cq,
     3           km,kh,dzh,itype,n)

c     ustar0=ustar
      tstar0=tstar

      if ((ilong.eq.iprint).and.(jlat.eq.jprint)) then
        iter=-1
        call output(u,v,t,q,e,lscale,z,zhat,dzh,
     2              km,kh,gm,gh,cm,ch,cq,z0m,z0h,z0q,
     3              ustar,tstar,qstar,lmonin,tgrnd,qgrnd,
     4              utop,vtop,ttop,qtop,
     5              dtime,bgrid,ilong,jlat,iter,itype,n)
      endif

c ----------------------------------------------------------------------

      do 300 iter=1,itmax
        call tqeqns_sta(u,v,t,q,kh,dz,dzh,
     2            ch,cq,tgrnd,qgrnd,ttop,qtop,n)
        call uveqns_sta(u,v,z,km,dz,dzh,
     2            ustar,cm,utop,vtop,coriol,ilong,jlat,n)

        call tcheck(t,tgrnd,n)
        call tcheck(q,qgrnd,n)
        call ucheck(u,v,z,ustar,lmonin,z0m,hemi,psi0,psi1,n)

        call elevl2(e,u,v,t,km,kh,lscale,dzh,ustar,n)

        do 200 i=1,n-1
          u(i)=w*usave(i)+(1.-w)*u(i)
          v(i)=w*vsave(i)+(1.-w)*v(i)
          t(i)=w*tsave(i)+(1.-w)*t(i)
          q(i)=w*qsave(i)+(1.-w)*q(i)
          e(i)=w*esave(i)+(1.-w)*e(i)
          usave(i)=u(i)
          vsave(i)=v(i)
          tsave(i)=t(i)
          qsave(i)=q(i)
          esave(i)=e(i)
 200    continue

        call getk_old(km,kh,gm,gh,u,v,t,e,lscale,z,zhat,dzh,itype,n)
        call getl1(e,zhat,dzh,lscale,n)
        call stars(ustar,tstar,qstar,lmonin,tgrnd,qgrnd,
     2             u,v,t,q,z,z0m,z0h,z0q,cm,ch,cq,
     3             km,kh,dzh,itype,n)

c       test=abs(ustar-ustar0)/(ustar+ustar0)
        test=abs(tstar-tstar0)/abs(tstar+tstar0)
        if (test.lt.tol) go to 400
c       ustar0=ustar
        tstar0=tstar

        if ((ilong.eq.iprint).and.(jlat.eq.jprint)) then
          call output(u,v,t,q,e,lscale,z,zhat,dzh,
     2                km,kh,gm,gh,cm,ch,cq,z0m,z0h,z0q,
     3                ustar,tstar,qstar,lmonin,tgrnd,qgrnd,
     4                utop,vtop,ttop,qtop,
     5                dtime,bgrid,ilong,jlat,iter,itype,n)
        endif

300   continue

c     write (99,9900) ilong,jlat,itype,test
400   continue

        if ((ilong.eq.iprint).and.(jlat.eq.jprint)) then
        call output(u,v,t,q,e,lscale,z,zhat,dzh,
     2              km,kh,gm,gh,cm,ch,cq,z0m,z0h,z0q,
     3              ustar,tstar,qstar,lmonin,tgrnd,qgrnd,
     4              utop,vtop,ttop,qtop,
     5              dtime,bgrid,ilong,jlat,iter,itype,n)
      endif

c     call check1(u,n,ilong,jlat,1)
c     call check1(v,n,ilong,jlat,2)
c     call check1(t,n,ilong,jlat,3)
c     call check1(q,n,ilong,jlat,4)
c     call check1(ustar,1,ilong,jlat,0)

      return
1000  format (1h ,i3,10(1x,1pe11.4))
2000  format (1h ,/,1x,'iter = ',i3,/)
3000  format (1h ,'test = ',1pe11.4,/)
8000  format (1h ,9(1pe11.4,1x),1x,1pe10.3,2x,1pe10.3)
9000  format (1h )
9900  format (1h ,'i = ',i2,2x,'j = ',i2,2x,'itype = ',i2,2x,
     2            'test = ',1pe11.4)
      end

      subroutine tcheck(t,tgrnd,n)
c ----------------------------------------------------------------------
c This routine makes sure that the temperature remains within
c  reasonable bounds during the initialization process. (Sometimes the
c  the computed temperature iterated out in left field someplace,
c  *way* outside any reasonable range.) This routine keeps the temp
c  between the maximum and minimum of the boundary temperatures.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,k,l,o-z)
      dimension t(n)
      dimension imax(20),imin(20)

      if (tgrnd.lt.t(n)) then
        tmin=tgrnd
        tmax=t(n)
        else
        tmin=t(n)
        tmax=tgrnd
      endif
      nmin=0
      nmax=0

      do 100 i=1,n-1
        if (t(i).lt.tmin) then
          nmin=nmin+1
          imin(nmin)=i
        endif
        if (t(i).gt.tmax) then
          nmax=nmax+1
          imax(nmax)=i
        endif
100   continue

      dt=0.01*(tmax-tmin)
      if (nmin.gt.0) then
        do 200 i=1,nmin
          t(imin(i))=tmin+float(i)*dt
200     continue
      endif

      if (nmax.gt.0) then
        do 400 i=1,nmax
          t(imax(i))=tmax-float(i)*dt
400     continue
      endif

      return
      end

      subroutine ucheck(u,v,z,ustar,lmonin,z0m,hemi,psi0,psi1,n)
c ----------------------------------------------------------------------
c This routine makes sure that the winds remain within reasonable
c  bounds during the initialization process. (Sometimes the computed
c  wind speed iterated out in left field someplace, *way* outside
c  any reasonable range.) Tests and corrects both direction and
c  magnitude of the wind rotation with altitude. Tests the total wind
c  speed via comparison to similarity theory. Note that it works from
c  the top down so that it can assume that at level (i), level (i+1)
c  displays reasonable behavior.
c ----------------------------------------------------------------------
      implicit real*8 (a-h,k,l,o-z)
      parameter (epslon=1.e-20)
      parameter (pi=3.141592654,radian=pi/180.)
      parameter (psistb=15.*radian,psiuns=5.*radian)
      parameter (gammau=19.3,gammas=4.8)
      parameter (kappa=0.40)
      dimension u(n),v(n),z(n)

      if (lmonin.ge.0.) then
        psilim=psistb
        else
        psilim=psiuns
      endif

c --------------------------------------------------------------------
c First, check the rotation of the wind vector:

      do 100 i=n-1,1,-1
        if (u(i).eq.0.) u(i)=epslon
        if (v(i).eq.0.) v(i)=epslon
        if ((u(i).gt.0.).and.(v(i).gt.0.)) then
          psiu=atan2(v(i),u(i))
        endif
        if ((u(i).lt.0.).and.(v(i).gt.0.)) then
          psiu=atan2(abs(u(i)),v(i))+0.5*pi
        endif
        if ((u(i).lt.0.).and.(v(i).lt.0.)) then
          psiu=atan2(abs(v(i)),abs(u(i)))+pi
        endif
        if ((u(i).gt.0.).and.(v(i).lt.0.)) then
          psiu=atan2(u(i),abs(v(i)))+1.5*pi
        endif
        psirot=psiu-psi0
c --------------------------------------------------------------------
c Next, check the magnitude of the wind. If the gradient is incorrect,
c  set the wind magnitude to that given by similarity theory:

        utotal=sqrt(u(i)*u(i)+v(i)*v(i))
        utest =sqrt(u(i+1)*u(i+1)+v(i+1)*v(i+1))
        if (utotal.gt.utest) then
          zbyl=z(i)/lmonin
          z0byl=z0m/lmonin
          if (lmonin.gt.0.) then
            dpsim=-gammas*(zbyl-z0byl)
            else
            x  = (1.-gammau* zbyl)**0.25
            x0 = (1.-gammau*z0byl)**0.25
            dpsim=log((1.+x )*(1.+x )*(1.+x *x )/
     2               ((1.+x0)*(1.+x0)*(1.+x0*x0)))-
     3            2.*(atan(x)-atan(x0))
          endif
          utotal=(ustar/kappa)*(log(z(i)/z0m)-dpsim)
          if (utotal.gt.utest) utotal=0.95*utest
          u(i)=utotal*cos(psiu)
          v(i)=utotal*sin(psiu)
        endif

        if (hemi*psirot.lt.0.) then
          angle=psi0+psi1*float(n-i)/float(n-1)
          u(i)=utotal*cos(angle)
          v(i)=utotal*sin(angle)
          go to 100
        endif

        if (hemi*psirot.gt.psilim) then
          angle=psi0+hemi*psilim*float(n-i)/float(n-1)
          u(i)=utotal*cos(angle)
          v(i)=utotal*sin(angle)
        endif

100   continue
c --------------------------------------------------------------------

      return
      end


      subroutine check1(a,n,ilong,jlat,id)
      implicit real*8 (a-h,o-z)
      dimension a(n)
c ----------------------------------------------------------------------
c checks for NaN'S and INF'S in real 1-D arrays.
c input:
c a      - real array
c n      - dimension of a
c ilong  - longitude identifier
c jlat   - latitude identifier
c id     - integer id
c output:
c prints to standard output array element number and id
c ----------------------------------------------------------------------
      character*16 str

      do 100 i=1,n
        write(str,'(e16.8)')a(i)
        k=index(str,'N')+index(str,'n')
        if (k.ne.0) then
          write (99,1000) ilong,jlat,i,id,a(i)
          if (id.lt.100) stop 'check1'
        endif
  100 continue

      return
 1000 format (1h ,'check1:  ilong = ',6x,i3,/,1x,
     2            '         jlat  = ',6x,i3,/,1x,
     3            '         i     = ',6x,i3,/,1x,
     4            '         id    = ',6x,i3,/,1x,
     5            '         value = ',1pe11.4,/)
      end

      subroutine output(u,v,t,q,e,lscale,z,zhat,dzh,
     2                  km,kh,gm,gh,cm,ch,cq,z0m,z0h,z0q,
     3                  ustar,tstar,qstar,lmonin,tgrnd,qgrnd,
     4                  utop,vtop,ttop,qtop,
     5                  dtime,bgrid,ilong,jlat,iter,itype,n)
      implicit real*8 (a-h,k,o-z)
      parameter (epslon=1.e-40,degree=180./3.141592654,kappa=0.40)
      parameter (sigmat=0.95,gammam=19.3,gammah=11.6)
      parameter (betam=4.8,betah=8./sigmat)
      parameter (grav=9.81)
      real*8 lscale(n-1),lmonin
      dimension u(n),v(n),t(n),q(n),e(n-1),z(n),zhat(n-1),dzh(n-1)
      dimension km(n-1),kh(n-1),gm(n-1),gh(n-1)
c ----------------------------------------------------------------------
c Output for diagnostic purposes:
      write (99,5000) ilong,jlat,itype,dtime,iter,ustar,tstar,qstar,
     2                lmonin,tgrnd,qgrnd,cm,ch,cq,z0m,z0h,z0q,
     3                utop,vtop,ttop,qtop,
     4                bgrid
      write (99,1000)
      psitop=atan2(v(n),u(n)+epslon)*degree
      if (psitop.lt.0.) psitop=psitop+360.
      do 100 i=1,n
        psi=atan2(v(i),u(i)+epslon)*degree
        if (psi.lt.0.) psi=psi+360.
        psi=psi-psitop
        utotal=sqrt(u(i)*u(i)+v(i)*v(i))
        call simil(utest,ttest,qtest,z(i),ustar,tstar,qstar,
     2             z0m,z0h,z0q,lmonin,tgrnd,qgrnd)
        write (99,3000) i,z(i),u(i),v(i),psi,utotal,utest,
     2                    t(i),ttest,q(i),qtest
100   continue
      write (99,9000)
      write (99,2000)
      do 200 i=1,n-1
        utot1=u(i+1)*u(i+1)+v(i+1)*v(i+1)
        utot2=u(i)  *u(i)  +v(i)  *v(i)
        if (utot1.ge.utot2) then
          sign=1.
          else
          sign=-1.
        endif
        shear =sqrt((u(i+1)-u(i))*(u(i+1)-u(i))
     2             +(v(i+1)-v(i))*(v(i+1)-v(i)))/dzh(i)
        tgrad =(t(i+1)-t(i))/dzh(i)
        qgrad =(q(i+1)-q(i))/dzh(i)
        egrad =(e(i+1)-e(i))/dzh(i)
        uflux=km(i)*shear*sign
        hflux=kh(i)*tgrad
        qflux=kh(i)*qgrad
        bvfrq2=grav*log(t(i+1)/t(i))/dzh(i)+1.e-12
        shear2=((u(i+1)-u(i))/dzh(i))**2+
     2         ((v(i+1)-v(i))/dzh(i))**2+1.e-8
        rich=bvfrq2/shear2
        write (99,3000) i,zhat(i),e(i),lscale(i),km(i),kh(i),
     2                    gm(i),gh(i),uflux,hflux,qflux,rich
200   continue
      write (99,9000)
      write (99,7000)
      do 300 i=1,n-1
        utot1=sqrt(u(i+1)*u(i+1)+v(i+1)*v(i+1))
        utot2=sqrt(  u(i)*  u(i)+  v(i)*  v(i))
        dudz=(utot1-utot2)/dzh(i)
        dtdz=(t(i+1)-t(i))/dzh(i)
        dqdz=(q(i+1)-q(i))/dzh(i)
        if (lmonin.lt.0.) then
          phim = 1./((1.-gammam*zhat(i)/lmonin)**0.25)
          phih = sigmat/sqrt(1.-gammah*zhat(i)/lmonin)
          else
          phim = 1.+betam*zhat(i)/lmonin
          phih = sigmat*(1.+betah*zhat(i)/lmonin)
        endif
        dudzs=ustar*phim/(kappa*zhat(i))
        dtdzs=tstar*phih/(kappa*zhat(i))
        dqdzs=qstar*phih/(kappa*zhat(i))
        uratio=dudzs/(dudz+1.e-40)
        tratio=dtdzs/(dtdz+1.e-40)
        qratio=dqdzs/(dqdz+1.e-40)

        dbydzh=1./dzh(i)
        shear2=dbydzh*dbydzh*((u(i+1)-u(i))**2+
     2                        (v(i+1)-v(i))**2)
        tgradl=dbydzh*log(t(i+1)/t(i))
        prod=km(i)*shear2-grav*kh(i)*tgradl

        write (99,3000) i,zhat(i),dudz,dudzs,dtdz,dtdzs,dqdz,dqdzs,
     2                            uratio,tratio,qratio,prod
300   continue
      write (99,9000)
      write (99,9000)
c ----------------------------------------------------------------------
      return
1000  format (1h ,'   i','      z     ','      U     ','      V     ',
     2                   '     psi    ','    Utotal  ','    Utest   ',
     2                   '      T     ','    Ttest   ','      Q     ',
     2                   '    Qtest   ',/)
2000  format (1h ,'   i','     zhat   ','      E     ','      L     ',
     2                   '      Km    ','      Kh    ',
     3                   '      Gm    ','      Gh    ','    uflux   ',
     4                   '    hflux   ','    qflux   ','    rich #  ',/)
3000  format (1h ,i4,11(1x,1pe11.4))
5000  format (1h ,'i      = ',9x,i2,/,1x,
     2            'j      = ',9x,i2,/,1x,
     2            'itype  = ',9x,i2,/,1x,
     2            'dtime  = ',1pe11.4,/,1x,
     3            'iter   = ',7x,i4,/,1x,
     4            'ustar  = ',1pe11.4,/,1x,
     5            'tstar  = ',1pe11.4,/,1x,
     5            'qstar  = ',1pe11.4,/,1x,
     5            'lmonin = ',1pe11.4,/,1x,
     5            'tgrnd  = ',1pe11.4,/,1x,
     5            'qgrnd  = ',1pe11.4,/,1x,
     6            'cm     = ',1pe11.4,/,1x,
     6            'ch     = ',1pe11.4,/,1x,
     6            'cq     = ',1pe11.4,/,1x,
     6            'z0m    = ',1pe11.4,/,1x,
     7            'z0h    = ',1pe11.4,/,1x,
     8            'z0q    = ',1pe11.4,/,1x,
     6            'utop   = ',1pe11.4,/,1x,
     6            'vtop   = ',1pe11.4,/,1x,
     7            'ttop   = ',1pe11.4,/,1x,
     8            'qtop   = ',1pe11.4,/,1x,
     8            'bgrid  = ',1pe11.4,/)
7000  format (1h ,'   i','     zhat   ','     dudz   ','   dudz sim ',
     2                   '     dtdz   ','   dtdz sim ','     dqdz   ',
     3                   '   dqdz sim ','    uratio  ','    tratio  ',
     4                   '    qratio  ','  production',/)
9000  format (1h )
      end

