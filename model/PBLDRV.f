      SUBROUTINE PBL(I,J,ITYPE,PTYPE)
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
      USE CONSTANT, only :  rgas,grav
      USE E001M12_COM
     &     , only : IM,JM,LM, t,q,u,v,p,ptop,ls1,psf
      USE DYNAMICS, only : pmid,pk,pedn
      USE GEOM, only : idij,idjj,kmaxj
      USE PBLCOM, ustar_type=>ustar
      USE SOCPBL, gravx=>grav
     &     ,uij=>u,vij=>v,tij=>t,qij=>q,eij=>e
     &     ,dpdxrij=>dpdxr,dpdyrij=>dpdyr
     &     ,dpdxr0ij=>dpdxr0,dpdyr0ij=>dpdyr0
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: I,J  !@var I,J grid point
      INTEGER, INTENT(IN) :: ITYPE  !@var ITYPE surface type
      REAL*8, INTENT(IN) :: PTYPE  !@var PTYPE percent surface type
C
C        ocean and ocean ice are treated as rough surfaces
C        roughness lengths from Brutsaert for rough surfaces
C

      REAL*8 KMSURF,KHSURF
      LOGICAL POLE

      REAL*8 ZS1,TGV,TKV,QG,HEMI,DTSURF
      REAL*8 US,VS,WS,TSV,QS,PSI,DBL,KM,KH,USTAR,PPBL,
     2               CM,CH,CQ,UG,VG,WG,ZMIX, TS

      COMMON /PBLPAR/ZS1,TGV,TKV,QG,HEMI,DTSURF,POLE

      COMMON /PBLOUT/US,VS,WS,TSV,QS,PSI,DBL,KM,KH,PPBL,
     2               UG,VG,WG,ZMIX

      REAL*8, PARAMETER ::  EPSLON=1.D-20,radian=3.141592654/180.
      REAL*8 Z0M,ztop,zpbl,pl1,tl1,pl,tl,tbar,thbar,zpbl1,coriol
      REAL*8 ttop,qtop,tgrnd,qgrnd,utop,vtop,z0h,z0q,ufluxs,vfluxs
     *     ,tfluxs,qfluxs,psitop,psisrf
      INTEGER LDC,jvpo,L,k

      INTEGER, DIMENSION(IM) :: IDI,IDJ
      INTEGER :: KMAX

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
        USTAR=USTAR_TYPE(I,J,ITYPE)
        DBL=0.3*USTAR/OMEGA2
        if (dbl.le.ztop) then
C THE VERTICAL LEVEL FOR WHICH WG IS COMPUTED IS THE FIRST:
          dbl=ztop
          L=1
          ELSE
          if (dbl.gt.3000.) dbl=3000.
C FIND THE VERTICAL LEVEL NEXT HIGHER THAN DBL AND COMPUTE WG THERE:
          zpbl=ztop
          pl1=pmid(1,i,j)         ! pij*sig(1)+ptop
          tl1=t(i,j,1)*pk(1,i,j)  !expbyk(pl1)
          do l=2,ls1
            pl=pmid(l,i,j)        !pij*sig(l)+ptop
            tl=t(i,j,l)*pk(l,i,j) !expbyk(pl)
            tbar=thbar(tl1,tl)
            zpbl=zpbl+(rgas/grav)*tbar*(pl1-pl)/pl1
            if (zpbl.ge.dbl) go to 200
            pl1=pl
            tl1=tl
          end do
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
        LDC=DCLEV(I,J)
        IF (LDC.EQ.0) LDC=1
        if (ldc.eq.1) then
          dbl=ztop
          l=1
          else
          zpbl=ztop
          pl1=pmid(1,i,j)          !pij*sig(1)+ptop
          tl1=t(i,j,1)*pk(1,i,j)   !expbyk(pl1)
          zpbl1=ztop
          do l=2,ldc
            pl=pmid(l,i,j)         !pij*sig(l)+ptop
            tl=t(i,j,l)*pk(l,i,j)  !expbyk(pl)
            tbar=thbar(tl1,tl)
            zpbl=zpbl+(rgas/grav)*tbar*(pl1-pl)/pl1
            if (zpbl.ge.3000.) then
              zpbl=zpbl1
              go to 400
            endif
            pl1=pl
            tl1=tl
            zpbl1=zpbl
          end do
400       continue
          l=l-1
          dbl=zpbl
        endif
C**********************************************************************
      ENDIF

C *********************************************************************
      ppbl=pedn(l,i,j)      ! sige(l)*pij+ptop
c      if (l.gt.ls1) ppbl=sige(l)*psfmpt+ptop
c      phi=radian*(float(j-1)*180./float(jm-1)-90.)
      coriol=sin(radian*(float(j-1)*180./float(jm-1)-90.))*omega2
      ttop=tkv
      qtop=q(i,j,1)
      tgrnd=tgv
      qgrnd=qg

      KMAX=KMAXJ(J)
      do k=1,kmax
         idi(k)=idij(k,i,j)
         idj(k)=idjj(k,j)
      enddo

      if (pole) then
        jvpo = idj(1)
        utop=.25*(u(1,jvpo,1)-u(iq2,jvpo,1)
     2          -(v(iq1,jvpo,1)-v(iq3,jvpo,1))*hemi)
        vtop=.25*(v(1,jvpo,1)-v(iq2,jvpo,1)
     2          +(u(iq1,jvpo,1)-u(iq3,jvpo,1))*hemi)
        ug  =.25*(u(1,jvpo,l)-u(iq2,jvpo,l)
     2          -(v(iq1,jvpo,l)-v(iq3,jvpo,l))*hemi)
        vg  =.25*(v(1,jvpo,l)-v(iq2,jvpo,l)
     2          +(u(iq1,jvpo,l)-u(iq3,jvpo,l))*hemi)
      else
c        utop=.25*(u(im1,j,1)+u(i,j,1)+u(im1,j+1,1)+u(i,j+1,1))
c        vtop=.25*(v(im1,j,1)+v(i,j,1)+v(im1,j+1,1)+v(i,j+1,1))
c        ug  =.25*(u(im1,j,l)+u(i,j,l)+u(im1,j+1,l)+u(i,j+1,l))
c        vg  =.25*(v(im1,j,l)+v(i,j,l)+v(im1,j+1,l)+v(i,j+1,l))
        utop=.25*(u(idi(1),idj(1),1)+u(idi(2),idj(2),1)+
     &            u(idi(3),idj(3),1)+u(idi(4),idj(4),1))
        vtop=.25*(v(idi(1),idj(1),1)+v(idi(2),idj(2),1)+
     &            v(idi(3),idj(3),1)+v(idi(4),idj(4),1))
        ug  =.25*(u(idi(1),idj(1),l)+u(idi(2),idj(2),l)+
     &            u(idi(3),idj(3),l)+u(idi(4),idj(4),l))
        vg  =.25*(v(idi(1),idj(1),l)+v(idi(2),idj(2),l)+
     &            v(idi(3),idj(3),l)+v(idi(4),idj(4),l))
      endif

      uij(:)=uabl(:,i,j,itype)
      vij(:)=vabl(:,i,j,itype)
      tij(:)=tabl(:,i,j,itype)
      qij(:)=qabl(:,i,j,itype)
      eij(1:n-1)=eabl(1:n-1,i,j,itype)

      cm=cmgs(i,j,itype)
      ch=chgs(i,j,itype)
      cq=cqgs(i,j,itype)

      dpdxrij  = dpdxr(i,j)
      dpdyrij  = dpdyr(i,j)
      dpdxr0ij = dpdxr0(i,j)
      dpdyr0ij = dpdyr0(i,j)

      call advanc(us,vs,tsv,qs,kmsurf,khsurf,ustar,ug,vg,cm,ch,cq,
     2            z0m,z0h,z0q,coriol,utop,vtop,ttop,qtop,tgrnd,
     3            qgrnd,zgs,ztop,zmix,dtsurf,ufluxs,vfluxs,
     4            tfluxs,qfluxs,i,j,itype)

      uabl(:,i,j,itype)=uij(:)
      vabl(:,i,j,itype)=vij(:)
      tabl(:,i,j,itype)=tij(:)
      qabl(:,i,j,itype)=qij(:)
      eabl(1:n-1,i,j,itype)=eij(1:n-1)

      cmgs(i,j,itype)=cm
      chgs(i,j,itype)=ch
      cqgs(i,j,itype)=cq
      ipbl(i,j,itype)=1

      ws    =sqrt(us*us+vs*vs)
      wg    =sqrt(ug*ug+vg*vg)
      km    =kmsurf
      kh    =khsurf
      psitop=atan2(vg,ug+epslon)
      psisrf=atan2(vs,us+epslon)
      psi   =psisrf-psitop
      ustar_type(i,j,itype)=ustar
C ******************************************************************
      TS=TSV ! /(1.+QS*RVX) ! rvx=0 for now
      WSAVG(I,J)=WSAVG(I,J)+WS*PTYPE
      TSAVG(I,J)=TSAVG(I,J)+TS*PTYPE
      if(itype.ne.4) QSAVG(I,J)=QSAVG(I,J)+QS*PTYPE
      USAVG(I,J)=USAVG(I,J)+US*PTYPE
      VSAVG(I,J)=VSAVG(I,J)+VS*PTYPE
      TAUAVG(I,J)=TAUAVG(I,J)+CM*WS*WS*PTYPE

      RETURN
      END SUBROUTINE PBL

      subroutine init_pbl(inipbl)
c -------------------------------------------------------------
c These routines include the array ipbl which indicates if the
c  computation for a particular ITYPE was done last time step.
c Sets up the initialization of wind, temperature, and moisture
c  fields in the boundary layer. The initial values of these
c  fields are obtained by solving the static equations of the
c  Level 2 model. This is used when starting from a restart
c  file that does not have this data stored.
c -------------------------------------------------------------
      USE CONSTANT, only : lhe,lhs,tf
      USE E001M12_COM
      USE PBLCOM, ustar_type=>ustar
      USE SOCPBL, only : npbl=>n,zgs,bgrid,inits,ccoeff0
     & ,  uinit=>u,vinit=>v,tinit=>t,qinit=>q,einit=>e
     &     ,dpdxrij=>dpdxr,dpdyrij=>dpdyr
     &     ,dpdxr0ij=>dpdxr0,dpdyr0ij=>dpdyr0
      USE DYNAMICS, only : pmid,pk,pedn,pek
      USE OCEAN, only : tocean
      USE SEAICE_COM, only : rsi
      USE GHYCOM, only : tearth
      USE LANDICE_COM, only : tlandi
      USE FLUXES, only : gtemp
      USE FILEMANAGER

      IMPLICIT NONE

!@var inipbl whether to init prog vars
      logical, intent(in) :: inipbl
!@var iu_CDN unit number for roughness length input file
      integer :: iu_CDN

      real*8, parameter :: ohmega=7.292e-5,rvx=0.

      integer :: ilong  !@var ilong  longitude identifier
      integer :: jlat   !@var jlat  latitude identifier
      real*8 tgvdat(im,jm,4)

      integer :: itype  !@var itype surface type
      integer i,j,iter,imax,im1,jvpo,lpbl !@var i,j,iter loop variable
      real*8 pland,pwater,plice,psoil,poice,pocean,pi,radian,
     *     ztop,elhx,coriol,tgrnd,pij,ps,psk,qgrnd,hemi
     *     ,utop,vtop,qtop,ttop,zgrnd,cm,ch,cq,ustar
      real*8 qsat

C things to be done regardless of inipbl
      call getunit("CDN",iu_CDN,.TRUE.,.true.)
      call readt (iu_CDN,0,roughl,im*jm,roughl,1)
      close (iu_CDN)

      call ccoeff0
      call getb(zgs,ztop,bgrid)

      if(.not.inipbl) return

      call pgrads1   !   added 6/19/00
      do j=1,jm
        do i=1,im
          pland=fland(i,j)
          pwater=1.-pland
          plice=flice(i,j)
          psoil=fearth(i,j)
          poice=rsi(i,j)*pwater
          pocean=pwater-poice
          tgvdat(i,j,1)=tocean(1,i,j) + TF
          if (pocean.le.0.) tgvdat(i,j,1)=0.
          tgvdat(i,j,2)=gtemp(1,2,i,j)+ TF
          if (poice.le.0.)  tgvdat(i,j,2)=0.
          tgvdat(i,j,3)=tlandi(1,i,j) + TF
          if (plice.le.0.)  tgvdat(i,j,3)=0.
          tgvdat(i,j,4)=tearth(i,j)   + TF
          if (psoil.le.0.)  tgvdat(i,j,4)=0.
        end do
      end do

      pi=dacos(-1.d0)
      radian=pi/180.

      do itype=1,4
        if ((itype.eq.1).or.(itype.eq.4)) then
          elhx=lhe
          else
          elhx=lhs
        endif
        do j=1,jm
          if ((j.eq.1).or.(j.eq.jm)) then
            imax=1
            else
            imax=im
          endif
          jlat=j
       coriol=2.*sin(radian*(float(j-1)*180./float(jm-1)-90.))*ohmega

          im1=im
          do i=1,imax
            tgrnd=tgvdat(i,j,itype)
            if (tgrnd.eq.0.) then
              ipbl(i,j,itype)=0
              im1=i
              go to 200
            endif
            ilong=i
            pij=p(i,j)
            ps=pedn(1,i,j)    !pij+ptop
            psk=pek(1,i,j)    !expbyk(ps)
            qgrnd=qsat(tgrnd,elhx,ps)

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

            dpdxrij  = dpdxr(i,j)
            dpdyrij  = dpdyr(i,j)
            dpdxr0ij = dpdxr0(i,j)
            dpdyr0ij = dpdyr0(i,j)

            call inits(tgrnd,qgrnd,zgrnd,zgs,ztop,utop,vtop,
     2                 ttop,qtop,coriol,cm,ch,cq,bgrid,ustar,
     3                 ilong,jlat,itype)
            cmgs(i,j,itype)=cm
            chgs(i,j,itype)=ch
            cqgs(i,j,itype)=cq

            do lpbl=1,npbl
              uabl(lpbl,i,j,itype)=uinit(lpbl)
              vabl(lpbl,i,j,itype)=vinit(lpbl)
              tabl(lpbl,i,j,itype)=tinit(lpbl)
              qabl(lpbl,i,j,itype)=qinit(lpbl)
            end do

            do lpbl=1,npbl-1
              eabl(lpbl,i,j,itype)=einit(lpbl)
            end do

            ipbl(i,j,itype)=1
            ustar_type(i,j,itype)=ustar
            im1=i

 200      end do
        end do
c     write (99,1000) itype
      end do

      return
 1000 format (1x,//,1x,'completed initialization, itype = ',i2,//)
      end subroutine init_pbl

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
      USE PBLCOM, only : npbl=>n,uabl,vabl,tabl,qabl,eabl,cmgs,chgs,cqgs
     *     ,ipbl,ustar_type=>ustar
      IMPLICIT NONE
      integer i,j,iter,lpbl,imax  !@var i,j,iter,lpbl loop variable

      do j=1,jm
        if ((j.eq.1).or.(j.eq.jm)) then
          imax=1
          else
          imax=im
        endif
        do i=1,imax

c ******* itype=1: Ocean

          if (ipbl(i,j,1).eq.0) then
            if (ipbl(i,j,2).eq.1) then
              do lpbl=1,npbl-1
                 uabl(lpbl,i,j,1)=uabl(lpbl,i,j,2)
                 vabl(lpbl,i,j,1)=vabl(lpbl,i,j,2)
                 tabl(lpbl,i,j,1)=tabl(lpbl,i,j,2)
                 qabl(lpbl,i,j,1)=qabl(lpbl,i,j,2)
                 eabl(lpbl,i,j,1)=eabl(lpbl,i,j,2)
              end do
              uabl(npbl,i,j,1)=uabl(npbl,i,j,2)
              vabl(npbl,i,j,1)=vabl(npbl,i,j,2)
              tabl(npbl,i,j,1)=tabl(npbl,i,j,2)
              qabl(npbl,i,j,1)=qabl(npbl,i,j,2)
              cmgs(i,j,1)=cmgs(i,j,2)
              chgs(i,j,1)=chgs(i,j,2)
              cqgs(i,j,1)=cqgs(i,j,2)
              ustar_type(i,j,1)=ustar_type(i,j,2)
            endif
          endif

c ******* itype=2: Ocean ice

          if (ipbl(i,j,2).eq.0) then
            if (ipbl(i,j,1).eq.1) then
              do lpbl=1,npbl-1
                 uabl(lpbl,i,j,2)=uabl(lpbl,i,j,1)
                 vabl(lpbl,i,j,2)=vabl(lpbl,i,j,1)
                 tabl(lpbl,i,j,2)=tabl(lpbl,i,j,1)
                 qabl(lpbl,i,j,2)=qabl(lpbl,i,j,1)
                 eabl(lpbl,i,j,2)=eabl(lpbl,i,j,1)
              end do
              uabl(npbl,i,j,2)=uabl(npbl,i,j,1)
              vabl(npbl,i,j,2)=vabl(npbl,i,j,1)
              tabl(npbl,i,j,2)=tabl(npbl,i,j,1)
              qabl(npbl,i,j,2)=qabl(npbl,i,j,1)
              cmgs(i,j,2)=cmgs(i,j,1)
              chgs(i,j,2)=chgs(i,j,1)
              cqgs(i,j,2)=cqgs(i,j,1)
              ustar_type(i,j,2)=ustar_type(i,j,1)
            endif
          endif

c ******* itype=3: Land ice

          if (ipbl(i,j,3).eq.0) then
            if (ipbl(i,j,4).eq.1) then
              do lpbl=1,npbl-1
                 uabl(lpbl,i,j,3)=uabl(lpbl,i,j,4)
                 vabl(lpbl,i,j,3)=vabl(lpbl,i,j,4)
                 tabl(lpbl,i,j,3)=tabl(lpbl,i,j,4)
                 qabl(lpbl,i,j,3)=qabl(lpbl,i,j,4)
                 eabl(lpbl,i,j,3)=eabl(lpbl,i,j,4)
              end do
              uabl(npbl,i,j,3)=uabl(npbl,i,j,4)
              vabl(npbl,i,j,3)=vabl(npbl,i,j,4)
              tabl(npbl,i,j,3)=tabl(npbl,i,j,4)
              qabl(npbl,i,j,3)=qabl(npbl,i,j,4)
              cmgs(i,j,3)=cmgs(i,j,4)
              chgs(i,j,3)=chgs(i,j,4)
              cqgs(i,j,3)=cqgs(i,j,4)
              ustar_type(i,j,3)=ustar_type(i,j,4)
            endif
          endif

c ******* itype=4: Land

          if (ipbl(i,j,4).eq.0) then
            if (ipbl(i,j,3).eq.1) then
              do lpbl=1,npbl-1
                 uabl(lpbl,i,j,4)=uabl(lpbl,i,j,3)
                 vabl(lpbl,i,j,4)=vabl(lpbl,i,j,3)
                 tabl(lpbl,i,j,4)=tabl(lpbl,i,j,3)
                 qabl(lpbl,i,j,4)=qabl(lpbl,i,j,3)
                 eabl(lpbl,i,j,4)=eabl(lpbl,i,j,3)
              end do
              uabl(npbl,i,j,4)=uabl(npbl,i,j,3)
              vabl(npbl,i,j,4)=vabl(npbl,i,j,3)
              tabl(npbl,i,j,4)=tabl(npbl,i,j,3)
              qabl(npbl,i,j,4)=qabl(npbl,i,j,3)
              cmgs(i,j,4)=cmgs(i,j,3)
              chgs(i,j,4)=chgs(i,j,3)
              cqgs(i,j,4)=cqgs(i,j,3)
              ustar_type(i,j,4)=ustar_type(i,j,3)
            endif
          endif

       end do
      end do

      return
      end subroutine loadbl

      subroutine pgrads1
      USE CONSTANT, only : rgas
      USE E001M12_COM, only : im,jm,t,p,sig,zatmo
      USE GEOM, only : dyp,dxp
      USE PBLCOM, only : dpdxr,dpdyr,phi,dpdxr0,dpdyr0,iq1,iq2,iq3
      USE DYNAMICS, only : pmid,pk,pedn
      IMPLICIT NONE

      integer i,j,iter  !@var i,j,iter loop variable
      real*8 p1k,t1,pij,rho1,dpx,dpy,dypsp,dypnp,rhojm,p1
      integer index1,index2
c      real*8 expbyk

c     for gcm main level 1:
      call geopot

      do j=2,jm-1
        do i=1,im

          pij=p(i,j)
          p1=100.*pmid(1,i,j)      !(pij*sig(1)+ptop)
          t1=t(i,j,1)*pk(1,i,j)    !expbyk(pij*sig(1)+ptop)
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
     2               (ZATMO(index1,j)-ZATMO(index2,j))/dxp(j)

          dpy=100.*(p(i,j+1)-p(i,j-1))

          dpdyr0(i,j)=0.5*dpy/(dyp(j)*rho1)+0.5*
     2               (ZATMO(i,j+1)-ZATMO(i,j-1))/dyp(j)

       end do
      end do

      i=1
      j=1
      pij=p(i,j)
      p1=100.*pmid(1,i,j)            !(pij*sig(1)+ptop)
      t1=t(i,j,1)*pk(1,i,j)          !expbyk(pij*sig(1)+ptop)
      rho1 =p1/(rgas*t1)
      dypsp=2.*dyp(1)

      j=jm
      pij=p(i,j)
      p1=100.*pmid(1,i,j)            !(pij*sig(1)+ptop)
      t1=t(i,j,1)*pk(1,i,j)          !expbyk(pij*sig(1)+ptop)
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
     4            0.125*(ZATMO(iq1  ,   2)-ZATMO(iq3  ,   2)+
     5                   ZATMO(iq1+1,   2)-ZATMO(iq3+1,   2))/dypsp
      dpdyr0(1, 1)=0.25*(p(    1,   2)-p(iq2  ,   2)+
     2                  p(    2,   2)-p(iq2+1,   2))*100./
     3                 (dypsp*rho1)+
     4            0.125*(ZATMO(    1,   2)-ZATMO(iq2  ,   2)+
     5                   ZATMO(    2,   2)-ZATMO(iq2+1,   2))/dypsp



      dpdxr0(1,jm)=0.25*(p(iq1  ,jm-1)-p(iq3  ,jm-1)+
     2                  p(iq1+1,jm-1)-p(iq3+1,jm-1))*100./
     3                 (dypnp*rhojm)+
     4            0.125*(ZATMO(iq1  ,jm-1)-ZATMO(iq3  ,jm-1)+
     5                   ZATMO(iq1+1,jm-1)-ZATMO(iq3+1,jm-1))/dypnp
      dpdyr0(1,jm)=0.25*(p(iq2  ,jm-1)-p(    1,jm-1)+
     2                  p(iq2+1,jm-1)-p(    2,jm-1))*100./
     3                 (dypnp*rhojm)+
     4            0.125*(ZATMO(iq2  ,jm-1)-ZATMO(    1,jm-1)+
     5                   ZATMO(iq2+1,jm-1)-ZATMO(    2,jm-1))/dypnp

      return
      end subroutine pgrads1

      subroutine geopot
      USE CONSTANT, only : rgas,grav
      USE E001M12_COM, only : im,jm,t,p,dsig,zatmo
      USE PBLCOM, only : phi
      USE SOCPBL, only : zgs
      USE DYNAMICS, only : pmid,pk,pedn
      IMPLICIT NONE

c      real*8 expbyk

      integer i,j,iter  !@var i,j,iter loop variable
      real*8 p1,p1k,pij,t1,z1

c     note: ZATMO(I,J) is the geopotential height (9.81*zatm)
c
c     for GCM main level 1:
      do j=2,jm-1
        do i=1,im
          pij=p(i,j)
          p1=pmid(1,i,j)     !pij*sig(1)+ptop
          p1k=pk(1,i,j)      !expbyk(p1)
          t1=t(i,j,1)*p1k
          z1=zgs+0.5*dsig(1)*rgas*t1*pij/(p1*grav)
          phi(i,j)=grav*z1+ZATMO(I,J)
        end do
      end do
      i=1
      j=1
      pij=p(i,j)
      p1=pmid(1,i,j)          !pij*sig(1)+ptop
      p1k=pk(1,i,j)           !expbyk(p1)
      t1=t(i,j,1)*p1k
      z1=zgs+0.5*dsig(1)*rgas*t1*pij/(p1*grav)
      phi(i,j)=grav*z1+ZATMO(I,J)

      j=jm
      pij=p(i,j)
      p1=pmid(1,i,j)          !pij*sig(1)+ptop
      p1k=pk(1,i,j)           !expbyk(p1)
      t1=t(i,j,1)*p1k
      z1=zgs+0.5*dsig(1)*rgas*t1*pij/(p1*grav)
      phi(i,j)=grav*z1+ZATMO(I,J)
      return
      end subroutine geopot

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
      USE CONSTANT, only : rgas,grav
      USE E001M12_COM, only : sige,psf,ptop,psfmpt
      IMPLICIT NONE

      REAL*8, INTENT(IN) :: ZGS
      REAL*8, INTENT(OUT) :: ZTOP,BGRID
      real*8 theta,z1,x

      theta=269.0727251
      z1=zgs+0.5*(1.-sige(2))*psfmpt*rgas*theta/(grav*psf)
      x=z1/100.
      ztop=z1
      bgrid=0.177427*x**4 - 1.0504*x**3 + 2.34169*x**2 -
     2      2.4772*x + 1.44509
      return
      end subroutine getb

      SUBROUTINE CHECKPBL(SUBR)
!@sum  CHECKPBL Checks whether PBL data are reasonable
!@auth Original Development Team
!@ver  1.0
      USE E001M12_COM, only : im,jm
      USE PBLCOM, only : wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg
     *     ,ustar
      IMPLICIT NONE

!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

C**** Check for NaN/INF in boundary layer data
      CALL CHECK3(wsavg,IM,JM,1,SUBR,'wsavg')
      CALL CHECK3(tsavg,IM,JM,1,SUBR,'tsavg')
      CALL CHECK3(qsavg,IM,JM,1,SUBR,'qsavg')
      CALL CHECK3(dclev,IM,JM,1,SUBR,'dclev')
      CALL CHECK3(usavg,IM,JM,1,SUBR,'usavg')
      CALL CHECK3(vsavg,IM,JM,1,SUBR,'vsavg')
      CALL CHECK3(tauavg,IM,JM,1,SUBR,'tauavg')
      CALL CHECK3(ustar,IM,JM,4,SUBR,'ustar')

      END SUBROUTINE CHECKPBL

