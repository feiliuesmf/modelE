#include "rundeck_opts.h"

      module PBL_DRV
#ifdef TRACERS_ON
      use tracer_com, only : ntm
#endif
      implicit none
ccc   save

c     input data:
!@var uocean,vocean ocean/ice velocities for use in drag calulation
      real*8 uocean,vocean

!@var evap_max maximal evaporation from unsaturated soil
!@var  fr_sat fraction of saturated soil
      real*8 :: evap_max,fr_sat
      real*8 :: psurf, trhr0
      common/pbl_loc/evap_max,fr_sat,uocean,vocean,psurf,trhr0

!$OMP  THREADPRIVATE (/pbl_loc/)



#ifdef TRACERS_ON
C**** Tracer input/output common block for PBL
!@var trtop,trs tracer mass ratio in level 1/surface
!@var trsfac, trconstflx factors in surface flux boundary cond.
!@var ntx number of tracers that need pbl calculation
!@var ntix index array to map local tracer number to global
      real*8, dimension(ntm) :: trtop,trs,trsfac,trconstflx
      integer ntx
      integer, dimension(ntm) :: ntix
#ifdef TRACERS_WATER
!@var tr_evap_max maximum amount of tracer available in ground reservoir
      real*8, dimension(ntm) :: tr_evap_max
#endif
#ifdef TRACERS_DRYDEP
!@var dep_vel turbulent deposition velocity = 1/bulk sfc. res. (m/s)
!@var gs_vel gravitational settling velocity (m/s)
      real*8, dimension(ntm) :: dep_vel,gs_vel
#endif
#ifdef TRACERS_AEROSOLS_Koch
      real*8 :: DMS_flux,ss1_flux,ss2_flux
#endif
      common /trspec/trtop,trs,trsfac,trconstflx
#ifdef TRACERS_WATER
     *     ,tr_evap_max
#endif
#ifdef TRACERS_DRYDEP
     *     ,dep_vel,gs_vel
#endif
#ifdef TRACERS_AEROSOLS_Koch
     *     ,DMS_flux,ss1_flux,ss2_flux
#endif

!$OMP  THREADPRIVATE (/trspec/)
#endif

      contains

      SUBROUTINE PBL(I,J,ITYPE,PTYPE)
!@sum  PBL calculate pbl profiles for each surface type
!@auth Greg. Hartke/Ye Cheng
!@ver  1.0

C    input: ZS1,TGV,TKV,QG_SAT,qg_aver,HEMI,DTSURF,POLE,UOCEAN,VOCEAN
C    output:US,VS,WS,WSM,WSH,TSV,QSRF,PSI,DBL,KMS,KHS,KQS,PPBL
C          ,UG,VG,WG,W2_1

      USE CONSTANT, only :  rgas,grav,omega2,deltx,teeny
      USE MODEL_COM
     &     , only : IM,JM,LM, t,q,u,v,p,ptop,ls1,psf,itime
      USE GEOM, only : idij,idjj,kmaxj,rapj,cosiv,siniv,sinp
      USE DYNAMICS, only : pmid,pk,pedn,pek
     &    ,DPDX_BY_RHO,DPDY_BY_RHO,DPDX_BY_RHO_0,DPDY_BY_RHO_0
      USE CLOUDS_COM, only : ddm1
      USE SOCPBL, only : npbl=>n
     &     ,dpdxr,dpdyr
     &     ,dpdxr0,dpdyr0
     &     ,advanc                      ! subroutine
     &     ,zgs,DTSURF                  ! global
     &     ,ZS1,TGV,TKV,QG_SAT,qg_aver,HEMI,POLE    ! rest local
     &     ,US,VS,WS,WSM,WSH,TSV,QSRF,PSI,DBL,KMS,KHS,KQS,PPBL
     &     ,UG,VG,WG,mdf
     &     ,ustar,cm,ch,cq,z0m,z0h,z0q,w2_1
#ifdef TRACERS_ON
     *     ,tr
#ifdef TRACERS_AEROSOLS_Koch
      USE AEROSOL_SOURCES, only: PBLH
#endif
#endif

      USE PBLCOM
      use QUSDEF, only : mz
      use SOMTQ_COM, only : tmom
 
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: I,J  !@var I,J grid point
      INTEGER, INTENT(IN) :: ITYPE  !@var ITYPE surface type
      REAL*8, INTENT(IN) :: PTYPE  !@var PTYPE percent surface type

      REAL*8, parameter :: dbl_max=3000., dbl_max_stable=500. ! meters
      real*8, parameter :: S1byG1=.57735d0
      REAL*8 Ts,ts_guess

#ifdef TRACERS_ON
      integer nx
#endif
c
      REAL*8 ztop,zpbl,pl1,tl1,pl,tl,tbar,thbar,zpbl1,coriol
      REAL*8 ttop,qtop,tgrndv,qgrnd,qgrnd_sat,utop,vtop,ufluxs,vfluxs
     *     ,tfluxs,qfluxs,psitop,psisrf
      INTEGER LDC,L,k

c**** special threadprivate common block (compaq compiler stupidity)
      real*8, dimension(npbl) :: upbl,vpbl,tpbl,qpbl
      real*8, dimension(npbl-1) :: epbl
      common/pbluvtq/upbl,vpbl,tpbl,qpbl,epbl
!$OMP  THREADPRIVATE (/pbluvtq/)
C**** end special threadprivate common block

C        ocean and ocean ice are treated as rough surfaces
C        roughness lengths from Brutsaert for rough surfaces

      IF (ITYPE.GT.2) THEN
        Z0M=30./(10.**ROUGHL(I,J))
      ENDIF
      ztop=zgs+zs1  ! zs1 is calculated before pbl is called
      IF (TKV.EQ.TGV) TGV=1.0001d0*TGV

      ! FIND THE PBL HEIGHT IN METERS (DBL) AND THE CORRESPONDING
      ! GCM LAYER (L) AT WHICH TO COMPUTE UG AND VG.
      ! LDC IS THE LAYER TO WHICH DRY CONVECTION/TURBULENCE MIXES

        IF (TKV.GE.TGV) THEN
          ! ATMOSPHERE IS STABLE WITH RESPECT TO THE GROUND
          ! DETERMINE VERTICAL LEVEL CORRESPONDING TO HEIGHT OF PBL:
          ! WHEN ATMOSPHERE IS STABLE, CAN COMPUTE DBL BUT DO NOT
          ! KNOW THE INDEX OF THE LAYER.
          ustar=ustar_pbl(i,j,itype)
          DBL=min(0.3d0*USTAR/OMEGA2,dbl_max_stable)
          if (dbl.le.ztop) then
            dbl=ztop
            L=1
          else
            ! FIND THE VERTICAL LEVEL NEXT HIGHER THAN DBL AND
            ! COMPUTE Ug and Vg THERE:
            zpbl=ztop
            pl1=pmid(1,i,j)         ! pij*sig(1)+ptop
            tl1=t(i,j,1)*(1.+deltx*q(i,j,1))*pk(1,i,j)
            do l=2,ls1
              pl=pmid(l,i,j)        !pij*sig(l)+ptop
              tl=t(i,j,l)*(1.+deltx*q(i,j,l))*pk(l,i,j) !virtual,absolute
              tbar=thbar(tl1,tl)
              zpbl=zpbl-(rgas/grav)*tbar*(pl-pl1)/(pl1+pl)*2.
              if (zpbl.ge.dbl) exit
              pl1=pl
              tl1=tl
            end do
          endif

      ELSE
        ! ATMOSPHERE IS UNSTABLE WITH RESPECT TO THE GROUND
        ! LDC IS THE LEVEL TO WHICH DRYCNV/ATURB MIXES.
        ! FIND DBL FROM LDC.  IF BOUNDARY
        ! LAYER HEIGHT IS LESS THAN DBL_MAX, ASSIGN LDC TO L, OTHERWISE
        ! MUST FIND INDEX FOR NEXT MODEL LAYER ABOVE 3 KM:

        LDC=max(int(DCLEV(I,J)+.5d0),1)
        IF (LDC.EQ.0) LDC=1
        if (ldc.eq.1) then
          dbl=ztop
          l=1
        else
          zpbl=ztop
          pl1=pmid(1,i,j)                             ! pij*sig(1)+ptop
          tl1=t(i,j,1)*(1.+deltx*q(i,j,1))*pk(1,i,j)  ! expbyk(pl1)
          zpbl1=ztop
          do l=2,ldc
            pl=pmid(l,i,j)                            ! pij*sig(l)+ptop
            tl=t(i,j,l)*(1.+deltx*q(i,j,l))*pk(l,i,j) ! expbyk(pl)
            tbar=thbar(tl1,tl)
            zpbl=zpbl-(rgas/grav)*tbar*(pl-pl1)/(pl1+pl)*2.
            if (zpbl.ge.dbl_max) then
              zpbl=zpbl1
              exit
            endif
            pl1=pl
            tl1=tl
            zpbl1=zpbl
          end do
          l=min(l,ldc)
          dbl=zpbl
        endif

      ENDIF

#if (defined TRACERS_DRYDEP) && (defined TRACERS_AEROSOLS_Koch)
      IF (ITYPE.eq.1) PBLH(I,J)=dbl
#endif
      ppbl=pedn(l,i,j)
      coriol=sinp(j)*omega2
      ttop=tkv
      qtop=q(i,j,1)
      tgrndv=tgv
      qgrnd_sat=qg_sat
      qgrnd=qg_aver

      utop=0. ; vtop=0. ;  ug=0. ; vg=0.
      ! pole and hemi are determined before pbl is called
      if (pole) then
        do k=1,kmaxj(j)
          utop = utop + rapj(k,j)*(u(idij(k,i,j),idjj(k,j),1)*cosiv(k) -
     2                        hemi*v(idij(k,i,j),idjj(k,j),1)*siniv(k))
          vtop = vtop + rapj(k,j)*(v(idij(k,i,j),idjj(k,j),1)*cosiv(k) +
     2                        hemi*u(idij(k,i,j),idjj(k,j),1)*siniv(k))
          ug   = ug   + rapj(k,j)*(u(idij(k,i,j),idjj(k,j),L)*cosiv(k) -
     2                        hemi*v(idij(k,i,j),idjj(k,j),L)*siniv(k))
          vg   = vg   + rapj(k,j)*(v(idij(k,i,j),idjj(k,j),L)*cosiv(k) +
     2                        hemi*u(idij(k,i,j),idjj(k,j),L)*siniv(k))
        end do
      else
        do k=1,kmaxj(j)
          utop = utop + u(idij(k,i,j),idjj(k,j),1)*rapj(k,j)
          vtop = vtop + v(idij(k,i,j),idjj(k,j),1)*rapj(k,j)
          ug   = ug   + u(idij(k,i,j),idjj(k,j),L)*rapj(k,j)
          vg   = vg   + v(idij(k,i,j),idjj(k,j),L)*rapj(k,j)
        end do
      endif

      upbl(:)=uabl(:,i,j,itype)
      vpbl(:)=vabl(:,i,j,itype)
      tpbl(:)=tabl(:,i,j,itype)
      qpbl(:)=qabl(:,i,j,itype)
      epbl(1:npbl-1)=eabl(1:npbl-1,i,j,itype)
#ifdef TRACERS_ON
      do nx=1,ntx
        tr(:,nx)=trabl(:,ntix(nx),i,j,itype)
      end do
#endif

      cm=cmgs(i,j,itype)
      ch=chgs(i,j,itype)
      cq=cqgs(i,j,itype)
      dpdxr  = DPDX_BY_RHO(i,j)
      dpdyr  = DPDY_BY_RHO(i,j)
      dpdxr0 = DPDX_BY_RHO_0(i,j)
      dpdyr0 = DPDY_BY_RHO_0(i,j)

      ts_guess = (t(i,j,1)-tmom(mz,i,j,1)*S1byG1)*pek(1,i,j)
     2          *(1+q(i,j,1)*deltx)
      mdf = ddm1(i,j)

      call advanc(
     3     coriol,utop,vtop,ttop,qtop,tgrndv,
     &     qgrnd,qgrnd_sat,evap_max,fr_sat,
#if defined(TRACERS_ON)
     *     trs,trtop,trsfac,trconstflx,ntx,ntix,
#if defined(TRACERS_WATER)
     *     tr_evap_max,
#endif
#if defined(TRACERS_DRYDEP)
     *     dep_vel,gs_vel,
#endif
#ifdef TRACERS_AEROSOLS_Koch
     *     DMS_flux,ss1_flux,ss2_flux,
#endif
#endif
     4     psurf,trhr0,ztop,dtsurf,ufluxs,vfluxs,tfluxs,qfluxs,
     5     uocean,vocean,ts_guess,i,j,itype)

      uabl(:,i,j,itype)=upbl(:)
      vabl(:,i,j,itype)=vpbl(:)
      tabl(:,i,j,itype)=tpbl(:)
      qabl(:,i,j,itype)=qpbl(:)
      eabl(1:npbl-1,i,j,itype)=epbl(1:npbl-1)
#ifdef TRACERS_ON
      do nx=1,ntx
        trabl(:,ntix(nx),i,j,itype)=tr(:,nx)
      end do
#endif

      cmgs(i,j,itype)=cm
      chgs(i,j,itype)=ch
      cqgs(i,j,itype)=cq
      ipbl(i,j,itype)=1  ! ipbl is used in subroutine init_pbl

      ws    =wsm
      wg    =sqrt(ug*ug+vg*vg)

      psitop=atan2(vg,ug+teeny)
      psisrf=atan2(vs,us+teeny)
      psi   =psisrf-psitop
      ustar_pbl(i,j,itype)=ustar
C ******************************************************************
      TS=TSV/(1.+QSRF*deltx)
      if ( ts.lt.152d0 .or. ts.gt.423d0 ) then
        write(6,*) 'PBL: Ts bad at',i,j,' itype',itype,ts
        if (ts.gt.1d3) call stop_model("PBL: Ts out of range",255)
        if (ts.lt.50d0) call stop_model("PBL: Ts out of range",255)
      end if
      WSAVG(I,J)=WSAVG(I,J)+WS*PTYPE
      TSAVG(I,J)=TSAVG(I,J)+TS*PTYPE
      if(itype.ne.4) QSAVG(I,J)=QSAVG(I,J)+QSRF*PTYPE
      USAVG(I,J)=USAVG(I,J)+US*PTYPE
      VSAVG(I,J)=VSAVG(I,J)+VS*PTYPE
      TAUAVG(I,J)=TAUAVG(I,J)+CM*WS*WS*PTYPE

      uflux(I,J)=uflux(I,J)+ufluxs*PTYPE
      vflux(I,J)=vflux(I,J)+vfluxs*PTYPE
      tflux(I,J)=tflux(I,J)+tfluxs*PTYPE
      qflux(I,J)=qflux(I,J)+qfluxs*PTYPE

      tgvAVG(I,J)=tgvAVG(I,J)+tgv*PTYPE
      qgAVG(I,J)=qgAVG(I,J)+qgrnd*PTYPE
      w2_l1(I,J)=w2_l1(I,J)+w2_1*PTYPE

      RETURN
      END SUBROUTINE PBL

      end module PBL_DRV

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
      USE FILEMANAGER
      USE PARAM
      USE CONSTANT, only : lhe,lhs,tf,omega2,deltx
      USE MODEL_COM
      USE GEOM, only : idij,idjj,imaxj,kmaxj,rapj,cosiv,siniv,sinp
      USE SOCPBL, only : npbl=>n,zgs,inits,ccoeff0,XCDpbl
     &     ,dpdxr,dpdyr,dpdxr0,dpdyr0
      USE PBLCOM
      USE DOMAIN_DECOMP, only : GRID, GET, READT_PARALLEL
      USE DOMAIN_DECOMP, only : HALO_UPDATE,CHECKSUM,NORTH
      USE DYNAMICS, only : pmid,pk,pedn,pek
     &    ,DPDX_BY_RHO,DPDY_BY_RHO,DPDX_BY_RHO_0,DPDY_BY_RHO_0
      USE SEAICE_COM, only : rsi,snowi
      USE FLUXES, only : gtemp

      IMPLICIT NONE

C**** ignore ocean currents for initialisation.
      real*8, parameter :: uocean=0.,vocean=0.
!@var inipbl whether to init prog vars
      logical, intent(in) :: inipbl
!@var iu_CDN unit number for roughness length input file
      integer :: iu_CDN
      integer :: ilong  !@var ilong  longitude identifier
      integer :: jlat   !@var jlat  latitude identifier
      real*8, dimension(im,GRID%J_STRT_HALO:GRID%J_STOP_HALO,4) ::
     *                                                      tgvdat

      integer :: itype  !@var itype surface type
      integer i,j,k,iter,lpbl !@var i,j,k,iter loop variable
      real*8 pland,pwater,plice,psoil,poice,pocean,
     *     ztop,elhx,coriol,tgrndv,pij,ps,psk,qgrnd
     *     ,utop,vtop,qtop,ttop,zgrnd,cm,ch,cq,ustar
      real*8 qsat

c**** special threadprivate common block (compaq compiler stupidity)
      real*8, dimension(npbl) :: upbl,vpbl,tpbl,qpbl
      real*8, dimension(npbl-1) :: epbl
      common/pbluvtq/upbl,vpbl,tpbl,qpbl,epbl
!$OMP  THREADPRIVATE (/pbluvtq/)
C**** end special threadprivate common block

      integer :: J_1, J_0
      integer :: J_1H, J_0H

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     *               J_STRT=J_0,       J_STOP=J_1)


C things to be done regardless of inipbl
      call openunit("CDN",iu_CDN,.TRUE.,.true.)
      CALL READT_PARALLEL(grid,iu_CDN,NAMEUNIT(iu_CDN),0,roughl,1)
      call closeunit(iu_CDN)
      call sync_param( 'XCDpbl', XCDpbl )

      do j=J_0,J_1
        do i=1,im
C**** fix roughness length for ocean ice that turned to land ice
          if (snowi(i,j).lt.-1.and.flice(i,j).gt.0) roughl(i,j)=1.84d0
          if (fland(i,j).gt.0.and.roughl(i,j).eq.0) then
            print*,"Roughness length not defined for i,j",i,j
     *           ,roughl(i,j),fland(i,j),flice(i,j)
            roughl(i,j)=roughl(10,40)
          end if
        end do
      end do

      call ccoeff0
      call getztop(zgs,ztop)

      if(.not.inipbl) return

      do j=J_0,J_1
      do i=1,im
        pland=fland(i,j)
        pwater=1.-pland
        plice=flice(i,j)
        psoil=fearth(i,j)
        poice=rsi(i,j)*pwater
        pocean=pwater-poice
        if (pocean.le.0.) then
          tgvdat(i,j,1)=0.
        else
          tgvdat(i,j,1)=gtemp(1,1,i,j)+TF
        end if
        if (poice.le.0.) then
          tgvdat(i,j,2)=0.
        else
          tgvdat(i,j,2)=gtemp(1,2,i,j)+TF
        end if
        if (plice.le.0.) then
          tgvdat(i,j,3)=0.
        else
          tgvdat(i,j,3)=gtemp(1,3,i,j)+TF
        end if
        if (psoil.le.0.) then
          tgvdat(i,j,4)=0.
        else
          tgvdat(i,j,4)=gtemp(1,4,i,j)+TF
        end if
      end do
      end do

      do itype=1,4
        if ((itype.eq.1).or.(itype.eq.4)) then
          elhx=lhe
        else
          elhx=lhs
        endif
C**** HALO UPDATES OF u AND v FOR DISTRIBUTED PARALLELIZATION 
        call CHECKSUM   (grid, u, __LINE__, __FILE__)
        call HALO_UPDATE(grid, u, from=NORTH)
        call CHECKSUM   (grid, v, __LINE__, __FILE__)
        call HALO_UPDATE(grid, v, from=NORTH)
        do j=J_0,J_1
          jlat=j
          coriol=sinp(j)*omega2

          do i=1,imaxj(j)
            tgrndv=tgvdat(i,j,itype)
            if (tgrndv.eq.0.) then
              ipbl(i,j,itype)=0
              go to 200
            endif
            ilong=i
            pij=p(i,j)
            ps=pedn(1,i,j)    !pij+ptop
            psk=pek(1,i,j)    !expbyk(ps)
            qgrnd=qsat(tgrndv,elhx,ps)

            utop = 0. ;  vtop = 0.
            if (j.eq.1) then
c ******************************************************************
c           At the south pole:
              do k=1,kmaxj(j)
                utop = utop + (u(idij(k,i,j),idjj(k,j),1)*cosiv(k) +
     2                    v(idij(k,i,j),idjj(k,j),1)*siniv(k))*rapj(k,j)
                vtop = vtop + (v(idij(k,i,j),idjj(k,j),1)*cosiv(k) -
     2                    u(idij(k,i,j),idjj(k,j),1)*siniv(k))*rapj(k,j)
              end do
c ******************************************************************

            else if (j.eq.jm) then
c ******************************************************************
c     At the north pole:
              do k=1,kmaxj(j)
                utop = utop + (u(idij(k,i,j),idjj(k,j),1)*cosiv(k) -
     2                    v(idij(k,i,j),idjj(k,j),1)*siniv(k))*rapj(k,j)
                vtop = vtop + (v(idij(k,i,j),idjj(k,j),1)*cosiv(k) +
     2                    u(idij(k,i,j),idjj(k,j),1)*siniv(k))*rapj(k,j)
              end do
c ******************************************************************

            else
c ******************************************************************
c     Away from the poles:
              do k=1,kmaxj(j)
                utop = utop + u(idij(k,i,j),idjj(k,j),1)*rapj(k,j)
                vtop = vtop + v(idij(k,i,j),idjj(k,j),1)*rapj(k,j)
              end do
c ******************************************************************
            endif

            qtop=q(i,j,1)
            ttop=t(i,j,1)*(1.+qtop*deltx)*psk

            zgrnd=.1d0 ! formal initialization
            if (itype.gt.2) zgrnd=30./(10.**roughl(i,j))

            dpdxr  = DPDX_BY_RHO(i,j)
            dpdyr  = DPDY_BY_RHO(i,j)
            dpdxr0 = DPDX_BY_RHO_0(i,j)
            dpdyr0 = DPDY_BY_RHO_0(i,j)

            call inits(tgrndv,qgrnd,zgrnd,zgs,ztop,utop,vtop,
     2                 ttop,qtop,coriol,cm,ch,cq,ustar,
     3                 uocean,vocean,ilong,jlat,itype)
            cmgs(i,j,itype)=cm
            chgs(i,j,itype)=ch
            cqgs(i,j,itype)=cq

            do lpbl=1,npbl
              uabl(lpbl,i,j,itype)=upbl(lpbl)
              vabl(lpbl,i,j,itype)=vpbl(lpbl)
              tabl(lpbl,i,j,itype)=tpbl(lpbl)
              qabl(lpbl,i,j,itype)=qpbl(lpbl)
            end do

            do lpbl=1,npbl-1
              eabl(lpbl,i,j,itype)=epbl(lpbl)
            end do

            ipbl(i,j,itype)=1
            ustar_pbl(i,j,itype)=ustar

 200      end do
        end do
      end do

      return
 1000 format (1x,//,1x,'completed initialization, itype = ',i2,//)
      end subroutine init_pbl

      subroutine loadbl
!@sum loadbl initiallise boundary layer calc each surface time step
!@auth Ye Cheng
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
      USE MODEL_COM
      USE GEOM, only : imaxj
      USE DOMAIN_DECOMP, only : GRID, GET
      USE PBLCOM, only : npbl,uabl,vabl,tabl,qabl,eabl,cmgs,chgs,cqgs
     *     ,ipbl,ustar_pbl,wsavg,tsavg,qsavg,usavg,vsavg,tauavg
     &     ,uflux,vflux,tflux,qflux,tgvavg,qgavg,w2_l1
#ifdef TRACERS_ON
     *     ,trabl
#endif
      IMPLICIT NONE
      integer i,j,iter,lpbl  !@var i,j,iter,lpbl loop variable

      integer :: J_1, J_0
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

      do j=J_0,J_1
      do i=1,imaxj(j)

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
#ifdef TRACERS_ON
              trabl(:,:,i,j,1)=trabl(:,:,i,j,2)
#endif
              cmgs(i,j,1)=cmgs(i,j,2)
              chgs(i,j,1)=chgs(i,j,2)
              cqgs(i,j,1)=cqgs(i,j,2)
              ustar_pbl(i,j,1)=ustar_pbl(i,j,2)
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
#ifdef TRACERS_ON
              trabl(:,:,i,j,2)=trabl(:,:,i,j,1)
#endif
              cmgs(i,j,2)=cmgs(i,j,1)
              chgs(i,j,2)=chgs(i,j,1)
              cqgs(i,j,2)=cqgs(i,j,1)
              ustar_pbl(i,j,2)=ustar_pbl(i,j,1)
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
#ifdef TRACERS_ON
              trabl(:,:,i,j,3)=trabl(:,:,i,j,4)
#endif
              cmgs(i,j,3)=cmgs(i,j,4)
              chgs(i,j,3)=chgs(i,j,4)
              cqgs(i,j,3)=cqgs(i,j,4)
              ustar_pbl(i,j,3)=ustar_pbl(i,j,4)
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
#ifdef TRACERS_ON
              trabl(:,:,i,j,4)=trabl(:,:,i,j,3)
#endif
              cmgs(i,j,4)=cmgs(i,j,3)
              chgs(i,j,4)=chgs(i,j,3)
              cqgs(i,j,4)=cqgs(i,j,3)
              ustar_pbl(i,j,4)=ustar_pbl(i,j,3)
            endif
          endif

C**** initialise some pbl common variables
          WSAVG(I,J)=0.
          TSAVG(I,J)=0.
          QSAVG(I,J)=0.
          USAVG(I,J)=0.
          VSAVG(I,J)=0.
          TAUAVG(I,J)=0.
          TGVAVG(I,J)=0.
          QGAVG(I,J)=0.
          w2_l1(I,J)=0.

          uflux(I,J)=0.
          vflux(I,J)=0.
          tflux(I,J)=0.
          qflux(I,J)=0.

      end do
      end do

      return
      end subroutine loadbl

      subroutine getztop(zgs,ztop)
!@sum  getztop computes the value of ztop which is the height in meters
!@+  of the first GCM layer from the surface.
!@+  This subroutine only needs to be called when the BL fields require
!@+  initialization.
!@+  This form for z1 = zgs + zs1 (in terms of GCM parameters) yields an
!@+  average value for zs1. The quantity theta was computed on the
!@+  assumption of zs1=200 m from the original 9-layer model (actually
!@+  was misconstrued as z1 = 200m when it should have been zs1 = 200m)
!@+  and is then applied to all vertical resolutions.
!@auth Greg. Hartke/Ye Cheng
!@var zgs The height of the surface layer.
!@var ztop The height of the top of the BL simulation domain.
!@+   Corresponds to averaged height of the middle of first model layer.

      USE CONSTANT, only : rgas,grav
      USE MODEL_COM, only : sige,psf,psfmpt
      IMPLICIT NONE

      REAL*8, INTENT(IN) :: ZGS
      REAL*8, INTENT(OUT) :: ZTOP
      real*8, parameter :: theta=269.0727251d0

      ztop=zgs+0.5d0*(1.-sige(2))*psfmpt*rgas*theta/(grav*psf)

      return
      end subroutine getztop

      SUBROUTINE CHECKPBL(SUBR)
!@sum  CHECKPBL Checks whether PBL data are reasonable
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : im,jm
      USE DOMAIN_DECOMP, only : GRID, GET
      USE PBLCOM, only : wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg
     *     ,ustar_pbl,uflux,vflux,tflux,qflux,tgvavg,qgavg,w2_l1
      IMPLICIT NONE

!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

      integer :: I_1, I_0, J_1, J_0, Ilen, Jlen
C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, I_STRT=I_0, I_STOP=I_1,
     *               J_STRT=J_0, J_STOP=J_1)

      Ilen = I_1-I_0+1
      Jlen = J_1-J_0+1

C**** Check for NaN/INF in boundary layer data
      CALL CHECK3(wsavg(I_0:I_1,J_0:J_1),Ilen,Jlen,1,SUBR,'wsavg')
      CALL CHECK3(tsavg(I_0:I_1,J_0:J_1),Ilen,Jlen,1,SUBR,'tsavg')
      CALL CHECK3(qsavg(I_0:I_1,J_0:J_1),Ilen,Jlen,1,SUBR,'qsavg')
      CALL CHECK3(dclev(I_0:I_1,J_0:J_1),Ilen,Jlen,1,SUBR,'dclev')
      CALL CHECK3(usavg(I_0:I_1,J_0:J_1),Ilen,Jlen,1,SUBR,'usavg')
      CALL CHECK3(vsavg(I_0:I_1,J_0:J_1),Ilen,Jlen,1,SUBR,'vsavg')
      CALL CHECK3(tauavg(I_0:I_1,J_0:J_1),Ilen,Jlen,1,SUBR,'tauavg')
      CALL CHECK3(ustar_pbl(I_0:I_1,J_0:J_1,1:4),Ilen,Jlen,4,SUBR,
     *           'ustar')

      CALL CHECK3(uflux(I_0:I_1,J_0:J_1),Ilen,Jlen,1,SUBR,'uflux')
      CALL CHECK3(vflux(I_0:I_1,J_0:J_1),Ilen,Jlen,1,SUBR,'vflux')
      CALL CHECK3(tflux(I_0:I_1,J_0:J_1),Ilen,Jlen,1,SUBR,'tflux')
      CALL CHECK3(qflux(I_0:I_1,J_0:J_1),Ilen,Jlen,1,SUBR,'qflux')

      CALL CHECK3(tgvavg(I_0:I_1,J_0:J_1),Ilen,Jlen,1,SUBR,'tgvavg')
      CALL CHECK3(qgavg(I_0:I_1,J_0:J_1),Ilen,Jlen,1,SUBR,'qgavg')
      CALL CHECK3(w2_l1(I_0:I_1,J_0:J_1),Ilen,Jlen,1,SUBR,'w2_l1')

      END SUBROUTINE CHECKPBL

