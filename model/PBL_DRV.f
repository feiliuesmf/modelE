#include "rundeck_opts.h"

      SUBROUTINE PBL(I,J,ITYPE,PTYPE)
!@sum  PBL calculate pbl profiles for each surface type
!@auth Greg. Hartke/Ye Cheng
!@ver  1.0

C    input: ZS1,TGV,TKV,QG,HEMI,DTSURF,POLE
C    output:US,VS,WS,WSH,WSQ,TSV,QS,PSI,DBL,KMS,KHS,KQS,PPBL
C          ,UG,VG,WG,ZMIX

      USE CONSTANT, only :  rgas,grav,omega2,deltx,teeny
      USE MODEL_COM
     &     , only : IM,JM,LM, t,q,u,v,p,ptop,ls1,psf,itime
      USE DYNAMICS, only : pmid,pk,pedn
     &    ,DPDX_BY_RHO,DPDY_BY_RHO,DPDX_BY_RHO_0,DPDY_BY_RHO_0
      USE GEOM, only : idij,idjj,kmaxj,rapj,cosiv,siniv,sinp
      USE PBLCOM, ustar_type=>ustar
      USE SOCPBL, only : uij=>u,vij=>v,tij=>t,qij=>q,eij=>e
     &     ,dpdxrij=>dpdxr,dpdyrij=>dpdyr
     &     ,dpdxr0ij=>dpdxr0,dpdyr0ij=>dpdyr0
     &     ,zgs,advanc
     &     ,ZS1,TGV,TKV,QG,HEMI,DTSURF,POLE
     &     ,US,VS,WS,WSH,WSQ,TSV,QS,PSI,DBL,KMS,KHS,KQS,PPBL
     &     ,UG,VG,WG,ZMIX
     &     ,ustar,cm,ch,cq,z0m,z0h,z0q
#ifdef TRACERS_ON
     &     ,trij=>tr
      USE TRACER_COM, only : ntm,needtrs,itime_tr0
#endif
      IMPLICIT NONE
     
      INTEGER, INTENT(IN) :: I,J  !@var I,J grid point
      INTEGER, INTENT(IN) :: ITYPE  !@var ITYPE surface type
      REAL*8, INTENT(IN) :: PTYPE  !@var PTYPE percent surface type

      REAL*8 Ts

#ifdef TRACERS_ON
C**** Tracer input/output common block
!@var trsfac, trconstflx factors in surface flux boundary cond.
!@var ntx number of tracers that need pbl calculation
      real*8, dimension(ntm) :: trtop,trs,trsfac,trconstflx
      integer itr,n,ntx
      common /trspec/trtop,trs,trsfac,trconstflx,ntx
#endif

      REAL*8 ztop,zpbl,pl1,tl1,pl,tl,tbar,thbar,zpbl1,coriol
      REAL*8 ttop,qtop,tgrnd,qgrnd,utop,vtop,ufluxs,vfluxs
     *     ,tfluxs,qfluxs,psitop,psisrf
      INTEGER LDC,L,k

C        ocean and ocean ice are treated as rough surfaces
C        roughness lengths from Brutsaert for rough surfaces

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
          ! pk(1,i,j) = expbyk(pl1)
          tl1=t(i,j,1)*(1.+deltx*q(i,j,1))*pk(1,i,j)
          do l=2,ls1
            pl=pmid(l,i,j)        !pij*sig(l)+ptop
            tl=t(i,j,l)*(1.+deltx*q(i,j,l))*pk(l,i,j) !virtual,absolute
            tbar=thbar(tl1,tl)
            zpbl=zpbl-(rgas/grav)*tbar*(pl-pl1)/(pl1+pl)*2.
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
        LDC=nint(DCLEV(I,J))
        IF (LDC.EQ.0) LDC=1
        if (ldc.eq.1) then
          dbl=ztop
          l=1
          else
          zpbl=ztop
          pl1=pmid(1,i,j)          !pij*sig(1)+ptop
          tl1=t(i,j,1)*(1.+deltx*q(i,j,1))*pk(1,i,j)   !expbyk(pl1)
          zpbl1=ztop
          do l=2,ldc
            pl=pmid(l,i,j)         !pij*sig(l)+ptop
            tl=t(i,j,l)*(1.+deltx*q(i,j,l))*pk(l,i,j)  !expbyk(pl)
            tbar=thbar(tl1,tl)
            zpbl=zpbl-(rgas/grav)*tbar*(pl-pl1)/(pl1+pl)*2.
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
      coriol=sinp(j)*omega2
      ttop=tkv
      qtop=q(i,j,1)
      tgrnd=tgv
      qgrnd=qg

      utop=0. ; vtop=0. ;  ug=0. ; vg=0.
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

      uij(:)=uabl(:,i,j,itype)
      vij(:)=vabl(:,i,j,itype)
      tij(:)=tabl(:,i,j,itype)
      qij(:)=qabl(:,i,j,itype)
      eij(1:npbl-1)=eabl(1:npbl-1,i,j,itype)
#ifdef TRACERS_ON
      n=0
      do itr=1,ntm
        if (itime_tr0(itr).le.itime.and.needtrs(itr)) then
          n=n+1
          trij(:,n)=trabl(:,itr,i,j,itype)
        end if
      end do
#endif

      cm=cmgs(i,j,itype)
      ch=chgs(i,j,itype)
      cq=cqgs(i,j,itype)

      dpdxrij  = DPDX_BY_RHO(i,j)
      dpdyrij  = DPDY_BY_RHO(i,j)
      dpdxr0ij = DPDX_BY_RHO_0(i,j)
      dpdyr0ij = DPDY_BY_RHO_0(i,j)

c     write(67,1003) "p-gradients: ",dpdxrij,dpdyrij,dpdxr0ij,dpdyr0ij
c1003 format(a,4(1pe14.4))
      call advanc(
     3     coriol,utop,vtop,ttop,qtop,tgrnd,qgrnd,
#ifdef TRACERS_ON
     *     trs,trtop,trsfac,trconstflx,ntx,
#endif
     4     ztop,dtsurf,ufluxs,vfluxs,tfluxs,qfluxs,i,j,itype)

      uabl(:,i,j,itype)=uij(:)
      vabl(:,i,j,itype)=vij(:)
      tabl(:,i,j,itype)=tij(:)
      qabl(:,i,j,itype)=qij(:)
      eabl(1:npbl-1,i,j,itype)=eij(1:npbl-1)
#ifdef TRACERS_ON
      n=0
      do itr=1,ntm
        if (itime_tr0(itr).le.itime.and.needtrs(itr)) then
          n=n+1
          trabl(:,itr,i,j,itype)=trij(:,n)
        end if
      end do
#endif

      cmgs(i,j,itype)=cm
      chgs(i,j,itype)=ch
      cqgs(i,j,itype)=cq
      ipbl(i,j,itype)=1

      ws    =sqrt(us*us+vs*vs)
      wg    =sqrt(ug*ug+vg*vg)

      psitop=atan2(vg,ug+teeny)
      psisrf=atan2(vs,us+teeny)
      psi   =psisrf-psitop
      ustar_type(i,j,itype)=ustar
C ******************************************************************
      TS=TSV/(1.+QS*deltx)
      WSAVG(I,J)=WSAVG(I,J)+WS*PTYPE
      TSAVG(I,J)=TSAVG(I,J)+TS*PTYPE
      if(itype.ne.4) QSAVG(I,J)=QSAVG(I,J)+QS*PTYPE
      USAVG(I,J)=USAVG(I,J)+US*PTYPE
      VSAVG(I,J)=VSAVG(I,J)+VS*PTYPE
      TAUAVG(I,J)=TAUAVG(I,J)+CM*WS*WS*PTYPE

      uflux(I,J)=uflux(I,J)+ufluxs*PTYPE
      vflux(I,J)=vflux(I,J)+vfluxs*PTYPE
      tflux(I,J)=tflux(I,J)+tfluxs*PTYPE
      qflux(I,J)=qflux(I,J)+qfluxs*PTYPE

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
      USE CONSTANT, only : lhe,lhs,tf,omega2,deltx
      USE MODEL_COM
      USE GEOM, only : idij,idjj,imaxj,kmaxj,rapj,cosiv,siniv,sinp
      USE PBLCOM, ustar_type=>ustar
      USE SOCPBL, only : npbl=>n,zgs,bgrid,inits,ccoeff0
     & ,  uinit=>u,vinit=>v,tinit=>t,qinit=>q,einit=>e
     &     ,dpdxrij=>dpdxr,dpdyrij=>dpdyr
     &     ,dpdxr0ij=>dpdxr0,dpdyr0ij=>dpdyr0
      USE DYNAMICS, only : pmid,pk,pedn,pek
     &    ,DPDX_BY_RHO,DPDY_BY_RHO,DPDX_BY_RHO_0,DPDY_BY_RHO_0
      USE SEAICE_COM, only : rsi
      USE FLUXES, only : gtemp
      USE FILEMANAGER

      IMPLICIT NONE

!@var inipbl whether to init prog vars
      logical, intent(in) :: inipbl
!@var iu_CDN unit number for roughness length input file
      integer :: iu_CDN
      integer :: ilong  !@var ilong  longitude identifier
      integer :: jlat   !@var jlat  latitude identifier
      real*8 tgvdat(im,jm,4)

      integer :: itype  !@var itype surface type
      integer i,j,k,iter,lpbl !@var i,j,k,iter loop variable
      real*8 pland,pwater,plice,psoil,poice,pocean,
     *     ztop,elhx,coriol,tgrnd,pij,ps,psk,qgrnd
     *     ,utop,vtop,qtop,ttop,zgrnd,cm,ch,cq,ustar
      real*8 qsat

C things to be done regardless of inipbl
      call openunit("CDN",iu_CDN,.TRUE.,.true.)
      call readt (iu_CDN,0,roughl,im*jm,roughl,1)
      call closeunit(iu_CDN)

      do j=1,jm
        do i=1,im
          if (fland(i,j).gt.0.and.roughl(i,j).eq.0) then
            print*,"Roughness length not defined for i,j",i,j
     *           ,roughl(i,j),fland(i,j),flice(i,j)
            roughl(i,j)=roughl(10,40)
          end if
        end do
      end do

      call ccoeff0
      call getb(zgs,ztop,bgrid)

      if(.not.inipbl) return

      do j=1,jm
      do i=1,im
        pland=fland(i,j)
        pwater=1.-pland
        plice=flice(i,j)
        psoil=fearth(i,j)
        poice=rsi(i,j)*pwater
        pocean=pwater-poice
        tgvdat(i,j,1)=gtemp(1,1,i,j)+TF
        if (pocean.le.0.) tgvdat(i,j,1)=0.
        tgvdat(i,j,2)=gtemp(1,2,i,j)+TF
        if (poice.le.0.)  tgvdat(i,j,2)=0.
        tgvdat(i,j,3)=gtemp(1,3,i,j)+TF
        if (plice.le.0.)  tgvdat(i,j,3)=0.
        tgvdat(i,j,4)=gtemp(1,4,i,j)+TF
        if (psoil.le.0.)  tgvdat(i,j,4)=0.
      end do
      end do

      do itype=1,4
        if ((itype.eq.1).or.(itype.eq.4)) then
          elhx=lhe
        else
          elhx=lhs
        endif
        do j=1,jm
          jlat=j
          coriol=sinp(j)*omega2

          do i=1,imaxj(j)
            tgrnd=tgvdat(i,j,itype)
            if (tgrnd.eq.0.) then
              ipbl(i,j,itype)=0
              go to 200
            endif
            ilong=i
            pij=p(i,j)
            ps=pedn(1,i,j)    !pij+ptop
            psk=pek(1,i,j)    !expbyk(ps)
            qgrnd=qsat(tgrnd,elhx,ps)

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
            if (itype.gt.2) zgrnd=30./(10.**roughl(i,j))

            dpdxrij  = DPDX_BY_RHO(i,j)
            dpdyrij  = DPDY_BY_RHO(i,j)
            dpdxr0ij = DPDX_BY_RHO_0(i,j)
            dpdyr0ij = DPDY_BY_RHO_0(i,j)

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
      USE MODEL_COM
      USE GEOM, only : imaxj
      USE PBLCOM, only : npbl,uabl,vabl,tabl,qabl,eabl,cmgs,chgs,cqgs
     *     ,ipbl,ustar_type=>ustar
      IMPLICIT NONE
      integer i,j,iter,lpbl  !@var i,j,iter,lpbl loop variable

      do j=1,jm
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
      USE MODEL_COM, only : sige,psf,ptop,psfmpt
      IMPLICIT NONE

      REAL*8, INTENT(IN) :: ZGS
      REAL*8, INTENT(OUT) :: ZTOP,BGRID
      real*8 theta,z1,x

      theta=269.0727251d0
      z1=zgs+0.5*(1.-sige(2))*psfmpt*rgas*theta/(grav*psf)
      x=z1/100.
      ztop=z1
      bgrid=(((0.177427d0*x - 1.0504d0)*x + 2.34169d0)*x -
     2      2.4772d0)*x + 1.44509d0
      return
      end subroutine getb

      SUBROUTINE CHECKPBL(SUBR)
!@sum  CHECKPBL Checks whether PBL data are reasonable
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : im,jm
      USE PBLCOM, only : wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg
     *     ,ustar,uflux,vflux,tflux,qflux
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

      CALL CHECK3(uflux,IM,JM,1,SUBR,'uflux')
      CALL CHECK3(vflux,IM,JM,1,SUBR,'vflux')
      CALL CHECK3(tflux,IM,JM,1,SUBR,'tflux')
      CALL CHECK3(qflux,IM,JM,1,SUBR,'qflux')

      END SUBROUTINE CHECKPBL

