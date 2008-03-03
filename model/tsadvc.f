#include "hycom_mpi_hacks.h"
      subroutine tsadvc(m,n,mm,nn,k1m,k1n)
c
c --- hycom version 1.0 -- cyclic in j
      USE HYCOM_SCALARS, only : wts1,onemm,wts2,delt1,lp,nstep,itest
     &     ,jtest,diagno,temdff,theta
      USE HYCOM_ARRAYS_GLOB

      USE HYCOM_DIM, only : ii, jj, kk, idm, jdm, kdm,
     &     J_0, J_1, I_0H, J_0H, ogrid,
     &     ip, iu, iv, iq, ifp, ilp, isp, ifq, ilq, isq,
     &     ifu, ilu, isu, ifv, ilv, isv

      use hycom_arrays_glob_renamer
      USE DOMAIN_DECOMP, only : get,AM_I_ROOT,HALO_UPDATE,NORTH,SOUTH,
     &                          GLOBALMIN,GLOBALMAX,ESMF_BCAST

      implicit none
c
!!      include 'dimensions.h'
      include 'dimension2.h'
!!      include 'common_blocks.h'
c
      real smin,smax,tmin,tmax,sminn,smaxx,tminn,tmaxx,posdef,flxdiv
     .    ,offset,factor,q,pold,pmid,pnew,snew,tnew,val
      real uflxn(idm,jdm,kdm),vflxn(idm,jdm,kdm),sign(idm,jdm,kdm),
     .     pn(idm,jdm,kdm+1)
      integer kp
      real sigocn,hfharm
      external sigocn,hfharm
c
      call scatter_hycom_arrays

      do 81 k=1,kk
      km=k+mm
      kn=k+nn
      kp=min(k+1,kk)
c
      smin=999.
      tmin=999.
      smax=-999.
      tmax=-999.
      sminn=999.
      tminn=999.
      smaxx=-999.
      tmaxx=-999.
c
c --- --------------------------------------
c --- advection of thermodynamic variable(s)
c --- --------------------------------------
c
      posdef=100.
c
      call cpy_p_par(dp_loc(I_0H,J_0H,km))
      call cpy_p_par(dp_loc(I_0H,J_0H,kn))
      CALL HALO_UPDATE(ogrid,vflx_loc, FROM=NORTH)

c
c$OMP PARALLEL DO PRIVATE(jb,pold,pmid,flxdiv,offset)
c$OMP+ SCHEDULE(STATIC,jchunk)
      do 49 j=J_0,J_1
      jb= PERIODIC_INDEX(j+1, jj)
      do 49 l=1,isp(j)
      do 49 i=ifp(j,l),ilp(j,l)
c
c --- time smoothing of thermodynamic variable(s) (part 1)
      pold=max(0.,dpold_loc(i,j,k))
      pmid=max(0.,dp_loc(i,j,km))
      temp_loc(i,j,km)=temp_loc(i,j,km)*(wts1*pmid+onemm)+
     .             temp_loc(i,j,kn)* wts2*pold
      saln_loc(i,j,km)=saln_loc(i,j,km)*(wts1*pmid+onemm)+
     .             saln_loc(i,j,kn)* wts2*pold
c
c --- before calling 'advem', make sure (a) mass fluxes are consistent
c --- with layer thickness change, and (b) all fields are positive-definite
      flxdiv=(uflx_loc(i+1,j,k)-uflx_loc(i,j,k)
     .       +vflx_loc(i,jb ,k)-vflx_loc(i,j,k))*delt1*scp2i_loc(i,j)
      util2_loc(i,j)=.5*(dpold_loc(i,j,k)+dp_loc(i,j,kn)-flxdiv)
      util1_loc(i,j)=.5*(dpold_loc(i,j,k)+dp_loc(i,j,kn)+flxdiv)
      offset=min(0.,util1_loc(i,j),util2_loc(i,j))
      util2_loc(i,j)=util2_loc(i,j)-offset
      util1_loc(i,j)=util1_loc(i,j)-offset
c
      temp_loc(i,j,kn)=temp_loc(i,j,kn)+posdef
      smin=min(smin,saln_loc(i,j,kn))
      smax=max(smax,saln_loc(i,j,kn))
      tmin=min(tmin,temp_loc(i,j,kn))
      tmax=max(tmax,temp_loc(i,j,kn))
 49   continue
c$OMP END PARALLEL DO
c
      call GLOBALMIN(ogrid,smin,val); smin=val;
      call GLOBALMAX(ogrid,smax,val); smax=val;
      call GLOBALMIN(ogrid,tmin,val); tmin=val;
      call GLOBALMAX(ogrid,tmax,val); tmax=val;

      call gather_hycom_arrays
      if (AM_I_ROOT()) then

      if (tmin.lt.0. .or. smin.lt.0.) then
      do 490 j=1,jj
      do 490 l=1,isp(j)
      do 490 i=ifp(j,l),ilp(j,l)
      if ((tmin.lt.0. and. temp(i,j,kn).eq.tmin) .or.
     .    (smin.lt.0. and. saln(i,j,kn).eq.smin)) then
        write (lp,101) nstep,i,j,k,' neg. temp/saln bfore advem call ',
     .  temp(i,j,kn)-posdef,saln(i,j,kn)
 101  format (i9,' i,j,k =',2i5,i3,a,2f8.2)
        itest=i
        jtest=j
      end if
 490  continue
!     call stencl(kn)
      end if
c
      call advem(2,temp(1,1,kn),uflx(1,1,k),vflx(1,1,k),scp2,scp2i,
     .           delt1,util1,util2)
c
      call advem(2,saln(1,1,kn),uflx(1,1,k),vflx(1,1,k),scp2,scp2i,
     .           delt1,util1,util2)

      endif ! AM_I_ROOT
      call scatter_hycom_arrays
c
c$OMP PARALLEL DO PRIVATE(pold,pmid,pnew) SCHEDULE(STATIC,jchunk)
      do 46 j=J_0,J_1
      do 46 l=1,isp(j)
      do 46 i=ifp(j,l),ilp(j,l)
      temp_loc(i,j,kn)=temp_loc(i,j,kn)-posdef
      sminn=min(sminn,saln_loc(i,j,kn))
      smaxx=max(smaxx,saln_loc(i,j,kn))
      tminn=min(tminn,temp_loc(i,j,kn))
      tmaxx=max(tmaxx,temp_loc(i,j,kn))
c
c --- time smoothing of thickness field
      pold=max(0.,dpold_loc(i,j,k))
      pmid=max(0.,dp_loc(i,j,km))
      pnew=max(0.,dp_loc(i,j,kn))
      dp_loc(i,j,km)=pmid*wts1+(pold+pnew)*wts2
c --- time smoothing of thermodynamic variable(s) (part 2)
      pmid=max(0.,dp_loc(i,j,km))
      temp_loc(i,j,km)=(temp_loc(i,j,km)+temp_loc(i,j,kn)*wts2*pnew)/
     .   (pmid+onemm)
      saln_loc(i,j,km)=(saln_loc(i,j,km)+saln_loc(i,j,kn)*wts2*pnew)/
     .   (pmid+onemm)
      th3d_loc(i,j,km)=sigocn(temp_loc(i,j,km),saln_loc(i,j,km))
c
c --- build up time integral of mass field variables
      pmid=max(0.,dp_loc(i,j,km))
      dpav_loc (i,j,k)=dpav_loc (i,j,k)+pmid
      temav_loc(i,j,k)=temav_loc(i,j,k)+temp_loc(i,j,km)*pmid
      salav_loc(i,j,k)=salav_loc(i,j,k)+saln_loc(i,j,km)*pmid
      th3av_loc(i,j,k)=th3av_loc(i,j,k)+th3d_loc(i,j,km)*pmid
 46   continue
c$OMP END PARALLEL DO
c
      call GLOBALMIN(ogrid,sminn,val); sminn=val;
      call GLOBALMAX(ogrid,smaxx,val); smaxx=val;
      call GLOBALMIN(ogrid,tminn,val); tminn=val;
      call GLOBALMAX(ogrid,tmaxx,val); tmaxx=val;

      call gather_hycom_arrays
      if (AM_I_ROOT()) then

      if (tminn+posdef.lt.0. .or. sminn.lt.0.) then
      do 492 j=1,jj
      do 492 l=1,isp(j)
      do 492 i=ifp(j,l),ilp(j,l)
      if (temp(i,j,kn).eq.tminn .or. saln(i,j,kn).eq.sminn)
     .  write (lp,101) nstep,i,j,k,' neg. temp/saln after advem call ',
     .  temp(i,j,kn),saln(i,j,kn)
 492  continue
      end if
c
cdiag if (itest.gt.0.and.jtest.gt.0)
cdiag.write (lp,'(i9,2i5,i3,'' th,s,dp after advection  '',2f9.3,f8.2)')
cdiag.nstep,itest,jtest,k,temp(itest,jtest,kn),saln(itest,jtest,kn),
cdiag.dp(itest,jtest,kn)/onem
c
      if (diagno)
     . write (lp,'(i9,i3,'' min/max of s after advection:'',4f7.2)')
     . nstep,k,sminn,smaxx
c
c --- --------------------------------------
c --- diffusion of thermodynamic variable(s)
c --- --------------------------------------
c
      endif ! AM_I_ROOT
      call scatter_hycom_arrays

      call cpy_p_par(temp_loc(I_0H,J_0H,kn))
      call cpy_p_par(saln_loc(I_0H,J_0H,kn))
      call cpy_p_par(th3d_loc(I_0H,J_0H,kn))

      !south halo updates on saln, temp, dp
      CALL HALO_UPDATE(ogrid,saln_loc, FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(ogrid,temp_loc, FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(ogrid,  dp_loc, FROM=SOUTH+NORTH)

c
c$OMP PARALLEL DO PRIVATE(ja,factor) SCHEDULE(STATIC,jchunk)
      do 145 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
c
      do 144 l=1,isu(j)
      do 144 i=ifu(j,l),ilu(j,l)
      factor=scuy_loc(i,j)*2.*hfharm(max(dp_loc(i-1,j,kn),onemm)
     .                          ,max(dp_loc(i  ,j,kn),onemm))
      uflux_loc (i,j)=factor*(temp_loc(i-1,j,kn)-temp_loc(i,j,kn))
 144  uflux2_loc(i,j)=factor*(saln_loc(i-1,j,kn)-saln_loc(i,j,kn))
c
      do 145 l=1,isv(j)
      do 145 i=ifv(j,l),ilv(j,l)
      factor=scvx_loc(i,j)*2.*hfharm(max(dp_loc(i,ja ,kn),onemm)
     .                          ,max(dp_loc(i,j  ,kn),onemm))
      vflux_loc (i,j)=factor*(temp_loc(i,ja ,kn)-temp_loc(i,j,kn))
 145  vflux2_loc(i,j)=factor*(saln_loc(i,ja ,kn)-saln_loc(i,j,kn))
c$OMP END PARALLEL DO
c
      CALL HALO_UPDATE(ogrid,vflux_loc, FROM=NORTH)
      CALL HALO_UPDATE(ogrid,vflux2_loc, FROM=NORTH)

c$OMP PARALLEL DO PRIVATE(jb,factor) SCHEDULE(STATIC,jchunk)
      do 146 j=J_0,J_1
      jb = PERIODIC_INDEX(j+1, jj)
      do 146 l=1,isp(j)
      do 146 i=ifp(j,l),ilp(j,l)
      factor=-temdff*delt1/(scp2_loc(i,j)*max(dp_loc(i,j,kn),onemm))
      util1_loc(i,j)=(uflux_loc (i+1,j)-uflux_loc (i,j)
     .           +vflux_loc (i,jb )-vflux_loc (i,j))*factor
      util2_loc(i,j)=(uflux2_loc(i+1,j)-uflux2_loc(i,j)
     .           +vflux2_loc(i,jb )-vflux2_loc(i,j))*factor
      temp_loc(i,j,kn)=temp_loc(i,j,kn)+util1_loc(i,j)
      saln_loc(i,j,kn)=saln_loc(i,j,kn)+util2_loc(i,j)
      th3d_loc(i,j,kn)=sigocn(temp_loc(i,j,kn),saln_loc(i,j,kn))
c
cdiag if (i.eq.itest.and.j.eq.jtest)
cdiag. write (lp,100) nstep,i,j,k,'t,s,dt,ds,dsigdt,dsigds,cabbl =',
cdiag. temp(i,j,kn),saln(i,j,kn),util1(i,j),util2(i,j),
cdiag. dsigdt(temp(i,j,kn),saln(i,j,kn))*util1(i,j),
cdiag. dsigds(temp(i,j,kn),saln(i,j,kn))*util2(i,j),cabbl(i,j,k)
c
 146  continue
c$OMP END PARALLEL DO
c
cdiag if (itest.gt.0.and.jtest.gt.0)
cdiag.write (lp,'(i9,2i5,i3,'' t,s,dp after isopyc.mix.'',2f9.3,f8.2)')
cdiag.nstep,itest,jtest,k,temp(itest,jtest,kn),saln(itest,jtest,kn),
cdiag.dp(itest,jtest,kn)/onem
c
      call cpy_p_par(temp_loc(I_0H,J_0H,km))
      call cpy_p_par(temp_loc(I_0H,J_0H,kn))
      call cpy_p_par(saln_loc(I_0H,J_0H,km))
      call cpy_p_par(saln_loc(I_0H,J_0H,kn))
      call cpy_p_par(th3d_loc(I_0H,J_0H,km))
      call cpy_p_par(th3d_loc(I_0H,J_0H,kn))

 81   continue

      call gather_hycom_arrays
      if (AM_I_ROOT()) then
c
c --- convert mass fluxes to density coord. prior to time integration
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call reflux(uflx ,vflx ,th3d(1,1,k1m),p,
     .            uflxn,vflxn,sign,pn,theta,kdm,kdm)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      endif ! AM_I_ROOT
      call scatter_hycom_arrays
      call esmf_bcast(ogrid,uflxn)
      call esmf_bcast(ogrid,vflxn)

c --- activate this loop if -reflux- is   n o t   called
ccc      do k=1,kk
ccc      do j=1,jj
ccc      do i=1,ii
ccc      uflxn(i,j,k)=uflx(i,j,k)
ccc      vflxn(i,j,k)=vflx(i,j,k)
ccc      end do
ccc      end do
ccc      end do
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c$OMP PARALLEL DO PRIVATE(ja) SCHEDULE(STATIC,jchunk)
      do 155 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
      do 155 k=1,kk
c
      do 154 l=1,isu(j)
      do 154 i=ifu(j,l),ilu(j,l)
      uflxav_loc(i,j,k)=uflxav_loc(i,j,k)+uflxn(i,j,k)		!  uflx time integral
     .   *.5*min(nstep,2)
 154  continue
c
      do 155 l=1,isv(j)
      do 155 i=ifv(j,l),ilv(j,l)
      vflxav_loc(i,j,k)=vflxav_loc(i,j,k)+vflxn(i,j,k)		!  vflx time integral
     .   *.5*min(nstep,2)
 155  continue
c$OMP END PARALLEL DO
c
      call gather_hycom_arrays

      return
      end
c
c
c  Revision history:
c
c> June 1995 - eliminated setting of salinity in massless layers (loop 46)
c>             (this is now done in mxlayr.f)
c> Aug. 1995 - added array -cabbl- to transmit cabbeling info to -diapfl-
c> Aug. 1995 - omitted t/s/dp time smoothing in case of abrupt mxlayr.thk.change
c> Sep. 1995 - increased temdff if mixed layer occupies >90% of column
