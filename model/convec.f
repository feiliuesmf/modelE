#include "hycom_mpi_hacks.h"
      subroutine convec(m,n,mm,nn,k1m,k1n)
c
c --- hycom version 0.9.1 -- cyclic and noncyclic b.c. combined
c
c --- convective adjustment. either -th3d- and -thstar- can be used as
c --- indicators of static instability. switch between the two by
c --- activating/deactivating lines containing -kappaf-
c
      USE HYCOM_DIM
      USE HYCOM_SCALARS
      USE HYCOM_ARRAYS
      USE DOMAIN_DECOMP, only : HALO_UPDATE,SOUTH,NORTH

      implicit none
      integer,intent(IN) :: m,n,mm,nn,k1m,k1n
      integer i,j,k,l,ja,kn,kp,kbase,kmax
      real q1,q2,sigup,uup,vup,siglo,ulo,vlo,tem,sal,thet,trc(ntrcr),
     .     homog,delp(kdm),dens(kdm),star(kdm),ttem(kdm),ssal(kdm),
     .     trac(kdm,ntrcr),pres(kdm+1),
     .     totem,tosal,tndcyt,tndcys		!  col.integrals (diag.use only)
      real sigocn,kappaf
      external sigocn,kappaf
c
c --- ---------------------
c --- convective adjustment
c --- ---------------------
c
 103  format (i9,2i5,a/(33x,i3,2f8.3,f8.3,f8.2,f8.1))
cdiag if (itest.gt.0 .and. jtest.gt.0) write (lp,103) nstep,itest,jtest,
cdiag.  '  entering convec:  temp    saln    dens    thkns    dpth',
cdiag.  (k,temp(itest,jtest,k+nn),saln(itest,jtest,k+nn),
cdiag.  th3d(itest,jtest,k+nn),dp(itest,jtest,k+nn)/onem,
cdiag.  p(itest,jtest,k+1)/onem,k=1,kk)
c
      CALL HALO_UPDATE(ogrid,th3d  ,FROM=SOUTH)

c$OMP PARALLEL DO PRIVATE(ja,kn,uup,ulo,q1,q2,vup,vlo)
c$OMP+ SCHEDULE(STATIC,jchunk)
      do 26 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
c
c --- convection of u
c
      do 16 k=2,kk
      kn=k+nn
c
      do 16 l=1,isu(j)
      do 16 i=ifu(j,l),ilu(j,l)
      if (th3d(i  ,j,kn).le.th3d(i  ,j,kn-1) .or.
     .    th3d(i-1,j,kn).le.th3d(i-1,j,kn-1)) then
        uup=u(i,j,kn-1)
        ulo=u(i,j,kn  )
        q1=max(dpu(i,j,kn-1),epsil)
        q2=max(dpu(i,j,kn  ),epsil)
        u(i,j,kn  )=(q1*u(i,j,kn-1)+q2*u(i,j,kn))/(q1+q2)
        u(i,j,kn-1)=u(i,j,kn)
cdiag   if (i.eq.itest .and. j.eq.jtest) write (lp,100) nstep,i,j,1,k,
cdiag.    '  upr,lwr,final u:',uup,ulo,u(i,j,kn),q2/(q1+q2)
      end if
 16   continue
c
c --- convection of v
c
      do 26 k=2,kk
      kn=k+nn
c
      do 26 l=1,isv(j)
      do 26 i=ifv(j,l),ilv(j,l)
      if (th3d(i,j  ,kn).le.th3d(i,j  ,kn-1) .or.
     .    th3d(i,ja ,kn).le.th3d(i,ja ,kn-1)) then
        vup=v(i,j,kn-1)
        vlo=v(i,j,kn  )
        q1=max(dpv(i,j,kn-1),epsil)
        q2=max(dpv(i,j,kn  ),epsil)
        v(i,j,kn  )=(q1*v(i,j,kn-1)+q2*v(i,j,kn))/(q1+q2)
        v(i,j,kn-1)=v(i,j,kn)
cdiag   if (i.eq.itest .and. j.eq.jtest) write (lp,100) nstep,i,j,1,k,
cdiag.    '  upr,lwr,final v:',vup,vlo,v(i,j,kn),q2/(q1+q2)
      end if
 26   continue
c$OMP END PARALLEL DO
c
c --- convection of thermodynamic variables and tracer
c
c$OMP PARALLEL DO SCHEDULE(STATIC,jchunk)
c$OMP. PRIVATE(kn,kmax,totem,tosal,homog,sigup,siglo,q1,q2,tem,sal,
c$OMP. trc,thet,kbase,tndcyt,tndcys,ttem,ssal,dens,star,delp,trac,pres)
      do 1 j=J_0,J_1
      do 1 l=1,isp(j)
      do 1 i=ifp(j,l),ilp(j,l)
c
c --- extract single column from 3-d mesh
c
      pres(1)=p(i,j,1)
      do 12 k=1,kk
      kn=k+nn
      ttem(k)=temp(i,j,kn)
      ssal(k)=saln(i,j,kn)
      dens(k)=th3d(i,j,kn)
      star(k)=dens(k)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc     .   +kappaf(ttem(k),ssal(k),pres(k))
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (dotrcr) trac(k,:)=tracer(i,j,k,:)
      delp(k)=dp(i,j,kn)
 12   pres(k+1)=pres(k)+delp(k)
c
      totem=0.
      tosal=0.
      do k=1,kk
        kn=k+nn
        totem=totem+ttem(k)*delp(k)
        tosal=tosal+ssal(k)*delp(k)
      end do
c
      do 4 kbase=1,kk-1
c
      homog=delp(kbase)
      kmax=kbase
      do 6 k=kbase+1,kk
      if (star(k).gt.star(kbase)) exit
      kmax=k
      sigup=dens(kbase)
      siglo=dens(k    )
      q1=max(homog  ,epsil)
      q2=max(delp(k),epsil)
      tem=(q1*ttem(kbase)+q2*ttem(k))/(q1+q2)
      sal=(q1*ssal(kbase)+q2*ssal(k))/(q1+q2)
      thet=sigocn(tem,sal)
      homog=homog+delp(k)
      ttem(kbase)=tem
      ssal(kbase)=sal
      dens(kbase)=thet
      star(kbase)=thet
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc     .   +kappaf(tem,sal,pres(kbase))
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (dotrcr) then
        trc(:)=(q1*trac(kbase,:)+q2*trac(k,:))/(q1+q2)
        trac(kbase,:)=trc(:)
      end if
c
cdiag if (i.eq.itest .and. j.eq.jtest) write (lp,100) nstep,i,j,kbase,
cdiag. k,'  upr,lwr,final dens:',sigup,siglo,
cdiag.  dens(k),q2/(q1+q2)
 100    format (i9,2i5,2i3,a,3f8.3,f5.2)
c
 6    continue
c
      do 7 kp=kbase+1,kmax
      ttem(kp)=tem
      ssal(kp)=sal
      dens(kp)=thet
      star(kp)=thet
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc     .   +kappaf(tem,sal,pres(kp))
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (dotrcr) trac(kp,:)=trc(:)
 7    continue
 4    continue
c
      tndcyt=-totem
      tndcys=-tosal
      totem=10.*pbot(i,j)
      tosal=35.*pbot(i,j)
      do k=1,kk
        tndcyt=tndcyt+ttem(k)*delp(k)
        tndcys=tndcys+ssal(k)*delp(k)
      end do
      if (abs(tndcyt).gt.acurcy*totem)
     .  write (lp,'(2i5,a,1p,2e16.8,e9.1)') i,j,
     .  '  convec - bad temp.intgl.',totem,tndcyt,tndcyt/totem
      if (abs(tndcys).gt.acurcy*tosal)
     .  write (lp,'(2i4,a,1p,2e16.8,e9.1)') i,j,
     .  '  convec - bad saln.intgl.',tosal,tndcys,tndcys/tosal
c
c --- put 1-d column back into 3-d grid
c
      do 8 k=1,kk
      kn=k+nn
      temp(i,j,kn)=ttem(k)
      saln(i,j,kn)=ssal(k)
      th3d(i,j,kn)=dens(k)
      if (dotrcr) tracer(i,j,k,:)=trac(k,:)
 8    continue
c
 1    continue
c$OMP END PARALLEL DO
c
cdiag if (itest.gt.0 .and. jtest.gt.0) write (lp,103) nstep,itest,jtest,
cdiag.  '  exiting  convec:  temp    saln    dens    thkns    dpth',
cdiag.  (k,temp(itest,jtest,k+nn),saln(itest,jtest,k+nn),
cdiag.  th3d(itest,jtest,k+nn),dp(itest,jtest,k+nn)/onem,
cdiag.  p(itest,jtest,k+1)/onem,k=1,kk)
c
      return
      end
c
c
c> Revision history:
c>
c> Oct. 1999 - convection of u and v added
c> May  2000 - switched to single-column structure
c> June 2000 - modified j-1,j+1 to accomodate both channel & closed basin b.c.
c> June 2000 - tightened up kbase-defining logic after loop 10
c> Apr. 2001 - eliminated stmt_funcs.h
c> Nov. 2001 - eliminated multiple passes of -kbase- through k space
c> Nov. 2001 - replaced -th3d- by -thstar- as indicator of static instability
c> Mar. 2002 - allowed mixed layer tracer to fill homogenized column
c> Sep. 2003 - recoded loop 6 (eliminating 'go to 7') to placate ibm compiler
c> Sep. 2003 - added logical switch to enable/disable tracer redistribution
c> Feb. 2005 - added multiple tracer capability
