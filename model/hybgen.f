      subroutine hybgen(m,n,mm,nn,k1m,k1n)
c
c --- hycom version 0.9.10
c --- this version allows switching between T/S and rho/S conservation
c
      implicit none
c
c --- ---------------------
c --- hybrid grid generator (coordinate restoration exclusively by "dilution")
c --- ---------------------
c
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
c
      real sigref,delp,dp0,dp0abv,dpsum,zinteg,tinteg,sinteg,uvintg,
     .     uvscl,siga,sigb,za,zb,phi,plo,pa,pb,dsgdt,dsgds,scalt,scals,
     .     tup,sup,tem,sal,thhat,s_hat,t_hat,p_hat,q,q1,q2,slak,
     .     torho,totem,tosal,tndcyr,tndcyt,tndcys,displ(kdm+1)
      real targt(kdm+1),dens(kdm),ttem(kdm),ssal(kdm),pres(kdm+1),
     .     uold(kdm),vold(kdm),pold(kdm+1),pnew(kdm+1),dplist(kdm),
     .     trac(kdm),trcint
      logical abort,tscnsv,vrbos
      data tscnsv/.true./		! if true, go with T/S conservation
      data abort/.false./
      real sigocn,tofsig,dsigdt,dsigds,cushn
      external sigocn,tofsig,dsigdt,dsigds,cushn
      integer lyr,k1,kp,iunit,lpunit,ko,khyb
     .        ,ntot2,ntot3,nwrk2,nwrk3
     .        ,ntot2d(jdm),ntot3d(jdm),nwrk2d(jdm),nwrk3d(jdm)
      common/nmpipe/iunit,lpunit
      character text*20
      data uvscl/0.02/			!  2 cm/s
      data scalt,scals/30.,10./
      parameter (slak=.5/86400.)	! intfc nudging time scale: 2 days
ccc   parameter (slak=1./86400.)	! intfc nudging time scale: 1 day
ccc   parameter (slak=1.e6)		! intfc nudging time scale: 1 microsec
c
      data (dplist(k),k=1,kdm)/
 
     .   20.0, 5.0, 7.6, 9.8,11.6,13.0,14.0,14.6,
     .   14.9,15.0,15.0,15.0,15.0,15.0,15.0,15.0/

ccc     .    5.0, 8.8,12.2,15.2,17.8,20.0,21.9,23.6,
ccc     .   25.1,26.4,27.5,28.4,29.1,29.6,29.9,30.0/	!  350.5

ccc     .    5.0, 8.8,12.3,15.5,18.4,21.0,23.3,25.3,
ccc     .   27.0,28.4,29.5,30.4,31.1,31.6,31.9,32.0/	!  371.5

ccc     .    5.0, 8.9,12.5,15.8,18.8,21.5,23.9,26.0,
ccc     .   27.8,29.3,30.5,31.4,32.1,32.6,32.9,33.0/	!  382.0

ccc     .    5.0, 9.3,13.2,16.7,19.8,22.5,24.9,27.0,
ccc     .   28.8,30.3,31.5,32.4,33.1,33.6,33.9,34.0/	!  396.0

ccc     .    5.0, 9.6,13.8,17.6,21.0,24.0,26.6,28.8,
ccc     .   30.6,32.0,33.1,33.9,34.4,34.7,34.9,35.0/	!  415.0

ccc     .    5.0, 9.6,13.8,17.6,21.0,24.0,26.6,28.8,
ccc     .   30.7,32.3,33.6,34.6,35.3,35.7,35.9,36.0/	!  420.5

ccc     .    5.0,10.0,14.6,18.8,22.6,26.0,29.0,31.6,
ccc     .   33.8,35.6,37.0,38.0,38.6,38.9,39.0,39.0/    !  457.5

ccc     .   20.0, 5.0, 7.6, 9.8,11.6,13.0,14.0,14.6,14.9,15.0,
ccc     .   15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0/	!  275.5 total

ccc     .   20.0, 5.0, 8.1,10.9,13.5,15.9,18.1,20.1,21.9,23.5,
ccc     .   24.9,26.1,27.1,27.9,28.5,29.0,29.4,29.7,29.9,30.0/  !  439.5

ccc     .    5.0, 8.3,11.2,13.8,16.2,18.4,20.4,22.2,23.8,25.2,
ccc     .   26.4,27.4,28.2,28.9,29.5,30.0,30.4,30.7,30.9,31.0/  !  457.9

ccc     .    5.0, 8.5,11.6,14.4,17.0,19.4,21.6,23.6,25.4,27.0,
ccc     .   28.4,29.6,30.6,31.4,32.0,32.4,32.7,32.9,33.0,33.0/  !  489.5

ccc     .    5.0, 8.6,11.9,14.9,17.7,20.3,22.7,24.9,26.9,28.7,
ccc     .   30.3,31.7,32.9,33.9,34.7,35.3,35.7,35.9,36.0,36.0/  !  524.0

ccc     .    5.0, 8.8,12.3,15.5,18.4,21.0,23.4,25.6,27.6,29.4,
ccc     .   31.0,32.4,33.6,34.6,35.4,36.0,36.4,36.7,36.9,37.0/  !  537.0

ccc     .    5.0, 8.8,12.3,15.5,18.4,21.1,23.6,25.9,28.0,29.9,
ccc     .   31.6,33.1,34.4,35.5,36.4,37.1,37.6,37.9,38.0,38.0/  !  548.1

ccc     .    5.0, 9.0,12.6,15.9,18.9,21.7,24.3,26.7,28.9,30.9,
ccc     .   32.7,34.3,35.7,36.9,37.9,38.7,39.3,39.7,39.9,40.0/  !  569.0

ccc     .    5.0, 9.2,13.0,16.4,19.5,22.4,25.1,27.6,29.9,32.0,
ccc     .   33.9,35.6,37.1,38.4,39.5,40.4,41.1,41.6,41.9,42.0/  !  591.6

ccc     .    5.0, 9.5,13.6,17.4,20.9,24.1,27.0,29.6,31.9,33.9,
ccc     .   35.7,37.3,38.7,39.9,40.9,41.7,42.3,42.7,42.9,43.0/  !  618.0

ccc     .    5.0, 9.5,13.6,17.4,20.9,24.1,27.0,29.6,31.9,34.0,
ccc     .   35.9,37.6,39.1,40.4,41.5,42.4,43.1,43.6,43.9,44.0/  !  624.5

ccc     .    5.0, 9.7,14.0,17.9,21.5,24.8,27.8,30.5,32.9,35.0,
ccc     .   36.9,38.6,40.1,41.4,42.5,43.4,44.1,44.6,44.9,45.0/  !  640.6

ccc     .    5.0, 9.8,14.2,18.3,22.1,25.6,28.8,31.7,34.3,36.6,
ccc     .   38.6,40.3,41.7,42.9,43.9,44.7,45.3,45.7,45.9,46.0/  !  661.4

ccc     .    5.0,10.0,14.6,18.8,22.7,26.3,29.6,32.6,35.3,37.7,
ccc     .   39.8,41.6,43.1,44.3,45.2,45.9,46.4,46.7,46.9,47.0/  !  679.5

ccc     .    5.0,10.0,14.6,18.9,22.9,26.6,30.0,33.1,35.9,38.4,
ccc     .   40.6,42.5,44.1,45.4,46.4,47.1,47.6,47.9,48.0,48.0/  !  693.0

ccc     .    5.0,10.0,14.6,18.9,22.9,26.6,30.0,33.1,35.9,38.4,
ccc     .   40.6,42.5,44.1,45.4,46.5,47.4,48.1,48.6,48.9,49.0/  !  696.5

ccc     .    5.0,10.1,14.9,19.4,23.6,27.5,31.1,34.4,37.4,40.1,
ccc     .   42.5,44.6,46.4,47.9,49.1,50.0,50.6,50.9,51.0,51.0/  !  727.5

ccc     .    5.0,10.3,15.2,19.8,24.1,28.1,31.8,35.2,38.3,41.1,
ccc     .   43.6,45.8,47.7,49.3,50.6,51.6,52.3,52.7,52.9,53.0/  !  748.4

ccc     .    5.0,11.0,16.6,21.8,26.6,31.0,35.0,38.6,41.8,44.6,
ccc     .   47.0,49.0,50.6,51.8,52.7,53.3,53.7,53.9,54.0,54.0/  !  792.0

ccc     .    5.0,11.0,16.6,21.8,26.6,31.0,35.0,38.6,41.8,44.6,
ccc     .   47.0,49.0,50.6,51.9,52.9,53.7,54.3,54.7,54.9,55.0/  !  796.0

ccc     .    5.0,11.2,17.0,22.4,27.4,32.0,36.2,40.0,43.4,46.4,
ccc     .   49.0,51.2,53.0,54.4,55.4,56.1,56.6,56.9,57.0,57.0/  !  827.6

ccc     .    5.0,11.2,17.0,22.4,27.4,32.0,36.2,40.0,43.4,46.4,
ccc     .   49.0,51.2,53.0,54.4,55.5,56.4,57.1,57.6,57.9,58.0/  !  831.1

ccc     .    5.0,11.2,17.0,22.4,27.4,32.0,36.2,40.0,43.4,46.4,
ccc     .   49.1,51.5,53.6,55.4,56.9,58.1,59.0,59.6,59.9,60.0/  !  844.1

ccc     .    5.0,11.5,17.6,23.3,28.6,33.5,38.0,42.1,45.8,49.1,
ccc     .   52.0,54.5,56.6,58.3,59.6,60.6,61.3,61.7,61.9,62.0/  !  883.0

ccc     .    5.0,11.5,17.6,23.3,28.6,33.5,38.0,42.1,45.8,49.1,
ccc     .   52.0,54.5,56.6,58.4,59.9,61.1,62.0,62.6,62.9,63.0/  !  887.5

ccc     .    5.0,12.0,18.6,24.8,30.6,36.0,41.0,45.6,49.8,53.6,
ccc     .   57.0,60.0,62.6,64.8,66.6,68.0,69.0,69.6,69.9,70.0/  !  974.5
c
      sigref=1000.*thref
c
cdiag if (itest.gt.0 .and. jtest.gt.0) then
cdiag   write (lp,103) nstep,itest,jtest,
cdiag.  '  entering hybgen:  temp    saln    dens    thkns    dpth',
cdiag.  (k,temp(itest,jtest,k+nn),saln(itest,jtest,k+nn),sigref*
cdiag.  (th3d(itest,jtest,k+nn)+thbase),dp(itest,jtest,k+nn)/onem,
cdiag.  p(itest,jtest,k+1)/onem,k=1,kk)
cdiag   write (lp,106) nstep,itest,jtest,
cdiag.  '  entering hybgen:  dpthu      u    dpthv      v',
cdiag.  (k,pu(itest,jtest,k+1)/onem,u(itest,jtest,k+nn),
cdiag.     pv(itest,jtest,k+1)/onem,v(itest,jtest,k+nn),k=1,kk)
cdiag end if
 103  format (i9,2i5,a/(33x,i3,2f8.3,f8.3,f8.2,f8.1))
 106  format (i9,2i5,a/(33x,i3,2(f8.1,f8.3)))
c
c$OMP PARALLEL DO
      do 19 j=1,jj
      do 19 k=1,kk
      do 19 l=1,isp(j)
      do 19 i=ifp(j,l),ilp(j,l)
 19   p(i,j,k+1)=p(i,j,k)+dp(i,j,k+nn)
c$OMP END PARALLEL DO
c
      abort=.false.
c$OMP PARALLEL DO PRIVATE(torho,totem,tosal,kp,q,q1,q2,tem,sal,dsgdt,
c$OMP+ dsgds,t_hat,s_hat,thhat,tup,sup,p_hat,tndcyr,tndcyt,tndcys,kn,
c$OMP+ dens,ttem,ssal,pres,targt,dp0,dp0abv,dpsum,k1,tinteg,sinteg,phi,
c$OMP+ plo,pa,pb,ntot2,ntot3,nwrk2,nwrk3,trac,trcint,displ,vrbos,khyb)
c$OMP+ SHARED(abort)
      do 12 j=1,jj
      ntot2=0
      ntot3=0
      nwrk2=0
      nwrk3=0
c
      do 2 l=1,isp(j)
      do 2 i=ifp(j,l),ilp(j,l)
c
c --- extract t,s,rho column from 3-d grid 
c
      pres(1)=p(i,j,1)
      do 3 k=1,kk
      kn=k+nn
      dens(k)=th3d(i,j,kn)
      ttem(k)=temp(i,j,kn)
      ssal(k)=saln(i,j,kn)
      pres(k+1)=pres(k)+dp(i,j,kn)
 3    targt(k)=theta(k)
c
      vrbos=.false.
cdiag if (i.eq.itest .and. j.eq.jtest) vrbos=.true.
      if (vrbos) then
        write (lp,99) nstep,i,j,'      o l d   p r o f i l e :'
        do k=1,kk,10
        write (lp,100) (pres(k1)/onem,k1=k,min(kk+1,k+10))
        write (lp,101) (sigref*(dens(k1)+thbase),k1=k,min(kk,k+9))
        write (lp,102) (ttem(k1),k1=k,min(kk,k+9))
        write (lp,102) (ssal(k1),k1=k,min(kk,k+9))
        end do
      end if
 99   format (i9,2i5,a)
 100  format (11f7.1)
 101  format (4x,10f7.2)
 102  format (4x,10f7.2)
c
      torho=0.
      totem=0.
      tosal=0.
      do k=1,kk
        torho=torho+dens(k)*(pres(k+1)-pres(k))
        totem=totem+ttem(k)*(pres(k+1)-pres(k))
        tosal=tosal+ssal(k)*(pres(k+1)-pres(k))
      end do
c
      kp=1
      do 4 k=2,kk
      if (pres(k).lt.pres(kk+1)-onecm) then
        kp=k
      else
c
c --- fill massless layers on sea floor with data from above
c
        if (vrbos)
     .    write (lp,'(i9,2i5,a,i3,a,i3/(10f8.4))') nstep,i,j,
     .    '  absorb layer',k,' in',kp,(pres(k1+1)-pres(k1),k1=kp+1,kk)
c
        q1=max(epsil,pres(k  )-pres(kp))
        q2=max(   0.,pres(k+1)-pres(k ))
        q=q1/(q1+q2)
        if (q.lt.0. .or. q.gt.1.) then
          write (lp,*) 'i,j,q1,q2,q=',i,j,q1,q2,q
          abort=.true.
        end if
        ttem(kp)=ttem(kp)*q+ttem(k)*(1.-q)
        ssal(kp)=ssal(kp)*q+ssal(k)*(1.-q)
        dens(kp)=dens(kp)*q+dens(k)*(1.-q)
        ttem(k)=ttem(kp)
        ssal(k)=max(ssal(kp),salmin(k))
        dens(k)=dens(kp)
      end if
 4    continue
c
      do 11 k=kp+1,kk
 11   pres(k)=pres(kk+1)
c
      do 16 k=1,kk
 16   dpold(i,j,k)=pres(k+1)-pres(k)
c
      do 10 k=3,kk
c
c --- is layer touching sea floor to light?
      if (pres(k  ).lt.pres(kk+1)-onecm .and.
     .    pres(k+1).eq.pres(kk+1) .and. dens(k-1).lt.dens(k) .and.
     .    dens(k).lt.theta(k)-.001*sigjmp) then
c
c --- water in layer k is too light. split layer into 2 sublayers, one
c --- having the desired density and one matching the density of layer k-1.
c --- t & s in sublayers are obtained by unmixing in direction perpendicular
c --- to isopycnals in t/s space. combine upper sublayer with layer k-1.
c
        tem=ttem(k)
        sal=ssal(k)
ccc        scalt=abs(ttem(k-1)-tem)
ccc        scals=abs(ssal(k-1)-sal)
ccc        if (scalt+scals.eq.0.) go to 10
        dsgdt=dsigdt(tem,sal)*scalt
        dsgds=dsigds(tem,sal)*scals
        q=1./(dsgdt*dsgdt+dsgds*dsgds)
c
c --- set properties in lower sublayer:
        s_hat=sal+(theta(k)-dens(k))*q*dsgds*scals
        if (tscnsv) then
          t_hat=tem+(theta(k)-dens(k))*q*dsgdt*scalt
          thhat=sigocn(t_hat,s_hat)-thbase
        else
          thhat=theta(k)
          t_hat=tofsig(thhat+thbase,s_hat)
        end if
c
c --- set properties in upper sublayer:
        sup=sal+(dens(k-1)-dens(k))*q*dsgds*scals
        if (tscnsv) then
          tup=tem+(dens(k-1)-dens(k))*q*dsgdt*scalt
          if (dsgdt*scalt+dsgds*scals.gt.0.) then
            q=(s_hat-sal)/(s_hat-sup)
          else
            q=(t_hat-tem)/(t_hat-tup)
          end if
        else
          q=(thhat-dens(k))/max(epsil,thhat-dens(k-1))
        end if
        if (q.ge.0. .and. q.le.1.) then
          p_hat=pres(k)*(1.-q)+pres(k+1)*q
c
          if (vrbos) then
            tup=tofsig(dens(k-1)+thbase,sup)
            write (lp,'(i9,2i5,i3,a,3f7.3,f8.2)') nstep,i,j,k,
     .      '  t,s,th,dp in upper sblyr:',tup,sup,
     .      sigref*sigocn(tup,sup),(p_hat-pres(k))/onem
            write (lp,'(22x,a,3f7.3,f8.2)') 
     .      '  t,s,th,dp in lower sblyr:',t_hat,s_hat,
     .      sigref*(thhat+thbase),(pres(k+1)-p_hat)/onem
            write (lp,'(22x,a,1p,2e11.3)') '  scalt,scals =',scalt,scals
          end if
c
c --- combine upper sublayer with layer k-1
          q=(p_hat-pres(k))/max(p_hat-pres(k-1),epsil)
          if (q.lt.0. .or. q.gt.1.) then
            write (lp,*) 'i,j,k,p_hat,pres(k),q=',
     .                    i,j,k,p_hat,pres(k),q
            abort=.true.
          end if
          ssal(k-1)=sup    *q+ssal(k-1)*(1.-q)
          if (tscnsv) then
            ttem(k-1)=tup*q+ttem(k-1)*(1.-q)
            dens(k-1)=sigocn(ttem(k-1),ssal(k-1))-thbase
          else
            ttem(k-1)=tofsig(dens(k-1)+thbase,ssal(k-1))
          end if
c
          if (vrbos)
     .    write (lp,'(22x,a,2f7.3,0p,f8.2)') '  old/new th(k-1):',
     .    sigref*(dens(k-1)+thbase),sigref*sigocn(ttem(k-1),ssal(k-1)),q
c
          ttem(k)=t_hat
          ssal(k)=s_hat
          dens(k)=thhat
          pres(k)=p_hat
        end if
        nwrk2=nwrk2+1
      end if
 10   continue
      ntot2=ntot2+1
c
c --- try to restore isopycnic conditions by moving layer interfaces
c
      khyb=kk
      dpsum=0.
      do 8 k=1,kk
      ntot3=ntot3+1
c
c --- set lower limits for layer thknss (dp0) and depth of lower intfc (dpsum)
c
      dp0abv=dp0
      dp0=dplist(k)*onem
c
c --- optional: reduce spacing of z layers near equator, but hide transition
c --- in a subtropical latitude band where z layers are least likely to exist
      if (k.gt.1) dp0=dp0*max(.6,min(1.,(abs(latij(i,j,3))+5.)*.04))
c
c --- shrink layer thickness in shallow regions, mimicking sigma coordinate
      dp0=min(dp0,5.*pbot(i,j)/float(kk))
      dpsum=dpsum+dp0
c
c --- maintain constant thickness in layer 1
      if (k.eq.1) then
        p_hat=dp0
        if (p_hat.gt.pres(2)) then
c --- layer 1 is too thin. entrain water from layers below
          p_hat=min(p_hat,pres(2)+
     .          max(onecm,10.*slak*delt1*(p_hat-pres(2))))
          go to 5
        else if (p_hat.lt.pres(2)) then
c --- layer 1 is too thick. expell layer 1 water into layer 2
          p_hat=max(p_hat,pres(2)+
     .          min(-onecm,10.*slak*delt1*(p_hat-pres(2))))
c
          if (vrbos) write (lp,105)
     .     i,j,k,'lower intfc',pres(k+1)/onem,'=>',p_hat/onem
 105      format (2i5,i3,2(3x,a,f9.3))
c
          q=(pres(2)-p_hat)/max(pres(3)-p_hat,epsil)
          if (q.lt.0. .or. q.gt.1.) then
            write (lp,*) 'i,j,k,pres(2),p_hat,q=',
     .                    i,j,k,pres(2),p_hat,q
            abort=.true.
          end if
          ssal(2)=ssal(2)*(1.-q)+ssal(1)*q
          if (tscnsv) then
            ttem(2)=ttem(2)*(1.-q)+ttem(1)*q
            dens(2)=sigocn(ttem(2),ssal(2))-thbase
          else
            dens(2)=dens(2)*(1.-q)+dens(1)*q
            ttem(2)=tofsig(dens(2)+thbase,ssal(2))
          end if
          pres(2)=p_hat
          nwrk3=nwrk3+1
        end if
        khyb=1				!  layer 1 is in hybrid domain
        go to 8
      end if				!  k = 1
c
c --- are we dealing with a near-massless layer on the sea floor?
      if (pres(k).eq.pres(kk+1)) dens(k)=max(targt(k),dens(k))
      if (pres(k).gt.pres(kk+1)-onecm) go to 8
c
c --- is lower intfc too close to the surface?
      p_hat=dpsum
      if (k.lt.kk .and. p_hat.gt.pres(k+1)) then
        p_hat=min(p_hat,pres(k+1)+
     .        max(onecm,slak*delt1*(p_hat-pres(k+1))))
        khyb=k				!  layer k is in hybrid domain
        go to 5
      end if
c
c --- is density noticeably different from target value?
      if (abs(dens(k)-targt(k)).lt..1*sigjmp) go to 8
c
      if (dens(k).le.targt(k)) go to 7		!  layer too light
c
c --- water in layer k is too  d e n s e . dilute with water from layer k-1
c                              ^^^^^^^^^
      if (k.eq.2) go to 6			!  don't touch layer 1
      q=(targt(k)-dens(k))/max(targt(k)-dens(k-1),sigjmp)
      p_hat=pres(k)*(1.-q)+pres(k+1)*q
c
c --- maintain minimum layer thickess of layer k-1
      p_hat=pres(k-1)+cushn(p_hat-pres(k-1),dp0abv)
      p_hat=min(p_hat,.5*(pres(k-1)+pres(k+1)))
c
      if (vrbos) write (lp,105)
     . i,j,k,'upper intfc',pres(k)/onem,'=>',p_hat/onem
c
      if (p_hat.lt.pres(k)) then
c
c --- upper intfc moves up. entrain layer k-1 water into layer k
c
        p_hat=max(p_hat,pres(k-1),pres(k)+
     .        min(-onecm,slak*delt1*(p_hat-pres(k))))
        if (k-1.le.khyb) then			!  use plm
          displ(1)=0.
          displ(2)=0.
          displ(3)=p_hat-pres(k)
          displ(4)=0.
          if (vrbos)
     .    write (lp,'(2i5,i3,a)') i,j,k,'  entrain from layer above'
          call plmad3(pres(k-2),displ,ssal(k-2),ssal(k-2),vrbos)
          if (tscnsv) then
            call plmad3(pres(k-2),displ,ttem(k-2),ttem(k-2),vrbos)
            dens(k-1)=sigocn(ttem(k-1),ssal(k-1))-thbase
            dens(k  )=sigocn(ttem(k  ),ssal(k  ))-thbase
          else
            call plmad3(pres(k-2),displ,dens(k-2),dens(k-2),vrbos)
            ttem(k-1)=tofsig(dens(k-1)+thbase,ssal(k-1))
            ttem(k  )=tofsig(dens(k  )+thbase,ssal(k  ))
          end if
        else					!  use pcm
          q=(pres(k)-p_hat)/max(pres(k+1)-p_hat,epsil)
          if (q.lt.0. .or. q.gt.1.) then
            write (lp,*) 'i,j,k,p_hat,pres(k),q=',
     .                    i,j,k,p_hat,pres(k),q
            abort=.true.
          end if
          ssal(k)=ssal(k)*(1.-q)+ssal(k-1)*q
          if (tscnsv) then
            ttem(k)=ttem(k)*(1.-q)+ttem(k-1)*q
            dens(k)=sigocn(ttem(k),ssal(k))-thbase
          else
            dens(k)=dens(k)*(1.-q)+dens(k-1)*q
            ttem(k)=tofsig(dens(k)+thbase,ssal(k))
          end if
        end if
        pres(k)=p_hat
        nwrk3=nwrk3+1
c
      else if (p_hat.gt.pres(k)) then		!  p_hat > pres(k)
c
c --- layer k-1 is too thin for allowing upper intfc to move up.  instead,
c --- move upper interface down and entrain layer k water into layer k-1
c
        p_hat=min(p_hat,pres(k+1),pres(k)+
     .        max(onecm,slak*delt1*(p_hat-pres(k))))
        khyb=k-1			!  layer k-1 is in hybrid domain
        if (k.lt.kk) then		!  use plm
          if (vrbos)
     .    write (lp,'(2i5,i3,a)') i,j,k,'  detrain into layer above'
          displ(1)=0.
          displ(2)=p_hat-pres(k)
          displ(3)=0.
          displ(4)=0.
          call plmad3(pres(k-1),displ,ssal(k-1),ssal(k-1),vrbos)
          if (tscnsv) then
            call plmad3(pres(k-1),displ,ttem(k-1),ttem(k-1),vrbos)
            dens(k-1)=sigocn(ttem(k-1),ssal(k-1))-thbase
            dens(k  )=sigocn(ttem(k  ),ssal(k  ))-thbase
          else
            call plmad3(pres(k-1),displ,dens(k-1),dens(k-1),vrbos)
            ttem(k-1)=tofsig(dens(k-1)+thbase,ssal(k-1))
            ttem(k  )=tofsig(dens(k  )+thbase,ssal(k  ))
          end if
        else					!  use pcm
          q=(p_hat-pres(k))/max(p_hat-pres(k-1),epsil)
          if (q.lt.0. .or. q.gt.1.) then
            write (lp,*) 'i,j,k,p_hat,pres(k),q=',
     .                    i,j,k,p_hat,pres(k),q
            abort=.true.
          end if
          ssal(k-1)=ssal(k-1)*(1.-q)+ssal(k)*q
          if (tscnsv) then
            ttem(k-1)=ttem(k-1)*(1.-q)+ttem(k)*q
            dens(k-1)=sigocn(ttem(k-1),ssal(k-1))-thbase
          else
            dens(k-1)=dens(k-1)*(1.-q)+dens(k)*q
            ttem(k-1)=tofsig(dens(k-1)+thbase,ssal(k-1))
          end if
        end if
        pres(k)=p_hat
        nwrk3=nwrk3+1
      end if
c
c --- do we need to inflate layer k by lowering  l o w e r  interface?
c
 6    p_hat=pres(k)+dplist(2)*onem
      if (k.lt.kk .and. pres(k+1).lt.pres(kk+1)-onemm .and.
     .    pres(k+1).lt.p_hat) then
        p_hat=min(p_hat,pres(k)+
     .        max(onecm,slak*delt1*(p_hat-pres(k))))
        khyb=k				!  layer k is in hybrid domain
        go to 5
      end if
      go to 8
c
c --- water in layer k is too  l i g h t . dilute with water from layer k+1
c                              ^^^^^^^^^
 7    if (k.ge.kk .or. pres(k+1).gt.pres(kk+1)-onemm) go to 8
c
      q=(dens(k)-targt(k))/max(dens(k+1)-targt(k),sigjmp)
      p_hat=pres(k+1)*(1.-q)+pres(k)*q
c
c --- curtail downward growth of layers (esp. lowest hybrid layer)
      p_hat=max(pres(k+1),min(p_hat,.5*(pres(k)+pres(k+2))))
      p_hat=min(p_hat,pres(k+1)+
     .      max(onecm,slak*delt1*(p_hat-pres(k+1))))
c
 5    p_hat=min(p_hat,pres(k+2))
      if (p_hat.gt.pres(k+1)+onemm) then
c
        if (vrbos) write (lp,105)
     .   i,j,k,'lower intfc',pres(k+1)/onem,'=>',p_hat/onem
c
        if (k.le.khyb .and. k.lt.kk-1) then	!  use plm
          if (vrbos)
     .    write (lp,'(2i5,i3,a)') i,j,k,'  entrain from layer below'
          displ(1)=0.
          displ(2)=p_hat-pres(k+1)
          displ(3)=0.
          displ(4)=0.
          call plmad3(pres(k),displ,ssal(k),ssal(k),vrbos)
          if (tscnsv) then
            call plmad3(pres(k),displ,ttem(k),ttem(k),vrbos)
            dens(k  )=sigocn(ttem(k  ),ssal(k  ))-thbase
            dens(k+1)=sigocn(ttem(k+1),ssal(k+1))-thbase
          else
            call plmad3(pres(k),displ,dens(k),dens(k),vrbos)
            ttem(k  )=tofsig(dens(k  )+thbase,ssal(k  ))
            ttem(k+1)=tofsig(dens(k+1)+thbase,ssal(k+1))
          end if
        else					!  use pcm
          q=(p_hat-pres(k+1))/max(p_hat-pres(k),epsil)
          if (q.lt.0. .or. q.gt.1.) then
            write (lp,*) 'i,j,k,p_hat,pres(k+1),q=',
     .                    i,j,k,p_hat,pres(k+1),q
            abort=.true.
          end if
          ssal(k)=ssal(k)*(1.-q)+ssal(k+1)*q
          if (tscnsv) then
            ttem(k)=ttem(k)*(1.-q)+ttem(k+1)*q
            dens(k)=sigocn(ttem(k),ssal(k))-thbase
          else
            dens(k)=dens(k)*(1.-q)+dens(k+1)*q
            ttem(k)=tofsig(dens(k)+thbase,ssal(k))
          end if
        end if
        pres(k+1)=p_hat
        nwrk3=nwrk3+1
      end if
 8    continue
c
      tndcyr=-torho
      tndcyt=-totem
      tndcys=-tosal
      do k=1,kk
        tndcyr=tndcyr+dens(k)*(pres(k+1)-pres(k))
        tndcyt=tndcyt+ttem(k)*(pres(k+1)-pres(k))
        tndcys=tndcys+ssal(k)*(pres(k+1)-pres(k))
      end do
      if (tscnsv) then
        if (abs(tndcyt).gt.acurcy*10.*pres(kk+1))
     .   write (lp,104) i,j,'  hybgen - bad temp.intgl.:',totem,
     .    tndcyt,tndcyt/(10.*pres(kk+1))
      else
        if (abs(tndcyr).gt.acurcy*thbase*pres(kk+1))
     .   write (lp,104) i,j,'  hybgen - bad dens.intgl.:',torho,
     .    tndcyr,tndcyr/(thbase*pres(kk+1))
      end if
      if (abs(tndcys).gt.acurcy*35.*pres(kk+1))
     . write (lp,104) i,j,'  hybgen - bad saln.intgl.:',tosal,
     .  tndcys,tndcys/(35.*pres(kk+1))
 104  format (2i5,a,1p,2e15.7,e9.1)
c
      if (vrbos) then
        write (lp,99) nstep,i,j,'      n e w   p r o f i l e :'
        do k=1,kk,10
        write (lp,100) (pres(k1)/onem,k1=k,min(kk+1,k+10))
        write (lp,101) (sigref*(dens(k1)+thbase),k1=k,min(kk,k+9))
        write (lp,102) (ttem(k1),k1=k,min(kk,k+9))
        write (lp,102) (ssal(k1),k1=k,min(kk,k+9))
        end do
      end if
c
c --- put 1-d column back into 3-d grid
c
      do 2 k=1,kk
      kn=k+nn
      th3d(i,j,kn)=dens(k)
      temp(i,j,kn)=ttem(k)
      saln(i,j,kn)=ssal(k)
      p(i,j,k+1)=pres(k+1)
      dp(i,j,kn)=pres(k+1)-pres(k)
      diaflx(i,j,k)=diaflx(i,j,k)+(dp(i,j,kn)-dpold(i,j,k))	!  diapyc.flux
 2    continue
c
      ntot2d(j)=ntot2
      ntot3d(j)=ntot3
      nwrk2d(j)=nwrk2
      nwrk3d(j)=nwrk3
 12   continue
c$OMP END PARALLEL DO
      if (abort) stop '(error in hybgen -- q out of bounds)'
c
c$OMP PARALLEL DO
      do 1 j=1,jj
      do 1 k=1,kk
      do 1 l=1,isp(j)
      do 1 i=ifp(j,l),ilp(j,l)
 1    p(i,j,k+1)=p(i,j,k)+dpold(i,j,k)
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ja)
      do 88 j=1,jj
      ja=mod(j-2+jj,jj)+1
      do 88 k=2,kk+1
c
      do 881 l=1,isu(j)
      do 881 i=ifu(j,l),ilu(j,l)
 881  pu(i,j,k)=min(depthu(i,j),.5*(p(i,j,k)+p(i-1,j,k)))
c
      do 882 l=1,isv(j)
      do 882 i=ifv(j,l),ilv(j,l)
 882  pv(i,j,k)=min(depthv(i,j),.5*(p(i,j,k)+p(i,ja ,k)))
 88   continue
c$OMP END PARALLEL DO
c
      call dpudpv(nn)
c
c$OMP PARALLEL DO PRIVATE(kn,pold,pnew,displ,trac,uold,vold,vrbos)
      do 13 j=1,jj
c
      if (trcout) then
c
c --- evaluate effect of regridding on tracer field(s)
c
        do 20 l=1,isp(j)
        do 20 i=ifp(j,l),ilp(j,l)
c
        vrbos=.false.
cdiag   if (i.eq.itest .and. j.eq.jtest) vrbos=.true.
c
        pold(1)=p(i,j,1)
        pnew(1)=p(i,j,1)
c
        do 21 k=1,kk
        pold(k+1)=pold(k)+dpold(i,j,k)
        pnew(k+1)=p(i,j,k+1)
        displ(k+1)=pnew(k+1)-pold(k+1)
 21     trac(k)=tracer(i,j,k)
        displ(   1)=0.
        displ(kk+1)=0.
c
        call plmadv(kk,pold,displ,trac,trac,vrbos)
c
        do 20 k=1,kk
 20     tracer(i,j,k)=trac(k)
c
      end if				!  trcout
c
c --- evaluate effect of regridding on -u-
c
      do 14 l=1,isu(j)
      do 14 i=ifu(j,l),ilu(j,l)
c
      vrbos=.false.
cdiag if (i.eq.itest .and. j.eq.jtest) vrbos=.true.
c
      pold(1)=p(i,j,1)
      pnew(1)=p(i,j,1)
c
      do 15 k=1,kk
      kn=k+nn
      uold(k)=u(i,j,kn)
      pold(k+1)=pu(i,j,k+1)
      pnew(k+1)=pnew(k)+dpu(i,j,kn)
      pu(i,j,k+1)=pnew(k+1)
 15   displ(k+1)=pnew(k+1)-pold(k+1)
      displ(   1)=0.
      displ(kk+1)=0.
c
      call plmadv(kk,pold,displ,uold,uold,vrbos)
c
      do 14 k=1,kk
 14   u(i,j,k+nn)=uold(k)
c
c --- evaluate effect of regridding on -v-
c
      do 24 l=1,isv(j)
      do 24 i=ifv(j,l),ilv(j,l)
c
      vrbos=.false.
cdiag if (i.eq.itest .and. j.eq.jtest) vrbos=.true.
c
      pold(1)=p(i,j,1)
      pnew(1)=p(i,j,1)
c
      do 25 k=1,kk
      kn=k+nn
      vold(k)=v(i,j,kn)
      pold(k+1)=pv(i,j,k+1)
      pnew(k+1)=pnew(k)+dpv(i,j,kn)
      pv(i,j,k+1)=pnew(k+1)
 25   displ(k+1)=pnew(k+1)-pold(k+1)
      displ(   1)=0.
      displ(kk+1)=0.
c
      call plmadv(kk,pold,displ,vold,vold,vrbos)
c
      do 24 k=1,kk
 24   v(i,j,k+nn)=vold(k)
c
 13   continue
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO
      do 9 j=1,jj
      do 9 k=1,kk
      do 9 l=1,isp(j)
      do 9 i=ifp(j,l),ilp(j,l)
 9    p(i,j,k+1)=p(i,j,k)+dp(i,j,k+nn)
c$OMP END PARALLEL DO
c
      call dpudpv(nn)
c
cdiag if (itest.gt.0 .and. jtest.gt.0) then
cdiag   write (lp,103) nstep,itest,jtest,
cdiag.  '  exiting  hybgen:  temp    saln    dens    thkns    dpth',
cdiag.  (k,temp(itest,jtest,k+nn),saln(itest,jtest,k+nn),sigref*
cdiag.  (th3d(itest,jtest,k+nn)+thbase),dp(itest,jtest,k+nn)/onem,
cdiag.  p(itest,jtest,k+1)/onem,k=1,kk)
cdiag   write (lp,106) nstep,itest,jtest,
cdiag.  '  exiting  hybgen:  dpthu      u    dpthv      v',
cdiag.  (k,pu(itest,jtest,k+1)/onem,u(itest,jtest,k+nn),
cdiag.     pv(itest,jtest,k+1)/onem,v(itest,jtest,k+nn),k=1,kk)
cdiag end if
c
      if (mod(time+.0001,1.).lt..0002) then
        nwrk2=0
        nwrk3=0
        ntot2=0
        ntot3=0
        do j=1,jdm
          nwrk2=nwrk2+nwrk2d(j)
          nwrk3=nwrk3+nwrk3d(j)
          ntot2=ntot2+ntot2d(j)
          ntot3=ntot3+ntot3d(j)
        end do
        write (lp,'(a,f6.1,a,i9,a)') 'hybgen - grid restoration at',
     .   100.*float(nwrk3)/float(ntot3),' per cent of',ntot3,' points'
        write (lp,'(a,f6.1,a,i9,a)') 'hybgen - new bottom layer at',
     .   100.*float(nwrk2)/float(ntot2),' per cent of',ntot2,' points'
      end if
      return
      end
c
c
      subroutine plmad3(x,dx,y,ynew,diagno)
c
c --- advection by piecewise linear method
c --- this is a special version taylored to nmax = 3
c                                           ^^^^^^^^
c --- input variables:
c --- y(nmax)    - function values at cell midpoints
c --- x(nmax+1)  - cell boundaries
c --- dx(nmax+1) - displacement of cell boundaries during 1 time step
c
c --- output variables:
c --- ynew(nmax) - function values after advection (overwriting of -y- allowed)
c
      implicit none
      integer nmax,ndim,n
      parameter (nmax=3,ndim=3)
      real x(nmax+1),dx(nmax+1),y(nmax),ynew(nmax),total,tndcy,onemu
     .    ,ytmp(ndim),slope(ndim),yleft,yrigh,wgt,slab,dxnew,acurcy
      logical diagno
      data acurcy/1.e-11/,onemu/.098/
      if (ndim.lt.nmax) stop '(ndim too small in plmadv)'
c
      total=0.
      do 3 n=1,nmax
      if (x(n).gt.x(n+1)+onemu) then
        write (*,'(a,4f9.0)') 'error: pressure inversion in plmadv',x
        stop '(plmadv)'
      end if
c
      ytmp(n)=y(n)*(x(n+1)-x(n))
      ynew(n)=y(n)
 3    total=total+ytmp(n)
c
c --- get -y- slope at point -2-
c
      if (y(2).le.min(y(1),y(3)) .or.
     .    y(2).ge.max(y(1),y(3))) then
        slope(2)=0.
      else if ((y(3)-y(1))*(y(1)+y(3)-2.*y(2)).gt.0.) then
        slope(2)=y(2)-y(1)
      else
        slope(2)=y(3)-y(2)
      end if
c
c --- now transport -y- across cell boundaries
c
      yleft=y(2)-.5*slope(2)
      yrigh=y(2)+.5*slope(2)
      if (dx(2).gt.0.) then
c --- velocity at left edge is negative. export slab to cell on left
        wgt=.5*dx(2)/(x(3)-x(2))
        slab=(yrigh*wgt+yleft*(1.-wgt))*dx(2)
        ytmp(2)=ytmp(2)-slab
        ytmp(1)=ytmp(1)+slab
      end if
      if (dx(3).lt.0.) then
c --- velocity at right edge is positive. export slab to cell on right
        wgt=-.5*dx(3)/(x(3)-x(2))
        slab=-(yleft*wgt+yrigh*(1.-wgt))*dx(3)
        ytmp(2)=ytmp(2)-slab
        ytmp(3)=ytmp(3)+slab
      end if
c
      if (diagno) write (*,100) 'plmadv in: ',x,y
 100  format (a,4f9.0,3f9.4)
c
      tndcy=-total
      do 5 n=1,nmax
      dxnew=x(n+1)-x(n)+dx(n+1)-dx(n)
      if (dxnew.ne.0.) ynew(n)=ytmp(n)/dxnew
 5    tndcy=tndcy+ynew(n)*dxnew
      if (abs(tndcy).gt.acurcy*abs(total)) write (*,'(a,1p,2e11.3)')
     .  'plmadv - bad intgl.:',total,tndcy
c
      if (diagno) write (*,100) 'plmadv out:',dx,ynew
c
      return
      end
c
c
      subroutine plmadv(nmax,x,dx,y,ynew,diagno)
c
c --- advection by piecewise linear method
c --- this version is customized for recursive updating of cell boundaries
c                                    ^^^^^^^^^
c --- input variables:
c --- y(nmax)    - function values at cell midpoints
c --- x(nmax+1)  - cell boundaries
c --- dx(nmax+1) - displacement of cell boundaries during 1 time step
c
c --- output variables:
c --- ynew(nmax) - function values after advection (overwriting of -y- allowed)
c
      implicit none
      integer nmax,ndim,n,na,nb,m
      parameter (ndim=25)
      real x(nmax+1),dx(nmax+1),y(nmax),ynew(nmax),
     .     xtmp(ndim+1),ytmp(ndim),slope(ndim),yold(ndim),
     .     yleft,yrigh,wgt,slab,dxnew,acurcy,total,tndcy,onemu,uvscl
      logical diagno
      data acurcy/1.e-11/,onemu/.098/,uvscl/.02/
      if (ndim.lt.nmax) stop '(ndim too small in plmadv)'
c
      total=0.
      do 3 n=1,nmax
      if (x(n).gt.x(n+1)+onemu) then
        write (*,'(a/(8f9.0))') 'error: x not monotonic in plmadv:',x
        stop
      end if
c
      xtmp(n)=x(n)
      ytmp(n)=y(n)*(x(n+1)-x(n))
      yold(n)=y(n)
      ynew(n)=y(n)
 3    total=total+ytmp(n)
      xtmp(nmax+1)=x(nmax+1)
c
      do 4 n=1,nmax
      na=max(   1,n-1)			!  non-cyclic
      nb=min(nmax,n+1)			!  non-cyclic
      if (xtmp(n)+dx(n).lt.xtmp(na )-onemu .or.
     .    xtmp(n)+dx(n).gt.xtmp(n+1)+onemu) then
        write (*,'(a,i3/(i3,3f15.2))') 'error: x+dx out of bounds, n =',
     .   n,(m,x(m),dx(m),xtmp(m),m=1,n-1),
     .      n,x(n),dx(n),xtmp(n)+dx(n),n+1,x(n+1)
        stop
      end if
c
      if (xtmp(n+1).gt.xtmp(n)) then
c
c --- get -y- slope at point -n-
c
        if (y(n).le.min(y(na),y(nb)) .or.
     .      y(n).ge.max(y(na),y(nb))) then
          slope(n)=0.
        else if ((y(nb)-y(na))*(y(na)+y(nb)-2.*y(n)).gt.0.) then
          slope(n)=y(n)-y(na)
        else
          slope(n)=y(nb)-y(n)
        end if
c
c --- now transport -y- across cell boundaries
c
        yleft=y(n)-.5*slope(n)
        yrigh=y(n)+.5*slope(n)
        if (dx(n).gt.0.) then
c --- velocity at left edge is negative. export slab to cell on left
          wgt=.5*dx(n)/(xtmp(n+1)-xtmp(n))
          slab=(yrigh*wgt+yleft*(1.-wgt))*dx(n)
          ytmp(n )=ytmp(n )-slab
          ytmp(na)=ytmp(na)+slab
        end if
        if (dx(n+1).lt.0.) then
c --- velocity at right edge is positive. export slab to cell on right
          wgt=-.5*dx(n+1)/(xtmp(n+1)-xtmp(n))
          slab=-(yleft*wgt+yrigh*(1.-wgt))*dx(n+1)
          ytmp(n )=ytmp(n )-slab
          ytmp(nb)=ytmp(nb)+slab
        end if
      end if
c
c --- before proceding to next cell, update cell boundary
 4    xtmp(n)=xtmp(n)+dx(n)
c
      tndcy=-total
      do 5 n=1,nmax
      dxnew=xtmp(n+1)-xtmp(n)
      if (dxnew.ne.0.) ynew(n)=ytmp(n)/dxnew
 5    tndcy=tndcy+ynew(n)*dxnew
      if (abs(tndcy).gt.acurcy*max(abs(total),uvscl*x(nmax+1))) then
        write (*,'(a,1p,2e11.3)') 'plmadv - bad intgl.:',total,tndcy
        write (*,100) 'plmadv:',
     .   (n,x(n),yold(n),dx(n),xtmp(n),ynew(n),n=1,nmax)
      end if
c
      if (diagno) write (*,100) 'plmadv:',
     .   (n,x(n),yold(n),dx(n),xtmp(n),ynew(n),n=1,nmax)
 100  format (a/(i3,f14.1,f14.7,2f14.1,f14.7))
c
      return
      end
c
c
      real function cushn(delp,dp0)
c
c --- c u s h i o n   function (from Bleck & Benjamin, 1993):
c
c                 (x + x1 - 2)**2
c --- cushn = 1 + ---------------    where  x = delp/dp0
c                  4 ( x1 - 1)
c
c --- if delp >>  0, -cushn- returns -delp-
c --- if delp <<  0, -cushn- returns -dp0-
c
      real delp,dp0,qq,x1,factor
ccc   parameter (x1=4.)			!  used in Bleck&Benjamin 1993
ccc   parameter (x1=6.)			!  used in Bleck 2001
      parameter (x1=8.)
c
      parameter (factor=.25/(x1-1.))
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc   qq=max(-1.,min(2.,delp/(2.*dp0)))
ccc   cushn=dp0*(1.+athird*(qq+1.)**2)
ccc  .            *max(1.,delp/(2.*dp0))
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc   qq=max(-.75,min(1.25,delp/(4.*dp0)))
ccc   cushn=dp0*(1.+(qq+.75)**2)
ccc  .            *max(1.,delp/(5.*dp0))
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc   qq=max(-1.,min(1.5,delp/(4.*dp0)))
ccc   cushn=dp0*(1+.8*(qq+1.)**2)
ccc  .            *max(1.,delp/(6.*dp0))
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      qq=max(2.-x1,min(x1,delp/dp0))
      cushn=dp0*(1.+factor*(qq+x1-2.)**2) * max(1.,delp/(x1*dp0))
      return
      end
c
c
c> Revision history:
c>
c> May  2001 - added u/v regridding
c> June 2001 - added interlayer ("diapycnal") mass flux diagnostics
c> Aug. 2001 - added -kn- to list of private variables in loop 13
c> Aug. 2001 - various refinements, incl. multi-layer ingest
c> Dec. 2001 - softened enforcement of minimum layer thickness constraint
c> Feb. 2002 - restricted ingest of multiple layers from below
c> Mar. 2002 - added passive tracer
c> Apr. 2002 - fixed -diaflx- diagnostics by defining -dpold-
c> June 2002 - changed 'slak' logic - now based on gradual restoring
c> June 2002 - curtail downward growth of layers (esp. lowest hybrid layer)
c> May  2003 - introduced variable -dpsum- to make -p_hat- k-dependent
c> Aug. 2003 - added option to maintain 10 cm min.thickness throughout column
c> Sep. 2003 - added logical switch to enable/disable tracer redistribution
c> Sep. 2003 - provided choice between T/S and rho/S conservation
c> Nov. 2003 - replaced dp0 in cushn call by dp0 appropriate for layer k-1
c> Nov. 2003 - accelerated relaxation time for 1st interface (slak x 10)
c> Nov. 2004 - allowed -dplist- values to shrink near equator
c> Nov. 2004 - conversion to piecewise linear advection scheme (PLM)
c> Nov. 2004 - extended PLM advection to velocity field
