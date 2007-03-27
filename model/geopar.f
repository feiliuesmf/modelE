      subroutine geopar
c
c --- set up model parameters related to geography
c
c --- hycom version 0.9 -- cyclic in j
css   USE GEOM, only : dxyp
c
      implicit none
      include 'dimensions.h'
      include 'dimension2.h'
      include 'common_blocks.h'
      include 'cpl.h'
      include 'a2o.h'
      include 'kprf_arrays.h'
c
      real realat,sphdis,glufac,zero
      integer idim,jdim,length,iz,jz,nt
      character util(idm*jdm+14)*2,preambl(5)*79
      real*4 real4(idm,jdm),lat4(idm,jdm,4),lon4(idm,jdm,4)
c --- 'glufac' = regional viscosity enhancement factor
      data glufac/3./
      data zero/0./
c
c --- read basin depth array
      write (lp,'(2a)') ' reading bathymetry file from ',flnmdep
      open (unit=9,file=flnmdep,form='unformatted',status='old'
     .     ,action='read')
      read (9) iz,jz
      if (iz.ne.idm .or. jz.ne.jdm) then
        write (lp,'(2(a,2i5))') 'depth file dimensions',iz,jz,
     .   '  should be',idm,jdm
        stop '(geopar)'
      end if
      rewind (9)
      read (9) iz,jz,((real4(i,j),i=1,iz),j=1,jz)
      close (unit=9)
c$OMP PARALLEL DO
      do 9 j=1,jj
      do 9 i=1,ii
 9    depths(i,j)=real4(i,j)
c$OMP END PARALLEL DO
c
cccc$OMP PARALLEL DO
c     do 7 j=1,jj
c     do 7 i=1,ii
c7    if (depths(i,j).gt.0.) depths(i,j)=max(botmin.,depths(i,j))
cccc$OMP END PARALLEL DO
c
c --- reset the Denmark Strait - done in advance
c
c     depths(69,168)=798.
c     depths(69,169)=798.
c     depths(70,167)=798.
c     depths(70,168)=798.
c     depths(71,167)=798.
c
c --- reset the Iceland-Faeroes ridge - done in advance
c     do i=72,74
c     depths(i,175)=798.
c     end do
c
      write (lp,*) 'shown below: bottom depth'
      call zebra(depths,idm,ii1,jj)
c
c --- determine do-loop limits for u,v,p,q points
      call bigrid(depths)
ccc      do 3 i=1,ii1
ccc 3    write (lp,'('' i='',i3,'' jfp,jlp='',7(1x,2i5))') i,
ccc     . (jfp(i,l),jlp(i,l),l=1,jsp(i))
ccc      do 5 j=1,jj
ccc 5    write (lp,'('' j='',i3,'' ifp,ilp='',7(1x,2i5))') j,
ccc     . (ifp(j,l),ilp(j,l),l=1,isp(j))
c
c --- smooth bottom topography (optional)
ccc      call psmoo(depths)
c
c     call prtmsk(ip,depths,util1,idm,ii1,jj,0.,1.,
c    .     'bottom depth (m)')
c
      write (lp,'(2a)') 'read lat/lon from ',flnmlat
      open (33,file=flnmlat,form='unformatted',status='old')
      read (33) iz,jz
      if (iz.ne.idm .or. jz.ne.jdm) then
        write (lp,'(2(a,2i5))') 'error - idm,jdm =',iz,jz,
     .   ' in lat/lon file should be',idm,jdm
        stop '(geopar)'
      end if
      rewind 33
      read (33) iz,jz,lat4,lon4
      close(33)
c
c$OMP PARALLEL DO PRIVATE(n)
      do 8 j=1,jj
      do 8 n=1,4
      do 8 i=1,ii
      latij(i,j,n)=lat4(i,j,n)
 8    lonij(i,j,n)=lon4(i,j,n)
c$OMP end PARALLEL DO
c
      write (lp,*) 'shown below: latitude of vorticity points'
      call zebra(latij(1,1,4),idm,ii,jj)
      write (lp,*) 'shown below: longitude of vorticity points'
      call zebra(lonij(1,1,4),idm,ii,jj)
c
c --- define coriolis parameter and grid size
c$OMP PARALLEL DO PRIVATE(ja,jb)
      do 56 j=1,jj
      ja=mod(j-2+jj,jj)+1
      jb=mod(j     ,jj)+1
      do 56 i=1,ii
c
      corio(i,j)=sind(latij(i,j,4))*4.*pi/86164.        !  86400 * 365 / 366
c
      scpy(i,j)=sphdis(latij(i,j ,2),lonij(i,j ,2),
     .                 latij(i,jb,2),lonij(i,jb,2))
c
      if (i.lt.ii) then
      scpx(i,j)=sphdis(latij(i  ,j,1),lonij(i  ,j,1),
     .                 latij(i+1,j,1),lonij(i+1,j,1))
      scp2(i,j)=scpx(i,j)*scpy(i,j)
      scp2i(i,j)=1./scp2(i,j)
      end if
c
      scuy(i,j)=sphdis(latij(i,j ,4),lonij(i,j ,4),
     .                 latij(i,jb,4),lonij(i,jb,4))
c
      if (i.gt.1) then
      scux(i,j)=sphdis(latij(i  ,j,3),lonij(i  ,j,3),
     .                 latij(i-1,j,3),lonij(i-1,j,3))
      scu2(i,j)=scux(i,j)*scuy(i,j)
      scuxi(i,j)=1./scux(i,j)
      end if
c
      scvy(i,j)=sphdis(latij(i,j ,3),lonij(i,j ,3),
     .                 latij(i,ja,3),lonij(i,ja,3))
c
      if (i.lt.ii) then
      scvx(i,j)=sphdis(latij(i  ,j,4),lonij(i  ,j,4),
     .                 latij(i+1,j,4),lonij(i+1,j,4))
      scv2(i,j)=scvx(i,j)*scvy(i,j)
      scvyi(i,j)=1./scvy(i,j)
      end if
c
      scqy(i,j)=sphdis(latij(i,j ,1),lonij(i,j ,1),
     .                 latij(i,ja,1),lonij(i,ja,1))
c
      if (i.gt.1) then
      scqx(i,j)=sphdis(latij(i  ,j,2),lonij(i  ,j,2),
     .                 latij(i-1,j,2),lonij(i-1,j,2))
      scq2(i,j)=scqx(i,j)*scqy(i,j)
      scq2i(i,j)=1./scq2(i,j)
      end if
c
      if (i.eq.71.and.j.eq.91) write(*,*)i,j,' scp=',scp2(i,j)
      if (i.eq.65.and.j.eq.136) write(*,*)i,j,' scp=',scp2(i,j)
 56   continue
c$OMP END PARALLEL DO
c
      if (beropn .and. scu2(ipacs,jpac).ne.scu2(iatln,jatl))
     .  write(*,'(a,6f13.5)') ' chk WRONG scu2'
     . ,scu2(ipacs,jpac),scu2(iatln,jatl)
c
      if ( latij(ipacs,jpac,3).ne.latij(iatls,jatl,3)
     . .or.latij(ipacn,jpac,3).ne.latij(iatln,jatl,3)
     . .or.lonij(ipacs,jpac,3).ne.lonij(iatls,jatl,3)
     . .or.lonij(ipacn,jpac,3).ne.lonij(iatln,jatl,3))
     .  write(*,'(a,8f9.2)') ' chk WRONG lat/lon '
     .,latij(ipacs,jpac,3),latij(iatls,jatl,3)
     .,latij(ipacn,jpac,3),latij(iatln,jatl,3)
     .,lonij(ipacs,jpac,3),lonij(iatls,jatl,3)
     .,lonij(ipacn,jpac,3),lonij(iatln,jatl,3)
c
      write(*,'(a,2f8.2)') 'lat/lon of Bering Strait:'
     .   ,latij(ipacs,jpac,3),lonij(ipacs,jpac,3)
c     write (lp,'('' shown below: coriolis parameter'')')
c     call zebra(corio,idm,ii,jj)
c     write (lp,'('' shown below: grid cell size'')')
c     call zebra(scp2,idm,ii,jj)
c
      area=0.
      avgbot=0.
c
      do 57 j=1,jj
      do 57 l=1,isp(j)
      do 57 i=ifp(j,l),ilp(j,l)
      avgbot=avgbot+depths(i,j)*scp2(i,j)
 57   area=area+scp2(i,j)
      avgbot=avgbot/area
      write (lp,100) avgbot,area
 100  format(' mean basin depth (m) and area (10^6 km^2):',f9.1,
     .       -12p,f9.1)
c
c --- initialize some arrays
c
css   if (nstep0.eq.0) then
      write (lp,*) 'laying out arrays in memory ...'
c$OMP PARALLEL DO
      do 209 j=1,jj
      do 209 i=1,ii
      p(i,j,1)=huge
      if (ip(i,j).eq.1) p(i,j,1)=zero
      pu(i,j,1)=huge
      pv(i,j,1)=huge
      pbot(i,j)=huge
      ubavg(i,j,1)=huge
      ubavg(i,j,2)=huge
      ubavg(i,j,3)=huge
      vbavg(i,j,1)=huge
      vbavg(i,j,2)=huge
      vbavg(i,j,3)=huge
      utotm(i,j)=huge
      vtotm(i,j)=huge
      utotn(i,j)=huge
      vtotn(i,j)=huge
      uflux (i,j)=huge
      vflux (i,j)=huge
      uflux1(i,j)=huge
      vflux1(i,j)=huge
      uflux2(i,j)=huge
      vflux2(i,j)=huge
      uflux3(i,j)=huge
      vflux3(i,j)=huge
      uja(i,j)=huge
      ujb(i,j)=huge
      via(i,j)=huge
      vib(i,j)=huge
      pgfx(i,j)=huge
      pgfy(i,j)=huge
      depthu(i,j)=huge
      depthv(i,j)=huge
      tprime(i,j)=huge
c
      srfhgt(i,j)=zero
      dpmixl(i,j)=zero
      oice(i,j)=zero
      taux(i,j)=zero
      tauy(i,j)=zero
      oflxa2o(i,j)=zero
      osalt(i,j)=zero
      oemnp(i,j)=zero
      ustar(i,j)=zero
      sswflx(i,j)=zero
c
      ubavav(i,j)=zero
      vbavav(i,j)=zero
      pbavav(i,j)=zero
      sfhtav(i,j)=zero
      dpmxav(i,j)=zero
      oiceav(i,j)=zero
      eminpav(i,j)=zero
      surflav(i,j)=zero
      sflxav(i,j)=zero
      brineav(i,j)=zero
c
      do 209 k=1,kk
      u  (i,j,k   )=huge
      u  (i,j,k+kk)=huge
      v  (i,j,k   )=huge
      v  (i,j,k+kk)=huge
      uflx(i,j,k)=huge
      vflx(i,j,k)=huge
      ufxcum(i,j,k)=huge
      vfxcum(i,j,k)=huge
      dpinit(i,j,k)=huge
      dpold (i,j,k)=huge
      dp (i,j,k   )=huge
      dp (i,j,k+kk)=huge
      dpu(i,j,k   )=huge
      dpu(i,j,k+kk)=huge
      dpv(i,j,k   )=huge
      dpv(i,j,k+kk)=huge
      p (i,j,k+1)=huge
      pu(i,j,k+1)=huge
      pv(i,j,k+1)=huge
c
      th3d(i,j,k+kk)=huge
      th3d(i,j,k   )=huge
      thermb(i,j,k+kk)=huge
      thermb(i,j,k   )=huge
      thstar(i,j,k)=huge
      do nt=1,ntrcr
        tracer(i,j,k,nt)=zero
      end do
      uav(i,j,k)=zero
      vav(i,j,k)=zero
      dpuav(i,j,k)=zero
      dpvav(i,j,k)=zero
      dpav (i,j,k)=zero
      temav(i,j,k)=zero
      salav(i,j,k)=zero
      th3av(i,j,k)=zero
      uflxav(i,j,k)=zero
      vflxav(i,j,k)=zero
      diaflx(i,j,k)=zero
 209  continue
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ja)
      do 210 j=1,jj
      ja=mod(j-2+jj,jj)+1
      do 210 l=1,isq(j)
      do 210 i=ifq(j,l),ilq(j,l)
      pbot(i  ,j  )=0.
      pbot(i-1,j  )=0.
      pbot(i  ,ja )=0.
      pbot(i-1,ja )=0.
      p(i  ,j  ,1)=0.
      p(i-1,j  ,1)=0.
      p(i  ,ja ,1)=0.
      p(i-1,ja ,1)=0.
      do 210 k=1,kk
      dp(i  ,j  ,k   )=0.
      dp(i  ,j  ,k+kk)=0.
      dp(i-1,j  ,k   )=0.
      dp(i-1,j  ,k+kk)=0.
      dp(i  ,ja ,k   )=0.
      dp(i  ,ja ,k+kk)=0.
      dp(i-1,ja ,k   )=0.
 210  dp(i-1,ja ,k+kk)=0.
c$OMP END PARALLEL DO
c
c --- initialize  u,ubavg,utotm,uflx,uflux,uflux2/3,uja,ujb  at points
c --- located upstream and downstream (in i direction) of p points.
c --- initialize  depthu,dpu,utotn,pgfx  upstream and downstream of p points
c --- as well as at lateral neighbors of interior u points.
c
c$OMP PARALLEL DO PRIVATE(ja,jb)
      do 156 j=1,jj
      ja=mod(j-2+jj,jj)+1
      jb=mod(j     ,jj)+1
      do 156 l=1,isu(j)
      do 156 i=ifu(j,l),ilu(j,l)
      pu(i,j,1)=0.
c
      depthu(i,ja)=0.
      utotn (i,ja)=0.
      pgfx  (i,ja)=0.
c
      depthu(i,jb)=0.
      utotn (i,jb)=0.
      pgfx  (i,jb)=0.
c
      do 156 k=1,kk
      dpu(i,ja,k   )=0.
      dpu(i,ja,k+kk)=0.
c
      dpu(i,jb,k   )=0.
      dpu(i,jb,k+kk)=0.
 156  continue
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO
      do 158 j=1,jj
      do 158 l=1,isp(j)
      do 158 i=ifp(j,l),ilp(j,l)+1
      depthu(i,j)=0.
      utotn (i,j)=0.
      pgfx  (i,j)=0.
      ubavg(i,j,1)=0.
      ubavg(i,j,2)=0.
      ubavg(i,j,3)=0.
      utotm (i,j)=0.
      uflux (i,j)=0.
      uflux2(i,j)=0.
      uflux3(i,j)=0.
      uja(i,j)=0.
      ujb(i,j)=0.
c
      do 158 k=1,kk
      dpu(i,j,k   )=0.
      dpu(i,j,k+kk)=0.
      uflx(i,j,k)=0.
      ufxcum(i,j,k)=0.
      u(i,j,k   )=0.
 158  u(i,j,k+kk)=0.
c$OMP END PARALLEL DO
c
c --- initialize  v,vbavg,vtotm,vflx,vflux,vflux2/3,via,vib  at points
c --- located upstream and downstream (in j direction) of p points.
c --- initialize  depthv,dpv,vtotn,pgfy  upstream and downstream of p points
c --- as well as at lateral neighbors of interior v points.
c
      do 166 i=1,ii1
      ia=mod(i-2+ii,ii)+1
      ib=i+1
      do 166 l=1,jsv(i)
      do 166 j=jfv(i,l),jlv(i,l)
      pv(i,j,1)=0.
c
      depthv(ia,j)=0.
      vtotn (ia,j)=0.
      pgfy  (ia,j)=0.
c
      depthv(ib,j)=0.
      vtotn (ib,j)=0.
      pgfy  (ib,j)=0.
c
      do 166 k=1,kk
      dpv(ia,j,k   )=0.
      dpv(ia,j,k+kk)=0.
c
      dpv(ib,j,k   )=0.
      dpv(ib,j,k+kk)=0.
 166  continue
c
      do 168 i=1,ii1
      do 168 l=1,jsp(i)
      do 168 jp=jfp(i,l),jlp(i,l)+1
      j=mod(jp-1+jj,jj)+1
      depthv(i,j)=0.
      vtotn (i,j)=0.
      pgfy  (i,j)=0.
      vbavg(i,j,1)=0.
      vbavg(i,j,2)=0.
      vbavg(i,j,3)=0.
      vtotm (i,j)=0.
      vflux (i,j)=0.
      vflux2(i,j)=0.
      vflux3(i,j)=0.
      via(i,j)=0.
      vib(i,j)=0.
c
      do 168 k=1,kk
      dpv(i,j,k   )=0.
      dpv(i,j,k+kk)=0.
      vflx(i,j,k)=0.
      vfxcum(i,j,k)=0.
      v(i,j,k   )=0.
 168  v(i,j,k+kk)=0.
      write (lp,*) '... array layout completed'
css   endif                    ! end of nstep=0
c
c --- set 'glue' to values > 1 in regions where extra viscosity is needed
c
c$OMP PARALLEL DO PRIVATE(ia,ib,ja,jb)
      do 154 j=1,jj
      ja=mod(j-2+jj,jj)+1
      jb=mod(j     ,jj)+1
      do 154 i=1,ii
      ia=mod(i-2+ii,ii)+1
      ib=mod(i     ,ii)+1
      glue(i,j)=1.
c --- add glue in coastal areas
      if (depths(i,j).lt.200.
     .    .or. depths(ia,j).lt.200. .or. depths(i,ja).lt.200.
     .    .or. depths(ib,j).lt.200. .or. depths(i,jb).lt.200.)
     .    glue(i,j)=max(glue(i,j),glufac)
      if (depths(i,j).lt.100.
     .    .or. depths(ia,j).lt.100. .or. depths(i,ja).lt.100.
     .    .or. depths(ib,j).lt.100. .or. depths(i,jb).lt.100.)
     .    glue(i,j)=max(glue(i,j),2.*glufac)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- add glue to mediterranean:
      if (jdm.eq.256) then              !  1.4 deg, equator at i=64:
        if ((i.ge.  29 .and. i.le.  41 .and. j.le.  26)
     . .or. (i .ge. 36 .and. i.le.  37 .and. j.ge. 253))
     .        glue(i,j)=glufac
c
      else if (jdm.eq.180) then         !  2.0 deg, equator at i=115:
        if (i.ge.  91 .and. i.le.  98 .and. j.le.  25)
     .        glue(i,j)=glufac
      else
        write (lp,*) 'unable to determine location of Medterranean'
        stop '(geopar)'
      end if
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 154  continue
c$OMP END PARALLEL DO
c
c --- read in all weights
c
c$OMP PARALLEL DO
      do 11 ja=1,jjo
      do 11 ia=1,iio
11    ocellsz(ia,ja)=scp2(ia,ja)
c$OMP END PARALLEL DO
c
c     open(301,file='agcmgdsz.8bin',status='unknown',form='unformatted')
c     write(301)dxyp
c     close(301)
c
      if (iio*jjo*((nwgta2o*2 +1)*4+nwgta2o *8).ne.10249200 .or.
     .    iio*jjo*((nwgta2o2*2+1)*4+nwgta2o2*8).ne.10249200 .or.
     .    iia*jja*((nwgto2a*2 +1)*4+nwgto2a *8).ne.2079936 ) then
        print *,' wrong coupler=',iio*jjo*((nwgta2o*2+1)*4+nwgta2o*8)
     .  ,iio*jjo*((nwgta2o2*2+1)*4+nwgta2o2*8)
     .  ,iia*jja*((nwgto2a*2 +1)*4+nwgto2a *8)
        stop ' wrong size in flxa2o'
      endif
c
      open(14,file=flnma2o,form='unformatted',status='old',
     .  access='direct',recl=iio*jjo*((nwgta2o*2+1)*4+nwgta2o*8))
      read(14,rec=1) ilista2o,jlista2o,wlista2o,nlista2o
      close(14)
c
      open(15,file=flnma2o_tau,form='unformatted',status='old',
     .  access='direct',recl=iio*jjo*((nwgta2o2*2+1)*4+nwgta2o2*8))
      read(15,rec=1) itaua2o,jtaua2o,wtaua2o,ntaua2o
      close(15)
c
      open(16,file=flnmo2a,form='unformatted',status='old',
     .  access='direct',recl=iia*jja*((nwgto2a*2+1)*4+nwgto2a*8))
      read(16,rec=1) ilisto2a,jlisto2a,wlisto2a,nlisto2a
      close(16)
c
      open(17,file=flnmo2a_e,form='unformatted',status='old',
     .  access='direct',recl=iia*jja*((nwgto2a*2+1)*4+nwgto2a*8))
      read(17,rec=1) ilisto2a_e,jlisto2a_e,wlisto2a_e,nlisto2a_e
      close(17)
c
      open(18,file=flnmo2a_n,form='unformatted',status='old',
     .  access='direct',recl=iia*jja*((nwgto2a*2+1)*4+nwgto2a*8))
      read(18,rec=1) ilisto2a_n,jlisto2a_n,wlisto2a_n,nlisto2a_n
      close(18)
c
      open(22,file=flnmcoso,form='unformatted',status='old')
      read(22) iz,jz,coso,sino
      close(22)
      if (iz.ne.idm .or. jz.ne.jdm) then
        print *,' iz,jz=',iz,jz
        stop '(wrong iz/jz in cososino.8bin)'
      endif
c
c --- 1:8: NAT, SAT, NIN, SIN, NPA, SPA, ARC, SO, MED
      open (34,file=flnmbas,form='formatted',status='old')
      do n=1,2
      read(34,*)
      read(34,'(90i1)') ((msk(i,j),j=(n-1)*jj/2+1,n*jj/2),i=1,ii)
      enddo 
      close(34)
c
      return
      end
c
      function sphdis(x1,y1,x2,y2)
c --- dist.(m) between 2 points on sphere, lat/lon (x1,y1) and lat/lon (x2,y2)
      implicit none
      real x1,y1,x2,y2,sphdis,ang,radius,radian
      data radius/6375.e3/,radian/57.2957795/
c
      ang=mod(y2-y1+540.,360.)-180.
      sphdis=radius*acos(min(1.,cosd(90.-x1)*cosd(90.-x2)
     .                         +sind(90.-x1)*sind(90.-x2)*cosd(ang)))
      if (sphdis.eq.0.) 
     .  sphdis=radius*sqrt((x2-x1)**2+(ang*cosd(.5*(x1+x2)))**2)/radian
cdiag if (sphdis.eq.0.) write (*,'(a,2f8.3,2x,2f8.3)')
cdiag.  'warning - zero distance between lat/lon points',x1,y1,x2,y2
      sphdis=max(sphdis,1.)
      return
      end
c
c> Revision history
c>
c> Mar. 2000 - conversion to SI units
c> May  2000 - changed loop 56 from i/j to j/i to improve memory layout
c> Oct. 2000 - added code to compute 'glue'
c> Apr. 2001 - eliminated stmt_funcs.h
c> Dec. 2001 - added clause to assure p=0 in single-width channels (loop 209)
