      subroutine overtn(mm)
c
c --- hycom version 0.9
c
c --- diagnose meridional overturning (in indiv. basins) and heat/salt flux.
c --- n o t e : this version reapportions mass fluxes using  r e f l u x
c
      USE MODEL_COM, only :
     *  itime,iyear1,nday,jdendofm,jyear,jmon,jday,jdate,jhour,aMON
     * ,xlabel,lrunid
c
      implicit none
#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"
c
      integer nbasin
      parameter (nbasin=3)                !  3 basins (Atl.,Ind.,Pac.)
      character text*48,flnm*60,ocean(nbasin+1)*8
      real flux(kdm,idm,nbasin+1),heatfl(idm,nbasin+1),sunda(kdm+1),
     .     uflxnw(idm,jdm,kdm),vflxnw(idm,jdm,kdm),signw(idm,jdm,kdm),
     .     pnw(idm,jdm,kdm+1),work(kdm,idm),scale,thru,reviof,x,x1,x3,
     .     thrufl,suheat(kdm+1),
     .     saltfl(idm,nbasin+1),susalt(kdm+1)
      integer iub,num,imin,imax,indoi,indoj1,indoj2,nio,nb
     .       ,no,imx,imn,it,io,keep,iof,idrk1,jdrk1,idrk2,jdrk2
     .       ,imno,imxo,ioa,iob,lgth,maskk(kdm,idm)
c
      common /basins/ iub(nbasin+1,idm,jdm),num(idm,nbasin+1),
     .       imin(nbasin+1),imax(nbasin+1)
      save   /basins/
c
      data ocean/'ATLANTIC','  INDIAN',' PACIFIC','  GLOBAL'/
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- specify segment of -u- points representing Indonesian throughflow
c     data indoi/118/,indoj1,indoj2/54,70/                !  2.0 deg (ieq=115)
c     data idrk1,jdrk1/148,147/,idrk2,jdrk2/158,149/        !  2.0 deg (ieq=115)
      data indoi/129/,indoj1,indoj2/53,70/                !  2.0 deg (ieq=122)
      data idrk1,jdrk1/162,147/,idrk2,jdrk2/172,149/        !  2.0 deg (ieq=122)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      iof(i)=2*i-2
      reviof(x)=.5*(x+2.)
c
c --- read basin masks from disk
      nio=12
      write (lp,'(a/9x,a)') 'looking for basinmask data in',flnmbas
      open (unit=nio,file=flnmbas,status='old',form='unformatted',
     .      err=37)
c
c --- basinmask file exists - no need to recreate it
      write (lp,*) '..... file found'
      read (nio) iub,imin,imax,num
      close (unit=nio)
      go to 38
c
 37   write (lp,*) '..... file not found. create it.'
      call bsnmsk(nbasin,ocean,indoi,indoj1,indoj2)
c
c --- save basin masks on disk
      open (unit=nio,file=flnmbas,status='unknown',form='unformatted')
      write (nio) iub,imin,imax,num
      close (unit=nio)
c
 38   continue                        
c
      no=19
      do 14 lgth=60,1,-1
      if (flnmovt(lgth:lgth).eq.'/') go to 28
 14   continue
      write (lp,*) 'overtn  --  cannot find slash in',flnmovt
      stop
c
      call getdte(Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)

 28   write(flnm,'(a3,i4.4,2a)') amon,Jyear,'.ovtn',xlabel(1:lrunid)
      write (text,'(i11.11)') int(time+.001E+00)
      text(1:5)='ovtn.'
      open (unit=no,file=flnm,status='unknown')
      write (no,'(a11,'' (reflux) - + - + - + - + - + - +'')') text
c
c --- diagnose indonesian throughflow
c
      sunda(kk+1)=0.
      suheat(kk+1)=0.
      susalt(kk+1)=0.
      do 26 k=1,kk
      km=k+mm
      sunda(k)=0.
      susalt(k)=0.
      suheat(k)=0.
      do 35 j=indoj1,indoj2
      sunda(k)=sunda(k)+uflx(indoi,j,k)
      susalt(k)=susalt(k)+uflx(indoi,j,k)
     .                  *(saln(indoi,j,km)+saln(indoi-1,j,km))
 35   suheat(k)=suheat(k)+uflx(indoi,j,k)
     .                  *(temp(indoi,j,km)+temp(indoi-1,j,km))
      sunda(kk+1)=sunda(kk+1)+sunda(k)
      susalt(kk+1)=susalt(kk+1)+susalt(k)
 26   suheat(kk+1)=suheat(kk+1)+suheat(k)
c
      do 30 k=1,kk+1
      sunda(k)=sunda(k)/onem * 1.e-6
      susalt(k)=.5*susalt(k)/(-35.*onem) * 1.e-6        ! S-ward > 0
 30   suheat(k)=.5*suheat(k)*spcifh/g * 1.e-15                ! S-ward > 0
c
c --- transform hybrid mass fluxes to isopycnal fluxes
c
      call reflux(uflx  ,vflx  ,th3d(1,1,1+mm),p,
     .            uflxnw,vflxnw,signw,pnw,theta,kdm,kdm)
c
      do 15 nb=1,nbasin+1
c
c --- integrate meridional mass fluxes in zonal direction
      do 1 k=1,kk
      do 151 i=1,ii
 151  flux(k,i,nb)=0.
      do 1 j=1,jj
      do 1 i=imin(nb),imax(nb)
      if (iub(nb,i,j).eq.1) then
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- optional: exclude mediterranean
        if (jdm.eq.180) then                !  2.0 deg, equator at i=115
          if ((i.ge.  91 .and. i.le.  98 .and. j.le.  18)
     .   .or. (i .ge. 95 .and. i.le.  96 .and. j .ge.178)) go to 1
        else if (jdm.eq.256) then        !  1.4 deg, equator at i=64
          if ((i.ge.  29 .and. i.le.  41 .and. j.le.  26)
     .   .or. (i .ge. 36 .and. i.le.  37 .and. j .ge.253)) go to 1
        else
          write (lp,'(a)') 'j-dimension error in overtn.f'
          stop '(jdm)'
        end if
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        flux(k,i,nb)=flux(k,i,nb)-uflxnw(i,j,k)/onem * 1.e-6        ! => Sv
      end if
 1    continue
c
c --- integrate meridional heat/salt fluxes vertically and in zonal direction
      do 111 i=1,ii
      saltfl(i,nb)=0.
 111  heatfl(i,nb)=0.
c
      do 112 i=imin(nb),imax(nb)
      do 11 k=1,kk
      km=k+mm
      do 11 j=1,jj
      if (iub(nb,i,j).eq.1) then
        saltfl(i,nb)=saltfl(i,nb)+uflx(i,j,k)
     .                         *(saln(i,j,km)+saln(i-1,j,km))
        heatfl(i,nb)=heatfl(i,nb)+uflx(i,j,k)
     .                         *(temp(i,j,km)+temp(i-1,j,km))
      end if
 11   continue
      saltfl(i,nb)=-.5*saltfl(i,nb)/(-35.*onem) * 1.e-6                ! N-ward > 0
 112  heatfl(i,nb)=-.5*heatfl(i,nb)*spcifh/g * 1.e-15                ! N-ward > 0
c
c --- add mass fluxes vertically to produce stream function values
      do 2 k=1,kk
      do 221 i=1,ii
 221  if (k.gt.1) flux(k,i,nb)=flux(k,i,nb)+flux(k-1,i,nb)
c
c --- subtract out portion due to indonesian throughflow
      if (nbasin.gt.2) then
       if (nb.eq.2) then
        do 29 i=indoi,imax(nb)
 29     flux(k,i,nb)=flux(k,i,nb)+sunda(k)                !  Indian
       else if (nb.eq.3) then
        do 39 i=indoi+1,imax(nb)
 39     flux(k,i,nb)=flux(k,i,nb)-sunda(k)                !  Pacific
       end if
      end if
c
 2    continue
c
c --- subtract heat transport through indonesian passage
      if (nbasin.gt.2) then
       if (nb.eq.2) then
        do 41 i=indoi,imax(nb)
        saltfl(i,nb)=saltfl(i,nb)+susalt(kk+1)                !  Indian
 41     heatfl(i,nb)=heatfl(i,nb)+suheat(kk+1)                !  Indian
       else if (nb.eq.3) then
        do 51 i=indoi+1,imax(nb)
        saltfl(i,nb)=saltfl(i,nb)-susalt(kk+1)                !  Pacific
 51     heatfl(i,nb)=heatfl(i,nb)-suheat(kk+1)                !  Pacific
       end if
      end if
c
c --- contract mass and heat flux arrays in i direction
c
      imn=imin(nb)
      imx=imax(nb)
      thru=indoi
c
ccc      do 4 it=1,int(log(float(ii1))/log(2.)-4.5)
      do 4 it=1,2
      thru=reviof(thru)
c
      imno=imn
      imxo=imx
      imn=reviof(float(imn))
      imx=reviof(float(imx))
c
      do 25 i=imn,imx
      io=iof(i)
      ioa=max(io-1,imno)
      iob=min(io+1,imxo)
      io=max(ioa,min(io,iob))
      saltfl(i,nb)=.5*saltfl(io,nb)
     .   +.25*(saltfl(ioa,nb)+saltfl(iob,nb))
      heatfl(i,nb)=.5*heatfl(io,nb)
     .   +.25*(heatfl(ioa,nb)+heatfl(iob,nb))
      do 25 k=1,kk
 25   flux(k,i,nb)=.5*flux(k,io,nb)
     .   +.25*(flux(k,ioa,nb)+flux(k,iob,nb))
c
 4    continue
c
      write (lp,100) ocean(nb),' northward heat flux (petawatts):',
     .  (heatfl(i,nb),i=imn,imx)
      write (lp,100) ocean(nb),' northward salt flux (Sv):',
     .  (saltfl(i,nb),i=imn,imx)
 100  format (2a/(11f7.3))
c
      if (nb.eq.1) then
c
c --- find scale factor for printing mass flux streamfunction
        scale=0.
        do 3 i=imn,imx
        do 3 k=1,kk
 3      scale=max(scale,2.*abs(flux(k,i,nb)))
        scale=10.**(int(log10(scale))-2)
      end if
c
      do 33 i=1,ii
      do 33 k=1,kdm
 33   maskk(k,i)=1
c
      write (text,101) 'merid.overturn.(',scale,') -- ',ocean(nb)
 101  format (a,1p,e7.1,2a)
      call prtmsk(maskk,flux(1,1,nb),work,kdm,kdm,imx,0.,1./scale,
     .   text)
c
      call zebra(flux(1,1,nb),kdm,kdm,imx)
c
      if (nbasin.gt.2 .and. nb.eq.nbasin+1) then
        write (lp,'(a,f9.4/(1x,19i4))')
     .  'indonesian throughflow by layer (.1 sv); passage at  i =',
     .  thru,(int(sunda(k)/scale+sign(.5,sunda(k))),k=1,kk)
        write (lp,'(i5,a)') int(sunda(kk+1)/scale+sign(.5,sunda(kk+1))),
     .  '  total'
      end if
c
c --- save everything in a special file
      keep=lp
      lp=no
c
      write (lp,100) ocean(nb),' northward heat flux (petawatts):',
     .  (heatfl(i,nb),i=1,imx)
      write (lp,100) ocean(nb),' northward salt flux (Sv):',
     .  (saltfl(i,nb),i=1,imx)
      write (text,101) 'merid.overturn.(',scale,') -- ',ocean(nb)
      call prtmsk(maskk,flux(1,1,nb),work,kdm,kdm,imx,0.,1./scale,
     .   text)
c
      if (nbasin.gt.2 .and. nb.eq.nbasin+1) then
        write (lp,'(a,f7.2/(1x,19i4))')
     .  'indonesian throughflow by layer (.1 sv); passage at  i =',
     .  thru,(int(sunda(k)/scale+sign(.5,sunda(k))),k=1,kk)
      end if
c
      lp=keep
 15   continue
c
c --- compute throughflow through various passages
c
      x1=thrufl(idrk1,jdrk1,idrk2,jdrk2,'(Drake Passage)')
      x3=thrufl(indoi,indoj1,indoi,indoj2,'(Indonesia)')
c
      keep=lp
      lp=no
      write (lp,'(f7.1,3(f7.1,2x,a))') sunda(kk+1),
     .  x3,'(Indonesia)',x1,'(Drake Passage)'
      lp=keep
c
      close (unit=no)
      write (lp,'(2a)') 'saved the above in  ',flnm
      return
      end
c
c
      real function thrufl(iaa,jaa,ibb,jbb,text)
c
c --- compute thruflow through section (iaa,jaa) - (ibb,jbb)
c --- (it is recommended to put both end points on land)
c
      implicit none
#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"
c
      real flo1,flo2
      integer iaa,jaa,ibb,jbb,iam1,iap1,jam1,jap1,ibm1,ibp1,jbm1,jbp1
      character text*(*)
c
ccc      write (lp,'(a,6i5)')
ccc     .   'thrufl called with iaa,jaa,ibb,jbb =',
ccc     .                       iaa,jaa,ibb,jbb
c --- reverse points if iaa > ibb
      if (iaa.gt.ibb) then
        ia=ibb
        ib=iaa
        ja=jbb
        jb=jaa
      else
        ia=iaa
        ib=ibb
        ja=jaa
        jb=jbb
      end if
c
      flo1=0.
      flo2=0.
c                                !                                 x
c                                !                                  x
      if (jb.ge.ja) then         !  cross section orientation:       x
c                                !                                    x
c                                !                                     x
      jap1=mod(ja,jj)+1
      do 2 i=ia,ib
      if (iv(i,jap1).eq.1) flo1=flo1
     .+vbavg(i,jap1,1)*min(depths(i,ja),depths(i,jap1))*scvx(i,jap1)
 2    continue
      ibm1=ib-1
      do 3 j=ja,jb
      if (iu(ib  ,j).eq.1) flo1=flo1
     .-ubavg(ib  ,j,1)*min(depths(ib,j),depths(ibm1,j))*scuy(ib  ,j)
 3    continue
c
      jbm1=mod(jb-2+jj,jj)+1
      do 4 i=ia,ib
      if (iv(i,jb  ).eq.1) flo2=flo2
     .+vbavg(i,jb  ,1)*min(depths(i,jb),depths(i,jbm1))*scvx(i,jb  )
 4    continue
      iap1=ia+1
      do 5 j=ja,jb
      if (iu(iap1,j).eq.1) flo2=flo2
     .-ubavg(iap1,j,1)*min(depths(ia,j),depths(iap1,j))*scuy(iap1,j)
 5    continue
c                                !                                     x
c                                !                                    x
      else       !  (jb < ja)    !  cross section orientation:       x
c                                !                                  x
c                                !                                 x
      jbp1=mod(jb,jj)+1
      do 6 i=ia,ib
      if (iv(i,jbp1).eq.1) flo1=flo1
     .+vbavg(i,jbp1,1)*min(depths(i,jb),depths(i,jbp1))*scvx(i,jbp1)
 6    continue
      iap1=ia+1
      do 7 j=jb,ja
      if (iu(iap1,j).eq.1) flo1=flo1
     .+ubavg(iap1,j,1)*min(depths(ia,j),depths(iap1,j))*scuy(iap1,j)
 7    continue
c
      jam1=mod(ja-2+jj,jj)+1
      do 8 i=ia,ib
      if (iv(i,ja  ).eq.1) flo2=flo2
     .+vbavg(i,ja  ,1)*min(depths(i,ja),depths(i,jam1))*scvx(i,ja  )
 8    continue
      ibm1=ib-1
      do 9 j=jb,ja
      if (iu(ib  ,j).eq.1) flo2=flo2
     .+ubavg(ib  ,j,1)*min(depths(ib,j),depths(ibm1,j))*scuy(ib  ,j)
 9    continue
c
      end if
c --- convert to sverdrups (ubavg,vbavg in m/s):
      flo1=flo1*sign(1,ibb-iaa)*1.e-6
      flo2=flo2*sign(1,ibb-iaa)*1.e-6
c
      write (lp,'(2f6.1,2(a,2i5),a,2x,a)') flo1,flo2,
     .   '  transport between (',iaa,jaa,') and (',ibb,jbb,')',text
c
      thrufl=.5*(flo1+flo2)
      return
      end
c
c
      subroutine bsnmsk(nbasn,ocean,indoi,indoj1,indoj2)
c
c --- generate -iu- masks for individual ocean basins
c
      implicit none
#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"
c
      integer nbasin,nbasn
      parameter (nbasin=3)
      character fmt*12,char2*2,ocean(nbasin+1)*8
      integer mask(nbasin+1,idm,jdm),iseed(nbasin),jseed(nbasin),
     .        ncnt,jsec,jfrst,jlast,jc,nb,mb,
     .        iub,num,imin,imax,indoi,indoj1,indoj2
      data fmt/'(i4,1x,75i1)'/
c
      common /basins/ iub(nbasin+1,idm,jdm),num(idm,nbasin+1),
     .       imin(nbasin+1),imax(nbasin+1)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- define one ocean point at northern tip of each basin
ccc      data iseed/1,49,1/,jseed/1,65,137/        !  3 basins  1.4deg
      data iseed/71,104,77/,jseed/164,32,91/        !  3 basins  2.0deg (ieq=115)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      if (nbasin.ne.nbasn) then
        write (lp,*) 'basin number mismatch in overtn.f and bsnmsk.f'
        stop
      end if
c
      do 211 j=1,jj
      do 21 i=1,ii
      do 21 nb=1,nbasin+1
      num(i,nb)=0
      iub(nb,i,j)=0
 21   mask(nb,i,j)=0
      do 211 i=1,ii
 211  mask(nbasin+1,i,j)=ip(i,j)
c
c --- place 'seed' near northern tip of each basin. allow seed to propagate
c --- south. ocean basins end where seeds from different basins make contact
c
      do 27 nb=1,nbasin
      mask(nb,iseed(nb),jseed(nb))=1
      write (lp,'(2a,2i6)') ocean(nb),' seed placed at',
     .   iseed(nb),jseed(nb)
      if (depths(iseed(nb),jseed(nb)).le.0.) then
        write (lp,'(2a,2i5)') ocean(nb),' seed landlocked:',
     .    iseed(nb),jseed(nb)
        stop '(overtn initlzn)'
      end if
      imin(nb)=ii
 27   imax(nb)=ii
      imin(nbasin+1)=ii
      imax(nbasin+1)=0
c
c --- propagate 'seed' southward and laterally through each basin.
c
      do 8 i=2,ii1
      do 110 nb=1,nbasin
c
      ncnt=0
      do 9 jc=-jj+1,jj
      j=mod(jc-1+jj,jj)+1
      ja=mod(j-2+jj,jj)+1
c --- block indonesian passage
      if (nbasin.gt.2 .and. nb.eq.3 .and. i.eq.indoi
     .            .and. j.ge.indoj1 .and. j.le.indoj2) go to 9
      ncnt=ncnt-mask(nb,i,j)
      if (depths(i,j).gt.0.)
     .  mask(nb,i,j)=max(mask(nb,i,j),mask(nb,i-1,j),mask(nb,i,ja))
      ncnt=ncnt+mask(nb,i,j)
 9    continue
c
      do 10 jc=jj,-jj+1,-1
      j=mod(jc-1+jj,jj)+1
      jb=mod(j     ,jj)+1
c --- block indonesian passage
      if (nbasin.gt.2 .and. nb.eq.3 .and. i.eq.indoi
     .            .and. j.ge.indoj1 .and. j.le.indoj2) go to 10
      ncnt=ncnt-mask(nb,i,j)
      if (depths(i,j).gt.0.)
     .  mask(nb,i,j)=max(mask(nb,i,j),mask(nb,i-1,j),mask(nb,i,jb))
      ncnt=ncnt+mask(nb,i,j)
 10   continue
      write (lp,'(i5,a,i6,7x,a)') ncnt,' pts added in row',i,ocean(nb)
 110  continue
c
c --- individual basins end at the latitude where they start overlapping
c
      do 8 j=1,jj
      do 8 nb=1,nbasin
      do 8 mb=nb+1,nbasin
      if (mask(nb,i,j)+mask(mb,i,j).gt.1) then
        if (i.lt.imax(nb)) then
          write (lp,'(2a,i5)') ocean(nb),' basin ends at i =',i
          imax(nb)=i
        end if
        if (i.lt.imax(mb)) then
          write (lp,'(2a,i5)') ocean(mb),' basin ends at i =',i
          imax(mb)=i
        end if
      end if
 8    continue
c
      do 5 nb=1,nbasin
c
      do 17 j=1,jj
      do 17 i=imax(nb),ii
 17   mask(nb,i,j)=0
c
c --- do northward sweep to catch remaining parts of each ocean basin
c
      do 22 i=imax(nb)-2,1,-1
c
      do 19 j=1,jj
      ja=mod(j-2+jj,jj)+1
c --- block indonesian passage
      if (nbasin.gt.2 .and. nb.eq.2 .and. i.eq.indoi-1
     .            .and. j.ge.indoj1 .and. j.le.indoj2) go to 19
      if (depths(i,j).gt.0.)
     .  mask(nb,i,j)=max(mask(nb,i,j),mask(nb,i+1,j),mask(nb,i,ja))
 19   continue
c
      do 20 j=jj,1,-1
      jb=mod(j     ,jj)+1
c --- block indonesian passage
      if (nbasin.gt.2 .and. nb.eq.2 .and. i.eq.indoi-1
     .            .and. j.ge.indoj1 .and. j.le.indoj2) go to 20
      if (depths(i,j).gt.0.)
     .  mask(nb,i,j)=max(mask(nb,i,j),mask(nb,i+1,j),mask(nb,i,jb))
 20   continue
 22   continue
c
c --- do another southward sweep
c
      do 18 i=2,imax(nb)-1
c
      do 23 j=1,jj
      ja=mod(j-2+jj,jj)+1
c --- block indonesian passage
      if (nbasin.gt.2 .and. nb.eq.3 .and. i.eq.indoi
     .            .and. j.ge.indoj1 .and. j.le.indoj2) go to 23
      if (depths(i,j).gt.0.)
     .  mask(nb,i,j)=max(mask(nb,i,j),mask(nb,i-1,j),mask(nb,i,ja))
 23   continue
c
      do 24 j=jj,1,-1
      jb=mod(j     ,jj)+1
c --- block indonesian passage
      if (nbasin.gt.2 .and. nb.eq.3 .and. i.eq.indoi
     .            .and. j.ge.indoj1 .and. j.le.indoj2) go to 24
      if (depths(i,j).gt.0.)
     .  mask(nb,i,j)=max(mask(nb,i,j),mask(nb,i-1,j),mask(nb,i,jb))
 24   continue
c
      do 18 j=1,jj
      iub(nb,i,j)=iu(i,j)*max(mask(nb,i,j),mask(nb,i-1,j))
 18   continue
      do 5 j=1,jj
 5    iub(nb,imax(nb),j)=iu(imax(nb),j)*mask(nb,imax(nb)-1,j)
c
      do 171 j=1,jj
      do 171 i=1,ii
 171  iub(nbasin+1,i,j)=iu(i,j)
c
      do 34 nb=1,nbasin+1
      do 13 i=1,ii
      num(i,nb)=0
      do 131 j=1,jj
 131  if (iub(nb,i,j).eq.1) num(i,nb)=num(i,nb)+iub(nb,i,j)
      if (num(i,nb).gt.0) then
        imin(nb)=min(i,imin(nb))
        imax(nb)=max(i,imax(nb))
ccc        write (lp,'(a,i9,a,i5)') ocean(nb),num(i,nb),' points in row',i
      end if
 13   continue
 34   write (lp,'(2a,2i7)') ocean(nb),' basin limits:',imin(nb),imax(nb)
c
c --- print out basin mask arrays
c --- data are written in strips 75 points wide
      do 12 nb=1,nbasin+1
      jsec=(jj-1)/75
c
      do 16 jfrst=0,75*jsec,75
      jlast=min(jj,jfrst+75)
      write (char2,'(i2)') jlast-jfrst
      fmt(8:9)=char2
      write (lp,'(a,'' p mask, cols'',i5,'' --'',i5)') ocean(nb),
     .   jfrst+1,jlast 
 16   write (lp,fmt) (i,(10*mask(nb,i,j),j=jfrst+1,jlast),i=1,ii1)
c
      do 36 jfrst=0,75*jsec,75
      jlast=min(jj,jfrst+75)
      write (char2,'(i2)') jlast-jfrst
      fmt(8:9)=char2
      write (lp,'(a,'' u mask, cols'',i5,'' --'',i5)') ocean(nb),
     .   jfrst+1,jlast 
 36   write (lp,fmt) (i,(10*iub(nb,i,j),j=jfrst+1,jlast),i=1,ii1)
 12   continue
c
      return
      end
c
c
c> Revision history:
c
c> May 2003 - fixed thruflo function
