      subroutine archiv(n,nn)
c
c --- write archive file for time level n to flnm ( b i n a r y  hycom fmt)
c
      USE MODEL_COM, only :
     *  itime,iyear1,nday,jdendofm,jyear,jmon,jday,jdate,jhour,aMON
     * ,xlabel,lrunid
c
      implicit none
#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"
#include "cpl.h"
c
      integer no,nop,length
      real factor,vol,tts2,temavg,sst
     .    ,icearea,icevol,icearean,icevoln,iceareas,icevols
      character flnm*20,intvl*5
      real*4 real4(idm,jdm)
     .   ,time4,watcum4,empcum4,thref4,thbase4,theta4(kdm)
      integer*4 length4,idm4,jdm4,kdm4,nstep4
c
      call getdte(Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)
c --- check if ogcm date matches agcm date
      if (abs((itime+1.)/nday-time).gt.1.e-5) then
        write(*,*) 'mismatching archive date in agcm/ogcm=',
     .     (itime+1.)/nday,time
        stop 'mismatching archive date'
      endif
c
      write(flnm,'(a3,i4.4,2a)') amon,Jyear,'.out',xlabel(1:lrunid)
c
      write (lp,*) 'shown below: sea surface height'
      call zebra(util1,idm,ii1,jj)
      write (lp,*) 'shown below: sea surface temperature'
      call zebra(temp(1,1,1+nn),idm,ii1,jj)
c
      write (intvl,'(i4.4,a)') jdate,' ' 
      nop=13
      no=4096
      length=((idm*jdm+no+15)/no)*no           ! f90 on compac for real*4
c
      write (lp,'(a/9x,a)') 'storing history data in',flnm
c
      open (unit=nop,file=flnm,status='unknown',form='unformatted',
     .      access='direct',recl=length)
      no=1
      length4=length
      idm4=idm
      jdm4=jdm
      kdm4=kdm
      nstep4=nstep
      time4=time
      thbase4=thbase
      do k=1,kk
        theta4(k)=theta(k)
      end do
      write (nop,rec=no) length4,idm4,jdm4,kdm4,nstep4,time4
     .      ,thbase4,theta4
c
      no=no+1
      call r8tor4(ubavg(1,1,n),real4)
      write (nop,rec=no) 'ubavg       ',1,real4
      no=no+1
      call r8tor4(vbavg(1,1,n),real4)
      write (nop,rec=no) 'vbavg       ',1,real4
ccc   no=no+1
ccc   call r8tor4(pbavg(1,1,n),real4)
ccc   write (nop,rec=no) 'pbavg       ',1,real4
      no=no+1
      call r8tor4(util1,real4)
      write (nop,rec=no) 'srfht       ',1,real4
      no=no+1
c
      call r8tor4(dpmixl(1,1),real4)
c
      write (nop,rec=no) 'mix_dpth(m) ',1,real4
      no=no+1
      call r8tor4(oice,real4)
      write (nop,rec=no) 'icecover(%) ',1,real4
c
      do 75 k=1,kk
      kn=k+nn
      no=no+1
      call r8tor4(u(1,1,kn),real4)
      write (nop,rec=no) 'u           ',k,real4
      no=no+1
      call r8tor4(v(1,1,kn),real4)
      write (nop,rec=no) 'v           ',k,real4
      no=no+1
      call r8tor4(dp(1,1,kn),real4)
      write (nop,rec=no) 'dp          ',k,real4
      no=no+1
      call r8tor4(temp(1,1,kn),real4)
      write (nop,rec=no) 'temp        ',k,real4
      no=no+1
      call r8tor4(saln(1,1,kn),real4)
      write (nop,rec=no) 'saln        ',k,real4
      no=no+1
      call r8tor4(th3d(1,1,kn),real4)
      write (nop,rec=no) 'th3d        ',k,real4
      no=no+1
      call r8tor4(tracer(1,1,k),real4)
      write (nop,rec=no) 'tracer      ',k,real4
 75   continue
c
c --- output time-averaged fields
c
      factor=baclin/(jdate*86400.)
c
c$OMP PARALLEL DO
      do 55 j=1,jj
      do 55 l=1,isp(j)
      do 55 i=ifp(j,l),ilp(j,l)
      eminpav(i,j)=eminpav(i,j)*factor
      surflav(i,j)=surflav(i,j)*factor
      tauxav(i,j)=tauxav(i,j)*factor
      tauyav(i,j)=tauyav(i,j)*factor
 55   continue
c$OMP END PARALLEL DO
c
      do 58 k=1,kk
c
c$OMP PARALLEL DO
      do 59 j=1,jj
      do 591 l=1,isu(j)
      do 591 i=ifu(j,l),ilu(j,l)
      if (dpuav(i,j,k).gt.0) uav(i,j,k)=uav(i,j,k)/dpuav(i,j,k)
 591  uflxav(i,j,k)=uflxav(i,j,k)*baclin
      do 592 l=1,isv(j)
      do 592 i=ifv(j,l),ilv(j,l)
      if (dpvav(i,j,k).gt.0) vav(i,j,k)=vav(i,j,k)/dpvav(i,j,k)
 592  vflxav(i,j,k)=vflxav(i,j,k)*baclin
      do 59 l=1,isp(j)
      do 59 i=ifp(j,l),ilp(j,l)
      if (dpav(i,j,k).gt.0.) then
        temav(i,j,k)=temav(i,j,k)/dpav(i,j,k)
        salav(i,j,k)=salav(i,j,k)/dpav(i,j,k)
        th3av(i,j,k)=th3av(i,j,k)/dpav(i,j,k)
      end if
      dpav(i,j,k)=dpav(i,j,k)*factor
c
      diaflx(i,j,k)=-diaflx(i,j,k)/(2.*onem)
c --- convert diapycnal thickness changes into actual interface fluxes
      if (k.gt.1) diaflx(i,j,k)=diaflx(i,j,k)+diaflx(i,j,k-1)
 59   continue
c$OMP END PARALLEL DO
ccc      write (lp,'(a,i3)') 'shown below: N.Atl. diaflx, bottm of layer',k
ccc      call zebra(diaflx(1,int(.8*jdm),k),idm,idm/3,idm/3)
c
      no=no+1
      call r8tor4(uflxav(1,1,k),real4)
      write (nop,rec=no) 'uflxav_'//intvl,k,real4
      no=no+1
      call r8tor4(vflxav(1,1,k),real4)
      write (nop,rec=no) 'vflxav_'//intvl,k,real4
      no=no+1
      call r8tor4(diaflx(1,1,k),real4)
      write (nop,rec=no) 'diaflx_'//intvl,k,real4
 58   continue
c
c$OMP PARALLEL DO
      do 56 j=1,jj
      do 561 l=1,isu(j)
      do 561 i=ifu(j,l),ilu(j,l)
 561  ubavav(i,j)=ubavav(i,j)*factor
      do 562 l=1,isv(j)
      do 562 i=ifv(j,l),ilv(j,l)
 562  vbavav(i,j)=vbavav(i,j)*factor
      do 56 l=1,isp(j)
      do 56 i=ifp(j,l),ilp(j,l)
      pbavav(i,j)=pbavav(i,j)*factor
 56   montav(i,j)=montav(i,j)*factor+thref*pbavav(i,j)
c$OMP END PARALLEL DO
c
      write (lp,'(3a)') 'shown below: ',intvl,'- day SSH average'
      call zebra(montav,idm,ii1,jj)
      write (lp,'(3a)') 'shown below: ',intvl,'- day SST average'
      call zebra(temav,idm,ii1,jj)
c
      no=no+1
      call r8tor4(ubavav,real4)
      write (nop,rec=no) 'ubavav_'//intvl,1,real4
      no=no+1
      call r8tor4(vbavav,real4)
      write (nop,rec=no) 'vbavav_'//intvl,1,real4
ccc   no=no+1
ccc   call r8tor4(pbavav,real4)
ccc   write (nop,rec=no) 'pbavav_'//intvl,1,real4
      no=no+1
      call r8tor4(montav,real4)
      write (nop,rec=no) 'sfhtav_'//intvl,1,real4
c
      do 57 k=1,kk
      no=no+1
      call r8tor4(uav(1,1,k),real4)
      write (nop,rec=no) '   uav_'//intvl,k,real4
      no=no+1
      call r8tor4(vav(1,1,k),real4)
      write (nop,rec=no) '   vav_'//intvl,k,real4
      no=no+1
      call r8tor4(dpav(1,1,k),real4)
      write (nop,rec=no) '  dpav_'//intvl,k,real4
      no=no+1
      call r8tor4(temav(1,1,k),real4)
      write (nop,rec=no) ' temav_'//intvl,k,real4
      no=no+1
      call r8tor4(salav(1,1,k),real4)
      write (nop,rec=no) ' salav_'//intvl,k,real4
      no=no+1
      call r8tor4(th3av(1,1,k),real4)
      write (nop,rec=no) 'th3dav_'//intvl,k,real4
 57   continue
c
c --- time-averaged surface fluxes:
      no=no+1
      call r8tor4(eminpav,real4)
      write (nop,rec=no) 'eminpav'//intvl,1,real4
      no=no+1
      call r8tor4(surflav,real4)
      write (nop,rec=no) 'surflav'//intvl,1,real4
      no=no+1
      call r8tor4(tauxav,real4)
      write (nop,rec=no) ' tauxav'//intvl,1,real4
      no=no+1
      call r8tor4(tauyav,real4)
      write (nop,rec=no) ' tauyav'//intvl,1,real4
c
      close (unit=nop)
c
c$OMP PARALLEL DO
      do 60 j=1,jj
c
      do 601 i=1,ii
      eminpav(i,j)=0.
      surflav(i,j)=0.
      tauxav(i,j)=0.
      tauyav(i,j)=0.
c
      ubavav(i,j)=0.
      vbavav(i,j)=0.
      pbavav(i,j)=0.
 601  montav(i,j)=0.
c
      do 60 k=1,kk
      do 602 i=1,ii
      uav(i,j,k)=0.
      vav(i,j,k)=0.
      dpuav(i,j,k)=0.
      dpvav(i,j,k)=0.
      dpav (i,j,k)=0.
      temav(i,j,k)=0.
      salav(i,j,k)=0.
      th3av(i,j,k)=0.
      uflxav(i,j,k)=0.
      vflxav(i,j,k)=0.
 602  diaflx(i,j,k)=0.
 60   continue
c$OMP END PARALLEL DO
c
      return
      end
c
      subroutine r8tor4(real8,real4)
c
      implicit none
#include "dimensions.h"
#include "dimension2.h"
c
      real real8(idm,jdm)
      real*4 real4(idm,jdm)
c
c$OMP PARALLEL DO
      do 1 j=1,jdm
      do 1 i=1,idm
 1    real4(i,j)=real8(i,j)
c$OMP END PARALLEL DO
c
      return
      end
c
c
c> Revision history:
c>
c> June 2001 - corrected sign error in diaflx
