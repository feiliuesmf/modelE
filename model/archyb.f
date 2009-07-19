!#include "hycom_mpi_hacks.h"
#include "rundeck_opts.h"
      subroutine archiv(n,nn)
c
c --- write archive file for time level n to flnm ( b i n a r y  hycom fmt)
c
      USE MODEL_COM, only :
     *  itime,iyear1,nday,jdendofm,jyear,jmon,jday,jdate,jhour,aMON
     * ,xlabel,lrunid
      USE HYCOM_SCALARS, only : nstep,time,lp,theta,huge,baclin,onem
     &     ,thref
      USE HYCOM_DIM_GLOB, only : ii1,jj,JDM,kk,isp,ifp,ilp,ntrcr,isu
     &     ,ifu,ilu,isv,ifv,ilv,ii,idm,kdm
      USE HYCOM_ARRAYS_GLOB
#ifdef TRACERS_OceanBiology
      USE obio_com, only :  tracav
#endif
c
      implicit none
      integer i,j,k,l,n,nn,kn
c
      integer no,nop,length,nt
      real factor,vol,tts2,temavg,sst,dpsmo(idm,jdm,kdm)
     .    ,icearea,icevol,icearean,icevoln,iceareas,icevols
      character flnm*40,intvl*3
      character what*12
      real*4 real4(idm,jdm)
     .   ,time4,watcum4,empcum4,thref4,theta4(kdm),unused
      integer*4 length4,idm4,jdm4,kdm4,nstep4
      integer*4 irecl ! specific record lenth, machine dependent
      logical, parameter :: smooth = .true.     ! smooth fields before saving
      data unused/0./
c
      call getdte(Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)
c --- check if ogcm date matches agcm date
      if (nstep.eq.1) then
        write(flnm,'(a3,i4.4,2a)') amon,0,'.out',xlabel(1:lrunid)
      elseif (abs((itime+1.)/nday-time).gt.1.e-5) then
        write(*,*) 'mismatching archive date in agcm/ogcm=',
     .     (itime+1.)/nday,time
        stop 'mismatching archive date'
      else
        write(flnm,'(a3,i4.4,2a)') amon,Jyear,'.out',xlabel(1:lrunid)
      endif
c
c     write (lp,*) 'shown below: sea surface height'
c     call zebra(srfhgt,idm,ii1,jj)
      write (lp,*) 'shown below: sea surface temperature'
      call zebra(temp(1,1,1+nn),idm,ii1,jj)
c
      if (jdate.lt.100) then
        write (intvl,'(i3.3)') jdate
      else
        stop 'jdate >100'
      endif
c
      nop=13
      no=4096 
      inquire (iolength=irecl)  real4(1,1) ! length of an unformatted real*4  
                                           ! irecl=1 on COMPAQ, irecl=4 on SGI
      length=((irecl*idm*JDM+no+15)/no)*no
c
      write (lp,'(a/9x,a)') 'storing history data in',flnm
c
      open (unit=nop,file=flnm,status='unknown',form='unformatted',
     .      access='direct',recl=length)
      no=1
      length4=length
      idm4=idm
      jdm4=JDM
      kdm4=kdm
      nstep4=nstep
      time4=time
      do k=1,kk
        theta4(k)=theta(k)
      end do
      write (nop,rec=no) length4,idm4,jdm4,kdm4,nstep4,time4
     .      ,unused,theta4
c
      if (smooth) then
c
cc$OMP PARALLEL DO SCHEDULE(STATIC,jchunk)
        do j=1,jj
         do k=1,kk
          do l=1,isp(j)
           do i=ifp(j,l),ilp(j,l)
            p(i,j,k+1)=p(i,j,k)+dp(i,j,k+nn)
           end do
          end do
        end do
       end do
cc$OMP END PARALLEL DO
c
       do k=2,kk
        call psmo1(p(1,1,k),pbot)
       end do
c
      dpsmo=huge
cc$OMP PARALLEL DO SCHEDULE(STATIC,jchunk)
       do j=1,jj
        do k=1,kk
         do l=1,isp(j)
          do i=ifp(j,l),ilp(j,l)
            dpsmo(i,j,k)=p(i,j,k+1)-p(i,j,k)
          end do
         end do
        end do
       end do
cc$OMP END PARALLEL DO
c
      end if                            !  smooth
c
      no=no+1
      call r8tor4(ubavg(1,1,n),real4)
      if (smooth) call usmoo4(real4)
      write (nop,rec=no) 'ubavg       ',0,real4
      write (lp,100)     'ubavg       ',0,no
      no=no+1
      call r8tor4(vbavg(1,1,n),real4)
      if (smooth) call vsmoo4(real4)
      write (nop,rec=no) 'vbavg       ',0,real4
      write (lp,100)     'vbavg       ',0,no
      no=no+1
      call r8tor4(srfhgt,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'srfhgt      ',0,real4
      write (lp,100)     'srfhgt      ',0,no
      no=no+1
      call r8tor4(dpmixl,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'mix_dpth    ',0,real4
      write (lp,100)     'mix_dpth    ',0,no
      no=no+1
      call r8tor4(oice,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'icecover(%) ',0,real4
      write (lp,100)     'iceconver(%)',0,no
c
      do 75 k=1,kk
      kn=k+nn
      no=no+1
      call r8tor4(u(1,1,kn),real4)
      if (smooth) call usmoo4(real4)
      write (nop,rec=no) 'u           ',k,real4
      write (lp,100)     'u           ',k,no
      no=no+1
      call r8tor4(v(1,1,kn),real4)
      if (smooth) call vsmoo4(real4)
      write (nop,rec=no) 'v           ',k,real4
      write (lp,100)     'v           ',k,no
      no=no+1
      call r8tor4(dpsmo(1,1,k),real4)
      write (nop,rec=no) 'dp          ',k,real4
      write (lp,100)     'dp          ',k,no
      no=no+1
      call r8tor4(temp(1,1,kn),real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'temp        ',k,real4
      write (lp,100)     'temp        ',k,no
      no=no+1
      call r8tor4(saln(1,1,kn),real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'saln        ',k,real4
      write (lp,100)     'saln        ',k,no
      no=no+1
      call r8tor4(th3d(1,1,kn),real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'th3d        ',k,real4
      write (lp,100)     'th3d        ',k,no
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- ifort miscomputes 'no' in the following loop:
ccc      do nt=1,ntrcr
ccc        no=no+1
ccc        call r8tor4(tracer(1,1,k,nt),real4)
ccc        if (smooth) call psmoo4(real4)
ccc        write (what,'(a6,i2,4x)') 'tracer',nt
ccc        write (nop,rec=no) what,k,real4
ccc      write (lp,100)       what,k,no
ccc      end do
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- temporary: store diffusivity as tracer 1
      tracer(:,:,:,1)=diadff*1.e4
      write (what,'(a12)') 'diapyc.diffu'
c
c --- code around compiler glitch:
      do nt=1,ntrcr
        call r8tor4(tracer(1,1,k,nt),real4)
        if (smooth) call psmoo4(real4)
        write (what,'(a6,i2,4x)') 'tracer',nt
        write (nop,rec=no+nt) what,k,real4
      write (lp,100)       what,k,no+nt
      end do
      no=no+ntrcr
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 75   continue
c
c --- output time-averaged fields
c
      factor=baclin/(jdate*86400.)
c
cc$OMP PARALLEL DO SCHEDULE(STATIC,jchunk)
      do 55 j=1,jj
      do 55 l=1,isp(j)
      do 55 i=ifp(j,l),ilp(j,l)
      eminpav(i,j)=eminpav(i,j)*factor
      surflav(i,j)=surflav(i,j)*factor
      salflav(i,j)=salflav(i,j)*factor
      brineav(i,j)=brineav(i,j)*factor
      tauxav(i,j)=tauxav(i,j)*factor
      tauyav(i,j)=tauyav(i,j)*factor
 55   continue
cc$OMP END PARALLEL DO
c
      do 58 k=1,kk
c
cc$OMP PARALLEL DO SCHEDULE(STATIC,jchunk)
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
#ifdef TRACERS_OceanBiology
        do nt=1,ntrcr
           tracav(i,j,k,nt)=tracav(i,j,k,nt)/dpav(i,j,k)
        enddo
#endif
      end if
      dpav(i,j,k)=dpav(i,j,k)*factor
c
      diaflx(i,j,k)=-diaflx(i,j,k)/(2.*onem)
c --- convert diapycnal thickness changes into actual interface fluxes
      if (k.gt.1) diaflx(i,j,k)=diaflx(i,j,k)+diaflx(i,j,k-1)
 59   continue
cc$OMP END PARALLEL DO
ccc      write (lp,'(a,i3)') 'shown below: N.Atl. diaflx, bottm of layer',k
ccc      call zebra(diaflx(1,int(.8*jdm),k),idm,idm/3,idm/3)
c
      no=no+1
      call r8tor4(uflxav(1,1,k),real4)
      write (nop,rec=no) '  uflxav_'//intvl,k,real4
      write (lp,100)     '  uflxav_'//intvl,k,no
      no=no+1
      call r8tor4(vflxav(1,1,k),real4)
      write (nop,rec=no) '  vflxav_'//intvl,k,real4
      write (lp,100)     '  vflxav_'//intvl,k,no
      no=no+1
      call r8tor4(diaflx(1,1,k),real4)
      write (nop,rec=no) '  diaflx_'//intvl,k,real4
      write (lp,100)     '  diaflx_'//intvl,k,no
 58   continue
c
cc$OMP PARALLEL DO SCHEDULE(STATIC,jchunk)
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
      sfhtav(i,j)=sfhtav(i,j)*factor+thref*pbavav(i,j)
      dpmxav(i,j)=dpmxav(i,j)*factor
 56   oiceav(i,j)=oiceav(i,j)*factor
cc$OMP END PARALLEL DO
c
c     write (lp,'(3a)') 'shown below: ',intvl,'- day SSH average'
c     call zebra(sfhtav,idm,ii1,jj)
      write (lp,'(3a)') 'shown below: ',intvl,'- day SST average'
      call zebra(temav,idm,ii1,jj)
c
      no=no+1
      call r8tor4(ubavav,real4)
      write (nop,rec=no) '  ubavav_'//intvl,0,real4
      write (lp,100)     '  ubavav_'//intvl,0,no
      no=no+1
      call r8tor4(vbavav,real4)
      write (nop,rec=no) '  vbavav_'//intvl,0,real4
      write (lp,100)     '  vbavav_'//intvl,0,no
      no=no+1
      call r8tor4(sfhtav,real4)
      write (nop,rec=no) '  sfhtav_'//intvl,0,real4
      write (lp,100)     '  sfhtav_'//intvl,0,no
      no=no+1
      call r8tor4(dpmxav,real4)
      write (nop,rec=no) '  dpmxav_'//intvl,0,real4
      write (lp,100)     '  dpmxav_'//intvl,0,no
      no=no+1
      call r8tor4(oiceav,real4)
      write (nop,rec=no) '  oiceav_'//intvl,1,real4
      write (lp,100)     '  oiceav_'//intvl,0,no
c
      do 57 k=1,kk
      no=no+1
      call r8tor4(uav(1,1,k),real4)
      write (nop,rec=no) '     uav_'//intvl,k,real4
      write (lp,100)     '     uav_'//intvl,k,no
      no=no+1
      call r8tor4(vav(1,1,k),real4)
      write (nop,rec=no) '     vav_'//intvl,k,real4
      write (lp,100)     '     vav_'//intvl,k,no
      no=no+1
      call r8tor4(dpav(1,1,k),real4)
      write (nop,rec=no) '    dpav_'//intvl,k,real4
      write (lp,100)     '    dpav_'//intvl,k,no
      no=no+1
      call r8tor4(temav(1,1,k),real4)
      write (nop,rec=no) '   temav_'//intvl,k,real4
      write (lp,100)     '   temav_'//intvl,k,no
      no=no+1
      call r8tor4(salav(1,1,k),real4)
      write (nop,rec=no) '   salav_'//intvl,k,real4
      write (lp,100)     '   salav_'//intvl,k,no
      no=no+1
      call r8tor4(th3av(1,1,k),real4)
      write (nop,rec=no) '  th3dav_'//intvl,k,real4
      write (lp,100)     '  th3dav_'//intvl,k,no
 57   continue
c
c --- time-averaged surface fluxes:
      no=no+1
      call r8tor4(eminpav,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) ' emp m/s '//intvl,0,real4
      write (lp,100)     ' emp m/s '//intvl,0,no
      no=no+1
      call r8tor4(surflav,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'hflx W/m2'//intvl,0,real4
      write (lp,100)     'hflx W/m2'//intvl,0,no
      no=no+1
      call r8tor4(salflav,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'sfx g/m2s'//intvl,0,real4
      write (lp,100)     'sfx g/m2s'//intvl,0,no
      no=no+1
      call r8tor4(brineav,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'brn g/m2s'//intvl,0,real4
      write (lp,100)     'brn g/m2s'//intvl,0,no
      no=no+1
      call r8tor4(tauxav,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'taux N/m2'//intvl,0,real4
      write (lp,100)     'taux N/m2'//intvl,0,no
      no=no+1
      call r8tor4(tauyav,real4)
      if (smooth) call psmoo4(real4)
      write (nop,rec=no) 'tauy N/m2'//intvl,0,real4
      write (lp,100)     'tauy N/m2'//intvl,0,no
c
      close (unit=nop)
 100  format (9x,a,' (layer',i3,') archived as record',i5)
      write (lp,*) no,' records archived'
c


#ifdef TRACERS_OceanBiology
      call obio_archyb(nn,dpav,temav,salav,th3av,tracav,dpmxav,oiceav)
#endif

cc$OMP PARALLEL DO SCHEDULE(STATIC,jchunk)
      do 60 j=1,jj
c
      do 601 i=1,ii
      eminpav(i,j)=0.
      surflav(i,j)=0.
       salflav(i,j)=0.
      brineav(i,j)=0.
      tauxav(i,j)=0.
      tauyav(i,j)=0.
c
      ubavav(i,j)=0.
      vbavav(i,j)=0.
      pbavav(i,j)=0.
      dpmxav(i,j)=0.
      sfhtav(i,j)=0.
 601  oiceav(i,j)=0.
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
cc$OMP END PARALLEL DO
c
      return
      end
c
      subroutine r8tor4(real8,real4)
c
      USE HYCOM_DIM_GLOB
      implicit none
      integer i,j
c
      real real8(idm,jdm)
      real*4 real4(idm,jdm)
c
cc$OMP PARALLEL DO SCHEDULE(STATIC,jchunk)
      do 1 j=1,jdm
      do 1 i=1,idm
 1    real4(i,j)=real8(i,j)
cc$OMP END PARALLEL DO
c
      return
      end
c
c
c> Revision history:
c>
c> June 2001 - corrected sign error in diaflx
c> Feb. 2005 - added multiple tracer capability
c> June 2005 - reordered output fields, added time-av'ged mxlyr & ice temp/th
