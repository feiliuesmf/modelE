      subroutine OCEANS
c
c --- ------------------------------
c --- MICOM-based hybrid ocean model
c
c ---     v e r s i o n    0.9
c --- ------------------------------
c
c                            Program Notes:
c
c                    Online Tracers in MICOM and HYCOM
c                    ---------------------------------
c
c 1. The number of tracers is set by 'ntrcr' in dimensions.h.
c
c 2. Source/sink functions for each tracer must be supplied by the user.
c In one current implementation, ntrcr=2, and there is a call in hycom.f
c to a routine called 'addtrc' which continuously adds tracer at 2
c "pipeline outfall" points on the sea floor. In another implementation,
c ntrcr=21, and hycom calls a biochemical submodel encompassing 21
c biological and chemical tracers.
c
c 3. Tracer transport takes place once every 'trcfrq' time steps; in
c other words, the transport time step is trcfrq*baclin. The transport
c equation is solved in flux form and is based on a time-integrated form
c of the mass continuity equation in which interlayer mass exchange is
c determined as a residual based on cumulative layer thickness change
c and time-integrated horizontal mass flux divergence.
c
c 4. Tracer transport in  n o t  done in leapfrog mode. In the present
c implementation, 'trcfrq' must be even. Only mass fluxes at even time
c steps are used in building up the mass flux time integrals. This is
c to avoid leapfrog-related inconsistencies between layer thickness
c tendencies and horizontal flux fields which would generate spurious
c diapycnal mass fluxes. The single remaining inconsistency is the one
c caused by time smoothing of the thickness field. It spawns small
c diapyncal fluxes (and hence interlayer tracer leakage) even in the
c absence of physically based vertical mixing processes.
c
c 5. Transport is handled by subroutine trcadv.f which has 3 entries:
c
c    -  tradv0 is called at the beginning of each 'long' time interval.
c       It initializes the mass flux time integrals and stores the
c       initial layer thickness.
c
c    -  tradv1 is called every even-numbered 'regular' model time step.
c       It builds up the mass flux time integrals.
c
c    -  tradv2 is called at the end of each 'long' time interval. It
c       diagnoses the diapycnal mass flux components and carries out the
c       actual 3-D tracer transport, using a combination of FCT (flux
c       corrected transport)in the horizontal and PPM (piecewise parabolic
c       method) in the vertical.
c
c 6. Solving transport equations in flux form makes it impossible to
c rigorously conserve total tracer when layer thickness is allowed to go
c to zero. Conservation errors are currently eliminated by applying a
c global correction factor to each tracer field following each transport
c step. With this corrective device in place, tracer is rigorously
c conserved in the model.
c
c 7. Upon exit from tradv2, tracer fields are valid at the current time
c level 'n' and hence can be paired with the current layer thickness
c dp(:,:,k+nn). This temporal consistency presents a window of
c opportunity, which remains open until the next call to tradv0, for
c performing vertical mixing operations on the tracer fields. The logical
c variable 'dotrcr' is set to 'true' during this time window and can be
c used to direct T/S mixing routines to also mix tracer. For numerical
c consistency, the time step for tracer mixing in those routines should
c be set to trcfrq*baclin. Since convective adjustment and Kraus-Turner
c mixed layer routines do not solve a diffusion equation with an explicit
c time step, mathematical consistency of tracer treatment by convec and
c mxlayr is not obvious. These routines operate on the layer thickness
c field every 'regular' time step. Their diffusive effect on the tracer
c fields is therefore captured mainly by diapycnal fluxes spawned by
c layer thickness trends which result from vertical mixing.
c
c 8. Since layer thickness in HYCOM's z-coordinate subdomain does not
c reflect the action of vertical mixing processes the way it does in
c MICOM, vertical mixing by convec and mxlayr is likely to be
c underestimated in HYCOM. This problem is alleviated by using
c vertical mixing schemes like KPP (with time step trcfrq*baclin).
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      USE FLUXES, only : e0,prec,eprec,evapor,flowo,eflowo,dmua,dmva
     . ,erunosi,runosi,srunosi,runpsi,srunpsi,dmui,dmvi,dmsi,dhsi,dssi
     . ,gtemp,sss,mlhc,ogeoza,uosurf,vosurf,MELTI,EMELTI,SMELTI
     . ,gmelt,egmelt,solar
      USE SEAICE_COM, only : rsi,msi
      USE SEAICE, only : fsss,tfrez
      USE GEOM, only : dxyp
      USE MODEL_COM, only : focean
      USE CONSTANT, only : lhm,shi,shw
      USE MODEL_COM, only: dtsrc
     *  ,itime,iyear1,nday,jdendofm,jyear,jmon,jday,jdate,jhour,aMON
c
      implicit none
c
#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"
#include "cpl.h"
#include "a2o.h"
#include "kprf_scalars.h"
#include "kprf_arrays.h"
c
      real sum,coord,x,x1,totl,sumice,fusion,saldif,sofsig,tf
     .    ,sigocn,kappaf,check,apehyc,pechg_hyc_bolus
     .    ,hyc_pechg1,hyc_pechg2,q,sum1,sum2,dpini(kdm)
     .    ,thkchg,flxdiv
      integer jj1,no,index,nflip,mo0,mo1,mo2,mo3,rename,iatest,jatest
     .       ,OMP_GET_NUM_THREADS,io,jo,ipa(iia,jja),nsub
      external rename
      logical master,slave,diag_ape
      character util(idm*jdm+14)*2,charac(20)*1,string*20,
     .          flnm*60
      data charac/'1','2','3','4','5','6','7','8','9','0',
     .            'A','B','C','D','E','F','G','H','I','J'/
      data iatest,jatest/9,7/
c     data iatest,jatest/33,40/           ! Iceland
      data nflip/0/
c
      real cnuity_time,tsadvc_time,momtum_time,barotp_time,trcadv_time,
     .     thermf_time,enloan_time,mxlayr_time,hybgen_time,
     .     agcm_time,ogcm_time
      integer after,before,rate,bfogcm
c
      real osst(idm,jdm),osss(idm,jdm),osiav(idm,jdm)
     . ,oogeoza(idm,jdm),usf(idm,jdm),vsf(idm,jdm)
     . ,utila(iia,jja)

#include "state_eqn.h"
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- initiate named-pipe comparison utility
c --- (see pipe.f for instructions)
ccc      inquire(file='master',exist=master)
ccc      inquire(file='slave',exist=slave)
ccc      if ((master .and. slave) .or. (.not.master .and. .not.slave))
ccc     .   stop '(master/slave ambiguity)'
ccc      if (master) call pipe_init(.true.)
ccc      if (slave) call pipe_init(.false.)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      call getdte(Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)
cd    write(*,'(a,i8,7i5,a)')'chk =',
cd   .    Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon
c
      if (mod(jhour,nhr).eq.0.and.mod(itime,nday/24).eq.0) then
c$OMP PARALLEL DO
        do 28 ia=1,iia
        do 28 ja=1,jja
          ataux(ia,ja)=0.
          atauy(ia,ja)=0.
        aflxa2o(ia,ja)=0.
          aemnp(ia,ja)=0.
          asalt(ia,ja)=0.
          aice(ia,ja)=0.
         austar(ia,ja)=0.
         aswflx(ia,ja)=0.
 28     continue
c$OMP END PARALLEL DO
      endif
c
c$OMP PARALLEL DO
c --- accumulate agcm fields over nhr 
      do 29 ia=1,iia
      do 29 ja=1,jja
      ipa(ia,ja)=0
      if (focean(ia,ja).eq.0.) goto 29
      ipa(ia,ja)=1
c --- accumulate 
      aemnp(ia,ja)=aemnp(ia,ja)                               ! kg/m*m => m/s
     .+((prec(ia,ja)-evapor(ia,ja,1))*(1.-rsi(ia,ja))         ! open water
     .+(flowo(ia,ja)+gmelt(ia,ja)+melti(ia,ja))/(dxyp(ja)*focean(ia,ja))!ocn/ice
     .+(runosi(ia,ja)+runpsi(ia,ja))*rsi(ia,ja))*thref        ! ice
     .                                /(3600.*real(nhr))
      aflxa2o(ia,ja)=aflxa2o(ia,ja)                           ! J/m*m => W/m*m
     . +((e0(ia,ja,1)+eprec(ia,ja))*(1.-rsi(ia,ja))           ! ocean water
     . +(eflowo(ia,ja)+egmelt(ia,ja)+emelti(ia,ja))
     .                              /(dxyp(ja)*focean(ia,ja)) ! ocn or ice
css  . +(erunosi(ia,ja)+erunpsi(ia,ja))*rsi(ia,ja))           ! ice
     . + erunosi(ia,ja)*rsi(ia,ja))                           ! ice
     .                                 /(3600.*real(nhr))
      asalt(ia,ja)=asalt(ia,ja)                            ! kg/m*m/sec salt
     .+((srunosi(ia,ja)+srunpsi(ia,ja))*rsi(ia,ja)         !
     .   +smelti(ia,ja)/(dxyp(ja)*focean(ia,ja)))          !
     .                             /(3600.*real(nhr))
       aice(ia,ja)= aice(ia,ja) + rsi(ia,ja)*dtsrc/(real(nhr)*3600.)
c --- dmua on B-grid, dmui on C-grid; Nick aug04
      ataux(ia,ja)=ataux(ia,ja)+(dmua(ia,ja,1)+dmui(ia,ja))    ! scaled by rsi 
     .                                     /(3600.*real(nhr))  ! kg/ms => N/m*m
      atauy(ia,ja)=atauy(ia,ja)+(dmva(ia,ja,1)+dmvi(ia,ja))
     .                                     /(3600.*real(nhr))  ! kg/ms => N/m*m
      austar(ia,ja)=austar(ia,ja)+(
     . sqrt(sqrt((dmua(ia,ja,1)+dmui(ia,ja))**2
     .          +(dmva(ia,ja,1)+dmvi(ia,ja))**2)/dtsrc*thref)) ! sqrt(T/r)=>m/s
     .                               *dtsrc/(real(nhr)*3600.)
      aswflx(ia,ja)=aswflx(ia,ja)+(solar(1,ia,ja)*(1.-rsi(ia,ja))!J/m*m=>W/m*m
     .                            +solar(3,ia,ja)*    rsi(ia,ja))
     .                                 /(3600.*real(nhr))
 29   continue
c$OMP END PARALLEL DO
c
      nsavea=nsavea+1
      if (mod(jhour,nhr).gt.0.or.mod(itime,nday/24).eq.0) return
      if (nsavea*24/nday.ne.nhr) then
        print *,'nonmatching b.c. accumulation periods: agcm/ogcm:',
     .  nsavea*24/nday,nhr
        stop 'agcm/ogcm accumulating periods differ'
      end if
      nsavea=0
c
      call veca2o(ataux,atauy,taux,tauy)          !wind stress
      call flxa2o(aice,oice)                      !ice coverage
      call flxa2o(aflxa2o,oflxa2o)                !heatflux everywhere
      call flxa2o(asalt,osalt)                    !saltflux from SI
      call flxa2o(aemnp,oemnp)                    !E - P everywhere
      call flxa2o(austar,ustar)                   !friction velocity
      call flxa2o(aswflx,sswflx)                  ! shortwave flux
c
      call system_clock(before)
      call system_clock(count_rate=rate)
      bfogcm=before
      agcm_time = real(bfogcm-afogcm)/real(rate)
      write (lp,99009) agcm_time,' sec for AGCM bfor ocn step ',nstep
99009 format (f12.3,a,i8)
c
      if (nstep.eq.0 .or.nstep.eq.nstep0) then
c
c$OMP PARALLEL SHARED(mo0)
c$    mo0=OMP_GET_NUM_THREADS()
c$OMP END PARALLEL
c$    call OMP_SET_DYNAMIC(.false.)
ccc   call OMP_SET_NUM_THREADS(max(mo0,16))
c$OMP PARALLEL SHARED(mo1)
ccc   mo1=OMP_GET_NUM_THREADS()
c$OMP END PARALLEL
c     write (lp,'(2(a,i5))') ' hycom thread count',mo0
c     write (lp,'(2(a,i5))') ' hycom thread count',mo0,' changed to',mo1
ccc     . ,' =======  number of threads:',mo1,' ======='
c
      lp=6
c
c --- compute eqn.of state check values
c
      check=36.719718               ! T:[-2:30],S:[18:38]
      check=36.876506               ! T:[-2:32],S:[16:38]
      if (abs(sigocn(4.,35.)-check).gt..0001) then
       write (lp,'(/2(a,f11.6))') 'error -- sigocn(t=4,s=35) should be',
     . check,', not',sigocn(4.,35.)
css     stop
      end if
c --- reference: t= 0.0, s=34.0, p=0 bar,kap(4.5,34.5,1.e7)=  0.11954594
      check=.11954594
      check=0.
      if (abs(kappaf(4.5,34.5,1.e7)-check).gt..00001) then
        write (lp,'(/a,2(a,f12.8))') 'error: kappa(4.5,34.5,10^7)',
     .  '  should be',check,', not',kappaf(4.5,34.5,1.e7)
        stop
      end if
c
      write (lp,109) thkdff,temdff,veldff,viscos,diapyc,vertmx
 109  format (' turb. flux parameters:',1p/
     .  ' thkdff,temdff,veldff =',3e9.2/
     .  ' viscos,diapyc,vertmx =',3e9.2)
      write (lp,'(a,1p,e11.3)') 'thbase =',thbase
c
c --- 'lstep' = number of barotropic time steps per baroclinic time step.
c --- lstep   m u s t   be even.
c
      lstep=baclin/batrop
      lstep=2*((lstep+1)/2)
      dlt=baclin/lstep
      nstepi=real(nhr)*3600./baclin + .0001
      write (lp,'(i4,'' barotropic steps per baroclinic time step'')')
     .  lstep
      write (lp,'(''ogcm exchange w. agcm every step,hr'',2i5)')
     .  nstepi,nhr
c
c --- set up parameters defining the geographic environment
c
css   call geopar                ! moved to agcm for zero start or below
css   call inicon                ! moved to agcm
      if (mxlkpp) then
        call inikpp                 ! kpp
        if (mxlgis) stop 'wrong kprf: kpp=gis=true'
        print *,'chk: mxlkpp=true'
      elseif (mxlgis) then
        call inigis                 ! giss
        print *,'chk: mxlgis=true'
      else
        stop 'wrong kprf: kpp=gis=false'
      endif
c
      trcfrq=int(43200./baclin+.0001)
      write (lp,'(a,2i4)') 'trcfrq set to',trcfrq
      if (trcout) dotrcr=.true.
c
      watcum=0.
      empcum=0.
c
      if (jyear-iyear1.le.20) then
        write (lp,'(''starting date in ogcm/agcm '',2i5,'' hr '',2i12
     .   /'' iyear1/jyear='',2i5)')
     .  int((nstep0+nstepi-1)*baclin)/(3600*24),itime/nday
     . ,int((nstep0+nstepi-1)*baclin)/3600,itime*24/nday,iyear1,jyear
      else
        write (lp,'(''starting date in ogcm/agcm '',2i12
     .   /'' iyear1/jyear='',2i5)') 
     .  int((nstep0+nstepi-1)*baclin)/(3600*24),itime/nday
     . ,iyear1,jyear,(nstep0+nstepi-1)*baclin/3600.,itime*24/nday
      endif
c
c     if(abs((nstep0+nstepi-1)*baclin/3600.-itime*24./nday).gt.1.e-5)
      if(abs(int((nstep0+nstepi-1)*baclin/3600)-itime*24/nday).gt.0)
     .                                                        then
        write (lp,'(a,f16.8,i10,f16.8)')'mismatch date found '
     . ,int((nstep0+nstepi-1)*baclin/3600),itime*24/nday
     . ,(nstep0+nstepi-1)*baclin/3600.-itime*24/nday
        write (lp,'(/(a,i9,a,i9,a,f7.1,a))')'chk model start step',nstep0
     . ,' changed to: '
     . ,int((itime*24/nday+(iyear1-jyear)*8760)*3600/baclin)-nstepi+1
     . ,' (day '
     .,(int((itime*24/nday+(iyear1-jyear)*8760)*3600/baclin)-nstepi+1)
     .  *baclin/86400,')'
        nstep0=int((itime*24/nday+(iyear1-jyear)*8760)*3600/baclin)
     .                                         -nstepi+3600/baclin
        nstep=nstep0
        time0=nstep0*baclin/86400
        tavini=time0
c
      do 19 k=1,kk
 19   dpini(k)=dp(itest,jtest,k)+dp(itest,jtest,k+kk)
c
      else
        write (lp,'(/(a,i9))') 'chk model starts at steps',nstep0
      endif
c
c$OMP PARALLEL DO
      do j=1,jj
      do i=1,ii
      osst(i,j)=0.
      osss(i,j)=0.
      osiav(i,j)=0.
      enddo
      enddo
c$OMP END PARALLEL DO
c
      endif       ! if (nstep.eq.0 .or.nstep.eq.nstep0) 
c
      if (nstepi.le.0) stop 'wrong nstepi'
c
c     print *,' shown below oice'
c     call zebra(oice,iio,iio,jjo)
c     print *,' shown below aflxa2o'
c     call zebra(aflxa2o,iia,iia,jja)
c
c$OMP PARALLEL DO PRIVATE(jb)
      do 202 j=1,jj
      jb=mod(j,jj)+1
      do 202 l=1,isp(j)
      do 202 i=ifp(j,l),ilp(j,l)
      osst(i,j)=0.
      osss(i,j)=0.
      osiav(i,j)=0.
c
      hekman(i,j)=ustar(i,j)*(cekman*4.0)/                ! kpp
     &           (abs(corio(i,j  ))+abs(corio(i+1,j  ))+
     &            abs(corio(i,jb ))+abs(corio(i+1,jb )))
 202  continue
c$OMP END PARALLEL DO
c
c --- ---------------------
c --- sub loop starts here
c --- ---------------------
c
c --- letter 'm' refers to mid-time level (example: dp(i,j,km) )
c --- letter 'n' refers to old and new time level
c
      nsaveo=0
      do 15 nsub=1,nstepi
      m=mod(nstep  ,2)+1
      n=mod(nstep+1,2)+1
      mm=(m-1)*kk
      nn=(n-1)*kk
      k1m=1+mm
      k1n=1+nn
      nstep=nstep+1
      time=time0+(nstep-nstep0)*baclin/86400.
c
      diagno=.false.
      if (JDendOfM(jmon).eq.jday.and.Jhour.eq.24.and.nsub.eq.nstepi) 
     .                                    diagno=.true. ! end of month
c
css   if (nstep.eq.1) diagno=.true.    ! initial condition
      diag_ape=diagno
c
      trcadv_time = 0.0
      if (dotrcr) then
        call system_clock(before)      ! time elapsed since last system_clock
c --- initialization of tracer transport arrays (incl. dpinit):
        call tradv0(m,mm)
cdiag print *,'past tradv0'
        dotrcr=.false.
c
c --- inject tracer into ogcm
c$OMP PARALLEL DO
        do 12 j=1,jj
        do 12 l=1,isp(j)
        do 12 i=ifp(j,l),ilp(j,l)
 12     tracer(i,j,1,1)=1.              !  surface ventilation tracer
c$OMP END PARALLEL DO
        call system_clock(after)        ! time elapsed since last system_clock
        trcadv_time = real(after-before)/real(rate)
      end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- available potential energy diagnostics:
      if (diag_ape)
     . write (501,'(f9.1,a,1p,e15.7)') time,'  APE (J):',
     .  hyc_pechg1(dp(1,1,k1n),th3d(1,1,k1n),31)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 103  format (f9.1,a,-12p,f9.3,' TW')
c
      call system_clock(before)
      call cnuity(m,n,mm,nn,k1m,k1n)
      call system_clock(after)
      cnuity_time = real(after-before)/real(rate)
c     call sstbud(0,'  initialization',temp(1,1,k1n))
c
ccc      write (string,'(a12,i8)') 'cnuity, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
      if (trcout) then
        before = after
c --- long time step tracer advection: build up mass flux time integral
        if (n.eq.oddev) call tradv1(n,nn)
cdiag print *,'past tradv1'
c
        if (mod(nstep,trcfrq).eq.0) then
          dotrcr=.true.         !  start tracer advection/turb.mixing cycle
          write (lp,'(a)') 'start tracer advection/turb.mixing cycle'
          before = after
c --- tracer transport:
          call tradv2(n,nn)
        end if
        trcadv_time = trcadv_time + real(after-before)/real(rate)
      end if
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (diag_ape) q=hyc_pechg1(dp(1,1,k1m),th3d(1,1,k1m),32)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      before = after
      call tsadvc(m,n,mm,nn,k1m,k1n)
      call system_clock(after)
      tsadvc_time = real(after-before)/real(rate)
c     call sstbud(1,' horiz.advection',temp(1,1,k1n))
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (diag_ape)
     . write (501,103) time,'  APE change due to time smoothing: ',
     .  hyc_pechg2(dp(1,1,k1m),th3d(1,1,k1m),32)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      write (string,'(a12,i8)') 'tsadvc, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
      before = after
      call momtum(m,n,mm,nn,k1m,k1n)
      call system_clock(after)
      momtum_time = real(after-before)/real(rate)
c
ccc      write (string,'(a12,i8)') 'momtum, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
      before = after
      call barotp(m,n,mm,nn,k1m,k1n)
      call system_clock(after)
      barotp_time = real(after-before)/real(rate)
c
ccc      write (string,'(a12,i8)') 'barotp, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
c     before = after
c     call convec(m,n,mm,nn,k1m,k1n)
c     call system_clock(after)
c     convec_time = real(after-before)/real(rate)
c     call sstbud(2,' convec.adjustmt',temp(1,1,k1n))
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (diag_ape)
     . write (501,103) time,'  APE change due to mass adjustment:',
     .  hyc_pechg2(dp(1,1,k1n),th3d(1,1,k1n),31)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
ccc      write (string,'(a12,i8)') 'convec, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
c     before = after
c     call diapfl(m,n,mm,nn,k1m,k1n)
c     call system_clock(after)
c     diapfl_time = real(after-before)/real(rate)
c     call sstbud(3,'   diapyc.mixing',temp(1,1,k1n))
c
ccc      write (string,'(a12,i8)') 'diapfl, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
      before = after
      call thermf(m,n,mm,nn,k1m,k1n)
      call system_clock(after)
      thermf_time = real(after-before)/real(rate)
c
ccc      write (string,'(a12,i8)') 'thermf, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
      before = after
      call eice(m,n,mm,nn,k1m,k1n)
      call system_clock(after)
      enloan_time = real(after-before)/real(rate)
c
      before = after
c     call mxlayr(m,n,mm,nn,k1m,k1n)
      if (nstep*baclin.gt.24.*3600) call mxkprf(m,n,mm,nn,k1m,k1n) !aft 24 hr
c
      call system_clock(after)
      mxlayr_time = real(after-before)/real(rate)
c     call sstbud(4,'  air-sea fluxes',tprime)
c     call sstbud(5,'     entrainment',temp(1,1,k1n))

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (diag_ape)
     . write (501,103) time,'  APE change due to thermal forcing:',
     .  hyc_pechg2(dp(1,1,k1n),th3d(1,1,k1n),31)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
ccc      write (string,'(a12,i8)') 'mxlayr, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
      before = after
      call hybgen(m,n,mm,nn,k1m,k1n)
      call system_clock(after)
      hybgen_time = real(after-before)/real(rate)
c
ccc      write (string,'(a12,i8)') 'hybgrd, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (diag_ape)
     . write (501,103) time,'  APE change due to vert.regridding:',
     .  hyc_pechg2(dp(1,1,k1n),th3d(1,1,k1n),31)
c     call sstbud(6,' vert.regridding',temp(1,1,k1n))
c
      if (dotrcr)
     .  write (lp,'(a)') 'tracer advection/turb.mixing cycle completed'
c
c ---------------------------------------------------------------------------
c
c --- output and diagnostic calculations
c
c ---------------------------------------------------------------------------
c
c --- after 1 day, check terms in time-integrated continuity eqn
c
      if (abs(time-tavini-1.).lt..001) then
      i=itest
      j=jtest
      jb=mod(j,jj)+1
      write (lp,'(i9,2i5,a/a)') nstep,i,j,
     . '  time-integrated continuity eqn diagnostics:',
     .  '     thknss_tndcy  horiz.flxdiv   vert.flxdiv      residuum'
      do k=1,kk
      km=k+mm
      kn=k+nn
      thkchg=dp(i,j,km)+dp(i,j,kn)-dpini(k)
      flxdiv=(uflxav(i+1,j,k)-uflxav(i,j,k)
     .       +vflxav(i,jb ,k)-vflxav(i,j,k))*scp2i(i,j)*delt1
      write (lp,102) k,thkchg,flxdiv,-diaflx(i,j,k),
     .  thkchg+flxdiv-diaflx(i,j,k)
 102  format (i3,4f14.1)
      end do
      end if
c
      if (diagno .or. mod(time+.0001,1.).lt..0002) then    !  once a day
c
c --- make line printer plot of mean layer thickness
        index=-1
        nflip=mod(nflip+1,2)
        do 705 k=1+(kk-1)*nflip,kk-(kk-1)*nflip,1-2*nflip
ccc     if (k.eq.kk-(kk-1)*nflip) index=1
        sum=0.
        totl=0.
        do 706 j=1,jj
        do 706 l=1,isp(j)
        do 706 i=ifp(j,l),ilp(j,l)
        totl=totl+oice(i,j)*scp2(i,j)
 706    sum=sum+dp(i,j,k+nn)*scp2(i,j)
        call linout(sum/(area*onem),charac(k),index)
 705    index=0
c --- add ice extend (%) to plot
        call linout(100.*totl/area,'I',1)
c
c --- diagnose mean sea surface height
        sum=0.
        sumice=0.
        do 704 j=1,jj
        do 704 l=1,isp(j)
        do 704 i=ifp(j,l),ilp(j,l)
c --- compute sea surface height (m)
        srfhgt(i,j)=(montg(i,j,1)+thref*pbavg(i,j,m))/g
        sumice=sumice+oice(i,j)*scp2(i,j)
 704    sum=sum+srfhgt(i,j)*scp2(i,j)
        write (lp,'(i9,'' mean sea srf.hgt. (mm):'',f9.2,f12.0)')nstep,
     .  sum*1.e3/area,sumice*1.e-6
c
c --- find the largest distance from a tencm layer from bottom
        do 707 j=1,jj
        do 707 l=1,isp(j)
        do 707 i=ifp(j,l),ilp(j,l)
        do 709 k=1,kk
        if (dp(i,j,k+nn).lt.tencm) then
          util2(i,j)=(p(i,j,kk+1)-p(i,j,k))/onem
          goto 708
        endif
 709    continue
 708    continue
 707    continue
c       call findmx(ip,util2,idm,ii,jj,'lowest')

c       call findmx(ipa,aflxa2o,iia,iia,jja,'aflxa2o')
c       call findmx(ip,osst,ii,ii,jj,'osst')
c       call findmx(ip,osss,ii,ii,jj,'osss')
c       call findmx(ip,usf,ii,ii,jj,'u_ocn')
c       call findmx(ip,vsf,ii,ii,jj,'v_ocn')
c
c       call findmx(ipa,asst,iia,iia,jja,'asst')
c       call findmx(ipa,sss,iia,iia,jja,'osss')
c       call findmx(ipa,uosurf,iia,iia,jja,'uosurf')
c       call findmx(ipa,vosurf,iia,iia,jja,'vosurf')
c
      end if                      ! once a day
c
      ogcm_time =
     .    cnuity_time
     .   +tsadvc_time
     .   +trcadv_time
     .   +momtum_time
     .   +barotp_time
     .   +thermf_time
     .   +enloan_time
     .   +mxlayr_time
     .   +hybgen_time

      write (lp,99009) ogcm_time,' sec for OGCM   at ocn step ',nstep
      write (lp,'(a/(5(4x,a,i5)))') 'timing (msec) by routine:'
     .   ,'cnuity',int(1000.*cnuity_time)
     .   ,'tsadvc',int(1000.*tsadvc_time)
     .   ,'trcadv',int(1000.*trcadv_time)
     .   ,'momtum',int(1000.*momtum_time)
     .   ,'barotp',int(1000.*barotp_time)
     .   ,'thermf',int(1000.*thermf_time)
     .   ,'enloan',int(1000.*enloan_time)
     .   ,'mxlayr',int(1000.*mxlayr_time)
     .   ,'hybgen',int(1000.*hybgen_time)
c
      if (mod(nstep,5).eq.0) call flush(lp)
c
      if (.not.diagno) go to 23
c
      write (lp,100) nstep,int((time+.001)/365.),mod(time+.001,365.)
 100  format (' ocn time step',i9,4x,'y e a r',i6,4x,'d a y',f9.2)
c
c --- output to history file
c
      call archiv(n,nn)
c --- diagnose meridional overturning and heat flux
c$OMP PARALLEL DO
        do 3 j=1,jj
        do 3 k=1,kk
        do 3 l=1,isp(j)
        do 3 i=ifp(j,l),ilp(j,l)
 3      p(i,j,k+1)=p(i,j,k)+dp(i,j,k+mm)
c$OMP END PARALLEL DO
css   call overtn(mm)
c
      write (lp,105) nstep
 105  format (' step',i9,' -- archiving completed --')
c
c --- test for cyclic-in-j vs. noncyclic domain
c --- (closed-basin coditions are indicated by a land barrier at j=jj)
c
      jj1=jj-1
      do i=1,ii1
        if (ip(i,jj).gt.0) then
          jj1=jj
        end if
      end do
c
c --- output to line printer
c
ccc      call prtmsk(ip,srfhgt,util3,idm,ii1,jj1,0.,1./(thref*onecm),
ccc     .     'sea surface height (cm)')
ccc      call prtmsk(ip,dpmixl(1,1),util3,idm,ii1,jj1,0.,1./onem,
ccc     .     'mixed layer depth (m)')
      call prtmsk(ip,temp(1,1,k1n),util3,idm,ii1,jj1,0.,10.,
     .     'mix.layer temp. (.1 deg)')
ccc      call prtmsk(ip,saln(1,1,k1n),util3,idm,ii1,jj1,35.,100.,
ccc     .     'mx.lay. salin. (.01 mil)')
      call prtmsk(ip,oice,util3,idm,ii1,jj1,0.,100.,
     .     'ice coverage (cm)')
        call prtmsk(ipa,asst,util3,iia,iia,jja,0.,1.,'asst ')
      do 77 j=1,jj1
      do 77 l=1,isu(j)
      do 77 i=ifu(j,l),ilu(j,l)
 77   util1(i,j)=u(i,j,k1n)+ubavg(i,j,n)
ccc      call prtmsk(iu,util1,util3,idm,ii1,jj1,0.,1000.,
ccc     .     'u vel. (mm/s), layer 1')
      do 78 i=1,ii1
      do 78 l=1,jsv(i)
      do 78 j=jfv(i,l),jlv(i,l)
 78   util2(i,j)=v(i,j,k1n)+vbavg(i,j,n)
ccc      call prtmsk(iv,util2,util3,idm,ii1,jj1,0.,1000.,
ccc     .     'v vel. (mm/s), layer 1')
ccc      call prtmsk(iu,ubavg(1,1,n),util3,idm,ii1,jj1,0.,1000.,
ccc     .     'barotrop. u vel. (mm/s)')
ccc      call prtmsk(iv,vbavg(1,1,n),util3,idm,ii1,jj1,0.,1000.,
ccc     .     'barotrop. v vel. (mm/s)')
 23   continue

c --- accumulate fields for agcm
c$OMP PARALLEL DO
      do 201 j=1,jj
      do 201 l=1,isp(j)
      do 201 i=ifp(j,l),ilp(j,l)
      osst(i,j)=osst(i,j)+temp(i,j,k1n)*baclin/(3600.*real(nhr))
      osss(i,j)=osss(i,j)+saln(i,j,k1n)*baclin/(3600.*real(nhr))
      osiav(i,j)=osiav(i,j)+odmsi(i,j)*baclin*dtsrc/(3600.*real(nhr)) !kg/m2=>kg*.5*hr/m2
      omlhc(i,j)=spcifh*max(dp(i,j,k1n)/onem,thkmin)/thref  ! J/(m2*C)
      oogeoza(i,j)=(montg(i,j,1)+thref*pbavg(i,j,m))*g/(thref*onem) ! m^2/s^2
 201  continue
c$OMP END PARALLEL DO
c
      nsaveo=nsaveo+1
 15   continue
      if (nsaveo*baclin.ne.nhr*3600) then
            print *, ' ogcm saved over hr=',nsaveo*baclin/3600.
            stop ' stop: ogcm saved over hr'
      end if
      nsaveo=0
c
c$OMP PARALLEL DO
      do 87 j=1,jj
      do 86 i=1,ii
      usf(i,j)=0.
      vsf(i,j)=0.
 86   continue
      do 87 l=1,isu(j)
      do 87 i=ifu(j,l),ilu(j,l)
 87   usf(i,j)=u(i,j,k1n)
c$OMP END PARALLEL DO
c
      do 88 i=1,ii
      do 88 l=1,jsv(i)
      do 88 j=jfv(i,l),jlv(i,l)
 88   vsf(i,j)=v(i,j,k1n)
c
      call ssto2a(osst,asst)
      call ssto2a(osss,sss)
css   call iceo2a(omlhc,mlhc)
      call ssto2a(oogeoza,ogeoza)
      call ssto2a(osiav,utila)                 !kg/m*m per agcm time step
      call veco2a(usf,vsf,uosurf,vosurf)
c
c       call findmx(ip,osst,ii,ii,jj,'osst')
c       call findmx(ip,osss,ii,ii,jj,'osss')
c       call findmx(ip,usf,ii,ii,jj,'u_ocn')
c       call findmx(ip,vsf,ii,ii,jj,'v_ocn')
c
c       call findmx(ipa,asst,iia,iia,jja,'asst')
c       call findmx(ipa,sss,iia,iia,jja,'osss')
c       call findmx(ipa,uosurf,iia,iia,jja,'uosurf')
c       call findmx(ipa,vosurf,iia,iia,jja,'vosurf')
c
c$OMP PARALLEL DO PRIVATE(tf)
      do 204 ia=1,iia
      do 204 ja=1,jja
      if (focean(ia,ja).gt.0.) then
        gtemp(1,1,ia,ja)=asst(ia,ja)
        tf=tfrez(sss(ia,ja),0.)
        dmsi(1,ia,ja)=utila(ia,ja)                        !kg/m2 per agcm step
        dhsi(1,ia,ja)=utila(ia,ja)                        !J/m2 per agcm step
     .        *(tf*shi-lhm*(1-0.001*fsss*sss(ia,ja)))
        dssi(1,ia,ja)=1.d-3*dmsi(1,ia,ja)*sss(ia,ja)*fsss !kg/m2 per agcm step
c --- evenly distribute new ice over open water and sea ice
        if (aice(ia,ja).gt.1.e-3) then
          dhsi(2,ia,ja)=dhsi(1,ia,ja)
          dmsi(2,ia,ja)=dmsi(1,ia,ja)
          dssi(2,ia,ja)=dssi(1,ia,ja)
        endif
      endif
 204  continue
c$OMP END PARALLEL DO
      call system_clock(afogcm)
      if (abs(time-(itime+1.)/nday).gt..01) then
        write(*,'(a,2f10.3)') 'stop: agcm/ogcm date after ogcm'
     .    ,(itime+1.)/nday,time
        stop ' stop: dates in agcm/ogcm do not match'
      end if
c
      return
      end
c
c> Revision history:
c>
c> May  1997 - removed statement "theta(1)=-thbase" after loop 14
c> June 1997 - added loop 60 to fix bug in vertical summation of -diaflx-
c> Mar. 2000 - conversion to SI units
c> June 2000 - added timing logic
c> Aug. 2000 - added code for changing number of OpenMP threads
c> Apr. 2001 - eliminated stmt_funcs.h
c> June 2001 - corrected sign error in diaflx
c> July 2001 - replaced archiving statements by 'call archiv'
c> Oct  2004 - map ice mass to agcm, and then calculate E
