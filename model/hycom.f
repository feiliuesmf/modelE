#include "hycom_mpi_hacks.h"
#include "rundeck_opts.h"
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
      USE DOMAIN_DECOMP, only: AM_I_ROOT, HALO_UPDATE, NORTH,
     &                         haveLatitude, GLOBALSUM
      USE HYCOM_ATM !, only : gather_atm, scatter_atm
!      USE FLUXES, only : e0,prec,eprec,evapor,flowo,eflowo,dmua,dmva
!     . ,erunosi,runosi,srunosi,runpsi,srunpsi,dmui,dmvi,dmsi,dhsi,dssi
!     . ,gtemp,sss,mlhc,ogeoza,uosurf,vosurf,MELTI,EMELTI,SMELTI
!     . ,gmelt,egmelt,solar,gtempr
#ifdef TRACERS_GASEXCH_Natassa
      USE FLUXES, only : TRGASEX,GTRACER

      USE TRACER_COM, only : ntm    !tracers involved in air-sea gas exch

      USE TRACER_GASEXCH_COM, only : atracflx,atrac,tracflx
#endif

      !USE SEAICE_COM, only : rsi,msi
      USE SEAICE, only : fsss,tfrez ! number, function - ok
      USE GEOM, only : dxyp ! ok
      !USE MODEL_COM, only : focean
      USE CONSTANT, only : lhm,shi,shw
      USE MODEL_COM, only: dtsrc
     *  ,itime,iyear1,nday,jdendofm,jyear,jmon,jday,jdate,jhour,aMON
#ifdef TRACERS_OceanBiology 
      USE obio_dim
      USE obio_forc, only:    avgq,awind,owind,asolz,osolz
      USE obio_com,  only:    gcmax,pCO2,dobio

      USE PBLCOM, only : wsavg 
      USE RAD_COM,   only: COSZ1
     .           ,FSRDIR,SRVISSURF,FSRDIF,DIRNIR,DIFNIR
#endif 
#ifdef OBIO_RAD_coupling
      USE obio_forc, only:    avisdir,avisdif,anirdir,anirdif
     .                       ,ovisdir,ovisdif,onirdir,onirdif
#ifdef CHL_from_OBIO
      USE FLUXES, only: chl
      USE obio_com, only: tot_chlo
#endif
#endif
      USE HYCOM_DIM_GLOB
      USE HYCOM_DIM, only : aJ_0, aJ_1, aJ_0H, aJ_1H,
     &                       J_0,  J_1,  J_0H,  J_1H,
     &      ogrid, isp_loc => isp, ifp_loc => ifp, ilp_loc => ilp
      USE HYCOM_SCALARS
      USE HYCOM_ARRAYS_GLOB
      USE hycom_arrays_glob_renamer
c
      USE KPRF_ARRAYS
      USE HYCOM_CPLER
      implicit none
c
!!#include "dimensions.h"
#include "dimension2.h"
!!#include "common_blocks.h"
!!#include "cpl.h"
!!! afogcm,nsavea should be initialized properly !
      integer :: afogcm=0,nsavea=0,nsaveo
!!#include "a2o.h"
#include "kprf_scalars.h"
!!#include "kprf_arrays.h"
c
      real sum,coord,x,x1,totl,sumice,fusion,saldif,sofsig,tf
     .    ,sigocn,kappaf,chk_rho,chk_kap,apehyc,pechg_hyc_bolus
     .    ,hyc_pechg1,hyc_pechg2,q,sum1,sum2,dpini(kdm)
     .    ,thkchg,flxdiv,den_str
      real totlj(J_0H:J_1H), sumj(J_0H:J_1H), sumicej(J_0H:J_1H)
      logical doThis
      integer jj1,no,index,nflip,mo0,mo1,mo2,mo3,rename,iatest,jatest
     .       ,OMP_GET_NUM_THREADS,io,jo,ipa(iia,jja),nsub
      integer ipa_loc(iia,aJ_0H:aJ_1H)
#ifdef TRACERS_GASEXCH_Natassa
      integer nt
#endif
#ifdef TRACERS_OceanBiology
      integer ihr,ichan,hour_of_day,day_of_month,iyear
      integer bef,aft                   !  bio routine timing variables
#endif
      external rename
      logical master,slave,diag_ape
      character util(idm*jdm+14)*2,charac(20)*1,string*20,
     .          flnm*60
      data charac/'1','2','3','4','5','6','7','8','9','0',
     .            'A','B','C','D','E','F','G','H','I','J'/
c     data iatest,jatest/26,7/            ! Weddell Sea
c     data iatest,jatest/33,40/           ! Iceland
      data iatest,jatest/31,41/           ! Iceland
      data nflip/0/
c
      real cnuity_time,tsadvc_time,momtum_time,barotp_time,trcadv_time,
     .     thermf_time,enloan_time,mxlayr_time,hybgen_time,
     .     agcm_time,ogcm_time
#ifdef TRACERS_OceanBiology
      real*4 ocnbio_time,ocnbio_total_time,ocnbio_avg_time
#endif
      integer after,before,rate,bfogcm
c
      real osst(idm,jdm),osss(idm,jdm),osiav(idm,jdm)
     . ,oogeoza(idm,jdm),usf(idm,jdm),vsf(idm,jdm)
     . ,utila(iia,jja)
      real osst_loc(idm,J_0H:J_1H),osss_loc(idm,J_0H:J_1H),
     &    osiav_loc(idm,J_0H:J_1H), oogeoza_loc(idm,J_0H:J_1H)
#ifdef TRACERS_GASEXCH_Natassa
     . ,otrac(idm,jdm,ntm)
#endif

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
      ! move to global atm grid

      call getdte(Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)

cdiag write(*,'(a,i8,7i5,a)')'chk =',
cdiag.    Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon
c
      if (mod(jhour,nhr).eq.0.and.mod(itime,nday/24).eq.0) then
c$OMP PARALLEL DO SCHEDULE(STATIC,jchunk)
        do 28 ia=1,iia
        do 28 ja=aJ_0,aJ_1
          ataux_loc(ia,ja)=0.
          atauy_loc(ia,ja)=0.
        aflxa2o_loc(ia,ja)=0.
          aemnp_loc(ia,ja)=0.
          asalt_loc(ia,ja)=0.
           aice_loc(ia,ja)=0.
         austar_loc(ia,ja)=0.
         aswflx_loc(ia,ja)=0.
#ifdef TRACERS_GASEXCH_Natassa
        do nt=1,ntm
        atracflx(ia,ja,nt)=0.
        enddo
#endif
#ifdef TRACERS_OceanBiology
          awind(ia,ja)=0.
          asolz(ia,ja)=0.
#endif
#ifdef OBIO_RAD_coupling
          avisdir(ia,ja)=0.
          avisdif(ia,ja)=0.
          anirdir(ia,ja)=0.
          anirdif(ia,ja)=0.
#endif
 28     continue
c$OMP END PARALLEL DO
      endif
c
c$OMP PARALLEL DO SCHEDULE(STATIC,jchunk)
c --- accumulate agcm fields over nhr 
      do 29 ia=1,iia
      do 29 ja=aJ_0,aJ_1
      ipa_loc(ia,ja)=0
      if (focean_loc(ia,ja).eq.0.) goto 29
      ipa_loc(ia,ja)=1
c --- accumulate
      aemnp_loc(ia,ja)=aemnp_loc(ia,ja)                               ! kg/m*m => m/s
     .+((prec_loc(ia,ja)-evapor_loc(ia,ja,1))*(1.-rsi_loc(ia,ja))         ! open water
     .+(flowo_loc(ia,ja)+gmelt_loc(ia,ja)+melti_loc(ia,ja))/
     .                       (dxyp(ja)*focean_loc(ia,ja))!ocn/ice
     .+(runosi_loc(ia,ja)+runpsi_loc(ia,ja))*rsi_loc(ia,ja))*thref        ! ice
     .                                /(3600.*real(nhr))
      aflxa2o_loc(ia,ja)=aflxa2o_loc(ia,ja)                           ! J/m*m => W/m*m
     . +((e0_loc(ia,ja,1)+eprec_loc(ia,ja))*(1.-rsi_loc(ia,ja))           ! ocean water
     . +(eflowo_loc(ia,ja)+egmelt_loc(ia,ja)+emelti_loc(ia,ja))
     .                              /(dxyp(ja)*focean_loc(ia,ja)) ! ocn or ice
css  . +(erunosi(ia,ja)+erunpsi(ia,ja))*rsi(ia,ja))           ! ice
     . + erunosi_loc(ia,ja)*rsi_loc(ia,ja))                           ! ice
     .                                 /(3600.*real(nhr))
      asalt_loc(ia,ja)=asalt_loc(ia,ja)                            ! kg/m*m/sec salt
     .+((srunosi_loc(ia,ja)+srunpsi_loc(ia,ja))*rsi_loc(ia,ja)         !
     .   +smelti_loc(ia,ja)/(dxyp(ja)*focean_loc(ia,ja)))          !
     .                             /(3600.*real(nhr))
       aice_loc(ia,ja)= aice_loc(ia,ja) + rsi_loc(ia,ja)*
     .                                  dtsrc/(real(nhr)*3600.)
c --- dmua on B-grid, dmui on C-grid; Nick aug04
      ataux_loc(ia,ja)=ataux_loc(ia,ja)+(dmua_loc(ia,ja,1)+
     .                                         dmui_loc(ia,ja))    ! scaled by rsi 
     .                                     /(3600.*real(nhr))  ! kg/ms => N/m*m
      atauy_loc(ia,ja)=atauy_loc(ia,ja)+(dmva_loc(ia,ja,1)+
     .                                         dmvi_loc(ia,ja))
     .                                     /(3600.*real(nhr))  ! kg/ms => N/m*m
      austar_loc(ia,ja)=austar_loc(ia,ja)+(
     . sqrt(sqrt((dmua_loc(ia,ja,1)+dmui_loc(ia,ja))**2
     .          +(dmva_loc(ia,ja,1)+dmvi_loc(ia,ja))**2)/dtsrc*thref)) ! sqrt(T/r)=>m/s
     .                               *dtsrc/(real(nhr)*3600.)
      aswflx_loc(ia,ja)=aswflx_loc(ia,ja)+(solar_loc(1,ia,ja)*
     .                               (1.-rsi_loc(ia,ja))!J/m*m=>W/m*m
     .                         +solar_loc(3,ia,ja)*    rsi_loc(ia,ja))
     .                                 /(3600.*real(nhr))
#ifdef TRACERS_GASEXCH_Natassa
      do nt=1,ntm
      atracflx(ia,ja,nt)= atracflx(ia,ja,nt)
     .                 + TRGASEX(nt,1,ia,ja)      ! in mol/m2/s
     .                 * dtsrc/(real(nhr)*3600.)
      enddo
#endif
#ifdef TRACERS_OceanBiology
      asolz(ia,ja)=asolz(ia,ja)                            !
     .            +COSZ1(ia,ja)*dtsrc/(3600.*real(nhr))    !
      awind(ia,ja)=awind(ia,ja)                            !
     .            +wsavg(ia,ja)*dtsrc/(3600.*real(nhr))    !
#endif
#ifdef OBIO_RAD_coupling
      avisdir(ia,ja)=avisdir(ia,ja)                        !
     . +FSRDIR(ia,ja)*SRVISSURF(ia,ja)*dtsrc/(3600.*real(nhr))    !
      avisdif(ia,ja)=avisdif(ia,ja)                        !
     .            +FSRDIF(ia,ja)*dtsrc/(3600.*real(nhr))    !
      anirdir(ia,ja)=anirdir(ia,ja)                        !
     .            +DIRNIR(ia,ja)*dtsrc/(3600.*real(nhr))    !
      anirdif(ia,ja)=anirdif(ia,ja)                        !
     .            +DIFNIR(ia,ja)*dtsrc/(3600.*real(nhr))    !

!     write(*,*)'hycom, fsrdir:',nstep,ia,ja,
!    .     FSRDIR(ia,ja),SRVISSURF(ia,ja),avisdir(ia,ja)

#endif
 29   continue
c$OMP END PARALLEL DO
c
c     ia=iatest
c     ja=jatest
c     write(402,'(a,2i3/13es12.3 )') 'emnp=',ia,ja
c    .,aemnp(ia,ja)
c    .,(prec(ia,ja)-evapor(ia,ja,1))*(1.-rsi(ia,ja))         ! open water
c    .,(flowo(ia,ja)+gmelt(ia,ja)+melti(ia,ja))/(dxyp(ja)*focean(ia,ja))!ocn/ice
c    .,(runosi(ia,ja)+runpsi(ia,ja))*rsi(ia,ja)
c    .  ,prec(ia,ja),evapor(ia,ja,1),rsi(ia,ja)
c    .  ,flowo(ia,ja),gmelt(ia,ja),melti(ia,ja),dxyp(ja),focean(ia,ja)
c    .  ,runosi(ia,ja),runpsi(ia,ja)
c
      call gather_atm
      call gather2_atm

      nsavea=nsavea+1

      if (mod(jhour,nhr).gt.0.or.mod(itime,nday/24).eq.0) goto 9666

      if (nsavea*24/nday.ne.nhr) then
        write(*,'(a,4i4,i8)')
     .  'nonmatching b.c. accumulation periods: agcm/ogcm:',
     .  nsavea*24/nday,nhr,nsavea,nday,itime
        stop 'agcm/ogcm accumulating periods differ'
      end if
      nsavea=0
c
      if (AM_I_ROOT()) then

      call veca2o(ataux,atauy,taux,tauy)          !wind stress
      call flxa2o(aice,oice)                      !ice coverage
      call flxa2o(aflxa2o,oflxa2o)                !heatflux everywhere
      call flxa2o(asalt,osalt)                    !saltflux from SI
      call flxa2o(aemnp,oemnp)                    !E - P everywhere
      call flxa2o(austar,ustar)                   !friction velocity
      call flxa2o(aswflx,sswflx)                  ! shortwave flux
#ifdef TRACERS_GASEXCH_Natassa
      do nt=1,ntm
      call flxa2o(atracflx(:,:,nt),tracflx(:,:,nt)) !tracer flux
      enddo
      do j=1,jj
      do l=1,isp(j)
      do i=ifp(j,l),ilp(j,l)
!     write(*,'(a,3i5,e12.4)')'hycom, TRGASEX2: ',
!    .  nstep,i,j,tracflx(i,j,ntm)
      enddo
      enddo
      enddo
#endif
#ifdef TRACERS_OceanBiology
      call flxa2o(asolz,osolz)
      call flxa2o(awind,owind)
#endif
#ifdef OBIO_RAD_coupling
      call flxa2o(avisdir,ovisdir)
      call flxa2o(avisdif,ovisdif)
      call flxa2o(anirdir,onirdir)
      call flxa2o(anirdif,onirdif)
#endif

      endif ! AM_I_ROOT

      call scatter_kprf_arrays
      call scatter_hycom_arrays
c
      call system_clock(before)
      call system_clock(count_rate=rate)
      bfogcm=before
      agcm_time = real(bfogcm-afogcm)/real(rate)

      if (AM_I_ROOT()) then
      write (lp,99009) agcm_time,' sec for AGCM bfor ocn step ',nstep
99009 format (f12.3,a,i8)
      endif ! AM_I_ROOT
c
!hack hack !!! copied from inicon:
cddd      if (nstep0.eq.0) then
cddd        delt1=baclin
cddd      else
cddd        delt1=baclin+baclin
cddd      endif

      if (nstep.eq.0 .or. nstep.eq.nstep0) then

c
c$OMP PARALLEL SHARED(mo0)
c$    mo0=OMP_GET_NUM_THREADS()
c$OMP END PARALLEL
c$    call OMP_SET_DYNAMIC(.false.)
c$    write (lp,'(2(a,i5))') ' hycom thread count',mo0
c
ccc$  call OMP_SET_NUM_THREADS(max(mo0,16))
ccc$  OMP PARALLEL SHARED(mo1)
ccc$  mo1=OMP_GET_NUM_THREADS()
ccc$  OMP END PARALLEL
ccc$  write (lp,'(2(a,i5))') ' hycom thread count',mo0,' changed to',mo1
ccc$     . ,' =======  number of threads:',mo1,' ======='
      if (AM_I_ROOT())
     &  write (lp,'(2(a,i5))') ' hycom thread count:',mo0
      jchunk=50
c     if (mo0.eq.4) jchunk=50           !  jdm=180
c     if (mo0.eq.6) jchunk=32           !  jdm=180
c     if (mo0.eq.8) jchunk=23           !  jdm=180
      if (jchunk.lt.0) then
        if (AM_I_ROOT()) write (lp,*) 'you forgot to define jchunk'
        stop   ! proper mpi_abort
      end if
c
      lp=6
c
c --- compute eqn.of state check values
c
      if (pref.eq.2.e7) then
        chk_rho=36.719718               ! T:[-2:30],S:[18:38]
        chk_rho=36.876506               ! T:[-2:32],S:[16:38]
c ---   chk_kap: t= 1.0, s=34.5, p=0 bar, kap_t(4.5,34.5,1.e7)= 0.08728276
        chk_kap=0.08728276
      elseif (pref.eq.0.) then
        chk_rho=27.786223
c ---   chk_kap: t= 1.0, s=34.5, p=0 pascal, kap_t(4.5,34.5,1.e7)=-0.09080282
        chk_kap=-0.09080282
      else
        stop 'wrong pref'    ! proper mpi_abort
      endif
c
      if (abs(sigocn(4.,35.)-chk_rho).gt..0001) then
      if (AM_I_ROOT())
     & write (lp,'(/2(a,f11.6))') 'error -- sigocn(t=4,s=35) should be',
     . chk_rho,', not',sigocn(4.,35.)
css     stop
      end if
      if (abs(kappaf(4.5,34.5,1.e7)-chk_kap).gt..00001) then
      if (AM_I_ROOT())
     &  write (lp,'(/a,2(a,f12.8))') 'error: kappa(4.5,34.5,10^7)',
     .  '  should be',chk_kap,', not',kappaf(4.5,34.5,1.e7)
      stop
      end if
c
      if (AM_I_ROOT())
     & write (lp,109) thkdff,temdff,veldff,viscos,diapyc,vertmx
 109  format (' turb. flux parameters:',1p/
     .  ' thkdff,temdff,veldff =',3e9.2/
     .  ' viscos,diapyc,vertmx =',3e9.2)
c
c --- 'lstep' = number of barotropic time steps per baroclinic time step.
c --- lstep   m u s t   be even.
c
      lstep=baclin/batrop
      lstep=2*((lstep+1)/2)
      dlt=baclin/lstep
      nstepi=real(nhr)*3600./baclin + .0001

      if (AM_I_ROOT()) then
      write (lp,'(i4,'' barotropic steps per baroclinic time step'')')
     .  lstep
      write (lp,*) "time0/nstep0=",time0,nstep0
      write (lp,'(''ogcm exchange w. agcm every step,hr'',2i5)')
     .  nstepi,nhr
      endif ! AM_I_ROOT
c
c --- set up parameters defining the geographic environment
c
css   call geopar                ! moved to agcm for zero start or below
css   call inicon                ! moved to agcm
      if (mxlkpp) then
        call inikpp                 ! kpp
        if (mxlgis) stop 'wrong kprf: kpp=gis=true'
        print *,'chk: mxlkpp=true'
        write(*,'(a,2f9.3)')' chk qkpar=',
     .  1./akpar(itest,jtest,1),1./(betabl(2)*onem)
      elseif (mxlgis) then
        call inigis                 ! giss
        print *,'chk: mxlgis=true'
      else
        stop 'wrong kprf: kpp=gis=false'
      endif
c
      mixfrq=int(43200./baclin+.0001)
      trcfrq=int(43200./baclin+.0001)
      if (AM_I_ROOT())
     &  write (lp,'(a,2i4)') 'trcfrq set to',trcfrq

      if (trcout) dotrcr=.true.
c
      watcum=0.
      empcum=0.
c
      if (jyear-iyear1.le.20) then
        if (AM_I_ROOT())
     .  write (lp,'(''starting date in ogcm/agcm '',2i5,'' hr '',2i12
     .   /'' iyear1/jyear='',2i5)')
     .  int((nstep0+nstepi-1)*baclin)/(3600*24),itime/nday
     . ,int((nstep0+nstepi-1)*baclin)/3600,itime*24/nday,iyear1,jyear
      else
        if (AM_I_ROOT())
     .  write (lp,'(''starting date in ogcm/agcm '',2i12
     .   /'' iyear1/jyear='',2i5)') 
     .  int((nstep0+nstepi-1)*baclin)/(3600*24),itime/nday
     . ,iyear1,jyear,(nstep0+nstepi-1)*baclin/3600.,itime*24/nday
      endif
c
c     if(abs((nstep0+nstepi-1)*baclin/3600.-itime*24./nday).gt.1.e-5)
      if(abs(int((nstep0+nstepi-1)*baclin/3600)-itime*24/nday).gt.0)
     .                                                        then
        if (AM_I_ROOT()) then
        write (lp,'(a,f16.8,i10,f16.8)')'mismatch date found '
     . ,int((nstep0+nstepi-1)*baclin/3600),itime*24/nday
     . ,(nstep0+nstepi-1)*baclin/3600.-itime*24/nday
        write (lp,'(/(a,i9,a,i9,a,f7.1,a))')'chk model start step'
     &       ,nstep0
     . ,' changed to: '
     . ,int((itime*24/nday+(iyear1-jyear)*8760)*3600/baclin)-nstepi+1
     . ,' (day '
     .,(int((itime*24/nday+(iyear1-jyear)*8760)*3600/baclin)-nstepi+1)
     .  *baclin/86400,')'

        end if   ! AM_I_ROOT

        nstep0=int((itime*24/nday+(iyear1-jyear)*8760)*3600/baclin)
     .                                         -nstepi+3600/baclin

        nstep=nstep0
        time0=nstep0*baclin/86400
        tavini=time0
c
! dpini is used to compute thkchg which is used only for printing
      if (haveLatitude(ogrid, J=jtest)) then
        do 19 k=1,kk
 19     dpini(k)=dp_loc(itest,jtest,k)+dp_loc(itest,jtest,k+kk)
      end if
c
      else
        write (lp,'(/(a,i9))') 'chk model starts at steps',nstep
      endif

c
      endif       ! if (nstep=0 or nstep=nstep0) 
c
      if (nstepi.le.0) stop 'wrong nstepi'
c
c     print *,' shown below oice'
c     call zebra(oice,iio,iio,jjo)
c     print *,' shown below aflxa2o'
c     call zebra(aflxa2o,iia,iia,jja)
c
      osst=0.  ;   osst_loc=0.;
      osss=0.  ;   osss_loc=0.;
      osiav=0. ;  osiav_loc=0.;

      CALL HALO_UPDATE(ogrid,corio_loc,     FROM=NORTH)
      hekman_loc=0.;
c
c$OMP PARALLEL DO PRIVATE(jb)
      do 202 j=J_0,J_1
      jb = PERIODIC_INDEX(j+1, jj)
      do 202 l=1,isp_loc(j)
      do 202 i=ifp_loc(j,l),ilp_loc(j,l)
c
      hekman_loc(i,j)=ustar_loc(i,j)*(cekman*4.0)/                ! kpp
     &           (abs(corio_loc(i,j  ))+abs(corio_loc(i+1,j  ))+
     &            abs(corio_loc(i,jb ))+abs(corio_loc(i+1,jb )))
#ifdef TRACERS_GASEXCH_Natassa
      do nt=1,ntm
      otrac(i,j,nt)=0.
      enddo
#endif
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
c     open(202,file='ip.array',form='unformatted',status='unknown')
c     write(202) ip,iu,iv,ifp,ilp,isp
c     close(202)
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
      if (nstep.eq.1) diagno=.true.
      diag_ape=.false.
c
      trcadv_time = 0.0
!---------------------below not implemented------------------
      if (dotrcr) then
        call system_clock(before)      ! time elapsed since last system_clock

c --- initialization of tracer transport arrays (incl. dpinit):
        call gather_hycom_arrays
        if (AM_I_ROOT()) then
          call tradv0(m,mm)
        endif
        call scatter_hycom_arrays
cdiag print *,'past tradv0'
        dotrcr=.false.
c
c --- inject tracer into ogcm
c     write(*,*) 'nstep=',nstep
c     call findmx(ip,temp(1,1,k1n),idm,ii,jj,'sst')
c     call findmx(ip,saln(1,1,k1n),idm,ii,jj,'sss')
c     if (trcout) call cfcflx(tracer,p,temp(1,1,k1n),saln(1,1,k1n)
c    .           ,latij(1,1,3),scp2,baclin*trcfrq)
#if defined(TRACERS_HYCOM_Ventilation)
       call gather_hycom_arrays
       if (AM_I_ROOT()) then
c$OMP PARALLEL DO
       do 12 j=1,jj
       do 12 l=1,isp(j)
       do 12 i=ifp(j,l),ilp(j,l)
12     tracer(i,j,1,1)=1.              !  surface ventilation tracer
c$OMP END PARALLEL DO
      endif
      call scatter_hycom_arrays
#endif
        call system_clock(after)        ! time elapsed since last system_clock
        trcadv_time = real(after-before)/real(rate)
      end if ! dotrcr  above if block is done only once

c
#ifdef TRACERS_OceanBiology
      call system_clock(before)      ! time elapsed since last system_clock
      trcout = .true.

      before = after
      ocnbio_time=0.0


      if (dobio) then
cdiag  do k=1,kdm
cdiag  km=k+mm
cdiag   if (k.eq.1)
cdiag.   write(lp,'(a,3i5,a,a)')'bfre obio ',nstep,itest,jtest
cdiag.   ,'    dp       nitr      ammo      sili      iron  '
cdiag.   ,'    diat      chlo      cyan      cocco     herb'
cdiag   write(lp,'(22x,i3,10(1x,es9.2))')
cdiag.     k,dp(itest,jtest,km)/onem,
cdiag.     tracer(itest,jtest,k,1),tracer(itest,jtest,k,2),
cdiag.     tracer(itest,jtest,k,3),tracer(itest,jtest,k,4),
cdiag.     tracer(itest,jtest,k,5),tracer(itest,jtest,k,6),
cdiag.     tracer(itest,jtest,k,7),tracer(itest,jtest,k,8),
cdiag.     tracer(itest,jtest,k,9)
cdiag  enddo
cdiag  call obio_limits('bfre obio_model')

         call obio_model(mm)

cdiag  do k=1,kdm
cdiag  km=k+mm
cdiag   if (k.eq.1)
cdiag.   write(lp,'(a,3i5,a,a)')'aftr obio ',nstep,itest,jtest
cdiag.   ,'    dp       nitr      ammo      sili      iron  '
cdiag.   ,'    diat      chlo      cyan      cocco     herb'
cdiag   write(lp,'(22x,i3,10(1x,es9.2))')
cdiag.     k,dp(itest,jtest,km)/onem,
cdiag.     tracer(itest,jtest,k,1),tracer(itest,jtest,k,2),
cdiag.     tracer(itest,jtest,k,3),tracer(itest,jtest,k,4),
cdiag.     tracer(itest,jtest,k,5),tracer(itest,jtest,k,6),
cdiag.     tracer(itest,jtest,k,7),tracer(itest,jtest,k,8),
cdiag.     tracer(itest,jtest,k,9)
cdiag  enddo
cdiag  call obio_limits('aftr obio_model')

      endif


      call system_clock(after)
      ocnbio_time = real(after-before)/real(rate)

#endif
!---------------------above not implemented------------------

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- available potential energy diagnostics:
      if (diag_ape) then
        call gather_hycom_arrays  !mkb check if all arrays are in
        if (AM_I_ROOT())
     .    write (501,'(f9.1,a,1p,e15.7)') time,'  APE (J):',
     .    hyc_pechg1(dp(1,1,k1n),th3d(1,1,k1n),31)
      end if ! diag_ape
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 103  format (f9.1,a,-12p,f9.3,' TW')
c
cddd!!! just checking...
cddd      call gather_hycom_arrays
cddd
cddd      if (AM_I_ROOT()) then ! work on global grids here

      call system_clock(before)

      call cnuity(m,n,mm,nn,k1m,k1n)

      call gather_hycom_arrays
      if (AM_I_ROOT()) then ! work on global grids here
        call hycom_arrays_checksum
      end if ! AM_I_ROOT

      call system_clock(after)
      cnuity_time = real(after-before)/real(rate)

c     call sstbud(0,'  initialization',temp(1,1,k1n))
c
ccc      write (string,'(a12,i8)') 'cnuity, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
!---------------------below not implemented------------------
      if (trcout) then
        before = after
c --- long time step tracer advection: build up mass flux time integral

        if (n.eq.oddev) then
          call gather_hycom_arrays
          if (AM_I_ROOT()) then
            call tradv1(n,nn)
          endif
          call scatter_hycom_arrays

        endif
cdiag print *,'past tradv1'
c
        if (mod(nstep,trcfrq).eq.0) then
          dotrcr=.true.         !  start tracer advection/turb.mixing cycle
          write (lp,'(a)') 'start tracer advection/turb.mixing cycle'
          before = after
c --- tracer transport:
          call gather_hycom_arrays
          if (AM_I_ROOT()) then
            call tradv2(n,nn)
          endif
          call scatter_hycom_arrays

        end if
        trcadv_time = trcadv_time + real(after-before)/real(rate)

#ifdef TRACERS_OceanBiology
cdiag  do k=1,kdm
cdiag  km=k+mm
cdiag   if (k.eq.1)
cdiag.   write(lp,'(a,3i5,a,a)')'aftr tadv ',nstep,itest,jtest
cdiag.   ,'    dp       nitr      ammo      sili      iron  '
cdiag.   ,'    diat      chlo      cyan      cocco     herb'
cdiag   write(lp,'(22x,i3,10(1x,es9.2))')
cdiag.     k,dp(itest,jtest,km)/onem,
cdiag.     tracer(itest,jtest,k,1),tracer(itest,jtest,k,2),
cdiag.     tracer(itest,jtest,k,3),tracer(itest,jtest,k,4),
cdiag.     tracer(itest,jtest,k,5),tracer(itest,jtest,k,6),
cdiag.     tracer(itest,jtest,k,7),tracer(itest,jtest,k,8),
cdiag.     tracer(itest,jtest,k,9)
cdiag  enddo
cdiag  call obio_limits('aftr trcadv')
#endif

      end if !trcout
!---------------------above not implemented------------------
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (diag_ape) then
        call gather_hycom_arrays
        if (AM_I_ROOT())
     .    q=hyc_pechg1(dp(1,1,k1m),th3d(1,1,k1m),32)
      endif  ! diag_ape
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      before = after

      call tsadvc(m,n,mm,nn,k1m,k1n)

      call system_clock(after)
      tsadvc_time = real(after-before)/real(rate)
c
c     call sstbud(1,' horiz.advection',temp(1,1,k1n))
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (diag_ape) then
        call gather_hycom_arrays
        if (AM_I_ROOT())
     .    write (501,103) time,'  APE change due to time smoothing: ',
     .    hyc_pechg2(dp(1,1,k1m),th3d(1,1,k1m),32)
      end if ! diag_ape
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
      if (diag_ape) then
        call gather_hycom_arrays
        if (AM_I_ROOT())
     .  write (501,103) time,'  APE change due to mass adjustment:',
     .  hyc_pechg2(dp(1,1,k1n),th3d(1,1,k1n),31)
      end if ! diag_ape
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
ccc      write (string,'(a12,i8)') 'convec, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
c     before = after

      call diapfl(m,n,mm,nn,k1m,k1n)

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

      if (nstep*baclin.gt.12.*3600) then
         call mxkprf(m,n,mm,nn,k1m,k1n) !aft 12 hr
      end if
c
      call system_clock(after)
      mxlayr_time = real(after-before)/real(rate)

c     call sstbud(4,'  air-sea fluxes',tprime)
c     call sstbud(5,'     entrainment',temp(1,1,k1n))

#ifdef TRACERS_OceanBiology
cdiag  do k=1,kdm
cdiag  km=k+mm
cdiag  if (k.eq.1)
cdiag.   write(lp,'(a,3i5,a,a)')'aftr kpp  ',nstep,itest,jtest
cdiag.   ,'    dp       nitr      ammo      sili      iron  '
cdiag.   ,'    diat      chlo      cyan      cocco     herb'
cdiag   write(lp,'(22x,i3,10(1x,es9.2))')
cdiag.     k,dp(itest,jtest,km)/onem,
cdiag.     tracer(itest,jtest,k,1),tracer(itest,jtest,k,2),
cdiag.     tracer(itest,jtest,k,3),tracer(itest,jtest,k,4),
cdiag.     tracer(itest,jtest,k,5),tracer(itest,jtest,k,6),
cdiag.     tracer(itest,jtest,k,7),tracer(itest,jtest,k,8),
cdiag.     tracer(itest,jtest,k,9)
cdiag  enddo
cdiag  call obio_limits('aftr kpp')
#endif
      if (AM_I_ROOT()) then

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (diag_ape)
     . write (501,103) time,'  APE change due to thermal forcing:',
     .  hyc_pechg2(dp(1,1,k1n),th3d(1,1,k1n),31)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
ccc      write (string,'(a12,i8)') 'mxlayr, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
      end if ! AM_I_ROOT
      before = after

      call hybgen(m,n,mm,nn,k1m,k1n)

      call system_clock(after)
      hybgen_time = real(after-before)/real(rate)
c
ccc      write (string,'(a12,i8)') 'hybgrd, step',nstep
ccc      call comparall(m,n,mm,nn,string)

#ifdef TRACERS_OceanBiology
cdiag  do k=1,kdm
cdiag  km=k+mm
cdiag   if (k.eq.1)
cdiag.   write(lp,'(a,3i5,a,a)')'aftr hbge ',nstep,itest,jtest
cdiag.   ,'    dp       nitr      ammo      sili      iron  '
cdiag.   ,'    diat      chlo      cyan      cocco     herb'
cdiag   write(lp,'(22x,i3,10(1x,es9.2))')
cdiag.     k,dp(itest,jtest,km)/onem,
cdiag.     tracer(itest,jtest,k,1),tracer(itest,jtest,k,2),
cdiag.     tracer(itest,jtest,k,3),tracer(itest,jtest,k,4),
cdiag.     tracer(itest,jtest,k,5),tracer(itest,jtest,k,6),
cdiag.     tracer(itest,jtest,k,7),tracer(itest,jtest,k,8),
cdiag.     tracer(itest,jtest,k,9)
cdiag  enddo
cdiag  call obio_limits('aftr hybgen')
#endif

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (AM_I_ROOT()) then
      if (diag_ape)
     . write (501,103) time,'  APE change due to vert.regridding:',
     .  hyc_pechg2(dp(1,1,k1n),th3d(1,1,k1n),31)
c     call sstbud(6,' vert.regridding',temp(1,1,k1n))
c
      if (dotrcr)
     .  write (lp,'(a)') 'tracer advection/turb.mixing cycle completed'
      end if ! AM_I_ROOT
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
        CALL HALO_UPDATE(ogrid, vflxav_loc,  FROM=NORTH)
        i=itest
        j=jtest
        if (haveLatitude(ogrid, J=j)) then
          jb = PERIODIC_INDEX(j+1, jj)
          write (lp,'(i9,2i5,a/a)') nstep,i,j,
     . '  time-integrated continuity eqn diagnostics:',
     .  '     thknss_tndcy  horiz.flxdiv   vert.flxdiv      residuum'
          do k=1,kk
            km=k+mm
            kn=k+nn
            thkchg=dp_loc(i,j,km)+dp_loc(i,j,kn)-dpini(k)
            flxdiv=(uflxav_loc(i+1,j,k)-uflxav_loc(i,j,k)
     .      +vflxav_loc(i,jb ,k)-vflxav_loc(i,j,k))*scp2i_loc(i,j)*delt1
            write (lp,102) k,thkchg,flxdiv,-diaflx_loc(i,j,k),
     .      thkchg+flxdiv-diaflx_loc(i,j,k)
 102  format (i3,4f14.1)
          end do
        end if
      end if
c
!------------------------------------------------------------
      if (diagno .or. mod(time+.0001,1.).lt..0002) then    !  once a day
c
c --- make line printer plot of mean layer thickness
!----note: this section is not used for any thing, so we are
!----      skipping this computation.
        doThis = .false.

        if( doThis ) then
          index=-1
          nflip=mod(nflip+1,2)
          do 705 k=1+(kk-1)*nflip,kk-(kk-1)*nflip,1-2*nflip
ccc       if (k.eq.kk-(kk-1)*nflip) index=1
          sum=0.
          totl=0.
          do 706 j=1,jj
          do 706 l=1,isp(j)
          do 706 i=ifp(j,l),ilp(j,l)
          totl=totl+oice(i,j)*scp2(i,j)
 706      sum=sum+dp(i,j,k+nn)*scp2(i,j)
css       call linout(sum/(area*onem),charac(k),index)
 705      index=0
c ---     add ice extend (%) to plot
css       call linout(100.*totl/area,'I',1)
        end if ! dothis
c
c --- diagnose mean sea surface height
        sumj(:)=0.; sumicej(:)=0.;
        do 704 j=J_0,J_1
        do 704 l=1,isp_loc(j)
        do 704 i=ifp_loc(j,l),ilp_loc(j,l)
c --- compute sea surface height (m)
        srfhgt_loc(i,j)=(montg_loc(i,j,1)+thref*pbavg_loc(i,j,m))/g
        sumicej(j)=sumicej(j)+oice_loc(i,j)*scp2_loc(i,j)
 704    sumj(j)=sumj(j)+srfhgt_loc(i,j)*scp2_loc(i,j)
         call GLOBALSUM(ogrid,sumicej,sumice, all=.true.)
         call GLOBALSUM(ogrid,sumj,sum, all=.true.)

        if (AM_I_ROOT())
     .  write (lp,'(i9,'' mean sea srf.hgt. (mm):'',f9.2,f12.0)')nstep,
     .  sum*1.e3/area,sumice*1.e-6
c
c --- find the largest distance from a tencm layer from bottom
        do 707 j=J_0,J_1
        do 707 l=1,isp_loc(j)
        do 707 i=ifp_loc(j,l),ilp_loc(j,l)
        do 709 k=1,kk
        if (dp_loc(i,j,k+nn).lt.tencm) then
          util2_loc(i,j)=(p_loc(i,j,kk+1)-p_loc(i,j,k))/onem
          goto 708
        endif
 709    continue
 708    continue
 707    continue
c
c       call findmx(ip,util2,idm,ii,jj,'lowest')
c       call findmx(ip,osst,ii,ii,jj,'osst')
c       call findmx(ip,osss,ii,ii,jj,'osss')
c       call findmx(iu,usf,ii,ii,jj,'u_ocn')
c       call findmx(iv,vsf,ii,ii,jj,'v_ocn')
c
c       call findmx(ipa,aflxa2o,iia,iia,jja,'aflxa2o')
c       call findmx(ipa,asst,iia,iia,jja,'asst')
c       call findmx(ipa,sss,iia,iia,jja,'osss')
c       call findmx(ipa,uosurf,iia,iia,jja,'uosurf')
c       call findmx(ipa,vosurf,iia,iia,jja,'vosurf')
c
      end if                      ! once a day
!------------------------------------------------------------
c
c     write (*,'(2i5,2f8.1,a,f9.3,i4/(5(5f9.2,3x,5f9.2/)))')
c    . itest,jtest,latij(itest,jtest,3),lonij(itest,jtest,3),
c    .  ' output tx,ty,ht,ep,t,s,ice,dp,day=',time,nstep
c    . ,((taux(i,j)*100.,j=jtest-2,jtest+2)
c    . , (tauy(i,j)*100.,j=jtest-2,jtest+2),i=itest-2,itest+2)
c    . ,((oflxa2o(i,j),j=jtest-2,jtest+2)
c    . , (oemnp(i,j)*1.e6,j=jtest-2,jtest+2),i=itest-2,itest+2)
c    . ,((temp(i,j,k1n),j=jtest-2,jtest+2)
c    . , (saln(i,j,k1n),j=jtest-2,jtest+2),i=itest-2,itest+2)
c    . ,((oice(i,j),j=jtest-2,jtest+2)
c    . , (dpmixl(i,j)/onem,j=jtest-2,jtest+2),i=itest-2,itest+2)
c
c     write(401,'(i4,f7.3,2f6.1,f6.2,3f8.1)')
c    .  nstep,time,latij(itest,jtest,3),lonij(itest,jtest,3)
c    . ,oice(i,j),dpmixl(i,j)/onem,oflxa2o(i,j),oemnp(i,j)*1.e6
c     do k=1,kk
c     kn=k+nn
c     write(401,'(i3,4f8.2)')
c    . k,p(i,j,k+1)/onem,temp(i,j,kn),saln(i,j,kn),th3d(i,j,kn)
c     enddo
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

      if (AM_I_ROOT()) then
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
      end if  ! AM_I_ROOT

#ifdef TRACERS_OceanBiology
      if (dobio .or. diagno) call obio_trint
#endif

c
      if (.not.diagno) go to 23
c
      if (AM_I_ROOT())
     .  write (lp,100) nstep,int((time+.001)/365.),mod(time+.001,365.)
 100  format (' ocn time step',i9,4x,'y e a r',i6,4x,'d a y',f9.2)
c
c --- output to history file
c
      call gather_before_archive()

      if (AM_I_ROOT()) then
        call archiv(n,nn)

c --- diagnose meridional overturning and heat flux
c$OMP PARALLEL DO
        do 3 j=1,jj
        do 3 k=1,kk
        do 3 l=1,isp(j)
        do 3 i=ifp(j,l),ilp(j,l)
 3      p(i,j,k+1)=p(i,j,k)+dp(i,j,k+mm)
c$OMP END PARALLEL DO
c     call overtn(mm)
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

      endif  !  AM_I_ROOT
      call set_data_after_archiv()

 23   continue

c --- accumulate fields for agcm
c$OMP PARALLEL DO

      do 201 j=J_0,J_1
      do 201 l=1,isp_loc(j)
      do 201 i=ifp_loc(j,l),ilp_loc(j,l)
      osst_loc(i,j)=osst_loc(i,j)+temp_loc(i,j,k1n)*baclin/(3600.*
     &                        real(nhr))
      osss_loc(i,j)=osss_loc(i,j)+saln_loc(i,j,k1n)*baclin/(3600.*
     &                        real(nhr))
      osiav_loc(i,j)=osiav_loc(i,j)+odmsi_loc(i,j)*baclin*dtsrc/(3600.*
     &                        real(nhr)) !kg/m2=>kg*.5*hr/m2
      omlhc_loc(i,j)=spcifh*max(dp_loc(i,j,k1n)/onem,thkmin)/thref  ! J/(m2*C)
      oogeoza_loc(i,j)=(montg_loc(i,j,1)+thref*pbavg_loc(i,j,m))*
     &                        g/(thref*onem) ! m^2/s^2

#ifdef TRACERS_GASEXCH_Natassa
      !here we define the tracer that participates in the
      !gas exchange flux.
#ifdef TRACERS_GASEXCH_CFC_Natassa
      do nt=1,ntm
      otrac(i,j,nt)=otrac(i,j,nt)
     .            +tracer(i,j,1,nt)
     .            *baclin/(3600.*real(nhr))
      enddo
#endif
#ifdef TRACERS_GASEXCH_CO2_Natassa
      !in this case, the tracer is not really active ocean tracer, because it is being
      !handled by the bgcm and it is only at the surface
      do nt=1,ntm
      otrac(i,j,nt)=otrac(i,j,nt)
     .             +pCO2(i,j)                 !pCO2 is in ppmv(uatm)
     .             *baclin/(3600.*real(nhr))
      enddo
#endif
#endif

 201  continue
c$OMP END PARALLEL DO
c
      nsaveo=nsaveo+1

      delt1=baclin+baclin

 15   continue

      if (nsaveo*baclin.ne.nhr*3600) then
            print *, ' ogcm saved over hr=',nsaveo*baclin/3600.
            stop ' stop: ogcm saved over hr'
      end if
      nsaveo=0

      call gather6hycom(osst,osss,osiav,oogeoza,
     &         osst_loc,osss_loc,osiav_loc,oogeoza_loc)
      call gather_hycom_arrays    !mkb 4
      call gather_kprf_arrays

      if (AM_I_ROOT()) then ! work on global grids here
c
c$OMP PARALLEL DO
      do 88 j=1,jj
      do 88 i=1,ii
      usf(i,j)=u(i,j,k1n)+ubavg(i,j,n)
 88   vsf(i,j)=v(i,j,k1n)+vbavg(i,j,n)
c$OMP END PARALLEL DO
c
c     call findmx(ip,osst,ii,ii,jj,'osst')
c     call findmx(ip,osss,ii,ii,jj,'osss')
c     call findmx(iu,usf,ii,ii,jj,'u_ocn')
c     call findmx(iv,vsf,ii,ii,jj,'v_ocn')
c
      call ssto2a(osst,asst)
      call tempro2a(osst,atempr)
      call ssto2a(osss,sss)
css   call iceo2a(omlhc,mlhc)
      call ssto2a(oogeoza,ogeoza)
      call ssto2a(osiav,utila)                 !kg/m*m per agcm time step
      call veco2a(usf,vsf,uosurf,vosurf)
#ifdef TRACERS_GASEXCH_Natassa
      do nt=1,ntm
      call ssto2a(otrac(:,:,nt),atrac(:,:,nt))
      enddo
#endif
#ifdef CHL_from_OBIO
      call ssto2a(tot_chlo,chl)
#endif
c
c     call findmx(ipa,asst,iia,iia,jja,'asst')
c     call findmx(ipa,sss,iia,iia,jja,'osss')
c     call findmx(ipa,uosurf,iia,iia,jja,'uosurf')
c     call findmx(ipa,vosurf,iia,iia,jja,'vosurf')
c
cdiag if (time.le.1)
cdiag.write (*,'(2i5,a,f9.3,i4/(5(5f9.2,3x,5f9.2/)))')
cdiag. iatest,jatest,' output raw uo,vo,ice,focean day=',time,nstep
cdiag. ,((uosurf(i,j)*100.,i=iatest-2,iatest+2)
cdiag. , (vosurf(i,j)*100.,i=iatest-2,iatest+2),j=jatest+2,jatest-2,-1)
cdiag. ,((aice(i,j)*100.,  i=iatest-2,iatest+2)
cdiag. , (focean(i,j)*100.,i=iatest-2,iatest+2),j=jatest+2,jatest-2,-1)
c
c$OMP PARALLEL DO PRIVATE(tf)
      do 204 ia=1,iia
      do 204 ja=1,jja
      if (focean(ia,ja).gt.0.) then
        gtemp(1,1,ia,ja)=asst(ia,ja)
        gtempr(1,ia,ja)=atempr(ia,ja)
#ifdef TRACERS_GASEXCH_Natassa
        do nt=1,ntm
        GTRACER(nt,1,ia,ja)=atrac(ia,ja,nt)
!       if (ia.eq.iatest.and.ja.eq.jatest)
!       if (nstep.eq.25)
!    .  write(6,'(a,3i5,e12.4)')'hycom, atrac at nstep,i,j=',
!    .      nstep,ia,ja,atrac(ia,ja,nt)
        enddo
#endif
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
#ifdef TRACERS_GASEXCH_Natassa
        do nt=1,ntm
        GTRACER(nt,1,ia,ja)=atrac(ia,ja,nt)
        if (ia.eq.iatest.and.ja.eq.jatest)
!    .  write(6,'(a,3i5,e12.4)')
     .  write(6,*)
     .      'hycom, atrac at nstep,i,j=',
     .      nstep,iatest,jatest,atrac(iatest,jatest,nt)
        enddo
#endif
      endif
 204  continue
c$OMP END PARALLEL DO
c --- enhance flow through Denmark Strait
      do j=39,43
      do i=j-10,j-8
      den_str=.08/(2.**(abs(j-41)+abs(i-(j-9))))
      uosurf(i,j)=uosurf(i,j)-den_str
      vosurf(i,j)=vosurf(i,j)-den_str
cdiag if (time.le.1) write(*,'(a,2i4,f8.2)') 'enhanced i,j=',i,j,den_str
      enddo
      enddo
c
cdiag if (time.le.1)
cdiag.write (*,'(2i5,a,f9.3,i4/(5(5f9.2,3x,5f9.2/)))')
cdiag. iatest,jatest,' output enh. uo,vo,ice,focean day=',time,nstep
cdiag. ,((uosurf(i,j)*100.,i=iatest-2,iatest+2)
cdiag. , (vosurf(i,j)*100.,i=iatest-2,iatest+2),j=jatest+2,jatest-2,-1)
cdiag. ,((aice(i,j)*100.,  i=iatest-2,iatest+2)
cdiag. , (focean(i,j)*100.,i=iatest-2,iatest+2),j=jatest+2,jatest-2,-1)
c     call prtmsk(ipa,uosurf,util3(31,31),iia,10,15,0.,1000.,'uo(mm/s)')
c     call prtmsk(ipa,vosurf,util3(31,31),iia,10,15,0.,1000.,'vo(mm/s)')
c
      call system_clock(afogcm)
      if (abs(time-(itime+1.)/nday).gt..01) then
        write(*,'(a,2f10.3)') 'stop: agcm/ogcm date after ogcm'
     .    ,(itime+1.)/nday,time
        stop ' stop: dates in agcm/ogcm do not match'
      end if
c

      endif ! AM_I_ROOT
 9666 continue
      call scatter_atm

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
c------------------------------------------------------------------
      subroutine gather2_atm

      USE HYCOM_ATM !, only : gather_atm, scatter_atm
      USE DOMAIN_DECOMP, ONLY: GRID, PACK_DATA

      call pack_data( grid,  ataux_loc,   ataux )
      call pack_data( grid,  atauy_loc,   atauy )
      call pack_data( grid,  aflxa2o_loc, aflxa2o )
      call pack_data( grid,  aemnp_loc,   aemnp )
      call pack_data( grid,  aice_loc,    aice )
      call pack_data( grid,  asalt_loc,   asalt )
      call pack_data( grid,  austar_loc,  austar )
      call pack_data( grid,  aswflx_loc,  aswflx )

      end subroutine gather2_atm
c------------------------------------------------------------------
      subroutine gather_before_archive
      USE HYCOM_ARRAYS_GLOB
      use hycom_arrays_glob_renamer
      USE HYCOM_DIM, only : ogrid
      USE DOMAIN_DECOMP, ONLY: PACK_DATA

      call pack_data( ogrid,  u_loc, u )
      call pack_data( ogrid,  v_loc, v )
      call pack_data( ogrid,  dp_loc, dp )
      call pack_data( ogrid,  p_loc, p )
      call pack_data( ogrid,  temp_loc, temp )
      call pack_data( ogrid,  saln_loc, saln )
      call pack_data( ogrid,  th3d_loc, th3d )
      call pack_data( ogrid,  dpmixl_loc, dpmixl )
      call pack_data( ogrid,  srfhgt_loc, srfhgt )
      call pack_data( ogrid,  ubavg_loc, ubavg )
      call pack_data( ogrid,  vbavg_loc, vbavg )
      call pack_data( ogrid,  uav_loc, uav )
      call pack_data( ogrid,  vav_loc, vav )
      call pack_data( ogrid,  dpuav_loc, dpuav )
      call pack_data( ogrid,  dpvav_loc, dpvav )
      call pack_data( ogrid,  temav_loc, temav )
      call pack_data( ogrid,  salav_loc, salav )
      call pack_data( ogrid,  th3av_loc, th3av )
      call pack_data( ogrid,  dpav_loc, dpav )
      call pack_data( ogrid,  ubavav_loc, ubavav )
      call pack_data( ogrid,  vbavav_loc, vbavav )
      call pack_data( ogrid,  pbavav_loc, pbavav )
      call pack_data( ogrid,  sfhtav_loc, sfhtav )
      call pack_data( ogrid,  uflxav_loc, uflxav )
      call pack_data( ogrid,  vflxav_loc, vflxav )
      call pack_data( ogrid,  diaflx_loc, diaflx )
      call pack_data( ogrid,  sflxav_loc, sflxav )
      call pack_data( ogrid,  brineav_loc, brineav )
      call pack_data( ogrid,  eminpav_loc, eminpav )
      call pack_data( ogrid,  surflav_loc, surflav )
      call pack_data( ogrid,  dpmxav_loc, dpmxav )
      call pack_data( ogrid,  oiceav_loc, oiceav )
      call pack_data( ogrid,  pbot_loc, pbot )
      call pack_data( ogrid,  tracer_loc, tracer )
      call pack_data( ogrid,  oice_loc, oice )
      call pack_data( ogrid,  util1_loc, util1 )

      end subroutine gather_before_archive
c------------------------------------------------------------------
      subroutine set_data_after_archiv

      USE HYCOM_DIM, only : ii, J_0,  J_1, kk, ogrid 
      USE HYCOM_ARRAYS_GLOB, only : p_glob => p, util1_glob => util1,
     &                                           util2_glob => util2
      USE HYCOM_ARRAYS
      USE DOMAIN_DECOMP, ONLY: UNPACK_DATA

      implicit none
      integer :: i, j, k

      call unpack_data( ogrid,  p_glob, p )
      call unpack_data( ogrid,  util1_glob, util1 )
      call unpack_data( ogrid,  util2_glob, util2 )

      do 60 j=J_0, J_1
c
      do 601 i=1,ii
      eminpav(i,j)=0.
      surflav(i,j)=0.
       sflxav(i,j)=0.
      brineav(i,j)=0.
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

      end subroutine set_data_after_archiv
c------------------------------------------------------------------
      subroutine gather6hycom(osst,osss,osiav,oogeoza,
     &         osst_loc,osss_loc,osiav_loc,oogeoza_loc)

      USE HYCOM_DIM, only : idm, jdm, J_0H,  J_1H, ogrid
      USE HYCOM_ARRAYS_GLOB, only : omlhc
      USE hycom_arrays_glob_renamer, only : omlhc_loc
      USE DOMAIN_DECOMP, ONLY: GRID, PACK_DATA
      implicit none
      real osst(idm,jdm),osss(idm,jdm),osiav(idm,jdm),oogeoza(idm,jdm)
      real osst_loc(idm,J_0H:J_1H),osss_loc(idm,J_0H:J_1H),
     &     osiav_loc(idm,J_0H:J_1H),oogeoza_loc(idm,J_0H:J_1H)
     
      call pack_data( ogrid,  osst_loc,      osst    )
      call pack_data( ogrid,  osss_loc,      osss    )
      call pack_data( ogrid,  osiav_loc,     osiav   )
      call pack_data( ogrid,  oogeoza_loc,   oogeoza )
      call pack_data( ogrid,  omlhc_loc,     omlhc   )

      end subroutine gather6hycom
