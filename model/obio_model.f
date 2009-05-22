#include "rundeck_opts.h"

#ifdef OBIO_ON_GARYocean
      subroutine obio_model
#else
      subroutine obio_model(nn,mm)
#endif

!@sum  OBIO_MODEL is the main ocean bio-geo-chem routine 
!@auth Natassa Romanou/Watson Gregg
!@ver  1.0

      USE obio_dim
      USE obio_incom
      USE obio_forc, only: solz,tirrq,Ed,Es
     .                    ,rmud,atmFe,avgq,ihra,sunz
     .                    ,wind
     .                    ,atmFe_all
     .                    ,owind,osolz
     .                    ,alk
     .                    ,tirrq3d
#ifdef OBIO_RAD_coupling
     .                    ,eda_frac,esa_frac
     .                    ,ovisdir,ovisdif,onirdir,onirdif
     .                    ,ovisdir_ij,ovisdif_ij,onirdir_ij,onirdif_ij
#else
     .                    ,Eda,Esa,Eda2,Esa2
#endif
      USE obio_com,  only: gcmax,day_of_month,hour_of_day
     .                    ,temp1d,dp1d,obio_P,det,car,avgq1d
     .                    ,ihra_ij,gcmax1d,atmFe_ij,covice_ij
     .                    ,P_tend,D_tend,C_tend,saln1d
     .                    ,pCO2,pCO2_ij,p1d,wsdet
     .                    ,rhs,alk1d
     .                    ,tzoo,tfac,rmuplsr,rikd,bn,wshc,Fescav
     .                    ,tzoo2d,tfac3d,rmuplsr3d,rikd3d
     .                    ,bn3d,obio_wsd2d,obio_wsh2d,wshc3d,Fescav3d 
     .                    ,acdom,pp2_1d,pp2tot_day
     .                    ,tot_chlo,acdom3d,pnoice
#ifdef OBIO_ON_GARYocean
     .                    ,itest,jtest,obio_deltat,nstep0
     .                    ,tracer =>tracer_loc        
     .                    ,tracer_glob => tracer
      USE ODIAG, only : ij_pCO2,ij_dic,ij_nitr,ij_diat
     .                 ,ij_amm,ij_sil,ij_chlo,ij_cyan,ij_cocc,ij_herb
     .                 ,ij_doc,ij_iron,ij_alk
      USE ODIAG, only : oij=>oij_loc
#endif

      USE MODEL_COM, only: JMON,jhour,nday,jdate,jday
     . ,itime,iyear1,jdendofm,jyear,aMON
     . ,xlabel,lrunid

      USE FILEMANAGER, only: openunit,closeunit

#ifdef TRACERS_GASEXCH_ocean_CO2
      USE TRACER_COM, only : ntm    !tracers involved in air-sea gas exch
      USE TRACER_GASEXCH_COM, only : tracflx,tracflx1d
#endif


#ifdef OBIO_ON_GARYocean
      USE MODEL_COM,  only : nstep=>itime,itimei
      USE GEOM,       only : LON_DG,LAT_DG
      USE CONSTANT,   only : grav
      USE OCEANR_DIM, only : ogrid
      USE OCEANRES,   only : kdm=>lmo,dzo
      USE OFLUXES,    only : oRSI,oAPRESS
      USE OCEAN,      only : ZOE=>ZE,g0m,s0m,mo,dxypo,focean,lmm
     .                      ,trmo,txmo,tymo,tzmo
      USE KPP_COM,    only : kpl,kpl_glob
      USE OCN_TRACER_COM,    only : obio_tr_mm
#else
      USE hycom_dim
      USE hycom_arrays, only: tracer,dpinit,temp,saln,oice
     .                            ,p,dpmixl,latij,lonij
      USE hycom_scalars, only: trcout,nstep,onem,nstep0
     .                        ,time,lp,itest,jtest
#endif

      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT,pack_data,unpack_data


      implicit none

      integer i,j,k,l,km,nn,mm 

      integer ihr,ichan,iyear,nt,ihr0,lgth,kmax
      integer iu_pco2,ll,iu_tend
      real    tot,dummy(6),dummy1
      real    rod(nlt),ros(nlt)
#ifdef OBIO_ON_GARYocean
      Real*8,External   :: VOLGSP
      real*8 temgs,g,s,temgsp,pres
      real*8 time,dtr,ftr,rho_water
      real*8 trmo_unit_factor(kdm,ntrac)
      integer i_0,i_1,j_0,j_1
#endif
      character string*80
      character jstring*3

      logical vrbos,noon,errcon

!--------------------------------------------------------
      diagno_bio=.false.
#ifdef OBIO_ON_GARYocean
      if (JDendOfM(jmon).eq.jday.and.Jhour.eq.12) then
          if (mod(nstep,2).eq.0)
     .    diagno_bio=.true. ! end of month,mid-day
      endif
#else
      if (JDendOfM(jmon).eq.jday.and.Jhour.eq.12)
     .    diagno_bio=.true. ! end of month,mid-day
#endif

!Cold initialization

#ifdef OBIO_ON_GARYocean
      print*, 'itimei,nstep,nstep0 =',
     .         itimei,nstep,nstep0

      if (nstep.eq.itimei) then
      print*, 'COLD INITIALIZATION....'
      nstep0=0
#else
      if (nstep.eq.1) then
        trcout=.true.
#endif

      if (AM_I_ROOT()) write(*,'(a)')'BIO:Ocean Biology starts ....'

      call obio_init
      !tracer array initialization.
      !note: we do not initialize obio_P,det and car
#ifdef OBIO_ON_GARYocean
      call obio_bioinit_g
#else
      call obio_bioinit(nn)
#endif
      endif   !if nstep=1 or nstep=itimei



!Warm initialization

#ifdef OBIO_ON_GARYocean
      if (nstep0 .gt. itimei .and. nstep.eq.nstep0) then
         call obio_init

         print*,'WARM INITIALIZATION'
         call obio_trint

      endif !for restart only
#else
       if (nstep0 .gt. 0 .and. nstep.eq.nstep0+1) then
         write(*,'(a)')'For restart runs.....'
         write(*,'(a,2i9,f10.3)')
     .            'nstep0,nstep,time=',nstep0,nstep,time
         call obio_init

         print*,'WARM INITIALIZATION'
         call obio_trint
       endif !for restart only
#endif

#ifdef OBIO_ON_GARYocean
      time = float(nstep)
      i_0=ogrid%I_STRT
      i_1=ogrid%I_STOP
      j_0=ogrid%J_STRT
      j_1=ogrid%J_STOP
#endif

#ifdef OBIO_ON_GARYocean
      write(*,'(a,2i5,2e12.4)')'TEST POINT at: ',itest,jtest,
     .      LAT_DG(jtest,2),LON_DG(itest,2)
#else
      write(*,'(a,2i5,2e12.4)')'TEST POINT at: ',itest,jtest,
     .      latij(itest,jtest,3),lonij(itest,jtest,3)
#endif

!--------------------------------------------------------

       day_of_month=jdate
       hour_of_day =jhour

       !if (time.kmod(nstep,24*30).eq.0) then
       !   day_of_month=1
       !else
       !   day_of_month=day_of_month+1
       !endif

       !if (mod(nstep,24).eq.0) then
       !  hour_of_day=1
       !else
       !  hour_of_day=hour_of_day+1
       !endif
 
      if (AM_I_ROOT()) then
       write(*,'(a,i15,1x,f9.3,2x,3i5)')
     .    'BIO: nstep,time,day_of_month,hour_of_day,jday=',
     .    nstep,time,day_of_month,hour_of_day,jday
      endif

         ihr0 = int(hour_of_day/2)

      if (diagno_bio) then
#ifdef OBIO_ON_GARYocean
        write(string,'(a3,i4.4,2a)') amon,Jyear,'.',xlabel(1:lrunid)
#else
!       do 16 lgth=60,1,-1
!         if (flnmovt(lgth:lgth).ne.' ') go to 17
! 16    continue
!       stop '(bioinit:flnmovt)'
! 17    write (string,'(i6.6)') int(time+.001)
!!      open(16,file=flnmovt(1:lgth)//'pco2_out.'//string(1:6)
!!    .        ,status='unknown')
        if (nstep.eq.1) then
         write(string,'(a3,i4.4,2a)') amon,0,'.',xlabel(1:lrunid)
       elseif (abs((itime+1.)/nday-time).gt.1.e-5) then
         write(*,*) 'mismatching archive date in agcm/ogcm=',
     .      (itime+1.)/nday,time
         stop 'mismatching archive date'
       else
         write(string,'(a3,i4.4,2a)') amon,Jyear,'.',xlabel(1:lrunid)
       endif
#endif

!       print*,' '
!       print*, 'BIO: saving in pco2 and tend files'
!       print*,' '
!       jstring='xxx'
!#ifdef OBIO_ON_GARYocean
!       if(j_1.lt.100)write(jstring(1:2),'(i2)')j_1
!       if(j_1.ge.100)write(jstring(1:3),'(i3)')j_1
!#else
!       if(j_1h.lt.100)write(jstring(1:2),'(i2)')j_1h
!       if(j_1h.ge.100)write(jstring(1:3),'(i3)')j_1h
!#endif
!      print*, 'pco2.'//jstring//string
!      call openunit('pco2.'//jstring//string,iu_pco2)
!      call openunit('tend.'//jstring//string,iu_tend)
!
      endif  !diagno_bio

#ifdef OBIO_ON_GARYocean
      write(*,'(/,a,2i5,2e12.4)')'obio_model, test point=',
     .      itest,jtest,LON_DG(itest,1),LAT_DG(jtest,1)
#endif


c$OMP PARALLEL DO PRIVATE(km,iyear,kmax,vrbos,errcon,tot,noon,rod,ros)
c$OMP. SHARED(hour_of_day,day_of_month,JMON)

#ifdef OBIO_ON_GARYocean
       do 1000 j=j_0,j_1
       do 1000 i=i_0,i_1
       IF(FOCEAN(I,J).gt.0.) THEN
#else
       do 1000 j=j_0h,j_1h             !1,jj
       do 1000 l=1,isp(j)
       do 1000 i=ifp(j,l),ilp(j,l)
#endif

cdiag  write(*,'(/,a,i5,2i4)')'obio_model, step,i,j=',nstep,i,j

       vrbos=.false.
       if (i.eq.itest.and.j.eq.jtest) vrbos=.true.

       !surface forcing atmFe test case
       !atmFe(i,j)=0.

       !fill in reduced rank arrays
       ihra_ij=ihra(i,j)
       atmFe_ij=atmFe(i,j)
#ifdef OBIO_ON_GARYocean
       covice_ij=oRSI(i,j)
#else
       !!covice_ij=covice(i,j)  !for standalone hycom
       covice_ij=oice(i,j)      !for modelE-hycom
#endif

       pCO2_ij=pCO2(i,j)

#ifdef OBIO_ON_GARYocean
       pres = oAPRESS(i,j)    !surface atm. pressure
       do k=1,lmm(i,j)
         s=S0M(I,J,k)/(MO(I,J,k)*DXYPO(J))
!!!!     temp1d(k)=TEMGS(g,s)           !potential temperature
         temp1d(k)=TEMGSP(g,s,pres)     !in situ   temperature
          saln1d(k)=s*1000.             !convert to psu (eg. ocean mean salinity=35psu)
           dp1d(k)=dzo(k)               !thickenss of each layer in meters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! array tracer has units [mole]/m3. convert to/from trmo with units kg as follows:
! trmo = tracer * MB*1e-3/rho_water     * mo * dxypo
! note that occassionally tracer is in milimoles and othertimes in micromoles
! trmo_unit_factor below has units m^3 and depends on which 
! layer you are: trmo = tracer * factor,  tracer=trmo/factor
! txmo,tzmo,tymo are all computed based on trmo
! this should be done at every timestep except when COLD INITIALIZATION
! because then trmo=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         rho_water = 1d0/VOLGSP(g,s,pres)

         if(vrbos.and.k.eq.1)write(*,'(a,4e12.4)')
     .             'obio_model,t,s,p,rho= '
     .             ,temp1d(k),saln1d(k),dp1d(k),rho_water

         do nt=1,ntrac
           trmo_unit_factor(k,nt) = 1d-3*1d-3*obio_tr_mm(nt)        ! milimoles/m3=> kg/m3
     .                            *  MO(I,J,k)*DXYPO(J)/rho_water   ! kg/m3 => kg
           if (nt.eq.4.or.nt.eq.13)
     .         trmo_unit_factor(k,nt) = 1d-3*trmo_unit_factor(k,nt) !nanomoles/m3 => kg
           if (nt.ge.5.and.nt.le.9)
     .         trmo_unit_factor(k,nt) =  1d-3*1d-3 
     .                                *  MO(I,J,k)*DXYPO(J)/rho_water ! mg/m3 => kg
           if (nt.eq.ntrac)    !factor for alkalinity
     .         trmo_unit_factor(k,nt) = 1d-3*1d-3*1d-3*obio_tr_mm(nt) ! umol/kg=micro-mol/kg=> kg,trac/kg,air
     .                                *  MO(I,J,k)*DXYPO(J)           ! kg,trac/kg,air=> kg,trac

           if (nstep0 .gt. itimei) then
              tracer(i,j,k,nt) = trmo(i,j,k,nt) / trmo_unit_factor(k,nt)
           endif
         enddo
#else
       do k=1,kdm
        km=k+mm
         temp1d(k)=temp(i,j,km)
          saln1d(k)=saln(i,j,km)
           dp1d(k)=dpinit(i,j,k)/onem
#endif
            avgq1d(k)=avgq(i,j,k)
             gcmax1d(k)=gcmax(i,j,k)
              tirrq(k)=tirrq3d(i,j,k)
#ifndef TRACERS_Alkalinity
              alk1d(k)=alk(i,j,k)
#endif
              !----daysetbio/daysetrad arrays----!
              tzoo=tzoo2d(i,j)
              tfac(k)=tfac3d(i,j,k)
              do nt=1,nchl
                rmuplsr(k,nt)=rmuplsr3d(i,j,k,nt)
                rikd(k,nt)=rikd3d(i,j,k,nt)
                obio_wsd(nt)=obio_wsd2d(i,j,nt)
                obio_wsh(nt)=obio_wsh2d(i,j,nt)
              enddo
              bn(k)=bn3d(i,j,k)
              wshc(k)=wshc3d(i,j,k)
              Fescav(k)=Fescav3d(i,j,k)
              do nt=1,nlt
                acdom(k,nt)=acdom3d(i,j,k,nt)
              enddo
              !----daysetbio arrays----!
              do nt=1,ntyp+n_inert
             obio_P(k,nt)=tracer(i,j,k,nt)
            enddo
           do nt=1,ndet
          det(k,nt)=tracer(i,j,k,ntyp+n_inert+nt)
         enddo
        do nt=1,ncar
       car(k,nt)=tracer(i,j,k,ntyp+n_inert+ndet+nt)
       enddo
#ifdef TRACERS_Alkalinity
       do nt=1,nalk
       alk1d(k)=tracer(i,j,k,ntyp+n_inert+ndet+ncar+nt)
       enddo
#endif
       enddo  !k=1,kdm or lmm

       p1d(1)=0.
       do k=2,kdm+1
          p1d(k)=p1d(k-1)+dp1d(k-1)    !in meters
       enddo

#ifdef OBIO_ON_GARYocean
      kmax = lmm(i,j)
cdiag if(vrbos)write(*,'(a,4i5)')'nstep,i,j,kmax= ',
cdiag.         nstep,i,j,kmax

#else

! --- find number of non-massless layers in column
!     do kmax=kdm,1,-1
!       if (dp1d(kmax).gt.0.) exit
!       if (kmax.eq.1) then
!         write (*,*) 'error - zero column depth'
!         return
!       end if
!     end do
      do kmax=1,kdm
        if (dp1d(kmax).lt.1.e-2) exit
      end do
      kmax=kmax-1
cdiag write(*,'(a,4i5)')'nstep,i,j,kmax= ',nstep,i,j,kmax

      !ensure that all massless layers have same concentration
      !as last mass one
      do k=kmax+1,kdm
       do nt=1,ntyp+n_inert
        if(obio_P(kmax,nt).le.0.)obio_P(kmax,nt)=1.e-8
        obio_P(k,nt)=obio_P(kmax,nt)
       enddo
       do nt=1,ndet
        det(k,nt)=   det(kmax,nt)
       enddo
       do nt=1,ncar
        car(k,nt)=   car(kmax,nt)
       enddo
      enddo
#endif


#ifdef OBIO_RAD_coupling
       ovisdir_ij=ovisdir(i,j)
       ovisdif_ij=ovisdif(i,j)
       onirdir_ij=onirdir(i,j)
       onirdif_ij=onirdif(i,j)

       if (vrbos)write(*,'(/,a,3i5,4e12.4)')
     .    'obio_model, radiation: ',
     .    nstep,i,j,ovisdir_ij,ovisdif_ij,onirdir_ij,onirdif_ij
#else
       !OASIM spectral irradiance data just above the surface
       !Eda and Esa fields have ihr=1:12, ie every 2hrs
!      do ihr=1,nhn
!       do ichan=1,nlt
!        Eda2(ichan,ihr)=Eda(i,j,ichan,ihr,l0)*w0
!    .                  +Eda(i,j,ichan,ihr,l1)*w1
!    .                  +Eda(i,j,ichan,ihr,l2)*w2
!    .                  +Eda(i,j,ichan,ihr,l3)*w3
!        Esa2(ichan,ihr)=Esa(i,j,ichan,ihr,l0)*w0
!    .                  +Esa(i,j,ichan,ihr,l1)*w1
!    .                  +Esa(i,j,ichan,ihr,l2)*w2
!    .                  +Esa(i,j,ichan,ihr,l3)*w3
!       enddo
!      enddo
#endif

       !solz is read inside hycom.f and forfun.f
!      do ihr=1,12
!       solz2(ihr)=solz_all(i,j,ihr,l0)*w0
!    .            +solz_all(i,j,ihr,l1)*w1
!    .            +solz_all(i,j,ihr,l2)*w2
!    .            +solz_all(i,j,ihr,l3)*w3

!       sunz2(ihr)=acos(solz2(ihr))*rad   !in degs
!      enddo

!        solz=solz2(ihr0)     !because of bi-hourly values
!        sunz=sunz2(ihr0)     !in degs
         solz=osolz(i,j)      !read instead from modelE
         sunz=acos(solz)*rad  !in degs


!      wind=wndspd(i,j,l0)*w0+wndspd(i,j,l1)*w1
!    .     +wndspd(i,j,l2)*w2+wndspd(i,j,l3)*w3

       wind=owind(i,j)

!      atmFe_ij=atmFe_all(i,j,l0)*w0 + atmFe_all(i,j,l1)*w1
!    .         +atmFe_all(i,j,l2)*w2 + atmFe_all(i,j,l3)*w3

       !atmospheric deposition iron * solubility
       atmFe_ij=atmFe_all(i,j,JMON)*solFe

       if (vrbos) then
         write(*,'(/,a,3i5,4e12.4)')'obio_model, forcing: ',
     .   nstep,i,j,solz,sunz,wind,atmFe_ij
       endif

#ifdef TRACERS_GASEXCH_ocean_CO2
       do nt=1,ntm
          tracflx1d(nt) = tracflx(i,j,nt)
!         write(*,'(/,a,3i5,2e12.4)')'obio_model, tracflx:',
!    .        nstep,i,j,tracflx(i,j,nt),tracflx1d(nt)
       enddo
#endif

       !------------------------------------------------------------
       !at the beginning of each day only
#ifdef OBIO_ON_GARYocean
       if (hour_of_day.le.1) then    
#else
       if (hour_of_day.eq.1) then
#endif

         if (day_of_month.eq.1)ihra_ij=1
          call obio_daysetrad(vrbos,i,j)
          ihra_ij = 0
          call obio_daysetbio(vrbos,i,j)

             !fill in the 3d arrays to keep in memory for rest of the day
             !when daysetbio is not called again
             tzoo2d(i,j)=tzoo
             do  k=1,kdm
               tfac3d(i,j,k)=tfac(k)
             do nt=1,nchl
               rmuplsr3d(i,j,k,nt)=rmuplsr(k,nt)
               rikd3d(i,j,k,nt)=rikd(k,nt)
               obio_wsd2d(i,j,nt)=obio_wsd(nt)
               obio_wsh2d(i,j,nt)=obio_wsh(nt)
             enddo
               bn3d(i,j,k)=bn(k)
               wshc3d(i,j,k)=wshc(k)
               Fescav3d(i,j,k)=Fescav(k)
             do nt=1,nlt
               acdom3d(i,j,k,nt)=acdom(k,nt)
             enddo
             enddo

         if (day_of_month.eq.1)ihra_ij=1


cdiag    if (vrbos) then
cdiag     write (*,103) nstep,i,j,
cdiag.' aftrsetrad   dpth     dp       nitr    ',
cdiag.'   ammo     sili     iron',
cdiag.    (k,p1d(k+1),dp1d(k),obio_P(k,1),obio_P(k,2),
cdiag.                        obio_P(k,3),obio_P(k,4),k=1,kdm)

cdiag     write (*,104) nstep,i,j,
cdiag.' aftrsetrad   dpth     dp       diat       chlo    ',
cdiag.' cyan       cocc     herb',
cdiag.    (k,p1d(k+1),dp1d(k),obio_P(k,5),obio_P(k,6),obio_P(k,7),
cdiag.                        obio_P(k,8),obio_P(k,9),k=1,kdm)
cdiag    endif
 103     format(i9,2i5,a,a/(25x,i3,6(1x,es9.2)))
 104     format(i9,2i5,a,a/(25x,i3,7(1x,es9.2)))

       endif   !end of calculations for the beginning of day


       !------------------------------------------------------------
#ifndef OBIO_RAD_coupling
       !Eda and Esa OASIM data is given every 2hrs
       if (mod(hour_of_day,2) .eq. 0) then
       !only every 2 hrs
#endif

         do ichan = 1,nlt
           Ed(ichan) = 0.0
           Es(ichan) = 0.0
         enddo
         rmud = 0.0
         iyear=2001


         !compute rod and ros only here. not ocean albedo.
         !ocean albedo is computed in ALBEDO.f
         !have to have hygr =  .true. 
         call obio_ocalbedo(wind,solz,dummy,dummy,dummy1,
     .                      rod,ros,.true.,vrbos,i,j)


cdiag    if (vrbos)
cdiag.     write(*,105)nstep,
cdiag.  ' channel, surf refl dir, surf refl diff',
cdiag.                        (k,rod(k),ros(k),k=1,nlt)
 105  format(i9,a/(18x,i3,2(1x,es9.2)))

         !only call obio_sfcirr for points in light
         tot = 0.0
         do ichan = 1,nlt
#ifdef OBIO_RAD_coupling
          if (ichan .le. 18) then     !visible + uv
              Ed(ichan) = ovisdir_ij * eda_frac(ichan)
              Es(ichan) = ovisdif_ij * esa_frac(ichan)
          else               !nir
              Ed(ichan) = onirdir_ij * eda_frac(ichan)
              Es(ichan) = onirdif_ij * esa_frac(ichan)
          endif
          tot = tot + Ed(ichan)+Es(ichan)

cdiag if (nstep.eq.12)
cdiag.write(*,'(/,a,4i5,3e12.4)')'obio_model, tirrq: ',
cdiag.           nstep,i,j,ichan,
cdiag.           ovisdir_ij,eda_frac(ichan),Ed(ichan)

#else
          Ed(ichan) = Eda2(ichan,ihr0)
          Es(ichan) = Esa2(ichan,ihr0)
          tot = tot + Eda2(ichan,ihr0)+Esa2(ichan,ihr0)

cdiag if (nstep.eq.12)
cdiag.write(*,'(/,a,4i5,3e12.4)')'obio_model, tirrq: ',
cdiag.           nstep,i,j,ichan,Ed(ichan)

#endif
         enddo  !ichan
         noon=.false.
         if (hour_of_day.eq.12)then
          if (i.eq.itest.and.j.eq.jtest)noon=.true.
         endif
         if (tot .ge. 0.1) call obio_sfcirr(noon,rod,ros,vrbos)
      
      !check
      if (vrbos) then
        write(*,'(a,3i9)')
     .       'obio_model: counter days,   i,j,nstep=',i,j,nstep
        write(*,'(3(a,i9))')'hour of day=',hour_of_day,
     .                       ', day of month=',day_of_month,
     .                       ', ihr0=',ihr0

cdiag   write(*,'(a)')
cdiag.'    k     dp          u            v         temp         saln'
cdiag   do k=1,kdm
cdiag   write(*,'(i5,5e12.4)')
cdiag   write(*,*)
cdiag.  k,dp1d(k),u(i,j,k+mm),v(i,j,k+mm),
cdiag.    temp(i,j,k+mm),saln(i,j,k+mm)
cdiag   enddo

        write(*,'(a)')
     .'    k     P(1)      P(2)         P(3)       P(4)         P(5) '
        do k=1,kdm
        write(*,'(i5,7e12.4)')
!       write(*,*)
     .   k,obio_P(k,1),obio_P(k,2),obio_P(k,3),obio_P(k,4),
     .   obio_P(k,5),obio_P(k,6),obio_P(k,7)
        enddo

        write(*,'(a)')
     .'    k     P(8)      P(9)         P(11)      P(12)        P(13)'
        do k=1,kdm
        write(*,'(i5,7e12.4)')
!       write(*,*)
     .   k,obio_P(k,8),obio_P(k,9),det(k,1),det(k,2),
     .   det(k,3),car(k,1),car(k,2)
        enddo

        write(*,'(2a)')
cdiag.'    Ed          Es          solz         sunz',
cdiag.'       atmFe       wind'
!       write(*,'(6e12.4)')
cdiag   write(*,*)
cdiag.   Ed(ichan),Es(ichan),solz,sunz,atmFe_ij,wind
       endif

cdiag    if (vrbos)
cdiag.     write(*,106)nstep,
cdiag.    '      channel, dir dwn irr, diff dwn irr,  tot',
cdiag.                  (ichan,Ed(ichan),Es(ichan),
cdiag.                  tot,ichan=1,nlt)
 106  format(i9,a/(18x,i3,3(1x,es9.2)))


         !this part decomposes light into gmao bands-- dont need it
         !         m = indext2   !use the indicator of "past"
         !         call gmaoirr(ihr,nwater,Edg,Esg)
         !         call oasimbands(idmo,ihr,nwater,Edg,Esg)

         do k=1,kdm
          tirrq(k) = 0.0
         enddo

         if (tot .ge. 0.1) call obio_edeu(kmax,vrbos,i,j)

cdiag    if (vrbos)
cdiag.     write(*,107)nstep,
cdiag.       '           k   avgq    tirrq',
cdiag.                 (k,avgq1d(k),tirrq(k),k=1,kdm)
 107  format(i9,a/(18x,i3,2(1x,es9.2)))


         if (tot .ge. 0.1) ihra_ij = ihra_ij + 1

#ifndef OBIO_RAD_coupling
       endif   !every even hour of the day
#endif

       !------------------------------------------------------------
       !compute tendency terms on the m level
cdiag     if (vrbos) then
cdiag     write (*,103) nstep,i,j,
cdiag.' bfreptend  dpth       dp       nitr    ',
cdiag.'     ammo     sili     iron',
cdiag.    (k,p(i,j,k+1)/onem,dp1d(k),obio_P(k,1),obio_P(k,2),
cdiag.                       obio_P(k,3),obio_P(k,4),k=1,kdm)

cdiag     write (*,104) nstep,i,j,
cdiag.' bfreptend  dpth     diat       chlo    ',
cdiag.' cyan       cocc     herb',
cdiag.    (k,p(i,j,k+1)/onem,obio_P(k,5),obio_P(k,6),obio_P(k,7),
cdiag.                       obio_P(k,8),obio_P(k,9),k=1,kdm)

cdiag     write (*,104) nstep,i,j,
cdiag.' bfreptend  dpth     nitr_det  sili_det iron_det   doc   ',
cdiag.'    dic ',
cdiag.    (k,p(i,j,k+1)/onem,det(k,1),det(k,2),det(k,3),
cdiag.                       car(k,1),car(k,2),k=1,kdm)
cdiag     endif
       !------------------------------------------------------------

cdiag  if (vrbos)write(*,*)'bfre obio_ptend: ',
cdiag.     nstep,(k,tirrq(k),k=1,kmax)

       call obio_ptend(vrbos,kmax,i,j)

       !------------------------------------------------------------
cdiag  if (vrbos)then
cdiag   write(*,108)nstep,' aftrptend dpth      dp     P_tend(1:9)',
cdiag.    (k,p(i,j,k+1)/onem,dp1d(k),P_tend(k,1),P_tend(k,2),
cdiag.     P_tend(k,3),P_tend(k,4),P_tend(k,5),P_tend(k,6),
cdiag.     P_tend(k,7),P_tend(k,8),P_tend(k,9),k=1,kdm)
cdiag   write(*,109)nstep,
cdiag.    ' aftrptend dpth      dp     D_tend(1:3) and C_tend(1:2)',
cdiag.    (k,p(i,j,k+1)/onem,dp1d(k),D_tend(k,1),D_tend(k,2),
cdiag.     D_tend(k,3),C_tend(k,1),C_tend(k,2),k=1,kdm)
cdiag  endif
 108  format(i6,a/(12x,i3,11(1x,es9.2)))
 109  format(i6,a/(12x,i3,7(1x,es9.2)))

       !------------------------------------------------------------

#ifdef OBIO_ON_GARYocean
       !update biology to new time level
       !also do phyto sinking and detrital settling here
       !MUST CALL sinksettl BEFORE update
       call obio_sinksettl(vrbos,kmax,errcon,i,j)
       call obio_update(vrbos,kmax,i,j)
#else
       !update biology from m to n level
       !also do phyto sinking and detrital settling here
       !MUST CALL sinksettl AFTER update
       call obio_update(vrbos,kmax,i,j)
       call obio_sinksettl(vrbos,kmax,errcon,i,j)
       if (errcon) then
          write (*,'(a,2i5)') 'error update at i,j =',i,j
          do k=1,kdm
          write (*,'(i5,3e12.4)')k,dp1d(k),p1d(k),wsdet(k,1)
          enddo
          stop
       endif
#endif

cdiag     if (vrbos) then
cdiag     write (*,*)'     '
cdiag     write (*,103) nstep,i,j,
cdiag.' aftrupdate dpth     dp        nitr      ammo      sili',
cdiag.'      iron',
cdiag.  (k,p(i,j,k+1)/onem,dp1d(k),obio_P(k,1),obio_P(k,2),
cdiag.                     obio_P(k,3),obio_P(k,4),k=1,kdm)

cdiag     write (*,104) nstep,i,j,
cdiag.' aftrupdate dpth     diat      chlo        cyan',
cdiag.'      cocc      herb',
cdiag.  (k,p(i,j,k+1)/onem,obio_P(k,5),obio_P(k,6),obio_P(k,7),
cdiag.                       obio_P(k,8),obio_P(k,9),k=1,kdm)

cdiag     write (*,104) nstep,i,j,
cdiag.' aftrupdate dpth     nitr_det  sili_det    iron_det',
cdiag.'      doc       dic ',
cdiag.    (k,p(i,j,k+1)/onem,det(k,1),det(k,2),det(k,3),
cdiag.                       car(k,1),car(k,2),k=1,kdm)
cdiag     endif


       !------------------------------------------------------------

      if (vrbos) then
       print*, 'OBIO TENDENCIES, 1-16, 1,7'
       do k=1,1
        write(*,'(16(e9.2,1x))')((rhs(k,nt,ll),ll=1,16),nt=1,7)
         print*, 'OBIO TENDENCIES, 1-16, 8,14'
        write(*,'(16(e9.2,1x))')((rhs(k,nt,ll),ll=1,16),nt=8,ntrac-1)
       enddo
      endif

!      if (diagno_bio) then
!        do k=1,kdm
!#ifdef OBIO_ON_GARYocean
!        write(iu_tend,'(4i5,6e12.4)')
!     .      nstep,i,j,k,p1d(k+1),rhs(k,14,5),rhs(k,14,10)
!     .                 ,rhs(k,14,14),rhs(k,14,15),rhs(k,14,16)
!#else
!        write(iu_tend,'(4i5,6e12.4)')
!     .      nstep,i,j,k,p(i,j,k+1)/onem,rhs(k,14,5),rhs(k,14,10)
!     .                 ,rhs(k,14,14),rhs(k,14,15),rhs(k,14,16)
!#endif
!        enddo
!      endif

       !------------------------------------------------------------
       !update 3d tracer array
       do k=1,kmax
        do nt=1,ntyp+n_inert
         tracer(i,j,k,nt)=obio_P(k,nt)
        enddo
        do nt=1,ndet
         tracer(i,j,k,ntyp+n_inert+nt)=det(k,nt)
        enddo
        do nt=1,ncar
         tracer(i,j,k,ntyp+n_inert+ndet+nt)=car(k,nt)
        enddo
#ifdef TRACERS_Alkalinity
        do nt=1,nalk
         tracer(i,j,k,ntyp+n_inert+ndet+ncar+nt)=alk1d(k)
        enddo
#endif
        !update avgq and gcmax arrays
        avgq(i,j,k)=avgq1d(k)
        gcmax(i,j,k)=gcmax1d(k)
        tirrq3d(i,j,k)=tirrq(k)
       enddo !k

!maybe we need to smooth dic here....

#ifdef OBIO_ON_GARYocean
      !update trmo etc arrays
       do k=1,kmax
       do nt=1,ntrac

          dtr = tracer(i,j,k,nt) * trmo_unit_factor(k,nt) 
     .        - trmo(i,j,k,nt)

        if (dtr.lt.0) then
          ftr = -dtr/trmo(i,j,k,nt)
          txmo(i,j,k,nt)=txmo(i,j,k,nt)*(1.-ftr)
          tymo(i,j,k,nt)=tymo(i,j,k,nt)*(1.-ftr)
          tzmo(i,j,k,nt)=tzmo(i,j,k,nt)*(1.-ftr)       
        endif

        trmo(i,j,k,nt)=trmo(i,j,k,nt)+dtr

       enddo
       enddo
#endif

       ihra(i,j)=ihra_ij

       !compute total chlorophyl at surface layer
       tot_chlo(i,j)=0.
       do nt=1,nchl
          tot_chlo(i,j)=tot_chlo(i,j)+obio_P(1,nnut+nt)
       enddo
       if (vrbos) then
          write(*,'(/,a,3i5,e12.4)')
     .       'obio_model, tot_chlo= ',nstep,i,j,tot_chlo(i,j)
       endif

       !compute total primary production per day
        if (hour_of_day.eq.1) then
          pp2tot_day(i,j)=0.
        else
          do nt=1,nchl
          do k=1,kdm
            if (p1d(kdm+1).gt.200.)     !total depth > 200m
     .        pp2tot_day(i,j)=pp2tot_day(i,j)+pp2_1d(k,nt)
          enddo
          enddo
       endif
cdiag  if (vrbos) then
cdiag  do nt=1,nchl
cdiag  do k=1,kdm
cdiag    write(*,'(a,7i5,3e12.4)')'PRIMARY PRODUCTIVITY ',
cdiag.       nstep,day_of_month,hour_of_day,i,j,k,nt,
cdiag.       p1d(kdm+1),pp2_1d(k,nt),pp2tot_day(i,j)
cdiag  enddo
cdiag  enddo
cdiag  endif

       !update pCO2 array
       pCO2(i,j)=pCO2_ij


!diagnostics
#ifdef OBIO_ON_GARYocean
       OIJ(I,J,IJ_nitr) = OIJ(I,J,IJ_nitr) + tracer(i,j,1,1) ! surf ocean nitrates
       OIJ(I,J,IJ_amm) = OIJ(I,J,IJ_amm) + tracer(i,j,1,2) ! surf ocean nitrates
       OIJ(I,J,IJ_sil) = OIJ(I,J,IJ_sil) + tracer(i,j,1,3) ! surf ocean nitrates
       OIJ(I,J,IJ_iron) = OIJ(I,J,IJ_iron) + tracer(i,j,1,4) ! surf ocean nitrates

       OIJ(I,J,IJ_diat) = OIJ(I,J,IJ_diat) + tracer(i,j,1,5) ! surf ocean diatoms
       OIJ(I,J,IJ_chlo) = OIJ(I,J,IJ_chlo) + tracer(i,j,1,6) ! surf ocean diatoms
       OIJ(I,J,IJ_cyan) = OIJ(I,J,IJ_cyan) + tracer(i,j,1,7) ! surf ocean diatoms
       OIJ(I,J,IJ_cocc) = OIJ(I,J,IJ_cocc) + tracer(i,j,1,8) ! surf ocean diatoms
       OIJ(I,J,IJ_herb) = OIJ(I,J,IJ_herb) + tracer(i,j,1,9) ! surf ocean diatoms

       OIJ(I,J,IJ_doc) = OIJ(I,J,IJ_doc) + tracer(i,j,1,14) ! surf ocean doc
       OIJ(I,J,IJ_dic) = OIJ(I,J,IJ_dic) + tracer(i,j,1,15) ! surf ocean dic
       OIJ(I,J,IJ_pCO2) = OIJ(I,J,IJ_pCO2) + pCO2(i,j)*(1.-oRSI(i,j)) ! surf ocean pco2

#ifdef TRACERS_Alkalinity
       OIJ(I,J,IJ_alk) = OIJ(I,J,IJ_alk) + tracer(i,j,1,16)    ! surf ocean alkalinity
#else
       OIJ(I,J,IJ_alk) = OIJ(I,J,IJ_alk) + alk(i,j,1)          ! surf ocean alkalinity
#endif
#endif

!       if (diagno_bio) then
!#ifdef OBIO_ON_GARYocean
!       if (kpl(i,j).gt.0) then    !??/why kpl(2,90)=0???
!#endif
!         write(iu_pco2,'(3i7,22e12.4)')
!     .     nstep,i,j
!     .    ,temp1d(1),saln1d(1),dp1d(1),car(1,2),pCO2(i,j)
!     .    ,covice_ij,obio_P(1,1:ntyp),det(1,1:ndet),car(1,1)
!#ifdef OBIO_ON_GARYocean
!     .    ,p1d(kpl(i,j))
!#else
!     .    ,dpmixl(i,j,mm)/onem
!#endif
!     .    ,alk1d(1),pp2tot_day(i,j)
!#ifdef OBIO_ON_GARYocean
!       endif   !kpl>0
!#endif
!       endif   !diagno_bio
 
#ifdef OBIO_ON_GARYocean
      endif   !if focean>0
#endif

 1000 continue
c$OMP END PARALLEL DO

!      if (diagno_bio) then
!        call closeunit(iu_pco2)
!        call closeunit(iu_tend)
!      endif     ! diagno_bio


#ifdef OBIO_ON_GARYocean
      !gather_tracer
      call pack_data( ogrid,  tracer, tracer_glob )
#endif


      call obio_trint

      return

      end subroutine obio_model
