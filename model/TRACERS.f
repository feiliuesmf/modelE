#include "rundeck_opts.h"

#ifdef TRACERS_ON
!@sum  TRACERS: generic tracer routines used for all tracers
!@+    Routines included: 
!@+      Generic diags: set_generic_tracer_diags
!@+      Apply previously set sources: apply_tracer_sources 
!@+      Radioactive Decay: tdecay 
!@+      Gravitaional Settling: trgrav
!@+      Check routine: checktr
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0

      SUBROUTINE set_generic_tracer_diags
!@sum init_tracer initializes trace gas attributes and diagnostics
!@auth J. Lerner
!@calls sync_param, SET_TCON
      use DAGCOM, only: ia_src,ia_12hr,ir_log2
      USE MODEL_COM, only: dtsrc,nisurf
      USE TRACER_COM
      USE TRACER_DIAG_COM
      USE CONSTANT, only: mair,sday
      USE PARAM
      implicit none
      integer :: l,k,n
      character*20 sum_unit(ntm),inst_unit(ntm)   ! for conservation
      character*10 CMR
      logical :: qcon(KTCON-1), T=.TRUE. , F=.FALSE.
!@var to_volume_MixRat: For printout of tracer concentration
!@+   to_volume_MixRat=1: printout is in Volume Mixing Ratio
!@+   to_volume_MixRat=0: printout is in Mass Mixing Ratio
      INTEGER :: to_volume_MixRat=0 
      character*50 :: unit_string
C**** Get itime_tr0 from rundeck if it exists
      call sync_param("itime_tr0",itime_tr0,ntm)
C**** Get to_volume_MixRat from rundecks if it exists
      call sync_param("to_volume_MixRat",to_volume_MixRat)

C**** Get factor to convert from mass mixing ratio to volume mr
      if (to_volume_MixRat .eq.1) then
        MMR_to_VMR(:) = mair/tr_mm(:)
        cmr = ' V/V air'
      else
        MMR_to_VMR(:) = 1.d0
        cmr = ' kg/kg air'
      endif
C****
C**** TAJLN(J,L,KQ,N)  (SUM OVER LONGITUDE AND TIME OF)
C****
C**** jlq_power Exponent associated with a physical process 
C****      (for printing tracers).   (ntm_power+jlq_power)
C**** jls_power Exponent associated with a source/sink (for printing)
C**** Defaults for JLN
      scale_jlq(:) = 1./DTsrc
      jgrid_jlq(:) = 1
      ia_jlq(:) = ia_src

C**** Tracer concentration
      do n=1,ntm
        k = 1        ! <<<<< Be sure to do this
        jlnt_conc = k
        sname_jln(k,n) = trim(trname(n))//'_CONCENTRATION' 
        lname_jln(k,n) = trim(trname(n))//' CONCENTRATION' 
        jlq_power(k) = 0.
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),cmr)
        scale_jlq(k) = 1.d0
        scale_jln(n) = MMR_to_VMR(n)

C**** Physical processes affecting tracers
C****   F (TOTAL NORTHWARD TRANSPORT OF TRACER MASS)  (kg)
        k = k + 1 
        jlnt_nt_tot = k
        sname_jln(k,n) = 'tr_nt_tot_'//trname(n)
        lname_jln(k,n) = 'TOTAL NORTHWARD TRANSPORT OF '//
     &     trim(trname(n))//' MASS'
        jlq_power(k) = 11.
        jgrid_jlq(k) = 2
C****   STM/SM (MEAN MERIDIONAL N.T. OF TRACER MASS)  (kg)
        k = k + 1 
        jlnt_nt_mm = k
        sname_jln(k,n) = 'tr_nt_mm_'//trname(n)
        lname_jln(k,n) = 'NORTHWARD TRANS. OF '//
     &     trim(trname(n))//' MASS BY MERIDIONAL CIRC.'
        jlq_power(k) = 10.
        jgrid_jlq(k) = 2
C****   F (TOTAL VERTICAL TRANSPORT OF TRACER MASS)  (kg)
        k = k + 1
        jlnt_vt_tot = k
        sname_jln(k,n) = 'tr_vt_tot_'//trname(n)
        lname_jln(k,n) = 'TOTAL VERTICAL TRANSPORT OF '//
     &     trim(trname(n))//' MASS'
        jlq_power(k) = 11.
C****   STM/SM (MEAN MERIDIONAL V.T. OF TRACER MASS)  (kg)
        k = k + 1
        jlnt_vt_mm = k
        sname_jln(k,n) = 'tr_vt_mm_'//trname(n)
        lname_jln(k,n) = 'VERTICAL TRANS. OF '//
     &     trim(trname(n))//' MASS BY MERIDIONAL CIRC.'
        jlq_power(k) = 10.
C****   TMBAR-TM (CHANGE OF TRACER MASS BY MOIST CONVEC)(kg)
        k = k + 1
        jlnt_mc = k
        sname_jln(k,n) = 'tr_mc_'//trname(n)
        lname_jln(k,n) = 'CHANGE OF '//
     &     trim(trname(n))//' MASS BY MOIST CONVECTION'
        jlq_power(k) = 10.
C****   TMBAR-TM (CHANGE OF TRACER MASS BY Large-scale CONDENSE)  (kg)
        k = k + 1
        jlnt_lscond = k
        sname_jln(k,n) = 'tr_lscond'//trname(n)
        lname_jln(k,n) ='CHANGE OF '//
     &     trim(trname(n))//' MASS BY LARGE-SCALE CONDENSE'
        jlq_power(k) = 10.
C****   TMBAR-TM (CHANGE OF TRACER MASS BY DRY CONVEC)  (kg)
        k = k + 1
        jlnt_turb = k
        sname_jln(k,n) = 'tr_turb_'//trname(n)
        lname_jln(k,n) = 'CHANGE OF '//
     &     trim(trname(n))//' MASS BY TURBULENCE/DRY CONVECTION'
        jlq_power(k) = 10.
      end do

      if (k.gt. ktajl) then
        write (6,*) 
     &   'tjl_defs: Increase ktajl=',ktajl,' to at least ',k
        stop 'ktajl too small'
      end if

C**** Construct UNITS string for output
      do n=1,ntm
      do k=2,ktajl
      units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),' kg/s')
      end do; end do

C**** CONTENTS OF TAIJLN(I,J,LM,N)  (SUM OVER TIME OF)
C****        TML (M*M * KG TRACER/KG AIR)
C**** Set defaults that are true for all tracers and layers
      ia_ijt    = ia_src
      ir_ijt(:) = ir_log2   !n
      ijtc_power(:) = ntm_power(:)+1   !n for concentration
      ijtm_power(:) = ntm_power(:)+4   !n for integrated mass
C**** Tracer concentrations (AIJLN)
      do n=1,ntm
      do l=1,lm
        write(sname_ijt(l,n),'(a,i2.2)') trim(TRNAME(n))//'_L_',l 
        write(lname_ijt(l,n),'(a,i2)')   trim(TRNAME(n))//' L ',l
        units_ijt(l,n) = unit_string(ijtc_power(n),cmr)
        scale_ijt(l,n) = MMR_to_VMR(n)*10.**(-ijtc_power(n))
      end do
      end do

C**** AIJN
C****     1  TM (SUM OVER ALL LAYERS) (M*M * KG TRACER/KG AIR)
C****     2  TRS (SURFACE TRACER CONC.) (M*M * KG TRACER/KG AIR)
C****     3  TM (SUM OVER ALL LAYERS) (M*M * KG TRACER)
      do n=1,ntm
C**** Summation of mass over all layers
      k = 1        ! <<<<< Be sure to do this
      tij_mass = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_Total_Mass'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//' Total Mass'
        units_tij(k,n) = unit_string(ijtm_power(n),' kg/m^2')
        scale_tij(k,n) = 10.**(-ijtm_power(n))
C**** Average concentration over layers
      k = k+1
      tij_conc = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_Average'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//' Average'
        units_tij(k,n) = unit_string(ijtc_power(n),cmr)
        scale_tij(k,n) = MMR_to_VMR(n)*10.**(-ijtc_power(n))
C**** Surface concentration
      k = k+1
      tij_surf = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_At_Surface'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//' At Surface'
        units_tij(k,n) = unit_string(ijtc_power(n),cmr)
        scale_tij(k,n)=MMR_to_VMR(n)*10.**(-ijtc_power(n))/dble(NIsurf)
#ifdef TRACERS_WATER
C**** Tracers in precipitation
      k = k+1
      tij_prec = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_in_prcp'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' in Precip'
        units_tij(k,n)=unit_string(ijtc_power(n),cmr)
        scale_tij(k,n)=10.**(-ijtc_power(n))
C**** Tracers in evaporation
      k = k+1
      tij_evap = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_in_evap'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' in Evaporation'
        units_tij(k,n)=unit_string(ijtc_power(n),cmr)
        scale_tij(k,n)=10.**(-ijtc_power(n))
C**** Tracers in river runoff
      k = k+1
      tij_rvr = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_in_rvr'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' in River Outflow'
        units_tij(k,n)=unit_string(ijtc_power(n),cmr)
        scale_tij(k,n)=10.**(-ijtc_power(n))
C**** Tracers in sea ice
      k = k+1
      tij_seaice = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_in_ice'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' in Sea Ice'
        units_tij(k,n)=unit_string(ijtc_power(n),cmr)
        scale_tij(k,n)=10.**(-ijtc_power(n))
C**** Tracers conc. in ground component (ie. water or ice surfaces)
      k = k+1
      tij_grnd = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_at_Grnd'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' at Ground'
        units_tij(k,n)=unit_string(ijtc_power(n),cmr)
        scale_tij(k,n)=10.**(-ijtc_power(n))
#endif
      end do

      if (k .gt. ktaij) then
        write (6,*) 
     &   'tij_defs: Increase ktaij=',ktaij,' to at least ',k
        stop 'ktaij too small'
      end if

      end SUBROUTINE set_generic_tracer_diags


      SUBROUTINE apply_tracer_source(dtsurf)
!@sum apply_tracer_source adds non-interactive surface sources to tracers
!@auth Jean Lerner/Gavin Schmidt
      USE MODEL_COM, only : jm
      USE GEOM, only : imaxj
      USE TRACER_COM, only : ntm,trm,trmom,ntsurfsrc
      USE QUSDEF, only : mz,mzz
      USE FLUXES, only : trsource,trflux1,trsrfflx
      USE TRACER_DIAG_COM, only : taijs,tajls,ijts_source,jls_source
     *     ,itcon_surf
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: dtsurf
      INTEGER n,ns,naij,najl,j,i
      REAL*8, DIMENSION(JM) :: dtracer

C**** This is tracer independent coding designed to work for all
C**** surface sources.
C**** Note that tracer flux is added to first layer either implicitly
C**** in ATURB or explcitly in 'apply_fluxes_to_atm' call in SURFACE.

      do n=1,ntm

        trflux1(:,:,n) = 0.
C**** Non-interactive sources
        do ns=1,ntsurfsrc(n)
C**** diagnostics
          naij = ijts_source(ns,n)
          taijs(:,:,naij) = taijs(:,:,naij) + trsource(:,:,ns,n)*dtsurf
          najl = jls_source(ns,n)
          do j=1,jm
            tajls(j,1,najl) = tajls(j,1,najl)+sum(trsource(:,j,ns,n))
     *           *dtsurf
            dtracer(j)=0.
            do i=1,imaxj(j)
              dtracer(j)=dtracer(j)+trsource(i,j,ns,n)*dtsurf
            end do
          end do
          call DIAGTCB(dtracer,itcon_surf(ns,n),n)
C**** trflux1 is total flux into first layer
          trflux1(:,:,n) = trflux1(:,:,n)+trsource(:,:,ns,n)
        end do
C**** Interactive sources
C**** diagnostics
c        naij = ijts_source(ns,n)  ????
c        taijs(:,:,naij) = taijs(:,:,naij) + trsrfflx(:,:,n)*dtsurf
c        najl = jls_source(ns,n)   ????
c        do j=1,jm
c          tajls(j,1,najl) = tajls(j,1,najl)+sum(trsrfflx(:,j,n))
c     *         *dtsurf
c        end do
c        call DIAGTCA(itcon_surf(ns,n),n)  ????

C**** Accumulate interactive sources as well
        trflux1(:,:,n) = trflux1(:,:,n)+trsrfflx(:,:,n)

C**** modify vertical moments        
        trmom( mz,:,:,1,n) = trmom( mz,:,:,1,n)-1.5*trflux1(:,:,n)
     *       *dtsurf
        trmom(mzz,:,:,1,n) = trmom(mzz,:,:,1,n)+0.5*trflux1(:,:,n)
     *       *dtsurf
      end do
C****
      RETURN
      END SUBROUTINE apply_tracer_source


      SUBROUTINE TDECAY
!@sum TDECAY decays radioactive tracers every source time step
!@auth Gavin Schmidt/Jean Lerner
      USE MODEL_COM, only : im,jm,lm,itime,dtsrc
      USE TRACER_COM, only : ntm,trm,trmom,trdecy,itime_tr0
#ifdef TRACERS_WATER
     *     ,trwm
#endif
      USE TRACER_DIAG_COM, only : tajls,jls_decay,itcon_decay
      IMPLICIT NONE
      real*8, save, dimension(ntm) :: expdec = 1.
      real*8, dimension(im,jm,lm) :: told
      integer, save :: ifirst=1
      integer n,najl,j,l

      if (ifirst.eq.1) then               
        do n=1,ntm
          if (trdecy(n).gt.0.0) expdec(n)=exp(-trdecy(n)*dtsrc)
        end do
        ifirst = 0
      end if

      do n=1,ntm
        if (trdecy(n).gt.0. .and. itime.ge.itime_tr0(n)) then
C**** Atmospheric decay
          told(:,:,:)=trm(:,:,:,n)
#ifdef TRACERS_WATER
     *               +trwm(:,:,:,n)
          trwm(:,:,:,n)=expdec(n)*trwm(:,:,:,n)
#endif
          trm(:,:,:,n)=expdec(n)*trm(:,:,:,n)
          trmom(:,:,:,:,n)=expdec(n)*trmom(:,:,:,:,n)
          najl = jls_decay(n)
          do l=1,lm
          do j=1,jm
          tajls(j,l,najl)=tajls(j,l,najl)+sum(trm(:,j,l,n)
#ifdef TRACERS_WATER
     *               +trwm(:,j,l,n)
#endif
     *               -told(:,j,l))
          enddo 
          enddo
          call DIAGTCA(itcon_decay(n),n)
        end if
      end do
C****
      return
      end subroutine tdecay


      SUBROUTINE TRGRAV
!@sum TRGRAV gravitationally settles particular tracers 
!@auth Gavin Schmidt/Reha Cakmur
      USE CONSTANT, only : visc_air,grav
      USE MODEL_COM, only : im,jm,lm,itime,dtsrc,zatmo
      USE GEOM, only : imaxj
      USE SOMTQ_COM, only : mz,mzz,mzx,myz
      USE TRACER_COM, only : ntm,trm,trmom,itime_tr0,trradius,trpdens
      USE TRACER_DIAG_COM, only : tajls,jls_grav,itcon_grav
      USE FLUXES, only : trgrdep
      USE DYNAMICS, only : gz
      IMPLICIT NONE
      real*8, save, dimension(ntm) :: stokevdt = 0.
      integer, save :: ifirst=1
      real*8, dimension(im,jm,lm) :: told
      real*8 fgrfluxd,fgrfluxu
      integer n,najl,i,j,l

      if (ifirst.eq.1) then               
C**** Calculate settling velocity based on Stokes' Law using particle
C**** density and effective radius
        do n=1,ntm
          if (trradius(n).gt.0.0) stokevdt(n)=dtsrc*2.*grav*trpdens(n)
     *         *trradius(n)**2/(9.*visc_air)
        end do
        ifirst = 0
      end if

      do n=1,ntm
        if (trradius(n).gt.0. .and. itime.ge.itime_tr0(n)) then
C**** Gravitional settling 
          do l=1,lm
          do j=1,jm
          do i=1,imaxj(j)
            told(i,j,l)=trm(i,j,l,n)
C**** Calculate height differences using geopotential
            if (l.eq.1) then   ! layer 1 calc
              fgrfluxd=stokevdt(n)*grav/(gz(i,j,l)-zatmo(i,j))
              trgrdep(i,j,n)=fgrfluxd*trm(i,j,l,n)
            else               ! above layer 1
              fgrfluxd=stokevdt(n)*grav/(gz(i,j,l)-gz(i,j,l-1))
            end if
            if (l.lt.lm) then  ! below top layer
              fgrfluxu=stokevdt(n)*grav/(gz(i,j,l+1)-gz(i,j,l))
            else               ! top layer
              fgrfluxu=0.
            end if
            trm(i,j,l,n)=trm(i,j,l  ,n)*(1.-fgrfluxd)
     *                 + trm(i,j,l+1,n)*    fgrfluxu
            trmom(mz ,i,j,l,n)=trmom(mz ,i,j,l,n)*(1.-fgrfluxd)
            trmom(mzz,i,j,l,n)=trmom(mzz,i,j,l,n)*(1.-fgrfluxd)
            trmom(mzx,i,j,l,n)=trmom(mzx,i,j,l,n)*(1.-fgrfluxd)
            trmom(myz,i,j,l,n)=trmom(myz,i,j,l,n)*(1.-fgrfluxd)
          end do
          end do
          end do

          najl = jls_grav(n)
          do l=1,lm
          do j=1,jm
          tajls(j,l,najl)=tajls(j,l,najl)+sum(trm(:,j,l,n)-told(:,j,l))
          enddo 
          enddo
          call DIAGTCA(itcon_grav(n),n)
        end if
      end do
C****
      return
      end subroutine trgrav


      SUBROUTINE read_monthly_sources(iu,jdlast,tlca,tlcb,data)
!@sum Read in monthly sources and interpolate to current day
!@+   Calling routine must have the two lines:
!@+      data jdlast /0/
!@+      save jdlast
!@+   Input: iu, the fileUnit#; jdlast
!@+   Output: interpolated data array
!@auth Jean Lerner and others
      USE MODEL_COM, only: jday,im,jm,idofm=>JDmidOfM
      implicit none
      real*8 frac
      real*8 tlca(im,jm),tlcb(im,jm),data(im,jm)
      integer imon,iu,jdlast

      if (jdlast.EQ.0) then ! NEED TO READ IN FIRST MONTH OF DATA
        imon=1          ! imon=January
        if (jday.le.16)  then ! JDAY in Jan 1-15, first month is Dec
          call readt(iu,0,tlca,im*jm,tlca,12)
          rewind iu
        else            ! JDAY is in Jan 16 to Dec 16, get first month
  120     imon=imon+1
          if (jday.gt.idofm(imon) .AND. imon.le.12) go to 120
          call readt(iu,0,tlca,im*jm,tlca,imon-1)
          if (imon.eq.13)  rewind iu
        end if
      else              ! Do we need to read in second month?
        if (jday.ne.jdlast+1) then ! Check that data is read in daily
          if (jday.ne.1 .OR. jdlast.ne.365) then
            write(6,*)
     *      'Incorrect values in Tracer Source:JDAY,JDLAST=',JDAY,JDLAST
            stop
          end if
          imon=imon-12  ! New year
          go to 130
        end if
        if (jday.le.idofm(imon)) go to 130
        imon=imon+1     ! read in new month of data
        tlca(:,:) = tlcb(:,:)
        if (imon.eq.13) rewind iu
      end if
      call readt(iu,0,tlcb,im*jm,tlcb,1)
  130 continue
c**** Interpolate two months of data to current day
      frac = float(idofm(imon)-jday)/(idofm(imon)-idofm(imon-1))
      data(:,:) = tlca(:,:)*frac + tlcb(:,:)*(1.-frac)
      return
      end subroutine read_monthly_sources


      subroutine checktr(subr)
!@sum  CHECKTR Checks whether tracer variables are reasonable
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : teeny
      USE MODEL_COM, only : ls1,im,jm,lm,q,wm
      USE DYNAMICS, only : am
      USE GEOM, only : dxyp,imaxj
      USE TRACER_COM
      IMPLICIT NONE
      LOGICAL QCHECKO
      INTEGER I,J,L,N, imax,jmax,lmax
      REAL*8 relerr, errmax
!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

      do n=1,ntm
        CALL CHECK3(trm(1,1,1,n),IM,JM,LM,SUBR,trname(n)(1:3))
        CALL CHECK3(trmom(1,1,1,1,n),NMOM,IM,JM*LM,SUBR,
     *       'X'//trname(n)(1:2))
#ifdef TRACERS_WATER
        CALL CHECK3(trwm(1,1,1,n),IM,JM,LM,SUBR,'W'//trname(n)(1:2))
#endif

C**** check whether air mass is conserved
        
        if (trname(n).eq.'Air') then
          errmax = 0. ; lmax=1 ; imax=1 ; jmax=1 
          do l=1,lm
          do j=1,jm
          do i=1,imaxj(j)
            relerr=abs(trm(i,j,l,n)-am(l,i,j)*dxyp(j))/
     *           (am(l,i,j)*dxyp(j))
            if (relerr.gt.errmax) then
              lmax=l ; imax=i ; jmax=j ; errmax=relerr
            end if
          end do
          end do
          end do
          print*,"Relative error in air mass after ",subr,":",imax
     *         ,jmax,lmax,errmax,trm(imax,jmax,lmax,n),am(lmax,imax
     *         ,jmax)*dxyp(jmax)
        end if
      
#ifdef TRACERS_WATER
        if (trname(n).eq.'Water') then
          errmax = 0. ; lmax=1 ; imax=1 ; jmax=1 
          do l=1,lm
          do j=1,jm
          do i=1,imaxj(j)
            relerr=max(abs(trm(i,j,l,n)-q(i,j,l)*am(l,i,j)*dxyp(j))/
     *           (q(i,j,l)*am(l,i,j)*dxyp(j)+teeny), abs(trwm(i,j,l,n)
     *           -wm(i,j,l)*am(l,i,j)*dxyp(j))/max(trwm(i,j,l,n),wm(i,j
     *           ,l)*am(l,i,j)*dxyp(j)))
            if (relerr.gt.errmax .and. trwm(i,j,l,n).gt.1.) then
              lmax=l ; imax=i ; jmax=j ; errmax=relerr
            end if
          end do
          end do
          end do
          print*,"Relative error in water mass after ",subr,":",imax
     *         ,jmax,lmax,errmax,trm(imax,jmax,lmax,n),trwm(imax,jmax
     *         ,lmax,n),q(imax,jmax,lmax)*am(lmax,imax,jmax)*dxyp(jmax)
     *         ,wm(imax,jmax,lmax)*am(lmax,imax,jmax)*dxyp(jmax)
        end if
#endif
      end do

      return
      end subroutine checktr

#endif

      SUBROUTINE io_tracer(kunit,iaction,ioerr)
!@sum  io_tracer reads and writes tracer variables to file
!@auth Jean Lerner
!@ver  1.0
      USE MODEL_COM, only: ioread,iowrite,irsfic,irerun,lhead
      USE TRACER_COM
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
#ifdef TRACERS_WATER
      CHARACTER*80 :: HEADER, MODULE_HEADER = "TRACERW01"

      write (MODULE_HEADER(lhead+1:80),'(a,i2,a,a,i1,a,i2,a,i2,a)')
     *     'R8 TRM(im,jm,lm,',NTM,')',
     *     ',TRmom(',NMOM,',im,jm,lm,',NTM,'),trwm(im,jm,lm,',NTM,')' 
#else
      CHARACTER*80 :: HEADER, MODULE_HEADER = "TRACER01"

      write (MODULE_HEADER(lhead+1:80),'(a,i2,a,a,i1,a,i2,a)')
     *           'R8 TRM(im,jm,lm,',NTM,')',
     *  ',TRmom(',NMOM,',im,jm,lm,',NTM,')'

#endif

      SELECT CASE (IACTION)

      CASE (:IOWRITE) ! output to end-of-month restart file
        WRITE (kunit,err=10) MODULE_HEADER,TRM,TRmom
#ifdef TRACERS_WATER
     *     ,TRWM
#endif
      CASE (IOREAD:)          ! input from restart file
        SELECT CASE (IACTION)
        CASE (IRSFIC)   ! initial conditions
          READ (kunit)
        CASE (ioread,irerun) ! restarts
          READ (kunit,err=10) HEADER,TRM,TRmom
#ifdef TRACERS_WATER
     *       ,TRWM
#endif
          IF (HEADER(1:lhead).ne.MODULE_HEADER(1:lhead)) THEN
            PRINT*,"Discrepancy in module version",HEADER,MODULE_HEADER
            GO TO 10
          END IF
        END SELECT
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_tracer
