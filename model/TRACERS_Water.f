#include "rundeck_opts.h"

#ifdef TRACERS_ON
!@sum  TRACERS_Water: tracer-dependent routines for water mass tracers
!@+    Routines included:
!@+      Those that MUST EXIST for all tracers: 
!@+        Diagnostic specs: init_tracer
!@+        Tracer initialisation + sources: tracer_ic, set_tracer_source
!@+        Entry points: daily_tracer
!@+      Those that are unique for specific tracers: 
!@+   
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0

      SUBROUTINE init_tracer
!@sum init_tracer initializes water tracer attributes and diagnostics
!@auth J. Lerner
!@calls sync_param, SET_TCON
      use DAGCOM, only: ia_src,ia_12hr,ir_log2
      USE MODEL_COM, only: dtsrc,nisurf
      USE TRACER_COM
      USE TRACER_DIAG_COM
      USE PARAM
      implicit none
      integer :: l,k,n
      character*20 sum_unit(ntm),inst_unit(ntm)   ! for conservation
      logical :: qcon(KTCON-1), qsum(KTCON-1), T=.TRUE. , F=.FALSE.
      character*50 :: unit_string

C**** Set some diags that are the same regardless
      call set_generic_tracer_diags

C**** Tracer sources and sinks
C**** Defaults for jls (sources, sinks, etc.)
C**** This needs to be 'hand coded' depending on circumstances
      do k=1,ktajls  ! max number of sources and sinks
        jgrid_jls(k) = 1
        ia_jls(k) = ia_src
        scale_jls(k) = 1./DTsrc
      end do
      jls_index(:) = 0
        
      k = 0
      n = n_Water
        k = k + 1 
        jls_source(1,n)=k
        sname_jls(k) = 'Surface_source_'//trname(n)
        lname_jls(k) = 'EVAPORATION OF '//trname(n)
        jls_index(k) = n
        jls_ltop(k) = 1
        jls_power(k) = 0.
        units_jls(k) = unit_string(jls_power(k),' kg/s')

      if (k.gt. ktajls) then
        write (6,*) 
     &   'tjl_defs: Increase ktajls=',ktajls,' to at least ',k
        stop 'ktajls too small'
      end if

C**** TAIJN: fill in if necessary

c      if (k .gt. ktaij) then
c        write (6,*) 
c     &   'tij_defs: Increase ktaij=',ktaij,' to at least ',k
c        stop 'ktaij too small'
c      end if

C**** Tracer sources and sinks
C**** Defaults for ijts (sources, sinks, etc.)
C**** This needs to be 'hand coded' depending on circumstances
      k = 1
        ijts_source(1,n_Water) = k
        ijts_index(k) = n_Water
        ia_ijts(k) = ia_src   
        lname_ijts(k) = 'Water Evaporative SOURCE'
        sname_ijts(k) = 'water_evap'
        ijts_power(k) = 0.
        units_ijts(k) = unit_string(ijts_power(k),' kg/s*m^2')
        scale_ijts(k) = 10.**(-ijts_power(k))/DTsrc

      if (k .gt. ktaijs) then
        write (6,*)'ijt_defs: Increase ktaijs=',ktaijs,' to at least ',k
        stop 'ktaijs too small'
      end if

C**** Initialize conservation diagnostics
C**** To add a new conservation diagnostic:
C****       Set up a QCON, and call SET_TCON to allocate array numbers,
C****       set up scales, titles, etc. 
C**** QCON denotes when the conservation diags should be accumulated
C**** QSUM says whether that diag is to be used in summation (if the
C****      routine DIAGCTB is used, this must be false).
C**** 1:NPTS+1 ==> INST,  DYN,   COND,   RAD,   PREC,   LAND,  SURF,
C****            FILTER,STRDG/OCEAN, DAILY, OCEAN1, OCEAN2,
C**** First 12 are standard for all tracers and GCM
      QCON=(/ t,                                           !instant.
     *        T,  T,  F,  F,  T,  T,  T,  T,  F,  F,  F,   !2-12 (npts)
     *        F,  F,  F,  F,  F,  F,  F,  F,  F,  F,       !13-22
     *        F,  F,  F,  F,  F,  F,  F,  F,  F,  F  /)    !21-ktcon-1
      QSUM=(/ f,                                           !instant.
     *        T,  T,  F,  F,  T,  T,  T,  T,  F,  F,  F,   !2-12 (npts)
     *        F,  F,  F,  F,  F,  F,  F,  F,  F,  F,       !13-22
     *        F,  F,  F,  F,  F,  F,  F,  F,  F,  F  /)    !21-ktcon-1
      do n=1,ntm
        kt_power_inst(n)   = ntm_power(n)+2
        kt_power_change(n) = ntm_power(n)-4
        scale_inst(n)   = 10d0**(-kt_power_inst(n))
        scale_change(n) = 10d0**(-kt_power_change(n))
        inst_unit(n) = unit_string(kt_power_inst(n),  ' kg/m^2)')
         sum_unit(n) = unit_string(kt_power_change(n),' kg/s/m^2)')
      end do
      N = n_Water
      itcon_mc(n)=13
      qcon(itcon_mc(n))=.true.  ; conpts(1) = 'MOIST CONV'
      qsum(itcon_mc(n)) = .false.
      itcon_ss(n)=14
      qcon(itcon_ss(n))=.true.  ; conpts(2) = 'LS COND'
      qsum(itcon_ss(n)) = .false.
      CALL SET_TCON(QCON,TRNAME(N),QSUM,inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)

C**** print out total tracer diagnostic array size
      WRITE (6,'(A14,2I8)') "KTACC=",KTACC

      return
      end subroutine init_tracer

      SUBROUTINE tracer_IC
!@sum tracer_IC initializes tracers when they are first switched on
!@auth Jean Lerner/Gavin Schmidt
      USE CONSTANT, only : rhow
      USE MODEL_COM, only: itime,im,jm,lm,q,wm,flice,fearth
      USE SOMTQ_COM, only : qmom
      USE GEOM, only: dxyp,bydxyp
      USE DYNAMICS, only: am,byam  ! Air mass of each box (kg/m^2)
      USE PBLCOM, only : npbl,trabl,qabl
      USE LANDICE_COM, only : trlndi,trsnowli,snowli,trli0
      USE LANDICE, only : ace1li,ace2li
      USE SEAICE_COM, only : rsi,msi,snowi,trsi,trsi0,ssi
      USE SEAICE, only : xsi,ace1i
      USE LAKES_COM, only : trlake,mwl,mldlk,flake
      USE GHYCOM, only : trbare,trvege,trsnowbv,wbare,wvege,snowbv,afb
      USE TRACER_COM, only : ntm,trm,trmom,itime_tr0,trname,needtrs,trwm
     *     ,trw0
      USE FLUXES, only : gtracer
      IMPLICIT NONE
      INTEGER i,n,l,j,ipbl,it
      REAL*8 conv

      do n=1,ntm

      if (itime.eq.itime_tr0(n)) then
      select case (trname(n)) 
          
        case default            
c          write(6,*) 'In TRACER_IC:',trname(n),' does not exist '
          stop "TRACER_IC"
  
        case ('Water')
          do l=1,lm
          do j=1,jm
            trm(:,j,l,n) =  q(:,j,l)*am(l,:,j)*dxyp(j)
            trwm(:,j,l,n)= wm(:,j,l)*am(l,:,j)*dxyp(j)
C**** for gradients defined on air mass
            do i=1,im
              trmom(:,i,j,l,n) = qmom(:,i,j,l)*am(l,i,j)*dxyp(j)
            end do
C**** for gradients defined on water mass (should be an option?)
c           trmom(:,:,j,l,n) = 0.
          end do; enddo
          do i=2,im
            trm(i, 1,:,n) = trm(1, 1,:,n) !poles
            trm(i,jm,:,n) = trm(1,jm,:,n) !poles
            trwm(i, 1,:,n)= trwm(1, 1,:,n) !poles
            trwm(i,jm,:,n)= trwm(1,jm,:,n) !poles
            trmom(:,i, 1,:,n)=0.
            trmom(:,i,jm,:,n)=0.
          enddo

          do j=1,jm
          do i=1,im
C**** lakes
            if (flake(i,j).gt.0) then
              trlake(n,1,i,j)=trw0(n)*mldlk(i,j)*rhow*flake(i,j)*dxyp(j)
              trlake(n,2,i,j)=trw0(n)*mwl(i,j)-trlake(n,1,i,j)
              gtracer(n,1,i,j)=trw0(n)
            end if
c**** ice 
            if (rsi(i,j).gt.0) then
              trsi(n,1,i,j)=trsi0(n)*
     *             (xsi(1)*(snowi(i,j)+ace1i)-ssi(1,i,j))
              trsi(n,2,i,j)=trsi0(n)*
     *             (xsi(2)*(snowi(i,j)+ace1i)-ssi(2,i,j))
              trsi(n,3,i,j)=trsi0(n)*(xsi(3)*msi(i,j)-ssi(3,i,j))
              trsi(n,4,i,j)=trsi0(n)*(xsi(4)*msi(i,j)-ssi(4,i,j))
              gtracer(n,2,i,j)=trsi0(n)
            else
              gtracer(n,2,i,j)=0.
            end if
c**** landice
            if (flice(i,j).gt.0) then
              trlndi(n,i,j)=trli0(n)*(ace1li+ace2li)
              trsnowli(n,i,j)=trli0(n)*snowli(i,j)
              gtracer(n,3,i,j)=trli0(n)
            else
              gtracer(n,3,i,j)=0.
            end if
c**** earth
            if (fearth(i,j).gt.0) then
              conv=rhow  ! convert from m to kg/m^2
              trbare  (n,:,i,j)=trw0(n)*wbare (:,i,j)*conv
              trvege  (n,:,i,j)=trw0(n)*wvege (:,i,j)*conv
              trsnowbv(n,1,i,j)=trw0(n)*snowbv(1,i,j)*conv
              trsnowbv(n,2,i,j)=trw0(n)*snowbv(2,i,j)*conv
              gtracer (n,4,i,j)=trw0(n)
            else
              gtracer(n,4,i,j)=0.
            end if
          end do
          end do

      end select

C**** Initialise pbl profile if necessary
      if (needtrs(n)) then
        do it=1,4
        do j=1,jm
        do ipbl=1,npbl
          trabl(ipbl,n,:,j,it) = trw0(n)*qabl(ipbl,:,j,it) 
        end do
        end do
        end do
      end if

      write(6,*) ' Tracer ',trname(n),' initialized at itime=',itime
      end if
      end do
#ifdef TRACERS_WATER
C**** Initialise ocean tracers if necessary
      call tracer_ic_ocean
#endif
C****
      end subroutine tracer_IC

      subroutine daily_tracer(iact)
!@sum daily_tracer is called once a day for tracers
!@auth Jean Lerner
C**** Note this routine must always exist (but can be a dummy routine)
      USE TRACER_COM, only : ntm,trname
      IMPLICIT NONE
      INTEGER n,iact

C**** Initialize tracers here to allow for tracers that 'turn on'
C****  at any time
      call tracer_IC

      return
C****
      end subroutine daily_tracer

      SUBROUTINE set_tracer_source
!@sum tracer_source calculates non-interactive sources for tracers
!@auth Jean Lerner/Gavin Schmidt
      USE MODEL_COM, only: itime
      USE TRACER_COM
      USE FLUXES, only : trsource
      implicit none
      integer :: n

C**** All sources are saved as kg/s
      do n=1,ntm

      if (itime.lt.itime_tr0(n)) cycle
      select case (trname(n))

      case default
!     write(6,*) ' Sources for ',trname(n),' tracer are not in this routine'
C****
C**** No non-interactive surface sources of Water
C****
      case ('Water')
        trsource(:,:,:,n)=0

      end select
      end do

C****
      END SUBROUTINE set_tracer_source

#endif

#ifdef TRACERS_WATER
C---SUBROUTINES FOR TRACER WET DEPOSITION-------------------------------------

      SUBROUTINE GET_COND_FACTOR(L,N,WMXTR,FCLOUD,FQ0,fq)
!@sum  GET_COND_FACTOR calculation of condensate fraction for tracers
!@+    within or below convective or large-scale clouds. Gas 
!@+    condensation uses Henry's Law if not freezing.
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436TdsM23 CLOUDCHCC and CLOUDCHEM subroutines)
c
C**** GLOBAL parameters and variables:
      USE CLOUDS, only: PL, TL, NTIX
      USE TRACER_COM, only: tr_RKD,tr_DHD,nWATER,nGAS,nPART,tr_wd_TYPE
      USE CONSTANT, only: TF, BYGASC, MAIR
c      
      IMPLICIT NONE
c      
C**** Local parameters and variables and arguments:
c
!@param BY298K unknown meaning for now (assumed= 1./298K)
!@var Ppas pressure at current altitude (in Pascal=kg/s2/m) 
!@var TFAC exponential coeffiecient of tracer condensation temperature
!@+   dependence (mole/joule)
!@var FCLOUD fraction of cloud available for tracer condensation
!@var SSFAC dummy variable (assumed units= kg water?)
!@var FQ            fraction of tracer that goes into condensate
!@var FQ0 [default] fraction of tracer that goes into condensate
!@var L index for altitude loop
!@var N index for tracer number loop
!@var WMXTR mixing ratio of water available for tracer condensation ( )?
!@var RKD dummy variable (= tr_RKD*EXP[ ])
      REAL*8, PARAMETER :: BY298K=3.3557D-3
      REAL*8 Ppas, tfac, ssfac, RKD
      REAL*8,  INTENT(IN) :: fq0, FCLOUD, WMXTR
      REAL*8,  INTENT(OUT):: fq
      INTEGER, INTENT(IN) :: L, N
c
C**** CALCULATE the fraction of tracer mass that becomes condensate:
c
      SELECT CASE(tr_wd_TYPE(NTIX(N)))                 
        CASE(nGAS)                                  ! gas tracer
          fq = 0.D0                                   ! frozen case
          IF(TL(L).ge.TF) THEN                        ! if not frozen then: 
            Ppas = PL(L)*1.D2                         ! pressure to pascals
            tfac = (1.D0/TL(L) - BY298K)*BYGASC
            IF(tr_DHD(NTIX(N)).ne.0.D0) THEN
              RKD=tr_RKD(NTIX(N))*DEXP(-tr_DHD(NTIX(N))*tfac)
            ELSE  
              RKD=tr_RKD(NTIX(N))
            END IF
c           clwc=WMXTR*MAIR*1.D-3*Ppas*BYGASC/(TL(L)*FCLOUD)
c           ssfac=RKD*GASC*TL(L)*clwc   ! Henry's Law
            ssfac=RKD*WMXTR*MAIR*1.D-3*Ppas/FCLOUD
            fq=ssfac / (1.D0 + ssfac)
          END IF
        CASE(nWATER)                                ! water tracer
          fq = FQ0                                  
        CASE(nPART)                                 ! particulate tracer
          fq = 0.D0                                   ! temporarily zero.
c NOTE 1: Dorothy has some code that will be put here to condense 
c aerosols. GSF 1/4/02
c
c NOTE 2:  Really, any aerosol 'formation', (meaning the production of
c an aerosol tracer due to cloud chemistry, or any flux among tracers),
c should be done elsewhere, like in a chemistry section of the model.
c But if it is impossible to supply that section with the variables
c needed from the wet deposition code, then the aerosol formation code
c should probably go here... If you add this, please make appropriate
c changes in the subroutine's name/summary above. GSF 1/4/02. 
        CASE DEFAULT                                ! error
          STOP 'tr_wd_TYPE(NTIX(N)) out of range in SCAVENGE_TRACER'            
      END SELECT
c      
      RETURN
      END SUBROUTINE GET_COND_FACTOR


      SUBROUTINE GET_PREC_FACTOR(N,BELOW_CLOUD,CM,FCLD,FQ0,fq)
!@sum  GET_PREC_FACTOR calculation of the precipitation scavenging fraction
!@+    for tracers WITHIN large scale clouds. Current version uses the
!@+    first order removal rate based on [Giorgi and Chameides, 1986], for
!@+    gaseous and particulate tracers. 
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436TdsM23 RAINOUT subroutine)
c
C**** GLOBAL parameters and variables:
      USE CLOUDS, only: NTIX
      USE MODEL_COM, only: dtsrc
      USE TRACER_COM, only: nWATER, nGAS, nPART, tr_wd_TYPE
c
      IMPLICIT NONE
c      
C**** Local parameters and variables and arguments:
!@var FQ            tracer fraction scavenged into precipitation
!@var FQ0 [default] tracer fraction scavenged into precipitation
!@var FCLD cloud fraction
!@var N index for tracer number loop
!@var CM conversion rate for large cloud water content
!@var BELOW_CLOUD logical- is the current level below cloud?
      LOGICAL, INTENT(IN) :: BELOW_CLOUD
      INTEGER, INTENT(IN) :: N
      REAL*8,  INTENT(IN) :: FQ0, FCLD, CM
      REAL*8,  INTENT(OUT):: FQ
c
      SELECT CASE(tr_wd_TYPE(NTIX(N)))                 
        CASE(nGAS)                                ! gas
          IF(BELOW_CLOUD) THEN
            fq = 0.D0
          ELSE
c           minus preserves FPRT sign convention in LSCOND
            fq = -(CM*DTsrc*(EXP(-CM*DTsrc)- 1.0)*FCLD)
          END IF
        CASE(nWATER)                              ! water/original method
          fq = FQ0                               
        CASE(nPART)                               ! aerosols
          IF(BELOW_CLOUD) THEN
            fq = 0.D0
          ELSE
c           minus preserves FPRT sign convention in LSCOND
            fq = -(CM*DTsrc*(EXP(-CM*DTsrc)- 1.0)*FCLD)
          END IF
        CASE DEFAULT                              ! error
          STOP 'tr_wd_TYPE(NTIX(N)) out of range in GET_FPRT'                              
      END SELECT
c
      RETURN
      END SUBROUTINE GET_PREC_FACTOR


      SUBROUTINE GET_WASH_FACTOR(N,b_beta_DT,PREC,fq)
!@sum  GET_WASH_FACTOR calculation of the fraction of tracer 
!@+    scavanged by precipitation below convective clouds ("washout"). 
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436TdsM23 CWASH and WASH_EVAP routines)
c
C**** GLOBAL parameters and variables:
      USE TRACER_COM, only: nWATER, nGAS, nPART, tr_wd_TYPE
      USE CLOUDS, only: NTIX
c
      IMPLICIT NONE
c 
C**** Local parameters and variables and arguments:
!@var FQ fraction of tracer scavenged by below-cloud precipitation
!@param rc_wash aerosol washout rate constant (mm-1)
!@var PREC precipitation amount from layer above for washout (mm)
!@var b_beta_DT precipitating grid box fraction from lowest percipitating
!@+   layer. The name was chosen to correspond to Koch et al. p. 23,802.
!@var N index for tracer number loop
      INTEGER, INTENT(IN) :: N
      REAL*8, INTENT(OUT):: FQ
      REAL*8, INTENT(IN) :: PREC, b_beta_DT
      REAL*8, PARAMETER :: rc_wash = 1.D-1
C
      SELECT CASE(tr_wd_TYPE(NTIX(N)))                 
        CASE(nGAS)                                ! gas
          fq = 0.D0
        CASE(nWATER)                              ! water/original method
          fq = 0.D0                              
        CASE(nPART)                               ! aerosols
          fq = b_beta_DT*(DEXP(-PREC*rc_wash)-1.)
        CASE DEFAULT                              ! error
          STOP 'tr_wd_TYPE(NTIX(N)) out of range in WASHOUT_TRACER'                              
      END SELECT   
c      
      RETURN
      END SUBROUTINE GET_WASH_FACTOR


      SUBROUTINE GET_EVAP_FACTOR(N,FQ0,fq)
!@sum  GET_EVAP_FACTOR calculation of the evaporation fraction for tracers.
!@auth Dorothy Koch (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on CB436TdsM23 EVAPD and WASH_EVAP routines)
c
C**** GLOBAL parameters and variables:
      USE CLOUDS, only: NTIX
      USE TRACER_COM, only: tr_evap_fact, tr_wd_TYPE
c
      IMPLICIT NONE
c      
C**** Local parameters and variables and arguments:
!@var FQ            fraction of tracer evaporated
!@var FQ0 [default] fraction of tracer evaporated
!@var N index for tracer number loop
      INTEGER, INTENT(IN) :: N
      REAL*8,  INTENT(OUT):: FQ
      REAL*8,  INTENT(IN) :: FQ0
c
      fq=FQ0*tr_evap_fact(tr_wd_TYPE(NTIX(N)))
c
      RETURN
      END SUBROUTINE GET_EVAP_FACTOR 
#endif
