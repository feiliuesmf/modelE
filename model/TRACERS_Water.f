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
      logical :: qcon(KTCON-1), T=.TRUE. , F=.FALSE.
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

C**** TAIJN

      if (k .gt. ktaij) then
        write (6,*) 
     &   'tij_defs: Increase ktaij=',ktaij,' to at least ',k
        stop 'ktaij too small'
      end if

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
C**** 1:NPTS+1 ==> INST,  DYN,   COND,   RAD,   PREC,   LAND,  SURF,
C****            FILTER,STRDG/OCEAN, DAILY, OCEAN1, OCEAN2,
C**** First 12 are standard for all tracers and GCM
      QCON=(/ t,                                           !instant.
     *        T,  T,  F,  T,  T,  T,  T,  T,  F,  F,  F,   !2-12 (npts)
     *        F,  F,  F,  F,  F,  F,  F,  F,  F,  F,       !13-22
     *        F,  F,  F,  F,  F,  F,  F,  F,  F,  F  /)   !21-ktcon-1
      do n=1,ntm
        kt_power_inst(n)   = ntm_power(n)+2
        kt_power_change(n) = ntm_power(n)-4
        scale_inst(n)   = 10d0**(-kt_power_inst(n))
        scale_change(n) = 10d0**(-kt_power_change(n))
        inst_unit(n) = unit_string(kt_power_inst(n),  ' kg/m^2)')
         sum_unit(n) = unit_string(kt_power_change(n),' kg/s/m^2)')
      end do
      N = n_Water
      CALL SET_TCON(QCON,TRNAME(N),inst_unit(n),
     *     sum_unit(n),scale_inst(n),scale_change(n), N,CONPTs)

C**** print out total tracer diagnostic array size
      WRITE (6,'(A14,2I8)') "KTACC=",KTACC

      return
      end subroutine init_tracer

      SUBROUTINE tracer_IC
!@sum tracer_IC initializes tracers when they are first switched on
!@auth Jean Lerner/Gavin Schmidt
      USE CONSTANT, only : rhow
      USE MODEL_COM, only: itime,im,jm,lm,q,wm,flice
      USE SOMTQ_COM, only : qmom
      USE GEOM, only: dxyp,bydxyp
      USE DYNAMICS, only: am,byam  ! Air mass of each box (kg/m^2)
      USE PBLCOM, only : npbl,trabl,qabl
      USE LANDICE_COM, only : trlndi,trsnowli,snowli,trli0
      USE LANDICE, only : ace1li,ace2li
      USE SEAICE_COM, only : rsi,msi,snowi,trsi,trsi0
      USE SEAICE, only : xsi,ace1i
      USE LAKES_COM, only : trlake,mwl,mldlk,flake
      USE TRACER_COM, only : ntm,trm,trmom,itime_tr0,trname,needtrs,trwm
     *     ,trw0
      IMPLICIT NONE
      INTEGER i,n,l,j,ipbl,it

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
C**** lakes
          do j=1,jm
          do i=1,im
            if (flake(i,j).gt.0) THEN
              trlake(:,1,i,j)=trw0(:)*mldlk(i,j)*rhow*flake(i,j)*dxyp(j)
              trlake(:,2,i,j)=trw0(:)*mwl(i,j)-trlake(:,1,i,j)
            end if
          end do
          end do
C**** ice 
          do j=1,jm
          do i=1,im
            if (rsi(i,j).gt.0) then
              TRSI(:,1,I,J)=TRSI0(:)*XSI(1)*(SNOWI(I,J)+ACE1I)
              TRSI(:,2,I,J)=TRSI0(:)*XSI(2)*(SNOWI(I,J)+ACE1I)
              TRSI(:,3,I,J)=TRSI0(:)*XSI(3)*MSI(I,J)
              TRSI(:,4,I,J)=TRSI0(:)*XSI(4)*MSI(I,J)
            end if
          end do
          end do
C**** landice
          do j=1,jm
          do i=1,im
            if (flice(i,j).gt.0) then
              TRLNDI(:,I,J)=TRLI0(:)*(ACE1LI+ACE2LI)
              TRSNOWLI(:,I,J)=TRLI0(:)*SNOWLI(I,J)
            end if
          end do
          end do

      end select

C**** Initialise pbl profile if necessary
      if (needtrs(n)) then
        do it=1,4
        do j=1,jm
        do ipbl=1,npbl
          trabl(ipbl,n,:,j,it) = qabl(ipbl,:,j,it) 
        end do
        end do
        end do
      end if

      write(6,*) ' Tracer ',trname(n),' initialized at itime=',itime
      end if
      end do
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
      USE FLUXES, only : trsource,tot_trsource
      implicit none
      integer :: ns,n

C**** All sources are saved as kg/s
      do n=1,ntm

      if (itime.lt.itime_tr0(n)) cycle
      select case (trname(n))

      case default
!     write(6,*) ' Sources for ',trname(n),' tracer are not in this routine'
C****
C**** No un-interactive surface sources of Water
C****
      case ('Water')
        trsource(:,:,:,n)=0

      end select

C**** tot_trsource is required for PBL calculations
      tot_trsource(:,:,n) = 0.
      do ns=1,ntsurfsrc(n)
        tot_trsource(:,:,n) = tot_trsource(:,:,n)+trsource(:,:,ns,n)
      end do

      end do
C****
      END SUBROUTINE set_tracer_source

#endif
