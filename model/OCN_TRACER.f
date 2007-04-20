#include "rundeck_opts.h"

#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
!@sum  OCN_TRACER: tracer-dependent routines for GISS Ocean tracers
!@+    Routines included:
!@+      Those that MUST EXIST for all tracers:
!@+        Diagnostic specs: init_tracer_ocean (if TRACERS_OCEAN true)
!@+        Tracer initialisation + sources: tracer_ic_ocean
!@+
!@auth Jean Lerner/Gavin Schmidt
!@ver  1.0

      subroutine tracer_ic_ocean
!@sum tracer_ic_ocean initialise ocean tracers
!@auth Gavin Schmidt
!@ver 1.0
      USE GEOM, only : dxyp
      USE MODEL_COM, only: itime
      USE TRACER_COM, only : itime_tr0, ntm, trname, trw0
      USE OCEAN, only : im,jm,lmo,focean,dxypo,mo,lmm,imaxj
#ifdef TRACERS_OCEAN
     *     ,trmo,txmo,tymo,tzmo,s0m,sxmo,symo,szmo,oc_tracer_mean
     *     ,trmo_glob,s0m_glob,mo_glob
#endif
      USE SEAICE, only : xsi,lmi
      USE STRAITS, only : nmst,msist,ssist
#ifdef TRACERS_OCEAN
     *     ,lmst,ist,jst,xst,yst,mmst,s0mst,sxmst,szmst,trmst,txmst
     *     ,tzmst
#endif
#ifdef TRACERS_WATER
     *     ,trsist
      USE FLUXES, only : gtracer
#endif
      USE FILEMANAGER, only : openunit,closeunit
      USE DOMAIN_DECOMP, only : grid, get, haveLatitude, GLOBALSUM,
     *     ESMF_BCAST
      USE DOMAIN_DECOMP, only : AM_I_ROOT, pack_data, unpack_data
      USE OCEAN, only : scatter_ocean, gather_ocean

      IMPLICIT NONE
      integer n,i,j,l,k,nst,i1,j1,i2,j2
      real*8 t01,t02
#ifdef TRACERS_SPECIAL_O18
      integer iu_O18ic,ip1,im1
      character*80 title
      real*4, dimension(im,jm,lmo) :: t0m4,tzm4
      real*8 fwm,trsum,wsum,tratio,afac,frac_tr
#endif
      real*8 :: OTRACJ(GRID%J_STRT:GRID%J_STOP)
      INTEGER :: J_0S, J_1S, J_0, J_1, J_0H, J_1H
      CALL GET(grid, J_STRT_SKP = J_0S, J_STOP_SKP = J_1S, J_STRT = J_0,
     *     J_STOP = J_1, J_STRT_HALO = J_0H, J_STOP_HALO = J_1H)

C**** Note that only sea ice related arrays are initialised if
C**** only TRACERS_WATER is true.

      do n=1,ntm

        if (itime.eq.itime_tr0(n)) then
        select case (trname(n))

        case default
#ifdef TRACERS_OCEAN
          trmo(:,:,:,n)=0.
          txmo(:,:,:,n)=0.
          tymo(:,:,:,n)=0.
          tzmo(:,:,:,n)=0.
#endif
#ifdef TRACERS_WATER
          do j=J_0,J_1
          do i=1,im
            if (focean(i,j).gt.0) gtracer(n,1,i,j)=0.
          end do
          end do
#endif
        case ('Water', 'H2O18', 'HDO')
#if (defined TRACERS_OCEAN) && (defined TRACERS_SPECIAL_O18)
C**** Open ic file for isotope tracers
          call openunit("H2O18ic",iu_O18ic,.true.,.true.)
#endif

#ifdef TRACERS_OCEAN
#ifndef TRACERS_SPECIAL_O18
C**** main ocean variabiles and gradients
          do j=J_0,J_1
          do i=1,im
            do l=1,lmm(i,j)
              trmo(i,j,l,n)=trw0(n)*(mo(i,j,l)*dxypo(j)-s0m(i,j,l))
              txmo(i,j,l,n)=-trw0(n)*sxmo(i,j,l)
              tymo(i,j,l,n)=-trw0(n)*symo(i,j,l)
              tzmo(i,j,l,n)=-trw0(n)*szmo(i,j,l)
            end do
          end do
          end do
#else
C**** read in initial conditions for isotopes
C**** search through for correct tracer (since there is no guarantee
C**** that they will be in same order as tracer indices).
C**** data are now in 'per mil'
          rewind (iu_O18ic)
 10       read  (iu_O18ic,err=800,end=810) title,t0m4,tzm4
          if (index(title,trim(trname(n))).eq.0) goto 10
          write (6,*) 'Read from H2O18ic: ',title
          call closeunit(iu_O18ic)

C**** Turn per mil data into mass ratios (using current standard)
          t0m4(:,:,:)=(t0m4(:,:,:)*1d-3+1.)*trw0(n)
          tzm4(:,:,:)=0.    ! tzm4(:,:,:)*1d-3*trw0(n) corrupted?
C****
          do l=1,lmo
            txmo(:,J_0:J_1,l,n) = 0.
            tymo(:,J_0:J_1,l,n) = 0.
C**** Define East-West horizontal gradients
            im1=im-1
            i=im
            do j=J_0S,J_1S
              do ip1=1,im
                if (lmm(i,j).ge.l) then
                  if (lmm(im1,j).ge.l) then
                    if (lmm(ip1,j).ge.l) then
                      txmo(i,j,l,n)=2.5d-1*(t0m4(ip1,j,l)-t0m4(im1,j,l))
                    else
                      txmo(i,j,l,n)=  5d-1*(t0m4(  i,j,l)-t0m4(im1,j,l))
                    end if
                  else
                    if (lmm(ip1,j).ge.l)
     *                   txmo(i,j,l,n)=5d-1*(t0m4(ip1,j,l)-t0m4(i,j,l))
                  end if
                end if
                im1=i
                i=ip1
              end do
            end do
C**** Define North-South horizontal gradients
            do j=J_0S,J_1S
              do i=1,im
                if (lmm(i,j).ge.l)  then
                  if (lmm(i,j-1).ge.l)  then
                    if (lmm(i,j+1).ge.l)  then
                      tymo(i,j,l,n)=2.5d-1*(t0m4(i,j+1,l)-t0m4(i,j-1,l))
                    else
                      tymo(i,j,l,n)=  5d-1*(t0m4(i,  j,l)-t0m4(i,j-1,l))
                    end if
                  else
                    if (lmm(i,j+1).ge.l)
     *                   tymo(i,j,l,n)=5d-1*(t0m4(i,j+1,l)-t0m4(i,j,l))
                  end if
                end if
              end do
            end do
C**** Multiply ratios by freshwater mass
            do j=J_0,J_1
            do i=1,im
              if (l.le.lmm(i,j)) then
                fwm = mo(i,j,l)*dxypo(j)-s0m(i,j,l)
                trmo(i,j,l,n)=t0m4(i,j,l)*fwm
                txmo(i,j,l,n)=txmo(i,j,l,n)*fwm-sxmo(i,j,l)*t0m4(i,j,l)
                tymo(i,j,l,n)=tymo(i,j,l,n)*fwm-symo(i,j,l)*t0m4(i,j,l)
                tzmo(i,j,l,n)=tzm4(i,j,l)  *fwm-szmo(i,j,l)*t0m4(i,j,l)
              else
                trmo(i,j,l,n)=0.
                txmo(i,j,l,n)=0.
                tymo(i,j,l,n)=0.
                tzmo(i,j,l,n)=0.
              end if
            end do
            end do
          end do

C**** Initiallise strait values based on adjacent ocean boxes
          call gather_ocean(1)  ! mo,g0m,gx-zmo,s0m,sx-zmo,trmo,tx-zmo

          if(am_I_root()) then
          do nst=1,nmst
#ifdef TRACERS_OCEAN
            i1=ist(nst,1)
            j1=jst(nst,1)
            i2=ist(nst,2)
            j2=jst(nst,2)
            do l=1,lmst(nst)
              t01=trmo_glob(i1,j1,l,n)/(mo_glob(i1,j1,l)*dxypo(j1)
     *             -s0m_glob(i1,j1,l))
              t02=trmo_glob(i2,j2,l,n)/(mo_glob(i2,j2,l)*dxypo(j2)
     *             -s0m_glob(i2,j2,l))
              trmst(l,nst,n) = 5d-1*(mmst(l,nst)-s0mst(l,nst))*(t01+t02)
              txmst(l,nst,n) = -5d-1*sxmst(l,nst)*(t01+t02)
              tzmst(l,nst,n) = -5d-1*szmst(l,nst)*(t01+t02)
            end do
            do l=lmst(nst)+1,lmo
              trmst(l,nst,n) = 0.
              txmst(l,nst,n) = 0.
              tzmst(l,nst,n) = 0.
            end do
#endif
#ifdef TRACERS_WATER
            trsist(n,1:2,nst) = trw0(n)*(msist(1,nst)*xsi(1:2)
     *           -ssist(1:2,nst))
            trsist(n,3:lmi,nst)=trw0(n)*(msist(2,nst)*xsi(3:lmi)
     *           -ssist(3:lmi,nst))
#endif
          end do

        end if

        call bcast_straits(.false.) ! bcst tracers

C**** Balance tracers so that average concentration is TRW0 
C**** or oc_tracer_mean
          
          CALL CONSERV_OTR(OTRACJ,N)
          OTRACJ(:)=OTRACJ(:)*DXYP(:)

          CALL GLOBALSUM(grid, OTRACJ, trsum, ALL=.true.)

          if (AM_I_ROOT()) then
            if (oc_tracer_mean(n).ne.-999.) then
              tratio=trw0(n)*(oc_tracer_mean(n)*1d-3+1.)
            else
              tratio=trw0(n)
            end if
            
            if (trname(n).eq.'Water') then
              wsum = trsum
              frac_tr = -999.
            else
              frac_tr = tratio/(trsum/wsum)
              write(6,*) "Average oceanic tracer concentration ",
     *             trname(n),(trsum/(wsum*trw0(n))-1d0)*1d3,frac_tr
            end if
          end if
          call ESMF_BCAST(grid,  frac_tr)
         
          if (frac_tr.ne.-999.) then    
            do l=1,lmo
              do j=J_0,J_1
                do i=1,imaxj(j)
                  trmo(i,j,l,n) = trmo(i,j,l,n) * frac_tr
                  if (frac_tr.lt.1.) then
                    TXMO(I,J,L,n)=frac_tr*TXMO(I,J,L,n)
                    TYMO(I,J,L,n)=frac_tr*TYMO(I,J,L,n)
                    TZMO(I,J,L,n)=frac_tr*TZMO(I,J,L,n)
                  end if
                end do
              end do
            end do
C**** straits
            if (am_i_root()) then
              DO NST=1,NMST
                DO L=1,LMST(NST)
                  TRMST(L,NST,N)=frac_tr*TRMST(L,NST,N)
                  if (frac_tr.lt.1) then
                    TXMST(L,NST,N)=frac_tr*TXMST(L,NST,N)
                    TZMST(L,NST,N)=frac_tr*TZMST(L,NST,N)
                  end if
                END DO
              END DO
            end if
            CALL ESMF_BCAST(grid, TRMST)
            CALL ESMF_BCAST(grid, TXMST)
            CALL ESMF_BCAST(grid, TZMST)

C**** Check
            CALL CONSERV_OTR(OTRACJ,N)
            OTRACJ(:)=OTRACJ(:)*DXYP(:)
            CALL GLOBALSUM(grid, OTRACJ, trsum, ALL=.true.)

            if (AM_I_ROOT()) then
              tratio=(trsum/(wsum*trw0(n))-1.)*1000.
              write(6,*) "New ocean tracer mean: ",tratio
     *             ,oc_tracer_mean(n)
            end if
          end if
#endif
#ifdef TRACERS_WATER
          do j=J_0,J_1
            do i=1,im
              if (focean(i,j).gt.0) gtracer(n,1,i,j)=trmo(i,j,1,n)/
     *              (mo(i,j,1)*dxypo(j)-s0m(i,j,1))
            end do
          end do
#endif
#endif
C****
        end select
        write(6,*) trname(n)," tracer initialised in ocean"
        end if
      end do

      return
 800  write(6,*) "Error reading input file H2O18ic"
      call stop_model('stopped in OCN_TRACER.f',255)
 810  write(6,*) "Tracer ",trname(n)," not found in file H2O18ic"
      call stop_model('stopped in OCN_TRACER.f',255)
      end subroutine tracer_ic_ocean
#endif

#ifdef TRACERS_OCEAN
      subroutine init_tracer_ocean
!@sum Initialise ocean tracer diagnostics (none so far)
      return
      end subroutine init_tracer_ocean
#endif

#ifdef TRACERS_OCEAN
      SUBROUTINE OC_TDECAY
!@sum OC_TDECAY decays radioactive tracers in ocean
!@auth Gavin Schmidt/Jean Lerner
      USE MODEL_COM, only : itime,dtsrc
      USE TRACER_COM, only : ntm,trdecay,itime_tr0
      USE OCEAN, only : trmo,txmo,tymo,tzmo
      IMPLICIT NONE
      real*8, save, dimension(ntm) :: expdec = 1.
      logical, save :: ifirst=.true.
      integer n

      if (ifirst) then               
        do n=1,ntm
          if (trdecay(n).gt.0.0) expdec(n)=exp(-trdecay(n)*dtsrc)
        end do
        ifirst = .false.
      end if

      do n=1,ntm
        if (trdecay(n).gt.0. .and. itime.ge.itime_tr0(n)) then
C**** Oceanic decay
          trmo(:,:,:,n)   = expdec(n)*trmo(:,:,:,n)
          txmo(:,:,:,n)   = expdec(n)*txmo(:,:,:,n)
          tymo(:,:,:,n)   = expdec(n)*tymo(:,:,:,n)
          tzmo(:,:,:,n)   = expdec(n)*tzmo(:,:,:,n)
        end if
      end do
C****
      return
      end subroutine oc_tdecay
#endif

      SUBROUTINE conserv_OTR(OTR,ITR)
!@sum  conserv_OTR calculates zonal ocean tracer on atmos grid
!@auth Gavin Schmidt
!@ver  1.0
      USE GEOM, only : bydxyp
      USE OCEAN, only : im,jm,fim,imaxj,focean,mo,lmm,trmo
      USE STRAITS, only : nmst,jst,trmst
      USE DOMAIN_DECOMP, only : GET, GRID
      IMPLICIT NONE
!@var OSALT zonal ocean tracer (kg/m^2)
      REAL*8, DIMENSION(GRID%J_STRT:GRID%J_STOP) :: OTR
      INTEGER I,J,L,N,ITR

      INTEGER :: J_0, J_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      CALL GET(grid, J_STRT=J_0,    J_STOP=J_1,
     &  HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,HAVE_NORTH_POLE=HAVE_NORTH_POLE)

      DO J=J_0,J_1
        OTR(J)=0
        DO I=1,IMAXJ(J)
          DO L=1,LMM(I,J)
            OTR(J) = OTR(J) + TRMO(I,J,L,ITR)*FOCEAN(I,J)*BYDXYP(J)
          END DO
        END DO
      END DO
      if (HAVE_SOUTH_POLE) OTR(1) =FIM*OTR(1)
      if (HAVE_NORTH_POLE) OTR(JM)=FIM*OTR(JM)
C**** include straits variables
      DO N=1,NMST
        J=JST(N,1)
        if(j.lt.j_0 .or. j.gt.j_1) cycle
        OTR(J)=OTR(J)+SUM(TRMST(:,N,ITR))*BYDXYP(J)
      END DO
C****
      RETURN
      END SUBROUTINE conserv_OTR
