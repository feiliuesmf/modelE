#include "rundeck_opts.h"

#ifdef TRACERS_WATER
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
      USE MODEL_COM, only: itime
      USE OCEAN, only : im,jm,lmo,focean,dxypo,mo,lmm
#ifdef TRACERS_OCEAN
     *     ,trmo,txmo,tymo,tzmo,s0m,sxmo,symo,szmo
#endif
      USE SEAICE, only : xsi,lmi
      USE STRAITS, only : nmst,msist,ssist,trsist
#ifdef TRACERS_OCEAN
     *     ,lmst,ist,jst,xst,yst,mmst,s0mst,sxmst,szmst,trmst,txmst
     *     ,tzmst
#endif
      USE TRACER_COM, only : itime_tr0, ntm, trname, trw0
      USE FLUXES, only : gtracer
      IMPLICIT NONE
      integer n,i,j,l,k,nst,i1,j1,i2,j2
      real*8 t01,t02

C**** Note that only sea ice related arrays are initialised if
C**** only TRACERS_WATER is true.

      do n=1,ntm

        if (itime.eq.itime_tr0(n)) then
        select case (trname(n)) 
          
        case default
c          write(6,*) 'In TRACER_IC:',trname(n),' does not exist '
          stop "TRACER_IC_OCEAN"
          
        case ('Water')
          
C**** main ocean variabiles and gradients
          do j=1,jm
          do i=1,im
            if (focean(i,j).gt.0) then
#ifdef TRACERS_OCEAN
              do l=1,lmm(i,j)
                trmo(i,j,l,n)=trw0(n)*(mo(i,j,l)*dxypo(j)-s0m(i,j,l))
                txmo(i,j,l,n)=-trw0(n)*sxmo(i,j,l)
                tymo(i,j,l,n)=-trw0(n)*symo(i,j,l)
                tzmo(i,j,l,n)=-trw0(n)*szmo(i,j,l)
              end do
#endif
              gtracer(n,1,i,j)=trw0(n)
            end if
#ifdef TRACERS_OCEAN
            do l=lmm(i,j)+1,lmo
              trmo(i,j,l,n)=0. ; txmo(i,j,l,n)=0.
              tymo(i,j,l,n)=0. ; tzmo(i,j,l,n)=0.
            end do
#endif
          end do
          end do
           
C**** Initiallise strait values based on adjacent ocean boxes
          do nst=1,nmst
#ifdef TRACERS_OCEAN
            i1=ist(nst,1)
            j1=jst(nst,1)
            i2=ist(nst,2)
            j2=jst(nst,2)
            do l=1,lmst(nst)
              t01=(trmo(i1,j1,l,n)+xst(nst,1)*txmo(i1,j1,l,n)
     *             +yst(nst,1)*tymo(i1,j1,l,n))/(mo(i1,j1,l)
     *             *dxypo(j1)-s0m(i1,j1,l))
              t02=(trmo(i2,j2,l,n)+xst(nst,2)*txmo(i2,j2,l,n)
     *             +yst(nst,2)*tymo(i2,j2,l,n))/(mo(i2,j2,l)
     *             *dxypo(j2)-s0m(i2,j2,l))
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
            trsist(n,1:2,nst) = trw0(n)*(msist(1,nst)*xsi(1:2)
     *           -ssist(1:2,nst))
            trsist(n,3:lmi,nst)=trw0(n)*(msist(2,nst)*xsi(3:lmi)
     *           -ssist(3:lmi,nst))
          end do
C**** 
        end select
        write(6,*) trname(n)," tracer initialised in ocean"
        end if
      end do

      return
      end subroutine tracer_ic_ocean
#endif

#ifdef TRACERS_OCEAN
      subroutine init_tracer_ocean
!@sum Initialise ocean tracer diagnostics (none so far)
      return
      end subroutine init_tracer_ocean
#endif
