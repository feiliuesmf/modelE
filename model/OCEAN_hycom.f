#include "rundeck_opts.h"
      SUBROUTINE init_OCEAN(iniOCEAN,istart)
      USE MODEL_COM, only : im,jm,focean
      USE FLUXES, only : gtemp
#ifdef TRACERS_GASEXCH_Natassa
     .                  ,GTRACER,TRGASEX

      USE TRACER_COM, only : ntm    !tracers involved in air-sea gas exch

      USE TRACER_GASEXCH_COM, only : atrac
#endif
      integer istart
#include "dimensions.h"
#include "dimension2.h"
#include "cpl.h"
css   if (istart.eq.2 .or. nstep0.eq.0) call geopar
      call inicon
c
c --- increase temp by 2 deg
c     do 21 j=1,jj
c     do 21 l=1,isp(j)
c     do 21 i=ifp(j,l),ilp(j,l)
c       if (latij(i,j,3).lt.-65..and.lonij(i,j,3).le.5.) then !lat[-90:90],lon[0:360]
c         p(i,j,1)=0.
c         do k=1,kk
c         p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
c         temp(i,j,   k)=temp(i,j,   k)+2.*(min(p(i,j,k+1),200.*onem)-
c    .    min(p(i,j,k),200.*onem))/max(dp(i,j,k),onemm)
c         temp(i,j,kk+k)=temp(i,j,kk+k)+2.*(min(p(i,j,k+1),200.*onem)-
c    .    min(p(i,j,k),200.*onem))/max(dp(i,j,k),onemm)
c         enddo
c       endif
c21   continue
c
      if (istart.gt.2) then               !istart=2 has done this in inirfn
      DO J=1,JM
      DO I=1,IM
        IF (FOCEAN(I,J).gt.0.) THEN
          GTEMP(1,1,I,J)=asst(I,J)
#ifdef TRACERS_GASEXCH_Natassa
        do nt=1,ntm
        GTRACER(nt,1,I,J)=atrac(I,J,nt)
        enddo
#endif
        END IF
      END DO
      END DO
      endif
c
      END SUBROUTINE init_OCEAN
c
      SUBROUTINE DUMMY_OCN
!@sum  DUMMY necessary entry points for non-dynamic/non-deep oceans
!@auth Gavin Schmidt
!@ver  1.0
css   ENTRY ODYNAM
      ENTRY ODIFS
      ENTRY io_ocdiag
      ENTRY init_ODEEP
      ENTRY reset_ODIAG
      ENTRY diag_OCEAN
      entry OSTRUC(QTCHNG)
      entry OCLIM(end_of_day)
      entry OSOURC (ROICE,SMSI,TGW,WTRO,OTDT,RUN0,F0DT,F2DT,RVRRUN
     *           ,RVRERUN,EVAPO,EVAPI,TFW,RUN4O,ERUN4O,RUN4I,ERUN4I
     *           ,ENRGFO,ACEFO,ACE2F,ENRGFI)
css   entry daily_OCEAN(end_of_day)
      entry PRECIP_OC
      entry GROUND_OC
      entry DIAGCO (M)
      entry io_oda(kunit,it,iaction,ioerr)
css   entry io_ocean(iu_GIC,ioread,ioerr)
css   entry CHECKO(SUBR)
c
      ENTRY ADVSI_DIAG
      entry alloc_ocean
c --- not calling ice dynamics
css      ENTRY DYNSI
css      ENTRY ADVSI
css      ENTRY io_icedyn
css      ENTRY io_icdiag
css      ENTRY init_icedyn
css      ENTRY reset_icdiag
css      ENTRY diag_ICEDYN
c
      RETURN
      END SUBROUTINE DUMMY_OCN
c
      SUBROUTINE io_ocean(kunit,iaction,ioerr)
!@sum  io_ocean outputs ocean related fields for restart
!@ver  1.0       
      USE MODEL_COM, only : ioread,iowrite,irsficno,irsfic
     *     ,irsficnt,irerun,lhead
      USE FLUXES, only : sss,ogeoza,uosurf,vosurf,dmsi,dhsi,dssi
#ifdef TRACERS_GASEXCH_Natassa
      USE TRACER_GASEXCH_COM, only : atrac
#endif
      IMPLICIT NONE
#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"
#include "cpl.h"
c
      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "OCDYN01"
#ifdef TRACERS_OCEAN
!@var TRHEADER Character string label for individual records
      CHARACTER*80 :: TRHEADER, TRMODULE_HEADER = "TROCDYN02"
c
      write (TRMODULE_HEADER(lhead+1:80),'(a13,i3,a1,i3,a)')
     *     'R8 dim(im,jm,',LMO,',',NTM,'):TRMO,TX,TY,TZ'
#endif
c
css   write (MODULE_HEADER(lhead+1:80),'(a13,i2,a)') 'R8 dim(im,jm,',
css  *   LMO,'):M,U,V,G0,GX,GY,GZ,S0,SX,SY,SZ, OGZ,OGZSV'
c
      write(*,'(a,i9,f9.0)')'chk ocean write at nstep/day=',nstep,time
      write (MODULE_HEADER(lhead+1:80),'(a,i8,f8.1,a)')
     . 'u,v,dp,t,s,th,tb,ub,vb,pb,pb,psi,thk,mxl,uf,vf,df,tcr3+o18+a8'
      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
css     WRITE (kunit,err=10) MODULE_HEADER,MO,UO,VO,G0M,GXMO,GYMO,GZMO
css  *     ,S0M,SXMO,SYMO,SZMO,OGEOZ,OGEOZ_SV
css#ifdef TRACERS_OCEAN
css       WRITE (kunit,err=10) TRMODULE_HEADER,tracer
css#endif
        WRITE (kunit,err=10) MODULE_HEADER,nstep,time
     . ,u,v,dp,temp,saln,th3d,ubavg,vbavg,pbavg,pbot,psikk,thkk
     . ,dpmixl,uflxav,vflxav,diaflx,tracer,dpinit,oddev
     . ,uav,vav,dpuav,dpvav,dpav,temav,salav,th3av,ubavav,vbavav
     . ,pbavav,sfhtav,eminpav,surflav,sflxav,brineav,dpmxav,oiceav
     . ,asst,sss,ogeoza,uosurf,vosurf,dhsi,dmsi,dssi         ! agcm grid
     . ,scpx,scux,scvx,scqx,scpy,scuy,scvy,scqy,scp2,scu2,scv2,scq2
#ifdef TRACERS_GASEXCH_Natassa
     . ,atrac
#endif
      CASE (IOREAD:)            ! input from restart file
        SELECT CASE (IACTION)
          CASE (IRSFICNO)   ! initial conditions (no ocean data)
            READ (kunit)
          CASE (ioread,irerun,irsfic) ! restarts
css         READ (kunit,err=10) HEADER,MO,UO,VO,G0M,GXMO,GYMO,GZMO,S0M
css  *           ,SXMO,SYMO,SZMO,OGEOZ,OGEOZ_SV
c
            print *,' calling geopar'
            call geopar
            READ (kunit,err=10) HEADER,nstep0,time0
     . ,u,v,dp,temp,saln,th3d,ubavg,vbavg,pbavg,pbot,psikk,thkk
     . ,dpmixl,uflxav,vflxav,diaflx,tracer,dpinit,oddev
     . ,uav,vav,dpuav,dpvav,dpav,temav,salav,th3av,ubavav,vbavav
     . ,pbavav,sfhtav,eminpav,surflav,sflxav,brineav,dpmxav,oiceav
     . ,asst,sss,ogeoza,uosurf,vosurf,dhsi,dmsi,dssi         ! agcm grid
     . ,scpx,scux,scvx,scqx,scpy,scuy,scvy,scqy,scp2,scu2,scv2,scq2
#ifdef TRACERS_GASEXCH_Natassa
     . ,atrac
#endif
      nstep0=time0*86400./baclin+.0001
      write(*,'(a,i9,f9.0)')'chk ocean read at nstep/day=',nstep0,time0
      nstep=nstep0
      time=time0
c
            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",HEADER
     *             ,MODULE_HEADER
              GO TO 10
            END IF
#ifdef TRACERS_OCEAN
            READ (kunit,err=10) TRHEADER,TRMO,TXMO,TYMO,TZMO
            IF (TRHEADER(1:LHEAD).NE.TRMODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",TRHEADER
     *             ,TRMODULE_HEADER
              GO TO 10
            END IF
#endif
          CASE (irsficnt) ! restarts (never any tracer data)
css         READ (kunit,err=10) HEADER,MO,UO,VO,G0M,GXMO,GYMO,GZMO,S0M
css  *           ,SXMO,SYMO,SZMO,OGEOZ,OGEOZ_SV
c
            call geopar
            READ (kunit,err=10) HEADER,nstep0,time0
     . ,u,v,dp,temp,saln,th3d,ubavg,vbavg,pbavg,pbot,psikk,thkk
     . ,dpmixl,uflxav,vflxav,diaflx,tracer,dpinit,oddev
     . ,uav,vav,dpuav,dpvav,dpav,temav,salav,th3av,ubavav,vbavav
     . ,pbavav,sfhtav,eminpav,surflav,sflxav,brineav,dpmxav,oiceav
     . ,asst,sss,ogeoza,uosurf,vosurf,dhsi,dmsi,dssi         ! agcm grid
     . ,scpx,scux,scvx,scqx,scpy,scuy,scvy,scqy,scp2,scu2,scv2,scq2
      nstep0=time0*86400./baclin+.0001
      write(*,'(a,i9,f9.0)')'chk ocean read at nstep/day=',nstep0,time0
            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",HEADER
     *             ,MODULE_HEADER
              GO TO 10
            END IF
          END SELECT
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
C****
      END SUBROUTINE io_ocean
c
      SUBROUTINE CHECKO(SUBR)
!@sum  CHECKO Checks whether Ocean are reasonable
!@ver  1.0
      USE MODEL_COM, only : im,jm
      USE FLUXES, only : gtemp
      USE MODEL_COM, only : focean

      IMPLICIT NONE
      integer i,j

!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

      print *,'SUBR=',SUBR
      write(*,'(10f7.2)') ((gtemp(1,1,i,j),i=1,10),j=15,20)
      write(*,'(a)') 'focean'
      write(*,'(10f7.2)') ((focean(i,j),i=1,10),j=15,20)
      END SUBROUTINE CHECKO
c

      SUBROUTINE daily_OCEAN
!@sum  daily_OCEAN performs the daily tasks for the ocean module
!@auth Original Development Team
!@ver  1.0
C****
      RETURN
      END SUBROUTINE daily_OCEAN
c
css   REAL*8 FUNCTION TFREZS (SIN)
C****
C**** TFREZS calculates the freezing temperature of sea water as a
C**** function of salinity.
C**** The reference for this function is:
C**** N.P. Fofonoff and R.C. Millard Jr., 1983.  Algorithms for
C**** Computation of Fundamental Properties of Seawater.  UNESCO
C**** Technical Papers in Marine Science, volume 44.
C****
C**** Input: SIN (1) = salinity (kg NaCl/kg sea water), from .004 to .04
C****
C**** Output: TFREZS (C) = freezing temperature of sea water
C****
css   IMPLICIT NONE
css    REAL*8, INTENT(IN) :: SIN
css      REAL*8 :: A01 = -.0575d0, A02 = -2.154996D-4, A03 =1.710523D-3
css      REAL*8 S,S32
C****
css      S   = SIN*1.D3
css      S32 = S*DSQRT(S)
css      TFREZS = (A01 + A02*S)*S + A03*S32
css      RETURN
css      END
c
      subroutine gather_odiags
C     nothing to gather - ocean prescribed
      return
      end subroutine gather_odiags
