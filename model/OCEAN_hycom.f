#include "rundeck_opts.h"
      SUBROUTINE init_OCEAN(iniOCEAN,istart)
      USE DOMAIN_DECOMP, only: AM_I_ROOT,ESMF_BCAST
      USE HYCOM_ATM, only : gather_atm,scatter_atm, focean,gtemp,gtempr,
     &     asst,atempr,im,jm
#ifdef TRACERS_GASEXCH_Natassa
     .    ,GTRACER
#endif
!!      USE MODEL_COM, only : im,jm,focean
!!      USE FLUXES, only : gtemp
#ifdef TRACERS_GASEXCH_Natassa
      USE FLUXES, only : TRGASEX !,GTRACER

      USE TRACER_COM, only : ntm    !tracers involved in air-sea gas exch

      USE TRACER_GASEXCH_COM, only : atrac
#endif

#ifdef TRACERS_OceanBiology
      USE obio_forc, only : avgq,tirrq3d,ihra
      USE obio_com,  only : gcmax
#endif
      USE HYCOM_DIM
      USE HYCOM_SCALARS, only : delt1, salmin
     &    , nstep0, nstep, time0, time
      USE HYCOM_ARRAYS_GLOB, only: scatter_hycom_arrays
      implicit none

      logical, intent(in) :: iniOCEAN
      integer, intent(in) :: istart
!!#include "dimensions.h"
!!#include "dimension2.h"
!!#include "cpl.h"
      integer i,j

      ! move to global atm grid
      call gather_atm

             call geopar(iniOCEAN)

      if (AM_I_ROOT()) then ! work on global grids here
      
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

      if (istart.gt.2) then               !istart=2 has done this inicon
      DO J=1,JM
      DO I=1,IM
        IF (FOCEAN(I,J).gt.0.) THEN
          GTEMP(1,1,I,J)=asst(I,J)
          gtempr(1,I,J)=atempr(I,J)
#ifdef TRACERS_GASEXCH_Natassa
        GTRACER(1:ntm,1,I,J)=atrac(I,J,1:ntm)
#endif
        END IF
      END DO
      END DO
      endif

!!! I guess I had a good reason for commenting this out... 
!!! (probably should done in inicon) IA

!!!! the following is already done in inicon
!!      if (istart.gt.2) then               !istart=2 has done this in inirfn
!!      DO J=1,JM
!!      DO I=1,IM
!!        IF (FOCEAN(I,J).gt.0.) THEN
!!          GTEMP(1,1,I,J)=asst(I,J)
!!! GTEMPR ??
!!#ifdef TRACERS_GASEXCH_Natassa
!!        do nt=1,ntm
!!        GTRACER(nt,1,I,J)=atrac(I,J,nt)
!!        enddo
!!#endif
!!        END IF
!!      END DO
 !!     END DO
!!      endif

      endif ! AM_I_ROOT
      call scatter_atm
      call scatter_hycom_arrays

!!! hack needed for serial inicon
      CALL ESMF_BCAST(ogrid, delt1 )

      CALL ESMF_BCAST(ogrid, salmin )
      CALL ESMF_BCAST(ogrid, nstep0 )
      CALL ESMF_BCAST(ogrid, nstep )
      CALL ESMF_BCAST(ogrid, time0 )
      CALL ESMF_BCAST(ogrid, time )


c
      END SUBROUTINE init_OCEAN
c
      SUBROUTINE DUMMY_OCN
!@sum  DUMMY necessary entry points for non-dynamic/non-deep oceans
!@auth Gavin Schmidt
!@ver  1.0
css   ENTRY ODYNAM
      !! fix later: implicit none

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
!!      entry alloc_ocean
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
      USE DOMAIN_DECOMP, only: AM_I_ROOT, pack_data, unpack_data,
     &     ESMF_BCAST
      USE HYCOM_ATM, only : gather_atm,scatter_atm,
     &     sss,ogeoza,uosurf,vosurf,dmsi,dhsi,dssi,asst,atempr
      USE MODEL_COM, only : ioread,iowrite,irsficno,irsfic
     *     ,irsficnt,irerun,lhead
!!      USE FLUXES, only : sss,ogeoza,uosurf,vosurf,dmsi,dhsi,dssi
      USE HYCOM_DIM_GLOB, only : kk,iia,jja,kdm,idm,jdm
      USE HYCOM_DIM, only : ogrid
      USE HYCOM_SCALARS, only : nstep,time,oddev,nstep0,time0,baclin
     &     ,onem
#ifdef TRACERS_GASEXCH_Natassa
      USE TRACER_GASEXCH_COM, only : atrac
#endif

#ifdef TRACERS_OceanBiology
      USE obio_forc, only : avgq,tirrq3d,ihra
      USE obio_com,  only : gcmax
#endif
      USE HYCOM_ARRAYS_GLOB
      IMPLICIT NONE
!!#include "dimensions.h"
!!#include "dimension2.h"
!!#include "common_blocks.h"
!!#include "cpl.h"
      !!! hack
      !!!real asst
c
      INTEGER, intent(in) :: kunit   !@var kunit unit number of read/write
      INTEGER, intent(in) :: iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "OCDYN01"
#if defined(TRACERS_GASEXCH_Natassa) || defined(TRACERS_OceanBiology)
      integer i,j,k
!@var TRNHEADER Character string label for individual records
      CHARACTER*80 :: TRNHEADER, TRNMODULE_HEADER = "TRGASEX-OBIO"
#endif
#ifdef TRACERS_OceanBiology
      real, allocatable :: avgq_glob(:,:,:),tirrq3d_glob(:,:,:),
     &     gcmax_glob(:,:,:)
      integer, allocatable :: ihra_glob(:,:)
#endif

#ifdef TRACERS_OCEAN
!@var TRHEADER Character string label for individual records
      CHARACTER*80 :: TRHEADER, TRMODULE_HEADER = "TROCDYN02"
c
      write (TRMODULE_HEADER(lhead+1:80),'(a13,i3,a1,i3,a)')
     *     'R8 dim(im,jm,',LMO,',',NTM,'):TRMO,TX,TY,TZ'
#endif

#ifdef TRACERS_OceanBiology
      if (AM_I_ROOT()) then
        allocate( avgq_glob(idm,jdm,kdm),tirrq3d_glob(idm,jdm,kdm),
     &       ihra_glob(idm,jdm), gcmax_glob(idm,jdm,kdm) )
      endif
      call pack_data(ogrid, avgq, avgq_glob)
      call pack_data(ogrid, tirrq3d, tirrq3d_glob)
      call pack_data(ogrid, ihra, ihra_glob)
      call pack_data(ogrid, gcmax, gcmax_glob)
#endif


      ! move to global atm grid
      call gather_atm
!!!!!!call scatter_tracer        ! caution: temporary fix
                                 ! to allow Romanou to run with
                                 ! global tracer
      call gather_hycom_arrays   !mkb Jun  6

      if (AM_I_ROOT()) then ! work on global grids here

c
css   write (MODULE_HEADER(lhead+1:80),'(a13,i2,a)') 'R8 dim(im,jm,',
css  *   LMO,'):M,U,V,G0,GX,GY,GZ,S0,SX,SY,SZ, OGZ,OGZSV'
c
      write(*,'(a,i9,f9.0)')'chk ocean write at nstep/day=',nstep,time
      write (MODULE_HEADER(lhead+1:80),'(a,i8,f8.1,a)')
     . 'u,v,dp,t,s,th,tb,ub,vb,pb,pb,psi,thk,mxl,uf,vf,df,tcr3+o18+a8'

#if defined(TRACERS_GASEXCH_Natassa) && defined(TRACERS_OceanBiology)
      write(*,'(a,i9,f9.0)')'chk GASEXCH write at nstep/day=',nstep,time
      write (TRNMODULE_HEADER(lhead+1:80),'(a29)')
     *     'atrac,avgq,gcmax,tirrq3d,ihra'
#else
#ifdef TRACERS_GASEXCH_Natassa
      write(*,'(a,i9,f9.0)')'chk GASEXCH write at nstep/day=',nstep,time
      write (TRNMODULE_HEADER(lhead+1:80),'(a5)')
     *     'atrac'
#endif
#ifdef TRACERS_OceanBiology
      write(*,'(a,i9,f9.0)')'chk GASEXCH write at nstep/day=',nstep,time
      write (TRNMODULE_HEADER(lhead+1:80),'(a23)')
     *     'avgq,gcmax,tirrq3d,ihra'
#endif
#endif

      SELECT CASE (IACTION)
c---------------------------------------------------------------------------------
      CASE (:IOWRITE)            ! output to standard restart file
css     WRITE (kunit,err=10) MODULE_HEADER,MO,UO,VO,G0M,GXMO,GYMO,GZMO
css  *     ,S0M,SXMO,SYMO,SZMO,OGEOZ,OGEOZ_SV
css#ifdef TRACERS_OCEAN
css       WRITE (kunit,err=10) TRMODULE_HEADER,tracer
css#endif
        WRITE (kunit,err=10) MODULE_HEADER,nstep,time
     . ,u,v,dp,temp,saln,th3d,thermb,ubavg,vbavg,pbavg,pbot,psikk,thkk
     . ,dpmixl,uflxav,vflxav,diaflx,tracer,dpinit,oddev
     . ,uav,vav,dpuav,dpvav,dpav,temav,salav,th3av,ubavav,vbavav
     . ,pbavav,sfhtav,eminpav,surflav,sflxav,brineav,dpmxav,oiceav
     . ,asst,atempr,sss,ogeoza,uosurf,vosurf,dhsi,dmsi,dssi         ! agcm grid

#if defined(TRACERS_GASEXCH_Natassa) && defined(TRACERS_OceanBiology)
      WRITE (kunit,err=10) TRNMODULE_HEADER,nstep,time
     . ,atrac,avgq_glob,gcmax_glob,tirrq3d_glob,ihra_glob
      i=100
      j=100
      do k=1,kdm
      write(*,'(a,i2,7(e12.4,1x),i3)') ' tst1a k=',k,
     .    dp(i,j,k)/onem,temp(i,j,k),avgq_glob(i,j,k),gcmax_glob(i,j,k),
     .    tracer(i,j,k,1),tracer(i,j,k,15),tirrq3d_glob(i,j,k),
     &       ihra_glob(i,j)
      enddo
#else
#ifdef TRACERS_GASEXCH_Natassa
      WRITE (kunit,err=10) TRNMODULE_HEADER,nstep,time
     . ,atrac
      i=100
      j=100
      do k=1,kdm
      write(*,'(a,i2,7(e12.4,1x),i3)') ' tst1a k=',k,
     .    dp(i,j,k)/onem,temp(i,j,k),avgq_glob(i,j,k),gcmax_glob(i,j,k),
     .    tracer(i,j,k,1),tracer(i,j,k,15),tirrq3d_glob(i,j,k),
     &       ihra_glob(i,j)
      enddo
#endif
#ifdef TRACERS_OceanBiology
      WRITE (kunit,err=10) TRNMODULE_HEADER,nstep,time
     . ,avgq_glob,gcmax_glob,tirrq3d_glob,ihra_glob
      i=100
      j=100
      do k=1,kdm
      write(*,'(a,i2,7(e12.4,1x),i3)') ' tst1a k=',k,
     .    dp(i,j,k)/onem,temp(i,j,k),avgq_glob(i,j,k),gcmax_glob(i,j,k),
     .    tracer(i,j,k,1),tracer(i,j,k,15),tirrq3d_glob(i,j,k),
     &       ihra_glob(i,j)
      enddo
#endif
#endif

c---------------------------------------------------------------------------------
      CASE (IOREAD:)            ! input from restart file
        SELECT CASE (IACTION)
c       --------------------------------------------------------------------------
          CASE (IRSFICNO)   ! initial conditions (no ocean data)
            READ (kunit)
c       --------------------------------------------------------------------------
          CASE (ioread,irerun,irsfic) ! restarts
css         READ (kunit,err=10) HEADER,MO,UO,VO,G0M,GXMO,GYMO,GZMO,S0M
css  *           ,SXMO,SYMO,SZMO,OGEOZ,OGEOZ_SV
c
            !!call geopar
            READ (kunit,err=10) HEADER,nstep0,time0
     . ,u,v,dp,temp,saln,th3d,thermb,ubavg,vbavg,pbavg,pbot,psikk,thkk
     . ,dpmixl,uflxav,vflxav,diaflx,tracer,dpinit,oddev
     . ,uav,vav,dpuav,dpvav,dpav,temav,salav,th3av,ubavav,vbavav
     . ,pbavav,sfhtav,eminpav,surflav,sflxav,brineav,dpmxav,oiceav
     . ,asst,atempr,sss,ogeoza,uosurf,vosurf,dhsi,dmsi,dssi         ! agcm grid

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

#if defined(TRACERS_GASEXCH_Natassa) && defined(TRACERS_OceanBiology)
      READ (kunit,err=10) TRNHEADER,nstep0,time0
     . ,atrac,avgq_glob,gcmax_glob,tirrq3d_glob,ihra_glob
      write(*,'(a,i9,f9.0)')'chk GASEXCH read at nstep/day=',nstep0,time0
      i=100
      j=100
      do k=1,kdm
      write(*,'(a,i2,7(e12.4,1x),i3)') ' tst1b k=',k,
     .    dp(i,j,k)/onem,temp(i,j,k),avgq_glob(i,j,k),gcmax_glob(i,j,k),
     .    tracer(i,j,k,1),tracer(i,j,k,15),tirrq3d_glob(i,j,k),
     &       ihra_glob(i,j)
      enddo
            IF (TRNHEADER(1:LHEAD).NE.TRNMODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",TRNHEADER
     .             ,TRNMODULE_HEADER
              GO TO 10
            END IF
#else
#ifdef TRACERS_GASEXCH_Natassa
      READ (kunit,err=10) TRNHEADER,nstep0,time0
     . ,atrac
      i=100
      j=100
      do k=1,kdm
      write(*,'(a,i2,7(e12.4,1x),i3)') ' tst1b k=',k,
     .    dp(i,j,k)/onem,temp(i,j,k),avgq_glob(i,j,k),gcmax_glob(i,j,k),
     .    tracer(i,j,k,1),tracer(i,j,k,15),tirrq3d_glob(i,j,k),
     &       ihra_glob(i,j)
      enddo
            IF (TRNHEADER(1:LHEAD).NE.TRNMODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",TRNHEADER
     .             ,TRNMODULE_HEADER
              GO TO 10
            END IF
#endif
#ifdef TRACERS_OceanBiology
      READ (kunit,err=10) TRNHEADER,nstep0,time0
     . ,avgq_glob,gcmax_glob,tirrq3d_glob,ihra_glob
      i=100
      j=100
      do k=1,kdm
      write(*,'(a,i2,7(e12.4,1x),i3)') ' tst1b k=',k,
     .    dp(i,j,k)/onem,temp(i,j,k),avgq_glob(i,j,k),gcmax_glob(i,j,k),
     .    tracer(i,j,k,1),tracer(i,j,k,15),tirrq3d_glob(i,j,k),
     &       ihra_glob(i,j)
      enddo
            IF (TRNHEADER(1:LHEAD).NE.TRNMODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",TRNHEADER
     .             ,TRNMODULE_HEADER
              GO TO 10
            END IF
#endif
#endif

#ifdef TRACERS_OCEAN
            READ (kunit,err=10) TRHEADER,TRMO,TXMO,TYMO,TZMO
            IF (TRHEADER(1:LHEAD).NE.TRMODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",TRHEADER
     *             ,TRMODULE_HEADER
              GO TO 10
            END IF
#endif

c       --------------------------------------------------------------------------
          CASE (irsficnt) ! restarts (never any tracer data)
css         READ (kunit,err=10) HEADER,MO,UO,VO,G0M,GXMO,GYMO,GZMO,S0M
css  *           ,SXMO,SYMO,SZMO,OGEOZ,OGEOZ_SV
c
            print*,'restarts (never any tracer data -irsficnt)'
            !!call geopar
            READ (kunit,err=10) HEADER,nstep0,time0
     . ,u,v,dp,temp,saln,th3d,thermb,ubavg,vbavg,pbavg,pbot,psikk,thkk
     . ,dpmixl,uflxav,vflxav,diaflx,tracer(:,:,:,1),dpinit,oddev
     . ,uav,vav,dpuav,dpvav,dpav,temav,salav,th3av,ubavav,vbavav
     . ,pbavav,sfhtav,eminpav,surflav,sflxav,brineav,dpmxav,oiceav
     . ,asst,atempr,sss,ogeoza,uosurf,vosurf,dhsi,dmsi,dssi         ! agcm grid

      nstep0=time0*86400./baclin+.0001
      write(*,'(a,i9,f9.0)')'chk ocean read at nstep/day=',nstep0,time0

            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",HEADER
     *             ,MODULE_HEADER
              GO TO 10
            END IF

      go to 222
#if defined(TRACERS_GASEXCH_Natassa) && defined(TRACERS_OceanBiology)
      READ (kunit,err=10) TRNHEADER,nstep0,time0
     . ,atrac,avgq_glob,gcmax_glob,tirrq3d_glob,ihra_glob
      write(*,'(a,i9,f9.0)')
     &     'chk GASEXCH read at nstep/day=',nstep0,time0
      i=100
      j=100
      do k=1,kdm
      write(*,'(a,i2,7(e12.4,1x),i3)') ' tst2 k=',k,
     .    dp(i,j,k)/onem,temp(i,j,k),avgq_glob(i,j,k),gcmax_glob(i,j,k),
     .    tracer(i,j,k,1),tracer(i,j,k,15),tirrq3d_glob(i,j,k),
     &       ihra_glob(i,j)
      enddo
      write(*,*)'atrac at (36,23) =',atrac(36,23,1)
            IF (TRNHEADER(1:LHEAD).NE.TRNMODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",TRNHEADER
     .             ,TRNMODULE_HEADER
              GO TO 10
            END IF
#else
#ifdef TRACERS_GASEXCH_Natassa
      READ (kunit,err=10) TRNHEADER,nstep0,time0
     . ,atrac
      i=100
      j=100
      do k=1,kdm
      write(*,'(a,i2,7(e12.4,1x),i3)') ' tst2 k=',k,
     .    dp(i,j,k)/onem,temp(i,j,k),avgq_glob(i,j,k),gcmax_glob(i,j,k),
     .    tracer(i,j,k,1),tracer(i,j,k,15),tirrq3d_glob(i,j,k),
     &       ihra_glob(i,j)
      enddo
            IF (TRNHEADER(1:LHEAD).NE.TRNMODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",TRNHEADER
     .             ,TRNMODULE_HEADER
              GO TO 10
            END IF
#endif
#ifdef TRACERS_OceanBiology
      READ (kunit,err=10) TRNHEADER,nstep0,time0
     . ,avgq_glob,gcmax_glob,tirrq3d_glob,ihra_glob
      i=100
      j=100
      do k=1,kdm
      write(*,'(a,i2,7(e12.4,1x),i3)') ' tst2 k=',k,
     .    dp(i,j,k)/onem,temp(i,j,k),avgq_glob(i,j,k),gcmax_glob(i,j,k),
     .    tracer(i,j,k,1),tracer(i,j,k,15),tirrq3d_glob(i,j,k),
     &       ihra_glob(i,j)
      enddo
            IF (TRNHEADER(1:LHEAD).NE.TRNMODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",TRNHEADER
     .             ,TRNMODULE_HEADER
              GO TO 10
            END IF
#endif
#endif

 222  continue

          END SELECT
      END SELECT

      endif ! AM_I_ROOT
      call scatter_atm
      CALL ESMF_BCAST(ogrid, nstep0 )
      CALL ESMF_BCAST(ogrid, time0 )

#ifdef TRACERS_OceanBiology
      call unpack_data(ogrid, avgq_glob, avgq)
      call unpack_data(ogrid, tirrq3d_glob, tirrq3d)
      call unpack_data(ogrid, ihra_glob, ihra)
      call unpack_data(ogrid, gcmax_glob, gcmax)
      if (AM_I_ROOT()) then
        deallocate( avgq_glob,tirrq3d_glob,
     &       ihra_glob, gcmax_glob )
      endif
#endif

      RETURN
 10   IOERR=1
      call scatter_atm
      ! why do we need return after error?
      call stop_model("error i/o in io_ocean",255)
      RETURN
C****
      END SUBROUTINE io_ocean
c
      SUBROUTINE CHECKO(SUBR)
!@sum  CHECKO Checks whether Ocean are reasonable
!@ver  1.0
!!      USE MODEL_COM, only : im,jm
!!      USE FLUXES, only : gtemp
!!      USE MODEL_COM, only : focean
      USE DOMAIN_DECOMP, only: AM_I_ROOT
      USE HYCOM_ATM, only : gather_atm, focean,gtemp,im,jm
      IMPLICIT NONE
      integer i,j

!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

      call gather_atm
      if (AM_I_ROOT()) then

      print *,'SUBR=',SUBR
      write(*,'(10f7.2)') ((gtemp(1,1,i,j),i=1,10),j=15,20)
      write(*,'(a)') 'focean'
      write(*,'(10f7.2)') ((focean(i,j),i=1,10),j=15,20)

      endif ! AM_I_ROOT
      ! no need to sctter since nothing changed

      END SUBROUTINE CHECKO
c

      SUBROUTINE daily_OCEAN
!@sum  daily_OCEAN performs the daily tasks for the ocean module
!@auth Original Development Team
!@ver  1.0
C****
      implicit none
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
      implicit none
      return
      end subroutine gather_odiags


      subroutine alloc_ocean
      USE DOMAIN_DECOMP, only: AM_I_ROOT
      USE HYCOM_DIM, only : init_hycom_grid, alloc_hycom_dim
      USE HYCOM_ARRAYS, only : alloc_hycom_arrays

      USE HYCOM_DIM_GLOB, only : alloc_hycom_dim_glob
      USE HYCOM_ARRAYS_GLOB, only : alloc_hycom_arrays_glob

      USE KPRF_ARRAYS, only : alloc_kprf_arrays, alloc_kprf_arrays_local
      USE HYCOM_ATM, only : alloc_hycom_atm

#ifdef TRACERS_GASEXCH_Natassa
      USE TRACER_GASEXCH_COM, only: alloc_tracer_gasexch_com
#endif

#ifdef TRACERS_OceanBiology
      USE obio_forc, only: alloc_obio_forc
      USE obio_com,  only: alloc_obio_com
#endif

      implicit none
      

      ! seems like this is ok place to create ocean grid since nobody
      ! uses it before this call...
      call init_hycom_grid

      

      call alloc_hycom_atm

      call alloc_hycom_dim
      call alloc_hycom_arrays

      call alloc_hycom_dim_glob
      call alloc_hycom_arrays_glob

      call alloc_kprf_arrays
      call alloc_kprf_arrays_local

#ifdef TRACERS_GASEXCH_Natassa
      call alloc_tracer_gasexch_com
#endif
#ifdef TRACERS_OceanBiology
      call alloc_obio_forc
      call alloc_obio_com
#endif

      !!call reset_hycom_arrays

      !if (AM_I_ROOT()) then
 !!!       call geopar
      !endif


      end subroutine alloc_ocean

#ifdef THIS_PART_IS_NOT_READY
      subroutine reset_hycom_arrays
      USE HYCOM_DIM_GLOB, only : ii,jj,kk
      USE HYCOM_ARRAYS_GLOB
      USE KPRF_ARRAYS, only : sswflx
      USE HYCOM_SCALARS, only : lp,huge
      implicit none

      integer i,j,k,ja,jb,ia
      real :: zero = 0.
      

      write (*,*) 'laying out arrays in memory ...'
c$OMP PARALLEL DO SCHEDULE(STATIC,jchunk)
      do 209 j=1,jj
      do 209 i=1,ii
      p(i,j,:)=huge
      pv(i,j,1)=huge
      pbot(i,j)=huge
      ubavg(i,j,:)=huge
      vbavg(i,j,:)=huge
       utotm(i,j)=huge
      vtotm(i,j)=huge
      utotn(i,j)=huge
      vtotn(i,j)=huge
      uflux (i,j)=huge
      vflux (i,j)=huge
      uflux1(i,j)=huge
      vflux1(i,j)=huge
      uflux2(i,j)=huge
      vflux2(i,j)=huge
      uflux3(i,j)=huge
      vflux3(i,j)=huge
      uja(i,j)=huge
      ujb(i,j)=huge
      via(i,j)=huge
      vib(i,j)=huge
      pgfx(i,j)=huge
      pgfy(i,j)=huge
      depthu(i,j)=huge
      depthv(i,j)=huge
      tprime(i,j)=huge
c
      srfhgt(i,j)=zero
      dpmixl(i,j)=zero
      oice(i,j)=zero
      taux(i,j)=zero
      tauy(i,j)=zero
      oflxa2o(i,j)=zero
      osalt(i,j)=zero
      oemnp(i,j)=zero
      ustar(i,j)=zero
      sswflx(i,j)=zero
c
      ubavav(i,j)=zero
      vbavav(i,j)=zero
      pbavav(i,j)=zero
      sfhtav(i,j)=zero
      dpmxav(i,j)=zero
      oiceav(i,j)=zero
      eminpav(i,j)=zero
      surflav(i,j)=zero
      sflxav(i,j)=zero
      brineav(i,j)=zero
c
      u  (i,j,:   )=huge
      v  (i,j,:   )=huge
      uflx(i,j,:)=huge
      vflx(i,j,:)=huge
      ufxcum(i,j,:)=huge
      vfxcum(i,j,:)=huge
      dpinit(i,j,:)=huge
      dpold (i,j,:)=huge
      dp (i,j,:   )=huge
      dpu(i,j,k   )=huge
      dpv(i,j,k   )=huge
      p (i,j,:)=huge
      pu(i,j,:)=huge
      pv(i,j,:)=huge
c
      th3d(i,j,:)=huge
      thermb(i,j,:)=huge
      thstar(i,j,:)=huge
!      do nt=1,ntrcr
        tracer(i,j,:,:)=zero
!      end do
      uav(i,j,:)=zero
      vav(i,j,:)=zero
      dpuav(i,j,:)=zero
      dpvav(i,j,:)=zero
      dpav (i,j,:)=zero
      temav(i,j,:)=zero
      salav(i,j,:)=zero
      th3av(i,j,:)=zero
      uflxav(i,j,:)=zero
      vflxav(i,j,:)=zero
      diaflx(i,j,:)=zero
 209  continue
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(ja) SCHEDULE(STATIC,jchunk)
      do 210 j=1,jj
      !!ja=mod(j-2+jj,jj)+1
      !!do 210 l=1,isq(j)
      !!do 210 i=ifq(j,l),ilq(j,l)
      do i=1,ii
      pbot(i  ,j  )=0.
      !pbot(i-1,j  )=0.
      !pbot(i  ,ja )=0.
      !pbot(i-1,ja )=0.
      p(i  ,j  ,1)=0.
      !p(i-1,j  ,1)=0.
      !p(i  ,ja ,1)=0.
      !p(i-1,ja ,1)=0.
      do 210 k=1,kk
      dp(i  ,j  ,:   )=0.
      dp(i  ,j  ,k+kk)=0.
      !dp(i-1,j  ,k   )=0.
      !dp(i-1,j  ,k+kk)=0.
      !dp(i  ,ja ,k   )=0.
      !dp(i  ,ja ,k+kk)=0.
      !dp(i-1,ja ,k   )=0.
 !210  !dp(i-1,ja ,k+kk)=0.
 210  continue
c$OMP END PARALLEL DO
c
c --- initialize  u,ubavg,utotm,uflx,uflux,uflux2/3,uja,ujb  at points
c --- located upstream and downstream (in i direction) of p points.
c --- initialize  depthu,dpu,utotn,pgfx  upstream and downstream of p points
c --- as well as at lateral neighbors of interior u points.
c
c$OMP PARALLEL DO PRIVATE(ja,jb) SCHEDULE(STATIC,jchunk)
      do 156 j=1,jj
      do 156 i=1,ii !ifu(j,l),ilu(j,l)
      pu(i,j,:)=0.
c
      depthu(i,j)=0.
      utotn (i,j)=0.
      pgfx  (i,j)=0.
      dpu(i,j,:   )=0.
c
 156  continue
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO SCHEDULE(STATIC,jchunk)
      do 158 j=1,jj
      !do 158 l=1,isp(j)
      do 158 i=1,ii !ifp(j,l),ilp(j,l)+1
      depthu(i,j)=0.
      utotn (i,j)=0.
      pgfx  (i,j)=0.
      ubavg(i,j,1)=0.
      ubavg(i,j,2)=0.
      ubavg(i,j,3)=0.
      utotm (i,j)=0.
      uflux (i,j)=0.
      uflux2(i,j)=0.
      uflux3(i,j)=0.
      uja(i,j)=0.
      ujb(i,j)=0.
c
      dpu(:,:,:)=0.
      uflx(:,:,:)=0.
      ufxcum(:,:,:)=0.
      u(:,:,:)=0.
c$OMP END PARALLEL DO
c
c --- initialize  v,vbavg,vtotm,vflx,vflux,vflux2/3,via,vib  at points
c --- located upstream and downstream (in j direction) of p points.
c --- initialize  depthv,dpv,vtotn,pgfy  upstream and downstream of p points
c --- as well as at lateral neighbors of interior v points.
c
      pv(:,:,:)=0.
c
      depthv(:,:)=0.
      vtotn (:,:)=0.
      pgfy  (:,:)=0.
c
      depthv(:,:)=0.
      vtotn (:,:)=0.
      pgfy  (:,:)=0.
c
      dpv(:,:,:   )=0.
c
      depthv(:,:)=0.
      vtotn (:,:)=0.
      pgfy  (:,:)=0.
      vbavg(:,:,:)=0.
      vtotm (:,:)=0.
      vflux (:,:)=0.
      vflux2(:,:)=0.
      vflux3(:,:)=0.
      via(:,:)=0.
      vib(:,:)=0.
c
      dpv(:,:,:   )=0.
      vflx(:,:,:)=0.
      vfxcum(:,:,:)=0.
      v(:,:,:   )=0.
      write (*,*) '... array layout completed'
css   endif                    ! end of nstep=0

      end subroutine reset_hycom_arrays
#endif
