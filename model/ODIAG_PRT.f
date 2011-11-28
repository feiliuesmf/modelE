#include "rundeck_opts.h"

      subroutine diag_ocean_prep
      use oceanr_dim, only : grid=>ogrid
      implicit none
      if(.not. grid%have_domain) return
      call basin_prep
      call oijl_prep
      return
      end subroutine diag_ocean_prep

      SUBROUTINE diag_OCEAN
#ifndef NEW_IO
!@sum  diag_OCEAN prints out diagnostics for ocean
!@auth Gavin Schmidt/Gary Russell
C**** Note this is an incorporation and modification of the stand alone
C**** ocean diagnostic programs from Gary. All diagnostics are on the
C**** ocean grid.
      USE CONSTANT, only : undef,teeny
      USE MODEL_COM, only : xlabel,lrunid,jmon0,jyear0,idacc,jdate0
     *     ,amon0,jdate,amon,jyear
      USE OCEAN, only : im,jm,lmo,dts,imaxj,lmm,ze
      USE DIAG_COM, only : qdiag,zoc_pout=>zoc,zoc1_pout=>zoc1
      USE MDIAG_COM, only : acc_period
      USE ODIAG
      USE FILEMANAGER, only : openunit
      IMPLICIT NONE
      INTEGER I,J,L,N,NOL(LMO),JEQ,JDLAT,KXLB
      REAL*8 DLON

      zoc_pout(1:lmo) = zoc(1:lmo)
      zoc1_pout(1:lmo+1) = zoc1(1:lmo+1)

C**** Calculate latitudes
      JEQ = JM/2
      JDLAT = NINT(180./(JM-1))
      DO J=1,JM
        FLAT(J,1)  = JDLAT*(J-JEQ-0.5)  ! primary grid
      END DO
      FLAT(1,2)=-90.
      DO J=1,JM-1
        FLAT(J+1,2)  = JDLAT*(J-JEQ)  ! secondary grid  (shifted by 1)
      END DO
C**** Calculate longitudes
      DLON=360./REAL(IM,KIND=8)
      DO I=1,IM
        FLON(I,1)  = -180.+DLON*(I-0.5)  ! primary grid
      END DO
      DO I=1,IM
        FLON(I,2)  = -180.+DLON*I        ! secondary grid
      END DO

C**** determine label to be added to all titles
      KXLB = INDEX(XLABEL(1:11),'(')-1
      IF(KXLB.le.0) KXLB = 10
      XLB = ' '
      XLB(1:13)=acc_period(1:3)//' '//acc_period(4:12)
      XLB = TRIM(XLB)//" "//XLABEL(1:KXLB)

C**** Open output files
      IF(QDIAG) then
        call open_ij(trim(acc_period)//'.oij'//XLABEL(1:LRUNID),im,jm)
        call open_jl(trim(acc_period)//'.ojl'//XLABEL(1:LRUNID),jm,lmo,0
     *       ,flat)
        call open_il(trim(acc_period)//'.oil'//XLABEL(1:LRUNID),im,lmo,0
     *       )
        call openunit(trim(acc_period)//'.otj'//XLABEL(1:LRUNID),iu_otj,
     *       .FALSE.,.FALSE.)
      END IF
C**** Call diagnostic calculating routines
      CALL OIJOUT  ! lat-lon diags
      CALL OSFOUT  ! overturning streamfunction diags
      CALL OJLOUT  ! lat-height and lon-height sections and means
      WRITE (6,*)
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
      CALL OTJOUT  ! transports
      CALL STROUT  ! strait diags
C**** Miscellaneous diags for print out

C****
C**** vertical mean diagnostics
C****
      WRITE (6,*)
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
C**** Calculate number of points in average
      NOL = 0.
      DO J=1,JM
      DO I=1,IMAXJ(J)
      DO L=1,LMM(I,J)
        NOL(L)=NOL(L)+1
      END DO
      END DO
      END DO

      WRITE(6,*) " Ocean Mean quantities:"
      WRITE(6,*) " Level   Rho0       Temp    Salinity"
      DO L=1,LMO
        WRITE(6,'(2X,I4,10F10.3)') L,OL(L,L_RHO)/(NOL(L)*IDACC(1)),
     *           OL(L,L_TEMP)/(NOL(L)*IDACC(1)),
     *       1d3*OL(L,L_SALT)/(NOL(L)*IDACC(1))
      END DO
      call sys_flush(6)
C****
      IF (QDIAG) THEN
        call close_ij
        call close_jl
        call close_il
      END IF
C****
      RETURN
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
#endif
      END SUBROUTINE diag_OCEAN

#ifndef NEW_IO
      SUBROUTINE OIJOUT
!@sum  OIJOUT prints out lat-lon diagnostics for ocean
!@auth Gavin Schmidt/Gary Russell
      USE CONSTANT, only : undef,teeny,rhows
      USE MODEL_COM, only : xlabel,lrunid,jmon0,jyear0,idacc,jdate0
     *     ,amon0,jdate,amon,jyear
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : ntm,trw0,trname,ntrocn,n_water,n_obio
#endif
      USE OCEAN, only: im,jm,lmo,focean,dxypo,dts,imaxj,lmm,ze,dxvo,dypo
      USE DIAG_COM, only : qdiag
      USE MDIAG_COM, only :
     &     sname_strlen,units_strlen,lname_strlen
      USE STRAITS, only : nmst,wist,dist,lmst,name_st

      USE OCEAN, only : oDLAT_DG, oLAT_DG, oLON_DG

      USE ODIAG
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : to_per_mil
#endif
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,JM) :: Q,SFIJ,QS,QS700
      REAL*8, DIMENSION(IM,JM,LMO) :: Q3
c for now we are assuming that igrid=jgrid in arguments to pout_ij
      INTEGER I,J,K,L,N,KB,IJGRID,IP1,k1,KK
      INTEGER :: LMINEF=1, LMAXEF=LMO, KVMF(3) = (/ 3, 6, 9/),
     *           LMINMF=1, LMAXMF=1,   KCMF(3) = (/ 3, 6, 9/),
     *           LMINSF=1, LMAXSF=LMO, KVDC(3) = (/ 3, 6, 9/)
     *     ,KCMFfull(13) = (/1,2,3,4,5,6,7,8,9,10,11,12,13/)
      REAL*8 GOS,SOS,FAC,FACST,GSMAX,GSMIN,CKMIN,CKMAX,ACMIN,ACMAX
     *     ,TSUM,TEMGS,QJ(JM),QSUM,MQ,DLON,byiacc,volgs
      CHARACTER NAME(KOLNST)*40,TITLE*80
      CHARACTER(len=lname_strlen) :: lname
      CHARACTER(len=sname_strlen) :: sname
      CHARACTER(len=units_strlen) :: units
      character*50 :: unit_string
      character(len=2), dimension(lmo) :: levstr
      real*8 :: scale_jm2

      QJ=0.
      QSUM=0.
      IJGRID=1
      Q=0.
      Q3=0.

C**** lat/lon diagnostics
      WRITE (6,*)
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)

      IF(QDIAG) then

      DO L=1,LMO
        WRITE(LEVSTR(L),'(I2.2)') L
      ENDDO

C****
C**** Ocean Potential Temperature (C)
C****
      K=IJL_PTM
      LNAME=LNAME_OIJL(K)
      UNITS=UNITS_OIJL(K)
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      TITLE(51:80)=XLB
C**** Loop over layers
      DO L=1,LMO
      DO J=1,JM
      DO I=1,IMAXJ(J)
        Q(I,J) = UNDEF
        IF(FOCEAN(I,J).gt..5 .and. OIJL(I,J,L,IJL_MO).gt.0.)  THEN
          GOS = OIJL(I,J,L,IJL_G0M) / (OIJL(I,J,L,IJL_MO)*DXYPO(J))
          SOS = OIJL(I,J,L,IJL_S0M) / (OIJL(I,J,L,IJL_MO)*DXYPO(J))
          Q(I,J) = TEMGS(GOS,SOS)
        END IF
      END DO
      END DO
      Q(2:IM,JM)=Q(1,JM)
      Q(2:IM,1)=Q(1,1)
      WRITE (LNAME(40:47),'(A5,I3)') 'Level',L
      WRITE (TITLE(40:47),'(A5,I3)') 'Level',L
      SNAME = 'oc_temp_L'//LEVSTR(L)
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO
C****
C**** Ocean Salinity (per mil)
C****
      K=IJL_S0M
      LNAME=LNAME_OIJL(K)
      UNITS=UNITS_OIJL(K)
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      TITLE(51:80)=XLB
C**** Loop over layers
      DO L=1,LMO
      DO J=1,JM
      DO I=1,IMAXJ(J)
        Q(I,J) = UNDEF
        IF(FOCEAN(I,J).gt..5 .and. OIJL(I,J,L,IJL_MO).gt.0.)
     *       Q(I,J) = OIJL(I,J,L,K)*SCALE_OIJL(K) /
     &       (OIJL(I,J,L,IJL_MO)*DXYPO(J))
      END DO
      END DO
      Q(2:IM,JM)=Q(1,JM)
      Q(2:IM,1)=Q(1,1)
      WRITE (LNAME(40:47),'(A5,I3)') 'Level',L
      WRITE (TITLE(40:47),'(A5,I3)') 'Level',L
      SNAME = 'oc_salt_L'//LEVSTR(L)
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO
C****
C**** Ocean Potential Density (kg/m^3) (w.r.t. 0m)
C****
      K=IJL_PDM
      LNAME=LNAME_OIJL(K)
      UNITS=UNITS_OIJL(K)
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      TITLE(51:80)=XLB
C**** Loop over layers
      DO L=1,LMO
      DO J=1,JM
      DO I=1,IMAXJ(J)
        Q(I,J) = UNDEF
        IF(FOCEAN(I,J).gt..5 .and. OIJL(I,J,L,IJL_MO).gt.0.)  THEN
          GOS = OIJL(I,J,L,IJL_G0M) / (OIJL(I,J,L,IJL_MO)*DXYPO(J))
          SOS = OIJL(I,J,L,IJL_S0M) / (OIJL(I,J,L,IJL_MO)*DXYPO(J))
          Q(I,J) = 1./VOLGS(GOS,SOS) - 1000.
        END IF
      END DO
      END DO
      Q(2:IM,JM)=Q(1,JM)
      Q(2:IM,1)=Q(1,1)
      WRITE (LNAME(40:47),'(A5,I3)') 'Level',L
      WRITE (TITLE(40:47),'(A5,I3)') 'Level',L
      SNAME = 'oc_pot_den_L'//LEVSTR(L)
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO
#ifdef TRACERS_OCEAN
C****
C**** Ocean Tracers
C****
      DO N=1,NTM
      LNAME="OCEAN "//trname(n)
      if (to_per_mil(n).gt.0) THEN
        UNITS="per mil"
      ELSE
        UNITS=unit_string(ntrocn(n),'kg/kg')
      END IF
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      TITLE(51:80)=XLB
C**** Loop over layers
      DO L=1,LMO
      DO J=1,JM
      DO I=1,IMAXJ(J)
        Q(I,J) = UNDEF
        IF(FOCEAN(I,J).gt..5 .and. OIJL(I,J,L,IJL_MO).gt.0.) THEN
          if (to_per_mil(n).gt.0 .and. TOIJL(I,J,L,TOIJL_CONC
     *         ,n_water).gt.0) THEN
            Q(I,J)=1d3*(TOIJL(I,J,L,TOIJL_CONC,N)/(TOIJL(I,J,L
     *           ,TOIJL_CONC,n_water)*trw0(n))-1.)
c          Q(I,J)=1d3*(TOIJL(I,J,L,TOIJL_CONC,N)/((OIJL(I,J,L,IJL_MO)
c     *         *DXYPO(J)-OIJL(I,J,L,IJL_S0M))*trw0(n))-1.)
          else
            Q(I,J)=10.**(-ntrocn(n))*TOIJL(I,J,L,TOIJL_CONC,N)/
     *           (OIJL(I,J,L,IJL_MO)*DXYPO(J))
          end if
        END IF
      END DO
      END DO
      Q(2:IM,JM)=Q(1,JM)
      Q(2:IM,1)=Q(1,1)
      WRITE (LNAME(40:47),'(A5,I3)') 'Level',L
      WRITE (TITLE(40:47),'(A5,I3)') 'Level',L
      SNAME='oc_'//trim(trname(n))//'_L'//LEVSTR(L)
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO
      END DO
#endif
C****
C**** East-West or North-South Velocities (cm/s)
C****
      K=IJL_MFU
      DO LMINMF=1,LMO
        LNAME=LNAME_OIJL(K)
        UNITS=UNITS_OIJL(K)
        LMAXMF=LMINMF
      Q = 0.
      DO J=1,JM
        I=IM
        DO IP1=1,IMAXJ(J)
          MQ = 0.
          DO L=LMINMF,LMAXMF
            MQ = MQ + (OIJL(I,J,L,IJL_MO)+OIJL(IP1,J,L,IJL_MO))
            Q(I,J) = Q(I,J) + OIJL(I,J,L,K)
          END DO
          Q(I,J) = SCALE_OIJL(K) * Q(I,J) / (.5*MQ*DYPO(J)+teeny)
          I=IP1
        END DO
      END DO
      Q(2:IM,JM)=Q(1,JM)
      Q(2:IM,1)=Q(1,1)
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      IF(LMINMF.eq.LMAXMF) THEN
        WRITE (TITLE(39:41),'(I3)') LMINMF
        WRITE (LNAME(39:41),'(I3)') LMINMF
      ELSEIF(LMINMF.lt.LMAXMF) THEN
        WRITE (TITLE(39:46),'(I3,A2,I3)') LMINMF," -",LMAXMF
        WRITE (LNAME(39:46),'(I3,A2,I3)') LMINMF," -",LMAXMF
      END IF
      SNAME="uvel_L"//LEVSTR(LMINMF)
      TITLE(51:80)=XLB
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO

      K=IJL_MFV
      DO LMINMF=1,LMO
        LNAME=LNAME_OIJL(K)
        UNITS=UNITS_OIJL(K)
        LMAXMF=LMINMF
      Q = 0.
      DO J=1,JM-1
        DO I=1,IMAXJ(J)
          MQ = 0.
          DO L=LMINMF,LMAXMF
            MQ = MQ + (OIJL(I,J,L,IJL_MO)+OIJL(I,J+1,L,IJL_MO))
            Q(I,J) = Q(I,J) + OIJL(I,J,L,K)
          END DO
          Q(I,J) =  SCALE_OIJL(K) * Q(I,J) / (.5*MQ*DXVO(J)+teeny)
        END DO
      END DO
      Q(2:IM,1)=Q(1,1)
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      IF(LMINMF.eq.LMAXMF) THEN
        WRITE (TITLE(39:41),'(I3)') LMINMF
        WRITE (LNAME(39:41),'(I3)') LMINMF
      ELSEIF(LMINMF.lt.LMAXMF) THEN
        WRITE (TITLE(39:46),'(I3,A2,I3)') LMINMF," -",LMAXMF
        WRITE (LNAME(39:46),'(I3,A2,I3)') LMINMF," -",LMAXMF
      END IF
      SNAME="vvel_L"//LEVSTR(LMINMF)
      TITLE(51:80)=XLB
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,2,2)
      END DO
C****
C**** Vertical Velocity (cm/s)
C****
      K=IJL_MFW
      DO L=1,LMO-1
c     IF(KVMF(K).le.0)  GO TO 370
c        L =KVMF(K)
        LNAME=LNAME_OIJL(K)
        UNITS=UNITS_OIJL(K)
        SNAME='vert_vel_L'//LEVSTR(L)
        DO J=1,JM
          DO I=1,IMAXJ(J)
            Q(I,J)=OIJL(I,J,L,K)*SCALE_OIJL(K)/(IDACC(1)*DXYPO(J))
          END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        WRITE (TITLE(40:47),'(A5,I3)') "Level",L
        WRITE (LNAME(40:47),'(A5,I3)') "Level",L
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO
C****
C**** Vertical Mass Flux^2 (kg2/m4)
C****
      K=IJL_MFW2
      DO L=1,LMO-1
c     IF(KVMF(K).le.0)  GO TO 370
c        L =KVMF(K)
        LNAME=LNAME_OIJL(K)
        UNITS=UNITS_OIJL(K)
        SNAME='vert_ms_flx_sq_L'//LEVSTR(L)
        DO J=1,JM
          DO I=1,IMAXJ(J)
            Q(I,J)=OIJL(I,J,L,K)*SCALE_OIJL(K)/
     *           (IDACC(1)*DXYPO(J)*DXYPO(J))
          END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        WRITE (TITLE(40:47),'(A5,I3)') "Level",L
        WRITE (LNAME(40:47),'(A5,I3)') "Level",L
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO
C****
C**** East-West or North-South Mass fluxes (kg/s)
C****
      K=IJL_MFU
      DO LMINMF=1,LMO
      LNAME="EAST-WEST MASS FLUX"
      UNITS="10^9 kg/s"
        LMAXMF=LMINMF
      Q = 0.
      DO J=1,JM
        I=IM
        DO IP1=1,IMAXJ(J)
          DO L=LMINMF,LMAXMF
            Q(I,J) = Q(I,J) + OIJL(I,J,L,K)
          END DO
          Q(I,J) = 1d-9* Q(I,J) / (IDACC(1)*DTS+teeny)
          I=IP1
        END DO
      END DO
      Q(2:IM,JM)=Q(1,JM)
      Q(2:IM,1)=Q(1,1)
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      IF(LMINMF.eq.LMAXMF) THEN
        WRITE (TITLE(40:47),'(A5,I3)') "Level",LMINMF ! anl
        WRITE (LNAME(40:47),'(A5,I3)') "Level",LMINMF ! anl
      ELSEIF(LMINMF.lt.LMAXMF) THEN
        WRITE (TITLE(39:46),'(I3,A2,I3)') LMINMF," -",LMAXMF
        WRITE (LNAME(39:46),'(I3,A2,I3)') LMINMF," -",LMAXMF
      END IF
      SNAME="umfl_L"//LEVSTR(LMINMF)
      TITLE(51:80)=XLB
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO

      K=IJL_MFV
      DO LMINMF=1,LMO
      LNAME="NORTH-SOUTH MASS FLUX"
      UNITS="10^9 kg/s"
        LMAXMF=LMINMF
      Q = 0.
      DO J=1,JM-1
        DO I=1,IMAXJ(J)
          DO L=LMINMF,LMAXMF
            Q(I,J) = Q(I,J) + OIJL(I,J,L,K)
          END DO
          Q(I,J) = 1d-9* Q(I,J) / (IDACC(1)*DTS+teeny)
        END DO
      END DO
      Q(2:IM,1)=Q(1,1)
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      IF(LMINMF.eq.LMAXMF) THEN
        WRITE (TITLE(40:47),'(A5,I3)') "Level",LMINMF ! anl
        WRITE (LNAME(40:47),'(A5,I3)') "Level",LMINMF ! anl
      ELSEIF(LMINMF.lt.LMAXMF) THEN
        WRITE (TITLE(39:46),'(I3,A2,I3)') LMINMF," -",LMAXMF
        WRITE (LNAME(39:46),'(I3,A2,I3)') LMINMF," -",LMAXMF
      END IF
      SNAME="vmfl_L"//LEVSTR(LMINMF)
c        WRITE (TITLE(40:47),'(A5,I3)') "Level",L ! anl
c        WRITE (LNAME(40:47),'(A5,I3)') "Level",L ! anl
      TITLE(51:80)=XLB
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,2,2)
      END DO
C****
C**** Vertical Mass Flux (kg/s)
C****
      DO L=1,LMO-1
c     IF(KVMF(K).le.0)  GO TO 370
c        L =KVMF(K)
        LNAME="DOWNWARD VERTICAL MASS FLUX"
        UNITS="10^9 kg/s"
        SNAME='vert_mfl_L'//LEVSTR(L)
        DO J=1,JM
          DO I=1,IMAXJ(J)
            Q(I,J) = 1d-9*OIJL(I,J,L,IJL_MFW) / real(IDACC(1)*DTS)
          END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        WRITE (TITLE(40:47),'(A5,I3)') "Level",L
        WRITE (LNAME(40:47),'(A5,I3)') "Level",L
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO
C****
C**** East-West or North-South Heat Flux (10^15 W)
C****
      DO K=IJL_GFLX,IJL_GFLX+1
c      IF(.not.QL(K))  GO TO 440
        LNAME=LNAME_OIJL(K)
        UNITS=UNITS_OIJL(K)
        IF (K.eq.IJL_GFLX) THEN
          SNAME="ew_hflx"
        ELSE
          SNAME="ns_hflx"
        END IF
        Q = 0.
        DO J=1,JM
        DO I=1,IMAXJ(J)
          DO L=LMINEF,LMAXEF
            Q(I,J) = Q(I,J) + OIJL(I,J,L,K)
          END DO
          Q(I,J) = SCALE_OIJL(K)*Q(I,J)/IDACC(1)
        END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        IF(LMINEF.eq.LMAXEF) WRITE (TITLE(39:41),'(I3)') LMINEF
        IF(LMINEF.lt.LMAXEF) WRITE (TITLE(39:46),'(I3,A2,I3)') LMINEF,
     *       " -",LMAXEF
        IF(LMINEF.eq.LMAXEF) WRITE (LNAME(39:41),'(I3)') LMINEF
        IF(LMINEF.lt.LMAXEF) WRITE (LNAME(39:46),'(I3,A2,I3)') LMINEF,
     *       " -",LMAXEF
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO
C****
C**** East-West or North-South Salt Flux (10^6 kg/s)
C****
      DO K=IJL_SFLX,IJL_SFLX+1
c      IF(.not.QL(K))  GO TO 540
        LNAME=LNAME_OIJL(K)
        UNITS=UNITS_OIJL(K)
        IF (K.eq.IJL_SFLX) THEN
          SNAME="ew_sflx"
        ELSE
          SNAME="ns_sflx"
        END IF
        Q = 0.
        DO J=1,JM
        DO I=1,IMAXJ(J)
          DO L=LMINSF,LMAXSF
            Q(I,J) = Q(I,J) + OIJL(I,J,L,K)
          END DO
          Q(I,J) = SCALE_OIJL(K)*Q(I,J)/IDACC(1)
        END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        IF(LMINSF.eq.LMAXSF) WRITE (TITLE(39:41),'(I3)') LMINSF
        IF(LMINSF.lt.LMAXSF) WRITE (TITLE(39:46),'(I3,A2,I3)') LMINSF,
     *       " -",LMAXSF
        IF(LMINSF.eq.LMAXSF) WRITE (LNAME(39:41),'(I3)') LMINSF
        IF(LMINSF.lt.LMAXSF) WRITE (LNAME(39:46),'(I3,A2,I3)') LMINSF,
     *       " -",LMAXSF
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO
C****
C**** Gent-McWilliams fluxes ([J or kg]/s or [J kg]/m2/s for vert)
C****
      DO KK=0,2
        K=KK+IJL_GGMFL
        DO L=1,lmo 
          LNAME=LNAME_OIJL(K)
          UNITS=UNITS_OIJL(K)
          SELECT CASE (KK)
          CASE (0)      ! E-W fluxes
            SNAME='gm_ew_hflx_L'//LEVSTR(L)
            DO J=1,JM
              DO I=1,IMAXJ(J)
                Q(I,J) = SCALE_OIJL(K)*OIJL(I,J,L,K)/IDACC(1)
              END DO
            END DO
          CASE (1)  ! N-S fluxes
            SNAME='gm_ns_hflx_L'//LEVSTR(L)
            DO J=1,JM
              DO I=1,IMAXJ(J)
                Q(I,J) = SCALE_OIJL(K)*OIJL(I,J,L,K)/IDACC(1)
              END DO
            END DO
          CASE (2)    !  Vertical fluxes
            SNAME='gm_vt_hflx_L'//LEVSTR(L)
            DO J=1,JM
              DO I=1,IMAXJ(J)
                Q(I,J) = SCALE_OIJL(K)*OIJL(I,J,L,K)/(IDACC(1)*DXYPO(J))
              END DO
            END DO
          END SELECT
          Q(2:IM,JM)=Q(1,JM)
          Q(2:IM,1)=Q(1,1)
          TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
          WRITE (TITLE(40:47),'(A5,I3)') "Level",L
          WRITE (LNAME(40:47),'(A5,I3)') "Level",L
          TITLE(51:80)=XLB
          CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
        END DO
      END DO
C****
C**** Gent-McWilliams Salt Fluxes
C****
      DO KK=0,2
        K=KK+IJL_SGMFL
        DO L=1,lmo
          LNAME=LNAME_OIJL(K)
          UNITS=UNITS_OIJL(K)
          SELECT CASE (KK)
          CASE (0)      ! E-W fluxes
            SNAME='gm_ew_sflx_L'//LEVSTR(L)
            DO J=1,JM
              DO I=1,IMAXJ(J)
                Q(I,J) = SCALE_OIJL(K)*OIJL(I,J,L,K)/IDACC(1)
              END DO
            END DO
          CASE (1)  ! N-S fluxes
            SNAME='gm_ns_sflx_L'//LEVSTR(L)
            DO J=1,JM
              DO I=1,IMAXJ(J)
                Q(I,J) = SCALE_OIJL(K)*OIJL(I,J,L,K)/IDACC(1)
              END DO
            END DO
          CASE (2)    !  Vertical fluxes
            SNAME='gm_vt_sflx_L'//LEVSTR(L)
            DO J=1,JM
              DO I=1,IMAXJ(J)
                Q(I,J) = SCALE_OIJL(K)*OIJL(I,J,L,K)/(IDACC(1)*DXYPO(J))
              END DO
            END DO
          END SELECT

          Q(2:IM,JM)=Q(1,JM)
          Q(2:IM,1)=Q(1,1)
          TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
          WRITE (TITLE(40:47),'(A5,I3)') "Level",L
          WRITE (LNAME(40:47),'(A5,I3)') "Level",L
          TITLE(51:80)=XLB
          CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
        END DO
      END DO
#ifdef TRACERS_OCEAN
C****
C**** Gent-McWilliams Tracer Fluxes
C****
      do n=1,ntm
      DO KK=0,2
        DO L=1,lmo
          SELECT CASE (KK)
          CASE (0)      ! E-W fluxes
            LNAME="GM/EDDY E-W FLUX "//trname(n)
            UNITS=unit_string(ntrocn(n),'kg/s')
            SNAME="oc_gm_ewflx"//trim(trname(n))//"_L"//LEVSTR(L)
          CASE (1)  ! N-S fluxes
            LNAME="GM/EDDY N-S FLUX "//trname(n)
            UNITS=unit_string(ntrocn(n),'kg/s')
            SNAME="oc_gm_nstflx"//trim(trname(n))//"_L"//LEVSTR(L)
          CASE (2)    !  Vertical fluxes
            LNAME="GM/EDDY VERT. FLUX "//trname(n)
            UNITS=unit_string(ntrocn(n)-6,'kg/m^2 s')
            SNAME="gm_vt_tflx"//trim(trname(n))//"_L"//LEVSTR(L)
          END SELECT
          if (KK.EQ.3) THEN     ! vert fluxes, scale/divide by area

          DO J=1,JM
            DO I=1,IMAXJ(J)
              IF(FOCEAN(I,J).gt..5 .and. OIJL(I,J,L,IJL_MO).gt.0.)
     *             THEN
                if (TOIJL(I,J,L,TOIJL_CONC,N).gt.0) Q(I,J)=10.**
     *                (6-ntrocn(n))*TOIJL(I,J,L,KK+TOIJL_GMFL,N)/
     *                (IDACC(1)*DTS*DXYPO(J))
              ENDIF
            END DO
          END DO

          ELSE
          
          DO J=1,JM
            DO I=1,IMAXJ(J)
              IF(FOCEAN(I,J).gt..5 .and. OIJL(I,J,L,IJL_MO).gt.0.)
     *             THEN
                if (TOIJL(I,J,L,TOIJL_CONC,N).gt.0) Q(I,J)=10.**
     *               (-ntrocn(n))*TOIJL(I,J,L,KK+TOIJL_GMFL,N)/(IDACC(1)
     *               *DTS)
              ENDIF
            END DO
          END DO

          END IF
          Q(2:IM,JM)=Q(1,JM)
          Q(2:IM,1)=Q(1,1)
          TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
          WRITE (TITLE(40:47),'(A5,I3)') "Level",L
          WRITE (LNAME(40:47),'(A5,I3)') "Level",L
          TITLE(51:80)=XLB
          CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
        END DO
      END DO
      enddo
#endif
C****
C**** Vertical Diffusion Coefficients (cm/s)
C****
      K=IJL_KVM
      DO L=1,lmo-1
        LNAME=LNAME_OIJL(K)
        UNITS=UNITS_OIJL(K)
        SNAME="kvm_L"//LEVSTR(L)
        Q=UNDEF
        DO J=1,JM
        DO I=1,IMAXJ(J)
          IF (OIJL(I,J,L+1,IJL_MO).gt.0) Q(I,J)=OIJL(I,J,L,K)*
     &         SCALE_OIJL(K)/IDACC(1)
        END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        WRITE (TITLE(40:47),'(A5,I3)') "Level",L
        WRITE (LNAME(40:47),'(A5,I3)') "Level",L
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,2,2)
      END DO

      K=IJL_KVG
      DO L=1,lmo-1
        LNAME=LNAME_OIJL(K)
        UNITS=UNITS_OIJL(K)
        SNAME="kvg_L"//LEVSTR(L)
        Q=UNDEF
        DO J=1,JM
        DO I=1,IMAXJ(J)
          IF (OIJL(I,J,L+1,IJL_MO).gt.0) Q(I,J)=OIJL(I,J,L,K)*
     &         SCALE_OIJL(K)/IDACC(1)
        END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        WRITE (TITLE(40:47),'(A5,I3)') "Level",L
        WRITE (LNAME(40:47),'(A5,I3)') "Level",L
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,2,2)
      END DO

      K=IJL_WGFL
      DO L=1,lmo-1
        LNAME=LNAME_OIJL(K)
        UNITS=UNITS_OIJL(K)
        SNAME="wgfl_L"//LEVSTR(L)
        Q=UNDEF
        DO J=1,JM
          DO I=1,IMAXJ(J)
            IF (OIJL(I,J,L+1,IJL_MO).gt.0) Q(I,J)=OIJL(I,J,L,K)
     &           *SCALE_OIJL(K)/(IDACC(1)*dxypo(j))
          END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        WRITE (TITLE(40:47),'(A5,I3)') "Level",L
        WRITE (LNAME(40:47),'(A5,I3)') "Level",L
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,2,2)
      END DO

      K=IJL_WSFL
      DO L=1,lmo-1
        LNAME=LNAME_OIJL(K)
        UNITS=UNITS_OIJL(K)
        SNAME="wsfl_L"//LEVSTR(L)
        Q=UNDEF
        DO J=1,JM
        DO I=1,IMAXJ(J)
          IF (OIJL(I,J,L+1,IJL_MO).gt.0) Q(I,J)=OIJL(I,J,L,K)
     &         *SCALE_OIJL(K)/(IDACC(1)*dxypo(j))
        END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        WRITE (TITLE(40:47),'(A5,I3)') "Level",L
        WRITE (LNAME(40:47),'(A5,I3)') "Level",L
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,2,2)
      END DO

#ifdef TRACERS_OCEAN
      do n=1,ntm
      DO L=1,lmo-1
        LNAME="VERT. DIFF. FLUX "//trname(n)
        UNITS=unit_string(ntrocn(n)-6,'kg/m^2 s')
        SNAME="wtfltr"//trim(trname(n))//"_L"//LEVSTR(L)
        Q=UNDEF
        DO J=1,JM
        DO I=1,IMAXJ(J)
          IF((OIJL(I,J,L+1,IJL_MO).gt.0).and.(TOIJL(I,J,L,TOIJL_CONC
     *         ,N).gt.0)) Q(I,J)=10.**(6-ntrocn(n))*TOIJL(I,J,L
     *         ,TOIJL_wtfl,N)/(IDACC(1)*DTS*DXYPO(J))
        END DO
        END DO
        Q(2:IM,JM)=Q(1,JM)
        Q(2:IM,1)=Q(1,1)
        TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
        WRITE (TITLE(40:47),'(A5,I3)') "Level",L
        WRITE (LNAME(40:47),'(A5,I3)') "Level",L
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,2,2)
      END DO
      enddo
#endif

C****
C**** Simple scaled OIJ diagnostics
C****
      DO K=1,KOIJ
        if(trim(sname_oij(k)).eq.'unused') cycle
        byiacc=1./(IDACC(IA_OIJ(K))+teeny)
        if(igrid_oij(k).eq.1 .and. jgrid_oij(k).eq.1) then
          Q=UNDEF
          DO J=1,JM
            DO I=1,IMAXJ(J)
              IF (FOCEAN(I,J).gt.0.5)
     *             Q(I,J)=SCALE_OIJ(K)*OIJ(I,J,K)*byiacc
            END DO
          END DO
          Q(2:IM,JM)=Q(1,JM)
          Q(2:IM,1)=Q(1,1)
        else ! horizontal fluxes
          Q=SCALE_OIJ(K)*OIJ(:,:,K)*byiacc
        endif
        lname=lname_oij(k)
        TITLE=trim(LNAME)//" ("//trim(UNITS_OIJ(K))//") "
        TITLE(51:80)=XLB
        CALL POUT_IJ(TITLE,SNAME_OIJ(K),LNAME_OIJ(K),UNITS_OIJ(K),Q,QJ
     *       ,QSUM,IGRID_OIJ(K),JGRID_OIJ(K))

      END DO

      END IF

C****
C**** Calculate Salt Stream Function and write it
C****
c      FAC   = -1d-6/(IDACC(1)*DTS)
c      FACST = -.5E-6/(IDACC(1)*DTS)
c      CALL STRMIJ (OIJL(1,1,1,IJL_SFLX),FAC,OLNST(1,1,LN_SFLX),FACST
c     *     ,SFIJS)
C**** Subtract .035 times the Mass Stream Function
c      DO J=1,JM
c        DO I=1,IM
c        IF (SFIJS(I,J).ne.UNDEF) SFIJS(I,J) = SFIJS(I,J)-35.*SFIJM(I,J)
c        END DO
c      END DO
c      TITLE = NAME(2)//'  Run '//XLABEL(1:6)
c      WRITE(TITLE(63:80),'(A6,I4)') JMON0,JYEAR0

C**** Output Key diagnostics: Gulf Stream, ACC, Kuroshio
      SFIJ = OIJ(:,:,IJ_SF)/idacc(1)
      WRITE(6,'(A)') " Key horizontal mass stream function diags:"

      GSMAX=0 ; GSMIN=100.
      CKMAX=0 ; CKMIN=100.
      ACMAX=0 ; ACMIN=100.
      DO J=2,JM-1
        DO I=1,IM
          if (oLAT_DG(j,2).ge.24 .and. oLAT_DG(j,2).le.38 .and. oLON_DG
     $         (i,2).ge.-80 .and. oLON_DG(i,2).le.-60) then
            IF (SFIJ(I,J).GT.GSMAX) GSMAX=SFIJ(I,J)
            IF (SFIJ(I,J).LT.GSMIN) GSMIN=SFIJ(I,J)
          end if
          if (oLAT_DG(j,2).ge.24 .and. oLAT_DG(j,2).le.38 .and. oLON_DG
     $         (i,2).ge.135 .and. oLON_DG(i,2).le.155) then
            IF (SFIJ(I,J).GT.CKMAX) CKMAX=SFIJ(I,J)
            IF (SFIJ(I,J).LT.CKMIN) CKMIN=SFIJ(I,J)
          end if
          if (oLAT_DG(j,2).ge.-72 .and. oLAT_DG(j,2).le.-54
     $         .and. oLON_DG(i,2).ge.-65-0.5*oDLAT_DG .and. oLON_DG(i,2)
     $         .le. -65+0.5*oDLAT_DG) then
            IF (SFIJ(I,J).GT.ACMAX) ACMAX=SFIJ(I,J)
            IF (SFIJ(I,J).LT.ACMIN) ACMIN=SFIJ(I,J)
          end if
        END DO
      END DO

      WRITE(6,'(a,F10.3)') " Gulf Stream (MAX in (24-38N,60-80W)):",
     *     GSMAX-GSMIN
      WRITE(6,'(a,F10.3)') " Kuroshio  (MAX in (24-38N,135-150E)):",
     *     CKMAX-CKMIN
      WRITE(6,'(a,F10.3)') " ACC                 (Drakes Passage):",
     *     ACMAX-ACMIN
C****
      IF (QDIAG) THEN
C****
C**** Ocean Heat Content (J/m^2)
C****
      K=IJL_G0M
      LNAME=LNAME_OIJL(K)
      UNITS='10^6 J/m^2' ! UNITS_OIJL(K) is J/kg and SCALE_OIJL(K) is 1
      scale_jm2 = 1d-6
      TITLE=TRIM(LNAME)//" ("//TRIM(UNITS)//")"
      TITLE(51:80)=XLB
C**** Loop over layers
      QS=0. ; QS700=0.
      DO L=1,LMO
      DO J=1,JM
      DO I=1,IMAXJ(J)
        Q(I,J) = UNDEF
        IF(FOCEAN(I,J).gt..5 .and. OIJL(I,J,L,IJL_MO).gt.0.)  THEN
          Q(I,J) = SCALE_JM2*OIJL(I,J,L,K) / (DXYPO(J)*IDACC(1))
          QS(I,J) = QS(I,J)+OIJL(I,J,L,K) / DXYPO(J)
          IF (ZE(L-1).LE.700) QS700(I,J)=QS700(I,J)+OIJL(I,J,L,K)
     *         /DXYPO(J)
        END IF
      END DO
      END DO
      Q(2:IM,JM)=Q(1,JM)
      Q(2:IM,1)=Q(1,1)
      WRITE (LNAME(40:47),'(A5,I3)') 'Level',L
      WRITE (TITLE(40:47),'(A5,I3)') 'Level',L
      SNAME="oc_ht_L"//LEVSTR(L)
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,Q,QJ,QSUM,IJGRID,IJGRID)
      END DO
C**** sum to 700m
      QS700(2:IM,JM)=QS700(1,JM)
      QS700(2:IM,1)=QS700(1,1)
      QS700 = SCALE_JM2*QS700/REAL(IDACC(1))
      WRITE (LNAME(40:47),'(A8)') 'SUM 700m'
      WRITE (TITLE(40:47),'(A8)') 'SUM 700m'
      SNAME="oc_ht_700m"
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,QS700,QJ,QSUM,IJGRID,IJGRID)
C**** total sum
      QS(2:IM,JM)=QS(1,JM)
      QS(2:IM,1)=QS(1,1)
      QS = SCALE_JM2*QS/REAL(IDACC(1))
      WRITE (LNAME(40:47),'(A8)') '  SUM  '
      WRITE (TITLE(40:47),'(A8)') '  SUM  '
      SNAME="oc_ht_total"
      CALL POUT_IJ(TITLE,SNAME,LNAME,UNITS,QS,QJ,QSUM,IJGRID,IJGRID)
      END IF
C****
      RETURN
      END SUBROUTINE OIJOUT

      SUBROUTINE OSFOUT
!@sum  OSFOUT prints out streamfunction diagnostics for ocean
!@auth Gavin Schmidt/Gary Russell
      USE CONSTANT, only : undef
      USE MODEL_COM, only : xlabel,lrunid,jmon0,jyear0,idacc,jdate0
     *     ,amon0,jdate,amon,jyear
      USE OCEAN, only : im,jm,lmo,dxypo,imaxj,ze
      USE DIAG_COM, only : qdiag
      USE MDIAG_COM, only :
     &     sname_strlen,units_strlen,lname_strlen

      USE OCEAN, only : oDLAT_DG, oLAT_DG

      USE ODIAG
      IMPLICIT NONE
      REAL*8, DIMENSION(JM+3,LMO+1) :: XJL
      INTEGER I,J,K,L,NS,N,KB,IP1,LP1,JEQ
      REAL*8 FAC,FACST
      CHARACTER TITLE*80
      CHARACTER(len=lname_strlen) :: lname
      CHARACTER(len=sname_strlen) :: sname
      CHARACTER(len=units_strlen) :: units

      TITLE(51:80)=XLB
      LNAME=""
      XJL=undef
C****
C**** Write the Mass Stream Function
C****
      TITLE(9:50)=" Mass Stream Function (Sv)"
      DO KB=1,4
        TITLE(1:8)=TRIM(BASIN(KB))
        lname(1:31)=TITLE(1:31)
        sname='sf_'//trim(BASIN(KB))
        units='Sv'
        XJL(2:JM,1:LMO+1)=SFM(1:JM-1,0:LMO,KB)
        where(xjl.ne.undef) xjl=xjl/idacc(1)
        LP1=LMO+1
        IF (QDIAG) CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,2,LP1,XJL
     *       ,ZOC1,"Latitude","Depth (m)")
C****
C**** Write data to PRinT file
C****
c        WRITE (6,921) 'Southern Hemisphere',(NINT(FLAT(J,2)),J=1,JEQ)
c        DO L=0,LMO
c          WRITE (6,922) ZOC1(L+1),(NINT(SFM(J,L,KB)),J=1,JEQ)
c        END DO
c        WRITE (6,921) 'Northern Hemisphere',(NINT(FLAT(J,2)),J=JEQ,
c    *     JM-1)
c        DO L=0,LMO
c          WRITE (6,922) ZOC1(L+1),(NINT(SFM(J,L,KB)),J=JEQ,JM-1)
c        END DO
      END DO

C**** Output Key diagnostics
      do L=1,LMO-1
        do J=2,JM-1
C**** North Atl. + North Pac. overturning
          if (oLAT_DG(j+1,2).gt.48-0.5*oDLAT_DG .and. oLAT_DG(j+1,2)
     *         .lt. 48+0.5*oDLAT_DG .and. 0.5*(ZE(L)+ZE(L-1)).le.900 
     *         .and. 0.5*(ZE(L)+ZE(L+1)).gt.900) then 
             WRITE(6,'(A46,F6.2)') 
     *            " North Atlantic overturning: 900m 48N: ",
     &           SFM(j,l,1)/idacc(1)
             WRITE(6,'(A46,F6.2)') 
     *            " North Pacific overturning:  900m 48N: ",
     &            SFM(j,l,2)/idacc(1)
          end if
C**** AABW production
          if (oLAT_DG(j+1,2).ge.-52-0.5*oDLAT_DG .and. oLAT_DG(j+1,2)
     *         .lt. -52+0.5*oDLAT_DG .and. 0.5*(ZE(L)+ZE(L-1)).le.3000
     *         .and. 0.5*(ZE(L)+ZE(L+1)).gt.3000) then
             WRITE(6,'(A46,F6.2)') " Antarctic Bottom Water production:"
     *            //" 3000m 52S: ",SFM(j,l,4)/idacc(1)
          end if
        end do
      end do

C****
C**** Calculate Salt Stream Function and write it
C****
c      FAC   = -1d-6/(IDACC(1)*NDYNO*DTO)
c      FACST = -.5E-6/(IDACC(1)*NDYNO*DTO)
c     CALL STRMJL(OIJL(1,1,1,IJL_SFLX),FAC,OLNST(1,1,LN_SFLX),FACST,SFS)
cC**** Subtract .035 times the Mass Stream Function
c      TITLE(9:50)=" Salt Stream Function (Sv)"
c      DO KB=1,4
c        DO L=0,LMO-1
c        DO J=1,JM-1
c          IF (SFS(J,L,KB).ne.UNDEF) SFS(J,L,KB) = SFS(J,L,KB) -
c     *         35.*SFM(J,L,KB)
c        END DO
c        END DO
c        TITLE(1:8)=TRIM(BASIN(KB))
c        lname(1:31)=TITLE(1:31)
c        sname='salt_sf_'//trim(BASIN(KB))
c        units='Sv'
c        XJL(2:JM,1:LMO+1)=SFS(1:JM-1,0:LMO,KB)
c        IF (QDIAG) CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,2,LMO+1,XJL
c     *       ,ZOC1,"Latitude","Depth (m)")
c      END DO
C****
      RETURN
  920 FORMAT ('1',A80/)
  921 FORMAT ('0',A20 / '0',6X,23I5 / 8X,23('-----'))
  922 FORMAT (F6.0,2X,23I5)

      END SUBROUTINE OSFOUT
#endif

      SUBROUTINE STROUT
!@sum  STROUT prints out strait diagnostics for ocean
!@auth Gavin Schmidt/Gary Russell
      USE CONSTANT, only : undef,teeny
      USE MODEL_COM, only : xlabel,lrunid,jmon0,jyear0,idacc,jdate0
     *     ,amon0,jdate,amon,jyear
      USE OCEAN, only : im,jm,lmo,dts,ze
      USE STRAITS, only : nmst,wist,dist,lmst,name_st
      USE ODIAG
      IMPLICIT NONE
      REAL*8, DIMENSION(LMO,NMST) :: AS
      INTEGER I,J,K,L,NS,N,KB,IP1
     *     ,LMSTMIN(LMO),SUMORMN(KOLNST)
      REAL*8 MQ
      CHARACTER :: TITLE*80
C****
C**** Strait diagnostics
C****
      TITLE=""
      TITLE(51:80)=XLB
      SUMORMN(:)=1
      SUMORMN(LN_KVM)=2
      SUMORMN(LN_ICFL)=0
      DO N=1,KOLNST
        if (n.eq.LN_MFLX .or. n.eq.LN_GFLX.or. n.eq.LN_SFLX.or. n.eq
     *       .LN_KVM.or. n.eq.LN_ICFL) THEN
          AS = 0.
          TITLE(1:40) = TRIM(LNAME_OLNST(N))//' ('//
     &         TRIM(UNITS_OLNST(N))//')'
          DO NS=1,NMST
            DO L=1,LMST(NS)
              AS(L,NS) = OLNST(L,NS,N)*SCALE_OLNST(N)/IDACC(IA_OLNST(N))
            END DO
          END DO
          IF (N.eq.LN_ICFL) THEN
c            LMSTMIN(:)=1
c            CALL STABLE (LMSTMIN,name_st,AS,TITLE,SUMORMN(N))
          ELSE
            CALL STABLE (LMST,name_st,AS,TITLE,SUMORMN(N))
          END IF
        END IF
      END DO
C****
      RETURN
      END SUBROUTINE STROUT

      SUBROUTINE STABLE (LMST,STRAIT,AX,TITLE,SUMORMN)
!@sum STABLE produces a layer by strait table on the line printer.
C****
C**** Input:
!@var LMST = number of layers for each strait
!@var STRAIT = names of straits
!@var AX = two dimensional input array
!@var TITLE = title of run
!@var SUMORMN flag for whether sum or mean or neither ar printed
C****
      USE OCEAN, only : lmo
      USE STRAITS, only : nmst
      IMPLICIT NONE
      INTEGER, PARAMETER :: LMAXST=6
      INTEGER, INTENT(IN), DIMENSION(NMST) :: LMST
      INTEGER, INTENT(IN) :: SUMORMN
      REAL*8, DIMENSION(LMO,NMST), INTENT(IN) :: AX
      CHARACTER*20, INTENT(IN), DIMENSION(NMST) :: STRAIT
      CHARACTER*80, INTENT(IN) :: TITLE
      INTEGER, SAVE :: LINECT = 0.
      REAL*8 SUML
      INTEGER NS,L,LMAX
C****
C**** Produce leading title lines
C****
      LINECT = LINECT + 19
      IF(LINECT.LE.63)  THEN
        WRITE (6,900) TITLE
      ELSE
        LINECT = 8
        WRITE (6,901) TITLE
      ENDIF
      WRITE (6,902) ('------',L=1,LMAXST),(L,L=1,LMAXST),
     *              ('------',L=1,LMAXST)
C****
C**** Calculate sums and print the table
C****
      DO NS=1,NMST
        LMAX = LMST(NS)
        SUML = 0.
        DO L=1,LMAX
          SUML = SUML + AX(L,NS)
        END DO
        SELECT CASE (SUMORMN)
        CASE (0)                ! no sum or mean
          WRITE (6,'(1X,A20,9X,10I6)') STRAIT(NS),NINT(AX(1:LMAX,NS))
        CASE (1)
          WRITE (6,'(1X,A20,F8.1,1X,10I6)') STRAIT(NS),SUML
     *         ,NINT(AX(1:LMAX,NS))
        CASE (2)
          SUML=SUML/(LMAX-1)
          WRITE (6,'(1X,A20,F8.1,1X,10I6)') STRAIT(NS),SUML
     *         ,NINT(AX(1:LMAX,NS))
        END SELECT
      END DO
C****
      RETURN
  900 FORMAT ('0'//1X,A80)
  901 FORMAT ('1',A80)
  902 FORMAT ('0',29('-'),6A6 / ' Strait',15X,'Sum/Mean',6I6 /
     *        ' ',29('-'),6A6)
  904 FORMAT (1X,A20,F8.1,1X,10I6)
      END SUBROUTINE STABLE

      SUBROUTINE STRMJL (OIJL,FAC,OLNST,FACST,SF)
!@sum  STRMJL calculates the latitude by layer stream function
!@auth G. Russell/G. Schmidt
C****
C**** Input:
!@var OIJL  = west-east and south-north tracer fluxes (kg/s)
!@var FAC   = global scaling factor
!@var OLNST = strait mass flux (kg/s)
!@var FACST = global scaling factor for straits
C**** Output:
!@var    SF = stream function (kg/s)
C****
      USE CONSTANT, only : undef
      USE OCEAN, only : im,jm,lmo
      USE STRAITS, only : nmst
      USE ODIAG, only : kbasin
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : pack_dataj,am_i_root
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: FAC,FACST
      REAL*8, INTENT(IN) ::
     &     OIJL(IM,grid%j_strt_halo:grid%j_stop_halo,LMO,2)
      REAL*8, INTENT(IN), DIMENSION(LMO,NMST) :: OLNST
      REAL*8, INTENT(OUT), DIMENSION(JM,0:LMO,4) :: SF
      REAL*8, DIMENSION(4) :: SUMB
      INTEGER :: KB,I,J,L,K
      REAL*8 SFx(grid%j_strt_halo:grid%j_stop_halo,0:LMO,4)
      REAL*8 :: SF_str(JM,0:LMO,4)
C****
C**** Zero out the Stream Function at the ocean bottom
C****
      SFX(:,:,:) = 0.
C****
C**** Integrate the Stream Function upwards from the ocean bottom
C****
      DO L=LMO-1,0,-1
      DO J=grid%j_strt,min(grid%j_stop,jm-1)
        DO KB=1,4
          SFX(J,L,KB) = SFX(J,L+1,KB)
        END DO
        DO I=1,IM
          KB = KBASIN(I,J+1)
          IF(KB.gt.0) SFX(J,L,KB) = SFX(J,L,KB) + OIJL(I,J,L+1,2)*FAC
        END DO
      END DO
      END DO

      call pack_dataj(grid,sfx,sf)

      if(am_i_root()) then
C****
C**** Add strait flow to the Stream Function
C****
      sf_str = 0.
      DO L=LMO-1,0,-1
        CALL STRMJL_STRAITS(L,SF_str,OLNST,FACST)
        sf_str(:,l,:) = sf_str(:,l,:) + sf_str(:,l+1,:)
      END DO
      sf = sf + sf_str

C****
C**** Calculate global Stream Function by summing it over 3 oceans
C****
      SF(:,:,4) = SF(:,:,1) + SF(:,:,2) + SF(:,:,3) + SF(:,:,4)
C****
C**** Mask streamfunction so that topography is put to skip value
C****
      DO J=1,JM-1
        SUMB = 0
        DO I=1,IM
          IF (KBASIN(I,J).gt.0) SUMB(KBASIN(I,J))=1.
        END DO
        SUMB(4)=SUMB(1)+SUMB(2)+SUMB(3)+SUMB(4)
        DO K=1,4
          IF (SUMB(K).eq.0) THEN
            SF(J,0:LMO,K) = UNDEF   ! UNDEF areas of no integral
          ELSE
            L=0
  620       L=L+1       ! UNDEF below bottom values
            IF (L.LE.LMO) THEN
              IF (SF(J,L,K).eq.0 .and. SF(J,L-1,K).eq.0) THEN
                SF(J,L:LMO,K) = UNDEF
              ELSE
                GOTO 620
              END IF
            END IF
          END IF
        END DO
      END DO
      endif ! am_i_root
C****
      RETURN
      END

      SUBROUTINE OBASIN
!@sum  OBASIN Read in KBASIN: 0=continent,1=Atlantic,2=Pacific,3=Indian
!@auth G. Russell
      USE OCEAN, only : IM,JM,focean
      USE ODIAG, only : kbasin
      USE FILEMANAGER
      IMPLICIT NONE
      CHARACTER TITLE*72, CBASIN(IM,JM)
      CHARACTER*6 :: FILEIN="KBASIN"
      INTEGER J,I,iu_KB,I72
C****
C**** read in basin data
      call openunit(FILEIN,iu_KB,.false.,.true.)

      READ  (iu_KB,900) TITLE
      WRITE (6,*) 'Read on unit ',iu_KB,': ',TITLE
      READ  (iu_KB,900)
      DO I72=1,1+(IM-1)/72
        DO J=JM,1,-1
          READ (iu_KB,901) (CBASIN(I,J),I=72*(I72-1)+1,MIN(IM,I72*72))
        END DO
      END DO
      call closeunit(iu_KB)

      DO J=1,JM
      DO I=1,IM
        SELECT CASE (CBASIN(I,J))
        CASE DEFAULT
          KBASIN(I,J) = 0
          IF (FOCEAN(I,J).gt.0) WRITE(6,*)
     *         "Warning: Ocean box not defined in KBASIN ",i,j
        CASE ('A')
          KBASIN(I,J) = 1
        CASE ('P')
          KBASIN(I,J) = 2
        CASE ('I')
          KBASIN(I,J) = 3
        CASE ('G')
          KBASIN(I,J) = 4
        END SELECT
      END DO
      END DO
C****
      RETURN
  900 FORMAT (A72)
  901 FORMAT (72A1)
      END SUBROUTINE OBASIN

      SUBROUTINE STRMIJ (MFU,FAC,OLNST,FACST,SF)
!@sum STRMIJ calculates the latitude by longitude stream function for a
!@+   given quantity.
C****
C**** Input:
!@var MFU  = west-east and south-north tracer fluxes (kg/s)
!@var   FAC = global scaling factor
!@var OLNST = strait mass flux (kg/s)
!@var FACST = global scaling factor for straits
C**** Output:
!@var    SF = stream function (kg/s)
C****
      Use OCEAN, Only: IM,JM,LMO, oDLAT_DG,oFJEQ=>FJEQ
      USE STRAITS, only : nmst,lmst
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: FAC,FACST
      REAL*8, INTENT(IN), DIMENSION(IM,JM) :: MFU
      REAL*8, INTENT(IN), DIMENSION(LMO,NMST) :: OLNST
      REAL*8, INTENT(OUT), DIMENSION(IM,JM) :: SF
      Integer*4 :: I,J,L, NSUM
      REAL*8 TSUM
C****
C**** Integrate up from South Pole
C****
      SF=0
      DO J=2,JM-1
        DO I=1,IM
          SF(I,J) = SF(I,J-1) + MFU(I,J)*FAC
        END DO
C****
C**** Add strait flow from to the Stream Function
C****
        CALL STRMIJ_STRAITS(J,SF,OLNST,FACST)
      EndDo  !  End of Do J=2,JM-1

C****
C**** Recalibrate SF to be 0 over middle of North America
C**** Include UO cells whose centers reside in 110-90 W, 32-40 N
C****
      NSUM = 0
      TSUM = 0
      Do J = Ceiling(oFJEQ+32/oDLAT_DG), Floor(oFJEQ+40/oDLAT_DG)
      Do I = Ceiling(70*IM/360.), Floor(90*IM/360.)
        TSUM = TSUM + SF(I,J)
        NSUM = NSUM + 1  ;  EndDo  ;  EndDo
      SF(:,:) = SF(:,:) - TSUM/NSUM
      RETURN
      END

#ifndef NEW_IO
      SUBROUTINE OJLOUT
!@sum OJLOUT calculates basin means, lat. and long. sections
!@+   and advective tracer fluxes
!@auth Gavin Schmidt/Gary Russell
      USE CONSTANT, only : undef,teeny
      USE MODEL_COM, only : idacc
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : ntm,trw0,trname,ntrocn,n_water
#endif
      USE OCEAN, only : im,jm,lmo,ze,imaxj,focean,dypo,dts,dxvo
     *     ,dxypo, oDLAT_DG, oDLON_DG
      USE DIAG_COM, only : qdiag
      USE MDIAG_COM, only :
     &     sname_strlen,units_strlen,lname_strlen
      USE ODIAG
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : to_per_mil
#endif
      IMPLICIT NONE
      REAL*8, DIMENSION(JM+3,LMO+1) :: XJL
      REAL*8 XB0(JM,LMO,NBAS),X0(IM,LMO+1),XS(IM,LMO+1),
     *     XBG(JM,LMO,NBAS),XBS(JM,LMO,NBAS),XG(IM,LMO+1)
#ifdef TRACERS_OCEAN
      REAL*8 XBT(JM,LMO,NBAS,NTM),XT(IM,LMO+1,NTM),XBTW(JM,LMO,NBAS)
#endif
      CHARACTER TITLE*80,EW*1,NS*1
      CHARACTER(len=lname_strlen) :: lname
      CHARACTER(len=sname_strlen) :: sname
      CHARACTER(len=units_strlen) :: units
      CHARACTER LABI*16,LABJ*16
      character*50 :: unit_string
      INTEGER I,J,L,KB,ISEC,II,ILON,JLAT,N,I1
      REAL*8 GOS,SOS,TEMGS,ASUM(IM),GSUM,ZONAL(LMO)

      IF (QDIAG) then
      XJL=0.

      XB0 = OJL(:,:,:,JL_M)
      XBG = OJL(:,:,:,JL_PT)
      XBS = OJL(:,:,:,JL_S)

      TITLE(51:80) = XLB
C****
C**** Quantites that are calculated
C****
C**** Basin and Zonal averages of tracers
C****
#ifdef TRACERS_OCEAN
      XBT = 0.
      DO L=1,LMO
        DO J= 1,JM
          DO I=1,IMAXJ(J)
            KB=KBASIN(I,J)
            IF(FOCEAN(I,J).gt..5)  THEN
              XBT(J,L,KB,:) = XBT(J,L,KB,:) + TOIJL(I,J,L,TOIJL_CONC,:)
            END IF
          END DO
        END DO
      END DO
      XBT(:,:,4,:)= SUM(XBT(:,:,1:4,:),DIM=3)
      if (n_water.ne.0) XBTW(:,:,:) =  XBT(:,:,:,n_water)
#endif

      DO KB=1,4
        DO J=1,JM
          DO L=1,LMO
            IF (XB0(J,L,KB).ne.0) THEN
#ifdef TRACERS_OCEAN
              do n=1,ntm
              if (to_per_mil(n).gt.0) then
                if (XBTW(j,l,kb).gt.0) then
                  XBT(j,l,kb,n)= 1d3*(XBT(j,l,kb,n)/
     *                 (XBTW(j,l,kb)*trw0(n))-1.)
                else
                  XBT(j,l,kb,n)=undef
                end if
c                XBT(j,l,kb,n)= 1d3*(XBT(j,l,kb,n)/
c     *               ((XB0(J,L,KB)-XBS(J,L,KB))*trw0(n))-1.)
              else
                XBT(j,l,kb,n)= 10.**(-ntrocn(n))*XBT(j,l,kb,n)/
     *               (XB0(J,L,KB))
              end if
              end do
#endif
              XBG(J,L,KB) = XBG(J,L,KB)/XB0(J,L,KB)
              XBS(J,L,KB) = XBS(J,L,KB)/XB0(J,L,KB)
            ELSE
              XBG(J,L,KB)=undef
              XBS(J,L,KB)=undef
#ifdef TRACERS_OCEAN
              XBT(j,l,kb,:)= undef
#endif
            END IF
          END DO
        END DO
        if (KB.le.3) then
          TITLE(21:50)=" "//TRIM(BASIN(KB))//" Basin Average"
          SNAME(5:30)="_"//BASIN(KB)(1:3)
        else
          TITLE(21:50)=" Zonal Average"
          SNAME(5:30)="_Zonal"
        end if
        LNAME(21:50)=TITLE(21:50)
        TITLE(1:20)="Temperature (C)"
        LNAME(1:20)="Temperature"
        SNAME(1:4)="Temp"
        UNITS="C"
        XJL(1:JM,1:LMO)=XBG(1:JM,1:LMO,KB)
        CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,1,LMO,XJL,ZOC
     *       ,"Latitude","Depth (m)")
        TITLE(1:20)="Salinity  (ppt)"
        LNAME(1:20)="Salinity"
        SNAME(1:4)="Salt"
        UNITS="ppt"
        XJL(1:JM,1:LMO)=XBS(1:JM,1:LMO,KB)
        CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,1,LMO,XJL,ZOC
     *       ,"Latitude","Depth (m)")
#ifdef TRACERS_OCEAN
        DO N=1,NTM
          if (to_per_mil(n).gt.0) then
            UNITS="permil"
          else
            UNITS=unit_string(ntrocn(n),'kg/kg')
          end if
          TITLE(1:20)=trim(trname(n))//" ("//trim(UNITS)//")"
          LNAME(1:20)=trim(trname(n))
          SNAME(1:4)=trname(n)(1:4)
          XJL(1:JM,1:LMO)=XBT(1:JM,1:LMO,KB,N)
          CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,1,LMO,XJL,ZOC
     *         ,"Latitude","Depth (m)")
        END DO
#endif
      END DO
C****
C**** Longitudinal sections
C****
      DO ISEC=1,NSEC
        XB0 = undef  ; XBS = undef ;  XBG = undef
#ifdef TRACERS_OCEAN
        XBT = undef
#endif
        IF (SEC_LON(ISEC).ne.-1) THEN
          I=nint( .5*im + SEC_LON(ISEC)/odlon_dg)
          IF (I.le.IM/2) THEN
            EW="W"
          ELSE
            EW="E"
          END IF
          DO L=1,LMO
            DO J= 1,JM
              IF (OIJL(I,J,L,IJL_MO).ne.0) THEN
                GOS = OIJL(I,J,L,IJL_G0M)/(OIJL(I,J,L,IJL_MO)*DXYPO(J))
                SOS = OIJL(I,J,L,IJL_S0M)/(OIJL(I,J,L,IJL_MO)*DXYPO(J))
                XBG(J,L,1)= TEMGS(GOS,SOS)
                XBS(J,L,1)= 1d3*SOS
#ifdef TRACERS_OCEAN
                do n=1,ntm
                  if (to_per_mil(n).gt.0) then
                    XBT(j,l,1,n)=1d3*(TOIJL(I,J,L,TOIJL_CONC,n)/
     *              (TOIJL(I,J,L,TOIJL_CONC,n_water)*trw0(n))-1.)
c                    XBT(j,l,1,n)=1d3*(TOIJL(I,J,L,TOIJL_CONC,n)/
c     *              ((OIJL(I,J,L,IJL_MO)*DXYPO(J)-OIJL(I,J,L,IJL_S0M))
c     *                   *trw0(n))-1.)
                  else
                    XBT(j,l,1,n)= 10.**(-ntrocn(n))*TOIJL(I,J,L
     *                   ,TOIJL_CONC,n)/(OIJL(I,J,L,IJL_MO)*DXYPO(J))
                  end if
                end do
#endif
              ELSE
                XBG(J,L,1)=undef ; XBS(J,L,1)=undef
#ifdef TRACERS_OCEAN
                XBT(J,L,1,:)=undef
#endif
              END IF
            END DO
          END DO
          ILON=ABS(SEC_LON(ISEC)-odlon_dg*0.5)
          WRITE(LABI,'(I3,A1)') ILON,EW
          TITLE(1:50)="Temperature Section         (C)"
          WRITE(TITLE(21:25),'(I3,A1)') ILON,EW
          LNAME(1:25)=TITLE(1:25)
          SNAME="Temp_"//adjustl(labi)
          UNITS="C"
          XJL(1:JM,1:LMO)=XBG(1:JM,1:LMO,1)
          CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,1,LMO,XJL,ZOC
     *         ,"Latitude","Depth (m)")
          TITLE(1:50)="Salinity Section            (ppt)"
          WRITE(TITLE(21:25),'(I3,A1)') ILON,EW
          LNAME(1:25)=TITLE(1:25)
          SNAME="Salt_"//adjustl(labi)
          UNITS="ppt"
          XJL(1:JM,1:LMO)=XBS(1:JM,1:LMO,1)
          CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,1,LMO,XJL,ZOC
     *         ,"Latitude","Depth (m)")
#ifdef TRACERS_OCEAN
          DO N=1,NTM
          if (to_per_mil(n).gt.0) then
            UNITS="permil"
          else
            UNITS=unit_string(ntrocn(n),'kg/kg')
          end if
          TITLE(1:50)=trname(n)//" Section          ("//TRIM(UNITS)//")"
          WRITE(TITLE(21:25),'(I3,A1)') ILON,EW
          LNAME(1:25)=TITLE(1:25)
          SNAME=trim(trname(n))//"_"//adjustl(labi)
          XJL(1:JM,1:LMO)=XBT(1:JM,1:LMO,1,N)
          CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,1,LMO,XJL,ZOC
     *         ,"Latitude","Depth (m)")
          END DO
#endif

C**** EW fluxes of mass, heat and salt
C**** Note that ocean fluxes are defined 1:JM-1 and need to be moved up
C**** one place for compatibility with AGCM output format routines
          DO L=1,LMO
            DO J= 1,JM-1
C**** GM fluxes are also saved, so add the GM heat and salt fluxes here
              IF (OIJL(I,J,L,IJL_MFU).ne.0) THEN
                XB0(J,L,1)= 1d0*OIJL(I,J,L,IJL_MFU)/(IDACC(1)*DTS
     *               *DYPO(J) *(ZE(L)-ZE(L-1)))
                XBG(J,L,1)= 1d-6*(OIJL(I,J,L,IJL_GFLX)+OIJL(I,J,L
     *               ,IJL_GGMFL))/(IDACC(1)*DTS*DYPO(J)*(ZE(L)-ZE(L-1)))
                XBS(J,L,1)= 1d-6*(OIJL(I,J,L,IJL_SFLX)+OIJL(I,J,L
     *               ,IJL_SGMFL))/(IDACC(1)*DTS*DYPO(J)*(ZE(L)-ZE(L-1)))
              ELSE
                XB0(J,L,1)=undef
                XBG(J,L,1)=undef
                XBS(J,L,1)=undef
              END IF
            END DO
          END DO
          ILON=ABS(SEC_LON(ISEC))
          WRITE(LABI,'(I3,A1)') ILON,EW
          TITLE(1:50)=" EW Mass Flux            (kg/s m^2)"
          WRITE(TITLE(21:25),'(I3,A1)') ILON,EW
          LNAME(1:25)=TITLE(1:25)
          SNAME="EWmflx_"//adjustl(labi)
          UNITS="kg/s m^2"
          XJL(2:JM,1:LMO)=XB0(1:JM-1,1:LMO,1)
          CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,2,LMO,XJL,ZOC
     *         ,"Latitude","Depth (m)")
          TITLE(1:50)=" EW Heat Flux            (10^6 J/s m^2)"
          WRITE(TITLE(21:25),'(I3,A1)') ILON,EW
          LNAME(1:25)=TITLE(1:25)
          SNAME="EWhflx_"//adjustl(labi)
          UNITS="10^6 J/s m^2"
          XJL(2:JM,1:LMO)=XBG(1:JM-1,1:LMO,1)
          CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,2,LMO,XJL,ZOC
     *         ,"Latitude","Depth (m)")
          TITLE(1:50)=" EW Salt Flux            (10^6 kg/s m^2)"
          WRITE(TITLE(21:25),'(I3,A1)') ILON,EW
          LNAME(1:25)=TITLE(1:25)
          SNAME="EWsflx_"//adjustl(labi)
          UNITS="10^6 kg/s m^2"
          XJL(2:JM,1:LMO)=XBS(1:JM-1,1:LMO,1)
          CALL POUT_JL(TITLE,LNAME,SNAME,UNITS,2,LMO,XJL,ZOC
     *         ,"Latitude","Depth (m)")
        END IF
      END DO
C****
C**** Latitudinal sections
C****
C**** Define starting point for wrap (~20W) (to avoid splitting pacific)
      I1 = 1+200./odlon_dg

      ASUM=0. ; GSUM=0. ; ZONAL=0.
      DO ISEC=1,NSEC
        X0=undef ; XS=undef ; XG=undef
#ifdef TRACERS_OCEAN
        XT=undef
#endif
        IF (SEC_LAT(ISEC).ne.-1) THEN
          IF (SEC_LAT(ISEC).le.0) THEN
            NS="S"
          ELSE
            NS="N"
          END IF
          J=nint(.5*jm+SEC_LAT(ISEC)/odlat_dg )
          DO L=1,LMO
            DO II= 1,IM
              I = II+I1-1
              IF (I.gt.IM) I=I-IM
              IF (OIJL(I,J,L,IJL_MO).ne.0) THEN
                GOS = OIJL(I,J,L,IJL_G0M)/(OIJL(I,J,L,IJL_MO)*DXYPO(J))
                SOS = OIJL(I,J,L,IJL_S0M)/(OIJL(I,J,L,IJL_MO)*DXYPO(J))
                XG(II,L)= TEMGS(GOS,SOS)
                XS(II,L)= 1d3*SOS
#ifdef TRACERS_OCEAN
                do n=1,ntm
                  if (to_per_mil(n).gt.0.and.TOIJL(I,J,L,TOIJL_CONC
     *                 ,n_water).gt.0) then
                    XT(II,l,n)= 1d3*(TOIJL(I,J,L,TOIJL_CONC,n)/
     *              (TOIJL(I,J,L,TOIJL_CONC,n_water)*trw0(n))-1.)
c                    XT(II,l,n)= 1d3*(TOIJL(I,J,L,TOIJL_CONC,n)/
c     *              ((OIJL(I,J,L,IJL_MO)*DXYPO(J)-OIJL(I,J,L,IJL_S0M))
c     *                   *trw0(n))-1.)
                  else
                    XT(II,l,n)= 10.**(-ntrocn(n))*TOIJL(I,J,L,TOIJL_CONC
     *                   ,n)/(OIJL(I,J,L,IJL_MO)*DXYPO(J))
                  end if
                end do
#endif
              ELSE
                XG(II,L)=undef ; XS(II,L)=undef
#ifdef TRACERS_OCEAN
                XT(II,L,:)=undef
#endif
              END IF
            END DO
          END DO
          JLAT=ABS(SEC_LAT(ISEC)-0.5*odlat_dg)
          WRITE(LABJ,'(I3,A1)') JLAT,NS
          TITLE(1:50)="Temperature Section              (C)"
          WRITE(TITLE(21:25),'(I3,A1)') JLAT,NS
          LNAME(1:25)=TITLE(1:25)
          SNAME="Temp_"//adjustl(labj)
          UNITS="C"
          CALL POUT_IL(TITLE,sname,lname,units,I1,1,LMO,XG
     *         ,ZOC,"Longitude","Depth (m)",ASUM,GSUM,ZONAL)
          TITLE(1:50)="Salinity Section                 (ppt)"
          WRITE(TITLE(21:25),'(I3,A1)') JLAT,NS
          LNAME(1:25)=TITLE(1:25)
          SNAME="Salt_"//adjustl(labj)
          UNITS="ppt"
          CALL POUT_IL(TITLE,sname,lname,units,I1,1,LMO,XS
     *         ,ZOC,"Longitude","Depth (m)",ASUM,GSUM,ZONAL)
#ifdef TRACERS_OCEAN
          DO N=1,NTM
          if (to_per_mil(n).gt.0) then
            UNITS="permil"
          else
            UNITS=unit_string(ntrocn(n),'kg/kg')
          end if
          TITLE(1:50)=trname(n)//" Section          ("//TRIM(UNITS)//")"
          WRITE(TITLE(21:25),'(I3,A1)') JLAT,NS
          LNAME(1:25)=TITLE(1:25)
          SNAME=trim(trname(n))//"_"//adjustl(labj)
          CALL POUT_IL(TITLE,sname,lname,units,I1,1,LMO,XT(1,1,N)
     *         ,ZOC,"Longitude","Depth (m)",ASUM,GSUM,ZONAL)
          END DO
#endif
C**** Fluxes
          DO L=1,LMO
            DO II= 1,IM
              I = II+I1-1
              IF (I.gt.IM) I=I-IM
              IF (OIJL(I,J,L,IJL_MFV).ne.0) THEN
                X0(II,L)= OIJL(I,J,L,IJL_MFV)/(IDACC(1)*DTS
     *               *DXVO(J) *(ZE(L)-ZE(L-1)))
                XG(II,L)= 1d-6*(OIJL(I,J,L,IJL_GFLX+1)+OIJL(I,J,L
     *               ,IJL_GGMFL+1))/(IDACC(1)*DTS*DXVO(J)*(ZE(L)-ZE(L-1)
     *               ))
                XS(II,L)= 1d-6*(OIJL(I,J,L,IJL_SFLX+1)+OIJL(I,J,L
     *               ,IJL_SGMFL+1))/(IDACC(1)*DTS*DXVO(J)*(ZE(L)-ZE(L-1)
     *               ))
              ELSE
                X0(II,L)=undef ; XG(II,L)=undef ; XS(II,L)=undef
              END IF
            END DO
          END DO
          JLAT=ABS(SEC_LAT(ISEC))
          WRITE(LABJ,'(I3,A1)') JLAT,NS
          TITLE(1:50)=" NS Mass Flux            (kg/s m^2)"
          WRITE(TITLE(21:25),'(I3,A1)') JLAT,NS
          LNAME(1:25)=TITLE(1:25)
          SNAME="NSmflx_"//adjustl(labj)
          UNITS="kg/s m^2"
          CALL POUT_IL(TITLE,sname,lname,units,I1,2,LMO,X0
     *         ,ZOC,"Longitude","Depth (m)",ASUM,GSUM,ZONAL)
          TITLE(1:50)=" NS Heat Flux            (10^6  J/s m^2)"
          WRITE(TITLE(21:25),'(I3,A1)') JLAT,NS
          LNAME(1:25)=TITLE(1:25)
          SNAME="NShflx_"//adjustl(labj)
          UNITS="10^6 W/m^2"
          CALL POUT_IL(TITLE,sname,lname,units,I1,2,LMO,XG
     *         ,ZOC,"Longitude","Depth (m)",ASUM,GSUM,ZONAL)
          TITLE(1:50)=" NS Salt Flux            (10^6 kg/s m^2)"
          WRITE(TITLE(21:25),'(I3,A1)') JLAT,NS
          LNAME(1:25)=TITLE(1:25)
          SNAME="NSsflx_"//adjustl(labj)
          UNITS="10^6 kg/s m^2"
          CALL POUT_IL(TITLE,sname,lname,units,I1,2,LMO,XS
     *         ,ZOC,"Longitude","Depth (m)",ASUM,GSUM,ZONAL)
        END IF
      END DO
      END IF
C****
      RETURN
      END SUBROUTINE OJLOUT
#endif

      subroutine oijl_prep
c
c Convert oijl accumulations into the desired units
c
      use ocean, only : im,jm,lmo,lmm,imaxj,focean,dxypo,dxvo,dypo,dts
      use odiag, only : koijl,oijl_out,oijl=>oijl_loc
     &     ,ijl_mo,ijl_mou,ijl_mov,ijl_g0m,ijl_s0m,ijl_ptm,ijl_pdm
     &     ,ijl_mfu,ijl_mfv,ijl_mfw,ijl_mfw2,ijl_ggmfl,ijl_sgmfl
     &     ,ijl_wgfl,ijl_wsfl,ijl_kvm,ijl_kvg,ijl_gflx,ijl_sflx
     &     ,oij=>oij_loc,ij_sf,olnst,ln_mflx
#ifdef OCN_Mesoscales
     &     ,ijl_ueddy,ijl_veddy,ijl_n2
#endif
#ifdef TRACERS_OCEAN
      use odiag, only :
     &     ktoijlx,toijl_out,divbya_toijl,kn_toijl,toijl_loc,toijl_conc
      USE OCN_TRACER_COM, only : trw0,n_Water,to_per_mil
#endif
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : am_i_root,halo_update,south
     &     ,pack_data,unpack_data ! for horz stream function
      implicit none
      integer i,j,l,k,kk,n
      real*8 mass,gos,sos,temgs,volgs,fac,facst
      integer :: j_0,j_1,j_0s,j_1s
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) :: mfu
      real*8, dimension(:,:), allocatable :: mfu_glob,sf_glob

      j_0 = grid%j_strt
      j_1 = grid%j_stop
      j_0s = grid%j_strt_skp
      j_1s = grid%j_stop_skp

      oijl_out(:,:,:,:) = 0.

c
c Cell-centered quantities. Some conversions to per square meter
c
      do l=1,lmo
      do j=j_0,j_1
      do i=1,imaxj(j)
        mass = oijl(i,j,l,ijl_mo)*dxypo(j)
        if(focean(i,j).le..5 .or. mass.le.0.) cycle
        oijl_out(i,j,l,ijl_mo) = mass
        oijl_out(i,j,l,ijl_g0m) = oijl(i,j,l,ijl_g0m)
        oijl_out(i,j,l,ijl_s0m) = oijl(i,j,l,ijl_s0m)

#ifdef OCN_Mesoscales
        oijl_out(i,j,l,ijl_n2) = oijl(i,j,l,ijl_n2)
        oijl_out(i,j,l,ijl_ueddy) = oijl(i,j,l,ijl_ueddy)
        oijl_out(i,j,l,ijl_veddy) = oijl(i,j,l,ijl_veddy)
#endif

c
c compute potential temperature and potential density
c
        gos = oijl(i,j,l,ijl_g0m) / mass
        sos = oijl(i,j,l,ijl_s0m) / mass
        oijl_out(i,j,l,ijl_ptm) = mass*temgs(gos,sos)
        oijl_out(i,j,l,ijl_pdm) = mass*(1d0/volgs(gos,sos)-1000d0)
      enddo
      enddo
      enddo

c
c Vertical fluxes.  Some conversions to per square meter
c
      do l=1,lmo-1
      do j=j_0,j_1
      do i=1,imaxj(j)
        oijl_out(i,j,l,ijl_mfw) = oijl(i,j,l,ijl_mfw)/dxypo(j)
        oijl_out(i,j,l,ijl_mfw2) = oijl(i,j,l,ijl_mfw2)/(dxypo(j)**2)
        oijl_out(i,j,l,ijl_ggmfl+2) = oijl(i,j,l,ijl_ggmfl+2)/dxypo(j)
        oijl_out(i,j,l,ijl_sgmfl+2) = oijl(i,j,l,ijl_sgmfl+2)/dxypo(j)
        oijl_out(i,j,l,ijl_wgfl) = oijl(i,j,l,ijl_wgfl)/dxypo(j)
        oijl_out(i,j,l,ijl_wsfl) = oijl(i,j,l,ijl_wsfl)/dxypo(j)
        oijl_out(i,j,l,ijl_kvm) = oijl(i,j,l,ijl_kvm)
        oijl_out(i,j,l,ijl_kvg) = oijl(i,j,l,ijl_kvg)
        oijl_out(i,j,l,ijl_gflx+2) = oijl(i,j,l,ijl_gflx+2)/dxypo(j)
        oijl_out(i,j,l,ijl_sflx+2) = oijl(i,j,l,ijl_sflx+2)/dxypo(j)
      enddo
      enddo
      enddo

c
c Horizontal fluxes.  Some conversions to per meter
c
      do k=1,koijl
        call halo_update(grid,oijl(:,:,:,k),from=south)
      enddo
      do l=1,lmo
        do j=j_0s,j_1s
        do i=1,im
          oijl_out(i,j,l,ijl_mfu) = oijl(i,j,l,ijl_mfu)
          oijl_out(i,j,l,ijl_gflx) = oijl(i,j,l,ijl_gflx)
          oijl_out(i,j,l,ijl_sflx) = oijl(i,j,l,ijl_sflx)
          oijl_out(i,j,l,ijl_ggmfl) = oijl(i,j,l,ijl_ggmfl)
          oijl_out(i,j,l,ijl_sgmfl) = oijl(i,j,l,ijl_sgmfl)
        enddo
        do i=1,im-1
          oijl_out(i,j,l,ijl_mou) =
     &         .5*(oijl(i,j,l,ijl_mo)+oijl(i+1,j,l,ijl_mo))*dypo(j)
        enddo
        i=im
          oijl_out(i,j,l,ijl_mou) =
     &         .5*(oijl(i,j,l,ijl_mo)+oijl(1,j,l,ijl_mo))*dypo(j)
        enddo ! j
        do j=max(2,j_0),j_1
        do i=1,im
          oijl_out(i,j,l,ijl_mfv) = oijl(i,j-1,l,ijl_mfv)
          oijl_out(i,j,l,ijl_mov) =
     &         .5*(oijl(i,j,l,ijl_mo)+oijl(i,j-1,l,ijl_mo))*dxvo(j-1)
          oijl_out(i,j,l,ijl_gflx+1) = oijl(i,j-1,l,ijl_gflx+1)
          oijl_out(i,j,l,ijl_sflx+1) = oijl(i,j-1,l,ijl_sflx+1)
          oijl_out(i,j,l,ijl_ggmfl+1) = oijl(i,j-1,l,ijl_ggmfl+1)
          oijl_out(i,j,l,ijl_sgmfl+1) = oijl(i,j-1,l,ijl_sgmfl+1)
        enddo
        enddo ! j
      enddo

C****
C**** Calculate Horizontal Mass Stream Function
C****
      do j=j_0s,j_1s
      do i=1,im
        mfu(i,j) = 0.
        do l=1,lmm(i,j)
          mfu(i,j) = mfu(i,j) + oijl(i,j,l,ijl_mfu)
        enddo
      enddo
      enddo
      if(am_i_root()) then
        allocate(mfu_glob(im,jm),sf_glob(im,jm))
      else
        allocate(mfu_glob(1,1),sf_glob(1,1))
      endif
      call pack_data(grid,mfu,mfu_glob)
      if(am_i_root()) then
        FAC   = -1d-9/dts
        FACST = -1d-9/dts
        CALL STRMIJ(MFU_GLOB,FAC,OLNST(1,1,LN_MFLX),FACST,SF_GLOB)
      endif
      call unpack_data(grid,sf_glob,oij(:,:,ij_sf))
      deallocate(mfu_glob,sf_glob)

#ifdef TRACERS_OCEAN
C****
C**** Tracers
C****
      toijl_out(:,:,:,:) = 0.
      kk = 1
      toijl_out(:,:,:,kk) = oijl_out(:,:,:,ijl_mo)
      do kk=2,ktoijlx
        k = kn_toijl(1,kk)
        n = kn_toijl(2,kk)
        if(k.le.0 .or. n.le.0) cycle
        toijl_out(:,:,:,kk) = toijl_loc(:,:,:,k,n)
        if(divbya_toijl(kk)) then
          do l=1,lmo; do j=j_0,j_1
            toijl_out(:,j,l,kk) = toijl_out(:,j,l,kk)/dxypo(j)
          enddo; enddo
        endif
        if(to_per_mil(n)>0 .and. n.ne.n_Water) then
          toijl_out(:,:,:,kk) = 1d3*(toijl_out(:,:,:,kk)/trw0(n)
     &         -toijl_loc(:,:,:,TOIJL_conc,n_water))
        endif
      enddo
#endif

      return
      end subroutine oijl_prep

      subroutine basin_prep
c
c Calculate zonal sums for ocean basins
c
      use ocean, only : im,jm,lmo,imaxj,focean,dxypo
      USE OCEAN, only : dts
      USE STRAITS, only : nmst
      use odiag, only : oijl=>oijl_loc,kbasin
     &     ,ijl_mo,ijl_g0m,ijl_s0m,ijl_mfu,ijl_mfv,ijl_gflx,ijl_sflx
     &     ,ijl_ggmfl,ijl_sgmfl
     &     ,ojl,ojl_out,nbas,nqty,jl_m,jl_pt,jl_s,jl_sf
     &     ,otj,otjcomp,otj_out
     &     ,olnst,ln_mflx,ln_gflx,ln_sflx
     &     ,sfm
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : pack_dataj,am_i_root
      implicit none
      integer i,j,l,kb
      real*8 gos,sos,temgs
      integer :: j_0,j_1, j_0h,j_1h
      real*8 :: ojlx(grid%j_strt_halo:grid%j_stop_halo,lmo,nbas,4)
      INTEGER K,KQ,NOIJL,NOIJLGM,NOLNST,NS,KK,KC
      REAL*8 SOIJL,SOIJLGM,FAC,FACST
      REAL*8 SCALEM(3),SCALES(3),SOLNST(NMST),MV(4),MT(4)
      INTEGER IS(4)
      REAL*8 ::
     &     X(grid%j_strt_halo:grid%j_stop_halo,4,3)
     &    ,XCOMP(grid%j_strt_halo:grid%j_stop_halo,4,3,3)

      j_0 = grid%j_strt
      j_1 = grid%j_stop
      j_0h = grid%j_strt_halo
      j_1h = grid%j_stop_halo

      ojlx = 0.
      do l=1,lmo
      do j=j_0,j_1
      do i=1,imaxj(j)
        if(focean(i,j).lt..5) cycle
        kb=kbasin(i,j)
        ojlx(j,l,kb,jl_m) = ojlx(j,l,kb,jl_m) +oijl(i,j,l,ijl_mo)
        ojlx(j,l,kb,jl_pt) = ojlx(j,l,kb,jl_pt) +oijl(i,j,l,ijl_g0m)
        ojlx(j,l,kb,jl_s) = ojlx(j,l,kb,jl_s) +oijl(i,j,l,ijl_s0m)
      enddo
      enddo
      enddo

      ojlx(:,:,4,:) = sum(ojlx(:,:,1:4,:),dim=3)

      do kb=1,4
      do l=1,lmo
      do j=j_0,j_1
        ojlx(j,l,kb,jl_m) = ojlx(j,l,kb,jl_m)*dxypo(j)
        if(ojlx(j,l,kb,jl_m).eq.0.) cycle
        gos = ojlx(j,l,kb,jl_pt)/ojlx(j,l,kb,jl_m)
        sos = ojlx(j,l,kb,jl_s )/ojlx(j,l,kb,jl_m)
        ojlx(j,l,kb,jl_pt) = temgs(gos,sos)*ojlx(j,l,kb,jl_m)
        ojlx(j,l,kb,jl_s) = ojlx(j,l,kb,jl_s)*1d3
      enddo
      enddo
      enddo

      call pack_dataj(grid,ojlx,ojl)

C****
C**** Calculate Mass Stream Function
C****
      FAC   = -1d-9/DTS
      FACST = -1d-9/DTS
      CALL STRMJL (OIJL(1,j_0h,1,IJL_MFU),FAC,
     &     OLNST(1,1,LN_MFLX),FACST,SFM)
      if(am_i_root()) then
        ojl(1,:,:,jl_sf) = 0.
        ojl(2:jm,:,:,jl_sf) = sfm(1:jm-1,1:lmo,:)
      endif

c copy to j,l,n array
      if(am_i_root()) ojl_out = reshape(ojl,shape(ojl_out))


c
c Transports
c

C****
C****  N          Contents of OIJL(I,J,L,N)       Scaling factor
C****  -          -------------------------       --------------
C**** IJL_MV      South-north mass flux (kg/s)    2/IDACC(1)*NDYNO
C**** IJL_GFLX+1  South-north heat flux (J/s)     1/IDACC(1)*DTS
C**** IJL_SFLX+1  South-north salt flux (kg/s)    1/IDACC(1)*DTS
C****
C****  N       Contents of OLNST(L,N)             Scaling factor
C****  -       ----------------------             --------------
C**** LN_MFLX  Strait mass flux (kg/s)            1/IDACC(1)*DTS
C**** LN_GFLX  Strait heat flux (J/s)            .5/IDACC(1)*DTS
C**** LN_SFLX  Strait salt flux (kg/s)           .5/IDACC(1)*DTS
C****
      SCALEM(1) = 1.D- 9/(DTS)
      SCALEM(2) = 1.D-15/(DTS)
      SCALEM(3) = 1.D- 6/(DTS)
      SCALES(1) = 1.D- 9/(DTS)
      SCALES(2) = .5D-15/(DTS)
      SCALES(3) = .5D- 6/(DTS)
C****
      X = 0. ; XCOMP =0.
C****
C**** Loop over quantites
C****
      DO KQ=1,3
C****
C**** Calculate the transport entering each domain from the
C**** southern boundary
C****
      SELECT CASE (KQ)
      CASE (1)
        NOIJL=IJL_MFV
        NOIJLGM=0
      CASE (2)
        NOIJL=IJL_GFLX+1
        NOIJLGM=IJL_GGMFL+1
      CASE (3)
        NOIJL=IJL_SFLX+1
        NOIJLGM=IJL_SGMFL+1
      END SELECT

      DO J=J_0,min(J_1,jm-1)
        DO I=1,IM
          IF(OIJL(I,J,1,IJL_MFV).eq.0) CYCLE
          KB = KBASIN(I,J+1)
          SOIJL = 0.
          SOIJLGM = 0.
          DO L=1,LMO
            SOIJL = SOIJL + OIJL(I,J,L,NOIJL)
            IF (NOIJLGM.gt.0) SOIJLGM = SOIJLGM + OIJL(I,J,L,NOIJLGM)
          END DO
          X(J,KB,KQ) = X(J,KB,KQ) + (SOIJL+SOIJLGM)*SCALEM(KQ)
          IF (KQ.ne.1) THEN
c        X(J,KB,5) = X(J,KB,5) + SOIJLGM*SCALEM(KQ)
c        X(J, 4,5) = X(J, 4,5) + SOIJLGM*SCALEM(KQ)
            XCOMP(J,KB,2,KQ) = XCOMP(J,KB,2,KQ) + SOIJLGM*SCALEM(KQ)
          END IF
        END DO
        X(J,4,KQ) = SUM(X(J,1:4,KQ))
        XCOMP(J,4,2,KQ) = SUM(XCOMP(J,1:4,2,KQ))
      END DO
C****
C**** Accumulate northward transports by overturning and by
C**** horizontal gyre
C**** 
      DO J=J_0,min(J_1,jm-1)
        DO L=1,LMO
          IS=0 ; MV=0 ; MT=0
          DO I=1,IM
            IF(OIJL(I,J,L,IJL_MFV).eq.0)  CYCLE
            KB = KBASIN(I,J+1)
            IS(KB) = IS(KB) + 1
            MV(KB) = MV(KB) + OIJL(I,J,L,IJL_MFV)
            IF (NOIJL.gt.0) MT(KB) = MT(KB) + OIJL(I,J,L,NOIJL)/OIJL(I,J
     *           ,L,IJL_MFV)
          END DO
          IS(4) = SUM(IS(1:4))
          MV(4) = SUM(MV(1:4))
          IF (NOIJL.gt.0) MT(4) = SUM(MT(1:4))
          DO KB=1,4
            IF(IS(KB).eq.0) CYCLE
c     X(J,KB,4) = X(J,KB,4) + SCALEM(2)*MV(KB)*MT(KB)/IS(KB)
            IF (KQ.eq.1) THEN
              XCOMP(J,KB,1,KQ) = XCOMP(J,KB,1,KQ) +SCALEM(KQ)*MV(KB)
            ELSE
              XCOMP(J,KB,1,KQ) = XCOMP(J,KB,1,KQ) +SCALEM(KQ)*MV(KB)
     *             *MT(KB)/IS(KB)
            END IF
          END DO
        END DO
        DO KB=1,4
c     X(J,KB,6) = X(J,KB,2) - X(J,KB,4) - X(J,KB,5)
          XCOMP(J,KB,3,KQ)=
     &         X(J,KB,KQ)-XCOMP(J,KB,1,KQ)-XCOMP(J,KB,2,KQ)
        END DO
      END DO

      END DO ! end loop over kq

      call pack_dataj(grid,x,otj(1:jm,:,:))
      call pack_dataj(grid,xcomp,otjcomp(1:jm,:,:,:))

      if(am_i_root()) then
      DO KQ=1,3
C****
C**** Calculate transport through straits from latitude to another
C**** within the same basin
C****
      SELECT CASE (KQ)
      CASE (1)
        NOLNST=LN_MFLX
      CASE (2)
        NOLNST=LN_GFLX
      CASE (3)
        NOLNST=LN_SFLX
      END SELECT

      DO NS=1,NMST
        SOLNST(NS) = 0.
        DO L=1,LMO
          SOLNST(NS) = SOLNST(NS) + OLNST(L,NS,NOLNST)
        END DO
      END DO

      CALL OTJ_STRAITS(OTJ,SOLNST,SCALES(KQ),KQ)

C**** Fill in south polar value
      OTJ(0,1:4,KQ)  = OTJ(1,1:4,KQ)

      END DO ! end loop over kq

C**** Replace salt by salt - .035 times mass
      DO KB=1,4
        DO J=0,JM
          OTJ(J,KB,3) = OTJ(J,KB,3) - .035d3*OTJ(J,KB,1)
          OTJCOMP(J,KB,:,3) = OTJCOMP(J,KB,:,3)-.035d3*OTJCOMP(J,KB,:,1)
        END DO
      END DO

C**** combine OTJ and OTJCOMP into a single output array
      otj_out(1,:) = 0.
      kk = 0
      do kq=1,3
      do kb=1,4
        kk = kk + 1
        otj_out(2:jm,kk) = otj(1:jm-1,kb,kq)
        do kc=1,3
          kk = kk + 1
          otj_out(2:jm,kk) = otjcomp(1:jm-1,kb,kc,kq)
        enddo
      enddo
      enddo

      endif ! am_i_root

      return
      end subroutine basin_prep

#ifndef NEW_IO
      SUBROUTINE OTJOUT
!@sum OTJOUT print vertically integrated basin and zonal
!@+   northward transports
!@auth Gavin Schmidt/Gary Russell
      USE MODEL_COM, only : idacc
      USE OCEAN, only : im,jm,lmo
      USE DIAG_COM, only : qdiag
      USE ODIAG
      IMPLICIT NONE
      CHARACTER TITLE*80, YAXIS(3)*8
      REAL*8 VLAT(0:JM)
      INTEGER, PARAMETER :: INC=1+(JM-1)/24
      CHARACTER*50, DIMENSION(11) ::  OTJNAME=(/
     *     'Northward Transport of Mass (10^9 kg/sec)        ',
     *     'North. Trans. of Potential Enthalpy (10^15 W)    ',
     *     'North. Trans. of Salt - .035*Mass (10^6 kg/sec)  ',
     *     'North. Trans. of Heat in Atl. Ocean (10^15 W)    ',
     *     'North. Trans. of Heat in Pac. Ocean (10^15 W)    ',
     *     'North. Trans. of Heat in Indian Ocean (10^15 W)  ',
     *     'North. Trans. of Heat in Global Ocean (10^15 W)  ',
     *     'North. Trans. of Salt in Atl. Ocean (10^6 kg/s)  ',
     *     'North. Trans. of Salt in Pac. Ocean (10^6 kg/s)  ',
     *     'North. Trans. of Salt in Indian Ocean (10^6 kg/s)',
     *     'North. Trans. of Salt in Global Ocean (10^6 kg/s)'/)
      DATA YAXIS /'  Mass  ', 'Enthalpy', '  Salt  '/
      INTEGER J,KQ,KB,IDLAT
      REAL*8 FJEQ
      REAL*8 X(0:JM,4,3),XCOMP(0:JM,4,3,3)

      X = OTJ/idacc(1)
      XCOMP = OTJCOMP/idacc(1)

      IDLAT= NINT(180./(JM-1))
      FJEQ = JM/2.
      DO J=1,JM-1
        VLAT( J) = IDLAT*(J-FJEQ)
      END DO
      VLAT( 0) = -90.
      VLAT(JM) =  90.
C****
C**** Write titles and data to disk.
C****
      TITLE(51:80)=XLB
      DO KQ=1,3
        TITLE(1:50)=OTJNAME(KQ)
        WRITE (6,907) TITLE(1:72)
C**** print out truncated series to PRT file
        WRITE (6,903) (NINT(VLAT(J)),J=JM,0,-INC)
        DO KB=1,4
          WRITE (6,906) BASIN(KB),(X(J,KB,KQ),J=JM,0,-INC)
        END DO
        WRITE(6,904)
        IF (QDIAG) THEN
          WRITE (iu_otj,*) TITLE
          WRITE (iu_otj,*) 'Latitude'
          WRITE (iu_otj,*) YAXIS(KQ)
          WRITE (iu_otj,*)
     *         ' Lat   Atlantic   Pacific    Indian     Global'
          WRITE (iu_otj,971) (VLAT(J),(X(J,KB,KQ),KB=1,4),J=0,JM)
          WRITE (iu_otj,*)
        END IF
      END DO
C****
C**** Write titles and data to disk for northward transports
C**** by components
C****
      DO KQ=2,3
        WRITE(6,908)
        DO KB=1,4
          TITLE(1:50)=OTJNAME(3+KB+4*(KQ-2))
          WRITE (6,907) TITLE(1:72)
C**** print out truncated series to PRT file
          WRITE(6,903) (NINT(VLAT(J)),J=JM,0,-INC)
          WRITE(6,906)BASIN(KB)(1:3)//" Overtrn",
     &         (XCOMP(J,KB,1,KQ),J=JM,0,-INC)
          WRITE(6,906)BASIN(KB)(1:3)//" GM flx ",
     &         (XCOMP(J,KB,2,KQ),J=JM,0,-INC)
          WRITE(6,906)BASIN(KB)(1:3)//" Hor gyr",
     &         (XCOMP(J,KB,3,KQ),J=JM,0,-INC)
          WRITE(6,906)BASIN(KB)(1:3)//" Total  ",
     &         (X(J,KB,KQ),J=JM,0,-INC)
          WRITE(6,904)
          IF (QDIAG) THEN
            WRITE (iu_otj,*) TITLE
            WRITE (iu_otj,*) 'Latitude'
            WRITE (iu_otj,*) YAXIS(2)
            WRITE (iu_otj,*)
     *           ' Lat      Total    Overturn    GM_flux    Hor_Gyre'
            WRITE (iu_otj,976)(VLAT(J),X(J,KB,KQ),
     &           XCOMP(J,KB,1:3,KQ),J=0,JM)
            WRITE (iu_otj,*)
          END IF
        END DO
      END DO
      RETURN
C****
  970 FORMAT (A50,A6,2X,2A4)
  971 FORMAT (F6.0,4F10.3)
  972 FORMAT (A,1P,3E15.4)
  973 FORMAT (A,2I6)
  976 FORMAT (F6.0,4F10.3)
  903 FORMAT (' ',131('-')/,' Latitude  ',24I5)
  904 FORMAT (' ',131('-'))
  906 FORMAT (' ',A11,24F5.1)
  907 FORMAT ('0',A)
  908 FORMAT ('1')

      END SUBROUTINE OTJOUT
#endif

      Subroutine STRMJL_STRAITS (L,SF,OLNST,FACST)
C****
C**** Add strait flow to JxL Strean Function
C****
C**** Input:   L = ocean model layer index
C****      OLNST = ocean strait mass flux (kg/s)
C****      FACST = global scaling factor for ocean straits
C**** Output: SF = JxL stream function (kg/s)
C****
      Use OCEAN,   Only: JM,LMO
      Use STRAITS, Only: NMST, IST,JST
      Use ODIAG,   Only: KBASIN
      Implicit None

      Integer*4,Intent(In) :: L
      Real*8,Intent(InOut) :: SF(JM,0:LMO,4)
      Real*8,Intent(In)    :: OLNST(LMO,NMST),FACST

C**** Local variables
      Integer*4 N, I1,J1,K1, I2,J2,K2

      Do 30 N=1,NMST
      I1 = IST(N,1)       ;  I2 = IST(N,2)
      J1 = JST(N,1)       ;  J2 = JST(N,2)
      K1 = KBASIN(I1,J1)  ;  K2 = KBASIN(I1,J1)
      If (J2 - J1) 10,30,20

C**** JST(N,2) < JST(N,1)
   10 If (K1 == K2)  Then
         SF(J2:J1-1,L,K1) = SF(J2:J1-1,L,K1) - OLNST(L+1,N)*FACST
      Else
         SF(J2:J1-1,L,K1) = SF(J2:J1-1,L,K1) - OLNST(L+1,N)*FACST*.5
         SF(J2:J1-1,L,K2) = SF(J2:J1-1,L,K2) - OLNST(L+1,N)*FACST*.5
      EndIf
      GoTo 30

C**** JST(N,2) > JST(N,1)
   20 If (K1 == K2)  Then
         SF(J1:J2-1,L,K1) = SF(J1:J2-1,L,K1) + OLNST(L+1,N)*FACST
      Else
         SF(J1:J2-1,L,K1) = SF(J1:J2-1,L,K1) + OLNST(L+1,N)*FACST*.5
         SF(J1:J2-1,L,K2) = SF(J1:J2-1,L,K2) + OLNST(L+1,N)*FACST*.5
      EndIf
   30 Continue
      Return
      EndSubroutine STRMJL_STRAITS

      Subroutine STRMIJ_STRAITS (J,SF,OLNST,FACST)
C****
C**** Add strait flow to IxJ Strean Function
C****
C**** Input:   J = ocean model latitude index
C****      OLNST = ocean strait mass flux (kg/s)
C****      FACST = global scaling factor for ocean straits
C**** Output: SF = IxJ stream function (kg/s)
C****
      Use OCEAN,   Only: IM,JM,LMO
      Use STRAITS, Only: NMST,LMST, IST,JST
      Implicit None

      Integer*4,Intent(In) :: J
      Real*8,Intent(InOut) :: SF(IM,JM)
      Real*8,Intent(In)    :: OLNST(LMO,NMST),FACST

C**** Local variables
      Integer*4 N,LM, I1,J1, I2,J2

      Do 40 N=1,NMST
      I1 = IST(N,1)       ;  I2 = IST(N,2)
      J1 = JST(N,1)       ;  J2 = JST(N,2)
      LM = LMST(N)
      If (J2 - J1) 10,20,30

C**** JST(N,2) < JST(N,1)
   10 If (J==J1 .or. J==J2)  Then
         SF(I1:I2-1,J) = SF(I1:I2-1,J) +
     *                   Sum(OLNST(1:LM,N))*FACST*.5/(J1-J2)
         GoTo 40  ;  EndIf
      If (J2 < J .and. J < J1)
     *   SF(I1:I2-1,J) = SF(I1:I2-1,J) +
     *                   Sum(OLNST(1:LM,N))*FACST/(J1-J2)
      GoTo 40

C**** JST(N,2) = JST(N,1)
   20 If (J==J1)
     *   SF(I1:I2-1,J) = SF(I1:I2-1,J) + Sum(OLNST(1:LM,N))*FACST
      GoTo 40

C**** JST(N,2) > JST(N,1)
   30 If (J==J1 .or. J==J2)  Then
         SF(I1:I2-1,J) = SF(I1:I2-1,J) +
     *                   Sum(OLNST(1:LM,N))*FACST*.5/(J2-J1)
         GoTo 40  ;  EndIf
      If (J1 < J .and. J < J2)
     *   SF(I1:I2-1,J) = SF(I1:I2-1,J) +
     *                   Sum(OLNST(1:LM,N))*FACST/(J2-J1)
   40 Continue
      Return
      EndSubroutine STRMIJ_STRAITS

      Subroutine OTJ_STRAITS (X,SOLNST,SCALE,KQ)
C****
C**** Calculate transport through straits from one latitude to another
C****
C**** Input: SOLNST = vertically integrated flux through each strait
C****         SCALE = strait scaling factor
C****            KQ = 1: mass flux; 2: heat flux; 3: salt flux
C**** Output:     X = northward flux as a function of J and basin
C****
      Use OCEAN,   Only: JM
      Use STRAITS, Only: NMST, IST,JST
      Use ODIAG,   Only: KBASIN
      Implicit None

      Real*8,Intent(InOut) :: X(0:JM,4,3)
      Real*8,Intent(In)    :: SOLNST(NMST),SCALE
      Integer*4,Intent(In) :: KQ

C**** Local variables
      Integer*4 N, I1,J1,K1, I2,J2,K2

      Do 30 N=1,NMST
      I1 = IST(N,1)       ;  I2 = IST(N,2)
      J1 = JST(N,1)       ;  J2 = JST(N,2)
      K1 = KBASIN(I1,J1)  ;  K2 = KBASIN(I2,J2)
      If (J2 - J1) 10,30,20

C**** JST(N,2) < JST(N,1)
   10 If (K1 == K2)  Then
        X(J2:J1-1,K1,KQ) = X(J2:J1-1,K1,KQ) - SOLNST(N)*SCALE
      Else
        X(J2:J1-1,K1,KQ) = X(J2:J1-1,K1,KQ) - SOLNST(N)*SCALE*.5
        X(J2:J1-1,K2,KQ) = X(J2:J1-1,K2,KQ) - SOLNST(N)*SCALE*.5
      Endif
        X(J2:J1-1, 4,KQ) = X(J2:J1-1, 4,KQ) - SOLNST(N)*SCALE
      GoTo 30

C**** JST(N,2) > JST(N,1)
   20 If (K1 == K2)  Then
        X(J1:J2-1,K1,KQ) = X(J1:J2-1,K1,KQ) + SOLNST(N)*SCALE
      Else
        X(J1:J2-1,K1,KQ) = X(J1:J2-1,K1,KQ) + SOLNST(N)*SCALE*.5
        X(J1:J2-1,K2,KQ) = X(J1:J2-1,K2,KQ) + SOLNST(N)*SCALE*.5
      Endif
        X(J1:J2-1, 4,KQ) = X(J1:J2-1, 4,KQ) + SOLNST(N)*SCALE
   30 Continue
      Return
      EndSubroutine OTJ_STRAITS
