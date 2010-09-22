#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_ON
#undef TRACERS_WATER
#endif
!@sum  LANDICE_DRV contains drivers for LANDICE related routines
!@auth Gavin Schmidt
!@ver  1.0
!@cont init_LI,PRECIP_LI,GROUND_LI,daily_LI

      SUBROUTINE init_LI(istart)
!@sum  init_ice initialises landice arrays
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : edpery,sday,lhm,tf
      USE FILEMANAGER
      USE MODEL_COM, only : im,jm,flice,focean,dtsrc
#ifdef SCM
     *                      ,I_TARG,J_TARG
#endif
      USE GEOM, only : axyp,imaxj,lat2d
      USE LANDICE, only: ace1li,ace2li,glmelt_on,glmelt_fac_nh
     *     ,glmelt_fac_sh,fwarea_sh,fwarea_nh,accpda,accpdg,eaccpda
     *     ,eaccpdg,snmin,dmicbimp,micbimp,deicbimp,eicbimp
#ifdef SCM
      USE SCMCOM, only : iu_scm_prt,SCM_SURFACE_FLAG,ATSKIN
#endif
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
     *     ,traccpda,traccpdg
#endif
#endif               /* TNL: inserted */
      USE LANDICE_COM, only : tlandi,snowli,loc_glm,mdwnimp,edwnimp
#ifdef TRACERS_WATER
     *     ,trsnowli,trlndi,trdwnimp
#endif
      USE FLUXES, only : gtemp,gmelt,egmelt,gtempr
#ifdef TRACERS_WATER
     *     ,gtracer
#ifdef TRACERS_OCEAN
     *     ,trgmelt
#endif
#endif
#ifdef TRACERS_OCEAN
#ifdef TRACERS_OCEAN_INDEP
#else
      USE TRACER_COM, only :  trglac,trw0
#endif
#endif
      USE DIAG_COM, only : npts,icon_MLI,icon_HLI,title_con,conpt0
     *     ,icon_MICB,icon_HICB
      USE PARAM
      USE DOMAIN_DECOMP_ATM, only : GRID,GET, GLOBALSUM, READT_PARALLEL
     &     ,AM_I_ROOT
      IMPLICIT NONE
      LOGICAL :: QCON(NPTS), T=.TRUE. , F=.FALSE.
      INTEGER, INTENT(IN) :: istart
      REAL*8 FAC_SH,FAC_NH
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     *                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     FWAREA_s,FWAREA_n,
     &     r8mask ! dummy array to read in 0/1 mask
      LOGICAL :: do_glmelt = .false.
      INTEGER I,J,N, iu_GL, I72
      INTEGER :: I_0,I_1, J_0,J_1
      CHARACTER*1, DIMENSION(IM,JM) :: CGLM   ! global array
      CHARACTER*72 :: TITLE

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** set GTEMP array for landice
      DO J=J_0,J_1
        DO I=I_0,I_1
          IF (FLICE(I,J).gt.0) THEN
            GTEMP(1:2,3,I,J)=TLANDI(1:2,I,J)
            GTEMPR(3,I,J)   =TLANDI(1,I,J)+TF
#ifdef SCM
            if (I.eq.I_TARG.and.J.eq.J_TARG) then
                if (SCM_SURFACE_FLAG.ge.1) then
                    GTEMP(1,3,I,J) = ATSKIN
                    GTEMP(2,3,I,J) = ATSKIN
                    GTEMPR(3,I,J) = ATSKIN + TF
                endif
            endif
#endif
#ifdef TRACERS_WATER
            if (istart.ge.9) then
            IF (SNOWLI(I,J).gt.SNMIN) THEN
              GTRACER(:,3,I,J)=TRSNOWLI(:,I,J)/SNOWLI(I,J)
            ELSE
              GTRACER(:,3,I,J)=TRLNDI(:,I,J)/(ACE1LI+ACE2LI)
            END IF
            end if
#endif
          END IF
        END DO
      END DO

C**** Calculate (fixed) iceberg melt terms from Antarctica and Greenland
      call sync_param("glmelt_on",glmelt_on)
      call sync_param("glmelt_fac_nh",glmelt_fac_nh)
      call sync_param("glmelt_fac_sh",glmelt_fac_sh)

C**** Note these parameters are highly resolution dependent!
C**** Around Antarctica, fresh water is added from 78S to 62S
C**** and from 0-60W and 135W to 165E
C**** Around Greenland the freshwater is added to both the east and
C**** west coasts in the one grid box closest to the coast which is
C**** determined by the direction in the river flow file.
C**** This information is now read in from the GLMELT file.

C**** Read in GLMELT file to distribute glacial melt (CGLM is global)
      IF (JM.gt.24) THEN ! for finer that old 8x10
         do_glmelt=.true.

#ifdef CUBED_SPHERE
c read a binary version of GLMELT
      call openunit("GLMELT",iu_GL,.true.,.true.)
      call readt_parallel(grid,iu_gl,nameunit(iu_gl),r8mask,1)
      loc_glm(i_0:i_1,j_0:j_1) = r8mask(i_0:i_1,j_0:j_1).eq.1d0
#else
c every PE reads the ascii version of GLMELT
      call openunit("GLMELT",iu_GL,.false.,.true.)
      READ  (iu_GL,'(A72)') TITLE
      IF(AM_I_ROOT())
     &     WRITE (6,*) 'Read on unit ',iu_GL,': ',TITLE
      READ  (iu_GL,*)
C**** assumes a 72-column width slab - will need adjusting for CS
      DO I72=1,1+(IM-1)/72
        DO J=JM,1,-1
          READ (iu_GL,'(72A1)')
     *         (CGLM(I,J),I=72*(I72-1)+1,MIN(IM,I72*72))
        END DO
      END DO
      DO J=J_0,J_1
        DO I=I_0,I_1
          LOC_GLM(I,J)=CGLM(I,J).eq."1" 
        ENDDO
      ENDDO
C****
#endif
      call closeunit(iu_GL)

C**** Calculate hemispheric areas (weighted by area and landmask)
C**** This could be extended by a different number in GLMELT to
C**** give ice sheet (rather than hemispheric) dependent areas
      FWAREA_s(:,:)=0
      FWAREA_n(:,:)=0
      DO J=J_0,J_1
        DO I=I_0,I_1
          IF (LOC_GLM(I,J) .and. FOCEAN(I,J).gt.0 ) THEN
            IF(LAT2D(I,J).LT.0.) THEN
              FWAREA_s(I,J)=AXYP(I,J)*FOCEAN(I,J)
            ELSE
              FWAREA_n(I,J)=AXYP(I,J)*FOCEAN(I,J)
            ENDIF
          END IF
        END DO
      END DO
      CALL GLOBALSUM(grid, FWAREA_s, FWAREA_SH, all=.true.)
      CALL GLOBALSUM(grid, FWAREA_n, FWAREA_NH, all=.true.)

      END IF

C**** Intialise gmelt fluxes
      GMELT = 0. ; EGMELT = 0.
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
      TRGMELT = 0.
#endif
#endif /* TNL: inserted */

      if (do_glmelt .and. glmelt_on > 0) then
C****  Note that water goes in with as if it is ice at 0 deg.
C****  Possibly this should be a function of in-situ freezing temp?
C****
C**** gmelt_fac_nh/sh should be used to match accumulation in control
C**** run so that global mean sea level is nearly constant. We can't
C**** get perfect values, but the drift should hopefully be less
C**** than 1 mm/year. This will minimise shock once interactive values
C**** are calculated.
C****
C**** Initial average mass fluxes for Antarctica/Greenland

      if (istart.lt.8 .or. ACCPDA.eq.0) then
C**** Initiallise total mass/energy fluxes (only at start of run)
C**** The net accumulation from IPCC2 report is 2016x10**12 kg/year
C**** for Antarctica and for Greenland it is 316x10**12 kg/year
        ACCPDA = glmelt_fac_sh*2016d12 ! kg/year
        ACCPDG = glmelt_fac_nh*316d12 ! kg/year
        EACCPDA = -LHM*ACCPDA ; EACCPDG = -LHM*ACCPDG ! J/year
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
        TRACCPDA(:) = trglac(:)*ACCPDA ! kg/year
        TRACCPDG(:) = trglac(:)*ACCPDG ! kg/year
#endif
#endif /* TNL: inserted */

C**** initiallise implicit accumulators (note that these arrays are
C**** not used until at least one full year has passed)
        MDWNIMP=0.  ;  EDWNIMP=0.
        MICBIMP=0.  ;  EICBIMP=0.
        dMICBIMP=0. ;  dEICBIMP=0.
#ifdef TRACERS_WATER
        TRDWNIMP=0.
c        TRICBIMP=0.
c        dTRICBIMP=0.
#endif
      end if

! accumulation (kg per source time step) per water column
      FAC_SH=DTsrc/(EDPERY*SDAY*FWAREA_SH)
      FAC_NH=DTsrc/(EDPERY*SDAY*FWAREA_NH)

C**** Set GL MELT arrays
#ifndef SCM
      DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
          IF (LOC_GLM(I,J) .and. lat2D(I,J).lt.0 ) THEN ! SH
             GMELT(I,J) =  ACCPDA*FAC_SH*AXYP(I,J)*FOCEAN(I,J) ! kg
            EGMELT(I,J) = EACCPDA*FAC_SH*AXYP(I,J)*FOCEAN(I,J) ! J
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
            TRGMELT(:,I,J)= TRACCPDA(:)*FAC_SH*AXYP(I,J)*FOCEAN(I,J)  ! kg
#endif
#endif /* TNL: inserted */
          END IF
          IF (LOC_GLM(I,J) .and. lat2D(I,J).gt.0 ) THEN ! NH
             GMELT(I,J) =  ACCPDG*FAC_NH*AXYP(I,J)*FOCEAN(I,J) ! kg
            EGMELT(I,J) = EACCPDG*FAC_NH*AXYP(I,J)*FOCEAN(I,J) ! J
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
            TRGMELT(:,I,J) = TRACCPDG(:)*FAC_NH*AXYP(I,J)*FOCEAN(I,J) ! kg
#endif
#endif /* TNL: inserted */
          END IF
        END DO
      END DO
#endif

      end if

C**** Set conservation diagnostics for land ice mass, energy
      QCON=(/ F, F, F, T, T, F, F, F, T, F, F/)
      CALL SET_CON(QCON,CONPT0,"LAND ICE","(KG/M^2)        ",
     *     "(10**-9 KG/SM^2)",1d0,1d9,icon_MLI)
      QCON=(/ F, F, F, T, T, F, F, F, T, F, F/)
      CALL SET_CON(QCON,CONPT0,"LNDI ENR","(10**6 J/M^2)   ",
     *     "(10**-3 W/M^2)  ",1d-6,1d3,icon_HLI)

C**** Set conservation diagnostics for implicit iceberg mass, energy
      QCON=(/ F, F, F, F, F, F, F, T, T, F, F/)
      CALL SET_CON(QCON,CONPT0,"ICEBERG ","(KG/M^2)        ",
     *     "(10**-7 KG/SM^2)",1d0,1d7,icon_MICB)
      QCON=(/ F, F, F, F, F, F, F, T, T, F, F/)
      CALL SET_CON(QCON,CONPT0,"ICBRG EN","(10**6 J/M^2)   ",
     *     "(10**-1 W/M^2)  ",1d-6,1d1,icon_HICB)
C****
      END SUBROUTINE init_LI

      SUBROUTINE PRECIP_LI
!@sum  PRECIP_LI driver for applying precipitation to land ice fraction
!@auth Original Development team
!@ver  1.0
!@calls LANDICE:PRECLI
      USE MODEL_COM, only : im,jm,flice,itlandi
#ifdef SCM
     *                      ,I_TARG,J_TARG
#endif
      USE CONSTANT, only : tf
      USE GEOM, only : imaxj,axyp,byaxyp
      USE FLUXES, only : runoli,prec,eprec,gtemp,gtempr
#ifdef SCM
      USE SCMCOM, only : iu_scm_prt,SCM_SURFACE_FLAG,ATSKIN
#endif
#ifdef TRACERS_WATER
     *     ,trunoli,trprec,gtracer
#endif
      USE LANDICE, only: ace1li,ace2li,precli,snmin
      USE LANDICE_COM, only : snowli,tlandi,mdwnimp,edwnimp
#ifdef TRACERS_WATER
     *     ,trsnowli,trlndi,ntm,trdwnimp
#endif
      USE DIAG_COM, only : aij=>aij_loc,jreg,ij_f0li,ij_f1li,ij_erun2
     *     ,ij_runli,j_run,j_implh,j_implm
      USE DOMAIN_DECOMP_ATM, only : GRID,GET
      IMPLICIT NONE

      REAL*8 SNOW,TG1,TG2,PRCP,ENRGP,EDIFS,DIFS,ERUN2,RUN0,PLICE,DXYPIJ
#ifdef TRACERS_WATER
!@var TRSNOW tracer amount in snow (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRSNOW
!@var TRLI tracer amount in land ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRLI
!@var TRPRCP tracer amount in precip (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRPRCP
!@var TRUN0 tracer runoff from ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRUN0
!@var TRDIFS implicit tracer flux at base of ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRDIFS
#endif
C**** Get useful grid parameters
      INTEGER :: I, J, JR
      INTEGER :: J_0, J_1, J_0H, J_1H ,I_0,I_1

      CALL GET(GRID,J_STRT=J_0      , J_STOP=J_1      ,
     &              J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
      DXYPIJ=AXYP(I,J)
      PLICE=FLICE(I,J)
      PRCP=PREC(I,J)
      JR=JREG(I,J)
      RUNOLI(I,J)=0
#ifdef TRACERS_WATER
      TRUNOLI(:,I,J)=0.
#endif
      IF (PLICE.gt.0 .and. PRCP.gt.0) THEN

        ENRGP=EPREC(I,J)      ! energy of precipitation
        SNOW=SNOWLI(I,J)
        TG1=TLANDI(1,I,J)
        TG2=TLANDI(2,I,J)
#ifdef TRACERS_WATER
        TRLI(:)=TRLNDI(:,I,J)
        TRSNOW(:)=TRSNOWLI(:,I,J)
        TRPRCP(:)=TRPREC(:,I,J)*BYAXYP(I,J)
#endif
        AIJ(I,J,IJ_F0LI)=AIJ(I,J,IJ_F0LI)+ENRGP*PLICE

        CALL PRECLI(SNOW,TG1,TG2,PRCP,ENRGP,
#ifdef TRACERS_WATER
     *       TRSNOW,TRLI,TRPRCP,TRDIFS,TRUN0,
#endif
     *       EDIFS,DIFS,ERUN2,RUN0)

C**** RESAVE PROGNOSTIC QUANTITIES AND FLUXES
        SNOWLI(I,J)=SNOW
        TLANDI(1,I,J)=TG1
        TLANDI(2,I,J)=TG2
        RUNOLI(I,J)  =RUN0
        GTEMP(1:2,3,I,J)=TLANDI(1:2,I,J)
        GTEMPR(3,I,J)   =TLANDI(1,I,J)+TF
#ifdef SCM
        if (I.eq.I_TARG.and.J.eq.J_TARG) then
            if (SCM_SURFACE_FLAG.ge.1) then
                GTEMP(1,3,I,J) = ATSKIN
                GTEMP(2,3,I,J) = ATSKIN
                GTEMPR(3,I,J) = ATSKIN + TF
            endif
        endif
#endif
C**** accumulate implicit fluxes for setting ocean balance
        MDWNIMP(I,J)=MDWNIMP(I,J)+DIFS *PLICE*DXYPIJ
        EDWNIMP(I,J)=EDWNIMP(I,J)+ERUN2*PLICE*DXYPIJ
#ifdef TRACERS_WATER
        TRLNDI(:,I,J)=TRLI(:)
        TRSNOWLI(:,I,J)=TRSNOW(:)
        TRUNOLI(:,I,J)=TRUN0(:)
        TRDWNIMP(:,I,J)=TRDWNIMP(:,I,J)+TRDIFS(:)*PLICE*DXYPIJ
        IF (SNOW.gt.SNMIN) THEN
          GTRACER(:,3,I,J)=TRSNOW(:)/SNOW
        ELSE
          GTRACER(:,3,I,J)=TRLI(:)/(ACE1LI+ACE2LI)
        END IF
#endif
C**** ACCUMULATE DIAGNOSTICS
c       CALL INC_AJ(I,J,ITLANDI,J_TYPE,       PLICE)
        CALL INC_AJ(I,J,ITLANDI,J_RUN ,  RUN0*PLICE)
C       CALL INC_AJ(I,J,ITLANDI,J_ERUN, ERUN0*PLICE) ! (Tg=0)
        CALL INC_AJ(I,J,ITLANDI,J_IMPLM, DIFS*PLICE)
        CALL INC_AJ(I,J,ITLANDI,J_IMPLH,ERUN2*PLICE)
        CALL INC_AREG(I,J,JR,J_RUN ,  RUN0*PLICE)
c       CALL INC_AREG(I,J,JR,J_ERUN, ERUN0*PLICE) ! (Tg=0)
        CALL INC_AREG(I,J,JR,J_IMPLM, DIFS*PLICE)
        CALL INC_AREG(I,J,JR,J_IMPLH,ERUN2*PLICE)
        AIJ(I,J,IJ_F1LI) =AIJ(I,J,IJ_F1LI) +EDIFS*PLICE
        AIJ(I,J,IJ_ERUN2)=AIJ(I,J,IJ_ERUN2)+ERUN2*PLICE
        AIJ(I,J,IJ_RUNLI)=AIJ(I,J,IJ_RUNLI)+RUN0 *PLICE
      END IF
      END DO
      END DO

      END SUBROUTINE PRECIP_LI

      SUBROUTINE GROUND_LI
!@sum  GROUND_LI driver for applying surface fluxes to land ice fraction
!@auth Original Development team
!@ver  1.0
!@calls LANDICE:LNDICE
      USE CONSTANT, only : tf,sday,edpery
      USE MODEL_COM, only : im,jm,flice,itlandi,itocean,itoice,dtsrc
#ifdef SCM
     *                      ,I_TARG,J_TARG
#endif
      USE GEOM, only : imaxj,axyp,byaxyp
      USE LANDICE, only : lndice,ace1li,ace2li,snmin,micbimp
     *     ,eicbimp,accpda,accpdg,eaccpda,eaccpdg
      USE SEAICE_COM, only : rsi
      USE SEAICE, only : rhos
#ifdef SCM
      USE SCMCOM, only : iu_scm_prt,SCM_SURFACE_FLAG,ATSKIN
#endif
      USE DIAG_COM, only : aij=>aij_loc,jreg,ij_runli,ij_f1li,ij_erun2
     *     ,j_wtr1,j_ace1,j_wtr2,j_ace2,j_snow,j_run
     *     ,j_implh,j_implm,j_rsnow,ij_rsnw,ij_rsit,ij_snow,ij_f0oc
     *     ,j_rvrd,j_ervr,ij_mrvr,ij_ervr,ij_zsnow,ij_fwoc,ij_li
      USE LANDICE_COM, only : snowli,tlandi,mdwnimp,edwnimp
#ifdef TRACERS_WATER
     *     ,ntm,trsnowli,trlndi,trdwnimp  !,tricbimp,traccpda,traccpdg
#endif
      USE FLUXES, only : e0,e1,evapor,gtemp,runoli,gmelt,egmelt,gtempr
#ifdef TRACERS_WATER
     *     ,trunoli,trevapor,gtracer
#ifdef TRACERS_DRYDEP
     *     ,trdrydep
#endif
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
     *     ,trgmelt
#endif   /* TNL: inserted */
      USE TRDIAG_COM, only : taijn=>taijn_loc , tij_rvr
#endif
#endif
      USE DOMAIN_DECOMP_ATM, only : GRID,GET
      USE TimerPackage_mod, only: startTimer => start
      USE TimerPackage_mod, only: stopTimer => stop
      IMPLICIT NONE

      REAL*8 SNOW,TG1,TG2,F0DT,F1DT,EVAP,EDIFS,DIFS,RUN0,PLICE,DXYPIJ
     *     ,SCOVLI
#ifdef TRACERS_WATER
!@var TRSNOW tracer amount in snow (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRSNOW
!@var TRLI tracer amount in land ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRLI
!@var TREVAP tracer amount in evaporation (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TREVAP
!@var TRUN0 tracer runoff from ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRUN0
!@var TRDIFS implicit tracer flux at base of ice (kg/m^2)
      REAL*8, DIMENSION(NTM) :: TRDIFS
#endif
      INTEGER I,J,JR,N
      INTEGER :: J_0,J_1, J_0H, J_1H ,I_0,I_1

      call startTimer('GROUND_LI()')
      CALL GET(GRID,J_STRT=J_0      ,J_STOP=J_1
     &             ,J_STRT_HALO=J_0H,J_STOP_HALO=J_1H)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
      DXYPIJ=AXYP(I,J)
      PLICE=FLICE(I,J)
      JR=JREG(I,J)
      RUNOLI(I,J)=0.
      IF (PLICE.gt.0) THEN

        SNOW=SNOWLI(I,J)
        TG1=TLANDI(1,I,J)
        TG2=TLANDI(2,I,J)
        F0DT=E0(I,J,3)
        F1DT=E1(I,J,3)
        EVAP=EVAPOR(I,J,3)
#ifdef TRACERS_WATER
        TRLI(:)=TRLNDI(:,I,J)
        TRSNOW(:)=TRSNOWLI(:,I,J)
        TREVAP(:)=TREVAPOR(:,3,I,J)
#ifdef TRACERS_DRYDEP
     *       -trdrydep(:,3,i,j)
#endif
#endif
        CALL LNDICE(SNOW,TG1,TG2,F0DT,F1DT,EVAP,
#ifdef TRACERS_WATER
     *     TRSNOW,TRLI,TREVAP,TRDIFS,TRUN0,
#endif
     *     EDIFS,DIFS,RUN0)
C**** RESAVE PROGNOSTIC QUANTITIES AND FLUXES
        SNOWLI(I,J)=SNOW
        TLANDI(1,I,J)=TG1
        TLANDI(2,I,J)=TG2
        RUNOLI(I,J) = RUN0
        GTEMP(1:2,3,I,J)=TLANDI(1:2,I,J)
        GTEMPR(3,I,J)   =TLANDI(1,I,J)+TF
#ifdef SCM
        if (I.eq.I_TARG.and.J.eq.J_TARG) then
            if (SCM_SURFACE_FLAG.ge.1) then
                GTEMP(1,3,I,J) = ATSKIN
                GTEMP(2,3,I,J) = ATSKIN
                GTEMPR(3,I,J) = ATSKIN + TF
            endif
        endif
#endif
C**** accumulate implicit fluxes for setting ocean balance
        MDWNIMP(I,J)=MDWNIMP(I,J)+DIFS *PLICE*DXYPIJ
        EDWNIMP(I,J)=EDWNIMP(I,J)+EDIFS*PLICE*DXYPIJ
#ifdef TRACERS_WATER
        TRLNDI(:,I,J)=TRLI(:)
        TRSNOWLI(:,I,J)=TRSNOW(:)
        TRUNOLI(:,I,J)=TRUN0(:)
        TRDWNIMP(:,I,J)=TRDWNIMP(:,I,J)+TRDIFS(:)*PLICE*DXYPIJ
        IF (SNOW.gt.SNMIN) THEN
          GTRACER(:,3,I,J)=TRSNOW(:)/SNOW
        ELSE
          GTRACER(:,3,I,J)=TRLI(:)/(ACE1LI+ACE2LI)
        END IF
#endif
C**** ACCUMULATE DIAGNOSTICS
        SCOVLI=0
        IF (SNOWLI(I,J).GT.0.) SCOVLI=PLICE
        CALL INC_AJ(I,J,ITLANDI,J_RSNOW,SCOVLI)
        CALL INC_AREG(I,J,JR,J_RSNOW,SCOVLI)
        AIJ(I,J,IJ_RSNW)=AIJ(I,J,IJ_RSNW)+SCOVLI
        AIJ(I,J,IJ_SNOW)=AIJ(I,J,IJ_SNOW)+SNOW*PLICE
        AIJ(I,J,IJ_RSIT)=AIJ(I,J,IJ_RSIT)+PLICE
        AIJ(I,J,IJ_LI)  =AIJ(I,J,IJ_LI)  +PLICE
        AIJ(I,J,IJ_ZSNOW)=AIJ(I,J,IJ_ZSNOW)+PLICE*SNOW/RHOS

        CALL INC_AJ(I,J,ITLANDI,J_RUN,   RUN0*PLICE)
        CALL INC_AJ(I,J,ITLANDI,J_SNOW,  SNOW*PLICE)
        CALL INC_AJ(I,J,ITLANDI,J_ACE1,ACE1LI*PLICE)
        CALL INC_AJ(I,J,ITLANDI,J_ACE2,ACE2LI*PLICE)
        CALL INC_AJ(I,J,ITLANDI,J_IMPLH,EDIFS*PLICE)
        CALL INC_AJ(I,J,ITLANDI,J_IMPLM, DIFS*PLICE)

        CALL INC_AREG(I,J,JR,J_RUN  ,  RUN0*PLICE)
        CALL INC_AREG(I,J,JR,J_SNOW ,  SNOW*PLICE)
        CALL INC_AREG(I,J,JR,J_ACE1 ,ACE1LI*PLICE)
        CALL INC_AREG(I,J,JR,J_ACE2 ,ACE2LI*PLICE)
        CALL INC_AREG(I,J,JR,J_IMPLH, EDIFS*PLICE)
        CALL INC_AREG(I,J,JR,J_IMPLM,  DIFS*PLICE)

        AIJ(I,J,IJ_F1LI) =AIJ(I,J,IJ_F1LI) +(EDIFS+F1DT)*PLICE
        AIJ(I,J,IJ_RUNLI)=AIJ(I,J,IJ_RUNLI)+RUN0 *PLICE
        AIJ(I,J,IJ_ERUN2)=AIJ(I,J,IJ_ERUN2)+EDIFS*PLICE
      END IF

C**** Accumulate diagnostics related to iceberg flux here also
C**** remove source time step icb flux 
      MICBIMP(:) = MICBIMP(:)-(/ ACCPDA,  ACCPDG /)*DTsrc/(EDPERY*SDAY)
      EICBIMP(:) = EICBIMP(:)-(/EACCPDA, EACCPDG /)*DTsrc/(EDPERY*SDAY)

      CALL INC_AJ(I,J,ITOCEAN,J_RVRD,(1.-RSI(I,J))* GMELT(I,J)
     *     *BYAXYP(I,J))
      CALL INC_AJ(I,J,ITOCEAN,J_ERVR,(1.-RSI(I,J))*EGMELT(I,J)
     *     *BYAXYP(I,J))
      CALL INC_AJ(I,J,ITOICE,J_RVRD,RSI(I,J)* GMELT(I,J)*BYAXYP(I,J))
      CALL INC_AJ(I,J,ITOICE,J_ERVR,RSI(I,J)*EGMELT(I,J)*BYAXYP(I,J))
      AIJ(I,J,IJ_F0OC) = AIJ(I,J,IJ_F0OC)+EGMELT(I,J)*BYAXYP(I,J)
      AIJ(I,J,IJ_FWOC) = AIJ(I,J,IJ_FWOC)+ GMELT(I,J)*BYAXYP(I,J)

      CALL INC_AREG(I,J,JR,J_RVRD, GMELT(I,J)*BYAXYP(I,J))
      CALL INC_AREG(I,J,JR,J_ERVR,EGMELT(I,J)*BYAXYP(I,J))
      AIJ(I,J,IJ_MRVR)=AIJ(I,J,IJ_MRVR) +  GMELT(I,J)
      AIJ(I,J,IJ_ERVR)=AIJ(I,J,IJ_ERVR) + EGMELT(I,J)
#ifdef TRACERS_WATER  /* TNL: inserted */
c      DO N=1,NTM
c        TRICBIMP(N,:) = TRICBIMP(N,:)-(/TRACCPDA(:), TRACCPDG(:)/)
c     *       *DTsrc/(EDPERY*SDAY)
c      END DO
#ifdef TRACERS_OCEAN
      TAIJN(I,J,TIJ_RVR,:)=TAIJN(I,J,TIJ_RVR,:)+ TRGMELT(:,I,J)
     *     *BYAXYP(I,J)
#endif
#endif  /* TNL: inserted */
C****
      END DO
      END DO

      call stopTimer('GROUND_LI()')
      END SUBROUTINE GROUND_LI

      SUBROUTINE conserv_MLI(ICE)
!@sum  conserv_MLI calculates total amount of snow and ice in land ice
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm,fim,flice
      USE GEOM, only : imaxj
      USE LANDICE_COM, only : snowli
      USE LANDICE, only : lndice,ace1li,ace2li
      USE DOMAIN_DECOMP_ATM, only : GRID,GET
      IMPLICIT NONE
!@var ICE total land ice snow and ice mass (kg/m^2)
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) :: ICE
      INTEGER I,J
      INTEGER :: J_0,J_1 ,I_0,I_1
      LOGICAl :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1,
     &         HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        ICE(I,J)=FLICE(I,J)*(ACE1LI+ACE2LI+SNOWLI(I,J))
      END DO
      END DO
      IF(HAVE_SOUTH_POLE) ICE(2:im,1) =ICE(1,1)
      IF(HAVE_NORTH_POLE) ICE(2:im,JM)=ICE(1,JM)

      RETURN
C****
      END SUBROUTINE conserv_MLI

      SUBROUTINE conserv_HLI(EICE)
!@sum  conserv_HLI calculates total land ice energy
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : shi,lhm
      USE MODEL_COM, only : im,jm,fim,flice
      USE GEOM, only : imaxj
      USE LANDICE_COM, only : snowli,tlandi
      USE LANDICE, only : ace1li,ace2li
      USE DOMAIN_DECOMP_ATM, only : GRID,GET
      IMPLICIT NONE
!@var EICE total land ice energy (J/m^2)
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) :: EICE
      INTEGER I,J
      INTEGER :: J_0,J_1 ,I_0,I_1
      LOGICAl :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1,
     &         HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        EICE(I,J)=FLICE(I,J)*((TLANDI(1,I,J)*SHI-LHM)*(ACE1LI
     *         +SNOWLI(I,J))+(TLANDI(2,I,J)*SHI-LHM)*ACE2LI)
      END DO
      END DO
      IF(HAVE_SOUTH_POLE) EICE(2:im,1) =EICE(1,1)
      IF(HAVE_NORTH_POLE) EICE(2:im,JM)=EICE(1,JM)
      RETURN
C****
      END SUBROUTINE conserv_HLI

      SUBROUTINE conserv_MICB(MICB)
!@sum  conserv_MICB calculates total amount of implicit iceberg
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm
      USE GEOM, only : imaxj,axyp
      USE LANDICE, only : micbimp
      USE DOMAIN_DECOMP_ATM, only : GRID,GET
      IMPLICIT NONE
!@var MICB implicit mass of iceberg (kg/m^2)
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) :: MICB
      INTEGER I,J
      INTEGER :: J_0,J_1 ,I_0,I_1
      LOGICAl :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1,
     &         HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** since this is just diagnostic of hemispheric mean for now, we put
C**** everything in one spot per hemisphere.
C**** (how should we do this in CUBED_SPHERE?)
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        MICB(I,J) = 0.
        IF (I.EQ.1 .and. J.eq.2)    MICB(I,J)=micbimp(1)/AXYP(I,J)
        IF (I.EQ.1 .and. J.eq.JM-1) MICB(I,J)=micbimp(2)/AXYP(I,J)
      END DO
      END DO
      IF(HAVE_SOUTH_POLE) MICB(2:im,1) =MICB(1,1)
      IF(HAVE_NORTH_POLE) MICB(2:im,JM)=MICB(1,JM)

      RETURN
C****
      END SUBROUTINE conserv_MICB

      SUBROUTINE conserv_HICB(EICB)
!@sum  conserv_HICB calculates implicit iceberg energy
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : im,jm
      USE GEOM, only : imaxj,axyp
      USE LANDICE, only : eicbimp
      USE DOMAIN_DECOMP_ATM, only : GRID,GET
      IMPLICIT NONE
!@var EICB impicit iceberg energy (J/m^2)
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) :: EICB
      INTEGER I,J
      INTEGER :: J_0,J_1 ,I_0,I_1
      LOGICAl :: HAVE_SOUTH_POLE,HAVE_NORTH_POLE

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1,
     &         HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** since this is just diagnostic of hemispheric mean for now, we put
C**** everything in one spot per hemisphere.
C**** (how should we do this in CUBED_SPHERE?)
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        EICB(I,J) = 0.
        IF (I.EQ.1 .and. J.eq.2)    EICB(I,J)=eicbimp(1)/AXYP(I,J)
        IF (I.EQ.1 .and. J.eq.JM-1) EICB(I,J)=eicbimp(2)/AXYP(I,J)
      END DO
      END DO
      IF(HAVE_SOUTH_POLE) EICB(2:im,1) =EICB(1,1)
      IF(HAVE_NORTH_POLE) EICB(2:im,JM)=EICB(1,JM)
      RETURN
C****
      END SUBROUTINE conserv_HICB

      SUBROUTINE daily_LI
!@sum  daily_ice does daily landice things
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : edpery,sday,lhm,shi
      USE MODEL_COM, only : im,jm,flice,focean,dtsrc,jday,jyear
     *     ,itime,itimei,nday,JDperY
      USE GEOM, only : axyp,imaxj,lat2d
      USE LANDICE, only: ace1li,ace2li,glmelt_on,glmelt_fac_nh
     *     ,glmelt_fac_sh,fwarea_sh,fwarea_nh,accpda,accpdg,eaccpda
     *     ,eaccpdg,micbimp,eicbimp,dmicbimp,deicbimp
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
     *     ,traccpda,traccpdg   ! tricbimp,dtricbimp
#endif
      USE TRDIAG_COM, only : to_per_mil
      USE TRACER_COM, only :  trw0, nWater, trname
#endif    /* TNL: inserted */
      USE LANDICE_COM, only : tlandi,snowli,mdwnimp,edwnimp,loc_glm
#ifdef TRACERS_WATER
     *     ,ntm,trsnowli,trlndi,trdwnimp
#endif
      USE FLUXES, only : gmelt,egmelt
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
     *     ,trgmelt
#endif
#endif    /* TNL: inserted */
      USE PARAM
      USE DOMAIN_DECOMP_ATM, only : GRID, GET, GLOBALSUM, AM_I_ROOT,
     *     ESMF_BCAST
      IMPLICIT NONE
!@var gm_relax Glacial Melt relaxation parameter (1/year)
      REAL*8, PARAMETER :: gm_relax = 0.1d0  ! 10 year relaxation

      REAL*8 mdwnimp_SH,mdwnimp_NH,edwnimp_SH,edwnimp_NH,FAC_SH,FAC_NH
#ifdef TRACERS_WATER
      REAL*8 trdwnimp_SH(NTM),trdwnimp_NH(NTM)
#endif
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO) :: mask_s,arr_s,arr_n
      INTEGER :: J_0,J_1,I_0,I_1,I,J,ITM

      CALL GET(GRID,J_STRT=J_0,J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** Keep track of accumlated implicit iceberg amount
C**** hemispheric masking (why isn't this in GEOM?)
      do j=j_0,j_1; do i=i_0,i_1
        if(lat2d(i,j).lt.0.) then
          mask_s(i,j) = 1.
        else
          mask_s(i,j) = 0.
        endif
      enddo; enddo

C**** For total cumulative accumulation, pre-remove accumulation up to
C**** the day before

      MICBIMP(:) = MICBIMP(:) - dMICBIMP(:) 
      EICBIMP(:) = EICBIMP(:) - dEICBIMP(:) 
c#ifdef TRACERS_WATER
c      TRICBIMP(:,:) = TRCIBIMP(:,:) - dTRICBIMP(:,:) 
c#endif

C**** Calculate and save mass/energy/tracer accumulation to date for
C**** this year
! mass and energy (kg, J)
      do j=j_0,j_1; do i=i_0,i_1
        arr_s(i,j) = mdwnimp(i,j)*mask_s(i,j)
        arr_n(i,j) = mdwnimp(i,j)*(1.-mask_s(i,j))
      enddo; enddo
      CALL GLOBALSUM(grid, arr_s, mdwnimp_SH ,ALL=.TRUE.)
      CALL GLOBALSUM(grid, arr_n, mdwnimp_NH ,ALL=.TRUE.)
      do j=j_0,j_1; do i=i_0,i_1
        arr_s(i,j) = edwnimp(i,j)*mask_s(i,j)
        arr_n(i,j) = edwnimp(i,j)*(1.-mask_s(i,j))
      enddo; enddo
      CALL GLOBALSUM(grid, arr_s, edwnimp_SH ,ALL=.TRUE.)
      CALL GLOBALSUM(grid, arr_n, edwnimp_NH ,ALL=.TRUE.)

      dMICBIMP(:)=(/ mdwnimp_SH, mdwnimp_NH /)
      dEICBIMP(:)=(/ edwnimp_SH, edwnimp_NH /)
      
#ifdef TRACERS_WATER
      DO ITM=1,NTM
        do j=j_0,j_1; do i=i_0,i_1
          arr_s(i,j) = trdwnimp(itm,i,j)*mask_s(i,j)
          arr_n(i,j) = trdwnimp(itm,i,j)*(1.-mask_s(i,j))
        enddo; enddo
        CALL GLOBALSUM(grid, arr_s, trdwnimp_SH(itm) ,ALL=.TRUE.)
        CALL GLOBALSUM(grid, arr_n, trdwnimp_NH(itm) ,ALL=.TRUE.)
c        dTRICBIMP(ITM,:)=(/ trdwnimp_SH(itm), trdwnimp_NH(itm) /)
      END DO
#endif

C**** For total cumulative accumulation, add back accumulation up to today
      MICBIMP(:) = MICBIMP(:) + dMICBIMP(:) 
      EICBIMP(:) = EICBIMP(:) + dEICBIMP(:) 
c#ifdef TRACERS_WATER
c      TRICBIMP(:,:) = TRCIBIMP(:,:) + dTRICBIMP(:,:) 
c#endif

C**** Every year update gmelt factors in order to balance downward
C**** implicit fluxes. If this is used, then ice sheets/snow are FORCED to
C**** be in balance. This may not be appropriate for transient runs but
C**** we aren't getting that right anyway.

      IF (JDAY.eq.1) THEN   ! Jan 1. only

! only adjust after at least one full year
        IF (itime.ge.itimei+JDperY*nday .and. glmelt_on==1) THEN

C*** prevent iceberg sucking
          if(mdwnimp_NH.lt.0) then
            write(99,*) "Limiting NH icesheet replacement mass/energy",
     *           mdwnimp_NH,edwnimp_NH
            mdwnimp_NH=0.
            edwnimp_NH=0.
#ifdef TRACERS_WATER
            trdwnimp_NH=0.
#endif
          endif

          if(mdwnimp_SH.lt.0) then
            write(99,*) "Limiting SH icesheet replacement mass/energy",
     *           mdwnimp_SH,edwnimp_SH
            mdwnimp_SH=0.
            edwnimp_SH=0.
#ifdef TRACERS_WATER
            trdwnimp_SH=0.
#endif
          endif

#ifdef TRACERS_WATER
          DO ITM=1,NTM
            if (to_per_mil(itm).gt.0) then
              if (mdwnimp_NH .gt. 0) print*,'Tracers(NH) (permil):',itm
     *             ,(trdwnimp_NH(ITM)/mdwnimp_NH/trw0(itm)-1.)*1000
     *             ,trdwnimp_NH(ITM),mdwnimp_NH
              if (mdwnimp_SH .gt. 0)  print*,'Tracers(SH) (permil):',itm
     *             ,(trdwnimp_SH(ITM)/mdwnimp_SH/trw0(itm)-1.)*1000
     *             ,trdwnimp_SH(ITM),mdwnimp_SH
            else
              print*,'Tracers(NH) (kg):',itm,trdwnimp_NH(ITM),mdwnimp_NH
              print*,'Tracers(SH) (kg):',itm,trdwnimp_SH(ITM),mdwnimp_SH
            end if
          END DO
#endif

C**** adjust hemispheric mean glacial melt amounts (only on root processor)
      if (AM_I_ROOT()) THEN
          write(6,*) "Adjusting glacial melt: ", jday,jyear
          write(6,*) "Mass (before): ",accpda,accpdg,mdwnimp_SH
     *         ,mdwnimp_NH
          write(6,*) "Temp (before): ",(eaccpda/accpda+lhm)/shi
     *         ,(eaccpdg/accpdg+lhm)/shi
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
          do itm=1,ntm
            if (to_per_mil(itm).gt.0) then
              write(6,*) trim(trname(itm))," (before) (permil)",1000
     *             *(traccpda(itm)/accpda/trw0(itm)-1.)
            else
              write(6,*) trim(trname(itm))," (before) (kg)",traccpda(itm
     *             ),accpda
            end if
          end do
#endif
#endif   /* TNL: inserted */
           accpda =  accpda + gm_relax*(mdwnimp_SH -  accpda)
           accpdg =  accpdg + gm_relax*(mdwnimp_NH -  accpdg)
          eaccpda = eaccpda + gm_relax*(edwnimp_SH - eaccpda)
          eaccpdg = eaccpdg + gm_relax*(edwnimp_NH - eaccpdg)
          write(6,*) "Mass (after): ",accpda,accpdg
          write(6,*) "Temp (after): ",(eaccpda/accpda+lhm)/shi,
     *         (eaccpdg/accpdg+lhm)/shi
        END IF
        call ESMF_BCAST(grid,  ACCPDA)
        call ESMF_BCAST(grid,  ACCPDG)
        call ESMF_BCAST(grid, EACCPDA)
        call ESMF_BCAST(grid, EACCPDG)

#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
        if (AM_I_ROOT()) THEN
          traccpda(:)=traccpda(:)+gm_relax*(trdwnimp_SH(:)-traccpda(:))
          traccpdg(:)=traccpdg(:)+gm_relax*(trdwnimp_NH(:)-traccpdg(:))

          do itm=1,ntm
            if (to_per_mil(itm).gt.0) then
              write(6,*) trim(trname(itm))," (after) (permil)",1000
     *             *(traccpda(itm)/accpda/trw0(itm)-1.)
            else
              write(6,*) trim(trname(itm))," (after) (kg)",traccpda(itm)
     *             ,accpda
            end if
          end do
        ENDIF
        call ESMF_BCAST(grid, TRACCPDA)
        call ESMF_BCAST(grid, TRACCPDG)
#endif
#endif   /* TNL: inserted */
! accumulation (kg per source time step) per water column
      FAC_SH=DTsrc/(EDPERY*SDAY*FWAREA_SH)
      FAC_NH=DTsrc/(EDPERY*SDAY*FWAREA_NH)

C**** Set GL MELT arrays
      DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
          IF (LOC_GLM(I,J) .and. lat2D(I,J).lt.0 ) THEN ! SH
              GMELT(I,J)  =  ACCPDA*FAC_SH*AXYP(I,J)*FOCEAN(I,J)  ! kg
              EGMELT(I,J) = EACCPDA*FAC_SH*AXYP(I,J)*FOCEAN(I,J) ! J
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
              TRGMELT(:,I,J)=TRACCPDA(:)*FAC_SH*AXYP(I,J)*FOCEAN(I,J)  ! kg
#endif
#endif    /* TNL: inserted */
          END IF
          IF (LOC_GLM(I,J) .and. lat2D(I,J).gt.0 ) THEN ! NH
             GMELT(I,J) =  ACCPDG*FAC_NH*AXYP(I,J)*FOCEAN(I,J) ! kg
            EGMELT(I,J) = EACCPDG*FAC_NH*AXYP(I,J)*FOCEAN(I,J) ! J
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
            TRGMELT(:,I,J) = TRACCPDG(:)*FAC_NH*AXYP(I,J)*FOCEAN(I,J)  ! kg
#endif
#endif    /* TNL: inserted */
          END IF
        END DO
      END DO

      END IF  ! run not in its first year

C**** reset implicit accumulators
      MDWNIMP=0.
      EDWNIMP=0.
#ifdef TRACERS_WATER
      TRDWNIMP=0.
#endif

      END IF  ! start of new year

      END SUBROUTINE daily_LI

      SUBROUTINE CHECKLI (SUBR)
!@sum  CHECKLI checks whether the land ice variables are reasonable.
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0 (based on LB265)
      USE CONSTANT, only : teeny
      USE MODEL_COM, only : im,jm,qcheck,flice
      USE DOMAIN_DECOMP_ATM, only : HALO_UPDATE, GET, GRID
      USE GEOM, only : imaxj
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm, trname
#endif
      USE LANDICE, only : ace1li, ace2li
      USE LANDICE_COM, only : tlandi,snowli,mdwnimp,edwnimp
#ifdef TRACERS_WATER
     *     ,trsnowli,trlndi,trdwnimp
#endif
      IMPLICIT NONE
      INTEGER :: J_0,J_1,J_0H,J_1H,J_0S,J_1S,I_0H,I_1H,I_0,I_1,njpol
      INTEGER I,J,N !@var I,J loop variables
      CHARACTER*6, INTENT(IN) :: SUBR
      LOGICAL QCHECKL
#ifdef TRACERS_WATER
      integer :: imax,jmax
      real*8 relerr,errmax
#endif
      CALL GET(grid, J_STRT=J_0,      J_STOP=J_1,
     *               J_STRT_HALO=J_0H,J_STOP_HALO=J_1H,
     &               J_STRT_SKP=J_0S, J_STOP_SKP=J_1S)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      njpol = grid%J_STRT_SKP-grid%J_STRT

C**** Check for NaN/INF in land ice data 
      CALL CHECK3C(tlandi(:,I_0:I_1,J_0:J_1),2,I_0,I_1,J_0,J_1,NJPOL,
     &     SUBR,'tli')
      CALL CHECK3B(snowli(I_0:I_1,J_0:J_1) ,I_0,I_1,J_0,J_1,NJPOL,1,
     &     SUBR,'sli')
      CALL CHECK3B(mdwnimp(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &     SUBR,'mdw')
      CALL CHECK3B(edwnimp(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &     SUBR,'edw')

      QCHECKL = .FALSE.
#ifdef TRACERS_WATER
      do n=1,ntm
C**** Check conservation of water tracers in land ice
        if (trname(n).eq.'Water') then
          errmax = 0. ; imax=I_0 ; jmax=J_0
          do j=J_0, J_1
          do i=I_0,imaxj(j)
            if (flice(i,j).gt.0) then
              relerr=max(
     *             abs(trlndi(n,i,j)/(ace1li+ace2li)-1.)
     *             ,abs((trsnowli(n,i,j)-snowli(i,j))/(snowli(i,j)+teeny
     *             )),abs((trdwnimp(n,i,j)-mdwnimp(i,j))/(mdwnimp(i,j)
     *             +teeny)))
            else
              relerr=abs(trdwnimp(n,i,j)-mdwnimp(i,j))/(mdwnimp(i,j)
     *             +teeny)
            end if
            if (relerr.gt.errmax) then
              imax=i ; jmax=j ; errmax=relerr
            end if
          end do
          end do
          print*,"Relative error in land ice mass after ",trim(subr),":"
     *         ,imax,jmax,errmax,trlndi(n,imax,jmax),ace1li+ace2li
     *         ,trsnowli(n,imax,jmax),snowli(imax,jmax),trdwnimp(n,imax
     *         ,jmax),mdwnimp(imax,jmax)
        end if
      end do
#endif

      IF (QCHECKL)
     &call stop_model('CHECKLI: Land Ice variables out of bounds',255)
      RETURN
C****
      END SUBROUTINE CHECKLI
