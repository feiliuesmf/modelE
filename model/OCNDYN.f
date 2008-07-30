#include "rundeck_opts.h"

      SUBROUTINE OCEANS
!@sum  OCEANS integrates ocean source terms and dynamics
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : rhows,grav
      USE MODEL_COM, only : idacc,modd5s,msurf
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : t_qlimit,ntm
      USE OCEAN, only : trmo,txmo,tymo,tzmo
C?*** For serial GM/straits computations, pack data into global arrays
      USE OCEAN, only : trmo_glob,txmo_glob,tymo_glob,tzmo_glob
#endif
      USE OCEAN, only : im,jm,lmo,ndyno,mo,g0m,gxmo,gymo,gzmo
     *     ,s0m,sxmo,symo,szmo,dts,dtofs,dto,dtolf,mdyno,msgso
     *     ,ogeoz,ogeoz_sv,opbot,ze,lmm,imaxj, UO,VONP,IVNP !,VOSP,IVSP
     *     ,OBottom_drag,OCoastal_drag,focean
C?*** For serial GM/straits computations, pack data into global arrays
      USE OCEAN, only : S0M_glob,SXMO_glob,SYMO_glob,SZMO_glob
      USE OCEAN, only : G0M_glob,GXMO_glob,GYMO_glob,GZMO_glob

      USE ODIAG, only : oijl=>oijl_loc,oij=>oij_loc
     *     ,ijl_mo,ijl_g0m,ijl_s0m,  ijl_gflx,ijl_sflx
     *     ,ijl_mfu,ijl_mfv,ijl_mfw, ijl_ggmfl,ijl_sgmfl,ij_ssh,ij_pb
#ifdef TRACERS_OCEAN
     *     ,toijl=>toijl_loc, toijl_conc,toijl_tflx,toijl_gmfl
#endif
      USE OCEAN_DYN, only : mmi,smu,smv,smw
      USE DOMAIN_DECOMP, only : grid,get, haveLatitude

C?*** For serial ODIF/GM/straits computations:
      USE DOMAIN_DECOMP, only : AM_I_ROOT, pack_data, unpack_data
      USE OCEAN, only : scatter_ocean, gather_ocean
      USE OCEAN, only : scatter_ocean_straits, gather_ocean_straits

      IMPLICIT NONE
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) ::
     *     MM0,MM1,MMX,UM0,VM0,UM1,VM1
      INTEGER NS,I,J,L,mnow,n

c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0H
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1, J_STRT_HALO = J_0H)

C**** initiallise work arrays
      MM0=0 ; MM1=0 ; MMX=0 ; UM0=0 ; VM0=0 ; UM1=0 ; VM1=0

C****
C**** Integrate Ocean Dynamics terms
C****
      OGEOZ_SV(:,:)=OGEOZ(:,:)

C**** Apply surface fluxes to ocean
      CALL GROUND_OC
         CALL CHECKO('GRNDOC')

C**** Apply ice/ocean and air/ocean stress to ocean
      CALL OSTRES
         CALL CHECKO('OSTRES')
         CALL TIMER (MNOW,MSURF)
         IF (MODD5S == 0) CALL DIAGCA (11)

C**** Calculate vertical diffusion
      CALL OCONV
         CALL CHECKO('OCONV ')

C**** Apply bottom and coastal drags
      if (OBottom_drag  == 1) CALL OBDRAG
      if (OCoastal_drag == 1) CALL OCOAST
         CALL TIMER (MNOW,MSGSO)

      IDACC(11) = IDACC(11) + 1

C****
C**** Integrate Ocean Dynamics
C****
C**** initialize summed mass fluxes
!$OMP PARALLEL DO  PRIVATE(L)
      DO L=1,LMO
        SMU(:,:,L) = 0. ; SMV(:,:,L) = 0. ; SMW(:,:,L) = 0.
      END DO
!$OMP END PARALLEL DO

      CALL OVTOM (MMI,UM0,VM0)

      CALL OPGF0

      NS=NDYNO
C**** Initial Forward step,  QMX = QM0 + DT*F(Q0)
      CALL OFLUX  (NS,MMI,.FALSE.)
      CALL OADVM (MM1,MMI,DTOFS)
C     CALL STADVM(MM1,MUST,DTOFS,.False.)
      CALL OADVV (UM1,VM1,UM0,VM0,DTOFS)
      CALL OPGF  (UM1,VM1,DTOFS)
C     CALL STPGF (MUST1,MUST,DTOFS)
      CALL OMTOV (MM1,UM1,VM1)
C**** Initial Backward step,  QM1 = QM0 + DT*F(Q0)
      CALL OFLUX  (NS,MMI,.FALSE.)
      CALL OADVM (MM1,MMI,DTO)
C     CALL STADVM(MM1,MUST,DTO,.False.)
      CALL OADVV (UM1,VM1,UM0,VM0,DTO)
      CALL OPGF  (UM1,VM1,DTO)
C     CALL STPGF (MUST1,MUST,DTO)
      CALL OMTOV (MM1,UM1,VM1)
C**** First even leap frog step,  Q2 = Q0 + 2*DT*F(Q1)
      CALL OFLUX  (NS,MMI,.TRUE.)
      CALL OADVM (MM0,MMI,DTOLF)
C     CALL STADVM(MM0,MUST1,DTOLF,.True.)
      CALL OADVV (UM0,VM0,UM0,VM0,DTOLF)
      CALL OPGF  (UM0,VM0,DTOLF)
C     CALL STPGF (MUST,MUST,DTOLF)
      CALL OMTOV (MM0,UM0,VM0)
      NS=NS-1
C**** Odd leap frog step,  Q3 = Q1 + 2*DT*F(Q2)
  420 Continue
      CALL OFLUX  (NS,MM1,.FALSE.)
      CALL OADVM (MM1,MM1,DTOLF)
C     CALL STADVM(MM1,MUST,DTOLF,.False.)
      CALL OADVV (UM1,VM1,UM1,VM1,DTOLF)
      CALL OPGF  (UM1,VM1,DTOLF)
C     CALL STPGF (MUST1,MUST1,DTOLF)
      CALL OMTOV (MM1,UM1,VM1)
      NS=NS-1
C**** Even leap frog step,  Q4 = Q2 + 2*DT*F(Q3)
      CALL OFLUX  (NS,MM0,.TRUE.)
      CALL OADVM (MM0,MM0,DTOLF)
C     CALL STADVM(MM0,MUST1,DTOLF,.True.)
      CALL OADVV (UM0,VM0,UM0,VM0,DTOLF)
      CALL OPGF  (UM0,VM0,DTOLF)
C     CALL STPGF (MUST,MUST,DTOLF)
      CALL OMTOV (MM0,UM0,VM0)
C**** Check for end of leap frog time scheme
      NS=NS-1
      IF(NS.GT.1)  GO TO 420
C     if(j_0 ==  1) UO(IVSP,1 ,:) = VOSP(:) ! not needed if Mod(IM,4)=0
      if(j_1 == JM) UO(IVNP,JM,:) = VONP(:) ! not needed if Mod(IM,4)=0
!$OMP PARALLEL DO  PRIVATE(L)
      DO L=1,LMO
        OIJL(:,:,L,IJL_MFU) = OIJL(:,:,L,IJL_MFU) + SMU(:,:,L)
        OIJL(:,:,L,IJL_MFV) = OIJL(:,:,L,IJL_MFV) + SMV(:,:,L)
        OIJL(:,:,L,IJL_MFW) = OIJL(:,:,L,IJL_MFW) + SMW(:,:,L)
      END DO
!$OMP END PARALLEL DO

C**** Advection of Potential Enthalpy and Salt
      CALL OADVT (G0M,GXMO,GYMO,GZMO,DTOLF,.FALSE.
     *        ,OIJL(1,J_0H,1,IJL_GFLX))
      CALL OADVT (S0M,SXMO,SYMO,SZMO,DTOLF,.TRUE.
     *        ,OIJL(1,J_0H,1,IJL_SFLX))

#ifdef TRACERS_OCEAN
      DO N=1,NTM
        CALL OADVT(TRMO(1,J_0H,1,N),TXMO(1,J_0H,1,N)
     *       ,TYMO(1,J_0H,1,N),TZMO(1,J_0H,1,N),DTOLF,t_qlimit(n)
     *       ,TOIJL(1,J_0H,1,TOIJL_TFLX,N))
      END DO
#endif
        CALL CHECKO ('OADVT ')
!$OMP PARALLEL DO PRIVATE(L)
        DO L=1,LMO
          OIJL(:,:,L,IJL_MO)  = OIJL(:,:,L,IJL_MO) +  MO(:,:,L)
          OIJL(:,:,L,IJL_G0M) = OIJL(:,:,L,IJL_G0M) + G0M(:,:,L)
          OIJL(:,:,L,IJL_S0M) = OIJL(:,:,L,IJL_S0M) + S0M(:,:,L)
        END DO
!$OMP END PARALLEL DO
        DO J=J_0,J_1
          DO I=1,IMAXJ(J)
            IF (FOCEAN(I,J).gt.0) THEN
              OIJ(I,J,IJ_SSH) = OIJ(I,J,IJ_SSH) + OGEOZ(I,J)
              OIJ(I,J,IJ_PB)  = OIJ(I,J,IJ_PB)  +
     +             (OPBOT(I,J)-ZE(LMM(I,J))*RHOWS*GRAV)
            END IF
          END DO
        END DO

#ifdef TRACERS_OCEAN
        DO N=1,NTM
!$OMP PARALLEL DO PRIVATE(L)
          DO L=1,LMO
            TOIJL(:,:,L,TOIJL_CONC,N)=TOIJL(:,:,L,TOIJL_CONC,N)
     *           +TRMO(:,:,L,N)
          END DO
!$OMP END PARALLEL DO
        END DO
#endif
        CALL TIMER (MNOW,MDYNO)

        IF (MODD5S == 0) CALL DIAGCA (12)

C**** Apply Wajowicz horizontal diffusion to UO and VO ocean currents
      CALL ODIFF(DTS)
      CALL CHECKO ('ODIFF0')

C**** Apply GM + Redi tracer fluxes
      CALL GMKDIF
      CALL GMFEXP(G0M,GXMO,GYMO,GZMO,.FALSE.,OIJL(1,J_0H,1,IJL_GGMFL))
      CALL GMFEXP(S0M,SXMO,SYMO,SZMO,.TRUE. ,OIJL(1,J_0H,1,IJL_SGMFL))
#ifdef TRACERS_OCEAN
      DO N = 1,NTM
        CALL GMFEXP(TRMO(1,J_0H,1,N),TXMO(1,J_0H,1,N),TYMO(1,J_0H,1,N),
     *    TZMO(1,J_0H,1,N),t_qlimit(n),TOIJL(1,J_0H,1,TOIJL_GMFL,N))
      END DO
#endif
      CALL CHECKO ('GMDIFF')
      CALL TIMER (MNOW,MSGSO)

c????          Non-parallelized parts : straits
c     straits: mo, G0M,Gx-zMO,S0M,Sx-zMO,TRMO,Tx-zMO,opress (ocean)

      call gather_ocean_straits (2)

      IF(AM_I_ROOT()) THEN
C****
C**** Acceleration and advection of tracers through ocean straits
C****
        CALL STPGF
        CALL STCONV
        CALL STBDRA
          CALL CHECKO_serial ('STBDRA')
        CALL STADV
          CALL CHECKO_serial ('STADV0')
      END IF

      call scatter_ocean_straits (2)
      call BCAST_straits (0)
c????      end of serialized part
        CALL CHECKO ('STADV ')
C**** remove STADVI since it is not really consistent with ICEDYN
c      CALL STADVI
c        CALL CHECKO ('STADVI')
#ifdef TRACERS_OCEAN
      CALL OC_TDECAY

#ifdef TRACERS_AGE_OCEAN
      CALL OCN_TR_AGE
#endif
#endif
        CALL TIMER (MNOW,MSGSO)
      CALL TOC2SST
      RETURN

      END SUBROUTINE OCEANS

      SUBROUTINE init_OCEAN(iniOCEAN,istart)
!@sum init_OCEAN initiallises ocean variables
!@auth Original Development Team
!@ver  1.0
      USE FILEMANAGER, only : openunit,closeunit
      USE TIMINGS, only : timing,ntimeacc
      USE PARAM
      USE CONSTANT, only : twopi,radius,by3,grav,rhow
      USE MODEL_COM, only : dtsrc,kocean
      USE OCEAN, only : im,jm,lmo,focean,lmm
     *     ,lmu,lmv,hatmo,hocean,ze,mo,g0m,gxmo,gymo,gzmo,s0m,sxmo
     *     ,symo,szmo,uo,vo,dxypo,ogeoz,dts,dtolf,dto,dtofs,mdyno,msgso
     *     ,ndyno,imaxj,ogeoz_sv,bydts,lmo_min,j1o
     *     ,OBottom_drag,OCoastal_drag,oc_salt_mean
#ifdef TRACERS_OCEAN
     *     ,oc_tracer_mean,ntm
#endif
      USE OCEANRES, only : dZO
      USE OCFUNC, only : vgsp,tgsp,hgsp,agsp,bgsp,cgs
      USE SW2OCEAN, only : init_solar
      USE FLUXES, only : ogeoza, uosurf, vosurf

      USE DOMAIN_DECOMP, only : grid,get

      IMPLICIT NONE
      INTEGER I,J,L,N,iu_OIC,iu_OFTAB,IP1,IM1,LMIJ,I1,J1,I2,J2
     *     ,iu_TOPO,II,JJ,flagij
      REAL*4, DIMENSION(IM,JM,LMO):: MO4,G0M4,S0M4,GZM4,SZM4
      CHARACTER*80 TITLE
      REAL*8 FJEQ,SM,SG0,SGZ,SS0,SSZ
      LOGICAL, INTENT(IN) :: iniOCEAN
      INTEGER, INTENT(IN) :: istart
      LOGICAL :: iniStraits
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S
      LOGICAL :: HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1
     *      ,J_STRT_SKP  = J_0S, J_STOP_SKP  = J_1S
     *      ,HAVE_NORTH_POLE = HAVE_NORTH_POLE)

C****
C**** Check that KOCEAN is set correctly
C****
      IF (KOCEAN == 0) THEN
        call stop_model(
     &       "Must have KOCEAN > 0 for interactive ocean runs",255)
      END IF
C****
C**** Select drag options
C****
      call sync_param("OBottom_drag",OBottom_drag)
      call sync_param("OCoastal_drag",OCoastal_drag)

C**** define initial condition options for global mean
      call sync_param("oc_salt_mean",oc_salt_mean)
#ifdef TRACERS_OCEAN
      call sync_param("oc_tracer_mean",oc_tracer_mean,ntm)
#endif

C****
C**** set up time steps from atmospheric model
C****
      call sync_param("DTO",DTO)

      DTS=DTSRC
      BYDTS=1d0/DTS
      NDYNO=2*NINT(.5*DTS/DTO)
      DTO=DTS/NDYNO
      DTOLF=2.*DTO
      DTOFS=2.*DTO*BY3
C**** Set up timing indexes
      CALL SET_TIMER(" OCEAN DYNAM",MDYNO)
      CALL SET_TIMER(" OCEAN PHYS.",MSGSO)
C****
C**** Arrays needed each ocean model run
C****
      CALL GEOMO

C**** Calculate ZE 

      ZE(0) = 0d0
      DO L = 1,LMO
        ZE(L) = ZE(L-1) + dZO(L) 
      END DO

C**** Read in table function for specific volume
      CALL openunit("OFTAB",iu_OFTAB,.TRUE.,.TRUE.)
      READ  (iu_OFTAB) TITLE,VGSP
      WRITE (6,*) 'Read from unit ',iu_OFTAB,': ',TITLE
      READ  (iu_OFTAB) TITLE,TGSP
      WRITE (6,*) 'Read from unit ',iu_OFTAB,': ',TITLE
      READ  (iu_OFTAB) TITLE,CGS
      WRITE (6,*) 'Read from unit ',iu_OFTAB,': ',TITLE
      READ  (iu_OFTAB) TITLE,HGSP
      WRITE (6,*) 'Read from unit ',iu_OFTAB,': ',TITLE
      READ  (iu_OFTAB) TITLE,AGSP
      WRITE (6,*) 'Read from unit ',iu_OFTAB,': ',TITLE
      READ  (iu_OFTAB) TITLE,BGSP
      WRITE (6,*) 'Read from unit ',iu_OFTAB,': ',TITLE
      call closeunit(iu_OFTAB)

C**** READ IN LANDMASKS AND TOPOGRAPHIC DATA
      call openunit("TOPO_OC",iu_TOPO,.true.,.true.)
      CALL READT (iu_TOPO,0,FOCEAN,IM*JM,FOCEAN,1) ! Ocean fraction
      CALL READT (iu_TOPO,0,HATMO ,IM*JM,HATMO ,4) ! Atmo. Topography
      CALL READT (iu_TOPO,0,HOCEAN,IM*JM,HOCEAN,1) ! Ocean depths
      call closeunit(iu_TOPO)

C**** Calculate J1O = least J with some ocean
      Do 130 J=1,JM
  130 If (Sum(FOCEAN(:,J)) > 0)  GoTo 140
  140 J1O = J
      write(6,*) "Minimum J with some ocean:",J1O
C**** Fix to 4 for temporary consistency
      J1O=4
      write(6,*) "Fixed J1O:",J1O

C**** Calculate LMM and modify HOCEAN
c      call sync_param("LMO_min",LMO_min)
      DO 170 J=1,JM    ! global arrays (fixed from now on)
      DO 170 I=1,IM
      LMM(I,J) =0
      IF(FOCEAN(I,J).LE.0.)  GO TO 170
      DO 150 L=LMO_min,LMO-1
  150 IF(HATMO(I,J)+HOCEAN(I,J) .le. 5d-1*(ZE(L)+ZE(L+1)))  GO TO 160
C     L=LMO
  160 LMM(I,J)=L
      HOCEAN(I,J) = -HATMO(I,J) + ZE(L)
  170 CONTINUE
C**** Calculate LMU
      I=IM
      DO 180 J=1,JM    ! global arrays (fixed from now on)
      DO 180 IP1=1,IM
      LMU(I,J) = MIN(LMM(I,J),LMM(IP1,J))
  180 I=IP1
C**** Calculate LMV
      DO 190 J=1,JM-1  ! global arrays (fixed from now on)
      DO 190 I=1,IM
 190  LMV(I,J) = MIN(LMM(I,J),LMM(I,J+1))
C****
      IF(iniOCEAN) THEN
C**** Initialize a run from ocean initial conditions
C???? For starters, let all processes read the IC
      CALL openunit("OIC",iu_OIC,.TRUE.,.TRUE.)
      READ  (iu_OIC,ERR=820) TITLE,MO4,G0M4,GZM4,S0M4,SZM4
      call closeunit(iu_OIC)
      WRITE (6,*) 'OIC read from unit ',iu_OIC,': ',TITLE
C**** Calculate layer mass from column mass and check for consistency
      DO 313 J=J_0,J_1
      DO 313 I=1,IM
      LMIJ=LMM(I,J)
      DO 311 L=1,LMIJ
  311 MO(I,J,L) = MO4(I,J,L)
      DO 312 L=LMIJ+1,LMO
  312 MO(I,J,L) = 0.
C**** if there is a problem try nearest neighbour
      IF((LMM(I,J).GT.0).AND.(ABS(MO(I,J,1)/ZE(1)-1d3).GT.5d1)) THEN
        WRITE (6,931) I,J,LMIJ,MO(I,J,1),ZE(1)
        II=0 ; JJ=0 ; flagij=0

        do j1=0,2
          do i1=0,4
            if(flagij.eq.0) then
              if(j1.eq.0) jj=j
              if(j1.eq.1) jj=j-1
              if(j1.eq.2) jj=j+1
              
              if(i1.eq.0) ii=i
              if(i1.eq.1) ii=i-1
              if(i1.eq.2) ii=i+1
              if(i1.eq.3) ii=i-2
              if(i1.eq.4) ii=i+2
              if(i1.eq.5) ii=i-3
              if(i1.eq.6) ii=i+3
              if(ii.gt.im) ii=ii-im
              if(ii.lt.1) ii=ii+im
              if(jj.gt.jm) then
                jj=jm
                if(ii.le.im/2) ii=ii+im/2
                if(ii.gt.im/2) ii=ii-im/2
              endif
              if(jj.lt.1) then
                jj=1
                if(ii.le.im/2) ii=ii+im/2
                if(ii.gt.im/2) ii=ii-im/2
              endif
              IF ((MO4(II,JJ,1).gt.0) .and. (LMM(II,JJ).ge.LMM(I,J))) 
     *             flagij=1
            endif
          enddo
        enddo
        IF (flagij.ne.0) THEN
          MO(I,J,1:LMM(I,J))=MO4(II,JJ,1:LMM(I,J))
          G0M4(I,J,1:LMM(I,J))=G0M4(II,JJ,1:LMM(I,J))
          S0M4(I,J,1:LMM(I,J))=S0M4(II,JJ,1:LMM(I,J))
          GZM4(I,J,1:LMM(I,J))=GZM4(II,JJ,1:LMM(I,J))
          SZM4(I,J,1:LMM(I,J))=SZM4(II,JJ,1:LMM(I,J))
          WRITE (6,*) "Inconsistency at ",I,J,"fixed from :",II,JJ
        END IF
      END IF
  313 CONTINUE
C**** Initialize velocity field to zero
      UO=0
      VO=0
C**** Define mean value of mass, potential heat, and salinity at poles
      DO 370 L=1,LMO
      if(HAVE_NORTH_POLE) then ! average polar ocean fields
        J=JM
        SM  = 0.
        SG0 = 0.
        SGZ = 0.
        SS0 = 0.
        SSZ = 0.
        DO I=1,IM
          SM  = SM  +   MO(I,J,L)
          SG0 = SG0 + G0M4(I,J,L)
          SGZ = SGZ + GZM4(I,J,L)
          SS0 = SS0 + S0M4(I,J,L)
          SSZ = SSZ + SZM4(I,J,L)
        end do
        DO I=1,IM
          MO(I,J,L)   = SM /IM
          G0M4(I,J,L) = SG0/IM
          GZM4(I,J,L) = SGZ/IM
          S0M4(I,J,L) = SS0/IM
          SZM4(I,J,L) = SSZ/IM
        end do
      end if
C**** Define East-West horizontal gradients
      GXMO=0
      GYMO=0
      SXMO=0
      SYMO=0
      IM1=IM-1
      I=IM
      DO 345 J=J_0S,J_1S
      DO 345 IP1=1,IM
      IF(LMM(I  ,J).LT.L)  GO TO 344
      IF(LMM(IM1,J).GE.L)  GO TO 342
      IF(LMM(IP1,J).LT.L)  GO TO 344
      GXMO(I,J,L) = .5*(G0M4(IP1,J,L)-G0M4(I,J,L))
      SXMO(I,J,L) = .5*(S0M4(IP1,J,L)-S0M4(I,J,L))
      GO TO 344
  342 IF(LMM(IP1,J).GE.L)  GO TO 343
      GXMO(I,J,L) = .5*(G0M4(I,J,L)-G0M4(IM1,J,L))
      SXMO(I,J,L) = .5*(S0M4(I,J,L)-S0M4(IM1,J,L))
      GO TO 344
  343 GXMO(I,J,L) = .25*(G0M4(IP1,J,L)-G0M4(IM1,J,L))
      SXMO(I,J,L) = .25*(S0M4(IP1,J,L)-S0M4(IM1,J,L))
  344 IM1=I
  345 I=IP1
C**** Define North-South horizontal gradients
      DO 354 J=J_0S,J_1S
      DO 354 I=1,IM
      IF(LMM(I,J  ).LT.L)  GO TO 354
      IF(LMM(I,J-1).GE.L)  GO TO 352
      IF(LMM(I,J+1).LT.L)  GO TO 354
      GYMO(I,J,L) = .5*(G0M4(I,J+1,L)-G0M4(I,J,L))
      SYMO(I,J,L) = .5*(S0M4(I,J+1,L)-S0M4(I,J,L))
      GO TO 354
  352 IF(LMM(I,J+1).GE.L)  GO TO 353
      GYMO(I,J,L) = .5*(G0M4(I,J,L)-G0M4(I,J-1,L))
      SYMO(I,J,L) = .5*(S0M4(I,J,L)-S0M4(I,J-1,L))
      GO TO 354
  353 GYMO(I,J,L) = .25*(G0M4(I,J+1,L)-G0M4(I,J-1,L))
      SYMO(I,J,L) = .25*(S0M4(I,J+1,L)-S0M4(I,J-1,L))
  354 CONTINUE
C**** Multiply specific quantities by mass
      DO 360 J=J_0,J_1
      DO 360 I=1,IM
      G0M(I,J,L)  = G0M4(I,J,L)*(MO(I,J,L)*DXYPO(J))
      GXMO(I,J,L) = GXMO(I,J,L)*(MO(I,J,L)*DXYPO(J))
      GYMO(I,J,L) = GYMO(I,J,L)*(MO(I,J,L)*DXYPO(J))
      GZMO(I,J,L) = GZM4(I,J,L)*(MO(I,J,L)*DXYPO(J))
      S0M(I,J,L)  = S0M4(I,J,L)*(MO(I,J,L)*DXYPO(J))
      SXMO(I,J,L) = SXMO(I,J,L)*(MO(I,J,L)*DXYPO(J))
      SYMO(I,J,L) = SYMO(I,J,L)*(MO(I,J,L)*DXYPO(J))
  360 SZMO(I,J,L) = SZM4(I,J,L)*(MO(I,J,L)*DXYPO(J))
  370 CONTINUE
C**** Initiallise geopotential field (needed by KPP)
      OGEOZ = 0.
      OGEOZ_SV = 0.

      END IF

C**** Extend ocean data to added layers at bottom if necessary
      if (istart.gt.0 .and. lmo_min .gt. 1) then
        do j=j_0s,j_1s
        do i=1,im
          if (lmm(i,j) == lmo_min .and. MO(i,j,lmo_min) == 0.) then
            do l=2,lmo_min
              if (MO(i,j,l) == 0.) then
                MO(i,j,l)   = (ZE(L)-ZE(L-1))*RHOW*
     *                     (1.+S0M(i,j,l-1)/(MO(i,j,l-1)*DXYPO(J)))
                G0M(i,j,l)  = G0M(i,j,l-1) *(MO(i,j,l)/MO(i,j,l-1))
                GXMO(i,j,l) = GXMO(i,j,l-1)*(MO(i,j,l)/MO(i,j,l-1))
                GYMO(i,j,l) = GYMO(i,j,l-1)*(MO(i,j,l)/MO(i,j,l-1))
                GZMO(i,j,l) = 0.
                S0M(i,j,l)  = S0M(i,j,l-1) *(MO(i,j,l)/MO(i,j,l-1))
                SXMO(i,j,l) = SXMO(i,j,l-1)*(MO(i,j,l)/MO(i,j,l-1))
                SYMO(i,j,l) = SYMO(i,j,l-1)*(MO(i,j,l)/MO(i,j,l-1))
                SZMO(i,j,l) = 0.
              end if
            end do
          end if
        end do
        end do
      end if

C**** zero out unphysical values (that might have come from a
C**** restart file with different topography)
      if (istart.lt.9) then
        DO J=J_0S,J_1S
          DO I=1,IMAXJ(J)
            VO(I,J,LMV(I,J)+1:LMO)=0.
          END DO
        END DO
        DO J=J_0,J_1
          DO I=1,IMAXJ(J)
            UO(I,J,LMU(I,J)+1:LMO)=0.
            MO(I,J,LMM(I,J)+1:LMO)=0.
          END DO
        END DO
      end if

C**** Initialize straits arrays
      iniStraits=iniOCEAN.or.(istart.lt.9)
      call init_STRAITS(iniStraits)

C**** Adjust global mean salinity if required (only at first start up)
      if (istart.gt.0 .and. istart.lt.9 .and. oc_salt_mean.ne.-999.)
     *     call adjust_mean_salt

C**** Initialize solar radiation penetration arrays
      call init_solar

C**** Initialize ocean diagnostics
      call init_ODIAG

C**** Initialize KPP mixing scheme
      call kmixinit(ZE)

#ifdef TRACERS_OCEAN
C**** Set diagnostics for ocean tracers
      call init_tracer_ocean
#endif

C**** Initialize some info passed to atmsophere
      uosurf=0 ; vosurf=0. ; ogeoza=0.

C**** Set atmospheric surface variables
      IF (ISTART.gt.0) CALL TOC2SST

      RETURN
C**** Terminate because of improper start up
  820 call stop_model('init_OCEAN: Error reading ocean IC',255)
C****
  931 FORMAT ('0Inconsistency between LMM and M:',3I4,2F10.1)

      RETURN
C****
      END SUBROUTINE init_OCEAN

      SUBROUTINE daily_OCEAN(end_of_day)
!@sum  daily_OCEAN performs the daily tasks for the ocean module
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : sday
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: end_of_day

C**** Only do this at end of the day
      IF (end_of_day) THEN

C**** Add glacial melt from Antarctica and Greenland
        CALL GLMELT(SDAY)

C**** set gtemp arrays for ocean
        CALL TOC2SST
      END IF
C****
      RETURN
      END SUBROUTINE daily_OCEAN

      SUBROUTINE io_ocean(kunit,iaction,ioerr)
!@sum  io_ocean is a driver to ocean related io routines
!@auth Gavin Schmidt
!@ver  1.0
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR

      call io_ocdyn  (kunit,iaction,ioerr)
      call io_straits(kunit,iaction,ioerr)

      RETURN
C****
      END SUBROUTINE io_ocean

      SUBROUTINE io_ocdyn(kunit,iaction,ioerr)
!@sum  io_ocdyn reads and writes ocean arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,irsficno,irsfic
     *     ,irsficnt,irerun,lhead
      USE OCEAN

      USE DOMAIN_DECOMP, only : AM_I_ROOT

      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "OCDYN01"
#ifdef TRACERS_OCEAN
!@var TRHEADER Character string label for individual records
      CHARACTER*80 :: TRHEADER, TRMODULE_HEADER = "TROCDYN02"

      write (TRMODULE_HEADER(lhead+1:80),'(a13,i3,a1,i3,a)')
     *     'R8 dim(im,jm,',LMO,',',NTM,'):TRMO,TX,TY,TZ'
#endif

      write (MODULE_HEADER(lhead+1:80),'(a13,i2,a)') 'R8 dim(im,jm,',
     *   LMO,'):M,U,V,G0,GX,GY,GZ,S0,SX,SY,SZ, OGZ,OGZSV'

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        call gather_ocean(0) ! mo,uo,vo,g0m,gx-z,s0m,sx-z,ogz's,tr
        if(AM_I_ROOT()) then
          WRITE (kunit,err=10) MODULE_HEADER,MO_glob,UO_glob,VO_glob
     *     ,G0M_glob,GXMO_glob,GYMO_glob,GZMO_glob
     *     ,S0M_glob,SXMO_glob,SYMO_glob,SZMO_glob
     *     ,OGEOZ_glob,OGEOZ_SV_glob
#ifdef TRACERS_OCEAN
          WRITE (kunit,err=10) TRMODULE_HEADER
     *     ,TRMO_glob,TXMO_glob,TYMO_glob,TZMO_glob
#endif
        end if
      CASE (IOREAD:)            ! input from restart file
        SELECT CASE (IACTION)
          CASE (IRSFICNO)   ! initial conditions (no ocean data)
            if(AM_I_ROOT()) READ (kunit)
          CASE (ioread,irerun,irsfic) ! restarts
           if(AM_I_ROOT()) then
            READ (kunit,err=10) HEADER,MO_glob,UO_glob,VO_glob
     *        ,G0M_glob,GXMO_glob,GYMO_glob,GZMO_glob
     *        ,S0M_glob,SXMO_glob,SYMO_glob,SZMO_glob
     *        ,OGEOZ_glob,OGEOZ_SV_glob
            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",HEADER
     *             ,MODULE_HEADER
              GO TO 10
            END IF
#ifdef TRACERS_OCEAN
            READ (kunit,err=10) TRHEADER
     *        ,TRMO_glob,TXMO_glob,TYMO_glob,TZMO_glob
            IF (TRHEADER(1:LHEAD).NE.TRMODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",TRHEADER
     *             ,TRMODULE_HEADER
              GO TO 10
            END IF
#endif
           end if
           call scatter_ocean (0)  ! mo,uo,vo,g0m,x-z,s0m,x-z,ogz's,tr
          CASE (irsficnt) ! restarts (never any tracer data)
           if(AM_I_ROOT()) then
            READ (kunit,err=10) HEADER,MO_glob,UO_glob,VO_glob
     *        ,G0M_glob,GXMO_glob,GYMO_glob,GZMO_glob
     *        ,S0M_glob,SXMO_glob,SYMO_glob,SZMO_glob
     *        ,OGEOZ_glob,OGEOZ_SV_glob
            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
              PRINT*,"Discrepancy in module version ",HEADER
     *             ,MODULE_HEADER
              GO TO 10
            END IF
           end if
           call scatter_ocean (-1) ! mo,uo,vo,g0m,gx-z,s0m,sx-z,ogz's
           return
        END SELECT
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
C****
      END SUBROUTINE io_ocdyn

      SUBROUTINE CHECKO_serial(SUBR)
!@sum  CHECKO Checks whether Ocean variables are reasonable (serial version)
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : byrt3,teeny
      USE MODEL_COM, only : qcheck
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : ntm, trname, t_qlimit
#endif
      USE OCEAN, only : im,jm,lmo,dxypo,focean,imaxj, lmm, mo=>mo_glob
     *  ,g0m=>g0m_glob, gxmo=>gxmo_glob,gymo=>gymo_glob,gzmo=>gzmo_glob
     *  ,s0m=>s0m_glob, sxmo=>sxmo_glob,symo=>symo_glob,szmo=>szmo_glob
     *  ,uo=>uo_glob, vo=>vo_glob
#ifdef TRACERS_OCEAN
     *  ,trmo=>trmo_glob, txmo=>txmo_glob, tymo=>tymo_glob,
     *   tzmo=>tzmo_glob
#endif
      IMPLICIT NONE
      REAL*8 SALIM,GO1,SO1,relerr,errmax,temgs
      LOGICAL QCHECKO
      INTEGER I,J,L,n,imax,jmax,lmax
!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

C**** Check for NaN/INF in ocean data
      IF (QCHECK) THEN
      CALL CHECK3B(MO  ,IM,1,JM,JM,LMO,SUBR,'mo ')
      CALL CHECK3B(G0M ,IM,1,JM,JM,LMO,SUBR,'g0m')
      CALL CHECK3B(GXMO,IM,1,JM,JM,LMO,SUBR,'gxm')
      CALL CHECK3B(GYMO,IM,1,JM,JM,LMO,SUBR,'gym')
      CALL CHECK3B(GZMO,IM,1,JM,JM,LMO,SUBR,'gzm')
      CALL CHECK3B(S0M ,IM,1,JM,JM,LMO,SUBR,'s0m')
      CALL CHECK3B(SXMO,IM,1,JM,JM,LMO,SUBR,'sxm')
      CALL CHECK3B(SYMO,IM,1,JM,JM,LMO,SUBR,'sym')
      CALL CHECK3B(SZMO,IM,1,JM,JM,LMO,SUBR,'szm')
      CALL CHECK3B(UO  ,IM,1,JM,JM,LMO,SUBR,'uo ')
      CALL CHECK3B(VO  ,IM,1,JM,JM,LMO,SUBR,'vo ')
#ifdef TRACERS_OCEAN
      CALL CHECK4B(TRMO,IM,1,JM,JM,LMO,NTM,SUBR,'tzm')
#endif

C**** Check for variables out of bounds
      QCHECKO=.FALSE.
      DO J=1,JM
      DO I=1,IMAXJ(J)
        IF(FOCEAN(I,J).gt.0.) THEN
C**** Check potential specific enthalpy/salinity
          DO L=1,LMM(I,J)
          GO1 = G0M(I,J,L)/(MO(I,J,L)*DXYPO(J))
          SO1 = S0M(I,J,L)/(MO(I,J,L)*DXYPO(J))
          IF(GO1.lt.-10000. .or. GO1.gt.200000.) THEN
            WRITE (6,*) 'After ',SUBR,': I,J,L,GO=',I,J,L,GO1,TEMGS(GO1
     *           ,SO1)
            IF (GO1.lt.-20000. .or. GO1.gt.200000.) QCHECKO=.TRUE.
          END IF
          IF(SO1.gt.0.045 .or. SO1.lt.0.) THEN
            WRITE (6,*) 'After ',SUBR,': I,J,L,SO=',I,J,L,1d3*SO1
            IF ((SO1.gt.0.05 .or. SO1.lt.0.) .and. .not. (I == 47.and.
     *           J == 30)) QCHECKO=.TRUE.
          END IF
          END DO
C**** Check all ocean currents
          DO L = 1,LMO
            IF(ABS(UO(I,J,L)).gt.7. .or. ABS(VO(I,J,L)).gt.5) THEN
              WRITE (6,*) 'After ',SUBR,': I,J,L,UO,VO=',I,J,L,UO(I,J,L)
     *             ,VO(I,J,L)
              QCHECKO=.TRUE.
            END IF
          END DO
C**** Check first layer ocean mass
          IF(MO(I,J,1).lt.2000. .or. MO(I,J,1).gt.20000.) THEN
!!          IF (I == 47.and.(J == 33.or.J == 34)) GOTO 230 ! not Caspian
            WRITE (6,*) 'After ',SUBR,': I,J,MO=',I,J,MO(I,J,1)
!!          QCHECKO=.TRUE.
          END IF
C**** Check ocean salinity in each eighth box for the first layer
 230      SALIM = .045
!!        IF(JM == 46 .and. I == 47 .and. J == 30) GOTO 240 ! not Persian Gulf   !.048
          IF(.5*ABS(SXMO(I,J,1))+.5*ABS(SYMO(I,J,1))+BYRT3*ABS(SZMO(I,J
     *         ,1)).lt.S0M(I,J,1) .and.(.5*ABS(SXMO(I,J,1))+.5
     *         *ABS(SYMO(I,J,1))+BYRT3*ABS(SZMO(I,J,1))+S0M(I,J,1))
     *         /(MO(I,J,1)*DXYPO(J)).lt.SALIM)  GO TO 240
          WRITE (6,*) 'After ',SUBR,': I,J,S0,SX,SY,SZ=',I,J,
     *         1d3*S0M(I,J,1)/(MO(I,J,1)*DXYPO(J)),1d3*SXMO(I,J,1)/(MO(I
     *         ,J,1)*DXYPO(J)),1d3*SYMO(I,J,1)/(MO(I,J,1)*DXYPO(J)),1d3
     *         *SZMO(I,J,1)/(MO(I,J,1)*DXYPO(J))
!!        QCHECKO=.TRUE.
 240      CONTINUE
        END IF
      END DO
      END DO

#ifdef TRACERS_OCEAN
      do n=1,ntm
C**** Check for negative tracers
        if (t_qlimit(n)) then
        do l=1,lmo
        do j=1,jm
        do i=1,imaxj(j)
          if (l.le.lmm(i,j)) then
            if (trmo(i,j,l,n).lt.0) then
              write(6,*) "Neg Tracer in ocean after ",subr,i,j,l,
     *             trname(n),trmo(i,j,l,n)
              QCHECKO=.true.
            end if
          end if
        end do
        end do
        end do
        end if
C**** Check conservation of water tracers in ocean
        if (trname(n) == 'Water') then
          errmax = 0. ; imax=1 ; jmax=1 ; lmax=1
          do l=1,lmo
          do j=1,jm
          do i=1,imaxj(j)
            if (l.le.lmm(i,j)) then
              relerr=max(
     *             abs(trmo(i,j,l,n)-mo(i,j,l)*dxypo(j)+s0m(i,j,l)),
     *             abs(txmo(i,j,l,n)+sxmo(i,j,l)),
     *             abs(tymo(i,j,l,n)+symo(i,j,l)),
     *             abs(tzmo(i,j,l,n)+szmo(i,j,l)))/
     *             (mo(i,j,l)*dxypo(j)-s0m(i,j,l))
            if (relerr.gt.errmax) then
              imax=i ; jmax=j ; lmax=l ; errmax=relerr
            end if
            end if
          end do
          end do
          end do
          print*,"Relative error in ocean fresh water mass after ",subr
     *         ,":",imax,jmax,lmax,errmax,trmo(imax,jmax,lmax,n),mo(imax
     *         ,jmax,lmax)*dxypo(jmax)-s0m(imax,jmax,lmax),txmo(imax
     *         ,jmax,lmax,n),-sxmo(imax,jmax,lmax),tymo(imax,jmax,lmax,n
     *         ),-symo(imax,jmax ,lmax),tzmo(imax,jmax,lmax,n),
     *         -szmo(imax,jmax,lmax)
        end if
      end do
#endif

      IF (QCHECKO)
     &     call stop_model("QCHECKO: Ocean Variables out of bounds",255)

      END IF
C****
      CALL CHECKOST(SUBR)
C****
      END SUBROUTINE CHECKO_serial

      SUBROUTINE CHECKO(SUBR)
!@sum  CHECKO Checks whether Ocean variables are reasonable (parallel version)
!@auth Original Development Team
!@ver  1.0
      USE CONSTANT, only : byrt3,teeny
      USE MODEL_COM, only : qcheck
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : ntm, trname, t_qlimit
#endif
      USE OCEAN
      USE DOMAIN_DECOMP, only : grid, GET, AM_I_ROOT
      IMPLICIT NONE
      REAL*8 SALIM,GO1,SO1,relerr,errmax,temgs
      LOGICAL QCHECKO
      INTEGER I,J,L,n,imax,jmax,lmax
!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

c**** Extract domain decomposition info
      INTEGER :: J_0S, J_0, J_1, J_0H, J_1H, JM_loc
      CALL GET(grid, J_STRT_SKP = J_0S, J_STRT = J_0, J_STOP = J_1,
     *   J_STRT_HALO = J_0H, J_STOP_HALO = J_1H)

C**** Check for NaN/INF in ocean data
      IF (QCHECK) THEN
      JM_loc = J_1H - J_0H + 1
      CALL CHECK3B(MO  ,IM,J_0H,J_1H,JM,LMO,SUBR,'mo ')
      CALL CHECK3B(G0M ,IM,J_0H,J_1H,JM,LMO,SUBR,'g0m')
      CALL CHECK3B(GXMO,IM,J_0H,J_1H,JM,LMO,SUBR,'gxm')
      CALL CHECK3B(GYMO,IM,J_0H,J_1H,JM,LMO,SUBR,'gym')
      CALL CHECK3B(GZMO,IM,J_0H,J_1H,JM,LMO,SUBR,'gzm')
      CALL CHECK3B(S0M ,IM,J_0H,J_1H,JM,LMO,SUBR,'s0m')
      CALL CHECK3B(SXMO,IM,J_0H,J_1H,JM,LMO,SUBR,'sxm')
      CALL CHECK3B(SYMO,IM,J_0H,J_1H,JM,LMO,SUBR,'sym')
      CALL CHECK3B(SZMO,IM,J_0H,J_1H,JM,LMO,SUBR,'szm')
      CALL CHECK3B(UO  ,IM,J_0H,J_1H,JM,LMO,SUBR,'uo ')
      CALL CHECK3B(VO  ,IM,J_0H,J_1H,JM,LMO,SUBR,'vo ')
#ifdef TRACERS_OCEAN
      CALL CHECK4B(TRMO,IM,J_0H,J_1H,JM,LMO,NTM,SUBR,'trm')
      CALL CHECK4B(TXMO,IM,J_0H,J_1H,JM,LMO,NTM,SUBR,'txm')
      CALL CHECK4B(TYMO,IM,J_0H,J_1H,JM,LMO,NTM,SUBR,'tym')
      CALL CHECK4B(TZMO,IM,J_0H,J_1H,JM,LMO,NTM,SUBR,'tzm')
#endif

C**** Check for variables out of bounds
      QCHECKO=.FALSE.
      DO J=J_0,J_1
      DO I=1,IMAXJ(J)
        IF(FOCEAN(I,J).gt.0.) THEN
C**** Check potential specific enthalpy/salinity
          DO L=1,LMM(I,J)
          GO1 = G0M(I,J,L)/(MO(I,J,L)*DXYPO(J))
          SO1 = S0M(I,J,L)/(MO(I,J,L)*DXYPO(J))
          IF(GO1.lt.-10000. .or. GO1.gt.200000.) THEN
            WRITE (6,*) 'After ',SUBR,': I,J,L,GO=',I,J,L,GO1,TEMGS(GO1
     *           ,SO1)
            IF (GO1.lt.-20000. .or. GO1.gt.200000.) QCHECKO=.TRUE.
          END IF
          IF(SO1.gt.0.045 .or. SO1.lt.0.) THEN
            WRITE (6,*) 'After ',SUBR,': I,J,L,SO=',I,J,L,1d3*SO1
            IF ((SO1.gt.0.05 .or. SO1.lt.0.) .and. .not. (I == 47.and.
     *           J == 30)) QCHECKO=.TRUE.
          END IF
          END DO
C**** Check all ocean currents
          DO L = 1,LMO
            IF(ABS(UO(I,J,L)).gt.7. .or. ABS(VO(I,J,L)).gt.5) THEN
              WRITE (6,*) 'After ',SUBR,': I,J,L,UO,VO=',I,J,L,UO(I,J,L)
     *             ,VO(I,J,L)
              QCHECKO=.TRUE.
            END IF
          END DO
C**** Check first layer ocean mass
          IF(MO(I,J,1).lt.2000. .or. MO(I,J,1).gt.20000.) THEN
!!          IF (I == 47.and.(J == 33.or.J == 34)) GOTO 230 ! not Caspian
            WRITE (6,*) 'After ',SUBR,': I,J,MO=',I,J,MO(I,J,1)
!!          QCHECKO=.TRUE.
          END IF
C**** Check ocean salinity in each eighth box for the first layer
 230      SALIM = .045
!!        IF(JM == 46 .and. I == 47 .and. J == 30) GOTO 240 ! not Persian Gulf   !.048
          IF(.5*ABS(SXMO(I,J,1))+.5*ABS(SYMO(I,J,1))+BYRT3*ABS(SZMO(I,J
     *         ,1)).lt.S0M(I,J,1) .and.(.5*ABS(SXMO(I,J,1))+.5
     *         *ABS(SYMO(I,J,1))+BYRT3*ABS(SZMO(I,J,1))+S0M(I,J,1))
     *         /(MO(I,J,1)*DXYPO(J)).lt.SALIM)  GO TO 240
          WRITE (6,*) 'After ',SUBR,': I,J,S0,SX,SY,SZ=',I,J,
     *         1d3*S0M(I,J,1)/(MO(I,J,1)*DXYPO(J)),1d3*SXMO(I,J,1)/(MO(I
     *         ,J,1)*DXYPO(J)),1d3*SYMO(I,J,1)/(MO(I,J,1)*DXYPO(J)),1d3
     *         *SZMO(I,J,1)/(MO(I,J,1)*DXYPO(J))
!!        QCHECKO=.TRUE.
 240      CONTINUE
        END IF
      END DO
      END DO

#ifdef TRACERS_OCEAN
      do n=1,ntm
C**** Check for negative tracers
        if (t_qlimit(n)) then
        do l=1,lmo
        do j=j_0,j_1
        do i=1,imaxj(j)
          if (l.le.lmm(i,j)) then
            if (trmo(i,j,l,n).lt.0) then
              write(6,*) "Neg Tracer in ocean after ",subr,i,j,l,
     *             trname(n),trmo(i,j,l,n)
              QCHECKO=.true.
            end if
          end if
        end do
        end do
        end do
        end if
C**** Check conservation of water tracers in ocean
        if (trname(n) == 'Water') then
          errmax = 0. ; imax=1 ; jmax=1 ; lmax=1
          do l=1,lmo
          do j=j_0,j_1
          do i=1,imaxj(j)
            if (l.le.lmm(i,j)) then
              relerr=max(
     *             abs(trmo(i,j,l,n)-mo(i,j,l)*dxypo(j)+s0m(i,j,l)),
     *             abs(txmo(i,j,l,n)+sxmo(i,j,l)),
     *             abs(tymo(i,j,l,n)+symo(i,j,l)),
     *             abs(tzmo(i,j,l,n)+szmo(i,j,l)))/
     *             (mo(i,j,l)*dxypo(j)-s0m(i,j,l))
            if (relerr.gt.errmax) then
              imax=i ; jmax=j ; lmax=l ; errmax=relerr
            end if
            end if
          end do
          end do
          end do
          print*,"Relative error in ocean fresh water mass after ",subr
     *         ,":",imax,jmax,lmax,errmax,trmo(imax,jmax,lmax,n),mo(imax
     *         ,jmax,lmax)*dxypo(jmax)-s0m(imax,jmax,lmax),txmo(imax
     *         ,jmax,lmax,n),-sxmo(imax,jmax,lmax),tymo(imax,jmax,lmax,n
     *         ),-symo(imax,jmax ,lmax),tzmo(imax,jmax,lmax,n),
     *         -szmo(imax,jmax,lmax)
        end if
      end do
#endif

      IF (QCHECKO)
     &     call stop_model("QCHECKO: Ocean Variables out of bounds",255)

      END IF
C****
      IF (AM_I_ROOT()) CALL CHECKOST(SUBR)
C****
      END SUBROUTINE CHECKO

      SUBROUTINE conserv_OKE(OKE)
!@sum  conserv_OKE calculates zonal ocean kinetic energy
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : shw
      USE OCEAN, only : im,jm,lmo,ivnp,fim,imaxj,focean,mo,uo,vo,lmm
      USE DOMAIN_DECOMP, only : GRID, GET, SOUTH, HALO_UPDATE
      IMPLICIT NONE
!@var OKE zonal ocean kinetic energy per unit area (J/m**2)
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: OKE
      INTEGER I,J,L,IP1

      INTEGER :: J_0S, J_1S
      LOGICAL :: HAVE_NORTH_POLE

      CALL GET(grid, J_STRT_SKP=J_0S,    J_STOP_SKP=J_1S,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)

      CALL HALO_UPDATE(grid, VO, FROM=SOUTH)
      OKE=0.
      DO J=J_0S,J_1S
        I=IM
        DO IP1=1,IM
          DO L=1,LMM(I,J)
            OKE(J) = OKE(J)+((MO(I,J,L)+MO(IP1,J,L))*UO(I,J,L)*UO(I,J,L)
     *      + MO(I,J,L)*(VO(I,J-1,L)*VO(I,J-1,L) + VO(I,J,L)*VO(I,J,L)))
          END DO
          I=IP1
        END DO
        OKE(J)=OKE(J)*0.25
      END DO
      if(HAVE_NORTH_POLE) then
        DO L=1,LMM(1,JM)
          OKE(JM) = OKE(JM) +
     +      MO(1,JM,L) * (1.5*IM*(UO(IM,JM,L)**2 + UO(IVNP,JM,L)**2) +
     +                    Sum(VO(:,JM-1,L)*VO(:,JM-1,L)) )
        END DO
        OKE(JM)= OKE(JM)*0.25
      end if
C****
      RETURN
      END SUBROUTINE conserv_OKE

      SUBROUTINE conserv_OCE(OCEANE)
!@sum  conserv_OCE calculates zonal ocean potential enthalpy(atmos grid)
!@auth Gavin Schmidt
!@ver  1.0
      USE GEOM, only : bydxyp
      USE OCEAN, only : im,jm,fim,imaxj,focean,g0m,lmm
      USE STRAITS, only : nmst,jst,g0mst
      USE DOMAIN_DECOMP, only : GRID, GET
      IMPLICIT NONE
!@var OCEANE zonal ocean potential enthalpy (J/m^2)
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: OCEANE
      INTEGER I,J,L,N

      INTEGER :: J_0, J_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      CALL GET(grid, J_STRT=J_0,    J_STOP=J_1,
     &  HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,HAVE_NORTH_POLE=HAVE_NORTH_POLE)

      DO J=J_0,J_1
        OCEANE(J)=0
        DO I=1,IMAXJ(J)
          DO L=1,LMM(I,J)
            OCEANE(J) = OCEANE(J) + G0M(I,J,L)*FOCEAN(I,J)*BYDXYP(J)
          END DO
        END DO
      END DO
      if(HAVE_SOUTH_POLE) OCEANE(1) =FIM*OCEANE(1)
      if(HAVE_NORTH_POLE) OCEANE(JM)=FIM*OCEANE(JM)
C**** include straits variables
      DO N=1,NMST
        J=JST(N,1)
       if(j.lt.j_0 .or. j.gt.j_1) cycle
        OCEANE(J)=OCEANE(J)+SUM(G0MST(:,N))*BYDXYP(J)
      END DO
C****
      RETURN
      END SUBROUTINE conserv_OCE

      SUBROUTINE conserv_OMS(OMASS)
!@sum  conserv_OMS calculates zonal ocean mass (on atmos grid)
!@auth Gavin Schmidt
!@ver  1.0
      USE GEOM, only : bydxyp
      USE OCEAN, only : im,jm,fim,imaxj,focean,mo,g0m,lmm,dxypo
      USE STRAITS, only : nmst,jst,mmst
      USE DOMAIN_DECOMP, only : GRID, GET, pack_data
      IMPLICIT NONE
!@var OMASS zonal ocean mass (kg/m^2)
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: OMASS
      REAL*8, DIMENSION(JM) :: OMSSV
      COMMON /OCCONS/OMSSV
      INTEGER I,J,L,N

      INTEGER :: J_0, J_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      CALL GET(grid, J_STRT=J_0,    J_STOP=J_1,
     &  HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,HAVE_NORTH_POLE=HAVE_NORTH_POLE)

      DO J=J_0,J_1
        OMASS(J)=0
        DO I=1,IMAXJ(J)
          DO L=1,LMM(I,J)
            OMASS(J) = OMASS(J) + MO(I,J,L)*FOCEAN(I,J)
          END DO
        END DO
        OMASS(J)=OMASS(J)*DXYPO(J)*BYDXYP(J)
      END DO
      if (HAVE_SOUTH_POLE) OMASS(1) =FIM*OMASS(1)
      if (HAVE_NORTH_POLE) OMASS(JM)=FIM*OMASS(JM)
C**** include straits variables
      DO N=1,NMST
        J=JST(N,1)
        if(j.lt.j_0 .or. j.gt.j_1) cycle
        OMASS(J)=OMASS(J)+SUM(MMST(:,N))*BYDXYP(J)
      END DO
C**** save mass for AM calculation
      OMSSV(J_0:J_1) = OMASS(J_0:J_1)
C****
      RETURN
      END SUBROUTINE conserv_OMS

      SUBROUTINE conserv_OAM(OAM)
!@sum  conserv_OAM calculates zonal ocean angular momentum
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : radius,omega
      USE OCEAN, only : im,jm,fim,imaxj,focean,mo,uo,cosm,cosq,lmu
      USE DOMAIN_DECOMP, only : GET, GRID
      IMPLICIT NONE
!@var OAM ocean angular momentum divided by area (kg/s)
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: OAM
      REAL*8, DIMENSION(JM) :: OMSSV
      COMMON /OCCONS/OMSSV
      INTEGER I,J,L,IP1
      REAL*8 UMILx2

      INTEGER :: J_0S, J_1S
      LOGICAL :: HAVE_NORTH_POLE

      CALL GET(grid, J_STRT_SKP=J_0S,    J_STOP_SKP=J_1S,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      OAM=0
      DO J=J_0S,J_1S
        UMILx2 = 0.
        I=IM
        DO IP1=1,IM
          DO L=1,LMU(I,J)
            UMILx2 = UMILx2 + UO(I,J,L)*(MO(I,J,L)+MO(IP1,J,L))
          END DO
          I=IP1
        END DO
c       OAM(J) = UMIL*COSPO(J) + OMSSV(J)*RADIUS*OMEGA*(COSVO(J-1)
c    *       *COSVO(J-1)+COSVO(J)*COSVO(J))
c       OAM(J)=0.5*RADIUS*OAM(J)
        OAM(J) = RADIUS*(OMSSV(J)*RADIUS*OMEGA*COSQ(J) +
     +                   .5*UMILx2*COSM(J))
      END DO

      if (HAVE_NORTH_POLE)
     *  OAM(JM) = RADIUS*OMSSV(JM)*RADIUS*OMEGA*COSQ(JM)
c     DO L=1,LMU(1,JM)
c       UMIL = UMIL + UO(1,JM,L)*MO(1,JM,L)
c     END DO
c     OAM(JM) = UMIL*COSPO(JM)*IM*2. +OMSSV(JM)*RADIUS*OMEGA
c    *     *COSVO(JM-1)*COSVO(JM-1)
c     OAM(JM)=0.5*RADIUS*OAM(JM)
C****
      RETURN
      END SUBROUTINE conserv_OAM

      SUBROUTINE conserv_OSL(OSALT)
!@sum  conserv_OSL calculates zonal ocean salt on atmos grid
!@auth Gavin Schmidt
!@ver  1.0
      USE GEOM, only : bydxyp
      USE OCEAN, only : im,jm,fim,imaxj,focean,mo,s0m,lmm
      USE STRAITS, only : nmst,jst,s0mst
      USE DOMAIN_DECOMP, only : GET, GRID
      IMPLICIT NONE
!@var OSALT zonal ocean salt (kg/m^2)
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: OSALT
      INTEGER I,J,L,N

      INTEGER :: J_0, J_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      CALL GET(grid, J_STRT=J_0,    J_STOP=J_1,
     &  HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,HAVE_NORTH_POLE=HAVE_NORTH_POLE)

      DO J=J_0,J_1
        OSALT(J)=0
        DO I=1,IMAXJ(J)
          DO L=1,LMM(I,J)
            OSALT(J) = OSALT(J) + S0M(I,J,L)*FOCEAN(I,J)*BYDXYP(J)
          END DO
        END DO
      END DO
      if (HAVE_SOUTH_POLE) OSALT(1) =FIM*OSALT(1)
      if (HAVE_NORTH_POLE) OSALT(JM)=FIM*OSALT(JM)
C**** include straits variables
      DO N=1,NMST
        J=JST(N,1)
        if(j.lt.j_0 .or. j.gt.j_1) cycle
        OSALT(J)=OSALT(J)+SUM(S0MST(:,N))*BYDXYP(J)
      END DO
C****
      RETURN
      END SUBROUTINE conserv_OSL

      Subroutine OVtoM (MM,UM,VM)
C****
C**** OVtoM converts density and velocity in concentration units
C**** to mass and momentum in mass units
C**** Input:  MO (kg/m^2), UO (m/s), VO (m/s)
C**** Define: VOSP from UO(IVSP,1), VONP from UO(IVNP,JM)
C**** Output: MM (kg), UM (kg*m/s), VM (kg*m/s)
C****
      Use OCEAN, Only: IM,JM,LMO, IVSP,IVNP, DXYP=>DXYPO, COSU,SINU,
     *                 COSI=>COSIC,SINI=>SINIC, MO,UO,VO, VONP
      Use DOMAIN_DECOMP, Only: GRID, HALO_UPDATE, NORTH
      Implicit None
      Real*8,Intent(Out),
     *  Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) :: MM,UM,VM
      Integer*4 I,J,L, J1,JN,J1P,JNP,JNQ,JNR
      Logical*4 QSP,QNP
C****
C**** Extract domain decomposition band parameters
C**** ASSUME THAT J1 NEVER EQUALS 2 AND JN NEVER EQUALS JM-1
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      J1P = Max(J1,2)       !    2      5     JM-3   Exclude SP
      JNP = Min(JN,JM-1)    !    4      8     JM-1   Exclude NP
      JNQ = Min(JN,JM-2)    !    4      8     JM-2   Exclude NP,NP-1
      JNR = Min(JN+1,JM-1)  !    5      9     JM-1   Halo, Exclude NP
C     QSP = J1==1           !    T      F      F
      QNP = JN==JM          !    F      F      T
C****
      Call HALO_UPDATE (GRID, MO, FROM=NORTH)
C****
!$OMP ParallelDo   Private (I,J,L)
      DO 50 L=1,LMO
C**** Define VOSP and VONP
C     If (QSP)  VOSP(L) = UO(IVSP,1 ,L)
      If (QNP)  VONP(L) = UO(IVNP,JM,L)
C**** Convert density to mass
      DO 10 J=J1,JNR
   10 MM(:,J,L) = MO(:,J,L)*DXYP(J)
C     If (QSP)  MM(1,1 ,L) = MO(1,1 ,L)*DXYP(1)
      If (QNP)  MM(1,JM,L) = MO(1,JM,L)*DXYP(JM)
C**** Convert U velocity to U momentum
      DO 30 J=J1P,JNP
      DO 20 I=1,IM-1
   20 UM(I ,J,L) = UO(I ,J,L)*(MM(I,J,L)+MM(I+1,J,L))*.5
   30 UM(IM,J,L) = UO(IM,J,L)*(MM(1,J,L)+MM(IM ,J,L))*.5
C**** Convert V velocity to V momentum
      DO 40 J=J1P,JNQ
   40 VM(:,J,L) = VO(:,J,L)*(MM(:,J,L)+MM(:,J+1,L))*.5
C     If (QSP)  VM(:, 1  ,L) = VO(:, 1  ,L)*(MM(:, 2  ,L)+MM(1, 1,L))*.5
      If (QNP)  VM(:,JM-1,L) = VO(:,JM-1,L)*(MM(:,JM-1,L)+MM(1,JM,L))*.5
C**** Convert U and V velocity to U and V momentum at poles
C     If (QSP)  UM(IM,1 ,L) = UO(IM,1 ,L)*MM(1,1 ,L)
      If (QNP)  UM(IM,JM,L) = UO(IM,JM,L)*MM(1,JM,L)
C     If (QSP)  UM(IVSP,1 ,L) =   VOSP(L)*MM(1,1 ,L)
      If (QNP)  UM(IVNP,JM,L) =   VONP(L)*MM(1,JM,L)
C**** Define MO, UO and VO at poles
C     If (QSP)  Then
C       MO(:,1 ,L) = MO(1 ,1 ,L)
C       UO(:,1 ,L) = UO(IM,1 ,L)*COSU(:) - VOSP(L)*SINU(:)
C       VO(:,0 ,L) = VOSP(L)*COSI(:) + UO(IM,1 ,L)*SINI(:)  ;  EndIf
      If (QNP)  Then
        MO(:,JM,L) = MO(1 ,JM,L)
        UO(:,JM,L) = UO(IM,JM,L)*COSU(:) + VONP(L)*SINU(:)
        VO(:,JM,L) = VONP(L)*COSI(:) - UO(IM,JM,L)*SINI(:)  ;  EndIf
   50 Continue
      Return
      EndSubroutine OVtoM

      Subroutine OMtoV (MM,UM,VM)
C****
C**** OMtov converts mass and momentum in mass units
C**** to density and velocity in concentration units
C**** Input:  MM (kg), UM (kg*m/s), VM (kg*m/s)
C**** Output: MO (kg/m^2), UO (m/s), VO (m/s)
C**** Define: VOSP from UM(IVSP,1)/MM, VONP from UM(IVNP,JM)/MM
C****
      Use OCEAN, Only: IM,JM,LMO,IVSP,IVNP,
     *                  LMOM=>LMM, LMOU=>LMU, LMOV=>LMV,
     *                  COSU,SINU, COSI=>COSIC,SINI=>SINIC,
     *                  zDXYP=>BYDXYPO, MO,UO,VO, VONP
      Use DOMAIN_DECOMP, Only: GRID, HALO_UPDATE, NORTH
      Implicit None
      Real*8, !!! Intent(IN), (except for HALO_UPDATEs)
     *  Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) :: MM,UM,VM
      Integer*4 I,J,L, J1,JN,J1P,JNP,JNQ
      Logical*4 QSP,QNP
C****
C**** Extract domain decomposition band parameters
C**** ASSUME THAT J1 NEVER EQUALS 2 AND JN NEVER EQUALS JM-1
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      J1P = Max(J1,2)       !    2      5     JM-3   Exclude SP
      JNP = Min(JN,JM-1)    !    4      8     JM-1   Exclude NP
      JNQ = Min(JN,JM-2)    !    4      8     JM-2   Exclude NP,NP-1
C     QSP = J1==1           !    T      F      F
      QNP = JN==JM          !    F      F      T
C****
      Call HALO_UPDATE (GRID, MM, FROM=NORTH)
C**** Convert mass to density
!$OMP ParallelDo   Private (I,J,L)
      Do 10 J=J1P,JNP
      Do 10 I=1,IM
      Do 10 L=1,LMOM(I,J)
   10 MO(I,J,L) = MM(I,J,L)*zDXYP(J)
C     If (QSP)  Then
C       Do L=1,LMOM(1,1)
C  11   MO(:,1,L) = MM(1,1,L)*zDXYP(1)  ;  EndIf
      If (QNP)  Then
        Do 12 L=1,LMOM(1,JM)
   12   MO(:,JM,L) = MM(1,JM,L)*zDXYP(JM)  ;  EndIf
C**** Convert U momentum to U velocity
!$OMP ParallelDo   Private (I,J,L)
      Do 30 J=J1P,JNP
      Do 20 I=1,IM-1
      Do 20 L=1,LMOU(I,J)
   20 UO(I,J,L) = UM(I,J,L)*2/(MM(I,J,L)+MM(I+1,J,L))
      Do 30 L=1,LMOU(IM,J)
   30 UO(IM,J,L) = UM(IM,J,L)*2/(MM(IM,J,L)+MM(1,J,L))
C**** Convert V momentum to V velocity
!$OMP ParallelDo   Private (I,J,L)
      Do 40 J=J1P,JNQ
      Do 40 I=1,IM
      Do 40 L=1,LMOV(I,J)
   40 VO(I,J,L) = VM(I,J,L)*2/(MM(I,J,L)+MM(I,J+1,L))
C     If (QSP)  Then
C       Do 41 I=1,IM
C       Do 41 L=1,LMOV(I,1)
C  41   VO(I,1,L) = VM(I,1,L)*2/(MM(1,1,L)+MM(I,2,L))  ;  EndIf
      If (QNP)  Then
        Do 42 I=1,IM
        Do 42 L=1,LMOV(I,JM-1)
   42   VO(I,JM-1,L) = VM(I,JM-1,L)*2/(MM(1,JM,L)+MM(I,JM-1,L))  ; EndIf
C**** Convert momentum to velocity at south pole, define UO and VO
C     If (QSP)  Then
C       Do 50 L=1,LMOM(1,1)
C       UO(IM,1,L) = UM(IM,1,L)/MM(1,1,L)
C       VOSP(L)  = UM(IVSP,1,L)/MM(1,1,L)
C       UO(:,1,L) = UO(IM,1,L)*COSU(:) - VOSP(L)*SINU(:)
C       VO(:,0,L) = VOSP(L)*COSI(:) + UO(IM,1,L)*SINI(:)
C  50   Continue  ;  EndIf
C**** Convert momentum to velocity at north pole, define UO and VO
      If (QNP)  Then
        Do 60 L=1,LMOM(1,JM)
        UO(IM,JM,L) = UM(IM,JM,L)/MM(1,JM,L)
        VONP(L)   = UM(IVNP,JM,L)/MM(1,JM,L)
        UO(:,JM,L) = UO(IM,JM,L)*COSU(:) + VONP(L)*SINU(:)
        VO(:,JM,L) = VONP(L)*COSI(:) - UO(IM,JM,L)*SINI(:)
   60   Continue  ;  EndIf
      Return
      EndSubroutine OMtoV

      Subroutine OFLUX (NS,MM,QEVEN)
C****
C**** OFLUX calculates the fluid fluxes
C**** Input:  M (kg/m^2), U (m/s), V (m/s)
C**** Output: MU (kg/s) = DY * U * M
C****         MV (kg/s) = DX * V * M
C****         MW (kg/s) = MW(L-1) + CONV-CONVs*dZO/ZOE + MM-MMs*dZO/ZOE
C****       CONV (kg/s) = MU-MU + MV-MV
C****
      Use OCEAN, Only: IM,JM,LMO, LMOM=>LMM, ZOE=>ZE,dZO, DTO,
     *                 DXYP=>DXYPO,DYP=>DYPO,DXV=>DXVO, MO,UO,VO
      Use OCEAN_DYN, Only: MU,MV,MW, SMU,SMV,SMW, CONV
      Use DOMAIN_DECOMP, Only: GRID, HALO_UPDATE, NORTH,SOUTH
      Implicit None
      Integer*4,Intent(In) :: NS     !  decrementing leap-frog step
      Logical*4,Intent(In) :: QEVEN  !  whether called from even step
      Real*8   ,Intent(In),
     *  Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) :: MM
C****
      Real*8    zNSxDTO, MVS,dMVS(IM),dMVSm, MVN,dMVN(IM),dMVNm,
     *          CONVs,MMs
      Integer*4 I,J,L,LM,IMAX, J1,JN,J1P,J1H,JNP,JNQ
      Logical*4 QSP,QNP
C****
C**** Extract domain decomposition band parameters
C**** ASSUME THAT J1 NEVER EQUALS 2 AND JN NEVER EQUALS JM-1
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      J1P = Max(J1,2)       !    2      5     JM-3   Exclude SP
      J1H = Max(J1-1,1)     !    1      4     JM-4   Halo minimum
      JNP = Min(JN,JM-1)    !    4      8     JM-1   Exclude NP
      JNQ = Min(JN,JM-2)    !    4      8     JM-2   Exclude NP,NP-1
C     QSP = J1==1           !    T      F      F
      QNP = JN==JM          !    F      F      T
C****
      Call HALO_UPDATE (GRID, MO, FROM=NORTH+SOUTH)
      Call HALO_UPDATE (GRID, VO, FROM=SOUTH)
      zNSxDTO = 1 / (NS*DTO)
C****
C**** Compute fluid fluxes for the C grid
C****
C**** Smooth the West-East velocity near the poles
!$OMP ParallelDo   Private (J,L)
      Do 110 L=1,LMO
      Do 110 J=J1P,JNP
  110 MU(:,J,L) = UO(:,J,L)
      Call OPFIL (MU)
C**** Compute MU, the West-East mass flux, at non-polar points
!$OMP ParallelDo   Private (I,J,L, MVS,dMVS,dMVSm, MVN,dMVN,dMVNm)
      DO 430 L=1,LMO
      DO 130 J=J1P,JNP
      DO 120 I=1,IM-1
  120 MU(I ,J,L) = .5*DYP(J)*MU(I ,J,L)*(MO(I,J,L)+MO(I+1,J,L))
  130 MU(IM,J,L) = .5*DYP(J)*MU(IM,J,L)*(MO(1,J,L)+MO(IM ,J,L))
C**** Compute MV, the South-North mass flux
      DO 210 J=J1H,JNP
  210 MV(:,J,L) = .5*DXV(J)*VO(:,J,L)*(MO(:,J,L)+MO(:,J+1,L))
C****
C**** Compute MU so that CONV is identical for each polar triangle
C****
C     If (QSP) then
C       MVS = Sum(MV(:,1,L)) / IM
C       dMVS(1) = 0
C       Do 310 I=2,IM
C 310   dMVS(I) = dMVS(I-1) + (MV(I,1,L)-MVS)
C       dMVSm   = Sum(dMVS(:)) / IM
C       MU(:,1,L) = dMVSm - dMVS(:)  ;  EndIf
      If (QNP) then
        MVN = Sum(MV(:,JM-1,L)) / IM
        dMVN(1) = 0
        Do 320 I=2,IM
  320   dMVN(I) = dMVN(I-1) + (MV(I,JM-1,L)-MVN)
        dMVNm   = Sum(dMVN(:)) / IM
        MU(:,JM,L) = dMVN(:) - dMVNm  ;  EndIf
C****
C**** Compute horizontal fluid convergence at non-polar points
C****
      Do 420 J=J1P,JNP
      Do 410 I=2,IM
  410 CONV(I,J,L) = MU(I-1,J,L)-MU(I,J,L) + (MV(I,J-1,L)-MV(I,J,L))
  420 CONV(1,J,L) = MU(IM ,J,L)-MU(1,J,L) + (MV(1,J-1,L)-MV(1,J,L))
C     If (QSP)  CONV(1,1 ,L) = - MVS
      If (QNP)  CONV(1,JM,L) =   MVN
  430 Continue
C****
C**** Compute vertically integrated column convergence and mass
C****
!$OMP ParallelDo   Private (I,J,L, IMAX,LM, CONVs,MMs)
      Do 630 J=J1,JN
      IMAX=IM  ;  If(J==1.or.J==JM) IMAX=1
      Do 630 I=1,IMAX
      LM=1
      If (LMOM(I,J) <= 1)  GoTo 620
      LM=LMOM(I,J)
      CONVs = Sum(CONV(I,J,1:LM)) / ZOE(LM)
      MMs   = Sum(  MM(I,J,1:LM)) / ZOE(LM)
C****
C**** Compute MW, the downward fluid flux
C****
      MW(I,J,1) = CONV(I,J,1) - CONVs*dZO(1) +
     +           (  MM(I,J,1) -   MMs*dZO(1)) * zNSxDTO
      Do 610 L=2,LM-1
  610 MW(I,J,L) = CONV(I,J,L) - CONVs*dZO(L) + MW(I,J,L-1) +
     +           (  MM(I,J,L) -   MMs*dZO(L)) * zNSxDTO
  620 MW(I,J,LM:LMO-1) = 0
  630 Continue
C****
C**** Sum mass fluxes to be used for advection of tracers
C****
      If (.not.QEVEN)  Return
!$OMP ParallelDo   Private (L)
      Do 710 L=1,LMO
      SMU(:,J1 :JN ,L) = SMU(:,J1 :JN ,L) + MU(:,J1 :JN ,L)
      SMV(:,J1H:JNP,L) = SMV(:,J1H:JNP,L) + MV(:,J1H:JNP,L)
      If (L==LMO)  GoTo 710
      SMW(:,J1P:JNP,L) = SMW(:,J1P:JNP,L) + MW(:,J1P:JNP,L)
  710 Continue
C     If (QSP)  SMW(1,1 ,:) = SMW(1,1 ,:) + MW(1,1 ,:)
      If (QNP)  SMW(1,JM,:) = SMW(1,JM,:) + MW(1,JM,:)
      Return
      EndSubroutine OFLUX

      Subroutine OPFIL (X)
C****
C**** OPFIL smoothes X in the zonal direction by reducing coefficients
C**** of its Fourier series for high wave numbers near the poles.
C**** The ocean polar filter here works with AVR4X5LD.Z12.gas1 .
C****
      Use CONSTANT, Only: TWOPI
      Use OCEAN, Only: IM,JM,LMO,J1O, DLON, DXP=>DXPO, DYP=>DYPO,
     *  JMPF=>J40S  !  greatest J in SH where polar filter is applied
      Use FILEMANAGER, Only: OPENUNIT, CLOSEUNIT
      Use DOMAIN_DECOMP, Only: GRID
      Implicit None
C****
      Real*8,Intent(InOut) ::
     *  X(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO)
      Integer*4,Parameter :: IMz2=IM/2  !  highest possible wave number
      Integer*4,Save ::
     *  JMEX,    !  exclude unfiltered latitudes, store J=JM-1 in JMEX
     *  NBASM,   !  maximum number of ocean basins at any latitude
     *  INDM,    !  number of entries in array REDUCO
     *  NMIN(JM) !  minimum wave number for Fourier smoothing
      Real*8,Save :: SMOOTH(IMz2,JM)  !  multiplies Fourier coefficient
C****
      Integer*2,Save,Allocatable,Dimension(:,:) ::
     *  NBAS    !  number of ocean basins for given layer and latitude
      Integer*2,Save,Allocatable,Dimension(:,:,:) ::
     *  IMINm1, !  western most cell in basin minus 1
     *  IWIDm1  !  number of cells in basin minus 1
      Integer*4,Save,Allocatable,Dimension(:,:) ::
     *  INDEX   !  index to concatenated matricies
      Real*4,Save,Allocatable,Dimension(:) ::
     *  REDUCO  !  concatenation of reduction matricies
C****
      Character*80 TITLE
      Integer*4 I,I1,INDX,IWm2, J,JA,JX, K,L, N,NB, IU_AVR
      Real*8 AN(0:IMz2), !  Fourier cosine coefficients
     *       BN(0:IMz2), !  Fourier sine coefficients
     *          Y(IM*2), !  original copy of X that wraps around IDL
     *       DRAT, REDUC
C****
      Integer*4,Save :: IFIRST=1, J1,JN,J1P,JNP
      If (IFIRST == 0)  GoTo 100
      IFIRST = 0
      JMEX = 2*JMPF-1
C**** Extract domain decomposition band parameters
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      J1P = Max(J1,2)       !    2      5     JM-3   Exclude SP
      JNP = Min(JN,JM-1)    !    4      8     JM-1   Exclude NP
C**** Calculate SMOOTHing factor for longitudes without land cells
      Do 30 J=J1O,JM-1
      DRAT = DXP(J)/DYP(3)
      Do 20 N=IMz2,1,-1
      SMOOTH(N,J) = IMz2*DRAT/N
   20 If (SMOOTH(N,J) >= 1)  GoTo 30
C     N = 0
   30 NMIN(J) = N+1
C**** Read in reduction contribution matrices from disk
      Call OPENUNIT ('AVR',IU_AVR,.True.,.True.)
      Read (IU_AVR) TITLE,NBASM,INDM
      Allocate (IMINm1(NBASM,LMO,J1O:JMEX), NBAS(LMO,J1O:JMEX),
     *          IWIDm1(NBASM,LMO,J1O:JMEX), INDEX(IM,2:JMPF),
     *          REDUCO(INDM))
      Read (IU_AVR) TITLE,NBAS,IMINm1,IWIDm1,INDEX,REDUCO
      Call CLOSEUNIT (IU_AVR)
      Write (6,*) 'Read from unit',IU_AVR,': ',TITLE
 100  CONTINUE
C****
C**** Loop over J and L.  JX = eXclude unfiltered latitudes
C****                     JA = Absolute latitude
C****
!$OMP ParallelDo   Private (I,I1,INDX,IWm2, J,JA,JX, K,L, N,NB,
!$OMP*                      AN,BN, REDUC,Y)
      Do 410 J=Max(J1O,J1P),JNP
      JX=J  ;  If(J > JMPF) JX=J+2*JMPF-JM
      JA=J  ;  If(J > JMPF) JA=JM+1-J
      If (JA > JMPF)  GoTo 410  !  skip latitudes J=JMPF+1,JM-JMPF
      Do 400 L=1,LMO
      If (IWIDm1(1,L,JX) >= IM)  GoTo 300
C****
C**** Land cells exist at this latitude and layer, loop over ocean
C**** basins.
C****
      Do 270 NB=1,NBAS(L,JX)
      I1   = IMINm1(NB,L,JX) + 1
      IWm2 = IWIDm1(NB,L,JX) - 1
      INDX = INDEX(IWm2+1,JA)
      If (I1+IWm2 > IM)  GoTo 200
C**** Ocean basin does not wrap around the IDL.
C**** Copy X to temporary array Y and filter X in place.
      Do 110 I=I1,I1+IWm2
  110 Y(I) = X(I,J,L)
      Do 140 I=I1,I1+IWm2
      REDUC = 0
      Do 130 K=I1,I1+IWm2
      INDX=INDX+1
  130 REDUC = REDUC + REDUCO(INDX)*Y(K)
  140 X(I,J,L) = X(I,J,L) - REDUC
      GoTo 270
C**** Ocean basin wraps around the IDL.
C**** Copy X to temporary array Y and filter X in place.
  200 Do 210 I=I1,IM
  210 Y(I) = X(I,J,L)
      Do 220 I=IM+1,I1+IWm2
  220 Y(I) = X(I-IM,J,L)
      Do 240 I=I1,IM
      REDUC = 0
      Do 230 K=I1,I1+IWm2
      INDX=INDX+1
  230 REDUC = REDUC + REDUCO(INDX)*Y(K)
  240 X(I,J,L) = X(I,J,L) - REDUC
      Do 260 I=1,I1+IWm2-IM
      REDUC = 0
      Do 250 K=I1,I1+IWm2
      INDX=INDX+1
  250 REDUC = REDUC + REDUCO(INDX)*Y(K)
  260 X(I,J,L) = X(I,J,L) - REDUC
  270 Continue
      GoTo 400
C****
C**** No land cells at this latitude and layer,
C**** perform standard polar filter
C****
  300 Call FFT (X(1,J,L),AN,BN)
      Do 310 N=NMIN(J),IMz2-1
      AN(N) = AN(N)*SMOOTH(N,J)
  310 BN(N) = BN(N)*SMOOTH(N,J)
      AN(IMz2) = AN(IMz2)*SMOOTH(IMz2,J)
      Call FFTI (AN,BN,X(1,J,L))
C****
  400 Continue
  410 Continue
      Return
      EndSubroutine OPFIL

      Subroutine OADVM (MM2,MM0,DT)
C****
C**** OADVM calculates the updated mass
C**** Input:  MM0 (kg), DT (s), CONV (kg/s), MW (kg/s)
C**** Output: MM2 (kg) = MM0 + DT*[CONV + MW(L-1) - MW(L)]
C****
      Use OCEAN, Only: IM,JM,LMO, LMOM=>LMM
      Use OCEAN_DYN, Only: MW, CONV
      Use DOMAIN_DECOMP, Only: GRID
      Implicit None
      Real*8,Intent(Out):: MM2(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO)
      Real*8,Intent(In) :: MM0(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO)
     *                    ,DT
      Integer*4 I,J,L,LM,IMAX, J1,JN
C****
C**** Extract domain decomposition band parameters
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
C****
C**** Compute the new mass MM2
C****
!$OMP ParallelDo   Private (I,J, IMAX,LM)
      Do 20 J=J1,JN
      IMAX=IM  ;  If(J==1.or.J==JM) IMAX=1
      Do 20 I=1,IMAX
      If (LMOM(I,J)==0)  GoTo 20
      If (LMOM(I,J)==1)  Then
        MM2(I,J,1) = MM0(I,J,1) + DT*CONV(I,J,1)
        GoTo 20  ;  EndIf
      LM = LMOM(I,J)
      MM2(I,J,1) = MM0(I,J,1) + DT*(CONV(I,J,1)-MW(I,J,1))
      Do 10 L=2,LM-1
   10 MM2(I,J,L) = MM0(I,J,L) + DT*(CONV(I,J,L) + MW(I,J,L-1)-MW(I,J,L))
      MM2(I,J,LM) = MM0(I,J,LM) + DT*(CONV(I,J,LM) + MW(I,J,LM-1))
   20 Continue
      Return
      EndSubroutine OADVM

      Subroutine OADVV (UM2,VM2,UM0,VM0,DT1)
C****
C**** OADVV advects oceanic momentum (with coriolis force)
C**** Input:  MO (kg/m^2), UO (m/s), VO (m/s) = from odd solution
C****         MU (kg/s), MV (kg/s), MW (kg/s) = fluid fluxes
C**** Output: UM2 (kg*m/s) = UM0 + DT*(MU*U-MU*U + MV*U-MV*U + M*CM*V)
C****         VM2 (kg*m/s) = VM0 + DT*(MU*V-MU*V + MV*V-MV*V - M*CM*U)
C****
      Use CONSTANT, Only: RADIUS,OMEGA
      Use OCEAN, Only: IM,JM,LMO, IVSP,IVNP, DXV=>DXVO,
     *                 COSU,SINU, SINxY,TANxY, MO,UO,VO
      Use OCEAN_DYN, Only: MU,MV,MW
      Use DOMAIN_DECOMP, Only: GRID, HALO_UPDATE, NORTH,SOUTH
      Implicit None
      Real*8,Intent(Out),
     *  Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) :: UM2,VM2
      Real*8,Intent(In),
     *  Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) :: UM0,VM0
      Real*8,Intent(In) :: DT1
C****
      Real*8 DT2,DT4,DT6,DT12
      Real*8,Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) ::
     *       DUM,DVM
      Real*8,Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     *       UMVc,UMVw,UMVe
      Real*8 UMU,VMV,VMUc,VMUs,VMUn, FLUX, dUMSP,dVMSP,dUMNP,dVMNP,
     *       USINYJ(IM),UTANUJ(IM)
      Integer*4 I,J,L, J1,JN,J1P,J1H,JNP,JNR
      Logical*4 QSP,QNP
C****
C**** Extract domain decomposition band parameters
C**** ASSUME THAT J1 NEVER EQUALS 2 AND JN NEVER EQUALS JM-1
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      J1P = Max(J1,2)       !    2      5     JM-3   Exclude SP
      J1H = Max(J1-1,1)     !    1      4     JM-4   Halo minimum
      JNP = Min(JN,JM-1)    !    4      8     JM-1   Exclude NP
      JNR = Min(JN+1,JM-1)  !    5      9     JM-1   Halo, Exclude NP
C     QSP = J1==1           !    T      F      F
      QNP = JN==JM          !    F      F      T
C****
C     Call HALO_UPDATE (GRID, MO, FROM=NORTH)        !  copied in OFLUX
      Call HALO_UPDATE (GRID, UO, FROM=NORTH+SOUTH)
      Call HALO_UPDATE (GRID, VO, FROM=NORTH)
C     Call HALO_UPDATE (GRID, VO, FROM=SOUTH)        !  copied in OFLUX
      Call HALO_UPDATE (GRID, MU, FROM=NORTH)
      Call HALO_UPDATE (GRID, MV, FROM=NORTH)
C     Call HALO_UPDATE (GRID, MV, FROM=SOUTH)        !  recalc in OFLUX
      Call HALO_UPDATE (GRID, MW, FROM=NORTH)
C****
      DT2  = DT1/2
      DT4  = DT1/4
      DT6  = DT1/6
      Dt12 = DT1/12
C****
C**** Horizontal advection of momentum
C****
!$OMP ParallelDo   Private (I,J,L, UMU, UMVc,UMVw,UMVe, USINYJ,
!$OMP*                             VMV, VMUc,VMUs,VMUn, UTANUJ)
      Do 480 L=1,LMO
C**** Zero out momentum changes
      DUM(:,:,L) = 0
      DVM(:,:,L) = 0
C**** Contribution of eastward flux to U momentum
      Do 310 J=J1,JN
      UMU = DT4*(UO(1,J,L)+UO(IM,J,L))*(MU(1,J,L)+MU(IM,J,L))
      DUM(1 ,J,L) = DUM(1 ,J,L) + UMU
      DUM(IM,J,L) = DUM(IM,J,L) - UMU
      Do 310 I=2,IM
      UMU = DT4*(UO(I,J,L)+UO(I-1,J,L))*(MU(I,J,L)+MU(I-1,J,L))
      DUM(I  ,J,L) = DUM(I  ,J,L) + UMU
  310 DUM(I-1,J,L) = DUM(I-1,J,L) - UMU
C**** Contribution of northward center to U momentum
      Do 320 J=J1H,JNP
      UMVc(IM,J) = DT6*(UO(IM,J,L)+UO(IM,J+1,L))*(MV(IM,J,L)+MV(1,J,L))
      DUM(IM,J  ,L) = DUM(IM,J  ,L) - UMVc(IM,J)
      DUM(IM,J+1,L) = DUM(IM,J+1,L) + UMVc(IM,J)
      Do 320 I=1,IM-1
      UMVc(I,J) = DT6*(UO(I,J,L)+UO(I,J+1,L))*(MV(I,J,L)+MV(I+1,J,L))
      DUM(I,J  ,L) = DUM(I,J  ,L) - UMVc(I,J)
  320 DUM(I,J+1,L) = DUM(I,J+1,L) + UMVc(I,J)
C**** Contribution of northward corner fluxes to U momentum
      Do 330 J=J1H,JNP
      UMVe(1,J) = DT12*(UO(IM,J,L)+UO(1 ,J+1,L))*MV(1,J,L)
      UMVw(1,J) = DT12*(UO(1 ,J,L)+UO(IM,J+1,L))*MV(1,J,L)
      DUM(1 ,J  ,L) = DUM(1 ,J  ,L) - UMVw(1,J)
      DUM(IM,J  ,L) = DUM(IM,J  ,L) - UMVe(1,J)
      DUM(1 ,J+1,L) = DUM(1 ,J+1,L) + UMVe(1,J)
      DUM(IM,J+1,L) = DUM(IM,J+1,L) + UMVw(1,J)
      Do 330 I=2,IM
      UMVe(I,J) = DT12*(UO(I-1,J,L)+UO(I  ,J+1,L))*MV(I,J,L)
      UMVw(I,J) = DT12*(UO(I  ,J,L)+UO(I-1,J+1,L))*MV(I,J,L)
      DUM(I  ,J  ,L) = DUM(I  ,J  ,L) - UMVw(I,J)
      DUM(I-1,J  ,L) = DUM(I-1,J  ,L) - UMVe(I,J)
      DUM(I  ,J+1,L) = DUM(I  ,J+1,L) + UMVe(I,J)
  330 DUM(I-1,J+1,L) = DUM(I-1,J+1,L) + UMVw(I,J)
C**** Contribution of eastward center flux to V momentum
      Do 340 J=J1,JNP
      VMUc = DT6*(VO(IM,J,L)+VO(1,J,L))*(MU(IM,J,L)+MU(IM,J+1,L))
      DVM(IM,J,L) = DVM(IM,J,L) - VMUc
      DVM(1 ,J,L) = DVM(1 ,J,L) + VMUc
      Do 340 I=1,IM-1
      VMUc = DT6*(VO(I,J,L)+VO(I+1,J,L))*(MU(I,J,L)+MU(I,J+1,L))
      DVM(I  ,J,L) = DVM(I  ,J,L) - VMUc
  340 DVM(I+1,J,L) = DVM(I+1,J,L) + VMUc
C**** Contribution of eastward corner and northward fluxes to V momentum
C     If (QSP)  Then
C       J = 1
C       VMV  = DT4 *(VO(IM,1,L)+V0(IM,0,L))*MV(IM,1,L)
C       VMUn = DT12*(V0(IM,0,L)+VO(1 ,1,L))*MU(IM,1,L)
C       VMUs = DT12*(VO(IM,1,L)+V0(1 ,0,L))*MU(IM,1,L)
C       DVM(IM,1,L) = DVM(IM,1,L) - VMUs + VMV
C       DVM(1 ,1,L) = DVM(1 ,1,L) + VMUn
C       Do 350 I=1,IM-1
C       VMV  = DT4* (VO(I,1,L)+V0(I  ,0,L))*MV(I,1,L)
C       VMUn = DT12*(V0(I,0,L)+VO(I+1,1,L))*MU(I,1,L)
C       VMUs = DT12*(VO(I,1,L)+V0(I+1,0,L))*MU(I,1,L)
C       DVM(I  ,1,L) = DVM(I  ,1,L) - VMUs + VMV
C 350   DVM(I+1,1,L) = DVM(I+1,1,L) + VMUn  ;  EndIf
      Do 360 J=J1P,JNR
      VMV  = DT4 *(VO(IM,J  ,L)+VO(IM,J-1,L))*(MV(IM,J,L)+MV(IM,J-1,L))
      VMUn = DT12*(VO(IM,J-1,L)+VO(1 ,J  ,L))* MU(IM,J,L)
      VMUs = DT12*(VO(IM,J  ,L)+VO(1 ,J-1,L))* MU(IM,J,L)
      DVM(IM,J  ,L) = DVM(IM,J  ,L) - VMUs + VMV
      DVM(IM,J-1,L) = DVM(IM,J-1,L) - VMUn - VMV
      DVM(1 ,J  ,L) = DVM(1 ,J  ,L) + VMUn
      DVM(1 ,J-1,L) = DVM(1 ,J-1,L) + VMUs
      Do 360 I=1,IM-1
      VMV  = DT4* (VO(I,J  ,L)+VO(I  ,J-1,L))*(MV(I,J,L)+MV(I,J-1,L))
      VMUn = DT12*(VO(I,J-1,L)+VO(I+1,J  ,L))* MU(I,J,L)
      VMUs = DT12*(VO(I,J  ,L)+VO(I+1,J-1,L))* MU(I,J,L)
      DVM(I  ,J  ,L) = DVM(I  ,J  ,L) - VMUs + VMV
      DVM(I  ,J-1,L) = DVM(I  ,J-1,L) - VMUn - VMV
      DVM(I+1,J  ,L) = DVM(I+1,J  ,L) + VMUn
  360 DVM(I+1,J-1,L) = DVM(I+1,J-1,L) + VMUs
      If (QNP)  Then
C       J = JM
        VMV  = DT4 *(VO(IM,JM  ,L)+VO(IM,JM-1,L))*MV(IM,JM-1,L)
        VMUn = DT12*(VO(IM,JM-1,L)+VO(1 ,JM  ,L))*MU(IM,JM,L)
        VMUs = DT12*(VO(IM,JM  ,L)+VO(1 ,JM-1,L))*MU(IM,JM,L)
        DVM(IM,JM-1,L) = DVM(IM,JM-1,L) - VMUn - VMV
        DVM(1 ,JM-1,L) = DVM(1 ,JM-1,L) + VMUs
        Do 370 I=1,IM-1
        VMV  = DT4* (VO(I,JM  ,L)+VO(I  ,JM-1,L))*MV(I,JM-1,L)
        VMUn = DT12*(VO(I,JM-1,L)+VO(I+1,JM  ,L))*MU(I,JM,L)
        VMUs = DT12*(VO(I,JM  ,L)+VO(I+1,JM-1,L))*MU(I,JM,L)
        DVM(I  ,JM-1,L) = DVM(I  ,JM-1,L) - VMUn - VMV
  370   DVM(I+1,JM-1,L) = DVM(I+1,JM-1,L) + VMUs  ;  EndIf
C****
C**** Coriolis force and metric term
C****
C**** U component
C     If (QSP)  Then
C       Do 410 I=1,IM-1
C 410   dUM(I,1,L) = dUM(I,1,L) +
C    +    DT2*OMEGA*RADIUS*SINxY(1)*(MV(I,1,L)+MV(I+1,1,L)) +
C    +    TANxY(1)*(UMVw(I,1)+UMVc(I,1)+UMVe(I+1,1))*.5
C       dUM(IM,1,L) = dUM(IM,1,L) +
C    +    DT2*OMEGA*RADIUS*SINxY(1)*(MV(IM,1,L)+MV(1,1,L)) +
C    +    TANxY(1)*(UMVw(IM,1)+UMVc(IM,1)+UMVe(1,1))*.5  !  EndIf
      Do 420 J=J1P,JNP
      dUM(IM,J,L) = dUM(IM,J,L) + DT2*OMEGA*RADIUS*SINxY(J)*
     *    (MV(IM,J-1,L)+MV(1,J-1,L)+MV(IM,J,L)+MV(1,J,L)) +
     +  TANxY(J)*(UMVe(IM,J-1)+UMVc(IM,J-1)+UMVw(1,J-1) +
     +            UMVw(IM,J  )+UMVc(IM,J  )+UMVe(1,J  ))*.5
      Do 420 I=1,IM-1
  420 dUM(I,J,L) = dUM(I,J,L) + DT2*OMEGA*RADIUS*SINxY(J)*
     *    (MV(I,J-1,L)+MV(I+1,J-1,L)+MV(I,J,L)+MV(I+1,J,L)) +
     +  TANxY(J)*(UMVe(I,J-1)+UMVc(I,J-1)+UMVw(I+1,J-1) +
     +            UMVw(I,J  )+UMVc(I,J  )+UMVe(I+1,J  ))*.5
      If (QNP)  Then
        Do 430 I=1,IM-1
  430   dUM(I,JM,L) = dUM(I,JM,L) +
     +    DT2*OMEGA*RADIUS*SINxY(JM)*(MV(I,JM-1,L)+MV(I+1,JM-1,L)) +
     +    TANxY(JM)*(UMVe(I,JM-1)+UMVc(I,JM-1)+UMVw(I+1,JM-1))*.5
        dUM(IM,JM,L) = dUM(IM,JM,L) +
     +    DT2*OMEGA*RADIUS*SINxY(JM)*(MV(IM,JM-1,L)+MV(1,JM-1,L)) +
     +    TANxY(JM)*(UMVe(IM,JM-1)+UMVc(IM,JM-1)+UMVw(1,JM-1))*.5
        EndIf
C**** V component
      Do 450 J=J1,JNP
      Do 440 I=1,IM
      USINYJ(I) =  UO(I,J,L)*SINxY(J) + UO(I,J+1,L)*SINxY(J+1)
  440 UTANUJ(I) = (UO(I,J,L)*TANxY(J) + UO(I,J+1,L)*TANxY(J+1))*
     *            (UO(I,J,L)+UO(I,J+1,L))
      dVM(1,J,L) = dVM(1,J,L) - DT4*DXV(J)*(MO(1,J,L)+MO(1,J+1,L))*
     *  (OMEGA*RADIUS*(USINYJ(IM)+USINYJ(1)) +
     +  .25*(UTANUJ(IM)+UTANUJ(1)) - (TANxY(J)+TANxY(J+1))*
     *    (UO(IM,J,L)-UO(1,J,L))*(UO(IM,J+1,L)-UO(1,J+1,L))/12)
      Do 450 I=2,IM
  450 dVM(I,J,L) = dVM(I,J,L) - DT4*DXV(J)*(MO(I,J,L)+MO(I,J+1,L))*
     *  (OMEGA*RADIUS*(USINYJ(I-1)+USINYJ(I)) +
     +  .25*(UTANUJ(I-1)+UTANUJ(I)) - (TANxY(J)+TANxY(J+1))*
     *    (UO(I-1,J,L)-UO(I,J,L))*(UO(I-1,J+1,L)-UO(I,J+1,L))/12)
  480 Continue
C****
C**** Vertical advection of momentum
C****
C**** U component
C$OMP ParallelDo   Private (I,J,L, FLUX)
      Do 560 J=J1,JN
      If (J==1)   GoTo 520
      If (J==JM)  GoTo 540
      Do 510 L=1,LMO-1
      FLUX = DT4*(MW(IM,J,L)+MW(1,J,L))*(UO(IM,J,L)+UO(IM,J,L+1))
      dUM(IM,J,L)   = dUM(IM,J,L)   - FLUX
      dUM(IM,J,L+1) = dUM(IM,J,L+1) + FLUX
      Do 510 I=1,IM-1
      FLUX = DT4*(MW(I,J,L)+MW(I+1,J,L))*(UO(I,J,L)+UO(I,J,L+1))
      dUM(I,J,L)   = dUM(I,J,L)   - FLUX
  510 dUM(I,J,L+1) = dUM(I,J,L+1) + FLUX
      GoTo 560
  520 Do 530 L=1,LMO-1
      Do 530 I=1,IM
      FLUX = DT2*MW(1,1,L)*(UO(I,1,L)+UO(I,1,L+1))
      dUM(I,1,L)   = dUM(I,1,L)   - FLUX
  530 dUM(I,1,L+1) = dUM(I,1,L+1) + FLUX
      GoTo 560
  540 Do 550 L=1,LMO-1
      Do 550 I=1,IM
      FLUX = DT2*MW(1,JM,L)*(UO(I,JM,L)+UO(I,JM,L+1))
      dUM(I,JM,L)   = dUM(I,JM,L)   - FLUX
  550 dUM(I,JM,L+1) = dUM(I,JM,L+1) + FLUX
  560 Continue
C**** V component
C$OMP ParallelDo   Private (I,J,L, FLUX)
      Do 660 J=J1,JNP
      If (J==1)  GoTo 620
      If (J==JM-1)  GoTo 640
      Do 610 L=1,LMO-1
      Do 610 I=1,IM
      FLUX = DT4*(MW(I,J,L)+MW(I,J+1,L))*(VO(I,J,L)+VO(I,J,L+1))
      dVM(I,J,L)   = dVM(I,J,L)   - FLUX
  610 dVM(I,J,L+1) = dVM(I,J,L+1) + FLUX
      GoTo 660
  620 Do 630 L=1,LMO-1
      Do 630 I=1,IM
      FLUX = DT4*(MW(1,1,L)+MW(I,2,L))*(VO(I,1,L)+VO(I,1,L+1))
      dVM(I,1,L)   = dVM(I,1,L)   - FLUX
  630 dVM(I,1,L+1) = dVM(I,1,L+1) + FLUX
      GoTo 660
  640 Do 650 L=1,LMO-1
      Do 650 I=1,IM
      FLUX = DT4*(MW(I,JM-1,L)+MW(1,JM,L))*(VO(I,JM-1,L)+VO(I,JM-1,L+1))
      dVM(I,JM-1,L)   = dVM(I,JM-1,L)   - FLUX
  650 dVM(I,JM-1,L+1) = dVM(I,JM-1,L+1) + FLUX
  660 Continue
C****
C**** Add changes to momentum
C****
C$OMP ParallelDo   Private (I,J,L, dUMSP,dVMSP,dUMNP,dVMNP)
      Do 730 L=1,LMO
C**** Calculate dUM and dVM at poles, then update UM and VM
C     If (QSP)  Then
C       dUMSP =   Sum (dUM(:,1,L)*COSU(:)) * 2/IM
C       dVMSP = - Sum (dUM(:,1,L)*SINU(:)) * 2/IM
C       UM2(IVSP,1,L) = UM0(IVSP,1,L) + dVMSP
C       UM2(IM  ,1,L) = UM0(IM  ,1,L) + dUMSP  ;  EndIf
      If (QNP)  Then
        dUMNP = Sum (dUM(:,JM,L)*COSU(:)) * 2/IM
        dVMNP = Sum (dUM(:,JM,L)*SINU(:)) * 2/IM
        UM2(IM  ,JM,L) = UM0(IM  ,JM,L) + dUMNP
        UM2(IVNP,JM,L) = UM0(IVNP,JM,L) + dVMNP  ;  EndIf
C**** Update UM and VM away from poles
      Do 710 J=J1P,JNP
  710 UM2(:,J,L) = UM0(:,J,L) + dUM(:,J,L)
      Do 720 J=J1,JNP
  720 VM2(:,J,L) = VM0(:,J,L) + dVM(:,J,L)
  730 Continue
      Return
      EndSubroutine OADVV

      Subroutine OPGF (UM,VM,DT1)
C****
C**** OPGF adds the pressure gradient force to the momentum
C**** Input: G0M (J), GZM, S0M (kg), SZM, DT (s), MO (kg/m^2)
C**** Output: UM (kg*m/s) = UM - DT*(DH*D(P)+MO*D(PHI))*DYP
C****         VM (kg*m/s) = VM - DT*(DH*D(P)+MO*D(PHI))*DXV
C****
      Use CONSTANT, Only: GRAV
      Use OCEAN, Only: IM,JM,LMO, IVNP,IVSP,
     *                 LMOM=>LMM,LMOU=>LMU,LMOV=>LMV,
     *                 MO, G0M,GZM=>GZMO, S0M,SZM=>SZMO,
     *                 FOCEAN, mZSOLID=>HOCEAN, DXPGF,DYPGF,
     *                 COSI=>COSIC,SINI=>SINIC, OPRESS,OGEOZ,OPBOT
      Use OCEAN_DYN, Only: VBAR,dH, GUP,GDN, SUP,SDN
      Use DOMAIN_DECOMP, Only: GRID, HALO_UPDATE, NORTH
      Implicit None
      Real*8,Parameter :: z12eH=.28867513d0  !  z12eH = 1/SQRT(12)
      Real*8,Intent(Out),
     *  Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) :: UM,VM
      Real*8,Intent(In) :: DT1
      Real*8,External   :: VOLGSP
C****
      Real*8,Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) ::
     *          dUM,dVM, P,ZG
      Real*8    PUP(LMO),PDN(LMO), PE,ZGE,VUP,VDN,dZGdP,
     *          dUMNP(LMO),dVMNP(LMO), DT2
      Integer*4 I,J,L,IMAX, J1,JN,J1P,JNP,JNH
      Logical*4 QSP,QNP
C****
C**** Extract domain decomposition band parameters
C**** ASSUME THAT J1 NEVER EQUALS 2 AND JN NEVER EQUALS JM-1
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      J1P = Max(J1,2)       !    2      5     JM-3   Exclude SP
      JNP = Min(JN,JM-1)    !    4      8     JM-1   Exclude NP
      JNH = Min(JN+1,JM)    !    5      9     JM     Halo maximum
C     QSP = J1==1           !    T      F      F
      QNP = JN==JM          !    F      F      T
C****
      Call HALO_UPDATE (GRID,OPRESS,FROM=NORTH)
C     Call HALO_UPDATE (GRID,  MO, FROM=NORTH)  !  copied in OFLUX
C     Call HALO_UPDATE (GRID, GUP, FROM=NORTH)  !  recalc in OPGF0
C     Call HALO_UPDATE (GRID, GDN, FROM=NORTH)  !  recalc in OPGF0
C     Call HALO_UPDATE (GRID, SUP, FROM=NORTH)  !  recalc in OPGF0
C     Call HALO_UPDATE (GRID, SDN, FROM=NORTH)  !  recalc in OPGF0
      DT2 = DT1/2
C****
C**** Calculate the mass weighted pressure P (Pa),
C**** geopotential ZG (m^2/s^2), and layer thickness dH (m)
C****
C$OMP ParallelDo   Private (I,J,L,IMAX, PUP,PDN,PE, ZGE,VUP,VDN,dZGdP)
      Do 130 J=J1,JNH
      IMAX=IM  ;  If(J==1.or.J==JM) IMAX=1
      Do 130 I=1,IMAX
      If (FOCEAN(I,J) == 0)  GoTo 130
C**** Calculate pressure by integrating from the top down
      PE = OPRESS(I,J)
      Do 110 L=1,LMOM(I,J)
      P(I,J,L) = PE     + MO(I,J,L)*GRAV*.5
      PUP(L) = P(I,J,L) - MO(I,J,L)*GRAV*z12eH
      PDN(L) = P(I,J,L) + MO(I,J,L)*GRAV*z12eH
  110 PE     = PE       + MO(I,J,L)*GRAV
C**** save bottom pressure diagnostic
      OPBOT(I,J)=PE
C**** Calculate geopotential by integrating from the bottom up
      ZGE = - mZSOLID(I,J)*GRAV
      Do 120 L=LMOM(I,J),1,-1
      VUP = VOLGSP (GUP(I,J,L),SUP(I,J,L),PUP(L))
      VDN = VOLGSP (GDN(I,J,L),SDN(I,J,L),PDN(L))
      dZGdP = VUP*(.5-z12eH) + VDN*(.5+z12eH)
      VBAR(I,J,L) = (VUP + VDN)*.5
      ZG(I,J,L) = MO(I,J,L)*GRAV*.5*dZGdP + ZGE
      dH(I,J,L) = MO(I,J,L)*VBAR(I,J,L)
  120 ZGE = ZGE + dH(I,J,L)*GRAV
      OGEOZ(I,J) = ZGE
  130 Continue
C**** Copy VBAR and DH to all longitudes at north pole
      If (QNP) Then
        Do 140 L=1,LMOM(1,JM)
        VBAR(2:IM,JM,L) = VBAR(1,JM,L)
  140     DH(2:IM,JM,L) =   DH(1,JM,L)
      EndIf
C****
C**** Calculate smoothed East-West Pressure Gradient Force
C****
C$OMP ParallelDo   Private (I,J,L)
      Do 220 J=J1P,JNP
      Do 210 I=1,IM-1
      Do 210 L=1,LMOU(I,J)
  210 dUM(I,J,L) = DT2*DYPGF(J)*
     *             ((dH(I,J,L)+dH(I+1,J,L))*( P(I,J,L)- P(I+1,J,L)) +
     +              (MO(I,J,L)+MO(I+1,J,L))*(ZG(I,J,L)-ZG(I+1,J,L)))
      Do 220 L=1,LMOU(IM,J)
  220 dUM(IM,J,L) = DT2*DYPGF(J)*
     *              ((dH(IM,J,L)+dH(1,J,L))*( P(IM,J,L)- P(1,J,L)) +
     +               (MO(IM,J,L)+MO(1,J,L))*(ZG(IM,J,L)-ZG(1,J,L)))
      Call OPFIL (dUM)
C****
C**** Calculate North-South Pressure Gradient Force
C****
C$OMP ParallelDo   Private (I,J,L)
      Do 340 J=J1,JNP
      If (J==JM-1)  GoTo 320
      Do 310 I=1,IM
      Do 310 L=1,LMOV(I,J)
  310 dVM(I,J,L) = DT2*DXPGF(J)*
     *  ((dH(I,J,L)+dH(I,J+1,L))*( P(I,J,L)- P(I,J+1,L)) +
     +   (MO(I,J,L)+MO(I,J+1,L))*(ZG(I,J,L)-ZG(I,J+1,L)))
      GoTo 340
  320 Do 330 I=1,IM
      Do 330 L=1,LMOV(I,JM-1)
  330 dVM(I,JM-1,L) = DT2*DXPGF(JM-1)*
     *  ((dH(I,JM-1,L)+dH(1,JM,L))*( P(I,JM-1,L)- P(1,JM,L)) +
     +   (MO(I,JM-1,L)+MO(1,JM,L))*(ZG(I,JM-1,L)-ZG(1,JM,L)))
  340 Continue
C****
C**** Pressure Gradient Force at north pole
C****
      If (QNP)  Then
        dUMNP(:) = 0  ;  dVMNP(:) = 0
        Do 510 I=1,IM
        Do 510 L=1,LMOV(I,JM-1)
        dUMNP(L) = dUMNP(L) - SINI(I)*dVM(I,JM-1,L)
  510   dVMNP(L) = dVMNP(L) + COSI(I)*dVM(I,JM-1,L)
        Do 520 L=1,LMOM(1,JM)
        dUMNP(L) = dUMNP(L)*4*DXPGF(JM) / (IM*DXPGF(JM-1))
  520   dVMNP(L) = dVMNP(L)*4*DXPGF(JM) / (IM*DXPGF(JM-1))  ;  EndIf
C****
C**** Add pressure gradient force to momentum
C****
C**** Update UM away from poles
C$OMP ParallelDo   Private (I,J,L)
      Do 630 J=J1,JNP
      Do 620 I=1,IM
      Do 610 L=1,LMOU(I,J)
  610 UM(I,J,L) = UM(I,J,L) + dUM(I,J,L)
C**** Update VM away from poles
      Do 620 L=1,LMOV(I,J)
  620 VM(I,J,L) = VM(I,J,L) + dVM(I,J,L)
  630 Continue
C**** UM and VM at poles
      If (QNP)  Then
        Do 650 L=1,LMOM(1,JM)
        UM(IM  ,JM,L) = UM(IM  ,JM,L) + dUMNP(L)
  650   UM(IVNP,JM,L) = UM(IVNP,JM,L) + dVMNP(L)  ;  EndIf
      Return
      End Subroutine OPGF

      Subroutine OPGF0
C****
C**** OPGF0 calculates GUP, GDN, SUP and SDN inside each ocean cell,
C**** which will be kept constant during the dynamics of momentum
C****
      Use OCEAN, Only: IM,JM,LMO, LMOM=>LMM,
     *                 G0M,GZM=>GZMO, S0M,SZM=>SZMO
      Use OCEAN_DYN, Only: MMI, GUP,GDN, SUP,SDN
      Use DOMAIN_DECOMP, Only: GRID, HALO_UPDATE, NORTH
      Implicit None
      Real*8,Parameter :: z12eH=.28867513d0  !  z12eH = 1/SQRT(12)
      Integer*4 I,J,L,IMAX, J1,JN,JNH
C****
C**** Extract domain decomposition band parameters
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      JNH = Min(JN+1,JM)    !    5      9     JM     Halo maximum
C****
      Call HALO_UPDATE (GRID, MMI, FROM=NORTH)
      Call HALO_UPDATE (GRID, G0M, FROM=NORTH)
      Call HALO_UPDATE (GRID, GZM, FROM=NORTH)
      Call HALO_UPDATE (GRID, S0M, FROM=NORTH)
      Call HALO_UPDATE (GRID, SZM, FROM=NORTH)
C****
C$OMP ParallelDo   Private (I,J,L,IMAX)
      Do 10 J=J1,JNH
      IMAX=IM  ;  If(J==1.or.J==JM) IMAX=1
      Do 10 I=1,IMAX
      Do 10 L=1,LMOM(I,J)
      GUP(I,J,L) = (G0M(I,J,L) - 2*z12eH*GZM(I,J,L)) / MMI(I,J,L)
      GDN(I,J,L) = (G0M(I,J,L) + 2*z12eH*GZM(I,J,L)) / MMI(I,J,L)
      SUP(I,J,L) = (S0M(I,J,L) - 2*z12eH*SZM(I,J,L)) / MMI(I,J,L)
   10 SDN(I,J,L) = (S0M(I,J,L) + 2*z12eH*SZM(I,J,L)) / MMI(I,J,L)
      Return
      End Subroutine OPGF0

      SUBROUTINE OADVT (RM,RX,RY,RZ,DT,QLIMIT, OIJL)
!@sum  OADVT advects tracers using the linear upstream scheme.
!@auth Gary Russell
!@ver  1.0
C****
C**** Input:  MB (kg) = mass before advection
C****          DT (s) = time step
C****       MU (kg/s) = west to east mass flux
C****       MV (kg/s) = south to north mass flux
C****       MW (kg/s) = downward vertical mass flux
C****          QLIMIT = whether slope limitations should be used
C**** Output: RM (kg) = tracer mass
C****   RX,RY,RZ (kg) = first moments of tracer mass
C****       OIJL (kg) = diagnostic accumulation of tracer mass flux
C****
      USE OCEAN, only : im,jm,lmo
      USE OCEAN_DYN, only : mb=>mmi,smu,smv,smw
      use domain_decomp, only : grid, get

      IMPLICIT NONE
      REAL*8, INTENT(INOUT),     DIMENSION
     *     (IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) :: RM,RX,RY,RZ
      REAL*8, INTENT(INOUT),
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO,3) :: OIJL
      REAL*8,
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) :: MA
      INTEGER I,J,L, J_0H
      LOGICAL, INTENT(IN) :: QLIMIT
      REAL*8, INTENT(IN) :: DT

      logical :: HAVE_NORTH_POLE  ! ,HAVE_SOUTH_POLE
         call get (grid, J_STRT_HALO=J_0H)
         call get (grid, HAVE_NORTH_POLE=HAVE_NORTH_POLE)
C        call get (grid, HAVE_SOUTH_POLE=HAVE_SOUTH_POLE)

C****
C**** Load mass after advection from mass before advection
C****
      MA = MB
C****
C**** Advect the tracer using the Slopes Scheme
C****
      CALL OADVTX (RM,RX,RY,RZ,MA,SMU,5d-1*DT,QLIMIT,OIJL(1,J_0H,1,1))
      CALL OADVTY (RM,RX,RY,RZ,MA,SMV,     DT,QLIMIT,OIJL(1,J_0H,1,2))
      CALL OADVTZ (RM,RX,RY,RZ,MA,SMW,     DT,QLIMIT,OIJL(1,J_0H,1,3))
      CALL OADVTX (RM,RX,RY,RZ,MA,SMU,5d-1*DT,QLIMIT,OIJL(1,J_0H,1,1))
C**** Fill in values at the poles
      if (HAVE_NORTH_POLE) then
      DO 20 L=1,LMO
      DO 20 I=1,IM
      RM(I,JM,L) = RM(1,JM,L)
   20 RZ(I,JM,L) = RZ(1,JM,L)
      end if
C     if (HAVE_SOUTH_POLE) then
C     DO 30 L=1,LMO
C     DO 30 I=1,IM
C     RM(I, 1,L) = RM(IM,1,L)
C  30 RZ(I, 1,L) = RZ(IM,1,L)
C     end if

      RETURN
      END SUBROUTINE OADVT

      SUBROUTINE OADVTX (RM,RX,RY,RZ,MO,MU,DT,QLIMIT,OIJL)
!@sum  OADVTX advects tracer in x direction using linear upstream scheme
!@auth Gary Russell
!@ver  1.0
C**** If QLIMIT is true, the gradients are
C**** limited to prevent the mean tracer from becoming negative.
C****
C**** Input: DT (s) = time step
C****     MU (kg/s) = west to east mass flux
C****        QLIMIT = whether slope limitations should be used
C**** Input and Output: RM (kg) = tracer mass
C****             RX,RY,RZ (kg) = first moments of tracer mass
C****                    M (kg) = ocean mass
C****
      USE OCEAN, only : im,jm,lmo,lmm,focean
      use domain_decomp, only : grid, get

      IMPLICIT NONE
      REAL*8, INTENT(INOUT),
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) ::
     *  RM,RX,RY,RZ, OIJL, MO
      REAL*8, INTENT(IN),
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) :: MU
      LOGICAL*4, INTENT(IN) :: QLIMIT
      REAL*8, INTENT(IN) :: DT
      REAL*8, DIMENSION(IM) :: AM,A,FM,FX,FY,FZ
      REAL*8 RXY
      INTEGER I,J,L,IM1,IP1,ICKERR

      INTEGER :: J_0S,J_1S

      CALL GET(grid, J_STRT_SKP=J_0S,    J_STOP_SKP=J_1S)

C**** Loop over layers and latitudes
      ICKERR=0
!$OMP PARALLEL DO  PRIVATE(I,J,L,IP1,IM1,A,AM,FM,FX,FY,FZ,RXY)
!$OMP&             REDUCTION(+:ICKERR)
      DO 320 L=1,LMO
      DO 320 J=J_0S,J_1S
C****
C**** Calculate FM (kg), FX (kg**2), FY (kg) and FZ (kg)
C****
      I=IM
      DO 140 IP1=1,IM
      AM(I) = DT*MU(I,J,L)
      IF(AM(I)) 110,120,130
C**** Ocean mass flux is negative
  110 A(I)  = AM(I)/MO(IP1,J,L)
      IF(A(I).LT.-1d0)  WRITE (6,*) 'A<-1:',I,J,L,A(I),MO(IP1,J,L)
      FM(I) = A(I)*(RM(IP1,J,L)-(1d0+A(I))*RX(IP1,J,L))
      FX(I) = AM(I)*(A(I)*A(I)*RX(IP1,J,L)-3d0*FM(I))
      FY(I) = A(I)*RY(IP1,J,L)
      FZ(I) = A(I)*RZ(IP1,J,L)
      GO TO 140
C**** Ocean mass flux is zero
  120 A(I)  = 0.
      FM(I) = 0.
      FX(I) = 0.
      FY(I) = 0.
      FZ(I) = 0.
      GO TO 140
C**** Ocean mass flux is positive
  130 A(I)  = AM(I)/MO(I,J,L)
      IF(A(I).GT.1d0)  WRITE (6,*) 'A>1:',I,J,L,A(I),MO(I,J,L)
      FM(I) = A(I)*(RM(I,J,L)+(1d0-A(I))*RX(I,J,L))
      FX(I) = AM(I)*(A(I)*A(I)*RX(I,J,L)-3d0*FM(I))
      FY(I) = A(I)*RY(I,J,L)
      FZ(I) = A(I)*RZ(I,J,L)
  140 I=IP1
C****
C**** Modify the tracer moments so that the tracer mass in each
C**** division is non-negative
C****
      IF(.NOT.QLIMIT)  GO TO 300
      IM1=IM
      DO 290 I=1,IM
      IF(A(IM1).GE.0.)  GO TO 240
C**** Water is leaving through the left edge: 2 or 3 divisions
      IF(FM(IM1).LE.0.)  GO TO 210
C**** Left most division is negative, RML = -FM(I-1) < 0: Case 2 or 4
      RX(I,J,L) = RM(I,J,L)/(1d0+A(IM1))
      FM(IM1) = 0.
      FX(IM1) = AM(IM1)*A(IM1)*A(IM1)*RX(I,J,L)
      IF(A(I).LE.0.)  GO TO 290
      FM(I) = A(I)*(RM(I,J,L)+(1d0-A(I))*RX(I,J,L))
      FX(I) = AM(I)*(A(I)*A(I)*RX(I,J,L)-3d0*FM(I))
      GO TO 290
C**** Left most division is non-negative, RML = -FM(I-() > 0:
C**** Case 1, 3 or 5
  210 IF(A(I).LE.0.)  GO TO 230
C**** Water is leaving through the right edge: 3 divisions
      IF(FM(I).GE.0.)  GO TO 290
C**** Right most division is negative, RMR = FM(I) < 0: Case 3 or 5
  220 RX(I,J,L) = -RM(I,J,L)/(1d0-A(I))
      FM(I) = 0.
      FX(I) = AM(I)*A(I)*A(I)*RX(I,J,L)
      FM(IM1) = A(IM1)*(RM(I,J,L)-(1d0+A(IM1))*RX(I,J,L))
      FX(IM1) = AM(IM1)*(A(IM1)*A(IM1)*RX(I,J,L)-3d0*FM(IM1))
      GO TO 290
C**** No water is leaving through the right edge: 2 divisions
  230 IF(RM(I,J,L)+FM(IM1).GE.0.)  GO TO 290
C**** Right most division is negative, RMR = RM(I,J)+FM(I-1) < 0: Case 3
      RX(I,J,L) = RM(I,J,L)/A(IM1)
      FM(IM1) = -RM(I,J,L)
      FX(IM1) = AM(IM1)*(A(IM1)+3d0)*RM(I,J,L)
      GO TO 290
C**** No water is leaving through the left edge: 1 or 2 divisions
  240 IF(A(I).LE.0.)  GO TO 290
C**** Water is leaving through the right edge: 2 divisions
      IF(FM(I).GE.0.)  GO TO 250
C**** Right most division is negative, RMR = FM(I) < 0: Case 3
      RX(I,J,L) = -RM(I,J,L)/(1d0-A(I))
      FM(I) = 0.
      FX(I) = AM(I)*A(I)*A(I)*RX(I,J,L)
      GO TO 290
C**** Right most division is non-negative, RMR = FM(I) > 0: Case 1 or 2
  250 IF(RM(I,J,L)-FM(I).GE.0.)  GO TO 290
C**** Left most division is negative, RML = RM(I,J)-FM(I) < 0: Case 2
      RX(I,J,L) = RM(I,J,L)/A(I)
      FM(I) = RM(I,J,L)
      FX(I) = AM(I)*(A(I)-3d0)*RM(I,J,L)
C****
  290 IM1=I
C****
C**** Calculate new tracer mass and first moments of tracer mass
C****
  300 IM1=IM
      DO 310 I=1,IM
      IF(L.GT.LMM(I,J))  GO TO 310
      RM(I,J,L) = RM(I,J,L) +  (FM(IM1)-FM(I))
      RX(I,J,L) = (RX(I,J,L)*MO(I,J,L) + (FX(IM1)-FX(I))
     *  + 3d0*((AM(IM1)+AM(I))*RM(I,J,L)-MO(I,J,L)*(FM(IM1)+FM(I))))
     *  / (MO(I,J,L)+AM(IM1)-AM(I))
      RY(I,J,L) = RY(I,J,L) + (FY(IM1)-FY(I))
      RZ(I,J,L) = RZ(I,J,L) + (FZ(IM1)-FZ(I))
C****
      if ( QLIMIT ) then ! limit tracer gradients
        RXY = abs(RX(I,J,L)) + abs(RY(I,J,L))
        if ( RXY > RM(I,J,L) ) then
          RX(I,J,L) = RX(I,J,L)*( RM(I,J,L)/(RXY + tiny(RXY)) )
          RY(I,J,L) = RY(I,J,L)*( RM(I,J,L)/(RXY + tiny(RXY)) )
        end if
        if ( abs(RZ(I,J,L)) > RM(I,J,L) )
     *       RZ(I,J,L) = sign(RM(I,J,L), RZ(I,J,L)+0d0)
      end if
C****
      MO(I,J,L) = MO(I,J,L) +  AM(IM1)-AM(I)
         IF(MO(I,J,L).LE.0.)                ICKERR=ICKERR+1
         IF(QLIMIT .AND. RM(I,J,L).LT.0.)   ICKERR=ICKERR+1
         OIJL(I,J,L) = OIJL(I,J,L) + FM(I)
  310 IM1=I
  320 CONTINUE
!$OMP END PARALLEL DO

C**** IF NO ERROR HAS OCCURRED - RETURN, ELSE STOP
      IF(ICKERR == 0)  RETURN
      DO 440 L=1,LMO
      DO 440 J=J_0S,J_1S
      DO 440 I=1,IM
         IF(FOCEAN(I,J).gt.0 .and. MO(I,J,L).LE.0.)  GO TO 800
         IF(QLIMIT .AND. RM(I,J,L).LT.0.) GO TO 810
  440 CONTINUE
      WRITE(6,*) 'ERROR CHECK INCONSISTENCY: OADVTX ',ICKERR
      call stop_model("OADVTX",255)

  800 WRITE (6,*) 'MO<0 in OADVTX:',I,J,L,MO(I,J,L)
  810 WRITE (6,*) 'RM in OADVTX:',I,J,L,RM(I,J,L)
c      WRITE (6,*) 'A=',(I,A(I),I=1,IM)
      call stop_model("OADVTX",255)
      END SUBROUTINE OADVTX

      subroutine oadvty (rm,rx,ry,rz,mo,mv,dt,qlimit,oijl)
!@sum  OADVTY advection driver for y-direction
!@auth Gary Russell, modified by T.Clune, R. Ruedy
c****
c**** oadvty advects tracers in the south to north direction using the
c**** linear upstream scheme.  if qlimit is true, the moments are
c**** limited to prevent the mean tracer from becoming negative.
c****
c**** input:
c****     mv (kg/s) = north-south mo flux, positive northward
c****        qlimit = whether slope limitations should be used
c****
c**** input/output:
c****     rm     (kg) = tracer mass
c****   rx,ry,rz (kg) = 1st moments of tracer mass
c****     mo     (kg) = ocean mass
c****
      use DOMAIN_DECOMP, only : grid, get, halo_update
      use DOMAIN_DECOMP, only : halo_update_column
      use DOMAIN_DECOMP, only : NORTH, SOUTH, AM_I_ROOT
      use OCEAN, only : im,jm,lmo,lmm,focean
      implicit none
      REAL*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &                  rm, rx,ry,rz, mo,mv, oijl
      real*8, intent(in) :: dt
      logical, intent(in) ::  qlimit

      REAL*8, dimension(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) ::
     &        BM, fm,fx,fy,fz
      integer :: i,j,L,ierr,ICKERR, err_loc(3)
      REAL*8, DIMENSION(LMO) ::
     &     m_np,rm_np,rzm_np ! ,m_sp,rm_sp,rzm_sp

c****Get relevant local distributed parameters
      INTEGER J_0,J_1,J_0H,J_1H,J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      CALL GET(grid, J_STRT = J_0,
     &               J_STOP = J_1, J_STOP_SKP=J_1S,
     &               J_STRT_HALO = J_0H,
     &               J_STOP_HALO = J_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

c**** loop over layers
      ICKERR=0
!$OMP  PARALLEL DO PRIVATE(I,J,L) DEFAULT(SHARED)
!!! !$OMP* SHARED(JM,im)
      do L=1,lmo

c****   fill in and save polar values
c       if (HAVE_SOUTH_POLE) then
c         m_sp(L) = mo(1,1,L)
c         mo(2:im,1,L) = mo(1,1,L)
c         rm_sp(L) = rm(1,1,L)
c         rm(2:im,1,L) = rm(1,1,L)
c         rzm_sp(L)  = rz(1,1,L)
c         rz(2:im,1,L)  = rz(1,1,L)
c       end if                       !SOUTH POLE

        if (HAVE_NORTH_POLE) then
          m_np(L) = mo(1,jm,L)
          mo(2:im,jm,L) = mo(1,jm,L)
          rm_np(L) = rm(1,jm,l)
          rm(2:im,jm,L) = rm(1,jm,L)
          rzm_np(L)  = rz(1,jm,l)
          rz(2:im,jm,L)  = rz(1,jm,L)
        end if                       !NORTH POLE
c****
c****   convert flux to water mass moving north (if >0)
        do j=j_0,j_1
        do i=1,im
          if(focean(i,j).gt.0. .and. L.le.lmm(i,j)) then
            bm(i,j,L)=dt*mv(i,j,L)
          else
            bm(i,j,L)=0.
            mo(i,j,L)=0.
          end if
        end do
        end do

c****   POLES: set horiz. moments to zero
        IF (HAVE_SOUTH_POLE) THEN
          rx(:,1,L) = 0. ; ry(:,1,L) = 0.
        end if
        IF (HAVE_NORTH_POLE) THEN
          bm(:,jm,L) = 0.
          rx(:,jm,L) = 0. ; ry(:,jm,L) = 0.
        END IF
      end do   ! loop over layers
!$OMP END PARALLEL DO

c****
c**** call 1-d advection routine
c****
        call advec_lin_1D_custom(
     &     rm(1,j_0h,1), rx(1,j_0h,1),ry(1,j_0h,1),rz(1,j_0h,1),
     &     fm(1,j_0h,1), fx(1,j_0h,1),fy(1,j_0h,1),fz(1,j_0h,1),
     &     mo(1,j_0h,1), bm(1,j_0h,1),
     &     qlimit,ierr,err_loc)

        if (ierr.gt.0) then
          write(6,*) "Error in oadvty: i,j,l=",err_loc
          if (ierr == 2) then
ccc         write(0,*) "Error in qlimit: abs(b) > 1"
ccc         call stop_model('Error in qlimit: abs(b) > 1',11)
            ICKERR=ICKERR+1
          endif
        end if
! horizontal moments are zero at pole
c       IF (HAVE_SOUTH_POLE) then
c          rx(:,1, :) = 0  ; ry(:,1, :) = 0
c       end if
        IF (HAVE_NORTH_POLE) then
           rx(:,jm,:) = 0. ; ry(:,jm,:) = 0.
        end if

!$OMP  PARALLEL DO PRIVATE(I,L) DEFAULT(SHARED)
!$OMP&             REDUCTION(+:ICKERR)
!!! !$OMP* SHARED(JM,IM)
      do l=1,lmo

c****   average and update polar boxes
c       if (HAVE_SOUTH_POLE) then
c         mo(:,1 ,l) = (m_sp(l) + sum(mo(:,1 ,l)-m_sp(l)))/im
c         rm(:,1 ,l) = (rm_sp(l) + sum(rm(:,1 ,l)-rm_sp(l)))/im
c         rz(:,1 ,l) = (rzm_sp(l) + sum(rz(:,1 ,l)-rzm_sp(l) ))/im
c       end if   !SOUTH POLE

        if (HAVE_NORTH_POLE .and. L.le.lmm(1,jm)) then
          mo(1,jm,l) = m_np(l)   + sum(bm(:,jm-1,l))/im
          rm(1,jm,l) = rm_np(l)  + sum(fm(:,jm-1,l))/im
          rz(1,jm,l) = rzm_np(l) + sum(fz(:,jm-1,l))/im
          if(mo(1,jm,l)<0. .or. (qlimit.and.rm(1,jm,l)<0.)) then
            ICKERR=ICKERR+1
            write(0,*) 'oadvty: mo or salt<0 at North Pole layer',
     *      L,mo(1,jm,L),rm(1,jm,L)
          endif
          mo(2:im,jm,l)=mo(1,jm,l)
          rm(2:im,jm,l)=rm(1,jm,l)
          rz(2:im,jm,l)=rz(1,jm,l)
        end if  !NORTH POLE

      enddo ! end loop over levels
!$OMP  END PARALLEL DO
c
c**** sum into oijl
c
!$OMP  PARALLEL DO PRIVATE(J,L)
      do l=1,lmo ; do j=J_0,J_1S
          oijl(:,j,l)  = oijl(:,j,l) + fm(:,j,l)
      enddo ; enddo
!$OMP  END PARALLEL DO
C
      IF(ICKERR.GT.0)  CALL stop_model('Stopped in oadvty',11)
C
      return
c****
      end subroutine oadvty

      subroutine advec_lin_1D_custom(s,sx,sy,sz, f,fx,fy,fz, mass,dm,
     *     qlimit,ierr, err_loc)
!@sum  advec_lin_1d_custom is a parallel variant of adv1d but for a
!@sum  linear upstream scheme.
!@auth T. Clune, R. Ruedy
c--------------------------------------------------------------
c adv1d advects tracers in j-direction using the lin ups scheme
c--------------------------------------------------------------
      use ocean, only: im,jm,lmo,lmm,focean
      USE DOMAIN_DECOMP, only: grid, GET
      USE DOMAIN_DECOMP, only: NORTH, SOUTH
      USE DOMAIN_DECOMP, only: HALO_UPDATE, HALO_UPDATE_COLUMN
      USE DOMAIN_DECOMP, only: CHECKSUM
      implicit none
      !
!@var s mean tracer amount (kg or J)
!@var sx,sy,sz lin tracer moments (kg or J)
!@var f tracer flux (diagnostic output) (kg or J)
!@var fx,fy,fz tracer moment flux (diagnostic output) (kg or J)
!@var mass mass field (kg)
!@var dm mass flux (kg)
!@var qlimit true if negative tracer is to be avoided
!@var ierr, nerr error codes
      logical, intent(in) :: qlimit
      REAL*8, dimension(im, grid%j_strt_halo:grid%j_stop_halo,lmo)
     *  :: s,sx,sy,sz, f,fx,fy,fz, mass,dm
      integer :: n,np1,nm1,nn,ns
      integer,intent(out) :: ierr,err_loc(3) ! 3 dimensions
      REAL*8 :: fracm,frac1,mnew
      INTEGER :: i
      INTEGER :: l
      INTEGER :: J_0, J_1, J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE

      ierr=0

      CALL GET(grid, J_STRT = J_0, J_STOP=J_1,
     & J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     & HAVE_SOUTH_POLE=HAVE_SOUTH_POLE)

      CALL HALO_UPDATE(grid, mass, FROM=NORTH)
      CALL HALO_UPDATE(grid, s, FROM=NORTH)
      CALL HALO_UPDATE(grid, sx, FROM=NORTH)
      CALL HALO_UPDATE(grid, sy, FROM=NORTH)
      CALL HALO_UPDATE(grid, sz, FROM=NORTH)

      DO l=1,lmo
        Do i=1, im
          Call calc_tracer_mass_flux() ! f from s,sy,mass,dm
        End Do                         ! temp. store dm/mass in fy
      end do

      CALL HALO_UPDATE(grid, fy, FROM=NORTH+SOUTH) ! still holding dm/m
      ! Limit fluxes to maintain positive mean values?
      If (qlimit) Then

        Call HALO_UPDATE(grid, f, FROM=NORTH+SOUTH)
        DO l=1,lmo
          Do i=1, im
            Call apply_limiter() ! adjusting f,sy
          End Do
        end do
        Call HALO_UPDATE(grid, f, FROM=NORTH+SOUTH)
      End If
c--------------------------------------------------------------------
         ! calculate tracer fluxes of slopes fx,fy,fz
c--------------------------------------------------------------------
      DO l=1,lmo
        Do i=1, im
          Call tracer_slopes()
        end do
      end do
c-------------------------------------------------------------------
c update tracer mass, moments of tracer mass, air mass distribution
c-------------------------------------------------------------------
      CALL HALO_UPDATE(grid, f,  FROM=SOUTH)
      CALL HALO_UPDATE(grid, dm, FROM=SOUTH)
      CALL HALO_UPDATE(grid, fx, FROM=SOUTH)
      CALL HALO_UPDATE(grid, fy, FROM=SOUTH)
      CALL HALO_UPDATE(grid, fz, FROM=SOUTH)

      DO l=1,lmo
        DO i = 1, im
          Call update_tracer_mass() ! s also: sx,sy,sz, mass
        end do
      enddo

      return

      Contains

      Integer Function NeighborByFlux(n, dm)
        Integer, Intent(In) :: n
        Real*8, Intent(In) :: dm

        Integer :: nn

        If (dm < 0) Then ! air mass flux is negative
          nn=n+1
        else ! air mass flux is positive
          nn=n
        endif
        NeighborByFlux = nn
      End Function NeighborByFlux

      Function FluxFraction(dm) Result(frac)
        Real*8, Intent(In) :: dm
        Real*8 :: frac

        If (dm < 0 ) Then
          frac = +1.
        Else ! Flux non negative
          frac = -1.
        End If
      End Function FluxFraction

      Function MassFraction(dm, mass) Result (fracm)
        Real*8, Intent(In) :: dm
        Real*8, Intent(In) :: mass
        Real*8 :: fracm

        If (mass > 0.0d0) Then
          fracm = dm / mass
        Else
          fracm = 0.d0
        End If
      End Function MassFraction

      Subroutine calc_tracer_mass_flux()

      Do n = J_0, J_1

        nn = NeighborByFlux(n, dm(i,n,l))
        fracm = MassFraction(dm(i,n,l), mass(i,nn,l))

        frac1 = fracm + FluxFraction(dm(i,n,l))

        f(i,n,l)=fracm*(s(i,nn,l)-frac1*sy(i,nn,l))
      ! temporary storage of fracm in fy, to be used below
        fy(i,n,l)=fracm
      !
      enddo
      End Subroutine calc_tracer_mass_flux

      Subroutine tracer_slopes()

      Do n = J_0, J_1

        nn = NeighborByFlux(n, dm(i,n,l))
      ! retrieving fracm, which was stored in fy
        fracm=fy(i,n,l)
      !
        fy(i,n,l)=
     &     dm(i,n,l)*(fracm*fracm*sy(i,nn,l)-3.*f(i,n,l))
      ! cross moments
        fx(i,n,l)=fracm*sx(i,nn,l)
        fz(i,n,l)=fracm*sz(i,nn,l)
      enddo
      End Subroutine tracer_slopes

      Subroutine apply_limiter
      REAL*8 :: an, anm1, fn, fnm1, sn, syn

c     If (HAVE_SOUTH_POLE) Then
c        n = J_0
c        an = fy(i,n,l) ! reading fracm which was stored in fx
c        anm1 = 0
c        fn = f(i,n,l)
c        fnm1 = 0
c        sn = s(i,n,l)
c        syn = sy(i,n,l)
c        call limitq_lin(anm1,an,fnm1,fn,sn,syn)
c        f(i,n,l) = fn
c        sy(i,n,l) = syn
c     End If

      DO n = J_0S, J_1S+1

         an = fy(i,n,l) ! reading fracm which was stored in fy
         anm1 = fy(i,n-1,l)
         fn = f(i,n,l)
         fnm1 = f(i,n-1,l)
         sn = s(i,n,l)
         syn = sy(i,n,l)
         call limitq_lin(anm1,an,fnm1,fn,sn,syn)
         f(i,n,l) = fn
         f(i,n-1,l) = fnm1
         sy(i,n,l) = syn

      enddo
      End Subroutine apply_limiter

      Subroutine update_tracer_mass() ! non-polar only
      REAL*8 :: sxy

      Do n=J_0S,J_1S

         if(focean(i,n).le.0. .or. l.gt.lmm(i,n)) cycle
         mnew=mass(i,n,l) + dm(i,n-1,l)-dm(i,n,l)

         s(i,n,l)=s(i,n,l) + (f(i,n-1,l)-f(i,n,l))

         sy(i,n,l)=( sy(i,n,l)*mass(i,n,l) + (fy(i,n-1,l)-fy(i,n,l))
     &              + 3.*( (dm(i,n-1,l)+dm(i,n,l))*s(i,n,l)
     &                     -mass(i,n,l)*(f(i,n-1,l)+f(i,n,l)) ) )/mnew
      ! cross moments
         sx(i,n,l) = sx(i,n,l) + (fx(i,n-1,l)-fx(i,n,l))
         sz(i,n,l) = sz(i,n,l) + (fz(i,n-1,l)-fz(i,n,l))
      !
         if (qlimit) then ! limit tracer gradients
           sxy = abs(sx(i,n,l)) + abs(sy(i,n,l))
           if ( sxy > s(i,n,l) ) then
             sx(i,n,l) = sx(i,n,l)*( s(i,n,l)/(sxy + tiny(sxy)) )
             sy(i,n,l) = sy(i,n,l)*( s(i,n,l)/(sxy + tiny(sxy)) )
           end if
           if ( abs(sz(i,n,l)) > s(i,n,l) )
     *       sz(i,n,l) = sign(s(i,n,l),sz(i,n,l)+0.d0)
         end if
c------------------------------------------------------------------
         mass(i,n,l) = mnew
         if(mass(i,n,l).le.0. .or. (qlimit.and.s(i,n,l).lt.0.)) then
            ierr=2
            err_loc=(/ i, n, l /)
            write(0,*) 'oadvty: mo or salt<0 at',err_loc,mnew,s(i,n,l)
            return
         endif
c-----------------------------------------------------------------

      enddo
      End Subroutine Update_Tracer_Mass

      end subroutine advec_lin_1D_custom

      subroutine limitq_lin(anm1,an,fnm1,fn,sn,sx)
!@sum  limitq adjusts moments to maintain non-neg. tracer means/fluxes
!@auth G. Russell, modified by Maxwell Kelley
        implicit none
        REAL*8 :: anm1,an,fnm1,fn,sn,sx
c local variables
        REAL*8 :: sl,sc,sr, frl,frl1, frr,frr1, gamma,g13ab,
     &       fr,fr1, fsign,su,sd
c****
c**** modify the tracer moments so that the tracer mass in each
c**** division is non-negative
c****
c**** no water leaving the box
        if(anm1.ge.0. .and. an.le.0.) return
c**** water is leaving through both the left and right edges
        if(anm1.lt.0. .and. an.gt.0.) then
           sl = -fnm1
           sr = +fn
c**** all divisions are non-negative
           if(sl.ge.0. .and. sr.ge.0.) return
c**** at least one division is negative
           frl = anm1
           frl1 = frl+1.
           frr = an
           frr1 = frr-1.
           if(sl.lt.0.) then           ! leftmost division
              sx = sn/frl1
              sr = frr*(sn-frr1*sx)
              sl = 0.
           else                        ! rightmost division
              sx = sn/frr1
              sl = -frl*(sn-frl1*sx)
              sr = 0.
           endif
           fnm1 = -sl
           fn   = +sr
        else
c**** water is leaving only through one edge
           if(an.gt.0.)  then ! right edge
              fr=an
              sd=fn
              fsign=-1.
           else                  ! left edge
              fr=anm1
              sd=-fnm1
              fsign=1.
           endif
           su = sn-sd
           if(sd.ge.0. .and. su.ge.0.) return
           fr1=fr+fsign
           if(sd.lt.0.)  then
c**** downstream division is negative
              sx = sn/fr1
              su = sn
           else
c**** upstream division is negative
              sx = sn/fr
              su = 0.
           endif
           sd = sn - su
           if(an.gt.0.) then
              fn=sd
           else
              fnm1=-sd
           endif
        endif
        return
      end subroutine limitq_lin

      SUBROUTINE OADVTZ (RM,RX,RY,RZ,MO,MW,DT,QLIMIT,OIJL)
C****
C**** OADVTZ advects tracers in the vertical direction using the
C**** linear upstream scheme.  If QLIMIT is true, the gradients are
C**** limited to prevent the mean tracer from becoming negative.
C****
C**** Input: DT (s) = time step
C****     MW (kg/s) = downward vertical mass flux
C****        QLIMIT = whether slope limitations should be used
C**** Input and Output: RM (kg) = tracer mass
C****             RX,RY,RZ (kg) = first moments of tracer mass
C****                    M (kg) = ocean mass
C****
      USE OCEAN, only : im,jm,lmo,lmm,focean
      use domain_decomp, only : grid, get

      IMPLICIT NONE
      REAL*8, INTENT(INOUT),
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) ::
     *  RM,RX,RY,RZ, OIJL, MO
      REAL*8, INTENT(IN),
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO-1) :: MW
      LOGICAL*4, INTENT(IN) :: QLIMIT
      REAL*8, INTENT(IN) :: DT
      REAL*8, DIMENSION(0:LMO) :: CM,C,FM,FX,FY,FZ
      INTEGER I,J,L,LMIJ,ICKERR,IMIN,IMAX
      REAL*8 SBMN,SFMN,SFZN,RXY

      INTEGER :: J_0,J_1

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

C****
C**** Loop over latitudes and longitudes
      ICKERR=0
!$OMP PARALLEL DO  PRIVATE(I,IMIN,IMAX,J,L,LMIJ, C,CM, FM,FX,FY,FZ,RXY)
!$OMP&             REDUCTION(+:ICKERR)
      DO J=J_0,J_1
        IMIN=1
        IMAX=IM
        IF (J == 1) IMIN=IM
        IF (J == JM) IMAX=1
        DO I=IMIN,IMAX
      CM(0) = 0.
       C(0) = 0.
      FM(0) = 0.
      FX(0) = 0.
      FY(0) = 0.
      FZ(0) = 0.
      LMIJ=LMM(I,J)
      IF(LMIJ.LE.1)  GO TO 330
      CM(LMIJ) = 0.
       C(LMIJ) = 0.
      FM(LMIJ) = 0.
      FX(LMIJ) = 0.
      FY(LMIJ) = 0.
      FZ(LMIJ) = 0.
C****
C**** Calculate FM (kg), FX (kg), FY (kg) and FZ (kg**2)
C****
      DO 120 L=1,LMIJ-1
      CM(L) = DT*MW(I,J,L)
      IF(CM(L).LT.0.)  GO TO 110
C**** Ocean mass flux is positive
      C(L)  = CM(L)/MO(I,J,L)
      IF(C(L).GT.1d0)  WRITE (6,*) 'C>1:',I,L,C(L),MO(I,J,L)
      FM(L) = C(L)*(RM(I,J,L)+(1d0-C(L))*RZ(I,J,L))
      FX(L) = C(L)*RX(I,J,L)
      FY(L) = C(L)*RY(I,J,L)
      FZ(L) = CM(L)*(C(L)*C(L)*RZ(I,J,L)-3d0*FM(L))
      GO TO 120
C**** Ocean mass flux is negative
  110 C(L)  = CM(L)/MO(I,J,L+1)
      IF(C(L).LT.-1d0)  WRITE (6,*) 'C<-1:',I,L,C(L),MO(I,J,L+1)
      FM(L) = C(L)*(RM(I,J,L+1)-(1d0+C(L))*RZ(I,J,L+1))
      FX(L) = C(L)*RX(I,J,L+1)
      FY(L) = C(L)*RY(I,J,L+1)
      FZ(L) = CM(L)*(C(L)*C(L)*RZ(I,J,L+1)-3d0*FM(L))
  120 CONTINUE
C****
C**** Modify the tracer moments so that the tracer mass in each
C**** division is non-negative
C****
      IF(.NOT.QLIMIT)  GO TO 300
      DO 290 L=1,LMIJ
      IF(C(L-1).GE.0.)  GO TO 240
C**** Water is leaving through the bottom edge: 2 or 3 divisions
      IF(FM(L-1).LE.0.)  GO TO 210
C**** Bottom most division is negative, RMB = -FM(L-1) < 0: Case 2 or 4
      RZ(I,J,L) = RM(I,J,L)/(1d0+C(L-1))
      FM(L-1) = 0.
      FZ(L-1) = CM(L-1)*C(L-1)*C(L-1)*RZ(I,J,L)
      IF(C(L).LE.0.)  GO TO 290
      FM(L) = C(L)*(RM(I,J,L)+(1d0-C(L))*RZ(I,J,L))
      FZ(L) = CM(L)*(C(L)*C(L)*RZ(I,J,L)-3d0*FM(L))
      GO TO 290
C**** Bottom most division is non-negative, RMB = -FM(L-1) > 0:
C**** Case 1, 3 or 5
  210 IF(C(L).LE.0.)  GO TO 230
C**** Water is leaving through the top edge: 3 divisions
      IF(FM(L).GE.0.)  GO TO 290
C**** Top most division is negative, RMT = FM(L) < 0: Case 3 or 5
      RZ(I,J,L) = -RM(I,J,L)/(1d0-C(L))
      FM(L) = 0.
      FZ(L) = CM(L)*C(L)*C(L)*RZ(I,J,L)
      FM(L-1) = C(L-1)*(RM(I,J,L)-(1d0+C(L-1))*RZ(I,J,L))
      FZ(L-1) = CM(L-1)*(C(L-1)*C(L-1)*RZ(I,J,L)-3d0*FM(L-1))
      GO TO 290
C**** No water is leaving through the top edge: 2 divisions
  230 IF(RM(I,J,L)+FM(L-1).GE.0.)  GO TO 290
C**** Top most division is negative, RMT = RM(I,J,L)+FM(L-1) < 0: Case 3
      RZ(I,J,L) = RM(I,J,L)/C(L-1)
      FM(L-1) = -RM(I,J,L)
      FZ(L-1) = CM(L-1)*(C(L-1)+3d0)*RM(I,J,L)
      GO TO 290
C**** No water is leaving through the bottom edge: 1 or 2 divisions
  240 IF(C(L).LE.0.)  GO TO 290
C**** Water is leaving through the top edge: 2 divisions
      IF(FM(L).GE.0.)  GO TO 250
C**** Top most division is negative, RMT = FM(L) < 0: Case 3
      RZ(I,J,L) = -RM(I,J,L)/(1d0-C(L))
      FM(L) = 0.
      FZ(L) = CM(L)*C(L)*C(L)*RZ(I,J,L)
      GO TO 290
C**** Top most division is non-negative, RMT = FM(L) > 0: Case 1 or 2
  250 IF(RM(I,J,L)-FM(L).GE.0.)  GO TO 290
C**** Bottom most division is negative, RMB = RM(I,J,L)-FM(L) < 0: Cas 2
      RZ(I,J,L) = RM(I,J,L)/C(L)
      FM(L) = RM(I,J,L)
      FZ(L) = CM(L)*(C(L)-3d0)*RM(I,J,L)
C****
  290 CONTINUE
C****
C**** Calculate new tracer mass and first moments of tracer mass
C****
  300 DO 310 L=1,LMIJ
      RM(I,J,L) = RM(I,J,L) + (FM(L-1)-FM(L))
      RX(I,J,L) = RX(I,J,L) + (FX(L-1)-FX(L))
      RY(I,J,L) = RY(I,J,L) + (FY(L-1)-FY(L))
      RZ(I,J,L) = (RZ(I,J,L)*MO(I,J,L) + (FZ(L-1)-FZ(L))
     *  + 3d0*((CM(L-1)+CM(L))*RM(I,J,L)-MO(I,J,L)*(FM(L-1)+FM(L))))
     *  / (MO(I,J,L)+CM(L-1)-CM(L))
C****
      if ( QLIMIT ) then ! limit tracer gradients
        RXY = abs(RX(I,J,L)) + abs(RY(I,J,L))
        if ( RXY > RM(I,J,L) ) then
          RX(I,J,L) = RX(I,J,L)*( RM(I,J,L)/(RXY + tiny(RXY)) )
          RY(I,J,L) = RY(I,J,L)*( RM(I,J,L)/(RXY + tiny(RXY)) )
        end if
        if ( abs(RZ(I,J,L)) > RM(I,J,L) )
     *       RZ(I,J,L) = sign(RM(I,J,L), RZ(I,J,L)+0d0)
      end if
C****
      MO(I,J,L) = MO(I,J,L) +  CM(L-1)-CM(L)
         IF(MO(I,J,L).LE.0.)              ICKERR=ICKERR+1
         IF(QLIMIT.AND.RM(I,J,L).LT.0.)   ICKERR=ICKERR+1
  310 CONTINUE
         DO 320 L=1,LMIJ-1
  320    OIJL(I,J,L) = OIJL(I,J,L) + FM(L)
  330 CONTINUE
      END DO
      END DO
!$OMP END PARALLEL DO

C**** IF NO ERROR HAS OCCURRED - RETURN, ELSE STOP
      IF(ICKERR == 0)  RETURN
      DO J=J_0,J_1
        IMIN=1
        IMAX=IM
        IF (J == 1) IMIN=IM
        IF (J == JM) IMAX=1
        DO I=IMIN,IMAX
        LMIJ=LMM(I,J)
        DO L=1,LMIJ
          IF(FOCEAN(I,J).gt.0 .and. MO(I,J,L).LE.0.)  GO TO 800
          IF(QLIMIT .AND. RM(I,J,L).LT.0.) GO TO 810
        END DO
      END DO
      END DO
      WRITE(6,*) 'ERROR CHECK INCONSISTENCY: OADVTZ ',ICKERR
      call stop_model("OAVDTZ",255)

  800 WRITE (6,*) 'MO<0 in OADVTZ:',I,J,L,MO(I,J,L)
  810 WRITE (6,*) 'RM in OADVTZ:',I,J,L,RM(I,J,L)
c      WRITE (6,*) 'C=',(L,C(L),L=0,LMIJ)
      call stop_model("OADVTZ",255)
      END SUBROUTINE OADVTZ

      Subroutine OBDRAG
!@sum  OBDRAG exerts a drag on the Ocean Model's bottom layer
!@auth Gary Russell
!@ver  1.0
      Use OCEAN, Only: im,jm,lmo,IVNP,J1O, mo,uo,vo, lmu,lmv, dts,
     *                 COSI=>COSIC,SINI=>SINIC
      use domain_decomp, only : grid, get, halo_update, north, south

      Implicit None
      REAL*8, PARAMETER :: BDRAGX=1d0, SDRAGX=1d-1
      REAL*8,DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO):: UT,VT
      Integer*4 I,J,L,Ip1,Im1
      REAL*8 WSQ

      INTEGER :: J_0, J_0S,J_1S  ; logical :: have_north_pole

      CALL GET(grid, J_STRT=J_0, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     *               have_north_pole=have_north_pole)
C****
C**** UO = UO*(1-x/y)  is approximated by  UO*y/(y+x)  for stability
C****
      call halo_update (grid, mo, FROM=NORTH)
      call halo_update (grid, uo, FROM=NORTH)
      call halo_update (grid, vo, FROM=SOUTH)
C**** Save UO,VO into UT,VT which will be unchanged
!$OMP ParallelDo   Private (L)
      DO 10 L=1,LMO
        UT(:,:,L) = UO(:,:,L)
 10     VT(:,:,L) = VO(:,:,L)
!$OMP EndParallelDo
C****
C**** Reduce West-East ocean current
C****
C**** Bottom drag in the interior
!$OMP ParallelDo   Private (I,J,L,Ip1, WSQ)
      DO 120 J=max(J1O,J_0),J_1S
      I=IM
      DO 110 IP1=1,IM
      IF(LMU(I,J) <= 0)  GO TO 110
      L=LMU(I,J)
      WSQ = UT(I,J,L)*UT(I,J,L) + 1d-20 +
     *  .25*(VT(I,J  ,L)*VT(I,J  ,L) + VT(IP1,J  ,L)*VT(IP1,J  ,L)
     *     + VT(I,J-1,L)*VT(I,J-1,L) + VT(IP1,J-1,L)*VT(IP1,J-1,L))
      UO(I,J,L) = UO(I,J,L) * (MO(I,J,L)+MO(IP1,J,L)) /
     *           (MO(I,J,L)+MO(IP1,J,L) + DTS*BDRAGX*SQRT(WSQ)*2d0)
  110 I=IP1
  120 CONTINUE
C**** Bottom drag at the poles
C     IF(LMU(1,1 or JM) <= 0)  GO TO
      if (have_north_pole) then
        L=LMU(1,JM)
        WSQ = UT(IM,JM,L)**2 + UT(IVNP,JM,L)**2 + 1d-20
        UO(IM,JM,L) = UO(IM,JM,L) * MO(1,JM,L) /
     /                           (MO(1,JM,L) + DTS*BDRAGX*Sqrt(WSQ))
        UO(IVNP,JM,L) = UO(IVNP,JM,L) * MO(1,JM,L) /
     /                           (MO(1,JM,L) + DTS*BDRAGX*Sqrt(WSQ))
      end if
C****
C**** Reduce South-North ocean current
C****
!$OMP ParallelDo   Private (I,J,L,Im1, WSQ)
      Do 240 J=max(J1O,J_0),J_1S
      If (J==JM-1)  GoTo 220
C**** Bottom drag away from north pole
      IM1=IM
      DO 210 I=1,IM
      IF(LMV(I,J) <= 0)  GO TO 210
      L=LMV(I,J)
      WSQ = VT(I,J,L)*VT(I,J,L) + 1d-20 +
     *  .25*(UT(IM1,J+1,L)*UT(IM1,J+1,L) + UT(I,J+1,L)*UT(I,J+1,L)
     *     + UT(IM1,J  ,L)*UT(IM1,J  ,L) + UT(I,J  ,L)*UT(I,J  ,L))
      VO(I,J,L) = VO(I,J,L) * (MO(I,J,L)+MO(I,J+1,L)) /
     *           (MO(I,J,L)+MO(I,J+1,L) + DTS*BDRAGX*SQRT(WSQ)*2d0)
  210 IM1=I
      GoTo 240
C**** Bottom drag near north pole
  220 Im1=IM
      Do 230 I=1,IM
      If (LMV(I,JM-1) <= 0)  GoTo 230
      L = LMV(I,JM-1)
      WSQ = VT(I,JM-1,L)*VT(I,JM-1,L) + 1d-20 +
     +   .5*(UT(IM,JM,L)*COSI(I) + UT(IVNP,JM,L)*SINI(I))**2 +
     +  .25*(UT(Im1,JM-1,L)*UT(Im1,JM-1,L) + UT(I,JM-1,L)*UT(I,JM-1,L))
      VO(I,JM-1,L) = VO(I,JM-1,L) * (MO(I,JM-1,L)+MO(1,JM,L)) /
     *              (MO(I,JM-1,L)+MO(1,JM,L) + DTS*BDRAGX*SQRT(WSQ)*2)
  230 Im1=I
  240 Continue
      RETURN
C****
      END Subroutine OBDRAG

      SUBROUTINE OCOAST
!@sum OCOAST reduces the horizontal perpendicular gradients of tracers
!@sum in coastline ocean grid boxes
!@auth Gary Russell
!@ver  1.0
      USE CONSTANT, only : sday
      USE OCEAN, only : im,jm,dts,lmm,gxmo,gymo,sxmo,symo
#ifdef TRACERS_OCEAN
     *     ,ntm,txmo,tymo
#endif
      use domain_decomp, only : grid, get

      IMPLICIT NONE
      INTEGER I,IM1,IP1,J,LMIN,L,N
      REAL*8 REDUCE

      integer :: J_0S, J_1S

      call get (grid, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S)


      REDUCE = 1d0 - DTS/(SDAY*2d1)
C**** Reduce West-East gradient of tracers
      IM1=IM-1
      I=IM
      DO 120 J=J_0S,J_1S
      DO 120 IP1=1,IM
      LMIN = MIN(LMM(IM1,J),LMM(IP1,J)) + 1
      DO 110 L=LMIN,LMM(I,J)
        GXMO(I,J,L) = GXMO(I,J,L)*REDUCE
        SXMO(I,J,L) = SXMO(I,J,L)*REDUCE
#ifdef TRACERS_OCEAN
        DO N = 1,NTM
          TXMO(I,J,L,N) = TXMO(I,J,L,N) *REDUCE
        END DO
#endif
 110  CONTINUE
      IM1=I
  120 I=IP1
C**** Reduce South-North gradient of tracers
      DO 220 J=J_0S,J_1S
      DO 220 I=1,IM
      LMIN = MIN(LMM(I,J-1),LMM(I,J+1)) + 1
      DO 210 L=LMIN,LMM(I,J)
        GYMO(I,J,L) = GYMO(I,J,L)*REDUCE
        SYMO(I,J,L) = SYMO(I,J,L)*REDUCE
#ifdef TRACERS_OCEAN
        DO N = 1,NTM
          TYMO(I,J,L,N) = TYMO(I,J,L,N) *REDUCE
        END DO
#endif
 210  CONTINUE
  220 CONTINUE
      RETURN
      END SUBROUTINE OCOAST

      SUBROUTINE OSTRES
!@sum OSTRES applies the atmospheric surface stress over open ocean
!@sum and the sea ice stress to the layer 1 ocean velocities
!@auth Gary Russell
!@ver  1.0
      USE OCEAN, only : im,jm,ivnp,uo,vo,mo,dxyso,dxyno,dxyvo,
     *                  lmu,lmv,cosic,sinic,ratoc
      USE FLUXES, only : dmua,dmva,dmui,dmvi
      use domain_decomp, only : grid, get, halo_update, north

      IMPLICIT NONE
      INTEGER I,J,IP1
C****
C**** All stress now defined for whole box, not just ocn or ice fraction
C**** FLUXCB  DMUA(1)  U momentum downward into open ocean (kg/m*s)
C****         DMVA(1)  V momentum downward into open ocean (kg/m*s)
C****         DMUA(2,JM,1)  polar atmo. mass slowed to zero (kg/m**2)
C****         DMUI     U momentum downward from sea ice (kg/m*s)
C****         DMVI     V momentum downward from sea ice (kg/m*s)

      integer :: J_0, J_1, J_0S, J_1S  ; logical :: have_north_pole

      call get (grid, J_STRT=J_0, J_STOP=J_1,
     *                J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     *                have_north_pole=have_north_pole)

C**** Scale stresses for ocean area
      DO J=J_0,J_1
        DO I=1,IM
          DMUA(I,J,1)=RATOC(J)*DMUA(I,J,1)
          DMVA(I,J,1)=RATOC(J)*DMVA(I,J,1)
          DMUI(I,J)=RATOC(J)*DMUI(I,J)
          DMVI(I,J)=RATOC(J)*DMVI(I,J)
        END DO
      END DO
C****
C**** Surface stress is applied to U component
C****
      I=IM
      DO J=J_0S,J_1S
      DO IP1=1,IM
        IF(LMU(I,J).gt.0.)  UO(I,J,1) = UO(I,J,1) +
     *       (DMUA(I,J,1) + DMUA(IP1,J,1) + 2d0*DMUI(I,J)) /
     *       (  MO(I,J,1) +   MO(IP1,J,1))
        I=IP1
      END DO
      END DO
      if (have_north_pole) then
        UO(IM  ,JM,1) = UO(IM  ,JM,1) + DMUA(1,JM,1)/MO(1,JM,1)
        UO(IVNP,JM,1) = UO(IVNP,JM,1) + DMVA(1,JM,1)/MO(1,JM,1)
      end if
C****
C**** Surface stress is applied to V component
C****
      call halo_update(grid, dmva, from=north)
      call halo_update(grid,   mo, from=north)
      DO J=J_0S,min(J_1S,JM-2)
      DO I=1,IM
        IF(LMV(I,J).GT.0.)  VO(I,J,1) = VO(I,J,1) +
     *       (DMVA(I,J  ,1)*DXYNO(J) + DMVA(I,J+1,1)*DXYSO(J+1)
     *       +2d0*DMVI(I,J)*DXYVO(J))
     * / (MO(I,J,1)*DXYNO(J) + MO(I,J+1,1)*DXYSO(J+1))
      END DO
      END DO
C**** Surface stress is applied to V component at the North Pole
      if (have_north_pole) then
      DO I=1,IM
        VO(I,JM-1,1) = VO(I,JM-1,1) +
     *       (DMVA(I,JM-1,1)*DXYNO(JM-1)+
     *       (DMVA(1,JM,1)*COSIC(I) - DMUA(1,JM,1)*SINIC(I))*DXYSO(JM)
     *       + 2d0*DMVI(I,JM-1)*DXYVO(JM-1)) /
     *  (MO(I,JM-1,1)*DXYNO(JM-1) + MO(I,JM,1)*DXYSO(JM))
      END DO
      end if
      RETURN
      END SUBROUTINE OSTRES

      SUBROUTINE GROUND_OC
!@sum  GROUND_OC adds vertical fluxes into the ocean
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : grav
      USE GEOM, only : dxyp,bydxyp
      USE OCEAN, only : im,jm,mo,g0m,s0m,focean,gzmo,imaxj,dxypo,bydxypo
     *     ,lmo,lmm,ratoc,rocat,opress
#ifdef TRACERS_OCEAN
     *     ,trmo,ntm
#endif
      USE FLUXES, only : solar,e0,evapor,dmsi,dhsi,dssi,runosi,erunosi
     *     ,flowo,eflowo,srunosi,apress,melti,emelti,smelti
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
     *     ,trflowo,trevapor,trunosi,trmelti
#ifdef TRACERS_DRYDEP
     *     ,trdrydep
#endif
#endif
     *     ,dtrsi
#endif
#ifdef TRACERS_GASEXCH_Natassa
      USE FLUXES, only : TRGASEX
#endif
      USE SEAICE_COM, only : rsi
      use domain_decomp, only : grid, get

      IMPLICIT NONE
      INTEGER I,J
      REAL*8 DXYPJ,BYDXYPJ,RUNO,RUNI,ERUNO,ERUNI,SROX(2),G0ML(LMO)
     *     ,MO1,SO1,ROICE,DMOO,DMOI,DEOO,DEOI,GZML(LMO),SRUNO,SRUNI,DSOO
     *     ,DSOI,POCEAN,POICE
#ifdef TRACERS_OCEAN
      REAL*8, DIMENSION(NTM) :: TRUNO,TRUNI,DTROO,DTROI,TRO1
#endif

      integer ::  J_1, J_0S
      logical :: have_south_pole, have_north_pole

      call get (grid, J_STOP=J_1, J_STRT_SKP=J_0S,
     *                have_north_pole=have_north_pole,
     *                have_south_pole=have_south_pole)

C****
C**** Add surface source of fresh water and heat
C****
      DO J=J_0S,J_1
        DXYPJ=DXYPO(J)
        BYDXYPJ=BYDXYPO(J)
      DO I=1,IMAXJ(J)
      IF(FOCEAN(I,J).gt.0.) THEN
        ROICE = RSI(I,J)
        POCEAN=FOCEAN(I,J)*(1.-ROICE)
        POICE =FOCEAN(I,J)*ROICE
        DXYPJ=DXYPJ*FOCEAN(I,J)     ! adjust areas for completeness
        BYDXYPJ=BYDXYPJ/FOCEAN(I,J) ! no change to results
C**** set mass & energy fluxes (incl. river/sea ice runoff + basal flux)
        RUNO = (FLOWO(I,J)+ MELTI(I,J))/(DXYPJ)-
     *                                   RATOC(J)*EVAPOR(I,J,1)
        RUNI = (FLOWO(I,J)+ MELTI(I,J))/(DXYPJ)+
     *                                   RATOC(J)*RUNOSI(I,J)
        ERUNO=(EFLOWO(I,J)+EMELTI(I,J))/(DXYPJ)+
     *                                   RATOC(J)*E0(I,J,1)
        ERUNI=(EFLOWO(I,J)+EMELTI(I,J))/(DXYPJ)+
     *                                   RATOC(J)*ERUNOSI(I,J)
        SRUNO=SMELTI(I,J)/(DXYPJ)
        SRUNI=SMELTI(I,J)/(DXYPJ)+RATOC(J)*SRUNOSI(I,J)
        G0ML(:) =  G0M(I,J,:)
        GZML(:) = GZMO(I,J,:)
        SROX(1)=SOLAR(1,I,J)*RATOC(J) ! open water
        SROX(2)=SOLAR(3,I,J)*RATOC(J) ! through ice
        MO1 = MO(I,J,1)
        SO1 = S0M(I,J,1)
#ifdef TRACERS_OCEAN
        TRO1(:) = TRMO(I,J,1,:)
#ifdef TRACERS_WATER
        TRUNO(:)=(TRFLOWO(:,I,J)+TRMELTI(:,I,J))/(DXYPJ)-
     *       RATOC(J)*TREVAPOR(:,1,I,J)
#ifdef TRACERS_DRYDEP
     *       + RATOC(J)*trdrydep(:,1,i,j)
#endif
        TRUNI(:)=(TRFLOWO(:,I,J)+TRMELTI(:,I,J))/(DXYPJ)+
     *       RATOC(J)*TRUNOSI(:,I,J)
#else
        TRUNO(:)=0. ; TRUNI(:)=0.
#endif

#ifdef TRACERS_GASEXCH_Natassa
        !note that TRGASEX is positive down i.e. same sign as
        !TRUNO, while TREVAPOR is positive up.
        TRUNO(:)=RATOC(J)*TRGASEX(:,1,I,J)
#endif

#endif

        CALL OSOURC(ROICE,MO1,G0ML,GZML,SO1,DXYPJ,BYDXYPJ,LMM(I,J),RUNO
     *         ,RUNI,ERUNO,ERUNI,SRUNO,SRUNI,SROX,
#ifdef TRACERS_OCEAN
     *         TRO1,TRUNO,TRUNI,DTROO,DTROI,
#endif
     *         DMOO,DEOO,DMOI,DEOI,DSOO,DSOI)

C**** update ocean variables
          MO(I,J,1) = MO1
         S0M(I,J,1) = SO1
         G0M(I,J,:) = G0ML(:)
        GZMO(I,J,:) = GZML(:)

C**** Store mass and energy fluxes for formation of sea ice
        DMSI(1,I,J)=DMOO*ROCAT(J)
        DMSI(2,I,J)=DMOI*ROCAT(J)
        DHSI(1,I,J)=DEOO*ROCAT(J)
        DHSI(2,I,J)=DEOI*ROCAT(J)
        DSSI(1,I,J)=DSOO*ROCAT(J)
        DSSI(2,I,J)=DSOI*ROCAT(J)
#ifdef TRACERS_OCEAN
        TRMO(I,J,1,:) = TRO1(:)
        DTRSI(:,1,I,J)=DTROO(:)*ROCAT(J)
        DTRSI(:,2,I,J)=DTROI(:)*ROCAT(J)
#endif

C**** Calculate pressure anomaly at ocean surface (and scale for areas)
C**** Updated using latest sea ice (this ensures that total column mass
C**** is consistent for OGEOZ calculation).
        OPRESS(I,J) = RATOC(J)*(APRESS(I,J))+GRAV*(
     *       (1.-RSI(I,J))*DMSI(1,I,J) + RSI(I,J)*DMSI(2,I,J))

        END IF
      END DO
      END DO

      if(have_south_pole) OPRESS(2:IM,1)  = OPRESS(1,1)
      if(have_north_pole) OPRESS(2:IM,JM) = OPRESS(1,JM)
C****
      RETURN
      END SUBROUTINE GROUND_OC

      SUBROUTINE OSOURC (ROICE,MO,G0ML,GZML,S0M,DXYPJ,BYDXYPJ,LMIJ,RUNO
     *     ,RUNI,ERUNO,ERUNI,SRUNO,SRUNI,SROX,
#ifdef TRACERS_OCEAN
     *     TROM,TRUNO,TRUNI,DTROO,DTROI,
#endif
     *     DMOO,DEOO,DMOI,DEOI,DSOO,DSOI)
!@sum  OSOURC applies fluxes to ocean in ice-covered and ice-free areas
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : ntm,trname
#endif
      USE SW2OCEAN, only : lsrpd,fsr,fsrz
      USE SEAICE, only : fsss, Ei
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: ROICE,DXYPJ,BYDXYPJ,RUNO,RUNI,ERUNO,ERUNI
     *     ,SROX(2),SRUNO,SRUNI
      INTEGER, INTENT(IN) :: LMIJ
      REAL*8, INTENT(INOUT) :: MO,G0ML(LSRPD),GZML(LSRPD),S0M
      REAL*8, INTENT(OUT) :: DMOO,DMOI,DEOO,DEOI,DSOO,DSOI
#ifdef TRACERS_OCEAN
      REAL*8, DIMENSION(NTM), INTENT(INOUT) :: TROM
      REAL*8, DIMENSION(NTM), INTENT(IN) :: TRUNO,TRUNI
      REAL*8, DIMENSION(NTM), INTENT(OUT) :: DTROO,DTROI
      REAL*8, DIMENSION(NTM) :: TMOO,TMOI,FRAC
#ifdef TRACERS_SPECIAL_O18
      REAL*8 fracls
#endif
#endif
      REAL*8 MOO,GOO,GMOO,GMOI,MOI,GOI,SMOO,SMOI,SOO,SOI,GFOO,GFOI,TFOO
     *     ,TFOI,SIOO,SIOI
      REAL*8 GFREZS,TFREZS,TSOL
      INTEGER L,LSR,N

      DMOO=0. ; DEOO=0. ; DMOI=0. ; DEOI=0. ; DSOO=0. ; DSOI=0.
#ifdef TRACERS_OCEAN
      DTROI(:) = 0. ; DTROO(:) = 0.
      do n=1,ntm
#ifdef TRACERS_SPECIAL_O18
        FRAC(n)=fracls(n)
#else
        FRAC(n)=1.
#endif
      end do
#endif

      LSR = MIN(LSRPD,LMIJ)
C****
C**** Open Ocean
C****
      MOO  = MO + RUNO
      GMOO = G0ML(1)*BYDXYPJ + ERUNO
      SMOO = S0M*BYDXYPJ + SRUNO
#ifdef TRACERS_OCEAN
      TMOO(:) = TROM(:)*BYDXYPJ+TRUNO(:)
#endif
      IF (ROICE.lt.1d0) THEN
C**** Remove insolation from layer 1 that goes to lower layers
      IF (LSR.gt.1) GMOO = GMOO - SROX(1)*FSR(2)

      GOO  = GMOO/MOO
      SOO  = SMOO/MOO
      GFOO = GFREZS(SOO)
      IF(GOO.lt.GFOO) THEN
C**** Open ocean is below freezing, calculate
C**** DMOO = mass of ocean that freezes over open fraction from
C**** GOO*MOO = GFOO*(MOO-DMOO) + Ei(TFOO,SIOO*1d3)*DMOO
        TFOO = TFREZS(SOO)
        SIOO = FSSS*SOO
        DMOO = MOO*(GOO-GFOO)/(Ei(TFOO,SIOO*1d3)-GFOO)
        DEOO = Ei(TFOO,SIOO*1d3)*DMOO
        DSOO = SIOO*DMOO
#ifdef TRACERS_OCEAN
        DTROO(:) = TMOO(:)*FRAC(:)*(DMOO-DSOO)/(MOO-SMOO)
#endif
      END IF
      END IF
C****
C**** Ocean underneath the ice
C****
      MOI  = MO + RUNI
      GMOI = G0ML(1)*BYDXYPJ + ERUNI
      SMOI = S0M*BYDXYPJ + SRUNI
#ifdef TRACERS_OCEAN
      TMOI(:) = TROM(:)*BYDXYPJ+TRUNI(:)
#endif
      IF(ROICE.gt.0.) THEN
C**** Remove insolation from layer 1 that goes to lower layers
        IF (LSR.gt.1) GMOI = GMOI - SROX(2)*FSR(2)

        GOI  = GMOI/MOI
        SOI  = SMOI/MOI
        GFOI = GFREZS(SOI)
        IF(GOI.LT.GFOI) THEN
C**** Ocean underneath the ice is below freezing, calculate
C**** DMOI = mass of ocean that freezes under sea ice fraction from
C**** GOI*MOI = GFOI*(MOI-DMOI) + Ei(TFOI,SIOI*1d3)*DMOI
          TFOI = TFREZS(SOI)
          SIOI = FSSS*SOI
          DMOI = MOI*(GOI-GFOI)/(Ei(TFOI,SIOI*1d3)-GFOI)
          DEOI = Ei(TFOI,SIOI*1d3)*DMOI
          DSOI = SIOI*DMOI
#ifdef TRACERS_OCEAN
          DTROI(:) = TMOI(:)*FRAC(:)*(DMOI-DSOI)/(MOI-SMOI)
#endif
        END IF
      END IF
C**** Update first layer variables
      MO     =  (MOI-DMOI)*ROICE + (1.-ROICE)*( MOO-DMOO)
      G0ML(1)=((GMOI-DEOI)*ROICE + (1.-ROICE)*(GMOO-DEOO))*DXYPJ
      S0M    =((SMOI-DSOI)*ROICE + (1.-ROICE)*(SMOO-DSOO))*DXYPJ
#ifdef TRACERS_OCEAN
      TROM(:)=((TMOI(:)-DTROI(:))*ROICE + (1.-ROICE)*(TMOO(:)-DTROO(:)))
     *     *DXYPJ
#endif
C**** add insolation to lower layers
      TSOL=(SROX(1)*(1.-ROICE)+SROX(2)*ROICE)*DXYPJ
      DO L=2,LSR-1
        G0ML(L)=G0ML(L)+TSOL*(FSR(L)-FSR(L+1))
        GZML(L)=GZML(L)+TSOL*FSRZ(L)
      END DO
      G0ML(LSR) = G0ML(LSR) + TSOL*FSR (LSR)
      GZML(LSR) = GZML(LSR) + TSOL*FSRZ(LSR)
C****
      RETURN
      END SUBROUTINE OSOURC

      SUBROUTINE PRECIP_OC
!@sum  PRECIP_OC driver for applying precipitation to ocean fraction
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE GEOM, only : dxyp
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : ntm
#ifdef TRACERS_WATER
      USE FLUXES, only : trpreca=>trprec,trunpsia=>trunpsi
#endif
#endif
      USE FLUXES, only : runpsia=>runpsi,srunpsia=>srunpsi,preca=>prec
     *     ,epreca=>eprec, erunpsia=>erunpsi 
      USE OCEAN, only : im,jm,mo,g0m,s0m,bydxypo,focean,imaxj
#ifdef TRACERS_OCEAN
     *     ,trmo,dxypo
#endif
      USE SEAICE_COM, only : rsia=>rsi
      USE DOMAIN_DECOMP, only : grid,get

      IMPLICIT NONE
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *   PREC,EPREC,RUNPSI,RSI,SRUNPSI,ERUNPSI
#ifdef TRACERS_OCEAN
      REAL*8, DIMENSION(NTM,IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *   trprec,trunpsi
#endif
      INTEGER I,J

      integer :: J_0, J_1

      call get (grid, J_STRT=J_0, J_STOP=J_1)

C**** save surface variables before any fluxes are added
      CALL KVINIT

C**** Convert fluxes on atmospheric grid to oceanic grid
C**** build in enough code to allow a different ocean grid.
C**** Since the geometry differs on B and C grids, some processing
C**** of fluxes is necessary anyway
      DO J=J_0,J_1
        DO I=1,IMAXJ(J)
          PREC   (I,J)=PRECA   (I,J)*DXYP(J)*BYDXYPO(J)  ! kg/m^2
          EPREC  (I,J)=EPRECA  (I,J)*DXYP(J)             ! J
          RUNPSI (I,J)=RUNPSIA (I,J)*DXYP(J)*BYDXYPO(J)  ! kg/m^2
          SRUNPSI(I,J)=SRUNPSIA(I,J)*DXYP(J)             ! kg
          ERUNPSI(I,J)=ERUNPSIA(I,J)*DXYP(J)             ! J
          RSI    (I,J)=RSIA    (I,J)
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
          TRPREC(:,I,J)=TRPRECA(:,I,J)                   ! kg
          TRUNPSI(:,I,J)=TRUNPSIA(:,I,J)*DXYP(J)         ! kg
#else
          TRPREC(:,I,J)=0.  ; TRUNPSI(:,I,J)=0.
#endif
#endif
        END DO
      END DO
C****
      DO J=J_0,J_1
        DO I=1,IMAXJ(J)
          IF(FOCEAN(I,J).gt.0. .and. PREC(I,J).gt.0.)  THEN
            MO (I,J,1)= MO(I,J,1) + ((1d0-RSI(I,J))*PREC(I,J) +
     *           RSI(I,J)*RUNPSI(I,J))*FOCEAN(I,J)
            G0M(I,J,1)=G0M(I,J,1)+ ((1d0-RSI(I,J))*EPREC(I,J) +
     *           RSI(I,J)*ERUNPSI(I,J))*FOCEAN(I,J)
            S0M(I,J,1)=S0M(I,J,1) + RSI(I,J)*SRUNPSI(I,J)*FOCEAN(I,J)
#ifdef TRACERS_OCEAN
            TRMO(I,J,1,:)=TRMO(I,J,1,:)+((1d0-RSI(I,J))*TRPREC(:,I,J)
     *             +RSI(I,J)*TRUNPSI(:,I,J))*FOCEAN(I,J)
#endif
          END IF
        END DO
      END DO

C**** Convert ocean surface temp to atmospheric SST array
      CALL TOC2SST

      RETURN
      END SUBROUTINE PRECIP_OC

      SUBROUTINE ODIFF (DTDIFF)
C???? ESMF-exception - ODIFF currently works with global arrays
!@sum  ODIFF applies Wasjowicz horizontal viscosity to velocities
!@auth Gavin Schmidt
!@ver  1.0
C****
C**** ODIFF calculates horizontal Wasjowicz viscosity terms in momentum
C**** equations implicitly using ADI method and assumes no slip/free
C**** slip conditions at the side. K_h (m^2/s) may vary spatially
C**** based on Munk length though must remain isotropic.
C**** (If longitudinal variation is wanted just make K arrays K(I,J))
C**** FSLIP = 0 implies no slip conditions, = 1 implies free slip
C**** Mass variation is included
C****
      USE CONSTANT, only :twopi,rhows,omega,radius
      USE OCEAN, only : im,jm,lmo,mo,uo,vo,
     *  IVNP,UONP,VONP, COSU,SINU, COSI=>COSIC,SINI=>SINIC,
     *  cospo,cosvo,rlat,lmu,lmv,dxpo,dypo,dxvo,dyvo,dxyvo,dxypo,bydxypo
      USE OCEAN_DYN, only : dh
      USE TRIDIAG_MOD, only : tridiag, tridiag_new

      USE DOMAIN_DECOMP, ONLY : grid, GET, AM_I_ROOT
      USE DOMAIN_DECOMP, ONLY : HALO_UPDATE, NORTH, SOUTH, ESMF_BCAST
      USE DOMAIN_DECOMP, ONLY : haveLatitude

      IMPLICIT NONE
      REAL*8, PARAMETER :: AKHMIN=1.5d8, FSLIP=0.

      REAL*8, DIMENSION(grid%j_strt_halo:grid%j_stop_halo) ::
     *     KYPXP,KXPYV,KYVXV,KXVYP
      REAL*8, SAVE, ALLOCATABLE, DIMENSION(:) ::
     *     BYDXYV, KHP,KHV,TANP,TANV,BYDXV,BYDXP,BYDYV,BYDYP
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,2) ::
     *      DUDX,DUDY,DVDX,DVDY
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo) ::
     *      FUX,FUY,FVX,FVY,BYMU,BYMV

      INTEGER, SAVE :: IFIRST = 1
C**** Local variables

      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::
     *     AU, BU, CU, RU, UU
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::
     *     AV, BV, CV, RV, UV
      REAL*8, SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::
     *      UXA,UXB,UXC,UYA,UYB,UYC,VXA,VXB,VXC,VYA,VYB,VYC
      REAL*8, SAVE, DIMENSION(LMO) :: UYPB
      REAL*8, SAVE, DIMENSION(IM,LMO) :: UYPA

      REAL*8, INTENT(IN) :: DTDIFF
      REAL*8, SAVE :: BYDXYPJM
      REAL*8 DSV,DSP,VLAT,DLAT,DT2,DTU,DTV,VX,VY,VT,UT,UX,UY
      INTEGER I,J,L,IP1,IM1,II

!     domain decomposition
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0H, J_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

c**** Extract domain decomposition info
      CALL GET(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

C****

      IF(IFIRST.ne.0)  THEN
      IFIRST = 0

        allocate( BYDXYV(grid%j_strt_halo:grid%j_stop_halo) )
        allocate( KHP   (grid%j_strt_halo:grid%j_stop_halo) )
        allocate( KHV   (grid%j_strt_halo:grid%j_stop_halo) )
        allocate( TANP  (grid%j_strt_halo:grid%j_stop_halo) )
        allocate( TANV  (grid%j_strt_halo:grid%j_stop_halo) )
        allocate( BYDXV (grid%j_strt_halo:grid%j_stop_halo) )
        allocate( BYDXP (grid%j_strt_halo:grid%j_stop_halo) )
        allocate( BYDYV (grid%j_strt_halo:grid%j_stop_halo) )
        allocate( BYDYP (grid%j_strt_halo:grid%j_stop_halo) )

        allocate( UXA(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( UXB(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( UXC(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( UYA(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( UYB(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( UYC(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( VXA(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( VXB(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( VXC(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( VYA(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( VYB(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( VYC(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )


      DO J=J_0,J_1S
C**** Calculate KH = rho_0 BETA* L_Munk^3 where DX=L_Munk
c      KHP(J)=2d0*RHOWS*OMEGA*COSP(J)*(DXP(J)**3)/RADIUS ! tracer lat
c      KHV(J)=2d0*RHOWS*OMEGA*COSV(J)*(DXV(J)**3)/RADIUS ! (v vel pts)
C**** Calculate KH=rho_0 BETA (sqrt(3) L_Munk/pi)^3, L_Munk=min(DX,DY)
        DSP=MIN(DXPO(J),DYPO(J))*2.*SQRT(3.)/TWOPI  ! tracer lat
        DSV=MIN(DXVO(J),DYVO(J))*2.*SQRT(3.)/TWOPI  ! v vel pts
        KHP(J)=2d0*RHOWS*OMEGA*COSPO(J)*(DSP**3)/RADIUS ! tracer lat
        KHV(J)=2d0*RHOWS*OMEGA*COSVO(J)*(DSV**3)/RADIUS ! (v vel pts)
        KHP(J)=MAX(KHP(J),AKHMIN)
        KHV(J)=MAX(KHV(J),AKHMIN)
        BYDXYV(J)=1D0/DXYVO(J)
        BYDXV(J)=1D0/DXVO(J)
        BYDXP(J)=1D0/DXPO(J)
        BYDYV(J)=1D0/DYVO(J)
        BYDYP(J)=1D0/DYPO(J)
        KYPXP(J)=KHP(J)*DYPO(J)*BYDXP(J)
        KXPYV(J)=KHV(J)*DXPO(J)*BYDYV(J)
        KYVXV(J)=KHV(J)*DYVO(J)*BYDXV(J)
        KXVYP(J)=KHP(J)*DXVO(J)*BYDYP(J)
C**** Discretisation errors need TANP/V to be defined like this
        DLAT = TWOPI*NINT(360d0/(JM-1))/720d0
        TANP(J)=TAN(RLAT(J))*TAN(0.5*DLAT)/(RADIUS*0.5*DLAT)
        VLAT = DLAT*(J+0.5-0.5*(1+JM))
        TANV(J)=TAN(VLAT)*SIN(DLAT)/(DLAT*RADIUS)
      END DO
      !make halo_update if 2 is not present
      !barrier synchoronization is required before sending the message,
      !j=2 is computed before the message is sent
      if( HAVE_SOUTH_POLE ) then
        KHV(1)=KHV(2)
        BYDXV(1)=1D0/DXVO(1)
        BYDYV(1)=1D0/DYVO(1)
        BYDYP(1)=1D0/DYPO(1)
      endif
      if( HAVE_NORTH_POLE ) then
        BYDXP(JM)=1D0/DXPO(JM)
        BYDXYPJM=1D0/(DXYPO(JM)*IM)
        TANP(JM) = 0
      endif

      CALL HALO_UPDATE(grid,TANP (grid%j_strt_halo:grid%j_stop_halo) ,
     *                 FROM=NORTH+SOUTH)
      CALL HALO_UPDATE(grid,BYDXP(grid%j_strt_halo:grid%j_stop_halo) ,
     *                 FROM=NORTH)
      CALL HALO_UPDATE(grid,TANV (grid%j_strt_halo:grid%j_stop_halo) ,
     *                 FROM=SOUTH)
      CALL HALO_UPDATE(grid,KHP (grid%j_strt_halo:grid%j_stop_halo) ,
     *                 FROM=NORTH)
      CALL HALO_UPDATE(grid,KXVYP (grid%j_strt_halo:grid%j_stop_halo) ,
     *                 FROM=SOUTH)
      CALL HALO_UPDATE(grid,DYPO (grid%j_strt_halo:grid%j_stop_halo) ,
     *                 FROM=SOUTH)

C****
C**** Calculate operators fixed in time for U and V equations
C****
      UXA=0. ; UXB=0. ; UXC=0. ; UYA=0. ; UYB=0. ; UYC=0.
      VXA=0. ; VXB=0. ; VXC=0. ; VYA=0. ; VYB=0. ; VYC=0.
      UYPB=0.; UYPA=0.

      DO L=1,LMO
C**** Calculate flux operators
C**** i.e DUDX(1) and DUDX(2) are the coefficients of u1 and u2 for
C**** calculating centered difference K_h du/dx
C**** including metric terms in y derivatives
        DUDX=0.
        DUDY=0.
        DVDX=0.
        DVDY=0.
        DO J=max(2,J_0S-1),J_1S
          I=IM
          DO IP1=1,IM
            IF (L.LE.LMU(IP1,J)) DUDX(IP1,J,1) = KYPXP(J)
            IF (L.LE.LMU(I,J)) THEN
              DUDX(IP1,J,2) = -KYPXP(J)
              IF (L.LE.LMU(I,J+1)) THEN
                DUDY(I,J,1) =  KXPYV(J)*(1. +0.5*TANV(J)*DYVO(J))
                DUDY(I,J,2) = -KXPYV(J)*(1. -0.5*TANV(J)*DYVO(J))
              ELSE
                DUDY(I,J,2) = -(1.-FSLIP)*2d0*KXPYV(J)
              END IF
            ELSE
              IF (L.LE.LMU(I,J+1)) DUDY(I,J,1) = (1.-FSLIP)*2d0*KXPYV(J)
            END IF
            IF (L.LE.LMV(I,J+1)) DVDY(I,J+1,1) = KXVYP(J)*(1. +
     *           0.5*TANP(J)*DYPO(J))
            IF (L.LE.LMV(I,J)) THEN
              DVDY(I,J+1,2) = -KXVYP(J)*(1.-0.5*TANP(J)*DYPO(J))
              IF (L.LE.LMV(IP1,J)) THEN
                DVDX(I,J,1) =  KYVXV(J)
                DVDX(I,J,2) = -KYVXV(J)
              ELSE
                DVDX(I,J,2) = -(1.-FSLIP)*2d0*KYVXV(J)
              END IF
            ELSE
              IF (L.LE.LMV(IP1,J)) DVDX(I,J,1) = (1.-FSLIP)*2d0*KYVXV(J)
            END IF
            I=IP1
          END DO
        END DO
        CALL HALO_UPDATE(grid, DUDY(:,grid%j_strt_halo:grid%j_stop_halo,
     *     :), FROM=SOUTH)
C****
C**** Combine to form tri-diagonal operators including first metric term
C****
        DO J=J_0S,J_1S
          IM1=IM
          DO I=1,IM
            IF(L.LE.LMU(IM1,J)) THEN
              UXA(IM1,J,L) = -DUDX(IM1,J,2)                  *BYDXYPO(J)
              UXB(IM1,J,L) = (DUDX(I  ,J,2) -DUDX(IM1,J  ,1))*BYDXYPO(J)
              UXC(IM1,J,L) =  DUDX(I  ,J,1)                  *BYDXYPO(J)
              UYA(IM1,J,L) = -DUDY(IM1,J-1,2)                *BYDXYPO(J)
     *                     + 0.5*TANP(J)*KHP(J)*BYDYP(J)
              UYB(IM1,J,L) =(DUDY(IM1,J  ,2)-DUDY(IM1,J-1,1))*BYDXYPO(J)
     *                     + TANP(J)*TANP(J)*KHP(J)
              UYC(IM1,J,L) = DUDY(IM1,J  ,1)                 *BYDXYPO(J)
     *                     - 0.5*TANP(J)*KHP(J)*BYDYP(J)
            END IF
            IF (L.LE.LMV(I,J)) THEN
              VXA(I,J,L) = -DVDX(IM1,J,2)                 *BYDXYV(J)
              VXB(I,J,L) = (DVDX(I  ,J,2) - DVDX(IM1,J,1))*BYDXYV(J)
              VXC(I,J,L) =  DVDX(I  ,J,1)                 *BYDXYV(J)
              VYA(I,J,L) = -DVDY(I,J  ,2)                 *BYDXYV(J)
     *                   + 0.5*TANV(J)*KHV(J)*BYDYV(J)
              VYB(I,J,L) = (DVDY(I,J+1,2) - DVDY(I  ,J,1))*BYDXYV(J)
     *                   + TANV(J)*TANV(J)*KHV(J)
              VYC(I,J,L) =  DVDY(I,J+1,1)                 *BYDXYV(J)
     *                   - 0.5*TANV(J)*KHV(J)*BYDYV(J)
            END IF
            IM1=I
          END DO
        END DO
C**** At North Pole
        !following uses jm and jm-1 data, if both are not on same
        !processor (need to check this), halo_date is required.
        IF(HAVE_NORTH_POLE) THEN
          IF(L.LE.LMU(1,JM)) THEN
            DO I=1,IM
              UYPB(L) = UYPB(L) - DUDY(I,JM-1,1)
              UYPA(I,L) = -DUDY(I,JM-1,2)*BYDXYPJM
            END DO
            UYPB(L) = UYPB(L)*BYDXYPJM
          END IF
        END IF
      END DO
      END IF
!**  IFIRST block ends here

C**** End of initialization from first call to ODIFF
C****
C**** Solve diffusion equations semi implicitly
C****
      DT2=DTDIFF*5d-1     ! half time-step
C**** Store North Pole velocity components, they will not be changed
      if(have_north_pole) then
        UONP(:) = UO(IM  ,JM,:)
        VONP(:) = UO(IVNP,JM,:)
      end if

!     also may later need HALO for DXPO(N), DXVO(S)

      CALL HALO_UPDATE(grid,DH,
     *                 FROM=NORTH)
      CALL HALO_UPDATE(grid,MO,
     *                 FROM=NORTH)

      CALL HALO_UPDATE(grid,UO (:,grid%j_strt_halo:grid%j_stop_halo,:) ,
     *                 FROM=NORTH)
      CALL HALO_UPDATE(grid,UO (:,grid%j_strt_halo:grid%j_stop_halo,:) ,
     *                 FROM=SOUTH)
      CALL HALO_UPDATE(grid,VO (:,grid%j_strt_halo:grid%j_stop_halo,:) ,
     *                 FROM=NORTH)
      CALL HALO_UPDATE(grid,VO (:,grid%j_strt_halo:grid%j_stop_halo,:) ,
     *                 FROM=SOUTH)

C****
!$OMP PARALLEL DEFAULT(NONE),
!$OMP&  PRIVATE(AU,AV, BYMU,BYMV,BU,BV, CU,CV, DTU, DTV,
!$OMP&          FUX,FUY,FVX,FVY, I,IP1,IM1, J, L, RU,RV,
!$OMP&          UU,UV,UT,UY,UX, VT,VY,VX),
!$OMP&  SHARED(have_north_pole, UO, VO, UONP, VONP, COSU,SINU,COSI,SINI,
!$OMP&         J_0, J_0S, J_1S, LMU, MO, LMV, TANP, TANV,
!$OMP&         BYDYP, BYDXV, BYDXP, BYDYV, KHV, KHP, grid, DT2, DH,
!$OMP&         UXA, UXB, UXC, UYA, UYB, UYC, VXA, VXB, VXC,
!$OMP&         VYA, VYB, VYC, DYVO, DXVO, DXPO, DYPO, BYDXYPO, BYDXYV)

      allocate( AU(IM,grid%j_strt_halo:grid%j_stop_halo) )
      allocate( BU(IM,grid%j_strt_halo:grid%j_stop_halo) )
      allocate( CU(IM,grid%j_strt_halo:grid%j_stop_halo) )
      allocate( RU(IM,grid%j_strt_halo:grid%j_stop_halo) )
      allocate( UU(IM,grid%j_strt_halo:grid%j_stop_halo) )

      allocate( AV(IM,grid%j_strt_halo:grid%j_stop_halo) )
      allocate( BV(IM,grid%j_strt_halo:grid%j_stop_halo) )
      allocate( CV(IM,grid%j_strt_halo:grid%j_stop_halo) )
      allocate( RV(IM,grid%j_strt_halo:grid%j_stop_halo) )
      allocate( UV(IM,grid%j_strt_halo:grid%j_stop_halo) )

!$OMP DO
      DO L=1,LMO
C**** Calculate rotating polar velocities from UONP and VONP
      if(have_north_pole) then
        UO(:,JM,L) = UONP(L)*COSU(:) + VONP(L)*SINU(:)
        VO(:,JM,L) = VONP(L)*COSI(:) - UONP(L)*SINI(:)
      end if
C**** Save (0.5*) mass reciprical for velocity points
      DO J=J_0S,J_1S
        I=IM
        DO IP1=1,IM
          IF (L.LE.LMU(I,J)) BYMU(I,J) = 1./(MO(I,J,L)+MO(IP1,J,L))
          IF (L.LE.LMV(I,J)) BYMV(I,J) = 1./(MO(I,J,L)+MO(I,J+1,L))
          I=IP1
        END DO
      END DO
      if( HAVE_NORTH_POLE ) then
        IF (L.LE.LMU(1,JM)) BYMU(1,JM) = 1./MO(1,JM,L)
      endif
C**** Calculate Wasjowicz boundary terms
C**** Need dv/dy,tv,dv/dx for u equation, du/dy,tu,du/dx for v equation
      FUX=0             ! flux in U equation at the x_+ boundary
      FUY=0             ! flux in U equation at the y_+ boundary
      FVX=0             ! flux in V equation at the x_+ boundary
      FVY=0             ! flux in V equation at the y_+ boundary
      DO J=J_0, J_1S
        IM1=IM-1
        DO I=1,IM
          UT=0          ! mean u*tan on x_+ boundary for V equation
          UY=0          ! mean du/dx on y_+ boundary for V equation
          UX=0          ! mean du/dy on x_+ boundary for V equation
          IF (L.LE.LMU(I  ,J+1)) THEN
            UT=     UO(I  ,J+1,L)*TANP(J+1)
            UX=     UO(I  ,J+1,L)
            UY=     UO(I  ,J+1,L)
          END IF
          IF (L.LE.LMU(I  ,J  )) THEN
            UT=UT + UO(I  ,J  ,L)*TANP(J  )
            UY=UY - UO(I  ,J  ,L)
          END IF
          IF (L.LE.LMU(IM1,J+1)) UX=UX-    UO(IM1,J+1,L)
          UT=0.5*UT
          UX=UX*BYDXP(J+1)
          UY=UY*BYDYV(J)
C****
          VT=0          ! mean v*tan on x_+ boundary for U equation
          VX=0          ! mean dv/dx on y_+ boundary for U equation
          VY=0          ! mean dv/dy on x_+ boundary for U equation
          IF (L.LE.LMV(I  ,J  )) THEN
            VT=     VO(I  ,J  ,L)*TANV(J  )
            VX=     VO(I  ,J  ,L)
            VY=     VO(I  ,J  ,L)
          END IF
          IF (J.GT.1) THEN
            IF  (L.LE.LMV(I  ,J-1)) THEN
            VT=VT + VO(I  ,J-1,L)*TANV(J-1)
            VY=VY - VO(I  ,J-1,L)
            END IF
          END IF
          IF (L.LE.LMV(IM1,J  )) VX=VX - VO(IM1,J  ,L)
          VT=0.5*VT
          VY=VY*BYDYP(J)
          VX=VX*BYDXV(J)
C**** Calculate fluxes (including FSLIP condition)
          IF (FSLIP == 1.) THEN
            IF (L.LE.LMV(I,J) .AND. L.LE.LMV(IM1,J))
     *           FUY(IM1,J)=KHV(J)*VX
            IF (L.LE.LMU(I,J) .AND. L.LE.LMU(I,J+1))
     *           FVX(I  ,J)=KHV(J)*(UY + UT)
          ELSE
            FUY(IM1,J)=KHV(J)*VX
            FVX(I  ,J)=KHV(J)*(UY + UT)
          END IF
          FUX(IM1,J)=KHP(J)*(VY + VT)
          IF (J.LT.JM-1) FVY(I,J)=KHP(J+1)*UX
          IM1=I
        END DO
      END DO

      CALL HALO_UPDATE(grid,FUY (:,grid%j_strt_halo:grid%j_stop_halo) ,
     *                 FROM=SOUTH)
      CALL HALO_UPDATE(grid,FVY (:,grid%j_strt_halo:grid%j_stop_halo),
     *                 FROM=SOUTH)

C**** Calculate tridiagonal matrix for first semi-implicit step (in x)
C**** Minor complication due to cyclic nature of boundary condition
      AU=0. ; BU=0. ; CU=0. ; RU=0.
      AV=0. ; BV=0. ; CV=0. ; RV=0.

      DO J=J_0S,J_1S
        IM1=IM-1
        I=IM
        DO IP1=1,IM
          BU(I,J) = 1d0
          BV(I,J) = 1d0
          IF (L.LE.LMU(I,J)) THEN
            DTU = DT2*(DH(I,J,L)+DH(IP1,J,L))*BYMU(I,J)
            IF (I.gt.1 ) AU(I,J) =        - DTU*UXA(I,J,L)
                         BU(I,J) = BU(I,J) - DTU*UXB(I,J,L)
            IF (I.lt.IM) CU(I,J) =        - DTU*UXC(I,J,L)
            RU(I,J) = UO(I,J,L) + DTU*(UYA(I,J,L)*UO(I,J-1,L)
     *           +UYB(I,J,L)*UO(I,J,L) + UYC(I,J,L)*UO(I,J+1,L))
C**** Make properly tridiagonal by making explicit cyclic terms
            IF (I == 1 ) RU(I,J)=RU(I,J) + DTU*UXA(I,J,L)*UO(IM,J,L)
            IF (I == IM) RU(I,J)=RU(I,J) + DTU*UXC(I,J,L)*UO(1,J,L)
C**** Add Wasjowicz cross-terms to RU + second metric term
            RU(I,J) = RU(I,J) + DTU*((DYPO(J)*(FUX(IM1,J) - FUX(I,J))
     *           + DXVO(J)*FUY(I,J) - DXVO(J-1)*FUY(I,J-1))*BYDXYPO(J)
     *           - 0.5*(TANV(J-1)*FUY(I,J-1) + TANV(J)*FUY(I,J)))
          END IF
          IF (L.LE.LMV(I,J)) THEN
            DTV = DT2*(DH(I,J,L)+DH(I,J+1,L))*BYMV(I,J)
            IF (I.gt.1 ) AV(I,J) =        - DTV*VXA(I,J,L)
                         BV(I,J) = BV(I,J) - DTV*VXB(I,J,L)
            IF (I.lt.IM) CV(I,J) =        - DTV*VXC(I,J,L)
            RV(I,J) = VO(I,J,L) + DTV*(VYA(I,J,L)*VO(I,J-1,L)
     *           +VYB(I,J,L)*VO(I,J,L) + VYC(I,J,L)*VO(I,J+1,L))
C**** Make properly tridiagonal by making explicit cyclic terms
            IF (I == 1 ) RV(I,J)=RV(I,J) + DTV*VXA(I,J,L)*VO(IM,J,L)
            IF (I == IM) RV(I,J)=RV(I,J) + DTV*VXC(I,J,L)*VO(1,J,L)
C**** Add Wasjowicz cross-terms to RV + second metric term
            RV(I,J) = RV(I,J) + DTV*((DYVO(J)*(FVX(I,J) - FVX(IM1,J))
     *           + DXPO(J)*FVY(I,J-1) - DXPO(J+1)*FVY(I,J))*BYDXYV(J)
     *           + 0.5*(TANP(J-1)*FVY(I,J-1) + TANP(J)*FVY(I,J)))
          END IF
          IM1=I
          I=IP1
        END DO
      END DO
C**** At North Pole (no metric terms)
c     BU(IIP) = 1d0
c     BV(IIP) = 1d0
c     IF (L.LE.LMU(1,JM)) THEN
c     DTU = DT2*DH(1,JM,L)*BYMU(1,JM)
c       RU(IIP) = 0.
c       DO I=1,IM       ! include Wasjowicz cross-terms at North Pole
c         RU(IIP) = RU(IIP) + DTU*(UYPA(I,L)*UO(I,JM-1,L)
c    *                      - DXVO(JM-1)*FUY(I,JM-1)*BYDXYPJM)
c       END DO
c     END IF
C**** Call tridiagonal solver
      DO J = J_0S, J_1S
        CALL TRIDIAG(AU(:,J), BU(:,J), CU(:,J), RU(:,J), UO(:,J,L), IM)
        CALL TRIDIAG(AV(:,J), BV(:,J), CV(:,J), RV(:,J), VO(:,J,L), IM)
      END DO

C****
C**** Now do semi implicit solution in y
C**** Recalculate rotating polar velocities from UONP and VONP
C     UO(:,JM,L) = UONP(L)*COSU(:) + VONP(L)*SINU(:)
C     VO(:,JM,L) = VONP(L)*COSI(:) - UONP(L)*SINI(:)

      CALL HALO_UPDATE(grid,UO (:,grid%j_strt_halo:grid%j_stop_halo,L) ,
     *                 FROM=NORTH)
      CALL HALO_UPDATE(grid,UO (:,grid%j_strt_halo:grid%j_stop_halo,L) ,
     *                 FROM=SOUTH)
      CALL HALO_UPDATE(grid,VO (:,grid%j_strt_halo:grid%j_stop_halo,L) ,
     *                 FROM=NORTH)
      CALL HALO_UPDATE(grid,VO (:,grid%j_strt_halo:grid%j_stop_halo,L) ,
     *                 FROM=SOUTH)

C**** Calc. cross-term fluxes + second metric term (at half time step)
C**** Need dv/dy,tv,dv/dx for u equation, du/dy,tu,du/dx for v equation
      FUX=0             ! flux in U equation at the x_+ boundary
      FUY=0             ! flux in U equation at the y_+ boundary
      FVX=0             ! flux in V equation at the x_+ boundary
      FVY=0             ! flux in V equation at the y_+ boundary
      DO J=J_0,J_1S
        IM1=IM-1
        DO I=1,IM
          UT=0         ! mean u*tan on x_+ boundary for V equation
          UY=0         ! mean du/dx on y_+ boundary for V equation
          UX=0         ! mean du/dy on x_+ boundary for V equation
          IF (L.LE.LMU(I  ,J+1)) THEN
            UT=     UO(I  ,J+1,L)*TANP(J+1)
            UX=     UO(I  ,J+1,L)
            UY=     UO(I  ,J+1,L)
          END IF
          IF (L.LE.LMU(I  ,J  )) THEN
            UT=UT + UO(I  ,J  ,L)*TANP(J  )
            UY=UY - UO(I  ,J  ,L)
          END IF
          IF (L.LE.LMU(IM1,J+1)) UX=UX-    UO(IM1,J+1,L)
          UT=0.5*UT
          UX=UX*BYDXP(J+1)
          UY=UY*BYDYV(J)
C****
          VT=0         ! mean v*tan on x_+ boundary for U equation
          VX=0         ! mean dv/dx on y_+ boundary for U equation
          VY=0         ! mean dv/dy on x_+ boundary for U equation
          IF (L.LE.LMV(I  ,J  )) THEN
            VT=     VO(I  ,J  ,L)*TANV(J  )
            VX=     VO(I  ,J  ,L)
            VY=     VO(I  ,J  ,L)
          END IF
          IF (J.GT.1) THEN
            IF (L.LE.LMV(I  ,J-1)) THEN
            VT=VT + VO(I  ,J-1,L)*TANV(J-1)
            VY=VY - VO(I  ,J-1,L)
          END IF
          END IF
          IF (L.LE.LMV(IM1,J  )) VX=VX - VO(IM1,J  ,L)
          VT=0.5*VT
          VY=VY*BYDYP(J)
          VX=VX*BYDXV(J)
C**** Calculate fluxes (including FSLIP condition)
          IF (FSLIP == 1.) THEN
            IF (L.LE.LMV(I,J) .AND. L.LE.LMV(IM1,J))
     *           FUY(IM1,J)=KHV(J)* VX
            IF (L.LE.LMU(I,J) .AND. L.LE.LMU(I,J+1))
     *           FVX(I,J)=KHV(J)*(UY + UT)
          ELSE
            FUY(IM1,J)=KHV(J)* VX
            FVX(I  ,J)=KHV(J)*(UY + UT)
          END IF
          FUX(IM1,J)=KHP(J)*(VY + VT)
          IF (J.LT.JM-1) FVY(I,J)=KHP(J+1)*UX
          IM1=I
        END DO
      END DO

      CALL HALO_UPDATE(grid,FUY (:,grid%j_strt_halo:grid%j_stop_halo) ,
     *                 FROM=SOUTH)
      CALL HALO_UPDATE(grid,FVY (:,grid%j_strt_halo:grid%j_stop_halo) ,
     *                 FROM=SOUTH)

C**** Calculate tridiagonal matrix for second semi-implicit step (in y)
C**** Minor complication due to singular nature of polar box
      AU=0. ; BU=0. ; CU=0. ; RU=0.; UU=0
      AV=0. ; BV=0. ; CV=0. ; RV=0.; UV=0.

      DO IP1=1,IM
        DO J=J_0S,J_1S
          !put following later into a subroutine
          if(ip1.eq.1) then
            if(J.eq.2) then
              IM1=IM-1; I=IM;
            elseif(J.eq.3) then
              IM1=IM; I=IP1;
            else
              IM1=IP1; I=IP1;
            endif
          endif
          if(ip1.gt.1) then
            if(j.eq.2) then
              IM1=IP1-1; I=IP1-1;
            else
              IM1=IP1; I=IP1;
            endif
          endif

          BU(I,J) = 1d0
          BV(I,J) = 1d0
          IF (L.LE.LMU(I,J)) THEN
            DTU = DT2*(DH(I,J,L)+DH(IP1,J,L))*BYMU(I,J)
            AU(I,J) =        - DTU*UYA(I,J,L)
            BU(I,J) = BU(I,J) - DTU*UYB(I,J,L)
            IF (J.lt.JM-1) CU(I,J) =     - DTU*UYC(I,J,L)
            RU(I,J) = UO(I,J,L) + DTU*(UXA(I,J,L)*UO(IM1,J,L)
     *           +UXB(I,J,L)*UO(I,J,L) + UXC(I,J,L)*UO(IP1,J,L))
c**** Make properly tridiagonal by making explicit polar terms
!mkt  UO(1,JM,L) changed to UO(I,JM,L)
            IF (J == JM-1) RU(I,J)=RU(I,J) + DTU*UYC(I,J,L)*UO(I,JM,L)
C**** Add Wasjowicz cross-terms to RU + second metric term
            RU(I,J) = RU(I,J) + DTU*((DYPO(J)*(FUX(IM1,J) - FUX(I,J))
     *           + DXVO(J)*FUY(I,J) - DXVO(J-1)*FUY(I,J-1))*BYDXYPO(J)
     *           - 0.5*(TANV(J-1)*FUY(I,J-1) + TANV(J)*FUY(I,J)))
          END IF
          IF (L.LE.LMV(I,J)) THEN
            DTV = DT2*(DH(I,J,L)+DH(I,J+1,L))*BYMV(I,J)
            AV(I,J) =        - DTV*VYA(I,J,L)
            BV(I,J) = BV(I,J) - DTV*VYB(I,J,L)
            IF (J.lt.JM-1) CV(I,J) =     - DTV*VYC(I,J,L)
            RV(I,J) = VO(I,J,L) + DTV*(VXA(I,J,L)*VO(IM1,J,L)
     *           +VXB(I,J,L)*VO(I,J,L) + VXC(I,J,L)*VO(IP1,J,L))
c**** Make properly tridiagonal by making explicit polar terms
            IF (J == JM-1) RV(I,J)=RV(I,J) + DTV*VYC(I,J,L)*VO(I,JM,L)
C**** Add Wasjowicz cross-terms to RV + second metric term
            RV(I,J) = RV(I,J) + DTV*((DYVO(J)*(FVX(I,J) - FVX(IM1,J))
     *           + DXPO(J)*FVY(I,J-1) - DXPO(J+1)*FVY(I,J))*BYDXYV(J)
     *           + 0.5*(TANP(J-1)*FVY(I,J-1) + TANP(J)*FVY(I,J)))
          END IF
          IM1=I
          I=IP1
        END DO
      END DO
C**** At North Pole (do partly explicitly) no metric terms
c     BU(IIP) = 1d0
c     BV(IIP) = 1d0
c     IF (L.LE.LMU(1,JM)) THEN
c       DTU = DT2*DH(1,JM,L)*BYMU(1,JM)
c       BU(IIP) = BU(IIP) - DTU*UYPB(L)
c       RU(IIP) = UO(1,JM,L)
c       DO I=1,IM       ! include Wasjowicz cross-terms at North Pole
c         RU(IIP)= RU(IIP) + DTU*(UYPA(I,L)*UO(I,JM-1,L)
c    *         - DXVO(JM-1)*FUY(I,JM-1)*BYDXYPJM)
c       END DO
c     END IF
C**** Call tridiagonal solver

      CALL TRIDIAG_new(AU, BU, CU, RU, UU, grid,J_LOWER=2,J_UPPER=JM-1)
      CALL TRIDIAG_new(AV, BV, CV, RV, UV, grid,J_LOWER=2,J_UPPER=JM-1)

      UO(:,J_0S:J_1S,L)=UU(:,J_0S:J_1S)
      VO(:,J_0S:J_1S,L)=UV(:,J_0S:J_1S)

C****
      END DO
!$OMP END DO
      deallocate(AU, BU, CU, RU, UU)
      deallocate(AV, BV, CV, RV, UV)
!$OMP END PARALLEL

C**** Restore unchanged UONP and VONP into prognostic locations in UO
      if(have_north_pole) then
        UO(IM  ,JM,:) = UONP(:)
        UO(IVNP,JM,:) = VONP(:)
      end if
C****
      RETURN
      END SUBROUTINE ODIFF


      SUBROUTINE TOC2SST
!@sum  TOC2SST convert ocean surface variables into atmospheric sst
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : tf
#if (defined TRACERS_OCEAN) || (defined TRACERS_WATER)
      USE OCN_TRACER_COM, only : trw0, ntm
#endif
      USE OCEAN, only : im,jm,IVSP,IVNP, imaxj,dxypo,lmm,
     *     COSU,SINU, COSI=>COSIC,SINI=>SINIC,
     *     focean,g0m,s0m,mo,ogeoz,uo,vo,ogeoz_sv
#ifdef TRACERS_OCEAN
     *     ,trmo
#endif
      USE FLUXES, only : gtemp, sss, mlhc, ogeoza, uosurf, vosurf,
     *      gtempr
#ifdef TRACERS_ON
     *     ,gtracer
#endif
#ifdef TRACERS_GASEXCH_Natassa
     *     ,trgasex
#endif

      use domain_decomp, only : grid, get

      IMPLICIT NONE
      INTEGER I,J
      REAL*8 TEMGS,shcgs,GO,SO,GO2,SO2,TO
      integer :: j_0,j_1,n
      logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      call get (grid, j_strt=j_0, j_stop=j_1,
     * HAVE_SOUTH_POLE=HAVE_SOUTH_POLE, HAVE_NORTH_POLE=HAVE_NORTH_POLE)

C****
C**** Note that currently everything is on same grid
C****
      DO J=J_0,J_1
        DO I=1,IMAXJ(J)
          IF (FOCEAN(I,J).gt.0.) THEN
            GO= G0M(I,J,1)/(MO(I,J,1)*DXYPO(J))
            SO= S0M(I,J,1)/(MO(I,J,1)*DXYPO(J))
            TO= TEMGS(GO,SO)
            GTEMP(1,1,I,J)= TO
            GTEMPR(1,I,J) = TO+TF
            SSS(I,J) = 1d3*SO
            MLHC(I,J)= MO(I,J,1)*SHCGS(GO,SO)
            IF (LMM(I,J).gt.1) THEN
              GO2= G0M(I,J,2)/(MO(I,J,2)*DXYPO(J))
              SO2= S0M(I,J,2)/(MO(I,J,2)*DXYPO(J))
              TO= TEMGS(GO2,SO2)
            END IF
            GTEMP(2,1,I,J)= TO
   ! atmospheric grid Ocean height
            OGEOZA(I,J)=0.5*(OGEOZ(I,J)+OGEOZ_SV(I,J))
            UOSURF(I,J)=UO(I,J,1)
            VOSURF(I,J)=VO(I,J,1)

#ifdef TRACERS_WATER
#ifdef TRACERS_OCEAN
            GTRACER(:,1,I,J)=TRMO(I,J,1,:)/(MO(I,J,1)*DXYPO(J)-
     *           S0M(I,J,1))
#else
            GTRACER(:,1,I,J)=trw0(:)
#endif
#endif

#ifdef TRACERS_GASEXCH_Natassa
            GTRACER(:,1,I,J)=TRMO(I,J,1,:)/(MO(I,J,1)*DXYPO(J))
#endif
          END IF

        END DO
      END DO
C**** do poles
      if (HAVE_NORTH_POLE) then
      IF (FOCEAN(1,JM).gt.0) THEN
        UOSURF(1,JM) = UO(IM,JM,1)*COSU(1) + UO(IVNP,JM,1)*SINU(1)
        VOSURF(1,JM) = UO(IVNP,JM,1)*COSI(1) - UO(IM,JM,1)*SINI(1)
        DO I=2,IM
          GTEMP(:,1,I,JM)=GTEMP(:,1,1,JM)
          GTEMPR(1,I,JM) =GTEMPR(1,1,JM) 
          SSS(I,JM)=SSS(1,JM)
          MLHC(I,JM)=MLHC(1,JM)
          UOSURF(I,JM) = UO(IM,JM,1)*COSU(I) + UO(IVNP,JM,1)*SINU(I)
          VOSURF(I,JM) = UO(IVNP,JM,1)*COSI(I) - UO(IM,JM,1)*SINI(I)
          OGEOZA(I,JM)=OGEOZA(1,JM)
#if (defined TRACERS_WATER) || (defined TRACERS_GASEXCH_Natassa)
          GTRACER(:,1,I,JM)=GTRACER(:,1,1,JM)
#endif
        END DO
      END IF
      end if

      if (HAVE_SOUTH_POLE) then
      IF (FOCEAN(1,1).gt.0) THEN
        DO I=2,IM
          GTEMP(:,1,I,1)=GTEMP(:,1,1,1)
          GTEMPR(1,I,1) =GTEMPR(1,1,1)
          SSS(I,1)=SSS(1,1)
          MLHC(I,1)=MLHC(1,1)
          UOSURF(I,1) = UO(IM,1,1)*COSU(1) - UO(IVSP,1,1)*SINU(1)
          VOSURF(I,0) = UO(IVSP,1,1)*COSI(I) - UO(IM,1,1)*SINI(I)
          OGEOZA(I,1)=OGEOZA(1,1)
#if (defined TRACERS_WATER) || (defined TRACERS_GASEXCH_Natassa)
          GTRACER(:,1,I,1)=GTRACER(:,1,1,1)
#endif
        END DO
      END IF
      end if

      RETURN
C****
      END SUBROUTINE TOC2SST

      SUBROUTINE io_oda(kunit,it,iaction,ioerr)
!@sum  io_oda dummy routine for consistency with uncoupled model
!@auth Gavin Schmidt
!@ver  1.0
      RETURN
      END SUBROUTINE io_oda

      SUBROUTINE ADVSI_DIAG
!@sum ADVSI_DIAG dummy routine for consistency with qflux model
      RETURN
      END SUBROUTINE ADVSI_DIAG

      SUBROUTINE AT2OT(FIELDA,FIELDO,NF,QCONSERV)
!@sum  AT2OT interpolates Atm Tracer grid to Ocean Tracer grid
!@auth Gavin Schmidt
!@ver  1.0
c GISS-ESMF EXCEPTIONAL CASE - AT2OT not needed yet - nothing done yet
      USE MODEL_COM, only : ima=>im,jma=>jm
      USE OCEAN, only : imo=>im,jmo=>jm,imaxj,ratoc
      use domain_decomp, only : grid,get

      IMPLICIT NONE
!@var QCONSERV true if integrated field must be conserved
      LOGICAL, INTENT(IN) :: QCONSERV
!@var N number of fields
      INTEGER, INTENT(IN) :: NF
!@var FIELDA array on atmospheric tracer grid
      REAL*8, INTENT(IN), DIMENSION(NF,IMA,JMA) :: FIELDA
!@var FIELDO array on oceanic tracer grid
      REAL*8, INTENT(OUT), DIMENSION(NF,IMO,JMO) :: FIELDO
      INTEGER I,J

C**** currently no need for interpolation,
C**** just scaling due to area differences for fluxes
      IF (QCONSERV) THEN
        DO J=1,JMO
          DO I=1,IMAXJ(J)
            FIELDO(:,I,J) = FIELDA(:,I,J)*RATOC(J)
          END DO
        END DO
        DO I=2,IMO
          FIELDO(:,I,JMO)=FIELDO(:,1,JMO)
          FIELDO(:,I,1  )=FIELDO(:,1,  1)
        END DO
      ELSE
        FIELDO = FIELDA
      END IF
C****
      RETURN
      END SUBROUTINE AT2OT

      SUBROUTINE OT2AT(FIELDO,FIELDA,NF,QCONSERV)
c GISS-ESMF EXCEPTIONAL CASE - OT2AT not needed yet - nothing done yet
!@sum  OT2AT interpolates Ocean Tracer grid to Atm Tracer grid
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ima=>im,jma=>jm
      USE GEOM, only : imaxj
      USE OCEAN, only : imo=>im,jmo=>jm,rocat
      IMPLICIT NONE
!@var QCONSERV true if integrated field must be conserved
      LOGICAL, INTENT(IN) :: QCONSERV
!@var N number of fields
      INTEGER, INTENT(IN) :: NF
!@var FIELDA array on atmospheric tracer grid
      REAL*8, INTENT(OUT), DIMENSION(NF,IMA,JMA) :: FIELDA
!@var FIELDO array on oceanic tracer grid
      REAL*8, INTENT(IN), DIMENSION(NF,IMO,JMO) :: FIELDO
      INTEGER I,J

C**** currently no need for interpolation,
C**** just scaling due to area differences for fluxes
      IF (QCONSERV) THEN
        DO J=1,JMA
          DO I=1,IMAXJ(J)
            FIELDA(:,I,J) = FIELDO(:,I,J)*ROCAT(J)
          END DO
        END DO
        DO I=2,IMA
          FIELDA(:,I,JMA)=FIELDA(:,1,JMA)
          FIELDA(:,I,1  )=FIELDA(:,1,  1)
        END DO
      ELSE
        FIELDA = FIELDO
      END IF
C****
      RETURN
      END SUBROUTINE OT2AT

      SUBROUTINE AT2OV(FIELDA,FIELDO,NF,QCONSERV,QU)
!@sum  AT2OV interpolates Atm Tracer grid to Ocean Velocity grid
!@auth Gavin Schmidt
!@ver  1.0
c GISS-ESMF EXCEPTIONAL CASE - AT2OV not needed yet - nothing done yet
      USE MODEL_COM, only : ima=>im,jma=>jm
      USE GEOM, only : imaxj
      USE OCEAN, only : imo=>im,jmo=>jm,ramvn,ramvs,ratoc
      IMPLICIT NONE
!@var QCONSERV true if integrated field must be conserved
!@var QU true if u-velocity pts. are wanted (false for v velocity pts.)
      LOGICAL, INTENT(IN) :: QCONSERV, QU
!@var N number of fields
      INTEGER, INTENT(IN) :: NF
!@var FIELDA array on atmospheric tracer grid
      REAL*8, INTENT(IN), DIMENSION(NF,IMA,JMA) :: FIELDA
!@var FIELDO array on oceanic tracer grid
      REAL*8, INTENT(OUT), DIMENSION(NF,IMO,JMO) :: FIELDO
      INTEGER I,J,IP1

C**** Ocean velocities are on C grid
      IF (QU) THEN  ! interpolate onto U points
        IF (QCONSERV) THEN
          DO J=1,JMO-1
            I=IMO
            DO IP1=1,IMO
              FIELDO(:,I,J)=0.5*RATOC(J)*(FIELDA(:,I,J)+FIELDA(:,IP1,J))
              I=IP1
            END DO
          END DO
C**** do poles
          DO I=1,IMO
            FIELDO(:,I,JMO) = FIELDA(:,1,JMO)*RATOC(JMO)
            FIELDO(:,I,  1) = FIELDA(:,1,  1)*RATOC(JMO)
          END DO
        ELSE   ! no area weighting
          DO J=1,JMO-1
            I=IMO
            DO IP1=1,IMO
              FIELDO(:,I,J) = 0.5*(FIELDA(:,I,J)+FIELDA(:,IP1,J))
              I=IP1
            END DO
          END DO
C**** do poles
          DO I=1,IMO
            FIELDO(:,I,JMO) = FIELDO(:,1,JMO)
            FIELDO(:,I,  1) = FIELDO(:,1,  1)
          END DO
        END IF
      ELSE                      ! interpolate onto V points
        IF (QCONSERV) THEN
          DO J=1,JMO-1
            DO I=1,IMO
              FIELDO(:,I,J) = RATOC(J)*RAMVN(J)*FIELDA(:,I,J)+
     *             RATOC(J+1)*RAMVS(J+1)*FIELDA(:,I,J+1)
            END DO
          END DO
        ELSE
          DO J=1,JMO-1
            DO I=1,IMO
              FIELDO(:,I,J) = 0.5*(FIELDA(:,I,J)+FIELDA(:,I,J+1))
            END DO
          END DO
        END IF
        FIELDO(:,:,JMO) = 0.
      END IF
C****
      RETURN
      END SUBROUTINE AT2OV

      SUBROUTINE GLMELT(DT)
!@sum  GLMELT adds glacial melt around Greenland and Antarctica to ocean
!@auth Sukeshi Sheth/Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : dtsrc
      USE OCEAN, only : imaxj,jm,g0m,mo,ze,focean,dxypo,lmm
#ifdef TRACERS_OCEAN
     *     ,trmo
#endif
      USE FLUXES, only : gmelt,egmelt
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
     *     ,trgmelt
#endif
#endif               /* TNL: inserted */
      USE OCEANRES, only : maxgl
      use domain_decomp, only : grid,get

      IMPLICIT NONE
      REAL*8, INTENT(IN) :: DT  !@var DT timestep for GLMELT call
      REAL*8 DZ
      INTEGER I,J,L

      integer :: j_0,j_1

      call get (grid, J_STRT=j_0, J_STOP=j_1)

      DO L=1,MAXGL
C**** divide over depth and scale for time step
        DO J=j_0,j_1
          DO I=1,IMAXJ(J)
            IF (FOCEAN(I,J).GT.0. .and. GMELT(I,J).gt.0) THEN
              DZ=DT*(ZE(L)-ZE(L-1))/(DTsrc*ZE(MIN(MAXGL,LMM(I,J))))
              MO(I,J,L) = MO(I,J,L)+GMELT(I,J)*DZ/(DXYPO(J)*FOCEAN(I,J))
              G0M(I,J,L)=G0M(I,J,L)+EGMELT(I,J)*DZ
#ifdef TRACERS_WATER  /* TNL: inserted */
#ifdef TRACERS_OCEAN
              TRMO(I,J,L,:)=TRMO(I,J,L,:)+TRGMELT(:,I,J)*DZ
#endif
#endif               /* TNL: inserted */
            END IF
          END DO
        END DO
      END DO
C****
      RETURN
      END SUBROUTINE GLMELT


      SUBROUTINE ADJUST_MEAN_SALT
!@sum  ADJUST_MEAN_SALT sets the global mean salinity in the ocean
!@auth Gavin Schmidt
!@ver  1.0
      USE GEOM, only : dxyp,imaxj
      USE CONSTANT, only : grav
      USE OCEAN, only : im,jm,lmo,focean,lmm,mo,s0m,sxmo
     *     ,symo,szmo,dxypo,oc_salt_mean,g0m
      USE STRAITS, only : s0mst,sxmst,szmst,nmst,lmst,g0mst,mmst,dist
     *     ,wist
      use DOMAIN_DECOMP, only: grid, GLOBALSUM, get, AM_I_ROOT,
     *     ESMF_BCAST
      IMPLICIT NONE
      REAL*8 :: totalSalt,totalMass

      REAL*8, DIMENSION(grid%J_STRT_HALO:grid%J_STOP_HALO) ::OSALTJ
     *     ,OMASSJ
      REAL*8 mean_S,frac_inc,T_ORIG,T_NEW,temgsp,shcgs,pres,g,s
      INTEGER I,J,L,N,J_0,J_1

      CALL GET(grid, J_STRT = J_0, J_STOP = J_1)

      call conserv_OSL(OSALTJ)
      call conserv_OMS(OMASSJ)
      OSALTJ(:)=OSALTJ(:)*DXYP(:)
      OMASSJ(:)=OMASSJ(:)*DXYP(:)

      CALL GLOBALSUM(grid, OSALTJ, totalSalt, ALL=.true.)
      CALL GLOBALSUM(grid, OMASSJ, totalMass, ALL=.true.)

      if (AM_I_ROOT()) then
        mean_S=1000*totalSalt/totalMass  ! psu
        frac_inc=oc_salt_mean/mean_S
        write(6,*) "Changing ocean salinity: ",mean_S,frac_inc
      end if
      call ESMF_BCAST(grid,  frac_inc)

C**** adjust open ocean salinity
      DO J=J_0,J_1
        DO I=1,IMAXJ(J)
          PRES=0
          DO L=1,LMM(I,J)
            PRES=PRES+MO(I,J,L)*GRAV*.5
            G=G0M(I,J,L)/(MO(I,J,L)*DXYPO(J))
            S=S0M(I,J,L)/(MO(I,J,L)*DXYPO(J))
            T_ORIG=TEMGSP (G,S,PRES)
            S0M(I,J,L)=frac_inc*S0M(I,J,L)
            if (frac_inc.lt.1.) then
              SXMO(I,J,L)=frac_inc*SXMO(I,J,L)
              SYMO(I,J,L)=frac_inc*SYMO(I,J,L)
              SZMO(I,J,L)=frac_inc*SZMO(I,J,L)
            end if
            S=S0M(I,J,L)/(MO(I,J,L)*DXYPO(J))
            T_NEW=TEMGSP (G,S,PRES)
C**** approximately adjust enthalpy to restore temperature
            G0M(I,J,L)=G0M(I,J,L)+(T_ORIG-T_NEW)*SHCGS(G,S)*MO(I,J,L)
     *           *DXYPO(J)
            G=G0M(I,J,L)/(MO(I,J,L)*DXYPO(J))
            IF (L.lt.LMM(I,J)) PRES=PRES+MO(I,J,L)*GRAV*.5
          END DO
        END DO
      END DO

C**** adjust strait salinity
      if (am_I_root()) then
        DO N=1,NMST
          PRES=0
          DO L=1,LMST(N)
            PRES=PRES+MMST(L,N)*GRAV*.5/(DIST(N)*WIST(N))
            G=G0MST(L,N)/MMST(L,N)
            S=S0MST(L,N)/MMST(L,N)
            T_ORIG=TEMGSP (G,S,PRES)
            S0MST(L,N)=frac_inc*S0MST(L,N)
            if (frac_inc.lt.1) then
              SXMST(L,N)=frac_inc*SXMST(L,N)
              SZMST(L,N)=frac_inc*SZMST(L,N)
            end if
            S=S0MST(L,N)/MMST(L,N)
            T_NEW=TEMGSP (G,S,PRES)
C**** approximately adjust enthalpy to restore temperature
            G0MST(L,N)=G0MST(L,N)+(T_ORIG-T_NEW)*SHCGS(G,S)*MMST(L,N)
            G=G0MST(L,N)/MMST(L,N)
            IF (L.lt.LMST(N)) PRES=PRES+MMST(L,N)*GRAV*.5/(DIST(N)
     *           *WIST(N))
          END DO
        END DO
      end if
      CALL ESMF_BCAST(grid, S0MST)
      CALL ESMF_BCAST(grid, SXMST)
      CALL ESMF_BCAST(grid, SZMST)

C**** Check
      call conserv_OSL(OSALTJ)
      OSALTJ(:)=OSALTJ(:)*DXYP(:)
      CALL GLOBALSUM(grid, OSALTJ, totalSalt, ALL=.true.)

      if (AM_I_ROOT()) then
        mean_S=1000*totalSalt/totalMass  ! psu
        write(6,*) "New ocean salinity: ",mean_S,oc_salt_mean
      end if

      RETURN
      END SUBROUTINE ADJUST_MEAN_SALT
