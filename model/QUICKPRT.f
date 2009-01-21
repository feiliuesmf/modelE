#include "rundeck_opts.h"

      MODULE WORKJK
      USE MODEL_COM, ONLY : JM,LM
!!!
!Replaces:
!!!      COMMON/WORKJK/DPJK,DPHEM,DPGLOB

      REAL*8, DIMENSION(JM,LM,2):: DPJK
      REAL*8, DIMENSION(2,LM,2) :: DPHEM
      REAL*8, DIMENSION(LM,2) :: DPGLOB

      END MODULE WORKJK


!------------------------------------------------

      BLOCK DATA BDWP
C****
C**** TITLES FOR SUBROUTINE DIAG7
C****
      COMMON/D7COM/LNAME,SNAME,UNITS
      CHARACTER LNAME(12)*50,SNAME(12)*30,UNITS(12)*50
      DATA LNAME/
     1'WAVE POWER FOR U NEAR 850 MB AND EQUATOR  ',
     2'WAVE POWER FOR V NEAR 850 MB AND EQUATOR  ',
     3'WAVE POWER FOR U NEAR 300 MB AND EQUATOR  ',
     4'WAVE POWER FOR V NEAR 300 MB AND EQUATOR  ',
     5'WAVE POWER FOR U NEAR 50 MB AND EQUATOR   ',
     6'WAVE POWER FOR V NEAR 50 MB AND EQUATOR   ',
     7'WAVE POWER FOR PHI AT 922 MB AND 50 DEG N.',
     8'WAVE POWER FOR PHI AT 700 MB AND 50 DEG N.',
     9'WAVE POWER FOR PHI AT 500 MB AND 50 DEG N.',
     A'WAVE POWER FOR PHI AT 300 MB AND 50 DEG N.',
     B'WAVE POWER FOR PHI AT 100 MB AND 50 DEG N.',
     C'WAVE POWER FOR PHI AT 10 MB AND 50 DEG N. '/
!      .........1.........2.........3.........4.........5.........6
      DATA UNITS/
     1'DAY*(m/s)^2   ','DAY*(m/s)^2   ','10 DAY*(m/s)^2',
     4'DAY*(m/s)^2   ','10 DAY*(m/s)^2','DAY*(m/s)^2   ',
     7'10**3 DAY*m^2 ','10**3 DAY*m^2 ','10**3 DAY*m^2 ',
     A'10**3 DAY*m^2 ','10**4 DAY*m^2 ','10**4 DAY*m^2 '/
      DATA SNAME/
     1'WPU850EQU'   ,'WPV850EQU'   ,'WPU300EQU'   ,'WPV300EQU'   ,
     5'WPU50EQU'    ,'WPV50EQU'    ,'WPPHI922_50N','WPPHI700_50N',
     9'WPPHI500_50N','WPPHI300_50N','WPPHI100_50N','WPPHI10_50N' /

      END BLOCK DATA BDWP

!------------------------------------------------

      MODULE BDJ
!@sum  stores information for outputting composite zonal diagnostics
!@auth M. Kelley
      IMPLICIT NONE
      SAVE
!@param nj_out number of j-format output fields = 11
      integer, parameter :: nj_out=11
!@var units string containing output field units
      CHARACTER(LEN=50), DIMENSION(nj_out) :: UNITS_J_O
!@var lname string describing output field
      CHARACTER(LEN=50), DIMENSION(nj_out) :: LNAME_J_O
!@var sname string referencing output field in self-desc. output file
      CHARACTER(LEN=30), DIMENSION(nj_out) :: NAME_J_O
!@var stitle short title for print out
      CHARACTER(LEN=16), DIMENSION(nj_out) :: STITLE_J_O
!@var INUM_J_O,IDEN_J_O numerator and denominator for calculated J diags
      INTEGER, DIMENSION(nj_out) :: INUM_J_O, IDEN_J_O
!@var SCALE_J_O scale for calculated J diags
      REAL*8, DIMENSION(nj_out) :: SCALE_J_O

      END MODULE BDJ

!------------------------------------------------

      MODULE BDjkjl
!@sum  stores information for outputting lat-sigma/pressure diagnostics
!@auth M. Kelley
      IMPLICIT NONE
      SAVE
!@param names of derived jk/jl output fields
      INTEGER :: jl_rad_cool,jk_dudt_econv,jl_nt_lh_e,jl_vt_lh_e,
     *  jk_psi_cp,jk_dudt_epdiv,jk_stdev_dp,
     *  jk_dtempdt_econv,jl_phi_amp_wave1,jl_phi_phase_wave1,
     *  jl_epflx_div,jk_vt_dse_e,jk_vt_lh_eddy,jk_vt_se_eddy,
     *  jk_tot_vt_se,jk_psi_tem,jk_epflx_v,
     *  jk_nt_eqgpv,jk_dyn_conv_eddy_geop,jk_nt_sheat_e,
     *  jk_dyn_conv_dse,jk_seke,jk_eke,
     *  jk_nt_dse_se,jk_nt_dse_e,jk_tot_nt_dse,
     *  jk_nt_lh_e,jk_nt_see,jk_tot_nt_se,
     *  jk_nt_am_stand_eddy,jk_nt_am_eddy,jk_tot_nt_am,
     *  jk_we_flx_nor,jk_we_flx_div,jk_refr_ind_wave1,
     *  jk_del_qgpv,jk_nt_lh_se,jk_wstar,jk_vstar,
     *  jl_mcdrgpm10,jl_mcdrgpm40,jl_mcdrgpm20,jl_sumdrg

      END MODULE BDjkjl

!------------------------------------------------

      MODULE BDIJ
!@sum  stores information for outputting lon-lat diagnostics
!@auth M. Kelley
      use MODEL_COM, only : IM,JM
      use DIAG_COM
      IMPLICIT NONE
      SAVE

!@param nij_o total number of diagnostic ij-format fields
      integer nij_o

!@var ij_xxx non single-aij diagnostic names
      INTEGER :: ij_topo, ij_jet, ij_wsmn, ij_jetdir, ij_wsdir, ij_grow,
     *  ij_netrdp, ij_albp, ij_albg, ij_albv, ij_ntdsese, ij_ntdsete,
     *  ij_fland, ij_dzt1, ij_albgv, ij_colh2o, ij_msu2,ij_msu3,ij_msu4,
     *  ij_Tatm, ij_RTSE, ij_HWV, ij_PVS

!@var SENTDSE stand.eddy northw. transport of dry static energy * 16
!@var TENTDSE trans.eddy northw. transport of dry static energy * 16
!@var TMSU2-4 MSU channel 2-4 temperatures (C)
      REAL*8, DIMENSION(IM,JM) :: SENTDSE,TENTDSE, TMSU2,TMSU3,TMSU4

      contains

      function mark (val,ibar,undef)
!@sum  mark selects a character (color) based on value and color bar
!@auth R. Ruedy
!@ver  1.0
      real*8 val,undef
      integer ibar,n
      character*1 mark

      if (val .eq. undef) then
        mark=' '
      else
      select case (ibar)
      case (ib_pct)                                ! 0.....100 %
        n = 2.5 + val
        if (val .ge. 20.) n=23
        if (val .le.  0.) n= 1
        mark = cbar(ib_pct)(n:n)
      case (ib_pos)                                ! 0++++++++++
        n = 2.5 + val
c          non-unif scaling: (currently not used)
c          if (n .gt. 13) n = (n+123)/10
        if (n .gt. 38) n=38
        if (n .lt. 1 .or. val .le. 0.) n= 1
        mark = cbar(ib_pos)(n:n)
      case (ib_npp)                                ! ---0+++++++
        n = 11.5 + val
        if (n .gt. 38.) n=38
        if (n .lt.  1 ) n= 1
        mark = cbar(ib_npp)(n:n)
      case (ib_nnp)                                ! -------0+++
        n = 28.5 + val
        if (n .gt. 38.) n=38
        if (n .lt.  1 ) n= 1
        mark = cbar(ib_nnp)(n:n)
      case (ib_hyb)                      ! hybrid: multiple scales
        n = 2.5 + val
        if (n .gt. 28) n=(n+263)/10
        if (n .gt. 35) n=(n+180)/6
        if (n .gt. 37) n=37
        if (val .le.  0.) n=1
        mark = cbar(ib_hyb)(n:n)
      case (ib_ntr)                !tracers       ! ---0+++++++
        if (val.lt.0.) then
          n = 11.5-LOG(-val)/LOG(2.)
          if (n .le.  0) n= 1
          if (n .gt. 11) n=11
        else if (val.eq.0.) then
          n = 11
        else
          n = 11.5+LOG( val)/LOG(2.)
          if (n .lt. 11) n=11
          if (n .ge. 38) n=38
        end if
        mark = cbar(ib_npp)(n:n)                  ! use ib_npp
      end select
      end if

      return
      end function mark

      function ib_of_legnd (leg)
!@sum  ib_of_legnd finds the 'colorbar' for the given legend
!@auth R. Ruedy
!@ver  1.0
      integer ib_of_legnd, leg

      ib_of_legnd = ib_pos
      if (legend(leg)(7:8) .eq. ',Z') ib_of_legnd = ib_nnp
      if (legend(leg)(7:8) .eq. ',9') ib_of_legnd = ib_npp
      if (index(legend(leg)(21:40),'-') .gt. 0) ib_of_legnd = ib_hyb
      if (index(legend(leg),'100 ') .gt. 0) ib_of_legnd = ib_pct
      if (legend(leg)(1:4) .eq. '9=-5') ib_of_legnd = ib_ntr

      return
      end function ib_of_legnd

      END MODULE BDIJ

!--------------------------------------------------------


!====================
!      SUBROUTINE DIAGKS
      module DIAGKS
C****
C**** THIS ROUTINE PRODUCES A SUMMARY OF KEY NUMBERS CALCULATED IN
C**** OTHER DIAGNOSTIC ROUTINES
C****
C**** CONTENTS OF KEYNR
C****
C****  N
C****
C****  1 MONTH
C****  2 TOTAL CLOUD COVER (PERCENT)
C****  3 SNOW AND ICE COVERAGE OF GLOBE (PERCENT)
C****  4 SNOW AND ICE COVERAGE OF NORTHERN HEMISPHERE (PERCENT)
C****  5 SNOW COVER--NORTHERN HEMSIPHERE (PERCENT)
C****  6 ICE COVER--NORTHERN HEMISPHERE (PERCENT)
C****  7 PLANETARY ALBEDO (PERCENT)
C****  8 SOLAR RADIATION ABSORBED BY ATMOSPHERE (WT/M**2)
C****  9 SOLAR RADIATION ABSORBED BY PLANET (WT/M**2)
C**** 10 NET HEAT AT GROUND (WT/M**2)
Cobso    ANGULAR MOMENTUM PER UNIT AREA (10**10 J*SEC/M**2)
Cobso    EVAPORATION (.1 MM/DAY)
C**** 11 PRECIPITATION (.1 MM/DAY)
C**** 12 SENSIBLE HEAT FLUX INTO GROUND (ABS.VALUE)
C**** 13 LATENT HEAT FLUX INTO GROUND (ABS.VALUE)
C**** 14 MEAN GROUND TEMPERATURE (DEGREES K)
C**** 15 MEAN GLOBAL ATMOSPHERIC TEMPERATURE (DEGREES K)
C**** 16 MERID. TEMPERATURE GRADIENT (N.HEMISPHERE)
C**** 17 MERID. TEMPERATURE GRADIENT (S.HEMISPHERE)
C**** 18 MEAN TROPOSPHERIC EKE-NORTHERN HEMISPHERE
C**** 19 MEAN TROPOSPHERIC EKE-SOUTHERN HEMISPHERE
C**** 20 MEAN TROPOSPHERIC ZKE-NORTHERN HEMISPHERE
C**** 21 MEAN TROPOSPHERIC ZKE-SOUTHERN HEMISPHERE
C**** 22 MEAN TROPOSPHERIC EPE-NORTHERN HEMISPHERE
C**** 23 MEAN TROPOSPHERIC ZPE-NORTHERN HEMISPHERE
C**** 24 MEAN EDDY KINETIC ENERGY AT EQUATOR
C**** 25 MAX. MEAN EDDY KINETIC ENERGY IN MID NORTH LATITUDES
C**** 26 MAX. ZONAL WIND (U COMPONENT) IN TROPOSPHERE (NH), M/SEC
C**** 27 LATITUDE CORRESPONDING TO 26
C**** 28 MAX. ZONAL WIND (U COMPONENT) IN TROPOSPHERE (SH), M/SEC
C**** 29 LATITUDE CORRESPONDING TO 28
C**** 30-32: 31 IS LARGEST VALUE OF STREAM FUNCTION, POSITIVE OR
C****    NEGATIVE; 30 AND 32 ARE THE MAGNITUDES OF THE LARGEST VALUES OF
C****    OPPOSITE SIGN TO THE NORTH AND SOUTH RESPECTIVELY
C**** 33-42 REFER TO NORTHERN HEMISPHERE ONLY
C**** 33 MAX.NORTHWARD TRANS. OF DRY STATIC ENERGY BY STANDING EDDIES
C**** 34 MAX.NORTHWARD TRANS. OF DRY STATIC ENERGY BY EDDIES
C**** 35 MAX. TOTAL NORTH. TRANS. OF DRY STATIC ENERGY
C**** 36 MAX.NORTHWARD TRANS. OF STATIC ENERGY BY EDDIES
C**** 37 MAX.TOTAL NORTH. TRANS. OF STATIC ENERGY
C**** 38 LATITUDE CORRESPONDING TO 37
C**** 39 MAX. NORTH. TRANS. OF ANGULAR MOMENTUM BY STANDING EDDIES
C**** 40 MAX. NORTH. TRANS. OF ANGULAR MOMENTUM BY EDDIES
C**** 41 MAX. TOTAL NORTH. TRANS. OF ANGULAR MOMENTUM
C**** 42 LATITUDE CORRESPONDING TO 41
C****
      USE CONSTANT, only : twopi
      USE MODEL_COM, only : jm,lm,jeq, JHOUR,JHOUR0,
     &     JDATE,JDATE0,JMON,JMON0,AMON,AMON0,JYEAR,JYEAR0,
     &     Itime,ItimeI,Itime0,XLABEL,AMONTH,nday,pmidl00
      USE DIAG_COM, only : keyct,keynr,ned,nkeynr
      USE PARAM
      IMPLICIT NONE
      PRIVATE
      SAVE

      public KEYDJ,KEYJKT,KEYJKJ,KEYJLS,KEYJKE,KEYJKN,KEYIJ
     &     ,KEYD4,DIAGKN

      !REAL*8, DIMENSION(JM) :: FLAT
      REAL*8, DIMENSION(JM,LM) :: FKEY
      !REAL*8, DIMENSION(JM) :: ASUM
      !REAL*8, DIMENSION(2) :: HSUM
      !INTEGER, DIMENSION(2*NED) :: IK

      INTEGER ::
     &     I,I35,I70,IEND,ISIGN,
     &     J,J60,JMAX,JNDEX,JSTART,
     &     K,KEYMAX,KNDEX,
     &     LL,LMAX,LNLM,LNM,LSLM,LSM

      REAL*8 ::
     &     A,BIG,CEPT,CHECK,
     &     HN,HS,
     &     SAVE,DAYS,TEQ,TNOR,TSOU,
     &     UNLM,UNM,USLM,USM,X60

      contains
C****
C**** ENTRIES CALLED FROM DIAGJ
C****
!      ENTRY KEYDJ (N,FGLOB,FNH)
      subroutine KEYDJ(name,FGLOB,FNH)
      character*20 name
      real*8 FGLOB,FNH

      SELECT CASE ( name )
      CASE ('J_totcld')           ; KEYNR( 2,keyct) = NINT(FGLOB)
      CASE ('J_snow_cover')       ; KEYNR( 5,keyct) = NINT(FNH)
      CASE ('J_ocn_lak_ice_frac') ; KEYNR( 6,keyct) = NINT(FNH)
      CASE ('plan_alb')           ; KEYNR( 7,keyct) = NINT(10.*fglob)
      CASE ('J_sw_abs_atm')       ; KEYNR( 8,keyct) = NINT(fglob)
      CASE ('J_net_rad_p0')       ; KEYNR( 9,keyct) = NINT(fglob)
      CASE ('J_nt_ht_z0')         ; KEYNR(10,keyct) = NINT(fglob)
      CASE ('J_prec')             ; KEYNR(11,keyct) = NINT(10.*fglob)
      CASE ('J_snsht_flx')        ; KEYNR(12,keyct) = NINT(-fglob)
      CASE ('J_evht_flx')         ; KEYNR(13,keyct) = NINT(-fglob)
      CASE ('J_tg1')              ; KEYNR(14,keyct) = NINT(.1*fglob)
!!!   CASE ('J_tair')             ; KEYNR(15,keyct) = NINT(.1*fglob)
      end select
      RETURN
      end subroutine KEYDJ
C****
C**** ENTRIES CALLED FROM DIAGJL VIA JLMAP OR FROM DIAGJK VIA JKMAP
C****
!      ENTRY KEYJKT (GSUM,ASUM)
      subroutine KEYJKT (GSUM,ASUM)
      real*8 GSUM
      REAL*8, DIMENSION(JM) :: ASUM
c     routine emptied - add content here
      RETURN
      end subroutine KEYJKT
C****
!      ENTRY KEYJKJ (L,FLAT)

      subroutine KEYJKJ (L,FLAT)
      integer L
      REAL*8, DIMENSION(JM) :: FLAT
c     routine emptied - add content here
      RETURN
      end subroutine KEYJKJ
C****
!      ENTRY KEYJLS (L,FLAT)
      subroutine KEYJLS (L,FLAT)
      integer L
      REAL*8, DIMENSION(JM) :: FLAT
C**** STREAM FUNCTION
      DO 290 J=2,JM
  290 FKEY(J,L)=FLAT(J)
      IF (L.NE.1) RETURN
  300 SAVE=0.
      HS=0.
      HN=0.
      DO 310 K=1,LM
      DO 310 I=2,JM
      CHECK=ABS(FKEY(I,K))
      IF (CHECK.LT.SAVE) GO TO 310
      SAVE=CHECK
      JNDEX=I
      KNDEX=K
  310 CONTINUE
      SAVE=FKEY(JNDEX,KNDEX)
      ISIGN=1
      IF (SAVE.GT.0.0) ISIGN=-1
      IF (JNDEX.LT.4) GO TO 325
      IEND=JNDEX-1
      DO 320 K=1,LM
      DO 320 I=2,IEND
      CHECK=FKEY(I,K)*ISIGN
  320 IF (CHECK.GT.HS)HS=CHECK
  325 CONTINUE
      IF (JNDEX.GT.(JM-2))GO TO 335
      JSTART=JNDEX+1
      DO 330 K=1,LM
      DO 330 I=JSTART,JM
      CHECK=FKEY(I,K)*ISIGN
  330 IF (CHECK.GT.HN)HN=CHECK
  335 CONTINUE
      KEYNR(30,KEYCT)=ABS(HN)+0.5
      KEYNR(31,KEYCT)=NINT(SAVE)
      KEYNR(32,KEYCT)=ABS(HS)+0.5
      RETURN
      end subroutine KEYJLS
C****
!      ENTRY KEYJKE (NT,HSUM,ASUM)
      subroutine KEYJKE (NT,HSUM,ASUM)
      integer NT
      REAL*8, DIMENSION(2) :: HSUM
      REAL*8, DIMENSION(JM) :: ASUM
C**** EDDY AND ZONAL KINETIC ENERGY
      IF (NT.EQ.19) GO TO 450
      KEYNR(18,KEYCT)=NINT(HSUM(2))
      KEYNR(19,KEYCT)=NINT(HSUM(1))
      KEYNR(20,KEYCT)=KEYNR(20,KEYCT)-NINT(HSUM(2))
      KEYNR(21,KEYCT)=KEYNR(21,KEYCT)-NINT(HSUM(1))
      KEYNR(24,KEYCT)=NINT(ASUM(JEQ))

      BIG=-99999.
      I35=2.+(JM-1.)*125./180.
      I70=2.+(JM-1.)*160./180.
      DO 440 I=I35,I70
      IF (ASUM(I).LT.BIG) GO TO 440
      BIG=ASUM(I)
  440 CONTINUE
      KEYNR(25,KEYCT)=NINT(BIG)
      RETURN
  450 KEYNR(20,KEYCT)=KEYNR(20,KEYCT)+NINT(HSUM(2))
      KEYNR(21,KEYCT)=KEYNR(21,KEYCT)+NINT(HSUM(1))
      RETURN
      end subroutine KEYJKE
C****
!      ENTRY KEYJKN (NT,ASUM,SUMFAC)
      subroutine KEYJKN (NT,ASUM,SUMFAC)
      integer NT
      REAL*8, DIMENSION(JM) :: ASUM
      REAL*8 SUMFAC
c     routine emptied - add content here
      RETURN
      end subroutine KEYJKN
C****
C**** ENTRY CALLED FROM DIAGIJ
C****
!      ENTRY KEYIJ(PISG,PISN)
      subroutine KEYIJ(PISG,PISN)
      REAL*8 PISG,PISN
      KEYNR(3,KEYCT)=NINT(PISG)
      KEYNR(4,KEYCT)=NINT(PISN)
      RETURN
      end subroutine KEYIJ
C****
C**** ENTRY CALLED FROM DIAG4
C****
!      ENTRY KEYD4 (IK)
      subroutine KEYD4 (IK)
      INTEGER, DIMENSION(2*NED) :: IK
      KEYNR(22,KEYCT)=(IK(10)+IK(20)+5)/10
      KEYNR(23,KEYCT)=(IK(8)+IK(18)+5)/10
      RETURN
      end subroutine KEYD4
C****
!      ENTRY DIAGKN
      subroutine DIAGKN
C**** PRINTS THE TABLE OF KEY NUMBERS
C****
      DAYS=(Itime-Itime0)/FLOAT(nday)
      KEYNR(1,KEYCT)=JMON0
      IF (Itime.eq.ItimeI+1) KEYNR(1,KEYCT)=0
      IF (KEYCT.GE.2) THEN
        if (KEYNR(1,KEYCT-1).EQ.JMON0) KEYCT=KEYCT-1
      ENDIF
      WRITE(6,901) XLABEL
      WRITE(6,910) JYEAR0,AMON0,JDATE0,JHOUR0,
     *  JYEAR,AMON,JDATE,JHOUR,ITIME,DAYS
      WRITE(6,902)
      DO 810 K=1,KEYCT
      IF (KEYNR(1,K).EQ.1) WRITE (6,905)
  810 WRITE(6,905) AMONTH(KEYNR(1,K)),(KEYNR(I,K),I=2,42)
      WRITE (6,915)
CB       DO 815 I=1,42
CB815    FKEYDS(I)=KEYNR(I,KEYCT)
      KEYCT=KEYCT+1
      KEYMAX=49
      IF (KEYNR(1,1).NE.0) KEYMAX=48
      IF (KEYCT.LE.KEYMAX) RETURN
C**** ROLL UP KEY NUMBERS 1 YEAR AT A TIME
      DO 820 K=1,36
      DO 820 I=1,NKEYNR
  820 KEYNR(I,K)=KEYNR(I,K+KEYMAX-36)
      DO 880 K=37,50
      DO 880 I=1,NKEYNR
  880 KEYNR(I,K)=0
      KEYCT=37
      RETURN
  901 FORMAT('1',A)
  902 FORMAT ('0',7X,'SN+IC NH NH AL AB NT NT PR        T   T-OF-ATM  EK
     *E   ZKE           EKE   JET-STREAMS STREAM-FN NOR-TRAN NOR-TRAN NO
     *RTH-TRANS'/
     *         5X,'CL GL    SN OI BE BY RD HT EC SN LAT OF  GL  GRAD ---
     *-- ----- EPE ZPE ------ NORTH SOUTH --------- DRY-STAT STAT-ENR ZO
     *N MOMENTM'/
     *      5X,'CV OB NH CV CV DO AT P0 Z0 IP HT  HT GD  OB NH SH NH ','
     *SH NH SH  NH  NH EQ  ML VL LT VL LT NH MAX SH SE ED TL ED TL LT SE
     * ED TL LT'/)
  905 FORMAT (1X,A3,4I3,I2,I4,5I3,I4,I3,I4,6I3,2I4,I3,I4,5I3,I4,11I3)
  910 FORMAT ('0',15X,'From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,
     *  '  Model-Time:',I9,5X,'Dif:',F7.2,' Days')
  915 FORMAT('0')
      end subroutine DIAGKN

      END MODULE DIAGKS


!==================



      MODULE DIAG_SERIAL
      USE MODEL_COM, ONLY : IM, JM
      USE DOMAIN_DECOMP_ATM, only : grid, DIST_GRID, AM_I_ROOT
      USE DIAGKS

      PRIVATE
      PUBLIC :: PRINT_DIAGS
      PUBLIC :: JLMAP
      PUBLIC :: MAPTXT
      PUBLIC :: IJ_avg

!ESMF: These globalsum routines are private to this module and execute
!      serially in a single processor.
      INTERFACE GLOBALSUM
        MODULE PROCEDURE GLOBALSUM_J
        MODULE PROCEDURE GLOBALSUM_JK
      END INTERFACE

      !REAL*8 :: FLAND_glob(IM,JM)
      !REAL*8 :: FEARTH_glob(IM,JM)
      REAL*8 :: FOCEAN_glob(IM,JM)
      REAL*8 :: FLICE_glob(IM,JM)
      REAL*8 :: ZATMO_glob(IM,JM)

      CONTAINS




      SUBROUTINE GLOBALSUM_J(grd_dum, garr, gsum,
     &                       hsum, istag, iskip, all)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      REAL*8,            INTENT(IN) :: garr(grd_dum%jm_world)
      REAL*8,            INTENT(OUT):: gsum
      REAL*8, OPTIONAL,  INTENT(OUT):: hsum(2)
      INTEGER,OPTIONAL,  INTENT(IN) :: istag
      INTEGER,OPTIONAL,  INTENT(IN) :: iskip
      LOGICAL,OPTIONAL,  INTENT(IN) :: all

      INTEGER :: IM, JM, J, ier
      LOGICAL :: istag_, iskip_


      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD

      istag_ = .false.
      If (Present(istag)) Then
        If (istag == 1) istag_ = .true.
      End If

      iskip_ = .false.
      If (Present(iskip)) Then
        If (iskip == 1) iskip_ = .true.
      End If

      If (istag_) then
        gsum = sum(garr(2:JM),1)
      ElseIf (iskip_) then
        gsum = sum(garr(2:JM-1),1)
      Else
        gsum = sum(garr(1:JM),1)
      EndIf
      If (Present(hsum)) then
        If (istag_) then
          hsum(1)   = Sum( garr(2     :JM/2),1   )
          hsum(2)   = Sum( garr(2+JM/2:JM  ),1   )
          hsum(1)   = hsum(1) + 0.5*garr(1+JM/2)
          hsum(2)   = hsum(2) + 0.5*garr(1+JM/2)
        Else
          hsum(1)   = Sum( garr(1     :JM/2),1   )
          hsum(2)   = Sum( garr(1+JM/2:JM  ),1   )
        EndIf
      EndIf

      END SUBROUTINE GLOBALSUM_J


      SUBROUTINE GLOBALSUM_JK(grd_dum, garr, gsum, hsum, istag, all)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      REAL*8,            INTENT(IN) :: garr(              : ,:)
!     REAL*8,            INTENT(IN) :: garr(grd_dum%jm_world,:)
      REAL*8,            INTENT(OUT):: gsum(size(garr,2))
      REAL*8, OPTIONAL,  INTENT(OUT):: hsum(2,size(garr,2))
      INTEGER,OPTIONAL,  INTENT(IN) :: istag
      LOGICAL,OPTIONAL,  INTENT(IN) :: all

      INTEGER :: k
      INTEGER :: ier
      INTEGER :: IM, JM
      LOGICAL :: istag_

      IM   = grd_dum%IM_WORLD
      JM   = grd_dum%JM_WORLD
      istag_ = .false.
      If (Present(istag)) Then
        If (istag == 1) istag_ = .true.
      End If

      If (istag_) then
        gsum = sum(garr(2:JM,:),1)
      Else
        gsum = sum(garr(1:JM,:),1)
      EndIf
      If (Present(hsum)) then
        If (istag_) then
          hsum(1,:)   = Sum( garr(2     :JM/2,:),1   )
          hsum(2,:)   = Sum( garr(2+JM/2:JM  ,:),1   )
          hsum(1,:)   = hsum(1,:) + 0.5*garr(1+JM/2,:)
          hsum(2,:)   = hsum(2,:) + 0.5*garr(1+JM/2,:)
        Else
          hsum(1,:)   = Sum( garr(1     :JM/2,:),1   )
          hsum(2,:)   = Sum( garr(1+JM/2:JM  ,:),1   )
        EndIf
      EndIf

      END SUBROUTINE GLOBALSUM_JK


!------------------------------------------------


      subroutine print_diags(partial)
!@sum print_diag prints out binary and ascii diag output.
!@auth  Original Development Team
      USE MODEL_COM, only : itime,itimeI
      USE DIAG_COM, only : kdiag,keynr,keyct,isccp_diags
      IMPLICIT NONE
!@var partial : accum period is complete (if =0) or partial (if =1)
      INTEGER, INTENT(IN) :: partial

#ifdef SCM
c     write(0,*) 'SCM no diags   print_diags'
      return
#endif

      CALL DIAG_GATHER

      IF (AM_I_ROOT()) THEN

      call certain_wind_averages

      IF (KDIAG(1).LT.9) CALL DIAGJ
      IF (KDIAG(2).LT.9) CALL DIAGJK
      IF (KDIAG(10).LT.9) CALL DIAGIL
      IF (KDIAG(7).LT.9) CALL DIAG7P
      IF (KDIAG(3).LT.9) CALL DIAGIJ
      IF (KDIAG(9).LT.9) CALL DIAGCP
      IF (KDIAG(5).LT.9) CALL DIAG5P
      IF (partial.eq.0 .and. KDIAG(6).LT.9) CALL DIAGDD  ! full period
      IF (KDIAG(13).LT.9) CALL DIAGDH
      IF (KDIAG(4).LT.9) CALL DIAG4
      IF (KDIAG(11).LT.9) CALL diag_RIVER
      IF (KDIAG(12).LT.9) CALL diag_OCEAN
      IF (KDIAG(12).LT.9) CALL diag_ICEDYN
      IF (isccp_diags.eq.1) CALL diag_ISCCP
      IF (partial.eq.0 .or. Itime.LE.ItimeI+1) THEN  ! full period or IC
        CALL DIAGKN
      ELSE                      ! RESET THE UNUSED KEYNUMBERS TO ZERO
        KEYNR(1:42,KEYCT)=0
      END IF
#ifdef TRACERS_ON
      IF (KDIAG(8).LT.9) then
        CALL DIAGJLT
        CALL DIAGIJT
        CALL DIAGTCP
      end if
#endif
      END IF ! AM_I_ROOT

      CALL DIAG_SCATTER

      return
      end subroutine print_diags


      SUBROUTINE J_TITLES
!@sum  J_TITLES calculated zonal diagnostics
!@auth M. Kelley/G. Schmidt
!@ver  1.0
      USE DIAG_COM, only : j_srincp0,j_srnfp0,j_plavis,j_planir,j_srnfg
     *     ,j_srincg,j_albvis,j_albnir,j_srrvis,j_srrnir,j_sravis
     *     ,j_sranir,j_clddep,j_pcldmc
      USE BDJ
      IMPLICIT NONE
      INTEGER :: K
C**** These information are for J zonal budget calulated diagnostics
C**** Note that we assume that they are all ratios of two existing
C**** records.
c
      k = 0
c
      k = k + 1
      name_j_o(k) = 'plan_alb'
      lname_j_o(k) = ' TOTAL PLANETARY ALBEDO'
      units_j_o(k) = '%'
      stitle_j_o(k)= ' PLANETARY ALBDO'
      inum_j_o(k)  = J_SRNFP0
      iden_j_o(k)  = J_SRINCP0
      scale_j_o(k) = 100.
c
      k = k + 1
      name_j_o(k) = 'plan_alb_vis'
      lname_j_o(k) = 'PLANETARY ALBEDO IN VISUAL'
      units_j_o(k) = '%'
      stitle_j_o(k)= ' PLAN ALB VISUAL'
      inum_j_o(k)  = J_PLAVIS
      iden_j_o(k)  = J_SRINCP0
      scale_j_o(k) = 100.
c
      k = k + 1
      name_j_o(k) = 'plan_alb_nir'
      lname_j_o(k) = 'PLANETARY ALBEDO IN NEAR IR'
      units_j_o(k) = '%'
      stitle_j_o(k)= ' PLAN ALB NEARIR'
      inum_j_o(k)  = J_PLANIR
      iden_j_o(k)  = J_SRINCP0
      scale_j_o(k) = 100.
c
      k = k + 1
      name_j_o(k) = 'surf_alb'
      lname_j_o(k) = 'GROUND ALBEDO'
      units_j_o(k) = '%'
      stitle_j_o(k)= ' SURFACE G ALBDO'
      inum_j_o(k)  = J_SRNFG
      iden_j_o(k)  = J_SRINCG
      scale_j_o(k) = 100.
c
      k = k + 1
      name_j_o(k) = 'surf_alb_vis'
      lname_j_o(k) = 'GROUND ALBEDO IN VISUAL'
      units_j_o(k) = '%'
      stitle_j_o(k)= ' SURF ALB VISUAL'
      inum_j_o(k)  = J_ALBVIS
      iden_j_o(k)  = J_SRINCP0
      scale_j_o(k) = 100.
c
      k = k + 1
      name_j_o(k) = 'surf_alb_nir'
      lname_j_o(k) = 'GROUND ALBEDO IN NEAR IR'
      units_j_o(k) = '%'
      stitle_j_o(k)= ' SURF ALB NEARIR'
      inum_j_o(k)  = J_ALBNIR
      iden_j_o(k)  = J_SRINCP0
      scale_j_o(k) = 100.
c
      k = k + 1
      name_j_o(k) = 'atm_alb_vis'
      lname_j_o(k) = 'ATMOSPHERIC ALBEDO IN VISUAL'
      units_j_o(k) = '%'
      stitle_j_o(k)= '0ATMO ALB VISUAL'
      inum_j_o(k)  = J_SRRVIS
      iden_j_o(k)  = J_SRINCP0
      scale_j_o(k) = 100.
c
      k = k + 1
      name_j_o(k) = 'atm_alb_nir'
      lname_j_o(k) = 'ATMOSPHERIC ALBEDO IN NEAR IR'
      units_j_o(k) = '%'
      stitle_j_o(k)= ' ATMO ALB NEARIR'
      inum_j_o(k)  = J_SRRNIR
      iden_j_o(k)  = J_SRINCP0
      scale_j_o(k) = 100.
c
      k = k + 1
      name_j_o(k) = 'atm_abs_vis'
      lname_j_o(k) = 'ATMOSPHERIC ABSORPTION IN VISUAL'
      units_j_o(k) = 'W/m**2'
      stitle_j_o(k)= ' ATMO ABS VISUAL'
      inum_j_o(k)  = J_SRAVIS
      iden_j_o(k)  = J_SRINCP0
      scale_j_o(k) = 100.
c
      k = k + 1
      name_j_o(k) = 'atm_abs_nir'
      lname_j_o(k) = 'ATMOSPHERIC ABSORPTION IN NEAR IR'
      units_j_o(k) = 'W/m**2'
      stitle_j_o(k)= ' ATMO ABS NEARIR'
      inum_j_o(k)  = J_SRANIR
      iden_j_o(k)  = J_SRINCP0
      scale_j_o(k) = 100.
c
      k = k + 1
      name_j_o(k) = 'mc_clddp'
      lname_j_o(k) = 'MOIST CONVECTIVE CLOUD DEPTH'
      units_j_o(k) = 'mb'
      stitle_j_o(k)= ' MC CLD DPTH(MB)'
      inum_j_o(k)  = J_CLDDEP
      iden_j_o(k)  = J_PCLDMC
      scale_j_o(k) = 1.

      RETURN
      END SUBROUTINE J_TITLES

      SUBROUTINE DIAGJ
c     routine emptied - add content here
      END SUBROUTINE DIAGJ



      SUBROUTINE JKJL_TITLEX
!@sum  JKJL_TITLEX titles etc for composite jl, jk output
!@auth G. Schmidt/M. Kelley/J. Lerner
!@ver  1.0
      use filemanager
      USE CONSTANT, only : sday,bygrav,sha,lhe
      USE MODEL_COM, only : byim,DTsrc,fim
      USE BDjkjl
      USE DIAG_COM
      IMPLICIT NONE
      INTEGER :: K,kk,iu_Ijk
      LOGICAL qIjk,Ql(KAJLx),Qk(KAJKx)
      character*80 line
c
c derived JL-arrays
c
      k = kajl
c
      k = k + 1
      jl_rad_cool = k                         ; jgrid_jl(k) = 1
      sname_jl(k) = 'rad_cool'
      lname_jl(k) = 'TOTAL RADIATION COOLING RATE'
      units_jl(k) = 'W/(m^2*mb)'
      ia_jl(k) = ia_rad
      pow_jl(k) = -2
      k = k + 1
      jl_phi_amp_wave1 = k                    ; jgrid_jl(k) = 1
      sname_jl(k) = 'phi_amp_wave1'
      lname_jl(k) ='AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 1'
      units_jl(k) = 'METERS'
      k = k + 1
      sname_jl(k) = 'phi_amp_wave2'           ; jgrid_jl(k) = 1
      lname_jl(k) ='AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 2'
      units_jl(k) = 'METERS'
      k = k + 1
      sname_jl(k) = 'phi_amp_wave3'           ; jgrid_jl(k) = 1
      lname_jl(k) ='AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 3'
      units_jl(k) = 'METERS'
      k = k + 1
      sname_jl(k) = 'phi_amp_wave4'           ; jgrid_jl(k) = 1
      lname_jl(k) ='AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 4'
      units_jl(k) = 'METERS'
      k = k + 1
      jl_phi_phase_wave1 = k                  ; jgrid_jl(k) = 1
      sname_jl(k) = 'phi_phase_wave1'
      lname_jl(k) = 'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 1'
      units_jl(k) = 'DEG WEST LONG'
      k = k + 1
      sname_jl(k) = 'phi_phase_wave2'         ; jgrid_jl(k) = 1
      lname_jl(k) = 'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 2'
      units_jl(k) = 'DEG WEST LONG'
      k = k + 1
      sname_jl(k) = 'phi_phase_wave3'         ; jgrid_jl(k) = 1
      lname_jl(k) = 'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 3'
      units_jl(k) = 'DEG WEST LONG'
      k = k + 1
      sname_jl(k) = 'phi_phase_wave4'         ; jgrid_jl(k) = 1
      lname_jl(k) = 'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 4'
      units_jl(k) = 'DEG WEST LONG'
      k = k + 1
      jl_epflx_div = k                        ; jgrid_jl(k) = 1
      sname_jl(k) = 'epflx_div'
      lname_jl(k) = 'DIVERGENCE OF THE ELIASSEN-PALM FLUX'
      units_jl(k) = 'm/s^2'
      scale_jl(k) = 1.
      pow_jl(k) = -6
      k = k + 1
      jl_mcdrgpm10 = k                        ; jgrid_jl(k) = 2
      sname_jl(k) = 'dudt_mcdrgpm10'
      lname_jl(k) = 'DU/DT BY STRAT. MC DRAG  C=+/-10R'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      k = k + 1
      jl_mcdrgpm40 = k                        ; jgrid_jl(k) = 2
      sname_jl(k) = 'dudt_mcdrgpm40' !AJL24+25
      lname_jl(k) = 'DU/DT BY STRAT. MC DRAG  C=+/-40R'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      k = k + 1
      jl_mcdrgpm20 = k                        ; jgrid_jl(k) = 2
      sname_jl(k) = 'dudt_mcdrgpm20' !AJL26+27
      lname_jl(k) = 'DU/DT BY STRAT MC DRAG C=+/-20R'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      k = k + 1
      jl_sumdrg = k                           ; jgrid_jl(k) = 2
      sname_jl(k) = 'dudt_sumdrg' !AJL(18+20-27)
      lname_jl(k) = 'ZONAL WIND CHANGE BY MTN+DEFORM+SHR+MC DRAG'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      k = k + 1
      jl_nt_lh_e = k
      sname_jl(k) = 'nt_lh_eddy'        ; jgrid_jl(k) = 2
      lname_jl(k) = 'N. TRANSPORT OF LATENT HEAT BY EDDIES (QDYN)'
      units_jl(k) = 'W/mb'
      scale_jl(k) = 100.*bygrav*xwon*lhe*fim/DTsrc
      pow_jl(k) = 10
      ia_jl(k) = ia_src
      k = k + 1
      jl_vt_lh_e = k
      sname_jl(k) = 'vt_lh_eddy1'        ; jgrid_jl(k) = 1
      lname_jl(k) = 'V. TRANSPORT OF LATENT HEAT BY EDDIES'
      units_jl(k) = 'W/m^2'
      scale_jl(k) = 100.*bygrav*xwon*lhe*byim/DTsrc
      pow_jl(k) = 0
      ia_jl(k) = ia_src

c Check the count
      if (k .gt. KAJLx) then
        write (*,*) 'Increase KAJLx=',KAJLx,' to at least ',k
        call stop_model('JL_TITLES: KAJLx too small',255)
      end if

      inquire(file='Ijk',exist=qIjk)
      if(.not.qIjk) then
         call openunit('Ijk',iu_Ijk,.false.,.false.)
         write(iu_Ijk, FMT='(a)') 'list of JL-fields'
         do kk = 1,k
           write(iu_Ijk, '(i3,1x,a)') kk,lname_jl(kk)
         end do
      else if(kdiag(2).gt.0) then
         Ql=.false.
         call openunit('Ijk',iu_Ijk,.false.,.true.)
         read (iu_Ijk,'(a)',end=20) line
   10    read (iu_Ijk,'(a)',end=20) line
         if(line(1:1).eq.'l') go to 20
         read(line,'(i3)') kk
         Ql(kk)=.true.
         go to 10
   20    continue
         do kk=1,KAJLx
           if(.not.Ql(kk)) sname_jl(kk)='skip'
         end do
       end if
c
c derived JK-arrays
c
      k = KAJK
c
      k = k + 1
      jk_dudt_econv = k                       ; jgrid_jk(k) = jgrid_u
      sname_jk(k) = 'dudt_eddy_conv'
      lname_jk(k) = 'DU/DT BY EDDY CONVERGENCE (CP)'
      units_jk(k) = '10**-6 m/s^2'
      scale_jk(k)= 1.D6
      k = k + 1
      jk_psi_cp = k                           ; jgrid_jk(k) = 2
      sname_jk(k) = 'psi_cp'
      lname_jk(k) = 'STREAM FUNCTION (CP)'
      units_jk(k) = 'kg/s' !'10**9 KILOGRAMS/SECOND'
      scale_jk(k) = 100.*BYGRAV
      ia_jk(k) = ia_dga
      pow_jk(k) = 9
      k = k + 1
      jk_dudt_epdiv = k                       ; jgrid_jk(k) = jgrid_u
      sname_jk(k) = 'dudt_epdiv'
      lname_jk(k) = 'DU/DT BY ELIASSEN-PALM DIVERGENCE (CP)'
      units_jk(k) = 'm/s^2'
      pow_jk(k) = -6
      k = k + 1
      jk_stdev_dp = k                         ; jgrid_jk(k) = 2
      sname_jk(k) = 'stdev_dp'
      lname_jk(k) = 'STANDARD DEVIATION OF PRESSURE DIFFERENCES'
      units_jk(k) = 'MB'
      scale_jk(k) = 1.
      k = k + 1
      jk_dtempdt_econv = k                    ; jgrid_jk(k) = 1
      sname_jk(k) = 'dtempdt_eddy_conv'
      lname_jk(k) = 'DTEMP/DT BY EDDY CONVERGENCE (CP)'
      units_jk(k) = 'K/DAY'
      scale_jk(k) = SDAY
      pow_jk(k) = -1
      k = k + 1
      jk_vt_dse_e = k                         ; jgrid_jk(k) = 1
      sname_jk(k) = 'vt_dse_e'
      lname_jk(k) = 'VERT. TRANS. OF DRY STATIC ENERGY BY EDDIES (CP)'
      units_jk(k) = 'W/m^2'
      scale_jk(k) = -100.*BYGRAV*BYIM
      pow_jk(k) = -1
      ia_jk(k) = ia_dga
      k = k + 1
      jk_vt_lh_eddy = k                       ; jgrid_jk(k) = 1
      sname_jk(k) = 'vt_lh_eddy'
      lname_jk(k) = 'VERTICAL TRANSPORT OF LATENT HEAT BY EDDIES (CP)'
      units_jk(k) = 'W/m^2'
      scale_jk(k) = -100.*BYGRAV*BYIM*LHE
      ia_jk(k) = ia_dga
      k = k + 1
      jk_vt_se_eddy = k                       ; jgrid_jk(k) = 1
      sname_jk(k) = 'vt_se_eddy'
      lname_jk(k) ='VERTICAL TRANSPORT OF STATIC ENERGY BY EDDIES (CP)'
      units_jk(k) = 'W/m^2'
      scale_jk(k) = -100.*BYGRAV*BYIM
      ia_jk(k) = ia_dga
      k = k + 1
      jk_tot_vt_se = k                        ; jgrid_jk(k) = 1
      sname_jk(k) = 'tot_vt_se'
      lname_jk(k) =
     &    'TOTAL LARGE SCALE VERT. TRANS. OF STATIC ENRG (CP)'
      units_jk(k) = 'W/m^2'
      scale_jk(k) = -100.*BYGRAV*BYIM
      pow_jk(k) = 1
      ia_jk(k) = ia_dga
      k = k + 1
      jk_psi_tem = k                          ; jgrid_jk(k) = 2
      sname_jk(k) = 'psi_tem'
      lname_jk(k) = 'TRANSFORMED STREAM FUNCTION (CP)'
      units_jk(k) = 'kg/s'
      scale_jk(k) = 100.*BYGRAV*XWON
      ia_jk(k) = ia_dga
      pow_jk(k) = 9
      k = k + 1
      jk_epflx_v = k                          ; jgrid_jk(k) = 1
      sname_jk(k) = 'epflx_vert_cp'
      lname_jk(k) = 'VERTICAL ELIASSEN-PALM FLUX (CP)'
      units_jk(k) = 'm^2/s^2'
      scale_jk(k) = -100.*BYGRAV*BYIM
      ia_jk(k) = ia_dga
      pow_jk(k) = -2
      k = k + 1
      jk_nt_eqgpv = k                         ; jgrid_jk(k) = 1
      sname_jk(k) = 'nt_eddy_qgpv'
      lname_jk(k) = 'NORTH. TRANS. OF EDDY Q-G POT. VORTICITY'
      units_jk(k) = '10**-6 m/s^2'
      scale_jk(k) = 1.D6
      k = k + 1
      jk_dyn_conv_eddy_geop = k               ; jgrid_jk(k) = 1
      sname_jk(k) = 'dyn_conv_eddy_geop'
      lname_jk(k) = 'DYNAMIC CONVERGENCE OF EDDY GEOPOTENTIAL'
      units_jk(k) = 'W/(m^2*mb)'
      scale_jk(k) = 1d2*BYGRAV
      pow_jk(k) = -4
      k = k + 1
      jk_nt_sheat_e = k                       ; jgrid_jk(k) = 2
      sname_jk(k) = 'nt_sheat_eddy'
      lname_jk(k) = 'NORTH. TRANS. OF SENSIBLE HEAT BY EDDIES'
      units_jk(k) = 'W/mb'
      scale_jk(k) = SHA*XWON*FIM*1d2*BYGRAV
      pow_jk(k) = 11
      k = k + 1
      jk_dyn_conv_dse = k                     ; jgrid_jk(k) = 1
      sname_jk(k) = 'dyn_conv_dse'
      lname_jk(k) = 'DYNAMIC CONVERGENCE OF DRY STATIC ENERGY'
      units_jk(k) = 'W/(m^2*mb)'
      scale_jk(k) = 1d2*BYGRAV
      pow_jk(k) = -2
      k = k + 1
      jk_seke = k                             ; jgrid_jk(k) = jgrid_ke
      sname_jk(k) = 'stand_eddy_ke'
      lname_jk(k) = 'STANDING EDDY KINETIC ENERGY'
      units_jk(k) = 'm^2/s^2'
      scale_jk(k) = .5
      k = k + 1
      jk_eke = k                              ; jgrid_jk(k) = jgrid_ke
      sname_jk(k) = 'eddy_ke'
      lname_jk(k) = 'EDDY KINETIC ENERGY'
      units_jk(k) = 'm^2/s^2'
      scale_jk(k) = .5
      k = k + 1
      jk_nt_dse_se = k                        ; jgrid_jk(k) = 2
      sname_jk(k) = 'nt_dse_stand_eddy'
      lname_jk(k) = 'NOR. TRANS. OF DRY STAT. ENERGY BY STAND. EDDIES'
      units_jk(k) = 'W/mb'
      scale_jk(k) = XWON*FIM*1d2*BYGRAV
      pow_jk(k) = 11
      k = k + 1
      jk_nt_dse_e = k                         ; jgrid_jk(k) = 2
      sname_jk(k) = 'nt_dse_eddy'
      lname_jk(k) = 'NORTH. TRANS. OF DRY STATIC ENERGY BY EDDIES'
      units_jk(k) = 'W/mb'
      scale_jk(k) = XWON*FIM*1d2*BYGRAV
      pow_jk(k) = 11
      k = k + 1
      jk_tot_nt_dse = k                       ; jgrid_jk(k) = 2
      sname_jk(k) = 'tot_nt_dse'
      lname_jk(k) = 'TOTAL NORTH. TRANSPORT OF DRY STATIC ENERGY'
      units_jk(k) = 'W/mb'
      scale_jk(k) = XWON*FIM*1d2*BYGRAV
      pow_jk(k) = 12
      k = k + 1
      jk_nt_lh_e = k                          ; jgrid_jk(k) = 2
      sname_jk(k) = 'nt_lh_e'
      lname_jk(k) = 'NORTHWARD TRANSPORT OF LATENT HEAT BY EDDIES'
      units_jk(k) = 'W/mb'
      scale_jk(k) = lhe*XWON*FIM*1d2*BYGRAV
      pow_jk(k) = 10
      k = k + 1
      jk_nt_lh_se = k
      sname_jk(k) = 'nt_lh_stand_eddy'        ; jgrid_jk(k) = 2
      lname_jk(k) = 'N. TRANSPORT OF LATENT HEAT BY STAND. EDDIES'
      units_jk(k) = 'W/mb'
      scale_jk(k) = lhe*XWON*FIM*1d2*BYGRAV
      pow_jk(k) = 9
      k = k + 1
      jk_nt_see = k                           ; jgrid_jk(k) = 2
      sname_jk(k) = 'nt_se_eddy'
      lname_jk(k) = 'NORTH.TRANSPORT OF STATIC ENERGY BY EDDIES'
      units_jk(k) = 'W/mb'
      scale_jk(k) = XWON*FIM*1d2*BYGRAV
      pow_jk(k) = 11
      k = k + 1
      jk_tot_nt_se = k                        ; jgrid_jk(k) = 2
      sname_jk(k) = 'tot_nt_se'
      lname_jk(k) = 'TOTAL NORTHWARD TRANSPORT OF STATIC ENERGY'
      units_jk(k) = 'W/mb'
      scale_jk(k) = XWON*FIM*1d2*BYGRAV
      pow_jk(k) = 12
      k = k + 1
      jk_nt_am_stand_eddy = k                 ; jgrid_jk(k) = 2
      sname_jk(k) = 'nt_u_stand_eddy'
      lname_jk(k) = 'NORTH. TRANS. ZONAL MOM. BY STAND. EDDIES'
      units_jk(k) = 'm^2/s^2'
      scale_jk(k) = 1.
      k = k + 1
      jk_nt_am_eddy = k                       ; jgrid_jk(k) = 2
      sname_jk(k) = 'nt_u_eddy'
      lname_jk(k) = 'NORTH. TRANS. ZONAL MOM. BY EDDIES'
      units_jk(k) = 'm^2/s^2'
      scale_jk(k) = 1.
      k = k + 1
      jk_tot_nt_am = k                        ; jgrid_jk(k) = 2
      sname_jk(k) = 'tot_nt_u'
      lname_jk(k) = 'TOTAL NORTH. TRANS. ZONAL MOM.'
      units_jk(k) = 'm^2/s^2'
      scale_jk(k) = 1.
      pow_jk(k) = 1
      k = k + 1
      jk_we_flx_nor = k                       ; jgrid_jk(k) = 2
      sname_jk(k) = 'we_flx_nor'
      lname_jk(k) = 'NORTHWARD WAVE ENERGY FLUX'
c      units_jk(k) = '10**11 JOULES/METER/UNIT SIGMA'
      units_jk(k) = 'm^2/s^2'
      scale_jk(k) = .25
      k = k + 1
      jk_we_flx_div = k                       ; jgrid_jk(k) = 1
      sname_jk(k) = 'we_flx_div'
      lname_jk(k) = 'DIVERGENCE OF THE WAVE ENERGY FLUX'
      units_jk(k) = 'm/s^2'
      scale_jk(k) = 1
      pow_jk(k) = -6
      k = k + 1
      jk_refr_ind_wave1 = k  !!!!! Refraction Inicies must be in order
      jgrid_jk(k) = 2
      sname_jk(k) = 'refr_ind_wave1'
      lname_jk(k) = 'REFRACTION INDEX FOR WAVE NUMBER 1'
      units_jk(k) = '10**-8 m^-2'
      k = k + 1
      sname_jk(k) = 'refr_ind_wave2'
      lname_jk(k) = 'REFRACTION INDEX FOR WAVE NUMBER 2'
      units_jk(k) = '10**-8 m^-2'
      k = k + 1
      sname_jk(k) = 'refr_ind_wave3'
      lname_jk(k) = 'REFRACTION INDEX FOR WAVE NUMBER 3'
      units_jk(k) = '10**-8 m^-2'
      k = k + 1
      sname_jk(k) = 'refr_ind_wave6'
      lname_jk(k) = 'REFRACTION INDEX FOR WAVE NUMBER 6'
      units_jk(k) = '10**-8 m^-2'
      k = k + 1
      sname_jk(k) = 'refr_ind_wave9'
      lname_jk(k) = 'REFRACTION INDEX FOR WAVE NUMBER 9'
      units_jk(k) = '10**-8 m^-2'
      k = k + 1
      jk_del_qgpv = k                         ; jgrid_jk(k) = 2
      sname_jk(k) = 'del_qgpv'
      lname_jk(k) = 'Q-G POT. VORTICITY CHANGE OVER LATITUDES'
      units_jk(k) = '1/(m*s)'
      scale_jk(k) = 1.
      pow_jk(k) = -12
      k = k + 1
      jk_wstar = k                            ; jgrid_jk(k) = 1
      sname_jk(k) = 'wstar'
      lname_jk(k) = 'W*    RESIDUAL VERTICAL VELOCITY'
      units_jk(k) = 'mb/s'
      pow_jk(k) = -5
      scale_jk(k) = 1
      k = k + 1
      jk_vstar = k                            ; jgrid_jk(k) = 2
      sname_jk(k) = 'vstar'
      lname_jk(k) = 'V* = V - D(V''TH''/DTHDP)/DP'
      units_jk(k) = 'm/s'
      pow_jk(k) = -2
c Check the count
      if (k .gt. KAJKx) then
        write (6,*) 'Increase KAJKx=',KAJKx,' to at least ',k
        call stop_model('JK_TITLES: KAJKx too small',255)
      end if

      if(.not.qIjk) then
         write(iu_Ijk, FMT='(a)') 'list of JK-fields'
         do kk = 1,k
           write (iu_Ijk, '(i3,1x,a)') kk,lname_jk(kk)
         end do
         call closeunit(iu_Ijk)
      else if(kdiag(2).gt.0) then
         Qk=.false.
   30    read (iu_Ijk,'(a)',end=40) line
         read(line,'(i3)') kk
         Qk(kk)=.true.
         go to 30
   40    continue
         do kk=1,KAJKx
           if(.not.Qk(kk)) sname_jk(kk)='skip'
         end do
         call closeunit(iu_Ijk)
       end if

      RETURN
      END SUBROUTINE JKJL_TITLEX

      SUBROUTINE DIAGJK
c     routine emptied - add content here
      END SUBROUTINE DIAGJK


      SUBROUTINE JKMAP(LNAME,SNAME,UNITS,POW10P,
     &     PM,AX,SCALET,SCALEJ,SCALEK,KMAX,JWT,J1,
     *  ARQX,SCALER,SCALJR,SCALLR)
      USE MODEL_COM, only :
     &     jm,lm
      USE DIAG_COM, only : LM_REQ
      IMPLICIT NONE

!@var units string containing output field units
      CHARACTER(LEN=50) :: UNITS
!@var lname string describing output field
      CHARACTER(LEN=50) :: LNAME
!@var sname string referencing output field
      CHARACTER(LEN=30) :: SNAME

      INTEGER :: J1,JWT,KMAX
      REAL*8 :: SCALET,SCALER,PRTFAC
      INTEGER :: POW10P
      REAL*8, DIMENSION(JM,LM) :: AX
      REAL*8, DIMENSION(JM,LM_REQ) :: ARQX
      REAL*8, DIMENSION(JM) :: SCALEJ,SCALJR
      REAL*8, DIMENSION(LM) :: SCALEK
      REAL*8, DIMENSION(LM_REQ) :: SCALLR
      REAL*8, DIMENSION(LM+LM_REQ) :: PM
c     routine emptied - add content here
      END SUBROUTINE JKMAP

      SUBROUTINE JLMAP(LNAME,SNAME,UNITS,POW10P,
     &     PL,AX,SCALET,SCALEJ,SCALEL,LMAX,JWT,J1,
     *  ARQX,SCALER,SCALJR,SCALLR)
      USE MODEL_COM, only :
     &     jm,lm
      USE DIAG_COM, only : LM_REQ
      IMPLICIT NONE

!@var units string containing output field units
      CHARACTER(LEN=50) :: UNITS
!@var lname string describing output field
      CHARACTER(LEN=50) :: LNAME
!@var sname string referencing output field
      CHARACTER(LEN=30) :: SNAME

      INTEGER :: J1,JWT,LMAX
      REAL*8 :: SCALET,SCALER,PRTFAC
      INTEGER :: POW10P
      REAL*8, DIMENSION(JM,LMAX) :: AX
      REAL*8, DIMENSION(JM,LM_REQ) :: ARQX
      REAL*8, DIMENSION(JM) :: SCALEJ,SCALJR
      REAL*8, DIMENSION(LM) :: SCALEL
      REAL*8, DIMENSION(LM_REQ) :: SCALLR
      REAL*8, DIMENSION(:) :: PL
c     routine emptied - add content here         
      END SUBROUTINE JLMAP

      SUBROUTINE JLVMAP(LNAME,SNAME,UNITS,POW10P,
     &     PL,AX,SCALET,SCALEJ,SCALEL,LMAX,JWT,J1,VWT)
      USE MODEL_COM, only :
     &     jm,lm
      USE DIAG_COM, only : LM_REQ
     *     ,jmby2
      IMPLICIT NONE

!@var units string containing output field units
      CHARACTER(LEN=50) :: UNITS
!@var lname string describing output field
      CHARACTER(LEN=50) :: LNAME
!@var sname string referencing output field
      CHARACTER(LEN=30) :: SNAME

      INTEGER :: J1,JWT,LMAX
      REAL*8 :: SCALET,PRTFAC
      INTEGER :: POW10P
      REAL*8, DIMENSION(JM,LM) :: AX,VWT
      REAL*8, DIMENSION(JM) :: SCALEJ,SCALJR
      REAL*8, DIMENSION(LM) :: SCALEL
      REAL*8, DIMENSION(LM+LM_REQ) :: PL
c     routine emptied - add content here
      END SUBROUTINE JLVMAP

      SUBROUTINE DIAGIL
c     routine emptied - add content here
      END SUBROUTINE DIAGIL


      SUBROUTINE ILMAP (sname,lname,unit,PL,AX,SCALEL,LMAX,JWT
     *     ,ISHIFT)
      USE MODEL_COM, only : im,jm,lm
      IMPLICIT NONE
      character(len=20), intent(in) :: sname,unit
      character(len=80), intent(in) :: lname
      REAL*8, DIMENSION(LM), INTENT(IN) :: PL,SCALEL
      REAL*8, DIMENSION(IM,LM), INTENT(IN) :: AX
      INTEGER, INTENT(IN) :: JWT,ISHIFT
      INTEGER :: LMAX
c     routine emptied - add content here
      END SUBROUTINE ILMAP

      SUBROUTINE DIAG7P
C****
C**** THIS ENTRY PRINTS THE TABLES
C****
      USE MODEL_COM, only :
     &     im,IDACC,JDATE,JDATE0,AMON,AMON0,JYEAR,JYEAR0,XLABEL,lrunid
      USE DIAG_COM, only : qdiag,ia_12hr,ia_inst,
     &     nwav_dag,wave,Max12HR_sequ,Min12HR_sequ,acc_period
      IMPLICIT NONE

      INTEGER, DIMENSION(44) :: IPOWER
      REAL*8, DIMENSION(120) :: POWER
      REAL*8, DIMENSION(43) :: XPOWER
      REAL*8, DIMENSION(13) :: FPE
!     Arrays for pdE
      REAL*8, DIMENSION(43+1,NWAV_DAG+1) :: FPOWER
      REAL*8, DIMENSION(41,2) :: period_e
      REAL*8, DIMENSION(nwav_dag) :: xnwav
!@var COMP_WAVE complex form of WAVE. Correct arg. to subr. MEM
      COMPLEX*16, DIMENSION(Max12HR_sequ) :: COMP_WAVE
      CHARACTER XLB*14,CLAT*16,CPRES*16,CBLANK*16,TITLEO*80
      DATA CLAT/'PERIOD EASTWARD'/,CPRES/'N'/,CBLANK/' '/

      INTEGER, PARAMETER :: MMAX=12,NUAMAX=120,NUBMAX=15

      COMMON/D7COM/LNAME,SNAME,UNITS
      CHARACTER TITLE(12)*66,LNAME(12)*50,SNAME(12)*30,UNITS(12)*50

      REAL*8, DIMENSION(12) :: SCALET
      DATA SCALET/1.,1., .1,1., .1,1., 4*1.D-3,1.D-4,1.D-5/

      REAL*8 :: BYIA12,PNU,POWX,VAR

      INTEGER ::
     &     IDACC9,K,KPAGE,KQ,KTABLE,
     &     M,MMAXP1,N,NMAX,NS,NUA,NX
      INTEGER :: ic

      NMAX=NWAV_DAG
      IDACC9=IDACC(ia_12hr)
      IF (IDACC9.LE.MMAX) RETURN
C**** PATCH NEEDED IF SEVERAL RESTART FILES WERE ACCUMULATED
      IF (IDACC(ia_inst).LE.1) GO TO 320
      IDACC9=Min12HR_sequ           ! in case a February was included
      BYIA12=1./IDACC(ia_inst)
      WAVE(:,:,:,:)=WAVE(:,:,:,:)*BYIA12
  320 CONTINUE
      IF (IDACC9.GT.Max12HR_sequ) IDACC9=Max12HR_sequ

C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
C**** NOTE: there are 41 periods.  Fpower is padded for POUT
C**** Save inverse period (day) so coordinate is monotonic
      IF(QDIAG) then
        do k=1,41
          period_e(41+1-k,1) = (k-25)  !/60.
          period_e(41+1-k,2) = (k-17)  !/60.
        end do
          call open_wp(trim(acc_period)//'.wp'//XLABEL(1:LRUNID)
     *     ,41,NMAX,0,period_e)
      do n=1,nmax
        xnwav(n) = n
      end do
      XLB=' '//acc_period(1:3)//' '//acc_period(4:12)
      fpower = -1.E30
      end if
C****
C**** OUTPUT WAVE POWER AT THE EQUATOR
C****
      MMAXP1=MMAX+1
      DO 400 KPAGE=1,2
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
      DO 390 KTABLE=1,3
      KQ=3*(KPAGE-1)+KTABLE
      TITLE(KQ)=TRIM(LNAME(KQ))//" ("//TRIM(UNITS(KQ))//") "
      WRITE (6,901) TITLE(KQ)
      DO 380 NX=1,NMAX
      N=NMAX+1-NX
      do ic=1,Max12HR_sequ
        comp_wave(ic)=cmplx( WAVE(1,ic,N,KQ) , WAVE(2,ic,N,KQ) )
      end do
      CALL MEM (COMP_WAVE,IDACC9,MMAX,NUAMAX,NUBMAX,POWER,FPE,
     *  VAR,PNU)
      POWX=.5*POWER(1)
      DO 330 NUA=2,27
  330 POWX=POWX+POWER(NUA)
      XPOWER(1)=SCALET(KQ)*POWX/26.5
      POWX=0.
      DO 340 NUA=28,34
  340 POWX=POWX+POWER(NUA)
      XPOWER(2)=SCALET(KQ)*POWX/7.
      XPOWER(3)=SCALET(KQ)*(POWER(35)+POWER(36)+POWER(37)+POWER(38))/4.
      XPOWER(4)=SCALET(KQ)*(POWER(39)+POWER(40))/2.
      DO 350 NUA=41,76
  350 XPOWER(NUA-36)=SCALET(KQ)*POWER(NUA)
      POWX=.5*POWER(1)
      DO 360 NUA=77,120
  360 POWX=POWX+POWER(NUA)
      XPOWER(41)=SCALET(KQ)*POWX/44.5
      XPOWER(42)=10.*SCALET(KQ)*VAR
      XPOWER(43)=1000.*SCALET(KQ)*(VAR-PNU)
      DO 370 NS=1,43
      IPOWER(NS)=XPOWER(NS)+.5
  370 CONTINUE
        DO NS=1,41
          FPOWER(41+1-NS,N)=XPOWER(NS)
        END DO
          FPOWER(42,N)=XPOWER(42)
          FPOWER(43,N)=XPOWER(43)
  380 WRITE (6,902) N,(IPOWER(NS),NS=1,43)
      WRITE (6,903) (FPE(M),M=1,MMAXP1)
      TITLEO=TITLE(KQ)//'*60day'//XLB
      IF(QDIAG) CALL POUT_WP(TITLEO,LNAME(KQ),SNAME(KQ),UNITS(KQ),
     *     1,NMAX,FPOWER,xnwav,CLAT,CPRES)
  390 CONTINUE
  400 CONTINUE
C****
C**** OUTPUT WAVE POWER AT 50 DEG NORTH
C****
      DO 500 KPAGE=3,4
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
      DO 490 KTABLE=1,3
      KQ=3*(KPAGE-1)+KTABLE
      TITLE(KQ)=TRIM(LNAME(KQ))//" ("//TRIM(UNITS(KQ))//") "
  410 WRITE (6,911) TITLE(KQ)
      DO 480 NX=1,NMAX
      N=NMAX+1-NX
      do ic=1,Max12HR_sequ
        comp_wave(ic)=cmplx( WAVE(1,ic,N,KQ) , WAVE(2,ic,N,KQ) )
      end do
      CALL MEM (COMP_WAVE,IDACC9,MMAX,NUAMAX,NUBMAX,POWER,FPE,
     *  VAR,PNU)
      DO 420 M=1,MMAXP1
  420 FPE(M)=1000.*SCALET(KQ)*FPE(M)
      POWX=.5*POWER(1)
      DO 430 NUA=2,45
  430 POWX=POWX+POWER(NUA)
      XPOWER(1)=SCALET(KQ)*POWX/44.5
      DO 440 NUA=46,81
  440 XPOWER(NUA-44)=SCALET(KQ)*POWER(NUA)
      XPOWER(38)=SCALET(KQ)*(POWER(82)+POWER(83))/2.
      XPOWER(39)=SCALET(KQ)*(POWER(84)+POWER(85)+POWER(86)+POWER(87))/4.
      POWX=0.
      DO 450 NUA=88,94
  450 POWX=POWX+POWER(NUA)
      XPOWER(40)=SCALET(KQ)*POWX/7.
      POWX=.5*POWER(1)
      DO 460 NUA=95,120
  460 POWX=POWX+POWER(NUA)
      XPOWER(41)=SCALET(KQ)*POWX/26.5
      XPOWER(42)=10.*SCALET(KQ)*VAR
      XPOWER(43)=1000.*SCALET(KQ)*(VAR-PNU)
      DO 470 NS=1,43
      IPOWER(NS)=XPOWER(NS)+.5
  470 CONTINUE
        DO NS=1,41
          FPOWER(41+1-NS,N)=XPOWER(NS)
        END DO
          FPOWER(42,N)=XPOWER(42)
          FPOWER(43,N)=XPOWER(43)
  480 WRITE (6,902) N,(IPOWER(NS),NS=1,43)
      WRITE (6,903) (FPE(M),M=1,MMAXP1)
      TITLEO=TITLE(KQ)//'*60day'//XLB
      IF(QDIAG) CALL POUT_WP(TITLEO,LNAME(KQ),SNAME(KQ),UNITS(KQ),
     *     2,NMAX,FPOWER,xnwav,CLAT,CPRES)
  490 CONTINUE
  500 CONTINUE
      if(qdiag) call close_wp
      RETURN
C****
  901 FORMAT ('0',30X,A64,8X,'*1/60 (1/DAY)'/'   PERIOD EASTWARD--',
     *   35('---')/' N    -2      *-3   -3.3      -4       -5    -6   -7
     *.5  -10-12-15-20-30-60    60 30 20 15 12 10    7.5    6     5
     *   4*   VAR ERR'/'   --',40('---'))
  902 FORMAT (I2,41I3,I4,I4)
  903 FORMAT ('   --',40('---')/(1X,13F10.4))
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
  911 FORMAT ('0',30X,A64,8X,'*1/60 (1/DAY)'/'   PERIOD EASTWARD--',
     *  35('---')/               ' N   *-4       -5    -6   -7.5  -10-12
     *-15-20-30-60    60 30 20 15 12 10    7.5    6     5        4
     * 3.3    3*       2    VAR ERR'/'   --',40('---'))
      END SUBROUTINE DIAG7P


      SUBROUTINE MEM (SERIES,ITM,MMAX,NUAMAX,NUBMAX,POWER,FPE,VAR,PNU)
      USE CONSTANT, only : pi
      IMPLICIT NONE
      DIMENSION C(1800),S(1800),B1(62),B2(62),A(12),AA(11),P(13)
      DIMENSION SERIES(*),POWER(*),FPE(*)
      REAL*8 ARG,PP,POWERX,P,C,S,POWER,FPE
      COMPLEX*16 CI,CSUM,A,AA,B1,B2,ANOM,ADEN
      COMPLEX*16 SERIES
      REAL*8 :: PNU,VAR
      INTEGER ::
     &     I,ITM,L,M,MMAX,MMAXP1,NU,NUA,
     &     NUAMAX,NUB,NUBMAX,NUMAX,NUTM
      CI=CMPLX(0.D0,1.D0)
      MMAXP1=MMAX+1
C**COSINE AND SINE FUNCTION
      NUMAX=NUAMAX*NUBMAX
      DO 20 NU=1,NUMAX
      ARG=2.0*PI*FLOAT(NU)/FLOAT(NUMAX)
      C(NU)=DCOS(ARG)
   20 S(NU)=DSIN(ARG)
   50 PP=0.0
      DO 60 I=1,ITM
   60 PP=PP+SERIES(I)*CONJG(SERIES(I))
      P(1)=PP/FLOAT(ITM)
      VAR=P(1)
      M=1
      B1(1)=SERIES(1)
      B2(ITM-1)=SERIES(ITM)
      DO 70 I=2,ITM-1
      B1(I)=SERIES(I)
   70 B2(I-1)=SERIES(I)
      GO TO 80
  100 DO 110 I=1,M
  110 AA(I)=A(I)
      M=M+1
      DO 120 I=1,ITM-M
      B1(I)=B1(I)-CONJG(AA(M-1))*B2(I)
  120 B2(I)=B2(I+1)-AA(M-1)*B1(I+1)
   80 ANOM=CMPLX(0.D0,0.D0)
      ADEN=CMPLX(0.D0,0.D0)
      DO 90 I=1,ITM-M
      ANOM=ANOM+CONJG(B1(I))*B2(I)
   90 ADEN=ADEN+B1(I)*CONJG(B1(I))+B2(I)*CONJG(B2(I))
      A(M)=(ANOM+ANOM)/ADEN
      P(M+1)=P(M)*(1.0-CONJG(A(M))*A(M))
      IF (M.EQ.1) GO TO 100
  130 CONTINUE
      DO 140 I=1,M-1
  140 A(I)=AA(I)-A(M)*CONJG(AA(M-I))
      IF (M.LT.MMAX) GO TO 100
C**FINAL PREDICTION ERROR
      DO 150 M=1,MMAXP1
  150 FPE(M)=P(M)*FLOAT(ITM+M-1)/FLOAT(ITM-M+1)
      DO 180 NUA=1,NUAMAX
      POWERX=0.
C**FREQUENCY BAND AVERAGE
      DO 170 NUB=1,NUBMAX
      NU=NUB+NUA*NUBMAX+(NUMAX-3*NUBMAX-1)/2
      CSUM=1.
      DO 160 M=1,MMAX
      NUTM=MOD(NU*M-1,NUMAX)+1
  160 CSUM=CSUM-A(M)*(C(NUTM)-CI*S(NUTM))
  170 POWERX=POWERX+P(MMAXP1)/(CSUM*CONJG(CSUM))
      POWER(NUA)=.5*POWERX/FLOAT(NUBMAX)
  180 CONTINUE
      PNU=0.0
      DO 210 L=1,NUAMAX
  210 PNU=PNU+POWER(L)
      PNU=PNU/(.5*NUAMAX)
      RETURN
      END SUBROUTINE MEM


      SUBROUTINE IJ_TITLEX
!@sum  IJ_TITLEX defines name,lname,units for composite ij output
!@+    the remaining attributes (value,wt,grid,range) are set in ij_MAPk
!@auth G. Schmidt/M. Kelley
!@ver  1.0
      USE DIAG_COM
      USE BDIJ
      IMPLICIT NONE
      INTEGER :: k,k1
c
      k = kaij
c
      k = k + 1

      ij_topo = k
      name_ij(k) = 'topog'
      lname_ij(k) = 'TOPOGRAPHY'
      units_ij(k) = 'METERS'

      k = k + 1
      ij_fland = k
      name_ij(k) = 'frac_land'
      lname_ij(k) = 'LAND COVERAGE'
      units_ij(k) = '%'

      k = k + 1
      ij_jet = k
      name_ij(k) = 'jet_speed'
      lname_ij(k) = 'JET SPEED'
      units_ij(k) = 'm/s'

      k = k + 1
      ij_wsmn = k
      name_ij(k) = 'rt_usmn2_vsmn2'
      lname_ij(k) = 'SURF WIND SPEED FROM Uav,Vav'
      units_ij(k) = 'm/s'

      k = k + 1
      ij_jetdir = k
      name_ij(k) = 'jet_dir'
      lname_ij(k) = 'JET DIRECTION'
      units_ij(k) = 'CW NOR'

      k = k + 1
      ij_wsdir = k
      name_ij(k) = 'srf_wind_dir'
      lname_ij(k) = 'SURFACE WIND DIRECTION'
      units_ij(k) = 'CW NOR'

      k = k + 1
      ij_netrdp = k
      name_ij(k) = 'net_rad_planet'
      lname_ij(k) = 'NET RAD. OF PLANET'
      units_ij(k) = 'W/m^2'

      k = k + 1
      ij_albp = k
      name_ij(k) = 'plan_alb'
      lname_ij(k) = 'PLANETARY ALBEDO'
      units_ij(k) = '%'

      k = k + 1
      ij_albg = k
      name_ij(k) = 'grnd_alb'
      lname_ij(k) = 'GROUND ALBEDO'
      units_ij(k) = '%'

      k = k + 1
      ij_albv = k
      name_ij(k) = 'vis_alb'
      lname_ij(k) = 'VISUAL ALBEDO'
      units_ij(k) = '%'

      k = k + 1
      ij_albgv = k
      name_ij(k) = 'grnd_alb_vis'
      lname_ij(k) = 'GROUND ALBEDO IN VISUAL RANGE'
      units_ij(k) = '%'

      k = k + 1
      ij_ntdsese = k
      name_ij(k) = 'stand_eddy_nt_dse'
      lname_ij(k) = 'NT DRY STAT ENR BY ST ED' ! NORTHWD TRANSP
      units_ij(k) = 'E14 WT'

      k = k + 1
      ij_ntdsete = k
      name_ij(k) = 'trans_eddy_nt_dse'
      lname_ij(k) = 'NT DRY STAT ENR BY TR ED' ! NORTHWD TRANSP
      units_ij(k) = 'E14 WT'

      ij_dzt1 = k+1
      do k1 = 1,kgz_max-1
        name_ij(k+k1) = 'dztemp_'//trim(pmname(k1))//
     *    '-'//trim(pmname(k1+1))
        lname_ij(k+k1) = 'THICKNESS TEMP '//trim(pmname(k1))//
     *    '-'//pmname(k1+1)
        units_ij(k+k1) = 'C'
      end do
      k = k + kgz_max -1

      k = k + 1
      ij_grow = k
      name_ij(k) = 'grow_seas'
      lname_ij(k) = 'GROWING SEASON'
      units_ij(k) = 'days'

      k = k + 1
      ij_colh2o = k
      name_ij(k) = 'pcol_h2o'
      lname_ij(k) = 'PRECIPITABLE WATER'
      units_ij(k) = 'cm'

      k = k + 1
      ij_msu2 = k
      name_ij(k) = 'Tmsu_ch2'
      lname_ij(k) = 'MSU-channel 2 TEMPERATURE'
      units_ij(k) = 'C'

      k = k + 1
      ij_msu3 = k
      name_ij(k) = 'Tmsu_ch3'
      lname_ij(k) = 'MSU-channel 3 TEMPERATURE'
      units_ij(k) = 'C'

      k = k + 1
      ij_msu4 = k
      name_ij(k) = 'Tmsu_ch4'
      lname_ij(k) = 'MSU-channel 4 TEMPERATURE'
      units_ij(k) = 'C'

      k = k + 1
      ij_Tatm = k
      name_ij(k) = 'Tatm'
      lname_ij(k) = 'ATMOSPHERIC TEMPERATURE'
      units_ij(k) = 'C'

      K = K+1
      IJ_RTSE = K
       NAME_IJ(K) = 'RTSE'
      LNAME_IJ(K) = 'THERMAL RADIATION EMITTED by SURFACE'
      UNITS_IJ(K) = 'W/m^2'

      K = K+1
      IJ_HWV = K
       NAME_IJ(K) = 'HWV'
      LNAME_IJ(K) = 'LATENT HEAT FLUX'
      UNITS_IJ(K) = 'W/m^2'

      K = K+1
      IJ_PVS = K
       NAME_IJ(K) = 'PVS'
      LNAME_IJ(K) = 'SURFACE VAPOR PRESSURE'
      UNITS_IJ(K) = 'mb'

c Check the count
      if (k .gt. kaijx) then
        write (6,*) 'Increase kaijx=',kaijx,' to at least ',k
        call stop_model('IJ_TITLES: kaijx too small',255)
      end if

      do k1 = k+1,kaijx
        write(name_ij(k1),'(a3,i3.3)') 'AIJ',k1
        lname_ij(k1) = 'unused'
        units_ij(k1) = 'unused'
      end do

      return

      END SUBROUTINE IJ_TITLEX


      subroutine IJ_MAPk (k,smap,smapj,gm,igrid,jgrid,irange,
     &     name,lname,units)
!@sum IJ_MAPk returns the map data and related terms for the k-th field
!+    (l)name/units are set in DEFACC/IJ_TITLEX but may be altered here
      USE CONSTANT, only : grav,rgas,sday,twopi,sha,kapa,bygrav,tf,undef
     *     ,teeny
      USE MODEL_COM, only : im,jm,fim,jeq,byim,DTsrc,ptop,IDACC,
     &     JHOUR,JHOUR0,JDATE,JDATE0,AMON,AMON0,JYEAR,JYEAR0,
     &     NDAY,Itime,Itime0,XLABEL,LRUNID
      USE DIAG_COM
      USE BDIJ

      IMPLICIT NONE

      REAL*8, DIMENSION(IM,JM) :: anum,adenom,smap
      REAL*8, DIMENSION(JM) :: smapj
      integer, intent(in) :: k
      integer i,j,l,k1,k2,iwt,igrid,jgrid,irange,n1,n2
      character(len=30) name,units
      character(len=80) lname
      real*8 :: gm,nh,sh, off, byiacc, scalek, an2Zan1
!@var  isumz,isumg = 1 or 2 if zon,glob sums or means are appropriate
      integer isumz,isumg

      isumz = 2 ; isumg = 2  !  default: in most cases MEANS are needed
      if (k.eq.ij_dsev) isumz = 1

c**** Find & scale the numerators and find the appropriate denominators
c****
      adenom = 1.                                             ! default
      anum = 0.

c**** the standard cases: aij(.,.,k) or aij(.,.,k)/aij(.,.,k1)
      if (k .le. kaij) then
        name = name_ij(k) ; lname = lname_ij(k) ; units = units_ij(k)
        iwt = iw_ij(k) ; jgrid = jgrid_ij(k) ; irange = ir_ij(k)
        igrid = igrid_ij(k)
c**** offsets ("  + " or "  - " in lname_ij, i.e. 2 blanks,+|-,1 blank)
        off = 0.
        k1 = index(lname_ij(k),'  - ')
        if (k1 .le. 0) k1 = index(lname_ij(k),'  + ')
        if (k1 .gt. 0) then
          if (index(lname_ij(k),'  + TF ') .gt. 0) then
            off = TF                       ! should do in accum-phase ??
          else if (index(lname_ij(k),'  - PTOP') .gt. 0)  then
            off = -ptop
          end if
          lname(k1:80) = ' '
        end if
        byiacc = 1./(idacc(ia_ij(k))+teeny)
        do j=1,jm
        do i=1,im
          anum(i,j) = aij(i,j,k)*(scale_ij(k)*byiacc) - off
        end do
        end do

c**** ratios (the denominators)
        k1 = index(lname_ij(k),' x ')
        if (k1 .gt. 0 .and. qdiag_ratios) then
          if (index(lname_ij(k),' x POPOCN') .gt. 0) then
            do j=1,jm      ! open ocean only (incl. open lake)
            do i=1,im
c              adenom(i,j) = 1-fland_glob(i,j) - aij(i,j,ij_rsoi)
c     *             /(idacc(ia_ij(ij_rsoi))+teeny)
              adenom(i,j) =  wt_ij(i,j,2)+aij(i,j,ij_lk)
     *             /(idacc(ia_ij(ij_lk))+teeny) - aij(i,j,ij_rsoi)
     *             /(idacc(ia_ij(ij_rsoi))+teeny)
            end do
            end do
          else if (index(lname_ij(k),' x POCEAN') .gt. 0) then
            do j=1,jm      ! full ocean box (no lake)
            do i=1,im
              adenom(i,j) = wt_ij(i,j,2) ! focean_glob
            end do
            end do
          else if (index(lname_ij(k),' x POICE') .gt. 0) then
            do j=1,jm      ! ice-covered only
            do i=1,im
              adenom(i,j)=aij(i,j,ij_rsoi)/(idacc(ia_ij(ij_rsoi))+teeny)
            end do
            end do
          else if (index(lname_ij(k),' x PLICE') .gt. 0) then
            do j=1,jm      ! land ice-covered only
            do i=1,im
              adenom(i,j)=aij(i,j,ij_li)/(idacc(ia_ij(ij_li))+teeny)
            end do
            end do
          else if (index(lname_ij(k),' x PSOIL') .gt. 0) then
            do j=1,jm      ! earth only (fland - flake - flice)
            do i=1,im
              adenom(i,j)=wt_ij(i,j,iw_soil)
            end do
            end do
          else if (index(lname_ij(k),' x TOTAL CLOUD') .gt. 0) then
            do j=1,jm
            do i=1,im
              adenom(i,j)=aij(i,j,ij_cldcv)/(idacc(ia_ij(ij_cldcv))
     *             +teeny)
            end do
            end do
          else if (index(lname_ij(k),' x TAU>1 CLOUD') .gt. 0) then
            do j=1,jm
            do i=1,im
              adenom(i,j)=aij(i,j,ij_cldcv1)/(idacc(ia_ij(ij_cldcv1))
     *             +teeny)
            end do
            end do
          else if (index(lname_ij(k),' x WATER CLOUD') .gt. 0) then
            do j=1,jm
            do i=1,im
              adenom(i,j)=aij(i,j,ij_wtrcld)/(idacc(ia_ij(ij_wtrcld))
     *             +teeny)
            end do
            end do
          else if (index(lname_ij(k),' x ICE CLOUD') .gt. 0) then
            do j=1,jm
            do i=1,im
              adenom(i,j)=aij(i,j,ij_icecld)/(idacc(ia_ij(ij_icecld))
     *             +teeny)
            end do
            end do
          else if (index(lname_ij(k),' x CLRSKY') .gt. 0) then
            do j=1,jm
            do i=1,im
              adenom(i,j)=1.-aij(i,j,ij_cldcv)/(idacc(ia_ij(ij_cldcv))
     *             +teeny)
            end do
            end do
          else if (index(lname_ij(k),' x TOTAL ISCCP') .gt. 0) then
            do j=1,jm
            do i=1,im
              adenom(i,j)=aij(i,j,ij_tcldi)/(idacc(ia_ij(ij_tcldi))
     *             +teeny)
            end do
            end do
          else if (index(lname_ij(k),' x SUNLIT ISCCP') .gt. 0) then
            do j=1,jm
            do i=1,im
              adenom(i,j)=aij(i,j,ij_scldi)/(idacc(ia_ij(ij_scldi))
     *             +teeny)
            end do
            end do
          else if (index(lname_ij(k),' x P1000') .gt. 0) then
            do j=1,jm
            do i=1,im
              adenom(i,j)=aij(i,j,ij_p1000)/(idacc(ia_ij(ij_p1000))
     *             +teeny)
            end do
            end do
          else if (index(lname_ij(k),' x P925') .gt. 0) then
            do j=1,jm
            do i=1,im
              adenom(i,j)=aij(i,j,ij_p925)/(idacc(ia_ij(ij_p925))
     *             +teeny)
            end do
            end do
          else if (index(lname_ij(k),' x P850') .gt. 0) then
            do j=1,jm
            do i=1,im
              adenom(i,j)=aij(i,j,ij_p850)/(idacc(ia_ij(ij_p850))
     *             +teeny)
            end do
            end do
          else if (index(lname_ij(k),' x P700') .gt. 0) then
            do j=1,jm
            do i=1,im
              adenom(i,j)=aij(i,j,ij_p700)/(idacc(ia_ij(ij_p700))
     *             +teeny)
            end do
            end do
          else if (index(lname_ij(k),' x P600') .gt. 0) then
            do j=1,jm
            do i=1,im
              adenom(i,j)=aij(i,j,ij_p600)/(idacc(ia_ij(ij_p600))
     *             +teeny)
            end do
            end do
          else if (index(lname_ij(k),' x P500') .gt. 0) then
            do j=1,jm
            do i=1,im
              adenom(i,j)=aij(i,j,ij_p500)/(idacc(ia_ij(ij_p500))
     *             +teeny)
            end do
            end do
          else if (index(lname_ij(k),' x LKICE') .gt. 0) then
            do j=1,jm
            do i=1,im
              adenom(i,j)=aij(i,j,ij_lkice)/(idacc(ia_ij(ij_lkice))
     *             +teeny)
            end do
            end do
          end if
          lname(k1:80) = ' '
        end if
        go to 100
      end if

c**** compound quantities defined with their attributes (k > kaij)
c****
      iwt = iw_all ; igrid = 1; jgrid = 1 ; irange = ir_pct   ! defaults
      name  = name_ij(k)
      lname = lname_ij(k) ; units = units_ij(k)

c**** time independent arrays
      if      (k.eq.ij_topo)  then
        anum = zatmo_glob*bygrav    ; irange = ir_0_3550

      else if (k.eq.ij_fland) then
        anum = 100.*wt_ij(:,:,iw_land)

c**** vectors: magnitude
      else if (k.eq.ij_jet.or.k.eq.ij_wsmn) then
          igrid = 2
          jgrid = 2 ;  n1 = ij_ujet ; n2 = ij_vjet ; irange = ir_0_71
        if (k.eq.ij_wsmn) then
          igrid = 1
          jgrid = 1 ;  n1 = ij_us   ; n2 = ij_vs   ; irange = ir_0_18
        end if
        byiacc=1./(idacc(ia_ij(n1))+teeny)
        do j=1,jm
        do i=1,im
          anum(i,j)=sqrt(aij(i,j,n1)**2+aij(i,j,n2)**2)*
     *      (scale_ij(n1)*byiacc)
        end do
        end do

c**** vectors: direction clockwise north (-180 -> 180)
      else if (k.eq.ij_jetdir .or. k.eq.ij_wsdir) then
        irange = ir_angl
        if (k.eq.ij_jetdir) then
          igrid = 2 ; jgrid = 2 ; n1 = ij_ujet ; n2 = ij_vjet
        else if (k.eq.ij_wsdir)  then
          igrid = 1 ; jgrid = 1 ; n1 = ij_us   ; n2 = ij_vs
        end if
        do j=2,jm
        do i=1,im
          anum(i,j)=360.*atan2(aij(i,j,n1)+teeny,aij(i,j,n2)+teeny)
     *         /twopi
        end do
        end do

c**** linear combinations (sums, differences, etc)
      else if (k.eq.ij_netrdp) then
        irange = ir_m530_190
        byiacc = 1./(idacc(ia_ij(ij_trnfp0))+teeny)
        do j=1,jm
        do i=1,im
          anum(i,j)=(aij(i,j,ij_trnfp0)+aij(i,j,ij_srnfp0))*byiacc
        end do
        end do

c**** ratios of lin. comb.: albedos from net radiation
      else if (k.eq.ij_albp .or. k.eq.ij_albg) then
        n1=ij_srnfp0
        n2=ij_srincp0
        if (k.eq.ij_albg) then
          n2=ij_srincg
          n1=ij_srnfg
        end if
        an2Zan1=idacc(ia_ij(n2))/(idacc(ia_ij(n1))+teeny)
        do j=1,jm
        do i=1,im
          adenom(i,j)=aij(i,j,n2)
          anum(i,j)=100.*(adenom(i,j)-aij(i,j,n1)*an2Zan1)
        end do
        end do

c**** ratios: albedos from reflected radiation
      else if (k.eq.ij_albv .or. k.eq.ij_albgv) then
        n1=ij_srref
        n2=ij_srincp0
        if (k.eq.ij_albgv) then
          n1=ij_srvis
          n2=ij_srincp0
        end if
        an2Zan1=idacc(ia_ij(n2))/(idacc(ia_ij(n1))+teeny)
        do j=1,jm
        do i=1,im
          anum(i,j)=100.*aij(i,j,n1)*an2Zan1
          adenom(i,j)=aij(i,j,n2)
        end do
        end do

c**** precomputed fields: northward tranports by eddies, Tmsu
      else if (k.eq.ij_ntdsese) then                   ! standing eddies
        byiacc=1./(idacc(ia_ij(ij_dsev))+teeny)   ; irange = ir_m95_265
        anum=SENTDSE*(byiacc*scale_ij(ij_dsev))  ;  igrid = 2; jgrid = 2
        isumz = 1 ; isumg = 2

      else if (k.eq.ij_ntdsete) then                  ! transient eddies
        byiacc=1./(idacc(ia_ij(ij_dsev))+teeny)   ; irange = ir_m1_3
        anum=TENTDSE*(byiacc*scale_ij(ij_dsev))  ;  igrid = 2; jgrid = 2
        isumz = 1 ; isumg = 2

      else if (k.eq.ij_msu2) then                   ! T_msu_ch2
        anum=tmsu2  ; igrid = 2; jgrid = 2 ; irange = ir_m80_28

      else if (k.eq.ij_msu3) then                   ! T_msu_ch3
        anum=tmsu3  ; igrid = 2; jgrid = 2 ; irange = ir_m80_28

      else if (k.eq.ij_msu4) then                   ! T_msu_ch4
        anum=tmsu4  ; igrid = 2; jgrid = 2 ; irange = ir_m80_28

c**** group of kgz_max-1 thickness temperatures (from heights)
      else if (k.ge.ij_dzt1 .and. k.le.ij_dzt1+kgz_max-2) then
        byiacc = 1./(idacc(ia_ij(ij_phi1k))+teeny) ; irange = ir_m80_28
        k1 = k-ij_dzt1+1  ; k2 = ij_phi1k + k1
        scalek = 1./(rgas*log(pmb(k1)/pmb(k1+1)))
        off = (ght(k1+1)-ght(k1)) * grav
        do j=1,jm
        do i=1,im
          anum(i,j)=((aij(i,j,k2)-aij(i,j,k2-1))*byiacc + off)*scalek-tf
        end do
        end do

c**** length of growing season   (not quite right ???)
      else if (k.eq.ij_grow) then
        byiacc = 1./(idacc(ia_inst)+teeny) ; irange = ir_0_180
        do j=1,jm
        do i=1,im
          anum(i,j)=(tsfrez(i,j,tf_last)-tsfrez(i,j,tf_day1))*byiacc
        end do
        end do

c**** precipitable water
      else if (k.eq.ij_colh2o) then
        igrid = 2; jgrid = 2; irange = ir_ij(ij_prec)
        byiacc = .1*100.*bygrav/(idacc(ia_dga)+teeny)
        anum = 0.
        do l=1,lm
        do j=1,jm
        do i=1,im
          anum(i,j) = anum(i,j) + aijk(i,j,l,ijk_q)
        end do
        end do
        end do
        anum = anum*byiacc

c**** column atmospheric temperature
      else if (k.eq.ij_tatm) then
        do j=2,jm
        do i=1,im
          anum(i,j) = sum(aijk(i,j,1:lm,ijk_t))/
     /                sum(aijk(i,j,1:lm,ijk_dp)) - TF
        end do
        end do

C**** Thermal Radiation Emitted by Surface (W/m^2)
      elseif (K == IJ_RTSE) then
        ANUM(:,:) = (AIJ(:,:,IJ_TRSUP) - AIJ(:,:,IJ_TRSDN)) /
     /              IDACC(IA_IJ(IJ_TRSUP))

C**** Water Vapor (latent) Heat flux (W/m^2)
      elseif (K == IJ_HWV) then
        ANUM(:,:) = AIJ(:,:,IJ_EVAP) * 2500000 /
     /              (IDACC(IA_IJ(IJ_EVAP)) * DTsrc)

C**** Surface Vapor Pressure (mb)
      elseif (K == IJ_PVS) then
        ANUM(:,:) = (AIJ(:,:,IJ_PRES) / IDACC(IA_IJ(IJ_PRES)) + PTOP) *
     *              (AIJ(:,:,IJ_QS  ) / IDACC(IA_IJ(IJ_QS  )))

      else  ! should not happen
        write (6,*) 'no field defined for ij_index',k
        call stop_model('ij_mapk: undefined extra ij_field',255)
      end if

c**** Find final field and zonal, global, and hemispheric means
  100 call ij_avg (anum,adenom,wt_ij(1,1,iwt),jgrid,isumz,isumg,  ! in
     *             smap,smapj,gm,nh,sh)                    ! out

c**** fill in some key numbers
      if (k .eq. IJ_RSIT) call keyij(gm,nh)

      return

      end subroutine ij_mapk


      subroutine ij_avg (anum,aden,wtij,jgrid,isumz,isumg,
     *                   smap,smapj,gm,nh,sh)
!@sum ij_avg finds num/den and various averages from num and den
!@auth R.Ruedy
!@ver  1.0
      USE MODEL_COM, only :  im,jm

      IMPLICIT NONE

      real*8, dimension(im,jm) :: anum,aden,wtij,smap
      real*8, dimension(jm) :: smapj
      real*8, dimension(2) :: znumh,zdenh
      real*8  gm,nh,sh, sumj,wt

      integer k,i,j,jgrid,isumz,isumg
c     routine emptied - add content here

      return
      end subroutine ij_avg


      SUBROUTINE DIAGIJ
c     routine emptied - add content here
      END SUBROUTINE DIAGIJ


      subroutine maptxt (smap,smapj,gm,irange,title,line,kcol,nlines)
!@sum  maptxt prints a maplet onto 1/3 of a virtual half-page line(1-51)
!@auth R.Ruedy
!@ver  1.0
      use constant, only : undef
      use model_com, only : im,jm
      use diag_com, only : inci,incj
      use bdij

      IMPLICIT NONE

      real*8, dimension(im,jm) :: smap
      real*8, dimension(jm) :: smapj
      real*8  gm,val,zmax,zmin
      integer irange,kcol,k,k1,ifm,nz1,nz2,n1,n2,ibar,i,j
     *  ,nlines
      character(len=133), dimension(53) :: line
      character(len=40) :: title

c**** find first and last column for zonal means nz1,nz2 and maps n1,n2
      nz1 = 2 + (kcol-1)*(9+im/inci) ; nz2 = nz1 + 4
      n1  = nz2 + 2  ;  n2 = n1 + im/inci - 1
c**** pick color bar and format for zonal mean according to the range
      ibar = ib_of_legnd(irange)
      zmax = 999.9 ; zmin = -99.9
      if (kcol.gt.1) then
        zmax = 99999.9 ; zmin = -9999.9
        nz1 = nz1 - 2
      end if
      ifm = 2
      do j=1,jm
        if (smapj(j).eq.undef) cycle
        if (smapj(j).gt.zmax .or. smapj(j).lt.zmin) ifm = 1
      end do

c**** title on line 1
      line(1)(n1-4:n2) = title(1:40) ; line(1)(1:1) = '0'  ! line feed
c**** use 2 lines for each shown latitude because of overstrike
      k=0
      do j=jm,1,-incj
        k = k+2
c**** zonal mean
        val=0.
        if (smapj(j) .ne. undef) val = smapj(j)
        if (kcol.eq.1) then
          if (ifm.eq.1) write(line(k)(nz1:nz2),'(i5)') nint(min(1d5,max(
     *         -1d5,val)))
          if (ifm.eq.2) write(line(k)(nz1:nz2),'(f5.1)') val
        else
          if (ifm.eq.1) write(line(k)(nz1:nz2),'(i7)') nint(min(1d5,max(
     *         -1d5,val)))
          if (ifm.eq.2) write(line(k)(nz1:nz2),'(f7.1)') val
        end if
c**** mark selected longitudes for each selected latitude
        k1 = n1-1
        do i=1,im,inci
          k1 = k1+1
          val = undef
          if (smap(i,j).ne.undef) val = smap(i,j)*fac_legnd(irange)
c**** for angles, change range from -180->180 to 0-360 (fac=.1)
          if (irange.eq.ir_angl .and. val.lt.-.5) val = val+36.
          if (irange.eq.ir_angl .and. val.le..1) val = .1
          line(k)(k1:k1) = mark(val,ibar,undef)
          if (wt_ij(i,j,iw_land).gt..5) line(k+1)(k1:k1)=line(k)(k1:k1)
          line(k+1)(1:1) = '+'             !  overstrike
        end do
      end do
c**** below map, show global mean and mark each quarter with a '+'
      k = k+2
      if (gm .ge. zmin .and. gm .le. zmax) then
        if (kcol.eq.1) write(line(k)(nz1:nz2),'(f5.1)') gm
        if (kcol.gt.1) write(line(k)(nz1:nz2),'(f7.1)') gm
      else
        if (gm.eq.undef) then
          if (kcol.eq.1) write(line(k)(nz1:nz2),'(a)') "Undef"
          if (kcol.gt.1) write(line(k)(nz1:nz2),'(a)') "  Undef"
        else
          if (kcol.eq.1) write(line(k)(nz1:nz2),'(i5)') nint(min(1d5
     *         ,max(-1d5,gm)))
          if (kcol.gt.1) write(line(k)(nz1:nz2),'(i7)') nint(min(1d5
     *         ,max(-1d5,gm)))
        end if
      end if
      do k1 = n1,n2,im/(4*inci)
        line(k)(k1:k1) = '+'
      end do
c**** last line: legend (there is a little more space in col 1 and 2)
      if(kcol.lt.3) n2 = n1 + 39
      line(k+1)(n1:n2) = legend(irange)
      nlines = k+1

      return
      end subroutine maptxt


      SUBROUTINE IJMAP (title,smap,smapj,jgrid)
C**** Print out full-page digital maps
      USE MODEL_COM, only :
     &     im,jm
      IMPLICIT NONE
      CHARACTER*48 TITLE
      REAL*8, DIMENSION(IM,JM) :: SMAP
      REAL*8, DIMENSION(JM) :: SMAPJ
      INTEGER :: I,J,jgrid
c     routine emptied - add content here
      END SUBROUTINE IJMAP


      subroutine set_ijout (nmaplets,nmaps,Iord,Qk,iu_Iij)
!@sum set_ijout either lists or sets the fields to be processed
!@auth Reto A. Ruedy
!@ver  1.0
      USE DIAG_COM
      USE BDIJ
      use filemanager

      IMPLICIT NONE
      character*80 line
      logical Qk(kaijx),Qktmp(kaijx)
      INTEGER Iord(kaijx+10),nmaplets,nmaps,iu_Iij,k,
     *   n,kmap(3)

c**** Just list what's available - then do same for ijk-fields
      if (kdiag(3) .eq. 0) then
        Qktmp = Qk
        call openunit('Iij',iu_Iij,.false.,.false.)
        write (iu_Iij,'(a)') 'List of fields shown as maplets'
        do n=1,nmaplets
          k = Iord(n)
          Qktmp(k) = .false.
          if (k.le.0) then
             write (iu_Iij,'(i3,1x,a)') k, '  blank maplet'
          else
             write (iu_Iij,'(i3,1x,a)') k,lname_ij(k)
          end if
        end do
        write (iu_Iij,'(a)') 'List of fields shown as 1-pg maps'
        do n=nmaplets+1,nmaplets+nmaps
          k = Iord(n)
          Qktmp(k) = .false.
          if (k.le.0) then
             cycle
          else
             write (iu_Iij,'(i3,1x,a)') k,lname_ij(k)
          end if
        end do
        write (iu_Iij,'(a)') 'List of other fields in binary output'
        do k=1,kaijx
          if (.not.Qktmp(k)) cycle
          write (iu_Iij,'(i3,1x,a)') k,lname_ij(k)
        end do
        kdiag(3)=9
        CALL IJKMAP (iu_Iij)
        kdiag(3)=0
        return
      end if

c**** Redefine nmaplets,nmaps,Iord,Qk if  kdiag(3) > 0
      call openunit('Iij',iu_Iij,.false.,.true.)

      nmaplets = 0 ; nmaps = 0 ; Iord = 0 ; Qk = .false.

      kmap = 0 ; n=0 ; k=0
   10 read (iu_Iij,'(a)',end=20) line
      if (line(1:1) .eq. 'l') go to 20
      if (line(1:1) .eq. 'L') then
        n=n+1
        go to 10
      end if
      k = k+1
      read(line,'(i3)') Iord(k)
      if (Iord(k).gt.0) Qk(Iord(k)) = .true.
      kmap(n) = kmap(n) + 1
      go to 10

   20 nmaplets = kmap(1) ; nmaps = kmap(2)
      if (.not.qdiag .or. kdiag(3).eq.1) call closeunit(iu_Iij)
      return
      end subroutine set_ijout


      SUBROUTINE DIAGCP
c     routine emptied - add content here
      END SUBROUTINE DIAGCP


      SUBROUTINE DIAG5P
c     routine emptied - add content here
      END SUBROUTINE DIAG5P


      SUBROUTINE DIAGDD
!@sum  DIAGDD prints out diurnal cycle diagnostics
!@auth G. Russell
!@ver  1.0
      USE MODEL_COM, only :
     &     idacc,JDATE,JDATE0,AMON,AMON0,JYEAR,JYEAR0,XLABEL,LRUNID,NDAY
      USE DIAG_COM, only :   kdiag,qdiag,acc_period,units_dd,ndiupt,
     &     adiurn,ijdd,namdd,ndiuvar,hr_in_day,scale_dd,lname_dd,name_dd
     *     ,ia_12hr
      IMPLICIT NONE

      REAL*8, DIMENSION(HR_IN_DAY+1) :: XHOUR
      INTEGER, DIMENSION(HR_IN_DAY+1) :: MHOUR
      REAL*8 :: AVE,AVED,AVEN,BYIDAC
      INTEGER :: I,IH,IREGF,IREGL,IS,K,KP,KQ,KR,NDAYS,KF,KNDIU,KR1,KR2
      CHARACTER*16, DIMENSION(NDIUVAR) :: UNITSO,LNAMEO,SNAMEO
      REAL*8, DIMENSION(HR_IN_DAY+1,NDIUVAR) :: FHOUR
      CHARACTER :: CPOUT*2
C****
      NDAYS=IDACC(ia_12hr)/2
      IF (NDAYS.LE.0) RETURN
      BYIDAC=24./(NDAY*NDAYS)
C****
      IREGF=1
      IREGL=NDIUPT-KDIAG(6)       ! kd6=KDIAG(6)>0: skip last kd6 points
      IF (KDIAG(6).LT.0.AND.KDIAG(6).GE.-NDIUPT) IREGF=-KDIAG(6)
      IF (KDIAG(6).LT.0) IREGL=IREGF       ! kd6<0: show only point -kd6
C**** for netcdf limits, loop in steps of 2000
      KNDIU=0
      DO KQ=1,NDIUVAR
        IF (LNAME_DD(KQ) == "unused") CYCLE
        KNDIU=KNDIU+1
      END DO
      DO KF=1,1+(KNDIU*(IREGL-IREGF+1)-1)/2000
C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      KR1=IREGF+(KF-1)*INT(2000/KNDIU)
      KR2=MIN(IREGL,IREGF+KF*INT(2000/KNDIU)-1)
      IF (QDIAG) THEN
        CPOUT=""
        IF (KNDIU*(IREGL-IREGF+1)/2000 > 1) THEN ! more than one file
          IF (KF <= 9) THEN
            WRITE(CPOUT(1:1),'(I1)') KF
          ELSE
            WRITE(CPOUT(1:2),'(I2)') KF
          END IF
        END IF
        call open_diurn (trim(acc_period)//'.diurn'//trim(cpout)
     *      //XLABEL(1:LRUNID),hr_in_day,KNDIU,KR1,KR2)
      END IF
C**** LOOP OVER EACH BLOCK OF DIAGS
      DO KR=KR1,KR2
        WRITE (6,901) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
        WRITE (6,903) NAMDD(KR),IJDD(1,KR),IJDD(2,KR),(I,I=1,HR_IN_DAY)
C**** KP packs the quantities for postprocessing (skipping unused)
        KP = 0
        DO KQ=1,NDIUVAR
          IF (MOD(KQ-1,5).eq.0) WRITE(6,*)
          IF (LNAME_DD(KQ).eq."unused") CYCLE
          KP = KP+1
          SELECT CASE (NAME_DD(KQ))
          CASE DEFAULT
C**** NORMAL QUANTITIES
            AVE=0.
            DO IH=1,HR_IN_DAY
              AVE=AVE+ADIURN(KQ,KR,IH)
              XHOUR(IH)=ADIURN(KQ,KR,IH)*SCALE_DD(KQ)*BYIDAC
            END DO
            XHOUR(HR_IN_DAY+1)=AVE/FLOAT(HR_IN_DAY)*SCALE_DD(KQ)*BYIDAC
C**** RATIO OF TWO QUANTITIES
          CASE ('LDC')
            AVEN=0.
            AVED=0.
            DO IH=1,HR_IN_DAY
              AVEN=AVEN+ADIURN(KQ,KR,IH)
              AVED=AVED+ADIURN(KQ-1,KR,IH)
              XHOUR(IH)=ADIURN(KQ,KR,IH)*SCALE_DD(KQ)/
     *             (ADIURN(KQ-1,KR,IH)+1D-20)
            END DO
            XHOUR(HR_IN_DAY+1)=AVEN*SCALE_DD(KQ)/(AVED+1D-20)
          END SELECT
          DO IS=1,HR_IN_DAY+1
            FHOUR(IS,KP)=XHOUR(IS)
            MHOUR(IS)=NINT(XHOUR(IS))
          END DO
          WRITE (6,904) LNAME_DD(KQ),MHOUR
          SNAMEO(KP)=NAME_DD(KQ)(1:16)
          LNAMEO(KP)=LNAME_DD(KQ)(1:16)
          UNITSO(KP)=UNITS_DD(KQ)(1:16)
        END DO
        IF (QDIAG) CALL POUT_DIURN(SNAMEO,LNAMEO,UNITSO,FHOUR,
     *       NAMDD(KR),IJDD(1,KR),IJDD(2,KR),HR_IN_DAY,KP)
      END DO
      IF (QDIAG) call close_diurn
      END DO

      RETURN
C****
  901 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
  903 FORMAT ('0',A4,I2,',',I2,' ',I2,23I5,'  AVE')
  904 FORMAT (A8,25I5)
      END SUBROUTINE DIAGDD

      SUBROUTINE DIAGDH
!@sum  DIAGDH prints out hourly diurnal cycle diagnostics
!@+       It uses the same quantities as DIAGHH and shares some arrays
!@+       When radiation is not called every hour this will not average
!@+       exactly to same numbers as in DIAGDD.
!@auth J. Lerner
!@ver  1.0
#ifndef NO_HDIURN
      USE MODEL_COM, only :   JDendOfM,JMON,NDAY,
     &     idacc,JDATE,JDATE0,AMON,AMON0,JYEAR,JYEAR0,XLABEL,LRUNID
      USE DIAG_COM, only :   kdiag,qdiag,acc_period,units_dd,hr_in_month
     *     ,hdiurn,ijdd,namdd,ndiuvar,hr_in_day,scale_dd,lname_dd
     *     ,name_dd,ia_12hr,NDIUPT
      IMPLICIT NONE
      REAL*8, DIMENSION(HR_IN_MONTH) :: XHOUR
      INTEGER, DIMENSION(HR_IN_MONTH) :: MHOUR
      INTEGER :: I,IH,IH0,IREGF,IREGL,IS,JD,jdayofm,K,KP,KQ,KR,NDAYS,KF,
     &     KNDIU,KR1,KR2
      CHARACTER*16, DIMENSION(NDIUVAR) :: UNITSO,LNAMEO,SNAMEO
      REAL*8, DIMENSION(HR_IN_MONTH,NDIUVAR) :: FHOUR
      CHARACTER :: CPOUT*2
C****
      NDAYS=IDACC(ia_12hr)/2
      IF (NDAYS.LE.0) RETURN
C****
C**** KP packs the quantities for postprocessing (skipping unused)
      jdayofM = JDendOfM(jmon)-JDendOfM(jmon-1)
      IREGF=1
      IREGL=NDIUPT-KDIAG(13)      ! kd13=KDIAG(13)>0: skip last kd13 pts
      IF (KDIAG(13).LT.0.AND.KDIAG(13).GE.-NDIUPT) IREGF=-KDIAG(13)
      IF (KDIAG(13).LT.0) IREGL=IREGF       ! kd13<0: show only pt -kd13
C**** for netcdf limits, loop in steps of 2000
      KNDIU=0
      DO KQ=1,NDIUVAR
        IF (LNAME_DD(KQ) == "unused") CYCLE
        KNDIU=KNDIU+1
      END DO
      DO KF=1,1+(KNDIU*(IREGL-IREGF+1)-1)/2000
      KR1=IREGF+(KF-1)*INT(2000/KNDIU)
      KR2=MIN(IREGL,IREGF+KF*INT(2000/KNDIU)-1)
      IF (QDIAG) THEN
        CPOUT=""
        IF (KNDIU*(IREGL-IREGF+1)/2000 > 1) THEN ! more than one file
          IF (KF <= 9) THEN
            WRITE(CPOUT(1:1),'(I1)') KF
          ELSE
            WRITE(CPOUT(1:2),'(I2)') KF
          END IF
        END IF
C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
        call open_hdiurn (trim(acc_period)//'.hdiurn'//trim(cpout)
     &       //XLABEL(1:LRUNID),hr_in_month,KNDIU,KR1,KR2)
      END IF
C**** LOOP OVER EACH BLOCK OF DIAGS
      DO KR=KR1,KR2
        WRITE (6,901) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
        WRITE (6,903)NAMDD(KR),IJDD(1,KR),IJDD(2,KR),(I,I=1,HR_IN_DAY)
C**** KP packs the quantities for postprocessing (skipping unused)
        KP = 0
        DO KQ=1,NDIUVAR
          IF (MOD(KQ-1,5).eq.0) WRITE(6,*)
          IF (LNAME_DD(KQ).eq."unused") CYCLE
          KP = KP+1
          SELECT CASE (NAME_DD(KQ))
          CASE DEFAULT
C**** NORMAL QUANTITIES
            DO IH=1,HR_IN_MONTH
              XHOUR(IH)=HDIURN(KQ,KR,IH)*SCALE_DD(KQ)*(24./NDAY)
            END DO
C**** RATIO OF TWO QUANTITIES
          CASE ('LDC')
            DO IH=1,HR_IN_MONTH
              XHOUR(IH)=HDIURN(KQ,KR,IH)*SCALE_DD(KQ)/
     *             (HDIURN(KQ-1,KR,IH)+1D-20)
            END DO
          END SELECT
          DO IS=1,HR_IN_MONTH
            FHOUR(IS,KP)=XHOUR(IS)
            MHOUR(IS)=NINT(XHOUR(IS))
          END DO
          ih0 = 1
          do jd = 1,jdayofm
            WRITE (6,904) LNAME_DD(KQ),(MHOUR(i),i=ih0,ih0+23),jd
            ih0 = ih0+24
          end do
          SNAMEO(KP)=NAME_DD(KQ)(1:16)
          LNAMEO(KP)=LNAME_DD(KQ)(1:16)
          UNITSO(KP)=UNITS_DD(KQ)(1:16)
        END DO
        IF (QDIAG) CALL POUT_HDIURN(SNAMEO,LNAMEO,UNITSO,FHOUR,
     *     NAMDD(KR),IJDD(1,KR),IJDD(2,KR),HR_IN_MONTH,KP)
      END DO
      IF (QDIAG) call close_hdiurn
      END DO
#endif
      RETURN
C****
  901 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
  903 FORMAT ('0',A4,I2,',',I2,' ',I2,23I5,'  Day')
  904 FORMAT (A8,24I5,I5)
      END SUBROUTINE DIAGDH


      SUBROUTINE DIAG4
!@sum  DIAG4 prints out a time history of the energy diagnostics
!@auth G. Russell
!@ver  1.0
      USE CONSTANT, only :
     &     grav,rgas,bygrav
      USE MODEL_COM, only :
     &     im,jm,lm,fim,
     &     IDACC,JHOUR,JHOUR0,JDATE,JDATE0,AMON,AMON0,
     &     JYEAR,JYEAR0,NDA4,NDAY,Itime0,XLABEL,istrat
      USE DIAG_COM, only :
     &     energy,ned,nehist,hist_days,xwon,ia_inst,ia_d4a
      IMPLICIT NONE

      REAL*8, DIMENSION(2) :: FAC
      REAL*8, DIMENSION(NED) :: SCALET
      REAL*8, DIMENSION(2*NED) :: SUME
      INTEGER, DIMENSION(2*NED) :: IK
      REAL*8, DIMENSION(NEHIST,HIST_DAYS+1) :: EHIST

      INTEGER ::
     &     I,IDACC5,ItimeX,IDAYX,IDAYXM,K,K0,KS,KN,KSPHER
      REAL*8 :: TOFDYX

      IDACC5=IDACC(ia_d4a)
      IF (IDACC5.LE.0) RETURN
      IF (IDACC(ia_inst).LT.1) IDACC(ia_inst)=1
      SCALET(1)=100.D-18*BYGRAV
      SCALET(2)=SCALET(1)
      SCALET(3)=SCALET(1)
      SCALET(4)=SCALET(1)
c     SCALET(5)=.5*SCALET(1)
      SCALET(5)=SCALET(1)
      SCALET(6)=SCALET(5)
      SCALET(7)=SCALET(1)*RGAS
      SCALET(8)=SCALET(7)
      SCALET(9)=SCALET(7)
      SCALET(10)=SCALET(7)
      DO K=1,NED
        SCALET(K)=XWON*SCALET(K)/IDACC(ia_inst)
      END DO
C****
      DO K0=1,MIN(1+ISTRAT,2)
        WRITE (6,901) XLABEL
        IF (K0.eq.1) THEN
          FAC(1) = 1.
          FAC(2) = 10.  ! a factor of 10 for LOW STRAT
          WRITE (6,902) JYEAR0,AMON0,JDATE0,JHOUR0,JYEAR,AMON,JDATE
     *         ,JHOUR
          WRITE (6,903)
        ELSE
          FAC(1) = 10.  ! 10 goes from 10^18 to 10^17
          FAC(2) = 100. ! another factor of 10 for HIGH STRAT
          WRITE (6,906) JYEAR0,AMON0,JDATE0,JHOUR0,JYEAR,AMON,JDATE
     *         ,JHOUR
          WRITE (6,907)
        END IF
        SUME(:)=0.
        DO I=1,IDACC5
          ItimeX=Itime0+I*NDA4-1
          IDAYX=1+ItimeX/NDAY
          IDAYXM=MOD(IDAYX,100000)
          TOFDYX=MOD(ItimeX,NDAY)*24./NDAY
          DO KSPHER=1,2
            DO K=1,NED
              KS=K+(KSPHER-1)*NED
              KN=KS+(K0-1)*2*NED
              IF (KN.le.NEHIST) THEN
                EHIST(KN,I)=ENERGY(KN,I)*SCALET(K)*FAC(KSPHER)
                IK(KS)=EHIST(KN,I)+.5
                SUME(KS)=SUME(KS)+ENERGY(KN,I)
              ELSE
                IK(KS)=-999
              END IF
            END DO
          END DO
          WRITE (6,904) IDAYXM,TOFDYX,IK
        END DO
        DO KSPHER=1,2
          DO K=1,NED
            KS=K+(KSPHER-1)*NED
            KN=KS+(K0-1)*2*NED
            IF (KN.le.NEHIST) THEN
              EHIST(KN,HIST_DAYS+1)=SUME(KS)*SCALET(K)*FAC(KSPHER)
     *             /IDACC5
              IK(KS)=EHIST(KN,HIST_DAYS+1)+.5
            ELSE
              IK(KS)=-999
            END IF
          END DO
        END DO
        WRITE (6,905) IK
        IF (K0.eq.1) CALL KEYD4 (IK)
      END DO
      RETURN
C****
  901 FORMAT ('1',A)
  902 FORMAT ('0** ENERGY HISTORY **      From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,
     *  '    UNITS OF 10**18 JOULES')
  903 FORMAT ('0',15X,21('-'),' TROPOSPHERE ',22('-'),5X,21('-'),
     *  '  LOW STRAT. * 10 ',17('-')/8X,2(11X,'ZKE',8X,'EKE',7X,
     *     'SEKE',9X,
     * 'ZPE',10X,'EPE')/3X,'DAY  HOUR     SH   NH    SH   NH     1    2
     *    SH    NH     SH    NH      SH   NH    SH   NH    SH   NH     S
     *H    NH     SH    NH'/1X,132('='))
  904 FORMAT (I6,F6.1,1X,3(I6,I5),2(I7,I6),2X,3(I6,I5),2(I7,I6))
  905 FORMAT (1X,132('=')/8X,'MEAN ',3(I6,I5),2(I7,I6),2X,3(I6,I5),
     *  2(I7,I6))
  906 FORMAT ('0** ENERGY HISTORY **      From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,
     *  '    UNITS OF 10**17 JOULES')
  907 FORMAT ('0',15X,19('-'),' MID STRATOSPHERE ',19('-'),5X,18('-'),
     *  ' HIGH STRAT. * 10  ',19('-')/8X,2(11X,'ZKE',8X,'EKE',7X,
     *   'NHKE',9X,
     *  'ZPE',10X,'EPE')/3X,'DAY  HOUR     SH   NH    SH   NH    1    2
     *    SH    NH     SH    NH      SH   NH    SH   NH    1    2      S
     *H    NH     SH    NH'/1X,132('='))
  920 FORMAT (1X)
      END SUBROUTINE DIAG4


      subroutine KEYVSUMS (QUANT,GSUM,HSUM,ASUM,SUMFAC)
      USE MODEL_COM, only : jm
      implicit none
!@var quant string designating the quantity for which to save keynrs
      CHARACTER(LEN=30) :: QUANT
      REAL*8, DIMENSION(JM) :: ASUM
      REAL*8, DIMENSION(2) :: HSUM
      REAL*8 :: GSUM,SUMFAC
      if(quant.eq.'temp') CALL KEYJKT (GSUM,ASUM)
      if(quant.eq.'eddy_ke') CALL KEYJKE (18,HSUM,ASUM)
      if(quant.eq.'tot_ke') CALL KEYJKE (19,HSUM,ASUM)
      if(quant.eq.'nt_dse_stand_eddy') CALL KEYJKN (33,ASUM,SUMFAC)
      if(quant.eq.'nt_dse_eddy') CALL KEYJKN (34,ASUM,SUMFAC)
      if(quant.eq.'tot_nt_dse') CALL KEYJKN (35,ASUM,SUMFAC)
!!!   if(quant.eq.'nt_lh_e') CALL KEYJKN (??,ASUM,SUMFAC)
!!!   if(quant.eq.'tot_nt_lh') CALL KEYJKN (??,ASUM,SUMFAC)
      if(quant.eq.'nt_se_eddy') CALL KEYJKN (36,ASUM,SUMFAC)
      if(quant.eq.'tot_nt_se') CALL KEYJKN (37,ASUM,SUMFAC)
!!!   if(quant.eq.'tot_nt_ke') CALL KEYJKN (??,ASUM,SUMFAC)
      if(quant.eq.'nt_u_stand_eddy') CALL KEYJKN (39,ASUM,SUMFAC)
      if(quant.eq.'nt_u_eddy') CALL KEYJKN (40,ASUM,SUMFAC)
      if(quant.eq.'tot_nt_u') CALL KEYJKN (41,ASUM,SUMFAC)
      RETURN
      end subroutine keyvsums


      subroutine keynrl(quant,l,flat)
      USE MODEL_COM, only : jm
      implicit none
      integer :: l
      REAL*8, DIMENSION(JM) :: FLAT
!@var quant string designating the quantity for which to save keynrs
      CHARACTER(LEN=30) :: QUANT
      if(quant.eq.'u') CALL KEYJKJ (L,FLAT)
      if(quant.eq.'psi_cp') CALL KEYJLS (L,FLAT)
      return
      end subroutine keynrl




      SUBROUTINE IJK_TITLES
!@sum  IJK_TITLES extra titles for composite ijk output
!@auth G. Schmidt/M. Kelley
!@ver  1.0
      USE CONSTANT, only : bygrav
      USE DIAG_COM, only : kaijk,kaijkx,
     *   units_ijk,name_ijk,lname_ijk,scale_ijk  !  diag_com
      IMPLICIT NONE
      INTEGER :: K
c
      k = kaijk
c
      k = k + 1
      name_ijk(k) = 'z'
      lname_ijk(k) = 'HEIGHT'
      units_ijk(k) = 'm'
      scale_ijk(k) = BYGRAV

      return
      END SUBROUTINE IJK_TITLES


      SUBROUTINE IJKMAP (iu_Iij)
!@sum  IJKMAP output 3-D constant pressure output fields
!@auth G. Schmidt
!@ver  1.0
C**** Note that since many IJK diags are weighted w.r.t pressure, all
C**** diagnostics must be divided by the accumulated pressure
C**** All titles/names etc. implicitly assume that this will be done.
C**** IJL diags are done separately
      USE CONSTANT, only : grav,sha,undef
      USE MODEL_COM, only : im,jm,lm,pmidl00,XLABEL,LRUNID,idacc
      USE DIAG_COM, only : kdiag,jgrid_ijk,
     &     aijk,acc_period,ijk_u,ijk_v,ijk_t,ijk_q,ijk_dp,ijk_dse
     *     ,scale_ijk,off_ijk,name_ijk,lname_ijk,units_ijk,kaijk,kaijkx
     *     ,ijl_cf,ijk_w,ia_rad,ia_dga
     *     ,ia_src,lh_diags
     *     ,ijl_llh,ijl_mctlh,ijl_mcdlh,ijl_mcslh
#ifdef CLD_AER_CDNC
     *    ,ijl_rewm,ijl_rews,ijl_cdwm,ijl_cdws,ijl_cwwm,ijl_cwws
     *    ,ijl_reim,ijl_reis,ijl_cdim,ijl_cdis,ijl_cwim,ijl_cwis
#endif
      use filemanager
      IMPLICIT NONE

      CHARACTER XLB*24,TITLEX*56
      CHARACTER*80 TITLEL(LM)
      REAL*8 SMAP(IM,JM,LM),SMAPJK(JM,LM),SMAPK(LM)
      REAL*8 flat,dp
      CHARACTER*8 CPRESS(LM)
      INTEGER i,j,l,kxlb,ni,k,iu_Iij
      logical, dimension (kaijkx) :: Qk

C****
C**** INITIALIZE CERTAIN QUANTITIES
C****
      call ijk_titles

      Qk = .true.
      do k=1,kaijkx
        if (lname_ijk(k).eq.'unused') Qk(k) = .false.
      end do
      if (kdiag(3).eq.9) then
         write (iu_Iij,'(a)') 'list of 3-d fields'
         do k=1,kaijkx
           if (lname_ijk(k).ne.'unused')
     *        write (iu_Iij,'(i3,1x,a)') k,lname_ijk(k)
         end do
         call closeunit(iu_Iij)
         return
      else if (kdiag(3).gt.1) then
         Qk = .false.
   10    read (iu_Iij,'(i3)',end=20) k
         Qk(k) = .true.
         go to 10
   20    continue
         call closeunit(iu_Iij)
      end if

C**** OPEN PLOTTABLE OUTPUT FILE
      call open_ijk(trim(acc_period)//'.ijk'//XLABEL(1:LRUNID),im,jm,lm)
      KXLB = INDEX(XLABEL(1:11),'(')-1
      IF(KXLB.le.0) KXLB = 10
      XLB = ' '
      XLB(1:13)=acc_period(1:3)//' '//acc_period(4:12)
      XLB(15:14+KXLB) = XLABEL(1:KXLB)
C****
C**** Complete 3D-field titles
C****
      DO L=1,LM
        WRITE(CPRESS(L),'(F8.3)') pmidl00(l)
      END DO

C**** Select fields
      DO K=1,KAIJKx
        if (.not.Qk(k).or.k.eq.ijk_dp.or.k.eq.ijl_cf) cycle
        SMAP(:,:,:) = UNDEF
        SMAPJK(:,:) = UNDEF
        SMAPK(:)    = UNDEF
        if (jgrid_ijk(k).eq.1) then
          TITLEX = lname_ijk(k)(1:17)//"   at        mb ("//
     *         trim(units_ijk(k))//")"
        else
          TITLEX = lname_ijk(k)(1:17)//"   at        mb ("//
     *         trim(units_ijk(k))//", UV grid)"
        end if
        IF (name_ijk(K).eq.'z') THEN ! special compound case
          DO L=1,LM
            DO J=2,JM
              NI = 0
              FLAT = 0.
              DO I=1,IM
                DP=AIJK(I,J,L,IJK_DP)
                IF(DP.GT.0.) THEN
                  SMAP(I,J,L) = SCALE_IJK(k)*(AIJK(I,J,L,IJK_DSE)-
     *                 SHA*AIJK(I,J,L,IJK_T))/DP
                  FLAT = FLAT+SMAP(I,J,L)
                  NI = NI+1
                END IF
              END DO
              IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
            END DO
            WRITE(TITLEX(23:30),'(A)') CPRESS(L)
            TITLEL(L) = TITLEX//XLB
          END DO
        ELSEIF (jgrid_ijk(k).eq.1 .or. name_ijk(k).eq."p") THEN ! no dp weight
          DO L=1,LM
            DO J=1,JM
              NI = 0
              FLAT = 0.
              DO I=1,IM
                SMAP(I,J,L) = SCALE_IJK(k)*AIJK(I,J,L,K)/IDACC(ia_dga)
                FLAT = FLAT+SMAP(I,J,L)
                NI = NI+1
              END DO
              IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
            END DO
            WRITE(TITLEX(23:30),'(A)') CPRESS(L)
            TITLEL(L) = TITLEX//XLB
          END DO
        ELSE    !  simple b-grid cases
          DO L=1,LM
            DO J=2,JM
              NI = 0
              FLAT = 0.
              DO I=1,IM
                DP=AIJK(I,J,L,IJK_DP)
                IF(DP.GT.0.) THEN
                  SMAP(I,J,L)=SCALE_IJK(K)*AIJK(I,J,L,K)/DP+OFF_IJK(K)
                  FLAT = FLAT+SMAP(I,J,L)
                  NI = NI+1
                END IF
              END DO
              IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
            END DO
            WRITE(TITLEX(23:30),'(A)') CPRESS(L)
            TITLEL(L) = TITLEX//XLB
          END DO
        END IF
        CALL POUT_IJK(TITLEL,name_ijk(k),lname_ijk(k),units_ijk(k)
     *       ,SMAP,SMAPJK,SMAPK,jgrid_ijk(k))
      END DO
C****
      call close_ijk
C****
C**** ijl output
C****
      call open_ijl(trim(acc_period)//'.ijl'//XLABEL(1:LRUNID),im,jm,lm)

      k=ijl_cf
      TITLEX = lname_ijk(k)(1:17)//"   at  Level    ("//
     *     trim(units_ijk(k))//")"
      SMAP(:,:,:) = UNDEF
      SMAPJK(:,:) = UNDEF
      SMAPK(:)    = UNDEF
      DO L=1,LM
        DO J=1,JM
          NI = 0
          FLAT = 0.
          DO I=1,IM
            SMAP(I,J,L)=SCALE_IJK(K)*AIJK(I,J,L,K)/IDACC(ia_rad)
            FLAT = FLAT+SMAP(I,J,L)
            NI = NI+1
          END DO
          IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
        END DO
        WRITE(TITLEX(31:33),'(I3)') L
        TITLEL(L) = TITLEX//XLB
      END DO
      CALL POUT_IJL(TITLEL,name_ijk(k),lname_ijk(k),units_ijk(k)
     *     ,SMAP,SMAPJK,SMAPK,jgrid_ijk(k))
C*** Begin *** 3-D latent heating diags ***
      if (lh_diags.eq.1) then
       k=ijl_llh
       TITLEX = lname_ijk(k)(1:35)//" at (        "//
     *     trim(units_ijk(k))//")"
       SMAP(:,:,:) = UNDEF
       SMAPJK(:,:) = UNDEF
       SMAPK(:)    = UNDEF
       DO L=1,LM
         DO J=1,JM
          NI = 0
          FLAT = 0.
          DO I=1,IM
            SMAP(I,J,L)=SCALE_IJK(K)*AIJK(I,J,L,K)/IDACC(ia_src)
            FLAT = FLAT+SMAP(I,J,L)
            NI = NI+1
          END DO
          IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
         END DO
         WRITE(TITLEX(41:47),'(A)') CPRESS(L)
         TITLEL(L) = TITLEX//XLB
       END DO
       CALL POUT_IJL(TITLEL,name_ijk(k),lname_ijk(k),units_ijk(k)
     *     ,SMAP,SMAPJK,SMAPK,jgrid_ijk(k))

       k=ijl_mctlh
       TITLEX = lname_ijk(k)(1:30)//" at (        "//
     *     trim(units_ijk(k))//")"
       SMAP(:,:,:) = UNDEF
       SMAPJK(:,:) = UNDEF
       SMAPK(:)    = UNDEF
       DO L=1,LM
         DO J=1,JM
          NI = 0
          FLAT = 0.
          DO I=1,IM
            SMAP(I,J,L)=SCALE_IJK(K)*AIJK(I,J,L,K)/IDACC(ia_src)
            FLAT = FLAT+SMAP(I,J,L)
            NI = NI+1
          END DO
          IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
         END DO
         WRITE(TITLEX(36:42),'(A)') CPRESS(L)
         TITLEL(L) = TITLEX//XLB
       END DO
       CALL POUT_IJL(TITLEL,name_ijk(k),lname_ijk(k),units_ijk(k)
     *     ,SMAP,SMAPJK,SMAPK,jgrid_ijk(k))

       k=ijl_mcdlh
       TITLEX = lname_ijk(k)(1:30)//" at (        "//
     *     trim(units_ijk(k))//")"
       SMAP(:,:,:) = UNDEF
       SMAPJK(:,:) = UNDEF
       SMAPK(:)    = UNDEF
       DO L=1,LM
         DO J=1,JM
          NI = 0
          FLAT = 0.
          DO I=1,IM
            SMAP(I,J,L)=SCALE_IJK(K)*AIJK(I,J,L,K)/IDACC(ia_src)
            FLAT = FLAT+SMAP(I,J,L)
            NI = NI+1
          END DO
          IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
         END DO
         WRITE(TITLEX(36:42),'(A)') CPRESS(L)
         TITLEL(L) = TITLEX//XLB
       END DO
       CALL POUT_IJL(TITLEL,name_ijk(k),lname_ijk(k),units_ijk(k)
     *     ,SMAP,SMAPJK,SMAPK,jgrid_ijk(k))

       k=ijl_mcslh
       TITLEX = lname_ijk(k)(1:30)//" at (        "//
     *     trim(units_ijk(k))//")"
       SMAP(:,:,:) = UNDEF
       SMAPJK(:,:) = UNDEF
       SMAPK(:)    = UNDEF
       DO L=1,LM
         DO J=1,JM
          NI = 0
          FLAT = 0.
          DO I=1,IM
            SMAP(I,J,L)=SCALE_IJK(K)*AIJK(I,J,L,K)/IDACC(ia_src)
            FLAT = FLAT+SMAP(I,J,L)
            NI = NI+1
          END DO
          IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
         END DO
         WRITE(TITLEX(36:42),'(A)') CPRESS(L)
         TITLEL(L) = TITLEX//XLB
       END DO
       CALL POUT_IJL(TITLEL,name_ijk(k),lname_ijk(k),units_ijk(k)
     *     ,SMAP,SMAPJK,SMAPK,jgrid_ijk(k))

       endif
C*** End 3-D latent heating ***

#ifdef CLD_AER_CDNC
C** Cold cloud part
      k=ijl_reim
      TITLEX = lname_ijk(k)(1:17)//"   at  Level    ("//
     *     trim(units_ijk(k))//")"
      SMAP(:,:,:) = UNDEF
      SMAPJK(:,:) = UNDEF
      SMAPK(:)    = UNDEF
      DO L=1,LM
        DO J=1,JM
          NI = 0
          FLAT = 0.
          DO I=1,IM
            SMAP(I,J,L)=SCALE_IJK(K)*AIJK(I,J,L,K)/IDACC(ia_src)
            FLAT = FLAT+SMAP(I,J,L)
            NI = NI+1
          END DO
          IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
        END DO
        WRITE(TITLEX(31:33),'(I3)') L
        TITLEL(L) = TITLEX//XLB
      END DO
      CALL POUT_IJL(TITLEL,name_ijk(k),lname_ijk(k),units_ijk(k)
     *     ,SMAP,SMAPJK,SMAPK,jgrid_ijk(k))

      k=ijl_reis
      TITLEX = lname_ijk(k)(1:17)//"   at  Level    ("//
     *     trim(units_ijk(k))//")"
      SMAP(:,:,:) = UNDEF
      SMAPJK(:,:) = UNDEF
      SMAPK(:)    = UNDEF
      DO L=1,LM
        DO J=1,JM
          NI = 0
          FLAT = 0.
          DO I=1,IM
            SMAP(I,J,L)=SCALE_IJK(K)*AIJK(I,J,L,K)/IDACC(ia_src)
            FLAT = FLAT+SMAP(I,J,L)
            NI = NI+1
          END DO
          IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
        END DO
        WRITE(TITLEX(31:33),'(I3)') L
        TITLEL(L) = TITLEX//XLB
      END DO
      CALL POUT_IJL(TITLEL,name_ijk(k),lname_ijk(k),units_ijk(k)
     *     ,SMAP,SMAPJK,SMAPK,jgrid_ijk(k))

      k=ijl_cdis
      TITLEX = lname_ijk(k)(1:17)//"   at  Level    ("//
     *     trim(units_ijk(k))//")"
      SMAP(:,:,:) = UNDEF
      SMAPJK(:,:) = UNDEF
      SMAPK(:)    = UNDEF
      DO L=1,LM
        DO J=1,JM
          NI = 0
          FLAT = 0.
          DO I=1,IM
            SMAP(I,J,L)=SCALE_IJK(K)*AIJK(I,J,L,K)/IDACC(ia_src)
            FLAT = FLAT+SMAP(I,J,L)
            NI = NI+1
          END DO
          IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
        END DO
        WRITE(TITLEX(31:33),'(I3)') L
        TITLEL(L) = TITLEX//XLB
      END DO
      CALL POUT_IJL(TITLEL,name_ijk(k),lname_ijk(k),units_ijk(k)
     *     ,SMAP,SMAPJK,SMAPK,jgrid_ijk(k))

      k=ijl_cdim
      TITLEX = lname_ijk(k)(1:17)//"   at  Level    ("//
     *     trim(units_ijk(k))//")"
      SMAP(:,:,:) = UNDEF
      SMAPJK(:,:) = UNDEF
      SMAPK(:)    = UNDEF
      DO L=1,LM
        DO J=1,JM
          NI = 0
          FLAT = 0.
          DO I=1,IM
            SMAP(I,J,L)=SCALE_IJK(K)*AIJK(I,J,L,K)/IDACC(ia_src)
            FLAT = FLAT+SMAP(I,J,L)
            NI = NI+1
          END DO
          IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
        END DO
        WRITE(TITLEX(31:33),'(I3)') L
        TITLEL(L) = TITLEX//XLB
      END DO
      CALL POUT_IJL(TITLEL,name_ijk(k),lname_ijk(k),units_ijk(k)
     *     ,SMAP,SMAPJK,SMAPK,jgrid_ijk(k))

      k=ijl_cwim
      TITLEX = lname_ijk(k)(1:17)//"   at  Level    ("//
     *     trim(units_ijk(k))//")"
      SMAP(:,:,:) = UNDEF
      SMAPJK(:,:) = UNDEF
      SMAPK(:)    = UNDEF
      DO L=1,LM
        DO J=1,JM
          NI = 0
          FLAT = 0.
          DO I=1,IM
            SMAP(I,J,L)=SCALE_IJK(K)*AIJK(I,J,L,K)/IDACC(ia_src)
            FLAT = FLAT+SMAP(I,J,L)
            NI = NI+1
          END DO
          IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
        END DO
        WRITE(TITLEX(31:33),'(I3)') L
        TITLEL(L) = TITLEX//XLB
      END DO
      CALL POUT_IJL(TITLEL,name_ijk(k),lname_ijk(k),units_ijk(k)
     *     ,SMAP,SMAPJK,SMAPK,jgrid_ijk(k))

      k=ijl_cwis
      TITLEX = lname_ijk(k)(1:17)//"   at  Level    ("//
     *     trim(units_ijk(k))//")"
      SMAP(:,:,:) = UNDEF
      SMAPJK(:,:) = UNDEF
      SMAPK(:)    = UNDEF
      DO L=1,LM
        DO J=1,JM
          NI = 0
          FLAT = 0.
          DO I=1,IM
            SMAP(I,J,L)=SCALE_IJK(K)*AIJK(I,J,L,K)/IDACC(ia_src)
            FLAT = FLAT+SMAP(I,J,L)
            NI = NI+1
          END DO
          IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
        END DO
        WRITE(TITLEX(31:33),'(I3)') L
        TITLEL(L) = TITLEX//XLB
      END DO
      CALL POUT_IJL(TITLEL,name_ijk(k),lname_ijk(k),units_ijk(k)
     *     ,SMAP,SMAPJK,SMAPK,jgrid_ijk(k))

C*** Warm cloud part
      k=ijl_rewm
      TITLEX = lname_ijk(k)(1:17)//"   at  Level    ("//
     *     trim(units_ijk(k))//")"
      SMAP(:,:,:) = UNDEF
      SMAPJK(:,:) = UNDEF
      SMAPK(:)    = UNDEF
      DO L=1,LM
        DO J=1,JM
          NI = 0
          FLAT = 0.
          DO I=1,IM
            SMAP(I,J,L)=SCALE_IJK(K)*AIJK(I,J,L,K)/IDACC(ia_src)
c      if (AIJK(I,J,L,K).gt.5.d0)
c    * write(6,*)"Reff",AIJK(I,J,L,K),I,J,L
            FLAT = FLAT+SMAP(I,J,L)
            NI = NI+1
          END DO
          IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
        END DO
        WRITE(TITLEX(31:33),'(I3)') L
        TITLEL(L) = TITLEX//XLB
      END DO
      CALL POUT_IJL(TITLEL,name_ijk(k),lname_ijk(k),units_ijk(k)
     *     ,SMAP,SMAPJK,SMAPK,jgrid_ijk(k))

      k=ijl_rews
      TITLEX = lname_ijk(k)(1:17)//"   at  Level    ("//
     *     trim(units_ijk(k))//")"
      SMAP(:,:,:) = UNDEF
      SMAPJK(:,:) = UNDEF
      SMAPK(:)    = UNDEF
      DO L=1,LM
        DO J=1,JM
          NI = 0
          FLAT = 0.
          DO I=1,IM
            SMAP(I,J,L)=SCALE_IJK(K)*AIJK(I,J,L,K)/IDACC(ia_src)
            FLAT = FLAT+SMAP(I,J,L)
            NI = NI+1
          END DO
          IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
        END DO
        WRITE(TITLEX(31:33),'(I3)') L
        TITLEL(L) = TITLEX//XLB
      END DO
      CALL POUT_IJL(TITLEL,name_ijk(k),lname_ijk(k),units_ijk(k)
     *     ,SMAP,SMAPJK,SMAPK,jgrid_ijk(k))

      k=ijl_cdws
      TITLEX = lname_ijk(k)(1:17)//"   at  Level    ("//
     *     trim(units_ijk(k))//")"
      SMAP(:,:,:) = UNDEF
      SMAPJK(:,:) = UNDEF
      SMAPK(:)    = UNDEF
      DO L=1,LM
        DO J=1,JM
          NI = 0
          FLAT = 0.
          DO I=1,IM
            SMAP(I,J,L)=SCALE_IJK(K)*AIJK(I,J,L,K)/IDACC(ia_src)
            FLAT = FLAT+SMAP(I,J,L)
            NI = NI+1
          END DO
          IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
        END DO
        WRITE(TITLEX(31:33),'(I3)') L
        TITLEL(L) = TITLEX//XLB
      END DO
      CALL POUT_IJL(TITLEL,name_ijk(k),lname_ijk(k),units_ijk(k)
     *     ,SMAP,SMAPJK,SMAPK,jgrid_ijk(k))

      k=ijl_cdwm
      TITLEX = lname_ijk(k)(1:17)//"   at  Level    ("//
     *     trim(units_ijk(k))//")"
      SMAP(:,:,:) = UNDEF
      SMAPJK(:,:) = UNDEF
      SMAPK(:)    = UNDEF
      DO L=1,LM
        DO J=1,JM
          NI = 0
          FLAT = 0.
          DO I=1,IM
            SMAP(I,J,L)=SCALE_IJK(K)*AIJK(I,J,L,K)/IDACC(ia_src)
            FLAT = FLAT+SMAP(I,J,L)
            NI = NI+1
          END DO
          IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
        END DO
        WRITE(TITLEX(31:33),'(I3)') L
        TITLEL(L) = TITLEX//XLB
      END DO
      CALL POUT_IJL(TITLEL,name_ijk(k),lname_ijk(k),units_ijk(k)
     *     ,SMAP,SMAPJK,SMAPK,jgrid_ijk(k))

      k=ijl_cwwm
      TITLEX = lname_ijk(k)(1:17)//"   at  Level    ("//
     *     trim(units_ijk(k))//")"
      SMAP(:,:,:) = UNDEF
      SMAPJK(:,:) = UNDEF
      SMAPK(:)    = UNDEF
      DO L=1,LM
        DO J=1,JM
          NI = 0
          FLAT = 0.
          DO I=1,IM
            SMAP(I,J,L)=SCALE_IJK(K)*AIJK(I,J,L,K)/IDACC(ia_src)
            FLAT = FLAT+SMAP(I,J,L)
            NI = NI+1
          END DO
          IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
        END DO
        WRITE(TITLEX(31:33),'(I3)') L
        TITLEL(L) = TITLEX//XLB
      END DO
      CALL POUT_IJL(TITLEL,name_ijk(k),lname_ijk(k),units_ijk(k)
     *     ,SMAP,SMAPJK,SMAPK,jgrid_ijk(k))

      k=ijl_cwws
      TITLEX = lname_ijk(k)(1:17)//"   at  Level    ("//
     *     trim(units_ijk(k))//")"
      SMAP(:,:,:) = UNDEF
      SMAPJK(:,:) = UNDEF
      SMAPK(:)    = UNDEF
      DO L=1,LM
        DO J=1,JM
          NI = 0
          FLAT = 0.
          DO I=1,IM
            SMAP(I,J,L)=SCALE_IJK(K)*AIJK(I,J,L,K)/IDACC(ia_src)
            FLAT = FLAT+SMAP(I,J,L)
            NI = NI+1
          END DO
          IF (NI.GT.0) SMAPJK(J,L) = FLAT/NI
        END DO
        WRITE(TITLEX(31:33),'(I3)') L
        TITLEL(L) = TITLEX//XLB
      END DO
      CALL POUT_IJL(TITLEL,name_ijk(k),lname_ijk(k),units_ijk(k)
     *     ,SMAP,SMAPJK,SMAPK,jgrid_ijk(k))
#endif
C****
      call close_ijl
C****
      RETURN
      END SUBROUTINE IJKMAP

      function NINTlimit( x )
      real*8 x
      integer NINTlimit
      real*8 y
      y = min (  2147483647.d0, x )
      y = max ( -2147483647.d0, y )
      NINTlimit = NINT( y )
      return
      end function NINTlimit

      subroutine diag_isccp
c     routine emptied - add content here
      end subroutine diag_isccp

      subroutine diag_msu
!@sum diag_msu computes MSU channel 2,3,4 temperatures as weighted means
!@auth Reto A Ruedy (input file created by Makiko Sato)
      use filemanager
      USE CONSTANT
      USE DIAG_COM
      USE MODEL_COM
      USE BDIJ

      implicit none

      integer, parameter :: nmsu=200 , ncols=4
      real*8 plbmsu(nmsu),wmsu(ncols,nmsu) ;  save plbmsu,wmsu
      real*8 tlmsu(nmsu),tmsu(ncols,im,jm)
      real*8 plb(0:lm+2),tlb(0:lm+2),tlm(lm)

      integer i,j,l,n,ip1,iu_msu  ;           integer, save :: ifirst=1
      real*8  ts,pland,dp

      if(ifirst.eq.1) then
c**** read in the weights file
        call openunit('MSU_wts',iu_msu,.false.,.true.)
        do n=1,4
          read(iu_msu,*)
        end do

        do l=1,nmsu
          read(iu_msu,*) plbmsu(l),(wmsu(n,l),n=1,ncols)
        end do

        call closeunit(iu_msu)
        ifirst=0
      end if

c**** Collect temperatures and pressures (on the secondary grid)
      do i=2,im
        aij(i,1,ij_ts)=aij(1,1,ij_ts)
        aij(i,jm,ij_ts)=aij(1,jm,ij_ts)
      end do
      do j=2,jm
      i=IM
      do ip1=1,im
        ts=.25*(aij(i,j,ij_ts)+aij(i,j-1,ij_ts)+
     +      aij(ip1,j,ij_ts)+aij(ip1,j-1,ij_ts))/idacc(ia_src)
        pland = .25*(wt_ij(i,j,iw_land) + wt_ij(i,j-1,iw_land)
     &       + wt_ij(ip1,j,iw_land) + wt_ij(ip1,j-1,iw_land))
        plb(lm+1)=pmtop
        do l=lm,1,-1
          dp=aijk(i,j,l,ijk_dp)
          plb(l)=plb(l+1)+dp/idacc(ia_dga)
          tlm(l)=ts
          if(dp.gt.0.) tlm(l)=aijk(i,j,l,ijk_t)/dp - tf
        end do
c**** find edge temperatures (assume continuity and given means)
        tlb(0)=ts ; plb(0)=plbmsu(1) ; tlb(1)=ts
        do l=1,lm
           tlb(l+1)=2*tlm(l)-tlb(l)
        end do
        tlb(lm+2)=tlb(lm+1) ; plb(lm+2)=0.
        call vntrp1 (lm+2,plb,tlb, nmsu-1,plbmsu,tlmsu)
c**** find MSU channel 2,3,4 temperatures
        tmsu(:,i,j)=0.
        do l=1,nmsu-1
          tmsu(:,i,j)=tmsu(:,i,j)+tlmsu(l)*wmsu(:,l)
        end do
        tmsu2(i,j)=(1-pland)*tmsu(1,i,j)+pland*tmsu(2,i,j)
        tmsu3(i,j)=tmsu(3,i,j)
        tmsu4(i,j)=tmsu(4,i,j)
        i=ip1
      end do
      end do

      return
      end subroutine diag_msu

      SUBROUTINE VNTRP1 (KM,P,AIN,  LMA,PE,AOUT)
C**** Vertically interpolates a 1-D array
C**** Input:       KM = number of input pressure levels
C****            P(K) = input pressure levels (mb)
C****          AIN(K) = input quantity at level P(K)
C****             LMA = number of vertical layers of output grid
C****           PE(L) = output pressure levels (mb) (edges of layers)
C**** Output: AOUT(L) = output quantity: mean between PE(L-1) & PE(L)
C****
      implicit none
      integer, intent(in) :: km,lma
      REAL*8, intent(in)  :: P(0:KM),AIN(0:KM),    PE(0:LMA)
      REAL*8, intent(out) :: AOUT(LMA)

      integer k,k1,l
      real*8 pdn,adn,pup,aup,psum,asum

C****
      PDN = PE(0)
      ADN = AIN(0)
      K=1
C**** Ignore input levels below ground level pe(0)=p(0)
      IF(P(1).GT.PE(0)) THEN
         DO K1=2,KM
         K=K1
         IF(P(K).LT.PE(0)) THEN  ! interpolate to ground level
           ADN=AIN(K)+(AIN(K-1)-AIN(K))*(PDN-P(K))/(P(K-1)-P(K))
           GO TO 300
         END IF
         END DO
         STOP 'VNTRP1 - error - should not get here'
      END IF
C**** Integrate - connecting input data by straight lines
  300 DO 330 L=1,LMA
      ASUM = 0.
      PSUM = 0.
      PUP = PE(L)
  310 IF(P(K).le.PUP)  GO TO 320
      PSUM = PSUM + (PDN-P(K))
      ASUM = ASUM + (PDN-P(K))*(ADN+AIN(K))/2.
      PDN  = P(K)
      ADN  = AIN(K)
      K=K+1
      IF(K.LE.KM) GO TO 310
      stop 'VNTRP1 - should not happen'
C****
  320 AUP  = AIN(K) + (ADN-AIN(K))*(PUP-P(K))/(PDN-P(K))
      PSUM = PSUM + (PDN-PUP)
      ASUM = ASUM + (PDN-PUP)*(ADN+AUP)/2.
      AOUT(L) = ASUM/PSUM
      PDN = PUP
  330 ADN = AUP
C****
      RETURN
      END subroutine vntrp1

      SUBROUTINE DIAG_GATHER
      USE MODEL_COM, only : IM, FOCEAN, FLICE, ZATMO
cddd      USE LAKES_COM, only : FLAKE
cddd      USE GHY_COM, only : FEARTH
cddd#ifdef USE_ENT
cddd      use ent_com, only : entcells
cddd      use ent_mod, only : ent_get_exports
cddd#else
cddd      USE VEG_COM,   only : vdata
cddd#endif
      USE DIAG_COM, only : AIJ,  AIJ_loc, AJ,   AJ_loc, AREGJ,
     *     AREGJ_loc, AJK,  AJK_loc, AIJK, AIJK_loc,
     *     ASJL, ASJL_loc, AJL,  AJL_loc , CONSRV, CONSRV_loc, TSFREZ,
     *     TSFREZ_loc, WT_IJ
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, PACK_DATA, PACK_DATAj !, GET
      USE DOMAIN_DECOMP_ATM, ONLY : CHECKSUMj,CHECKSUM,am_i_root
      USE CONSTANT, only : NaN
      IMPLICIT NONE
cddd      INTEGER :: J_0, J_1, J_0H, J_1H
cddd      REAL*8, ALLOCATABLE :: tmp(:,:)
cddd#ifdef USE_ENT
cddd      REAL*8, ALLOCATABLE :: fract_vege(:,:)
cddd      INTEGER i,j
cddd#endif

      if(am_i_root()) call alloc_diag_com_glob

      call Gather_Diagnostics()
c      CALL PACK_DATAj(GRID, AJ_loc,  AJ)
c      CALL PACK_DATA(GRID, AREGJ_loc,  AREGJ)
c      CALL PACK_DATAj(GRID, APJ_loc, APJ)
c      CALL PACK_DATAj(GRID, AJK_loc, AJK)
c      CALL PACK_DATA (GRID, AIJ_loc, AIJ)
c      CALL PACK_DATA (GRID, AIJK_loc, AIJK)
c      CALL PACK_DATAj(GRID, ASJL_loc, ASJL)
c      CALL PACK_DATAj(GRID, AJL_loc,  AJL)
c      CALL PACK_DATAj(GRID, CONSRV_loc,  CONSRV)
c      CALL PACK_DATA (GRID, TSFREZ_loc,  TSFREZ)

#ifdef TRACERS_ON
      call gather_trdiag
#endif

! Now the external arrays
      !!CALL PACK_DATA(GRID, fland, fland_glob)
      !!CALL PACK_DATA(GRID, fearth, fearth_glob)
      CALL PACK_DATA(GRID, focean, focean_glob)
      CALL PACK_DATA(GRID, flice, flice_glob)
      CALL PACK_DATA(GRID, zatmo, zatmo_glob)

cddd      CALL GET(GRID, J_STRT=J_0, J_STOP=J_1,
cddd     &     J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
cddd      ALLOCATE(tmp(IM, J_0H:J_1H))
cddd
cddd      wt_ij(:,:,1) = 1.
      wt_ij(:,:,:) = 1.  ! NaN
      !!CALL PACK_DATA(GRID, focean, wt_ij(:,:,2))
      !!CALL PACK_DATA(GRID, flake,  wt_ij(:,:,3))  ! not correct
      !!CALL PACK_DATA(GRID, flice,  wt_ij(:,:,4))
      !!CALL PACK_DATA(GRID, fearth, wt_ij(:,:,5))  ! not correct

cddd#ifdef USE_ENT
cddd      ALLOCATE(fract_vege(IM, J_0H:J_1H))
cddd      call ent_get_exports( entcells(1:IM,J_0:J_1),
cddd     &           fraction_of_vegetated_soil=fract_vege(1:IM,J_0:J_1) )
cddd      tmp(:,J_0:J_1) = fearth(:,J_0:J_1) * (1.d0-fract_vege(:,J_0:J_1))
cddd      CALL PACK_DATA(GRID, tmp, wt_ij(:,:,6))
cddd      tmp(:,J_0:J_1) = fearth(:,J_0:J_1) * fract_vege(:,J_0:J_1)
cddd      CALL PACK_DATA(GRID, tmp, wt_ij(:,:,7))
cddd      DEALLOCATE(fract_vege)
cddd#else
cddd      tmp(:,J_0:J_1) = fearth(:,J_0:J_1) *
cddd     &     (vdata(:,J_0:J_1,1)+vdata(:,J_0:J_1,10))
cddd      CALL PACK_DATA(GRID, tmp, wt_ij(:,:,6))
cddd      tmp(:,J_0:J_1) = fearth(:,J_0:J_1) *
cddd     &     (1.-(vdata(:,J_0:J_1,1)+vdata(:,J_0:J_1,10)))
cddd      CALL PACK_DATA(GRID, tmp, wt_ij(:,:,7))
cddd#endif
cddd      DEALLOCATE(tmp)

      call gather_odiags ()

      END SUBROUTINE DIAG_GATHER

      SUBROUTINE DIAG_SCATTER
      USE DIAG_COM, only : AIJ, AIJ_loc, AJ,  AJ_loc, AREGJ, AREGJ_loc,
     *     AJK, AJK_loc, AIJK, AIJK_loc, ASJL, ASJL_loc,
     *     AJL,  AJL_loc, CONSRV, CONSRV_loc, TSFREZ, TSFREZ_loc
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, UNPACK_DATA, UNPACK_DATAj
      USE DOMAIN_DECOMP_ATM, ONLY : am_i_root
      IMPLICIT NONE

#ifndef CUBE_GRID
      CALL UNPACK_DATAj(GRID, AJ,  AJ_loc)
      CALL UNPACK_DATA(GRID, AREGJ,  AREGJ_loc)
#endif
      CALL UNPACK_DATAj(GRID, AJK, AJK_loc)
      CALL UNPACK_DATA (GRID, AIJ, AIJ_loc)
      CALL UNPACK_DATA (GRID, AIJK, AIJK_loc)
      CALL UNPACK_DATAj(GRID, ASJL, ASJL_loc)
      CALL UNPACK_DATAj(GRID, AJL,  AJL_loc)
      CALL UNPACK_DATAj(GRID, CONSRV,  CONSRV_loc)
      CALL UNPACK_DATA (GRID, TSFREZ,  TSFREZ_loc)

#ifdef TRACERS_ON
      call scatter_trdiag
#endif
      if(am_i_root()) call dealloc_diag_com_glob

      END SUBROUTINE DIAG_SCATTER

      subroutine certain_wind_averages
      use diag_com, only : im,jm,lm,ail,ajl
     &     ,jl_uepac,jl_vepac,jl_wepac,jl_uwpac,jl_vwpac,jl_wwpac
     &     ,il_u,il_v,il_w
      implicit none
      INTEGER, PARAMETER ::
     *     I150E = IM*(180+150)/360+1, ! WEST EDGE OF 150 EAST
     *     I110W = IM*(180-110)/360+1, ! WEST EDGE OF 110 WEST
     *     I135W = IM*(180-135)/360+1  ! WEST EDGE OF 135 WEST
      integer :: i,j,l
C****
C**** CERTAIN HORIZONTAL WIND AVERAGES
C****
      DO L=1,LM
      DO J=1,JM
        ajl(j,l,jl_uepac) = 0d0
        ajl(j,l,jl_vepac) = 0d0
        DO I=I135W,I110W        ! EAST PACIFIC
          AJL(J,L,JL_UEPAC)=AJL(J,L,JL_UEPAC)+AIL(I,J,L,IL_U)
          AJL(J,L,JL_VEPAC)=AJL(J,L,JL_VEPAC)+AIL(I,J,L,IL_V)
        END DO
        ajl(j,l,jl_uepac) = ajl(j,l,jl_uepac)/(1+I110W-I135W)
        ajl(j,l,jl_vepac) = ajl(j,l,jl_vepac)/(1+I110W-I135W)
        ajl(j,l,jl_uwpac) = 0d0
        ajl(j,l,jl_vwpac) = 0d0
        DO I=I150E,IM           ! WEST PACIFIC
          AJL(J,L,JL_UWPAC)=AJL(J,L,JL_UWPAC)+AIL(I,J,L,IL_U)
          AJL(J,L,JL_VWPAC)=AJL(J,L,JL_VWPAC)+AIL(I,J,L,IL_V)
        END DO
        ajl(j,l,jl_uwpac) = ajl(j,l,jl_uwpac)/(1+IM-I150E)
        ajl(j,l,jl_vwpac) = ajl(j,l,jl_vwpac)/(1+IM-I150E)
      END DO
      END DO
C****
C**** CERTAIN VERTICAL WIND AVERAGES
C****
      DO L=1,LM-1
      DO J=1,JM
        ajl(j,l,jl_wepac) = 0d0
        DO I=I135W,I110W        ! EAST PACIFIC
          AJL(J,L,JL_WEPAC)=AJL(J,L,JL_WEPAC)+AIL(I,J,L,IL_W)
        END DO
        ajl(j,l,jl_wepac) = ajl(j,l,jl_wepac)/(1+I110W-I135W)
        ajl(j,l,jl_wwpac) = 0d0
        DO I=I150E,IM           ! WEST PACIFIC
          AJL(J,L,JL_WWPAC)=AJL(J,L,JL_WWPAC)+AIL(I,J,L,IL_W)
        END DO
        ajl(j,l,jl_wwpac) = ajl(j,l,jl_wwpac)/(1+IM-I150E)
      END DO
      END DO
      return
      end subroutine certain_wind_averages

      END MODULE DIAG_SERIAL
