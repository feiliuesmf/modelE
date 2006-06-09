#include "rundeck_opts.h"
      MODULE WORKJK
      USE MODEL_COM, ONLY : JM,LM
      REAL*8, DIMENSION(JM,LM,2):: DPJK
      REAL*8, DIMENSION(2,LM,2) :: DPHEM
      REAL*8, DIMENSION(LM,2) :: DPGLOB
      END MODULE WORKJK
      BLOCK DATA BDWP
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
      MODULE BDJ
      IMPLICIT NONE
      SAVE
      integer, parameter :: nj_out=11
      CHARACTER(LEN=50), DIMENSION(nj_out) :: UNITS_J_O
      CHARACTER(LEN=50), DIMENSION(nj_out) :: LNAME_J_O
      CHARACTER(LEN=30), DIMENSION(nj_out) :: NAME_J_O
      CHARACTER(LEN=16), DIMENSION(nj_out) :: STITLE_J_O
      INTEGER, DIMENSION(nj_out) :: INUM_J_O, IDEN_J_O
      REAL*8, DIMENSION(nj_out) :: SCALE_J_O
      END MODULE BDJ
      MODULE BDjkjl
      IMPLICIT NONE
      SAVE
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
      MODULE BDIJ
      use MODEL_COM, only : IM,JM
      use DIAG_COM
      IMPLICIT NONE
      SAVE
      integer nij_o
      INTEGER :: ij_topo, ij_jet, ij_wsmn, ij_jetdir, ij_wsdir, ij_grow,
     *  ij_netrdp, ij_albp, ij_albg, ij_albv, ij_ntdsese, ij_ntdsete,
     *  ij_fland, ij_dzt1, ij_albgv, ij_colh2o !,ij_msu2,ij_msu3,ij_msu4
      REAL*8, DIMENSION(IM,JM) :: SENTDSE,TENTDSE, TMSU2,TMSU3,TMSU4
      contains
      function mark (val,ibar,undef)
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
      module DIAGKS
      USE CONSTANT, only :
     &     twopi
      USE MODEL_COM, only :
     &     jm,lm,jeq, JHOUR,JHOUR0,
     &     JDATE,JDATE0,JMON,JMON0,AMON,AMON0,JYEAR,JYEAR0,
     &     Itime,ItimeI,Itime0,PTOP,SIG,XLABEL,
     &     PSFMPT,AMONTH,nday
      USE GEOM, only :
     &     DLAT,DXYP,LAT_DG
      USE DIAG_COM, only :
     &     keyct,keynr,ned,nkeynr
      USE PARAM
      IMPLICIT NONE
      PRIVATE
      SAVE
      public KEYDJ,KEYDJA,KEYJKT,KEYJKJ,KEYJLS,KEYJKE,KEYJKN,KEYIJ
     &     ,KEYD4,DIAGKN
      INTEGER*4 :: NDEX(42) = (/
     *     1,2,31,32,3,   4,5,6,7,8,   9,10,11,12,13,
     *             14,15,16,17,18,  19,20,21,22,23,  24,25,26,27,28,
     *             29,30,33,34,35,  36,37,38,39,40,  41,42/)
      REAL*8, DIMENSION(JM,LM) :: FKEY
      INTEGER ::
     &     I,I35,I70,IEND,ISIGN,
     &     J,J60,JMAX,JNDEX,JSTART,
     &     K,KEYMAX,KNDEX,
     &     LL,LMAX,LNLM,LNM,LSLM,LSM,
     &     NTDIF
      REAL*8 ::
     &     A,BIG,CEPT,CHECK,
     &     HN,HS,
     &     SAVE,DAYS,TEQ,TNOR,TSOU,
     &     UNLM,UNM,USLM,USM,X60
      contains
      subroutine KEYDJ(N,FGLOB,FNH)
      integer N
      real*8 FGLOB,FNH
      end subroutine KEYDJ
      subroutine KEYDJA (FGLOB)
      real*8 FGLOB
      end subroutine KEYDJA
      subroutine KEYJKT (GSUM,ASUM)
      real*8 GSUM
      REAL*8, DIMENSION(JM) :: ASUM
      end subroutine KEYJKT
      subroutine KEYJKJ (L,FLAT)
      integer L
      REAL*8, DIMENSION(JM) :: FLAT
      end subroutine KEYJKJ
      subroutine KEYJLS (L,FLAT)
      integer L
      REAL*8, DIMENSION(JM) :: FLAT
      end subroutine KEYJLS
      subroutine KEYJKE (NT,HSUM,ASUM)
      integer NT
      REAL*8, DIMENSION(2) :: HSUM
      REAL*8, DIMENSION(JM) :: ASUM
      end subroutine KEYJKE
      subroutine KEYJKN (NT,ASUM,SUMFAC)
      integer NT
      REAL*8, DIMENSION(JM) :: ASUM
      REAL*8 SUMFAC
      end subroutine KEYJKN
      subroutine KEYIJ(PISG,PISN)
      REAL*8 PISG,PISN
      end subroutine KEYIJ
      subroutine KEYD4 (IK)
      INTEGER, DIMENSION(2*NED) :: IK
      end subroutine KEYD4
      subroutine DIAGKN
      end subroutine DIAGKN
      END MODULE DIAGKS
      MODULE DIAG_SERIAL
      USE MODEL_COM, ONLY : IM, JM
      USE DOMAIN_DECOMP, only : grid, DIST_GRID, AM_I_ROOT
      USE DIAGKS
      PRIVATE
      PUBLIC :: PRINT_DIAGS
      PUBLIC :: JLMAP
      PUBLIC :: MAPTXT
      PUBLIC :: IJ_avg
      INTERFACE GLOBALSUM
        MODULE PROCEDURE GLOBALSUM_J
        MODULE PROCEDURE GLOBALSUM_JK
      END INTERFACE
      REAL*8 :: FLAND_glob(IM,JM)
      REAL*8 :: FEARTH_glob(IM,JM)
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
      END SUBROUTINE GLOBALSUM_J
      SUBROUTINE GLOBALSUM_JK(grd_dum, garr, gsum, hsum, istag, all)
      IMPLICIT NONE
      TYPE (DIST_GRID),  INTENT(IN) :: grd_dum
      REAL*8,            INTENT(IN) :: garr(              : ,:)
      REAL*8,            INTENT(OUT):: gsum(size(garr,2))
      REAL*8, OPTIONAL,  INTENT(OUT):: hsum(2,size(garr,2))
      INTEGER,OPTIONAL,  INTENT(IN) :: istag
      LOGICAL,OPTIONAL,  INTENT(IN) :: all
      INTEGER :: k
      INTEGER :: ier
      INTEGER :: IM, JM
      LOGICAL :: istag_
      END SUBROUTINE GLOBALSUM_JK
      subroutine print_diags(partial)
      USE MODEL_COM, only : itime,itimeI
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: partial
      end subroutine print_diags
      SUBROUTINE J_TITLES
      END SUBROUTINE J_TITLES
      SUBROUTINE DIAGJ
      END SUBROUTINE DIAGJ
      SUBROUTINE JKJL_TITLEX
      END SUBROUTINE JKJL_TITLEX
      SUBROUTINE DIAGJK
      END SUBROUTINE DIAGJK
      SUBROUTINE JKMAP(LNAME,SNAME,UNITS,POW10P,
     &     PM,AX,SCALET,SCALEJ,SCALEK,KMAX,JWT,J1,
     *  ARQX,SCALER,SCALJR,SCALLR)
      USE MODEL_COM, only :
     &     jm,lm,JDATE,JDATE0,JMON0,JMON,AMON0,AMON,JYEAR,JYEAR0,XLABEL
      IMPLICIT NONE
      CHARACTER(LEN=50) :: UNITS,UNITS_WITH_SCALE
      CHARACTER(LEN=50) :: LNAME
      CHARACTER(LEN=30) :: SNAME
      CHARACTER(LEN=64) :: TITLE
      INTEGER, DIMENSION(JM) :: MLAT
      REAL*8, DIMENSION(JM) :: FLAT,ASUM
      REAL*8, DIMENSION(2) :: AHEM,AHEML
      INTEGER :: J1,JWT,KMAX
      REAL*8 :: SCALET,SCALER,PRTFAC
      INTEGER :: POW10P
      REAL*8, DIMENSION(JM,LM) :: AX
      REAL*8, DIMENSION(JM,1) :: ARQX
      REAL*8, DIMENSION(JM) :: SCALEJ,SCALJR
      REAL*8, DIMENSION(LM) :: SCALEK
      REAL*8, DIMENSION(1) :: SCALLR
      REAL*8, DIMENSION(1) :: PM
      REAL*8, DIMENSION(JM,LM) :: CX
      CHARACTER*4 DASH,WORD(4)
      DATA DASH/'----'/,WORD/'SUM','MEAN',' ','.1*'/
      INTEGER :: IWORD,J,JHEMI,K,L ,ksx,klmax
      REAL*8 :: AGLOB,FGLOB,FLATJ,G1,H1,H2,SUMFAC
      CHARACTER XLB*16,CLAT*16,CPRES*16,CBLANK*16,TITLEO*80,TPOW*8
      DATA CLAT/'LATITUDE'/,CPRES/'PRESSURE (MB)'/,CBLANK/' '/
      CHARACTER*200 :: out_line
      END SUBROUTINE JKMAP
      SUBROUTINE JLMAP(LNAME,SNAME,UNITS,POW10P,
     &     PL,AX,SCALET,SCALEJ,SCALEL,LMAX,JWT,J1,
     *  ARQX,SCALER,SCALJR,SCALLR)
      USE MODEL_COM, only :
     &     jm,lm,DSIG,JDATE,JDATE0,AMON,AMON0,JYEAR,JYEAR0,SIGE,XLABEL
      IMPLICIT NONE
      CHARACTER(LEN=50) :: UNITS,UNITS_WITH_SCALE
      CHARACTER(LEN=50) :: LNAME
      CHARACTER(LEN=30) :: SNAME
      CHARACTER(LEN=64) :: TITLE
      REAL*8, DIMENSION(JM) :: FLAT,ASUM
      REAL*8, DIMENSION(2) :: FHEM,HSUM
      INTEGER :: J1,JWT,LMAX
      REAL*8 :: SCALET,SCALER,PRTFAC
      INTEGER :: POW10P
      REAL*8, DIMENSION(JM,LM) :: AX
      REAL*8, DIMENSION(JM,1) :: ARQX
      REAL*8, DIMENSION(JM) :: SCALEJ,SCALJR
      REAL*8, DIMENSION(LM) :: SCALEL
      REAL*8, DIMENSION(1) :: SCALLR
      REAL*8, DIMENSION(:) :: PL
      CHARACTER*4 DASH,WORD(4)
      DATA DASH/'----'/,WORD/'SUM','MEAN',' ','.1*'/
      INTEGER :: IWORD,J,JH,JHEMI,K,L  ,ksx,klmax
      REAL*8 :: FGLOB,GSUM,SDSIG,SUMFAC
      CHARACTER XLB*16,CLAT*16,CPRES*16,CBLANK*16,TITLEO*80,TPOW*8
      DATA CLAT/'LATITUDE'/,CPRES/'PRESSURE (MB)'/,CBLANK/' '/
      END SUBROUTINE JLMAP
      SUBROUTINE JLVMAP(LNAME,SNAME,UNITS,POW10P,
     &     PL,AX,SCALET,SCALEJ,SCALEL,LMAX,JWT,J1,VWT)
      USE MODEL_COM, only :
     &     jm,lm,DSIG,JDATE,JDATE0,AMON,AMON0,JYEAR,JYEAR0,SIGE,XLABEL
      IMPLICIT NONE
      CHARACTER(LEN=50) :: UNITS,UNITS_WITH_SCALE
      CHARACTER(LEN=50) :: LNAME
      CHARACTER(LEN=30) :: SNAME
      CHARACTER(LEN=64) :: TITLE
      REAL*8, DIMENSION(JM) :: FLAT,ASUM,SVWTJ
      REAL*8, DIMENSION(2) :: FHEM,HSUM,HVWT
      REAL*8, DIMENSION(2) :: HM
      REAL*8               :: GM
      INTEGER :: J1,JWT,LMAX
      REAL*8 :: SCALET,PRTFAC
      INTEGER :: POW10P
      REAL*8, DIMENSION(JM,LM) :: AX,VWT
      REAL*8, DIMENSION(JM) :: SCALEJ,SCALJR
      REAL*8, DIMENSION(LM) :: SCALEL
      REAL*8, DIMENSION(1) :: PL
      CHARACTER*4 DASH,WORD(4)
      DATA DASH/'----'/,WORD/'SUM','MEAN',' ','.1*'/
      INTEGER :: IWORD,J,JHEMI,K,L
      REAL*8 :: FGLOB,GSUM,SUMFAC,GVWT
      CHARACTER XLB*16,CLAT*16,CPRES*16,CBLANK*16,TITLEO*80,TPOW*8
      DATA CLAT/'LATITUDE'/,CPRES/'PRESSURE (MB)'/,CBLANK/' '/
      END SUBROUTINE JLVMAP
      SUBROUTINE DIAGIL
      END SUBROUTINE DIAGIL
      SUBROUTINE ILMAP (sname,lname,unit,PL,AX,SCALEL,LMAX,JWT
     *     ,ISHIFT)
      USE MODEL_COM, only : im,jm,lm,dsig,jdate,jdate0,amon,amon0,jyear
     *     ,jyear0,sige,xlabel
      IMPLICIT NONE
      CHARACTER XLB*80,CWORD*8
      character(len=20), intent(in) :: sname,unit
      character(len=80), intent(in) :: lname
      CHARACTER*64 :: TITLE
      CHARACTER*4, PARAMETER :: DASH='----'
      CHARACTER*4, DIMENSION(2), PARAMETER :: WORD=(/'SUM ','MEAN'/)
      CHARACTER*16, PARAMETER :: CBLANK=' ', CLAT='LONGITUDE',
     *     CPRES='PRESSURE (MB)'
      REAL*8, DIMENSION(LM), INTENT(IN) :: PL,SCALEL
      REAL*8, DIMENSION(IM,LM), INTENT(IN) :: AX
      REAL*8 :: XIL(IM,LM),ZONAL(LM) ! used for post-proc
      REAL*8 :: FGLOB,FLON,GSUM,SDSIG
      REAL*8, DIMENSION(IM) :: ASUM
      INTEGER, DIMENSION(IM) :: MLON
      INTEGER, INTENT(IN) :: JWT,ISHIFT
      INTEGER :: I,L,LMAX
      END SUBROUTINE ILMAP
      SUBROUTINE DIAG7P
      END SUBROUTINE DIAG7P
      SUBROUTINE MEM (SERIES,ITM,MMAX,NUAMAX,NUBMAX,POWER,FPE,VAR,PNU)
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
      END SUBROUTINE MEM
      SUBROUTINE IJ_TITLEX
      END SUBROUTINE IJ_TITLEX
      subroutine IJ_MAPk (k,smap,smapj,gm,igrid,jgrid,irange,
     &     name,lname,units)
      USE MODEL_COM, only :
     &     im,jm,fim,jeq,byim,DTsrc,ptop,
     &     IDACC,
     &     JHOUR,JHOUR0,JDATE,JDATE0,AMON,AMON0,JYEAR,JYEAR0,
     &     NDAY,Itime,Itime0,XLABEL,LRUNID
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,JM) :: anum,adenom,smap
      REAL*8, DIMENSION(JM) :: smapj
      integer, intent(in) :: k
      integer i,j,l,k1,k2,iwt,igrid,jgrid,irange,n1,n2
      character(len=30) name,units
      character(len=80) lname
      real*8 :: gm,nh,sh, off, byiacc, scalek
      integer isumz,isumg
      end subroutine ij_mapk
      subroutine ij_avg (anum,aden,wtij,jgrid,isumz,isumg,
     *                   smap,smapj,gm,nh,sh)
      USE MODEL_COM, only :  im,jm,fim,jeq
      IMPLICIT NONE
      real*8, dimension(im,jm) :: anum,aden,wtij,smap
      real*8, dimension(jm) :: smapj
      real*8, dimension(2) :: znumh,zdenh
      real*8  gm,nh,sh, sumj,wt
      real*8, dimension(jm) :: sumjA, wtA
      integer k,i,j,jgrid,isumz,isumg,jhemi
      end subroutine ij_avg
      SUBROUTINE DIAGIJ
      END SUBROUTINE DIAGIJ
      subroutine maptxt (smap,smapj,gm,irange,title,line,kcol,nlines)
      use model_com, only : im,jm
      IMPLICIT NONE
      real*8, dimension(im,jm) :: smap
      real*8, dimension(jm) :: smapj
      real*8  gm,val,zmax,zmin
      integer irange,kcol,k,k1,ifm,nz1,nz2,n1,n2,ibar,i,j
     *  ,nlines
      character(len=133), dimension(53) :: line
      character(len=40) :: title
      end subroutine maptxt
      SUBROUTINE IJMAP (title,smap,smapj,jgrid)
      USE MODEL_COM, only :
     &     im,jm,NDAY,JHOUR,JHOUR0,JDATE,JDATE0,AMON,AMON0,
     &     JYEAR,JYEAR0,Itime,Itime0,XLABEL,lrunid
      IMPLICIT NONE
      CHARACTER*48 TITLE
      CHARACTER(LEN=3), DIMENSION(IM) :: LINE
      CHARACTER(LEN=9) :: AVG
      REAL*8, DIMENSION(IM,JM) :: SMAP
      REAL*8, DIMENSION(JM) :: SMAPJ
      REAL*8 :: DAYS
      INTEGER :: I,J,jgrid
      END SUBROUTINE IJMAP
      subroutine set_ijout (nmaplets,nmaps,Iord,Qk,iu_Iij)
      IMPLICIT NONE
      character*80 line
      logical Qk(1),Qktmp(1)
      INTEGER Iord(1),nmaplets,nmaps,iu_Iij,k,
     *   n,kmap(3)
      end subroutine set_ijout
      SUBROUTINE DIAGCP
      END SUBROUTINE DIAGCP
      SUBROUTINE DIAG5P
      END SUBROUTINE DIAG5P
      SUBROUTINE DIAGDD
      END SUBROUTINE DIAGDD
      SUBROUTINE DIAGDH
      END SUBROUTINE DIAGDH
      SUBROUTINE DIAG4
      END SUBROUTINE DIAG4
      subroutine KEYVSUMS (QUANT,GSUM,HSUM,ASUM,SUMFAC)
      USE MODEL_COM, only : jm
      implicit none
      CHARACTER(LEN=30) :: QUANT
      REAL*8, DIMENSION(JM) :: ASUM
      REAL*8, DIMENSION(2) :: HSUM
      REAL*8 :: GSUM,SUMFAC
      end subroutine keyvsums
      subroutine keynrl(quant,l,flat)
      USE MODEL_COM, only : jm
      implicit none
      integer :: l
      REAL*8, DIMENSION(JM) :: FLAT
      CHARACTER(LEN=30) :: QUANT
      end subroutine keynrl
      SUBROUTINE IJK_TITLES
      END SUBROUTINE IJK_TITLES
      SUBROUTINE IJKMAP (iu_Iij)
      USE MODEL_COM, only :
     &     im,jm,lm,
     &     PTOP,SIG,PSFMPT,XLABEL,LRUNID,idacc
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
      end subroutine diag_isccp
      subroutine diag_msu
      end subroutine diag_msu
      SUBROUTINE VNTRP1 (KM,P,AIN,  LMA,PE,AOUT)
      implicit none
      integer, intent(in) :: km,lma
      REAL*8, intent(in)  :: P(0:KM),AIN(0:KM),    PE(0:LMA)
      REAL*8, intent(out) :: AOUT(LMA)
      integer k,k1,l
      real*8 pdn,adn,pup,aup,psum,asum
      END subroutine vntrp1
      SUBROUTINE DIAG_GATHER
      END SUBROUTINE DIAG_GATHER
      SUBROUTINE DIAG_SCATTER
      END SUBROUTINE DIAG_SCATTER
      END MODULE DIAG_SERIAL
