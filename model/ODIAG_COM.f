#include "rundeck_opts.h"

      MODULE ODIAG
!@sum  ODIAG ocean diagnostic arrays (incl. dynamic sea ice)
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
#ifdef TRACERS_OCEAN
      USE TRACER_COM, only : ntm
#endif
      USE OCEAN, only : im,jm,lmo
      USE STRAITS, only : nmst
      USE DAGCOM, only : npts  ! needed for conservation diags
      IMPLICIT NONE
      SAVE
      INTEGER, PARAMETER :: KOIJ=17,KOIJL=22,KOL=6,KOLNST=8,
     *     KACCO=IM*JM*KOIJ + IM*JM*LMO*KOIJL + LMO*KOL + LMO*NMST
     *     *KOLNST
!@var OIJ   lat-lon ocean diagnostics (on ocean grid)
!@var OIJL  3-dimensional ocean diagnostics
!@var OL    vertical ocean diagnostics
!@var OLNST strait diagnostics
      REAL*8, DIMENSION(IM,JM,KOIJ)  :: OIJ
      REAL*8, DIMENSION(IM,JM,LMO,KOIJL) :: OIJL
      REAL*8, DIMENSION(LMO,KOL)   :: OL
      REAL*8, DIMENSION(LMO,NMST,KOLNST):: OLNST
!@var IJ_xxx Names for OIJ diagnostics
      INTEGER IJ_USI,IJ_VSI,IJ_DMUI,IJ_DMVI,IJ_PICE,IJ_HBL,IJ_BO
     *     ,IJ_BOSOL,IJ_USTAR,IJ_MUSI,IJ_MVSI,IJ_HUSI,IJ_HVSI,IJ_SUSI
     *     ,IJ_SVSI,IJ_RSI,IJ_SSH
!@var lname_oij Long names for OIJ diagnostics
      CHARACTER*50, DIMENSION(KOIJ) :: LNAME_OIJ
!@var sname_oij Short names for OIJ diagnostics
      CHARACTER*30, DIMENSION(KOIJ) :: SNAME_OIJ
!@var units_oij Units for OIJ diagnostics
      CHARACTER*50, DIMENSION(KOIJ) :: UNITS_OIJ
!@var ia_oij IDACC numbers for OIJ diagnostics
      INTEGER, DIMENSION(KOIJ) :: IA_OIJ
!@var scale_oij scales for OIJ diagnostics
      REAL*8, DIMENSION(KOIJ) :: SCALE_OIJ
!@var ijgrid_oij Grid descriptor for OIJ diagnostics
       INTEGER, DIMENSION(KOIJ) :: IJGRID_OIJ

!@var IJL_xxx Names for OIJL diagnostics
      INTEGER IJL_MO,IJL_G0M,IJL_S0M,IJL_GFLX,IJL_SFLX,IJL_MFU,IJL_MFV
     *     ,IJL_MFW,IJL_GGMFL,IJL_SGMFL,IJL_KVM,IJL_KVG,IJL_WGFL
     *     ,IJL_WSFL
!@var LN_xxx Names for OLNST diagnostics
      INTEGER LN_KVM,LN_KVG,LN_WGFL,LN_WSFL,LN_MFLX,LN_GFLX,LN_SFLX
     *     ,LN_ICFL
!@var L_xxx Names for OL diagnostics
      INTEGER L_RHO,L_TEMP,L_SALT
!@var icon_xx indexes for conservation quantities
      INTEGER icon_OCE,icon_OKE,icon_OAM,icon_OMS,icon_OSL
!@var kbasin integer index of which basin a particular ocean point is in
      INTEGER, DIMENSION(IM,JM) :: KBASIN
!@var XLB label for diagnostic titles
      CHARACTER XLB*30
!@var FLAT latitude values on primary and secondary ocean grids
      REAL*8, DIMENSION(JM,2) :: FLAT
!@var FLON longitude values on primary and secondary ocean grids
      REAL*8, DIMENSION(IM,2) :: FLON
!@var FLEV depth values on primary and secondary ocean grids
      REAL*8, DIMENSION(0:LMO,2) :: FLEV
!@var iu_otj unit number for ascii output of ocean transports
      INTEGER iu_otj
!@var BASIN names of ocean basins for diag output
      CHARACTER*16, DIMENSION(4) :: BASIN=
     *     (/"Atlantic","Pacific ","Indian  ","Global  "/)
C****

#ifdef TRACERS_OCEAN
!@var KTOIJL number of 3-dimensional ocean tracer diagnostics
      INTEGER, PARAMETER :: KTOIJL=10
!@var TOIJL  3-dimensional ocean tracer diagnostics
      REAL*8, DIMENSION(IM,JM,LMO,KTOIJL,NTM)  :: TOIJL
!@var toijl_xxx indices for TOIJL diags
      INTEGER, PARAMETER :: toijl_conc=1,toijl_tflx=2,toijl_gmfl=6
!@var KTOIJ number of lat/lon ocean tracer diagnostics
      INTEGER, PARAMETER :: KTOIJ=2
!@var TOIJ  lat/lon ocean tracer diagnostics
      REAL*8, DIMENSION(IM,JM,KTOIJ,NTM)  :: TOIJ
!@var toij_xxx indices for TOIJL diags
      INTEGER :: toij_tusi,toij_tvsi
!@var lname_toij Long names for TOIJ diagnostics
      CHARACTER*50, DIMENSION(KTOIJ) :: LNAME_TOIJ
!@var sname_toij Short names for TOIJ diagnostics
      CHARACTER*30, DIMENSION(KTOIJ) :: SNAME_TOIJ
!@var units_toij Units for TOIJ diagnostics
      CHARACTER*50, DIMENSION(KTOIJ) :: UNITS_TOIJ
!@var ia_toij IDACC numbers for TOIJ diagnostics
      INTEGER, DIMENSION(KTOIJ) :: IA_TOIJ
!@var scale_toij scales for TOIJ diagnostics
      REAL*8, DIMENSION(KTOIJ) :: SCALE_TOIJ
!@var ijgrid_toij Grid descriptor for TOIJ diagnostics
       INTEGER, DIMENSION(KTOIJ) :: IJGRID_TOIJ
!@var TLNST strait diagnostics
      REAL*8, DIMENSION(LMO,NMST,KOLNST,NTM):: TLNST
#endif
      END MODULE ODIAG

      SUBROUTINE io_ocdiag(kunit,it,iaction,ioerr)
!@sum  io_ocdiag reads and writes ocean diagnostic arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,iowrite_mon,iowrite_single
     *     ,irsfic,irerun,ioread_single,lhead
      USE ODIAG
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "OCDIAG01"
!@var it input/ouput value of hour
      INTEGER, INTENT(INOUT) :: it
!@var OIJ4,OIJL4,OLNST4,OL4 dummy arrays for reading diag. files
      REAL*4, DIMENSION(IM,JM,KOIJ)  :: OIJ4
      REAL*4, DIMENSION(IM,JM,LMO,KOIJL) :: OIJL4
      REAL*4, DIMENSION(LMO,KOL)   :: OL4
      REAL*4, DIMENSION(LMO,NMST,KOLNST):: OLNST4
#ifdef TRACERS_OCEAN
      REAL*4, DIMENSION(IM,JM,LMO,KTOIJL,NTM)  :: TOIJL4
      REAL*4, DIMENSION(LMO,NMST,KOLNST,NTM):: TLNST4
      REAL*4, DIMENSION(IM,JM,KTOIJ,NTM)  :: TOIJ4
!@var TR_HEADER Character string label for individual tracer records
      CHARACTER*80 :: TR_HEADER, TR_MODULE_HEADER = "TROCDIAG01"

      write(TR_MODULE_HEADER(lhead+1:80)
     *     ,'(a19,i2,a1,i2,a14,i2,a1,i2,a9,i4,a4)')
     *     'R8 Toijl(im,jm,lmo,',ktoijl,',',ntm,'), Toij(im,jm,',
     *     ktoij,',',ntm,'), Tlnst(',LMO*NMST*KOLNST*NTM,'),it'
#endif
      write(MODULE_HEADER(lhead+1:80),'(a13,i2,a13,i2,a1,  i2,a5,i2,
     *  a1,i2,a8,i4,a)') 'R8 Oij(im,jm,',koij,'),Oijl(im,jm,',lmo,',',
     *  koijl,'),Ol(',lmo,   ',',kol,'),OLNST(',LMO*NMST*KOLNST,'),it'

      SELECT CASE (IACTION)
      CASE (IOWRITE,IOWRITE_MON)  ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,OIJ,OIJL,OL,OLNST,it
#ifdef TRACERS_OCEAN
        WRITE (kunit,err=10) TR_MODULE_HEADER,TOIJL,TOIJ,TLNST,it
#endif
      CASE (IOWRITE_SINGLE)    ! output to acc file
        MODULE_HEADER(LHEAD+1:LHEAD+2) = 'R4'
        WRITE (kunit,err=10) MODULE_HEADER,REAL(OIJ,KIND=4),
     *        REAL(OIJL,KIND=4),REAL(OL,KIND=4),REAL(OLNST,KIND=4),it
#ifdef TRACERS_OCEAN
        TR_MODULE_HEADER(LHEAD+1:LHEAD+2) = 'R4'
        WRITE (kunit,err=10) TR_MODULE_HEADER,REAL(TOIJL,KIND=4),
     *       REAL(TOIJ,KIND=4),REAL(TLNST,KIND=4),it
#endif
      CASE (IOREAD:)            ! input from restart file
        SELECT CASE (IACTION)
        CASE (IRSFIC)           ! initial conditions
        CASE (ioread_single)    ! accumulate diagnostic files
          READ (kunit,err=10) HEADER,OIJ4,OIJL4,OL4,OLNST4,it
C**** accumulate diagnostics
          OIJ=OIJ+OIJ4
          OIJL=OIJL+OIJL4
          OL=OL+OL4
          OLNST=OLNST+OLNST4
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER
     *           ,MODULE_HEADER
            GO TO 10
          END IF
#ifdef TRACERS_OCEAN
          READ (kunit,err=10) TR_HEADER,TOIJL4,TOIJ4,TLNST4,it
C**** accumulate diagnostics
          TOIJL=TOIJL+TOIJL4
          TOIJ =TOIJ +TOIJ4
          TLNST=TLNST+TLNST4
          IF (TR_HEADER(1:LHEAD).NE.TR_MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",TR_HEADER
     *           ,TR_MODULE_HEADER
            GO TO 10
          END IF
#endif
        CASE (ioread,irerun)    ! restarts
          READ (kunit,err=10) HEADER,OIJ,OIJL,OL,OLNST,it
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version",HEADER
     *           ,MODULE_HEADER
            GO TO 10
          END IF
#ifdef TRACERS_OCEAN
          READ (kunit,err=10) TR_HEADER,TOIJL,TLNST,it
          IF (TR_HEADER(1:LHEAD).NE.TR_MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",TR_HEADER
     *           ,TR_MODULE_HEADER
            GO TO 10
          END IF
#endif
        END SELECT
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
C****
      END SUBROUTINE io_ocdiag

      SUBROUTINE DIAGCO (M)
!@sum  DIAGCO Keeps track of the ocean conservation properties
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE ODIAG, only : icon_OCE,icon_OKE,icon_OMS,icon_OSL,icon_OAM
      IMPLICIT NONE
!@var M index denoting from where DIAGCO is called (see DIAGCA)
      INTEGER, INTENT(IN) :: M
      REAL*8, EXTERNAL :: conserv_OCE,conserv_OKE,conserv_OMS
     *     ,conserv_OSL,conserv_OAM

C**** OCEAN MASS
      CALL conserv_DIAG(M,conserv_OMS,icon_OMS)

C**** OCEAN ANGULAR MOMENTUM
      CALL conserv_DIAG(M,conserv_OAM,icon_OAM)

C**** OCEAN KINETIC ENERGY
      CALL conserv_DIAG(M,conserv_OKE,icon_OKE)

C**** OCEAN POTENTIAL ENTHALPY
      CALL conserv_DIAG(M,conserv_OCE,icon_OCE)

C**** OCEAN SALT
      CALL conserv_DIAG(M,conserv_OSL,icon_OSL)

C****
      RETURN
      END SUBROUTINE DIAGCO

      SUBROUTINE init_ODIAG
!@sum  init_ODIAG initialises ocean diagnostics
!@auth Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only : bygrav
      USE MODEL_COM, only : dtsrc
      USE DAGCOM, only : ia_src
      USE ODIAG
      IMPLICIT NONE
      LOGICAL :: QCON(NPTS), T = .TRUE. , F = .FALSE.
      INTEGER k

C**** Set names for OIJL diagnostics
      IJL_MO=1  ; IJL_G0M=2  ; IJL_S0M=3   ; IJL_MFU=4    ; IJL_MFV=5
      IJL_MFW=6 ; IJL_KVM=13; IJL_KVG=14 ; IJL_WGFL=15 ; IJL_WSFL=16
C**** These flux diagnostics need 3 spots each
      IJL_GFLX=7 ; IJL_SFLX=10 ; IJL_GGMFL=17 ; IJL_SGMFL=20
C**** Set names for OLNST diagnostics
      LN_KVM=1  ; LN_KVG=2  ; LN_WGFL=3 ; LN_WSFL=4
      LN_MFLX=5 ; LN_GFLX=6 ; LN_SFLX=7 ; LN_ICFL=8
C**** Set names for OL diagnostics
      L_RHO=1   ; L_TEMP=2  ; L_SALT=3

C**** set properties for OIJ diagnostics
      k=0

      k=k+1
      IJ_USI=k
      lname_oij(k)="Sea ice EW velocity x POICEU"
      sname_oij(k)="oij_usi"
      units_oij(k)="m/s"
      ia_oij(k)=ia_src
      scale_oij(k)=1.
      ijgrid_oij(k)=2

      k=k+1
      IJ_VSI=k
      lname_oij(k)="Sea ice NS velocity x POICEV"
      sname_oij(k)="oij_vsi"
      units_oij(k)="m/s"
      ia_oij(k)=ia_src
      scale_oij(k)=1.
      ijgrid_oij(k)=2

      k=k+1
      IJ_DMUI=k
      lname_oij(k)="Ice-ocean EW stress"
      sname_oij(k)="oij_dmui"
      units_oij(k)="kg/m s^2"
      ia_oij(k)=ia_src
      scale_oij(k)=1./dtsrc
      ijgrid_oij(k)=1

      k=k+1
      IJ_DMVI=k
      lname_oij(k)="Ice-ocean NS stress"
      sname_oij(k)="oij_dmvi"
      units_oij(k)="kg/m s^2"
      ia_oij(k)=ia_src
      scale_oij(k)=1./dtsrc
      ijgrid_oij(k)=1

      k=k+1
      IJ_PICE=k
      lname_oij(k)="Sea ice internal pressure x POICE"
      sname_oij(k)="oij_psi"
      units_oij(k)="10^3 kg/m s^2"
      ia_oij(k)=ia_src
      scale_oij(k)=1d-3
      ijgrid_oij(k)=1

      k=k+1
      IJ_MUSI=k
      lname_oij(k)="Sea ice NS mass flux"
      sname_oij(k)="oij_musi"
      units_oij(k)="kg/s"
      ia_oij(k)=ia_src
      scale_oij(k)=1./dtsrc
      ijgrid_oij(k)=2

      k=k+1
      IJ_MVSI=k
      lname_oij(k)="Sea ice EW mass flux"
      sname_oij(k)="oij_mvsi"
      units_oij(k)="kg/s"
      ia_oij(k)=ia_src
      scale_oij(k)=1./dtsrc
      ijgrid_oij(k)=2

      k=k+1
      IJ_HUSI=k
      lname_oij(k)="Sea ice NS heat flux"
      sname_oij(k)="oij_husi"
      units_oij(k)="W"
      ia_oij(k)=ia_src
      scale_oij(k)=1./dtsrc
      ijgrid_oij(k)=2

      k=k+1
      IJ_HVSI=k
      lname_oij(k)="Sea ice EW heat flux"
      sname_oij(k)="oij_hvsi"
      units_oij(k)="W"
      ia_oij(k)=ia_src
      scale_oij(k)=1./dtsrc
      ijgrid_oij(k)=2

      k=k+1
      IJ_SUSI=k
      lname_oij(k)="Sea ice NS salt flux"
      sname_oij(k)="oij_susi"
      units_oij(k)="kg/s"
      ia_oij(k)=ia_src
      scale_oij(k)=1./dtsrc
      ijgrid_oij(k)=2

      k=k+1
      IJ_SVSI=k
      lname_oij(k)="Sea ice EW salt flux"
      sname_oij(k)="oij_svsi"
      units_oij(k)="kg/s"
      ia_oij(k)=ia_src
      scale_oij(k)=1./dtsrc
      ijgrid_oij(k)=2

      k=k+1
      IJ_HBL=k
      lname_oij(k)="Ocean Boundary layer depth (KPP) x PO4"
      sname_oij(k)="oij_hbl"
      units_oij(k)="m"
      ia_oij(k)=ia_src
      scale_oij(k)=0.25
      ijgrid_oij(k)=1

      k=k+1
      IJ_BO=k
      lname_oij(k)="Surface buoyancy forcing (KPP) x PO4"
      sname_oij(k)="oij_bo"
      units_oij(k)="10^-7 m^2/s^3"
      ia_oij(k)=ia_src
      scale_oij(k)=0.25*1d7
      ijgrid_oij(k)=1

      k=k+1
      IJ_BOSOL=k
      lname_oij(k)="Surface solar buoyancy flux x PO4"
      sname_oij(k)="oij_bosol"
      units_oij(k)="10^-7 m^2/s^3"
      ia_oij(k)=ia_src
      scale_oij(k)=0.25*1d7
      ijgrid_oij(k)=1

      k=k+1
      IJ_USTAR=k
      lname_oij(k)="Surface friction speed x PO4"
      sname_oij(k)="oij_ustar"
      units_oij(k)="m/s"
      ia_oij(k)=ia_src
      scale_oij(k)=0.25
      ijgrid_oij(k)=1

      k=k+1
      IJ_SSH=k
      lname_oij(k)="Ocean surface height"
      sname_oij(k)="oij_ssh"
      units_oij(k)="m"
      ia_oij(k)=ia_src
      scale_oij(k)=bygrav
      ijgrid_oij(k)=1

      k=k+1
      IJ_RSI=k
      lname_oij(k)="Ocean ice fraction"
      sname_oij(k)="oij_rsi"
      units_oij(k)="%"
      ia_oij(k)=ia_src
      scale_oij(k)=100.
      ijgrid_oij(k)=1

      if (k.gt.KOIJ) then
        write(6,*) "Too many OIJ diagnostics: increase KOIJ to at least"
     *       ,k
        stop "OIJ diagnostic error"
      end if

#ifdef TRACERS_OCEAN
C**** simple tracer diags same description for all tracers
C**** set properties for TOIJ diagnostics
      k=0
     
      k=k+1
      TOIJ_TUSI=k
      lname_toij(k)="Sea ice NS tracer flux"
      sname_toij(k)="toij_tusi"
      units_toij(k)="kg/s"
      ia_toij(k)=ia_src
      scale_toij(k)=1./dtsrc
      ijgrid_toij(k)=2

      k=k+1
      TOIJ_TVSI=k
      lname_toij(k)="Sea ice EW tracer flux"
      sname_toij(k)="toij_tvsi"
      units_toij(k)="kg/s"
      ia_toij(k)=ia_src
      scale_toij(k)=1./dtsrc
      ijgrid_toij(k)=2

      if (k.gt.KTOIJ) then
        write(6,*) "Too many TOIJ diags: increase KTOIJ to at least"
     *       ,k
        stop "TOIJ diagnostic error"
      end if
#endif

C**** Set up oceanic component conservation diagnostics
C**** Oceanic mass
      QCON=(/ F, F, F, T, F, T, F, T, T, T, T/)
      CALL SET_CON(QCON,"OCN MASS","(10**2 KG/M^2) ",
     *     "(10**-8 KG/SM^2)",1d-2,1d8,icon_OMS)
C**** Oceanic angular momentum
      QCON=(/ F, F, F, T, F, T, F, T, T, T, T/)
      CALL SET_CON(QCON,"OCN AM  ","(10**12 JS/M^2)",
     *     "(10**2 J/M^2)   ",1d-12,1d-2,icon_OAM)
C**** Oceanic kinetic energy
      QCON=(/ F, F, F, T, F, T, F, T, T, T, T/)
      CALL SET_CON(QCON,"OCEAN KE","(J/M^2)        ",
     *     "(10**-6 W/M^2)  ",1d0,1d6,icon_OKE)
C**** Oceanic potential enthalpy (heat)
      QCON=(/ F, F, F, T, F, T, F, T, T, T, T/)
      CALL SET_CON(QCON,"OCN HEAT","(10**6 J/M^2)  ",
     *     "(10**-2 KG/SM^2)",1d-6,1d2,icon_OCE)
C**** Oceanic salt mass
      QCON=(/ F, F, F, T, F, T, F, T, T, T, T/)
      CALL SET_CON(QCON,"OCN SALT","(10 KG/M^2)    ",
     *     "(10**-9 KG/SM^2)",1d-1,1d9,icon_OSL)
C**** Initialise ocean basins
      CALL OBASIN

C****
      RETURN
      END SUBROUTINE init_ODIAG

      SUBROUTINE reset_odiag(isum)
!@sum  reset_odiag zeros out ocean diagnostics if needed
!@auth G. Schmidt
!@ver  1.0
      USE ODIAG, only : oij,oijl,ol,olnst
#ifdef TRACERS_OCEAN
     *     ,toijl,toij,tlnst
#endif
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: isum  ! needed for plug-play compatibility

      OIJ=0. ; OIJL=0. ; OL=0. ; OLNST=0. 

#ifdef TRACERS_OCEAN
      TOIJL=0. ; TLNST = 0. ; TOIJ=0.
#endif
      return
      END SUBROUTINE reset_odiag

