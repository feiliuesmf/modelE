      MODULE ODIAG
!@sum  ODIAG ocean diagnostic arrays (incl. dynamic sea ice)
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE OCEAN, only : im,jm,lmo
      USE STRAITS, only : nmst
      USE DAGCOM, only : npts  ! needed for conservation diags
      IMPLICIT NONE
      SAVE
      INTEGER, PARAMETER :: KOIJ=11,KOIJL=22,KOL=6,KOLNST=8,
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
     *     ,IJ_BOSOL,IJ_USTAR,IJ_MUSI,IJ_MVSI
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
C****
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

      write(MODULE_HEADER(lhead+1:80),'(a13,i2,a13,i2,a1,  i2,a5,i2,
     *  a1,i2,a8,i4,a)') 'R8 Oij(im,jm,',koij,'),Oijl(im,jm,',lmo,',',
     *  koijl,'),Ol(',lmo,   ',',kol,'),OLNST(',LMO*NMST*KOLNST,'),it'
   
      SELECT CASE (IACTION)
      CASE (IOWRITE,IOWRITE_MON)  ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,OIJ,OIJL,OL,OLNST,it
      CASE (IOWRITE_SINGLE)    ! output to acc file
        MODULE_HEADER(LHEAD+1:LHEAD+2) = 'R4'
        WRITE (kunit,err=10) MODULE_HEADER,SNGL(OIJ),SNGL(OIJL),SNGL(OL)
     *     ,SNGL(OLNST),it
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
            PRINT*,"Discrepancy in module version",HEADER
     *           ,MODULE_HEADER
            GO TO 10
          END IF
        CASE (ioread,irerun)    ! restarts
          READ (kunit,err=10) HEADER,OIJ,OIJL,OL,OLNST,it
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version",HEADER
     *           ,MODULE_HEADER
            GO TO 10
          END IF
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
!@var M index denoting from where DIAGCO is called
      INTEGER, INTENT(IN) :: M
C****
C**** THE PARAMETER M INDICATES WHEN DIAGCO IS BEING CALLED
C**** M=1  INITIALIZE CURRENT QUANTITY
C****   2  AFTER DYNAMICS
C****   3  AFTER CONDENSATION
C****   4  AFTER RADIATION
C****   5  AFTER PRECIPITATION
C****   6  AFTER LAND SURFACE (INCL. RIVER RUNOFF)
C****   7  AFTER FULL SURFACE INTERACTION
C****   8  AFTER FILTER
C****   9  AFTER STRATOSPHERIC DRAG
C****  10  AFTER OCEAN DYNAMICS
C****  11  AFTER OCEAN SUB-GRIDSCALE PHYS
C****  12  AFTER DAILY
C****
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
      USE ODIAG
      IMPLICIT NONE
      LOGICAL :: QCON(NPTS), T = .TRUE. , F = .FALSE.

C**** Set names for OIJ diagnostics
      IJ_USI=1  ; IJ_VSI=2  ; IJ_DMUI=3  ; IJ_DMVI=4
      IJ_PICE=5 ; IJ_MUSI=6 ; IJ_MVSI=7
      IJ_HBL=8   ; IJ_BO=9  ; IJ_BOSOL=10; IJ_USTAR=11
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

C**** Set up oceanic component conservation diagnostics
C**** Oceanic mass
      QCON=(/ F, F, F, T, F, T, F, T, T, T, T/)
      CALL SET_CON(QCON,"OCN MASS","(10^2 KG/M**2)  ",
     *     "(10^-8 KG/S/M^2)",1d-2,1d8,icon_OMS)
C**** Oceanic angular momentum
      QCON=(/ F, F, F, T, F, T, F, T, T, T, T/)
      CALL SET_CON(QCON,"OCN ANGM","(10^12 J S/M**2)",
     *     "(10^2 W/M^2)    ",1d-12,1d-2,icon_OAM)
C**** Oceanic kinetic energy
      QCON=(/ F, F, F, T, F, T, F, T, T, T, T/)
      CALL SET_CON(QCON,"OCEAN KE","(J/M**2)        ",
     *     "(10^-6 W/M^2)   ",1d0,1d6,icon_OKE)
C**** Oceanic potential enthalpy (heat)
      QCON=(/ F, F, F, T, F, T, F, T, T, T, T/)
      CALL SET_CON(QCON,"OCN HEAT","(10^6 J/M**2)   ",
     *     "(10^-2 KG/S/M^2)",1d-6,1d2,icon_OCE)
C**** Oceanic salt mass
      QCON=(/ F, F, F, T, F, T, F, T, T, T, T/)
      CALL SET_CON(QCON,"OCN SALT","(10 KG/M**2)    ",
     *     "(10^-8 KG/S/M^2)",1d-1,1d8,icon_OSL)
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
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: isum

      OIJ=0 ; OIJL=0 ; OL=0 ; OLNST=0

      return
      END SUBROUTINE reset_odiag

