      MODULE STRAITS
!@sum  STRAITS ocean strait related variables 
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE OCEAN, only : lmo
      USE SEAICE, only : lmi
      IMPLICIT NONE
      SAVE
C**** These values are highly resolution dependent
C****
C****     Strait           From         To       LM    Width
C****     ------           ----         --       --    -----
C****  1  Fury & Hecla   19,42 ES    20,40 WN     2    20000
C****  2  Nares          22,43 EN    24,44 WS     5     5000
C****  3  Gibraltar      35,32 EN    37,33 WS     5    25000
C****  4  English        36,36 EN    37,37 WS     2    35000
C****  5  Kattegat       38,38 EN    40,38 WS     2    60000
C****  6  Bosporous      42,33 EN    43,34 WS     2     6000
C****  7  Red Sea        44,29 ES    45,28 WN     6   250000
C****  8  Bab al Mandab  45,28 ES    46,27 WN     6    25000
C****  9  Hormuz         47,30 ES    49,29 WN     2    50000
C**** 10  Malacca        56,25 EN    58,24 WS     3    50000
C**** 11  Korea          62,32 EN    63,33 WS     4   170000
C**** 12  Soya-kaikyo    64,34 EN    65,35 WS     2    40000
C****
      INTEGER, PARAMETER :: NMST=12  !@param NMST no. of ocean straits
!@var MMST mass of water in strait (kg)
!@var MUST mass flux of water in strait (kg/s)
!@var G0MST,GXMST,GYMST pot. enthalpy of water in strait (+ moments) (J)
!@var S0MST,SXMST,SYMST salinity of water in strait (+ moments) (kg)
      REAL*8, DIMENSION(LMO,NMST) :: MMST,MUST,G0MST,GXMST,GZMST,S0MST
     *     ,SXMST,SZMST 

!@var WIST width of strait (m)
!@var DIST distance along strait (m)
!@var DISTPG distance between centre points of adjoining ocean boxes (m)
      REAL*8, DIMENSION(NMST) :: DIST,DISTPG,
     *     WIST = (/  2d4,   5d3, 2.5d4, 3.5d4,   6d4, 6d3,
     *              2.5d5, 2.5d4,   5d4,   5d4, 1.7d5, 4d4/)

!@var XST,YST 
      REAL*8, DIMENSION(NMST,2) :: XST = RESHAPE( (/
     *     0d0  ,6d-1, 6d-1,  1d0,  1d0, 0d0, 6d-1, 6d-1,  1d0,6d-1,
     *     0d0  ,6d-1,  0d0,-8d-1,-8d-1, 0d0,-6d-1,-8d-1,-8d-1,-1d0,
     *     -7d-1,-1d0,-6d-1,-1d0 /), (/NMST,2/) ),
     *     YST = RESHAPE( (/
     *     -1d0,8d-1, 8d-1,  0d0,  0d0, 1d0,-8d-1,-8d-1, 0d0,-8d-1,
     *      1d0,8d-1,  1d0,-6d-1,-6d-1,-1d0,-8d-1,-6d-1,6d-1,  0d0,
     *     7d-1, 0d0,-8d-1,  0d0 /), (/NMST,2/) )

!@var IST,JST i,j coordinates of ends of straits
      INTEGER, DIMENSION(NMST,2) ::
     *     IST = RESHAPE( (/
     *     19, 22, 35, 36, 38, 42, 44, 45, 47, 56, 62, 64,
     *     20, 24, 37, 37, 40, 43, 45, 46, 49, 58, 63, 65/),
     *     (/NMST,2/) ),
     *     JST = RESHAPE( (/
     *     42, 43, 32, 36, 38, 33, 29, 28, 30, 25, 32, 34,
     *     40, 44, 33, 37, 38, 34, 28, 27, 29, 24, 33, 35/),
     *     (/NMST,2/) )

!@var LMST no. of levels in strait
      INTEGER, DIMENSION(NMST) :: LMST = (/ 
     *     2,  5,  5,  2,  2,  2,  6,  6,  2,  3,  4,  2/)

!@var RSIST Sea ice fraction in strait
!@var RSIXST Center of sea ice in strait (m) 
!@var MSIST Mass of ice within strait (kg/m**2)
!@var HSIST Enthalphy of ice within strait (J/m**2)
      REAL*8, DIMENSION(NMST) :: RSIST,RSIXST
      REAL*8, DIMENSION(2,NMST) :: MSIST
      REAL*8, DIMENSION(LMI,NMST) :: HSIST

      END MODULE STRAITS

      SUBROUTINE io_straits(kunit,iaction,ioerr)
!@sum  io_straits reads and writes ocean straits arrays to file
!@auth Gavin Schmidt
!@ver  1.0
      USE MODEL_COM, only : ioread,iowrite,irsfic,irerun
      USE STRAITS
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*8 :: HEADER, MODULE_HEADER = "OCSTR01"

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER,MUST,G0MST,GXMST,GZMST,S0MST
     *       ,SXMST,SZMST,RSIST,RSIXST,MSIST,HSIST
      CASE (IOREAD:)            ! input from restart file
        SELECT CASE (IACTION)
        CASE (IRSFIC)           ! initial conditions
        CASE (ioread,irerun)    ! restarts
          READ (kunit,err=10) HEADER,MUST,G0MST,GXMST,GZMST,S0MST
     *         ,SXMST,SZMST,RSIST,RSIXST,MSIST,HSIST
          IF (HEADER.NE.MODULE_HEADER) THEN
            PRINT*,"Discrepancy in module version ",HEADER
     *           ,MODULE_HEADER
            GO TO 10
          END IF
        END SELECT
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
C****
      END SUBROUTINE io_straits
