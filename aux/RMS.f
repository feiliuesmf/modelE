!@sum RMS calculates RMS differences from observations and stores results
!@auth Gary Russell/Gavin Schmidt
!@ver 1.0
C**** Currently this is only good for 72x46 model, should be more general.
C**** If you want to add a target, increase NDIAG and LMAX, enter
C**** the TITLE for the diagnostic, and set up a releavnt OBS line that
C**** describes the observations.
      USE MODEL_COM, only : im,jm,fearth
      USE GEOM, only : geom_b
      USE FILEMANAGER, only : openunit,closeunit
      IMPLICIT NONE
!@var LMAX number of RMS calculations
!@var NDIAG number of target diagnostics
      INTEGER, PARAMETER :: LMAX=33, NDIAG=25
C****
      CHARACTER OBS(LMAX)*63, RMS(0:LMAX,2)*140, FILEIN*80, TITLE1*80,
     *     TITLE2*80, RUN*20, DATE*8, MONFILE(2)*80
      REAL*4 TOTAL(20,2),DATA1(IM,JM),DATA2(IM,JM)
      CHARACTER*7, DIMENSION(2) :: MONTH =(/'January','   July'/)
      CHARACTER*6 :: BLANK ='      '
      DATA RMS /' ',LMAX*' ', ' ',LMAX*' '/
!@var NREC array of positions of targets in diagnostic file      
      INTEGER, DIMENSION(NDIAG) :: NREC
!@var TITLE array of target diagnostics to be found in output
      CHARACTER*45, DIMENSION(NDIAG) :: TITLE = (/
     *     'SNOW COVERAGE (%)                        ',
     *     'PRECIPITATION (mm/day)                   ',
     *     'TOTAL CLOUD COVER (%)                    ',
     *     'NET SOLAR RADIATION, SURF (W/m^2)        ',
     *     'INCIDENT SOLAR RADIATION, SURF (W/m^2)   ',
     *     'SURFACE AIR TEMPERATURE (C)              ',
     *     'U COMPONENT OF SURFACE AIR WIND (m/s)    ',
     *     'V COMPONENT OF SURFACE AIR WIND (m/s)    ',
     *     'LOW LEVEL CLOUDINESS (ISCCP) (%)         ',
     *     'MIDDLE LEVEL CLOUDINESS (ISCCP) (%)      ',
     *     'HIGH LEVEL CLOUDINESS (ISCCP) (%)        ',
     *     'SURFACE AIR SPECIFIC HUMIDITY (10^-4 g/g)',
     *     'SEA LEVEL PRESSURE (mb-1000)             ',
     *     'CLOUD TOP PRESSURE (ISCCP) (mb)          ',
     *     'PLANETARY ALBEDO (%)                     ',
     *     'GROUND ALBEDO IN VISUAL RANGE (%)        ',
     *     'DIURNAL SURF AIR TEMP (C)                ',
     *     'NET SOLAR RADIATION, TOA (W/m^2)         ',
     *     'NET THERMAL RADIATION, TOA (W/m^2)       ',
     *     'TEMPERATURE AT 850mb (C)                 ',
     *     'TEMPERATURE AT 500mb (C)                 ',
     *     'TEMPERATURE AT 300mb (C)                 ',
     *     'SPECIFIC HUMIDITY AT 850mb (g/kg)        ',
     *     'SPECIFIC HUMIDITY AT 500mb (g/kg)        ',
     *     '500 mb HEIGHT (m-5600)                   '/)
C****
C****  Diag  Obsr Domain  N Obseration file     Mon B Wt Fac   Off
C****  ----  ---- ------  - ---------------     --- - -- ----  ---
      DATA OBS /         
     1'PREC  Lega NHGrnd  2 PREC72X46.LEGATES   1 7   99  1    0    ',
     2'PREC  Shea NHGrnd  2 PREC72X46.SHEA      1 7   99  1    0    ',
     3'PREC  GPCP Glob60  2 PREC72X46.GPCP      1 7   99  1    0    ',
     4'SNOWC Robn NHGrnd  1 SNOWC72X46.ROBINSON 1 7 Y 10  1    0    ',
C****		                                            
cNU  x'TCT   ISCC Glob60 xx TCT72X46.ISCCP      1 7    5  1    0    ',
     5'T300  NCEP Global 22 T300.4X5.NCEP       1 7   10  1   273.16',
     6'T500  NCEP NHGrnd 21 T50072X46.NCEP      1 7   15  1   273.16',
     7'T850  NCEP NHGrnd 20 T85072X46.NCEP      1 7   20  1   273.16',
     8'TAS   Lega NHGrnd  6 TSA72X46.LEGATES    1 7   30  1    0    ',
     9'TAS   Shea NHGrnd  6 TSA72X46.SHEA       1 7   30  1    0    ',
     O'TAS   May  NHGrnd  6 TSA72X46.MAY        1 7   35  1    0    ',
     1'TAS   Oort NHGrnd  6 OORT4X5.TSURF       1 7   15  1    0    ',
     2'TASDV May  NHGrnd 17 TSADV72X46.MAY      1 7   25  1    0    ',
cNU  x'TASSD May  NHGrnd xx TSASD72X46.MAY      1 7 Y 10  1    0    ',
C****		                                            
     4'RTEP  ERBE Global 19 OTRP72X46.ERBE      1 7   10 -1.   0    ',
     5'RTEP  ISCC Global 19 OTRP72X46.ISCCP     1 3    5 -1.   0    ',
     6'RSAP  ERBE Global 18 ASRP72X46.ERBE      1 7   10  1    0    ',
     7'RSAP  ISCC Global 18 ASRP72X46.ISCCP     1 3    5  1    0    ',
     8'RSIS  Bish Global  5 BISH4X5             1 7   12  1    0    ',
     9'RSIS  ISCC Glob60  5 ISRS72X46.ISCCP     1 3    5  1    0    ',
     O'ALBP  ERBE Global 15 ALBP72X46.ERBE      1 7   20  1    0    ',
     1'ALBP  ISCC Global 15 ALBP72X46.ISCCP     1 3   10  1    0    ',
     2'ALBS  ISCC Glob60 16 ALBS72X46.ISCCP     1 3 Y 10  1    0    ',
C****		                                            
     3'CLDL  ISCC Glob60  9 CLDL72X46.ISCCP     1 7   10  1    0    ',
     4'CLDH  ISCC Glob60 11 CLDH72X46.ISCCP     1 7   10  1    0    ',
     5'CLDT  ISCC Glob60  3 CLDT72X46.ISCCP     1 7 Y 15  1    0    ',
C****		                                            
     6'PCT   ISCC Glob60 14 PCT72X46.ISCCP      1 7    1  1    0    ',
     7'PSL   Shea Global 13 PSL72X46.SHEA       1 7   15  1    0    ',
     8'PSL   Oort Global 13 OORT4X5.SLP         1 7   10  1    0    ',
     9'Z500  NCEP NHGrnd 25 Z50072X46.NCEP      1 7    1  1   5600. ',
     O'Z500  Oort NHGrnd 25 OORT4X5.H500        1 7 Y .5  1    0    ',
C****		                                            
     1'Q500  NCEP NHGrnd 24 Q50072X46.NCEP      1 7   15  1    0    ',
     2'Q850  NCEP NHGrnd 23 Q85072X46.NCEP      1 7   20  1    0    ',
     3'QS    Oort NHGrnd 12 OORT4X5.QSURF       1 7 Y 10 .1    0    ',
C****		                                            
     4'UVS   Oort NHGrnd  7 OORT4X5.UVSURF      1 7    5  1    0    '/
C****
      INTEGER IARGC,NARGS,M,N,L,NRUN,IREC,K,KM,KMIN,KMAX,NMOD,iu_TOPO
      REAL*4 WRMS,RMSDIF,RRMS,XRMS
C****
      NARGS = IARGC()
      IF(NARGS.lt.3)  GO TO 800
      CALL GEOM_B
C****
C**** Read in FEARTH
C****
      call openunit("TOPO",iu_TOPO,.true.,.true.)
      CALL READT (iu_TOPO,0,FEARTH,IM*JM,FEARTH,3) ! Earth frac. (no LI)
      call closeunit(iu_TOPO)
C****
C**** Use first file to find target diagnostics
C****
      CALL GETARG (1,RUN)
      NRUN = LEN_TRIM(RUN)
      CALL GETARG (2,MONFILE(1))  ! January file
      CALL GETARG (3,MONFILE(2))  ! July file

      OPEN(3,FILE=MONFILE(1),STATUS="OLD",FORM="UNFORMATTED")

      DO N=1,NDIAG
        REWIND (3)
        IREC=1
 10     READ(3,END=20) TITLE1
        IF (INDEX(TITLE1,TITLE(N)).gt.0) THEN
          NREC(N)=IREC
          CYCLE
        ELSE
          IREC=IREC+1
          GOTO 10
        END IF
 20     WRITE (0,*) ' Target diagnostic not found: ',N,TITLE(N)
        WRITE (0,*) ' RMS computation will not be complete.'
        NREC(N)=-1
      END DO

      CLOSE (3)
C****
C**** Calculate new RMS statistics and fill into RMS array
C****

      DO M=1,2
        KMAX = LEN_TRIM(RMS(0,M)) + 11
        KMIN = KMAX-10
        RMS(0,M)(KMIN:KMAX) = RUN(1:NRUN)

        OPEN(3,FILE=MONFILE(M),STATUS="OLD",FORM="UNFORMATTED")
          
C**** loop over calculations
        DO L=1,LMAX
C**** find appropriate target
          READ (OBS(L)(19:21),*) NMOD    ! target index
          REWIND(3)
          IF (NREC(NMOD).gt.-1) THEN
            DO N=1,NREC(NMOD)-1
              READ(3) 
            END DO
            IF(OBS(L)(1:3).eq.'UVS')  THEN
              READ(3) TITLE1,DATA1
              READ(3) TITLE2,DATA2
            ELSE
              READ(3) TITLE1,DATA1
            END IF
            
            RRMS = RMSDIF (M,OBS(L),DATA1,DATA2)
            
            WRITE (RMS(L,M)(KMIN:KMAX),921) RRMS
          ELSE
            WRITE (RMS(L,M)(KMIN:KMAX),'(a)') "*****" 
          END IF
        END DO
        CLOSE (3)
      END DO
C****
C**** Calculate weighted RMS (printed as final line)
C****
      DO M=1,2
        KM = LEN_TRIM(RMS(0,M))/7
        DO K=1,KM
          TOTAL(K,M) = 0.
          DO L=1,LMAX
            if (RMS(L,M)(K*7-6:K*7-6).ne."I" .and. RMS(L,M)(K*7-6:K*7-6)
     *           .ne."*") READ (RMS(L,M)(K*7-6:K*7),921) XRMS 
            READ (OBS(L)(47:49),*) WRMS
            TOTAL(K,M) = TOTAL(K,M) + XRMS*WRMS
          END DO
        END DO
      END DO
C****
C**** Write RMS file to database
C****
      OPEN (3,FILE='/u/cmrun/modelE/RMS/'//trim(RUN),STATUS="UNKNOWN"
     *     ,FORM="FORMATTED") 
      DO M=1,2
        KMAX = LEN_TRIM(RMS(0,M))
        KMIN = KMAX-6
        CALL DATE_AND_TIME (DATE)
        WRITE (3,940) 'RMS    RMS differences for various runs    ' //
     *       DATE(1:4) // '/' // DATE(5:6) // '/' // DATE(7:8)
        WRITE (3,940)
        WRITE (3,940) MONTH(M)
        WRITE (3,940)
        WRITE (3,901) 'Diag  Obsr Domain',RMS(0,M)(1:KMAX)
        WRITE (3,941) ('   ----',K=1,KMAX/7)
        DO L=1,LMAX
          WRITE (3,901) OBS(L)(1:18),RMS(L,M)(1:KMAX)
          IF(OBS(L)(46:46).eq.'Y')  WRITE (3,901)
        END DO
        WRITE (3,942) ( '  -----',K=1,KMAX/7)
        WRITE (3,943) (NINT(TOTAL(K,M)),K=1,KMAX/7)
        IF(M.eq.1)  WRITE (3,944) CHAR(12) !  write form feed
      END DO
      CLOSE (3)
      STOP
C****
  800 WRITE (0,*)
     *' Usage: RMS run jan_diag_file jul_diag_file           2001/10/03'
      WRITE (0,*)
     *' Calculate RMS statistics for JAN and JUL diagnostics from RUN'
      WRITE (0,*)
     *' The DIAG_FILES should be the lat-lon diagnostics for the '
      WRITE (0,*)
     *' appropriate period, for January and July'
      WRITE (0,*)
     *' Example:  RMS E001 JAN1950-1955.ijE001 JUL1950-1955.ijE001'
      GO TO 999
  801 WRITE (0,*) ' Error opening datafile: ',FILEIN(1:40)
      STOP 801
C****
  901 FORMAT (A17,2X,A)
  921 FORMAT (F7.2)
  940 FORMAT (A)
  941 FORMAT ('----  ---- ------',1X,20A7)
  942 FORMAT (19X,20A7)
  943 FORMAT (19X,20I7)
  944 FORMAT (A1,$)
  999 END

      FUNCTION RMSDIF (M,OBS,UM,VM)
!@sum RMSDIF calculates the RMS difference between a model data record
!@+ and an observation data record
!@auth Gary Russell/Gavin Schmidt
C**** Input: M = 1 for January, 2 for July
C****      UM,VM data arrays
      USE CONSTANT, only : SKIP=>undef
      USE MODEL_COM, only : im,jm,fearth
      USE GEOM,only : dxyp
      IMPLICIT NONE
      REAL*4, PARAMETER :: SKIPOBS=-999999.
      CHARACTER OBS*64, FILEIN*80, TITLE*80
      REAL*4 UM(IM,JM),VM(IM,JM), UO(IM,JM),VO(IM,JM),FAC,OFFSET,W,Q
     *     ,RMSDIF,DLAT
      INTEGER N,NOBS,I,J,J60,M
C**** factors and offsets to transform model output to commensurate units	
      READ(OBS(56:61),*) OFFSET
      READ(OBS(51:54),*) FAC
C****
C**** Open model and observation data files
C****
      FILEIN = '/raid4/OBS/' // OBS(22:41)
      OPEN (2,FILE=FILEIN,FORM='UNFORMATTED',STATUS='OLD',ERR=801)
C****
C**** Skip through unused records
C****
      READ (OBS(40+2*M:40+2*M),*) NOBS
      DO N=1,NOBS-1
        READ (2)
      END DO
C****
C**** Determine kind of RMS
C****
      IF(OBS( 1: 3).eq.'UVS'   )  GO TO 400
      IF(OBS(12:17).eq.'NHGrnd')  GO TO 300
      IF(OBS(12:17).eq.'Glob60')  GO TO 200
C     IF(OBS(12:17).eq.'Global')  GO TO 100
C****
C**** Calculate global RMS
C****
c  100 READ (1,ERR=810) TITLE,UM
 100  READ (2,ERR=810) TITLE,UO
      W = 0.
      Q = 0.
      DO J=1,JM
      DO I=1,IM
        IF(UM(I,J).le.SKIP .or. UO(I,J).eq.SKIPOBS) CYCLE
        W = W + DXYP(J)
        Q = Q + DXYP(J)*(UM(I,J)*FAC+OFFSET-UO(I,J))**2
      END DO
      END DO
      GO TO 500
C****
C**** Calculate global RMS from -60 to +60
C****
c  200 READ (1,ERR=810) TITLE,UM
 200  READ (2,ERR=810) TITLE,UO
      DLAT = .5*NINT(360./JM)
      J60  = 1 + NINT(JM/2 - 60./DLAT)
      W = 0.
      Q = 0.
      DO J=J60,1+JM-J60
      DO I=1,IM
        IF(UM(I,J).le.SKIP .or. UO(I,J).eq.SKIPOBS) CYCLE
        W = W + DXYP(J)
        Q = Q + DXYP(J)*(UM(I,J)*FAC+OFFSET-UO(I,J))**2
      END DO
      END DO
      GO TO 500
C****
C**** Calculate RMS over northern hemisphere ground
C****
c  300 READ (1,ERR=810) TITLE,UM
 300  READ (2,ERR=810) TITLE,UO
      W = 0.
      Q = 0.
      DO J=JM/2+1,JM
      DO I=1,IM
        IF(UM(I,J).le.SKIP .or. UO(I,J).eq.SKIPOBS) CYCLE
        W = W + FEARTH(I,J)*DXYP(J)
        Q = Q + FEARTH(I,J)*DXYP(J)*(UM(I,J)*FAC+OFFSET-UO(I,J))**2
      END DO
      END DO
      GO TO 500
C****
C**** Calculate vector RMS over northern hemisphere ground
C****
c  400 READ (1,ERR=810) TITLE,UM
c      READ (1,ERR=810) TITLE,VM
 400  READ (2,ERR=810) TITLE,UO,VO
      W = 0.
      Q = 0.
      DO J=JM/2+1,JM
      DO I=1,IM
        IF(UM(I,J).le.SKIP .or. UO(I,J).eq.SKIPOBS .or.
     *       VM(I,J).le.SKIP .or. VO(I,J).eq.SKIPOBS) CYCLE
        W = W + FEARTH(I,J)*DXYP(J)
        Q = Q + FEARTH(I,J)*DXYP(J)*((UM(I,J)*FAC+OFFSET-UO(I,J))**2 +
     +       (VM(I,J)*FAC+OFFSET-VO(I,J))**2)
      END DO
      END DO
C****
C**** Final calculation of RMS and close datafiles
C****
  500 RMSDIF = SQRT(Q/(W+1.e-20))
      CLOSE (2)
      RETURN
C****
  801 WRITE (0,*) ' Error opening datafile: ',FILEIN(1:40)
      STOP 801
  810 WRITE (0,*) ' Error reading datafile: ',OBS
      STOP 810
      END




