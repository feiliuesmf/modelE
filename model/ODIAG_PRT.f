      SUBROUTINE diag_OCEAN
!@sum  diag_OCEAN prints out diagnostics for ocean
!@auth Gavin Schmidt/Gary Russell
!@ver  1.0 
C**** Note this is an incorporation and modification of the stand alone
C**** ocean diagnostic programs from Gary. All diagnostics are on the
C**** ocean grid.
      USE CONSTANT, only : skip=>undef
      USE MODEL_COM, only : xlabel,jmon0,jyear0,idacc,jdate0,amon0
     *     ,jdate,amon,jyear 
      USE OCEAN, only : im,jm,lmo,focean,dxyp=>dxypo,ndyno,dts,dto
     *     ,imaxj,lmm,ze
      USE ODIAG
      USE STRAITS, only : nmst,wist,dist,lmst,name_st
      IMPLICIT NONE
      REAL*4, DIMENSION(IM,JM) :: Q,SFIJM,SFIJS
      REAL*4, DIMENSION(IM,JM,LMO) :: Q3
      REAL*8, DIMENSION(LMO,NMST) :: AS
      REAL*4, DIMENSION(JM-1,0:LMO,0:4) :: SFM,SFS
      INTEGER I,J,K,L,LMINEF,LMAXEF,KVMF(LMO),LMAXNF,LMINMF,NS,N,KB
     *     ,LMINSF,LMAXSF,KVDC(LMO),KCMF(LMO),LMAXMF,NOL(LMO)
     *     ,LMSTMIN(LMO),SUMORMN(LMO),JEQ,JDLAT,JLAT(JM)
      REAL*4 R4LEV(0:LMO)
      REAL*8 GOS,SOS,FAC,FACST,GSMAX,GSMIN,CKMIN,CKMAX,ACMIN,ACMAX
     *     ,SCALEO(10),TSUM,TEMGS,scalej(jm),scalel(lmo),AJLTMP(JM,LMO)
      CHARACTER NAME(KOLNST)*50,OUTMON,OUTYR,NAMESF(8)*50,TITLE*80
     *     ,NAMEL(20)*80,sname*30,units*50
      LOGICAL QKV
C**** depth/latitude diagnostics by basin
      DATA NAMESF/
     1 'ATLANTIC MASS STREAM FUNCTION',
     2 'PACIFIC MASS STREAM FUNCTION',
     3 'INDIAN MASS STREAM FUNCTION',
     4 'GLOBAL MASS STREAM FUNCTION',
     5 'ATLANT SALT STREAM FUNC - .035*MASS',
     6 'PACIFC SALT STREAM FUNC - .035*MASS',
     7 'INDIAN SALT STREAM FUNC - .035*MASS',
     8 'GLOBAL SALT STREAM FUNC - .035*MASS'/

C**** Calculate latitudes
      JEQ = JM/2
      JDLAT = NINT(180./(JM-1))
      DO J=1,JM-1
        JLAT(J)  = JDLAT*(J-JEQ)
      END DO
      R4LEV = ZE
C****
C**** lat/lon diagnostics
C****
      WRITE (6,*)
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)

C**** Ocean Surface Height (cm)
c      IF(.not.QOSH)  GO TO 220
c      DO 211 I=1,IM*(JM-1)+1
c      Q(I,1) = SKIP
c      IF(FOCEAN(I,1).le..5)  GO TO 211
c      Q(I,1) = 1.E2*(AIJ29(I,1)/IDACC(1))
C****     *            -  (AIJ16(I,1)/IDACC(1))/916.6)
c  211 CONTINUE
c      CALL WRITED (NAME(1),XLABEL,OUTMON,OUTYR)

C**** Loop over layers
      DO L=1,LMO
C****
C**** Ocean Potential Temperature (C)
C****
      DO J=1,JM
      DO I=1,IM
        Q(I,J) = SKIP
        IF(FOCEAN(I,J).gt..5 .and. OIJL(I,J,L,IJL_MO).le.0.)  THEN
          GOS = OIJL(I,J,L,IJL_G0M) / (OIJL(I,J,L,IJL_MO)*DXYP(J))
          SOS = OIJL(I,J,L,IJL_S0M) / (OIJL(I,J,L,IJL_MO)*DXYP(J))
          Q(I,J) = TEMGS(GOS,SOS)
        END IF
      END DO
      END DO
c      IF (L.gt.1) WRITE (NAME(2)(8:15),'(A5,I3)') 'Level',L
c      CALL WRITED (NAME(2),XLABEL,OUTMON,OUTYR)
C****
C**** Ocean Salinity (per mil)
C****
      DO J=1,JM
      DO I=1,IM
        Q(I,J) = SKIP
        IF(FOCEAN(I,J).gt..5 .and. OIJL(I,J,L,IJL_MO).le.0.) 
     *       Q(I,J) = 1.E3*OIJL(I,J,L,IJL_S0M) / (OIJL(I,J,L,IJL_MO)
     *       *DXYP(J))
      END DO
      END DO
c      IF (L.gt.1) WRITE (NAME(3)(8:15),'(A5,I3)') 'Level',L
c      CALL WRITED (NAME(3),XLABEL,OUTMON,OUTYR)
      END DO
C****
C**** Ocean Column and Sea Ice Mass (kg/m**2)
C****
c  240 IF(.not.QOCM)  GO TO 300
c      IDACC5 = IDACC(5)
c      IF(IDACC5.le.0)  IDACC5 = IDACC(2)
c      DO 242 I=1,IM*(JM-1)+1
c      Q(I,1) = SKIP
c      IF(FOCEAN(I,1).le..5)  GO TO 242
c      Q(I,1) = 0.
c      DO 241 L=1,LMO
c  241 Q(I,1) = Q(I,1) + OIJL(I,1,L,IJL_MO)
c      Q(I,1) = Q(I,1)/IDACC5 + AIJ16(I,1)/IDACC(1)
c  242 CONTINUE
c      CALL WRITED (NAME(4),XLABEL,OUTMON,OUTYR)
C****
C**** East-West or North-South Mass Flux (10^9 kg/s)
C****
  300 DO K=2,3
c      IF(.not.QL(K))  GO TO 340
        Q = 0.
        DO I=1,IM*(JM-1)
        DO L=LMINMF,LMAXMF
          Q(I,1) = Q(I,1) + OIJL(I,1,L,K)
        END DO
        Q(I,1) = 2.E-9*Q(I,1) / (IDACC(1)*NDYNO)
      END DO
      IF(LMINMF.eq.LMAXMF) WRITE (NAMEL(K)(41:43),'(I3)') LMINMF
      IF(LMINMF.lt.LMAXMF) WRITE (NAMEL(K)(41:49),'(2I4)') LMINMF
     *     ,LMAXMF
c      CALL WRITED (NAMEL(K),XLABEL,OUTMON,OUTYR)
c  340 CONTINUE
      END DO
C****
C**** Vertical Mass Flux (10^-2 kg/s*m)
C****
      DO K=1,LMO
c     IF(KVMF(K).le.0)  GO TO 370
        L =KVMF(K)
        DO J=1,JM
          DO I=1,IMAXJ(J)
            Q(I,J) = 2.E2*OIJL(I,J,L,IJL_MFW) / (IDACC(1)*NDYNO*DXYP(J))
          END DO
        END DO
        WRITE (NAMEL(4)(43:45),'(I3)') L
c      CALL WRITED (NAMEL(4),XLABEL,OUTMON,OUTYR)
      END DO
C****
C**** East-West or North-South Heat Flux (10^15 W)
C****
      DO K=6,7
c      IF(.not.QL(K))  GO TO 440
        Q = 0.
        DO I=1,IM*(JM-1)
          DO L=LMINEF,LMAXEF
            Q(I,1) = Q(I,1) + OIJL(I,1,L,K)
          END DO
          Q(I,1) = 1.E-15*Q(I,1) / (IDACC(1)*DTS)
        END DO
        IF(LMINEF.eq.LMAXEF) WRITE (NAMEL(K)(39:41),'(I3)') LMINEF
        IF(LMINEF.lt.LMAXEF) WRITE (NAMEL(K)(39:47),'(2I4)') LMINEF
     *       ,LMAXEF
c      CALL WRITED (NAMEL(K),XLABEL,OUTMON,OUTYR)
      END DO
C****
C**** East-West or North-South Salt Flux (10^6 kg/s)
C****
      DO K=10,11
c      IF(.not.QL(K))  GO TO 540
        Q = 0.
        DO I=1,IM*(JM-1)
          DO L=LMINSF,LMAXSF
            Q(I,1) = Q(I,1) + OIJL(I,1,L,K)
          END DO
          Q(I,1) = 1.E-6*Q(I,1) / (IDACC(1)*DTS)
        END DO
        IF(LMINSF.eq.LMAXSF) WRITE (NAMEL(K)(41:43),'(I3)') LMINSF
        IF(LMINSF.lt.LMAXSF) WRITE (NAMEL(K)(41:49),'(2I4)') LMINSF
     *       ,LMAXSF
c     CALL WRITED (NAMEL(K),XLABEL,OUTMON,OUTYR)
      END DO
C****
C**** Gent Mcwilliams coefficient (10^-2 kg/s*m)
C****
      DO K=1,LMO
c      IF(KCMF(K).le.0)  GO TO 620
        L =KCMF(K)
        DO J=1,JM
          DO I=1,IMAXJ(J)
            Q(I,J) = 1.E2*OIJL(I,J,L,IJL_GGMFL)/(IDACC(1)*DTS*DXYP(J))
          END DO
        END DO
        WRITE (NAMEL(13)(44:46),'(I3)') L
c      CALL WRITED (NAMEL(13),XLABEL,OUTMON,OUTYR)
      END DO
C****
C**** Vertical Diffusion Coefficient (cm/s)
C****
      DO K=1,LMO
c      IF(KVDC(K).le.0)  GO TO 720
        L =KVDC(K)
        DO I=IM+1,IM*(JM-1)+1
          Q(I,1) = 1.E4*.00097*.00097*OIJL(I,1,L,IJL_KVM)/(IDACC(1)*4.)
        END DO
        WRITE (NAMEL(14)(45:47),'(I3)') L
c      CALL WRITED (NAMEL(14),XLABEL,OUTMON,OUTYR)
      END DO

      IF (QKV) THEN
        DO L=1,LMO
          DO I=1,IM*JM
            Q3(I,1,L)=1E4*.00097**2*OIJL(I,1,L,IJL_KVM)/(IDACC(1)*4.)
          END DO
        END DO
        TITLE = NAMEL(15)//'  Run '//XLABEL//'  '//OUTMON//OUTYR
c        WRITE(2) TITLE,Q3
c       DO L=1,LMO
c         DO I=1,IM*JM
c           Q3(I,1,L)=1.E4*.00097*.00097*OIJL(I,1,L,IJL_KVG)/(IDACC(1)*4.)
c         END DO
c       END DO
c       TITLE = NAMEL(16)//'  Run '//XLABEL//'  '//OUTMON//OUTYR
c       WRITE(2) TITLE,Q3
      END IF
C****
C**** Calculate Horizontal Mass Stream Function and write it
C****
      NAME(1)= 'HORIZONTAL MASS TRANSPORT STREAMFUNCTION (Sv)'
      NAME(2)= 'HORIZ. SALT STREAMFUNCTION (10**6 kg/s)-.035*MASS'
      FAC   = -2.E-9/(IDACC(1)*NDYNO)
      FACST = -1.E-9/(IDACC(1)*NDYNO*DTO)
      CALL STRMIJ (OIJL(1,1,1,IJL_MFU),FAC,OLNST(1,1,LN_MFLX),FACST
     *     ,SFIJM)
      TITLE = NAME(1)//'  Run '//XLABEL(1:6)
      WRITE(TITLE(63:80),'(A6,I4)') JMON0,JYEAR0
c      WRITE(2) TITLE,SFIJM
c      WRITE(6,'(A80)') TITLE

c      call ijmap (title,sfijm,skip,1)

C****
C**** Calculate Salt Stream Function and write it
C****
      FAC   = -1.E-6/(IDACC(1)*NDYNO*DTO)
      FACST = -.5E-6/(IDACC(1)*NDYNO*DTO)
      CALL STRMIJ (OIJL(1,1,1,IJL_SFLX),FAC,OLNST(1,1,LN_SFLX),FACST
     *     ,SFIJS)
C**** Subtract .035 times the Mass Stream Function
      DO J=1,JM
        DO I=1,IM
          IF (SFIJS(I,J).ne.SKIP) SFIJS(I,J) = SFIJS(I,J)-35.*SFIJM(I,J)
        END DO
      END DO
      TITLE = NAME(2)//'  Run '//XLABEL(1:6)
      WRITE(TITLE(63:80),'(A6,I4)') JMON0,JYEAR0
c      WRITE(2) TITLE,SFIJS
C**** Output Key diagnostics: Gulf Stream, ACC, Kuroshio
      WRITE(6,'(A17)') "Key horizontal mass stream function diags:"
      GSMAX=0
      GSMIN=100.
      DO J=30,33
        DO I=20,24
          IF (SFIJM(I,J).GT.GSMAX) GSMAX=SFIJM(I,J)
          IF (SFIJM(I,J).LT.GSMIN) GSMIN=SFIJM(I,J)
        END DO
      END DO
      WRITE(6,'(a,F10.3)') "Gulf Stream (MAX in (20-24,30-33)):", GSMAX
     *     -GSMIN
      CKMAX=0
      CKMIN=100.
      DO J=30,33
        DO I=63,67
          IF (SFIJM(I,J).GT.CKMAX) CKMAX=SFIJM(I,J)
          IF (SFIJM(I,J).LT.CKMIN) CKMIN=SFIJM(I,J)
        END DO
      END DO
      WRITE(6,'(a,F10.3)') "Kuroshio    (MAX in (63-67,30-33)):", CKMAX
     *     -CKMIN
      ACMAX=0
      ACMIN=100.
      I=23
      DO J=6,10
        IF (SFIJM(I,J).GT.ACMAX) ACMAX=SFIJM(I,J)
        IF (SFIJM(I,J).LT.ACMIN) ACMIN=SFIJM(I,J)
      END DO
      WRITE(6,'(a,F10.3)') "ACC               (Drakes Passage):", ACMAX
     *     -ACMIN
C****
C**** Calculate Mass Stream Function and write it
C****
      FAC   = -2.E-9/(IDACC(1)*NDYNO)
      FACST = -1.E-9/(IDACC(1)*NDYNO*DTO)
      CALL STRMJL (OIJL(1,1,1,IJL_MFU),FAC,OLNST(1,1,LN_MFLX),FACST,SFM)
      DO KB=1,4
c        CALL WRITED (KB,NAMESF(KB)  ,SFM(1,0,KB))
        sname='sf_by_basin'
        units='Sv'
        scalel=1.
        scalej=1.
C****
C**** Write data to PRinT file
C****
c      TITLE = TRIM(NAMESF(KB))//' ('//trim(units)//')'
      PRINT*,NAMESF(KB),units
      call sys_flush(6)
c      WRITE (6,920) TITLE
      WRITE (6,921) 'Southern Hemisphere',(JLAT(J),J=1,JEQ)
      DO 220 L=0,LMO
  220 WRITE (6,922) R4LEV(L),(NINT(SFM(J,L,KB)),J=1,JEQ)
      WRITE (6,921) 'Northern Hemisphere',(JLAT(J),J=JEQ,JM-1)
      DO 240 L=0,LMO
  240 WRITE (6,922) R4LEV(L),(NINT(SFM(J,L,KB)),J=JEQ,JM-1)
  920 FORMAT ('1',A80/)
  921 FORMAT ('0',A20 / '0',6X,23I5 / 8X,23('-----'))
  922 FORMAT (F6.0,2X,23I5)

      END DO
C****
C**** Calculate Salt Stream Function and write it
C****
      FAC   = -1.E-6/(IDACC(1)*NDYNO*DTO)
      FACST = -.5E-6/(IDACC(1)*NDYNO*DTO)
      CALL STRMJL(OIJL(1,1,1,IJL_SFLX),FAC,OLNST(1,1,LN_SFLX),FACST,SFS)
C**** Subtract .035 times the Mass Stream Function
      DO KB=1,4
        DO L=0,LMO-1
        DO J=1,JM-1
          IF (SFS(J,L,KB).ne.SKIP) SFS(J,L,KB) = SFS(J,L,KB) -
     *         35.*SFM(J,L,KB)
        END DO
        END DO
c        CALL WRITED (KB,NAMESF(KB+4),SFS(1,0,KB))
      END DO
C**** 
C**** vertical diagnostics
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
        WRITE(6,'(2X,I4,10F10.3)') L,OL(L,L_RHO)/(4.*NOL(L)*IDACC(1)),
     *       OL(L,L_TEMP)/(4.*NOL(L)*IDACC(1)),
     *       1d3*OL(L,L_SALT)/(4.*NOL(L)*IDACC(1))
      END DO
C****
C**** Strait diagnostics
C****
      WRITE (TITLE(50:80),'(A10,2X,A4,I4,8X)') XLABEL(1:10),AMON0,JYEAR0
      SUMORMN(:)=1
      NAME(LN_MFLX)='Strait Transport of Mass (10^6 kg/s)'
      SCALEO(LN_MFLX) = 1.D- 6 / DTS
      NAME(LN_GFLX)='Strait Trans of Potential Enthalpy (10^11 W)'
      SCALEO(LN_GFLX) = .5D-11 / DTS
      NAME(LN_SFLX)='Strait Transport of Salt (10^5 kg/s)'
      SCALEO(LN_SFLX) = .5D- 5 / DTS
      NAME(LN_KVM)= 'Vertical Diffusion Coefficient (m^2/s)'
      SCALEO(LN_KVM) = .5
      SUMORMN(LN_KVM)=2
      NAME(LN_ICFL)='Sea Ice Mass Flux (10^4 kg/s)'
      SCALEO(LN_ICFL) = .5D-4 / DTS
      SUMORMN(LN_ICFL)=0
      DO N=1,KOLNST
        if (n.eq.LN_MFLX .or. n.eq.LN_GFLX.or. n.eq.LN_SFLX.or. n.eq
     *       .LN_KVM.or. n.eq.LN_ICFL) THEN 
          AS = 0.
          TITLE(1:50) = NAME(N)
          DO NS=1,NMST
            DO L=1,LMST(NS)
              AS(L,NS) = OLNST(L,NS,N)*SCALEO(N)/IDACC(1)
            END DO
          END DO
          IF (N.eq.LN_ICFL) THEN
            LMSTMIN(:)=1
            CALL STABLE (LMSTMIN,name_st,AS,TITLE,SUMORMN(N))
          ELSE
            CALL STABLE (LMST,name_st,AS,TITLE,SUMORMN(N))
          END IF
          call sys_flush(6)
        END IF
      END DO
C****
      RETURN
      END SUBROUTINE diag_OCEAN

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
      END

      SUBROUTINE STRMJL (OIJL,FAC,OLNST,FACST,SF)
!@sum  STRMJL calculates the latitude by layer stream function
!@auth G. Russell/G. Schmidt
!@ver  1.0 (from OSFJL064)
C****
C**** Input: 
!@var OIJL  = west-east and south-north tracer fluxes (kg/s)
!@var FAC   = global scaling factor
!@var OLNST = strait mass flux (kg/s)
!@var FACST = global scaling factor for straits
C**** Output:  
!@var    SF = stream function (kg/s)
C****
c      USE CONSTANT, only : skip=>undef
      USE OCEAN, only : im,jm,lmo
      USE STRAITS, only : nmst
      USE ODIAG, only : kbasin
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: FAC,FACST
      REAL*8, INTENT(IN), DIMENSION(IM,JM,LMO,2) :: OIJL
      REAL*8, INTENT(IN), DIMENSION(LMO,NMST) :: OLNST
      REAL*4, INTENT(OUT), DIMENSION(JM-1,0:LMO,0:4) :: SF
      REAL*4, DIMENSION(4) :: SUMB
      REAL*8 :: SKIP = 0.
      INTEGER :: KB,I,J,L,K
C****
C**** Zero out the Stream Function at the ocean bottom
C****
      SF(:,:,:) = 0.
C****
C**** Integrate the Stream Function upwards from the ocean bottom
C****
      DO L=LMO-1,0,-1
      DO J=1,JM-1
        DO KB=1,4
          SF(J,L,KB) = SF(J,L+1,KB)
        END DO
        DO I=1,IM
          KB = KBASIN(I,J+1)
          SF(J,L,KB) = SF(J,L,KB) + OIJL(I,J,L+1,2)*FAC
        END DO
      END DO
C****
C**** Add south-north flow between basins
C****
C**** Add Bering Strait flow into Atlantic from Pacific
C      DO 310 J=1,38
C      SF(J,L,1) = SF(J,L,1) + OIJL(3,39,L+1,2)*FAC
C  310 SF(J,L,2) = SF(J,L,2) - OIJL(3,39,L+1,2)*FAC
C**** Add Indonesian Straits flow into Pacific from Indian
C      DO 320 J=1,20
C      SF(J,L,2) = SF(J,L,2) + OIJL(59,21,L+1,2)*FAC
C      SF(J,L,2) = SF(J,L,2) + OIJL(60,21,L+1,2)*FAC
C      SF(J,L,3) = SF(J,L,3) - OIJL(59,21,L+1,2)*FAC
C  320 SF(J,L,3) = SF(J,L,3) - OIJL(60,21,L+1,2)*FAC
C      DO 330 J=1,19
C      SF(J,L,2) = SF(J,L,2) + OIJL(62,20,L+1,2)*FAC
C  330 SF(J,L,3) = SF(J,L,3) - OIJL(62,20,L+1,2)*FAC
C****
C**** Add west-east flow between basins
C****
C**** Add Drake Passage flow into Atlantic from Pacific
C      FLOW = 0.
C      DO 410 J=8,1,-1
C      FLOW = FLOW + OIJL(23,J+1,L+1,1)*FAC
C      SF(J,L,1) = SF(J,L,1) + FLOW
C  410 SF(J,L,2) = SF(J,L,2) - FLOW
C**** Add flow south of Africa into Indian from Atlantic
C      FLOW = 0.
C      DO 420 J=13,1,-1
C      FLOW = FLOW + OIJL(40,J+1,L+1,1)*FAC
C      SF(J,L,3) = SF(J,L,3) + FLOW
C  420 SF(J,L,1) = SF(J,L,1) - FLOW
C**** Add flow south of Australia into Pacific from Indian
C      FLOW = 0.
C      DO 430 J=12,1,-1
C      FLOW = FLOW + OIJL(65,J+1,L+1,1)*FAC
C      SF(J,L,2) = SF(J,L,2) + FLOW
C  430 SF(J,L,3) = SF(J,L,3) - FLOW
C****
C**** Add strait flow from to the Stream Function
C****
C**** Fury & Hecla: (19,42) to (20,40)
      SF(41,L,1) = SF(41,L,1) - OLNST(L+1,1)*FACST
      SF(40,L,1) = SF(40,L,1) - OLNST(L+1,1)*FACST
C**** Nares: (22,43) to (24,44)
      SF(43,L,1) = SF(43,L,1) + OLNST(L+1,2)*FACST
C**** Gibralter: (35,32) to (37,33)
      SF(32,L,1) = SF(32,L,1) + OLNST(L+1,3)*FACST
C**** English: (36,36) to (37,37)
      SF(36,L,1) = SF(36,L,1) + OLNST(L+1,4)*FACST
C**** Bosporous: (42,33) to (43,34)
      SF(33,L,1) = SF(33,L,1) + OLNST(L+1,6)*FACST
C**** Red Sea: (44,29) to (45,28)
      SF(28,L,3) = SF(28,L,3) - OLNST(L+1,7)*FACST
C**** Bab-al-Mandab: (45,28) to (46,27)
      SF(27,L,3) = SF(27,L,3) - OLNST(L+1,8)*FACST
C**** Hormuz: (47,30) to (49,29)
      SF(29,L,3) = SF(29,L,3) - OLNST(L+1,9)*FACST
C**** Korea: (62,32) to (63,33)
      SF(32,L,2) = SF(32,L,2) + OLNST(L+1,11)*FACST
C**** Soya: (64,34) to (65,35)
      SF(34,L,2) = SF(34,L,2) + OLNST(L+1,12)*FACST
C**** Malacca: (56,25) to (58,24), from Indian Ocean to Pacific Ocean
      SF(24,L,2) = SF(24,L,2) - OLNST(L+1,10)*FACST  !*.5
c      DO 510 J=1,23
c      SF(J,L,2) = SF(J,L,2) - OLNST(L+1,10)*FACST*.5
c  510 SF(J,L,3) = SF(J,L,3) + OLNST(L+1,10)*FACST*.5
      END DO
C****
C**** Calculate global Stream Function by summing it over 3 oceans
C****
      SF(:,:,4) = SF(:,:,1) + SF(:,:,2) + SF(:,:,3) + SF(:,:,4)
C**** 
C**** Mask streamfunction so that topography is put to skip value
C****
      DO J=1,JM
        SUMB = 0
        DO I=1,IM
          SUMB(KBASIN(I,J))=1.
        END DO
        SUMB(4)=SUMB(1)+SUMB(2)+SUMB(3)+SUMB(4)
        DO K=1,4
          IF (SUMB(K).eq.0) THEN
            SF(J,0:LMO,K) = SKIP   ! SKIP areas of no integral
          ELSE
            L=0
  620       L=L+1       ! SKIP below bottom values
            IF (SF(J,L,K).eq.0 .and. SF(J,L-1,K).eq.0) THEN
              SF(J,L:LMO,K) = SKIP
            ELSE
              GOTO 620
            END IF
          END IF
        END DO
      END DO
C****
      RETURN
      END

      SUBROUTINE OBASIN
!@sum  OBASIN Read in KBASIN: 0=continent, 1=Atlantic, 2=Pacific, 3=Indian
!@auth G. Russell
!@ver  1.0
      USE OCEAN, only : IM,JM
      USE ODIAG, only : kbasin
      USE FILEMANAGER
      IMPLICIT NONE
      CHARACTER TITLE*72, CBASIN(IM,JM)
      INTEGER J,I,iu_KB
C****
C**** read in basin data
c      OPEN  (3,FILE='/u/gavin/gissgcm/data/KB4X512.OCN')
      call openunit('KBASIN',iu_KB,.false.,.true.)
      
      READ  (iu_KB,900) TITLE
      WRITE (6,*) 'Read on unit ',iu_KB,': ',TITLE
      READ  (iu_KB,900)
      DO J=JM,1,-1
        READ  (iu_KB,901) (CBASIN(I,J),I=1,72)
      END DO
      call closeunit(iu_KB)

      DO J=1,JM
      DO I=1,IM
        SELECT CASE (CBASIN(I,J))
        CASE DEFAULT
          KBASIN(I,J) = 0
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


      SUBROUTINE STRMIJ (OIJL,FAC,OLNST,FACST,SF)
!@sum STRMIJ calculates the latitude by longitude stream function for a
!@+   given quantity.
C****
C**** Input: 
!@var OIJL  = west-east and south-north tracer fluxes (kg/s)
!@var   FAC = global scaling factor
!@var OLNST = strait mass flux (kg/s)
!@var FACST = global scaling factor for straits
C**** Output:  
!@var    SF = stream function (kg/s)
C****
      USE CONSTANT, only : skip=>undef
      USE OCEAN, only : im,jm,lmo
      USE STRAITS, only : nmst
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: FAC,FACST
      REAL*8, INTENT(IN), DIMENSION(IM,JM,LMO,2) :: OIJL
      REAL*8, INTENT(IN), DIMENSION(LMO,NMST) :: OLNST
      REAL*4, INTENT(OUT), DIMENSION(IM,JM) :: SF
      REAL*4, DIMENSION(4) :: SUMB
      INTEGER :: KB,I,J,L
      REAL*4 TSUM
C****
C**** Integrate up from South Pole
C****
      SF=0
      DO J=2,JM-1
        DO I=1,IM
          SF(I,J) = SF(I,J-1)
          DO L=1,LMO
            SF(I,J) = SF(I,J) + OIJL(I,J,L,1)*FAC
          END DO
        END DO
C****
C**** Add strait flow from to the Stream Function
C****
      DO L=1,LMO
C**** Fury & Hecla: (19,42) to (20,40)
      IF (J.eq.41) SF(19,41) = SF(19,41) + OLNST(L,1)*FACST
C**** Nares: (22,43) to (24,44)
      IF (J.eq.44) SF(22,44) = SF(22,44) + OLNST(L,2)*FACST
      IF (J.eq.44) SF(23,44) = SF(23,44) + OLNST(L,2)*FACST
C**** Gibrater: (35,32) to (37,33)
      IF (J.eq.33) SF(35,33) = SF(35,33) + OLNST(L,3)*FACST
      IF (J.eq.33) SF(36,33) = SF(36,33) + OLNST(L,3)*FACST
C**** Engish: (36,36) to (37,37)
      IF (J.eq.37) SF(36,37) = SF(36,37) + OLNST(L,4)*FACST
C**** Bosporous: (42,33) to (43,34)
      IF (J.eq.34) SF(42,34) = SF(42,34) + OLNST(L,6)*FACST
C**** Red Sea: (44,29) to (45,28)
      IF (J.eq.29) SF(44,29) = SF(44,29) + OLNST(L,7)*FACST
C**** Bab-al-Mandab: (45,28) to (46,27)
      IF (J.eq.28) SF(45,28) = SF(45,28) + OLNST(L,8)*FACST
C**** Hormuz: (47,30) to (49,29)
      IF (J.eq.30) SF(47,30) = SF(47,30) + OLNST(L,9)*FACST
      IF (J.eq.30) SF(48,30) = SF(48,30) + OLNST(L,9)*FACST
C**** Korea: (62,32) to (63,33)
      IF (J.eq.33) SF(62,33) = SF(62,33) + OLNST(L,11)*FACST
C**** Soya: (64,34) to (65,35)
      IF (J.eq.35) SF(64,35) = SF(64,35) + OLNST(L,12)*FACST
C**** Malacca: (56,25) to (58,24), 
      IF (J.eq.25) SF(56,25) = SF(56,25) + OLNST(L,10)*FACST
      IF (J.eq.25) SF(57,25) = SF(57,25) + OLNST(L,10)*FACST
      END DO
      END DO
C**** Correct SF for mean E-W drift (SF over topography --> 0)
      TSUM=0 
      DO I=14,17   ! (mid N. America as example)
        TSUM=TSUM+SF(I,34)
      END DO
      TSUM=TSUM/4.
      DO J=1,JM-1
        DO I=1,IM
          IF (SF(I,J).ne.SKIP) SF(I,J)=SF(I,J)-TSUM
        END DO
      END DO
      RETURN
      END
