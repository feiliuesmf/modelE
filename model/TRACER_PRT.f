!@sum TRACER_PRT Tracer diagnostics.  These routines are generic for
!@+    all tracers
!@+   TRACEA and DIAGTCA are called throughout the day
!@+   The other routines are for printing
!@auth Jean Lerner (with ideas stolen from G.Schmidt, R. Ruedy, etc.)

      SUBROUTINE TRACEA
!@sum TRACEA accumulates tracer concentration diagnostics (IJL, JL)
!@auth J.Lerner
!@ver  1.0
      USE MODEL_COM, only: im,jm,lm,itime
      USE GEOM, only: imaxj
      USE SOMTQ_COM, only: mz  
      USE TRACER_COM 
      USE TRACER_DIAG_COM 
      USE DYNAMICS, only: am,byam 
      implicit none

      integer i,j,l,k,n
      real*8 tsum,asum

C**** Compute AM: kg air/sq meter
      call calc_am(lm)
C****
C**** Accumulate concentraton for all tracers
C****
      do 600 n=1,ntm
      IF (itime.lt.itime_tr0(n)) cycle
C**** Latitude-longitude by layer concentration
      do l=1,lm
        taijln(:,:,l,n) = taijln(:,:,l,n) + trm(:,:,l,n)*byam(l,:,:)
      end do
C**** Average concentration; surface concentration; total mass
      l = 1
      do j=1,jm
      do i=1,im
        tsum = sum(trm(i,j,:,n))  !sum over l
        asum = sum(am(:,i,j))     !sum over l
        taijkn(i,j,ijkn_mass,n) = taijkn(i,j,ijkn_mass,n)+tsum  !MASS 
        taijkn(i,j,ijkn_conc,n) = taijkn(i,j,ijkn_conc,n)+tsum/asum
        taijkn(i,j,ijkn_surf,n) = taijkn(i,j,ijkn_surf,n)
     *     + max((trm(i,j,l,n)-trmom(mz,i,j,l,n))*byam(l,i,j) ,0.)
      enddo; enddo
C**** Zonal mean concentration
      do l=1,lm
      do j=1,jm
        tsum = sum(trm(1:imaxj(j),j,l,n)) !sum over i
        asum = sum(am(l,1:imaxj(j),j))    !sum over i
        tajln(j,l,jlnt_conc,n) = tajln(j,l,jlnt_conc,n)+tsum/asum
      enddo; enddo
  600 continue
      return
      end SUBROUTINE TRACEA


      SUBROUTINE DIAGTCA (M,NT)
!@sum  DIAGCAT Keeps track of the conservation properties of tracers
!@auth Gary Russell/Gavin Schmidt/Jean Lerner
!@ver  1.0
      USE MODEL_COM, only: im,jm,lm,fim
      USE TRACER_COM
      USE TRACER_DIAG_COM, only: tconsrv,nofmt,title_tcon
      USE GEOM, only: imaxj
      IMPLICIT NONE

!@var M index denoting which process changed the tracer
      INTEGER, INTENT(IN) :: m
!@var NT index denoting tracer number
      INTEGER, INTENT(IN) :: nt
!@var TOTAL amount of conserved quantity at this time
      REAL*8, DIMENSION(JM) :: total
      INTEGER :: i,j,l,nm,ni
      REAL*8 sstm,stm,mnow
C****
C**** THE PARAMETER M INDICATES WHEN DIAGCA IS BEING CALLED
C**** M=1,2...12:  See DIAGCA in DIAG.f
C****   13+ AFTER Sources and Sinks
C****
C**** NOFMT contains the indexes of the TCONSRV array where each
C**** change is to be stored for each quantity. If NOFMT(M,NT)=0,
C**** no calculation is done.
C**** NOFMT(1,NT) is the index for the instantaneous value.

      if (nofmt(m,nt).gt.0) then
C**** Calculate current value TOTAL
        do j=1,jm
          sstm = 0.
          do l=1,lm
            stm = 0.
            do i=1,imaxj(j)
              stm = stm+trm(i,j,l,nt)
            end do
            sstm = sstm+stm
          end do
          total(j) = sstm
        end do
        total(1) = fim*total(1)
        total(jm)= fim*total(jm)
        nm=nofmt(m,nt)
        ni=nofmt(1,nt)
c**** Accumulate difference from last time in TCONSRV(NM)
        if (m.gt.1) then
          tconsrv(:,nm,nt) = 
     &    tconsrv(:,nm,nt)+(total(:)-tconsrv(:,ni,nt))  !do 1,jm
        end if
C**** Save current value in TCONSRV(NI)
        tconsrv(:,ni,nt)=total(:)   !do 1,jm
      end if
      return
      end subroutine diagtca


      SUBROUTINE DIAGTCP
!@sum  DIAGCP produces tables of the conservation diagnostics
!@auth Gary Russell/Gavin Schmidt
!@ver  1.0
      USE CONSTANT, only: teeny, twopi
      USE MODEL_COM, only:
     &     jm,fim,idacc,jhour,jhour0,jdate,jdate0,amon,amon0,
     &     jyear,jyear0,nday,jeq,itime,itime0,xlabel
      USE GEOM, only:
     &     areag,dlon,dxyp,LAT_DG
      USE TRACER_COM, only: ntm ,itime_tr0
      USE tracer_DIAG_COM, only:
     &     TCONSRV,ktcon,scale_TCON,title_TCON,nsum_TCON,ia_TCON,nofmt
      USE DAGCOM, only: inc=>incj,xwon,kdiag
      IMPLICIT NONE

      INTEGER, DIMENSION(JM) :: MAREA
      DOUBLE PRECISION, DIMENSION(KTCON) :: FGLOB
      DOUBLE PRECISION, DIMENSION(2,KTCON) :: FHEM
      DOUBLE PRECISION, DIMENSION(JM+3,KTCON,NTM) :: CNSLAT
      CHARACTER*4, PARAMETER :: HEMIS(2) = (/' SH ',' NH '/),
     *     DASH = ('----')
      INTEGER :: j,jhemi,jnh,jp1,jpm,jsh,jx,n,k,KTCON_max
      DOUBLE PRECISION :: aglob,ahem,feq,fnh,fsh,days

      if (kdiag(8).ge.2) return
C**** CALCULATE SCALING FACTORS
      IF (IDACC(12).LT.1) IDACC(12)=1
C**** Calculate areas
      DO J=1,JM
        MAREA(J)=1.D-10*XWON*FIM*DXYP(J)+.5
      END DO
      AGLOB=1.D-10*AREAG*XWON
      AHEM=1.D-10*(.5*AREAG)*XWON
C****
C**** Outer loop over tracers
C****
      DO 900 N=1,NTM
      IF (itime.LT.itime_tr0(N)) cycle
C**** CALCULATE SUM OF CHANGES
C**** LOOP BACKWARDS SO THAT INITIALIZATION IS DONE BEFORE SUMMATION!
      DO J=1,JM
        DO K=KTCON,1,-1
          IF (NSUM_TCON(K,N).eq.0) THEN
            TCONSRV(J,K,N)=0.
          ELSEIF (NSUM_TCON(K,N).gt.0) THEN
            TCONSRV(J,NSUM_TCON(K,N),N)=
     *      TCONSRV(J,NSUM_TCON(K,N),N)+TCONSRV(J,K,N)
     *           *SCALE_TCON(K,N)*IDACC(12)/(IDACC(IA_TCON(K))+teeny)
          KTCON_max = NSUM_TCON(K,N)
          END IF
        END DO
      END DO
C**** CALCULATE ALL OTHER CONSERVED QUANTITIES ON TRACER GRID
      DO K=1,KTCON_max
        FGLOB(K)=0.
        FHEM(1,K)=0.
        FHEM(2,K)=0.
        DO JSH=1,JEQ-1
          JNH=1+JM-JSH
          FSH=TCONSRV(JSH,K,N)*SCALE_TCON(K,N)/(IDACC(IA_TCON(K))+teeny)
          FNH=TCONSRV(JNH,K,N)*SCALE_TCON(K,N)/(IDACC(IA_TCON(K))+teeny)
          FGLOB (K)=FGLOB (K)+(FSH+FNH)
          FHEM(1,K)=FHEM(1,K)+FSH
          FHEM(2,K)=FHEM(2,K)+FNH
          CNSLAT(JSH,K,N)=FSH/(FIM*DXYP(JSH))
          CNSLAT(JNH,K,N)=FNH/(FIM*DXYP(JNH))
        END DO
        FGLOB (K)=FGLOB (K)/AREAG
        FHEM(1,K)=FHEM(1,K)/(.5*AREAG)
        FHEM(2,K)=FHEM(2,K)/(.5*AREAG)
        CNSLAT(JM+1,K,N)=FHEM(1,K)
        CNSLAT(JM+2,K,N)=FHEM(2,K)
        CNSLAT(JM+3,K,N)=FGLOB (K)
      END DO
C**** LOOP OVER HEMISPHERES
      DAYS=(Itime-Itime0)/DFLOAT(nday)
      DO JHEMI=2,1,-1
        WRITE (6,'(''1'',A)') XLABEL
        WRITE (6,902) JYEAR0,AMON0,JDATE0,JHOUR0,
     *       JYEAR,AMON,JDATE,JHOUR,ITIME,DAYS
        JP1=1+(JHEMI-1)*(JEQ-1)
        JPM=JHEMI*(JEQ-1)
C**** PRODUCE TABLES 
        WRITE (6,'(''0'')')
        WRITE (6,903) (DASH,J=JP1,JPM,INC)
        WRITE (6,904) HEMIS(JHEMI),(NINT(LAT_DG(JX,1)),JX=JPM,JP1,-INC)
        WRITE (6,903) (DASH,J=JP1,JPM,INC)
        DO K=1,KTCON_max
        WRITE (6,905) TITLE_TCON(K,N),FGLOB(K),FHEM(JHEMI,K),
     *         (NINT(CNSLAT(JX,K,N)),JX=JPM,JP1,-INC)
        END DO
        WRITE (6,906) AGLOB,AHEM,(MAREA(JX),JX=JPM,JP1,-INC)
      END DO
  900 CONTINUE
      RETURN
C****
  902 FORMAT ('0Conservation Quantities       From:',
     *  I6,A6,I2,',  Hr',I3,  6X,  'To:',I6,A6,I2,', Hr',I3,
     *  '  Model-Time:',I9,5X,'Dif:',F7.2,' Days')
  903 FORMAT (1X,28('--'),13(A4,'--'))
  904 FORMAT (41X,'GLOBAL',A8,2X,13I6)
  905 FORMAT (A38,2F9.2,1X,13I6)
  906 FORMAT ('0AREA (10^10 M^2)',F30.1,F9.1,1X,13I6)
      END SUBROUTINE DIAGTCP


      MODULE BDjlt
!@sum  stores info for outputting lat-sig/pressure diags for tracers
!@auth J. Lerner
      IMPLICIT NONE
      SAVE
!@param names of derived JLt/JLs tracer output fields
      INTEGER jlnt_nt_eddy,jlnt_vt_eddy  
      END MODULE BDjlt


      SUBROUTINE JLt_TITLEX
!@sum JLt_TITLEX sets up titles, etc. for composit JL output (tracers)
!@auth J. Lerner
!@ver  1.0
      USE TRACER_COM
      USE TRACER_DIAG_COM
      USE BDjlt
      IMPLICIT NONE
      character*50 :: unit_string
      integer k,n,kpmax

!@var jlnt_xx Names for TAJL diagnostics
      do n=1,ntm
      k = ktajl +1
        jlnt_nt_eddy = k
        sname_jln(k,n) = 'tr_nt_eddy_'//trname(n)
        lname_jln(k,n) = 'NORTHWARD TRANS. OF '//
     &     trim(trname(n))//' MASS BY EDDIES'
        jlq_power(k) = 10.
        jgrid_jlq(k) = 2
      k = k+1
        jlnt_vt_eddy = k
        sname_jln(k,n) = 'tr_vt_eddy_'//trname(n)
        lname_jln(k,n) = 'VERTICAL TRANS. OF '//
     &     trim(trname(n))//' MASS BY EDDIES'
        jlq_power(k) = 10.
      end do

      kpmax = k

      if (k .gt. ktajlx) then
        write (6,*) 'JLt_TITLEX: Increase ktajlx=',ktajlx,
     &         ' to at least ',k
        stop 'ktajlx too small'
      end if
C**** Construct UNITS string for output
      do n=1,ntm
      do k=ktajl+1,kpmax
      units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),' kg/s')
      end do; end do

      END SUBROUTINE JLt_TITLEX


      SUBROUTINE DIAGJLT
!@sum DIAGJLT controls calls to the pressure/lat print routine for 
!@+      tracers
!@auth J. Lerner
!@calls open_jl, JLt_TITLEX, JLMAP_t
C****
C**** THIS ROUTINE PRODUCES LATITUDE BY LAYER TABLES OF TRACERS
C****
      USE GEOM, only: BYDXYP,dxyp,LAT_DG
      USE DAGCOM, only: linect,plm,acc_period,qdiag,lm_req
      USE TRACER_COM
      USE TRACER_DIAG_COM
      USE MODEL_COM, only: jm,lm,fim,itime,idacc,xlabel,LRUNID,psfmpt
     &   ,sige,ptop
      USE BDJLT
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(JM) :: ONESPO
      DOUBLE PRECISION, DIMENSION(JM+LM) :: ONES
      DOUBLE PRECISION, DIMENSION(JM,LM) :: A
      DOUBLE PRECISION, DIMENSION(LM) :: PM
      DOUBLE PRECISION :: scalet
      INTEGER :: J,L,N,K,jtpow

C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      IF(QDIAG) call open_jl(trim(acc_period)//'.jlt'//XLABEL(1:LRUNID)
     *  ,jm,lm,0,lat_dg)
!?   *  ,jm,lm,lm_req,lat_dg)
C****
C**** INITIALIZE CERTAIN QUANTITIES
C****
      call JLt_TITLEX

      do l=1,lm
        pm(l)=psfmpt*sige(l+1)+ptop
      end do
      onespo(:)  = 1.d0
      onespo(1)  = fim
      onespo(jm) = fim
      ones(:) = 1.d0
      linect = 65
C****
C**** LOOP OVER TRACE GASES
C****
      DO 400 N=1,NTM
      IF (itime.LT.itime_tr0(N)) cycle
C****
C**** TRACER CONCENTRATION
C****
      k = jlnt_conc
      scalet = scale_jln(n)*scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *  plm,tajln(1,1,k,n),scalet,bydxyp,ones,lm,2,jgrid_jlq(k))
C****
C**** NORTHWARD TRANSPORTS: Total and eddies
C****
      a(:,:) = 0.
      k = jlnt_nt_tot
      do 205 l=1,lm
      do 205 j=1,jm-1
  205 a(j+1,l) = tajln(j,l,k,n)
      scalet = scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *  plm,a,scalet,ones,ones,lm,1,jgrid_jlq(k))
      do 210 l=1,lm
      do 210 j=1,jm-1
  210 a(j+1,l) = tajln(j,l,k,n) -tajln(j,l,jlnt_nt_mm,n)
      k = jlnt_nt_eddy
      scalet = scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *  plm,a,scalet,ones,ones,lm,1,jgrid_jlq(k))
C****
C**** VERTICAL TRANSPORTS: Total and eddies
C****
      k = jlnt_vt_tot
      scalet = scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *  pm,tajln(1,1,k,n),scalet,onespo,ones,lm-1,1,jgrid_jlq(k))
      a(:,:) = 0.
      do 260 l=1,lm-1
      do 260 j=1,jm
  260 a(j,l) = tajln(j,l,k,n)-tajln(j,l,jlnt_vt_mm,n)
      k = jlnt_vt_eddy
      scalet = scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *  pm,a,scalet,onespo,ones,lm-1,1,jgrid_jlq(k))
C****
C**** Convective Processes
C****
C**** MOIST CONVECTION 
      k = jlnt_mc
      scalet = scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *  plm,tajln(1,1,k,n),scalet,ones,ones,lm,1,Jgrid_jlq(k))
C**** TURBULENCE (or Dry Convection)
      k = jlnt_turb
      scalet = scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *  plm,tajln(1,1,k,n),scalet,ones,ones,lm,1,jgrid_jlq(k))
C**** LARGE-SCALE CONDENSATION
      k = jlnt_lscond
      scalet = scale_jlq(k)/idacc(ia_jlq(k))
      jtpow = ntm_power(n)+jlq_power(k)
      scalet = scalet*10.**(-jtpow)
      CALL JLMAP_t (lname_jln(k,n),sname_jln(k,n),units_jln(k,n),
     *  plm,tajln(1,1,k,n),scalet,ones,ones,lm,1,jgrid_jlq(k))
  400 CONTINUE
C****
C**** Sources and sinks
C****
      do k=1,ktajls
        n = jls_source(k)
        IF (n.eq.0 .or. itime.LT.itime_tr0(n))  cycle
        scalet = scale_jls(k)/idacc(ia_jls(k))
        scalet = scalet*10.**(-jls_power(k))
        CALL JLMAP_t (lname_jls(k),sname_jls(k),units_jls(k),
     *    plm,tajls(1,1,k),scalet,ones,ones,jls_ltop(k),1,jgrid_jls(k))
      end do

      if (qdiag) call close_jl
      RETURN
      END SUBROUTINE DIAGJLT


      SUBROUTINE JLMAP_t (LNAME,SNAME,UNITS,
     &     PL,AX,SCALET,SCALEJ,SCALEL,LMAX,JWT,JG)
C****
C**** THIS ROUTINE PRODUCES LAYER BY LATITUDE TABLES ON THE LINE
C**** PRINTER.  THE INTERIOR NUMBERS OF THE TABLE ARE CALCULATED AS
C****               AX * SCALET * SCALEJ * SCALEL.
C**** WHEN JWT=1, THE INSIDE NUMBERS ARE NOT AREA WEIGHTED AND THE
C****    HEMISPHERIC AND GLOBAL NUMBERS ARE SUMMATIONS.
C**** WHEN JWT=2, ALL NUMBERS ARE PER UNIT AREA.
C**** JG INDICATES PRIMARY OR SECONDARY GRID.
C**** THE BOTTOM LINE IS CALCULATED AS THE SUMMATION OF DSIG TIMES THE
C**** NUMBERS ABOVE 
C****
      USE DAGCOM, only: QDIAG,acc_period,inc=>incj,linect,jmby2,LM_REQ
      USE RESOLUTION, only: dsig,sige
      USE MODEL_COM, only:
     &     jm,lm,JDATE,JDATE0,AMON,AMON0,JYEAR,JYEAR0,XLABEL
      USE GEOM, only:
     &     WTJ,JRANGE_HEMI,LAT_DG
      IMPLICIT NONE

!@var units string containing output field units
      CHARACTER(LEN=50) :: UNITS
!@var lname string describing output field
      CHARACTER(LEN=50) :: LNAME
!@var sname string referencing output field
      CHARACTER(LEN=30) :: SNAME
!@var title string, formed as concatentation of lname//units
      CHARACTER(LEN=64) :: TITLE

      INTEGER :: JG,JWT,LMAX
      DOUBLE PRECISION :: SCALET
      DOUBLE PRECISION, DIMENSION(JM,LM) :: AX
      DOUBLE PRECISION, DIMENSION(JM) :: SCALEJ
      DOUBLE PRECISION, DIMENSION(LM) :: SCALEL,PL

      CHARACTER*4 DASH,WORD(4)
      DATA DASH/'----'/,WORD/'SUM','MEAN',' ','.1*'/

      INTEGER, DIMENSION(JM) :: MLAT
      DOUBLE PRECISION, DIMENSION(JM) :: ASUM
      DOUBLE PRECISION, DIMENSION(2) :: FHEM,HSUM
      INTEGER :: IWORD,J,JHEMI,K,L
      DOUBLE PRECISION :: FGLOB,FLATJ,GSUM,SDSIG

!?    REAL*8, DIMENSION(JM+3,LM+LM_REQ+1) :: XJL ! for binary output
      REAL*8, DIMENSION(JM+3,LM+1) :: XJL ! for binary output
      CHARACTER XLB*16,CLAT*16,CPRES*16,CBLANK*16,TITLEO*80
      DATA CLAT/'LATITUDE'/,CPRES/'PRESSURE (MB)'/,CBLANK/' '/

C form title string
      title = trim(lname)//' ('//trim(units)//')'

C****
C**** WRITE XLABEL ON THE TOP OF EACH OUTPUT PAGE
C****
      LINECT = LINECT+LMAX+7
      IF(LINECT.LE.60) GO TO 20
      WRITE (6,907) XLABEL(1:105),JDATE0,AMON0,JYEAR0,JDATE,AMON,JYEAR
      LINECT = LMAX+8
   20 CONTINUE
      WRITE (6,901) TITLE,(DASH,J=JG,JM,INC)
      WRITE (6,904) WORD(JWT),(NINT(LAT_DG(J,JG)),J=JM,JG,-INC)
      WRITE (6,905) (DASH,J=JG,JM,INC)
C****
C**** CALCULATE TABLE NUMBERS AND WRITE THEM TO THE LINE PRINTER
C****
         XJL(:,:) = -1.E30
      SDSIG = 1.-SIGE(LMAX+1)
      ASUM(:) = 0.
      HSUM(:) = 0.
      GSUM = 0.
      DO 240 L=LMAX,1,-1
      FGLOB = 0.
      DO 230 JHEMI=1,2
      FHEM(JHEMI) = 0.
      DO 220 J=JRANGE_HEMI(1,JHEMI,JG),JRANGE_HEMI(2,JHEMI,JG)
      FLATJ = AX(J,L)*SCALET*SCALEJ(J)*SCALEL(L)
         XJL(J,L) = FLATJ
      MLAT(J) = NINT(FLATJ)
      IF (JWT.EQ.1) THEN     
        ASUM(J) = ASUM(J)+FLATJ  !!!!most
      ELSE
        ASUM(J) = ASUM(J)+FLATJ*DSIG(L)/SDSIG  !!!! concentration
      ENDIF
      FHEM(JHEMI) = FHEM(JHEMI)+FLATJ*WTJ(J,JWT,JG)
  220 CONTINUE
  230 FGLOB = FGLOB+FHEM(JHEMI)/JWT
         XJL(JM+3,L)=FHEM(1)   ! SOUTHERN HEM
         XJL(JM+2,L)=FHEM(2)   ! NORTHERN HEM
         XJL(JM+1,L)=FGLOB     ! GLOBAL
      WRITE (6,902) PL(L),FGLOB,FHEM(2),FHEM(1),(MLAT(J),J=JM,JG,-INC)
      IF (JWT.EQ.1) THEN
        HSUM(1) = HSUM(1)+FHEM(1)
        HSUM(2) = HSUM(2)+FHEM(2)
        GSUM    = GSUM   +FGLOB
      ELSE
        HSUM(1) = HSUM(1)+FHEM(1)*DSIG(L)/SDSIG
        HSUM(2) = HSUM(2)+FHEM(2)*DSIG(L)/SDSIG
        GSUM    = GSUM   +FGLOB  *DSIG(L)/SDSIG
      ENDIF
  240 CONTINUE
      WRITE (6,905) (DASH,J=JG,JM,INC)
      ASUM(jmby2+1) = ASUM(jmby2+1)/JG
         DO 180 J=JG,JM
! 180    XJL(J   ,LM+LM_REQ+1)=ASUM(J)
!        XJL(JM+3,LM+LM_REQ+1)=HSUM(1)   ! SOUTHERN HEM
!        XJL(JM+2,LM+LM_REQ+1)=HSUM(2)   ! NORTHERN HEM
!        XJL(JM+1,LM+LM_REQ+1)=GSUM      ! GLOBAL
  180    XJL(J   ,LM+1)=ASUM(J)
         XJL(JM+3,LM+1)=HSUM(1)   ! SOUTHERN HEM
         XJL(JM+2,LM+1)=HSUM(2)   ! NORTHERN HEM
         XJL(JM+1,LM+1)=GSUM      ! GLOBAL
         XLB=' '//acc_period(1:3)//' '//acc_period(4:12)//'  '
         TITLEO=TITLE//XLB
         IF(QDIAG) CALL POUT_JL(TITLEO,LNAME,SNAME,UNITS,
     *        JG,LMAX,XJL,PL,CLAT,CPRES)
      IF (LMAX.EQ.1) THEN
         LINECT = LINECT-1
         RETURN
      ENDIF
      IF (LMAX.EQ.1) RETURN
      WRITE (6,903) WORD(JWT),GSUM,HSUM(2),HSUM(1),
     *    (NINT(ASUM(J)),J=JM,JG,-INC)
      RETURN
C****
  901 FORMAT ('0',30X,A64/2X,32('-'),24A4)
  902 FORMAT (1X,F8.3,3F8.1,1X,24I4)
  903 FORMAT (1X,A6,2X,3F8.1,1X,24I4)
  904 FORMAT ('  P(MB)    ',A4,' G     NH      SH  ',24I4)
  905 FORMAT (2X,32('-'),24A4)
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
      END


      SUBROUTINE DIAGIJt
!@sum  DIAGIJt produces lat-lon fields as maplets (6/page) or full-page
!@+    digital maps, and binary (netcdf etc) files (if qdiag=true)
!@auth Jean Lerner (adapted from work of G. Russell,M. Kelley,R. Ruedy)
!@ver   1.0
      USE DAGCOM                         !kdiag
      USE CONSTANT, only: teeny
      USE MODEL_COM, only:
     &     im,jm,lm,
     &     JHOUR,JHOUR0,JDATE,JDATE0,AMON,AMON0,JYEAR,JYEAR0,
     &     NDAY,Itime,Itime0,XLABEL,LRUNID
      USE TRACER_COM
      USE TRACER_DIAG_COM

      IMPLICIT NONE

      integer, parameter :: ktmax = ktaijln+ktaijkn+ktaijs
!@var Qk: if Qk(k)=.true. field k still has to be processed
      logical, dimension (ktmax) :: Qk
!@var Iord: Index array, fields are processed in order Iord(k), k=1,2,..
!@+     only important for fields 1->nmaplets+nmaps (appear in printout)
!@+     Iord(k)=0 indicates that a blank space replaces a maplet
      INTEGER Iord(ktmax+10),nmaplets,nmaps ! 10 extra blank maplets
      INTEGER kmaplets,nt(ktmax+10),ijtype(ktmax+10)
      REAL*8, DIMENSION(IM,JM) :: SMAP
      REAL*8, DIMENSION(JM) :: SMAPJ
      CHARACTER xlb*32,title*48,lname*80,name*30,units*30
!@var LINE virtual half page (with room for overstrikes)
      CHARACTER*133 LINE(53)
      INTEGER ::  I,J,K,kx,L,M,N,kcolmn,nlines,jgrid,irange
      DOUBLE PRECISION :: DAYS,gm

      if (kdiag(8).ge.1) return
C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      IF(QDIAG)call open_ij(trim(acc_period)//'.ijt'//XLABEL(1:LRUNID) 
     *     ,im,jm)

C**** INITIALIZE CERTAIN QUANTITIES
C**** standard printout
!     nmaplets = ktmax   ! ktaijln+ktaijkn+ktaijs
      nmaps = 0

C**** Fill in maplet indices for tracer concentrations
      k = 0
      do n=1,ntm
      if (itime.lt.itime_tr0(n)) cycle
      do l=1,lm
        k = k+1
        iord(k) = l
        nt(k) = n
        ijtype(k) = 1
      end do
C**** Fill in maplet indices for tracer sums and means
      do l=1,ktaijk
        k = k+1
        iord(k) = l
        nt(k) = n
        ijtype(k) = 2
      end do
      end do
C**** Fill in maplet indices for sources and sinks
      do kx=1,ktaijs
      n = ijts_source(kx)
      if (itime.lt.itime_tr0(n)) cycle
        k = k+1
        iord(k) = kx
        nt(k) = kx
        ijtype(k) = 3
      end do

      nmaplets = k
c**** always skip unused fields
      Qk = .true.
      do k=1,ktmax
        if (index(lname_ij(k),'unused').gt.0) Qk(k) = .false.
      end do

      xlb=acc_period(1:3)//' '//acc_period(4:12)//' '//XLABEL(1:LRUNID)
C****
      DAYS=(Itime-Itime0)/DFLOAT(nday)
C**** Collect the appropriate weight-arrays in WT_IJ
        wt_ij(:,:,:) = 1.

C**** Fill in the undefined pole box duplicates
      do i=2,im
        taijln(i,1,:,:) = taijln(1,1,:,:)
        taijln(i,jm,:,:) = taijln(1,jm,:,:)
        taijkn(i,1,:,:) = taijkn(1,1,:,:)
        taijkn(i,jm,:,:) = taijkn(1,jm,:,:)
        taijs (i,1,:) = taijs(1,1,:)
        taijs (i,jm,:) = taijs(1,jm,:)
      end do

C**** Print out 6-map pages
      do n=1,nmaplets
        if (mod(n-1,6) .eq. 0) then
c**** print header lines
          write (6,'("1",a)') xlabel
          write (6,902) jyear0,amon0,jdate0,jhour0,
     *      jyear,amon,jdate,jhour,itime,days
        end if
        kcolmn = 1 + mod(n-1,3)
        if (kcolmn .eq. 1) line=' '
c**** Find, then display the appropriate array
        if (Iord(n) .gt. 0 .and. Qk(n)) then
          call ijt_mapk (ijtype(n),nt(n),
     &         Iord(n),smap,smapj,gm,jgrid,irange,name,lname,units)
          title=trim(lname)//' ('//trim(units)//')'
          call maptxt(smap,smapj,gm,irange,title,line,kcolmn,nlines)
          if(qdiag) call pout_ij(title//xlb,name,lname,units,
     *                            smap,smapj,gm,jgrid)
          Qk(n) = .false.
        end if
c**** copy virtual half-page to paper if appropriate
        if (kcolmn.eq.3 .or. n.eq.nmaplets) then
          do k=1,nlines
            write (6,'(a133)') line(k)
          end do
        end if
      end do

      if (.not.qdiag) RETURN
C**** produce binary files of remaining fields if appropriate
      do k=1,ktmax
        if (Qk(k)) then
          call ijt_mapk (ijtype(n),nt(n),
     &         Iord(n),smap,smapj,gm,jgrid,irange,name,lname,units)
          title=trim(lname)//' ('//trim(units)//')'
          call pout_ij(title//xlb,name,lname,units,smap,smapj,gm,jgrid)
        end if
      end do
      if(qdiag) call close_ij

      RETURN
C****
  902 FORMAT ('0',15X,'From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,'  Model-Time:',I9,5X,
     *  'Dif:',F7.2,' Days')
      END SUBROUTINE DIAGIJt


      subroutine IJt_MAPk(nmap,n,
     &      k,smap,smapj,gm,jgrid,irange,name,lname,units)
!@sum ijt_MAPk returns the map data and related terms for the k-th field
!@+   for tracers and tracer sources/sinks
      USE DAGCOM
      USE CONSTANT, only: teeny
      USE GEOM, only: dxyp
      USE MODEL_COM, only:im,jm, IDACC
      use TRACER_DIAG_COM

      IMPLICIT NONE

      REAL*8, DIMENSION(IM,JM) :: anum,adenom,smap
      REAL*8, DIMENSION(JM) :: smapj
      integer i,j,k,iwt,jgrid,irange,n,nmap
      character(len=30) name,units
      character(len=80) lname
      real*8 :: gm,nh,sh, off, byiacc
!@var isumz,isumg = 1 or 2 if zon,glob sums or means are appropriate
      integer isumz,isumg

      isumz = 2 ; isumg = 2  !  default: in most cases MEANS are needed
      iwt = 1 ; jgrid = 1

      adenom = 1.
c**** tracer amounts
      if (nmap.eq.1) then
        name = sname_ijt(k,n) 
        lname = lname_ijt(k,n) 
        units = units_ijt(k,n)
        irange = ir_ijt(n)
        byiacc = 1./(idacc(ia_ijt)+teeny)
        do j=1,jm
        do i=1,im
          anum(i,j)=taijln(i,j,k,n)*byiacc*scale_ijt(k,n)/dxyp(j)
        end do
        end do
c**** tracer sums and means
      else if (nmap.eq.2) then
        name = sname_ijkn(k,n) 
        lname = lname_ijkn(k,n) 
        units = units_ijkn(k,n)
        irange = ir_ijt(n)
        byiacc = 1./(idacc(ia_ijt)+teeny)
        do j=1,jm
        do i=1,im
          anum(i,j)=taijkn(i,j,k,n)*byiacc*scale_ijkn(k,n)/dxyp(j)
        end do
        end do
c**** Sources and sinks
      else if (nmap.eq.3) then
        name = sname_ijts(k) 
        lname = lname_ijts(k) 
        units = units_ijts(k)
        irange = ir_ijt(1)
        byiacc = 1./(idacc(ia_ijts(k))+teeny)  
        do j=1,jm
        do i=1,im
          anum(i,j)=taijs(i,j,k)*byiacc*scale_ijts(k)/dxyp(j)
        end do
        end do
C**** PROBLEM
      else  ! should not happen
        write (6,*) 'no field defined for ijt_index',k
        stop 'ijt_mapk: undefined extra ij_field for tracers'
      end if

c**** Find final field and zonal, global, and hemispheric means
  100 continue
      call ij_avg (anum,adenom,wt_ij(1,1,iwt),jgrid,isumz,isumg, ! in
     *             smap,smapj,gm,nh,sh)                    ! out
      return
      end subroutine ijt_mapk


      SUBROUTINE io_trdiag(kunit,it,iaction,ioerr)
!@sum  io_trdiag reads and writes tracer diagnostics arrays to file
!@auth Jean Lerner
!@ver  1.0
      USE MODEL_COM, only: ioread,iowrite,iowrite_mon,iowrite_single
     *     ,irsfic,irerun,ioread_single,lhead
      USE TRACER_DIAG_COM, only: TACC,ktacc
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "TRdiag01"
!@var it input/ouput value of hour
      INTEGER, INTENT(INOUT) :: it
!@var TACC4(KTACC) dummy array for reading diagnostics files
      REAL*4 TACC4(KTACC)

      write (MODULE_HEADER(lhead+1:80),'(a,i8,a)')
     *   'R8 TACC(',ktacc,'),it'

      SELECT CASE (IACTION)
      CASE (IOWRITE,IOWRITE_MON) ! output to standard restart file
        WRITE (kunit,err=10) MODULE_HEADER, TACC,it
      CASE (IOWRITE_SINGLE)    ! output to acc file
        MODULE_HEADER(LHEAD+1:LHEAD+2) = 'R4'
        WRITE (kunit,err=10) MODULE_HEADER, SNGL(TACC),it
      CASE (IOREAD:)          ! input from restart file
        SELECT CASE (IACTION)
        CASE (IRSFIC)           ! initial conditions
        CASE (ioread_single)    ! accumulate diagnostic files
          READ (kunit,err=10) HEADER,TACC4
          TACC = TACC+TACC4     ! accumulate diagnostics
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version",HEADER
     *           ,MODULE_HEADER
            GO TO 10
          END IF
        CASE (ioread,irerun)    ! restarts
          READ (kunit,err=10) HEADER, TACC,it
          IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version",HEADER
     *           ,MODULE_HEADER
            GO TO 10
          end IF
        END SELECT
      END SELECT

      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_trdiag



