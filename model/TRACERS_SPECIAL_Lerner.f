#include "rundeck_opts.h"
#ifdef TRACERS_ON
#ifdef TRACERS_SPECIAL_Lerner

      MODULE PRATHER_CHEM_COM
!@sum Variables for chemical tracer routines that were provided by
!@+    Michael Prather.  The chemistry is parameterized as frequencies
C**** Stratospheric chemistry: consolidated commons
C**** These variables are used by both ozone and strat chem routines
c     EQUIVALENCE (SUT,TXM),(SVT,TYM),(SWT,TZM),(SUU,TXXM),(SVV,TYYM),
c    *  (SWW,TZZM),(SUV,TXYM),(SUW,TZXM),(SVW,TYZM),  (STT,TRM)
      USE MODEL_COM, only: jm,lm,ls1
!@var NSTRTC=# of strat chem layers (counting top down)
      integer, parameter :: nstrtc=lm-ls1
      integer jlatmd(jm)
      real*8 p0l(lm+1)

      contains
      subroutine set_prather_constants
      USE MODEL_COM, only: jm,lm,psfmpt,sige,ptop
      implicit none
      real*8 yedge(jm+1),yedge1,yedgen,xlatmd
      integer j,jxxx,l,lr

C---calculate nearest latitude to std lat's
      yedge1=-90.
      yedgen= 90.
      do j=1,jm+1
        yedge(j)=yedge1*float(jm+1-j)/float(jm) +
     *           yedgen*float(j-1)/float(jm)
      end do
      do j=1,jm
        xlatmd = 0.5*(yedge(j)+yedge(j+1))
        jxxx = xlatmd/10. + 10.
        jlatmd(j) = min(18,max(1,jxxx))
      end do
C---just calculate average P(mbar) at edge of each level
      do lr=1,lm+1
        p0l(lr) = sige(lm+2-lr)*psfmpt+ptop
      end do
      return

      end subroutine set_prather_constants
      end MODULE PRATHER_CHEM_COM


      MODULE TRACER_MPchem_COM
!@sum Variables for Prather's Stratospheric chemistry loss model
      USE MODEL_COM, only: jm,lm
      USE TRACER_COM, only: ntm
!@var n_MPtable_max:  Number of tracers that will use the frequency 
!@+    tables and share strat chem code
      integer, parameter :: n_MPtable_max=3
!@var n_MPtable: Index for tracers that use the frequency tables
      integer, dimension(ntm) :: n_MPtable
!@var tcscale: Scale factor for frequency tables
      real*8, dimension(n_MPtable_max) :: tcscale
!@param lz_tscparm Number of heights in stratchem tables
      integer, parameter :: lz_tscparm=20
!@var tscparm: Contains mean loss prequency in grid box
      real*4 tscparm(lz_tscparm,18,12,n_MPtable_max)
!@var TLtrm,TLtzm,TLtzzm: loss freq and moments of loss freq from tables
      real*8, dimension (jm,lm,n_MPtable_max) :: tltrm,tltzm,tltzzm
      end MODULE TRACER_MPchem_COM


      SUBROUTINE STRATCHEM_SETUP(nsc,tname)
C**** Prather stratospheric chemistry
      USE FILEMANAGER, only: openunit,closeunit
      USE TRACER_MPchem_COM, only: tscparm,lz_tscparm
      USE PRATHER_CHEM_COM, only: set_prather_constants
      implicit none
!     nsc = n_MPtable(n)
      integer nsc,j,k,m,iu
      character*80 titlch
      character*8 tname
      character*16 filein

      if (nsc.eq.0) stop ' NSC=0 in STRATCHEM_SETUP'
      filein = trim(tname)//'_TABLE'
      call openunit(filein,iu,.false.,.true.)
      read (iu,'(a)')   titlch
      write(6,'(1x,a)') titlch
      do m=1,12
        do j=1,18
          read(iu,'(20x,6e10.3/(8e10.3))')
     *          (tscparm(k,j,m,nsc),k=lz_tscparm,1,-1)
        end do
      end do
      call closeunit(iu)
      write(6,'(2A)') ' STRATCHEM TABLES READ for ',tname

      call set_prather_constants

      return
      end SUBROUTINE STRATCHEM_SETUP


      SUBROUTINE Strat_chem_Prather(n,ns)
!@sum Strat_chem_Prather calculates stratospheric chemistry for 
!@+        N2O, CFC and CH4
!@auth Michael Prather (J.Lerner adopted code)
!@var nsc index for chemical tracer number (to access tables)
!@var T0L is the amount 'lost' to chemistry this time step
!@var facbb: APPLY AN AD-HOC FACTOR TO BRING CH4 INTO BALANCE
      USE TRACER_COM
      USE TRACER_DIAG_COM, only : tajls,jls_3Dsource,itcon_3Dsrc
      USE MODEL_COM, only: im,jm,lm,ls1,dtsrc
      USE GEOM, only: imaxj
      USE DYNAMICS, only: ltropo
      USE FLUXES, only: tr3Dsource
      USE QUSDEF, only : mz,mzz
      USE TRACER_MPchem_COM, only: tltrm,tltzm,tltzzm,n_MPtable,tcscale
      implicit none
      integer i,j,l,lr,n,ns,najl,nsc
      real*8, parameter :: by3=1./3.d0,by7=1./7.d0
      real*8 facbb,told(im,jm,lm)
      real*8 f0l,f1l,f2l,t0l,f1l,f2l,g0l,g1l,g2l,t0l,t1l,t2l

      nsc = n_MPtable(n)
      facbb = 1.
      if (trname(n).eq.'CH4') facbb = (40.d0/25.73d0)*(40.d0/35.177d0)
C-----STRATOSPHERIC LOSS occurs in top NSTRTC layers (ls1+1 to lm)
C-----ALLOW FOR ONLY ONE SET OF TSCPARM's FOR STRATOSPHERIC LOSS,
C-----uses TCSCALE for different tracers to scale loss
C-----uses S.O.M. formulation for vertical losses
C-----NOTE that TLTRM(J,LR,N) stored from top (=LM) down
      
      told(:,:,:) = trm(:,:,:,n)
      do 150 l=lm,ls1+1,-1
      lr = lm+1-l
      do 140 j=1,jm
  !     IF (L.LT.LMSTRT(J)) GO TO 140
C-----TSCPARM->TLtrm contains mean loss freq in grid box:
        f0l = max(TLtrm(j,lr,nsc),0d0)
        if (f0l.le.0.) go to 140
        f1l = tltzm(j,lr,nsc)  ! FOM of loss freq from tables
        f2l = tltzzm(j,lr,nsc)  ! SOM of loss freq from tables
        do 130 i=1,imaxj(j)
!!?????   if (l.le.ltropo(i,j)) go to 130
          if (trm(i,j,l,n).le.0.) go to 130
C------Couple the moments of the loss freq with moments of the tracer:
C       dSo  = dt*( So*Lo + Sz*Lz/3 + Szz*Lzz/5)
C       dSz  = dt*( Sz*Lo + So*Lz + 2(Sz*Lzz + Szz*Lz)/5)
C       dSzz = dt*( Szz*Lo + SoLzz +2Sz*Lz/3 + 2Szz*Lzz/7)
          t0l = f0l*trm(i,j,l,n) + 
     *           f1l*trmom(mz,i,j,l,n)*by3+f2l*trmom(mzz,i,j,l,n)*.2d0
          if (t0l.lt.0.) go to 130
          t1l = f0l*trmom( mz,i,j,l,n) + f1l*trm(i,j,l,n) +
     *      0.4*(f2l*trmom(mz,i,j,l,n)    +f1l*trmom(mzz,i,j,l,n))
          t2l = f0l*trmom(mzz,i,j,l,n) + f2l*trm(i,j,l,n) +
     *      2.0*(f1l*trmom(mz,i,j,l,n)*by3+f2l*trmom(mzz,i,j,l,n)*by7)
C---calculate e-folding of tracer mass (trm) & scale moment losses
C---  to this change (T0L):
          g0l = t0l/trm(i,j,l,n)
          g1l = t1l/t0l
          g2l = t2l/t0l
          t0l = (1.0 - exp(-g0l*dtsrc*tcscale(nsc)))*trm(i,j,l,n)
          t0l = t0l*facbb  ! APPLY AN AD-HOC FACTOR
          trm(i,j,l,n) = trm(i,j,l,n) - t0l
          trmom( mz,i,j,l,n) = trmom( mz,i,j,l,n) - t0l*g1l
          trmom(mzz,i,j,l,n) = trmom(mzz,i,j,l,n) - t0l*g2l
  130     CONTINUE
  140     CONTINUE
  150   CONTINUE
      najl = jls_3Dsource(ns,n)
      do l=1,lm
      do j=1,jm
      do i=1,imaxj(j)
        tajls(j,l,najl) = tajls(j,l,najl)+(trm(i,j,l,n)-told(i,j,l))
      end do
      end do
      end do
      call DIAGTCA(itcon_3Dsrc(ns,n),n)

      return
      END SUBROUTINE Strat_chem_Prather


      SUBROUTINE STRTL
C**** This is called at the beginning of each month
C**** Prather strat chem
      USE MODEL_COM, only: jm,lm,jmon
      USE PRATHER_CHEM_COM, only: nstrtc,jlatmd,p0l
      USE TRACER_MPchem_COM, only: tscparm,n_MPtable_max,
     *    tltrm,tltzm,tltzzm,lz_tscparm
      implicit none
C-----------------------------------------------------------------------
C---monthly set up of chemical loss parameters
      real*8 strt0l(lm),strt1l(lm),strt2l(lm),strtx(lz_tscparm)
      real*8 ps(lz_tscparm+6),f(lz_tscparm+5)
C-----------------------------------------------------------------------
C--tscparm(lz_tscparm,18,12,N) defined for 18 lats (85S, 75S, ...85N) 
C --                   & 12 months
C----  do NOT interpolate, just pick nearest latitude
C---assume given MONTH = month #, NTM=# tracers, JM=#lats, etc.
      integer n,j,jj,k,lr

      DO 800 N=1,n_MPtable_max
        DO 700 J=1,JM
          JJ = JLATMD(J)
          DO 500 K=1,lz_tscparm
            STRTX(K) = tscparm(K,JJ,jmon,N)
500       CONTINUE
          CALL STRT2M(STRTX,lz_tscparm,STRT0L,STRT1L,STRT2L,P0L,NSTRTC
     *      ,ps,f)
C----store loss freq & moments in TLtrm/tzm/tzzm for exact
C---- CTM layers LM down
          DO 600 LR=1,NSTRTC
            TLtrm(J,LR,N) = STRT0L(LR)
            TLtzm(J,LR,N) = STRT1L(LR)
            TLtzzm(J,LR,N) = STRT2L(LR)
600       CONTINUE
700     CONTINUE
800   CONTINUE
      RETURN
      END SUBROUTINE STRTL

      SUBROUTINE Trop_chem_CH4(n,ns)
!@sum Trop_chem_CH4 calculates tropospheric chemistry for CH4
!@+     by applying a pre-determined chemical loss rate
!@auth Jean Lerner
      USE TRACER_COM
      USE DYNAMICS, only: ltropo
      USE MODEL_COM, only: im,jm,lm,ls1,byim,jyear,jhour,jday,itime
      USE GEOM, only: imaxj
      USE FLUXES, only: tr3Dsource
      USE FILEMANAGER, only: openunit,closeunit
      implicit none
      integer n,ns,i,j,l,ifile
      REAL*4 frqlos(im,jm,ls1),taux
      real*8 tauy,tune
      parameter (tune = 445./501.)
      logical, save :: ifirst=.true.
      save frqlos,tauy

C**** Read chemical loss rate dataset (5-day frequency)
      if (mod(jday,5).gt.0 .and. .not.ifirst) go to 550
      if (jhour.ne.0 .and. .not.ifirst) go to 550
      ifirst = .false.
      call openunit('CH4_TROP_FRQ',ifile,.true.,.true.)
  510 continue
      read (ifile,end=515) taux,frqlos
      tauy = nint(taux)+(jyear-1950)*8760.
      IF (itime+ 60..gt.tauy+120.) go to 510
      IF (itime+180..le.tauy+120.) then
        write(6,*)' PROBLEM MATCHING itime ON FLUX FILE',TAUX,TAUY,JYEAR
        STOP 'PROBLEM MATCHING itime ON FLUX FILE in Trop_chem_CH4'
      end if
      rewind (ifile)
      go to 518
C**** FOR END OF YEAR, USE FIRST RECORD
  515 continue
      rewind (ifile)
      read (ifile,end=515) taux,frqlos
      rewind (ifile)
      tauy = nint(taux)+(jyear-1950)*8760.
  518 continue
      call closeunit(ifile)
      WRITE(6,'(2A,2F10.0,2I10)')
     * ' *** Chemical Loss Rates in Trop_chem_CH4 read for',
     * ' taux,tauy,itime,jyear=', taux,tauy,itime,jyear
C**** AVERAGE POLES
      do l=1,ls1
        frqlos(1, 1,l) = sum(frqlos(:, 1,l))*byim
        frqlos(1,jm,l) = sum(frqlos(:,jm,l))*byim
      end do
C**** APPLY AN AD-HOC FACTOR TO BRING INTO BALANCE
      frqlos(:,:,:) = frqlos(:,:,:)*tune
C**** Apply the chemistry
  550 continue
      do j=1,jm
      do i=1,imaxj(j)
        do l = 1,ltropo(i,j)
          if (l.le.ls1)
     *     tr3Dsource(i,j,l,ns,n) = -frqlos(i,j,l)*trm(i,j,l,n)
        end do
      end do
      end do
!     IF (jhour.eq.0.) write(6,'(a,i10,F10.0)')
!    *  ' Tropo chem for CH4 performed',itime,TAUY
      return
      END SUBROUTINE Trop_chem_CH4


      MODULE LINOZ_CHEM_COM
C**** linoz with sol variability 
!@sum Variables for linoz chemistry.  Original code was provided by
!@+    Michael Prather and Chris McLinden.  
C**** Ozone tracer; linoz chemistry
C     NCTABLE=# of linoz tables (includes solar UV)
C     n_O3=tracer number for linoz O3
!@dbparam dsol describes portion of solar cycle being modeled for linoz
!@+      +1.0 = solar max, 0.0 = neutral, -1.0 = solar min
      USE MODEL_COM, only: im,jm,lm,ptop,psfmpt,sige,dtsrc
      USE TRACER_COM, only: ntm,tr_mm
      USE PRATHER_CHEM_COM, only: set_prather_constants
      USE CONSTANT, only: mair
      implicit none
!@param lz_tlparm Number of heights in linoz tables
      integer, PARAMETER :: lz_tlparm=25,nctable=8
C****         25 heights, 18 lats, 12 months, 7 parameters
      real*8 TLPARM(lz_tlparm,18,12,nctable),
     *  TLT0M(JM,LM,nctable),TLTZM(JM,LM,nctable),TLTZZM(JM,LM,nctable)
      real*8 tmmvv(ntm),dsol
!@var lbc Top layer for ozone boundary conditions in troposphere
      integer lbc
      real*8 dtchem

      contains

      SUBROUTINE LINOZ_SETUP(n_O3)
C**** Needed for linoz chemistry
      USE FILEMANAGER, only: openunit,closeunit
      implicit none
      integer iu,j,k,l,m,n,n_O3
      character*80 titlch

      call openunit('LINOZ_TABLE',iu,.false.,.true.)
      read (iu,'(a)')   titlch
      write(6,'(1x,a)') titlch
      do n=1,nctable
        read (iu,'(a)')   titlch
        write(6,'(1x,a)') titlch
        do m=1,12
          do j=1,18
            read(iu,'(20x,6e10.3/(8e10.3))') 
     *           (tlparm(k,j,m,n),k=lz_tlparm,1,-1)
          end do
        end do
      end do
      write(6,'(a)') ' linoz tables read'
      call closeunit(iu)

C**** Calculate level for tropophere ozone chem
      do l=1,lm
        if (psfmpt*sige(l+1)+ptop .le. 900.)then
          lbc = l
          write(6,'(a,i3,f8.1)') 
     *     ' Top layer for tropo O3 chem is ',lbc, psfmpt*sige(l)+ptop
          exit
        end if
      end do
C**** Get other useful things
      call set_prather_constants
      tmmvv(n_O3) = mair/tr_mm(n_O3)
      dtchem = dtsrc
      return

      end SUBROUTINE LINOZ_SETUP
      end MODULE LINOZ_CHEM_COM


      SUBROUTINE Trop_chem_O3(n,ns)
c
c-----------------------------------------------------------------------
c   Lower boundary conditions for Linearized Strat. Chem.
c-----------------------------------------------------------------------
c
      USE DYNAMICS, only: am   ! Air mass of each box (kg/m^2)
      USE TRACER_COM
      USE MODEL_COM, only: jm,itime
      USE GEOM, only: imaxj,dxyp
      USE TRACER_DIAG_COM, only : tajls,jls_3Dsource,itcon_3Dsrc
      USE LINOZ_CHEM_COM, only: lbc,dtchem,tmmvv
      USE FLUXES, only: tr3Dsource
      implicit none
      real*8, parameter :: taubc=172800.d0  !in seconds
      real*8, parameter :: tmrbc=20.0d-9    !mixing ratio
      integer i,j,l,n,najl,ns
      real*8 coeff,ratio,T0Mold,dmass,sdmass,scalmom

c lower boundary condition : relax each species to tmrbc(n)
c using a lifetime of taubc(n) for lowest lbc levels

      coeff = 1.0-exp(-dtchem/taubc)
      ratio = tmrbc / tmmvv(n)
      najl = jls_3Dsource(ns,n)
        do l=1,lbc
        do j=1,jm
        sdmass = 0.
        do i=1,imaxj(j)
          T0Mold=trm(i,j,l,n)
          dmass = (ratio*am(l,i,j)*dxyp(j)-trm(i,j,l,n))*coeff
          if(trm(i,j,l,n) + dmass .lt.0.) then
               write(6,'(a,3i4,2f15.2,i9)') 
     *           ' Negative tracer in Trop_chem_O3',
     *                i,j,l,trm(i,j,l,n),dmass,itime
               dmass = -trm(i,j,l,n)
               trm(i,j,l,n) = 0.d0
            else
               trm(i,j,l,n) = trm(i,j,l,n) + dmass
          end if
!         tr3Dsource(i,j,l,ns,n) = (T0Mold + dmass)/dtchem
          sdmass = sdmass+dmass
c scale remaining moments by fractional change in total tracer mass
        if (dmass.lt.0.d0) then
          scalmom = trm(i,j,l,n)/T0Mold
          trmom(1:nmom,I,J,L,n) = trmom(1:nmom,I,J,L,n) * scalmom
        end if
        enddo
          tajls(j,l,najl) = tajls(j,l,najl) + sdmass
        enddo
        enddo
      call DIAGTCA(itcon_3Dsrc(ns,n),n)
      RETURN
      END SUBROUTINE Trop_chem_O3


      SUBROUTINE Strat_chem_O3(n,ns)
c
c-----------------------------------------------------------------------
c   CHEM3 applies linearized chemistry based on tables from
c    PRATMO model using climatological T, O3, time of year
c-----------------------------------------------------------------------
c  stratospheric chem occurs in top NSTRTC layers of CTM
c  TLT0M(J,LR,N) is stored LR from top (=LM) down (=LM+1-NCSTRT)
c
c Stratospheric Chemistry Tables for O3:
c ======================================
c   7 tables, each a function of month (12), latitude
c   (18, -85 to 85 in 10 deg. increments) and altitude
c   (25, z*=10-58 km in 2 km increments).
c  1- ozone (Logan climatology), v/v
c  2- Temperature climatology, K
c  3- Column ozone climatology, Logan ozone integrated above box, DU
c  4- ozone (P-L) for climatological ozone, v/v/s
c  5- d(P-L) / dO3, 1/s
c  6- d(P-L) / dT, v/v/s/K
c  7- d(P-L) / d(column O3), v/v/s/DU
c  8- d(P-L) / d(sol.flx.) (-) 
c
c dsol is variable which describes portion of solar cycle being modeled
c +1.0 = solar max, 0.0 = neutral, -1.0 = solar min
c can go beyond +1,-1 for more extreme situations
c (but remember this is a linear approximation)
c use dsol=0.0 for 'standard linoz' runs

      USE MODEL_COM, only: itime,im,jm,lm,t
      USE DYNAMICS, only: am,ltropo   ! Air mass of each box (kg/m^2)
      USE GEOM, only: imaxj,dxyp
      USE TRACER_COM
      USE DYNAMICS, only: pk
      USE TRACER_DIAG_COM, only : tajls,jls_3Dsource,itcon_3Dsrc
      USE PRATHER_CHEM_COM, only: nstrtc
      USE LINOZ_CHEM_COM, only: dtchem,tmmvv,tlT0M,TLTZM,TLTZZM,dsol
      USE FLUXES, only: tr3Dsource
      implicit none
      real*8 dero3,scalmom,dcolo3(im,jm,lm),
     &  pmltot,dertmp,dtmp,derco3,colo3(im,jm,lm),dco3,sso3,
     &  climo3,climpml,dersol
      real*8 dmass,T0Mold
      integer i,j,l,lr,lsmin,n,ns,najl

      lsmin   = lm+1 - nstrtc
      najl = jls_3Dsource(ns,n)
c start at top layer and continue to lowest layer for strat. chem
      DO 330 L = LM,LSMIN,-1
        LR = LM+1-L
        DO 320 J = 1,JM
!         IF (L .LT. LMSTRT(J))  GOTO 320
          do 310 i=1,imaxj(j)
            if (l.le.ltropo(i,j)) go to 310  ! stay above tropopause
            if (trm(i,j,l,n).le.0.d0) goto 310

c calculate ozone column above box (and save)
c   dcolo3 = ozone column (in DU) in given layer
c   colo3 =  ozone column above layer + half of column in layer
            if (l.eq.lm) then           !top model layer
              dcolo3(i,j,l) = trm(i,j,l,n) / dxyp(j) *
     &          6.022d23/(tr_mm(n)*1d-3)/ 2.687d16 * 1d-4
              colo3(i,j,l) = dcolo3(i,j,l)*0.5
            else
              dcolo3(i,j,l) = trm(i,j,l,n)/ dxyp(J) *
     &          6.022d23/(tr_mm(n)*1d-3)/ 2.687d16 * 1d-4
              colo3(i,j,l) = colo3(i,j,l+1) +
     &          (dcolo3(i,j,l)+dcolo3(i,j,l+1))*0.5
            endif

c ****** O3 Chemistry  ******
c store tracer mass before chemistry
            T0Mold=trm(i,j,l,n)
c climatological P-L:
            climpml=tlT0M(j,lr,4)/tmmvv(n)*am(l,i,j)*dxyp(j)
c local ozone feedback:
            dero3=tlT0M(j,lr,5)
            climo3=tlT0M(j,lr,1)/tmmvv(n)*am(l,i,j)*dxyp(j)
c column ozone feedback:
            derco3=tlT0M(j,lr,7)/tmmvv(n)*am(l,i,j)*dxyp(j)
            dco3=(colo3(i,j,l)-tlT0M(j,lr,3))
c temperature feedback: T is potential temp, need to convert
            dertmp=tlT0M(j,lr,6)/tmmvv(n)*am(l,i,j)*dxyp(j)
            dtmp=(t(i,j,l)*PK(I,J,L)-tlT0M(j,lr,2))
c define sol.flux. derivative and convert from mixing ratio to mass
            dersol = tlT0M(j,lr,8)/tmmvv(n)*am(l,i,j)*dxyp(j)
c calulate steady-state ozone:
c           sso3=climo3 - (climpml+dco3*derco3+dtmp*dertmp)/dero3
            sso3=climo3 -
     *        (climpml+dco3*derco3+dtmp*dertmp+dsol*dersol)/dero3
c change in ozone mass due to chemistry:
            dmass=(sso3-T0Mold)*(1.0-exp(dero3*dtchem))
c update ozone mass
            if (T0Mold+dmass.lt.0.) then
               write(6,'(a,3i4,2f15.2,i9)') 
     *           ' Negative tracer in Strat_chem_O3',
     *                i,j,l,T0Mold,dmass,itime
               dmass = -T0Mold
               trm(i,j,l,n) = 0.d0
            else
               trm(i,j,l,n) = T0Mold + dmass
            end if
!           tr3Dsource(i,j,l,ns,n) = (T0Mold + dmass)/dtchem
            tajls(j,l,najl) = tajls(j,l,najl) + dmass

c scale remaining moments by fractional change in total tracer mass
        if (dmass.lt.0.d0) then
          scalmom = trm(i,j,l,n)/T0Mold
          trmom(1:nmom,I,J,L,n) = trmom(1:nmom,I,J,L,n) * scalmom
        end if
  310 continue
  320 continue
  330 continue
      call DIAGTCA(itcon_3Dsrc(ns,n),n)
      return
      end SUBROUTINE Strat_chem_O3


c-----------------------------------------------------------------------
      SUBROUTINE linoz_STRATL
c-----------------------------------------------------------------------
c-------- monthly fixup of chemistry PARAM'S
c
      USE PRATHER_CHEM_COM, only: jlatmd,p0l,NSTRTC
      USE LINOZ_CHEM_COM, only: nctable,lz_tlparm,TLPARM,
     *    tlT0M,TLTZM,TLTZZM
      USE MODEL_COM, only: jmon,jm,lm
      implicit none
      real*8  STRT0L(LM),STRT1L(LM),STRT2L(LM),STRTX(lz_tlparm)
      real*8 ps(lz_tlparm+6),f(lz_tlparm+5)
      integer j,jj,k,lr,n

c-------- TLPARM(25,18,12,N) defined for -----------------------------
c lz_tlparm  25 layers from 58 km to 10 km by 2 km intervals
c            18 LATS (85S, 75S, ...85N)
c            12 months
c            N tables = NCTABLE
c-------- skip interpolating, pick nearest latitude --------------------

      DO N = 1,NCTABLE
      DO J = 1,JM
        JJ = JLATMD(J)
        DO K = 1,lz_tlparm
c nearest latitude, single month
          STRTX(K) = TLPARM(K,JJ,jmon,N)
        ENDDO

c-------- stratospheric chem occurs in top NSTRTC layers ---------------
c
        CALL STRT2M(STRTX,lz_tlparm,STRT0L,STRT1L,STRT2L,P0L,NSTRTC
     *      ,ps,f)
c
c-------store loss freq/yields & moments in TLT0M/TLTZM/TLTZZM
c-------             for exact CTM layers LM down
          DO LR = 1,NSTRTC
            TLT0M(J,LR,N) = STRT0L(LR)
            TLTZM(J,LR,N) = STRT1L(LR)
            TLTZZM(J,LR,N) = STRT2L(LR)
          ENDDO
        ENDDO   ! J
      ENDDO     ! N
      return
      END SUBROUTINE linoz_STRATL


      SUBROUTINE STRT2M (STRTX,NX,STRT0L,STRT1L,STRT2L,P0L,NSTRT,ps,f)
c-----------------------------------------------------------------------
      implicit none
      integer  NX,NSTRT,nl,l,k
      real*8   P0L(*),STRT0L(*),STRT1L(*),STRT2L(*),STRTX(*)
      real*8    P1,P2,F0,F1,F2,PS(*),F(*)
      real*8    XPSD,XPSLM1,XPSL
c-----------------------------------------------------------------------
c       set up std z* atmosphere: p = 1000 * 10**(-z*/16 km)
c       assume that stratospheric chemical parameters always start at
cc       52 km (N=27) scan downward from 52 km to 14 km (NX=20) by 2 km
c       58 km (N=30) scan downward from 58 km to 10 km (NX=25) by 2 km
c       intervals, constant >58km
c-------- N.B. F(@30km) assumed to be constant from 29-31 km (by mass)
c
      NL = NX+5  !??? in old start code nl=nx+7.....
      XPSD       = 10.D0 **(-0.125D0)      !Z=2 km
      XPSLM1     = 1000.D0                 !Z=0 km, P=1000 mb
      PS(1)      = 1000.D0
      DO L = 2,NL
        XPSL     = XPSLM1 *XPSD
        PS(L)    = 0.5D0 *(XPSLM1 +XPSL) !Z=2*I km, P=1000*10**(-Z/16)
        XPSLM1   = XPSL
      ENDDO
      PS(NL+1)   = 0.D0
      DO L = 1,NL-NX
        F(L)     = 0.D0
      ENDDO
c-------- K=1 is at the top of atmosphere ------------------------------
      DO K = 1,NX
        F(NL+1-K)= STRTX(K)
      ENDDO
      DO K = 1,NSTRT
        P1       = P0L(K+1)
        P2       = P0L(K)
        CALL SOMLFQ(P1,P2,F0,F1,F2,PS,F,NL)
        STRT0L(K)= F0
        STRT1L(K)= F1
        STRT2L(K)= F2
      ENDDO
      RETURN
      END SUBROUTINE STRT2M


      SUBROUTINE SOMLFQ(P1,P2,F0,F1,F2,PS,F,NL)
c-----------------------------------------------------------------------
c--- calculate loss freq moments from a set of loss freq's at std z*
c--  given a CTM model interval pressure range: P1 > P2 (decreasing up)
c-----  the pressure levels BETWEEN z* values are:
c                      PS(i) > PS(i+1) bounds z*(i)
c----- NL:  z* levels, ==> PS(NL+1) = 0
c-----            (extrapolate chemical loss to top)
c      Z1 = 16.D0*LOG10(1000.D0/P1)
c      Z2 = 16.D0*LOG10(1000.D0/P2)
c
c----- The MOMENTS for a square-wave or 'bar':
c-----             F(x)=f0  b<=x<=c, =0.0 else
c---    S0 =   f0 (x)                      [from x=b to x=c]
c---    S1 = 3 f0 (x^2 - x)                [from x=b to x=c]
c---    S2 = 5 f0 (2x^3 - 3x^2 + x)        [from x=b to x=c]
c-----------------------------------------------------------------------
      implicit none
      integer  NL
      real*8   P1,P2,F0,F1,F2,PS(NL+1),F(NL),sgnf0
      integer  I
      real*8   XB,XC,PC,PB,THIRD
      parameter (THIRD = 1.D0/3.D0)
c-----------------------------------------------------------------------
      F0     = 0.D0
      F1     = 0.D0
      F2     = 0.D0
      DO I = 1,NL
        PC   = MIN(P1,PS(I))
        PB   = MAX(P2,PS(I+1))
        IF (PC .GT. PB)  THEN
C------ have condition:  P1>=PC > PB>=P2, 0<=XB < XC<=1 --------------
          XC = (PC-P2)/(P1-P2)
          XB = (PB-P2)/(P1-P2)
c
c------ assume that the loss freq, F, is constant over interval [XB,XC],
c------ F0: (c-b),  F1: 6((c2-c)-(b2-b)), F2: 5((2c3-3c2+c)-(2b3-3b2+b))
c------ calculate its contribution to the moments in the interval [0,1]
c
          F0 = F0 +F(I) *(XC-XB)
          F1 = F1 +F(I) *3.D0 *((XC*XC-XC) - (XB*XB-XB))
          F2 = F2 +F(I) *5.D0 *
     &         ((XC+XC-1.D0)*(XC*XC-XC) - (XB+XB-1.D0)*(XB*XB-XB))
        ENDIF
      ENDDO
c
c-------- RESTRAIN moments: force monotonicity & positive at min end pt
c
c -=-=- cam: tables can be + or -
      if (f0.ne.0.0) then
        sgnf0=f0 / abs(f0)
      else
        sgnf0=1.0
      endif
      f0=abs(f0)

c      F0 = MAX(F0, 0.D0)
      IF (F2 .GT. 0.D0)  THEN
c
c-------- do not allow reversal of curvature: F2 > 0 -------------------
        F2   = MIN(F2, ABS(F1)*THIRD, 5.D-1*F0)
        IF (F1 .LT .0.D0)  THEN
          F1 = MAX(-(F0+F2), F1)
         ELSE
          F1 = MIN(+(F0+F2), F1)
        ENDIF
      ELSE
c
c-------- F2 < 0 = curved down at ends, allow if F1 < F0 ---------------
        F1  = MIN(F0,MAX(-F0,F1))
        F2  = MAX(F2,(ABS(F1)-F0),(-ABS(F1)*THIRD))
      ENDIF
c
c -=-=- cae: apply sign
      f0=sgnf0 * f0
      f1=sgnf0 * f1
      f2=sgnf0 * f2
      RETURN
      END SUBROUTINE SOMLFQ


      MODULE CH4_SOURCES
      USE TRACER_COM
!@var CH4_src CH4 surface sources and sinks (kg/s)
      integer, parameter :: nch4src=14
      real*8 CH4_src(im,jm,nch4src)
      END MODULE CH4_SOURCES


      subroutine read_CH4_sources(nt,iact)
!@sum reads in CH4 sources and sinks
!@auth Jean Lerner
C****
C**** There are 3 monthly sources and 11 annual sources
C**** Annual sources are read in at start and re-start of run only
C**** Monthly sources are interpolated each day
      USE MODEL_COM, only: itime,JDperY,im,jm,jday,focean,fearth
      USE LAKES_COM, only: flake
      USE CONSTANT, only: sday
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      USE TRACER_COM, only: itime_tr0,trname
      use CH4_SOURCES, only: src=>ch4_src,nsrc=>nch4src
      implicit none
      character*80 title
!@var adj Factors that tune the total amount of individual sources
      real*8 adj(nsrc)
      data adj/1.3847,1.0285,3.904,1.659,1.233,1.194,0.999,
!    *  7.2154, 3.7247d-5,3.1399d-4,5.4838d-5,   !Model II prime
     *  7.2154, 3.5997d-5,17.330d-4,5.3558d-5,
     *  0.4369,0.7533,0.9818/
!@var nanns,nmons: number of annual and monthly input files
      integer, parameter :: nanns=11,nmons=3
      integer ann_units(nanns-3),mon_units(nmons)
      character*12 :: ann_files(nanns-3) =
     *  (/'CH4_ANIMALS ','CH4_COALMINE','CH4_GASLEAK ','CH4_GASVENT ',
     *    'CH4_CITYDUMP','CH4_SOIL_ABS','CH4_TERMITES','CH4_COALBURN'/)
      logical :: ann_bins(nanns-3)=(/.true.,.true.,.true.,.true.,
     *        .true.,.true.,.true.,.true./)
      character*8 :: mon_files(nmons) =
     *   (/'CH4_BURN','CH4_RICE','CH4_WETL'/)
      real*8 adj_wet(jm)
      data adj_wet/15*0.6585,16*1.761,15*0.6585/ !zonal adj FOR WETLANDS
      data iwet/14/  !!! position of wetlands array in src
      logical :: mon_bins(nmons)=(/.true.,.true.,.true./)
      real*8 tlca(im,jm,nmons),tlcb(im,jm,nmons)  ! for monthly sources
      integer i,j,nt,iact,iu,k,jdlast,iwet
      logical :: ifirst=.true.
      integer :: jdlast=0
      save ifirst,jdlast,tlca,tlcb

      if (itime.lt.itime_tr0(nt)) return
C****
C**** Annual Sources and sinks
C**** Apply adjustment factors to bring sources into balance
C**** Annual sources are in KG C/M2/Y
C**** Sources need to be kg/m^2 s; convert /year to /s
C****
      if (ifirst) then
        k = 0
        call openunits(ann_files,ann_units,ann_bins,nanns-3)
        do iu = ann_units(1),ann_units(nanns-3)
          k = k+1
          call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          src(:,:,k) = src(:,:,k)*adj(k)/(sday*JDperY)
        end do
        call closeunits(ann_units,nanns-3)
        ! 3 miscellaneous sources
          k = k+1
          src(:,:,k) = focean(:,:)*adj(k)/(sday*JDperY)
          k = k+1
          src(:,:,k) =  flake(:,:)*adj(k)/(sday*JDperY)
          k = k+1
          src(:,:,k) = fearth(:,:)*adj(k)/(sday*JDperY)
      endif
C****
C**** Monthly sources are interpolated to the current day
C****
C**** Also, Apply adjustment factors to bring sources into balance
C**** Monthly sources are in KG C/M2/S => src in kg/m^2 s
C****
      write(6,*) ' In read_CH4_sources iact,ifirst',iact,ifirst
      if (iact.eq.1 .and. .not.ifirst) return
      ifirst = .false.
      call openunits(mon_files,mon_units,mon_bins,nmons)
      call sys_flush(6)
      j = 0
      do k = nanns+1,nsrc
        j = j+1
        write(6,*) j,mon_units(j)
        call sys_flush(6)
        call read_monthly_sources(mon_units(j),jdlast,
     *    tlca(1,1,j),tlcb(1,1,j),src(1,1,k))
        src(:,:,k) = src(:,:,k)*adj(k)
      end do
      jdlast = jday
      call closeunits(mon_units,nmons)
      write(6,*) trname(nt),'Sources interpolated to current day'
C****
C**** Zonal adjustment for combined wetlands and tundra
C****
      do j=1,jm
        src(:,j,iwet) = src(:,j,iwet)*adj_wet(j)
      end do
      return
      end subroutine read_CH4_sources


      MODULE CO2_SOURCES
      USE TRACER_COM
!@var co2_src C02 surface sources and sinks (kg/s)
      integer, parameter :: nco2src=6
      real*8 co2_src(im,jm,nco2src)
      END MODULE CO2_SOURCES


      subroutine read_CO2_sources(nt,iact)
!@sum reads in CO2 sources and sinks
!@auth Jean Lerner
C****
C**** There are two monthly sources and 4 annual sources
C**** Annual sources are read in at start and re-start of run only
C**** Monthly sources are interpolated each day
      USE MODEL_COM, only: itime,jday,JDperY,im,jm
      USE CONSTANT, only: sday
      USE FILEMANAGER, only: openunit,closeunit, openunits,closeunits
      USE TRACER_COM, only: itime_tr0,trname
      use CO2_SOURCES, only: src=>co2_src,nsrc=>nco2src
      implicit none
      character*80 title
!@var adj Factors that tune the total amount of individual sources
      real*8 adj(nsrc)
      data adj/3.81d0,3.67d0,3.67d0,19.54d0,3.67d0,6.42d0/
!@var nanns,nmons: number of annual and monthly input files
      integer, parameter :: nanns=4,nmons=2
      integer ann_units(nanns),mon_units(nmons)
      character*12 :: ann_files(nanns) =
     *  (/'CO2_FOS_FUEL','CO2_FERT    ','CO2_REGROWTH','CO2_LAND_USE'/)
      logical :: ann_bins(nanns)=(/.true.,.true.,.true.,.true./)
      character*9 :: mon_files(nmons) = (/'CO2_VEG  ','CO2_OCEAN'/)
      logical :: mon_bins(nmons)=(/.true.,.true./)
      real*8 tlca(im,jm,nmons),tlcb(im,jm,nmons)  ! for monthly sources
      integer i,j,nt,iact,iu,k,jdlast
      logical :: ifirst=.true.
      integer :: jdlast=0
      save ifirst,jdlast,tlca,tlcb

      if (itime.lt.itime_tr0(nt)) return
C****
C**** Annual Sources and sink
C**** Apply adjustment factors to bring sources into balance
C**** Annual sources are in KG C/M2/Y
C**** Sources need to be kg/m^2 s; convert /year to /s
C****
      if (ifirst) then
        call openunits(ann_files,ann_units,ann_bins,nanns)
        k = 0
        do iu = ann_units(1),ann_units(nanns)
          k = k+1
          call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          src(:,:,k) = src(:,:,k)*adj(k)/(sday*JDperY)
        end do
        call closeunits(ann_units,nanns)
      endif
C****
C**** Monthly sources are interpolated to the current day
C****
C**** Also, Apply adjustment factors to bring sources into balance
C**** Monthly sources are in KG C/M2/S => src in kg/m^2 s
      write(6,*) ' In read_CO2_sources iact,ifirst',iact,ifirst
      if (iact.eq.1 .and. .not.ifirst) return
      ifirst = .false.
      call openunits(mon_files,mon_units,mon_bins,nmons)
      j = 0
      do k=nanns+1,nsrc
        j = j+1
        call read_monthly_sources(mon_units(j),jdlast,
     *    tlca(1,1,j),tlcb(1,1,j),src(1,1,k))
        src(:,:,k) = src(:,:,k)*adj(k)
      end do
      jdlast = jday
      call closeunits(mon_units,nmons)
      write(6,*) trname(nt),'Sources interpolated to current day'
C****
      return
      end subroutine read_CO2_sources

      SUBROUTINE get_14CO2_IC(CO2IJL)
!@sum GET_14CO2_IC Calculates initial distribution for 14CO2 tracer
!@auth J.Lerner (modified from program by G. Russell)
C**** NOTE: tracer is supposed to start on 10/16
C**** October 1963 14CO2 Concentrations for GCM  2/26/99
C**** 2/2/2: generalized code for modelE
C****        
      USE MODEL_COM, ONLY: im,jm,lm,ls1,sige,psf,ptop
      USE DYNAMICS, only: pedn
      USE FILEMANAGER, only: openunit,closeunit
      IMPLICIT NONE
      PARAMETER (kmw=60)
      REAL*4 CO2W(37,0:30) ! ,AMASS(IM,JM)
      real*8 p(0:60),CO2JK(JM,0:kmw),CO2IJL(IM,JM,LM)
      CHARACTER*80 TITLE
      integer i,j,jw,k,l,n,iu_in,iu_out,kmw
      real*8 pup,cup,pdn,cdn,psum,csum,psurf,ptrop,sig,w,zk !,stratm

C****
C**** Read in CO2 concentrations from workshop
C****
c     OPEN (1,FILE='workshop.14co2',STATUS='OLD')
      call openunit('14CO2_IC_DATA',iu_in,.false.)
      DO 110 N=1,10
  110 READ (iu_in,911)
      READ (iu_in,911) ((CO2W(J,K),J=1,13),K=0,30)
      READ (iu_in,911)
      READ (iu_in,911)
      READ (iu_in,912) ((CO2W(J,K),J=14,25),K=0,30)
      READ (iu_in,912)
      READ (iu_in,912)
      READ (iu_in,912) ((CO2W(J,K),J=26,37),K=0,30)
      call closeunit(iu_in)
C**** Calculate workshop pressure levels (Pa)
      DO K=0,kmw
        ZK = 2.*K
        P(K) = 100000.*10**(-ZK/16.)
      end do
C****
C**** Interpolate workshop CO2W to CO2JK on GCM latitudes
C****
      DO K=31,kmw       ! above 17 Pa set all equal
        CO2JK( 1,K) = 250.
        CO2JK(JM,K) = 250.
        DO J=2,JM-1
          CO2JK(J,K) = 250.
      end do; end do
      DO K=0,30
        CO2JK( 1,K) = CO2W( 1,K)
        CO2JK(JM,K) = CO2W(37,K)
        DO J=2,JM-1
          W = 1. + (J-1)*36./(JM-1)
          JW=W
          CO2JK(J,K) = CO2W(JW,K)*(JW+1-W) + CO2W(JW+1,K)*(W-JW)
      end do; end do
C****
C**** Read in GCM atmospheric mass for October
C****
c     OPEN (2,FILE='AIJX011.O',FORM='UNFORMATTED',STATUS='OLD')
   !  call openunit('OCTAirMass',iu_in)
   !  READ  (iu_in)
   !  READ  (iu_in) TITLE,AMASS  !! in kg/m**2
   !  call closeunit(iu_in)
   !  WRITE (6,*) 'Read: ',TITLE
C****
C**** Interpoate CO2J to CO2IJL on GCM grid boxes conserving vertical
C**** means
C****
C**** psf, ptrop, pdn ..... in pascals (mb*100)
      ptrop = ptop*100.
    ! STRATM = 101.9368   ! (kg/m**2)(10 mb x 100)/grav
      DO 440 J=1,JM
      DO 440 I=1,IM
    ! PDN = (AMASS(I,J)*SIGE(1)+STRATM)*GRAV
      PDN = pedn(1,i,j)*100.
      psurf = pdn
      CDN = CO2JK(J,0)
      K=1
      DO 430 L=1,LM
      PSUM = 0.
      CSUM = 0.
      if (l.eq.ls1) psurf = psf*100.
      PUP  =  SIGE(L+1)*(psurf-ptrop)+ptrop
  410 IF(P(K).LE.PUP)  GO TO 420
      PSUM = PSUM +  PDN-P(K)
      CSUM = CSUM + (PDN-P(K))*(CDN+CO2JK(J,K))/2.
      PDN  = P(K)
      CDN  = CO2JK(J,K)
      K=K+1
      if (k.gt.kmw) stop ' Please increase kmw in get_14CO2_IC'
      GO TO 410
C****
  420 CUP  = CO2JK(J,K) + (CO2JK(J,K-1)-CO2JK(J,K))*(PUP-P(K))/
     /       (P(K-1)-P(K))
      PSUM = PSUM +  PDN-PUP
      CSUM = CSUM + (PDN-PUP)*(CDN+CUP)/2.
      CO2IJL(I,J,L) = CSUM/PSUM
      PDN = PUP
  430 CDN = CUP
  440 CONTINUE
C****
C**** Scale data to proper units (10**-18 kg 14CO2/kg air)
C****
      CO2IJL(:,:,:) = CO2IJL(:,:,:)*4.82d0*(14.+16.*2.)/29.029d0
      RETURN
C****
  911 FORMAT (5X,13F5.0)
  912 FORMAT (5X,12F5.0)
      END SUBROUTINE get_14CO2_IC

#endif
#endif

