#include "rundeck_opts.h"
      SUBROUTINE cheminit
!@sum cheminit initialize model chemistry
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on cheminit0C5_M23p & ds4p_chem_init_M23.f)
!@calls jplrts,phtlst,inphot,wave,reactn
c
C**** GLOBAL parameters and variables:
C
      USE FILEMANAGER, only: openunit,closeunit
      USE MODEL_COM, only: Itime, ItimeI, IM, LM
      USE TRCHEM_Shindell_COM, only:nc,ny,numfam,JPPJ,nn,ks,nps,nds,
     &                          ndnr,kps,kds,kpnr,kdnr,nnr,kss,nr,npnr,
     &                          nr2,nr3,nmm,nhet,prnls,prnrts,prnchg,
     &                          lprn,jprn,iprn,ay,nss,pHOx,pOx,pNOx,
     &                          yCH3O2,yC2O3,yROR,yXO2,yAldehyde,yNO3,
     &                          yRXPAR,yXO2N
#ifdef SHINDELL_STRAT_CHEM
     &                          ,pCLOx,pCLx,pOClOx,pBrOx,yCl2,yCl2O2
#endif
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
!@var iu_data temporary unit number
!@var title read in from file
!@var i,l loop dummy
      integer iu_data,i,l
      character*80 title
C
C     Read chem diagnostics parameters and molecule names
C     from MOLEC file:
C
      call openunit('MOLEC',iu_data,.false.,.true.)
      read(iu_data,100)prnls,prnrts,prnchg,lprn,jprn,iprn
 100  format(/3(50x,l1/),3(50x,i8/))
      read(iu_data,110)ay
#ifdef SHINDELL_STRAT_CHEM
 110  format(5(///10(a8)),(///3(a8)))
#else
 110  format(3(///10(a8)),(///5(a8)))
#endif
      call closeunit(iu_data)
C
C     Read JPL chemical reactions/rates from unit JPLRX:
C
      call jplrts
C
C     Read photolysis parameters and reactions from unit JPLPH :

      call phtlst
C
c     fastj initialization routine
      call inphot
c
c     Set up wavelengths 200 -730 nm & O2 and O3 cross sections
      call wave
C
c     Set up arrays of reaction numbers involving each molecule.
      call reactn
C
C     Initialize a few (IM,JM,LM) arrays, first hour only:
      IF(Itime.eq.ItimeI) THEN
        pHOx     =1.d0
        pOx      =1.d0
        pNOx     =1.d0
        yCH3O2   =1.d5
        yC2O3    =0.d0
        yROR     =0.d0
        yXO2     =0.d0
        yAldehyde=0.d0
        yNO3     =0.d0
        yXO2N    =0.d0
        yRXPAR   =0.d0
#ifdef SHINDELL_STRAT_CHEM
        pClOx    =1.d0
        pClx     =0.d0
        pOClOx   =0.d0
        pBrOx    =1.d0
        yCl2     =0.d0
        yCl2O2   =0.d0
#endif
      END IF
C
      return
      END SUBROUTINE cheminit


c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE jplrts
!@sum jplrts read/set up chemical reaction rates from JPL
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on cheminit0C5_M23p & ds4p_chem_init_M23)
!@calls lstnum
c
C**** GLOBAL parameters and variables:
C
      USE FILEMANAGER, only: openunit,closeunit
      USE TRCHEM_Shindell_COM, only: nr,nr2,nr3,nmm,nhet,pe,ea,nst,ro,
     &                          r1,sn,sb,nn,nnr,nc
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
C
!@var ate temporary reactants names array
!@var i,ii,j dummy loop variable
!@var iu_data temporary unit number
      CHARACTER*8, DIMENSION(4) :: ate
      INTEGER i,ii,j,iu_data
C
C     Read in the number of each type of reaction:
      call openunit('JPLRX',iu_data,.false.,.true.)
      read(iu_data,124)nr,nr2,nr3,nmm,nhet
      write(6,*)
      write(6,*) 'Chemical reactions used in the model: '
C
      do i=1,nr               ! >>> begin loop over total reactions <<<
       if(i.le.nr-nhet) then !non-hetero
         if(i.le.nr2) then   !mono or bi
           if(i.gt.nr2-nmm) then ! read monomolecular reactions
             if(i.eq.nr2-nmm+1)read(iu_data,22)ate
             read(iu_data,16)ate,pe(i),ea(i),nst(i-nr2+nmm)
             write(6,30) i,ate(1),' + ',ate(2),
     *       ' --> ',ate(3),' + ',ate(4)
           else                  ! read bimolecular reactions
   5         read(iu_data,16)ate,pe(i),ea(i)
             write(6,30) i,ate(1),' + ',ate(2),
     *       ' --> ',ate(3),' + ',ate(4)
           end if
         else                    ! read trimolecular reactions
 20        if(i.eq.nr2+1)read(iu_data,22)ate
           ii=i-nr2
           read(iu_data,21)ate,ro(ii),sn(ii),r1(ii),sb(ii)
           write(6,30) i,ate(1),' + ',ate(2),
     *     ' --> ',ate(3),' + ',ate(4)
         end if
       else                     ! read heterogeneous reactions
         if(i.eq.nr-(nhet-1))read(iu_data,22)ate
         read(iu_data,31)ate
         write(6,30) i,ate(1),' + ',ate(2),
     *   ' --> ',ate(3),' + ',ate(4)
       end if ! (i.le.nr-nhet)
c
      do j=1,2
        call lstnum(ate(j),nn(j,i))
        call lstnum(ate(j+2),nnr(j,i))
      end do
      end do                ! >>> end loop over total reactions <<<
C
 124  format(///5(/43x,i3)///)
  27  format(/(30x,i2))
  21  format(4x,a8,1x,a8,3x,a8,1x,a8,e8.2,f5.2,e9.2,f4.1)
  22  format(/10x,4a8/)
  25  format(//32x,2f7.1,i6)
  16  format(4x,a8,1x,a8,3x,a8,1x,a8,e8.2,f8.0,i4)
  31  format(4x,a8,1x,a8,3x,a8,1x,a8)
  30  format(1x,i3,2x,a8,a3,a8,a5,a8,a3,a8)
      call closeunit(iu_data)
      return
      end SUBROUTINE jplrts
C
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE lstnum(at,ks)
!@sum lstnum find molecule number in param list of molecules
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on cheminit0C5_M23p & ds4p_chem_init_M23)
c
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: nc,ay
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
C
!@var at local copy of species name
!@var ks local variable to be passed back to jplrts nnr or nn array.
!@var j dummy loop variable
      INTEGER j
      INTEGER,     INTENT(OUT) :: ks
      CHARACTER*8, INTENT(IN)  :: at
C     This code could probably be cleaned up in some kind of while-loop:
      j=0
   1  j=j+1
      if(j.gt.nc)goto2
      if(at.ne.ay(j))goto1
      ks=j
      goto3
   2  ks=nc+1
C
   3  continue
      return
      end SUBROUTINE lstnum
C
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE phtlst
!@sum phtlst read Photolysis Reactions and parameters
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on cheminit0C5_M23p & ds4p_chem_init_M23)
!@calls lstnum
c
C**** GLOBAL parameters and variables:
C
      USE FILEMANAGER, only: openunit,closeunit
      USE TRCHEM_Shindell_COM, only: nss,ks,kss,nc,JPPJ
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
C
!@var ate species name
!@var nabs,al,nll,nhu,o2up,o3up currently read from JPLPH, but not used
!@var i,j dummy loop variables
!@var iu_data temporary unit number
      CHARACTER*8, DIMENSION(3) :: ate
      INTEGER nabs,nll,nhu,i,j,iu_data
      REAL*8  al,o2up,o3up
C
C     Read in photolysis parameters:
      call openunit('JPLPH',iu_data,.false.,.true.)
      read(iu_data,121)nss,nabs,al,nll,nhu,o2up,o3up
C     check on the number of photolysis reactions:
      IF(nss.ne.JPPJ)
     & call stop_model('WARNING: nss.ne.JPPJ, check # photo rxns',255)
C
c     assign ks and kss gas numbers of photolysis reactants from list
      write(6,*) ' '
      write(6,*) 'Photolysis reactions used in the model: '
      do i=1,JPPJ
        read(iu_data,112)ate
        write(6,172) i,ate(1),' + hv   --> ',ate(2),' + ',ate(3)
        call lstnum(ate(1),ks(i))
        do j=2,3
           call lstnum(ate(j),kss(j-1,i))
        end do
      end do
 121  format(//2(45x,i2/),43x,f4.2/44x,i3/45x,i2/2(40x,e7.1/))
 112  format(4x,a8,3x,a8,1x,a8)
 172  format(1x,i2,2x,a8,a12,a8,a3,a8)
      call closeunit(iu_data)
C
      return
      end SUBROUTINE phtlst
C
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE reactn
!@sum reactn read chemical and photochemical reaction lists
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on cheminit0C5_M23p & ds4p_chem_init_M23)
!@calls guide,printls
c
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: nps,nds,kps,kds,nn,nnr,nr,ks,kss,
     &                          JPPJ,npnr,ndnr,kpnr,kdnr,prnls
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
C     (none)
C
c     Chemical reaction lists
      call guide(npnr,ndnr,kpnr,kdnr,nn,nnr,2,nr)
c     Photolysis reaction lists
      call guide( nps, nds, kps, kds,ks,kss,1,JPPJ)
C
C     print out some diagnostics:
      if(prnls) call printls
      return
      end SUBROUTINE reactn
C
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE guide(npr,ndr,kpr,kdr,xx,nnn,ns,nre)
!@sum guide read chemical and photochemical reaction lists
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on cheminit0C5_M23p & ds4p_chem_init_M23)
!@calls calcls
C
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: p_1,p_2,p_3,p_4
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
C
!@var xx   = either nn  or ks   from reactn sub
!@var nnn  = either nnr or kss  from reactn sub
!@var kpr  = either kps or kpnr from reactn sub
!@var kdr  = either kds or kdnr from reactn sub
!@var npr  = either nps or npnr from reactn sub
!@var ndr  = either nds or ndnr from reactn sub
!@var ns   = either 1   or    2 from reactn sub
!@var nre number of reactions
      INTEGER,  DIMENSION(p_4)     :: kpr, kdr
      INTEGER,  DIMENSION(p_3)     :: npr, ndr
      INTEGER, DIMENSION(p_1,p_2)  :: xx, nnn
      INTEGER ns, nre
C
c     Chemical and photolytic destruction
      call calcls(xx,ns,nnn,2,ndr,kdr,nre)
c     Chemical and photolytic production
      call calcls(nnn,2,xx,ns,npr,kpr,nre)
      return
      end SUBROUTINE guide
C
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE calcls(nn,ns,nnn,nns,ndr,kdr,nre)
!@sum calcls Set up reaction lists for calculated gases (1 to ny)
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on cheminit0C5_M23p & ds4p_chem_init_M23)
C
C Note: code could use more modernizations (e.g. the goto 5's):
C
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: ny, numfam, p_2, p_3, p_4, nfam,
     &                               prnls
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
C
!@var kdr  = either kdr or kpr from guide sub
!@var ndr  = either ndr or npr from guide sub
!@var nre  number of reactions
!@var nns number of partic_ on opposite side of reaction
!@var ns   = either ns or   2 from guide sub
!@var nn   = either xx or nnn from guide sub
!@var nnn  = either xx or nnn from guide sub
!@var ii,k,j,i,ij,i2,newfam,ifam dummy variables
      INTEGER, DIMENSION(p_4)    :: kdr
      INTEGER, DIMENSION(p_3)    :: ndr
      INTEGER nre, nns, ns, k, j, i, ij, i2, newfam, ifam, ii
      INTEGER, DIMENSION(ns,p_2) :: nn 
      INTEGER, DIMENSION(nns,p_2):: nnn

      k=1
      nfam(numfam+1)=ny+1
      do 5 j=1,numfam      !families, list only interfamily reactions
        newfam=0
        kdr(j)=k
        do 5 i=1,nre       ! 1 to # chem or phot reactions
          do 5 ij=1,ns    !ns # partic (prod & chem dest=2,phot dest=1)
c           check if molecule # nn() is element of family j
            if(nn(ij,i).ge.nfam(j).and.nn(ij,i).lt.nfam(j+1))then
c             check if reaction is intrafamily
              do i2=1,nns  ! nns # partic on opposite side of reac.
                if(nnn(i2,i).ge.nfam(j).and.nnn(i2,i).lt.nfam(j+1))goto5
              enddo
c             don't write same reaction twice
              if(k.eq.1)goto2
              if(ndr(k-1).eq.i.and.newfam.ne.0)goto5
 2            ndr(k)=i
              k=k+1
              newfam=1
            endif
 5    continue
C
      do 10 j=numfam+1,nfam(1)-1     ! individual non-family molecules
        kdr(j)=k
        do 10 i=1,nre                ! 1 to # chem or phot reactions
          do 10 ij=1,ns  !ns # partic (prod & chem dest=2,phot dest=1)
            if(nn(ij,i).ne.j)goto10    ! nn is mol # of participant
            ndr(k)=i
            k=k+1
  10  continue
C
      do 100 j=nfam(1),ny !indiv family mols.,list only intrafamily
        do ii=1,numfam-1
          if(j.lt.nfam(ii+1))then
            ifam=ii
            goto 110
          endif
        enddo
        ifam=numfam
 110    kdr(j)=k
        do 100 i=1,nre          ! 1 to # chem or phot reactions
          do 100 ij=1,ns  !ns # partic (prod & chem dest=2,phot dest=1)
            if(nn(ij,i).ne.j)goto100       ! nn is mol # of participant
c           check that reaction is intrafamily
            do i2=1,nns  ! nns # participants on opposite side of reac.
              if(nnn(i2,i).ge.nfam(ifam).and.nnn(i2,i).lt.nfam(ifam+1))
     &        then
                 ndr(k)=i
                 k=k+1
                 goto100
              endif
            enddo
 100  continue
      kdr(ny+1)=k
      if(prnls)write(*,*) 'nn array size :',k

      return
      end SUBROUTINE calcls
C
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE printls
!@sum printls print out some chemistry diagnostics (reaction lists)
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on cheminit0C5_M23p & ds4p_chem_init_M23)
C
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: kpnr,npnr,kdnr,ndnr,kps,nps,
     &                         ny,nn,nnr,ks,kss,ay,kds,nds,nc
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
C
!@var ireac,igas,ichange,ii dummy variables
      INTEGER ireac,igas,ichange,ii
C
c     print reaction lists
      write(6,*)
      write(6,*) '______________ CHEMICAL PRODUCTION _______________'
      ireac=0
      do igas=1,ny
        write(6,*)
        write(6,10) ay(igas)
        ichange=kpnr(igas+1)-kpnr(igas)
        if(ichange.ge.1) then
          do ii=1,ichange
            ireac=ireac+1
            if (nnr(2,npnr(ireac)).gt.nc) then
               write(6,20) ' Reaction # ',npnr(ireac),' produces ',
     *         ay(nnr(1,npnr(ireac))),' and  ','X'
            else
               write(6,20) ' Reaction # ',npnr(ireac),' produces ',
     *         ay(nnr(1,npnr(ireac))),' and  ',ay(nnr(2,npnr(ireac)))
            end if
          enddo
        end if
      end do
      write(6,*)
      write(6,*) '______________ CHEMICAL DESTRUCTION _______________'
      ireac=0
      do igas=1,ny
        write(6,*)
        write(6,10) ay(igas)
        ichange=kdnr(igas+1)-kdnr(igas)
        if(ichange.ge.1) then
          do ii=1,ichange
            ireac=ireac+1
            write(6,20) ' Reaction # ',ndnr(ireac),' destroys ',
     *      ay(nn(1,ndnr(ireac))),' and  ',ay(nn(2,ndnr(ireac)))
          enddo
        end if
      end do
      write(6,*)
      write(6,*) '______________ PHOTOLYTIC PRODUCTION _______________'
      ireac=0
      do igas=1,ny
        write(6,*)
        write(6,10) ay(igas)
        ichange=kps(igas+1)-kps(igas)
        if(ichange.ge.1) then
          do ii=1,ichange
            ireac=ireac+1
            write(6,20) ' Reaction # ',nps(ireac),' produces ',
     *      ay(kss(1,nps(ireac))),' and  ', ay(kss(2,nps(ireac)))
          enddo
        end if
      end do
      write(6,*)
      write(6,*) '______________ PHOTOLYTIC DESTRUCTION _______________'
      ireac=0
      do igas=1,ny
        write(6,*)
        write(6,10) ay(igas)
        ichange=kds(igas+1)-kds(igas)
        if(ichange.ge.1) then
          do ii=1,ichange
            ireac=ireac+1
            write(6,30) ' Reaction # ',nds(ireac),' destroys ',
     *      ay(ks(nds(ireac)))
          enddo
        end if
      end do
      return
  10  format(1x,a8)
  20  format(a12,i4,a10,a8,a6,a8)
  30  format(a12,i4,a10,a8)
C
      end SUBROUTINE printls
C
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE wave
!@sum wave Set up Wavelengths 200-730 nm, and O2 & O3 Cross Sections
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on cheminit0C5_M23p & ds4p_chem_init_M23)
C
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: n_bnd3, wlt, sech, sO3, sO2
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
C
!@var lw dummy variable
      INTEGER lw
C
      do lw=1,n_bnd3
          wlt(lw)=195.d0+5.d0*real(lw)
c         sech(2,x)=O3 cross sections, sech(1,x)=O2 cross sections
          sech(2,lw)=sO3(lw)
          if(lw.le.9) then
            sech(1,lw)=sO2(lw)
          end if
      end do
      return
      end SUBROUTINE wave
C
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE inphot
!@sum inphot initialise photolysis rate data, called directly from the
!@+   cinit routine in ASAD. Currently use to read the JPL spectral data
!@+   and standard O3 and T profiles and to set the appropriate reaction
!@+   index.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on cheminit0C5_M23p & ds4p_chem_init_M23)
!@calls RD_TJPL,RD_PROF
C
C**** GLOBAL parameters and variables:
C
      USE FILEMANAGER, only: openunit,closeunit
      USE TRCHEM_Shindell_COM, only: jfacta, jlabel,jppj
#ifdef SHINDELL_STRAT_CHEM
     &                  ,MXFASTJ,MIEDX2,title_aer_pf,NAA
#endif
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
C
!@var ipr Photolysis reaction counter
!@var cline dummmy text
!@var i dummy loop variable
!@var iu_data temporary unit number
!@var temp1 temp variable to read in jfacta
!@var temp2 temp variable to read in jlabel
      integer iu_data, ipr, i
      character*120 cline
      real*8 temp1
      character*7 temp2
c
c Reread the ratj_GISS.d file to map photolysis rate to reaction
c Read in quantum yield jfacta and fastj label jlabel
      ipr=0
      call openunit('RATJ',iu_data,.false.,.true.)
 10   read(iu_data,'(a)',err=20) cline
      if(cline(2:5).eq.'9999') then
         go to 20
      elseif(cline(1:1).eq.'#' .or. cline(5:5).eq.'$') then
         go to 10
      else
         ipr=ipr+1
         backspace iu_data
         read(iu_data,'(78x,f5.1,2x,a7)',err=20) temp1,temp2
         jfacta(ipr) = temp1
         jlabel(ipr) = temp2
         jfacta(ipr)=jfacta(ipr)*1.d-2 
         go to 10
      endif
 20   call closeunit(iu_data)
      if(ipr.ne.jppj) then
         write(6,1000) ipr,jppj
         call stop_model('problem with # photolysis reactions',255)
      endif
 1000 format(' Error: ',i3,' photolysis labels but ',i3,' reactions')
c
c Print details
      write(6,1100) ipr
      write(6,1200) (i, jlabel(i), jfacta(i),i=1,ipr)
c
 1100 format(' Fast-J Photolysis Scheme: considering ',i2,' reactions')
 1200 format(3x,10(3(i2,': ',a7,' (Q.Y. ',f5.3,') '),/,3x))
c
c Read in JPL spectral data set
      CALL openunit('SPECFJ',iu_data,.false.,.true.)
      CALL RD_TJPL(iu_data)
      call closeunit(iu_data)
c
c Read in T & O3 climatology
      CALL openunit('ATMFJ',iu_data,.false.,.true.)
      CALL RD_PROF(iu_data)
      call closeunit(iu_data)
c
#ifdef SHINDELL_STRAT_CHEM
c  Ensure all aerosol types are valid selections
      do i=1,MXFASTJ
        write(6,1001) MIEDX2(i),title_aer_pf(MIEDX2(i))
        if(MIEDX2(i).gt.NAA.or.MIEDX2(i).le.0) then
          write(6,1201) MIEDX2(i),NAA
          call stop_model('Problem in inphot.',255)
        endif
      enddo
 1001 format('Using Aerosol type: ',i2,1x,a)
 1201 format('Aerosol type ',i1,' unsuitable; supplied values must be',
     $       ' between 1 and ',i1)
#endif
      return
      end SUBROUTINE inphot
C
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE RD_TJPL(NJ1)
!@sum RD_TJPL Read wavelength bins, solar fluxes, Rayleigh parameters,
!@+   T-dependent cross sections and Rayleigh/aerosol scattering phase
!@+   functions with temperature dependences. Current data originates
!@+   from JPL'97.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on cheminit0C5_M23p & ds4p_chem_init_M23)
C
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: NJVAL,WBIN,WL,NWWW,FL,QRAYL,QBC,
     &                          Q1D,TQQ,QQQ,NAA,QAAFASTJ,NK,WAAFASTJ,
     &                          PAA,zpdep,npdep,jpdep,lpdep,NS,TITLE0,
     &                          NW1,NW2,TITLEJ,JPPJ,jind,jlabel,jfacta,
     &                          title_aer_pf,QO2,QO3
#ifdef SHINDELL_STRAT_CHEM
     &                          ,SSA,RAA,NP
#endif
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
C
!@var NJ1 local copy of unit number to read
!@var i,j,k,iw dummy loop variables
!@var jj dummy variable
!@var nQQQ minus no. additional J-values from X-sects (O2,O3P,O3D+NQQQ)
!@var NJVAL2 temporary test for NJVAL= its constant value...
      INTEGER, INTENT(IN) :: NJ1
      INTEGER i,j,k,iw,jj,nqqq,NJVAL2
C
      TQQ = 0.d0

C------read in-------spectral data------------------------------------
      READ(NJ1,'(A)') TITLE0
      WRITE(6,'(1X,A)') TITLE0
      READ(NJ1,'(10X,4I5)') NJVAL2,NWWW,NW1,NW2
      IF(NJVAL.ne.NJVAL2) THEN
        WRITE(6,*)'NJVAL (constant)= ',NJVAL,' but it is ', NJVAL2,
     &  'when read in from SPECFJ file.  Please reconcile.'
        call stop_model('NJVAL problem in RD_TJPL',255)
      END IF
      NQQQ = NJVAL-3
      READ(NJ1,102) (WBIN(IW),IW=1,NWWW)
      READ(NJ1,102) (WBIN(IW+1),IW=1,NWWW)
      READ(NJ1,102) (WL(IW),IW=1,NWWW)
      READ(NJ1,102) (FL(IW),IW=1,NWWW)
      READ(NJ1,102) (QRAYL(IW),IW=1,NWWW)
      READ(NJ1,102) (QBC(IW),IW=1,NWWW)   !From Loiusse et al[JGR,96]
c
C-Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields(each at 3 temps)
      DO K=1,3
        READ(NJ1,103) TITLEJ(K,1),TQQ(K,1), (QO2(IW,K),IW=1,NWWW)
      ENDDO
      DO K=1,3
        READ(NJ1,103) TITLEJ(K,2),TQQ(K,2), (QO3(IW,K),IW=1,NWWW)
      ENDDO
      DO K=1,3
        READ(NJ1,103) TITLEJ(K,3),TQQ(K,3), (Q1D(IW,K),IW=1,NWWW)
      ENDDO
      do k=1,3
        write(6,200) titlej(1,k),(tqq(i,k),i=1,3)
      enddo
c
C---Read remaining species:  X-sections at 2 T's
      DO J=1,NQQQ
        READ(NJ1,103) TITLEJ(1,J+3),TQQ(1,J+3),(QQQ(IW,1,J),IW=1,NWWW)
        READ(NJ1,103) TITLEJ(2,J+3),TQQ(2,J+3),(QQQ(IW,2,J),IW=1,NWWW)
        write(6,200) titlej(1,j+3),(tqq(i,j+3),i=1,2)
      ENDDO
      READ(NJ1,'(A)') TITLE0
c
c---Pressure dependencies
#ifdef SHINDELL_STRAT_CHEM
      npdep=0
#else
      read(NJ1,104) npdep
      do k=1,npdep
        read(NJ1,105) lpdep(k),(zpdep(iw,k),iw=1,nwww)
        write(6,201)  lpdep(k),(zpdep(iw,k),iw=1,nwww)
      enddo
      read(NJ1,'(A)') TITLE0
#endif
c
c---Zero index arrays
      jind=0
      jpdep=0
c
C---Set mapping index
      do j=1,NJVAL
        do k=1,jppj
          if(jlabel(k).eq.titlej(1,j)) jind(k)=j
        enddo
#ifndef SHINDELL_STRAT_CHEM
        do k=1,npdep
          if(lpdep(k).eq.titlej(1,j)) jpdep(j)=k
        enddo
#endif
      enddo
      do k=1,jppj
        if(jfacta(k).eq.0.d0)
     &             write(6,*) 'Not using photolysis reaction ',k
        if(jind(k).eq.0) then
          if(jfacta(k).eq.0.d0) then
            jind(k)=1
          else
            write(6,*) 'Which J-rate for photolysis reaction ',k,' ?'
            call stop_model('J-rate problem in RD_TJPL',255)
          endif
        endif
      enddo
c
C---Read aerosol phase functions:
#ifdef SHINDELL_STRAT_CHEM
      read(NJ1,'(A10,I5,/)') TITLE0,NAA
      if(NAA.gt.NP)then 
        write(6,350) NAA
        call stop_model('NAA too large in RD_TJPL',255)
      endif
      
      NK=4        ! Fix number of wavelengths at 4
      do j=1,NAA
        read(NJ1,110) title_aer_pf(j)
        do k=1,NK
          read(NJ1,*) WAAFASTJ(k,j),QAAFASTJ(k,j),RAA(k,j),SSA(k,j),
     &                (PAA(i,k,j),i=1,8)
        enddo
      enddo
#else
      READ(NJ1,'(A10,I5)') TITLE0,NAA
      write(6,*)'Title0 is',Title0
      write(6,*)'NAA is',NAA

      DO J=1,NAA
        READ(NJ1,'(A5,I3,I2,14F5.0)') title_aer_pf(J),JJ,NK,
     $            (WAAFASTJ(K,J),QAAFASTJ(K,J),K=1,NK)
        DO K=1,NK
          READ(NJ1,'(8X,8F9.5)') (PAA(I,K,J), I=1,8)
        ENDDO
      ENDDO
#endif   
      write(6,*) 'Aerosol phase functions & wavelengths'
      DO J=1,NAA
        write(6,'(1x,A8,I2,A,9F8.1)')
     $             title_aer_pf(J),J,'  wavel=',(WAAFASTJ(K,J),K=1,NK)
        write(6,'(9x,I2,A,9F8.4)') J,'  Qext =',(QAAFASTJ(K,J),K=1,NK)
      ENDDO   
C--------
  101 FORMAT(8E10.3)
#ifdef SHINDELL_STRAT_CHEM
  102 FORMAT((10X,6E10.3)/(10X,6E10.3)/(10X,6E10.3))
  103 FORMAT(A7,F3.0,6E10.3/(10X,6E10.3)/(10X,6E10.3))
#else
  102 FORMAT(10X,7E10.3)
  103 FORMAT(A7,F3.0,7E10.3)
#endif
  104 FORMAT(13x,i2)
  105 FORMAT(A7,3x,7E10.3)
  110 format(3x,a20)
  200 format(1x,' x-sect:',a10,3(3x,f6.2))
  201 format(1x,' pr.dep:',a10,7(1pE10.3))
  350 format(' Too many phase functions supplied; increase NP to ',i2)
      RETURN
      END SUBROUTINE RD_TJPL
C
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE rd_prof(nj2)
!@sum rd_prof input T & O3 reference profiles, define Black Carbon prof.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@ver  1.0 (based on cheminit0C5_M23p & ds4p_chem_init_M23)
C
C**** GLOBAL parameters and variables:
C
      USE TRCHEM_Shindell_COM, only: TITLE0,
#ifdef SHINDELL_STRAT_CHEM
     & TREF2, OREF2, BREF2, ZZHT
#else 
     & TREF, OREF, BREF
#endif
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
C
!@var nj2 local unit number
!@var ia,i,m,l,lat,mon,ntlats,ntmons,n216 local dummy variables
      INTEGER, INTENT(IN) :: nj2
      integer ia, i, m, l, lat, mon, ntlats, ntmons, n216
#ifdef SHINDELL_STRAT_CHEM
      REAL*8 ofac, ofak
#endif
C
      READ(NJ2,'(A)') TITLE0
      WRITE(6,'(1X,A)') TITLE0
      READ(NJ2,'(2I5)') NTLATS,NTMONS
      WRITE(6,1000) NTLATS,NTMONS
      N216 = MIN0(216, NTLATS*NTMONS)
      DO IA=1,N216
        READ(NJ2,'(1X,I3,3X,I2)') LAT, MON
        M = MIN(12, MAX(1, MON))
        L = MIN(18, MAX(1, (LAT+95)/10))
#ifdef SHINDELL_STRAT_CHEM
        READ(NJ2,'(3X,11F7.1)') (TREF2(I,L,M), I=1,41)
        READ(NJ2,'(3X,11F7.4)') (OREF2(I,L,M), I=1,31)
#else
        READ(NJ2,'(3X,11F7.1)') (TREF(I,L,M), I=1,41)
        READ(NJ2,'(3X,11F7.4)') (OREF(I,L,M), I=1,31)
#endif
      ENDDO
C      
#ifdef SHINDELL_STRAT_CHEM      
c
c  Extend climatology to 100 km
      ofac=exp(-2.d5/ZZHT)
      do i=32,51
        ofak=ofac**(i-31)
        do m=1,ntmons
          do l=1,ntlats
            oref2(i,l,m)=oref2(31,l,m)*ofak
          enddo
        enddo
      enddo
      do l=1,ntlats
        do m=1,ntmons
          do i=42,51
            tref2(i,l,m)=tref2(41,l,m)
          enddo
        enddo
      enddo
c
c Approximate Black Carbon up to 10 km; surface 200 ng/m3  (Liousse et al)
c Scale: 1 ng/m3 = 1.0d-15 g/cm3 (1.0d-11 g/m2/cm as BREF is in cm))
      do i=1,6;  BREF2(i) =10.d0*1.0d-11; end do
      do i=7,51; BREF2(i) =0.d0         ; end do
#else
      do i=1,11 ; BREF(i) =10.d0*1.0d-11; end do
      do i=12,41; BREF(i) =0.d0         ; end do
#endif
      return
 1000 format(1x,'Data: ',i3,' Lats x ',i2,' Months')
C
      end SUBROUTINE rd_prof
