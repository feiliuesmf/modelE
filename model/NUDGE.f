c *****************************************************************
      MODULE NUDGE_COM
!@sum  NUDGE_COM contains all the nudging related variables
!@auth 
!@ver  
      USE MODEL_COM, only : im,jm,lm
c      USE DOMAIN_DECOMP, only : grid
      IMPLICIT NONE
      SAVE
!@param  nlevnc vertical levels of NCEP data  
      INTEGER, PARAMETER :: nlevnc =17
!@var  U1, V1 NCEP wind at prior ncep timestep (m/s)
!@var  U2, V2 NCEP wind at the following ncep timestep (m/s)
      REAL*4, DIMENSION(IM,JM,LM) :: u1,v1,u2,v2
      REAL*4, DIMENSION(IM,JM-1,nlevnc) :: UN1,VN1,UN2,VN2
      REAL*4, DIMENSION(nlevnc) :: pl
!@var netcdf integer
      INTEGER :: ncidu,ncidv,uid,vid,plid
      INTEGER :: step_rea=1,zirk=0
!@var tau nuding time interpoltation
      REAL*8 :: tau
!@param  anudgeu anudgev relaxation constant
      REAL*8 :: anudgeu = 0., anudgev = 0.

      END MODULE NUDGE_COM

c------------------------------------------------------------------
      SUBROUTINE NUDGE(UGCM,VGCM,DTSTEP)
c******************************************************************
!@sum  Nudging of the horizontal wind comonents to reanalysis data sets
!@auth Susanne Bauer
!@ver  
      USE MODEL_COM, only : im,jm,lm,plbot
      USE DOMAIN_DECOMP, only : grid
      USE NUDGE_COM, only : u1,v1,u2,v2,tau,anudgeu,anudgev,pl,nlevnc
      IMPLICIT NONE
      REAL*8,DIMENSION(IM,JM,LM) ::  ugcm,vgcm
c     LOCAL
      INTEGER i,j,l
      REAL*8   alphau,alphav,a,dtstep

c - - - - - - - - - - - - - - - - - - - - - - - - - |
c   x_gcm = a * x_gcm + ( 1 - a ) * x_reanalys      |
C----------------------------------------------------
             alphau=dtstep * anudgeu  !alphau !anudgeu
             alphav=dtstep * anudgev  !alphav !nudgev

         do l=1,lm
           if (plbot(l+1).lt.pl(nlevnc)) cycle
         do j=2,jm  ! Please pay attention j starts at 2
         do i=1,im
            a=(1.-tau)*u1(i,j,l)+tau*u2(i,j,l)         !time interpolation
            ugcm(i,j,l) = (ugcm(i,j,l)+ (a * alphau))/ (1+alphau) !nudging
            a=(1.-tau)*v1(i,j,l)+tau*v2(i,j,l)         !time interpolation
            vgcm(i,j,l) = (vgcm(i,j,l)+ (a * alphav))/ (1+alphav) !nudging
         enddo
         enddo
         enddo

      RETURN
      END SUBROUTINE NUDGE

c------------------------------------------------------------------
      SUBROUTINE NUDGE_PREP
c******************************************************************
!@sum  Nudging of the horizontal wind comonents to reanalysis data sets
!@auth Susanne
!@ver  
      USE MODEL_COM, only: im,jm,lm,jhour,jday,itime,nday,jyear,iyear1
      USE NUDGE_COM
      IMPLICIT NONE
      include 'netcdf.inc'
      INTEGER i
      integer start(4),count(4),status
      character(len=3) :: nstr1,nstr2
C-----------------------------------------------------------------------

      if (jday.eq.1.and.mod(itime,nday).eq.0) then
            zirk = jyear - iyear1
      endif 
      if (step_rea.eq.1460.) then
            zirk = jyear - iyear1 + 1
c      endif
c -----------------------------------------------------------------
c   Opening of the files to be read
c -----------------------------------------------------------------

        if (zirk.lt.10) then
           write(nstr1,'(I1)') zirk
        elseif (zirk.lt.100) then
           write(nstr1,'(I2)') zirk
        endif   

        if (zirk+1.lt.10) then
           write(nstr2,'(I1)') zirk +1
        elseif (zirk+1.lt.100) then
           write(nstr2,'(I2)') zirk + 1
        endif   

      print*, '  I N    N U D G E : OPEN NF FILES','  u',nstr2,'.nc'

            status=NF_CLOSE('u'//trim(nstr1)//'.nc',NCNOWRIT,ncidu)
            status=NF_CLOSE('v'//trim(nstr1)//'.nc',NCNOWRIT,ncidv)

            status=NF_OPEN('u'//trim(nstr2)//'.nc',NCNOWRIT,ncidu)
            status=NF_OPEN('v'//trim(nstr2)//'.nc',NCNOWRIT,ncidv)

            status=NF_INQ_VARID(ncidu,'level',plid)
            status=NF_INQ_VARID(ncidu,'uwnd',uid)
            status=NF_INQ_VARID(ncidv,'vwnd',vid)

        endif  ! step = 1460
C-----------------------------------------------------------------------
      if (mod(itime,nday/4).eq.0.) then
           vn1(:,:,:) = vn2(:,:,:)
           un1(:,:,:) = un2(:,:,:)
           step_rea = INT( (((jday - 1) * 24) + jhour)/6) + 2
            if (step_rea.eq.1461.) step_rea = 1
           call read_ana(step_rea)
      endif
C-----------------------------------------------------------------------
      tau  = mod(itime,nday/4)/float(nday/4)

      call vinterana2mod(un1,nlevnc,pl,u1)
      call vinterana2mod(vn1,nlevnc,pl,v1)
      call vinterana2mod(un2,nlevnc,pl,u2)
      call vinterana2mod(vn2,nlevnc,pl,v2)

      RETURN
      END SUBROUTINE NUDGE_PREP

   
c -----------------------------------------------------------------
      SUBROUTINE READ_ANA(timestep)
c******************************************************************
!@sum  read in analysis data sets
!@auth Susanne Bauer
!@ver  
      USE MODEL_COM, only : im,jm,lm
      USE NUDGE_COM
      USE DOMAIN_DECOMP, only : grid
      IMPLICIT NONE
      include 'netcdf.inc'

      integer timestep,l
      integer start(4),count(4),status

c -----------------------------------------------------------------
c   read u, v, and pl
c -----------------------------------------------------------------
      status=NF_GET_VARA_REAL(ncidu,plid,1,nlevnc,pl)
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=timestep

      count(1)=im
      count(2)=jm-1
      count(3)=nlevnc
      count(4)=1

c  u
      status=NF_GET_VARA_REAL(ncidu,uid,start,count,un2)
c  v
      status=NF_GET_VARA_REAL(ncidv,vid,start,count,vn2)

      return
      end subroutine read_ana

c -----------------------------------------------------------------
      subroutine vinterana2mod(varo,lmo,po,varn)
c**********************************************************
!@sum  vertical interpolation
!@auth Susanne Bauer
!@ver  
      USE MODEL_COM, only : im,jm,lm
      USE DYNAMICS, only : PMID  ! Pressure at mid point of box (mb)
      USE DOMAIN_DECOMP, only : grid
       IMPLICIT NONE
c 
c ==============

       INTEGER lmo ! vertical dimensions of input
       INTEGER  i,j,lo,ln

        real po(lmo)! pressure levels in millibars of the input
        REAL varo(im,jm-1,lmo)!  Variable on the old grid (input)
        REAL varn(im,jm,lm) !  Variable on the new grid (output)
        real coef,dp1,dp2
 
        do j= 2,jm              ! Please pay attention j starts at 2
          do i=1,im
          do ln=1,lm

            lo=1
            if (pmid(ln,i,j).ge.po(1))then
              varn(i,j,ln) =  varo(i,j-1,1)
            else if (pmid(ln,i,j).le.po(lmo)) then
              varn(i,j,ln) =  varo(i,j-1,lmo)
            else
 10           if ( (pmid(ln,i,j).le.po(lo)).and.
     &             (pmid(ln,i,j).gt.po(lo+1)) )then
c                coef=(pmid(ln,i,j)-po(lo))/(po(lo+1)-po(lo))
c                varn(i,j,ln)=varo(i,j-1,lo)
c     &               +coef*(varo(i,j-1,lo+1)-varo(i,j-1,lo))

                dp1 = (-pmid(ln,i,j)+po(lo))
                dp2 = (pmid(ln,i,j)-po(lo+1))

                varn(i,j,ln)=((varo(i,j-1,lo) * dp1)
     &               + (varo(i,j-1,lo+1) *dp2)) / (dp1 + dp2)
              else
                lo=lo+1
                if (lo.le.lmo) goto 10
              end if
            endif

          enddo
        enddo
        enddo

      return
      end  subroutine vinterana2mod
c------------------------------------------------------------------
      SUBROUTINE NUDGE_INIT
c******************************************************************
!@sum  Initialization for Nudging
!@auth Susanne
!@ver  
      USE MODEL_COM, only : im,jm,lm,jhour,jday,itime,nday,jyear,iyear1
      USE NUDGE_COM
      USE PARAM
      IMPLICIT NONE
      include 'netcdf.inc'

      integer start(4),count(4),status
      character(len=3) :: nstr1
C**** Rundeck parameters:
      call sync_param( "ANUDGEU", ANUDGEU )
      call sync_param( "ANUDGEV", ANUDGEV )
c -----------------------------------------------------------------
c   Initialisation of the files to be read
c -----------------------------------------------------------------
           step_rea = INT( (((jday - 1) * 24) + jhour)/6) + 1
           print*,'READING REANALYSIS INIT ',jhour, jday, step_rea
C-----------------------------------------------------------------------

      if (jday.eq.1.and.mod(itime,nday).eq.0) then
            zirk = jyear - iyear1
      endif 
      if (step_rea.eq.1460.) then
            zirk = jyear - iyear1 + 1
      endif
c -----------------------------------------------------------------
c   Opening of the files to be read
c -----------------------------------------------------------------

        if (zirk.lt.10) then
           write(nstr1,'(I1)') zirk
        elseif (zirk.lt.100) then
           write(nstr1,'(I2)') zirk
        endif   

            status=NF_OPEN('u'//trim(nstr1)//'.nc',NCNOWRIT,ncidu)
            status=NF_OPEN('v'//trim(nstr1)//'.nc',NCNOWRIT,ncidv)

            status=NF_INQ_VARID(ncidu,'level',plid)
            status=NF_INQ_VARID(ncidu,'uwnd',uid)
            status=NF_INQ_VARID(ncidv,'vwnd',vid) 

C------------------------------------------------------------------
c -----------------------------------------------------------------
c   read u, v, and pl
c -----------------------------------------------------------------
      status=NF_GET_VARA_REAL(ncidu,plid,1,nlevnc,pl)
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=step_rea

      count(1)=im
      count(2)=jm-1
      count(3)=nlevnc
      count(4)=1

      status=NF_GET_VARA_REAL(ncidu,uid,start,count,un1)
      status=NF_GET_VARA_REAL(ncidv,vid,start,count,vn1)

      start(4)=step_rea+1
      status=NF_GET_VARA_REAL(ncidu,uid,start,count,un2)
      status=NF_GET_VARA_REAL(ncidv,vid,start,count,vn2)

      return
      end subroutine nudge_init


c -----------------------------------------------------------------
