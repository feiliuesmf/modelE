c *****************************************************************
      MODULE NUDGE_COM
!@sum  NUDGE_COM contains all the nudging related variables
!@auth 
!@ver  
      USE MODEL_COM, only : im,jm,lm
      IMPLICIT NONE
      SAVE
!@var  U1, V1 NCEP U wind at prior ncep timestep (m/s)
      REAL*4, DIMENSION(IM,JM,LM) :: U1,V1
!@var  U2, V2 NCEP U wind at the following ncep timestep (m/s)
      REAL*4, DIMENSION(IM,JM,LM) :: U2,V2
!@var netcdf integer
      INTEGER :: ncidu,ncidv,uid,vid,plid,step_rea=1
!@var tau nuding time interpoltation
      INTEGER :: tau
!@param  nlevnc vertical levels of NCEP data  
      INTEGER, PARAMETER :: nlevnc =17
!@param  anudgeu anudgev relaxation constant
      REAL*8 :: anudgeu = 0., anudgev = 0.

      END MODULE NUDGE_COM

c------------------------------------------------------------------
      SUBROUTINE NUDGE(UGCM,VGCM,DTSTEP)
c******************************************************************
!@sum  Nudging of the horizontal wind comonents to reanalysis data sets
!@auth Susanne
!@ver  
      USE MODEL_COM, only : im,jm,lm
      USE NUDGE_COM, only : U1,V1,U2,V2,tau,anudgeu,anudgev
      IMPLICIT NONE

      REAL*8 UGCM(IM,JM,LM),VGCM(IM,JM,LM)

c     LOCAL
      INTEGER i,j,l
      REAL  alphau,alphav,a,dtstep


c - - - - - - - - - - - - - - - - - - - - - - - - - |
c   x_gcm = a * x_gcm + ( 1 - a ) * x_reanalys      |
C----------------------------------------------------
             alphau=dtstep * anudgeu !alphau !anudgeu
             alphav=dtstep * anudgev !alphav !nudgev

         do l=1,lm
         do j=2,jm
         do i=1,im
            a=(1.-tau)*u1(i,j,l)+tau*u2(i,j,l)         !time interpolation
            ugcm(i,j,l) =(ugcm(i,j,l)+ a * alphau)/ (1+alphau) !nudging
            a=(1.-tau)*v1(i,j,l)+tau*v2(i,j,l)         !time interpolation
            vgcm(i,j,l) =(vgcm(i,j,l)+ a * alphav)/ (1+alphav) !nudging
         enddo
         enddo
         enddo

      RETURN
      END SUBROUTINE NUDGE

c------------------------------------------------------------------
      SUBROUTINE NUDGE_PREP
c******************************************************************
c(UNDG,VNDG)
!@sum  Nudging of the horizontal wind comonents to reanalysis data sets
!@auth Susanne
!@ver  
      USE MODEL_COM, only : im,jm,lm,jhour
      USE NUDGE_COM, only : U1,U2,V1,V2,tau,step_rea
      IMPLICIT NONE
      INTEGER i
C-----------------------------------------------------------------------
      if (jhour.eq.0.or.jhour.eq.6.or.jhour.eq.12.or.jhour.eq.18) then
           v1(:,:,:) = v2(:,:,:)
           u1(:,:,:) = u2(:,:,:)
           step_rea=step_rea+1
           print*,'READING REANALYSIS AT ',jhour,step_rea
           call read_ana(u2,v2,step_rea)
      endif

C-----------------------------------------------------------------------
      tau = jhour
      if(tau.gt.6) tau = tau - 6.
      if(tau.gt.6) tau = tau - 6.
      if(tau.gt.6) tau = tau - 6.
        tau=tau/6.  !time interpolation
      if(tau.eq.1.) tau = 0.  

      RETURN
      END SUBROUTINE NUDGE_PREP

   
c -----------------------------------------------------------------
      SUBROUTINE READ_ANA(UN,VN,timestep)
c******************************************************************
!@sum  read in analysis data sets
!@auth Susanne Bauer
!@ver  
      USE MODEL_COM, only : im,jm,lm
      USE NUDGE_COM, only : nlevnc,ncidu,ncidv,uid,vid,plid
      IMPLICIT NONE
      include 'netcdf.inc'

      integer timestep,l

      REAL*4 U(IM,JM-1,nlevnc),V(IM,JM-1,nlevnc),PL(nlevnc)
      REAL*4 UN(IM,JM,LM),VN(IM,JM,LM)

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
      status=NF_GET_VARA_REAL(ncidu,uid,start,count,u)
c  v
      status=NF_GET_VARA_REAL(ncidv,vid,start,count,v)

      call vinterana2mod(u,nlevnc,pl,un)
      call vinterana2mod(v,nlevnc,pl,vn)

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
       IMPLICIT NONE
c 
c ==============

       INTEGER lmo ! vertical dimensions of input
       INTEGER  i,j,lo,ln

        real po(lmo)! pressure levels in millibars of the input
        REAL varo(im,jm-1,lmo)!  Variable on the old grid (input)
        REAL varn(im,jm  ,lm) !  Variable on the new grid (output)
        real coef


        do i=1,im
        do j=2,jm
           do ln=1,lm
              if (pmid(ln,i,j).ge.po(1))then
                 varn(i,j,ln) =  varo(i,j-1,1)
              else if (pmid(ln,i,j).le.po(lmo)) then
                 varn(i,j,ln) =  varo(i,j-1,lmo)
              else
                 do lo=1,lmo-1 
                    if ( (pmid(ln,i,j).le.po(lo)).and.
     &                 (pmid(ln,i,j).gt.po(lo+1)) )then
                       coef=(pmid(ln,i,j)-po(lo))
     &                 /(po(lo+1)-po(lo))
                       varn(i,j,ln)=varo(i,j-1,lo)
     &                 +coef*(varo(i,j-1,lo+1)-varo(i,j-1,lo))
                    end if
                 enddo           
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
      USE MODEL_COM, only : im,jm,lm
      USE NUDGE_COM
      USE PARAM
      IMPLICIT NONE
      include 'netcdf.inc'

      REAL*4 U(IM,JM-1,nlevnc),V(IM,JM-1,nlevnc),PL(nlevnc)

      integer start(4),count(4),status
C**** Rundeck parameters:
      call sync_param( "ANUDGEU", ANUDGEU )
      call sync_param( "ANUDGEV", ANUDGEV )
c -----------------------------------------------------------------
c   Initialisation of the files to be read
c -----------------------------------------------------------------
            status=NF_OPEN('u.nc',NCNOWRIT,ncidu)
            status=NF_INQ_VARID(ncidu,'uwnd',uid)
            status=NF_INQ_VARID(ncidu,'level',plid)

            status=NF_OPEN('v.nc',NCNOWRIT,ncidv)
            status=NF_INQ_VARID(ncidv,'vwnd',vid)

c -----------------------------------------------------------------
c   read u, v, and pl
c -----------------------------------------------------------------
      status=NF_GET_VARA_REAL(ncidu,plid,1,nlevnc,pl)
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=1

      count(1)=im
      count(2)=jm-1
      count(3)=nlevnc
      count(4)=1

c  u
      status=NF_GET_VARA_REAL(ncidu,uid,start,count,u)
c  v
      status=NF_GET_VARA_REAL(ncidv,vid,start,count,v)

      call vinterana2mod(u,nlevnc,pl,u2)
      call vinterana2mod(v,nlevnc,pl,v2)

      return
      end subroutine nudge_init


c -----------------------------------------------------------------
