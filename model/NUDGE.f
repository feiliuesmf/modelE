c *****************************************************************
      MODULE NUDGE_COM
!@sum  NUDGE_COM contains all the nudging related variables
!@auth 
!@ver  
      USE MODEL_COM, only : im,jm,lm
      USE DOMAIN_DECOMP, only : grid
      IMPLICIT NONE
      SAVE
!@var  U1, V1 NCEP U wind at prior ncep timestep (m/s)
      REAL*4, DIMENSION(IM,JM,LM) :: U1,V1
!@var  U2, V2 NCEP U wind at the following ncep timestep (m/s)
      REAL*4, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) :: U2,V2
!@var netcdf integer
      INTEGER :: ncidu,ncidv,uid,vid,plid,step_rea=1,zirk=0
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
      USE DOMAIN_DECOMP, only : grid
      USE NUDGE_COM, only : U1,V1,U2,V2,tau,anudgeu,anudgev
      IMPLICIT NONE

      REAL*8 DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) :: UGCM, VGCM

c     LOCAL
      INTEGER i,j,l
      REAL  alphau,alphav,a,dtstep

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0S, J_1S, J_0SG, J_1SG

C****
C**** Extract useful local domain parameters from "grid"
C****
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      J_0 = grid%J_STRT
      J_1 = grid%J_STOP
      J_0S = grid%J_STRT_SKP
      J_1S = grid%J_STOP_SKP
      J_0SG = grid%J_STRT_STGR
      J_1SG = grid%J_STOP_STGR


c - - - - - - - - - - - - - - - - - - - - - - - - - |
c   x_gcm = a * x_gcm + ( 1 - a ) * x_reanalys      |
C----------------------------------------------------
             alphau=0.0001  !in sync_para in modelE
             alphav=0.0001
             alphau=dtstep * anudgeu !alphau !anudgeu
             alphav=dtstep * anudgev !alphav !nudgev
c             print*, ' A L P H A ', alphau, alphav,anudgeu,anudgev
         do l=1,lm
         do j= J_0SG, J_1SG
         do i=1,im
            a=(1.-tau)*u1(i,j,l)+tau*u2(i,j,l)         !time interpolation
c            ugcm(i,j,l)=alphau*ugcm(i,j,l)+(1-alphau)*a 
            ugcm(i,j,l) =(ugcm(i,j,l)+ a * alphau)/ (1+alphau) !nudging
            a=(1.-tau)*v1(i,j,l)+tau*v2(i,j,l)         !time interpolation
c            vgcm(i,j,l)=alphav*vgcm(i,j,l)+(1-alphav)*a 
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
      USE MODEL_COM, only : im,jm,lm,jhour,jday
      USE NUDGE_COM
      IMPLICIT NONE
      include 'netcdf.inc'
      INTEGER i
      integer start(4),count(4),status
C-----------------------------------------------------------------------
      if (jday.eq.1.and.jhour.eq.0) then
            zirk = zirk + 1
c -----------------------------------------------------------------
c   Opening of the files to be read
c -----------------------------------------------------------------
      print*, '  I N    N U D G E : OPEN NF FILES', jday,zirk
      if(zirk.eq.2) then
            status=NF_CLOSE('u.nc',NCNOWRIT,ncidu)
            status=NF_CLOSE('v.nc',NCNOWRIT,ncidv)

            status=NF_OPEN('u1.nc',NCNOWRIT,ncidu)
            status=NF_OPEN('v1.nc',NCNOWRIT,ncidv)
      endif     
      if(zirk.eq.3) then
            status=NF_CLOSE('u1.nc',NCNOWRIT,ncidu)
            status=NF_CLOSE('v1.nc',NCNOWRIT,ncidv)

            status=NF_OPEN('u2.nc',NCNOWRIT,ncidu)
            status=NF_OPEN('v2.nc',NCNOWRIT,ncidv)
      endif     
      if(zirk.eq.4) then
            status=NF_CLOSE('u2.nc',NCNOWRIT,ncidu)
            status=NF_CLOSE('v2.nc',NCNOWRIT,ncidv)

            status=NF_OPEN('u3.nc',NCNOWRIT,ncidu)
            status=NF_OPEN('v3.nc',NCNOWRIT,ncidv)
      endif     
      if(zirk.eq.5) then
            status=NF_CLOSE('u3.nc',NCNOWRIT,ncidu)
            status=NF_CLOSE('v3.nc',NCNOWRIT,ncidv)

            status=NF_OPEN('u4.nc',NCNOWRIT,ncidu)
            status=NF_OPEN('v4.nc',NCNOWRIT,ncidv)
      endif     
      if(zirk.eq.6) then
            status=NF_CLOSE('u4.nc',NCNOWRIT,ncidu)
            status=NF_CLOSE('v4.nc',NCNOWRIT,ncidv)

            status=NF_OPEN('u5.nc',NCNOWRIT,ncidu)
            status=NF_OPEN('v5.nc',NCNOWRIT,ncidv)
      endif     
      if(zirk.eq.7) then
            status=NF_CLOSE('u5.nc',NCNOWRIT,ncidu)
            status=NF_CLOSE('v5.nc',NCNOWRIT,ncidv)

            status=NF_OPEN('u6.nc',NCNOWRIT,ncidu)
            status=NF_OPEN('v6.nc',NCNOWRIT,ncidv)
      endif     
      if(zirk.eq.8) then
            status=NF_CLOSE('u6.nc',NCNOWRIT,ncidu)
            status=NF_CLOSE('v6.nc',NCNOWRIT,ncidv)

            status=NF_OPEN('u7.nc',NCNOWRIT,ncidu)
            status=NF_OPEN('v7.nc',NCNOWRIT,ncidv)
      endif     

      if(zirk.eq.9) then
            status=NF_CLOSE('u7.nc',NCNOWRIT,ncidu)
            status=NF_CLOSE('v7.nc',NCNOWRIT,ncidv)

            status=NF_OPEN('u8.nc',NCNOWRIT,ncidu)
            status=NF_OPEN('v8.nc',NCNOWRIT,ncidv)
      endif     

      if(zirk.eq.10) then
            status=NF_CLOSE('u8.nc',NCNOWRIT,ncidu)
            status=NF_CLOSE('v8.nc',NCNOWRIT,ncidv)

            status=NF_OPEN('u9.nc',NCNOWRIT,ncidu)
            status=NF_OPEN('v9.nc',NCNOWRIT,ncidv)
      endif     

          step_rea = 0
            status=NF_INQ_VARID(ncidu,'uwnd',uid)
            status=NF_INQ_VARID(ncidu,'level',plid)
            status=NF_INQ_VARID(ncidv,'vwnd',vid)
      endif

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
      USE DOMAIN_DECOMP, only : grid
      IMPLICIT NONE
      include 'netcdf.inc'

      integer timestep,l

      REAL*4 U(IM,JM-1,nlevnc),V(IM,JM-1,nlevnc),PL(nlevnc)
      REAL*4 DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) :: UN, VN

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
      USE DOMAIN_DECOMP, only : grid
       IMPLICIT NONE
c 
c ==============

       INTEGER lmo ! vertical dimensions of input
       INTEGER  i,j,lo,ln

        real po(lmo)! pressure levels in millibars of the input
        REAL varo(im,jm-1,lmo)!  Variable on the old grid (input)
        REAL varn(im,grid%J_STRT_HALO:grid%J_STOP_HALO,lm) !  Variable on the new grid (output)
        real coef

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0S, J_1S, J_0SG, J_1SG

C****
C**** Extract useful local domain parameters from "grid"
C****
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      J_0 = grid%J_STRT
      J_1 = grid%J_STOP
      J_0S = grid%J_STRT_SKP
      J_1S = grid%J_STOP_SKP
      J_0SG = grid%J_STRT_STGR
      J_1SG = grid%J_STOP_STGR


        do i=1,im
        do j= J_OSG, J_1SG
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
      print*, '  I N    N U D G E : OPEN NF FILES'
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
