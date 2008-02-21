c *****************************************************************
      MODULE NUDGE_COM
!@sum  NUDGE_COM contains all the nudging related variables
!@auth
!@ver
      USE MODEL_COM, only : im,jm,lm
!      USE DOMAIN_DECOMP, only : grid
      IMPLICIT NONE
      SAVE
!@param  nlevnc vertical levels of NCEP data
      INTEGER, PARAMETER :: nlevnc =17
!@var  U1, V1 NCEP wind at prior ncep timestep (m/s)
!@var  U2, V2 NCEP wind at the following ncep timestep (m/s)
      REAL*4, ALLOCATABLE, DIMENSION(:,:,:) :: u1,v1,u2,v2
      REAL*4, ALLOCATABLE, DIMENSION(:,:,:) :: UN1,VN1,UN2,VN2
      REAL*4, DIMENSION(nlevnc) :: pl
      REAL*8, DIMENSION(IM,JM,LM) :: u18,v18,u28,v28
      REAL*8, DIMENSION(IM,JM,LM) :: UN18,VN18,UN28,VN28
      REAL*8, DIMENSION(LM) :: pl8
!@var netcdf integer
      INTEGER :: ncidu,ncidv,uid,vid,plid
      INTEGER :: step_rea=1,zirk=0
!@var tau nudging time interpolation
      REAL*8 :: tau
!@param  anudgeu anudgev relaxation constant
      REAL*8 :: anudgeu = 0., anudgev = 0.

      END MODULE NUDGE_COM

      SUBROUTINE ALLOC_NUDGE(grid)

      USE DOMAIN_DECOMP, ONLY: DIST_GRID
      USE NUDGE_COM
!      USE MODEL_COM, only : im,jm,lm

      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H
      INTEGER :: IER

      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO
      ALLOCATE ( u1(IM,J_0H:J_1H,LM),
     $           v1(IM,J_0H:J_1H,LM),
     $           u2(IM,J_0H:J_1H,LM),
     $           v2(IM,J_0H:J_1H,LM),
     $   STAT = IER)
      ALLOCATE ( un1(IM,J_0H:J_1H,nlevnc),
     $           vn1(IM,J_0H:J_1H,nlevnc),
     $           un2(IM,J_0H:J_1H,nlevnc),
     $           vn2(IM,J_0H:J_1H,nlevnc),
     $   STAT = IER)

      RETURN
      END SUBROUTINE ALLOC_NUDGE

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
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     &     UGCM, VGCM
c     LOCAL
      INTEGER i,j,l
      REAL*8  alphau,alphav,a,dtstep

      INTEGER :: J_0SG, J_1SG

C****
C**** Extract useful local domain parameters from "grid"
C****
      J_0SG = grid%J_STRT_STGR
      J_1SG = grid%J_STOP_STGR
c - - - - - - - - - - - - - - - - - - - - - - - - - |
c   x_gcm = a * x_gcm + ( 1 - a ) * x_reanalys      |
C----------------------------------------------------
             alphau=dtstep * anudgeu  !alphau !anudgeu
             alphav=dtstep * anudgev  !alphav !nudgev

         do l=1,lm
           if (plbot(l+1).lt.pl(nlevnc)) cycle
         do j=J_0SG,J_1SG  ! Please pay attention j starts at 2
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
      USE DOMAIN_DECOMP, only: grid,am_i_root
      USE MODEL_COM, only: im,jm,lm,jhour,jday,itime,nday,jyear,iyear1
      USE NUDGE_COM
      IMPLICIT NONE
      include 'netcdf.inc'
      INTEGER i
      INTEGER J_0,J_1S
      integer start(4),count(4),status
      character(len=3) :: nstr1,nstr2
C-----------------------------------------------------------------------
      J_0=grid%J_STRT
      J_1S=grid%J_STOP_SKP

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

        if(am_i_root()) then
            status=NF_CLOSE('u'//trim(nstr1)//'.nc',NCNOWRIT,ncidu)
            status=NF_CLOSE('v'//trim(nstr1)//'.nc',NCNOWRIT,ncidv)

            status=NF_OPEN('u'//trim(nstr2)//'.nc',NCNOWRIT,ncidu)
            status=NF_OPEN('v'//trim(nstr2)//'.nc',NCNOWRIT,ncidv)

            status=NF_INQ_VARID(ncidu,'level',plid)
            status=NF_INQ_VARID(ncidu,'uwnd',uid)
            status=NF_INQ_VARID(ncidv,'vwnd',vid)
          endif 

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
      USE DOMAIN_DECOMP, only : grid, unpack_data, am_i_root, esmf_bcast
      IMPLICIT NONE
      include 'netcdf.inc'

      REAL*4, DIMENSION(IM,JM,nlevnc) :: scratch_glob
      real*8 scratch(IM,grid%j_strt_halo:grid%j_stop_halo,nlevnc)
      integer timestep,l

      integer start(4),count(4),status
      integer J_0,J_1S

      J_0=grid%J_STRT
      J_1S=grid%J_STOP_SKP

c -----------------------------------------------------------------
c   read u, v, and pl
c -----------------------------------------------------------------
      if(am_i_root()) then
         status=NF_GET_VARA_REAL(ncidu,plid,1,nlevnc,pl)
         pl8=0.d0
         pl8(1:nlevnc)=dble(pl(1:nlevnc))
      endif
      call esmf_bcast(grid,pl8)
      pl(1:nlevnc)=sngl(pl8(1:nlevnc))
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=timestep

      count(1)=im
      count(2)=jm-1
      count(3)=nlevnc
      count(4)=1

c  u
      scratch_glob=0.
      if(am_i_root()) then
         status=NF_GET_VARA_REAL(ncidu,uid,start,count,
     &                        scratch_glob(:,1:JM-1,:))
      endif
      call unpack_data(grid, dble(scratch_glob), scratch)
      un2=sngl(scratch)
c  v

      if(am_i_root()) then
         status=NF_GET_VARA_REAL(ncidv,vid,start,count,
     &                        scratch_glob(:,1:JM-1,:))
      endif
      call unpack_data(grid, dble(scratch_glob), scratch)
      vn2=sngl(scratch)

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
      USE DOMAIN_DECOMP, only : grid, HALO_UPDATE, SOUTH
       IMPLICIT NONE
c
c ==============

       INTEGER lmo ! vertical dimensions of input
       INTEGER  i,j,lo,ln

        real po(lmo)! pressure levels in millibars of the input
        REAL*8 scratch(im,grid%J_STRT_HALO:grid%J_STOP_HALO,lmo)
        REAL varo(im,grid%J_STRT_HALO:grid%J_STOP_HALO,lmo)!  Variable on the old grid (input)
        REAL varn(im,grid%J_STRT_HALO:grid%J_STOP_HALO,lm) !  Variable on the new grid (output)
        real coef,dp1,dp2
        integer J_0SG,J_1SG

        J_0SG = grid%J_STRT_STGR
        J_1SG = grid%J_STOP_STGR
        scratch = dble(varo)
        CALL HALO_UPDATE(grid,scratch,FROM=SOUTH)
        varo = sngl(scratch)
        do j= J_0SG, J_1SG             ! Please pay attention j starts at 2
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
      USE DOMAIN_DECOMP, only: grid, unpack_data, am_i_root, esmf_bcast
      USE MODEL_COM, only : im,jm,lm,jhour,jday,itime,nday,jyear,iyear1
     &                      ,itime,itimei
      USE NUDGE_COM
      USE PARAM
      IMPLICIT NONE
      include 'netcdf.inc'

      REAL*4, DIMENSION(IM,JM,nlevnc) :: scratch_glob
      real*8 scratch(IM,grid%j_strt_halo:grid%j_stop_halo,nlevnc)
      integer start(4),count(4),status
      character(len=3) :: nstr1

C**** Rundeck parameters:
      call sync_param( "ANUDGEU", ANUDGEU )
      call sync_param( "ANUDGEV", ANUDGEV )

C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      IF(Itime == ItimeI) THEN  
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

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

C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      ENDIF                      
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

c -----------------------------------------------------------------
c   Opening of the files to be read
c -----------------------------------------------------------------

        if (zirk.lt.10) then
           write(nstr1,'(I1)') zirk
        elseif (zirk.lt.100) then
           write(nstr1,'(I2)') zirk
        endif

      if(am_i_root()) then
            status=NF_OPEN('u'//trim(nstr1)//'.nc',NCNOWRIT,ncidu)
            status=NF_OPEN('v'//trim(nstr1)//'.nc',NCNOWRIT,ncidv)

            status=NF_INQ_VARID(ncidu,'level',plid)
            status=NF_INQ_VARID(ncidu,'uwnd',uid)
            status=NF_INQ_VARID(ncidv,'vwnd',vid)
      endif

C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      IF(Itime == ItimeI) THEN   
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

C------------------------------------------------------------------
c -----------------------------------------------------------------
c   read u, v, and pl
c -----------------------------------------------------------------
      if(am_i_root()) then
         status=NF_GET_VARA_REAL(ncidu,plid,1,nlevnc,pl)
         pl8=0.d0
         pl8(1:nlevnc)=dble(pl(1:nlevnc))
      endif
      call esmf_bcast(grid,pl8)
      pl(1:nlevnc)=sngl(pl8(1:nlevnc))

      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=step_rea

      count(1)=im
      count(2)=jm-1
      count(3)=nlevnc
      count(4)=1

      if(am_i_root()) then
      scratch_glob=0.
         status=NF_GET_VARA_REAL(ncidu,uid,start,count,
     &                        scratch_glob(:,1:JM-1,:))
      endif
      call unpack_data(grid, dble(scratch_glob), scratch)
      un1=sngl(scratch)
      if(am_i_root()) then
         status=NF_GET_VARA_REAL(ncidv,vid,start,count,
     &                        scratch_glob(:,1:JM-1,:))
      endif
      call unpack_data(grid, dble(scratch_glob), scratch)
      vn1=sngl(scratch)

      start(4)=step_rea+1
      if(am_i_root()) then
         status=NF_GET_VARA_REAL(ncidu,uid,start,count,
     &                        scratch_glob(:,1:JM-1,:))
      endif
      call unpack_data(grid, dble(scratch_glob), scratch)
      un2=sngl(scratch)
      if(am_i_root()) then
         status=NF_GET_VARA_REAL(ncidv,vid,start,count,
     &                        scratch_glob(:,1:JM-1,:))
      endif
      call unpack_data(grid, dble(scratch_glob), scratch)
      vn2=sngl(scratch)

C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      ENDIF 
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      return
      end subroutine nudge_init


c -----------------------------------------------------------------


      subroutine io_nudge(kunit,iaction,ioerr)
!@sum  io_nudge reads and writes met-nuding arrays to file
!@auth G. Faluvegi
!@ver  1.0
! note that the nudging was never updated for MPI/ESMF, so for the
! moment neither is this routine. This was just put in to try and
! fix an irreproducibility issue when nudging is on.
      use model_com, only : ioread,iowrite,irerun
      use model_com, only : im,jm
      use domain_decomp, only : grid, am_i_root
      use domain_decomp, only : pack_data, unpack_data, esmf_bcast
      use nudge_com
      implicit none
      
      real*8 scratch(IM,grid%j_strt_halo:grid%j_stop_halo,LM)
      integer kunit   !@var kunit unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
!@var ioerr 1 (or -1) if there is (or is not) an error in i/o
      integer, intent(inout) :: ioerr
C      INTEGER :: J_1H, J_0H

C      J_0H = grid%J_STRT_HALO
C      J_1H = grid%J_STOP_HALO

      select case (iaction)

      case (:iowrite)            ! output to standard restart file
      
        call pack_data(grid,dble(u1),u18)
        call pack_data(grid,dble(u2),u28)
        call pack_data(grid,dble(v1),v18)
        call pack_data(grid,dble(v2),v28)
        un18(:,:,:)=0.d0
        call pack_data(grid,dble(un1),un18(:,:,1:nlevnc))
        un28(:,:,:)=0.d0
        call pack_data(grid,dble(un2),un28(:,:,1:nlevnc))
        vn18(:,:,:)=0.d0
        call pack_data(grid,dble(vn1),vn18(:,:,1:nlevnc))
        vn28(:,:,:)=0.d0
        call pack_data(grid,dble(vn2),vn28(:,:,1:nlevnc))
        pl8(:)=0.d0
        pl8(1:nlevnc)=dble(pl(1:nlevnc))

        if(am_i_root()) write (kunit,err=10) u18,u28,v18,v28,un18,un28,
     &      vn18,vn28,pl8,zirk,step_rea,ncidu,ncidv,plid,uid,vid,tau

      case (ioread,irerun)      ! input from restart file (not for irsfic)

        if(am_i_root()) read (kunit,err=10)  u18,u28,v18,v28,un18,un28,
     &      vn18,vn28,pl8,zirk,step_rea,ncidu,ncidv,plid,uid,vid,tau

        call unpack_data(grid, u18, scratch)
        u1(:,:,:)=sngl(scratch)
        call unpack_data(grid, u28, scratch)
        u2(:,:,:)=sngl(scratch)
        call unpack_data(grid, v18, scratch)
        v1(:,:,:)=sngl(scratch)
        call unpack_data(grid, v28, scratch)
        v2(:,:,:)=sngl(scratch)
        call unpack_data(grid, un18, scratch)
        un1(:,:,:)=sngl(scratch(:,:,1:nlevnc))
        call unpack_data(grid, un28, scratch)
        un2(:,:,:)=sngl(scratch(:,:,1:nlevnc))
        call unpack_data(grid, vn18, scratch)
        vn1(:,:,:)=sngl(scratch(:,:,1:nlevnc))
        call unpack_data(grid, vn28, scratch)
        vn2(:,:,:)=sngl(scratch(:,:,1:nlevnc))
        call esmf_bcast(grid,pl8)
        pl(1:nlevnc)=sngl(pl8(1:nlevnc))
        call esmf_bcast(grid,zirk)
        call esmf_bcast(grid,step_rea)
        call esmf_bcast(grid,ncidu)
        call esmf_bcast(grid,ncidv)
        call esmf_bcast(grid,plid)
        call esmf_bcast(grid,uid)
        call esmf_bcast(grid,vid)
        call esmf_bcast(grid,tau)

      end select

      return
 10   ioerr=1

      return
      end subroutine io_nudge
