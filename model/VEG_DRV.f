#include "rundeck_opts.h"

      module veg_drv
!@sum veg_drv contains variables and routines for vegetation driver
!@auth I. Alienov, N. Kiang

      use model_com, only : im,jm
      implicit none
      private
      save

      public init_vegetation,reset_veg_to_defaults
     &     ,veg_set_cell, veg_save_cell,upd_gh

      real*8,public :: cosday,sinday

      contains



      subroutine init_vegetation(redogh,istart)
!@sum initializes vegetation
      integer, intent(in) :: istart
      logical, intent(in) :: redogh

      call read_veg_data(redogh,istart)
      call upd_gh
      end subroutine init_vegetation


      subroutine read_veg_data(redogh,istart)
!@sum reads vegetation arrays and rundeck parameters
      use filemanager
      use param
      use DOMAIN_DECOMP, only : GRID, GET, READT_PARALLEL
      use vegetation, only : cond_scheme,vegCO2X_off,crops_yr
      use veg_com
      use model_com, only : fearth,jyear

      implicit none

      integer, intent(in) :: istart
      logical, intent(in) :: redogh

      INTEGER :: J_1, J_0, J_1H, J_0H
      INTEGER :: J_0S, J_1S, J_0STG, J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      integer iu_veg
      integer i, j, k
      logical veg_data_missing

      real*8, allocatable :: veg_c4(:,:)
      real*8 :: vc4
      integer :: read_c4_grass = 0

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

c**** read rundeck parameters
      call sync_param( "cond_scheme", cond_scheme)  !nyk 5/1/03
      call sync_param( "vegCO2X_off", vegCO2X_off)  !nyk 3/2/04
      call sync_param( "crops_yr", crops_yr)
      call sync_param( "read_c4_grass", read_c4_grass)

c**** read land surface parameters or use defaults
      call openunit("VEG",iu_VEG,.true.,.true.)
      do k=1,10                 !  11 ????
        CALL READT_PARALLEL
     *    (grid,iu_VEG,NAMEUNIT(iu_VEG),0,vdata(:,:,K),1)
      end do
c**** zero-out vdata(11) until it is properly read in
      vdata(:,:,11) = 0.
      vdata(:,:,12) = 0.
      call closeunit(iu_VEG)
c**** add data on c4 grass
      if ( read_c4_grass >= 1 ) then
        print *,"Adding c4 grass data"
        call openunit("VEG_C4",iu_VEG,.true.,.true.)
        allocate( veg_c4(im,J_0H:J_1H) )
        CALL READT_PARALLEL
     *       (grid,iu_VEG,NAMEUNIT(iu_VEG),0,veg_c4(:,:),1)
        do j=J_0,J_1
          do i=1,im
            if (fearth(i,j).gt.0) then
              ! normalize to earth fraction
              vc4 = veg_c4(i,j) / fearth(i,j)
              ! add c4 grass up to max amount of current grass
              vdata(i,j,12) = min( vdata(i,j,3), vc4 )
              vdata(i,j,3) = vdata(i,j,3) - vdata(i,j,12)
            end if
          enddo
        enddo      
        deallocate( veg_c4 )
        call closeunit(iu_VEG)
      endif
C**** Update vegetation file if necessary (i.e. crops_yr =0 or >0)
      if(crops_yr.eq.0) call updveg(jyear,.false.)
      if(crops_yr.gt.0) call updveg(crops_yr,.false.)

      if (istart.le.2 .or. redogh) then ! initial. foliage arrays (adf)
        Cint(:,:)=0.0127D0      ! internal CO2
        Qfol(:,:)=3.D-6         ! surface mixing ratio
        cnc_ij(:,:) = 0.d0
      end if
      if (istart.le.0) return   ! avoid reading unneeded files

     !!! spgsn=.1d0

c**** check whether ground hydrology data exist at this point.
      veg_data_missing = .false.
      do j=J_0,J_1
        do i=1,im
          if (fearth(i,j).gt.0) then
            if ( sum(vdata(i,j,1:12)).eq.0 ) then
              print *,"No vegetation data: i,j=",i,j,vdata(i,j,1:12)
              veg_data_missing = .true.
            end if
          end if
        enddo
      enddo
      if ( veg_data_missing ) then
        write(6,*) 'Vegetation data is missing at some pts'
        write(6,*) 'If you have a non-standard land mask, please'
        write(6,*) 'consider using extended GH data and rfs file.'
        call stop_model(
     &       'Vegetation data is missing at some cells',255)
      endif

      end subroutine read_veg_data


      subroutine upd_gh
!@sum initializes (or re-initializes) the vegetation data
      use constant, only : twopi,one
      use DOMAIN_DECOMP, only : GRID, GET, READT_PARALLEL
      use model_com, only : fearth,jeq
      use veg_com !, only : vdata,Cint,Qfol
      use ghy_com, only : ngm
      use ghy_com, only : dz_ij
!      use surf_albedo, only: albvnh  !nyk !!! not used? i.a.

      implicit none

!---  local vars
      integer i, j
      real*8 dif,frdn,frup,pearth,phase,scs0,scsim,scsre,sfv,sla0
      real*8 almass0, almassre, almassim  !nyk
!      real*8 sla0f, slimf, slref !nyk
      real*8 slim,slre,svh,z
      real*8 snm,snf  ! temporary sums (adf)
      real*8 cwc_sum
      integer iv, l
      real*8 fv,dz(ngm)
      integer n
!---  parameters for different types of vegetation
!---          tundr  grass  shrub  trees  decid evrgr  rainf crops
      real*8, parameter :: alamax(11) =
     $     (/ 1.5d0, 2.0d0, 2.5d0, 4.0d0, 6.0d0,10.0d0,8.0d0,4.5d0
     &     ,0.d0, 0.d0, 2.d0 /)
      real*8, parameter :: alamin(11) =
     $     (/ 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 8.0d0,6.0d0,1.0d0
     &     ,0.d0, 0.d0, 1.d0 /)

      real*8, parameter :: aroot(11) =
     $     (/ 12.5d0, 0.9d0, 0.8d0,0.25d0,0.25d0,0.25d0,1.1d0,0.9d0
     &     ,0.d0, 0.d0, 0.9d0 /)
      real*8, parameter :: broot(11) =
     $     (/  1.0d0, 0.9d0, 0.4d0,2.00d0,2.00d0,2.00d0,0.4d0,0.9d0
     &     ,0.d0, 0.d0, 0.9d0 /)
      real*8, parameter :: rsar(11) =
     $     (/100d0, 100d0, 200d0, 200d0, 200d0, 300d0,250d0, 125d0
     &     ,0.d0, 0.d0, 100d0 /)
      real*8, parameter :: vhght(11) =
     $     (/0.1d0, 1.5d0,   5d0,  15d0,  20d0,  30d0, 25d0,1.75d0
     &     ,0.d0, 0.d0, 1.5d0 /)
! Mean canopy nitrogen (nmv; g/m2[leaf]) and Rubisco factors (nfv) for each
! vegetation type (adf)
      real*8, parameter :: nmv(11) =
     $     (/1.6d0,0.82d0,2.38d0,1.03d0,1.25d0,2.9d0,2.7d0,2.50d0
     &     ,0.d0, 0.d0, 0.82d0 /)
      real*8, parameter :: nfv(11) =
     $     (/1.4d0,1.5d0 ,1.3d0 ,1.3d0 ,1.5d0 ,0.9d0,1.1d0,1.3d0
     &     ,0.d0, 0.d0, 1.5d0 /)
      integer, parameter :: laday(11) =
     $     (/ 196,  196,  196,  196,  196,  196,  196,  196
     &     ,0, 0, 196 /)
      real*8, parameter :: can_w_coef(11) =
     &     (/ 1.d-4, 1.d-4, 1.d-4, 1.d-4, 1.d-4, 1.d-4, 1.d-4, 1.d-4
     &     ,0.d0, 0.d0, 1.d-4 /)

! Specific leaf areas (sleafa, kg[C]/m2) (adf, nyk)
! Values below 1/(m2/kg) to get kg/m2 for multiplying.
! Sources: White, M.A., et.al. (2000), Earth Interactions, 4:1-85.
!          Leonardos,E.D.,et.al.(2003), Physiologia Plantarum, 117:521+.
!               From winter wheat grown at 20 C (20 m2/kg[dry mass])
!                                      and  5 C (13 m2/kg[dry mass])
!          Francesco Tubiello, personal communication, crop 18-20 m2/kg.
      real*8, parameter :: sleafa(11) =
     $     (/ 1./30.5d0,1./49.0d0,1./30.5d0,1./40.5d0,1./32.0d0,
     $        1./8.2d0,1./32.0d0,1./18.0d0, 0.d0, 0.d0, 1./49.0d0 /) 

c****             tundr grass shrub trees decid evrgr rainf crops
c****
c**** laday(veg type, lat belt) = day of peak lai
c old peak lai:  2nd line is for latitudes < 23.5 deg
c****    1  temperate latitudes
c****    2  non-temperate latitudes
c     data  laday/ 196,  196,  196,  196,  196,  196,  105,  196/
c     data  laday/ 196,  288,  288,  288,  288,  196,  105,  288/
c****
c**** contents of ala(k,i,j),  lai coefficients
c****   1  average leaf area index
c****   2  real amplitude of leaf area index
c****   3  imaginary amplitude of leaf area index
c****
c**** contents of almass(i,j),  leaf mass
c****      leaf mass = ala*sleafa = kg[C]/ground area
c****
c**** contents of acs(k,i,j),  cs coefficients
c****   1  average stomatal conductance
c****   2  real amplitude of stomatal conductance
c****   3  imaginary amplitude of stomatal conductance
c****
!@dbparam ghy_default_data if == 1 reset all GHY data to defaults
!@+ (do not read it from files)
      integer :: ghy_default_data = 0
      integer :: northsouth  !1=south, 2=north hemisphere

      INTEGER :: I_0, I_1, J_1, J_0, J_1H, J_0H
      INTEGER :: J_0S, J_1S, J_0STG, J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

!      entry upd_gh ! need to redo if vdata changes
C****
C**** Extract parameters from "grid" in case we entered here
C****
      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
c****
c**** set the global arrays  ala, acs, afb, afr, anm, anf
c****
      ala(:,:,:)=0.
      alaf(:,:,:,:)=0.          !nyk lai by veg functional type
      alaif(:,:,:)=0.
      acs(:,:,:)=0.
      avh(:,:)=0.
      afb(:,:)=0.
      afr(:,:,:)=0.
      acs(1,:,:)=.01d0
      anm(:,:)=0.0D0            ! Global mean canopy N array (adf)
      anf(:,:)=0.0D0            ! Global Rubisco factors (adf)
      almass(:,:,:)=0.0D0       ! Global leaf mass at a time, nyk
!rar  aalbveg(:,:)=0.08D0 ! no need, it is set in daily_earth
      can_w_capacity(:,:) = 0.d0

      do j=J_0,J_1
        if(j.le.jm/2) then
          northsouth=1.d0  !southern hemisphere
        else
          northsouth=2.d0  !northern hemisphere
        end if
        do i=1,im
          pearth=fearth(i,j)
          afb(i,j)=vdata(i,j,1)+vdata(i,j,10)
          if(afb(i,j).gt..999) afb(i,j)=1.
          if(pearth.le.0..or.afb(i,j).ge.1.) cycle
c**** calculate lai, cs coefficicents
          sfv=0.d0
          sla0=0.     !For calculating sla
          slre=0.     !For calculating sla
          slim=0.     !For calculating sla
          !sla0f=0.     !For calculating alaf
          !slref=0.     !For calculating alaf
          !slimf=0.     !For calculating alaf
          almass0=0.  !nyk For specific leaf area
          almassre=0. !nyk For specific leaf area
          almassim=0. !nyk For specific leaf area
          scs0=0.
          scsre=0.
          scsim=0.
          svh=0.
          snm=0. ! adf
          snf=0. ! adf
          cwc_sum = 0.d0
          do iv=1,11
            if ( iv==9 .or. iv==10 ) cycle
            phase=twopi*laday(iv)/365.
            if(j.lt.jeq) phase=phase+twopi/2.
            fv=vdata(i,j,iv+1)
            sfv=sfv+fv
            svh=svh+fv*vhght(iv)
            snm=snm+fv*nmv(iv) ! adf
            snf=snf+fv*nfv(iv) ! adf
            dif=(alamax(iv) - alamin(iv))
            sla0=sla0+fv*(alamax(iv) + alamin(iv))
            slre=slre+fv*dif*cos(phase)
            slim=slim+fv*dif*sin(phase)
            !nyk-------------
            !alaf and alaif
            alaf(1,iv,i,j) = 0.5*(alamax(iv) + alamin(iv))
            alaf(2,iv,i,j) = 0.5*dif*cos(phase)
            alaf(3,iv,i,j) = 0.5*dif*sin(phase)
            alaif(iv,i,j)= alaf(1,iv,i,j)+
     $           cosday*alaf(2,iv,i,j)+sinday*alaf(3,iv,i,j) !save ij
           !nyk-------------
            !almaxmin = almaxmin + sleafa(iv)*fv*(alamax(iv)-alamin(iv))
            almass0 = almass0 + sleafa(iv)*fv*(alamax(iv) - alamin(iv))
            !almass0 = almass0 + sleafa(iv)*fv*(alamax(iv) + alamin(iv))
            !almassre = almassre + sleafa(iv)*fv*dif*cos(phase)
            !almassim = almassim + sleafa(iv)*fv*dif*sin(phase)
            !----------------
            scs0=scs0+fv*(alamax(iv) + alamin(iv))/rsar(iv)
            scsre=scsre+fv*dif*cos(phase)/rsar(iv)
            scsim=scsim+fv*dif*sin(phase)/rsar(iv)
            cwc_sum = cwc_sum + fv*can_w_coef(iv)
          end do
          ala(1,i,j)=.5/sfv*sla0
          ala(2,i,j)=.5/sfv*slre
          ala(3,i,j)=.5/sfv*slim

          acs(1,i,j)=.5/sfv*scs0
          acs(2,i,j)=.5/sfv*scsre
          acs(3,i,j)=.5/sfv*scsim
          avh(i,j)=svh/sfv
          can_w_capacity(i,j) = cwc_sum/sfv
          anm(i,j)=snm/sfv ! adf
          anf(i,j)=snf/sfv ! adf
          almass(1,i,j) = almass0   !This just computes total growth for
          almass(2,i,j) = 0.   ! year via difference between max and min
          almass(3,i,j) = 0.
          !almass(1,i,j)= 0.5/sfv*almass0 !nyk
          !almass(2,i,j)= 0.5/sfv*almassre !nyk
          !almass(3,i,j)= 0.5/sfv*almassim !nyk

c**** calculate root fraction afr averaged over vegetation types
          do n=1,ngm
            dz(n)=dz_ij(i,j,n)
            if(dz(n).le.0.) go to 320
          end do
 320      n=n-1
          do iv=1,11
            if ( iv==9 .or. iv==10 ) cycle
            fv=vdata(i,j,iv+1)
            z=0.
            frup=0.
            do l=1,n
              z=z+dz(l)
              frdn=aroot(iv)*z**broot(iv)
              frdn=min(frdn,one)
              if(l.eq.n)frdn=1.
              afr(l,i,j) = afr(l,i,j) + fv*(frdn-frup)
              frup=frdn
            end do
          end do
          do l=1,n
            afr(l,i,j) = afr(l,i,j)/(1.-afb(i,j))
          end do
        end do
      end do

      ! write(98,*) afr

      return
      end subroutine upd_gh


      subroutine reset_veg_to_defaults( reset_prognostic )
      use veg_com, only: vdata
      use DOMAIN_DECOMP, only : GRID, GET
      logical, intent(in) :: reset_prognostic
      integer i,j

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0S, J_1S, J_0STG, J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      do j=J_0,J_1
      do i=1,im

      vdata(i,j,1:12)= (/  0.00000000d+00,  0.00000000d+00,
     &     0.00000000d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.62451148d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.37548852d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.00000000d+00 /)
      enddo
      enddo

      end subroutine reset_veg_to_defaults


      subroutine veg_set_cell (i0,j0,ghy_data_only)
!@sum resets the vegetation module to a new cell i0,j0
      use constant, only : gasc
      use ghy_com, only : ngm
      use sle001, only : fr,snowm,ws,shc,shw,pres,ts
      use vegetation, only : alaie,rs,nm,nf,alai,vh
     &     ,alait,vfraction
     &     ,fdir,parinc,vegalbedo,sbeta,Ci,Qf ! added by adf, nyk
     &     ,Ca     ! nyk
     &     ,cond_scheme,vegCO2X_off  !nyk
     &     ,cnc
      use rad_com, only : cosz1, CO2X, CO2ppm  !nyk
     &    ,FSRDIR,SRVISSURF  !adf, nyk
      use veg_com, only :
     &     Cint,Qfol           ! added by adf
     &     ,cnc_ij
     &     ,aalbveg    ! nyk
     &     ,afr,avh,anm,anf,ala,alaif,alaf,vdata,acs


      implicit none

      integer, intent(in) :: i0,j0
!@var ghy_data_only set only ghy data (no call to vegetatin expected)
      logical, optional :: ghy_data_only
      real*8, parameter :: spgsn=.1d0
      integer l
      real*8 aa
!      real*8 aalbveg0, sfv  !nyk
!      integer northsouth,iv  !nyk
!      real*8 alaic,vh,shtpr,alai  ! adf
      real*8 alaic

      ! call stop_model("veg_set_cell called", 255)

c**** fr: root fraction in layer l  (1=fr(1)+fr(2)+...+fr(n))
      do l=1,ngm
        fr(l)=afr(l,i0,j0)
      end do
c**** vh: vegetation height
      vh=avh(i0,j0)
      nm=anm(i0,j0) ! mean canopy nitrogen (g/m2[leaf]) (adf)
      nf=anf(i0,j0) ! canopy nitrogen factor (adf)
      snowm=vh*spgsn

c**** alai: leaf area index
      alai=ala(1,i0,j0)+cosday*ala(2,i0,j0)+sinday*ala(3,i0,j0)
      alai=max(alai,1.d0)
      !lai by functional type, not normalized by cover fraction!
      alaif(:,i0,j0) = 0.d0
      alait(:) = 0.d0
      vfraction(:) = 0.d0
      do L=1,11
        if ( L==9 .or. L==10 ) cycle
        alaif(L,i0,j0)=alaf(1,L,i0,j0)+
     &       cosday*alaf(2,L,i0,j0)+sinday*alaf(3,L,i0,j0) !save ij
        alait(L) = alaif(L,i0,j0)  !for VEGETATION.f
        vfraction(L)=vdata(i0,j0,L+1) !for VEGETATION.f
      end do

      alaic=5.0
      alaie=alaic*(1.-exp(-alai/alaic))
c**** rs: minimum stomatal resistance
      rs=alai/(acs(1,i0,j0)+cosday*acs(2,i0,j0)+sinday*acs(3,i0,j0))
c???  cnc=alai/rs   redefined before being used (qsbal,cond)

      !---------------------------------------------------------
      !nyk vegetation albedo.  Only really updated daily, but have to
      !get it initialized somewhere after ALBVNH is calculated.
      !ALBVNH is unfortunately set *after* ground hydr is initialized.
      !albvnh(9,6,2)=albvnh(1+8veg,6bands,2hemi), band 1 is VIS.
!      aalbveg0=0.d0               !nyk
!      sfv=0.d0
!      if (j0.le.jm/2) then
!
!       northsouth=1.d0         !southern hemisphere
!      else
!        northsouth=2.d0         !northern hemisphere
!      end if
!      do iv=1,8
!        aalbveg0 = aalbveg0 + fv*(ALBVNH(iv+1,1,northsouth))
!        sfv = sfv + fv
!      end do
!      aalbveg(i0,j0) = aalbveg0/sfv
!      write (99,*) 'aalbveg ghinij', aalbveg(i0,j0) !nyk
      !---------------------------------------------------------

      ws(0,2)=.0001d0*alai
!!!      ws(0,2)=can_w_capacity(i0,j0)*alai
c shc(0,2) is the heat capacity of the canopy
      aa=ala(1,i0,j0)
      shc(0,2)=(.010d0+.002d0*aa+.001d0*aa**2)*shw

! This is a hack to prevent the program from using uninitialized
! radiation arrays at the beginning of the run.
! If program returns at this point, it sets only the data needed for
! ground hydrology but not the vegetation which is ok if veg_conductance
! is not called.
      if ( present(ghy_data_only) ) then
        if ( ghy_data_only ) return
      endif

!----------------------------------------------------------------------!
      if (cond_scheme.eq.2) then  !new conductance scheme
! Sine of solar elevation (rad).
        sbeta=cosz1(i0,j0)
! adf Fraction of solar radiation at ground that is direct beam.
        fdir=FSRDIR(i0,j0)
! nyk Calculate incident PAR, photosynthetically active radiation.
!     *SRVISSURF is from SRDVIS:
!     = incident visible solar radiation (dir+dif) on the surface
!     = estimated as 53% of total solar flux density, so wavelength
!       range is UV through ~760 or 770 nm cutoff (not strict)
!     *SRDVIS is normalized to solar zenith = 0.
!     *From integrating solar flux density (TOA) over solar spectrum,
!       PAR(400-700 nm) is ~ 82% of the flux density of SRDVIS.
        parinc=0.82*SRVISSURF(i0,j0)*sbeta
! nyk Get vegetation grid albedo, temporary until canopy scheme in place
        vegalbedo = aalbveg(i0,j0)
! Internal foliage CO2 concentration (mol/m3).
        Ci=Cint(i0,j0)
! Foliage surface mixing ratio (kg/kg).
        Qf=Qfol(i0,j0)
        CNC = cnc_ij(i0,j0)
! Atmospheric CO2 concentration at land surface (mol/m3)
        if (vegCO2X_off==0) then  !Default
           Ca = CO2X*CO2ppm*(1.0D-06)*pres*100.0/gasc/ts
        else
           Ca = CO2ppm*(1.0D-06)*pres*100.0/gasc/ts
        endif
!        write(99,*) "vegCO2X_off,CO2X,CO2ppm,Ca,pres,gasc,ts",
!     %       vegCO2X_off,CO2X,CO2ppm,Ca,pres,gasc,ts
!        stop
      end if
!----------------------------------------------------------------------!
      return
      end subroutine veg_set_cell


      subroutine veg_save_cell(i,j)
      !Save vegetation arrays for cell.
      use vegetation, only: cond_scheme, Ci, Qf, CNC
      use veg_com, only:  Cint, Qfol, cnc_ij

      integer, intent(in):: i,j

      if (cond_scheme.eq.2) then      !new conductance scheme only
        Cint(i,j)=Ci                  ! Save last value
        Qfol(i,j)=Qf                  ! Save last value
        cnc_ij(i,j)=CNC
      end if
      end subroutine veg_save_cell



      end module veg_drv

!***********************************************************************

      subroutine updveg (year,reset_veg)
!@sum  reads appropriate crops data and updates the vegetation file
!@auth R. Ruedy
!@ver  1.0
      USE FILEMANAGER
      USE MODEL_COM, only : im,jm
      USE DOMAIN_DECOMP, only : GRID, GET
      use veg_com, only : vdata
      USE GEOM, only : imaxj
      use veg_drv, only : upd_gh
      implicit none
      integer, intent(in) :: year
      logical, intent(in) :: reset_veg

      real*8 wt,crops         ! temporary vars ! ,crop(im,jm)

      integer :: year1,year2,year_old=-1, iu, i,j,k
      real*8, Allocatable :: crop1(:,:), crop2(:,:), vdata0(:,:,:)
      save   year1,year2,year_old,vdata0,crop1,crop2,iu      ! to limit i/o

      character*80 title
      real*4 crop4(im,jm)

      INTEGER :: I_0, I_1, J_1, J_0, J_0H, J_1H
      INTEGER :: J_0S, J_1S, J_0STG, J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      CALL GET(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               J_STRT_HALO =J_0H,   J_STOP_HALO =J_1H,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

C**** check whether update is needed
      if (year.eq.year_old) return

C**** first iteration actions:
      if (year_old.lt.0) then
C****     check whether a no-crops vege-file was used
        do j=J_0S,J_1S
        do i=1,im
           if(vdata(i,j,9).gt.0.)
     *     call stop_model('updveg: use no_crops_VEG_file',255)
        end do
        end do
C****     open and read input file
        call openunit('CROPS',iu,.true.,.true.)
        read(iu) title,crop4

        read(title,*) year1
        Allocate(crop1(IM,J_0H:J_1H), crop2(IM,J_0H:J_1H),
     &           VDATA0(IM,J_0H:J_1H,12))
        crop1(:,J_0:J_1)=crop4(:,J_0:J_1) 
        crop2(:,J_0:J_1)=crop4(:,J_0:J_1)
        year2=year1
        if (year1.ge.year)          year2=year+1
C****     save orig. (no-crop) vdata to preserve restart-independence
        vdata0 = vdata
      end if

      wt=0.
      do while (year2.lt.year)
         year1 = year2 ; crop1 = crop2
         read (iu,end=10) title,crop4
         read(title,*) year2
         crop2(:,J_0:J_1) = crop4(:,J_0:J_1)
      end do
      wt = (year-year1)/(real(year2-year1,kind=8))
   10 continue
      write(6,*) 'Using crops data from year',year1+wt*(year2-year1)

C**** Modify the vegetation fractions
      do j=J_0,J_1
      do i=1,imaxj(j)
         if (crop1(i,j).ge.0.) then
            crops = crop1(i,j) + wt*(crop2(i,j)-crop1(i,j))
            do k=1,12
              vdata(i,j,k) = vdata0(i,j,k)*(1.-crops)
            end do
              vdata(i,j,9) = crops
         end if
      end do
      end do
      if(reset_veg) call upd_gh

      year_old = year

      return
      end subroutine updveg



