#include "rundeck_opts.h"
      MODULE ATM_COM
!@sum  ATM_COM Main atmospheric variables
!@auth Original Development Team
      USE RESOLUTION, only : LM
#ifdef USE_ESMF
      USE ESMF_MOD, only: ESMF_Clock
#endif
      IMPLICIT NONE
      SAVE

#ifdef USE_ESMF
! different components can define their own clocks later
      Type (ESMF_CLOCK) :: atmclock
#endif

!@var LM_REQ Extra number of radiative equilibrium layers
      INTEGER, PARAMETER :: LM_REQ=3
!@var REQ_FAC/REQ_FAC_M factors for REQ layer pressures
      REAL*8, PARAMETER, DIMENSION(LM_REQ-1) ::
     *     REQ_FAC=(/ .5d0, .2d0 /)               ! edge
      REAL*8, PARAMETER, DIMENSION(LM_REQ) ::
     *     REQ_FAC_M=(/ .75d0, .35d0, .1d0 /),    ! mid-points
     *     REQ_FAC_D=(/ .5d0,  .3d0,  .2d0 /)     ! delta

!@var PL00, PMIDL00, PDSIGL00, AML00 press (mb), mid-pressure (mb),
!@+        mass (kg/m2) for mean profile
!@var PEDNL00 edge pressure for mean profile (mb)
      REAL*8, DIMENSION(LM+LM_REQ) ::
     &     PL00, PMIDL00, PDSIGL00, AML00, BYAML00
      REAL*8, DIMENSION(LM+LM_REQ+1) :: PEDNL00

!**** Model control parameters:

!@var ij_debug: if i > 0, print out some extra info on bad ij box
      integer, dimension(2) :: ij_debug = (/ 0 , 1 /)

!**** Diagnostic control parameters
!@dbparam Kradia if -1 save data for, if 1|2 do   inst|adj forcing run
      integer :: Kradia=0,iu_rad

!**** Main atmospheric prognostic variables
!@var U,V east-west, and north-south velocities (m/s)
!@var T potential temperature (referenced to 1 mb) (K)
!@var Q specific humidity (kg water vapor/kg air)
!@var WM cloud liquid water amount (kg water/kg air)
#ifdef BLK_2MOM
!@var WMICE cloud ice amount (kg water/kg air)
#endif
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: U
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: V
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: T
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: Q
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: WM
#ifdef BLK_2MOM
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: WMICE
#endif
!@var P surface pressure (hecto-Pascals - PTOP)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: P

      real*8, parameter :: temperature_istart1=250. ! not used

!**** Boundary condition arrays:
!@var ZATMO: surface elevation (m)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: ZATMO

C**** Some helpful arrays (arrays should be L first)
!@var  PLIJ  Surface pressure: P(I,J) or PSF-PTOP (mb)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PLIJ
!@var  PDSIG  Surface pressure * DSIG(L) (mb)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PDSIG
!@var  AM  Air mass of each box (kg/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: AM     ! PLIJ*DSIG(L)*100/grav
!@var  BYAM  1/Air mass (m^2/kg)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: BYAM
!@var  PMID  Pressure at mid point of box (mb)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PMID    ! SIG(L)*PLIJ+PTOP
!@var  PK   PMID**KAPA
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PK
!@var  PEUP  Pressure at lower edge of box (incl. surface) (mb)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PEDN  ! SIGE(L)*PLIJ+PTOP
!@var  PEK  PEUP**KAPA
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PEK
!@var  SQRTP  square root of P (used in diagnostics)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: SQRTP
#ifdef etc_subdd
!@var  TTROPO  Temperature at mid point of tropopause level, extra subdaily
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: TTROPO
#endif
!@var  PTROPO  Pressure at mid point of tropopause level (mb)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: PTROPO
!@var  LTROPO  Tropopause layer
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: LTROPO

C**** module should own dynam variables used by other routines
!@var PTOLD pressure at beginning of dynamic time step (for clouds)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)    :: PTOLD
!@var SD_CLOUDS vert. integrated horizontal convergence (for clouds)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: SD_CLOUDS
!@var GZ geopotential height (for Clouds and Diagnostics)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: GZ
!@var DPDX_BY_RHO,DPDY_BY_RHO (pressure gradients)/density at L=1
      REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: DPDX_BY_RHO,DPDY_BY_RHO
!@var DPDX_BY_RHO_0,DPDY_BY_RHO_0 surface (pressure gradients)/density
      REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: DPDX_BY_RHO_0
      REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: DPDY_BY_RHO_0

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PHI

!@var PUA,PVA,SDA,PS save PU,PV,SD,P for hourly tracer advection
!@var MB Air mass array for tracers (before advection)
!@var MA Air mass array for tracers (updated during advection)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: PUA,PVA,SDA,MB,MA
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: PS

!@var DKE change in KE due to dissipation (SURF/DC/MC) (m^2/s^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: DKE
!@var KEA KE on the A grid (m^2/s^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: KEA ! ke on A grid
!@var UALIJ,VALIJ U,V on the A grid (m/s)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: UALIJ,VALIJ
!@var WSAVE vertical velocity (m/s)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: WSAVE

!@var SRFP actual surface pressure (hecto-Pascals)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: SRFP

      TARGET :: SRFP

      END MODULE ATM_COM

      SUBROUTINE ALLOC_ATM_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
      USE CONSTANT, only : GRAV
      USE FILEMANAGER
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID,HALO_UPDATE,READT_PARALLEL
     &     ,hassouthpole,hasnorthpole
      USE RESOLUTION, ONLY : IM,JM,LM,PSFMPT
      USE ATM_COM, ONLY : temperature_istart1
      USE ATM_COM, ONLY : ZATMO,P,U,V,T,Q,WM
#ifdef BLK_2MOM
     *  ,WMICE
#endif
      USE ATM_COM, ONLY : 
     &     PLIJ,PDSIG,AM,BYAM,PMID,PK,
     &     PEDN,PEK,SD_CLOUDS,GZ,PHI,
     &     PUA,PVA,SDA,MB,MA,DKE,KEA,
     &     UALIJ,VALIJ,WSAVE,SRFP,
     &     SQRTP,PTROPO,LTROPO,PS,PTOLD,
#ifdef etc_subdd
     &     TTROPO,
#endif
     &     DPDX_BY_RHO,DPDY_BY_RHO,DPDX_BY_RHO_0,DPDY_BY_RHO_0
#ifdef CUBED_SPHERE
       use GEOM, only : geom_cs
#else
       USE GEOM, only : geom_b
#endif
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid
      INTEGER :: iu_TOPO
      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: I, J, I_0, I_1, J_1, J_0
      INTEGER :: IER

      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO

      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      J_0 = grid%J_STRT
      J_1 = grid%J_STOP

C****
C**** CALCULATE SPHERICAL GEOMETRY
C****
#ifdef CUBED_SPHERE
      call geom_cs
#else
      CALL GEOM_B
#endif

      ALLOCATE(ZATMO(I_0H:I_1H,J_0H:J_1H), STAT = IER)

      ALLOCATE(P(I_0H:I_1H,J_0H:J_1H), STAT = IER)

      ALLOCATE(U(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)

      ALLOCATE(V(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)

      ALLOCATE(T(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)

      ALLOCATE(Q(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)

      ALLOCATE(WM(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)
#ifdef BLK_2MOM
      ALLOCATE(WMICE(I_0H:I_1H,J_0H:J_1H,LM), STAT = IER)
#endif


      U(:,:,:)=0.
      V(:,:,:)=0.
      T(:,:,:)=temperature_istart1  ! will be changed to pot.temp later
      Q(:,:,:)=3.D-6
      P(:,:)=PSFMPT
      WM    (:,:,:)=0.
#ifdef BLK_2MOM
      WMICE (:,:,:)=0.
#endif

      call openunit("TOPO",iu_TOPO,.true.,.true.)
      CALL READT_PARALLEL(grid,iu_TOPO,NAMEUNIT(iu_TOPO),ZATMO ,5) ! Topography
      ZATMO(I_0:I_1,J_0:J_1) = ZATMO(I_0:I_1,J_0:J_1)*GRAV   ! Geopotential
      CALL HALO_UPDATE(GRID, ZATMO)
C**** Check polar uniformity
      if(hassouthpole(grid)) then
        do i=2,im
          if (zatmo(i,1).ne.zatmo(1,1)) then
            print*,"Polar topography not uniform, corrected",i,1
     *           ,zatmo(i,1),zatmo(1,1)
            zatmo(i,1)=zatmo(1,1)
          end if
        end do
      end if
      if(hasnorthpole(grid)) then
        do i=2,im
          if (zatmo(i,jm).ne.zatmo(1,jm)) then
            print*,"Polar topography not uniform, corrected",i,jm
     *           ,zatmo(i,jm),zatmo(1,jm)
            zatmo(i,jm)=zatmo(1,jm)
          end if
        end do
      end if
      call closeunit(iu_TOPO)

      ! K-I-J arrays
      ALLOCATE ( PLIJ(LM,I_0H:I_1H,J_0H:J_1H), 
     $           PDSIG(LM,I_0H:I_1H,J_0H:J_1H),
     $             AM(LM,I_0H:I_1H,J_0H:J_1H),  
     $           BYAM(LM,I_0H:I_1H,J_0H:J_1H),
     $           PMID(LM,I_0H:I_1H,J_0H:J_1H),    
     $           PK(LM,I_0H:I_1H,J_0H:J_1H),
     $           PEDN(LM+1,I_0H:I_1H,J_0H:J_1H), 
     $            PEK(LM+1,I_0H:I_1H,J_0H:J_1H),
     $   STAT = IER)
      ! I-J-K arrays
      ALLOCATE( SD_CLOUDS(I_0H:I_1H,J_0H:J_1H,LM),  
     $                 GZ(I_0H:I_1H,J_0H:J_1H,LM), 
     $                PHI(I_0H:I_1H,J_0H:J_1H,LM), 
     $                PUA(I_0H:I_1H,J_0H:J_1H,LM), 
     $                PVA(I_0H:I_1H,J_0H:J_1H,LM), 
     $                SDA(I_0H:I_1H,J_0H:J_1H,LM),  
     $                MB(I_0H:I_1H,J_0H:J_1H,LM), 
     $                MA(I_0H:I_1H,J_0H:J_1H,LM), 
     $                DKE(I_0H:I_1H,J_0H:J_1H,LM), 
     $                KEA(I_0H:I_1H,J_0H:J_1H,LM), 
     $                UALIJ(LM,I_0H:I_1H,J_0H:J_1H), 
     $                VALIJ(LM,I_0H:I_1H,J_0H:J_1H), 
     $              WSAVE(I_0H:I_1H,J_0H:J_1H,LM-1), 
     $   STAT = IER)

      ! I-J arrays
      ALLOCATE(  SQRTP(I_0H:I_1H,J_0H:J_1H), 
     $          SRFP(I_0H:I_1H,J_0H:J_1H),
     $          PTROPO(I_0H:I_1H,J_0H:J_1H),
     $          LTROPO(I_0H:I_1H,J_0H:J_1H),  
#ifdef etc_subdd
     $          TTROPO(I_0H:I_1H,J_0H:J_1H),   ! extra subdaily
#endif
     $          PTOLD(I_0H:I_1H,J_0H:J_1H),
     $          DPDX_BY_RHO(I_0H:I_1H,J_0H:J_1H), 
     $          DPDY_BY_RHO(I_0H:I_1H,J_0H:J_1H),
     $          DPDX_BY_RHO_0(I_0H:I_1H,J_0H:J_1H), 
     $          DPDY_BY_RHO_0(I_0H:I_1H,J_0H:J_1H),
     $          PS(I_0H:I_1H,J_0H:J_1H),
     $   STAT = IER)

! correct or wrong, but being static all arrays were initialized
! to zero by default. They have to be initialized to something now 
! to avoid floating point exceptions...
      DPDX_BY_RHO(I_0H:I_1H,J_0H:J_1H) = 0.d0
      DPDY_BY_RHO(I_0H:I_1H,J_0H:J_1H) = 0.d0
      DPDX_BY_RHO_0(I_0H:I_1H,J_0H:J_1H) = 0.d0
      DPDY_BY_RHO_0(I_0H:I_1H,J_0H:J_1H) = 0.d0

      SD_CLOUDS(I_0H:I_1H,J_0H:J_1H,1:LM) = 0.d0

      END SUBROUTINE ALLOC_ATM_COM

      SUBROUTINE io_atm(kunit,iaction,ioerr)
!@sum  io_model reads and writes model variables to file
!@auth Gavin Schmidt
      USE ATM_COM
      USE RESOLUTION, only : IM,JM,LM
      USE MODEL_COM, only : IOWRITE,IOREAD,LHEAD
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, only: PACK_DATA,UNPACK_DATA,AM_I_ROOT
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "MODEL01"
!@var U_glob Work array for parallel I/O
!@var V_glob Work array for parallel I/O
!@var T_glob Work array for parallel I/O
!@var Q_glob Work array for parallel I/O
!@var WM_glob Work array for parallel I/O
#ifdef BLK_2MOM
!@var WMICE_glob Work array for parallel I/O
#endif
!@var P_glob Work array for parallel I/O
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: !(IM,JM,LM)
     &     U_glob,V_glob,T_glob,Q_glob,WM_glob
#ifdef BLK_2MOM
     *,WMICE_glob
#endif
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: P_glob!(IM,JM)
      integer :: img, jmg, lmg

      MODULE_HEADER(lhead+1:80) = 'R8 dim(im,jm,lm):u,v,t, p(im,jm),'//
     *  ' dim(im,jm,lm):q,MliqW'

      if(am_i_root()) then
         img = IM
         jmg = JM
         lmg = LM
      else
         img = 1
         jmg = 1
         lmg = 1
      end if
      allocate(U_glob(img,jmg,lmg))
      allocate(V_glob(img,jmg,lmg))
      allocate(T_glob(img,jmg,lmg))
      allocate(P_glob(IM,JM))
      allocate(Q_glob(img,jmg,lmg))
      allocate(WM_glob(img,jmg,lmg))
#ifdef BLK_2MOM
      allocate(WMICE_glob(img,jmg,lmg))
#endif

      SELECT CASE (IACTION)
      CASE (:IOWRITE) ! output to end-of-month restart file
        CALL PACK_DATA(grid, U, U_GLOB)
        CALL PACK_DATA(grid, V, V_GLOB)
        CALL PACK_DATA(grid, T, T_GLOB)
        CALL PACK_DATA(grid, Q, Q_GLOB)
        CALL PACK_DATA(grid, WM, WM_GLOB)
#ifdef BLK_2MOM
        CALL PACK_DATA(grid, WMICE, WMICE_GLOB)
#endif
        CALL PACK_DATA(grid, P, P_GLOB)
        IF (AM_I_ROOT())
     &    WRITE (kunit,err=10) MODULE_HEADER,U_glob,V_glob,T_glob,
     &                         P_glob,Q_glob,WM_glob
#ifdef BLK_2MOM
     &   ,WMICE_glob
#endif
      CASE (IOREAD:)          ! input from restart file
        if ( AM_I_ROOT() ) then
          READ (kunit,err=10) HEADER,U_glob,V_glob,T_glob,
     &                           P_glob,Q_glob,WM_glob
#ifdef BLK_2MOM
     &   ,WMICE_glob
#endif
          IF (HEADER(1:LHEAD).ne.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
          END IF
        end if
        CALL UNPACK_DATA(grid, U_GLOB, U)
        CALL UNPACK_DATA(grid, V_GLOB, V)
        CALL UNPACK_DATA(grid, T_GLOB, T)
        CALL UNPACK_DATA(grid, Q_GLOB, Q)
        CALL UNPACK_DATA(grid, WM_GLOB, WM)
#ifdef BLK_2MOM
        CALL UNPACK_DATA(grid, WMICE_GLOB, WMICE)
#endif
        CALL UNPACK_DATA(grid, P_GLOB, P)
      END SELECT
      call freemem
      RETURN
 10   IOERR=1
      call freemem
      RETURN
      contains
      subroutine freemem
        deallocate(U_glob)
        deallocate(V_glob)
        deallocate(T_glob)
        deallocate(P_glob)
        deallocate(Q_glob)
        deallocate(WM_glob)
#ifdef BLK_2MOM
        deallocate(WMICE_glob)
#endif
      end subroutine freemem
      END SUBROUTINE io_atm

#ifdef NEW_IO
ccc was not sure where to dump these routines ... IA
      module conserv_diags
      implicit none

      contains
      subroutine declare_conserv_diags( grid, fid, name_dims )
      use domain_decomp_atm, only : dist_grid, get
      use pario, only : defvar
      implicit none
      type (dist_grid), intent(in) :: grid
      integer ::  fid
      character(len=*) :: name_dims
      integer :: i_0h, i_1h, j_0h, j_1h
      integer :: ier
      real*8, allocatable :: buf(:,:)
      call get( grid, j_strt_halo=j_0h, j_stop_halo=j_1h,
     &     i_strt_halo=i_0h, i_stop_halo=i_1h )
      allocate( buf(i_0h:i_1h,j_0h:j_1h), stat=ier)
      call defvar(grid, fid, buf, trim(name_dims))     
      deallocate( buf )
      end subroutine declare_conserv_diags

      subroutine dump_conserv_diags( grid, fid, name, conserv )
      use domain_decomp_atm, only : dist_grid, get
      use pario, only : write_dist_data
      implicit none
      type (dist_grid), intent(in) :: grid
      integer ::  fid
      character(len=*) :: name
      external :: conserv
      integer :: i_0h, i_1h, j_0h, j_1h
      integer :: ier
      real*8, allocatable :: buf(:,:)
      call get( grid, j_strt_halo=j_0h, j_stop_halo=j_1h,
     &     i_strt_halo=i_0h, i_stop_halo=i_1h )
      allocate( buf(i_0h:i_1h,j_0h:j_1h), stat=ier)
      call conserv(buf)
      call write_dist_data(grid,fid,trim(name),buf)
      deallocate( buf )
      end subroutine dump_conserv_diags

      end module conserv_diags

      subroutine def_rsf_atm(fid)
!@sum  def_rsf_model defines U,V,T,P,Q,WM array structure in restart files
!@auth M. Kelley
!@ver  beta
      use atm_com
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      use conserv_diags
      implicit none
      integer fid   !@var fid file id
      character(len=20) :: ijlstr
      ijlstr='(dist_im,dist_jm,lm)'
      call defvar(grid,fid,u,'u'//ijlstr)
      call defvar(grid,fid,v,'v'//ijlstr)
      call defvar(grid,fid,t,'t'//ijlstr)
      call defvar(grid,fid,q,'q'//ijlstr)
      call defvar(grid,fid,wm,'wm'//ijlstr)
      call defvar(grid,fid,p,'p(dist_im,dist_jm)')
#ifdef BLK_2MOM
      call defvar(grid,fid,wmice,'wmice'//ijlstr)
#endif
      call declare_conserv_diags( grid, fid, 'watmo(dist_im,dist_jm)' )
      call declare_conserv_diags( grid, fid, 'ekatmo(dist_im,dist_jm)' )
      call declare_conserv_diags( grid, fid, 'epatmo(dist_im,dist_jm)' )
      call declare_conserv_diags( grid, fid, 'ewatmo(dist_im,dist_jm)' )
      return
      end subroutine def_rsf_atm

      subroutine new_io_atm(fid,iaction)
!@sum  new_io_model read/write U,V,T,P,Q,WM arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : iowrite,ioread
      use atm_com
      use domain_decomp_atm, only: grid
      use pario, only : write_dist_data,read_dist_data
      use conserv_diags
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      external conserv_WM, conserv_KE, conserv_PE, conserv_EWM
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid, fid, 'u', u)
        call write_dist_data(grid, fid, 'v', v)
        call write_dist_data(grid, fid, 't', t)
        call write_dist_data(grid, fid, 'p', p)
        call write_dist_data(grid, fid, 'q', q)
        call write_dist_data(grid, fid, 'wm', wm)
#ifdef BLK_2MOM
        call write_dist_data(grid, fid, 'wmice', wmice)
#endif
        call dump_conserv_diags( grid, fid, 'watmo', conserv_WM )
        call dump_conserv_diags( grid, fid, 'ekatmo', conserv_KE )
        call dump_conserv_diags( grid, fid, 'epatmo', conserv_PE )
        call dump_conserv_diags( grid, fid, 'ewatmo', conserv_EWM  )
      case (ioread)             ! input from restart file
        call read_dist_data(grid, fid, 'u', u)
        call read_dist_data(grid, fid, 'v', v)
        call read_dist_data(grid, fid, 't', t)
        call read_dist_data(grid, fid, 'p', p)
        call read_dist_data(grid, fid, 'q', q)
        call read_dist_data(grid, fid, 'wm', wm)
#ifdef BLK_2MOM
        call read_dist_data(grid, fid, 'wmice', wmice)
#endif
      end select
      return
      end subroutine new_io_atm
#endif /* NEW_IO */
