      module lsm
      use ent_mod
      use sle001, only : ngm, imt, nlsn, LS_NFRAC, advnc
      implicit none

      integer, parameter :: im=72, jm=46
      integer, parameter :: j_0h=1, j_1h=jm

      type t_lsm_state
        real*8, pointer :: w_ij(:,:,:,:)
        real*8, pointer :: ht_ij(:,:,:,:)
        integer, pointer :: nsn_ij    (:, :, :)     
        real*8, pointer :: dzsn_ij   (:, :, :, :)
        real*8, pointer :: wsn_ij    (:, :, :, :)
        real*8, pointer :: hsn_ij    (:, :, :, :)
        real*8, pointer :: fr_snow_ij(:, :, :)         
        real*8, pointer :: Qf_ij(:,:)
        type (entcelltype_public), pointer :: entcells(:,:)
      end type t_lsm_state

      type t_lsm_bc
        real*8, pointer ::
     &       fearth(:,:),
     &       top_index_ij(:,:),                 
     &       top_dev_ij(:,:),                   
     &       dz_ij(:,:,:),                   
     &       q_ij(:,:,:,:),              
     &       qk_ij(:,:,:,:),             
     &       sl_ij(:,:)        
      end type t_lsm_bc

      contains

      subroutine lsm_init_state( s, fearth, jday, jyear )
      type(t_lsm_state) :: s
      real*8 :: fearth(:,:)
      integer jday, jyear
      integer i,j

      allocate(
     &     s%w_ij(0:ngm,ls_nfrac,im,j_0h:j_1h),
     &     s%ht_ij(0:ngm,ls_nfrac,im,j_0h:j_1h),
     &     s%nsn_ij(     2,im,j_0h:j_1h),
     &     s%dzsn_ij(nlsn,2,im,j_0h:j_1h),
     &     s%wsn_ij(nlsn,2,im,j_0h:j_1h),
     &     s%hsn_ij(nlsn,2,im,j_0h:j_1h),
     &     s%fr_snow_ij(2,im,j_0h:j_1h),
     &     s%Qf_ij(im,j_0h:j_1h) )
      
      read(952)
     &       s%w_ij(0:ngm,1:ls_nfrac,:,:),         
     &       s%ht_ij(0:ngm,1:ls_nfrac,:,:),        
     &       s%nsn_ij    (1:2, :,:),              
     &       s%dzsn_ij   (1:nlsn, 1:2, :,:),      
     &       s%wsn_ij    (1:nlsn, 1:2, :,:),      
     &       s%hsn_ij    (1:nlsn, 1:2, :,:),      
     &       s%fr_snow_ij(1:2, :,:),
     &       s%Qf_ij(:,:)

      allocate( s%entcells(im,j_0h:j_1h) )
      call ent_cell_nullify( s%entcells )

      do j=1,jm
        do i=1,im
          if ( fearth(i,j) > 0.d0 )
     &         call ent_cell_construct( s%entcells(i,j) )
        enddo
      enddo

      call set_vegetation_data( s%entcells,
     &     IM, JM, 1, IM, 1, JM, jday, jyear )

      end subroutine lsm_init_state


      subroutine lsm_init_bc( bc )
      type (t_lsm_bc) :: bc

      allocate(
     &     bc%fearth      (im, j_0h:j_1h),
     &     bc%top_index_ij(im, j_0h:j_1h),
     &     bc%top_dev_ij  (im, j_0h:j_1h),
     &     bc%dz_ij       (im, j_0h:j_1h, ngm),
     &     bc%q_ij        (im, j_0h:j_1h, imt, ngm),
     &     bc%qk_ij       (im, j_0h:j_1h, imt, ngm),
     &     bc%sl_ij       (im, j_0h:j_1h) )

        read(953)
     &       bc%fearth(:,:),
     &       bc%top_index_ij(:,:),                 
     &       bc%top_dev_ij(:,:),                   
     &       bc%dz_ij(:,:,1:ngm),                   
     &       bc%q_ij(:,:,1:imt,1:ngm),              
     &       bc%qk_ij(:,:,1:imt,1:ngm),             
     &       bc%sl_ij(:,:)        

      end subroutine lsm_init_bc


      subroutine get_forcings(
     &         force_Ca              ,
     &         force_cos_zen_angle   ,
     &         force_vis_rad         ,
     &         force_direct_vis_rad  ,
     &         force_prec_ms         ,
     &         force_eprec_w         ,
     &         force_sprec_ms        ,
     &         force_seprec_w        ,
     &         force_srheat          ,
     &         force_trheat          ,
     &         force_ts              ,
     &         force_qs              ,
     &         force_ps              ,
     &         force_rhosrf          ,
     &         force_cdh             ,
     &         force_qm1             ,
     &         force_ws              ,
     &         force_pbl_args_ws0    ,
     &         force_pbl_args_tprime ,
     &         force_pbl_args_qprime )
      real*8
     &         force_Ca (1:im,1:jm)             ,
     &         force_cos_zen_angle (1:im,1:jm)  ,
     &         force_vis_rad (1:im,1:jm)        ,
     &         force_direct_vis_rad (1:im,1:jm) ,
     &         force_prec_ms (1:im,1:jm)        ,
     &         force_eprec_w (1:im,1:jm)        ,
     &         force_sprec_ms (1:im,1:jm)       ,
     &         force_seprec_w (1:im,1:jm)       ,
     &         force_srheat (1:im,1:jm)         ,
     &         force_trheat (1:im,1:jm)         ,
     &         force_ts (1:im,1:jm)             ,
     &         force_qs (1:im,1:jm)             ,
     &         force_ps (1:im,1:jm)             ,
     &         force_rhosrf (1:im,1:jm)         ,
     &         force_cdh (1:im,1:jm)            ,
     &         force_qm1 (1:im,1:jm)            ,
     &         force_ws (1:im,1:jm)             ,
     &         force_pbl_args_ws0 (1:im,1:jm)   ,
     &         force_pbl_args_tprime (1:im,1:jm),
     &         force_pbl_args_qprime (1:im,1:jm)

      print *,"reading forcings"

      read(951)
     &         force_Ca (1:im,1:jm)             ,
     &         force_cos_zen_angle (1:im,1:jm)  ,
     &         force_vis_rad (1:im,1:jm)        ,
     &         force_direct_vis_rad (1:im,1:jm) ,
     &         force_prec_ms (1:im,1:jm)        ,
     &         force_eprec_w (1:im,1:jm)        ,
     &         force_sprec_ms (1:im,1:jm)       ,
     &         force_seprec_w (1:im,1:jm)       ,
     &         force_srheat (1:im,1:jm)         ,
     &         force_trheat (1:im,1:jm)         ,
     &         force_ts (1:im,1:jm)             ,
     &         force_qs (1:im,1:jm)             ,
     &         force_ps (1:im,1:jm)             ,
     &         force_rhosrf (1:im,1:jm)         ,
     &         force_cdh (1:im,1:jm)            ,
     &         force_qm1 (1:im,1:jm)            ,
     &         force_ws (1:im,1:jm)             ,
     &         force_pbl_args_ws0 (1:im,1:jm)   ,
     &         force_pbl_args_tprime (1:im,1:jm),
     &         force_pbl_args_qprime (1:im,1:jm)

      print *,"forcings ok"

      end subroutine get_forcings


      subroutine lsm_initialize( dt_in )
      use sle001, only : hl0, dt
      real*8 dt_in

      call ent_initialize(
     &     do_soilresp=.true.
     &     ,do_phenology=.false.
     &     ,do_patchdynamics=.false.
     &     )

      dt = dt_in
      call hl0

      end subroutine lsm_initialize


      subroutine lsm_run( s, bc, jday, jyear )
      use sle001, only : tp,tbcs
      type(t_lsm_state) :: s
      type (t_lsm_bc) :: bc
      integer, intent(in) :: jday, jyear
      real*8 fb,fv
      integer i,j

      real*8 ::
     &     force_Ca (1:im,1:jm)             ,
     &     force_cos_zen_angle (1:im,1:jm)  ,
     &     force_vis_rad (1:im,1:jm)        ,
     &     force_direct_vis_rad (1:im,1:jm) ,
     &     force_prec_ms (1:im,1:jm)        ,
     &     force_eprec_w (1:im,1:jm)        ,
     &     force_sprec_ms (1:im,1:jm)       ,
     &     force_seprec_w (1:im,1:jm)       ,
     &     force_srheat (1:im,1:jm)         ,
     &     force_trheat (1:im,1:jm)         ,
     &     force_ts (1:im,1:jm)             ,
     &     force_qs (1:im,1:jm)             ,
     &     force_ps (1:im,1:jm)             ,
     &     force_rhosrf (1:im,1:jm)         ,
     &     force_cdh (1:im,1:jm)            ,
     &     force_qm1 (1:im,1:jm)            ,
     &     force_ws (1:im,1:jm)             ,
     &     force_pbl_args_ws0 (1:im,1:jm)   ,
     &     force_pbl_args_tprime (1:im,1:jm),
     &     force_pbl_args_qprime (1:im,1:jm)

          call get_forcings(
     &         force_Ca (1:im,1:jm)             ,
     &         force_cos_zen_angle (1:im,1:jm)  ,
     &         force_vis_rad (1:im,1:jm)        ,
     &         force_direct_vis_rad (1:im,1:jm) ,
     &         force_prec_ms (1:im,1:jm)        ,
     &         force_eprec_w (1:im,1:jm)        ,
     &         force_sprec_ms (1:im,1:jm)       ,
     &         force_seprec_w (1:im,1:jm)       ,
     &         force_srheat (1:im,1:jm)         ,
     &         force_trheat (1:im,1:jm)         ,
     &         force_ts (1:im,1:jm)             ,
     &         force_qs (1:im,1:jm)             ,
     &         force_ps (1:im,1:jm)             ,
     &         force_rhosrf (1:im,1:jm)         ,
     &         force_cdh (1:im,1:jm)            ,
     &         force_qm1 (1:im,1:jm)            ,
     &         force_ws (1:im,1:jm)             ,
     &         force_pbl_args_ws0 (1:im,1:jm)   ,
     &         force_pbl_args_tprime (1:im,1:jm),
     &         force_pbl_args_qprime (1:im,1:jm) )


      call update_vegetation_data( s%entcells,
     &     im, jm, 1, im, 1, jm, jday, jyear )

      loop_j: do j = 1,jm
        loop_i: do i = 1,im

          if ( bc%fearth(i,j) <= 0.d0 ) cycle loop_i


          ! really fb, fv are not needed for Ent, but just in case...
          call ent_get_exports( s%entcells(i,j),
     &         fraction_of_vegetated_soil=fv
     &         )
          fb = 1.d0 - fv

          call advnc(
!-------------- Ent specific
     &         s%entcells(i,j), force_Ca(i,j),
     &         force_cos_zen_angle(i,j), force_vis_rad(i,j),
     &         force_direct_vis_rad(i,j),
     &         s%Qf_ij(i,j),
!-------------- old vegetation scheme (not implemented at the moment)
!     &         vegcell,
!-------------- prognostic vars
     &         s%w_ij(0:ngm,1:LS_NFRAC,i,j),         
     &         s%ht_ij(0:ngm,1:LS_NFRAC,i,j),        
     &         s%nsn_ij    (1:2, i, j),              
     &         s%dzsn_ij   (1:nlsn, 1:2, i, j),      
     &         s%wsn_ij    (1:nlsn, 1:2, i, j),      
     &         s%hsn_ij    (1:nlsn, 1:2, i, j),      
     &         s%fr_snow_ij(1:2, i, j),          
!-------------- BC's    
     &         bc%top_index_ij(i, j),                 
     &         bc%top_dev_ij(i, j),                   
     &         bc%dz_ij(i,j,1:ngm),                   
     &         bc%q_ij(i,j,1:imt,1:ngm),              
     &         bc%qk_ij(i,j,1:imt,1:ngm),             
     &         bc%sl_ij(i,j),                         
     &         fb,                                 
     &         fv,                        
!-------------- forcings         
     &         force_prec_ms (i,j)        ,
     &         force_eprec_w (i,j)        ,
     &         force_sprec_ms (i,j)       ,
     &         force_seprec_w (i,j)       ,
     &         force_srheat (i,j)         ,
     &         force_trheat (i,j)         ,
     &         force_ts (i,j)             ,
     &         force_qs (i,j)             ,
     &         force_ps (i,j)             ,
     &         force_rhosrf (i,j)         ,
     &         force_cdh (i,j)            ,
     &         force_qm1 (i,j)            ,
     &         force_ws (i,j)             ,
     &         force_pbl_args_ws0 (i,j)   ,
     &         force_pbl_args_tprime (i,j),
     &         force_pbl_args_qprime (i,j) )

          write(936,*) i,j, tbcs, s%w_ij(0:ngm,1:LS_NFRAC,i,j),
     &     s%ht_ij(0:ngm,1:LS_NFRAC,i,j)



          print *,i,j,tp(1,1),tp(0,2)
        enddo loop_i
      enddo loop_j

      end subroutine lsm_run

      subroutine set_vegetation_data( entcells,
     &     im, jm, i0, i1, j0, j1, jday, year )
!@sum read standard GISS vegetation BC's and pass them to Ent for
!@+   initialization of Ent cells. Halo cells ignored, i.e.
!@+   entcells should be a slice without halo
      use ent_prescribed_drv, only : init_canopy_physical,prescr_vegdata
      use ent_prescr_veg, only : prescr_calc_shc,prescr_calcconst
      type(entcelltype_public), intent(out) :: entcells(I0:I1,J0:J1)
      integer, intent(in) :: im, jm, i0, i1, j0, j1, jday, year
      !Local variables
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: vegdata !cohort
      real*8, dimension(N_BANDS,N_COVERTYPES,I0:I1,J0:J1) :: albedodata !patch, NOTE:snow
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: laidata  !cohort
      real*8, dimension(N_COVERTYPES) :: hdata    !cohort
      real*8, dimension(N_COVERTYPES) :: nmdata    !cohort
      real*8, dimension(N_COVERTYPES,N_DEPTH) :: rootprofdata !Root fraction of veg type.
      real*8, dimension(N_COVERTYPES) :: popdata !Dummy population density:  0-bare soil, 1-vegetated
      real*8, dimension(N_COVERTYPES) :: dbhdata !Diameter at breast height for woody veg.(cm)
      real*8, dimension(N_COVERTYPES) :: craddata !Crown radius (m)
      real*8, dimension(N_COVERTYPES,N_BPOOLS,I0:I1,J0:J1) :: cpooldata !Carbon pools in individuals
      integer, dimension(N_COVERTYPES) :: soildata ! soil types 1-bright 2-dark
      real*8, dimension(N_SOIL_TEXTURES,I0:I1,J0:J1) :: soil_texture
      real*8, dimension(I0:I1,J0:J1) :: Ci_ini,CNC_ini,Tcan_ini,Qf_ini
      real*8, dimension(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS,
     &     I0:I1,J0:J1):: Tpool_ini
      !-----Local---------
      integer i,j
      real*8 heat_capacity

      call prescr_calcconst()

      call prescr_vegdata(jday, year, 
     &     IM,JM,I0,I1,J0,J1,vegdata,albedodata,laidata,hdata,nmdata,
     &     popdata,dbhdata,craddata,cpooldata,rootprofdata,
     &     soildata,soil_texture,Tpool_ini,.true.)
      call init_canopy_physical(I0, I1, J0, J1,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini)

      !Translate gridded data to Entdata structure
      !GISS data:  a patch per vegetation cover fraction, one cohort per patch
      call ent_cell_set(entcells, vegdata, popdata, laidata,
     &     hdata, dbhdata, craddata, cpooldata, nmdata, rootprofdata, 
     &     soildata, albedodata, soil_texture,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini, Tpool_ini)

      ! just in case, do nothing, just set heat capacities
      call ent_prescribe_vegupdate(entcells)

      end subroutine set_vegetation_data


      subroutine update_vegetation_data( entcells,
     &     im, jm, i0, i1, j0, j1, jday, year )
!@sum read standard GISS vegetation BC's and pass them to Ent for
!@+   initialization of Ent cells. Halo cells ignored, i.e.
!@+   entcells should be a slice without halo
      use ent_prescribed_drv, only:
     &     prescr_get_laidata,prescr_veg_albedodata
      !use ent_prescr_veg, only: prescr_get_laidata,prescr_veg_albedodata
      type(entcelltype_public), intent(out) :: entcells(I0:I1,J0:J1)
      integer, intent(in) :: im, jm, i0, i1, j0, j1, jday, year
      !Local variables
      real*8, dimension(N_BANDS,N_COVERTYPES,I0:I1,J0:J1) :: albedodata !patch, NOTE:snow
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: laidata  !cohort
      !-----Local---------
      integer i,j

      call prescr_get_laidata(jday,JM,I0,I1,J0,J1,laidata)
      call prescr_veg_albedodata(jday,JM,I0,I1,J0,J1,albedodata)

      call ent_prescribe_vegupdate(entcells
     &     ,laidata=laidata(2:2+N_PFT-1,:,:)
     &     ,albedodata=albedodata(:,2:2+N_PFT-1,:,:)
     &     )

      end subroutine update_vegetation_data

      end module lsm

      MODULE TRIDIAG_MOD
!@sum TRIDIAG_MOD contains subroutine TRIDIAG
      PRIVATE
      PUBLIC TRIDIAG

      contains

      SUBROUTINE TRIDIAG(A,B,C,R,U,N)
!@sum  TRIDIAG  solves a tridiagonal matrix equation (A,B,C)U=R
!@auth Numerical Recipes
!@ver  1.0
      IMPLICIT NONE
c      INTEGER, PARAMETER :: NMAX = 8000  !@var NMAX workspace
      INTEGER, INTENT(IN):: N         !@var N    dimension of arrays
      REAL*8, INTENT(IN) :: A(N)   !@var A    coefficients of u_i-1
      REAL*8, INTENT(IN) :: B(N)   !@var B    coefficients of u_i
      REAL*8, INTENT(IN) :: C(N)   !@var C    coefficients of u_i+1
      REAL*8, INTENT(IN) :: R(N)   !@var R    RHS vector
      REAL*8, INTENT(OUT):: U(N)   !@var U    solution vector
      REAL*8 :: BET                !@var BET  work variable
      REAL*8 :: GAM(Size(A))       !@var GAM  work array
      INTEGER :: J                 !@var J    loop variable

c      IF ( N > NMAX )
c     &     call stop_model("TRIDIAG: N > NMAX, increase NMAX",255)
      BET=B(1)
      IF (BET.eq.0) call stop_model("TRIDIAG: DENOMINATOR = ZERO",255)
      U(1)=R(1)/BET
      DO J=2,N
        GAM(J)=C(J-1)/BET
        BET=B(J)-A(J)*GAM(J)
        IF (BET.eq.0) call stop_model("TRIDIAG: DENOMINATOR = ZERO",255)
        U(J)=(R(J)-A(J)*U(J-1))/BET
      END DO
      DO J=N-1,1,-1
        U(J)=U(J)-GAM(J+1)*U(J+1)
      END DO
      RETURN
      END SUBROUTINE TRIDIAG

      end MODULE TRIDIAG_MOD


      program lsm_standalone
      use lsm
      use ent_mod
      implicit none
      type (t_lsm_bc) :: bc
      type(t_lsm_state) :: lsm_state
      integer jday, jyear
      integer :: itime,itime0=1,itime1=2 !48

      ! set jday, jyear as in default GCM setup
      jday = 335
      jyear = 1949

!need to initialize GHY and ENT here ...
      call lsm_initialize( 1800.d0 )

      call lsm_init_bc( bc )
      call lsm_init_state( lsm_state, bc%fearth, jday, jyear )

      ! main time loop
      do itime=itime0,itime1

        call lsm_run( lsm_state, bc, jday, jyear )

      enddo

      end
