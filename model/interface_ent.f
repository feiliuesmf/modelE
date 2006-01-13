

      module interface_ent



      contains



!from reth:

!    ws(0,2) - saturated amount of water in canopy (m)

!from retp:
!    shc(0,2)



      subroutine veg_conductance(
     &   cnc
     &       ,gpp
     &       ,trans_sw       !nyk
     &       ,betad          ! evaporation efficiency
     &       ,tp_can          ! canopy temperature C
     &       ,qv
     &       ,dts
     &       )
      real*8 cnc,gpp,trans_sw,betad,tp_can,qv,dts

      end subroutine veg_conductance


      ! dummy progs (so far..)

      subroutine reset_veg_to_defaults( reset_prognostic )
      logical, intent(in) :: reset_prognostic

      end subroutine reset_veg_to_defaults

      subroutine io_ent    (kunit,iaction,ioerr)
      integer kunit             !@var kunit unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
!@var ioerr 1 (or -1) if there is (or is not) an error in i/o
      integer, intent(inout) :: ioerr

      end subroutine io_ent


      subroutine init_module_ent(iniENT,grid,jday,dxyp)
      USE DOMAIN_DECOMP, ONLY : DIST_GRID
      logical :: iniENT
      TYPE (DIST_GRID), INTENT(IN) :: grid
      integer jday
      real*8 dxyp(:)
      

      end subroutine init_module_ent



      subroutine ent_get_value( i, j,
     &     canopy_holding_capacity,
     &     canopy_heat_capacity,
     &     fraction_of_vegetated_soil
     &     )
      integer, intent(in) :: i, j
      real*8, optional :: canopy_holding_capacity
      real*8, optional :: canopy_heat_capacity
      real*8, optional :: fraction_of_vegetated_soil
      !----------


      if ( present(canopy_holding_capacity) ) then
        !canopy_holding_capacity = .0001d0 * alai
      endif

      if ( present(canopy_heat_capacity) ) then
        !aa=ala(1,i0,j0)
        !canopy_heat_capacity=(.010d0+.002d0*aa+.001d0*aa**2)*shw
      endif

      if ( present(fraction_of_vegetated_soil) ) then
        !fraction_of_vegetated = sum( vdata(i,j,2:9 )
      endif


      end subroutine ent_get_value


      subroutine ent_reset_veg_to_defaults( reset_prognostic )
!@sum fills vegetation arrays with some default vegetation
!@+   for unusual experiments (paleo climate, different location
!@+   of continents, etc.) 
      logical reset_prognostic
      end subroutine ent_reset_veg_to_defaults

      subroutine ent_update_crops(jyear)
      integer, intent (in) :: jyear

      end subroutine ent_update_crops

      end module interface_ent
