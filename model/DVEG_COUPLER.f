      module DVEG_COUPLER
!@sum contains routines which provide interface between GCM and 
!@+   dynamic vegetation module
      implicit none

!!! so far this is just an example how to pass GCM data to Ent module

      contains

      subroutine step_dveg(dt)
!@sum extracts all the data necessary for dynamic vegetation time step,
!@+   passes them to the dynamic vegetation routine which performs
!@+   one time step, then copies updated values into GCM arrays if 
!@+   necessary
      implicit none
      use ENT_INTERFACE, only : ent_step
      use MODEL_COM, only :
!@var dt dynamic vegetation time step (s)
      real*8, intent(in) :: dt

      call ent_step(
     &     dt                   ! time step (s)
     &     )

      end subroutine step_dveg

      subroutine init_dveg
!@sum extracts all the data necessary for initialization of dynamic
!@+   vegetation module from GCM, then passes them to the initialization
!@+   routine
      use ENT, only : ent_init
      use GEOM : dxyp


      call ent_init(
     &     dxyp                 ! area of a grid cell (m^2)
     &     )

      end subroutine init_dveg

      end module DVEG_COUPLER
