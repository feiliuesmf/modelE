#include "rundeck_opts.h"

      module irrigate_crop

!@sum Contains routines to read in irrigation data and to compute the
!@sum potential irrigation rate
!@auth M.J. Puma, R. Ruedy

      implicit none
      private

      public irrigate_flux
      public init_irrigate
      public irrigate_extract

      SAVE

!@var Monthly irrigation arrays read in to compute daily irrigation
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: irr_month_0,irr_month_1
!@var iu_irrigate is unit number for irrigation climatologies
      integer :: iu_irrigate
!@dbparam irrig_cycl determines whether prescribed irrigation data
!@+       vary year to year:  0=no; 1=yes
      integer :: irrig_cycl


      contains

!-----------------------------------------------------------------------

      subroutine init_irrigate()
      use filemanager, only : openunit
      use param, only : sync_param
      USE DOMAIN_DECOMP_ATM, ONLY : grid, get
      implicit none

      integer :: I_0,I_1,J_0,J_1
      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER, i,j

      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO

      ALLOCATE( irr_month_0   ( I_0H:I_1H , J_0H:J_1H ),
     &          irr_month_1   ( I_0H:I_1H , J_0H:J_1H ),
     &          STAT=IER)

      call openunit("IRRIG",iu_irrigate,.true.,.true.)
      call sync_param("irrig_cycl",irrig_cycl)
      call GET(grid,I_STRT=I_0,I_STOP=I_1,J_STRT=J_0,J_STOP=J_1)

!**** Compute irrigation rate
      call irrigate_flux()

      do j=J_0,J_1
         do i=I_0,I_1
!****      Initialize actual water and energy values
           call irrigate_extract(i,j)
         enddo
      enddo

      end subroutine init_irrigate

!-----------------------------------------------------------------------

      subroutine irrigate_flux()
!@sum  Calculates daily irrigation from monthly irrigation data.  Routine
!@sum  sets imon, irr_month_0, and irr_month_1 depending on jday and
!@sum  jyear (cyclical case does not need jyear).
      use model_com, only : jday,Itime,jmon,JDmidOfM,itimei,JMperY,jyear
      USE DOMAIN_DECOMP_ATM, ONLY : grid, get, am_i_root
     &                          ,READT_PARALLEL,REWIND_PARALLEL
     &                          ,READ_PARALLEL, MREAD_PARALLEL
     &                          ,BACKSPACE_PARALLEL
      USE FILEMANAGER, only : NAMEUNIT
      USE FLUXES,only :irrig_water_pot
      USE CONSTANT, only : rhow,teeny,shw

      implicit none

!@var m_yr_mnth: title of the monthly irrigation data (YYYYMM)
       integer :: m_yr_mnth
!@var JDLAST julian day that irrigate_flux was last called
      integer, SAVE :: jdlast=0
!@var IMON counter to identify data for later of two months
!@var needed for interpolation (irr_month_1) ->
!@var IMON-1 is for irr_month_0
      integer, SAVE :: imon = 0
      integer :: i,j,m,lstmon
      integer :: I_0,I_1,J_0,J_1
      real*8 :: FRAC

      call GET(grid,I_STRT=I_0,I_STOP=I_1,J_STRT=J_0,J_STOP=J_1)

      if (irrig_cycl == 0) then
!**** Non-cyclical case: data vary from year to year
         if(jdlast == 0)then ! first call
!           Position within the irrigation file based on jyear and jmon
            lstmon = jyear*100 + jmon
            imon = jmon + 1
            if (jday < JDmidOFM(imon))then
               imon = jmon
               lstmon=lstmon-1
               if(jmon==1) lstmon = lstmon-100+12
            endif

            call READ_PARALLEL(m_yr_mnth, iu_irrigate)

            do while (m_yr_mnth < lstmon)
               call READ_PARALLEL(m_yr_mnth, iu_irrigate)
            enddo

            if(m_yr_mnth /= lstmon ) then
              call stop_model(
     &        'Read error: irrig. data not available for start date'
     &         ,255)
            endif

!           Read in irr_month_0
            call BACKSPACE_PARALLEL( iu_irrigate )
            call MREAD_PARALLEL(grid,iu_irrigate,
     &                 NAMEUNIT(iu_irrigate),m,irr_month_0)

            if(AM_I_ROOT())then
               write(6,*) 'Read in irrigate for month,',jmon,m
            endif

         else

           if (jday /= jdlast+1) then !have a new year: reset imon 13->1
             imon=imon-12
             go to 330
           end if

           if (jday <= JDmidOFM(imon)) go to 330

           imon=imon+1          ! read in new month of climatological data
           irr_month_0 = irr_month_1

         endif

!        Read in irr_month_1
         CALL MREAD_PARALLEL(grid,iu_irrigate,
     &          NAMEUNIT(iu_irrigate),m,irr_month_1)
         if(AM_I_ROOT())then
            WRITE(6,*) 'Read in irrigate for next month,',jmon,m
         endif

 330     jdlast=jday

!**** End of non-cyclical case

      else

!**** Cyclical case: data does not from year to year
         if (JDLAST == 0) then ! need to read in first month climatology
           IMON=1          ! IMON=January
           if (JDAY < 16)  then ! JDAY in Jan 1-15, first month is Dec
             print *, "iu_irrigate, IPOS =", iu_irrigate, 12
             CALL READT_PARALLEL(grid,iu_irrigate,NAMEUNIT(iu_irrigate),
     &                        irr_month_0,12)
             CALL REWIND_PARALLEL( iu_irrigate )
           else            ! JDAY is in Jan 16 to Dec 16, get first month
  620        IMON=IMON+1
             IF (JDAY > JDmidOFM(IMON) .and. IMON <= 12) GO TO 620
             CALL READT_PARALLEL
     &       (grid,iu_irrigate,NAMEUNIT(iu_irrigate),irr_month_0,IMON-1)
             if (IMON == 13)  CALL REWIND_PARALLEL( iu_irrigate )
           endif
         ELSE                      ! Do we need to read in second month?
           IF (JDAY /= JDLAST+1) THEN ! Check that data is read in daily
             IF (JDAY /= 1 .or. JDLAST /= 365) THEN
            WRITE (6,*) 'Incorrect values in irrigate,JDAY,JDLAST=',JDAY
     &           ,JDLAST
               call stop_model(
     &           'ERROR READING IN SETTING IRRIGATION CLIMATOLOGY',255)
             END IF
             IMON=IMON-12          ! New year
             GO TO 630
           END IF
           IF (JDAY <= JDmidOFM(IMON)) GO TO 630
           IMON=IMON+1          ! read in new month of climatological data
           irr_month_0 = irr_month_1
           IF (IMON == 13) CALL REWIND_PARALLEL( iu_irrigate )
         END IF
         CALL READT_PARALLEL(grid,iu_irrigate,NAMEUNIT(iu_irrigate),
     &                    irr_month_1,1)
 630     JDLAST=JDAY
      endif

!**** Interpolate daily irrigation depth to the current day
      FRAC = REAL(JDmidOFM(IMON)-JDAY,KIND=8)/
     &           (JDmidOFM(IMON)-JDmidOFM(IMON-1))
      if(am_I_root())
     &write(6,*) 'interp. irrigate to it,day,fr0',itime,jday,frac
!**** Compute irrigation rate
      do j=J_0,J_1
         do i=I_0,I_1
!****      Input irrigation values in m/s
           irrig_water_pot(i,j)= ( FRAC*irr_month_0(i,j)
     &                          + (1.d0-FRAC)*irr_month_1(i,j) )
!****      Check to make sure no negative irrigation rates
           if( irrig_water_pot(i,j) < 0.d0 ) then
              irrig_water_pot(i,j) = 0.d0
           endif

         enddo
      enddo

      end subroutine irrigate_flux

!-----------------------------------------------------------------------

      subroutine irrigate_extract(i,j)
      USE FLUXES, only :irrig_water_pot,irrig_water_act,irrig_energy_act
     &                 ,irrig_gw,irrig_gw_energy

      USE LAKES_COM, only : mwl,gml,tlake,mldlk,flake
      USE GEOM, only : imaxj,axyp
      USE CONSTANT, only : rhow,teeny,shw
      USE MODEL_COM, only : dtsrc
      USE GHY_COM, only : fearth
      USE sle001, only : tp

      implicit none
      integer,intent(in):: i,j

!@var Mass of irrigation water withdrawn from rivers/groundwater [kg]
      real*8 :: m_irr_pot
!@var Mass of water available for irrigation [kg]
      real*8 :: m_avail
!@var Temperature of abstracted irrigation water [deg C]
      real*8 :: T_irr, T_irr_g
!@var Vegetated and bare fractions of the grid cell
      real*8 :: fv, fb
!@param Minimum depth of lake for abstracting irrigation water [m]
      real*8,parameter :: hlake_min = 1.d0
!@param Conversion factor from m^3 water/(m^2 grid cell) to kg water
      real*8 :: m_to_kg
!@param Flag for external irrigation source (from deep aquifers - not modeled) 
      integer,parameter :: flag_irrig_grndwat = 1

      m_to_kg = rhow*axyp(i,j)

      call get_fb_fv( fb, fv, i, j )

!***  Irrigation based on potential values in m/s per area grid cell 
      if ( (irrig_water_pot(i,j) > teeny)  .and. 
     &     (fearth(i,j)          > teeny)  .and.
     &     (fv                   > teeny)) then

         m_irr_pot = irrig_water_pot(i,j)* m_to_kg *dtsrc

         if (flake(i,j) > 0.d0) then
            m_avail = mwl(i,j) - hlake_min*flake(i,j)*m_to_kg
            m_avail = max(m_avail, 0.d0)
            T_irr = tlake(i,j)
         else
            m_avail = mwl(i,j)
            T_irr = gml(i,j)/(mwl(i,j)*shw+teeny)
         endif
!        Check these limits !!!!
         T_irr = max(T_irr, 0.d0)
         T_irr_g = max(tp(1,2),0.d0)

!***     Set actual irrigation rates and update mwl and gml (if necessary)
         if (m_avail <= teeny) then

            if(flag_irrig_grndwat == 0) then
               irrig_water_act(i,j) = 0.d0
               irrig_energy_act(i,j) = 0.d0
               irrig_gw(i,j) = 0.d0
               irrig_gw_energy(i,j) = 0.d0
            else
               irrig_water_act(i,j) = irrig_water_pot(i,j)
               irrig_energy_act(i,j)= irrig_water_act(i,j)*shw*
     &                                 T_irr_g*rhow
               irrig_gw(i,j)        = irrig_water_act(i,j)
               irrig_gw_energy(i,j) = irrig_energy_act(i,j)
            endif

         elseif (m_avail >= m_irr_pot) then

            mwl(i,j) = mwl(i,j) - irrig_water_pot(i,j)*m_to_kg*dtsrc
            gml(i,j) = gml(i,j) - irrig_water_pot(i,j)
     &                               *m_to_kg*dtsrc*shw*T_irr
            irrig_water_act(i,j) = irrig_water_pot(i,j)
            irrig_energy_act(i,j)=irrig_water_act(i,j)*shw*T_irr*rhow
            irrig_gw(i,j) = 0.d0
            irrig_gw_energy(i,j) = 0.d0
            
         else !!! (m_avail < m_irr_pot)

            mwl(i,j) = mwl(i,j) - m_avail
            gml(i,j) = gml(i,j) - m_avail*shw*T_irr

            if(flag_irrig_grndwat == 0) then
               irrig_water_act(i,j) = m_avail / m_to_kg / dtsrc
               irrig_energy_act(i,j) = irrig_water_act(i,j)*shw*T_irr
     &                                 *rhow
               irrig_gw(i,j) = 0.d0
               irrig_gw_energy(i,j) = 0.d0
            else
               irrig_water_act(i,j) = irrig_water_pot(i,j)
               irrig_energy_act(i,j)=irrig_water_act(i,j)*shw*T_irr*rhow
               irrig_gw(i,j)=irrig_water_pot(i,j)-m_avail/m_to_kg/dtsrc
               irrig_gw_energy(i,j) = irrig_gw(i,j)*shw*T_irr_g*rhow
            endif

         endif

      else
         irrig_water_act(i,j) = 0.d0
         irrig_energy_act(i,j) = 0.d0
         irrig_gw(i,j) = 0.d0
         irrig_gw_energy(i,j) = 0.d0
      endif

!!!!! HACK HACK test test
      if( (irrig_water_act(i,j)-irrig_water_pot(i,j)) > 1.d-30 ) then
         write(6,*) 'error in irrigation: act> pot'
         write(6,*) 'i,j,irrig_water_act(i,j), irrig_water_pot(i,j):',
     &              i,j,irrig_water_act(i,j), irrig_water_pot(i,j)
         write(6,*) 'm_avail, m_irr_pot:', m_avail, m_irr_pot
         irrig_water_act(i,j)=irrig_water_pot(i,j)
      endif
!!!!!
      end subroutine irrigate_extract

!-----------------------------------------------------------------------

      end module irrigate_crop
