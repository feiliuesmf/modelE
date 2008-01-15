#include "rundeck_opts.h"
      MODULE COSMO_SOURCES

      USE TRACER_COM
      USE DOMAIN_DECOMP, only : GRID, get

      SAVE
!@var be7_ and be10_src_3d Source functions for 7Be & 10Be (kg tracer)/( (kg air)/m^2 s )
      real*8, allocatable, dimension(:,:,:) :: be7_src_3d, be10_src_3d
!@var be7_src_param global multiplier of be7_src_3d to match obs
      real*8 :: be7_src_param=1    !default value
      real*8, allocatable, dimension(:,:) :: BE7W_acc, BE7D_acc 
      INTEGER :: J_1H, J_0H, variable_phi

      END MODULE COSMO_SOURCES


      SUBROUTINE init_cosmo
!      IMPLICIT NONE
      USE COSMO_SOURCES, only : be7_src_3d, be10_src_3d,
     $     BE7W_acc, BE7D_acc, variable_phi
      USE PARAM, only : sync_param
      USE TRACER_COM
      USE DOMAIN_DECOMP, only : GRID, get
      INTEGER :: J_1H, J_0H

      call sync_param("variable_phi",variable_phi)

      CALL GET(grid,J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      
      ALLOCATE(be7_src_3d (1:im, J_0H:J_1H, 1:lm))
      ALLOCATE(be10_src_3d (1:im, J_0H:J_1H, 1:lm))
      ALLOCATE(BE7W_acc (1:im, J_0H:J_1H))
      ALLOCATE(BE7D_acc (1:im, J_0H:J_1H))
 
      if (variable_phi .eq. 0) call read_Be_source_noAlpha
      print*, "variable_phi = ", variable_phi

      if (variable_phi .eq. 1) call read_Be_source
      print*, "variable_phi = ", variable_phi
         
!      if (variable_phi .eq. 2) call update_annual_phi
!      print*, "variable_phi = ", variable_phi
      
      if (variable_phi .eq. 3) call update_daily_phi   
      print*, "variable_phi = ", variable_phi

      END SUBROUTINE init_cosmo
      
      SUBROUTINE read_Be_source_noAlpha
!@sum reads in cosmogenic Be7 source from appropriate "old" (no alpha particles) versions 
!@sum of the Beer production files
!@auth C Salyk
      USE CONSTANT, only : avog
      USE COSMO_SOURCES, only: be7_src_3d
      USE TRACER_COM
      USE GEOM, only: dxyp
      USE DOMAIN_DECOMP, only : GRID, get
      USE FILEMANAGER, only: openunit,closeunit
      IMPLICIT NONE
      integer, parameter :: layers=23
      real*8, dimension(jm,lm) :: ibe, ibm
      integer iuc, j,i,l
C**** constants used in file to make numbers neater
      real*8 :: tfacti, tfact2
      INTEGER :: J_1H, J_0H

      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

C**** Open source file
      call openunit('BE7_COSMO_NO_ALPHA', iuc, .false., .true.)
      read(iuc,*) tfacti
      read(iuc,*) ibe
C**** ibe has units atoms/g/s
      call closeunit(iuc)

C**** convert from atoms/g/s to (kg tracer)/ (kg air/m^2) /s
      do l=1,lm; do j=J_0H,J_1H; do i=1,im
        be7_src_3d(i,j,l)=ibe(j,l)*dxyp(j)*(tr_mm(n_Be7)*tfacti/avog)
      end do ; end do ; end do

      END SUBROUTINE read_Be_source_noAlpha



      SUBROUTINE read_Be_source
!@sum reads in cosmogenic Be7 source from appropriate file
!@auth C Salyk
      USE CONSTANT, only : avog
      USE COSMO_SOURCES, only: be7_src_3d, be10_src_3d
      USE TRACER_COM
      USE GEOM, only: dxyp
      USE DOMAIN_DECOMP, only : GRID, get
      USE FILEMANAGER, only: openunit,closeunit
      IMPLICIT NONE
      integer, parameter :: layers=23
      real*8, dimension(jm,lm) :: ibe, ibm
      real*8, dimension(jm,lm) :: ibe_10
      integer iuc, j,i,l
C**** constants used in file to make numbers neater
      real*8 :: tfacti, tfact2, tfacti_10
      INTEGER :: J_1H, J_0H

      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
c      ALLOCATE(ibe (1:im, J_0H:J_1H))
c      ALLOCATE(ibm (1:im, J_0H:J_1H))
c      ALLOCATE(ibe_10 (1:im, J_0H:J_1H))

C**** Open source file
      print*, "layers = ", layers
      print*, "about to open Be7_cosmo"
      call openunit('BE7_COSMO', iuc, .false., .true.)
      print*, "opened Be7_cosmo file"
      read(iuc,*) tfacti
      print*, "read tfacti"
      print*, "tfacti = ", tfacti
      print*, "reading ibe"
      read(iuc,*) ibe
      print*, "read ibe"
C**** ibe has units atoms/g/s
      call closeunit(iuc)
      print*, "closed be7 cosmo file"
      
C**** convert from atoms/g/s to (kg tracer)/ (kg air/m^2) /s
      print*, "converting"
      do l=1,lm; do j=J_0H,J_1H ; do i=1,im
         be7_src_3d(i,j,l)=ibe(j,l)*dxyp(j)*(tr_mm(n_Be7)*tfacti/avog)
      end do ; end do ; end do

C     repeat for Be10:
      print*, "about to open Be10_cosmo"
      call openunit('BE10_COSMO', iuc, .false., .true.)
      print*, "opened Be10_cosmo file"
      read(iuc,*) tfacti_10
      print*, "read tfacti_10"
      read(iuc,*) ibe_10
      print*, "read ibe_10"
C**** ibe has units atoms/g/s
      call closeunit(iuc)
      print*, "closed be10 cosmo file"
      
C**** convert from atoms/g/s to (kg tracer)/ (kg air/m^2) /s
      print*, "converting"
      do l=1,lm; do j=J_0H,J_1H ; do i=1,im
         be10_src_3d(i,j,l)=ibe_10(j,l)*dxyp(j)*(tr_mm(n_Be10)*tfacti_10
     $        /avog)
!         if (i .eq. 1) then
!            print*, "i,j,k =", i, j, l
!         end if
      end do ; end do ; end do
      print*, "finished converting"

      
      END SUBROUTINE read_Be_source

      
!      SUBROUTINE update_annual_phi
!@sum interpolates betw. diff. Be7 production values for each year to get correct Be7 
!@production corresponding to phi value for a given year.
!@auth C Field
!      phi_yr = model_yr + 0.5
!      do i = 1,nyrs
!         if (phi_yr .eq. year(i)) then
!            phi = phi_record(i)
!         end if
!      end do
      
!      do i = 1,nphi
!      if ((be7_phi(i) .lt. phi) .and.(be7_phi(i+1) .gt. phi)) then
            
!           delta_phi = phi - be7_phi(i)
            
!     slope=(be7_src_prod_alpha(:,:,:,i+1)-be7_src_prod_alpha(:,:,:,i)))/((be7_phi(i+1))-(be7_phi(i)))
      
!     new_prod = be7_src_prod_alpha(:,:,:,i)+(slope(:,:,:)*delta_phi)
            
            
!         else if (be7_phi(i) .eq. phi) then
!            new_prod = be7_src_prod_alpha(:,:,:,i)
!         end if
!      end do
      
      
!      END SUBROUTINE update_annual_phi


 
      SUBROUTINE update_daily_phi
!@sum to be used with doing model runs for Usoskin experiments, Jan. 2005 - Feb. 2005
!@sum reads in phi timeseries, Be7 production values and geomagnetic (Pc) values.
!@sum Production values are in atoms/g/s.
!@auth C Field
      USE FILEMANAGER, only: openunit,closeunit
      USE GEOM, only: dxyp
      USE MODEL_COM, only : jday, jmon, itime
      USE CONSTANT, only : avog
      USE TRACER_COM
      USE DOMAIN_DECOMP, only : GRID, get
      USE COSMO_SOURCES, only : be7_src_3d

      IMPLICIT NONE
      character title*80
      integer, parameter :: npress=23, npc=41, nphi=31
      integer, parameter :: nmonths=2, ndays=59
      real :: day(ndays), daily_phi(ndays)
      real :: phi_day
      real :: pc_list(npc), phi_list(nphi)
      real*8 :: be7_prod(1:npress,1:npc,1:nphi)
      real*8 :: be7_scr(1:npress,1:npc)
      real*8 pc_month
      real*4, dimension(im,jm,nmonths) :: pc_table
      real*8, dimension(lm) :: new_prod, new_prod1, new_prod2, scr_prod
       
      real :: phi_hi, phi_low, pc_hi, pc_low
      real :: delta_pc, delta_phi, delta_phi_1, delta_pc_3
      real :: slope(npress), slope_1(npress)
      real :: slope_2(npress), slope_3(npress), slope_4(npress) 
      integer :: iuc, i, j, k, m, n, iphi, ipc
      INTEGER :: J_1H, J_0H

      CALL GET(grid,J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

! files needed:
!     BE_DAILY_PROD = Be23_TEST_all_UsoskinLinear_v3.dat = contains production values for all vertical layers, all pc values and all phi values (units are atoms/g/s).
!     BE_PC_VALUES = Pc_output_JF_Data_TEST.dat = contains pc values for each I/J box and each month from 1955 to 2006 

      print*, "reading files for daily phi"
      call openunit('BE7_DAILY_PROD', iuc, .false., .true.)
      print*, "opened"
      read(iuc,*,end=20) (((be7_prod(i,j,k), i=1,npress), j=1,npc), k=1
     $     ,nphi)
 !     print*, be7_prod(1,1,1)
 20   call closeunit(iuc)

      print*, "reading files for SCR production"
      call openunit('BE7_SCR_PROD', iuc, .false., .true.)
      print*, "opened"
      read(iuc,*,end=40) ((be7_scr(i,j), i=1,npress), j=1,npc)
 !     print*, be7_scr(1,1)
 40   call closeunit(iuc)
      
      print*, "reading"
      call openunit('PC_LOOKUP', iuc, .true., .true.)
      do k=1,nmonths
         read(iuc) title,pc_table(:,:,k)
         print*,"read ",title
      end do
!      print*, "pc_table (1,1,1) = ", pc_table(1,1,1) 
      call closeunit(iuc)  
      
       print*, "reading PC list"
      call openunit('PC_LIST', iuc, .false., .true.)
         read(iuc,*) (pc_list(i), i=1,npc)
!         print*, "pc_list (1) = ", pc_list(1) 
      call closeunit(iuc)

       print*, "reading Phi list"
      call openunit('PHI_LIST', iuc, .false., .true.)
         read(iuc,*) (phi_list(i), i=1,nphi)
         do j = 1,nphi
        print*, "phi_list = ", phi_list(j)
        end do

      call closeunit(iuc)

      print*, "reading PHI_BY_DAY"
      call openunit('PHI_BY_DAY', iuc, .false., .true.)
      do i = 1,ndays
         read(iuc,*) day(i), daily_phi(i)
!         print*, "daily_phi (1) = ", daily_phi(1)
      end do
      call closeunit(iuc)  

      

c     for each day, calculate Be7 production based on that month's pc
c     value and that day's phi value      

      if (jday > ndays) then
        phi_day = daily_phi(ndays) ! ? default 
      else
        phi_day = daily_phi(jday)
      end if
        print*, "phi_day = ", phi_day

      do m = 2,nphi
         if ((phi_day .lt. phi_list(m)) .and. (phi_day .ge. phi_list(m-1
     $        ))) then
            phi_low = phi_list(m-1)
            phi_hi = phi_list(m)
            iphi = m-1
            exit 
         end if
      end do
!      print*, "phi values: ", phi_low, phi_hi, iphi
      

      do j = J_0H,J_1H
         do i = 1,IM
            pc_month =  pc_table(i,j,jmon)
            if ((i .eq. 10) .and. (j .eq. J_1H)) then
               print*, "pc_month = ", pc_month, jmon,i,j
            end if
            
            do n = 1,npc
               if ((pc_month .lt. pc_list(n)) .and. (pc_month
     $              .ge. pc_list(n-1))) then
                  pc_low = pc_list(n-1) 
                  pc_hi = pc_list(n) 
                  ipc = n-1
                  exit
               end if
            end do
            if ((i .eq. 10) .and. (j .eq. J_1H)) then
               print*, "pc values (n. pole): ", pc_month, pc_low, pc_hi,
     $              ipc
            end if
            

c     there are 4 possible cases:
c     1. no interpolation needed
c     2. interpolation needed over pc values
c     3. interpolation needed over phi values
c     4. interpolation needed over both pc and phi
c     all four scenarios are accounted for in the following code:

c     first interpolate between the 2 phis at pc_lo:
               delta_phi_1 = phi_day - phi_low
               
               slope_1 = (be7_prod(:,ipc,iphi+1) - be7_prod(:,ipc,iphi)
     $              )/(phi_hi - phi_low)  
               
               new_prod1(:) = be7_prod(:,ipc,iphi) + (slope_1
     $              *delta_phi_1)
               if ((i .eq. 1) .and. (j .eq. 1)) then
                  print*, "new_prod1 = ", new_prod1(1)
               end if

c     next interpolate betw. the 2 phis at pc_hi:
               slope_2 = (be7_prod(:,ipc+1,iphi+1) - be7_prod(:,ipc+1
     $              ,iphi))/(phi_hi - phi_low)
               
               new_prod2(:) = be7_prod(:,ipc+1,iphi) + (slope_2
     $              *delta_phi_1)
               if ((i .eq. 1) .and. (j .eq. 1)) then
                  print*, "new_prod2 = ", new_prod2(1)
               end if
               

c     now interpolate betw. the 2 interpolated values that you've just
c     created:
               delta_pc_3 = pc_month - pc_low
               
               slope_3 = (new_prod2(:) - new_prod1(:))/(pc_hi -
     $              pc_low)
               
               new_prod(:) = new_prod1(:) + (slope_3*delta_pc_3)
               if ((i .eq. 1) .and. (j .eq. 1)) then
                  print*, "new_prod = ", new_prod(1)
               end if

! calculate Be7 production from SCR if JDAY = Jan.20:
c      if (jday .eq. 20) then
c         slope_4 = (be7_scr(:,ipc+1) - be7_scr(:,ipc))/(pc_hi - pc_low
c     $        )
c         scr_prod(:) = be7_scr(:,ipc) + (slope_4*delta_pc_3)
c         
c         new_prod = new_prod + scr_prod
c      end if
            

C**** convert from atoms/g/s to (kg tracer)/ (kg air/m^2) /s
            do k=1,npress
!               print*, "converting atoms/g/s to kg/kg"
              be7_src_3d(i,j,k)=new_prod(k)*dxyp(j)*(tr_mm(n_Be7)
     $              /avog)
              if ((i .eq. 10) .and. (j .eq. 46) .and. (k .eq. 1)) then
                 print*, "DXYP = ", dxyp(j)
                 print*, "tr_mm(n_Be7) = ", tr_mm(n_Be7)
                 print*, " "
              end if
              if ((i .eq. 37) .and. (j .eq. J_1H)) then
                 print*, "itime = ", itime, "be7_src_3d (n. pole) = ",
     $                be7_src_3d(37,J_1H,k)
              end if
           end do 
            
         end do
      end do
   
      print*, "closing update_daily_phi"

      END SUBROUTINE update_daily_phi
      









      



