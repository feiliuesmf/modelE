#include "rundeck_opts.h"
      MODULE COSMO_SOURCES
      
      USE TRACER_COM

      real*8 be7_src_3d(im,jm,lm)

      END MODULE COSMO_SOURCES


      SUBROUTINE read_Be_source
!@sum reads in cosmogenic Be7 source from appropriate file
!@auth C Salyk 
      USE COSMO_SOURCES, only: be7_src_3d
      USE TRACER_COM
      USE DYNAMICS, only: am
      USE GEOM, only: dxyp
      USE FILEMANAGER, only: openunit,closeunit   
      IMPLICIT NONE
      SAVE


      real*8, dimension(jm,lm) :: ibe, ibm
      integer iuc, j,i,l
      real*8, dimension(im,jm,lm) :: temp_src_3d

!constants used in file to make numbers neater
      real*8 :: tfacti, tfact2
      real*8, parameter :: avo=6.02d23

**** Open source file
        call openunit('BE7_COSMO', iuc, .false.)
        read(iuc,*) tfacti
        read(iuc,*) ibe
**** ibe has units atoms/g/s
        read(iuc,*) tfact2
        read(iuc,*) ibm
        call closeunit(iuc)

        do l=1,lm; do j=1,jm; do i=1,im
****  convert from atoms/g/s to kg/s
          temp_src_3d(i,j,l)=ibe(j,l)*am(l,i,j)*dxyp(j)*tr_mm(n_Be7)/avo
        end do; end do; end do

**** loop over longitude
        do i=1,im
        be7_src_3d(:,:,:)=temp_src_3d(:,:,:) * tfacti
        end do

**** only for diagnostic purposes
!        open(unit=21, file='temp.out', form='formatted')
!        write(21,900) tfacti
! 900    format (E8.3)
!        write(21,901) tfact2
! 901    format(E10.5) 
!        write(21,902) ibe
! 902    format(23I5)
!       close(21)

        END SUBROUTINE read_Be_source



