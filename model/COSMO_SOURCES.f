#include "rundeck_opts.h"
      MODULE COSMO_SOURCES
      
      USE TRACER_COM

      real*8 be7_src_3d(im,jm,lm)
      real*8 :: be7_src_param=1    !default value

      END MODULE COSMO_SOURCES


      SUBROUTINE read_Be_source
!@sum reads in cosmogenic Be7 source from appropriate file
!@auth C Salyk 
      USE COSMO_SOURCES, only: be7_src_3d
      USE TRACER_COM
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
        call closeunit(iuc)
 
       do l=1,lm; do j=1,jm; do i=1,im
****  convert from atoms/g/s 
          temp_src_3d(i,j,l)=ibe(j,l)*dxyp(j)*tr_mm(n_Be7)/avo
        end do; end do; end do

**** loop over longitude
        do i=1,im
        be7_src_3d(:,:,:)=temp_src_3d(:,:,:) * tfacti
        end do

        END SUBROUTINE read_Be_source



