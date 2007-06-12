#include "rundeck_opts.h"
      MODULE COSMO_SOURCES

      USE TRACER_COM
      SAVE
!@var be7_src_3d Source function for 7Be (kg tracer)/( (kg air)/m^2 s )
      real*8 be7_src_3d(im,jm,lm)
!@var be7_src_param global multiplier of be7_src_3d to match obs
      real*8 :: be7_src_param=1    !default value

      END MODULE COSMO_SOURCES


      SUBROUTINE read_Be_source
!@sum reads in cosmogenic Be7 source from appropriate file
!@auth C Salyk
      USE CONSTANT, only : avog
      USE COSMO_SOURCES, only: be7_src_3d
      USE TRACER_COM
      USE GEOM, only: dxyp
      USE FILEMANAGER, only: openunit,closeunit
      IMPLICIT NONE
      real*8, dimension(jm,lm) :: ibe, ibm
      integer iuc, j,i,l

C**** constants used in file to make numbers neater
      real*8 :: tfacti, tfact2

C**** Open source file
      call openunit('BE7_COSMO', iuc, .false., .true.)
      read(iuc,*) tfacti
      read(iuc,*) ibe
C**** ibe has units atoms/g/s
      call closeunit(iuc)

C**** convert from atoms/g/s to (kg tracer)/ (kg air/m^2) /s
      do l=1,lm; do j=1,jm ; do i=1,im
        be7_src_3d(i,j,l)=ibe(j,l)*dxyp(j)*(tr_mm(n_Be7)*tfacti/avog)
      end do ; end do ; end do

      END SUBROUTINE read_Be_source



