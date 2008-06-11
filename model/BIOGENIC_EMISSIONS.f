#include "rundeck_opts.h"

      module biogenic_emis
      integer, parameter :: npolynb=20, ntype=16, nvegtype=74
!@dbparam base_isopreneX factor to tune the base isoprene emissions
!@+ globally when #defined BIOGENIC_EMISSIONS
      real*8 :: base_isopreneX=1.d0
      real*8, allocatable, dimension(:,:,:) :: baseisop
      real*8 isopcoeff(npolynb)
      end module biogenic_emis


      subroutine alloc_biogenic_emis(grid)
!@SUM  To alllocate arrays whose sizes now need to be determined
!@+    at run-time
!@auth G.Faluvegi
!@ver  1.0
      use domain_decomp, only : dist_grid, get
      use biogenic_emis, only:  baseisop,nvegtype
      use model_com, only: im

      IMPLICIT NONE

      type (dist_grid), intent(in) :: grid
      integer :: J_1H, J_0H

      call get( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
      allocate( baseisop(IM,J_0H:J_1H,nvegtype) )

      end subroutine alloc_biogenic_emis


      subroutine isoprene_emission(i,j,itype,emisop)

      use biogenic_emis
      use ghy_com, only: canopy_temp_ij  
               !canopy_temp_ij = array(im,jm) of canopy temperatures (C)
               !Same temperature for all veg types over grid cell.
      use rad_com, only : cosz1,cfrac
      use geom, only : bydxyp
      use pbl_drv, only : t_pbl_args
      use constant, only : rgas
      use tracers_drydep, only : xylai,ijreg,ijuse
      use constant, only : tf
    
      implicit none
 
      type (t_pbl_args) :: pbl_args
      integer :: inveg
      integer, intent(in) :: i,j,itype  
      real*8 :: tlai,embio,clight,biofit,tcorr,sfdens
      real*8, dimension(i,j) :: tmmpk
      real*8, intent(out) :: emisop

      select case(itype)

      case(1:3)      ! ocean, ocean ice, landice:

        emisop=0.d0  ! no emissions

      case(4)        ! land:

! Temperature of canpoy 

        ! estimated surface density
        sfdens=100.d0*pbl_args%psurf/(rgas*pbl_args%tgv) 
        tmmpk(i,j) = canopy_temp_ij(i,j) + tf

        emisop=0.d0
        tlai=0.d0

        do inveg=1,ijreg(i,j)
          tlai=tlai+xylai(i,j,inveg)*baseisop(i,j,inveg)
        end do

! Light correction

        if ((cosz1(i,j) > 0.).and.(tlai > 0.)) then ! Only calculate
        ! for grid cell with sunlight and isoprene-emitting vegetation

          embio=0.d0
          do inveg=1,ijreg(i,j)
            if (xylai(i,j,inveg)*baseisop(i,j,inveg) > 0.) then
              clight=biofit(isopcoeff,xylai(i,j,inveg),
     &               cosz1(i,j),cfrac(i,j))
              embio=embio+baseisop(i,j,inveg)*
     &              clight*ijuse(i,j,inveg)*1.d-3
            endif
          end do

! Temperature correction

          if(tmmpk(i,j) > tf) THEN
            emisop=tcorr(tmmpk(i,j))*embio
          else
            emisop=0.d0
          endif

        endif  

! emisop = kg C emitted from grid cell per second
! 2D interative isoprene source. Convert to units kg C/m2/s:

        emisop=emisop*bydxyp(j)

      end select

      return                                                          
      end                              



      subroutine rdisopcf                                              
!@sum These polynomial coefficients normally should be read
!@+  in from the file 'isoprene.coef'. I am hardcoding them
!@+  here, as the quickest way to make sure this is ESMF-
!@+  compliant, as I don't suspect we will commit this code.

      use biogenic_emis

      implicit none

! polynomial coefficients for biogenic isoprene emissions:
! From Guenther.

      if(npolynb /= 20)call stop_model('npolynb problem',255)

      isopcoeff(1) = -1.86E-01; isopcoeff(2) =  2.19E+00
      isopcoeff(3) =  2.12E+00; isopcoeff(4) = -2.43E-01
      isopcoeff(5) = -4.72E+00; isopcoeff(6) =  1.05E+01
      isopcoeff(7) =  3.77E-01; isopcoeff(8) = -3.05E+00
      isopcoeff(9) =  3.05E-01; isopcoeff(10)=  5.16E-01
      isopcoeff(11)=  2.72E+00; isopcoeff(12)= -4.82E+00
      isopcoeff(13)= -9.59E-04; isopcoeff(14)= -1.29E+00
      isopcoeff(15)=  9.37E-01; isopcoeff(16)= -3.31E-01
      isopcoeff(17)=  1.07E+00; isopcoeff(18)= -7.59E-02
      isopcoeff(19)= -3.01E-01; isopcoeff(20)= -4.07E-01

      return                                                            
      end                       


      subroutine rdisobase                                               
!@sum These baseline emissions factors normally should be read
!@+  in from the file 'isopemis.table'. I am hardcoding them
!@+  here, as the quickest way to make sure this is ESMF-
!@+  compliant, as I don't suspect we will commit this code.
!@+  Units are atoms C cm^-2 leaf s^-1
!@+  Construct the base emission for each grid box                         
!@+  Output is baseisop in kg C cm^-2 * surface area of cell (cm^2)
!@+  emitted in 1 hour time step                                           

      use model_com, only: im,jm
      use biogenic_emis
      use tracers_drydep, only : ijreg,ijland
      use constant, only   : avog
      use geom, only : dxyp,imaxj
      use domain_decomp, only : get, grid

      implicit none

      integer:: i,j,k,j_0,j_1
      real*8, dimension(nvegtype) ::  converT
      real*8 :: factor                       

      CALL GET(grid, J_STRT=J_0, J_STOP=J_1)

      if(nvegtype /= 74) call stop_model('nvegtype problem',255)

      convert( 0+1)=0.00E+00; convert( 1+1)=2.79E+12
      convert( 2+1)=8.71E+11; convert( 3+1)=2.79E+12
      convert( 4+1)=2.79E+12; convert( 5+1)=2.79E+12
      convert( 6+1)=2.79E+12; convert( 7+1)=2.79E+12
      convert( 8+1)=3.34E+12; convert( 9+1)=2.79E+12
      convert(10+1)=2.79E+12; convert(11+1)=2.79E+12
      convert(12+1)=2.79E+12; convert(13+1)=2.79E+12
      convert(14+1)=2.79E+12; convert(15+1)=2.79E+12
      convert(16+1)=2.79E+12; convert(17+1)=2.79E+12
      convert(18+1)=2.79E+12; convert(19+1)=2.79E+12
      convert(20+1)=1.67E+12; convert(21+1)=1.67E+12
      convert(22+1)=1.67E+12; convert(23+1)=1.39E+12
      convert(24+1)=3.34E+12; convert(25+1)=6.27E+12
      convert(26+1)=6.27E+12; convert(27+1)=3.34E+12
      convert(28+1)=0.93E+12; convert(29+1)=0.93E+12
      convert(30+1)=8.71E+11; convert(31+1)=8.71E+11
      convert(32+1)=2.49E+12; convert(33+1)=1.39E+12
      convert(34+1)=2.79E+12; convert(35+1)=2.79E+12
      convert(36+1)=8.71E+11; convert(37+1)=8.71E+11
      convert(38+1)=8.71E+11; convert(39+1)=8.71E+11
      convert(40+1)=0.93E+12; convert(41+1)=1.39E+12
      convert(42+1)=8.71E+11; convert(43+1)=2.79E+12
      convert(44+1)=1.67E+12; convert(45+1)=1.67E+12
      convert(46+1)=3.34E+12; convert(47+1)=2.79E+12
      convert(48+1)=9.41E+12; convert(49+1)=2.79E+12
      convert(50+1)=3.34E+12; convert(51+1)=3.34E+12
      convert(52+1)=2.79E+12; convert(53+1)=2.79E+12
      convert(54+1)=3.34E+12; convert(55+1)=2.79E+12
      convert(56+1)=4.18E+12; convert(57+1)=4.18E+12
      convert(58+1)=1.39E+12; convert(59+1)=3.34E+12
      convert(60+1)=3.14E+12; convert(61+1)=3.34E+12
      convert(62+1)=1.39E+12; convert(63+1)=3.34E+12
      convert(64+1)=8.71E+11; convert(65+1)=2.79E+12
      convert(66+1)=2.79E+12; convert(67+1)=2.79E+12
      convert(68+1)=2.79E+12; convert(69+1)=3.34E+12
      convert(70+1)=2.79E+12; convert(71+1)=3.34E+12
      convert(72+1)=8.71E+11; convert(73+1)=2.79E+12

! Isoprene is traced in terms of equivalent C atoms.
! Compute the baseline ISOPRENE emissions, which depend on veg type   
! 12.d-3 is the carbon mol wt. in kg/mole.
! 1.d4 because the convert data from file is /1000: greg says "huh?"

      factor = 12.d-3*1.d4/avog      

      do J=J_0,J_1
        do I=1,imaxj(J)
          do k=1,ijreg(i,j)
            baseisop(i,j,K) = 
     &      convert(ijland(i,j,k)+1)*factor*dxyp(j)*base_isopreneX
          enddo                                                       
        enddo                                                          
      enddo                                                             
                                                                        
      return                                                            
      end                                                               
                      

! Local Air temperature correction for isoprene emissions:
! From Guenther et al. 1992
	
      real*8 function tcorr(temp)

      implicit none

      real*8 :: r,ct1,ct2,t1,t3,temp

      r=8.314d0
      ct1=9.5d4
      ct2=2.3d5
      t1=3.03e2
      t3=3.14e2

      tcorr =    exp( ct1/(r*t1*temp)*(temp-t1) ) /
     &   (1.D0 + exp( ct2/(r*t1*temp)*(temp-t3) ))

      return
      end
